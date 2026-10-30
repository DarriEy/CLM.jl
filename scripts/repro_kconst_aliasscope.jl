# Minimal reproduction of the #365 CNFUN clamp miscompile pattern.
# Needs only KernelAbstractions. Run: julia --project=. scripts/repro_kconst_aliasscope.jl
# Expected: (1) 0.0 with KA >= 0.9.45 or Julia 1.10; ~10 with Julia 1.12.6 + KA 0.9.42
# (LLVM 18.1.7). (2) is the upstream reduction (JuliaLang/julia#63129); it is 0.0 on
# 1.12.6 and fails on 1.13 -- shown for reference.
using KernelAbstractions
using Base.Experimental: @aliasscope
println("Julia ", VERSION, "  KA ", pkgversion(KernelAbstractions))

# (1) The CNFUN pattern: two @view slices of one scratch array, an in-place clamp,
#     then a read-back accumulation. `lim` is the only @Const argument.
@kernel function clamp_kernel!(S, out, @Const(lim), n)
    p = @index(Global)
    @inbounds begin
        a = @view S[:, :, 1]; m = @view S[:, :, 2]
        acc = zero(eltype(S))
        for j in 1:n
            if a[p, j] + m[p, j] > lim[p, j]
                d = a[p, j] + m[p, j] - lim[p, j]
                t = a[p, j]
                a[p, j] -= d * a[p, j] / (a[p, j] + m[p, j])
                m[p, j] -= d * m[p, j] / (t + m[p, j])
            end
            acc += m[p, j]          # must read the CLAMPED value
        end
        out[p] = acc
    end
end
function clamp_ref!(S, out, lim, n)
    for p in axes(S, 1)
        acc = 0.0
        for j in 1:n
            a, m = S[p, j, 1], S[p, j, 2]
            if a + m > lim[p, j]
                d = a + m - lim[p, j]
                S[p, j, 1] = a - d * a / (a + m)
                S[p, j, 2] = m - d * m / (a + m)
            end
            acc += S[p, j, 2]
        end
        out[p] = acc
    end
    out
end

# (2) Upstream reduction (JuliaLang/julia#63129): no Const loads at all.
function fwd!(B, a, n)
    for l in axes(B, 1)
        @aliasscope begin
            @inbounds for i in 2:n
                B[l, i] -= a[i] * B[l, i-1]
            end
        end
    end
    B
end
function fwd_plain!(B, a, n)
    for l in axes(B, 1), i in 2:n
        B[l, i] -= a[i] * B[l, i-1]
    end
    B
end

# (data: plain rand(); the error is large for any draw where the clamp fires)
np, n = 8, 20
S = rand(np, n, 2); lim = 0.5 .* rand(np, n)
o1 = zeros(np); o2 = zeros(np)
S1 = copy(S); clamp_kernel!(CPU())(S1, o1, lim, n; ndrange = np); KernelAbstractions.synchronize(CPU())
clamp_ref!(copy(S), o2, lim, n)
println("(1) CNFUN clamp pattern, KA @Const kernel : max|err| = ", maximum(abs.(o1 .- o2)))
B = rand(1, 10); a = rand(10)
println("(2) @aliasscope recurrence, no Const     : max|err| = ",
        maximum(abs.(fwd!(copy(B), a, 10) .- fwd_plain!(copy(B), a, 10))))
