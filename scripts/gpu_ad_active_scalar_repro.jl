# Minimal reproduction of the CUDA whole-driver reverse-AD blocker, and a check of
# CLM's workaround.
#
# A KernelAbstractions kernel takes a Float64 scalar argument and is launched on
# CUDA inside a function that Enzyme reverse-differentiates (runtime activity, as in
# CLM's compositional_reverse!). Enzyme annotates the scalar Active and
# KernelAbstractions' Enzyme extension then errors with
#     "Active kernel arguments not supported on GPU"
# whether the scalar is read from a global Ref, read from a Const argument, or
# passed through `ignore_derivatives`. Only a literal that LLVM constant-folds gets
# through. The same code differentiates correctly on the CPU backend.
#
#   variant :globalref  scalar read from a global Ref (like PHASE_CHANGE_MASS_K[])
#   variant :constaux   scalar read from a Const NamedTuple argument (like aux.dtime)
#   variant :ignored    scalars passed through EnzymeCore.ignore_derivatives
#   variant :literal    scalars written as literals at the call site
#   variant :clm        launched through CLM._kcall! (the launch-level rule)
#
# Observed on Julia 1.12.6, Enzyme 0.13.216, EnzymeCore 0.8.26,
# KernelAbstractions 0.9.45, CUDA.jl 5.11.3 (RTX 5070 Ti, sm_120).
#
# Usage: julia --project=scripts scripts/gpu_ad_active_scalar_repro.jl [variant ...]
# Exit status 0 iff the :clm variant matches the exact gradient on CUDA.
using CLM, Enzyme, KernelAbstractions, CUDA, Printf
const EnzymeCore = Enzyme.EnzymeCore
const KA = KernelAbstractions

const KREF = Ref(0.25)

@kernel function axpy2_kernel!(y, @Const(x), a, b)
    i = @index(Global)
    @inbounds y[i] = a * x[i] + b * x[i]^2
end

launch!(y, x, a, b) = (axpy2_kernel!(KA.get_backend(y))(y, x, a, b; ndrange = length(y));
                       KA.synchronize(KA.get_backend(y)); nothing)

phase_globalref!(st, aux) = launch!(st.y, st.x, convert(Float64, aux.a), KREF[])
phase_constaux!(st, aux)  = launch!(st.y, st.x, convert(Float64, aux.a), aux.b)
phase_ignored!(st, aux)   = launch!(st.y, st.x,
    EnzymeCore.ignore_derivatives(convert(Float64, aux.a)),
    EnzymeCore.ignore_derivatives(KREF[]))
phase_literal!(st, aux)   = launch!(st.y, st.x, 1.5, 0.25)
phase_clm!(st, aux) = CLM._kcall!(axpy2_kernel!, KA.get_backend(st.y), st.y, st.x,
    convert(Float64, aux.a), KREF[]; ndrange = length(st.y))

const VARIANTS = (; globalref = phase_globalref!, constaux = phase_constaux!,
                    ignored = phase_ignored!, literal = phase_literal!, clm = phase_clm!)

function run(variant, AT)
    f = VARIANTS[variant]
    x = AT(collect(1.0:8.0)); y = AT(zeros(8))
    st = (; x, y); dst = (; x = AT(zeros(8)), y = AT(ones(8)))
    aux = (; a = 1.5, b = 0.25)
    Enzyme.API.strictAliasing!(false)
    mode = Enzyme.set_runtime_activity(Enzyme.Reverse)
    Enzyme.autodiff(mode, f, Enzyme.Const, Enzyme.Duplicated(st, dst), Enzyme.Const(aux))
    xe = collect(1.0:8.0)
    err = maximum(abs, Array(dst.x) .- (1.5 .+ 2 .* 0.25 .* xe))
    # y is overwritten, so its incoming cotangent must be consumed (zeroed).
    err = max(err, maximum(abs, Array(dst.y)))
    @printf("%-10s %-5s max|grad - exact| = %.3e\n", variant, AT === Array ? "CPU" : "CUDA", err)
    return err
end

function main()
    variants = isempty(ARGS) ? collect(keys(VARIANTS)) : Symbol.(ARGS)
    println("Julia ", VERSION, ", Enzyme ", pkgversion(Enzyme), ", KernelAbstractions ",
            pkgversion(KA), ", CUDA.jl ", pkgversion(CUDA))
    ok = true
    for v in variants, AT in (Array, CuArray)
        try
            err = run(v, AT)
            v === :clm && (ok &= err == 0)
        catch e
            @printf("%-10s %-5s ERROR: %s\n", v, AT === Array ? "CPU" : "CUDA",
                    first(split(sprint(showerror, e), '\n')))
            v === :clm && (ok = false)
        end
    end
    println(ok ? "PASS" : "FAIL")
    return ok ? 0 : 1
end

exit(main())
