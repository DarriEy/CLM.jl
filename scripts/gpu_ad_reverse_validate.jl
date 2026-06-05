# ==========================================================================
# gpu_ad_reverse_validate.jl — REVERSE-mode (Enzyme) AD-through-kernel validation.
#
# Phase C bring-up. This harness validates that Enzyme reverse-mode AD
# differentiates cleanly THROUGH the project's KernelAbstractions kernels, for the
# representative kernel shapes the driver is built from:
#
#   1. elementwise per-column kernel          (compute_forc_q!)
#   2. loop-carried batched linear solve      (tridiagonal_multi!, in-kernel Thomas)
#   3. atomic patch→column scatter            (_scatter_add!, Atomix.@atomic path)
#   4. grouped device-view struct kernel       (the immutable struct-of-arrays bundle
#                                              used by every field-heavy driver loop)
#
# TWO RUN MODES — the SAME probes, parameterized by a device-transfer function:
#
#   • CPU  (the default here, on Apple-Silicon dev machines): probes run on the KA
#     CPU backend (plain Arrays). Validation = CROSS-AD: reverse-mode gradient must
#     equal the forward-mode (ForwardDiff) gradient. This PROVES the kernels are
#     reverse-differentiable; it is fully validatable WITHOUT a GPU.
#
#   • CUDA (when an NVIDIA box is available): each probe ALSO runs on CuArrays and
#     adds a PARITY check — device-reverse gradient vs CPU-reverse gradient. Reverse
#     through the device kernel launch is the Enzyme-over-CUDA path; Enzyme does NOT
#     support the Apple Metal GPU (it hangs compiling the adjoint), so the device leg
#     is CUDA-only and was NOT exercised on the Metal dev hardware.
#
# ALSO NOTE: a full-driver Enzyme pass is separately blocked on a Julia-1.12 codegen
# bug ("instruction does not dominate all uses"); these single-kernel probes are
# deliberately small to sidestep it. Revisit the full driver on a Julia LTS.
#
# LESSONS BAKED IN (cost real debugging during bring-up):
#   • Scalar-RETURNING differentiated fns hit Enzyme's "Duplicated Returns not yet
#     handled" — use the IN-PLACE form: write into a Duplicated output array, seed its
#     shadow = 1, read the input shadow as the vjp (= ∂(sum out)/∂input).
#   • Grouped device-view structs bundle ACTIVE arrays (state) with CONST arrays
#     (params/indices) in one immutable. Enzyme's static activity analysis rejects that
#     ("Constant memory is stored to a differentiable variable") → run under
#     Enzyme.set_runtime_activity(Reverse). Used for ALL probes (harmless when unneeded).
#   • A device-view struct that may hold BOTH Dual (AD) and Float64 (param) arrays needs
#     a SEPARATE type param per field — a single {V} can't unify Dual vs Float64 (this
#     also bites the ForwardDiff cross-check, not just Enzyme). See _ProbeStruct{Vx,Vy}.
#   • _scatter_add!'s Atomix.@atomic path (taken for plain Float/Int arrays, i.e. the
#     reverse-mode arrays) reverse-differentiates correctly on CPU — no special handling.
#
# RUN:
#   CPU  (validatable here):  julia --project=scripts scripts/gpu_ad_reverse_validate.jl
#   CUDA (device parity):     julia --project=scripts scripts/gpu_ad_reverse_validate.jl
#                             (auto-adds the device leg when CUDA is functional)
# ==========================================================================

using CLM
using Printf
using Random
using ForwardDiff
import Enzyme
import KernelAbstractions as KA
using KernelAbstractions: @kernel, @index, @Const
const Adapt = CLM.Adapt   # CLM depends on Adapt; reuse it (not a direct scripts dep)

# Reverse-mode with runtime activity — REQUIRED for the grouped-struct probe (mixed
# active/const fields), harmless for the others.
const REV = Enzyme.set_runtime_activity(Enzyme.Reverse)

maxabs(x) = maximum(abs.(x); init = 0.0)
const CROSS_TOL  = 1e-9   # reverse-vs-forward gradient, same precision
const PARITY_TOL = 1e-9   # device-reverse vs CPU-reverse (CUDA carries Float64)

# --- optional CUDA backend (device-reverse leg only; CPU leg always runs) --------
function _cuda_or_nothing()
    try
        @eval using CUDA
        mod = Base.invokelatest(getfield, @__MODULE__, :CUDA)
        Base.invokelatest(() -> mod.functional()) ? mod : nothing
    catch
        nothing
    end
end

# ==========================================================================
# Probe 4 needs its kernel + struct at top level. {Vx,Vy}: x may be Dual (active)
# while y stays Float64 (const) — separate params so both unify (the device-view-
# struct template; mirrors CfEn*/Swm*/Pcb* in src/).
# ==========================================================================
struct _ProbeStruct{Vx,Vy}
    x::Vx
    y::Vy
end
Adapt.@adapt_structure _ProbeStruct

@kernel function _probe_struct_kernel!(out, @Const(inb))
    c = @index(Global)
    @inbounds out[c] = inb.x[c] * inb.y[c] + inb.x[c]^2
end
function _probe_struct_run!(out, x, y)
    inb = _ProbeStruct(x, y)
    backend = CLM._kernel_backend(out)
    _probe_struct_kernel!(backend)(out, inb; ndrange = length(out))
    KA.synchronize(backend)
    nothing
end

@kernel function _probe_scatter_kernel!(col_out, @Const(val), @Const(p2c))
    p = @index(Global)
    @inbounds CLM._scatter_add!(col_out, p2c[p], val[p])
end
function _probe_scatter_run!(col_out, val, p2c)
    backend = CLM._kernel_backend(col_out)
    fill!(col_out, zero(eltype(col_out)))
    _probe_scatter_kernel!(backend)(col_out, val, p2c; ndrange = length(val))
    KA.synchronize(backend)
    nothing
end

# ==========================================================================
# Each probe returns (name, cross, parity, ok). `dev` moves arrays to the device
# (identity on CPU). `cu` is the CUDA module or nothing — when present we also run
# the device-reverse leg and report parity(device-rev vs cpu-rev).
# The reverse vjp uses the IN-PLACE form: Duplicated(out, ones) seeds the cotangent,
# the input shadow accumulates ∂(sum out)/∂input == ForwardDiff.gradient(sum∘f).
# ==========================================================================

# -- Probe 1: compute_forc_q! — elementwise per-column kernel -----------------
function probe_forc_q(cu)
    rng = MersenneTwister(1); nc = 32
    gc = collect(1:nc)
    vp = 500.0 .+ 600.0 .* rand(rng, nc)
    pb = 80000.0 .+ 10000.0 .* rand(rng, nc)
    runner(out, vpx, gcx, pbx) = (CLM.compute_forc_q!(out, gcx, vpx, pbx); nothing)

    rev(dev) = begin
        out = dev(zeros(nc)); dout = dev(ones(nc))
        vpx = dev(copy(vp)); dvp = dev(zeros(nc))
        Enzyme.autodiff(REV, runner, Enzyme.Const,
            Enzyme.Duplicated(out, dout), Enzyme.Duplicated(vpx, dvp),
            Enzyme.Const(dev(gc)), Enzyme.Const(dev(pb)))
        Array(dvp)
    end
    fwd = ForwardDiff.gradient(v -> (o = zeros(eltype(v), nc); CLM.compute_forc_q!(o, gc, v, pb); sum(o)), vp)
    _finish("compute_forc_q! (∂/∂vp)", rev, fwd, cu)
end

# -- Probe 2: tridiagonal_multi! — loop-carried batched Thomas solve ----------
function probe_tridiagonal(cu)
    rng = MersenneTwister(2); nc = 6; nl = 8
    a = zeros(nc, nl); b = zeros(nc, nl); c = zeros(nc, nl); r = zeros(nc, nl)
    for col in 1:nc, j in 1:nl
        a[col, j] = j == 1  ? 0.0 : -0.5 - 0.1 * rand(rng)
        c[col, j] = j == nl ? 0.0 : -0.5 - 0.1 * rand(rng)
        b[col, j] = 2.0 + rand(rng)          # diagonally dominant → well-conditioned
        r[col, j] = rand(rng)
    end
    jtop = fill(1, nc); mask = trues(nc)
    runner(u, rx, ax, bx, cx, jt, mk) = (CLM.tridiagonal_multi!(u, ax, bx, cx, rx, jt, mk, nc, nl); nothing)

    rev(dev) = begin
        u = dev(zeros(nc, nl)); du = dev(ones(nc, nl))
        rx = dev(copy(r)); dr = dev(zeros(nc, nl))
        Enzyme.autodiff(REV, runner, Enzyme.Const,
            Enzyme.Duplicated(u, du), Enzyme.Duplicated(rx, dr),
            Enzyme.Const(dev(a)), Enzyme.Const(dev(b)), Enzyme.Const(dev(c)),
            Enzyme.Const(dev(jtop)), Enzyme.Const(dev(mask)))
        Array(dr)
    end
    fwd = ForwardDiff.gradient(rv -> begin
        uu = zeros(eltype(rv), nc, nl)
        CLM.tridiagonal_multi!(uu, a, b, c, reshape(rv, nc, nl), jtop, mask, nc, nl)
        sum(uu)
    end, vec(r))
    _finish("tridiagonal_multi! (∂/∂r)", dev -> vec(rev(dev)), fwd, cu)
end

# -- Probe 3: _scatter_add! — atomic patch→column scatter ---------------------
function probe_scatter(cu)
    rng = MersenneTwister(3); np = 24; ncol = 5
    p2c = rand(rng, 1:ncol, np)
    val = rand(rng, np)

    rev(dev) = begin
        co = dev(zeros(ncol)); dco = dev(ones(ncol))
        vx = dev(copy(val)); dv = dev(zeros(np))
        Enzyme.autodiff(REV, _probe_scatter_run!, Enzyme.Const,
            Enzyme.Duplicated(co, dco), Enzyme.Duplicated(vx, dv), Enzyme.Const(dev(p2c)))
        Array(dv)
    end
    fwd = ForwardDiff.gradient(v -> (c = zeros(eltype(v), ncol); _probe_scatter_run!(c, v, p2c); sum(c)), val)
    _finish("_scatter_add! (∂/∂val, atomic)", rev, fwd, cu)
end

# -- Probe 4: grouped device-view struct kernel -------------------------------
function probe_struct(cu)
    rng = MersenneTwister(4); nc = 16
    x = rand(rng, nc); y = rand(rng, nc)

    rev(dev) = begin
        out = dev(zeros(nc)); dout = dev(ones(nc))
        xx = dev(copy(x)); dx = dev(zeros(nc))
        Enzyme.autodiff(REV, _probe_struct_run!, Enzyme.Const,
            Enzyme.Duplicated(out, dout), Enzyme.Duplicated(xx, dx), Enzyme.Const(dev(y)))
        Array(dx)
    end
    fwd = ForwardDiff.gradient(v -> (o = zeros(eltype(v), nc); _probe_struct_run!(o, v, y); sum(o)), x)
    _finish("grouped device-view struct (∂/∂x)", rev, fwd, cu)
end

# Run the CPU-reverse leg (always) + the CUDA-reverse leg (if available); compare
# cross-AD (reverse vs forward) and, on CUDA, parity (device-rev vs cpu-rev).
function _finish(name, rev, fwd, cu)
    cpu = rev(identity)
    cross = maxabs(cpu .- fwd)
    parity = NaN
    if cu !== nothing
        cu_xfer(x) = Base.invokelatest(() -> cu.cu(x))
        dev = rev(cu_xfer)
        parity = maxabs(dev .- cpu)
    end
    ok = cross < CROSS_TOL && (isnan(parity) || parity < PARITY_TOL)
    (name, cross, parity, ok)
end

function main()
    println("=" ^ 78)
    println("REVERSE-mode (Enzyme) AD-through-kernel validation")
    cu = _cuda_or_nothing()
    if cu === nothing
        println("Mode: CPU (KA CPU backend) — validating reverse-vs-forward cross-AD.")
        println("      No functional CUDA; device-reverse parity leg SKIPPED (Enzyme has no")
        println("      Metal support — the device leg is CUDA-only, run it on an NVIDIA box).")
    else
        println("Mode: CPU + CUDA — cross-AD on CPU AND device-reverse-vs-cpu-reverse parity.")
    end
    println("=" ^ 78)

    probes = (probe_forc_q, probe_tridiagonal, probe_scatter, probe_struct)
    results = Any[]
    for p in probes
        try
            push!(results, p(cu))
        catch err
            push!(results, ("$(p) ERRORED: $(typeof(err))", NaN, NaN, false))
            @warn "reverse-AD probe threw" probe=p exception=(err, catch_backtrace())
        end
    end

    nfail = 0
    for (nm, cross, parity, ok) in results
        ptxt = isnan(parity) ? "   n/a    " : @sprintf("%.3e", parity)
        @printf("  [%s] %-36s  cross(rev-vs-fwd)=%.3e  parity(dev-vs-cpu)=%s\n",
                ok ? "PASS" : "FAIL", nm, cross, ptxt)
        ok || (nfail += 1)
    end
    println()
    if nfail == 0
        println("  ✓ Reverse-mode AD differentiates through all kernel shapes (CPU).",
                cu === nothing ? " Device parity pending CUDA." : " Device parity ✓.")
    else
        println("  ✗ $nfail probe(s) failed — see warnings above.")
    end
    return nfail == 0 ? 0 : 1
end

exit(main())
