# ==========================================================================
# gpu_ad_reverse_validate.jl — REVERSE-mode (Enzyme) AD-on-GPU scaffolding.
#
# ⚠️  UNVALIDATED SCAFFOLDING (2026-06-01). This script has NOT been run — the
#     development machine is Apple Silicon (Metal), and Enzyme does NOT
#     differentiate the Metal device kernel launch (no Apple-GPU Enzyme support;
#     it hangs compiling the adjoint). Reverse-mode AD through a KernelAbstractions
#     kernel ON the GPU is an Enzyme-over-CUDA path. This file encodes the intended
#     structure so it can be run + debugged on a CUDA (NVIDIA) box; expect to
#     iterate. It is written to mirror the FORWARD harness gpu_ad_validate.jl.
#
# ALSO GATED: full-driver Enzyme is blocked on a Julia-1.12 codegen bug
#     ("instruction does not dominate all uses"); these single-kernel probes are
#     deliberately small to sidestep that. Revisit on a Julia LTS with the fix.
#
# WHAT IT CHECKS (when run on CUDA):
#   For a kernel out = f(x), seed the output cotangent d̄out = 1 and accumulate the
#   input gradient d̄x via Enzyme.Reverse through the on-device kernel launch, then:
#     1. PARITY   — device-reverse d̄x vs CPU-reverse d̄x (does reverse-AD give the
#                   SAME answer on the GPU as the CPU?). PASS/FAIL on this.
#     2. CROSS-AD — device-reverse d̄x vs FORWARD-mode partials (for element-wise
#                   kernels with d̄out=1, the vjp equals the forward partial).
#
# RUN (on a CUDA machine):
#   julia --project=scripts scripts/gpu_ad_reverse_validate.jl
# ==========================================================================

using CLM
using Printf
using Random
using ForwardDiff
import Enzyme
import KernelAbstractions as KA

# --- CUDA-only backend detection (reverse-AD on GPU is the Enzyme+CUDA path) ----
# Deliberately NOT reusing gpu_backends.jl's Metal/AMD fallbacks: Enzyme reverse
# through the device launch is only expected to work on CUDA.
function _cuda_or_nothing()
    try
        @eval using CUDA
        mod = Base.invokelatest(getfield, @__MODULE__, :CUDA)
        Base.invokelatest(() -> mod.functional()) ? mod : nothing
    catch
        nothing
    end
end

maxabs(x) = maximum(abs.(x); init = 0.0)
relerr(a, b) = maxabs(a .- b) / max(maxabs(b), eps(Float64))
const PARITY_TOL = 1e-9      # CUDA carries Float64 — device-rev vs CPU-rev should be tight
const CROSS_TOL  = 1e-9      # vjp(d̄out=1) vs forward partial, same precision

# ---------------------------------------------------------------------------
# Test 1: compute_forc_q!(out, gridcell, vp, pbot) — element-wise, no solver.
# vjp w.r.t. vp with d̄out=1 equals d(out_i)/d(vp_i) (forward partial), since the
# kernel is element-wise. This is the smallest end-to-end reverse-AD-on-GPU probe.
# ---------------------------------------------------------------------------
function reverse_probe_forc_q(cu)
    rng = MersenneTwister(1); nc = 32
    gc  = collect(1:nc)
    vp  = 500.0 .+ 600.0 .* rand(rng, nc)
    pb  = 80000.0 .+ 10000.0 .* rand(rng, nc)

    # f!(out, vp): launch the kernel writing out from vp (pbot/gc are Const).
    runner(out, vpx, gcx, pbx) = (CLM.compute_forc_q!(out, gcx, vpx, pbx); nothing)

    # ---- CPU reverse (reference) ----
    out_c  = zeros(nc);  dout_c = ones(nc)         # seed cotangent = 1
    vp_c   = copy(vp);   dvp_c  = zeros(nc)
    Enzyme.autodiff(Enzyme.Reverse, runner,
        Enzyme.Const,
        Enzyme.Duplicated(out_c, dout_c),
        Enzyme.Duplicated(vp_c, dvp_c),
        Enzyme.Const(gc), Enzyme.Const(pb))

    # ---- device reverse (the actual question) ----
    out_d  = cu(zeros(nc));  dout_d = cu(ones(nc))
    vp_d   = cu(copy(vp));   dvp_d  = cu(zeros(nc))
    gc_d   = cu(gc);         pb_d   = cu(pb)
    Enzyme.autodiff(Enzyme.Reverse, runner,
        Enzyme.Const,
        Enzyme.Duplicated(out_d, dout_d),
        Enzyme.Duplicated(vp_d, dvp_d),
        Enzyme.Const(gc_d), Enzyme.Const(pb_d))

    # ---- forward-mode cross-check (partial of each out_i w.r.t its vp_i) ----
    DT = ForwardDiff.Dual{Nothing,Float64,1}
    vpf = DT.(vp); for i in eachindex(vpf); vpf[i] = ForwardDiff.Dual(vp[i], 1.0); end
    outf = zeros(DT, nc); CLM.compute_forc_q!(outf, gc, vpf, DT.(pb))
    fwd = ForwardDiff.partials.(outf, 1)

    parity = maxabs(dvp_c .- Array(dvp_d))
    cross  = maxabs(Array(dvp_d) .- fwd)
    ok = parity < PARITY_TOL && cross < CROSS_TOL
    return ("compute_forc_q! (vjp wrt vp)", parity, cross, ok)
end

# ---------------------------------------------------------------------------
# Add further probes here as reverse-AD-on-CUDA is brought up, e.g.:
#   - tridiagonal_multi! (batched Thomas solve — loop-carried; tests reverse
#     through an in-kernel sequential sweep)
#   - a soil_temperature rhs/mat kernel (struct-arg; tests reverse through the
#     grouped device-view-struct pattern)
# Each follows the same shape: Duplicated(primal, shadow) for differentiated
# arrays, Const for indices/params, seed the output shadow, compare device-rev
# vs CPU-rev (parity) and vs forward partials (cross-AD).
# ---------------------------------------------------------------------------

function main()
    println("=" ^ 70)
    println("REVERSE-mode (Enzyme) AD-on-GPU validation  [UNVALIDATED SCAFFOLDING]")
    println("=" ^ 70)
    cu_mod = _cuda_or_nothing()
    if cu_mod === nothing
        println("""
  No functional CUDA backend found. Reverse-mode AD on GPU requires NVIDIA CUDA
  (Enzyme does not support the Apple/Metal GPU). This script is scaffolding only
  and has NOT been validated; run it on a CUDA box (and a Julia LTS that resolves
  the 1.12 'instruction does not dominate all uses' codegen bug). Nothing to do
  here — exiting cleanly.""")
        return 0
    end
    cu(x) = Base.invokelatest(() -> cu_mod.cu(x))
    @printf("  Backend: CUDA   (working precision: Float64)\n\n")

    results = Any[]
    for probe in (reverse_probe_forc_q,)
        try
            push!(results, probe(cu))
        catch err
            push!(results, ("(probe errored: $(typeof(err)))", NaN, NaN, false))
            @warn "reverse-AD probe threw — expected during bring-up" exception=err
        end
    end

    nfail = 0
    for (nm, parity, cross, ok) in results
        @printf("  [%s] %-28s  parity(dev-rev vs cpu-rev)=%.3e  cross(vs fwd)=%.3e\n",
                ok ? "PASS" : "FAIL", nm, parity, cross)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  Reverse-mode AD matches on CUDA ✓ (validate before trusting)" :
                         "  Reverse-mode AD NOT matching — debug on CUDA hardware.")
    return nfail == 0 ? 0 : 1
end

exit(main())
