# ==========================================================================
# gpu_validate.jl — real-hardware GPU validation for CLM.jl's kernels.
#
# Auto-detects an available GPU backend and runs each Phase-4 kernel on a CPU
# `Array` (the correctness reference) AND on a device array of the SAME inputs at
# the SAME precision, comparing the two to a tolerance. This is a PARITY test:
# "does the device execute the kernel the same way the CPU does?", run at the
# precision the backend actually supports.
#
# SUPPORTED BACKENDS (auto-detected, first functional one wins)
#   CUDA     — Float64   (NVIDIA)
#   AMDGPU   — Float64   (AMD ROCm)
#   oneAPI   — Float64   (Intel; Data-Center GPUs do f64, integrated may not)
#   Metal    — Float32   (Apple Silicon; the hardware has NO Float64 at all)
#
# Float64 backends are checked at a tight tolerance; Float32-only backends
# (Metal) are checked at a relaxed same-precision tolerance — the CPU reference
# is itself run in Float32 so we compare like with like (this is a parity test,
# not an accuracy-vs-Fortran test).
#
# RUN
#   # auto-detect, using the script environment (which carries the GPU packages):
#   julia --project=scripts scripts/gpu_validate.jl
#   # the core package's own env works too if the backend pkg is installed:
#   julia --project=. -e 'using CUDA' -e 'include("scripts/gpu_validate.jl")'
#
# WHAT IS VALIDATED
# -----------------
# Per-kernel CPU-vs-device parity for the Phase-4 kernels currently in the
# validation set (forc_q, batched tridiagonal Thomas, batched pentadiagonal band
# GE). The full clm_drv! timestep is NOT yet 100% kernelized, so a full-driver GPU
# run is still deferred. As more loops are kernelized, add them to KERNEL_TESTS.
# ==========================================================================

using CLM
using Printf
using Random

# Backend auto-detection (detect_backend / _dev) is shared with gpu_ad_validate.jl.
include(joinpath(@__DIR__, "gpu_backends.jl"))

maxabsdiff(a, b) = maximum(abs.(Array(a) .- Array(b)); init = 0.0)

# Per-precision tolerances. Float64 backends must match the CPU almost exactly;
# Float32-only backends (Metal) accumulate FMA/ordering differences through the
# solver sweeps, so they use a relaxed — but still failure-catching — bound.
function _tol(::Type{FT}, which::Symbol) where {FT}
    if FT === Float32
        which === :forcq ? 1f-5 :
        which === :tri   ? 5f-3 :
        which === :band  ? 5f-3 : 1f-3
    else
        which === :forcq ? 1e-12 :
        which === :tri   ? 1e-10 :
        which === :band  ? 1e-9  : 1e-10
    end
end

# --- the kernel parity tests ------------------------------------------------
# Each takes (to, FT) and returns (name, max_abs_diff, ok). Inputs are built at FT
# on the CPU; the reference is the CPU kernel at FT; the device run uses the SAME
# helper on device arrays. Every per-thread-indexed array (incl. integer index and
# Bool mask arrays) is moved to the device — a kernel cannot index a host array.

function test_forc_q(to, ::Type{FT}) where {FT}
    rng = MersenneTwister(1)
    nc = 64
    gridcell = collect(1:nc)
    vp  = FT.(500 .+ 600 .* rand(rng, nc))
    pbot = FT.(80000 .+ 10000 .* rand(rng, nc))
    out_cpu = zeros(FT, nc)
    CLM.compute_forc_q!(out_cpu, gridcell, vp, pbot)
    out_dev = to(zeros(FT, nc))
    CLM.compute_forc_q!(out_dev, _dev(to, gridcell), _dev(to, vp), _dev(to, pbot))
    d = maxabsdiff(out_cpu, out_dev)
    return ("compute_forc_q!", d, d < _tol(FT, :forcq))
end

function test_tridiagonal_multi(to, ::Type{FT}) where {FT}
    rng = MersenneTwister(2)
    ncols, nlevs = 48, 20
    a = FT.(-0.5 .* ones(ncols, nlevs))
    c = FT.(-0.5 .* ones(ncols, nlevs))
    b = FT.(2 .+ 0.5 .* rand(rng, ncols, nlevs))
    r = FT.(rand(rng, ncols, nlevs))
    a[:, 1] .= 0
    c[:, nlevs] .= 0
    jtop = ones(Int, ncols)
    mask = trues(ncols)
    u_cpu = zeros(FT, ncols, nlevs)
    CLM.tridiagonal_multi!(u_cpu, a, b, c, r, jtop, mask, ncols, nlevs)
    u_dev = to(zeros(FT, ncols, nlevs))
    CLM.tridiagonal_multi!(u_dev, _dev(to, a), _dev(to, b), _dev(to, c), _dev(to, r),
                           _dev(to, jtop), _dev(to, mask), ncols, nlevs)
    d = maxabsdiff(u_cpu, u_dev)
    return ("tridiagonal_multi! (batched Thomas)", d, d < _tol(FT, :tri))
end

function test_batched_band(to, ::Type{FT}) where {FT}
    rng = MersenneTwister(3)
    nc = 32
    kl = ku = 2
    nband = 2kl + 1          # 5
    nlevsno = 5
    nlev = 25                # tvector second dim (snow+soil padded)
    jtop = ones(Int, nc)               # 1-based fortran top
    jbot = fill(15, nc)                # 15 active soil layers
    mask = trues(nc)
    bmatrix = zeros(FT, nc, nband, nlev)
    for c in 1:nc, lev in 1:nlev
        bmatrix[c, kl + 1, lev] = FT(4) + rand(rng, FT)  # diagonal band
        bmatrix[c, kl,     lev] = FT(-0.7)               # sub
        bmatrix[c, kl + 2, lev] = FT(-0.7)               # super
        bmatrix[c, 1,      lev] = FT(-0.3)               # sub2
        bmatrix[c, nband,  lev] = FT(-0.3)               # super2
    end
    rvector = FT.(rand(rng, nc, nlev))
    t_cpu = zeros(FT, nc, nlev)
    t_dev = to(zeros(FT, nc, nlev))
    CLM.batched_band_solve!(t_cpu, bmatrix, rvector, jtop, jbot, mask, kl, ku, nlevsno)
    CLM.batched_band_solve!(t_dev, _dev(to, bmatrix), _dev(to, rvector),
                            _dev(to, jtop), _dev(to, jbot), _dev(to, mask), kl, ku, nlevsno)
    d = maxabsdiff(t_cpu, t_dev)
    return ("batched_band_solve! (pentadiagonal GE)", d, d < _tol(FT, :band))
end

const KERNEL_TESTS = [test_forc_q, test_tridiagonal_multi, test_batched_band]

# --- driver -----------------------------------------------------------------
# `backend` is detected at TOP LEVEL (see bottom of file) so the `using <Backend>`
# has advanced the world age before main() runs — otherwise every device method
# (size, get_backend, …) is "too new to be called from this world context".
function main(backend)
    println("=" ^ 70)
    println("GPU VALIDATION for CLM.jl kernels (CPU reference vs device)")
    println("=" ^ 70)

    if backend === nothing
        println("""
  No GPU backend detected (CUDA / AMDGPU / oneAPI / Metal).
  Nothing to validate here — the kernels are exercised on the KA CPU backend by
  the test suite. Install a backend package and re-run on real GPU hardware.
""")
        return 0
    end

    name, to, FT = backend
    @printf("  Detected backend: %s   (working precision: %s)\n", name, FT)
    if FT === Float32
        println("  NOTE: Float32-only GPU — relaxed same-precision parity tolerance.")
    end
    println()

    npass = 0
    nfail = 0
    for t in KERNEL_TESTS
        knm, d, ok = try
            t(to, FT)
        catch e
            (string(t) * "  (ERROR: " * sprint(showerror, e) * ")", NaN, false)
        end
        status = ok ? "PASS" : "FAIL"
        @printf("  [%s] %-42s  max|dev-cpu| = %.3e\n", status, knm, d)
        ok ? (npass += 1) : (nfail += 1)
    end

    println()
    @printf("  %d passed, %d failed\n", npass, nfail)
    if nfail == 0
        println("\n  ALL KERNELS MATCH CPU ON $(name) ($(FT)) ✓")
        println("  (Full clm_drv! GPU run is deferred until all per-column loops are kernelized.)")
    else
        println("\n  SOME KERNELS DIVERGED — investigate before trusting the GPU path.")
    end
    return nfail == 0 ? 0 : 1
end

# Detect (and load) the backend as its own top-level statement, THEN run main in
# the advanced world age where the backend's device methods are callable.
const BACKEND = detect_backend()
exit(main(BACKEND))
