# ==========================================================================
# gpu_validate.jl — real-hardware GPU validation for CLM.jl's kernels.
#
# This machine (Apple Silicon) has only a Metal GPU, which is Float32-only, while
# CLM runs in Float64 — so this script CANNOT run meaningfully here. It is written
# to run unchanged on a CUDA or AMD (ROCm) box: it auto-detects a Float64-capable
# backend, and if none is present it prints a notice and exits 0 (no-op).
#
#   # on a CUDA box:
#   julia --project=. -e 'using CUDA' -e 'include("scripts/gpu_validate.jl")'
#   # or simply (auto-detect):
#   julia --project=. scripts/gpu_validate.jl
#
# WHAT IT VALIDATES
# -----------------
# Per-kernel CPU-vs-device parity for the Phase-4 kernels: each kernel is run on a
# CPU `Array` (the correctness reference) and on a device array of the SAME inputs,
# and the results are compared to a tight tolerance. This is the honest scope today:
# the full clm_drv! timestep is NOT yet 100% kernelized (the ~156 element-wise ops +
# the 4 iterative solvers + the 2 batched linear solvers are; many other per-column
# loops still index scalars and would trip GPU scalar-indexing), so a full-driver GPU
# run is deferred. As more loops are kernelized, add them to KERNEL_TESTS below.
#
# The kernels themselves are backend-agnostic (they dispatch on the output array's
# backend via `_launch!`), so "running on GPU" here means: construct the inputs as
# device arrays and call the SAME helper the driver calls.
# ==========================================================================

using CLM
using Printf
using Random

# --- backend auto-detection -------------------------------------------------
# Returns (name, to_device) where to_device(::Array) -> device array, or nothing.
function detect_backend()
    # CUDA
    try
        @eval using CUDA
        if Base.invokelatest(CUDA.functional)
            return ("CUDA", x -> Base.invokelatest(CUDA.cu, x))   # cu() keeps Float64
        end
    catch
    end
    # AMDGPU (ROCm)
    try
        @eval using AMDGPU
        if Base.invokelatest(AMDGPU.functional)
            return ("AMDGPU", x -> Base.invokelatest(AMDGPU.ROCArray, x))
        end
    catch
    end
    return nothing
end

# Device arrays disallow BitArray and require Bool masks as Vector{Bool}.
_dev(to, x::BitArray) = to(collect(Bool, x))
_dev(to, x::AbstractArray) = to(x)

maxabsdiff(a, b) = maximum(abs.(Array(a) .- Array(b)); init = 0.0)

# --- the kernel parity tests ------------------------------------------------
# Each returns (name, max_abs_diff, ok). Inputs are built on CPU; the reference is
# the CPU kernel; the device run uses the SAME helper on device arrays.

function test_forc_q(to)
    rng = MersenneTwister(1)
    nc = 64
    gridcell = collect(1:nc)
    vp = 500.0 .+ 600.0 .* rand(rng, nc)
    pbot = 80000.0 .+ 10000.0 .* rand(rng, nc)
    # CPU reference
    out_cpu = zeros(nc)
    CLM.compute_forc_q!(out_cpu, gridcell, vp, pbot)
    # device
    out_dev = to(zeros(nc))
    CLM.compute_forc_q!(out_dev, _dev(to, gridcell), _dev(to, vp), _dev(to, pbot))
    d = maxabsdiff(out_cpu, out_dev)
    return ("compute_forc_q!", d, d < 1e-12)
end

function test_tridiagonal_multi(to)
    rng = MersenneTwister(2)
    ncols, nlevs = 48, 20
    # Diagonally dominant tridiagonal systems, one per column.
    a = -0.5 .* ones(ncols, nlevs)
    c = -0.5 .* ones(ncols, nlevs)
    b = 2.0 .+ 0.5 .* rand(rng, ncols, nlevs)
    r = rand(rng, ncols, nlevs)
    a[:, 1] .= 0.0
    c[:, nlevs] .= 0.0
    jtop = ones(Int, ncols)
    mask = trues(ncols)
    u_cpu = zeros(ncols, nlevs)
    CLM.tridiagonal_multi!(u_cpu, a, b, c, r, jtop, mask, ncols, nlevs)
    u_dev = to(zeros(ncols, nlevs))
    CLM.tridiagonal_multi!(u_dev, _dev(to, a), _dev(to, b), _dev(to, c), _dev(to, r),
                           jtop, _dev(to, mask), ncols, nlevs)
    d = maxabsdiff(u_cpu, u_dev)
    return ("tridiagonal_multi! (batched Thomas)", d, d < 1e-10)
end

function test_batched_band(to)
    rng = MersenneTwister(3)
    nc = 32
    kl = ku = 2
    nband = 2kl + 1          # 5
    nlevsno = 5
    nlev = 25                # tvector second dim (snow+soil padded)
    # Each column: full soil column jtop=1..jbot (no snow layers -> snl=0 path).
    jtop = ones(Int, nc)               # 1-based fortran top
    jbot = fill(15, nc)                # 15 active soil layers
    mask = trues(nc)
    # bmatrix[c, band, level]: diagonally dominant pentadiagonal.
    bmatrix = zeros(nc, nband, nlev)
    for c in 1:nc, lev in 1:nlev
        bmatrix[c, kl + 1, lev] = 4.0 + rand(rng)      # diagonal band
        bmatrix[c, kl,     lev] = -0.7                 # sub
        bmatrix[c, kl + 2, lev] = -0.7                 # super
        bmatrix[c, 1,      lev] = -0.3                 # sub2
        bmatrix[c, nband,  lev] = -0.3                 # super2
    end
    rvector = rand(rng, nc, nlev)
    t_cpu = zeros(nc, nlev)
    t_dev = to(zeros(nc, nlev))
    CLM.batched_band_solve!(t_cpu, bmatrix, rvector, jtop, jbot, mask, kl, ku, nlevsno)
    CLM.batched_band_solve!(t_dev, _dev(to, bmatrix), _dev(to, rvector), jtop, jbot,
                            _dev(to, mask), kl, ku, nlevsno)
    d = maxabsdiff(t_cpu, t_dev)
    return ("batched_band_solve! (pentadiagonal GE)", d, d < 1e-9)
end

const KERNEL_TESTS = [test_forc_q, test_tridiagonal_multi, test_batched_band]

# --- driver -----------------------------------------------------------------
function main()
    println("=" ^ 70)
    println("GPU VALIDATION for CLM.jl kernels (CPU reference vs device)")
    println("=" ^ 70)

    backend = detect_backend()
    if backend === nothing
        println("""
  No Float64-capable GPU backend detected (CUDA / AMDGPU).
  This is expected on Apple Silicon (Metal is Float32-only) and on CPU-only CI.
  Nothing to validate here — the kernels are exercised on the KA CPU backend by
  the test suite. Run this script on a CUDA or ROCm box to validate on real GPU.
""")
        return 0
    end

    name, to = backend
    @printf("  Detected backend: %s\n\n", name)

    npass = 0
    nfail = 0
    for t in KERNEL_TESTS
        knm, d, ok = try
            t(to)
        catch e
            (string(t), NaN, false)
        end
        status = ok ? "PASS" : "FAIL"
        @printf("  [%s] %-42s  max|dev-cpu| = %.3e\n", status, knm, d)
        ok ? (npass += 1) : (nfail += 1)
    end

    println()
    @printf("  %d passed, %d failed\n", npass, nfail)
    if nfail == 0
        println("\n  ALL KERNELS MATCH CPU ON $(name) ✓")
        println("  (Full clm_drv! GPU run is deferred until all per-column loops are kernelized.)")
    else
        println("\n  SOME KERNELS DIVERGED — investigate before trusting the GPU path.")
    end
    return nfail == 0 ? 0 : 1
end

exit(main())
