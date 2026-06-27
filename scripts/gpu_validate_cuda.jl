#!/usr/bin/env julia
# ==========================================================================
# gpu_validate_cuda.jl — one-command NVIDIA CUDA backend validation.
#
# The CLM.jl GPU backends are validated at 0.0 parity on Apple Metal; the CUDA
# (and AMDGPU) backends are CPU-proxy-validated only — never exercised on real
# silicon. This script closes that gap: run it on a CUDA box (a rented cloud GPU,
# or the JuliaGPU Buildkite `juliagpu` queue) to assert the CUDA backend is
# functional and matches the CPU reference across the kernel + AD parity suites.
#
#   julia --project=. scripts/gpu_validate_cuda.jl
#
# Exit code 0 = all CUDA validations passed; 1 = a divergence/failure; 2 = no
# functional CUDA device (so CI can distinguish "no GPU" from "GPU mismatch").
#
# What it runs:
#   1. CUDA functional check + device report (name, compute capability, FP64).
#   2. The backend registry round-trip — clm_set_backend(:cuda) + clm_device_array.
#   3. scripts/gpu_validate.jl     — the full kernel-parity suite (backend-agnostic;
#      detect_backend() selects CUDA on a CUDA box), CUDA arrays vs CPU.
#   4. scripts/gpu_ad_validate.jl  — forward-mode AD parity on the device.
#
# (3) and (4) are the proven per-module validation harnesses — this driver just
# asserts CUDA, points them at the device, and aggregates a single verdict, so
# the rented box is exactly `git clone && julia --project=. scripts/gpu_validate_cuda.jl`.
# ==========================================================================

using CLM
using Printf

# --- 1. CUDA functional? ---------------------------------------------------
# `using CUDA` binds in a newer world age than this body, so reach the loaded
# module's methods through invokelatest (same pattern as scripts/gpu_backends.jl).
const _CUDA = try
    @eval using CUDA
    Base.invokelatest(getfield, @__MODULE__, :CUDA)
catch e
    @warn "CUDA.jl could not be loaded (install with Pkg.add(\"CUDA\") on a GPU box)" exception=(e, catch_backtrace())
    nothing
end

if _CUDA === nothing || !Base.invokelatest(() -> _CUDA.functional())
    println("""
    ✗ CUDA is NOT functional in this environment.
      This script must run on a host with an NVIDIA GPU + driver + CUDA.jl
      (a rented cloud GPU — RunPod/Vast.ai/Lambda — or the JuliaGPU Buildkite
      `juliagpu` queue). Nothing to validate here.""")
    exit(2)
end

# --- 2. Device report ------------------------------------------------------
Base.invokelatest() do
    dev = _CUDA.device()
    @printf("  CUDA device       : %s\n", _CUDA.name(dev))
    @printf("  compute capability: %s\n", string(_CUDA.capability(dev)))
    try; @printf("  CUDA runtime      : %s\n", string(_CUDA.runtime_version())); catch; end
    try; @printf("  CUDA driver       : %s\n", string(_CUDA.driver_version()));  catch; end
    # CLM.jl is Float64 throughout; consumer cards do FP64 (slowly) — fine for
    # correctness parity, which is what this validates.
    println("  precision under test: Float64")
end
println()

# --- 3. Backend registry round-trip (the clm_set_backend abstraction) ------
function registry_check()
    try
        CLM.clm_set_backend(:cuda; require_functional = true)
        x = collect(1.0:8.0)
        d = CLM.clm_device_array(x)            # generalizes adapt(CuArray, x)
        back = Array(d)
        CLM.clm_set_backend(:cpu)
        return (back == x, back == x ? "" : "device round-trip mismatch")
    catch e
        try; CLM.clm_set_backend(:cpu); catch; end
        return (false, sprint(showerror, e))
    end
end

# --- 4. Run the proven backend-agnostic suites as subprocesses --------------
# Each in a fresh process so `detect_backend()`'s `using CUDA` advances the world
# age cleanly before the kernels run.
function run_suite(script)
    jl   = Base.julia_cmd()
    proj = abspath(joinpath(@__DIR__, ".."))
    path = joinpath(@__DIR__, script)
    buf  = IOBuffer()
    ok = try
        run(pipeline(`$jl --project=$proj $path`; stdout = buf, stderr = buf))
        true
    catch
        false
    end
    return (ok, String(take!(buf)))
end

results = Tuple{String,Bool}[]

println("=== [1/3] backend registry: clm_set_backend(:cuda) + device round-trip ===")
(rok, rmsg) = registry_check()
println(rok ? "  [PASS] clm_set_backend(:cuda) + clm_device_array round-trip" :
              "  [FAIL] $rmsg")
push!(results, ("backend-registry", rok))

for (i, s) in enumerate(("gpu_validate.jl", "gpu_ad_validate.jl"))
    println("\n=== [$(i+1)/3] $s on CUDA ===")
    (ok, out) = run_suite(s)
    for ln in split(out, '\n')
        occursin(r"\b(passed|PASS|FAIL|MATCH|DIVERGED|ERROR)\b|✓", ln) && println("  ", strip(ln))
    end
    push!(results, (s, ok))
end

# --- verdict ---------------------------------------------------------------
println("\n", "="^64)
for (nm, ok) in results
    @printf("  %-22s %s\n", nm, ok ? "PASS" : "FAIL")
end
npass = count(r -> r[2], results); n = length(results)
@printf("\n  %d/%d CUDA validation suites passed\n", npass, n)
if npass == n
    devname = Base.invokelatest(() -> _CUDA.name(_CUDA.device()))
    println("\n  ✓ CUDA BACKEND VALIDATED vs CPU on $devname")
else
    println("\n  ✗ some CUDA validations diverged/failed — see suite output above")
end
exit(npass == n ? 0 : 1)
