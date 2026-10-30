#!/usr/bin/env julia
# ==========================================================================
# gpu_validate_amdgpu.jl — one-command AMD ROCm (AMDGPU.jl) backend validation.
#
# The AMDGPU twin of scripts/gpu_validate_cuda.jl. CLM.jl's :amdgpu backend
# (ext/CLMAMDGPUExt.jl) has NOT yet run on AMD silicon; this is the script to
# run on a ROCm box to change that. See docs/GPU_AMDGPU.md.
#
#   julia --project=scripts -e 'using Pkg; Pkg.develop(path="."); Pkg.add("AMDGPU"); Pkg.instantiate()'
#   julia --project=scripts scripts/gpu_validate_amdgpu.jl
#
# Exit code 0 = all AMDGPU validations passed; 1 = a divergence/failure; 2 = no
# functional AMDGPU device (so CI can distinguish "no GPU" from "GPU mismatch").
#
# What it runs (each suite in a fresh process, pinned to AMDGPU via
# CLM_GPU_BACKEND=AMDGPU so scripts/gpu_backends.jl cannot pick another GPU):
#   1. AMDGPU functional check + device report.
#   2. Backend registry round-trip — clm_set_backend(:amdgpu) + clm_device_array.
#   3. test/test_amdgpu_backend.jl — single-column clm_drv! device vs CPU.
#   4. scripts/gpu_validate.jl — core kernel-parity suite.
#   5. scripts/gpu_ad_validate.jl — forward-mode AD parity on the device.
#   6. scripts/gpu_validate_clmdrv_e2e.jl — full biogeophysics/hydrology timestep.
#   7. scripts/gpu_validate_clmdrv_cn_e2e.jl — full timestep with CN enabled.
#   8. scripts/gpu_validate_clmdrv_realinit_e2e.jl — multi-step real NetCDF init
#      (requires a discoverable CLM parameter file).
# ==========================================================================

using CLM
using Printf

const _AMDGPU = try
    @eval using AMDGPU
    Base.invokelatest(getfield, @__MODULE__, :AMDGPU)
catch e
    @warn "AMDGPU.jl could not be loaded (Pkg.add(\"AMDGPU\") on a ROCm box)" exception=(e, catch_backtrace())
    nothing
end

if _AMDGPU === nothing || !Base.invokelatest(() -> _AMDGPU.functional())
    println("""
    ✗ AMDGPU is NOT functional in this environment.
      This script must run on a host with an AMD GPU + ROCm + AMDGPU.jl.
      Nothing to validate here.""")
    exit(2)
end

# --- device report ---------------------------------------------------------
Base.invokelatest() do
    @printf("  AMDGPU device     : %s\n", string(_AMDGPU.device()))
    try; @printf("  AMDGPU.jl         : %s\n", string(pkgversion(_AMDGPU))); catch; end
    println("  precision under test: Float64")
end
println()

function registry_check()
    try
        CLM.clm_set_backend(:amdgpu; require_functional = true)
        x = collect(1.0:8.0)
        d = CLM.clm_device_array(x)            # adapt(ROCArray, x)
        ok = Base.invokelatest(() -> d isa _AMDGPU.ROCArray{Float64}) && Array(d) == x
        CLM.clm_set_backend(:cpu)
        return (ok, ok ? "" : "device round-trip mismatch / wrong array type")
    catch e
        try; CLM.clm_set_backend(:cpu); catch; end
        return (false, sprint(showerror, e))
    end
end

function run_suite(path, proj)
    jl  = Base.julia_cmd()
    buf = IOBuffer()
    cmd = addenv(`$jl --project=$proj $path`, "CLM_GPU_BACKEND" => "AMDGPU")
    ok = try
        run(pipeline(cmd; stdout = buf, stderr = buf)); true
    catch
        false
    end
    return (ok, String(take!(buf)))
end

const SCRIPTS = abspath(@__DIR__)
const SUITES = (
    joinpath(SCRIPTS, "..", "test", "test_amdgpu_backend.jl"),
    joinpath(SCRIPTS, "gpu_validate.jl"),
    joinpath(SCRIPTS, "gpu_ad_validate.jl"),
    joinpath(SCRIPTS, "gpu_validate_clmdrv_e2e.jl"),
    joinpath(SCRIPTS, "gpu_validate_clmdrv_cn_e2e.jl"),
    joinpath(SCRIPTS, "gpu_validate_clmdrv_realinit_e2e.jl"),
)
const NVALIDATIONS = 1 + length(SUITES)

results = Tuple{String,Bool}[]
println("=== [1/$NVALIDATIONS] backend registry: clm_set_backend(:amdgpu) + device round-trip ===")
(rok, rmsg) = registry_check()
println(rok ? "  [PASS] clm_set_backend(:amdgpu) + clm_device_array round-trip" : "  [FAIL] $rmsg")
push!(results, ("backend-registry", rok))

for (i, s) in enumerate(SUITES)
    nm = basename(s)
    println("\n=== [$(i+1)/$NVALIDATIONS] $nm on AMDGPU ===")
    (ok, out) = run_suite(s, SCRIPTS)
    if ok
        for ln in split(out, '\n')
            occursin(r"\b(passed|PASS|FAIL|MATCH|DIVERGED|ERROR|amdgpu)\b|✓|Test Summary", ln) &&
                println("  ", strip(ln))
        end
    else
        print(out)    # keep the full kernel/IR diagnosis on failure
    end
    push!(results, (nm, ok))
end

println("\n", "="^64)
for (nm, ok) in results
    @printf("  %-40s %s\n", nm, ok ? "PASS" : "FAIL")
end
npass = count(r -> r[2], results); n = length(results)
@printf("\n  %d/%d AMDGPU validation suites passed\n", npass, n)
println(npass == n ? "\n  ✓ AMDGPU BACKEND VALIDATED vs CPU" :
                     "\n  ✗ some AMDGPU validations diverged/failed — see suite output above")
exit(npass == n ? 0 : 1)
