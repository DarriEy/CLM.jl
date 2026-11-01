# ==========================================================================
# test_amdgpu_backend.jl — the :amdgpu (AMD ROCm) backend extension.
#
# AMDGPU is a WEAKDEP: nothing here may make it a hard dependency, and the
# default suite must stay green on a machine that has never heard of AMDGPU.
# The file therefore has three tiers, each gated on what is actually present:
#
#   1. AMDGPU.jl not installed in the active environment stack
#        -> everything below is @test_skip (the default on CI and laptops).
#   2. AMDGPU.jl installed but no functional ROCm device (AMDGPU.functional()
#      is false — e.g. a CUDA or Metal box with AMDGPU.jl added)
#        -> the CLMAMDGPUExt extension must load and register :amdgpu; the
#           registry hooks resolve to ROCArray / ROCBackend; the device run is
#           @test_skip.
#   3. A functional AMD GPU
#        -> one short whole-timestep `clm_drv!` on the device (Float64 ROCArray
#           state) compared against the same step on the CPU. It uses the
#           small driver fixture from the GPU harnesses
#           (ng=2, nl=3, nc=4, np=6), not ng=nl=nc=np=1: on CPU the 1x1x1x1
#           fixture leaves t_grnd / t_veg / eflx_sh_tot all NaN, so a
#           single-column comparison would only check the soil-water fields.
#           The CUDA whole-driver gate uses this same 4-column fixture.
#
# Tier 3 has NOT been executed on AMD hardware as of this writing (no ROCm
# device was reachable); see docs/GPU_AMDGPU.md for what is and is not verified.
# ==========================================================================

using Test
using CLM
using Printf
import KernelAbstractions as KA

const _AMDGPU_INSTALLED = Base.find_package("AMDGPU") !== nothing

# Load AMDGPU at top level (its own statement) so the extension's __init__ has
# run and the new world age is visible to everything below.
_AMDGPU_INSTALLED && @eval import AMDGPU

const _AMDGPU_FUNCTIONAL = _AMDGPU_INSTALLED && Base.invokelatest(() -> AMDGPU.functional())

# Small driver fixture shared with the GPU harnesses (scripts/clmdrv_make_data.jl).
# Only included when a device run will actually happen.
if _AMDGPU_FUNCTIONAL
    include(joinpath(@__DIR__, "..", "scripts", "clmdrv_make_data.jl"))
end

# Device arrays reject BitArray; masks become Vector{Bool} before adaptation
# (same host-prep step as scripts/gpu_adapt.jl, Float64 path).
struct _AMDHostPrep end
CLM.Adapt.adapt_storage(::_AMDHostPrep, x::AbstractArray{Bool}) = collect(Bool, x)
CLM.Adapt.adapt_storage(::_AMDHostPrep, x::AbstractArray) = x
_to_amd(x) = CLM.clm_device_array(CLM.Adapt.adapt(_AMDHostPrep(), x); backend = :amdgpu)

function _finite_reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        n += 1
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) /
                   (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return (m, n)
end

const _AMD_DRV_ARGS = (true, 1.0, 0.0, 0.0, 0.4091, false, false, "20260101", false)

_amd_run_drv!(inst, filt, filt_ia, bounds, config, ps) =
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds, _AMD_DRV_ARGS...;
                 nstep = 1, is_first_step = false, is_beg_curr_day = false,
                 is_beg_curr_year = false, photosyns = ps)

@testset "AMDGPU (ROCm) backend extension" begin
    if !_AMDGPU_INSTALLED
        @test_skip "AMDGPU.jl not installed in this environment — :amdgpu extension not exercised"
    else
        ext = Base.get_extension(CLM, :CLMAMDGPUExt)
        @test ext !== nothing                                   # weakdep → extension loaded
        # The device-side `_smooth_f64() = false` overlay must actually be installed
        # (an earlier guard looked for the macro in the wrong module and skipped it).
        @test ext !== nothing && ext._SMOOTH_OVERRIDE_INSTALLED
        @test CLM.clm_backend_registered(:amdgpu)
        @test CLM.clm_backend_functional(:amdgpu) == _AMDGPU_FUNCTIONAL
        # Hooks resolve to the AMDGPU types without touching a device.
        @test CLM.clm_ka_backend(backend = :amdgpu) isa AMDGPU.ROCBackend
        @test CLM.clm_float_type(backend = :amdgpu) == Float64 || haskey(ENV, "CLM_FLOAT")
        # Selection without a device is allowed (plumbing); requiring one is not.
        @test CLM.clm_set_backend(:amdgpu) == :amdgpu
        CLM.clm_set_backend(:cpu)
        if !_AMDGPU_FUNCTIONAL
            @test_throws ErrorException CLM.clm_set_backend(:amdgpu; require_functional = true)
            @test CLM.clm_backend() == :cpu
        end
    end
end

@testset "AMDGPU one-step clm_drv! vs CPU (skips without a ROCm device)" begin
    if !_AMDGPU_FUNCTIONAL
        @test_skip "no functional AMD GPU (AMDGPU.functional() == false or AMDGPU.jl absent)"
    else
        @test CLM.clm_set_backend(:amdgpu; require_functional = true) == :amdgpu
        try
            # CPU reference (Float64).
            instH, bounds, filtH, filt_iaH, configH, psH =
                make_driver_data_physical()
            instH.canopystate.frac_veg_nosno_alb_patch .= 1
            _amd_run_drv!(instH, filtH, filt_iaH, bounds, configH, psH)

            # Identical fresh state on the device (Float64 ROCArray, no demotion).
            instB, boundsB, filtB, filt_iaB, configB, psB =
                make_driver_data_physical()
            instB.canopystate.frac_veg_nosno_alb_patch .= 1
            inst_d = _to_amd(instB)
            @test inst_d.temperature.t_soisno_col isa AMDGPU.ROCArray
            @test eltype(inst_d.temperature.t_soisno_col) == Float64
            @test !(KA.get_backend(inst_d.temperature.t_soisno_col) isa KA.CPU)
            _amd_run_drv!(inst_d, _to_amd(filtB), _to_amd(filt_iaB), boundsB, configB,
                          _to_amd(psB))
            AMDGPU.synchronize()

            H = instH; D = inst_d
            wsH = H.water.waterstatebulk_inst.ws; wsD = D.water.waterstatebulk_inst.ws
            checks = (
                ("t_grnd_col",    H.temperature.t_grnd_col,         D.temperature.t_grnd_col),
                ("t_soisno_col",  H.temperature.t_soisno_col,       D.temperature.t_soisno_col),
                ("t_veg_patch",   H.temperature.t_veg_patch,        D.temperature.t_veg_patch),
                ("h2osoi_liq",    wsH.h2osoi_liq_col,               wsD.h2osoi_liq_col),
                ("h2osoi_ice",    wsH.h2osoi_ice_col,               wsD.h2osoi_ice_col),
                ("eflx_sh_tot",   H.energyflux.eflx_sh_tot_patch,   D.energyflux.eflx_sh_tot_patch),
                ("z0mg_col",      H.frictionvel.z0mg_col,           D.frictionvel.z0mg_col),
            )
            nfinite = 0
            for (nm, a, b) in checks
                d, n = _finite_reldiff(a, b)
                nfinite += n
                n == 0 && continue
                @printf("  amdgpu %-14s rel|dev-cpu| = %.3e (%d finite)\n", nm, d, n)
                @test d < 1e-6        # Float64 on both sides (CUDA reached ≤5e-8 here)
            end
            @test nfinite > 0         # a NaN-only fixture would prove nothing
        finally
            CLM.clm_set_backend(:cpu)
        end
    end
end
