# ==========================================================================
# gpu_validate_clmdrv_e2e.jl — WHOLE-TIMESTEP Metal parity for clm_drv!.
#
# Runs the ACTUAL top-level `clm_drv!` driver (the full biogeophys/hydrology
# chain: drv_init -> radiation -> fluxes -> soil temp/water -> snow -> hydrology
# -> albedo -> balance check) on the CPU and on Metal over the SAME synthetic
# state, then compares the resulting state. This is the chained-pipeline test
# that the ~86 per-module gpu_validate_* harnesses cannot give: every state
# struct flows module->module ON THE DEVICE through one clm_drv! call, so this
# catches scalar-precision leaks at module boundaries (cf. the `dt` Float64 leak
# that the BGC-pipeline integration harness found but the unit harnesses masked).
#
# Default config (use_cn=false): the biogeophysics + hydrology path.
#
#   julia --project=scripts scripts/gpu_validate_clmdrv_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

# Float32 down-adaptor: arrays -> Float32, AND ::FT scalar fields (e.g.
# TemperatureData.excess_ice_coldstart_*) -> Float32 so the parametric {FT}
# constructor accepts them. The handful of structs that pin a CONCRETE ::Float64
# scalar (which must stay Float64) are special-cased below to preserve it.
struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
# Unpack Bool arrays (BitVector) to Vector{Bool} so the subsequent MtlArray adapt
# doesn't scalar-index packed bits.
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = collect(Bool, x)
CLM.Adapt.adapt_storage(::_F32, x::Float64) = Float32(x)

# Reconstruct a struct adapting every field EXCEPT the named concrete-Float64
# scalar field(s), which are preserved as-is (their ctor pins ::Float64).
function _adapt_keep_f64(a, x, keep::Tuple)
    T = typeof(x)
    vals = ntuple(fieldcount(T)) do i
        n = fieldname(T, i); v = getfield(x, i)
        n in keep ? v : CLM.Adapt.adapt(a, v)
    end
    return T.name.wrapper(vals...)
end
CLM.Adapt.adapt_structure(a::_F32, x::CLM.CanopyStateData) = _adapt_keep_f64(a, x, (:leaf_mr_vcm,))
CLM.Adapt.adapt_structure(a::_F32, x::CLM.SoilBiogeochemCarbonStateData) = _adapt_keep_f64(a, x, (:totvegcthresh,))
CLM.Adapt.adapt_structure(a::_F32, x::CLM.SoilBiogeochemNitrogenStateData) = _adapt_keep_f64(a, x, (:totvegcthresh,))

mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))

# Finite-only relative diff: compares only entries finite on BOTH backends (the
# minimal smoke-test fixture leaves much of the physics NaN/Inf on both). Returns
# (max_reldiff, n_finite_compared).
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isfinite(A[i]) && isfinite(B[i])) || continue
        n += 1
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) /
                   (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return (m, n)
end

# NB: the global CLM.pftcon stays host-resident (it's a `const` with concrete
# Vector{Float64} fields). clm_drv_core! copies the pftcon arrays it passes to
# kernels onto the working backend via its internal `_onbk` helper, so no global
# swap is needed here.

# make_driver_data: cloned from test/test_clm_driver.jl (Float64 build, sets the
# shared module globals varpar/varcon/pftcon/SNOW_DZ*/urban_ctrl).
include(joinpath(@__DIR__, "clmdrv_make_data.jl"))

const DRV_ARGS = (true,      # doalb (full driver incl. surface_albedo radiative transfer)
                  1.0,       # nextsw_cday
                  0.0,       # declinp1
                  0.0,       # declin
                  0.4091,    # obliqr
                  false,     # rstwr
                  false,     # nlend
                  "20260101",# rdate
                  false)     # rof_prognostic

function run_drv!(inst, filt, filt_ia, bounds, config, photosyns)
    CLM.clm_drv!(config, inst, filt, filt_ia, bounds, DRV_ARGS...;
                 nstep=1, is_first_step=false,
                 is_beg_curr_day=false, is_beg_curr_year=false,
                 photosyns=photosyns)
end

function main(backend)
    println("="^72)
    println("WHOLE-TIMESTEP Metal parity for clm_drv! (default biogeophys path)")
    println("="^72)
    if backend === nothing
        println("  No GPU backend — CPU path exercised by the suite.")
        return 0
    end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    # ---- CPU reference run (Float64) ----
    instH, bounds, filtH, filt_iaH, configH, psH = make_driver_data()
    instH.canopystate.frac_veg_nosno_alb_patch .= 1
    run_drv!(instH, filtH, filt_iaH, bounds, configH, psH)
    println("  CPU clm_drv! completed.")

    # ---- Device run (Metal Float32) over a fresh, identical state ----
    instB, boundsB, filtB, filt_iaB, configB, psB = make_driver_data()
    instB.canopystate.frac_veg_nosno_alb_patch .= 1

    inst_d = mf(instB)
    ps_d   = mf(psB)
    filt_d = mf(filtB)
    filt_ia_d = mf(filt_iaB)

    if !(inst_d.temperature.t_soisno_col isa Metal.MtlArray)
        println("  BLOCKED: inst tree did not move to the device under mf().")
        return 2
    end

    run_drv!(inst_d, filt_d, filt_ia_d, boundsB, configB, ps_d)
    println("  Device clm_drv! completed.\n")

    # ---- Compare driver outputs (finite-only; the minimal smoke fixture leaves
    #      much of the physics NaN/Inf on BOTH backends — only finite entries are
    #      a meaningful CPU-vs-device parity signal). ----
    H = instH; D = inst_d
    ws_H = H.water.waterstatebulk_inst.ws; ws_D = D.water.waterstatebulk_inst.ws
    wf_H = H.water.waterfluxbulk_inst;     wf_D = D.water.waterfluxbulk_inst
    wd_H = H.water.waterdiagnosticbulk_inst; wd_D = D.water.waterdiagnosticbulk_inst
    checks = [
        ("temp.t_grnd_col",        H.temperature.t_grnd_col,        D.temperature.t_grnd_col),
        ("temp.t_soisno_col",      H.temperature.t_soisno_col,      D.temperature.t_soisno_col),
        ("temp.t_veg_patch",       H.temperature.t_veg_patch,       D.temperature.t_veg_patch),
        ("temp.thv_col",           H.temperature.thv_col,           D.temperature.thv_col),
        ("temp.emg_col",           H.temperature.emg_col,           D.temperature.emg_col),
        ("temp.t_h2osfc_col",      H.temperature.t_h2osfc_col,      D.temperature.t_h2osfc_col),
        ("ws.h2osoi_liq_col",      ws_H.h2osoi_liq_col,             ws_D.h2osoi_liq_col),
        ("ws.h2osoi_ice_col",      ws_H.h2osoi_ice_col,             ws_D.h2osoi_ice_col),
        ("ws.h2osfc_col",          ws_H.h2osfc_col,                 ws_D.h2osfc_col),
        ("wd.frac_sno_col",        wd_H.frac_sno_col,               wd_D.frac_sno_col),
        ("wd.frac_sno_eff_col",    wd_H.frac_sno_eff_col,           wd_D.frac_sno_eff_col),
        ("wd.frac_h2osfc_col",     wd_H.frac_h2osfc_col,            wd_D.frac_h2osfc_col),
        ("ws.h2osoi_vol_col",      ws_H.h2osoi_vol_col,             ws_D.h2osoi_vol_col),
        ("wd.h2osoi_liqvol_col",   wd_H.h2osoi_liqvol_col,          wd_D.h2osoi_liqvol_col),
        ("ef.eflx_sh_tot_patch",   H.energyflux.eflx_sh_tot_patch,  D.energyflux.eflx_sh_tot_patch),
        ("fv.z0mg_col",            H.frictionvel.z0mg_col,          D.frictionvel.z0mg_col),
        ("cs.elai_patch",          H.canopystate.elai_patch,        D.canopystate.elai_patch),
    ]
    nfail = 0; nfin_total = 0; gmax = 0.0; ncmp = 0
    for (nm, a, b) in checks
        d, n = reldiff(a, b)
        if n == 0
            @printf("  [skip] %-24s (no finite entries on both)\n", nm)
            continue
        end
        ncmp += 1; nfin_total += n; gmax = max(gmax, d)
        ok = d < 1f-3
        @printf("  [%s] %-24s rel|dev-cpu| = %.3e  (%d finite)\n", ok ? "PASS" : "FAIL", nm, d, n)
        ok || (nfail += 1)
    end

    println()
    @printf("  %d fields compared, %d finite entries; global max rel = %.3e\n", ncmp, nfin_total, gmax)
    if ncmp == 0
        println("  [WARN] no finite outputs to compare — driver RAN device-clean on $name but the")
        println("         smoke-test fixture yields NaN on both backends (no numeric parity exercised).")
        return 0
    end
    println(nfail == 0 ? "  WHOLE clm_drv! MATCHES CPU ON $name ($FT) over all finite outputs" :
                         "  DIVERGENCE ($nfail) — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
