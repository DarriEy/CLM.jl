# ==========================================================================
# gpu_validate_bareground_e2e.jl — end-to-end Metal parity for the WHOLE
# bareground_fluxes! driver (Phase 1 init + Phase 2 Monin-Obukhov init +
# Phase 3 stability iteration (friction_velocity! + stability kernel) +
# Phase 4 post-iteration fluxes / diagnostics).
#
# Mirrors test/test_bareground_fluxes.jl's single-patch smoke setup at Float32,
# runs bareground_fluxes! on the CPU, adapts every state struct (and the forcing
# arrays the kernels touch) to Metal, runs the SAME call on the device, and
# compares the mutated outputs with reldiff.
#
# CRITICAL: every required scalar param / input is set to a real CLM default so
# the CPU reference is FINITE — reldiff silently PASSES when both sides are NaN,
# so a false PASS from unset NaN inputs is guarded against here (we assert the
# CPU reference fields are finite before trusting parity).
#
#   julia --project=scripts scripts/gpu_validate_bareground_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end

# Build the single-patch bareground instance at precision FT (mirrors the test).
function build(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    np = 1; nc = 1; ng = 1; nl = 1
    nlevgrnd = CLM.varpar.nlevgrnd
    nlevsno  = CLM.varpar.nlevsno

    canopystate = CLM.CanopyStateData{FT}(); CLM.canopystate_init!(canopystate, np)
    canopystate.displa_patch[1] = 0.0

    energyflux = CLM.EnergyFluxData{FT}(); CLM.energyflux_init!(energyflux, np, nc, nl, ng)
    energyflux.htvp_col[1] = CLM.HVAP
    energyflux.btran_patch[1] = 0.0
    energyflux.cgrnds_patch[1] = 0.0; energyflux.cgrndl_patch[1] = 0.0; energyflux.cgrnd_patch[1] = 0.0
    energyflux.dhsdt_canopy_patch[1] = 0.0; energyflux.eflx_sh_stem_patch[1] = 0.0
    for j in 1:nlevgrnd; energyflux.rresis_patch[1, j] = 0.0; end

    frictionvel = CLM.FrictionVelocityData{FT}(); CLM.frictionvel_init!(frictionvel, np, nc)
    frictionvel.zetamaxstable = 0.5; frictionvel.zsno = 0.00085; frictionvel.zlnd = 0.000775
    frictionvel.z0mg_col[1] = 0.01; frictionvel.z0hg_col[1] = 0.01; frictionvel.z0qg_col[1] = 0.01
    frictionvel.z0mg_patch[1] = 0.01; frictionvel.z0hg_patch[1] = 0.01; frictionvel.z0qg_patch[1] = 0.01
    frictionvel.z0mv_patch[1] = 0.0; frictionvel.z0hv_patch[1] = 0.0; frictionvel.z0qv_patch[1] = 0.0
    frictionvel.forc_hgt_u_patch[1] = 30.0; frictionvel.forc_hgt_t_patch[1] = 30.0
    frictionvel.forc_hgt_q_patch[1] = 30.0
    frictionvel.ustar_patch[1] = 0.5; frictionvel.um_patch[1] = 5.0
    frictionvel.obu_patch[1] = -100.0; frictionvel.zeta_patch[1] = -0.1
    frictionvel.ram1_patch[1] = 50.0; frictionvel.num_iter_patch[1] = 0.0
    frictionvel.kbm1_patch[1] = 0.0

    temperature = CLM.TemperatureData{FT}(); CLM.temperature_init!(temperature, np, nc, nl, ng)
    temperature.t_veg_patch[1] = 290.0; temperature.thm_patch[1] = 290.0
    temperature.t_grnd_col[1] = 288.0; temperature.t_h2osfc_col[1] = 288.0
    temperature.thv_col[1] = 291.0; temperature.beta_col[1] = 1.0
    temperature.t_ref2m_patch[1] = 290.0; temperature.t_ref2m_r_patch[1] = 290.0
    for j in 1:size(temperature.t_soisno_col, 2); temperature.t_soisno_col[1, j] = 288.0; end

    soilstate = CLM.SoilStateData{FT}(); CLM.soilstate_init!(soilstate, np, nc)
    soilstate.soilbeta_col[1] = 0.8; soilstate.soilresis_col[1] = 100.0
    for j in 1:nlevgrnd
        soilstate.watsat_col[1, j] = 0.45
        soilstate.rootr_patch[1, j] = 0.0
    end

    waterfluxbulk = CLM.WaterFluxBulkData{FT}(); CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)
    waterfluxbulk.wf.qflx_tran_veg_patch[1] = 0.0; waterfluxbulk.wf.qflx_evap_veg_patch[1] = 0.0
    waterfluxbulk.wf.qflx_evap_soi_patch[1] = 0.0; waterfluxbulk.wf.qflx_evap_tot_patch[1] = 0.0

    waterstatebulk = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
    waterstatebulk.ws.h2osoi_liq_col[1, 1 + nlevsno] = 50.0
    waterstatebulk.ws.h2osoi_ice_col[1, 1 + nlevsno] = 0.0

    waterdiagbulk = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
    waterdiagbulk.qg_col[1] = 0.005; waterdiagbulk.qg_snow_col[1] = 0.005
    waterdiagbulk.qg_soil_col[1] = 0.005; waterdiagbulk.qg_h2osfc_col[1] = 0.005
    waterdiagbulk.dqgdT_col[1] = 0.0003; waterdiagbulk.q_ref2m_patch[1] = 0.008
    waterdiagbulk.rh_ref2m_patch[1] = 50.0; waterdiagbulk.rh_ref2m_r_patch[1] = 50.0

    photosyns = CLM.PhotosynthesisData{FT}(); CLM.photosynthesis_data_init!(photosyns, np)
    photosyns.rssun_patch[1] = 200.0; photosyns.rssha_patch[1] = 300.0

    patch_data = CLM.PatchData{FT}(); CLM.patch_init!(patch_data, np)
    patch_data.column[1] = 1; patch_data.gridcell[1] = 1; patch_data.landunit[1] = 1
    patch_data.itype[1] = 1; patch_data.active[1] = true

    col_data = CLM.ColumnData{FT}(); CLM.column_init!(col_data, nc)
    col_data.snl[1] = 0; col_data.zii[1] = 1000.0; col_data.dz[1, 1 + nlevsno] = 0.1

    lun_data = CLM.LandunitData{FT}(); CLM.landunit_init!(lun_data, nl)
    lun_data.itype[1] = CLM.ISTSOIL

    S = (; canopystate, energyflux, frictionvel, temperature, soilstate,
           waterfluxbulk, waterstatebulk, waterdiagbulk, photosyns, patch_data,
           col_data, lun_data)
    mask = falses(np); mask[1] = true
    forc = (; q = FT[0.008], pbot = FT[101325.0], th = FT[290.0], rho = FT[1.2],
              t = FT[288.0], u_grc = FT[3.0], v_grc = FT[1.0],
              hgt_t = FT[10.0], hgt_u = FT[10.0], hgt_q = FT[10.0],
              # CH4 conductance output: a per-patch FT array so it stays device-
              # resident under adapt (use_lch4=false, so values are unused but the
              # kernel still needs an isbits device buffer to compile on Metal).
              grnd_ch4 = zeros(FT, np))
    return (; np, nc, S, mask, forc)
end

run_bgf!(S, mask, forc) = CLM.bareground_fluxes!(
    S.canopystate, S.energyflux, S.frictionvel, S.temperature, S.soilstate,
    S.waterfluxbulk, S.waterstatebulk, S.waterdiagbulk, S.photosyns,
    S.patch_data, S.col_data, S.lun_data, mask, 1:1,
    forc.q, forc.pbot, forc.th, forc.rho, forc.t,
    forc.u_grc, forc.v_grc, forc.hgt_t, forc.hgt_u, forc.hgt_q;
    use_lch4 = false, z0param_method = "ZengWang2007",
    grnd_ch4_cond_patch = forc.grnd_ch4)

function setconfig!()
    CLM.soil_resistance_read_nl!(soil_resis_method = CLM.SOIL_RESIS_LEEPIELKE_1992)
    CLM.bareground_fluxes_read_params!(a_coef = 0.5, a_exp = 1.0, wind_min = 1.0)
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for bareground_fluxes! (whole driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    setconfig!()
    H = build(FT)   # CPU reference
    setconfig!()
    B = build(FT)   # device copy

    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    Sd = map(ad, B.S)
    forc_d = (; (k => ad(getfield(B.forc, k)) for k in keys(B.forc))...)
    mask_d = B.mask   # host-only mask (kept on host by design)

    if !(Sd.temperature.t_veg_patch isa Metal.MtlArray)
        println("  BLOCKED: a bareground state struct did not move to the device under adapt.")
        return 2
    end

    run_bgf!(H.S, H.mask, H.forc)     # CPU reference
    run_bgf!(Sd, mask_d, forc_d)      # device

    # Guard against a false PASS: the CPU reference fields must be FINITE.
    refchecks = [
        ("eflx_sh_grnd", H.S.energyflux.eflx_sh_grnd_patch),
        ("taux",         H.S.energyflux.taux_patch),
        ("cgrnd",        H.S.energyflux.cgrnd_patch),
        ("t_ref2m",      H.S.temperature.t_ref2m_patch),
        ("rh_ref2m",     H.S.waterdiagbulk.rh_ref2m_patch),
        ("ustar",        H.S.frictionvel.ustar_patch),
        ("qflx_evap_soi",H.S.waterfluxbulk.wf.qflx_evap_soi_patch),
    ]
    for (nm, a) in refchecks
        if !all(isfinite, Array(a))
            @printf("  BLOCKED: CPU reference field %s is NOT finite — parity would be a false PASS.\n", nm)
            return 2
        end
    end

    checks = [
        ("eflx_sh_grnd",  H.S.energyflux.eflx_sh_grnd_patch,  Sd.energyflux.eflx_sh_grnd_patch),
        ("eflx_sh_tot",   H.S.energyflux.eflx_sh_tot_patch,   Sd.energyflux.eflx_sh_tot_patch),
        ("eflx_sh_snow",  H.S.energyflux.eflx_sh_snow_patch,  Sd.energyflux.eflx_sh_snow_patch),
        ("eflx_sh_soil",  H.S.energyflux.eflx_sh_soil_patch,  Sd.energyflux.eflx_sh_soil_patch),
        ("eflx_sh_h2osfc",H.S.energyflux.eflx_sh_h2osfc_patch,Sd.energyflux.eflx_sh_h2osfc_patch),
        ("taux",          H.S.energyflux.taux_patch,          Sd.energyflux.taux_patch),
        ("tauy",          H.S.energyflux.tauy_patch,          Sd.energyflux.tauy_patch),
        ("cgrnds",        H.S.energyflux.cgrnds_patch,        Sd.energyflux.cgrnds_patch),
        ("cgrndl",        H.S.energyflux.cgrndl_patch,        Sd.energyflux.cgrndl_patch),
        ("cgrnd",         H.S.energyflux.cgrnd_patch,         Sd.energyflux.cgrnd_patch),
        ("ustar",         H.S.frictionvel.ustar_patch,        Sd.frictionvel.ustar_patch),
        ("ram1",          H.S.frictionvel.ram1_patch,         Sd.frictionvel.ram1_patch),
        ("kbm1",          H.S.frictionvel.kbm1_patch,         Sd.frictionvel.kbm1_patch),
        ("num_iter",      H.S.frictionvel.num_iter_patch,     Sd.frictionvel.num_iter_patch),
        ("z0hg_col",      H.S.frictionvel.z0hg_col,           Sd.frictionvel.z0hg_col),
        ("t_ref2m",       H.S.temperature.t_ref2m_patch,      Sd.temperature.t_ref2m_patch),
        ("t_ref2m_r",     H.S.temperature.t_ref2m_r_patch,    Sd.temperature.t_ref2m_r_patch),
        ("t_veg",         H.S.temperature.t_veg_patch,        Sd.temperature.t_veg_patch),
        ("q_ref2m",       H.S.waterdiagbulk.q_ref2m_patch,    Sd.waterdiagbulk.q_ref2m_patch),
        ("rh_ref2m",      H.S.waterdiagbulk.rh_ref2m_patch,   Sd.waterdiagbulk.rh_ref2m_patch),
        ("rh_ref2m_r",    H.S.waterdiagbulk.rh_ref2m_r_patch, Sd.waterdiagbulk.rh_ref2m_r_patch),
        ("qflx_evap_soi", H.S.waterfluxbulk.wf.qflx_evap_soi_patch, Sd.waterfluxbulk.wf.qflx_evap_soi_patch),
        ("qflx_evap_tot", H.S.waterfluxbulk.wf.qflx_evap_tot_patch, Sd.waterfluxbulk.wf.qflx_evap_tot_patch),
        ("qflx_ev_snow",  H.S.waterfluxbulk.qflx_ev_snow_patch, Sd.waterfluxbulk.qflx_ev_snow_patch),
        ("rssun",         H.S.photosyns.rssun_patch,          Sd.photosyns.rssun_patch),
    ]
    nfail = 0; worst = 0.0
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-2
        worst = max(worst, d)
        @printf("  [%s] %-15s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    @printf("\n  worst reldiff = %.3e\n", worst)
    println(nfail == 0 ? "  WHOLE bareground_fluxes! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
