# ==========================================================================
# gpu_validate_canopy_e2e.jl — end-to-end Metal parity for the WHOLE
# canopy_fluxes! driver (Newton iteration + all surrounding kernels).
#
# Mirrors test/test_canopy_fluxes.jl's single-patch smoke setup at Float32,
# runs canopy_fluxes! on the CPU, adapts every state struct (and the forcing /
# PFT-parameter arrays the kernels touch) to Metal, runs the SAME call on the
# device, and compares the mutated outputs. The full photosynthesis! solver is
# NOT exercised (nrad/parsun left empty, as in the test), so the inner Ci
# solver is out of scope here — this validates the canopy energy-balance /
# Monin-Obukhov Newton path and its surrounding loops end to end.
#
#   julia --project=scripts scripts/gpu_validate_canopy_e2e.jl
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

# Build the single-patch canopy instance at precision FT (mirrors the test).
function build(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    np = 1; nc = 1; ng = 1; nl = 1

    canopystate = CLM.CanopyStateData{FT}(); CLM.canopystate_init!(canopystate, np)
    canopystate.elai_patch[1] = 2.0;  canopystate.esai_patch[1] = 0.5
    canopystate.laisun_patch[1] = 1.2; canopystate.laisha_patch[1] = 0.8
    canopystate.displa_patch[1] = 5.0; canopystate.htop_patch[1] = 10.0
    canopystate.frac_veg_nosno_patch[1] = 1; canopystate.dleaf_patch[1] = 0.04
    canopystate.stem_biomass_patch[1] = 0.0; canopystate.leaf_biomass_patch[1] = 0.0

    energyflux = CLM.EnergyFluxData{FT}(); CLM.energyflux_init!(energyflux, np, nc, nl, ng)
    energyflux.btran_patch[1] = 0.5; energyflux.bsun_patch[1] = 0.5; energyflux.bsha_patch[1] = 0.5
    energyflux.htvp_col[1] = CLM.HVAP
    energyflux.cgrnds_patch[1] = 0.0; energyflux.cgrndl_patch[1] = 0.0; energyflux.cgrnd_patch[1] = 0.0
    energyflux.dhsdt_canopy_patch[1] = 0.0; energyflux.eflx_sh_stem_patch[1] = 0.0

    frictionvel = CLM.FrictionVelocityData{FT}(); CLM.frictionvel_init!(frictionvel, np, nc)
    frictionvel.zetamaxstable = 0.5; frictionvel.zsno = 0.00085; frictionvel.zlnd = 0.000775
    frictionvel.z0mv_patch[1] = 0.5; frictionvel.z0hv_patch[1] = 0.5; frictionvel.z0qv_patch[1] = 0.5
    frictionvel.z0mg_col[1] = 0.01; frictionvel.z0hg_col[1] = 0.01; frictionvel.z0qg_col[1] = 0.01
    frictionvel.forc_hgt_u_patch[1] = 30.0; frictionvel.forc_hgt_t_patch[1] = 30.0
    frictionvel.forc_hgt_q_patch[1] = 30.0
    frictionvel.ustar_patch[1] = 0.5; frictionvel.um_patch[1] = 5.0; frictionvel.uaf_patch[1] = 3.0
    frictionvel.taf_patch[1] = 290.0; frictionvel.qaf_patch[1] = 0.008
    frictionvel.obu_patch[1] = -100.0; frictionvel.zeta_patch[1] = -0.1
    frictionvel.vpd_patch[1] = 1.0; frictionvel.rb1_patch[1] = 50.0
    frictionvel.ram1_patch[1] = 50.0; frictionvel.num_iter_patch[1] = 0.0

    temperature = CLM.TemperatureData{FT}(); CLM.temperature_init!(temperature, np, nc, nl, ng)
    temperature.t_veg_patch[1] = 290.0; temperature.t_stem_patch[1] = 289.0
    temperature.t_skin_patch[1] = 290.0; temperature.thm_patch[1] = 290.0
    temperature.t_grnd_col[1] = 288.0; temperature.t_h2osfc_col[1] = 288.0
    temperature.thv_col[1] = 291.0; temperature.emv_patch[1] = 0.97; temperature.emg_col[1] = 0.96
    temperature.t_ref2m_patch[1] = 290.0; temperature.t_ref2m_r_patch[1] = 290.0
    for j in 1:size(temperature.t_soisno_col, 2); temperature.t_soisno_col[1, j] = 288.0; end

    solarabs = CLM.SolarAbsorbedData{FT}(); CLM.solarabs_init!(solarabs, np, nl)
    solarabs.sabv_patch[1] = 150.0

    soilstate = CLM.SoilStateData{FT}(); CLM.soilstate_init!(soilstate, nc, nl)
    soilstate.soilbeta_col[1] = 0.8; soilstate.soilresis_col[1] = 100.0
    for j in 1:CLM.varpar.nlevgrnd; soilstate.watsat_col[1, j] = 0.45; end

    waterfluxbulk = CLM.WaterFluxBulkData{FT}(); CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)
    waterfluxbulk.wf.qflx_tran_veg_patch[1] = 0.0; waterfluxbulk.wf.qflx_evap_veg_patch[1] = 0.0
    waterfluxbulk.wf.qflx_evap_soi_patch[1] = 0.0

    waterstatebulk = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
    waterstatebulk.ws.liqcan_patch[1] = 0.1; waterstatebulk.ws.snocan_patch[1] = 0.0

    waterdiagbulk = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
    waterdiagbulk.fwet_patch[1] = 0.1; waterdiagbulk.fdry_patch[1] = 0.8
    waterdiagbulk.frac_sno_eff_col[1] = 0.0; waterdiagbulk.frac_h2osfc_col[1] = 0.0
    waterdiagbulk.snow_depth_col[1] = 0.0; waterdiagbulk.qg_col[1] = 0.005
    waterdiagbulk.qg_snow_col[1] = 0.005; waterdiagbulk.qg_soil_col[1] = 0.005
    waterdiagbulk.qg_h2osfc_col[1] = 0.005; waterdiagbulk.dqgdT_col[1] = 0.0003
    waterdiagbulk.rh_af_patch[1] = 0.6

    photosyns = CLM.PhotosynthesisData{FT}(); CLM.photosynthesis_data_init!(photosyns, np)
    photosyns.rssun_patch[1] = 200.0; photosyns.rssha_patch[1] = 300.0

    patch_data = CLM.PatchData{FT}(); CLM.patch_init!(patch_data, np)
    patch_data.column[1] = 1; patch_data.gridcell[1] = 1; patch_data.landunit[1] = 1
    patch_data.itype[1] = 1; patch_data.active[1] = true

    col_data = CLM.ColumnData{FT}(); CLM.column_init!(col_data, nc); col_data.snl[1] = 0
    gridcell_data = CLM.GridcellData{FT}(); CLM.gridcell_init!(gridcell_data, ng)

    S = (; canopystate, energyflux, frictionvel, temperature, solarabs, soilstate,
           waterfluxbulk, waterstatebulk, waterdiagbulk, photosyns, patch_data,
           col_data, gridcell_data)
    mask = falses(np); mask[1] = true
    forc = (; lwrad = FT[350.0], q = FT[0.008], pbot = FT[101325.0], th = FT[290.0],
              rho = FT[1.2], t = FT[288.0], u_grc = FT[3.0], v_grc = FT[1.0],
              pco2 = FT[40.0], po2 = FT[21000.0], hgt_t = FT[10.0], hgt_u = FT[10.0],
              hgt_q = FT[10.0], dayl = FT[43200.0], max_dayl = FT[50000.0],
              downreg = FT[1.0], leafn = FT[1.0])
    mp = CLM.MXPFT + 1
    pft = (; dleaf = fill(FT(0.04), mp), z0v_Cr = fill(FT(0.35), mp),
             z0v_Cs = fill(FT(0.003), mp), z0v_c = fill(FT(0.25), mp),
             z0v_cw = fill(FT(2.0), mp), z0v_LAImax = fill(FT(8.0), mp),
             grnd_ch4 = fill(FT(0.0), np))
    return (; np, nc, S, mask, forc, pft, dtime = FT(1800.0))
end

run_cf!(S, mask, forc, pft, dtime) = CLM.canopy_fluxes!(
    S.canopystate, S.energyflux, S.frictionvel, S.temperature, S.solarabs,
    S.soilstate, S.waterfluxbulk, S.waterstatebulk, S.waterdiagbulk, S.photosyns,
    S.patch_data, S.col_data, S.gridcell_data, mask, 1:1, 1:1,
    forc.lwrad, forc.q, forc.pbot, forc.th, forc.rho, forc.t,
    forc.u_grc, forc.v_grc, forc.pco2, forc.po2, forc.hgt_t, forc.hgt_u, forc.hgt_q,
    forc.dayl, forc.max_dayl, forc.downreg, forc.leafn, dtime;
    dleaf_pft = pft.dleaf, z0v_Cr_pft = pft.z0v_Cr, z0v_Cs_pft = pft.z0v_Cs,
    z0v_c_pft = pft.z0v_c, z0v_cw_pft = pft.z0v_cw, z0v_LAImax_pft = pft.z0v_LAImax,
    grnd_ch4_cond_patch = pft.grnd_ch4)

function setconfig!()
    CLM.soil_resistance_read_nl!(soil_resis_method = CLM.SOIL_RESIS_LEEPIELKE_1992)
    CLM.canopy_fluxes_read_nml!(use_undercanopy_stability = false,
        use_biomass_heat_storage = false, itmax_canopy_fluxes = 40)
    CLM.canopy_fluxes_read_params!()
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for canopy_fluxes! (whole driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    setconfig!()
    B = build(FT)
    setconfig!()
    H = build(FT)   # independent host copy (CPU reference)

    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    Sd = map(ad, B.S)                       # adapt every state struct
    forc_d = (; (k => ad(getfield(B.forc, k)) for k in keys(B.forc))...)
    pft_d  = (; (k => ad(getfield(B.pft,  k)) for k in keys(B.pft))...)
    mask_d = B.mask   # host-only mask (kept on host by design)

    if !(Sd.temperature.t_veg_patch isa Metal.MtlArray)
        println("  BLOCKED: a canopy state struct did not move to the device under adapt.")
        return 2
    end

    run_cf!(H.S, H.mask, H.forc, H.pft, H.dtime)   # CPU reference
    run_cf!(Sd, mask_d, forc_d, pft_d, B.dtime)    # device

    checks = [
        ("t_veg",        H.S.temperature.t_veg_patch,       Sd.temperature.t_veg_patch),
        ("t_ref2m",      H.S.temperature.t_ref2m_patch,     Sd.temperature.t_ref2m_patch),
        ("t_skin",       H.S.temperature.t_skin_patch,      Sd.temperature.t_skin_patch),
        ("eflx_sh_veg",  H.S.energyflux.eflx_sh_veg_patch,  Sd.energyflux.eflx_sh_veg_patch),
        ("eflx_sh_grnd", H.S.energyflux.eflx_sh_grnd_patch, Sd.energyflux.eflx_sh_grnd_patch),
        ("taux",         H.S.energyflux.taux_patch,         Sd.energyflux.taux_patch),
        ("tauy",         H.S.energyflux.tauy_patch,         Sd.energyflux.tauy_patch),
        ("dlrad",        H.S.energyflux.dlrad_patch,        Sd.energyflux.dlrad_patch),
        ("ulrad",        H.S.energyflux.ulrad_patch,        Sd.energyflux.ulrad_patch),
        ("cgrnd",        H.S.energyflux.cgrnd_patch,        Sd.energyflux.cgrnd_patch),
        ("ustar",        H.S.frictionvel.ustar_patch,       Sd.frictionvel.ustar_patch),
        ("obu",          H.S.frictionvel.obu_patch,         Sd.frictionvel.obu_patch),
        ("um",           H.S.frictionvel.um_patch,          Sd.frictionvel.um_patch),
        ("qflx_evap_veg",H.S.waterfluxbulk.wf.qflx_evap_veg_patch, Sd.waterfluxbulk.wf.qflx_evap_veg_patch),
        ("liqcan",       H.S.waterstatebulk.ws.liqcan_patch, Sd.waterstatebulk.ws.liqcan_patch),
        ("num_iter",     H.S.frictionvel.num_iter_patch,    Sd.frictionvel.num_iter_patch),
    ]
    nfail = 0
    for (nm, a, b) in checks
        d = reldiff(a, b); ok = d < 1f-2
        @printf("  [%s] %-15s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  WHOLE canopy_fluxes! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
