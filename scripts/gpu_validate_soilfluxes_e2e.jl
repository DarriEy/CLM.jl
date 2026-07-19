# ==========================================================================
# gpu_validate_soilfluxes_e2e.jl — end-to-end GPU parity for the WHOLE
# soil_fluxes! driver (every per-(c)/per-(p) loop + the errsoi p2c scatter).
#
# Builds a small Float32 instance (one non-urban soil column / patch) mirroring
# test/test_soil_fluxes.jl's smoke setup, runs soil_fluxes! on the CPU, adapts
# every state struct (+ masks/forcing) to the GPU, runs the SAME call on the device,
# and compares the mutated outputs field-by-field. This exercises the flux
# correction, evap partition + limiting, ground-heat/total-flux/lwrad loops, the
# errsoi layer summation, and the patch->column errsoi scatter together.
#
#   julia --project=scripts scripts/gpu_validate_soilfluxes_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

# reldiff: NaN-aware (both-NaN agrees; one-sided NaN flags divergence). Also asserts
# the CPU reference is FINITE so a both-NaN false PASS can't slip through.
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(A[i] - B[i]) / (1.0 + max(abs(A[i]), abs(B[i]))))
    end
    return m
end
cpu_has_finite(a) = any(isfinite, Array(a))

function build(::Type{FT}) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    np = 1; nc = 1; ng = 1; nl = 1
    nlevsno  = CLM.varpar.nlevsno
    nlev_soisno = nlevsno + CLM.varpar.nlevmaxurbgrnd

    energyflux = CLM.EnergyFluxData{FT}(); CLM.energyflux_init!(energyflux, np, nc, nl, ng)
    temperature = CLM.TemperatureData{FT}(); CLM.temperature_init!(temperature, np, nc, nl, ng)
    canopystate = CLM.CanopyStateData{FT}(); CLM.canopystate_init!(canopystate, np)
    waterstatebulk = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
    waterdiagbulk = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
    waterfluxbulk = CLM.WaterFluxBulkData{FT}(); CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)
    solarabs = CLM.SolarAbsorbedData{FT}(); CLM.solarabs_init!(solarabs, np, nl)

    patch_data = CLM.PatchData{FT}(); CLM.patch_init!(patch_data, np)
    patch_data.column[1] = 1; patch_data.gridcell[1] = 1; patch_data.landunit[1] = 1
    patch_data.itype[1] = 1; patch_data.active[1] = true; patch_data.wtcol[1] = 1.0

    col_data = CLM.ColumnData{FT}(); CLM.column_init!(col_data, nc)
    col_data.snl[1] = 0; col_data.itype[1] = 1; col_data.landunit[1] = 1
    col_data.patchi[1] = 1; col_data.patchf[1] = 1

    lun_data = CLM.LandunitData{FT}(); CLM.landunit_init!(lun_data, nl)
    lun_data.itype[1] = CLM.ISTSOIL; lun_data.urbpoi[1] = false

    temperature.t_grnd_col[1] = 280.0; temperature.t_h2osfc_col[1] = 280.0
    temperature.t_h2osfc_bef_col[1] = 279.5; temperature.emg_col[1] = 0.96
    temperature.xmf_col[1] = 0.0; temperature.xmf_h2osfc_col[1] = 0.0
    temperature.c_h2osfc_col[1] = 4.188e6; temperature.t_skin_patch[1] = 280.0
    for j in 1:nlev_soisno
        temperature.t_ssbef_col[1, j] = 279.5
        temperature.t_soisno_col[1, j] = 280.0
        temperature.fact_col[1, j] = 1.0e7
    end

    for j in 1:size(waterstatebulk.ws.h2osoi_liq_col, 2)
        waterstatebulk.ws.h2osoi_liq_col[1, j] = 20.0
        waterstatebulk.ws.h2osoi_ice_col[1, j] = 5.0
    end

    waterdiagbulk.frac_sno_eff_col[1] = 0.0; waterdiagbulk.frac_h2osfc_col[1] = 0.0

    energyflux.cgrnds_patch[1] = 10.0; energyflux.cgrndl_patch[1] = 0.001
    energyflux.htvp_col[1] = CLM.HVAP
    energyflux.eflx_sh_grnd_patch[1] = 50.0; energyflux.eflx_sh_veg_patch[1] = 10.0
    energyflux.eflx_sh_stem_patch[1] = 2.0
    energyflux.dlrad_patch[1] = 300.0; energyflux.ulrad_patch[1] = 350.0
    energyflux.eflx_h2osfc_to_snow_col[1] = 0.0; energyflux.eflx_building_heat_errsoi_col[1] = 0.0
    energyflux.eflx_lwrad_net_patch[1] = 0.0; energyflux.eflx_lwrad_out_patch[1] = 0.0
    energyflux.eflx_wasteheat_patch[1] = 0.0; energyflux.eflx_heat_from_ac_patch[1] = 0.0
    energyflux.eflx_traffic_patch[1] = 0.0; energyflux.eflx_ventilation_patch[1] = 0.0
    energyflux.eflx_lwrad_net_u_patch[1] = 0.0; energyflux.eflx_lwrad_out_u_patch[1] = 0.0

    waterfluxbulk.wf.qflx_evap_soi_patch[1] = 1.0e-5
    waterfluxbulk.wf.qflx_evap_veg_patch[1] = 2.0e-5
    waterfluxbulk.wf.qflx_tran_veg_patch[1] = 1.0e-5
    waterfluxbulk.wf.qflx_evap_tot_patch[1] = 3.0e-5
    waterfluxbulk.wf.qflx_evap_can_patch[1] = 1.0e-5
    waterfluxbulk.qflx_ev_snow_patch[1] = 5.0e-6
    waterfluxbulk.qflx_ev_soil_patch[1] = 5.0e-6
    waterfluxbulk.qflx_ev_h2osfc_patch[1] = 0.0

    canopystate.frac_veg_nosno_patch[1] = 0
    solarabs.sabg_soil_patch[1] = 100.0; solarabs.sabg_snow_patch[1] = 0.0; solarabs.sabg_patch[1] = 100.0

    S = (; energyflux, temperature, canopystate, waterstatebulk, waterdiagbulk,
           waterfluxbulk, solarabs, patch_data, col_data, lun_data)
    return (; np, nc, nl, S,
              forc_lwrad = FT[300.0], dtime = FT(1800.0))
end

run_sf!(S, m_c, m_p, m_u, forc_lwrad, dtime) = CLM.soil_fluxes!(
    S.energyflux, S.temperature, S.canopystate, S.waterstatebulk, S.waterdiagbulk,
    S.waterfluxbulk, S.solarabs, S.patch_data, S.col_data, S.lun_data,
    m_c, m_p, m_u, 1:length(m_c), 1:length(m_p), forc_lwrad, dtime)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for soil_fluxes! (whole driver)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)               # CPU reference
    B = build(FT)               # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    Sd = map(ad, B.S)
    dmask(m) = to(collect(Bool, m))

    if !(Sd.waterdiagbulk.frac_sno_eff_col isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    m_c = trues(H.nc); m_p = trues(H.np); m_u = falses(H.np)

    run_sf!(H.S, m_c, m_p, m_u, H.forc_lwrad, H.dtime)                       # CPU
    run_sf!(Sd, dmask(m_c), dmask(m_p), dmask(m_u), to(B.forc_lwrad), B.dtime)  # device

    checks = [
        ("eflx_sh_grnd",   H.S.energyflux.eflx_sh_grnd_patch,   Sd.energyflux.eflx_sh_grnd_patch),
        ("eflx_sh_tot",    H.S.energyflux.eflx_sh_tot_patch,    Sd.energyflux.eflx_sh_tot_patch),
        ("eflx_soil_grnd", H.S.energyflux.eflx_soil_grnd_patch, Sd.energyflux.eflx_soil_grnd_patch),
        ("eflx_soil_grnd_r",H.S.energyflux.eflx_soil_grnd_r_patch,Sd.energyflux.eflx_soil_grnd_r_patch),
        ("eflx_lh_tot",    H.S.energyflux.eflx_lh_tot_patch,    Sd.energyflux.eflx_lh_tot_patch),
        ("eflx_lh_grnd",   H.S.energyflux.eflx_lh_grnd_patch,   Sd.energyflux.eflx_lh_grnd_patch),
        ("eflx_lh_vegt",   H.S.energyflux.eflx_lh_vegt_patch,   Sd.energyflux.eflx_lh_vegt_patch),
        ("eflx_lwrad_out", H.S.energyflux.eflx_lwrad_out_patch, Sd.energyflux.eflx_lwrad_out_patch),
        ("eflx_lwrad_net", H.S.energyflux.eflx_lwrad_net_patch, Sd.energyflux.eflx_lwrad_net_patch),
        ("errsoi_patch",   H.S.energyflux.errsoi_patch,         Sd.energyflux.errsoi_patch),
        ("errsoi_col",     H.S.energyflux.errsoi_col,           Sd.energyflux.errsoi_col),
        ("qflx_evap_soi",  H.S.waterfluxbulk.wf.qflx_evap_soi_patch, Sd.waterfluxbulk.wf.qflx_evap_soi_patch),
        ("qflx_evap_tot",  H.S.waterfluxbulk.wf.qflx_evap_tot_patch, Sd.waterfluxbulk.wf.qflx_evap_tot_patch),
        ("qflx_liqevap",   H.S.waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch,
                           Sd.waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch),
        ("qflx_solidevap", H.S.waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch,
                           Sd.waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch),
        ("ev_snow",        H.S.waterfluxbulk.qflx_ev_snow_patch, Sd.waterfluxbulk.qflx_ev_snow_patch),
        ("t_skin",         H.S.temperature.t_skin_patch,        Sd.temperature.t_skin_patch),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-16s CPU reference is all-NaN/Inf — skipping (no parity signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-16s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  WHOLE soil_fluxes! MATCHES CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
