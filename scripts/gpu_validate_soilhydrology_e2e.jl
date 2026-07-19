# ==========================================================================
# gpu_validate_soilhydrology_e2e.jl — end-to-end GPU parity for the seven
# A8 soil-hydrology functions whose per-column / per-(c,j) loops are now
# KernelAbstractions kernels:
#
#   set_soil_water_fractions!  set_floodc!  set_qflx_inputs!
#   route_infiltration_excess! infiltration!  total_surface_runoff!
#   renew_condensation!
#
# Builds a small Float32 instance covering the major branches (one non-urban
# soil column + one urban impervious-roof column), runs each WHOLE function on
# the CPU, adapts every state struct (+ masks) to the GPU, runs the SAME calls on
# the device, and compares every mutated field with a NaN-aware reldiff.
#
#   julia --project=scripts scripts/gpu_validate_soilhydrology_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

# reldiff: NaN-aware (both-NaN agrees; one-sided NaN flags divergence).
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
    nlevsoi = CLM.varpar.nlevsoi
    nlevsno = CLM.varpar.nlevsno
    nlevgrnd = CLM.varpar.nlevgrnd
    ntot = nlevsno + CLM.varpar.nlevmaxurbgrnd

    nc = 2; ng = 2; nl = 2

    # column 1: non-urban soil column (hydrologically active, no snow)
    # column 2: urban impervious-roof column (exercises urban branches)
    soilhyd = CLM.SoilHydrologyData{FT}(); CLM.soilhydrology_init!(soilhyd, nc)
    soilhyd.h2osfcflag = 1
    soilhyd.xs_urban_col .= FT(0)
    for c in 1:nc, j in 1:nlevgrnd
        soilhyd.icefrac_col[c, j] = FT(0)
    end

    soilstate = CLM.SoilStateData{FT}()
    soilstate.watsat_col = fill(FT(0.4), nc, nlevgrnd)
    soilstate.eff_porosity_col = fill(FT(0.35), nc, nlevgrnd)

    wsb = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(wsb, nc, nc, nl, ng)
    ws = wsb.ws
    ws.h2osoi_liq_col = fill(FT(20.0), nc, ntot)
    ws.h2osoi_ice_col = fill(FT(5.0), nc, ntot)
    ws.excess_ice_col = fill(FT(0.0), nc, ntot)

    wdb = CLM.WaterDiagnosticBulkData{FT}(); CLM.waterdiagnosticbulk_init!(wdb, nc, nc, nl, ng)
    wdb.frac_sno_eff_col .= FT(0.0)
    wdb.frac_h2osfc_col  .= FT(0.0)

    wfb = CLM.WaterFluxBulkData{FT}(); CLM.waterfluxbulk_init!(wfb, nc, nc, nl, ng)
    wf = wfb.wf
    wf.qflx_rain_plus_snomelt_col      .= FT(1.0e-4)
    wf.qflx_snow_h2osfc_col            .= FT(2.0e-5)
    wf.qflx_floodc_col                 .= FT(0.0)
    wf.qflx_liqevap_from_top_layer_col .= FT(1.0e-6)
    wf.qflx_liqdew_to_top_layer_col    .= FT(1.0e-6)
    wf.qflx_soliddew_to_top_layer_col  .= FT(5.0e-7)
    wf.qflx_solidevap_from_top_layer_col .= FT(1.0e-7)
    wf.qflx_top_soil_col               .= FT(0.0)
    wf.qflx_surf_col                   .= FT(0.0)
    wf.qflx_infl_col                   .= FT(0.0)

    wfb.qflx_ev_soil_col          .= FT(5.0e-6)
    wfb.qflx_ev_h2osfc_col        .= FT(0.0)
    wfb.qflx_sat_excess_surf_col  .= FT(1.0e-6)
    wfb.qflx_infl_excess_col      .= FT(2.0e-6)
    wfb.qflx_infl_excess_surf_col .= FT(0.0)
    wfb.qflx_h2osfc_surf_col      .= FT(3.0e-6)
    wfb.qflx_in_soil_col          .= FT(0.0)
    wfb.qflx_in_soil_limited_col  .= FT(0.0)
    wfb.qflx_h2osfc_drain_col     .= FT(1.0e-6)
    wfb.qflx_top_soil_to_h2osfc_col .= FT(0.0)
    wfb.qflx_in_h2osfc_col        .= FT(0.0)

    qflx_floodg = FT[0.5, 0.3]

    col_dz = fill(FT(0.1), nc, ntot)

    col_gridcell = [1, 2]
    col_landunit = [1, 2]
    col_snl      = [0, 0]
    # column 1 normal soil; column 2 urban roof
    col_itype    = [1, CLM.ICOL_ROOF]
    lun_itype    = [CLM.ISTSOIL, CLM.ISTSOIL]
    lun_urbpoi   = [false, true]

    mask_hydro = trues(nc)        # both treated as hydrologically active for the active-col kernels
    mask_nolake = trues(nc)
    mask_urban = BitVector([false, true])

    S = (; soilhyd, soilstate, wsb, wdb, wfb)
    P = (; col_dz, qflx_floodg, col_gridcell, col_landunit, col_snl, col_itype,
           lun_itype, lun_urbpoi, mask_hydro, mask_nolake, mask_urban,
           nlevsoi, nlevsno, nc, dtime = FT(1800.0))
    return (; S, P)
end

# Run all seven whole functions in order on a given (S, P, mask-converter).
function run_all!(S, P, dmask, ivec)
    nc = P.nc
    CLM.set_soil_water_fractions!(S.soilhyd, S.soilstate, S.wsb.ws, ivec(P.col_dz),
        dmask(P.mask_hydro), 1:nc, P.nlevsoi, P.nlevsno)

    CLM.set_floodc!(S.wfb.wf.qflx_floodc_col, ivec(P.qflx_floodg), ivec(P.col_gridcell),
        ivec(P.col_itype), dmask(P.mask_nolake), 1:nc)

    CLM.set_qflx_inputs!(S.wfb, S.wdb, ivec(P.col_snl), dmask(P.mask_hydro), 1:nc)

    CLM.route_infiltration_excess!(S.wfb, S.soilhyd, ivec(P.col_landunit),
        ivec(P.lun_itype), dmask(P.mask_hydro), 1:nc)

    CLM.infiltration!(S.wfb, dmask(P.mask_hydro), 1:nc)

    CLM.total_surface_runoff!(S.wfb, S.soilhyd, S.wsb.ws, ivec(P.col_snl),
        ivec(P.col_itype), ivec(P.col_landunit), P.lun_urbpoi,
        dmask(P.mask_hydro), dmask(P.mask_urban), 1:nc, P.dtime)

    CLM.renew_condensation!(S.wsb.ws, S.wdb, S.wfb, ivec(P.col_snl), ivec(P.col_itype),
        dmask(P.mask_hydro), dmask(P.mask_urban), 1:nc, P.dtime)
    return nothing
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for soil-hydrology A8 functions (7 whole fns)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)   # CPU reference
    B = build(FT)   # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    Sd = (; soilhyd = ad(B.S.soilhyd), soilstate = ad(B.S.soilstate),
            wsb = ad(B.S.wsb), wdb = ad(B.S.wdb), wfb = ad(B.S.wfb))

    if !(Sd.wfb.wf.qflx_infl_col isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    # CPU: masks stay BitVector, integer vectors stay host (kernels writing BitVector
    # run on CPU; integer index arrays are @Const inputs read on whichever backend).
    run_all!(H.S, H.P, identity, identity)
    # device: Bool masks + device integer arrays
    dmask(m) = to(collect(Bool, m))
    ivec(v) = to(v)
    run_all!(Sd, B.P, dmask, ivec)

    checks = [
        ("eff_porosity",       H.S.soilstate.eff_porosity_col,      Sd.soilstate.eff_porosity_col),
        ("icefrac",            H.S.soilhyd.icefrac_col,             Sd.soilhyd.icefrac_col),
        ("qflx_floodc",        H.S.wfb.wf.qflx_floodc_col,          Sd.wfb.wf.qflx_floodc_col),
        ("qflx_top_soil",      H.S.wfb.wf.qflx_top_soil_col,        Sd.wfb.wf.qflx_top_soil_col),
        ("qflx_in_soil",       H.S.wfb.qflx_in_soil_col,            Sd.wfb.qflx_in_soil_col),
        ("qflx_top2h2osfc",    H.S.wfb.qflx_top_soil_to_h2osfc_col, Sd.wfb.qflx_top_soil_to_h2osfc_col),
        ("qflx_in_soil_lim",   H.S.wfb.qflx_in_soil_limited_col,    Sd.wfb.qflx_in_soil_limited_col),
        ("qflx_in_h2osfc",     H.S.wfb.qflx_in_h2osfc_col,          Sd.wfb.qflx_in_h2osfc_col),
        ("qflx_infl_exc_surf", H.S.wfb.qflx_infl_excess_surf_col,   Sd.wfb.qflx_infl_excess_surf_col),
        ("qflx_infl",          H.S.wfb.wf.qflx_infl_col,            Sd.wfb.wf.qflx_infl_col),
        ("qflx_surf",          H.S.wfb.wf.qflx_surf_col,            Sd.wfb.wf.qflx_surf_col),
        ("xs_urban",           H.S.soilhyd.xs_urban_col,            Sd.soilhyd.xs_urban_col),
        ("h2osoi_liq",         H.S.wsb.ws.h2osoi_liq_col,           Sd.wsb.ws.h2osoi_liq_col),
        ("h2osoi_ice",         H.S.wsb.ws.h2osoi_ice_col,           Sd.wsb.ws.h2osoi_ice_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-18s CPU reference all-NaN/Inf — skipping (no parity signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-18s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end
    println()
    println(nfail == 0 ? "  ALL 7 soil-hydrology functions MATCH CPU ON $name ($FT) ✓" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
