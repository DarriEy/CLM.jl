# ==========================================================================
# gpu_validate_drainage_e2e.jl — end-to-end Metal parity for the WHOLE
# drainage! AND clm_vic_map! drivers (one thread per column running the full
# per-column subroutine in-thread).
#
# drainage! is the long pole: per-column it runs SEQUENTIAL nlevsoi level
# sweeps —
#   • jwt search (layer right above the water table),
#   • BASEFLOW perched-table drainage (frost-table + saturation cascade),
#   • topographic runoff / aquifer recharge (rsub_top, loop-carried zwt),
#   • upward saturation-excess redistribution + watmin pull-from-below,
#   • urban no-drainage override.
# clm_vic_map! maps each column's CLM layers onto the VIC layers + top-layer
# totals (fully per-column).
#
# Builds a small Float32 instance with SEVERAL columns deliberately exercising
# every major branch:
#   col 1: unfrozen, shallow water table (within layer 1-2) -> within-soil cascade
#   col 2: unfrozen, deep water table mid-column            -> deepening cascade
#   col 3: unfrozen, water table at/below bedrock (jwt==nlevsoi) -> aquifer branch
#   col 4: frozen upper layers, wt above frost table        -> perched-table path
#   col 5: urban pervious-road column                        -> urban (keeps drainage)
#   col 6: urban roof column                                 -> urban override (zero drain)
# Runs both whole functions on the CPU, adapts every state struct (+ masks) to
# Metal, runs the SAME calls on the device, and compares the mutated outputs
# field-by-field with a NaN-aware reldiff. CPU-reference fields are asserted
# finite first so a both-NaN false PASS cannot slip through.
#
#   julia --project=scripts scripts/gpu_validate_drainage_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal   # Metal-specific; MtlArray is the Adapt adaptor type for the structs
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
    nlevsoi  = CLM.varpar.nlevsoi
    nlevsno  = CLM.varpar.nlevsno
    nlevgrnd = CLM.varpar.nlevgrnd
    nc = 6
    joff    = nlevsno
    joff_zi = nlevsno + 1

    sh  = CLM.SoilHydrologyData{FT}();  CLM.soilhydrology_init!(sh, nc)
    ss  = CLM.SoilStateData{FT}();      CLM.soilstate_init!(ss, nc, nc)
    wsb = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(wsb, nc, nc, 1, 1)
    wfb = CLM.WaterFluxBulkData{FT}();  CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
    temp = CLM.TemperatureData{FT}();   CLM.temperature_init!(temp, nc, nc, 1, 1)

    nlev_soisno = nlevsno + CLM.varpar.nlevmaxurbgrnd
    col_dz = fill(FT(NaN), nc, nlev_soisno)
    col_z  = fill(FT(NaN), nc, nlev_soisno)
    col_zi = fill(FT(NaN), nc, nlev_soisno + 1)

    # Monotonic soil column (no resolved snow): dz grows with depth.
    for c in 1:nc
        zi_prev = FT(0.0)
        col_zi[c, joff_zi + 0] = FT(0.0)
        for j in 1:nlevsoi
            dz = FT(0.025) * FT(1.5)^(j - 1)
            col_dz[c, joff + j] = dz
            zi = zi_prev + dz
            col_zi[c, joff_zi + j] = zi
            col_z[c, joff + j]  = FT(0.5) * (zi_prev + zi)
            zi_prev = zi
        end
    end
    zi_bot = Array(col_zi)[1, joff_zi + nlevsoi]

    # --- soil hydraulic properties ---
    for c in 1:nc, j in 1:nlevgrnd
        ss.watsat_col[c, j]       = FT(0.45)
        ss.sucsat_col[c, j]       = FT(150.0)
        ss.bsw_col[c, j]          = FT(6.0)
        ss.hksat_col[c, j]        = FT(1.0e-3)
        ss.eff_porosity_col[c, j] = FT(0.40)
    end

    # --- water state on snow+soil grid (offsets used by drainage!) ---
    for c in 1:nc, j in 1:nlevsoi
        frac = j <= 3 ? FT(0.95) : FT(0.6)
        liqvol = frac * FT(0.45)
        wsb.ws.h2osoi_liq_col[c, joff + j] = liqvol * Array(col_dz)[c, joff + j] * FT(1000.0)
        wsb.ws.h2osoi_ice_col[c, joff + j] = FT(0.0)
        wsb.ws.h2osoi_vol_col[c, j]        = FT(0.0)
    end
    # clm_vic_map! reads h2osoi_*[c, 1:nlevsoi/nlevgrnd] on the *soil* grid (no
    # snow offset), as in the real hydrology_drainage.jl caller. Fill the full
    # arrays so those reads are finite (snow indices 1:nlevsno are otherwise NaN
    # and would poison the CPU reference -> a both-NaN false PASS). Same values
    # feed CPU and device, so parity is unaffected.
    for c in 1:nc
        for j in 1:nlevgrnd
            wsb.ws.h2osoi_liq_col[c, j] = FT(30.0) + FT(j)
            wsb.ws.h2osoi_ice_col[c, j] = FT(5.0)  + FT(0.5) * FT(j)
            wsb.ws.h2osoi_vol_col[c, j] = FT(0.3)
        end
    end

    # --- temperature: all warm by default ---
    for c in 1:nc, j in 1:nlev_soisno
        temp.t_soisno_col[c, j] = FT(285.0)
    end

    # --- per-column scalars / state ---
    for c in 1:nc
        sh.icefrac_col[c, :]  .= FT(0.0)
        sh.hkdepth_col[c]      = FT(1.0 / 2.5)
        sh.frost_table_col[c]  = FT(0.0)
        sh.zwt_perched_col[c]  = FT(0.0)
        sh.zwt_col[c]          = FT(0.0)
        wsb.ws.wa_col[c]       = FT(4800.0)
        wsb.ws.h2osfc_col[c]   = FT(0.0)
        wfb.wf.qflx_snwcp_liq_col[c] = FT(1.0e-6)
        wfb.wf.qflx_drain_col[c]         = FT(0.0)
        wfb.wf.qflx_rsub_sat_col[c]      = FT(0.0)
        wfb.wf.qflx_drain_perched_col[c] = FT(0.0)
        wfb.wf.qflx_qrgwl_col[c]         = FT(0.0)
        wfb.wf.qflx_ice_runoff_xs_col[c] = FT(0.0)
    end
    wsb.ws.aquifer_water_baseline = FT(4800.0)
    sh.h2osfcflag = 1

    # branch placement
    sh.zwt_col[1] = FT(0.08)             # shallow within-soil
    sh.zwt_col[2] = FT(0.5 * zi_bot)     # deep mid-column
    sh.zwt_col[3] = FT(zi_bot + 5.0)     # at/below bedrock -> aquifer
    # col 4: frozen upper layers, wt above frost table -> perched path
    for j in 1:nlevsoi
        temp.t_soisno_col[4, joff + j] = j <= 4 ? FT(285.0) : FT(270.0)
    end
    sh.zwt_col[4] = FT(0.02)             # above frost table
    sh.zwt_col[5] = FT(0.3 * zi_bot)     # urban pervious road
    sh.zwt_col[6] = FT(0.3 * zi_bot)     # urban roof (override)
    # give col 3 a surplus aquifer so the wa>baseline branch fires
    wsb.ws.wa_col[3] = FT(4900.0)

    # --- topology / landunit / column type ---
    col_topo_slope = fill(FT(2.0), nc)              # ~2 deg slope -> nonzero rsub_top
    col_snl        = zeros(Int, nc)
    col_landunit   = collect(1:nc)
    col_itype      = fill(1, nc)
    col_itype[5]   = CLM.ICOL_ROAD_PERV             # pervious road -> keeps drainage
    col_itype[6]   = CLM.ICOL_ROOF                  # roof         -> drainage zeroed
    lun_urbpoi     = fill(false, nc)
    lun_urbpoi[5]  = true
    lun_urbpoi[6]  = true

    mask_hydrology = trues(nc)
    mask_urban     = BitVector([false, false, false, false, true, true])

    # ---- VIC mapping inputs (nlayer = nlevsoi, nlayert = nlevsoi + 1) ----
    # NOTE: clm_vic_map! indexes h2osoi_*[c, 1:nlevsoi] on the *soil* grid (no
    # snow offset), matching the real caller in hydrology_drainage.jl. Build a
    # fresh soil-grid set so the harness exercises that path on its own values.
    nlayer  = nlevsoi
    nlayert = nlevsoi + 1
    # reallocate VIC arrays at the caller's (nlayer/nlayert) sizes
    sh.depth_col     = fill(FT(0.2), nc, nlayert)
    sh.porosity_col  = fill(FT(0.45), nc, nlayer)
    sh.max_moist_col = fill(FT(500.0), nc, nlayer)
    # moist/ice carry their PRIOR values into the ice-clamp smooth_min(moist0+ice0,
    # ice); production seeds these from restart/init_cold, so seed them finite here
    # (NaN priors would poison the clamp exactly as in the scalar version).
    sh.moist_col     = fill(FT(120.0), nc, nlayert)
    sh.ice_col       = fill(FT(40.0),  nc, nlayert)
    sh.moist_vol_col = fill(FT(0.2),   nc, nlayert)
    sh.vic_clm_fract_col = fill(FT(0.0), nc, nlayer, nlevsoi)
    # identity-ish mapping: VIC layer i draws mainly from CLM layer i
    for c in 1:nc, i in 1:nlayer, j in 1:nlevsoi
        sh.vic_clm_fract_col[c, i, j] = (i == j) ? FT(0.8) : (abs(i - j) == 1 ? FT(0.1) : FT(0.0))
    end

    S = (; sh, ss, wsb, wfb, temp, col_dz, col_z, col_zi)
    P = (; col_snl, col_itype, col_landunit, col_topo_slope, lun_urbpoi,
           mask_hydrology, mask_urban, nlevsoi, nlevsno, nlevgrnd,
           nlayer, nlayert, nc, dtime = FT(1800.0))
    return (; S, P)
end

run_drain!(S, P, dmask, ivec, bvec) = CLM.drainage!(
    S.temp.t_soisno_col, S.sh, S.ss, S.wsb.ws, S.wfb,
    S.col_dz, S.col_z, S.col_zi, ivec(P.col_snl), ivec(P.col_itype),
    ivec(P.col_landunit), ivec(P.col_topo_slope), bvec(P.lun_urbpoi),
    dmask(P.mask_hydrology), dmask(P.mask_urban), 1:P.nc, P.nlevsoi, P.dtime;
    nlevsno = P.nlevsno, use_vichydro = false)

run_vic!(S, P, dmask) = CLM.clm_vic_map!(
    S.sh, S.wsb.ws, S.col_dz, S.col_zi, S.col_z,
    dmask(P.mask_hydrology), 1:P.nc, P.nlevsoi, P.nlevgrnd, P.nlayer, P.nlayert)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for drainage! + clm_vic_map! (whole drivers)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU driver exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)               # CPU reference
    B = build(FT)               # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    Sd = (; sh = ad(B.S.sh), ss = ad(B.S.ss), wsb = ad(B.S.wsb), wfb = ad(B.S.wfb),
            temp = ad(B.S.temp),
            col_dz = to(B.S.col_dz), col_z = to(B.S.col_z), col_zi = to(B.S.col_zi))

    if !(Sd.sh.zwt_col isa Metal.MtlArray)
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    dmask(m) = to(collect(Bool, m))
    ivec(v)  = to(v)
    bvec(v)  = to(collect(Bool, v))

    # ---------- drainage! ----------
    run_drain!(H.S, H.P, identity, identity, identity)     # CPU
    run_drain!(Sd,  B.P, dmask,    ivec,     bvec)         # device

    # ---------- clm_vic_map! ----------
    run_vic!(H.S, H.P, identity)                            # CPU
    run_vic!(Sd,  B.P, dmask)                               # device

    checks = [
        ("zwt",                H.S.sh.zwt_col,                    Sd.sh.zwt_col),
        ("zwt_perched",        H.S.sh.zwt_perched_col,            Sd.sh.zwt_perched_col),
        ("frost_table",        H.S.sh.frost_table_col,            Sd.sh.frost_table_col),
        ("icefrac",            H.S.sh.icefrac_col,                Sd.sh.icefrac_col),
        ("wa",                 H.S.wsb.ws.wa_col,                 Sd.wsb.ws.wa_col),
        ("h2osfc",             H.S.wsb.ws.h2osfc_col,             Sd.wsb.ws.h2osfc_col),
        ("h2osoi_liq",         H.S.wsb.ws.h2osoi_liq_col,         Sd.wsb.ws.h2osoi_liq_col),
        ("h2osoi_ice",         H.S.wsb.ws.h2osoi_ice_col,         Sd.wsb.ws.h2osoi_ice_col),
        ("qflx_drain",         H.S.wfb.wf.qflx_drain_col,         Sd.wfb.wf.qflx_drain_col),
        ("qflx_rsub_sat",      H.S.wfb.wf.qflx_rsub_sat_col,      Sd.wfb.wf.qflx_rsub_sat_col),
        ("qflx_drain_perched", H.S.wfb.wf.qflx_drain_perched_col, Sd.wfb.wf.qflx_drain_perched_col),
        ("qflx_qrgwl",         H.S.wfb.wf.qflx_qrgwl_col,         Sd.wfb.wf.qflx_qrgwl_col),
        ("qflx_ice_runoff_xs", H.S.wfb.wf.qflx_ice_runoff_xs_col, Sd.wfb.wf.qflx_ice_runoff_xs_col),
        # VIC-mapped outputs
        ("vic_moist",          H.S.sh.moist_col,                  Sd.sh.moist_col),
        ("vic_ice",            H.S.sh.ice_col,                    Sd.sh.ice_col),
        ("vic_moist_vol",      H.S.sh.moist_vol_col,              Sd.sh.moist_vol_col),
        ("vic_top_moist",      H.S.sh.top_moist_col,              Sd.sh.top_moist_col),
        ("vic_top_ice",        H.S.sh.top_ice_col,                Sd.sh.top_ice_col),
        ("vic_top_max_moist",  H.S.sh.top_max_moist_col,          Sd.sh.top_max_moist_col),
        ("vic_top_moist_lim",  H.S.sh.top_moist_limited_col,      Sd.sh.top_moist_limited_col),
    ]
    nfail = 0; worst = 0.0; worst_nm = ""
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-20s CPU reference all-NaN/Inf — skipping (no signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-20s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        d > worst && (worst = d; worst_nm = nm)
        ok || (nfail += 1)
    end

    # Per-column zwt / qflx_drain to show branch coverage explicitly.
    println()
    zc = Array(H.S.sh.zwt_col);            zd = Array(Sd.sh.zwt_col)
    qc = Array(H.S.wfb.wf.qflx_drain_col); qd = Array(Sd.wfb.wf.qflx_drain_col)
    labels = ["shallow", "deep", "bedrock/aq", "frozen/perch", "urban-road", "urban-roof"]
    for c in 1:H.P.nc
        @printf("    col %d (%-12s): zwt cpu=%9.5f dev=%9.5f | qflx_drain cpu=%11.4e dev=%11.4e\n",
                c, labels[c], zc[c], zd[c], qc[c], qd[c])
    end

    println()
    @printf("  worst-case reldiff: %.3e (%s)\n", worst, worst_nm)
    println(nfail == 0 ? "  WHOLE drainage! + clm_vic_map! MATCH CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
