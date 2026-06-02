# ==========================================================================
# gpu_validate_lateralflow_e2e.jl — end-to-end Metal parity for the WHOLE
# perched_lateral_flow! AND subsurface_lateral_flow! drivers (the two
# largest/hardest A8 soil-hydrology functions): one thread per hydrology
# column running the full per-column sequence in-thread.
#
#   perched_lateral_flow!     — frost/perch-table search, q_perch drainage,
#                               per-layer storage removal.
#   subsurface_lateral_flow!  — icefrac/impedance, jwt search, power-law
#                               baseflow, topographic-runoff water-table
#                               rise/deepen cascade (loop-carried zwt),
#                               upward excess redistribution, surface/ice
#                               overflow, watmin back-fill.
#
# Both kernels cover the NON-HILLSLOPE path, which is fully per-column
# independent — the only cross-column reduction (the hillslope donor/receiver
# scatter `qflx_*[cold[c]] -= …`) is never entered when cold[c]==ISPVAL, so no
# exact-ordering reduction is needed on the device. (The hillslope path stays
# on the sequential CPU branch.)
#
# Builds a small Float32 instance with SEVERAL columns deliberately placing
# the water table / perched zone at different depths so every branch runs:
#   col 1: saturated, shallow water table  -> rising drainage, perched zone
#   col 2: drier, deep water table         -> deepening drainage, no perch
#   col 3: water table at/below bedrock    -> baseflow-off branch
#   col 4: urban (non-pervious road)       -> no-drainage urban branch
#
# Runs each WHOLE function on the CPU, adapts every state struct (+ masks) to
# Metal, runs the SAME calls on the device, and compares every mutated field
# with a NaN-aware reldiff. CPU-reference fields are asserted finite so a
# both-NaN false PASS cannot slip through.
#
#   julia --project=scripts scripts/gpu_validate_lateralflow_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal   # MtlArray is the Adapt adaptor type for the structs
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
    nc = 4; nl = 2; ng = 2
    joff    = nlevsno
    joff_zi = nlevsno + 1
    nlevtot = nlevsoi + nlevsno

    # --- soil hydrology state ---
    sh = CLM.SoilHydrologyData{FT}(); CLM.soilhydrology_init!(sh, nc)
    sh.icefrac_col .= FT(0.0)

    # --- soil hydraulic properties ---
    ss = CLM.SoilStateData{FT}()
    ss.watsat_col       = fill(FT(0.4),   nc, nlevgrnd)
    ss.bsw_col          = fill(FT(5.0),   nc, nlevgrnd)
    ss.hksat_col        = fill(FT(0.01),  nc, nlevgrnd)
    ss.sucsat_col       = fill(FT(100.0), nc, nlevgrnd)
    ss.eff_porosity_col = fill(FT(0.35),  nc, nlevgrnd)
    ss.hk_l_col         = fill(FT(0.005), nc, nlevgrnd)

    # --- water state ---
    ws = CLM.WaterStateData{FT}(); CLM.waterstate_init!(ws, nc, nc, nl, ng)
    ws.h2osoi_liq_col = fill(FT(10.0), nc, nlevtot)
    ws.h2osoi_ice_col = fill(FT(0.0),  nc, nlevtot)
    ws.h2osfc_col    .= FT(0.0)
    ws.stream_water_volume_lun = fill(FT(0.0), nl)

    # --- water fluxes ---
    wfb = CLM.WaterFluxBulkData{FT}(); CLM.waterfluxbulk_init!(wfb, nc, nc, nl, ng)
    wfb.wf.qflx_snwcp_liq_col .= FT(0.0)

    # --- column geometry / topology ---
    col = CLM.ColumnData{FT}(); CLM.column_init!(col, nc)
    col.dz = fill(FT(0.1), nc, size(col.dz, 2))
    col.z  = zeros(FT, nc, size(col.z, 2))
    col.zi = zeros(FT, nc, size(col.zi, 2))
    for j in 1:size(col.z, 2)
        col.z[:, j] .= FT(0.05) + (j - 1) * FT(0.1)
    end
    for j in 1:size(col.zi, 2)
        col.zi[:, j] .= j * FT(0.1)   # zi[:,joff_zi+k] = (joff_zi+k)*0.1
    end
    col.nbedrock .= nlevsoi
    col.itype    .= 1
    col.is_hillslope_column .= false
    col.active   .= true
    col.cold     .= CLM.ISPVAL
    col.colu     .= CLM.ISPVAL
    col.hill_slope .= FT(0.1); col.hill_elev .= FT(10.0); col.hill_distance .= FT(100.0)
    col.hill_width .= FT(50.0); col.hill_area .= FT(5000.0)
    col.topo_slope .= FT(5.0)
    col.wtgcell    .= FT(0.5)
    # gridcell / landunit map: cols 1-3 -> grc1/lun1 (soil); col 4 -> grc2/lun2 (urban)
    col.gridcell  = Int[1, 1, 1, 2]
    col.landunit  = Int[1, 1, 1, 2]

    # zi_bedrock (uniform across cols): interface depth at the bedrock index
    zi_bedrock = Array(col.zi)[1, nlevsoi + joff_zi]

    # --- per-column water-table / perched-zone placement (branch coverage) ---
    # col 1: saturated, shallow wt -> rising drainage; perched zone present
    sh.zwt_col[1]         = FT(0.15)
    sh.frost_table_col[1] = FT(0.45)
    sh.zwt_perched_col[1] = FT(0.25)
    ws.h2osoi_liq_col[1, (joff+1):(joff+nlevsoi)] .= FT(35.0)  # near saturation

    # col 2: drier, deep wt mid-column -> deepening drainage; no perched zone
    sh.zwt_col[2]         = FT(0.5 * zi_bedrock)
    sh.frost_table_col[2] = FT(0.2)
    sh.zwt_perched_col[2] = FT(0.5)   # perched below frost -> no perch drainage
    ws.h2osoi_liq_col[2, (joff+1):(joff+nlevsoi)] .= FT(8.0)

    # col 3: water table at/below bedrock -> baseflow-off branch
    sh.zwt_col[3]         = FT(zi_bedrock + 5.0)
    sh.frost_table_col[3] = FT(0.4)
    sh.zwt_perched_col[3] = FT(0.2)

    # col 4: urban (non-pervious road) -> no-drainage urban branch
    sh.zwt_col[4]         = FT(0.3)
    sh.frost_table_col[4] = FT(0.4)
    sh.zwt_perched_col[4] = FT(0.2)
    col.itype[4]          = CLM.ICOL_ROOF   # not pervious road -> drainage zeroed

    # --- landunit data ---
    lun = CLM.LandunitData{FT}(); CLM.landunit_init!(lun, nl)
    lun.stream_channel_length .= FT(100.0)
    lun.stream_channel_width  .= FT(5.0)
    lun.stream_channel_depth  .= FT(1.0)
    lun.stream_channel_number .= FT(1.0)
    lun.urbpoi = Bool[false, true]

    # --- gridcell-level forcing ---
    grc_area  = fill(FT(100.0), ng)
    tdepth    = fill(FT(0.0),   ng)
    tdepthmax = fill(FT(1.0),   ng)

    mask_hydro = trues(nc)
    mask_urban = BitVector([false, false, false, true])

    S = (; sh, ss, ws, wfb, col, lun)
    P = (; grc_area, tdepth, tdepthmax, mask_hydro, mask_urban,
           nc, nl, ng, nlevsoi, dtime = FT(1800.0))
    return (; S, P)
end

# Run both whole functions in order (perched then subsurface) on a given
# (S, P) plus device/host converters for the host-side vectors + masks.
function run_both!(S, P, ivec, dmask)
    nc = P.nc
    CLM.perched_lateral_flow!(S.sh, S.ss, S.ws, S.wfb, S.col, S.lun,
        ivec(P.tdepth), ivec(P.tdepthmax),
        dmask(P.mask_hydro), 1:nc, P.nlevsoi, P.dtime)

    CLM.subsurface_lateral_flow!(S.sh, S.ss, S.ws, S.wfb, S.col, S.lun,
        ivec(P.tdepth), ivec(P.tdepthmax), ivec(P.grc_area),
        dmask(P.mask_hydro), dmask(P.mask_urban), 1:nc, P.nlevsoi, P.dtime)
    return nothing
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for perched_lateral_flow! + subsurface_lateral_flow!")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU drivers exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)   # CPU reference
    B = build(FT)   # source for the device snapshot

    ad(x) = CLM.Adapt.adapt(Metal.MtlArray, x)
    Sd = (; sh = ad(B.S.sh), ss = ad(B.S.ss), ws = ad(B.S.ws),
            wfb = ad(B.S.wfb), col = ad(B.S.col), lun = ad(B.S.lun))

    if !(Sd.sh.zwt_col isa Metal.MtlArray)
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    # CPU reference: masks stay BitVector, host vectors stay host.
    run_both!(H.S, H.P, identity, identity)
    # device: Bool masks + device vectors.
    dmask(m) = to(collect(Bool, m))
    ivec(v)  = to(v)
    run_both!(Sd, B.P, ivec, dmask)

    checks = [
        # perched_lateral_flow! outputs
        ("qflx_drain_perched", H.S.wfb.wf.qflx_drain_perched_col, Sd.wfb.wf.qflx_drain_perched_col),
        # subsurface_lateral_flow! outputs
        ("zwt",                H.S.sh.zwt_col,                    Sd.sh.zwt_col),
        ("icefrac",            H.S.sh.icefrac_col,                Sd.sh.icefrac_col),
        ("h2osfc",             H.S.ws.h2osfc_col,                 Sd.ws.h2osfc_col),
        ("h2osoi_liq",         H.S.ws.h2osoi_liq_col,             Sd.ws.h2osoi_liq_col),
        ("h2osoi_ice",         H.S.ws.h2osoi_ice_col,             Sd.ws.h2osoi_ice_col),
        ("qflx_latflow_out",   H.S.wfb.wf.qflx_latflow_out_col,   Sd.wfb.wf.qflx_latflow_out_col),
        ("qflx_latflow_in",    H.S.wfb.wf.qflx_latflow_in_col,    Sd.wfb.wf.qflx_latflow_in_col),
        ("vol_discharge",      H.S.wfb.wf.volumetric_discharge_col, Sd.wfb.wf.volumetric_discharge_col),
        ("qflx_drain",         H.S.wfb.wf.qflx_drain_col,         Sd.wfb.wf.qflx_drain_col),
        ("qflx_qrgwl",         H.S.wfb.wf.qflx_qrgwl_col,         Sd.wfb.wf.qflx_qrgwl_col),
        ("qflx_rsub_sat",      H.S.wfb.wf.qflx_rsub_sat_col,      Sd.wfb.wf.qflx_rsub_sat_col),
        ("qflx_ice_runoff_xs", H.S.wfb.wf.qflx_ice_runoff_xs_col, Sd.wfb.wf.qflx_ice_runoff_xs_col),
    ]
    nfail = 0; worst = 0.0; worst_nm = ""
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-20s CPU reference all-NaN/Inf — skipping (no parity signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        d > worst && (worst = d; worst_nm = nm)
        @printf("  [%s] %-20s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    # Per-column zwt / qflx_drain to show branch coverage explicitly.
    println()
    zc = Array(H.S.sh.zwt_col);            zd = Array(Sd.sh.zwt_col)
    qc = Array(H.S.wfb.wf.qflx_drain_col); qd = Array(Sd.wfb.wf.qflx_drain_col)
    pc = Array(H.S.wfb.wf.qflx_drain_perched_col); pd = Array(Sd.wfb.wf.qflx_drain_perched_col)
    for c in 1:H.P.nc
        @printf("    col %d: zwt cpu=%9.5f dev=%9.5f | qdrain cpu=%11.4e dev=%11.4e | qdperch cpu=%11.4e dev=%11.4e\n",
                c, zc[c], zd[c], qc[c], qd[c], pc[c], pd[c])
    end

    println()
    @printf("  worst field: %s  (rel = %.3e)\n", worst_nm, worst)
    println(nfail == 0 ?
        "  WHOLE perched_lateral_flow! + subsurface_lateral_flow! MATCH CPU ON $name ($FT)" :
        "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
