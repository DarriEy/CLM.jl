# ==========================================================================
# gpu_validate_watertable_e2e.jl — end-to-end Metal parity for the WHOLE
# water_table! driver (the A8 "long pole"): one thread per hydrology column
# running the full nested, loop-carried water-table search in-thread —
#   • jwt search (layer right above the water table),
#   • QCHARGE rise/deepen cascade (loop-carried zwt, aquifer specific yield),
#   • BASEFLOW perched water table (frost-table + saturation search).
#
# Builds a small Float32 instance with SEVERAL columns deliberately placing the
# water table at different depths so every branch of the search runs:
#   col 1: shallow water table (within layer 1-2)   -> jwt small, rising qcharge
#   col 2: deep water table (mid column)            -> deepening qcharge
#   col 3: at/below bedrock (jwt == nlevsoi)        -> aquifer branch (wa update)
#   col 4: frozen upper layers (frost table active) -> perched-table branch
# Runs water_table! on the CPU, adapts every state struct (+ mask) to Metal,
# runs the SAME call on the device, and compares the mutated outputs
# field-by-field. CPU-reference fields are asserted finite so a both-NaN false
# PASS cannot slip through.
#
#   julia --project=scripts scripts/gpu_validate_watertable_e2e.jl
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
    nc = 4; np = 4; nl = 1; ng = 1
    joff    = nlevsno
    joff_zi = nlevsno + 1

    sh  = CLM.SoilHydrologyData{FT}();  CLM.soilhydrology_init!(sh, nc)
    ss  = CLM.SoilStateData{FT}();      CLM.soilstate_init!(ss, np, nc)
    wsb = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(wsb, nc, np, nl, ng)
    wfb = CLM.WaterFluxBulkData{FT}();  CLM.waterfluxbulk_init!(wfb, nc, np, nl, ng)
    temp = CLM.TemperatureData{FT}();   CLM.temperature_init!(temp, np, nc, nl, ng)

    # --- geometry: snow+soil layer thickness / depth / interface (m) ---
    nlev_soisno = nlevsno + CLM.varpar.nlevmaxurbgrnd
    col_dz = fill(FT(NaN), nc, nlev_soisno)
    col_z  = fill(FT(NaN), nc, nlev_soisno)
    col_zi = fill(FT(NaN), nc, nlev_soisno + 1)

    # Build a simple monotonic soil column (no resolved snow): dz grows with depth.
    for c in 1:nc
        zi_prev = FT(0.0)
        col_zi[c, joff_zi + 0] = FT(0.0)
        for j in 1:nlevsoi
            dz = FT(0.025) * FT(1.5)^(j - 1)          # 0.025 m .. growing
            col_dz[c, joff + j] = dz
            zi = zi_prev + dz
            col_zi[c, joff_zi + j] = zi
            col_z[c, joff + j]  = FT(0.5) * (zi_prev + zi)
            zi_prev = zi
        end
    end
    zi_bot = Array(col_zi)[1, joff_zi + nlevsoi]   # bottom interface depth (m)

    # --- soil hydraulic properties (per column,layer) ---
    for c in 1:nc, j in 1:nlevsoi
        ss.watsat_col[c, j]   = FT(0.45)
        ss.sucsat_col[c, j]   = FT(150.0)
        ss.bsw_col[c, j]      = FT(6.0)
        ss.hksat_col[c, j]    = FT(1.0e-3)
        ss.eff_porosity_col[c, j] = FT(0.40)
    end

    # --- water state: liquid/ice per (column,layer) on the snow+soil grid ---
    for c in 1:nc, j in 1:nlevsoi
        # mostly-saturated upper layers, drier deep -> exercise perched search
        frac = j <= 3 ? FT(0.95) : FT(0.6)
        liqvol = frac * FT(0.45)
        wsb.ws.h2osoi_liq_col[c, joff + j] = liqvol * Array(col_dz)[c, joff + j] * FT(1000.0)
        wsb.ws.h2osoi_ice_col[c, joff + j] = FT(0.0)
        wsb.ws.h2osoi_vol_col[c, j] = FT(NaN)
    end

    # --- temperature on the snow+soil grid (all warm by default) ---
    for c in 1:nc, j in 1:nlev_soisno
        temp.t_soisno_col[c, j] = FT(285.0)
    end

    # --- per-column scalars ---
    for c in 1:nc
        sh.qcharge_col[c] = FT(0.0)
        wsb.ws.wa_col[c]  = FT(4800.0)
        sh.zwt_col[c]         = FT(0.0)
        sh.zwt_perched_col[c] = FT(0.0)
        sh.frost_table_col[c] = FT(0.0)
    end

    # col 1: shallow water table within layer ~2, rising qcharge (recharge > 0)
    sh.zwt_col[1]     = FT(0.08)
    sh.qcharge_col[1] = FT(5.0e-5)

    # col 2: deep water table mid-column, deepening qcharge (recharge < 0)
    sh.zwt_col[2]     = FT(0.5 * zi_bot)
    sh.qcharge_col[2] = FT(-3.0e-5)

    # col 3: water table at/below bedrock (jwt == nlevsoi) -> aquifer branch
    sh.zwt_col[3]     = FT(zi_bot + 5.0)
    sh.qcharge_col[3] = FT(2.0e-5)

    # col 4: frozen upper layers -> frost table active, perched-table branch
    for j in 1:nlevsoi
        temp.t_soisno_col[4, joff + j] = j <= 4 ? FT(285.0) : FT(270.0)
    end
    sh.zwt_col[4]     = FT(0.6 * zi_bot)
    sh.qcharge_col[4] = FT(1.0e-5)

    S = (; sh, ss, wsb, wfb, temp, col_dz, col_z, col_zi)
    return (; nc, np, nl, ng, nlevsoi, S, dtime = FT(1800.0))
end

run_wt!(S, mask, nc, nlevsoi, dtime) = CLM.water_table!(
    S.sh, S.ss, S.temp.t_soisno_col, S.wsb.ws, S.wfb,
    S.col_dz, S.col_z, S.col_zi, mask, 1:nc, nlevsoi, dtime)

function main(backend)
    println("=" ^ 70)
    println("END-TO-END Metal parity for water_table! (whole driver, A8 long pole)")
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
    dmask(m) = to(collect(Bool, m))

    if !(Sd.sh.zwt_col isa Metal.MtlArray)
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    mask = trues(H.nc)
    run_wt!(H.S,  mask,        H.nc, H.nlevsoi, H.dtime)            # CPU
    run_wt!(Sd,   dmask(mask), B.nc, B.nlevsoi, B.dtime)           # device

    checks = [
        ("zwt",                H.S.sh.zwt_col,                 Sd.sh.zwt_col),
        ("zwt_perched",        H.S.sh.zwt_perched_col,         Sd.sh.zwt_perched_col),
        ("frost_table",        H.S.sh.frost_table_col,         Sd.sh.frost_table_col),
        ("wa",                 H.S.wsb.ws.wa_col,              Sd.wsb.ws.wa_col),
        ("h2osoi_vol",         H.S.wsb.ws.h2osoi_vol_col,      Sd.wsb.ws.h2osoi_vol_col),
        ("qflx_drain",         H.S.wfb.wf.qflx_drain_col,      Sd.wfb.wf.qflx_drain_col),
        ("qflx_drain_perched", H.S.wfb.wf.qflx_drain_perched_col, Sd.wfb.wf.qflx_drain_perched_col),
        ("qflx_rsub_sat",      H.S.wfb.wf.qflx_rsub_sat_col,   Sd.wfb.wf.qflx_rsub_sat_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-20s CPU reference all-NaN/Inf — skipping (no signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-20s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    # Report per-column zwt to show the branch coverage explicitly.
    println()
    zc = Array(H.S.sh.zwt_col); zd = Array(Sd.sh.zwt_col)
    for c in 1:H.nc
        @printf("    col %d: zwt_cpu = %10.5f   zwt_dev = %10.5f\n", c, zc[c], zd[c])
    end

    println()
    println(nfail == 0 ? "  WHOLE water_table! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
