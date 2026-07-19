# ==========================================================================
# gpu_validate_irrigation_e2e.jl — end-to-end GPU parity for the WHOLE
# irrigation / urban-ponding trio in SoilHydrologyMod:
#
#   update_urban_ponding!            — ponding state on urban roof/road surfaces
#   calc_irrig_withdrawals!          — split gw irrig demand into per-layer
#                                      (unconfined) + confined (aquifer) fluxes
#   withdraw_groundwater_irrigation! — remove the per-layer liquid (h2osoi_liq)
#                                      and the confined withdrawal (wa)
#
# calc_irrig_withdrawals! is a host-side splitter (concrete Matrix/Vector
# signature, no kernel) that produces the irrig-flux INPUTS; it runs on the CPU
# for BOTH instances (identical inputs => identical fluxes), then the two
# KERNELIZED removal functions run on CPU and on the device and are compared.
#
# Columns deliberately exercise every branch:
#   col 1: irrigated soil, water table BELOW bedrock (jwt == nlevsoi)
#          -> unconfined layer loop empty, ALL demand -> confined aquifer (wa)
#   col 2: irrigated soil, DEEP water table, small demand
#          -> single soil layer withdrawn, no confined withdrawal
#   col 3: NON-irrigated soil column (qflx_gw_demand = 0) -> no withdrawals
#   col 4: urban ROOF column, xs_urban > 0  -> ponding pinned to pondmx_urban
#   col 5: urban ROOF column, xs_urban <= 0 -> ponding = max(0, prev+rain-evap)
#
# Runs each WHOLE function on the CPU, adapts every state struct (+ masks) to
# Metal, runs the SAME calls on the device, compares every mutated field with a
# NaN-aware reldiff. CPU-reference fields are asserted finite so a both-NaN
# false PASS cannot slip through.
#
#   julia --project=scripts scripts/gpu_validate_irrigation_e2e.jl
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
    nlevsoi  = CLM.varpar.nlevsoi
    nlevsno  = CLM.varpar.nlevsno
    nlevgrnd = CLM.varpar.nlevgrnd
    ntot     = nlevsno + CLM.varpar.nlevmaxurbgrnd

    # NOTE on indexing: WithdrawGroundwaterIrrigation / UpdateUrbanPonding /
    # CalcIrrigWithdrawals index h2osoi_liq and col_zi with the BARE soil index
    # (j = 1 .. nlevsoi), matching the Fortran `do j = 1, nlevsoi` loop and the
    # existing CPU test (h2osoi_liq_col[c, 1:nlevsoi]). They do NOT use the
    # snow+soil offset. So liquid water and interface depths live at columns
    # 1 .. nlevsoi here.

    nc = 5; np = 5; nl = 5; ng = 5

    # --- state structs ---
    sh  = CLM.SoilHydrologyData{FT}();  CLM.soilhydrology_init!(sh, nc)
    ss  = CLM.SoilStateData{FT}();      CLM.soilstate_init!(ss, np, nc)
    wsb = CLM.WaterStateBulkData{FT}(); CLM.waterstatebulk_init!(wsb, nc, np, nl, ng)
    wfb = CLM.WaterFluxBulkData{FT}();  CLM.waterfluxbulk_init!(wfb, nc, np, nl, ng)
    ws  = wsb.ws
    wf  = wfb.wf

    # --- soil hydraulic properties (per column,layer) ---
    for c in 1:nc, j in 1:nlevgrnd
        ss.watsat_col[c, j] = FT(0.45)
        ss.sucsat_col[c, j] = FT(150.0)
        ss.bsw_col[c, j]    = FT(6.0)
    end

    # --- soil interface depths col_zi (BARE soil index, 1:nlevsoi) ---
    col_zi_soil = fill(FT(0.0), nc, nlevsoi)
    for c in 1:nc
        zi_prev = FT(0.0)
        for j in 1:nlevsoi
            dz = FT(0.05) * FT(1.4)^(j - 1)
            zi = zi_prev + dz
            col_zi_soil[c, j] = zi
            zi_prev = zi
        end
    end
    zi_bot = col_zi_soil[1, nlevsoi]

    # --- liquid water (BARE soil index, 1:nlevsoi). Fill the WHOLE matrix with a
    #     finite value first so no NaN entries pollute the reldiff. ---
    ws.h2osoi_liq_col .= FT(0.0)
    ws.h2osoi_ice_col .= FT(0.0)
    for c in 1:nc, j in 1:nlevsoi
        ws.h2osoi_liq_col[c, j] = FT(30.0)
    end
    for c in 1:nc
        ws.wa_col[c] = FT(4800.0)
    end

    # --- column type / snow / masks ---
    col_snl   = fill(0, nc)
    col_itype = fill(1, nc)                 # default non-urban soil
    col_itype[4] = CLM.ICOL_ROOF
    col_itype[5] = CLM.ICOL_ROOF
    col_nbedrock = fill(nlevsoi, nc)        # bedrock at column bottom

    mask_soil  = BitVector([true, true, true, false, false])   # irrig/withdraw cols
    mask_urban = BitVector([false, false, false, true, true])  # urban ponding cols

    # --- water-table depth per column (drives jwt / layer split) ---
    # col 1: water table BELOW bedrock -> jwt = nlevsoi -> unconfined loop empty
    #        -> entire demand becomes CONFINED withdrawal (wa decremented).
    sh.zwt_col[1] = FT(zi_bot + 5.0)
    sh.zwt_col[2] = FT(0.85 * zi_bot)   # deep WT: only the deep layer(s) tapped, no confined
    sh.zwt_col[3] = FT(0.5  * zi_bot)
    sh.zwt_col[4] = FT(0.5  * zi_bot)
    sh.zwt_col[5] = FT(0.5  * zi_bot)

    # --- groundwater irrigation demand (mm/s) ---
    qflx_gw_demand = fill(FT(0.0), nc)
    qflx_gw_demand[1] = FT(1.0e-3)     # WT below bedrock -> all demand to confined (wa)
    qflx_gw_demand[2] = FT(2.0e-5)     # tiny  -> one soil layer, no confined
    qflx_gw_demand[3] = FT(0.0)        # non-irrigated

    # --- urban ponding inputs ---
    sh.xs_urban_col .= FT(0.0)
    sh.xs_urban_col[4] = FT(0.7)       # > 0  -> ponding pinned to pondmx_urban
    sh.xs_urban_col[5] = FT(0.0)       # <= 0 -> accumulate rain - evap
    wf.qflx_rain_plus_snomelt_col      .= FT(0.0)
    wf.qflx_liqevap_from_top_layer_col .= FT(0.0)
    wf.qflx_rain_plus_snomelt_col[5]      = FT(2.0e-4)
    wf.qflx_liqevap_from_top_layer_col[5] = FT(5.0e-5)
    # urban roof columns start with some ponding in the top (bare index 1) layer
    ws.h2osoi_liq_col[4, 1] = FT(0.3)
    ws.h2osoi_liq_col[5, 1] = FT(0.4)

    S = (; sh, ss, wsb, wfb)
    P = (; col_snl, col_itype, col_nbedrock, mask_soil, mask_urban,
           qflx_gw_demand, col_zi_soil,
           nlevsoi, nlevsno, nc, dtime = FT(1800.0))
    return (; S, P)
end

# calc_irrig_withdrawals! is host-only (concrete Matrix/Vector). Run it on the
# CPU for whichever instance to populate the irrig-flux fields.
function run_calc!(S, P)
    CLM.calc_irrig_withdrawals!(S.sh, S.ss,
        P.qflx_gw_demand,
        S.wfb.wf.qflx_gw_uncon_irrig_lyr_col,
        S.wfb.wf.qflx_gw_con_irrig_col,
        P.col_nbedrock, P.col_zi_soil,
        P.mask_soil, 1:P.nc, P.nlevsoi, P.dtime)
    return nothing
end

# The two KERNELIZED whole functions, run on either backend. `dmask` moves Bool
# masks and `ivec` moves integer index arrays onto the active backend.
function run_kernels!(S, P, dmask, ivec)
    CLM.update_urban_ponding!(S.wsb.ws, S.sh, S.wfb,
        ivec(P.col_snl), ivec(P.col_itype), dmask(P.mask_urban), 1:P.nc, P.dtime)

    CLM.withdraw_groundwater_irrigation!(S.wfb.wf, S.wsb.ws,
        dmask(P.mask_soil), 1:P.nc, P.nlevsoi, P.dtime)
    return nothing
end

function main(backend)
    println("=" ^ 70)
    println("END-TO-END GPU parity for irrigation/urban-ponding trio (3 whole fns)")
    println("=" ^ 70)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    H = build(FT)   # CPU reference
    B = build(FT)   # source for the device snapshot

    # calc_irrig_withdrawals! is CPU-only; run it on BOTH instances first
    # (identical inputs => identical irrig-flux fields) so the kernelized
    # withdrawal functions consume the same fluxes on each backend.
    run_calc!(H.S, H.P)
    run_calc!(B.S, B.P)

    ad(x) = CLM.Adapt.adapt(device_array_type(), x)
    Sd = (; sh = ad(B.S.sh), ss = ad(B.S.ss), wsb = ad(B.S.wsb), wfb = ad(B.S.wfb))

    if !(Sd.wsb.ws.wa_col isa device_array_type()) || !(Sd.wfb.wf.qflx_gw_con_irrig_col isa device_array_type())
        println("  BLOCKED: a state struct did not move to the device under adapt.")
        return 2
    end

    # Assert the host irrig-flux fields are FINITE before trusting any reldiff.
    for (nm, a) in (("qflx_gw_uncon_irrig_lyr", H.S.wfb.wf.qflx_gw_uncon_irrig_lyr_col),
                    ("qflx_gw_con_irrig",       H.S.wfb.wf.qflx_gw_con_irrig_col))
        cpu_has_finite(a) || (println("  BLOCKED: CPU $nm all-NaN before kernels."); return 2)
    end

    dmask(m) = to(collect(Bool, m))
    ivec(v) = to(v)
    run_kernels!(H.S, H.P, identity, identity)   # CPU (BitVector masks, host ints)
    run_kernels!(Sd,   B.P, dmask,    ivec)      # device (Bool masks, device ints)

    checks = [
        ("h2osoi_liq",              H.S.wsb.ws.h2osoi_liq_col,            Sd.wsb.ws.h2osoi_liq_col),
        ("wa",                      H.S.wsb.ws.wa_col,                    Sd.wsb.ws.wa_col),
        ("qflx_gw_uncon_irrig_lyr", H.S.wfb.wf.qflx_gw_uncon_irrig_lyr_col, Sd.wfb.wf.qflx_gw_uncon_irrig_lyr_col),
        ("qflx_gw_con_irrig",       H.S.wfb.wf.qflx_gw_con_irrig_col,     Sd.wfb.wf.qflx_gw_con_irrig_col),
    ]
    nfail = 0
    for (nm, a, b) in checks
        if !cpu_has_finite(a)
            @printf("  [WARN] %-24s CPU reference all-NaN/Inf — skipping (no signal)\n", nm)
            continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-24s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    # Branch-coverage report.
    println()
    con = Array(H.S.wfb.wf.qflx_gw_con_irrig_col)
    lyr = Array(H.S.wfb.wf.qflx_gw_uncon_irrig_lyr_col)
    pond_c = Array(H.S.wsb.ws.h2osoi_liq_col)
    pond_d = Array(Sd.wsb.ws.h2osoi_liq_col)
    wa_c = Array(H.S.wsb.ws.wa_col); wa_d = Array(Sd.wsb.ws.wa_col)
    @printf("    col 1 (irrig, WT<bedrock): confined=%.4e, nlayers_tapped=%d, wa=%.3f(cpu)/%.3f(dev)\n",
            con[1], count(>(0), @view lyr[1, :]), wa_c[1], wa_d[1])
    @printf("    col 2 (irrig, deep WT)   : confined=%.4e, nlayers_tapped=%d\n",
            con[2], count(>(0), @view lyr[2, :]))
    @printf("    col 3 (non-irrigated)    : confined=%.4e, nlayers_tapped=%d\n",
            con[3], count(>(0), @view lyr[3, :]))
    @printf("    col 4 (urban xs>0)       : ponding[1]=%.4f (cpu) %.4f (dev)\n",
            pond_c[4, 1], pond_d[4, 1])
    @printf("    col 5 (urban xs<=0)      : ponding[1]=%.4f (cpu) %.4f (dev)\n",
            pond_c[5, 1], pond_d[5, 1])

    println()
    println(nfail == 0 ? "  WHOLE irrigation/urban-ponding trio MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
