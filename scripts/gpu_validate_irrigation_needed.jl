# ==========================================================================
# gpu_validate_irrigation_needed.jl — host-vs-Metal parity for the kernelized
# IrrigationMod science functions (feature-gated: crop + irrigate):
#   calc_irrigation_needed!  (patch-check + layer accumulation + deficit + rate,
#                             incl. the volr-limit path: c2g_irrig!/irrig_volr_ratio!/
#                             irrig_deficit_limited!)
#   calc_irrigation_fluxes!  (calc_bulk_withdrawals! patch kernel + p2c_irrig! +
#                             calc_application_fluxes! drip/sprinkler)
#   calc_total_gw_uncon_irrig!  (per-column layer sum — tested standalone)
#
# Each whole function runs on the CPU (Float64) and on Metal (Float32); every
# mutated field is compared with a NaN-aware reldiff.
#
#   julia --project=scripts scripts/gpu_validate_irrigation_needed.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
mf(x) = mf(Metal.MtlArray, x)

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0; n = 0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + abs(Float64(A[i])))); n += 1
    end
    return m, n
end

# Build an irrigation fixture (nc columns / np patches) that produces a real deficit.
# `rof` toggles the volr river-volume-limit path.
function build(::Type{FT}; rof::Bool) where {FT}
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5); CLM.varcon_init!()
    nc, np, nlevsoi, ng = 2, 2, 5, 2
    irrig = CLM.IrrigationData{FT}()
    CLM.irrigation_init_allocate!(irrig, np, nc, nlevsoi)
    irrig.dtime = 1800
    irrig.irrig_nsteps_per_day = CLM.calc_irrig_nsteps_per_day(irrig.params.irrig_length, 1800)
    irrig.relsat_wilting_point_col .= FT(0.1)
    irrig.relsat_target_col .= FT(0.9)
    irrig.irrig_method_patch .= CLM.IRRIG_METHOD_DRIP
    irrig.params.limit_irrigation_if_rof_enabled = rof

    col = CLM.ColumnData{FT}()
    col.landunit = collect(1:nc); col.gridcell = collect(1:nc); col.wtgcell = ones(FT, nc)
    col.z = fill(FT(0.05), nc, nlevsoi); col.dz = fill(FT(0.1), nc, nlevsoi)
    col.nbedrock = fill(nlevsoi, nc)

    grc = CLM.GridcellData{FT}(); grc.londeg = zeros(FT, ng); grc.area = fill(FT(1.0e4), ng)

    pch = CLM.PatchData{FT}()
    pch.column = collect(1:np); pch.gridcell = collect(1:np); pch.landunit = collect(1:np)
    pch.itype = fill(CLM.nc3irrig, np); pch.wtcol = ones(FT, np); pch.active = trues(np)

    npft = CLM.MXPFT + 1
    pftcon_irrigated = zeros(FT, npft); pftcon_irrigated[CLM.nc3irrig] = FT(1.0)

    elai = fill(FT(2.0), np)
    t_soisno = fill(FT(CLM.TFRZ + 10.0), nc, nlevsoi)
    eff_porosity = fill(FT(0.4), nc, nlevsoi)
    h2osoi_liq = fill(FT(0.05 * CLM.DENH2O * 0.1), nc, nlevsoi)  # dry → deficit
    volr = fill(FT(50.0), ng)                                     # some river water
    lts = fill(irrig.params.irrig_start_time - irrig.dtime + 1, np)
    mask = trues(np)
    return (; irrig, col, grc, pch, pftcon_irrigated, elai, t_soisno, eff_porosity,
            h2osoi_liq, volr, lts, mask, nc, np, nlevsoi, ng)
end

# ---- run calc_irrigation_needed! then calc_irrigation_fluxes!, on a given backend --
function run_pipeline!(f, adapt; rof::Bool)
    irr  = adapt(f.irrig); col = adapt(f.col); grc = adapt(f.grc); pch = adapt(f.pch)
    pfti = adapt(f.pftcon_irrigated); elai = adapt(f.elai); tsn = adapt(f.t_soisno)
    effp = adapt(f.eff_porosity); liq = adapt(f.h2osoi_liq); volr = adapt(f.volr)
    lts  = adapt(f.lts); mask = adapt(f.mask)
    CLM.calc_irrigation_needed!(irr, elai, tsn, effp, liq, volr, rof, pfti, lts,
        col, grc, pch, mask, 1:f.nc, 1:f.np, 1:f.ng, f.nlevsoi)
    ss = CLM.SoilStateData(); CLM.soilstate_init!(ss, f.np, f.nc)
    sh = CLM.SoilHydrologyData(); CLM.soilhydrology_init!(sh, f.nc)
    wfb = CLM.WaterFluxBulkData(); CLM.waterfluxbulk_init!(wfb, f.nc, f.np, 1, 1)
    wfb.wf.qflx_sfc_irrig_col .= 0
    ss = adapt(ss); sh = adapt(sh); wfb = adapt(wfb)
    CLM.calc_irrigation_fluxes!(irr, sh, ss, wfb, col, pch, adapt(f.mask), adapt(f.mask),
        1:f.nc, 1:f.np, f.nlevsoi, 1800.0)
    return (; irr, wf = wfb.wf)
end

function main(backend)
    println("="^64); println("Irrigation (IrrigationMod science) — host vs Metal"); println("="^64)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, _ = backend; @printf("  Backend: %s\n\n", name)
    id(x) = x
    nfail = 0
    for rof in (false, true)
        @printf("  --- scenario: rof=%s (volr-limit %s) ---\n", rof, rof ? "ON" : "off")
        h = run_pipeline!(build(Float64; rof=rof), id;  rof=rof)
        d = run_pipeline!(build(Float32; rof=rof), mf;  rof=rof)
        checks = [("sfc_irrig_rate",   h.irr.sfc_irrig_rate_patch,   d.irr.sfc_irrig_rate_patch),
                  ("irrig_rate_demand",h.irr.irrig_rate_demand_patch,d.irr.irrig_rate_demand_patch),
                  ("n_irrig_steps",    h.irr.n_irrig_steps_left_patch,d.irr.n_irrig_steps_left_patch),
                  ("qflx_sfc_irrig_col",h.wf.qflx_sfc_irrig_col,     d.wf.qflx_sfc_irrig_col),
                  ("qflx_irrig_drip",  h.wf.qflx_irrig_drip_patch,   d.wf.qflx_irrig_drip_patch)]
        for (nm, a, b) in checks
            r, n = reldiff(a, b); ok = r < 1f-3
            @printf("    [%s] %-20s rel=%.2e over %d\n", ok ? "PASS" : "FAIL", nm, r, n)
            ok || (nfail += 1)
        end
    end

    # calc_total_gw_uncon_irrig! standalone (per-column layer sum)
    println("  --- calc_total_gw_uncon_irrig! (layer sum) ---")
    function run_gw(::Type{FT}, adapt) where {FT}
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        nc, nlevsoi = 3, 5
        wf = CLM.WaterFluxData{FT}()
        wf.qflx_gw_uncon_irrig_lyr_col = FT[FT(0.1)*(c+j) for c in 1:nc, j in 1:nlevsoi]
        wf.qflx_gw_uncon_irrig_col = fill(FT(NaN), nc)
        wf = adapt(wf); mask = adapt(trues(nc))
        CLM.calc_total_gw_uncon_irrig!(wf, mask, 1:nc, nlevsoi)
        return wf.qflx_gw_uncon_irrig_col
    end
    r, n = reldiff(run_gw(Float64, id), run_gw(Float32, mf))
    ok = r < 1f-3; @printf("    [%s] qflx_gw_uncon_irrig_col rel=%.2e over %d\n", ok ? "PASS" : "FAIL", r, n); ok || (nfail += 1)

    println()
    println(nfail == 0 ? "  IRRIGATION science kernels MATCH host on $name" : "  DIVERGENCE ($nfail)")
    return nfail == 0 ? 0 : 1
end
exit(main(detect_backend()))
