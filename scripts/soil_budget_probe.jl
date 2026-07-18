# ==========================================================================
# soil_budget_probe.jl — Julia-side TOP-SOIL-LAYER water budget on the Bow
# FATES parity case, to isolate the process behind the top-layer over-drying
# documented in PR #247 (h2ovol1 J 0.11–0.14 vs F 0.32–0.34, persistent from
# day 1). Reuses the exact build()/step loop of scripts/fates_fortran_parity.jl
# so the run is identical to the FATES parity harness; adds a per-step column-1
# water-budget dump.
#
#   julia +1.12 --project=. scripts/soil_budget_probe.jl [nsteps]
# writes scripts/soil_budget_julia.txt
# ==========================================================================
include(joinpath(@__DIR__, "fates_fortran_parity.jl"))

using Printf, Dates, NCDatasets

const NBP = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 480
const OUTP = get(ENV, "SOIL_BUDGET_OUT", joinpath(@__DIR__, "soil_budget_julia.txt"))

# local build allowing use_bedrock override (to test the harness config mismatch)
function build_ub(use_bedrock::Bool)
    inst, bounds, filt, tm = _C.clm_initialize!(; fsurdat=FSURDAT, paramfile=FPARAM,
        use_fates=true, start_date=START, dtime=Int(DTIME), use_bedrock=use_bedrock)
    for g in 1:bounds.endg
        isfinite(inst.gridcell.latdeg[g]) || (inst.gridcell.latdeg[g] = inst.gridcell.lat[g]*180/π)
        isfinite(inst.gridcell.londeg[g]) || (inst.gridcell.londeg[g] = inst.gridcell.lon[g]*180/π)
    end
    let col = inst.column, lun = inst.landunit
        for c in 1:bounds.endc
            if lun.itype[col.landunit[c]] == _C.ISTSOIL
                col.is_fates[c] = true
            end
        end
    end
    config = _C.CLMDriverConfig(use_fates=true)
    return inst, bounds, filt, _C.clump_filter_inactive_and_active, config, tm
end

function budget_run()
    # USE_BEDROCK=true  reproduces the PR #247 over-drying (shallow Bow bedrock →
    #                   zwt clamped 8.6→2.28 m → ~10 mm/day spurious drainage);
    # USE_BEDROCK=false MATCHES the Fortran case (use_bedrock=.false.) → zwt 8.6 m,
    #                   qflx_drain≡0, top layer tracks the Fortran ground truth.
    ub = get(ENV, "USE_BEDROCK", "false") == "true"
    @printf("USE_BEDROCK = %s\n", ub)
    inst, bounds, filt, filt_ia, config, _tm = build_ub(ub)
    ng, nc = bounds.endg, bounds.endc

    # the FATES soil column (first is_fates ISTSOIL col) — the one FATES site 1 maps to
    col = inst.column; lun = inst.landunit
    cfates = 0
    for c in 1:nc
        if lun.itype[col.landunit[c]] == _C.ISTSOIL && col.is_fates[c]
            cfates = c; break
        end
    end
    cfates == 0 && error("no FATES soil column found")
    joff = _C.varpar.nlevsno   # snow/soil combined index offset

    ss  = inst.soilstate
    wsb = inst.water.waterstatebulk_inst
    wfb = inst.water.waterfluxbulk_inst
    temp = inst.temperature

    fr = _C.ForcingReader()
    _C.forcing_reader_init!(fr, joinpath(FORCDIR, "clmforc.$(year(START)).nc"))
    let tf = joinpath(FORCDIR, "topo_forcing.nc")
        if isfile(tf)
            ds = NCDataset(tf, "r")
            if haskey(ds, "TOPO")
                ft = Float64(ds["TOPO"][1])
                for gg in 1:ng; inst.atm2lnd.forc_topo_grc[gg] = ft; end
                for cc in 1:nc; inst.topo.topo_col[cc] = ft; end
            end
            close(ds)
        end
    end

    io = open(OUTP, "w")
    sh = inst.soilhydrology
    cols = ["nstep","tod","day",
            "hliq1","hliq2","hliq3","hice1","hvol1","hvol2","hvol3",
            "t_soi1","watsat1","watfc1","effporo1","soilresis","dsl","soilbeta",
            "q_ev_soil","q_liqevap_top","q_infl","q_in_soil","q_surf","q_drain",
            "q_tran_veg","q_root1","q_root2","q_root3",
            "zwt","zwt_perched","wa","q_rsub_sat","q_drain_perched","hkdepth","topo_slope","frost_table"]
    println(io, "# soil budget julia  col=", cfates, " joff=", joff)
    println(io, join(cols, " "))

    g(v, i) = (v isa AbstractVector && 1 <= i <= length(v)) ? Float64(v[i]) : NaN
    m(v, i, j) = (v isa AbstractMatrix && 1 <= i <= size(v,1) && 1 <= j <= size(v,2)) ? Float64(v[i,j]) : NaN

    for i in 0:NBP
        cur = START + Second(round(Int, max(i - 1, 0) * DTIME))
        tod = hour(cur)*3600 + minute(cur)*60 + second(cur)
        calday = dayofyear(cur) + tod / _C.SECSPDAY
        (declin, _e1)   = _C.compute_orbital(calday)
        (declinm1, _e2) = _C.compute_orbital(calday - DTIME / _C.SECSPDAY)
        _C.init_daylength!(inst.gridcell, declin, declinm1, _C.ORB_OBLIQR_DEFAULT, 1:ng)
        nextsw = calday + DTIME / _C.SECSPDAY

        _C.read_forcing_step!(fr, inst.atm2lnd, cur, ng, nc)
        _C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)

        _C.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            _C.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=i, is_first_step=(i <= 1),
            is_beg_curr_day=(i > 0 && tod == 0), is_end_curr_day=false, is_beg_curr_year=false,
            dtime=DTIME, mon=month(cur), day=day(cur), secs=tod,
            jday=dayofyear(cur), photosyns=inst.photosyns)

        c = cfates
        vals = [i, tod, dayofyear(cur),
            m(wsb.ws.h2osoi_liq_col, c, 1+joff), m(wsb.ws.h2osoi_liq_col, c, 2+joff), m(wsb.ws.h2osoi_liq_col, c, 3+joff),
            m(wsb.ws.h2osoi_ice_col, c, 1+joff),
            m(wsb.ws.h2osoi_vol_col, c, 1), m(wsb.ws.h2osoi_vol_col, c, 2), m(wsb.ws.h2osoi_vol_col, c, 3),
            m(temp.t_soisno_col, c, 1+joff), m(ss.watsat_col, c, 1), m(ss.watfc_col, c, 1), m(ss.eff_porosity_col, c, 1),
            g(ss.soilresis_col, c), g(ss.dsl_col, c), g(ss.soilbeta_col, c),
            g(wfb.qflx_ev_soil_col, c), g(wfb.wf.qflx_liqevap_from_top_layer_col, c),
            g(wfb.wf.qflx_infl_col, c), g(wfb.qflx_in_soil_col, c), g(wfb.wf.qflx_surf_col, c), g(wfb.wf.qflx_drain_col, c),
            g(wfb.wf.qflx_tran_veg_col, c),
            m(wfb.qflx_rootsoi_col, c, 1), m(wfb.qflx_rootsoi_col, c, 2), m(wfb.qflx_rootsoi_col, c, 3),
            g(sh.zwt_col, c), g(sh.zwt_perched_col, c), g(wsb.ws.wa_col, c),
            g(wfb.wf.qflx_rsub_sat_col, c), g(wfb.wf.qflx_drain_perched_col, c),
            g(sh.hkdepth_col, c), g(col.topo_slope, c), g(sh.frost_table_col, c)]
        println(io, join(vals, " "))
    end
    close(io)
    _C.forcing_reader_close!(fr)
    @printf("wrote %s (%d steps, col=%d)\n", OUTP, NBP, cfates)
end

budget_run()
