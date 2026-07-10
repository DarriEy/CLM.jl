#!/usr/bin/env julia
# ==========================================================================
# fates_dieback_probe.jl — diagnose the yr2→3 stand-clearing die-back seen in the
# 10-year fates_longhorizon spin-up (carbon 10,627 → ~900, elai 3.1 → 0.11).
#
# Reuses build() + the forcing loop from fates_longhorizon.jl and, each day in the
# crash window, decomposes the mortality by FATES mode (background / carbon-starvation
# / hydraulic-failure / freezing / senescence / age-senescence / damage) + tracks
# disturbance rate, leaf carbon, and the largest cohort — to tell whether the crash is
# physical self-thinning (carbon-starvation dominated) or an over-aggressive mode/bug.
#
#   FATES_NDAYS=800 julia +1.12 --project=. scripts/fates_dieback_probe.jl
# ==========================================================================
using CLM, Printf, Dates
const _C = CLM
include(joinpath(@__DIR__, "fates_longhorizon.jl"))   # build(), helpers, forcing (guarded exit)

_c12 = _C.carbon12_element
_organs = (_C.leaf_organ, _C.fnrt_organ, _C.sapw_organ, _C.store_organ, _C.struct_organ)
coh_carbon(cc) = sum(_C.GetState(cc.prt, o, _c12) for o in _organs)   # kgC/indiv

# Per-day mortality "flux" by mode = Σ_cohorts n·(mort_rate/yr)/365·cohort_carbon.
# Absolute units are arbitrary-but-consistent; the MODE RATIO + time pattern are the signal.
function mort_decomp(site)
    m = zeros(7)  # b, c, h, fr, s, as, dg
    leafC = 0.0; storeC = 0.0; ncoh = 0
    cp = site.oldest_patch
    while cp !== nothing
        cc = cp.tallest
        while cc !== nothing
            ncoh += 1
            cm = coh_carbon(cc); w = cc.n / 365.0 * cm
            for (k, v) in enumerate((cc.bmort, cc.cmort, cc.hmort, cc.frmort, cc.smort, cc.asmort, cc.dgmort))
                isfinite(v) && (m[k] += v * w)
            end
            leafC  += _C.GetState(cc.prt, _C.leaf_organ,  _c12) * cc.n
            storeC += _C.GetState(cc.prt, _C.store_organ, _c12) * cc.n
            cc = cc.shorter
        end
        cp = cp.younger
    end
    dist = try; sum(site.disturbance_rates); catch; NaN; end
    return m, leafC, storeC, dist, ncoh
end

# Site phenology status — is the leaf-off cold-driven (cstatus) or drought (dstatus)?
function phen_state(site)
    t24 = 0.0; a = 0.0; cp = site.oldest_patch
    while cp !== nothing
        t24 += _C.GetMean(cp.tveg24) * cp.area; a += cp.area; cp = cp.younger
    end
    tC = a > 0 ? t24/a - 273.15 : NaN
    ds = try; join(string.(Int.(site.dstatus[1:min(8,length(site.dstatus))])), ""); catch; "?"; end
    return (cstat=Int(site.cstatus), dstat=ds, tC=tC,
            nchill=try Int(site.nchilldays) catch; -1 end,
            gdd=try site.grow_deg_days catch; NaN end)
end

# Largest-dbh cohort state (connects to the maxdbh-freeze lead).
function biggest(site)
    md = 0.0; cc_big = nothing; cp = site.oldest_patch
    while cp !== nothing
        cc = cp.tallest
        while cc !== nothing
            (isfinite(cc.dbh) && cc.dbh > md) && (md = cc.dbh; cc_big = cc)
            cc = cc.shorter
        end
        cp = cp.younger
    end
    cc_big === nothing && return (dbh=0.0, n=0.0, status=0, cmort=NaN, npp=NaN, leaf=NaN)
    return (dbh=cc_big.dbh, n=cc_big.n, status=cc_big.status_coh, cmort=cc_big.cmort,
            npp=cc_big.npp_acc_hold, leaf=_C.GetState(cc_big.prt, _C.leaf_organ, _c12))
end

function main()
    ndays = parse(Int, get(ENV, "FATES_NDAYS", "800"))
    lo = parse(Int, get(ENV, "CRASH_LO", "690")); hi = parse(Int, get(ENV, "CRASH_HI", "770"))
    println("="^96); println("  FATES die-back probe — mortality-mode decomposition through the yr2→3 crash"); println("="^96)
    inst, fates, config, bounds, filt, filt_ia = build()
    site = fates.sites[1]; photosyns = inst.photosyns
    dtime = 1800.0; steps_per_day = Int(round(86400/dtime)); nsteps = steps_per_day*ndays
    forcing = get(ENV, "FATES_FORCING", "$DATA/domain_Aripuana_Amazon/data/forcing/CLM_input/clmforc.2004.nc")
    fyr = 2004
    fr = _C.ForcingReader(); _C.forcing_reader_init!(fr, forcing); fr.interp_time = true
    start_date = DateTime(fyr,1,1)
    @printf("  cold start carbon=%.4g  (crash window days %d–%d; total %d)\n", total_site_carbon(site), lo, hi, ndays)
    @printf("  %5s %10s %7s %9s %9s | %8s %8s %8s %8s %8s %8s | %8s | %6s %7s %6s\n",
        "day","carbon","elai","leafC","storeC","bmort","cmort","hmort","frmort","smort/as","dgmort","disturb","bigDBH","bigN","bigSt")
    println("  " * "-"^94)
    daycount = 0
    for i in 1:nsteps
        step_start = start_date + Second((i-1)*Int(dtime))
        is_beg = (Dates.hour(step_start)==0 && Dates.minute(step_start)==0 && Dates.second(step_start)==0)
        yr_off = Dates.year(step_start) - fyr
        forcing_time = yr_off == 0 ? step_start : step_start - Dates.Year(yr_off)
        _C.read_forcing_step!(fr, inst.atm2lnd, forcing_time, 1, 1;
            gridcell_latdeg=inst.gridcell.latdeg, gridcell_londeg=inst.gridcell.londeg, dtime=Int(dtime))
        inst.atm2lnd.forc_topo_grc[1]=200.0; inst.topo.topo_col[1]=200.0
        _C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        sod = Dates.hour(step_start)*3600 + Dates.minute(step_start)*60
        calday = Dates.dayofyear(step_start) + sod/86400.0
        (declin,_e) = _C.compute_orbital(calday); nextsw_cday = calday + Int(dtime)/86400.0
        try
            _C.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin, 0.4091,
                false, false, "20260101", false; nstep=i, is_first_step=(i==1), is_beg_curr_day=is_beg,
                is_end_curr_day=false, is_beg_curr_year=false, dtime=dtime, mon=Dates.month(step_start),
                day=Dates.day(step_start), secs=sod, jday=Dates.dayofyear(step_start), photosyns=photosyns)
        catch e
            @printf("  ✗ ERROR day %d: %s\n", daycount, first(split(sprint(showerror,e),"\n"))); break
        end
        if is_beg && i>1
            daycount = (i-1) ÷ steps_per_day
            # print in the crash window, else every 60 days for context
            if lo <= daycount <= hi || daycount % 30 == 0
                m, leafC, storeC, dist, _ = mort_decomp(site)
                c = total_site_carbon(site); el = inst.canopystate.elai_patch[2]; ph = phen_state(site)
                tvraw = try; Float64(inst.fates.bc_in[1].t_veg_pa[1]) - 273.15; catch; NaN; end
                tvcanopy = try; Float64(inst.temperature.t_veg_patch[2]) - 273.15; catch; NaN; end
                @printf("  %5d %9.4g %6.3f | t_veg_pa=%6.1f t_veg_patch=%6.1f tveg24=%5.1f | cstat=%d nchill=%3d | cmort=%6.3g\n",
                    daycount, c, el, tvraw, tvcanopy, ph.tC, ph.cstat, ph.nchill, m[2])
            end
        end
    end
    _C.forcing_reader_close!(fr)
    println("\n  legend: bmort=background cmort=carbon-starvation hmort=hydraulic frmort=freezing smort/as=senescence dgmort=damage")
    return 0
end
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
