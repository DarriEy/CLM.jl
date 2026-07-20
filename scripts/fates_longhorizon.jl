#!/usr/bin/env julia
# ==========================================================================
# fates_longhorizon.jl — long-horizon FATES stability + demography validation.
#
# The deepest existing FATES live test (test_fates_spinup.jl) runs 15 days. This
# runs a carbon-only FATES clm_drv! loop for a FULL YEAR+ (FATES_NDAYS, default
# 365 @ dtime=1800s), sampling the demographic trajectory daily, to validate that
# over demographic timescales the model:
#   - never errors / goes non-finite,
#   - conserves carbon EVERY day (daily TotalBalanceCheck holds — a leak accumulates
#     over a year and would show up),
#   - the census stays physical every day (n>=0, dbh>0, patch areas sum to AREA),
#   - the demography ADVANCES and stays bounded (recruitment/mortality/growth/
#     disturbance happen; biomass approaches a quasi-steady band, no blow-up/collapse).
#
#   FATES_NDAYS=365 julia +1.12 --project=. scripts/fates_longhorizon.jl
# ==========================================================================
# NB: `Base.include(@__MODULE__, ...)`, not bare `include`. Several of these
# scripts are loaded by their tests into a fresh `Module(:X)`, which does NOT
# bind a bare `include` — that form throws UndefVarError there.
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "test", "testdata.jl"))

using CLM, Printf, Dates
const _C = CLM
const DATA = symfluence_data_root()

# ---- census / carbon helpers (mirror test_fates_live_modes.jl / test_fates_spinup.jl) ----
function census(inst)
    site = inst.fates.sites[1]
    npatch = 0; ncoh = 0; totarea = 0.0; bad = String[]
    cp = site.oldest_patch
    while cp !== nothing
        npatch += 1
        isfinite(cp.area) || push!(bad, "patch.area"); totarea += cp.area
        cc = cp.tallest
        while cc !== nothing
            ncoh += 1
            (isfinite(cc.n) && cc.n >= 0.0) || push!(bad, "cohort.n")
            (isfinite(cc.dbh) && cc.dbh > 0.0) || push!(bad, "cohort.dbh")
            cc = cc.shorter
        end
        cp = cp.younger
    end
    return (npatch=npatch, ncoh=ncoh, totarea=totarea, bad=bad)
end
function total_site_carbon(site)
    tot = 0.0; cp = site.oldest_patch
    while cp !== nothing
        cc = cp.tallest
        while cc !== nothing
            for org in (_C.leaf_organ, _C.fnrt_organ, _C.sapw_organ, _C.store_organ, _C.struct_organ)
                tot += _C.GetState(cc.prt, org, _C.carbon12_element) * cc.n
            end
            cc = cc.shorter
        end
        cp = cp.younger
    end
    return tot
end
function max_dbh(site)
    md = 0.0; cp = site.oldest_patch
    while cp !== nothing
        cc = cp.tallest
        while cc !== nothing; md = max(md, cc.dbh); cc = cc.shorter; end
        cp = cp.younger
    end
    return md
end

# Per-day carbon budget + cohort leaf/stem structure (why NPP<0?). *_acc_hold hold the
# previous day's accumulated GPP/NPP/resp [kgC/indiv/day]; scale by n [/m2] for /m2/day.
function carbon_budget(site)
    gpp=0.0; npp=0.0; resp=0.0; lc=0.0; sc=0.0; stc=0.0; fc=0.0; storc=0.0
    println("    coh  pft     n/m2    dbh  status  treelai   c_area   leafC   sapwC  structC  npp/d[gC/m2]")
    ic=0; cp=site.oldest_patch
    while cp!==nothing
        cc=cp.tallest
        while cc!==nothing
            ic+=1
            gh = isfinite(cc.gpp_acc_hold) ? cc.gpp_acc_hold : 0.0
            nh = isfinite(cc.npp_acc_hold) ? cc.npp_acc_hold : 0.0
            rh = isfinite(cc.resp_acc_hold) ? cc.resp_acc_hold : 0.0
            gpp+=gh*cc.n; npp+=nh*cc.n; resp+=rh*cc.n
            l=_C.GetState(cc.prt,_C.leaf_organ,_C.carbon12_element)
            s=_C.GetState(cc.prt,_C.sapw_organ,_C.carbon12_element)
            st=_C.GetState(cc.prt,_C.struct_organ,_C.carbon12_element)
            f=_C.GetState(cc.prt,_C.fnrt_organ,_C.carbon12_element)
            sto=_C.GetState(cc.prt,_C.store_organ,_C.carbon12_element)
            lc+=l*cc.n; sc+=s*cc.n; stc+=st*cc.n; fc+=f*cc.n; storc+=sto*cc.n
            if ic<=8
                @printf("    %3d  %3d  %7.4f  %5.2f  %5d  %7.3f  %7.4f  %6.1f  %6.1f  %7.1f  %8.3f\n",
                    ic, cc.pft, cc.n, cc.dbh, cc.status_coh, cc.treelai, cc.c_area,
                    l*cc.n*1000, s*cc.n*1000, st*cc.n*1000, nh*cc.n*1000)
            end
            cc=cc.shorter
        end
        cp=cp.younger
    end
    @printf("    TOTAL/m2/day: GPP=%.4g NPP=%.4g RESP=%.4g gC/m2/day  |  pools gC/m2: leaf=%.1f sapw=%.1f struct=%.1f fnrt=%.1f store=%.1f\n",
        gpp*1000, npp*1000, resp*1000, lc*1000, sc*1000, stc*1000, fc*1000, storc*1000)
    return nothing
end

# use_bedrock. This harness INHERITED the driver default, so #252 (which made that
# default conditional -- `use_fates => false`, previously `true`) silently flipped it
# under us. At #252 that flip IS what makes multi-year runs here die at day 74 on a
# 400 W/m2 longwave imbalance (pinning it back to `true` there rescues the run) -- but
# on current main pinning no longer rescues, so a SECOND defect followed it. See
# fates_fortran_parity/README.md D4. `nothing` = inherit the (correct, CTSM-matching)
# conditional default; FATES_USE_BEDROCK=1/0 pins it, for attributing that failure.
_use_bedrock() = (v = get(ENV, "FATES_USE_BEDROCK", "");
                  v == "1" ? true : v == "0" ? false : nothing)

function build()
    fsurdat = get(ENV, "FATES_FSURDAT",
        "$DATA/domain_Aripuana_Amazon/settings/CLM/parameters/surfdata_clm.nc")
    paramfile = get(ENV, "FATES_PARAMFILE",
        "$DATA/domain_Aripuana_Amazon/settings/CLM/parameters/clm5_params.nc")
    fyr = parse(Int, get(ENV, "FATES_YEAR", "2004"))
    # Optional climate-appropriate PFT screening (fixed-biogeography). Default OFF ->
    # the all-PFT NBG cold start (the boom-bust baseline).
    #
    # STATUS (D4, 2026-07-20): #197's 4-year evidence for this screen NO LONGER
    # REPRODUCES -- both arms die before day 100 on a fatal longwave imbalance (see
    # fates_fortran_parity/README.md D4). Note also that D3/#274's suggestion that the
    # screen compensates for the running-mean bug is DISPROVED: that bug was inert at
    # this harness's dtime=1800. The screen is neither validated nor retired; treat any
    # boom-bust claim from this script as unverified until the multi-year blocker is fixed.
    #
    # FATES_BIOGEOG selects a
    # named screen: "drop_cold_deciduous" (or "1") seeds only non-cold-deciduous PFTs
    # (evergreen + drought-deciduous) — the tropical-appropriate set for Aripuana, so
    # no cold-deciduous cohort is ever seeded/recruited to drive the die-back cycle.
    biogeog_screen = :none
    let bg = get(ENV, "FATES_BIOGEOG", "")
        if bg == "1" || bg == "drop_cold_deciduous"
            biogeog_screen = :drop_cold_deciduous
        elseif bg != "" && bg != "0" && bg != "none"
            biogeog_screen = Symbol(bg)
        end
    end
    # REAL cold-start. clm_initialize! reads surfdata/params, runs the CLM init_cold
    # routines and — crucially — computes a CONSISTENT initial surface albedo, so the
    # step-1 surface-energy balance produces FINITE t_grnd/t_veg. That is the piece the
    # earlier hand-rolled scaffold could never prime: the cold-start albedo<-hydrology<-
    # energy circular dependency leaves fabd/fabi=NaN, and sabv=fabd*solad=NaN*0=NaN
    # NaNs the whole energy chain → t_veg=NaN → FATES photosynthesis endruns on
    # imaginary roots the first daylit step. With use_fates=true, clm_initialize! ALSO
    # bootstraps + cold-starts one carbon-only FATES site per CLM column onto inst.fates
    # (real FATES default PFT params). This is the only setup that hands FATES a valid
    # leaf temperature and hence non-NaN photosynthesis.
    inst, bounds, filt, _tm = _C.clm_initialize!(; fsurdat=fsurdat, paramfile=paramfile,
        use_fates=true, start_date=DateTime(fyr,1,1), dtime=1800,
        fates_biogeog_screen=biogeog_screen, use_bedrock=_use_bedrock())
    biogeog_screen == :none || @printf("  fixed-biogeog screen: %s (climate-appropriate PFTs only)\n", biogeog_screen)
    # coszen closure reads grc.lat/lon (radians); ensure latdeg/londeg exist for the
    # forcing reader's time-interpolation of FSDS (fills from lat/lon if surfdata omits).
    for g in 1:bounds.endg
        isfinite(inst.gridcell.latdeg[g]) || (inst.gridcell.latdeg[g] = inst.gridcell.lat[g]*180/π)
        isfinite(inst.gridcell.londeg[g]) || (inst.gridcell.londeg[g] = inst.gridcell.lon[g]*180/π)
    end
    # clm_initialize! attaches inst.fates + cold-starts the site cohorts, but nothing in
    # the codebase flags WHICH HLM columns are FATES. Every driver block that couples
    # FATES to the biogeophysics (radiation → sunlit PAR, photosynthesis → GPP) is gated
    # on col.is_fates[c]; without it, demography runs (via config.use_fates) but the
    # canopy never photosynthesizes (solar_zenith_flag stays false, parsun=0, GPP=0).
    # Flag the natural-vegetated soil column(s) — the s-th is_fates column maps 1:1 onto
    # FATES site s (fates_veg_patches / FatesSunShadeFracs).
    let col=inst.column, lun=inst.landunit, nflag=0
        for c in 1:bounds.endc
            if lun.itype[col.landunit[c]] == _C.ISTSOIL
                col.is_fates[c] = true; nflag += 1
            end
        end
        @printf("  flagged %d FATES column(s); fates.nsites=%d\n", nflag, inst.fates.nsites)
    end
    # The generic surfdata cold-start applies a COLD deep-soil profile (~272 K, partly
    # frozen) — unphysical for a tropical Amazon site and it freezes deep roots out of
    # the water-uptake calc. Warm the FATES column's whole soil column + ground to a
    # tropical ~299 K so the entire root zone is thawed. (opt out: FATES_WARMSOIL=0)
    if get(ENV,"FATES_WARMSOIL","1")=="1"
        col=inst.column; temp=inst.temperature; joff=_C.varpar.nlevsno; ngr=_C.varpar.nlevgrnd
        for c in 1:bounds.endc
            col.is_fates[c] || continue
            temp.t_grnd_col[c]=299.0
            for j in 1:ngr; temp.t_soisno_col[c, joff+j]=299.0; end
        end
    end
    # DIAGNOSTIC-ONLY escape hatch. On current `main` this site trips a FATAL
    # longwave energy-balance error (errlon ~ -410 W/m2 at p=7) around day 74-98,
    # so a multi-YEAR demographic trajectory cannot be measured at all. That is a
    # real, unfixed defect and it is NOT the subject of this harness. Setting
    # FATES_SOFT_BALANCE=1 degrades every hard balance error to a warning so the
    # demography can be observed past it.
    #
    # A run with this set has a KNOWN, UNCLOSED energy balance: its absolute
    # carbon numbers are NOT trustworthy. Use it only to read the SHAPE of the
    # demographic trajectory (boom-bust vs. bounded), never as a parity or
    # conservation result. Default OFF -> byte-identical, still fatal.
    if get(ENV,"FATES_SOFT_BALANCE","0")=="1"
        inst.balcheck.hard_error = false
        println("  [DIAG] FATES_SOFT_BALANCE=1: balance hard errors degraded to warnings.")
        println("         Absolute carbon is NOT trustworthy in this run — trajectory SHAPE only.")
    end
    config = _C.CLMDriverConfig(use_fates=true)
    filt_ia = _C.clump_filter_inactive_and_active
    return inst, inst.fates, config, bounds, filt, filt_ia
end

function set_forcing!(inst, sunfrac)
    a=inst.atm2lnd; gg=1; T=290.0; pbot=9.9e4; q=0.008
    th=T*(1.0e5/pbot)^0.286; rho=pbot/(287.058*T*(1.0+0.61*q)); vp=q*pbot/(0.622+0.378*q)
    a.forc_t_not_downscaled_grc[gg]=T; a.forc_th_not_downscaled_grc[gg]=th
    a.forc_pbot_not_downscaled_grc[gg]=pbot; a.forc_q_not_downscaled_grc[gg]=q
    a.forc_rho_not_downscaled_grc[gg]=rho; a.forc_lwrad_not_downscaled_grc[gg]=350.0
    a.forc_rain_not_downscaled_grc[gg]=0.0; a.forc_snow_not_downscaled_grc[gg]=0.0
    a.forc_u_grc[gg]=2.0; a.forc_v_grc[gg]=0.0; a.forc_hgt_grc[gg]=30.0; a.forc_hgt_u_grc[gg]=30.0
    a.forc_hgt_t_grc[gg]=30.0; a.forc_hgt_q_grc[gg]=30.0; a.forc_vp_grc[gg]=vp
    a.forc_pco2_grc[gg]=367.0e-6*pbot; a.forc_po2_grc[gg]=0.209*pbot; a.forc_topo_grc[gg]=150.0
    for b in 1:size(a.forc_solad_not_downscaled_grc,2)
        a.forc_solad_not_downscaled_grc[gg,b]=(b==1 ? 300.0 : 250.0)*sunfrac
        a.forc_solai_grc[gg,b]=(b==1 ? 100.0 : 80.0)*sunfrac
    end
    inst.topo.topo_col[1]=150.0
    return nothing
end

function main()
    ndays = parse(Int, get(ENV, "FATES_NDAYS", "365"))
    println("="^70); println("  FATES long-horizon stability + demography — $ndays days"); println("="^70)
    inst, fates, config, bounds, filt, filt_ia = build()
    site = fates.sites[1]; photosyns = inst.photosyns
    dtime = 1800.0; steps_per_day = Int(round(86400/dtime)); nsteps = steps_per_day*ndays
    # Real met forcing (default Aripuana tropical 2004): hourly TBOT/FSDS/PREC/…, model
    # steps half-hourly (interp). read_forcing_step! also sets forc_pco2/po2 from pbot.
    forcing = get(ENV, "FATES_FORCING", "$DATA/domain_Aripuana_Amazon/data/forcing/CLM_input/clmforc.2004.nc")
    fyr = parse(Int, get(ENV, "FATES_YEAR", "2004"))
    fr = _C.ForcingReader(); _C.forcing_reader_init!(fr, forcing); fr.interp_time = true
    start_date = DateTime(fyr, 1, 1); daycount = 0
    println("  forcing: $forcing  (lat -7.51, lon -59.63, tropical)\n")

    c0 = total_site_carbon(site); cen0 = census(inst)
    @printf("  cold start: ncoh=%d npatch=%d carbon=%.4g dbh=%.4f\n\n", cen0.ncoh, cen0.npatch, c0, max_dbh(site))
    @printf("  %6s %6s %7s %12s %9s %9s %8s %8s %8s\n", "day", "ncoh", "npatch", "carbon", "maxdbh", "elai", "GPP/d", "NPP/d", "bal")
    println("  " * "-"^62)

    no_error=true; errmsg=""; nan_days=String[]; days_advanced=0; day_bal_ok=0
    min_ncoh=cen0.ncoh; max_ncoh=cen0.ncoh; max_npatch=cen0.npatch
    min_c=c0; max_c=c0; births=0; deaths=0; prev_ncoh=cen0.ncoh
    acc_gpp=0.0; acc_npp=0.0; day_gpp=0.0; day_npp=0.0   # per-day flux [gC/m2/day]
    t0 = time()
    for i in 1:nsteps
        step_start = start_date + Second((i-1)*Int(dtime))
        is_beg = (Dates.hour(step_start)==0 && Dates.minute(step_start)==0 && Dates.second(step_start)==0)
        # Multi-year forcing wrap: the forcing file holds ONE year (fyr=2004). For
        # runs >365 days the model date rolls into fyr+1… which the file lacks, and
        # the reader (interp_time=true → _find_bracket_time) would CLAMP to the last
        # 2004 timestep — freezing stale Dec-31 forcing for the entire 2nd year on.
        # Map the model instant back onto the forcing year (same month/day/hour, so a
        # full realistic diurnal+annual cycle repeats each year). yr_off=0 in year 1
        # → forcing_time==step_start (byte-identical to the old single-year path).
        yr_off = Dates.year(step_start) - fyr
        forcing_time = yr_off == 0 ? step_start : step_start - Dates.Year(yr_off)
        _C.read_forcing_step!(fr, inst.atm2lnd, forcing_time, 1, 1;
            gridcell_latdeg=inst.gridcell.latdeg, gridcell_londeg=inst.gridcell.londeg, dtime=Int(dtime))
        inst.atm2lnd.forc_topo_grc[1]=200.0; inst.topo.topo_col[1]=200.0   # matched → no lapse correction
        _C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        sod = Dates.hour(step_start)*3600 + Dates.minute(step_start)*60 + Dates.second(step_start)
        calday = Dates.dayofyear(step_start) + sod/86400.0
        (declin, _e) = _C.compute_orbital(calday); nextsw_cday = calday + Int(dtime)/86400.0
        if get(ENV,"FATES_TVEG","")=="1" && i<=6
            a=inst.atm2lnd
            @printf("  [TVEG i=%d %02d:%02d] PRE t_grnd=%.3f t_veg[2]=%.3f t_soisno[c,6]=%.3f | forc_t=%.2f lw=%.2f solad=%.2f rho=%.4f q=%.5f pbot=%.1f\n",
                i, Dates.hour(step_start), Dates.minute(step_start),
                inst.temperature.t_grnd_col[1], inst.temperature.t_veg_patch[2],
                inst.temperature.t_soisno_col[1,6],
                a.forc_t_downscaled_col[1], a.forc_lwrad_downscaled_col[1],
                a.forc_solad_downscaled_col[1,1], a.forc_rho_downscaled_col[1],
                a.forc_q_downscaled_col[1], a.forc_pbot_downscaled_col[1])
        end
        try
            _C.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin, 0.4091,
                false, false, "20260101", false; nstep=i, is_first_step=(i==1), is_beg_curr_day=is_beg,
                is_end_curr_day=false, is_beg_curr_year=false, dtime=dtime, mon=Dates.month(step_start),
                day=Dates.day(step_start), secs=sod, jday=Dates.dayofyear(step_start), photosyns=photosyns)
        catch e
            no_error=false; errmsg=sprint(showerror, e); @printf("  ✗ ERROR at step %d (day %d): %s\n", i, daycount, first(split(errmsg,"\n"))); break
        end
        # accumulate stand GPP/NPP this step: gpp_tstep [kgC/indiv/step] × n [/m2] × 1000 → gC/m2
        let cp=site.oldest_patch
            while cp!==nothing
                cc=cp.tallest
                while cc!==nothing
                    isfinite(cc.gpp_tstep) && (acc_gpp += cc.gpp_tstep*cc.n*1000)
                    isfinite(cc.npp_tstep) && (acc_npp += cc.npp_tstep*cc.n*1000)
                    cc=cc.shorter
                end
                cp=cp.younger
            end
        end
        if get(ENV,"FATES_PROBE","")=="1" && i in (24,36)
            solad = inst.atm2lnd.forc_solad_downscaled_col[1,:]
            gmax=0.0; parmax=0.0; szf=false; btmax=0.0; vcmax=0.0; anetmax=-1e9; cp=site.oldest_patch
            while cp!==nothing
                (isdefined(cp,:solar_zenith_flag)) && (szf = szf || cp.solar_zenith_flag)
                try; parmax = max(parmax, maximum(x->isfinite(x) ? x : 0.0, cp.ed_parsun_z)); catch; end
                try; btmax = max(btmax, maximum(x->isfinite(x) ? x : 0.0, cp.btran_ft)); catch; end
                cc=cp.tallest
                while cc!==nothing
                    isfinite(cc.gpp_tstep)&&(gmax=max(gmax,cc.gpp_tstep))
                    isfinite(cc.vcmax25top)&&(vcmax=max(vcmax,cc.vcmax25top))
                    (isdefined(cc,:ts_net_uptake) && !isempty(cc.ts_net_uptake)) && (anetmax=max(anetmax, maximum(x->isfinite(x) ? x : -1e9, cc.ts_net_uptake)))
                    cc=cc.shorter
                end
                cp=cp.younger
            end
            @printf("  [PROBE i=%d %02d:00] coszen=%.3f solad=%s szen=%s parsun=%.3g btran=%.4f vc25=%.1f anet=%.4g gpp=%.5g\n",
                    i, Dates.hour(step_start), inst.surfalb.coszen_col[1], string(round.(solad,digits=1)), szf, parmax, btmax, vcmax, anetmax, gmax)
            bc=inst.fates.bc_in[1]
            @printf("        FATES soil bc: smp_sl=%s  h2o_liqvol=%s  eff_por=%s  tempk=%s\n",
                string(round.(bc.smp_sl,digits=0)), string(round.(bc.h2o_liqvol_sl,digits=3)),
                string(round.(bc.eff_porosity_sl,digits=3)), string(round.(bc.tempk_sl,digits=1)))
            ep=_C.EDPftvarcon_inst[]
            @printf("        PFT7 smpsc=%.0f smpso=%.0f\n", ep.smpsc[7], ep.smpso[7])
        end
        if is_beg && i>1
            # ELAPSED days since the start of the run, NOT day-of-year. This was
            # `Dates.dayofyear(step_start) - 1`, which RESETS every 1 January — so on
            # any run past 365 days every day-indexed number the script reports (the
            # sampled rows, `nan_days`, the "ran N/M days" summary and the s/yr-sim
            # rate) silently restarted from 0, and a year-2 failure printed as a
            # two-digit day indistinguishable from a year-1 one. Harmless while
            # multi-year FATES was impossible; actively misleading now that it isn't.
            daycount = Dates.value(Date(step_start) - Date(start_date))
            day_gpp=acc_gpp; day_npp=acc_npp; acc_gpp=0.0; acc_npp=0.0   # close out the day
            days_advanced += 1
            bal_ok = _C.TotalBalanceCheck(site, -1) === nothing
            bal_ok && (day_bal_ok += 1)
            cen = census(inst); c = total_site_carbon(site); md = max_dbh(site)
            isempty(cen.bad) || push!(nan_days, "day$daycount:$(cen.bad)")
            isfinite(c) || push!(nan_days, "day$daycount:carbon")
            dn = cen.ncoh - prev_ncoh; dn>0 && (births += dn); dn<0 && (deaths += -dn); prev_ncoh = cen.ncoh
            min_ncoh=min(min_ncoh,cen.ncoh); max_ncoh=max(max_ncoh,cen.ncoh); max_npatch=max(max_npatch,cen.npatch)
            min_c=min(min_c,c); max_c=max(max_c,c)
            if daycount % max(1, ndays ÷ 24) == 0 || daycount == ndays  # ~24 sampled rows
                @printf("  %6d %6d %7d %12.5g %9.4f %9.4f %8.3f %8.3f %8s\n", daycount, cen.ncoh, cen.npatch, c, md,
                        inst.canopystate.elai_patch[2], day_gpp, day_npp, bal_ok ? "ok" : "FAIL")
            end
        end
    end
    _C.forcing_reader_close!(fr)
    dt = time()-t0
    cenF = census(inst); cF = total_site_carbon(site)
    println("\n  " * "-"^62)
    @printf("  ran %d/%d days in %.1f s (%.2f s/yr-sim)\n", daycount, ndays, dt, dt*365/max(daycount,1))
    @printf("  carbon:  cold=%.4g  final=%.4g  min=%.4g  max=%.4g\n", c0, cF, min_c, max_c)
    @printf("  ncoh:    range [%d, %d]   npatch max=%d   births=%d deaths=%d\n", min_ncoh, max_ncoh, max_npatch, births, deaths)
    @printf("  daily balance held: %d/%d days\n", day_bal_ok, days_advanced)
    if get(ENV,"FATES_BUDGET","")=="1"
        println("\n  --- carbon budget (final day) ---"); carbon_budget(site)
    end
    # ---- validation verdict ----
    checks = [
        ("no error over horizon",              no_error),
        ("no non-finite day",                  isempty(nan_days)),
        ("carbon conserved EVERY day",         day_bal_ok == days_advanced && days_advanced >= ndays-1),
        ("cohorts persisted (ncoh>=1)",        min_ncoh >= 1),
        ("demography advanced (births+deaths)",births + deaths > 0),
        ("disturbance spawned patch (npatch>=2)", max_npatch >= 2),
        ("carbon bounded (no blow-up/collapse)", isfinite(cF) && cF > 0 && max_c < 50*max(c0,1e-6)),
        ("final balance holds",                _C.TotalBalanceCheck(site,-1) === nothing),
    ]
    println()
    nfail = 0
    for (nm, ok) in checks; ok || (nfail += 1); @printf("  [%s] %s\n", ok ? "PASS" : "FAIL", nm); end
    println(nfail == 0 ? "\n  ★ FATES stable + demographically active over $ndays days" : "\n  ✗ $nfail check(s) failed")
    return nfail
end
# Run as a script; when `include`d (e.g. by the Metal FATES harness reusing build()),
# just define the functions.
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
