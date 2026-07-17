#!/usr/bin/env julia
# ==========================================================================
# fates_multisite_validation.jl — multi-site FATES equilibrium/stability
# validation with climate-appropriate PFT screening (PR #197).
#
# fates_longhorizon.jl proved ONE tropical site (Aripuana) reaches a bounded
# carbon band under the :drop_cold_deciduous screen. This harness generalizes
# that to a SITE REGISTRY spanning contrasting biomes, each seeded with the
# climatically-correct FATES PFT set (either a named screen or an explicit
# per-FATES-PFT area vector), and scores each stand for a *sensible multi-year
# equilibrium*:
#   - carbon reaches a bounded band (not boom-bust, not monotone collapse to 0),
#   - LAI / maxdbh / cohort count stay physically reasonable for that biome,
#   - carbon mass is conserved EVERY day (daily TotalBalanceCheck holds),
#   - the census stays physical every day (n>=0, dbh>0, area sums to AREA).
#
# FATES carries MODULE-GLOBAL state (param Refs, numpft, the identity HLM<->FATES
# crosswalk the fixed-biogeog path installs), so — like fates_longhorizon.jl —
# each site MUST run in its OWN process. This script runs ONE site per invocation
# (chosen by SITE=... env or ARGS[1]); a companion driver launches the sites as
# parallel background processes and aggregates the per-site summary lines.
#
#   SITE=aripuana_screened FATES_NDAYS=1460 julia +1.12 --project=. \
#       scripts/fates_multisite_validation.jl
#
# Each run appends one summary line to $FATES_MS_RESULTS (default scratchpad).
# ==========================================================================
using CLM, Printf, Dates
const _C = CLM
const DATA = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data"

# ---- FATES default PFT phenology (data/fates/fates_params_default.cdl) --------
# idx name                                 evergrn season_decid stress_decid
#  1  broadleaf_evergreen_tropical_tree       1        0            0
#  2  needleleaf_evergreen_extratrop_tree     1        0            0
#  3  needleleaf_colddecid_extratrop_tree     0        1            0
#  4  broadleaf_evergreen_extratrop_tree      1        0            0
#  5  broadleaf_hydrodecid_tropical_tree      0        0            1
#  6  broadleaf_colddecid_extratrop_tree      0        1            0
#  7  broadleaf_evergreen_extratrop_shrub     1        0            0
#  8  broadleaf_hydrodecid_extratrop_shrub    0        0            1
#  9  broadleaf_colddecid_extratrop_shrub     0        1            0
# 10  broadleaf_evergreen_arctic_shrub        1        0            0
# 11  broadleaf_colddecid_arctic_shrub        0        1            0
# 12  arctic_c3_grass                         0        1            0
# 13  cool_c3_grass                           0        0            1
# 14  c4_grass                                0        0            1
#
# Explicit climate-appropriate area vectors (equal presence among the chosen set;
# 0 excludes the PFT from seeding + recruitment). Boreal/temperate sites CANNOT
# use the :drop_cold_deciduous screen (that keeps warm/tropical PFTs and drops the
# cold-deciduous trees that ARE appropriate there), so they pass an explicit set.
const NPFT = 14
_areafrac(idxs) = [i in idxs ? 1.0 : 0.0 for i in 1:NPFT]
# Boreal tree set: needleleaf evergreen (2) + needleleaf colddecid (3) +
#                  broadleaf colddecid tree (6).
const AREAFRAC_BOREAL   = _areafrac((2, 3, 6))
# Temperate set: broadleaf colddecid tree (6, dominant) + needleleaf evergreen
#                (2) + broadleaf evergreen extratrop (4).
const AREAFRAC_TEMPERATE = _areafrac((2, 4, 6))

# ---- Site registry -----------------------------------------------------------
# soil_tk: cold-start applies a generic ~272 K partly-frozen deep profile; we warm
# the FATES column's soil+ground to a climate-appropriate GROWING-SEASON value so
# frozen-root exclusion doesn't spuriously zero water uptake (real leaf temp still
# comes from the energy balance on the real forcing). Tropical ~299, boreal ~281,
# temperate ~286 K.
# soil_h2o: OPTIONAL cold-start root-zone LIQUID moisture prime (fraction of watsat),
# mirroring soil_tk. The generic surfdata cold start loads the FATES column with
# 0.75*watsat of water, but at the ~272 K cold-start temperature that water is stored
# entirely as ICE; the harness then warms the soil (soil_tk) but never thaws the ice.
# Over a wrapped single-year of forcing the ice melts and the meltwater DRAINS away, so
# by growing season the boreal root zone is bone-dry (smp << wilting) -> btran=0 -> GPP=0
# -> the stand starves on maintenance respiration. A real spun-up boreal soil holds
# field-capacity liquid water through the growing season, so this is a cold-start /
# short-spin artifact, not real physics. `soil_h2o` re-asserts that growing-season
# liquid state at cold start for cold/under-moist sites (see build_site).
# `nothing` = no prime (tropical/temperate already reach healthy equilibria — leave them
# byte-identical). Bounded by watsat, so it can never exceed pore space.
struct SiteCfg
    key::String; label::String; domain::String; fyr::Int
    screen::Symbol; areafrac::Union{Nothing,Vector{Float64}}
    soil_tk::Float64; band_lo::Float64; band_hi::Float64
    soil_h2o::Union{Nothing,Float64}
end
# convenience ctor keeping the pre-existing positional call sites working (no prime)
SiteCfg(key,label,domain,fyr,screen,areafrac,soil_tk,band_lo,band_hi) =
    SiteCfg(key,label,domain,fyr,screen,areafrac,soil_tk,band_lo,band_hi,nothing)
const SITES = Dict(
    "aripuana_screened" => SiteCfg("aripuana_screened",
        "Aripuana Amazon (tropical) — SCREENED :drop_cold_deciduous",
        "domain_Aripuana_Amazon", 2004, :drop_cold_deciduous, nothing, 299.0, 800.0, 4000.0),
    "aripuana_baseline" => SiteCfg("aripuana_baseline",
        "Aripuana Amazon (tropical) — BASELINE all-14-PFT (boom-bust)",
        "domain_Aripuana_Amazon", 2004, :none, nothing, 299.0, 0.0, Inf),
    "krycklan_screened" => SiteCfg("krycklan_screened",
        "Krycklan (boreal Sweden) — SCREENED boreal tree set {2,3,6}",
        "domain_Boreal_Krycklan_Sweden", 2010, :none, AREAFRAC_BOREAL, 281.0, 200.0, 6000.0, 0.75),
    "krycklan_baseline" => SiteCfg("krycklan_baseline",
        "Krycklan (boreal Sweden) — BASELINE all-14-PFT",
        "domain_Boreal_Krycklan_Sweden", 2010, :none, nothing, 281.0, 0.0, Inf),
    "hubbardbrook_screened" => SiteCfg("hubbardbrook_screened",
        "Hubbard Brook (temperate USA) — SCREENED temperate set {2,4,6}",
        "domain_Temperate_HubbardBrook_USA", 2015, :none, AREAFRAC_TEMPERATE, 286.0, 200.0, 6000.0),
    "hubbardbrook_baseline" => SiteCfg("hubbardbrook_baseline",
        "Hubbard Brook (temperate USA) — BASELINE all-14-PFT",
        "domain_Temperate_HubbardBrook_USA", 2015, :none, nothing, 286.0, 0.0, Inf),
)

# ---- census / carbon helpers (mirror fates_longhorizon.jl) -------------------
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

# ---- growing-season root-zone moisture prescription --------------------------
# Prescribe a spun-up boreal soil-moisture boundary condition on the FATES column's
# root zone: set each (thawed) soil layer to `frac`*watsat of LIQUID water, zeroing its
# ice. This is NOT a water-conservation claim — it is a prescribed-soil-moisture forcing
# (a standard land-model technique for isolating vegetation dynamics from water-balance
# spin-up), applied only to sites that carry a soil_h2o. The short single-year wrapped
# forcing here cannot reproduce the multi-decade boreal water balance: the cold-start
# ice load drains below field capacity before leaf-on, so without this the boreal root
# zone is desert-dry (smp << wilting) all growing season -> btran=0 -> GPP=0 -> collapse.
# A real spun-up boreal forest's thawed active layer holds ~field-capacity water through
# the growing season. `thawed_only=true` (the daily re-assertion) touches only layers
# above freezing, so it never fights the freeze/thaw physics; frozen layers are correctly
# left as non-uptake layers. Bounded by watsat -> never exceeds pore space. Carbon
# conservation (the harness's actual invariant) is unaffected.
# Effective moisture-prime fraction: the per-site cfg.soil_h2o, overridable for tuning
# via FATES_SOIL_H2O (returns `nothing` unchanged so unprimed sites stay unprimed).
_soil_h2o(cfg::SiteCfg) = cfg.soil_h2o === nothing ? nothing :
    something(tryparse(Float64, get(ENV, "FATES_SOIL_H2O", "")), cfg.soil_h2o)

function prime_fates_soil_moisture!(inst, bounds, frac::Float64; thawed_only::Bool)
    col=inst.column; ss=inst.soilstate; temp=inst.temperature
    ws=inst.water.waterstatebulk_inst.ws
    joff=_C.varpar.nlevsno; nsoi=_C.varpar.nlevsoi; np=0
    for c in 1:bounds.endc
        col.is_fates[c] || continue
        np += 1
        for j in 1:nsoi
            ws_j = ss.watsat_col[c, j]
            (isfinite(ws_j) && ws_j > 0.0) || continue
            thawed_only && temp.t_soisno_col[c, joff+j] <= _C.TFRZ && continue
            vol = frac * ws_j
            ws.h2osoi_vol_col[c, j]      = vol
            ws.h2osoi_liq_col[c, joff+j] = vol * col.dz[c, joff+j] * _C.DENH2O
            ws.h2osoi_ice_col[c, joff+j] = 0.0
        end
    end
    return np
end

# ---- per-site build (parameterized clone of fates_longhorizon.build) ---------
function build_site(cfg::SiteCfg)
    dom = "$DATA/$(cfg.domain)"
    fsurdat  = "$dom/settings/CLM/parameters/surfdata_clm.nc"
    paramfile = "$dom/settings/CLM/parameters/clm5_params.nc"
    inst, bounds, filt, _tm = _C.clm_initialize!(; fsurdat=fsurdat, paramfile=paramfile,
        use_fates=true, start_date=DateTime(cfg.fyr,1,1), dtime=1800,
        fates_biogeog_screen=cfg.screen, fates_pft_areafrac=cfg.areafrac)
    if cfg.areafrac !== nothing
        @printf("  fixed-biogeog: explicit PFT areafrac, PFTs=%s\n",
                string(findall(>(0.0), cfg.areafrac)))
    elseif cfg.screen != :none
        @printf("  fixed-biogeog screen: %s\n", cfg.screen)
    else
        @printf("  baseline: all-%d-PFT NBG cold start\n", NPFT)
    end
    for g in 1:bounds.endg
        isfinite(inst.gridcell.latdeg[g]) || (inst.gridcell.latdeg[g] = inst.gridcell.lat[g]*180/π)
        isfinite(inst.gridcell.londeg[g]) || (inst.gridcell.londeg[g] = inst.gridcell.lon[g]*180/π)
    end
    let col=inst.column, lun=inst.landunit, nflag=0
        for c in 1:bounds.endc
            if lun.itype[col.landunit[c]] == _C.ISTSOIL
                col.is_fates[c] = true; nflag += 1
            end
        end
        @printf("  flagged %d FATES column(s); fates.nsites=%d; lat=%.3f\n",
                nflag, inst.fates.nsites, inst.gridcell.latdeg[1])
    end
    # Climate-appropriate growing-season soil warm-up (see SiteCfg.soil_tk note).
    let col=inst.column, temp=inst.temperature, joff=_C.varpar.nlevsno, ngr=_C.varpar.nlevgrnd
        for c in 1:bounds.endc
            col.is_fates[c] || continue
            temp.t_grnd_col[c]=cfg.soil_tk
            for j in 1:ngr; temp.t_soisno_col[c, joff+j]=cfg.soil_tk; end
        end
    end
    # Climate-appropriate growing-season root-zone LIQUID moisture prescription (see the
    # SiteCfg.soil_h2o note + prime_fates_soil_moisture!). Applied once here at cold start;
    # run_site re-asserts it each growing-season day for THAWED layers, because the short
    # wrapped forcing drains the root zone below field capacity long before the growing
    # season otherwise. Gated on cfg.soil_h2o !== nothing so tropical/temperate stay
    # byte-identical.
    if _soil_h2o(cfg) !== nothing
        np = prime_fates_soil_moisture!(inst, bounds, _soil_h2o(cfg); thawed_only=false)
        @printf("  soil moisture prime: %d FATES col(s) root zone set to %.2f*watsat LIQUID\n", np, _soil_h2o(cfg))
    end
    config = _C.CLMDriverConfig(use_fates=true)
    filt_ia = _C.clump_filter_inactive_and_active
    return inst, inst.fates, config, bounds, filt, filt_ia
end

function run_site(cfg::SiteCfg, ndays::Int)
    println("="^74); println("  FATES multi-site validation — $(cfg.label)"); println("  horizon: $ndays days (~$(round(ndays/365,digits=2)) yr)"); println("="^74)
    inst, fates, config, bounds, filt, filt_ia = build_site(cfg)
    site = fates.sites[1]; photosyns = inst.photosyns
    dtime = 1800.0; steps_per_day = Int(round(86400/dtime)); nsteps = steps_per_day*ndays
    forcing = "$DATA/$(cfg.domain)/data/forcing/CLM_input/clmforc.$(cfg.fyr).nc"
    fr = _C.ForcingReader(); _C.forcing_reader_init!(fr, forcing); fr.interp_time = true
    start_date = DateTime(cfg.fyr, 1, 1); daycount = 0
    println("  forcing: clmforc.$(cfg.fyr).nc (single-year, wrapped)\n")

    c0 = total_site_carbon(site); cen0 = census(inst)
    @printf("  cold start: ncoh=%d npatch=%d carbon=%.4g dbh=%.4f\n\n", cen0.ncoh, cen0.npatch, c0, max_dbh(site))
    @printf("  %6s %6s %7s %12s %9s %9s %8s %8s %8s\n", "day", "ncoh", "npatch", "carbon", "maxdbh", "elai", "GPP/d", "NPP/d", "bal")
    println("  " * "-"^66)

    no_error=true; errmsg=""; nan_days=String[]; days_advanced=0; day_bal_ok=0
    min_ncoh=cen0.ncoh; max_ncoh=cen0.ncoh; max_npatch=cen0.npatch
    min_c=c0; max_c=c0; births=0; deaths=0; prev_ncoh=cen0.ncoh
    acc_gpp=0.0; acc_npp=0.0; day_gpp=0.0; day_npp=0.0
    elai_min=Inf; elai_max=-Inf; maxdbh_final=0.0
    carbon_series = Float64[]   # per-day carbon, for equilibrium analysis
    t0 = time()
    for i in 1:nsteps
        step_start = start_date + Second((i-1)*Int(dtime))
        is_beg = (Dates.hour(step_start)==0 && Dates.minute(step_start)==0 && Dates.second(step_start)==0)
        yr_off = Dates.year(step_start) - cfg.fyr
        forcing_time = yr_off == 0 ? step_start : step_start - Dates.Year(yr_off)
        _C.read_forcing_step!(fr, inst.atm2lnd, forcing_time, 1, 1;
            gridcell_latdeg=inst.gridcell.latdeg, gridcell_londeg=inst.gridcell.londeg, dtime=Int(dtime))
        inst.atm2lnd.forc_topo_grc[1]=200.0; inst.topo.topo_col[1]=200.0
        _C.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
        sod = Dates.hour(step_start)*3600 + Dates.minute(step_start)*60 + Dates.second(step_start)
        calday = Dates.dayofyear(step_start) + sod/86400.0
        (declin, _e) = _C.compute_orbital(calday); nextsw_cday = calday + Int(dtime)/86400.0
        # Re-assert the prescribed growing-season root-zone moisture for THAWED layers
        # once per day (see prime_fates_soil_moisture!). Only for sites carrying soil_h2o.
        if is_beg && _soil_h2o(cfg) !== nothing
            prime_fates_soil_moisture!(inst, bounds, _soil_h2o(cfg); thawed_only=true)
        end
        try
            _C.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw_cday, declin, declin, 0.4091,
                false, false, "20260101", false; nstep=i, is_first_step=(i==1), is_beg_curr_day=is_beg,
                is_end_curr_day=false, is_beg_curr_year=false, dtime=dtime, mon=Dates.month(step_start),
                day=Dates.day(step_start), secs=sod, jday=Dates.dayofyear(step_start), photosyns=photosyns)
        catch e
            no_error=false; errmsg=sprint(showerror, e)
            @printf("  x ERROR at step %d (day %d): %s\n", i, daycount, first(split(errmsg,"\n"))); break
        end
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
        # Optional soil-bc / btran diagnostic (FATES_PROBE=1). Mirrors the probe in
        # fates_longhorizon.jl: at midday on selected days it prints btran, the packed
        # FATES soil bc (smp_sl / h2o_liqvol_sl / eff_porosity_sl / tempk_sl) plus the
        # raw column liquid/ice split, to distinguish a genuinely DRY root zone from
        # water that is present but locked as ICE by the cold-start freeze.
        if get(ENV,"FATES_PROBE","")=="1" && Dates.hour(step_start)==12 && Dates.minute(step_start)==0 &&
                (daycount <= 8 || daycount % 15 == 0)
            btmax=0.0; gmax=0.0; vcmax=0.0; cp=site.oldest_patch
            while cp!==nothing
                try; btmax = max(btmax, maximum(x->isfinite(x) ? x : 0.0, cp.btran_ft)); catch; end
                cc=cp.tallest
                while cc!==nothing
                    isfinite(cc.gpp_tstep)&&(gmax=max(gmax,cc.gpp_tstep))
                    isfinite(cc.vcmax25top)&&(vcmax=max(vcmax,cc.vcmax25top))
                    cc=cc.shorter
                end
                cp=cp.younger
            end
            bc=inst.fates.bc_in[1]
            wsb=inst.water.waterstatebulk_inst; joff=_C.varpar.nlevsno
            fc = findfirst(inst.column.is_fates[1:bounds.endc]); fc===nothing && (fc=1)
            liq = [wsb.ws.h2osoi_liq_col[fc, joff+j] for j in 1:5]
            ice = [wsb.ws.h2osoi_ice_col[fc, joff+j] for j in 1:5]
            @printf("  [PROBE day%d %02d:00] btran=%.4f vc25=%.1f gpp=%.5g elai=%.3f tveg=%.2f\n",
                    daycount, Dates.hour(step_start), btmax, vcmax, gmax,
                    inst.canopystate.elai_patch[2], inst.temperature.t_veg_patch[2])
            @printf("        smp_sl(1:5)=%s  liqvol(1:5)=%s  eff_por(1:5)=%s  tempk(1:5)=%s\n",
                    string(round.(bc.smp_sl[1:5],digits=0)), string(round.(bc.h2o_liqvol_sl[1:5],digits=3)),
                    string(round.(bc.eff_porosity_sl[1:5],digits=3)), string(round.(bc.tempk_sl[1:5],digits=1)))
            @printf("        h2osoi_liq(1:5)=%s  h2osoi_ice(1:5)=%s [kg/m2]\n",
                    string(round.(liq,digits=2)), string(round.(ice,digits=2)))
            ep=_C.EDPftvarcon_inst[]
            @printf("        PFT2 smpsc=%.0f smpso=%.0f | PFT3 smpsc=%.0f smpso=%.0f\n",
                    ep.smpsc[2], ep.smpso[2], ep.smpsc[3], ep.smpso[3])
        end
        if is_beg && i>1
            daycount = (Dates.value(step_start - start_date) ÷ 86400_000)  # cumulative day index
            day_gpp=acc_gpp; day_npp=acc_npp; acc_gpp=0.0; acc_npp=0.0
            days_advanced += 1
            bal_ok = _C.TotalBalanceCheck(site, -1) === nothing
            bal_ok && (day_bal_ok += 1)
            cen = census(inst); c = total_site_carbon(site); md = max_dbh(site)
            el = inst.canopystate.elai_patch[2]
            isempty(cen.bad) || push!(nan_days, "day$daycount:$(cen.bad)")
            isfinite(c) || push!(nan_days, "day$daycount:carbon")
            dn = cen.ncoh - prev_ncoh; dn>0 && (births += dn); dn<0 && (deaths += -dn); prev_ncoh = cen.ncoh
            min_ncoh=min(min_ncoh,cen.ncoh); max_ncoh=max(max_ncoh,cen.ncoh); max_npatch=max(max_npatch,cen.npatch)
            min_c=min(min_c,c); max_c=max(max_c,c)
            isfinite(el) && (elai_min=min(elai_min,el); elai_max=max(elai_max,el))
            maxdbh_final=md
            push!(carbon_series, c)
            if daycount % max(1, ndays ÷ 30) == 0 || daycount == ndays
                @printf("  %6d %6d %7d %12.5g %9.4f %9.4f %8.3f %8.3f %8s\n", daycount, cen.ncoh, cen.npatch, c, md,
                        el, day_gpp, day_npp, bal_ok ? "ok" : "FAIL")
            end
        end
    end
    _C.forcing_reader_close!(fr)
    dt = time()-t0
    cenF = census(inst); cF = total_site_carbon(site)

    # ---- equilibrium analysis (second-half band + trend) --------------------
    n = length(carbon_series)
    half = n >= 4 ? carbon_series[(n÷2+1):end] : carbon_series
    band_lo = isempty(half) ? NaN : minimum(half)
    band_hi = isempty(half) ? NaN : maximum(half)
    band_ratio = (isfinite(band_lo) && band_lo > 0) ? band_hi/band_lo : Inf
    peak = isempty(carbon_series) ? NaN : maximum(carbon_series)
    early = n >= 1 ? carbon_series[max(1, n÷10)] : c0        # ~10% in
    # collapse: final <20% of peak AND declining.  boom-bust: rose >2.5x then fell >2.5x.
    collapsed = isfinite(cF) && isfinite(peak) && peak>0 && cF < 0.2*peak
    boom_bust = isfinite(peak) && isfinite(early) && early>0 && peak > 2.5*early && cF < peak/2.5
    in_band = isfinite(cF) && cF >= cfg.band_lo && cF <= cfg.band_hi

    println("\n  " * "-"^66)
    @printf("  ran %d/%d days in %.1f s (%.2f s/yr-sim)\n", daycount, ndays, dt, dt*365/max(daycount,1))
    @printf("  carbon:  cold=%.4g  final=%.4g  min=%.4g  max=%.4g  peak=%.4g\n", c0, cF, min_c, max_c, peak)
    @printf("  2nd-half band: [%.4g, %.4g]  (hi/lo=%.2f)\n", band_lo, band_hi, band_ratio)
    @printf("  ncoh:    range [%d, %d]   npatch max=%d   births=%d deaths=%d\n", min_ncoh, max_ncoh, max_npatch, births, deaths)
    @printf("  elai:    range [%.3f, %.3f]   maxdbh final=%.3f cm\n", elai_min, elai_max, maxdbh_final)
    @printf("  daily balance held: %d/%d days\n", day_bal_ok, days_advanced)

    checks = [
        ("no error over horizon",                 no_error),
        ("no non-finite day",                     isempty(nan_days)),
        ("carbon conserved EVERY day",            day_bal_ok == days_advanced && days_advanced >= ndays-1),
        ("cohorts persisted (ncoh>=1)",           min_ncoh >= 1),
        ("demography advanced (births+deaths)",   births + deaths > 0),
        ("carbon finite & positive",              isfinite(cF) && cF > 0),
        ("no collapse to ~0 (final>=20% peak)",   !collapsed),
        ("no boom-bust (rise&fall >2.5x)",        !boom_bust),
        ("final in expected biome band",          in_band),
        ("final balance holds",                   _C.TotalBalanceCheck(site,-1) === nothing),
    ]
    println()
    nfail = 0
    for (nm, ok) in checks; ok || (nfail += 1); @printf("  [%s] %s\n", ok ? "PASS" : "FAIL", nm); end
    println(nfail == 0 ? "\n  * $(cfg.key): healthy multi-year equilibrium ($ndays days)" :
                         "\n  x $(cfg.key): $nfail check(s) failed")

    # ---- append machine-readable summary line -------------------------------
    resfile = get(ENV, "FATES_MS_RESULTS",
        "$(get(ENV,"TMPDIR","/tmp"))fates_multisite_results.csv")
    open(resfile, "a") do io
        @printf(io, "%s\t%s\t%d\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.2f\t%d\t%d\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\n",
            cfg.key, cfg.label, daycount, c0, cF, min_c, max_c, peak, band_ratio,
            min_ncoh, max_ncoh, isfinite(elai_min) ? round(Int,elai_min*1000)/1000 : -1,
            isfinite(elai_max) ? round(Int,elai_max*1000)/1000 : -1, maxdbh_final,
            day_bal_ok, days_advanced, collapsed ? 1 : 0, boom_bust ? 1 : 0,
            in_band ? 1 : 0, nfail)
    end
    println("  summary -> $resfile")
    return nfail
end

function main()
    site_key = length(ARGS) >= 1 ? ARGS[1] : get(ENV, "SITE", "aripuana_screened")
    haskey(SITES, site_key) || (println("unknown SITE=$site_key; known: $(sort(collect(keys(SITES))))"); return 2)
    ndays = parse(Int, get(ENV, "FATES_NDAYS", "1460"))  # default 4 yr
    return run_site(SITES[site_key], ndays)
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
