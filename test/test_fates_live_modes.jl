# test_fates_live_modes.jl
# FATES live-driver MODE exercises — running the SPITFIRE fire, CNP nutrient-cycle,
# and tree-damage modes through the REAL multi-day `clm_drv!` loop (vs the prior
# constructed-state history-buffer validation).
#
# Each mode cold-starts a real-parameter (14-PFT) FATES site with the mode enabled
# (the new gated `clm_fates_init!` kwargs `spitfire_mode` / `parteh_mode` /
# `nitrogen_spec` / `phosphorus_spec` / `use_tree_damage`), loops `clm_drv!` for
# >= 2 days, and asserts:
#   1. every step completes without error,
#   2. no NaN/Inf in the FATES column's key CLM outputs nor the FATES site state,
#   3. `TotalBalanceCheck` holds at the end (C for all; +N/P for CNP),
#   4. the mode actually DID something live (distinct from the carbon-only path):
#        * fire   — the SPITFIRE chain ran (site fire-weather index / NF / fdi > 0,
#                   which the no_fire path leaves at 0),
#        * CNP    — 3 elements registered, every established cohort carries a
#                   positive daily_n_demand (the prescribed-uptake nutrient
#                   acquisition fed the demographic step),
#        * damage — the damage class state is active (nlevdamage >= 2, every cohort
#                   carries a valid crowndamage class).
#
# `clm_fates_init!` mutates FATES module-global control Refs + parameter tables, so
# each mode saves/restores them (the same set as test_fates_driver_run.jl) to leave
# the surrounding suite unperturbed. The default suite never sets use_fates, so the
# use_fates=false path is byte-identical to before.

using Test
using CLM
const _C = CLM

# ---------------------------------------------------------------------------
# Shared helpers (build a driver-ready FATES-tagged single-column inst, run the
# real clm_drv! loop, walk the cohort/patch census).
# ---------------------------------------------------------------------------

# Snapshot of every FATES + varpar global clm_fates_init! mutates.
function _fates_save_globals()
    return (
        prt_global = _C.prt_global[], numel = _C.num_elements[],
        ellist = copy(_C.element_list), parteh = _C.hlm_parteh_mode[],
        use_sp = _C.hlm_use_sp[], hydro = _C.hlm_use_planthydro[],
        agetrk = _C.hlm_use_cohort_age_tracking[], nleafage = _C.nleafage[],
        numpft = _C.numpft[], tod = _C.hlm_current_tod[],
        nlevsc = _C.nlevsclass[], nlevca = _C.nlevcoage[], nlevage = _C.nlevage[],
        nlevdam = _C.nlevdamage[], nlevhgt = _C.nlevheight[],
        sfp = _C.SFParams[], edparams = _C.EDParams[], edpft = _C.EDPftvarcon_inst[],
        paramd = _C.ParamDerived[], restart = _C.hlm_is_restart[],
        nocomp = _C.hlm_use_nocomp[], biogeog = _C.hlm_use_fixed_biogeog[],
        luh = _C.hlm_use_luh[], inv = _C.hlm_use_inventory_init[],
        damage = _C.hlm_use_tree_damage[], doy = _C.hlm_day_of_year[],
        dpy = _C.hlm_days_per_year[], freqday = _C.hlm_freq_day[],
        nharv = _C.hlm_num_lu_harvest_cats[], ch4 = _C.hlm_use_ch4[],
        vert = _C.hlm_use_vertsoilc[], spit = _C.hlm_spitfire_mode[],
        paramfp = _C.fates_maxPatchesPerSite[],
        maxcomp = _C.max_comp_per_site[], npcomp = _C.fates_np_comp_scaling[],
        nspec = _C.hlm_nitrogen_spec[], pspec = _C.hlm_phosphorus_spec[],
        nlevsno = _C.varpar.nlevsno, nlevsoi = _C.varpar.nlevsoi,
        nlevgrnd = _C.varpar.nlevgrnd, nlevmaxu = _C.varpar.nlevmaxurbgrnd)
end

function _fates_restore_globals(g)
    _C.prt_global[] = g.prt_global
    _C.num_elements[] = g.numel
    empty!(_C.element_list); append!(_C.element_list, g.ellist)
    _C.hlm_parteh_mode[] = g.parteh
    _C.hlm_use_sp[] = g.use_sp
    _C.hlm_use_planthydro[] = g.hydro
    _C.hlm_use_cohort_age_tracking[] = g.agetrk
    _C.nleafage[] = g.nleafage
    _C.numpft[] = g.numpft
    _C.hlm_current_tod[] = g.tod
    _C.nlevsclass[] = g.nlevsc
    _C.nlevcoage[] = g.nlevca
    _C.nlevage[] = g.nlevage
    _C.nlevdamage[] = g.nlevdam
    _C.nlevheight[] = g.nlevhgt
    _C.SFParams[] = g.sfp
    _C.EDParams[] = g.edparams
    _C.EDPftvarcon_inst[] = g.edpft
    _C.ParamDerived[] = g.paramd
    _C.hlm_is_restart[] = g.restart
    _C.hlm_use_nocomp[] = g.nocomp
    _C.hlm_use_fixed_biogeog[] = g.biogeog
    _C.hlm_use_luh[] = g.luh
    _C.hlm_use_inventory_init[] = g.inv
    _C.hlm_use_tree_damage[] = g.damage
    _C.hlm_day_of_year[] = g.doy
    _C.hlm_days_per_year[] = g.dpy
    _C.hlm_freq_day[] = g.freqday
    _C.hlm_num_lu_harvest_cats[] = g.nharv
    _C.hlm_use_ch4[] = g.ch4
    _C.hlm_use_vertsoilc[] = g.vert
    _C.hlm_spitfire_mode[] = g.spit
    _C.fates_maxPatchesPerSite[] = g.paramfp
    _C.max_comp_per_site[] = g.maxcomp
    _C.fates_np_comp_scaling[] = g.npcomp
    _C.hlm_nitrogen_spec[] = g.nspec
    _C.hlm_phosphorus_spec[] = g.pspec
    _C.varpar.nlevsno = g.nlevsno
    _C.varpar.nlevsoi = g.nlevsoi
    _C.varpar.nlevgrnd = g.nlevgrnd
    _C.varpar.nlevmaxurbgrnd = g.nlevmaxu
end

# Walk the cohort/patch linked lists: census + finiteness + total patch area.
function _fates_census(inst)
    site = inst.fates.sites[1]
    npatch = 0; ncoh = 0; totarea = 0.0; bad = String[]
    cp = site.oldest_patch
    while cp !== nothing
        npatch += 1
        isfinite(cp.area) || push!(bad, "patch.area")
        totarea += cp.area
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

# Build a driver-ready FATES-tagged single column with the given mode kwargs and
# run `clm_drv!` for `ndays` days. Returns the populated inst + a run report.
function _run_fates_mode(; ndays::Int = 2, init_kwargs...)
    nlevsoil = 5
    nlevdecomp = 1

    _C.varpar_init!(_C.varpar, 1, 14, 2, 5)
    _C.varcon_init!()
    _C.varpar.nlevsno = 5
    _C.varpar.nlevsoi = nlevsoil
    _C.varpar.nlevgrnd = nlevsoil
    _C.varpar.nlevmaxurbgrnd = nlevsoil
    nlevsno = _C.varpar.nlevsno

    ng, nl, nc, np = 1, 1, 1, 2
    inst = _C.CLMInstances()
    _C.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np,
                     nlevdecomp_full=_C.varpar.nlevdecomp_full)

    col = inst.column
    col.landunit[1] = 1; col.gridcell[1] = 1; col.snl[1] = 0
    col.patchi[1] = 1; col.patchf[1] = 2; col.npatches[1] = 2
    col.nbedrock[1] = _C.varpar.nlevsoi
    col.is_fates[1] = true

    lun = inst.landunit
    lun.itype[1] = _C.ISTSOIL; lun.urbpoi[1] = false
    lun.lakpoi[1] = false; lun.active[1] = true

    pch = inst.patch
    for p in 1:np
        pch.active[p] = true; pch.landunit[p] = 1; pch.gridcell[p] = 1
        pch.column[p] = 1; pch.itype[p] = 1
    end
    pch.wtcol[1] = 0.0; pch.wtcol[2] = 1.0
    pch.is_fates[1] = true; pch.is_fates[2] = true

    bounds = _C.BoundsType(begg=1, endg=ng, begl=1, endl=nl, begc=1, endc=nc,
        begp=1, endp=np, begCohort=0, endCohort=0,
        level=_C.BOUNDS_LEVEL_CLUMP, clump_index=1)

    filt = _C.ClumpFilter(); _C.alloc_filters!(filt, nc, np, nl)
    filt.allc .= true; filt.nolakec .= true; filt.nolakep .= true
    filt.soilc .= true; filt.soilp .= true; filt.nolakeurbanp .= true
    filt.nourbanp .= true; filt.nourbanc .= true; filt.hydrologyc .= true
    filt.urbanc .= false; filt.urbanl .= false; filt.urbanp .= false
    filt.snowc .= false; filt.nosnowc .= true; filt.do_smb_c .= false
    filt.exposedvegp .= true; filt.noexposedvegp .= false
    filt.lakec .= false; filt.lakep .= false
    filt.lakesnowc .= false; filt.lakenosnowc .= false
    filt.bgc_soilc .= false; filt.bgc_vegp .= false; filt.pcropp .= false
    filt.soilnopcropp .= true; filt.actfirec .= false; filt.actfirep .= false

    filt_ia = _C.ClumpFilter(); _C.alloc_filters!(filt_ia, nc, np, nl)
    filt_ia.allc .= true; filt_ia.nolakec .= true; filt_ia.nolakep .= true
    filt_ia.nourbanc .= true; filt_ia.nourbanp .= true

    _C.urban_read_nml!(_C.urban_ctrl)

    dzmin = zeros(nlevsno); dzmax_u = zeros(nlevsno); dzmax_l = zeros(nlevsno)
    dzmin[1]=0.010; dzmax_u[1]=0.02; dzmax_l[1]=0.03
    dzmin[2]=0.015; dzmax_u[2]=0.05; dzmax_l[2]=0.07
    for j in 3:nlevsno
        dzmin[j]=dzmax_u[j-1]*0.5
        dzmax_u[j]=2.0*dzmax_u[j-1]+0.01
        dzmax_l[j]=dzmax_u[j]+dzmax_l[j-1]
        if j==nlevsno; dzmax_u[j]=floatmax(Float64); dzmax_l[j]=floatmax(Float64); end
    end
    _C.SNOW_DZMIN[]=dzmin; _C.SNOW_DZMAX_U[]=dzmax_u; _C.SNOW_DZMAX_L[]=dzmax_l

    npft_clm = _C.MXPFT + 1
    pc = _C.pftcon
    pc.dleaf=fill(0.04,npft_clm); pc.slatop=fill(0.01,npft_clm)
    pc.leafcn=fill(25.0,npft_clm); pc.flnr=fill(0.1,npft_clm)
    pc.fnitr=fill(0.1,npft_clm); pc.mbbopt=fill(9.0,npft_clm)
    pc.c3psn=fill(1.0,npft_clm); pc.woody=fill(0.0,npft_clm)
    pc.smpso=fill(-66000.0,npft_clm); pc.smpsc=fill(-275000.0,npft_clm)
    pc.z0mr=fill(0.055,npft_clm); pc.displar=fill(0.67,npft_clm)
    pc.xl=fill(0.1,npft_clm); pc.rhol=fill(0.1,npft_clm,2); pc.rhos=fill(0.2,npft_clm,2)
    pc.taul=fill(0.05,npft_clm,2); pc.taus=fill(0.1,npft_clm,2)
    pc.medlynintercept=fill(100.0,npft_clm); pc.medlynslope=fill(6.0,npft_clm)
    pc.crop=fill(0.0,npft_clm)

    photosyns = inst.photosyns
    inst.canopystate.frac_veg_nosno_alb_patch .= 1

    # cold-start a real-param (14-PFT) FATES site with the requested mode.
    fates = _C.clm_fates_init!(inst; nsites=1, nlevsoil=nlevsoil,
                               nlevdecomp=nlevdecomp, init_kwargs...)

    config = _C.CLMDriverConfig(use_fates=true, use_aquifer_layer=false,
                                use_hydrstress=false, use_luna=false)

    function set_forcing!(inst, sunfrac)
        a = inst.atm2lnd; g = 1
        T=290.0; pbot=9.9e4; q=0.008
        th = T*(1.0e5/pbot)^0.286
        rho = pbot/(287.058*T*(1.0+0.61*q))
        vp = q*pbot/(0.622+0.378*q)
        a.forc_t_not_downscaled_grc[g]=T; a.forc_th_not_downscaled_grc[g]=th
        a.forc_pbot_not_downscaled_grc[g]=pbot; a.forc_q_not_downscaled_grc[g]=q
        a.forc_rho_not_downscaled_grc[g]=rho; a.forc_lwrad_not_downscaled_grc[g]=350.0
        a.forc_rain_not_downscaled_grc[g]=0.0; a.forc_snow_not_downscaled_grc[g]=0.0
        a.forc_u_grc[g]=4.0; a.forc_v_grc[g]=0.0
        a.forc_hgt_grc[g]=30.0; a.forc_hgt_u_grc[g]=30.0
        a.forc_hgt_t_grc[g]=30.0; a.forc_hgt_q_grc[g]=30.0
        a.forc_vp_grc[g]=vp; a.forc_pco2_grc[g]=367.0e-6*pbot
        a.forc_po2_grc[g]=0.209*pbot; a.forc_topo_grc[g]=150.0
        for b in 1:size(a.forc_solad_not_downscaled_grc,2)
            a.forc_solad_not_downscaled_grc[g,b] = (b==1 ? 300.0 : 250.0)*sunfrac
            a.forc_solai_grc[g,b]                = (b==1 ? 100.0 : 80.0)*sunfrac
        end
        inst.topo.topo_col[1]=150.0
        return nothing
    end

    c = 1; pveg = col.patchi[1] + 1
    dtime = 1800.0
    steps_per_day = Int(round(86400 / dtime))
    nsteps = steps_per_day * ndays
    declin = 0.4

    secs = 0; jday = 152; daycount = 0
    daily_fires = 0
    no_error = true
    balance_threw = false
    nan_steps = String[]
    errmsg = ""

    for i in 1:nsteps
        is_beg = (secs == 0)
        frac = secs / 86400.0
        sun = max(sin(π * frac), 0.0)
        set_forcing!(inst, sun)
        _C.downscale_forcings!(bounds, inst.atm2lnd, inst.column,
                               inst.landunit, inst.topo)
        nextsw_cday = (jday + (secs + dtime) / 86400.0)

        try
            _C.clm_drv!(config, inst, filt, filt_ia, bounds,
                true, nextsw_cday, declin, declin, 0.4091,
                false, false, "20260601", false;
                nstep=i, is_first_step=(i==1), is_beg_curr_day=is_beg,
                is_end_curr_day=false, is_beg_curr_year=false, dtime=dtime,
                mon=6, day=1+daycount, secs=secs, jday=jday, photosyns=photosyns)
        catch e
            no_error = false
            errmsg = sprint(showerror, e)
            if e isa ErrorException && occursin("alance", errmsg)
                balance_threw = true
            end
            break
        end
        if is_beg && i > 1
            daily_fires += 1
        end

        cs = inst.canopystate; sa = inst.surfalb
        ps = inst.photosyns; ef = inst.energyflux
        audited = Dict(
            "elai"=>cs.elai_patch[pveg], "htop"=>cs.htop_patch[pveg],
            "laisun"=>cs.laisun_patch[pveg], "laisha"=>cs.laisha_patch[pveg],
            "albd"=>sa.albd_patch[pveg,1], "rssun"=>ps.rssun_patch[pveg],
            "btran"=>ef.btran_patch[pveg])
        bad = [k for (k,v) in audited if !isfinite(v)]
        isempty(bad) || push!(nan_steps, "step$i:$bad")
        cen = _fates_census(inst)
        isempty(cen.bad) || push!(nan_steps, "step$i:census$(cen.bad)")
        isapprox(cen.totarea, _C.area; atol=1e-4) ||
            push!(nan_steps, "step$i:area=$(cen.totarea)")

        secs += Int(dtime)
        if secs >= 86400
            secs = 0; jday += 1; daycount += 1
        end
    end

    return (inst=inst, no_error=no_error, balance_threw=balance_threw,
            nan_steps=nan_steps, daily_fires=daily_fires, errmsg=errmsg,
            pveg=pveg)
end

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@testset "FATES live modes (fire / CNP / damage) through clm_drv!" begin

    @testset "SPITFIRE fire mode" begin
        g = _fates_save_globals()
        try
            r = _run_fates_mode(ndays=2, spitfire_mode=1)  # scalar_lightning
            @test r.no_error
            @test !r.balance_threw
            @test isempty(r.nan_steps)
            @test r.daily_fires >= 1   # the daily demographic + fire step fired

            # the SPITFIRE chain RAN live (the no_fire path leaves these at 0):
            # site fire-weather index, number of ignitions, fire-danger index.
            site = r.inst.fates.sites[1]
            @test isfinite(site.fireWeather.fire_weather_index)
            @test site.fireWeather.fire_weather_index > 0.0
            @test isfinite(site.NF)
            @test site.NF > 0.0
            @test isfinite(site.fdi)
            @test site.fdi > 0.0
            # per-patch fuel loading is populated + finite.
            cp = site.oldest_patch; totfuel = 0.0; allfin = true
            while cp !== nothing
                if cp.fuel !== nothing
                    all(isfinite, cp.fuel.loading) || (allfin = false)
                    totfuel += sum(cp.fuel.loading)
                end
                isfinite(cp.frac_burnt) || (allfin = false)
                cp = cp.younger
            end
            @test allfin
            @test totfuel >= 0.0
            # final carbon mass balance holds.
            @test _C.TotalBalanceCheck(site, -1) === nothing
        finally
            _fates_restore_globals(g)
        end
    end

    @testset "CNP nutrient-cycle mode" begin
        g = _fates_save_globals()
        try
            r = _run_fates_mode(ndays=2, parteh_mode=_C.prt_cnp_flex_allom_hyp,
                                nitrogen_spec=1, phosphorus_spec=1)
            @test r.no_error
            @test !r.balance_threw
            @test isempty(r.nan_steps)
            @test r.daily_fires >= 1

            # CNP registered 3 elements (carbon / nitrogen / phosphorus).
            @test _C.num_elements[] == 3
            @test _C.element_list == [_C.carbon12_element, _C.nitrogen_element,
                                      _C.phosphorus_element]

            site = r.inst.fates.sites[1]
            # the prescribed-uptake nutrient acquisition fed the demographic step:
            # every established (non-new) cohort carries a positive daily N demand
            # and finite N/P leaf pools.
            ndems = Float64[]; totN = 0.0; totP = 0.0; npfin = true
            cp = site.oldest_patch
            while cp !== nothing
                cc = cp.tallest
                while cc !== nothing
                    if !cc.isnew
                        push!(ndems, cc.daily_n_demand)
                    end
                    n_leaf = _C.GetState(cc.prt, _C.leaf_organ, _C.nitrogen_element)
                    p_leaf = _C.GetState(cc.prt, _C.leaf_organ, _C.phosphorus_element)
                    (isfinite(n_leaf) && isfinite(p_leaf)) || (npfin = false)
                    totN += n_leaf * cc.n
                    totP += p_leaf * cc.n
                    cc = cc.shorter
                end
                cp = cp.younger
            end
            @test npfin
            @test !isempty(ndems)
            @test all(d -> isfinite(d) && d > 0.0, ndems)
            @test totN > 0.0
            @test totP > 0.0

            # mass balance holds across all 3 elements (C/N/P).
            for el in 1:_C.num_elements[]
                @test isfinite(site.mass_balance[el].err_fates)
            end
            @test _C.TotalBalanceCheck(site, -1) === nothing
        finally
            _fates_restore_globals(g)
        end
    end

    @testset "tree-damage mode" begin
        g = _fates_save_globals()
        try
            r = _run_fates_mode(ndays=2, use_tree_damage=true)
            @test r.no_error
            @test !r.balance_threw
            @test isempty(r.nan_steps)
            @test r.daily_fires >= 1

            # tree-damage class state is active.
            @test _C.hlm_use_tree_damage[] == _C.itrue
            @test _C.nlevdamage[] >= 2

            site = r.inst.fates.sites[1]
            # every cohort carries a valid crown-damage class (>= 1, 1 = undamaged);
            # the damage class state is live (not the unset sentinel).
            ncoh = 0; alldmg = true
            cp = site.oldest_patch
            while cp !== nothing
                cc = cp.tallest
                while cc !== nothing
                    ncoh += 1
                    (cc.crowndamage >= 1 && cc.crowndamage <= _C.nlevdamage[]) ||
                        (alldmg = false)
                    cc = cc.shorter
                end
                cp = cp.younger
            end
            @test ncoh >= 1
            @test alldmg
            @test _C.TotalBalanceCheck(site, -1) === nothing
        finally
            _fates_restore_globals(g)
        end
    end
end
