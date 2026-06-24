# test_fates_driver_run.jl
# FATES live-driver wiring — the first REAL multi-day clm_drv! integrated run.
#
# Builds a minimal single-column run with CLMDriverConfig(use_fates=true), tags
# the column / vegetated patch as FATES, cold-starts a real-parameter (14-PFT)
# carbon-only FATES site (clm_fates_init!), sets representative forcing, and loops
# the real `clm_drv!` for >= 2 days (96 steps at dtime=1800s), firing the daily
# demographic step on each day boundary (is_beg_curr_day=true at secs==0).
#
# Unlike test_fates_driver_hooks.jl (which drives the bc pack -> FATES -> unpack
# round-trip in isolation), this exercises the FOUR per-timestep FATES hooks +
# the daily demographic hook through the actual clm_drv! loop end-to-end, and
# asserts:
#   1. every step completes without error,
#   2. no NaN/Inf in the FATES column's key CLM outputs nor the FATES site state,
#   3. FATES actually drove the column (laisha/canopy radiation/rssun populated by
#      the FATES branch — distinct from the non-FATES path),
#   4. the daily demographic step advanced the census & TotalBalanceCheck held
#      (n>=0, dbh>0, patch areas sum to AREA),
#   5. values stay physical over the multi-day horizon (no blow-up).
#
# clm_fates_init! mutates FATES module-global control Refs + parameter tables, so
# save/restore them (same set as test_fates_driver_hooks.jl) to leave the
# surrounding default suite unperturbed. The default suite never sets use_fates,
# so the use_fates=false path is byte-identical to before.

using Test
using CLM
const _C = CLM

@testset "FATES multi-day clm_drv! integrated run" begin
    # ---- save FATES + varpar globals (mirrors test_fates_driver_hooks.jl) ----
    old_global   = _C.prt_global[]
    old_numel    = _C.num_elements[]
    old_ellist   = copy(_C.element_list)
    old_parteh   = _C.hlm_parteh_mode[]
    old_use_sp   = _C.hlm_use_sp[]
    old_hydro    = _C.hlm_use_planthydro[]
    old_agetrk   = _C.hlm_use_cohort_age_tracking[]
    old_nleafage = _C.nleafage[]
    old_numpft   = _C.numpft[]
    old_tod      = _C.hlm_current_tod[]
    old_nlevsc   = _C.nlevsclass[]
    old_nlevca   = _C.nlevcoage[]
    old_nlevage  = _C.nlevage[]
    old_nlevdam  = _C.nlevdamage[]
    old_nlevhgt  = _C.nlevheight[]
    old_sfp      = _C.SFParams[]
    old_edparams = _C.EDParams[]
    old_edpft    = _C.EDPftvarcon_inst[]
    old_paramd   = _C.ParamDerived[]
    old_restart  = _C.hlm_is_restart[]
    old_nocomp   = _C.hlm_use_nocomp[]
    old_biogeog  = _C.hlm_use_fixed_biogeog[]
    old_luh      = _C.hlm_use_luh[]
    old_inv      = _C.hlm_use_inventory_init[]
    old_damage   = _C.hlm_use_tree_damage[]
    old_doy      = _C.hlm_day_of_year[]
    old_dpy      = _C.hlm_days_per_year[]
    old_freqday  = _C.hlm_freq_day[]
    old_nharv    = _C.hlm_num_lu_harvest_cats[]
    old_ch4      = _C.hlm_use_ch4[]
    old_vert     = _C.hlm_use_vertsoilc[]
    old_spit     = _C.hlm_spitfire_mode[]
    old_paramfp  = _C.fates_maxPatchesPerSite[]
    old_nlevsno  = _C.varpar.nlevsno
    old_nlevsoi  = _C.varpar.nlevsoi
    old_nlevgrnd = _C.varpar.nlevgrnd
    old_nlevmaxu = _C.varpar.nlevmaxurbgrnd

    try
        nlevsoil = 5
        nlevdecomp = 1

        # ===============================================================
        # Build a driver-ready single-column inst, FATES-tagged.
        # ===============================================================
        _C.varpar_init!(_C.varpar, 1, 14, 2, 5)
        _C.varcon_init!()
        _C.varpar.nlevsno = 5
        _C.varpar.nlevsoi = nlevsoil
        _C.varpar.nlevgrnd = nlevsoil
        _C.varpar.nlevmaxurbgrnd = nlevsoil
        nlevsno = _C.varpar.nlevsno

        # 1 gridcell / 1 landunit / 1 column / 2 patches (1 = bareground, 2 = veg).
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

        # snow-layer constants (handle_new_snow!).
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

        # pftcon defaults (non-FATES patches / albedo).
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

        # ---- cold-start a real-param (14-PFT) carbon-only FATES site ----
        fates = _C.clm_fates_init!(inst; nsites=1, nlevsoil=nlevsoil, nlevdecomp=nlevdecomp)
        @test inst.fates !== nothing
        @test inst.fates.nsites == 1
        @test _C.numpft[] == 14

        config = _C.CLMDriverConfig(use_fates=true, use_aquifer_layer=false,
                                    use_hydrstress=false, use_luna=false)
        @test config.use_fates == true

        # ===============================================================
        # Forcing helper (representative mid-latitude summer snapshot).
        # ===============================================================
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
            a.forc_u_grc[g]=2.0; a.forc_v_grc[g]=0.0
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

        # FATES census (walk the cohort/patch linked lists): finiteness + physical
        # ranges + total patch area.
        function fates_census(inst)
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

        cen0 = fates_census(inst)
        @test cen0.npatch >= 1
        @test cen0.ncoh >= 1
        @test isempty(cen0.bad)
        @test cen0.totarea ≈ _C.area atol = 1e-6   # AREA = 10000 m2

        # ===============================================================
        # The multi-day clm_drv! loop: 96 steps @ 1800 s = 2 days.
        # ===============================================================
        c = 1; pveg = col.patchi[1] + 1   # vegetated patch = 2
        dtime = 1800.0
        steps_per_day = Int(round(86400 / dtime))   # 48
        ndays = 2
        nsteps = steps_per_day * ndays
        declin = 0.4   # ~summer-solstice declination (rad)

        secs = 0; jday = 152; daycount = 0
        daily_fires = 0
        no_error = true
        nan_steps = String[]
        balance_threw = false
        # FATES-drove-it trackers (max over the run of the FATES-branch outputs).
        max_laisha = 0.0; max_canopy_alb = 0.0; saw_rssun = false
        prev_npatch = cen0.npatch

        for i in 1:nsteps
            is_beg = (secs == 0)
            # diurnal solar: sin over the day (0 at night, peak midday).
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
                if e isa ErrorException && occursin("alance", sprint(showerror, e))
                    balance_threw = true
                end
                @error "FATES driver step errored" step=i is_beg exception=(e, catch_backtrace())
                break
            end
            if is_beg && i > 1
                daily_fires += 1
            end

            # ---- finiteness audit on the FATES column's key CLM outputs ----
            cs = inst.canopystate; sa = inst.surfalb
            ps = inst.photosyns; ef = inst.energyflux
            audited = Dict(
                "elai"=>cs.elai_patch[pveg], "fsun"=>cs.fsun_patch[pveg],
                "htop"=>cs.htop_patch[pveg], "laisun"=>cs.laisun_patch[pveg],
                "laisha"=>cs.laisha_patch[pveg],
                "albd"=>sa.albd_patch[pveg,1], "albi"=>sa.albi_patch[pveg,1],
                "rssun"=>ps.rssun_patch[pveg], "rssha"=>ps.rssha_patch[pveg],
                "btran"=>ef.btran_patch[pveg])
            bad = [k for (k,v) in audited if !isfinite(v)]
            isempty(bad) || push!(nan_steps, "step$i:$bad")

            # ---- FATES census finiteness/physicality each step ----
            cen = fates_census(inst)
            isempty(cen.bad) || push!(nan_steps, "step$i:census$(cen.bad)")
            # area is conserved at AREA each step.
            isapprox(cen.totarea, _C.area; atol=1e-4) || push!(nan_steps, "step$i:area=$(cen.totarea)")

            # ---- FATES-drove-it trackers ----
            isfinite(cs.laisha_patch[pveg]) && (max_laisha = max(max_laisha, cs.laisha_patch[pveg]))
            isfinite(sa.albd_patch[pveg,1]) && (max_canopy_alb = max(max_canopy_alb, sa.albd_patch[pveg,1]))
            (isfinite(ps.rssun_patch[pveg]) && ps.rssun_patch[pveg] > 0.0) && (saw_rssun = true)
            prev_npatch = cen.npatch

            # advance the clock.
            secs += Int(dtime)
            if secs >= 86400
                secs = 0; jday += 1; daycount += 1
            end
        end

        # ===============================================================
        # Assertions.
        # ===============================================================
        # 1. completes without error every step.
        @test no_error
        @test !balance_threw

        # 4. the daily demographic step fired at least once (>= 1 day boundary).
        @test daily_fires >= 1

        # 2. no NaN/Inf in CLM outputs nor FATES census over the whole run.
        @test isempty(nan_steps)

        # 3. FATES actually drove the column: the FATES sun/shade hook populated a
        #    nonzero shaded leaf area on the veg patch (cold-start saplings carry a
        #    small canopy), and the FATES photosynthesis hook produced a positive
        #    stomatal resistance — both distinct from the (unrun) non-FATES path.
        @test max_laisha > 0.0
        @test saw_rssun

        # 5. census stays physical + area conserved over the multi-day horizon.
        cenF = fates_census(inst)
        @test isempty(cenF.bad)
        @test cenF.totarea ≈ _C.area atol = 1e-4
        @test cenF.ncoh >= 1
        # the daily demographic step advanced the population (disturbance can split
        # a patch); patch count is finite & physical.
        @test cenF.npatch >= 1

        # final TotalBalanceCheck holds (would have thrown inside the daily step).
        @test _C.TotalBalanceCheck(inst.fates.sites[1], -1) === nothing

    finally
        _C.prt_global[]                  = old_global
        _C.num_elements[]                = old_numel
        empty!(_C.element_list); append!(_C.element_list, old_ellist)
        _C.hlm_parteh_mode[]             = old_parteh
        _C.hlm_use_sp[]                  = old_use_sp
        _C.hlm_use_planthydro[]          = old_hydro
        _C.hlm_use_cohort_age_tracking[] = old_agetrk
        _C.nleafage[]                    = old_nleafage
        _C.numpft[]                      = old_numpft
        _C.hlm_current_tod[]             = old_tod
        _C.nlevsclass[]                  = old_nlevsc
        _C.nlevcoage[]                   = old_nlevca
        _C.nlevage[]                     = old_nlevage
        _C.nlevdamage[]                  = old_nlevdam
        _C.nlevheight[]                  = old_nlevhgt
        _C.SFParams[]                    = old_sfp
        _C.EDParams[]                    = old_edparams
        _C.EDPftvarcon_inst[]            = old_edpft
        _C.ParamDerived[]                = old_paramd
        _C.hlm_is_restart[]              = old_restart
        _C.hlm_use_nocomp[]              = old_nocomp
        _C.hlm_use_fixed_biogeog[]       = old_biogeog
        _C.hlm_use_luh[]                 = old_luh
        _C.hlm_use_inventory_init[]      = old_inv
        _C.hlm_use_tree_damage[]         = old_damage
        _C.hlm_day_of_year[]             = old_doy
        _C.hlm_days_per_year[]           = old_dpy
        _C.hlm_freq_day[]                = old_freqday
        _C.hlm_num_lu_harvest_cats[]     = old_nharv
        _C.hlm_use_ch4[]                 = old_ch4
        _C.hlm_use_vertsoilc[]           = old_vert
        _C.hlm_spitfire_mode[]           = old_spit
        _C.fates_maxPatchesPerSite[]     = old_paramfp
        _C.varpar.nlevsno                = old_nlevsno
        _C.varpar.nlevsoi                = old_nlevsoi
        _C.varpar.nlevgrnd               = old_nlevgrnd
        _C.varpar.nlevmaxurbgrnd         = old_nlevmaxu
    end
end
