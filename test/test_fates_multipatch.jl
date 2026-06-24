# test_fates_multipatch.jl
# FATES live-driver wiring — multi-veg-patch / multi-site coupling generalization.
#
# The single-veg-patch MVP mapped a FATES column to exactly ONE vegetated patch via
# `p = col.patchi[c] + 1`. The Fortran host instead walks ALL the column's vegetated
# patches (`p = ifp + col.patchi(c)`, ifp = 1..site%youngest_patch%patchno) and
# rebuilds the HLM patch weights after dynamics (`setFilters`). This test exercises
# the generalized Julia coupling:
#
#   PART A (multi-veg-patch, deterministic): cold-start a single FATES site, then
#     splice a SECOND vegetated patch into its age-ordered linked list (halving the
#     area, mirroring a disturbance split). Assert:
#       * fates_veg_patches walks BOTH patches → (ifp,p) = (1,patchi+1),(2,patchi+2);
#       * update_hlm_dynamics! + the bc pack/unpack touch BOTH patch slots (the bc
#         *_pa arrays carry per-patch values, not just slot 1);
#       * fates_set_filters! sets pch.wtcol for BOTH veg patches from
#         canopy_fraction_pa (+ the bare-ground remainder), and the column weights
#         sum to ~1; both veg patches become active.
#
#   PART B (multi-site): cold-start nsites=2 on two FATES-tagged columns and loop the
#     real clm_drv! for >= 1 day. Assert both sites run finite end-to-end (no NaN in
#     either column's CLM outputs nor either site's FATES census) and the daily step
#     advances both.
#
# clm_fates_init! mutates FATES module-global control Refs + parameter tables, so
# save/restore them (same set as test_fates_driver_run.jl).

using Test
using CLM
const _C = CLM

@testset "FATES multi-veg-patch / multi-site coupling" begin
    # ---- save FATES + varpar globals (mirrors test_fates_driver_run.jl) ----
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
    old_varctl_fates = _C.varctl.use_fates

    try
        # allocate the FATES-only patch fields (is_veg/is_bareground/wt_ed) in
        # patch_init! so the setFilters rebuild can populate them (a real FATES run
        # sets this; the default suite leaves it false → fields stay empty).
        _C.varctl.use_fates = true

        nlevsoil = 5
        nlevdecomp = 1

        # =================================================================
        # PART A — multi-veg-patch ifp walk + setFilters weight rebuild.
        #
        # A focused unit test of the two new mapping primitives, exercised on a
        # FATES site whose patch linked list is grown to TWO vegetated patches (a
        # disturbance-style split). The site's canopy_fraction_pa is set directly
        # (the cold-start cohort crown areas are NaN until a radiation pass runs, so
        # we don't rely on update_hlm_dynamics! here — that path is covered live by
        # PART B). This isolates: (a) the oldest->younger ifp walk → (ifp,p) pairs,
        # (b) the patchi+ifp HLM mapping, (c) the fates_set_filters! weight rebuild.
        # =================================================================
        @testset "multi-veg-patch ifp walk + setFilters" begin
            _C.varpar_init!(_C.varpar, 1, 14, 2, 5)
            _C.varcon_init!()
            _C.varpar.nlevsno = 5
            _C.varpar.nlevsoi = nlevsoil
            _C.varpar.nlevgrnd = nlevsoil
            _C.varpar.nlevmaxurbgrnd = nlevsoil

            # 1 gridcell / 1 landunit / 1 column / 3 patches:
            #   patch 1 = bareground (col.patchi[c]+0),
            #   patches 2,3 = two vegetated FATES patches (+1, +2).
            ng, nl, nc, np = 1, 1, 1, 3
            inst = _C.CLMInstances()
            _C.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np,
                             nlevdecomp_full=_C.varpar.nlevdecomp_full)

            col = inst.column
            col.landunit[1] = 1; col.gridcell[1] = 1; col.snl[1] = 0
            col.patchi[1] = 1; col.patchf[1] = 3; col.npatches[1] = 3
            col.nbedrock[1] = _C.varpar.nlevsoi
            col.is_fates[1] = true

            lun = inst.landunit
            lun.itype[1] = _C.ISTSOIL; lun.urbpoi[1] = false
            lun.lakpoi[1] = false; lun.active[1] = true

            pch = inst.patch
            for p in 1:np
                pch.active[p] = true; pch.landunit[p] = 1; pch.gridcell[p] = 1
                pch.column[p] = 1; pch.itype[p] = 1; pch.is_fates[p] = true
            end

            fates = _C.clm_fates_init!(inst; nsites=1, nlevsoil=nlevsoil,
                                       nlevdecomp=nlevdecomp)
            site = fates.sites[1]

            # ---- before split: exactly ONE vegetated patch (NBG cold start) ----
            walk0 = _C.fates_veg_patches(site, 1, col)
            @test length(walk0) == 1
            @test walk0[1] == (1, col.patchi[1] + 1)   # (ifp=1, p=patchi+1)

            # ---- grow the patch list to TWO vegetated patches (split) ----
            # Halve the existing patch's area and link a sibling veg patch as the
            # younger patch. (We test the mapping/weight logic, so the sibling needs
            # no cohorts — its canopy_fraction is supplied below.)
            p1 = site.oldest_patch
            @test p1 !== nothing
            half = p1.area / 2
            p1.area = half

            p2 = _C.fates_patch_type()
            _C.Create!(p2, 0.0, half, p1.land_use_label, p1.nocomp_pft_label,
                       _C.num_swb, _C.numpft[], site.nlevsoil, _C.hlm_current_tod[],
                       _C.ed_params().regeneration_model)
            p2.older = p1; p2.younger = nothing
            p1.younger = p2
            site.youngest_patch = p2
            _C.set_patchno(site)

            # both patches are vegetated (non-bareground); patchno of youngest = 2.
            @test p1.nocomp_pft_label != _C.nocomp_bareground
            @test p2.nocomp_pft_label != _C.nocomp_bareground
            @test site.youngest_patch.patchno == 2

            # ---- the ifp walk now yields BOTH vegetated patches ----
            walk = _C.fates_veg_patches(site, 1, col)
            @test length(walk) == 2
            @test walk[1] == (1, col.patchi[1] + 1)
            @test walk[2] == (2, col.patchi[1] + 2)

            # ---- drive the bc *_pa arrays per patch (distinct values) + unpack ----
            bc_out = fates.bc_out[1]
            bc_out.htop_pa[1] = 12.0; bc_out.htop_pa[2] = 6.0
            bc_out.hbot_pa[1] = 0.2;  bc_out.hbot_pa[2] = 0.1
            bc_out.elai_pa[1] = 3.0;  bc_out.elai_pa[2] = 1.5
            bc_out.esai_pa[1] = 0.5;  bc_out.esai_pa[2] = 0.25
            bc_out.z0m_pa[1]  = 1.0;  bc_out.z0m_pa[2]  = 0.5
            bc_out.displa_pa[1] = 8.0; bc_out.displa_pa[2] = 4.0
            bc_out.dleaf_pa[1]  = 0.04; bc_out.dleaf_pa[2] = 0.04
            # canopy fractions: each veg patch covers 0.4 of the site (bg = 0.2).
            bc_out.canopy_fraction_pa[1] = 0.4
            bc_out.canopy_fraction_pa[2] = 0.4

            for (ifp, p) in walk
                _C.fates_unpack_bcout_canopy_structure!(inst; s=1, c=1, p=p, ifp=ifp)
            end
            cs = inst.canopystate
            # the two veg patches got DISTINCT slots (not both written to patchi+1).
            @test cs.htop_patch[col.patchi[1] + 1] == bc_out.htop_pa[1]   # 12.0
            @test cs.htop_patch[col.patchi[1] + 2] == bc_out.htop_pa[2]   # 6.0
            @test cs.elai_patch[col.patchi[1] + 1] == bc_out.elai_pa[1]   # 3.0
            @test cs.elai_patch[col.patchi[1] + 2] == bc_out.elai_pa[2]   # 1.5

            # ---- setFilters rebuild: HLM patch weights from canopy fractions ----
            _C.fates_set_filters!(inst; c=1, s=1)
            bg = col.patchi[1]                 # bareground HLM patch
            pv1 = col.patchi[1] + 1
            pv2 = col.patchi[1] + 2
            sum_cf = bc_out.canopy_fraction_pa[1] + bc_out.canopy_fraction_pa[2]  # 0.8
            # bareground weight = 1 - sum(canopy fractions); veg weights = fractions.
            @test pch.wtcol[bg]  ≈ max(0.0, 1.0 - sum_cf)   # 0.2
            @test pch.wtcol[pv1] ≈ bc_out.canopy_fraction_pa[1]   # 0.4
            @test pch.wtcol[pv2] ≈ bc_out.canopy_fraction_pa[2]   # 0.4
            # the column patch weights (bareground + the two veg patches) sum to 1.
            @test pch.wtcol[bg] + pch.wtcol[pv1] + pch.wtcol[pv2] ≈ 1.0 atol = 1e-12
            # both veg patches flagged is_veg + active; bareground flagged.
            @test pch.is_veg[pv1] && pch.is_veg[pv2]
            @test pch.is_bareground[bg]
            @test pch.active[pv1] && pch.active[pv2]
            @test pch.wt_ed[pv1] ≈ bc_out.canopy_fraction_pa[1]
            @test pch.wt_ed[pv2] ≈ bc_out.canopy_fraction_pa[2]

            # ---- a MERGE (fewer patches) zeroes + deactivates the freed slot ----
            # Re-link to a single veg patch (p1 only) and rebuild filters: the now
            # unused veg slot (patchi+2) must be zeroed + deactivated.
            p1.younger = nothing
            site.youngest_patch = p1
            _C.set_patchno(site)
            @test length(_C.fates_veg_patches(site, 1, col)) == 1
            _C.fates_set_filters!(inst; c=1, s=1)
            @test pch.wtcol[pv1] ≈ bc_out.canopy_fraction_pa[1]   # still 0.4
            @test pch.wtcol[pv2] == 0.0                            # freed slot
            @test pch.active[pv2] == false
            @test pch.wtcol[bg] ≈ 1.0 - bc_out.canopy_fraction_pa[1]  # 0.6
        end

        # =================================================================
        # PART B — multi-site (nsites=2) clm_drv! integrated run, >= 1 day.
        # =================================================================
        @testset "multi-site clm_drv! run finite" begin
            _C.varpar_init!(_C.varpar, 1, 14, 2, 5)
            _C.varcon_init!()
            _C.varpar.nlevsno = 5
            _C.varpar.nlevsoi = nlevsoil
            _C.varpar.nlevgrnd = nlevsoil
            _C.varpar.nlevmaxurbgrnd = nlevsoil
            nlevsno = _C.varpar.nlevsno

            # 2 gridcells / 2 landunits / 2 columns / 6 patches.
            # Each column reserves bareground (patchi+0) + TWO veg slots (patchi+1,+2)
            # so an organic disturbance split (one cold-start patch → two over a day)
            # has an HLM slot to land in.
            #   col 1: patches 1(bg),2(veg),3(veg);  col 2: patches 4(bg),5,6(veg).
            ng, nl, nc, np = 2, 2, 2, 6
            inst = _C.CLMInstances()
            _C.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np,
                             nlevdecomp_full=_C.varpar.nlevdecomp_full)

            col = inst.column
            for (cc, (g, pi, pf)) in enumerate(((1,1,3),(2,4,6)))
                col.landunit[cc] = cc; col.gridcell[cc] = g; col.snl[cc] = 0
                col.patchi[cc] = pi; col.patchf[cc] = pf; col.npatches[cc] = 3
                col.nbedrock[cc] = _C.varpar.nlevsoi
                col.is_fates[cc] = true
            end

            lun = inst.landunit
            for l in 1:nl
                lun.itype[l] = _C.ISTSOIL; lun.urbpoi[l] = false
                lun.lakpoi[l] = false; lun.active[l] = true
            end

            pch = inst.patch
            # patch -> (gridcell, column): 1,2,3 -> (1,1);  4,5,6 -> (2,2).
            pmap = ((1,1),(1,1),(1,1),(2,2),(2,2),(2,2))
            for p in 1:np
                g, cc = pmap[p]
                pch.active[p] = true; pch.landunit[p] = cc; pch.gridcell[p] = g
                pch.column[p] = cc; pch.itype[p] = 1; pch.is_fates[p] = true
                # bareground slot (patchi) weight 0, first veg slot (patchi+1) weight 1.
                pch.wtcol[p] = (p == col.patchi[cc] + 1 ? 1.0 : 0.0)
            end

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

            # ---- cold-start TWO real-param (14-PFT) carbon-only FATES sites ----
            fates = _C.clm_fates_init!(inst; nsites=2, nlevsoil=nlevsoil,
                                       nlevdecomp=nlevdecomp)
            @test inst.fates !== nothing
            @test inst.fates.nsites == 2

            config = _C.CLMDriverConfig(use_fates=true, use_aquifer_layer=false,
                                        use_hydrstress=false, use_luna=false)

            function set_forcing!(inst, sunfrac)
                a = inst.atm2lnd
                for g in 1:ng
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
                end
                for cc in 1:nc; inst.topo.topo_col[cc]=150.0; end
                return nothing
            end

            function fates_census(site)
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

            # both sites cold-start finite.
            for s in 1:2
                cen = fates_census(fates.sites[s])
                @test cen.npatch >= 1
                @test cen.ncoh >= 1
                @test isempty(cen.bad)
                @test cen.totarea ≈ _C.area atol = 1e-6
            end

            dtime = 1800.0
            steps_per_day = Int(round(86400 / dtime))
            nsteps = steps_per_day + 1   # just past 1 day boundary
            declin = 0.4
            secs = 0; jday = 152; daycount = 0
            no_error = true; nan_steps = String[]; daily_fires = 0
            # veg patch (patchi+1) per column.
            vegp = (col.patchi[1] + 1, col.patchi[2] + 1)
            max_npatch = (0, 0)   # max FATES patch count seen per site
            weight_bad = String[]

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
                    @error "multi-site FATES step errored" step=i exception=(e, catch_backtrace())
                    break
                end
                (is_beg && i > 1) && (daily_fires += 1)

                # finiteness audit on BOTH FATES columns' key CLM outputs + census.
                cs = inst.canopystate; sa = inst.surfalb
                ps = inst.photosyns; ef = inst.energyflux
                for cc in 1:2
                    pveg = vegp[cc]
                    vals = (cs.elai_patch[pveg], cs.laisha_patch[pveg],
                            sa.albd_patch[pveg,1], ps.rssun_patch[pveg],
                            ef.btran_patch[pveg])
                    all(isfinite, vals) || push!(nan_steps, "step$i:col$cc")
                end
                for s in 1:2
                    cen = fates_census(fates.sites[s])
                    isempty(cen.bad) || push!(nan_steps, "step$i:site$(s)census")
                    isapprox(cen.totarea, _C.area; atol=1e-4) ||
                        push!(nan_steps, "step$i:site$(s)area")
                end

                # track max FATES patch count + verify the HLM weight rebuild keeps
                # each FATES column's patch weights (over its reserved range) summing
                # to ~1 after the daily setFilters rebuild.
                for s in 1:2
                    cc = s
                    yp = fates.sites[s].youngest_patch
                    npq = yp === nothing ? 0 : yp.patchno
                    max_npatch = (s == 1 ? (max(max_npatch[1], npq), max_npatch[2]) :
                                          (max_npatch[1], max(max_npatch[2], npq)))
                end
                if is_beg && i > 1
                    # after a daily rebuild, the column patch weights sum to ~1.
                    for cc in 1:2
                        wsum = 0.0
                        for p in col.patchi[cc]:col.patchf[cc]; wsum += pch.wtcol[p]; end
                        isapprox(wsum, 1.0; atol=1e-6) ||
                            push!(weight_bad, "col$cc:wsum=$wsum")
                    end
                end

                secs += Int(dtime)
                if secs >= 86400; secs = 0; jday += 1; daycount += 1; end
            end

            @test no_error
            @test isempty(nan_steps)
            @test isempty(weight_bad)   # column weights sum to 1 post-rebuild.
            @test daily_fires >= 1    # the daily step fired on the day boundary.

            # both sites still physical + area-conserved, TotalBalanceCheck holds.
            for s in 1:2
                cen = fates_census(fates.sites[s])
                @test isempty(cen.bad)
                @test cen.totarea ≈ _C.area atol = 1e-4
                @test cen.ncoh >= 1
                @test _C.TotalBalanceCheck(fates.sites[s], -1) === nothing
            end
        end

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
        _C.varctl.use_fates              = old_varctl_fates
    end
end
