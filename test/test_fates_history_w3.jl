# test_fates_history_w3.jl
# FATES history buffer-fills, Wave 3 (B18-followup H-series, FINAL wave):
#   * update_history_nutrflux!     — per-cohort CNP flux diagnostics (N/P uptake,
#     efflux, demand, fixation, excess-resp) -> *_si and *_scpf bins.
#   * the FIRE-dimensioned groups   — dyn1 site-level SPITFIRE scalars + patch-area
#     weighted fuel/intensity (NESTEROV/FDI/ROS/FUEL_*/SUM_FUEL/...), and dyn2's
#     fire-by-age + per-fuel-class diagnostics.
#   * the TREE-DAMAGE cross-tab     — dyn2 cohort + canopy/understory + site-level
#     *_si_cdpf cross-tab (damage x size x pft) and dyn1 CROWNAREA_*_DAMAGE.
#   * update_history_hydraulics!   — si_hydr site water + root-weighted soil means,
#     plus the size x pft co_hydr cohort-mean / flux scpf diagnostics.
#
# These groups read state (CNP fluxes, a `fuel` object, fire-weather + fire-mort
# arrays, damaged cohorts, si_hydr/co_hydr) that is unset in the default carbon-only
# cold-start. Per the validated W2 approach we CONSTRUCT representative state on the
# cohorts/patches/site, run the relevant update_history_*, then assert the ih_*
# handles populate FINITE / physical values with CORRECT binning (the cohort's flux
# lands in its size x pft / damage-class bin) and that aggregation/conservation holds.
#
# clm_fates_init! mutates FATES module-global control Refs; save/restore them so the
# surrounding default suite stays unperturbed.

using Test
using CLM

@testset "FATES history W3 (nutrflux + fire + damage + hydraulics)" begin
    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_elpos    = copy(CLM.element_pos)
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_sp   = CLM.hlm_use_sp[]
    old_st3      = CLM.hlm_use_ed_st3[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_agetrk   = CLM.hlm_use_cohort_age_tracking[]
    old_nleafage = CLM.nleafage[]
    old_numpft   = CLM.numpft[]
    old_tod      = CLM.hlm_current_tod[]
    old_nlevsc   = CLM.nlevsclass[]
    old_nlevca   = CLM.nlevcoage[]
    old_nlevage  = CLM.nlevage[]
    old_nlevdam  = CLM.nlevdamage[]
    old_nlevht   = CLM.nlevheight[]
    old_sfp      = CLM.SFParams[]
    old_restart  = CLM.hlm_is_restart[]
    old_nocomp   = CLM.hlm_use_nocomp[]
    old_biogeog  = CLM.hlm_use_fixed_biogeog[]
    old_luh      = CLM.hlm_use_luh[]
    old_inv      = CLM.hlm_use_inventory_init[]
    old_damage   = CLM.hlm_use_tree_damage[]
    old_doy      = CLM.hlm_day_of_year[]
    old_dpy      = CLM.hlm_days_per_year[]
    old_freqday  = CLM.hlm_freq_day[]
    old_nharv    = CLM.hlm_num_lu_harvest_cats[]
    old_ch4      = CLM.hlm_use_ch4[]
    old_vert     = CLM.hlm_use_vertsoilc[]
    old_spit     = CLM.hlm_spitfire_mode[]
    old_hdyn     = CLM.hlm_hist_level_dynam[]
    old_hhfr     = CLM.hlm_hist_level_hifrq[]
    old_paramfp  = CLM.fates_maxPatchesPerSite[]
    old_modelday = CLM.hlm_model_day[]
    old_nlevsno  = CLM.varpar.nlevsno
    old_nlevsoi  = CLM.varpar.nlevsoi
    old_nlevgrnd = CLM.varpar.nlevgrnd
    old_nlevmaxu = CLM.varpar.nlevmaxurbgrnd

    try
        nlevsoil = 3
        nlevdecomp = 1
        ng, nl, nc, np = 1, 1, 1, 2
        CLM.varpar.nlevsno = 5
        CLM.varpar.nlevsoi = nlevsoil
        CLM.varpar.nlevgrnd = nlevsoil
        CLM.varpar.nlevmaxurbgrnd = nlevsoil

        inst = CLM.CLMInstances()
        CLM.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np)

        col = inst.column
        col.gridcell[1] = 1
        col.patchi[1]   = 1
        col.npatches[1] = 2
        col.is_fates[1] = true
        inst.patch.column[1] = 1
        inst.patch.column[2] = 1
        inst.patch.is_fates[1] = true
        inst.patch.is_fates[2] = true

        cveg = 1

        # ---- cold-started carbon-only FATES site at level-2 history ----
        fates = CLM.clm_fates_init!(inst; nsites = 1, nlevsoil = nlevsoil,
                                    nlevdecomp = nlevdecomp, hist_dimlevel = 2)
        @test inst.fates !== nothing
        hist = inst.fates.hist
        @test CLM.hlm_hist_level_dynam[] == 2
        @test CLM.hlm_hist_level_hifrq[] == 2

        # representative soil-column state for the daily step
        ss   = inst.soilstate
        temp = inst.temperature
        wdb  = inst.water.waterdiagnosticbulk_inst
        joff = CLM.varpar.nlevsno
        for j in 1:nlevsoil
            temp.t_soisno_col[cveg, joff + j] = 290.0
            wdb.h2osoi_liqvol_col[cveg, j]    = 0.30
            ss.watsat_col[cveg, j]            = 0.45
            ss.eff_porosity_col[cveg, j]      = 0.45
            ss.sucsat_col[cveg, j]            = 100.0
            ss.bsw_col[cveg, j]               = 6.0
        end

        site = inst.fates.sites[1]
        io   = site.h_gid
        hbuf(sym) = hist.hvars[hist.ih[sym]].r81d[io]
        hrow(sym) = @view hist.hvars[hist.ih[sym]].r82d[io, :]

        # advance one daily step so cohorts are non-new + class indices set
        CLM.fates_daily_dynamics_step!(inst; nlevsoil=nlevsoil,
                                       nlevdecomp=nlevdecomp, use_fates_bgc=false)

        # gather the live cohorts (with prt) + their class indices for binning checks
        coh = CLM.fates_cohort_type[]
        pp = site.oldest_patch
        while pp !== nothing
            cc = pp.shortest
            while cc !== nothing
                cc.prt !== nothing && push!(coh, cc)
                cc = cc.taller
            end
            pp = pp.younger
        end
        @test !isempty(coh)

        # ==============================================================
        # 1. NUTRFLUX — construct CNP flux state, run update_history_nutrflux!
        # ==============================================================
        # The cold start is carbon-only (num_elements==1), so the N/P diagnostics are
        # `:has_n`/`:has_p`-gated off in the registry. Splice nitrogen+phosphorus into
        # element_list and REBUILD the history interface so the N/P (+ scpf) handles
        # register, then set representative daily fluxes on each cohort.
        CLM.num_elements[] = 3
        empty!(CLM.element_list)
        append!(CLM.element_list, [CLM.carbon12_element, CLM.nitrogen_element, CLM.phosphorus_element])

        histn = CLM._build_fates_history_interface(1)
        @test haskey(histn.ih, :ih_nh4uptake_si)
        @test haskey(histn.ih, :ih_nh4uptake_scpf)
        @test haskey(histn.ih, :ih_excess_resp_si)
        ion = site.h_gid
        hbufn(sym) = histn.hvars[histn.ih[sym]].r81d[ion]
        hrown(sym) = @view histn.hvars[histn.ih[sym]].r82d[ion, :]

        for cc in coh
            cc.isnew            = false
            cc.resp_excess      = 2.0e-7
            cc.daily_nh4_uptake = 3.0e-7
            cc.daily_no3_uptake = 1.0e-7
            cc.sym_nfix_daily   = 5.0e-8
            cc.daily_n_efflux   = 2.0e-8
            cc.daily_n_demand   = 6.0e-7
            cc.daily_p_gain     = 4.0e-8
            cc.daily_p_efflux   = 1.0e-8
            cc.daily_p_demand   = 9.0e-8
        end

        CLM.update_history_nutrflux!(histn, site)

        # excess-resp accumulates with uconv = n*ha_per_m2*days_per_sec
        sumn = sum(c.n for c in coh)
        @test isfinite(hbufn(:ih_excess_resp_si))
        @test hbufn(:ih_excess_resp_si) > 0.0
        @test hbufn(:ih_excess_resp_si) ≈ 2.0e-7 * sumn * CLM.ha_per_m2 * CLM.days_per_sec rtol=1e-10
        # N/P site-level uptake/demand finite + positive
        for sym in (:ih_nh4uptake_si, :ih_no3uptake_si, :ih_nfix_si, :ih_ndemand_si,
                    :ih_puptake_si, :ih_pdemand_si)
            @test isfinite(hbufn(sym)) && hbufn(sym) > 0.0
        end
        # ammonium uptake equals Σ n * flux * conv
        @test hbufn(:ih_nh4uptake_si) ≈ 3.0e-7 * sumn * CLM.ha_per_m2 * CLM.days_per_sec rtol=1e-10
        # scpf disaggregation: the cohort's flux lands in its size x pft bin, and the
        # site total equals the sum over scpf bins (conservation)
        nh4_scpf = hrown(:ih_nh4uptake_scpf)
        @test all(isfinite, nh4_scpf)
        @test nh4_scpf[coh[1].size_by_pft_class] > 0.0
        @test sum(nh4_scpf) ≈ hbufn(:ih_nh4uptake_si) rtol=1e-10

        # restore carbon-only element registry for the remaining groups
        CLM.num_elements[] = 1
        empty!(CLM.element_list); push!(CLM.element_list, CLM.carbon12_element)

        # ==============================================================
        # 2. FIRE — construct a fuel object + fire-weather + per-patch fire state
        #    + site fire-mortality arrays, then run dyn1 + dyn2.
        # ==============================================================
        site.fireWeather = CLM.nesterov_index()
        site.fireWeather.fire_weather_index  = 12.5
        site.fireWeather.effective_windspeed = 60.0    # m/min
        site.fdi          = 0.4
        site.NF_successful = 0.02

        pp = site.oldest_patch
        while pp !== nothing
            pp.ros_front  = 1.2          # m/min
            pp.fi         = 150.0        # kW/m
            pp.tfc_ros    = 0.05         # kgC/m2/day
            pp.frac_burnt = 0.1          # frac/day
            if pp.fuel === nothing
                pp.fuel = CLM.fuel_type()
                CLM.init_fuel!(pp.fuel)
            end
            f = pp.fuel
            f.bulk_density_notrunks     = 12.0
            f.average_moisture_notrunks = 0.15
            f.SAV_notrunks              = 60.0
            f.MEF_notrunks              = 0.25
            f.non_trunk_loading         = 0.8
            for k in 1:CLM.num_fuel_classes
                f.frac_loading[k]       = 1.0 / CLM.num_fuel_classes
                f.effective_moisture[k] = 0.12
                f.frac_burnt[k]         = 0.3
            end
            pp = pp.younger
        end

        # set site fire-mortality scpf arrays (allocated by zero_site!) so the dyn2
        # crown/cambial fire-mort by scpf is exercised
        if !isempty(site.fmort_rate_crown)
            fill!(site.fmort_rate_crown, 1.0)
            fill!(site.fmort_rate_cambial, 0.5)
            fill!(site.fmort_rate_canopy, 2.0)
            fill!(site.fmort_rate_ustory, 1.0)
        end

        CLM.update_history_dyn1!(hist, 1, fates.nsites, fates.sites, fates.bc_in)

        # -- site fire scalars finite + physical --
        @test hbuf(:ih_nesterov_fire_danger_si) ≈ 12.5
        @test hbuf(:ih_fire_fdi_si) ≈ 0.4
        @test hbuf(:ih_effect_wspeed_si) ≈ 60.0 / CLM.sec_per_min rtol=1e-12
        @test hbuf(:ih_fire_nignitions_si) ≈ 0.02 / CLM.m2_per_km2 / CLM.sec_per_day rtol=1e-12
        # -- patch-area-weighted fuel/intensity (positive, finite) --
        for sym in (:ih_spitfire_ros_si, :ih_tfc_ros_si, :ih_fire_intensity_si,
                    :ih_fire_area_si, :ih_fire_fuel_bulkd_si, :ih_fire_fuel_eff_moist_si,
                    :ih_fire_fuel_sav_si, :ih_fire_fuel_mef_si, :ih_sum_fuel_si,
                    :ih_fire_intensity_area_product_si)
            @test isfinite(hbuf(sym)) && hbuf(sym) > 0.0
        end
        # fire intensity = Σ fi*area/AREA*J_per_kJ ; patches sum to total area
        site_area_frac = 0.0
        pp = site.oldest_patch
        while pp !== nothing; site_area_frac += pp.area * CLM.AREA_INV; pp = pp.younger; end
        @test hbuf(:ih_fire_intensity_si) ≈ 150.0 * site_area_frac * CLM.J_per_kJ rtol=1e-10
        @test hbuf(:ih_fire_fuel_bulkd_si) ≈ 12.0 * site_area_frac rtol=1e-10

        CLM.update_history_dyn2!(hist, 1, fates.nsites, fates.sites, fates.bc_in)

        # -- fire-by-age: a populated age bin --
        area_burnt_age  = hrow(:ih_area_burnt_si_age)
        fire_intens_age = hrow(:ih_fire_intensity_si_age)
        fire_fuel_age   = hrow(:ih_fire_sum_fuel_si_age)
        @test all(isfinite, area_burnt_age)
        @test sum(area_burnt_age) > 0.0
        @test sum(fire_intens_age) > 0.0
        @test sum(fire_fuel_age) > 0.0
        # -- per-fuel-class diagnostics: every fuel class carries finite values --
        litter_moist_fuel = hrow(:ih_litter_moisture_si_fuel)
        fuel_amount_fuel  = hrow(:ih_fuel_amount_si_fuel)
        burnt_frac_fuel   = hrow(:ih_burnt_frac_litter_si_fuel)
        @test all(isfinite, litter_moist_fuel)
        @test all(>(0.0), litter_moist_fuel[1:CLM.num_fuel_classes])
        @test all(>(0.0), fuel_amount_fuel[1:CLM.num_fuel_classes])
        @test all(>(0.0), burnt_frac_fuel[1:CLM.num_fuel_classes])
        # total fuel amount over classes == sum_fuel (frac_loading sums to 1)
        @test sum(fuel_amount_fuel) ≈ 0.8 * site_area_frac rtol=1e-10
        # crown/cambial fire mortality by scpf is populated
        @test sum(hrow(:ih_crownfiremort_si_scpf)) > 0.0
        @test sum(hrow(:ih_cambialfiremort_si_scpf)) > 0.0

        # ==============================================================
        # 3. TREE-DAMAGE cross-tab — rebuild the history interface with
        #    hlm_use_tree_damage on (registers the *_cdpf / CROWNAREA_*_DAMAGE
        #    handles), construct damaged cohorts + site damage arrays, run dyn2/dyn1.
        # ==============================================================
        CLM.hlm_use_tree_damage[] = CLM.itrue
        histd = CLM._build_fates_history_interface(1)
        @test haskey(histd.ih, :ih_nplant_si_cdpf)
        @test haskey(histd.ih, :ih_mortality_si_cdpf)
        @test haskey(histd.ih, :ih_crownarea_canopy_damage_si)

        nld  = CLM.nlevdamage[]
        nsc  = CLM.nlevsclass[]
        npft = CLM.numpft[]
        @test nld >= 1

        # allocate + zero the site damage arrays (EDInit gates these on the flag)
        site.term_nindivs_canopy_damage = zeros(nld, nsc, npft)
        site.term_nindivs_ustory_damage = zeros(nld, nsc, npft)
        site.imort_rate_damage          = zeros(nld, nsc, npft)
        site.fmort_rate_canopy_damage   = zeros(nld, nsc, npft)
        site.fmort_rate_ustory_damage   = zeros(nld, nsc, npft)
        # put a known mortality count into damage class 1, size 1, pft 1
        site.term_nindivs_canopy_damage[1, 1, 1] = 4.0
        site.fmort_rate_ustory_damage[1, 1, 1]   = 2.0
        site.crownarea_canopy_damage = 3.0
        site.crownarea_ustory_damage = 1.5

        # mark each cohort as damaged (crowndamage class > 1 where >1 exists)
        dmgclass = nld > 1 ? 2 : 1
        for cc in coh
            cc.crowndamage = dmgclass
            cc.dgmort = 0.01
            cc.cmort  = isfinite(cc.cmort)  ? cc.cmort  : 0.0
        end

        iod = site.h_gid
        hrowd(sym) = @view histd.hvars[histd.ih[sym]].r82d[iod, :]
        hbufd(sym) = histd.hvars[histd.ih[sym]].r81d[iod]

        CLM.update_history_dyn1!(histd, 1, fates.nsites, fates.sites, fates.bc_in)
        # CROWNAREA_*_DAMAGE (dyn1, gated on tree-damage)
        @test hbufd(:ih_crownarea_canopy_damage_si) ≈ 3.0 * CLM.days_per_year / CLM.m2_per_ha rtol=1e-10
        @test hbufd(:ih_crownarea_ustory_damage_si) ≈ 1.5 * CLM.days_per_year / CLM.m2_per_ha rtol=1e-10

        CLM.update_history_dyn2!(histd, 1, fates.nsites, fates.sites, fates.bc_in)

        # cohort-level cdpf: nplant lands in the cohort's damage x size x pft bin
        nplant_cdpf = hrowd(:ih_nplant_si_cdpf)
        @test all(isfinite, nplant_cdpf)
        @test sum(nplant_cdpf) > 0.0
        icdpf1 = CLM.get_cdamagesizepft_class_index(coh[1].dbh, coh[1].crowndamage, coh[1].pft)
        @test nplant_cdpf[icdpf1] > 0.0
        # damage mortality (m11) by scpf + cdpf populated
        @test sum(hrowd(:ih_m11_si_scpf)) > 0.0
        @test sum(hrowd(:ih_m11_si_cdpf)) > 0.0
        # site-level cross-tab: the seeded term/fmort counts show up in cdpf mortality
        mort_cdpf = hrowd(:ih_mortality_si_cdpf)
        @test all(isfinite, mort_cdpf)
        @test sum(mort_cdpf) > 0.0
        # the site-array injection (class 1,1,1) contributes to the (icdam=1,scls=1,pft=1) cdpf bin
        site_icdpf = (1 - 1) * nsc + 1 + (1 - 1) * nsc * nld
        @test mort_cdpf[site_icdpf] > 0.0

        CLM.hlm_use_tree_damage[] = CLM.ifalse

        # ==============================================================
        # 4. HYDRAULICS — construct si_hydr + per-cohort co_hydr, enable planthydro,
        #    run update_history_hydraulics!.
        # ==============================================================
        CLM.hlm_use_planthydro[] = CLM.itrue
        bcin = fates.bc_in[1]
        bcin.nlevsoil = nlevsoil
        bcin.h2o_liqvol_sl = fill(0.30, nlevsoil)
        bcin.watsat_sl     = fill(0.45, nlevsoil)
        bcin.dz_sisl       = fill(0.1, nlevsoil)

        # one rhizosphere layer mapping to all soil layers
        sh = CLM.ed_site_hydr_type()
        sh.nlevrhiz = 1
        sh.map_r2s  = reshape([1, nlevsoil], 1, 2)
        sh.dz_rhiz  = [0.1 * nlevsoil]
        sh.l_aroot_layer = [5.0]
        sh.rs1      = [1.0e-4]
        wrf = CLM.wrf_type_cch()
        CLM.set_wrf_param!(wrf, [0.45, -1.0e-4, 6.0])  # th_sat, psi_sat[MPa], beta
        sh.wrf_soil = CLM.WRFType[wrf]
        sh.h2oveg            = 1.5
        sh.h2oveg_hydro_err  = 0.01
        sh.rootuptake_sl     = fill(1.0e-6, nlevsoil)
        nsc  = CLM.nlevsclass[]
        npft = CLM.numpft[]
        sh.sapflow_scpf       = fill(2.0, nsc, npft)
        sh.rootuptake0_scpf   = fill(1.0, nsc, npft)
        sh.rootuptake10_scpf  = fill(1.0, nsc, npft)
        sh.rootuptake50_scpf  = fill(1.0, nsc, npft)
        sh.rootuptake100_scpf = fill(1.0, nsc, npft)
        site.si_hydr = sh

        # per-cohort co_hydr state (1 rhiz layer)
        for cc in coh
            cc.isnew = false
            ch = CLM.ed_cohort_hydr_type()
            ch.v_aroot_layer = [1.0]
            ch.th_aroot  = [0.4]
            ch.psi_aroot = [-0.5]
            ch.ftc_aroot = [0.8]
            ch.th_troot  = 0.42
            ch.psi_troot = -0.4
            ch.ftc_troot = 0.85
            ch.th_ag  = [0.45, 0.43]   # (1)=leaf, (2)=stem
            ch.psi_ag = [-0.3, -0.35]
            ch.ftc_ag = [0.9, 0.88]
            ch.btran  = 0.7
            ch.qtop   = 1.0e-6
            ch.errh2o = 1.0e-9
            ch.iterh1 = 3.0
            ch.iterh2 = 2.0
            cc.co_hydr = ch
        end

        # rebuild a planthydro-aware history interface so the hydro groups register
        histh = CLM._build_fates_history_interface(1)
        @test haskey(histh.ih, :ih_h2oveg_si)
        @test haskey(histh.ih, :ih_btran_scpf)
        ioh = site.h_gid
        hbufh(sym) = histh.hvars[histh.ih[sym]].r81d[ioh]
        hrowh(sym) = @view histh.hvars[histh.ih[sym]].r82d[ioh, :]

        CLM.update_history_hydraulics!(histh, 1, fates.nsites, fates.sites,
                                       fates.bc_in, 1800.0)

        # -- site-level: water storage + total root uptake + root-weighted soil means --
        @test hbufh(:ih_h2oveg_si) ≈ 1.5
        @test hbufh(:ih_h2oveg_hydro_err_si) ≈ 0.01
        @test hbufh(:ih_rootuptake_si) ≈ nlevsoil * 1.0e-6 rtol=1e-12
        @test hbufh(:ih_rootwgt_soilvwc_si) ≈ 0.30 rtol=1e-10        # uniform column
        @test hbufh(:ih_rootwgt_soilvwcsat_si) ≈ 0.45 rtol=1e-10
        @test isfinite(hbufh(:ih_rootwgt_soilmatpot_si))
        @test hbufh(:ih_rootwgt_soilmatpot_si) < 0.0                 # unsaturated -> negative Pa

        # -- size x pft cohort means: btran in (0,1], lands in cohort's scpf bin --
        btran_scpf = hrowh(:ih_btran_scpf)
        @test all(isfinite, btran_scpf)
        @test btran_scpf[coh[1].size_by_pft_class] > 0.0
        # n-weighted mean btran == 0.7 in any occupied bin (all cohorts set to 0.7)
        sb = coh[1].size_by_pft_class
        @test btran_scpf[sb] ≈ 0.7 rtol=1e-8
        # leaf/stem water potentials [Pa] finite + negative; leaf uses ag(1)
        lwp = hrowh(:ih_lwp_scpf); swp = hrowh(:ih_swp_scpf)
        @test all(isfinite, lwp) && all(isfinite, swp)
        @test lwp[sb] ≈ -0.3 * CLM.pa_per_mpa rtol=1e-8
        @test swp[sb] ≈ -0.35 * CLM.pa_per_mpa rtol=1e-8
        # conductivity fractions in (0,1]
        lflc = hrowh(:ih_lflc_scpf)
        @test lflc[sb] ≈ 0.9 rtol=1e-8
        # transpiration flux scpf finite + positive (qtop * n/Σn / dt)
        tran = hrowh(:ih_tran_scpf)
        @test all(isfinite, tran)
        @test tran[sb] > 0.0
        # per-soil-layer + sapflow scpf populated (the FATES soil hist dim may be
        # narrower than nlevsoil; each registered layer carries the per-layer uptake)
        rusl = hrowh(:ih_rootuptake_sl)
        nsl_io = length(rusl)
        @test all(isfinite, rusl)
        @test sum(rusl) ≈ min(nsl_io, nlevsoil) * 1.0e-6 rtol=1e-10
        @test hrowh(:ih_sapflow_scpf)[sb] ≈ 2.0 * CLM.ha_per_m2 rtol=1e-10

        CLM.hlm_use_planthydro[] = CLM.ifalse

    finally
        CLM.prt_global[]                  = old_global
        CLM.num_elements[]                = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        copyto!(CLM.element_pos, old_elpos)
        CLM.hlm_parteh_mode[]             = old_parteh
        CLM.hlm_use_sp[]                  = old_use_sp
        CLM.hlm_use_ed_st3[]              = old_st3
        CLM.hlm_use_planthydro[]          = old_hydro
        CLM.hlm_use_cohort_age_tracking[] = old_agetrk
        CLM.nleafage[]                    = old_nleafage
        CLM.numpft[]                      = old_numpft
        CLM.hlm_current_tod[]             = old_tod
        CLM.nlevsclass[]                  = old_nlevsc
        CLM.nlevcoage[]                   = old_nlevca
        CLM.nlevage[]                     = old_nlevage
        CLM.nlevdamage[]                  = old_nlevdam
        CLM.nlevheight[]                  = old_nlevht
        CLM.SFParams[]                    = old_sfp
        CLM.hlm_is_restart[]              = old_restart
        CLM.hlm_use_nocomp[]              = old_nocomp
        CLM.hlm_use_fixed_biogeog[]       = old_biogeog
        CLM.hlm_use_luh[]                 = old_luh
        CLM.hlm_use_inventory_init[]      = old_inv
        CLM.hlm_use_tree_damage[]         = old_damage
        CLM.hlm_day_of_year[]             = old_doy
        CLM.hlm_days_per_year[]           = old_dpy
        CLM.hlm_freq_day[]                = old_freqday
        CLM.hlm_num_lu_harvest_cats[]     = old_nharv
        CLM.hlm_use_ch4[]                 = old_ch4
        CLM.hlm_use_vertsoilc[]           = old_vert
        CLM.hlm_spitfire_mode[]           = old_spit
        CLM.hlm_hist_level_dynam[]        = old_hdyn
        CLM.hlm_hist_level_hifrq[]        = old_hhfr
        CLM.fates_maxPatchesPerSite[]     = old_paramfp
        CLM.hlm_model_day[]               = old_modelday
        CLM.varpar.nlevsno                = old_nlevsno
        CLM.varpar.nlevsoi                = old_nlevsoi
        CLM.varpar.nlevgrnd               = old_nlevgrnd
        CLM.varpar.nlevmaxurbgrnd         = old_nlevmaxu
    end
end
