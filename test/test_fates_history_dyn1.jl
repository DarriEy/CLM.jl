# test_fates_history_dyn1.jl
# FATES history buffer-fills, Wave 1 (B18-followup H-series):
#   * the history interface is instantiated + Init'd inside clm_fates_init! (W6
#     deferred wiring) and attached as inst.fates.hist,
#   * update_history_dyn1! (daily) is driven through the W5 daily-dynamics hook
#     (fates_daily_dynamics_step!), filling the dynamics-group buffers from the
#     freshly-advanced site/patch/cohort state,
#   * update_history_hifrq1! (per-timestep) is driven through fates_hifrq_history_step!
#     after setting representative per-step cohort carbon fluxes.
#
# Asserts that a sample of dyn1 handles (a biomass pool, an NPP, a mortality flux,
# an LAI) and hifrq1 handles (gpp/npp/nep + organ MR + conductance) hold FINITE,
# physical values (>= 0 where expected; not the ignore value), and that a variable
# in a deferred (not-yet-filled) fill-group stays at its flush value.
#
# clm_fates_init! mutates FATES module-global control Refs + parameter tables;
# save/restore them so the surrounding default suite stays unperturbed.

using Test
using CLM

@testset "FATES history W1 (dyn1 core + hifrq1)" begin
    # ---- save FATES globals (same set as the W5 daily-dynamics test) ----
    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
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
        pveg = col.patchi[cveg] + 1

        # ---- cold-started carbon-only FATES site (builds inst.fates.hist) ----
        fates = CLM.clm_fates_init!(inst; nsites = 1,
                                    nlevsoil = nlevsoil, nlevdecomp = nlevdecomp)

        # ==============================================================
        # 1. The history interface is wired + Init'd in the run.
        # ==============================================================
        @test inst.fates !== nothing
        @test inst.fates.hist !== nothing
        hist = inst.fates.hist
        @test hist.num_history_vars_ > 0
        @test length(hist.hvars) == hist.num_history_vars_
        # site h_gid mapped to its history-io slot
        @test inst.fates.sites[1].h_gid == 1
        # both update groups active (set in clm_fates_init!)
        @test CLM.hlm_hist_level_dynam[] > 0
        @test CLM.hlm_hist_level_hifrq[] > 0

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

        # handy buffer accessor by ih_* handle
        hbuf(sym) = hist.hvars[hist.ih[sym]].r81d[io]

        # ==============================================================
        # 2. Drive the DAILY step -> update_history_dyn1! fills the buffers.
        # ==============================================================
        CLM.fates_daily_dynamics_step!(inst; nlevsoil=nlevsoil,
                                       nlevdecomp=nlevdecomp, use_fates_bgc=false)

        # -- counts --
        @test hbuf(:ih_npatches_si) >= 1.0
        @test hbuf(:ih_ncohorts_si) >= 1.0

        # -- biomass pools (finite, non-negative) --
        for sym in (:ih_totvegc_si, :ih_leafc_si, :ih_fnrtc_si, :ih_sapwc_si,
                    :ih_storec_si, :ih_agb_si, :ih_bdead_si, :ih_balive_si,
                    :ih_canopy_biomass_si)
            v = hbuf(sym)
            @test isfinite(v)
            @test v >= 0.0
        end
        # the seeded cohort carries real carbon — at least one pool is positive
        @test hbuf(:ih_totvegc_si) > 0.0

        # -- LAI / area aggregates --
        for sym in (:ih_lai_si, :ih_elai_si, :ih_area_plant_si, :ih_area_trees_si,
                    :ih_trimming_si)
            v = hbuf(sym)
            @test isfinite(v)
            @test v >= 0.0
        end

        # -- mortality / flux diagnostics (finite; zero-or-positive) --
        for sym in (:ih_canopy_mortality_carbonflux_si,
                    :ih_understory_mortality_carbonflux_si,
                    :ih_demotion_carbonflux_si, :ih_promotion_carbonflux_si,
                    :ih_litter_in_si, :ih_litter_out_si, :ih_seed_bank_si,
                    :ih_cbal_err_fates_si, :ih_fates_fraction_si)
            @test isfinite(hbuf(sym))
        end
        # fates fraction is set to exactly 1 on a fates column
        @test hbuf(:ih_fates_fraction_si) == 1.0

        # -- site scalars copied through --
        @test isfinite(hbuf(:ih_gdd_si))
        @test isfinite(hbuf(:ih_canopy_spread_si))

        # -- a DEFERRED (dyn2 / per-scpf) variable stays at its flush value --
        # update_history_dyn2! is still stubbed, so its 2D buffers are untouched.
        if hist.ih[:ih_totvegc_scpf] > 0
            hv2 = hist.hvars[hist.ih[:ih_totvegc_scpf]]
            @test all(hv2.r82d[io, :] .== hv2.flushval)
        end

        # ==============================================================
        # 3. Drive the HIFRQ (per-timestep) fill. Set representative finite
        #    per-step cohort carbon fluxes first (these are NaN until the
        #    photosynthesis driver runs; the daily step has cleared isnew).
        # ==============================================================
        ncoh_set = 0
        pp = site.oldest_patch
        while pp !== nothing
            # The carbon-only cold start leaves total_canopy_area unset (NaN) on a
            # patch whose canopy hasn't been summarized; the in-solve photosynthesis
            # path would set it. Give it a representative finite value so the
            # vegetated-area normalizer (and conductance/tveg weighting) is exercised.
            if !isfinite(pp.total_canopy_area)
                pp.total_canopy_area = pp.area
            end
            # patch conductances + per-band radiation error so the conductance /
            # rad-error diagnostics get a finite, non-ignore value.
            pp.c_stomata = 0.2
            pp.c_lblayer = 2.0
            if !isempty(pp.rad_error)
                pp.rad_error[CLM.ivis] = 1.0e-3
                pp.rad_error[CLM.inir] = 2.0e-3
            end
            cc = pp.shortest
            while cc !== nothing
                cc.isnew            = false
                cc.gpp_tstep        = 1.0e-6
                cc.npp_tstep        = 6.0e-7
                cc.resp_tstep       = 4.0e-7
                cc.resp_g_tstep     = 1.0e-7
                cc.resp_m           = 3.0e-7
                cc.resp_m_unreduced = 3.0e-7
                cc.rdark            = 1.0e-8
                cc.froot_mr         = 1.0e-8
                cc.livecroot_mr     = 1.0e-8
                cc.livestem_mr      = 1.0e-8
                ncoh_set += 1
                cc = cc.taller
            end
            pp = pp.younger
        end
        @test ncoh_set >= 1

        # per-step veg temperature into bc_in (read by the TVEG diagnostic)
        bcin = inst.fates.bc_in[1]
        if !isempty(bcin.t_veg_pa)
            fill!(bcin.t_veg_pa, 290.0)
        end

        dtime = 1800.0
        CLM.fates_hifrq_history_step!(inst; dt_tstep=dtime)

        # hifrq core: finite, physically-signed
        @test isfinite(hbuf(:ih_gpp_si)) && hbuf(:ih_gpp_si) > 0.0
        @test isfinite(hbuf(:ih_npp_si)) && hbuf(:ih_npp_si) > 0.0
        @test isfinite(hbuf(:ih_nep_si))                       # = NPP - HR (HR=0 here)
        @test isfinite(hbuf(:ih_hr_si))  && hbuf(:ih_hr_si) >= 0.0
        @test isfinite(hbuf(:ih_aresp_si)) && hbuf(:ih_aresp_si) > 0.0
        @test isfinite(hbuf(:ih_growth_resp_si)) && hbuf(:ih_growth_resp_si) >= 0.0
        @test isfinite(hbuf(:ih_maint_resp_si))  && hbuf(:ih_maint_resp_si)  >= 0.0
        # organ maintenance respiration
        for sym in (:ih_leaf_mr_si, :ih_froot_mr_si, :ih_livecroot_mr_si, :ih_livestem_mr_si)
            @test isfinite(hbuf(sym)) && hbuf(sym) >= 0.0
        end
        # canopy/understory partition
        @test isfinite(hbuf(:ih_gpp_canopy_si)) || isfinite(hbuf(:ih_gpp_understory_si))
        # conductances are vegetated-area weighted finite values (not the ignore val)
        @test isfinite(hbuf(:ih_c_stomata_si))
        @test hbuf(:ih_c_stomata_si) != CLM.hlm_hio_ignore_val[]
        @test isfinite(hbuf(:ih_c_lblayer_si))
        @test hbuf(:ih_c_lblayer_si) != CLM.hlm_hio_ignore_val[]
        # NEP == NPP (no heterotrophic respiration in carbon-only) to within fp
        @test hbuf(:ih_nep_si) ≈ hbuf(:ih_npp_si) atol = 1e-20

        # radiation-error diagnostics got the solver-weighted value (not ignore)
        @test isfinite(hbuf(:ih_vis_rad_err_si))
        @test hbuf(:ih_vis_rad_err_si) != CLM.hlm_hio_ignore_val[]

    finally
        CLM.prt_global[]                  = old_global
        CLM.num_elements[]                = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
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
        CLM.varpar.nlevsno                = old_nlevsno
        CLM.varpar.nlevsoi                = old_nlevsoi
        CLM.varpar.nlevgrnd               = old_nlevgrnd
        CLM.varpar.nlevmaxurbgrnd         = old_nlevmaxu
    end
end
