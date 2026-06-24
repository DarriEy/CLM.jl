# test_fates_history_dyn2.jl
# FATES history buffer-fills, Wave 2 (B18-followup H-series):
#   * update_history_dyn2! (daily, level-2) — the size/PFT/age-class DISAGGREGATED
#     daily diagnostics: biomass/NPP/mortality/number-density/basal-area/crown-area
#     binned into the *_si_pft / *_si_scls / *_si_scpf / *_si_scag / *_si_age handles,
#     plus the mortality m1..m10 series and CWD-by-cwdsc.
#   * update_history_hifrq2! (per-timestep, level-2) — the per-canopy-layer × PFT ×
#     leaf-layer PAR/LAI profiles + the respiration disaggregation by scpf/scls.
#
# Both level-2 updates only run when the history dim-level is > 1, so we cold-start
# with hist_dimlevel=2 (the registry then allocates the dynam1/hifrq1-gated 2D
# vars). Asserts that a sample of the newly-filled DISAGGREGATED handles populate
# FINITE / physical values AND that the cohort's mass lands in the correct
# size/pft class bin (the class-binning is verified against the cohort's own
# computed size_class / size_by_pft_class). Confirms a W2b-deferred group
# (fire/damage/landuse-dimensioned) stays at its flush value.
#
# clm_fates_init! mutates FATES module-global control Refs + parameter tables;
# save/restore them so the surrounding default suite stays unperturbed.

using Test
using CLM

@testset "FATES history W2 (dyn2 + hifrq2)" begin
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
        @test inst.fates.hist !== nothing
        hist = inst.fates.hist
        @test CLM.hlm_hist_level_dynam[] == 2
        @test CLM.hlm_hist_level_hifrq[] == 2

        # level-2 (dynam1/hifrq1-gated) 2D vars must now be registered
        for sym in (:ih_biomass_si_scls, :ih_nplant_si_scpf, :ih_biomass_si_pft,
                    :ih_gpp_si_age, :ih_m1_si_scpf, :ih_ar_si_scpf, :ih_parsun_z_si_cnlf)
            @test haskey(hist.ih, sym) && hist.ih[sym] > 0
        end

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

        # 2D buffer accessor by ih_* handle (returns the io_si row)
        hrow(sym) = @view hist.hvars[hist.ih[sym]].r82d[io, :]

        # ==============================================================
        # 1. Drive the DAILY step -> update_history_dyn2! fills the level-2 buffers.
        # ==============================================================
        CLM.fates_daily_dynamics_step!(inst; nlevsoil=nlevsoil,
                                       nlevdecomp=nlevdecomp, use_fates_bgc=false)

        # gather the live cohorts' class indices so we can verify the binning.
        coh_scls = Int[]; coh_scpf = Int[]; coh_pft = Int[]
        total_n  = 0.0
        pp = site.oldest_patch
        while pp !== nothing
            cc = pp.shortest
            while cc !== nothing
                if cc.prt !== nothing
                    push!(coh_scls, cc.size_class)
                    push!(coh_scpf, cc.size_by_pft_class)
                    push!(coh_pft,  cc.pft)
                    total_n += cc.n
                end
                cc = cc.taller
            end
            pp = pp.younger
        end
        @test !isempty(coh_scls)

        # -- *_si_scls biomass: finite, non-negative, with mass in the cohort's bin --
        bio_scls = hrow(:ih_biomass_si_scls)
        @test all(isfinite, bio_scls)
        @test all(>=(0.0), bio_scls)
        @test bio_scls[coh_scls[1]] > 0.0          # CORRECT size-class binning
        # nothing should have landed in a bin no cohort occupies
        for k in eachindex(bio_scls)
            if !(k in coh_scls)
                @test bio_scls[k] == 0.0
            end
        end

        # -- *_si_scpf (size x pft) number density: lands in the joint bin --
        np_scpf = hrow(:ih_nplant_si_scpf)
        @test all(isfinite, np_scpf)
        @test np_scpf[coh_scpf[1]] > 0.0           # CORRECT size x pft binning
        # total number density (per ha->per m2) conserved across the scpf bins
        @test sum(np_scpf) ≈ total_n / CLM.m2_per_ha rtol = 1e-10

        # -- *_si_pft aggregate biomass: finite, in the cohort's pft column --
        bio_pft = hrow(:ih_biomass_si_pft)
        @test all(isfinite, bio_pft)
        @test bio_pft[coh_pft[1]] > 0.0
        # nindivs by pft also conserved
        nind_pft = hrow(:ih_nindivs_si_pft)
        @test sum(nind_pft) ≈ total_n * CLM.AREA_INV rtol = 1e-10

        # -- *_si_age: at least one age bin carries area + biomass --
        area_age = hrow(:ih_area_si_age)
        @test all(isfinite, area_age)
        @test sum(area_age) > 0.0
        bio_age = hrow(:ih_biomass_si_age)
        @test all(isfinite, bio_age)
        @test sum(bio_age) > 0.0

        # -- a mortality m*_si_scpf handle is finite + non-negative --
        m1_scpf = hrow(:ih_m1_si_scpf)
        @test all(isfinite, m1_scpf)
        @test all(>=(0.0), m1_scpf)

        # -- carbon pool scpf state (filled in the element loop) --
        totvegc_scpf = hrow(:ih_totvegc_scpf)
        @test all(isfinite, totvegc_scpf)
        @test totvegc_scpf[coh_scpf[1]] > 0.0

        # -- a W2b-DEFERRED group stays UNFILLED --
        # update_history_dyn2! zeros every group_dyna_complx slot at site-zero, but
        # the fire-dimensioned fills are deferred (cpatch.fuel is nothing / FI unset
        # in carbon-only), so the fire-intensity-by-age handle sits at exactly 0.0
        # (zeroed, never accumulated) — distinct from the filled handles above.
        fire_int = hrow(:ih_fire_intensity_si_age)
        @test all(==(0.0), fire_int)

        # ==============================================================
        # 2. Drive the HIFRQ level-2 fill. Set representative finite per-step
        #    cohort fluxes + patch radiation profiles first.
        # ==============================================================
        # The carbon-only cold start does NOT summarize the canopy (radiation runs
        # only in the full driver): canopy_area_profile is unallocated, nleaf is
        # unset, cohort nv == 0. Build a representative 1-canopy-layer profile here
        # (set nleaf -> ReAllocateDynamics! -> fill the per-bin radiation/LAI) so the
        # cnlf / cnlfpft / can profile fills are exercised end-to-end.
        nleaf_test = 2
        pp = site.oldest_patch
        while pp !== nothing
            if !isfinite(pp.total_canopy_area)
                pp.total_canopy_area = pp.area
            end
            pp.c_stomata = 0.2
            pp.c_lblayer = 2.0
            pp.ncl_p = 1
            for ip in 1:CLM.numpft[]
                pp.nleaf[1, ip] = nleaf_test
            end
            CLM.ReAllocateDynamics!(pp)
            for ip in 1:CLM.numpft[], il in 1:nleaf_test
                pp.canopy_area_profile[1, ip, il] = 0.5
                pp.ed_parsun_z[1, ip, il] = 100.0
                pp.ed_parsha_z[1, ip, il] = 20.0
                pp.elai_profile[1, ip, il] = 0.4
                pp.f_sun[1, ip, il] = 0.6
                pp.parprof_pft_dir_z[1, ip, il] = 80.0
                pp.parprof_pft_dif_z[1, ip, il] = 15.0
            end
            cc = pp.shortest
            while cc !== nothing
                cc.isnew        = false
                cc.nv           = nleaf_test
                cc.gpp_tstep    = 1.0e-6
                cc.npp_tstep    = 6.0e-7
                cc.resp_tstep   = 4.0e-7
                cc.resp_g_tstep = 1.0e-7
                cc.resp_m       = 3.0e-7
                cc.rdark        = 1.0e-8
                cc.froot_mr     = 1.0e-8
                cc.livecroot_mr = 1.0e-8
                cc.livestem_mr  = 1.0e-8
                if !isempty(cc.ts_net_uptake)
                    fill!(cc.ts_net_uptake, 1.0e-9)
                end
                cc = cc.taller
            end
            pp = pp.younger
        end

        dtime = 1800.0
        CLM.fates_hifrq_history_step!(inst; dt_tstep=dtime)

        # -- respiration by size x pft (ih_ar_si_scpf): finite, lands in cohort bin --
        ar_scpf = hrow(:ih_ar_si_scpf)
        @test all(isfinite, ar_scpf)
        @test all(>=(0.0), ar_scpf)
        @test ar_scpf[coh_scpf[1]] > 0.0           # CORRECT respiration binning

        # -- canopy-layer x leaf radiation profile (ih_parsun_z_si_cnlf) --
        parsun_cnlf = hrow(:ih_parsun_z_si_cnlf)
        # at least one bin is a finite, positive, non-ignore PAR value
        finite_pos = filter(v -> isfinite(v) && v != CLM.hlm_hio_ignore_val[] && v > 0.0,
                            collect(parsun_cnlf))
        @test !isempty(finite_pos)

        # -- LAI sunlit profile (ih_laisun_z_si_cnlf) populated --
        laisun_cnlf = hrow(:ih_laisun_z_si_cnlf)
        finite_lai = filter(v -> isfinite(v) && v != CLM.hlm_hio_ignore_val[] && v > 0.0,
                            collect(laisun_cnlf))
        @test !isempty(finite_lai)

        # -- GPP by age (ih_gpp_si_age) finite + a populated age bin --
        gpp_age = hrow(:ih_gpp_si_age)
        @test all(isfinite, gpp_age)
        @test sum(gpp_age) > 0.0

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
        CLM.hlm_model_day[]               = old_modelday
        CLM.varpar.nlevsno                = old_nlevsno
        CLM.varpar.nlevsoi                = old_nlevsoi
        CLM.varpar.nlevgrnd               = old_nlevgrnd
        CLM.varpar.nlevmaxurbgrnd         = old_nlevmaxu
    end
end
