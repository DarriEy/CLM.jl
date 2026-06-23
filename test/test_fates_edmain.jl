# test_fates_edmain.jl
# Tests for FATES Batch 16 (Tier F): EDMainMod — the DAILY ecosystem-dynamics
# DRIVER (ed_ecosystem_dynamics / ed_integrate_state_variables / ed_update_site /
# TotalBalanceCheck / bypass_dynamics).
#
# NOTE: three sibling modules (FatesNormanRadMod / FatesPlantRespPhotosynthMod /
# FatesInventoryInitMod) are ported in PARALLEL in this batch and are NOT present
# in this worktree, so the full ed_ecosystem_dynamics / ed_update_site paths that
# reach photosynthesis / radiation / inventory cannot be RUN here. These tests
# therefore exercise only the orchestration logic that is self-contained:
#   * bypass_dynamics — sets cohort fluxes/rates/flags to trivial values.
#   * TotalBalanceCheck — passes on a balanced site, FLAGS (aborts) an imbalance,
#     updates old_stock / err_fates only on the final (-1) call, and is entirely
#     SKIPPED in SP mode.
#   * SiteMassStock round-trip used by the balance check.
#
# Setup mirrors test_fates_edpatchdynamics.jl (single woody carbon PFT).

using Test
using CLM

# --- single woody evergreen carbon PFT (same allometry as the cohort suite) ---
function _setup_edmain_pft!()
    npft = 1
    p = CLM.prt_params
    CLM.allocate_prt_params!(p, npft, CLM.num_organ_types, 1)

    p.c2b          .= 2.0
    p.wood_density .= 0.6
    p.slatop       .= 0.012
    p.slamax       .= 0.020
    p.allom_agb_frac .= 0.6
    p.allom_dbh_maxheight .= 90.0
    p.allom_la_per_sa_int .= 0.8e3
    p.allom_la_per_sa_slp .= 0.0
    p.allom_sai_scaler    .= 0.1
    p.allom_l2fr          .= 1.0
    p.cushion             .= 1.0
    p.allom_blca_expnt_diff      .= 0.0
    p.allom_d2ca_coefficient_min .= 0.3
    p.allom_d2ca_coefficient_max .= 0.6
    p.allom_h2cd1 .= 0.5
    p.allom_h2cd2 .= 1.0

    p.allom_dmode .= 1
    p.woody       .= 1
    p.allom_cmode .= 1
    p.allom_smode .= 1
    p.allom_fmode .= 1
    p.allom_stmode .= 1
    p.allom_hmode .= 1
    p.allom_amode .= 1
    p.allom_lmode .= 1

    p.allom_d2h1 .= 0.64
    p.allom_d2h2 .= 0.37
    p.allom_d2h3 .= -0.034

    p.allom_agb1 .= 0.06896; p.allom_agb2 .= 0.572; p.allom_agb3 .= 1.94; p.allom_agb4 .= 0.931
    p.allom_d2bl1 .= 0.07; p.allom_d2bl2 .= 1.3; p.allom_d2bl3 .= 0.55

    p.fnrt_prof_mode .= 1
    p.fnrt_prof_a    .= 0.976
    p.fnrt_prof_b    .= 0.0

    p.stress_decid .= 0
    p.season_decid .= 0
    p.evergreen    .= CLM.itrue
    p.leaf_stor_priority .= 0.8
    p.leaf_long          .= 1.5
    p.leaf_long_ustory   .= 1.5

    p.seed_alloc         .= 0.1
    p.seed_alloc_mature  .= 0.0
    p.dbh_repro_threshold .= 1000.0
    p.repro_alloc_a      .= 0.0
    p.repro_alloc_b      .= 0.0

    pd = CLM.param_derived_type()
    pd.branch_frac = fill(0.25, npft)
    pd.jmax25top = fill(85.0, npft, 1)
    pd.tpu25top  = fill(8.0, npft, 1)
    pd.kp25top   = fill(0.6, npft, 1)
    CLM.ParamDerived[] = pd

    ev = CLM.EDPftvarcon_type()
    ev.vcmax25top = fill(50.0, npft, 1)
    ev.damage_recovery_scalar = fill(0.0, npft)
    ev.lf_flab = fill(0.25, npft); ev.lf_fcel = fill(0.50, npft); ev.lf_flig = fill(0.25, npft)
    ev.fr_flab = fill(0.25, npft); ev.fr_fcel = fill(0.50, npft); ev.fr_flig = fill(0.25, npft)
    ev.allom_frbstor_repro = fill(0.0, npft)
    ev.landusechange_frac_burned   = fill(0.0, npft)
    ev.landusechange_frac_exported = fill(0.0, npft)
    CLM.EDPftvarcon_inst[] = ev

    edp = CLM.ed_params_type()
    edp.regeneration_model = CLM.default_regeneration
    edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
    edp.ED_val_history_sizeclass_bin_edges  = [0.0, 5.0, 10.0, 20.0, 50.0]
    edp.ED_val_history_coageclass_bin_edges = [0.0, 3.0, 6.0, 9.0, 12.0]
    edp.ED_val_history_ageclass_bin_edges   = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    edp.ED_val_cohort_size_fusion_tol = 0.06
    edp.ED_val_cohort_age_fusion_tol  = 0.08
    edp.max_cohort_per_patch          = 100
    edp.ED_val_patch_fusion_tol       = 0.05
    edp.fates_mortality_disturbance_fraction = 1.0
    edp.ED_val_understorey_death      = 0.55983
    edp.logging_coll_under_frac       = 0.55983
    edp.maxpatches_by_landuse         = fill(10, CLM.n_landuse_cats)
    nlv = CLM.nlevleaf
    edp.dinc_vai   = fill(1.0, nlv)
    edp.dlower_vai = collect(1.0:1.0:Float64(nlv))
    CLM.EDParams[] = edp

    CLM.nlevsclass[] = length(edp.ED_val_history_sizeclass_bin_edges)
    CLM.nlevcoage[]  = length(edp.ED_val_history_coageclass_bin_edges)

    sfp = CLM.sf_params_type()
    sfp.SF_val_CWD_frac = [0.045, 0.075, 0.21, 0.67]
    CLM.SFParams[] = sfp

    return npft
end

function _seed_prt(ipft::Int, dbh0::Float64)
    prt = CLM.callom_prt_vartypes()
    CLM.InitPRTVartype!(prt)
    canopy_trim = 1.0; crowndamage = 1; ef = 1.0
    l2fr = CLM.prt_params.allom_l2fr[ipft]
    tgt_leaf, _    = CLM.bleaf(dbh0, ipft, crowndamage, canopy_trim, ef)
    tgt_fnrt, _    = CLM.bfineroot(dbh0, ipft, canopy_trim, l2fr, ef)
    _, tgt_sapw, _ = CLM.bsap_allom(dbh0, ipft, crowndamage, canopy_trim, ef)
    tgt_store, _   = CLM.bstore_allom(dbh0, ipft, crowndamage, canopy_trim)
    tgt_agw, _     = CLM.bagw_allom(dbh0, ipft, crowndamage, ef)
    tgt_bgw, _     = CLM.bbgw_allom(dbh0, ipft, ef)
    tgt_struct, _  = CLM.bdead_allom(tgt_agw, tgt_bgw, tgt_sapw, ipft)
    CLM.SetState!(prt, CLM.leaf_organ,   CLM.carbon12_element, tgt_leaf)
    CLM.SetState!(prt, CLM.fnrt_organ,   CLM.carbon12_element, tgt_fnrt)
    CLM.SetState!(prt, CLM.sapw_organ,   CLM.carbon12_element, tgt_sapw)
    CLM.SetState!(prt, CLM.store_organ,  CLM.carbon12_element, tgt_store)
    CLM.SetState!(prt, CLM.struct_organ, CLM.carbon12_element, tgt_struct)
    CLM.SetState!(prt, CLM.repro_organ,  CLM.carbon12_element, 0.0)
    return prt
end

function _build_site(npft::Int, nlevsoil::Int)
    site = CLM.ed_site_type()
    nmt = CLM.n_term_mort_types
    nsc = CLM.nlevsclass[]

    site.term_nindivs_canopy    = zeros(nmt, nsc, npft)
    site.term_nindivs_ustory    = zeros(nmt, nsc, npft)
    site.term_carbonflux_canopy = zeros(nmt, npft)
    site.term_carbonflux_ustory = zeros(nmt, npft)
    site.term_abg_flux          = zeros(nsc, npft)
    site.growthflux_fusion      = zeros(nsc, npft)
    site.spread                 = 1.0

    site.nlevsoil     = nlevsoil
    site.zi_soil      = collect(0.0:0.1:0.1 * nlevsoil)
    site.rootfrac_scr = zeros(nlevsoil)

    site.flux_diags = CLM.site_fluxdiags_type()
    site.flux_diags.elem = [CLM.elem_diag_type() for _ in 1:CLM.num_elements[]]
    for el in 1:CLM.num_elements[]
        site.flux_diags.elem[el].surf_fine_litter_input = zeros(npft)
        site.flux_diags.elem[el].root_litter_input      = zeros(npft)
    end

    site.mass_balance  = [CLM.site_massbal_type() for _ in 1:CLM.num_elements[]]
    site.iflux_balance = [CLM.site_ifluxbal_type() for _ in 1:CLM.num_elements[]]
    site.area_by_age   = zeros(length(CLM.ed_params().ED_val_history_ageclass_bin_edges))

    site.lat = 50.0
    site.lon = 250.0
    return site
end

function _make_patch(area::Float64, npft::Int, nlevsoil::Int; land_use_label=CLM.primaryland)
    patch = CLM.fates_patch_type()
    CLM.Create!(patch, 0.0, area, land_use_label, CLM.fates_unset_int, CLM.num_swb,
                npft, nlevsoil, 0, CLM.ed_params().regeneration_model)
    for el in 1:CLM.num_elements[]
        CLM.InitConditions!(patch.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
    patch.canopy_layer_tlai = zeros(Float64, CLM.nclmax)
    return patch
end

function _add_cohort!(site, patch, ipft, dbh0, nn; clayer=1)
    prt = _seed_prt(ipft, dbh0)
    h, _ = CLM.h_allom(dbh0, ipft)
    CLM.create_cohort(site, patch, ipft, nn, h, 0.0, dbh0, prt, 1.0, 1.0, 1.0,
                      CLM.leaves_on, 0, 1.0, 0.0, clayer, 1, site.spread)
    c = patch.tallest
    while c !== nothing; c.isnew = false; c = c.shorter; end
    return patch.tallest
end

@testset "FATES Batch 16: EDMainMod" begin

    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_sp   = CLM.hlm_use_sp[]
    old_st3      = CLM.hlm_use_ed_st3[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_agetrk   = CLM.hlm_use_cohort_age_tracking[]
    old_nocomp   = CLM.hlm_use_nocomp[]
    old_luh      = CLM.hlm_use_luh[]
    old_fixbio   = CLM.hlm_use_fixed_biogeog[]
    old_freqday  = CLM.hlm_freq_day[]
    old_dpy      = CLM.hlm_days_per_year[]
    old_tod      = CLM.hlm_current_tod[]
    old_numpft   = CLM.numpft[]
    old_nharv    = CLM.hlm_num_lu_harvest_cats[]
    old_nleafage = CLM.nleafage[]
    old_nlevsc   = CLM.nlevsclass[]
    old_nlevca   = CLM.nlevcoage[]
    old_sfp      = CLM.SFParams[]

    try
        npft = _setup_edmain_pft!()
        ipft = 1
        CLM.InitPRTGlobalAllometricCarbon!()

        CLM.num_elements[] = 1
        empty!(CLM.element_list)
        push!(CLM.element_list, CLM.carbon12_element)

        CLM.hlm_parteh_mode[]             = CLM.prt_carbon_allom_hyp
        CLM.hlm_use_sp[]                  = CLM.ifalse
        CLM.hlm_use_ed_st3[]              = CLM.ifalse
        CLM.hlm_use_planthydro[]          = CLM.ifalse
        CLM.hlm_use_cohort_age_tracking[] = CLM.ifalse
        CLM.hlm_use_nocomp[]              = CLM.ifalse
        CLM.hlm_use_luh[]                 = CLM.ifalse
        CLM.hlm_use_fixed_biogeog[]       = CLM.ifalse
        CLM.hlm_freq_day[]                = 1.0 / 365.0
        CLM.hlm_days_per_year[]           = 365
        CLM.hlm_current_tod[]             = 0
        CLM.numpft[]                      = npft
        CLM.hlm_num_lu_harvest_cats[]     = 1
        CLM.nleafage[]                    = 1

        nlevsoil = 3

        # ==============================================================
        # bypass_dynamics: sets trivial cohort fluxes / rates / flags.
        # ==============================================================
        bsite = _build_site(npft, nlevsoil)
        bp = _make_patch(CLM.area, npft, nlevsoil)
        bp.older = nothing; bp.younger = nothing
        bsite.oldest_patch = bp; bsite.youngest_patch = bp
        bc = _add_cohort!(bsite, bp, ipft, 12.0, 0.03; clayer=1)

        # Pre-load non-trivial accumulators / rates / flags so we can see them reset.
        bc.isnew    = true
        bc.npp_acc  = 2.0; bc.gpp_acc = 3.0; bc.resp_acc = 1.0
        bc.bmort = 9.0; bc.hmort = 9.0; bc.cmort = 9.0; bc.frmort = 9.0
        bc.smort = 9.0; bc.asmort = 9.0; bc.dgmort = 9.0
        bc.dndt = 9.0; bc.dhdt = 9.0; bc.ddbhdt = 9.0

        CLM.bypass_dynamics(bsite)

        @test bc.isnew == false
        # _hold = _acc * days_per_year, captured BEFORE the _acc are zeroed.
        @test bc.npp_acc_hold  ≈ 2.0 * 365
        @test bc.gpp_acc_hold  ≈ 3.0 * 365
        @test bc.resp_acc_hold ≈ 1.0 * 365
        @test bc.npp_acc  == 0.0
        @test bc.gpp_acc  == 0.0
        @test bc.resp_acc == 0.0
        @test bc.bmort == 0.0 && bc.hmort == 0.0 && bc.cmort == 0.0
        @test bc.frmort == 0.0 && bc.smort == 0.0 && bc.asmort == 0.0 && bc.dgmort == 0.0
        @test bc.dndt == 0.0 && bc.dhdt == 0.0 && bc.ddbhdt == 0.0

        # ==============================================================
        # TotalBalanceCheck: passes on a balanced site (no net flux,
        # no stock change), and updates old_stock/err_fates on the
        # final (-1) call only.
        # ==============================================================
        tsite = _build_site(npft, nlevsoil)
        tp = _make_patch(CLM.area, npft, nlevsoil)
        tp.older = nothing; tp.younger = nothing
        tsite.oldest_patch = tp; tsite.youngest_patch = tp
        _add_cohort!(tsite, tp, ipft, 15.0, 0.02; clayer=1)

        # All mass-balance fluxes start at zero; prime old_stock to the current
        # stock so change_in_stock == 0 -> no error regardless of fluxes.
        total0, _, _, _ = CLM.SiteMassStock(tsite, 1)
        @test total0 > 0.0
        tsite.mass_balance[1].old_stock = total0

        # A non-final call should NOT touch old_stock / err_fates and must pass.
        tsite.mass_balance[1].err_fates = -123.0
        CLM.TotalBalanceCheck(tsite, 2)
        @test tsite.mass_balance[1].old_stock == total0
        @test tsite.mass_balance[1].err_fates == -123.0   # untouched (non-final)

        # The final (-1) call records old_stock = current total + err_fates.
        CLM.TotalBalanceCheck(tsite, -1)
        @test tsite.mass_balance[1].old_stock ≈ total0
        # net_flux (all zero) - change_in_stock (zero) == 0
        @test tsite.mass_balance[1].err_fates ≈ 0.0 atol = 1e-12

        # ==============================================================
        # TotalBalanceCheck: a positive uncompensated stock change with no
        # matching flux must be FLAGGED (aborts via fates_endrun -> error).
        # ==============================================================
        usite = _build_site(npft, nlevsoil)
        up = _make_patch(CLM.area, npft, nlevsoil)
        up.older = nothing; up.younger = nothing
        usite.oldest_patch = up; usite.youngest_patch = up
        _add_cohort!(usite, up, ipft, 15.0, 0.02; clayer=1)

        # old_stock = 0 while the live biomass is large -> change_in_stock huge
        # and positive, but flux_in - flux_out == 0 -> large error_frac.
        usite.mass_balance[1].old_stock = 0.0
        @test_throws ErrorException CLM.TotalBalanceCheck(usite, 5)

        # ==============================================================
        # SP mode short-circuits TotalBalanceCheck entirely: even a wildly
        # imbalanced site does NOT throw.
        # ==============================================================
        CLM.hlm_use_sp[] = CLM.itrue
        usite.mass_balance[1].old_stock = 0.0   # still imbalanced
        CLM.TotalBalanceCheck(usite, 5)         # no throw in SP mode
        @test true
        CLM.hlm_use_sp[] = CLM.ifalse

    finally
        CLM.prt_global[]                  = old_global
        CLM.num_elements[]                = old_numel
        empty!(CLM.element_list)
        append!(CLM.element_list, old_ellist)
        CLM.hlm_parteh_mode[]             = old_parteh
        CLM.hlm_use_sp[]                  = old_use_sp
        CLM.hlm_use_ed_st3[]              = old_st3
        CLM.hlm_use_planthydro[]          = old_hydro
        CLM.hlm_use_cohort_age_tracking[] = old_agetrk
        CLM.hlm_use_nocomp[]              = old_nocomp
        CLM.hlm_use_luh[]                 = old_luh
        CLM.hlm_use_fixed_biogeog[]       = old_fixbio
        CLM.hlm_freq_day[]                = old_freqday
        CLM.hlm_days_per_year[]           = old_dpy
        CLM.hlm_current_tod[]             = old_tod
        CLM.numpft[]                      = old_numpft
        CLM.hlm_num_lu_harvest_cats[]     = old_nharv
        CLM.nleafage[]                    = old_nleafage
        CLM.nlevsclass[]                  = old_nlevsc
        CLM.nlevcoage[]                   = old_nlevca
        CLM.SFParams[]                    = old_sfp
    end
end
