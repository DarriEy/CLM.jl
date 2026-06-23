# test_fates_edpatchdynamics.jl
# Tests for FATES Batch 14 (Tier F): EDPatchDynamicsMod — the patch DISTURBANCE
# ENGINE that creates / splits / fuses / terminates patches on the site's
# age-ordered `older`/`younger` patch linked list, and moves cohorts + litter
# between patches during disturbance.
#
# Strategy (mirrors test_fates_edcohortdynamics.jl's setup):
#   * Build a synthetic single-PFT (woody, carbon-only) prt_params table, install
#     param_derived, the FATES PFT table, and ed_params (regeneration + bins +
#     fusion tolerances + max_cohort_per_patch + maxpatches_by_landuse +
#     understorey-death + mortality-disturbance fraction + patch-fusion tol).
#   * Register the carbon-only PARTEH hypothesis; flags = carbon, no SP, no plant
#     hydro, no cohort-age-tracking, no nocomp, no luh, no fixed-biogeog, nleafage=1.
#   * Build a minimal ed_site_type (termination/imort/fmort/growth tallies,
#     flux_diags, soil layering + root scratch, mass_balance, disturbance_rates).
#   * Build patches with Create! (litter, fuel, running means) so the disturbance
#     paths have real pools.
#   * Exercise (invariant-based, hand-computed where possible):
#       - set_patchno / countPatches / check_patch_area invariants.
#       - GetPseudoPatchAge / InsertPatch ordering.
#       - patch_pft_size_profile binning.
#       - DistributeSeeds mass conservation.
#       - split_patch: area + cohort-number + carbon conservation across old/new.
#       - TransLitterNewPatch: litter mass conservation (no burn).
#       - mortality_litter_fluxes / fire_litter_fluxes: killed biomass == sum of
#         litter destinations (+ burn flux).
#       - spawn_patches (treefall): total patch area == site area; cohort number
#         conserved (donor + new) for understory/grass.
#       - fuse_2_patches / fuse_patches: area conservation, relinking.
#       - terminate_patches: a tiny patch is fused away, list relinked.

using Test
using CLM

# --- single woody evergreen carbon PFT (same allometry as the cohort suite) ---
function _setup_edpatch_pft!()
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
    # Patch-dynamics-specific PFT params (mortality litter shunt + LUC fluxes).
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

function _seed_prt_pd(ipft::Int, dbh0::Float64)
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

# Minimal ed_site_type with all the arrays the patch-dynamics paths touch.
function _build_patch_site(npft::Int, nlevsoil::Int)
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

    # impact / fire mortality tallies (spawn survivorship paths).
    site.imort_rate       = zeros(nsc, npft)
    site.imort_carbonflux = zeros(npft)
    site.imort_abg_flux   = zeros(nsc, npft)
    site.fmort_rate_canopy = zeros(nsc, npft)
    site.fmort_rate_ustory = zeros(nsc, npft)
    site.fmort_rate_cambial = zeros(nsc, npft)
    site.fmort_rate_crown   = zeros(nsc, npft)
    site.fmort_carbonflux_canopy = zeros(npft)
    site.fmort_carbonflux_ustory = zeros(npft)
    site.fmort_abg_flux          = zeros(nsc, npft)

    site.nlevsoil     = nlevsoil
    site.zi_soil      = collect(0.0:0.1:0.1 * nlevsoil)
    site.rootfrac_scr = zeros(nlevsoil)

    site.flux_diags = CLM.site_fluxdiags_type()
    site.flux_diags.elem = [CLM.elem_diag_type() for _ in 1:CLM.num_elements[]]
    for el in 1:CLM.num_elements[]
        site.flux_diags.elem[el].surf_fine_litter_input = zeros(npft)
        site.flux_diags.elem[el].root_litter_input      = zeros(npft)
    end

    site.mass_balance = [CLM.site_massbal_type() for _ in 1:CLM.num_elements[]]

    site.disturbance_rates = zeros(CLM.N_DIST_TYPES, CLM.n_landuse_cats, CLM.n_landuse_cats)
    site.landuse_transition_matrix = zeros(CLM.n_landuse_cats, CLM.n_landuse_cats)
    site.landuse_vector_gt_min = fill(false, CLM.n_landuse_cats)
    site.min_allowed_landuse_fraction = 0.0
    site.area_bareground = 0.0
    site.area_PFT = zeros(npft, CLM.n_landuse_cats)

    return site
end

# Create a patch via Create! so it has litter, fuel, and running means.
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

# Add a cohort directly to a patch (carbon-seeded), returning the cohort.
function _add_patch_cohort!(site, patch, ipft, dbh0, nn; clayer=1)
    prt = _seed_prt_pd(ipft, dbh0)
    h, _ = CLM.h_allom(dbh0, ipft)
    CLM.create_cohort(site, patch, ipft, nn, h, 0.0, dbh0, prt, 1.0, 1.0, 1.0,
                      CLM.leaves_on, 0, 1.0, 0.0, clayer, 1, site.spread)
    # the just-inserted cohort is somewhere in the list; mark non-new.
    c = patch.tallest
    while c !== nothing; c.isnew = false; c = c.shorter; end
    return patch.tallest
end

# Walk oldest->youngest counting patches and summing area.
function _site_patch_area(site)
    tot = 0.0
    p = site.oldest_patch
    while p !== nothing; tot += p.area; p = p.younger; end
    return tot
end
function _site_npatches(site)
    n = 0
    p = site.oldest_patch
    while p !== nothing; n += 1; p = p.younger; end
    return n
end

# Sum litter mass [kg] on a patch (= density * area), carbon element.
function _patch_litter_kg(patch)
    l = patch.litter[1]
    dens = sum(l.ag_cwd) + sum(l.bg_cwd) + sum(l.leaf_fines) +
           sum(l.root_fines) + sum(l.seed) + sum(l.seed_germ)
    return dens * patch.area
end

# Number-weighted live carbon [kg] across all cohorts of a patch.
function _patch_live_kg(patch)
    tot = 0.0
    c = patch.tallest
    while c !== nothing
        for org in (CLM.leaf_organ, CLM.fnrt_organ, CLM.sapw_organ,
                    CLM.store_organ, CLM.struct_organ, CLM.repro_organ)
            tot += c.n * CLM.GetState(c.prt, org, CLM.carbon12_element)
        end
        c = c.shorter
    end
    return tot
end

# Total cohort number across a patch list.
function _patch_n(patch)
    tot = 0.0
    c = patch.tallest
    while c !== nothing; tot += c.n; c = c.shorter; end
    return tot
end

# A minimal bc_in for the litter-flux / spawn paths: only max_rooting_depth_index_col
# and the (empty) harvest fields are read on the treefall/fire/no-LUH path.
function _dummy_bc()
    bc = CLM.bc_in_type()
    bc.max_rooting_depth_index_col = 3
    bc.hlm_harvest_rates    = Float64[]
    bc.hlm_harvest_catnames = String[]
    bc.hlm_harvest_units    = CLM.fates_unset_int
    bc.site_area            = CLM.area
    return bc
end

@testset "FATES Batch 14: EDPatchDynamicsMod" begin

    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_sp   = CLM.hlm_use_sp[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_agetrk   = CLM.hlm_use_cohort_age_tracking[]
    old_nocomp   = CLM.hlm_use_nocomp[]
    old_luh      = CLM.hlm_use_luh[]
    old_fixbio   = CLM.hlm_use_fixed_biogeog[]
    old_freqday  = CLM.hlm_freq_day[]
    old_tod      = CLM.hlm_current_tod[]
    old_numpft   = CLM.numpft[]
    old_nharv    = CLM.hlm_num_lu_harvest_cats[]
    old_nleafage = CLM.nleafage[]
    old_nlevsc   = CLM.nlevsclass[]
    old_nlevca   = CLM.nlevcoage[]
    old_sfp      = CLM.SFParams[]

    try
        npft = _setup_edpatch_pft!()
        ipft = 1
        CLM.InitPRTGlobalAllometricCarbon!()

        CLM.num_elements[] = 1
        empty!(CLM.element_list)
        push!(CLM.element_list, CLM.carbon12_element)

        CLM.hlm_parteh_mode[]             = CLM.prt_carbon_allom_hyp
        CLM.hlm_use_sp[]                  = CLM.ifalse
        CLM.hlm_use_planthydro[]          = CLM.ifalse
        CLM.hlm_use_cohort_age_tracking[] = CLM.ifalse
        CLM.hlm_use_nocomp[]              = CLM.ifalse
        CLM.hlm_use_luh[]                 = CLM.ifalse
        CLM.hlm_use_fixed_biogeog[]       = CLM.ifalse
        CLM.hlm_freq_day[]                = 1.0 / 365.0
        CLM.hlm_current_tod[]             = 0
        CLM.numpft[]                      = npft
        CLM.hlm_num_lu_harvest_cats[]     = 1
        CLM.nleafage[]                    = 1

        nlevsoil = 3
        site_area = CLM.area   # 10000.0

        # ------------------------------------------------------------------
        # Build a 3-patch site (oldest -> youngest) summing to the site area.
        # ------------------------------------------------------------------
        site = _build_patch_site(npft, nlevsoil)
        pA = _make_patch(6000.0, npft, nlevsoil)
        pB = _make_patch(3000.0, npft, nlevsoil)
        pC = _make_patch(1000.0, npft, nlevsoil)
        # link oldest(pA) <-> pB <-> youngest(pC)
        pA.older = nothing; pA.younger = pB
        pB.older = pA;      pB.younger = pC
        pC.older = pB;      pC.younger = nothing
        site.oldest_patch   = pA
        site.youngest_patch = pC

        # --- set_patchno / countPatches / check_patch_area invariants ---
        CLM.set_patchno(site)
        @test pA.patchno == 1 && pB.patchno == 2 && pC.patchno == 3
        @test CLM.countPatches(1, [site]) == 3
        @test _site_patch_area(site) ≈ site_area
        # check_patch_area should not throw (areas sum exactly).
        CLM.check_patch_area(site)
        @test _site_patch_area(site) ≈ site_area

        # tiny precision drift gets folded onto the largest patch.
        pA.area = 6000.0 - 1.0e-9
        CLM.check_patch_area(site)
        @test _site_patch_area(site) ≈ site_area rtol=1e-15
        @test pA.area ≈ 6000.0   # largest patch absorbed the drift

        # ------------------------------------------------------------------
        # GetPseudoPatchAge / InsertPatch ordering.
        # ------------------------------------------------------------------
        # All three patches share land use + (unset) nocomp label, so pseudo-age
        # is monotone in age. Set ages so the new patch (age 0) is youngest.
        pA.age = 10.0; pB.age = 5.0; pC.age = 1.0
        @test CLM.GetPseudoPatchAge(pA) > CLM.GetPseudoPatchAge(pC)
        pNew = _make_patch(0.0, npft, nlevsoil)   # area 0 just for ordering test
        pNew.age = 0.0
        CLM.InsertPatch(site, pNew)
        # age 0 is the youngest -> inserted at the youngest head.
        @test site.youngest_patch === pNew
        @test pNew.older === pC
        # remove pNew again (restore the 3-patch site) and re-fix the area.
        site.youngest_patch = pC
        pC.younger = nothing
        pNew.older = nothing

        # ------------------------------------------------------------------
        # patch_pft_size_profile binning. dbh=10 -> bin 2 ((5,20] edge), dbh=30
        # -> bin 3 ((20,50] edge). One cohort each on pB.
        # ------------------------------------------------------------------
        _add_patch_cohort!(site, pB, ipft, 10.0, 0.01)
        _add_patch_cohort!(site, pB, ipft, 30.0, 0.005)
        CLM.patch_pft_size_profile(pB)
        @test sum(pB.pft_agb_profile) > 0.0
        @test pB.pft_agb_profile[ipft, 2] > 0.0   # dbh 10 lands in bin 2
        @test pB.pft_agb_profile[ipft, 3] > 0.0   # dbh 30 lands in bin 3
        @test pB.pft_agb_profile[ipft, 1] == 0.0  # nothing in (0,5]

        # ------------------------------------------------------------------
        # DistributeSeeds: total mass into the seed pools == input mass.
        # seed density added per patch = seed_mass/site_area, mass = density*area,
        # summed over patches -> seed_mass * (sum area)/site_area == seed_mass.
        # ------------------------------------------------------------------
        seed_before = 0.0
        p = site.oldest_patch
        while p !== nothing; seed_before += sum(p.litter[1].seed) * p.area; p = p.younger; end
        seed_mass = 4.2
        CLM.DistributeSeeds(site, seed_mass, 1, ipft)
        seed_after = 0.0
        p = site.oldest_patch
        while p !== nothing; seed_after += sum(p.litter[1].seed) * p.area; p = p.younger; end
        @test (seed_after - seed_before) ≈ seed_mass rtol=1e-12

        # ------------------------------------------------------------------
        # split_patch: split pA (keep 0.7), conserving area + cohort number +
        # live carbon across the old + new patches.
        # ------------------------------------------------------------------
        _add_patch_cohort!(site, pA, ipft, 20.0, 0.02)
        _add_patch_cohort!(site, pA, ipft, 8.0,  0.05)
        area_pre   = pA.area
        n_pre      = _patch_n(pA)
        livec_pre  = _patch_live_kg(pA)
        litt_pre   = _patch_litter_kg(pA)

        new_patch = CLM.fates_patch_type()
        frac_keep = 0.7
        CLM.split_patch(site, pA, new_patch, frac_keep)

        @test pA.area ≈ area_pre * frac_keep rtol=1e-12
        @test new_patch.area ≈ area_pre * (1.0 - frac_keep) rtol=1e-12
        @test (pA.area + new_patch.area) ≈ area_pre rtol=1e-12
        # cohort number conserved (donor keeps frac, new gets 1-frac).
        @test (_patch_n(pA) + _patch_n(new_patch)) ≈ n_pre rtol=1e-12
        # live carbon: donor scaled by frac, new by (1-frac) -> total preserved.
        @test (_patch_live_kg(pA) + _patch_live_kg(new_patch)) ≈ livec_pre rtol=1e-10
        # litter mass conserved across the split (no burn in split_patch).
        @test (_patch_litter_kg(pA) + _patch_litter_kg(new_patch)) ≈ litt_pre rtol=1e-10

        # ------------------------------------------------------------------
        # TransLitterNewPatch in isolation: seed donor with litter, move part of
        # its area into a fresh patch; total litter mass (no fire) is conserved.
        # ------------------------------------------------------------------
        donor = _make_patch(2000.0, npft, nlevsoil)
        # give the donor some litter density.
        donor.litter[1].ag_cwd     .= 0.30
        donor.litter[1].leaf_fines .= 0.10
        donor.litter[1].seed       .= 0.05
        fill!(donor.fuel.frac_burnt, 0.0)   # no burn for a clean conservation test
        areadis = 800.0
        recv = _make_patch(areadis, npft, nlevsoil)
        for el in 1:CLM.num_elements[]
            CLM.InitConditions!(recv.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        end
        donor_litt_before = _patch_litter_kg(donor)   # over full donor area
        CLM.TransLitterNewPatch(site, donor, recv, areadis)
        # donor mass now over its remaining area (it has not shrunk yet here),
        # so account litter on (donor area - areadis) + recv area.
        donor_remaining = donor.area - areadis
        donor_after = (sum(donor.litter[1].ag_cwd) + sum(donor.litter[1].bg_cwd) +
                       sum(donor.litter[1].leaf_fines) + sum(donor.litter[1].root_fines) +
                       sum(donor.litter[1].seed) + sum(donor.litter[1].seed_germ)) * donor_remaining
        recv_after = _patch_litter_kg(recv)
        @test (donor_after + recv_after) ≈ donor_litt_before rtol=1e-10
        # with existing_litt_localization == 1, all donatable litter goes to recv.
        @test recv_after > 0.0

        # ------------------------------------------------------------------
        # mortality_litter_fluxes: killed biomass deposited into litter equals
        # the dead-plant mass (carbon element, no seeds/burn for trees). Use an
        # understory cohort so num_dead = understorey_death * n * areadis/area.
        # ------------------------------------------------------------------
        dpatch = _make_patch(2000.0, npft, nlevsoil)
        ucohort = _add_patch_cohort!(site, dpatch, ipft, 12.0, 0.04; clayer=2)  # understory
        # zero the receiving + donor litter and the site flux tallies.
        recv2 = _make_patch(500.0, npft, nlevsoil)
        for el in 1:CLM.num_elements[]
            CLM.InitConditions!(recv2.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            CLM.InitConditions!(dpatch.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        end
        site.mass_balance[1].burn_flux_to_atm = 0.0

        areadis2 = 500.0
        # expected number killed (understory woody).
        num_dead = CLM.ed_params().ED_val_understorey_death * ucohort.n *
                   (areadis2 / dpatch.area)
        # the dead plant's per-plant carbon mass (allom_frbstor_repro == 0 -> all
        # store goes to root litter, none to a separate seed shunt that would be
        # excluded; with frbstor_repro=0 seed_mass=0, so all dead C -> litter).
        deadc = num_dead * (
            CLM.GetState(ucohort.prt, CLM.leaf_organ,   CLM.carbon12_element) +
            CLM.GetState(ucohort.prt, CLM.fnrt_organ,   CLM.carbon12_element) +
            CLM.GetState(ucohort.prt, CLM.sapw_organ,   CLM.carbon12_element) +
            CLM.GetState(ucohort.prt, CLM.store_organ,  CLM.carbon12_element) +
            CLM.GetState(ucohort.prt, CLM.struct_organ, CLM.carbon12_element) +
            CLM.GetState(ucohort.prt, CLM.repro_organ,  CLM.carbon12_element))

        CLM.mortality_litter_fluxes(site, dpatch, recv2, areadis2, _dummy_bc())

        litt_dest = _patch_litter_kg(recv2) +
            (sum(dpatch.litter[1].ag_cwd) + sum(dpatch.litter[1].bg_cwd) +
             sum(dpatch.litter[1].leaf_fines) + sum(dpatch.litter[1].root_fines) +
             sum(dpatch.litter[1].seed) + sum(dpatch.litter[1].seed_germ)) *
            (dpatch.area - areadis2)
        @test litt_dest ≈ deadc rtol=1e-9
        @test site.mass_balance[1].burn_flux_to_atm == 0.0  # no burn for treefall

        # ------------------------------------------------------------------
        # fire_litter_fluxes: killed + partially-burned biomass == litter + burn
        # flux. fraction_crown_burned = 0 and frac_burnt = 0 -> pure litter, so
        # the killed biomass lands fully in litter and burn flux stays 0.
        # ------------------------------------------------------------------
        fpatch = _make_patch(2000.0, npft, nlevsoil)
        fpatch.fire = CLM.itrue
        fcohort = _add_patch_cohort!(site, fpatch, ipft, 15.0, 0.03; clayer=1)
        fcohort.fire_mort = 0.5
        fcohort.fraction_crown_burned = 0.0
        fill!(fpatch.fuel.frac_burnt, 0.0)
        recv3 = _make_patch(400.0, npft, nlevsoil)
        for el in 1:CLM.num_elements[]
            CLM.InitConditions!(recv3.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            CLM.InitConditions!(fpatch.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        end
        site.mass_balance[1].burn_flux_to_atm = 0.0

        areadis3 = 400.0
        num_dead_fire = fcohort.fire_mort * fcohort.n * (areadis3 / fpatch.area)
        deadc_fire = num_dead_fire * (
            CLM.GetState(fcohort.prt, CLM.leaf_organ,   CLM.carbon12_element) +
            CLM.GetState(fcohort.prt, CLM.fnrt_organ,   CLM.carbon12_element) +
            CLM.GetState(fcohort.prt, CLM.sapw_organ,   CLM.carbon12_element) +
            CLM.GetState(fcohort.prt, CLM.store_organ,  CLM.carbon12_element) +
            CLM.GetState(fcohort.prt, CLM.struct_organ, CLM.carbon12_element) +
            CLM.GetState(fcohort.prt, CLM.repro_organ,  CLM.carbon12_element))

        CLM.fire_litter_fluxes(site, fpatch, recv3, areadis3, _dummy_bc())

        litt_dest_fire = _patch_litter_kg(recv3) +
            (sum(fpatch.litter[1].ag_cwd) + sum(fpatch.litter[1].bg_cwd) +
             sum(fpatch.litter[1].leaf_fines) + sum(fpatch.litter[1].root_fines)) *
            (fpatch.area - areadis3)
        @test (litt_dest_fire + site.mass_balance[1].burn_flux_to_atm) ≈ deadc_fire rtol=1e-9
        @test site.mass_balance[1].burn_flux_to_atm ≈ 0.0 atol=1e-12

        # ------------------------------------------------------------------
        # fuse_2_patches: fuse a small donor into a recipient; area conserved,
        # donor relinked out of the list.
        # ------------------------------------------------------------------
        fsite = _build_patch_site(npft, nlevsoil)
        q1 = _make_patch(7000.0, npft, nlevsoil); q1.age = 8.0
        q2 = _make_patch(3000.0, npft, nlevsoil); q2.age = 2.0
        q1.older = nothing; q1.younger = q2
        q2.older = q1;      q2.younger = nothing
        fsite.oldest_patch = q1; fsite.youngest_patch = q2
        CLM.set_patchno(fsite)
        area_sum_pre = q1.area + q2.area
        # fuse donor q1 into recipient q2.
        CLM.fuse_2_patches(fsite, q1, q2)
        @test q2.area ≈ area_sum_pre rtol=1e-12
        @test _site_npatches(fsite) == 1
        @test fsite.oldest_patch === q2
        @test fsite.youngest_patch === q2
        @test q2.older === nothing && q2.younger === nothing

        # ------------------------------------------------------------------
        # fuse_patches: two patches with near-identical (empty) biomass profiles
        # fuse into one; total area conserved.
        # ------------------------------------------------------------------
        gsite = _build_patch_site(npft, nlevsoil)
        r1 = _make_patch(5000.0, npft, nlevsoil); r1.age = 6.0
        r2 = _make_patch(5000.0, npft, nlevsoil); r2.age = 4.0
        r1.older = nothing; r1.younger = r2
        r2.older = r1;      r2.younger = nothing
        gsite.oldest_patch = r1; gsite.youngest_patch = r2
        CLM.set_patchno(gsite)
        CLM.fuse_patches(gsite, _dummy_bc())
        @test _site_patch_area(gsite) ≈ site_area rtol=1e-10
        @test _site_npatches(gsite) == 1   # identical empty profiles -> fused

        # ------------------------------------------------------------------
        # terminate_patches: a sub-min_patch_area patch is fused into a neighbor.
        # Three patches; make the OLDEST tiny (not the youngest, which is spared).
        # ------------------------------------------------------------------
        tsite = _build_patch_site(npft, nlevsoil)
        t1 = _make_patch(CLM.min_patch_area * 0.5, npft, nlevsoil); t1.age = 9.0  # tiny, oldest
        t2 = _make_patch(6000.0, npft, nlevsoil); t2.age = 5.0
        t3 = _make_patch(4000.0 - CLM.min_patch_area * 0.5, npft, nlevsoil); t3.age = 1.0
        t1.older = nothing; t1.younger = t2
        t2.older = t1;      t2.younger = t3
        t3.older = t2;      t3.younger = nothing
        tsite.oldest_patch = t1; tsite.youngest_patch = t3
        CLM.set_patchno(tsite)
        npatch_pre = _site_npatches(tsite)
        area_pre_t = _site_patch_area(tsite)
        @test area_pre_t ≈ site_area rtol=1e-10
        CLM.terminate_patches(tsite, _dummy_bc())
        @test _site_npatches(tsite) == npatch_pre - 1    # tiny patch fused away
        @test _site_patch_area(tsite) ≈ site_area rtol=1e-10

        # ------------------------------------------------------------------
        # spawn_patches (treefall): give a single patch a treefall disturbance
        # rate, spawn, and confirm total patch area still equals the site area.
        # ------------------------------------------------------------------
        ssite = _build_patch_site(npft, nlevsoil)
        sp1 = _make_patch(site_area, npft, nlevsoil); sp1.age = 3.0
        sp1.older = nothing; sp1.younger = nothing
        ssite.oldest_patch = sp1; ssite.youngest_patch = sp1
        CLM.set_patchno(ssite)
        # an understory cohort so survivorship preserves some plants -> new patch.
        _add_patch_cohort!(ssite, sp1, ipft, 10.0, 0.02; clayer=2)
        # set a treefall disturbance rate directly (bypassing disturbance_rates).
        fill!(sp1.disturbance_rates, 0.0)
        fill!(sp1.landuse_transition_rates, 0.0)
        sp1.disturbance_rates[CLM.dtype_ifall] = 0.1
        sp1.frac_burnt = 0.0

        CLM.spawn_patches(ssite, _dummy_bc())

        @test _site_patch_area(ssite) ≈ site_area rtol=1e-8   # area conserved
        @test _site_npatches(ssite) >= 1
        # patchno bookkeeping was refreshed by spawn -> check_patch_area passed.
        CLM.check_patch_area(ssite)
        @test _site_patch_area(ssite) ≈ site_area rtol=1e-8

    finally
        CLM.prt_global[]                  = old_global
        CLM.num_elements[]                = old_numel
        empty!(CLM.element_list)
        append!(CLM.element_list, old_ellist)
        CLM.hlm_parteh_mode[]             = old_parteh
        CLM.hlm_use_sp[]                  = old_use_sp
        CLM.hlm_use_planthydro[]          = old_hydro
        CLM.hlm_use_cohort_age_tracking[] = old_agetrk
        CLM.hlm_use_nocomp[]              = old_nocomp
        CLM.hlm_use_luh[]                 = old_luh
        CLM.hlm_use_fixed_biogeog[]       = old_fixbio
        CLM.hlm_freq_day[]                = old_freqday
        CLM.hlm_current_tod[]             = old_tod
        CLM.numpft[]                      = old_numpft
        CLM.hlm_num_lu_harvest_cats[]     = old_nharv
        CLM.nleafage[]                    = old_nleafage
        CLM.nlevsclass[]                  = old_nlevsc
        CLM.nlevcoage[]                   = old_nlevca
        CLM.SFParams[]                    = old_sfp
    end
end
