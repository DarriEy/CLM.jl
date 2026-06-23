# test_fates_edcanopystructure.jl
# Tests for FATES Batch 15 (Tier F): EDCanopyStructureMod — the canopy STRUCTURE
# ENGINE that arranges cohorts into discrete canopy layers (overstory/understory)
# so each layer's crown area fits the patch area (Perfect Plasticity
# Approximation), and builds the per-patch canopy-layer x PFT x leaf-layer LAI/SAI
# profiles.
#
# Strategy (mirrors test_fates_edpatchdynamics.jl's setup):
#   * Build a synthetic single-PFT (woody, carbon-only) prt_params table, install
#     param_derived, the FATES PFT table, and ed_params (bins + fusion tolerances +
#     comp-exclusion + canopy-closure-thresh + dinc/dlower_vai + radiation_model).
#   * Register the carbon-only PARTEH hypothesis; flags = carbon, no SP, no plant
#     hydro, no cohort-age-tracking, no nocomp, nleafage=1, ohcomp-excln stochastic.
#   * Build a minimal ed_site_type (demotion/promotion rate tallies, termination
#     tallies, flux_diags, soil layering, mass_balance, disturbance_rates).
#   * Build a patch with Create! (litter/fuel/running means) and add cohorts.
#   * Exercise (invariant-based, hand-computed where possible):
#       - NumPotentialCanopyLayers / CanopyLayerArea / calc_areaindex helpers.
#       - canopy_spread! produces a sane [0,1] scaling factor.
#       - canopy_structure! organizes an over-full set of cohorts into layers so
#         each layer's total crown area <= patch area (the key PPA invariant), and
#         ncl_p is set; total cohort number + crown area conserved across the move.
#       - DemoteFromLayer! / PromoteIntoLayer! conserve cohort number + crown area
#         while moving area between layers.
#       - canopy_summarization! + leaf_area_profile! build LAI profiles whose
#         crown-fraction-weighted sum reproduces the per-patch LAI.

using Test
using CLM

# --- single woody evergreen carbon PFT (same allometry as the cohort suite) ---
function _setup_edcanopy_pft!()
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

    # leaf-N vertical scaling (used by tree_lai's decay coefficient).
    p.leafn_vert_scaler_coeff1 .= 0.00963
    p.leafn_vert_scaler_coeff2 .= 2.43

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
    # Canopy-structure-specific PFT params for update_hlm_dynamics!.
    ev.z0mr    = fill(0.055, npft)
    ev.displar = fill(0.67, npft)
    ev.dleaf   = fill(0.04, npft)
    # Wood-product pool split (read by UpdateHarvestC!; fluxes are zero here).
    ev.harvest_pprod10       = fill(1.0, npft)
    ev.landusechange_pprod10 = fill(1.0, npft)
    CLM.EDPftvarcon_inst[] = ev

    edp = CLM.ed_params_type()
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
    # Canopy-structure parameters:
    edp.ED_val_comp_excln             = 3.0   # stochastic demotion (>= 0)
    edp.ED_val_canopy_closure_thresh  = 0.8
    edp.radiation_model               = 1     # NOT twostr_solver (skip rad elements)
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

# Minimal ed_site_type with all the arrays the canopy-structure paths touch.
function _build_canopy_site(npft::Int, nlevsoil::Int)
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
    site.lat                    = 0.0
    site.lon                    = 0.0
    site.snow_depth             = 0.0

    # demotion/promotion tallies (indexed by size class).
    site.demotion_rate   = zeros(nsc)
    site.promotion_rate  = zeros(nsc)
    site.demotion_carbonflux  = 0.0
    site.promotion_carbonflux = 0.0

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

    return site
end

# Create a patch via Create! so it has litter, fuel, and running means.
function _make_canopy_patch(area::Float64, npft::Int, nlevsoil::Int; land_use_label=CLM.primaryland)
    patch = CLM.fates_patch_type()
    CLM.Create!(patch, 0.0, area, land_use_label, CLM.fates_unset_int, CLM.num_swb,
                npft, nlevsoil, 0, CLM.ed_params().regeneration_model)
    for el in 1:CLM.num_elements[]
        CLM.InitConditions!(patch.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
    patch.canopy_layer_tlai = zeros(Float64, CLM.nclmax)
    patch.total_canopy_area = 0.0
    patch.ncl_p             = 1
    patch.nocomp_pft_label  = CLM.fates_unset_int
    return patch
end

# Add a cohort directly to a patch (carbon-seeded), returning the cohort.
function _add_canopy_cohort!(site, patch, ipft, dbh0, nn; clayer=1)
    prt = _seed_prt(ipft, dbh0)
    h, _ = CLM.h_allom(dbh0, ipft)
    CLM.create_cohort(site, patch, ipft, nn, h, 0.0, dbh0, prt, 1.0, 1.0, 1.0,
                      CLM.leaves_on, 0, 1.0, 0.0, clayer, 1, site.spread)
    c = patch.tallest
    while c !== nothing; c.isnew = false; c = c.shorter; end
    return patch.tallest
end

# Total cohort number across a patch's cohort list.
function _coh_n(patch)
    tot = 0.0
    c = patch.tallest
    while c !== nothing; tot += c.n; c = c.shorter; end
    return tot
end

# Total crown area across a patch's cohort list (both layers).
function _coh_carea(patch)
    tot = 0.0
    c = patch.tallest
    while c !== nothing; tot += c.c_area; c = c.shorter; end
    return tot
end

# Count cohorts in a given canopy layer.
function _coh_count_layer(patch, lyr)
    n = 0
    c = patch.tallest
    while c !== nothing
        c.canopy_layer == lyr && (n += 1)
        c = c.shorter
    end
    return n
end

@testset "FATES Batch 15: EDCanopyStructureMod" begin

    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_sp   = CLM.hlm_use_sp[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_agetrk   = CLM.hlm_use_cohort_age_tracking[]
    old_numpft   = CLM.numpft[]
    old_nleafage = CLM.nleafage[]
    old_nlevsc   = CLM.nlevsclass[]
    old_nlevca   = CLM.nlevcoage[]
    old_sfp      = CLM.SFParams[]

    try
        npft = _setup_edcanopy_pft!()
        ipft = 1
        CLM.InitPRTGlobalAllometricCarbon!()

        CLM.num_elements[] = 1
        empty!(CLM.element_list)
        push!(CLM.element_list, CLM.carbon12_element)
        old_elpos = CLM.element_pos[CLM.carbon12_element]
        CLM.element_pos[CLM.carbon12_element] = 1   # carbon is mass_balance[1]

        CLM.hlm_parteh_mode[]             = CLM.prt_carbon_allom_hyp
        CLM.hlm_use_sp[]                  = CLM.ifalse
        CLM.hlm_use_planthydro[]          = CLM.ifalse
        CLM.hlm_use_cohort_age_tracking[] = CLM.ifalse
        CLM.numpft[]                      = npft
        CLM.nleafage[]                    = 1

        nlevsoil  = 3
        site_area = CLM.area   # 10000.0

        # ==============================================================
        # 1. NumPotentialCanopyLayers / CanopyLayerArea helper correctness.
        # ==============================================================
        site = _build_canopy_site(npft, nlevsoil)
        patch = _make_canopy_patch(site_area, npft, nlevsoil)
        site.oldest_patch   = patch
        site.youngest_patch = patch
        patch.older = nothing; patch.younger = nothing

        # One small cohort in layer 1: only one layer, area below patch area.
        _add_canopy_cohort!(site, patch, ipft, 10.0, 1.0e-4; clayer=1)
        @test CLM.NumPotentialCanopyLayers(patch, site.spread) == 1
        a1 = CLM.CanopyLayerArea(patch, site.spread, 1)
        @test a1 ≈ patch.tallest.c_area
        @test a1 < patch.area
        # Layer 2 is empty -> zero area.
        @test CLM.CanopyLayerArea(patch, site.spread, 2) == 0.0

        # Add a cohort flagged into layer 2 -> two layers reported.
        _add_canopy_cohort!(site, patch, ipft, 8.0, 1.0e-4; clayer=2)
        @test CLM.NumPotentialCanopyLayers(patch, site.spread) == 2

        # ==============================================================
        # 2. canopy_spread! -> a sane [0,1] scaling factor.
        # ==============================================================
        site.spread = 0.5
        CLM.canopy_spread!(site)
        @test 0.0 <= site.spread <= 1.0
        # Tiny crown area is far below closure -> spread increases by inc (0.05).
        @test site.spread ≈ 0.55

        # ==============================================================
        # 3. canopy_structure! organizes an OVER-FULL set of cohorts into
        #    layers so each layer's total crown area <= patch area, and
        #    conserves total cohort number + crown area through the move.
        # ==============================================================
        site2 = _build_canopy_site(npft, nlevsoil)
        site2.spread = 1.0
        patch2 = _make_canopy_patch(site_area, npft, nlevsoil)
        site2.oldest_patch   = patch2
        site2.youngest_patch = patch2
        patch2.older = nothing; patch2.younger = nothing

        # Several large cohorts all initially in layer 1 -> total crown area
        # (~97 m2/plant at dbh 50, x large n) far exceeds the patch area,
        # forcing demotion into layer 2.
        for (dbh0, nn) in [(50.0, 40.0), (45.0, 40.0), (40.0, 40.0), (35.0, 40.0)]
            _add_canopy_cohort!(site2, patch2, ipft, dbh0, nn; clayer=1)
        end

        n_before     = _coh_n(patch2)
        carea_before = CLM.CanopyLayerArea(patch2, site2.spread, 1)  # recomputes c_area
        carea_total_before = _coh_carea(patch2)
        @test carea_before > patch2.area   # genuinely over-full to begin with

        CLM.canopy_structure!(site2, nothing)

        # KEY INVARIANT: each occupied layer's crown area <= patch area (tol).
        z = CLM.NumPotentialCanopyLayers(patch2, site2.spread)
        @test patch2.ncl_p == z
        @test z >= 2   # the over-full top layer must have spilled into layer 2
        for ly in 1:z
            la = CLM.CanopyLayerArea(patch2, site2.spread, ly)
            @test la - patch2.area <= 1.0e-6
        end
        # Top layer should be (essentially) full.
        @test CLM.CanopyLayerArea(patch2, site2.spread, 1) ≈ patch2.area rtol=1e-6

        # Conservation: total number + total crown area unchanged by the
        # demotion/splitting (no termination happened: all cohorts are robust).
        @test _coh_n(patch2) ≈ n_before rtol=1e-9
        @test _coh_carea(patch2) ≈ carea_total_before rtol=1e-6

        # ==============================================================
        # 4. DemoteFromLayer! conserves cohort number + crown area while
        #    moving area from layer 1 into layer 2.
        # ==============================================================
        site3 = _build_canopy_site(npft, nlevsoil)
        site3.spread = 1.0
        patch3 = _make_canopy_patch(site_area, npft, nlevsoil)
        site3.oldest_patch = patch3; site3.youngest_patch = patch3
        for (dbh0, nn) in [(50.0, 50.0), (44.0, 50.0), (38.0, 50.0)]
            _add_canopy_cohort!(site3, patch3, ipft, dbh0, nn; clayer=1)
        end
        n_pre     = _coh_n(patch3)
        carea_pre = (CLM.CanopyLayerArea(patch3, site3.spread, 1); _coh_carea(patch3))
        layer1_pre = CLM.CanopyLayerArea(patch3, site3.spread, 1)
        @test layer1_pre > patch3.area

        CLM.DemoteFromLayer!(site3, patch3, 1, nothing)

        # Layer 1 trimmed to the patch area.
        @test CLM.CanopyLayerArea(patch3, site3.spread, 1) ≈ patch3.area rtol=1e-6
        # The excess landed in layer 2.
        @test CLM.CanopyLayerArea(patch3, site3.spread, 2) > 0.0
        # Number + total crown area conserved across the split/demotion.
        @test _coh_n(patch3) ≈ n_pre rtol=1e-9
        @test _coh_carea(patch3) ≈ carea_pre rtol=1e-6

        # ==============================================================
        # 5. PromoteIntoLayer! conserves number + crown area, filling an
        #    under-full top layer from the layer below.
        # ==============================================================
        site4 = _build_canopy_site(npft, nlevsoil)
        site4.spread = 1.0
        patch4 = _make_canopy_patch(site_area, npft, nlevsoil)
        site4.oldest_patch = patch4; site4.youngest_patch = patch4
        # Top layer under-full (~4850 m2 < 10000); understory (~12500 m2) has
        # plenty to promote upward to fill the gap.
        _add_canopy_cohort!(site4, patch4, ipft, 50.0, 50.0; clayer=1)  # partial top
        _add_canopy_cohort!(site4, patch4, ipft, 45.0, 80.0; clayer=2)  # understory
        _add_canopy_cohort!(site4, patch4, ipft, 40.0, 80.0; clayer=2)

        layer1_under = CLM.CanopyLayerArea(patch4, site4.spread, 1)
        @test layer1_under < patch4.area
        n_pre4     = _coh_n(patch4)
        carea_pre4 = _coh_carea(patch4)

        CLM.PromoteIntoLayer!(site4, patch4, 1)

        # Top layer brought up to (essentially) full.
        @test CLM.CanopyLayerArea(patch4, site4.spread, 1) ≈ patch4.area rtol=1e-6
        @test _coh_n(patch4) ≈ n_pre4 rtol=1e-9
        @test _coh_carea(patch4) ≈ carea_pre4 rtol=1e-6

        # ==============================================================
        # 6. canopy_summarization! + leaf_area_profile! build LAI profiles
        #    whose crown-fraction-weighted sum reproduces patch LAI; and
        #    calc_areaindex integrates them consistently.
        # ==============================================================
        site5 = _build_canopy_site(npft, nlevsoil)
        site5.spread = 1.0
        patch5 = _make_canopy_patch(site_area, npft, nlevsoil)
        patch5.nocomp_pft_label = CLM.fates_unset_int  # not bareground
        site5.oldest_patch = patch5; site5.youngest_patch = patch5
        # A single, robust top-layer cohort that fits inside the patch area.
        _add_canopy_cohort!(site5, patch5, ipft, 50.0, 0.02; clayer=1)
        # Organize layers first (sets ncl_p, c_area, total_canopy_area downstream).
        CLM.canopy_structure!(site5, nothing)

        bc_in = [CLM.bc_in_type()]
        CLM.canopy_summarization!(1, [site5], bc_in)

        @test patch5.total_canopy_area > 0.0
        @test patch5.total_canopy_area <= patch5.area + 1.0e-6
        @test patch5.ncl_p >= 1

        # The leaf-area profile must sum (crown-fraction weighted) back to the
        # patch total LAI computed directly from the cohort tree LAI.
        # patch tlai = sum_cohort treelai * c_area / total_canopy_area.
        direct_tlai = 0.0
        c = patch5.tallest
        while c !== nothing
            direct_tlai += c.treelai * c.c_area / patch5.total_canopy_area
            c = c.shorter
        end

        # Integrated profile total: sum over cl,ft,iv of
        #   canopy_area_profile * tlai_profile  (this is exactly calc_areaindex
        #   before the ai_min floor — with a single full top layer it equals the
        #   direct per-patch LAI).
        prof_tlai = 0.0
        for cl in 1:patch5.ncl_p
            for ft in 1:CLM.numpft[]
                for iv in 1:patch5.nrad[cl, ft]
                    prof_tlai += patch5.canopy_area_profile[cl, ft, iv] *
                                 patch5.tlai_profile[cl, ft, iv]
                end
            end
        end
        @test prof_tlai ≈ direct_tlai rtol=1e-8
        @test prof_tlai > 0.0

        # calc_areaindex applies a 0.1 floor then matches the integration.
        ai_tlai = CLM.calc_areaindex(patch5, "tlai")
        @test ai_tlai ≈ max(0.1, prof_tlai) rtol=1e-8
        @test CLM.calc_areaindex(patch5, "elai") <= ai_tlai + 1e-12  # ELAI <= TLAI (snow-free here so equal)
        @test CLM.calc_areaindex(patch5, "tsai") >= 0.1              # at least the floor
        @test_throws Exception CLM.calc_areaindex(patch5, "bogus")

        # ==============================================================
        # 7. UpdatePatchLAI! / UpdateCohortLAI! produce positive cohort LAI
        #    and a sensible vegetation-layer count.
        # ==============================================================
        c = patch5.tallest
        @test c.treelai > 0.0
        @test c.treesai >= 0.0
        @test c.nv >= 1
        # nleaf for the occupied (cl=1, ft=1) bin matches the cohort NV.
        @test patch5.nleaf[1, ipft] >= c.nv

        # ==============================================================
        # 8. update_hlm_dynamics! packs sane boundary-condition outputs.
        # ==============================================================
        bc_out = [CLM.bc_out_type()]
        np = 1   # one vegetated patch
        bc_out[1].canopy_fraction_pa    = zeros(np)
        bc_out[1].dleaf_pa              = zeros(np)
        bc_out[1].z0m_pa                = zeros(np)
        bc_out[1].displa_pa             = zeros(np)
        bc_out[1].htop_pa               = zeros(np)
        bc_out[1].hbot_pa               = zeros(np)
        bc_out[1].nocomp_pft_label_pa   = zeros(Int, np)
        bc_out[1].elai_pa               = zeros(np)
        bc_out[1].esai_pa               = zeros(np)
        bc_out[1].tlai_pa               = zeros(np)
        bc_out[1].tsai_pa               = zeros(np)
        bc_out[1].frac_veg_nosno_alb_pa = zeros(np)
        # UpdateHarvestC! writes the wood-product flux fields (zero flux here).
        bc_out[1].hrv_deadstemc_to_prod10c  = 0.0
        bc_out[1].hrv_deadstemc_to_prod100c = 0.0

        CLM.update_hlm_dynamics!(1, [site5], [1], bc_out)

        @test bc_out[1].htop_pa[1] ≈ patch5.tallest.height
        @test bc_out[1].z0m_pa[1] > 0.0
        @test bc_out[1].displa_pa[1] > 0.0
        @test bc_out[1].dleaf_pa[1] ≈ 0.04 rtol=1e-9   # single PFT -> its dleaf
        @test 0.0 <= bc_out[1].canopy_fraction_pa[1] <= 1.0
        @test bc_out[1].frac_veg_nosno_alb_pa[1] == 1.0
        @test bc_out[1].tlai_pa[1] >= bc_out[1].elai_pa[1] - 1e-12

    finally
        CLM.prt_global[]                  = old_global
        CLM.num_elements[]                = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        CLM.hlm_parteh_mode[]             = old_parteh
        CLM.hlm_use_sp[]                  = old_use_sp
        CLM.hlm_use_planthydro[]          = old_hydro
        CLM.hlm_use_cohort_age_tracking[] = old_agetrk
        CLM.numpft[]                      = old_numpft
        CLM.nleafage[]                    = old_nleafage
        CLM.nlevsclass[]                  = old_nlevsc
        CLM.nlevcoage[]                   = old_nlevca
        CLM.SFParams[]                    = old_sfp
        if @isdefined(old_elpos); CLM.element_pos[CLM.carbon12_element] = old_elpos; end
    end
end
