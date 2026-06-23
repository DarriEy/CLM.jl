# test_fates_inventoryinit.jl
# Tests for FATES Batch 16 (Tier F): FatesInventoryInitMod — inventory-based site
# initialization.  This module reads forest-inventory "type 1" PSS (patch) + CSS
# (cohort) records and assembles them into a FATES site's age-ordered patch list /
# height-ordered cohort lists, reusing the Batch 12–15 demographic constructors.
#
# Strategy (mirrors test_fates_edcohortdynamics.jl / test_fates_edpatchdynamics.jl):
#   * Build a synthetic single-PFT (woody, evergreen, carbon-only) prt_params table,
#     param_derived, the FATES PFT table, ed_params (regeneration + bins + fusion
#     tolerances). Register the carbon-only PARTEH hypothesis (no SP / hydro / age).
#   * Build a minimal-but-complete ed_site_type (the arrays fuse_cohorts /
#     fuse_patches / count_cohorts touch).
#   * Exercise the PSS/CSS parsing + site assembly from IN-MEMORY line records:
#       - assess_inventory_sites / count_inventory_sites: descriptor parsing.
#       - set_inventory_patch_type1!: one patch record -> age/area/age_class.
#       - initialize_site_by_inventory!: end-to-end. Assert
#           * correct patch count + per-patch areas (fraction*AREA),
#           * patches linked youngest..oldest in age order (set_patchno renumber),
#           * correct cohort count + PFT/dbh/density assignment,
#           * cohorts height-sorted within a patch,
#           * total area bookkeeping (sum of patch areas == sum of fractions*AREA).
#       - the pft==0 special case (one cohort per PFT).
#       - write_inventory_type1: round-trips a header + one line per patch/cohort.
#
# The fusion tolerances are set TIGHT so the distinct test patches/cohorts do NOT
# fuse away, letting us assert exact counts.

using Test
using CLM

# --- single woody evergreen carbon PFT (same allometry as the cohort suite) ---
function _setup_inv_pft!(npft::Int)
    p = CLM.prt_params
    CLM.allocate_prt_params!(p, npft, CLM.num_organ_types, 1)  # 1 leaf age class

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

    # Evergreen: no seasonal/stress deciduous -> the "fully flushed" branch.
    p.stress_decid .= 0
    p.season_decid .= 0
    p.evergreen    .= CLM.itrue
    p.phen_fnrt_drop_fraction .= 0.0
    p.phen_stem_drop_fraction .= 0.0
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
    CLM.EDPftvarcon_inst[] = ev

    edp = CLM.ed_params_type()
    edp.regeneration_model = CLM.default_regeneration
    edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
    edp.ED_val_history_sizeclass_bin_edges  = [0.0, 5.0, 10.0, 20.0, 50.0]
    edp.ED_val_history_coageclass_bin_edges = [0.0, 3.0, 6.0, 9.0, 12.0]
    edp.ED_val_history_ageclass_bin_edges   = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    # TIGHT fusion tolerances so distinct test cohorts/patches survive.
    edp.ED_val_cohort_size_fusion_tol = 1.0e-6
    edp.ED_val_cohort_age_fusion_tol  = 1.0e-6
    edp.max_cohort_per_patch          = 100
    edp.ED_val_patch_fusion_tol       = 1.0e-6
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

# Minimal ed_site_type with all the arrays the patch/cohort-fusion paths touch.
function _build_inv_site(npft::Int, nlevsoil::Int)
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

    # Phenology state used by set_inventory_cohort_type1! (evergreen branch).
    site.cstatus      = CLM.phen_cstat_notcold
    site.elong_factor = fill(1.0, CLM.maxpft)

    return site
end

# Minimal bc_in for the fuse / litter-flux paths (mirrors the patch-dynamics
# suite's _dummy_bc): only a few fields are read on the no-LUH/no-harvest path.
function _dummy_bc()
    bc = CLM.bc_in_type()
    bc.max_rooting_depth_index_col = 3
    bc.hlm_harvest_rates    = Float64[]
    bc.hlm_harvest_catnames = String[]
    bc.hlm_harvest_units    = CLM.fates_unset_int
    bc.site_area            = CLM.area
    return bc
end

# Walk youngest -> oldest collecting patches.
function _patch_list(site)
    ps = CLM.fates_patch_type[]
    p = site.youngest_patch
    while p !== nothing
        push!(ps, p)
        p = p.older
    end
    return ps
end

# Walk tallest -> shortest collecting cohorts of a patch.
function _cohort_list(patch)
    cs = CLM.fates_cohort_type[]
    c = patch.tallest
    while c !== nothing
        push!(cs, c)
        c = c.shorter
    end
    return cs
end

@testset "FATES Batch 16: FatesInventoryInitMod" begin

    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_sp   = CLM.hlm_use_sp[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_agetrk   = CLM.hlm_use_cohort_age_tracking[]
    old_nleafage = CLM.nleafage[]
    old_numpft   = CLM.numpft[]
    old_tod      = CLM.hlm_current_tod[]
    old_nlevsc   = CLM.nlevsclass[]
    old_nlevca   = CLM.nlevcoage[]
    old_sfp      = CLM.SFParams[]

    try
        npft = 1
        _setup_inv_pft!(npft)
        CLM.InitPRTGlobalAllometricCarbon!()

        # Carbon-only element registry.
        CLM.num_elements[] = 1
        empty!(CLM.element_list)
        push!(CLM.element_list, CLM.carbon12_element)

        CLM.hlm_parteh_mode[]             = CLM.prt_carbon_allom_hyp
        CLM.hlm_use_sp[]                  = CLM.ifalse
        CLM.hlm_use_planthydro[]          = CLM.ifalse
        CLM.hlm_use_cohort_age_tracking[] = CLM.ifalse
        CLM.nleafage[]                    = 1
        CLM.numpft[]                      = npft
        CLM.hlm_current_tod[]             = 0

        nlevsoil = 3

        # =================================================================
        # 1. Descriptor parsing: count + assess inventory sites.
        # =================================================================
        sitelist = [
            "type lat lon pss css",
            "1  45.0  250.0  site_a.pss  site_a.css",
        ]
        @test CLM.count_inventory_sites(sitelist) == 1
        fmt, pss, css, lat, lon = CLM.assess_inventory_sites(sitelist, 1)
        @test fmt == [1]
        @test pss == ["site_a.pss"]
        @test css == ["site_a.css"]
        @test lat == [45.0]
        @test lon == [250.0]

        # Negative longitude is wrapped to [0,360).
        _, _, _, _, lon2 = CLM.assess_inventory_sites(
            ["hdr", "1 0.0 -110.0 a.pss a.css"], 1)
        @test lon2[1] ≈ 250.0

        # =================================================================
        # 2. set_inventory_patch_type1!: one patch record -> age/area.
        # =================================================================
        np = CLM.fates_patch_type()
        CLM.Create!(np, 0.0, 0.0, CLM.primaryland, CLM.fates_unset_int, CLM.num_swb,
                    npft, nlevsoil, 0, CLM.ed_params().regeneration_model)
        pcur = CLM.LineCursor(["time patch trk age area",
                               "1900.0 patchX 2 12.5 0.25"])
        CLM.skip_header!(pcur)
        pname, pios = CLM.set_inventory_patch_type1!(np, pcur, 1)
        @test pios == 0
        @test pname == "patchX"
        @test np.age == 12.5
        @test np.area ≈ 0.25 * CLM.area     # fraction -> m2 of notional AREA
        @test np.age_class == CLM.get_age_class_index(12.5)

        # =================================================================
        # 3. Full site assembly from in-memory PSS/CSS (the core path).
        #    Two patches of distinct ages (0.5 yr young, 30 yr old); the young
        #    patch has 2 cohorts of different dbh, the old patch has 1.
        # =================================================================
        pss_lines = [
            "time patch trk age area",
            "2000.0 young 1  0.5  0.30",
            "2000.0 old   2  30.0 0.70",
        ]
        css_lines = [
            "time patch dbh height pft nplant",
            "2000.0 young 10.0 -1.0 1 0.20",
            "2000.0 young 25.0 -1.0 1 0.05",
            "2000.0 old   40.0 -1.0 1 0.02",
        ]

        site = _build_inv_site(npft, nlevsoil)
        site.lat = 45.0
        site.lon = 250.0
        bc_in = _dummy_bc()

        total_cohorts = CLM.initialize_site_by_inventory!(site, bc_in, 1,
                                                          pss_lines, css_lines)

        patches = _patch_list(site)
        @test length(patches) == 2

        # Patches ordered youngest -> oldest by age.
        @test patches[1].age == 0.5
        @test patches[2].age == 30.0
        # Renumbered 1..N by set_patchno (youngest first).
        @test patches[1].patchno == 1
        @test patches[2].patchno == 2
        # Linked-list head/tail well-formed.
        @test site.youngest_patch === patches[1]
        @test site.oldest_patch   === patches[2]
        @test site.youngest_patch.younger === nothing
        @test site.oldest_patch.older === nothing

        # Per-patch areas = fraction * AREA.
        @test patches[1].area ≈ 0.30 * CLM.area
        @test patches[2].area ≈ 0.70 * CLM.area
        # Total area bookkeeping: areas sum to (sum of fractions)*AREA == AREA.
        @test (patches[1].area + patches[2].area) ≈ CLM.area

        # Cohort counts.
        young_cohorts = _cohort_list(patches[1])
        old_cohorts   = _cohort_list(patches[2])
        @test length(young_cohorts) == 2
        @test length(old_cohorts) == 1
        @test total_cohorts == 3

        # PFT assignment.
        @test all(c.pft == 1 for c in young_cohorts)
        @test old_cohorts[1].pft == 1

        # dbh assignment (cohorts are height-sorted tallest->shortest; bigger dbh
        # => taller for this allometry).
        @test young_cohorts[1].dbh ≈ 25.0
        @test young_cohorts[2].dbh ≈ 10.0
        @test old_cohorts[1].dbh ≈ 40.0

        # Cohorts sorted by size (descending height) within a patch.
        yheights = [c.height for c in young_cohorts]
        @test issorted(yheights; rev=true)

        # Density assignment: n = nplant * patch.area (absolute m2), per the quirk.
        @test young_cohorts[1].n ≈ 0.05 * patches[1].area
        @test young_cohorts[2].n ≈ 0.20 * patches[1].area
        @test old_cohorts[1].n   ≈ 0.02 * patches[2].area

        # =================================================================
        # 4. pft == 0 SPECIAL CASE: one cohort per PFT (here npft=1 -> 1 cohort,
        #    but n is divided by ncohorts_to_create == numpft).
        # =================================================================
        pss0 = ["time patch trk age area",
                "2000.0 p0 1 5.0 1.0"]
        css0 = ["time patch dbh height pft nplant",
                "2000.0 p0 12.0 -1.0 0 0.10"]
        site0 = _build_inv_site(npft, nlevsoil)
        site0.lat = 45.0; site0.lon = 250.0
        tc0 = CLM.initialize_site_by_inventory!(site0, _dummy_bc(), 1, pss0, css0)
        p0 = _patch_list(site0)[1]
        c0 = _cohort_list(p0)
        @test length(c0) == npft            # one cohort per PFT
        @test tc0 == npft
        @test c0[1].pft == 1
        # n = nplant * area / numpft.
        @test c0[1].n ≈ 0.10 * p0.area / npft

        # =================================================================
        # 5. write_inventory_type1: round-trip the assembled site to lines.
        # =================================================================
        pss_name, css_name, pss_out, css_out = CLM.write_inventory_type1(site)
        @test startswith(pss_name, "pss_out_")
        @test startswith(css_name, "css_out_")
        @test pss_out[1] == "time patch trk age area"
        @test css_out[1] == "time patch dbh height pft nplant"
        # header + one line per patch / cohort.
        @test length(pss_out) == 1 + length(patches)
        @test length(css_out) == 1 + total_cohorts

        # =================================================================
        # 6. Public entry: initialize_sites_by_inventory! via in-memory providers.
        # =================================================================
        pss_map = Dict("site_a.pss" => pss_lines)
        css_map = Dict("site_a.css" => css_lines)
        site2 = _build_inv_site(npft, nlevsoil)
        site2.lat = 45.0; site2.lon = 250.0
        CLM.initialize_sites_by_inventory!(1, [site2], [_dummy_bc()];
            sitelist_lines = ["type lat lon pss css",
                              "1 45.0 250.0 site_a.pss site_a.css"],
            pss_provider = p -> pss_map[p],
            css_provider = c -> css_map[c])
        @test length(_patch_list(site2)) == 2

        # A model site too far from any inventory site aborts.
        site_far = _build_inv_site(npft, nlevsoil)
        site_far.lat = 0.0; site_far.lon = 0.0
        @test_throws Exception CLM.initialize_sites_by_inventory!(1, [site_far], [nothing];
            sitelist_lines = ["type lat lon pss css",
                              "1 45.0 250.0 site_a.pss site_a.css"],
            pss_provider = p -> pss_map[p],
            css_provider = c -> css_map[c])

    finally
        CLM.prt_global[]                   = old_global
        CLM.num_elements[]                 = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        CLM.hlm_parteh_mode[]              = old_parteh
        CLM.hlm_use_sp[]                   = old_use_sp
        CLM.hlm_use_planthydro[]           = old_hydro
        CLM.hlm_use_cohort_age_tracking[]  = old_agetrk
        CLM.nleafage[]                     = old_nleafage
        CLM.numpft[]                       = old_numpft
        CLM.hlm_current_tod[]              = old_tod
        CLM.nlevsclass[]                   = old_nlevsc
        CLM.nlevcoage[]                    = old_nlevca
        CLM.SFParams[]                     = old_sfp
    end
end
