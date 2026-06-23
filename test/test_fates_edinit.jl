# test_fates_edinit.jl
# Tests for FATES Batch 17 (Tier F): EDInitMod — FATES cold-start state
# initialization.  This is the DEFAULT (non-inventory) cold-start that allocates
# + zeroes the site-level arrays (init_site_vars! / zero_site!), sets the initial
# site properties (set_site_properties!), builds the starting patch(es)
# (init_patches!) and seeds each with one small cohort per PFT (init_cohorts!).
#
# Strategy (mirrors test_fates_inventoryinit.jl):
#   * Build a synthetic 2-PFT (woody, evergreen, carbon-only) prt_params table,
#     param_derived, EDPftvarcon (initd / hgt_min / vcmax25top / hlm_pft_map),
#     ed_params (regeneration + bins + TIGHT fusion tols). Register the carbon-only
#     PARTEH hypothesis.  Set the HLM globals for the default cold-start path:
#     no inventory / no nocomp / no fixed-biogeog / no LUH / no SP / no hydro.
#   * init_site_vars! -> the size-x-pft / soil / flux-diag / seed arrays are
#     allocated and the (Nesterov) fire-weather object is attached.
#   * zero_site! -> the accumulators / fire / disturbance / mass-balance / flux
#     diagnostics are fully zeroed, phenology/elong unset, spread 0.
#   * set_site_properties! -> sane cold-start phenology dates/status/GDD, water +
#     veg-temp memory, recruit l2fr, use_this_pft all on, min landuse fraction.
#   * init_patches! (1 site) -> exactly ONE primaryland patch with the full site
#     area; init_cohorts! seeds it with one cohort per PFT at the prescribed
#     initial small size, height-sorted, with the mass-balance old_stock set.
#
# The fusion tolerances are set TIGHT so the two distinct-PFT cohorts do NOT fuse.

using Test
using CLM

# --- two woody evergreen carbon PFTs (same allometry as the cohort/inv suites) ---
function _setup_edinit_pft!(npft::Int)
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

    # Evergreen: the "fully flushed" cohort branch.
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
    # Cold-start init knobs read by init_cohorts!:
    ev.initd   = fill(0.2, npft)   # positive => initial recruit density
    ev.hgt_min = fill(1.25, npft)  # sapling height (m)
    ev.hlm_pft_map = [Float64(i == j) for i in 1:npft, j in 1:npft]  # identity HLM->FATES map
    CLM.EDPftvarcon_inst[] = ev

    edp = CLM.ed_params_type()
    edp.regeneration_model = CLM.default_regeneration
    edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
    edp.ED_val_history_sizeclass_bin_edges  = [0.0, 5.0, 10.0, 20.0, 50.0]
    edp.ED_val_history_coageclass_bin_edges = [0.0, 3.0, 6.0, 9.0, 12.0]
    edp.ED_val_history_ageclass_bin_edges   = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    edp.ED_val_cohort_size_fusion_tol = 1.0e-6
    edp.ED_val_cohort_age_fusion_tol  = 1.0e-6
    edp.max_cohort_per_patch          = 100
    edp.ED_val_patch_fusion_tol       = 1.0e-6
    edp.fates_mortality_disturbance_fraction = 1.0
    edp.ED_val_understorey_death      = 0.55983
    edp.logging_coll_under_frac       = 0.55983
    edp.maxpatches_by_landuse         = fill(10, CLM.n_landuse_cats)
    edp.max_nocomp_pfts_by_landuse    = fill(npft, CLM.n_landuse_cats)
    edp.crop_lu_pft_vector            = fill(1, CLM.n_landuse_cats)
    nlv = CLM.nlevleaf
    edp.dinc_vai   = fill(1.0, nlv)
    edp.dlower_vai = collect(1.0:1.0:Float64(nlv))
    CLM.EDParams[] = edp

    CLM.nlevsclass[] = length(edp.ED_val_history_sizeclass_bin_edges)
    CLM.nlevcoage[]  = length(edp.ED_val_history_coageclass_bin_edges)
    CLM.nlevage[]    = length(edp.ED_val_history_ageclass_bin_edges)
    CLM.nlevdamage[] = length(edp.ED_val_history_damage_bin_edges)

    sfp = CLM.sf_params_type()
    sfp.SF_val_CWD_frac = [0.045, 0.075, 0.21, 0.67]
    CLM.SFParams[] = sfp

    return npft
end

# A bc_in with the soil layering init_site_vars! reads + the fuse/litter paths use.
function _edinit_bc(nlevsoil::Int)
    bc = CLM.bc_in_type()
    CLM.allocate_bcin!(bc; npatches=1, nlevsoil=nlevsoil, nlevdecomp=nlevsoil)
    bc.nlevsoil = nlevsoil
    # Sane soil layering for the static-array copy in init_site_vars!.
    bc.zi_sisl  .= collect(0.0:0.1:0.1 * nlevsoil)
    bc.dz_sisl  .= fill(0.1, nlevsoil)
    bc.z_sisl   .= collect(0.05:0.1:0.1 * nlevsoil)
    bc.max_rooting_depth_index_col = 3
    bc.hlm_harvest_rates    = Float64[]
    bc.hlm_harvest_catnames = String[]
    bc.hlm_harvest_units    = CLM.fates_unset_int
    bc.site_area            = CLM.area
    return bc
end

# Walk youngest -> oldest collecting patches.
function _edinit_patch_list(site)
    ps = CLM.fates_patch_type[]
    p = site.youngest_patch
    while p !== nothing
        push!(ps, p)
        p = p.older
    end
    return ps
end

# Walk tallest -> shortest collecting cohorts of a patch.
function _edinit_cohort_list(patch)
    cs = CLM.fates_cohort_type[]
    c = patch.tallest
    while c !== nothing
        push!(cs, c)
        c = c.shorter
    end
    return cs
end

@testset "FATES Batch 17: EDInitMod" begin

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
    old_nlevage  = CLM.nlevage[]
    old_nlevdam  = CLM.nlevdamage[]
    old_sfp      = CLM.SFParams[]
    old_restart  = CLM.hlm_is_restart[]
    old_nocomp   = CLM.hlm_use_nocomp[]
    old_biogeog  = CLM.hlm_use_fixed_biogeog[]
    old_luh      = CLM.hlm_use_luh[]
    old_inv      = CLM.hlm_use_inventory_init[]
    old_damage   = CLM.hlm_use_tree_damage[]
    old_doy      = CLM.hlm_day_of_year[]

    try
        npft = 2
        _setup_edinit_pft!(npft)
        CLM.InitPRTGlobalAllometricCarbon!()

        # Carbon-only element registry.
        CLM.num_elements[] = 1
        empty!(CLM.element_list)
        push!(CLM.element_list, CLM.carbon12_element)

        # Default cold-start path globals.
        CLM.hlm_parteh_mode[]             = CLM.prt_carbon_allom_hyp
        CLM.hlm_use_sp[]                  = CLM.ifalse
        CLM.hlm_use_planthydro[]          = CLM.ifalse
        CLM.hlm_use_cohort_age_tracking[] = CLM.ifalse
        CLM.hlm_use_nocomp[]              = CLM.ifalse
        CLM.hlm_use_fixed_biogeog[]       = CLM.ifalse
        CLM.hlm_use_luh[]                 = CLM.ifalse
        CLM.hlm_use_inventory_init[]      = CLM.ifalse
        CLM.hlm_use_tree_damage[]         = CLM.ifalse
        CLM.hlm_is_restart[]              = CLM.ifalse
        CLM.nleafage[]                    = 1
        CLM.numpft[]                      = npft
        CLM.hlm_current_tod[]             = 0
        CLM.hlm_day_of_year[]             = 1

        nlevsoil = 3
        bc = _edinit_bc(nlevsoil)

        # =================================================================
        # 1. init_site_vars!: arrays allocated, fire-weather attached.
        # =================================================================
        site = CLM.ed_site_type()
        CLM.init_site_vars!(site, bc)

        nsc = CLM.nlevsclass[]
        nmt = CLM.n_term_mort_types
        @test size(site.term_nindivs_canopy) == (nmt, nsc, npft)
        @test size(site.imort_rate) == (nsc, npft)
        @test size(site.area_PFT) == (npft, CLM.n_landuse_cats)
        @test length(site.use_this_pft) == npft
        @test length(site.mass_balance) == CLM.num_elements[]
        @test length(site.iflux_balance) == CLM.num_elements[]
        @test length(site.flux_diags.elem) == CLM.num_elements[]
        @test site.nlevsoil == nlevsoil
        @test length(site.zi_soil) == nlevsoil + 1
        @test length(site.dz_soil) == nlevsoil
        # Static soil arrays copied from the boundary condition.
        @test site.dz_soil == bc.dz_sisl
        @test site.z_soil  == bc.z_sisl
        # Fire weather object attached (Nesterov).
        @test site.fireWeather isa CLM.nesterov_index

        # =================================================================
        # 2. zero_site!: accumulators zeroed, phenology/elong unset, spread 0.
        # =================================================================
        CLM.zero_site!(site)
        @test site.oldest_patch === nothing
        @test site.youngest_patch === nothing
        @test site.cstatus == CLM.fates_unset_int
        @test all(site.dstatus[1:npft] .== CLM.fates_unset_int)
        @test all(isnan, site.elong_factor[1:npft])
        @test site.spread == 0.0
        @test site.NF == 0.0
        @test site.NF_successful == 0.0
        @test all(site.term_nindivs_canopy .== 0.0)
        @test all(site.imort_rate .== 0.0)
        @test all(site.area_PFT .== 0.0)
        @test site.area_bareground == 0.0
        @test site.transition_landuse_from_off_to_on == false
        @test site.mass_balance[1].old_stock == 0.0

        # =================================================================
        # 3. set_site_properties!: sane cold-start initial values.
        # =================================================================
        sites = [site]
        bc_in = [bc]
        CLM.set_site_properties!(1, sites, bc_in)

        @test site.nchilldays == 0
        @test site.ncolddays  == 0
        @test site.phen_model_date == 0
        @test site.grow_deg_days == 30.0
        # cleafondate = cleafon(100) - day_of_year(1).
        @test site.cleafondate  == 100 - 1
        @test site.cleafoffdate == 300 - 1
        @test site.cstatus == CLM.phen_cstat_notcold
        @test all(site.dstatus[1:npft] .== CLM.phen_dstat_moiston)
        @test all(site.elong_factor[1:npft] .== 1.0)
        # Water + veg-temp memory.
        @test all(site.liqvol_memory[1:CLM.numWaterMem, 1:npft] .== 0.5)
        @test all(site.smp_memory[1:CLM.numWaterMem, 1:npft] .== 0.0)
        @test all(site.vegtemp_memory[1:CLM.num_vegtemp_mem] .== 0.0)
        @test site.ema_npp == -9999.0
        # No fixed-biogeog => every PFT is used.
        @test all(site.use_this_pft[1:npft] .== CLM.itrue)
        # Recruit l2fr seeded from allom_l2fr.
        @test all(all(site.rec_l2fr[ft, :] .== CLM.prt_params.allom_l2fr[ft]) for ft in 1:npft)
        # min allowed landuse fraction (non-nocomp branch).
        @test site.min_allowed_landuse_fraction ≈ CLM.min_patch_area_forced / CLM.area

        # A restart leaves the values untouched (no-op).
        site_r = CLM.ed_site_type()
        CLM.init_site_vars!(site_r, bc)
        CLM.zero_site!(site_r)
        CLM.hlm_is_restart[] = CLM.itrue
        CLM.set_site_properties!(1, [site_r], [bc])
        @test site_r.cstatus == CLM.fates_unset_int   # untouched
        @test site_r.grow_deg_days |> isnan
        CLM.hlm_is_restart[] = CLM.ifalse

        # =================================================================
        # 4. init_patches!: one primaryland patch with the full site area,
        #    seeded with one cohort per PFT (init_cohorts!).
        # =================================================================
        CLM.init_patches!(1, sites, bc_in)

        patches = _edinit_patch_list(site)
        @test length(patches) == 1                       # default: 1 patch
        p1 = patches[1]
        @test p1.land_use_label == CLM.primaryland
        @test p1.patchno == 1
        @test site.youngest_patch === p1
        @test site.oldest_patch === p1

        # The single patch holds the whole site area.
        @test p1.area ≈ CLM.area
        # Patch areas sum to the site area (total-area accounting).
        total = sum(p.area for p in patches)
        @test total ≈ CLM.area

        # init_cohorts! seeded one cohort per PFT at the prescribed small size.
        cohorts = _edinit_cohort_list(p1)
        @test length(cohorts) == npft
        @test sort([c.pft for c in cohorts]) == collect(1:npft)
        # Each cohort starts at the prescribed initial density n = initd * area.
        for c in cohorts
            @test c.n ≈ CLM.EDPftvarcon_inst[].initd[c.pft] * p1.area
            # Small initial size: dbh from hgt_min via allometry => positive, modest.
            @test c.dbh > 0.0
            @test c.height > 0.0
            @test c.height ≈ CLM.EDPftvarcon_inst[].hgt_min[c.pft] rtol=1e-6
        end

        # Cohorts are height-sorted (tallest -> shortest) within the patch.
        heights = [c.height for c in cohorts]
        @test issorted(heights; rev=true)

        # Mass-balance stock initialized to the live + litter + seed stock.
        @test site.mass_balance[1].old_stock > 0.0
        # Patch fire variables zeroed for the first timestep.
        @test p1.livegrass == 0.0
        @test p1.fire == 0
        @test all(p1.scorch_ht .== 0.0)

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
        CLM.nlevage[]                      = old_nlevage
        CLM.nlevdamage[]                   = old_nlevdam
        CLM.SFParams[]                     = old_sfp
        CLM.hlm_is_restart[]               = old_restart
        CLM.hlm_use_nocomp[]               = old_nocomp
        CLM.hlm_use_fixed_biogeog[]        = old_biogeog
        CLM.hlm_use_luh[]                  = old_luh
        CLM.hlm_use_inventory_init[]       = old_inv
        CLM.hlm_use_tree_damage[]          = old_damage
        CLM.hlm_day_of_year[]              = old_doy
    end
end
