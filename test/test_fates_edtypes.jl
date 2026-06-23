# test_fates_edtypes.jl
# Tests for FATES Batch 10 (Tier F): EDTypesMod — the SITE type, the TOP of the
# FATES type system. ed_site_type owns the age-ordered patch linked list
# (oldest_patch/youngest_patch), the site plant-hydraulics object, the
# fire-weather object, the per-element mass-balance + integrated-flux-balance
# arrays, the flux diagnostics, the full phenology state, soil layering, and the
# termination/recruitment/demotion/promotion/disturbance accumulators.
#
# Strategy:
#   * Construct an ed_site_type with a small synthetic dimension set; assert
#     scalar field init, fixed-size array shapes (sized by maxpft / num_vegtemp_mem
#     / numWaterMem / N_DIST_TYPES / n_landuse_cats), and that the helper objects
#     (resources_management, flux_diags) were default-constructed.
#   * Allocate + zero the helper types via the bang routines (ZeroMassBalState!,
#     ZeroMassBalFlux!, ZeroFluxDiags!).
#   * Build a 3-patch oldest_patch/youngest_patch age-ordered linked list (each
#     patch carrying a cohort or two), and assert the age ordering traverses
#     oldest -> youngest and the per-land-use / tree-grass-area site routines.

using Test
using CLM

@testset "FATES Batch 10: EDTypesMod (ed_site_type)" begin

    # ----------------------------------------------------------------------
    # Module-level constants owned by EDTypesMod.
    # ----------------------------------------------------------------------
    @test CLM.init_recruit_trim == 0.8
    @test CLM.n_rad_stream_types == 2
    @test CLM.idirect == 1 && CLM.idiffuse == 2
    @test CLM.do_fates_salinity == false
    @test CLM.area == 10000.0
    @test CLM.area_inv == 1.0e-4
    @test CLM.numWaterMem == 10
    @test CLM.num_vegtemp_mem == 10
    @test CLM.phen_cstat_nevercold == 0
    @test CLM.phen_dstat_pshed == 4
    @test CLM.HEIGHTMAX == 30.0
    @test CLM.min_nppatch == CLM.min_npm2 * CLM.min_patch_area
    @test CLM.homogenize_seed_pfts == false

    # ----------------------------------------------------------------------
    # Helper types default-construct + zero.
    # ----------------------------------------------------------------------
    mb = CLM.site_massbal_type()
    mb.old_stock = 9.0; mb.gpp_acc = 1.0
    fill!(mb.wood_product_harvest, 2.0)
    @test length(mb.wood_product_harvest) == CLM.maxpft
    @test length(mb.wood_product_landusechange) == CLM.maxpft
    CLM.ZeroMassBalState!(mb)
    @test mb.old_stock == 0.0 && mb.err_fates == 0.0
    @test mb.gpp_acc == 1.0  # flux not touched by state-zero
    CLM.ZeroMassBalFlux!(mb)
    @test mb.gpp_acc == 0.0
    @test all(mb.wood_product_harvest .== 0.0)

    # site_fluxdiags_type with a per-element elem array.
    fd = CLM.site_fluxdiags_type()
    @test length(fd.elem) == 0
    e1 = CLM.elem_diag_type()
    @test length(e1.cwd_ag_input) == CLM.ncwd
    @test length(e1.cwd_bg_input) == CLM.ncwd
    e1.surf_fine_litter_input = zeros(3)
    e1.root_litter_input = zeros(3)
    push!(fd.elem, e1)
    fd.npp = 5.0; fd.nh4_uptake = 2.0
    fd.elem[1].cwd_ag_input[1] = 7.0
    fd.elem[1].burned_liveveg = 3.0
    # ZeroFluxDiags! loops over num_elements[]; set it to 1 for this test.
    old_numel = CLM.num_elements[]
    try
        CLM.num_elements[] = 1
        CLM.ZeroFluxDiags!(fd)
    finally
        CLM.num_elements[] = old_numel
    end
    @test fd.npp == 0.0 && fd.nh4_uptake == 0.0
    @test all(fd.elem[1].cwd_ag_input .== 0.0)
    @test fd.elem[1].burned_liveveg == 0.0
    @test all(fd.elem[1].surf_fine_litter_input .== 0.0)

    # ----------------------------------------------------------------------
    # Construct the site with a synthetic dimension set.
    # ----------------------------------------------------------------------
    site = CLM.ed_site_type()

    # Pointers default to nothing.
    @test site.oldest_patch === nothing
    @test site.youngest_patch === nothing
    @test site.si_hydr === nothing
    @test site.fireWeather === nothing

    # Default scalars / sentinels.
    @test isnan(site.lat) && isnan(site.lon)
    @test site.h_gid == CLM.fates_unset_int
    @test site.cstatus == CLM.fates_unset_int
    @test site.nlevsoil == CLM.fates_unset_int
    @test site.transition_landuse_from_off_to_on == false

    # Default-constructed helper objects.
    @test site.resources_management isa CLM.ed_resources_management_type
    @test site.flux_diags isa CLM.site_fluxdiags_type

    # Fixed-size array shapes.
    @test length(site.dstatus) == CLM.maxpft
    @test length(site.dleafondate) == CLM.maxpft
    @test length(site.elong_factor) == CLM.maxpft
    @test length(site.recruitment_rate) == CLM.maxpft
    @test length(site.vegtemp_memory) == CLM.num_vegtemp_mem
    @test size(site.liqvol_memory) == (CLM.numWaterMem, CLM.maxpft)
    @test size(site.smp_memory) == (CLM.numWaterMem, CLM.maxpft)
    @test size(site.disturbance_rates) == (CLM.N_DIST_TYPES, CLM.n_landuse_cats, CLM.n_landuse_cats)
    @test size(site.landuse_transition_matrix) == (CLM.n_landuse_cats, CLM.n_landuse_cats)

    # Allocatable arrays default to "not allocated" (zero-sized).
    @test isempty(site.area_by_age)
    @test isempty(site.seed_in) && isempty(site.seed_out)
    @test isempty(site.term_carbonflux_canopy)
    @test isempty(site.fmort_rate_canopy_damage)

    # Attach the per-element mass-balance / iflux-balance arrays.
    site.mass_balance = [CLM.site_massbal_type()]
    site.iflux_balance = [CLM.site_ifluxbal_type()]
    @test length(site.mass_balance) == 1
    @test length(site.iflux_balance) == 1

    # ----------------------------------------------------------------------
    # Build a 3-patch oldest -> youngest linked list, each with a cohort or two.
    # Use Create! for realistic patch construction. PFT1 woody, PFT2 grass.
    # ----------------------------------------------------------------------
    old_numpft  = CLM.numpft[]
    old_ne      = CLM.num_elements[]
    old_ellist  = copy(CLM.element_list)
    old_woody   = copy(CLM.prt_params.woody)
    try
        npft = 2
        CLM.numpft[] = npft
        CLM.num_elements[] = 1
        empty!(CLM.element_list); push!(CLM.element_list, CLM.carbon12_element)
        CLM.prt_params.woody = [CLM.itrue, CLM.ifalse]

        num_swb = CLM.num_swb
        nlsoil  = 4
        regen   = CLM.default_regeneration

        # ages: oldest (200, secondary) -> mid (50, primary) -> youngest (5, secondary)
        pOld = CLM.fates_patch_type()
        CLM.Create!(pOld, 200.0, 3000.0, CLM.secondaryland, CLM.fates_unset_int,
                    num_swb, npft, nlsoil, 0, regen)
        pMid = CLM.fates_patch_type()
        CLM.Create!(pMid, 50.0, 4000.0, CLM.primaryland, CLM.fates_unset_int,
                    num_swb, npft, nlsoil, 0, regen)
        pYng = CLM.fates_patch_type()
        CLM.Create!(pYng, 5.0, 3000.0, CLM.secondaryland, CLM.fates_unset_int,
                    num_swb, npft, nlsoil, 0, regen)

        # Cohorts: oldest patch gets 2 (woody + grass), mid + youngest get 1 each.
        cO1 = CLM.fates_cohort_type(height = 25.0, pft = 1, c_area = 1500.0)
        cO2 = CLM.fates_cohort_type(height = 12.0, pft = 2, c_area = 800.0)
        cO1.shorter = cO2; cO2.taller = cO1
        pOld.tallest = cO1; pOld.shortest = cO2

        cM1 = CLM.fates_cohort_type(height = 18.0, pft = 1, c_area = 2000.0)
        pMid.tallest = cM1; pMid.shortest = cM1

        cY1 = CLM.fates_cohort_type(height = 8.0, pft = 2, c_area = 1000.0)
        pYng.tallest = cY1; pYng.shortest = cY1

        # Wire the age-ordered patch list: oldest <-> mid <-> youngest.
        pOld.younger = pMid; pMid.older = pOld
        pMid.younger = pYng; pYng.older = pMid
        site.oldest_patch = pOld
        site.youngest_patch = pYng

        # --- Age ordering traverses oldest -> youngest ----------------------
        @test site.oldest_patch === pOld
        @test site.oldest_patch.younger === pMid
        @test site.oldest_patch.younger.younger === pYng
        @test site.oldest_patch.younger.younger.younger === nothing
        @test site.youngest_patch === pYng
        @test site.youngest_patch.older === pMid
        @test site.youngest_patch.older.older === pOld
        # Walking youngest -> oldest gives ascending age.
        @test site.youngest_patch.age <
              site.youngest_patch.older.age <
              site.youngest_patch.older.older.age

        # Count patches by traversal.
        n = 0; p = site.oldest_patch
        while p !== nothing; n += 1; p = p.younger; end
        @test n == 3

        # --- get_current_landuse_statevector --------------------------------
        # secondaryland area = 3000 + 3000 = 6000; primaryland = 4000.
        sv = CLM.get_current_landuse_statevector(site)
        @test length(sv) == CLM.n_landuse_cats
        @test sv[CLM.secondaryland] ≈ 6000.0 / CLM.area
        @test sv[CLM.primaryland]   ≈ 4000.0 / CLM.area
        @test sv[CLM.cropland] == 0.0

        # --- get_secondary_young_fraction -----------------------------------
        # secondary patches: pOld age 200 (>= threshold 94 -> old),
        # pYng age 5 (< threshold -> young). young/(young+old) = 3000/6000 = 0.5.
        f = CLM.get_secondary_young_fraction(site)
        @test f ≈ 0.5

        # --- CalculateTreeGrassAreaSite! ------------------------------------
        # Each patch's tree/grass area is capped at its own patch area, then
        # summed and normalized by AREA; fractions are sane.
        tf, gf, bf = CLM.CalculateTreeGrassAreaSite!(site)
        @test tf >= 0.0 && gf >= 0.0
        @test gf <= 1.0 - tf + 1e-12
        @test bf ≈ 1.0 - tf - gf
        @test tf + gf + bf ≈ 1.0

        # --- no secondary -> -1 sentinel ------------------------------------
        site2 = CLM.ed_site_type()
        pP = CLM.fates_patch_type()
        CLM.Create!(pP, 10.0, 5000.0, CLM.primaryland, CLM.fates_unset_int,
                    num_swb, npft, nlsoil, 0, regen)
        site2.oldest_patch = pP; site2.youngest_patch = pP
        @test CLM.get_secondary_young_fraction(site2) == -1.0

    finally
        CLM.numpft[] = old_numpft
        CLM.num_elements[] = old_ne
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        CLM.prt_params.woody = old_woody
    end
end
