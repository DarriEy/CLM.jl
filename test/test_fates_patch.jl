# test_fates_patch.jl
# Tests for FATES Batch 9 (Tier F): FatesPatchMod — the PATCH type. A FATES patch
# is a collection of cohorts sharing a disturbance history/age; it owns the cohort
# linked list, the per-element litter pools, the SPITFIRE fuel object, the
# two-stream radiation object, the canopy/radiation profile arrays, and the
# per-patch running means.
#
# Strategy:
#   * Set the FATES interface dimensions (numpft) and the PARTEH element registry
#     (carbon-only: num_elements = 1, element_list = [carbon12_element]).
#   * Configure prt_params.woody so the tree/grass-area helpers can dispatch.
#   * Create a patch via Create! with a synthetic dimension set, then assert
#     scalar init, fixed-size array shapes, and that the litter/fuel/twostr/
#     running-mean objects were constructed.
#   * Allocate the dynamic profile arrays (ReAllocateDynamics!) and assert shapes.
#   * Attach a 3-cohort tallest/shortest height-ordered linked list and assert the
#     ordering + the tree/grass-area accounting (a CheckPatchArea-style invariant).
#   * Run ZeroValues! and assert the accountable fields are zeroed.

using Test
using CLM

@testset "FATES Batch 9: FatesPatchMod" begin

    # Preserve and restore module-global state mutated by this suite.
    old_numpft   = CLM.numpft[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_woody    = nothing

    try
        npft = 2
        CLM.numpft[]       = npft
        CLM.num_elements[] = 1
        empty!(CLM.element_list)
        push!(CLM.element_list, CLM.carbon12_element)

        # prt_params.woody: PFT1 woody (tree), PFT2 non-woody (grass).
        p = CLM.prt_params
        old_woody = copy(p.woody)
        p.woody = [CLM.itrue, CLM.ifalse]

        num_swb     = CLM.num_swb
        num_levsoil = 4
        current_tod = 0
        regen       = CLM.default_regeneration  # not TRS -> no seedling means

        # ----------------------------------------------------------------------
        # Create a patch.
        # ----------------------------------------------------------------------
        patch = CLM.fates_patch_type()
        CLM.Create!(patch, 5.0, 1000.0, CLM.primaryland, CLM.fates_unset_int,
                    num_swb, npft, num_levsoil, current_tod, regen)

        # --- Scalar field initialization ------------------------------------
        @test patch.age == 5.0
        @test patch.area == 1000.0
        @test patch.age_class == 1
        @test patch.ncl_p == 1
        @test patch.land_use_label == CLM.primaryland
        @test patch.changed_landuse_this_ts == false
        # primaryland -> age_since_anthro_disturbance is the unset sentinel.
        @test patch.age_since_anthro_disturbance == CLM.fates_unset_r8

        # Pointers nulled by Create!/NanValues!.
        @test patch.tallest === nothing
        @test patch.shortest === nothing
        @test patch.older === nothing
        @test patch.younger === nothing

        # --- Fixed-size array shapes ----------------------------------------
        @test length(patch.tr_soil_dir) == num_swb
        @test length(patch.fab) == num_swb
        @test length(patch.fragmentation_scaler) == num_levsoil
        @test size(patch.pft_agb_profile) == (CLM.maxpft, CLM.N_DBH_BINS)
        @test length(patch.canopy_layer_tlai) == CLM.nclmax
        @test size(patch.canopy_mask) == (CLM.nclmax, CLM.maxpft)
        @test length(patch.btran_ft) == CLM.maxpft
        @test length(patch.disturbance_rates) == CLM.N_DIST_TYPES
        @test length(patch.landuse_transition_rates) == CLM.n_landuse_cats
        @test length(patch.scorch_ht) == CLM.maxpft

        # tr_soil_* set to 1.0 in Create!.
        @test all(patch.tr_soil_dir .== 1.0)
        @test all(patch.tr_soil_dif .== 1.0)

        # ZeroValues!-driven fields.
        @test all(patch.canopy_layer_tlai .== 0.0)
        @test patch.total_tree_area == 0.0
        @test patch.total_grass_area == 0.0
        @test patch.zstar == 0.0
        @test all(patch.disturbance_rates .== 0.0)
        @test patch.frac_burnt == 0.0
        @test patch.livegrass == 0.0

        # --- Contained objects constructed ----------------------------------
        # Litter: one per element, with pool arrays sized by pft/soil.
        @test length(patch.litter) == CLM.num_elements[]
        @test patch.litter[1].element_id == CLM.carbon12_element
        @test length(patch.litter[1].seed) == npft
        @test size(patch.litter[1].bg_cwd) == (CLM.ncwd, num_levsoil)

        # Fuel object constructed and zeroed.
        @test patch.fuel !== nothing
        @test patch.fuel.non_trunk_loading == 0.0

        # Two-stream object present but scattering elements not yet allocated.
        @test patch.twostr isa CLM.twostream_type
        @test isempty(patch.twostr.scelg)

        # Running means constructed (non-TRS: only the veg-temperature means).
        @test patch.tveg24 !== nothing
        @test patch.tveg_lpa !== nothing
        @test patch.tveg_longterm !== nothing
        @test isempty(patch.sdlng_mdd)       # TRS-only, not allocated here
        @test isempty(patch.sdlng_emerg_smp)
        # 24-hr mean seeded to the default veg temperature.
        @test patch.tveg24.l_mean ≈ 15.0 + CLM.t_water_freeze_k_1atm

        # --- Dynamic profile arrays: not yet allocated ----------------------
        @test isempty(patch.elai_profile)
        @test isempty(patch.fabd_sun_z)

        # ----------------------------------------------------------------------
        # Allocate the dynamic profile arrays.
        # ----------------------------------------------------------------------
        patch.ncl_p = CLM.nclmax
        fill!(patch.nleaf, 0)
        patch.nleaf[1, 1] = 6   # max leaf layers across (cl, pft)
        CLM.ReAllocateDynamics!(patch)

        nveg_expected = 6 + 1   # ReAllocateDynamics! adds a buffer of 1
        @test size(patch.elai_profile) == (CLM.nclmax, npft, nveg_expected)
        @test size(patch.tlai_profile) == (CLM.nclmax, npft, nveg_expected)
        @test size(patch.f_sun) == (CLM.nclmax, npft, nveg_expected)
        @test size(patch.nrmlzd_parprof_pft_dir_z) ==
              (CLM.num_rad_stream_types, CLM.nclmax, npft, nveg_expected)

        # ZeroDynamics! / NanDynamics! operate on the now-allocated arrays.
        CLM.ZeroDynamics!(patch)
        @test all(patch.elai_profile .== 0.0)
        @test all(patch.f_sun .== 0.0)
        CLM.NanDynamics!(patch)
        @test all(isnan, patch.tlai_profile)

        # ----------------------------------------------------------------------
        # Attach a 3-cohort tallest/shortest height-ordered linked list.
        # tA (tallest, woody) <-> tB (mid, grass) <-> tC (shortest, woody)
        # ----------------------------------------------------------------------
        tA = CLM.fates_cohort_type(height = 30.0, pft = 1, c_area = 100.0)
        tB = CLM.fates_cohort_type(height = 20.0, pft = 2, c_area = 200.0)
        tC = CLM.fates_cohort_type(height = 10.0, pft = 1, c_area = 50.0)

        tA.shorter = tB; tB.taller = tA
        tB.shorter = tC; tC.taller = tB

        patch.tallest  = tA
        patch.shortest = tC

        # Ordering: walk down (tallest -> shorter) gives descending heights.
        @test patch.tallest === tA
        @test patch.tallest.shorter === tB
        @test patch.tallest.shorter.shorter === tC
        @test patch.tallest.shorter.shorter.shorter === nothing
        @test patch.tallest.height >
              patch.tallest.shorter.height >
              patch.tallest.shorter.shorter.height
        # Walk up (shortest -> taller) gives ascending heights.
        @test patch.shortest === tC
        @test patch.shortest.taller === tB
        @test patch.shortest.taller.taller === tA
        @test patch.shortest.height <
              patch.shortest.taller.height <
              patch.shortest.taller.taller.height

        # ----------------------------------------------------------------------
        # Tree/grass area accounting (CheckPatchArea-style invariant):
        #   tree = woody cohorts (PFT1: tA + tC = 150), grass = PFT2 (tB = 200),
        #   each capped at the patch area (1000).
        # ----------------------------------------------------------------------
        CLM.UpdateTreeGrassArea!(patch)
        @test patch.total_tree_area == min(150.0, patch.area)
        @test patch.total_grass_area == min(200.0, patch.area)
        # Invariant: the partitioned areas never exceed the patch area.
        @test patch.total_tree_area <= patch.area
        @test patch.total_grass_area <= patch.area

        # ----------------------------------------------------------------------
        # ZeroValues! resets the accountable scalar/array fields.
        # ----------------------------------------------------------------------
        patch.frac_burnt = 0.7
        patch.livegrass  = 3.0
        fill!(patch.disturbance_rates, 0.5)
        CLM.ZeroValues!(patch)
        @test patch.frac_burnt == 0.0
        @test patch.livegrass == 0.0
        @test all(patch.disturbance_rates .== 0.0)
        @test patch.total_tree_area == 0.0
        @test patch.zstar == 0.0

    finally
        CLM.numpft[]       = old_numpft
        CLM.num_elements[] = old_numel
        empty!(CLM.element_list)
        append!(CLM.element_list, old_ellist)
        if old_woody !== nothing
            CLM.prt_params.woody = old_woody
        end
    end
end
