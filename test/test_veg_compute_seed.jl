@testset "VegComputeSeed" begin

    # ----------------------------------------------------------------
    # Helper: create a minimal PftconType with relevant fields populated
    # ----------------------------------------------------------------
    function make_pftcon(; npft=15)
        p = CLM.PftconType()
        n = npft + 1  # Julia 1-based: index 1 = Fortran PFT 0

        p.c3psn     = zeros(n)
        p.evergreen = zeros(n)
        p.woody     = zeros(n)
        p.leafcn    = fill(25.0, n)
        p.deadwdcn  = fill(100.0, n)

        # PFT 0 (noveg): bare ground, index 1
        # PFT 1 (index 2): needleleaf evergreen temperate tree (c3, evergreen, woody)
        p.c3psn[2]     = 1.0
        p.evergreen[2] = 1.0
        p.woody[2]     = 1.0
        p.leafcn[2]    = 30.0
        p.deadwdcn[2]  = 120.0

        # PFT 4 (index 5): broadleaf evergreen tropical tree (c3, evergreen, woody)
        p.c3psn[5]     = 1.0
        p.evergreen[5] = 1.0
        p.woody[5]     = 1.0
        p.leafcn[5]    = 40.0
        p.deadwdcn[5]  = 200.0

        # PFT 7 (index 8): broadleaf deciduous temperate tree (c3, deciduous, woody)
        p.c3psn[8]     = 1.0
        p.evergreen[8] = 0.0
        p.woody[8]     = 1.0
        p.leafcn[8]    = 25.0
        p.deadwdcn[8]  = 100.0

        # PFT 12 (index 13): c3 arctic grass (c3, deciduous, non-woody)
        p.c3psn[13]     = 1.0
        p.evergreen[13] = 0.0
        p.woody[13]     = 0.0
        p.leafcn[13]    = 20.0

        # PFT 14 (index 15): c4 grass (c4, deciduous, non-woody)
        p.c3psn[15]     = 0.0
        p.evergreen[15] = 0.0
        p.woody[15]     = 0.0
        p.leafcn[15]    = 22.0

        return p
    end

    # ----------------------------------------------------------------
    # Tests for leaf_proportions
    # ----------------------------------------------------------------
    @testset "leaf_proportions" begin
        pc = make_pftcon()

        @testset "evergreen PFT with ignore_current_state" begin
            # PFT 1 is evergreen: should put all in leaf
            pleaf, pstor, pxfer = CLM.leaf_proportions(true, 1, 5.0, 3.0, 2.0, pc)
            @test pleaf  == 1.0
            @test pstor  == 0.0
            @test pxfer  == 0.0
        end

        @testset "deciduous PFT with ignore_current_state" begin
            # PFT 7 is deciduous: should put all in storage
            pleaf, pstor, pxfer = CLM.leaf_proportions(true, 7, 5.0, 3.0, 2.0, pc)
            @test pleaf  == 0.0
            @test pstor  == 1.0
            @test pxfer  == 0.0
        end

        @testset "zero total leaf mass uses defaults (evergreen)" begin
            pleaf, pstor, pxfer = CLM.leaf_proportions(false, 1, 0.0, 0.0, 0.0, pc)
            @test pleaf  == 1.0
            @test pstor  == 0.0
            @test pxfer  == 0.0
        end

        @testset "zero total leaf mass uses defaults (deciduous)" begin
            pleaf, pstor, pxfer = CLM.leaf_proportions(false, 7, 0.0, 0.0, 0.0, pc)
            @test pleaf  == 0.0
            @test pstor  == 1.0
            @test pxfer  == 0.0
        end

        @testset "non-zero leaf state, don't ignore" begin
            # PFT 7 (deciduous) with actual leaf state
            pleaf, pstor, pxfer = CLM.leaf_proportions(false, 7, 6.0, 3.0, 1.0, pc)
            tot = 6.0 + 3.0 + 1.0
            @test pleaf  ≈ 6.0 / tot
            @test pstor  ≈ 3.0 / tot
            @test pxfer  ≈ 1.0 / tot
            @test pleaf + pstor + pxfer ≈ 1.0
        end

        @testset "all in storage, non-zero" begin
            pleaf, pstor, pxfer = CLM.leaf_proportions(false, 12, 0.0, 5.0, 0.0, pc)
            @test pleaf  == 0.0
            @test pstor  == 1.0
            @test pxfer  == 0.0
        end

        @testset "all in xfer, non-zero" begin
            pleaf, pstor, pxfer = CLM.leaf_proportions(false, 12, 0.0, 0.0, 4.0, pc)
            @test pleaf  == 0.0
            @test pstor  == 0.0
            @test pxfer  == 1.0
        end
    end

    # ----------------------------------------------------------------
    # Tests for species_type_multiplier
    # ----------------------------------------------------------------
    @testset "species_type_multiplier" begin
        pc = make_pftcon()

        @testset "C12 species always returns 1.0" begin
            @test CLM.species_type_multiplier(CLM.CN_SPECIES_C12, 1, CLM.COMPONENT_LEAF, pc) == 1.0
            @test CLM.species_type_multiplier(CLM.CN_SPECIES_C12, 14, CLM.COMPONENT_LEAF, pc) == 1.0
            @test CLM.species_type_multiplier(CLM.CN_SPECIES_C12, 7, CLM.COMPONENT_DEADWOOD, pc) == 1.0
        end

        @testset "C13 species: C3 plant uses C3_R2" begin
            # PFT 1 is C3 (c3psn = 1.0)
            mult = CLM.species_type_multiplier(CLM.CN_SPECIES_C13, 1, CLM.COMPONENT_LEAF, pc)
            @test mult ≈ CLM.C3_R2
        end

        @testset "C13 species: C4 plant uses C4_R2" begin
            # PFT 14 is C4 (c3psn = 0.0)
            mult = CLM.species_type_multiplier(CLM.CN_SPECIES_C13, 14, CLM.COMPONENT_LEAF, pc)
            @test mult ≈ CLM.C4_R2
        end

        @testset "C14 species returns C14RATIO" begin
            mult = CLM.species_type_multiplier(CLM.CN_SPECIES_C14, 1, CLM.COMPONENT_LEAF, pc)
            @test mult == CLM.C14RATIO
        end

        @testset "N species: leaf component uses 1/leafcn" begin
            # PFT 1: leafcn = 30.0
            mult = CLM.species_type_multiplier(CLM.CN_SPECIES_N, 1, CLM.COMPONENT_LEAF, pc)
            @test mult ≈ 1.0 / 30.0
        end

        @testset "N species: deadwood component uses 1/deadwdcn" begin
            # PFT 7: deadwdcn = 100.0
            mult = CLM.species_type_multiplier(CLM.CN_SPECIES_N, 7, CLM.COMPONENT_DEADWOOD, pc)
            @test mult ≈ 1.0 / 100.0
        end

        @testset "unknown species throws error" begin
            @test_throws ErrorException CLM.species_type_multiplier(99, 1, CLM.COMPONENT_LEAF, pc)
        end

        @testset "N species with unknown component throws error" begin
            @test_throws ErrorException CLM.species_type_multiplier(CLM.CN_SPECIES_N, 1, 99, pc)
        end
    end

    # ----------------------------------------------------------------
    # Tests for compute_seed_amounts!
    # ----------------------------------------------------------------
    @testset "compute_seed_amounts!" begin
        pc = make_pftcon()

        # Helper to build test data for np patches
        function make_seed_data(; np=4)
            patch = CLM.PatchData()
            patch.itype = zeros(Int, np)

            mask   = trues(np)
            bounds = 1:np

            compute_here          = fill(true, np)
            ignore_current_state  = fill(false, np)
            leaf_patch            = zeros(np)
            leaf_storage_patch    = zeros(np)
            leaf_xfer_patch       = zeros(np)
            seed_leaf_patch       = fill(-999.0, np)
            seed_leaf_storage_patch = fill(-999.0, np)
            seed_leaf_xfer_patch  = fill(-999.0, np)
            seed_deadstem_patch   = fill(-999.0, np)

            return (patch=patch, mask=mask, bounds=bounds,
                    compute_here=compute_here,
                    ignore_current_state=ignore_current_state,
                    leaf_patch=leaf_patch,
                    leaf_storage_patch=leaf_storage_patch,
                    leaf_xfer_patch=leaf_xfer_patch,
                    seed_leaf_patch=seed_leaf_patch,
                    seed_leaf_storage_patch=seed_leaf_storage_patch,
                    seed_leaf_xfer_patch=seed_leaf_xfer_patch,
                    seed_deadstem_patch=seed_deadstem_patch)
        end

        @testset "noveg patch: all seeds zero" begin
            d = make_seed_data(np=1)
            d.patch.itype[1] = 0  # noveg

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C12,
                leafc_seed=1.0, deadstemc_seed=0.5,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            @test d.seed_leaf_patch[1]         == 0.0
            @test d.seed_leaf_storage_patch[1] == 0.0
            @test d.seed_leaf_xfer_patch[1]    == 0.0
            @test d.seed_deadstem_patch[1]     == 0.0
        end

        @testset "compute_here false: output not modified" begin
            d = make_seed_data(np=1)
            d.patch.itype[1] = 1
            d.compute_here[1] = false

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C12,
                leafc_seed=1.0, deadstemc_seed=0.5,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            # Sentinel value should be unchanged
            @test d.seed_leaf_patch[1]         == -999.0
            @test d.seed_leaf_storage_patch[1] == -999.0
            @test d.seed_leaf_xfer_patch[1]    == -999.0
            @test d.seed_deadstem_patch[1]     == -999.0
        end

        @testset "mask false: output not modified" begin
            d = make_seed_data(np=1)
            d.patch.itype[1] = 1
            d.mask[1] = false

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C12,
                leafc_seed=1.0, deadstemc_seed=0.5,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            @test d.seed_leaf_patch[1]         == -999.0
            @test d.seed_leaf_storage_patch[1] == -999.0
        end

        @testset "C12 woody evergreen: seed goes to leaf, deadstem nonzero" begin
            # PFT 1: evergreen woody c3 tree
            d = make_seed_data(np=1)
            d.patch.itype[1] = 1
            d.ignore_current_state[1] = true  # force default proportions

            leafc_seed    = 2.0
            deadstemc_seed = 1.0

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C12,
                leafc_seed=leafc_seed, deadstemc_seed=deadstemc_seed,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            # C12 multiplier = 1.0; evergreen -> all in leaf
            @test d.seed_leaf_patch[1]         ≈ leafc_seed * 1.0
            @test d.seed_leaf_storage_patch[1] ≈ 0.0
            @test d.seed_leaf_xfer_patch[1]    ≈ 0.0
            @test d.seed_deadstem_patch[1]     ≈ deadstemc_seed * 1.0
        end

        @testset "C12 woody deciduous: seed goes to storage" begin
            # PFT 7: deciduous woody c3 tree
            d = make_seed_data(np=1)
            d.patch.itype[1] = 7
            d.ignore_current_state[1] = true

            leafc_seed    = 3.0
            deadstemc_seed = 1.5

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C12,
                leafc_seed=leafc_seed, deadstemc_seed=deadstemc_seed,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            # Deciduous -> all in storage
            @test d.seed_leaf_patch[1]         ≈ 0.0
            @test d.seed_leaf_storage_patch[1] ≈ leafc_seed * 1.0
            @test d.seed_leaf_xfer_patch[1]    ≈ 0.0
            @test d.seed_deadstem_patch[1]     ≈ deadstemc_seed * 1.0
        end

        @testset "C12 non-woody grass: no deadstem seed" begin
            # PFT 12: c3 grass, deciduous, non-woody
            d = make_seed_data(np=1)
            d.patch.itype[1] = 12
            d.ignore_current_state[1] = true

            leafc_seed    = 2.0
            deadstemc_seed = 1.0

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C12,
                leafc_seed=leafc_seed, deadstemc_seed=deadstemc_seed,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            # Non-woody: deadstem = 0
            @test d.seed_leaf_patch[1]         ≈ 0.0
            @test d.seed_leaf_storage_patch[1] ≈ leafc_seed * 1.0  # deciduous -> storage
            @test d.seed_leaf_xfer_patch[1]    ≈ 0.0
            @test d.seed_deadstem_patch[1]     ≈ 0.0
        end

        @testset "C13 species: C3 multiplier" begin
            # PFT 1: c3, evergreen, woody
            d = make_seed_data(np=1)
            d.patch.itype[1] = 1
            d.ignore_current_state[1] = true

            leafc_seed    = 2.0
            deadstemc_seed = 1.0

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C13,
                leafc_seed=leafc_seed, deadstemc_seed=deadstemc_seed,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            # C13 + C3 -> multiplier = C3_R2; evergreen -> all in leaf
            @test d.seed_leaf_patch[1]         ≈ leafc_seed * CLM.C3_R2
            @test d.seed_leaf_storage_patch[1] ≈ 0.0
            @test d.seed_deadstem_patch[1]     ≈ deadstemc_seed * CLM.C3_R2
        end

        @testset "C13 species: C4 multiplier" begin
            # PFT 14: c4 grass, deciduous, non-woody
            d = make_seed_data(np=1)
            d.patch.itype[1] = 14
            d.ignore_current_state[1] = true

            leafc_seed    = 2.0
            deadstemc_seed = 1.0

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C13,
                leafc_seed=leafc_seed, deadstemc_seed=deadstemc_seed,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            # C13 + C4 -> multiplier = C4_R2; deciduous -> storage; non-woody -> no deadstem
            @test d.seed_leaf_patch[1]         ≈ 0.0
            @test d.seed_leaf_storage_patch[1] ≈ leafc_seed * CLM.C4_R2
            @test d.seed_deadstem_patch[1]     ≈ 0.0
        end

        @testset "C14 species" begin
            # PFT 4: c3 evergreen woody tree
            d = make_seed_data(np=1)
            d.patch.itype[1] = 4
            d.ignore_current_state[1] = true

            leafc_seed    = 2.0
            deadstemc_seed = 1.0

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C14,
                leafc_seed=leafc_seed, deadstemc_seed=deadstemc_seed,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            @test d.seed_leaf_patch[1]     ≈ leafc_seed * CLM.C14RATIO
            @test d.seed_deadstem_patch[1] ≈ deadstemc_seed * CLM.C14RATIO
        end

        @testset "N species: woody tree" begin
            # PFT 7: deciduous woody tree, leafcn=25, deadwdcn=100
            d = make_seed_data(np=1)
            d.patch.itype[1] = 7
            d.ignore_current_state[1] = true

            leafc_seed    = 2.0
            deadstemc_seed = 1.0

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_N,
                leafc_seed=leafc_seed, deadstemc_seed=deadstemc_seed,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            # N leaf multiplier = 1/leafcn = 1/25; deciduous -> storage
            @test d.seed_leaf_patch[1]         ≈ 0.0
            @test d.seed_leaf_storage_patch[1] ≈ leafc_seed * (1.0 / 25.0)
            @test d.seed_leaf_xfer_patch[1]    ≈ 0.0
            # N deadstem multiplier = 1/deadwdcn = 1/100
            @test d.seed_deadstem_patch[1]     ≈ deadstemc_seed * (1.0 / 100.0)
        end

        @testset "N species: non-woody grass" begin
            # PFT 12: c3 grass, leafcn=20, non-woody -> no deadstem
            d = make_seed_data(np=1)
            d.patch.itype[1] = 12
            d.ignore_current_state[1] = true

            leafc_seed    = 3.0
            deadstemc_seed = 2.0

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_N,
                leafc_seed=leafc_seed, deadstemc_seed=deadstemc_seed,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            @test d.seed_leaf_storage_patch[1] ≈ leafc_seed * (1.0 / 20.0)
            @test d.seed_deadstem_patch[1]     ≈ 0.0
        end

        @testset "current state proportions used when not ignored" begin
            # PFT 7: deciduous woody tree
            d = make_seed_data(np=1)
            d.patch.itype[1] = 7
            d.ignore_current_state[1] = false
            d.leaf_patch[1]         = 6.0
            d.leaf_storage_patch[1] = 3.0
            d.leaf_xfer_patch[1]    = 1.0

            leafc_seed    = 10.0
            deadstemc_seed = 5.0

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C12,
                leafc_seed=leafc_seed, deadstemc_seed=deadstemc_seed,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            tot = 6.0 + 3.0 + 1.0
            @test d.seed_leaf_patch[1]         ≈ leafc_seed * 6.0 / tot
            @test d.seed_leaf_storage_patch[1] ≈ leafc_seed * 3.0 / tot
            @test d.seed_leaf_xfer_patch[1]    ≈ leafc_seed * 1.0 / tot
            @test d.seed_deadstem_patch[1]     ≈ deadstemc_seed
            # Verify total leaf seed
            total_seed = d.seed_leaf_patch[1] + d.seed_leaf_storage_patch[1] + d.seed_leaf_xfer_patch[1]
            @test total_seed ≈ leafc_seed
        end

        @testset "multiple patches with mixed types" begin
            np = 4
            d = make_seed_data(np=np)
            # Patch 1: noveg (type 0)
            # Patch 2: evergreen woody (type 1)
            # Patch 3: deciduous woody (type 7)
            # Patch 4: c4 grass (type 14)
            d.patch.itype .= [0, 1, 7, 14]
            d.ignore_current_state .= true

            leafc_seed    = 2.0
            deadstemc_seed = 1.0

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C12,
                leafc_seed=leafc_seed, deadstemc_seed=deadstemc_seed,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            # Patch 1 (noveg): all zero
            @test d.seed_leaf_patch[1]         == 0.0
            @test d.seed_leaf_storage_patch[1] == 0.0
            @test d.seed_leaf_xfer_patch[1]    == 0.0
            @test d.seed_deadstem_patch[1]     == 0.0

            # Patch 2 (evergreen woody): leaf = seed, deadstem = seed
            @test d.seed_leaf_patch[2]         ≈ leafc_seed
            @test d.seed_leaf_storage_patch[2] ≈ 0.0
            @test d.seed_deadstem_patch[2]     ≈ deadstemc_seed

            # Patch 3 (deciduous woody): storage = seed, deadstem = seed
            @test d.seed_leaf_patch[3]         ≈ 0.0
            @test d.seed_leaf_storage_patch[3] ≈ leafc_seed
            @test d.seed_deadstem_patch[3]     ≈ deadstemc_seed

            # Patch 4 (c4 grass, deciduous non-woody): storage = seed, deadstem = 0
            @test d.seed_leaf_patch[4]         ≈ 0.0
            @test d.seed_leaf_storage_patch[4] ≈ leafc_seed
            @test d.seed_deadstem_patch[4]     ≈ 0.0
        end

        @testset "partial mask: only selected patches modified" begin
            np = 3
            d = make_seed_data(np=np)
            d.patch.itype .= [1, 7, 12]
            d.ignore_current_state .= true

            # Only patches 1 and 3 in mask
            d.mask[2] = false

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C12,
                leafc_seed=2.0, deadstemc_seed=1.0,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            # Patch 1 (mask=true): modified
            @test d.seed_leaf_patch[1] ≈ 2.0  # evergreen -> leaf

            # Patch 2 (mask=false): sentinel preserved
            @test d.seed_leaf_patch[2] == -999.0

            # Patch 3 (mask=true): modified
            @test d.seed_leaf_storage_patch[3] ≈ 2.0  # deciduous grass -> storage
        end

        @testset "seed conservation: leaf + storage + xfer = total leaf seed" begin
            # Using actual current state proportions
            d = make_seed_data(np=1)
            d.patch.itype[1] = 4  # evergreen woody
            d.ignore_current_state[1] = false
            d.leaf_patch[1]         = 3.5
            d.leaf_storage_patch[1] = 1.2
            d.leaf_xfer_patch[1]    = 0.3

            leafc_seed = 5.0

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C12,
                leafc_seed=leafc_seed, deadstemc_seed=0.0,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            total = d.seed_leaf_patch[1] + d.seed_leaf_storage_patch[1] + d.seed_leaf_xfer_patch[1]
            @test total ≈ leafc_seed
        end

        @testset "C13 conservation for C3 tree with proportions" begin
            d = make_seed_data(np=1)
            d.patch.itype[1] = 1  # c3 evergreen woody tree
            d.ignore_current_state[1] = false
            d.leaf_patch[1]         = 4.0
            d.leaf_storage_patch[1] = 2.0
            d.leaf_xfer_patch[1]    = 4.0

            leafc_seed = 6.0

            CLM.compute_seed_amounts!(d.mask, d.bounds, d.patch, pc;
                species=CLM.CN_SPECIES_C13,
                leafc_seed=leafc_seed, deadstemc_seed=0.0,
                leaf_patch=d.leaf_patch,
                leaf_storage_patch=d.leaf_storage_patch,
                leaf_xfer_patch=d.leaf_xfer_patch,
                compute_here_patch=d.compute_here,
                ignore_current_state_patch=d.ignore_current_state,
                seed_leaf_patch=d.seed_leaf_patch,
                seed_leaf_storage_patch=d.seed_leaf_storage_patch,
                seed_leaf_xfer_patch=d.seed_leaf_xfer_patch,
                seed_deadstem_patch=d.seed_deadstem_patch)

            total = d.seed_leaf_patch[1] + d.seed_leaf_storage_patch[1] + d.seed_leaf_xfer_patch[1]
            @test total ≈ leafc_seed * CLM.C3_R2

            # Verify individual proportions
            tot_leaf = 4.0 + 2.0 + 4.0
            @test d.seed_leaf_patch[1]         ≈ leafc_seed * CLM.C3_R2 * 4.0 / tot_leaf
            @test d.seed_leaf_storage_patch[1] ≈ leafc_seed * CLM.C3_R2 * 2.0 / tot_leaf
            @test d.seed_leaf_xfer_patch[1]    ≈ leafc_seed * CLM.C3_R2 * 4.0 / tot_leaf
        end
    end
end
