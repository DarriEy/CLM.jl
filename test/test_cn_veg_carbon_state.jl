@testset "CNVegCarbonStateData" begin

    @testset "default construction" begin
        cs = CLM.CNVegCarbonStateData()
        @test cs.species == 0
        @test length(cs.leafc_patch) == 0
        @test length(cs.frootc_patch) == 0
        @test length(cs.rootc_col) == 0
        @test length(cs.seedc_grc) == 0
        @test size(cs.reproductivec_patch) == (0, 0)
        @test length(cs.matrix_cap_leafc_patch) == 0
        @test length(cs.leafc0_patch) == 0
        @test cs.dribble_crophrv_xsmrpool_2atm == false
    end

    @testset "cnveg_carbon_state_init! basic" begin
        np = 10
        nc = 5
        ng = 2
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, np, nc, ng)

        # Patch-level 1D
        @test length(cs.leafc_patch) == np
        @test all(isnan, cs.leafc_patch)
        @test length(cs.leafc_storage_patch) == np
        @test length(cs.frootc_patch) == np
        @test length(cs.livestemc_patch) == np
        @test length(cs.deadstemc_patch) == np
        @test length(cs.livecrootc_patch) == np
        @test length(cs.deadcrootc_patch) == np
        @test length(cs.gresp_storage_patch) == np
        @test length(cs.gresp_xfer_patch) == np
        @test length(cs.cpool_patch) == np
        @test length(cs.xsmrpool_patch) == np
        @test length(cs.xsmrpool_loss_patch) == np
        @test length(cs.ctrunc_patch) == np
        @test length(cs.woodc_patch) == np
        @test length(cs.leafcmax_patch) == np
        @test length(cs.cropseedc_deficit_patch) == np

        # Summary
        @test length(cs.dispvegc_patch) == np
        @test length(cs.storvegc_patch) == np
        @test length(cs.totc_patch) == np
        @test length(cs.totvegc_patch) == np

        # Patch-level 2D
        @test size(cs.reproductivec_patch) == (np, CLM.NREPR)
        @test all(isnan, cs.reproductivec_patch)

        # Column-level
        @test length(cs.rootc_col) == nc
        @test all(isnan, cs.rootc_col)
        @test length(cs.leafc_col) == nc
        @test length(cs.deadstemc_col) == nc
        @test length(cs.fuelc_col) == nc
        @test length(cs.fuelc_crop_col) == nc
        @test length(cs.totvegc_col) == nc
        @test length(cs.totc_p2c_col) == nc

        # Gridcell-level
        @test length(cs.seedc_grc) == ng
        @test all(isnan, cs.seedc_grc)

        # Matrix CN fields should be empty (use_matrixcn=false)
        @test length(cs.matrix_cap_leafc_patch) == 0
        @test length(cs.leafc0_patch) == 0
        @test length(cs.leafc_SASUsave_patch) == 0
        @test length(cs.matrix_calloc_leaf_acc_patch) == 0
    end

    @testset "cnveg_carbon_state_init! with use_matrixcn" begin
        np = 6
        nc = 3
        ng = 1
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, np, nc, ng; use_matrixcn=true)

        # Matrix capacity fields
        @test length(cs.matrix_cap_leafc_patch) == np
        @test all(isnan, cs.matrix_cap_leafc_patch)
        @test length(cs.matrix_cap_frootc_patch) == np
        @test length(cs.matrix_cap_livestemc_patch) == np
        @test length(cs.matrix_cap_deadstemc_patch) == np
        @test length(cs.matrix_cap_livecrootc_patch) == np
        @test length(cs.matrix_cap_deadcrootc_patch) == np
        @test length(cs.matrix_cap_reproc_patch) == np

        # SASU initial pools
        @test length(cs.leafc0_patch) == np
        @test length(cs.frootc0_patch) == np
        @test length(cs.reproc0_patch) == np

        # SASU save fields
        @test length(cs.leafc_SASUsave_patch) == np
        @test length(cs.grainc_SASUsave_patch) == np

        # Matrix accumulation fields
        @test length(cs.matrix_calloc_leaf_acc_patch) == np
        @test length(cs.matrix_calloc_grain_acc_patch) == np
        @test length(cs.matrix_ctransfer_leafst_to_leafxf_acc_patch) == np
        @test length(cs.matrix_cturnover_leaf_acc_patch) == np
        @test length(cs.matrix_cturnover_grainxf_acc_patch) == np
    end

    @testset "cnveg_carbon_state_set_values!" begin
        np = 4
        nc = 2
        ng = 1
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, np, nc, ng)

        mask_patch = BitVector([true, false, true, false])
        mask_col = BitVector([true, false])

        CLM.cnveg_carbon_state_set_values!(cs, mask_patch, 0.0, mask_col, 0.0)

        # Patches 1,3 should be zero; 2,4 should still be NaN
        @test cs.leafc_patch[1] == 0.0
        @test isnan(cs.leafc_patch[2])
        @test cs.leafc_patch[3] == 0.0
        @test isnan(cs.leafc_patch[4])
        @test cs.cpool_patch[1] == 0.0
        @test cs.woodc_patch[3] == 0.0

        # Columns: 1 should be zero, 2 NaN
        @test cs.rootc_col[1] == 0.0
        @test isnan(cs.rootc_col[2])
        @test cs.leafc_col[1] == 0.0
        @test cs.totvegc_col[1] == 0.0
    end

    @testset "cnveg_carbon_state_set_values! with matrixcn" begin
        np = 3
        nc = 2
        ng = 1
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, np, nc, ng; use_matrixcn=true)

        mask_patch = BitVector([true, true, true])
        mask_col = BitVector([true, true])

        CLM.cnveg_carbon_state_set_values!(cs, mask_patch, 5.0, mask_col, 3.0;
                                            use_matrixcn=true)

        @test cs.matrix_cap_leafc_patch[1] == 5.0
        @test cs.leafc0_patch[2] == 5.0
        @test cs.matrix_calloc_leaf_acc_patch[3] == 5.0
        @test cs.matrix_cturnover_leaf_acc_patch[1] == 5.0
        @test cs.rootc_col[1] == 3.0
    end

    @testset "cnveg_carbon_state_zero_dwt!" begin
        np = 5
        nc = 3
        ng = 1
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, np, nc, ng)

        # Set some values first
        cs.dispvegc_patch .= 10.0
        cs.storvegc_patch .= 20.0
        cs.totc_patch .= 30.0

        CLM.cnveg_carbon_state_zero_dwt!(cs, 1:np)

        @test all(x -> x == 0.0, cs.dispvegc_patch)
        @test all(x -> x == 0.0, cs.storvegc_patch)
        @test all(x -> x == 0.0, cs.totc_patch)
    end

    @testset "cnveg_carbon_state_summary!" begin
        np = 3
        nc = 2
        ng = 1
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, np, nc, ng)

        # Set all pools to known values
        for p in 1:np
            cs.leafc_patch[p]              = 1.0
            cs.frootc_patch[p]             = 2.0
            cs.livestemc_patch[p]          = 3.0
            cs.deadstemc_patch[p]          = 4.0
            cs.livecrootc_patch[p]         = 5.0
            cs.deadcrootc_patch[p]         = 6.0
            cs.cpool_patch[p]              = 0.5
            cs.leafc_storage_patch[p]      = 0.1
            cs.frootc_storage_patch[p]     = 0.1
            cs.livestemc_storage_patch[p]  = 0.1
            cs.deadstemc_storage_patch[p]  = 0.1
            cs.livecrootc_storage_patch[p] = 0.1
            cs.deadcrootc_storage_patch[p] = 0.1
            cs.leafc_xfer_patch[p]         = 0.1
            cs.frootc_xfer_patch[p]        = 0.1
            cs.livestemc_xfer_patch[p]     = 0.1
            cs.deadstemc_xfer_patch[p]     = 0.1
            cs.livecrootc_xfer_patch[p]    = 0.1
            cs.deadcrootc_xfer_patch[p]    = 0.1
            cs.gresp_storage_patch[p]      = 0.1
            cs.gresp_xfer_patch[p]         = 0.1
            cs.xsmrpool_patch[p]           = 0.2
            cs.ctrunc_patch[p]             = 0.3
        end

        mask = BitVector([true, true, true])
        CLM.cnveg_carbon_state_summary!(cs, mask, 1:np)

        # dispvegc = 1+2+3+4+5+6 = 21
        @test cs.dispvegc_patch[1] ≈ 21.0

        # storvegc = 0.5 + 6*0.1 + 6*0.1 + 0.1 + 0.1 = 0.5 + 0.6 + 0.6 + 0.2 = 1.9
        @test cs.storvegc_patch[1] ≈ 1.9

        # totvegc = dispvegc + storvegc = 21 + 1.9 = 22.9
        @test cs.totvegc_patch[1] ≈ 22.9

        # totc = totvegc + xsmrpool + ctrunc = 22.9 + 0.2 + 0.3 = 23.4
        @test cs.totc_patch[1] ≈ 23.4

        # woodc = deadstemc + livestemc + deadcrootc + livecrootc = 4+3+6+5 = 18
        @test cs.woodc_patch[1] ≈ 18.0
    end

    @testset "cnveg_carbon_state_init_cold! basic" begin
        np = 4
        nc = 2
        ng = 1
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, np, nc, ng)

        CLM.cnveg_carbon_state_init_cold!(cs, 1:np)

        for p in 1:np
            @test cs.leafcmax_patch[p] == 0.0
            @test cs.leafc_patch[p] == 0.0
            @test cs.leafc_storage_patch[p] ≈ CLM.INITIAL_VEGC
            @test cs.frootc_patch[p] == 0.0
            @test cs.frootc_storage_patch[p] ≈ CLM.INITIAL_VEGC
            @test cs.leafc_xfer_patch[p] == 0.0
            @test cs.livestemc_patch[p] == 0.0
            @test cs.deadstemc_patch[p] == 0.0
            @test cs.livecrootc_patch[p] == 0.0
            @test cs.deadcrootc_patch[p] == 0.0
            @test cs.gresp_storage_patch[p] == 0.0
            @test cs.cpool_patch[p] == 0.0
            @test cs.xsmrpool_patch[p] == 0.0
            @test cs.ctrunc_patch[p] == 0.0
            @test cs.woodc_patch[p] == 0.0
            @test cs.totc_patch[p] == 0.0
        end
    end

    @testset "cnveg_carbon_state_init_cold! with ratio" begin
        np = 2
        nc = 1
        ng = 1
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, np, nc, ng)

        ratio = 0.5
        CLM.cnveg_carbon_state_init_cold!(cs, 1:np; ratio=ratio)

        @test cs.leafc_storage_patch[1] ≈ CLM.INITIAL_VEGC * ratio
        @test cs.frootc_storage_patch[1] ≈ CLM.INITIAL_VEGC * ratio
    end

    @testset "cnveg_carbon_state_init_cold! with matrixcn" begin
        np = 3
        nc = 2
        ng = 1
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, np, nc, ng; use_matrixcn=true)

        CLM.cnveg_carbon_state_init_cold!(cs, 1:np; use_matrixcn=true)

        for p in 1:np
            @test cs.leafc0_patch[p] ≈ 1.0e-30
            @test cs.frootc0_patch[p] ≈ 1.0e-30
            @test cs.reproc0_patch[p] ≈ 1.0e-30
            @test cs.leafc_SASUsave_patch[p] == 0.0
            @test cs.grainc_SASUsave_patch[p] == 0.0
            @test cs.matrix_calloc_leaf_acc_patch[p] == 0.0
            @test cs.matrix_cturnover_leaf_acc_patch[p] == 0.0
            @test cs.matrix_cap_leafc_storage_patch[p] ≈ CLM.INITIAL_VEGC
            @test cs.matrix_cap_frootc_storage_patch[p] ≈ CLM.INITIAL_VEGC
        end
    end

    @testset "stub functions run without error" begin
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, 5, 3, 1)

        @test CLM.cnveg_carbon_state_init_history!(cs, 1:5, 1:3) === nothing
        @test CLM.cnveg_carbon_state_restart!(cs, 1:5, 1:3) === nothing
        @test CLM.cnveg_carbon_state_dynamic_patch_adjustments!(cs, 1:5) === nothing
    end

    @testset "field mutability" begin
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, 3, 2, 1)

        cs.leafc_patch[1] = 42.0
        @test cs.leafc_patch[1] == 42.0

        cs.rootc_col[2] = 99.0
        @test cs.rootc_col[2] == 99.0

        cs.seedc_grc[1] = 7.5
        @test cs.seedc_grc[1] == 7.5

        cs.reproductivec_patch[1, 1] = 3.14
        @test cs.reproductivec_patch[1, 1] == 3.14
    end

    @testset "re-init overwrites previous state" begin
        cs = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs, 3, 2, 1)
        cs.leafc_patch[1] = 999.0

        CLM.cnveg_carbon_state_init!(cs, 7, 4, 2)
        @test length(cs.leafc_patch) == 7
        @test all(isnan, cs.leafc_patch)
        @test length(cs.rootc_col) == 4
        @test length(cs.seedc_grc) == 2
    end

    @testset "module-level constants" begin
        @test CLM.SPINUP_FACTOR_DEADWOOD_DEFAULT == 1.0
        @test CLM.SPINUP_FACTOR_AD == 10.0
        @test CLM.INITIAL_VEGC == 20.0
    end

end
