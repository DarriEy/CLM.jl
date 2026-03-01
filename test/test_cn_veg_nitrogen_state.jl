@testset "CNVegNitrogenStateData" begin

    @testset "default construction" begin
        ns = CLM.CNVegNitrogenStateData()
        @test length(ns.leafn_patch) == 0
        @test length(ns.frootn_patch) == 0
        @test length(ns.totvegn_col) == 0
        @test length(ns.seedn_grc) == 0
        @test size(ns.reproductiven_patch) == (0, 0)
        @test length(ns.matrix_cap_leafn_patch) == 0
        @test length(ns.leafn0_patch) == 0
    end

    @testset "cnveg_nitrogen_state_init! basic" begin
        np = 10
        nc = 5
        ng = 2
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, np, nc, ng)

        # Patch-level 1D
        @test length(ns.leafn_patch) == np
        @test all(isnan, ns.leafn_patch)
        @test length(ns.leafn_storage_patch) == np
        @test length(ns.frootn_patch) == np
        @test length(ns.livestemn_patch) == np
        @test length(ns.deadstemn_patch) == np
        @test length(ns.livecrootn_patch) == np
        @test length(ns.deadcrootn_patch) == np
        @test length(ns.retransn_patch) == np
        @test length(ns.npool_patch) == np
        @test length(ns.ntrunc_patch) == np
        @test length(ns.cropseedn_deficit_patch) == np
        @test length(ns.leafn_storage_xfer_acc_patch) == np
        @test length(ns.storage_ndemand_patch) == np

        # Summary
        @test length(ns.dispvegn_patch) == np
        @test length(ns.storvegn_patch) == np
        @test length(ns.totvegn_patch) == np
        @test length(ns.totn_patch) == np

        # Patch-level 2D
        @test size(ns.reproductiven_patch) == (np, CLM.NREPR)
        @test all(isnan, ns.reproductiven_patch)

        # Column-level
        @test length(ns.totvegn_col) == nc
        @test all(isnan, ns.totvegn_col)
        @test length(ns.totn_p2c_col) == nc

        # Gridcell-level
        @test length(ns.seedn_grc) == ng
        @test all(isnan, ns.seedn_grc)

        # Matrix CN fields should be empty (use_matrixcn=false)
        @test length(ns.matrix_cap_leafn_patch) == 0
        @test length(ns.leafn0_patch) == 0
        @test length(ns.leafn_SASUsave_patch) == 0
        @test length(ns.matrix_nalloc_leaf_acc_patch) == 0
    end

    @testset "cnveg_nitrogen_state_init! with use_matrixcn" begin
        np = 6
        nc = 3
        ng = 1
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, np, nc, ng; use_matrixcn=true)

        # Matrix capacity fields
        @test length(ns.matrix_cap_leafn_patch) == np
        @test all(isnan, ns.matrix_cap_leafn_patch)
        @test length(ns.matrix_cap_frootn_patch) == np
        @test length(ns.matrix_cap_livestemn_patch) == np
        @test length(ns.matrix_cap_deadstemn_patch) == np
        @test length(ns.matrix_cap_livecrootn_patch) == np
        @test length(ns.matrix_cap_deadcrootn_patch) == np
        @test length(ns.matrix_cap_repron_patch) == np

        # SASU initial pools
        @test length(ns.leafn0_patch) == np
        @test length(ns.frootn0_patch) == np
        @test length(ns.repron0_patch) == np
        @test length(ns.retransn0_patch) == np

        # SASU save fields
        @test length(ns.leafn_SASUsave_patch) == np
        @test length(ns.grainn_SASUsave_patch) == np

        # Matrix accumulation fields
        @test length(ns.matrix_nalloc_leaf_acc_patch) == np
        @test length(ns.matrix_nalloc_grain_acc_patch) == np
        @test length(ns.matrix_ntransfer_leafst_to_leafxf_acc_patch) == np
        @test length(ns.matrix_nturnover_leaf_acc_patch) == np
        @test length(ns.matrix_nturnover_grainxf_acc_patch) == np
        @test length(ns.matrix_nturnover_retransn_acc_patch) == np

        # Retranslocation transfer fields
        @test length(ns.matrix_ntransfer_retransn_to_leaf_acc_patch) == np
        @test length(ns.matrix_ntransfer_leaf_to_retransn_acc_patch) == np
    end

    @testset "cnveg_nitrogen_state_set_values!" begin
        np = 4
        nc = 2
        ng = 1
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, np, nc, ng)

        mask_patch = BitVector([true, false, true, false])
        mask_col = BitVector([true, false])

        CLM.cnveg_nitrogen_state_set_values!(ns, mask_patch, 0.0, mask_col, 0.0)

        # Patches 1,3 should be zero; 2,4 should still be NaN
        @test ns.leafn_patch[1] == 0.0
        @test isnan(ns.leafn_patch[2])
        @test ns.leafn_patch[3] == 0.0
        @test isnan(ns.leafn_patch[4])
        @test ns.npool_patch[1] == 0.0
        @test ns.retransn_patch[3] == 0.0

        # Columns: 1 should be zero, 2 NaN
        @test ns.totvegn_col[1] == 0.0
        @test isnan(ns.totvegn_col[2])
        @test ns.totn_p2c_col[1] == 0.0
    end

    @testset "cnveg_nitrogen_state_set_values! with matrixcn" begin
        np = 3
        nc = 2
        ng = 1
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, np, nc, ng; use_matrixcn=true)

        mask_patch = BitVector([true, true, true])
        mask_col = BitVector([true, true])

        CLM.cnveg_nitrogen_state_set_values!(ns, mask_patch, 5.0, mask_col, 3.0;
                                              use_matrixcn=true)

        @test ns.matrix_cap_leafn_patch[1] == 5.0
        @test ns.leafn0_patch[2] == 5.0
        @test ns.matrix_nalloc_leaf_acc_patch[3] == 5.0
        @test ns.matrix_nturnover_leaf_acc_patch[1] == 5.0
        @test ns.matrix_ntransfer_retransn_to_leaf_acc_patch[2] == 5.0
        @test ns.matrix_ntransfer_leaf_to_retransn_acc_patch[1] == 5.0
        @test ns.totvegn_col[1] == 3.0
    end

    @testset "cnveg_nitrogen_state_zero_dwt!" begin
        np = 5
        nc = 3
        ng = 1
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, np, nc, ng)

        # Set some values first
        ns.dispvegn_patch .= 10.0
        ns.storvegn_patch .= 20.0
        ns.totvegn_patch .= 30.0
        ns.totn_patch .= 40.0

        CLM.cnveg_nitrogen_state_zero_dwt!(ns, 1:np)

        @test all(x -> x == 0.0, ns.dispvegn_patch)
        @test all(x -> x == 0.0, ns.storvegn_patch)
        @test all(x -> x == 0.0, ns.totvegn_patch)
        @test all(x -> x == 0.0, ns.totn_patch)
    end

    @testset "cnveg_nitrogen_state_summary!" begin
        np = 3
        nc = 2
        ng = 1
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, np, nc, ng)

        # Set all pools to known values
        for p in 1:np
            ns.leafn_patch[p]              = 1.0
            ns.frootn_patch[p]             = 2.0
            ns.livestemn_patch[p]          = 3.0
            ns.deadstemn_patch[p]          = 4.0
            ns.livecrootn_patch[p]         = 5.0
            ns.deadcrootn_patch[p]         = 6.0
            ns.leafn_storage_patch[p]      = 0.1
            ns.frootn_storage_patch[p]     = 0.1
            ns.livestemn_storage_patch[p]  = 0.1
            ns.deadstemn_storage_patch[p]  = 0.1
            ns.livecrootn_storage_patch[p] = 0.1
            ns.deadcrootn_storage_patch[p] = 0.1
            ns.leafn_xfer_patch[p]         = 0.1
            ns.frootn_xfer_patch[p]        = 0.1
            ns.livestemn_xfer_patch[p]     = 0.1
            ns.deadstemn_xfer_patch[p]     = 0.1
            ns.livecrootn_xfer_patch[p]    = 0.1
            ns.deadcrootn_xfer_patch[p]    = 0.1
            ns.npool_patch[p]              = 0.5
            ns.retransn_patch[p]           = 0.3
            ns.ntrunc_patch[p]             = 0.2
        end

        mask = BitVector([true, true, true])
        CLM.cnveg_nitrogen_state_summary!(ns, mask, 1:np)

        # dispvegn = 1+2+3+4+5+6 = 21
        @test ns.dispvegn_patch[1] ≈ 21.0

        # storvegn = 6*0.1(storage) + 6*0.1(xfer) + 0.5(npool) + 0.3(retransn) = 2.0
        @test ns.storvegn_patch[1] ≈ 2.0

        # totvegn = dispvegn + storvegn = 21 + 2.0 = 23.0
        @test ns.totvegn_patch[1] ≈ 23.0

        # totn = totvegn + ntrunc = 23.0 + 0.2 = 23.2
        @test ns.totn_patch[1] ≈ 23.2
    end

    @testset "cnveg_nitrogen_state_init_cold! basic" begin
        np = 4
        nc = 2
        ng = 1
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, np, nc, ng)

        CLM.cnveg_nitrogen_state_init_cold!(ns, 1:np)

        for p in 1:np
            @test ns.leafn_patch[p] == 0.0
            @test ns.leafn_storage_patch[p] == 0.0
            @test ns.frootn_patch[p] == 0.0
            @test ns.frootn_storage_patch[p] == 0.0
            @test ns.leafn_xfer_patch[p] == 0.0
            @test ns.livestemn_patch[p] == 0.0
            @test ns.deadstemn_patch[p] == 0.0
            @test ns.livecrootn_patch[p] == 0.0
            @test ns.deadcrootn_patch[p] == 0.0
            @test ns.retransn_patch[p] == 0.0
            @test ns.npool_patch[p] == 0.0
            @test ns.ntrunc_patch[p] == 0.0
            @test ns.dispvegn_patch[p] == 0.0
            @test ns.totvegn_patch[p] == 0.0
            @test ns.totn_patch[p] == 0.0
        end
    end

    @testset "cnveg_nitrogen_state_init_cold! with matrixcn" begin
        np = 3
        nc = 2
        ng = 1
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, np, nc, ng; use_matrixcn=true)

        CLM.cnveg_nitrogen_state_init_cold!(ns, 1:np; use_matrixcn=true)

        for p in 1:np
            @test ns.leafn0_patch[p] ≈ 1.0e-30
            @test ns.frootn0_patch[p] ≈ 1.0e-30
            @test ns.repron0_patch[p] ≈ 1.0e-30
            @test ns.retransn0_patch[p] ≈ 1.0e-30
            @test ns.leafn_SASUsave_patch[p] == 0.0
            @test ns.grainn_SASUsave_patch[p] == 0.0
            @test ns.matrix_nalloc_leaf_acc_patch[p] == 0.0
            @test ns.matrix_nturnover_leaf_acc_patch[p] == 0.0
            @test ns.matrix_ntransfer_retransn_to_leaf_acc_patch[p] == 0.0
            @test ns.matrix_ntransfer_leaf_to_retransn_acc_patch[p] == 0.0
            @test ns.matrix_nturnover_retransn_acc_patch[p] == 0.0
            @test ns.matrix_cap_leafn_patch[p] == 0.0
            @test ns.matrix_cap_frootn_patch[p] == 0.0
        end
    end

    @testset "stub functions run without error" begin
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, 5, 3, 1)

        @test CLM.cnveg_nitrogen_state_init_history!(ns, 1:5, 1:3) === nothing
        @test CLM.cnveg_nitrogen_state_restart!(ns, 1:5, 1:3) === nothing
        @test CLM.cnveg_nitrogen_state_dynamic_patch_adjustments!(ns, 1:5) === nothing
    end

    @testset "field mutability" begin
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, 3, 2, 1)

        ns.leafn_patch[1] = 42.0
        @test ns.leafn_patch[1] == 42.0

        ns.totvegn_col[2] = 99.0
        @test ns.totvegn_col[2] == 99.0

        ns.seedn_grc[1] = 7.5
        @test ns.seedn_grc[1] == 7.5

        ns.reproductiven_patch[1, 1] = 3.14
        @test ns.reproductiven_patch[1, 1] == 3.14
    end

    @testset "re-init overwrites previous state" begin
        ns = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns, 3, 2, 1)
        ns.leafn_patch[1] = 999.0

        CLM.cnveg_nitrogen_state_init!(ns, 7, 4, 2)
        @test length(ns.leafn_patch) == 7
        @test all(isnan, ns.leafn_patch)
        @test length(ns.totvegn_col) == 4
        @test length(ns.seedn_grc) == 2
    end

end
