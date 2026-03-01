@testset "CNVegNitrogenFlux" begin

    @testset "default construction" begin
        nf = CLM.CNVegNitrogenFluxData()
        @test length(nf.m_leafn_to_litter_patch) == 0
        @test length(nf.fire_nloss_patch) == 0
        @test size(nf.reproductiven_xfer_to_reproductiven_patch) == (0, 0)
        @test size(nf.phenology_n_to_litr_n_col) == (0, 0, 0)
        @test nf.ileafst_to_ileafxf_ph == 0
        @test nf.ileaf_to_iout_gm == 0
    end

    @testset "cnveg_nitrogen_flux_init! basic" begin
        np = 10; nc = 5; ng = 2
        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, np, nc, ng)

        # Patch-level 1D
        @test length(nf.m_leafn_to_litter_patch) == np
        @test all(isnan, nf.m_leafn_to_litter_patch)
        @test length(nf.fire_nloss_patch) == np
        @test all(isnan, nf.fire_nloss_patch)
        @test length(nf.ndeploy_patch) == np
        @test length(nf.plant_ndemand_patch) == np

        # Patch 2D (nrepr)
        @test size(nf.reproductiven_xfer_to_reproductiven_patch) == (np, CLM.NREPR)
        @test size(nf.npool_to_reproductiven_patch) == (np, CLM.NREPR)

        # Patch 3D (perharv)
        @test size(nf.repr_grainn_to_food_perharv_patch) == (np, CLM.MXHARVESTS, CLM.NREPR)

        # Column-level
        @test length(nf.fire_nloss_col) == nc
        @test length(nf.crop_harvestn_to_cropprodn_col) == nc
        @test all(isnan, nf.fire_nloss_col)

        # Gridcell-level
        @test length(nf.dwt_seedn_to_leaf_grc) == ng
        @test length(nf.dwt_conv_nflux_grc) == ng

        # 3D column arrays
        @test size(nf.phenology_n_to_litr_n_col) == (nc, 1, 1)
        @test size(nf.m_decomp_npools_to_fire_vr_col) == (nc, 1, 1)

        # Matrix CN fields should be empty when not enabled
        @test length(nf.matrix_Ninput_patch) == 0
    end

    @testset "cnveg_nitrogen_flux_init! with use_matrixcn" begin
        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, 6, 3, 1; use_matrixcn=true)

        @test length(nf.matrix_Ninput_patch) == 6
        @test all(isnan, nf.matrix_Ninput_patch)
        @test size(nf.matrix_nalloc_patch) == (6, CLM.NVEGPOOL_NATVEG)
    end

    @testset "cnveg_nitrogen_flux_set_values!" begin
        np = 4; nc = 2; ng = 1
        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, np, nc, ng)

        mask_patch = BitVector([true, false, true, false])
        mask_col   = BitVector([true, false])

        CLM.cnveg_nitrogen_flux_set_values!(nf, mask_patch, 0.0, mask_col, 0.0)

        # Masked-in patches should be 0
        @test nf.m_leafn_to_litter_patch[1] == 0.0
        @test nf.ndeploy_patch[1] == 0.0
        @test nf.fire_nloss_patch[3] == 0.0

        # Masked-out patches should remain NaN
        @test isnan(nf.m_leafn_to_litter_patch[2])
        @test isnan(nf.ndeploy_patch[4])

        # Column-level
        @test nf.fire_nloss_col[1] == 0.0
        @test isnan(nf.fire_nloss_col[2])
    end

    @testset "cnveg_nitrogen_flux_set_values! with crop" begin
        np = 4; nc = 2; ng = 1
        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, np, nc, ng)

        mask_patch = BitVector([true, true, true, true])
        mask_col   = BitVector([true, true])

        CLM.cnveg_nitrogen_flux_set_values!(nf, mask_patch, 5.0, mask_col, 3.0;
                                             use_crop=true)

        @test nf.livestemn_to_litter_patch[1] == 5.0
        @test nf.reproductiven_xfer_to_reproductiven_patch[1,1] == 5.0
        @test nf.soyfixn_patch[2] == 5.0
        @test nf.fire_nloss_col[1] == 3.0
    end

    @testset "cnveg_nitrogen_flux_zero_dwt!" begin
        ng = 2; nc = 3
        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, 5, nc, ng; nlevdecomp_full=2)

        nf.dwt_seedn_to_leaf_grc .= 10.0
        nf.dwt_conv_nflux_grc .= 20.0
        nf.dwt_frootn_to_litr_n_col .= 5.0
        nf.dwt_livecrootn_to_cwdn_col .= 7.0

        CLM.cnveg_nitrogen_flux_zero_dwt!(nf, 1:ng, 1:nc; nlevdecomp_full=2)

        @test all(x -> x == 0.0, nf.dwt_seedn_to_leaf_grc)
        @test all(x -> x == 0.0, nf.dwt_conv_nflux_grc)
        @test all(x -> x == 0.0, nf.dwt_frootn_to_litr_n_col)
        @test all(x -> x == 0.0, nf.dwt_livecrootn_to_cwdn_col)
    end

    @testset "cnveg_nitrogen_flux_zero_gru!" begin
        ng = 3
        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, 5, 3, ng)

        nf.gru_conv_nflux_grc .= 42.0
        nf.gru_wood_productn_gain_grc .= 7.0
        CLM.cnveg_nitrogen_flux_zero_gru!(nf, 1:ng)

        @test all(x -> x == 0.0, nf.gru_conv_nflux_grc)
        @test all(x -> x == 0.0, nf.gru_wood_productn_gain_grc)
    end

    @testset "cnveg_nitrogen_flux_init_cold!" begin
        np = 4; nc = 2
        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, np, nc, 1; nlevdecomp_full=2, ndecomp_pools=2)

        CLM.cnveg_nitrogen_flux_init_cold!(nf, 1:np, 1:nc; nlevdecomp_full=2, ndecomp_pools=2)

        for p in 1:np
            @test nf.fert_counter_patch[p] == 0.0
            @test nf.fert_patch[p] == 0.0
            @test nf.soyfixn_patch[p] == 0.0
        end

        for c in 1:nc
            @test nf.dwt_frootn_to_litr_n_col[c,1,1] == 0.0
            @test nf.gru_n_to_litr_n_col[c,2,1] == 0.0
        end
    end

    @testset "cnveg_nitrogen_flux_summary!" begin
        np = 4; nc = 2; ng = 1
        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, np, nc, ng)

        mask_patch = BitVector([true, true, false, false])

        # Set all flux vars to 0
        mask_all_p = trues(np)
        mask_all_c = trues(nc)
        CLM.cnveg_nitrogen_flux_set_values!(nf, mask_all_p, 0.0, mask_all_c, 0.0)

        # Set specific values for patch 1
        nf.sminn_to_npool_patch[1]      = 1.0
        nf.retransn_to_npool_patch[1]   = 2.0
        nf.free_retransn_to_npool_patch[1] = 0.5

        nf.m_leafn_to_fire_patch[1]     = 0.1
        nf.m_frootn_to_fire_patch[1]    = 0.2

        CLM.cnveg_nitrogen_flux_summary!(nf, mask_patch, 1:np)

        # NDEPLOY = 1.0 + 2.0 + 0.5 = 3.5
        @test nf.ndeploy_patch[1] ≈ 3.5

        # Fire N loss should include the two fire fluxes we set
        @test nf.fire_nloss_patch[1] ≈ 0.3

        # Patch 2 should be all zeros
        @test nf.ndeploy_patch[2] ≈ 0.0

        # Masked-out patch 3 should remain 0 (from set_values)
        @test nf.ndeploy_patch[3] == 0.0
    end

    @testset "stub functions run without error" begin
        nf = CLM.CNVegNitrogenFluxData()
        @test CLM.cnveg_nitrogen_flux_init_history!(nf) === nothing
        @test CLM.cnveg_nitrogen_flux_restart!(nf) === nothing
    end

    @testset "field mutability" begin
        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, 3, 2, 1)

        nf.m_leafn_to_litter_patch[1] = 42.0
        @test nf.m_leafn_to_litter_patch[1] == 42.0

        nf.fire_nloss_col[2] = 99.0
        @test nf.fire_nloss_col[2] == 99.0

        nf.reproductiven_xfer_to_reproductiven_patch[1, 1] = 3.14
        @test nf.reproductiven_xfer_to_reproductiven_patch[1, 1] == 3.14
    end

    @testset "re-init overwrites previous state" begin
        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, 3, 2, 1)
        nf.m_leafn_to_litter_patch[1] = 999.0

        CLM.cnveg_nitrogen_flux_init!(nf, 7, 4, 2)
        @test length(nf.m_leafn_to_litter_patch) == 7
        @test all(isnan, nf.m_leafn_to_litter_patch)
    end

end
