@testset "CNVegCarbonFlux" begin

    @testset "default construction" begin
        cf = CLM.CNVegCarbonFluxData()
        @test cf.dribble_crophrv_xsmrpool_2atm == false
        @test length(cf.m_leafc_to_litter_patch) == 0
        @test length(cf.gpp_patch) == 0
        @test size(cf.reproductivec_xfer_to_reproductivec_patch) == (0, 0)
        @test size(cf.phenology_c_to_litr_c_col) == (0, 0, 0)
        @test cf.ileafst_to_ileafxf_ph == 0
    end

    @testset "cnveg_carbon_flux_init! basic" begin
        np = 10; nc = 5; ng = 2
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, nc, ng)

        # Patch-level 1D
        @test length(cf.m_leafc_to_litter_patch) == np
        @test all(isnan, cf.m_leafc_to_litter_patch)
        @test length(cf.gpp_patch) == np
        @test all(isnan, cf.gpp_patch)
        @test length(cf.fire_closs_patch) == np
        @test length(cf.npp_Nactive_patch) == np

        # Patch 2D (nrepr)
        @test size(cf.reproductivec_xfer_to_reproductivec_patch) == (np, CLM.NREPR)
        @test size(cf.reproductive_mr_patch) == (np, CLM.NREPR)
        @test size(cf.cpool_to_reproductivec_patch) == (np, CLM.NREPR)

        # Patch 3D (perharv)
        @test size(cf.repr_grainc_to_food_perharv_patch) == (np, CLM.MXHARVESTS, CLM.NREPR)

        # Column-level
        @test length(cf.sr_col) == nc
        @test length(cf.nep_col) == nc
        @test length(cf.cwdc_loss_col) == nc
        @test all(isnan, cf.sr_col)

        # Special init values
        @test all(x -> x == 0.0, cf.hrv_xsmrpool_to_atm_patch)
        @test all(x -> x == 0.0, cf.xsmrpool_to_atm_patch)
        @test all(x -> x == CLM.SPVAL, cf.lag_npp_col)

        # Gridcell-level
        @test length(cf.nbp_grc) == ng
        @test length(cf.nee_grc) == ng

        # 3D column arrays
        @test size(cf.phenology_c_to_litr_c_col) == (nc, 1, 1)
        @test size(cf.m_decomp_cpools_to_fire_vr_col) == (nc, 1, 1)

        # Matrix CN fields should be empty when not enabled
        @test length(cf.matrix_Cinput_patch) == 0
    end

    @testset "cnveg_carbon_flux_init! with use_matrixcn" begin
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, 6, 3, 1; use_matrixcn=true)

        @test length(cf.matrix_Cinput_patch) == 6
        @test all(isnan, cf.matrix_Cinput_patch)
        @test length(cf.matrix_C13input_patch) == 6
        @test length(cf.matrix_C14input_patch) == 6
        @test size(cf.matrix_alloc_patch) == (6, CLM.NVEGPOOL_NATVEG)
    end

    @testset "cnveg_carbon_flux_set_values!" begin
        np = 4; nc = 2; ng = 1
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, nc, ng)

        mask_patch = BitVector([true, false, true, false])
        mask_col   = BitVector([true, false])

        CLM.cnveg_carbon_flux_set_values!(cf, mask_patch, 0.0, mask_col, 0.0)

        # Masked-in patches should be 0
        @test cf.m_leafc_to_litter_patch[1] == 0.0
        @test cf.gpp_patch[1] == 0.0
        @test cf.fire_closs_patch[3] == 0.0
        @test cf.npp_Nactive_patch[1] == 0.0

        # Masked-out patches should remain NaN
        @test isnan(cf.m_leafc_to_litter_patch[2])
        @test isnan(cf.gpp_patch[4])

        # Column-level
        @test cf.sr_col[1] == 0.0
        @test isnan(cf.sr_col[2])
        @test cf.cwdc_loss_col[1] == 0.0

        # Reproductive 2D should remain NaN for non-crop mode
        @test isnan(cf.reproductive_mr_patch[2,1])
        @test cf.reproductive_mr_patch[1,1] == 0.0
    end

    @testset "cnveg_carbon_flux_set_values! with crop" begin
        np = 4; nc = 2; ng = 1
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, nc, ng)

        mask_patch = BitVector([true, true, true, true])
        mask_col   = BitVector([true, true])

        CLM.cnveg_carbon_flux_set_values!(cf, mask_patch, 5.0, mask_col, 3.0;
                                           use_crop=true)

        @test cf.xsmrpool_to_atm_patch[1] == 5.0
        @test cf.livestemc_to_litter_patch[2] == 5.0
        @test cf.reproductivec_xfer_to_reproductivec_patch[1,1] == 5.0
        @test cf.xsmrpool_to_atm_col[1] == 3.0
    end

    @testset "cnveg_carbon_flux_zero_dwt!" begin
        ng = 2; nc = 3
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, 5, nc, ng; nlevdecomp_full=2, ndecomp_pools=2)

        cf.dwt_seedc_to_leaf_grc .= 10.0
        cf.dwt_conv_cflux_grc .= 20.0
        cf.dwt_frootc_to_litr_c_col .= 5.0
        cf.dwt_livecrootc_to_cwdc_col .= 7.0

        CLM.cnveg_carbon_flux_zero_dwt!(cf, 1:ng, 1:nc;
                                         nlevdecomp_full=2, ndecomp_pools=2)

        @test all(x -> x == 0.0, cf.dwt_seedc_to_leaf_grc)
        @test all(x -> x == 0.0, cf.dwt_conv_cflux_grc)
        @test all(x -> x == 0.0, cf.dwt_frootc_to_litr_c_col)
        @test all(x -> x == 0.0, cf.dwt_livecrootc_to_cwdc_col)
    end

    @testset "cnveg_carbon_flux_zero_gru!" begin
        ng = 3
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, 5, 3, ng)

        cf.gru_conv_cflux_grc .= 42.0
        CLM.cnveg_carbon_flux_zero_gru!(cf, 1:ng)

        @test all(x -> x == 0.0, cf.gru_conv_cflux_grc)
    end

    @testset "cnveg_carbon_flux_init_cold!" begin
        np = 4; nc = 2
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, nc, 1; nlevdecomp_full=2, ndecomp_pools=2)

        CLM.cnveg_carbon_flux_init_cold!(cf, 1:np, 1:nc;
                                          nlevdecomp_full=2, ndecomp_pools=2)

        for p in 1:np
            @test cf.gpp_before_downreg_patch[p] == 0.0
            @test cf.gpp_patch[p] == 0.0
            @test cf.availc_patch[p] == 0.0
            @test cf.tempsum_npp_patch[p] == 0.0
            @test cf.annsum_npp_patch[p] == 0.0
        end

        for c in 1:nc
            @test cf.annsum_npp_col[c] == 0.0
            @test cf.dwt_frootc_to_litr_c_col[c,1,1] == 0.0
            @test cf.gru_c_to_litr_c_col[c,2,2] == 0.0
        end
    end

    @testset "cnveg_carbon_flux_summary!" begin
        np = 4; nc = 2; ng = 1
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, nc, ng)

        mask_patch = BitVector([true, true, false, false])

        # Set all flux vars to 0 first via set_values
        mask_all_p = trues(np)
        mask_all_c = trues(nc)
        CLM.cnveg_carbon_flux_set_values!(cf, mask_all_p, 0.0, mask_all_c, 0.0)

        # Set specific values for patch 1
        cf.leaf_mr_patch[1]     = 1.0
        cf.froot_mr_patch[1]    = 2.0
        cf.livestem_mr_patch[1] = 0.5
        cf.livecroot_mr_patch[1] = 0.3

        cf.cpool_leaf_gr_patch[1]     = 0.1
        cf.cpool_froot_gr_patch[1]    = 0.2
        cf.cpool_livestem_gr_patch[1] = 0.05
        cf.cpool_deadstem_gr_patch[1] = 0.04
        cf.cpool_livecroot_gr_patch[1] = 0.03
        cf.cpool_deadcroot_gr_patch[1] = 0.02

        cf.transfer_leaf_gr_patch[1] = 0.01
        cf.transfer_froot_gr_patch[1] = 0.01
        cf.transfer_livestem_gr_patch[1] = 0.01
        cf.transfer_deadstem_gr_patch[1] = 0.01
        cf.transfer_livecroot_gr_patch[1] = 0.01
        cf.transfer_deadcroot_gr_patch[1] = 0.01

        cf.cpool_leaf_storage_gr_patch[1] = 0.005
        cf.cpool_froot_storage_gr_patch[1] = 0.005
        cf.cpool_livestem_storage_gr_patch[1] = 0.005
        cf.cpool_deadstem_storage_gr_patch[1] = 0.005
        cf.cpool_livecroot_storage_gr_patch[1] = 0.005
        cf.cpool_deadcroot_storage_gr_patch[1] = 0.005

        cf.psnsun_to_cpool_patch[1]  = 10.0
        cf.psnshade_to_cpool_patch[1] = 5.0

        CLM.cnveg_carbon_flux_summary!(cf, mask_patch, 1:np)

        # MR = 1.0 + 2.0 + 0.5 + 0.3 = 3.8
        @test cf.mr_patch[1] ≈ 3.8

        # Current GR = 0.1 + 0.2 + 0.05 + 0.04 + 0.03 + 0.02 = 0.44
        @test cf.current_gr_patch[1] ≈ 0.44

        # Transfer GR = 6 * 0.01 = 0.06
        @test cf.transfer_gr_patch[1] ≈ 0.06

        # Storage GR = 6 * 0.005 = 0.03
        @test cf.storage_gr_patch[1] ≈ 0.03

        # GR = 0.44 + 0.06 + 0.03 = 0.53
        @test cf.gr_patch[1] ≈ 0.53

        # AR = MR + GR = 3.8 + 0.53 = 4.33
        @test cf.ar_patch[1] ≈ 4.33

        # GPP = 10 + 5 = 15
        @test cf.gpp_patch[1] ≈ 15.0

        # NPP = GPP - AR = 15 - 4.33 = 10.67
        @test cf.npp_patch[1] ≈ 10.67

        # Patch 2 should be all zeros (was set to 0)
        @test cf.gpp_patch[2] ≈ 0.0
        @test cf.npp_patch[2] ≈ 0.0

        # Masked-out patch 3 should remain 0 (from set_values)
        @test cf.gpp_patch[3] == 0.0
    end

    @testset "stub functions run without error" begin
        cf = CLM.CNVegCarbonFluxData()
        @test CLM.cnveg_carbon_flux_init_history!(cf) === nothing
        @test CLM.cnveg_carbon_flux_restart!(cf) === nothing
    end

    @testset "field mutability" begin
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, 3, 2, 1)

        cf.gpp_patch[1] = 42.0
        @test cf.gpp_patch[1] == 42.0

        cf.sr_col[2] = 99.0
        @test cf.sr_col[2] == 99.0

        cf.reproductive_mr_patch[1, 1] = 3.14
        @test cf.reproductive_mr_patch[1, 1] == 3.14
    end

    @testset "re-init overwrites previous state" begin
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, 3, 2, 1)
        cf.gpp_patch[1] = 999.0

        CLM.cnveg_carbon_flux_init!(cf, 7, 4, 2)
        @test length(cf.gpp_patch) == 7
        @test all(isnan, cf.gpp_patch)
    end

end
