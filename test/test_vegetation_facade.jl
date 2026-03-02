@testset "Vegetation Facade" begin

    # ----------------------------------------------------------------
    # Helper: create minimal CNVegetationData for testing
    # ----------------------------------------------------------------
    function make_veg_facade_data(; nc=4, np=6, ng=2, nlevdecomp=1,
                                    ndecomp_pools=7,
                                    ndecomp_cascade_transitions=5,
                                    nrepr=1)
        cfg = CLM.CNVegetationConfig(use_cn=true)
        veg = CLM.CNVegetationData(config=cfg)
        CLM.cn_vegetation_init!(veg, np, nc, ng;
                                 nlevdecomp=nlevdecomp,
                                 ndecomp_pools=ndecomp_pools,
                                 ndecomp_cascade_transitions=ndecomp_cascade_transitions,
                                 nrepr=nrepr)

        # Set up state values for testing
        mask_soilc = trues(nc)
        mask_soilp = trues(np)

        CLM.cnveg_carbon_state_set_values!(veg.cnveg_carbonstate_inst,
            mask_soilp, 10.0, mask_soilc, 10.0; nrepr=nrepr)
        veg.cnveg_carbonstate_inst.cpool_patch             .= 100.0
        veg.cnveg_carbonstate_inst.leafc_patch             .= 50.0
        veg.cnveg_carbonstate_inst.frootc_patch            .= 40.0
        veg.cnveg_carbonstate_inst.livecrootc_patch        .= 20.0
        veg.cnveg_carbonstate_inst.totvegc_col             .= 500.0

        CLM.cnveg_carbon_flux_set_values!(veg.cnveg_carbonflux_inst,
            mask_soilp, 0.0, mask_soilc, 0.0;
            nrepr=nrepr, nlevdecomp_full=nlevdecomp,
            ndecomp_pools=ndecomp_pools)
        veg.cnveg_carbonflux_inst.nbp_grc           .= 1.5
        veg.cnveg_carbonflux_inst.rr_patch           .= 0.01
        veg.cnveg_carbonflux_inst.annsum_npp_patch   .= 500.0
        veg.cnveg_carbonflux_inst.agnpp_patch        .= 0.001
        veg.cnveg_carbonflux_inst.bgnpp_patch        .= 0.0005

        CLM.cnveg_nitrogen_state_set_values!(veg.cnveg_nitrogenstate_inst,
            mask_soilp, 5.0, mask_soilc, 5.0; nrepr=nrepr)
        veg.cnveg_nitrogenstate_inst.leafn_patch .= 3.0

        CLM.cnveg_nitrogen_flux_set_values!(veg.cnveg_nitrogenflux_inst,
            mask_soilp, 0.0, mask_soilc, 0.0;
            nrepr=nrepr, nlevdecomp_full=nlevdecomp,
            ndecomp_pools=ndecomp_pools)

        veg.cnveg_state_inst.downreg_patch .= 0.1

        return (veg=veg, nc=nc, np=np, ng=ng,
                nlevdecomp=nlevdecomp,
                ndecomp_pools=ndecomp_pools,
                ndecomp_cascade_transitions=ndecomp_cascade_transitions,
                nrepr=nrepr,
                mask_soilc=mask_soilc, mask_soilp=mask_soilp)
    end

    # ================================================================
    # Test CNVegetationConfig defaults
    # ================================================================
    @testset "CNVegetationConfig defaults" begin
        cfg = CLM.CNVegetationConfig()
        @test cfg.use_cn == false
        @test cfg.use_c13 == false
        @test cfg.use_c14 == false
        @test cfg.use_cndv == false
        @test cfg.use_fates_bgc == false
        @test cfg.reseed_dead_plants == false
        @test cfg.dribble_crophrv_xsmrpool_2atm == false
        @test cfg.skip_steps == 0
    end

    # ================================================================
    # Test CNVegetationData construction
    # ================================================================
    @testset "CNVegetationData construction" begin
        veg = CLM.CNVegetationData()
        @test veg.config isa CLM.CNVegetationConfig
        @test veg.driver_config isa CLM.CNDriverConfig
        @test veg.cnveg_state_inst isa CLM.CNVegStateData
        @test veg.cnveg_carbonstate_inst isa CLM.CNVegCarbonStateData
        @test veg.cnveg_carbonflux_inst isa CLM.CNVegCarbonFluxData
        @test veg.cnveg_nitrogenstate_inst isa CLM.CNVegNitrogenStateData
        @test veg.cnveg_nitrogenflux_inst isa CLM.CNVegNitrogenFluxData
    end

    # ================================================================
    # Test cn_vegetation_init!
    # ================================================================
    @testset "cn_vegetation_init!" begin
        cfg = CLM.CNVegetationConfig(use_cn=true)
        veg = CLM.CNVegetationData(config=cfg)

        @test CLM.cn_vegetation_init!(veg, 6, 4, 2) === nothing

        # Verify sub-types were initialized
        @test length(veg.cnveg_state_inst.downreg_patch) == 6
        @test length(veg.cnveg_carbonstate_inst.leafc_patch) == 6
        @test length(veg.cnveg_nitrogenstate_inst.leafn_patch) == 6
        @test length(veg.cnveg_carbonflux_inst.rr_patch) == 6
        @test length(veg.cnveg_nitrogenflux_inst.plant_ndemand_patch) == 6

        # Driver config should be synced
        @test veg.driver_config.use_cn == true
    end

    @testset "cn_vegetation_init! without CN" begin
        cfg = CLM.CNVegetationConfig(use_cn=false)
        veg = CLM.CNVegetationData(config=cfg)

        @test CLM.cn_vegetation_init!(veg, 6, 4, 2) === nothing

        # State should be initialized (always)
        @test length(veg.cnveg_state_inst.downreg_patch) == 6

        # Carbon/nitrogen types should not be initialized (empty)
        @test length(veg.cnveg_carbonstate_inst.leafc_patch) == 0
    end

    @testset "cn_vegetation_init! with isotopes" begin
        cfg = CLM.CNVegetationConfig(use_cn=true, use_c13=true, use_c14=true)
        veg = CLM.CNVegetationData(config=cfg)

        @test CLM.cn_vegetation_init!(veg, 6, 4, 2) === nothing

        # Isotope types should also be initialized
        @test length(veg.c13_cnveg_carbonstate_inst.leafc_patch) == 6
        @test length(veg.c14_cnveg_carbonstate_inst.leafc_patch) == 6
        @test length(veg.c13_cnveg_carbonflux_inst.rr_patch) == 6
        @test length(veg.c14_cnveg_carbonflux_inst.rr_patch) == 6
    end

    # ================================================================
    # Test cn_vegetation_init2!
    # ================================================================
    @testset "cn_vegetation_init2!" begin
        d = make_veg_facade_data()
        @test CLM.cn_vegetation_init2!(d.veg; bounds=1:d.nc) === nothing
    end

    # ================================================================
    # Test cn_vegetation_init_each_timestep!
    # ================================================================
    @testset "cn_vegetation_init_each_timestep!" begin
        d = make_veg_facade_data()
        @test CLM.cn_vegetation_init_each_timestep!(d.veg;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np) === nothing
    end

    # ================================================================
    # Test getter functions with use_cn=true
    # ================================================================
    @testset "get_net_carbon_exchange_grc" begin
        d = make_veg_facade_data()
        result = CLM.get_net_carbon_exchange_grc(d.veg, 1:d.ng)
        @test length(result) == d.ng
        @test all(result .≈ -1.5)  # -nbp_grc
    end

    @testset "get_leafn_patch" begin
        d = make_veg_facade_data()
        result = CLM.get_leafn_patch(d.veg, 1:d.np)
        @test length(result) == d.np
        @test all(result .≈ 3.0)
    end

    @testset "get_downreg_patch" begin
        d = make_veg_facade_data()
        result = CLM.get_downreg_patch(d.veg, 1:d.np)
        @test length(result) == d.np
        @test all(result .≈ 0.1)
    end

    @testset "get_root_respiration_patch" begin
        d = make_veg_facade_data()
        result = CLM.get_root_respiration_patch(d.veg, 1:d.np)
        @test length(result) == d.np
        @test all(result .≈ 0.01)
    end

    @testset "get_annsum_npp_patch" begin
        d = make_veg_facade_data()
        result = CLM.get_annsum_npp_patch(d.veg, 1:d.np)
        @test length(result) == d.np
        @test all(result .≈ 500.0)
    end

    @testset "get_agnpp_patch" begin
        d = make_veg_facade_data()
        result = CLM.get_agnpp_patch(d.veg, 1:d.np)
        @test length(result) == d.np
        @test all(result .≈ 0.001)
    end

    @testset "get_bgnpp_patch" begin
        d = make_veg_facade_data()
        result = CLM.get_bgnpp_patch(d.veg, 1:d.np)
        @test length(result) == d.np
        @test all(result .≈ 0.0005)
    end

    @testset "get_totvegc_col" begin
        d = make_veg_facade_data()
        result = CLM.get_totvegc_col(d.veg, 1:d.nc)
        @test length(result) == d.nc
        @test all(result .≈ 500.0)
    end

    @testset "get_froot_carbon_patch with CN" begin
        d = make_veg_facade_data()
        result = CLM.get_froot_carbon_patch(d.veg, 1:d.np)
        @test length(result) == d.np
        @test all(result .≈ 40.0)  # frootc_patch
    end

    @testset "get_croot_carbon_patch with CN" begin
        d = make_veg_facade_data()
        result = CLM.get_croot_carbon_patch(d.veg, 1:d.np)
        @test length(result) == d.np
        @test all(result .≈ 20.0)  # livecrootc_patch
    end

    # ================================================================
    # Test getter functions with use_cn=false
    # ================================================================
    @testset "getters with use_cn=false return NaN or zero" begin
        cfg = CLM.CNVegetationConfig(use_cn=false)
        veg = CLM.CNVegetationData(config=cfg)
        CLM.cn_vegetation_init!(veg, 6, 4, 2)

        np = 6; nc = 4; ng = 2

        # Net carbon exchange returns zeros when not CN
        nce = CLM.get_net_carbon_exchange_grc(veg, 1:ng)
        @test all(nce .== 0.0)

        # Others return NaN when not CN
        @test all(isnan.(CLM.get_leafn_patch(veg, 1:np)))
        @test all(isnan.(CLM.get_downreg_patch(veg, 1:np)))
        @test all(isnan.(CLM.get_root_respiration_patch(veg, 1:np)))
        @test all(isnan.(CLM.get_annsum_npp_patch(veg, 1:np)))
        @test all(isnan.(CLM.get_agnpp_patch(veg, 1:np)))
        @test all(isnan.(CLM.get_bgnpp_patch(veg, 1:np)))
        @test all(isnan.(CLM.get_totvegc_col(veg, 1:nc)))
    end

    @testset "get_froot_carbon_patch without CN uses LAI" begin
        cfg = CLM.CNVegetationConfig(use_cn=false)
        veg = CLM.CNVegetationData(config=cfg)
        CLM.cn_vegetation_init!(veg, 3, 2, 1)

        tlai = [2.0, 3.0, 1.0]
        slatop = [0.05, 0.1]        # per PFT
        froot_leaf = [0.6, 0.8]     # per PFT
        ivt = [1, 2, 1]             # PFT index per patch

        result = CLM.get_froot_carbon_patch(veg, 1:3;
            tlai=tlai, slatop=slatop, froot_leaf=froot_leaf, ivt=ivt)

        # patch 1: 2.0/0.05 * 0.6 = 24.0
        @test result[1] ≈ 24.0
        # patch 2: 3.0/0.1 * 0.8 = 24.0
        @test result[2] ≈ 24.0
        # patch 3: 1.0/0.05 * 0.6 = 12.0
        @test result[3] ≈ 12.0
    end

    @testset "get_croot_carbon_patch without CN uses LAI" begin
        cfg = CLM.CNVegetationConfig(use_cn=false)
        veg = CLM.CNVegetationData(config=cfg)
        CLM.cn_vegetation_init!(veg, 3, 2, 1)

        tlai = [2.0, 3.0, 1.0]
        slatop = [0.05, 0.1]
        stem_leaf = [0.5, 0.4]
        croot_stem = [0.3, 0.2]
        ivt = [1, 2, 1]

        result = CLM.get_croot_carbon_patch(veg, 1:3;
            tlai=tlai, slatop=slatop, stem_leaf=stem_leaf,
            croot_stem=croot_stem, ivt=ivt)

        # patch 1: 2.0/0.05 * 0.5 * 0.3 = 6.0
        @test result[1] ≈ 6.0
        # patch 2: 3.0/0.1 * 0.4 * 0.2 = 2.4
        @test result[2] ≈ 2.4
        # patch 3: 1.0/0.05 * 0.5 * 0.3 = 3.0
        @test result[3] ≈ 3.0
    end

    # ================================================================
    # Test _sync_driver_config!
    # ================================================================
    @testset "_sync_driver_config!" begin
        cfg = CLM.CNVegetationConfig(
            use_cn=true, use_c13=true, use_c14=false,
            use_fun=true, use_matrixcn=true)
        veg = CLM.CNVegetationData(config=cfg)
        CLM._sync_driver_config!(veg)

        @test veg.driver_config.use_cn == true
        @test veg.driver_config.use_c13 == true
        @test veg.driver_config.use_c14 == false
        @test veg.driver_config.use_fun == true
        @test veg.driver_config.use_matrixcn == true
    end

    # ================================================================
    # Test cn_vegetation_balance_check! — skips during initial steps
    # ================================================================
    @testset "cn_vegetation_balance_check!" begin
        d = make_veg_facade_data()
        d.veg.config.skip_steps = 5

        # Soil biogeochem data for the call
        cs_soil = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs_soil, d.nc, d.ng, d.nlevdecomp, d.ndecomp_pools)
        cf_soil = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_soil, d.nc, d.nlevdecomp, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)
        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, d.nc, d.ng, d.nlevdecomp, d.ndecomp_pools)
        nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf_soil, d.nc, d.nlevdecomp, d.ndecomp_pools,
                                          d.ndecomp_cascade_transitions)

        # Should skip (nstep <= skip_steps)
        @test CLM.cn_vegetation_balance_check!(d.veg;
            mask_bgc_soilc=d.mask_soilc,
            bounds_col=1:d.nc,
            nstep_since_startup=3,
            soilbgc_cf=cf_soil,
            soilbgc_nf=nf_soil,
            soilbgc_cs=cs_soil,
            soilbgc_ns=ns_soil) === nothing

        # Should not skip (nstep > skip_steps)
        @test CLM.cn_vegetation_balance_check!(d.veg;
            mask_bgc_soilc=d.mask_soilc,
            bounds_col=1:d.nc,
            nstep_since_startup=10,
            soilbgc_cf=cf_soil,
            soilbgc_nf=nf_soil,
            soilbgc_cs=cs_soil,
            soilbgc_ns=ns_soil) === nothing
    end

    # ================================================================
    # Test cn_vegetation_end_of_timestep!
    # ================================================================
    @testset "cn_vegetation_end_of_timestep!" begin
        d = make_veg_facade_data()

        # Without CNDV, should just return
        @test CLM.cn_vegetation_end_of_timestep!(d.veg;
            bounds_patch=1:d.np) === nothing

        # With CNDV but not end of year
        d.veg.config.use_cndv = true
        @test CLM.cn_vegetation_end_of_timestep!(d.veg;
            bounds_patch=1:d.np,
            is_end_curr_year=false) === nothing

        # With CNDV and end of year (not first step)
        @test CLM.cn_vegetation_end_of_timestep!(d.veg;
            bounds_patch=1:d.np,
            is_end_curr_year=true,
            is_first_step=false) === nothing
    end

    # ================================================================
    # Test cn_vegetation_init_column_balance! and gridcell balance
    # ================================================================
    @testset "cn_vegetation_init_column_balance!" begin
        d = make_veg_facade_data()
        cs_soil = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs_soil, d.nc, d.ng, d.nlevdecomp, d.ndecomp_pools)
        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, d.nc, d.ng, d.nlevdecomp, d.ndecomp_pools)

        @test CLM.cn_vegetation_init_column_balance!(d.veg;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            soilbgc_cs=cs_soil,
            soilbgc_ns=ns_soil) === nothing
    end

    @testset "cn_vegetation_init_gridcell_balance!" begin
        d = make_veg_facade_data()
        cs_soil = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs_soil, d.nc, d.ng, d.nlevdecomp, d.ndecomp_pools)
        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, d.nc, d.ng, d.nlevdecomp, d.ndecomp_pools)

        @test CLM.cn_vegetation_init_gridcell_balance!(d.veg;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            soilbgc_cs=cs_soil,
            soilbgc_ns=ns_soil) === nothing
    end

    # ================================================================
    # Test ecosystem dynamics pre-drainage facade
    # ================================================================
    @testset "cn_vegetation_ecosystem_pre_drainage!" begin
        d = make_veg_facade_data()

        # Set up additional required data
        cs_soil = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs_soil, d.nc, d.ng, d.nlevdecomp, d.ndecomp_pools)
        CLM.soil_bgc_carbon_state_set_values!(cs_soil, d.mask_soilc, 100.0)

        cf_soil = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_soil, d.nc, d.nlevdecomp, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)
        CLM.soil_bgc_carbon_flux_set_values!(cf_soil, d.mask_soilc, 0.0)

        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, d.nc, d.ng, d.nlevdecomp, d.ndecomp_pools)
        CLM.soil_bgc_nitrogen_state_set_values!(ns_soil, d.mask_soilc, 5.0)

        nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf_soil, d.nc, d.nlevdecomp, d.ndecomp_pools,
                                          d.ndecomp_cascade_transitions)
        CLM.soil_bgc_nitrogen_flux_set_values!(nf_soil, d.mask_soilc, 0.0)

        soilbgc_st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(soilbgc_st, d.nc, d.np, d.nlevdecomp,
                                  d.ndecomp_cascade_transitions)

        cascade_donor_pool    = [1, 2, 3, 1, 2]
        cascade_receiver_pool = [2, 3, 0, 4, 4]
        patch_column = [1, 1, 2, 2, 3, 4]
        ivt = fill(1, d.np)
        woody = zeros(Float64, 80)
        harvdate = fill(999, d.np)
        col_is_fates = fill(false, d.nc)

        # Also set up veg carbon/nitrogen state properly for the driver
        CLM.cnveg_carbon_state_set_values!(d.veg.cnveg_carbonstate_inst,
            d.mask_soilp, 10.0, d.mask_soilc, 10.0; nrepr=d.nrepr)
        d.veg.cnveg_carbonstate_inst.cpool_patch .= 100.0
        d.veg.cnveg_carbonstate_inst.leafc_patch .= 50.0
        d.veg.cnveg_carbonstate_inst.frootc_patch .= 40.0
        d.veg.cnveg_carbonstate_inst.xsmrpool_patch .= 15.0

        initial_leafc = copy(d.veg.cnveg_carbonstate_inst.leafc_patch)

        @test CLM.cn_vegetation_ecosystem_pre_drainage!(d.veg;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=1,
            i_litr_max=3,
            i_cwd=4,
            patch_column=patch_column,
            ivt=ivt,
            woody=woody,
            harvdate=harvdate,
            col_is_fates=col_is_fates,
            cascade_donor_pool=cascade_donor_pool,
            cascade_receiver_pool=cascade_receiver_pool,
            dt=1800.0,
            soilbgc_cs=cs_soil,
            soilbgc_cf=cf_soil,
            soilbgc_ns=ns_soil,
            soilbgc_nf=nf_soil,
            soilbgc_state=soilbgc_st) === nothing

        # With zero fluxes (zeroed by the driver), states unchanged
        @test d.veg.cnveg_carbonstate_inst.leafc_patch ≈ initial_leafc
    end

    # ================================================================
    # Test ecosystem dynamics post-drainage facade
    # ================================================================
    @testset "cn_vegetation_ecosystem_post_drainage!" begin
        d = make_veg_facade_data()

        cs_soil = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs_soil, d.nc, d.ng, d.nlevdecomp, d.ndecomp_pools)

        cf_soil = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_soil, d.nc, d.nlevdecomp, d.ndecomp_pools,
                                        d.ndecomp_cascade_transitions)

        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, d.nc, d.ng, d.nlevdecomp, d.ndecomp_pools)
        CLM.soil_bgc_nitrogen_state_set_values!(ns_soil, d.mask_soilc, 5.0)

        nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf_soil, d.nc, d.nlevdecomp, d.ndecomp_pools,
                                          d.ndecomp_cascade_transitions)
        CLM.soil_bgc_nitrogen_flux_set_values!(nf_soil, d.mask_soilc, 0.0)

        initial_sminn = copy(ns_soil.sminn_vr_col)

        @test CLM.cn_vegetation_ecosystem_post_drainage!(d.veg;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            dt=1800.0,
            soilbgc_cs=cs_soil,
            soilbgc_cf=cf_soil,
            soilbgc_ns=ns_soil,
            soilbgc_nf=nf_soil) === nothing

        # With zero leaching fluxes, sminn should be unchanged
        @test ns_soil.sminn_vr_col ≈ initial_sminn
    end

end
