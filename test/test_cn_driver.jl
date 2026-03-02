@testset "CN Driver" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for testing CN driver
    # Uses set_values! functions to properly zero all flux fields.
    # ----------------------------------------------------------------
    function make_cn_driver_data(; nc=4, np=6, ng=2, nlevdecomp=1,
                                   ndecomp_pools=7,
                                   ndecomp_cascade_transitions=5,
                                   nrepr=1)
        i_litr_min = 1
        i_litr_max = 3
        i_cwd = 4
        dt = 1800.0  # 30 minutes

        # --- Masks ---
        mask_soilc = trues(nc)
        mask_soilp = trues(np)

        # --- CNVeg carbon state ---
        cs_veg = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs_veg, np, nc, ng; nrepr=nrepr)
        CLM.cnveg_carbon_state_set_values!(cs_veg, mask_soilp, 10.0,
                                            mask_soilc, 10.0;
                                            nrepr=nrepr)
        cs_veg.cpool_patch             .= 100.0
        cs_veg.leafc_patch             .= 50.0
        cs_veg.leafc_storage_patch     .= 10.0
        cs_veg.leafc_xfer_patch        .= 5.0
        cs_veg.frootc_patch            .= 40.0
        cs_veg.frootc_storage_patch    .= 8.0
        cs_veg.frootc_xfer_patch       .= 4.0
        cs_veg.livestemc_patch         .= 30.0
        cs_veg.livestemc_storage_patch .= 6.0
        cs_veg.livestemc_xfer_patch    .= 3.0
        cs_veg.deadstemc_patch         .= 200.0
        cs_veg.deadstemc_storage_patch .= 4.0
        cs_veg.deadstemc_xfer_patch    .= 2.0
        cs_veg.livecrootc_patch        .= 20.0
        cs_veg.livecrootc_storage_patch .= 4.0
        cs_veg.livecrootc_xfer_patch   .= 2.0
        cs_veg.deadcrootc_patch        .= 100.0
        cs_veg.deadcrootc_storage_patch .= 2.0
        cs_veg.deadcrootc_xfer_patch   .= 1.0
        cs_veg.xsmrpool_patch          .= 15.0
        cs_veg.xsmrpool_loss_patch     .= 0.0
        cs_veg.gresp_storage_patch     .= 5.0
        cs_veg.gresp_xfer_patch        .= 2.0
        cs_veg.seedc_grc               .= 1000.0
        cs_veg.cropseedc_deficit_patch .= -5.0
        for k in 1:nrepr
            cs_veg.reproductivec_patch[:, k]         .= 10.0
            cs_veg.reproductivec_storage_patch[:, k] .= 2.0
            cs_veg.reproductivec_xfer_patch[:, k]    .= 1.0
        end

        # --- CNVeg carbon flux ---
        cf_veg = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf_veg, np, nc, ng;
                                     nrepr=nrepr,
                                     nlevdecomp_full=nlevdecomp,
                                     ndecomp_pools=ndecomp_pools)
        # Zero all flux fields using the SetValues function
        CLM.cnveg_carbon_flux_set_values!(cf_veg, mask_soilp, 0.0,
                                           mask_soilc, 0.0;
                                           nrepr=nrepr,
                                           nlevdecomp_full=nlevdecomp,
                                           ndecomp_pools=ndecomp_pools)

        # --- CNVeg nitrogen state ---
        ns_veg = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns_veg, np, nc, ng; nrepr=nrepr)
        CLM.cnveg_nitrogen_state_set_values!(ns_veg, mask_soilp, 5.0,
                                              mask_soilc, 5.0;
                                              nrepr=nrepr)
        ns_veg.leafn_patch             .= 5.0
        ns_veg.leafn_storage_patch     .= 1.0
        ns_veg.leafn_xfer_patch        .= 0.5
        ns_veg.frootn_patch            .= 4.0
        ns_veg.frootn_storage_patch    .= 0.8
        ns_veg.frootn_xfer_patch       .= 0.4
        ns_veg.livestemn_patch         .= 3.0
        ns_veg.livestemn_storage_patch .= 0.6
        ns_veg.livestemn_xfer_patch    .= 0.3
        ns_veg.deadstemn_patch         .= 20.0
        ns_veg.deadstemn_storage_patch .= 0.4
        ns_veg.deadstemn_xfer_patch    .= 0.2
        ns_veg.livecrootn_patch        .= 2.0
        ns_veg.livecrootn_storage_patch .= 0.4
        ns_veg.livecrootn_xfer_patch   .= 0.2
        ns_veg.deadcrootn_patch        .= 10.0
        ns_veg.deadcrootn_storage_patch .= 0.2
        ns_veg.deadcrootn_xfer_patch   .= 0.1
        ns_veg.retransn_patch          .= 2.0
        ns_veg.npool_patch             .= 50.0

        # --- CNVeg nitrogen flux ---
        nf_veg = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf_veg, np, nc, ng;
                                       nrepr=nrepr,
                                       nlevdecomp_full=nlevdecomp,
                                       ndecomp_pools=ndecomp_pools,
                                       i_litr_max=i_litr_max)
        # Zero all flux fields
        CLM.cnveg_nitrogen_flux_set_values!(nf_veg, mask_soilp, 0.0,
                                             mask_soilc, 0.0;
                                             nrepr=nrepr,
                                             nlevdecomp_full=nlevdecomp,
                                             ndecomp_pools=ndecomp_pools,
                                             i_litr_max=i_litr_max)

        # --- Soil biogeochem carbon state ---
        cs_soil = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs_soil, nc, ng, nlevdecomp, ndecomp_pools)
        CLM.soil_bgc_carbon_state_set_values!(cs_soil, mask_soilc, 100.0)

        # --- Soil biogeochem carbon flux ---
        cf_soil = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_soil, nc, nlevdecomp, ndecomp_pools,
                                        ndecomp_cascade_transitions)
        CLM.soil_bgc_carbon_flux_set_values!(cf_soil, mask_soilc, 0.0)

        # --- Soil biogeochem nitrogen state ---
        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, nc, ng, nlevdecomp, ndecomp_pools)
        CLM.soil_bgc_nitrogen_state_set_values!(ns_soil, mask_soilc, 5.0)

        # --- Soil biogeochem nitrogen flux ---
        nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf_soil, nc, nlevdecomp, ndecomp_pools,
                                          ndecomp_cascade_transitions)
        CLM.soil_bgc_nitrogen_flux_set_values!(nf_soil, mask_soilc, 0.0)

        # --- Soil biogeochem state ---
        soilbgc_st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(soilbgc_st, nc, np, nlevdecomp,
                                  ndecomp_cascade_transitions)

        # --- Decomposition cascade ---
        cascade_donor_pool    = [1, 2, 3, 1, 2]
        cascade_receiver_pool = [2, 3, 0, 4, 4]  # 0 = terminal

        # --- Patch data ---
        patch_column = [1, 1, 2, 2, 3, 4]
        ivt          = fill(1, np)
        woody        = zeros(Float64, 80)
        woody[2]     = 1.0
        harvdate     = fill(999, np)
        col_is_fates = fill(false, nc)

        # --- Config ---
        config = CLM.CNDriverConfig()

        return (config=config,
                cs_veg=cs_veg, cf_veg=cf_veg,
                ns_veg=ns_veg, nf_veg=nf_veg,
                cs_soil=cs_soil, cf_soil=cf_soil,
                ns_soil=ns_soil, nf_soil=nf_soil,
                soilbgc_st=soilbgc_st,
                cascade_donor_pool=cascade_donor_pool,
                cascade_receiver_pool=cascade_receiver_pool,
                mask_soilc=mask_soilc, mask_soilp=mask_soilp,
                patch_column=patch_column, ivt=ivt, woody=woody,
                harvdate=harvdate, col_is_fates=col_is_fates,
                nc=nc, np=np, ng=ng, dt=dt,
                nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                ndecomp_cascade_transitions=ndecomp_cascade_transitions,
                nrepr=nrepr,
                i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd)
    end

    # ================================================================
    # Test CNDriverConfig construction
    # ================================================================
    @testset "CNDriverConfig defaults" begin
        config = CLM.CNDriverConfig()
        @test config.use_c13 == false
        @test config.use_c14 == false
        @test config.use_cn == false
        @test config.use_fun == false
        @test config.use_nitrif_denitrif == false
        @test config.use_matrixcn == false
        @test config.use_soil_matrixcn == false
        @test config.decomp_method == 1
    end

    # ================================================================
    # Test cn_driver_init! runs without error
    # ================================================================
    @testset "cn_driver_init!" begin
        config = CLM.CNDriverConfig()
        @test CLM.cn_driver_init!(config) === nothing

        config.use_cn = true
        @test CLM.cn_driver_init!(config) === nothing
    end

    # ================================================================
    # Test cn_driver_no_leaching! with zero fluxes — states unchanged
    # ================================================================
    @testset "cn_driver_no_leaching! zero fluxes preserve state" begin
        d = make_cn_driver_data()

        initial_leafc = copy(d.cs_veg.leafc_patch)
        initial_cpool = copy(d.cs_veg.cpool_patch)
        initial_leafn = copy(d.ns_veg.leafn_patch)
        initial_decomp_c = copy(d.cs_soil.decomp_cpools_vr_col)

        CLM.cn_driver_no_leaching!(d.config;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            patch_column=d.patch_column,
            ivt=d.ivt,
            woody=d.woody,
            harvdate=d.harvdate,
            col_is_fates=d.col_is_fates,
            cascade_donor_pool=d.cascade_donor_pool,
            cascade_receiver_pool=d.cascade_receiver_pool,
            dt=d.dt,
            cnveg_cs=d.cs_veg,
            cnveg_cf=d.cf_veg,
            cnveg_ns=d.ns_veg,
            cnveg_nf=d.nf_veg,
            soilbgc_cs=d.cs_soil,
            soilbgc_cf=d.cf_soil,
            soilbgc_ns=d.ns_soil,
            soilbgc_nf=d.nf_soil,
            soilbgc_state=d.soilbgc_st)

        # With zero fluxes, states should be unchanged
        @test d.cs_veg.leafc_patch ≈ initial_leafc
        @test d.cs_veg.cpool_patch ≈ initial_cpool
        @test d.ns_veg.leafn_patch ≈ initial_leafn
        @test d.cs_soil.decomp_cpools_vr_col ≈ initial_decomp_c
    end

    # ================================================================
    # Test cn_driver_no_leaching! zeros flux fields
    # ================================================================
    @testset "cn_driver_no_leaching! zeros flux fields" begin
        d = make_cn_driver_data()

        # Set nonzero values to verify they get zeroed
        d.cf_veg.psnsun_to_cpool_patch  .= 999.0
        d.cf_veg.psnshade_to_cpool_patch .= 999.0
        d.nf_veg.plant_ndemand_patch     .= 999.0
        d.nf_soil.ndep_to_sminn_col      .= 999.0
        d.nf_soil.nfix_to_sminn_col      .= 999.0
        d.nf_soil.fert_to_sminn_col      .= 999.0
        d.nf_soil.soyfixn_to_sminn_col   .= 999.0
        d.nf_soil.ffix_to_sminn_col      .= 999.0

        CLM.cn_driver_no_leaching!(d.config;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            patch_column=d.patch_column,
            ivt=d.ivt,
            woody=d.woody,
            harvdate=d.harvdate,
            col_is_fates=d.col_is_fates,
            cascade_donor_pool=d.cascade_donor_pool,
            cascade_receiver_pool=d.cascade_receiver_pool,
            dt=d.dt,
            cnveg_cs=d.cs_veg,
            cnveg_cf=d.cf_veg,
            cnveg_ns=d.ns_veg,
            cnveg_nf=d.nf_veg,
            soilbgc_cs=d.cs_soil,
            soilbgc_cf=d.cf_soil,
            soilbgc_ns=d.ns_soil,
            soilbgc_nf=d.nf_soil,
            soilbgc_state=d.soilbgc_st)

        # These fields should be zeroed by the driver
        @test all(d.cf_veg.psnsun_to_cpool_patch .== 0.0)
        @test all(d.cf_veg.psnshade_to_cpool_patch .== 0.0)
        @test all(d.nf_veg.plant_ndemand_patch .== 0.0)
        @test all(d.nf_soil.ndep_to_sminn_col .== 0.0)
        @test all(d.nf_soil.nfix_to_sminn_col .== 0.0)
        @test all(d.nf_soil.fert_to_sminn_col .== 0.0)
        @test all(d.nf_soil.soyfixn_to_sminn_col .== 0.0)
        @test all(d.nf_soil.ffix_to_sminn_col .== 0.0)
    end

    # ================================================================
    # Test cn_driver_no_leaching! with nonzero CStateUpdate0 flux
    # ================================================================
    @testset "cn_driver_no_leaching! CStateUpdate0 flux applied" begin
        d = make_cn_driver_data()

        # Set photosynthesis fluxes — CStateUpdate0 adds these to cpool
        psn_flux = 1.0e-6
        d.cf_veg.psnsun_to_cpool_patch  .= psn_flux
        d.cf_veg.psnshade_to_cpool_patch .= psn_flux

        initial_cpool = copy(d.cs_veg.cpool_patch)

        CLM.cn_driver_no_leaching!(d.config;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            patch_column=d.patch_column,
            ivt=d.ivt,
            woody=d.woody,
            harvdate=d.harvdate,
            col_is_fates=d.col_is_fates,
            cascade_donor_pool=d.cascade_donor_pool,
            cascade_receiver_pool=d.cascade_receiver_pool,
            dt=d.dt,
            cnveg_cs=d.cs_veg,
            cnveg_cf=d.cf_veg,
            cnveg_ns=d.ns_veg,
            cnveg_nf=d.nf_veg,
            soilbgc_cs=d.cs_soil,
            soilbgc_cf=d.cf_soil,
            soilbgc_ns=d.ns_soil,
            soilbgc_nf=d.nf_soil,
            soilbgc_state=d.soilbgc_st)

        # NOTE: The driver zeros psnsun/psnshade first, then calls CStateUpdate0
        # which uses the (now-zeroed) psn fluxes. So cpool should be unchanged
        # because the flux was zeroed before CStateUpdate0 applies it.
        # This verifies the correct ordering: zero first, then update.
        @test d.cs_veg.cpool_patch ≈ initial_cpool
    end

    # ================================================================
    # Test cn_driver_no_leaching! skips veg block when no veg patches
    # ================================================================
    @testset "cn_driver_no_leaching! no veg patches" begin
        d = make_cn_driver_data()

        # All veg patches masked off
        mask_no_veg = falses(d.np)

        initial_leafc = copy(d.cs_veg.leafc_patch)
        initial_decomp = copy(d.cs_soil.decomp_cpools_vr_col)

        CLM.cn_driver_no_leaching!(d.config;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=mask_no_veg,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            patch_column=d.patch_column,
            ivt=d.ivt,
            woody=d.woody,
            harvdate=d.harvdate,
            col_is_fates=d.col_is_fates,
            cascade_donor_pool=d.cascade_donor_pool,
            cascade_receiver_pool=d.cascade_receiver_pool,
            dt=d.dt,
            cnveg_cs=d.cs_veg,
            cnveg_cf=d.cf_veg,
            cnveg_ns=d.ns_veg,
            cnveg_nf=d.nf_veg,
            soilbgc_cs=d.cs_soil,
            soilbgc_cf=d.cf_soil,
            soilbgc_ns=d.ns_soil,
            soilbgc_nf=d.nf_soil,
            soilbgc_state=d.soilbgc_st)

        # Veg pools should be unchanged (veg block skipped)
        @test d.cs_veg.leafc_patch == initial_leafc
        # Soil pools should also be unchanged (no fluxes)
        @test d.cs_soil.decomp_cpools_vr_col == initial_decomp
    end

    # ================================================================
    # Test cn_driver_leaching! runs without error
    # ================================================================
    @testset "cn_driver_leaching! basic" begin
        d = make_cn_driver_data()

        CLM.cn_driver_leaching!(d.config;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt,
            soilbgc_ns=d.ns_soil,
            soilbgc_nf=d.nf_soil,
            cnveg_ns=d.ns_veg,
            cnveg_nf=d.nf_veg)

        # With zero fluxes, states should be unchanged
        @test all(d.ns_soil.sminn_vr_col .≈ 5.0)
    end

    # ================================================================
    # Test cn_driver_leaching! with leaching flux
    # ================================================================
    @testset "cn_driver_leaching! applies leaching" begin
        d = make_cn_driver_data()

        # Set nonzero leaching flux
        d.nf_soil.sminn_leached_vr_col .= 1.0e-6

        initial_sminn = copy(d.ns_soil.sminn_vr_col)

        CLM.cn_driver_leaching!(d.config;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt,
            soilbgc_ns=d.ns_soil,
            soilbgc_nf=d.nf_soil,
            cnveg_ns=d.ns_veg,
            cnveg_nf=d.nf_veg)

        # Leaching should reduce sminn
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = initial_sminn[c, j] - 1.0e-6 * d.dt
            @test d.ns_soil.sminn_vr_col[c, j] ≈ expected
        end
    end

    # ================================================================
    # Test cn_driver_summarize_states! runs without error
    # ================================================================
    @testset "cn_driver_summarize_states!" begin
        d = make_cn_driver_data()

        result = CLM.cn_driver_summarize_states!(d.config;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            cnveg_cs=d.cs_veg,
            cnveg_ns=d.ns_veg,
            soilbgc_cs=d.cs_soil,
            soilbgc_ns=d.ns_soil)

        @test result === nothing
    end

    # ================================================================
    # Test cn_driver_summarize_fluxes! runs without error
    # ================================================================
    @testset "cn_driver_summarize_fluxes!" begin
        d = make_cn_driver_data()

        result = CLM.cn_driver_summarize_fluxes!(d.config;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            cnveg_cf=d.cf_veg,
            cnveg_nf=d.nf_veg,
            soilbgc_cf=d.cf_soil,
            soilbgc_nf=d.nf_soil)

        @test result === nothing
    end

    # ================================================================
    # Test decomp method constants
    # ================================================================
    @testset "decomp method constants" begin
        @test CLM.CENTURY_DECOMP == 1
        @test CLM.MIMICS_DECOMP == 2
    end

end
