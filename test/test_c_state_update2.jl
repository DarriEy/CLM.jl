@testset "C State Update 2" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for testing CStateUpdate2
    # ----------------------------------------------------------------
    function make_cstate_update2_data(; nc=4, np=6, ng=2, nlevdecomp=1,
                                       ndecomp_pools=7, ndecomp_cascade_transitions=5,
                                       nrepr=1)
        i_litr_min = 1
        i_litr_max = 3
        i_cwd = 4
        dt = 1800.0  # 30 minutes

        # --- CNVeg carbon state ---
        cs_veg = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs_veg, np, nc, ng; nrepr=nrepr)
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
        cs_veg.gresp_storage_patch     .= 5.0
        cs_veg.gresp_xfer_patch        .= 2.0

        # --- CNVeg carbon flux ---
        cf_veg = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf_veg, np, nc, ng;
                                     nrepr=nrepr,
                                     nlevdecomp_full=nlevdecomp,
                                     ndecomp_pools=ndecomp_pools)

        # Gap mortality fluxes (CStateUpdate2)
        cf_veg.m_leafc_to_litter_patch              .= 1.0e-6
        cf_veg.m_frootc_to_litter_patch             .= 0.8e-6
        cf_veg.m_livestemc_to_litter_patch          .= 0.5e-6
        cf_veg.m_deadstemc_to_litter_patch          .= 0.3e-6
        cf_veg.m_livecrootc_to_litter_patch         .= 0.4e-6
        cf_veg.m_deadcrootc_to_litter_patch         .= 0.2e-6
        cf_veg.m_leafc_storage_to_litter_patch      .= 0.1e-6
        cf_veg.m_frootc_storage_to_litter_patch     .= 0.08e-6
        cf_veg.m_livestemc_storage_to_litter_patch  .= 0.05e-6
        cf_veg.m_deadstemc_storage_to_litter_patch  .= 0.03e-6
        cf_veg.m_livecrootc_storage_to_litter_patch .= 0.04e-6
        cf_veg.m_deadcrootc_storage_to_litter_patch .= 0.02e-6
        cf_veg.m_leafc_xfer_to_litter_patch         .= 0.05e-6
        cf_veg.m_frootc_xfer_to_litter_patch        .= 0.04e-6
        cf_veg.m_livestemc_xfer_to_litter_patch     .= 0.03e-6
        cf_veg.m_deadstemc_xfer_to_litter_patch     .= 0.02e-6
        cf_veg.m_livecrootc_xfer_to_litter_patch    .= 0.015e-6
        cf_veg.m_deadcrootc_xfer_to_litter_patch    .= 0.01e-6
        cf_veg.m_gresp_storage_to_litter_patch      .= 0.01e-6
        cf_veg.m_gresp_xfer_to_litter_patch         .= 0.005e-6
        cf_veg.gap_mortality_c_to_litr_c_col        .= 2.0e-6
        cf_veg.gap_mortality_c_to_cwdc_col          .= 1.0e-6

        # Harvest mortality fluxes (CStateUpdate2h)
        cf_veg.hrv_leafc_to_litter_patch            .= 0.5e-6
        cf_veg.hrv_frootc_to_litter_patch           .= 0.4e-6
        cf_veg.hrv_livestemc_to_litter_patch        .= 0.3e-6
        cf_veg.hrv_livecrootc_to_litter_patch       .= 0.2e-6
        cf_veg.hrv_deadcrootc_to_litter_patch       .= 0.1e-6
        cf_veg.wood_harvestc_patch                  .= 0.6e-6
        cf_veg.hrv_leafc_storage_to_litter_patch    .= 0.05e-6
        cf_veg.hrv_frootc_storage_to_litter_patch   .= 0.04e-6
        cf_veg.hrv_livestemc_storage_to_litter_patch .= 0.03e-6
        cf_veg.hrv_deadstemc_storage_to_litter_patch .= 0.02e-6
        cf_veg.hrv_livecrootc_storage_to_litter_patch .= 0.015e-6
        cf_veg.hrv_deadcrootc_storage_to_litter_patch .= 0.01e-6
        cf_veg.hrv_leafc_xfer_to_litter_patch       .= 0.025e-6
        cf_veg.hrv_frootc_xfer_to_litter_patch      .= 0.02e-6
        cf_veg.hrv_livestemc_xfer_to_litter_patch   .= 0.015e-6
        cf_veg.hrv_deadstemc_xfer_to_litter_patch   .= 0.01e-6
        cf_veg.hrv_livecrootc_xfer_to_litter_patch  .= 0.008e-6
        cf_veg.hrv_deadcrootc_xfer_to_litter_patch  .= 0.005e-6
        cf_veg.hrv_gresp_storage_to_litter_patch    .= 0.005e-6
        cf_veg.hrv_gresp_xfer_to_litter_patch       .= 0.003e-6
        cf_veg.hrv_xsmrpool_to_atm_patch            .= 0.01e-6
        cf_veg.harvest_c_to_litr_c_col              .= 1.5e-6
        cf_veg.harvest_c_to_cwdc_col                .= 0.8e-6

        # Gross unrepresented landcover change fluxes (CStateUpdate2g)
        cf_veg.gru_leafc_to_litter_patch            .= 0.3e-6
        cf_veg.gru_frootc_to_litter_patch           .= 0.25e-6
        cf_veg.gru_livestemc_to_atm_patch           .= 0.2e-6
        cf_veg.gru_deadstemc_to_atm_patch           .= 0.15e-6
        cf_veg.gru_wood_productc_gain_patch         .= 0.1e-6
        cf_veg.gru_livecrootc_to_litter_patch       .= 0.12e-6
        cf_veg.gru_deadcrootc_to_litter_patch       .= 0.08e-6
        cf_veg.gru_leafc_storage_to_atm_patch       .= 0.03e-6
        cf_veg.gru_frootc_storage_to_atm_patch      .= 0.025e-6
        cf_veg.gru_livestemc_storage_to_atm_patch   .= 0.02e-6
        cf_veg.gru_deadstemc_storage_to_atm_patch   .= 0.015e-6
        cf_veg.gru_livecrootc_storage_to_atm_patch  .= 0.012e-6
        cf_veg.gru_deadcrootc_storage_to_atm_patch  .= 0.008e-6
        cf_veg.gru_leafc_xfer_to_atm_patch          .= 0.015e-6
        cf_veg.gru_frootc_xfer_to_atm_patch         .= 0.012e-6
        cf_veg.gru_livestemc_xfer_to_atm_patch      .= 0.01e-6
        cf_veg.gru_deadstemc_xfer_to_atm_patch      .= 0.008e-6
        cf_veg.gru_livecrootc_xfer_to_atm_patch     .= 0.006e-6
        cf_veg.gru_deadcrootc_xfer_to_atm_patch     .= 0.004e-6
        cf_veg.gru_gresp_storage_to_atm_patch       .= 0.003e-6
        cf_veg.gru_gresp_xfer_to_atm_patch          .= 0.002e-6
        cf_veg.gru_xsmrpool_to_atm_patch            .= 0.005e-6
        cf_veg.gru_c_to_litr_c_col                  .= 1.0e-6
        cf_veg.gru_c_to_cwdc_col                    .= 0.5e-6

        # --- Soil biogeochem carbon state ---
        cs_soil = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs_soil, nc, ng, nlevdecomp, ndecomp_pools)
        cs_soil.decomp_cpools_vr_col .= 100.0

        # --- Masks ---
        mask_soilc = trues(nc)
        mask_soilp = trues(np)

        return (cs_veg=cs_veg, cf_veg=cf_veg, cs_soil=cs_soil,
                mask_soilc=mask_soilc, mask_soilp=mask_soilp,
                nc=nc, np=np, ng=ng, dt=dt,
                nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd)
    end

    # ================================================================
    # Test CStateUpdate2 — gap-phase mortality
    # ================================================================
    @testset "c_state_update2! gap mortality" begin
        d = make_cstate_update2_data()

        initial_decomp = copy(d.cs_soil.decomp_cpools_vr_col)
        initial_leafc = copy(d.cs_veg.leafc_patch)
        initial_gresp_storage = copy(d.cs_veg.gresp_storage_patch)
        initial_gresp_xfer = copy(d.cs_veg.gresp_xfer_patch)

        CLM.c_state_update2!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        # Litter pools should increase from gap mortality
        for c in 1:d.nc, j in 1:d.nlevdecomp, i in d.i_litr_min:d.i_litr_max
            expected = initial_decomp[c, j, i] + d.cf_veg.gap_mortality_c_to_litr_c_col[c, j, i] * d.dt
            @test d.cs_soil.decomp_cpools_vr_col[c, j, i] ≈ expected
        end

        # CWD pool should increase
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = initial_decomp[c, j, d.i_cwd] + d.cf_veg.gap_mortality_c_to_cwdc_col[c, j] * d.dt
            @test d.cs_soil.decomp_cpools_vr_col[c, j, d.i_cwd] ≈ expected
        end

        # Patch-level: leafc should decrease from gap mortality
        for p in 1:d.np
            expected = initial_leafc[p] - d.cf_veg.m_leafc_to_litter_patch[p] * d.dt
            @test d.cs_veg.leafc_patch[p] ≈ expected
        end

        # gresp_storage should decrease
        for p in 1:d.np
            expected = initial_gresp_storage[p] - d.cf_veg.m_gresp_storage_to_litter_patch[p] * d.dt
            @test d.cs_veg.gresp_storage_patch[p] ≈ expected
        end

        # gresp_xfer should decrease
        for p in 1:d.np
            expected = initial_gresp_xfer[p] - d.cf_veg.m_gresp_xfer_to_litter_patch[p] * d.dt
            @test d.cs_veg.gresp_xfer_patch[p] ≈ expected
        end
    end

    # ================================================================
    # Test CStateUpdate2 — storage and transfer pools
    # ================================================================
    @testset "c_state_update2! storage and transfer pools" begin
        d = make_cstate_update2_data()

        initial_leafc_storage = copy(d.cs_veg.leafc_storage_patch)
        initial_leafc_xfer = copy(d.cs_veg.leafc_xfer_patch)
        initial_deadcrootc_storage = copy(d.cs_veg.deadcrootc_storage_patch)
        initial_deadcrootc_xfer = copy(d.cs_veg.deadcrootc_xfer_patch)

        CLM.c_state_update2!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        for p in 1:d.np
            @test d.cs_veg.leafc_storage_patch[p] ≈
                  initial_leafc_storage[p] - d.cf_veg.m_leafc_storage_to_litter_patch[p] * d.dt
            @test d.cs_veg.leafc_xfer_patch[p] ≈
                  initial_leafc_xfer[p] - d.cf_veg.m_leafc_xfer_to_litter_patch[p] * d.dt
            @test d.cs_veg.deadcrootc_storage_patch[p] ≈
                  initial_deadcrootc_storage[p] - d.cf_veg.m_deadcrootc_storage_to_litter_patch[p] * d.dt
            @test d.cs_veg.deadcrootc_xfer_patch[p] ≈
                  initial_deadcrootc_xfer[p] - d.cf_veg.m_deadcrootc_xfer_to_litter_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test CStateUpdate2 — mask filtering
    # ================================================================
    @testset "c_state_update2! mask filtering" begin
        d = make_cstate_update2_data()

        # Mask out everything
        d.mask_soilc .= false
        d.mask_soilp .= false

        initial_decomp = copy(d.cs_soil.decomp_cpools_vr_col)
        initial_leafc = copy(d.cs_veg.leafc_patch)

        CLM.c_state_update2!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        @test d.cs_soil.decomp_cpools_vr_col == initial_decomp
        @test d.cs_veg.leafc_patch == initial_leafc
    end

    # ================================================================
    # Test CStateUpdate2h — harvest mortality
    # ================================================================
    @testset "c_state_update2h! harvest mortality" begin
        d = make_cstate_update2_data()

        initial_decomp = copy(d.cs_soil.decomp_cpools_vr_col)
        initial_leafc = copy(d.cs_veg.leafc_patch)
        initial_deadstemc = copy(d.cs_veg.deadstemc_patch)
        initial_xsmrpool = copy(d.cs_veg.xsmrpool_patch)
        initial_gresp_storage = copy(d.cs_veg.gresp_storage_patch)
        initial_gresp_xfer = copy(d.cs_veg.gresp_xfer_patch)

        CLM.c_state_update2h!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        # Litter pools should increase from harvest
        for c in 1:d.nc, j in 1:d.nlevdecomp, i in d.i_litr_min:d.i_litr_max
            expected = initial_decomp[c, j, i] + d.cf_veg.harvest_c_to_litr_c_col[c, j, i] * d.dt
            @test d.cs_soil.decomp_cpools_vr_col[c, j, i] ≈ expected
        end

        # CWD pool should increase
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = initial_decomp[c, j, d.i_cwd] + d.cf_veg.harvest_c_to_cwdc_col[c, j] * d.dt
            @test d.cs_soil.decomp_cpools_vr_col[c, j, d.i_cwd] ≈ expected
        end

        # leafc should decrease
        for p in 1:d.np
            @test d.cs_veg.leafc_patch[p] ≈ initial_leafc[p] - d.cf_veg.hrv_leafc_to_litter_patch[p] * d.dt
        end

        # deadstemc uses wood_harvestc_patch
        for p in 1:d.np
            @test d.cs_veg.deadstemc_patch[p] ≈ initial_deadstemc[p] - d.cf_veg.wood_harvestc_patch[p] * d.dt
        end

        # xsmrpool
        for p in 1:d.np
            @test d.cs_veg.xsmrpool_patch[p] ≈ initial_xsmrpool[p] - d.cf_veg.hrv_xsmrpool_to_atm_patch[p] * d.dt
        end

        # gresp_storage
        for p in 1:d.np
            @test d.cs_veg.gresp_storage_patch[p] ≈
                  initial_gresp_storage[p] - d.cf_veg.hrv_gresp_storage_to_litter_patch[p] * d.dt
        end

        # gresp_xfer
        for p in 1:d.np
            @test d.cs_veg.gresp_xfer_patch[p] ≈
                  initial_gresp_xfer[p] - d.cf_veg.hrv_gresp_xfer_to_litter_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test CStateUpdate2h — harvest storage and transfer pools
    # ================================================================
    @testset "c_state_update2h! harvest storage/transfer" begin
        d = make_cstate_update2_data()

        initial_leafc_storage = copy(d.cs_veg.leafc_storage_patch)
        initial_deadcrootc_xfer = copy(d.cs_veg.deadcrootc_xfer_patch)

        CLM.c_state_update2h!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        for p in 1:d.np
            @test d.cs_veg.leafc_storage_patch[p] ≈
                  initial_leafc_storage[p] - d.cf_veg.hrv_leafc_storage_to_litter_patch[p] * d.dt
            @test d.cs_veg.deadcrootc_xfer_patch[p] ≈
                  initial_deadcrootc_xfer[p] - d.cf_veg.hrv_deadcrootc_xfer_to_litter_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test CStateUpdate2g — gross unrepresented landcover change
    # ================================================================
    @testset "c_state_update2g! gross unrep LCC" begin
        d = make_cstate_update2_data()

        initial_decomp = copy(d.cs_soil.decomp_cpools_vr_col)
        initial_leafc = copy(d.cs_veg.leafc_patch)
        initial_deadstemc = copy(d.cs_veg.deadstemc_patch)
        initial_xsmrpool = copy(d.cs_veg.xsmrpool_patch)
        initial_gresp_storage = copy(d.cs_veg.gresp_storage_patch)
        initial_gresp_xfer = copy(d.cs_veg.gresp_xfer_patch)

        CLM.c_state_update2g!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        # Litter pools should increase
        for c in 1:d.nc, j in 1:d.nlevdecomp, i in d.i_litr_min:d.i_litr_max
            expected = initial_decomp[c, j, i] + d.cf_veg.gru_c_to_litr_c_col[c, j, i] * d.dt
            @test d.cs_soil.decomp_cpools_vr_col[c, j, i] ≈ expected
        end

        # CWD pool should increase
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = initial_decomp[c, j, d.i_cwd] + d.cf_veg.gru_c_to_cwdc_col[c, j] * d.dt
            @test d.cs_soil.decomp_cpools_vr_col[c, j, d.i_cwd] ≈ expected
        end

        # leafc should decrease
        for p in 1:d.np
            @test d.cs_veg.leafc_patch[p] ≈ initial_leafc[p] - d.cf_veg.gru_leafc_to_litter_patch[p] * d.dt
        end

        # deadstemc: two separate decrements (gru_deadstemc_to_atm + gru_wood_productc_gain)
        for p in 1:d.np
            expected = initial_deadstemc[p] -
                       d.cf_veg.gru_deadstemc_to_atm_patch[p] * d.dt -
                       d.cf_veg.gru_wood_productc_gain_patch[p] * d.dt
            @test d.cs_veg.deadstemc_patch[p] ≈ expected
        end

        # xsmrpool
        for p in 1:d.np
            @test d.cs_veg.xsmrpool_patch[p] ≈ initial_xsmrpool[p] - d.cf_veg.gru_xsmrpool_to_atm_patch[p] * d.dt
        end

        # gresp_storage
        for p in 1:d.np
            @test d.cs_veg.gresp_storage_patch[p] ≈
                  initial_gresp_storage[p] - d.cf_veg.gru_gresp_storage_to_atm_patch[p] * d.dt
        end

        # gresp_xfer
        for p in 1:d.np
            @test d.cs_veg.gresp_xfer_patch[p] ≈
                  initial_gresp_xfer[p] - d.cf_veg.gru_gresp_xfer_to_atm_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test CStateUpdate2g — storage and transfer pools
    # ================================================================
    @testset "c_state_update2g! gru storage/transfer" begin
        d = make_cstate_update2_data()

        initial_leafc_storage = copy(d.cs_veg.leafc_storage_patch)
        initial_leafc_xfer = copy(d.cs_veg.leafc_xfer_patch)
        initial_deadcrootc_storage = copy(d.cs_veg.deadcrootc_storage_patch)
        initial_deadcrootc_xfer = copy(d.cs_veg.deadcrootc_xfer_patch)

        CLM.c_state_update2g!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        for p in 1:d.np
            @test d.cs_veg.leafc_storage_patch[p] ≈
                  initial_leafc_storage[p] - d.cf_veg.gru_leafc_storage_to_atm_patch[p] * d.dt
            @test d.cs_veg.leafc_xfer_patch[p] ≈
                  initial_leafc_xfer[p] - d.cf_veg.gru_leafc_xfer_to_atm_patch[p] * d.dt
            @test d.cs_veg.deadcrootc_storage_patch[p] ≈
                  initial_deadcrootc_storage[p] - d.cf_veg.gru_deadcrootc_storage_to_atm_patch[p] * d.dt
            @test d.cs_veg.deadcrootc_xfer_patch[p] ≈
                  initial_deadcrootc_xfer[p] - d.cf_veg.gru_deadcrootc_xfer_to_atm_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test CStateUpdate2g — mask filtering
    # ================================================================
    @testset "c_state_update2g! mask filtering" begin
        d = make_cstate_update2_data()

        d.mask_soilc .= false
        d.mask_soilp .= false

        initial_decomp = copy(d.cs_soil.decomp_cpools_vr_col)
        initial_leafc = copy(d.cs_veg.leafc_patch)

        CLM.c_state_update2g!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        @test d.cs_soil.decomp_cpools_vr_col == initial_decomp
        @test d.cs_veg.leafc_patch == initial_leafc
    end

    # ================================================================
    # Test use_matrixcn skips patch veg updates
    # ================================================================
    @testset "c_state_update2! use_matrixcn skips veg" begin
        d = make_cstate_update2_data()

        initial_leafc = copy(d.cs_veg.leafc_patch)
        initial_gresp_storage = copy(d.cs_veg.gresp_storage_patch)

        CLM.c_state_update2!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            use_matrixcn=true,
            dt=d.dt)

        # With use_matrixcn, veg displayed/storage/xfer pools should NOT change
        @test d.cs_veg.leafc_patch == initial_leafc

        # But gresp_storage/xfer still change (outside the matrixcn guard)
        for p in 1:d.np
            @test d.cs_veg.gresp_storage_patch[p] ≈
                  initial_gresp_storage[p] - d.cf_veg.m_gresp_storage_to_litter_patch[p] * d.dt
        end
    end

end
