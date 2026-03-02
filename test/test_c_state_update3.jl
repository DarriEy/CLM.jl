@testset "C State Update 3" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for testing CStateUpdate3
    # ----------------------------------------------------------------
    function make_cstate_update3_data(; nc=4, np=6, ng=2, nlevdecomp=1,
                                       ndecomp_pools=7, nrepr=1)
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
        cs_veg.gresp_storage_patch     .= 5.0
        cs_veg.gresp_xfer_patch        .= 2.0

        # --- CNVeg carbon flux ---
        cf_veg = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf_veg, np, nc, ng;
                                     nrepr=nrepr,
                                     nlevdecomp_full=nlevdecomp,
                                     ndecomp_pools=ndecomp_pools)

        # Fire mortality fluxes — patch level to combustion
        cf_veg.m_leafc_to_fire_patch              .= 1.0e-6
        cf_veg.m_frootc_to_fire_patch             .= 0.8e-6
        cf_veg.m_livestemc_to_fire_patch          .= 0.5e-6
        cf_veg.m_deadstemc_to_fire_patch          .= 0.3e-6
        cf_veg.m_livecrootc_to_fire_patch         .= 0.4e-6
        cf_veg.m_deadcrootc_to_fire_patch         .= 0.2e-6
        cf_veg.m_leafc_storage_to_fire_patch      .= 0.1e-6
        cf_veg.m_frootc_storage_to_fire_patch     .= 0.08e-6
        cf_veg.m_livestemc_storage_to_fire_patch  .= 0.05e-6
        cf_veg.m_deadstemc_storage_to_fire_patch  .= 0.03e-6
        cf_veg.m_livecrootc_storage_to_fire_patch .= 0.04e-6
        cf_veg.m_deadcrootc_storage_to_fire_patch .= 0.02e-6
        cf_veg.m_leafc_xfer_to_fire_patch         .= 0.05e-6
        cf_veg.m_frootc_xfer_to_fire_patch        .= 0.04e-6
        cf_veg.m_livestemc_xfer_to_fire_patch     .= 0.03e-6
        cf_veg.m_deadstemc_xfer_to_fire_patch     .= 0.02e-6
        cf_veg.m_livecrootc_xfer_to_fire_patch    .= 0.015e-6
        cf_veg.m_deadcrootc_xfer_to_fire_patch    .= 0.01e-6
        cf_veg.m_gresp_storage_to_fire_patch      .= 0.01e-6
        cf_veg.m_gresp_xfer_to_fire_patch         .= 0.005e-6

        # Fire mortality fluxes — patch level to litter
        cf_veg.m_leafc_to_litter_fire_patch              .= 0.5e-6
        cf_veg.m_frootc_to_litter_fire_patch             .= 0.4e-6
        cf_veg.m_livestemc_to_litter_fire_patch          .= 0.3e-6
        cf_veg.m_deadstemc_to_litter_fire_patch          .= 0.15e-6
        cf_veg.m_livecrootc_to_litter_fire_patch         .= 0.2e-6
        cf_veg.m_deadcrootc_to_litter_fire_patch         .= 0.1e-6
        cf_veg.m_leafc_storage_to_litter_fire_patch      .= 0.05e-6
        cf_veg.m_frootc_storage_to_litter_fire_patch     .= 0.04e-6
        cf_veg.m_livestemc_storage_to_litter_fire_patch  .= 0.03e-6
        cf_veg.m_deadstemc_storage_to_litter_fire_patch  .= 0.015e-6
        cf_veg.m_livecrootc_storage_to_litter_fire_patch .= 0.02e-6
        cf_veg.m_deadcrootc_storage_to_litter_fire_patch .= 0.01e-6
        cf_veg.m_leafc_xfer_to_litter_fire_patch         .= 0.025e-6
        cf_veg.m_frootc_xfer_to_litter_fire_patch        .= 0.02e-6
        cf_veg.m_livestemc_xfer_to_litter_fire_patch     .= 0.015e-6
        cf_veg.m_deadstemc_xfer_to_litter_fire_patch     .= 0.01e-6
        cf_veg.m_livecrootc_xfer_to_litter_fire_patch    .= 0.008e-6
        cf_veg.m_deadcrootc_xfer_to_litter_fire_patch    .= 0.005e-6
        cf_veg.m_gresp_storage_to_litter_fire_patch      .= 0.005e-6
        cf_veg.m_gresp_xfer_to_litter_fire_patch         .= 0.003e-6

        # Live-to-dead fire transfers
        cf_veg.m_livestemc_to_deadstemc_fire_patch       .= 0.2e-6
        cf_veg.m_livecrootc_to_deadcrootc_fire_patch     .= 0.15e-6

        # Column-level fire fluxes
        cf_veg.fire_mortality_c_to_cwdc_col   .= 1.0e-6
        cf_veg.m_c_to_litr_fire_col           .= 2.0e-6
        cf_veg.m_decomp_cpools_to_fire_vr_col .= 0.5e-6

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
    # Test column-level soil pool updates from fire
    # ================================================================
    @testset "c_state_update3! column soil pools" begin
        d = make_cstate_update3_data()

        initial_decomp = copy(d.cs_soil.decomp_cpools_vr_col)

        CLM.c_state_update3!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        # CWD pool should increase from fire_mortality_c_to_cwdc
        # then decrease from m_decomp_cpools_to_fire_vr
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = initial_decomp[c, j, d.i_cwd] +
                       d.cf_veg.fire_mortality_c_to_cwdc_col[c, j] * d.dt -
                       d.cf_veg.m_decomp_cpools_to_fire_vr_col[c, j, d.i_cwd] * d.dt
            @test d.cs_soil.decomp_cpools_vr_col[c, j, d.i_cwd] ≈ expected
        end

        # Litter pools should increase from m_c_to_litr_fire
        # then decrease from m_decomp_cpools_to_fire_vr
        for c in 1:d.nc, j in 1:d.nlevdecomp, i in d.i_litr_min:d.i_litr_max
            expected = initial_decomp[c, j, i] +
                       d.cf_veg.m_c_to_litr_fire_col[c, j, i] * d.dt -
                       d.cf_veg.m_decomp_cpools_to_fire_vr_col[c, j, i] * d.dt
            @test d.cs_soil.decomp_cpools_vr_col[c, j, i] ≈ expected
        end

        # Non-litter, non-CWD pools should only decrease from fire losses
        for c in 1:d.nc, j in 1:d.nlevdecomp
            for l in (d.i_litr_max+1):(d.ndecomp_pools)
                l == d.i_cwd && continue
                expected = initial_decomp[c, j, l] -
                           d.cf_veg.m_decomp_cpools_to_fire_vr_col[c, j, l] * d.dt
                @test d.cs_soil.decomp_cpools_vr_col[c, j, l] ≈ expected
            end
        end
    end

    # ================================================================
    # Test patch-level gresp pools (always updated, outside matrixcn guard)
    # ================================================================
    @testset "c_state_update3! gresp pools" begin
        d = make_cstate_update3_data()

        initial_gresp_storage = copy(d.cs_veg.gresp_storage_patch)
        initial_gresp_xfer = copy(d.cs_veg.gresp_xfer_patch)

        CLM.c_state_update3!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        for p in 1:d.np
            expected_gs = initial_gresp_storage[p] -
                d.cf_veg.m_gresp_storage_to_fire_patch[p] * d.dt -
                d.cf_veg.m_gresp_storage_to_litter_fire_patch[p] * d.dt
            @test d.cs_veg.gresp_storage_patch[p] ≈ expected_gs

            expected_gx = initial_gresp_xfer[p] -
                d.cf_veg.m_gresp_xfer_to_fire_patch[p] * d.dt -
                d.cf_veg.m_gresp_xfer_to_litter_fire_patch[p] * d.dt
            @test d.cs_veg.gresp_xfer_patch[p] ≈ expected_gx
        end
    end

    # ================================================================
    # Test patch-level displayed pools (leafc, livestemc with live-to-dead)
    # ================================================================
    @testset "c_state_update3! displayed pools" begin
        d = make_cstate_update3_data()

        initial_leafc = copy(d.cs_veg.leafc_patch)
        initial_livestemc = copy(d.cs_veg.livestemc_patch)
        initial_deadstemc = copy(d.cs_veg.deadstemc_patch)
        initial_livecrootc = copy(d.cs_veg.livecrootc_patch)
        initial_deadcrootc = copy(d.cs_veg.deadcrootc_patch)

        CLM.c_state_update3!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        for p in 1:d.np
            # leafc: fire + litter_fire
            @test d.cs_veg.leafc_patch[p] ≈ initial_leafc[p] -
                d.cf_veg.m_leafc_to_fire_patch[p] * d.dt -
                d.cf_veg.m_leafc_to_litter_fire_patch[p] * d.dt

            # livestemc: fire + litter_fire + live-to-dead
            @test d.cs_veg.livestemc_patch[p] ≈ initial_livestemc[p] -
                d.cf_veg.m_livestemc_to_fire_patch[p] * d.dt -
                d.cf_veg.m_livestemc_to_litter_fire_patch[p] * d.dt -
                d.cf_veg.m_livestemc_to_deadstemc_fire_patch[p] * d.dt

            # deadstemc: fire + litter_fire - live-to-dead (gain)
            @test d.cs_veg.deadstemc_patch[p] ≈ initial_deadstemc[p] -
                d.cf_veg.m_deadstemc_to_fire_patch[p] * d.dt -
                d.cf_veg.m_deadstemc_to_litter_fire_patch[p] * d.dt +
                d.cf_veg.m_livestemc_to_deadstemc_fire_patch[p] * d.dt

            # livecrootc: fire + litter_fire + live-to-dead
            @test d.cs_veg.livecrootc_patch[p] ≈ initial_livecrootc[p] -
                d.cf_veg.m_livecrootc_to_fire_patch[p] * d.dt -
                d.cf_veg.m_livecrootc_to_litter_fire_patch[p] * d.dt -
                d.cf_veg.m_livecrootc_to_deadcrootc_fire_patch[p] * d.dt

            # deadcrootc: fire + litter_fire - live-to-dead (gain)
            @test d.cs_veg.deadcrootc_patch[p] ≈ initial_deadcrootc[p] -
                d.cf_veg.m_deadcrootc_to_fire_patch[p] * d.dt -
                d.cf_veg.m_deadcrootc_to_litter_fire_patch[p] * d.dt +
                d.cf_veg.m_livecrootc_to_deadcrootc_fire_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test patch-level storage and transfer pools
    # ================================================================
    @testset "c_state_update3! storage and transfer pools" begin
        d = make_cstate_update3_data()

        initial_leafc_storage = copy(d.cs_veg.leafc_storage_patch)
        initial_leafc_xfer = copy(d.cs_veg.leafc_xfer_patch)
        initial_deadcrootc_storage = copy(d.cs_veg.deadcrootc_storage_patch)
        initial_deadcrootc_xfer = copy(d.cs_veg.deadcrootc_xfer_patch)

        CLM.c_state_update3!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        for p in 1:d.np
            @test d.cs_veg.leafc_storage_patch[p] ≈
                  initial_leafc_storage[p] -
                  d.cf_veg.m_leafc_storage_to_fire_patch[p] * d.dt -
                  d.cf_veg.m_leafc_storage_to_litter_fire_patch[p] * d.dt
            @test d.cs_veg.leafc_xfer_patch[p] ≈
                  initial_leafc_xfer[p] -
                  d.cf_veg.m_leafc_xfer_to_fire_patch[p] * d.dt -
                  d.cf_veg.m_leafc_xfer_to_litter_fire_patch[p] * d.dt
            @test d.cs_veg.deadcrootc_storage_patch[p] ≈
                  initial_deadcrootc_storage[p] -
                  d.cf_veg.m_deadcrootc_storage_to_fire_patch[p] * d.dt -
                  d.cf_veg.m_deadcrootc_storage_to_litter_fire_patch[p] * d.dt
            @test d.cs_veg.deadcrootc_xfer_patch[p] ≈
                  initial_deadcrootc_xfer[p] -
                  d.cf_veg.m_deadcrootc_xfer_to_fire_patch[p] * d.dt -
                  d.cf_veg.m_deadcrootc_xfer_to_litter_fire_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test mask filtering — nothing should change when masks are false
    # ================================================================
    @testset "c_state_update3! mask filtering" begin
        d = make_cstate_update3_data()

        d.mask_soilc .= false
        d.mask_soilp .= false

        initial_decomp = copy(d.cs_soil.decomp_cpools_vr_col)
        initial_leafc = copy(d.cs_veg.leafc_patch)
        initial_gresp_storage = copy(d.cs_veg.gresp_storage_patch)

        CLM.c_state_update3!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        @test d.cs_soil.decomp_cpools_vr_col == initial_decomp
        @test d.cs_veg.leafc_patch == initial_leafc
        @test d.cs_veg.gresp_storage_patch == initial_gresp_storage
    end

    # ================================================================
    # Test use_matrixcn skips patch veg displayed/storage/xfer updates
    # ================================================================
    @testset "c_state_update3! use_matrixcn skips veg" begin
        d = make_cstate_update3_data()

        initial_leafc = copy(d.cs_veg.leafc_patch)
        initial_leafc_storage = copy(d.cs_veg.leafc_storage_patch)
        initial_leafc_xfer = copy(d.cs_veg.leafc_xfer_patch)
        initial_gresp_storage = copy(d.cs_veg.gresp_storage_patch)

        CLM.c_state_update3!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            use_matrixcn=true,
            dt=d.dt)

        # With use_matrixcn, veg displayed/storage/xfer pools should NOT change
        @test d.cs_veg.leafc_patch == initial_leafc
        @test d.cs_veg.leafc_storage_patch == initial_leafc_storage
        @test d.cs_veg.leafc_xfer_patch == initial_leafc_xfer

        # But gresp_storage/xfer still change (outside the matrixcn guard)
        for p in 1:d.np
            @test d.cs_veg.gresp_storage_patch[p] ≈
                  initial_gresp_storage[p] -
                  d.cf_veg.m_gresp_storage_to_fire_patch[p] * d.dt -
                  d.cf_veg.m_gresp_storage_to_litter_fire_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test use_soil_matrixcn skips soil pool updates
    # ================================================================
    @testset "c_state_update3! use_soil_matrixcn skips soil" begin
        d = make_cstate_update3_data()

        initial_decomp = copy(d.cs_soil.decomp_cpools_vr_col)

        CLM.c_state_update3!(d.cs_veg, d.cf_veg, d.cs_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            use_soil_matrixcn=true,
            dt=d.dt)

        # With use_soil_matrixcn, soil pools should NOT change
        @test d.cs_soil.decomp_cpools_vr_col == initial_decomp
    end

end
