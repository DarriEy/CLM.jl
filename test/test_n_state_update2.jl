@testset "N State Update 2" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for testing NStateUpdate2
    # ----------------------------------------------------------------
    function make_nstate_update2_data(; nc=4, np=6, ng=2, nlevdecomp=1,
                                       ndecomp_pools=7, nrepr=1)
        i_litr_min = 1
        i_litr_max = 3
        i_cwd = 4
        dt = 1800.0  # 30 minutes

        # --- CNVeg nitrogen state ---
        ns_veg = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns_veg, np, nc, ng; nrepr=nrepr)
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

        # --- CNVeg nitrogen flux ---
        nf_veg = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf_veg, np, nc, ng;
                                       nrepr=nrepr,
                                       nlevdecomp_full=nlevdecomp,
                                       ndecomp_pools=ndecomp_pools,
                                       i_litr_max=i_litr_max)

        # Gap mortality fluxes (NStateUpdate2)
        nf_veg.m_leafn_to_litter_patch              .= 1.0e-6
        nf_veg.m_frootn_to_litter_patch             .= 0.8e-6
        nf_veg.m_livestemn_to_litter_patch          .= 0.5e-6
        nf_veg.m_deadstemn_to_litter_patch          .= 0.3e-6
        nf_veg.m_livecrootn_to_litter_patch         .= 0.4e-6
        nf_veg.m_deadcrootn_to_litter_patch         .= 0.2e-6
        nf_veg.m_retransn_to_litter_patch           .= 0.15e-6
        nf_veg.m_leafn_storage_to_litter_patch      .= 0.1e-6
        nf_veg.m_frootn_storage_to_litter_patch     .= 0.08e-6
        nf_veg.m_livestemn_storage_to_litter_patch  .= 0.05e-6
        nf_veg.m_deadstemn_storage_to_litter_patch  .= 0.03e-6
        nf_veg.m_livecrootn_storage_to_litter_patch .= 0.04e-6
        nf_veg.m_deadcrootn_storage_to_litter_patch .= 0.02e-6
        nf_veg.m_leafn_xfer_to_litter_patch         .= 0.05e-6
        nf_veg.m_frootn_xfer_to_litter_patch        .= 0.04e-6
        nf_veg.m_livestemn_xfer_to_litter_patch     .= 0.03e-6
        nf_veg.m_deadstemn_xfer_to_litter_patch     .= 0.02e-6
        nf_veg.m_livecrootn_xfer_to_litter_patch    .= 0.015e-6
        nf_veg.m_deadcrootn_xfer_to_litter_patch    .= 0.01e-6
        nf_veg.gap_mortality_n_to_litr_n_col        .= 2.0e-6
        nf_veg.gap_mortality_n_to_cwdn_col          .= 1.0e-6

        # Harvest mortality fluxes (NStateUpdate2h)
        nf_veg.hrv_leafn_to_litter_patch             .= 0.5e-6
        nf_veg.hrv_frootn_to_litter_patch            .= 0.4e-6
        nf_veg.hrv_livestemn_to_litter_patch         .= 0.3e-6
        nf_veg.hrv_livecrootn_to_litter_patch        .= 0.2e-6
        nf_veg.hrv_deadcrootn_to_litter_patch        .= 0.1e-6
        nf_veg.wood_harvestn_patch                   .= 0.6e-6
        nf_veg.hrv_retransn_to_litter_patch          .= 0.12e-6
        nf_veg.hrv_leafn_storage_to_litter_patch     .= 0.05e-6
        nf_veg.hrv_frootn_storage_to_litter_patch    .= 0.04e-6
        nf_veg.hrv_livestemn_storage_to_litter_patch .= 0.03e-6
        nf_veg.hrv_deadstemn_storage_to_litter_patch .= 0.02e-6
        nf_veg.hrv_livecrootn_storage_to_litter_patch .= 0.015e-6
        nf_veg.hrv_deadcrootn_storage_to_litter_patch .= 0.01e-6
        nf_veg.hrv_leafn_xfer_to_litter_patch        .= 0.025e-6
        nf_veg.hrv_frootn_xfer_to_litter_patch       .= 0.02e-6
        nf_veg.hrv_livestemn_xfer_to_litter_patch    .= 0.015e-6
        nf_veg.hrv_deadstemn_xfer_to_litter_patch    .= 0.01e-6
        nf_veg.hrv_livecrootn_xfer_to_litter_patch   .= 0.008e-6
        nf_veg.hrv_deadcrootn_xfer_to_litter_patch   .= 0.005e-6
        nf_veg.harvest_n_to_litr_n_col               .= 1.5e-6
        nf_veg.harvest_n_to_cwdn_col                 .= 0.8e-6

        # Gross unrepresented landcover change fluxes (NStateUpdate2g)
        nf_veg.gru_leafn_to_litter_patch            .= 0.3e-6
        nf_veg.gru_frootn_to_litter_patch           .= 0.25e-6
        nf_veg.gru_livestemn_to_atm_patch           .= 0.2e-6
        nf_veg.gru_deadstemn_to_atm_patch           .= 0.15e-6
        nf_veg.gru_wood_productn_gain_patch         .= 0.1e-6
        nf_veg.gru_livecrootn_to_litter_patch       .= 0.12e-6
        nf_veg.gru_deadcrootn_to_litter_patch       .= 0.08e-6
        nf_veg.gru_retransn_to_litter_patch         .= 0.06e-6
        nf_veg.gru_leafn_storage_to_atm_patch       .= 0.03e-6
        nf_veg.gru_frootn_storage_to_atm_patch      .= 0.025e-6
        nf_veg.gru_livestemn_storage_to_atm_patch   .= 0.02e-6
        nf_veg.gru_deadstemn_storage_to_atm_patch   .= 0.015e-6
        nf_veg.gru_livecrootn_storage_to_atm_patch  .= 0.012e-6
        nf_veg.gru_deadcrootn_storage_to_atm_patch  .= 0.008e-6
        nf_veg.gru_leafn_xfer_to_atm_patch          .= 0.015e-6
        nf_veg.gru_frootn_xfer_to_atm_patch         .= 0.012e-6
        nf_veg.gru_livestemn_xfer_to_atm_patch      .= 0.01e-6
        nf_veg.gru_deadstemn_xfer_to_atm_patch      .= 0.008e-6
        nf_veg.gru_livecrootn_xfer_to_atm_patch     .= 0.006e-6
        nf_veg.gru_deadcrootn_xfer_to_atm_patch     .= 0.004e-6
        nf_veg.gru_n_to_litr_n_col                  .= 1.0e-6
        nf_veg.gru_n_to_cwdn_col                    .= 0.5e-6

        # --- Soil biogeochem nitrogen state ---
        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, nc, ng, nlevdecomp, ndecomp_pools)
        ns_soil.decomp_npools_vr_col .= 10.0

        # --- Masks ---
        mask_soilc = trues(nc)
        mask_soilp = trues(np)

        return (ns_veg=ns_veg, nf_veg=nf_veg, ns_soil=ns_soil,
                mask_soilc=mask_soilc, mask_soilp=mask_soilp,
                nc=nc, np=np, ng=ng, dt=dt,
                nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd)
    end

    # ================================================================
    # Test NStateUpdate2 — gap-phase mortality
    # ================================================================
    @testset "n_state_update2! gap mortality" begin
        d = make_nstate_update2_data()

        initial_decomp = copy(d.ns_soil.decomp_npools_vr_col)
        initial_leafn = copy(d.ns_veg.leafn_patch)
        initial_retransn = copy(d.ns_veg.retransn_patch)

        CLM.n_state_update2!(d.ns_veg, d.nf_veg, d.ns_soil;
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
            expected = initial_decomp[c, j, i] + d.nf_veg.gap_mortality_n_to_litr_n_col[c, j, i] * d.dt
            @test d.ns_soil.decomp_npools_vr_col[c, j, i] ≈ expected
        end

        # CWD pool should increase
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = initial_decomp[c, j, d.i_cwd] + d.nf_veg.gap_mortality_n_to_cwdn_col[c, j] * d.dt
            @test d.ns_soil.decomp_npools_vr_col[c, j, d.i_cwd] ≈ expected
        end

        # Patch-level: leafn should decrease from gap mortality
        for p in 1:d.np
            expected = initial_leafn[p] - d.nf_veg.m_leafn_to_litter_patch[p] * d.dt
            @test d.ns_veg.leafn_patch[p] ≈ expected
        end

        # retransn should decrease
        for p in 1:d.np
            expected = initial_retransn[p] - d.nf_veg.m_retransn_to_litter_patch[p] * d.dt
            @test d.ns_veg.retransn_patch[p] ≈ expected
        end
    end

    # ================================================================
    # Test NStateUpdate2 — storage and transfer pools
    # ================================================================
    @testset "n_state_update2! storage and transfer pools" begin
        d = make_nstate_update2_data()

        initial_leafn_storage = copy(d.ns_veg.leafn_storage_patch)
        initial_leafn_xfer = copy(d.ns_veg.leafn_xfer_patch)
        initial_deadcrootn_storage = copy(d.ns_veg.deadcrootn_storage_patch)
        initial_deadcrootn_xfer = copy(d.ns_veg.deadcrootn_xfer_patch)

        CLM.n_state_update2!(d.ns_veg, d.nf_veg, d.ns_soil;
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
            @test d.ns_veg.leafn_storage_patch[p] ≈
                  initial_leafn_storage[p] - d.nf_veg.m_leafn_storage_to_litter_patch[p] * d.dt
            @test d.ns_veg.leafn_xfer_patch[p] ≈
                  initial_leafn_xfer[p] - d.nf_veg.m_leafn_xfer_to_litter_patch[p] * d.dt
            @test d.ns_veg.deadcrootn_storage_patch[p] ≈
                  initial_deadcrootn_storage[p] - d.nf_veg.m_deadcrootn_storage_to_litter_patch[p] * d.dt
            @test d.ns_veg.deadcrootn_xfer_patch[p] ≈
                  initial_deadcrootn_xfer[p] - d.nf_veg.m_deadcrootn_xfer_to_litter_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test NStateUpdate2 — mask filtering
    # ================================================================
    @testset "n_state_update2! mask filtering" begin
        d = make_nstate_update2_data()

        # Mask out everything
        d.mask_soilc .= false
        d.mask_soilp .= false

        initial_decomp = copy(d.ns_soil.decomp_npools_vr_col)
        initial_leafn = copy(d.ns_veg.leafn_patch)

        CLM.n_state_update2!(d.ns_veg, d.nf_veg, d.ns_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        @test d.ns_soil.decomp_npools_vr_col == initial_decomp
        @test d.ns_veg.leafn_patch == initial_leafn
    end

    # ================================================================
    # Test NStateUpdate2h — harvest mortality
    # ================================================================
    @testset "n_state_update2h! harvest mortality" begin
        d = make_nstate_update2_data()

        initial_decomp = copy(d.ns_soil.decomp_npools_vr_col)
        initial_leafn = copy(d.ns_veg.leafn_patch)
        initial_deadstemn = copy(d.ns_veg.deadstemn_patch)
        initial_retransn = copy(d.ns_veg.retransn_patch)

        CLM.n_state_update2h!(d.ns_veg, d.nf_veg, d.ns_soil;
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
            expected = initial_decomp[c, j, i] + d.nf_veg.harvest_n_to_litr_n_col[c, j, i] * d.dt
            @test d.ns_soil.decomp_npools_vr_col[c, j, i] ≈ expected
        end

        # CWD pool should increase
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = initial_decomp[c, j, d.i_cwd] + d.nf_veg.harvest_n_to_cwdn_col[c, j] * d.dt
            @test d.ns_soil.decomp_npools_vr_col[c, j, d.i_cwd] ≈ expected
        end

        # leafn should decrease
        for p in 1:d.np
            @test d.ns_veg.leafn_patch[p] ≈ initial_leafn[p] - d.nf_veg.hrv_leafn_to_litter_patch[p] * d.dt
        end

        # deadstemn uses wood_harvestn_patch
        for p in 1:d.np
            @test d.ns_veg.deadstemn_patch[p] ≈ initial_deadstemn[p] - d.nf_veg.wood_harvestn_patch[p] * d.dt
        end

        # retransn should decrease
        for p in 1:d.np
            @test d.ns_veg.retransn_patch[p] ≈ initial_retransn[p] - d.nf_veg.hrv_retransn_to_litter_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test NStateUpdate2h — harvest storage and transfer pools
    # ================================================================
    @testset "n_state_update2h! harvest storage/transfer" begin
        d = make_nstate_update2_data()

        initial_leafn_storage = copy(d.ns_veg.leafn_storage_patch)
        initial_deadcrootn_xfer = copy(d.ns_veg.deadcrootn_xfer_patch)

        CLM.n_state_update2h!(d.ns_veg, d.nf_veg, d.ns_soil;
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
            @test d.ns_veg.leafn_storage_patch[p] ≈
                  initial_leafn_storage[p] - d.nf_veg.hrv_leafn_storage_to_litter_patch[p] * d.dt
            @test d.ns_veg.deadcrootn_xfer_patch[p] ≈
                  initial_deadcrootn_xfer[p] - d.nf_veg.hrv_deadcrootn_xfer_to_litter_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test NStateUpdate2g — gross unrepresented landcover change
    # ================================================================
    @testset "n_state_update2g! gross unrep LCC" begin
        d = make_nstate_update2_data()

        initial_decomp = copy(d.ns_soil.decomp_npools_vr_col)
        initial_leafn = copy(d.ns_veg.leafn_patch)
        initial_deadstemn = copy(d.ns_veg.deadstemn_patch)
        initial_retransn = copy(d.ns_veg.retransn_patch)

        CLM.n_state_update2g!(d.ns_veg, d.nf_veg, d.ns_soil;
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
            expected = initial_decomp[c, j, i] + d.nf_veg.gru_n_to_litr_n_col[c, j, i] * d.dt
            @test d.ns_soil.decomp_npools_vr_col[c, j, i] ≈ expected
        end

        # CWD pool should increase
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = initial_decomp[c, j, d.i_cwd] + d.nf_veg.gru_n_to_cwdn_col[c, j] * d.dt
            @test d.ns_soil.decomp_npools_vr_col[c, j, d.i_cwd] ≈ expected
        end

        # leafn should decrease
        for p in 1:d.np
            @test d.ns_veg.leafn_patch[p] ≈ initial_leafn[p] - d.nf_veg.gru_leafn_to_litter_patch[p] * d.dt
        end

        # deadstemn: two separate decrements (gru_deadstemn_to_atm + gru_wood_productn_gain)
        for p in 1:d.np
            expected = initial_deadstemn[p] -
                       d.nf_veg.gru_deadstemn_to_atm_patch[p] * d.dt -
                       d.nf_veg.gru_wood_productn_gain_patch[p] * d.dt
            @test d.ns_veg.deadstemn_patch[p] ≈ expected
        end

        # retransn should decrease
        for p in 1:d.np
            @test d.ns_veg.retransn_patch[p] ≈ initial_retransn[p] - d.nf_veg.gru_retransn_to_litter_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test NStateUpdate2g — storage and transfer pools
    # ================================================================
    @testset "n_state_update2g! gru storage/transfer" begin
        d = make_nstate_update2_data()

        initial_leafn_storage = copy(d.ns_veg.leafn_storage_patch)
        initial_leafn_xfer = copy(d.ns_veg.leafn_xfer_patch)
        initial_deadcrootn_storage = copy(d.ns_veg.deadcrootn_storage_patch)
        initial_deadcrootn_xfer = copy(d.ns_veg.deadcrootn_xfer_patch)

        CLM.n_state_update2g!(d.ns_veg, d.nf_veg, d.ns_soil;
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
            @test d.ns_veg.leafn_storage_patch[p] ≈
                  initial_leafn_storage[p] - d.nf_veg.gru_leafn_storage_to_atm_patch[p] * d.dt
            @test d.ns_veg.leafn_xfer_patch[p] ≈
                  initial_leafn_xfer[p] - d.nf_veg.gru_leafn_xfer_to_atm_patch[p] * d.dt
            @test d.ns_veg.deadcrootn_storage_patch[p] ≈
                  initial_deadcrootn_storage[p] - d.nf_veg.gru_deadcrootn_storage_to_atm_patch[p] * d.dt
            @test d.ns_veg.deadcrootn_xfer_patch[p] ≈
                  initial_deadcrootn_xfer[p] - d.nf_veg.gru_deadcrootn_xfer_to_atm_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test NStateUpdate2g — mask filtering
    # ================================================================
    @testset "n_state_update2g! mask filtering" begin
        d = make_nstate_update2_data()

        d.mask_soilc .= false
        d.mask_soilp .= false

        initial_decomp = copy(d.ns_soil.decomp_npools_vr_col)
        initial_leafn = copy(d.ns_veg.leafn_patch)

        CLM.n_state_update2g!(d.ns_veg, d.nf_veg, d.ns_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        @test d.ns_soil.decomp_npools_vr_col == initial_decomp
        @test d.ns_veg.leafn_patch == initial_leafn
    end

    # ================================================================
    # Test use_matrixcn skips patch veg updates
    # ================================================================
    @testset "n_state_update2! use_matrixcn skips veg" begin
        d = make_nstate_update2_data()

        initial_leafn = copy(d.ns_veg.leafn_patch)
        initial_retransn = copy(d.ns_veg.retransn_patch)

        CLM.n_state_update2!(d.ns_veg, d.nf_veg, d.ns_soil;
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
        @test d.ns_veg.leafn_patch == initial_leafn
        @test d.ns_veg.retransn_patch == initial_retransn
    end

end
