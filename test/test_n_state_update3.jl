@testset "N State Update 3" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for testing NStateUpdate3
    # ----------------------------------------------------------------
    function make_nstate_update3_data(; nc=4, np=6, ng=2, nlevdecomp=1,
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

        # Fire mortality fluxes — patch level to combustion
        nf_veg.m_leafn_to_fire_patch              .= 1.0e-6
        nf_veg.m_frootn_to_fire_patch             .= 0.8e-6
        nf_veg.m_livestemn_to_fire_patch          .= 0.5e-6
        nf_veg.m_deadstemn_to_fire_patch          .= 0.3e-6
        nf_veg.m_livecrootn_to_fire_patch         .= 0.4e-6
        nf_veg.m_deadcrootn_to_fire_patch         .= 0.2e-6
        nf_veg.m_leafn_storage_to_fire_patch      .= 0.1e-6
        nf_veg.m_frootn_storage_to_fire_patch     .= 0.08e-6
        nf_veg.m_livestemn_storage_to_fire_patch  .= 0.05e-6
        nf_veg.m_deadstemn_storage_to_fire_patch  .= 0.03e-6
        nf_veg.m_livecrootn_storage_to_fire_patch .= 0.04e-6
        nf_veg.m_deadcrootn_storage_to_fire_patch .= 0.02e-6
        nf_veg.m_leafn_xfer_to_fire_patch         .= 0.05e-6
        nf_veg.m_frootn_xfer_to_fire_patch        .= 0.04e-6
        nf_veg.m_livestemn_xfer_to_fire_patch     .= 0.03e-6
        nf_veg.m_deadstemn_xfer_to_fire_patch     .= 0.02e-6
        nf_veg.m_livecrootn_xfer_to_fire_patch    .= 0.015e-6
        nf_veg.m_deadcrootn_xfer_to_fire_patch    .= 0.01e-6
        nf_veg.m_retransn_to_fire_patch           .= 0.12e-6

        # Fire mortality fluxes — patch level to litter
        nf_veg.m_leafn_to_litter_fire_patch              .= 0.5e-6
        nf_veg.m_frootn_to_litter_fire_patch             .= 0.4e-6
        nf_veg.m_livestemn_to_litter_fire_patch          .= 0.3e-6
        nf_veg.m_deadstemn_to_litter_fire_patch          .= 0.15e-6
        nf_veg.m_livecrootn_to_litter_fire_patch         .= 0.2e-6
        nf_veg.m_deadcrootn_to_litter_fire_patch         .= 0.1e-6
        nf_veg.m_leafn_storage_to_litter_fire_patch      .= 0.05e-6
        nf_veg.m_frootn_storage_to_litter_fire_patch     .= 0.04e-6
        nf_veg.m_livestemn_storage_to_litter_fire_patch  .= 0.03e-6
        nf_veg.m_deadstemn_storage_to_litter_fire_patch  .= 0.015e-6
        nf_veg.m_livecrootn_storage_to_litter_fire_patch .= 0.02e-6
        nf_veg.m_deadcrootn_storage_to_litter_fire_patch .= 0.01e-6
        nf_veg.m_leafn_xfer_to_litter_fire_patch         .= 0.025e-6
        nf_veg.m_frootn_xfer_to_litter_fire_patch        .= 0.02e-6
        nf_veg.m_livestemn_xfer_to_litter_fire_patch     .= 0.015e-6
        nf_veg.m_deadstemn_xfer_to_litter_fire_patch     .= 0.01e-6
        nf_veg.m_livecrootn_xfer_to_litter_fire_patch    .= 0.008e-6
        nf_veg.m_deadcrootn_xfer_to_litter_fire_patch    .= 0.005e-6
        nf_veg.m_retransn_to_litter_fire_patch           .= 0.06e-6

        # Live-to-dead fire transfers
        nf_veg.m_livestemn_to_deadstemn_fire_patch       .= 0.2e-6
        nf_veg.m_livecrootn_to_deadcrootn_fire_patch     .= 0.15e-6

        # Column-level fire fluxes
        nf_veg.fire_mortality_n_to_cwdn_col   .= 1.0e-6
        nf_veg.m_n_to_litr_fire_col           .= 2.0e-6
        nf_veg.m_decomp_npools_to_fire_vr_col .= 0.5e-6

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
    # Test column-level soil pool updates from fire
    # ================================================================
    @testset "n_state_update3! column soil pools" begin
        d = make_nstate_update3_data()

        initial_decomp = copy(d.ns_soil.decomp_npools_vr_col)

        CLM.n_state_update3!(d.ns_veg, d.nf_veg, d.ns_soil;
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

        # CWD pool should increase from fire_mortality_n_to_cwdn
        # then decrease from m_decomp_npools_to_fire_vr
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected = initial_decomp[c, j, d.i_cwd] +
                       d.nf_veg.fire_mortality_n_to_cwdn_col[c, j] * d.dt -
                       d.nf_veg.m_decomp_npools_to_fire_vr_col[c, j, d.i_cwd] * d.dt
            @test d.ns_soil.decomp_npools_vr_col[c, j, d.i_cwd] ≈ expected
        end

        # Litter pools should increase from m_n_to_litr_fire
        # then decrease from m_decomp_npools_to_fire_vr
        for c in 1:d.nc, j in 1:d.nlevdecomp, i in d.i_litr_min:d.i_litr_max
            expected = initial_decomp[c, j, i] +
                       d.nf_veg.m_n_to_litr_fire_col[c, j, i] * d.dt -
                       d.nf_veg.m_decomp_npools_to_fire_vr_col[c, j, i] * d.dt
            @test d.ns_soil.decomp_npools_vr_col[c, j, i] ≈ expected
        end

        # Non-litter, non-CWD pools should only decrease from fire losses
        for c in 1:d.nc, j in 1:d.nlevdecomp
            for l in (d.i_litr_max+1):(d.ndecomp_pools)
                l == d.i_cwd && continue
                expected = initial_decomp[c, j, l] -
                           d.nf_veg.m_decomp_npools_to_fire_vr_col[c, j, l] * d.dt
                @test d.ns_soil.decomp_npools_vr_col[c, j, l] ≈ expected
            end
        end
    end

    # ================================================================
    # Test patch-level displayed pools (with live-to-dead transfers)
    # ================================================================
    @testset "n_state_update3! displayed pools" begin
        d = make_nstate_update3_data()

        initial_leafn = copy(d.ns_veg.leafn_patch)
        initial_livestemn = copy(d.ns_veg.livestemn_patch)
        initial_deadstemn = copy(d.ns_veg.deadstemn_patch)
        initial_livecrootn = copy(d.ns_veg.livecrootn_patch)
        initial_deadcrootn = copy(d.ns_veg.deadcrootn_patch)
        initial_retransn = copy(d.ns_veg.retransn_patch)

        CLM.n_state_update3!(d.ns_veg, d.nf_veg, d.ns_soil;
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
            # leafn: fire + litter_fire
            @test d.ns_veg.leafn_patch[p] ≈ initial_leafn[p] -
                d.nf_veg.m_leafn_to_fire_patch[p] * d.dt -
                d.nf_veg.m_leafn_to_litter_fire_patch[p] * d.dt

            # livestemn: fire + litter_fire + live-to-dead
            @test d.ns_veg.livestemn_patch[p] ≈ initial_livestemn[p] -
                d.nf_veg.m_livestemn_to_fire_patch[p] * d.dt -
                d.nf_veg.m_livestemn_to_litter_fire_patch[p] * d.dt -
                d.nf_veg.m_livestemn_to_deadstemn_fire_patch[p] * d.dt

            # deadstemn: fire + litter_fire - live-to-dead (gain)
            @test d.ns_veg.deadstemn_patch[p] ≈ initial_deadstemn[p] -
                d.nf_veg.m_deadstemn_to_fire_patch[p] * d.dt -
                d.nf_veg.m_deadstemn_to_litter_fire_patch[p] * d.dt +
                d.nf_veg.m_livestemn_to_deadstemn_fire_patch[p] * d.dt

            # livecrootn: fire + litter_fire + live-to-dead
            @test d.ns_veg.livecrootn_patch[p] ≈ initial_livecrootn[p] -
                d.nf_veg.m_livecrootn_to_fire_patch[p] * d.dt -
                d.nf_veg.m_livecrootn_to_litter_fire_patch[p] * d.dt -
                d.nf_veg.m_livecrootn_to_deadcrootn_fire_patch[p] * d.dt

            # deadcrootn: fire + litter_fire - live-to-dead (gain)
            @test d.ns_veg.deadcrootn_patch[p] ≈ initial_deadcrootn[p] -
                d.nf_veg.m_deadcrootn_to_fire_patch[p] * d.dt -
                d.nf_veg.m_deadcrootn_to_litter_fire_patch[p] * d.dt +
                d.nf_veg.m_livecrootn_to_deadcrootn_fire_patch[p] * d.dt

            # retransn: fire + litter_fire
            @test d.ns_veg.retransn_patch[p] ≈ initial_retransn[p] -
                d.nf_veg.m_retransn_to_fire_patch[p] * d.dt -
                d.nf_veg.m_retransn_to_litter_fire_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test patch-level storage and transfer pools
    # ================================================================
    @testset "n_state_update3! storage and transfer pools" begin
        d = make_nstate_update3_data()

        initial_leafn_storage = copy(d.ns_veg.leafn_storage_patch)
        initial_leafn_xfer = copy(d.ns_veg.leafn_xfer_patch)
        initial_deadcrootn_storage = copy(d.ns_veg.deadcrootn_storage_patch)
        initial_deadcrootn_xfer = copy(d.ns_veg.deadcrootn_xfer_patch)

        CLM.n_state_update3!(d.ns_veg, d.nf_veg, d.ns_soil;
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
            @test d.ns_veg.leafn_storage_patch[p] ≈
                  initial_leafn_storage[p] -
                  d.nf_veg.m_leafn_storage_to_fire_patch[p] * d.dt -
                  d.nf_veg.m_leafn_storage_to_litter_fire_patch[p] * d.dt
            @test d.ns_veg.leafn_xfer_patch[p] ≈
                  initial_leafn_xfer[p] -
                  d.nf_veg.m_leafn_xfer_to_fire_patch[p] * d.dt -
                  d.nf_veg.m_leafn_xfer_to_litter_fire_patch[p] * d.dt
            @test d.ns_veg.deadcrootn_storage_patch[p] ≈
                  initial_deadcrootn_storage[p] -
                  d.nf_veg.m_deadcrootn_storage_to_fire_patch[p] * d.dt -
                  d.nf_veg.m_deadcrootn_storage_to_litter_fire_patch[p] * d.dt
            @test d.ns_veg.deadcrootn_xfer_patch[p] ≈
                  initial_deadcrootn_xfer[p] -
                  d.nf_veg.m_deadcrootn_xfer_to_fire_patch[p] * d.dt -
                  d.nf_veg.m_deadcrootn_xfer_to_litter_fire_patch[p] * d.dt
        end
    end

    # ================================================================
    # Test mask filtering — nothing should change when masks are false
    # ================================================================
    @testset "n_state_update3! mask filtering" begin
        d = make_nstate_update3_data()

        d.mask_soilc .= false
        d.mask_soilp .= false

        initial_decomp = copy(d.ns_soil.decomp_npools_vr_col)
        initial_leafn = copy(d.ns_veg.leafn_patch)
        initial_retransn = copy(d.ns_veg.retransn_patch)

        CLM.n_state_update3!(d.ns_veg, d.nf_veg, d.ns_soil;
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

        @test d.ns_soil.decomp_npools_vr_col == initial_decomp
        @test d.ns_veg.leafn_patch == initial_leafn
        @test d.ns_veg.retransn_patch == initial_retransn
    end

    # ================================================================
    # Test use_matrixcn skips patch veg updates
    # ================================================================
    @testset "n_state_update3! use_matrixcn skips veg" begin
        d = make_nstate_update3_data()

        initial_leafn = copy(d.ns_veg.leafn_patch)
        initial_leafn_storage = copy(d.ns_veg.leafn_storage_patch)
        initial_leafn_xfer = copy(d.ns_veg.leafn_xfer_patch)
        initial_retransn = copy(d.ns_veg.retransn_patch)

        CLM.n_state_update3!(d.ns_veg, d.nf_veg, d.ns_soil;
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

        # With use_matrixcn, veg displayed/storage/xfer/retrans pools should NOT change
        @test d.ns_veg.leafn_patch == initial_leafn
        @test d.ns_veg.leafn_storage_patch == initial_leafn_storage
        @test d.ns_veg.leafn_xfer_patch == initial_leafn_xfer
        @test d.ns_veg.retransn_patch == initial_retransn
    end

    # ================================================================
    # Test use_soil_matrixcn skips soil pool updates
    # ================================================================
    @testset "n_state_update3! use_soil_matrixcn skips soil" begin
        d = make_nstate_update3_data()

        initial_decomp = copy(d.ns_soil.decomp_npools_vr_col)

        CLM.n_state_update3!(d.ns_veg, d.nf_veg, d.ns_soil;
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
        @test d.ns_soil.decomp_npools_vr_col == initial_decomp
    end

    # ================================================================
    # Test NStateUpdateLeaching — without nitrif_denitrif
    # ================================================================
    @testset "n_state_update_leaching! without nitrif_denitrif" begin
        nc = 4
        nlevdecomp = 2
        ndecomp_pools = 7
        dt = 1800.0

        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, nc, 2, nlevdecomp, ndecomp_pools)
        ns_soil.sminn_vr_col .= 10.0

        nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf_soil, nc, nlevdecomp, ndecomp_pools, 1)
        nf_soil.sminn_leached_vr_col .= 1.0e-6

        mask_soilc = trues(nc)
        initial_sminn = copy(ns_soil.sminn_vr_col)

        CLM.n_state_update_leaching!(ns_soil, nf_soil;
            mask_soilc=mask_soilc,
            bounds_col=1:nc,
            nlevdecomp=nlevdecomp,
            use_nitrif_denitrif=false,
            dt=dt)

        for c in 1:nc, j in 1:nlevdecomp
            expected = initial_sminn[c, j] - nf_soil.sminn_leached_vr_col[c, j] * dt
            @test ns_soil.sminn_vr_col[c, j] ≈ expected
        end
    end

    # ================================================================
    # Test NStateUpdateLeaching — with nitrif_denitrif
    # ================================================================
    @testset "n_state_update_leaching! with nitrif_denitrif" begin
        nc = 4
        nlevdecomp = 2
        ndecomp_pools = 7
        dt = 1800.0

        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, nc, 2, nlevdecomp, ndecomp_pools)
        ns_soil.smin_no3_vr_col .= 5.0
        ns_soil.smin_nh4_vr_col .= 3.0
        ns_soil.sminn_vr_col    .= 8.0

        nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf_soil, nc, nlevdecomp, ndecomp_pools, 1)
        nf_soil.smin_no3_leached_vr_col .= 0.5e-6
        nf_soil.smin_no3_runoff_vr_col  .= 0.3e-6

        mask_soilc = trues(nc)
        initial_no3 = copy(ns_soil.smin_no3_vr_col)
        initial_nh4 = copy(ns_soil.smin_nh4_vr_col)

        CLM.n_state_update_leaching!(ns_soil, nf_soil;
            mask_soilc=mask_soilc,
            bounds_col=1:nc,
            nlevdecomp=nlevdecomp,
            use_nitrif_denitrif=true,
            dt=dt)

        for c in 1:nc, j in 1:nlevdecomp
            expected_no3 = max(initial_no3[c, j] -
                (nf_soil.smin_no3_leached_vr_col[c, j] + nf_soil.smin_no3_runoff_vr_col[c, j]) * dt, 0.0)
            @test ns_soil.smin_no3_vr_col[c, j] ≈ expected_no3
            @test ns_soil.sminn_vr_col[c, j] ≈ expected_no3 + initial_nh4[c, j]
        end
    end

    # ================================================================
    # Test NStateUpdateLeaching — max(0) clamp
    # ================================================================
    @testset "n_state_update_leaching! no3 clamped at zero" begin
        nc = 2
        nlevdecomp = 1
        ndecomp_pools = 7
        dt = 1800.0

        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, nc, 2, nlevdecomp, ndecomp_pools)
        ns_soil.smin_no3_vr_col .= 1.0e-6  # very small
        ns_soil.smin_nh4_vr_col .= 2.0

        nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf_soil, nc, nlevdecomp, ndecomp_pools, 1)
        nf_soil.smin_no3_leached_vr_col .= 1.0  # large leaching to force clamp
        nf_soil.smin_no3_runoff_vr_col  .= 1.0

        mask_soilc = trues(nc)

        CLM.n_state_update_leaching!(ns_soil, nf_soil;
            mask_soilc=mask_soilc,
            bounds_col=1:nc,
            nlevdecomp=nlevdecomp,
            use_nitrif_denitrif=true,
            dt=dt)

        for c in 1:nc, j in 1:nlevdecomp
            @test ns_soil.smin_no3_vr_col[c, j] == 0.0
            @test ns_soil.sminn_vr_col[c, j] ≈ 2.0  # nh4 only
        end
    end

    # ================================================================
    # Test NStateUpdateLeaching — mask filtering
    # ================================================================
    @testset "n_state_update_leaching! mask filtering" begin
        nc = 4
        nlevdecomp = 1
        ndecomp_pools = 7
        dt = 1800.0

        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, nc, 2, nlevdecomp, ndecomp_pools)
        ns_soil.sminn_vr_col .= 10.0

        nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf_soil, nc, nlevdecomp, ndecomp_pools, 1)
        nf_soil.sminn_leached_vr_col .= 1.0e-6

        mask_soilc = falses(nc)
        initial_sminn = copy(ns_soil.sminn_vr_col)

        CLM.n_state_update_leaching!(ns_soil, nf_soil;
            mask_soilc=mask_soilc,
            bounds_col=1:nc,
            nlevdecomp=nlevdecomp,
            use_nitrif_denitrif=false,
            dt=dt)

        @test ns_soil.sminn_vr_col == initial_sminn
    end

end
