@testset "N State Update 1" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for testing
    # ----------------------------------------------------------------
    function make_nstate_update_data(; nc=4, np=6, ng=2, nlevdecomp=1,
                                      ndecomp_pools=7, nrepr=1)
        i_litr_min = 1
        i_litr_max = 3
        i_cwd = 4
        dt = 1800.0  # 30 minutes

        # --- CNVeg nitrogen state ---
        ns_veg = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns_veg, np, nc, ng; nrepr=nrepr)
        # Initialize pools to known values
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
        ns_veg.npool_patch             .= 10.0
        ns_veg.seedn_grc               .= 100.0
        ns_veg.cropseedn_deficit_patch .= -0.5
        for k in 1:nrepr
            ns_veg.reproductiven_patch[:, k]         .= 1.0
            ns_veg.reproductiven_storage_patch[:, k] .= 0.2
            ns_veg.reproductiven_xfer_patch[:, k]    .= 0.1
        end

        # --- CNVeg nitrogen flux ---
        nf_veg = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf_veg, np, nc, ng;
                                       nrepr=nrepr,
                                       nlevdecomp_full=nlevdecomp,
                                       ndecomp_pools=ndecomp_pools,
                                       i_litr_max=i_litr_max)
        # Set all relevant fluxes to small nonzero values
        nf_veg.leafn_xfer_to_leafn_patch        .= 2.0e-7
        nf_veg.frootn_xfer_to_frootn_patch      .= 1.0e-7
        nf_veg.livestemn_xfer_to_livestemn_patch .= 1.5e-7
        nf_veg.deadstemn_xfer_to_deadstemn_patch .= 1.0e-7
        nf_veg.livecrootn_xfer_to_livecrootn_patch .= 0.8e-7
        nf_veg.deadcrootn_xfer_to_deadcrootn_patch .= 0.5e-7
        nf_veg.leafn_to_litter_patch            .= 3.0e-7
        nf_veg.frootn_to_litter_patch           .= 2.0e-7
        nf_veg.leafn_to_retransn_patch          .= 1.0e-7
        nf_veg.livestemn_to_deadstemn_patch     .= 1.0e-8
        nf_veg.livestemn_to_retransn_patch      .= 0.5e-8
        nf_veg.livecrootn_to_deadcrootn_patch   .= 1.0e-8
        nf_veg.livecrootn_to_retransn_patch     .= 0.5e-8
        nf_veg.frootn_to_retransn_patch         .= 0.0
        nf_veg.livestemn_to_litter_patch        .= 0.0
        nf_veg.livestemn_to_biofueln_patch      .= 0.0
        nf_veg.livestemn_to_removedresiduen_patch .= 0.0
        nf_veg.leafn_to_biofueln_patch          .= 0.0
        nf_veg.leafn_to_removedresiduen_patch   .= 0.0
        nf_veg.crop_seedn_to_leaf_patch         .= 0.0
        for k in 1:nrepr
            nf_veg.repr_grainn_to_food_patch[:, k]         .= 0.0
            nf_veg.repr_grainn_to_seed_patch[:, k]         .= 0.0
            nf_veg.repr_structuren_to_cropprod_patch[:, k] .= 0.0
            nf_veg.repr_structuren_to_litter_patch[:, k]   .= 0.0
            nf_veg.reproductiven_xfer_to_reproductiven_patch[:, k] .= 0.0
            nf_veg.reproductiven_storage_to_xfer_patch[:, k] .= 0.0
            nf_veg.npool_to_reproductiven_patch[:, k]      .= 0.0
            nf_veg.npool_to_reproductiven_storage_patch[:, k] .= 0.0
        end
        nf_veg.sminn_to_npool_patch             .= 5.0e-7
        nf_veg.retransn_to_npool_patch          .= 2.0e-7
        nf_veg.free_retransn_to_npool_patch     .= 0.0
        nf_veg.npool_to_leafn_patch             .= 5.0e-7
        nf_veg.npool_to_leafn_storage_patch     .= 2.0e-7
        nf_veg.npool_to_frootn_patch            .= 3.0e-7
        nf_veg.npool_to_frootn_storage_patch    .= 1.0e-7
        nf_veg.npool_to_livestemn_patch         .= 2.0e-7
        nf_veg.npool_to_livestemn_storage_patch .= 1.0e-7
        nf_veg.npool_to_deadstemn_patch         .= 1.5e-7
        nf_veg.npool_to_deadstemn_storage_patch .= 0.5e-7
        nf_veg.npool_to_livecrootn_patch        .= 1.0e-7
        nf_veg.npool_to_livecrootn_storage_patch .= 0.5e-7
        nf_veg.npool_to_deadcrootn_patch        .= 0.5e-7
        nf_veg.npool_to_deadcrootn_storage_patch .= 0.3e-7
        nf_veg.leafn_storage_to_xfer_patch      .= 0.0
        nf_veg.frootn_storage_to_xfer_patch     .= 0.0
        nf_veg.livestemn_storage_to_xfer_patch  .= 0.0
        nf_veg.deadstemn_storage_to_xfer_patch  .= 0.0
        nf_veg.livecrootn_storage_to_xfer_patch .= 0.0
        nf_veg.deadcrootn_storage_to_xfer_patch .= 0.0

        # DynPatch-specific fluxes
        nf_veg.dwt_frootn_to_litr_n_col        .= 0.0
        nf_veg.dwt_livecrootn_to_cwdn_col      .= 0.0
        nf_veg.dwt_deadcrootn_to_cwdn_col      .= 0.0
        nf_veg.dwt_seedn_to_leaf_grc            .= 0.0
        nf_veg.dwt_seedn_to_deadstem_grc        .= 0.0
        nf_veg.phenology_n_to_litr_n_col        .= 0.0

        # --- Soil biogeochem nitrogen flux ---
        nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf_soil, nc, nlevdecomp, ndecomp_pools, 5)
        nf_soil.decomp_npools_sourcesink_col .= 0.0

        # --- Soil biogeochem nitrogen state ---
        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, nc, ng, nlevdecomp, ndecomp_pools)
        ns_soil.decomp_npools_vr_col .= 10.0

        # --- Masks ---
        mask_soilc = trues(nc)
        mask_soilp = trues(np)
        mask_soilc_with_inactive = trues(nc)

        # --- Patch data ---
        ivt          = fill(1, np)  # non-woody, non-crop by default
        woody        = zeros(Float64, 80)
        woody[2]     = 1.0  # PFT 2 is woody
        col_is_fates = fill(false, nc)

        return (ns_veg=ns_veg, nf_veg=nf_veg, nf_soil=nf_soil, ns_soil=ns_soil,
                mask_soilc=mask_soilc, mask_soilp=mask_soilp,
                mask_soilc_with_inactive=mask_soilc_with_inactive,
                ivt=ivt, woody=woody, col_is_fates=col_is_fates,
                nc=nc, np=np, ng=ng, dt=dt,
                nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                nrepr=nrepr, i_litr_min=i_litr_min, i_litr_max=i_litr_max,
                i_cwd=i_cwd)
    end

    # ================================================================
    # Test NStateUpdateDynPatch
    # ================================================================
    @testset "n_state_update_dyn_patch!" begin
        d = make_nstate_update_data()

        # Set some dynamic patch fluxes
        for i in d.i_litr_min:d.i_litr_max
            d.nf_veg.dwt_frootn_to_litr_n_col[:, :, i] .= 1.0e-6
        end
        d.nf_veg.dwt_livecrootn_to_cwdn_col .= 2.0e-6
        d.nf_veg.dwt_deadcrootn_to_cwdn_col .= 3.0e-6
        d.nf_veg.dwt_seedn_to_leaf_grc       .= 0.5e-6
        d.nf_veg.dwt_seedn_to_deadstem_grc   .= 0.3e-6

        initial_seedn = copy(d.ns_veg.seedn_grc)
        initial_decomp = copy(d.ns_soil.decomp_npools_vr_col)

        CLM.n_state_update_dyn_patch!(d.ns_veg, d.nf_veg, d.ns_soil;
            mask_soilc_with_inactive=d.mask_soilc_with_inactive,
            bounds_col=1:d.nc,
            bounds_grc=1:d.ng,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            dt=d.dt)

        # Litter pools should increase
        for c in 1:d.nc, j in 1:d.nlevdecomp, i in d.i_litr_min:d.i_litr_max
            @test d.ns_soil.decomp_npools_vr_col[c, j, i] > initial_decomp[c, j, i]
        end

        # CWD pool should increase from livecroot + deadcroot
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected_delta = (2.0e-6 + 3.0e-6) * d.dt
            @test d.ns_soil.decomp_npools_vr_col[c, j, d.i_cwd] ≈
                  initial_decomp[c, j, d.i_cwd] + expected_delta
        end

        # Seed pool should decrease
        for g in 1:d.ng
            expected = initial_seedn[g] - (0.5e-6 + 0.3e-6) * d.dt
            @test d.ns_veg.seedn_grc[g] ≈ expected
        end
    end

    # ================================================================
    # Test NStateUpdate1 — non-woody, non-crop patch
    # ================================================================
    @testset "n_state_update1! non-woody non-crop" begin
        d = make_nstate_update_data()

        # Set phenology source-sink
        d.nf_veg.phenology_n_to_litr_n_col .= 1.0e-7

        initial_npool  = copy(d.ns_veg.npool_patch)
        initial_leafn  = copy(d.ns_veg.leafn_patch)
        initial_frootn = copy(d.ns_veg.frootn_patch)

        CLM.n_state_update1!(d.ns_veg, d.nf_veg, d.nf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            ivt=d.ivt,
            woody=d.woody,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        # npool should change: gains from sminn and retrans, loses to leaf/froot alloc
        # Net: +sminn(5e-7) +retrans(2e-7) -leafn(5e-7) -leafn_stor(2e-7) -frootn(3e-7) -frootn_stor(1e-7) = -4e-7
        for p in 1:d.np
            net = (d.nf_veg.sminn_to_npool_patch[p] +
                   d.nf_veg.retransn_to_npool_patch[p] -
                   d.nf_veg.npool_to_leafn_patch[p] -
                   d.nf_veg.npool_to_leafn_storage_patch[p] -
                   d.nf_veg.npool_to_frootn_patch[p] -
                   d.nf_veg.npool_to_frootn_storage_patch[p]) * d.dt
            @test d.ns_veg.npool_patch[p] ≈ initial_npool[p] + net
        end

        # leafn should increase: xfer growth + alloc - litterfall - retrans
        for p in 1:d.np
            net = (d.nf_veg.leafn_xfer_to_leafn_patch[p] +
                   d.nf_veg.npool_to_leafn_patch[p] -
                   d.nf_veg.leafn_to_litter_patch[p] -
                   d.nf_veg.leafn_to_retransn_patch[p]) * d.dt
            @test d.ns_veg.leafn_patch[p] ≈ initial_leafn[p] + net
        end

        # Soil decomp sourcesink should be set for litter pools
        for c in 1:d.nc, j in 1:d.nlevdecomp, i in d.i_litr_min:d.i_litr_max
            expected = d.nf_veg.phenology_n_to_litr_n_col[c, j, i] * d.dt
            @test d.nf_soil.decomp_npools_sourcesink_col[c, j, i] ≈ expected atol=1e-15
        end

        # CWD sourcesink should be zero
        for c in 1:d.nc, j in 1:d.nlevdecomp
            @test d.nf_soil.decomp_npools_sourcesink_col[c, j, d.i_cwd] ≈ 0.0 atol=1e-15
        end
    end

    # ================================================================
    # Test NStateUpdate1 — woody patch
    # ================================================================
    @testset "n_state_update1! woody patch" begin
        d = make_nstate_update_data()

        # Make all patches woody (PFT 2)
        d.ivt .= 2

        initial_livestemn  = copy(d.ns_veg.livestemn_patch)
        initial_deadstemn  = copy(d.ns_veg.deadstemn_patch)
        initial_retransn   = copy(d.ns_veg.retransn_patch)
        initial_npool      = copy(d.ns_veg.npool_patch)

        CLM.n_state_update1!(d.ns_veg, d.nf_veg, d.nf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            ivt=d.ivt,
            woody=d.woody,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        # For woody patches: livestemn should change (xfer + alloc - turnover - retrans)
        for p in 1:d.np
            xfer_delta    = d.nf_veg.livestemn_xfer_to_livestemn_patch[p] * d.dt
            alloc_delta   = d.nf_veg.npool_to_livestemn_patch[p] * d.dt
            turnover      = d.nf_veg.livestemn_to_deadstemn_patch[p] * d.dt
            retrans       = d.nf_veg.livestemn_to_retransn_patch[p] * d.dt
            expected = initial_livestemn[p] + xfer_delta + alloc_delta - turnover - retrans
            @test d.ns_veg.livestemn_patch[p] ≈ expected
        end

        # deadstemn should gain from xfer, alloc, and live-to-dead turnover
        for p in 1:d.np
            xfer_delta  = d.nf_veg.deadstemn_xfer_to_deadstemn_patch[p] * d.dt
            alloc_delta = d.nf_veg.npool_to_deadstemn_patch[p] * d.dt
            turnover    = d.nf_veg.livestemn_to_deadstemn_patch[p] * d.dt
            expected = initial_deadstemn[p] + xfer_delta + alloc_delta + turnover
            @test d.ns_veg.deadstemn_patch[p] ≈ expected
        end

        # npool should decrease more for woody (extra alloc to stem/croot pools)
        for p in 1:d.np
            @test d.ns_veg.npool_patch[p] < initial_npool[p]
        end
    end

    # ================================================================
    # Test NStateUpdate1 — FATES column (skip phenology input)
    # ================================================================
    @testset "n_state_update1! FATES column" begin
        d = make_nstate_update_data()

        # Mark column 1 as FATES
        d.col_is_fates[1] = true

        # Set phenology flux for non-FATES columns
        d.nf_veg.phenology_n_to_litr_n_col .= 1.0e-7

        CLM.n_state_update1!(d.ns_veg, d.nf_veg, d.nf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            ivt=d.ivt,
            woody=d.woody,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        # Non-FATES column 2 should have phenology input
        for j in 1:d.nlevdecomp, i in d.i_litr_min:d.i_litr_max
            @test d.nf_soil.decomp_npools_sourcesink_col[2, j, i] ≈ 1.0e-7 * d.dt atol=1e-15
        end
    end

    # ================================================================
    # Test mask filtering: masked columns/patches are skipped
    # ================================================================
    @testset "n_state_update1! mask filtering" begin
        d = make_nstate_update_data()

        # Mask out all columns and patches
        d.mask_soilc .= false
        d.mask_soilp .= false

        initial_npool = copy(d.ns_veg.npool_patch)
        initial_leafn = copy(d.ns_veg.leafn_patch)

        CLM.n_state_update1!(d.ns_veg, d.nf_veg, d.nf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            ivt=d.ivt,
            woody=d.woody,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        # Nothing should change when all masks are false
        @test d.ns_veg.npool_patch == initial_npool
        @test d.ns_veg.leafn_patch == initial_leafn
    end

    # ================================================================
    # Test storage-to-xfer transfer for woody patches
    # ================================================================
    @testset "n_state_update1! storage to xfer" begin
        d = make_nstate_update_data()

        # Make patches woody
        d.ivt .= 2

        # Set storage-to-xfer fluxes
        d.nf_veg.leafn_storage_to_xfer_patch      .= 1.0e-7
        d.nf_veg.frootn_storage_to_xfer_patch     .= 0.8e-7
        d.nf_veg.livestemn_storage_to_xfer_patch  .= 0.5e-7
        d.nf_veg.deadstemn_storage_to_xfer_patch  .= 0.3e-7
        d.nf_veg.livecrootn_storage_to_xfer_patch .= 0.2e-7
        d.nf_veg.deadcrootn_storage_to_xfer_patch .= 0.1e-7

        initial_leafn_storage = copy(d.ns_veg.leafn_storage_patch)
        initial_leafn_xfer   = copy(d.ns_veg.leafn_xfer_patch)

        CLM.n_state_update1!(d.ns_veg, d.nf_veg, d.nf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            ivt=d.ivt,
            woody=d.woody,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        for p in 1:d.np
            # leafn_storage should decrease by storage_to_xfer and increase from alloc
            alloc_delta = d.nf_veg.npool_to_leafn_storage_patch[p] * d.dt
            xfer_delta  = d.nf_veg.leafn_storage_to_xfer_patch[p] * d.dt
            expected_storage = initial_leafn_storage[p] + alloc_delta - xfer_delta
            @test d.ns_veg.leafn_storage_patch[p] ≈ expected_storage

            # leafn_xfer should decrease by xfer_to_leafn and increase from storage_to_xfer
            xfer_to_leaf = d.nf_veg.leafn_xfer_to_leafn_patch[p] * d.dt
            stor_to_xfer = d.nf_veg.leafn_storage_to_xfer_patch[p] * d.dt
            expected_xfer = initial_leafn_xfer[p] - xfer_to_leaf + stor_to_xfer
            @test d.ns_veg.leafn_xfer_patch[p] ≈ expected_xfer
        end
    end

    # ================================================================
    # Test retranslocation and npool accounting
    # ================================================================
    @testset "n_state_update1! retransn and npool" begin
        d = make_nstate_update_data()

        initial_retransn = copy(d.ns_veg.retransn_patch)
        initial_npool    = copy(d.ns_veg.npool_patch)

        CLM.n_state_update1!(d.ns_veg, d.nf_veg, d.nf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            ivt=d.ivt,
            woody=d.woody,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        for p in 1:d.np
            # retransn: gains from leafn_to_retransn, loses to retransn_to_npool
            retrans_gain = d.nf_veg.leafn_to_retransn_patch[p] * d.dt
            retrans_loss = d.nf_veg.retransn_to_npool_patch[p] * d.dt
            expected = initial_retransn[p] + retrans_gain - retrans_loss
            @test d.ns_veg.retransn_patch[p] ≈ expected
        end
    end

end
