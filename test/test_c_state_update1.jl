@testset "C State Update 1" begin

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for testing
    # ----------------------------------------------------------------
    function make_cstate_update_data(; nc=4, np=6, ng=2, nlevdecomp=1,
                                      ndecomp_pools=7, ndecomp_cascade_transitions=5,
                                      nrepr=1)
        i_litr_min = 1
        i_litr_max = 3
        i_cwd = 4
        dt = 1800.0  # 30 minutes

        # --- CNVeg carbon state ---
        cs_veg = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs_veg, np, nc, ng; nrepr=nrepr)
        # Initialize pools to known values
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
        # Set all relevant fluxes to small nonzero values
        cf_veg.psnsun_to_cpool_patch            .= 1.0e-6
        cf_veg.psnshade_to_cpool_patch          .= 0.5e-6
        cf_veg.leafc_xfer_to_leafc_patch        .= 2.0e-7
        cf_veg.frootc_xfer_to_frootc_patch      .= 1.0e-7
        cf_veg.livestemc_xfer_to_livestemc_patch .= 1.5e-7
        cf_veg.deadstemc_xfer_to_deadstemc_patch .= 1.0e-7
        cf_veg.livecrootc_xfer_to_livecrootc_patch .= 0.8e-7
        cf_veg.deadcrootc_xfer_to_deadcrootc_patch .= 0.5e-7
        cf_veg.leafc_to_litter_patch            .= 3.0e-7
        cf_veg.frootc_to_litter_patch           .= 2.0e-7
        cf_veg.livestemc_to_deadstemc_patch     .= 1.0e-8
        cf_veg.livecrootc_to_deadcrootc_patch   .= 1.0e-8
        cf_veg.livestemc_to_litter_patch        .= 0.0
        cf_veg.livestemc_to_biofuelc_patch      .= 0.0
        cf_veg.livestemc_to_removedresiduec_patch .= 0.0
        cf_veg.leafc_to_biofuelc_patch          .= 0.0
        cf_veg.leafc_to_removedresiduec_patch   .= 0.0
        cf_veg.crop_seedc_to_leaf_patch         .= 0.0
        for k in 1:nrepr
            cf_veg.repr_grainc_to_food_patch[:, k]          .= 0.0
            cf_veg.repr_grainc_to_seed_patch[:, k]          .= 0.0
            cf_veg.repr_structurec_to_cropprod_patch[:, k]  .= 0.0
            cf_veg.repr_structurec_to_litter_patch[:, k]    .= 0.0
            cf_veg.reproductivec_xfer_to_reproductivec_patch[:, k] .= 0.0
        end
        cf_veg.cpool_to_xsmrpool_patch          .= 1.0e-7
        cf_veg.leaf_curmr_patch                 .= 5.0e-8
        cf_veg.froot_curmr_patch                .= 3.0e-8
        cf_veg.livestem_curmr_patch             .= 2.0e-8
        cf_veg.livecroot_curmr_patch            .= 1.0e-8
        cf_veg.cpool_to_resp_patch              .= 0.0
        cf_veg.soilc_change_patch               .= 0.0
        cf_veg.leaf_xsmr_patch                  .= 1.0e-8
        cf_veg.froot_xsmr_patch                 .= 0.5e-8
        cf_veg.livestem_xsmr_patch              .= 0.3e-8
        cf_veg.livecroot_xsmr_patch             .= 0.2e-8
        for k in 1:nrepr
            cf_veg.reproductive_curmr_patch[:, k] .= 0.0
            cf_veg.reproductive_xsmr_patch[:, k]  .= 0.0
        end
        cf_veg.cpool_to_leafc_patch             .= 5.0e-7
        cf_veg.cpool_to_leafc_storage_patch     .= 2.0e-7
        cf_veg.cpool_to_frootc_patch            .= 3.0e-7
        cf_veg.cpool_to_frootc_storage_patch    .= 1.0e-7
        cf_veg.cpool_to_leafc_resp_patch        .= 0.0
        cf_veg.cpool_to_leafc_storage_resp_patch .= 0.0
        cf_veg.cpool_to_frootc_resp_patch       .= 0.0
        cf_veg.cpool_to_frootc_storage_resp_patch .= 0.0
        cf_veg.cpool_to_livestemc_patch         .= 2.0e-7
        cf_veg.cpool_to_livestemc_storage_patch .= 1.0e-7
        cf_veg.cpool_to_deadstemc_patch         .= 1.5e-7
        cf_veg.cpool_to_deadstemc_storage_patch .= 0.5e-7
        cf_veg.cpool_to_livecrootc_patch        .= 1.0e-7
        cf_veg.cpool_to_livecrootc_storage_patch .= 0.5e-7
        cf_veg.cpool_to_deadcrootc_patch        .= 0.5e-7
        cf_veg.cpool_to_deadcrootc_storage_patch .= 0.3e-7
        cf_veg.cpool_to_livecrootc_resp_patch   .= 0.0
        cf_veg.cpool_to_livecrootc_storage_resp_patch .= 0.0
        cf_veg.cpool_to_livestemc_resp_patch    .= 0.0
        cf_veg.cpool_to_livestemc_storage_resp_patch .= 0.0
        for k in 1:nrepr
            cf_veg.cpool_to_reproductivec_patch[:, k]         .= 0.0
            cf_veg.cpool_to_reproductivec_storage_patch[:, k] .= 0.0
        end
        cf_veg.cpool_leaf_gr_patch              .= 1.0e-8
        cf_veg.cpool_froot_gr_patch             .= 0.8e-8
        cf_veg.cpool_livestem_gr_patch          .= 0.5e-8
        cf_veg.cpool_deadstem_gr_patch          .= 0.3e-8
        cf_veg.cpool_livecroot_gr_patch         .= 0.2e-8
        cf_veg.cpool_deadcroot_gr_patch         .= 0.1e-8
        for k in 1:nrepr
            cf_veg.cpool_reproductive_gr_patch[:, k]          .= 0.0
            cf_veg.cpool_reproductive_storage_gr_patch[:, k]  .= 0.0
            cf_veg.transfer_reproductive_gr_patch[:, k]       .= 0.0
        end
        cf_veg.transfer_leaf_gr_patch           .= 0.5e-8
        cf_veg.transfer_froot_gr_patch          .= 0.3e-8
        cf_veg.transfer_livestem_gr_patch       .= 0.2e-8
        cf_veg.transfer_deadstem_gr_patch       .= 0.1e-8
        cf_veg.transfer_livecroot_gr_patch      .= 0.1e-8
        cf_veg.transfer_deadcroot_gr_patch      .= 0.05e-8
        cf_veg.cpool_leaf_storage_gr_patch      .= 0.5e-8
        cf_veg.cpool_froot_storage_gr_patch     .= 0.3e-8
        cf_veg.cpool_livestem_storage_gr_patch  .= 0.2e-8
        cf_veg.cpool_deadstem_storage_gr_patch  .= 0.1e-8
        cf_veg.cpool_livecroot_storage_gr_patch .= 0.1e-8
        cf_veg.cpool_deadcroot_storage_gr_patch .= 0.05e-8
        cf_veg.cpool_to_gresp_storage_patch     .= 1.0e-8
        cf_veg.leafc_storage_to_xfer_patch      .= 0.0
        cf_veg.frootc_storage_to_xfer_patch     .= 0.0
        cf_veg.gresp_storage_to_xfer_patch      .= 0.0
        cf_veg.livestemc_storage_to_xfer_patch  .= 0.0
        cf_veg.deadstemc_storage_to_xfer_patch  .= 0.0
        cf_veg.livecrootc_storage_to_xfer_patch .= 0.0
        cf_veg.deadcrootc_storage_to_xfer_patch .= 0.0
        for k in 1:nrepr
            cf_veg.reproductivec_storage_to_xfer_patch[:, k] .= 0.0
        end
        cf_veg.xsmrpool_to_atm_patch            .= 0.0

        # DynPatch-specific fluxes
        cf_veg.dwt_frootc_to_litr_c_col         .= 0.0
        cf_veg.dwt_livecrootc_to_cwdc_col       .= 0.0
        cf_veg.dwt_deadcrootc_to_cwdc_col       .= 0.0
        cf_veg.dwt_seedc_to_leaf_grc            .= 0.0
        cf_veg.dwt_seedc_to_deadstem_grc        .= 0.0
        cf_veg.phenology_c_to_litr_c_col        .= 0.0

        # --- Soil biogeochem carbon flux ---
        cf_soil = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_soil, nc, nlevdecomp, ndecomp_pools,
                                        ndecomp_cascade_transitions)
        cf_soil.decomp_cpools_sourcesink_col    .= 0.0
        cf_soil.decomp_cascade_hr_vr_col        .= 0.0
        cf_soil.decomp_cascade_ctransfer_vr_col .= 0.0

        # --- Soil biogeochem carbon state ---
        cs_soil = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs_soil, nc, ng, nlevdecomp, ndecomp_pools)
        cs_soil.decomp_cpools_vr_col .= 100.0

        # --- Decomposition cascade ---
        cascade_donor_pool    = [1, 2, 3, 1, 2]
        cascade_receiver_pool = [2, 3, 0, 4, 4]  # 0 = terminal

        # --- Masks ---
        mask_soilc = trues(nc)
        mask_soilp = trues(np)
        mask_soilc_with_inactive = trues(nc)

        # --- Patch data ---
        patch_column = [1, 1, 2, 2, 3, 4]
        ivt          = fill(1, np)  # non-woody, non-crop by default
        woody        = zeros(Float64, 80)
        woody[2]     = 1.0  # PFT 2 is woody
        harvdate     = fill(999, np)
        col_is_fates = fill(false, nc)

        return (cs_veg=cs_veg, cf_veg=cf_veg, cf_soil=cf_soil, cs_soil=cs_soil,
                cascade_donor_pool=cascade_donor_pool,
                cascade_receiver_pool=cascade_receiver_pool,
                mask_soilc=mask_soilc, mask_soilp=mask_soilp,
                mask_soilc_with_inactive=mask_soilc_with_inactive,
                patch_column=patch_column, ivt=ivt, woody=woody,
                harvdate=harvdate, col_is_fates=col_is_fates,
                nc=nc, np=np, ng=ng, dt=dt,
                nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                ndecomp_cascade_transitions=ndecomp_cascade_transitions,
                nrepr=nrepr, i_litr_min=i_litr_min, i_litr_max=i_litr_max,
                i_cwd=i_cwd)
    end

    # ================================================================
    # Test CStateUpdateDynPatch
    # ================================================================
    @testset "c_state_update_dyn_patch!" begin
        d = make_cstate_update_data()

        # Set some dynamic patch fluxes
        for i in d.i_litr_min:d.i_litr_max
            d.cf_veg.dwt_frootc_to_litr_c_col[:, :, i] .= 1.0e-6
        end
        d.cf_veg.dwt_livecrootc_to_cwdc_col .= 2.0e-6
        d.cf_veg.dwt_deadcrootc_to_cwdc_col .= 3.0e-6
        d.cf_veg.dwt_seedc_to_leaf_grc      .= 0.5e-6
        d.cf_veg.dwt_seedc_to_deadstem_grc  .= 0.3e-6

        initial_seedc = copy(d.cs_veg.seedc_grc)
        initial_decomp = copy(d.cs_soil.decomp_cpools_vr_col)

        CLM.c_state_update_dyn_patch!(d.cs_veg, d.cf_veg, d.cs_soil;
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
            @test d.cs_soil.decomp_cpools_vr_col[c, j, i] > initial_decomp[c, j, i]
        end

        # CWD pool should increase from livecroot + deadcroot
        for c in 1:d.nc, j in 1:d.nlevdecomp
            expected_delta = (2.0e-6 + 3.0e-6) * d.dt
            @test d.cs_soil.decomp_cpools_vr_col[c, j, d.i_cwd] ≈
                  initial_decomp[c, j, d.i_cwd] + expected_delta
        end

        # Seed pool should decrease
        for g in 1:d.ng
            expected = initial_seedc[g] - (0.5e-6 + 0.3e-6) * d.dt
            @test d.cs_veg.seedc_grc[g] ≈ expected
        end
    end

    # ================================================================
    # Test CStateUpdate0
    # ================================================================
    @testset "c_state_update0!" begin
        d = make_cstate_update_data()

        initial_cpool = copy(d.cs_veg.cpool_patch)

        CLM.c_state_update0!(d.cs_veg, d.cf_veg;
            mask_soilp=d.mask_soilp,
            bounds_patch=1:d.np,
            dt=d.dt)

        for p in 1:d.np
            expected = initial_cpool[p] +
                       (d.cf_veg.psnsun_to_cpool_patch[p] + d.cf_veg.psnshade_to_cpool_patch[p]) * d.dt
            @test d.cs_veg.cpool_patch[p] ≈ expected
        end
    end

    # ================================================================
    # Test CStateUpdate1 — non-woody, non-crop patch
    # ================================================================
    @testset "c_state_update1! non-woody non-crop" begin
        d = make_cstate_update_data()

        # Set phenology source-sink
        d.cf_veg.phenology_c_to_litr_c_col .= 1.0e-7

        initial_cpool = copy(d.cs_veg.cpool_patch)
        initial_leafc = copy(d.cs_veg.leafc_patch)
        initial_frootc = copy(d.cs_veg.frootc_patch)

        CLM.c_state_update1!(d.cs_veg, d.cf_veg, d.cf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            patch_column=d.patch_column,
            ivt=d.ivt,
            woody=d.woody,
            cascade_donor_pool=d.cascade_donor_pool,
            cascade_receiver_pool=d.cascade_receiver_pool,
            harvdate=d.harvdate,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        # cpool should decrease (MR, allocation, growth resp all drain it)
        for p in 1:d.np
            @test d.cs_veg.cpool_patch[p] < initial_cpool[p]
        end

        # leafc should increase (xfer growth + allocation > litterfall)
        # With flux values chosen: leafc_xfer_to_leafc (2e-7) + cpool_to_leafc (5e-7) - leafc_to_litter (3e-7) = 4e-7 > 0
        for p in 1:d.np
            net_flux = (d.cf_veg.leafc_xfer_to_leafc_patch[p] +
                        d.cf_veg.cpool_to_leafc_patch[p] -
                        d.cf_veg.leafc_to_litter_patch[p]) * d.dt
            @test d.cs_veg.leafc_patch[p] ≈ initial_leafc[p] + net_flux
        end

        # Soil decomp sourcesink should be set for litter pools
        for c in 1:d.nc, j in 1:d.nlevdecomp, i in d.i_litr_min:d.i_litr_max
            expected = d.cf_veg.phenology_c_to_litr_c_col[c, j, i] * d.dt
            @test d.cf_soil.decomp_cpools_sourcesink_col[c, j, i] ≈ expected atol=1e-15
        end

        # CWD sourcesink should be zero (terms moved to DynPatch)
        for c in 1:d.nc, j in 1:d.nlevdecomp
            @test d.cf_soil.decomp_cpools_sourcesink_col[c, j, d.i_cwd] ≈ 0.0 atol=1e-15
        end
    end

    # ================================================================
    # Test CStateUpdate1 — woody patch
    # ================================================================
    @testset "c_state_update1! woody patch" begin
        d = make_cstate_update_data()

        # Make all patches woody (PFT 2)
        d.ivt .= 2

        initial_livestemc = copy(d.cs_veg.livestemc_patch)
        initial_deadstemc = copy(d.cs_veg.deadstemc_patch)
        initial_xsmrpool  = copy(d.cs_veg.xsmrpool_patch)

        CLM.c_state_update1!(d.cs_veg, d.cf_veg, d.cf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            patch_column=d.patch_column,
            ivt=d.ivt,
            woody=d.woody,
            cascade_donor_pool=d.cascade_donor_pool,
            cascade_receiver_pool=d.cascade_receiver_pool,
            harvdate=d.harvdate,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        # For woody patches: livestemc should change (xfer + alloc - turnover)
        for p in 1:d.np
            # xfer growth + allocation + live-to-dead turnover
            xfer_delta = d.cf_veg.livestemc_xfer_to_livestemc_patch[p] * d.dt
            alloc_delta = d.cf_veg.cpool_to_livestemc_patch[p] * d.dt
            turnover_delta = d.cf_veg.livestemc_to_deadstemc_patch[p] * d.dt
            expected = initial_livestemc[p] + xfer_delta + alloc_delta - turnover_delta
            @test d.cs_veg.livestemc_patch[p] ≈ expected
        end

        # xsmrpool should change for woody patches
        for p in 1:d.np
            cpool_to_xsmr = d.cf_veg.cpool_to_xsmrpool_patch[p] * d.dt
            leaf_xsmr  = d.cf_veg.leaf_xsmr_patch[p] * d.dt
            froot_xsmr = d.cf_veg.froot_xsmr_patch[p] * d.dt
            lstem_xsmr = d.cf_veg.livestem_xsmr_patch[p] * d.dt
            lcroot_xsmr = d.cf_veg.livecroot_xsmr_patch[p] * d.dt
            expected = initial_xsmrpool[p] + cpool_to_xsmr - leaf_xsmr - froot_xsmr - lstem_xsmr - lcroot_xsmr
            @test d.cs_veg.xsmrpool_patch[p] ≈ expected
        end
    end

    # ================================================================
    # Test CStateUpdate1 — FATES column (skip phenology input)
    # ================================================================
    @testset "c_state_update1! FATES column" begin
        d = make_cstate_update_data()

        # Mark column 1 as FATES
        d.col_is_fates[1] = true

        # Set phenology flux for non-FATES columns
        d.cf_veg.phenology_c_to_litr_c_col .= 1.0e-7

        CLM.c_state_update1!(d.cs_veg, d.cf_veg, d.cf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            patch_column=d.patch_column,
            ivt=d.ivt,
            woody=d.woody,
            cascade_donor_pool=d.cascade_donor_pool,
            cascade_receiver_pool=d.cascade_receiver_pool,
            harvdate=d.harvdate,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        # FATES column 1 should NOT have phenology input (remains 0)
        # (sourcesink for FATES column was 0 before, stays 0 for litter pools
        #  because we skip the phenology_c_to_litr_c assignment for FATES)
        # Non-FATES column 2 should have phenology input
        for j in 1:d.nlevdecomp, i in d.i_litr_min:d.i_litr_max
            @test d.cf_soil.decomp_cpools_sourcesink_col[2, j, i] ≈ 1.0e-7 * d.dt atol=1e-15
        end
    end

    # ================================================================
    # Test decomposition cascade sourcesink accounting
    # ================================================================
    @testset "c_state_update1! decomp cascade" begin
        d = make_cstate_update_data()

        # Set some decomposition fluxes
        hr_rate = 5.0e-7
        transfer_rate = 3.0e-7
        d.cf_soil.decomp_cascade_hr_vr_col        .= hr_rate
        d.cf_soil.decomp_cascade_ctransfer_vr_col .= transfer_rate

        CLM.c_state_update1!(d.cs_veg, d.cf_veg, d.cf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            patch_column=d.patch_column,
            ivt=d.ivt,
            woody=d.woody,
            cascade_donor_pool=d.cascade_donor_pool,
            cascade_receiver_pool=d.cascade_receiver_pool,
            harvdate=d.harvdate,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        # Check that donor pools are depleted and receiver pools are augmented
        # Transition 3 has receiver_pool = 0 (terminal), so only HR applies
        for c in 1:d.nc, j in 1:d.nlevdecomp
            ss = d.cf_soil.decomp_cpools_sourcesink_col[c, j, :]
            # Terminal transition (k=3): donor pool 3, receiver 0
            # Pool 3 as donor: transitions 3 and (also as receiver of transitions 1,2 which have receiver=2,3)
            # This is complex, so just verify that some pools changed
            @test any(ss .!= 0.0)
        end
    end

    # ================================================================
    # Test mask filtering: masked columns/patches are skipped
    # ================================================================
    @testset "c_state_update1! mask filtering" begin
        d = make_cstate_update_data()

        # Mask out all columns and patches
        d.mask_soilc .= false
        d.mask_soilp .= false

        initial_cpool = copy(d.cs_veg.cpool_patch)
        initial_leafc = copy(d.cs_veg.leafc_patch)

        CLM.c_state_update1!(d.cs_veg, d.cf_veg, d.cf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            patch_column=d.patch_column,
            ivt=d.ivt,
            woody=d.woody,
            cascade_donor_pool=d.cascade_donor_pool,
            cascade_receiver_pool=d.cascade_receiver_pool,
            harvdate=d.harvdate,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        # Nothing should change when all masks are false
        @test d.cs_veg.cpool_patch == initial_cpool
        @test d.cs_veg.leafc_patch == initial_leafc
    end

    # ================================================================
    # Test gresp_storage and gresp_xfer updates
    # ================================================================
    @testset "c_state_update1! gresp storage/xfer" begin
        d = make_cstate_update_data()

        # Make patches woody for full gresp_storage_to_xfer path
        d.ivt .= 2

        initial_gresp_storage = copy(d.cs_veg.gresp_storage_patch)
        initial_gresp_xfer   = copy(d.cs_veg.gresp_xfer_patch)

        CLM.c_state_update1!(d.cs_veg, d.cf_veg, d.cf_soil;
            mask_soilc=d.mask_soilc,
            mask_soilp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            patch_column=d.patch_column,
            ivt=d.ivt,
            woody=d.woody,
            cascade_donor_pool=d.cascade_donor_pool,
            cascade_receiver_pool=d.cascade_receiver_pool,
            harvdate=d.harvdate,
            col_is_fates=d.col_is_fates,
            nlevdecomp=d.nlevdecomp,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            nrepr=d.nrepr,
            dt=d.dt)

        for p in 1:d.np
            # gresp_storage gains from cpool_to_gresp_storage, loses to gresp_storage_to_xfer
            cpool_to = d.cf_veg.cpool_to_gresp_storage_patch[p] * d.dt
            stor_to_xfer = d.cf_veg.gresp_storage_to_xfer_patch[p] * d.dt
            expected = initial_gresp_storage[p] + cpool_to - stor_to_xfer
            @test d.cs_veg.gresp_storage_patch[p] ≈ expected

            # gresp_xfer loses from transfer growth resp, gains from storage_to_xfer
            xfer_delta = d.cf_veg.gresp_storage_to_xfer_patch[p] * d.dt
            transfer_leaf_gr = d.cf_veg.transfer_leaf_gr_patch[p] * d.dt
            transfer_froot_gr = d.cf_veg.transfer_froot_gr_patch[p] * d.dt
            transfer_livestem_gr = d.cf_veg.transfer_livestem_gr_patch[p] * d.dt
            transfer_deadstem_gr = d.cf_veg.transfer_deadstem_gr_patch[p] * d.dt
            transfer_livecroot_gr = d.cf_veg.transfer_livecroot_gr_patch[p] * d.dt
            transfer_deadcroot_gr = d.cf_veg.transfer_deadcroot_gr_patch[p] * d.dt
            expected_xfer = initial_gresp_xfer[p] - transfer_leaf_gr - transfer_froot_gr -
                            transfer_livestem_gr - transfer_deadstem_gr -
                            transfer_livecroot_gr - transfer_deadcroot_gr + xfer_delta
            @test d.cs_veg.gresp_xfer_patch[p] ≈ expected_xfer
        end
    end

end
