@testset "CN Products Module" begin

    # ====================================================================
    # Helper: build minimal subgrid data for ng gridcells, np patches
    # ====================================================================
    function make_subgrid(; ng=2, np=4)
        # Map patches evenly to gridcells
        # For ng=1: all patches -> gridcell 1
        # For ng=2: first half -> g=1, second half -> g=2
        patch_gc = Int[]
        if ng == 1
            patch_gc = fill(1, np)
        else
            half = div(np, 2)
            patch_gc = vcat(fill(1, half), fill(2, np - half))
        end
        patches_per_gc = div(np, ng)

        # PFT types: use types 0..np-1 (0-based Fortran indices)
        patch_ivt = collect(0:np-1)

        # Build a PatchData with the minimum required fields for p2g_1d!
        pch = CLM.PatchData()
        pch.gridcell = patch_gc
        pch.itype    = patch_ivt
        pch.active   = fill(true, np)
        pch.wtgcell  = fill(1.0 / patches_per_gc, np)  # equal weight within each gridcell
        pch.wtcol    = fill(1.0, np)
        pch.wtlunit  = fill(1.0, np)
        pch.column   = collect(1:np)
        pch.landunit = fill(1, np)             # all patches in landunit 1

        # LandunitData: single landunit of type ISTSOIL
        lun = CLM.LandunitData()
        lun.itype   = [CLM.ISTSOIL]
        lun.active  = [true]
        lun.wtgcell = [1.0]
        lun.gridcell = [1]
        lun.urbpoi  = [false]
        lun.canyon_hwr = [0.0]

        # ColumnData: one column per patch (simplest mapping)
        col = CLM.ColumnData()
        col.active   = fill(true, np)
        col.wtgcell  = fill(1.0, np)
        col.wtlunit  = fill(1.0, np)
        col.landunit = fill(1, np)
        col.gridcell = patch_gc
        col.itype    = fill(1, np)

        # Bounds
        bounds = CLM.BoundsType(
            begg = 1, endg = ng,
            begl = 1, endl = 1,
            begc = 1, endc = np,
            begp = 1, endp = np,
        )

        return pch, col, lun, bounds
    end

    # ====================================================================
    # Helper: build PFT parameters for np PFT types
    # ====================================================================
    function make_pft_params(np)
        pprod10     = fill(0.2, np)
        pprod100    = fill(0.3, np)
        pprodharv10 = fill(0.4, np)
        return pprod10, pprod100, pprodharv10
    end

    # ====================================================================
    # Test: cn_products_full_init!
    # ====================================================================
    @testset "cn_products_full_init!" begin
        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, 3, 5)

        @test length(prod.cropprod1_grc) == 3
        @test length(prod.prod10_grc) == 3
        @test length(prod.prod100_grc) == 3
        @test length(prod.tot_woodprod_grc) == 3
        @test length(prod.product_loss_grc) == 3
        @test length(prod.dwt_prod10_gain_grc) == 3
        @test length(prod.gru_prod10_gain_patch) == 5
        @test length(prod.hrv_deadstem_to_prod10_patch) == 5
        @test length(prod.crop_harvest_to_cropprod1_patch) == 5

        # All should be zero
        @test all(prod.cropprod1_grc .== 0.0)
        @test all(prod.prod10_grc .== 0.0)
        @test all(prod.prod100_grc .== 0.0)
        @test all(prod.gru_prod10_gain_patch .== 0.0)
        @test all(prod.cropprod1_loss_grc .== 0.0)
    end

    # ====================================================================
    # Test: cn_products_set_values!
    # ====================================================================
    @testset "cn_products_set_values!" begin
        ng = 2; np = 4
        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, ng, np)

        # Set some non-zero values
        prod.dwt_prod10_gain_grc .= 99.0
        prod.dwt_prod100_gain_grc .= 99.0
        prod.dwt_cropprod1_gain_grc .= 99.0
        prod.crop_harvest_to_cropprod1_grc .= 99.0
        prod.hrv_deadstem_to_prod10_grc .= 99.0
        prod.hrv_deadstem_to_prod100_grc .= 99.0

        bounds = CLM.BoundsType(begg=1, endg=ng, begp=1, endp=np)

        CLM.cn_products_set_values!(prod, bounds, 0.0)

        @test all(prod.dwt_prod10_gain_grc .== 0.0)
        @test all(prod.dwt_prod100_gain_grc .== 0.0)
        @test all(prod.dwt_cropprod1_gain_grc .== 0.0)
        @test all(prod.crop_harvest_to_cropprod1_grc .== 0.0)
        @test all(prod.hrv_deadstem_to_prod10_grc .== 0.0)
        @test all(prod.hrv_deadstem_to_prod100_grc .== 0.0)
    end

    # ====================================================================
    # Test: cn_products_compute_product_summary! (decay + gain update)
    # ====================================================================
    @testset "cn_products_compute_product_summary!" begin
        ng = 2; np = 4
        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, ng, np)
        bounds = CLM.BoundsType(begg=1, endg=ng, begp=1, endp=np)

        # Set initial pool states
        prod.cropprod1_grc[1] = 100.0
        prod.prod10_grc[1]    = 200.0
        prod.prod100_grc[1]   = 300.0
        prod.cropprod1_grc[2] = 50.0
        prod.prod10_grc[2]    = 100.0
        prod.prod100_grc[2]   = 150.0

        # Set some gain fluxes (g/m2/s)
        prod.dwt_cropprod1_gain_grc[1]        = 1.0e-5
        prod.dwt_prod10_gain_grc[1]           = 2.0e-5
        prod.dwt_prod100_gain_grc[1]          = 3.0e-5
        prod.gru_prod10_gain_grc[1]           = 0.5e-5
        prod.gru_prod100_gain_grc[1]          = 0.8e-5
        prod.crop_harvest_to_cropprod1_grc[1] = 0.2e-5
        prod.hrv_deadstem_to_prod10_grc[1]    = 0.3e-5
        prod.hrv_deadstem_to_prod100_grc[1]   = 0.4e-5

        dt = 1800.0  # 30-minute timestep

        # Expected loss rates before update
        kprod1   = 7.2e-8
        kprod10  = 7.2e-9
        kprod100 = 7.2e-10

        expected_crop1_loss_g1 = 100.0 * kprod1
        expected_p10_loss_g1   = 200.0 * kprod10
        expected_p100_loss_g1  = 300.0 * kprod100

        # Run the routine
        CLM.cn_products_compute_product_summary!(prod, bounds, dt)

        # Check loss fluxes were computed
        @test prod.cropprod1_loss_grc[1] ≈ expected_crop1_loss_g1
        @test prod.prod10_loss_grc[1]    ≈ expected_p10_loss_g1
        @test prod.prod100_loss_grc[1]   ≈ expected_p100_loss_g1

        # Check updated pool states for gridcell 1
        expected_crop1 = 100.0 +
            1.0e-5 * dt +     # dwt gain
            0.2e-5 * dt -     # harvest gain
            expected_crop1_loss_g1 * dt  # loss

        expected_p10 = 200.0 +
            2.0e-5 * dt +     # dwt gain
            0.5e-5 * dt +     # gru gain
            0.3e-5 * dt -     # harvest gain
            expected_p10_loss_g1 * dt    # loss

        expected_p100 = 300.0 +
            3.0e-5 * dt +     # dwt gain
            0.8e-5 * dt +     # gru gain
            0.4e-5 * dt -     # harvest gain
            expected_p100_loss_g1 * dt   # loss

        @test prod.cropprod1_grc[1] ≈ expected_crop1
        @test prod.prod10_grc[1]    ≈ expected_p10
        @test prod.prod100_grc[1]   ≈ expected_p100

        # Check gridcell 2 (no gain fluxes, only losses)
        expected_crop1_loss_g2 = 50.0 * kprod1
        expected_p10_loss_g2   = 100.0 * kprod10
        expected_p100_loss_g2  = 150.0 * kprod100

        @test prod.cropprod1_grc[2] ≈ 50.0 - expected_crop1_loss_g2 * dt
        @test prod.prod10_grc[2]    ≈ 100.0 - expected_p10_loss_g2 * dt
        @test prod.prod100_grc[2]   ≈ 150.0 - expected_p100_loss_g2 * dt
    end

    # ====================================================================
    # Test: cn_products_compute_summary! (summary variables)
    # ====================================================================
    @testset "cn_products_compute_summary!" begin
        ng = 2; np = 4
        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, ng, np)
        bounds = CLM.BoundsType(begg=1, endg=ng, begp=1, endp=np)

        prod.prod10_grc[1]  = 100.0
        prod.prod100_grc[1] = 200.0

        prod.cropprod1_loss_grc[1] = 0.01
        prod.prod10_loss_grc[1]    = 0.02
        prod.prod100_loss_grc[1]   = 0.03

        prod.dwt_prod10_gain_grc[1]  = 1.0e-5
        prod.dwt_prod100_gain_grc[1] = 2.0e-5
        prod.gru_prod10_gain_grc[1]  = 3.0e-5
        prod.gru_prod100_gain_grc[1] = 4.0e-5

        CLM.cn_products_compute_summary!(prod, bounds)

        # tot_woodprod = prod10 + prod100
        @test prod.tot_woodprod_grc[1] ≈ 300.0

        # tot_woodprod_loss = prod10_loss + prod100_loss
        @test prod.tot_woodprod_loss_grc[1] ≈ 0.05

        # product_loss = cropprod1_loss + prod10_loss + prod100_loss
        @test prod.product_loss_grc[1] ≈ 0.06

        # dwt_woodprod_gain = dwt_prod10_gain + dwt_prod100_gain
        @test prod.dwt_woodprod_gain_grc[1] ≈ 3.0e-5

        # gru_woodprod_gain = gru_prod10_gain + gru_prod100_gain
        @test prod.gru_woodprod_gain_grc[1] ≈ 7.0e-5

        # Gridcell 2 should be all zeros
        @test prod.tot_woodprod_grc[2] ≈ 0.0
        @test prod.tot_woodprod_loss_grc[2] ≈ 0.0
        @test prod.product_loss_grc[2] ≈ 0.0
    end

    # ====================================================================
    # Test: cn_products_partition_wood_fluxes! (PFT partitioning logic)
    # ====================================================================
    @testset "cn_products_partition_wood_fluxes!" begin
        ng = 2; np = 4
        pch, col, lun, bounds = make_subgrid(ng=ng, np=np)
        pprod10, pprod100, pprodharv10 = make_pft_params(np)

        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, ng, np)

        # All patches are in the soil filter
        filter_soilp = collect(1:np)
        num_soilp = np

        # Input fluxes
        gru_wood  = [0.1, 0.2, 0.3, 0.4]     # g/m2/s
        wood_harv = [1.0, 2.0, 3.0, 4.0]      # g/m2/s
        dwt_wood  = zeros(np)                   # no dwt fluxes for this test

        CLM.cn_products_partition_wood_fluxes!(prod, bounds,
            num_soilp, filter_soilp,
            dwt_wood, gru_wood, wood_harv;
            pprod10 = pprod10,
            pprod100 = pprod100,
            pprodharv10 = pprodharv10,
            patch_gridcell = pch.gridcell,
            patch_itype = pch.itype,
            pch = pch, col = col, lun = lun)

        # Check gross unrepresented partitioning
        # pprod10/(pprod10+pprod100) = 0.2/0.5 = 0.4
        # pprod100/(pprod10+pprod100) = 0.3/0.5 = 0.6
        for p in 1:np
            @test prod.gru_prod10_gain_patch[p]  ≈ gru_wood[p] * 0.4
            @test prod.gru_prod100_gain_patch[p] ≈ gru_wood[p] * 0.6
        end

        # Check harvest partitioning
        # pprodharv10 = 0.4 for all PFTs
        for p in 1:np
            @test prod.hrv_deadstem_to_prod10_patch[p]  ≈ wood_harv[p] * 0.4
            @test prod.hrv_deadstem_to_prod100_patch[p] ≈ wood_harv[p] * 0.6
        end
    end

    # ====================================================================
    # Test: cn_products_partition_crop_fluxes!
    # ====================================================================
    @testset "cn_products_partition_crop_fluxes!" begin
        ng = 2; np = 4
        pch, col, lun, bounds = make_subgrid(ng=ng, np=np)

        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, ng, np)

        filter_soilp = collect(1:np)
        num_soilp = np

        dwt_crop_gain  = [0.001, 0.002, 0.003, 0.004]
        crop_harvest   = [0.01, 0.02, 0.03, 0.04]

        CLM.cn_products_partition_crop_fluxes!(prod, bounds,
            num_soilp, filter_soilp,
            dwt_crop_gain, crop_harvest;
            patch_gridcell = pch.gridcell,
            pch = pch, col = col, lun = lun)

        # crop_harvest_to_cropprod1_patch should match input directly
        for p in 1:np
            @test prod.crop_harvest_to_cropprod1_patch[p] ≈ crop_harvest[p]
        end

        # dwt_cropprod1_gain_grc: patches 1,2 -> g=1, patches 3,4 -> g=2
        @test prod.dwt_cropprod1_gain_grc[1] ≈ 0.001 + 0.002
        @test prod.dwt_cropprod1_gain_grc[2] ≈ 0.003 + 0.004
    end

    # ====================================================================
    # Test: Zero pprod_tot path (avoid divide-by-zero)
    # ====================================================================
    @testset "Zero pprod_tot path" begin
        ng = 1; np = 2
        pch, col, lun, bounds = make_subgrid(ng=ng, np=np)

        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, ng, np)

        filter_soilp = collect(1:np)
        num_soilp = np

        # Zero pprod10 and pprod100 => pprod_tot = 0
        pprod10     = [0.0, 0.0]
        pprod100    = [0.0, 0.0]
        pprodharv10 = [0.5, 0.5]

        gru_wood  = [1.0, 2.0]
        wood_harv = [0.5, 0.5]
        dwt_wood  = [0.0, 0.0]   # must be zero since pprod_tot = 0

        CLM.cn_products_partition_wood_fluxes!(prod, bounds,
            num_soilp, filter_soilp,
            dwt_wood, gru_wood, wood_harv;
            pprod10 = pprod10,
            pprod100 = pprod100,
            pprodharv10 = pprodharv10,
            patch_gridcell = pch.gridcell,
            patch_itype = pch.itype,
            pch = pch, col = col, lun = lun)

        # When pprod_tot = 0, fractions are 0 => no gru goes to product pools
        @test prod.gru_prod10_gain_patch[1]  ≈ 0.0
        @test prod.gru_prod100_gain_patch[1] ≈ 0.0
        @test prod.gru_prod10_gain_patch[2]  ≈ 0.0
        @test prod.gru_prod100_gain_patch[2] ≈ 0.0
    end

    # ====================================================================
    # Test: dwt_wood > 0 with pprod_tot == 0 raises error
    # ====================================================================
    @testset "dwt_wood > 0 with pprod_tot = 0 raises error" begin
        ng = 1; np = 2
        pch, col, lun, bounds = make_subgrid(ng=ng, np=np)

        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, ng, np)

        filter_soilp = collect(1:np)
        num_soilp = np

        pprod10     = [0.0, 0.0]
        pprod100    = [0.0, 0.0]
        pprodharv10 = [0.5, 0.5]

        gru_wood  = [0.0, 0.0]
        wood_harv = [0.0, 0.0]
        dwt_wood  = [0.5, 0.0]   # positive with pprod_tot == 0 -> error

        @test_throws ErrorException CLM.cn_products_partition_wood_fluxes!(
            prod, bounds,
            num_soilp, filter_soilp,
            dwt_wood, gru_wood, wood_harv;
            pprod10 = pprod10,
            pprod100 = pprod100,
            pprodharv10 = pprodharv10,
            patch_gridcell = pch.gridcell,
            patch_itype = pch.itype,
            pch = pch, col = col, lun = lun)
    end

    # ====================================================================
    # Test: Full integration (update + summary)
    # ====================================================================
    @testset "Full integration: update + product summary + summary" begin
        ng = 2; np = 4
        pch, col, lun, bounds = make_subgrid(ng=ng, np=np)
        pprod10, pprod100, pprodharv10 = make_pft_params(np)

        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, ng, np)

        # Set initial pool states
        prod.cropprod1_grc[1] = 50.0
        prod.prod10_grc[1]    = 100.0
        prod.prod100_grc[1]   = 200.0
        prod.cropprod1_grc[2] = 25.0
        prod.prod10_grc[2]    = 50.0
        prod.prod100_grc[2]   = 100.0

        filter_soilp = collect(1:np)
        num_soilp = np

        # Input fluxes
        dwt_wood_gain  = zeros(np)
        gru_wood_gain  = [0.1, 0.2, 0.3, 0.4]
        wood_harvest   = [1.0, 2.0, 3.0, 4.0]
        dwt_crop_gain  = zeros(np)
        crop_harvest   = [0.01, 0.02, 0.03, 0.04]

        dt = 1800.0

        # Step 1: zero gain flux accumulators
        CLM.cn_products_set_values!(prod, bounds, 0.0)

        # Step 2: partition fluxes
        CLM.cn_products_update!(prod, bounds,
            num_soilp, filter_soilp,
            dwt_wood_gain, gru_wood_gain, wood_harvest,
            dwt_crop_gain, crop_harvest;
            pprod10 = pprod10,
            pprod100 = pprod100,
            pprodharv10 = pprodharv10,
            patch_gridcell = pch.gridcell,
            patch_itype = pch.itype,
            pch = pch, col = col, lun = lun)

        # Step 3: compute pool updates
        CLM.cn_products_compute_product_summary!(prod, bounds, dt)

        # Step 4: compute summaries
        CLM.cn_products_compute_summary!(prod, bounds)

        # Verify summaries are consistent
        @test prod.tot_woodprod_grc[1] ≈ prod.prod10_grc[1] + prod.prod100_grc[1]
        @test prod.tot_woodprod_grc[2] ≈ prod.prod10_grc[2] + prod.prod100_grc[2]

        @test prod.tot_woodprod_loss_grc[1] ≈ prod.prod10_loss_grc[1] + prod.prod100_loss_grc[1]

        @test prod.product_loss_grc[1] ≈ (prod.cropprod1_loss_grc[1] +
                                           prod.prod10_loss_grc[1] +
                                           prod.prod100_loss_grc[1])

        @test prod.dwt_woodprod_gain_grc[1] ≈ prod.dwt_prod10_gain_grc[1] + prod.dwt_prod100_gain_grc[1]
        @test prod.gru_woodprod_gain_grc[1] ≈ prod.gru_prod10_gain_grc[1] + prod.gru_prod100_gain_grc[1]

        # Pool states should be positive (gains > losses for these inputs)
        @test prod.cropprod1_grc[1] > 0.0
        @test prod.prod10_grc[1] > 0.0
        @test prod.prod100_grc[1] > 0.0
        @test prod.cropprod1_grc[2] > 0.0
        @test prod.prod10_grc[2] > 0.0
        @test prod.prod100_grc[2] > 0.0
    end

    # ====================================================================
    # Test: Decay rate correctness (90% loss in pool-specific time)
    # ====================================================================
    @testset "Decay rates approximate 90% loss over pool lifespan" begin
        ng = 1; np = 1
        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, ng, np)
        bounds = CLM.BoundsType(begg=1, endg=ng, begp=1, endp=np)

        # Test 1-year pool: kprod1 = 7.2e-8 /s
        # Over 1 year: remaining = (1 - k*dt)^(N_steps)
        kprod1 = 7.2e-8
        dt = 1800.0                  # 30-min timestep
        secs_per_year = 365.25 * 86400.0
        nsteps = round(Int, secs_per_year / dt)

        remaining = (1.0 - kprod1 * dt) ^ nsteps
        # Should be approximately 0.1 (i.e., 90% lost)
        @test remaining < 0.15   # generous tolerance
        @test remaining > 0.05

        # Test 10-year pool
        kprod10 = 7.2e-9
        nsteps_10 = round(Int, 10.0 * secs_per_year / dt)
        remaining_10 = (1.0 - kprod10 * dt) ^ nsteps_10
        @test remaining_10 < 0.15
        @test remaining_10 > 0.05

        # Test 100-year pool
        kprod100 = 7.2e-10
        nsteps_100 = round(Int, 100.0 * secs_per_year / dt)
        remaining_100 = (1.0 - kprod100 * dt) ^ nsteps_100
        @test remaining_100 < 0.15
        @test remaining_100 > 0.05
    end

    # ====================================================================
    # Test: Multiple timesteps with no gains -> pure decay
    # ====================================================================
    @testset "Pure decay over multiple timesteps" begin
        ng = 1; np = 1
        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, ng, np)
        bounds = CLM.BoundsType(begg=1, endg=ng, begp=1, endp=np)

        initial_crop = 1000.0
        initial_p10  = 2000.0
        initial_p100 = 3000.0

        prod.cropprod1_grc[1] = initial_crop
        prod.prod10_grc[1]    = initial_p10
        prod.prod100_grc[1]   = initial_p100

        dt = 3600.0  # 1 hour
        nsteps = 100

        for _ in 1:nsteps
            CLM.cn_products_compute_product_summary!(prod, bounds, dt)
        end

        # All pools should have decreased
        @test prod.cropprod1_grc[1] < initial_crop
        @test prod.prod10_grc[1]    < initial_p10
        @test prod.prod100_grc[1]   < initial_p100

        # 1-year pool decays fastest
        crop_frac = prod.cropprod1_grc[1] / initial_crop
        p10_frac  = prod.prod10_grc[1]    / initial_p10
        p100_frac = prod.prod100_grc[1]   / initial_p100

        @test crop_frac < p10_frac < p100_frac

        # All should still be positive
        @test prod.cropprod1_grc[1] > 0.0
        @test prod.prod10_grc[1]    > 0.0
        @test prod.prod100_grc[1]   > 0.0
    end

    # ====================================================================
    # Test: Mass conservation -- gains - losses = change in pool
    # ====================================================================
    @testset "Mass conservation check" begin
        ng = 1; np = 2
        pch, col, lun, bounds = make_subgrid(ng=ng, np=np)
        pprod10, pprod100, pprodharv10 = make_pft_params(np)

        prod = CLM.CNProductsFullData()
        CLM.cn_products_full_init!(prod, ng, np)

        # Set initial states
        prod.cropprod1_grc[1] = 500.0
        prod.prod10_grc[1]    = 1000.0
        prod.prod100_grc[1]   = 2000.0

        filter_soilp = collect(1:np)
        num_soilp = np
        dt = 1800.0

        dwt_wood_gain  = zeros(np)
        gru_wood_gain  = [0.5, 0.5]
        wood_harvest   = [2.0, 2.0]
        dwt_crop_gain  = zeros(np)
        crop_harvest   = [0.05, 0.05]

        # Record state before
        crop1_before = prod.cropprod1_grc[1]
        p10_before   = prod.prod10_grc[1]
        p100_before  = prod.prod100_grc[1]

        CLM.cn_products_set_values!(prod, bounds, 0.0)

        CLM.cn_products_update!(prod, bounds,
            num_soilp, filter_soilp,
            dwt_wood_gain, gru_wood_gain, wood_harvest,
            dwt_crop_gain, crop_harvest;
            pprod10 = pprod10,
            pprod100 = pprod100,
            pprodharv10 = pprodharv10,
            patch_gridcell = pch.gridcell,
            patch_itype = pch.itype,
            pch = pch, col = col, lun = lun)

        CLM.cn_products_compute_product_summary!(prod, bounds, dt)

        # Check mass conservation for 1-yr crop pool
        crop_gain = (prod.dwt_cropprod1_gain_grc[1] +
                     prod.crop_harvest_to_cropprod1_grc[1]) * dt
        crop_loss = prod.cropprod1_loss_grc[1] * dt
        crop_change = prod.cropprod1_grc[1] - crop1_before
        @test crop_change ≈ crop_gain - crop_loss

        # Check mass conservation for 10-yr pool
        p10_gain = (prod.dwt_prod10_gain_grc[1] +
                    prod.gru_prod10_gain_grc[1] +
                    prod.hrv_deadstem_to_prod10_grc[1]) * dt
        p10_loss = prod.prod10_loss_grc[1] * dt
        p10_change = prod.prod10_grc[1] - p10_before
        @test p10_change ≈ p10_gain - p10_loss

        # Check mass conservation for 100-yr pool
        p100_gain = (prod.dwt_prod100_gain_grc[1] +
                     prod.gru_prod100_gain_grc[1] +
                     prod.hrv_deadstem_to_prod100_grc[1]) * dt
        p100_loss = prod.prod100_loss_grc[1] * dt
        p100_change = prod.prod100_grc[1] - p100_before
        @test p100_change ≈ p100_gain - p100_loss
    end

end  # @testset "CN Products Module"
