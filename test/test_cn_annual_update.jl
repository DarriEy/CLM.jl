@testset "CN Annual Update" begin

    # Ensure varpar is initialized (needed by column_init!)
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    # =====================================================================
    # Helper: set up a small domain with nc columns and np patches
    # =====================================================================
    function make_annual_update_data(;
            nc::Int = 2,
            np::Int = 4,
            dt::Float64 = 1800.0,
            days_per_year::Float64 = 365.0,
            # per-column: initial annsum_counter values
            annsum_counter::Vector{Float64} = fill(0.0, nc),
            # per-column: is_fates flags
            is_fates::Vector{Bool} = fill(false, nc),
            # per-patch: which column each patch belongs to
            patch_column::Vector{Int} = Int[],
            # per-patch: tempsum_potential_gpp values
            tempsum_potential_gpp::Vector{Float64} = fill(0.0, np),
            tempmax_retransn::Vector{Float64} = fill(0.0, np),
            tempavg_t2m::Vector{Float64} = fill(0.0, np),
            tempsum_npp::Vector{Float64} = fill(0.0, np),
            tempsum_litfall::Vector{Float64} = fill(0.0, np))

        # Column data
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.is_fates[1:nc] .= is_fates

        # Patch data
        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        if !isempty(patch_column)
            patch.column[1:np] .= patch_column
        end
        # Set active and wtcol for p2c averaging
        patch.active[1:np] .= true
        # Equal weight within each column for simplicity
        patch.wtcol[1:np] .= 0.0
        for p in 1:np
            c = patch.column[p]
            if 1 <= c <= nc
                patch.wtcol[p] = 1.0 / count(==(c), patch_column)
            end
        end

        # Set patchi / patchf for p2c_1d_filter!
        for c in 1:nc
            plist = findall(==(c), patch_column)
            if !isempty(plist)
                col.patchi[c] = first(plist)
                col.patchf[c] = last(plist)
            else
                col.patchi[c] = 1
                col.patchf[c] = 0  # empty range
            end
        end

        # CN veg state
        cnveg_state = CLM.CNVegStateData()
        CLM.cnveg_state_init!(cnveg_state, np, nc)
        cnveg_state.annsum_counter_col[1:nc]          .= annsum_counter
        cnveg_state.tempsum_potential_gpp_patch[1:np]  .= tempsum_potential_gpp
        cnveg_state.annsum_potential_gpp_patch[1:np]   .= 0.0
        cnveg_state.tempmax_retransn_patch[1:np]       .= tempmax_retransn
        cnveg_state.annmax_retransn_patch[1:np]        .= 0.0
        cnveg_state.tempavg_t2m_patch[1:np]            .= tempavg_t2m
        cnveg_state.annavg_t2m_patch[1:np]             .= 0.0
        cnveg_state.annavg_t2m_col[1:nc]               .= 0.0

        # CN veg carbon flux
        cnveg_cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cnveg_cf, np, nc, 1)
        cnveg_cf.tempsum_npp_patch[1:np]     .= tempsum_npp
        cnveg_cf.annsum_npp_patch[1:np]      .= 0.0
        cnveg_cf.tempsum_litfall_patch[1:np]  .= tempsum_litfall
        cnveg_cf.annsum_litfall_patch[1:np]   .= 0.0
        cnveg_cf.annsum_npp_col[1:nc]         .= 0.0

        # Masks: all columns and patches active
        mask_c = trues(nc)
        mask_p = trues(np)

        bounds_c = 1:nc
        bounds_p = 1:np

        return col, patch, cnveg_state, cnveg_cf, mask_c, mask_p, bounds_c, bounds_p
    end

    # =====================================================================
    # Test 1: counter increments but no end-of-year trigger
    # =====================================================================
    @testset "Counter increments, no end-of-year" begin
        dt = 1800.0
        days_per_year = 365.0
        nc = 1; np = 1

        col, patch, vs, cf, mask_c, mask_p, bc, bp = make_annual_update_data(
            nc=nc, np=np,
            annsum_counter = [0.0],
            patch_column = [1],
            tempsum_potential_gpp = [100.0],
            tempsum_npp = [50.0])

        CLM.cn_annual_update!(mask_c, mask_p, bc, bp, col, patch, vs, cf;
                              dt=dt, days_per_year=days_per_year)

        # Counter should have been incremented by dt
        @test vs.annsum_counter_col[1] ≈ dt

        # Annual accumulators should NOT have been updated (still zero from init)
        @test vs.annsum_potential_gpp_patch[1] ≈ 0.0
        @test cf.annsum_npp_patch[1] ≈ 0.0

        # Temporary accumulators should be unchanged
        @test vs.tempsum_potential_gpp_patch[1] ≈ 100.0
        @test cf.tempsum_npp_patch[1] ≈ 50.0
    end

    # =====================================================================
    # Test 2: end-of-year triggers annual update
    # =====================================================================
    @testset "End-of-year triggers annual update" begin
        dt = 1800.0
        days_per_year = 365.0
        secspyear = days_per_year * CLM.SECSPDAY
        nc = 1; np = 2

        # Set counter just below threshold so adding dt pushes it over
        initial_counter = secspyear - dt + 1.0

        col, patch, vs, cf, mask_c, mask_p, bc, bp = make_annual_update_data(
            nc=nc, np=np,
            annsum_counter = [initial_counter],
            patch_column = [1, 1],
            tempsum_potential_gpp = [100.0, 200.0],
            tempmax_retransn = [5.0, 10.0],
            tempavg_t2m = [280.0, 290.0],
            tempsum_npp = [50.0, 75.0],
            tempsum_litfall = [20.0, 30.0])

        CLM.cn_annual_update!(mask_c, mask_p, bc, bp, col, patch, vs, cf;
                              dt=dt, days_per_year=days_per_year)

        # Counter should have been reset to zero
        @test vs.annsum_counter_col[1] ≈ 0.0

        # Annual accumulators should be set from temporary values
        @test vs.annsum_potential_gpp_patch[1] ≈ 100.0
        @test vs.annsum_potential_gpp_patch[2] ≈ 200.0

        @test vs.annmax_retransn_patch[1] ≈ 5.0
        @test vs.annmax_retransn_patch[2] ≈ 10.0

        @test vs.annavg_t2m_patch[1] ≈ 280.0
        @test vs.annavg_t2m_patch[2] ≈ 290.0

        # NPP and litfall are multiplied by dt when promoted
        @test cf.annsum_npp_patch[1] ≈ 50.0 * dt
        @test cf.annsum_npp_patch[2] ≈ 75.0 * dt

        @test cf.annsum_litfall_patch[1] ≈ 20.0 * dt
        @test cf.annsum_litfall_patch[2] ≈ 30.0 * dt

        # Temporary accumulators should be reset to zero
        @test vs.tempsum_potential_gpp_patch[1] ≈ 0.0
        @test vs.tempsum_potential_gpp_patch[2] ≈ 0.0
        @test vs.tempmax_retransn_patch[1] ≈ 0.0
        @test vs.tempmax_retransn_patch[2] ≈ 0.0
        @test vs.tempavg_t2m_patch[1] ≈ 0.0
        @test vs.tempavg_t2m_patch[2] ≈ 0.0
        @test cf.tempsum_npp_patch[1] ≈ 0.0
        @test cf.tempsum_npp_patch[2] ≈ 0.0
        @test cf.tempsum_litfall_patch[1] ≈ 0.0
        @test cf.tempsum_litfall_patch[2] ≈ 0.0

        # Column-level p2c averaging: equal-weight patches
        # annsum_npp_col should be average of the two patches
        expected_npp_col = 0.5 * (50.0 * dt) + 0.5 * (75.0 * dt)
        @test cf.annsum_npp_col[1] ≈ expected_npp_col

        # annavg_t2m_col should be average of the two patches
        expected_t2m_col = 0.5 * 280.0 + 0.5 * 290.0
        @test vs.annavg_t2m_col[1] ≈ expected_t2m_col
    end

    # =====================================================================
    # Test 3: FATES columns are skipped
    # =====================================================================
    @testset "FATES columns skipped" begin
        dt = 1800.0
        days_per_year = 365.0
        secspyear = days_per_year * CLM.SECSPDAY
        nc = 2; np = 2

        initial_counter = secspyear - dt + 1.0

        col, patch, vs, cf, mask_c, mask_p, bc, bp = make_annual_update_data(
            nc=nc, np=np,
            annsum_counter = [initial_counter, initial_counter],
            is_fates = [true, false],
            patch_column = [1, 2],
            tempsum_potential_gpp = [100.0, 200.0],
            tempsum_npp = [50.0, 75.0],
            tempsum_litfall = [20.0, 30.0],
            tempmax_retransn = [5.0, 10.0],
            tempavg_t2m = [280.0, 290.0])

        CLM.cn_annual_update!(mask_c, mask_p, bc, bp, col, patch, vs, cf;
                              dt=dt, days_per_year=days_per_year)

        # Column 1 (FATES): counter should NOT have been incremented
        @test vs.annsum_counter_col[1] ≈ initial_counter

        # Patch 1 (FATES column): should NOT have been updated
        @test vs.annsum_potential_gpp_patch[1] ≈ 0.0
        @test vs.tempsum_potential_gpp_patch[1] ≈ 100.0

        # Column 2 (non-FATES): counter should be reset
        @test vs.annsum_counter_col[2] ≈ 0.0

        # Patch 2 (non-FATES column): should have been updated
        @test vs.annsum_potential_gpp_patch[2] ≈ 200.0
        @test vs.tempsum_potential_gpp_patch[2] ≈ 0.0
        @test cf.annsum_npp_patch[2] ≈ 75.0 * dt
        @test cf.tempsum_npp_patch[2] ≈ 0.0
    end

    # =====================================================================
    # Test 4: masked-out columns/patches are skipped
    # =====================================================================
    @testset "Masked-out columns/patches skipped" begin
        dt = 1800.0
        days_per_year = 365.0
        secspyear = days_per_year * CLM.SECSPDAY
        nc = 2; np = 2

        initial_counter = secspyear - dt + 1.0

        col, patch, vs, cf, mask_c, mask_p, bc, bp = make_annual_update_data(
            nc=nc, np=np,
            annsum_counter = [initial_counter, initial_counter],
            patch_column = [1, 2],
            tempsum_potential_gpp = [100.0, 200.0],
            tempsum_npp = [50.0, 75.0],
            tempsum_litfall = [20.0, 30.0],
            tempmax_retransn = [5.0, 10.0],
            tempavg_t2m = [280.0, 290.0])

        # Mask out column 1 and patch 1
        mask_c[1] = false
        mask_p[1] = false

        CLM.cn_annual_update!(mask_c, mask_p, bc, bp, col, patch, vs, cf;
                              dt=dt, days_per_year=days_per_year)

        # Column 1: counter unchanged (masked out)
        @test vs.annsum_counter_col[1] ≈ initial_counter

        # Patch 1: not updated (masked out)
        @test vs.annsum_potential_gpp_patch[1] ≈ 0.0
        @test vs.tempsum_potential_gpp_patch[1] ≈ 100.0

        # Column 2: counter reset (end of year)
        @test vs.annsum_counter_col[2] ≈ 0.0

        # Patch 2: updated
        @test vs.annsum_potential_gpp_patch[2] ≈ 200.0
        @test vs.tempsum_potential_gpp_patch[2] ≈ 0.0
    end

    # =====================================================================
    # Test 5: mixed columns - one at end-of-year, one not
    # =====================================================================
    @testset "Mixed columns - partial end-of-year" begin
        dt = 1800.0
        days_per_year = 365.0
        secspyear = days_per_year * CLM.SECSPDAY
        nc = 2; np = 2

        # Column 1: about to reach end-of-year
        # Column 2: far from end-of-year
        col1_counter = secspyear - dt + 1.0
        col2_counter = 1000.0

        col, patch, vs, cf, mask_c, mask_p, bc, bp = make_annual_update_data(
            nc=nc, np=np,
            annsum_counter = [col1_counter, col2_counter],
            patch_column = [1, 2],
            tempsum_potential_gpp = [100.0, 200.0],
            tempsum_npp = [50.0, 75.0],
            tempsum_litfall = [20.0, 30.0],
            tempmax_retransn = [5.0, 10.0],
            tempavg_t2m = [280.0, 290.0])

        CLM.cn_annual_update!(mask_c, mask_p, bc, bp, col, patch, vs, cf;
                              dt=dt, days_per_year=days_per_year)

        # Column 1: end of year reached, counter reset
        @test vs.annsum_counter_col[1] ≈ 0.0
        # Patch 1: annual values updated
        @test vs.annsum_potential_gpp_patch[1] ≈ 100.0
        @test vs.tempsum_potential_gpp_patch[1] ≈ 0.0
        @test cf.annsum_npp_patch[1] ≈ 50.0 * dt

        # Column 2: not end of year, counter incremented
        @test vs.annsum_counter_col[2] ≈ col2_counter + dt
        # Patch 2: annual values NOT updated
        @test vs.annsum_potential_gpp_patch[2] ≈ 0.0
        @test vs.tempsum_potential_gpp_patch[2] ≈ 200.0
        @test cf.annsum_npp_patch[2] ≈ 0.0
    end

    # =====================================================================
    # Test 6: no BGC veg patches (empty patch loop)
    # =====================================================================
    @testset "No BGC veg patches" begin
        dt = 1800.0
        days_per_year = 365.0
        secspyear = days_per_year * CLM.SECSPDAY
        nc = 1; np = 1

        initial_counter = secspyear - dt + 1.0

        col, patch, vs, cf, mask_c, mask_p, bc, bp = make_annual_update_data(
            nc=nc, np=np,
            annsum_counter = [initial_counter],
            patch_column = [1],
            tempsum_potential_gpp = [100.0],
            tempsum_npp = [50.0],
            tempsum_litfall = [20.0])

        # Mask out all patches
        mask_p[1] = false

        CLM.cn_annual_update!(mask_c, mask_p, bc, bp, col, patch, vs, cf;
                              dt=dt, days_per_year=days_per_year)

        # Column counter should still be reset (column loop is independent)
        @test vs.annsum_counter_col[1] ≈ 0.0

        # But patch-level accumulators should not have changed
        @test vs.annsum_potential_gpp_patch[1] ≈ 0.0
        @test vs.tempsum_potential_gpp_patch[1] ≈ 100.0
        @test cf.annsum_npp_patch[1] ≈ 0.0
        @test cf.tempsum_npp_patch[1] ≈ 50.0
    end

    # =====================================================================
    # Test 7: multiple patches per column with different weights
    # =====================================================================
    @testset "p2c averaging with unequal patch weights" begin
        dt = 1800.0
        days_per_year = 365.0
        secspyear = days_per_year * CLM.SECSPDAY
        nc = 1; np = 3

        initial_counter = secspyear - dt + 1.0

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.is_fates[1] = false
        col.patchi[1] = 1
        col.patchf[1] = 3

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.column[1:3] .= 1
        patch.active[1:3] .= true
        # Unequal weights: 0.5, 0.3, 0.2
        patch.wtcol[1] = 0.5
        patch.wtcol[2] = 0.3
        patch.wtcol[3] = 0.2

        vs = CLM.CNVegStateData()
        CLM.cnveg_state_init!(vs, np, nc)
        vs.annsum_counter_col[1] = initial_counter
        vs.tempsum_potential_gpp_patch[1:3] .= [100.0, 200.0, 300.0]
        vs.annsum_potential_gpp_patch[1:3]  .= 0.0
        vs.tempmax_retransn_patch[1:3]      .= [1.0, 2.0, 3.0]
        vs.annmax_retransn_patch[1:3]       .= 0.0
        vs.tempavg_t2m_patch[1:3]           .= [270.0, 280.0, 290.0]
        vs.annavg_t2m_patch[1:3]            .= 0.0
        vs.annavg_t2m_col[1]                = 0.0

        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, nc, 1)
        cf.tempsum_npp_patch[1:3]     .= [10.0, 20.0, 30.0]
        cf.annsum_npp_patch[1:3]      .= 0.0
        cf.tempsum_litfall_patch[1:3]  .= [5.0, 10.0, 15.0]
        cf.annsum_litfall_patch[1:3]   .= 0.0
        cf.annsum_npp_col[1]           = 0.0

        mask_c = trues(nc)
        mask_p = trues(np)

        CLM.cn_annual_update!(mask_c, mask_p, 1:nc, 1:np, col, patch, vs, cf;
                              dt=dt, days_per_year=days_per_year)

        # Check patch-level updates
        @test vs.annsum_potential_gpp_patch[1] ≈ 100.0
        @test vs.annsum_potential_gpp_patch[2] ≈ 200.0
        @test vs.annsum_potential_gpp_patch[3] ≈ 300.0

        @test cf.annsum_npp_patch[1] ≈ 10.0 * dt
        @test cf.annsum_npp_patch[2] ≈ 20.0 * dt
        @test cf.annsum_npp_patch[3] ≈ 30.0 * dt

        # Column-level p2c: weighted average
        expected_npp_col = 0.5 * (10.0 * dt) + 0.3 * (20.0 * dt) + 0.2 * (30.0 * dt)
        @test cf.annsum_npp_col[1] ≈ expected_npp_col

        expected_t2m_col = 0.5 * 270.0 + 0.3 * 280.0 + 0.2 * 290.0
        @test vs.annavg_t2m_col[1] ≈ expected_t2m_col
    end

    # =====================================================================
    # Test 8: exact boundary - counter exactly equals secspyear
    # =====================================================================
    @testset "Counter exactly at secspyear boundary" begin
        dt = 1800.0
        days_per_year = 365.0
        secspyear = days_per_year * CLM.SECSPDAY
        nc = 1; np = 1

        # Counter is exactly secspyear - dt, so after adding dt it equals secspyear
        initial_counter = secspyear - dt

        col, patch, vs, cf, mask_c, mask_p, bc, bp = make_annual_update_data(
            nc=nc, np=np,
            annsum_counter = [initial_counter],
            patch_column = [1],
            tempsum_potential_gpp = [42.0],
            tempsum_npp = [21.0],
            tempsum_litfall = [7.0],
            tempmax_retransn = [3.0],
            tempavg_t2m = [285.0])

        CLM.cn_annual_update!(mask_c, mask_p, bc, bp, col, patch, vs, cf;
                              dt=dt, days_per_year=days_per_year)

        # Should trigger end-of-year (>= comparison)
        @test vs.annsum_counter_col[1] ≈ 0.0
        @test vs.annsum_potential_gpp_patch[1] ≈ 42.0
        @test vs.annmax_retransn_patch[1] ≈ 3.0
        @test vs.annavg_t2m_patch[1] ≈ 285.0
        @test cf.annsum_npp_patch[1] ≈ 21.0 * dt
        @test cf.annsum_litfall_patch[1] ≈ 7.0 * dt
    end

    # =====================================================================
    # Test 9: dt=0 produces no change
    # =====================================================================
    @testset "Zero dt produces no change" begin
        days_per_year = 365.0
        nc = 1; np = 1

        col, patch, vs, cf, mask_c, mask_p, bc, bp = make_annual_update_data(
            nc=nc, np=np,
            annsum_counter = [1000.0],
            patch_column = [1],
            tempsum_potential_gpp = [100.0],
            tempsum_npp = [50.0],
            tempsum_litfall = [20.0])

        CLM.cn_annual_update!(mask_c, mask_p, bc, bp, col, patch, vs, cf;
                              dt=0.0, days_per_year=days_per_year)

        # Counter should not have changed
        @test vs.annsum_counter_col[1] ≈ 1000.0
        # No end-of-year: patch accumulators unchanged
        @test vs.tempsum_potential_gpp_patch[1] ≈ 100.0
        @test vs.annsum_potential_gpp_patch[1] ≈ 0.0
    end

end
