@testset "dyn_patch_state_updater" begin

    tol = 1.0e-12

    # ----------------------------------------------------------------------
    # Synthetic subgrid: 4 patches, each on its own column, all in one gcell.
    #
    #   p=1  shrinking   : 0.40 -> 0.20
    #   p=2  growing     : 0.10 -> 0.30
    #   p=3  initiating  : 0.00 -> 0.25  (grows from zero area)
    #   p=4  unchanged   : 0.25 -> 0.25
    #
    # Column weights mirror the patch weights (one patch per column), and we
    # use a single gridcell so "mass per unit area gridcell" sums directly.
    # ----------------------------------------------------------------------
    npatch = 4
    ncol   = 4

    bounds = CLM.BoundsType(
        begg = 1, endg = 1,
        begl = 1, endl = 1,
        begc = 1, endc = ncol,
        begp = 1, endp = npatch,
        level = CLM.BOUNDS_LEVEL_PROC,
    )

    # Old / new patch weights on the gridcell.
    pwt_old = [0.40, 0.10, 0.00, 0.25]
    pwt_new = [0.20, 0.30, 0.25, 0.25]

    # Each patch sits on its own column; column gridcell weight == patch weight.
    pch = CLM.PatchData{Float64}()
    pch.column   = collect(1:npatch)
    pch.itype    = [7, 13, 1, 0]    # arbitrary PFT types (0-based; bare ground = 0)
    pch.wtgcell  = copy(pwt_old)

    col = CLM.ColumnData{Float64}()
    col.wtgcell  = copy(pwt_old)

    # Filter that includes all patches (incl. ones that just became inactive).
    filterp = collect(1:npatch)

    # ----------------------------------------------------------------------
    # set_old_weights! / set_new_weights!
    # ----------------------------------------------------------------------
    updater = CLM.PatchStateUpdater(bounds)
    @test length(updater.pwtgcell_old) == npatch
    @test length(updater.cwtgcell_old) == ncol

    CLM.set_old_weights!(updater, bounds, pch, col)
    @test updater.pwtgcell_old == pwt_old
    @test updater.cwtgcell_old == pwt_old

    # Apply the new weights to patch & column, then snapshot.
    pch.wtgcell .= pwt_new
    CLM.set_new_weights!(updater, bounds, pch)
    @test updater.pwtgcell_new == pwt_new
    @test updater.dwt ≈ (pwt_new .- pwt_old)

    # growing fractions: valid only for growing patches (p=2, p=3)
    @test updater.growing_old_fraction[2] ≈ pwt_old[2] / pwt_new[2]
    @test updater.growing_new_fraction[2] ≈ (pwt_new[2] - pwt_old[2]) / pwt_new[2]
    @test updater.growing_old_fraction[3] ≈ 0.0 / pwt_new[3]   # initiating -> 0
    @test updater.growing_new_fraction[3] ≈ 1.0                # all new area
    # safe defaults for non-growing patches
    @test updater.growing_old_fraction[1] == 1.0
    @test updater.growing_new_fraction[1] == 0.0
    @test updater.growing_old_fraction[4] == 1.0
    @test updater.growing_new_fraction[4] == 0.0

    # ----------------------------------------------------------------------
    # Query helpers
    # ----------------------------------------------------------------------
    @testset "query helpers" begin
        owz  = CLM.old_weight_was_zero(updater, bounds)
        grew = CLM.patch_grew(updater, bounds)
        init = CLM.patch_initiating(updater, bounds)

        @test owz  == BitVector([false, false, true,  false])
        @test grew == BitVector([false, true,  true,  false])
        @test init == BitVector([false, false, true,  false])
    end

    # ----------------------------------------------------------------------
    # update_patch_state! — mass conservation (no seed).
    #
    # Total state*area before == total after + (-flux_out) accumulated for
    # shrinking patches.  We track flux as mass-per-unit-area-gridcell so the
    # accumulated flux is directly comparable to var*pwt.
    # ----------------------------------------------------------------------
    @testset "update_patch_state! conservation (no seed)" begin
        var0 = [10.0, 20.0, 0.0, 5.0]   # per-area state; p=3 starts empty (new area)
        var  = copy(var0)

        flux_grc = zeros(Float64, npatch)
        flux_col = zeros(Float64, npatch)

        CLM.update_patch_state!(updater, bounds, pch, filterp, var;
                                flux_out_col_area = flux_col,
                                flux_out_grc_area = flux_grc)

        # --- per-patch behavior ---
        # p=1 shrinking: per-area state unchanged, flux out = var*dwt (negative)
        @test var[1] ≈ var0[1]
        @test flux_grc[1] ≈ var0[1] * (pwt_new[1] - pwt_old[1])     # negative
        @test flux_grc[1] < 0.0
        # flux_out_col_area uses OLD column weight
        @test flux_col[1] ≈ var0[1] * ((pwt_new[1] - pwt_old[1]) / pwt_old[1])

        # p=2 growing (no seed): per-area diluted by growing_old_fraction
        @test var[2] ≈ var0[2] * (pwt_old[2] / pwt_new[2])
        @test flux_grc[2] == 0.0       # growing patches accumulate no flux out

        # p=3 initiating, no seed: stays at 0
        @test var[3] == 0.0
        @test flux_grc[3] == 0.0

        # p=4 unchanged: nothing happens
        @test var[4] ≈ var0[4]
        @test flux_grc[4] == 0.0

        # --- mass conservation (gridcell-area basis) ---
        # mass_before = sum var0*pwt_old ; mass_after = sum var*pwt_new
        # The flux_out_grc_area accounts EXACTLY for the difference (no seed).
        mass_before = sum(var0 .* pwt_old)
        mass_after  = sum(var  .* pwt_new)
        total_flux  = sum(flux_grc)     # negative: mass that left shrinking patches
        @test isapprox(mass_after - total_flux, mass_before; atol = tol)
    end

    # ----------------------------------------------------------------------
    # update_patch_state! — with seed + seed_addition.
    # Growing/initiating patches receive seed over their new-area fraction;
    # seed_addition accumulates seed*dwt (gridcell-area basis).
    # ----------------------------------------------------------------------
    @testset "update_patch_state! with seed" begin
        # Reset state on the updater-consistent patches.
        var0 = [10.0, 20.0, 0.0, 5.0]
        var  = copy(var0)

        seed     = [0.0, 100.0, 50.0, 0.0]      # seed amount (per new area)
        seed_add = zeros(Float64, npatch)
        flux_grc = zeros(Float64, npatch)

        CLM.update_patch_state!(updater, bounds, pch, filterp, var;
                                flux_out_grc_area = flux_grc,
                                seed = seed, seed_addition = seed_add)

        # p=2 growing with seed
        @test var[2] ≈ var0[2] * (pwt_old[2] / pwt_new[2]) +
                       seed[2] * ((pwt_new[2] - pwt_old[2]) / pwt_new[2])
        @test seed_add[2] ≈ seed[2] * (pwt_new[2] - pwt_old[2])

        # p=3 initiating: var becomes seed*growing_new_fraction (=seed*1.0)
        @test var[3] ≈ seed[3] * 1.0
        @test seed_add[3] ≈ seed[3] * (pwt_new[3] - pwt_old[3])

        # Conservation WITH seeding: mass_after = mass_before - total_flux + seed_added
        mass_before = sum(var0 .* pwt_old)
        mass_after  = sum(var  .* pwt_new)
        total_flux  = sum(flux_grc)
        seed_total  = sum(seed_add)
        @test isapprox(mass_after - total_flux, mass_before + seed_total; atol = tol)
    end

    # seed_addition without seed must error
    @testset "seed_addition requires seed" begin
        var = [1.0, 1.0, 1.0, 1.0]
        sadd = zeros(Float64, npatch)
        @test_throws ErrorException CLM.update_patch_state!(
            updater, bounds, pch, filterp, var; seed_addition = sadd)
    end

    # ----------------------------------------------------------------------
    # update_patch_state_partition_flux_by_type!
    # Partition the (negative) flux out of shrinking patches by PFT type.
    # ----------------------------------------------------------------------
    @testset "update_patch_state_partition_flux_by_type!" begin
        # flux1 fraction indexed by PFT type (Fortran 0:mxpft -> length MXPFT+1)
        frac = fill(0.25, CLM.MXPFT + 1)
        frac[pch.itype[1] + 1] = 0.75    # p=1 (the only shrinking patch) type -> 0.75

        var0 = [10.0, 20.0, 0.0, 5.0]
        var  = copy(var0)
        flux1 = zeros(Float64, npatch)
        flux2 = zeros(Float64, npatch)

        CLM.update_patch_state_partition_flux_by_type!(
            updater, bounds, pch, filterp, frac, var, flux1, flux2)

        # Total flux out of p=1 (gridcell-area basis):
        total_p1 = var0[1] * (pwt_new[1] - pwt_old[1])   # negative
        @test flux1[1] ≈ total_p1 * 0.75
        @test flux2[1] ≈ total_p1 * (1.0 - 0.75)
        # flux1 + flux2 reconstructs the full flux out for every patch (only
        # shrinking patches have non-zero flux out).
        for p in 1:npatch
            expected = (pwt_new[p] < pwt_old[p]) ? var0[p] * (pwt_new[p] - pwt_old[p]) : 0.0
            @test isapprox(flux1[p] + flux2[p], expected; atol = tol)
        end
        # Non-shrinking patches contribute nothing.
        @test flux1[2] == 0.0 && flux2[2] == 0.0
        @test flux1[3] == 0.0 && flux2[3] == 0.0
        @test flux1[4] == 0.0 && flux2[4] == 0.0

        # Wrong-length fraction array must error.
        @test_throws AssertionError CLM.update_patch_state_partition_flux_by_type!(
            updater, bounds, pch, filterp, fill(0.5, 3), copy(var0), flux1, flux2)
    end

end
