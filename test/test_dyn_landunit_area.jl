@testset "dyn_landunit_area" begin

    # Needed so column_init!/patch_init! can size snow-layer arrays.
    CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)
    CLM.varcon_init!()

    tol = 1.0e-14

    # ----------------------------------------------------------------------
    # update_landunit_weights_one_gcell!
    # ----------------------------------------------------------------------
    @testset "update_landunit_weights_one_gcell!" begin
        # Helper: fresh zero weight vector of length MAX_LUNIT.
        zerow() = zeros(Float64, CLM.MAX_LUNIT)

        # --- Already sums to 1: untouched ----------------------------------
        w = zerow()
        w[CLM.ISTSOIL] = 0.6
        w[CLM.ISTCROP] = 0.4
        before = copy(w)
        CLM.update_landunit_weights_one_gcell!(w)
        @test w == before
        @test abs(sum(w) - 1.0) <= tol

        # --- Sums to < 1: ISTSOIL absorbs the deficit ----------------------
        w = zerow()
        w[CLM.ISTSOIL] = 0.2
        w[CLM.ISTCROP] = 0.3   # sum = 0.5
        CLM.update_landunit_weights_one_gcell!(w)
        @test abs(sum(w) - 1.0) <= tol
        @test w[CLM.ISTSOIL] ≈ 0.7        # 0.2 + (1 - 0.5)
        @test w[CLM.ISTCROP] ≈ 0.3        # unchanged

        # --- Sums to < 1 with zero soil: soil grows from 0 -----------------
        w = zerow()
        w[CLM.ISTDLAK] = 0.3
        CLM.update_landunit_weights_one_gcell!(w)
        @test abs(sum(w) - 1.0) <= tol
        @test w[CLM.ISTSOIL] ≈ 0.7
        @test w[CLM.ISTDLAK] ≈ 0.3

        # --- Sums to > 1: soil decreased FIRST (priority head) -------------
        w = zerow()
        w[CLM.ISTSOIL] = 0.5
        w[CLM.ISTCROP] = 0.4
        w[CLM.ISTDLAK] = 0.3   # sum = 1.2, excess = 0.2
        CLM.update_landunit_weights_one_gcell!(w)
        @test abs(sum(w) - 1.0) <= tol
        @test w[CLM.ISTSOIL] ≈ 0.3        # 0.5 - 0.2 excess
        @test w[CLM.ISTCROP] ≈ 0.4        # untouched (excess absorbed by soil)
        @test w[CLM.ISTDLAK] ≈ 0.3        # untouched

        # --- Overflow: first priority landunit clamps at 0, spills to next -
        # soil=0.1, crop=0.4, lake=1.0 -> sum=1.5, excess=0.5.
        # decrease soil first: 0.1-0.5 -> clamp 0 (removed 0.1), excess now 0.4
        # decrease crop next:  0.4-0.4 -> 0.0   (removed 0.4), excess now 0.0
        w = zerow()
        w[CLM.ISTSOIL] = 0.1
        w[CLM.ISTCROP] = 0.4
        w[CLM.ISTDLAK] = 1.0
        CLM.update_landunit_weights_one_gcell!(w)
        @test abs(sum(w) - 1.0) <= tol
        @test isapprox(w[CLM.ISTSOIL], 0.0; atol=1e-14)
        @test isapprox(w[CLM.ISTCROP], 0.0; atol=1e-14)
        @test w[CLM.ISTDLAK] ≈ 1.0        # never decreased; lake last in order

        # --- istice (type 4, ISTICE) is excluded from decrease_order -------
        # It is never reduced; soil/crop/etc absorb the excess instead.
        w = zerow()
        w[CLM.ISTSOIL] = 0.5
        w[CLM.ISTICE]  = 0.8   # sum = 1.3, excess = 0.3
        CLM.update_landunit_weights_one_gcell!(w)
        @test abs(sum(w) - 1.0) <= tol
        @test w[CLM.ISTICE] ≈ 0.8         # istice preserved exactly
        @test w[CLM.ISTSOIL] ≈ 0.2        # soil absorbed the 0.3 excess

        # --- Decrease order: urban md before hd before tbd, then wet, lake -
        # Put weight only on a single non-soil landunit plus enough to overflow,
        # confirming order by checking which gets decreased.
        w = zerow()
        w[CLM.ISTURB_MD]  = 0.6
        w[CLM.ISTURB_HD]  = 0.6   # sum = 1.2, excess = 0.2; soil=0 first (clamp 0)
        CLM.update_landunit_weights_one_gcell!(w)
        @test abs(sum(w) - 1.0) <= tol
        # soil clamps 0 (no change to excess since it was already 0),
        # crop clamps 0, then ISTURB_MD decreased: 0.6 - 0.2 = 0.4
        @test w[CLM.ISTURB_MD] ≈ 0.4
        @test w[CLM.ISTURB_HD] ≈ 0.6      # untouched (later in order)

        # --- Wrong length asserts -----------------------------------------
        @test_throws AssertionError CLM.update_landunit_weights_one_gcell!(
            zeros(Float64, CLM.MAX_LUNIT - 1))
    end

    # ----------------------------------------------------------------------
    # update_landunit_weights! over a small synthetic gridcell
    # ----------------------------------------------------------------------
    @testset "update_landunit_weights!" begin
        # One gridcell, 3 landunits: soil, crop, lake. Weights don't sum to 1.
        ng, nl = 1, 3
        grc = CLM.GridcellData()
        lun = CLM.LandunitData()
        CLM.gridcell_init!(grc, ng)
        CLM.landunit_init!(lun, nl)

        # landunit 1 = soil, 2 = crop, 3 = lake on gridcell 1
        lun.gridcell .= 1
        lun.itype[1] = CLM.ISTSOIL
        lun.itype[2] = CLM.ISTCROP
        lun.itype[3] = CLM.ISTDLAK
        grc.landunit_indices[CLM.ISTSOIL, 1] = 1
        grc.landunit_indices[CLM.ISTCROP, 1] = 2
        grc.landunit_indices[CLM.ISTDLAK, 1] = 3

        # Deficit case: 0.3 + 0.2 + 0.1 = 0.6 -> soil grows to absorb 0.4.
        lun.wtgcell[1] = 0.3
        lun.wtgcell[2] = 0.2
        lun.wtgcell[3] = 0.1

        bounds = CLM.BoundsType(begg=1, endg=1, begl=1, endl=3,
                                level=CLM.BOUNDS_LEVEL_PROC)
        CLM.update_landunit_weights!(bounds, grc, lun)

        @test abs(sum(lun.wtgcell) - 1.0) <= tol
        @test lun.wtgcell[1] ≈ 0.7   # soil absorbed the residual
        @test lun.wtgcell[2] ≈ 0.2
        @test lun.wtgcell[3] ≈ 0.1

        # Excess case: 0.5 + 0.4 + 0.3 = 1.2 -> soil decreased first.
        lun.wtgcell[1] = 0.5
        lun.wtgcell[2] = 0.4
        lun.wtgcell[3] = 0.3
        CLM.update_landunit_weights!(bounds, grc, lun)
        @test abs(sum(lun.wtgcell) - 1.0) <= tol
        @test lun.wtgcell[1] ≈ 0.3
        @test lun.wtgcell[2] ≈ 0.4
        @test lun.wtgcell[3] ≈ 0.3
    end

    # ----------------------------------------------------------------------
    # PriorWeights snapshot
    # ----------------------------------------------------------------------
    @testset "PriorWeights" begin
        # 4 patches, 3 columns.
        np, nc = 4, 3
        bounds = CLM.BoundsType(begp=1, endp=np, begc=1, endc=nc,
                                level=CLM.BOUNDS_LEVEL_PROC)

        prior = CLM.PriorWeights(bounds)
        @test length(prior.pwtgcell) == np
        @test length(prior.cactive) == nc

        pch = CLM.PatchData()
        col = CLM.ColumnData()
        CLM.patch_init!(pch, np)
        CLM.column_init!(col, nc)

        pch.wtgcell .= [0.1, 0.2, 0.3, 0.4]
        col.active  .= [true, false, true]

        CLM.set_prior_weights!(prior, bounds, pch, col)
        @test prior.pwtgcell == [0.1, 0.2, 0.3, 0.4]
        @test collect(prior.cactive) == [true, false, true]

        # Mutate the originals; the snapshot must NOT change.
        pch.wtgcell .= [0.9, 0.9, 0.9, 0.9]
        col.active  .= [false, true, false]
        @test prior.pwtgcell == [0.1, 0.2, 0.3, 0.4]
        @test collect(prior.cactive) == [true, false, true]

        # Constructor rejects non-PROC bounds.
        @test_throws AssertionError CLM.PriorWeights(
            CLM.BoundsType(begp=1, endp=np, begc=1, endc=nc,
                           level=CLM.BOUNDS_LEVEL_CLUMP))
    end
end
