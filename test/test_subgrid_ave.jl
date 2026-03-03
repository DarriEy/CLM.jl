@testset "SubgridAve" begin

    # =====================================================================
    # Helper: build a minimal subgrid hierarchy for testing
    # =====================================================================
    # Layout:
    #   1 gridcell  →  2 landunits  →  3 columns  →  4 patches
    #
    #   gridcell 1
    #   ├── landunit 1 (soil, wtgcell=0.6)
    #   │   ├── column 1 (wtlunit=0.4, wtgcell=0.24)
    #   │   │   └── patch 1 (wtcol=1.0, wtlunit=0.4, wtgcell=0.24)
    #   │   └── column 2 (wtlunit=0.6, wtgcell=0.36)
    #   │       └── patch 2 (wtcol=1.0, wtlunit=0.6, wtgcell=0.36)
    #   └── landunit 2 (crop, wtgcell=0.4)
    #       └── column 3 (wtlunit=1.0, wtgcell=0.4)
    #           ├── patch 3 (wtcol=0.7, wtlunit=0.7, wtgcell=0.28)
    #           └── patch 4 (wtcol=0.3, wtlunit=0.3, wtgcell=0.12)
    function make_test_hierarchy()
        ngrc = 1; nlun = 2; ncol = 3; npch = 4

        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ngrc)
        grc.active[1] = true

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nlun)
        lun.itype[1] = CLM.ISTSOIL
        lun.itype[2] = CLM.ISTCROP
        lun.urbpoi[1] = false
        lun.urbpoi[2] = false
        lun.active[1] = true
        lun.active[2] = true
        lun.gridcell[1] = 1
        lun.gridcell[2] = 1
        lun.wtgcell[1] = 0.6
        lun.wtgcell[2] = 0.4
        lun.coli[1] = 1
        lun.colf[1] = 2
        lun.ncolumns[1] = 2
        lun.coli[2] = 3
        lun.colf[2] = 3
        lun.ncolumns[2] = 1
        lun.patchi[1] = 1
        lun.patchf[1] = 2
        lun.npatches[1] = 2
        lun.patchi[2] = 3
        lun.patchf[2] = 4
        lun.npatches[2] = 2

        col = CLM.ColumnData()
        CLM.column_init!(col, ncol)
        col.landunit[1] = 1; col.landunit[2] = 1; col.landunit[3] = 2
        col.gridcell[1] = 1; col.gridcell[2] = 1; col.gridcell[3] = 1
        col.active[1] = true; col.active[2] = true; col.active[3] = true
        col.itype[1] = 1; col.itype[2] = 1; col.itype[3] = 1
        col.wtlunit[1] = 0.4; col.wtlunit[2] = 0.6; col.wtlunit[3] = 1.0
        col.wtgcell[1] = 0.24; col.wtgcell[2] = 0.36; col.wtgcell[3] = 0.4
        col.patchi[1] = 1; col.patchi[2] = 2; col.patchi[3] = 3
        col.patchf[1] = 1; col.patchf[2] = 2; col.patchf[3] = 4
        col.npatches[1] = 1; col.npatches[2] = 1; col.npatches[3] = 2

        pch = CLM.PatchData()
        CLM.patch_init!(pch, npch)
        pch.column[1] = 1; pch.column[2] = 2; pch.column[3] = 3; pch.column[4] = 3
        pch.landunit[1] = 1; pch.landunit[2] = 1; pch.landunit[3] = 2; pch.landunit[4] = 2
        pch.gridcell[1] = 1; pch.gridcell[2] = 1; pch.gridcell[3] = 1; pch.gridcell[4] = 1
        pch.active[1] = true; pch.active[2] = true; pch.active[3] = true; pch.active[4] = true
        pch.wtcol[1] = 1.0; pch.wtcol[2] = 1.0; pch.wtcol[3] = 0.7; pch.wtcol[4] = 0.3
        pch.wtlunit[1] = 0.4; pch.wtlunit[2] = 0.6; pch.wtlunit[3] = 0.7; pch.wtlunit[4] = 0.3
        pch.wtgcell[1] = 0.24; pch.wtgcell[2] = 0.36; pch.wtgcell[3] = 0.28; pch.wtgcell[4] = 0.12

        bounds = CLM.BoundsType(begg=1, endg=ngrc, begl=1, endl=nlun,
                                begc=1, endc=ncol, begp=1, endp=npch)

        return grc, lun, col, pch, bounds
    end

    # =====================================================================
    # create_scale_l2g_lookup
    # =====================================================================
    @testset "create_scale_l2g_lookup" begin
        # unity: all 1.0
        lu = CLM.create_scale_l2g_lookup("unity")
        @test length(lu) == CLM.MAX_LUNIT
        @test all(x -> x == 1.0, lu)

        # natveg: only soil
        lu = CLM.create_scale_l2g_lookup("natveg")
        @test lu[CLM.ISTSOIL] == 1.0
        @test lu[CLM.ISTCROP] == CLM.SPVAL

        # veg: soil + crop
        lu = CLM.create_scale_l2g_lookup("veg")
        @test lu[CLM.ISTSOIL] == 1.0
        @test lu[CLM.ISTCROP] == 1.0
        @test lu[CLM.ISTICE] == CLM.SPVAL

        # ice: only ice
        lu = CLM.create_scale_l2g_lookup("ice")
        @test lu[CLM.ISTICE] == 1.0
        @test lu[CLM.ISTSOIL] == CLM.SPVAL

        # nonurb: all except urban
        lu = CLM.create_scale_l2g_lookup("nonurb")
        @test lu[CLM.ISTSOIL] == 1.0
        @test lu[CLM.ISTCROP] == 1.0
        for i in CLM.ISTURB_MIN:CLM.ISTURB_MAX
            @test lu[i] == CLM.SPVAL
        end

        # lake: only lake
        lu = CLM.create_scale_l2g_lookup("lake")
        @test lu[CLM.ISTDLAK] == 1.0
        @test lu[CLM.ISTSOIL] == CLM.SPVAL

        # veg_plus_lake: soil + crop + lake
        lu = CLM.create_scale_l2g_lookup("veg_plus_lake")
        @test lu[CLM.ISTSOIL] == 1.0
        @test lu[CLM.ISTCROP] == 1.0
        @test lu[CLM.ISTDLAK] == 1.0
        @test lu[CLM.ISTICE] == CLM.SPVAL

        # unsupported → error
        @test_throws ErrorException CLM.create_scale_l2g_lookup("bogus")
    end

    # =====================================================================
    # set_c2l_scale!
    # =====================================================================
    @testset "set_c2l_scale!" begin
        _, lun, col, _, bounds = make_test_hierarchy()

        # unity
        scale = zeros(bounds.endc)
        CLM.set_c2l_scale!(scale, bounds, "unity", col, lun)
        @test all(x -> x == 1.0, scale)

        # unsupported → error
        @test_throws ErrorException CLM.set_c2l_scale!(scale, bounds, "bogus", col, lun)
    end

    @testset "set_c2l_scale! urbanf" begin
        # Build a minimal urban hierarchy: 1 gridcell, 1 urban landunit, 5 urban columns
        nlun = 1; ncol = 5
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nlun)
        lun.itype[1] = CLM.ISTURB_HD
        lun.urbpoi[1] = true
        lun.canyon_hwr[1] = 2.0

        col = CLM.ColumnData()
        CLM.column_init!(col, ncol)
        for c in 1:ncol
            col.landunit[c] = 1
        end
        col.itype[1] = CLM.ICOL_ROOF
        col.itype[2] = CLM.ICOL_SUNWALL
        col.itype[3] = CLM.ICOL_SHADEWALL
        col.itype[4] = CLM.ICOL_ROAD_IMPERV
        col.itype[5] = CLM.ICOL_ROAD_PERV

        bounds = CLM.BoundsType(begc=1, endc=ncol)
        scale = zeros(ncol)

        # urbanf
        CLM.set_c2l_scale!(scale, bounds, "urbanf", col, lun)
        @test scale[1] == 1.0                  # roof
        @test scale[2] ≈ 3.0 * 2.0             # sunwall: 3*canyon_hwr
        @test scale[3] ≈ 3.0 * 2.0             # shadewall
        @test scale[4] == 3.0                   # road_imperv
        @test scale[5] == 3.0                   # road_perv

        # urbans
        CLM.set_c2l_scale!(scale, bounds, "urbans", col, lun)
        @test scale[1] == 1.0                   # roof
        @test scale[2] ≈ (3.0 * 2.0) / (2.0 * 2.0 + 1.0)  # sunwall
        @test scale[3] ≈ (3.0 * 2.0) / (2.0 * 2.0 + 1.0)  # shadewall
        @test scale[4] ≈ 3.0 / (2.0 * 2.0 + 1.0)           # road_imperv
        @test scale[5] ≈ 3.0 / (2.0 * 2.0 + 1.0)           # road_perv
    end

    # =====================================================================
    # p2c_1d!
    # =====================================================================
    @testset "p2c_1d! basic averaging" begin
        _, _, col, pch, bounds = make_test_hierarchy()

        # Patch values: 10, 20, 30, 40
        parr = [10.0, 20.0, 30.0, 40.0]
        carr = zeros(bounds.endc)

        CLM.p2c_1d!(carr, parr, bounds, "unity", pch)

        # Column 1 has patch 1 (wtcol=1.0): 10/1.0 = 10
        @test carr[1] ≈ 10.0
        # Column 2 has patch 2 (wtcol=1.0): 20/1.0 = 20
        @test carr[2] ≈ 20.0
        # Column 3 has patches 3,4 (wtcol=0.7,0.3): (30*0.7+40*0.3)/(0.7+0.3)=33
        @test carr[3] ≈ 33.0
    end

    @testset "p2c_1d! with spval" begin
        _, _, col, pch, bounds = make_test_hierarchy()

        # Patch 1 has spval → column 1 should remain spval
        parr = [CLM.SPVAL, 20.0, 30.0, 40.0]
        carr = zeros(bounds.endc)

        CLM.p2c_1d!(carr, parr, bounds, "unity", pch)

        @test carr[1] == CLM.SPVAL
        @test carr[2] ≈ 20.0
        @test carr[3] ≈ 33.0
    end

    @testset "p2c_1d! inactive patches skipped" begin
        _, _, col, pch, bounds = make_test_hierarchy()

        pch.active[4] = false  # deactivate patch 4

        parr = [10.0, 20.0, 30.0, 40.0]
        carr = zeros(bounds.endc)

        CLM.p2c_1d!(carr, parr, bounds, "unity", pch)

        # Column 3 now only has patch 3 (wtcol=0.7): 30*0.7/0.7 = 30
        @test carr[3] ≈ 30.0
    end

    # =====================================================================
    # p2c_2d!
    # =====================================================================
    @testset "p2c_2d! basic averaging" begin
        _, _, col, pch, bounds = make_test_hierarchy()
        num2d = 2

        parr = [10.0 100.0; 20.0 200.0; 30.0 300.0; 40.0 400.0]
        carr = zeros(bounds.endc, num2d)

        CLM.p2c_2d!(carr, parr, bounds, num2d, "unity", pch)

        @test carr[1, 1] ≈ 10.0
        @test carr[1, 2] ≈ 100.0
        @test carr[3, 1] ≈ 33.0
        @test carr[3, 2] ≈ 330.0
    end

    # =====================================================================
    # p2c_1d_filter!
    # =====================================================================
    @testset "p2c_1d_filter!" begin
        _, _, col, pch, bounds = make_test_hierarchy()

        patcharr = [10.0, 20.0, 30.0, 40.0]
        colarr = fill(999.0, bounds.endc)
        mask_c = BitVector([true, false, true])

        CLM.p2c_1d_filter!(colarr, patcharr, mask_c, col, pch)

        # Column 1 (masked in): patch 1, wtcol=1.0 → 10
        @test colarr[1] ≈ 10.0
        # Column 2 (masked out): untouched → 999
        @test colarr[2] == 999.0
        # Column 3 (masked in): patches 3,4, wtcol=0.7,0.3 → 30*0.7+40*0.3=33
        @test colarr[3] ≈ 33.0
    end

    # =====================================================================
    # p2c_2d_filter!
    # =====================================================================
    @testset "p2c_2d_filter!" begin
        _, _, col, pch, bounds = make_test_hierarchy()
        lev = 2

        patcharr = [10.0 100.0; 20.0 200.0; 30.0 300.0; 40.0 400.0]
        colarr = fill(999.0, bounds.endc, lev)
        mask_c = BitVector([true, false, true])

        CLM.p2c_2d_filter!(colarr, patcharr, lev, mask_c, col, pch)

        @test colarr[1, 1] ≈ 10.0
        @test colarr[1, 2] ≈ 100.0
        @test colarr[2, 1] == 999.0  # masked out
        @test colarr[3, 1] ≈ 33.0
        @test colarr[3, 2] ≈ 330.0
    end

    # =====================================================================
    # c2l_1d!
    # =====================================================================
    @testset "c2l_1d! basic averaging" begin
        _, lun, col, _, bounds = make_test_hierarchy()

        # Column values: 10, 20, 30
        carr = [10.0, 20.0, 30.0]
        larr = zeros(bounds.endl)

        CLM.c2l_1d!(larr, carr, bounds, "unity", col, lun)

        # Landunit 1: cols 1,2 with wtlunit 0.4, 0.6
        # (10*0.4 + 20*0.6) / (0.4+0.6) = 16
        @test larr[1] ≈ 16.0
        # Landunit 2: col 3 with wtlunit 1.0
        # 30*1.0 / 1.0 = 30
        @test larr[2] ≈ 30.0
    end

    @testset "c2l_1d! with spval input" begin
        _, lun, col, _, bounds = make_test_hierarchy()

        # Column 1 has spval
        carr = [CLM.SPVAL, 20.0, 30.0]
        larr = zeros(bounds.endl)

        CLM.c2l_1d!(larr, carr, bounds, "unity", col, lun)

        # Landunit 1: only col 2 contributes → 20*0.6/0.6 = 20
        @test larr[1] ≈ 20.0
        @test larr[2] ≈ 30.0
    end

    @testset "c2l_1d! with inactive columns" begin
        _, lun, col, _, bounds = make_test_hierarchy()

        col.active[1] = false  # deactivate column 1

        carr = [10.0, 20.0, 30.0]
        larr = zeros(bounds.endl)

        CLM.c2l_1d!(larr, carr, bounds, "unity", col, lun)

        # Landunit 1: only col 2 (active) → 20*0.6/0.6 = 20
        @test larr[1] ≈ 20.0
    end

    @testset "c2l_1d! include_inactive" begin
        _, lun, col, _, bounds = make_test_hierarchy()

        col.active[1] = false

        carr = [10.0, 20.0, 30.0]
        larr = zeros(bounds.endl)

        CLM.c2l_1d!(larr, carr, bounds, "unity", col, lun; include_inactive=true)

        # Both cols contribute: (10*0.4 + 20*0.6) / 1.0 = 16
        @test larr[1] ≈ 16.0
    end

    # =====================================================================
    # c2l_2d!
    # =====================================================================
    @testset "c2l_2d! basic averaging" begin
        _, lun, col, _, bounds = make_test_hierarchy()
        num2d = 2

        carr = [10.0 100.0; 20.0 200.0; 30.0 300.0]
        larr = zeros(bounds.endl, num2d)

        CLM.c2l_2d!(larr, carr, bounds, num2d, "unity", col, lun)

        @test larr[1, 1] ≈ 16.0
        @test larr[1, 2] ≈ 160.0
        @test larr[2, 1] ≈ 30.0
        @test larr[2, 2] ≈ 300.0
    end

    # =====================================================================
    # p2l_1d!
    # =====================================================================
    @testset "p2l_1d! basic averaging" begin
        _, lun, col, pch, bounds = make_test_hierarchy()

        # Patch values: 10, 20, 30, 40
        parr = [10.0, 20.0, 30.0, 40.0]
        larr = zeros(bounds.endl)

        CLM.p2l_1d!(larr, parr, bounds, "unity", "unity", pch, col, lun)

        # Landunit 1: patches 1,2 with wtlunit 0.4, 0.6
        # scale_p2c=1, scale_c2l=1
        # (10*1*1*0.4 + 20*1*1*0.6) / (0.4+0.6) = 16
        @test larr[1] ≈ 16.0
        # Landunit 2: patches 3,4 with wtlunit 0.7, 0.3
        # (30*1*1*0.7 + 40*1*1*0.3) / (0.7+0.3) = 33
        @test larr[2] ≈ 33.0
    end

    # =====================================================================
    # p2l_2d!
    # =====================================================================
    @testset "p2l_2d! basic averaging" begin
        _, lun, col, pch, bounds = make_test_hierarchy()
        num2d = 2

        parr = [10.0 100.0; 20.0 200.0; 30.0 300.0; 40.0 400.0]
        larr = zeros(bounds.endl, num2d)

        CLM.p2l_2d!(larr, parr, bounds, num2d, "unity", "unity", pch, col, lun)

        @test larr[1, 1] ≈ 16.0
        @test larr[1, 2] ≈ 160.0
        @test larr[2, 1] ≈ 33.0
        @test larr[2, 2] ≈ 330.0
    end

    # =====================================================================
    # l2g_1d!
    # =====================================================================
    @testset "l2g_1d! basic averaging" begin
        _, lun, _, _, bounds = make_test_hierarchy()

        larr = [5.0, 15.0]
        garr = zeros(bounds.endg)

        CLM.l2g_1d!(garr, larr, bounds, "unity", lun)

        # Gridcell 1: lun 1 (wt=0.6), lun 2 (wt=0.4)
        # (5*0.6 + 15*0.4) / (0.6+0.4) = 9
        @test garr[1] ≈ 9.0
    end

    @testset "l2g_1d! with natveg scale" begin
        _, lun, _, _, bounds = make_test_hierarchy()

        larr = [5.0, 15.0]
        garr = zeros(bounds.endg)

        CLM.l2g_1d!(garr, larr, bounds, "natveg", lun)

        # Only soil landunit (1) contributes: 5*1*0.6/0.6 = 5
        @test garr[1] ≈ 5.0
    end

    @testset "l2g_1d! inactive landunit skipped" begin
        _, lun, _, _, bounds = make_test_hierarchy()

        lun.active[2] = false

        larr = [5.0, 15.0]
        garr = zeros(bounds.endg)

        CLM.l2g_1d!(garr, larr, bounds, "unity", lun)

        # Only landunit 1 active: 5*0.6/0.6 = 5
        @test garr[1] ≈ 5.0
    end

    # =====================================================================
    # l2g_2d!
    # =====================================================================
    @testset "l2g_2d! basic averaging" begin
        _, lun, _, _, bounds = make_test_hierarchy()
        num2d = 2

        larr = [5.0 50.0; 15.0 150.0]
        garr = zeros(bounds.endg, num2d)

        CLM.l2g_2d!(garr, larr, bounds, num2d, "unity", lun)

        @test garr[1, 1] ≈ 9.0
        @test garr[1, 2] ≈ 90.0
    end

    # =====================================================================
    # c2g_1d!
    # =====================================================================
    @testset "c2g_1d! basic averaging" begin
        _, lun, col, _, bounds = make_test_hierarchy()

        carr = [10.0, 20.0, 30.0]
        garr = zeros(bounds.endg)

        CLM.c2g_1d!(garr, carr, bounds, "unity", "unity", col, lun)

        # Gridcell 1: cols 1,2,3 with wtgcell 0.24, 0.36, 0.4
        # (10*0.24 + 20*0.36 + 30*0.4) / (0.24+0.36+0.4) = (2.4+7.2+12)/1.0 = 21.6
        @test garr[1] ≈ 21.6
    end

    # =====================================================================
    # c2g_2d!
    # =====================================================================
    @testset "c2g_2d! basic averaging" begin
        _, lun, col, _, bounds = make_test_hierarchy()
        num2d = 2

        carr = [10.0 100.0; 20.0 200.0; 30.0 300.0]
        garr = zeros(bounds.endg, num2d)

        CLM.c2g_2d!(garr, carr, bounds, num2d, "unity", "unity", col, lun)

        @test garr[1, 1] ≈ 21.6
        @test garr[1, 2] ≈ 216.0
    end

    # =====================================================================
    # p2g_1d!
    # =====================================================================
    @testset "p2g_1d! basic averaging" begin
        _, lun, col, pch, bounds = make_test_hierarchy()

        parr = [10.0, 20.0, 30.0, 40.0]
        garr = zeros(bounds.endg)

        CLM.p2g_1d!(garr, parr, bounds, "unity", "unity", "unity", pch, col, lun)

        # Gridcell 1: patches 1,2,3,4 with wtgcell 0.24,0.36,0.28,0.12
        # (10*0.24 + 20*0.36 + 30*0.28 + 40*0.12) / (0.24+0.36+0.28+0.12)
        # = (2.4+7.2+8.4+4.8) / 1.0 = 22.8
        @test garr[1] ≈ 22.8
    end

    # =====================================================================
    # p2g_2d!
    # =====================================================================
    @testset "p2g_2d! basic averaging" begin
        _, lun, col, pch, bounds = make_test_hierarchy()
        num2d = 2

        parr = [10.0 100.0; 20.0 200.0; 30.0 300.0; 40.0 400.0]
        garr = zeros(bounds.endg, num2d)

        CLM.p2g_2d!(garr, parr, bounds, num2d, "unity", "unity", "unity", pch, col, lun)

        @test garr[1, 1] ≈ 22.8
        @test garr[1, 2] ≈ 228.0
    end

    # =====================================================================
    # All-spval input
    # =====================================================================
    @testset "all spval input → output stays spval" begin
        _, lun, col, pch, bounds = make_test_hierarchy()

        parr = fill(CLM.SPVAL, 4)
        carr = zeros(bounds.endc)
        CLM.p2c_1d!(carr, parr, bounds, "unity", pch)
        @test all(x -> x == CLM.SPVAL, carr)

        larr = zeros(bounds.endl)
        CLM.p2l_1d!(larr, parr, bounds, "unity", "unity", pch, col, lun)
        @test all(x -> x == CLM.SPVAL, larr)

        garr = zeros(bounds.endg)
        CLM.p2g_1d!(garr, parr, bounds, "unity", "unity", "unity", pch, col, lun)
        @test all(x -> x == CLM.SPVAL, garr)
    end

    # =====================================================================
    # Zero-weight patches
    # =====================================================================
    @testset "zero-weight patches skipped" begin
        _, lun, col, pch, bounds = make_test_hierarchy()

        pch.wtcol[1] = 0.0
        pch.wtlunit[1] = 0.0
        pch.wtgcell[1] = 0.0

        parr = [10.0, 20.0, 30.0, 40.0]
        carr = zeros(bounds.endc)

        CLM.p2c_1d!(carr, parr, bounds, "unity", pch)

        # Column 1: patch 1 has wtcol=0 → no contribution → spval
        @test carr[1] == CLM.SPVAL
        @test carr[2] ≈ 20.0
    end

    # =====================================================================
    # build_scale_l2g!
    # =====================================================================
    @testset "build_scale_l2g!" begin
        _, lun, _, _, bounds = make_test_hierarchy()

        scale_l2g = zeros(bounds.endl)
        CLM.build_scale_l2g!(scale_l2g, bounds, "unity", lun)
        @test all(x -> x == 1.0, scale_l2g)

        CLM.build_scale_l2g!(scale_l2g, bounds, "natveg", lun)
        @test scale_l2g[1] == 1.0                 # ISTSOIL
        @test scale_l2g[2] == CLM.SPVAL           # ISTCROP → excluded

        CLM.build_scale_l2g!(scale_l2g, bounds, "veg", lun)
        @test scale_l2g[1] == 1.0                 # ISTSOIL
        @test scale_l2g[2] == 1.0                 # ISTCROP → included
    end

end
