@testset "Subgrid Build" begin
    using Dates

    # ---- Test 1: add_landunit!, add_column!, add_patch! basic counters ----
    @testset "add_landunit!, add_column!, add_patch!" begin
        nlevsno_orig = CLM.varpar.nlevsno
        nlevmaxurbgrnd_orig = CLM.varpar.nlevmaxurbgrnd

        # Minimal setup: 1 gridcell, 2 landunits, 3 columns, 4 patches
        ng, nl, nc, np = 1, 2, 3, 4
        CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)
        CLM.varcon_init!()

        grc = CLM.GridcellData()
        lun = CLM.LandunitData()
        col = CLM.ColumnData()
        pch = CLM.PatchData()
        CLM.gridcell_init!(grc, ng)
        CLM.landunit_init!(lun, nl)
        CLM.column_init!(col, nc)
        CLM.patch_init!(pch, np)

        li = Ref(0)
        ci = Ref(0)
        pi = Ref(0)

        # Landunit 1 (soil), 1 column, 2 patches
        CLM.add_landunit!(lun, li, 1, CLM.ISTSOIL, 0.7)
        CLM.add_column!(col, lun, ci, li[], 1, 1.0)
        CLM.add_patch!(pch, col, lun, pi, ci[], 0, 0.8)
        CLM.add_patch!(pch, col, lun, pi, ci[], 1, 0.2)

        # Landunit 2 (lake), 2 columns, 2 patches
        CLM.add_landunit!(lun, li, 1, CLM.ISTDLAK, 0.3)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTDLAK, 0.5)
        CLM.add_patch!(pch, col, lun, pi, ci[], CLM.noveg, 1.0)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTDLAK, 0.5)
        CLM.add_patch!(pch, col, lun, pi, ci[], CLM.noveg, 1.0)

        @test li[] == 2
        @test ci[] == 3
        @test pi[] == 4
        @test lun.itype[1] == CLM.ISTSOIL
        @test lun.itype[2] == CLM.ISTDLAK
        @test lun.wtgcell[1] == 0.7
        @test lun.wtgcell[2] == 0.3
        @test col.landunit[1] == 1
        @test col.landunit[2] == 2
        @test col.landunit[3] == 2
        @test pch.column[1] == 1
        @test pch.column[2] == 1
        @test pch.column[3] == 2
        @test pch.column[4] == 3
    end

    # ---- Test 2: clm_ptrs_compdown! and clm_ptrs_check! ----
    @testset "clm_ptrs_compdown! and clm_ptrs_check!" begin
        ng, nl, nc, np = 1, 2, 3, 5
        CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)

        grc = CLM.GridcellData()
        lun = CLM.LandunitData()
        col = CLM.ColumnData()
        pch = CLM.PatchData()
        CLM.gridcell_init!(grc, ng)
        CLM.landunit_init!(lun, nl)
        CLM.column_init!(col, nc)
        CLM.patch_init!(pch, np)

        li = Ref(0)
        ci = Ref(0)
        pi = Ref(0)

        # Soil: 1 column, 3 patches
        CLM.add_landunit!(lun, li, 1, CLM.ISTSOIL, 0.7)
        CLM.add_column!(col, lun, ci, li[], 1, 1.0)
        CLM.add_patch!(pch, col, lun, pi, ci[], 0, 0.5)
        CLM.add_patch!(pch, col, lun, pi, ci[], 1, 0.3)
        CLM.add_patch!(pch, col, lun, pi, ci[], 2, 0.2)

        # Lake: 2 columns, 1 patch each
        CLM.add_landunit!(lun, li, 1, CLM.ISTDLAK, 0.3)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTDLAK, 0.5)
        CLM.add_patch!(pch, col, lun, pi, ci[], CLM.noveg, 1.0)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTDLAK, 0.5)
        CLM.add_patch!(pch, col, lun, pi, ci[], CLM.noveg, 1.0)

        bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl,
                                  begc=1, endc=nc, begp=1, endp=np,
                                  level=0, clump_index=1)

        # Fill down-pointers
        CLM.clm_ptrs_compdown!(bounds, grc, lun, col, pch)

        # Check: col.patchi/patchf
        @test col.patchi[1] == 1
        @test col.patchf[1] == 3
        @test col.npatches[1] == 3
        @test col.patchi[2] == 4
        @test col.patchf[2] == 4

        # Check: lun.coli/colf
        @test lun.coli[1] == 1
        @test lun.colf[1] == 1
        @test lun.coli[2] == 2
        @test lun.colf[2] == 3

        # Check: landunit_indices
        @test grc.landunit_indices[CLM.ISTSOIL, 1] == 1
        @test grc.landunit_indices[CLM.ISTDLAK, 1] == 2

        # Validate hierarchy
        CLM.clm_ptrs_check!(bounds, grc, lun, col, pch)
    end

    # ---- Test 3: compute_higher_order_weights! ----
    @testset "compute_higher_order_weights!" begin
        ng, nl, nc, np = 1, 1, 1, 2
        CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)

        grc = CLM.GridcellData()
        lun = CLM.LandunitData()
        col = CLM.ColumnData()
        pch = CLM.PatchData()
        CLM.gridcell_init!(grc, ng)
        CLM.landunit_init!(lun, nl)
        CLM.column_init!(col, nc)
        CLM.patch_init!(pch, np)

        li = Ref(0)
        ci = Ref(0)
        pi = Ref(0)

        CLM.add_landunit!(lun, li, 1, CLM.ISTSOIL, 0.6)
        CLM.add_column!(col, lun, ci, li[], 1, 1.0)
        CLM.add_patch!(pch, col, lun, pi, ci[], 0, 0.7)
        CLM.add_patch!(pch, col, lun, pi, ci[], 1, 0.3)

        bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl,
                                  begc=1, endc=nc, begp=1, endp=np,
                                  level=0, clump_index=1)

        CLM.compute_higher_order_weights!(bounds, col, lun, pch)

        @test col.wtgcell[1] ≈ 1.0 * 0.6   # col.wtlunit * lun.wtgcell
        @test pch.wtlunit[1] ≈ 0.7 * 1.0    # pch.wtcol * col.wtlunit
        @test pch.wtgcell[1] ≈ 0.7 * 0.6    # pch.wtcol * col.wtgcell
        @test pch.wtgcell[2] ≈ 0.3 * 0.6
    end

    # ---- Test 4: set_active! ----
    @testset "set_active!" begin
        ng, nl, nc, np = 1, 2, 2, 2
        CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)

        lun = CLM.LandunitData()
        col = CLM.ColumnData()
        pch = CLM.PatchData()
        CLM.landunit_init!(lun, nl)
        CLM.column_init!(col, nc)
        CLM.patch_init!(pch, np)

        li = Ref(0)
        ci = Ref(0)
        pi = Ref(0)

        # Soil with weight > 0
        CLM.add_landunit!(lun, li, 1, CLM.ISTSOIL, 0.8)
        CLM.add_column!(col, lun, ci, li[], 1, 1.0)
        CLM.add_patch!(pch, col, lun, pi, ci[], 0, 1.0)

        # Wetland with weight 0
        CLM.add_landunit!(lun, li, 1, CLM.ISTWET, 0.0)
        CLM.add_column!(col, lun, ci, li[], CLM.ISTWET, 1.0)
        CLM.add_patch!(pch, col, lun, pi, ci[], CLM.noveg, 1.0)

        bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl,
                                  begc=1, endc=nc, begp=1, endp=np,
                                  level=0, clump_index=1)

        old_all_active = CLM.varctl.all_active
        CLM.varctl.all_active = false
        CLM.set_active!(bounds, lun, col, pch)

        # Soil is always active
        @test lun.active[1] == true
        @test col.active[1] == true
        @test pch.active[1] == true

        # Wetland with zero weight: not active
        @test lun.active[2] == false
        @test col.active[2] == false
        @test pch.active[2] == false

        CLM.varctl.all_active = old_all_active
    end

    # ---- Test 5: TimeManager ----
    @testset "TimeManager" begin
        tm = CLM.TimeManager(
            start_date=DateTime(2000, 1, 1),
            current_date=DateTime(2000, 1, 1),
            dtime=1800
        )
        @test CLM.get_nstep(tm) == 0
        @test CLM.get_curr_calday(tm) ≈ 1.0

        CLM.advance_timestep!(tm)
        @test CLM.get_nstep(tm) == 1
        (yr, mon, d, tod) = CLM.get_curr_date(tm)
        @test yr == 2000
        @test mon == 1
        @test d == 1
        @test tod == 1800
        @test CLM.is_beg_curr_day(tm) == false

        # Advance to midnight next day
        for _ in 1:47
            CLM.advance_timestep!(tm)
        end
        @test CLM.is_beg_curr_day(tm) == true
        @test CLM.is_beg_curr_year(tm) == false
    end

    # ---- Test 6: surfrd_utils ----
    @testset "surfrd_utils" begin
        # renormalize!
        arr = [0.3 0.2 0.5; 0.0 0.0 0.0]
        CLM.renormalize!(arr; target=1.0)
        @test sum(arr[1, :]) ≈ 1.0
        @test arr[2, 1] ≈ 0.0  # zero row stays zero
    end

    # ---- Test 7: SurfaceInputData count_subgrid_elements ----
    @testset "count_subgrid_elements" begin
        ng = 1
        CLM.varpar_init!(CLM.varpar, 1, 15, 2, 5)
        CLM.varcon_init!()
        old_create_crop = CLM.varctl.create_crop_landunit
        old_rzwu = CLM.varctl.run_zero_weight_urban
        CLM.varctl.create_crop_landunit = false
        CLM.varctl.run_zero_weight_urban = false

        surf = CLM.SurfaceInputData()
        surf.wt_lunit = zeros(ng, CLM.MAX_LUNIT)
        surf.wt_lunit[1, CLM.ISTSOIL] = 0.8
        surf.wt_lunit[1, CLM.ISTDLAK] = 0.2
        surf.wt_nat_patch = zeros(ng, 2)
        surf.wt_nat_patch[1, 1] = 0.6
        surf.wt_nat_patch[1, 2] = 0.4
        surf.wt_cft = zeros(ng, 0)
        surf.wt_glc_mec = zeros(ng, 0)

        (nl, nc, np) = CLM.count_subgrid_elements(surf, ng)

        # Soil: 1 landunit, 1 column, 2 patches
        # Lake: 1 landunit, 1 column, 1 patch (always created)
        @test nl == 2
        @test nc == 2
        @test np == 3

        CLM.varctl.create_crop_landunit = old_create_crop
        CLM.varctl.run_zero_weight_urban = old_rzwu
    end
end
