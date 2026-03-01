@testset "Filters" begin

    # Helper: create a BoundsType for a single-clump domain
    function make_bounds(; ngridcells, nlandunits, ncols, npatches)
        CLM.BoundsType(
            begg=1, endg=ngridcells,
            begl=1, endl=nlandunits,
            begc=1, endc=ncols,
            begp=1, endp=npatches,
            begCohort=0, endCohort=0,
            level=CLM.BOUNDS_LEVEL_CLUMP,
            clump_index=1
        )
    end

    @testset "alloc_filters!" begin
        filt = CLM.ClumpFilter()
        CLM.alloc_filters!(filt, 5, 8, 3)

        # Column masks have correct size
        @test length(filt.allc) == 5
        @test length(filt.lakec) == 5
        @test length(filt.soilc) == 5
        @test length(filt.icec) == 5

        # Patch masks have correct size
        @test length(filt.soilp) == 8
        @test length(filt.lakep) == 8
        @test length(filt.nolakeurbanp) == 8

        # Landunit masks have correct size
        @test length(filt.urbanl) == 3
        @test length(filt.nourbanl) == 3

        # All masks start as false
        @test count(filt.allc) == 0
        @test count(filt.soilp) == 0
        @test count(filt.urbanl) == 0
    end

    @testset "alloc_all_filters!" begin
        CLM.alloc_all_filters!(4, 6, 2)
        @test length(CLM.clump_filter.allc) == 4
        @test length(CLM.clump_filter_inactive_and_active.allc) == 4
        @test length(CLM.clump_filter.soilp) == 6
        @test length(CLM.clump_filter.urbanl) == 2
    end

    @testset "set_filters_one_group! basic domain" begin
        # Initialize varpar for ColumnData
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        # Domain: 1 gridcell, 3 landunits (soil, lake, urban), 5 columns, 5 patches
        ngridcells = 1
        nlandunits = 3
        ncols = 5
        npatches = 5

        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ngridcells)
        grc.active[1] = true

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nlandunits)
        lun.gridcell .= 1
        lun.itype .= [CLM.ISTSOIL, CLM.ISTDLAK, CLM.ISTURB_TBD]
        lun.lakpoi .= [false, true, false]
        lun.urbpoi .= [false, false, true]
        lun.active .= true

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, ncols)
        col_data.landunit .= [1, 1, 2, 3, 3]
        col_data.gridcell .= 1
        col_data.itype .= [CLM.ISTSOIL, CLM.ISTSOIL, CLM.ISTDLAK, CLM.ICOL_ROOF, CLM.ICOL_ROAD_PERV]
        col_data.active .= true
        col_data.hydrologically_active .= [true, true, false, false, true]

        pch = CLM.PatchData()
        CLM.patch_init!(pch, npatches)
        pch.landunit .= [1, 1, 2, 3, 3]
        pch.column .= [1, 2, 3, 4, 5]
        pch.gridcell .= 1
        pch.itype .= [1, 2, 0, 0, 0]
        pch.active .= true

        bounds = make_bounds(ngridcells=1, nlandunits=3, ncols=5, npatches=5)

        # Save and restore varctl flags
        saved_use_cn = CLM.varctl.use_cn
        saved_use_fates = CLM.varctl.use_fates
        saved_use_fates_bgc = CLM.varctl.use_fates_bgc
        CLM.varctl.use_cn = false
        CLM.varctl.use_fates = false
        CLM.varctl.use_fates_bgc = false

        filt = CLM.ClumpFilter()
        CLM.alloc_filters!(filt, ncols, npatches, nlandunits)
        CLM.set_filters_one_group!(filt, bounds, false,
                                    col_data, lun, pch, grc)

        # Column filters
        @test filt.allc == BitVector([1,1,1,1,1])
        @test filt.lakec == BitVector([0,0,1,0,0])
        @test filt.nolakec == BitVector([1,1,0,1,1])
        @test filt.soilc == BitVector([1,1,0,0,0])
        @test filt.bgc_soilc == BitVector([0,0,0,0,0])  # use_cn=false
        @test filt.hydrologyc == BitVector([1,1,0,0,1])
        @test filt.urbanc == BitVector([0,0,0,1,1])
        @test filt.nourbanc == BitVector([1,1,1,0,0])
        @test filt.icec == BitVector([0,0,0,0,0])

        # Patch filters
        @test filt.lakep == BitVector([0,0,1,0,0])
        @test filt.nolakep == BitVector([1,1,0,1,1])
        @test filt.nolakeurbanp == BitVector([1,1,0,0,0])
        @test filt.soilp == BitVector([1,1,0,0,0])
        @test filt.bgc_vegp == BitVector([0,0,0,0,0])  # use_cn=false
        @test filt.urbanp == BitVector([0,0,0,1,1])
        @test filt.nourbanp == BitVector([1,1,1,0,0])

        # Crop filters (use_fates=false, all PFT types < NPCROPMIN)
        @test filt.pcropp == BitVector([0,0,0,0,0])
        @test filt.soilnopcropp == BitVector([1,1,0,0,0])

        # Landunit filters
        @test filt.urbanl == BitVector([0,0,1])
        @test filt.nourbanl == BitVector([1,1,0])

        # Restore flags
        CLM.varctl.use_cn = saved_use_cn
        CLM.varctl.use_fates = saved_use_fates
        CLM.varctl.use_fates_bgc = saved_use_fates_bgc
    end

    @testset "set_filters_one_group! active/inactive" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        ngridcells = 1
        nlandunits = 2
        ncols = 3
        npatches = 3

        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ngridcells)
        grc.active[1] = true

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nlandunits)
        lun.gridcell .= 1
        lun.itype .= [CLM.ISTSOIL, CLM.ISTSOIL]
        lun.lakpoi .= false
        lun.urbpoi .= false
        lun.active .= [true, false]

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, ncols)
        col_data.landunit .= [1, 1, 2]
        col_data.gridcell .= 1
        col_data.active .= [true, true, false]
        col_data.hydrologically_active .= [true, true, true]

        pch = CLM.PatchData()
        CLM.patch_init!(pch, npatches)
        pch.landunit .= [1, 1, 2]
        pch.column .= [1, 2, 3]
        pch.gridcell .= 1
        pch.itype .= [1, 2, 3]
        pch.active .= [true, true, false]

        bounds = make_bounds(ngridcells=1, nlandunits=2, ncols=3, npatches=3)

        saved_use_fates = CLM.varctl.use_fates
        CLM.varctl.use_fates = false

        # Active only (include_inactive=false)
        filt_active = CLM.ClumpFilter()
        CLM.alloc_filters!(filt_active, ncols, npatches, nlandunits)
        CLM.set_filters_one_group!(filt_active, bounds, false,
                                    col_data, lun, pch, grc)

        @test filt_active.allc == BitVector([1,1,0])
        @test filt_active.soilc == BitVector([1,1,0])
        @test filt_active.soilp == BitVector([1,1,0])
        @test filt_active.nourbanl == BitVector([1,0])

        # Include inactive (include_inactive=true)
        filt_all = CLM.ClumpFilter()
        CLM.alloc_filters!(filt_all, ncols, npatches, nlandunits)
        CLM.set_filters_one_group!(filt_all, bounds, true,
                                    col_data, lun, pch, grc)

        @test filt_all.allc == BitVector([1,1,1])
        @test filt_all.soilc == BitVector([1,1,1])
        @test filt_all.soilp == BitVector([1,1,1])
        @test filt_all.nourbanl == BitVector([1,1])

        CLM.varctl.use_fates = saved_use_fates
    end

    @testset "set_filters_one_group! BGC flags" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        ngridcells = 1
        nlandunits = 2
        ncols = 2
        npatches = 2

        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ngridcells)
        grc.active[1] = true

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nlandunits)
        lun.gridcell .= 1
        lun.itype .= [CLM.ISTSOIL, CLM.ISTDLAK]
        lun.lakpoi .= [false, true]
        lun.urbpoi .= false
        lun.active .= true

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, ncols)
        col_data.landunit .= [1, 2]
        col_data.gridcell .= 1
        col_data.active .= true
        col_data.hydrologically_active .= true

        pch = CLM.PatchData()
        CLM.patch_init!(pch, npatches)
        pch.landunit .= [1, 2]
        pch.column .= [1, 2]
        pch.gridcell .= 1
        pch.itype .= [1, 0]
        pch.active .= true

        bounds = make_bounds(ngridcells=1, nlandunits=2, ncols=2, npatches=2)

        # Test with use_cn=true
        saved_use_cn = CLM.varctl.use_cn
        saved_use_fates = CLM.varctl.use_fates
        saved_use_fates_bgc = CLM.varctl.use_fates_bgc
        CLM.varctl.use_cn = true
        CLM.varctl.use_fates = false
        CLM.varctl.use_fates_bgc = false

        filt = CLM.ClumpFilter()
        CLM.alloc_filters!(filt, ncols, npatches, nlandunits)
        CLM.set_filters_one_group!(filt, bounds, false,
                                    col_data, lun, pch, grc)

        @test filt.bgc_soilc == BitVector([1,0])  # soil only, not lake
        @test filt.bgc_vegp == BitVector([1,0])    # soil only, not lake

        # Test with use_fates_bgc=true (no use_cn)
        CLM.varctl.use_cn = false
        CLM.varctl.use_fates_bgc = true

        filt2 = CLM.ClumpFilter()
        CLM.alloc_filters!(filt2, ncols, npatches, nlandunits)
        CLM.set_filters_one_group!(filt2, bounds, false,
                                    col_data, lun, pch, grc)

        @test filt2.bgc_soilc == BitVector([1,0])  # fates_bgc triggers bgc_soilc
        @test filt2.bgc_vegp == BitVector([0,0])    # bgc_vegp requires use_cn

        CLM.varctl.use_cn = saved_use_cn
        CLM.varctl.use_fates = saved_use_fates
        CLM.varctl.use_fates_bgc = saved_use_fates_bgc
    end

    @testset "crop filters" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        ngridcells = 1
        nlandunits = 2
        ncols = 3
        npatches = 3

        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ngridcells)
        grc.active[1] = true

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nlandunits)
        lun.gridcell .= 1
        lun.itype .= [CLM.ISTSOIL, CLM.ISTCROP]
        lun.lakpoi .= false
        lun.urbpoi .= false
        lun.active .= true

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, ncols)
        col_data.landunit .= [1, 2, 2]
        col_data.gridcell .= 1
        col_data.active .= true
        col_data.hydrologically_active .= true

        pch = CLM.PatchData()
        CLM.patch_init!(pch, npatches)
        pch.landunit .= [1, 2, 2]
        pch.column .= [1, 2, 3]
        pch.gridcell .= 1
        pch.itype .= [1, 15, 17]  # natural veg, generic crop, prognostic crop
        pch.active .= true

        bounds = make_bounds(ngridcells=1, nlandunits=2, ncols=3, npatches=3)

        saved_use_fates = CLM.varctl.use_fates
        CLM.varctl.use_fates = false

        filt = CLM.ClumpFilter()
        CLM.alloc_filters!(filt, ncols, npatches, nlandunits)
        CLM.set_filters_one_group!(filt, bounds, false,
                                    col_data, lun, pch, grc)

        # patch 1: itype=1 < 17, landunit=soil → soilnopcropp
        # patch 2: itype=15 < 17, landunit=crop → soilnopcropp
        # patch 3: itype=17 >= 17 → pcropp
        @test filt.pcropp == BitVector([0, 0, 1])
        @test filt.soilnopcropp == BitVector([1, 1, 0])

        # Both soil and crop landunits count as soil
        @test filt.soilc == BitVector([1, 1, 1])
        @test filt.soilp == BitVector([1, 1, 1])

        # Test with use_fates=true: crop filters should be empty
        CLM.varctl.use_fates = true

        filt2 = CLM.ClumpFilter()
        CLM.alloc_filters!(filt2, ncols, npatches, nlandunits)
        CLM.set_filters_one_group!(filt2, bounds, false,
                                    col_data, lun, pch, grc)

        @test filt2.pcropp == BitVector([0, 0, 0])
        @test filt2.soilnopcropp == BitVector([0, 0, 0])

        CLM.varctl.use_fates = saved_use_fates
    end

    @testset "ice and SMB filters" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        ngridcells = 1
        nlandunits = 2
        ncols = 2
        npatches = 2

        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ngridcells)
        grc.active[1] = true

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nlandunits)
        lun.gridcell .= 1
        lun.itype .= [CLM.ISTSOIL, CLM.ISTICE]
        lun.lakpoi .= false
        lun.urbpoi .= false
        lun.active .= true

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, ncols)
        col_data.landunit .= [1, 2]
        col_data.gridcell .= 1
        col_data.active .= true
        col_data.hydrologically_active .= [true, true]

        pch = CLM.PatchData()
        CLM.patch_init!(pch, npatches)
        pch.landunit .= [1, 2]
        pch.column .= [1, 2]
        pch.gridcell .= 1
        pch.itype .= [1, 0]
        pch.active .= true

        bounds = make_bounds(ngridcells=1, nlandunits=2, ncols=2, npatches=2)

        saved_use_fates = CLM.varctl.use_fates
        CLM.varctl.use_fates = false

        filt = CLM.ClumpFilter()
        CLM.alloc_filters!(filt, ncols, npatches, nlandunits)

        # With melt_replaced_by_ice = [true] for gridcell 1
        CLM.set_filters_one_group!(filt, bounds, false,
                                    col_data, lun, pch, grc;
                                    melt_replaced_by_ice=[true])

        @test filt.icec == BitVector([0, 1])
        # SMB: both soil and ice qualify when melt_replaced_by_ice is true
        @test filt.do_smb_c == BitVector([1, 1])

        # Without melt_replaced_by_ice
        filt2 = CLM.ClumpFilter()
        CLM.alloc_filters!(filt2, ncols, npatches, nlandunits)
        CLM.set_filters_one_group!(filt2, bounds, false,
                                    col_data, lun, pch, grc)

        @test filt2.icec == BitVector([0, 1])
        @test filt2.do_smb_c == BitVector([0, 0])  # empty when no melt vector

        # With melt_replaced_by_ice = [false]
        filt3 = CLM.ClumpFilter()
        CLM.alloc_filters!(filt3, ncols, npatches, nlandunits)
        CLM.set_filters_one_group!(filt3, bounds, false,
                                    col_data, lun, pch, grc;
                                    melt_replaced_by_ice=[false])

        @test filt3.do_smb_c == BitVector([0, 0])

        CLM.varctl.use_fates = saved_use_fates
    end

    @testset "set_exposedvegp_filter!" begin
        npatches = 5
        ncols = 5
        nlandunits = 3

        filt = CLM.ClumpFilter()
        CLM.alloc_filters!(filt, ncols, npatches, nlandunits)

        # Simulate: patches 1,2 are nolakeurban; 3,4,5 are lake/urban
        filt.nolakeurbanp = BitVector([1, 1, 0, 0, 0])

        bounds = make_bounds(ngridcells=1, nlandunits=3, ncols=5, npatches=5)

        frac_veg_nosno = [1, 0, 1, 0, 0]
        CLM.set_exposedvegp_filter!(filt, bounds, frac_veg_nosno)

        # Patch 1: nolakeurban + frac>0 → exposed
        # Patch 2: nolakeurban + frac=0 → noexposed
        # Patches 3-5: not nolakeurban → neither
        @test filt.exposedvegp == BitVector([1,0,0,0,0])
        @test filt.noexposedvegp == BitVector([0,1,0,0,0])

        # Re-call with different frac_veg_nosno to verify masks are cleared
        frac_veg_nosno2 = [0, 1, 0, 0, 0]
        CLM.set_exposedvegp_filter!(filt, bounds, frac_veg_nosno2)

        @test filt.exposedvegp == BitVector([0,1,0,0,0])
        @test filt.noexposedvegp == BitVector([1,0,0,0,0])
    end

    @testset "set_filters! sets both filter groups" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        ngridcells = 1
        nlandunits = 1
        ncols = 2
        npatches = 2

        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ngridcells)
        grc.active[1] = true

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nlandunits)
        lun.gridcell .= 1
        lun.itype .= [CLM.ISTSOIL]
        lun.lakpoi .= false
        lun.urbpoi .= false
        lun.active .= true

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, ncols)
        col_data.landunit .= [1, 1]
        col_data.gridcell .= 1
        col_data.active .= [true, false]  # column 2 inactive
        col_data.hydrologically_active .= true

        pch = CLM.PatchData()
        CLM.patch_init!(pch, npatches)
        pch.landunit .= [1, 1]
        pch.column .= [1, 2]
        pch.gridcell .= 1
        pch.itype .= [1, 2]
        pch.active .= [true, false]

        CLM.alloc_all_filters!(ncols, npatches, nlandunits)

        bounds = make_bounds(ngridcells=1, nlandunits=1, ncols=2, npatches=2)

        saved_use_fates = CLM.varctl.use_fates
        CLM.varctl.use_fates = false

        CLM.set_filters!(bounds, col_data, lun, pch, grc)

        # Active-only filter: column 2 excluded
        @test CLM.clump_filter.allc == BitVector([1,0])
        @test CLM.clump_filter.soilc == BitVector([1,0])

        # Inactive+active filter: column 2 included
        @test CLM.clump_filter_inactive_and_active.allc == BitVector([1,1])
        @test CLM.clump_filter_inactive_and_active.soilc == BitVector([1,1])

        CLM.varctl.use_fates = saved_use_fates
    end

    @testset "NPCROPMIN constant" begin
        @test CLM.NPCROPMIN == 17
    end
end
