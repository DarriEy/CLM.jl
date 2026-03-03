@testset "TopoData" begin

    # Helper to create minimal ColumnData and LandunitData for testing
    function make_test_col_lun(ncols::Int, nlun::Int)
        col = CLM.ColumnData()
        CLM.column_init!(col, ncols)

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nlun)

        return col, lun
    end

    @testset "TopoData default construction" begin
        t = CLM.TopoData()
        @test length(t.topo_col) == 0
        @test length(t.needs_downscaling_col) == 0
    end

    @testset "topo_init_allocate!" begin
        t = CLM.TopoData()
        CLM.topo_init_allocate!(t, 5)
        @test length(t.topo_col) == 5
        @test all(isnan, t.topo_col)
        @test length(t.needs_downscaling_col) == 5
        @test all(x -> x == false, t.needs_downscaling_col)
    end

    @testset "topo_init_history! is no-op" begin
        t = CLM.TopoData()
        CLM.topo_init_allocate!(t, 3)
        @test CLM.topo_init_history!(t) === nothing
    end

    @testset "topo_init_cold! non-ice columns" begin
        ncols = 3
        nlun = 1
        col, lun = make_test_col_lun(ncols, nlun)

        # Set up: all columns on landunit 1 which is soil type
        lun.itype[1] = CLM.ISTSOIL
        for c in 1:ncols
            col.landunit[c] = 1
            col.gridcell[c] = 1
            col.itype[c] = 1
            col.is_hillslope_column[c] = false
        end

        t = CLM.TopoData()
        CLM.topo_init_allocate!(t, ncols)
        CLM.topo_init_cold!(t, ncols, col, lun)

        @test all(x -> x == 0.0, t.topo_col)
        @test all(x -> x == false, t.needs_downscaling_col)
    end

    @testset "topo_init_cold! ice columns with topo_glc_mec" begin
        ncols = 2
        nlun = 1
        col, lun = make_test_col_lun(ncols, nlun)

        lun.itype[1] = CLM.ISTICE
        for c in 1:ncols
            col.landunit[c] = 1
            col.gridcell[c] = 1
            col.itype[c] = c  # ice class 1 and 2
            col.is_hillslope_column[c] = false
        end

        # topo_glc_mec: gridcells x ice_classes
        topo_glc_mec = [100.0 200.0; 300.0 400.0]  # 2 gridcells x 2 ice classes

        t = CLM.TopoData()
        CLM.topo_init_allocate!(t, ncols)
        CLM.topo_init_cold!(t, ncols, col, lun; topo_glc_mec = topo_glc_mec)

        @test t.topo_col[1] == 100.0  # gridcell 1, ice class 1
        @test t.topo_col[2] == 200.0  # gridcell 1, ice class 2
        @test t.needs_downscaling_col[1] == true
        @test t.needs_downscaling_col[2] == true
    end

    @testset "topo_init_cold! ice columns without topo_glc_mec" begin
        ncols = 1
        nlun = 1
        col, lun = make_test_col_lun(ncols, nlun)

        lun.itype[1] = CLM.ISTICE
        col.landunit[1] = 1
        col.gridcell[1] = 1
        col.itype[1] = 1
        col.is_hillslope_column[1] = false

        t = CLM.TopoData()
        CLM.topo_init_allocate!(t, ncols)
        CLM.topo_init_cold!(t, ncols, col, lun)

        @test t.topo_col[1] == 0.0
        @test t.needs_downscaling_col[1] == true
    end

    @testset "topo_init_cold! hillslope columns with downscaling" begin
        ncols = 2
        nlun = 1
        col, lun = make_test_col_lun(ncols, nlun)

        lun.itype[1] = CLM.ISTSOIL
        for c in 1:ncols
            col.landunit[c] = 1
            col.gridcell[c] = 1
            col.itype[c] = 1
        end
        col.is_hillslope_column[1] = true
        col.is_hillslope_column[2] = false
        col.hill_elev[1] = 500.0

        t = CLM.TopoData()
        CLM.topo_init_allocate!(t, ncols)
        CLM.topo_init_cold!(t, ncols, col, lun;
                            use_hillslope = true,
                            downscale_hillslope_meteorology = true)

        @test t.topo_col[1] == 500.0
        @test t.needs_downscaling_col[1] == true
        @test t.topo_col[2] == 0.0
        @test t.needs_downscaling_col[2] == false
    end

    @testset "topo_init! combines allocate + cold start" begin
        ncols = 3
        nlun = 1
        col, lun = make_test_col_lun(ncols, nlun)

        lun.itype[1] = CLM.ISTSOIL
        for c in 1:ncols
            col.landunit[c] = 1
            col.gridcell[c] = 1
            col.itype[c] = 1
            col.is_hillslope_column[c] = false
        end

        t = CLM.TopoData()
        CLM.topo_init!(t, ncols, col, lun)

        @test length(t.topo_col) == ncols
        @test all(x -> x == 0.0, t.topo_col)
        @test all(x -> x == false, t.needs_downscaling_col)
    end

    @testset "topo_restart! is no-op" begin
        t = CLM.TopoData()
        @test CLM.topo_restart!(t) === nothing
    end

    @testset "topo_update! basic — no downscaling, no hillslope" begin
        ncols = 3
        nlun = 1
        col, lun = make_test_col_lun(ncols, nlun)

        lun.itype[1] = CLM.ISTSOIL
        lun.coli[1] = 1
        lun.colf[1] = ncols
        for c in 1:ncols
            col.landunit[c] = 1
            col.gridcell[c] = 1
            col.itype[c] = 1
            col.is_hillslope_column[c] = false
        end

        t = CLM.TopoData()
        CLM.topo_init!(t, ncols, col, lun)

        atm_topo = [1500.0]  # one gridcell
        CLM.topo_update!(t, ncols, col, lun, 1:1, atm_topo)

        # All columns should get atmosphere topo since no downscaling needed
        @test all(x -> x == 1500.0, t.topo_col)
        @test all(x -> x == false, t.needs_downscaling_col)
    end

    @testset "topo_update! with ice_cols_need_downscaling callback" begin
        ncols = 3
        nlun = 1
        col, lun = make_test_col_lun(ncols, nlun)

        lun.itype[1] = CLM.ISTSOIL
        lun.coli[1] = 1
        lun.colf[1] = ncols
        for c in 1:ncols
            col.landunit[c] = 1
            col.gridcell[c] = 1
            col.itype[c] = 1
            col.is_hillslope_column[c] = false
        end

        t = CLM.TopoData()
        CLM.topo_init!(t, ncols, col, lun)

        # Callback that marks column 2 as needing downscaling
        ice_cb! = (mask_ice, needs_ds) -> begin
            needs_ds[2] = true
        end

        # Callback that sets topo for downscaled columns
        glc2lnd_cb! = (topo_col, needs_ds) -> begin
            for c in eachindex(topo_col)
                if needs_ds[c]
                    topo_col[c] = 2000.0
                end
            end
        end

        mask_ice = BitVector([false, true, false])
        atm_topo = [1500.0]

        CLM.topo_update!(t, ncols, col, lun, 1:1, atm_topo;
                         mask_ice = mask_ice,
                         ice_cols_need_downscaling! = ice_cb!,
                         update_glc2lnd_topo! = glc2lnd_cb!)

        @test t.topo_col[1] == 1500.0   # not downscaled
        @test t.topo_col[2] == 2000.0   # set by glc2lnd callback
        @test t.topo_col[3] == 1500.0   # not downscaled
        @test t.needs_downscaling_col[2] == true
    end

    @testset "topo_update! with hillslope meteorology downscaling" begin
        ncols = 3
        nlun = 1
        col, lun = make_test_col_lun(ncols, nlun)

        lun.itype[1] = CLM.ISTSOIL
        lun.coli[1] = 1
        lun.colf[1] = ncols
        for c in 1:ncols
            col.landunit[c] = 1
            col.gridcell[c] = 1
            col.itype[c] = 1
            col.is_hillslope_column[c] = true
            col.hill_area[c] = 100.0
        end
        col.hill_elev[1] = 100.0
        col.hill_elev[2] = 200.0
        col.hill_elev[3] = 300.0

        t = CLM.TopoData()
        CLM.topo_init!(t, ncols, col, lun;
                        use_hillslope = true,
                        downscale_hillslope_meteorology = true)

        atm_topo = [1000.0]

        CLM.topo_update!(t, ncols, col, lun, 1:1, atm_topo;
                         use_hillslope = true,
                         downscale_hillslope_meteorology = true)

        # Mean hillslope elevation = (100*100 + 200*100 + 300*100) / (100+100+100) = 200
        # Column 1: atm_topo + (100 - 200) = 1000 - 100 = 900
        # Column 2: atm_topo + (200 - 200) = 1000 + 0 = 1000
        # Column 3: atm_topo + (300 - 200) = 1000 + 100 = 1100
        @test t.topo_col[1] ≈ 900.0
        @test t.topo_col[2] ≈ 1000.0
        @test t.topo_col[3] ≈ 1100.0
        @test all(t.needs_downscaling_col)
    end

    @testset "topo_downscale_filter" begin
        t = CLM.TopoData()
        CLM.topo_init_allocate!(t, 4)
        t.needs_downscaling_col[1] = false
        t.needs_downscaling_col[2] = true
        t.needs_downscaling_col[3] = false
        t.needs_downscaling_col[4] = true

        mask = CLM.topo_downscale_filter(t, 4)
        @test mask isa BitVector
        @test mask == BitVector([false, true, false, true])
    end

    @testset "topo_clean!" begin
        t = CLM.TopoData()
        CLM.topo_init_allocate!(t, 5)
        @test length(t.topo_col) == 5

        CLM.topo_clean!(t)
        @test length(t.topo_col) == 0
        @test length(t.needs_downscaling_col) == 0
    end

    @testset "field mutability" begin
        t = CLM.TopoData()
        CLM.topo_init_allocate!(t, 3)

        t.topo_col[1] = 500.0
        @test t.topo_col[1] == 500.0

        t.needs_downscaling_col[2] = true
        @test t.needs_downscaling_col[2] == true
    end

    @testset "re-init overwrites previous state" begin
        ncols = 3
        nlun = 1
        col, lun = make_test_col_lun(ncols, nlun)

        lun.itype[1] = CLM.ISTSOIL
        for c in 1:ncols
            col.landunit[c] = 1
            col.gridcell[c] = 1
            col.itype[c] = 1
            col.is_hillslope_column[c] = false
        end

        t = CLM.TopoData()
        CLM.topo_init!(t, ncols, col, lun)
        t.topo_col[1] = 999.0

        CLM.topo_init!(t, ncols, col, lun)
        @test t.topo_col[1] == 0.0
    end

    @testset "topo_update! with update_glc_classes callback" begin
        ncols = 2
        nlun = 1
        col, lun = make_test_col_lun(ncols, nlun)

        lun.itype[1] = CLM.ISTSOIL
        lun.coli[1] = 1
        lun.colf[1] = ncols
        for c in 1:ncols
            col.landunit[c] = 1
            col.gridcell[c] = 1
            col.itype[c] = 1
            col.is_hillslope_column[c] = false
        end

        t = CLM.TopoData()
        CLM.topo_init!(t, ncols, col, lun)

        classes_updated = Ref(false)
        glc_classes_cb! = (topo_col) -> begin
            classes_updated[] = true
        end

        atm_topo = [1500.0]
        CLM.topo_update!(t, ncols, col, lun, 1:1, atm_topo;
                         update_glc_classes! = glc_classes_cb!)

        @test classes_updated[] == true
    end

end
