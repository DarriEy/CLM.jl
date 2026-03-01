@testset "ColumnData" begin

    @testset "default construction" begin
        col = CLM.ColumnData()
        @test length(col.gridcell) == 0
        @test length(col.landunit) == 0
        @test length(col.wtgcell) == 0
        @test length(col.active) == 0
        @test length(col.snl) == 0
        @test size(col.dz) == (0, 0)
        @test size(col.z) == (0, 0)
        @test size(col.zi) == (0, 0)
        @test length(col.hill_elev) == 0
        @test length(col.urbpoi) == 0
    end

    @testset "column_init!" begin
        # Initialize varpar so dimension parameters are available
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        nc = 10
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        nlevsno = CLM.varpar.nlevsno
        nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd
        nlevlak = CLM.varpar.nlevlak

        # --- Check sizes: g/l/c/p hierarchy ---
        @test length(col.gridcell) == nc
        @test length(col.wtgcell) == nc
        @test length(col.landunit) == nc
        @test length(col.wtlunit) == nc
        @test length(col.patchi) == nc
        @test length(col.patchf) == nc
        @test length(col.npatches) == nc

        # --- Check sizes: topological mapping ---
        @test length(col.itype) == nc
        @test length(col.lun_itype) == nc
        @test length(col.active) == nc
        @test length(col.type_is_dynamic) == nc
        @test length(col.is_fates) == nc

        # --- Check sizes: topography ---
        @test length(col.micro_sigma) == nc
        @test length(col.topo_slope) == nc
        @test length(col.topo_std) == nc

        # --- Check sizes: vertical levels ---
        @test length(col.snl) == nc
        @test size(col.dz) == (nc, nlevsno + nlevmaxurbgrnd)
        @test size(col.z) == (nc, nlevsno + nlevmaxurbgrnd)
        @test size(col.zi) == (nc, nlevsno + nlevmaxurbgrnd + 1)
        @test length(col.zii) == nc
        @test size(col.dz_lake) == (nc, nlevlak)
        @test size(col.z_lake) == (nc, nlevlak)
        @test length(col.lakedepth) == nc
        @test length(col.nbedrock) == nc

        # --- Check sizes: hillslope hydrology ---
        @test length(col.col_ndx) == nc
        @test length(col.colu) == nc
        @test length(col.cold) == nc
        @test length(col.hillslope_ndx) == nc
        @test length(col.hill_elev) == nc
        @test length(col.hill_slope) == nc
        @test length(col.hill_area) == nc
        @test length(col.hill_width) == nc
        @test length(col.hill_distance) == nc
        @test length(col.hill_aspect) == nc

        # --- Check sizes: other characteristics ---
        @test length(col.is_hillslope_column) == nc
        @test length(col.hydrologically_active) == nc
        @test length(col.urbpoi) == nc
        @test size(col.levgrnd_class) == (nc, nlevmaxurbgrnd)

        # --- Check integer sentinel values ---
        @test all(==(CLM.ISPVAL), col.gridcell)
        @test all(==(CLM.ISPVAL), col.landunit)
        @test all(==(CLM.ISPVAL), col.patchi)
        @test all(==(CLM.ISPVAL), col.patchf)
        @test all(==(CLM.ISPVAL), col.npatches)
        @test all(==(CLM.ISPVAL), col.itype)
        @test all(==(CLM.ISPVAL), col.lun_itype)
        @test all(==(CLM.ISPVAL), col.snl)
        @test all(==(CLM.ISPVAL), col.nbedrock)
        @test all(==(CLM.ISPVAL), col.col_ndx)
        @test all(==(CLM.ISPVAL), col.colu)
        @test all(==(CLM.ISPVAL), col.cold)
        @test all(==(CLM.ISPVAL), col.hillslope_ndx)
        @test all(==(CLM.ISPVAL), col.levgrnd_class)

        # --- Check NaN initialization for real fields ---
        @test all(isnan, col.wtgcell)
        @test all(isnan, col.wtlunit)
        @test all(isnan, col.micro_sigma)
        @test all(isnan, col.topo_slope)
        @test all(isnan, col.topo_std)
        @test all(isnan, col.dz)
        @test all(isnan, col.z)
        @test all(isnan, col.zi)
        @test all(isnan, col.zii)
        @test all(isnan, col.dz_lake)
        @test all(isnan, col.z_lake)

        # --- Check SPVAL initialization ---
        @test all(==(CLM.SPVAL), col.lakedepth)
        @test all(==(CLM.SPVAL), col.hill_elev)
        @test all(==(CLM.SPVAL), col.hill_slope)
        @test all(==(CLM.SPVAL), col.hill_area)
        @test all(==(CLM.SPVAL), col.hill_width)
        @test all(==(CLM.SPVAL), col.hill_distance)
        @test all(==(CLM.SPVAL), col.hill_aspect)

        # --- Check boolean defaults ---
        @test all(.!col.active)
        @test all(.!col.type_is_dynamic)
        @test all(.!col.is_fates)
        @test all(.!col.is_hillslope_column)
        @test all(.!col.hydrologically_active)
        @test all(.!col.urbpoi)
    end

    @testset "column_clean!" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        col = CLM.ColumnData()
        CLM.column_init!(col, 10)
        CLM.column_clean!(col)

        # All vectors should be empty
        @test length(col.gridcell) == 0
        @test length(col.wtgcell) == 0
        @test length(col.landunit) == 0
        @test length(col.wtlunit) == 0
        @test length(col.patchi) == 0
        @test length(col.patchf) == 0
        @test length(col.npatches) == 0
        @test length(col.itype) == 0
        @test length(col.lun_itype) == 0
        @test length(col.active) == 0
        @test length(col.is_fates) == 0
        @test length(col.type_is_dynamic) == 0
        @test length(col.snl) == 0
        @test size(col.dz) == (0, 0)
        @test size(col.z) == (0, 0)
        @test size(col.zi) == (0, 0)
        @test length(col.zii) == 0
        @test length(col.lakedepth) == 0
        @test size(col.dz_lake) == (0, 0)
        @test size(col.z_lake) == (0, 0)
        @test length(col.micro_sigma) == 0
        @test length(col.topo_slope) == 0
        @test length(col.topo_std) == 0
        @test length(col.nbedrock) == 0
        @test size(col.levgrnd_class) == (0, 0)
        @test length(col.is_hillslope_column) == 0
        @test length(col.hydrologically_active) == 0
        @test length(col.col_ndx) == 0
        @test length(col.colu) == 0
        @test length(col.cold) == 0
        @test length(col.hillslope_ndx) == 0
        @test length(col.hill_elev) == 0
        @test length(col.hill_slope) == 0
        @test length(col.hill_area) == 0
        @test length(col.hill_width) == 0
        @test length(col.hill_distance) == 0
        @test length(col.hill_aspect) == 0
        @test length(col.urbpoi) == 0
    end

    @testset "field mutability" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        col = CLM.ColumnData()
        CLM.column_init!(col, 3)

        # Write to fields and verify
        col.gridcell[1] = 42
        col.wtgcell[2] = 0.75
        col.landunit[1] = 10
        col.wtlunit[1] = 0.5
        col.patchi[1] = 1
        col.patchf[1] = 8
        col.npatches[1] = 8
        col.itype[1] = CLM.ISTSOIL
        col.lun_itype[1] = CLM.ISTSOIL
        col.active[1] = true
        col.type_is_dynamic[1] = true
        col.is_fates[2] = true
        col.snl[1] = -3
        col.dz[1, 1] = 0.05
        col.z[1, 1] = 0.025
        col.zi[1, 1] = 0.0
        col.zii[1] = 100.0
        col.lakedepth[2] = 50.0
        col.dz_lake[1, 1] = 1.0
        col.z_lake[1, 1] = 0.5
        col.nbedrock[1] = 15
        col.micro_sigma[1] = 0.5
        col.topo_slope[1] = 0.01
        col.topo_std[1] = 50.0
        col.col_ndx[1] = 1
        col.colu[1] = 2
        col.cold[1] = 3
        col.hillslope_ndx[1] = 1
        col.hill_elev[1] = 100.0
        col.hill_slope[1] = 0.05
        col.hill_area[1] = 1000.0
        col.hill_width[1] = 50.0
        col.hill_distance[1] = 200.0
        col.hill_aspect[1] = 1.57
        col.is_hillslope_column[1] = true
        col.hydrologically_active[1] = true
        col.urbpoi[3] = true
        col.levgrnd_class[1, 1] = 1

        @test col.gridcell[1] == 42
        @test col.wtgcell[2] == 0.75
        @test col.landunit[1] == 10
        @test col.wtlunit[1] == 0.5
        @test col.patchi[1] == 1
        @test col.patchf[1] == 8
        @test col.npatches[1] == 8
        @test col.itype[1] == CLM.ISTSOIL
        @test col.lun_itype[1] == CLM.ISTSOIL
        @test col.active[1] == true
        @test col.type_is_dynamic[1] == true
        @test col.is_fates[2] == true
        @test col.snl[1] == -3
        @test col.dz[1, 1] == 0.05
        @test col.z[1, 1] == 0.025
        @test col.zi[1, 1] == 0.0
        @test col.zii[1] == 100.0
        @test col.lakedepth[2] == 50.0
        @test col.dz_lake[1, 1] == 1.0
        @test col.z_lake[1, 1] == 0.5
        @test col.nbedrock[1] == 15
        @test col.micro_sigma[1] == 0.5
        @test col.topo_slope[1] == 0.01
        @test col.topo_std[1] == 50.0
        @test col.col_ndx[1] == 1
        @test col.colu[1] == 2
        @test col.cold[1] == 3
        @test col.hillslope_ndx[1] == 1
        @test col.hill_elev[1] == 100.0
        @test col.hill_slope[1] == 0.05
        @test col.hill_area[1] == 1000.0
        @test col.hill_width[1] == 50.0
        @test col.hill_distance[1] == 200.0
        @test col.hill_aspect[1] == 1.57
        @test col.is_hillslope_column[1] == true
        @test col.hydrologically_active[1] == true
        @test col.urbpoi[3] == true
        @test col.levgrnd_class[1, 1] == 1
    end

    @testset "column_update_itype!" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        col = CLM.ColumnData()
        CLM.column_init!(col, 3)

        # Set up a dynamic column on a soil landunit
        col.type_is_dynamic[1] = true
        col.itype[1] = CLM.ISTSOIL
        col.lun_itype[1] = CLM.ISTSOIL

        # Update should succeed and set hydrologically_active
        CLM.column_update_itype!(col, 1, CLM.ISTCROP)
        @test col.itype[1] == CLM.ISTCROP
        @test col.hydrologically_active[1] == true

        # Set up a dynamic column on an urban landunit
        col.type_is_dynamic[2] = true
        col.lun_itype[2] = CLM.ISTURB_HD

        # Update to impervious road — not hydrologically active
        CLM.column_update_itype!(col, 2, CLM.ICOL_ROAD_IMPERV)
        @test col.itype[2] == CLM.ICOL_ROAD_IMPERV
        @test col.hydrologically_active[2] == false

        # Update to pervious road — hydrologically active
        CLM.column_update_itype!(col, 2, CLM.ICOL_ROAD_PERV)
        @test col.itype[2] == CLM.ICOL_ROAD_PERV
        @test col.hydrologically_active[2] == true

        # Attempt to update a non-dynamic column should throw an error
        col.type_is_dynamic[3] = false
        @test_throws ErrorException CLM.column_update_itype!(col, 3, CLM.ISTSOIL)
    end

    @testset "re-init overwrites previous state" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        col = CLM.ColumnData()
        CLM.column_init!(col, 3)
        col.wtgcell[1] = 999.0

        # Re-init with different size
        CLM.column_init!(col, 7)
        @test length(col.wtgcell) == 7
        @test all(isnan, col.wtgcell)  # old value gone
    end

end
