@testset "Hillslope Hydrology Module" begin
    # ------------------------------------------------------------------
    # Tests for HillslopeHydrologyMod port.
    # Verifies:
    #   1. hillslope_properties_init!
    #   2. check_aquifer_layer!
    #   3. hillslope_soil_thickness_profile! (set_lowland_upland)
    #   4. hillslope_soil_thickness_profile_linear!
    #   5. hillslope_set_lowland_upland_pfts!
    #   6. hillslope_pft_from_file!
    #   7. hillslope_stream_outflow!
    #   8. hillslope_update_stream_water!
    # ------------------------------------------------------------------

    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    CLM.varcon_init!()
    nlevsoi  = CLM.varpar.nlevsoi
    nlevgrnd = CLM.varpar.nlevgrnd
    nlevsno  = CLM.varpar.nlevsno

    # ------------------------------------------------------------------
    # 1. hillslope_properties_init!
    # ------------------------------------------------------------------
    @testset "hillslope_properties_init!" begin
        CLM.hillslope_properties_init!(
            pft_distribution_method = "Standard",
            soil_profile_method = "Uniform")
        @test CLM.hillslope_config[].pft_distribution_method == CLM.PFT_STANDARD
        @test CLM.hillslope_config[].soil_profile_method == CLM.SOIL_PROFILE_UNIFORM

        CLM.hillslope_properties_init!(
            pft_distribution_method = "FromFile",
            soil_profile_method = "SetLowlandUpland")
        @test CLM.hillslope_config[].pft_distribution_method == CLM.PFT_FROM_FILE
        @test CLM.hillslope_config[].soil_profile_method == CLM.SOIL_PROFILE_SET_LOWLAND_UPLAND

        CLM.hillslope_properties_init!(
            pft_distribution_method = "PftLowlandUpland",
            soil_profile_method = "Linear")
        @test CLM.hillslope_config[].pft_distribution_method == CLM.PFT_LOWLAND_UPLAND
        @test CLM.hillslope_config[].soil_profile_method == CLM.SOIL_PROFILE_LINEAR

        # Reset to default
        CLM.hillslope_properties_init!()

        # Invalid method should throw
        @test_throws ErrorException CLM.hillslope_properties_init!(
            pft_distribution_method = "InvalidMethod")
        @test_throws ErrorException CLM.hillslope_properties_init!(
            soil_profile_method = "InvalidMethod")
    end

    # ------------------------------------------------------------------
    # 2. check_aquifer_layer!
    # ------------------------------------------------------------------
    @testset "check_aquifer_layer!" begin
        # Should not throw when one or both are false
        CLM.check_aquifer_layer!(false, false)
        CLM.check_aquifer_layer!(true, false)
        CLM.check_aquifer_layer!(false, true)

        # Should throw when both are true
        @test_throws ErrorException CLM.check_aquifer_layer!(true, true)
    end

    # ------------------------------------------------------------------
    # 3. hillslope_soil_thickness_profile! (set_lowland_upland)
    # ------------------------------------------------------------------
    @testset "hillslope_soil_thickness_profile! set_lowland_upland" begin
        nc = 3
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        col.is_hillslope_column .= true
        col.active .= true

        # Column 1: upland (has downhill neighbor)
        col.cold[1] = 2
        # Column 2: lowland (no downhill neighbor)
        col.cold[2] = CLM.ISPVAL
        # Column 3: upland
        col.cold[3] = 2

        # Set nbedrock to initial value
        col.nbedrock .= nlevsoi

        CLM.hillslope_soil_thickness_profile!(
            col, 1:nc,
            soil_profile_method = CLM.SOIL_PROFILE_SET_LOWLAND_UPLAND,
            nlevsoi = nlevsoi,
            soil_depth_lowland = 8.0,
            soil_depth_upland = 2.0)

        # Upland columns should have shallower bedrock than lowland
        @test col.nbedrock[1] < col.nbedrock[2]  # upland < lowland
        @test col.nbedrock[3] == col.nbedrock[1]  # both upland
        @test col.nbedrock[2] > 0  # lowland nbedrock set
        @test col.nbedrock[1] > 0  # upland nbedrock set
    end

    # ------------------------------------------------------------------
    # 4. hillslope_soil_thickness_profile_linear!
    # ------------------------------------------------------------------
    @testset "hillslope_soil_thickness_profile_linear!" begin
        nc = 4
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        col.is_hillslope_column .= true
        col.active .= true
        col.nbedrock .= nlevsoi

        # Set hill_distance: linearly increasing
        col.hill_distance[1] = 0.0    # closest to stream (lowland)
        col.hill_distance[2] = 100.0
        col.hill_distance[3] = 200.0
        col.hill_distance[4] = 300.0  # farthest (upland)

        CLM.hillslope_soil_thickness_profile_linear!(
            col, 1:nc,
            nlevsoi = nlevsoi,
            soil_depth_lowland = 8.0,
            soil_depth_upland = 2.0)

        # Lowland (dist=0) should have deepest soil, upland (dist=300) shallowest
        @test col.nbedrock[1] >= col.nbedrock[4]
        # Middle columns should be in between
        @test col.nbedrock[2] >= col.nbedrock[3]
        @test col.nbedrock[3] >= col.nbedrock[4]
    end

    # ------------------------------------------------------------------
    # 5. hillslope_set_lowland_upland_pfts!
    # ------------------------------------------------------------------
    @testset "hillslope_set_lowland_upland_pfts!" begin
        nc = 2
        np = 2

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)

        col.is_hillslope_column .= true
        col.patchi .= [1, 2]
        col.patchf .= [1, 2]

        # Column 1: lowland (no downhill neighbor)
        col.cold[1] = CLM.ISPVAL
        # Column 2: upland
        col.cold[2] = 1

        pch.itype .= [0, 0]
        pch.mxy .= [0, 0]

        lowland_ivt = 7   # broadleaf deciduous tree
        upland_ivt = 13   # c3 non-arctic grass

        CLM.hillslope_set_lowland_upland_pfts!(
            col, pch, 1:nc,
            lowland_ivt = lowland_ivt,
            upland_ivt = upland_ivt)

        @test pch.itype[1] == lowland_ivt  # lowland column
        @test pch.itype[2] == upland_ivt   # upland column
        @test pch.mxy[1] == lowland_ivt + 1
        @test pch.mxy[2] == upland_ivt + 1
    end

    # ------------------------------------------------------------------
    # 6. hillslope_pft_from_file!
    # ------------------------------------------------------------------
    @testset "hillslope_pft_from_file!" begin
        nc = 3
        np = 3

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        pch = CLM.PatchData()
        CLM.patch_init!(pch, np)

        col.is_hillslope_column .= [true, true, false]
        col.patchi .= [1, 2, 3]
        col.patchf .= [1, 2, 3]

        pch.itype .= [0, 0, 0]
        pch.mxy .= [0, 0, 0]

        col_pftndx = [5, 10, 0]

        CLM.hillslope_pft_from_file!(
            col, pch, 1:nc, col_pftndx)

        @test pch.itype[1] == 5   # hillslope column 1
        @test pch.itype[2] == 10  # hillslope column 2
        @test pch.itype[3] == 0   # non-hillslope, unchanged
    end

    # ------------------------------------------------------------------
    # 7. hillslope_stream_outflow!
    # ------------------------------------------------------------------
    @testset "hillslope_stream_outflow!" begin
        nl = 2
        dtime = 1800.0
        bounds_l = 1:nl

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)

        lun.itype .= [CLM.ISTSOIL, CLM.ISTDLAK]
        lun.active .= [true, true]

        # Set stream channel properties for landunit 1
        lun.stream_channel_length[1] = 1000.0  # m
        lun.stream_channel_width[1]  = 5.0     # m
        lun.stream_channel_depth[1]  = 1.0     # m
        lun.stream_channel_slope[1]  = 0.01    # m/m
        lun.stream_channel_number[1] = 2.0     # 2 channel reaches

        # Landunit 2: lake, no stream
        lun.stream_channel_length[2] = 0.0
        lun.stream_channel_width[2]  = 0.0

        stream_water_volume    = [500.0, 0.0]   # m3
        volumetric_streamflow  = [0.0, 0.0]

        CLM.hillslope_stream_outflow!(
            stream_water_volume, volumetric_streamflow,
            lun, bounds_l, dtime)

        # Landunit 1: should have positive streamflow
        @test volumetric_streamflow[1] > 0.0
        # Streamflow should not exceed available water / dtime
        @test volumetric_streamflow[1] <= stream_water_volume[1] / dtime

        # Landunit 2: lake, no streamflow
        @test volumetric_streamflow[2] == 0.0
    end

    @testset "hillslope_stream_outflow! zero water" begin
        nl = 1
        dtime = 1800.0

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)

        lun.itype .= [CLM.ISTSOIL]
        lun.active .= [true]
        lun.stream_channel_length[1] = 1000.0
        lun.stream_channel_width[1]  = 5.0
        lun.stream_channel_depth[1]  = 1.0
        lun.stream_channel_slope[1]  = 0.01
        lun.stream_channel_number[1] = 1.0

        stream_water_volume    = [0.0]
        volumetric_streamflow  = [0.0]

        CLM.hillslope_stream_outflow!(
            stream_water_volume, volumetric_streamflow,
            lun, 1:nl, dtime)

        # Zero water → zero or non-negative streamflow
        @test volumetric_streamflow[1] >= 0.0
    end

    @testset "hillslope_stream_outflow! overbank" begin
        nl = 1
        dtime = 1800.0

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)

        lun.itype .= [CLM.ISTSOIL]
        lun.active .= [true]
        lun.stream_channel_length[1] = 1000.0
        lun.stream_channel_width[1]  = 5.0
        lun.stream_channel_depth[1]  = 1.0
        lun.stream_channel_slope[1]  = 0.01
        lun.stream_channel_number[1] = 1.0

        # Large water volume → stream_depth > channel_depth → overbank
        stream_water_volume    = [50000.0]  # m3
        volumetric_streamflow  = [0.0]

        CLM.hillslope_stream_outflow!(
            stream_water_volume, volumetric_streamflow,
            lun, 1:nl, dtime)

        @test volumetric_streamflow[1] > 0.0
        @test volumetric_streamflow[1] <= stream_water_volume[1] / dtime
    end

    # ------------------------------------------------------------------
    # 8. hillslope_update_stream_water!
    # ------------------------------------------------------------------
    @testset "hillslope_update_stream_water!" begin
        nl = 1
        nc = 2
        ng = 1
        dtime = 1800.0

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ng)

        lun.itype .= [CLM.ISTSOIL]
        lun.active .= [true]
        lun.gridcell .= [1]
        lun.wtgcell .= [0.5]
        lun.coli .= [1]
        lun.colf .= [nc]
        lun.stream_channel_length[1] = 1000.0
        lun.stream_channel_width[1]  = 5.0
        lun.stream_channel_depth[1]  = 1.0
        lun.stream_channel_slope[1]  = 0.01
        lun.stream_channel_number[1] = 1.0

        grc.area[1] = 100.0  # km^2

        col.is_hillslope_column .= true
        col.active .= true
        col.gridcell .= [1, 1]
        col.landunit .= [1, 1]
        col.wtgcell .= [0.25, 0.25]

        stream_water_volume   = [1000.0]
        volumetric_streamflow = [0.1]     # m3/s
        stream_water_depth    = [0.0]

        qflx_drain         = [0.001, 0.002]    # mm/s
        qflx_drain_perched = [0.0005, 0.001]   # mm/s
        qflx_surf          = [0.0001, 0.0002]  # mm/s

        initial_volume = stream_water_volume[1]

        CLM.hillslope_update_stream_water!(
            stream_water_volume, volumetric_streamflow, stream_water_depth,
            qflx_drain, qflx_drain_perched, qflx_surf,
            col, lun, grc,
            1:nl, dtime)

        # Stream water volume should have changed
        @test stream_water_volume[1] != initial_volume
        # Stream water depth should be positive
        @test stream_water_depth[1] >= 0.0

        # Compute expected additions
        col1_area = grc.area[1] * 1.0e6 * col.wtgcell[1]
        col2_area = grc.area[1] * 1.0e6 * col.wtgcell[2]
        inflow_vol = ((qflx_drain[1] + qflx_drain_perched[1] + qflx_surf[1]) * 1.0e-3 * col1_area +
                      (qflx_drain[2] + qflx_drain_perched[2] + qflx_surf[2]) * 1.0e-3 * col2_area) * dtime
        expected_volume = initial_volume + inflow_vol - volumetric_streamflow[1] * dtime
        @test stream_water_volume[1] ≈ expected_volume
    end

    @testset "hillslope_update_stream_water! negative drainage" begin
        nl = 1
        nc = 1
        ng = 1
        dtime = 1800.0

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ng)

        lun.itype .= [CLM.ISTSOIL]
        lun.active .= [true]
        lun.gridcell .= [1]
        lun.wtgcell .= [1.0]
        lun.coli .= [1]
        lun.colf .= [1]
        lun.stream_channel_length[1] = 100.0
        lun.stream_channel_width[1]  = 2.0
        lun.stream_channel_depth[1]  = 0.5
        lun.stream_channel_slope[1]  = 0.01
        lun.stream_channel_number[1] = 1.0

        grc.area[1] = 1.0

        col.is_hillslope_column .= true
        col.active .= true
        col.gridcell .= [1]
        col.landunit .= [1]
        col.wtgcell .= [1.0]

        # Large streamflow, small volume → should go negative then be clamped
        stream_water_volume   = [10.0]
        volumetric_streamflow = [100.0]   # much larger than volume/dtime
        stream_water_depth    = [0.0]

        # Very small drainage
        qflx_drain         = [0.0]
        qflx_drain_perched = [0.0]
        qflx_surf          = [0.0]

        CLM.hillslope_update_stream_water!(
            stream_water_volume, volumetric_streamflow, stream_water_depth,
            qflx_drain, qflx_drain_perched, qflx_surf,
            col, lun, grc,
            1:nl, dtime)

        # Volume should be clamped to 0
        @test stream_water_volume[1] == 0.0
        # Streamflow adjusted for negative
        @test stream_water_depth[1] == 0.0
    end

    # ------------------------------------------------------------------
    # 9. set_hillslope_soil_thickness! (integration)
    # ------------------------------------------------------------------
    @testset "set_hillslope_soil_thickness!" begin
        nc = 2
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        col.is_hillslope_column .= true
        col.active .= true
        col.cold[1] = CLM.ISPVAL  # lowland
        col.cold[2] = 1           # upland
        col.nbedrock .= nlevsoi

        CLM.set_hillslope_soil_thickness!(
            col, CLM.LandunitData(), 1:nc,
            soil_profile_method = CLM.SOIL_PROFILE_SET_LOWLAND_UPLAND,
            nlevsoi = nlevsoi,
            soil_depth_lowland = 8.0,
            soil_depth_upland = 2.0)

        # Lowland deeper than upland
        @test col.nbedrock[1] > col.nbedrock[2]
    end

    # ------------------------------------------------------------------
    # 10. Constants check
    # ------------------------------------------------------------------
    @testset "Module constants" begin
        @test CLM.STREAMFLOW_MANNING == 0
        @test CLM.PFT_STANDARD == 0
        @test CLM.PFT_FROM_FILE == 1
        @test CLM.PFT_UNIFORM_DOMINANT_PFT == 2
        @test CLM.PFT_LOWLAND_DOMINANT_PFT == 3
        @test CLM.PFT_LOWLAND_UPLAND == 4
        @test CLM.SOIL_PROFILE_UNIFORM == 0
        @test CLM.SOIL_PROFILE_FROM_FILE == 1
        @test CLM.SOIL_PROFILE_SET_LOWLAND_UPLAND == 2
        @test CLM.SOIL_PROFILE_LINEAR == 3
    end
end
