@testset "LandunitData" begin

    @testset "default construction" begin
        lun = CLM.LandunitData()
        @test length(lun.gridcell) == 0
        @test length(lun.wtgcell) == 0
        @test length(lun.active) == 0
        @test length(lun.canyon_hwr) == 0
        @test length(lun.stream_channel_depth) == 0
    end

    @testset "landunit_init!" begin
        nl = 5
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)

        # Check sizes — hierarchy fields
        @test length(lun.gridcell) == nl
        @test length(lun.wtgcell) == nl
        @test length(lun.coli) == nl
        @test length(lun.colf) == nl
        @test length(lun.ncolumns) == nl
        @test length(lun.nhillslopes) == nl
        @test length(lun.patchi) == nl
        @test length(lun.patchf) == nl
        @test length(lun.npatches) == nl

        # Check sizes — topological mapping
        @test length(lun.itype) == nl
        @test length(lun.ifspecial) == nl
        @test length(lun.lakpoi) == nl
        @test length(lun.urbpoi) == nl
        @test length(lun.glcpoi) == nl
        @test length(lun.active) == nl

        # Check sizes — urban properties
        @test length(lun.canyon_hwr) == nl
        @test length(lun.wtroad_perv) == nl
        @test length(lun.wtlunit_roof) == nl
        @test length(lun.ht_roof) == nl
        @test length(lun.z_0_town) == nl
        @test length(lun.z_d_town) == nl

        # Check sizes — hillslope variables
        @test length(lun.stream_channel_depth) == nl
        @test length(lun.stream_channel_width) == nl
        @test length(lun.stream_channel_length) == nl
        @test length(lun.stream_channel_slope) == nl
        @test length(lun.stream_channel_number) == nl

        # Check integer sentinel values
        @test all(==(CLM.ISPVAL), lun.gridcell)
        @test all(==(CLM.ISPVAL), lun.coli)
        @test all(==(CLM.ISPVAL), lun.colf)
        @test all(==(CLM.ISPVAL), lun.ncolumns)
        @test all(==(CLM.ISPVAL), lun.nhillslopes)
        @test all(==(CLM.ISPVAL), lun.patchi)
        @test all(==(CLM.ISPVAL), lun.patchf)
        @test all(==(CLM.ISPVAL), lun.npatches)
        @test all(==(CLM.ISPVAL), lun.itype)

        # Check NaN initialization for real fields
        @test all(isnan, lun.wtgcell)
        @test all(isnan, lun.canyon_hwr)
        @test all(isnan, lun.wtroad_perv)
        @test all(isnan, lun.wtlunit_roof)
        @test all(isnan, lun.ht_roof)
        @test all(isnan, lun.z_0_town)
        @test all(isnan, lun.z_d_town)
        @test all(isnan, lun.stream_channel_depth)
        @test all(isnan, lun.stream_channel_width)
        @test all(isnan, lun.stream_channel_length)
        @test all(isnan, lun.stream_channel_slope)
        @test all(isnan, lun.stream_channel_number)

        # Check boolean defaults
        @test all(.!lun.ifspecial)
        @test all(.!lun.lakpoi)
        @test all(.!lun.urbpoi)
        @test all(.!lun.glcpoi)
    end

    @testset "landunit_clean!" begin
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 10)
        CLM.landunit_clean!(lun)

        @test length(lun.gridcell) == 0
        @test length(lun.wtgcell) == 0
        @test length(lun.coli) == 0
        @test length(lun.colf) == 0
        @test length(lun.ncolumns) == 0
        @test length(lun.nhillslopes) == 0
        @test length(lun.patchi) == 0
        @test length(lun.patchf) == 0
        @test length(lun.npatches) == 0
        @test length(lun.itype) == 0
        @test length(lun.ifspecial) == 0
        @test length(lun.lakpoi) == 0
        @test length(lun.urbpoi) == 0
        @test length(lun.glcpoi) == 0
        @test length(lun.active) == 0
        @test length(lun.canyon_hwr) == 0
        @test length(lun.wtroad_perv) == 0
        @test length(lun.wtlunit_roof) == 0
        @test length(lun.ht_roof) == 0
        @test length(lun.z_0_town) == 0
        @test length(lun.z_d_town) == 0
        @test length(lun.stream_channel_depth) == 0
        @test length(lun.stream_channel_width) == 0
        @test length(lun.stream_channel_length) == 0
        @test length(lun.stream_channel_slope) == 0
        @test length(lun.stream_channel_number) == 0
    end

    @testset "field mutability" begin
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 3)

        # Write to fields and verify
        lun.gridcell[1] = 42
        lun.wtgcell[2] = 0.75
        lun.coli[1] = 10
        lun.colf[1] = 15
        lun.ncolumns[1] = 6
        lun.nhillslopes[1] = 2
        lun.patchi[1] = 1
        lun.patchf[1] = 8
        lun.npatches[1] = 8
        lun.itype[1] = 1
        lun.ifspecial[2] = true
        lun.lakpoi[3] = true
        lun.urbpoi[1] = true
        lun.glcpoi[2] = true
        lun.active[1] = true
        lun.canyon_hwr[1] = 0.5
        lun.wtroad_perv[1] = 0.3
        lun.wtlunit_roof[1] = 0.4
        lun.ht_roof[1] = 10.0
        lun.z_0_town[1] = 0.1
        lun.z_d_town[1] = 5.0
        lun.stream_channel_depth[1] = 2.0
        lun.stream_channel_width[1] = 10.0
        lun.stream_channel_length[1] = 500.0
        lun.stream_channel_slope[1] = 0.01
        lun.stream_channel_number[1] = 3.0

        @test lun.gridcell[1] == 42
        @test lun.wtgcell[2] == 0.75
        @test lun.coli[1] == 10
        @test lun.colf[1] == 15
        @test lun.ncolumns[1] == 6
        @test lun.nhillslopes[1] == 2
        @test lun.patchi[1] == 1
        @test lun.patchf[1] == 8
        @test lun.npatches[1] == 8
        @test lun.itype[1] == 1
        @test lun.ifspecial[2] == true
        @test lun.lakpoi[3] == true
        @test lun.urbpoi[1] == true
        @test lun.glcpoi[2] == true
        @test lun.active[1] == true
        @test lun.canyon_hwr[1] == 0.5
        @test lun.wtroad_perv[1] == 0.3
        @test lun.wtlunit_roof[1] == 0.4
        @test lun.ht_roof[1] == 10.0
        @test lun.z_0_town[1] == 0.1
        @test lun.z_d_town[1] == 5.0
        @test lun.stream_channel_depth[1] == 2.0
        @test lun.stream_channel_width[1] == 10.0
        @test lun.stream_channel_length[1] == 500.0
        @test lun.stream_channel_slope[1] == 0.01
        @test lun.stream_channel_number[1] == 3.0
    end

    @testset "re-init overwrites previous state" begin
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 3)
        lun.wtgcell[1] = 999.0

        # Re-init with different size
        CLM.landunit_init!(lun, 7)
        @test length(lun.wtgcell) == 7
        @test all(isnan, lun.wtgcell)  # old value gone
    end

end
