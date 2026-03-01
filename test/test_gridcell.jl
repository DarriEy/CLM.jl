@testset "GridcellData" begin

    @testset "default construction" begin
        grc = CLM.GridcellData()
        @test length(grc.area) == 0
        @test length(grc.lat) == 0
        @test length(grc.active) == 0
        @test size(grc.landunit_indices) == (0, 0)
    end

    @testset "gridcell_init!" begin
        ng = 5
        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ng)

        # Check sizes
        @test length(grc.area) == ng
        @test length(grc.lat) == ng
        @test length(grc.lon) == ng
        @test length(grc.latdeg) == ng
        @test length(grc.londeg) == ng
        @test length(grc.active) == ng
        @test length(grc.nbedrock) == ng
        @test length(grc.max_dayl) == ng
        @test length(grc.dayl) == ng
        @test length(grc.prev_dayl) == ng
        @test size(grc.landunit_indices) == (CLM.MAX_LUNIT, ng)

        # Check NaN initialization for real fields
        @test all(isnan, grc.area)
        @test all(isnan, grc.lat)
        @test all(isnan, grc.lon)
        @test all(isnan, grc.latdeg)
        @test all(isnan, grc.londeg)
        @test all(isnan, grc.max_dayl)
        @test all(isnan, grc.dayl)
        @test all(isnan, grc.prev_dayl)

        # Check active defaults to true
        @test all(grc.active)

        # Check integer sentinel values
        @test all(==(CLM.ISPVAL), grc.nbedrock)
        @test all(==(CLM.ISPVAL), grc.landunit_indices)
    end

    @testset "gridcell_clean!" begin
        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, 10)
        CLM.gridcell_clean!(grc)

        @test length(grc.area) == 0
        @test length(grc.lat) == 0
        @test length(grc.lon) == 0
        @test length(grc.latdeg) == 0
        @test length(grc.londeg) == 0
        @test length(grc.active) == 0
        @test length(grc.nbedrock) == 0
        @test length(grc.max_dayl) == 0
        @test length(grc.dayl) == 0
        @test length(grc.prev_dayl) == 0
        @test size(grc.landunit_indices) == (0, 0)
    end

    @testset "field mutability" begin
        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, 3)

        # Write to fields and verify
        grc.area[1] = 100.0
        grc.lat[2] = 0.5
        grc.lon[3] = 1.2
        grc.latdeg[1] = 45.0
        grc.londeg[1] = -120.0
        grc.active[2] = false
        grc.nbedrock[3] = 5
        grc.max_dayl[1] = 43200.0
        grc.dayl[1] = 36000.0
        grc.prev_dayl[1] = 35000.0
        grc.landunit_indices[1, 1] = 42

        @test grc.area[1] == 100.0
        @test grc.lat[2] == 0.5
        @test grc.lon[3] == 1.2
        @test grc.latdeg[1] == 45.0
        @test grc.londeg[1] == -120.0
        @test grc.active[2] == false
        @test grc.nbedrock[3] == 5
        @test grc.max_dayl[1] == 43200.0
        @test grc.dayl[1] == 36000.0
        @test grc.prev_dayl[1] == 35000.0
        @test grc.landunit_indices[1, 1] == 42
    end

    @testset "re-init overwrites previous state" begin
        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, 3)
        grc.area[1] = 999.0

        # Re-init with different size
        CLM.gridcell_init!(grc, 7)
        @test length(grc.area) == 7
        @test all(isnan, grc.area)  # old value gone
    end

end
