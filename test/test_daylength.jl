@testset "Daylength" begin
    # --- daylength() pure function tests ---

    @testset "daylength basic" begin
        # Equator at equinox (decl=0) → 12 hours = 43200 s
        dayl = CLM.daylength(0.0, 0.0)
        @test dayl ≈ 43200.0 atol=1.0

        # Equator should give ~12 h regardless of declination
        dayl_summer = CLM.daylength(0.0, 0.3)
        @test dayl_summer ≈ 43200.0 atol=1.0
    end

    @testset "daylength latitude effects" begin
        # Summer in northern hemisphere: longer days at higher latitudes
        decl_summer = 0.4   # ~23° declination (near solstice)
        dayl_45 = CLM.daylength(deg2rad(45.0), decl_summer)
        dayl_60 = CLM.daylength(deg2rad(60.0), decl_summer)
        @test dayl_60 > dayl_45

        # Winter in northern hemisphere: shorter days at higher latitudes
        decl_winter = -0.4
        dayl_45w = CLM.daylength(deg2rad(45.0), decl_winter)
        dayl_60w = CLM.daylength(deg2rad(60.0), decl_winter)
        @test dayl_60w < dayl_45w
    end

    @testset "daylength symmetry" begin
        # Southern hemisphere summer should mirror northern hemisphere summer
        decl = 0.3
        dayl_north = CLM.daylength(deg2rad(45.0), decl)
        dayl_south = CLM.daylength(deg2rad(-45.0), -decl)
        @test dayl_north ≈ dayl_south atol=1e-10
    end

    @testset "daylength polar cases" begin
        # Near pole in summer: should approach 24 hours (86400 s)
        lat_near_pole = deg2rad(89.0)
        decl_summer = 0.4
        dayl = CLM.daylength(lat_near_pole, decl_summer)
        @test dayl > 80000.0  # more than ~22 hours

        # Near pole in winter: should approach 0 hours
        dayl_winter = CLM.daylength(lat_near_pole, -decl_summer)
        @test dayl_winter < 6000.0  # less than ~1.7 hours
    end

    @testset "daylength NaN for invalid inputs" begin
        # lat too far beyond pole
        @test isnan(CLM.daylength(2.0, 0.0))
        @test isnan(CLM.daylength(-2.0, 0.0))

        # decl at or beyond pole
        pole = π / 2.0
        @test isnan(CLM.daylength(0.0, pole))
        @test isnan(CLM.daylength(0.0, -pole))
    end

    @testset "daylength range" begin
        # Daylength should always be between 0 and 86400 for valid inputs
        for lat_deg in -85:5:85
            for decl in -0.4:0.1:0.4
                dayl = CLM.daylength(deg2rad(Float64(lat_deg)), Float64(decl))
                @test 0.0 <= dayl <= 86401.0  # slightly above 86400 possible due to secs_per_radian
            end
        end
    end

    # --- compute_max_daylength! tests ---

    @testset "compute_max_daylength!" begin
        n = 5
        lat = [deg2rad(-60.0), deg2rad(-30.0), 0.0, deg2rad(30.0), deg2rad(60.0)]
        obliquity = deg2rad(23.44)  # Earth's obliquity
        max_dayl = zeros(n)

        CLM.compute_max_daylength!(lat, obliquity, max_dayl, 1:n)

        # At equator, max daylength ≈ 12 hours
        @test max_dayl[3] ≈ 43200.0 atol=1.0

        # Higher latitudes have longer max daylength
        @test max_dayl[4] > max_dayl[3]
        @test max_dayl[5] > max_dayl[4]

        # Symmetric between hemispheres
        @test max_dayl[1] ≈ max_dayl[5] atol=1e-10
        @test max_dayl[2] ≈ max_dayl[4] atol=1e-10
    end

    # --- init_daylength! tests ---

    @testset "init_daylength!" begin
        n = 3
        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, n)

        grc.lat .= [deg2rad(0.0), deg2rad(45.0), deg2rad(-45.0)]

        declin = 0.2
        declinm1 = 0.18
        obliquity = deg2rad(23.44)

        CLM.init_daylength!(grc, declin, declinm1, obliquity, 1:n)

        # All fields should be set (no NaN)
        @test all(!isnan, grc.dayl[1:n])
        @test all(!isnan, grc.prev_dayl[1:n])
        @test all(!isnan, grc.max_dayl[1:n])

        # prev_dayl should differ from dayl (different declinations)
        @test grc.prev_dayl[2] != grc.dayl[2]

        # Verify dayl matches direct daylength call
        @test grc.dayl[1] ≈ CLM.daylength(grc.lat[1], declin)
        @test grc.dayl[2] ≈ CLM.daylength(grc.lat[2], declin)
        @test grc.dayl[3] ≈ CLM.daylength(grc.lat[3], declin)
    end

    # --- update_daylength! tests ---

    @testset "update_daylength!" begin
        n = 2
        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, n)

        grc.lat .= [deg2rad(30.0), deg2rad(50.0)]

        # First: initialize
        declin0 = 0.1
        declinm1 = 0.08
        obliquity = deg2rad(23.44)
        CLM.init_daylength!(grc, declin0, declinm1, obliquity, 1:n)

        saved_dayl = copy(grc.dayl)
        saved_prev = copy(grc.prev_dayl)

        # First step: should NOT update dayl or prev_dayl
        CLM.update_daylength!(grc, 0.15, obliquity, true, 1:n)
        @test grc.dayl ≈ saved_dayl
        @test grc.prev_dayl ≈ saved_prev

        # Second step: should update
        declin_new = 0.2
        CLM.update_daylength!(grc, declin_new, obliquity, false, 1:n)

        # prev_dayl should now equal old dayl
        @test grc.prev_dayl ≈ saved_dayl

        # dayl should be recomputed with new declination
        @test grc.dayl[1] ≈ CLM.daylength(grc.lat[1], declin_new)
        @test grc.dayl[2] ≈ CLM.daylength(grc.lat[2], declin_new)
    end
end
