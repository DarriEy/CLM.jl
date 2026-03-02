@testset "Ozone" begin

    # =====================================================================
    # Constants
    # =====================================================================
    @testset "Ozone Constants" begin
        @test CLM.STRESS_METHOD_LOMBARDOZZI2015 == 1
        @test CLM.STRESS_METHOD_FALK == 2
        @test CLM.O3_KO3 == 1.67
        @test CLM.O3_LAI_THRESH == 0.5
        @test CLM.O3_FLUX_THRESHOLD == 0.8
        @test CLM.O3_NEEDLELEAF_PHOTO_INT == 0.8390
        @test CLM.O3_NEEDLELEAF_PHOTO_SLOPE == 0.0
        @test CLM.O3_BROADLEAF_PHOTO_INT == 0.8752
        @test CLM.O3_BROADLEAF_PHOTO_SLOPE == 0.0
        @test CLM.O3_NONWOODY_PHOTO_INT == 0.8021
        @test CLM.O3_NONWOODY_PHOTO_SLOPE == -0.0009
        @test CLM.O3_NEEDLELEAF_COND_INT == 0.7823
        @test CLM.O3_NEEDLELEAF_COND_SLOPE == 0.0048
        @test CLM.O3_BROADLEAF_COND_INT == 0.9125
        @test CLM.O3_BROADLEAF_COND_SLOPE == 0.0
        @test CLM.O3_NONWOODY_COND_INT == 0.7511
        @test CLM.O3_NONWOODY_COND_SLOPE == 0.0
        @test CLM.O3_NEEDLELEAF_JMAX_INT == 1.0
        @test CLM.O3_NEEDLELEAF_JMAX_SLOPE == 0.0
        @test CLM.O3_BROADLEAF_JMAX_INT == 1.0
        @test CLM.O3_BROADLEAF_JMAX_SLOPE == -0.0037
        @test CLM.O3_NONWOODY_JMAX_INT == 1.0
        @test CLM.O3_NONWOODY_JMAX_SLOPE == 0.0
    end

    # =====================================================================
    # OzoneData init/clean
    # =====================================================================
    @testset "OzoneData init/clean" begin
        oz = CLM.OzoneData()
        CLM.ozone_init!(oz, 5; stress_method="stress_lombardozzi2015")

        @test oz.stress_method == CLM.STRESS_METHOD_LOMBARDOZZI2015
        @test length(oz.o3uptakesha_patch) == 5
        @test length(oz.o3uptakesun_patch) == 5
        @test length(oz.tlai_old_patch) == 5
        @test length(oz.o3coefvsha_patch) == 5
        @test length(oz.o3coefvsun_patch) == 5
        @test length(oz.o3coefgsha_patch) == 5
        @test length(oz.o3coefgsun_patch) == 5
        @test length(oz.o3coefjmaxsha_patch) == 5
        @test length(oz.o3coefjmaxsun_patch) == 5

        # Cold-start values
        @test all(oz.o3uptakesha_patch .== 0.0)
        @test all(oz.o3uptakesun_patch .== 0.0)
        @test all(oz.tlai_old_patch .== 0.0)
        @test all(oz.o3coefvsha_patch .== 1.0)
        @test all(oz.o3coefvsun_patch .== 1.0)
        @test all(oz.o3coefgsha_patch .== 1.0)
        @test all(oz.o3coefgsun_patch .== 1.0)
        @test all(oz.o3coefjmaxsha_patch .== 1.0)
        @test all(oz.o3coefjmaxsun_patch .== 1.0)

        CLM.ozone_clean!(oz)
        @test isempty(oz.o3uptakesha_patch)
        @test isempty(oz.o3uptakesun_patch)
        @test isempty(oz.tlai_old_patch)
        @test isempty(oz.o3coefvsha_patch)
    end

    @testset "OzoneData init Falk requires LUNA" begin
        oz = CLM.OzoneData()
        @test_throws ErrorException CLM.ozone_init!(oz, 5; stress_method="stress_falk", use_luna=false)
        # Should succeed with use_luna=true
        CLM.ozone_init!(oz, 5; stress_method="stress_falk", use_luna=true)
        @test oz.stress_method == CLM.STRESS_METHOD_FALK
    end

    @testset "OzoneData init unknown method" begin
        oz = CLM.OzoneData()
        @test_throws ErrorException CLM.ozone_init!(oz, 5; stress_method="unknown")
    end

    # =====================================================================
    # Ozone uptake — single point
    # =====================================================================
    @testset "calc_ozone_uptake_one_point basic" begin
        # PFT parameters (78+1 entries, using simple defaults)
        np = CLM.MXPFT + 1
        evergreen_pft = zeros(np)
        leaf_long_pft = fill(1.0, np)
        dtime = 1800.0  # 30 min timestep

        # Test with zero ozone forcing -> no uptake
        result = CLM.calc_ozone_uptake_one_point(
            forc_ozone = 0.0,
            forc_pbot  = 101325.0,
            forc_th    = 300.0,
            rs         = 200.0,
            rb         = 50.0,
            ram        = 100.0,
            tlai_val   = 2.0,
            tlai_old   = 2.0,
            pft_type   = 5,
            o3uptake   = 0.0,
            evergreen_pft = evergreen_pft,
            leaf_long_pft = leaf_long_pft,
            dtime      = dtime)
        @test result == 0.0

        # Test with LAI below threshold -> uptake reset to 0
        result = CLM.calc_ozone_uptake_one_point(
            forc_ozone = 50.0e-9,
            forc_pbot  = 101325.0,
            forc_th    = 300.0,
            rs         = 200.0,
            rb         = 50.0,
            ram        = 100.0,
            tlai_val   = 0.3,    # below O3_LAI_THRESH = 0.5
            tlai_old   = 0.3,
            pft_type   = 5,
            o3uptake   = 10.0,
            evergreen_pft = evergreen_pft,
            leaf_long_pft = leaf_long_pft,
            dtime      = dtime)
        @test result == 0.0
    end

    @testset "calc_ozone_uptake_one_point accumulation" begin
        np = CLM.MXPFT + 1
        evergreen_pft = zeros(np)
        leaf_long_pft = fill(1.0, np)
        dtime = 1800.0

        # Moderate ozone exposure with LAI above threshold -> positive uptake
        result = CLM.calc_ozone_uptake_one_point(
            forc_ozone = 50.0e-9,
            forc_pbot  = 101325.0,
            forc_th    = 300.0,
            rs         = 200.0,
            rb         = 50.0,
            ram        = 100.0,
            tlai_val   = 3.0,
            tlai_old   = 3.0,
            pft_type   = 5,
            o3uptake   = 0.0,
            evergreen_pft = evergreen_pft,
            leaf_long_pft = leaf_long_pft,
            dtime      = dtime)
        @test result >= 0.0
        @test isfinite(result)
    end

    @testset "calc_ozone_uptake_one_point healing" begin
        np = CLM.MXPFT + 1
        evergreen_pft = zeros(np)
        leaf_long_pft = fill(1.0, np)
        dtime = 1800.0

        # Growing leaves (tlai > tlai_old) should reduce damage via healing
        result_growing = CLM.calc_ozone_uptake_one_point(
            forc_ozone = 50.0e-9,
            forc_pbot  = 101325.0,
            forc_th    = 300.0,
            rs         = 200.0,
            rb         = 50.0,
            ram        = 100.0,
            tlai_val   = 4.0,
            tlai_old   = 2.0,
            pft_type   = 5,
            o3uptake   = 1.0,
            evergreen_pft = evergreen_pft,
            leaf_long_pft = leaf_long_pft,
            dtime      = dtime)

        result_stable = CLM.calc_ozone_uptake_one_point(
            forc_ozone = 50.0e-9,
            forc_pbot  = 101325.0,
            forc_th    = 300.0,
            rs         = 200.0,
            rb         = 50.0,
            ram        = 100.0,
            tlai_val   = 4.0,
            tlai_old   = 4.0,
            pft_type   = 5,
            o3uptake   = 1.0,
            evergreen_pft = evergreen_pft,
            leaf_long_pft = leaf_long_pft,
            dtime      = dtime)

        # Growing leaves should result in less uptake than stable leaves
        @test result_growing <= result_stable
    end

    @testset "calc_ozone_uptake_one_point evergreen decay" begin
        np = CLM.MXPFT + 1
        evergreen_pft = zeros(np)
        evergreen_pft[2] = 1.0  # PFT 2 is needleleaf_evergreen_boreal_tree
        leaf_long_pft = fill(1.0, np)
        leaf_long_pft[2] = 3.0  # 3-year leaf longevity
        dtime = 1800.0

        # Evergreen with existing uptake should have decay
        result_evergreen = CLM.calc_ozone_uptake_one_point(
            forc_ozone = 0.0,    # no new ozone
            forc_pbot  = 101325.0,
            forc_th    = 300.0,
            rs         = 200.0,
            rb         = 50.0,
            ram        = 100.0,
            tlai_val   = 3.0,
            tlai_old   = 3.0,
            pft_type   = 2,
            o3uptake   = 5.0,
            evergreen_pft = evergreen_pft,
            leaf_long_pft = leaf_long_pft,
            dtime      = dtime)

        # With no new ozone and evergreen decay, uptake should decrease
        @test result_evergreen < 5.0
        @test result_evergreen >= 0.0
    end

    # =====================================================================
    # Ozone stress — Lombardozzi 2015 one point
    # =====================================================================
    @testset "calc_ozone_stress_lombardozzi2015_one_point" begin
        np = CLM.MXPFT + 1
        woody_pft = zeros(np)
        # Set woody for broadleaf types (pft_type > 3)
        for i in 4:9
            woody_pft[i] = 1.0
        end

        # Zero uptake -> coefficients = 1
        o3coefv, o3coefg = CLM.calc_ozone_stress_lombardozzi2015_one_point(
            pft_type = 5, o3uptake = 0.0, woody_pft = woody_pft)
        @test o3coefv == 1.0
        @test o3coefg == 1.0

        # Needleleaf (pft_type <= 3)
        o3coefv, o3coefg = CLM.calc_ozone_stress_lombardozzi2015_one_point(
            pft_type = 2, o3uptake = 5.0, woody_pft = woody_pft)
        @test 0.0 <= o3coefv <= 1.0
        @test 0.0 <= o3coefg <= 1.0
        # needleleaf photo slope is 0, so o3coefv = needleleafPhotoInt
        @test o3coefv ≈ CLM.O3_NEEDLELEAF_PHOTO_INT
        # needleleaf cond slope is positive, so o3coefg > condInt
        @test o3coefg ≈ min(1.0, CLM.O3_NEEDLELEAF_COND_INT + CLM.O3_NEEDLELEAF_COND_SLOPE * 5.0)

        # Broadleaf (pft_type > 3, woody = 1)
        o3coefv, o3coefg = CLM.calc_ozone_stress_lombardozzi2015_one_point(
            pft_type = 5, o3uptake = 5.0, woody_pft = woody_pft)
        @test o3coefv ≈ CLM.O3_BROADLEAF_PHOTO_INT  # slope is 0
        @test o3coefg ≈ CLM.O3_BROADLEAF_COND_INT   # slope is 0

        # Nonwoody (pft_type > 3, woody = 0)
        o3coefv, o3coefg = CLM.calc_ozone_stress_lombardozzi2015_one_point(
            pft_type = 13, o3uptake = 10.0, woody_pft = woody_pft)
        expected_v = max(0.0, min(1.0, CLM.O3_NONWOODY_PHOTO_INT + CLM.O3_NONWOODY_PHOTO_SLOPE * 10.0))
        expected_g = max(0.0, min(1.0, CLM.O3_NONWOODY_COND_INT + CLM.O3_NONWOODY_COND_SLOPE * 10.0))
        @test o3coefv ≈ expected_v
        @test o3coefg ≈ expected_g
    end

    # =====================================================================
    # Ozone stress — Falk one point
    # =====================================================================
    @testset "calc_ozone_stress_falk_one_point" begin
        np = CLM.MXPFT + 1
        woody_pft = zeros(np)
        for i in 4:9
            woody_pft[i] = 1.0
        end

        # Zero uptake -> coefficient = 1
        result = CLM.calc_ozone_stress_falk_one_point(
            pft_type = 5, o3uptake = 0.0, o3coefjmax = 0.8, woody_pft = woody_pft)
        @test result == 1.0

        # Broadleaf (pft > 3, woody = 1) with uptake
        result = CLM.calc_ozone_stress_falk_one_point(
            pft_type = 5, o3uptake = 10.0, o3coefjmax = 1.0, woody_pft = woody_pft)
        expected = max(0.0, min(1.0, CLM.O3_BROADLEAF_JMAX_INT + CLM.O3_BROADLEAF_JMAX_SLOPE * 10.0))
        @test result ≈ expected

        # Needleleaf (pft <= 3) — slopes are 0
        result = CLM.calc_ozone_stress_falk_one_point(
            pft_type = 2, o3uptake = 10.0, o3coefjmax = 1.0, woody_pft = woody_pft)
        @test result ≈ 1.0  # needleleaf slope is 0, int is 1

        # Nonwoody (pft > 3, woody = 0) — slopes are 0
        result = CLM.calc_ozone_stress_falk_one_point(
            pft_type = 13, o3uptake = 10.0, o3coefjmax = 1.0, woody_pft = woody_pft)
        @test result ≈ 1.0  # nonwoody slope is 0, int is 1
    end

    # =====================================================================
    # Full calc_ozone_uptake! integration test
    # =====================================================================
    @testset "calc_ozone_uptake! integration" begin
        npatches = 3
        ncols = 2
        ngridcells = 1

        # Set up ozone data
        oz = CLM.OzoneData()
        CLM.ozone_init!(oz, npatches; stress_method="stress_lombardozzi2015")

        # Set up patch data
        pch = CLM.PatchData()
        pch.column   = [1, 1, 2]
        pch.gridcell = [1, 1, 1]
        pch.itype    = [2, 5, 13]   # needleleaf, broadleaf, nonwoody

        # PFT parameters
        np = CLM.MXPFT + 1
        evergreen_pft = zeros(np)
        evergreen_pft[2] = 1.0
        leaf_long_pft = fill(1.0, np)
        leaf_long_pft[2] = 3.0

        mask = BitVector([true, true, false])  # only first two exposed
        bounds = 1:3
        forc_pbot = [101325.0, 101325.0]
        forc_th   = [300.0, 300.0]
        rssun     = [200.0, 150.0, 300.0]
        rssha     = [250.0, 180.0, 350.0]
        rb_vec    = [50.0, 60.0, 70.0]
        ram_vec   = [100.0, 110.0, 120.0]
        tlai_vec  = [3.0, 4.0, 1.0]
        forc_o3   = [50.0e-9]

        CLM.calc_ozone_uptake!(oz, pch, mask, bounds,
            forc_pbot, forc_th, rssun, rssha, rb_vec, ram_vec,
            tlai_vec, forc_o3, evergreen_pft, leaf_long_pft, 1800.0)

        # First two patches should have uptake computed
        @test oz.o3uptakesha_patch[1] >= 0.0
        @test oz.o3uptakesun_patch[1] >= 0.0
        @test oz.o3uptakesha_patch[2] >= 0.0
        @test oz.o3uptakesun_patch[2] >= 0.0

        # Third patch (not exposed) should be unchanged (zero from init)
        @test oz.o3uptakesha_patch[3] == 0.0
        @test oz.o3uptakesun_patch[3] == 0.0

        # tlai_old should be updated for exposed patches
        @test oz.tlai_old_patch[1] == 3.0
        @test oz.tlai_old_patch[2] == 4.0
        @test oz.tlai_old_patch[3] == 0.0
    end

    # =====================================================================
    # Full calc_ozone_stress! integration — Lombardozzi2015
    # =====================================================================
    @testset "calc_ozone_stress! Lombardozzi2015 integration" begin
        npatches = 3
        oz = CLM.OzoneData()
        CLM.ozone_init!(oz, npatches; stress_method="stress_lombardozzi2015")

        # Set some uptake values
        oz.o3uptakesha_patch .= [5.0, 10.0, 0.0]
        oz.o3uptakesun_patch .= [3.0, 8.0, 0.0]

        pch = CLM.PatchData()
        pch.itype = [2, 13, 5]

        np = CLM.MXPFT + 1
        woody_pft = zeros(np)
        for i in 4:9
            woody_pft[i] = 1.0
        end

        mask_exposed   = BitVector([true, true, false])
        mask_noexposed = BitVector([false, false, true])
        bounds = 1:3

        CLM.calc_ozone_stress!(oz, mask_exposed, mask_noexposed, bounds, pch, woody_pft)

        # Exposed patches should have updated stress coefficients
        @test 0.0 <= oz.o3coefvsha_patch[1] <= 1.0
        @test 0.0 <= oz.o3coefvsun_patch[1] <= 1.0
        @test 0.0 <= oz.o3coefgsha_patch[1] <= 1.0
        @test 0.0 <= oz.o3coefgsun_patch[1] <= 1.0

        # Non-exposed patch should be reset to 1
        @test oz.o3coefvsha_patch[3] == 1.0
        @test oz.o3coefvsun_patch[3] == 1.0
        @test oz.o3coefgsha_patch[3] == 1.0
        @test oz.o3coefgsun_patch[3] == 1.0
    end

    # =====================================================================
    # Full calc_ozone_stress! integration — Falk
    # =====================================================================
    @testset "calc_ozone_stress! Falk integration" begin
        npatches = 3
        oz = CLM.OzoneData()
        CLM.ozone_init!(oz, npatches; stress_method="stress_falk", use_luna=true)

        oz.o3uptakesha_patch .= [5.0, 10.0, 0.0]
        oz.o3uptakesun_patch .= [3.0, 8.0, 0.0]

        pch = CLM.PatchData()
        pch.itype = [2, 5, 13]

        np = CLM.MXPFT + 1
        woody_pft = zeros(np)
        for i in 4:9
            woody_pft[i] = 1.0
        end

        mask_exposed   = BitVector([true, true, false])
        mask_noexposed = BitVector([false, false, true])
        bounds = 1:3

        # When is_time_to_run_luna=false, Falk stress should NOT run
        CLM.calc_ozone_stress!(oz, mask_exposed, mask_noexposed, bounds, pch, woody_pft;
                               is_time_to_run_luna=false)
        # All should still be 1.0 (initial cold-start value)
        @test oz.o3coefjmaxsha_patch[1] == 1.0
        @test oz.o3coefjmaxsun_patch[1] == 1.0

        # When is_time_to_run_luna=true, Falk stress should run
        CLM.calc_ozone_stress!(oz, mask_exposed, mask_noexposed, bounds, pch, woody_pft;
                               is_time_to_run_luna=true)

        # Exposed patches should have updated jmax coefficients
        @test 0.0 <= oz.o3coefjmaxsha_patch[1] <= 1.0
        @test 0.0 <= oz.o3coefjmaxsun_patch[1] <= 1.0
        @test 0.0 <= oz.o3coefjmaxsha_patch[2] <= 1.0

        # Broadleaf (pft 5) with uptake 10 should show reduction
        expected_sha = max(0.0, min(1.0, CLM.O3_BROADLEAF_JMAX_INT + CLM.O3_BROADLEAF_JMAX_SLOPE * 10.0))
        @test oz.o3coefjmaxsha_patch[2] ≈ expected_sha

        # Non-exposed patch should be reset to 1
        @test oz.o3coefjmaxsha_patch[3] == 1.0
        @test oz.o3coefjmaxsun_patch[3] == 1.0
    end

    # =====================================================================
    # o3 concentration conversion sanity check
    # =====================================================================
    @testset "o3 concentration conversion" begin
        # Verify the mol/mol -> nmol/m^3 conversion manually
        forc_ozone = 50.0e-9   # 50 ppb
        forc_pbot  = 101325.0  # Pa
        forc_th    = 300.0     # K

        o3concnmolm3 = forc_ozone * 1.0e9 * (forc_pbot / (forc_th * CLM.RGAS * 0.001))
        @test o3concnmolm3 > 0.0
        @test isfinite(o3concnmolm3)

        # Expected: ~50 * 101325 / (300 * 8.31446 * 0.001) = ~50 * 40605 ≈ 2.03e6
        @test o3concnmolm3 ≈ 50.0 * 101325.0 / (300.0 * CLM.RGAS * 0.001) rtol=1e-10
    end

end
