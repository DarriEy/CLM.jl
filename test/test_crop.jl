@testset "CropData" begin

    @testset "cphase constants" begin
        @test CLM.cphase_not_planted == 0.0
        @test CLM.cphase_planted     == 1.0
        @test CLM.cphase_leafemerge  == 2.0
        @test CLM.cphase_grainfill   == 3.0
        @test CLM.cphase_harvest     == 4.0
    end

    @testset "baset mapping constants" begin
        @test CLM.BASET_MAP_CONSTANT == "constant"
        @test CLM.BASET_MAP_LATVARY  == "varytropicsbylat"
    end

    @testset "default construction" begin
        cr = CLM.CropData()
        @test length(cr.nyrs_crop_active_patch) == 0
        @test length(cr.croplive_patch) == 0
        @test length(cr.fertnitro_patch) == 0
        @test length(cr.hui_patch) == 0
        @test size(cr.sdates_thisyr_patch) == (0, 0)
        @test size(cr.hdates_thisyr_patch) == (0, 0)
        @test cr.baset_mapping == CLM.BASET_MAP_CONSTANT
        @test cr.baset_latvary_intercept == 12.0
        @test cr.baset_latvary_slope == 0.4
    end

    @testset "crop_init! basic allocation" begin
        np = 10
        mxs = CLM.MXSOWINGS
        mxh = CLM.MXHARVESTS

        cr = CLM.CropData()
        CLM.crop_init!(cr, np)

        # Integer 1D
        @test length(cr.nyrs_crop_active_patch) == np
        @test all(x -> x == 0, cr.nyrs_crop_active_patch)
        @test length(cr.harvdate_patch) == np
        @test all(x -> x == typemax(Int), cr.harvdate_patch)
        @test length(cr.sowing_reason_patch) == np
        @test all(x -> x == -1, cr.sowing_reason_patch)
        @test length(cr.sowing_count) == np
        @test all(x -> x == 0, cr.sowing_count)
        @test length(cr.harvest_count) == np
        @test all(x -> x == 0, cr.harvest_count)

        # Boolean 1D
        @test length(cr.croplive_patch) == np
        @test !any(cr.croplive_patch)
        @test length(cr.sown_in_this_window) == np
        @test !any(cr.sown_in_this_window)

        # Real 1D — SPVAL initialized
        @test length(cr.fertnitro_patch) == np
        @test all(x -> x == CLM.SPVAL, cr.fertnitro_patch)
        @test length(cr.hui_patch) == np
        @test all(x -> x == CLM.SPVAL, cr.hui_patch)
        @test length(cr.gddaccum_patch) == np
        @test all(x -> x == CLM.SPVAL, cr.gddaccum_patch)
        @test length(cr.gddtsoi_patch) == np
        @test all(x -> x == CLM.SPVAL, cr.gddtsoi_patch)
        @test length(cr.latbaset_patch) == np
        @test all(x -> x == CLM.SPVAL, cr.latbaset_patch)

        # Real 1D — zero initialized
        @test length(cr.vf_patch) == np
        @test all(x -> x == 0.0, cr.vf_patch)

        # Real 1D — cphase_not_planted
        @test length(cr.cphase_patch) == np
        @test all(x -> x == CLM.cphase_not_planted, cr.cphase_patch)

        # GDD20 fields — SPVAL
        @test all(x -> x == CLM.SPVAL, cr.gdd20_baseline_patch)
        @test all(x -> x == CLM.SPVAL, cr.gdd20_season_start_patch)
        @test all(x -> x == CLM.SPVAL, cr.gdd20_season_end_patch)

        # Integer 2D (patch × mxsowings)
        @test size(cr.rx_swindow_starts_thisyr_patch) == (np, mxs)
        @test all(x -> x == -1, cr.rx_swindow_starts_thisyr_patch)
        @test size(cr.rx_swindow_ends_thisyr_patch) == (np, mxs)
        @test all(x -> x == -1, cr.rx_swindow_ends_thisyr_patch)

        # Real 2D (patch × mxsowings) — SPVAL
        @test size(cr.rx_cultivar_gdds_thisyr_patch) == (np, mxs)
        @test all(x -> x == CLM.SPVAL, cr.rx_cultivar_gdds_thisyr_patch)
        @test size(cr.sdates_thisyr_patch) == (np, mxs)
        @test all(x -> x == CLM.SPVAL, cr.sdates_thisyr_patch)
        @test size(cr.swindow_starts_thisyr_patch) == (np, mxs)
        @test all(x -> x == CLM.SPVAL, cr.swindow_starts_thisyr_patch)
        @test size(cr.swindow_ends_thisyr_patch) == (np, mxs)
        @test all(x -> x == CLM.SPVAL, cr.swindow_ends_thisyr_patch)

        # Real 2D (patch × mxharvests) — SPVAL
        @test size(cr.sdates_perharv_patch) == (np, mxh)
        @test all(x -> x == CLM.SPVAL, cr.sdates_perharv_patch)
        @test size(cr.syears_perharv_patch) == (np, mxh)
        @test all(x -> x == CLM.SPVAL, cr.syears_perharv_patch)
        @test size(cr.hdates_thisyr_patch) == (np, mxh)
        @test all(x -> x == CLM.SPVAL, cr.hdates_thisyr_patch)
        @test size(cr.gddaccum_thisyr_patch) == (np, mxh)
        @test all(x -> x == CLM.SPVAL, cr.gddaccum_thisyr_patch)
        @test size(cr.hui_thisyr_patch) == (np, mxh)
        @test all(x -> x == CLM.SPVAL, cr.hui_thisyr_patch)
        @test size(cr.sowing_reason_perharv_patch) == (np, mxh)
        @test all(x -> x == CLM.SPVAL, cr.sowing_reason_perharv_patch)
        @test size(cr.harvest_reason_thisyr_patch) == (np, mxh)
        @test all(x -> x == CLM.SPVAL, cr.harvest_reason_thisyr_patch)
    end

    @testset "crop_read_nml!" begin
        cr = CLM.CropData()

        # Default (constant)
        CLM.crop_read_nml!(cr)
        @test cr.baset_mapping == CLM.BASET_MAP_CONSTANT
        @test cr.baset_latvary_intercept == 12.0
        @test cr.baset_latvary_slope == 0.4

        # Latvary
        CLM.crop_read_nml!(cr; baset_mapping=CLM.BASET_MAP_LATVARY,
                           baset_latvary_intercept=10.0,
                           baset_latvary_slope=0.3)
        @test cr.baset_mapping == CLM.BASET_MAP_LATVARY
        @test cr.baset_latvary_intercept == 10.0
        @test cr.baset_latvary_slope == 0.3

        # Invalid mapping
        @test_throws ErrorException CLM.crop_read_nml!(cr; baset_mapping="invalid")
    end

    @testset "crop_init_cold! constant mapping" begin
        np = 5
        cr = CLM.CropData()
        CLM.crop_init!(cr, np)
        cr.baset_mapping = CLM.BASET_MAP_CONSTANT

        CLM.crop_init_cold!(cr, 1:np)

        # With constant mapping, latbaset should be NaN
        @test all(isnan, cr.latbaset_patch)
        # nyrs_crop_active should be 0
        @test all(x -> x == 0, cr.nyrs_crop_active_patch)
    end

    @testset "crop_init_cold! latvary mapping with crop patches" begin
        np = 3
        cr = CLM.CropData()
        CLM.crop_init!(cr, np)
        cr.baset_mapping = CLM.BASET_MAP_LATVARY
        cr.baset_latvary_intercept = 12.0
        cr.baset_latvary_slope = 0.4

        # Simulate: patches 1,3 are crop (landunit type == istcrop); patch 2 is not
        patch_landunit = [1, 2, 1]
        lun_itype = [15, 1]  # landunit 1 is crop (istcrop=15), landunit 2 is not
        istcrop_val = 15
        patch_gridcell = [1, 1, 2]
        patch_itype = [0, 1, 0]  # 0-based Fortran PFT indices
        fert_cft = [100.0 200.0; 50.0 150.0]  # (gridcell × ivt)
        pftcon_baset = [5.0, 8.0]
        grc_latdeg = [30.0, 60.0]

        CLM.crop_init_cold!(cr, 1:np;
                            patch_landunit=patch_landunit,
                            lun_itype=lun_itype,
                            istcrop=istcrop_val,
                            patch_gridcell=patch_gridcell,
                            patch_itype=patch_itype,
                            fert_cft=fert_cft,
                            pftcon_baset=pftcon_baset,
                            grc_latdeg=grc_latdeg)

        # Patch 1: crop, gridcell=1, ivt=1 => fertnitro = fert_cft[1,1] = 100.0
        @test cr.fertnitro_patch[1] == 100.0
        # Patch 1: latbaset = latbaset(5.0, 30.0, 12.0, 0.4)
        expected_lb1 = CLM.latbaset(5.0, 30.0, 12.0, 0.4)
        @test cr.latbaset_patch[1] == expected_lb1

        # Patch 2: not crop — fertnitro unchanged (SPVAL from init)
        @test cr.fertnitro_patch[2] == CLM.SPVAL

        # Patch 3: crop, gridcell=2, ivt=1 => fertnitro = fert_cft[2,1] = 50.0
        @test cr.fertnitro_patch[3] == 50.0
        expected_lb3 = CLM.latbaset(5.0, 60.0, 12.0, 0.4)
        @test cr.latbaset_patch[3] == expected_lb3
    end

    @testset "latbaset function" begin
        # At equator (lat=0): baset + intercept - min(intercept, 0) = baset + intercept
        @test CLM.latbaset(5.0, 0.0, 12.0, 0.4) ≈ 5.0 + 12.0

        # At maxlat = intercept/slope = 12.0/0.4 = 30.0:
        # baset + intercept - min(intercept, slope*30) = baset + 12 - 12 = baset
        @test CLM.latbaset(5.0, 30.0, 12.0, 0.4) ≈ 5.0

        # Beyond maxlat (lat=50): slope*50 = 20 > intercept=12
        # baset + 12 - min(12, 20) = baset + 12 - 12 = baset
        @test CLM.latbaset(5.0, 50.0, 12.0, 0.4) ≈ 5.0

        # Negative latitude should give same result as positive (abs(lat))
        @test CLM.latbaset(5.0, -30.0, 12.0, 0.4) ≈ CLM.latbaset(5.0, 30.0, 12.0, 0.4)

        # Intermediate latitude (lat=15): slope*15 = 6 < intercept=12
        # baset + 12 - 6 = baset + 6
        @test CLM.latbaset(5.0, 15.0, 12.0, 0.4) ≈ 11.0
    end

    @testset "crop_increment_year!" begin
        np = 4
        cr = CLM.CropData()
        CLM.crop_init!(cr, np)

        # Patches 1,3 are crop; 2,4 are not
        mask = BitVector([true, false, true, false])

        # Not at year boundary — no increment
        CLM.crop_increment_year!(cr, mask, 1:np; kmo=6, kda=15, mcsec=0, is_first_step=false)
        @test all(x -> x == 0, cr.nyrs_crop_active_patch)

        # At year boundary, not first step — increment crop patches
        CLM.crop_increment_year!(cr, mask, 1:np; kmo=1, kda=1, mcsec=0, is_first_step=false)
        @test cr.nyrs_crop_active_patch[1] == 1
        @test cr.nyrs_crop_active_patch[2] == 0  # not crop
        @test cr.nyrs_crop_active_patch[3] == 1
        @test cr.nyrs_crop_active_patch[4] == 0

        # At year boundary but first step — no increment
        CLM.crop_increment_year!(cr, mask, 1:np; kmo=1, kda=1, mcsec=0, is_first_step=true)
        @test cr.nyrs_crop_active_patch[1] == 1  # unchanged
        @test cr.nyrs_crop_active_patch[3] == 1

        # Second year boundary
        CLM.crop_increment_year!(cr, mask, 1:np; kmo=1, kda=1, mcsec=0, is_first_step=false)
        @test cr.nyrs_crop_active_patch[1] == 2
        @test cr.nyrs_crop_active_patch[3] == 2
    end

    @testset "stub functions run without error" begin
        cr = CLM.CropData()
        CLM.crop_init!(cr, 5)
        @test CLM.crop_init_history!(cr, 1:5) === nothing
        @test CLM.crop_init_acc_buffer!(cr) === nothing
        @test CLM.crop_init_acc_vars!(cr, 1:5) === nothing
        @test CLM.crop_restart!(cr, 1:5) === nothing
        @test CLM.crop_check_dates!() === nothing
        @test CLM.crop_update_acc_vars!(cr, 1:5, zeros(5), zeros(5, 2)) === nothing
    end

    @testset "field mutability" begin
        cr = CLM.CropData()
        CLM.crop_init!(cr, 3)

        cr.croplive_patch[1] = true
        @test cr.croplive_patch[1] == true

        cr.hui_patch[2] = 42.0
        @test cr.hui_patch[2] == 42.0

        cr.cphase_patch[3] = CLM.cphase_grainfill
        @test cr.cphase_patch[3] == CLM.cphase_grainfill

        cr.sdates_thisyr_patch[1, 1] = 120.0
        @test cr.sdates_thisyr_patch[1, 1] == 120.0
    end

    @testset "re-init overwrites previous state" begin
        cr = CLM.CropData()
        CLM.crop_init!(cr, 3)
        cr.hui_patch[1] = 999.0

        CLM.crop_init!(cr, 7)
        @test length(cr.hui_patch) == 7
        @test all(x -> x == CLM.SPVAL, cr.hui_patch)
        @test length(cr.croplive_patch) == 7
        @test !any(cr.croplive_patch)
        @test size(cr.sdates_thisyr_patch) == (7, CLM.MXSOWINGS)
    end

end
