@testset "Atm2LndData" begin

    @testset "Atm2LndParamsData default construction" begin
        p = CLM.Atm2LndParamsData()
        @test p.repartition_rain_snow == false
        @test p.glcmec_downscale_longwave == false
        @test isnan(p.lapse_rate)
        @test isnan(p.lapse_rate_longwave)
        @test isnan(p.longwave_downscaling_limit)
        @test isnan(p.precip_repartition_glc_all_snow_t)
        @test isnan(p.precip_repartition_glc_frac_rain_slope)
        @test isnan(p.precip_repartition_nonglc_all_snow_t)
        @test isnan(p.precip_repartition_nonglc_frac_rain_slope)
    end

    @testset "Atm2LndData default construction" begin
        a = CLM.Atm2LndData()
        @test length(a.forc_u_grc) == 0
        @test length(a.forc_t_downscaled_col) == 0
        @test length(a.fsd24_patch) == 0
        @test size(a.forc_solad_not_downscaled_grc) == (0, 0)
        @test size(a.forc_aer_grc) == (0, 0)
    end

    @testset "compute_ramp_params" begin
        # all_snow_t_c = -5, all_rain_t_c = 5 => slope = 1/10 = 0.1
        (snow_k, slope) = CLM.compute_ramp_params(-5.0, 5.0)
        @test snow_k ≈ -5.0 + CLM.TFRZ
        @test slope ≈ 0.1

        # all_snow_t_c = 0, all_rain_t_c = 2 => slope = 0.5
        (snow_k2, slope2) = CLM.compute_ramp_params(0.0, 2.0)
        @test snow_k2 ≈ CLM.TFRZ
        @test slope2 ≈ 0.5
    end

    @testset "atm2lnd_params_init! basic (no repartition, no downscale)" begin
        p = CLM.Atm2LndParamsData()
        CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.006)

        @test p.repartition_rain_snow == false
        @test p.glcmec_downscale_longwave == false
        @test p.lapse_rate ≈ 0.006
        @test isnan(p.lapse_rate_longwave)
        @test isnan(p.longwave_downscaling_limit)
        @test isnan(p.precip_repartition_glc_all_snow_t)
        @test isnan(p.precip_repartition_glc_frac_rain_slope)
        @test isnan(p.precip_repartition_nonglc_all_snow_t)
        @test isnan(p.precip_repartition_nonglc_frac_rain_slope)
    end

    @testset "atm2lnd_params_init! with longwave downscaling" begin
        p = CLM.Atm2LndParamsData()
        CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = true,
            lapse_rate = 0.006,
            lapse_rate_longwave = 0.032,
            longwave_downscaling_limit = 0.5)

        @test p.glcmec_downscale_longwave == true
        @test p.lapse_rate_longwave ≈ 0.032
        @test p.longwave_downscaling_limit ≈ 0.5
    end

    @testset "atm2lnd_params_init! longwave downscale error checks" begin
        p = CLM.Atm2LndParamsData()
        # Missing lapse_rate_longwave
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = true,
            lapse_rate = 0.006)

        # Missing longwave_downscaling_limit
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = true,
            lapse_rate = 0.006,
            lapse_rate_longwave = 0.032)

        # longwave_downscaling_limit out of range
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = true,
            lapse_rate = 0.006,
            lapse_rate_longwave = 0.032,
            longwave_downscaling_limit = 1.5)
    end

    @testset "atm2lnd_params_init! with repartitioning" begin
        p = CLM.Atm2LndParamsData()
        CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = true,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.006,
            precip_repartition_glc_all_snow_t = -2.0,
            precip_repartition_glc_all_rain_t = 0.0,
            precip_repartition_nonglc_all_snow_t = -5.0,
            precip_repartition_nonglc_all_rain_t = 2.5)

        @test p.repartition_rain_snow == true
        # glc: all_snow_t_k = -2.0 + TFRZ, slope = 1/(0-(-2)) = 0.5
        @test p.precip_repartition_glc_all_snow_t ≈ -2.0 + CLM.TFRZ
        @test p.precip_repartition_glc_frac_rain_slope ≈ 0.5
        # nonglc: all_snow_t_k = -5.0 + TFRZ, slope = 1/(2.5-(-5)) = 1/7.5
        @test p.precip_repartition_nonglc_all_snow_t ≈ -5.0 + CLM.TFRZ
        @test p.precip_repartition_nonglc_frac_rain_slope ≈ 1.0 / 7.5
    end

    @testset "atm2lnd_params_init! repartitioning error checks" begin
        p = CLM.Atm2LndParamsData()
        # Missing required param
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = true,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.006)

        # all_rain <= all_snow (glc)
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = true,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.006,
            precip_repartition_glc_all_snow_t = 0.0,
            precip_repartition_glc_all_rain_t = -1.0,
            precip_repartition_nonglc_all_snow_t = -5.0,
            precip_repartition_nonglc_all_rain_t = 2.5)

        # all_rain <= all_snow (nonglc)
        @test_throws ErrorException CLM.atm2lnd_params_init!(p;
            repartition_rain_snow = true,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.006,
            precip_repartition_glc_all_snow_t = -2.0,
            precip_repartition_glc_all_rain_t = 0.0,
            precip_repartition_nonglc_all_snow_t = 2.5,
            precip_repartition_nonglc_all_rain_t = 2.5)
    end

    @testset "atm2lnd_init! allocation" begin
        ng, nc, np = 4, 8, 12
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, ng, nc, np)

        # Gridcell-level vectors
        @test length(a.forc_u_grc) == ng
        @test all(x -> x == 0.0, a.forc_u_grc)
        @test length(a.forc_v_grc) == ng
        @test length(a.forc_wind_grc) == ng
        @test length(a.forc_hgt_grc) == ng
        @test length(a.forc_topo_grc) == ng
        @test length(a.forc_hgt_u_grc) == ng
        @test length(a.forc_hgt_t_grc) == ng
        @test length(a.forc_hgt_q_grc) == ng
        @test length(a.forc_vp_grc) == ng
        @test length(a.forc_pco2_grc) == ng
        @test length(a.forc_ndep_grc) == ng
        @test length(a.forc_pc13o2_grc) == ng
        @test length(a.forc_po2_grc) == ng
        @test length(a.forc_pch4_grc) == ng
        @test length(a.forc_o3_grc) == ng
        @test length(a.forc_solar_not_downscaled_grc) == ng

        # Gridcell-level not-downscaled
        @test length(a.forc_t_not_downscaled_grc) == ng
        @test length(a.forc_pbot_not_downscaled_grc) == ng
        @test length(a.forc_th_not_downscaled_grc) == ng
        @test length(a.forc_rho_not_downscaled_grc) == ng
        @test length(a.forc_lwrad_not_downscaled_grc) == ng

        # 2D gridcell arrays
        @test size(a.forc_solad_not_downscaled_grc) == (ng, CLM.NUMRAD)
        @test size(a.forc_solai_grc) == (ng, CLM.NUMRAD)
        @test size(a.forc_aer_grc) == (ng, 14)
        @test all(x -> x == 0.0, a.forc_aer_grc)

        # Column-level downscaled
        @test length(a.forc_t_downscaled_col) == nc
        @test all(x -> x == 0.0, a.forc_t_downscaled_col)
        @test length(a.forc_pbot_downscaled_col) == nc
        @test length(a.forc_th_downscaled_col) == nc
        @test length(a.forc_rho_downscaled_col) == nc
        @test length(a.forc_lwrad_downscaled_col) == nc
        @test size(a.forc_solad_downscaled_col) == (nc, CLM.NUMRAD)
        @test length(a.forc_solar_downscaled_col) == nc

        # Patch-level time-averaged
        @test length(a.fsd24_patch) == np
        @test all(isnan, a.fsd24_patch)
        @test length(a.fsd240_patch) == np
        @test all(isnan, a.fsd240_patch)
        @test length(a.fsi24_patch) == np
        @test all(isnan, a.fsi24_patch)
        @test length(a.fsi240_patch) == np
        @test all(isnan, a.fsi240_patch)
        @test length(a.t_mo_patch) == np
        @test all(isnan, a.t_mo_patch)
        @test length(a.t_mo_min_patch) == np
        @test all(x -> x == CLM.SPVAL, a.t_mo_min_patch)
    end

    @testset "atm2lnd_init_for_testing! default params" begin
        ng, nc, np = 3, 5, 7
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init_for_testing!(a, ng, nc, np)

        @test a.params.repartition_rain_snow == false
        @test a.params.glcmec_downscale_longwave == false
        @test a.params.lapse_rate ≈ 0.01
        @test length(a.forc_u_grc) == ng
    end

    @testset "atm2lnd_init_for_testing! with custom params" begin
        ng, nc, np = 2, 3, 4
        custom_params = CLM.Atm2LndParamsData(
            repartition_rain_snow = true,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.005,
            precip_repartition_glc_all_snow_t = 270.0,
            precip_repartition_glc_frac_rain_slope = 0.5)

        a = CLM.Atm2LndData()
        CLM.atm2lnd_init_for_testing!(a, ng, nc, np; params = custom_params)

        @test a.params.repartition_rain_snow == true
        @test a.params.lapse_rate ≈ 0.005
        @test a.params.precip_repartition_glc_all_snow_t ≈ 270.0
        @test length(a.forc_u_grc) == ng
    end

    @testset "atm2lnd_read_namelist!" begin
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, 3, 5, 7)
        CLM.atm2lnd_read_namelist!(a;
            repartition_rain_snow = false,
            glcmec_downscale_longwave = false,
            lapse_rate = 0.008)
        @test a.params.lapse_rate ≈ 0.008
        @test a.params.repartition_rain_snow == false
    end

    @testset "atm2lnd_update_acc_vars! basic" begin
        ng, nc, np = 2, 4, 6
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, ng, nc, np)

        # Set some forcing values
        a.forc_solad_not_downscaled_grc[1, 1] = 100.0
        a.forc_solad_not_downscaled_grc[2, 1] = 200.0
        a.forc_solai_grc[1, 1] = 50.0
        a.forc_solai_grc[2, 1] = 75.0

        # patches 1-3 on gridcell 1, patches 4-6 on gridcell 2
        patch_gridcell = [1, 1, 1, 2, 2, 2]
        patch_column   = [1, 2, 3, 4, 4, 4]

        CLM.atm2lnd_update_acc_vars!(a, 1:np, patch_gridcell, patch_column)

        # Check direct beam
        @test a.fsd24_patch[1] ≈ 100.0
        @test a.fsd240_patch[3] ≈ 100.0
        @test a.fsd24_patch[4] ≈ 200.0
        @test a.fsd240_patch[6] ≈ 200.0

        # Check diffuse
        @test a.fsi24_patch[1] ≈ 50.0
        @test a.fsi240_patch[2] ≈ 50.0
        @test a.fsi24_patch[5] ≈ 75.0
        @test a.fsi240_patch[6] ≈ 75.0
    end

    @testset "atm2lnd_clean!" begin
        ng, nc, np = 3, 5, 7
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, ng, nc, np)

        @test length(a.forc_u_grc) == ng

        CLM.atm2lnd_clean!(a)

        @test length(a.forc_u_grc) == 0
        @test length(a.forc_t_downscaled_col) == 0
        @test length(a.fsd24_patch) == 0
        @test size(a.forc_solad_not_downscaled_grc) == (0, 0)
        @test size(a.forc_aer_grc) == (0, 0)
        @test length(a.t_mo_patch) == 0
        @test length(a.t_mo_min_patch) == 0
    end

    @testset "stub functions run without error" begin
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, 3, 5, 7)
        @test CLM.atm2lnd_init_history!(a, 3, 5, 7) === nothing
        @test CLM.atm2lnd_init_acc_buffer!(a) === nothing
        @test CLM.atm2lnd_init_acc_vars!(a, 1:7) === nothing
        @test CLM.atm2lnd_restart!(a, 1:7) === nothing
    end

    @testset "field mutability" begin
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, 3, 5, 7)

        a.forc_u_grc[1] = 5.0
        @test a.forc_u_grc[1] == 5.0

        a.forc_t_downscaled_col[2] = 300.0
        @test a.forc_t_downscaled_col[2] == 300.0

        a.forc_aer_grc[1, 3] = 0.001
        @test a.forc_aer_grc[1, 3] == 0.001
    end

    @testset "re-init overwrites previous state" begin
        a = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a, 3, 5, 7)
        a.forc_u_grc[1] = 999.0

        CLM.atm2lnd_init!(a, 10, 20, 30)
        @test length(a.forc_u_grc) == 10
        @test all(x -> x == 0.0, a.forc_u_grc)
        @test length(a.forc_t_downscaled_col) == 20
        @test length(a.fsd24_patch) == 30
    end

end
