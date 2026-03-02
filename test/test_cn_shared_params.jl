@testset "CNSharedParamsData" begin

    @testset "default construction" begin
        p = CLM.CNSharedParamsData()
        @test p.Q10 == 1.5
        @test p.minpsi == -10.0
        @test p.maxpsi == -0.1
        @test p.rf_cwdl2 == 0.0
        @test p.tau_cwd == 0.0
        @test p.cwd_flig == 0.0
        @test p.froz_q10 == 1.5
        @test p.decomp_depth_efolding == 0.5
        @test p.mino2lim == 0.0
        @test p.organic_max == 0.0
        @test p.constrain_stress_deciduous_onset == false
        @test p.use_matrixcn == false
        @test p.use_fun == false
        @test p.nlev_soildecomp_standard == 5
        @test p.upper_soil_layer == -1
    end

    @testset "cn_shared_params_read_netcdf!" begin
        p = CLM.CNSharedParamsData()
        CLM.cn_shared_params_read_netcdf!(p;
            q10_mr=2.0,
            minpsi_hr=-8.0,
            maxpsi_hr=-0.05,
            rf_cwdl2=0.3,
            tau_cwd=3.33,
            cwd_flig=0.24,
            decomp_depth_efolding=1.0,
            froz_q10=1.2,
            mino2lim=0.2,
            organic_max=130.0)

        @test p.Q10 == 2.0
        @test p.minpsi == -8.0
        @test p.maxpsi == -0.05
        @test p.rf_cwdl2 == 0.3
        @test p.tau_cwd == 3.33
        @test p.cwd_flig == 0.24
        @test p.decomp_depth_efolding == 1.0
        @test p.froz_q10 == 1.2
        @test p.mino2lim == 0.2
        @test p.organic_max == 130.0
    end

    @testset "cn_shared_params_read_namelist!" begin
        p = CLM.CNSharedParamsData()

        # Default (false)
        CLM.cn_shared_params_read_namelist!(p)
        @test p.constrain_stress_deciduous_onset == false

        # Set to true
        CLM.cn_shared_params_read_namelist!(p; constrain_stress_deciduous_onset=true)
        @test p.constrain_stress_deciduous_onset == true
    end

    @testset "cn_shared_params_read! combined" begin
        p = CLM.CNSharedParamsData()
        CLM.cn_shared_params_read!(p;
            q10_mr=1.8,
            minpsi_hr=-5.0,
            maxpsi_hr=-0.02,
            rf_cwdl2=0.25,
            tau_cwd=2.5,
            cwd_flig=0.20,
            decomp_depth_efolding=0.8,
            froz_q10=1.1,
            mino2lim=0.15,
            organic_max=100.0,
            constrain_stress_deciduous_onset=true)

        # NetCDF params
        @test p.Q10 == 1.8
        @test p.minpsi == -5.0
        @test p.maxpsi == -0.02
        @test p.rf_cwdl2 == 0.25
        @test p.tau_cwd == 2.5
        @test p.cwd_flig == 0.20
        @test p.decomp_depth_efolding == 0.8
        @test p.froz_q10 == 1.1
        @test p.mino2lim == 0.15
        @test p.organic_max == 100.0

        # Namelist param
        @test p.constrain_stress_deciduous_onset == true
    end

    @testset "cn_shared_params_set_soil_depth!" begin
        p = CLM.CNSharedParamsData()
        @test p.upper_soil_layer == -1

        # Provide a mock depth-to-layer function
        mock_find_layer = depth -> begin
            # Simple mock: layer 2 contains depth 0.12
            if depth <= 0.05
                return 1
            elseif depth <= 0.15
                return 2
            else
                return 3
            end
        end

        CLM.cn_shared_params_set_soil_depth!(p, mock_find_layer)
        @test p.upper_soil_layer == 2  # 0.12 falls in layer 2

        # Test with a different mock
        mock_find_layer2 = depth -> 4
        CLM.cn_shared_params_set_soil_depth!(p, mock_find_layer2)
        @test p.upper_soil_layer == 4
    end

    @testset "module-level flags mutability" begin
        p = CLM.CNSharedParamsData()

        p.use_matrixcn = true
        @test p.use_matrixcn == true

        p.use_fun = true
        @test p.use_fun == true

        p.nlev_soildecomp_standard = 10
        @test p.nlev_soildecomp_standard == 10
    end

    @testset "field mutability" begin
        p = CLM.CNSharedParamsData()

        p.Q10 = 3.0
        @test p.Q10 == 3.0

        p.minpsi = -20.0
        @test p.minpsi == -20.0

        p.organic_max = 200.0
        @test p.organic_max == 200.0
    end

end
