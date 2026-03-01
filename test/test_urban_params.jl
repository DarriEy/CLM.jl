@testset "UrbanParamsData" begin

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    numrad  = CLM.NUMRAD
    nlevurb = CLM.varpar.nlevurb

    # ---------------------------------------------------------------
    # Urban control and namelist functions
    # ---------------------------------------------------------------
    @testset "UrbanControl defaults" begin
        ctrl = CLM.UrbanControl()
        @test ctrl.urban_hac == CLM.URBAN_HAC_OFF
        @test ctrl.urban_explicit_ac == true
        @test ctrl.urban_traffic == false
        @test ctrl.building_temp_method == CLM.BUILDING_TEMP_METHOD_PROG
        @test ctrl.read_namelist == false
    end

    @testset "is_simple/prog_build_temp before namelist read" begin
        saved = CLM.urban_ctrl.read_namelist
        CLM.urban_ctrl.read_namelist = false
        @test_throws ErrorException CLM.is_simple_build_temp()
        @test_throws ErrorException CLM.is_prog_build_temp()
        CLM.urban_ctrl.read_namelist = saved
    end

    @testset "urban_read_nml!" begin
        ctrl = CLM.UrbanControl()
        CLM.urban_read_nml!(ctrl; urban_hac=CLM.URBAN_HAC_ON,
                            building_temp_method=CLM.BUILDING_TEMP_METHOD_SIMPLE)
        @test ctrl.read_namelist == true
        @test ctrl.urban_hac == CLM.URBAN_HAC_ON
        @test ctrl.building_temp_method == CLM.BUILDING_TEMP_METHOD_SIMPLE

        # Traffic should throw
        @test_throws ErrorException CLM.urban_read_nml!(ctrl; urban_traffic=true)
    end

    @testset "is_simple/prog_build_temp after namelist" begin
        saved_method = CLM.urban_ctrl.building_temp_method
        saved_read   = CLM.urban_ctrl.read_namelist

        CLM.urban_ctrl.read_namelist = true
        CLM.urban_ctrl.building_temp_method = CLM.BUILDING_TEMP_METHOD_PROG
        @test CLM.is_prog_build_temp() == true
        @test CLM.is_simple_build_temp() == false

        CLM.urban_ctrl.building_temp_method = CLM.BUILDING_TEMP_METHOD_SIMPLE
        @test CLM.is_simple_build_temp() == true
        @test CLM.is_prog_build_temp() == false

        CLM.urban_ctrl.building_temp_method = saved_method
        CLM.urban_ctrl.read_namelist = saved_read
    end

    # ---------------------------------------------------------------
    # UrbanInputData
    # ---------------------------------------------------------------
    @testset "UrbanInputData default construction" begin
        ui = CLM.UrbanInputData()
        @test size(ui.canyon_hwr) == (0, 0)
        @test size(ui.alb_roof_dir) == (0, 0, 0)
        @test size(ui.tk_wall) == (0, 0, 0)
        @test size(ui.nlev_improad) == (0, 0)
    end

    @testset "urbinp_init!" begin
        ng = 5
        numurbl = CLM.NUMURBL
        ui = CLM.UrbanInputData()
        CLM.urbinp_init!(ui, ng, numurbl, numrad, nlevurb)

        @test size(ui.canyon_hwr) == (ng, numurbl)
        @test size(ui.wtlunit_roof) == (ng, numurbl)
        @test size(ui.em_roof) == (ng, numurbl)
        @test size(ui.alb_roof_dir) == (ng, numurbl, numrad)
        @test size(ui.alb_wall_dif) == (ng, numurbl, numrad)
        @test size(ui.tk_wall) == (ng, numurbl, nlevurb)
        @test size(ui.cv_improad) == (ng, numurbl, nlevurb)
        @test size(ui.thick_wall) == (ng, numurbl)
        @test size(ui.nlev_improad) == (ng, numurbl)
        @test size(ui.t_building_min) == (ng, numurbl)
        @test size(ui.ht_roof) == (ng, numurbl)

        # Check NaN initialization for floats
        @test all(isnan, ui.canyon_hwr)
        @test all(isnan, ui.alb_roof_dir)
        @test all(isnan, ui.tk_wall)

        # Check integer init for nlev_improad
        @test all(x -> x == 0, ui.nlev_improad)
    end

    @testset "urbinp_clean!" begin
        ui = CLM.UrbanInputData()
        CLM.urbinp_init!(ui, 3, 3, 2, 5)
        CLM.urbinp_clean!(ui)

        @test size(ui.canyon_hwr) == (0, 0)
        @test size(ui.alb_roof_dir) == (0, 0, 0)
        @test size(ui.tk_wall) == (0, 0, 0)
        @test size(ui.nlev_improad) == (0, 0)
    end

    # ---------------------------------------------------------------
    # UrbanParamsData
    # ---------------------------------------------------------------
    @testset "UrbanParamsData default construction" begin
        up = CLM.UrbanParamsData()
        @test length(up.wind_hgt_canyon) == 0
        @test length(up.em_roof) == 0
        @test length(up.vf_sr) == 0
        @test size(up.alb_roof_dir) == (0, 0)
        @test size(up.tk_wall) == (0, 0)
    end

    @testset "urbanparams_init!" begin
        nl = 6
        up = CLM.UrbanParamsData()
        CLM.urbanparams_init!(up, nl; nlevurb=nlevurb)

        # 1D fields
        @test length(up.wind_hgt_canyon) == nl
        @test length(up.em_roof) == nl
        @test length(up.em_improad) == nl
        @test length(up.em_perroad) == nl
        @test length(up.em_wall) == nl
        @test length(up.t_building_min) == nl
        @test length(up.thick_wall) == nl
        @test length(up.thick_roof) == nl
        @test length(up.nlev_improad) == nl
        @test length(up.vf_sr) == nl
        @test length(up.vf_wr) == nl
        @test length(up.vf_sw) == nl
        @test length(up.vf_rw) == nl
        @test length(up.vf_ww) == nl
        @test length(up.eflx_traffic_factor) == nl

        # 2D albedo fields
        @test size(up.alb_roof_dir) == (nl, numrad)
        @test size(up.alb_wall_dif) == (nl, numrad)

        # 2D thermal property fields
        @test size(up.tk_wall) == (nl, nlevurb)
        @test size(up.tk_roof) == (nl, nlevurb)
        @test size(up.cv_wall) == (nl, nlevurb)
        @test size(up.cv_roof) == (nl, nlevurb)
        @test size(up.tk_improad) == (nl, nlevurb)
        @test size(up.cv_improad) == (nl, nlevurb)

        # NaN initialization
        @test all(isnan, up.wind_hgt_canyon)
        @test all(isnan, up.em_roof)
        @test all(isnan, up.alb_roof_dir)
        @test all(isnan, up.tk_wall)
        @test all(isnan, up.vf_sr)
        @test all(isnan, up.eflx_traffic_factor)

        # Integer sentinel
        @test all(x -> x == typemax(Int), up.nlev_improad)
    end

    @testset "urbanparams_init! with nlevurb=0" begin
        up = CLM.UrbanParamsData()
        CLM.urbanparams_init!(up, 4; nlevurb=0)

        @test size(up.tk_wall) == (0, 0)
        @test size(up.tk_roof) == (0, 0)
        @test size(up.cv_wall) == (0, 0)
        @test size(up.cv_roof) == (0, 0)
        # tk_improad and cv_improad are always allocated
        @test size(up.tk_improad) == (4, 0)
        @test size(up.cv_improad) == (4, 0)
    end

    @testset "urbanparams_clean!" begin
        up = CLM.UrbanParamsData()
        CLM.urbanparams_init!(up, 5; nlevurb=nlevurb)
        CLM.urbanparams_clean!(up)

        @test length(up.wind_hgt_canyon) == 0
        @test length(up.em_roof) == 0
        @test length(up.vf_sr) == 0
        @test size(up.alb_roof_dir) == (0, 0)
        @test size(up.tk_wall) == (0, 0)
        @test length(up.nlev_improad) == 0
        @test length(up.eflx_traffic_factor) == 0
    end

    # ---------------------------------------------------------------
    # urbanparams_populate! — integration test with mock data
    # ---------------------------------------------------------------
    @testset "urbanparams_populate! view factors and aerodynamics" begin
        nl = 3   # 3 landunits: 2 urban, 1 non-urban
        ng = 2   # 2 gridcells
        numurbl = CLM.NUMURBL

        # Set up LandunitData
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)

        # Landunit 1: urban TBD (itype=7), gridcell 1
        lun.urbpoi[1]   = true
        lun.itype[1]    = CLM.ISTURB_TBD  # 7
        lun.gridcell[1] = 1

        # Landunit 2: urban HD (itype=8), gridcell 2
        lun.urbpoi[2]   = true
        lun.itype[2]    = CLM.ISTURB_HD  # 8
        lun.gridcell[2] = 2

        # Landunit 3: non-urban
        lun.urbpoi[3]   = false
        lun.itype[3]    = CLM.ISTSOIL
        lun.gridcell[3] = 1

        # Set up UrbanInputData
        ui = CLM.UrbanInputData()
        CLM.urbinp_init!(ui, ng, numurbl, numrad, nlevurb)

        # Fill input data with valid values for gridcell 1, density type 1 (TBD)
        dindx1 = 1  # ISTURB_TBD - ISTURB_MIN + 1 = 7 - 7 + 1 = 1
        ui.canyon_hwr[1, dindx1]      = 1.0
        ui.wtlunit_roof[1, dindx1]    = 0.5
        ui.wtroad_perv[1, dindx1]     = 0.4
        ui.em_roof[1, dindx1]         = 0.9
        ui.em_improad[1, dindx1]      = 0.95
        ui.em_perroad[1, dindx1]      = 0.95
        ui.em_wall[1, dindx1]         = 0.85
        ui.ht_roof[1, dindx1]         = 10.0
        ui.wind_hgt_canyon[1, dindx1]  = 5.0
        ui.thick_wall[1, dindx1]      = 0.15
        ui.thick_roof[1, dindx1]      = 0.15
        ui.nlev_improad[1, dindx1]    = 3
        ui.t_building_min[1, dindx1]  = 288.0
        for ib in 1:numrad
            ui.alb_roof_dir[1, dindx1, ib]    = 0.2
            ui.alb_roof_dif[1, dindx1, ib]    = 0.2
            ui.alb_improad_dir[1, dindx1, ib] = 0.1
            ui.alb_improad_dif[1, dindx1, ib] = 0.1
            ui.alb_perroad_dir[1, dindx1, ib] = 0.15
            ui.alb_perroad_dif[1, dindx1, ib] = 0.15
            ui.alb_wall_dir[1, dindx1, ib]    = 0.25
            ui.alb_wall_dif[1, dindx1, ib]    = 0.25
        end
        for k in 1:nlevurb
            ui.tk_wall[1, dindx1, k]    = 1.0
            ui.tk_roof[1, dindx1, k]    = 0.8
            ui.tk_improad[1, dindx1, k] = 1.2
            ui.cv_wall[1, dindx1, k]    = 1.5e6
            ui.cv_roof[1, dindx1, k]    = 1.2e6
            ui.cv_improad[1, dindx1, k] = 1.8e6
        end

        # Fill input data for gridcell 2, density type 2 (HD)
        dindx2 = 2  # ISTURB_HD - ISTURB_MIN + 1 = 8 - 7 + 1 = 2
        ui.canyon_hwr[2, dindx2]      = 2.0
        ui.wtlunit_roof[2, dindx2]    = 0.6
        ui.wtroad_perv[2, dindx2]     = 0.3
        ui.em_roof[2, dindx2]         = 0.92
        ui.em_improad[2, dindx2]      = 0.93
        ui.em_perroad[2, dindx2]      = 0.94
        ui.em_wall[2, dindx2]         = 0.88
        ui.ht_roof[2, dindx2]         = 20.0
        ui.wind_hgt_canyon[2, dindx2]  = 6.0
        ui.thick_wall[2, dindx2]      = 0.2
        ui.thick_roof[2, dindx2]      = 0.2
        ui.nlev_improad[2, dindx2]    = 4
        ui.t_building_min[2, dindx2]  = 293.0
        for ib in 1:numrad
            ui.alb_roof_dir[2, dindx2, ib]    = 0.3
            ui.alb_roof_dif[2, dindx2, ib]    = 0.3
            ui.alb_improad_dir[2, dindx2, ib] = 0.12
            ui.alb_improad_dif[2, dindx2, ib] = 0.12
            ui.alb_perroad_dir[2, dindx2, ib] = 0.18
            ui.alb_perroad_dif[2, dindx2, ib] = 0.18
            ui.alb_wall_dir[2, dindx2, ib]    = 0.28
            ui.alb_wall_dif[2, dindx2, ib]    = 0.28
        end
        for k in 1:nlevurb
            ui.tk_wall[2, dindx2, k]    = 1.1
            ui.tk_roof[2, dindx2, k]    = 0.9
            ui.tk_improad[2, dindx2, k] = 1.3
            ui.cv_wall[2, dindx2, k]    = 1.6e6
            ui.cv_roof[2, dindx2, k]    = 1.3e6
            ui.cv_improad[2, dindx2, k] = 1.9e6
        end

        # Set up UrbanParamsData and populate
        up = CLM.UrbanParamsData()
        CLM.urbanparams_init!(up, nl; nlevurb=nlevurb)
        CLM.urbanparams_populate!(up, lun, ui, 1:nl)

        # Check view factors for canyon_hwr = 1.0 (landunit 1)
        # vf_sr = sqrt(1^2 + 1) - 1 = sqrt(2) - 1 ≈ 0.4142
        @test up.vf_sr[1] ≈ sqrt(2.0) - 1.0 atol=1e-10
        @test up.vf_wr[1] ≈ 0.5 * (1.0 - up.vf_sr[1]) atol=1e-10

        # vf_sw = 0.5 * (1 + 1 - sqrt(2)) / 1 ≈ 0.2929
        expected_vf_sw = 0.5 * (1.0 + 1.0 - sqrt(2.0)) / 1.0
        @test up.vf_sw[1] ≈ expected_vf_sw atol=1e-10
        @test up.vf_rw[1] ≈ up.vf_sw[1] atol=1e-10
        @test up.vf_ww[1] ≈ 1.0 - up.vf_sw[1] - up.vf_rw[1] atol=1e-10

        # View factors should sum to 1
        @test up.vf_sr[1] + 2.0 * up.vf_wr[1] ≈ 1.0 atol=1e-06
        @test up.vf_sw[1] + up.vf_rw[1] + up.vf_ww[1] ≈ 1.0 atol=1e-06

        # View factors for canyon_hwr = 2.0 (landunit 2)
        @test up.vf_sr[2] ≈ sqrt(4.0 + 1.0) - 2.0 atol=1e-10
        @test up.vf_sr[2] + 2.0 * up.vf_wr[2] ≈ 1.0 atol=1e-06
        @test up.vf_sw[2] + up.vf_rw[2] + up.vf_ww[2] ≈ 1.0 atol=1e-06

        # Check emissivities were copied
        @test up.em_roof[1] == 0.9
        @test up.em_wall[1] == 0.85
        @test up.em_roof[2] == 0.92

        # Check albedos copied
        @test up.alb_roof_dir[1, 1] == 0.2
        @test up.alb_roof_dir[2, 1] == 0.3

        # Check thermal properties copied
        @test up.tk_wall[1, 1] == 1.0
        @test up.tk_roof[2, 1] == 0.9

        # Check landunit properties were set
        @test lun.canyon_hwr[1] == 1.0
        @test lun.canyon_hwr[2] == 2.0
        @test lun.ht_roof[1] == 10.0
        @test lun.ht_roof[2] == 20.0

        # Check traffic factor is 0 when urban_traffic=false (default)
        @test up.eflx_traffic_factor[1] == 0.0
        @test up.eflx_traffic_factor[2] == 0.0

        # Check non-urban landunit gets SPVAL
        @test up.eflx_traffic_factor[3] == CLM.SPVAL
        @test up.t_building_min[3] == CLM.SPVAL
        @test up.vf_sr[3] == CLM.SPVAL
        @test up.vf_wr[3] == CLM.SPVAL
        @test up.vf_sw[3] == CLM.SPVAL
        @test up.vf_rw[3] == CLM.SPVAL
        @test up.vf_ww[3] == CLM.SPVAL

        # Check z_d_town and z_0_town were computed (non-NaN)
        @test !isnan(lun.z_d_town[1])
        @test !isnan(lun.z_0_town[1])
        @test lun.z_d_town[1] > 0.0
        @test lun.z_0_town[1] > 0.0

        # With urban_hac=OFF, t_building_min should be overridden to 200.0
        @test up.t_building_min[1] == 200.0
        @test up.t_building_min[2] == 200.0
    end

    @testset "urbanparams_populate! with urban_traffic" begin
        nl = 1
        ng = 1
        numurbl = CLM.NUMURBL

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        lun.urbpoi[1]   = true
        lun.itype[1]    = CLM.ISTURB_TBD
        lun.gridcell[1] = 1

        ui = CLM.UrbanInputData()
        CLM.urbinp_init!(ui, ng, numurbl, numrad, nlevurb)
        ui.canyon_hwr[1, 1]      = 1.5
        ui.wtlunit_roof[1, 1]   = 0.5
        ui.wtroad_perv[1, 1]    = 0.4
        ui.em_roof[1, 1]        = 0.9
        ui.em_improad[1, 1]     = 0.95
        ui.em_perroad[1, 1]     = 0.95
        ui.em_wall[1, 1]        = 0.85
        ui.ht_roof[1, 1]        = 15.0
        ui.wind_hgt_canyon[1, 1] = 5.0
        ui.thick_wall[1, 1]     = 0.15
        ui.thick_roof[1, 1]     = 0.15
        ui.nlev_improad[1, 1]   = 3
        ui.t_building_min[1, 1] = 288.0
        for ib in 1:numrad
            ui.alb_roof_dir[1, 1, ib]    = 0.2
            ui.alb_roof_dif[1, 1, ib]    = 0.2
            ui.alb_improad_dir[1, 1, ib] = 0.1
            ui.alb_improad_dif[1, 1, ib] = 0.1
            ui.alb_perroad_dir[1, 1, ib] = 0.15
            ui.alb_perroad_dif[1, 1, ib] = 0.15
            ui.alb_wall_dir[1, 1, ib]    = 0.25
            ui.alb_wall_dif[1, 1, ib]    = 0.25
        end
        for k in 1:nlevurb
            ui.tk_wall[1, 1, k]    = 1.0
            ui.tk_roof[1, 1, k]    = 0.8
            ui.tk_improad[1, 1, k] = 1.2
            ui.cv_wall[1, 1, k]    = 1.5e6
            ui.cv_roof[1, 1, k]    = 1.2e6
            ui.cv_improad[1, 1, k] = 1.8e6
        end

        up = CLM.UrbanParamsData()
        CLM.urbanparams_init!(up, nl; nlevurb=nlevurb)
        CLM.urbanparams_populate!(up, lun, ui, 1:nl; urban_traffic=true)

        # Traffic factor = 3.6 * (1.5 - 0.5) + 1.0 = 4.6
        @test up.eflx_traffic_factor[1] ≈ 4.6 atol=1e-10
    end

    @testset "urbanparams_populate! with use_vancouver" begin
        nl = 1
        ng = 1
        numurbl = CLM.NUMURBL

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        lun.urbpoi[1]   = true
        lun.itype[1]    = CLM.ISTURB_TBD
        lun.gridcell[1] = 1

        ui = CLM.UrbanInputData()
        CLM.urbinp_init!(ui, ng, numurbl, numrad, nlevurb)
        ui.canyon_hwr[1, 1]      = 1.0
        ui.wtlunit_roof[1, 1]   = 0.5
        ui.wtroad_perv[1, 1]    = 0.4
        ui.em_roof[1, 1]        = 0.9
        ui.em_improad[1, 1]     = 0.95
        ui.em_perroad[1, 1]     = 0.95
        ui.em_wall[1, 1]        = 0.85
        ui.ht_roof[1, 1]        = 10.0
        ui.wind_hgt_canyon[1, 1] = 5.0
        ui.thick_wall[1, 1]     = 0.15
        ui.thick_roof[1, 1]     = 0.15
        ui.nlev_improad[1, 1]   = 3
        ui.t_building_min[1, 1] = 288.0
        for ib in 1:numrad
            ui.alb_roof_dir[1, 1, ib]    = 0.2
            ui.alb_roof_dif[1, 1, ib]    = 0.2
            ui.alb_improad_dir[1, 1, ib] = 0.1
            ui.alb_improad_dif[1, 1, ib] = 0.1
            ui.alb_perroad_dir[1, 1, ib] = 0.15
            ui.alb_perroad_dif[1, 1, ib] = 0.15
            ui.alb_wall_dir[1, 1, ib]    = 0.25
            ui.alb_wall_dif[1, 1, ib]    = 0.25
        end
        for k in 1:nlevurb
            ui.tk_wall[1, 1, k]    = 1.0
            ui.tk_roof[1, 1, k]    = 0.8
            ui.tk_improad[1, 1, k] = 1.2
            ui.cv_wall[1, 1, k]    = 1.5e6
            ui.cv_roof[1, 1, k]    = 1.2e6
            ui.cv_improad[1, 1, k] = 1.8e6
        end

        up = CLM.UrbanParamsData()
        CLM.urbanparams_init!(up, nl; nlevurb=nlevurb)
        CLM.urbanparams_populate!(up, lun, ui, 1:nl; use_vancouver=true)

        # Vancouver special values
        @test lun.z_d_town[1] == 3.5
        @test lun.z_0_town[1] == 0.35
        @test up.t_building_min[1] == 200.0
    end

    # ---------------------------------------------------------------
    # check_urban
    # ---------------------------------------------------------------
    @testset "check_urban valid data" begin
        ng = 2
        numurbl = CLM.NUMURBL
        ui = CLM.UrbanInputData()
        CLM.urbinp_init!(ui, ng, numurbl, numrad, nlevurb)

        # Fill with valid data
        for g in 1:ng, d in 1:numurbl
            ui.canyon_hwr[g, d]      = 1.0
            ui.wtlunit_roof[g, d]    = 0.5
            ui.wtroad_perv[g, d]     = 0.4
            ui.em_roof[g, d]         = 0.9
            ui.em_improad[g, d]      = 0.95
            ui.em_perroad[g, d]      = 0.95
            ui.em_wall[g, d]         = 0.85
            ui.ht_roof[g, d]         = 10.0
            ui.wind_hgt_canyon[g, d]  = 5.0
            ui.thick_wall[g, d]      = 0.15
            ui.thick_roof[g, d]      = 0.15
            ui.nlev_improad[g, d]    = 2
            ui.t_building_min[g, d]  = 288.0
            for ib in 1:numrad
                ui.alb_roof_dir[g, d, ib]    = 0.2
                ui.alb_roof_dif[g, d, ib]    = 0.2
                ui.alb_improad_dir[g, d, ib] = 0.1
                ui.alb_improad_dif[g, d, ib] = 0.1
                ui.alb_perroad_dir[g, d, ib] = 0.15
                ui.alb_perroad_dif[g, d, ib] = 0.15
                ui.alb_wall_dir[g, d, ib]    = 0.25
                ui.alb_wall_dif[g, d, ib]    = 0.25
            end
            for k in 1:nlevurb
                ui.tk_wall[g, d, k]    = 1.0
                ui.tk_roof[g, d, k]    = 0.8
                ui.tk_improad[g, d, k] = 1.2
                ui.cv_wall[g, d, k]    = 1.5e6
                ui.cv_roof[g, d, k]    = 1.2e6
                ui.cv_improad[g, d, k] = 1.8e6
            end
        end

        pcturb = fill(10.0, ng, numurbl)
        # Should not throw
        @test CLM.check_urban(ui, pcturb, 1:ng; caller="test") === nothing
    end

    @testset "check_urban invalid data" begin
        ng = 1
        numurbl = CLM.NUMURBL
        ui = CLM.UrbanInputData()
        CLM.urbinp_init!(ui, ng, numurbl, numrad, nlevurb)

        # Set canyon_hwr to 0 (invalid) for density type 1
        ui.canyon_hwr[1, 1] = 0.0

        pcturb = fill(10.0, ng, numurbl)
        @test_throws ErrorException CLM.check_urban(ui, pcturb, 1:ng; caller="test")
    end

    @testset "check_urban zero pct urban skips validation" begin
        ng = 1
        numurbl = CLM.NUMURBL
        ui = CLM.UrbanInputData()
        CLM.urbinp_init!(ui, ng, numurbl, numrad, nlevurb)

        # All input data is NaN (invalid), but pct urban is 0 — should pass
        pcturb = fill(0.0, ng, numurbl)
        @test CLM.check_urban(ui, pcturb, 1:ng; caller="test") === nothing
    end

    # ---------------------------------------------------------------
    # Field mutability
    # ---------------------------------------------------------------
    @testset "field mutability" begin
        up = CLM.UrbanParamsData()
        CLM.urbanparams_init!(up, 4; nlevurb=nlevurb)

        up.em_roof[2] = 0.91
        @test up.em_roof[2] == 0.91

        up.alb_roof_dir[1, 2] = 0.33
        @test up.alb_roof_dir[1, 2] == 0.33

        up.tk_wall[3, 1] = 2.5
        @test up.tk_wall[3, 1] == 2.5

        up.vf_sr[1] = 0.42
        @test up.vf_sr[1] == 0.42
    end

    @testset "re-init overwrites previous state" begin
        up = CLM.UrbanParamsData()
        CLM.urbanparams_init!(up, 3; nlevurb=nlevurb)
        up.em_roof[1] = 999.0

        CLM.urbanparams_init!(up, 6; nlevurb=nlevurb)
        @test length(up.em_roof) == 6
        @test all(isnan, up.em_roof)
    end

end
