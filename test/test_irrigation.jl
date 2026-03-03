@testset "Irrigation (IrrigationMod)" begin

    # Ensure varpar is initialized for column_init!, waterflux_init!, etc.
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    # -----------------------------------------------------------------------
    # Test constants
    # -----------------------------------------------------------------------
    @testset "Irrigation constants" begin
        @test CLM.WILTING_POINT_SMP == -150000.0
        @test CLM.M3_OVER_KM2_TO_MM == 1.0e-3
        @test CLM.IRRIG_METHOD_UNSET == 0
        @test CLM.IRRIG_METHOD_DRIP == 1
        @test CLM.IRRIG_METHOD_SPRINKLER == 2
    end

    # -----------------------------------------------------------------------
    # Test IrrigationParamsData construction
    # -----------------------------------------------------------------------
    @testset "IrrigationParamsData construction" begin
        params = CLM.IrrigationParamsData()
        @test params.irrig_min_lai == 0.0
        @test params.irrig_start_time == 21600
        @test params.irrig_length == 14400
        @test params.irrig_target_smp == -3400.0
        @test params.irrig_depth == 0.6
        @test params.irrig_threshold_fraction == 1.0
        @test params.limit_irrigation_if_rof_enabled == false
        @test params.use_groundwater_irrigation == false
        @test params.irrig_method_default == CLM.IRRIG_METHOD_DRIP

        # Custom construction
        params2 = CLM.IrrigationParamsData(
            irrig_min_lai = 0.5,
            irrig_target_smp = -5000.0,
            irrig_depth = 1.0
        )
        @test params2.irrig_min_lai == 0.5
        @test params2.irrig_target_smp == -5000.0
        @test params2.irrig_depth == 1.0
    end

    # -----------------------------------------------------------------------
    # Test IrrigationData construction and init/clean
    # -----------------------------------------------------------------------
    @testset "IrrigationData construction and init/clean" begin
        irrig = CLM.IrrigationData()
        @test irrig.dtime == 0
        @test irrig.irrig_nsteps_per_day == 0
        @test length(irrig.sfc_irrig_rate_patch) == 0
        @test length(irrig.irrig_method_patch) == 0
        @test size(irrig.relsat_wilting_point_col) == (0, 0)
        @test size(irrig.relsat_target_col) == (0, 0)

        # Test allocate
        np = 4
        nc = 2
        nlevsoi = 3
        CLM.irrigation_init_allocate!(irrig, np, nc, nlevsoi)

        @test length(irrig.qflx_irrig_demand_patch) == np
        @test all(isnan, irrig.qflx_irrig_demand_patch)
        @test size(irrig.relsat_wilting_point_col) == (nc, nlevsoi)
        @test all(isnan, irrig.relsat_wilting_point_col)
        @test size(irrig.relsat_target_col) == (nc, nlevsoi)
        @test length(irrig.sfc_irrig_rate_patch) == np
        @test length(irrig.irrig_rate_demand_patch) == np
        @test length(irrig.irrig_method_patch) == np
        @test all(x -> x == CLM.ISPVAL, irrig.irrig_method_patch)
        @test length(irrig.n_irrig_steps_left_patch) == np
        @test all(x -> x == 0, irrig.n_irrig_steps_left_patch)

        # Test clean
        CLM.irrigation_clean!(irrig)
        @test length(irrig.sfc_irrig_rate_patch) == 0
        @test length(irrig.irrig_method_patch) == 0
        @test size(irrig.relsat_wilting_point_col) == (0, 0)
        @test size(irrig.relsat_target_col) == (0, 0)
        @test length(irrig.qflx_irrig_demand_patch) == 0
    end

    # -----------------------------------------------------------------------
    # Test calc_irrig_nsteps_per_day
    # -----------------------------------------------------------------------
    @testset "calc_irrig_nsteps_per_day" begin
        # Exact multiple: 14400 sec / 1800 sec = 8 steps
        @test CLM.calc_irrig_nsteps_per_day(14400, 1800) == 8
        # Not exact multiple: 14400 / 3600 = 4
        @test CLM.calc_irrig_nsteps_per_day(14400, 3600) == 4
        # Round up: 14400 / 5000 = 2.88 -> 3
        @test CLM.calc_irrig_nsteps_per_day(14400, 5000) == 3
        # Edge: 1 sec / 1 sec = 1 step
        @test CLM.calc_irrig_nsteps_per_day(1, 1) == 1
        # Full day: 86400 / 1800 = 48
        @test CLM.calc_irrig_nsteps_per_day(86400, 1800) == 48
    end

    # -----------------------------------------------------------------------
    # Test relsat_to_h2osoi
    # -----------------------------------------------------------------------
    @testset "relsat_to_h2osoi" begin
        # Full saturation, 0.5 porosity, 0.1 m layer
        # vol_liq = 0.5 * 1.0 = 0.5
        # h2osoi = 0.5 * 1000 * 0.1 = 50.0 kg/m2
        h = CLM.relsat_to_h2osoi(1.0, 0.5, 0.1)
        @test h ≈ 50.0 atol=1e-10

        # Zero saturation
        h = CLM.relsat_to_h2osoi(0.0, 0.5, 0.1)
        @test h ≈ 0.0 atol=1e-10

        # Half saturation, 0.4 porosity, 0.2 m layer
        # vol_liq = 0.4 * 0.5 = 0.2
        # h2osoi = 0.2 * 1000 * 0.2 = 40.0
        h = CLM.relsat_to_h2osoi(0.5, 0.4, 0.2)
        @test h ≈ 40.0 atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test check_namelist_validity!
    # -----------------------------------------------------------------------
    @testset "check_namelist_validity!" begin
        # Valid params should not throw
        params = CLM.IrrigationParamsData(
            irrig_min_lai = 0.0,
            irrig_start_time = 21600,
            irrig_length = 14400,
            irrig_target_smp = -3400.0,
            irrig_depth = 0.6,
            irrig_threshold_fraction = 1.0,
            irrig_river_volume_threshold = 0.1,
            limit_irrigation_if_rof_enabled = false,
            use_groundwater_irrigation = false
        )
        @test CLM.check_namelist_validity!(params) === nothing

        # Negative irrig_min_lai should throw
        params_bad = CLM.IrrigationParamsData(irrig_min_lai = -1.0)
        @test_throws ErrorException CLM.check_namelist_validity!(params_bad)

        # irrig_start_time out of range
        params_bad2 = CLM.IrrigationParamsData(irrig_start_time = -1)
        @test_throws ErrorException CLM.check_namelist_validity!(params_bad2)

        params_bad3 = CLM.IrrigationParamsData(irrig_start_time = 86400)
        @test_throws ErrorException CLM.check_namelist_validity!(params_bad3)

        # irrig_length out of range
        params_bad4 = CLM.IrrigationParamsData(irrig_length = 0)
        @test_throws ErrorException CLM.check_namelist_validity!(params_bad4)

        # irrig_target_smp positive
        params_bad5 = CLM.IrrigationParamsData(irrig_target_smp = 100.0)
        @test_throws ErrorException CLM.check_namelist_validity!(params_bad5)

        # irrig_target_smp below wilting point
        params_bad6 = CLM.IrrigationParamsData(irrig_target_smp = -200000.0)
        @test_throws ErrorException CLM.check_namelist_validity!(params_bad6)

        # irrig_threshold_fraction out of range
        params_bad7 = CLM.IrrigationParamsData(irrig_threshold_fraction = -0.1)
        @test_throws ErrorException CLM.check_namelist_validity!(params_bad7)

        params_bad8 = CLM.IrrigationParamsData(irrig_threshold_fraction = 1.5)
        @test_throws ErrorException CLM.check_namelist_validity!(params_bad8)

        # use_groundwater without limit_irrigation
        params_bad9 = CLM.IrrigationParamsData(
            use_groundwater_irrigation = true,
            limit_irrigation_if_rof_enabled = false
        )
        @test_throws ErrorException CLM.check_namelist_validity!(params_bad9)

        # use_groundwater with aquifer layer
        params_gw = CLM.IrrigationParamsData(
            use_groundwater_irrigation = true,
            limit_irrigation_if_rof_enabled = true,
            irrig_river_volume_threshold = 0.1
        )
        @test_throws ErrorException CLM.check_namelist_validity!(params_gw, use_aquifer_layer=true)
    end

    # -----------------------------------------------------------------------
    # Test set_irrig_method!
    # -----------------------------------------------------------------------
    @testset "set_irrig_method!" begin
        irrig = CLM.IrrigationData()
        irrig.params = CLM.IrrigationParamsData(irrig_method_default = CLM.IRRIG_METHOD_DRIP)

        np = 3
        nc = 2
        ng = 1
        nlevsoi = 2
        CLM.irrigation_init_allocate!(irrig, np, nc, nlevsoi)

        # Setup patch data
        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        patch_data.gridcell .= [1, 1, 1]
        patch_data.itype .= [1, 2, 3]

        # Setup gridcell data
        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ng)

        # pftcon irrigated: type 1 and 2 are irrigated
        pftcon_irrigated = [1.0, 1.0, 0.0]

        # irrig_method_surface: [ng × max_pft_types]
        irrig_method_surface = fill(CLM.IRRIG_METHOD_UNSET, ng, 3)
        irrig_method_surface[1, 1] = CLM.IRRIG_METHOD_SPRINKLER
        irrig_method_surface[1, 2] = CLM.IRRIG_METHOD_UNSET

        bounds_p = 1:np

        CLM.set_irrig_method!(irrig, pftcon_irrigated, irrig_method_surface,
                              patch_data, grc, bounds_p)

        # Patch 1: irrigated, method = sprinkler (from surface data)
        @test irrig.irrig_method_patch[1] == CLM.IRRIG_METHOD_SPRINKLER
        # Patch 2: irrigated, method unset -> use default (drip)
        @test irrig.irrig_method_patch[2] == CLM.IRRIG_METHOD_DRIP
        # Patch 3: not irrigated -> default
        @test irrig.irrig_method_patch[3] == CLM.IRRIG_METHOD_DRIP
    end

    # -----------------------------------------------------------------------
    # Test p2c_irrig! — patch-to-column average
    # -----------------------------------------------------------------------
    @testset "p2c_irrig!" begin
        np = 4
        nc = 2
        bounds_c = 1:nc
        bounds_p = 1:np

        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        # Patches 1,2 on column 1; patches 3,4 on column 2
        patch_data.column .= [1, 1, 2, 2]
        patch_data.wtcol .= [0.6, 0.4, 0.3, 0.7]

        mask_soilc = trues(nc)
        patch_in = [10.0, 20.0, 30.0, 40.0]
        col_out = zeros(nc)

        CLM.p2c_irrig!(col_out, patch_in, patch_data, mask_soilc, bounds_c, bounds_p)

        # Column 1: 10*0.6 + 20*0.4 = 6 + 8 = 14.0
        @test col_out[1] ≈ 14.0 atol=1e-10
        # Column 2: 30*0.3 + 40*0.7 = 9 + 28 = 37.0
        @test col_out[2] ≈ 37.0 atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test c2g_irrig! — column-to-gridcell average
    # -----------------------------------------------------------------------
    @testset "c2g_irrig!" begin
        nc = 3
        ng = 2
        bounds_c = 1:nc
        bounds_g = 1:ng

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)
        col_data.gridcell .= [1, 1, 2]
        col_data.wtgcell .= [0.5, 0.5, 1.0]

        carr = [10.0, 20.0, 30.0]
        garr = zeros(ng)

        CLM.c2g_irrig!(garr, carr, col_data, bounds_c, bounds_g)

        # Gridcell 1: 10*0.5 + 20*0.5 = 5 + 10 = 15.0
        @test garr[1] ≈ 15.0 atol=1e-10
        # Gridcell 2: 30*1.0 = 30.0
        @test garr[2] ≈ 30.0 atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test calc_total_gw_uncon_irrig!
    # -----------------------------------------------------------------------
    @testset "calc_total_gw_uncon_irrig!" begin
        nc = 2
        nlevsoi = 3
        bounds_c = 1:nc
        mask_soilc = trues(nc)

        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, nc, 1, 1, 1)
        wf.qflx_gw_uncon_irrig_lyr_col .= 0.0
        wf.qflx_gw_uncon_irrig_col .= 0.0

        # Set layer values
        wf.qflx_gw_uncon_irrig_lyr_col[1, 1] = 1.0
        wf.qflx_gw_uncon_irrig_lyr_col[1, 2] = 2.0
        wf.qflx_gw_uncon_irrig_lyr_col[1, 3] = 3.0
        wf.qflx_gw_uncon_irrig_lyr_col[2, 1] = 0.5
        wf.qflx_gw_uncon_irrig_lyr_col[2, 2] = 1.5
        wf.qflx_gw_uncon_irrig_lyr_col[2, 3] = 0.0

        CLM.calc_total_gw_uncon_irrig!(wf, mask_soilc, bounds_c, nlevsoi)

        @test wf.qflx_gw_uncon_irrig_col[1] ≈ 6.0 atol=1e-10
        @test wf.qflx_gw_uncon_irrig_col[2] ≈ 2.0 atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test calc_application_fluxes! — drip/sprinkler routing
    # -----------------------------------------------------------------------
    @testset "calc_application_fluxes! drip/sprinkler" begin
        nc = 1
        np = 2
        nlevsoi = 2
        bounds_c = 1:nc
        bounds_p = 1:np

        irrig = CLM.IrrigationData()
        CLM.irrigation_init_allocate!(irrig, np, nc, nlevsoi)
        irrig.irrig_method_patch[1] = CLM.IRRIG_METHOD_DRIP
        irrig.irrig_method_patch[2] = CLM.IRRIG_METHOD_SPRINKLER

        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, nc, np, 1, 1)
        wf.qflx_sfc_irrig_col .= 0.0
        wf.qflx_gw_uncon_irrig_col .= 0.0
        wf.qflx_gw_con_irrig_col .= 0.0

        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        patch_data.column .= [1, 1]

        mask_soilc = trues(nc)
        mask_soilp = trues(np)

        qflx_sfc_irrig_bulk_patch = [5.0, 3.0]
        qflx_sfc_irrig_bulk_col = [8.0]
        qflx_gw_demand_bulk_patch = [0.0, 0.0]
        qflx_gw_demand_bulk_col = [0.0]

        CLM.calc_application_fluxes!(irrig, wf, true,
                                     qflx_sfc_irrig_bulk_patch,
                                     qflx_sfc_irrig_bulk_col,
                                     qflx_gw_demand_bulk_patch,
                                     qflx_gw_demand_bulk_col,
                                     patch_data, mask_soilc, mask_soilp,
                                     bounds_c, bounds_p)

        # Patch 1 is drip: total irrig = 5.0 + 0.0 = 5.0
        @test wf.qflx_irrig_drip_patch[1] ≈ 5.0 atol=1e-10
        @test wf.qflx_irrig_sprinkler_patch[1] ≈ 0.0 atol=1e-10

        # Patch 2 is sprinkler: total irrig = 3.0 + 0.0 = 3.0
        @test wf.qflx_irrig_drip_patch[2] ≈ 0.0 atol=1e-10
        @test wf.qflx_irrig_sprinkler_patch[2] ≈ 3.0 atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test calc_bulk_withdrawals! — basic bulk withdrawal
    # -----------------------------------------------------------------------
    @testset "calc_bulk_withdrawals! basic" begin
        nc = 1
        np = 2
        nlevsoi = 2
        bounds_c = 1:nc
        bounds_p = 1:np

        irrig = CLM.IrrigationData()
        CLM.irrigation_init_allocate!(irrig, np, nc, nlevsoi)
        irrig.params = CLM.IrrigationParamsData(use_groundwater_irrigation = false)

        # Set patch 1 to irrigate, patch 2 not irrigating
        irrig.n_irrig_steps_left_patch[1] = 3
        irrig.sfc_irrig_rate_patch[1] = 2.5
        irrig.irrig_rate_demand_patch[1] = 2.5
        irrig.n_irrig_steps_left_patch[2] = 0

        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, 1, 1)
        waterfluxbulk.wf.qflx_sfc_irrig_col .= 0.0

        soilhydrology = CLM.SoilHydrologyData()
        soilstate = CLM.SoilStateData()

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)

        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        patch_data.column .= [1, 1]
        patch_data.wtcol .= [0.6, 0.4]

        mask_soilc = trues(nc)
        mask_soilp = trues(np)

        qflx_sfc_bulk_p = zeros(np)
        qflx_gw_demand_bulk_p = zeros(np)
        qflx_gw_demand_bulk_c = zeros(nc)

        CLM.calc_bulk_withdrawals!(irrig, waterfluxbulk, soilhydrology, soilstate,
                                   col_data, patch_data, mask_soilc, mask_soilp,
                                   bounds_c, bounds_p, nlevsoi, 1800.0,
                                   qflx_sfc_bulk_p, qflx_gw_demand_bulk_p,
                                   qflx_gw_demand_bulk_c)

        # Patch 1 should have sfc_irrig = 2.5
        @test qflx_sfc_bulk_p[1] ≈ 2.5 atol=1e-10
        @test irrig.qflx_irrig_demand_patch[1] ≈ 2.5 atol=1e-10
        # GW demand = demand - sfc = 2.5 - 2.5 = 0
        @test qflx_gw_demand_bulk_p[1] ≈ 0.0 atol=1e-10

        # Patch 2 not irrigating
        @test qflx_sfc_bulk_p[2] ≈ 0.0 atol=1e-10
        @test irrig.qflx_irrig_demand_patch[2] ≈ 0.0 atol=1e-10

        # Steps left should be decremented for patch 1
        @test irrig.n_irrig_steps_left_patch[1] == 2
        @test irrig.n_irrig_steps_left_patch[2] == 0

        # Column sfc irrig = 2.5*0.6 + 0.0*0.4 = 1.5
        @test waterfluxbulk.wf.qflx_sfc_irrig_col[1] ≈ 1.5 atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test use_groundwater_irrigation accessor
    # -----------------------------------------------------------------------
    @testset "use_groundwater_irrigation accessor" begin
        irrig = CLM.IrrigationData()
        irrig.params.use_groundwater_irrigation = false
        @test CLM.use_groundwater_irrigation(irrig) == false

        irrig.params.use_groundwater_irrigation = true
        @test CLM.use_groundwater_irrigation(irrig) == true
    end

    # -----------------------------------------------------------------------
    # Test calc_deficit_volr_limited! — river volume limitation
    # -----------------------------------------------------------------------
    @testset "calc_deficit_volr_limited!" begin
        nc = 2
        ng = 1
        bounds_c = 1:nc
        bounds_g = 1:ng

        irrig = CLM.IrrigationData()
        irrig.params = CLM.IrrigationParamsData(irrig_river_volume_threshold = 0.1)

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)
        col_data.gridcell .= [1, 1]
        col_data.wtgcell .= [0.5, 0.5]

        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, ng)
        grc.area[1] = 100.0  # km^2

        # Deficit is 10 mm per column
        deficit = [10.0, 10.0]
        deficit_volr_limited = zeros(nc)
        check_for_irrig_col = trues(nc)

        # Large volr: no limitation
        volr = [1.0e9]  # very large
        CLM.calc_deficit_volr_limited!(irrig, deficit, volr, deficit_volr_limited,
                                       check_for_irrig_col, col_data, grc,
                                       bounds_c, bounds_g)
        @test deficit_volr_limited[1] ≈ 10.0 atol=1e-6
        @test deficit_volr_limited[2] ≈ 10.0 atol=1e-6

        # Zero volr: should limit to 0
        volr = [0.0]
        CLM.calc_deficit_volr_limited!(irrig, deficit, volr, deficit_volr_limited,
                                       check_for_irrig_col, col_data, grc,
                                       bounds_c, bounds_g)
        @test deficit_volr_limited[1] ≈ 0.0 atol=1e-10
        @test deficit_volr_limited[2] ≈ 0.0 atol=1e-10

        # Negative volr: should limit to 0
        volr = [-100.0]
        CLM.calc_deficit_volr_limited!(irrig, deficit, volr, deficit_volr_limited,
                                       check_for_irrig_col, col_data, grc,
                                       bounds_c, bounds_g)
        @test deficit_volr_limited[1] ≈ 0.0 atol=1e-10
        @test deficit_volr_limited[2] ≈ 0.0 atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Test irrigation_restart! stub
    # -----------------------------------------------------------------------
    @testset "irrigation_restart! stub" begin
        irrig = CLM.IrrigationData()
        @test CLM.irrigation_restart!(irrig) === nothing
    end

end
