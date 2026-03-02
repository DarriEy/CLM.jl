@testset "Canopy Hydrology" begin

    @testset "CanopyHydrologyParamsData defaults" begin
        p = CLM.CanopyHydrologyParamsData()
        @test p.liq_canopy_storage_scalar == 0.04
        @test p.snow_canopy_storage_scalar == 6.0
        @test p.snowcan_unload_temp_fact == 1.87e5
        @test p.snowcan_unload_wind_fact == 1.56e5
        @test p.interception_fraction == 1.0
        @test p.maximum_leaf_wetted_fraction == 1.0
    end

    @testset "CanopyHydrologyControl defaults" begin
        ctrl = CLM.CanopyHydrologyControl()
        @test ctrl.use_clm5_fpi == false
    end

    @testset "canopy_hydrology_read_nml!" begin
        CLM.canopy_hydrology_read_nml!(use_clm5_fpi=true)
        @test CLM.canopy_hydrology_ctrl.use_clm5_fpi == true
        CLM.canopy_hydrology_read_nml!(use_clm5_fpi=false)
        @test CLM.canopy_hydrology_ctrl.use_clm5_fpi == false
    end

    @testset "canopy_hydrology_read_params!" begin
        CLM.canopy_hydrology_read_params!(
            liq_canopy_storage_scalar=0.05,
            snow_canopy_storage_scalar=7.0,
            snowcan_unload_temp_fact=2.0e5,
            snowcan_unload_wind_fact=2.0e5,
            interception_fraction=0.8,
            maximum_leaf_wetted_fraction=0.9)
        @test CLM.canopy_hydrology_params.liq_canopy_storage_scalar == 0.05
        @test CLM.canopy_hydrology_params.snow_canopy_storage_scalar == 7.0
        @test CLM.canopy_hydrology_params.snowcan_unload_temp_fact == 2.0e5
        @test CLM.canopy_hydrology_params.snowcan_unload_wind_fact == 2.0e5
        @test CLM.canopy_hydrology_params.interception_fraction == 0.8
        @test CLM.canopy_hydrology_params.maximum_leaf_wetted_fraction == 0.9
        # Reset to defaults
        CLM.canopy_hydrology_read_params!()
    end

    @testset "sum_flux_top_of_canopy_inputs!" begin
        np = 2
        nc = 1

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.column[1] = 1
        patch.column[2] = 1

        mask_nolakep = BitVector([true, true])
        bounds_p = 1:np

        forc_rain = [3.0e-4]       # column-level rain
        forc_snow_col = [1.0e-4]   # column-level snow
        qflx_irrig_sprinkler = [0.5e-4, 1.0e-4]  # patch-level irrigation
        qflx_liq_above_canopy = zeros(np)
        forc_snow_patch = zeros(np)

        CLM.sum_flux_top_of_canopy_inputs!(
            patch, mask_nolakep, bounds_p,
            forc_rain, forc_snow_col, qflx_irrig_sprinkler,
            qflx_liq_above_canopy, forc_snow_patch)

        @test qflx_liq_above_canopy[1] ≈ 3.0e-4 + 0.5e-4
        @test qflx_liq_above_canopy[2] ≈ 3.0e-4 + 1.0e-4
        @test forc_snow_patch[1] ≈ 1.0e-4
        @test forc_snow_patch[2] ≈ 1.0e-4
    end

    @testset "bulk_flux_canopy_interception_and_throughfall! — vegetated" begin
        np = 1
        nc = 1

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.column[1] = 1

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.itype[1] = 1  # soil

        mask_nolakep = BitVector([true])
        bounds_p = 1:np

        frac_veg_nosno = [1]
        elai = [2.0]
        esai = [0.5]
        forc_snow = [2.0e-4]
        qflx_liq_above_canopy = [5.0e-4]

        qflx_through_snow = zeros(np)
        qflx_through_liq = zeros(np)
        qflx_intercepted_snow = zeros(np)
        qflx_intercepted_liq = zeros(np)
        check_point = fill(false, np)

        # Default fpi (not clm5)
        CLM.canopy_hydrology_ctrl.use_clm5_fpi = false

        CLM.bulk_flux_canopy_interception_and_throughfall!(
            patch, col, mask_nolakep, bounds_p,
            frac_veg_nosno, elai, esai, forc_snow, qflx_liq_above_canopy,
            qflx_through_snow, qflx_through_liq,
            qflx_intercepted_snow, qflx_intercepted_liq,
            check_point)

        @test check_point[1] == true

        # Verify interception coefficients
        lai_sai = elai[1] + esai[1]
        fpiliq = 0.25 * (1.0 - exp(-0.5 * lai_sai))
        fpisnow = 1.0 - exp(-0.5 * lai_sai)

        @test qflx_through_snow[1] ≈ forc_snow[1] * (1.0 - fpisnow)
        @test qflx_through_liq[1] ≈ qflx_liq_above_canopy[1] * (1.0 - fpiliq)
        @test qflx_intercepted_snow[1] ≈ forc_snow[1] * fpisnow
        @test qflx_intercepted_liq[1] ≈ qflx_liq_above_canopy[1] * fpiliq

        # Conservation: throughfall + interception = total input
        @test qflx_through_snow[1] + qflx_intercepted_snow[1] ≈ forc_snow[1]
        @test qflx_through_liq[1] + qflx_intercepted_liq[1] ≈ qflx_liq_above_canopy[1]
    end

    @testset "bulk_flux_canopy_interception_and_throughfall! — clm5 fpi" begin
        np = 1
        nc = 1

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.column[1] = 1

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.itype[1] = 1

        mask_nolakep = BitVector([true])
        bounds_p = 1:np

        frac_veg_nosno = [1]
        elai = [2.0]
        esai = [0.5]
        forc_snow = [2.0e-4]
        qflx_liq_above_canopy = [5.0e-4]

        qflx_through_snow = zeros(np)
        qflx_through_liq = zeros(np)
        qflx_intercepted_snow = zeros(np)
        qflx_intercepted_liq = zeros(np)
        check_point = fill(false, np)

        CLM.canopy_hydrology_ctrl.use_clm5_fpi = true
        CLM.canopy_hydrology_params.interception_fraction = 0.8

        CLM.bulk_flux_canopy_interception_and_throughfall!(
            patch, col, mask_nolakep, bounds_p,
            frac_veg_nosno, elai, esai, forc_snow, qflx_liq_above_canopy,
            qflx_through_snow, qflx_through_liq,
            qflx_intercepted_snow, qflx_intercepted_liq,
            check_point)

        lai_sai = elai[1] + esai[1]
        fpiliq_clm5 = 0.8 * tanh(lai_sai)
        @test qflx_intercepted_liq[1] ≈ qflx_liq_above_canopy[1] * fpiliq_clm5
        @test qflx_through_liq[1] ≈ qflx_liq_above_canopy[1] * (1.0 - fpiliq_clm5)

        # Reset
        CLM.canopy_hydrology_ctrl.use_clm5_fpi = false
        CLM.canopy_hydrology_read_params!()
    end

    @testset "bulk_flux_canopy_interception_and_throughfall! — no veg" begin
        np = 1
        nc = 1

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.column[1] = 1

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.itype[1] = 1  # soil

        mask_nolakep = BitVector([true])
        bounds_p = 1:np

        frac_veg_nosno = [0]
        elai = [0.0]
        esai = [0.0]
        forc_snow = [2.0e-4]
        qflx_liq_above_canopy = [5.0e-4]

        qflx_through_snow = zeros(np)
        qflx_through_liq = zeros(np)
        qflx_intercepted_snow = zeros(np)
        qflx_intercepted_liq = zeros(np)
        check_point = fill(false, np)

        CLM.bulk_flux_canopy_interception_and_throughfall!(
            patch, col, mask_nolakep, bounds_p,
            frac_veg_nosno, elai, esai, forc_snow, qflx_liq_above_canopy,
            qflx_through_snow, qflx_through_liq,
            qflx_intercepted_snow, qflx_intercepted_liq,
            check_point)

        @test check_point[1] == false
        @test qflx_intercepted_snow[1] == 0.0
        @test qflx_intercepted_liq[1] == 0.0
        @test qflx_through_snow[1] ≈ forc_snow[1]
        @test qflx_through_liq[1] ≈ qflx_liq_above_canopy[1]
    end

    @testset "bulk_flux_canopy_interception_and_throughfall! — sunwall" begin
        np = 1
        nc = 1

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.column[1] = 1

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.itype[1] = CLM.ICOL_SUNWALL

        mask_nolakep = BitVector([true])
        bounds_p = 1:np

        frac_veg_nosno = [0]
        elai = [0.0]
        esai = [0.0]
        forc_snow = [2.0e-4]
        qflx_liq_above_canopy = [5.0e-4]

        qflx_through_snow = zeros(np)
        qflx_through_liq = zeros(np)
        qflx_intercepted_snow = zeros(np)
        qflx_intercepted_liq = zeros(np)
        check_point = fill(false, np)

        CLM.bulk_flux_canopy_interception_and_throughfall!(
            patch, col, mask_nolakep, bounds_p,
            frac_veg_nosno, elai, esai, forc_snow, qflx_liq_above_canopy,
            qflx_through_snow, qflx_through_liq,
            qflx_intercepted_snow, qflx_intercepted_liq,
            check_point)

        @test qflx_through_snow[1] == 0.0
        @test qflx_through_liq[1] == 0.0
    end

    @testset "update_state_add_interception_to_canopy!" begin
        np = 1
        mask_soilp = BitVector([true])
        bounds_p = 1:np
        dtime = 1800.0

        qflx_intercepted_snow = [1.0e-4]
        qflx_intercepted_liq = [2.0e-4]
        snocan = [0.5]
        liqcan = [0.3]

        CLM.update_state_add_interception_to_canopy!(
            mask_soilp, bounds_p, dtime,
            qflx_intercepted_snow, qflx_intercepted_liq,
            snocan, liqcan)

        @test snocan[1] ≈ max(0.0, 0.5 + dtime * 1.0e-4)
        @test liqcan[1] ≈ max(0.0, 0.3 + dtime * 2.0e-4)
    end

    @testset "bulk_flux_canopy_excess!" begin
        np = 1
        mask_soilp = BitVector([true])
        bounds_p = 1:np
        dtime = 1800.0

        elai = [2.0]
        esai = [0.5]

        # Set canopy water above max
        liq_max = CLM.canopy_hydrology_params.liq_canopy_storage_scalar * (elai[1] + esai[1])
        sno_max = CLM.canopy_hydrology_params.snow_canopy_storage_scalar * (elai[1] + esai[1])

        snocan = [sno_max + 2.0]  # excess snow
        liqcan = [liq_max + 0.5]  # excess liquid

        check_point = [true]

        qflx_snocanfall = zeros(np)
        qflx_liqcanfall = zeros(np)

        CLM.bulk_flux_canopy_excess!(
            mask_soilp, bounds_p, dtime,
            elai, esai, snocan, liqcan, check_point,
            qflx_snocanfall, qflx_liqcanfall)

        @test qflx_liqcanfall[1] ≈ max((liqcan[1] - liq_max) / dtime, 0.0)
        @test qflx_snocanfall[1] ≈ max((snocan[1] - sno_max) / dtime, 0.0)
        @test qflx_liqcanfall[1] > 0.0
        @test qflx_snocanfall[1] > 0.0
    end

    @testset "bulk_flux_canopy_excess! — no excess" begin
        np = 1
        mask_soilp = BitVector([true])
        bounds_p = 1:np
        dtime = 1800.0

        elai = [2.0]
        esai = [0.5]
        snocan = [0.0]
        liqcan = [0.0]
        check_point = [true]

        qflx_snocanfall = zeros(np)
        qflx_liqcanfall = zeros(np)

        CLM.bulk_flux_canopy_excess!(
            mask_soilp, bounds_p, dtime,
            elai, esai, snocan, liqcan, check_point,
            qflx_snocanfall, qflx_liqcanfall)

        @test qflx_liqcanfall[1] == 0.0
        @test qflx_snocanfall[1] == 0.0
    end

    @testset "update_state_remove_canfall_from_canopy!" begin
        np = 1
        mask_soilp = BitVector([true])
        bounds_p = 1:np
        dtime = 1800.0

        qflx_liqcanfall = [1.0e-5]
        qflx_snocanfall = [2.0e-5]
        liqcan = [0.5]
        snocan = [1.0]

        liqcan_before = liqcan[1]
        snocan_before = snocan[1]

        CLM.update_state_remove_canfall_from_canopy!(
            mask_soilp, bounds_p, dtime,
            qflx_liqcanfall, qflx_snocanfall,
            liqcan, snocan)

        @test liqcan[1] ≈ liqcan_before - dtime * qflx_liqcanfall[1]
        @test snocan[1] ≈ snocan_before - dtime * qflx_snocanfall[1]
    end

    @testset "bulk_flux_snow_unloading! — warm with snow" begin
        np = 1
        nc = 1
        ng = 1

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.column[1] = 1
        patch.gridcell[1] = 1

        mask_soilp = BitVector([true])
        bounds_p = 1:np
        dtime = 1800.0

        frac_veg_nosno = [1]
        forc_t = [275.0]     # above 270.15 K threshold
        forc_wind = [5.0]
        snocan = [3.0]

        qflx_snotempunload = zeros(np)
        qflx_snowindunload = zeros(np)
        qflx_snow_unload = zeros(np)

        CLM.bulk_flux_snow_unloading!(
            patch, mask_soilp, bounds_p, dtime,
            frac_veg_nosno, forc_t, forc_wind, snocan,
            qflx_snotempunload, qflx_snowindunload, qflx_snow_unload)

        @test qflx_snotempunload[1] > 0.0
        @test qflx_snowindunload[1] > 0.0
        @test qflx_snow_unload[1] > 0.0
        @test qflx_snow_unload[1] <= snocan[1] / dtime
    end

    @testset "bulk_flux_snow_unloading! — no veg" begin
        np = 1

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.column[1] = 1
        patch.gridcell[1] = 1

        mask_soilp = BitVector([true])
        bounds_p = 1:np
        dtime = 1800.0

        frac_veg_nosno = [0]
        forc_t = [275.0]
        forc_wind = [5.0]
        snocan = [3.0]

        qflx_snotempunload = zeros(np)
        qflx_snowindunload = zeros(np)
        qflx_snow_unload = zeros(np)

        CLM.bulk_flux_snow_unloading!(
            patch, mask_soilp, bounds_p, dtime,
            frac_veg_nosno, forc_t, forc_wind, snocan,
            qflx_snotempunload, qflx_snowindunload, qflx_snow_unload)

        @test qflx_snotempunload[1] == 0.0
        @test qflx_snowindunload[1] == 0.0
        @test qflx_snow_unload[1] == 0.0
    end

    @testset "update_state_remove_snow_unloading! — partial" begin
        np = 1
        mask_soilp = BitVector([true])
        bounds_p = 1:np
        dtime = 1800.0

        snocan = [3.0]
        qflx_snow_unload = [0.001]  # will not empty canopy

        CLM.update_state_remove_snow_unloading!(
            mask_soilp, bounds_p, dtime,
            qflx_snow_unload, snocan)

        @test snocan[1] ≈ 3.0 - 0.001 * dtime
    end

    @testset "update_state_remove_snow_unloading! — exact emptying" begin
        np = 1
        mask_soilp = BitVector([true])
        bounds_p = 1:np
        dtime = 1800.0

        snocan_val = 3.0
        snocan = [snocan_val]
        qflx_snow_unload = [snocan_val / dtime]  # exact emptying

        CLM.update_state_remove_snow_unloading!(
            mask_soilp, bounds_p, dtime,
            qflx_snow_unload, snocan)

        @test snocan[1] == 0.0
    end

    @testset "sum_flux_fluxes_onto_ground!" begin
        np = 2
        nc = 1

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.column[1] = 1
        patch.column[2] = 1
        patch.wtcol[1] = 0.6
        patch.wtcol[2] = 0.4

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        mask_nolakep = BitVector([true, true])
        mask_nolakec = BitVector([true])
        bounds_p = 1:np
        bounds_c = 1:nc

        qflx_through_snow = [1.0e-4, 2.0e-4]
        qflx_snocanfall = [0.5e-4, 0.5e-4]
        qflx_snow_unload = [0.2e-4, 0.3e-4]
        qflx_through_liq = [3.0e-4, 4.0e-4]
        qflx_liqcanfall = [0.3e-4, 0.2e-4]
        qflx_irrig_drip = [0.1e-4, 0.1e-4]

        qflx_snow_grnd_col = zeros(nc)
        qflx_liq_grnd_col = zeros(nc)
        qflx_snow_h2osfc = zeros(nc)

        CLM.sum_flux_fluxes_onto_ground!(
            patch, col, mask_nolakep, mask_nolakec,
            bounds_p, bounds_c,
            qflx_through_snow, qflx_snocanfall, qflx_snow_unload,
            qflx_through_liq, qflx_liqcanfall, qflx_irrig_drip,
            qflx_snow_grnd_col, qflx_liq_grnd_col, qflx_snow_h2osfc)

        # Expected: weighted average of patch sums
        snow_p1 = qflx_through_snow[1] + qflx_snocanfall[1] + qflx_snow_unload[1]
        snow_p2 = qflx_through_snow[2] + qflx_snocanfall[2] + qflx_snow_unload[2]
        liq_p1 = qflx_through_liq[1] + qflx_liqcanfall[1] + qflx_irrig_drip[1]
        liq_p2 = qflx_through_liq[2] + qflx_liqcanfall[2] + qflx_irrig_drip[2]

        @test qflx_snow_grnd_col[1] ≈ snow_p1 * 0.6 + snow_p2 * 0.4
        @test qflx_liq_grnd_col[1] ≈ liq_p1 * 0.6 + liq_p2 * 0.4
        @test qflx_snow_h2osfc[1] == 0.0
    end

    @testset "bulk_diag_frac_wet! — wet canopy" begin
        np = 1
        mask_soilp = BitVector([true])
        bounds_p = 1:np

        frac_veg_nosno = [1]
        elai = [2.0]
        esai = [0.5]
        snocan = [0.5]
        liqcan = [0.3]

        fwet = zeros(np)
        fdry = zeros(np)
        fcansno = zeros(np)

        CLM.bulk_diag_frac_wet!(
            mask_soilp, bounds_p,
            frac_veg_nosno, elai, esai, snocan, liqcan,
            fwet, fdry, fcansno)

        h2ocan = snocan[1] + liqcan[1]
        vegt = 1.0 * (elai[1] + esai[1])
        expected_fwet = (h2ocan / (vegt * CLM.canopy_hydrology_params.liq_canopy_storage_scalar))^0.666666666666
        expected_fwet = min(expected_fwet, CLM.canopy_hydrology_params.maximum_leaf_wetted_fraction)
        expected_fcansno = (snocan[1] / (vegt * CLM.canopy_hydrology_params.snow_canopy_storage_scalar))^0.15
        expected_fcansno = min(expected_fcansno, 1.0)
        expected_fdry = (1.0 - expected_fwet) * elai[1] / (elai[1] + esai[1])

        @test fwet[1] ≈ expected_fwet
        @test fdry[1] ≈ expected_fdry
        @test fcansno[1] ≈ expected_fcansno
    end

    @testset "bulk_diag_frac_wet! — dry canopy" begin
        np = 1
        mask_soilp = BitVector([true])
        bounds_p = 1:np

        frac_veg_nosno = [1]
        elai = [2.0]
        esai = [0.5]
        snocan = [0.0]
        liqcan = [0.0]

        fwet = ones(np)
        fdry = ones(np)
        fcansno = ones(np)

        CLM.bulk_diag_frac_wet!(
            mask_soilp, bounds_p,
            frac_veg_nosno, elai, esai, snocan, liqcan,
            fwet, fdry, fcansno)

        @test fwet[1] == 0.0
        @test fcansno[1] == 0.0
        @test fdry[1] ≈ elai[1] / (elai[1] + esai[1])
    end

    @testset "bulk_diag_frac_wet! — no vegetation" begin
        np = 1
        mask_soilp = BitVector([true])
        bounds_p = 1:np

        frac_veg_nosno = [0]
        elai = [0.0]
        esai = [0.0]
        snocan = [0.0]
        liqcan = [0.0]

        fwet = ones(np)
        fdry = ones(np)
        fcansno = ones(np)

        CLM.bulk_diag_frac_wet!(
            mask_soilp, bounds_p,
            frac_veg_nosno, elai, esai, snocan, liqcan,
            fwet, fdry, fcansno)

        @test fwet[1] == 0.0
        @test fdry[1] == 0.0
    end

    @testset "tracer_flux_canopy_interception_and_throughfall!" begin
        np = 1
        mask_nolakep = BitVector([true])
        bounds_p = 1:np

        # Bulk values
        bulk_forc_snow = [2.0e-4]
        bulk_qflx_liq_above_canopy = [5.0e-4]
        bulk_qflx_through_snow = [1.0e-4]
        bulk_qflx_intercepted_snow = [1.0e-4]
        bulk_qflx_through_liq = [3.0e-4]
        bulk_qflx_intercepted_liq = [2.0e-4]

        # Tracer values (proportional)
        ratio_snow = 0.9
        ratio_liq = 0.8
        trac_forc_snow = [bulk_forc_snow[1] * ratio_snow]
        trac_qflx_liq_above_canopy = [bulk_qflx_liq_above_canopy[1] * ratio_liq]

        trac_qflx_through_snow = zeros(np)
        trac_qflx_intercepted_snow = zeros(np)
        trac_qflx_through_liq = zeros(np)
        trac_qflx_intercepted_liq = zeros(np)

        CLM.tracer_flux_canopy_interception_and_throughfall!(
            mask_nolakep, bounds_p,
            bulk_forc_snow, bulk_qflx_liq_above_canopy,
            bulk_qflx_through_snow, bulk_qflx_intercepted_snow,
            bulk_qflx_through_liq, bulk_qflx_intercepted_liq,
            trac_forc_snow, trac_qflx_liq_above_canopy,
            trac_qflx_through_snow, trac_qflx_intercepted_snow,
            trac_qflx_through_liq, trac_qflx_intercepted_liq)

        @test trac_qflx_through_snow[1] ≈ ratio_snow * bulk_qflx_through_snow[1]
        @test trac_qflx_intercepted_snow[1] ≈ ratio_snow * bulk_qflx_intercepted_snow[1]
        @test trac_qflx_through_liq[1] ≈ ratio_liq * bulk_qflx_through_liq[1]
        @test trac_qflx_intercepted_liq[1] ≈ ratio_liq * bulk_qflx_intercepted_liq[1]
    end

    @testset "tracer_flux_canopy_excess!" begin
        np = 1
        mask_soilp = BitVector([true])
        bounds_p = 1:np

        bulk_liqcan = [0.5]
        bulk_snocan = [1.0]
        bulk_qflx_liqcanfall = [0.001]
        bulk_qflx_snocanfall = [0.002]

        ratio_liq = 0.9
        ratio_sno = 0.8
        trac_liqcan = [bulk_liqcan[1] * ratio_liq]
        trac_snocan = [bulk_snocan[1] * ratio_sno]

        trac_qflx_liqcanfall = zeros(np)
        trac_qflx_snocanfall = zeros(np)

        CLM.tracer_flux_canopy_excess!(
            mask_soilp, bounds_p,
            bulk_liqcan, bulk_snocan,
            bulk_qflx_liqcanfall, bulk_qflx_snocanfall,
            trac_liqcan, trac_snocan,
            trac_qflx_liqcanfall, trac_qflx_snocanfall)

        @test trac_qflx_liqcanfall[1] ≈ ratio_liq * bulk_qflx_liqcanfall[1]
        @test trac_qflx_snocanfall[1] ≈ ratio_sno * bulk_qflx_snocanfall[1]
    end

    @testset "tracer_flux_snow_unloading!" begin
        np = 1
        mask_soilp = BitVector([true])
        bounds_p = 1:np

        bulk_snocan = [3.0]
        bulk_qflx_snow_unload = [0.001]
        ratio = 0.85
        trac_snocan = [bulk_snocan[1] * ratio]

        trac_qflx_snow_unload = zeros(np)

        CLM.tracer_flux_snow_unloading!(
            mask_soilp, bounds_p,
            bulk_snocan, bulk_qflx_snow_unload,
            trac_snocan, trac_qflx_snow_unload)

        @test trac_qflx_snow_unload[1] ≈ ratio * bulk_qflx_snow_unload[1]
    end

    @testset "Water conservation — full pipeline" begin
        # Test that total water is conserved through the interception pipeline
        np = 1
        nc = 1
        dtime = 1800.0

        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        patch.column[1] = 1
        patch.gridcell[1] = 1
        patch.wtcol[1] = 1.0

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.itype[1] = 1

        mask_soilp = BitVector([true])
        mask_nolakep = BitVector([true])
        mask_nolakec = BitVector([true])
        bounds_p = 1:np
        bounds_c = 1:nc

        # Input forcing
        forc_rain = [5.0e-4]
        forc_snow_col = [2.0e-4]
        qflx_irrig_sprinkler = [1.0e-4]
        qflx_irrig_drip = [0.5e-4]

        total_input_liq = forc_rain[1] + qflx_irrig_sprinkler[1] + qflx_irrig_drip[1]
        total_input_snow = forc_snow_col[1]

        # Step 1: Top-of-canopy
        qflx_liq_above = zeros(np)
        forc_snow_p = zeros(np)
        CLM.sum_flux_top_of_canopy_inputs!(
            patch, mask_nolakep, bounds_p,
            forc_rain, forc_snow_col, qflx_irrig_sprinkler,
            qflx_liq_above, forc_snow_p)

        # Step 2: Interception
        frac_veg_nosno = [1]
        elai = [2.0]
        esai = [0.5]

        qflx_through_snow = zeros(np)
        qflx_through_liq = zeros(np)
        qflx_intercepted_snow = zeros(np)
        qflx_intercepted_liq = zeros(np)
        check_point = fill(false, np)

        CLM.bulk_flux_canopy_interception_and_throughfall!(
            patch, col, mask_nolakep, bounds_p,
            frac_veg_nosno, elai, esai, forc_snow_p, qflx_liq_above,
            qflx_through_snow, qflx_through_liq,
            qflx_intercepted_snow, qflx_intercepted_liq,
            check_point)

        # Verify conservation at interception step
        @test qflx_through_snow[1] + qflx_intercepted_snow[1] ≈ forc_snow_p[1]
        @test qflx_through_liq[1] + qflx_intercepted_liq[1] ≈ qflx_liq_above[1]

        # Step 3: Add interception
        snocan = [0.0]
        liqcan = [0.0]
        CLM.update_state_add_interception_to_canopy!(
            mask_soilp, bounds_p, dtime,
            qflx_intercepted_snow, qflx_intercepted_liq,
            snocan, liqcan)

        # Step 4: Canopy excess
        qflx_snocanfall = zeros(np)
        qflx_liqcanfall = zeros(np)
        CLM.bulk_flux_canopy_excess!(
            mask_soilp, bounds_p, dtime,
            elai, esai, snocan, liqcan, check_point,
            qflx_snocanfall, qflx_liqcanfall)

        # Step 5: Remove canfall
        CLM.update_state_remove_canfall_from_canopy!(
            mask_soilp, bounds_p, dtime,
            qflx_liqcanfall, qflx_snocanfall,
            liqcan, snocan)

        # Step 6: Snow unloading
        forc_t = [275.0]
        forc_wind = [5.0]
        qflx_snotempunload = zeros(np)
        qflx_snowindunload = zeros(np)
        qflx_snow_unload = zeros(np)

        CLM.bulk_flux_snow_unloading!(
            patch, mask_soilp, bounds_p, dtime,
            frac_veg_nosno, forc_t, forc_wind, snocan,
            qflx_snotempunload, qflx_snowindunload, qflx_snow_unload)

        # Step 7: Remove snow unloading
        CLM.update_state_remove_snow_unloading!(
            mask_soilp, bounds_p, dtime,
            qflx_snow_unload, snocan)

        # Step 8: Sum onto ground
        qflx_snow_grnd_col = zeros(nc)
        qflx_liq_grnd_col = zeros(nc)
        qflx_snow_h2osfc = zeros(nc)

        CLM.sum_flux_fluxes_onto_ground!(
            patch, col, mask_nolakep, mask_nolakec,
            bounds_p, bounds_c,
            qflx_through_snow, qflx_snocanfall, qflx_snow_unload,
            qflx_through_liq, qflx_liqcanfall, qflx_irrig_drip,
            qflx_snow_grnd_col, qflx_liq_grnd_col, qflx_snow_h2osfc)

        # Conservation check: ground flux + canopy storage change = total input
        # canopy_storage_change = (snocan + liqcan) / dtime  (net rate of storage)
        # total_to_ground = qflx_snow_grnd_col + qflx_liq_grnd_col
        # total_input = forc_rain + forc_snow + qflx_irrig_sprinkler + qflx_irrig_drip
        total_to_ground = qflx_snow_grnd_col[1] + qflx_liq_grnd_col[1]
        canopy_storage_rate = (snocan[1] + liqcan[1]) / dtime

        @test total_to_ground + canopy_storage_rate ≈ total_input_liq + total_input_snow atol=1e-12
    end

end
