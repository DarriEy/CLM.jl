@testset "Lake Hydrology Module" begin
    # ------------------------------------------------------------------
    # Tests for LakeHydrologyMod port.
    # Verifies:
    #   1. sum_flux_fluxes_onto_ground!
    #   2. lake_sublimation_dew!
    #   3. lake_frac_sno_eff!
    #   4. lake_patch_to_col_fluxes!
    #   5. lake_soil_hydrology!
    #   6. lake_check_single_snow_layer!
    #   7. lake_snow_above_unfrozen!
    #   8. lake_snow_diagnostics!
    #   9. lake_water_balance!
    #   10. lake_update_h2osoi_vol!
    #   11. lake_top_layer_diagnostics!
    #   12. lake_hydrology! (orchestrator smoke test)
    # ------------------------------------------------------------------

    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsno  = CLM.varpar.nlevsno
    nlevsoi  = CLM.varpar.nlevsoi
    nlevgrnd = CLM.varpar.nlevgrnd
    nlevlak  = CLM.varpar.nlevlak
    nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd

    # ------------------------------------------------------------------
    # 1. sum_flux_fluxes_onto_ground!
    # ------------------------------------------------------------------
    @testset "sum_flux_fluxes_onto_ground!" begin
        nc = 3
        bounds = 1:nc
        forc_snow = [0.5, 1.0, 0.0]
        forc_rain = [2.0, 0.0, 1.5]
        qflx_snow_grnd = fill(NaN, nc)
        qflx_liq_grnd  = fill(NaN, nc)
        mask_lake = BitVector([true, true, false])

        CLM.sum_flux_fluxes_onto_ground!(
            forc_snow, forc_rain,
            qflx_snow_grnd, qflx_liq_grnd,
            mask_lake, bounds)

        @test qflx_snow_grnd[1] ≈ 0.5
        @test qflx_snow_grnd[2] ≈ 1.0
        @test isnan(qflx_snow_grnd[3])  # not a lake column

        @test qflx_liq_grnd[1] ≈ 2.0
        @test qflx_liq_grnd[2] ≈ 0.0
        @test isnan(qflx_liq_grnd[3])   # not a lake column
    end

    # ------------------------------------------------------------------
    # 2. lake_sublimation_dew!
    # ------------------------------------------------------------------
    @testset "lake_sublimation_dew!" begin
        # Test case: no snow layers, evaporation
        nc = 2
        np = 2
        nlevtot = nlevsno + nlevmaxurbgrnd
        dtime = 1800.0
        bounds_patch = 1:np

        qflx_liqevap   = fill(0.0, np)
        qflx_solidevap = fill(0.0, np)
        qflx_soliddew  = fill(0.0, np)
        qflx_liqdew    = fill(0.0, np)
        qflx_ev_snow   = fill(0.0, np)
        qflx_evap_soi  = [0.001, -0.002]  # p1: evap, p2: dew

        h2osoi_liq  = fill(10.0, nc, nlevtot)
        h2osoi_ice  = fill(5.0, nc, nlevtot)
        h2osno_no_layers = [2.0, 0.0]
        snow_depth  = [0.01, 0.0]
        frac_sno    = [0.5, 0.0]
        t_grnd      = [260.0, 260.0]   # below freezing
        t_soisno    = fill(260.0, nc, nlevtot)
        snl         = [0, 0]  # no snow layers
        patch_column = [1, 2]
        mask_lakep  = BitVector([true, true])

        CLM.lake_sublimation_dew!(
            qflx_liqevap, qflx_solidevap,
            qflx_soliddew, qflx_liqdew,
            qflx_ev_snow,
            qflx_evap_soi, h2osoi_liq, h2osoi_ice,
            h2osno_no_layers, snow_depth, frac_sno,
            t_grnd, t_soisno, snl,
            patch_column, mask_lakep, bounds_patch,
            dtime, nlevsno)

        # Patch 1: evaporation, no snow layers
        # solidevap = min(0.001, 2.0/1800) = 0.001 (since 2.0/1800 ≈ 0.00111 > 0.001)
        @test qflx_solidevap[1] ≈ 0.001
        @test qflx_liqevap[1] ≈ 0.0

        # Patch 2: dew, t_grnd < tfrz-0.1 → frost
        @test qflx_soliddew[2] ≈ 0.002
        @test qflx_liqdew[2] ≈ 0.0
        # Snow pack should have increased due to frost
        @test h2osno_no_layers[2] ≈ 0.002 * dtime
    end

    # ------------------------------------------------------------------
    # 3. lake_frac_sno_eff!
    # ------------------------------------------------------------------
    @testset "lake_frac_sno_eff!" begin
        nc = 3
        bounds = 1:nc
        frac_sno = [0.5, 0.0, 0.8]
        frac_sno_eff = fill(NaN, nc)
        mask_lake = BitVector([true, false, true])

        CLM.lake_frac_sno_eff!(frac_sno, frac_sno_eff, mask_lake, bounds)

        @test frac_sno_eff[1] ≈ 0.5
        @test isnan(frac_sno_eff[2])  # not a lake column
        @test frac_sno_eff[3] ≈ 0.8
    end

    # ------------------------------------------------------------------
    # 4. lake_patch_to_col_fluxes!
    # ------------------------------------------------------------------
    @testset "lake_patch_to_col_fluxes!" begin
        nc = 2
        np = 2
        bounds_patch = 1:np

        # Column-level (output)
        qflx_evap_tot_col      = fill(NaN, nc)
        qflx_liqevap_col       = fill(NaN, nc)
        qflx_liqdew_col        = fill(NaN, nc)
        qflx_soliddew_col      = fill(NaN, nc)
        qflx_solidevap_col     = fill(NaN, nc)
        qflx_ev_snow_col       = fill(NaN, nc)

        # Patch-level (input)
        qflx_evap_tot = [1.0, 2.0]
        qflx_liqevap  = [0.1, 0.2]
        qflx_liqdew   = [0.3, 0.4]
        qflx_soliddew = [0.5, 0.6]
        qflx_solidevap = [0.7, 0.8]
        qflx_ev_snow  = [0.9, 1.1]

        patch_column = [1, 2]
        mask_lakep = BitVector([true, true])

        CLM.lake_patch_to_col_fluxes!(
            qflx_evap_tot_col, qflx_liqevap_col,
            qflx_liqdew_col, qflx_soliddew_col,
            qflx_solidevap_col, qflx_ev_snow_col,
            qflx_evap_tot, qflx_liqevap,
            qflx_liqdew, qflx_soliddew,
            qflx_solidevap, qflx_ev_snow,
            patch_column, mask_lakep, bounds_patch)

        @test qflx_evap_tot_col[1] ≈ 1.0
        @test qflx_evap_tot_col[2] ≈ 2.0
        @test qflx_ev_snow_col[2] ≈ 1.1
    end

    # ------------------------------------------------------------------
    # 5. lake_soil_hydrology!
    # ------------------------------------------------------------------
    @testset "lake_soil_hydrology!" begin
        nc = 2
        bounds = 1:nc
        nlevtot = nlevsno + nlevmaxurbgrnd

        h2osoi_liq = fill(5.0, nc, nlevtot)
        h2osoi_ice = fill(1.0, nc, nlevtot)
        h2osoi_vol = fill(0.0, nc, nlevmaxurbgrnd)
        dz = fill(0.1, nc, nlevtot)
        watsat = fill(0.4, nc, nlevmaxurbgrnd)
        mask_lake = BitVector([true, true])

        # Make soil layer 1 unsaturated (low liq/ice)
        jj = 1 + nlevsno
        h2osoi_liq[1, jj] = 1.0   # small amount
        h2osoi_ice[1, jj] = 0.5
        # vol = 1.0/(0.1*1000) + 0.5/(0.1*917) ≈ 0.01 + 0.00545 = 0.01545 < 0.4

        CLM.lake_soil_hydrology!(
            h2osoi_liq, h2osoi_ice, h2osoi_vol,
            dz, watsat, mask_lake, bounds,
            nlevsoi, nlevsno)

        # Layer 1 was undersaturated, so should be filled to saturation
        expected_liq = (watsat[1, 1] * dz[1, jj] - h2osoi_ice[1, jj] / CLM.DENICE) * CLM.DENH2O
        @test h2osoi_liq[1, jj] ≈ expected_liq
    end

    # ------------------------------------------------------------------
    # 6. lake_check_single_snow_layer!
    # ------------------------------------------------------------------
    @testset "lake_check_single_snow_layer!" begin
        nc = 2
        np = 2
        nlevtot = nlevsno + nlevmaxurbgrnd
        dtime = 1800.0

        eflx_sh_tot    = fill(0.0, np)
        eflx_sh_grnd   = fill(0.0, np)
        eflx_soil_grnd = fill(0.0, np)
        eflx_gnet      = fill(10.0, np)
        eflx_grnd_lake = fill(0.0, np)

        t_soisno     = fill(CLM.TFRZ - 1.0, nc, nlevtot)
        h2osoi_ice   = fill(0.0, nc, nlevtot)
        h2osoi_liq   = fill(1.0, nc, nlevtot)
        snl          = [-1, 0]   # col 1 has one snow layer, col 2 has none
        h2osno_no_layers = [0.0, 0.0]
        h2osno_total = [2.0, 0.0]
        snow_depth   = [0.01, 0.0]
        qflx_snow_drain = fill(0.0, nc)
        patch_column = [1, 2]
        mask_lakep = BitVector([true, true])
        bounds_patch = 1:np

        # Set layer 0 (Julia index nlevsno) to have no ice → should remove layer
        j = nlevsno  # snow layer 0 Julia index
        h2osoi_ice[1, j] = 0.0
        h2osoi_liq[1, j] = 0.5
        t_soisno[1, j] = CLM.TFRZ + 1.0  # above freezing

        CLM.lake_check_single_snow_layer!(
            eflx_sh_tot, eflx_sh_grnd, eflx_soil_grnd,
            eflx_gnet, eflx_grnd_lake,
            t_soisno, h2osoi_ice, h2osoi_liq,
            snl, h2osno_no_layers, h2osno_total,
            snow_depth, qflx_snow_drain,
            patch_column, mask_lakep, bounds_patch,
            dtime, nlevsno)

        # Column 1: snow layer removed
        @test snl[1] == 0
        @test h2osno_no_layers[1] == 0.0
        @test snow_depth[1] == 0.0
        @test qflx_snow_drain[1] > 0.0

        # Column 2: no snow layer, eflx_grnd_lake = eflx_gnet
        @test eflx_grnd_lake[2] ≈ 10.0
    end

    # ------------------------------------------------------------------
    # 7. lake_snow_above_unfrozen!
    # ------------------------------------------------------------------
    @testset "lake_snow_above_unfrozen!" begin
        nc = 2
        nlevtot = nlevsno + nlevmaxurbgrnd
        dtime = 1800.0
        bounds = 1:nc

        t_lake = fill(CLM.TFRZ - 5.0, nc, nlevlak)
        lake_icefrac = fill(0.0, nc, nlevlak)
        t_soisno = fill(CLM.TFRZ - 2.0, nc, nlevtot)
        h2osoi_ice = fill(0.0, nc, nlevtot)
        h2osoi_liq = fill(0.0, nc, nlevtot)
        snl = [-1, 0]  # col 1 has snow, col 2 does not
        h2osno_no_layers = [0.0, 0.0]
        h2osno_total = [1.0, 0.0]
        snow_depth = [0.01, 0.0]
        qflx_snomelt = fill(0.0, nc)
        eflx_snomelt = fill(0.0, nc)
        qflx_snomelt_lyr = fill(0.0, nc, nlevsno)
        qflx_snow_drain = fill(0.0, nc)
        dz_lake = fill(1.0, nc, nlevlak)
        mask_lake = BitVector([true, true])

        # Make column 1 have unfrozen top lake layer with snow above
        t_lake[1, 1] = CLM.TFRZ + 5.0  # warm lake
        lake_icefrac[1, 1] = 0.0
        # Snow layer 0 (Julia index nlevsno)
        j = nlevsno
        h2osoi_ice[1, j] = 0.5  # some snow ice
        h2osoi_liq[1, j] = 0.1
        t_soisno[1, j] = CLM.TFRZ - 1.0

        CLM.lake_snow_above_unfrozen!(
            t_lake, lake_icefrac, t_soisno,
            h2osoi_ice, h2osoi_liq,
            snl, h2osno_no_layers, h2osno_total,
            snow_depth, qflx_snomelt, eflx_snomelt,
            qflx_snomelt_lyr, qflx_snow_drain,
            dz_lake, mask_lake, bounds,
            dtime, nlevsno)

        # Column 1: lake is warm enough to melt snow, snow should be removed
        @test snl[1] == 0
        @test snow_depth[1] == 0.0
        @test qflx_snomelt[1] > 0.0

        # Column 2: no snow, should be unchanged
        @test snl[2] == 0
        @test qflx_snomelt[2] == 0.0
    end

    # ------------------------------------------------------------------
    # 8. lake_snow_diagnostics!
    # ------------------------------------------------------------------
    @testset "lake_snow_diagnostics!" begin
        nc = 2
        nlevtot = nlevsno + nlevmaxurbgrnd
        bounds = 1:nc

        snowice = fill(NaN, nc)
        snowliq = fill(NaN, nc)
        t_sno_mul_mss = fill(NaN, nc)

        h2osoi_ice = fill(0.0, nc, nlevtot)
        h2osoi_liq = fill(0.0, nc, nlevtot)
        t_soisno = fill(CLM.TFRZ - 5.0, nc, nlevtot)
        snl = [-2, 0]  # col 1 has 2 snow layers
        mask_lake = BitVector([true, true])
        mask_lakesnow = BitVector([true, false])

        # Set snow layers for column 1
        # Snow layers: -1 → Julia index nlevsno-1, 0 → Julia index nlevsno
        j_m1 = nlevsno - 1
        j_0  = nlevsno
        h2osoi_ice[1, j_m1] = 3.0
        h2osoi_liq[1, j_m1] = 0.5
        h2osoi_ice[1, j_0]  = 2.0
        h2osoi_liq[1, j_0]  = 0.3

        CLM.lake_snow_diagnostics!(
            snowice, snowliq, t_sno_mul_mss,
            h2osoi_ice, h2osoi_liq, t_soisno,
            snl, mask_lake, mask_lakesnow,
            bounds, nlevsno)

        @test snowice[1] ≈ 5.0   # 3.0 + 2.0
        @test snowliq[1] ≈ 0.8   # 0.5 + 0.3
        @test snowice[2] ≈ 0.0   # no snow
        @test snowliq[2] ≈ 0.0
    end

    # ------------------------------------------------------------------
    # 9. lake_water_balance!
    # ------------------------------------------------------------------
    @testset "lake_water_balance!" begin
        nc = 1
        np = 1
        nlevtot = nlevsno + nlevmaxurbgrnd
        dtime = 1800.0
        bounds_patch = 1:np

        qflx_drain_perched = fill(NaN, nc)
        qflx_rsub_sat = fill(NaN, nc)
        qflx_infl = fill(NaN, nc)
        qflx_surf = fill(NaN, nc)
        qflx_drain = fill(NaN, nc)
        qflx_qrgwl = fill(NaN, nc)
        qflx_floodc = fill(NaN, nc)
        qflx_runoff = fill(NaN, nc)
        qflx_rain_plus_snomelt = fill(NaN, nc)
        qflx_top_soil = fill(NaN, nc)
        qflx_ice_runoff_snwcp = fill(NaN, nc)

        h2osoi_liq = fill(10.0, nc, nlevtot)
        h2osoi_ice = fill(1.0, nc, nlevtot)
        h2osoi_vol = fill(0.3, nc, nlevmaxurbgrnd)

        forc_rain = [2.0]
        forc_snow = [1.0]
        qflx_evap_tot = [0.5]
        qflx_snwcp_ice = [0.1]
        qflx_snwcp_discarded_ice = [0.01]
        qflx_snwcp_discarded_liq = [0.02]
        qflx_liq_grnd = [1.5]
        qflx_snow_drain = [0.3]
        qflx_floodg = [0.5]
        begwb = [100.0]
        endwb = [100.0]  # same as begwb → no change in water storage
        dz = fill(0.1, nc, nlevtot)
        patch_column = [1]
        patch_gridcell = [1]
        mask_lakep = BitVector([true])

        CLM.lake_water_balance!(
            qflx_drain_perched, qflx_rsub_sat, qflx_infl,
            qflx_surf, qflx_drain, qflx_qrgwl,
            qflx_floodc, qflx_runoff,
            qflx_rain_plus_snomelt, qflx_top_soil,
            qflx_ice_runoff_snwcp,
            h2osoi_liq, h2osoi_ice, h2osoi_vol,
            forc_rain, forc_snow, qflx_evap_tot,
            qflx_snwcp_ice, qflx_snwcp_discarded_ice,
            qflx_snwcp_discarded_liq,
            qflx_liq_grnd, qflx_snow_drain,
            qflx_floodg, begwb, endwb,
            dz, patch_column, patch_gridcell,
            mask_lakep, bounds_patch,
            dtime, nlevsno, nlevgrnd)

        @test qflx_drain_perched[1] == 0.0
        @test qflx_rsub_sat[1] == 0.0
        @test qflx_infl[1] == 0.0
        @test qflx_surf[1] == 0.0
        @test qflx_drain[1] == 0.0

        # qflx_qrgwl = rain + snow - evap - snwcp_ice - disc_ice - disc_liq - (endwb-begwb)/dtime + flood
        expected_qrgwl = 2.0 + 1.0 - 0.5 - 0.1 - 0.01 - 0.02 - 0.0 + 0.5
        @test qflx_qrgwl[1] ≈ expected_qrgwl

        @test qflx_runoff[1] ≈ 0.0 + expected_qrgwl  # drain + qrgwl
        @test qflx_rain_plus_snomelt[1] ≈ 1.5 + 0.3
        @test qflx_top_soil[1] ≈ 1.5 + 0.3
        @test qflx_ice_runoff_snwcp[1] ≈ 0.1
    end

    # ------------------------------------------------------------------
    # 10. lake_update_h2osoi_vol!
    # ------------------------------------------------------------------
    @testset "lake_update_h2osoi_vol!" begin
        nc = 1
        nlevtot = nlevsno + nlevmaxurbgrnd
        bounds = 1:nc

        h2osoi_liq = fill(10.0, nc, nlevtot)
        h2osoi_ice = fill(2.0, nc, nlevtot)
        h2osoi_vol = fill(NaN, nc, nlevmaxurbgrnd)
        dz = fill(0.1, nc, nlevtot)
        mask_lake = BitVector([true])

        CLM.lake_update_h2osoi_vol!(
            h2osoi_liq, h2osoi_ice, h2osoi_vol,
            dz, mask_lake, bounds,
            nlevgrnd, nlevsno)

        jj = 1 + nlevsno
        expected_vol = h2osoi_liq[1, jj] / (dz[1, jj] * CLM.DENH2O) +
                       h2osoi_ice[1, jj] / (dz[1, jj] * CLM.DENICE)
        @test h2osoi_vol[1, 1] ≈ expected_vol
    end

    # ------------------------------------------------------------------
    # 11. lake_top_layer_diagnostics!
    # ------------------------------------------------------------------
    @testset "lake_top_layer_diagnostics!" begin
        nc = 2
        nlevtot = nlevsno + nlevmaxurbgrnd
        bounds = 1:nc

        h2osno_top  = fill(NaN, nc)
        snw_rds     = fill(100.0, nc, nlevsno)
        snot_top    = fill(NaN, nc)
        dTdz_top    = fill(NaN, nc)
        snw_rds_top = fill(NaN, nc)
        sno_liq_top = fill(NaN, nc)

        h2osoi_ice = fill(0.0, nc, nlevtot)
        h2osoi_liq = fill(0.0, nc, nlevtot)
        snl = [-1, 0]  # col 1 has 1 snow layer

        mask_lakesnow   = BitVector([true, false])
        mask_lakenosnow = BitVector([false, true])

        # Set snow layer 0 for column 1
        j = nlevsno  # layer 0
        h2osoi_ice[1, j] = 3.0
        h2osoi_liq[1, j] = 0.5

        CLM.lake_top_layer_diagnostics!(
            h2osno_top, snw_rds, snot_top, dTdz_top,
            snw_rds_top, sno_liq_top,
            h2osoi_ice, h2osoi_liq, snl,
            mask_lakesnow, mask_lakenosnow,
            bounds, nlevsno)

        # Column 1: snow present
        @test h2osno_top[1] ≈ 3.5  # 3.0 + 0.5

        # Column 2: no snow — SPVAL
        @test h2osno_top[2] == 0.0
        @test snot_top[2] ≈ CLM.SPVAL
        @test dTdz_top[2] ≈ CLM.SPVAL
        @test snw_rds_top[2] ≈ CLM.SPVAL
        @test sno_liq_top[2] ≈ CLM.SPVAL
        @test snw_rds[2, 1] == 0.0
    end

    # ------------------------------------------------------------------
    # 12. lake_hydrology! (orchestrator smoke test)
    # ------------------------------------------------------------------
    @testset "lake_hydrology! smoke test" begin
        nc = 2
        np = 2
        nl = 1
        ng = 1
        nlevtot = nlevsno + nlevmaxurbgrnd
        dtime = 1800.0

        # Initialize data structures
        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)

        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nl, ng)

        lakestate = CLM.LakeStateData()
        CLM.lakestate_init!(lakestate, nc, np)

        soilstate = CLM.SoilStateData()
        CLM.soilstate_init!(soilstate, np, nc)

        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)

        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)

        waterbalancebulk = CLM.WaterBalanceData()
        CLM.waterbalance_init!(waterbalancebulk, nc, np, ng)

        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)

        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)

        # Setup: 2 lake columns, each with 1 patch
        col_data.landunit .= 1
        col_data.gridcell .= 1
        col_data.itype .= CLM.ISTDLAK
        col_data.snl .= 0  # no snow layers

        patch_data.column .= [1, 2]
        patch_data.gridcell .= [1, 1]
        patch_data.landunit .= [1, 1]

        # Initialize geometry
        for c in 1:nc
            for j in 1:nlevtot
                col_data.dz[c, j] = 0.1
                col_data.z[c, j]  = 0.05 + (j - 1) * 0.1
            end
            for j in 1:nlevlak
                col_data.dz_lake[c, j] = 1.0
                col_data.z_lake[c, j]  = 0.5 + (j - 1) * 1.0
            end
        end

        # Initialize temperatures
        temperature.t_soisno_col .= 270.0
        temperature.t_grnd_col .= 270.0
        temperature.t_lake_col .= 275.0
        temperature.t_sno_mul_mss_col .= 0.0
        temperature.snot_top_col .= CLM.SPVAL
        temperature.dTdz_top_col .= CLM.SPVAL

        # Initialize water state
        ws = waterstatebulk.ws
        ws.h2osoi_liq_col .= 10.0
        ws.h2osoi_ice_col .= 1.0
        ws.h2osoi_vol_col .= 0.3
        ws.h2osno_no_layers_col .= 0.0

        # Initialize water diagnostics
        waterdiagbulk.frac_sno_col .= 0.0
        waterdiagbulk.frac_sno_eff_col .= 0.0
        waterdiagbulk.snow_depth_col .= 0.0
        waterdiagbulk.snowice_col .= 0.0
        waterdiagbulk.snowliq_col .= 0.0
        waterdiagbulk.snw_rds_col .= 0.0
        waterdiagbulk.snw_rds_top_col .= CLM.SPVAL
        waterdiagbulk.h2osno_top_col .= 0.0
        waterdiagbulk.sno_liq_top_col .= CLM.SPVAL

        # Initialize energy fluxes
        energyflux.eflx_sh_tot_patch .= 0.0
        energyflux.eflx_sh_grnd_patch .= 0.0
        energyflux.eflx_soil_grnd_patch .= 0.0
        energyflux.eflx_gnet_patch .= 10.0
        energyflux.eflx_grnd_lake_patch .= 0.0
        energyflux.eflx_snomelt_col .= 0.0

        # Initialize lake state
        lakestate.lake_icefrac_col .= 0.0

        # Initialize soil state
        soilstate.watsat_col .= 0.4

        # Initialize water fluxes
        wf = waterfluxbulk.wf
        wf.qflx_evap_soi_patch .= 0.0001
        wf.qflx_evap_tot_patch .= 0.0001
        wf.qflx_ev_snow_patch .= 0.0
        wf.qflx_solidevap_from_top_layer_patch .= 0.0
        wf.qflx_liqevap_from_top_layer_patch .= 0.0
        wf.qflx_soliddew_to_top_layer_patch .= 0.0
        wf.qflx_liqdew_to_top_layer_patch .= 0.0
        wf.qflx_liq_grnd_col .= 0.0
        wf.qflx_snow_grnd_col .= 0.0
        wf.qflx_evap_tot_col .= 0.0
        wf.qflx_liqevap_from_top_layer_col .= 0.0
        wf.qflx_liqdew_to_top_layer_col .= 0.0
        wf.qflx_soliddew_to_top_layer_col .= 0.0
        wf.qflx_solidevap_from_top_layer_col .= 0.0
        wf.qflx_snomelt_col .= 0.0
        wf.qflx_snow_drain_col .= 0.0
        wf.qflx_snwcp_ice_col .= 0.0
        wf.qflx_snwcp_discarded_ice_col .= 0.0
        wf.qflx_snwcp_discarded_liq_col .= 0.0
        wf.qflx_drain_perched_col .= 0.0
        wf.qflx_rsub_sat_col .= 0.0
        wf.qflx_surf_col .= 0.0
        wf.qflx_drain_col .= 0.0
        wf.qflx_infl_col .= 0.0
        wf.qflx_qrgwl_col .= 0.0
        wf.qflx_runoff_col .= 0.0
        wf.qflx_floodc_col .= 0.0
        wf.qflx_rain_plus_snomelt_col .= 0.0
        wf.qflx_top_soil_col .= 0.0
        wf.qflx_ice_runoff_snwcp_col .= 0.0

        waterfluxbulk.qflx_ev_snow_col .= 0.0
        waterfluxbulk.qflx_snomelt_lyr_col .= 0.0

        # Water balance
        waterbalancebulk.begwb_col .= 100.0
        waterbalancebulk.endwb_col .= 0.0

        mask_lake  = BitVector([true, true])
        mask_lakep = BitVector([true, true])
        bounds_col = 1:nc
        bounds_patch = 1:np

        forc_rain  = fill(0.001, nc)
        forc_snow  = fill(0.0005, nc)
        qflx_floodg = [0.0]

        # Run the orchestrator
        CLM.lake_hydrology!(
            temperature, energyflux, lakestate, soilstate,
            waterstatebulk, waterdiagbulk, waterbalancebulk, waterfluxbulk,
            col_data, patch_data,
            mask_lake, mask_lakep,
            forc_rain, forc_snow, qflx_floodg,
            bounds_col, bounds_patch,
            dtime, nlevsno, nlevsoi, nlevgrnd)

        # Verify basic results
        # Precip fluxes should be copied to ground
        @test wf.qflx_snow_grnd_col[1] ≈ 0.0005
        @test wf.qflx_liq_grnd_col[1] ≈ 0.001

        # Drainage should be zero for lakes
        @test wf.qflx_drain_col[1] == 0.0
        @test wf.qflx_surf_col[1] == 0.0
        @test wf.qflx_infl_col[1] == 0.0

        # qflx_qrgwl should balance the water budget
        @test !isnan(wf.qflx_qrgwl_col[1])

        # eflx_grnd_lake should be set
        @test !isnan(energyflux.eflx_grnd_lake_patch[1])

        # Soil should be at/near saturation
        joff = nlevsno
        vol = ws.h2osoi_vol_col[1, 1]
        @test !isnan(vol)
    end
end
