@testset "Soil Fluxes" begin

    @testset "soil_fluxes! single non-urban patch smoke test" begin
        # Initialize varpar
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        np = 1   # 1 patch
        nc = 1   # 1 column
        ng = 1   # 1 gridcell
        nl = 1   # 1 landunit
        nlevgrnd = CLM.varpar.nlevgrnd
        nlevsno  = CLM.varpar.nlevsno
        nlevurb  = CLM.varpar.nlevurb
        nlev_soisno = nlevsno + CLM.varpar.nlevmaxurbgrnd

        dtime = 1800.0  # 30 min time step

        # --- Set up all required data structures ---
        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nl, ng)

        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)

        canopystate = CLM.CanopyStateData()
        CLM.canopystate_init!(canopystate, np)

        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)

        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)

        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)

        solarabs = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs, np, nl)

        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        patch_data.column[1] = 1
        patch_data.gridcell[1] = 1
        patch_data.landunit[1] = 1
        patch_data.itype[1] = 1
        patch_data.active[1] = true
        patch_data.wtcol[1] = 1.0

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)
        col_data.snl[1] = 0          # no snow layers
        col_data.itype[1] = 1        # soil column
        col_data.landunit[1] = 1
        col_data.patchi[1] = 1
        col_data.patchf[1] = 1

        lun_data = CLM.LandunitData()
        CLM.landunit_init!(lun_data, nl)
        lun_data.itype[1] = CLM.ISTSOIL
        lun_data.urbpoi[1] = false

        # Set temperature values
        temperature.t_grnd_col[1] = 280.0
        temperature.t_h2osfc_col[1] = 280.0
        temperature.t_h2osfc_bef_col[1] = 279.5
        temperature.emg_col[1] = 0.96
        temperature.xmf_col[1] = 0.0
        temperature.xmf_h2osfc_col[1] = 0.0
        temperature.c_h2osfc_col[1] = 4.188e6
        temperature.t_skin_patch[1] = 280.0
        # t_ssbef and t_soisno: set all levels
        for j in 1:nlev_soisno
            temperature.t_ssbef_col[1, j] = 279.5
            temperature.t_soisno_col[1, j] = 280.0
            temperature.fact_col[1, j] = 1.0e7  # large value to make errsoi term small
        end

        # Water state: some liquid, some ice in top soil layer
        for j in 1:size(waterstatebulk.ws.h2osoi_liq_col, 2)
            waterstatebulk.ws.h2osoi_liq_col[1, j] = 20.0
            waterstatebulk.ws.h2osoi_ice_col[1, j] = 5.0
        end

        # Water diagnostics
        waterdiagbulk.frac_sno_eff_col[1] = 0.0
        waterdiagbulk.frac_h2osfc_col[1] = 0.0

        # Energy flux pre-conditions (set by earlier routines)
        energyflux.cgrnds_patch[1] = 10.0    # dSH/dT
        energyflux.cgrndl_patch[1] = 0.001   # dLH/dT
        energyflux.htvp_col[1] = CLM.HVAP
        energyflux.eflx_sh_grnd_patch[1] = 50.0
        energyflux.eflx_sh_veg_patch[1] = 10.0
        energyflux.eflx_sh_stem_patch[1] = 2.0
        energyflux.dlrad_patch[1] = 300.0
        energyflux.ulrad_patch[1] = 350.0
        energyflux.eflx_h2osfc_to_snow_col[1] = 0.0
        energyflux.eflx_building_heat_errsoi_col[1] = 0.0
        energyflux.eflx_lwrad_net_patch[1] = 0.0
        energyflux.eflx_lwrad_out_patch[1] = 0.0
        energyflux.eflx_wasteheat_patch[1] = 0.0
        energyflux.eflx_heat_from_ac_patch[1] = 0.0
        energyflux.eflx_traffic_patch[1] = 0.0
        energyflux.eflx_ventilation_patch[1] = 0.0
        energyflux.eflx_lwrad_net_u_patch[1] = 0.0
        energyflux.eflx_lwrad_out_u_patch[1] = 0.0

        # Water flux pre-conditions
        waterfluxbulk.wf.qflx_evap_soi_patch[1] = 1.0e-5
        waterfluxbulk.wf.qflx_evap_veg_patch[1] = 2.0e-5
        waterfluxbulk.wf.qflx_tran_veg_patch[1] = 1.0e-5
        waterfluxbulk.wf.qflx_evap_tot_patch[1] = 3.0e-5
        waterfluxbulk.wf.qflx_evap_can_patch[1] = 1.0e-5
        waterfluxbulk.qflx_ev_snow_patch[1] = 5.0e-6
        waterfluxbulk.qflx_ev_soil_patch[1] = 5.0e-6
        waterfluxbulk.qflx_ev_h2osfc_patch[1] = 0.0

        # Canopy state
        canopystate.frac_veg_nosno_patch[1] = 0  # bare ground

        # Solar absorbed
        solarabs.sabg_soil_patch[1] = 100.0
        solarabs.sabg_snow_patch[1] = 0.0
        solarabs.sabg_patch[1] = 100.0

        # Atmospheric forcing
        forc_lwrad_col = [300.0]

        # Masks
        mask_nolakec = trues(nc)
        mask_nolakep = trues(np)
        mask_urbanp  = falses(np)

        bounds_col   = 1:nc
        bounds_patch = 1:np

        # Call soil_fluxes!
        CLM.soil_fluxes!(
            energyflux, temperature, canopystate,
            waterstatebulk, waterdiagbulk, waterfluxbulk, solarabs,
            patch_data, col_data, lun_data,
            mask_nolakec, mask_nolakep, mask_urbanp,
            bounds_col, bounds_patch,
            forc_lwrad_col, dtime)

        # --- Finiteness checks ---
        @test isfinite(energyflux.eflx_sh_grnd_patch[1])
        @test isfinite(energyflux.eflx_sh_tot_patch[1])
        @test isfinite(energyflux.eflx_soil_grnd_patch[1])
        @test isfinite(energyflux.eflx_soil_grnd_r_patch[1])
        @test isfinite(energyflux.eflx_lh_tot_patch[1])
        @test isfinite(energyflux.eflx_lh_tot_r_patch[1])
        @test isfinite(energyflux.eflx_sh_tot_r_patch[1])
        @test isfinite(energyflux.eflx_lwrad_out_patch[1])
        @test isfinite(energyflux.eflx_lwrad_net_patch[1])
        @test isfinite(energyflux.eflx_lwrad_net_r_patch[1])
        @test isfinite(energyflux.eflx_lwrad_out_r_patch[1])
        @test isfinite(energyflux.eflx_lh_vege_patch[1])
        @test isfinite(energyflux.eflx_lh_vegt_patch[1])
        @test isfinite(energyflux.eflx_lh_grnd_patch[1])
        @test isfinite(energyflux.errsoi_patch[1])
        @test isfinite(energyflux.errsoi_col[1])

        @test isfinite(waterfluxbulk.wf.qflx_evap_soi_patch[1])
        @test isfinite(waterfluxbulk.wf.qflx_evap_tot_patch[1])
        @test isfinite(waterfluxbulk.wf.qflx_evap_can_patch[1])
        @test isfinite(waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[1])
        @test isfinite(waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[1])

        @test isfinite(temperature.t_skin_patch[1])

        # --- Physical reasonableness checks ---

        # Total sensible heat = veg + ground + stem for non-urban
        @test energyflux.eflx_sh_tot_patch[1] ≈
            energyflux.eflx_sh_veg_patch[1] + energyflux.eflx_sh_grnd_patch[1] + energyflux.eflx_sh_stem_patch[1]

        # Total evaporation = vegetation + soil
        @test waterfluxbulk.wf.qflx_evap_tot_patch[1] ≈
            waterfluxbulk.wf.qflx_evap_veg_patch[1] + waterfluxbulk.wf.qflx_evap_soi_patch[1]

        # Canopy evaporation = veg evap - transpiration
        @test waterfluxbulk.wf.qflx_evap_can_patch[1] ≈
            waterfluxbulk.wf.qflx_evap_veg_patch[1] - waterfluxbulk.wf.qflx_tran_veg_patch[1]

        # Latent heat of transpiration
        @test energyflux.eflx_lh_vegt_patch[1] ≈
            waterfluxbulk.wf.qflx_tran_veg_patch[1] * CLM.HVAP

        # Rural soil heat flux should match total for istsoil
        @test energyflux.eflx_soil_grnd_r_patch[1] == energyflux.eflx_soil_grnd_patch[1]

        # Rural total latent heat should match total for istsoil
        @test energyflux.eflx_lh_tot_r_patch[1] == energyflux.eflx_lh_tot_patch[1]

        # Skin temperature should be positive and reasonable for bare ground
        @test temperature.t_skin_patch[1] > 200.0
        @test temperature.t_skin_patch[1] < 400.0

        # Column errsoi should equal patch errsoi (single patch, wtcol=1)
        @test energyflux.errsoi_col[1] ≈ energyflux.errsoi_patch[1]

        # Evaporation partitioning: liquid + solid should sum to total
        # (for positive evaporation, dew terms are zero)
        if waterfluxbulk.qflx_ev_snow_patch[1] >= 0.0
            @test waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[1] +
                waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[1] ≈
                waterfluxbulk.qflx_ev_snow_patch[1]
        end
    end

    @testset "soil_fluxes! with empty masks" begin
        # Test with no active columns/patches — should be a no-op
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        np = 1; nc = 1; ng = 1; nl = 1
        nlevsno = CLM.varpar.nlevsno

        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nl, ng)
        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)
        canopystate = CLM.CanopyStateData()
        CLM.canopystate_init!(canopystate, np)
        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)
        solarabs = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs, np, nl)
        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)
        lun_data = CLM.LandunitData()
        CLM.landunit_init!(lun_data, nl)

        # Empty masks
        mask_nolakec = falses(nc)
        mask_nolakep = falses(np)
        mask_urbanp  = falses(np)

        # This should run without error
        CLM.soil_fluxes!(
            energyflux, temperature, canopystate,
            waterstatebulk, waterdiagbulk, waterfluxbulk, solarabs,
            patch_data, col_data, lun_data,
            mask_nolakec, mask_nolakep, mask_urbanp,
            1:nc, 1:np,
            [300.0], 1800.0)

        @test true  # if we get here, no error
    end

    @testset "soil_fluxes! dew formation (negative evaporation, t_grnd < TFRZ)" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        np = 1; nc = 1; ng = 1; nl = 1
        nlevsno = CLM.varpar.nlevsno
        nlev_soisno = nlevsno + CLM.varpar.nlevmaxurbgrnd
        dtime = 1800.0

        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nl, ng)
        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)
        canopystate = CLM.CanopyStateData()
        CLM.canopystate_init!(canopystate, np)
        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)
        solarabs = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs, np, nl)

        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        patch_data.column[1] = 1
        patch_data.gridcell[1] = 1
        patch_data.landunit[1] = 1
        patch_data.itype[1] = 1
        patch_data.active[1] = true
        patch_data.wtcol[1] = 1.0

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)
        col_data.snl[1] = 0
        col_data.itype[1] = 1
        col_data.landunit[1] = 1
        col_data.patchi[1] = 1
        col_data.patchf[1] = 1

        lun_data = CLM.LandunitData()
        CLM.landunit_init!(lun_data, nl)
        lun_data.itype[1] = CLM.ISTSOIL
        lun_data.urbpoi[1] = false

        # Cold ground — below freezing
        temperature.t_grnd_col[1] = 265.0
        temperature.t_h2osfc_col[1] = 265.0
        temperature.t_h2osfc_bef_col[1] = 265.0
        temperature.emg_col[1] = 0.96
        temperature.xmf_col[1] = 0.0
        temperature.xmf_h2osfc_col[1] = 0.0
        temperature.c_h2osfc_col[1] = 4.188e6
        temperature.t_skin_patch[1] = 265.0
        for j in 1:nlev_soisno
            temperature.t_ssbef_col[1, j] = 265.0
            temperature.t_soisno_col[1, j] = 265.0
            temperature.fact_col[1, j] = 1.0e7
        end

        for j in 1:size(waterstatebulk.ws.h2osoi_liq_col, 2)
            waterstatebulk.ws.h2osoi_liq_col[1, j] = 10.0
            waterstatebulk.ws.h2osoi_ice_col[1, j] = 10.0
        end

        waterdiagbulk.frac_sno_eff_col[1] = 0.0
        waterdiagbulk.frac_h2osfc_col[1] = 0.0

        energyflux.cgrnds_patch[1] = 10.0
        energyflux.cgrndl_patch[1] = 0.001
        energyflux.htvp_col[1] = CLM.HSUB  # sublimation
        energyflux.eflx_sh_grnd_patch[1] = 30.0
        energyflux.eflx_sh_veg_patch[1] = 5.0
        energyflux.eflx_sh_stem_patch[1] = 1.0
        energyflux.dlrad_patch[1] = 250.0
        energyflux.ulrad_patch[1] = 280.0
        energyflux.eflx_h2osfc_to_snow_col[1] = 0.0
        energyflux.eflx_building_heat_errsoi_col[1] = 0.0
        energyflux.eflx_lwrad_net_patch[1] = 0.0
        energyflux.eflx_lwrad_out_patch[1] = 0.0

        canopystate.frac_veg_nosno_patch[1] = 0

        # Negative evaporation (dew/frost)
        waterfluxbulk.wf.qflx_evap_soi_patch[1] = -1.0e-5
        waterfluxbulk.wf.qflx_evap_veg_patch[1] = 0.0
        waterfluxbulk.wf.qflx_tran_veg_patch[1] = 0.0
        waterfluxbulk.qflx_ev_snow_patch[1] = -1.0e-5
        waterfluxbulk.qflx_ev_soil_patch[1] = -1.0e-5
        waterfluxbulk.qflx_ev_h2osfc_patch[1] = 0.0

        solarabs.sabg_soil_patch[1] = 50.0
        solarabs.sabg_snow_patch[1] = 0.0
        solarabs.sabg_patch[1] = 50.0

        CLM.soil_fluxes!(
            energyflux, temperature, canopystate,
            waterstatebulk, waterdiagbulk, waterfluxbulk, solarabs,
            patch_data, col_data, lun_data,
            trues(nc), trues(np), falses(np),
            1:nc, 1:np,
            [250.0], dtime)

        # Below freezing: dew should go to solid (frost)
        @test waterfluxbulk.wf.qflx_soliddew_to_top_layer_patch[1] > 0.0
        @test waterfluxbulk.wf.qflx_liqdew_to_top_layer_patch[1] == 0.0
        @test waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[1] == 0.0
        @test waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[1] == 0.0

        @test isfinite(energyflux.errsoi_patch[1])
        @test isfinite(energyflux.eflx_soil_grnd_patch[1])
    end

    @testset "soil_fluxes! with snow layers" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        np = 1; nc = 1; ng = 1; nl = 1
        nlevsno = CLM.varpar.nlevsno
        nlev_soisno = nlevsno + CLM.varpar.nlevmaxurbgrnd
        dtime = 1800.0

        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nl, ng)
        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)
        canopystate = CLM.CanopyStateData()
        CLM.canopystate_init!(canopystate, np)
        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)
        solarabs = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs, np, nl)

        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        patch_data.column[1] = 1
        patch_data.gridcell[1] = 1
        patch_data.landunit[1] = 1
        patch_data.itype[1] = 1
        patch_data.active[1] = true
        patch_data.wtcol[1] = 1.0

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)
        col_data.snl[1] = -3         # 3 snow layers
        col_data.itype[1] = 1
        col_data.landunit[1] = 1
        col_data.patchi[1] = 1
        col_data.patchf[1] = 1

        lun_data = CLM.LandunitData()
        CLM.landunit_init!(lun_data, nl)
        lun_data.itype[1] = CLM.ISTSOIL
        lun_data.urbpoi[1] = false

        temperature.t_grnd_col[1] = 270.0
        temperature.t_h2osfc_col[1] = 270.0
        temperature.t_h2osfc_bef_col[1] = 269.5
        temperature.emg_col[1] = 0.97
        temperature.xmf_col[1] = 0.0
        temperature.xmf_h2osfc_col[1] = 0.0
        temperature.c_h2osfc_col[1] = 4.188e6
        temperature.t_skin_patch[1] = 270.0
        for j in 1:nlev_soisno
            temperature.t_ssbef_col[1, j] = 269.5
            temperature.t_soisno_col[1, j] = 270.0
            temperature.fact_col[1, j] = 1.0e7
        end

        for j in 1:size(waterstatebulk.ws.h2osoi_liq_col, 2)
            waterstatebulk.ws.h2osoi_liq_col[1, j] = 5.0
            waterstatebulk.ws.h2osoi_ice_col[1, j] = 50.0
        end

        waterdiagbulk.frac_sno_eff_col[1] = 0.8
        waterdiagbulk.frac_h2osfc_col[1] = 0.0

        energyflux.cgrnds_patch[1] = 8.0
        energyflux.cgrndl_patch[1] = 0.0005
        energyflux.htvp_col[1] = CLM.HSUB
        energyflux.eflx_sh_grnd_patch[1] = 20.0
        energyflux.eflx_sh_veg_patch[1] = 3.0
        energyflux.eflx_sh_stem_patch[1] = 0.5
        energyflux.dlrad_patch[1] = 220.0
        energyflux.ulrad_patch[1] = 240.0
        energyflux.eflx_h2osfc_to_snow_col[1] = 0.0
        energyflux.eflx_building_heat_errsoi_col[1] = 0.0
        energyflux.eflx_lwrad_net_patch[1] = 0.0
        energyflux.eflx_lwrad_out_patch[1] = 0.0

        canopystate.frac_veg_nosno_patch[1] = 0

        waterfluxbulk.wf.qflx_evap_soi_patch[1] = 5.0e-6
        waterfluxbulk.wf.qflx_evap_veg_patch[1] = 0.0
        waterfluxbulk.wf.qflx_tran_veg_patch[1] = 0.0
        waterfluxbulk.qflx_ev_snow_patch[1] = 5.0e-6
        waterfluxbulk.qflx_ev_soil_patch[1] = 0.0
        waterfluxbulk.qflx_ev_h2osfc_patch[1] = 0.0

        solarabs.sabg_soil_patch[1] = 30.0
        solarabs.sabg_snow_patch[1] = 80.0
        solarabs.sabg_patch[1] = 110.0

        CLM.soil_fluxes!(
            energyflux, temperature, canopystate,
            waterstatebulk, waterdiagbulk, waterfluxbulk, solarabs,
            patch_data, col_data, lun_data,
            trues(nc), trues(np), falses(np),
            1:nc, 1:np,
            [280.0], dtime)

        # With snow layers (snl=-3), t_grnd0 should use snow fraction weighting
        @test isfinite(energyflux.eflx_soil_grnd_patch[1])
        @test isfinite(energyflux.errsoi_patch[1])
        @test isfinite(energyflux.eflx_lwrad_out_patch[1])
        @test isfinite(energyflux.eflx_lwrad_net_patch[1])
        @test isfinite(temperature.t_skin_patch[1])

        # Evaporation from snow layers should be constrained
        @test waterfluxbulk.qflx_ev_snow_patch[1] >= 0.0

        # Partitioning should be consistent
        @test isfinite(waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[1])
        @test isfinite(waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[1])
    end
end
