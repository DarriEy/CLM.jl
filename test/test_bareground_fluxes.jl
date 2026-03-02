@testset "Bare Ground Fluxes" begin

    @testset "BareGroundFluxesParamsData" begin
        bgf = CLM.BareGroundFluxesParamsData()
        @test bgf.a_coef == 0.0
        @test bgf.a_exp == 0.0
        @test bgf.wind_min == 0.0
    end

    @testset "bareground_fluxes_read_params!" begin
        CLM.bareground_fluxes_read_params!(
            a_coef=0.5, a_exp=1.0, wind_min=1.0)
        @test CLM.bareground_fluxes_params.a_coef == 0.5
        @test CLM.bareground_fluxes_params.a_exp == 1.0
        @test CLM.bareground_fluxes_params.wind_min == 1.0

        # Reset to defaults
        CLM.bareground_fluxes_read_params!(
            a_coef=0.0, a_exp=0.0, wind_min=0.0)
    end

    @testset "bareground_fluxes! single patch smoke test" begin
        # Initialize varpar for nlevgrnd
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        np = 1   # 1 patch
        nc = 1   # 1 column
        ng = 1   # 1 gridcell
        nl = 1   # 1 landunit
        nlevgrnd = CLM.varpar.nlevgrnd
        nlevsno  = CLM.varpar.nlevsno

        # Set parameters for ZengWang2007 method
        CLM.bareground_fluxes_read_params!(
            a_coef=0.5, a_exp=1.0, wind_min=1.0)

        # --- Set up all required data structures ---
        canopystate = CLM.CanopyStateData()
        CLM.canopystate_init!(canopystate, np)
        canopystate.displa_patch[1] = 0.0

        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nl, ng)
        energyflux.htvp_col[1] = CLM.HVAP
        energyflux.btran_patch[1] = 0.0
        energyflux.cgrnds_patch[1] = 0.0
        energyflux.cgrndl_patch[1] = 0.0
        energyflux.cgrnd_patch[1] = 0.0
        energyflux.dhsdt_canopy_patch[1] = 0.0
        energyflux.eflx_sh_stem_patch[1] = 0.0
        for j in 1:nlevgrnd
            energyflux.rresis_patch[1, j] = 0.0
        end

        frictionvel = CLM.FrictionVelocityData()
        CLM.frictionvel_init!(frictionvel, np, nc)
        frictionvel.zetamaxstable = 0.5
        frictionvel.zsno = 0.00085
        frictionvel.zlnd = 0.000775
        frictionvel.z0mg_col[1] = 0.01
        frictionvel.z0hg_col[1] = 0.01
        frictionvel.z0qg_col[1] = 0.01
        frictionvel.z0mg_patch[1] = 0.01
        frictionvel.z0hg_patch[1] = 0.01
        frictionvel.z0qg_patch[1] = 0.01
        frictionvel.z0mv_patch[1] = 0.0
        frictionvel.z0hv_patch[1] = 0.0
        frictionvel.z0qv_patch[1] = 0.0
        frictionvel.forc_hgt_u_patch[1] = 30.0
        frictionvel.forc_hgt_t_patch[1] = 30.0
        frictionvel.forc_hgt_q_patch[1] = 30.0
        frictionvel.ustar_patch[1] = 0.5
        frictionvel.um_patch[1] = 5.0
        frictionvel.obu_patch[1] = -100.0
        frictionvel.zeta_patch[1] = -0.1
        frictionvel.ram1_patch[1] = 50.0
        frictionvel.num_iter_patch[1] = 0.0
        frictionvel.kbm1_patch[1] = 0.0

        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)
        temperature.t_veg_patch[1] = 290.0
        temperature.thm_patch[1] = 290.0
        temperature.t_grnd_col[1] = 288.0
        temperature.t_h2osfc_col[1] = 288.0
        temperature.thv_col[1] = 291.0
        temperature.beta_col[1] = 1.0
        temperature.t_ref2m_patch[1] = 290.0
        temperature.t_ref2m_r_patch[1] = 290.0
        # t_soisno_col: combined snow+soil indexing
        for j in 1:size(temperature.t_soisno_col, 2)
            temperature.t_soisno_col[1, j] = 288.0
        end

        soilstate = CLM.SoilStateData()
        CLM.soilstate_init!(soilstate, np, nc)
        soilstate.soilbeta_col[1] = 0.8
        soilstate.soilresis_col[1] = 100.0
        for j in 1:nlevgrnd
            soilstate.watsat_col[1, j] = 0.45
            soilstate.rootr_patch[1, j] = 0.0
        end

        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)
        waterfluxbulk.wf.qflx_tran_veg_patch[1] = 0.0
        waterfluxbulk.wf.qflx_evap_veg_patch[1] = 0.0
        waterfluxbulk.wf.qflx_evap_soi_patch[1] = 0.0
        waterfluxbulk.wf.qflx_evap_tot_patch[1] = 0.0

        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
        # First soil layer: some liquid water, no ice
        waterstatebulk.ws.h2osoi_liq_col[1, 1 + nlevsno] = 50.0   # kg/m2
        waterstatebulk.ws.h2osoi_ice_col[1, 1 + nlevsno] = 0.0    # kg/m2

        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
        waterdiagbulk.qg_col[1] = 0.005
        waterdiagbulk.qg_snow_col[1] = 0.005
        waterdiagbulk.qg_soil_col[1] = 0.005
        waterdiagbulk.qg_h2osfc_col[1] = 0.005
        waterdiagbulk.dqgdT_col[1] = 0.0003
        waterdiagbulk.q_ref2m_patch[1] = 0.008
        waterdiagbulk.rh_ref2m_patch[1] = 50.0
        waterdiagbulk.rh_ref2m_r_patch[1] = 50.0

        photosyns = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(photosyns, np)
        photosyns.rssun_patch[1] = 200.0
        photosyns.rssha_patch[1] = 300.0

        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        patch_data.column[1] = 1
        patch_data.gridcell[1] = 1
        patch_data.landunit[1] = 1
        patch_data.itype[1] = 1
        patch_data.active[1] = true

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)
        col_data.snl[1] = 0
        col_data.zii[1] = 1000.0
        col_data.dz[1, 1 + nlevsno] = 0.1   # first soil layer thickness

        lun_data = CLM.LandunitData()
        CLM.landunit_init!(lun_data, nl)
        lun_data.itype[1] = CLM.ISTSOIL

        mask_noexposedveg = falses(np)
        mask_noexposedveg[1] = true

        bounds_patch = 1:np

        # Atmospheric forcing
        forc_q_col     = [0.008]
        forc_pbot_col  = [101325.0]
        forc_th_col    = [290.0]
        forc_rho_col   = [1.2]
        forc_t_col     = [288.0]
        forc_u_grc     = [3.0]
        forc_v_grc     = [1.0]
        forc_hgt_t_grc = [10.0]
        forc_hgt_u_grc = [10.0]
        forc_hgt_q_grc = [10.0]

        # Set up surface resistance method (Lee-Pielke 1992)
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_LEEPIELKE_1992)

        # Call bareground_fluxes!
        CLM.bareground_fluxes!(
            canopystate, energyflux, frictionvel, temperature,
            soilstate, waterfluxbulk, waterstatebulk, waterdiagbulk,
            photosyns, patch_data, col_data, lun_data,
            mask_noexposedveg, bounds_patch,
            forc_q_col, forc_pbot_col, forc_th_col,
            forc_rho_col, forc_t_col,
            forc_u_grc, forc_v_grc,
            forc_hgt_t_grc, forc_hgt_u_grc, forc_hgt_q_grc;
            use_lch4=false,
            z0param_method="ZengWang2007")

        # --- Finiteness checks ---
        @test isfinite(energyflux.eflx_sh_grnd_patch[1])
        @test isfinite(energyflux.eflx_sh_tot_patch[1])
        @test isfinite(energyflux.eflx_sh_snow_patch[1])
        @test isfinite(energyflux.eflx_sh_soil_patch[1])
        @test isfinite(energyflux.eflx_sh_h2osfc_patch[1])
        @test isfinite(energyflux.taux_patch[1])
        @test isfinite(energyflux.tauy_patch[1])
        @test isfinite(energyflux.cgrnds_patch[1])
        @test isfinite(energyflux.cgrndl_patch[1])
        @test isfinite(energyflux.cgrnd_patch[1])
        @test isfinite(temperature.t_ref2m_patch[1])
        @test isfinite(temperature.t_ref2m_r_patch[1])
        @test isfinite(waterdiagbulk.q_ref2m_patch[1])
        @test isfinite(waterdiagbulk.rh_ref2m_patch[1])
        @test isfinite(waterdiagbulk.rh_ref2m_r_patch[1])
        @test isfinite(frictionvel.ustar_patch[1])
        @test isfinite(frictionvel.ram1_patch[1])
        @test isfinite(frictionvel.kbm1_patch[1])
        @test isfinite(waterfluxbulk.wf.qflx_evap_soi_patch[1])
        @test isfinite(waterfluxbulk.wf.qflx_evap_tot_patch[1])
        @test isfinite(waterfluxbulk.qflx_ev_snow_patch[1])
        @test isfinite(waterfluxbulk.qflx_ev_soil_patch[1])
        @test isfinite(waterfluxbulk.qflx_ev_h2osfc_patch[1])

        # --- Physical reasonableness checks ---
        # 2m air temperature should be near ambient (~288-292 K)
        @test temperature.t_ref2m_patch[1] > 200.0
        @test temperature.t_ref2m_patch[1] < 400.0

        # Friction velocity should be positive
        @test frictionvel.ustar_patch[1] > 0.0

        # Aerodynamic resistance should be positive
        @test frictionvel.ram1_patch[1] > 0.0

        # Sensible heat derivative should be positive
        @test energyflux.cgrnds_patch[1] > 0.0

        # For bare ground: transpiration and vegetation evaporation should be zero
        @test waterfluxbulk.wf.qflx_tran_veg_patch[1] == 0.0
        @test waterfluxbulk.wf.qflx_evap_veg_patch[1] == 0.0

        # btran should be zero for bare ground
        @test energyflux.btran_patch[1] == 0.0

        # Displa should be zero (bare ground)
        @test canopystate.displa_patch[1] == 0.0

        # dlrad and ulrad should be zero (no canopy)
        @test energyflux.dlrad_patch[1] == 0.0
        @test energyflux.ulrad_patch[1] == 0.0

        # Total evaporation should equal soil evaporation for bare ground
        @test waterfluxbulk.wf.qflx_evap_tot_patch[1] == waterfluxbulk.wf.qflx_evap_soi_patch[1]

        # Total sensible heat should equal ground sensible heat for bare ground
        @test energyflux.eflx_sh_tot_patch[1] == energyflux.eflx_sh_grnd_patch[1]

        # Wind stress should have correct sign (opposing wind direction)
        @test energyflux.taux_patch[1] < 0.0  # opposing positive forc_u
        @test energyflux.tauy_patch[1] < 0.0  # opposing positive forc_v

        # Relative humidity should be between 0 and 100
        @test waterdiagbulk.rh_ref2m_patch[1] >= 0.0
        @test waterdiagbulk.rh_ref2m_patch[1] <= 100.0

        # For istsoil: rural values should match non-rural
        @test waterdiagbulk.rh_ref2m_r_patch[1] == waterdiagbulk.rh_ref2m_patch[1]
        @test temperature.t_ref2m_r_patch[1] == temperature.t_ref2m_patch[1]

        # Number of iterations should be 3 (niters)
        @test frictionvel.num_iter_patch[1] == 3.0

        # rootr and rresis should be zero for bare ground
        for j in 1:nlevgrnd
            @test soilstate.rootr_patch[1, j] == 0.0
            @test energyflux.rresis_patch[1, j] == 0.0
        end

        # Reset params to defaults
        CLM.bareground_fluxes_read_params!(
            a_coef=0.0, a_exp=0.0, wind_min=0.0)
    end

    @testset "bareground_fluxes! with empty mask" begin
        # Test with no active patches — should be a no-op
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        np = 1; nc = 1; ng = 1; nl = 1
        nlevgrnd = CLM.varpar.nlevgrnd
        nlevsno  = CLM.varpar.nlevsno

        CLM.bareground_fluxes_read_params!(
            a_coef=0.5, a_exp=1.0, wind_min=1.0)

        canopystate = CLM.CanopyStateData()
        CLM.canopystate_init!(canopystate, np)
        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nl, ng)
        frictionvel = CLM.FrictionVelocityData()
        CLM.frictionvel_init!(frictionvel, np, nc)
        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)
        soilstate = CLM.SoilStateData()
        CLM.soilstate_init!(soilstate, np, nc)
        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)
        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
        photosyns = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(photosyns, np)
        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)
        lun_data = CLM.LandunitData()
        CLM.landunit_init!(lun_data, nl)

        # Empty mask — no active patches
        mask_empty = falses(np)
        bounds_patch = 1:np

        # This should run without error
        CLM.bareground_fluxes!(
            canopystate, energyflux, frictionvel, temperature,
            soilstate, waterfluxbulk, waterstatebulk, waterdiagbulk,
            photosyns, patch_data, col_data, lun_data,
            mask_empty, bounds_patch,
            [0.008], [101325.0], [290.0], [1.2], [288.0],
            [3.0], [1.0], [10.0], [10.0], [10.0];
            use_lch4=false, z0param_method="ZengWang2007")

        @test true  # if we get here, no error
    end
end
