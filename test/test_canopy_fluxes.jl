@testset "Canopy Fluxes" begin

    @testset "Constants" begin
        @test CLM.BTRAN0 == 0.0
        @test CLM.ZII_CANOPY == 1000.0
        @test CLM.BETA_CANOPY == 1.0
        @test CLM.DELMAX_CANOPY == 1.0
        @test CLM.DLEMIN_CANOPY == 0.1
        @test CLM.DTMIN_CANOPY == 0.01
        @test CLM.ITMIN_CANOPY == 2
        @test CLM.RIA_CANOPY == 0.5
        @test CLM.ABOVE_CANOPY == 1
        @test CLM.BELOW_CANOPY == 2
        @test CLM.K_VERT_CANOPY == 0.1
        @test CLM.K_CYL_VOL_CANOPY == 1.0
        @test CLM.K_CYL_AREA_CANOPY == 1.0
        @test CLM.K_INTERNAL_CANOPY == 0.0
        @test CLM.MIN_STEM_DIAMETER == 0.05
    end

    @testset "CanopyFluxesParamsData" begin
        cf = CLM.CanopyFluxesParamsData()
        @test cf.lai_dl == 0.5
        @test cf.z_dl == 0.05
        @test cf.a_coef == 0.5
        @test cf.a_exp == 1.0
        @test cf.csoilc == 0.004
        @test cf.cv == 0.01
        @test cf.wind_min == 1.0
    end

    @testset "CanopyFluxesControl" begin
        ctrl = CLM.CanopyFluxesControl()
        @test ctrl.perchroot == false
        @test ctrl.perchroot_alt == false
        @test ctrl.use_undercanopy_stability == false
        @test ctrl.itmax_canopy_fluxes == 40
        @test ctrl.use_biomass_heat_storage == false
    end

    @testset "canopy_fluxes_read_nml!" begin
        # Reset to defaults first
        CLM.canopy_fluxes_read_nml!(
            use_undercanopy_stability=true,
            use_biomass_heat_storage=true,
            itmax_canopy_fluxes=50)
        @test CLM.canopy_fluxes_ctrl.use_undercanopy_stability == true
        @test CLM.canopy_fluxes_ctrl.use_biomass_heat_storage == true
        @test CLM.canopy_fluxes_ctrl.itmax_canopy_fluxes == 50

        # Error case
        @test_throws ErrorException CLM.canopy_fluxes_read_nml!(itmax_canopy_fluxes=0)

        # Reset back to defaults
        CLM.canopy_fluxes_read_nml!(
            use_undercanopy_stability=false,
            use_biomass_heat_storage=false,
            itmax_canopy_fluxes=40)
    end

    @testset "canopy_fluxes_read_params!" begin
        CLM.canopy_fluxes_read_params!(
            lai_dl=0.6, z_dl=0.06, a_coef=0.6,
            a_exp=1.1, csoilc=0.005, cv=0.02, wind_min=0.5)
        @test CLM.canopy_fluxes_params.lai_dl == 0.6
        @test CLM.canopy_fluxes_params.z_dl == 0.06
        @test CLM.canopy_fluxes_params.a_coef == 0.6
        @test CLM.canopy_fluxes_params.a_exp == 1.1
        @test CLM.canopy_fluxes_params.csoilc == 0.005
        @test CLM.canopy_fluxes_params.cv == 0.02
        @test CLM.canopy_fluxes_params.wind_min == 0.5

        # Reset to defaults
        CLM.canopy_fluxes_read_params!()
    end

    @testset "canopy_fluxes! single patch smoke test" begin
        # Initialize varpar for nlevgrnd
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

        np = 1   # 1 patch
        nc = 1   # 1 column
        ng = 1   # 1 gridcell
        nl = 1   # 1 landunit
        nlevgrnd = CLM.varpar.nlevgrnd
        nlevsno  = CLM.varpar.nlevsno
        nlevcan  = CLM.NLEVCAN
        nlevtot  = nlevsno + nlevgrnd

        # --- Set up all required data structures ---
        canopystate = CLM.CanopyStateData()
        CLM.canopystate_init!(canopystate, np)
        canopystate.elai_patch[1] = 2.0
        canopystate.esai_patch[1] = 0.5
        canopystate.laisun_patch[1] = 1.2
        canopystate.laisha_patch[1] = 0.8
        canopystate.displa_patch[1] = 5.0
        canopystate.htop_patch[1] = 10.0
        canopystate.frac_veg_nosno_patch[1] = 1
        canopystate.dleaf_patch[1] = 0.04
        canopystate.stem_biomass_patch[1] = 0.0
        canopystate.leaf_biomass_patch[1] = 0.0

        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nl, ng)
        energyflux.btran_patch[1] = 0.5
        energyflux.bsun_patch[1] = 0.5
        energyflux.bsha_patch[1] = 0.5
        energyflux.htvp_col[1] = CLM.HVAP
        energyflux.cgrnds_patch[1] = 0.0
        energyflux.cgrndl_patch[1] = 0.0
        energyflux.cgrnd_patch[1] = 0.0
        energyflux.dhsdt_canopy_patch[1] = 0.0
        energyflux.eflx_sh_stem_patch[1] = 0.0

        frictionvel = CLM.FrictionVelocityData()
        CLM.frictionvel_init!(frictionvel, np, nc)
        frictionvel.zetamaxstable = 0.5
        frictionvel.zsno = 0.00085
        frictionvel.zlnd = 0.000775
        frictionvel.z0mv_patch[1] = 0.5
        frictionvel.z0hv_patch[1] = 0.5
        frictionvel.z0qv_patch[1] = 0.5
        frictionvel.z0mg_col[1] = 0.01
        frictionvel.z0hg_col[1] = 0.01
        frictionvel.z0qg_col[1] = 0.01
        frictionvel.forc_hgt_u_patch[1] = 30.0
        frictionvel.forc_hgt_t_patch[1] = 30.0
        frictionvel.forc_hgt_q_patch[1] = 30.0
        frictionvel.ustar_patch[1] = 0.5
        frictionvel.um_patch[1] = 5.0
        frictionvel.uaf_patch[1] = 3.0
        frictionvel.taf_patch[1] = 290.0
        frictionvel.qaf_patch[1] = 0.008
        frictionvel.obu_patch[1] = -100.0
        frictionvel.zeta_patch[1] = -0.1
        frictionvel.vpd_patch[1] = 1.0
        frictionvel.rb1_patch[1] = 50.0
        frictionvel.ram1_patch[1] = 50.0
        frictionvel.num_iter_patch[1] = 0.0

        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)
        temperature.t_veg_patch[1] = 290.0
        temperature.t_stem_patch[1] = 289.0
        temperature.t_skin_patch[1] = 290.0
        temperature.thm_patch[1] = 290.0
        temperature.t_grnd_col[1] = 288.0
        temperature.t_h2osfc_col[1] = 288.0
        temperature.thv_col[1] = 291.0
        temperature.emv_patch[1] = 0.97
        temperature.emg_col[1] = 0.96
        temperature.t_ref2m_patch[1] = 290.0
        temperature.t_ref2m_r_patch[1] = 290.0
        # t_soisno_col: combined snow+soil indexing
        for j in 1:size(temperature.t_soisno_col, 2)
            temperature.t_soisno_col[1, j] = 288.0
        end

        solarabs = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs, np, nl)
        solarabs.sabv_patch[1] = 150.0

        soilstate = CLM.SoilStateData()
        CLM.soilstate_init!(soilstate, nc, nl)
        soilstate.soilbeta_col[1] = 0.8
        soilstate.soilresis_col[1] = 100.0
        for j in 1:nlevgrnd
            soilstate.watsat_col[1, j] = 0.45
        end

        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)
        waterfluxbulk.wf.qflx_tran_veg_patch[1] = 0.0
        waterfluxbulk.wf.qflx_evap_veg_patch[1] = 0.0
        waterfluxbulk.wf.qflx_evap_soi_patch[1] = 0.0

        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
        waterstatebulk.ws.liqcan_patch[1] = 0.1
        waterstatebulk.ws.snocan_patch[1] = 0.0

        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
        waterdiagbulk.fwet_patch[1] = 0.1
        waterdiagbulk.fdry_patch[1] = 0.8
        waterdiagbulk.frac_sno_eff_col[1] = 0.0
        waterdiagbulk.frac_h2osfc_col[1] = 0.0
        waterdiagbulk.snow_depth_col[1] = 0.0
        waterdiagbulk.qg_col[1] = 0.005
        waterdiagbulk.qg_snow_col[1] = 0.005
        waterdiagbulk.qg_soil_col[1] = 0.005
        waterdiagbulk.qg_h2osfc_col[1] = 0.005
        waterdiagbulk.dqgdT_col[1] = 0.0003
        waterdiagbulk.rh_af_patch[1] = 0.6
        waterdiagbulk.q_ref2m_patch[1] = 0.008
        waterdiagbulk.rh_ref2m_patch[1] = 50.0
        waterdiagbulk.rh_ref2m_r_patch[1] = 50.0
        waterdiagbulk.vpd_ref2m_patch[1] = 1000.0
        waterdiagbulk.iwue_ln_patch[1] = 0.0

        photosyns = CLM.PhotosynthesisData()
        CLM.photosynthesis_data_init!(photosyns, np)
        photosyns.rssun_patch[1] = 200.0
        photosyns.rssha_patch[1] = 300.0
        photosyns.fpsn_patch[1] = 5.0
        photosyns.psnsun_patch[1] = 3.0
        photosyns.psnsha_patch[1] = 2.0
        photosyns.gs_mol_sun_patch[1, 1] = 0.1
        photosyns.gs_mol_sha_patch[1, 1] = 0.05

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

        gridcell_data = CLM.GridcellData()
        CLM.gridcell_init!(gridcell_data, ng)

        mask_exposedveg = falses(np)
        mask_exposedveg[1] = true

        bounds_patch = 1:np
        bounds_col = 1:nc

        # Atmospheric forcing
        forc_lwrad_col = [350.0]
        forc_q_col     = [0.008]
        forc_pbot_col  = [101325.0]
        forc_th_col    = [290.0]
        forc_rho_col   = [1.2]
        forc_t_col     = [288.0]
        forc_u_grc     = [3.0]
        forc_v_grc     = [1.0]
        forc_pco2_grc  = [40.0]
        forc_po2_grc   = [21000.0]
        forc_hgt_t_grc = [10.0]
        forc_hgt_u_grc = [10.0]
        forc_hgt_q_grc = [10.0]
        dayl_grc       = [43200.0]
        max_dayl_grc   = [50000.0]
        downreg_patch  = [1.0]
        leafn_patch    = [1.0]
        dtime = 1800.0

        # Set up surface resistance method
        CLM.soil_resistance_read_nl!(soil_resis_method=CLM.SOIL_RESIS_LEEPIELKE_1992)

        # Reset control to defaults
        CLM.canopy_fluxes_read_nml!(
            use_undercanopy_stability=false,
            use_biomass_heat_storage=false,
            itmax_canopy_fluxes=40)
        CLM.canopy_fluxes_read_params!()

        # Call canopy_fluxes!
        CLM.canopy_fluxes!(
            canopystate, energyflux, frictionvel, temperature,
            solarabs, soilstate, waterfluxbulk, waterstatebulk,
            waterdiagbulk, photosyns, patch_data, col_data, gridcell_data,
            mask_exposedveg, bounds_patch, bounds_col,
            forc_lwrad_col, forc_q_col, forc_pbot_col, forc_th_col,
            forc_rho_col, forc_t_col, forc_u_grc, forc_v_grc,
            forc_pco2_grc, forc_po2_grc, forc_hgt_t_grc, forc_hgt_u_grc,
            forc_hgt_q_grc, dayl_grc, max_dayl_grc,
            downreg_patch, leafn_patch, dtime)

        # Basic sanity checks: outputs should be finite
        @test isfinite(temperature.t_veg_patch[1])
        @test isfinite(energyflux.eflx_sh_veg_patch[1])
        @test isfinite(energyflux.eflx_sh_grnd_patch[1])
        @test isfinite(energyflux.taux_patch[1])
        @test isfinite(energyflux.tauy_patch[1])
        @test isfinite(energyflux.dlrad_patch[1])
        @test isfinite(energyflux.ulrad_patch[1])
        @test isfinite(energyflux.cgrnds_patch[1])
        @test isfinite(energyflux.cgrndl_patch[1])
        @test isfinite(energyflux.cgrnd_patch[1])
        @test isfinite(temperature.t_ref2m_patch[1])
        @test isfinite(temperature.t_skin_patch[1])
        @test isfinite(frictionvel.ustar_patch[1])
        @test isfinite(frictionvel.obu_patch[1])

        # Physical reasonableness checks
        @test temperature.t_veg_patch[1] > 200.0  # not absurdly cold
        @test temperature.t_veg_patch[1] < 400.0  # not absurdly hot
        @test temperature.t_ref2m_patch[1] > 200.0
        @test temperature.t_ref2m_patch[1] < 400.0
        @test frictionvel.ustar_patch[1] > 0.0    # friction velocity positive
        @test energyflux.dlrad_patch[1] > 0.0     # downward longwave positive
        @test energyflux.ulrad_patch[1] > 0.0     # upward longwave positive

        # Reset surface resistance method
        CLM.soil_resistance_read_nl!()
    end
end
