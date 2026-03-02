@testset "Urban Fluxes" begin

    @testset "UrbanFluxesParamsData" begin
        ufp = CLM.UrbanFluxesParamsData()
        @test ufp.wind_min == 0.0
    end

    @testset "urban_fluxes_read_params!" begin
        CLM.urban_fluxes_read_params!(wind_min=1.0)
        @test CLM.urban_fluxes_params.wind_min == 1.0
        CLM.urban_fluxes_read_params!(wind_min=0.0)
    end

    @testset "simple_wasteheatfromac" begin
        # Save and set urban control
        old_hac = CLM.urban_ctrl.urban_hac
        old_read = CLM.urban_ctrl.read_namelist
        CLM.urban_ctrl.read_namelist = true

        # Test with wasteheat off
        CLM.urban_ctrl.urban_hac = CLM.URBAN_HAC_OFF
        (wh, hfac) = CLM.simple_wasteheatfromac(100.0, 50.0)
        @test wh == 0.0
        @test hfac == 0.0

        # Test with HAC on (no wasteheat)
        CLM.urban_ctrl.urban_hac = CLM.URBAN_HAC_ON
        (wh, hfac) = CLM.simple_wasteheatfromac(100.0, 50.0)
        @test wh == 0.0
        @test hfac == 100.0  # abs(eflx_urban_ac)

        # Test with wasteheat on
        CLM.urban_ctrl.urban_hac = CLM.URBAN_WASTEHEAT_ON
        (wh, hfac) = CLM.simple_wasteheatfromac(100.0, 50.0)
        @test wh ≈ CLM.AC_WASTEHEAT_FACTOR * 100.0 + CLM.HT_WASTEHEAT_FACTOR * 50.0
        @test hfac == 100.0  # abs(eflx_urban_ac)

        # Restore
        CLM.urban_ctrl.urban_hac = old_hac
        CLM.urban_ctrl.read_namelist = old_read
    end

    @testset "urban_fluxes! single landunit smoke test" begin
        # Initialize varpar
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        nlevgrnd = CLM.varpar.nlevgrnd
        nlevsno  = CLM.varpar.nlevsno
        nlevurb  = CLM.varpar.nlevurb

        # Configuration: 1 gridcell, 1 urban landunit, 5 urban columns, 5 patches
        ng = 1; nl = 1; nc = 5; np = 5

        # Set urban fluxes params
        CLM.urban_fluxes_read_params!(wind_min=1.0)

        # Set up urban control for simple building temp
        old_hac = CLM.urban_ctrl.urban_hac
        old_read = CLM.urban_ctrl.read_namelist
        old_btm = CLM.urban_ctrl.building_temp_method
        CLM.urban_ctrl.read_namelist = true
        CLM.urban_ctrl.building_temp_method = CLM.BUILDING_TEMP_METHOD_SIMPLE
        CLM.urban_ctrl.urban_hac = CLM.URBAN_HAC_OFF

        # --- Set up data structures ---
        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, nl, ng)
        for c in 1:nc
            energyflux.htvp_col[c] = CLM.HVAP
            energyflux.eflx_urban_ac_col[c] = 0.0
            energyflux.eflx_urban_heat_col[c] = 0.0
        end
        energyflux.eflx_urban_ac_lun[1] = 0.0
        energyflux.eflx_urban_heat_lun[1] = 0.0

        frictionvel = CLM.FrictionVelocityData()
        CLM.frictionvel_init!(frictionvel, np, nc)
        frictionvel.zetamaxstable = 0.5
        for p in 1:np
            frictionvel.forc_hgt_u_patch[p] = 30.0
            frictionvel.forc_hgt_t_patch[p] = 30.0
            frictionvel.forc_hgt_q_patch[p] = 30.0
            frictionvel.u10_clm_patch[p] = 5.0
        end

        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)
        temperature.taf_lun[1] = 290.0
        for c in 1:nc
            temperature.t_grnd_col[c] = 290.0
            for j in 1:size(temperature.t_soisno_col, 2)
                temperature.t_soisno_col[c, j] = 290.0
            end
        end

        soilstate = CLM.SoilStateData()
        CLM.soilstate_init!(soilstate, np, nc)
        soilstate.soilalpha_u_col[1:nc] .= 0.5
        for c in 1:nc
            for j in 1:nlevgrnd
                soilstate.rootr_road_perv_col[c, j] = (j <= 5 ? 0.2 : 0.0)
            end
        end

        urbanparams = CLM.UrbanParamsData()
        CLM.urbanparams_init!(urbanparams, nl; nlevurb=nlevurb, numrad=CLM.NUMRAD)
        urbanparams.wind_hgt_canyon[1] = 5.0
        urbanparams.eflx_traffic_factor[1] = 0.0

        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)

        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)
        for c in 1:nc
            waterstatebulk.ws.h2osoi_liq_col[c, 1] = 0.5
            waterstatebulk.ws.h2osoi_ice_col[c, 1] = 0.0
        end

        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)
        waterdiagbulk.qaf_lun[1] = 0.005
        for c in 1:nc
            waterdiagbulk.qg_col[c] = 0.005
            waterdiagbulk.snow_depth_col[c] = 0.0
            waterdiagbulk.frac_sno_col[c] = 0.0
            waterdiagbulk.dqgdT_col[c] = 0.0003
        end

        # --- Set up hierarchy ---
        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)
        # 5 patches, one per urban column
        for p in 1:np
            patch_data.column[p] = p
            patch_data.gridcell[p] = 1
            patch_data.landunit[p] = 1
            patch_data.active[p] = true
        end

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)
        # Column types: roof, sunwall, shadewall, road_imperv, road_perv
        col_data.itype[1] = CLM.ICOL_ROOF
        col_data.itype[2] = CLM.ICOL_SUNWALL
        col_data.itype[3] = CLM.ICOL_SHADEWALL
        col_data.itype[4] = CLM.ICOL_ROAD_IMPERV
        col_data.itype[5] = CLM.ICOL_ROAD_PERV
        for c in 1:nc
            col_data.landunit[c] = 1
            col_data.snl[c] = 0
        end

        lun_data = CLM.LandunitData()
        CLM.landunit_init!(lun_data, nl)
        lun_data.gridcell[1] = 1
        lun_data.itype[1] = CLM.ISTURB_MIN
        lun_data.urbpoi[1] = true
        lun_data.active[1] = true
        lun_data.patchi[1] = 1
        lun_data.patchf[1] = np
        lun_data.coli[1] = 1
        lun_data.colf[1] = nc
        lun_data.canyon_hwr[1] = 1.0
        lun_data.wtroad_perv[1] = 0.4
        lun_data.wtlunit_roof[1] = 0.5
        lun_data.ht_roof[1] = 10.0
        lun_data.z_0_town[1] = 0.5
        lun_data.z_d_town[1] = 5.0

        # Filters
        filter_nourbanl = Int[]
        filter_urbanl   = [1]
        filter_urbanc   = [1, 2, 3, 4, 5]
        filter_urbanp   = [1, 2, 3, 4, 5]

        # Atmospheric forcing
        forc_t_grc   = [288.0]
        forc_th_grc  = [290.0]
        forc_rho_grc = [1.2]
        forc_q_grc   = [0.008]
        forc_pbot_grc = [101325.0]
        forc_u_grc   = [3.0]
        forc_v_grc   = [1.0]

        # Call urban_fluxes!
        CLM.urban_fluxes!(
            energyflux, frictionvel, temperature, soilstate,
            urbanparams, waterfluxbulk, waterstatebulk, waterdiagbulk,
            patch_data, col_data, lun_data,
            0, filter_nourbanl,
            1, filter_urbanl,
            5, filter_urbanc,
            5, filter_urbanp,
            1:nl, 1:nc, 1:np,
            forc_t_grc, forc_th_grc, forc_rho_grc, forc_q_grc,
            forc_pbot_grc, forc_u_grc, forc_v_grc;
            dtime=3600.0, nstep=1)

        # --- Finiteness checks ---
        for p in 1:np
            @test isfinite(energyflux.eflx_sh_grnd_patch[p])
            @test isfinite(energyflux.eflx_sh_tot_patch[p])
            @test isfinite(energyflux.eflx_sh_tot_u_patch[p])
            @test isfinite(energyflux.taux_patch[p])
            @test isfinite(energyflux.tauy_patch[p])
            @test isfinite(energyflux.cgrnds_patch[p])
            @test isfinite(energyflux.cgrndl_patch[p])
            @test isfinite(energyflux.cgrnd_patch[p])
            @test isfinite(temperature.t_ref2m_patch[p])
            @test isfinite(waterdiagbulk.q_ref2m_patch[p])
            @test isfinite(waterdiagbulk.rh_ref2m_patch[p])
            @test isfinite(frictionvel.ram1_patch[p])
        end

        # --- Physical reasonableness checks ---
        # Urban canopy air temperature should be near ambient (200-400 K range)
        @test temperature.taf_lun[1] > 200.0
        @test temperature.taf_lun[1] < 400.0

        # 2m temperatures should be reasonable
        for p in 1:np
            @test temperature.t_ref2m_patch[p] > 200.0
            @test temperature.t_ref2m_patch[p] < 400.0
        end

        # Wind stress should oppose wind direction
        for p in 1:np
            @test energyflux.taux_patch[p] < 0.0
            @test energyflux.tauy_patch[p] < 0.0
        end

        # Relative humidity in [0, 100]
        for p in 1:np
            @test waterdiagbulk.rh_ref2m_patch[p] >= 0.0
            @test waterdiagbulk.rh_ref2m_patch[p] <= 100.0
        end

        # Aerodynamic resistance should be positive
        for p in 1:np
            @test frictionvel.ram1_patch[p] > 0.0
        end

        # dlrad and ulrad should be zero (no canopy for urban)
        for p in 1:np
            @test energyflux.dlrad_patch[p] == 0.0
            @test energyflux.ulrad_patch[p] == 0.0
        end

        # eflx_sh_snow/soil/h2osfc should be zero for urban
        for p in 1:np
            @test energyflux.eflx_sh_snow_patch[p] == 0.0
            @test energyflux.eflx_sh_soil_patch[p] == 0.0
            @test energyflux.eflx_sh_h2osfc_patch[p] == 0.0
        end

        # Total SH should equal ground SH for urban
        for p in 1:np
            @test energyflux.eflx_sh_tot_patch[p] == energyflux.eflx_sh_grnd_patch[p]
        end

        # Wall evaporation should be zero
        # patch 2 = sunwall, patch 3 = shadewall
        @test waterfluxbulk.wf.qflx_evap_soi_patch[2] == 0.0
        @test waterfluxbulk.wf.qflx_evap_soi_patch[3] == 0.0

        # Roots should only be nonzero for pervious road
        for j in 1:min(5, nlevgrnd)
            # patch 5 is pervious road
            @test soilstate.rootr_patch[5, j] ≈ soilstate.rootr_road_perv_col[5, j]
            # other patches should have zero roots
            @test soilstate.rootr_patch[1, j] == 0.0
            @test soilstate.rootr_patch[2, j] == 0.0
            @test soilstate.rootr_patch[3, j] == 0.0
            @test soilstate.rootr_patch[4, j] == 0.0
        end

        # Building temperature should be finite (simple method)
        @test isfinite(temperature.t_building_lun[1])
        @test temperature.t_building_lun[1] > 200.0
        @test temperature.t_building_lun[1] < 400.0

        # Derivative cgrnds should be positive for roof/road/wall
        for p in 1:np
            @test energyflux.cgrnds_patch[p] > 0.0
        end

        # Restore
        CLM.urban_ctrl.urban_hac = old_hac
        CLM.urban_ctrl.read_namelist = old_read
        CLM.urban_ctrl.building_temp_method = old_btm
        CLM.urban_fluxes_read_params!(wind_min=0.0)
    end

    @testset "calc_simple_internal_building_temp!" begin
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        nlevurb = CLM.varpar.nlevurb

        old_read = CLM.urban_ctrl.read_namelist
        old_btm  = CLM.urban_ctrl.building_temp_method
        CLM.urban_ctrl.read_namelist = true
        CLM.urban_ctrl.building_temp_method = CLM.BUILDING_TEMP_METHOD_SIMPLE

        nl = 1; nc = 3  # roof, sunwall, shadewall

        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, 3, nc, nl, 1)
        # Set inner wall/roof temperatures in t_soisno_col at nlevurb
        temperature.t_soisno_col[1, nlevurb] = 295.0   # roof inner
        temperature.t_soisno_col[2, nlevurb] = 293.0   # sunwall inner
        temperature.t_soisno_col[3, nlevurb] = 291.0   # shadewall inner

        col_data = CLM.ColumnData()
        CLM.column_init!(col_data, nc)
        col_data.itype[1] = CLM.ICOL_ROOF
        col_data.itype[2] = CLM.ICOL_SUNWALL
        col_data.itype[3] = CLM.ICOL_SHADEWALL
        col_data.landunit[1] = 1
        col_data.landunit[2] = 1
        col_data.landunit[3] = 1

        lun_data = CLM.LandunitData()
        CLM.landunit_init!(lun_data, nl)
        lun_data.ht_roof[1] = 10.0
        lun_data.canyon_hwr[1] = 1.0
        lun_data.wtlunit_roof[1] = 0.5

        CLM.calc_simple_internal_building_temp!(
            temperature, col_data, lun_data,
            3, [1, 2, 3], 1, [1])

        @test isfinite(temperature.t_building_lun[1])
        # Building temp should be a weighted average of inner temps
        @test temperature.t_building_lun[1] > 290.0
        @test temperature.t_building_lun[1] < 296.0

        CLM.urban_ctrl.read_namelist = old_read
        CLM.urban_ctrl.building_temp_method = old_btm
    end

end
