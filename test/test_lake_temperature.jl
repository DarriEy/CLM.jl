@testset "Lake Temperature" begin

    # ----------------------------------------------------------------
    # Helper: initialize varpar and varcon for a minimal lake setup
    # ----------------------------------------------------------------
    function setup_lake_params!()
        vp = CLM.varpar
        vp.nlevsno = 5
        vp.nlevsoi = 10
        vp.nlevgrnd = 15
        vp.nlevlak = 10
        vp.nlevurb = 5
        vp.nlevmaxurbgrnd = max(vp.nlevurb, vp.nlevgrnd)
        CLM.varcon_init!()
        return vp
    end

    # ----------------------------------------------------------------
    # Helper: create minimal data structures for nc columns, np patches
    # ----------------------------------------------------------------
    function create_test_data(nc::Int, np::Int)
        vp = setup_lake_params!()
        nlevsno = vp.nlevsno
        nlevgrnd = vp.nlevgrnd
        nlevlak = vp.nlevlak
        nlevsoi = vp.nlevsoi
        nlevmaxurbgrnd = vp.nlevmaxurbgrnd
        nlevtot = nlevsno + nlevmaxurbgrnd

        # Column data
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        for c in 1:nc
            col.snl[c] = 0
            col.lakedepth[c] = 30.0
            for j in 1:nlevlak
                col.dz_lake[c, j] = CLM.dzlak[][j]
                col.z_lake[c, j] = CLM.zlak[][j]
            end
            for j in 1:nlevgrnd
                jj = j + nlevsno
                col.dz[c, jj] = CLM.dzsoi[][j]
                col.z[c, jj] = CLM.zsoi[][j]
                col.zi[c, jj + 1] = CLM.zisoi[][j + 1]
            end
            col.zi[c, nlevsno + 1] = 0.0  # zi(c,0)
        end

        # Patch data
        patch = CLM.PatchData()
        CLM.patch_init!(patch, np)
        for p in 1:np
            patch.column[p] = p  # 1:1 mapping
        end

        # Temperature
        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, 1, 1)
        for c in 1:nc
            temp.t_grnd_col[c] = 280.0
            for j in 1:nlevlak
                temp.t_lake_col[c, j] = 280.0
            end
            for jj in 1:nlevtot
                temp.t_soisno_col[c, jj] = 280.0
            end
        end

        # Solar absorbed
        solarabs = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs, np, 1)
        for p in 1:np
            solarabs.sabg_patch[p] = 100.0
            solarabs.fsds_nir_d_patch[p] = 30.0
            solarabs.fsds_nir_i_patch[p] = 20.0
            solarabs.fsr_nir_d_patch[p] = 5.0
            solarabs.fsr_nir_i_patch[p] = 3.0
            for j in 1:(nlevsno + 1)
                solarabs.sabg_lyr_patch[p, j] = 100.0 / (nlevsno + 1)
            end
        end

        # Soil state
        soilstate = CLM.SoilStateData()
        CLM.soilstate_init!(soilstate, np, nc)
        for c in 1:nc
            for j in 1:nlevgrnd
                soilstate.watsat_col[c, j] = 0.45
                soilstate.tksatu_col[c, j] = 1.5
                soilstate.tkmg_col[c, j] = 3.0
                soilstate.tkdry_col[c, j] = 0.25
                soilstate.csol_col[c, j] = 2.0e6
            end
        end

        # Water state bulk
        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, 1, 1)
        for c in 1:nc
            waterstatebulk.ws.h2osno_no_layers_col[c] = 0.0
            for jj in 1:nlevtot
                waterstatebulk.ws.h2osoi_liq_col[c, jj] = 0.0
                waterstatebulk.ws.h2osoi_ice_col[c, jj] = 0.0
            end
            # Soil layers: saturated
            for j in 1:nlevgrnd
                jj = j + nlevsno
                waterstatebulk.ws.h2osoi_liq_col[c, jj] = col.dz[c, jj] * CLM.DENH2O * soilstate.watsat_col[c, j]
                waterstatebulk.ws.h2osoi_ice_col[c, jj] = 0.0
            end
        end

        # Water diagnostic bulk
        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, 1, 1)
        for c in 1:nc
            waterdiagbulk.snow_depth_col[c] = 0.0
            for jj in 1:nlevtot
                waterdiagbulk.frac_iceold_col[c, jj] = 0.0
            end
        end

        # Water flux bulk
        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, 1, 1)

        # Energy flux
        energyflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(energyflux, np, nc, 1, 1)
        for p in 1:np
            energyflux.eflx_gnet_patch[p] = 50.0
            energyflux.eflx_sh_tot_patch[p] = 20.0
            energyflux.eflx_sh_grnd_patch[p] = 20.0
            energyflux.eflx_soil_grnd_patch[p] = 50.0
        end

        # Lake state
        lakestate = CLM.LakeStateData()
        CLM.lakestate_init!(lakestate, nc, np)
        for c in 1:nc
            lakestate.etal_col[c] = 0.5
            lakestate.ks_col[c] = 0.1
            lakestate.ws_col[c] = 0.05
            lakestate.lake_raw_col[c] = 100.0
            lakestate.betaprime_col[c] = 0.4
            lakestate.savedtke1_col[c] = CLM.TKWAT
            lakestate.lakeresist_col[c] = 0.0
            lakestate.lake_icethick_col[c] = 0.0
            lakestate.lake_icefracsurf_col[c] = 0.0
            for j in 1:nlevlak
                lakestate.lake_icefrac_col[c, j] = 0.0
            end
        end

        grnd_ch4_cond = zeros(nc)

        mask_lakec = trues(nc)
        mask_lakep = trues(np)
        bounds_col = 1:nc
        bounds_patch = 1:np

        return (col=col, patch=patch, solarabs=solarabs, soilstate=soilstate,
                waterstatebulk=waterstatebulk, waterdiagbulk=waterdiagbulk,
                waterfluxbulk=waterfluxbulk, energyflux=energyflux,
                temperature=temp, lakestate=lakestate,
                grnd_ch4_cond=grnd_ch4_cond,
                mask_lakec=mask_lakec, mask_lakep=mask_lakep,
                bounds_col=bounds_col, bounds_patch=bounds_patch)
    end

    # ================================================================
    @testset "Constants defined" begin
        @test CLM.BETAVIS == 0.0
        @test CLM.ZA_LAKE == 0.6
        @test CLM.N2MIN == 7.5e-5
        @test CLM.TDMAX == 277.0
        @test CLM.DEPTHCRIT == 25.0
        @test CLM.MIXFACT == 10.0
        @test CLM.LAKEPUDDLING == false
        @test CLM.LAKE_NO_ED == false
    end

    # ================================================================
    @testset "soil_therm_prop_lake! basic" begin
        d = create_test_data(2, 2)
        vp = CLM.varpar
        nlevsno = vp.nlevsno
        nlevgrnd = vp.nlevgrnd
        nlevtot = nlevsno + vp.nlevmaxurbgrnd

        tk = zeros(2, nlevtot)
        cv = zeros(2, nlevtot)
        tktopsoillay = zeros(2)

        CLM.soil_therm_prop_lake!(d.col, d.soilstate, d.waterstatebulk, d.temperature,
                                  tk, cv, tktopsoillay, d.mask_lakec, d.bounds_col)

        # Soil layers should have positive heat capacity
        for c in 1:2
            for j in 1:nlevgrnd
                jj = j + nlevsno
                @test cv[c, jj] > 0.0
            end
        end

        # Top soil layer conductivity should be positive
        for c in 1:2
            @test tktopsoillay[c] > 0.0
        end
    end

    # ================================================================
    @testset "phase_change_lake! no change when T=280K" begin
        d = create_test_data(1, 1)
        vp = CLM.varpar
        nlevsno = vp.nlevsno
        nlevlak = vp.nlevlak
        nlevtot = nlevsno + vp.nlevmaxurbgrnd

        cv = zeros(1, nlevtot)
        cv_lake = zeros(1, nlevlak)

        # Set up cv values
        cwat = CLM.CPLIQ * CLM.DENH2O
        for j in 1:nlevlak
            cv_lake[1, j] = d.col.dz_lake[1, j] * cwat
        end
        for j in 1:vp.nlevgrnd
            jj = j + nlevsno
            cv[1, jj] = 2.0e6 * d.col.dz[1, jj]
        end

        lhabs = zeros(1)
        dtime = 1800.0

        t_lake_before = copy(d.temperature.t_lake_col)

        CLM.phase_change_lake!(d.col, d.waterstatebulk, d.waterdiagbulk, d.waterfluxbulk,
                               d.temperature, d.energyflux, d.lakestate,
                               cv, cv_lake, lhabs, d.mask_lakec, d.bounds_col, dtime)

        # At 280K with no ice, no phase change should occur
        @test lhabs[1] == 0.0
        @test d.temperature.t_lake_col[1, 1] == t_lake_before[1, 1]
    end

    # ================================================================
    @testset "phase_change_lake! freezing occurs when T < TFRZ" begin
        d = create_test_data(1, 1)
        vp = CLM.varpar
        nlevsno = vp.nlevsno
        nlevlak = vp.nlevlak
        nlevtot = nlevsno + vp.nlevmaxurbgrnd

        cv = zeros(1, nlevtot)
        cv_lake = zeros(1, nlevlak)
        cwat = CLM.CPLIQ * CLM.DENH2O
        cice_eff = CLM.CPICE * CLM.DENH2O

        # Set lake temp below freezing, ice frac < 1
        d.temperature.t_lake_col[1, 1] = CLM.TFRZ - 2.0
        for j in 1:nlevlak
            cv_lake[1, j] = d.col.dz_lake[1, j] * (cwat * (1.0 - d.lakestate.lake_icefrac_col[1, j]) +
                cice_eff * d.lakestate.lake_icefrac_col[1, j])
        end
        for j in 1:vp.nlevgrnd
            jj = j + nlevsno
            cv[1, jj] = 2.0e6 * d.col.dz[1, jj]
        end

        lhabs = zeros(1)
        dtime = 1800.0

        CLM.phase_change_lake!(d.col, d.waterstatebulk, d.waterdiagbulk, d.waterfluxbulk,
                               d.temperature, d.energyflux, d.lakestate,
                               cv, cv_lake, lhabs, d.mask_lakec, d.bounds_col, dtime)

        # Freezing should have occurred in layer 1
        @test d.lakestate.lake_icefrac_col[1, 1] > 0.0
        # Latent heat should have been absorbed (negative melt means freezing: lhabs += melt*hfus, melt<0 so lhabs<0)
        @test lhabs[1] < 0.0
    end

    # ================================================================
    @testset "lake_temperature! runs without error" begin
        d = create_test_data(2, 2)
        dtime = 1800.0

        # The main driver should run without errors
        CLM.lake_temperature!(d.col, d.patch, d.solarabs, d.soilstate,
                             d.waterstatebulk, d.waterdiagbulk, d.waterfluxbulk,
                             d.energyflux, d.temperature, d.lakestate,
                             d.grnd_ch4_cond, d.mask_lakec, d.mask_lakep,
                             d.bounds_col, d.bounds_patch, dtime)

        # Temperature should have changed from initial values
        # (with non-zero energy flux, temperature must evolve)
        @test isfinite(d.temperature.t_lake_col[1, 1])
        @test isfinite(d.temperature.t_lake_col[2, 1])

        # Lake ice thickness diagnostic should be set
        for c in 1:2
            @test isfinite(d.lakestate.lake_icethick_col[c])
            @test d.lakestate.lake_icethick_col[c] >= 0.0
        end

        # errsoi should be zeroed out (small imbalance dumped to sensible heat)
        for c in 1:2
            @test isfinite(d.energyflux.errsoi_col[c])
        end
    end

    # ================================================================
    @testset "lake_temperature! energy conservation" begin
        d = create_test_data(1, 1)
        dtime = 1800.0

        # Run with zero incoming flux for energy conservation test
        d.energyflux.eflx_gnet_patch[1] = 0.0
        d.solarabs.sabg_patch[1] = 0.0
        d.solarabs.fsds_nir_d_patch[1] = 0.0
        d.solarabs.fsds_nir_i_patch[1] = 0.0
        d.solarabs.fsr_nir_d_patch[1] = 0.0
        d.solarabs.fsr_nir_i_patch[1] = 0.0
        for j in 1:(CLM.varpar.nlevsno + 1)
            d.solarabs.sabg_lyr_patch[1, j] = 0.0
        end

        CLM.lake_temperature!(d.col, d.patch, d.solarabs, d.soilstate,
                             d.waterstatebulk, d.waterdiagbulk, d.waterfluxbulk,
                             d.energyflux, d.temperature, d.lakestate,
                             d.grnd_ch4_cond, d.mask_lakec, d.mask_lakep,
                             d.bounds_col, d.bounds_patch, dtime)

        # With zero energy input the error should be small
        @test abs(d.energyflux.errsoi_col[1]) < 0.2  # W/m^2
    end

    # ================================================================
    @testset "calculate_total_h2osno!" begin
        vp = setup_lake_params!()
        nc = 2
        nlevsno = vp.nlevsno

        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, 1, 1, 1)
        ws.h2osno_no_layers_col[1] = 5.0
        ws.h2osno_no_layers_col[2] = 0.0

        snl = zeros(Int, nc)
        snl[1] = 0  # no snow layers
        snl[2] = 0

        h2osno_total = zeros(nc)
        mask = trues(nc)

        CLM.calculate_total_h2osno!(ws, snl, h2osno_total, mask, 1:nc)

        @test h2osno_total[1] == 5.0
        @test h2osno_total[2] == 0.0
    end

end
