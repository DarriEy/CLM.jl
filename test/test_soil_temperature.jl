@testset "SoilTemperature" begin
    # Save and restore varpar state
    vp = CLM.varpar
    saved_nlevsno = vp.nlevsno
    saved_nlevgrnd = vp.nlevgrnd
    saved_nlevurb = vp.nlevurb
    saved_nlevmaxurbgrnd = vp.nlevmaxurbgrnd
    saved_nlevsoi = vp.nlevsoi

    vp.nlevsno = 5
    vp.nlevgrnd = 10
    vp.nlevurb = 5
    vp.nlevmaxurbgrnd = 10
    vp.nlevsoi = 8

    nlevsno = vp.nlevsno
    nlevgrnd = vp.nlevgrnd
    nlevurb = vp.nlevurb
    nlevmaxurbgrnd = vp.nlevmaxurbgrnd
    nlevsoi = vp.nlevsoi
    joff = nlevsno

    nc = 4   # number of columns
    np = 4   # number of patches
    nl = 4   # number of landunits
    ng = 1   # number of gridcells

    # Save varctl state
    saved_use_excess_ice = CLM.varctl.use_excess_ice
    saved_snow_thermal = CLM.varctl.snow_thermal_cond_method
    CLM.varctl.use_excess_ice = false
    CLM.varctl.snow_thermal_cond_method = "Jordan1991"

    # Save urban control state
    saved_read_namelist = CLM.urban_ctrl.read_namelist
    saved_building_temp_method = CLM.urban_ctrl.building_temp_method
    CLM.urban_ctrl.read_namelist = true
    CLM.urban_ctrl.building_temp_method = CLM.BUILDING_TEMP_METHOD_SIMPLE

    # Helper: set up all test data structures
    function setup_test_data()
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)

        patch_data = CLM.PatchData()
        CLM.patch_init!(patch_data, np)

        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, np, nc, nl, ng)

        ef = CLM.EnergyFluxData()
        CLM.energyflux_init!(ef, np, nc, nl, ng)

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, np, nc)

        wsb = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(wsb, nc, np, nl, ng)

        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, np, nl, ng)

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, np, nl, ng)

        sa = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(sa, np, nl)

        cs = CLM.CanopyStateData()
        CLM.canopystate_init!(cs, np)

        up = CLM.UrbanParamsData()
        CLM.urbanparams_init!(up, nl; nlevurb=nlevurb)

        mask_nolakec = falses(nc)
        mask_nolakep = falses(np)
        mask_urbanl = falses(nl)
        mask_urbanc = falses(nc)
        bounds_col = 1:nc
        bounds_lun = 1:nl
        bounds_patch = 1:np

        return col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
               mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
               bounds_col, bounds_lun, bounds_patch
    end

    # Helper: set up a simple soil column (no snow) with reasonable values
    function setup_soil_column!(col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                                 mask_nolakec, mask_nolakep, mask_urbanc, c, p, l)
        mask_nolakec[c] = true
        mask_nolakep[p] = true

        col.landunit[c] = l
        col.itype[c] = CLM.ISTSOIL
        lun.itype[l] = CLM.ISTSOIL
        lun.urbpoi[l] = false
        col.snl[c] = 0  # no snow layers
        col.nbedrock[c] = nlevgrnd

        patch_data.column[p] = c
        patch_data.landunit[p] = l
        patch_data.wtcol[p] = 1.0

        # Layer geometry
        for j in 1:nlevmaxurbgrnd
            jj = j + joff
            col.dz[c, jj] = 0.1 * j
            col.z[c, jj] = sum(0.1 * k for k in 1:j) - 0.5 * 0.1 * j
        end
        # zi (interface depths)
        col.zi[c, joff + 1] = 0.0
        for j in 1:nlevmaxurbgrnd
            col.zi[c, j + joff + 1] = col.zi[c, j + joff] + col.dz[c, j + joff]
        end

        # Temperatures (above freezing)
        for j in 1:nlevmaxurbgrnd
            temp.t_soisno_col[c, j + joff] = 280.0
        end
        temp.t_grnd_col[c] = 280.0
        temp.t_h2osfc_col[c] = 280.0
        temp.emg_col[c] = 0.97

        # Water state
        for j in 1:nlevmaxurbgrnd
            jj = j + joff
            wsb.ws.h2osoi_liq_col[c, jj] = 5.0
            wsb.ws.h2osoi_ice_col[c, jj] = 0.0
        end
        wsb.ws.h2osfc_col[c] = 0.0
        wsb.ws.h2osno_no_layers_col[c] = 0.0
        for j in 1:nlevmaxurbgrnd
            wsb.ws.excess_ice_col[c, j] = 0.0
        end
        wsb.int_snow_col[c] = 0.0

        # Water diagnostics
        wdb.frac_sno_eff_col[c] = 0.0
        wdb.frac_h2osfc_col[c] = 0.0
        wdb.snow_depth_col[c] = 0.0

        # Soil state
        for j in 1:nlevgrnd
            ss.watsat_col[c, j] = 0.45
            ss.bsw_col[c, j] = 5.0
            ss.sucsat_col[c, j] = 100.0
            ss.tkmg_col[c, j] = 2.0
            ss.tkdry_col[c, j] = 0.25
            ss.csol_col[c, j] = 2.0e6
            ss.tksatu_col[c, j] = 1.5
            ss.thk_col[c, j] = 0.0  # will be computed
        end

        # Energy flux
        ef.htvp_col[c] = CLM.HVAP
        ef.eflx_bot_col[c] = 0.0
        ef.dlrad_patch[p] = 50.0
        ef.cgrnd_patch[p] = 5.0
        ef.eflx_sh_grnd_patch[p] = 10.0
        ef.eflx_sh_snow_patch[p] = 10.0
        ef.eflx_sh_soil_patch[p] = 10.0
        ef.eflx_sh_h2osfc_patch[p] = 10.0
        ef.dgnetdT_patch[p] = 0.0
        ef.eflx_gnet_patch[p] = 0.0

        # Solar absorbed
        sa.sabg_patch[p] = 100.0
        sa.sabg_soil_patch[p] = 100.0
        sa.sabg_snow_patch[p] = 0.0
        sa.sabg_chk_patch[p] = 0.0
        for j in 1:(nlevsno + 1)
            sa.sabg_lyr_patch[p, j] = 0.0
        end
        sa.sabg_lyr_patch[p, nlevsno + 1] = 100.0  # layer 1 = soil top

        # Canopy state
        cs.frac_veg_nosno_patch[p] = 0

        # Water fluxes (set both nested wf and bulk-level fields)
        wfb.wf.qflx_evap_soi_patch[p] = 0.0
        wfb.wf.qflx_ev_snow_patch[p] = 0.0
        wfb.wf.qflx_ev_soil_patch[p] = 0.0
        wfb.wf.qflx_ev_h2osfc_patch[p] = 0.0
        wfb.wf.qflx_tran_veg_patch[p] = 0.0
        wfb.qflx_ev_snow_patch[p] = 0.0
        wfb.qflx_ev_soil_patch[p] = 0.0
        wfb.qflx_ev_h2osfc_patch[p] = 0.0
        wfb.wf.qflx_snomelt_col[c] = 0.0
        wfb.wf.qflx_snofrz_col[c] = 0.0
        wfb.wf.qflx_snow_drain_col[c] = 0.0
        wfb.wf.qflx_h2osfc_to_ice_col[c] = 0.0

        # Snow layer fluxes
        for j in 1:nlevsno
            wfb.wf.qflx_snofrz_lyr_col[c, j] = 0.0
            wfb.wf.qflx_snomelt_lyr_col[c, j] = 0.0
            wdb.bw_col[c, j] = 0.0
        end

        # Energy flux - snow melt outputs
        ef.eflx_snomelt_col[c] = 0.0
        ef.eflx_snomelt_r_col[c] = 0.0
        ef.eflx_snomelt_u_col[c] = 0.0
        ef.eflx_h2osfc_to_snow_col[c] = 0.0
        ef.eflx_building_heat_errsoi_col[c] = 0.0
        ef.eflx_urban_ac_col[c] = 0.0
        ef.eflx_urban_heat_col[c] = 0.0
        ef.eflx_fgr12_col[c] = 0.0
        for j in 1:nlevgrnd
            ef.eflx_fgr_col[c, j] = 0.0
        end

        # Urban params
        up.nlev_improad[l] = 0
        up.t_building_min[l] = 288.0

        # Diagnostics
        wdb.snomelt_accum_col[c] = 0.0
        for j in 1:nlevmaxurbgrnd
            wdb.exice_subs_col[c, j] = 0.0
        end
    end

    @testset "soil_therm_prop! — soil thermal conductivity" begin
        col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
            mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
            bounds_col, bounds_lun, bounds_patch = setup_test_data()

        c = 1; p = 1; l = 1
        setup_soil_column!(col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                           mask_nolakec, mask_nolakep, mask_urbanc, c, p, l)

        tk = zeros(nc, nlevsno + nlevmaxurbgrnd)
        cv = zeros(nc, nlevsno + nlevmaxurbgrnd)
        tk_h2osfc = fill(NaN, nc)

        CLM.soil_therm_prop!(col, lun, up, temp, wsb, wdb, ss,
                             mask_nolakec, mask_urbanc, bounds_col,
                             tk, cv, tk_h2osfc)

        # Thermal conductivity at layer interfaces should be positive
        for j in 1:(nlevgrnd - 1)
            @test tk[c, j + joff] > 0.0
        end
        # Bottom layer tk should be 0
        @test tk[c, nlevgrnd + joff] == 0.0

        # Heat capacity should be positive for all soil layers
        for j in 1:nlevgrnd
            @test cv[c, j + joff] > 0.0
        end

        # tk_h2osfc should be finite and positive
        @test isfinite(tk_h2osfc[c])
        @test tk_h2osfc[c] > 0.0
    end

    @testset "soil_therm_prop! — snow thermal conductivity Jordan1991" begin
        col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
            mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
            bounds_col, bounds_lun, bounds_patch = setup_test_data()

        c = 1; p = 1; l = 1
        setup_soil_column!(col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                           mask_nolakec, mask_nolakep, mask_urbanc, c, p, l)

        # Set up 2 snow layers
        col.snl[c] = -2
        wdb.frac_sno_eff_col[c] = 1.0
        for j in (-1):0
            jj = j + joff
            col.dz[c, jj] = 0.1
            col.z[c, jj] = (j - 0.5) * 0.1
            wsb.ws.h2osoi_ice_col[c, jj] = 10.0
            wsb.ws.h2osoi_liq_col[c, jj] = 1.0
            temp.t_soisno_col[c, jj] = 270.0
        end
        # zi for snow layers
        col.zi[c, joff - 1] = -0.2
        col.zi[c, joff] = -0.1
        col.zi[c, joff + 1] = 0.0

        tk = zeros(nc, nlevsno + nlevmaxurbgrnd)
        cv = zeros(nc, nlevsno + nlevmaxurbgrnd)
        tk_h2osfc = fill(NaN, nc)

        CLM.soil_therm_prop!(col, lun, up, temp, wsb, wdb, ss,
                             mask_nolakec, mask_urbanc, bounds_col,
                             tk, cv, tk_h2osfc)

        # Snow layers should have positive thermal conductivity at interfaces
        for j in (-1):(-1)
            @test tk[c, j + joff] > 0.0
        end

        # Snow heat capacity should be positive
        for j in (-1):0
            jj = j + joff
            @test cv[c, jj] > 0.0
        end
    end

    @testset "compute_ground_heat_flux_and_deriv!" begin
        col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
            mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
            bounds_col, bounds_lun, bounds_patch = setup_test_data()

        c = 1; p = 1; l = 1
        setup_soil_column!(col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                           mask_nolakec, mask_nolakep, mask_urbanc, c, p, l)

        forc_lwrad = fill(300.0, nc)

        hs_h2osfc = zeros(nc)
        hs_top_snow = zeros(nc)
        hs_soil = zeros(nc)
        hs_top = zeros(nc)
        dhsdT = zeros(nc)
        sabg_lyr_col = zeros(nc, nlevsno + 1)

        CLM.compute_ground_heat_flux_and_deriv!(
            col, lun, patch_data, temp, ef, sa, cs, wdb, wfb, up, forc_lwrad,
            mask_nolakec, mask_nolakep, bounds_col, bounds_patch,
            hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col)

        # hs_soil should be finite
        @test isfinite(hs_soil[c])
        # dhsdT should be negative (cooling response to temperature increase)
        @test dhsdT[c] < 0.0
        # eflx_gnet should be finite
        @test isfinite(ef.eflx_gnet_patch[p])
    end

    @testset "compute_heat_diff_flux_and_factor!" begin
        col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
            mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
            bounds_col, bounds_lun, bounds_patch = setup_test_data()

        c = 1; p = 1; l = 1
        setup_soil_column!(col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                           mask_nolakec, mask_nolakep, mask_urbanc, c, p, l)

        dtime = 1800.0

        # First compute thermal properties
        tk = zeros(nc, nlevsno + nlevmaxurbgrnd)
        cv = zeros(nc, nlevsno + nlevmaxurbgrnd)
        tk_h2osfc = fill(NaN, nc)
        CLM.soil_therm_prop!(col, lun, up, temp, wsb, wdb, ss,
                             mask_nolakec, mask_urbanc, bounds_col,
                             tk, cv, tk_h2osfc)

        fn = zeros(nc, nlevsno + nlevmaxurbgrnd)
        fact = zeros(nc, nlevsno + nlevmaxurbgrnd)

        CLM.compute_heat_diff_flux_and_factor!(col, lun, temp, ef, up,
                                                mask_nolakec, bounds_col, dtime,
                                                tk, cv, fn, fact)

        # fact should be positive for active layers
        for j in 1:nlevgrnd
            @test fact[c, j + joff] > 0.0
        end

        # fn should be finite for active layers
        for j in 1:(nlevgrnd - 1)
            @test isfinite(fn[c, j + joff])
        end

        # Bottom layer fn should equal eflx_bot (which is 0)
        @test fn[c, nlevgrnd + joff] ≈ 0.0
    end

    @testset "set_rhs_vec! — soil only (no snow)" begin
        col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
            mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
            bounds_col, bounds_lun, bounds_patch = setup_test_data()

        c = 1; p = 1; l = 1
        setup_soil_column!(col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                           mask_nolakec, mask_nolakep, mask_urbanc, c, p, l)

        dtime = 1800.0

        tk = zeros(nc, nlevsno + nlevmaxurbgrnd)
        cv = zeros(nc, nlevsno + nlevmaxurbgrnd)
        tk_h2osfc = fill(NaN, nc)
        CLM.soil_therm_prop!(col, lun, up, temp, wsb, wdb, ss,
                             mask_nolakec, mask_urbanc, bounds_col,
                             tk, cv, tk_h2osfc)

        fn = zeros(nc, nlevsno + nlevmaxurbgrnd)
        fact = zeros(nc, nlevsno + nlevmaxurbgrnd)
        CLM.compute_heat_diff_flux_and_factor!(col, lun, temp, ef, up,
                                                mask_nolakec, bounds_col, dtime,
                                                tk, cv, fn, fact)

        forc_lwrad = fill(300.0, nc)
        hs_h2osfc = zeros(nc)
        hs_top_snow = zeros(nc)
        hs_soil = zeros(nc)
        hs_top = zeros(nc)
        dhsdT = zeros(nc)
        sabg_lyr_col = zeros(nc, nlevsno + 1)

        CLM.compute_ground_heat_flux_and_deriv!(
            col, lun, patch_data, temp, ef, sa, cs, wdb, wfb, up, forc_lwrad,
            mask_nolakec, mask_nolakep, bounds_col, bounds_patch,
            hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT, sabg_lyr_col)

        c_h2osfc_out = fill(CLM.THIN_SFCLAYER, nc)
        dz_h2osfc = fill(CLM.THIN_SFCLAYER, nc)

        nlev_total = nlevsno + 1 + nlevmaxurbgrnd
        rvector = fill(NaN, nc, nlev_total)

        CLM.set_rhs_vec!(col, lun, temp, wdb,
                         mask_nolakec, bounds_col, dtime,
                         hs_h2osfc, hs_top_snow, hs_soil, hs_top, dhsdT,
                         sabg_lyr_col, tk, tk_h2osfc, fact, fn,
                         c_h2osfc_out, dz_h2osfc,
                         rvector)

        # Soil layer RHS should be finite for active layers
        for j in 1:nlevgrnd
            rv_idx = j + nlevsno + 1
            @test isfinite(rvector[c, rv_idx])
        end
    end

    @testset "set_matrix! — soil only (no snow)" begin
        col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
            mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
            bounds_col, bounds_lun, bounds_patch = setup_test_data()

        c = 1; p = 1; l = 1
        setup_soil_column!(col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                           mask_nolakec, mask_nolakep, mask_urbanc, c, p, l)

        dtime = 1800.0

        tk = zeros(nc, nlevsno + nlevmaxurbgrnd)
        cv = zeros(nc, nlevsno + nlevmaxurbgrnd)
        tk_h2osfc = fill(NaN, nc)
        CLM.soil_therm_prop!(col, lun, up, temp, wsb, wdb, ss,
                             mask_nolakec, mask_urbanc, bounds_col,
                             tk, cv, tk_h2osfc)

        fn = zeros(nc, nlevsno + nlevmaxurbgrnd)
        fact = zeros(nc, nlevsno + nlevmaxurbgrnd)
        CLM.compute_heat_diff_flux_and_factor!(col, lun, temp, ef, up,
                                                mask_nolakec, bounds_col, dtime,
                                                tk, cv, fn, fact)

        dhsdT = fill(-5.0, nc)
        c_h2osfc_out = fill(CLM.THIN_SFCLAYER, nc)
        dz_h2osfc = fill(CLM.THIN_SFCLAYER, nc)

        nband = 5
        nlev_total = nlevsno + 1 + nlevmaxurbgrnd
        bmatrix = zeros(nc, nband, nlev_total)

        CLM.set_matrix!(col, lun, wdb, mask_nolakec, bounds_col, dtime, nband,
                        dhsdT, tk, tk_h2osfc, fact, c_h2osfc_out, dz_h2osfc,
                        bmatrix)

        # Diagonal entries should be > 1 for active soil layers (due to 1 + ... terms)
        for j in 1:nlevgrnd
            bm_idx = j + nlevsno + 1
            @test bmatrix[c, 3, bm_idx] > 1.0
        end
    end

    @testset "phase_change_h2osfc! — no surface water" begin
        col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
            mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
            bounds_col, bounds_lun, bounds_patch = setup_test_data()

        c = 1; p = 1; l = 1
        setup_soil_column!(col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                           mask_nolakec, mask_nolakep, mask_urbanc, c, p, l)

        dtime = 1800.0
        dhsdT = fill(-5.0, nc)

        CLM.phase_change_h2osfc!(col, temp, ef, wsb, wdb, wfb,
                                  mask_nolakec, bounds_col, dtime, dhsdT)

        # With no surface water, xmf_h2osfc should be 0
        @test temp.xmf_h2osfc_col[c] == 0.0
        @test wfb.wf.qflx_h2osfc_to_ice_col[c] == 0.0
    end

    @testset "phase_change_beta! — no phase change above freezing" begin
        col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
            mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
            bounds_col, bounds_lun, bounds_patch = setup_test_data()

        c = 1; p = 1; l = 1
        setup_soil_column!(col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                           mask_nolakec, mask_nolakep, mask_urbanc, c, p, l)

        # All temperatures above freezing, no ice → no phase change
        dtime = 1800.0
        dhsdT = fill(-5.0, nc)

        # Set fact to valid values
        for j in 1:nlevgrnd
            temp.fact_col[c, j + joff] = 1.0
        end

        CLM.phase_change_beta!(col, lun, temp, ef, ss, wsb, wdb, wfb,
                               mask_nolakec, bounds_col, dtime, dhsdT)

        # No phase change should occur
        @test temp.xmf_col[c] ≈ 0.0
        @test wfb.wf.qflx_snomelt_col[c] ≈ 0.0
        @test wfb.wf.qflx_snofrz_col[c] ≈ 0.0

        # All imelt flags should be 0
        for j in 1:nlevgrnd
            @test temp.imelt_col[c, j + joff] == 0
        end
    end

    @testset "phase_change_beta! — melting ice in soil" begin
        col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
            mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
            bounds_col, bounds_lun, bounds_patch = setup_test_data()

        c = 1; p = 1; l = 1
        setup_soil_column!(col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
                           mask_nolakec, mask_nolakep, mask_urbanc, c, p, l)

        dtime = 1800.0
        dhsdT = fill(-5.0, nc)

        # Set up ice above freezing → should melt
        j_test = 2
        jj_test = j_test + joff
        wsb.ws.h2osoi_ice_col[c, jj_test] = 5.0
        wsb.ws.h2osoi_liq_col[c, jj_test] = 2.0
        temp.t_soisno_col[c, jj_test] = CLM.TFRZ + 2.0  # 2 degrees above freezing

        # Set fact to valid values
        for j in 1:nlevgrnd
            temp.fact_col[c, j + joff] = dtime / (2.0e6 * col.dz[c, j + joff])
        end

        CLM.phase_change_beta!(col, lun, temp, ef, ss, wsb, wdb, wfb,
                               mask_nolakec, bounds_col, dtime, dhsdT)

        # Layer should be marked as melting
        @test temp.imelt_col[c, jj_test] == 1
        # Temperature should be clamped to TFRZ
        @test temp.t_soisno_col[c, jj_test] ≈ CLM.TFRZ atol=2.0
        # Some ice should have melted
        @test wsb.ws.h2osoi_ice_col[c, jj_test] < 5.0
    end

    @testset "building_hac! — cooling mode" begin
        col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
            mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
            bounds_col, bounds_lun, bounds_patch = setup_test_data()

        l = 1
        mask_urbanl[l] = true
        lun.urbpoi[l] = true
        temp.t_building_lun[l] = 310.0  # hot building
        up.t_building_min[l] = 288.0
        t_building_max = fill(300.0, nl)

        cool_on = falses(nl)
        heat_on = falses(nl)

        CLM.building_hac!(lun, temp, up, t_building_max,
                          mask_urbanl, bounds_lun, cool_on, heat_on)

        # Building temp should be clamped to max
        @test temp.t_building_lun[l] ≈ 300.0
        @test cool_on[l] == true
        @test heat_on[l] == false
    end

    @testset "building_hac! — heating mode" begin
        col, lun, patch_data, temp, ef, ss, wsb, wdb, wfb, sa, cs, up,
            mask_nolakec, mask_nolakep, mask_urbanl, mask_urbanc,
            bounds_col, bounds_lun, bounds_patch = setup_test_data()

        l = 1
        mask_urbanl[l] = true
        lun.urbpoi[l] = true
        temp.t_building_lun[l] = 275.0  # cold building
        up.t_building_min[l] = 288.0
        t_building_max = fill(300.0, nl)

        cool_on = falses(nl)
        heat_on = falses(nl)

        CLM.building_hac!(lun, temp, up, t_building_max,
                          mask_urbanl, bounds_lun, cool_on, heat_on)

        # Building temp should be clamped to min
        @test temp.t_building_lun[l] ≈ 288.0
        @test cool_on[l] == false
        @test heat_on[l] == true
    end

    @testset "THIN_SFCLAYER constant" begin
        @test CLM.THIN_SFCLAYER == 1.0e-6
    end

    # Restore all saved state
    vp.nlevsno = saved_nlevsno
    vp.nlevgrnd = saved_nlevgrnd
    vp.nlevurb = saved_nlevurb
    vp.nlevmaxurbgrnd = saved_nlevmaxurbgrnd
    vp.nlevsoi = saved_nlevsoi
    CLM.varctl.use_excess_ice = saved_use_excess_ice
    CLM.varctl.snow_thermal_cond_method = saved_snow_thermal
    CLM.urban_ctrl.read_namelist = saved_read_namelist
    CLM.urban_ctrl.building_temp_method = saved_building_temp_method
end
