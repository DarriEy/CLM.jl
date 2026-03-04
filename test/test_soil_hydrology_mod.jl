@testset "Soil Hydrology Module" begin
    # ------------------------------------------------------------------
    # Tests for SoilHydrologyMod port.
    # Verifies:
    #   1. SoilHydrologyParams defaults
    #   2. SoilHydrologyConfig creation
    #   3. set_soil_water_fractions!
    #   4. set_floodc!
    #   5. set_qflx_inputs!
    #   6. infiltration!
    #   7. total_surface_runoff!
    #   8. update_urban_ponding!
    #   9. renew_condensation! basic check
    #   10. theta_based_water_table!
    #   11. calc_irrig_withdrawals!
    #   12. withdraw_groundwater_irrigation!
    # ------------------------------------------------------------------

    # Ensure varpar is initialized for dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevgrnd = CLM.varpar.nlevgrnd
    nlevsoi  = CLM.varpar.nlevsoi
    nlayer   = CLM.NLAYER
    nlayert  = CLM.varpar.nlayert

    # ------------------------------------------------------------------
    # 1. SoilHydrologyParams defaults
    # ------------------------------------------------------------------
    @testset "SoilHydrologyParams defaults" begin
        p = CLM.SoilHydrologyParams()
        @test p.aq_sp_yield_min == 0.01
        @test p.n_baseflow == 1.0
        @test p.perched_baseflow_scalar == 5.0e-4
        @test p.e_ice == 6.0
    end

    # ------------------------------------------------------------------
    # 2. SoilHydrologyConfig creation
    # ------------------------------------------------------------------
    @testset "SoilHydrologyConfig defaults" begin
        cfg = CLM.SoilHydrologyConfig()
        @test cfg.head_gradient_method == CLM.HEAD_GRADIENT_DARCY
        @test cfg.transmissivity_method == CLM.TRANSMISSIVITY_LAYERSUM
        @test cfg.baseflow_scalar == 1.0e-2
    end

    @testset "init_soil_hydrology_config" begin
        cfg = CLM.init_soil_hydrology_config(
            head_gradient_method = CLM.HEAD_GRADIENT_KINEMATIC,
            baseflow_scalar = 2.0e-2
        )
        @test cfg.head_gradient_method == CLM.HEAD_GRADIENT_KINEMATIC
        @test cfg.baseflow_scalar == 2.0e-2
        @test CLM.HEAD_GRADIENT_METHOD[] == CLM.HEAD_GRADIENT_KINEMATIC
        @test CLM.BASEFLOW_SCALAR[] == 2.0e-2

        # Reset to defaults
        CLM.init_soil_hydrology_config()
    end

    # ------------------------------------------------------------------
    # 3. set_soil_water_fractions!
    # ------------------------------------------------------------------
    @testset "set_soil_water_fractions!" begin
        nc = 3
        nlevsno_test = 12
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)

        ss = CLM.SoilStateData()
        ss.watsat_col = fill(0.4, nc, nlevgrnd)
        ss.eff_porosity_col = fill(0.35, nc, nlevgrnd)

        ws = CLM.WaterStateData()
        # Snow+soil arrays: nlevsno + nlevsoi columns
        ntot = nlevsno_test + nlevsoi
        ws.h2osoi_liq_col = fill(10.0, nc, ntot)
        ws.h2osoi_ice_col = fill(2.0, nc, ntot)
        ws.excess_ice_col = fill(0.0, nc, nlevsoi)

        col_dz = fill(0.1, nc, ntot)
        mask = trues(nc)

        CLM.set_soil_water_fractions!(sh, ss, ws, col_dz, mask, 1:nc, nlevsoi, nlevsno_test)

        # Check that icefrac was computed (not NaN) for active layers
        for c in 1:nc
            for j in 1:nlevsoi
                @test !isnan(sh.icefrac_col[c, j])
                @test sh.icefrac_col[c, j] >= 0.0
                @test sh.icefrac_col[c, j] <= 1.0
            end
        end

        # Check effective porosity was updated
        for c in 1:nc
            for j in 1:nlevsoi
                @test ss.eff_porosity_col[c, j] >= 0.01
                @test ss.eff_porosity_col[c, j] <= 0.4
            end
        end

        # Specific check: with h2osoi_ice=2.0, dz=0.1, DENICE=917
        # vol_ice = min(0.4, 2.0/(0.1*917)) = min(0.4, 0.0218) = 0.0218
        # eff_porosity = max(0.01, 0.4 - 0.0218) = 0.3782
        # icefrac = min(1.0, 0.0218/0.4) = 0.0545
        @test ss.eff_porosity_col[1, 1] ≈ max(0.01, 0.4 - 2.0 / (0.1 * CLM.DENICE)) atol=1e-10
        @test sh.icefrac_col[1, 1] ≈ min(1.0, (2.0 / (0.1 * CLM.DENICE)) / 0.4) atol=1e-10
    end

    # ------------------------------------------------------------------
    # 4. set_floodc!
    # ------------------------------------------------------------------
    @testset "set_floodc!" begin
        nc = 4
        qflx_floodc = zeros(nc)
        qflx_floodg = [0.5, 0.3]  # 2 gridcells
        col_gridcell = [1, 1, 2, 2]
        col_itype = [1, CLM.ICOL_SUNWALL, 1, CLM.ICOL_SHADEWALL]
        mask = trues(nc)

        CLM.set_floodc!(qflx_floodc, qflx_floodg, col_gridcell, col_itype, mask, 1:nc)

        @test qflx_floodc[1] ≈ 0.5  # normal column, gridcell 1
        @test qflx_floodc[2] ≈ 0.0  # sunwall -> 0
        @test qflx_floodc[3] ≈ 0.3  # normal column, gridcell 2
        @test qflx_floodc[4] ≈ 0.0  # shadewall -> 0
    end

    # ------------------------------------------------------------------
    # 5. infiltration!
    # ------------------------------------------------------------------
    @testset "infiltration!" begin
        nc = 3
        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)

        wfb.qflx_in_soil_limited_col .= [0.1, 0.2, 0.3]
        wfb.qflx_h2osfc_drain_col .= [0.01, 0.02, 0.03]

        mask = trues(nc)
        CLM.infiltration!(wfb, mask, 1:nc)

        @test wfb.wf.qflx_infl_col[1] ≈ 0.11
        @test wfb.wf.qflx_infl_col[2] ≈ 0.22
        @test wfb.wf.qflx_infl_col[3] ≈ 0.33
    end

    # ------------------------------------------------------------------
    # 6. theta_based_water_table!
    # ------------------------------------------------------------------
    @testset "theta_based_water_table!" begin
        nc = 2
        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)

        ss = CLM.SoilStateData()
        ss.watsat_col = fill(0.4, nc, nlevgrnd)

        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, nc, 1, 1)
        # Set soil moisture to be saturated everywhere
        for c in 1:nc, j in 1:nlevsoi
            ws.h2osoi_vol_col[c, j] = 0.4  # = watsat
        end

        # Column geometry
        col_dz = fill(0.1, nc, nlevsoi)
        col_z = zeros(nc, nlevsoi)
        col_zi = zeros(nc, nlevsoi)
        for j in 1:nlevsoi
            col_z[:, j] .= 0.05 + (j - 1) * 0.1
            col_zi[:, j] .= j * 0.1
        end

        col_nbedrock = fill(nlevsoi, nc)
        mask = trues(nc)

        # Make columns saturated: h2osoi_liq high, h2osoi_ice=0
        ws.h2osoi_liq_col = fill(0.4 * 0.1 * CLM.DENH2O, nc, nlevsoi)  # saturated
        ws.h2osoi_ice_col = fill(0.0, nc, nlevsoi)

        CLM.theta_based_water_table!(sh, ss, ws, col_dz, col_z, col_zi,
            col_nbedrock, mask, 1:nc, nlevsoi)

        # When fully saturated, water table should be at top (zi[c,1])
        for c in 1:nc
            @test sh.zwt_col[c] ≈ col_zi[c, 1]
        end

        # Make columns dry: h2osoi_liq very low
        ws.h2osoi_liq_col = fill(0.01, nc, nlevsoi)  # very dry
        ws.h2osoi_ice_col = fill(0.0, nc, nlevsoi)

        CLM.theta_based_water_table!(sh, ss, ws, col_dz, col_z, col_zi,
            col_nbedrock, mask, 1:nc, nlevsoi)

        # When dry, water table should be at bottom
        for c in 1:nc
            @test sh.zwt_col[c] ≈ col_zi[c, nlevsoi]
        end
    end

    # ------------------------------------------------------------------
    # 7. withdraw_groundwater_irrigation!
    # ------------------------------------------------------------------
    @testset "withdraw_groundwater_irrigation!" begin
        nc = 2
        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, nc, nc, 1, 1)
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, nc, 1, 1)

        # Set initial state
        ws.wa_col .= 5000.0
        ws.h2osoi_liq_col .= 10.0

        # Set irrigation fluxes
        wf.qflx_gw_uncon_irrig_lyr_col .= 0.0
        wf.qflx_gw_uncon_irrig_lyr_col[1, 1] = 0.01  # mm/s from layer 1
        wf.qflx_gw_con_irrig_col .= 0.0
        wf.qflx_gw_con_irrig_col[1] = 0.005  # mm/s confined

        mask = trues(nc)
        dtime = 1800.0  # 30 min

        CLM.withdraw_groundwater_irrigation!(wf, ws, mask, 1:nc, nlevsoi, dtime)

        # h2osoi_liq should be reduced in layer 1 for column 1
        @test ws.h2osoi_liq_col[1, 1] ≈ 10.0 - 0.01 * 1800.0
        # h2osoi_liq for column 2 unchanged
        @test ws.h2osoi_liq_col[2, 1] ≈ 10.0

        # wa should be reduced for column 1
        @test ws.wa_col[1] ≈ 5000.0 - 0.005 * 1800.0
        # wa for column 2 unchanged
        @test ws.wa_col[2] ≈ 5000.0
    end

    # ------------------------------------------------------------------
    # 8. renew_condensation! basic
    # ------------------------------------------------------------------
    @testset "renew_condensation! basic" begin
        nc = 2

        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, nc, 1, 1)
        ws.h2osoi_liq_col .= 5.0
        ws.h2osoi_ice_col .= 3.0

        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, nc, 1, 1)
        wdb.frac_h2osfc_col .= 0.0

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
        wfb.wf.qflx_liqdew_to_top_layer_col .= 0.001
        wfb.wf.qflx_soliddew_to_top_layer_col .= 0.0005
        wfb.wf.qflx_solidevap_from_top_layer_col .= 0.0001

        col_snl = fill(0, nc)  # no snow layers
        col_itype = fill(1, nc)
        mask_hydrology = trues(nc)
        mask_urban = falses(nc)
        dtime = 3600.0

        CLM.renew_condensation!(ws, wdb, wfb, col_snl, col_itype,
            mask_hydrology, mask_urban, 1:nc, dtime)

        # Top soil layer at Julia index nlevsno+1 (Fortran layer 1 when snl=0)
        jj = CLM.varpar.nlevsno + 1
        # h2osoi_liq should increase by liqdew * dtime
        @test ws.h2osoi_liq_col[1, jj] ≈ 5.0 + 0.001 * 3600.0
        # h2osoi_ice should increase by soliddew - solidevap
        expected_ice = 3.0 + 0.0005 * 3600.0 - 0.0001 * 3600.0
        @test ws.h2osoi_ice_col[1, jj] ≈ expected_ice
    end

    # ------------------------------------------------------------------
    # 9. Module-level constants
    # ------------------------------------------------------------------
    @testset "Constants" begin
        @test CLM.HEAD_GRADIENT_KINEMATIC == 0
        @test CLM.HEAD_GRADIENT_DARCY == 1
        @test CLM.TRANSMISSIVITY_UNIFORM == 0
        @test CLM.TRANSMISSIVITY_LAYERSUM == 1
        @test CLM.TOLERANCE_SOILHYDRO == 1.0e-12
    end

    # ------------------------------------------------------------------
    # 10. perched_lateral_flow! (non-hillslope path)
    # ------------------------------------------------------------------
    @testset "perched_lateral_flow! non-hillslope" begin
        nc = 2
        nl = 1
        ng = 1

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.frost_table_col .= 0.5   # frost table at 0.5 m
        sh.zwt_col .= 1.0           # water table below frost table
        sh.zwt_perched_col .= 0.3   # perched wt above frost table

        ss = CLM.SoilStateData()
        ss.watsat_col = fill(0.4, nc, nlevgrnd)
        ss.bsw_col = fill(5.0, nc, nlevgrnd)
        ss.hksat_col = fill(0.01, nc, nlevgrnd)
        ss.sucsat_col = fill(100.0, nc, nlevgrnd)
        ss.eff_porosity_col = fill(0.35, nc, nlevgrnd)

        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, nc, nl, ng)
        ws.h2osoi_liq_col = fill(10.0, nc, nlevsoi)
        ws.h2osoi_ice_col = fill(0.0, nc, nlevsoi)
        ws.stream_water_volume_lun = fill(0.0, nl)

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, nl, ng)

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.dz = fill(0.1, nc, size(col.dz, 2))
        col.z = zeros(nc, size(col.z, 2))
        col.zi = zeros(nc, size(col.zi, 2))
        for j in 1:size(col.z, 2)
            col.z[:, j] .= 0.05 + (j - 1) * 0.1
        end
        for j in 1:size(col.zi, 2)
            col.zi[:, j] .= j * 0.1
        end
        col.nbedrock .= nlevsoi
        col.landunit .= 1
        col.gridcell .= 1
        col.is_hillslope_column .= false
        col.active .= true
        col.cold .= CLM.ISPVAL
        col.hill_slope .= 0.1
        col.hill_elev .= 10.0
        col.hill_distance .= 100.0
        col.hill_width .= 50.0
        col.hill_area .= 5000.0
        col.topo_slope .= 5.0

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        lun.stream_channel_length .= 100.0
        lun.stream_channel_width .= 5.0
        lun.stream_channel_depth .= 1.0

        tdepth = fill(0.0, ng)
        tdepthmax = fill(1.0, ng)

        mask = trues(nc)

        CLM.perched_lateral_flow!(sh, ss, ws, wfb,
            col, lun, tdepth, tdepthmax,
            mask, 1:nc, nlevsoi, 1800.0)

        # Since frost_table > zwt_perched and non-hillslope, drainage should be computed
        for c in 1:nc
            @test !isnan(wfb.wf.qflx_drain_perched_col[c])
            @test isfinite(wfb.wf.qflx_drain_perched_col[c])
        end
    end

    # ------------------------------------------------------------------
    # 11. subsurface_lateral_flow! (non-hillslope path)
    # ------------------------------------------------------------------
    @testset "subsurface_lateral_flow! non-hillslope" begin
        nc = 2
        nl = 1
        ng = 1

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.zwt_col .= 0.5           # water table in middle
        sh.frost_table_col .= 1.0   # frost table below wt

        ss = CLM.SoilStateData()
        ss.watsat_col = fill(0.4, nc, nlevgrnd)
        ss.bsw_col = fill(5.0, nc, nlevgrnd)
        ss.hksat_col = fill(0.01, nc, nlevgrnd)
        ss.sucsat_col = fill(100.0, nc, nlevgrnd)
        ss.eff_porosity_col = fill(0.35, nc, nlevgrnd)
        ss.hk_l_col = fill(0.005, nc, nlevgrnd)

        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, nc, nl, ng)
        ws.h2osoi_liq_col = fill(10.0, nc, nlevsoi)
        ws.h2osoi_ice_col = fill(0.0, nc, nlevsoi)
        ws.h2osfc_col .= 0.0
        ws.stream_water_volume_lun = fill(0.0, nl)

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, nl, ng)
        wfb.wf.qflx_snwcp_liq_col .= 0.0

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.dz = fill(0.1, nc, size(col.dz, 2))
        col.z = zeros(nc, size(col.z, 2))
        col.zi = zeros(nc, size(col.zi, 2))
        for j in 1:size(col.z, 2)
            col.z[:, j] .= 0.05 + (j - 1) * 0.1
        end
        for j in 1:size(col.zi, 2)
            col.zi[:, j] .= j * 0.1
        end
        col.nbedrock .= nlevsoi
        col.itype .= 1
        col.landunit .= 1
        col.gridcell .= 1
        col.is_hillslope_column .= false
        col.active .= true
        col.cold .= CLM.ISPVAL
        col.colu .= CLM.ISPVAL
        col.hill_slope .= 0.1
        col.hill_elev .= 10.0
        col.hill_distance .= 100.0
        col.hill_width .= 50.0
        col.hill_area .= 5000.0
        col.topo_slope .= 5.0
        col.wtgcell .= 0.5

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nl)
        lun.stream_channel_length .= 100.0
        lun.stream_channel_width .= 5.0
        lun.stream_channel_depth .= 1.0
        lun.stream_channel_number .= 1.0
        lun.urbpoi .= false

        grc_area = fill(100.0, ng)
        tdepth = fill(0.0, ng)
        tdepthmax = fill(1.0, ng)

        mask_hydro = trues(nc)
        mask_urban = falses(nc)

        CLM.subsurface_lateral_flow!(sh, ss, ws, wfb,
            col, lun, tdepth, tdepthmax, grc_area,
            mask_hydro, mask_urban, 1:nc, nlevsoi, 1800.0)

        # Non-hillslope with zwt <= zi[nbedrock] should produce baseflow
        for c in 1:nc
            @test isfinite(wfb.wf.qflx_drain_col[c])
            @test isfinite(wfb.wf.qflx_latflow_out_col[c])
            @test wfb.wf.qflx_latflow_out_col[c] >= 0.0
        end
    end

end
