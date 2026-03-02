@testset "Hydrology Drainage Module" begin
    # ------------------------------------------------------------------
    # Tests for HydrologyDrainageMod port.
    # Verifies:
    #   1. compute_wetland_ice_hydrology!
    #   2. compute_total_runoff!
    #   3. compute_water_mass_non_lake! (stub)
    #   4. adjust_runoff_terms! (stub)
    #   5. hydrology_drainage! (orchestrator smoke test)
    # ------------------------------------------------------------------

    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    nlevsno  = CLM.varpar.nlevsno
    nlevsoi  = CLM.varpar.nlevsoi
    nlevgrnd = CLM.varpar.nlevgrnd
    nlevurb  = CLM.varpar.nlevurb
    nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd

    # ------------------------------------------------------------------
    # 1. compute_wetland_ice_hydrology!
    # ------------------------------------------------------------------
    @testset "compute_wetland_ice_hydrology!" begin
        nc = 4  # 4 columns: wetland, ice, urban non-perv, soil
        bounds = 1:nc

        # Allocate flux vectors
        qflx_drain         = fill(1.0, nc)
        qflx_drain_perched = fill(1.0, nc)
        qflx_surf          = fill(1.0, nc)
        qflx_infl          = fill(1.0, nc)
        qflx_qrgwl         = fill(0.0, nc)
        qflx_latflow_out   = fill(1.0, nc)
        qflx_rsub_sat      = fill(0.0, nc)
        qflx_evap_tot      = fill(0.5, nc)
        qflx_snwcp_ice     = fill(0.1, nc)
        qflx_snwcp_disc_ice = fill(0.01, nc)
        qflx_snwcp_disc_liq = fill(0.02, nc)
        qflx_ice_runoff_snwcp = fill(0.0, nc)

        forc_rain  = fill(2.0, nc)
        forc_snow  = fill(1.0, nc)
        # 2 gridcells, columns 1-2 map to gridcell 1, columns 3-4 to gridcell 2
        qflx_floodg = [0.5, 0.3]

        begwb = fill(100.0, nc)
        endwb = fill(110.0, nc)

        # Column-landunit-gridcell mapping
        col_landunit = [1, 2, 3, 4]
        col_gridcell = [1, 1, 2, 2]
        col_itype    = [1, 1, CLM.ICOL_ROOF, 1]  # col 3 is urban roof

        # Landunit types
        lun_itype = [CLM.ISTWET, CLM.ISTICE, CLM.ISTURB_TBD, CLM.ISTSOIL]
        lun_urbpoi = BitVector([false, false, true, false])

        mask_nolake = BitVector([true, true, true, true])
        dtime = 1800.0

        CLM.compute_wetland_ice_hydrology!(
            qflx_drain, qflx_drain_perched, qflx_surf, qflx_infl,
            qflx_qrgwl, qflx_latflow_out, qflx_rsub_sat,
            qflx_evap_tot, qflx_snwcp_ice,
            qflx_snwcp_disc_ice, qflx_snwcp_disc_liq,
            qflx_ice_runoff_snwcp,
            forc_rain, forc_snow, qflx_floodg,
            begwb, endwb,
            col_landunit, col_gridcell, col_itype,
            lun_itype, lun_urbpoi,
            mask_nolake, bounds, dtime)

        # Column 1 (wetland): all drainage/surf/infl zeroed
        @test qflx_drain[1] == 0.0
        @test qflx_drain_perched[1] == 0.0
        @test qflx_surf[1] == 0.0
        @test qflx_infl[1] == 0.0
        @test qflx_latflow_out[1] == 0.0
        # qflx_qrgwl = rain + snow + flood - evap - snwcp_ice - disc_ice - disc_liq - (endwb - begwb)/dtime
        expected_qrgwl1 = 2.0 + 1.0 + 0.5 - 0.5 - 0.1 - 0.01 - 0.02 - (110.0 - 100.0) / 1800.0
        @test qflx_qrgwl[1] ≈ expected_qrgwl1

        # Column 2 (ice): same pattern
        @test qflx_drain[2] == 0.0
        @test qflx_surf[2] == 0.0
        expected_qrgwl2 = 2.0 + 1.0 + 0.5 - 0.5 - 0.1 - 0.01 - 0.02 - (110.0 - 100.0) / 1800.0
        @test qflx_qrgwl[2] ≈ expected_qrgwl2

        # Column 3 (urban roof, not pervious): drain_perched=0, rsub_sat=SPVAL, infl=0
        @test qflx_drain_perched[3] == 0.0
        @test qflx_rsub_sat[3] ≈ CLM.SPVAL
        @test qflx_infl[3] == 0.0
        # drain should be unchanged (1.0) since this is urban, not wetland/ice
        @test qflx_drain[3] == 1.0

        # Column 4 (soil): no changes to drainage fluxes
        @test qflx_drain[4] == 1.0
        @test qflx_drain_perched[4] == 1.0
        @test qflx_surf[4] == 1.0
        @test qflx_infl[4] == 1.0

        # All columns: qflx_ice_runoff_snwcp = qflx_snwcp_ice
        @test qflx_ice_runoff_snwcp[1] ≈ 0.1
        @test qflx_ice_runoff_snwcp[2] ≈ 0.1
        @test qflx_ice_runoff_snwcp[3] ≈ 0.1
        @test qflx_ice_runoff_snwcp[4] ≈ 0.1
    end

    # ------------------------------------------------------------------
    # 2. compute_total_runoff!
    # ------------------------------------------------------------------
    @testset "compute_total_runoff!" begin
        nc = 3  # urban, soil, crop
        bounds = 1:nc

        qflx_runoff   = zeros(nc)
        qflx_runoff_u = zeros(nc)
        qflx_runoff_r = zeros(nc)

        qflx_drain         = [1.0, 2.0, 3.0]
        qflx_surf          = [0.5, 1.0, 0.5]
        qflx_qrgwl         = [0.0, 0.1, 0.0]
        qflx_drain_perched = [0.0, 0.5, 0.2]

        col_landunit = [1, 2, 3]
        lun_itype  = [CLM.ISTURB_TBD, CLM.ISTSOIL, CLM.ISTCROP]
        lun_urbpoi = BitVector([true, false, false])
        mask_nolake = BitVector([true, true, true])

        CLM.compute_total_runoff!(
            qflx_runoff, qflx_runoff_u, qflx_runoff_r,
            qflx_drain, qflx_surf, qflx_qrgwl, qflx_drain_perched,
            col_landunit, lun_itype, lun_urbpoi,
            mask_nolake, bounds)

        # Column 1 (urban): runoff = 1.0 + 0.5 + 0.0 + 0.0 = 1.5
        @test qflx_runoff[1] ≈ 1.5
        @test qflx_runoff_u[1] ≈ 1.5

        # Column 2 (soil): runoff = 2.0 + 1.0 + 0.1 + 0.5 = 3.6
        @test qflx_runoff[2] ≈ 3.6
        @test qflx_runoff_r[2] ≈ 3.6

        # Column 3 (crop): runoff = 3.0 + 0.5 + 0.0 + 0.2 = 3.7
        @test qflx_runoff[3] ≈ 3.7
        @test qflx_runoff_r[3] ≈ 3.7
    end

    # ------------------------------------------------------------------
    # 3. compute_water_mass_non_lake! (stub)
    # ------------------------------------------------------------------
    @testset "compute_water_mass_non_lake! (stub)" begin
        nc = 2
        np = 1; nl = 1; ng = 1

        waterstatebulk = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(waterstatebulk, nc, np, nl, ng)

        waterdiagbulk = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(waterdiagbulk, nc, np, nl, ng)

        endwb = fill(NaN, nc)
        mask_nolake = BitVector([true, true])
        bounds = 1:nc

        CLM.compute_water_mass_non_lake!(endwb, waterstatebulk, waterdiagbulk,
            mask_nolake, bounds)

        # Stub sets endwb to 0.0
        @test endwb[1] == 0.0
        @test endwb[2] == 0.0
    end

    # ------------------------------------------------------------------
    # 4. adjust_runoff_terms! (stub)
    # ------------------------------------------------------------------
    @testset "adjust_runoff_terms! (stub)" begin
        nc = 2
        np = 1; nl = 1; ng = 1

        waterfluxbulk = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(waterfluxbulk, nc, np, nl, ng)

        qflx_qrgwl = [1.0, 2.0]
        qflx_ice_runoff_snwcp = [0.5, 0.6]
        mask_do_smb = BitVector([true, true])
        bounds = 1:nc

        # Should run without error (stub — no-op)
        CLM.adjust_runoff_terms!(qflx_qrgwl, qflx_ice_runoff_snwcp,
            waterfluxbulk, mask_do_smb, bounds)

        # Values unchanged (stub is a no-op)
        @test qflx_qrgwl[1] == 1.0
        @test qflx_ice_runoff_snwcp[2] == 0.6
    end

    # ------------------------------------------------------------------
    # 5. hydrology_drainage! (orchestrator smoke test)
    # ------------------------------------------------------------------
    @testset "hydrology_drainage! smoke test" begin
        # Test the orchestrator with hydrology mask empty (skip drainage call)
        # to verify the volumetric soil water update, wetland hydrology, and
        # runoff computation work correctly together.
        nc = 2
        np = 1; nl = 1; ng = 1

        temperature = CLM.TemperatureData()
        CLM.temperature_init!(temperature, np, nc, nl, ng)

        soilhydrology = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(soilhydrology, nc)

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

        lun_data = CLM.LandunitData()
        CLM.landunit_init!(lun_data, nl)

        # Setup: 2 soil columns, 1 landunit, 1 gridcell
        col_data.landunit .= 1
        col_data.gridcell .= 1
        col_data.itype .= 1  # soil type
        col_data.topo_slope .= 0.01
        lun_data.itype .= CLM.ISTSOIL
        lun_data.urbpoi .= false

        # Initialize temperature
        temperature.t_soisno_col .= 280.0

        # Initialize water state arrays to reasonable values
        waterstatebulk.ws.h2osoi_liq_col .= 10.0
        waterstatebulk.ws.h2osoi_ice_col .= 1.0
        waterstatebulk.ws.h2osoi_vol_col .= 0.3
        waterstatebulk.ws.h2osfc_col .= 0.0
        waterstatebulk.ws.wa_col .= 5000.0

        # Initialize water flux arrays to zero
        wf = waterfluxbulk.wf
        wf.qflx_drain_col .= 0.5         # pre-set drainage from upstream call
        wf.qflx_drain_perched_col .= 0.1
        wf.qflx_surf_col .= 0.2
        wf.qflx_infl_col .= 0.0
        wf.qflx_qrgwl_col .= 0.0
        wf.qflx_latflow_out_col .= 0.0
        wf.qflx_rsub_sat_col .= 0.0
        wf.qflx_evap_tot_col .= 0.0
        wf.qflx_snwcp_ice_col .= 0.0
        wf.qflx_snwcp_discarded_ice_col .= 0.0
        wf.qflx_snwcp_discarded_liq_col .= 0.0
        wf.qflx_ice_runoff_snwcp_col .= 0.0
        wf.qflx_runoff_col .= 0.0
        wf.qflx_runoff_u_col .= 0.0
        wf.qflx_runoff_r_col .= 0.0

        # Water balance
        waterbalancebulk.begwb_col .= 0.0
        waterbalancebulk.endwb_col .= 0.0

        # Column geometry
        nlevtot = nlevsno + nlevmaxurbgrnd
        for c in 1:nc
            for j in 1:nlevtot
                col_data.dz[c, j] = 0.1
                col_data.z[c, j] = 0.05 + (j - 1) * 0.1
            end
        end

        # Use empty hydrology mask to skip drainage call
        mask_nolake    = BitVector([true, true])
        mask_hydrology = BitVector([false, false])  # skip drainage
        mask_urban     = BitVector([false, false])
        mask_do_smb    = BitVector([false, false])
        bounds = 1:nc
        dtime = 1800.0

        forc_rain  = fill(0.001, nc)
        forc_snow  = fill(0.0, nc)
        qflx_floodg = [0.0]

        CLM.hydrology_drainage!(
            temperature, soilhydrology, soilstate,
            waterstatebulk, waterdiagbulk, waterbalancebulk, waterfluxbulk,
            col_data, lun_data,
            mask_nolake, mask_hydrology, mask_urban, mask_do_smb,
            forc_rain, forc_snow, qflx_floodg,
            bounds, dtime, nlevsno, nlevsoi, nlevgrnd, nlevurb;
            use_vichydro=false,
            use_aquifer_layer=true,
            use_hillslope_routing=false)

        # h2osoi_vol should have been updated from h2osoi_liq/ice
        joff = nlevsno
        expected_vol = waterstatebulk.ws.h2osoi_liq_col[1, 1 + joff] /
                       (col_data.dz[1, 1 + joff] * CLM.DENH2O) +
                       waterstatebulk.ws.h2osoi_ice_col[1, 1 + joff] /
                       (col_data.dz[1, 1 + joff] * CLM.DENICE)
        @test waterstatebulk.ws.h2osoi_vol_col[1, 1] ≈ expected_vol

        # Runoff = drain + surf + qrgwl + drain_perched
        # For soil columns (not wetland/ice), fluxes unchanged: 0.5 + 0.2 + 0.0 + 0.1 = 0.8
        @test wf.qflx_runoff_col[1] ≈ 0.8
        @test wf.qflx_runoff_r_col[1] ≈ 0.8  # soil → rural runoff

        # ice_runoff_snwcp = snwcp_ice = 0.0
        @test wf.qflx_ice_runoff_snwcp_col[1] ≈ 0.0
    end
end
