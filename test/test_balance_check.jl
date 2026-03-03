@testset "Balance Check (BalanceCheckMod)" begin

    # -----------------------------------------------------------------------
    # Test BalanceCheckInit / Clean / GetSkipSteps
    # -----------------------------------------------------------------------
    @testset "BalanceCheckInit / Clean / GetSkipSteps" begin
        bc = CLM.BalanceCheckData()
        @test bc.skip_steps == -999

        # Before init, GetSkipSteps should error
        @test_throws ErrorException CLM.get_balance_check_skip_steps(bc)

        # Init with 1800s timestep: skip_size=3600 → nint(3600/1800)=2, max(2,2)+1=3
        CLM.balance_check_init!(bc, 1800.0)
        @test bc.skip_steps == 3

        @test CLM.get_balance_check_skip_steps(bc) == 3

        # Init with 600s timestep: nint(3600/600)=6, max(2,6)+1=7
        CLM.balance_check_init!(bc, 600.0)
        @test bc.skip_steps == 7

        # Init with very large timestep: nint(3600/7200)=1 → max(2,1)+1=3
        CLM.balance_check_init!(bc, 7200.0)
        @test bc.skip_steps == 3

        # Clean resets
        CLM.balance_check_clean!(bc)
        @test bc.skip_steps == -999
        @test_throws ErrorException CLM.get_balance_check_skip_steps(bc)
    end

    # -----------------------------------------------------------------------
    # Helper: create minimal test data
    # -----------------------------------------------------------------------
    function make_bc_test_data(; nc=4, np=4, nl=2, ng=1)
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        nlevsno = CLM.varpar.nlevsno
        nlevgrnd = CLM.varpar.nlevgrnd
        nlevtot = nlevsno + nlevgrnd

        # BalanceCheck state
        bc = CLM.BalanceCheckData()
        CLM.balance_check_init!(bc, 1800.0)

        # Column data
        col_data = CLM.ColumnData()
        col_data.gridcell = fill(1, nc)
        col_data.landunit = fill(1, nc)
        col_data.itype    = fill(1, nc)  # default soil
        col_data.active   = fill(true, nc)
        col_data.snl      = fill(0, nc)
        col_data.wtgcell  = fill(1.0 / nc, nc)
        col_data.zi       = zeros(nc, nlevtot + 1)

        # Landunit data
        lun_data = CLM.LandunitData()
        lun_data.itype   = fill(CLM.ISTSOIL, nl)
        lun_data.urbpoi  = falses(nl)
        lun_data.gridcell = fill(1, nl)

        # Patch data
        pat_data = CLM.PatchData()
        pat_data.column   = collect(1:np)
        pat_data.landunit = fill(1, np)
        pat_data.gridcell = fill(1, np)
        pat_data.active   = fill(true, np)

        # Gridcell data
        grc_data = CLM.GridcellData()
        grc_data.area = [1.0e6]  # 1 km^2

        # Masks
        mask_nolake = trues(nc)
        mask_lake = falses(nc)
        mask_allc = trues(nc)

        # Soil hydrology
        soilhyd = CLM.SoilHydrologyData()
        soilhyd.zwt_col = fill(5.0, nc)  # deep water table

        # Lake state
        lakestate = CLM.LakeStateData()

        # Water data (bulk only)
        water = CLM.WaterData()
        CLM.water_init!(water, nc, np, nl, ng)

        # Water balance
        wb = water.waterbalancebulk_inst
        wb.begwb_col  .= 100.0
        wb.endwb_col  .= 100.0
        wb.begwb_grc  .= 100.0
        wb.endwb_grc  .= 100.0
        wb.h2osno_old_col .= 0.0
        wb.errh2o_col .= 0.0
        wb.errh2osno_col .= 0.0
        wb.snow_sources_col .= 0.0
        wb.snow_sinks_col .= 0.0
        wb.wa_reset_nonconservation_gain_col .= 0.0

        # Water flux (bulk)
        wf = water.waterfluxbulk_inst.wf
        for f in fieldnames(typeof(wf))
            v = getfield(wf, f)
            if v isa Vector{Float64} && length(v) == nc
                v .= 0.0
            elseif v isa Vector{Float64} && length(v) == np
                v .= 0.0
            end
        end

        # Water diagnostic bulk
        wdb = water.waterdiagnosticbulk_inst
        wdb.frac_sno_eff_col .= 0.0
        wdb.frac_sno_col .= 0.0
        wdb.snow_depth_col .= 0.0
        wdb.qflx_prec_grnd_col .= 0.0

        # Water state (bulk)
        wstate = water.waterstatebulk_inst.ws
        wstate.aquifer_water_baseline = 0.0
        wstate.wa_col .= 0.0
        wstate.h2osno_no_layers_col .= 0.0
        wstate.stream_water_volume_lun .= 0.0

        # Energy flux
        eflux = CLM.EnergyFluxData()
        CLM.energyflux_init!(eflux, np, nc, nl, ng)
        eflux.dhsdt_canopy_patch .= 0.0
        eflux.eflx_lwrad_out_patch .= 300.0
        eflux.eflx_lwrad_net_patch .= 50.0
        eflux.eflx_sh_tot_patch .= 0.0
        eflux.eflx_lh_tot_patch .= 0.0
        eflux.eflx_soil_grnd_patch .= 0.0
        eflux.eflx_wasteheat_patch .= 0.0
        eflux.eflx_ventilation_patch .= 0.0
        eflux.eflx_heat_from_ac_patch .= 0.0
        eflux.eflx_traffic_patch .= 0.0
        eflux.errsoi_col .= 0.0
        eflux.errsol_patch .= 0.0
        eflux.errseb_patch .= 0.0
        eflux.errlon_patch .= 0.0
        eflux.netrad_patch .= 0.0

        # Solar absorbed
        solarabs = CLM.SolarAbsorbedData()
        CLM.solarabs_init!(solarabs, np, nl)
        solarabs.fsa_patch .= 100.0
        solarabs.fsr_patch .= 50.0
        solarabs.sabv_patch .= 40.0
        solarabs.sabg_patch .= 60.0
        solarabs.sabg_soil_patch .= 60.0
        solarabs.sabg_snow_patch .= 0.0
        solarabs.sabg_chk_patch .= 60.0

        # Canopy state
        canst = CLM.CanopyStateData()
        CLM.canopystate_init!(canst, np)
        canst.elai_patch .= 2.0
        canst.esai_patch .= 0.5

        # Surface albedo
        surfalb = CLM.SurfaceAlbedoData()
        CLM.surfalb_init!(surfalb, np, nc, ng)
        surfalb.fabd_patch .= 0.3
        surfalb.fabi_patch .= 0.1
        surfalb.albd_patch .= 0.2
        surfalb.albi_patch .= 0.15
        surfalb.ftdd_patch .= 0.5
        surfalb.ftid_patch .= 0.1
        surfalb.ftii_patch .= 0.4

        # Forcing arrays (from not-yet-ported types)
        forc_rain_col = fill(0.0, nc)
        forc_snow_col = fill(0.0, nc)
        forc_rain_grc = fill(0.0, ng)
        forc_snow_grc = fill(0.0, ng)
        forc_solad_col = fill(75.0, nc, 2)  # vis + nir direct
        forc_solai_grc = fill(0.0, ng, 2)   # vis + nir diffuse
        forc_lwrad_col = fill(250.0, nc)
        forc_flood_grc = fill(0.0, ng)
        qflx_ice_runoff_col = fill(0.0, nc)
        qflx_evap_tot_grc = fill(0.0, ng)
        qflx_surf_grc = fill(0.0, ng)
        qflx_qrgwl_grc = fill(0.0, ng)
        qflx_drain_grc = fill(0.0, ng)
        qflx_drain_perched_grc = fill(0.0, ng)
        qflx_ice_runoff_grc = fill(0.0, ng)
        qflx_sfc_irrig_grc = fill(0.0, ng)
        qflx_streamflow_grc = fill(0.0, ng)

        bounds_c = 1:nc
        bounds_p = 1:np
        bounds_l = 1:nl
        bounds_g = 1:ng
        dtime = 1800.0

        return (bc=bc, water=water, wb=wb, wf=wf, wdb=wdb, wstate=wstate,
                eflux=eflux, solarabs=solarabs, canst=canst, surfalb=surfalb,
                col_data=col_data, lun_data=lun_data, pat_data=pat_data,
                grc_data=grc_data, soilhyd=soilhyd, lakestate=lakestate,
                mask_nolake=mask_nolake, mask_lake=mask_lake, mask_allc=mask_allc,
                forc_rain_col=forc_rain_col, forc_snow_col=forc_snow_col,
                forc_rain_grc=forc_rain_grc, forc_snow_grc=forc_snow_grc,
                forc_solad_col=forc_solad_col, forc_solai_grc=forc_solai_grc,
                forc_lwrad_col=forc_lwrad_col, forc_flood_grc=forc_flood_grc,
                qflx_ice_runoff_col=qflx_ice_runoff_col,
                qflx_evap_tot_grc=qflx_evap_tot_grc,
                qflx_surf_grc=qflx_surf_grc, qflx_qrgwl_grc=qflx_qrgwl_grc,
                qflx_drain_grc=qflx_drain_grc,
                qflx_drain_perched_grc=qflx_drain_perched_grc,
                qflx_ice_runoff_grc=qflx_ice_runoff_grc,
                qflx_sfc_irrig_grc=qflx_sfc_irrig_grc,
                qflx_streamflow_grc=qflx_streamflow_grc,
                bounds_c=bounds_c, bounds_p=bounds_p,
                bounds_l=bounds_l, bounds_g=bounds_g, dtime=dtime)
    end

    # -----------------------------------------------------------------------
    # Test BeginWaterColumnBalance
    # -----------------------------------------------------------------------
    @testset "BeginWaterColumnBalance" begin
        d = make_bc_test_data()

        # Should run without error
        CLM.begin_water_column_balance!(d.water, d.soilhyd, d.lakestate,
            d.col_data, d.lun_data, d.mask_nolake, d.mask_lake, d.bounds_c;
            use_aquifer_layer=false)

        # begwb should be set (to 0.0 since stubs are used)
        wb = d.water.waterbalancebulk_inst
        for c in d.bounds_c
            @test wb.begwb_col[c] == 0.0
        end
    end

    # -----------------------------------------------------------------------
    # Test WaterGridcellBalance
    # -----------------------------------------------------------------------
    @testset "WaterGridcellBalance begwb" begin
        d = make_bc_test_data()

        dribble_liq = fill(0.0, last(d.bounds_g))
        dribble_ice = fill(0.0, last(d.bounds_g))

        CLM.water_gridcell_balance!(d.water, d.lakestate,
            d.col_data, d.lun_data, d.grc_data,
            d.mask_nolake, d.mask_lake,
            d.bounds_c, d.bounds_l, d.bounds_g,
            "begwb";
            use_aquifer_layer=false,
            use_hillslope_routing=false,
            qflx_liq_dynbal_left_to_dribble=dribble_liq,
            qflx_ice_dynbal_left_to_dribble=dribble_ice)

        wb = d.water.waterbalancebulk_inst
        # With stubs, begwb_grc should be 0.0
        @test wb.begwb_grc[1] == 0.0
    end

    @testset "WaterGridcellBalance endwb" begin
        d = make_bc_test_data()

        dribble_liq = fill(0.0, last(d.bounds_g))
        dribble_ice = fill(0.0, last(d.bounds_g))

        CLM.water_gridcell_balance!(d.water, d.lakestate,
            d.col_data, d.lun_data, d.grc_data,
            d.mask_nolake, d.mask_lake,
            d.bounds_c, d.bounds_l, d.bounds_g,
            "endwb";
            use_aquifer_layer=false,
            use_hillslope_routing=false,
            qflx_liq_dynbal_left_to_dribble=dribble_liq,
            qflx_ice_dynbal_left_to_dribble=dribble_ice)

        wb = d.water.waterbalancebulk_inst
        @test wb.endwb_grc[1] == 0.0
    end

    @testset "WaterGridcellBalance invalid flag" begin
        d = make_bc_test_data()

        @test_throws ErrorException CLM.water_gridcell_balance!(d.water, d.lakestate,
            d.col_data, d.lun_data, d.grc_data,
            d.mask_nolake, d.mask_lake,
            d.bounds_c, d.bounds_l, d.bounds_g,
            "invalid")
    end

    # -----------------------------------------------------------------------
    # Test BalanceCheck — zero fluxes (perfect balance)
    # -----------------------------------------------------------------------
    @testset "BalanceCheck perfect balance" begin
        d = make_bc_test_data()

        # Set begwb == endwb, all fluxes zero → perfect balance
        d.wb.begwb_col .= 100.0
        d.wb.endwb_col .= 100.0
        d.wb.begwb_grc .= 100.0
        d.wb.endwb_grc .= 100.0

        # Solar balance: fsa + fsr = forc_solad(1) + forc_solad(2) + forc_solai(1) + forc_solai(2)
        # 100 + 50 = 75 + 75 + 0 + 0 → errsol = 0
        d.solarabs.fsa_patch .= 100.0
        d.solarabs.fsr_patch .= 50.0
        d.forc_solad_col .= 75.0
        d.forc_solai_grc .= 0.0

        # Longwave: eflx_lwrad_out - eflx_lwrad_net - forc_lwrad = 0
        # 300 - 50 - 250 = 0
        d.eflux.eflx_lwrad_out_patch .= 300.0
        d.eflux.eflx_lwrad_net_patch .= 50.0
        d.forc_lwrad_col .= 250.0

        # Surface energy: sabv + sabg_chk + forc_lwrad - eflx_lwrad_out - sh - lh - soil_grnd - dhsdt = 0
        # 40 + 60 + 250 - 300 - 25 - 25 - 0 - 0 = 0
        d.solarabs.sabv_patch .= 40.0
        d.solarabs.sabg_chk_patch .= 60.0
        d.eflux.eflx_sh_tot_patch .= 25.0
        d.eflux.eflx_lh_tot_patch .= 25.0
        d.eflux.eflx_soil_grnd_patch .= 0.0
        d.eflux.dhsdt_canopy_patch .= 0.0

        # Should not error (balance is perfect)
        CLM.balance_check!(d.bc, d.water.waterfluxbulk_inst,
            d.water.waterstatebulk_inst, d.wb, d.wdb,
            d.eflux, d.solarabs, d.canst, d.surfalb,
            d.col_data, d.lun_data, d.pat_data, d.grc_data,
            d.mask_allc, d.bounds_c, d.bounds_p, d.bounds_g,
            1, 100, d.dtime;
            forc_rain_col=d.forc_rain_col,
            forc_snow_col=d.forc_snow_col,
            forc_rain_grc=d.forc_rain_grc,
            forc_snow_grc=d.forc_snow_grc,
            forc_solad_col=d.forc_solad_col,
            forc_solai_grc=d.forc_solai_grc,
            forc_lwrad_col=d.forc_lwrad_col,
            forc_flood_grc=d.forc_flood_grc,
            qflx_ice_runoff_col=d.qflx_ice_runoff_col,
            qflx_evap_tot_grc=d.qflx_evap_tot_grc,
            qflx_surf_grc=d.qflx_surf_grc,
            qflx_qrgwl_grc=d.qflx_qrgwl_grc,
            qflx_drain_grc=d.qflx_drain_grc,
            qflx_drain_perched_grc=d.qflx_drain_perched_grc,
            qflx_ice_runoff_grc=d.qflx_ice_runoff_grc,
            qflx_sfc_irrig_grc=d.qflx_sfc_irrig_grc,
            qflx_streamflow_grc=d.qflx_streamflow_grc)

        # Water errors should be zero
        for c in d.bounds_c
            @test d.wb.errh2o_col[c] == 0.0
        end

        # Energy errors should be zero
        for p in d.bounds_p
            @test d.eflux.errsol_patch[p] == 0.0
            @test d.eflux.errlon_patch[p] == 0.0
            @test abs(d.eflux.errseb_patch[p]) < 1.0e-10
        end
    end

    # -----------------------------------------------------------------------
    # Test BalanceCheck — water balance error detection
    # -----------------------------------------------------------------------
    @testset "BalanceCheck water error detection" begin
        d = make_bc_test_data()

        # Create a large water imbalance
        d.wb.begwb_col .= 100.0
        d.wb.endwb_col .= 200.0  # 100 mm gain with no fluxes
        d.wb.begwb_grc .= 100.0
        d.wb.endwb_grc .= 200.0

        # Make energy checks pass
        d.solarabs.fsa_patch .= 100.0
        d.solarabs.fsr_patch .= 50.0
        d.forc_solad_col .= 75.0
        d.forc_solai_grc .= 0.0
        d.eflux.eflx_lwrad_out_patch .= 300.0
        d.eflux.eflx_lwrad_net_patch .= 50.0
        d.forc_lwrad_col .= 250.0
        d.solarabs.sabv_patch .= 40.0
        d.solarabs.sabg_chk_patch .= 60.0
        d.eflux.eflx_sh_tot_patch .= 25.0
        d.eflux.eflx_lh_tot_patch .= 25.0

        # DAnstep > skip_steps triggers the error
        @test_throws ErrorException CLM.balance_check!(d.bc, d.water.waterfluxbulk_inst,
            d.water.waterstatebulk_inst, d.wb, d.wdb,
            d.eflux, d.solarabs, d.canst, d.surfalb,
            d.col_data, d.lun_data, d.pat_data, d.grc_data,
            d.mask_allc, d.bounds_c, d.bounds_p, d.bounds_g,
            1, 100, d.dtime;
            forc_rain_col=d.forc_rain_col,
            forc_snow_col=d.forc_snow_col,
            forc_rain_grc=d.forc_rain_grc,
            forc_snow_grc=d.forc_snow_grc,
            forc_solad_col=d.forc_solad_col,
            forc_solai_grc=d.forc_solai_grc,
            forc_lwrad_col=d.forc_lwrad_col,
            forc_flood_grc=d.forc_flood_grc,
            qflx_ice_runoff_col=d.qflx_ice_runoff_col,
            qflx_evap_tot_grc=d.qflx_evap_tot_grc,
            qflx_surf_grc=d.qflx_surf_grc,
            qflx_qrgwl_grc=d.qflx_qrgwl_grc,
            qflx_drain_grc=d.qflx_drain_grc,
            qflx_drain_perched_grc=d.qflx_drain_perched_grc,
            qflx_ice_runoff_grc=d.qflx_ice_runoff_grc,
            qflx_sfc_irrig_grc=d.qflx_sfc_irrig_grc,
            qflx_streamflow_grc=d.qflx_streamflow_grc)
    end

    # -----------------------------------------------------------------------
    # Test that urban wall columns get zero precipitation
    # -----------------------------------------------------------------------
    @testset "Urban wall zero precipitation" begin
        d = make_bc_test_data()

        # Set column 2 as sunwall, column 3 as shadewall
        d.col_data.itype[2] = CLM.ICOL_SUNWALL
        d.col_data.itype[3] = CLM.ICOL_SHADEWALL
        d.lun_data.itype[1] = CLM.ISTURB_HD
        d.lun_data.urbpoi[1] = true

        # Set some precipitation
        d.forc_rain_col .= 5.0
        d.forc_snow_col .= 2.0

        # Run balance check with DAnstep < skip_steps so no error is thrown
        CLM.balance_check!(d.bc, d.water.waterfluxbulk_inst,
            d.water.waterstatebulk_inst, d.wb, d.wdb,
            d.eflux, d.solarabs, d.canst, d.surfalb,
            d.col_data, d.lun_data, d.pat_data, d.grc_data,
            d.mask_allc, d.bounds_c, d.bounds_p, d.bounds_g,
            1, 1, d.dtime;  # DAnstep=1 < skip_steps=3
            forc_rain_col=d.forc_rain_col,
            forc_snow_col=d.forc_snow_col,
            forc_rain_grc=d.forc_rain_grc,
            forc_snow_grc=d.forc_snow_grc,
            forc_solad_col=d.forc_solad_col,
            forc_solai_grc=d.forc_solai_grc,
            forc_lwrad_col=d.forc_lwrad_col,
            forc_flood_grc=d.forc_flood_grc,
            qflx_ice_runoff_col=d.qflx_ice_runoff_col,
            qflx_evap_tot_grc=d.qflx_evap_tot_grc,
            qflx_surf_grc=d.qflx_surf_grc,
            qflx_qrgwl_grc=d.qflx_qrgwl_grc,
            qflx_drain_grc=d.qflx_drain_grc,
            qflx_drain_perched_grc=d.qflx_drain_perched_grc,
            qflx_ice_runoff_grc=d.qflx_ice_runoff_grc,
            qflx_sfc_irrig_grc=d.qflx_sfc_irrig_grc,
            qflx_streamflow_grc=d.qflx_streamflow_grc)

        # Wall columns should have zero water error (no precip, zero fluxes, begwb==endwb)
        # Column 2 (sunwall) and 3 (shadewall): errh2o should be 0
        @test d.wb.errh2o_col[2] ≈ 0.0 atol=1e-15
        @test d.wb.errh2o_col[3] ≈ 0.0 atol=1e-15

        # Non-wall columns should have non-zero errh2o due to precipitation with no storage change
        @test abs(d.wb.errh2o_col[1]) > 0.0
    end

    # -----------------------------------------------------------------------
    # Test BalanceCheck skip during startup
    # -----------------------------------------------------------------------
    @testset "BalanceCheck skip during startup" begin
        d = make_bc_test_data()

        # Create large water imbalance
        d.wb.begwb_col .= 100.0
        d.wb.endwb_col .= 200.0
        d.wb.begwb_grc .= 100.0
        d.wb.endwb_grc .= 200.0

        # DAnstep <= skip_steps → should warn but not error
        # skip_steps = 3, DAnstep = 2
        CLM.balance_check!(d.bc, d.water.waterfluxbulk_inst,
            d.water.waterstatebulk_inst, d.wb, d.wdb,
            d.eflux, d.solarabs, d.canst, d.surfalb,
            d.col_data, d.lun_data, d.pat_data, d.grc_data,
            d.mask_allc, d.bounds_c, d.bounds_p, d.bounds_g,
            1, 2, d.dtime;
            forc_rain_col=d.forc_rain_col,
            forc_snow_col=d.forc_snow_col,
            forc_rain_grc=d.forc_rain_grc,
            forc_snow_grc=d.forc_snow_grc,
            forc_solad_col=d.forc_solad_col,
            forc_solai_grc=d.forc_solai_grc,
            forc_lwrad_col=d.forc_lwrad_col,
            forc_flood_grc=d.forc_flood_grc,
            qflx_ice_runoff_col=d.qflx_ice_runoff_col,
            qflx_evap_tot_grc=d.qflx_evap_tot_grc,
            qflx_surf_grc=d.qflx_surf_grc,
            qflx_qrgwl_grc=d.qflx_qrgwl_grc,
            qflx_drain_grc=d.qflx_drain_grc,
            qflx_drain_perched_grc=d.qflx_drain_perched_grc,
            qflx_ice_runoff_grc=d.qflx_ice_runoff_grc,
            qflx_sfc_irrig_grc=d.qflx_sfc_irrig_grc,
            qflx_streamflow_grc=d.qflx_streamflow_grc)

        # Should not throw — just verify it ran
        @test true
    end

end
