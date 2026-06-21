# =========================================================================
# Soil lateral-flow + irrigation conservation tests
# =========================================================================
#
# Focused water-conservation / finite checks for the lateral-redistribution,
# infiltration-excess routing, and groundwater-irrigation withdrawal routines
# ported from `SoilHydrologyMod.F90`:
#
#   set_floodc!                        — gridcell flood flux -> non-lake columns
#   route_infiltration_excess!         — infiltration-excess runoff routing
#   withdraw_groundwater_irrigation!   — remove gw irrigation from aquifers
#   perched_lateral_flow!              — drainage from perched saturated zone
#   subsurface_lateral_flow!           — subsurface lateral flow / baseflow
#
# The existing test_soil_hydrology_mod.jl exercises these for finiteness; this
# file adds *mass-conservation* invariants:
#   * groundwater withdrawal removes exactly flux*dtime from storage
#   * infiltration-excess routing partitions the input flux without loss
#   * perched lateral flow: water removed from soil storage equals the
#     reported net drainage flux * dtime (per column water balance)
#   * non-hillslope lateral redistribution does not create water globally
# =========================================================================

@testset "Soil Lateral Flow & Irrigation (conservation)" begin
    # Ensure varpar is initialized for dimension parameters (nlevsoi etc.).
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    nlevsoi  = CLM.varpar.nlevsoi
    nlevgrnd = CLM.varpar.nlevgrnd
    nlevsno  = CLM.varpar.nlevsno
    nlevtot  = nlevsoi + nlevsno
    dtime    = 1800.0

    # Helper: build a simple uniform-soil ColumnData for `nc` non-hillslope
    # columns on a single gridcell/landunit, with monotone z/zi grids.
    function build_columns(nc)
        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        col.dz = fill(0.1, nc, size(col.dz, 2))
        col.z  = zeros(nc, size(col.z, 2))
        col.zi = zeros(nc, size(col.zi, 2))
        for j in 1:size(col.z, 2)
            col.z[:, j]  .= 0.05 + (j - 1) * 0.1
        end
        for j in 1:size(col.zi, 2)
            col.zi[:, j] .= j * 0.1
        end
        col.nbedrock .= nlevsoi
        col.itype    .= 1
        col.landunit .= 1
        col.gridcell .= 1
        col.is_hillslope_column .= false
        col.active   .= true
        col.cold     .= CLM.ISPVAL
        col.colu     .= CLM.ISPVAL
        col.hill_slope    .= 0.1
        col.hill_elev     .= 10.0
        col.hill_distance .= 100.0
        col.hill_width    .= 50.0
        col.hill_area     .= 5000.0
        col.topo_slope    .= 5.0
        col.wtgcell       .= 1.0 / nc
        return col
    end

    function build_soilstate(nc)
        ss = CLM.SoilStateData()
        ss.watsat_col       = fill(0.4, nc, nlevgrnd)
        ss.bsw_col          = fill(5.0, nc, nlevgrnd)
        ss.hksat_col        = fill(0.01, nc, nlevgrnd)
        ss.sucsat_col       = fill(100.0, nc, nlevgrnd)
        ss.eff_porosity_col = fill(0.35, nc, nlevgrnd)
        ss.hk_l_col         = fill(0.005, nc, nlevgrnd)
        return ss
    end

    # --------------------------------------------------------------------
    # 1. withdraw_groundwater_irrigation! — exact storage removal
    # --------------------------------------------------------------------
    @testset "groundwater irrigation withdrawal conserves" begin
        nc = 3
        wf = CLM.WaterFluxData()
        CLM.waterflux_init!(wf, nc, nc, 1, 1)
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, nc, 1, 1)

        ws.wa_col         .= 5000.0
        ws.h2osoi_liq_col .= 20.0

        wf.qflx_gw_uncon_irrig_lyr_col .= 0.0
        wf.qflx_gw_uncon_irrig_lyr_col[1, 1] = 0.01
        wf.qflx_gw_uncon_irrig_lyr_col[2, 2] = 0.004
        wf.qflx_gw_con_irrig_col .= 0.0
        wf.qflx_gw_con_irrig_col[1] = 0.005

        # Total water leaving storage should equal total flux * dtime.
        liq_before = sum(ws.h2osoi_liq_col)
        wa_before  = sum(ws.wa_col)
        flux_total = (sum(wf.qflx_gw_uncon_irrig_lyr_col) +
                      sum(wf.qflx_gw_con_irrig_col)) * dtime

        mask = trues(nc)
        CLM.withdraw_groundwater_irrigation!(wf, ws, mask, 1:nc, nlevsoi, dtime)

        liq_after = sum(ws.h2osoi_liq_col)
        wa_after  = sum(ws.wa_col)
        removed   = (liq_before - liq_after) + (wa_before - wa_after)

        @test removed ≈ flux_total atol = 1e-9
        # Per-layer/per-column exactness.
        @test ws.h2osoi_liq_col[1, 1] ≈ 20.0 - 0.01 * dtime
        @test ws.h2osoi_liq_col[2, 2] ≈ 20.0 - 0.004 * dtime
        @test ws.wa_col[1]            ≈ 5000.0 - 0.005 * dtime
        # Untouched columns unchanged.
        @test ws.h2osoi_liq_col[3, 1] ≈ 20.0
        @test ws.wa_col[3]            ≈ 5000.0
    end

    # --------------------------------------------------------------------
    # 2. set_floodc! — flood flux mapped to non-wall columns
    # --------------------------------------------------------------------
    @testset "set_floodc! maps gridcell flood flux" begin
        nc = 3
        ng = 2
        qflx_floodc  = zeros(nc)
        qflx_floodg  = [0.002, 0.007]
        col_gridcell = [1, 1, 2]
        # column 2 is a sunwall -> should get zero flood
        col_itype    = [1, CLM.ICOL_SUNWALL, 1]
        mask         = trues(nc)

        CLM.set_floodc!(qflx_floodc, qflx_floodg, col_gridcell, col_itype, mask, 1:nc)

        @test qflx_floodc[1] ≈ qflx_floodg[1]
        @test qflx_floodc[2] ≈ 0.0
        @test qflx_floodc[3] ≈ qflx_floodg[2]
    end

    # --------------------------------------------------------------------
    # 3. route_infiltration_excess! — partition conserves the input flux
    # --------------------------------------------------------------------
    @testset "route_infiltration_excess! conserves (h2osfc on/off)" begin
        nc = 2
        # vegetated (soil) columns so the istsoil/istcrop branch is taken
        col_landunit = [1, 1]
        lun_itype    = [CLM.ISTSOIL]

        for h2osfcflag in (1, 0)
            wfb = CLM.WaterFluxBulkData()
            CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
            wfb.qflx_in_soil_col            .= [0.03, 0.05]
            wfb.qflx_top_soil_to_h2osfc_col .= [0.01, 0.02]
            wfb.qflx_infl_excess_col        .= [0.004, 0.006]

            sh = CLM.SoilHydrologyData()
            CLM.soilhydrology_init!(sh, nc)
            sh.h2osfcflag = h2osfcflag

            mask = trues(nc)
            CLM.route_infiltration_excess!(wfb, sh, col_landunit, lun_itype,
                                           mask, 1:nc)

            for c in 1:nc
                # qflx_in_soil_limited always loses exactly the excess.
                @test wfb.qflx_in_soil_limited_col[c] ≈
                      wfb.qflx_in_soil_col[c] - wfb.qflx_infl_excess_col[c]

                # The excess plus the original h2osfc input must be fully
                # accounted for between the h2osfc pool and surface runoff,
                # with no creation/loss.
                routed = wfb.qflx_in_h2osfc_col[c] + wfb.qflx_infl_excess_surf_col[c]
                @test routed ≈ wfb.qflx_top_soil_to_h2osfc_col[c] +
                               wfb.qflx_infl_excess_col[c]

                if h2osfcflag != 0
                    @test wfb.qflx_infl_excess_surf_col[c] ≈ 0.0
                else
                    @test wfb.qflx_infl_excess_surf_col[c] ≈
                          wfb.qflx_infl_excess_col[c]
                end
            end
        end
    end

    # --------------------------------------------------------------------
    # 4. perched_lateral_flow! — per-column soil-storage balance
    # --------------------------------------------------------------------
    @testset "perched_lateral_flow! soil-storage balance" begin
        nc = 2
        col = build_columns(nc)
        ss  = build_soilstate(nc)

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.frost_table_col .= 0.5
        sh.zwt_col         .= 1.0
        sh.zwt_perched_col .= 0.3

        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, nc, 1, 1)
        ws.h2osoi_liq_col = fill(10.0, nc, nlevtot)
        ws.h2osoi_ice_col = fill(0.0, nc, nlevtot)
        ws.stream_water_volume_lun = fill(0.0, 1)

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 1)
        lun.stream_channel_length .= 100.0
        lun.stream_channel_width  .= 5.0
        lun.stream_channel_depth  .= 1.0

        tdepth    = fill(0.0, 1)
        tdepthmax = fill(1.0, 1)
        mask      = trues(nc)

        liq_before = copy(ws.h2osoi_liq_col)
        CLM.perched_lateral_flow!(sh, ss, ws, wfb, col, lun,
                                  tdepth, tdepthmax, mask, 1:nc, nlevsoi, dtime)

        for c in 1:nc
            qd = wfb.wf.qflx_drain_perched_col[c]
            @test isfinite(qd)
            # Non-hillslope draining column should lose water (qd >= 0).
            @test qd >= -1e-12
            # Water actually removed from soil storage must equal the
            # reported (residual-adjusted) drainage flux * dtime.
            removed = sum(liq_before[c, :]) - sum(ws.h2osoi_liq_col[c, :])
            @test removed ≈ qd * dtime atol = 1e-7
            @test removed >= -1e-9   # never adds water
        end
    end

    # --------------------------------------------------------------------
    # 5. subsurface_lateral_flow! — non-hillslope global non-creation
    # --------------------------------------------------------------------
    @testset "subsurface_lateral_flow! does not create water" begin
        nc = 2
        col = build_columns(nc)
        ss  = build_soilstate(nc)

        sh = CLM.SoilHydrologyData()
        CLM.soilhydrology_init!(sh, nc)
        sh.zwt_col         .= 0.5
        sh.frost_table_col .= 1.0

        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, nc, 1, 1)
        ws.h2osoi_liq_col = fill(10.0, nc, nlevtot)
        ws.h2osoi_ice_col = fill(0.0, nc, nlevtot)
        ws.h2osfc_col    .= 0.0
        ws.stream_water_volume_lun = fill(0.0, 1)

        wfb = CLM.WaterFluxBulkData()
        CLM.waterfluxbulk_init!(wfb, nc, nc, 1, 1)
        wfb.wf.qflx_snwcp_liq_col .= 0.0

        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, 1)
        lun.stream_channel_length .= 100.0
        lun.stream_channel_width  .= 5.0
        lun.stream_channel_depth  .= 1.0
        lun.stream_channel_number .= 1.0
        lun.urbpoi .= false

        grc_area  = fill(100.0, 1)
        tdepth    = fill(0.0, 1)
        tdepthmax = fill(1.0, 1)
        mask_hydro = trues(nc)
        mask_urban = falses(nc)

        liq_before = sum(ws.h2osoi_liq_col)
        CLM.subsurface_lateral_flow!(sh, ss, ws, wfb, col, lun,
                                     tdepth, tdepthmax, grc_area,
                                     mask_hydro, mask_urban, 1:nc, nlevsoi, dtime)
        liq_after = sum(ws.h2osoi_liq_col)

        for c in 1:nc
            @test isfinite(wfb.wf.qflx_drain_col[c])
            @test isfinite(wfb.wf.qflx_latflow_out_col[c])
            @test wfb.wf.qflx_latflow_out_col[c] >= -1e-12
        end
        # Baseflow leaves the column: total soil liquid may only decrease.
        @test liq_after <= liq_before + 1e-9
    end
end
