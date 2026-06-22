# Include the source file into the CLM module if not already present.
if !isdefined(CLM, :dyn_hwcontent_final!)
    Base.include(CLM, joinpath(@__DIR__, "..", "src", "infrastructure",
                               "dyn_cons_biogeophys.jl"))
end

@testset "dynConsBiogeophysMod" begin

    # Dimension parameters
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    CLM.varcon_init!()

    nlevsno        = CLM.varpar.nlevsno
    nlevgrnd       = CLM.varpar.nlevgrnd
    nlevurb        = CLM.varpar.nlevurb
    nlevmaxurbgrnd = CLM.varpar.nlevmaxurbgrnd
    nlevlak        = CLM.varpar.nlevlak
    nlevtot        = nlevsno + nlevmaxurbgrnd

    # -----------------------------------------------------------------------
    # Build a minimal single-gridcell, single-landunit (natural-veg soil),
    # `nc`-column world with known water/heat state. Columns are soil columns
    # (no urban, no lake, no snow layers); each maps 1:1 to its own landunit on
    # gridcell 1 by default unless `lun_types` overrides. All columns belong to
    # gridcell 1.
    # -----------------------------------------------------------------------
    function setup(nc; lun_types=nothing, wtgcell=nothing)
        bounds = CLM.BoundsType(begg=1, endg=1, begl=1, endl=nc, begc=1, endc=nc,
                                begp=1, endp=nc)

        col = CLM.ColumnData()
        CLM.column_init!(col, nc)
        lun = CLM.LandunitData()
        CLM.landunit_init!(lun, nc)
        grc = CLM.GridcellData()
        CLM.gridcell_init!(grc, 1)
        grc.active[1] = true

        # Weights: split gridcell evenly across the nc columns by default.
        wts = wtgcell !== nothing ? wtgcell : fill(1.0 / nc, nc)

        for c in 1:nc
            lt = lun_types !== nothing ? lun_types[c] : CLM.ISTSOIL
            # landunit c lives on gridcell 1
            lun.gridcell[c] = 1
            lun.itype[c]    = lt
            lun.active[c]   = true
            lun.wtgcell[c]  = wts[c]
            lun.urbpoi[c]   = false
            lun.lakpoi[c]   = (lt == CLM.ISTDLAK)
            lun.canyon_hwr[c] = 1.0

            # register this landunit in the gridcell landunit_indices table
            grc.landunit_indices[lt, 1] = c

            col.landunit[c]              = c
            col.gridcell[c]              = 1
            col.itype[c]                 = CLM.ISTSOIL   # soil column
            col.lun_itype[c]             = lt
            col.active[c]                = true
            col.hydrologically_active[c] = (lt == CLM.ISTSOIL)
            col.wtgcell[c]               = wts[c]
            col.wtlunit[c]               = 1.0
            col.snl[c]                   = 0
            for j in 1:nlevtot
                col.dz[c, j] = 0.1
            end
            for j in 1:nlevlak
                col.dz_lake[c, j] = 1.0
            end
        end

        # Water state
        ws = CLM.WaterStateData()
        CLM.waterstate_init!(ws, nc, nc, nc, 1)
        for c in 1:nc
            ws.h2osno_no_layers_col[c]   = 0.0
            ws.h2osfc_col[c]             = 0.0
            ws.wa_col[c]                 = CLM.AQUIFER_WATER_BASELINE
            ws.dynbal_baseline_liq_col[c] = 0.0
            ws.dynbal_baseline_ice_col[c] = 0.0
            for j in 1:nlevtot
                ws.h2osoi_liq_col[c, j] = 0.0
                ws.h2osoi_ice_col[c, j] = 0.0
                ws.excess_ice_col[c, j] = 0.0
            end
        end
        ws.aquifer_water_baseline = CLM.AQUIFER_WATER_BASELINE

        wdb = CLM.WaterDiagnosticBulkData()
        CLM.waterdiagnosticbulk_init!(wdb, nc, nc, nc, 1)
        for c in 1:nc
            wdb.total_plant_stored_h2o_col[c] = 0.0
        end

        wsb = CLM.WaterStateBulkData()
        CLM.waterstatebulk_init!(wsb, nc, nc, nc, 1)
        wsb.ws = ws

        temp = CLM.TemperatureData()
        CLM.temperature_init!(temp, nc, nc, nc, 1)
        for c in 1:nc
            temp.t_h2osfc_col[c] = CLM.TFRZ
            temp.dynbal_baseline_heat_col[c] = 0.0
            for j in 1:nlevtot
                temp.t_soisno_col[c, j] = CLM.TFRZ
            end
            for j in 1:nlevlak
                temp.t_lake_col[c, j] = CLM.TFRZ
            end
        end

        ss = CLM.SoilStateData()
        CLM.soilstate_init!(ss, nc, nc)
        for c in 1:nc
            for j in 1:nlevmaxurbgrnd
                if j <= nlevgrnd
                    ss.watsat_col[c, j] = 0.4
                    ss.csol_col[c, j]   = 2.0e6
                else
                    ss.watsat_col[c, j] = 0.0
                end
            end
        end

        ls = CLM.LakeStateData()
        CLM.lakestate_init!(ls, nc, nc)
        for c in 1:nc
            for j in 1:nlevlak
                ls.lake_icefrac_col[c, j] = 0.0
            end
        end

        up = CLM.UrbanParamsData()
        CLM.urbanparams_init!(up, nc; nlevurb=nlevurb)
        for l in 1:nc
            up.nlev_improad[l] = 2
            for j in 1:nlevurb
                up.cv_wall[l, j]    = 1.5e6
                up.cv_roof[l, j]    = 1.2e6
                up.cv_improad[l, j] = 1.8e6
            end
        end

        return bounds, col, lun, grc, ws, wdb, wsb, temp, ss, ls, up
    end

    # Build the non-lake / lake masks from landunit types.
    function masks(lun, col, nc)
        mask_nolakec = falses(nc)
        mask_lakec   = falses(nc)
        for c in 1:nc
            if lun.itype[col.landunit[c]] == CLM.ISTDLAK
                mask_lakec[c] = true
            else
                mask_nolakec[c] = true
            end
        end
        return mask_nolakec, mask_lakec
    end

    ctl = CLM.dyn_subgrid_control_init()           # default: zero-fluxes = false
    ctl_zero = CLM.dyn_subgrid_control_init(for_testing_zero_dynbal_fluxes=true)

    # ===================================================================
    @testset "state init" begin
        dynbal = CLM.dyn_cons_biogeophys_state_init(3)
        @test length(dynbal.liq1_grc) == 3
        @test all(dynbal.qflx_liq_dynbal_grc .== 0.0)
        @test all(dynbal.eflx_dynbal_grc .== 0.0)
    end

    # ===================================================================
    @testset "no weight/state change => zero dynbal fluxes" begin
        nc = 2
        bounds, col, lun, grc, ws, wdb, wsb, temp, ss, ls, up = setup(nc)
        mask_nolakec, mask_lakec = masks(lun, col, nc)

        # put some water + nonzero temperature into both columns
        for c in 1:nc
            for j in (nlevsno+1):(nlevsno+nlevgrnd)
                ws.h2osoi_liq_col[c, j] = 20.0
                ws.h2osoi_ice_col[c, j] = 5.0
                temp.t_soisno_col[c, j] = CLM.TFRZ + 8.0
            end
        end

        dynbal = CLM.dyn_cons_biogeophys_state_init(1)

        CLM.dyn_hwcontent_init!(dynbal, bounds, mask_nolakec, mask_lakec,
                                col, lun, up, ss, wsb, wdb, temp, ls)
        # NO change to weights or state
        CLM.dyn_hwcontent_final!(dynbal, bounds, mask_nolakec, mask_lakec,
                                 col, lun, up, ss, wsb, wdb, temp, ls, ctl)

        @test dynbal.liq2_grc[1] ≈ dynbal.liq1_grc[1] atol=1e-9
        @test dynbal.ice2_grc[1] ≈ dynbal.ice1_grc[1] atol=1e-9
        @test dynbal.heat2_grc[1] ≈ dynbal.heat1_grc[1] atol=1e-3
        @test dynbal.qflx_liq_dynbal_grc[1] ≈ 0.0 atol=1e-9
        @test dynbal.qflx_ice_dynbal_grc[1] ≈ 0.0 atol=1e-9
        @test dynbal.eflx_dynbal_grc[1] ≈ 0.0 atol=1e-3
    end

    # ===================================================================
    @testset "weight change => dynbal water flux = gridcell content delta" begin
        # Two soil columns with DIFFERENT water content. Shifting the gridcell
        # weight from the dry column to the wet column changes the gridcell-mean
        # water content; the dynbal flux must equal that change exactly.
        nc = 2
        bounds, col, lun, grc, ws, wdb, wsb, temp, ss, ls, up = setup(nc;
            wtgcell=[0.5, 0.5])
        mask_nolakec, mask_lakec = masks(lun, col, nc)

        # Column 1 = wet, column 2 = dry; both at a fixed temperature
        for j in (nlevsno+1):(nlevsno+nlevgrnd)
            ws.h2osoi_liq_col[1, j] = 30.0
            ws.h2osoi_liq_col[2, j] = 10.0
            ws.h2osoi_ice_col[1, j] = 4.0
            ws.h2osoi_ice_col[2, j] = 1.0
            temp.t_soisno_col[1, j] = CLM.TFRZ + 5.0
            temp.t_soisno_col[2, j] = CLM.TFRZ + 5.0
        end

        dynbal = CLM.dyn_cons_biogeophys_state_init(1)
        CLM.dyn_hwcontent_init!(dynbal, bounds, mask_nolakec, mask_lakec,
                                col, lun, up, ss, wsb, wdb, temp, ls)

        liq1 = dynbal.liq1_grc[1]
        ice1 = dynbal.ice1_grc[1]

        # Land-cover change: shift all weight onto the wet column 1.
        col.wtgcell[1] = 1.0; col.wtgcell[2] = 0.0
        lun.wtgcell[1] = 1.0; lun.wtgcell[2] = 0.0

        CLM.dyn_hwcontent_final!(dynbal, bounds, mask_nolakec, mask_lakec,
                                 col, lun, up, ss, wsb, wdb, temp, ls, ctl)

        liq2 = dynbal.liq2_grc[1]
        ice2 = dynbal.ice2_grc[1]

        # The gridcell content actually changed (sanity)
        @test liq2 > liq1
        @test ice2 > ice1

        # The dynbal flux equals the (after - before) gridcell content delta.
        @test dynbal.qflx_liq_dynbal_grc[1] ≈ (liq2 - liq1) atol=1e-9
        @test dynbal.qflx_ice_dynbal_grc[1] ≈ (ice2 - ice1) atol=1e-9

        # Cross-check the gridcell means against hand computation.
        # Column liquid content (soil only) = sum over nlevgrnd layers.
        col1_liq = 30.0 * nlevgrnd
        col2_liq = 10.0 * nlevgrnd
        @test liq1 ≈ 0.5 * col1_liq + 0.5 * col2_liq atol=1e-6
        @test liq2 ≈ col1_liq atol=1e-6
    end

    # ===================================================================
    @testset "weight change => dynbal energy flux tracks heat delta" begin
        nc = 2
        bounds, col, lun, grc, ws, wdb, wsb, temp, ss, ls, up = setup(nc;
            wtgcell=[0.5, 0.5])
        mask_nolakec, mask_lakec = masks(lun, col, nc)

        # Two columns with the SAME water but DIFFERENT temperature.
        for j in (nlevsno+1):(nlevsno+nlevgrnd)
            ws.h2osoi_liq_col[1, j] = 20.0
            ws.h2osoi_liq_col[2, j] = 20.0
            temp.t_soisno_col[1, j] = CLM.TFRZ + 12.0   # warm
            temp.t_soisno_col[2, j] = CLM.TFRZ + 2.0    # cool
        end

        dynbal = CLM.dyn_cons_biogeophys_state_init(1)
        CLM.dyn_hwcontent_init!(dynbal, bounds, mask_nolakec, mask_lakec,
                                col, lun, up, ss, wsb, wdb, temp, ls)

        # Shift all weight onto the warm column.
        col.wtgcell[1] = 1.0; col.wtgcell[2] = 0.0
        lun.wtgcell[1] = 1.0; lun.wtgcell[2] = 0.0

        CLM.dyn_hwcontent_final!(dynbal, bounds, mask_nolakec, mask_lakec,
                                 col, lun, up, ss, wsb, wdb, temp, ls, ctl)

        # Gridcell heat content rose (now all the warm column).
        @test dynbal.heat2_grc[1] > dynbal.heat1_grc[1]

        # Since water content is unchanged (delta_liq = 0), the heat adjustment
        # is a no-op and the energy flux equals the raw heat delta.
        @test dynbal.qflx_liq_dynbal_grc[1] ≈ 0.0 atol=1e-9
        @test dynbal.eflx_dynbal_grc[1] ≈
              (dynbal.heat2_grc[1] - dynbal.heat1_grc[1]) atol=1e-3
        @test dynbal.eflx_dynbal_grc[1] != 0.0
    end

    # ===================================================================
    @testset "get_for_testing_zero_dynbal_fluxes => all fluxes zeroed" begin
        nc = 2
        bounds, col, lun, grc, ws, wdb, wsb, temp, ss, ls, up = setup(nc;
            wtgcell=[0.5, 0.5])
        mask_nolakec, mask_lakec = masks(lun, col, nc)

        for j in (nlevsno+1):(nlevsno+nlevgrnd)
            ws.h2osoi_liq_col[1, j] = 30.0
            ws.h2osoi_liq_col[2, j] = 10.0
            temp.t_soisno_col[1, j] = CLM.TFRZ + 10.0
            temp.t_soisno_col[2, j] = CLM.TFRZ + 1.0
        end

        dynbal = CLM.dyn_cons_biogeophys_state_init(1)
        CLM.dyn_hwcontent_init!(dynbal, bounds, mask_nolakec, mask_lakec,
                                col, lun, up, ss, wsb, wdb, temp, ls)

        # Real land-cover change
        col.wtgcell[1] = 1.0; col.wtgcell[2] = 0.0
        lun.wtgcell[1] = 1.0; lun.wtgcell[2] = 0.0

        # ...but with zero-dynbal-fluxes requested
        CLM.dyn_hwcontent_final!(dynbal, bounds, mask_nolakec, mask_lakec,
                                 col, lun, up, ss, wsb, wdb, temp, ls, ctl_zero)

        # Content still differs, but every dynbal flux is forced to zero.
        @test dynbal.liq2_grc[1] != dynbal.liq1_grc[1]
        @test dynbal.qflx_liq_dynbal_grc[1] == 0.0
        @test dynbal.qflx_ice_dynbal_grc[1] == 0.0
        @test dynbal.eflx_dynbal_grc[1] == 0.0
    end

    # ===================================================================
    @testset "set_baselines: glacier baseline = ice-col minus natveg" begin
        # Gridcell 1 has a natural-veg soil landunit (col 1) and a glacier
        # landunit (col 2). The glacier water/heat baseline should be the
        # glacier column's soil content minus the natural-veg column content.
        nc = 2
        bounds, col, lun, grc, ws, wdb, wsb, temp, ss, ls, up = setup(nc;
            lun_types=[CLM.ISTSOIL, CLM.ISTICE], wtgcell=[0.6, 0.4])
        mask_nolakec, mask_lakec = masks(lun, col, nc)

        # natveg column (1) vs glacier column (2) soil liquid content
        for j in (nlevsno+1):(nlevsno+nlevgrnd)
            ws.h2osoi_liq_col[1, j] = 12.0   # natveg
            ws.h2osoi_liq_col[2, j] = 50.0   # glacier "virtual" soil
            ws.h2osoi_ice_col[1, j] = 3.0
            ws.h2osoi_ice_col[2, j] = 40.0
        end

        mask_icec = falses(nc); mask_icec[2] = true
        mask_lakec_bl = falses(nc)

        CLM.dyn_hwcontent_set_baselines!(bounds, mask_icec, mask_lakec_bl,
                                         col, lun, grc, up, ss, ls, wsb, wdb, temp;
                                         reset_all_baselines=true,
                                         reset_lake_baselines=false)

        # Expected: glacier baseline_liq = (glacier soil liq) - (natveg soil liq)
        natveg_liq  = 12.0 * nlevgrnd
        glacier_liq = 50.0 * nlevgrnd
        @test ws.dynbal_baseline_liq_col[2] ≈ (glacier_liq - natveg_liq) atol=1e-6

        natveg_ice  = 3.0 * nlevgrnd
        glacier_ice = 40.0 * nlevgrnd
        @test ws.dynbal_baseline_ice_col[2] ≈ (glacier_ice - natveg_ice) atol=1e-6

        # natveg column baseline untouched (still 0)
        @test ws.dynbal_baseline_liq_col[1] == 0.0
    end

end
