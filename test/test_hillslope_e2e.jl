# =========================================================================
# Hillslope hydrology — END-TO-END on a SYNTHETIC catena
# =========================================================================
#
# Exercises the full opt-in (`use_hillslope`) plumbing added on this branch:
#   surfdata reader (surfrd_hillslope!)  ->  count_subgrid_elements (hillslope)
#   ->  initGridCells! multi-column builder  ->  init_hillslope_columns!
#   ->  subsurface_lateral_flow! (hillslope branch)  ->  stream outflow/update.
#
# The geometry is a CLEARLY SYNTHETIC test fixture: 1 gridcell, 1 hillslope,
# 4 columns (upland -> lowland) with monotone elevation/distance and equal,
# mass-consistent areas summing to the landunit. There is NO CTSM-produced
# hillslope fsurdat in this environment and NO Fortran hillslope run to compare
# against, so this validates PORT SELF-CONSISTENCY (conservation) and matches
# an INDEPENDENT plain-scalar oracle of the ported math — it does NOT claim
# CTSM parity.
#
# Gates:
#   1. Builder produces a valid 4-column catena (cold/colu chain to the stream).
#   2. Inter-column lateral flux + volumetric discharge match an independent
#      scalar re-implementation of subsurface_lateral_flow! (rtol ~1e-8).
#   3. Water routing between columns conserves exactly (inflow == uphill outflow).
#   4. Stream-water budget closes (hillslope_update_stream_water!).
#   5. DEFAULT (use_hillslope=false) natveg build is single-column, unchanged.
# =========================================================================

@testset "Hillslope Hydrology END-TO-END (synthetic catena)" begin
    using NCDatasets

    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    CLM.varcon_init!()
    nlevsoi = CLM.varpar.nlevsoi
    nlevsno = CLM.varpar.nlevsno

    # ---------------------------------------------------------------------
    # Synthetic hillslope surfdata fixture (SYNTHETIC — not a real fsurdat)
    # ---------------------------------------------------------------------
    # 4-column catena: col 1 = uppermost/highest/farthest, col 4 = lowland
    # (nearest the stream). Downhill chain 1->2->3->4->stream.
    scratch = get(ENV, "CLM_SCRATCH", mktempdir())
    fixture = joinpath(scratch, "SYNTHETIC_hillslope_surfdata.nc")

    ncol      = 4
    dcol      = Int32[2, 3, 4, -999]          # downhill_column_index (-999 = to stream)
    hill_elev = [30.0, 20.0, 10.0, 1.0]       # m, monotone decreasing downhill
    hill_dist = [300.0, 200.0, 100.0, 10.0]   # m, monotone decreasing downhill
    hill_area = [2500.0, 2500.0, 2500.0, 2500.0]  # m2, sum = 10000
    hill_wid  = [50.0, 50.0, 50.0, 50.0]      # m
    hill_slp  = [0.10, 0.10, 0.10, 0.10]      # m/m

    isfile(fixture) && rm(fixture)
    NCDataset(fixture, "c") do ds
        defDim(ds, "gridcell", 1)
        defDim(ds, "nhillslope", 1)
        defDim(ds, "nmaxhillcol", ncol)
        defDim(ds, "natpft", 15)

        # --- base surfdata: one natural-veg gridcell, single C3-grass PFT ---
        defVar(ds, "PCT_NATVEG", Float64, ("gridcell",))[:] = [100.0]
        defVar(ds, "PCT_CROP",   Float64, ("gridcell",))[:] = [0.0]
        pnp = defVar(ds, "PCT_NAT_PFT", Float64, ("natpft", "gridcell"))
        natpft = zeros(15, 1)
        natpft[14, 1] = 100.0             # PFT index 13 (m=14): C3 non-arctic grass
        pnp[:, :] = natpft

        # --- hillslope geometry ---
        defVar(ds, "nhillcolumns", Int32, ("gridcell",))[:] = Int32[ncol]
        defVar(ds, "pct_hillslope", Float64, ("nhillslope", "gridcell"))[:, :] = reshape([100.0], 1, 1)
        defVar(ds, "hillslope_index",       Int32, ("nmaxhillcol", "gridcell"))[:, :] = reshape(Int32[1,1,1,1], ncol, 1)
        defVar(ds, "column_index",          Int32, ("nmaxhillcol", "gridcell"))[:, :] = reshape(Int32[1,2,3,4], ncol, 1)
        defVar(ds, "downhill_column_index", Int32, ("nmaxhillcol", "gridcell"))[:, :] = reshape(dcol, ncol, 1)
        defVar(ds, "hillslope_slope",     Float64, ("nmaxhillcol", "gridcell"))[:, :] = reshape(hill_slp, ncol, 1)
        defVar(ds, "hillslope_aspect",    Float64, ("nmaxhillcol", "gridcell"))[:, :] = reshape(zeros(ncol), ncol, 1)
        defVar(ds, "hillslope_area",      Float64, ("nmaxhillcol", "gridcell"))[:, :] = reshape(hill_area, ncol, 1)
        defVar(ds, "hillslope_distance",  Float64, ("nmaxhillcol", "gridcell"))[:, :] = reshape(hill_dist, ncol, 1)
        defVar(ds, "hillslope_width",     Float64, ("nmaxhillcol", "gridcell"))[:, :] = reshape(hill_wid, ncol, 1)
        defVar(ds, "hillslope_elevation", Float64, ("nmaxhillcol", "gridcell"))[:, :] = reshape(hill_elev, ncol, 1)

        # --- stream channel geometry (routing) ---
        defVar(ds, "hillslope_stream_depth", Float64, ("gridcell",))[:] = [1.0]
        defVar(ds, "hillslope_stream_width", Float64, ("gridcell",))[:] = [3.0]
        defVar(ds, "hillslope_stream_slope", Float64, ("gridcell",))[:] = [0.02]
    end

    # ---------------------------------------------------------------------
    # Save + restore global control flags / method Refs around the testset.
    # ---------------------------------------------------------------------
    saved = (
        use_hillslope         = CLM.varctl.use_hillslope,
        use_hillslope_routing = CLM.varctl.use_hillslope_routing,
        nhillslope            = CLM.varctl.nhillslope,
        max_columns_hillslope = CLM.varctl.max_columns_hillslope,
        head                  = CLM.HEAD_GRADIENT_METHOD[],
        trans                 = CLM.TRANSMISSIVITY_METHOD[],
    )

    # Helper: build the full subgrid hierarchy from `fixture` under the current
    # varctl flags. Returns the instance tuple.
    function build_subgrid()
        surf = CLM.SurfaceInputData()
        CLM.surfrd_get_data!(surf, 1, 1, fixture)
        ng = 1
        (nl, nc, np) = CLM.count_subgrid_elements(surf, ng)
        bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl,
                                begc=1, endc=nc, begp=1, endp=np,
                                level=CLM.BOUNDS_LEVEL_CLUMP, clump_index=1)
        grc = CLM.GridcellData(); CLM.gridcell_init!(grc, ng)
        lun = CLM.LandunitData(); CLM.landunit_init!(lun, nl)
        col = CLM.ColumnData();   CLM.column_init!(col, nc)
        pch = CLM.PatchData();    CLM.patch_init!(pch, np)
        CLM.initGridCells!(bounds, surf, grc, lun, col, pch; lat=50.0, lon=250.0, area=1.0)
        CLM.set_active!(bounds, lun, col, pch)
        CLM.init_hillslope_columns!(bounds, surf, grc, lun, col, pch)
        return (surf=surf, bounds=bounds, grc=grc, lun=lun, col=col, pch=pch,
                nl=nl, nc=nc, np=np)
    end

    try
        # -----------------------------------------------------------------
        # 5 (checked first): DEFAULT byte-identical single-column build.
        # -----------------------------------------------------------------
        @testset "default (use_hillslope=false) is single-column" begin
            CLM.varctl.use_hillslope = false
            CLM.varctl.use_hillslope_routing = false
            s = build_subgrid()
            # No hillslope read, no hillslope geometry retained.
            @test isempty(s.surf.nhillcolumns)
            # Exactly one natural-veg (ISTSOIL) column, weight 1.0, not hillslope.
            soil_cols = findall(==(CLM.ISTSOIL), s.lun.itype)
            @test length(soil_cols) == 1
            lsoil = soil_cols[1]
            ncol_soil = s.lun.colf[lsoil] - s.lun.coli[lsoil] + 1
            @test ncol_soil == 1
            c1 = s.lun.coli[lsoil]
            @test s.col.wtlunit[c1] == 1.0
            @test s.col.is_hillslope_column[c1] == false
            @test all(s.col.cold .== CLM.ISPVAL)
        end

        # -----------------------------------------------------------------
        # Hillslope ON: build the catena.
        # -----------------------------------------------------------------
        CLM.varctl.use_hillslope = true
        CLM.varctl.use_hillslope_routing = true
        CLM.varctl.nhillslope = 1
        CLM.varctl.max_columns_hillslope = ncol
        s = build_subgrid()
        bounds, grc, lun, col, pch = s.bounds, s.grc, s.lun, s.col, s.pch

        lsoil = findfirst(==(CLM.ISTSOIL), lun.itype)
        c0 = lun.coli[lsoil]                    # first catena column (uppermost)
        cs = c0 .+ (0:ncol-1)                    # the 4 catena columns, in order

        @testset "1. builder produced a valid catena" begin
            @test s.nc == ncol + 1               # 4 hillslope + 1 lake
            @test lun.colf[lsoil] - lun.coli[lsoil] + 1 == ncol
            @test all(col.is_hillslope_column[cs])
            # downhill chain: c0->c0+1->c0+2->c0+3->stream(ISPVAL)
            @test col.cold[cs[1]] == cs[2]
            @test col.cold[cs[2]] == cs[3]
            @test col.cold[cs[3]] == cs[4]
            @test col.cold[cs[4]] == CLM.ISPVAL
            # uphill chain is the mirror
            @test col.colu[cs[1]] == CLM.ISPVAL
            @test col.colu[cs[2]] == cs[1]
            @test col.colu[cs[3]] == cs[2]
            @test col.colu[cs[4]] == cs[3]
            # per-column geometry threaded from surfdata
            @test col.hill_elev[cs] ≈ hill_elev
            @test col.hill_distance[cs] ≈ hill_dist
            @test col.hill_area[cs] ≈ hill_area
            # weights recomputed from hillslope areas (equal-area -> 0.25 each, sum 1)
            @test sum(col.wtlunit[cs]) ≈ 1.0
            @test all(col.wtlunit[cs] .≈ 0.25)
            # stream channel geometry set from surfdata + geometry
            @test lun.stream_channel_width[lsoil] ≈ 3.0
            @test lun.stream_channel_depth[lsoil] ≈ 1.0
            @test lun.stream_channel_length[lsoil] > 0.0
            @test lun.stream_channel_number[lsoil] > 0.0
        end

        # -----------------------------------------------------------------
        # Overlay physical soil/water state on the catena (vertical grid +
        # soil hydraulics; the subgrid builder does not do initVertical).
        # -----------------------------------------------------------------
        nc = s.nc
        for j in 1:size(col.dz, 2)
            col.dz[:, j] .= 0.1
        end
        col.z  = zeros(nc, size(col.z, 2))
        col.zi = zeros(nc, size(col.zi, 2))
        for j in 1:size(col.z, 2);  col.z[:, j]  .= 0.05 + (j - 1) * 0.1;  end
        for j in 1:size(col.zi, 2); col.zi[:, j] .= j * 0.1;               end
        col.nbedrock .= nlevsoi
        col.topo_slope .= 5.0

        ss = CLM.SoilStateData()
        nlevgrnd = CLM.varpar.nlevgrnd
        ss.watsat_col       = fill(0.4,  nc, nlevgrnd)
        ss.bsw_col          = fill(5.0,  nc, nlevgrnd)
        ss.hksat_col        = fill(0.01, nc, nlevgrnd)
        ss.sucsat_col       = fill(100.0, nc, nlevgrnd)
        ss.eff_porosity_col = fill(0.35, nc, nlevgrnd)
        ss.hk_l_col         = fill(0.005, nc, nlevgrnd)

        sh = CLM.SoilHydrologyData(); CLM.soilhydrology_init!(sh, nc)
        sh.zwt_col         .= 0.5
        sh.frost_table_col .= 1.0

        nlevtot = nlevsoi + nlevsno
        ws = CLM.WaterStateData(); CLM.waterstate_init!(ws, nc, nc, s.nl, 1)
        ws.h2osoi_liq_col = fill(10.0, nc, nlevtot)
        ws.h2osoi_ice_col = fill(0.0,  nc, nlevtot)
        ws.h2osfc_col    .= 0.0
        ws.stream_water_volume_lun .= 0.0

        wfb = CLM.WaterFluxBulkData(); CLM.waterfluxbulk_init!(wfb, nc, nc, s.nl, 1)
        wfb.wf.qflx_snwcp_liq_col .= 0.0

        dtime = 1800.0
        tdepth    = fill(0.0, 1)
        tdepthmax = fill(1.0, 1)
        grc_area  = [grc.area[1]]
        mask_hydro = falses(nc); mask_hydro[cs] .= true    # hydrology on catena only
        mask_urban = falses(nc)

        # Clean, closed-form regime for the oracle: KINEMATIC head gradient
        # (head = hill_slope) + UNIFORM transmissivity.
        CLM.HEAD_GRADIENT_METHOD[]  = CLM.HEAD_GRADIENT_KINEMATIC
        CLM.TRANSMISSIVITY_METHOD[] = CLM.TRANSMISSIVITY_UNIFORM

        # Snapshot the flow INPUTS (the routine mutates zwt / h2osoi_ice / _liq),
        # so the oracle re-derives transmissivity from the same pre-flow state.
        zwt_pre = copy(sh.zwt_col)
        ice_pre = copy(ws.h2osoi_ice_col)

        liq_before = sum(ws.h2osoi_liq_col[cs, :])
        CLM.subsurface_lateral_flow!(sh, ss, ws, wfb, col, lun,
            tdepth, tdepthmax, grc_area, mask_hydro, mask_urban,
            1:nc, nlevsoi, dtime;
            use_hillslope_routing = true, nhillslope = 1)

        wf = wfb.wf

        # -----------------------------------------------------------------
        # Independent plain-scalar ORACLE of the hillslope lateral-flow math.
        # -----------------------------------------------------------------
        joff    = nlevsno
        joff_zi = nlevsno + 1
        eice    = CLM.soilhydrology_params.e_ice
        # jwt (layer above the water table) for each catena column
        function jwt_of(c)
            jw = nlevsoi
            for j in 1:nlevsoi
                if zwt_pre[c] <= col.zi[c, j + joff_zi]
                    jw = j - 1; break
                end
            end
            return jw
        end
        # ice impedance of layer j (matches the ported smooth_min primitive)
        function iceimp(c, j)
            vol_ice = CLM.smooth_min(ss.watsat_col[c, j],
                        ice_pre[c, j + joff] / (col.dz[c, j + joff] * CLM.DENICE);
                        k = CLM.SOIL_HYDRAULIC_K)
            icf = CLM.smooth_min(1.0, vol_ice / ss.watsat_col[c, j]; k = CLM.SOIL_HYDRAULIC_K)
            return 10.0^(-eice * icf)
        end
        # outflow volume (m3/s) from column c (UNIFORM transmissivity, KINEMATIC head)
        out_vol = Dict{Int,Float64}()
        for c in cs
            jw = jwt_of(c)
            hg = min(max(hill_slp[c - c0 + 1], -2.0), 2.0)   # KINEMATIC, capped
            zi_bed = col.zi[c, col.nbedrock[c] + joff_zi]
            transmis = 0.0
            if zwt_pre[c] <= zi_bed
                transmis = 1.0e-3 * iceimp(c, jw + 1) * ss.hksat_col[c, jw + 1] *
                           (zi_bed - zwt_pre[c])
            end
            out_vol[c] = transmis * col.hill_width[c] * hg
        end

        @testset "2. lateral flux matches independent oracle" begin
            for c in cs
                # qflx_latflow_out = 1e3 * out_vol / hill_area  (flux form)
                expect_out = 1.0e3 * out_vol[c] / col.hill_area[c]
                @test wf.qflx_latflow_out_col[c] ≈ expect_out rtol=1e-8
            end
            # volumetric discharge is exported only from the lowland column
            clow = cs[4]
            expect_disch = out_vol[clow] * (grc.area[1] * 1.0e6 * col.wtgcell[clow] / col.hill_area[clow])
            @test wf.volumetric_discharge_col[clow] ≈ expect_disch rtol=1e-8
            @test wf.volumetric_discharge_col[cs[1]] == 0.0    # non-lowland: no discharge
        end

        @testset "3. inter-column water routing conserves exactly" begin
            # uppermost column receives no inflow
            @test wf.qflx_latflow_in_col[cs[1]] == 0.0
            # every downhill column's inflow volume == its uphill's outflow volume
            for k in 1:ncol-1
                c   = cs[k]
                cd  = cs[k+1]
                in_vol_downhill  = wf.qflx_latflow_in_col[cd] * col.hill_area[cd] / 1.0e3
                out_vol_uphill   = wf.qflx_latflow_out_col[c] * col.hill_area[c]  / 1.0e3
                @test in_vol_downhill ≈ out_vol_uphill rtol=1e-10
                @test in_vol_downhill ≈ out_vol[c]     rtol=1e-8
            end
            # no water created: catena soil storage may only decrease
            liq_after = sum(ws.h2osoi_liq_col[cs, :])
            @test liq_after <= liq_before + 1e-9
        end

        # -----------------------------------------------------------------
        # 4. Stream water budget: drive the stream in/out-flow and confirm the
        #    volume bookkeeping closes against an independent recomputation.
        # -----------------------------------------------------------------
        @testset "4. stream-water budget closes" begin
            volumetric_streamflow = zeros(s.nl)
            stream_water_depth    = zeros(s.nl)
            ws.stream_water_volume_lun[lsoil] = 1000.0   # some standing stream water

            qflx_drain         = wf.qflx_drain_col
            qflx_drain_perched = wf.qflx_drain_perched_col
            qflx_surf          = zeros(nc)

            # 4a. discharge via Manning
            CLM.hillslope_stream_outflow!(ws.stream_water_volume_lun,
                volumetric_streamflow, lun, bounds.begl:bounds.endl, dtime)
            @test volumetric_streamflow[lsoil] >= 0.0
            @test volumetric_streamflow[lsoil] <= ws.stream_water_volume_lun[lsoil] / dtime + 1e-12

            # 4b. update stream volume — recompute the expected result independently
            vol0 = ws.stream_water_volume_lun[lsoil]
            inflow_vol = 0.0
            for c in cs
                col_area_m2 = grc.area[1] * 1.0e6 * col.wtgcell[c]
                inflow_vol += (qflx_drain[c] + qflx_drain_perched[c] + qflx_surf[c]) *
                              1.0e-3 * col_area_m2
            end
            inflow_vol *= dtime
            expected_vol = vol0 + inflow_vol - volumetric_streamflow[lsoil] * dtime

            CLM.hillslope_update_stream_water!(ws.stream_water_volume_lun,
                volumetric_streamflow, stream_water_depth,
                qflx_drain, qflx_drain_perched, qflx_surf,
                col, lun, grc, bounds.begl:bounds.endl, dtime)

            @test ws.stream_water_volume_lun[lsoil] ≈ expected_vol rtol=1e-10
            @test stream_water_depth[lsoil] >= 0.0
        end

    finally
        CLM.varctl.use_hillslope         = saved.use_hillslope
        CLM.varctl.use_hillslope_routing = saved.use_hillslope_routing
        CLM.varctl.nhillslope            = saved.nhillslope
        CLM.varctl.max_columns_hillslope = saved.max_columns_hillslope
        CLM.HEAD_GRADIENT_METHOD[]  = saved.head
        CLM.TRANSMISSIVITY_METHOD[] = saved.trans
        isfile(fixture) && rm(fixture; force=true)
    end
end
