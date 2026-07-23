# =========================================================================
# Hillslope hydrology — FULL DRIVER multi-column timestep (end-to-end)
# =========================================================================
#
# Exercises the deferred surface from PR #302: the full clm_run!/clm_drv!
# timestep on a MULTI-COLUMN hillslope catena (N soil columns per ISTSOIL
# landunit), plus the runtime stream-routing wire:
#
#   clm_initialize! (hillslope_file read + catena build)
#     -> clm_drv! -> hydrology_drainage! (subsurface_lateral_flow! +
#        hillslope_stream_outflow! + hillslope_update_stream_water!)
#     -> lnd2atm exclusion of hillslope cols from qflx_{surf,drain,drain_perched}_grc
#     -> hillslope_streamflow_to_grc! -> balance_check! (streamflow term).
#
# Two parts:
#   (A) Always-run, in-memory: unit-conservation of the new streamflow
#       aggregation + the stream-water budget on a hand-built catena. No
#       external assets, deterministic, safe under --check-bounds=yes.
#   (B) Asset-gated: a real clm_run! multi-column timestep on the Bow base
#       surfdata + a synthetic 4-column catena + generated forcing, asserting
#       finite state everywhere, water conservation (errh2o_grc), and the
#       DEFAULT (use_hillslope=false) single-column / zero-streamflow invariant.
# =========================================================================

include(joinpath(@__DIR__, "testdata.jl"))
include(joinpath(@__DIR__, "generate_forcing.jl"))

@testset "Hillslope FULL-DRIVER multi-column timestep (E2E)" begin

    # ---------------------------------------------------------------------
    # (A) In-memory conservation of the new streamflow aggregation.
    # ---------------------------------------------------------------------
    @testset "A. hillslope_streamflow_to_grc! + stream-budget conservation" begin
        # Minimal 1-gridcell / 1-landunit / 4-column catena.
        ncol = 4
        grc = CLM.GridcellData(); CLM.gridcell_init!(grc, 1); grc.area[1] = 100.0  # km2
        lun = CLM.LandunitData(); CLM.landunit_init!(lun, 1)
        lun.itype[1] = CLM.ISTSOIL
        lun.gridcell[1] = 1
        lun.active[1] = true
        lun.coli[1] = 1; lun.colf[1] = ncol
        lun.stream_channel_length[1] = 300.0
        lun.stream_channel_width[1]  = 3.0
        lun.stream_channel_depth[1]  = 1.0
        lun.stream_channel_slope[1]  = 0.02
        lun.stream_channel_number[1] = 1.0

        col = CLM.ColumnData(); CLM.column_init!(col, ncol)
        for c in 1:ncol
            col.gridcell[c] = 1; col.landunit[c] = 1
            col.is_hillslope_column[c] = true
            col.active[c] = true
            col.wtgcell[c] = 0.25 / 1.0          # each column 1/4 of the landunit (lun wt 1)
        end

        # State/flux vectors.
        swv  = zeros(1)                            # stream_water_volume_lun (m3)
        vsf  = zeros(1)                            # volumetric_streamflow_lun (m3/s)
        swd  = zeros(1)                            # stream_water_depth_lun (m)
        qdr  = fill(1.0e-4, ncol)                  # qflx_drain_col (mm/s)
        qdrp = fill(2.0e-5, ncol)                  # qflx_drain_perched_col (mm/s)
        qsrf = fill(5.0e-4, ncol)                  # qflx_surf_col (mm/s)
        dtime = 1800.0

        # One routing step: outflow (from current, empty channel → 0) then update.
        CLM.hillslope_stream_outflow!(swv, vsf, lun, 1:1, dtime)
        @test vsf[1] == 0.0                        # empty channel discharges nothing
        CLM.hillslope_update_stream_water!(swv, vsf, swd, qdr, qdrp, qsrf,
                                           col, lun, grc, 1:1, dtime)

        # Independent expected inflow: Σ_c (qdr+qdrp+qsrf)*1e-3*(area*1e6*wt)*dt.
        inflow = 0.0
        for c in 1:ncol
            a = grc.area[1] * 1.0e6 * col.wtgcell[c]
            inflow += (qdr[c] + qdrp[c] + qsrf[c]) * 1.0e-3 * a * dtime
        end
        @test isapprox(swv[1], inflow; rtol=1e-12)         # channel gained exactly inflow
        @test swd[1] ≈ swv[1] / 300.0 / 3.0                # depth bookkeeping

        # Second step: now the channel is non-empty → it discharges.
        CLM.hillslope_stream_outflow!(swv, vsf, lun, 1:1, dtime)
        @test vsf[1] > 0.0
        vol_before = swv[1]
        CLM.hillslope_update_stream_water!(swv, vsf, swd, qdr, qdrp, qsrf,
                                           col, lun, grc, 1:1, dtime)
        # Stream budget closes: ΔV = inflow - streamflow*dt (mass conservation).
        @test isapprox(swv[1] - vol_before, inflow - vsf[1]*dtime; rtol=1e-10)

        # hillslope_streamflow_to_grc!: m3/s → mm/s, summed (not weighted) over lun.
        qsf_grc = zeros(1)
        CLM.hillslope_streamflow_to_grc!(qsf_grc, vsf, lun, grc, 1:1, 1:1)
        @test qsf_grc[1] ≈ vsf[1] * 1.0e3 / (grc.area[1] * 1.0e6)
        @test qsf_grc[1] > 0.0

        # Inactive landunit contributes nothing.
        lun.active[1] = false
        CLM.hillslope_streamflow_to_grc!(qsf_grc, vsf, lun, grc, 1:1, 1:1)
        @test qsf_grc[1] == 0.0
        lun.active[1] = true
    end

    # ---------------------------------------------------------------------
    # (B) Full clm_run! multi-column timestep (asset-gated on Bow base).
    # ---------------------------------------------------------------------
    fsurdat, paramfile = bow_params()
    snowopt = snicar_optics(); snowage = snicar_aging()
    _assets_ok = isfile(fsurdat) && isfile(paramfile) && isfile(snowopt) && isfile(snowage)
    if !_assets_ok
        testdata_missing("hillslope full-driver e2e", fsurdat, paramfile, snowopt, snowage)
    end

    # Save/restore the global hillslope control flags so this test never leaks a
    # `use_hillslope=true` varctl into downstream tests (defensive — the default
    # sub-run below already resets them when it executes).
    _saved_uh  = CLM.varctl.use_hillslope
    _saved_uhr = CLM.varctl.use_hillslope_routing
    try
    if _assets_ok
    # Synthetic 4-column catena fixture on the Bow single gridcell (SYNTHETIC —
    # idealized upland→lowland chain; no CTSM hillslope fsurdat needed).
    hillf = tempname() * "_hillslope.nc"
    ncol = 4
    NCDataset(hillf, "c") do ds
        defDim(ds, "lsmlon", 1); defDim(ds, "lsmlat", 1)
        defDim(ds, "nhillslope", 1); defDim(ds, "nmaxhillcol", ncol)
        defVar(ds, "LONGXY", Float64, ("lsmlon","lsmlat"))[:,:] = reshape([250.0],1,1)
        defVar(ds, "LATIXY", Float64, ("lsmlon","lsmlat"))[:,:] = reshape([51.0],1,1)
        defVar(ds, "nhillcolumns", Int32, ("lsmlon","lsmlat"))[:,:] = reshape(Int32[ncol],1,1)
        defVar(ds, "pct_hillslope", Float64, ("lsmlon","lsmlat","nhillslope"))[:,:,:] = reshape([100.0],1,1,1)
        _l(v) = reshape(v, 1, 1, ncol)
        defVar(ds, "hillslope_index",       Int32,   ("lsmlon","lsmlat","nmaxhillcol"))[:,:,:] = _l(Int32[1,1,1,1])
        defVar(ds, "column_index",          Int32,   ("lsmlon","lsmlat","nmaxhillcol"))[:,:,:] = _l(Int32[1,2,3,4])
        defVar(ds, "downhill_column_index", Int32,   ("lsmlon","lsmlat","nmaxhillcol"))[:,:,:] = _l(Int32[-999,1,2,3])
        defVar(ds, "hillslope_slope",     Float64, ("lsmlon","lsmlat","nmaxhillcol"))[:,:,:] = _l(fill(0.1,ncol))
        defVar(ds, "hillslope_aspect",    Float64, ("lsmlon","lsmlat","nmaxhillcol"))[:,:,:] = _l(zeros(ncol))
        defVar(ds, "hillslope_area",      Float64, ("lsmlon","lsmlat","nmaxhillcol"))[:,:,:] = _l([2500.0,2250.0,2000.0,1750.0])
        defVar(ds, "hillslope_distance",  Float64, ("lsmlon","lsmlat","nmaxhillcol"))[:,:,:] = _l([10.0,100.0,200.0,300.0])
        defVar(ds, "hillslope_width",     Float64, ("lsmlon","lsmlat","nmaxhillcol"))[:,:,:] = _l(fill(50.0,ncol))
        defVar(ds, "hillslope_elevation", Float64, ("lsmlon","lsmlat","nmaxhillcol"))[:,:,:] = _l([1.0,10.0,20.0,30.0])
        defVar(ds, "hillslope_stream_depth", Float64, ("lsmlon","lsmlat"))[:,:] = reshape([1.0],1,1)
        defVar(ds, "hillslope_stream_width", Float64, ("lsmlon","lsmlat"))[:,:] = reshape([3.0],1,1)
        defVar(ds, "hillslope_stream_slope", Float64, ("lsmlon","lsmlat"))[:,:] = reshape([0.02],1,1)
    end

    forcing = tempname() * "_forcing.nc"
    generate_forcing(forcing; start_date=DateTime(2000,1,1),
                     end_date=DateTime(2000,1,2), dtime=1800)

    nsteps = 6
    start_date = DateTime(2000,1,1,0,0,0)
    end_date   = start_date + Second(nsteps*1800)

    # Capture per-step diagnostics.
    caps = NamedTuple[]
    probe = function(inst, tm)
        wb = inst.water.waterbalancebulk_inst
        ws = inst.water.waterstatebulk_inst.ws
        sh = inst.soilhydrology
        push!(caps, (errg = isempty(wb.errh2o_grc) ? NaN : maximum(abs.(wb.errh2o_grc)),
                     finite = all(isfinite, ws.h2osoi_liq_col) &&
                              all(isfinite, sh.zwt_col) &&
                              all(isfinite, ws.stream_water_volume_lun),
                     swv = copy(ws.stream_water_volume_lun)))
    end

    @testset "B. clm_run! builds + runs the 4-column catena, conserves" begin
        inst = CLM.clm_run!(; fsurdat=fsurdat, paramfile=paramfile, fforcing=forcing,
            fhistory=tempname()*".nc", start_date=start_date, end_date=end_date,
            dtime=1800, use_cn=false,
            use_hillslope=true, use_hillslope_routing=true, hillslope_file=hillf,
            gridcell_area_km2=12259.269, use_bedrock=false, use_aquifer_layer=false,
            fsnowoptics=snowopt, fsnowaging=snowage, verbose=false, step_probe=probe)

        col = inst.column
        hcols = findall(col.is_hillslope_column)
        # Catena built: 4 connected soil columns in one ISTSOIL landunit.
        @test length(hcols) == ncol
        @test col.cold[hcols] == [CLM.ISPVAL, hcols[1], hcols[2], hcols[3]]
        @test col.colu[hcols] == [hcols[2], hcols[3], hcols[4], CLM.ISPVAL]

        # Ran all steps with finite state everywhere.
        @test length(caps) == nsteps
        @test all(c.finite for c in caps)

        # Water conservation: gridcell balance closes each step (past the cold-start
        # transient the residual is ~machine precision). Bound generously but non-
        # vacuously — a broken multi-column balance is orders of magnitude larger.
        @test maximum(c.errg for c in caps) < 1.0e-4      # mm H2O
        @test maximum(caps[end].errg) < 1.0e-6            # steady-state near machine eps

        # Stream routing ran without producing garbage. The stream volume is
        # forcing-dependent (the synthetic Bow winter forcing can be frozen/dry ⇒
        # zero drainage ⇒ empty channel), so assert finite + non-negative here; the
        # nonzero-inflow channel fill is proven deterministically in part A (and on
        # the real Aripuanã tropical forcing the channel fills to ~5e5 m³).
        @test all(isfinite, caps[end].swv)
        @test caps[end].swv[1] >= 0.0
    end

    @testset "B. DEFAULT (use_hillslope=false) is single-column, zero streamflow" begin
        caps2 = NamedTuple[]
        probe2 = function(inst, tm)
            ws = inst.water.waterstatebulk_inst.ws
            push!(caps2, (nhill = count(inst.column.is_hillslope_column),
                          swv = copy(ws.stream_water_volume_lun),
                          finite = all(isfinite, ws.h2osoi_liq_col)))
        end
        inst = CLM.clm_run!(; fsurdat=fsurdat, paramfile=paramfile, fforcing=forcing,
            fhistory=tempname()*".nc", start_date=start_date, end_date=end_date,
            dtime=1800, use_cn=false, use_bedrock=false, use_aquifer_layer=false,
            fsnowoptics=snowopt, fsnowaging=snowage, verbose=false, step_probe=probe2)
        # No hillslope columns; single ISTSOIL column; no stream water accumulates.
        @test all(c.nhill == 0 for c in caps2)
        @test all(all(iszero, c.swv) || isempty(c.swv) for c in caps2)
        @test all(c.finite for c in caps2)
        soil_cols = findall(==(CLM.ISTSOIL), inst.landunit.itype)
        @test length(soil_cols) == 1
    end
    end  # if _assets_ok
    finally
        CLM.varctl.use_hillslope = _saved_uh
        CLM.varctl.use_hillslope_routing = _saved_uhr
    end
end
