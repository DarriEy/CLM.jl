# ==========================================================================
# Multi-gridcell atmospheric-forcing spatial interpolation tests
#
# Verifies:
#   1. A 2D (lon,lat) atm forcing grid is mapped to land gridcells by
#      nearest-neighbour: each land gridcell reads the forcing of its nearest
#      atm cell (distinct per cell, correct mapping).
#   2. A single-point forcing file still broadcasts to all land gridcells,
#      byte-identical to the historical single-gridcell behavior.
#   3. Time interpolation (closest-time selection) still works per-gridcell.
# ==========================================================================

using NCDatasets
using Dates

# ---- helpers to build tiny synthetic forcing files -----------------------

"""
Write a synthetic 2D (lon,lat,time) forcing file. TBOT (and the other vars)
are made spatially distinct so the nearest-neighbour mapping is unambiguous:
TBOT[i,j,t] = 200 + 10*i + 100*j + t  (so every atm cell has a unique value).
"""
function _write_grid_forcing(path; lons, lats, times, ref=DateTime(2000,1,1))
    nlon = length(lons); nlat = length(lats); ntime = length(times)
    ds = NCDataset(path, "c")
    defDim(ds, "lon", nlon); defDim(ds, "lat", nlat); defDim(ds, "time", ntime)
    defVar(ds, "lon", Float64, ("lon",))
    defVar(ds, "lat", Float64, ("lat",))
    defVar(ds, "time", Float64, ("time",);
           attrib = Dict("units" => "days since 2000-01-01", "calendar" => "noleap"))
    for v in ("TBOT","PSRF","WIND","FLDS","FSDS","PRECTmms","QBOT")
        defVar(ds, v, Float64, ("lon","lat","time"))
    end
    ds["lon"][:] = lons
    ds["lat"][:] = lats
    for t in 1:ntime
        ds["time"][t] = Dates.value(times[t] - ref) / (1000*86400)
    end
    for t in 1:ntime, j in 1:nlat, i in 1:nlon
        base = 200.0 + 10.0*i + 100.0*j + (t-1)
        ds["TBOT"][i,j,t]     = base                # distinct per cell & time
        ds["PSRF"][i,j,t]     = 85000.0 + 10.0*i + 100.0*j
        ds["WIND"][i,j,t]     = 1.0 + i + j
        ds["FLDS"][i,j,t]     = 200.0 + i + j
        ds["FSDS"][i,j,t]     = 0.0
        ds["PRECTmms"][i,j,t] = 0.0
        ds["QBOT"][i,j,t]     = 0.003
    end
    close(ds)
    return path
end

"""Write a single-point (time-only) forcing file (the legacy layout)."""
function _write_point_forcing(path; times, ref=DateTime(2000,1,1))
    ntime = length(times)
    ds = NCDataset(path, "c")
    defDim(ds, "time", ntime)
    defVar(ds, "time", Float64, ("time",);
           attrib = Dict("units" => "days since 2000-01-01", "calendar" => "noleap"))
    for v in ("TBOT","PSRF","WIND","FLDS","FSDS","PRECTmms","QBOT")
        defVar(ds, v, Float64, ("time",))
    end
    for t in 1:ntime
        ds["time"][t]     = Dates.value(times[t] - ref) / (1000*86400)
        ds["TBOT"][t]     = 270.0 + t
        ds["PSRF"][t]     = 85000.0
        ds["WIND"][t]     = 3.0
        ds["FLDS"][t]     = 250.0
        ds["FSDS"][t]     = 0.0
        ds["PRECTmms"][t] = 0.0
        ds["QBOT"][t]     = 0.003
    end
    close(ds)
    return path
end

@testset "Multi-gridcell forcing spatial interpolation" begin

    tmpdir = mktempdir()
    times = collect(DateTime(2000,1,1):Hour(1):DateTime(2000,1,1,3))  # 4 steps

    @testset "2D atm grid → nearest-neighbour per land gridcell" begin
        # 3 lon × 2 lat atm grid
        lons = [-120.0, -110.0, -100.0]
        lats = [40.0, 50.0]
        fpath = _write_grid_forcing(joinpath(tmpdir, "grid.nc"); lons=lons, lats=lats, times=times)

        fr = CLM.ForcingReader()
        CLM.forcing_reader_init!(fr, fpath)
        @test !fr.is_point
        @test fr.nlon == 3 && fr.nlat == 2
        @test fr.n_atm == 6

        # 4 land gridcells, each near a distinct atm cell:
        #   g1 ~ (lon=-120, lat=40) → atm (i=1,j=1) → aidx 1
        #   g2 ~ (lon=-100, lat=40) → atm (i=3,j=1) → aidx 3
        #   g3 ~ (lon=-110, lat=50) → atm (i=2,j=2) → aidx 5
        #   g4 ~ (lon=-100, lat=50) → atm (i=3,j=2) → aidx 6
        ng = 4
        gc_lat = [40.3, 39.7, 49.8, 50.2]
        gc_lon = [-119.8, -100.4, -110.1, -99.5]

        a2l = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a2l, ng, ng, ng)

        CLM.read_forcing_step!(fr, a2l, DateTime(2000,1,1,0), ng, ng;
                               gridcell_latdeg=gc_lat, gridcell_londeg=gc_lon)

        # The map is built and cached
        @test fr.map_g2atm == [1, 3, 5, 6]

        # TBOT[i,j,t=1] = 200 + 10i + 100j + 0
        expected = [200.0 + 10*1 + 100*1,    # g1: i1 j1
                    200.0 + 10*3 + 100*1,    # g2: i3 j1
                    200.0 + 10*2 + 100*2,    # g3: i2 j2
                    200.0 + 10*3 + 100*2]    # g4: i3 j2
        @test a2l.forc_t_not_downscaled_grc[1:ng] ≈ expected

        # Distinct per cell
        @test length(unique(a2l.forc_t_not_downscaled_grc[1:ng])) == ng

        # PSRF also mapped distinctly (no spatial collapse to cell 1)
        psrf_exp = [85000.0 + 10*1 + 100*1, 85000.0 + 10*3 + 100*1,
                    85000.0 + 10*2 + 100*2, 85000.0 + 10*3 + 100*2]
        @test a2l.forc_pbot_not_downscaled_grc[1:ng] ≈ psrf_exp

        CLM.forcing_reader_close!(fr)
    end

    @testset "single-point file broadcasts (byte-identical to legacy)" begin
        fpath = _write_point_forcing(joinpath(tmpdir, "point.nc"); times=times)

        fr = CLM.ForcingReader()
        CLM.forcing_reader_init!(fr, fpath)
        @test fr.is_point            # no horizontal coords → point
        @test fr.n_atm == 1

        ng = 3
        gc_lat = [40.0, 50.0, 60.0]
        gc_lon = [-120.0, -110.0, -100.0]

        a2l = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a2l, ng, ng, ng)
        CLM.read_forcing_step!(fr, a2l, DateTime(2000,1,1,0), ng, ng;
                               gridcell_latdeg=gc_lat, gridcell_londeg=gc_lon)

        # All gridcells get the SAME (broadcast) value at t=1: TBOT=271
        @test all(a2l.forc_t_not_downscaled_grc[1:ng] .== 271.0)
        @test all(a2l.forc_pbot_not_downscaled_grc[1:ng] .== 85000.0)

        # Byte-identical to calling WITHOUT coords (the historical signature)
        a2l_legacy = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a2l_legacy, ng, ng, ng)
        fr2 = CLM.ForcingReader()
        CLM.forcing_reader_init!(fr2, fpath)
        CLM.read_forcing_step!(fr2, a2l_legacy, DateTime(2000,1,1,0), ng, ng)
        @test a2l.forc_t_not_downscaled_grc      == a2l_legacy.forc_t_not_downscaled_grc
        @test a2l.forc_pbot_not_downscaled_grc   == a2l_legacy.forc_pbot_not_downscaled_grc
        @test a2l.forc_vp_grc                    == a2l_legacy.forc_vp_grc
        @test a2l.forc_solad_not_downscaled_grc  == a2l_legacy.forc_solad_not_downscaled_grc

        CLM.forcing_reader_close!(fr)
        CLM.forcing_reader_close!(fr2)
    end

    @testset "per-gridcell time interpolation (closest-time) still works" begin
        lons = [-120.0, -100.0]; lats = [45.0]
        fpath = _write_grid_forcing(joinpath(tmpdir, "grid_t.nc"); lons=lons, lats=lats, times=times)

        fr = CLM.ForcingReader()
        CLM.forcing_reader_init!(fr, fpath)
        ng = 2
        gc_lat = [45.0, 45.0]
        gc_lon = [-119.0, -101.0]   # g1→i1, g2→i2
        a2l = CLM.Atm2LndData()
        CLM.atm2lnd_init!(a2l, ng, ng, ng)

        # target nearest time t=3 (index 3 → t-1=2): TBOT = 200 + 10i + 100 + 2
        CLM.read_forcing_step!(fr, a2l, DateTime(2000,1,1,2,5), ng, ng;
                               gridcell_latdeg=gc_lat, gridcell_londeg=gc_lon)
        @test fr.time_index == 3
        @test a2l.forc_t_not_downscaled_grc[1] ≈ 200.0 + 10*1 + 100*1 + 2
        @test a2l.forc_t_not_downscaled_grc[2] ≈ 200.0 + 10*2 + 100*1 + 2
        # still distinct per cell
        @test a2l.forc_t_not_downscaled_grc[1] != a2l.forc_t_not_downscaled_grc[2]

        CLM.forcing_reader_close!(fr)
    end

    rm(tmpdir; recursive=true, force=true)
end
