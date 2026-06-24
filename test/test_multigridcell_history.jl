# Tests for the multi-gridcell 2D (lon×lat) history output mode in
# src/infrastructure/history_io.jl. Mirrors the I/O test style of
# test_history_io.jl: a lightweight stand-in `inst` carrying patch/column
# subgrid index+weight maps, written to NetCDF, re-read, and asserted.
#
# With `dov2xy=true` + a g→(i,j) grid map, fields must land on (lon,lat[,lev],
# time) dims with each gridcell's value at its (i,j) and non-active cells at the
# fill value. The 1D-gridcell-dim mode (no map) and single-gridcell path must be
# unchanged.

using Test
using NCDatasets

# Subgrid map stand-in (matches history_io's _hist_subgrid_to_gridcell_maps:
# needs .gridcell + .wtgcell). Defined at top level (not inside @testset).
Base.@kwdef mutable struct _MgMap
    gridcell::Vector{Int}     = Int[]
    wtgcell::Vector{Float64}  = Float64[]
end
Base.@kwdef mutable struct _MgInst
    g_vals::Vector{Float64} = Float64[]                     # gridcell-level 1D
    p_vals::Vector{Float64} = Float64[]                     # patch-level 1D
    g2d::Matrix{Float64}    = Matrix{Float64}(undef, 0, 0)  # gridcell 2D (ng,nlev)
    patch    = nothing
    column   = nothing
    landunit = nothing
end

@testset "multigridcell_history (2D lon×lat)" begin
    # A 2×3 lon×lat grid (nlon=2, nlat=3 → 6 cells). We populate only 4 of the
    # 6 cells with active gridcells; the other 2 must come back as fill value.
    #   g : (i,j)        lon/lat (degrees)
    #   1 : (1,1)        lon=10, lat=40
    #   2 : (2,1)        lon=20, lat=40
    #   3 : (1,2)        lon=10, lat=50
    #   4 : (2,3)        lon=20, lat=60
    # missing: (2,2)=(lon=20,lat=50) and (1,3)=(lon=10,lat=60).
    ng   = 4
    lon  = [10.0, 20.0, 10.0, 20.0]
    lat  = [40.0, 40.0, 50.0, 60.0]
    g2i  = [1, 2, 1, 2]
    g2j  = [1, 1, 2, 3]
    nlon, nlat = 2, 3
    FILL = Float32(CLM.SPVAL)

    @testset "derive g→(i,j) map from lon/lat" begin
        m = CLM.hist_grid_map_from_lonlat(lon, lat)
        @test m.nlon == nlon
        @test m.nlat == nlat
        @test m.lon1d ≈ [10.0, 20.0]
        @test m.lat1d ≈ [40.0, 50.0, 60.0]
        @test m.g2i == g2i
        @test m.g2j == g2j
    end

    # _hist_ngridcell reads inst.{patch,column,landunit}.gridcell — every inst
    # needs all three maps (gridcell-level fields infer ng from these).
    gmap_for_ng() = _MgMap(gridcell = collect(1:ng), wtgcell = ones(ng))

    @testset "1D gridcell field → (lon,lat,time) + fill" begin
        inst = _MgInst(g_vals = zeros(ng),
                       patch = gmap_for_ng(),
                       column = _MgMap(), landunit = _MgMap())
        tape = CLM.HistoryTape(dov2xy = true)
        CLM.hist_set_grid_map!(tape; lon = lon, lat = lat)
        @test tape.xy2d
        @test tape.nlon == nlon && tape.nlat == nlat

        CLM.hist_addfld!(tape, "GFLD", "u", i -> i.g_vals; level = "gridcell")

        gv = [11.0, 22.0, 33.0, 44.0]
        inst.g_vals .= gv
        CLM.hist_accumulate!(tape, inst)

        mktempdir() do dir
            fn = joinpath(dir, "mg_1d.nc")
            CLM.hist_write!(tape, fn; time = 1.0, inst = inst)
            NCDataset(fn, "r") do ds
                @test ds.dim["lon"] == nlon
                @test ds.dim["lat"] == nlat
                @test ds.dim["time"] == 1
                @test !haskey(ds.dim, "gridcell")
                @test Array(ds["lon"]) ≈ [10.0, 20.0]
                @test Array(ds["lat"]) ≈ [40.0, 50.0, 60.0]

                rd = Array(ds["GFLD"].var)          # raw (no _FillValue masking)
                @test size(rd) == (nlon, nlat, 1)
                grid = rd[:, :, 1]
                for g in 1:ng
                    @test grid[g2i[g], g2j[g]] ≈ gv[g]
                end
                @test grid[2, 2] == FILL
                @test grid[1, 3] == FILL
                @test ds["GFLD"].attrib["subgrid_level"] == "gridcell"
            end
        end
    end

    @testset "patch field area-weight aggregated → (lon,lat,time)" begin
        # 4 patches, one per gridcell (wt 1.0 each) so the gridcell aggregate is
        # the patch value itself — keeps the 2D-placement assertion crisp.
        inst = _MgInst(
            p_vals = zeros(ng),
            patch  = _MgMap(gridcell = [1, 2, 3, 4], wtgcell = ones(4)),
            column = _MgMap(gridcell = Int[],        wtgcell = Float64[]),
            landunit = _MgMap(gridcell = Int[],      wtgcell = Float64[]),
        )
        tape = CLM.HistoryTape(dov2xy = true)
        # explicit g2i/g2j + axis coords path
        CLM.hist_set_grid_map!(tape; g2i = g2i, g2j = g2j, nlon = nlon, nlat = nlat,
                               lon1d = [10.0, 20.0], lat1d = [40.0, 50.0, 60.0])
        CLM.hist_addfld!(tape, "PFLD", "u", i -> i.p_vals; level = "patch")

        pv = [5.0, 6.0, 7.0, 8.0]
        inst.p_vals .= pv
        CLM.hist_accumulate!(tape, inst)

        mktempdir() do dir
            fn = joinpath(dir, "mg_patch.nc")
            CLM.hist_write!(tape, fn; time = 2.0, inst = inst)
            NCDataset(fn, "r") do ds
                @test ds.dim["lon"] == nlon
                @test ds.dim["lat"] == nlat
                grid = Array(ds["PFLD"].var)[:, :, 1]
                for g in 1:ng
                    @test grid[g2i[g], g2j[g]] ≈ pv[g]
                end
                @test grid[2, 2] == FILL
                @test grid[1, 3] == FILL
            end
        end
    end

    @testset "levsoi field → (lon,lat,levsoi,time)" begin
        nlev = 2
        inst = _MgInst(g2d = zeros(ng, nlev),
                       patch = gmap_for_ng(),
                       column = _MgMap(), landunit = _MgMap())
        tape = CLM.HistoryTape(dov2xy = true)
        CLM.hist_set_grid_map!(tape; g2i = g2i, g2j = g2j, nlon = nlon, nlat = nlat)
        CLM.hist_addfld!(tape, "TSOI", "K", i -> i.g2d;
                         level = "gridcell", nlev = nlev, levdim = "levsoi")

        gvals = [281.0 301.0;
                 282.0 302.0;
                 283.0 303.0;
                 284.0 304.0]
        inst.g2d .= gvals
        CLM.hist_accumulate!(tape, inst)

        mktempdir() do dir
            fn = joinpath(dir, "mg_levsoi.nc")
            CLM.hist_write!(tape, fn; time = 3.0, inst = inst)
            NCDataset(fn, "r") do ds
                @test ds.dim["lon"] == nlon
                @test ds.dim["lat"] == nlat
                @test ds.dim["levsoi"] == nlev
                @test ds.dim["time"] == 1
                rd = Array(ds["TSOI"].var)         # raw (lon, lat, levsoi, time)
                @test size(rd) == (nlon, nlat, nlev, 1)
                for g in 1:ng, ell in 1:nlev
                    @test rd[g2i[g], g2j[g], ell, 1] ≈ gvals[g, ell]
                end
                @test rd[2, 2, 1, 1] == FILL
                @test rd[2, 2, 2, 1] == FILL
                @test rd[1, 3, 1, 1] == FILL
            end
        end
    end

    @testset "1D gridcell-dim mode unchanged when xy2d off" begin
        # dov2xy true but NO grid map → falls back to the flat `gridcell` dim
        # exactly as before this change.
        inst = _MgInst(g_vals = zeros(ng),
                       patch = gmap_for_ng(),
                       column = _MgMap(), landunit = _MgMap())
        tape = CLM.HistoryTape(dov2xy = true)   # xy2d defaults false
        @test tape.xy2d == false
        CLM.hist_addfld!(tape, "GFLD", "u", i -> i.g_vals; level = "gridcell")
        gv = [11.0, 22.0, 33.0, 44.0]
        inst.g_vals .= gv
        CLM.hist_accumulate!(tape, inst)
        mktempdir() do dir
            fn = joinpath(dir, "mg_flat.nc")
            CLM.hist_write!(tape, fn; time = 1.0, inst = inst)
            NCDataset(fn, "r") do ds
                @test ds.dim["gridcell"] == ng
                @test !haskey(ds.dim, "lon")
                @test !haskey(ds.dim, "lat")
                @test Array(ds["GFLD"])[:, 1] ≈ gv
            end
        end
    end

    @testset "single-gridcell path unchanged" begin
        # ng=1: xy2d on but a 1×1 grid still round-trips the single value.
        inst = _MgInst(g_vals = zeros(1),
                       patch = _MgMap(gridcell = [1], wtgcell = [1.0]),
                       column = _MgMap(), landunit = _MgMap())
        tape = CLM.HistoryTape(dov2xy = true)
        CLM.hist_set_grid_map!(tape; lon = [12.5], lat = [55.0])
        @test tape.nlon == 1 && tape.nlat == 1
        CLM.hist_addfld!(tape, "GFLD", "u", i -> i.g_vals; level = "gridcell")
        inst.g_vals .= [42.0]
        CLM.hist_accumulate!(tape, inst)
        mktempdir() do dir
            fn = joinpath(dir, "mg_single.nc")
            CLM.hist_write!(tape, fn; time = 0.0, inst = inst)
            NCDataset(fn, "r") do ds
                @test ds.dim["lon"] == 1
                @test ds.dim["lat"] == 1
                @test Array(ds["GFLD"])[1, 1, 1] ≈ 42.0
            end
        end
    end
end
