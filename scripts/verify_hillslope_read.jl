#!/usr/bin/env julia
# Verify that CLM.jl's hillslope init path (init_hillslope!) consumes the
# synthetic hillslope_file: read the netcdf variables (exact CTSM names/dims),
# map the active gridcell's catena onto a minimal 1-landunit/4-column fixture,
# and run init_hillslope!. Assert the resulting column geometry & weights.
#
# Usage: julia --project=<CLM.jl> verify_hillslope_read.jl <hillslope_file.nc>

using NCDatasets
using CLM
using CLM: init_hillslope!, ColumnData, LandunitData, GridcellData, ISTSOIL, ISPVAL

const F = length(ARGS) >= 1 ? ARGS[1] : "/private/tmp/claude-501/hillslope_ne3np4_synthetic_c.nc"

ds = NCDataset(F, "r")
nhill  = ds.dim["nhillslope"]
nmax   = ds.dim["nmaxhillcol"]
nhc    = Array{Int}(ds["nhillcolumns"][:])
gc     = argmax(nhc)                       # active gridcell (the one with columns)
ncol   = nhc[gc]
@assert ncol == nmax "test fixture expects a full catena"

# file layout is (extra, gridcell); pull the active gridcell's row
pcth   = Array{Float64}(ds["pct_hillslope"][gc, :])            # length nhill
hndx   = Array{Int}(ds["hillslope_index"][gc, :])              # length nmax
cndx   = Array{Int}(ds["column_index"][gc, :])
dndx   = Array{Int}(ds["downhill_column_index"][gc, :])
hslp   = Array{Float64}(ds["hillslope_slope"][gc, :])
hasp   = Array{Float64}(ds["hillslope_aspect"][gc, :])
harea  = Array{Float64}(ds["hillslope_area"][gc, :])
hdist  = Array{Float64}(ds["hillslope_distance"][gc, :])
hwid   = Array{Float64}(ds["hillslope_width"][gc, :])
helev  = Array{Float64}(ds["hillslope_elevation"][gc, :])
sdep   = Float64(ds["hillslope_stream_depth"][gc])
swid   = Float64(ds["hillslope_stream_width"][gc])
sslp   = Float64(ds["hillslope_stream_slope"][gc])
close(ds)

println("read active gridcell $gc: ncol=$ncol nhill=$nhill")
println("  column_index          = ", cndx)
println("  downhill_column_index = ", dndx)
println("  hillslope_area        = ", harea)

# ---- minimal fixture: 1 gridcell, 1 istsoil landunit, ncol columns
FT = Float64
grc = GridcellData(); grc.area = FT[1000.0]                       # km2 (nominal)
lun = LandunitData()
lun.itype   = Int[ISTSOIL]
lun.wtgcell = FT[1.0]
lun.coli    = Int[1]
lun.colf    = Int[ncol]
lun.gridcell = Int[1]
lun.stream_channel_depth  = fill(FT(0), 1)
lun.stream_channel_width  = fill(FT(0), 1)
lun.stream_channel_length = fill(FT(0), 1)
lun.stream_channel_slope  = fill(FT(0), 1)
lun.stream_channel_number = fill(FT(0), 1)

col = ColumnData()
col.wtlunit           = fill(FT(NaN), ncol)
col.colu              = fill(ISPVAL, ncol)
col.cold              = fill(ISPVAL, ncol)
col.hillslope_ndx     = fill(ISPVAL, ncol)
col.hill_elev         = fill(FT(0), ncol)
col.hill_slope        = fill(FT(0), ncol)
col.hill_area         = fill(FT(0), ncol)
col.hill_width        = fill(FT(0), ncol)
col.hill_distance     = fill(FT(0), ncol)
col.hill_aspect       = fill(FT(0), ncol)
col.is_hillslope_column = fill(false, ncol)

# kwargs indexed by landunit l=1 (rows) x column-in-hillslope ci (cols)
ncolumns_hillslope = Int[ncol]
pct_hillslope = reshape(pcth, 1, nhill)         # (nlandunit, nhillslope)
hill_ndx  = reshape(hndx, 1, nmax)
col_ndx   = reshape(cndx, 1, nmax)
col_dndx  = reshape(dndx, 1, nmax)
hill_slope  = reshape(hslp, 1, nmax)
hill_aspect = reshape(hasp, 1, nmax)
hill_area   = reshape(harea, 1, nmax)
hill_dist   = reshape(hdist, 1, nmax)
hill_width  = reshape(hwid, 1, nmax)
hill_elev   = reshape(helev, 1, nmax)

init_hillslope!(col, lun, grc, 1:1, 1:ncol;
    nhillslope = nhill, max_columns_hillslope = nmax,
    ncolumns_hillslope = ncolumns_hillslope,
    pct_hillslope = pct_hillslope,
    hill_ndx = hill_ndx, col_ndx_in = col_ndx, col_dndx = col_dndx,
    hill_slope_in = hill_slope, hill_aspect_in = hill_aspect,
    hill_area_in = hill_area, hill_dist_in = hill_dist,
    hill_width_in = hill_width, hill_elev_in = hill_elev,
    stream_channel_depth = fill(sdep, 1),
    stream_channel_width = fill(swid, 1),
    stream_channel_slope = fill(sslp, 1),
    use_hillslope_routing = true)

println("\n-- init_hillslope! results --")
println("  is_hillslope_column = ", col.is_hillslope_column)
println("  hillslope_ndx       = ", col.hillslope_ndx)
println("  cold (downhill)     = ", col.cold)
println("  colu (uphill)       = ", col.colu)
println("  hill_area           = ", col.hill_area)
println("  wtlunit             = ", col.wtlunit, "  sum=", sum(col.wtlunit))
println("  stream depth/width/slope = ", lun.stream_channel_depth[1], " ", lun.stream_channel_width[1], " ", lun.stream_channel_slope[1])
println("  stream channel number    = ", lun.stream_channel_number[1])
println("  stream channel length    = ", lun.stream_channel_length[1])

# ---- assertions
@assert all(col.is_hillslope_column) "all cols should be flagged hillslope"
@assert col.cold[1] == ISPVAL "lowland col1 must have no downhill (ISPVAL)"
@assert col.cold[2] == 1 && col.cold[3] == 2 && col.cold[4] == 3 "downhill chain must map ci->ci-1"
@assert isapprox(sum(col.wtlunit), 1.0; atol=1e-12) "column weights must sum to 1"
@assert col.hill_area == harea "areas must round-trip"
@assert lun.stream_channel_depth[1] == sdep "stream depth must round-trip"
println("\nPORT READER OK: init_hillslope! consumed the synthetic hillslope_file and produced a consistent catena.")
