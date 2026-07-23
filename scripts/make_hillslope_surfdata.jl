#!/usr/bin/env julia
# =============================================================================
# make_hillslope_surfdata.jl
#
# Generate a CTSM-conformant HILLSLOPE data file (`hillslope_file`) for parity
# testing of CLM.jl's hillslope hydrology against Fortran CTSM.
#
# IMPORTANT — what this file is:
#   * REAL BASE GRID / DOMAIN: longitude & latitude (LONGXY/LATIXY) are copied
#     verbatim from a real CTSM surface dataset, so the file passes CTSM's
#     `check_domain_attributes` against the same grid's mesh.
#   * SYNTHETIC GEOMETRY: the hillslope catena itself (1 hillslope, 4 columns,
#     upland->lowland) is an idealized, mass-consistent synthetic geometry. It
#     is NOT derived from a DEM. This is a legitimate CTSM-faithful *test* input
#     (idealized geometry on a real base), clearly labelled as synthetic.
#
# CTSM reads these variables in HillslopeHydrologyMod.F90::InitHillslope and
# surfrdMod.F90::surfrd_hillslope. All variables use dim1name=grlnd (=> the
# gridcell dimension is the *fastest*/last dim in CDL order; extra dim first),
# matching e.g. PCT_NAT_PFT(natpft, gridcell) in the base surfdata.
#
# Usage:
#   julia --project=<CLM.jl> make_hillslope_surfdata.jl <base_surfdata.nc> <out_hillslope.nc>
# =============================================================================

using NCDatasets
using Printf

const BASE = length(ARGS) >= 1 ? ARGS[1] :
    "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/surfdata_esmf/ctsm5.4.0/surfdata_ne3np4_hist_2000_78pfts_c251022.nc"
const OUT  = length(ARGS) >= 2 ? ARGS[2] :
    "/private/tmp/claude-501/hillslope_ne3np4_synthetic_c.nc"

# ---- synthetic catena definition (1 hillslope, 4 columns, lowland ci=1 -> upland ci=4)
const NHILL = 1     # number of representative hillslopes per landunit
const NMAXCOL = 4   # max columns per landunit (this catena uses all 4)

# per-column arrays are ordered ci = 1..NMAXCOL (ci=1 == lowland, adjacent to stream)
const COLUMN_INDEX    = Int32[1, 2, 3, 4]
const HILLSLOPE_INDEX = Int32[1, 1, 1, 1]                 # all belong to hillslope #1
const DOWNHILL_INDEX  = Int32[-999, 1, 2, 3]              # ci1 discharges to stream (-999)
const SEG_LEN         = 25.0                              # along-hill length of each column [m]
const HILL_DISTANCE   = Float64[0.0, 25.0, 50.0, 75.0]   # lower-edge distance from hill bottom [m]
const HILL_WIDTH      = Float64[100.0, 90.0, 80.0, 70.0] # lower-edge width [m]
const HILL_AREA       = HILL_WIDTH .* SEG_LEN             # column area [m2] = width*seg_len
const HILL_ELEV       = Float64[1.0, 6.0, 11.0, 16.0]    # mean elev rel. to gridcell mean [m], monotonic up
const HILL_SLOPE      = Float64[0.05, 0.08, 0.10, 0.12]  # along-hill slope [m/m]
const HILL_ASPECT     = fill(Float64(pi), NMAXCOL)       # azimuth [radians], south-facing (0..2pi)
const HILL_PFTNDX     = Int32[13, 13, 13, 13]            # natpft index (c3 grass); used only by FromFile/LowlandUpland

# stream channel props (only consumed when use_hillslope_routing=.true.)
const STREAM_DEPTH = 2.0    # [m]
const STREAM_WIDTH = 5.0    # [m]
const STREAM_SLOPE = 1.0e-3 # [m/m]

# ---- read base grid (support both 1D unstructured `gridcell` and 2D `lsmlat,lsmlon`).
# NCDatasets presents dims REVERSED from CDL; we build arrays in Julia-native
# (grid-fastest-first) order so the written CDL matches CTSM's convention
# (extra dim first, then grid dims): e.g. PCT_NAT_PFT(natpft, [lsmlat,] gridcell/lsmlon).
base = NCDataset(BASE, "r")
if haskey(base.dim, "gridcell")
    jgriddims  = ("gridcell",)                       # unstructured (ne3np4-style)
    ng = base.dim["gridcell"]
    jgridshape = (ng,)
    base_kind = "unstructured (gridcell)"
else
    @assert haskey(base.dim, "lsmlon") && haskey(base.dim, "lsmlat") "base must have gridcell or lsmlon/lsmlat"
    nlon = base.dim["lsmlon"]; nlat = base.dim["lsmlat"]
    jgriddims  = ("lsmlon", "lsmlat")                # 2D structured (single-point-style)
    ng = nlat * nlon
    jgridshape = (nlon, nlat)
    base_kind = "structured 2D (lsmlat,lsmlon)"
end
lon = Array{Float64}(base["LONGXY"][:])[:]         # flatten (Julia col-major)
lat = Array{Float64}(base["LATIXY"][:])[:]
pctnatveg = Array{Float64}(base["PCT_NATVEG"][:])[:]
close(base)
println("base grid kind: ", base_kind)

# choose the most-vegetated gridcell to host the catena (must be istsoil w/ wtgcell>0)
active = argmax(pctnatveg)
@printf("base grid: %d gridcells; active gridcell = %d (0-based %d) PCT_NATVEG=%.2f lon=%.3f lat=%.3f\n",
        ng, active, active-1, pctnatveg[active], lon[active], lat[active])

# ---- allocate full-grid arrays, grid-first (ng, extra) => Julia-native for writing
# (inactive gridcells => nhillcolumns=0, all-zero geometry)
nhillcolumns = zeros(Int32, ng)
pct_hillslope = zeros(Float64, ng, NHILL)            # (gridcell, nhillslope)
hillslope_index       = zeros(Int32, ng, NMAXCOL)    # (gridcell, nmaxhillcol)
column_index          = zeros(Int32, ng, NMAXCOL)
downhill_column_index = fill(Int32(-999), ng, NMAXCOL)
hillslope_slope       = zeros(Float64, ng, NMAXCOL)
hillslope_aspect      = zeros(Float64, ng, NMAXCOL)
hillslope_area        = zeros(Float64, ng, NMAXCOL)
hillslope_distance    = zeros(Float64, ng, NMAXCOL)
hillslope_width       = zeros(Float64, ng, NMAXCOL)
hillslope_elevation   = zeros(Float64, ng, NMAXCOL)
hillslope_pftndx      = zeros(Int32, ng, NMAXCOL)
hillslope_stream_depth = zeros(Float64, ng)
hillslope_stream_width = zeros(Float64, ng)
hillslope_stream_slope = zeros(Float64, ng)

# ---- write the synthetic catena into the active gridcell
nhillcolumns[active] = NMAXCOL
pct_hillslope[active, 1] = 100.0        # single hillslope occupies 100% of the hillslope fraction
hillslope_index[active, :]       .= HILLSLOPE_INDEX
column_index[active, :]          .= COLUMN_INDEX
downhill_column_index[active, :] .= DOWNHILL_INDEX
hillslope_slope[active, :]       .= HILL_SLOPE
hillslope_aspect[active, :]      .= HILL_ASPECT
hillslope_area[active, :]        .= HILL_AREA
hillslope_distance[active, :]    .= HILL_DISTANCE
hillslope_width[active, :]       .= HILL_WIDTH
hillslope_elevation[active, :]   .= HILL_ELEV
hillslope_pftndx[active, :]      .= HILL_PFTNDX
hillslope_stream_depth[active] = STREAM_DEPTH
hillslope_stream_width[active] = STREAM_WIDTH
hillslope_stream_slope[active] = STREAM_SLOPE

# ---- consistency checks (fail loudly)
@assert sum(nhillcolumns) == NMAXCOL "exactly one active gridcell expected"
@assert isapprox(sum(pct_hillslope[active, :]), 100.0) "pct_hillslope must sum to 100 over nhillslope"
@assert all(diff(HILL_ELEV) .> 0) "elevation must increase upland"
@assert all(diff(HILL_DISTANCE) .>= 0) "distance must be nondecreasing upland"
@assert DOWNHILL_INDEX[1] == -999 "lowland column must discharge to stream"
@assert all(HILL_AREA .> 0) "areas must be positive"
# weights that InitHillslope will recompute must sum to 1 over the columns
wts = (HILL_AREA ./ sum(HILL_AREA)) .* (100.0 * 0.01)
@assert isapprox(sum(wts), 1.0) "recomputed column weights must sum to 1"
@printf("catena OK: total hill area=%.1f m2, recomputed wtlunit sum=%.6f\n", sum(HILL_AREA), sum(wts))

# ---- write the netcdf hillslope_file
isfile(OUT) && rm(OUT)
ds = NCDataset(OUT, "c")
for (dn, dl) in zip(jgriddims, jgridshape); defDim(ds, dn, dl); end
defDim(ds, "nhillslope", NHILL)
defDim(ds, "nmaxhillcol", NMAXCOL)

ds.attrib["title"] = "Synthetic CTSM hillslope catena on a real base surface-data grid"
ds.attrib["note"]  = "REAL base grid/domain (LONGXY/LATIXY copied from base surfdata); " *
                     "SYNTHETIC idealized hillslope geometry (1 hillslope, 4 cols, upland->lowland). " *
                     "For CLM.jl<->Fortran CTSM hillslope-hydrology parity testing only."
ds.attrib["base_surfdata"] = BASE
ds.attrib["base_grid_kind"] = base_kind
ds.attrib["active_gridcell_1based"] = active
ds.attrib["created_by"] = "scripts/make_hillslope_surfdata.jl"

# reshape a flat (ng[, extra]) array to the base grid shape (jgridshape[, extra])
gridonly(v) = reshape(v, jgridshape)
withextra(v) = reshape(v, (jgridshape..., size(v, 2)))

function put(name, dims, data, atts=Dict{String,Any}())
    v = defVar(ds, name, eltype(data) <: Integer ? Int32 : Float64, dims)
    for (k, val) in atts; v.attrib[k] = val; end
    v[:] = data
end

gd  = jgriddims                    # grid dims (Julia order)
gde(extra) = (gd..., extra)        # grid dims + one extra (Julia order)

# domain coordinates (needed by check_domain_attributes)
put("LONGXY", gd, gridonly(lon), Dict("long_name"=>"longitude","units"=>"degrees east"))
put("LATIXY", gd, gridonly(lat), Dict("long_name"=>"latitude","units"=>"degrees north"))

# hillslope descriptors
put("nhillcolumns", gd, gridonly(nhillcolumns),
    Dict("long_name"=>"number of hillslope columns in each gridcell","units"=>"unitless"))
put("pct_hillslope", gde("nhillslope"), withextra(pct_hillslope),
    Dict("long_name"=>"percent of landunit occupied by each hillslope","units"=>"per cent"))
put("hillslope_index", gde("nmaxhillcol"), withextra(hillslope_index),
    Dict("long_name"=>"hillslope index of each column","units"=>"unitless"))
put("column_index", gde("nmaxhillcol"), withextra(column_index),
    Dict("long_name"=>"column index (1..nhillcolumns) of each column","units"=>"unitless"))
put("downhill_column_index", gde("nmaxhillcol"), withextra(downhill_column_index),
    Dict("long_name"=>"downhill (downstream) column index; -999 => discharges to stream","units"=>"unitless"))
put("hillslope_slope", gde("nmaxhillcol"), withextra(hillslope_slope),
    Dict("long_name"=>"along-hill slope of each column","units"=>"m/m"))
put("hillslope_aspect", gde("nmaxhillcol"), withextra(hillslope_aspect),
    Dict("long_name"=>"azimuth of each column","units"=>"radians"))
put("hillslope_area", gde("nmaxhillcol"), withextra(hillslope_area),
    Dict("long_name"=>"area of each column","units"=>"m2"))
put("hillslope_distance", gde("nmaxhillcol"), withextra(hillslope_distance),
    Dict("long_name"=>"distance of lower edge of column from hillslope bottom","units"=>"m"))
put("hillslope_width", gde("nmaxhillcol"), withextra(hillslope_width),
    Dict("long_name"=>"width of lower edge of each column","units"=>"m"))
put("hillslope_elevation", gde("nmaxhillcol"), withextra(hillslope_elevation),
    Dict("long_name"=>"mean elevation of column relative to gridcell mean","units"=>"m"))
put("hillslope_pftndx", gde("nmaxhillcol"), withextra(hillslope_pftndx),
    Dict("long_name"=>"pft index of each column (used by FromFile / PftLowlandUpland)","units"=>"unitless"))

# stream channel properties (used only when use_hillslope_routing=.true.)
put("hillslope_stream_depth", gd, gridonly(hillslope_stream_depth),
    Dict("long_name"=>"stream channel depth","units"=>"m"))
put("hillslope_stream_width", gd, gridonly(hillslope_stream_width),
    Dict("long_name"=>"stream channel width","units"=>"m"))
put("hillslope_stream_slope", gd, gridonly(hillslope_stream_slope),
    Dict("long_name"=>"stream channel slope","units"=>"m/m"))

close(ds)
@printf("wrote hillslope_file: %s\n", OUT)
