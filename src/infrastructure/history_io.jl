# ==========================================================================
# Minimal CLM-style history (h0) writer — a working subset of
# main/histFileMod.F90 (6100 lines; we do NOT port the full multi-tape
# framework). This provides a single time-averaged history tape that:
#   - registers fields with a getter into CLMInstances and an avgflag
#   - accumulates each field over an output interval, honoring the avgflag
#   - writes a CLM-style h0 NetCDF (time + spatial dim) with each field's
#     time-aggregated value, then resets the accumulators
#
# Averaging flags follow Fortran histFileMod (hist_update_hbuf_field):
#   "A" — time average (sum / nacs)
#   "I" — instantaneous (last sample)
#   "X" — maximum over the interval
#   "M" — minimum over the interval
#
# Subgrid handling (histFileMod hist_dov2xy / hist_type1d_pertape):
#   Each field carries a native 1d subgrid level (PATCH/COLUMN/LANDUNIT/
#   GRIDCELL). Fields of DIFFERENT levels coexist in one tape because each
#   level gets its OWN spatial NetCDF dim sized to that level's element count
#   (patch np, column nc, …). A field may additionally have a vertical (2d)
#   level — e.g. TSOI on levsoi — which gets its own level dim.
#
#   With `dov2xy=true` (Fortran hist_dov2xy) every field is instead area-
#   weight-aggregated up to the gridcell level (p2g / c2g / l2g via the wt*
#   weights in inst.patch/column/landunit), so the tape writes a single
#   shared gridcell dim — matching the default CLM gridded h0 output.
#
# This is a lighter, more generic complement to history_writer.jl (which
# does gridcell area-weighted daily aggregation). Here fields are written on
# their native subgrid level (or remapped to gridcell), NCDatasets style.
# ==========================================================================

"""
    HistSubgrid

Native 1d subgrid level of a history field, mirroring Fortran `type1d`
(`namep`/`namec`/`namel`/`nameg`). Determines the field's spatial dimension
and the weights used when remapping up to gridcell (`hist_dov2xy`).
"""
@enum HistSubgrid HIST_PATCH HIST_COLUMN HIST_LANDUNIT HIST_GRIDCELL

# String tag (used for dim names and NetCDF attributes) for each level.
_hist_subgrid_name(s::HistSubgrid) =
    s == HIST_PATCH     ? "patch"     :
    s == HIST_COLUMN    ? "column"    :
    s == HIST_LANDUNIT  ? "landunit"  :
                          "gridcell"

# Parse a user-facing string/symbol into a HistSubgrid.
function _hist_subgrid(level)::HistSubgrid
    level isa HistSubgrid && return level
    s = lowercase(String(level))
    s == "patch"     ? HIST_PATCH     :
    s == "column"    ? HIST_COLUMN    :
    s == "landunit"  ? HIST_LANDUNIT  :
    s == "gridcell"  ? HIST_GRIDCELL  :
    error("hist_addfld!: unknown subgrid level $(level); valid: patch/column/landunit/gridcell")
end

"""
    HistField

One registered history field: a `name`, `units`, a `getter` mapping
`CLMInstances` to the current per-element values, an averaging flag
(`"A"`, `"I"`, `"X"`, or `"M"`), and a native subgrid `level`.

The getter returns either an `AbstractVector` (1d field, one value per
subgrid element) or an `AbstractMatrix` of shape `(nelem, nlev)` (2d field
with a vertical dimension named `levdim`, length `nlev`).

`accum` is the running accumulator and `nacs` the per-slot sample count, both
stored flattened (length `nelem*nlev`, column-major over levels). `nelem` is
the spatial element count, fixed on the first accumulate (0 until then).
"""
mutable struct HistField
    name::String
    units::String
    long_name::String
    getter::Function
    avgflag::String
    level::HistSubgrid
    nlev::Int           # vertical levels (1 = pure 1d field)
    levdim::String      # name of the vertical dim (e.g. "levsoi"); "" if 1d
    accum::Vector{Float64}
    nacs::Vector{Int}
    nelem::Int
end

"""
    HistoryTape

A single time-averaged history tape: an ordered list of registered fields and
the total number of samples taken since the last write. Each field is written
on its own native subgrid level (a separate spatial NetCDF dim per level),
unless `dov2xy=true`, in which case all fields are area-weight aggregated to
the gridcell level and share one `gridcell` dim.
"""
Base.@kwdef mutable struct HistoryTape
    fields::Vector{HistField}   = HistField[]
    nsamples::Int               = 0
    dov2xy::Bool                = false   # true => area-weight aggregate all fields to gridcell

    # --- 2D (lon×lat) gridded output (histFileMod hist_dov2xy true-2D path) --
    # When `xy2d=true` AND `dov2xy=true` AND a g→(i,j) map is present, each
    # field's gridcell-aggregated value is written on `(lon, lat[, lev], time)`
    # dims instead of the flat `gridcell` dim, scattering gridcell `g` into grid
    # cell `(g2i[g], g2j[g])` and filling unmapped cells with the fill value.
    # `g2i`/`g2j` are 1-based lon/lat indices per gridcell (length ng); `nlon`/
    # `nlat` are the grid extents; `lon1d`/`lat1d` (optional) are the coordinate
    # values written for the lon/lat axes. With an empty map (or `xy2d=false`)
    # the tape falls back to the flat 1D `gridcell` dim (single-gridcell-safe).
    xy2d::Bool                  = false
    g2i::Vector{Int}            = Int[]
    g2j::Vector{Int}            = Int[]
    nlon::Int                   = 0
    nlat::Int                   = 0
    lon1d::Vector{Float64}      = Float64[]
    lat1d::Vector{Float64}      = Float64[]

    # --- CTSM time-slice / file metadata (histFileMod history_tape) ----------
    # `mfilt` = max time samples per history file (namelist hist_mfilt); when
    # `ntimes` reaches `mfilt` the file is "full" → rolls over to a new file.
    # `ndens` = output precision (namelist hist_ndens): 1 => double, 2 => single.
    # `ntimes` = time records written into the CURRENT file; `begtime` = start
    # of the current averaging interval (written as time_bounds(1)).
    mfilt::Int                  = 1
    ndens::Int                  = 2       # 2 = single (Float32), 1 = double (Float64)
    ntimes::Int                 = 0
    begtime::Float64            = 0.0
    # Set-level rollover bookkeeping: which output file (0-based) this tape is
    # currently writing into, and that file's path ("" = not yet opened). When
    # `ntimes` reaches `mfilt` the file is full and the next write rolls to a
    # new file (`file_index` + 1, `ntimes` reset to 0).
    file_index::Int             = 0
    cur_filename::String        = ""
end

const _HIST_VALID_AVGFLAGS = ("A", "I", "X", "M")

"""
    hist_addfld!(tape, name, units, getter; avgflag="A", long_name=name,
                 level="patch", nlev=1, levdim="levsoi")

Register a field on `tape`. `getter(inst::CLMInstances)` returns the current
values: an `AbstractVector` for a 1d field, or an `AbstractMatrix`
`(nelem, nlev)` for a 2d field. `avgflag` selects the time aggregation
(`"A"` average, `"I"` instantaneous, `"X"` max, `"M"` min). `level` is the
native subgrid level (`"patch"`/`"column"`/`"landunit"`/`"gridcell"`). For a
2d field pass `nlev>1` and a `levdim` name for the vertical dimension.
"""
function hist_addfld!(tape::HistoryTape, name::AbstractString, units::AbstractString,
                      getter::Function; avgflag::AbstractString="A",
                      long_name::AbstractString=name,
                      level="patch", nlev::Integer=1,
                      levdim::AbstractString="levsoi")
    avgflag in _HIST_VALID_AVGFLAGS ||
        error("hist_addfld!: unknown averaging flag $(avgflag); valid: $(_HIST_VALID_AVGFLAGS)")
    any(f -> f.name == name, tape.fields) &&
        error("hist_addfld!: field $(name) already registered")
    sg = _hist_subgrid(level)
    push!(tape.fields, HistField(String(name), String(units), String(long_name),
                                 getter, String(avgflag), sg,
                                 Int(nlev), nlev > 1 ? String(levdim) : "",
                                 Float64[], Int[], 0))
    return tape
end

# Initialize a field's accumulator (flattened nelem*nlev) to the natural
# identity for its avgflag, recording the spatial element count.
function _hist_reset_field!(fld::HistField, nelem::Int)
    fld.nelem = nelem
    n = nelem * fld.nlev
    fld.accum = fill(_hist_init_value(fld.avgflag), n)
    fld.nacs  = zeros(Int, n)
    return nothing
end

_hist_init_value(avgflag::String) =
    avgflag == "X" ? -Inf :
    avgflag == "M" ?  Inf :
    0.0  # "A" sum, "I" last

# Spatial element count of a field's raw getter output (rows for a matrix).
_hist_nelem(data::AbstractVector) = length(data)
_hist_nelem(data::AbstractMatrix) = size(data, 1)

"""
    hist_accumulate!(tape, inst)

Sample every registered field from `inst` and fold it into the accumulators
according to each field's avgflag. Values equal to `SPVAL` or non-finite are
skipped for that slot (no sample counted), matching Fortran behavior.
"""
function hist_accumulate!(tape::HistoryTape, inst)
    for fld in tape.fields
        data = fld.getter(inst)
        nelem = _hist_nelem(data)
        nelem == 0 && continue

        # A 2d field self-adopts its level count from the getter's matrix
        # (robust to varpar.nlevsoi changing between registration and use).
        if data isa AbstractMatrix
            nlev = size(data, 2)
            fld.nlev > 1 || error("hist_accumulate!: field $(fld.name) getter returned a matrix but the field is registered 1d (nlev=1)")
            fld.nlev = nlev
        end

        # Lazily (re)fix the field's spatial length / accumulator layout.
        if fld.nelem != nelem || length(fld.accum) != nelem * fld.nlev
            _hist_reset_field!(fld, nelem)
        end

        if data isa AbstractMatrix
            @inbounds for ell in 1:fld.nlev, k in 1:nelem
                _hist_fold!(fld, (ell - 1) * nelem + k, Float64(data[k, ell]))
            end
        else
            @inbounds for k in 1:nelem
                _hist_fold!(fld, k, Float64(data[k]))
            end
        end
    end
    tape.nsamples += 1
    return nothing
end

# Fold one value into a single accumulator slot per the field's avgflag.
@inline function _hist_fold!(fld::HistField, idx::Int, v::Float64)
    (isfinite(v) && v != SPVAL) || return nothing
    if fld.avgflag == "I"
        fld.accum[idx] = v
        fld.nacs[idx] = 1
    elseif fld.avgflag == "A"
        fld.accum[idx] += v
        fld.nacs[idx] += 1
    elseif fld.avgflag == "X"
        fld.accum[idx] = max(fld.accum[idx], v)
        fld.nacs[idx] = 1
    elseif fld.avgflag == "M"
        fld.accum[idx] = min(fld.accum[idx], v)
        fld.nacs[idx] = 1
    end
    return nothing
end

"""
    hist_field_value(fld) -> Vector{Float64} or Matrix{Float64}

Finalize a field's accumulated value: for `"A"` divide by the per-slot sample
count; otherwise return the running value. Slots that never received a valid
sample are set to `SPVAL`. A 1d field returns a `Vector` of length `nelem`; a
2d field returns a `Matrix` of shape `(nelem, nlev)`.
"""
function hist_field_value(fld::HistField)
    n = length(fld.accum)
    out = fill(SPVAL, n)
    @inbounds for k in 1:n
        if fld.nacs[k] > 0
            out[k] = fld.avgflag == "A" ? fld.accum[k] / fld.nacs[k] : fld.accum[k]
        end
    end
    return fld.nlev > 1 ? reshape(out, fld.nelem, fld.nlev) : out
end

# ----------------------------------------------------------------------
# Gridcell remapping (Fortran hist_dov2xy / hist_dov2xy=.true.):
# area-weight-aggregate a field from its native subgrid level up to the
# gridcell level using the wt* weights in inst.patch/column/landunit.
# Works for both 1d and 2d (per-level) fields.
# ----------------------------------------------------------------------

# Map level -> (index-into-gridcell vector, weight-to-gridcell vector).
function _hist_subgrid_to_gridcell_maps(fld::HistField, inst)
    if fld.level == HIST_PATCH
        return inst.patch.gridcell, inst.patch.wtgcell
    elseif fld.level == HIST_COLUMN
        return inst.column.gridcell, inst.column.wtgcell
    elseif fld.level == HIST_LANDUNIT
        return inst.landunit.gridcell, inst.landunit.wtgcell
    else  # HIST_GRIDCELL — identity
        return Int[], Float64[]
    end
end

"""
    hist_dov2xy(fld, vals, ng, inst) -> Vector or Matrix on the gridcell dim

Area-weight-aggregate a finalized field value `vals` (Vector or Matrix from
`hist_field_value`) from its native subgrid level up to `ng` gridcells, using
the field's wt-to-gridcell weights. Gridcells with no contributing weight are
set to `SPVAL`. Gridcell-level fields pass through (truncated/padded to `ng`).
"""
function hist_dov2xy(fld::HistField, vals, ng::Int, inst)
    idx, wt = _hist_subgrid_to_gridcell_maps(fld, inst)
    nlev = fld.nlev

    if fld.level == HIST_GRIDCELL
        out = fill(SPVAL, ng, nlev)
        src = nlev > 1 ? vals : reshape(vals, :, 1)
        @inbounds for ell in 1:nlev, g in 1:min(ng, size(src, 1))
            out[g, ell] = src[g, ell]
        end
        return nlev > 1 ? out : vec(out)
    end

    nelem = fld.nelem
    gval = zeros(Float64, ng, nlev)
    gwt  = zeros(Float64, ng, nlev)
    src  = nlev > 1 ? vals : reshape(vals, :, 1)

    @inbounds for ell in 1:nlev, e in 1:nelem
        (e <= length(idx) && e <= length(wt)) || continue
        g = idx[e]
        (g >= 1 && g <= ng) || continue
        v = src[e, ell]
        w = wt[e]
        (isfinite(v) && v != SPVAL && isfinite(w)) || continue
        gval[g, ell] += v * w
        gwt[g, ell]  += w
    end

    out = fill(SPVAL, ng, nlev)
    @inbounds for ell in 1:nlev, g in 1:ng
        if gwt[g, ell] > 0.0
            out[g, ell] = gval[g, ell] / gwt[g, ell]
        end
    end
    return nlev > 1 ? out : vec(out)
end

# ----------------------------------------------------------------------
# 2D (lon×lat) gridded output (Fortran hist_dov2xy true-2D path).
# A multi-gridcell run lays its `ng` gridcells onto a regular `(nlon, nlat)`
# grid via a g→(i,j) map. `hist_grid_map_from_lonlat` derives that map (and the
# 1D lon/lat axis coordinates) from per-gridcell lon/lat; `hist_set_grid_map!`
# stamps a tape with an explicit or derived map; `_hist_scatter_2d` scatters a
# gridcell-aggregated value onto the `(nlon, nlat[, nlev])` grid, filling
# unmapped cells with SPVAL.
# ----------------------------------------------------------------------

"""
    hist_grid_map_from_lonlat(lon, lat; atol=1e-6)
        -> (g2i, g2j, nlon, nlat, lon1d, lat1d)

Derive a gridcell → `(i, j)` lon/lat-index map from per-gridcell `lon`/`lat`
coordinate vectors (length `ng`). The unique sorted longitude values define the
`nlon` lon axis and the unique sorted latitude values the `nlat` lat axis; each
gridcell `g` maps to `(g2i[g], g2j[g])` = the index of its lon/lat on those
axes. `lon1d`/`lat1d` are the sorted unique coordinate values (the axis coords).

Values within `atol` of each other are treated as the same axis coordinate, so
floating point grid coordinates collapse cleanly to integer indices. Use this
when a multi-gridcell run carries only gridcell lon/lat (no explicit i/j map).
"""
function hist_grid_map_from_lonlat(lon::AbstractVector, lat::AbstractVector; atol::Real=1e-6)
    ng = length(lon)
    ng == length(lat) ||
        error("hist_grid_map_from_lonlat: lon/lat length mismatch ($(ng) vs $(length(lat)))")
    lon1d = _hist_unique_sorted(lon, atol)
    lat1d = _hist_unique_sorted(lat, atol)
    g2i = Vector{Int}(undef, ng)
    g2j = Vector{Int}(undef, ng)
    @inbounds for g in 1:ng
        g2i[g] = _hist_axis_index(lon1d, lon[g], atol)
        g2j[g] = _hist_axis_index(lat1d, lat[g], atol)
    end
    return (g2i = g2i, g2j = g2j, nlon = length(lon1d), nlat = length(lat1d),
            lon1d = lon1d, lat1d = lat1d)
end

# Sorted unique axis coordinates, collapsing values within `atol`.
function _hist_unique_sorted(v::AbstractVector, atol::Real)
    s = sort(Float64.(collect(v)))
    out = Float64[]
    for x in s
        if isempty(out) || abs(x - out[end]) > atol
            push!(out, x)
        end
    end
    return out
end

# 1-based index of `x` on a sorted axis `ax` (nearest within `atol`).
function _hist_axis_index(ax::Vector{Float64}, x::Real, atol::Real)
    best = 1
    bestd = Inf
    @inbounds for i in eachindex(ax)
        d = abs(ax[i] - Float64(x))
        if d < bestd
            bestd = d
            best = i
        end
    end
    return best
end

"""
    hist_set_grid_map!(tape; g2i=nothing, g2j=nothing, nlon=0, nlat=0,
                       lon=nothing, lat=nothing, lon1d=Float64[], lat1d=Float64[],
                       atol=1e-6) -> tape

Stamp `tape` with a gridcell → `(i, j)` 2D-output map and enable `xy2d`. Either
pass an explicit `g2i`/`g2j` (1-based) with `nlon`/`nlat` (and optional axis
coords `lon1d`/`lat1d`), or pass per-gridcell `lon`/`lat` and the map is derived
via `hist_grid_map_from_lonlat`. With `dov2xy=true` the tape then writes fields
on `(lon, lat[, lev], time)` dims.
"""
function hist_set_grid_map!(tape::HistoryTape;
                            g2i::Union{AbstractVector,Nothing}=nothing,
                            g2j::Union{AbstractVector,Nothing}=nothing,
                            nlon::Integer=0, nlat::Integer=0,
                            lon::Union{AbstractVector,Nothing}=nothing,
                            lat::Union{AbstractVector,Nothing}=nothing,
                            lon1d::AbstractVector=Float64[],
                            lat1d::AbstractVector=Float64[],
                            atol::Real=1e-6)
    if g2i === nothing || g2j === nothing
        (lon !== nothing && lat !== nothing) ||
            error("hist_set_grid_map!: pass explicit g2i/g2j or per-gridcell lon/lat")
        m = hist_grid_map_from_lonlat(lon, lat; atol=atol)
        tape.g2i = m.g2i; tape.g2j = m.g2j
        tape.nlon = m.nlon; tape.nlat = m.nlat
        tape.lon1d = m.lon1d; tape.lat1d = m.lat1d
    else
        length(g2i) == length(g2j) ||
            error("hist_set_grid_map!: g2i/g2j length mismatch")
        tape.g2i = Int.(collect(g2i)); tape.g2j = Int.(collect(g2j))
        tape.nlon = nlon > 0 ? Int(nlon) : (isempty(tape.g2i) ? 0 : maximum(tape.g2i))
        tape.nlat = nlat > 0 ? Int(nlat) : (isempty(tape.g2j) ? 0 : maximum(tape.g2j))
        tape.lon1d = Float64.(collect(lon1d))
        tape.lat1d = Float64.(collect(lat1d))
    end
    tape.xy2d = true
    return tape
end

# Whether the tape's 2D map is usable (xy2d on + a non-degenerate grid covering
# all `ng` gridcells). Single-gridcell / no-map tapes fall back to the 1D path.
function _hist_xy2d_active(tape::HistoryTape, ng::Int)
    tape.xy2d || return false
    (tape.nlon >= 1 && tape.nlat >= 1) || return false
    (length(tape.g2i) >= ng && length(tape.g2j) >= ng) || return false
    return true
end

"""
    _hist_scatter_2d(gvals, tape, ng) -> Array (nlon, nlat[, nlev])

Scatter a gridcell-aggregated value `gvals` (Vector length `ng`, or Matrix
`(ng, nlev)` from `hist_dov2xy`) onto the tape's `(nlon, nlat[, nlev])` grid via
`g→(i,j)`. Cells with no gridcell are filled with SPVAL.
"""
function _hist_scatter_2d(gvals, tape::HistoryTape, ng::Int)
    src = gvals isa AbstractMatrix ? gvals : reshape(gvals, :, 1)
    nlev = size(src, 2)
    out = fill(SPVAL, tape.nlon, tape.nlat, nlev)
    @inbounds for g in 1:ng
        i = tape.g2i[g]; j = tape.g2j[g]
        (i >= 1 && i <= tape.nlon && j >= 1 && j <= tape.nlat) || continue
        for ell in 1:nlev
            out[i, j, ell] = src[g, ell]
        end
    end
    return nlev > 1 ? out : reshape(out, tape.nlon, tape.nlat)
end

# Output NetCDF element type for a tape's `ndens` (Fortran hist_ndens):
# ndens==1 => ncd_double (Float64), else (ndens==2) => ncd_float (Float32).
_hist_ncprec(ndens::Integer) = ndens == 1 ? Float64 : Float32

# Compute the CTSM calendar date metadata (mcdate/mcsec/mdcur/mscur) for a
# model `time` expressed in the units of `time_units` ("<unit> since <date>").
# `mcdate` = yr*10000+mon*100+day, `mcsec` = seconds of the current day,
# `mdcur` = whole days since the base date, `mscur` = seconds within that day.
# Mirrors histFileMod get_curr_date / get_curr_time. Only the common "days
# since" / "seconds since" / "hours since" forms are decoded; anything else
# falls back to (time, 0) days/seconds without a calendar date.
function _hist_time_metadata(time::Real, time_units::AbstractString)
    parts = split(time_units)
    unit  = isempty(parts) ? "days" : lowercase(parts[1])
    factor = unit == "seconds" ? 1.0 / 86400.0 :
             unit == "hours"   ? 1.0 / 24.0    :
             unit == "minutes" ? 1.0 / 1440.0  :
             1.0   # "days" (default)
    days_since = Float64(time) * factor

    # Locate the base date string after "since".
    base = DateTime(2000, 1, 1)
    si = findfirst(==("since"), lowercase.(parts))
    if si !== nothing && si < length(parts)
        datestr = join(parts[(si+1):end], " ")
        try
            base = _hist_parse_basedate(datestr)
        catch
            base = DateTime(2000, 1, 1)
        end
    end

    whole_days  = floor(Int, days_since)
    frac_sec    = round(Int, (days_since - whole_days) * 86400.0)
    cur = base + Dates.Day(whole_days) + Dates.Second(frac_sec)
    yr  = Dates.year(cur)
    mon = Dates.month(cur)
    dy  = Dates.day(cur)
    mcsec  = (Dates.hour(cur) * 3600) + (Dates.minute(cur) * 60) + Dates.second(cur)
    mcdate = yr * 10000 + mon * 100 + dy
    mdcur  = whole_days
    mscur  = frac_sec
    return (mcdate = mcdate, mcsec = mcsec, mdcur = mdcur, mscur = mscur)
end

# Parse the "<date>[ <time>]" portion of a CF time-units string into DateTime.
function _hist_parse_basedate(s::AbstractString)
    s = strip(s)
    toks = split(s)
    d = Date(toks[1], dateformat"y-m-d")
    if length(toks) >= 2
        t = Time(toks[2], dateformat"H:M:S")
        return DateTime(d) + Dates.Hour(Dates.hour(t)) +
               Dates.Minute(Dates.minute(t)) + Dates.Second(Dates.second(t))
    end
    return DateTime(d)
end

# Number of gridcells inferred from the subgrid index vectors in `inst`.
function _hist_ngridcell(inst)
    ng = 0
    for v in (inst.column.gridcell, inst.patch.gridcell, inst.landunit.gridcell)
        if !isempty(v)
            m = maximum(v)
            m > ng && (ng = m)
        end
    end
    return ng
end

"""
    hist_write!(tape, filename; time=0.0, time_units="days since 2000-01-01",
                calendar="noleap", append=false, reset=true, inst=nothing,
                begtime=nothing)

Write the time-aggregated tape to a CLM-style h0 NetCDF `filename`. Defines a
`time` dimension (unlimited) and, for each subgrid level present, a spatial
dim sized to that level's element count (or, if `tape.dov2xy`, a single shared
`gridcell` dim). 2d fields additionally define their vertical level dim.
Appends one time record with each field's finalized value and (by default)
resets the accumulators for the next interval.

When `tape.dov2xy` is on AND the tape carries a 2D grid map (`tape.xy2d` with a
`g→(i,j)` map covering all gridcells, e.g. via `hist_set_grid_map!`), fields are
instead written on `(lon, lat[, lev], time)` dims: each gridcell `g`'s value is
scattered to grid cell `(g2i[g], g2j[g])` and cells with no gridcell are filled
with the `_FillValue` (SPVAL). With no map (or a single gridcell) it falls back
to the flat `gridcell` dim, so the single-gridcell path is unchanged.

Alongside `time` it writes the CTSM averaging-interval metadata (mirroring
`histFileMod`): `time_bounds(hist_interval,time)` (interval endpoints, with the
`time` variable carrying `bounds="time_bounds"`), and the calendar bookkeeping
`mcdate`/`mcsec`/`mdcur`/`mscur`. The interval start is `begtime` (defaults to
the tape's stored `tape.begtime`); the end is `time`. After the record the
tape's `begtime` advances to `time` so the next interval picks up where this
one ended.

Output precision follows the tape's `ndens` (Fortran hist_ndens): `ndens==1`
writes Float64, otherwise Float32. The `time`/`time_bounds` axes are always
Float64 (CTSM writes them at `ncprec`, but full precision avoids date-decode
loss).

`inst` is required when `tape.dov2xy=true` (supplies the subgrid weights for
gridcell aggregation). If `append=true` and the file exists, the record is
appended; otherwise the file is (re)created.
"""
function hist_write!(tape::HistoryTape, filename::AbstractString;
                     time::Real=0.0,
                     time_units::AbstractString="days since 2000-01-01",
                     calendar::AbstractString="noleap",
                     append::Bool=false,
                     reset::Bool=true,
                     inst=nothing,
                     begtime::Union{Real,Nothing}=nothing)
    isempty(tape.fields) && error("hist_write!: no fields registered")
    tape.nsamples == 0 && error("hist_write!: no samples accumulated")
    tape.dov2xy && inst === nothing &&
        error("hist_write!: dov2xy=true requires `inst` for gridcell aggregation")

    ng = tape.dov2xy ? _hist_ngridcell(inst) : 0
    prec = _hist_ncprec(tape.ndens)   # field output precision (Float32/Float64)
    btime = begtime === nothing ? tape.begtime : Float64(begtime)
    meta  = _hist_time_metadata(time, time_units)

    # True-2D gridded output only when dov2xy is on AND a usable (lon×lat) map
    # covers all ng gridcells; otherwise the 1D gridcell/native path (a single-
    # gridcell run with xy2d=false stays byte-identical to before).
    xy2d = tape.dov2xy && _hist_xy2d_active(tape, ng)

    # Finalize each field's value (and remap to gridcell if dov2xy), and
    # determine the spatial dim name + length it will be written on. In 2D
    # mode the value is scattered to (lon, lat[, lev]) and written on the
    # (lon, lat) dim pair instead of a single spatial dim.
    nf = length(tape.fields)
    finals  = Vector{Any}(undef, nf)
    dimname = Vector{String}(undef, nf)
    dimlen  = Vector{Int}(undef, nf)
    for (i, fld) in enumerate(tape.fields)
        v = hist_field_value(fld)
        if tape.dov2xy
            gv = hist_dov2xy(fld, v, ng, inst)
            if xy2d
                finals[i]  = _hist_scatter_2d(gv, tape, ng)
                dimname[i] = "lon"   # paired with "lat" in the 2D branch below
                dimlen[i]  = tape.nlon
            else
                finals[i]  = gv
                dimname[i] = "gridcell"
                dimlen[i]  = ng
            end
        else
            finals[i]  = v
            dimname[i] = _hist_subgrid_name(fld.level)
            dimlen[i]  = fld.nelem
        end
    end

    creating = !(append && isfile(filename))
    mode = creating ? "c" : "a"
    NCDataset(filename, mode) do ds
        if creating
            defDim(ds, "time", Inf)  # unlimited
            defDim(ds, "hist_interval", 2)  # CTSM time-bounds endpoint dim
            defVar(ds, "time", Float64, ("time",);
                   attrib = Dict("units" => time_units,
                                 "calendar" => calendar,
                                 "long_name" => "time",
                                 "bounds" => "time_bounds"))
            # Averaging-interval endpoints + calendar bookkeeping (histFileMod).
            defVar(ds, "time_bounds", Float64, ("hist_interval", "time");
                   attrib = Dict("long_name" => "history time interval endpoints",
                                 "units" => time_units,
                                 "calendar" => calendar))
            defVar(ds, "mcdate", Int32, ("time",);
                   attrib = Dict("long_name" => "current date (YYYYMMDD)"))
            defVar(ds, "mcsec", Int32, ("time",);
                   attrib = Dict("long_name" => "current seconds of current date",
                                 "units" => "s"))
            defVar(ds, "mdcur", Int32, ("time",);
                   attrib = Dict("long_name" => "current day (from base day)"))
            defVar(ds, "mscur", Int32, ("time",);
                   attrib = Dict("long_name" => "current seconds of current day"))
            seen = Set{String}()
            if xy2d
                # 2D grid: a shared (lon, lat) dim pair + optional axis coords.
                defDim(ds, "lon", tape.nlon)
                defDim(ds, "lat", tape.nlat)
                push!(seen, "lon"); push!(seen, "lat")
                if length(tape.lon1d) == tape.nlon
                    defVar(ds, "lon", Float64, ("lon",);
                           attrib = Dict("long_name" => "coordinate longitude",
                                         "units" => "degrees_east"))
                    ds["lon"][:] = tape.lon1d
                end
                if length(tape.lat1d) == tape.nlat
                    defVar(ds, "lat", Float64, ("lat",);
                           attrib = Dict("long_name" => "coordinate latitude",
                                         "units" => "degrees_north"))
                    ds["lat"][:] = tape.lat1d
                end
            else
                # One spatial dim per distinct level name.
                for i in 1:nf
                    if !(dimname[i] in seen)
                        defDim(ds, dimname[i], dimlen[i])
                        push!(seen, dimname[i])
                    end
                end
            end
            # Vertical level dims for 2d fields.
            for fld in tape.fields
                if fld.nlev > 1 && !(fld.levdim in seen)
                    defDim(ds, fld.levdim, fld.nlev)
                    push!(seen, fld.levdim)
                end
            end
            # Field variables.  1D mode: (space,time) / (space,lev,time).
            # 2D mode: (lon,lat,time) / (lon,lat,lev,time).
            for (i, fld) in enumerate(tape.fields)
                dims = if xy2d
                    fld.nlev > 1 ? ("lon", "lat", fld.levdim, "time") :
                                   ("lon", "lat", "time")
                else
                    fld.nlev > 1 ? (dimname[i], fld.levdim, "time") :
                                   (dimname[i], "time")
                end
                defVar(ds, fld.name, prec, dims;
                       attrib = Dict("long_name" => fld.long_name,
                                     "units" => fld.units,
                                     "cell_methods" => _hist_cell_method(fld.avgflag),
                                     "subgrid_level" =>
                                         tape.dov2xy ? "gridcell" : _hist_subgrid_name(fld.level),
                                     "_FillValue" => prec(SPVAL)))
            end
        end

        ti = ds.dim["time"] + 1
        ds["time"][ti] = Float64(time)
        ds["time_bounds"][:, ti] = Float64[btime, Float64(time)]
        ds["mcdate"][ti] = Int32(meta.mcdate)
        ds["mcsec"][ti]  = Int32(meta.mcsec)
        ds["mdcur"][ti]  = Int32(meta.mdcur)
        ds["mscur"][ti]  = Int32(meta.mscur)
        for (i, fld) in enumerate(tape.fields)
            if xy2d
                if fld.nlev > 1
                    ds[fld.name][:, :, :, ti] = finals[i]
                else
                    ds[fld.name][:, :, ti] = finals[i]
                end
            elseif fld.nlev > 1
                ds[fld.name][:, :, ti] = finals[i]
            else
                ds[fld.name][:, ti] = finals[i]
            end
        end
    end

    # Advance the averaging-interval start and bump the in-file record count.
    tape.begtime = Float64(time)
    tape.ntimes += 1

    if reset
        for fld in tape.fields
            if !isempty(fld.accum)
                _hist_reset_field!(fld, fld.nelem)
            end
        end
        tape.nsamples = 0
    end

    return filename
end

_hist_cell_method(avgflag::String) =
    avgflag == "A" ? "time: mean" :
    avgflag == "X" ? "time: maximum" :
    avgflag == "M" ? "time: minimum" :
    "time: point"  # "I"

# ----------------------------------------------------------------------
# Helpers for the built-in column 2d field (TSOI). Defined before the
# master field list so the `const` literal can reference them.
# ----------------------------------------------------------------------

# Number of soil layers for column 2d fields (TSOI/H2OSOI). Falls back to the
# standard CLM nlevsoi (10) if the global varpar block is unset.
function _hist_nlevsoi()
    n = varpar.nlevsoi
    return n > 0 ? n : 10
end

# Snow-layer offset into t_soisno_col: Fortran layer j maps to Julia j+nlevsno
# (snow levels -nlevsno+1:0 precede soil levels 1:nlevgrnd). Falls back to 5,
# the standard CLM nlevsno, if the global varpar block is unset.
function _hist_nlevsno()
    n = varpar.nlevsno
    return n > 0 ? n : 5
end

"""
    history_tsoi_col(inst) -> Matrix{Float64} (nc, nlevsoi)

Soil temperature on the soil layers (excludes snow layers), as a
`(ncolumn, nlevsoi)` matrix for a 2d column history field. `t_soisno_col` is
stored `(col, -nlevsno+1 : nlevmaxurbgrnd)`; soil layer j lives at Julia
column offset `nlevsno + j`.
"""
function history_tsoi_col(inst::CLMInstances)
    t = inst.temperature.t_soisno_col
    (t isa AbstractMatrix && size(t, 1) > 0) || return Matrix{Float64}(undef, 0, 0)
    nl = _hist_nlevsoi()
    nsno = _hist_nlevsno()
    nc = size(t, 1)
    out = Matrix{Float64}(undef, nc, nl)
    @inbounds for j in 1:nl, c in 1:nc
        col = nsno + j
        out[c, j] = col <= size(t, 2) ? Float64(t[c, col]) : SPVAL
    end
    return out
end

"""
    default_history_tape(; dov2xy=false) -> HistoryTape

A starter tape registering the key lnd2atm / energy / water diagnostics on
their native subgrid level (a mix of patch 1d, column 1d, and the column 2d
field `TSOI`). Useful as a quick-start h0 set; callers can `hist_addfld!` more
fields as needed. With `dov2xy=true` all fields are area-weight aggregated to
the gridcell level on write.
"""
function default_history_tape(; dov2xy::Bool=false)
    tape = HistoryTape(dov2xy=dov2xy)
    for spec in _MASTER_HIST_FIELDS
        hist_addfld!(tape, spec.name, spec.units, spec.getter;
                     avgflag=spec.avgflag, long_name=spec.long_name,
                     level=spec.level, nlev=spec.nlev, levdim=spec.levdim)
    end
    return tape
end

# ==========================================================================
# Multi-tape framework — mirrors histFileMod's master field list +
# htapes_build (fincl/fexcl per tape) + set_hist_filename (h0/h1/h2…).
#
# A `HistFieldSpec` is a registration in the *master* field list (analogous
# to Fortran's `allhistfldlist`): name, units, long_name, getter, the field's
# default avgflag, and `on_by_default` (whether it lands on a tape absent any
# fincl/fexcl, like Fortran's per-tape `actflag`).
# ==========================================================================

"""
    HistFieldSpec

One entry in the master history field list. `on_by_default=false` mirrors
Fortran's "inactive" master fields that only appear when named in `fincl`.
"""
struct HistFieldSpec
    name::String
    units::String
    long_name::String
    getter::Function
    avgflag::String
    on_by_default::Bool
    level::HistSubgrid
    nlev::Int
    levdim::String
end

"""
    hist_master_field(name, units, getter; avgflag="A", long_name=name,
                      on_by_default=true, level="patch", nlev=1,
                      levdim="levsoi") -> HistFieldSpec

Build a master-list field spec. Equivalent to Fortran `hist_addfld1d` /
`hist_addfld2d`'s registration into `allhistfldlist` (the per-tape selection
happens later in `htapes_build`). `level` is the native subgrid level; pass
`nlev>1` + `levdim` for a 2d (vertical-level) field.
"""
function hist_master_field(name::AbstractString, units::AbstractString, getter::Function;
                           avgflag::AbstractString="A", long_name::AbstractString=name,
                           on_by_default::Bool=true,
                           level="patch", nlev::Integer=1,
                           levdim::AbstractString="levsoi")
    avgflag in _HIST_VALID_AVGFLAGS ||
        error("hist_master_field: unknown averaging flag $(avgflag); valid: $(_HIST_VALID_AVGFLAGS)")
    return HistFieldSpec(String(name), String(units), String(long_name),
                         getter, String(avgflag), on_by_default,
                         _hist_subgrid(level), Int(nlev),
                         nlev > 1 ? String(levdim) : "")
end

# The default master field list (the set of always-on h0 diagnostics), a mix
# of patch 1d, column 1d, and a column 2d (vertical-level) field.
const _MASTER_HIST_FIELDS = HistFieldSpec[
    hist_master_field("EFLX_LH_TOT", "W/m^2",
                      inst -> inst.energyflux.eflx_lh_tot_patch;
                      level="patch", long_name="total latent heat flux"),
    hist_master_field("FSH", "W/m^2",
                      inst -> inst.energyflux.eflx_sh_tot_patch;
                      level="patch", long_name="sensible heat flux"),
    hist_master_field("FSA", "W/m^2",
                      inst -> inst.solarabs.fsa_patch;
                      level="patch", long_name="absorbed solar radiation"),
    hist_master_field("TSA", "K",
                      inst -> inst.temperature.t_ref2m_patch;
                      level="patch", long_name="2m air temperature"),
    hist_master_field("TG", "K",
                      inst -> inst.temperature.t_grnd_col;
                      level="column", long_name="ground temperature"),
    hist_master_field("QRUNOFF", "mm/s",
                      inst -> inst.water.waterfluxbulk_inst.wf.qflx_runoff_col;
                      level="column", long_name="total liquid runoff"),
    hist_master_field("QFLX_EVAP_TOT", "mm/s",
                      inst -> inst.water.waterfluxbulk_inst.wf.qflx_evap_tot_col;
                      level="column", long_name="total evaporation"),
    hist_master_field("H2OSNO", "mm",
                      history_h2osno_total_col;
                      level="column", long_name="snow water equivalent"),
    hist_master_field("TSOI", "K",
                      history_tsoi_col;
                      level="column", nlev=_hist_nlevsoi(), levdim="levsoi",
                      long_name="soil temperature"),
]

"""
    default_master_field_list() -> Vector{HistFieldSpec}

A copy of the built-in master history field list. Callers can append their own
`hist_master_field` specs before building a tape set.
"""
default_master_field_list() = copy(_MASTER_HIST_FIELDS)

# --- fincl/fexcl `name:flag` parsing (Fortran getname / getflag) ------------

"""
    hist_getname(entry) -> String

Name portion of a fincl/fexcl entry, dropping any `":flag"` suffix.
Mirrors Fortran `getname`.
"""
function hist_getname(entry::AbstractString)
    i = findfirst(==(':'), entry)
    return i === nothing ? String(strip(entry)) : String(strip(entry[1:prevind(entry, i)]))
end

"""
    hist_getflag(entry) -> String

Averaging-flag portion of a fincl/fexcl entry (after `":"`), or `""` if none.
Mirrors Fortran `getflag`.
"""
function hist_getflag(entry::AbstractString)
    i = findfirst(==(':'), entry)
    return i === nothing ? "" : String(strip(entry[nextind(entry, i):end]))
end

"""
    set_hist_filename(case, tape_index; ext_date="", compname="clm2",
                      file_index=0) -> String

Construct a CLM history filename. `tape_index` is 1-based in the code but the
filename uses `h(tape_index-1)`, so tape 1 → `…clm2.h0…`, tape 2 → `…h1…`.
Mirrors Fortran `set_hist_filename` (the date string is opaque here; pass it
via `ext_date`, e.g. `"2000-01"`). `file_index` (the mfilt rollover counter)
appends a `-NNNNN` suffix when `ext_date` is empty and `file_index>0`, so
successive files of one tape get distinct names. Always ends in `.nc`.
"""
function set_hist_filename(case::AbstractString, tape_index::Integer;
                           ext_date::AbstractString="", compname::AbstractString="clm2",
                           file_index::Integer=0)
    tape_index >= 1 || error("set_hist_filename: tape_index must be ≥ 1")
    hidx = tape_index - 1
    suffix = isempty(ext_date) ?
        (file_index > 0 ? "-" * lpad(file_index, 5, '0') : "") :
        "." * ext_date
    return string(case, ".", compname, ".h", hidx, suffix, ".nc")
end

# ==========================================================================
# HistoryTapeSet — the multi-tape container (h0/h1/h2…).
# ==========================================================================

"""
    HistoryTapeSet

A set of history tapes (h0, h1, …), each with its own field list and an
output frequency `nhtfrq` (Fortran semantics: `0` = monthly/event-driven,
`>0` = every N accumulate steps, `<0` = every |N| hours — interpreted by the
caller). Built from a master field list via `build_history_tapes`.

Fields:
- `tapes`   — the per-tape `HistoryTape`s (index 1 = h0).
- `nhtfrq`  — per-tape output frequency.
- `case`    — case id used for filenames.
- `compname`— component name (default `"clm2"`).
"""
Base.@kwdef mutable struct HistoryTapeSet
    tapes::Vector{HistoryTape} = HistoryTape[]
    nhtfrq::Vector{Int}        = Int[]
    case::String               = "clmrun"
    compname::String           = "clm2"
end

"""
    build_history_tapes(; master=default_master_field_list(),
                          fincl=[String[], …], fexcl=[String[], …],
                          nhtfrq=[0, …], avgflag_pertape=["", …],
                          case="clmrun", compname="clm2",
                          hist_empty_htapes=false) -> HistoryTapeSet

Construct a `HistoryTapeSet`, mirroring Fortran `htapes_build` /
`htapes_fieldlist`. For each tape `t`:

- A master field is added when it is named in `fincl[t]` (even if not
  on-by-default), OR — unless `hist_empty_htapes` — when it is on-by-default
  and not named in `fexcl[t]`.
- The avgflag is resolved with this precedence (highest first):
  1. a per-field `"NAME:FLAG"` override in that tape's `fincl[t]`,
  2. the tape-wide `avgflag_pertape[t]` (if non-empty),
  3. the field's master-list default avgflag.
- `fincl`/`fexcl` entries naming a field absent from `master` raise an error
  (Fortran validates fincl/fexcl against `allhistfldlist`).

`fincl`, `fexcl`, `nhtfrq`, and `avgflag_pertape` are per-tape lists; the
number of tapes is the longest of these (default 1 tape = h0). Shorter lists
are padded (empty include/exclude, `nhtfrq=0`, no tape-wide flag).

`hist_mfilt` (per tape) is the max number of time samples per history file
(Fortran namelist `hist_mfilt`, default 1 for h0). `hist_ndens` (per tape) is
the output precision (`hist_ndens`, default 2 = single/Float32; 1 = double).
"""
function build_history_tapes(; master::Vector{HistFieldSpec}=default_master_field_list(),
                               fincl::AbstractVector=Vector{String}[],
                               fexcl::AbstractVector=Vector{String}[],
                               nhtfrq::AbstractVector{<:Integer}=Int[],
                               avgflag_pertape::AbstractVector{<:AbstractString}=String[],
                               hist_mfilt::AbstractVector{<:Integer}=Int[],
                               hist_ndens::AbstractVector{<:Integer}=Int[],
                               case::AbstractString="clmrun",
                               compname::AbstractString="clm2",
                               hist_empty_htapes::Bool=false)
    ntapes = max(length(fincl), length(fexcl), length(nhtfrq),
                 length(avgflag_pertape), length(hist_mfilt),
                 length(hist_ndens), 1)

    _get(v, t, default) = t <= length(v) ? v[t] : default

    # Validate every fincl/fexcl name against the master list (Fortran does
    # this up front and endruns on an unknown name).
    masternames = Set(s.name for s in master)
    for t in 1:ntapes
        for entry in _get(fincl, t, String[])
            nm = hist_getname(entry)
            isempty(nm) && continue
            nm in masternames ||
                error("build_history_tapes: fincl field '$(nm)' (tape $(t)) not in master list")
        end
        for entry in _get(fexcl, t, String[])
            nm = hist_getname(entry)
            isempty(nm) && continue
            nm in masternames ||
                error("build_history_tapes: fexcl field '$(nm)' (tape $(t)) not in master list")
        end
    end

    tapes = HistoryTape[]
    freqs = Int[]
    for t in 1:ntapes
        tinc = _get(fincl, t, String[])
        texc = _get(fexcl, t, String[])
        tflag = String(_get(avgflag_pertape, t, ""))

        # name → fincl override flag (last occurrence wins, as in a scan)
        inc_flag = Dict{String,String}()
        for entry in tinc
            nm = hist_getname(entry)
            isempty(nm) && continue
            inc_flag[nm] = hist_getflag(entry)
        end
        exc_names = Set(hist_getname(e) for e in texc if !isempty(hist_getname(e)))

        # Default mfilt mirrors Fortran: 1 for the first tape, 30 for the rest.
        tape = HistoryTape(mfilt = Int(_get(hist_mfilt, t, t == 1 ? 1 : 30)),
                           ndens = Int(_get(hist_ndens, t, 2)))
        for spec in master
            included = haskey(inc_flag, spec.name)
            if included
                flag = inc_flag[spec.name]
            elseif !hist_empty_htapes && spec.on_by_default && !(spec.name in exc_names)
                flag = ""
            else
                continue
            end
            # avgflag precedence: per-field override > per-tape > field default
            avg = !isempty(flag) ? flag :
                  !isempty(tflag) ? tflag :
                  spec.avgflag
            hist_addfld!(tape, spec.name, spec.units, spec.getter;
                         avgflag=avg, long_name=spec.long_name,
                         level=spec.level, nlev=spec.nlev, levdim=spec.levdim)
        end
        push!(tapes, tape)
        push!(freqs, Int(_get(nhtfrq, t, 0)))
    end

    return HistoryTapeSet(tapes=tapes, nhtfrq=freqs,
                          case=String(case), compname=String(compname))
end

"""
    hist_accumulate!(set::HistoryTapeSet, inst)

Accumulate `inst` into every tape in the set.
"""
function hist_accumulate!(set::HistoryTapeSet, inst)
    for tape in set.tapes
        hist_accumulate!(tape, inst)
    end
    return nothing
end

"""
    hist_write!(set::HistoryTapeSet; ext_date="", reset=true, kwargs...) -> Vector{String}

Write every tape in the set, honoring each tape's `mfilt` (max time samples per
file). A tape keeps appending time records to its current file until that file
holds `mfilt` records; the next write then rolls over to a fresh file
(`set_hist_filename` with an incremented file index). Files land under the
current directory unless `dir` is given.

When `ext_date` is given it stamps the filename (`<case>.<compname>.h<t-1>.
<ext_date>.nc`); otherwise rolled-over files are distinguished by a trailing
`-NNNNN` file-index suffix (the first file has none). Tapes with no accumulated
samples are skipped (returns only the files actually written). Extra `kwargs`
(e.g. `time`, `time_units`, `calendar`) are forwarded to per-tape
`hist_write!`; `append`/`begtime` are managed here per the rollover state.
"""
function hist_write!(set::HistoryTapeSet; ext_date::AbstractString="",
                     dir::AbstractString="", reset::Bool=true, kwargs...)
    written = String[]
    for (t, tape) in enumerate(set.tapes)
        (isempty(tape.fields) || tape.nsamples == 0) && continue

        # Roll over to a new file if the current one is full (or none yet open).
        if isempty(tape.cur_filename) || tape.ntimes >= tape.mfilt
            if !isempty(tape.cur_filename)
                tape.file_index += 1   # advance only on an actual rollover
            end
            fn = set_hist_filename(set.case, t; ext_date=ext_date,
                                   compname=set.compname,
                                   file_index=tape.file_index)
            isempty(dir) || (fn = joinpath(dir, fn))
            tape.cur_filename = fn
            tape.ntimes = 0
        end

        # Append within the current file (the very first record creates it).
        appending = isfile(tape.cur_filename) && tape.ntimes > 0
        hist_write!(tape, tape.cur_filename; reset=reset, append=appending, kwargs...)
        push!(written, tape.cur_filename)
    end
    return unique(written)
end
