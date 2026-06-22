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
                calendar="noleap", append=false, reset=true, inst=nothing)

Write the time-aggregated tape to a CLM-style h0 NetCDF `filename`. Defines a
`time` dimension (unlimited) and, for each subgrid level present, a spatial
dim sized to that level's element count (or, if `tape.dov2xy`, a single shared
`gridcell` dim). 2d fields additionally define their vertical level dim.
Appends one time record with each field's finalized value and (by default)
resets the accumulators for the next interval.

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
                     inst=nothing)
    isempty(tape.fields) && error("hist_write!: no fields registered")
    tape.nsamples == 0 && error("hist_write!: no samples accumulated")
    tape.dov2xy && inst === nothing &&
        error("hist_write!: dov2xy=true requires `inst` for gridcell aggregation")

    ng = tape.dov2xy ? _hist_ngridcell(inst) : 0

    # Finalize each field's value (and remap to gridcell if dov2xy), and
    # determine the spatial dim name + length it will be written on.
    nf = length(tape.fields)
    finals  = Vector{Any}(undef, nf)
    dimname = Vector{String}(undef, nf)
    dimlen  = Vector{Int}(undef, nf)
    for (i, fld) in enumerate(tape.fields)
        v = hist_field_value(fld)
        if tape.dov2xy
            finals[i]  = hist_dov2xy(fld, v, ng, inst)
            dimname[i] = "gridcell"
            dimlen[i]  = ng
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
            defVar(ds, "time", Float64, ("time",);
                   attrib = Dict("units" => time_units,
                                 "calendar" => calendar,
                                 "long_name" => "time"))
            seen = Set{String}()
            # One spatial dim per distinct level name.
            for i in 1:nf
                if !(dimname[i] in seen)
                    defDim(ds, dimname[i], dimlen[i])
                    push!(seen, dimname[i])
                end
            end
            # Vertical level dims for 2d fields.
            for fld in tape.fields
                if fld.nlev > 1 && !(fld.levdim in seen)
                    defDim(ds, fld.levdim, fld.nlev)
                    push!(seen, fld.levdim)
                end
            end
            # Field variables: 1d -> (space,time); 2d -> (space,lev,time).
            for (i, fld) in enumerate(tape.fields)
                dims = fld.nlev > 1 ? (dimname[i], fld.levdim, "time") :
                                      (dimname[i], "time")
                defVar(ds, fld.name, Float64, dims;
                       attrib = Dict("long_name" => fld.long_name,
                                     "units" => fld.units,
                                     "cell_methods" => _hist_cell_method(fld.avgflag),
                                     "subgrid_level" =>
                                         tape.dov2xy ? "gridcell" : _hist_subgrid_name(fld.level),
                                     "_FillValue" => SPVAL))
            end
        end

        ti = ds.dim["time"] + 1
        ds["time"][ti] = Float64(time)
        for (i, fld) in enumerate(tape.fields)
            if fld.nlev > 1
                ds[fld.name][:, :, ti] = finals[i]
            else
                ds[fld.name][:, ti] = finals[i]
            end
        end
    end

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
    set_hist_filename(case, tape_index; ext_date="", compname="clm2") -> String

Construct a CLM history filename. `tape_index` is 1-based in the code but the
filename uses `h(tape_index-1)`, so tape 1 → `…clm2.h0…`, tape 2 → `…h1…`.
Mirrors Fortran `set_hist_filename` (the date string is opaque here; pass it
via `ext_date`, e.g. `"2000-01"`). Always ends in `.nc`.
"""
function set_hist_filename(case::AbstractString, tape_index::Integer;
                           ext_date::AbstractString="", compname::AbstractString="clm2")
    tape_index >= 1 || error("set_hist_filename: tape_index must be ≥ 1")
    hidx = tape_index - 1
    suffix = isempty(ext_date) ? "" : "." * ext_date
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
"""
function build_history_tapes(; master::Vector{HistFieldSpec}=default_master_field_list(),
                               fincl::AbstractVector=Vector{String}[],
                               fexcl::AbstractVector=Vector{String}[],
                               nhtfrq::AbstractVector{<:Integer}=Int[],
                               avgflag_pertape::AbstractVector{<:AbstractString}=String[],
                               case::AbstractString="clmrun",
                               compname::AbstractString="clm2",
                               hist_empty_htapes::Bool=false)
    ntapes = max(length(fincl), length(fexcl), length(nhtfrq),
                 length(avgflag_pertape), 1)

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

        tape = HistoryTape()
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

Write every tape in the set to its own `set_hist_filename`-named file
(`<case>.<compname>.h<t-1>[.ext_date].nc`) under the current directory unless
`dir` is given. Tapes with no accumulated samples are skipped (returns only
the files actually written). Extra `kwargs` (e.g. `time`, `time_units`,
`calendar`, `append`) are forwarded to per-tape `hist_write!`.
"""
function hist_write!(set::HistoryTapeSet; ext_date::AbstractString="",
                     dir::AbstractString="", reset::Bool=true, kwargs...)
    written = String[]
    for (t, tape) in enumerate(set.tapes)
        (isempty(tape.fields) || tape.nsamples == 0) && continue
        fn = set_hist_filename(set.case, t; ext_date=ext_date, compname=set.compname)
        isempty(dir) || (fn = joinpath(dir, fn))
        hist_write!(tape, fn; reset=reset, kwargs...)
        push!(written, fn)
    end
    return written
end
