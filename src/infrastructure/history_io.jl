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
# This is a lighter, more generic complement to history_writer.jl (which
# does gridcell area-weighted daily aggregation). Here fields are written on
# their native subgrid level, matching the surfdata.jl NCDatasets style.
# ==========================================================================

"""
    HistField

One registered history field: a `name`, `units`, a `getter` mapping
`CLMInstances` to an `AbstractVector` of values, and an averaging flag
(`"A"`, `"I"`, `"X"`, or `"M"`).

`accum` is the running accumulator and `nacs` the per-element sample count.
"""
mutable struct HistField
    name::String
    units::String
    long_name::String
    getter::Function
    avgflag::String
    accum::Vector{Float64}
    nacs::Vector{Int}
end

"""
    HistoryTape

A single time-averaged history tape: an ordered list of registered fields,
the spatial dimension name/length they share, and the total number of
samples taken since the last write.
"""
Base.@kwdef mutable struct HistoryTape
    fields::Vector{HistField}   = HistField[]
    nelem::Int                  = 0       # spatial dim length (0 = infer on first accumulate)
    dimname::String             = "subgrid"
    nsamples::Int               = 0
end

const _HIST_VALID_AVGFLAGS = ("A", "I", "X", "M")

"""
    hist_addfld!(tape, name, units, getter; avgflag="A", long_name=name)

Register a field on `tape`. `getter(inst::CLMInstances)::AbstractVector`
returns the current per-element values. `avgflag` selects the time
aggregation: `"A"` average, `"I"` instantaneous, `"X"` max, `"M"` min.
"""
function hist_addfld!(tape::HistoryTape, name::AbstractString, units::AbstractString,
                      getter::Function; avgflag::AbstractString="A",
                      long_name::AbstractString=name)
    avgflag in _HIST_VALID_AVGFLAGS ||
        error("hist_addfld!: unknown averaging flag $(avgflag); valid: $(_HIST_VALID_AVGFLAGS)")
    any(f -> f.name == name, tape.fields) &&
        error("hist_addfld!: field $(name) already registered")
    push!(tape.fields, HistField(String(name), String(units), String(long_name),
                                 getter, String(avgflag), Float64[], Int[]))
    return tape
end

# Initialize a field's accumulator to the natural identity for its avgflag.
function _hist_reset_field!(fld::HistField, n::Int)
    fld.accum = fill(_hist_init_value(fld.avgflag), n)
    fld.nacs  = zeros(Int, n)
    return nothing
end

_hist_init_value(avgflag::String) =
    avgflag == "X" ? -Inf :
    avgflag == "M" ?  Inf :
    0.0  # "A" sum, "I" last

"""
    hist_accumulate!(tape, inst)

Sample every registered field from `inst` and fold it into the accumulators
according to each field's avgflag. Values equal to `SPVAL` or non-finite are
skipped for that element (no sample counted), matching Fortran behavior.
"""
function hist_accumulate!(tape::HistoryTape, inst)
    for fld in tape.fields
        data = fld.getter(inst)
        n = length(data)
        n == 0 && continue

        # Lazily fix the tape's spatial length on first real sample.
        if tape.nelem == 0
            tape.nelem = n
        end
        if length(fld.accum) != n
            _hist_reset_field!(fld, n)
        end

        @inbounds for k in 1:n
            v = Float64(data[k])
            (isfinite(v) && v != SPVAL) || continue
            if fld.avgflag == "I"
                fld.accum[k] = v
                fld.nacs[k] = 1
            elseif fld.avgflag == "A"
                fld.accum[k] += v
                fld.nacs[k] += 1
            elseif fld.avgflag == "X"
                fld.accum[k] = max(fld.accum[k], v)
                fld.nacs[k] = 1
            elseif fld.avgflag == "M"
                fld.accum[k] = min(fld.accum[k], v)
                fld.nacs[k] = 1
            end
        end
    end
    tape.nsamples += 1
    return nothing
end

"""
    hist_field_value(fld) -> Vector{Float64}

Finalize a single field's accumulated value: for `"A"` divide by the sample
count; otherwise return the running value. Elements that never received a
valid sample are set to `SPVAL`.
"""
function hist_field_value(fld::HistField)
    n = length(fld.accum)
    out = fill(SPVAL, n)
    @inbounds for k in 1:n
        if fld.nacs[k] > 0
            out[k] = fld.avgflag == "A" ? fld.accum[k] / fld.nacs[k] : fld.accum[k]
        end
    end
    return out
end

"""
    hist_write!(tape, filename; time=0.0, time_units="days since 2000-01-01",
                calendar="noleap", append=false, reset=true)

Write the time-aggregated tape to a CLM-style h0 NetCDF `filename`. Defines a
`time` dimension (unlimited) and the tape's spatial dimension, appends one
time record with each field's finalized value, and (by default) resets the
accumulators for the next interval.

If `append=true` and the file exists, the record is appended to the existing
file; otherwise the file is (re)created.
"""
function hist_write!(tape::HistoryTape, filename::AbstractString;
                     time::Real=0.0,
                     time_units::AbstractString="days since 2000-01-01",
                     calendar::AbstractString="noleap",
                     append::Bool=false,
                     reset::Bool=true)
    isempty(tape.fields) && error("hist_write!: no fields registered")
    tape.nsamples == 0 && error("hist_write!: no samples accumulated")

    n = tape.nelem
    creating = !(append && isfile(filename))

    mode = creating ? "c" : "a"
    NCDataset(filename, mode) do ds
        if creating
            defDim(ds, tape.dimname, n)
            defDim(ds, "time", Inf)  # unlimited
            defVar(ds, "time", Float64, ("time",);
                   attrib = Dict("units" => time_units,
                                 "calendar" => calendar,
                                 "long_name" => "time"))
            for fld in tape.fields
                defVar(ds, fld.name, Float64, (tape.dimname, "time");
                       attrib = Dict("long_name" => fld.long_name,
                                     "units" => fld.units,
                                     "cell_methods" => _hist_cell_method(fld.avgflag),
                                     "_FillValue" => SPVAL))
            end
        end

        ti = ds.dim["time"] + 1
        ds["time"][ti] = Float64(time)
        for fld in tape.fields
            ds[fld.name][:, ti] = hist_field_value(fld)
        end
    end

    if reset
        for fld in tape.fields
            if !isempty(fld.accum)
                _hist_reset_field!(fld, length(fld.accum))
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

"""
    default_history_tape() -> HistoryTape

A starter tape registering the key lnd2atm / energy / water diagnostics on
their native subgrid level. Useful as a quick-start h0 set; callers can
`hist_addfld!` more fields as needed.
"""
function default_history_tape()
    tape = HistoryTape()
    for spec in _MASTER_HIST_FIELDS
        hist_addfld!(tape, spec.name, spec.units, spec.getter;
                     avgflag=spec.avgflag, long_name=spec.long_name)
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
end

"""
    hist_master_field(name, units, getter; avgflag="A", long_name=name,
                      on_by_default=true) -> HistFieldSpec

Build a master-list field spec. Equivalent to Fortran `hist_addfld1d`'s
registration into `allhistfldlist` (the per-tape selection happens later in
`htapes_build`).
"""
function hist_master_field(name::AbstractString, units::AbstractString, getter::Function;
                           avgflag::AbstractString="A", long_name::AbstractString=name,
                           on_by_default::Bool=true)
    avgflag in _HIST_VALID_AVGFLAGS ||
        error("hist_master_field: unknown averaging flag $(avgflag); valid: $(_HIST_VALID_AVGFLAGS)")
    return HistFieldSpec(String(name), String(units), String(long_name),
                         getter, String(avgflag), on_by_default)
end

# The default master field list (the set of always-on h0 diagnostics).
const _MASTER_HIST_FIELDS = HistFieldSpec[
    hist_master_field("EFLX_LH_TOT", "W/m^2",
                      inst -> inst.energyflux.eflx_lh_tot_patch;
                      long_name="total latent heat flux"),
    hist_master_field("FSH", "W/m^2",
                      inst -> inst.energyflux.eflx_sh_tot_patch;
                      long_name="sensible heat flux"),
    hist_master_field("FSA", "W/m^2",
                      inst -> inst.solarabs.fsa_patch;
                      long_name="absorbed solar radiation"),
    hist_master_field("TSA", "K",
                      inst -> inst.temperature.t_ref2m_patch;
                      long_name="2m air temperature"),
    hist_master_field("TG", "K",
                      inst -> inst.temperature.t_grnd_col;
                      long_name="ground temperature"),
    hist_master_field("QRUNOFF", "mm/s",
                      inst -> inst.water.waterfluxbulk_inst.wf.qflx_runoff_col;
                      long_name="total liquid runoff"),
    hist_master_field("QFLX_EVAP_TOT", "mm/s",
                      inst -> inst.water.waterfluxbulk_inst.wf.qflx_evap_tot_col;
                      long_name="total evaporation"),
    hist_master_field("H2OSNO", "mm",
                      history_h2osno_total_col;
                      long_name="snow water equivalent"),
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
                         avgflag=avg, long_name=spec.long_name)
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
