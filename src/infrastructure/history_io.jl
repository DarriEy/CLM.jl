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
    # Energy fluxes (patch level)
    hist_addfld!(tape, "EFLX_LH_TOT", "W/m^2",
                 inst -> inst.energyflux.eflx_lh_tot_patch;
                 long_name="total latent heat flux")
    hist_addfld!(tape, "FSH", "W/m^2",
                 inst -> inst.energyflux.eflx_sh_tot_patch;
                 long_name="sensible heat flux")
    hist_addfld!(tape, "FSA", "W/m^2",
                 inst -> inst.solarabs.fsa_patch;
                 long_name="absorbed solar radiation")
    hist_addfld!(tape, "TSA", "K",
                 inst -> inst.temperature.t_ref2m_patch;
                 long_name="2m air temperature")
    # Ground / state (column level)
    hist_addfld!(tape, "TG", "K",
                 inst -> inst.temperature.t_grnd_col;
                 long_name="ground temperature")
    hist_addfld!(tape, "QRUNOFF", "mm/s",
                 inst -> inst.water.waterfluxbulk_inst.wf.qflx_runoff_col;
                 long_name="total liquid runoff")
    hist_addfld!(tape, "QFLX_EVAP_TOT", "mm/s",
                 inst -> inst.water.waterfluxbulk_inst.wf.qflx_evap_tot_col;
                 long_name="total evaporation")
    hist_addfld!(tape, "H2OSNO", "mm",
                 history_h2osno_total_col;
                 long_name="snow water equivalent")
    return tape
end
