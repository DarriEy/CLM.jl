# ==========================================================================
# Ported from: src/main/accumulMod.F90
# Accumulation fields: time averages, running means, running accumulations
#
# Three types of accumulations:
#   - timeavg:  Average over time interval (valid only at end of interval)
#   - runmean:  Running mean over time interval (valid once simulation
#               exceeds the averaging interval)
#   - runaccum: Running accumulation (continuously accumulated; reset by
#               markreset_accum_field!)
# ==========================================================================

# --------------------------------------------------------------------------
# Accumulation type constants
# --------------------------------------------------------------------------

const ACCTYPE_TIMEAVG  = 1
const ACCTYPE_RUNMEAN  = 2
const ACCTYPE_RUNACCUM = 3

# --------------------------------------------------------------------------
# AccumField — per-field data structure
# --------------------------------------------------------------------------

"""
    AccumField

One accumulated field. Holds values, reset flags, step counters, and metadata.

Ported from `accum_field` in `accumulMod.F90`.
"""
Base.@kwdef mutable struct AccumField
    name::String = ""                                          # field name
    desc::String = ""                                          # field description
    units::String = ""                                         # field units
    acctype::Int = 0                                           # accumulation type (ACCTYPE_*)
    type1d::String = ""                                        # subgrid type: gridcell/landunit/column/pft
    type2d::String = ""                                        # level type (levsoi, numrad, etc.)
    beg1d::Int = 1                                             # subgrid beginning index
    end1d::Int = 0                                             # subgrid ending index
    numlev::Int = 1                                            # number of vertical levels
    active::Vector{Bool} = Bool[]                              # whether each point is active
    initval::Float64 = 0.0                                     # initial/reset value
    val::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)       # accumulated field [npts × numlev]
    reset::Matrix{Bool} = Matrix{Bool}(undef, 0, 0)           # whether field needs reset [npts × numlev]
    period::Int = 0                                            # accumulation period (in timesteps)
    scale_by_thickness::Bool = false                           # scale vertically interpolated variable by soil thickness
    old_name::String = ""                                      # previous name (for restart compatibility)
    nsteps::Matrix{Int} = Matrix{Int}(undef, 0, 0)            # steps accumulated since last reset [npts × numlev]
    ndays_reset_shifted::Matrix{Int} = Matrix{Int}(undef, 0, 0)  # days that reset has shifted timeavg out of sync
end

# --------------------------------------------------------------------------
# AccumManager — holds the collection of accumulation fields
# --------------------------------------------------------------------------

"""
    AccumManager

Container for all accumulation fields. Replaces the Fortran module-level
`accum` array and `naccflds` counter.
"""
Base.@kwdef mutable struct AccumManager
    fields::Vector{AccumField} = AccumField[]
end

# --------------------------------------------------------------------------
# acctype_to_string
# --------------------------------------------------------------------------

"""
    acctype_to_string(acctype) -> String

Return a string representation of an ACCTYPE parameter.

Ported from `acctype_to_string` in `accumulMod.F90`.
"""
function acctype_to_string(acctype::Int)::String
    if acctype == ACCTYPE_TIMEAVG
        return "timeavg"
    elseif acctype == ACCTYPE_RUNMEAN
        return "runmean"
    elseif acctype == ACCTYPE_RUNACCUM
        return "runaccum"
    else
        return "unknown"
    end
end

# --------------------------------------------------------------------------
# find_field — find field index by name
# --------------------------------------------------------------------------

"""
    find_field(mgr, field_name) -> Int

Find field index given field name. Errors if not found.

Ported from `find_field` in `accumulMod.F90`.
"""
function find_field(mgr::AccumManager, field_name::String)::Int
    for i in eachindex(mgr.fields)
        if mgr.fields[i].name == field_name
            return i
        end
    end
    error("Accumulation field '$(field_name)' not found")
end

# --------------------------------------------------------------------------
# init_accum_field! — initialize a new accumulation field
# --------------------------------------------------------------------------

"""
    init_accum_field!(mgr; name, units, desc, accum_type, accum_period,
                      numlev, subgrid_type, init_value, active, npts,
                      type2d="", scale_by_thickness=false, old_name="",
                      step_size=1800)

Initialize an accumulation field and add it to the manager.

- `accum_type`: one of `"timeavg"`, `"runmean"`, `"runaccum"`
- `accum_period`: if positive, period in timesteps; if negative, period in
  days (converted to timesteps using `step_size`)
- `active`: reference to the active mask for the subgrid type
- `npts`: number of spatial points
- `step_size`: model timestep in seconds (used when `accum_period < 0`)

Ported from `init_accum_field` in `accumulMod.F90`.
"""
function init_accum_field!(mgr::AccumManager;
                           name::String,
                           units::String,
                           desc::String,
                           accum_type::String,
                           accum_period::Int,
                           numlev::Int,
                           subgrid_type::String,
                           init_value::Float64,
                           active::Vector{Bool},
                           npts::Int,
                           type2d::String = "",
                           scale_by_thickness::Bool = false,
                           old_name::String = "",
                           step_size::Int = 1800)

    if numlev > 1
        if type2d == ""
            error("Field '$(name)': 2d accumulation fields require type2d")
        end
    end

    # Create new field
    af = AccumField()

    # Set accumulation type
    if accum_type == "timeavg"
        af.acctype = ACCTYPE_TIMEAVG
    elseif accum_type == "runmean"
        af.acctype = ACCTYPE_RUNMEAN
    elseif accum_type == "runaccum"
        af.acctype = ACCTYPE_RUNACCUM
    else
        error("init_accum_field: unknown accum_type '$(accum_type)'")
    end

    af.name = name
    af.units = units
    af.desc = desc
    af.initval = init_value
    af.old_name = old_name
    af.type1d = subgrid_type
    af.type2d = type2d
    af.scale_by_thickness = scale_by_thickness

    # Convert period from days to timesteps if negative
    period = accum_period
    if period < 0
        period = (-period) * round(Int, SECSPDAY) ÷ step_size
    end
    af.period = period

    af.beg1d = 1
    af.end1d = npts
    af.numlev = numlev
    af.active = active

    # Validate init_value for timeavg and runaccum
    if (af.acctype == ACCTYPE_TIMEAVG || af.acctype == ACCTYPE_RUNACCUM)
        if init_value != 0.0
            error("init_accum_field: init_value must be 0 for timeavg and runaccum fields (field '$(name)')")
        end
    end

    # Allocate and initialize
    af.val = fill(init_value, npts, numlev)
    af.reset = fill(false, npts, numlev)
    af.nsteps = fill(0, npts, numlev)
    af.ndays_reset_shifted = fill(0, npts, numlev)

    push!(mgr.fields, af)
    nothing
end

# --------------------------------------------------------------------------
# print_accum_fields
# --------------------------------------------------------------------------

"""
    print_accum_fields(mgr)

Diagnostic printout of accumulated fields.

Ported from `print_accum_fields` in `accumulMod.F90`.
"""
function print_accum_fields(mgr::AccumManager)
    println()
    println("Initializing variables for time accumulation .....")
    println(repeat("-", 60))
    println("Accumulated fields")
    println(" No Name     Units    Type    Period Inival Description")
    println(repeat("_", 71))
    for (nf, af) in enumerate(mgr.fields)
        if af.period != typemax(Int)
            println(" $(lpad(nf, 2)) $(rpad(af.name, 8)) $(rpad(af.units, 8)) " *
                    "$(rpad(acctype_to_string(af.acctype), 8)) $(lpad(af.period, 5)) " *
                    "$(lpad(round(Int, af.initval), 4)) $(af.desc)")
        else
            println(" $(lpad(nf, 2)) $(rpad(af.name, 8)) $(rpad(af.units, 8)) " *
                    "$(rpad(acctype_to_string(af.acctype), 8))  N.A. " *
                    "$(lpad(round(Int, af.initval), 4)) $(af.desc)")
        end
    end
    println(repeat("_", 71))
    println()
    println(repeat("-", 60))
    println("Successfully initialized variables for accumulation")
    println()
    nothing
end

# --------------------------------------------------------------------------
# Internal extract helpers
# --------------------------------------------------------------------------

"""
    extract_accum_field_basic!(af, level, nstep, field)

Extract values for one level — used for runmean and runaccum fields.
Simply copies the accumulated value.

Ported from `extract_accum_field_basic` in `accumulMod.F90`.
"""
function extract_accum_field_basic!(af::AccumField, level::Int, nstep::Int,
                                    field::AbstractVector{Float64})
    begi = af.beg1d
    endi = af.end1d
    @assert length(field) == endi - begi + 1

    for k in begi:endi
        kf = k - begi + 1
        field[kf] = af.val[k, level]
    end
    nothing
end

"""
    extract_accum_field_timeavg!(af, level, nstep, field)

Extract values for one level of a timeavg field. Returns the accumulated
value at the end of each averaging period; returns `SPVAL` at other times.

Ported from `extract_accum_field_timeavg` in `accumulMod.F90`.
"""
function extract_accum_field_timeavg!(af::AccumField, level::Int, nstep::Int,
                                      field::AbstractVector{Float64})
    begi = af.beg1d
    endi = af.end1d
    @assert length(field) == endi - begi + 1

    for k in begi:endi
        kf = k - begi + 1
        effective_nstep = nstep - af.ndays_reset_shifted[k, level]
        if mod(effective_nstep, af.period) == 0
            field[kf] = af.val[k, level]
        else
            field[kf] = SPVAL
        end
    end
    nothing
end

# --------------------------------------------------------------------------
# extract_accum_field! — public interface (single-level)
# --------------------------------------------------------------------------

"""
    extract_accum_field!(mgr, name, field::AbstractVector{Float64}, nstep)

Extract single-level accumulated field values.

Ported from `extract_accum_field_sl` in `accumulMod.F90`.
"""
function extract_accum_field!(mgr::AccumManager, name::String,
                              field::AbstractVector{Float64}, nstep::Int)
    nf = find_field(mgr, name)
    af = mgr.fields[nf]
    @assert af.numlev == 1

    if af.acctype == ACCTYPE_TIMEAVG
        extract_accum_field_timeavg!(af, 1, nstep, field)
    else
        extract_accum_field_basic!(af, 1, nstep, field)
    end
    nothing
end

# --------------------------------------------------------------------------
# extract_accum_field! — public interface (multi-level)
# --------------------------------------------------------------------------

"""
    extract_accum_field!(mgr, name, field::AbstractMatrix{Float64}, nstep)

Extract multi-level accumulated field values.

Ported from `extract_accum_field_ml` in `accumulMod.F90`.
"""
function extract_accum_field!(mgr::AccumManager, name::String,
                              field::AbstractMatrix{Float64}, nstep::Int)
    nf = find_field(mgr, name)
    af = mgr.fields[nf]
    @assert size(field, 2) == af.numlev

    for j in 1:af.numlev
        field_j = @view field[:, j]
        if af.acctype == ACCTYPE_TIMEAVG
            extract_accum_field_timeavg!(af, j, nstep, field_j)
        else
            extract_accum_field_basic!(af, j, nstep, field_j)
        end
    end
    nothing
end

# --------------------------------------------------------------------------
# Internal update helpers
# --------------------------------------------------------------------------

"""
    update_accum_field_timeavg!(af, level, nstep, field)

Update values for one level of a timeavg field:
1. Reset at start of new accumulation period
2. Accumulate field values
3. Normalize at end of period

Ported from `update_accum_field_timeavg` in `accumulMod.F90`.
"""
function update_accum_field_timeavg!(af::AccumField, level::Int, nstep::Int,
                                     field::AbstractVector{Float64})
    begi = af.beg1d
    endi = af.end1d
    @assert length(field) == endi - begi + 1

    # Phase 1: reset at start of new accumulation period
    for k in begi:endi
        effective_nstep = nstep - af.ndays_reset_shifted[k, level]
        time_to_reset = (mod(effective_nstep, af.period) == 1 || af.period == 1) && effective_nstep != 0
        if af.active[k] && (time_to_reset || af.reset[k, level])
            if af.reset[k, level] && !time_to_reset
                af.ndays_reset_shifted[k, level] += af.nsteps[k, level]
            end
            af.val[k, level] = af.initval
            af.nsteps[k, level] = 0
            af.reset[k, level] = false
        end
    end

    # Ignore reset requests that occurred when patch was inactive
    af.reset[begi:endi, level] .= false

    # Phase 2: accumulate
    for k in begi:endi
        if af.active[k]
            kf = k - begi + 1
            af.val[k, level] += field[kf]
            af.nsteps[k, level] += 1
        end
    end

    # Phase 3: normalize at end of period
    for k in begi:endi
        effective_nstep = nstep - af.ndays_reset_shifted[k, level]
        if af.active[k] && mod(effective_nstep, af.period) == 0
            af.val[k, level] /= af.nsteps[k, level]
        end
    end

    nothing
end

"""
    update_accum_field_runmean!(af, level, nstep, field)

Update values for one level of a runmean field. Computes a running mean
with the accumulation period as the window size. The running mean is valid
once nsteps reaches the accumulation period.

Note: unlike other accumulator types, runmean preserves reset requests
that occurred when the patch was inactive.

Ported from `update_accum_field_runmean` in `accumulMod.F90`.
"""
function update_accum_field_runmean!(af::AccumField, level::Int, nstep::Int,
                                     field::AbstractVector{Float64})
    begi = af.beg1d
    endi = af.end1d
    @assert length(field) == endi - begi + 1

    for k in begi:endi
        if af.active[k]
            if af.reset[k, level]
                af.nsteps[k, level] = 0
                af.val[k, level] = af.initval
                af.reset[k, level] = false
            end
            kf = k - begi + 1
            af.nsteps[k, level] += 1
            # Cap nsteps at period
            af.nsteps[k, level] = min(af.nsteps[k, level], af.period)
            accper = af.nsteps[k, level]
            af.val[k, level] = ((accper - 1) * af.val[k, level] + field[kf]) / accper
        end
    end

    nothing
end

"""
    update_accum_field_runaccum!(af, level, nstep, field)

Update values for one level of a runaccum field. Continuously accumulates,
clamped to [0, 99999]. Reset sets the value to zero.

Note: unlike other accumulator types, runaccum cannot reset AND update
in the same call.

Ported from `update_accum_field_runaccum` in `accumulMod.F90`.
"""
function update_accum_field_runaccum!(af::AccumField, level::Int, nstep::Int,
                                      field::AbstractVector{Float64})
    begi = af.beg1d
    endi = af.end1d
    @assert length(field) == endi - begi + 1

    for k in begi:endi
        if af.active[k]
            kf = k - begi + 1
            if af.reset[k, level]
                # Reset: cannot update in same call
                af.val[k, level] = 0.0
                af.nsteps[k, level] = trunc(Int, af.initval)
                af.reset[k, level] = false
            else
                af.val[k, level] = min(max(af.val[k, level] + field[kf], 0.0), 99999.0)
                af.nsteps[k, level] += 1
            end
        end
    end

    # Ignore reset requests that occurred when patch was inactive
    af.reset[begi:endi, level] .= false

    nothing
end

# --------------------------------------------------------------------------
# update_accum_field! — public interface (single-level)
# --------------------------------------------------------------------------

"""
    update_accum_field!(mgr, name, field::AbstractVector{Float64}, nstep)

Update single-level accumulation field with new values.

Values of `field` are ignored in inactive points.

Ported from `update_accum_field_sl` in `accumulMod.F90`.
"""
function update_accum_field!(mgr::AccumManager, name::String,
                             field::AbstractVector{Float64}, nstep::Int)
    nf = find_field(mgr, name)
    af = mgr.fields[nf]
    @assert af.numlev == 1

    if af.acctype == ACCTYPE_TIMEAVG
        update_accum_field_timeavg!(af, 1, nstep, field)
    elseif af.acctype == ACCTYPE_RUNMEAN
        update_accum_field_runmean!(af, 1, nstep, field)
    else
        update_accum_field_runaccum!(af, 1, nstep, field)
    end
    nothing
end

# --------------------------------------------------------------------------
# update_accum_field! — public interface (multi-level)
# --------------------------------------------------------------------------

"""
    update_accum_field!(mgr, name, field::AbstractMatrix{Float64}, nstep)

Update multi-level accumulation field with new values.

Values of `field` are ignored in inactive points.

Ported from `update_accum_field_ml` in `accumulMod.F90`.
"""
function update_accum_field!(mgr::AccumManager, name::String,
                             field::AbstractMatrix{Float64}, nstep::Int)
    nf = find_field(mgr, name)
    af = mgr.fields[nf]
    @assert size(field, 2) == af.numlev

    for j in 1:af.numlev
        field_j = @view field[:, j]
        if af.acctype == ACCTYPE_TIMEAVG
            update_accum_field_timeavg!(af, j, nstep, field_j)
        elseif af.acctype == ACCTYPE_RUNMEAN
            update_accum_field_runmean!(af, j, nstep, field_j)
        else
            update_accum_field_runaccum!(af, j, nstep, field_j)
        end
    end
    nothing
end

# --------------------------------------------------------------------------
# markreset_accum_field!
# --------------------------------------------------------------------------

"""
    markreset_accum_field!(mgr, name; kf=nothing, level=nothing)

Mark accumulator values as needing to be reset. The actual reset happens
in the next `update_accum_field!` call.

- No `kf`, no `level`: reset all points, all levels
- `kf` only: reset all levels for point `kf`
- `level` only: reset all points for level `level`
- Both: reset one point at one level

Ported from `markreset_accum_field` in `accumulMod.F90`.
"""
function markreset_accum_field!(mgr::AccumManager, name::String;
                                kf::Union{Int, Nothing} = nothing,
                                level::Union{Int, Nothing} = nothing)
    nf = find_field(mgr, name)
    af = mgr.fields[nf]
    begi = af.beg1d
    endi = af.end1d
    numlev = af.numlev

    if kf !== nothing
        k = kf + begi - 1
        @assert k >= begi && k <= endi
        if level !== nothing
            af.reset[k, level] = true
        else
            af.reset[k, 1:numlev] .= true
        end
    elseif level !== nothing
        af.reset[begi:endi, level] .= true
    else
        af.reset[begi:endi, 1:numlev] .= true
    end
    nothing
end

# --------------------------------------------------------------------------
# get_accum_reset
# --------------------------------------------------------------------------

"""
    get_accum_reset(mgr, name) -> Matrix{Bool}

Get a copy of the reset array for an accumulator field.

Ported from `get_accum_reset` in `accumulMod.F90`.
"""
function get_accum_reset(mgr::AccumManager, name::String)::Matrix{Bool}
    nf = find_field(mgr, name)
    return copy(mgr.fields[nf].reset)
end

# --------------------------------------------------------------------------
# accum_rest! — restart stub
# --------------------------------------------------------------------------

"""
    accum_rest!(mgr)

Stub for restart variable read/write. No-op in Julia port.

Ported from `accumulRest` in `accumulMod.F90`.
"""
function accum_rest!(mgr::AccumManager)
    nothing
end

# --------------------------------------------------------------------------
# clean_accum_fields!
# --------------------------------------------------------------------------

"""
    clean_accum_fields!(mgr)

Deallocate space and reset the accumulation fields list.

Ported from `clean_accum_fields` in `accumulMod.F90`.
"""
function clean_accum_fields!(mgr::AccumManager)
    empty!(mgr.fields)
    nothing
end
