# FatesGlobals.jl
# Julia port of FATES src/fates/main/FatesGlobals.F90
#
# Global data + logging used inside FATES: fates_log, fates_endrun, warning
# handling, and string-formatting helpers. Logging maps to Julia @warn/error.
# Deps: FatesConstantsMod (fates_r8).

# ---------------------------------------------------------------------------
# Module-level globals (Fortran module variables -> Refs so they are mutable)
# ---------------------------------------------------------------------------
const fates_log_ = Ref{Int}(6)                 # log unit (stdout-ish default)
const fates_global_verbose_ = Ref{Bool}(false)

# Warning handling: stop writing the same warning over and over.
const max_ids = 200       # maximum number of unique warning ids
const max_warnings = 100  # maximum warnings before we stop writing
const warning_override = false  # set true to suppress all of these warnings

# warn_counts/warn_active indexed 0:max_ids in Fortran -> store as Vector with a
# +1 offset accessor (index 0 -> slot 1).
const warn_counts = zeros(Int, max_ids + 1)
const warn_active = trues(max_ids + 1)

# =====================================================================================

function FatesGlobalsInit(log_unit::Integer, global_verbose::Bool)
    fates_log_[] = log_unit
    fates_global_verbose_[] = global_verbose
    # call TwoStreamLogInit(log_unit) -- no-op. The two-stream code IS ported
    # (TwoStreamMLPEMod.jl / FatesTwoStreamUtilsMod.jl); it just has no separate
    # log unit to initialize in this port.
    return nothing
end

# =====================================================================================

fates_log() = fates_log_[]

fates_global_verbose() = fates_global_verbose_[]

"""
    fates_endrun(msg)

Abort the model for abnormal termination (derived from CLM's endrun_vanilla()).
Maps `shr_sys_abort` to a Julia `error`.
"""
function fates_endrun(msg::AbstractString)
    error("ENDRUN: " * msg)
end

# =====================================================================================

"""
    FatesWarn(msg; index=0)

Emit a warning, but stop emitting once a given id has saturated `max_warnings`.
"""
function FatesWarn(msg::AbstractString; index::Integer=0)
    warning_override && return nothing  # exit early if warnings are off

    ind = Int(index)
    # warn_counts/warn_active are 0-based in Fortran; offset by +1.
    slot = ind + 1
    warn_counts[slot] += 1

    if warn_active[slot]
        @warn "FATESWARN: " * strip(I2S(ind)) * " m: " * msg
        if warn_counts[slot] > max_warnings
            warn_active[slot] = false
            @warn "FATESWARN: " * strip(I2S(ind)) * " has saturated messaging, no longer reporting"
        end
    end
    return nothing
end

# =====================================================================================

function FatesReportTotalWarnings()
    for ind in 0:max_ids
        slot = ind + 1
        if warn_counts[slot] > 0
            @warn "FATESWARN: " * strip(I2S(ind)) * " was triggered " *
                  strip(I2S(warn_counts[slot])) * " times"
        end
    end
    return nothing
end

# =====================================================================================
# String formatting helpers (Fortran write-to-string idioms)
# =====================================================================================

"""
    N2S(real_in) -> String

Format a real in a scientific-notation string (Fortran `write(str,'(E12.6)')`).
We avoid a Printf dependency; the value is the load-bearing part, not the exact
field width.
"""
function N2S(real_in::Real)
    return string(float(real_in))
end

"""
    I2S(int_in) -> String

Format an integer as a string, right-justified to width 15 (Fortran `'(I15)'`).
"""
function I2S(int_in::Integer)
    return lpad(string(int_in), 15)
end

"""
    A2S(reals_in) -> String

Concatenate an array of reals into a comma-separated string.
"""
function A2S(reals_in::AbstractVector{<:Real})
    str = ", "
    for i in 1:length(reals_in)
        str = strip(str) * ", " * N2S(reals_in[i])
    end
    return str
end
