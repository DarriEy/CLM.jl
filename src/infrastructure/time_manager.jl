# ==========================================================================
# Ported from: src/main/clm_time_manager.F90 (2086 lines → ~80 lines)
# Minimal time management using Julia Dates stdlib
# ==========================================================================

# --------------------------------------------------------------------------
# Calendar selector + average days-per-year
# Ported from clm_time_manager.F90: the module-level `calendar` variable, the
# NO_LEAP_C / GREGORIAN_C parameters, and `get_average_days_per_year()`.
#
# `get_average_days_per_year` gives the constant days-per-year used to convert
# per-year parameters (e.g. decomposition turnover times tau_* → rate constants
# k_*) into per-second values, so a parameter has one fixed value rather than a
# slightly different one on leap vs. non-leap years. Calendar-dependent:
#   NO_LEAP   → 365.0     (CTSM default; keeps the port byte-identical)
#   GREGORIAN → 365.2425  (average once the 4/100/400 leap rule is applied)
# --------------------------------------------------------------------------

const NO_LEAP_C   = "NO_LEAP"
const GREGORIAN_C = "GREGORIAN"

const DAYS_PER_YEAR_NOLEAP    = 365.0
const DAYS_PER_YEAR_GREGORIAN = 365.2425

"""
Module-level calendar in effect for average-days-per-year conversions, mirroring
`clm_time_manager.F90`'s `calendar` module variable. Defaults to `NO_LEAP_C`
(the CTSM default → 365.0 days/yr), so the CN decomposition rate constants stay
byte-identical to the historical port. Set to `GREGORIAN_C` for a Gregorian run.
"""
const CLM_CALENDAR = Ref{String}(NO_LEAP_C)

"""
    get_average_days_per_year(calendar::AbstractString = CLM_CALENDAR[]) -> Float64

Average number of days per year for the given calendar (mirrors Fortran
`get_average_days_per_year`): 365.0 for `NO_LEAP`, 365.2425 for `GREGORIAN`.
Errors on an unrecognized calendar. With no argument it reads the module-level
`CLM_CALENDAR[]` (defaults to `NO_LEAP` → 365.0).
"""
function get_average_days_per_year(calendar::AbstractString = CLM_CALENDAR[])
    cal = uppercase(String(calendar))
    if cal == NO_LEAP_C
        return DAYS_PER_YEAR_NOLEAP
    elseif cal == GREGORIAN_C
        return DAYS_PER_YEAR_GREGORIAN
    else
        error("get_average_days_per_year: unrecognized calendar = $(calendar)")
    end
end

"""
    TimeManager

Minimal time manager for CLM.jl. Tracks simulation time using Julia DateTime.
"""
Base.@kwdef mutable struct TimeManager
    start_date::DateTime    = DateTime(2000, 1, 1)
    current_date::DateTime  = DateTime(2000, 1, 1)
    dtime::Int              = 1800          # timestep in seconds
    nstep::Int              = 0             # current timestep number
    calendar::String        = "NO_LEAP"     # calendar type (only NO_LEAP supported)
end

"""
    advance_timestep!(tm::TimeManager)

Advance the time manager by one timestep.
"""
function advance_timestep!(tm::TimeManager)
    tm.nstep += 1
    tm.current_date += Second(tm.dtime)
    nothing
end

"""
    get_curr_date(tm::TimeManager) -> (yr, mon, day, tod)

Return current date as (year, month, day, time-of-day in seconds).
"""
function get_curr_date(tm::TimeManager)
    dt = tm.current_date
    yr = year(dt)
    mon = month(dt)
    d = day(dt)
    tod = hour(dt) * 3600 + minute(dt) * 60 + second(dt)
    return (yr, mon, d, tod)
end

"""
    get_curr_calday(tm::TimeManager) -> Float64

Return current calendar day (1.0 = Jan 1 00:00, etc.) for NO_LEAP calendar.
"""
function get_curr_calday(tm::TimeManager)
    dt = tm.current_date
    doy = dayofyear(dt)
    frac = (hour(dt) * 3600 + minute(dt) * 60 + second(dt)) / SECSPDAY
    return Float64(doy) + frac
end

"""
    is_beg_curr_day(tm::TimeManager) -> Bool

Return true if we are at the beginning of a day (time-of-day == 0).
"""
function is_beg_curr_day(tm::TimeManager)
    dt = tm.current_date
    return hour(dt) == 0 && minute(dt) == 0 && second(dt) == 0
end

"""
    is_end_curr_day(tm::TimeManager) -> Bool

Return true if the just-completed timestep ended at a day boundary.
"""
function is_end_curr_day(tm::TimeManager)
    return tm.nstep > 0 && is_beg_curr_day(tm)
end

"""
    is_beg_curr_year(tm::TimeManager) -> Bool

Return true if we are at Jan 1 00:00.
"""
function is_beg_curr_year(tm::TimeManager)
    dt = tm.current_date
    return month(dt) == 1 && day(dt) == 1 && is_beg_curr_day(tm)
end

"""
    get_nstep(tm::TimeManager) -> Int

Return current timestep number.
"""
get_nstep(tm::TimeManager) = tm.nstep

"""
    get_average_days_per_year(tm::TimeManager) -> Float64

Average days per year for the calendar configured on `tm` (see the
`AbstractString` method). NO_LEAP → 365.0, GREGORIAN → 365.2425.
"""
get_average_days_per_year(tm::TimeManager) = get_average_days_per_year(tm.calendar)
