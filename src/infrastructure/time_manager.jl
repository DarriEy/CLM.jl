# ==========================================================================
# Ported from: src/main/clm_time_manager.F90 (2086 lines → ~80 lines)
# Minimal time management using Julia Dates stdlib
# ==========================================================================

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

Return true if the next timestep would start a new day.
"""
function is_end_curr_day(tm::TimeManager)
    next_dt = tm.current_date + Second(tm.dtime)
    return day(next_dt) != day(tm.current_date) || month(next_dt) != month(tm.current_date)
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
