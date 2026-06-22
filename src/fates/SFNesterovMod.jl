# SFNesterovMod.jl
# Julia port of FATES src/fates/fire/SFNesterovMod.F90
#
# The Nesterov fire-weather-index concrete implementation. The Fortran defines a
# `nesterov_index` type that *extends* the abstract `fire_weather` type (ported
# in SFFireWeatherMod.jl), overriding the deferred `Init`/`UpdateIndex`
# procedures. The Nesterov index (Nesterov 1968, Eq. 5 of Thonicke et al. 2010)
# accumulates temperature x dewpoint-depression each dry day, and resets to zero
# on any day with precipitation above `min_precip_thresh`.
#
# Julia port: `nesterov_index <: fire_weather` (a mutable struct carrying the
# base fields `fire_weather_index` + `effective_windspeed` that the abstract
# type requires), with method overloads of `init_fire_weather!` (Fortran `Init`)
# and `update_index!` (Fortran `UpdateIndex`). The dewpoint / calc-index helpers
# port verbatim. `fates_r8` -> Float64. Deps: SFFireWeatherMod (fire_weather,
# init_fire_weather!, update_index!), FatesConstantsMod (dewpoint_a, dewpoint_b).

# threshold for precipitation above which to zero NI [mm/day]
const min_precip_thresh = 3.0

"""
    nesterov_index

Concrete [`fire_weather`](@ref) subtype implementing the Nesterov fire-weather
index. Carries the base fire-weather state (`fire_weather_index`,
`effective_windspeed`) required by the abstract type.
"""
Base.@kwdef mutable struct nesterov_index <: fire_weather
    fire_weather_index::Float64  = 0.0  # fire weather index (the Nesterov index)
    effective_windspeed::Float64 = 0.0  # effective wind speed [m/min]
end

"""
    init_fire_weather!(this::nesterov_index)

Initialize the Nesterov index class attributes (Fortran
`init_nesterov_fire_weather`): zero the fire-weather index and effective wind
speed.
"""
function init_fire_weather!(this::nesterov_index)
    this.fire_weather_index  = 0.0
    this.effective_windspeed = 0.0
    return this
end

"""
    update_index!(this::nesterov_index, temp_C, precip, rh, wind)

Update the Nesterov index (Fortran `update_nesterov_index`). If `precip`
(daily precipitation [mm]) exceeds [`min_precip_thresh`](@ref) the index is
rezeroed; otherwise the current day's Nesterov index (temperature times
dewpoint depression) is accumulated.

`temp_C` daily averaged temperature [deg C], `precip` daily precipitation [mm],
`rh` daily relative humidity [%], `wind` daily wind speed [m/min] (unused by the
index update but kept for interface compatibility).
"""
function update_index!(this::nesterov_index, temp_C::Float64, precip::Float64,
                       rh::Float64, wind::Float64)
    if precip > min_precip_thresh  # rezero NI if it rains
        this.fire_weather_index = 0.0
    else
        # Calculate dewpoint temperature
        t_dew = dewpoint(temp_C, rh)

        # Accumulate Nesterov index over fire season.
        this.fire_weather_index += calc_nesterov_index(temp_C, t_dew)
    end
    return this
end

"""
    calc_nesterov_index(temp_C, t_dew) -> Float64

Current day's Nesterov index for the given daily averaged temperature `temp_C`
[deg C] and dewpoint temperature `t_dew` [deg C]. Nesterov 1968, Eq. 5 of
Thonicke et al. 2010. Clamped at zero (cannot be negative).
"""
function calc_nesterov_index(temp_C::Float64, t_dew::Float64)
    ni = (temp_C - t_dew) * temp_C
    if ni < 0.0  # can't be negative
        ni = 0.0
    end
    return ni
end

"""
    dewpoint(temp_C, rh) -> Float64

Dewpoint temperature [deg C] from air temperature `temp_C` [deg C] and relative
humidity `rh` [%]. Equation 8 from Lawrence 2005
(https://doi.org/10.1175/BAMS-86-2-225).
"""
function dewpoint(temp_C::Float64, rh::Float64)
    yipsolon = log(max(1.0, rh) / 100.0) + (dewpoint_a * temp_C) / (dewpoint_b + temp_C)
    return (dewpoint_b * yipsolon) / (dewpoint_a - yipsolon)
end
