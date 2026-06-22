# SFFireWeatherMod.jl
# Julia port of FATES src/fates/fire/SFFireWeatherMod.F90
#
# SPITFIRE fire-weather base state. The Fortran defines an *abstract* type
# `fire_weather` (holding the fire-weather index + effective wind speed) with two
# deferred type-bound procedures (Init, UpdateIndex) that concrete subtypes must
# implement, plus a concrete base method UpdateEffectiveWindSpeed. The concrete
# fire-weather class (e.g. the Nesterov index) lives in a sibling module ported
# later.
#
# Fortran abstract type -> Julia `abstract type fire_weather`.
# Fields of the abstract type -> live in concrete subtypes (Julia abstract types
#   hold no fields), so we provide accessors and require subtypes to carry
#   `fire_weather_index` and `effective_windspeed` fields.
# Deferred procedures -> generic functions `init_fire_weather!` / `update_index!`
#   that error unless a subtype overloads them.
# Concrete base method -> `update_effective_windspeed!` (works on any subtype).
# `fates_r8` -> Float64. Deps: FatesConstantsMod (fates_r8).

"""
    fire_weather

Abstract base type for SPITFIRE fire-weather state. Concrete subtypes must
provide the fields `fire_weather_index::Float64` (fire weather index) and
`effective_windspeed::Float64` (effective wind speed corrected for tree/grass
cover [m/min]), and must implement [`init_fire_weather!`](@ref) and
[`update_index!`](@ref).
"""
abstract type fire_weather end

# Wind attenuation factors (Fortran `parameter`s in UpdateEffectiveWindSpeed).
const wind_atten_treed = 0.4   # wind attenuation factor for tree fraction
const wind_atten_grass = 0.6   # wind attenuation factor for grass fraction

"""
    init_fire_weather!(this::fire_weather)

Deferred initializer (Fortran `Init`). Concrete subtypes must overload this; the
abstract fallback errors.
"""
function init_fire_weather!(this::fire_weather)
    error("init_fire_weather! not implemented for $(typeof(this))")
end

"""
    update_index!(this::fire_weather, temp_C, precip, rh, wind)

Deferred fire-weather-index update (Fortran `UpdateIndex`). Concrete subtypes
(e.g. the Nesterov index) must overload this; the abstract fallback errors.
"""
function update_index!(this::fire_weather, temp_C::Float64, precip::Float64,
                       rh::Float64, wind::Float64)
    error("update_index! not implemented for $(typeof(this))")
end

"""
    update_effective_windspeed!(this::fire_weather, wind_speed, tree_fraction,
                                grass_fraction, bare_fraction)

Calculate the effective wind speed, corrected for tree/grass/bare cover.
Concrete base method (Fortran `UpdateEffectiveWindSpeed`), works on any subtype.

`wind_speed`, `tree_fraction`, `grass_fraction`, `bare_fraction` are the wind
speed [m/min] and the [0-1] cover fractions. Sets `this.effective_windspeed`.
"""
function update_effective_windspeed!(this::fire_weather, wind_speed::Float64,
                                     tree_fraction::Float64, grass_fraction::Float64,
                                     bare_fraction::Float64)
    this.effective_windspeed = wind_speed * (tree_fraction * wind_atten_treed +
        (grass_fraction + bare_fraction) * wind_atten_grass)
    return this
end

"""
    base_fire_weather

A minimal concrete subtype of [`fire_weather`](@ref) carrying just the base
state. Useful for exercising the base accessors / `update_effective_windspeed!`
before the index-specific subtypes (Nesterov, etc.) are ported.
"""
Base.@kwdef mutable struct base_fire_weather <: fire_weather
    fire_weather_index::Float64  = 0.0  # fire weather index
    effective_windspeed::Float64 = 0.0  # effective wind speed [m/min]
end

"""
    init_fire_weather!(this::base_fire_weather)

Initialize the base fire-weather state (zero index and effective wind speed).
"""
function init_fire_weather!(this::base_fire_weather)
    this.fire_weather_index  = 0.0
    this.effective_windspeed = 0.0
    return this
end
