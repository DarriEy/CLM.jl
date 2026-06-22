# FatesFuelMod.jl
# Julia port of FATES src/fates/fire/FatesFuelMod.F90
#
# The SPITFIRE `fuel_type`: holds per-fuel-class loading, accumulates litter into
# the fuel pools, computes total/non-trunk loadings, fuel moisture, average bulk
# density, surface-area-to-volume ratio, and fuel-moisture-of-extinction, plus
# fractional loadings. Only the non-trunk (1-h/10-h/100-h) fuel classes influence
# fire spread (Rothermel 1972; Wilson 1982; Pyne et al. 1996).
#
# Julia port: Fortran derived type `fuel_type` -> mutable struct `fuel_type`,
# with per-fuel-class arrays (sized `num_fuel_classes` from FatesFuelClassesMod)
# as `Vector{Float64}` (SoA). Type-bound procedures become bang-functions
# dispatching on `fuel_type`. The `fire_weather` argument uses the merged
# abstract type (SFFireWeatherMod); the Nesterov branch dispatches on
# `nesterov_index` (SFNesterovMod). `fates_r8` -> Float64.
#
# Deps: FatesFuelClassesMod (num_fuel_classes, fuel_classes + accessors),
# FatesConstantsMod (nearzero), SFFireWeatherMod (fire_weather), SFNesterovMod
# (nesterov_index), FatesGlobals (fates_log, fates_endrun). Standalone — NOT
# added to CLMInstances or any dual-copied struct.

"""
    fuel_type

SPITFIRE fuel state for a patch (Fortran `fuel_type`). Holds per-fuel-class
loading and derived quantities (fractional loading, fraction burnt, effective
moisture) plus the non-trunk aggregate diagnostics used by the fire-spread
calculation.

Fields (per-class arrays are length [`num_fuel_classes`](@ref)):
- `loading` — fuel loading of each fuel class [kgC/m2]
- `effective_moisture` — fuel effective moisture (moisture/MEF) of each class [m3/m3]
- `frac_loading` — fractional loading of each fuel class [0-1]
- `frac_burnt` — fraction of litter burnt by fire of each class [0-1]
- `non_trunk_loading` — total fuel loading excluding trunks [kgC/m2]
- `average_moisture_notrunks` — weighted-average fuel moisture across non-trunk classes [m3/m3]
- `bulk_density_notrunks` — weighted-average bulk density across non-trunk classes [kg/m3]
- `SAV_notrunks` — weighted-average surface-area-to-volume ratio across non-trunk classes [/cm]
- `MEF_notrunks` — weighted-average moisture of extinction across non-trunk classes [m3/m3]
"""
Base.@kwdef mutable struct fuel_type
    loading::Vector{Float64}            = zeros(Float64, num_fuel_classes)  # [kgC/m2]
    effective_moisture::Vector{Float64} = zeros(Float64, num_fuel_classes)  # [m3/m3]
    frac_loading::Vector{Float64}       = zeros(Float64, num_fuel_classes)  # [0-1]
    frac_burnt::Vector{Float64}         = zeros(Float64, num_fuel_classes)  # [0-1]
    non_trunk_loading::Float64          = 0.0  # total fuel loading excluding trunks [kgC/m2]
    average_moisture_notrunks::Float64  = 0.0  # weighted avg fuel moisture, non-trunks [m3/m3]
    bulk_density_notrunks::Float64      = 0.0  # weighted avg bulk density, non-trunks [kg/m3]
    SAV_notrunks::Float64               = 0.0  # weighted avg SAV, non-trunks [/cm]
    MEF_notrunks::Float64               = 0.0  # weighted avg moisture of extinction, non-trunks [m3/m3]
end

"""
    init_fuel!(this::fuel_type)

Initialize the fuel class (Fortran `Init`): just zero everything.
"""
function init_fuel!(this::fuel_type)
    this.loading[1:num_fuel_classes]            .= 0.0
    this.frac_loading[1:num_fuel_classes]       .= 0.0
    this.frac_burnt[1:num_fuel_classes]         .= 0.0
    this.effective_moisture[1:num_fuel_classes] .= 0.0
    this.non_trunk_loading         = 0.0
    this.average_moisture_notrunks = 0.0
    this.bulk_density_notrunks     = 0.0
    this.SAV_notrunks              = 0.0
    this.MEF_notrunks              = 0.0
    return this
end

"""
    fuse!(this::fuel_type, self_area, donor_area, donor_fuel::fuel_type)

Fuse the attributes of `this` fuel object with a `donor_fuel` (Fortran `Fuse`),
area-weighting by `self_area` and `donor_area` [m2].
"""
function fuse!(this::fuel_type, self_area::Float64, donor_area::Float64,
               donor_fuel::fuel_type)
    self_weight  = self_area / (donor_area + self_area)
    donor_weight = 1.0 - self_weight

    for i in 1:num_fuel_classes
        this.loading[i] = this.loading[i] * self_weight +
            donor_fuel.loading[i] * donor_weight
        this.frac_loading[i] = this.frac_loading[i] * self_weight +
            donor_fuel.frac_loading[i] * donor_weight
        this.frac_burnt[i] = this.frac_burnt[i] * self_weight +
            donor_fuel.frac_burnt[i] * donor_weight
        this.effective_moisture[i] = this.effective_moisture[i] * self_weight +
            donor_fuel.effective_moisture[i] * donor_weight
    end

    this.non_trunk_loading = this.non_trunk_loading * self_weight +
        donor_fuel.non_trunk_loading * donor_weight
    this.average_moisture_notrunks = this.average_moisture_notrunks * self_weight +
        donor_fuel.average_moisture_notrunks * donor_weight
    this.bulk_density_notrunks = this.bulk_density_notrunks * self_weight +
        donor_fuel.bulk_density_notrunks * donor_weight
    this.SAV_notrunks = this.SAV_notrunks * self_weight +
        donor_fuel.SAV_notrunks * donor_weight
    this.MEF_notrunks = this.MEF_notrunks * self_weight +
        donor_fuel.MEF_notrunks * donor_weight

    return this
end

"""
    update_loading!(this::fuel_type, leaf_litter, twig_litter,
                    small_branch_litter, large_branch_litter, trunk_litter,
                    live_grass)

Update the loading for each fuel class from input litter pools (Fortran
`UpdateLoading`). All inputs are [kgC/m2].
"""
function update_loading!(this::fuel_type, leaf_litter::Float64, twig_litter::Float64,
                         small_branch_litter::Float64, large_branch_litter::Float64,
                         trunk_litter::Float64, live_grass_litter::Float64)
    # NB: `live_grass_litter` rather than the Fortran `live_grass` — the latter
    # collides with the `live_grass(fuel_classes)` accessor function in Julia.
    this.loading[dead_leaves(fuel_classes)]    = leaf_litter
    this.loading[twigs(fuel_classes)]          = twig_litter
    this.loading[small_branches(fuel_classes)] = small_branch_litter
    this.loading[large_branches(fuel_classes)] = large_branch_litter
    this.loading[live_grass(fuel_classes)]     = live_grass_litter
    this.loading[trunks(fuel_classes)]         = trunk_litter
    return this
end

"""
    sum_loading!(this::fuel_type)

Sum up the loading, excluding trunks (Fortran `SumLoading`). Only the 1-h, 10-h
and 100-h fuel classes influence fire spread.
"""
function sum_loading!(this::fuel_type)
    this.non_trunk_loading = 0.0
    for i in 1:num_fuel_classes
        if i != trunks(fuel_classes)
            this.non_trunk_loading += this.loading[i]
        end
    end
    return this
end

"""
    calculate_fractional_loading!(this::fuel_type)

Calculate the fractional loading for the fuel (Fortran
`CalculateFractionalLoading`). Trunks get a fractional loading of zero; if there
is no non-trunk loading, all fractional loadings are zeroed.
"""
function calculate_fractional_loading!(this::fuel_type)
    # sum up loading just in case
    sum_loading!(this)

    if this.non_trunk_loading > nearzero
        for i in 1:num_fuel_classes
            if i != trunks(fuel_classes)
                this.frac_loading[i] = this.loading[i] / this.non_trunk_loading
            else
                this.frac_loading[i] = 0.0
            end
        end
    else
        this.frac_loading[1:num_fuel_classes] .= 0.0
        this.non_trunk_loading = 0.0
    end
    return this
end

"""
    update_fuel_moisture!(this::fuel_type, sav_fuel, drying_ratio,
                          fireWeatherClass::fire_weather)

Update the fuel moisture depending on which fire-weather class is in use (Fortran
`UpdateFuelMoisture`). `sav_fuel` is the per-class surface-area-to-volume ratio
[/cm], `drying_ratio` the drying ratio, and `fireWeatherClass` a
[`fire_weather`](@ref) subtype. Computes per-class effective moisture and the
non-trunk weighted-average moisture + moisture of extinction.
"""
function update_fuel_moisture!(this::fuel_type, sav_fuel::AbstractVector{Float64},
                               drying_ratio::Float64, fireWeatherClass::fire_weather)
    moisture               = zeros(Float64, num_fuel_classes)  # fuel moisture [m3/m3]
    moisture_of_extinction = zeros(Float64, num_fuel_classes)  # MEF [m3/m3]

    if this.non_trunk_loading + this.loading[trunks(fuel_classes)] > nearzero
        # calculate fuel moisture [m3/m3] for each fuel class depending on which
        # fire weather class is in use
        if fireWeatherClass isa nesterov_index
            calculate_fuel_moisture_nesterov(sav_fuel, drying_ratio,
                fireWeatherClass.fire_weather_index, moisture)
        else
            println(_fates_log_io(), "Unknown fire weather class selected.")
            println(_fates_log_io(), "Choose a different fire weather class or update this subroutine.")
            fates_endrun("FatesFuelMod: unknown fire weather class")
        end

        this.average_moisture_notrunks = 0.0
        this.MEF_notrunks = 0.0
        for i in 1:num_fuel_classes
            # calculate moisture of extinction and fuel effective moisture
            moisture_of_extinction[i] = moisture_of_extinction_fn(sav_fuel[i])
            this.effective_moisture[i] = moisture[i] / moisture_of_extinction[i]

            # average fuel moisture and MEF across all fuel types except trunks [m3/m3]
            if i != trunks(fuel_classes)
                this.average_moisture_notrunks += this.frac_loading[i] * moisture[i]
                this.MEF_notrunks += this.frac_loading[i] * moisture_of_extinction[i]
            end
        end
    else
        this.effective_moisture[1:num_fuel_classes] .= 0.0
        this.average_moisture_notrunks = 0.0
        this.MEF_notrunks = 0.0
    end
    return this
end

"""
    calculate_fuel_moisture_nesterov(sav_fuel, drying_ratio, NI, moisture)

Update the fuel `moisture` (in place, length [`num_fuel_classes`](@ref)) for the
Nesterov index (Fortran `CalculateFuelMoistureNesterov`). `sav_fuel` is the
per-class surface-area-to-volume ratio [/cm], `drying_ratio` the drying ratio,
and `NI` the Nesterov index. Live grass uses the twig (1-hour) SAV.
"""
function calculate_fuel_moisture_nesterov(sav_fuel::AbstractVector{Float64},
                                          drying_ratio::Float64, NI::Float64,
                                          moisture::AbstractVector{Float64})
    for i in 1:num_fuel_classes
        if i == live_grass(fuel_classes)
            # live grass moisture is a function of SAV and changes via Nesterov
            # Index along the same relationship as the 1 hour fuels; live grass
            # has same SAV as dead grass but retains more moisture here
            alpha_FMC = sav_fuel[twigs(fuel_classes)] / drying_ratio
        else
            alpha_FMC = sav_fuel[i] / drying_ratio
        end
        moisture[i] = exp(-1.0 * alpha_FMC * NI)
    end
    return moisture
end

"""
    moisture_of_extinction_fn(sav) -> Float64

Moisture of extinction [m3/m3] from the input surface-area-to-volume ratio `sav`
[/cm] (Fortran `MoistureOfExtinction`). Eq. 27 of Peterson and Ryan (1986).
Errors if `sav <= 0`.
"""
function moisture_of_extinction_fn(sav::Float64)
    MEF_a = 0.524
    MEF_b = 0.066

    if sav <= 0.0
        println(_fates_log_io(), "SAV cannot be negative - SAV")
        fates_endrun("FatesFuelMod: SAV cannot be negative")
        return 0.0
    else
        return MEF_a - MEF_b * log(sav)
    end
end

"""
    average_bulk_density_notrunks!(this::fuel_type, bulk_density)

Calculate the average bulk density excluding trunks (Fortran
`AverageBulkDensity_NoTrunks`). `bulk_density` is the per-class bulk density
[kg/m3]. If there is non-trunk loading it is the frac-loading-weighted average
over non-trunk classes; otherwise the plain mean over all classes.
"""
function average_bulk_density_notrunks!(this::fuel_type,
                                        bulk_density::AbstractVector{Float64})
    if this.non_trunk_loading > nearzero
        this.bulk_density_notrunks = 0.0
        for i in 1:num_fuel_classes
            # average bulk density across all fuel types except trunks
            if i != trunks(fuel_classes)
                this.bulk_density_notrunks += this.frac_loading[i] * bulk_density[i]
            end
        end
    else
        this.bulk_density_notrunks = sum(@view bulk_density[1:num_fuel_classes]) / num_fuel_classes
    end
    return this
end

"""
    average_sav_notrunks!(this::fuel_type, sav_fuel)

Calculate the average surface-area-to-volume ratio excluding trunks (Fortran
`AverageSAV_NoTrunks`). `sav_fuel` is the per-class SAV [/cm]. If there is
non-trunk loading it is the frac-loading-weighted average over non-trunk
classes; otherwise the plain mean over all classes.
"""
function average_sav_notrunks!(this::fuel_type, sav_fuel::AbstractVector{Float64})
    if this.non_trunk_loading > nearzero
        this.SAV_notrunks = 0.0
        for i in 1:num_fuel_classes
            # average SAV across all fuel types except trunks
            if i != trunks(fuel_classes)
                this.SAV_notrunks += this.frac_loading[i] * sav_fuel[i]
            end
        end
    else
        this.SAV_notrunks = sum(@view sav_fuel[1:num_fuel_classes]) / num_fuel_classes
    end
    return this
end

# Helper: resolve a writable IO from the FATES log unit (6 -> stdout).
_fates_log_io() = (fates_log() == 6 ? stdout : stderr)
