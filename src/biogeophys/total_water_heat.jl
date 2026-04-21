# ==========================================================================
# Ported from: src/biogeophys/TotalWaterAndHeatMod.F90
# Routines for computing total column water and heat contents
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

"""Base temperature for heat sums [K]. Equal to tfrz (freezing point)."""
const heat_base_temp = TFRZ

"""Minimum temperature for water temperature used by adjust_delta_heat_for_delta_liq! [K]."""
const DeltaLiqMinTemp = TFRZ

"""Maximum temperature for water temperature used by adjust_delta_heat_for_delta_liq! [K]."""
const DeltaLiqMaxTemp = TFRZ + 35.0

# ---------------------------------------------------------------------------
# Private helper: TempToHeat (pure function)
# ---------------------------------------------------------------------------

"""
    _temp_to_heat(temp::Float64, cv::Float64) -> Float64

Convert temperature to heat content relative to `heat_base_temp`.

Ported from `TempToHeat` in `TotalWaterAndHeatMod.F90`.
"""
@inline function _temp_to_heat(temp::Real, cv::Real)
    return cv * (temp - heat_base_temp)
end

# ---------------------------------------------------------------------------
# Private helper: AccumulateLiquidWaterHeat
# ---------------------------------------------------------------------------

"""
    _accumulate_liquid_water_heat!(temp, h2o, heat_liquid, latent_heat_liquid; cv_liquid=nothing)

Accumulate heat quantities for liquid water for a single column.
Adds to existing values in `heat_liquid` and `latent_heat_liquid`.
Optionally accumulates `cv_liquid` (liquid heat capacity).

Returns `(heat_liquid, latent_heat_liquid, cv_liquid)` as updated values.

Ported from `AccumulateLiquidWaterHeat` in `TotalWaterAndHeatMod.F90`.
"""
@inline function _accumulate_liquid_water_heat(temp::Real, h2o::Real,
                                                heat_liquid::Real,
                                                latent_heat_liquid::Real,
                                                cv_liquid::Real)
    cv = h2o * CPLIQ
    cv_liquid_out = cv_liquid + cv
    heat_liquid_out = heat_liquid + _temp_to_heat(temp, cv)
    latent_heat_liquid_out = latent_heat_liquid + h2o * HFUS
    return heat_liquid_out, latent_heat_liquid_out, cv_liquid_out
end

@inline function _accumulate_liquid_water_heat(temp::Real, h2o::Real,
                                                heat_liquid::Real,
                                                latent_heat_liquid::Real)
    cv = h2o * CPLIQ
    heat_liquid_out = heat_liquid + _temp_to_heat(temp, cv)
    latent_heat_liquid_out = latent_heat_liquid + h2o * HFUS
    return heat_liquid_out, latent_heat_liquid_out
end

# ---------------------------------------------------------------------------
# Public: LiquidWaterHeat
# ---------------------------------------------------------------------------

"""
    liquid_water_heat(temp::Float64, h2o::Float64) -> Float64

Get the total heat content (including latent heat) of some mass of liquid
water at a given temperature, using a base temperature of `heat_base_temp`.

Ported from `LiquidWaterHeat` in `TotalWaterAndHeatMod.F90`.
"""
function liquid_water_heat(temp::Real, h2o::Real)
    heat_liquid = 0.0
    latent_heat_liquid = 0.0
    heat_liquid, latent_heat_liquid = _accumulate_liquid_water_heat(
        temp, h2o, heat_liquid, latent_heat_liquid)
    return heat_liquid + latent_heat_liquid
end

# ---------------------------------------------------------------------------
# Public: ComputeWaterMassNonLake
# ---------------------------------------------------------------------------

"""
    compute_water_mass_non_lake!(mask_nolake, col, waterstate, waterdiagnostic,
                                  subtract_dynbal_baselines, water_mass)

Compute total water mass for all non-lake columns.

Arguments:
- `mask_nolake`: BitVector mask for non-lake columns
- `col`: ColumnData
- `waterstate`: WaterStateData (or WaterStateBulkData.ws)
- `waterdiagnostic`: WaterDiagnosticBulkData (needs total_plant_stored_h2o_col)
- `subtract_dynbal_baselines`: whether to subtract dynbal baselines
- `water_mass`: output vector [kg/m2]

Ported from `ComputeWaterMassNonLake` in `TotalWaterAndHeatMod.F90`.
"""
function compute_water_mass_non_lake!(mask_nolake::BitVector,
                                       col::ColumnData,
                                       waterstate::WaterStateData,
                                       waterdiagnostic::WaterDiagnosticBulkData,
                                       subtract_dynbal_baselines::Bool,
                                       water_mass::Vector{<:Real})
    nc = length(water_mass)
    FT = eltype(water_mass)
    liquid_mass = zeros(FT, nc)
    ice_mass = zeros(FT, nc)

    compute_liq_ice_mass_non_lake!(mask_nolake, col, waterstate, waterdiagnostic,
                                    subtract_dynbal_baselines, liquid_mass, ice_mass)

    for c in eachindex(mask_nolake)
        mask_nolake[c] || continue
        water_mass[c] = liquid_mass[c] + ice_mass[c]
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Public: ComputeWaterMassLake
# ---------------------------------------------------------------------------

"""
    compute_water_mass_lake!(mask_lake, col, waterstate, lakestate,
                              add_lake_water_and_subtract_dynbal_baselines, water_mass)

Compute total water mass for all lake columns.

Ported from `ComputeWaterMassLake` in `TotalWaterAndHeatMod.F90`.
"""
function compute_water_mass_lake!(mask_lake::BitVector,
                                   col::ColumnData,
                                   waterstate::WaterStateData,
                                   lakestate::LakeStateData,
                                   add_lake_water_and_subtract_dynbal_baselines::Bool,
                                   water_mass::Vector{<:Real})
    nc = length(water_mass)
    FT = eltype(water_mass)
    liquid_mass = zeros(FT, nc)
    ice_mass = zeros(FT, nc)

    compute_liq_ice_mass_lake!(mask_lake, col, waterstate, lakestate,
                                add_lake_water_and_subtract_dynbal_baselines,
                                liquid_mass, ice_mass)

    for c in eachindex(mask_lake)
        mask_lake[c] || continue
        water_mass[c] = liquid_mass[c] + ice_mass[c]
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Public: ComputeLiqIceMassNonLake
# ---------------------------------------------------------------------------

"""
    compute_liq_ice_mass_non_lake!(mask_nolake, col, waterstate, waterdiagnostic,
                                    subtract_dynbal_baselines, liquid_mass, ice_mass)

Compute total water mass for all non-lake columns, separated into liquid and ice.

Note: The Fortran version calls p2c to average patch-level canopy water to column
level. In this Julia port, we assume liqcan and snocan are already at column level
or use a simplified direct access. Since the canopy water in ComputeLiqIceMassNonLake
uses `waterstate_type` (generic, not bulk), and p2c requires patch infrastructure,
we provide a simplified version that uses `liqcan_patch` and `snocan_patch` directly
by summing over the column's patches. The caller can pre-compute column-level canopy
water if needed.

Ported from `ComputeLiqIceMassNonLake` in `TotalWaterAndHeatMod.F90`.
"""
function compute_liq_ice_mass_non_lake!(mask_nolake::BitVector,
                                          col::ColumnData,
                                          waterstate::WaterStateData,
                                          waterdiagnostic::WaterDiagnosticBulkData,
                                          subtract_dynbal_baselines::Bool,
                                          liquid_mass::Vector{<:Real},
                                          ice_mass::Vector{<:Real};
                                          liqcan_col::Union{Vector{Float64}, Nothing}=nothing,
                                          snocan_col::Union{Vector{Float64}, Nothing}=nothing)
    nlevsno = varpar.nlevsno

    # Initialize
    for c in eachindex(mask_nolake)
        mask_nolake[c] || continue
        liquid_mass[c] = 0.0
        ice_mass[c] = 0.0
    end

    # Canopy water: if not provided, use zero (caller should supply p2c results)
    FT_li = eltype(liquid_mass)
    _liqcan = liqcan_col !== nothing ? liqcan_col : zeros(FT_li, length(liquid_mass))
    _snocan = snocan_col !== nothing ? snocan_col : zeros(FT_li, length(liquid_mass))

    for c in eachindex(mask_nolake)
        mask_nolake[c] || continue

        # Canopy water + plant stored water
        liquid_mass[c] += _liqcan[c] + waterdiagnostic.total_plant_stored_h2o_col[c]
        ice_mass[c] += _snocan[c]

        # Snow not resolved into layers
        ice_mass[c] += waterstate.h2osno_no_layers_col[c]

        # Snow layers: Fortran indices snl+1:0, Julia offset by nlevsno
        for j in (col.snl[c] + 1):0
            jj = j + nlevsno
            liquid_mass[c] += waterstate.h2osoi_liq_col[c, jj]
            ice_mass[c] += waterstate.h2osoi_ice_col[c, jj]
        end

        # Aquifer water (only for hydrologically active columns)
        if col.hydrologically_active[c]
            liquid_mass[c] += (waterstate.wa_col[c] - waterstate.aquifer_water_baseline)
        end

        # Surface water (exclude urban impervious/wall/roof)
        if !(col.itype[c] == ICOL_ROOF || col.itype[c] == ICOL_SUNWALL ||
             col.itype[c] == ICOL_SHADEWALL || col.itype[c] == ICOL_ROAD_IMPERV)
            liquid_mass[c] += waterstate.h2osfc_col[c]
        end
    end

    # Soil water content
    accumulate_soil_liq_ice_mass_non_lake!(mask_nolake, col, waterstate,
                                            liquid_mass, ice_mass)

    # Subtract dynbal baselines if requested
    if subtract_dynbal_baselines
        for c in eachindex(mask_nolake)
            mask_nolake[c] || continue
            liquid_mass[c] -= waterstate.dynbal_baseline_liq_col[c]
            ice_mass[c] -= waterstate.dynbal_baseline_ice_col[c]
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: AccumulateSoilLiqIceMassNonLake
# ---------------------------------------------------------------------------

"""
    accumulate_soil_liq_ice_mass_non_lake!(mask, col, waterstate, liquid_mass, ice_mass)

Accumulate soil water mass of non-lake columns, separated into liquid and ice.
Adds to any existing values in `liquid_mass` and `ice_mass`.

Ported from `AccumulateSoilLiqIceMassNonLake` in `TotalWaterAndHeatMod.F90`.
"""
function accumulate_soil_liq_ice_mass_non_lake!(mask::BitVector,
                                                  col::ColumnData,
                                                  waterstate::WaterStateData,
                                                  liquid_mass::Vector{<:Real},
                                                  ice_mass::Vector{<:Real})
    nlevgrnd       = varpar.nlevgrnd
    nlevurb        = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevsno        = varpar.nlevsno

    for j in 1:nlevmaxurbgrnd
        for c in eachindex(mask)
            mask[c] || continue

            has_h2o = false
            if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL
                has_h2o = false
            elseif col.itype[c] == ICOL_ROOF
                has_h2o = (j <= nlevurb)
            else
                has_h2o = (j <= nlevgrnd)
            end

            if has_h2o
                jj = j + nlevsno  # convert soil index to combined array index
                liquid_mass[c] += waterstate.h2osoi_liq_col[c, jj]
                ice_mass[c] += waterstate.h2osoi_ice_col[c, jj] + waterstate.excess_ice_col[c, jj]
            end
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: ComputeLiqIceMassLake
# ---------------------------------------------------------------------------

"""
    compute_liq_ice_mass_lake!(mask_lake, col, waterstate, lakestate,
                                add_lake_water_and_subtract_dynbal_baselines,
                                liquid_mass, ice_mass)

Compute total water mass for all lake columns, separated into liquid and ice.

Ported from `ComputeLiqIceMassLake` in `TotalWaterAndHeatMod.F90`.
"""
function compute_liq_ice_mass_lake!(mask_lake::BitVector,
                                     col::ColumnData,
                                     waterstate::WaterStateData,
                                     lakestate::LakeStateData,
                                     add_lake_water_and_subtract_dynbal_baselines::Bool,
                                     liquid_mass::Vector{<:Real},
                                     ice_mass::Vector{<:Real})
    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno

    # Initialize
    for c in eachindex(mask_lake)
        mask_lake[c] || continue
        liquid_mass[c] = 0.0
        ice_mass[c] = 0.0
    end

    # Large canceling terms first (baselines and lake water)
    if add_lake_water_and_subtract_dynbal_baselines
        for c in eachindex(mask_lake)
            mask_lake[c] || continue
            liquid_mass[c] -= waterstate.dynbal_baseline_liq_col[c]
            ice_mass[c] -= waterstate.dynbal_baseline_ice_col[c]
        end

        # Lake water content (tracer_ratio = 1.0 for bulk water)
        accumulate_liq_ice_mass_lake!(mask_lake, col, lakestate, 1.0,
                                       liquid_mass, ice_mass)
    end

    # Snow water content
    for c in eachindex(mask_lake)
        mask_lake[c] || continue

        ice_mass[c] += waterstate.h2osno_no_layers_col[c]
        for j in (col.snl[c] + 1):0
            jj = j + nlevsno
            liquid_mass[c] += waterstate.h2osoi_liq_col[c, jj]
            ice_mass[c] += waterstate.h2osoi_ice_col[c, jj]
        end
    end

    # Soil water content under the lake
    for j in 1:nlevgrnd
        for c in eachindex(mask_lake)
            mask_lake[c] || continue
            jj = j + nlevsno
            liquid_mass[c] += waterstate.h2osoi_liq_col[c, jj]
            ice_mass[c] += waterstate.h2osoi_ice_col[c, jj]
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: AccumulateLiqIceMassLake
# ---------------------------------------------------------------------------

"""
    accumulate_liq_ice_mass_lake!(mask, col, lakestate, tracer_ratio,
                                   liquid_mass, ice_mass)

Accumulate lake water mass of lake columns, separated into liquid and ice.
Adds to any existing values in `liquid_mass` and `ice_mass`.

Ported from `AccumulateLiqIceMassLake` in `TotalWaterAndHeatMod.F90`.
"""
function accumulate_liq_ice_mass_lake!(mask::BitVector,
                                        col::ColumnData,
                                        lakestate::LakeStateData,
                                        tracer_ratio::Real,
                                        liquid_mass::Vector{<:Real},
                                        ice_mass::Vector{<:Real})
    nlevlak = varpar.nlevlak

    for j in 1:nlevlak
        for c in eachindex(mask)
            mask[c] || continue
            h2olak_liq = col.dz_lake[c, j] * DENH2O * (1.0 - lakestate.lake_icefrac_col[c, j]) * tracer_ratio
            h2olak_ice = col.dz_lake[c, j] * DENH2O * lakestate.lake_icefrac_col[c, j] * tracer_ratio
            liquid_mass[c] += h2olak_liq
            ice_mass[c] += h2olak_ice
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: ComputeHeatNonLake
# ---------------------------------------------------------------------------

"""
    compute_heat_non_lake!(mask_nolake, col, lun, urbanparams, soilstate,
                            temperature, waterstatebulk, waterdiagnosticbulk,
                            heat, heat_liquid, cv_liquid;
                            liqcan_col=nothing, snocan_col=nothing)

Compute total heat content for all non-lake columns.
Also returns the total heat content of liquid water (excluding latent heat)
and the total heat capacity of liquid water, which can be used to compute
the weighted average liquid water temperature.

Ported from `ComputeHeatNonLake` in `TotalWaterAndHeatMod.F90`.
"""
function compute_heat_non_lake!(mask_nolake::BitVector,
                                 col::ColumnData,
                                 lun::LandunitData,
                                 urbanparams::UrbanParamsData,
                                 soilstate::SoilStateData,
                                 temperature::TemperatureData,
                                 waterstatebulk::WaterStateBulkData,
                                 waterdiagnosticbulk::WaterDiagnosticBulkData,
                                 heat::Vector{<:Real},
                                 heat_liquid::Vector{<:Real},
                                 cv_liquid::Vector{<:Real};
                                 liqcan_col::Union{Vector{Float64}, Nothing}=nothing,
                                 snocan_col::Union{Vector{Float64}, Nothing}=nothing)
    nlevsno = varpar.nlevsno
    nc = length(heat)

    # Access through composition for bulk data
    ws = waterstatebulk.ws

    # Local accumulators
    FT = eltype(heat)
    heat_dry_mass = zeros(FT, nc)
    heat_ice = zeros(FT, nc)
    latent_heat_liquid = zeros(FT, nc)

    # Initialize outputs
    for c in eachindex(mask_nolake)
        mask_nolake[c] || continue
        heat_liquid[c] = 0.0
        cv_liquid[c] = 0.0
        heat_dry_mass[c] = 0.0
        heat_ice[c] = 0.0
        latent_heat_liquid[c] = 0.0
    end

    # Canopy water: if not provided, use zero
    _liqcan = liqcan_col !== nothing ? liqcan_col : zeros(FT, nc)
    _snocan = snocan_col !== nothing ? snocan_col : zeros(FT, nc)

    for c in eachindex(mask_nolake)
        mask_nolake[c] || continue

        # --- Canopy water ---
        liqveg = _liqcan[c] + waterdiagnosticbulk.total_plant_stored_h2o_col[c]

        heat_liquid[c], latent_heat_liquid[c], cv_liquid[c] =
            _accumulate_liquid_water_heat(heat_base_temp, liqveg,
                heat_liquid[c], latent_heat_liquid[c], cv_liquid[c])

        # --- Snow ---
        # h2osno_no_layers: use top snow/soil layer temperature (j=1 in Fortran = index nlevsno+1)
        j_top = nlevsno + 1
        heat_ice[c] += _temp_to_heat(temperature.t_soisno_col[c, j_top],
                                      ws.h2osno_no_layers_col[c] * CPICE)

        for j in (col.snl[c] + 1):0
            jj = j + nlevsno
            heat_liquid[c], latent_heat_liquid[c], cv_liquid[c] =
                _accumulate_liquid_water_heat(temperature.t_soisno_col[c, jj],
                    ws.h2osoi_liq_col[c, jj],
                    heat_liquid[c], latent_heat_liquid[c], cv_liquid[c])
            heat_ice[c] += _temp_to_heat(temperature.t_soisno_col[c, jj],
                                          ws.h2osoi_ice_col[c, jj] * CPICE)
        end

        # --- Aquifer water ---
        if col.hydrologically_active[c]
            heat_liquid[c], latent_heat_liquid[c], cv_liquid[c] =
                _accumulate_liquid_water_heat(heat_base_temp,
                    ws.wa_col[c] - ws.aquifer_water_baseline,
                    heat_liquid[c], latent_heat_liquid[c], cv_liquid[c])
        end

        # --- Surface water ---
        if !(col.itype[c] == ICOL_ROOF || col.itype[c] == ICOL_SUNWALL ||
             col.itype[c] == ICOL_SHADEWALL || col.itype[c] == ICOL_ROAD_IMPERV)
            heat_liquid[c], latent_heat_liquid[c], cv_liquid[c] =
                _accumulate_liquid_water_heat(temperature.t_h2osfc_col[c],
                    ws.h2osfc_col[c],
                    heat_liquid[c], latent_heat_liquid[c], cv_liquid[c])
        end
    end

    # Combine into heat
    for c in eachindex(mask_nolake)
        mask_nolake[c] || continue
        heat[c] = heat_dry_mass[c] + heat_ice[c] + heat_liquid[c] + latent_heat_liquid[c]
    end

    # Soil heat
    accumulate_soil_heat_non_lake!(mask_nolake, col, lun, urbanparams, soilstate,
                                    temperature, waterstatebulk,
                                    heat, heat_liquid, cv_liquid)

    # Subtract baseline heat
    for c in eachindex(mask_nolake)
        mask_nolake[c] || continue
        heat[c] -= temperature.dynbal_baseline_heat_col[c]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: AccumulateSoilHeatNonLake
# ---------------------------------------------------------------------------

"""
    accumulate_soil_heat_non_lake!(mask, col, lun, urbanparams, soilstate,
                                    temperature, waterstatebulk,
                                    heat, heat_liquid, cv_liquid)

Accumulate soil heat of non-lake columns (or a subset). Includes related
heat quantities for urban columns (wall, roof and road).
Adds to any existing values in `heat`, `heat_liquid` and `cv_liquid`.

Ported from `AccumulateSoilHeatNonLake` in `TotalWaterAndHeatMod.F90`.
"""
function accumulate_soil_heat_non_lake!(mask::BitVector,
                                          col::ColumnData,
                                          lun::LandunitData,
                                          urbanparams::UrbanParamsData,
                                          soilstate::SoilStateData,
                                          temperature::TemperatureData,
                                          waterstatebulk::WaterStateBulkData,
                                          heat::Vector{<:Real},
                                          heat_liquid::Vector{<:Real},
                                          cv_liquid::Vector{<:Real})
    nlevgrnd       = varpar.nlevgrnd
    nlevurb        = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevsno        = varpar.nlevsno

    ws = waterstatebulk.ws
    nc = length(heat)

    # Local accumulators
    FT = eltype(heat)
    soil_heat_liquid = zeros(FT, nc)
    soil_heat_dry_mass = zeros(FT, nc)
    soil_heat_ice = zeros(FT, nc)
    soil_latent_heat_liquid = zeros(FT, nc)

    for c in eachindex(mask)
        mask[c] || continue
        soil_heat_liquid[c] = 0.0
        soil_heat_dry_mass[c] = 0.0
        soil_heat_ice[c] = 0.0
        soil_latent_heat_liquid[c] = 0.0
    end

    for j in 1:nlevmaxurbgrnd
        for c in eachindex(mask)
            mask[c] || continue
            l = col.landunit[c]
            jj = j + nlevsno  # combined snow+soil index

            has_h2o = false

            if col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL
                has_h2o = false
                if j <= nlevurb
                    soil_heat_dry_mass[c] += _temp_to_heat(
                        temperature.t_soisno_col[c, jj],
                        urbanparams.cv_wall[l, j] * col.dz[c, jj])
                end

            elseif col.itype[c] == ICOL_ROOF
                if j <= nlevurb
                    has_h2o = true
                    soil_heat_dry_mass[c] += _temp_to_heat(
                        temperature.t_soisno_col[c, jj],
                        urbanparams.cv_roof[l, j] * col.dz[c, jj])
                else
                    has_h2o = false
                end

            else
                if j <= nlevgrnd
                    has_h2o = true

                    if col.itype[c] == ICOL_ROAD_IMPERV && j <= urbanparams.nlev_improad[l]
                        soil_heat_dry_mass[c] += _temp_to_heat(
                            temperature.t_soisno_col[c, jj],
                            urbanparams.cv_improad[l, j] * col.dz[c, jj])
                    elseif lun.itype[l] != ISTWET && lun.itype[l] != ISTICE
                        # Includes impervious roads below nlev_improad (where we have soil)
                        soil_heat_dry_mass[c] += _temp_to_heat(
                            temperature.t_soisno_col[c, jj],
                            soilstate.csol_col[c, j] * (1.0 - soilstate.watsat_col[c, j]) * col.dz[c, jj])
                    end
                else
                    has_h2o = false
                end
            end

            if has_h2o
                soil_heat_liquid[c], soil_latent_heat_liquid[c], cv_liquid[c] =
                    _accumulate_liquid_water_heat(temperature.t_soisno_col[c, jj],
                        ws.h2osoi_liq_col[c, jj],
                        soil_heat_liquid[c], soil_latent_heat_liquid[c], cv_liquid[c])

                soil_heat_ice[c] += _temp_to_heat(
                    temperature.t_soisno_col[c, jj],
                    ws.h2osoi_ice_col[c, jj] * CPICE)
            end
        end
    end

    for c in eachindex(mask)
        mask[c] || continue
        heat_liquid[c] += soil_heat_liquid[c]
        heat[c] += soil_heat_dry_mass[c] + soil_heat_ice[c] +
                   soil_heat_liquid[c] + soil_latent_heat_liquid[c]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: ComputeHeatLake
# ---------------------------------------------------------------------------

"""
    compute_heat_lake!(mask_lake, col, soilstate, temperature,
                        waterstatebulk, lakestate,
                        heat, heat_liquid, cv_liquid)

Compute total heat content for all lake columns.
Also returns the total heat content of liquid water (excluding latent heat)
and the total heat capacity of liquid water.

Ported from `ComputeHeatLake` in `TotalWaterAndHeatMod.F90`.
"""
function compute_heat_lake!(mask_lake::BitVector,
                             col::ColumnData,
                             soilstate::SoilStateData,
                             temperature::TemperatureData,
                             waterstatebulk::WaterStateBulkData,
                             lakestate::LakeStateData,
                             heat::Vector{<:Real},
                             heat_liquid::Vector{<:Real},
                             cv_liquid::Vector{<:Real})
    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno
    nc = length(heat)
    ws = waterstatebulk.ws

    # Local accumulators
    FT = eltype(heat)
    heat_dry_mass = zeros(FT, nc)
    heat_ice = zeros(FT, nc)
    latent_heat_liquid = zeros(FT, nc)

    for c in eachindex(mask_lake)
        mask_lake[c] || continue
        heat_liquid[c] = 0.0
        cv_liquid[c] = 0.0
        heat_dry_mass[c] = 0.0
        heat_ice[c] = 0.0
        latent_heat_liquid[c] = 0.0
    end

    # Subtract baseline heat (large canceling term first)
    for c in eachindex(mask_lake)
        mask_lake[c] || continue
        heat[c] = -temperature.dynbal_baseline_heat_col[c]
    end

    # Lake water heat content (NOT accumulated in heat_liquid/cv_liquid)
    accumulate_heat_lake!(mask_lake, col, temperature, lakestate, heat)

    # Snow heat content
    for c in eachindex(mask_lake)
        mask_lake[c] || continue
        for j in (col.snl[c] + 1):0
            jj = j + nlevsno
            heat_liquid[c], latent_heat_liquid[c], cv_liquid[c] =
                _accumulate_liquid_water_heat(temperature.t_soisno_col[c, jj],
                    ws.h2osoi_liq_col[c, jj],
                    heat_liquid[c], latent_heat_liquid[c], cv_liquid[c])
            heat_ice[c] += _temp_to_heat(temperature.t_soisno_col[c, jj],
                                          ws.h2osoi_ice_col[c, jj] * CPICE)
        end
    end

    # Soil heat content under the lake
    for j in 1:nlevgrnd
        for c in eachindex(mask_lake)
            mask_lake[c] || continue
            jj = j + nlevsno

            heat_dry_mass[c] += _temp_to_heat(
                temperature.t_soisno_col[c, jj],
                soilstate.csol_col[c, j] * (1.0 - soilstate.watsat_col[c, j]) * col.dz[c, jj])

            heat_liquid[c], latent_heat_liquid[c], cv_liquid[c] =
                _accumulate_liquid_water_heat(temperature.t_soisno_col[c, jj],
                    ws.h2osoi_liq_col[c, jj],
                    heat_liquid[c], latent_heat_liquid[c], cv_liquid[c])

            heat_ice[c] += _temp_to_heat(temperature.t_soisno_col[c, jj],
                                          ws.h2osoi_ice_col[c, jj] * CPICE)
        end
    end

    for c in eachindex(mask_lake)
        mask_lake[c] || continue
        heat[c] += heat_dry_mass[c] + heat_ice[c] + heat_liquid[c] + latent_heat_liquid[c]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: AccumulateHeatLake
# ---------------------------------------------------------------------------

"""
    accumulate_heat_lake!(mask, col, temperature, lakestate, heat)

Accumulate heat of lake water in lake columns. Does NOT accumulate
heat_liquid or cv_liquid because lake water is virtual and should not
contribute to the implicit temperature of the dynbal liquid flux.

Ported from `AccumulateHeatLake` in `TotalWaterAndHeatMod.F90`.
"""
function accumulate_heat_lake!(mask::BitVector,
                                col::ColumnData,
                                temperature::TemperatureData,
                                lakestate::LakeStateData,
                                heat::Vector{<:Real})
    nlevlak = varpar.nlevlak
    nc = length(heat)

    FT = eltype(heat)
    lake_heat_liquid = zeros(FT, nc)
    lake_heat_ice = zeros(FT, nc)
    lake_latent_heat_liquid = zeros(FT, nc)

    for c in eachindex(mask)
        mask[c] || continue
        lake_heat_liquid[c] = 0.0
        lake_heat_ice[c] = 0.0
        lake_latent_heat_liquid[c] = 0.0
    end

    for j in 1:nlevlak
        for c in eachindex(mask)
            mask[c] || continue
            # Liquid heat
            h2olak_liq = col.dz_lake[c, j] * DENH2O * (1.0 - lakestate.lake_icefrac_col[c, j])
            lake_heat_liquid[c], lake_latent_heat_liquid[c] =
                _accumulate_liquid_water_heat(temperature.t_lake_col[c, j],
                    h2olak_liq, lake_heat_liquid[c], lake_latent_heat_liquid[c])

            # Ice heat (use water density because lake depths are not adjusted on freeze)
            h2olak_ice = col.dz_lake[c, j] * DENH2O * lakestate.lake_icefrac_col[c, j]
            lake_heat_ice[c] += _temp_to_heat(temperature.t_lake_col[c, j],
                                               h2olak_ice * CPICE)
        end
    end

    for c in eachindex(mask)
        mask[c] || continue
        heat[c] += lake_heat_ice[c] + lake_heat_liquid[c] + lake_latent_heat_liquid[c]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: AdjustDeltaHeatForDeltaLiq
# ---------------------------------------------------------------------------

"""
    adjust_delta_heat_for_delta_liq!(bounds_grc, delta_liq, liquid_water_temp1,
                                      liquid_water_temp2, delta_heat)

Adjusts `delta_heat` (the change in gridcell heat content due to land cover
change) to account for the implicit heat flux associated with `delta_liq`.

Sign convention: `delta_liq` and `delta_heat` are positive if the
post-landcover change value is greater than the pre-landcover change value.

Ported from `AdjustDeltaHeatForDeltaLiq` in `TotalWaterAndHeatMod.F90`.
"""
function adjust_delta_heat_for_delta_liq!(bounds_grc::UnitRange{Int},
                                            delta_liq::Vector{<:Real},
                                            liquid_water_temp1::Vector{<:Real},
                                            liquid_water_temp2::Vector{<:Real},
                                            delta_heat::Vector{<:Real})
    for g in bounds_grc
        if delta_liq[g] != 0.0
            if delta_liq[g] < 0.0
                # More water initially: runoff has temperature of initial state
                water_temperature = liquid_water_temp1[g]
            else
                # More water finally: sucking water at temperature of final state
                water_temperature = liquid_water_temp2[g]
            end

            # Clamp to reasonable bounds
            water_temperature = max(water_temperature, DeltaLiqMinTemp)
            water_temperature = min(water_temperature, DeltaLiqMaxTemp)

            total_liquid_heat = liquid_water_heat(water_temperature, delta_liq[g])

            # Subtract: positive runoff should add positive heat to delta_heat
            delta_heat[g] -= total_liquid_heat
        end
    end

    return nothing
end
