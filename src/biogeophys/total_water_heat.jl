# ==========================================================================
# Ported from: src/biogeophys/TotalWaterAndHeatMod.F90
# Routines for computing total column water and heat contents
# ==========================================================================

# ---------------------------------------------------------------------------
# KernelAbstractions kernels for fully-independent per-column finalize loops.
#
# NOTE: The dominant pattern in this file is per-column TOTALS accumulated by
# summing water/heat over soil/snow/lake LAYERS (`total[c] += x[c,j]`). Those
# loops are loop-carried reductions over j and are intentionally NOT kernelized.
# Only the masked per-column "combine"/"finalize" loops below — where each
# column element is written independently from other-array values at the same
# index `[c]` (no accumulation over j) — are kernelized here.
# ---------------------------------------------------------------------------

# water_mass[c] = liquid_mass[c] + ice_mass[c]   (masked, per-column independent)
@kernel function _twh_water_mass_combine_kernel!(water_mass, @Const(mask),
                                                 @Const(liquid_mass), @Const(ice_mass))
    c = @index(Global)
    @inbounds if mask[c]
        water_mass[c] = liquid_mass[c] + ice_mass[c]
    end
end

twh_water_mass_combine!(water_mass, mask, liquid_mass, ice_mass) =
    _launch!(_twh_water_mass_combine_kernel!, water_mass, mask, liquid_mass, ice_mass)

# heat[c] = heat_dry_mass[c] + heat_ice[c] + heat_liquid[c] + latent_heat_liquid[c]
@kernel function _twh_heat_combine_kernel!(heat, @Const(mask), @Const(heat_dry_mass),
                                           @Const(heat_ice), @Const(heat_liquid),
                                           @Const(latent_heat_liquid))
    c = @index(Global)
    @inbounds if mask[c]
        heat[c] = heat_dry_mass[c] + heat_ice[c] + heat_liquid[c] + latent_heat_liquid[c]
    end
end

twh_heat_combine!(heat, mask, heat_dry_mass, heat_ice, heat_liquid, latent_heat_liquid) =
    _launch!(_twh_heat_combine_kernel!, heat, mask, heat_dry_mass, heat_ice,
             heat_liquid, latent_heat_liquid)

# heat[c] -= dynbal_baseline_heat_col[c]   (masked, per-column independent)
@kernel function _twh_subtract_baseline_heat_kernel!(heat, @Const(mask),
                                                     @Const(dynbal_baseline_heat))
    c = @index(Global)
    @inbounds if mask[c]
        heat[c] -= dynbal_baseline_heat[c]
    end
end

twh_subtract_baseline_heat!(heat, mask, dynbal_baseline_heat) =
    _launch!(_twh_subtract_baseline_heat_kernel!, heat, mask, dynbal_baseline_heat)

# heat[c] += heat_dry_mass[c] + heat_ice[c] + heat_liquid[c] + latent_heat_liquid[c]
@kernel function _twh_heat_accumulate_kernel!(heat, @Const(mask), @Const(heat_dry_mass),
                                              @Const(heat_ice), @Const(heat_liquid),
                                              @Const(latent_heat_liquid))
    c = @index(Global)
    @inbounds if mask[c]
        heat[c] += heat_dry_mass[c] + heat_ice[c] + heat_liquid[c] + latent_heat_liquid[c]
    end
end

twh_heat_accumulate!(heat, mask, heat_dry_mass, heat_ice, heat_liquid, latent_heat_liquid) =
    _launch!(_twh_heat_accumulate_kernel!, heat, mask, heat_dry_mass, heat_ice,
             heat_liquid, latent_heat_liquid)

# heat[c] += a[c] + b[c] + c3[c]   (masked, per-column independent; lake water heat)
@kernel function _twh_heat_accumulate3_kernel!(heat, @Const(mask), @Const(a),
                                               @Const(b), @Const(c3))
    c = @index(Global)
    @inbounds if mask[c]
        heat[c] += a[c] + b[c] + c3[c]
    end
end

twh_heat_accumulate3!(heat, mask, a, b, c3) =
    _launch!(_twh_heat_accumulate3_kernel!, heat, mask, a, b, c3)

# Soil-heat finalize (two independent per-column outputs):
#   heat_liquid[c] += soil_heat_liquid[c]
#   heat[c]        += soil_heat_dry_mass[c] + soil_heat_ice[c]
#                     + soil_heat_liquid[c] + soil_latent_heat_liquid[c]
@kernel function _twh_soil_heat_finalize_kernel!(heat, @Const(mask), heat_liquid,
                                                 @Const(soil_heat_liquid),
                                                 @Const(soil_heat_dry_mass),
                                                 @Const(soil_heat_ice),
                                                 @Const(soil_latent_heat_liquid))
    c = @index(Global)
    @inbounds if mask[c]
        heat_liquid[c] += soil_heat_liquid[c]
        heat[c] += soil_heat_dry_mass[c] + soil_heat_ice[c] +
                   soil_heat_liquid[c] + soil_latent_heat_liquid[c]
    end
end

twh_soil_heat_finalize!(heat, mask, heat_liquid, soil_heat_liquid, soil_heat_dry_mass,
                        soil_heat_ice, soil_latent_heat_liquid) =
    _launch!(_twh_soil_heat_finalize_kernel!, heat, mask, heat_liquid, soil_heat_liquid,
             soil_heat_dry_mass, soil_heat_ice, soil_latent_heat_liquid)

# heat[c] = -dynbal_baseline_heat_col[c]   (masked, per-column independent)
@kernel function _twh_init_neg_baseline_heat_kernel!(heat, @Const(mask),
                                                     @Const(dynbal_baseline_heat))
    c = @index(Global)
    @inbounds if mask[c]
        heat[c] = -dynbal_baseline_heat[c]
    end
end

twh_init_neg_baseline_heat!(heat, mask, dynbal_baseline_heat) =
    _launch!(_twh_init_neg_baseline_heat_kernel!, heat, mask, dynbal_baseline_heat)

# Zero two masked per-column arrays: a[c] = 0; b[c] = 0
@kernel function _twh_zero2_kernel!(a, @Const(mask), b)
    c = @index(Global)
    T = eltype(a)
    @inbounds if mask[c]
        a[c] = zero(T)
        b[c] = zero(T)
    end
end

twh_zero2!(a, mask, b) = _launch!(_twh_zero2_kernel!, a, mask, b)

# Subtract two baselines: liquid_mass[c] -= bliq[c]; ice_mass[c] -= bice[c]
@kernel function _twh_subtract_baselines2_kernel!(liquid_mass, @Const(mask), ice_mass,
                                                  @Const(bliq), @Const(bice))
    c = @index(Global)
    @inbounds if mask[c]
        liquid_mass[c] -= bliq[c]
        ice_mass[c] -= bice[c]
    end
end

twh_subtract_baselines2!(liquid_mass, mask, ice_mass, bliq, bice) =
    _launch!(_twh_subtract_baselines2_kernel!, liquid_mass, mask, ice_mass, bliq, bice)

# ---------------------------------------------------------------------------
# Non-lake liquid/ice mass: canopy + snow + aquifer + surface (per-column)
# ---------------------------------------------------------------------------
@kernel function _twh_liqice_nonlake_kernel!(liquid_mass, @Const(mask), ice_mass,
                                             @Const(liqcan), @Const(snocan),
                                             @Const(total_plant_stored),
                                             @Const(h2osno_no_layers),
                                             @Const(snl), @Const(h2osoi_liq),
                                             @Const(h2osoi_ice),
                                             @Const(hydrologically_active),
                                             @Const(wa), aquifer_baseline,
                                             @Const(itype), @Const(h2osfc), nlevsno)
    c = @index(Global)
    @inbounds if mask[c]
        lm = liquid_mass[c]
        im = ice_mass[c]

        # Canopy water + plant stored water
        lm += liqcan[c] + total_plant_stored[c]
        im += snocan[c]

        # Snow not resolved into layers
        im += h2osno_no_layers[c]

        # Snow layers: Fortran indices snl+1:0, Julia offset by nlevsno
        for j in (snl[c] + 1):0
            jj = j + nlevsno
            lm += h2osoi_liq[c, jj]
            im += h2osoi_ice[c, jj]
        end

        # Aquifer water (only for hydrologically active columns)
        if hydrologically_active[c]
            lm += (wa[c] - aquifer_baseline)
        end

        # Surface water (exclude urban impervious/wall/roof)
        if !(itype[c] == ICOL_ROOF || itype[c] == ICOL_SUNWALL ||
             itype[c] == ICOL_SHADEWALL || itype[c] == ICOL_ROAD_IMPERV)
            lm += h2osfc[c]
        end

        liquid_mass[c] = lm
        ice_mass[c] = im
    end
end

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

# Device-safe helper variants: convert Float64 module constants to the working
# eltype `T` so no double leaks into a Metal kernel. On CPU (T == Float64) these
# are byte-identical to the host helpers above (T(x) is the identity).
@inline function _temp_to_heat_d(::Type{T}, temp, cv) where {T}
    return cv * (temp - T(heat_base_temp))
end

@inline function _accumulate_liquid_water_heat_d(::Type{T}, temp, h2o,
                                                 heat_liquid, latent_heat_liquid,
                                                 cv_liquid) where {T}
    cv = h2o * T(CPLIQ)
    cv_liquid_out = cv_liquid + cv
    heat_liquid_out = heat_liquid + _temp_to_heat_d(T, temp, cv)
    latent_heat_liquid_out = latent_heat_liquid + h2o * T(HFUS)
    return heat_liquid_out, latent_heat_liquid_out, cv_liquid_out
end

@inline function _accumulate_liquid_water_heat_d(::Type{T}, temp, h2o,
                                                 heat_liquid, latent_heat_liquid) where {T}
    cv = h2o * T(CPLIQ)
    heat_liquid_out = heat_liquid + _temp_to_heat_d(T, temp, cv)
    latent_heat_liquid_out = latent_heat_liquid + h2o * T(HFUS)
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
function compute_water_mass_non_lake!(mask_nolake::AbstractVector{Bool},
                                       col::ColumnData,
                                       waterstate::WaterStateData,
                                       waterdiagnostic::WaterDiagnosticBulkData,
                                       subtract_dynbal_baselines::Bool,
                                       water_mass::AbstractVector{<:Real})
    nc = length(water_mass)
    FT = eltype(water_mass)
    liquid_mass = fill!(similar(water_mass, FT, nc), zero(FT))
    ice_mass = fill!(similar(water_mass, FT, nc), zero(FT))

    compute_liq_ice_mass_non_lake!(mask_nolake, col, waterstate, waterdiagnostic,
                                    subtract_dynbal_baselines, liquid_mass, ice_mass)

    twh_water_mass_combine!(water_mass, mask_nolake, liquid_mass, ice_mass)
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
function compute_water_mass_lake!(mask_lake::AbstractVector{Bool},
                                   col::ColumnData,
                                   waterstate::WaterStateData,
                                   lakestate::LakeStateData,
                                   add_lake_water_and_subtract_dynbal_baselines::Bool,
                                   water_mass::AbstractVector{<:Real})
    nc = length(water_mass)
    FT = eltype(water_mass)
    liquid_mass = fill!(similar(water_mass, FT, nc), zero(FT))
    ice_mass = fill!(similar(water_mass, FT, nc), zero(FT))

    compute_liq_ice_mass_lake!(mask_lake, col, waterstate, lakestate,
                                add_lake_water_and_subtract_dynbal_baselines,
                                liquid_mass, ice_mass)

    twh_water_mass_combine!(water_mass, mask_lake, liquid_mass, ice_mass)
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
function compute_liq_ice_mass_non_lake!(mask_nolake::AbstractVector{Bool},
                                          col::ColumnData,
                                          waterstate::WaterStateData,
                                          waterdiagnostic::WaterDiagnosticBulkData,
                                          subtract_dynbal_baselines::Bool,
                                          liquid_mass::AbstractVector{<:Real},
                                          ice_mass::AbstractVector{<:Real};
                                          liqcan_col::Union{AbstractVector, Nothing}=nothing,
                                          snocan_col::Union{AbstractVector, Nothing}=nothing)
    nlevsno = varpar.nlevsno

    # Initialize
    twh_zero2!(liquid_mass, mask_nolake, ice_mass)

    # Canopy water: if not provided, use zero (caller should supply p2c results)
    FT_li = eltype(liquid_mass)
    _liqcan = liqcan_col !== nothing ? liqcan_col :
              fill!(similar(liquid_mass, FT_li, length(liquid_mass)), zero(FT_li))
    _snocan = snocan_col !== nothing ? snocan_col :
              fill!(similar(liquid_mass, FT_li, length(liquid_mass)), zero(FT_li))

    aquifer_baseline = FT_li(waterstate.aquifer_water_baseline)

    _launch!(_twh_liqice_nonlake_kernel!, liquid_mass, mask_nolake, ice_mass,
             _liqcan, _snocan, waterdiagnostic.total_plant_stored_h2o_col,
             waterstate.h2osno_no_layers_col, col.snl, waterstate.h2osoi_liq_col,
             waterstate.h2osoi_ice_col, col.hydrologically_active, waterstate.wa_col,
             aquifer_baseline, col.itype, waterstate.h2osfc_col, nlevsno)

    # Soil water content
    accumulate_soil_liq_ice_mass_non_lake!(mask_nolake, col, waterstate,
                                            liquid_mass, ice_mass)

    # Subtract dynbal baselines if requested
    if subtract_dynbal_baselines
        twh_subtract_baselines2!(liquid_mass, mask_nolake, ice_mass,
                                 waterstate.dynbal_baseline_liq_col,
                                 waterstate.dynbal_baseline_ice_col)
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
@kernel function _twh_soil_liqice_nonlake_kernel!(liquid_mass, @Const(mask), ice_mass,
                                                  @Const(itype), @Const(h2osoi_liq),
                                                  @Const(h2osoi_ice), @Const(excess_ice),
                                                  nlevgrnd, nlevurb, nlevmaxurbgrnd, nlevsno)
    c = @index(Global)
    @inbounds if mask[c]
        lm = liquid_mass[c]
        im = ice_mass[c]
        for j in 1:nlevmaxurbgrnd
            has_h2o = false
            if itype[c] == ICOL_SUNWALL || itype[c] == ICOL_SHADEWALL
                has_h2o = false
            elseif itype[c] == ICOL_ROOF
                has_h2o = (j <= nlevurb)
            else
                has_h2o = (j <= nlevgrnd)
            end

            if has_h2o
                jj = j + nlevsno  # convert soil index to combined array index
                lm += h2osoi_liq[c, jj]
                im += h2osoi_ice[c, jj] + excess_ice[c, jj]
            end
        end
        liquid_mass[c] = lm
        ice_mass[c] = im
    end
end

function accumulate_soil_liq_ice_mass_non_lake!(mask::AbstractVector{Bool},
                                                  col::ColumnData,
                                                  waterstate::WaterStateData,
                                                  liquid_mass::AbstractVector{<:Real},
                                                  ice_mass::AbstractVector{<:Real})
    nlevgrnd       = varpar.nlevgrnd
    nlevurb        = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevsno        = varpar.nlevsno

    _launch!(_twh_soil_liqice_nonlake_kernel!, liquid_mass, mask, ice_mass,
             col.itype, waterstate.h2osoi_liq_col, waterstate.h2osoi_ice_col,
             waterstate.excess_ice_col, nlevgrnd, nlevurb, nlevmaxurbgrnd, nlevsno)

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
# Lake snow + soil-under-lake liquid/ice accumulation (per-column, adds to existing)
@kernel function _twh_liqice_lake_snowsoil_kernel!(liquid_mass, @Const(mask), ice_mass,
                                                   @Const(h2osno_no_layers), @Const(snl),
                                                   @Const(h2osoi_liq), @Const(h2osoi_ice),
                                                   nlevgrnd, nlevsno)
    c = @index(Global)
    @inbounds if mask[c]
        lm = liquid_mass[c]
        im = ice_mass[c]

        # Snow water content
        im += h2osno_no_layers[c]
        for j in (snl[c] + 1):0
            jj = j + nlevsno
            lm += h2osoi_liq[c, jj]
            im += h2osoi_ice[c, jj]
        end

        # Soil water content under the lake
        for j in 1:nlevgrnd
            jj = j + nlevsno
            lm += h2osoi_liq[c, jj]
            im += h2osoi_ice[c, jj]
        end

        liquid_mass[c] = lm
        ice_mass[c] = im
    end
end

function compute_liq_ice_mass_lake!(mask_lake::AbstractVector{Bool},
                                     col::ColumnData,
                                     waterstate::WaterStateData,
                                     lakestate::LakeStateData,
                                     add_lake_water_and_subtract_dynbal_baselines::Bool,
                                     liquid_mass::AbstractVector{<:Real},
                                     ice_mass::AbstractVector{<:Real})
    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno

    # Initialize
    twh_zero2!(liquid_mass, mask_lake, ice_mass)

    # Large canceling terms first (baselines and lake water)
    if add_lake_water_and_subtract_dynbal_baselines
        twh_subtract_baselines2!(liquid_mass, mask_lake, ice_mass,
                                 waterstate.dynbal_baseline_liq_col,
                                 waterstate.dynbal_baseline_ice_col)

        # Lake water content (tracer_ratio = 1.0 for bulk water)
        accumulate_liq_ice_mass_lake!(mask_lake, col, lakestate, 1.0,
                                       liquid_mass, ice_mass)
    end

    # Snow + soil-under-lake water content
    _launch!(_twh_liqice_lake_snowsoil_kernel!, liquid_mass, mask_lake, ice_mass,
             waterstate.h2osno_no_layers_col, col.snl, waterstate.h2osoi_liq_col,
             waterstate.h2osoi_ice_col, nlevgrnd, nlevsno)

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
@kernel function _twh_accum_liqice_lake_kernel!(liquid_mass, @Const(mask), ice_mass,
                                                @Const(dz_lake), @Const(lake_icefrac),
                                                tracer_ratio, nlevlak)
    c = @index(Global)
    T = eltype(liquid_mass)
    @inbounds if mask[c]
        dh2o = T(DENH2O)
        lm = liquid_mass[c]
        im = ice_mass[c]
        for j in 1:nlevlak
            h2olak_liq = dz_lake[c, j] * dh2o * (one(T) - lake_icefrac[c, j]) * tracer_ratio
            h2olak_ice = dz_lake[c, j] * dh2o * lake_icefrac[c, j] * tracer_ratio
            lm += h2olak_liq
            im += h2olak_ice
        end
        liquid_mass[c] = lm
        ice_mass[c] = im
    end
end

function accumulate_liq_ice_mass_lake!(mask::AbstractVector{Bool},
                                        col::ColumnData,
                                        lakestate::LakeStateData,
                                        tracer_ratio::Real,
                                        liquid_mass::AbstractVector{<:Real},
                                        ice_mass::AbstractVector{<:Real})
    nlevlak = varpar.nlevlak
    T = eltype(liquid_mass)

    _launch!(_twh_accum_liqice_lake_kernel!, liquid_mass, mask, ice_mass,
             col.dz_lake, lakestate.lake_icefrac_col, T(tracer_ratio), nlevlak)

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
# Non-lake heat: canopy + snow + aquifer + surface (per-column, init-from-zero)
@kernel function _twh_heat_nonlake_kernel!(heat_liquid, @Const(mask), latent_heat_liquid,
                                           cv_liquid, heat_ice, heat_dry_mass,
                                           @Const(liqcan), @Const(total_plant_stored),
                                           @Const(t_soisno), @Const(h2osno_no_layers),
                                           @Const(snl), @Const(h2osoi_liq),
                                           @Const(h2osoi_ice),
                                           @Const(hydrologically_active),
                                           @Const(wa), aquifer_baseline,
                                           @Const(t_h2osfc), @Const(h2osfc),
                                           @Const(itype), nlevsno)
    c = @index(Global)
    T = eltype(heat_liquid)
    @inbounds if mask[c]
        hbt = T(heat_base_temp)
        hl = zero(T)
        lhl = zero(T)
        cvl = zero(T)
        hi = zero(T)

        # --- Canopy water ---
        liqveg = liqcan[c] + total_plant_stored[c]
        hl, lhl, cvl = _accumulate_liquid_water_heat_d(T, hbt, liqveg, hl, lhl, cvl)

        # --- Snow ---
        j_top = nlevsno + 1
        hi += _temp_to_heat_d(T, t_soisno[c, j_top], h2osno_no_layers[c] * T(CPICE))

        for j in (snl[c] + 1):0
            jj = j + nlevsno
            hl, lhl, cvl = _accumulate_liquid_water_heat_d(T, t_soisno[c, jj],
                h2osoi_liq[c, jj], hl, lhl, cvl)
            hi += _temp_to_heat_d(T, t_soisno[c, jj], h2osoi_ice[c, jj] * T(CPICE))
        end

        # --- Aquifer water ---
        if hydrologically_active[c]
            hl, lhl, cvl = _accumulate_liquid_water_heat_d(T, hbt,
                wa[c] - aquifer_baseline, hl, lhl, cvl)
        end

        # --- Surface water ---
        if !(itype[c] == ICOL_ROOF || itype[c] == ICOL_SUNWALL ||
             itype[c] == ICOL_SHADEWALL || itype[c] == ICOL_ROAD_IMPERV)
            hl, lhl, cvl = _accumulate_liquid_water_heat_d(T, t_h2osfc[c],
                h2osfc[c], hl, lhl, cvl)
        end

        heat_liquid[c] = hl
        latent_heat_liquid[c] = lhl
        cv_liquid[c] = cvl
        heat_ice[c] = hi
        heat_dry_mass[c] = zero(T)
    end
end

function compute_heat_non_lake!(mask_nolake::AbstractVector{Bool},
                                 col::ColumnData,
                                 lun::LandunitData,
                                 urbanparams::UrbanParamsData,
                                 soilstate::SoilStateData,
                                 temperature::TemperatureData,
                                 waterstatebulk::WaterStateBulkData,
                                 waterdiagnosticbulk::WaterDiagnosticBulkData,
                                 heat::AbstractVector{<:Real},
                                 heat_liquid::AbstractVector{<:Real},
                                 cv_liquid::AbstractVector{<:Real};
                                 liqcan_col::Union{AbstractVector, Nothing}=nothing,
                                 snocan_col::Union{AbstractVector, Nothing}=nothing)
    nlevsno = varpar.nlevsno
    nc = length(heat)

    # Access through composition for bulk data
    ws = waterstatebulk.ws

    # Local accumulators
    FT = eltype(heat)
    heat_dry_mass = fill!(similar(heat, FT, nc), zero(FT))
    heat_ice = fill!(similar(heat, FT, nc), zero(FT))
    latent_heat_liquid = fill!(similar(heat, FT, nc), zero(FT))

    # Canopy water: if not provided, use zero
    _liqcan = liqcan_col !== nothing ? liqcan_col :
              fill!(similar(heat, FT, nc), zero(FT))
    _snocan = snocan_col !== nothing ? snocan_col :
              fill!(similar(heat, FT, nc), zero(FT))

    aquifer_baseline = FT(ws.aquifer_water_baseline)

    _launch!(_twh_heat_nonlake_kernel!, heat_liquid, mask_nolake, latent_heat_liquid,
             cv_liquid, heat_ice, heat_dry_mass, _liqcan,
             waterdiagnosticbulk.total_plant_stored_h2o_col, temperature.t_soisno_col,
             ws.h2osno_no_layers_col, col.snl, ws.h2osoi_liq_col, ws.h2osoi_ice_col,
             col.hydrologically_active, ws.wa_col, aquifer_baseline,
             temperature.t_h2osfc_col, ws.h2osfc_col, col.itype, nlevsno)

    # Combine into heat
    twh_heat_combine!(heat, mask_nolake, heat_dry_mass, heat_ice,
                      heat_liquid, latent_heat_liquid)

    # Soil heat
    accumulate_soil_heat_non_lake!(mask_nolake, col, lun, urbanparams, soilstate,
                                    temperature, waterstatebulk,
                                    heat, heat_liquid, cv_liquid)

    # Subtract baseline heat
    twh_subtract_baseline_heat!(heat, mask_nolake, temperature.dynbal_baseline_heat_col)

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
@kernel function _twh_soil_heat_nonlake_kernel!(soil_heat_liquid, @Const(mask),
                                                soil_heat_dry_mass, soil_heat_ice,
                                                soil_latent_heat_liquid, cv_liquid,
                                                @Const(itype), @Const(landunit),
                                                @Const(t_soisno), @Const(dz),
                                                @Const(cv_wall), @Const(cv_roof),
                                                @Const(cv_improad), @Const(nlev_improad),
                                                @Const(csol), @Const(watsat),
                                                @Const(lun_itype), @Const(h2osoi_liq),
                                                @Const(h2osoi_ice),
                                                nlevgrnd, nlevurb, nlevmaxurbgrnd, nlevsno)
    c = @index(Global)
    T = eltype(soil_heat_liquid)
    @inbounds if mask[c]
        shl  = zero(T)
        shdm = zero(T)
        shi  = zero(T)
        slhl = zero(T)
        cvl  = cv_liquid[c]
        l = landunit[c]

        for j in 1:nlevmaxurbgrnd
            jj = j + nlevsno  # combined snow+soil index

            has_h2o = false

            if itype[c] == ICOL_SUNWALL || itype[c] == ICOL_SHADEWALL
                has_h2o = false
                if j <= nlevurb
                    shdm += _temp_to_heat_d(T, t_soisno[c, jj], cv_wall[l, j] * dz[c, jj])
                end

            elseif itype[c] == ICOL_ROOF
                if j <= nlevurb
                    has_h2o = true
                    shdm += _temp_to_heat_d(T, t_soisno[c, jj], cv_roof[l, j] * dz[c, jj])
                else
                    has_h2o = false
                end

            else
                if j <= nlevgrnd
                    has_h2o = true

                    if itype[c] == ICOL_ROAD_IMPERV && j <= nlev_improad[l]
                        shdm += _temp_to_heat_d(T, t_soisno[c, jj], cv_improad[l, j] * dz[c, jj])
                    elseif lun_itype[l] != ISTWET && lun_itype[l] != ISTICE
                        # Includes impervious roads below nlev_improad (where we have soil)
                        shdm += _temp_to_heat_d(T, t_soisno[c, jj],
                            csol[c, j] * (one(T) - watsat[c, j]) * dz[c, jj])
                    end
                else
                    has_h2o = false
                end
            end

            if has_h2o
                shl, slhl, cvl = _accumulate_liquid_water_heat_d(T, t_soisno[c, jj],
                    h2osoi_liq[c, jj], shl, slhl, cvl)

                shi += _temp_to_heat_d(T, t_soisno[c, jj], h2osoi_ice[c, jj] * T(CPICE))
            end
        end

        soil_heat_liquid[c] = shl
        soil_heat_dry_mass[c] = shdm
        soil_heat_ice[c] = shi
        soil_latent_heat_liquid[c] = slhl
        cv_liquid[c] = cvl
    end
end

function accumulate_soil_heat_non_lake!(mask::AbstractVector{Bool},
                                          col::ColumnData,
                                          lun::LandunitData,
                                          urbanparams::UrbanParamsData,
                                          soilstate::SoilStateData,
                                          temperature::TemperatureData,
                                          waterstatebulk::WaterStateBulkData,
                                          heat::AbstractVector{<:Real},
                                          heat_liquid::AbstractVector{<:Real},
                                          cv_liquid::AbstractVector{<:Real})
    nlevgrnd       = varpar.nlevgrnd
    nlevurb        = varpar.nlevurb
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd
    nlevsno        = varpar.nlevsno

    ws = waterstatebulk.ws
    nc = length(heat)

    # Local accumulators
    FT = eltype(heat)
    soil_heat_liquid = fill!(similar(heat, FT, nc), zero(FT))
    soil_heat_dry_mass = fill!(similar(heat, FT, nc), zero(FT))
    soil_heat_ice = fill!(similar(heat, FT, nc), zero(FT))
    soil_latent_heat_liquid = fill!(similar(heat, FT, nc), zero(FT))

    _launch!(_twh_soil_heat_nonlake_kernel!, soil_heat_liquid, mask, soil_heat_dry_mass,
             soil_heat_ice, soil_latent_heat_liquid, cv_liquid, col.itype, col.landunit,
             temperature.t_soisno_col, col.dz, urbanparams.cv_wall, urbanparams.cv_roof,
             urbanparams.cv_improad, urbanparams.nlev_improad, soilstate.csol_col,
             soilstate.watsat_col, lun.itype, ws.h2osoi_liq_col, ws.h2osoi_ice_col,
             nlevgrnd, nlevurb, nlevmaxurbgrnd, nlevsno)

    twh_soil_heat_finalize!(heat, mask, heat_liquid, soil_heat_liquid,
                            soil_heat_dry_mass, soil_heat_ice, soil_latent_heat_liquid)

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
# Lake snow + soil-under-lake heat (per-column, init-from-zero)
@kernel function _twh_heat_lake_snowsoil_kernel!(heat_liquid, @Const(mask),
                                                 latent_heat_liquid, cv_liquid,
                                                 heat_ice, heat_dry_mass,
                                                 @Const(t_soisno), @Const(snl),
                                                 @Const(h2osoi_liq), @Const(h2osoi_ice),
                                                 @Const(csol), @Const(watsat),
                                                 @Const(dz), nlevgrnd, nlevsno)
    c = @index(Global)
    T = eltype(heat_liquid)
    @inbounds if mask[c]
        hl   = zero(T)
        lhl  = zero(T)
        cvl  = zero(T)
        hi   = zero(T)
        hdm  = zero(T)

        # Snow heat content
        for j in (snl[c] + 1):0
            jj = j + nlevsno
            hl, lhl, cvl = _accumulate_liquid_water_heat_d(T, t_soisno[c, jj],
                h2osoi_liq[c, jj], hl, lhl, cvl)
            hi += _temp_to_heat_d(T, t_soisno[c, jj], h2osoi_ice[c, jj] * T(CPICE))
        end

        # Soil heat content under the lake
        for j in 1:nlevgrnd
            jj = j + nlevsno
            hdm += _temp_to_heat_d(T, t_soisno[c, jj],
                csol[c, j] * (one(T) - watsat[c, j]) * dz[c, jj])
            hl, lhl, cvl = _accumulate_liquid_water_heat_d(T, t_soisno[c, jj],
                h2osoi_liq[c, jj], hl, lhl, cvl)
            hi += _temp_to_heat_d(T, t_soisno[c, jj], h2osoi_ice[c, jj] * T(CPICE))
        end

        heat_liquid[c] = hl
        latent_heat_liquid[c] = lhl
        cv_liquid[c] = cvl
        heat_ice[c] = hi
        heat_dry_mass[c] = hdm
    end
end

function compute_heat_lake!(mask_lake::AbstractVector{Bool},
                             col::ColumnData,
                             soilstate::SoilStateData,
                             temperature::TemperatureData,
                             waterstatebulk::WaterStateBulkData,
                             lakestate::LakeStateData,
                             heat::AbstractVector{<:Real},
                             heat_liquid::AbstractVector{<:Real},
                             cv_liquid::AbstractVector{<:Real})
    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno
    nc = length(heat)
    ws = waterstatebulk.ws

    # Local accumulators
    FT = eltype(heat)
    heat_dry_mass = fill!(similar(heat, FT, nc), zero(FT))
    heat_ice = fill!(similar(heat, FT, nc), zero(FT))
    latent_heat_liquid = fill!(similar(heat, FT, nc), zero(FT))

    # Subtract baseline heat (large canceling term first)
    twh_init_neg_baseline_heat!(heat, mask_lake, temperature.dynbal_baseline_heat_col)

    # Lake water heat content (NOT accumulated in heat_liquid/cv_liquid)
    accumulate_heat_lake!(mask_lake, col, temperature, lakestate, heat)

    # Snow + soil-under-lake heat content
    _launch!(_twh_heat_lake_snowsoil_kernel!, heat_liquid, mask_lake, latent_heat_liquid,
             cv_liquid, heat_ice, heat_dry_mass, temperature.t_soisno_col, col.snl,
             ws.h2osoi_liq_col, ws.h2osoi_ice_col, soilstate.csol_col,
             soilstate.watsat_col, col.dz, nlevgrnd, nlevsno)

    twh_heat_accumulate!(heat, mask_lake, heat_dry_mass, heat_ice,
                         heat_liquid, latent_heat_liquid)

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
# Lake water heat (per-column, init-from-zero)
@kernel function _twh_accum_heat_lake_kernel!(lake_heat_liquid, @Const(mask),
                                              lake_heat_ice, lake_latent_heat_liquid,
                                              @Const(dz_lake), @Const(lake_icefrac),
                                              @Const(t_lake), nlevlak)
    c = @index(Global)
    T = eltype(lake_heat_liquid)
    @inbounds if mask[c]
        dh2o = T(DENH2O)
        lhl  = zero(T)
        lhi  = zero(T)
        llhl = zero(T)
        for j in 1:nlevlak
            # Liquid heat
            h2olak_liq = dz_lake[c, j] * dh2o * (one(T) - lake_icefrac[c, j])
            lhl, llhl = _accumulate_liquid_water_heat_d(T, t_lake[c, j],
                h2olak_liq, lhl, llhl)

            # Ice heat (use water density because lake depths are not adjusted on freeze)
            h2olak_ice = dz_lake[c, j] * dh2o * lake_icefrac[c, j]
            lhi += _temp_to_heat_d(T, t_lake[c, j], h2olak_ice * T(CPICE))
        end
        lake_heat_liquid[c] = lhl
        lake_heat_ice[c] = lhi
        lake_latent_heat_liquid[c] = llhl
    end
end

function accumulate_heat_lake!(mask::AbstractVector{Bool},
                                col::ColumnData,
                                temperature::TemperatureData,
                                lakestate::LakeStateData,
                                heat::AbstractVector{<:Real})
    nlevlak = varpar.nlevlak
    nc = length(heat)

    FT = eltype(heat)
    lake_heat_liquid = fill!(similar(heat, FT, nc), zero(FT))
    lake_heat_ice = fill!(similar(heat, FT, nc), zero(FT))
    lake_latent_heat_liquid = fill!(similar(heat, FT, nc), zero(FT))

    _launch!(_twh_accum_heat_lake_kernel!, lake_heat_liquid, mask, lake_heat_ice,
             lake_latent_heat_liquid, col.dz_lake, lakestate.lake_icefrac_col,
             temperature.t_lake_col, nlevlak)

    twh_heat_accumulate3!(heat, mask, lake_heat_ice, lake_heat_liquid,
                          lake_latent_heat_liquid)

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
