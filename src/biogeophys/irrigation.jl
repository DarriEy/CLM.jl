# ==========================================================================
# Ported from: src/biogeophys/IrrigationMod.F90
# Irrigation module
#
# Calculates irrigation flux. Two main entry points:
#
#   - calc_irrigation_needed!: Determine whether and how much irrigation is
#     needed. Should be called once per timestep.
#
#   - calc_irrigation_fluxes!: Apply irrigation fluxes (withdrawal and
#     application). Should be called once per timestep.
#
# Public types:
#   IrrigationParamsData  — Irrigation parameters
#   IrrigationData        — Irrigation state/flux data
#
# Public functions:
#   irrigation_init!               — Initialize irrigation data
#   irrigation_init_for_testing!   — Initialize for unit testing
#   irrigation_init_allocate!      — Allocate arrays
#   irrigation_init_cold!          — Cold start initialization
#   irrigation_clean!              — Deallocate memory
#   calc_irrigation_needed!        — Determine irrigation need
#   calc_irrigation_fluxes!        — Apply irrigation fluxes
#   calc_bulk_withdrawals!         — Calculate bulk water withdrawals
#   calc_total_gw_uncon_irrig!     — Sum unconfined GW irrigation by layer
#   calc_application_fluxes!       — Set drip/sprinkler application fluxes
#   calc_deficit_volr_limited!     — Limit deficit by river volume
#   relsat_to_h2osoi               — Convert relative saturation to kg/m2
#   calc_irrig_nsteps_per_day      — Compute irrigation steps per day
#   point_needs_check_for_irrig    — Check if patch needs irrigation check
#   set_irrig_method!              — Set irrigation method per patch
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

# Soil matric potential at wilting point (mm)
const WILTING_POINT_SMP = -150000.0

# Conversion factor: m3/km2 to mm
const M3_OVER_KM2_TO_MM = 1.0e-3

# Irrigation methods
const IRRIG_METHOD_UNSET     = 0
const IRRIG_METHOD_DRIP      = 1
const IRRIG_METHOD_SPRINKLER = 2

# ---------------------------------------------------------------------------
# IrrigationParamsData — Irrigation parameters
# Ported from: irrigation_params_type in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    IrrigationParamsData

Irrigation parameters structure.

Ported from `irrigation_params_type` in `IrrigationMod.F90`.
"""
Base.@kwdef mutable struct IrrigationParamsData
    # Minimum LAI for irrigation
    irrig_min_lai::Float64 = 0.0
    # Time of day to check whether we need irrigation, seconds (0 = midnight)
    irrig_start_time::Int = 21600
    # Desired amount of time to irrigate per day (sec)
    irrig_length::Int = 14400
    # Target soil matric potential for irrigation (mm)
    irrig_target_smp::Float64 = -3400.0
    # Soil depth to which we measure for irrigation (m)
    irrig_depth::Float64 = 0.6
    # Threshold fraction for irrigation trigger (0 = wilting point, 1 = target)
    irrig_threshold_fraction::Float64 = 1.0
    # Threshold for river water volume below which irrigation is shut off (fraction)
    irrig_river_volume_threshold::Float64 = 0.1
    # Whether irrigation is limited based on river storage (when ROF is enabled)
    limit_irrigation_if_rof_enabled::Bool = false
    # Use groundwater supply for irrigation
    use_groundwater_irrigation::Bool = false
    # Default irrigation method (IRRIG_METHOD_DRIP or IRRIG_METHOD_SPRINKLER)
    irrig_method_default::Int = IRRIG_METHOD_DRIP
end

# ---------------------------------------------------------------------------
# IrrigationData — Irrigation state and flux data
# Ported from: irrigation_type in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    IrrigationData

Irrigation state and flux data structure.

Ported from `irrigation_type` in `IrrigationMod.F90`.
"""
Base.@kwdef mutable struct IrrigationData
    # Parameters
    params::IrrigationParamsData = IrrigationParamsData()
    # Land model time step (sec)
    dtime::Int = 0
    # Number of irrigation time steps per day
    irrig_nsteps_per_day::Int = 0
    # Relative saturation at wilting point [col, nlevsoi]
    relsat_wilting_point_col::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    # Relative saturation at irrigation target [col, nlevsoi]
    relsat_target_col::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    # Patch irrigation application method [patch]
    irrig_method_patch::Vector{Int} = Int[]
    # Current irrigation rate from surface water [mm/s] [patch]
    sfc_irrig_rate_patch::Vector{Float64} = Float64[]
    # Current irrigation rate demand, neglecting surface water source limitation [mm/s] [patch]
    irrig_rate_demand_patch::Vector{Float64} = Float64[]
    # Number of time steps for which we still need to irrigate today [patch]
    n_irrig_steps_left_patch::Vector{Int} = Int[]
    # Irrigation flux neglecting surface water source limitation [mm/s] [patch]
    qflx_irrig_demand_patch::Vector{Float64} = Float64[]
end

# ---------------------------------------------------------------------------
# irrigation_init_allocate! — Allocate irrigation data arrays
# Ported from: IrrigationInitAllocate in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    irrigation_init_allocate!(irrig, np, nc, nlevsoi)

Allocate and initialize all irrigation data arrays.

Ported from `IrrigationInitAllocate` in `IrrigationMod.F90`.
"""
function irrigation_init_allocate!(irrig::IrrigationData, np::Int, nc::Int, nlevsoi::Int)
    irrig.qflx_irrig_demand_patch    = fill(NaN, np)
    irrig.relsat_wilting_point_col   = fill(NaN, nc, nlevsoi)
    irrig.relsat_target_col          = fill(NaN, nc, nlevsoi)
    irrig.sfc_irrig_rate_patch       = fill(NaN, np)
    irrig.irrig_rate_demand_patch    = fill(NaN, np)
    irrig.irrig_method_patch         = fill(ISPVAL, np)
    irrig.n_irrig_steps_left_patch   = fill(0, np)
    return nothing
end

# ---------------------------------------------------------------------------
# irrigation_init_history! — Initialize history fields (stub)
# Ported from: IrrigationInitHistory in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    irrigation_init_history!(irrig, np)

Initialize irrigation history fields. Stub: to be filled when histFileMod is ported.

Ported from `IrrigationInitHistory` in `IrrigationMod.F90`.
"""
function irrigation_init_history!(irrig::IrrigationData, np::Int)
    irrig.qflx_irrig_demand_patch .= SPVAL
    return nothing
end

# ---------------------------------------------------------------------------
# irrigation_init_cold! — Cold start initialization
# Ported from: IrrigationInitCold in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    irrigation_init_cold!(irrig, soilstate, swrc, col_data, pftcon_irrigated,
                          irrig_method_surface, patch_data, grc,
                          bounds_c, bounds_p, nlevsoi, dtime)

Do cold-start initialization for irrigation data.

Ported from `IrrigationInitCold` in `IrrigationMod.F90`.
"""
function irrigation_init_cold!(
    irrig::IrrigationData,
    soilstate::SoilStateData,
    swrc::SoilWaterRetentionCurve,
    col_data::ColumnData,
    pftcon_irrigated::Vector{Float64},
    irrig_method_surface::Matrix{Int},
    patch_data::PatchData,
    grc::GridcellData,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    nlevsoi::Int,
    dtime::Int
)
    for j in 1:nlevsoi
        for c in bounds_c
            s_wp = soil_suction_inverse!(swrc, c, j, WILTING_POINT_SMP, soilstate)
            s_wp = min(s_wp, 1.0)
            s_wp = max(s_wp, 0.0)
            irrig.relsat_wilting_point_col[c, j] = s_wp

            s_tgt = soil_suction_inverse!(swrc, c, j, irrig.params.irrig_target_smp, soilstate)
            s_tgt = min(s_tgt, 1.0)
            s_tgt = max(s_tgt, 0.0)
            irrig.relsat_target_col[c, j] = s_tgt
        end
    end

    set_irrig_method!(irrig, pftcon_irrigated, irrig_method_surface,
                      patch_data, grc, bounds_p)

    irrig.dtime = dtime
    irrig.irrig_nsteps_per_day = calc_irrig_nsteps_per_day(irrig.params.irrig_length, dtime)

    for p in bounds_p
        irrig.qflx_irrig_demand_patch[p] = 0.0
    end

    return nothing
end

# ---------------------------------------------------------------------------
# irrigation_init! — Full initialization
# Ported from: IrrigationInit in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    irrigation_init!(irrig, params, soilstate, swrc, col_data,
                     pftcon_irrigated, irrig_method_surface,
                     patch_data, grc, bounds_c, bounds_p, nlevsoi, dtime)

Full initialization: allocate, history, and cold start.

Ported from `IrrigationInit` in `IrrigationMod.F90`.
"""
function irrigation_init!(
    irrig::IrrigationData,
    params::IrrigationParamsData,
    soilstate::SoilStateData,
    swrc::SoilWaterRetentionCurve,
    col_data::ColumnData,
    pftcon_irrigated::Vector{Float64},
    irrig_method_surface::Matrix{Int},
    patch_data::PatchData,
    grc::GridcellData,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    nlevsoi::Int,
    dtime::Int
)
    irrig.params = params
    np = length(bounds_p)
    nc = length(bounds_c)
    irrigation_init_allocate!(irrig, np, nc, nlevsoi)
    irrigation_init_history!(irrig, np)
    irrigation_init_cold!(irrig, soilstate, swrc, col_data, pftcon_irrigated,
                          irrig_method_surface, patch_data, grc,
                          bounds_c, bounds_p, nlevsoi, dtime)
    return nothing
end

# ---------------------------------------------------------------------------
# irrigation_init_for_testing! — Initialize for unit testing
# Ported from: InitForTesting in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    irrigation_init_for_testing!(irrig, params, dtime, relsat_wilting_point,
                                 relsat_target, pftcon_irrigated,
                                 irrig_method_surface, patch_data, grc,
                                 bounds_c, bounds_p, nlevsoi)

Initialize irrigation for unit testing with prescribed internal values.

Ported from `InitForTesting` in `IrrigationMod.F90`.
"""
function irrigation_init_for_testing!(
    irrig::IrrigationData,
    params::IrrigationParamsData,
    dtime::Int,
    relsat_wilting_point::Matrix{Float64},
    relsat_target::Matrix{Float64},
    pftcon_irrigated::Vector{Float64},
    irrig_method_surface::Matrix{Int},
    patch_data::PatchData,
    grc::GridcellData,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    nlevsoi::Int
)
    np = length(bounds_p)
    nc = length(bounds_c)
    irrigation_init_allocate!(irrig, np, nc, nlevsoi)
    irrig.params = params
    irrig.dtime = dtime
    set_irrig_method!(irrig, pftcon_irrigated, irrig_method_surface,
                      patch_data, grc, bounds_p)
    irrig.irrig_nsteps_per_day = calc_irrig_nsteps_per_day(params.irrig_length, dtime)
    irrig.relsat_wilting_point_col .= relsat_wilting_point
    irrig.relsat_target_col .= relsat_target
    return nothing
end

# ---------------------------------------------------------------------------
# irrigation_clean! — Deallocate memory
# Ported from: IrrigationClean in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    irrigation_clean!(irrig)

Deallocate all irrigation data arrays.

Ported from `IrrigationClean` in `IrrigationMod.F90`.
"""
function irrigation_clean!(irrig::IrrigationData)
    irrig.qflx_irrig_demand_patch    = Float64[]
    irrig.relsat_wilting_point_col   = Matrix{Float64}(undef, 0, 0)
    irrig.relsat_target_col          = Matrix{Float64}(undef, 0, 0)
    irrig.sfc_irrig_rate_patch       = Float64[]
    irrig.irrig_rate_demand_patch    = Float64[]
    irrig.irrig_method_patch         = Int[]
    irrig.n_irrig_steps_left_patch   = Int[]
    return nothing
end

# ---------------------------------------------------------------------------
# irrigation_restart! — Handle restart (stub)
# Ported from: Restart in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    irrigation_restart!(irrig)

Handle restart of irrigation variables. Stub: to be filled when restart I/O is ported.

Ported from `Restart` in `IrrigationMod.F90`.
"""
function irrigation_restart!(irrig::IrrigationData)
    return nothing  # Stub
end

# ---------------------------------------------------------------------------
# calc_irrig_nsteps_per_day — Compute number of irrigation steps per day
# Ported from: CalcIrrigNstepsPerDay in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    calc_irrig_nsteps_per_day(irrig_length, dtime) -> Int

Given irrig_length (sec) and dtime (sec), return number of irrigation steps per day
(rounded up).

Ported from `CalcIrrigNstepsPerDay` in `IrrigationMod.F90`.
"""
function calc_irrig_nsteps_per_day(irrig_length::Int, dtime::Int)
    return div(irrig_length + dtime - 1, dtime)  # round up
end

# ---------------------------------------------------------------------------
# set_irrig_method! — Set irrigation method per patch
# Ported from: SetIrrigMethod in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    set_irrig_method!(irrig, pftcon_irrigated, irrig_method_surface,
                      patch_data, grc, bounds_p)

Set `irrig.irrig_method_patch` based on surface dataset values.

Ported from `SetIrrigMethod` in `IrrigationMod.F90`.
"""
function set_irrig_method!(
    irrig::IrrigationData,
    pftcon_irrigated::Vector{Float64},
    irrig_method_surface::Matrix{Int},
    patch_data::PatchData,
    grc::GridcellData,
    bounds_p::UnitRange{Int}
)
    for p in bounds_p
        g = patch_data.gridcell[p]
        m = patch_data.itype[p]

        if m >= 1 && m <= size(irrig_method_surface, 2) &&
           pftcon_irrigated[m] == 1.0
            method_val = irrig_method_surface[g, m]
            if method_val == IRRIG_METHOD_UNSET
                irrig.irrig_method_patch[p] = irrig.params.irrig_method_default
            elseif method_val == IRRIG_METHOD_DRIP || method_val == IRRIG_METHOD_SPRINKLER
                irrig.irrig_method_patch[p] = method_val
            else
                error("Invalid irrigation method specified for patch $p: $method_val")
            end
        else
            irrig.irrig_method_patch[p] = irrig.params.irrig_method_default
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# relsat_to_h2osoi — Convert relative saturation to kg/m2 water
# Ported from: RelsatToH2osoi in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    relsat_to_h2osoi(relsat, eff_porosity, dz) -> Float64

Convert relative saturation to kg/m2 water for a single column and layer.

Ported from `RelsatToH2osoi` in `IrrigationMod.F90`.
"""
function relsat_to_h2osoi(relsat::Float64, eff_porosity::Float64, dz::Float64)
    vol_liq = eff_porosity * relsat
    return vol_liq * DENH2O * dz
end

# ---------------------------------------------------------------------------
# point_needs_check_for_irrig — Determine if a patch needs irrigation check
# Ported from: PointNeedsCheckForIrrig in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    point_needs_check_for_irrig(irrig, pft_type, elai, londeg,
                                 pftcon_irrigated, local_time_sec) -> Bool

Determine whether a given patch needs to be checked for irrigation now.

Ported from `PointNeedsCheckForIrrig` in `IrrigationMod.F90`.
"""
function point_needs_check_for_irrig(
    irrig::IrrigationData,
    pft_type::Int,
    elai::Float64,
    londeg::Float64,
    pftcon_irrigated::Vector{Float64},
    local_time_sec::Int
)
    if pftcon_irrigated[pft_type] == 1.0 && elai > irrig.params.irrig_min_lai
        # Compute seconds since irrig start time, accounting for local time
        # local_time_sec is seconds since midnight at this longitude
        # We compute seconds since (irrig_start_time - dtime)
        start_offset = irrig.params.irrig_start_time - irrig.dtime
        seconds_since = mod(local_time_sec - start_offset, ISECSPDAY)
        if seconds_since < 0
            seconds_since += ISECSPDAY
        end
        if seconds_since < irrig.dtime
            return true
        end
    end
    return false
end

# ---------------------------------------------------------------------------
# p2c_irrig! — Patch-to-column area-weighted average for irrigation
# Ported from: p2c in subgridAveMod.F90 (local version)
# ---------------------------------------------------------------------------

"""
    p2c_irrig!(col_out, patch_in, patch_data, mask_soilc, bounds_c, bounds_p)

Patch-to-column area-weighted average.

Ported from `p2c` in `subgridAveMod.F90`.
"""
function p2c_irrig!(
    col_out::Vector{Float64},
    patch_in::Vector{Float64},
    patch_data::PatchData,
    mask_soilc::BitVector,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int}
)
    for c in bounds_c
        mask_soilc[c] || continue
        col_out[c] = 0.0
    end
    for p in bounds_p
        c = patch_data.column[p]
        (c in bounds_c && mask_soilc[c]) || continue
        col_out[c] += patch_in[p] * patch_data.wtcol[p]
    end
    return nothing
end

# ---------------------------------------------------------------------------
# c2g_irrig! — Column-to-gridcell area-weighted average for irrigation
# Ported from: c2g in subgridAveMod.F90 (local version)
# ---------------------------------------------------------------------------

"""
    c2g_irrig!(garr, carr, col_data, bounds_c, bounds_g)

Column-to-gridcell area-weighted average (unity scaling).

Ported from `c2g` in `subgridAveMod.F90`.
"""
function c2g_irrig!(
    garr::Vector{Float64},
    carr::Vector{Float64},
    col_data::ColumnData,
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int}
)
    for g in bounds_g
        garr[g] = 0.0
    end
    for c in bounds_c
        g = col_data.gridcell[c]
        if g in bounds_g
            garr[g] += carr[c] * col_data.wtgcell[c]
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# calc_deficit_volr_limited! — Limit deficit by river volume
# Ported from: CalcDeficitVolrLimited in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    calc_deficit_volr_limited!(irrig, deficit, volr, deficit_volr_limited,
                               check_for_irrig_col, col_data, grc,
                               bounds_c, bounds_g)

Calculate deficit limited by river volume for each column.

Ported from `CalcDeficitVolrLimited` in `IrrigationMod.F90`.
"""
function calc_deficit_volr_limited!(
    irrig::IrrigationData,
    deficit::Vector{Float64},
    volr::Vector{Float64},
    deficit_volr_limited::Vector{Float64},
    check_for_irrig_col::BitVector,
    col_data::ColumnData,
    grc::GridcellData,
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int}
)
    ng = length(bounds_g)
    deficit_grc = zeros(Float64, length(volr))
    deficit_limited_ratio_grc = ones(Float64, length(volr))

    # Average deficit to gridcell level (unity scaling)
    c2g_irrig!(deficit_grc, deficit, col_data, bounds_c, bounds_g)

    for g in bounds_g
        if volr[g] > 0.0
            available_volr = volr[g] * (1.0 - irrig.params.irrig_river_volume_threshold)
            max_deficit_supported_by_volr = available_volr / grc.area[g] * M3_OVER_KM2_TO_MM
        else
            # Negative volr treated as 0
            max_deficit_supported_by_volr = 0.0
        end

        if deficit_grc[g] > max_deficit_supported_by_volr
            # Inadequate river storage, adjust irrigation demand
            deficit_limited_ratio_grc[g] = max_deficit_supported_by_volr / deficit_grc[g]
        else
            # Adequate river storage, no adjustment
            deficit_limited_ratio_grc[g] = 1.0
        end
    end

    deficit_volr_limited[bounds_c] .= 0.0
    for c in bounds_c
        check_for_irrig_col[c] || continue
        g = col_data.gridcell[c]
        deficit_volr_limited[c] = deficit[c] * deficit_limited_ratio_grc[g]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# calc_irrigation_needed! — Determine irrigation need
# Ported from: CalcIrrigationNeeded in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    calc_irrigation_needed!(irrig, elai, t_soisno, eff_porosity, h2osoi_liq,
                            volr, rof_prognostic, pftcon_irrigated,
                            local_time_sec_patch, col_data, grc,
                            patch_data, mask_exposedveg,
                            bounds_c, bounds_p, bounds_g, nlevsoi)

Calculate whether and how much irrigation is needed for each column.
Does NOT actually set the irrigation flux.

Ported from `CalcIrrigationNeeded` in `IrrigationMod.F90`.
"""
function calc_irrigation_needed!(
    irrig::IrrigationData,
    elai::Vector{Float64},
    t_soisno::Matrix{Float64},
    eff_porosity::Matrix{Float64},
    h2osoi_liq::Matrix{Float64},
    volr::Vector{Float64},
    rof_prognostic::Bool,
    pftcon_irrigated::Vector{Float64},
    local_time_sec_patch::Vector{Int},
    col_data::ColumnData,
    grc::GridcellData,
    patch_data::PatchData,
    mask_exposedveg::BitVector,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    bounds_g::UnitRange{Int},
    nlevsoi::Int
)
    nc = length(bounds_c)

    # Determine which patches need to be checked for irrigation
    check_for_irrig_patch = falses(length(elai))
    check_for_irrig_col = falses(length(check_for_irrig_patch) > 0 ?
                                  maximum(bounds_c) : 0)
    if length(check_for_irrig_col) < maximum(bounds_c)
        check_for_irrig_col = falses(maximum(bounds_c))
    end

    for p in bounds_p
        mask_exposedveg[p] || continue
        g = patch_data.gridcell[p]

        check_for_irrig_patch[p] = point_needs_check_for_irrig(
            irrig, patch_data.itype[p], elai[p], grc.londeg[g],
            pftcon_irrigated, local_time_sec_patch[p])

        if check_for_irrig_patch[p]
            c = patch_data.column[p]
            check_for_irrig_col[c] = true
        end
    end

    # Initialize accumulators for columns that need irrigation check
    h2osoi_liq_tot = zeros(Float64, maximum(bounds_c))
    h2osoi_liq_target_tot = zeros(Float64, maximum(bounds_c))
    h2osoi_liq_wilting_point_tot = zeros(Float64, maximum(bounds_c))
    reached_max_depth = falses(maximum(bounds_c))

    # Measure soil water and see if irrigation is needed
    for j in 1:nlevsoi
        for c in bounds_c
            check_for_irrig_col[c] || continue

            if !reached_max_depth[c]
                if col_data.z[c, j] > irrig.params.irrig_depth
                    reached_max_depth[c] = true
                elseif j > col_data.nbedrock[c]
                    reached_max_depth[c] = true
                elseif t_soisno[c, j] <= TFRZ
                    # Frozen level: don't look below
                    reached_max_depth[c] = true
                else
                    h2osoi_liq_tot[c] += h2osoi_liq[c, j]

                    h2osoi_liq_target = relsat_to_h2osoi(
                        irrig.relsat_target_col[c, j],
                        eff_porosity[c, j],
                        col_data.dz[c, j])
                    h2osoi_liq_target_tot[c] += h2osoi_liq_target

                    h2osoi_liq_wp = relsat_to_h2osoi(
                        irrig.relsat_wilting_point_col[c, j],
                        eff_porosity[c, j],
                        col_data.dz[c, j])
                    h2osoi_liq_wilting_point_tot[c] += h2osoi_liq_wp
                end
            end
        end
    end

    # Compute deficits
    deficit = zeros(Float64, maximum(bounds_c))
    for c in bounds_c
        check_for_irrig_col[c] || continue

        h2osoi_liq_at_threshold = h2osoi_liq_wilting_point_tot[c] +
            irrig.params.irrig_threshold_fraction *
            (h2osoi_liq_target_tot[c] - h2osoi_liq_wilting_point_tot[c])

        if h2osoi_liq_tot[c] < h2osoi_liq_at_threshold
            deficit[c] = h2osoi_liq_target_tot[c] - h2osoi_liq_tot[c]
            if deficit[c] < 0.0
                error("calc_irrigation_needed!: deficit < 0 at column $c. " *
                      "This implies irrigation target < irrigation threshold.")
            end
        else
            deficit[c] = 0.0
        end
    end

    # Limit deficits by available volr if desired
    limit_irrigation = irrig.params.limit_irrigation_if_rof_enabled && rof_prognostic
    deficit_volr_limited = zeros(Float64, maximum(bounds_c))

    if limit_irrigation
        calc_deficit_volr_limited!(irrig, deficit, volr, deficit_volr_limited,
                                   check_for_irrig_col, col_data, grc,
                                   bounds_c, bounds_g)
    else
        deficit_volr_limited[bounds_c] .= deficit[bounds_c]
    end

    # Convert deficits to irrigation rate
    for p in bounds_p
        mask_exposedveg[p] || continue
        c = patch_data.column[p]

        if check_for_irrig_patch[p]
            # Convert from mm to mm/sec
            irrig.sfc_irrig_rate_patch[p] = deficit_volr_limited[c] /
                (irrig.dtime * irrig.irrig_nsteps_per_day)
            irrig.irrig_rate_demand_patch[p] = deficit[c] /
                (irrig.dtime * irrig.irrig_nsteps_per_day)
            irrig.n_irrig_steps_left_patch[p] = irrig.irrig_nsteps_per_day
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# calc_bulk_withdrawals! — Calculate bulk water irrigation withdrawals
# Ported from: CalcBulkWithdrawals in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    calc_bulk_withdrawals!(irrig, waterfluxbulk, soilhydrology, soilstate,
                           col_data, patch_data,
                           mask_soilc, mask_soilp,
                           bounds_c, bounds_p, nlevsoi, dtime,
                           qflx_sfc_irrig_bulk_patch,
                           qflx_gw_demand_bulk_patch,
                           qflx_gw_demand_bulk_col)

Calculate irrigation withdrawals for bulk water.

Ported from `CalcBulkWithdrawals` in `IrrigationMod.F90`.
"""
function calc_bulk_withdrawals!(
    irrig::IrrigationData,
    waterfluxbulk::WaterFluxBulkData,
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    col_data::ColumnData,
    patch_data::PatchData,
    mask_soilc::BitVector,
    mask_soilp::BitVector,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    nlevsoi::Int,
    dtime::Float64,
    qflx_sfc_irrig_bulk_patch::Vector{Float64},
    qflx_gw_demand_bulk_patch::Vector{Float64},
    qflx_gw_demand_bulk_col::Vector{Float64}
)
    # Calculate per-patch surface irrigation and demand
    for p in bounds_p
        mask_soilp[p] || continue

        if irrig.n_irrig_steps_left_patch[p] > 0
            qflx_sfc_irrig_bulk_patch[p]    = irrig.sfc_irrig_rate_patch[p]
            irrig.qflx_irrig_demand_patch[p] = irrig.irrig_rate_demand_patch[p]
            qflx_gw_demand_bulk_patch[p]    = irrig.qflx_irrig_demand_patch[p] -
                                               qflx_sfc_irrig_bulk_patch[p]
            irrig.n_irrig_steps_left_patch[p] -= 1
        else
            qflx_sfc_irrig_bulk_patch[p]    = 0.0
            irrig.qflx_irrig_demand_patch[p] = 0.0
            qflx_gw_demand_bulk_patch[p]    = 0.0
        end
    end

    # Average patch surface irrigation to column
    p2c_irrig!(waterfluxbulk.wf.qflx_sfc_irrig_col, qflx_sfc_irrig_bulk_patch,
               patch_data, mask_soilc, bounds_c, bounds_p)

    # Average patch groundwater demand to column
    p2c_irrig!(qflx_gw_demand_bulk_col, qflx_gw_demand_bulk_patch,
               patch_data, mask_soilc, bounds_c, bounds_p)

    # If using groundwater irrigation, compute withdrawals
    if irrig.params.use_groundwater_irrigation
        calc_irrig_withdrawals!(
            soilhydrology, soilstate,
            qflx_gw_demand_bulk_col,
            waterfluxbulk.wf.qflx_gw_uncon_irrig_lyr_col,
            waterfluxbulk.wf.qflx_gw_con_irrig_col,
            col_data.nbedrock, col_data.z,
            mask_soilc, bounds_c, nlevsoi, dtime)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# calc_total_gw_uncon_irrig! — Sum unconfined GW irrigation by layer
# Ported from: CalcTotalGWUnconIrrig in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    calc_total_gw_uncon_irrig!(waterflux, mask_soilc, bounds_c, nlevsoi)

Calculate total irrigation withdrawal from unconfined aquifer (sum over layers).

Ported from `CalcTotalGWUnconIrrig` in `IrrigationMod.F90`.
"""
function calc_total_gw_uncon_irrig!(
    waterflux::WaterFluxData,
    mask_soilc::BitVector,
    bounds_c::UnitRange{Int},
    nlevsoi::Int
)
    for c in bounds_c
        mask_soilc[c] || continue
        waterflux.qflx_gw_uncon_irrig_col[c] = 0.0
    end
    for j in 1:nlevsoi
        for c in bounds_c
            mask_soilc[c] || continue
            waterflux.qflx_gw_uncon_irrig_col[c] +=
                waterflux.qflx_gw_uncon_irrig_lyr_col[c, j]
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# calc_application_fluxes! — Set drip/sprinkler application fluxes
# Ported from: CalcApplicationFluxes in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    calc_application_fluxes!(irrig, waterflux, is_bulk,
                             qflx_sfc_irrig_bulk_patch, qflx_sfc_irrig_bulk_col,
                             qflx_gw_demand_bulk_patch, qflx_gw_demand_bulk_col,
                             patch_data, mask_soilc, mask_soilp,
                             bounds_c, bounds_p)

Calculate irrigation application fluxes (drip / sprinkler).

Ported from `CalcApplicationFluxes` in `IrrigationMod.F90`.
"""
function calc_application_fluxes!(
    irrig::IrrigationData,
    waterflux::WaterFluxData,
    is_bulk::Bool,
    qflx_sfc_irrig_bulk_patch::Vector{Float64},
    qflx_sfc_irrig_bulk_col::Vector{Float64},
    qflx_gw_demand_bulk_patch::Vector{Float64},
    qflx_gw_demand_bulk_col::Vector{Float64},
    patch_data::PatchData,
    mask_soilc::BitVector,
    mask_soilp::BitVector,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int}
)
    # Compute total groundwater irrigation withdrawn per column
    qflx_gw_irrig_withdrawn_col = zeros(Float64, length(waterflux.qflx_gw_uncon_irrig_col))
    for c in bounds_c
        mask_soilc[c] || continue
        qflx_gw_irrig_withdrawn_col[c] =
            waterflux.qflx_gw_uncon_irrig_col[c] + waterflux.qflx_gw_con_irrig_col[c]
    end

    for p in bounds_p
        mask_soilp[p] || continue
        c = patch_data.column[p]

        # Compute surface irrigation for this patch
        if is_bulk
            qflx_sfc_irrig = qflx_sfc_irrig_bulk_patch[p]
        else
            if qflx_sfc_irrig_bulk_col[c] > 0.0
                qflx_sfc_irrig = waterflux.qflx_sfc_irrig_col[c] *
                    (qflx_sfc_irrig_bulk_patch[p] / qflx_sfc_irrig_bulk_col[c])
            else
                qflx_sfc_irrig = 0.0
            end
        end

        # Compute groundwater irrigation for this patch
        if qflx_gw_demand_bulk_col[c] > 0.0
            qflx_gw_irrig_applied = qflx_gw_irrig_withdrawn_col[c] *
                (qflx_gw_demand_bulk_patch[p] / qflx_gw_demand_bulk_col[c])
        else
            qflx_gw_irrig_applied = 0.0
        end

        qflx_irrig_tot = qflx_sfc_irrig + qflx_gw_irrig_applied

        # Set drip/sprinkler irrigation based on method
        waterflux.qflx_irrig_drip_patch[p]      = 0.0
        waterflux.qflx_irrig_sprinkler_patch[p] = 0.0

        if irrig.irrig_method_patch[p] == IRRIG_METHOD_DRIP
            waterflux.qflx_irrig_drip_patch[p] = qflx_irrig_tot
        elseif irrig.irrig_method_patch[p] == IRRIG_METHOD_SPRINKLER
            waterflux.qflx_irrig_sprinkler_patch[p] = qflx_irrig_tot
        else
            error("irrig_method_patch set to invalid value $(irrig.irrig_method_patch[p]) for patch $p")
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# calc_irrigation_fluxes! — Apply irrigation fluxes
# Ported from: CalcIrrigationFluxes in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    calc_irrigation_fluxes!(irrig, soilhydrology, soilstate,
                            waterfluxbulk, col_data, patch_data,
                            mask_soilc, mask_soilp,
                            bounds_c, bounds_p, nlevsoi, dtime)

Apply irrigation computed by `calc_irrigation_needed!` to set various fluxes.
Should be called once AND ONLY ONCE per time step.

Sets irrigation withdrawal and application fluxes in waterfluxbulk:
- qflx_sfc_irrig_col
- qflx_gw_uncon_irrig_lyr_col, qflx_gw_uncon_irrig_col
- qflx_gw_con_irrig_col
- qflx_irrig_drip_patch, qflx_irrig_sprinkler_patch

Ported from `CalcIrrigationFluxes` in `IrrigationMod.F90`.
"""
function calc_irrigation_fluxes!(
    irrig::IrrigationData,
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    waterfluxbulk::WaterFluxBulkData,
    col_data::ColumnData,
    patch_data::PatchData,
    mask_soilc::BitVector,
    mask_soilp::BitVector,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    nlevsoi::Int,
    dtime::Float64
)
    np = length(bounds_p)
    nc = length(bounds_c)

    qflx_sfc_irrig_bulk_patch  = zeros(Float64, length(irrig.sfc_irrig_rate_patch))
    qflx_gw_demand_bulk_patch  = zeros(Float64, length(irrig.sfc_irrig_rate_patch))
    qflx_gw_demand_bulk_col   = zeros(Float64, length(waterfluxbulk.wf.qflx_sfc_irrig_col))

    # Calculate bulk withdrawals
    calc_bulk_withdrawals!(irrig, waterfluxbulk, soilhydrology, soilstate,
                           col_data, patch_data, mask_soilc, mask_soilp,
                           bounds_c, bounds_p, nlevsoi, dtime,
                           qflx_sfc_irrig_bulk_patch,
                           qflx_gw_demand_bulk_patch,
                           qflx_gw_demand_bulk_col)

    # Note: Tracer withdrawals are skipped in this port (tracers not yet ported)

    # Sum up unconfined GW irrigation if using groundwater
    if irrig.params.use_groundwater_irrigation
        calc_total_gw_uncon_irrig!(waterfluxbulk.wf, mask_soilc, bounds_c, nlevsoi)
    end

    # Calculate application fluxes for bulk water
    calc_application_fluxes!(irrig, waterfluxbulk.wf, true,
                             qflx_sfc_irrig_bulk_patch,
                             waterfluxbulk.wf.qflx_sfc_irrig_col,
                             qflx_gw_demand_bulk_patch,
                             qflx_gw_demand_bulk_col,
                             patch_data, mask_soilc, mask_soilp,
                             bounds_c, bounds_p)

    return nothing
end

# ---------------------------------------------------------------------------
# use_groundwater_irrigation — Accessor
# Ported from: UseGroundwaterIrrigation in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    use_groundwater_irrigation(irrig) -> Bool

Returns true if groundwater irrigation is enabled.

Ported from `UseGroundwaterIrrigation` in `IrrigationMod.F90`.
"""
function use_groundwater_irrigation(irrig::IrrigationData)
    return irrig.params.use_groundwater_irrigation
end

# ---------------------------------------------------------------------------
# check_namelist_validity! — Validate parameters
# Ported from: CheckNamelistValidity in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    check_namelist_validity!(params; use_aquifer_layer=false)

Check for validity of irrigation parameters.

Ported from `CheckNamelistValidity` in `IrrigationMod.F90`.
"""
function check_namelist_validity!(params::IrrigationParamsData; use_aquifer_layer::Bool=false)
    if params.irrig_min_lai < 0.0
        error("irrig_min_lai must be >= 0, got $(params.irrig_min_lai)")
    end
    if params.irrig_start_time < 0 || params.irrig_start_time >= ISECSPDAY
        error("irrig_start_time must be >= 0 and < $ISECSPDAY, got $(params.irrig_start_time)")
    end
    if params.irrig_length <= 0 || params.irrig_length > ISECSPDAY
        error("irrig_length must be > 0 and <= $ISECSPDAY, got $(params.irrig_length)")
    end
    if params.irrig_target_smp >= 0.0
        error("irrig_target_smp must be negative, got $(params.irrig_target_smp)")
    end
    if params.irrig_target_smp < WILTING_POINT_SMP
        error("irrig_target_smp must be >= WILTING_POINT_SMP ($WILTING_POINT_SMP), got $(params.irrig_target_smp)")
    end
    if params.irrig_depth < 0.0
        error("irrig_depth must be > 0, got $(params.irrig_depth)")
    end
    if params.irrig_threshold_fraction < 0.0 || params.irrig_threshold_fraction > 1.0
        error("irrig_threshold_fraction must be between 0 and 1, got $(params.irrig_threshold_fraction)")
    end
    if params.limit_irrigation_if_rof_enabled
        if params.irrig_river_volume_threshold < 0.0 || params.irrig_river_volume_threshold > 1.0
            error("irrig_river_volume_threshold must be between 0 and 1, got $(params.irrig_river_volume_threshold)")
        end
    end
    if params.use_groundwater_irrigation && !params.limit_irrigation_if_rof_enabled
        error("use_groundwater_irrigation only makes sense if limit_irrigation_if_rof_enabled is set")
    end
    if use_aquifer_layer && params.use_groundwater_irrigation
        error("use_groundwater_irrigation and use_aquifer_layer may not be used simultaneously")
    end
    return nothing
end
