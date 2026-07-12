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
Base.@kwdef mutable struct IrrigationParamsData{FT<:Real}
    # Minimum LAI for irrigation
    irrig_min_lai::FT = 0.0
    # Time of day to check whether we need irrigation, seconds (0 = midnight)
    irrig_start_time::Int = 21600
    # Desired amount of time to irrigate per day (sec)
    irrig_length::Int = 14400
    # Target soil matric potential for irrigation (mm)
    irrig_target_smp::FT = -3400.0
    # Soil depth to which we measure for irrigation (m)
    irrig_depth::FT = 0.6
    # Threshold fraction for irrigation trigger (0 = wilting point, 1 = target)
    irrig_threshold_fraction::FT = 1.0
    # Threshold for river water volume below which irrigation is shut off (fraction)
    irrig_river_volume_threshold::FT = 0.1
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
Base.@kwdef mutable struct IrrigationData{FT<:Real,
                              V<:AbstractVector{FT},
                              M<:AbstractMatrix{FT},
                              VI<:AbstractVector{<:Integer}}
    # Parameters
    params::IrrigationParamsData = IrrigationParamsData()
    # Land model time step (sec)
    dtime::Int = 0
    # Number of irrigation time steps per day
    irrig_nsteps_per_day::Int = 0
    # Relative saturation at wilting point [col, nlevsoi]
    relsat_wilting_point_col::M = Matrix{Float64}(undef, 0, 0)
    # Relative saturation at irrigation target [col, nlevsoi]
    relsat_target_col::M = Matrix{Float64}(undef, 0, 0)
    # Patch irrigation application method [patch]
    irrig_method_patch::VI = Int[]
    # Current irrigation rate from surface water [mm/s] [patch]
    sfc_irrig_rate_patch::V = Float64[]
    # Current irrigation rate demand, neglecting surface water source limitation [mm/s] [patch]
    irrig_rate_demand_patch::V = Float64[]
    # Number of time steps for which we still need to irrigate today [patch]
    n_irrig_steps_left_patch::VI = Int[]
    # Irrigation flux neglecting surface water source limitation [mm/s] [patch]
    qflx_irrig_demand_patch::V = Float64[]
end

IrrigationData{FT}(; kwargs...) where {FT<:Real} =
    IrrigationData{FT, Vector{FT}, Matrix{FT}, Vector{Int}}(; kwargs...)
Adapt.@adapt_structure IrrigationData


# ---------------------------------------------------------------------------
# irrigation_init_allocate! — Allocate irrigation data arrays
# Ported from: IrrigationInitAllocate in IrrigationMod.F90
# ---------------------------------------------------------------------------

"""
    irrigation_init_allocate!(irrig, np, nc, nlevsoi)

Allocate and initialize all irrigation data arrays.

Ported from `IrrigationInitAllocate` in `IrrigationMod.F90`.
"""
function irrigation_init_allocate!(irrig::IrrigationData{FT}, np::Int, nc::Int, nlevsoi::Int) where {FT}
    irrig.qflx_irrig_demand_patch    = fill(FT(NaN), np)
    irrig.relsat_wilting_point_col   = fill(FT(NaN), nc, nlevsoi)
    irrig.relsat_target_col          = fill(FT(NaN), nc, nlevsoi)
    irrig.sfc_irrig_rate_patch       = fill(FT(NaN), np)
    irrig.irrig_rate_demand_patch    = fill(FT(NaN), np)
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

Initialize irrigation history fields (seeds the demand field to SPVAL).

The per-type "register my history fields" body is a no-op: history I/O IS ported
(`src/infrastructure/history_io.jl`), but CLM.jl declares history fields in a
central registry rather than per type.

Ported from `IrrigationInitHistory` in `IrrigationMod.F90`.
"""
function irrigation_init_history!(irrig::IrrigationData{FT}, np::Int) where {FT}
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
    pftcon_irrigated::Vector{<:Real},
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
    pftcon_irrigated::Vector{<:Real},
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
    relsat_wilting_point::Matrix{<:Real},
    relsat_target::Matrix{<:Real},
    pftcon_irrigated::Vector{<:Real},
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
function irrigation_clean!(irrig::IrrigationData{FT}) where {FT}
    irrig.qflx_irrig_demand_patch    = FT[]
    irrig.relsat_wilting_point_col   = Matrix{FT}(undef, 0, 0)
    irrig.relsat_target_col          = Matrix{FT}(undef, 0, 0)
    irrig.sfc_irrig_rate_patch       = FT[]
    irrig.irrig_rate_demand_patch    = FT[]
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

Handle restart of irrigation variables. No-op: restart I/O IS ported
(`src/infrastructure/restart_io.jl`), but restart variables are declared in a
central registry rather than per type. The irrigation counters are not currently
in that registry.

Ported from `Restart` in `IrrigationMod.F90`.
"""
function irrigation_restart!(irrig::IrrigationData{FT}) where {FT}
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
    pftcon_irrigated::Vector{<:Real},
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
function relsat_to_h2osoi(relsat::Real, eff_porosity::Real, dz::Real)
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
    elai::Real,
    londeg::Real,
    pftcon_irrigated::Vector{<:Real},
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
# Zero the masked columns (one thread per column), then atomic-scatter each patch's
# area-weighted contribution into its column (one thread per patch, `_scatter_add!`
# handles the many-patch→one-column race on the device; plain += on the KA CPU path).
@kernel function _p2c_irrig_zero_kernel!(col_out, @Const(mask_soilc), cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_soilc[c]
        col_out[c] = zero(eltype(col_out))
    end
end
@kernel function _p2c_irrig_scatter_kernel!(col_out, @Const(patch_in), @Const(pcol),
        @Const(wtcol), @Const(mask_soilc), cmin::Int, cmax::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        c = pcol[p]
        if cmin <= c <= cmax && mask_soilc[c]
            _scatter_add!(col_out, c, patch_in[p] * wtcol[p])
        end
    end
end

function p2c_irrig!(
    col_out::AbstractVector{<:Real},
    patch_in::AbstractVector{<:Real},
    patch_data::PatchData,
    mask_soilc::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int}
)
    FT = eltype(col_out)
    pcol  = _to_backend_like(col_out, FT, patch_data.column)
    wtcol = _to_backend_like(col_out, FT, patch_data.wtcol)
    _launch!(_p2c_irrig_zero_kernel!, col_out, mask_soilc, first(bounds_c), last(bounds_c))
    _launch!(_p2c_irrig_scatter_kernel!, col_out, patch_in, pcol, wtcol, mask_soilc,
             first(bounds_c), last(bounds_c), first(bounds_p), last(bounds_p);
             ndrange = length(patch_in))
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
# Zero the gridcells (one thread per gridcell), then atomic-scatter each column's
# area-weighted contribution into its gridcell (one thread per column).
@kernel function _c2g_irrig_zero_kernel!(garr, gmin::Int, gmax::Int)
    g = @index(Global)
    @inbounds if gmin <= g <= gmax
        garr[g] = zero(eltype(garr))
    end
end
@kernel function _c2g_irrig_scatter_kernel!(garr, @Const(carr), @Const(cgrid),
        @Const(wtgcell), gmin::Int, gmax::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax
        g = cgrid[c]
        if gmin <= g <= gmax
            _scatter_add!(garr, g, carr[c] * wtgcell[c])
        end
    end
end

function c2g_irrig!(
    garr::AbstractVector{<:Real},
    carr::AbstractVector{<:Real},
    col_data::ColumnData,
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int}
)
    FT = eltype(garr)
    cgrid   = _to_backend_like(garr, FT, col_data.gridcell)
    wtgcell = _to_backend_like(garr, FT, col_data.wtgcell)
    _launch!(_c2g_irrig_zero_kernel!, garr, first(bounds_g), last(bounds_g))
    _launch!(_c2g_irrig_scatter_kernel!, garr, carr, cgrid, wtgcell,
             first(bounds_g), last(bounds_g), first(bounds_c), last(bounds_c);
             ndrange = length(carr))
    return nothing
end

# ---------------------------------------------------------------------------
# calc_deficit_volr_limited! — Limit deficit by river volume
# Ported from: CalcDeficitVolrLimited in IrrigationMod.F90
# ---------------------------------------------------------------------------

# Per-gridcell deficit-limited ratio from available river volume (independent
# per gridcell). Gated on [gmin, gmax] so entries outside bounds_g keep their
# initialized value (1.0).
@kernel function _irrig_volr_ratio_kernel!(ratio_grc, @Const(volr), @Const(area),
                                           river_volume_threshold, @Const(deficit_grc_arg),
                                           m3_over_km2_to_mm, gmin::Int, gmax::Int)
    g = @index(Global)
    @inbounds if gmin <= g <= gmax
        T = eltype(ratio_grc)
        if volr[g] > zero(T)
            available_volr = volr[g] * (one(T) - river_volume_threshold)
            max_deficit_supported_by_volr = available_volr / area[g] * m3_over_km2_to_mm
        else
            max_deficit_supported_by_volr = zero(T)
        end
        ratio_grc[g] = deficit_grc_arg[g] > max_deficit_supported_by_volr ?
            max_deficit_supported_by_volr / deficit_grc_arg[g] : one(T)
    end
end

# river_volume_threshold is converted to the array element type by the caller so no
# Float64 reaches the kernel (Metal rejects it); M3_OVER_KM2_TO_MM likewise.
irrig_volr_ratio!(ratio_grc, volr, area, river_volume_threshold, deficit_grc_arg,
                  gmin::Int, gmax::Int) =
    _launch!(_irrig_volr_ratio_kernel!, ratio_grc, volr, area, river_volume_threshold,
             deficit_grc_arg, eltype(ratio_grc)(M3_OVER_KM2_TO_MM), gmin, gmax)

# Per-column deficit scaled by its gridcell's limited ratio (independent per
# column). Gated on the check_for_irrig mask; unmasked columns keep their
# pre-zeroed value.
@kernel function _irrig_deficit_limited_kernel!(deficit_volr_limited, @Const(mask),
                                                @Const(deficit), @Const(ratio_grc),
                                                @Const(gridcell), cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        deficit_volr_limited[c] = deficit[c] * ratio_grc[gridcell[c]]
    end
end

irrig_deficit_limited!(deficit_volr_limited, mask, deficit, ratio_grc, gridcell,
                       cmin::Int, cmax::Int) =
    _launch!(_irrig_deficit_limited_kernel!, deficit_volr_limited, mask, deficit,
             ratio_grc, gridcell, cmin, cmax)

"""
    calc_deficit_volr_limited!(irrig, deficit, volr, deficit_volr_limited,
                               check_for_irrig_col, col_data, grc,
                               bounds_c, bounds_g)

Calculate deficit limited by river volume for each column.

Ported from `CalcDeficitVolrLimited` in `IrrigationMod.F90`.
"""
function calc_deficit_volr_limited!(
    irrig::IrrigationData,
    deficit::AbstractVector{<:Real},
    volr::AbstractVector{<:Real},
    deficit_volr_limited::AbstractVector{<:Real},
    check_for_irrig_col::AbstractVector{Bool},
    col_data::ColumnData,
    grc::GridcellData,
    bounds_c::UnitRange{Int},
    bounds_g::UnitRange{Int}
)
    FT = eltype(deficit)
    deficit_grc = fill!(similar(deficit, FT, length(volr)), zero(FT))
    deficit_limited_ratio_grc = fill!(similar(deficit, FT, length(volr)), one(FT))

    # Average deficit to gridcell level (unity scaling)
    c2g_irrig!(deficit_grc, deficit, col_data, bounds_c, bounds_g)

    irrig_volr_ratio!(deficit_limited_ratio_grc, volr, grc.area,
                      FT(irrig.params.irrig_river_volume_threshold), deficit_grc,
                      first(bounds_g), last(bounds_g))

    fill!(view(deficit_volr_limited, bounds_c), zero(FT))
    irrig_deficit_limited!(deficit_volr_limited, check_for_irrig_col, deficit,
                           deficit_limited_ratio_grc, col_data.gridcell,
                           first(bounds_c), last(bounds_c))

    return nothing
end

# ---------------------------------------------------------------------------
# calc_irrigation_needed! — Determine irrigation need
# Ported from: CalcIrrigationNeeded in IrrigationMod.F90
# ---------------------------------------------------------------------------

# (1) Per-patch: inline point_needs_check_for_irrig, set the patch flag, and OR it
# into the column flag (multiple patches writing `true` to one column race-safely —
# an idempotent store, no read-modify-write). cfip/cfic must be pre-filled false.
@kernel function _irrig_need_patchcheck_kernel!(cfip, cfic, @Const(ivt), @Const(elai),
        @Const(pftcon_irrigated), @Const(local_time_sec), @Const(pcol), @Const(mask_ev),
        irrig_min_lai, irrig_start_time::Int, dtime_i::Int, isecspday::Int,
        pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_ev[p]
        need = false
        if pftcon_irrigated[ivt[p]] == one(eltype(pftcon_irrigated)) && elai[p] > irrig_min_lai
            start_offset = irrig_start_time - dtime_i
            seconds_since = mod(local_time_sec[p] - start_offset, isecspday)
            if seconds_since < 0
                seconds_since += isecspday
            end
            if seconds_since < dtime_i
                need = true
            end
        end
        cfip[p] = need
        if need
            cfic[pcol[p]] = true
        end
    end
end

# (2) Per-column: walk the soil layers with the reached_max_depth state machine,
# summing liquid / target / wilting-point water into LOCAL scalars (relsat_to_h2osoi
# inlined with T(DENH2O) so no Float64 reaches the kernel), then store once.
@kernel function _irrig_need_layer_kernel!(liq_tot, target_tot, wp_tot, @Const(cfic),
        @Const(z), @Const(nbedrock), @Const(t_soisno), @Const(h2osoi_liq),
        @Const(relsat_target), @Const(relsat_wp), @Const(eff_porosity), @Const(dz),
        irrig_depth, tfrz, denh2o, cmin::Int, cmax::Int, nlevsoi::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && cfic[c]
        T = eltype(liq_tot)
        lt = zero(T); tt = zero(T); wt = zero(T)
        reached = false
        for j in 1:nlevsoi
            if !reached
                if z[c, j] > irrig_depth
                    reached = true
                elseif j > nbedrock[c]
                    reached = true
                elseif t_soisno[c, j] <= tfrz
                    reached = true
                else
                    lt += h2osoi_liq[c, j]
                    tt += eff_porosity[c, j] * relsat_target[c, j] * denh2o * dz[c, j]
                    wt += eff_porosity[c, j] * relsat_wp[c, j] * denh2o * dz[c, j]
                end
            end
        end
        liq_tot[c] = lt; target_tot[c] = tt; wp_tot[c] = wt
    end
end

# (3) Per-column deficit (the host loop's deficit<0 error() is dropped: in this branch
# liq_tot < threshold <= target_tot, so target_tot - liq_tot > 0 by construction).
@kernel function _irrig_need_deficit_kernel!(deficit, @Const(cfic), @Const(liq_tot),
        @Const(target_tot), @Const(wp_tot), threshold_fraction, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && cfic[c]
        T = eltype(deficit)
        thresh = wp_tot[c] + threshold_fraction * (target_tot[c] - wp_tot[c])
        deficit[c] = liq_tot[c] < thresh ? (target_tot[c] - liq_tot[c]) : zero(T)
    end
end

# (4) Per-patch: convert the column deficit to an irrigation rate for checked patches.
@kernel function _irrig_need_rate_kernel!(sfc_rate, rate_demand, nsteps_left,
        @Const(cfip), @Const(deficit_vl), @Const(deficit), @Const(pcol), @Const(mask_ev),
        denom, nsteps::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_ev[p] && cfip[p]
        c = pcol[p]
        sfc_rate[p]    = deficit_vl[c] / denom
        rate_demand[p] = deficit[c] / denom
        nsteps_left[p] = nsteps
    end
end

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
    elai::AbstractVector{<:Real},
    t_soisno::AbstractMatrix{<:Real},
    eff_porosity::AbstractMatrix{<:Real},
    h2osoi_liq::AbstractMatrix{<:Real},
    volr::AbstractVector{<:Real},
    rof_prognostic::Bool,
    pftcon_irrigated::AbstractVector{<:Real},
    local_time_sec_patch::AbstractVector{<:Integer},
    col_data::ColumnData,
    grc::GridcellData,
    patch_data::PatchData,
    mask_exposedveg::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    bounds_g::UnitRange{Int},
    nlevsoi::Int
)
    FT   = eltype(h2osoi_liq)
    maxc = last(bounds_c)
    np   = length(elai)

    # Backend-matched temporaries (host Array on the CPU path, device array on the GPU).
    cfip = fill!(similar(elai, Bool, np), false)     # check_for_irrig_patch
    cfic = fill!(similar(elai, Bool, maxc), false)   # check_for_irrig_col
    liq_tot    = fill!(similar(elai, FT, maxc), zero(FT))
    target_tot = fill!(similar(elai, FT, maxc), zero(FT))
    wp_tot     = fill!(similar(elai, FT, maxc), zero(FT))
    deficit    = fill!(similar(elai, FT, maxc), zero(FT))
    deficit_vl = fill!(similar(elai, FT, maxc), zero(FT))

    # Index/topology arrays onto the state backend (Int-preserving where relevant).
    pcol = _to_backend_like(liq_tot, FT, patch_data.column)
    ivt  = _to_backend_like(liq_tot, FT, patch_data.itype)
    lts  = _to_backend_like(liq_tot, FT, local_time_sec_patch)
    pfti = _to_backend_like(liq_tot, FT, pftcon_irrigated)
    nbed = _to_backend_like(liq_tot, FT, col_data.nbedrock)

    pmin, pmax = first(bounds_p), last(bounds_p)
    cmin, cmax = first(bounds_c), last(bounds_c)

    # (1) which patches/columns need checking
    _launch!(_irrig_need_patchcheck_kernel!, cfip, cfic, ivt, elai, pfti, lts, pcol,
             mask_exposedveg, FT(irrig.params.irrig_min_lai),
             Int(irrig.params.irrig_start_time), Int(irrig.dtime), Int(ISECSPDAY),
             pmin, pmax)

    # (2) soil-water accumulation over layers
    _launch!(_irrig_need_layer_kernel!, liq_tot, target_tot, wp_tot, cfic,
             col_data.z, nbed, t_soisno, h2osoi_liq,
             irrig.relsat_target_col, irrig.relsat_wilting_point_col, eff_porosity,
             col_data.dz, FT(irrig.params.irrig_depth), FT(TFRZ), FT(DENH2O),
             cmin, cmax, nlevsoi)

    # (3) deficits
    _launch!(_irrig_need_deficit_kernel!, deficit, cfic, liq_tot, target_tot, wp_tot,
             FT(irrig.params.irrig_threshold_fraction), cmin, cmax)

    # (4) limit by available river volume if enabled, else copy through
    if irrig.params.limit_irrigation_if_rof_enabled && rof_prognostic
        calc_deficit_volr_limited!(irrig, deficit, volr, deficit_vl,
                                   cfic, col_data, grc, bounds_c, bounds_g)
    else
        copyto!(view(deficit_vl, bounds_c), view(deficit, bounds_c))
    end

    # (5) convert to irrigation rate for checked patches
    denom = FT(irrig.dtime) * FT(irrig.irrig_nsteps_per_day)
    _launch!(_irrig_need_rate_kernel!, irrig.sfc_irrig_rate_patch,
             irrig.irrig_rate_demand_patch, irrig.n_irrig_steps_left_patch,
             cfip, deficit_vl, deficit, pcol, mask_exposedveg,
             denom, Int(irrig.irrig_nsteps_per_day), pmin, pmax)

    return nothing
end

# ---------------------------------------------------------------------------
# calc_bulk_withdrawals! — Calculate bulk water irrigation withdrawals
# Ported from: CalcBulkWithdrawals in IrrigationMod.F90
# ---------------------------------------------------------------------------

# One thread per patch: set surface-irrig / demand fluxes and decrement the
# remaining-steps counter (every write is to the patch's own index → race-free,
# byte-identical to the host loop on the KA CPU backend).
@kernel function _bulk_withdrawals_patch_kernel!(qflx_sfc, qflx_demand, qflx_gw_demand,
        nsteps_left, @Const(sfc_rate), @Const(rate_demand), @Const(mask_soilp),
        pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        T = eltype(qflx_sfc)
        if nsteps_left[p] > 0
            qflx_sfc[p]       = sfc_rate[p]
            qflx_demand[p]    = rate_demand[p]
            qflx_gw_demand[p] = qflx_demand[p] - qflx_sfc[p]
            nsteps_left[p]   -= 1
        else
            qflx_sfc[p]       = zero(T)
            qflx_demand[p]    = zero(T)
            qflx_gw_demand[p] = zero(T)
        end
    end
end

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
    mask_soilc::AbstractVector{Bool},
    mask_soilp::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    nlevsoi::Int,
    dtime::Real,
    qflx_sfc_irrig_bulk_patch::AbstractVector{<:Real},
    qflx_gw_demand_bulk_patch::AbstractVector{<:Real},
    qflx_gw_demand_bulk_col::AbstractVector{<:Real}
)
    # Calculate per-patch surface irrigation and demand (one thread per patch)
    _launch!(_bulk_withdrawals_patch_kernel!, qflx_sfc_irrig_bulk_patch,
             irrig.qflx_irrig_demand_patch, qflx_gw_demand_bulk_patch,
             irrig.n_irrig_steps_left_patch, irrig.sfc_irrig_rate_patch,
             irrig.irrig_rate_demand_patch, mask_soilp,
             first(bounds_p), last(bounds_p))

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
# One thread per column: sum the per-layer unconfined GW irrigation over nlevsoi into
# a LOCAL scalar (the fixed-index `col[c] += lyr[c,j]`-in-loop form miscompiles on the
# KA CPU backend under --check-bounds), then store once. Byte-identical (same j order).
@kernel function _total_gw_uncon_kernel!(tot, @Const(lyr), @Const(mask_soilc),
        cmin::Int, cmax::Int, nlevsoi::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_soilc[c]
        s = zero(eltype(tot))
        for j in 1:nlevsoi
            s += lyr[c, j]
        end
        tot[c] = s
    end
end

function calc_total_gw_uncon_irrig!(
    waterflux::WaterFluxData,
    mask_soilc::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    nlevsoi::Int
)
    _launch!(_total_gw_uncon_kernel!, waterflux.qflx_gw_uncon_irrig_col,
             waterflux.qflx_gw_uncon_irrig_lyr_col, mask_soilc,
             first(bounds_c), last(bounds_c), nlevsoi)
    return nothing
end

# ---------------------------------------------------------------------------
# calc_application_fluxes! — Set drip/sprinkler application fluxes
# Ported from: CalcApplicationFluxes in IrrigationMod.F90
# ---------------------------------------------------------------------------

# Per-column total groundwater irrigation withdrawn (unconfined + confined).
@kernel function _app_gw_withdrawn_kernel!(withdrawn, @Const(uncon), @Const(con),
        @Const(mask_soilc), cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_soilc[c]
        withdrawn[c] = uncon[c] + con[c]
    end
end

# Per-patch drip/sprinkler routing. `is_bulk`/method ids are isbits scalars. The
# invalid-method error() of the host loop is dropped — set_irrig_method! guarantees
# the method is DRIP or SPRINKLER, so it never fires (byte-identical on valid input).
@kernel function _app_fluxes_patch_kernel!(drip, sprinkler, @Const(withdrawn),
        @Const(sfc_irrig_col), @Const(sfc_bulk_patch), @Const(sfc_bulk_col),
        @Const(gw_demand_patch), @Const(gw_demand_col), @Const(irrig_method),
        @Const(pcol), @Const(mask_soilp), is_bulk::Bool, drip_id::Int,
        sprinkler_id::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask_soilp[p]
        T = eltype(drip)
        c = pcol[p]
        if is_bulk
            qflx_sfc = sfc_bulk_patch[p]
        else
            qflx_sfc = sfc_bulk_col[c] > zero(T) ?
                sfc_irrig_col[c] * (sfc_bulk_patch[p] / sfc_bulk_col[c]) : zero(T)
        end
        qflx_gw = gw_demand_col[c] > zero(T) ?
            withdrawn[c] * (gw_demand_patch[p] / gw_demand_col[c]) : zero(T)
        qflx_tot = qflx_sfc + qflx_gw
        drip[p]      = zero(T)
        sprinkler[p] = zero(T)
        if irrig_method[p] == drip_id
            drip[p] = qflx_tot
        elseif irrig_method[p] == sprinkler_id
            sprinkler[p] = qflx_tot
        end
    end
end

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
    qflx_sfc_irrig_bulk_patch::AbstractVector{<:Real},
    qflx_sfc_irrig_bulk_col::AbstractVector{<:Real},
    qflx_gw_demand_bulk_patch::AbstractVector{<:Real},
    qflx_gw_demand_bulk_col::AbstractVector{<:Real},
    patch_data::PatchData,
    mask_soilc::AbstractVector{Bool},
    mask_soilp::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int}
)
    # Compute total groundwater irrigation withdrawn per column (backend-matched temp)
    FT = eltype(waterflux.qflx_gw_uncon_irrig_col)
    qflx_gw_irrig_withdrawn_col = fill!(
        similar(waterflux.qflx_gw_uncon_irrig_col, FT,
                length(waterflux.qflx_gw_uncon_irrig_col)), zero(FT))
    _launch!(_app_gw_withdrawn_kernel!, qflx_gw_irrig_withdrawn_col,
             waterflux.qflx_gw_uncon_irrig_col, waterflux.qflx_gw_con_irrig_col,
             mask_soilc, first(bounds_c), last(bounds_c))

    pcol = _to_backend_like(qflx_gw_irrig_withdrawn_col, FT, patch_data.column)
    _launch!(_app_fluxes_patch_kernel!, waterflux.qflx_irrig_drip_patch,
             waterflux.qflx_irrig_sprinkler_patch, qflx_gw_irrig_withdrawn_col,
             waterflux.qflx_sfc_irrig_col, qflx_sfc_irrig_bulk_patch,
             qflx_sfc_irrig_bulk_col, qflx_gw_demand_bulk_patch,
             qflx_gw_demand_bulk_col, irrig.irrig_method_patch, pcol, mask_soilp,
             is_bulk, IRRIG_METHOD_DRIP, IRRIG_METHOD_SPRINKLER,
             first(bounds_p), last(bounds_p);
             ndrange = length(waterflux.qflx_irrig_drip_patch))

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
    mask_soilc::AbstractVector{Bool},
    mask_soilp::AbstractVector{Bool},
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    nlevsoi::Int,
    dtime::Real
)
    np = length(bounds_p)
    nc = length(bounds_c)

    FT = eltype(irrig.sfc_irrig_rate_patch)
    # Backend-matched temporaries (host Array on CPU, device array on the GPU) so the
    # kernels launched inside calc_bulk_withdrawals!/calc_application_fluxes! don't mix
    # host and device operands.
    qflx_sfc_irrig_bulk_patch = fill!(similar(irrig.sfc_irrig_rate_patch, FT,
                                              length(irrig.sfc_irrig_rate_patch)), zero(FT))
    qflx_gw_demand_bulk_patch = fill!(similar(irrig.sfc_irrig_rate_patch, FT,
                                              length(irrig.sfc_irrig_rate_patch)), zero(FT))
    qflx_gw_demand_bulk_col   = fill!(similar(waterfluxbulk.wf.qflx_sfc_irrig_col, FT,
                                    length(waterfluxbulk.wf.qflx_sfc_irrig_col)), zero(FT))

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
