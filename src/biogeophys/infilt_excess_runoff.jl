# ==========================================================================
# Ported from: src/biogeophys/InfiltrationExcessRunoffMod.F90
# Compute infiltration excess runoff and related variables
#
# Public functions:
#   infiltration_excess_runoff!    -- Calculate surface runoff due to
#                                     infiltration excess
#   compute_qinmax_hksat!          -- Compute qinmax from hksat
#                                     (CLM default parameterization)
#
# Public types:
#   InfiltrationExcessRunoffData   -- State container for this module
#   InfiltrationExcessRunoffParams -- Module-level parameters (e_ice)
# ==========================================================================

# ---- Module-level parameters ----

"""
    InfiltrationExcessRunoffParams

Parameters read from parameter file for infiltration excess runoff.

Ported from `params_type` in `InfiltrationExcessRunoffMod.F90`.
"""
Base.@kwdef mutable struct InfiltrationExcessRunoffParams
    e_ice::Float64 = 6.0   # Soil ice impedance factor (unitless)
end

const infilt_excess_params = InfiltrationExcessRunoffParams()

# ---- Module-level constants ----

# For methods that don't generate infiltration excess runoff, we specify a huge
# qinmax value -- an effectively infinite max infiltration rate.
# 1e200 mm H2O/s is large enough to be effectively infinite, while not so
# large as to cause floating point overflows elsewhere.
const QINMAX_UNLIMITED = 1.0e200   # mm H2O/s

# Method flags for computing qinmax
const QINMAX_METHOD_NONE  = 0
const QINMAX_METHOD_HKSAT = 1

# ---- Data type ----

"""
    InfiltrationExcessRunoffData

Infiltration excess runoff state data. Holds `qinmax_col` (maximum
infiltration rate) and `qinmax_method` (method selector).

Ported from `infiltration_excess_runoff_type` in
`InfiltrationExcessRunoffMod.F90`.
"""
Base.@kwdef mutable struct InfiltrationExcessRunoffData{FT<:Real,
                                                        V<:AbstractVector{FT}}
    qinmax_col::V = Float64[]                  # col maximum infiltration rate (mm H2O/s)
    qinmax_method::Int = QINMAX_METHOD_HKSAT   # method for computing qinmax
end

# Single-type-parameter convenience constructor (mirrors the other state structs):
# `InfiltrationExcessRunoffData{FT}()` defaults the storage to `Vector{FT}`, so the
# AD dual-copy path (`wrapper{D}()`) and Float32 construction keep working.
InfiltrationExcessRunoffData{FT}(; kwargs...) where {FT<:Real} =
    InfiltrationExcessRunoffData{FT, Vector{FT}}(; kwargs...)

# Device movement: relocate the single array field; qinmax_method is a plain Int.
Adapt.@adapt_structure InfiltrationExcessRunoffData

# ---- Init / Allocate / Cold ----

"""
    infilt_excess_runoff_init!(ier::InfiltrationExcessRunoffData, nc::Int;
                                use_vichydro::Bool=false)

Allocate and initialize an `InfiltrationExcessRunoffData` instance for `nc`
columns. Sets `qinmax_method` based on `use_vichydro`.

Ported from `Init`, `InitAllocate`, and `InitCold` in
`InfiltrationExcessRunoffMod.F90`.
"""
function infilt_excess_runoff_init!(ier::InfiltrationExcessRunoffData{FT}, nc::Int;
                                     use_vichydro::Bool = false) where {FT}
    # InitAllocate
    ier.qinmax_col = fill(FT(NaN), nc)

    # InitCold
    if use_vichydro
        ier.qinmax_method = QINMAX_METHOD_NONE
    else
        ier.qinmax_method = QINMAX_METHOD_HKSAT
    end

    return nothing
end

"""
    infilt_excess_runoff_clean!(ier::InfiltrationExcessRunoffData)

Deallocate (reset to empty) all fields of an `InfiltrationExcessRunoffData`
instance.
"""
function infilt_excess_runoff_clean!(ier::InfiltrationExcessRunoffData)
    ier.qinmax_col = Float64[]
    return nothing
end

# =========================================================================
# Science routines
# =========================================================================

"""
    infiltration_excess_runoff!(ier, soilhydrology, soilstate,
                                 fsat_col, waterfluxbulk,
                                 waterdiagnosticbulk,
                                 mask_hydrology, bounds;
                                 params)

Calculate surface runoff due to infiltration excess.

Sets `waterfluxbulk.qflx_infl_excess_col` and `ier.qinmax_col`. These
are valid within the hydrology filter. Both give averages over the entire
column. However, qinmax is implicitly 0 over the fraction of the column
given by fsat, and qflx_infl_excess is implicitly 0 over both fsat and
frac_h2osfc.

# Arguments
- `ier::InfiltrationExcessRunoffData`  : infiltration excess state (output: qinmax_col)
- `soilhydrology::SoilHydrologyData`  : soil hydrology state (input: icefrac_col)
- `soilstate::SoilStateData`          : soil state (input: hksat_col)
- `fsat_col::Vector{<:Real}`         : fractional area with water table at surface (from SaturatedExcessRunoff)
- `waterfluxbulk::WaterFluxBulkData`  : water fluxes (input: qflx_in_soil_col; output: qflx_infl_excess_col)
- `waterdiagnosticbulk::WaterDiagnosticBulkData` : water diagnostics (input: frac_h2osfc_col)
- `mask_hydrology::BitVector`          : column-level hydrology mask
- `bounds::UnitRange{Int}`             : column bounds
- `params::InfiltrationExcessRunoffParams` : module parameters (keyword, default = infilt_excess_params)

Ported from `InfiltrationExcessRunoff` in `InfiltrationExcessRunoffMod.F90`.
"""
# ---- infiltration_excess_runoff! : per-column qinmax NONE-fill kernel ----
# Fill qinmax_on_unsaturated_area with the (effectively infinite) unlimited value
# for the QINMAX_METHOD_NONE path. One thread per column, mask-gated.
@kernel function _infexcess_qinmax_none_kernel!(qinmax_unsat, @Const(mask), qinmax_unlimited)
    c = @index(Global)
    @inbounds if mask[c]
        qinmax_unsat[c] = qinmax_unlimited
    end
end

# ---- infiltration_excess_runoff! : per-column column-average + excess kernel ----
# qinmax[c]           = (1 - fsat[c]) * qinmax_on_unsaturated_area[c]
# qflx_infl_excess[c] = max(0, qflx_in_soil[c] - (1 - frac_h2osfc[c]) * qinmax[c])
# Each column is fully independent (no loop-carried deps). Literals converted to
# the working eltype so no Float64 reaches a Float32-only backend (Metal).
@kernel function _infexcess_colavg_kernel!(qinmax, qflx_infl_excess, @Const(mask),
                                           @Const(fsat), @Const(qinmax_unsat),
                                           @Const(qflx_in_soil), @Const(frac_h2osfc))
    c = @index(Global)
    @inbounds if mask[c]
        one_ft = one(eltype(qinmax))
        zero_ft = zero(eltype(qinmax))
        qinmax[c] = (one_ft - fsat[c]) * qinmax_unsat[c]
        qflx_infl_excess[c] = max(zero_ft,
            qflx_in_soil[c] - (one_ft - frac_h2osfc[c]) * qinmax[c])
    end
end

function infiltration_excess_runoff!(
    ier::InfiltrationExcessRunoffData,
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    fsat_col::AbstractVector{<:Real},
    waterfluxbulk::WaterFluxBulkData,
    waterdiagnosticbulk::WaterDiagnosticBulkData,
    mask_hydrology,
    bounds::UnitRange{Int};
    params::InfiltrationExcessRunoffParams = infilt_excess_params
)
    # Aliases matching Fortran associate block
    qinmax           = ier.qinmax_col
    fsat             = fsat_col
    qflx_infl_excess = waterfluxbulk.qflx_infl_excess_col
    qflx_in_soil     = waterfluxbulk.qflx_in_soil_col
    frac_h2osfc      = waterdiagnosticbulk.frac_h2osfc_col

    # Validate method on host before launching (String/Int config resolved on host).
    if ier.qinmax_method != QINMAX_METHOD_NONE && ier.qinmax_method != QINMAX_METHOD_HKSAT
        error("InfiltrationExcessRunoff: Unrecognized qinmax_method: $(ier.qinmax_method)")
    end

    # Temporary workspace for qinmax on unsaturated area — device-resident, same
    # backend/eltype as the output qinmax (similar(), not zeros()).
    FT_ier = eltype(qinmax)
    qinmax_on_unsaturated_area = similar(qinmax, last(bounds))

    # --- Compute qinmax on unsaturated area ---
    if ier.qinmax_method == QINMAX_METHOD_NONE
        _launch!(_infexcess_qinmax_none_kernel!, qinmax_on_unsaturated_area,
                 mask_hydrology, FT_ier(QINMAX_UNLIMITED))
    else  # QINMAX_METHOD_HKSAT
        compute_qinmax_hksat!(soilhydrology, soilstate,
                              qinmax_on_unsaturated_area,
                              mask_hydrology, bounds;
                              params = params)
    end

    # --- Compute column-averaged qinmax and infiltration excess ---
    _launch!(_infexcess_colavg_kernel!, qinmax, qflx_infl_excess, mask_hydrology,
             fsat, qinmax_on_unsaturated_area, qflx_in_soil, frac_h2osfc)

    return nothing
end

"""
    compute_qinmax_hksat!(soilhydrology, soilstate,
                           qinmax_on_unsaturated_area,
                           mask_hydrology, bounds;
                           params)

Compute qinmax using the CLM default parameterization based on hksat.

For each column, computes:
    qinmax = min over layers 1:3 of
        10^(-e_ice * icefrac(c,j)) * hksat(c,j)

# Arguments
- `soilhydrology::SoilHydrologyData`    : input icefrac_col (ncols, nlevgrnd)
- `soilstate::SoilStateData`            : input hksat_col (ncols, nlevgrnd)
- `qinmax_on_unsaturated_area::Vector{<:Real}` : output
- `mask_hydrology::BitVector`            : column-level hydrology mask
- `bounds::UnitRange{Int}`               : column bounds
- `params::InfiltrationExcessRunoffParams` : module parameters (keyword)

Ported from `ComputeQinmaxHksat` in `InfiltrationExcessRunoffMod.F90`.
"""
# ---- compute_qinmax_hksat! : per-column min over top-3 layers kernel ----
# qinmax_unsat[c] = min over j in 1:3 of 10^(-e_ice * icefrac[c,j]) * hksat[c,j].
# One thread per column with an internal SEQUENTIAL j loop (a small min reduction;
# columns are independent). e_ice is passed as a working-eltype scalar and the
# base 10 / Inf literals are eltype-converted so no Float64 reaches a Float32-only
# backend (Metal); on Float64 these conversions are byte-identical.
@kernel function _qinmax_hksat_kernel!(qinmax_unsat, @Const(mask),
                                       @Const(icefrac), @Const(hksat), e_ice)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(qinmax_unsat)
        ten = T(10)
        qmin = T(Inf)
        for j in 1:3
            val = ten^(-e_ice * icefrac[c, j]) * hksat[c, j]
            if val < qmin
                qmin = val
            end
        end
        qinmax_unsat[c] = qmin
    end
end

function compute_qinmax_hksat!(
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    qinmax_on_unsaturated_area::AbstractVector{<:Real},
    mask_hydrology,
    bounds::UnitRange{Int};
    params::InfiltrationExcessRunoffParams = infilt_excess_params
)
    icefrac = soilhydrology.icefrac_col
    hksat   = soilstate.hksat_col

    FT = eltype(qinmax_on_unsaturated_area)
    _launch!(_qinmax_hksat_kernel!, qinmax_on_unsaturated_area, mask_hydrology,
             icefrac, hksat, FT(params.e_ice); ndrange = last(bounds))

    return nothing
end
