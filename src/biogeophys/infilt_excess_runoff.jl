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
Base.@kwdef mutable struct InfiltrationExcessRunoffData
    qinmax_col::Vector{Float64} = Float64[]   # col maximum infiltration rate (mm H2O/s)
    qinmax_method::Int = QINMAX_METHOD_HKSAT  # method for computing qinmax
end

# ---- Init / Allocate / Cold ----

"""
    infilt_excess_runoff_init!(ier::InfiltrationExcessRunoffData, nc::Int;
                                use_vichydro::Bool=false)

Allocate and initialize an `InfiltrationExcessRunoffData` instance for `nc`
columns. Sets `qinmax_method` based on `use_vichydro`.

Ported from `Init`, `InitAllocate`, and `InitCold` in
`InfiltrationExcessRunoffMod.F90`.
"""
function infilt_excess_runoff_init!(ier::InfiltrationExcessRunoffData, nc::Int;
                                     use_vichydro::Bool = false)
    # InitAllocate
    ier.qinmax_col = fill(NaN, nc)

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
- `fsat_col::Vector{Float64}`         : fractional area with water table at surface (from SaturatedExcessRunoff)
- `waterfluxbulk::WaterFluxBulkData`  : water fluxes (input: qflx_in_soil_col; output: qflx_infl_excess_col)
- `waterdiagnosticbulk::WaterDiagnosticBulkData` : water diagnostics (input: frac_h2osfc_col)
- `mask_hydrology::BitVector`          : column-level hydrology mask
- `bounds::UnitRange{Int}`             : column bounds
- `params::InfiltrationExcessRunoffParams` : module parameters (keyword, default = infilt_excess_params)

Ported from `InfiltrationExcessRunoff` in `InfiltrationExcessRunoffMod.F90`.
"""
function infiltration_excess_runoff!(
    ier::InfiltrationExcessRunoffData,
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    fsat_col::Vector{Float64},
    waterfluxbulk::WaterFluxBulkData,
    waterdiagnosticbulk::WaterDiagnosticBulkData,
    mask_hydrology::BitVector,
    bounds::UnitRange{Int};
    params::InfiltrationExcessRunoffParams = infilt_excess_params
)
    # Aliases matching Fortran associate block
    qinmax           = ier.qinmax_col
    fsat             = fsat_col
    qflx_infl_excess = waterfluxbulk.qflx_infl_excess_col
    qflx_in_soil     = waterfluxbulk.qflx_in_soil_col
    frac_h2osfc      = waterdiagnosticbulk.frac_h2osfc_col

    # Temporary workspace for qinmax on unsaturated area
    nc = length(bounds)
    qinmax_on_unsaturated_area = Vector{Float64}(undef, last(bounds))

    # --- Compute qinmax on unsaturated area ---
    if ier.qinmax_method == QINMAX_METHOD_NONE
        for c in bounds
            mask_hydrology[c] || continue
            qinmax_on_unsaturated_area[c] = QINMAX_UNLIMITED
        end
    elseif ier.qinmax_method == QINMAX_METHOD_HKSAT
        compute_qinmax_hksat!(soilhydrology, soilstate,
                              qinmax_on_unsaturated_area,
                              mask_hydrology, bounds;
                              params = params)
    else
        error("InfiltrationExcessRunoff: Unrecognized qinmax_method: $(ier.qinmax_method)")
    end

    # --- Compute column-averaged qinmax and infiltration excess ---
    for c in bounds
        mask_hydrology[c] || continue
        qinmax[c] = (1.0 - fsat[c]) * qinmax_on_unsaturated_area[c]
        qflx_infl_excess[c] = max(0.0,
            qflx_in_soil[c] - (1.0 - frac_h2osfc[c]) * qinmax[c])
    end

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
- `qinmax_on_unsaturated_area::Vector{Float64}` : output
- `mask_hydrology::BitVector`            : column-level hydrology mask
- `bounds::UnitRange{Int}`               : column bounds
- `params::InfiltrationExcessRunoffParams` : module parameters (keyword)

Ported from `ComputeQinmaxHksat` in `InfiltrationExcessRunoffMod.F90`.
"""
function compute_qinmax_hksat!(
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    qinmax_on_unsaturated_area::Vector{Float64},
    mask_hydrology::BitVector,
    bounds::UnitRange{Int};
    params::InfiltrationExcessRunoffParams = infilt_excess_params
)
    icefrac = soilhydrology.icefrac_col
    hksat   = soilstate.hksat_col

    for c in bounds
        mask_hydrology[c] || continue
        # Fortran: minval(10^(-e_ice * icefrac(c,1:3)) * hksat(c,1:3))
        qmin = Inf
        for j in 1:3
            val = 10.0^(-params.e_ice * icefrac[c, j]) * hksat[c, j]
            if val < qmin
                qmin = val
            end
        end
        qinmax_on_unsaturated_area[c] = qmin
    end

    return nothing
end
