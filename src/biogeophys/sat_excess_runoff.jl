# ==========================================================================
# Ported from: src/biogeophys/SaturatedExcessRunoffMod.F90
# Calculates surface runoff due to saturated surface excess.
#
# Also computes fsat (fraction of each column that is saturated), using
# either the TOPModel-based or VIC-based parameterization.
#
# Public functions:
#   saturated_excess_runoff_init!   -- Initialize SaturatedExcessRunoffData
#   saturated_excess_runoff!        -- Main science routine
#   compute_fsat_topmodel!          -- TOPModel fsat computation
#   compute_fsat_vic!               -- VIC fsat computation
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level parameter type and constants
# ---------------------------------------------------------------------------

"""
    SaturatedExcessRunoffParams

Module-level parameters for saturated excess runoff, read from the params file
in Fortran via `readParams`. In Julia, set directly.

Ported from `params_type` in `SaturatedExcessRunoffMod.F90`.
"""
Base.@kwdef mutable struct SaturatedExcessRunoffParams
    fff::Float64 = NaN   # Decay factor for fractional saturated area (1/m)
end

# Global parameter instance (mirrors Fortran params_inst)
const sat_excess_runoff_params = SaturatedExcessRunoffParams()

# Fsat method constants
const FSAT_METHOD_TOPMODEL = 1
const FSAT_METHOD_VIC      = 2

# ---------------------------------------------------------------------------
# Data type
# ---------------------------------------------------------------------------

"""
    SaturatedExcessRunoffData

Saturated excess runoff data structure. Holds fsat and fcov at the column
level, plus the fsat_method flag.

Ported from `saturated_excess_runoff_type` in `SaturatedExcessRunoffMod.F90`.
"""
Base.@kwdef mutable struct SaturatedExcessRunoffData{FT<:Real}
    fsat_col::Vector{FT} = Float64[]    # col fractional area with water table at surface
    fcov_col::Vector{FT} = Float64[]    # col fractional impermeable area
    fsat_method::Int = FSAT_METHOD_TOPMODEL  # method selector
end

# ---------------------------------------------------------------------------
# Initialization
# ---------------------------------------------------------------------------

"""
    saturated_excess_runoff_init!(ser::SaturatedExcessRunoffData, nc::Int;
                                  use_vichydro::Bool=false)

Allocate and initialize a `SaturatedExcessRunoffData` instance for `nc`
columns. Sets `fsat_method` based on the `use_vichydro` flag.

Combines `Init`, `InitAllocate`, `InitHistory`, and `InitCold` from the
Fortran source.
"""
function saturated_excess_runoff_init!(ser::SaturatedExcessRunoffData, nc::Int;
                                       use_vichydro::Bool=false)
    # InitAllocate: allocate with NaN (matching Fortran)
    ser.fsat_col = fill(NaN, nc)
    ser.fcov_col = fill(NaN, nc)

    # InitCold: set fsat_method based on use_vichydro
    if use_vichydro
        ser.fsat_method = FSAT_METHOD_VIC
    else
        ser.fsat_method = FSAT_METHOD_TOPMODEL
    end

    return nothing
end

"""
    saturated_excess_runoff_clean!(ser::SaturatedExcessRunoffData)

Deallocate (reset to empty) all fields.
"""
function saturated_excess_runoff_clean!(ser::SaturatedExcessRunoffData)
    ser.fsat_col = Float64[]
    ser.fcov_col = Float64[]
    return nothing
end

# ---------------------------------------------------------------------------
# Science routines
# ---------------------------------------------------------------------------

"""
    compute_fsat_topmodel!(mask_hydrology::BitVector, bounds_col::UnitRange{Int},
                           frost_table::Vector{<:Real},
                           zwt::Vector{<:Real},
                           zwt_perched::Vector{<:Real},
                           wtfact::Vector{<:Real},
                           fff::Float64,
                           fsat::Vector{<:Real})

Compute fsat using the TOPModel-based parameterization (CLM default).

If the frost table is between the perched water table and the main water
table, uses the perched water table depth; otherwise uses the main water
table depth.

Ported from `ComputeFsatTopmodel` in `SaturatedExcessRunoffMod.F90`.
"""
function compute_fsat_topmodel!(mask_hydrology::BitVector,
                                 bounds_col::UnitRange{Int},
                                 frost_table::Vector{<:Real},
                                 zwt::Vector{<:Real},
                                 zwt_perched::Vector{<:Real},
                                 wtfact::Vector{<:Real},
                                 fff::Real,
                                 fsat::Vector{<:Real})
    for c in bounds_col
        mask_hydrology[c] || continue

        if frost_table[c] > zwt_perched[c] && frost_table[c] <= zwt[c]
            # use perched water table to determine fsat (if present)
            fsat[c] = wtfact[c] * exp(-0.5 * fff * zwt_perched[c])
        else
            fsat[c] = wtfact[c] * exp(-0.5 * fff * zwt[c])
        end
    end

    return nothing
end

"""
    compute_fsat_vic!(mask_hydrology::BitVector, bounds_col::UnitRange{Int},
                      b_infil::Vector{<:Real},
                      top_max_moist::Vector{<:Real},
                      top_moist_limited::Vector{<:Real},
                      fsat::Vector{<:Real})

Compute fsat using the VIC-based parameterization.

Citation: Wood et al. 1992, "A land-surface hydrology parameterization with
subgrid variability for general circulation models", JGR 97(D3), 2717-2728.

fsat is equivalent to A in VIC papers.

Ported from `ComputeFsatVic` in `SaturatedExcessRunoffMod.F90`.
"""
function compute_fsat_vic!(mask_hydrology::BitVector,
                            bounds_col::UnitRange{Int},
                            b_infil::Vector{<:Real},
                            top_max_moist::Vector{<:Real},
                            top_moist_limited::Vector{<:Real},
                            fsat::Vector{<:Real})
    for c in bounds_col
        mask_hydrology[c] || continue

        ex = b_infil[c] / (1.0 + b_infil[c])
        # fsat is equivalent to A in VIC papers
        fsat[c] = 1.0 - (1.0 - top_moist_limited[c] / top_max_moist[c])^ex
    end

    return nothing
end

"""
    saturated_excess_runoff!(ser::SaturatedExcessRunoffData,
                              mask_hydrology::BitVector,
                              bounds_col::UnitRange{Int},
                              col::ColumnData,
                              lun::LandunitData,
                              soilhydrology_inst::SoilHydrologyData,
                              soilstate_inst::SoilStateData,
                              waterfluxbulk_inst::WaterFluxBulkData;
                              crop_fsat_equals_zero::Bool=false,
                              hillslope_fsat_equals_zero::Bool=false)

Calculate surface runoff due to saturated surface excess.

Sets `ser.fsat_col`, `ser.fcov_col`, and
`waterfluxbulk_inst.qflx_sat_excess_surf_col`.

Steps:
1. Compute fsat via TOPModel or VIC method
2. Optionally zero fsat for crop columns
3. Optionally zero fsat for upland hillslope columns
4. Compute qflx_sat_excess_surf = fsat * qflx_rain_plus_snomelt
5. Set fcov = fsat for history output
6. Add flood water flux to runoff for urban columns

Ported from `SaturatedExcessRunoff` in `SaturatedExcessRunoffMod.F90`.
"""
function saturated_excess_runoff!(ser::SaturatedExcessRunoffData,
                                   mask_hydrology::BitVector,
                                   bounds_col::UnitRange{Int},
                                   col::ColumnData,
                                   lun::LandunitData,
                                   soilhydrology_inst::SoilHydrologyData,
                                   soilstate_inst::SoilStateData,
                                   waterfluxbulk_inst::WaterFluxBulkData;
                                   crop_fsat_equals_zero::Bool=false,
                                   hillslope_fsat_equals_zero::Bool=false)

    # Aliases (matching Fortran associate block)
    fsat                   = ser.fsat_col
    fcov                   = ser.fcov_col
    qflx_sat_excess_surf   = waterfluxbulk_inst.qflx_sat_excess_surf_col
    qflx_floodc            = waterfluxbulk_inst.wf.qflx_floodc_col
    qflx_rain_plus_snomelt = waterfluxbulk_inst.wf.qflx_rain_plus_snomelt_col

    # -------------------------------------------------------------------
    # 1. Compute fsat
    # -------------------------------------------------------------------
    if ser.fsat_method == FSAT_METHOD_TOPMODEL
        compute_fsat_topmodel!(mask_hydrology, bounds_col,
            soilhydrology_inst.frost_table_col,
            soilhydrology_inst.zwt_col,
            soilhydrology_inst.zwt_perched_col,
            soilstate_inst.wtfact_col,
            sat_excess_runoff_params.fff,
            fsat)
    elseif ser.fsat_method == FSAT_METHOD_VIC
        compute_fsat_vic!(mask_hydrology, bounds_col,
            soilhydrology_inst.b_infil_col,
            soilhydrology_inst.top_max_moist_col,
            soilhydrology_inst.top_moist_limited_col,
            fsat)
    else
        error("saturated_excess_runoff!: Unrecognized fsat_method: $(ser.fsat_method)")
    end

    # -------------------------------------------------------------------
    # 2. Set fsat to zero for crop columns
    # -------------------------------------------------------------------
    if crop_fsat_equals_zero
        for c in bounds_col
            mask_hydrology[c] || continue
            l = col.landunit[c]
            if lun.itype[l] == ISTCROP
                fsat[c] = 0.0
            end
        end
    end

    # -------------------------------------------------------------------
    # 3. Set fsat to zero for upland hillslope columns
    # -------------------------------------------------------------------
    if hillslope_fsat_equals_zero
        for c in bounds_col
            mask_hydrology[c] || continue
            if col.is_hillslope_column[c] && col.active[c]
                # Set fsat to zero for upland columns (cold != ISPVAL)
                if col.cold[c] != ISPVAL
                    fsat[c] = 0.0
                end
            end
        end
    end

    # -------------------------------------------------------------------
    # 4. Compute qflx_sat_excess_surf and set fcov
    # -------------------------------------------------------------------
    for c in bounds_col
        mask_hydrology[c] || continue
        # only send fast runoff directly to streams
        qflx_sat_excess_surf[c] = fsat[c] * qflx_rain_plus_snomelt[c]
        # Set fcov just to have it on the history file
        fcov[c] = fsat[c]
    end

    # -------------------------------------------------------------------
    # 5. For urban columns, send flood water flux to runoff
    # -------------------------------------------------------------------
    for c in bounds_col
        mask_hydrology[c] || continue
        if col.urbpoi[c]
            # send flood water flux to runoff for all urban columns
            qflx_sat_excess_surf[c] = qflx_sat_excess_surf[c] + qflx_floodc[c]
        end
    end

    return nothing
end
