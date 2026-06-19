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
Base.@kwdef mutable struct SaturatedExcessRunoffData{FT<:Real,
                             V<:AbstractVector{FT}}
    fsat_col::V = Float64[]    # col fractional area with water table at surface
    fcov_col::V = Float64[]    # col fractional impermeable area
    fsat_method::Int = FSAT_METHOD_TOPMODEL  # method selector
end

# Convenience constructor: SaturatedExcessRunoffData{FT}() resolves the array
# type to Vector{FT} (mirrors the other Data structs; needed because @kwdef does
# not synthesize a single-parameter keyword constructor).
SaturatedExcessRunoffData{FT}(; kwargs...) where {FT<:Real} =
    SaturatedExcessRunoffData{FT, Vector{FT}}(; kwargs...)

# Device-movable: lets the whole struct be Adapt.adapt'd onto a GPU backend
# (fields become device arrays). The Int fsat_method is left untouched by adapt.
Adapt.@adapt_structure SaturatedExcessRunoffData

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
                           fff::Real,
                           fsat::Vector{<:Real})

Compute fsat using the TOPModel-based parameterization (CLM default). `fff` (the
TOPModel decay factor) is `::Real`, so a ForwardDiff `Dual` for it differentiates
through `fsat = wtfact·exp(-0.5·fff·zwt)` — fsat IS differentiable w.r.t. fff.

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

# ---------------------------------------------------------------------------
# GPU-kernelized whole-function path for saturated_excess_runoff!
#
# Every loop in this routine is fully independent per-column (no loop-carried
# deps, no inter-column coupling), so the entire routine fuses into ONE
# per-column kernel. fsat_method and the two host flags are resolved on the
# host and passed as scalars; the per-column arrays are grouped into immutable
# device-view bundles (Adapt.@adapt_structure) to stay under the Metal arg cap.
# Literals are converted to the working eltype (T(...)) so no Float64 reaches a
# Float32-only backend, while remaining byte-identical on Float64 CPU.
#
# These definitions live directly above saturated_excess_runoff! for
# merge-hygiene (disjoint line-region).
# ---------------------------------------------------------------------------

# Column-level hydrology / soilstate inputs the fsat computation reads.
Base.@kwdef struct _SERHydDV{V}
    frost_table_col::V
    zwt_col::V
    zwt_perched_col::V
    wtfact_col::V
    b_infil_col::V
    top_max_moist_col::V
    top_moist_limited_col::V
end
Adapt.@adapt_structure _SERHydDV

_ser_hyd_dv(sh, ss) = _SERHydDV(;
    frost_table_col       = sh.frost_table_col,
    zwt_col               = sh.zwt_col,
    zwt_perched_col       = sh.zwt_perched_col,
    wtfact_col            = ss.wtfact_col,
    b_infil_col           = sh.b_infil_col,
    top_max_moist_col     = sh.top_max_moist_col,
    top_moist_limited_col = sh.top_moist_limited_col)

# Column / landunit topology + flux arrays the branch logic and outputs touch.
Base.@kwdef struct _SERColDV{V,VB,VI}
    landunit::VI
    itype_lun::VI
    is_hillslope_column::VB
    active::VB
    cold::VI
    urbpoi::VB
    qflx_floodc_col::V
    qflx_rain_plus_snomelt_col::V
end
Adapt.@adapt_structure _SERColDV

_ser_col_dv(col, lun, wfb) = _SERColDV(;
    landunit                   = col.landunit,
    itype_lun                  = lun.itype,
    is_hillslope_column        = col.is_hillslope_column,
    active                     = col.active,
    cold                       = col.cold,
    urbpoi                     = col.urbpoi,
    qflx_floodc_col            = wfb.wf.qflx_floodc_col,
    qflx_rain_plus_snomelt_col = wfb.wf.qflx_rain_plus_snomelt_col)

# Whole-function per-column kernel. fsat_method/flags are Int/Bool host scalars.
@kernel function _sat_excess_runoff_kernel!(fsat, fcov, qflx_sat_excess_surf,
        @Const(mask), hyd, cdv, fff,
        fsat_method::Int, crop_zero::Bool, hillslope_zero::Bool,
        istcrop::Int, ispval::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = typeof(fff)

        # --- 1. Compute fsat (TOPModel or VIC) ---
        if fsat_method == FSAT_METHOD_TOPMODEL
            ft = hyd.frost_table_col[c]
            zw = hyd.zwt_col[c]
            zwp = hyd.zwt_perched_col[c]
            wt = hyd.wtfact_col[c]
            if ft > zwp && ft <= zw
                fs = wt * exp(-T(0.5) * fff * zwp)
            else
                fs = wt * exp(-T(0.5) * fff * zw)
            end
        else
            bi = hyd.b_infil_col[c]
            ex = bi / (one(T) + bi)
            fs = one(T) - (one(T) - hyd.top_moist_limited_col[c] / hyd.top_max_moist_col[c])^ex
        end

        # --- 2. Zero fsat for crop columns ---
        if crop_zero
            l = cdv.landunit[c]
            if cdv.itype_lun[l] == istcrop
                fs = zero(T)
            end
        end

        # --- 3. Zero fsat for upland hillslope columns ---
        if hillslope_zero
            if cdv.is_hillslope_column[c] && cdv.active[c]
                if cdv.cold[c] != ispval
                    fs = zero(T)
                end
            end
        end

        fsat[c] = fs

        # --- 4. qflx_sat_excess_surf and fcov ---
        q = fs * cdv.qflx_rain_plus_snomelt_col[c]
        fcov[c] = fs

        # --- 5. Urban flood water added to runoff ---
        if cdv.urbpoi[c]
            q = q + cdv.qflx_floodc_col[c]
        end
        qflx_sat_excess_surf[c] = q
    end
end

"""
    saturated_excess_runoff!(ser::SaturatedExcessRunoffData,
                              mask_hydrology::AbstractVector{Bool},
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
                                   mask_hydrology::AbstractVector{Bool},
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

    # Validate fsat_method on the host (preserves the scalar error path; the
    # device kernel only sees a valid method selector).
    if ser.fsat_method != FSAT_METHOD_TOPMODEL && ser.fsat_method != FSAT_METHOD_VIC
        error("saturated_excess_runoff!: Unrecognized fsat_method: $(ser.fsat_method)")
    end

    # Convert the host Float64 fff parameter to the working eltype so a
    # Float32-only backend (Metal) carries no Float64 (byte-identical on F64).
    fff = oftype(zero(eltype(fsat)), sat_excess_runoff_params.fff)

    hyd = _ser_hyd_dv(soilhydrology_inst, soilstate_inst)
    cdv = _ser_col_dv(col, lun, waterfluxbulk_inst)

    # Whole routine in one per-column kernel (all loops are independent).
    _launch!(_sat_excess_runoff_kernel!, fsat, fcov, qflx_sat_excess_surf,
        mask_hydrology, hyd, cdv, fff,
        ser.fsat_method, crop_fsat_equals_zero, hillslope_fsat_equals_zero,
        ISTCROP, ISPVAL; ndrange = length(mask_hydrology))

    return nothing
end
