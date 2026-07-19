# ==========================================================================
# Predefined calibration parameter sets for CLM.jl
#
# Each function returns a Vector{CalibrationParameter} with sensible
# defaults, bounds, and transforms for common calibration targets.
#
# Apply functions now set CalibrationOverrides fields so that Dual-typed
# parameter values flow through the computation graph during AD.
# For array-based params (baseflow_scalar, fff, ksat_scale), we set both
# the override AND the arrays, since these are used directly in hydrology
# without an override injection point.
# ==========================================================================

# --- Apply functions using CalibrationOverrides ---

function _apply_baseflow_scalar!(inst, val)
    # The override IS the live mechanism: clm_run! reads
    # `overrides.baseflow_scalar` into `init_soil_hydrology_config`, which sets
    # the `BASEFLOW_SCALAR[]` Ref that soil_hydrology.jl:3002 actually uses.
    # There was also a loop over `inst.column.baseflow_scalar` here -- a field
    # ColumnData does not have -- so this function threw FieldError on its FIRST
    # call and the parameter had never once been applied. It survived because
    # test_calibration.jl asserted only `.name` and never invoked `apply!`.
    inst.overrides.baseflow_scalar = val
    return nothing
end

function _apply_fff!(inst, val)
    # Same defect and same fix as _apply_baseflow_scalar!: the override is live
    # (clm_run! copies it into `sat_excess_runoff_params.fff`), while the loop
    # over the non-existent `inst.column.fff` made every call throw.
    inst.overrides.fff = val
    return nothing
end

function _apply_medlynslope!(inst, val)
    inst.overrides.medlyn_slope = val
    return nothing
end

function _apply_csoilc!(inst, val)
    inst.overrides.csoilc = val
    return nothing
end

function _apply_vcmax25_scale!(inst, val)
    inst.overrides.vcmax25_scale = val
    return nothing
end

function _apply_jmax25_scale!(inst, val)
    inst.overrides.jmax25top_sf = val
    return nothing
end

function _apply_ksat_scale!(inst, val)
    inst.overrides.ksat_scale = val
    if hasproperty(inst, :soilhydrology)
        for idx in eachindex(inst.soilhydrology.hksat_col)
            inst.soilhydrology.hksat_col[idx] *= val
        end
    end
    return nothing
end

# --- Parameter set constructors ---

"""
    default_hydrology_params()

Hydrology calibration parameters: baseflow scalar, surface runoff decay factor.
"""
function default_hydrology_params()
    return CalibrationParameter[
        CalibrationParameter(
            "baseflow_scalar", 0.001, _apply_baseflow_scalar!,
            (1e-6, 0.1), :log),
        CalibrationParameter(
            "fff", 0.5, _apply_fff!,
            (0.01, 10.0), :log),
    ]
end

"""
    default_canopy_params()

Canopy/surface exchange calibration parameters.
"""
function default_canopy_params()
    return CalibrationParameter[
        CalibrationParameter(
            "medlynslope", 6.0, _apply_medlynslope!,
            (1.0, 20.0), :identity),
        CalibrationParameter(
            "csoilc", 0.004, _apply_csoilc!,
            (0.001, 0.02), :log),
    ]
end

"""
    default_photosynthesis_params()

Photosynthesis calibration parameters: Vcmax and Jmax scaling.
"""
function default_photosynthesis_params()
    return CalibrationParameter[
        CalibrationParameter(
            "vcmax25_scale", 1.0, _apply_vcmax25_scale!,
            (0.3, 3.0), :identity),
        CalibrationParameter(
            "jmax25_scale", 1.0, _apply_jmax25_scale!,
            (0.3, 3.0), :identity),
    ]
end

"""
    default_soil_params()

Soil thermal/hydraulic calibration parameters.
"""
function default_soil_params()
    return CalibrationParameter[
        CalibrationParameter(
            "ksat_scale", 1.0, _apply_ksat_scale!,
            (0.1, 10.0), :log),
    ]
end

"""
    all_default_params()

Combine all default parameter sets into one vector.
"""
function all_default_params()
    return vcat(
        default_hydrology_params(),
        default_canopy_params(),
        default_photosynthesis_params(),
        default_soil_params(),
    )
end

"""
    default_theta(params::Vector{CalibrationParameter})

Return the default θ vector for a set of calibration parameters.
"""
function default_theta(params::Vector{CalibrationParameter})
    return [default_theta(p) for p in params]
end
