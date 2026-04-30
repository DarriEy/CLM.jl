# ==========================================================================
# CalibrationOverrides — lightweight parameter override struct for AD
#
# Carried by CLMInstances. NaN = "use default from const params".
# When calibrating, apply!() sets real (or Dual) values which physics
# functions pick up via isnan() checks at injection points.
# ==========================================================================

"""
    CalibrationOverrides{FT<:Real}

Lightweight struct holding calibration parameter overrides.
Fields default to `NaN`, meaning "use the default value from const params".
When a field is set to a non-NaN value, physics functions use it instead.

This allows Dual-typed parameter values to flow through the computation graph
during ForwardDiff AD, enabling true gradient computation.
"""
mutable struct CalibrationOverrides{FT<:Real}
    csoilc::FT            # canopy soil turbulence coeff
    jmax25top_sf::FT      # jmax25top scale factor
    vcmax25_scale::FT     # multiplicative scale on vcmax25top
    medlyn_slope::FT      # Medlyn stomatal slope override
    baseflow_scalar::FT   # baseflow scalar override
    fff::FT               # runoff decay factor
    ksat_scale::FT        # hksat multiplicative scale
    bsw_mult::FT          # Clapp-Hornberger b multiplier
    watsat_mult::FT       # porosity multiplier
    sucsat_mult::FT       # saturated suction multiplier

    # Inner constructor: all NaN
    function CalibrationOverrides{FT}() where {FT<:Real}
        new{FT}(FT(NaN), FT(NaN), FT(NaN), FT(NaN), FT(NaN), FT(NaN), FT(NaN),
                FT(NaN), FT(NaN), FT(NaN))
    end
end

# Convenience: CalibrationOverrides() defaults to Float64
CalibrationOverrides() = CalibrationOverrides{Float64}()

"""
    validate_overrides!(overrides::CalibrationOverrides)

Check that every field is either NaN (intentional sentinel = use default) or finite
(valid override value). Throws `ArgumentError` if any field is Inf or -Inf, which
would indicate numerical corruption rather than an intentional sentinel.
"""
function validate_overrides!(overrides::CalibrationOverrides)
    for fname in fieldnames(CalibrationOverrides)
        v = getfield(overrides, fname)
        if isinf(v)
            throw(ArgumentError(
                "CalibrationOverrides field `$fname` is $v — " *
                "expected NaN (use default) or a finite value"))
        end
    end
    return nothing
end

"""
    active_override_fields(overrides) -> Vector{Symbol}

Return names of fields that have been set to non-NaN (active) values.
Used by validate_overrides_post! to detect accidental NaN corruption.
"""
function active_override_fields(overrides::CalibrationOverrides)
    return [fname for fname in fieldnames(CalibrationOverrides)
            if !isnan(getfield(overrides, fname))]
end

"""
    validate_overrides_post!(overrides, active_fields)

Post-run validation: check that fields which were finite before the run
are still finite after. If a field that was set to a real value has become
NaN during the computation, this indicates numerical corruption.
"""
function validate_overrides_post!(overrides::CalibrationOverrides, active_fields::Vector{Symbol})
    for fname in active_fields
        v = getfield(overrides, fname)
        if isnan(v) || isinf(v)
            throw(ArgumentError(
                "CalibrationOverrides field `$fname` became $v during computation — " *
                "numerical corruption detected (was finite before run)"))
        end
    end
    return nothing
end
