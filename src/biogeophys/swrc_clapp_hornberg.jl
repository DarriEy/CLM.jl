# ==========================================================================
# Ported from: src/biogeophys/SoilWaterRetentionCurveClappHornberg1978Mod.F90
# Implementation of soil water retention curve using the Clapp-Hornberg 1978
# parameterizations.
# ==========================================================================

"""
    SoilWaterRetentionCurveClappHornberg1978 <: SoilWaterRetentionCurve

Concrete implementation of [`SoilWaterRetentionCurve`](@ref) using the
Clapp and Hornberger (1978) parameterizations for hydraulic conductivity
and soil suction.

Ported from `soil_water_retention_curve_clapp_hornberg_1978_type` in
`SoilWaterRetentionCurveClappHornberg1978Mod.F90`.
"""
struct SoilWaterRetentionCurveClappHornberg1978 <: SoilWaterRetentionCurve end

"""
    soil_hk!(swrc::SoilWaterRetentionCurveClappHornberg1978, c, j, s, imped, soilstate)

Compute hydraulic conductivity using the Clapp-Hornberg 1978 parameterization.

    hk = imped * hksat * s^(2*bsw + 3)
    dhkds = (2*bsw + 3) * hk / s

# Arguments
- `c::Int`: column index
- `j::Int`: level index
- `s::Float64`: relative saturation, [0, 1]
- `imped::Float64`: ice impedance
- `soilstate::SoilStateData`: soil state data (uses `hksat_col`, `bsw_col`)

# Returns
- `(hk, dhkds)`: hydraulic conductivity [mm/s] and its derivative w.r.t. s

Ported from `soil_hk` in `SoilWaterRetentionCurveClappHornberg1978Mod.F90`.
"""
function soil_hk!(swrc::SoilWaterRetentionCurveClappHornberg1978, c::Int, j::Int,
                  s::Float64, imped::Float64, soilstate::SoilStateData)
    hksat = soilstate.hksat_col[c, j]
    bsw   = soilstate.bsw_col[c, j]

    # compute hydraulic conductivity
    hk = imped * hksat * s^(2.0 * bsw + 3.0)

    # compute the derivative
    dhkds = (2.0 * bsw + 3.0) * hk / s

    return (hk, dhkds)
end

"""
    soil_suction!(swrc::SoilWaterRetentionCurveClappHornberg1978, c, j, s, soilstate)

Compute soil suction potential using the Clapp-Hornberg 1978 parameterization.

    smp = -sucsat * s^(-bsw)
    dsmpds = -bsw * smp / s

# Arguments
- `c::Int`: column index
- `j::Int`: level index
- `s::Float64`: relative saturation, [0, 1]
- `soilstate::SoilStateData`: soil state data (uses `bsw_col`, `sucsat_col`)

# Returns
- `(smp, dsmpds)`: soil suction [mm] (negative) and its derivative w.r.t. s

Ported from `soil_suction` in `SoilWaterRetentionCurveClappHornberg1978Mod.F90`.
"""
function soil_suction!(swrc::SoilWaterRetentionCurveClappHornberg1978, c::Int, j::Int,
                       s::Real, soilstate::SoilStateData)
    bsw    = soilstate.bsw_col[c, j]
    sucsat = soilstate.sucsat_col[c, j]

    # compute soil suction potential, negative
    smp = -sucsat * s^(-bsw)

    # compute derivative
    dsmpds = -bsw * smp / s

    return (smp, dsmpds)
end

"""
    soil_suction_inverse!(swrc::SoilWaterRetentionCurveClappHornberg1978, c, j, smp_target, soilstate)

Compute relative saturation at which soil suction is equal to a target value.
This is done by inverting the soil_suction equation to solve for s.

    s_target = (-smp_target / sucsat)^(-1/bsw)

# Arguments
- `c::Int`: column index
- `j::Int`: level index
- `smp_target::Float64`: target soil suction, negative [mm]
- `soilstate::SoilStateData`: soil state data (uses `bsw_col`, `sucsat_col`)

# Returns
- `s_target::Float64`: relative saturation at which smp = smp_target [0, 1]

Ported from `soil_suction_inverse` in `SoilWaterRetentionCurveClappHornberg1978Mod.F90`.
"""
function soil_suction_inverse!(swrc::SoilWaterRetentionCurveClappHornberg1978, c::Int, j::Int,
                               smp_target::Real, soilstate::SoilStateData)
    bsw    = soilstate.bsw_col[c, j]
    sucsat = soilstate.sucsat_col[c, j]

    s_target = (-smp_target / sucsat)^(-1.0 / bsw)

    return s_target
end
