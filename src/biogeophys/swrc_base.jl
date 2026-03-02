# ==========================================================================
# Ported from: src/biogeophys/SoilWaterRetentionCurveMod.F90
# Abstract base type for soil water retention curve functions.
# ==========================================================================

"""
    SoilWaterRetentionCurve

Abstract base type for soil water retention curve implementations.

Concrete subtypes must implement:
- [`soil_hk!`](@ref)
- [`soil_suction!`](@ref)
- [`soil_suction_inverse!`](@ref)

Ported from `soil_water_retention_curve_type` in `SoilWaterRetentionCurveMod.F90`.
"""
abstract type SoilWaterRetentionCurve end

"""
    soil_hk!(swrc::SoilWaterRetentionCurve, c, j, s, imped, soilstate;
             dhkds=nothing) -> (hk, dhkds_out)

Compute hydraulic conductivity.

# Arguments
- `swrc`: soil water retention curve instance
- `c::Int`: column index
- `j::Int`: level index
- `s::Float64`: relative saturation, [0, 1]
- `imped::Float64`: ice impedance
- `soilstate::SoilStateData`: soil state data

# Returns
- `hk::Float64`: hydraulic conductivity [mm/s]
- `dhkds_out::Union{Float64, Nothing}`: d[hk]/ds [mm/s], or `nothing` if not requested

Ported from `soil_hk_interface` in `SoilWaterRetentionCurveMod.F90`.
"""
function soil_hk!(swrc::SoilWaterRetentionCurve, c::Int, j::Int,
                  s::Float64, imped::Float64, soilstate::SoilStateData)
    error("soil_hk! not implemented for $(typeof(swrc))")
end

"""
    soil_suction!(swrc::SoilWaterRetentionCurve, c, j, s, soilstate) -> (smp, dsmpds_out)

Compute soil suction potential.

# Arguments
- `swrc`: soil water retention curve instance
- `c::Int`: column index
- `j::Int`: level index
- `s::Float64`: relative saturation, [0, 1]
- `soilstate::SoilStateData`: soil state data

# Returns
- `smp::Float64`: soil suction, negative [mm]
- `dsmpds::Union{Float64, Nothing}`: d[smp]/ds [mm], or `nothing` if not requested

Ported from `soil_suction_interface` in `SoilWaterRetentionCurveMod.F90`.
"""
function soil_suction!(swrc::SoilWaterRetentionCurve, c::Int, j::Int,
                       s::Float64, soilstate::SoilStateData)
    error("soil_suction! not implemented for $(typeof(swrc))")
end

"""
    soil_suction_inverse!(swrc::SoilWaterRetentionCurve, c, j, smp_target, soilstate) -> s_target

Compute relative saturation at which soil suction is equal to a target value.
This is done by inverting the soil_suction equation to solve for s.

# Arguments
- `swrc`: soil water retention curve instance
- `c::Int`: column index
- `j::Int`: level index
- `smp_target::Float64`: target soil suction, negative [mm]
- `soilstate::SoilStateData`: soil state data

# Returns
- `s_target::Float64`: relative saturation at which smp = smp_target [0, 1]

Ported from `soil_suction_inverse_interface` in `SoilWaterRetentionCurveMod.F90`.
"""
function soil_suction_inverse!(swrc::SoilWaterRetentionCurve, c::Int, j::Int,
                               smp_target::Float64, soilstate::SoilStateData)
    error("soil_suction_inverse! not implemented for $(typeof(swrc))")
end
