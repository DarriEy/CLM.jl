# ==========================================================================
# Ported from: src/biogeophys/DaylengthMod.F90
# Computes daylength, max daylength, and manages daylength updates
# ==========================================================================

# --- Constants local to this module ---
const SECS_PER_RADIAN = 13750.9871          # seconds per radian of hour-angle
const LAT_EPSILON     = 10.0 * eps(1.0)     # epsilon for latitudes near pole
const POLE            = RPI / 2.0           # π/2
const OFFSET_POLE     = POLE - LAT_EPSILON  # pole offset to avoid cos(lat)<0

"""
    daylength(lat, decl)

Compute daylength in seconds for a single grid cell.

`lat` and `decl` are in radians.  Returns `NaN` when preconditions are
violated (|lat| ≥ pole + ε, or |decl| ≥ pole).

Ported from elemental function `daylength` in `DaylengthMod.F90`.
"""
function daylength(lat::Real, decl::Real)::Real
    # lat must be less than π/2 within a small tolerance
    if abs(lat) >= (POLE + LAT_EPSILON)
        return NaN
    end

    # decl must be strictly less than π/2
    if abs(decl) >= POLE
        return NaN
    end

    # Ensure latitude isn't too close to pole
    my_lat = min(OFFSET_POLE, max(-OFFSET_POLE, lat))

    temp = -(sin(my_lat) * sin(decl)) / (cos(my_lat) * cos(decl))
    temp = min(1.0, max(-1.0, temp))
    return 2.0 * SECS_PER_RADIAN * acos(temp)
end

"""
    compute_max_daylength!(lat, obliquity, max_daylength, bounds)

Compute maximum daylength for each grid cell.

Ported from subroutine `ComputeMaxDaylength` in `DaylengthMod.F90`.
"""
function compute_max_daylength!(lat::Vector{<:Real},
                                obliquity::Real,
                                max_daylength::Vector{<:Real},
                                bounds::UnitRange{Int})
    for g in bounds
        max_decl = obliquity
        if lat[g] < 0.0
            max_decl = -max_decl
        end
        max_daylength[g] = daylength(lat[g], max_decl)
    end
    return nothing
end

"""
    init_daylength!(grc::GridcellData, declin, declinm1, obliquity, bounds)

Initialize daylength, previous daylength, and max daylength for all grid cells.

Should be called with `declin` set to the value for the first model time step,
and `declinm1` to the value for the previous time step.

Ported from subroutine `InitDaylength` in `DaylengthMod.F90`.
"""
function init_daylength!(grc::GridcellData,
                         declin::Real,
                         declinm1::Real,
                         obliquity::Real,
                         bounds::UnitRange{Int})
    for g in bounds
        grc.prev_dayl[g] = daylength(grc.lat[g], declinm1)
        grc.dayl[g]      = daylength(grc.lat[g], declin)
    end

    compute_max_daylength!(grc.lat, obliquity, grc.max_dayl, bounds)

    return nothing
end

"""
    update_daylength!(grc::GridcellData, declin, obliquity, is_first_step, bounds)

Update daylength, previous daylength, and max daylength for all grid cells.

Should be called exactly once per time step.  On the first step of a run
segment (`is_first_step == true`), dayl and prev_dayl are left as-is (they
were set during initialization).

Ported from subroutine `UpdateDaylength` in `DaylengthMod.F90`.
"""
function update_daylength!(grc::GridcellData,
                           declin::Real,
                           obliquity::Real,
                           is_first_step::Bool,
                           bounds::UnitRange{Int})
    if !is_first_step
        for g in bounds
            grc.prev_dayl[g] = grc.dayl[g]
            grc.dayl[g]      = daylength(grc.lat[g], declin)
        end
    end

    compute_max_daylength!(grc.lat, obliquity, grc.max_dayl, bounds)

    return nothing
end
