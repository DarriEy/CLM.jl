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
function daylength(lat::Real, decl::Real)
    # Work in the promoted element type so the body carries no Float64 literals
    # on a Float32-only backend (Metal); identity on Float64 (byte-identical).
    T = promote_type(typeof(lat), typeof(decl))
    lt   = T(lat)
    dcl  = T(decl)
    pole = T(POLE)
    eps_lat       = T(LAT_EPSILON)
    offset_pole   = T(OFFSET_POLE)
    secs_per_rad  = T(SECS_PER_RADIAN)

    # lat must be less than π/2 within a small tolerance
    if abs(lt) >= (pole + eps_lat)
        return T(NaN)
    end

    # decl must be strictly less than π/2
    if abs(dcl) >= pole
        return T(NaN)
    end

    # Ensure latitude isn't too close to pole
    my_lat = min(offset_pole, max(-offset_pole, lt))

    temp = -(sin(my_lat) * sin(dcl)) / (cos(my_lat) * cos(dcl))
    temp = min(one(T), max(-one(T), temp))
    return T(2) * secs_per_rad * acos(temp)
end

# --------------------------------------------------------------------------
# Kernel: maximum daylength per grid cell. `daylength` is eltype-generic and
# device-callable; `obliquity` is passed pre-converted to the working eltype so
# no Float64 reaches a Float32-only backend (Metal). One thread per grid cell;
# `lo` is the first index of `bounds` so threads 1:length(bounds) map to g.
# --------------------------------------------------------------------------
@kernel function _max_daylength_kernel!(max_daylength, @Const(lat), obliquity, lo::Int)
    i = @index(Global)
    @inbounds begin
        g = lo + i - 1
        T = eltype(max_daylength)
        max_decl = obliquity
        if lat[g] < zero(eltype(lat))
            max_decl = -max_decl
        end
        max_daylength[g] = daylength(lat[g], max_decl)
    end
end

"""
    compute_max_daylength!(lat, obliquity, max_daylength, bounds)

Compute maximum daylength for each grid cell. Backend-agnostic (CPU loop or
GPU kernel); one thread per grid cell.

Ported from subroutine `ComputeMaxDaylength` in `DaylengthMod.F90`.
"""
function compute_max_daylength!(lat::AbstractVector{<:Real},
                                obliquity::Real,
                                max_daylength::AbstractVector{<:Real},
                                bounds::UnitRange{Int})
    isempty(bounds) && return nothing
    obl = convert(eltype(max_daylength), obliquity)
    _launch!(_max_daylength_kernel!, max_daylength, lat, obl, first(bounds);
             ndrange = length(bounds))
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

# --------------------------------------------------------------------------
# Kernel: advance daylength one step. Each thread reads its OWN old dayl into
# prev_dayl, then overwrites dayl — fully per-grid-cell independent. `declin` is
# pre-converted to the working eltype. `lo` = first(bounds).
# --------------------------------------------------------------------------
@kernel function _update_dayl_kernel!(dayl, prev_dayl, @Const(lat), declin, lo::Int)
    i = @index(Global)
    @inbounds begin
        g = lo + i - 1
        prev_dayl[g] = dayl[g]
        dayl[g]      = daylength(lat[g], declin)
    end
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
    if !is_first_step && !isempty(bounds)
        decl = convert(eltype(grc.dayl), declin)
        _launch!(_update_dayl_kernel!, grc.dayl, grc.prev_dayl, grc.lat, decl,
                 first(bounds); ndrange = length(bounds))
    end

    compute_max_daylength!(grc.lat, obliquity, grc.max_dayl, bounds)

    return nothing
end
