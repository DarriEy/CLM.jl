# ==========================================================================
# Ported from: src/biogeophys/QSatMod.F90 (130 lines)
# Saturation vapor pressure and specific humidity calculations
# Uses Flatau et al. (1992) polynomial approximations
# ==========================================================================

# --- Polynomial coefficients for saturation vapor pressure over water (0-100°C) ---
const QSAT_A = (
    6.11213476,
    0.444007856,
    0.143064234e-1,
    0.264461437e-3,
    0.305903558e-5,
    0.196237241e-7,
    0.892344772e-10,
    -0.373208410e-12,
    0.209339997e-15,
)

# --- Coefficients for d(es)/d(T) over water ---
const QSAT_B = (
    0.444017302,
    0.286064092e-1,
    0.794683137e-3,
    0.121211669e-4,
    0.103354611e-6,
    0.404125005e-9,
    -0.788037859e-12,
    -0.114596802e-13,
    0.381294516e-16,
)

# --- Coefficients for saturation vapor pressure over ice (-75 to 0°C) ---
const QSAT_C = (
    6.11123516,
    0.503109514,
    0.188369801e-1,
    0.420547422e-3,
    0.614396778e-5,
    0.602780717e-7,
    0.387940929e-9,
    0.149436277e-11,
    0.262655803e-14,
)

# --- Coefficients for d(es)/d(T) over ice ---
const QSAT_D = (
    0.503277922,
    0.377289173e-1,
    0.126801703e-2,
    0.249468427e-4,
    0.313703411e-6,
    0.257180651e-8,
    0.133268878e-10,
    0.394116744e-13,
    0.498070196e-16,
)

"""
    qsat(T, p) -> (qs, es, dqsdT, desdT)

Compute saturation specific humidity, vapor pressure, and their derivatives.

# Arguments
- `T::Float64`: Temperature [K]
- `p::Float64`: Pressure [Pa]

# Returns
- `qs::Float64`: Saturation specific humidity [kg/kg]
- `es::Float64`: Saturation vapor pressure [Pa]
- `dqsdT::Float64`: d(qs)/d(T) [kg/kg/K]
- `desdT::Float64`: d(es)/d(T) [Pa/K]

Ported from QSat() in QSatMod.F90. Uses Flatau et al. (1992) polynomials.
"""
function qsat(T::Real, p::Real)
    # eltype-generic (R = working precision) so this lowers to valid Metal IR
    # under Float32; byte-identical to the Float64 literals on the CPU path.
    R = typeof(T)
    td = T - R(TFRZ)  # temperature in °C

    # Saturation vapor pressure: water polynomials (td >= 0) vs ice (td < 0).
    # smooth_ifelse is exact for Float64 (preserves Fortran parity bit-for-bit
    # in :auto mode) but smoothly blends the kink at the freezing point for
    # ForwardDiff.Dual, so the derivative is continuous across TFRZ for AD.
    es = smooth_ifelse(td,
        R(100.0) * _poly8(td, QSAT_A),    # over water
        R(100.0) * _poly8(td, QSAT_C))    # over ice
    desdT = smooth_ifelse(td,
        R(100.0) * _poly8(td, QSAT_B),    # over water
        R(100.0) * _poly8(td, QSAT_D))    # over ice

    # Specific humidity from vapor pressure
    # qs = 0.622 * es / (p - 0.378*es)
    # dqs/dT = 0.622 * desdT * p / (p - 0.378*es)^2   (quotient rule)
    vp = one(R) / (p - R(0.378) * es)
    vp1 = min(es * vp, one(R))
    qs = R(0.622) * vp1
    dqsdT = R(0.622) * p * desdT * vp * vp

    return (qs, es, dqsdT, desdT)
end

"""
    qsat_no_derivs(T, p) -> (qs, es)

Compute saturation specific humidity and vapor pressure without derivatives.
Slightly more efficient when derivatives are not needed.
"""
function qsat_no_derivs(T::Real, p::Real)
    R = typeof(T)
    td = T - R(TFRZ)

    # See qsat() above: exact branch for Float64, smooth blend across TFRZ for Dual.
    es = smooth_ifelse(td,
        R(100.0) * _poly8(td, QSAT_A),    # over water
        R(100.0) * _poly8(td, QSAT_C))    # over ice

    vp1 = min(es / (p - R(0.378) * es), one(R))
    qs = R(0.622) * vp1

    return (qs, es)
end

# Evaluate 8th-degree polynomial: c[1] + c[2]*x + c[3]*x^2 + ... + c[9]*x^8.
# Coefficients are converted to the eltype of x (R) so the result stays in the
# working precision — Float64 coefficients would otherwise force a Float64
# result (and break Metal under Float32). Byte-identical when x is Float64.
@inline function _poly8(x::Real, c::NTuple{9,<:Real})
    R = typeof(x)
    return R(c[1]) + x*(R(c[2]) + x*(R(c[3]) + x*(R(c[4]) + x*(R(c[5]) +
           x*(R(c[6]) + x*(R(c[7]) + x*(R(c[8]) + x*R(c[9]))))))))
end
