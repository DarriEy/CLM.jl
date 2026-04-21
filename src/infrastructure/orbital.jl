# ==========================================================================
# Ported from: src/share/shr_orb_mod.F90
# Solar orbital parameters — simplified analytical form
#
# Public functions:
#   compute_orbital  — Solar declination and obliquity from calendar day
#   shr_orb_decl     — Full orbital declination with eccentricity factor
# ==========================================================================

# Default orbital parameters (year 1990 CE, Berger 1978)
const ORB_ECCEN_DEFAULT  = 0.016724  # eccentricity
const ORB_OBLIQ_DEFAULT  = 23.4441   # obliquity [degrees]
const ORB_MVELP_DEFAULT  = 102.7     # moving vernal equinox longitude of perihelion [degrees]
const ORB_OBLIQR_DEFAULT = deg2rad(ORB_OBLIQ_DEFAULT)  # obliquity [radians]

"""
    compute_orbital(calday; eccen, obliqr, mvelpp, lambm0)
    -> (declin, eccf)

Compute solar declination angle and earth-sun distance factor for a
given calendar day using simplified Berger (1978) orbital parameters.

# Arguments
- `calday::Float64` — Calendar day (1.0 = Jan 1, fractional for sub-daily)

# Returns
- `declin::Float64` — Solar declination angle [radians]
- `eccf::Float64`   — Earth-sun distance factor (1/r²)

Ported from `shr_orb_decl` in `shr_orb_mod.F90`.
"""
function compute_orbital(calday::Real;
                         eccen::Real = ORB_ECCEN_DEFAULT,
                         obliqr::Real = ORB_OBLIQR_DEFAULT,
                         mvelpp::Real = deg2rad(ORB_MVELP_DEFAULT + 180.0),
                         lambm0::Real = NaN)
    dayspy = 365.0
    ve = 80.5  # vernal equinox calendar day

    # Compute lambm0 if not provided (mean longitude at vernal equinox)
    if isnan(lambm0)
        beta = sqrt(1.0 - eccen^2)
        lambm0 = 2.0 * (
            (0.5 * eccen + 0.125 * eccen^3) * (1.0 + beta) * sin(mvelpp)
            - 0.25 * eccen^2 * (0.5 + beta) * sin(2.0 * mvelpp)
            + 0.125 * eccen^3 * (1.0 / 3.0 + beta) * sin(3.0 * mvelpp)
        )
    end

    # Mean longitude at present day
    lambm = lambm0 + (calday - ve) * 2.0 * π / dayspy
    lmm = lambm - mvelpp

    # True longitude (Berger 1978, equation of center)
    sinl = sin(lmm)
    lamb = lambm + eccen * (
        2.0 * sinl +
        eccen * (1.25 * sin(2.0 * lmm) +
        eccen * (13.0 / 12.0 * sin(3.0 * lmm) - 0.25 * sinl))
    )

    # Inverse normalized sun-earth distance
    invrho = (1.0 + eccen * cos(lamb - mvelpp)) / (1.0 - eccen^2)

    # Solar declination
    declin = asin(sin(obliqr) * sin(lamb))

    # Earth-sun distance factor
    eccf = invrho^2

    return (declin, eccf)
end
