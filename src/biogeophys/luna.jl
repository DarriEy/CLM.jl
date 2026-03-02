# ==========================================================================
# Ported from: src/biogeophys/LunaMod.F90
# LUNA (Leaf Utilization of Nitrogen for Assimilation) model.
# Calculates photosynthetic capacities based on prescribed leaf nitrogen
# content, using the LUNA model developed by Chonggang Xu, Ashehad Ali
# and Rosie Fisher. Currently only works for C3 plants.
# See Xu et al 2012; Ali et al 2015a. Ecological Applications.
#
# Public functions:
#   update_photosynthesis_capacity!  — Update canopy nitrogen profile
#   acc24_climate_luna!              — Accumulate 24 hr climates
#   acc240_climate_luna!             — Accumulate 10 day climates
#   clear24_climate_luna!            — Clear 24 hr climates
#   is_time_to_run_luna              — Check if it is time to run LUNA
#
# Private functions:
#   nitrogen_allocation!             — Nitrogen partitioning via LUNA
#   nitrogen_investments!            — Nitrogen investment calculations
#   nue_ref                          — Reference NUE
#   nue_calc                         — Current NUE
#   vcmx_t_kattge                    — Vcmax temperature response (Kattge & Knorr 2007)
#   jmx_t_kattge                     — Jmax temperature response (Kattge & Knorr 2007)
#   vcmx_t_leuning                   — Vcmax temperature response (Leuning 2002)
#   jmx_t_leuning                    — Jmax temperature response (Leuning 2002)
#   resp_t_bernacchi                 — Respiration temperature response (Bernacchi 2001)
#   photosynthesis_luna!             — Photosynthesis for nitrogen allocation
#   quadratic_luna                   — Quadratic solver (LUNA-specific)
# ==========================================================================

# --- Module-level constants ---

const LUNA_Cv = 1.2e-5 * 3600.0             # conversion factor from umol CO2 to g carbon
const LUNA_Fc25 = 294.2                      # Fc25 = 6.22*47.3 see Rogers (2014) Photosynthesis Research
const LUNA_Fj25 = 1257.0                     # Fj25 = 8.06*156 see COSTE 2005 and Xu et al 2012
const LUNA_NUEr25 = 33.69                    # nitrogen use efficiency for respiration, see Xu et al 2012
const LUNA_Cb = 1.78                         # nitrogen use efficiency for chlorophyll for light capture, see Evans 1989
const LUNA_O2ref = 209460.0                  # ppm of O2 in the air
const LUNA_CO2ref = 380.0                    # reference CO2 concentration for calculation of reference NUE
const LUNA_forc_pbot_ref = 101325.0          # reference air pressure for calculation of reference NUE
const LUNA_Q10Enz = 2.0                      # Q10 value for enzyme decay rate
const LUNA_NMCp25 = 0.715                    # estimated by assuming 80% maintenance respiration is used for photosynthesis enzyme maintenance
const LUNA_Trange1 = 5.0                     # lower temperature limit (oC) for nitrogen optimization
const LUNA_Trange2 = 42.0                    # upper temperature limit (oC) for nitrogen optimization
const LUNA_SNC = 0.004                       # structural nitrogen concentration (g N g-1 dry mass carbon)
const LUNA_mp = 9.0                          # slope of stomatal conductance
const LUNA_PARLowLim = 200.0                 # minimum photosynthetically active radiation for nitrogen optimization

# ==========================================================================
# LUNA parameters type
# ==========================================================================

"""
    LunaParamsData

Parameters for the LUNA model (from `params_type` in Fortran).
"""
Base.@kwdef mutable struct LunaParamsData
    cp25_yr2000::Float64 = NaN               # CO2 compensation point at 25°C at present day O2 (mol/mol)
    kc25_coef::Float64 = NaN                 # Michaelis-Menten const. at 25°C for CO2 (unitless)
    ko25_coef::Float64 = NaN                 # Michaelis-Menten const. at 25°C for O2 (unitless)
    luna_theta_cj::Float64 = NaN             # LUNA empirical curvature parameter for ac, aj photosynthesis co-limitation
    enzyme_turnover_daily::Float64 = NaN     # Daily turnover rate for photosynthetic enzyme at 25oC
    relhExp::Float64 = NaN                   # Impact of relative humidity on electron transport rate
    minrelh::Float64 = NaN                   # Minimum relative humidity for nitrogen optimization (fraction)
    jmaxb0::Vector{Float64} = Float64[]      # Baseline proportion of nitrogen allocated for electron transport
    jmaxb1::Vector{Float64} = Float64[]      # Coefficient for electron transport rate response to light
    wc2wjb0::Vector{Float64} = Float64[]     # Baseline ratio of rubisco limited vs light limited rate (Wc:Wj)
end

"""
    luna_params_init!(lp, mxpft)

Allocate and initialize LUNA parameters with NaN.
"""
function luna_params_init!(lp::LunaParamsData, mxpft::Int)
    lp.jmaxb0  = fill(NaN, mxpft + 1)
    lp.jmaxb1  = fill(NaN, mxpft + 1)
    lp.wc2wjb0 = fill(NaN, mxpft + 1)
    return nothing
end

"""
    luna_params_clean!(lp)

Deallocate LUNA parameters.
"""
function luna_params_clean!(lp::LunaParamsData)
    lp.jmaxb0  = Float64[]
    lp.jmaxb1  = Float64[]
    lp.wc2wjb0 = Float64[]
    return nothing
end

# ==========================================================================
# Temperature response functions
# ==========================================================================

"""
    vcmx_t_kattge(tgrow, tleaf) -> Float64

Temperature response for Vcmax, with temperature acclimation (Kattge & Knorr 2007).
"""
function vcmx_t_kattge(tgrow::Float64, tleaf::Float64)
    TlimVcmx = 668.39 - 1.07 * min(max(tgrow, 11.0), 35.0)
    Vcmxf1 = 1.0 + exp((TlimVcmx * (25.0 + TFRZ) - 200000.0) / (RGAS * (25.0 + TFRZ)))
    Vcmxf2 = exp((72000.0 / (RGAS * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    Vcmxf3 = 1.0 + exp((TlimVcmx * (tleaf + TFRZ) - 200000.0) / (RGAS * (tleaf + TFRZ)))
    return Vcmxf1 * Vcmxf2 / Vcmxf3
end

"""
    jmx_t_kattge(tgrow, tleaf) -> Float64

Temperature response for Jmax, with temperature acclimation (Kattge & Knorr 2007).
"""
function jmx_t_kattge(tgrow::Float64, tleaf::Float64)
    TlimJmx = 659.7 - 0.75 * min(max(tgrow, 11.0), 35.0)
    Jmxf1 = 1.0 + exp((TlimJmx * (25.0 + TFRZ) - 200000.0) / (RGAS * (25.0 + TFRZ)))
    Jmxf2 = exp((50000.0 / (RGAS * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (tleaf + TFRZ)))
    Jmxf3 = 1.0 + exp((TlimJmx * (tleaf + TFRZ) - 200000.0) / (RGAS * (tleaf + TFRZ)))
    return Jmxf1 * Jmxf2 / Jmxf3
end

"""
    vcmx_t_leuning(tgrow, tleaf) -> Float64

Temperature response for Vcmax without temperature acclimation (Leuning 2002).
"""
function vcmx_t_leuning(tgrow::Float64, tleaf::Float64)
    TlimVcmx = 486.0
    Vcmxf1 = 1.0 + exp((TlimVcmx * (25.0 + TFRZ) - 149252.0) / (RGAS * (25.0 + TFRZ)))
    Vcmxf2 = exp((73637.0 / (RGAS * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    Vcmxf3 = 1.0 + exp((TlimVcmx * (tleaf + TFRZ) - 149252.0) / (RGAS * (tleaf + TFRZ)))
    return Vcmxf1 * Vcmxf2 / Vcmxf3
end

"""
    jmx_t_leuning(tgrow, tleaf) -> Float64

Temperature response for Jmax without temperature acclimation (Leuning 2002).
"""
function jmx_t_leuning(tgrow::Float64, tleaf::Float64)
    TlimJmx = 495.0
    Jmxf1 = 1.0 + exp((TlimJmx * (25.0 + TFRZ) - 152044.0) / (RGAS * (25.0 + TFRZ)))
    Jmxf2 = exp((50300.0 / (RGAS * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    Jmxf3 = 1.0 + exp((TlimJmx * (tleaf + TFRZ) - 152044.0) / (RGAS * (tleaf + TFRZ)))
    return Jmxf1 * Jmxf2 / Jmxf3
end

"""
    resp_t_bernacchi(tleaf) -> Float64

Temperature response for respiration (Bernacchi PCE 2001).
"""
function resp_t_bernacchi(tleaf::Float64)
    return exp(18.72 - 46.39 / (RGAS * 1.0e-3 * (tleaf + TFRZ)))
end

# ==========================================================================
# Quadratic solver (LUNA-specific, matches Fortran Quadratic subroutine)
# ==========================================================================

"""
    quadratic_luna(a, b, c) -> (r1, r2)

Solve a*x^2 + b*x + c = 0. Matches the Fortran LUNA Quadratic subroutine exactly.
"""
function quadratic_luna(a::Float64, b::Float64, c::Float64)
    r1 = 1.0e36
    r2 = 1.0e36

    if a == 0.0
        return (r1, r2)
    end

    if b >= 0.0
        q = -0.5 * (b + sqrt(b * b - 4.0 * a * c))
    else
        q = -0.5 * (b - sqrt(b * b - 4.0 * a * c))
    end

    r1 = q / a

    if q != 0.0
        r2 = c / q
    else
        r2 = 1.0e36
    end

    return (r1, r2)
end

# ==========================================================================
# Reference NUE calculation
# ==========================================================================

"""
    nue_ref(luna_params) -> (NUEjref, NUEcref, Kj2Kcref)

Calculate reference nitrogen use efficiency at 25°C and 380 ppm CO2.
"""
function nue_ref(luna_params::LunaParamsData)
    tgrow = 25.0
    tleaf = 25.0
    Fc = vcmx_t_kattge(tgrow, tleaf) * LUNA_Fc25
    Fj = jmx_t_kattge(tgrow, tleaf) * LUNA_Fj25
    CO2c = LUNA_CO2ref * LUNA_forc_pbot_ref * 1.0e-6  # Pa
    O2c  = LUNA_O2ref  * LUNA_forc_pbot_ref * 1.0e-6   # Pa
    k_c = luna_params.kc25_coef * exp((79430.0 / (RGAS * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    k_o = luna_params.ko25_coef * exp((36380.0 / (RGAS * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    c_p = luna_params.cp25_yr2000 * exp((37830.0 / (RGAS * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    awc = k_c * (1.0 + O2c / k_o)
    ci = 0.7 * CO2c
    Kj = max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p)
    Kc = max(ci - c_p, 0.0) / (ci + awc)
    NUEjref = Kj * Fj
    NUEcref = Kc * Fc
    Kj2Kcref = Kj / Kc
    return (NUEjref, NUEcref, Kj2Kcref)
end

# ==========================================================================
# Current NUE calculation
# ==========================================================================

"""
    nue_calc(O2a, ci, tgrow, tleaf, luna_params) -> (NUEj, NUEc, Kj2Kc)

Calculate nitrogen use efficiency under current environmental conditions.
"""
function nue_calc(O2a::Float64, ci::Float64, tgrow::Float64, tleaf::Float64,
                  luna_params::LunaParamsData)
    Fc = vcmx_t_kattge(tgrow, tleaf) * LUNA_Fc25
    Fj = jmx_t_kattge(tgrow, tleaf) * LUNA_Fj25
    k_c = luna_params.kc25_coef * exp((79430.0 / (RGAS * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    k_o = luna_params.ko25_coef * exp((36380.0 / (RGAS * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    c_p = luna_params.cp25_yr2000 * exp((37830.0 / (RGAS * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    awc = k_c * (1.0 + O2a / k_o)
    Kj = max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p)
    Kc = max(ci - c_p, 0.0) / (ci + awc)
    NUEj = Kj * Fj
    NUEc = Kc * Fc
    Kj2Kc = Kj / Kc
    return (NUEj, NUEc, Kj2Kc)
end

# ==========================================================================
# Photosynthesis for LUNA nitrogen allocation
# ==========================================================================

"""
    photosynthesis_luna!(forc_pbot, tleafd, relh, CO2a, O2a, rb, Vcmax, JmeanL, luna_params)
        -> (ci, Kc, Kj, A)

Solve photosynthesis equations for LUNA nitrogen allocation.
Uses Ball-Berry stomatal conductance and Farquhar model.
"""
function photosynthesis_luna!(forc_pbot::Float64, tleafd::Float64, relh::Float64,
                              CO2a::Float64, O2a::Float64, rb::Float64,
                              Vcmax::Float64, JmeanL::Float64,
                              luna_params::LunaParamsData)
    rsmax0 = 2.0 * 1.0e4
    bp = 2000.0
    tleaf = tleafd
    tleafk = tleaf + TFRZ
    aquad = 1.0
    relhc = max(luna_params.minrelh, relh)
    bbb = 1.0 / bp
    mbb = LUNA_mp
    CO2c = CO2a
    O2c = O2a
    ci = 0.7 * CO2c
    ciold = ci - 0.02
    cf = forc_pbot / (8.314 * tleafk) * 1.0e6
    gb_mol = cf / rb
    k_c = luna_params.kc25_coef * exp((79430.0 / (8.314 * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    k_o = luna_params.ko25_coef * exp((36380.0 / (8.314 * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    c_p = luna_params.cp25_yr2000 * exp((37830.0 / (8.314 * (25.0 + TFRZ))) * (1.0 - (TFRZ + 25.0) / (TFRZ + tleaf)))
    awc = k_c * (1.0 + O2c / k_o)

    # Initialize variables used across loop scopes
    gs_mol = bbb
    Kc_val = 0.0
    Wc = 0.0
    Wj = 0.0

    # Rubisco limitation iteration
    i = 1
    while abs(ci - ciold) > 0.01 && i < 100
        i += 1
        ciold = ci
        Kc_val = max(ci - c_p, 0.0) / (ci + awc)
        Wc = Kc_val * Vcmax
        gs_mol = bbb + mbb * Wc / CO2c * forc_pbot * relhc
        phi = forc_pbot * (1.37 * gs_mol + 1.6 * gb_mol) / (gb_mol * gs_mol)
        bquad = awc - CO2c + phi * Vcmax
        cquad = -(c_p * phi * Vcmax + awc * CO2c)
        r1, r2 = quadratic_luna(aquad, bquad, cquad)
        ci = max(r1, r2)
        if ci < 0.0
            ci = c_p + 0.5 * ciold
        end
    end

    Kj = max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p)
    Kc_val = max(ci - c_p, 0.0) / (ci + awc)
    Wc = Kc_val * Vcmax
    Wj = Kj * JmeanL
    ciold = ci - 0.02

    # Light limitation iteration (if Wj < Wc)
    if Wj < Wc
        i = 1
        while abs(ci - ciold) > 0.01 && i < 100
            i += 1
            ciold = ci
            gs_mol = bbb + mbb * Wj / CO2c * forc_pbot * relhc
            phi = forc_pbot * (1.37 * gs_mol + 1.6 * gb_mol) / (gb_mol * gs_mol)
            bquad = 2.0 * c_p - CO2c + phi * JmeanL / 4.0
            cquad = -(c_p * phi * JmeanL / 4.0 + 2.0 * c_p * CO2c)
            r1, r2 = quadratic_luna(aquad, bquad, cquad)
            ci = max(r1, r2)
            if ci < 0.0
                ci = c_p + 0.5 * ciold
            end
            Kj = max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p)
            Wj = Kj * JmeanL
        end
        Kj = max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p)
        Kc_val = max(ci - c_p, 0.0) / (ci + awc)
        Wc = Kc_val * Vcmax
        Wj = Kj * JmeanL
    end

    A = (1.0 - luna_params.luna_theta_cj) * max(Wc, Wj) + luna_params.luna_theta_cj * min(Wc, Wj)
    rs = cf / gs_mol
    rs = min(rsmax0, rs)

    return (ci, Kc_val, Kj, A)
end

# ==========================================================================
# Nitrogen investments calculation
# ==========================================================================

"""
    nitrogen_investments!(KcKjFlag, FNCa, Nlc, forc_pbot10, relh10, CO2a10, O2a10,
                          PARi10, PARimx10, rb10, hourpd, tair10, tleafd10, tleafn10,
                          Kj2Kc, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref,
                          NUEr, o3coefjmax, jmaxb0, wc2wjb0, Kc_in, Kj_in, ci_in,
                          luna_params) -> (Vcmax, Jmax, JmeanL, JmaxL, Net, Ncb, Nresp, PSN, RESP, Kc, Kj, ci)

Calculate nitrogen investment for electron transport, carboxylation, respiration
given nitrogen allocation in light capture [Nlc].
"""
function nitrogen_investments!(KcKjFlag::Int, FNCa::Float64, Nlc::Float64,
                               forc_pbot10::Float64, relh10::Float64,
                               CO2a10::Float64, O2a10::Float64,
                               PARi10::Float64, PARimx10::Float64,
                               rb10::Float64, hourpd::Float64,
                               tair10::Float64, tleafd10::Float64, tleafn10::Float64,
                               Kj2Kc::Float64, JmaxCoef::Float64,
                               Fc::Float64, Fj::Float64,
                               NUEc::Float64, NUEj::Float64,
                               NUEcref::Float64, NUEjref::Float64,
                               NUEr::Float64, o3coefjmax::Float64,
                               jmaxb0::Float64, wc2wjb0::Float64,
                               Kc_in::Float64, Kj_in::Float64, ci_in::Float64,
                               luna_params::LunaParamsData)
    leaf_mr_vcm = 0.015

    theta = 0.292 / (1.0 + 0.076 / (Nlc * LUNA_Cb))
    ELTRNabsorb = theta * PARi10
    Jmaxb0act = jmaxb0 * FNCa * Fj

    # o3coefjmax defaults to 1 unless ozone stress_method == 'stress_falk'
    Jmax = Jmaxb0act + JmaxCoef * ELTRNabsorb * o3coefjmax

    JmaxL = theta * PARimx10 / sqrt(1.0 + (theta * PARimx10 / Jmax)^2.0)
    NUEchg = (NUEc / NUEcref) * (NUEjref / NUEj)
    Wc2Wj = wc2wjb0 * (NUEchg^0.5)
    Vcmax = Wc2Wj * JmaxL * Kj2Kc
    JmeanL = theta * PARi10 / sqrt(1.0 + (ELTRNabsorb / Jmax)^2.0)

    Kc = Kc_in
    Kj = Kj_in
    ci = ci_in

    if KcKjFlag == 0  # update Kc, Kj and ci from photosynthesis
        ci, Kc, Kj, A = photosynthesis_luna!(forc_pbot10, tleafd10, relh10, CO2a10, O2a10, rb10, Vcmax, JmeanL, luna_params)
    else
        Wc = Kc * Vcmax
        Wj = Kj * JmeanL
        A = (1.0 - luna_params.luna_theta_cj) * max(Wc, Wj) + luna_params.luna_theta_cj * min(Wc, Wj)
    end

    PSN = LUNA_Cv * A * hourpd
    Vcmaxnight = vcmx_t_kattge(tair10, tleafn10) / vcmx_t_kattge(tair10, tleafd10) * Vcmax
    RESP = LUNA_Cv * leaf_mr_vcm * (Vcmax * hourpd + Vcmaxnight * (24.0 - hourpd))
    Net = Jmax / Fj
    Ncb = Vcmax / Fc
    Nresp = RESP / NUEr

    return (Vcmax, Jmax, JmeanL, JmaxL, Net, Ncb, Nresp, PSN, RESP, Kc, Kj, ci)
end

# ==========================================================================
# Nitrogen allocation optimization
# ==========================================================================

"""
    nitrogen_allocation!(FNCa, forc_pbot10, relh10, CO2a10, O2a10, PARi10, PARimx10,
                         rb10, hourpd, tair10, tleafd10, tleafn10,
                         jmaxb0, jmaxb1, wc2wjb0, PNlcold, PNetold,
                         PNrespold, PNcbold, dayl_factor, o3coefjmax, luna_params)
                         -> (PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt)

LUNA model nitrogen partitioning optimization.
"""
function nitrogen_allocation!(FNCa::Float64, forc_pbot10::Float64, relh10::Float64,
                              CO2a10::Float64, O2a10::Float64,
                              PARi10::Float64, PARimx10::Float64,
                              rb10::Float64, hourpd::Float64,
                              tair10::Float64, tleafd10::Float64, tleafn10::Float64,
                              jmaxb0_val::Float64, jmaxb1_val::Float64, wc2wjb0_val::Float64,
                              PNlcold::Float64, PNetold::Float64,
                              PNrespold::Float64, PNcbold::Float64,
                              dayl_factor::Float64, o3coefjmax::Float64,
                              luna_params::LunaParamsData)

    NUEjref, NUEcref, Kj2Kcref = nue_ref(luna_params)
    Nlc   = PNlcold * FNCa
    Net   = PNetold * FNCa
    Nresp = PNrespold * FNCa
    Ncb   = PNcbold * FNCa
    if Nlc > FNCa * 0.5
        Nlc = 0.5 * FNCa
    end
    chg_per_step = 0.02 * FNCa
    PNlc = PNlcold
    PNlcoldi = PNlcold - 0.001
    PARi10c  = max(LUNA_PARLowLim, PARi10)
    PARimx10c = max(LUNA_PARLowLim, PARimx10)
    increase_flag = 0
    jj = 1
    tleafd10c = min(max(tleafd10, LUNA_Trange1), LUNA_Trange2)
    tleafn10c = min(max(tleafn10, LUNA_Trange1), LUNA_Trange2)
    ci = 0.7 * CO2a10
    JmaxCoef = jmaxb1_val * dayl_factor * (1.0 - exp(-luna_params.relhExp * max(relh10 -
        luna_params.minrelh, 0.0) / (1.0 - luna_params.minrelh)))

    # Initialize Kc, Kj and loop variables for proper scoping
    Kc = 0.0
    Kj = 0.0
    Nstore = 0.0
    Npsntarget = 0.0
    PSN = 0.0

    while PNlcoldi != PNlc && jj < 100
        Fc = vcmx_t_kattge(tair10, tleafd10c) * LUNA_Fc25
        Fj = jmx_t_kattge(tair10, tleafd10c) * LUNA_Fj25
        NUEr = LUNA_Cv * LUNA_NUEr25 * (resp_t_bernacchi(tleafd10c) * hourpd +
            resp_t_bernacchi(tleafn10c) * (24.0 - hourpd))

        # Calculate NUE
        NUEj, NUEc, Kj2Kc = nue_calc(O2a10, ci, tair10, tleafd10c, luna_params)

        # Baseline nitrogen investments
        KcKjFlag = 0
        Vcmax, Jmax, JmeanL, JmaxL, Net, Ncb, Nresp, PSN, RESP, Kc, Kj, ci =
            nitrogen_investments!(KcKjFlag, FNCa, Nlc, forc_pbot10, relh10, CO2a10, O2a10,
                PARi10c, PARimx10c, rb10, hourpd, tair10, tleafd10c, tleafn10c,
                Kj2Kc, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, o3coefjmax,
                jmaxb0_val, wc2wjb0_val, Kc, Kj, ci, luna_params)

        Npsntarget = Nlc + Ncb + Net
        PNlcoldi = Nlc / FNCa
        Nstore = FNCa - Npsntarget - Nresp

        # Test increase of light capture nitrogen
        if Nstore > 0.0 && (increase_flag == 1 || jj == 1)
            Nlc2 = Nlc + chg_per_step
            if Nlc2 / FNCa > 0.95
                Nlc2 = 0.95 * FNCa
            end
            KcKjFlag = 1
            _, _, _, _, Net2, Ncb2, Nresp2, PSN2, RESP2, _, _, _ =
                nitrogen_investments!(KcKjFlag, FNCa, Nlc2, forc_pbot10, relh10, CO2a10, O2a10,
                    PARi10c, PARimx10c, rb10, hourpd, tair10, tleafd10c, tleafn10c,
                    Kj2Kc, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, o3coefjmax,
                    jmaxb0_val, wc2wjb0_val, Kc, Kj, ci, luna_params)

            Npsntarget2 = Nlc2 + Ncb2 + Net2
            Carboncost2 = (Npsntarget2 - Npsntarget) * LUNA_NMCp25 * LUNA_Cv *
                (resp_t_bernacchi(tleafd10c) * hourpd + resp_t_bernacchi(tleafn10c) * (24.0 - hourpd))
            Carbongain2 = PSN2 - PSN
            if Carbongain2 > Carboncost2 && (Npsntarget2 + Nresp2 < 0.95 * FNCa)
                Nlc = Nlc2
                Net = Net2
                Ncb = Ncb2
                Nstore = FNCa - Npsntarget2 - Nresp2
                if jj == 1
                    increase_flag = 1
                end
            end
        end

        # Test decrease of light capture nitrogen
        if increase_flag == 0
            if Nstore < 0.0
                Nlc1 = Nlc * 0.8  # bigger step of decrease if negative
            else
                Nlc1 = Nlc - chg_per_step
            end
            if Nlc1 < 0.05
                Nlc1 = 0.05
            end
            KcKjFlag = 1
            _, _, _, _, Net1, Ncb1, Nresp1, PSN1, RESP1, _, _, _ =
                nitrogen_investments!(KcKjFlag, FNCa, Nlc1, forc_pbot10, relh10, CO2a10, O2a10,
                    PARi10c, PARimx10c, rb10, hourpd, tair10, tleafd10c, tleafn10c,
                    Kj2Kc, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, o3coefjmax,
                    jmaxb0_val, wc2wjb0_val, Kc, Kj, ci, luna_params)

            Npsntarget1 = Nlc1 + Ncb1 + Net1
            Carboncost1 = (Npsntarget - Npsntarget1) * LUNA_NMCp25 * LUNA_Cv *
                (resp_t_bernacchi(tleafd10c) * hourpd + resp_t_bernacchi(tleafn10c) * (24.0 - hourpd))
            Carbongain1 = PSN - PSN1
            if (Carbongain1 < Carboncost1 && Nlc1 > 0.05) || (Npsntarget + Nresp) > 0.95 * FNCa
                Nlc = Nlc1
                Net = Net1
                Ncb = Ncb1
                Nstore = FNCa - Npsntarget1 - Nresp1
            end
        end

        PNlc = Nlc / FNCa
        jj += 1
    end

    PNlcopt    = Nlc / FNCa
    PNstoreopt = Nstore / FNCa
    PNcbopt    = Ncb / FNCa
    PNetopt    = Net / FNCa
    PNrespopt  = Nresp / FNCa

    return (PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt)
end

# ==========================================================================
# is_time_to_run_luna
# ==========================================================================

"""
    is_time_to_run_luna(is_end_curr_day) -> Bool

Check if it is time to run the LUNA module (end of current day).
"""
function is_time_to_run_luna(is_end_curr_day::Bool)
    return is_end_curr_day
end

# ==========================================================================
# Clear 24hr climate accumulations for LUNA
# ==========================================================================

"""
    clear24_climate_luna!(solarabs, photosyns, temperature, patchdata, mask_patch, bounds)

Zero out the 24 hr climate accumulations for the LUNA model.
"""
function clear24_climate_luna!(solarabs::SolarAbsorbedData,
                               photosyns::PhotosynthesisData,
                               temperature::TemperatureData,
                               patchdata::PatchData,
                               mask_patch::BitVector,
                               bounds::UnitRange{Int})
    for p in bounds
        mask_patch[p] || continue
        temperature.t_veg_day_patch[p]   = 0.0
        temperature.t_veg_night_patch[p] = 0.0
        solarabs.par24d_z_patch[p, :]   .= 0.0
        solarabs.par24x_z_patch[p, :]   .= 0.0
        photosyns.fpsn24_patch[p]        = 0.0
        temperature.nnightsteps_patch[p] = 0
        temperature.ndaysteps_patch[p]   = 0
    end
    return nothing
end

# ==========================================================================
# Accumulate 24hr climate for LUNA
# ==========================================================================

"""
    acc24_climate_luna!(canopystate, photosyns, surfalb, solarabs, temperature,
                        patchdata, mask_patch, bounds, dtime)

Accumulate the 24 hr climates for the LUNA model.
Called each timestep from CanopyFluxes.
"""
function acc24_climate_luna!(canopystate::CanopyStateData,
                             photosyns::PhotosynthesisData,
                             surfalb::SurfaceAlbedoData,
                             solarabs::SolarAbsorbedData,
                             temperature::TemperatureData,
                             patchdata::PatchData,
                             mask_patch::BitVector,
                             bounds::UnitRange{Int},
                             dtime::Float64)
    for p in bounds
        mask_patch[p] || continue

        # Check whether it is the first day
        if temperature.t_veg_day_patch[p] != SPVAL
            if solarabs.sabv_patch[p] > 0.0
                temperature.t_veg_day_patch[p] += temperature.t_veg_patch[p]
                temperature.ndaysteps_patch[p] += 1
            else
                temperature.t_veg_night_patch[p] += temperature.t_veg_patch[p]
                temperature.nnightsteps_patch[p] += 1
            end

            nrad_p = surfalb.nrad_patch[p]
            for z in 1:nrad_p
                # Average of sunlit and shaded leaves
                tlaii = canopystate.laisun_z_patch[p, z] + canopystate.laisha_z_patch[p, z]
                if tlaii > 0.0
                    # RF & GBB: Make LUNA predict sunlit fraction N fractionation
                    TRad = solarabs.parsun_z_patch[p, z]
                    solarabs.par24d_z_patch[p, z] += dtime * TRad
                    if TRad > solarabs.par24x_z_patch[p, z]
                        solarabs.par24x_z_patch[p, z] = TRad
                    end
                end
            end

            photosyns.fpsn24_patch[p] += dtime * photosyns.fpsn_patch[p]
        end
    end
    return nothing
end

# ==========================================================================
# Accumulate 10-day (240hr) running mean climate for LUNA
# ==========================================================================

"""
    acc240_climate_luna!(temperature, photosyns, surfalb, solarabs, waterdiag,
                         frictionvel, patchdata, mask_patch, bounds,
                         oair, cair, rb, rh, dtime)

Accumulate 10-day running mean climates for the LUNA model.
Called at end of day from CanopyFluxes.
"""
function acc240_climate_luna!(temperature::TemperatureData,
                              photosyns::PhotosynthesisData,
                              surfalb::SurfaceAlbedoData,
                              solarabs::SolarAbsorbedData,
                              waterdiag::WaterDiagnosticBulkData,
                              frictionvel::FrictionVelocityData,
                              patchdata::PatchData,
                              mask_patch::BitVector,
                              bounds::UnitRange{Int},
                              oair::Vector{Float64},
                              cair::Vector{Float64},
                              rb::Vector{Float64},
                              rh::Vector{Float64},
                              dtime::Float64)
    for p in bounds
        mask_patch[p] || continue

        if temperature.t_veg_day_patch[p] != SPVAL
            # Calculate 10-day running mean radiations
            ndaysteps_p = temperature.ndaysteps_patch[p]
            if ndaysteps_p > 0
                par24d_z_i = solarabs.par24d_z_patch[p, :] ./ (dtime * ndaysteps_p)
            else
                par24d_z_i = zeros(size(solarabs.par24d_z_patch, 2))
            end

            if solarabs.par240d_z_patch[p, 1] == SPVAL  # first day
                solarabs.par240x_z_patch[p, :] .= solarabs.par24x_z_patch[p, :]
                solarabs.par240d_z_patch[p, :] .= par24d_z_i
            else
                solarabs.par240x_z_patch[p, :] .= 0.9 .* solarabs.par240x_z_patch[p, :] .+ 0.1 .* solarabs.par24x_z_patch[p, :]
                solarabs.par240d_z_patch[p, :] .= 0.9 .* solarabs.par240d_z_patch[p, :] .+ 0.1 .* par24d_z_i
            end

            # 10-day running mean daytime temperature
            if ndaysteps_p > 0
                t_veg_dayi = temperature.t_veg_day_patch[p] / ndaysteps_p
            else
                nnightsteps_p = temperature.nnightsteps_patch[p]
                t_veg_dayi = temperature.t_veg_night_patch[p] / nnightsteps_p
            end
            if temperature.t_veg10_day_patch[p] == SPVAL
                temperature.t_veg10_day_patch[p] = t_veg_dayi
            end
            temperature.t_veg10_day_patch[p] = 0.9 * temperature.t_veg10_day_patch[p] + 0.1 * t_veg_dayi

            # 10-day running mean nighttime temperature
            nnightsteps_p = temperature.nnightsteps_patch[p]
            if nnightsteps_p > 0
                t_veg_nighti = temperature.t_veg_night_patch[p] / nnightsteps_p
            else
                t_veg_nighti = temperature.t_veg_day_patch[p] / ndaysteps_p
            end
            if temperature.t_veg10_night_patch[p] == SPVAL
                temperature.t_veg10_night_patch[p] = t_veg_nighti
            end
            temperature.t_veg10_night_patch[p] = 0.9 * temperature.t_veg10_night_patch[p] + 0.1 * t_veg_nighti

            # 10-day running mean relative humidity
            if waterdiag.rh10_af_patch[p] == SPVAL
                waterdiag.rh10_af_patch[p] = rh[p]
            end
            waterdiag.rh10_af_patch[p] = 0.9 * waterdiag.rh10_af_patch[p] + 0.1 * min(1.0, rh[p])

            # 10-day running mean boundary layer resistance
            if frictionvel.rb10_patch[p] == SPVAL
                frictionvel.rb10_patch[p] = rb[p]
            end
            frictionvel.rb10_patch[p] = 0.9 * frictionvel.rb10_patch[p] + 0.1 * rb[p]
        end
    end
    return nothing
end

# ==========================================================================
# Update photosynthetic capacity (main LUNA driver)
# ==========================================================================

"""
    update_photosynthesis_capacity!(photosyns, temperature, canopystate, surfalb, solarabs,
                                    waterdiag, frictionvel, patchdata, gridcell,
                                    mask_patch, bounds, dayl_factor,
                                    forc_pbot10, CO2_p240, O2_p240,
                                    c3psn_pft, slatop_pft, leafcn_pft, rhol_pft, taul_pft,
                                    o3coefjmax, luna_params, dtime, nlevcan)

Update Vcmax25 and Jmax25 profiles using the LUNA nitrogen allocation model.
"""
function update_photosynthesis_capacity!(photosyns::PhotosynthesisData,
                                         temperature::TemperatureData,
                                         canopystate::CanopyStateData,
                                         surfalb::SurfaceAlbedoData,
                                         solarabs::SolarAbsorbedData,
                                         waterdiag::WaterDiagnosticBulkData,
                                         frictionvel::FrictionVelocityData,
                                         patchdata::PatchData,
                                         gridcell::GridcellData,
                                         mask_patch::BitVector,
                                         bounds::UnitRange{Int},
                                         dayl_factor::Vector{Float64},
                                         forc_pbot10::Vector{Float64},
                                         CO2_p240::Vector{Float64},
                                         O2_p240::Vector{Float64},
                                         c3psn_pft::Vector{Float64},
                                         slatop_pft::Vector{Float64},
                                         leafcn_pft::Vector{Float64},
                                         rhol_pft::Matrix{Float64},
                                         taul_pft::Matrix{Float64},
                                         o3coefjmax::Vector{Float64},
                                         luna_params::LunaParamsData,
                                         dtime::Float64,
                                         nlevcan_val::Int)
    fnps = 0.15
    FNCa_z = zeros(nlevcan_val)

    for p in bounds
        mask_patch[p] || continue

        ft = patchdata.itype[p]
        g  = patchdata.gridcell[p]
        c  = patchdata.column[p]

        # Check whether it is the first day
        if temperature.t_veg_day_patch[p] != SPVAL
            # Get climate drivers
            CO2a10   = CO2_p240[p]
            O2a10    = O2_p240[p]
            hourpd   = gridcell.dayl[g] / 3600.0
            tleafd10 = temperature.t_veg10_day_patch[p] - TFRZ
            tleafn10 = temperature.t_veg10_night_patch[p] - TFRZ
            tleaf10  = (gridcell.dayl[g] * tleafd10 + (86400.0 - gridcell.dayl[g]) * tleafn10) / 86400.0
            tair10   = temperature.t_a10_patch[p] - TFRZ
            relh10   = min(1.0, waterdiag.rh10_af_patch[p])
            rb10v    = frictionvel.rb10_patch[p]

            # Enzyme turnover rate
            EnzTurnoverTFactor = LUNA_Q10Enz^(0.1 * (min(40.0, tleaf10) - 25.0))
            max_daily_pchg = EnzTurnoverTFactor * luna_params.enzyme_turnover_daily

            # Radiation absorption
            rabsorb = 1.0 - rhol_pft[ft, 1] - taul_pft[ft, 1]

            # Nitrogen allocation model
            if canopystate.tlai_patch[p] > 0.0 && photosyns.lnca_patch[p] > 0.0
                RadTop = solarabs.par240d_z_patch[p, 1] / rabsorb
                PARTop = RadTop * 4.6  # conversion from W/m2 to umol/m2/s

                if round(Int, c3psn_pft[ft]) == 1
                    if photosyns.fpsn24_patch[p] > 0.0  # only optimize if growth and C3
                        nrad_p = surfalb.nrad_patch[p]

                        # First pass: compute FNCa_z profile
                        for z in 1:nrad_p
                            if surfalb.tlai_z_patch[p, z] > 0.0
                                qabs = solarabs.par240d_z_patch[p, z] / rabsorb
                                PARi10 = qabs * 4.6
                            else
                                PARi10 = 0.01
                            end
                            relRad = PARi10 / PARTop
                            relCLNCa = 0.1802 * log(relRad) + 1.0  # see Ali et al 2015
                            relCLNCa = max(0.2, relCLNCa)
                            relCLNCa = min(1.0, relCLNCa)

                            SNCa = 1.0 / slatop_pft[ft] * LUNA_SNC
                            lnc_p = photosyns.lnca_patch[p]
                            if 0.9 * lnc_p > SNCa
                                FNCa_z[z] = relCLNCa * (lnc_p - SNCa)
                            else
                                FNCa_z[z] = relCLNCa * 0.1 * lnc_p
                            end
                        end

                        # Second pass: nitrogen allocation per layer
                        for z in 1:nrad_p
                            FNCa = FNCa_z[z]
                            if FNCa > 15.0  # boundary check
                                FNCa = 15.0
                                @warn "LUNA: leaf nitrogen content unrealistically high (>15.0 g N/m2 leaf) for patch=$p z=$z pft=$ft"
                            end

                            radmax2mean = solarabs.par240x_z_patch[p, z] / solarabs.par240d_z_patch[p, z]
                            if surfalb.tlai_z_patch[p, z] > 0.0
                                qabs = solarabs.par240d_z_patch[p, z] / rabsorb
                                PARi10 = qabs * 4.6
                            else
                                PARi10 = 0.01
                            end
                            PARimx10 = PARi10 * radmax2mean

                            # Nitrogen allocation model
                            PNlcold = photosyns.pnlc_z_patch[p, z]
                            PNetold = 0.0
                            PNrespold = 0.0
                            PNcbold = 0.0

                            PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt =
                                nitrogen_allocation!(FNCa, forc_pbot10[p], relh10, CO2a10, O2a10,
                                    PARi10, PARimx10, rb10v, hourpd, tair10, tleafd10, tleafn10,
                                    luna_params.jmaxb0[ft], luna_params.jmaxb1[ft], luna_params.wc2wjb0[ft],
                                    PNlcold, PNetold, PNrespold, PNcbold, dayl_factor[p],
                                    o3coefjmax[p], luna_params)

                            vcmx25_opt = PNcbopt * FNCa * LUNA_Fc25
                            jmx25_opt  = PNetopt * FNCa * LUNA_Fj25

                            # Constrained update of vcmx25
                            chg = vcmx25_opt - photosyns.vcmx25_z_patch[p, z]
                            chg_constrn = min(abs(chg), photosyns.vcmx25_z_patch[p, z] * max_daily_pchg)
                            photosyns.vcmx25_z_patch[p, z] += sign(chg) * chg_constrn
                            photosyns.vcmx25_z_last_valid_patch[p, z] = photosyns.vcmx25_z_patch[p, z]

                            # Constrained update of jmx25
                            chg = jmx25_opt - photosyns.jmx25_z_patch[p, z]
                            chg_constrn = min(abs(chg), photosyns.jmx25_z_patch[p, z] * max_daily_pchg)
                            photosyns.jmx25_z_patch[p, z] += sign(chg) * chg_constrn
                            photosyns.jmx25_z_last_valid_patch[p, z] = photosyns.jmx25_z_patch[p, z]

                            photosyns.pnlc_z_patch[p, z] = PNlcopt

                            if photosyns.enzs_z_patch[p, z] < 1.0
                                photosyns.enzs_z_patch[p, z] *= (1.0 + max_daily_pchg)
                            end

                            # Sanity checks
                            if isnan(photosyns.vcmx25_z_patch[p, z])
                                error("LUNA: Vcmx25 is NaN for patch=$p z=$z pft=$ft")
                            end
                            if photosyns.vcmx25_z_patch[p, z] > 1000.0 || photosyns.vcmx25_z_patch[p, z] < 0.0
                                @warn "LUNA: Vcmx25 unrealistic (>1000 or negative) for patch=$p z=$z pft=$ft, resetting to 50"
                                photosyns.vcmx25_z_patch[p, z] = 50.0
                            end
                            if isnan(photosyns.jmx25_z_patch[p, z])
                                error("LUNA: Jmx25 is NaN for patch=$p z=$z pft=$ft")
                            end
                            if photosyns.jmx25_z_patch[p, z] > 2000.0 || photosyns.jmx25_z_patch[p, z] < 0.0
                                @warn "LUNA: Jmx25 unrealistic (>2000 or negative) for patch=$p z=$z pft=$ft, resetting to 85"
                                photosyns.jmx25_z_patch[p, z] = 85.0
                            end
                        end  # z loop (growth)
                    else  # decay during drought or winter
                        max_daily_decay = min(0.5, 0.1 * max_daily_pchg)
                        nrad_p = surfalb.nrad_patch[p]
                        for z in 1:nrad_p
                            if photosyns.enzs_z_patch[p, z] > 0.5
                                photosyns.enzs_z_patch[p, z] *= (1.0 - max_daily_decay)
                                photosyns.jmx25_z_patch[p, z] *= (1.0 - max_daily_decay)
                                photosyns.vcmx25_z_patch[p, z] *= (1.0 - max_daily_decay)
                            end
                        end
                    end  # growth check
                end  # C3 check
            else  # no LAI or no LNC: use last valid values
                nrad_p = surfalb.nrad_patch[p]
                for z in 1:nrad_p
                    photosyns.jmx25_z_patch[p, z] = photosyns.jmx25_z_last_valid_patch[p, z]
                    photosyns.vcmx25_z_patch[p, z] = photosyns.vcmx25_z_last_valid_patch[p, z]
                end
            end  # LAI/LNC check
        end  # first day check
    end  # patch loop
    return nothing
end
