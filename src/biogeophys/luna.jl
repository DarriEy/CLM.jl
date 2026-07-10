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

# Isbits scalar view of LunaParamsData for device kernels. LunaParamsData's Vector fields
# can't cross to a GPU, but the LUNA optimization CHAIN reads only these scalar params — so
# the chain functions take this by-value struct instead. Field NAMES match LunaParamsData so
# the chain bodies read `luna_params.<field>` unchanged (only the type annotation changes).
struct _LunaScalars{T}
    cp25_yr2000::T
    kc25_coef::T
    ko25_coef::T
    luna_theta_cj::T
    enzyme_turnover_daily::T
    relhExp::T
    minrelh::T
end
@inline _LunaScalars(lp::LunaParamsData, ::Type{T}) where {T} = _LunaScalars{T}(
    T(lp.cp25_yr2000), T(lp.kc25_coef), T(lp.ko25_coef), T(lp.luna_theta_cj),
    T(lp.enzyme_turnover_daily), T(lp.relhExp), T(lp.minrelh))

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

# Module-level LUNA parameter store (mirrors `const params_inst = PhotoParamsData()`
# in photosynthesis.jl). Populated from the CLM params NetCDF by luna_read_params!.
const luna_params_inst = LunaParamsData()

"""
    luna_read_params!(lp, ds)

Read LUNA parameters (PFT vectors jmaxb0/jmaxb1/wc2wjb0 + scalars) from the open
CLM params NetCDF `ds` into `lp`. PFT vectors are 1-based length MXPFT+1; the file
is 0-based PFT so index p stores file row p (the natural-veg/crop CLM convention).
"""
function luna_read_params!(lp::LunaParamsData, ds)
    luna_params_init!(lp, MXPFT)   # allocate the PFT vectors (NaN)
    rs(name, default) = (haskey(ds, name) && !ismissing(ds[name][1])) ? Float64(ds[name][1]) : default
    # LunaMod readParams converts these from mol/mol to "Luna units" via ×1e5
    # (LunaMod.F90: params_inst%{cp25_yr2000,ko25_coef,kc25_coef} *= 1.e5_r8).
    # Without this, the reference NUE (Kc/Kj at 25°C) is wrong → vcmx25_opt ~0.67×.
    lp.cp25_yr2000          = rs("cp25_yr2000", 42.75e-6) * 1.0e5
    lp.kc25_coef            = rs("kc25_coef", 404.9e-6)   * 1.0e5
    lp.ko25_coef            = rs("ko25_coef", 278.4e-3)   * 1.0e5
    lp.luna_theta_cj        = rs("theta_cj_luna", rs("luna_theta_cj", 0.98))
    lp.enzyme_turnover_daily = rs("enzyme_turnover_daily", 0.0114)
    lp.relhExp              = rs("relhExp", 6.0686)
    lp.minrelh              = rs("minrelh", 0.65)
    for (name, dest) in (("jmaxb0", lp.jmaxb0), ("jmaxb1", lp.jmaxb1), ("wc2wjb0", lp.wc2wjb0))
        haskey(ds, name) || continue
        v = ds[name][:]
        @inbounds for p in 1:min(length(dest), length(v))
            dest[p] = ismissing(v[p]) ? dest[p] : Float64(v[p])
        end
    end
    return nothing
end

# ==========================================================================
# Temperature response functions
# ==========================================================================

"""
    vcmx_t_kattge(tgrow, tleaf) -> Float64

Temperature response for Vcmax, with temperature acclimation (Kattge & Knorr 2007).
"""
function vcmx_t_kattge(tgrow::Real, tleaf::Real)
    T = promote_type(typeof(tgrow), typeof(tleaf))
    tfrz = T(TFRZ); rgas = T(RGAS)
    TlimVcmx = T(668.39) - T(1.07) * smooth_min(smooth_max(tgrow, T(11.0)), T(35.0))
    Vcmxf1 = T(1.0) + exp((TlimVcmx * (T(25.0) + tfrz) - T(200000.0)) / (rgas * (T(25.0) + tfrz)))
    Vcmxf2 = exp((T(72000.0) / (rgas * (T(25.0) + tfrz))) * (T(1.0) - (tfrz + T(25.0)) / (tfrz + tleaf)))
    Vcmxf3 = T(1.0) + exp((TlimVcmx * (tleaf + tfrz) - T(200000.0)) / (rgas * (tleaf + tfrz)))
    return Vcmxf1 * Vcmxf2 / Vcmxf3
end

"""
    jmx_t_kattge(tgrow, tleaf) -> Float64

Temperature response for Jmax, with temperature acclimation (Kattge & Knorr 2007).
"""
function jmx_t_kattge(tgrow::Real, tleaf::Real)
    T = promote_type(typeof(tgrow), typeof(tleaf))
    tfrz = T(TFRZ); rgas = T(RGAS)
    TlimJmx = T(659.7) - T(0.75) * smooth_min(smooth_max(tgrow, T(11.0)), T(35.0))
    Jmxf1 = T(1.0) + exp((TlimJmx * (T(25.0) + tfrz) - T(200000.0)) / (rgas * (T(25.0) + tfrz)))
    Jmxf2 = exp((T(50000.0) / (rgas * (T(25.0) + tfrz))) * (T(1.0) - (tfrz + T(25.0)) / (tleaf + tfrz)))
    Jmxf3 = T(1.0) + exp((TlimJmx * (tleaf + tfrz) - T(200000.0)) / (rgas * (tleaf + tfrz)))
    return Jmxf1 * Jmxf2 / Jmxf3
end

"""
    vcmx_t_leuning(tgrow, tleaf) -> Float64

Temperature response for Vcmax without temperature acclimation (Leuning 2002).
"""
function vcmx_t_leuning(tgrow::Real, tleaf::Real)
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
function jmx_t_leuning(tgrow::Real, tleaf::Real)
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
function resp_t_bernacchi(tleaf::Real)
    # Bernacchi (2001) activation energy 46.39 is in kJ/mol, so RGAS (J/mol/K) must
    # be in kJ/mol/K here (RGAS*1e-3 == Fortran LunaMod's rgas[J/kmol]*1.e-6). Without
    # the 1e-3 the response is ~1.3e8 instead of 1.0 at the 25°C reference.
    T = typeof(tleaf)
    return exp(T(18.72) - T(46.39) / (T(RGAS) * T(1.0e-3) * (tleaf + T(TFRZ))))
end

# ==========================================================================
# Quadratic solver (LUNA-specific, matches Fortran Quadratic subroutine)
# ==========================================================================

"""
    quadratic_luna(a, b, c) -> (r1, r2)

Solve a*x^2 + b*x + c = 0. Matches the Fortran LUNA Quadratic subroutine exactly.
"""
function quadratic_luna(a::Real, b::Real, c::Real)
    T = promote_type(typeof(a), typeof(b), typeof(c))
    r1 = T(1.0e36)
    r2 = T(1.0e36)

    if a == zero(T)
        return (r1, r2)
    end

    if b >= zero(T)
        q = T(-0.5) * (b + sqrt(b * b - T(4.0) * a * c))
    else
        q = T(-0.5) * (b - sqrt(b * b - T(4.0) * a * c))
    end

    r1 = q / a

    if q != zero(T)
        r2 = c / q
    else
        r2 = T(1.0e36)
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
    Kj = smooth_max(ci - c_p, 0.0) / (4.0 * ci + 8.0 * c_p)
    Kc = smooth_max(ci - c_p, 0.0) / (ci + awc)
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
function nue_calc(O2a::Real, ci::Real, tgrow::Real, tleaf::Real,
                  luna_params::_LunaScalars)
    T = promote_type(typeof(O2a), typeof(ci), typeof(tgrow), typeof(tleaf))
    tfrz = T(TFRZ); rgas = T(RGAS)
    Fc = vcmx_t_kattge(tgrow, tleaf) * T(LUNA_Fc25)
    Fj = jmx_t_kattge(tgrow, tleaf) * T(LUNA_Fj25)
    k_c = T(luna_params.kc25_coef) * exp((T(79430.0) / (rgas * (T(25.0) + tfrz))) * (T(1.0) - (tfrz + T(25.0)) / (tfrz + tleaf)))
    k_o = T(luna_params.ko25_coef) * exp((T(36380.0) / (rgas * (T(25.0) + tfrz))) * (T(1.0) - (tfrz + T(25.0)) / (tfrz + tleaf)))
    c_p = T(luna_params.cp25_yr2000) * exp((T(37830.0) / (rgas * (T(25.0) + tfrz))) * (T(1.0) - (tfrz + T(25.0)) / (tfrz + tleaf)))
    awc = k_c * (T(1.0) + O2a / k_o)
    Kj = smooth_max(ci - c_p, zero(T)) / (T(4.0) * ci + T(8.0) * c_p)
    Kc = smooth_max(ci - c_p, zero(T)) / (ci + awc)
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
function photosynthesis_luna!(forc_pbot::Real, tleafd::Real, relh::Real,
                              CO2a::Real, O2a::Real, rb::Real,
                              Vcmax::Real, JmeanL::Real,
                              luna_params::_LunaScalars)
    T = promote_type(typeof(forc_pbot), typeof(tleafd), typeof(relh), typeof(CO2a),
                     typeof(O2a), typeof(rb), typeof(Vcmax), typeof(JmeanL))
    tfrz = T(TFRZ)
    rsmax0 = T(2.0) * T(1.0e4)
    bp = T(2000.0)
    tleaf = tleafd
    tleafk = tleaf + tfrz
    aquad = T(1.0)
    relhc = smooth_max(T(luna_params.minrelh), relh)
    bbb = T(1.0) / bp
    mbb = T(LUNA_mp)
    CO2c = CO2a
    O2c = O2a
    ci = T(0.7) * CO2c
    ciold = ci - T(0.02)
    cf = forc_pbot / (T(8.314) * tleafk) * T(1.0e6)
    gb_mol = cf / rb
    k_c = T(luna_params.kc25_coef) * exp((T(79430.0) / (T(8.314) * (T(25.0) + tfrz))) * (T(1.0) - (tfrz + T(25.0)) / (tfrz + tleaf)))
    k_o = T(luna_params.ko25_coef) * exp((T(36380.0) / (T(8.314) * (T(25.0) + tfrz))) * (T(1.0) - (tfrz + T(25.0)) / (tfrz + tleaf)))
    c_p = T(luna_params.cp25_yr2000) * exp((T(37830.0) / (T(8.314) * (T(25.0) + tfrz))) * (T(1.0) - (tfrz + T(25.0)) / (tfrz + tleaf)))
    awc = k_c * (T(1.0) + O2c / k_o)

    # Initialize variables used across loop scopes
    gs_mol = bbb
    Kc_val = zero(T)
    Wc = zero(T)
    Wj = zero(T)

    # Rubisco limitation iteration
    i = 1
    while abs(ci - ciold) > T(0.01) && i < 100
        i += 1
        ciold = ci
        Kc_val = smooth_max(ci - c_p, zero(T)) / (ci + awc)
        Wc = Kc_val * Vcmax
        gs_mol = bbb + mbb * Wc / CO2c * forc_pbot * relhc
        phi = forc_pbot * (T(1.37) * gs_mol + T(1.6) * gb_mol) / (gb_mol * gs_mol)
        bquad = awc - CO2c + phi * Vcmax
        cquad = -(c_p * phi * Vcmax + awc * CO2c)
        r1, r2 = quadratic_luna(aquad, bquad, cquad)
        ci = smooth_max(r1, r2)
        if ci < zero(T)
            ci = c_p + T(0.5) * ciold
        end
    end

    Kj = smooth_max(ci - c_p, zero(T)) / (T(4.0) * ci + T(8.0) * c_p)
    Kc_val = smooth_max(ci - c_p, zero(T)) / (ci + awc)
    Wc = Kc_val * Vcmax
    Wj = Kj * JmeanL
    ciold = ci - T(0.02)

    # Light limitation iteration (if Wj < Wc)
    if Wj < Wc
        i = 1
        while abs(ci - ciold) > T(0.01) && i < 100
            i += 1
            ciold = ci
            gs_mol = bbb + mbb * Wj / CO2c * forc_pbot * relhc
            phi = forc_pbot * (T(1.37) * gs_mol + T(1.6) * gb_mol) / (gb_mol * gs_mol)
            bquad = T(2.0) * c_p - CO2c + phi * JmeanL / T(4.0)
            cquad = -(c_p * phi * JmeanL / T(4.0) + T(2.0) * c_p * CO2c)
            r1, r2 = quadratic_luna(aquad, bquad, cquad)
            ci = smooth_max(r1, r2)
            if ci < zero(T)
                ci = c_p + T(0.5) * ciold
            end
            Kj = smooth_max(ci - c_p, zero(T)) / (T(4.0) * ci + T(8.0) * c_p)
            Wj = Kj * JmeanL
        end
        Kj = smooth_max(ci - c_p, zero(T)) / (T(4.0) * ci + T(8.0) * c_p)
        Kc_val = smooth_max(ci - c_p, zero(T)) / (ci + awc)
        Wc = Kc_val * Vcmax
        Wj = Kj * JmeanL
    end

    A = (T(1.0) - T(luna_params.luna_theta_cj)) * smooth_max(Wc, Wj) + T(luna_params.luna_theta_cj) * smooth_min(Wc, Wj)
    rs = cf / gs_mol
    rs = smooth_min(rsmax0, rs)

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
function nitrogen_investments!(KcKjFlag::Int, FNCa::Real, Nlc::Real,
                               forc_pbot10::Real, relh10::Real,
                               CO2a10::Real, O2a10::Real,
                               PARi10::Real, PARimx10::Real,
                               rb10::Real, hourpd::Real,
                               tair10::Real, tleafd10::Real, tleafn10::Real,
                               Kj2Kc::Real, JmaxCoef::Real,
                               Fc::Real, Fj::Real,
                               NUEc::Real, NUEj::Real,
                               NUEcref::Real, NUEjref::Real,
                               NUEr::Real, o3coefjmax::Real,
                               jmaxb0::Real, wc2wjb0::Real,
                               Kc_in::Real, Kj_in::Real, ci_in::Real,
                               luna_params::_LunaScalars)
    T = typeof(FNCa)
    leaf_mr_vcm = T(0.015)

    theta = T(0.292) / (T(1.0) + T(0.076) / (Nlc * T(LUNA_Cb)))
    ELTRNabsorb = theta * PARi10
    Jmaxb0act = jmaxb0 * FNCa * Fj

    # o3coefjmax defaults to 1 unless ozone stress_method == 'stress_falk'
    Jmax = Jmaxb0act + JmaxCoef * ELTRNabsorb * o3coefjmax

    JmaxL = theta * PARimx10 / sqrt(T(1.0) + (theta * PARimx10 / Jmax)^T(2.0))
    NUEchg = (NUEc / NUEcref) * (NUEjref / NUEj)
    Wc2Wj = wc2wjb0 * (NUEchg^T(0.5))
    Vcmax = Wc2Wj * JmaxL * Kj2Kc
    JmeanL = theta * PARi10 / sqrt(T(1.0) + (ELTRNabsorb / Jmax)^T(2.0))

    Kc = Kc_in
    Kj = Kj_in
    ci = ci_in

    if KcKjFlag == 0  # update Kc, Kj and ci from photosynthesis
        ci, Kc, Kj, A = photosynthesis_luna!(forc_pbot10, tleafd10, relh10, CO2a10, O2a10, rb10, Vcmax, JmeanL, luna_params)
    else
        Wc = Kc * Vcmax
        Wj = Kj * JmeanL
        A = (T(1.0) - T(luna_params.luna_theta_cj)) * smooth_max(Wc, Wj) + T(luna_params.luna_theta_cj) * smooth_min(Wc, Wj)
    end

    PSN = T(LUNA_Cv) * A * hourpd
    Vcmaxnight = vcmx_t_kattge(tair10, tleafn10) / vcmx_t_kattge(tair10, tleafd10) * Vcmax
    RESP = T(LUNA_Cv) * leaf_mr_vcm * (Vcmax * hourpd + Vcmaxnight * (T(24.0) - hourpd))
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
function nitrogen_allocation!(FNCa::Real, forc_pbot10::Real, relh10::Real,
                              CO2a10::Real, O2a10::Real,
                              PARi10::Real, PARimx10::Real,
                              rb10::Real, hourpd::Real,
                              tair10::Real, tleafd10::Real, tleafn10::Real,
                              jmaxb0_val::Real, jmaxb1_val::Real, wc2wjb0_val::Real,
                              PNlcold::Real, PNetold::Real,
                              PNrespold::Real, PNcbold::Real,
                              dayl_factor::Real, o3coefjmax::Real,
                              NUEjref::Real, NUEcref::Real,
                              luna_params::_LunaScalars)

    T = typeof(FNCa)
    Nlc   = PNlcold * FNCa
    Net   = PNetold * FNCa
    Nresp = PNrespold * FNCa
    Ncb   = PNcbold * FNCa
    if Nlc > FNCa * T(0.5)
        Nlc = T(0.5) * FNCa
    end
    chg_per_step = T(0.02) * FNCa
    PNlc = PNlcold
    PNlcoldi = PNlcold - T(0.001)
    PARi10c  = smooth_max(T(LUNA_PARLowLim), PARi10)
    PARimx10c = smooth_max(T(LUNA_PARLowLim), PARimx10)
    increase_flag = 0
    jj = 1
    tleafd10c = smooth_min(smooth_max(tleafd10, T(LUNA_Trange1)), T(LUNA_Trange2))
    tleafn10c = smooth_min(smooth_max(tleafn10, T(LUNA_Trange1)), T(LUNA_Trange2))
    ci = T(0.7) * CO2a10
    JmaxCoef = jmaxb1_val * dayl_factor * (T(1.0) - exp(-T(luna_params.relhExp) * smooth_max(relh10 -
        T(luna_params.minrelh), zero(T)) / (T(1.0) - T(luna_params.minrelh))))

    # Initialize Kc, Kj and loop variables for proper scoping
    Kc = zero(T)
    Kj = zero(T)
    Nstore = zero(T)
    Npsntarget = zero(T)
    PSN = zero(T)

    while PNlcoldi != PNlc && jj < 100
        Fc = vcmx_t_kattge(tair10, tleafd10c) * T(LUNA_Fc25)
        Fj = jmx_t_kattge(tair10, tleafd10c) * T(LUNA_Fj25)
        NUEr = T(LUNA_Cv) * T(LUNA_NUEr25) * (resp_t_bernacchi(tleafd10c) * hourpd +
            resp_t_bernacchi(tleafn10c) * (T(24.0) - hourpd))

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
        if Nstore > zero(T) && (increase_flag == 1 || jj == 1)
            Nlc2 = Nlc + chg_per_step
            if Nlc2 / FNCa > T(0.95)
                Nlc2 = T(0.95) * FNCa
            end
            KcKjFlag = 1
            _, _, _, _, Net2, Ncb2, Nresp2, PSN2, RESP2, _, _, _ =
                nitrogen_investments!(KcKjFlag, FNCa, Nlc2, forc_pbot10, relh10, CO2a10, O2a10,
                    PARi10c, PARimx10c, rb10, hourpd, tair10, tleafd10c, tleafn10c,
                    Kj2Kc, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, o3coefjmax,
                    jmaxb0_val, wc2wjb0_val, Kc, Kj, ci, luna_params)

            Npsntarget2 = Nlc2 + Ncb2 + Net2
            Carboncost2 = (Npsntarget2 - Npsntarget) * T(LUNA_NMCp25) * T(LUNA_Cv) *
                (resp_t_bernacchi(tleafd10c) * hourpd + resp_t_bernacchi(tleafn10c) * (T(24.0) - hourpd))
            Carbongain2 = PSN2 - PSN
            if Carbongain2 > Carboncost2 && (Npsntarget2 + Nresp2 < T(0.95) * FNCa)
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
            if Nstore < zero(T)
                Nlc1 = Nlc * T(0.8)  # bigger step of decrease if negative
            else
                Nlc1 = Nlc - chg_per_step
            end
            if Nlc1 < T(0.05)
                Nlc1 = T(0.05)
            end
            KcKjFlag = 1
            _, _, _, _, Net1, Ncb1, Nresp1, PSN1, RESP1, _, _, _ =
                nitrogen_investments!(KcKjFlag, FNCa, Nlc1, forc_pbot10, relh10, CO2a10, O2a10,
                    PARi10c, PARimx10c, rb10, hourpd, tair10, tleafd10c, tleafn10c,
                    Kj2Kc, JmaxCoef, Fc, Fj, NUEc, NUEj, NUEcref, NUEjref, NUEr, o3coefjmax,
                    jmaxb0_val, wc2wjb0_val, Kc, Kj, ci, luna_params)

            Npsntarget1 = Nlc1 + Ncb1 + Net1
            Carboncost1 = (Npsntarget - Npsntarget1) * T(LUNA_NMCp25) * T(LUNA_Cv) *
                (resp_t_bernacchi(tleafd10c) * hourpd + resp_t_bernacchi(tleafn10c) * (T(24.0) - hourpd))
            Carbongain1 = PSN - PSN1
            if (Carbongain1 < Carboncost1 && Nlc1 > T(0.05)) || (Npsntarget + Nresp) > T(0.95) * FNCa
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
# One thread per patch; zeroes the 24 hr accumulators (byte-identical to the host loop
# on the KA CPU backend). par24d_z/par24x_z rows zeroed by an in-thread z-loop.
@kernel function _luna_clear24_kernel!(t_veg_day, t_veg_night, nnightsteps, ndaysteps,
        fpsn24, par24d_z, par24x_z, @Const(mask), lo::Int, hi::Int, ncan::Int)
    p = @index(Global)
    @inbounds if lo <= p <= hi && mask[p]
        T = eltype(t_veg_day)
        t_veg_day[p]   = zero(T)
        t_veg_night[p] = zero(T)
        fpsn24[p]      = zero(T)
        nnightsteps[p] = 0
        ndaysteps[p]   = 0
        for z in 1:ncan
            par24d_z[p, z] = zero(T)
            par24x_z[p, z] = zero(T)
        end
    end
end

function clear24_climate_luna!(solarabs::SolarAbsorbedData,
                               photosyns::PhotosynthesisData,
                               temperature::TemperatureData,
                               patchdata::PatchData,
                               mask_patch::AbstractVector{Bool},
                               bounds::UnitRange{Int})
    isempty(bounds) && return nothing
    _launch!(_luna_clear24_kernel!, temperature.t_veg_day_patch, temperature.t_veg_night_patch,
        temperature.nnightsteps_patch, temperature.ndaysteps_patch, photosyns.fpsn24_patch,
        solarabs.par24d_z_patch, solarabs.par24x_z_patch, mask_patch,
        first(bounds), last(bounds), size(solarabs.par24d_z_patch, 2))
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
                             mask_patch::AbstractVector{Bool},
                             bounds::UnitRange{Int},
                             dtime::Real)
    isempty(bounds) && return nothing
    FT = eltype(temperature.t_veg_patch)
    _launch!(_luna_acc24_kernel!, temperature.t_veg_day_patch, temperature.t_veg_night_patch,
        temperature.ndaysteps_patch, temperature.nnightsteps_patch, temperature.t_veg_patch,
        solarabs.sabv_patch, solarabs.parsun_z_patch, solarabs.par24d_z_patch,
        solarabs.par24x_z_patch, photosyns.fpsn24_patch, photosyns.fpsn_patch,
        canopystate.laisun_z_patch, canopystate.laisha_z_patch, surfalb.nrad_patch,
        mask_patch, first(bounds), last(bounds), FT(dtime), FT(SPVAL))
    return nothing
end

# One thread per patch (byte-identical to the host loop on KA CPU). par24d_z[p,z]+= has
# z as the loop var → distinct element per iter (safe); fpsn24 accumulates once per patch.
@kernel function _luna_acc24_kernel!(t_veg_day, t_veg_night, ndaysteps, nnightsteps, t_veg,
        @Const(sabv), @Const(parsun_z), par24d_z, par24x_z, fpsn24, @Const(fpsn),
        @Const(laisun_z), @Const(laisha_z), @Const(nrad), @Const(mask),
        lo::Int, hi::Int, dtime, spval)
    p = @index(Global)
    @inbounds if lo <= p <= hi && mask[p]
        T = eltype(t_veg_day)
        if t_veg_day[p] != spval
            if sabv[p] > zero(T)
                t_veg_day[p] += t_veg[p]
                ndaysteps[p] += 1
            else
                t_veg_night[p] += t_veg[p]
                nnightsteps[p] += 1
            end

            nrad_p = nrad[p]
            for z in 1:nrad_p
                tlaii = laisun_z[p, z] + laisha_z[p, z]
                if tlaii > zero(T)
                    TRad = parsun_z[p, z]
                    par24d_z[p, z] += dtime * TRad
                    if TRad > par24x_z[p, z]
                        par24x_z[p, z] = TRad
                    end
                end
            end

            fpsn24[p] += dtime * fpsn[p]
        end
    end
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
                              mask_patch::AbstractVector{Bool},
                              bounds::UnitRange{Int},
                              oair::AbstractVector{<:Real},
                              cair::AbstractVector{<:Real},
                              rb::AbstractVector{<:Real},
                              rh::AbstractVector{<:Real},
                              dtime::Real)
    isempty(bounds) && return nothing
    FT = eltype(temperature.t_veg_day_patch)
    _launch!(_luna_acc240_kernel!, temperature.t_veg_day_patch, temperature.t_veg_night_patch,
        temperature.ndaysteps_patch, temperature.nnightsteps_patch,
        temperature.t_veg10_day_patch, temperature.t_veg10_night_patch,
        solarabs.par24d_z_patch, solarabs.par24x_z_patch,
        solarabs.par240d_z_patch, solarabs.par240x_z_patch,
        waterdiag.rh10_af_patch, frictionvel.rb10_patch, _to_backend_like(temperature.t_veg_day_patch, FT, rb),
        _to_backend_like(temperature.t_veg_day_patch, FT, rh),
        mask_patch, first(bounds), last(bounds), FT(dtime), FT(SPVAL),
        size(solarabs.par24d_z_patch, 2))
    return nothing
end

# One thread per patch (byte-identical to the host loop on KA CPU). The par24d_z_i slice
# temp is inlined per-z (used only at the same z); par240*_z[p,z] writes use z as the loop
# var → distinct element per iter.
@kernel function _luna_acc240_kernel!(t_veg_day, t_veg_night, ndaysteps, nnightsteps,
        t_veg10_day, t_veg10_night, @Const(par24d_z), @Const(par24x_z),
        par240d_z, par240x_z, rh10_af, rb10, @Const(rb), @Const(rh),
        @Const(mask), lo::Int, hi::Int, dtime, spval, ncan::Int)
    p = @index(Global)
    @inbounds if lo <= p <= hi && mask[p]
        T = eltype(t_veg_day)
        if t_veg_day[p] != spval
            ndaysteps_p = ndaysteps[p]
            first_day = par240d_z[p, 1] == spval
            for z in 1:ncan
                par24d_z_i = ndaysteps_p > 0 ? par24d_z[p, z] / (dtime * ndaysteps_p) : zero(T)
                if first_day
                    par240x_z[p, z] = par24x_z[p, z]
                    par240d_z[p, z] = par24d_z_i
                else
                    par240x_z[p, z] = T(0.9) * par240x_z[p, z] + T(0.1) * par24x_z[p, z]
                    par240d_z[p, z] = T(0.9) * par240d_z[p, z] + T(0.1) * par24d_z_i
                end
            end

            # 10-day running mean daytime temperature
            if ndaysteps_p > 0
                t_veg_dayi = t_veg_day[p] / ndaysteps_p
            else
                t_veg_dayi = t_veg_night[p] / nnightsteps[p]
            end
            if t_veg10_day[p] == spval
                t_veg10_day[p] = t_veg_dayi
            end
            t_veg10_day[p] = T(0.9) * t_veg10_day[p] + T(0.1) * t_veg_dayi

            # 10-day running mean nighttime temperature
            nnightsteps_p = nnightsteps[p]
            if nnightsteps_p > 0
                t_veg_nighti = t_veg_night[p] / nnightsteps_p
            else
                t_veg_nighti = t_veg_day[p] / ndaysteps_p
            end
            if t_veg10_night[p] == spval
                t_veg10_night[p] = t_veg_nighti
            end
            t_veg10_night[p] = T(0.9) * t_veg10_night[p] + T(0.1) * t_veg_nighti

            # 10-day running mean relative humidity
            if rh10_af[p] == spval
                rh10_af[p] = rh[p]
            end
            rh10_af[p] = T(0.9) * rh10_af[p] + T(0.1) * smooth_min(one(T), rh[p])

            # 10-day running mean boundary layer resistance
            if rb10[p] == spval
                rb10[p] = rb[p]
            end
            rb10[p] = T(0.9) * rb10[p] + T(0.1) * rb[p]
        end
    end
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
# ==========================================================================
# GPU kernelization of the LUNA end-of-day capacity update.
# One thread per patch runs the whole per-patch Rubisco-N optimization (the
# bounded nitrogen_allocation! fixed-point + the z-level passes). The chain it
# calls is already T-generic + allocation-free + device-callable; per-thread
# FNCa scratch is a (np x nlevcan) backend array (row p). The ~30 state/input
# arrays across 9 structs are grouped into device-view bundles to stay under
# Metal's arg limit; the sanity @warn/error checks become pure clamps (device-
# callable, byte-identical on any path where they don't fire). luna_params'
# Vector fields can't cross to a GPU, so its scalars ride in an isbits
# _LunaScalars and its 3 pft-vectors are moved onto the backend in _LunaPft.
# ==========================================================================
struct _LunaPsn{V,M}
    vcmx25_z::M; jmx25_z::M; vcmx25_z_last::M; jmx25_z_last::M
    pnlc_z::M; enzs_z::M; lnca::V; fpsn24::V
end
Adapt.@adapt_structure _LunaPsn

struct _LunaClim{V}
    t_veg_day::V; t_veg10_day::V; t_veg10_night::V; t_a10::V
end
Adapt.@adapt_structure _LunaClim

struct _LunaRad{V,VI,M}
    tlai::V; nrad::VI; tlai_z::M; par240d_z::M; par240x_z::M; rh10_af::V; rb10::V
end
Adapt.@adapt_structure _LunaRad

struct _LunaIdx{V,VI}
    itype::VI; gridcell::VI; dayl::V; dayl_factor::V
    forc_pbot10::V; CO2::V; O2::V; o3coefjmax::V
end
Adapt.@adapt_structure _LunaIdx

struct _LunaPft{V,M}
    c3psn::V; slatop::V; rhol::M; taul::M; jmaxb0::V; jmaxb1::V; wc2wjb0::V
end
Adapt.@adapt_structure _LunaPft

# Per-patch LUNA capacity update (device-callable; one thread owns patch p).
@inline function _luna_patch_body!(p, psn, clim, rad, idx, pft, lp, NUEjref, NUEcref, fnca)
    T = eltype(psn.vcmx25_z)
    @inbounds begin
        ft = idx.itype[p] + 1
        g  = idx.gridcell[p]

        if clim.t_veg_day[p] != T(SPVAL)
            CO2a10   = idx.CO2[p]
            O2a10    = idx.O2[p]
            hourpd   = idx.dayl[g] / T(3600.0)
            tleafd10 = clim.t_veg10_day[p] - T(TFRZ)
            tleafn10 = clim.t_veg10_night[p] - T(TFRZ)
            tleaf10  = (idx.dayl[g] * tleafd10 + (T(86400.0) - idx.dayl[g]) * tleafn10) / T(86400.0)
            tair10   = clim.t_a10[p] - T(TFRZ)
            relh10   = smooth_min(T(1.0), rad.rh10_af[p])
            rb10v    = rad.rb10[p]

            EnzTurnoverTFactor = T(LUNA_Q10Enz)^(T(0.1) * (smooth_min(T(40.0), tleaf10) - T(25.0)))
            max_daily_pchg = EnzTurnoverTFactor * lp.enzyme_turnover_daily

            rabsorb = T(1.0) - pft.rhol[ft, 1] - pft.taul[ft, 1]

            if rad.tlai[p] > T(0.0) && psn.lnca[p] > T(0.0)
                RadTop = rad.par240d_z[p, 1] / rabsorb
                PARTop = RadTop * T(4.6)

                if unsafe_trunc(Int, round(pft.c3psn[ft])) == 1   # round(Int,·) boxes on Metal
                    if psn.fpsn24[p] > T(0.0)
                        nrad_p = rad.nrad[p]

                        for z in 1:nrad_p
                            if rad.tlai_z[p, z] > T(0.0)
                                qabs = rad.par240d_z[p, z] / rabsorb
                                PARi10 = qabs * T(4.6)
                            else
                                PARi10 = T(0.01)
                            end
                            relRad = PARi10 / PARTop
                            relCLNCa = T(0.1802) * log(relRad) + T(1.0)
                            relCLNCa = smooth_max(T(0.2), relCLNCa)
                            relCLNCa = smooth_min(T(1.0), relCLNCa)

                            SNCa = T(1.0) / pft.slatop[ft] * T(LUNA_SNC)
                            lnc_p = psn.lnca[p]
                            if T(0.9) * lnc_p > SNCa
                                fnca[p, z] = relCLNCa * (lnc_p - SNCa)
                            else
                                fnca[p, z] = relCLNCa * T(0.1) * lnc_p
                            end
                        end

                        for z in 1:nrad_p
                            FNCa = fnca[p, z]
                            if FNCa > T(15.0)   # boundary clamp (was @warn)
                                FNCa = T(15.0)
                            end

                            radmax2mean = rad.par240x_z[p, z] / rad.par240d_z[p, z]
                            if rad.tlai_z[p, z] > T(0.0)
                                qabs = rad.par240d_z[p, z] / rabsorb
                                PARi10 = qabs * T(4.6)
                            else
                                PARi10 = T(0.01)
                            end
                            PARimx10 = PARi10 * radmax2mean

                            PNlcold = psn.pnlc_z[p, z]
                            PNetold = T(0.0)
                            PNrespold = T(0.0)
                            PNcbold = T(0.0)

                            PNstoreopt, PNlcopt, PNetopt, PNrespopt, PNcbopt =
                                nitrogen_allocation!(FNCa, idx.forc_pbot10[p], relh10, CO2a10, O2a10,
                                    PARi10, PARimx10, rb10v, hourpd, tair10, tleafd10, tleafn10,
                                    pft.jmaxb0[ft], pft.jmaxb1[ft], pft.wc2wjb0[ft],
                                    PNlcold, PNetold, PNrespold, PNcbold, idx.dayl_factor[p],
                                    idx.o3coefjmax[p], NUEjref, NUEcref, lp)

                            vcmx25_opt = PNcbopt * FNCa * T(LUNA_Fc25)
                            jmx25_opt  = PNetopt * FNCa * T(LUNA_Fj25)

                            chg = vcmx25_opt - psn.vcmx25_z[p, z]
                            chg_constrn = smooth_min(smooth_abs(chg), psn.vcmx25_z[p, z] * max_daily_pchg)
                            psn.vcmx25_z[p, z] += sign(chg) * chg_constrn
                            psn.vcmx25_z_last[p, z] = psn.vcmx25_z[p, z]

                            chg = jmx25_opt - psn.jmx25_z[p, z]
                            chg_constrn = smooth_min(smooth_abs(chg), psn.jmx25_z[p, z] * max_daily_pchg)
                            psn.jmx25_z[p, z] += sign(chg) * chg_constrn
                            psn.jmx25_z_last[p, z] = psn.jmx25_z[p, z]

                            psn.pnlc_z[p, z] = PNlcopt

                            if psn.enzs_z[p, z] < T(1.0)
                                psn.enzs_z[p, z] *= (T(1.0) + max_daily_pchg)
                            end

                            # Sanity clamps (were @warn/error; device-callable, byte-identical
                            # on any path where the abnormal branch does not fire).
                            if isnan(psn.vcmx25_z[p, z])
                                psn.vcmx25_z[p, z] = T(50.0)
                            end
                            if psn.vcmx25_z[p, z] > T(1000.0) || psn.vcmx25_z[p, z] < T(0.0)
                                psn.vcmx25_z[p, z] = T(50.0)
                            end
                            if isnan(psn.jmx25_z[p, z])
                                psn.jmx25_z[p, z] = T(85.0)
                            end
                            if psn.jmx25_z[p, z] > T(2000.0) || psn.jmx25_z[p, z] < T(0.0)
                                psn.jmx25_z[p, z] = T(85.0)
                            end
                        end  # z loop (growth)
                    else  # decay during drought or winter
                        max_daily_decay = smooth_min(T(0.5), T(0.1) * max_daily_pchg)
                        nrad_p = rad.nrad[p]
                        for z in 1:nrad_p
                            if psn.enzs_z[p, z] > T(0.5)
                                psn.enzs_z[p, z] *= (T(1.0) - max_daily_decay)
                                psn.jmx25_z[p, z] *= (T(1.0) - max_daily_decay)
                                psn.vcmx25_z[p, z] *= (T(1.0) - max_daily_decay)
                            end
                        end
                    end  # growth check
                end  # C3 check
            else  # no LAI or no LNC: use last valid values
                nrad_p = rad.nrad[p]
                for z in 1:nrad_p
                    psn.jmx25_z[p, z] = psn.jmx25_z_last[p, z]
                    psn.vcmx25_z[p, z] = psn.vcmx25_z_last[p, z]
                end
            end  # LAI/LNC check
        end  # first day check
    end
    return nothing
end

@kernel function _luna_capacity_kernel!(psn, clim, rad, idx, pft, lp, @Const(mask),
        lo::Int, hi::Int, NUEjref, NUEcref, fnca)
    p = @index(Global)
    @inbounds if lo <= p <= hi && mask[p]
        _luna_patch_body!(p, psn, clim, rad, idx, pft, lp, NUEjref, NUEcref, fnca)
    end
end

function update_photosynthesis_capacity!(photosyns::PhotosynthesisData,
                                         temperature::TemperatureData,
                                         canopystate::CanopyStateData,
                                         surfalb::SurfaceAlbedoData,
                                         solarabs::SolarAbsorbedData,
                                         waterdiag::WaterDiagnosticBulkData,
                                         frictionvel::FrictionVelocityData,
                                         patchdata::PatchData,
                                         gridcell::GridcellData,
                                         mask_patch::AbstractVector{Bool},
                                         bounds::UnitRange{Int},
                                         dayl_factor::AbstractVector{<:Real},
                                         forc_pbot10::AbstractVector{<:Real},
                                         CO2_p240::AbstractVector{<:Real},
                                         O2_p240::AbstractVector{<:Real},
                                         c3psn_pft::AbstractVector{<:Real},
                                         slatop_pft::AbstractVector{<:Real},
                                         leafcn_pft::AbstractVector{<:Real},
                                         rhol_pft::AbstractMatrix{<:Real},
                                         taul_pft::AbstractMatrix{<:Real},
                                         o3coefjmax::AbstractVector{<:Real},
                                         luna_params::LunaParamsData,
                                         dtime::Real,
                                         nlevcan_val::Int)
    # FT follows the STATE (photosyns), not the loose forcing args — on a device run the
    # state is Float32 while the forcing/pft args may arrive as host Float64 (moved to the
    # backend below); on host both are Float64. Wrong FT here -> MtlArray{Float64} (Metal-fatal).
    FT = eltype(photosyns.vcmx25_z_patch)
    # Reference NUE is patch-invariant (host-side); the chain reads only luna_params
    # SCALARS -> carry them in an isbits _LunaScalars (its Vector fields can't cross to a GPU).
    NUEjref, NUEcref, _Kj2Kcref = nue_ref(luna_params)
    lp = _LunaScalars(luna_params, FT)

    _ref = photosyns.vcmx25_z_patch
    psn = _LunaPsn(photosyns.vcmx25_z_patch, photosyns.jmx25_z_patch,
        photosyns.vcmx25_z_last_valid_patch, photosyns.jmx25_z_last_valid_patch,
        photosyns.pnlc_z_patch, photosyns.enzs_z_patch,
        photosyns.lnca_patch, photosyns.fpsn24_patch)
    clim = _LunaClim(temperature.t_veg_day_patch, temperature.t_veg10_day_patch,
        temperature.t_veg10_night_patch, temperature.t_a10_patch)
    rad = _LunaRad(canopystate.tlai_patch, surfalb.nrad_patch, surfalb.tlai_z_patch,
        solarabs.par240d_z_patch, solarabs.par240x_z_patch,
        waterdiag.rh10_af_patch, frictionvel.rb10_patch)
    # loose forcing + host pft params moved onto the state backend (no-op on host).
    idx = _LunaIdx(patchdata.itype, patchdata.gridcell, gridcell.dayl,
        _to_backend_like(_ref, FT, dayl_factor), _to_backend_like(_ref, FT, forc_pbot10),
        _to_backend_like(_ref, FT, CO2_p240), _to_backend_like(_ref, FT, O2_p240),
        _to_backend_like(_ref, FT, o3coefjmax))
    pft = _LunaPft(_to_backend_like(_ref, FT, c3psn_pft), _to_backend_like(_ref, FT, slatop_pft),
        _to_backend_like(_ref, FT, rhol_pft), _to_backend_like(_ref, FT, taul_pft),
        _to_backend_like(_ref, FT, luna_params.jmaxb0), _to_backend_like(_ref, FT, luna_params.jmaxb1),
        _to_backend_like(_ref, FT, luna_params.wc2wjb0))

    # Per-thread FNCa scratch (row p), allocated once on the state backend.
    np = size(photosyns.vcmx25_z_patch, 1)
    fnca = fill!(similar(_ref, FT, np, max(nlevcan_val, 1)), zero(FT))

    let be = _kernel_backend(_ref)
        lo = first(bounds); hi = last(bounds)
        if hi >= lo
            _luna_capacity_kernel!(be)(psn, clim, rad, idx, pft, lp, mask_patch,
                lo, hi, FT(NUEjref), FT(NUEcref), fnca; ndrange = hi)
            KA.synchronize(be)
        end
    end
    return nothing
end
