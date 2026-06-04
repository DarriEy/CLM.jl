# ==========================================================================
# Ported from: src/biogeochem/VOCEmissionMod.F90
# Volatile organic compound emission module (MEGAN v2.1)
#
# Calculates VOC emissions following MEGAN (Model of Emissions of Gases
# and Aerosols from Nature) v2.1 for 20 compound classes.
# E = epsilon * gamma * rho
# where epsilon = baseline emission factors, gamma = activity factors
# (light, temperature, leaf age, LAI, soil moisture, CO2), rho = 1.
#
# Public types:
#   MEGANFactors        — MEGAN compound class parameters (Agro, Amat, etc.)
#   MEGANCompound       — Single MEGAN compound descriptor
#   MEGANMechComp       — Mechanism compound mapping
#   VOCEmisData         — VOC emission state/flux data
#
# Public functions:
#   vocemis_init!        — Initialize VOC emission data arrays
#   vocemis_clean!       — Deallocate VOC emission data arrays
#   voc_emission!        — Main VOC emission driver
#   get_map_EF           — Mapped emission factor for isoprene
#   get_gamma_P          — Activity factor for PPFD
#   get_gamma_L          — Activity factor for LAI
#   get_gamma_SM         — Activity factor for soil moisture
#   get_gamma_T          — Activity factor for temperature
#   get_gamma_A          — Activity factor for leaf age
#   get_gamma_C          — Activity factor for CO2 (isoprene only)
# ==========================================================================

# ---------------------------------------------------------------------------
# MEGAN emission factors and compound class parameters
# Ported from MEGANFactorsMod.F90
# ---------------------------------------------------------------------------

"""
    MEGANFactors

MEGAN compound class parameters indexed by class number.
Contains activity factor parameters for leaf age, temperature,
and light dependence.
Ported from `MEGANFactorsMod` in CLM Fortran.
"""
# Reparametrized {FT,V} + @adapt_structure so the per-class factor vectors are
# device-movable (get_gamma_A/get_gamma_T index them inside the Metal kernel).
Base.@kwdef mutable struct MEGANFactors{FT<:Real, V<:AbstractVector{FT}}
    n_classes::Int = 20
    # Leaf age fractions (indexed by class number, 1:n_classes)
    Agro ::V = Float64[]  # growing leaves factor
    Amat ::V = Float64[]  # mature leaves factor
    Anew ::V = Float64[]  # new leaves factor
    Aold ::V = Float64[]  # old leaves factor
    # Temperature parameters
    betaT::V = Float64[]  # beta for T-dependent LIF
    ct1  ::V = Float64[]  # ct1 coefficient for T activity factor
    ct2  ::V = Float64[]  # ct2 coefficient for T activity factor
    Ceo  ::V = Float64[]  # Ceo coefficient for Eopt
    # Light dependence fraction
    LDF  ::V = Float64[]  # light-dependent fraction (0-1)
end

MEGANFactors{FT}(; kwargs...) where {FT<:Real} =
    MEGANFactors{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure MEGANFactors

"""
    MEGANCompound

Descriptor for a single MEGAN mega-compound.
Ported from `shr_megan_megcomp_t` in `shr_megan_mod`.
"""
Base.@kwdef mutable struct MEGANCompound{FT<:Real}
    name          ::String          = ""        # compound name (e.g. "isoprene")
    index         ::Int             = 0         # index in mega-compound list
    class_number  ::Int             = 0         # MEGAN class number (1-20)
    molec_weight  ::FT         = 0.0       # molecular weight [g/mol]
    coeff         ::FT         = 1.0       # coefficient for mapping to mechanism compound
    emis_factors  ::Vector{FT} = Float64[] # emission factors per PFT [ug m-2 h-1]
end

"""
    MEGANMechComp

Mechanism compound mapping: maps one or more MEGAN compounds
to a single chemical mechanism compound.
Ported from `shr_megan_mechcomps` in `shr_megan_mod`.
"""
Base.@kwdef mutable struct MEGANMechComp
    name          ::String       = ""   # mechanism compound name
    n_megan_comps ::Int          = 0    # number of MEGAN compounds mapped
    megan_indices ::Vector{Int}  = Int[] # indices into MEGANCompound array
end

# ---------------------------------------------------------------------------
# VOCEmisData — VOC emission state/flux data
# ---------------------------------------------------------------------------

"""
    VOCEmisData

VOC emission state and flux data.
Contains diagnostic gamma factors and fluxes at the patch level.
Ported from `vocemis_type` in `VOCEmissionMod.F90`.
"""
Base.@kwdef mutable struct VOCEmisData{FT<:Real,
                           V<:AbstractVector{FT},
                           M<:AbstractMatrix{FT}}
    # Diagnostic coefficients (patch level)
    Eopt_out_patch    ::V = Float64[]  # Eopt coefficient
    topt_out_patch    ::V = Float64[]  # topt coefficient
    alpha_out_patch   ::V = Float64[]  # alpha coefficient
    cp_out_patch      ::V = Float64[]  # cp coefficient

    # PAR diagnostics (patch level)
    paru_out_patch    ::V = Float64[]  # sunlit PAR [umol/m2/s]
    par24u_out_patch  ::V = Float64[]  # sunlit PAR 24hr avg
    par240u_out_patch ::V = Float64[]  # sunlit PAR 240hr avg
    para_out_patch    ::V = Float64[]  # shade PAR [umol/m2/s]
    par24a_out_patch  ::V = Float64[]  # shade PAR 24hr avg
    par240a_out_patch ::V = Float64[]  # shade PAR 240hr avg

    # Gamma diagnostics (patch level)
    gamma_out_patch   ::V = Float64[]  # total gamma
    gammaL_out_patch  ::V = Float64[]  # gamma for LAI
    gammaT_out_patch  ::V = Float64[]  # gamma for temperature
    gammaP_out_patch  ::V = Float64[]  # gamma for PPFD
    gammaA_out_patch  ::V = Float64[]  # gamma for leaf age
    gammaS_out_patch  ::V = Float64[]  # gamma for soil moisture
    gammaC_out_patch  ::V = Float64[]  # gamma for CO2

    # Fluxes (patch level)
    vocflx_tot_patch  ::V = Float64[]  # total VOC flux [moles/m2/sec]
    vocflx_patch      ::M = Matrix{Float64}(undef, 0, 0)  # per-mechanism-compound flux [moles/m2/sec]

    # Gridcell isoprene emission factors (6 × ngrc)
    efisop_grc        ::M = Matrix{Float64}(undef, 0, 0)  # [ug m-2 h-1]

    # Per-mega-compound flux output (n_megcomps × np) [kg/m2/sec] for history
    meg_flux_out      ::M = Matrix{Float64}(undef, 0, 0)
end

VOCEmisData{FT}(; kwargs...) where {FT<:Real} =
    VOCEmisData{FT, Vector{FT}, Matrix{FT}}(; kwargs...)
Adapt.@adapt_structure VOCEmisData


# ---------------------------------------------------------------------------
# vocemis_init! — Initialize VOC emission data
# ---------------------------------------------------------------------------

"""
    vocemis_init!(voc, np, ng, n_megcomps, n_mechcomps)

Allocate and initialize VOC emission data arrays.
Ported from `InitAllocate` in `VOCEmissionMod.F90`.
"""
function vocemis_init!(voc::VOCEmisData, np::Int, ng::Int,
                       n_megcomps::Int, n_mechcomps::Int)
    # Diagnostic coefficients
    voc.Eopt_out_patch    = fill(NaN, np)
    voc.topt_out_patch    = fill(NaN, np)
    voc.alpha_out_patch   = fill(NaN, np)
    voc.cp_out_patch      = fill(NaN, np)

    # PAR diagnostics
    voc.paru_out_patch    = fill(NaN, np)
    voc.par24u_out_patch  = fill(NaN, np)
    voc.par240u_out_patch = fill(NaN, np)
    voc.para_out_patch    = fill(NaN, np)
    voc.par24a_out_patch  = fill(NaN, np)
    voc.par240a_out_patch = fill(NaN, np)

    # Gamma diagnostics
    voc.gamma_out_patch   = fill(NaN, np)
    voc.gammaL_out_patch  = fill(NaN, np)
    voc.gammaT_out_patch  = fill(NaN, np)
    voc.gammaP_out_patch  = fill(NaN, np)
    voc.gammaA_out_patch  = fill(NaN, np)
    voc.gammaS_out_patch  = fill(NaN, np)
    voc.gammaC_out_patch  = fill(NaN, np)

    # Fluxes
    voc.vocflx_tot_patch  = fill(NaN, np)
    voc.vocflx_patch      = fill(NaN, np, max(n_mechcomps, 1))

    # Gridcell isoprene emission factors
    voc.efisop_grc        = fill(NaN, 6, ng)

    # Per-mega-compound flux output
    voc.meg_flux_out      = zeros(max(n_megcomps, 1), np)

    return nothing
end

# ---------------------------------------------------------------------------
# vocemis_clean! — Deallocate VOC emission data
# ---------------------------------------------------------------------------

"""
    vocemis_clean!(voc)

Deallocate VOC emission data arrays.
"""
function vocemis_clean!(voc::VOCEmisData{FT}) where {FT}
    voc.Eopt_out_patch    = FT[]
    voc.topt_out_patch    = FT[]
    voc.alpha_out_patch   = FT[]
    voc.cp_out_patch      = FT[]
    voc.paru_out_patch    = FT[]
    voc.par24u_out_patch  = FT[]
    voc.par240u_out_patch = FT[]
    voc.para_out_patch    = FT[]
    voc.par24a_out_patch  = FT[]
    voc.par240a_out_patch = FT[]
    voc.gamma_out_patch   = FT[]
    voc.gammaL_out_patch  = FT[]
    voc.gammaT_out_patch  = FT[]
    voc.gammaP_out_patch  = FT[]
    voc.gammaA_out_patch  = FT[]
    voc.gammaS_out_patch  = FT[]
    voc.gammaC_out_patch  = FT[]
    voc.vocflx_tot_patch  = FT[]
    voc.vocflx_patch      = Matrix{FT}(undef, 0, 0)
    voc.efisop_grc        = Matrix{FT}(undef, 0, 0)
    voc.meg_flux_out      = Matrix{FT}(undef, 0, 0)
    return nothing
end

# ---------------------------------------------------------------------------
# megan_factors_init! — Initialize default MEGAN factors
# ---------------------------------------------------------------------------

"""
    megan_factors_init!(mf, n_classes)

Initialize MEGAN factors with default values for `n_classes` compound classes.
Ported from `megan_factors_init` in `MEGANFactorsMod.F90`.
Default values follow MEGAN v2.1 (Guenther et al., 2012).
"""
function megan_factors_init!(mf::MEGANFactors{FT}, n_classes::Int=20) where {FT}
    mf.n_classes = n_classes
    mf.Agro  = ones(n_classes)
    mf.Amat  = ones(n_classes)
    mf.Anew  = ones(n_classes)
    mf.Aold  = ones(n_classes)
    mf.betaT = fill(0.10, n_classes)
    mf.ct1   = fill(95.0, n_classes)
    mf.ct2   = fill(230.0, n_classes)
    mf.Ceo   = fill(2.0, n_classes)
    mf.LDF   = fill(0.6, n_classes)
    return nothing
end

# ---------------------------------------------------------------------------
# get_map_EF — Mapped emission factor for isoprene
# ---------------------------------------------------------------------------

"""
    get_map_EF(ivt, g, efisop_grc; kwargs...)

Get mapped emission factor for isoprene using gridded values for 6 PFT
groups following Guenther et al. (2006).
Returns emission factor in [ug m-2 h-1].

Ported from `get_map_EF` in `VOCEmissionMod.F90`.
"""
function get_map_EF(ivt::Int, g::Int, efisop_grc::AbstractMatrix{<:Real};
                    ndllf_evr_tmp_tree::Int = 1,
                    ndllf_evr_brl_tree::Int = 2,
                    ndllf_dcd_brl_tree::Int = 3,
                    nbrdlf_evr_trp_tree::Int = 4,
                    nbrdlf_evr_tmp_tree::Int = 5,
                    nbrdlf_dcd_trp_tree::Int = 6,
                    nbrdlf_dcd_tmp_tree::Int = 7,
                    nbrdlf_dcd_brl_tree::Int = 8,
                    nbrdlf_evr_shrub::Int = 9,
                    nbrdlf_dcd_brl_shrub::Int = 11,
                    nc3_arctic_grass::Int = 12,
                    nc4_grass::Int = 14,
                    nc3crop::Int = 15)
    T = eltype(efisop_grc)
    ef = zero(T)

    if ivt == ndllf_evr_tmp_tree || ivt == ndllf_evr_brl_tree
        # fineleaf evergreen
        ef = efisop_grc[2, g]
    elseif ivt == ndllf_dcd_brl_tree
        # fineleaf deciduous
        ef = efisop_grc[3, g]
    elseif ivt >= nbrdlf_evr_trp_tree && ivt <= nbrdlf_dcd_brl_tree
        # broadleaf trees
        ef = efisop_grc[1, g]
    elseif ivt >= nbrdlf_evr_shrub && ivt <= nbrdlf_dcd_brl_shrub
        # shrubs
        ef = efisop_grc[4, g]
    elseif ivt >= nc3_arctic_grass && ivt <= nc4_grass
        # grass
        ef = efisop_grc[5, g]
    elseif ivt >= nc3crop
        # crops
        ef = efisop_grc[6, g]
    end

    return ef
end

# ---------------------------------------------------------------------------
# get_gamma_P — Activity factor for PPFD
# ---------------------------------------------------------------------------

"""
    get_gamma_P(par_sun, par24_sun, par240_sun,
                par_sha, par24_sha, par240_sha,
                fsun, fsun240, forc_solad240, forc_solai240, LDF_in)

Activity factor for PPFD (Guenther et al., 2006).
Returns (gamma_p, cp, alpha).

Ported from `get_gamma_P` in `VOCEmissionMod.F90`.
"""
function get_gamma_P(par_sun::Real, par24_sun::Real, par240_sun::Real,
                     par_sha::Real, par24_sha::Real, par240_sha::Real,
                     fsun::Real, fsun240::Real,
                     forc_solad240::Real, forc_solai240::Real,
                     LDF_in::Real)
    T = typeof(par_sun)
    # Empirical coefficients
    ca1 = T(0.004)
    ca2 = T(0.0005)
    ca3 = T(0.0468)
    par0_sun = T(200.0)    # std conditions for past 24 hrs [umol/m2/s]
    par0_shade = T(50.0)   # std conditions for past 24 hrs [umol/m2/s]
    alpha_fix = T(0.001)   # empirical coefficient
    cp_fix = T(1.21)       # empirical coefficient

    cp = zero(T)
    alpha = zero(T)
    gamma_p_LDF = zero(T)

    if fsun240 > zero(T) && fsun240 < one(T) && forc_solad240 > zero(T) && forc_solai240 > zero(T)
        # With alpha and cp calculated based on eq 6 and 7:
        # SUN:
        alpha = ca1 - ca2 * log(par240_sun)
        cp = ca3 * exp(ca2 * (par24_sun - par0_sun)) * par240_sun^T(0.6)
        gamma_p_LDF = fsun * (cp * alpha * par_sun *
                      (one(T) + alpha * alpha * par_sun * par_sun)^T(-0.5))
        # SHADE:
        alpha = ca1 - ca2 * log(par240_sha)
        cp = ca3 * exp(ca2 * (par_sha - par0_shade)) * par240_sha^T(0.6)
        gamma_p_LDF = gamma_p_LDF + (one(T) - fsun) *
                      (cp * alpha * par_sha * (one(T) + alpha * alpha * par_sha * par_sha)^T(-0.5))
    else
        # With fixed alpha and cp (from MEGAN User's Guide):
        alpha = alpha_fix
        cp = cp_fix
        # SUN: direct + diffuse
        gamma_p_LDF = fsun * (cp * alpha * par_sun *
                      (one(T) + alpha * alpha * par_sun * par_sun)^T(-0.5))
        # SHADE: diffuse
        gamma_p_LDF = gamma_p_LDF + (one(T) - fsun) *
                      (cp * alpha * par_sha * (one(T) + alpha * alpha * par_sha * par_sha)^T(-0.5))
    end

    # Total activity factor accounting for light-dependent fraction
    gamma_p = (one(T) - LDF_in) + LDF_in * gamma_p_LDF

    return (gamma_p, cp, alpha)
end

# ---------------------------------------------------------------------------
# get_gamma_L — Activity factor for LAI
# ---------------------------------------------------------------------------

"""
    get_gamma_L(fsun240, elai)

Activity factor for LAI (Guenther et al., 2006, eq 3).
Ported from `get_gamma_L` in `VOCEmissionMod.F90`.
"""
function get_gamma_L(fsun240::Real, elai::Real)
    T = typeof(elai)
    cce  = T(0.30)   # factor to set emissions to unity @ std (accumulated vars available)
    cce1 = T(0.24)   # same but for non-accumulated vars

    if fsun240 > zero(T) && fsun240 < T(1.0e30)
        return cce * elai
    else
        return cce1 * elai
    end
end

# ---------------------------------------------------------------------------
# get_gamma_SM — Activity factor for soil moisture
# ---------------------------------------------------------------------------

"""
    get_gamma_SM(btran)

Activity factor for soil moisture of isoprene (Wang et al., 2022, JAMES).
Based on eq. (11) in the paper.
Ported from `get_gamma_SM` in `VOCEmissionMod.F90`.
"""
function get_gamma_SM(btran::Real)
    T = typeof(btran)
    a1 = T(-7.4463)
    b1 = T(3.2552)
    btran_threshold = T(0.2)

    if btran >= one(T)
        return one(T)
    else
        return one(T) / (one(T) + b1 * exp(a1 * (btran - btran_threshold)))
    end
end

# ---------------------------------------------------------------------------
# get_gamma_T — Activity factor for temperature
# ---------------------------------------------------------------------------

"""
    get_gamma_T(t_veg240, t_veg24, t_veg, ct1_in, ct2_in,
                betaT_in, LDF_in, Ceo_in, ivt;
                tfrz, nbrdlf_dcd_brl_shrub, nc3_arctic_grass)

Activity factor for temperature.
Returns (gamma_t, Eopt, topt).

Includes Wang et al. (2024, GRL) and Wang et al. (2024, Nature Communications)
updates for boreal broadleaf deciduous shrub and Arctic C3 grass.

Ported from `get_gamma_T` in `VOCEmissionMod.F90`.
"""
function get_gamma_T(t_veg240::Real, t_veg24::Real, t_veg::Real,
                     ct1_in::Real, ct2_in::Real,
                     betaT_in::Real, LDF_in::Real, Ceo_in::Real,
                     ivt::Int;
                     tfrz::Real = TFRZ,
                     nbrdlf_dcd_brl_shrub::Int = 11,
                     nc3_arctic_grass::Int = 12)
    T = typeof(t_veg)
    tfrz = T(tfrz)   # default TFRZ is a Float64 global; keep arithmetic in T (Metal)
    # Empirical coefficients
    co1 = T(313.0)
    co2 = T(0.6)
    co4 = T(0.05)
    tstd0 = T(297.0)          # std temperature [K]
    topt_fix = T(317.0)       # std temperature [K]
    Eopt_fix = T(2.26)        # empirical coefficient
    ct3 = T(0.00831)          # empirical coefficient
    tstd = T(303.15)          # std temperature [K]
    bet = T(0.09)             # beta empirical coefficient [K-1]

    # Boreal-specific parameters
    std_act_energy_isopr = T(95.0)
    empirical_param_1 = T(9.49)
    empirical_param_2 = T(0.53)
    empirical_param_3 = T(0.12)
    empirical_param_4 = T(7.9)
    empirical_param_5 = T(0.217)
    bet_arc_c3_max = T(300.0)

    Eopt = zero(T)
    topt = zero(T)

    # Light dependent fraction (Guenther et al., 2006)
    if t_veg240 > zero(T) && t_veg240 < T(1.0e30)
        # topt and Eopt from eq 8 and 9:
        topt = co1 + co2 * (t_veg240 - tstd0)

        if ivt == nbrdlf_dcd_brl_shrub
            # boreal-deciduous-shrub (BEAR-oNS campaign willows)
            Eopt = empirical_param_4 * exp(empirical_param_5 * (t_veg24 - tfrz - T(24.0)))
        elseif ivt == nc3_arctic_grass
            # boreal-grass
            Eopt = exp(empirical_param_3 * (t_veg240 - tfrz - T(15.0)))
        else
            Eopt = Ceo_in * exp(co4 * (t_veg24 - tstd0)) * exp(co4 * (t_veg240 - tstd0))
        end
    else
        topt = topt_fix
        Eopt = Eopt_fix
    end

    x = ((one(T) / topt) - (one(T) / t_veg)) / ct3

    # For the boreal grass (BEAR-oNS campaign)
    if ivt == nc3_arctic_grass
        bet_arc_c3 = std_act_energy_isopr + empirical_param_1 *
                     exp(empirical_param_2 * (tfrz + T(15.0) - t_veg240))
        bet_arc_c3 = min(bet_arc_c3, bet_arc_c3_max)
        gamma_t_LDF = Eopt * exp(bet_arc_c3 * ((one(T) / (tfrz + T(30.0)) - one(T) / t_veg) / ct3))
    else
        gamma_t_LDF = Eopt * (ct2_in * exp(ct1_in * x) /
                      (ct2_in - ct1_in * (one(T) - exp(ct2_in * x))))
    end

    # Light independent fraction (of exp(beta T) form)
    gamma_t_LIF = exp(betaT_in * (t_veg - tstd))

    # Total activity factor
    gamma_t = (one(T) - LDF_in) * gamma_t_LIF + LDF_in * gamma_t_LDF

    return (gamma_t, Eopt, topt)
end

# ---------------------------------------------------------------------------
# get_gamma_A — Activity factor for leaf age
# ---------------------------------------------------------------------------

"""
    get_gamma_A(ivt, elai240, elai, nclass, mf;
                ndllf_dcd_brl_tree, nbrdlf_dcd_trp_tree)

Activity factor for leaf age (Guenther et al., 2006).
Evergreens get gamma_a=1.0; deciduous PFTs use new/growing/mature/old fractions.
Ported from `get_gamma_A` in `VOCEmissionMod.F90`.
"""
function get_gamma_A(ivt::Int, elai240::Real, elai::Real,
                     nclass::Int, mf;
                     ndllf_dcd_brl_tree::Int = 3,
                     nbrdlf_dcd_trp_tree::Int = 6)
    T = typeof(elai)
    # non-evergreen test
    if ivt == ndllf_dcd_brl_tree || ivt >= nbrdlf_dcd_trp_tree
        if elai240 > zero(T) && elai240 < T(1.0e30)
            elai_prev = T(2.0) * elai240 - elai  # accumulated average over last 10 days

            if elai_prev == elai
                fnew = zero(T)
                fgro = zero(T)
                fmat = one(T)
                fold = zero(T)
            elseif elai_prev > elai
                fnew = zero(T)
                fgro = zero(T)
                fmat = one(T) - (elai_prev - elai) / elai_prev
                fold = (elai_prev - elai) / elai_prev
            else  # elai_prev < elai
                fnew = one(T) - (elai_prev / elai)
                fgro = zero(T)
                fmat = elai_prev / elai
                fold = zero(T)
            end

            return fnew * mf.Anew[nclass] + fgro * mf.Agro[nclass] +
                   fmat * mf.Amat[nclass] + fold * mf.Aold[nclass]
        else
            return one(T)
        end
    else
        return one(T)
    end
end

# ---------------------------------------------------------------------------
# get_gamma_C — Activity factor for CO2 (isoprene only)
# ---------------------------------------------------------------------------

"""
    get_gamma_C(cisun, cisha, forc_pbot, fsun, co2_ppmv)

Activity factor for instantaneous CO2 changes (Heald et al., 2009).
Includes both short-term (intercellular CO2) and long-term (ambient CO2)
exposure effects.
Ported from `get_gamma_C` in `VOCEmissionMod.F90`.
"""
function get_gamma_C(cisun::Real, cisha::Real,
                     forc_pbot::Real, fsun::Real, co2_ppmv::Real)
    T = typeof(forc_pbot)
    # Long-term exposure coefficients
    Ismax_ca   = T(1.344)
    h_ca       = T(1.4614)
    Cstar_ca   = T(585.0)
    CiCa_ratio = T(0.7)

    # LONG-TERM EXPOSURE (based on ambient CO2, Ca)
    gamma_ca = Ismax_ca - (Ismax_ca * (CiCa_ratio * co2_ppmv)^h_ca) /
               (Cstar_ca^h_ca + (CiCa_ratio * co2_ppmv)^h_ca)

    # SHORT-TERM EXPOSURE (based on intercellular CO2, Ci)
    # Determine coefficients based on long-term CO2 growth environment
    if co2_ppmv < T(400.0)
        Ismax = T(1.072)
        h     = T(1.70)
        Cstar = T(1218.0)
    elseif co2_ppmv > T(400.0) && co2_ppmv < T(600.0)
        fint  = (co2_ppmv - T(400.0)) / T(200.0)
        Ismax = fint * T(1.036) + (one(T) - fint) * T(1.072)
        h     = fint * T(2.0125) + (one(T) - fint) * T(1.70)
        Cstar = fint * T(1150.0) + (one(T) - fint) * T(1218.0)
    elseif co2_ppmv > T(600.0) && co2_ppmv < T(800.0)
        fint  = (co2_ppmv - T(600.0)) / T(200.0)
        Ismax = fint * T(1.046) + (one(T) - fint) * T(1.036)
        h     = fint * T(1.5380) + (one(T) - fint) * T(2.0125)
        Cstar = fint * T(2025.0) + (one(T) - fint) * T(1150.0)
    else  # co2_ppmv >= 800.0
        Ismax = T(1.014)
        h     = T(2.861)
        Cstar = T(1525.0)
    end

    # Intercellular CO2 concentrations
    # Fortran uses (x .eq. x) as NaN check; Julia: !isnan(x)
    if !isnan(cisun) && !isnan(cisha) && forc_pbot > zero(T) && fsun > zero(T)
        ci = (fsun * cisun + (one(T) - fsun) * cisha) / forc_pbot * T(1.0e6)
        gamma_ci = Ismax - (Ismax * ci^h) / (Cstar^h + ci^h)
    elseif cisun > zero(T) && cisun < T(1.0e30) && forc_pbot > zero(T) && fsun == one(T)
        ci = cisun / forc_pbot * T(1.0e6)
        gamma_ci = Ismax - (Ismax * ci^h) / (Cstar^h + ci^h)
    elseif cisha > zero(T) && cisha < T(1.0e30) && forc_pbot > zero(T) && fsun == zero(T)
        ci = cisha / forc_pbot * T(1.0e6)
        gamma_ci = Ismax - (Ismax * ci^h) / (Cstar^h + ci^h)
    else
        gamma_ci = one(T)
    end

    return gamma_ci * gamma_ca
end

# ---------------------------------------------------------------------------
# voc_emission! — Main VOC emission driver
# ---------------------------------------------------------------------------

"""
    voc_emission!(voc, meg_compounds, mech_comps, mf,
                  patch, bounds_p, mask_soilp,
                  forc_solad_col, forc_solai_grc, forc_pbot_col, forc_pco2_grc,
                  fsd24_patch, fsd240_patch, fsi24_patch, fsi240_patch,
                  canopystate, photosyn, temperature, energyflux;
                  kwargs...)

Main VOC emission driver. Calculates MEGAN v2.1 VOC emissions for all
active soil patches.

Ported from `VOCEmission` subroutine in `VOCEmissionMod.F90`.
"""
# ==========================================================================
# voc_emission! device-view bundles (group the ~50 kernel arrays under Metal's
# ~31-arg limit). Immutable + Adapt-able; aliased to Fortran-named locals at the
# kernel top so the per-patch body stays verbatim.
# ==========================================================================
Base.@kwdef struct _VocOut{V,M}   # arrays the kernel writes (+ efisop_grc it reads)
    vocflx_patch::M; vocflx_tot::V; meg_flux_out::M; efisop_grc::M
    Eopt_out::V; topt_out::V; alpha_out::V; cp_out::V
    paru_out::V; par24u_out::V; par240u_out::V
    para_out::V; par24a_out::V; par240a_out::V
    gamma_out::V; gammaL_out::V; gammaT_out::V; gammaP_out::V
    gammaA_out::V; gammaS_out::V; gammaC_out::V
end
Adapt.@adapt_structure _VocOut

Base.@kwdef struct _VocMeg{V,VI,VB,M}   # flattened mega-compound metadata
    is_isoprene::VB; class_number::VI; coeff::V; molec_weight::V; emis_factors::M
end
Adapt.@adapt_structure _VocMeg

Base.@kwdef struct _VocMech{VI,MI}   # flattened mechanism-compound mapping
    n_megan::VI; megan_idx::MI
end
Adapt.@adapt_structure _VocMech

# MEGANFactors is a *mutable* struct (not isbits) so it cannot be a kernel arg;
# its per-class vectors are repacked into this immutable bundle (get_gamma_A/T
# read them on the device). Same field names as MEGANFactors so get_gamma_A's
# `mf.Anew[nclass]` etc. work unchanged.
Base.@kwdef struct _VocMF{V}
    Anew::V; Agro::V; Amat::V; Aold::V
    betaT::V; ct1::V; ct2::V; Ceo::V; LDF::V
end
Adapt.@adapt_structure _VocMF

Base.@kwdef struct _VocForce{V,M}   # forcing / canopy / temperature inputs
    forc_solad_col::M; forc_solai_grc::M; forc_pbot_col::V; forc_pco2_grc::V
    fsd24::V; fsd240::V; fsi24::V; fsi240::V
    fsun::V; fsun24::V; fsun240::V; elai::V; elai240::V
    cisun_z::M; cisha_z::M; t_veg::V; t_veg24::V; t_veg240::V; btran::V
end
Adapt.@adapt_structure _VocForce

# One thread per patch. Per-thread `vmeg[p, imeg]` holds the per-mega-compound
# moles flux (n_megcomps is a runtime value, so it lives in a device-resident
# scratch matrix rather than a thread-local array). The get_* helpers are called
# positionally → they use their DEFAULT PFT-index kwargs, which equal this
# function's defaults (the only configuration any caller uses).
@kernel function _voc_emission_kernel!(out, meg, mech, mf, force, vmeg,
        @Const(mask_soilp), @Const(p_gridcell), @Const(p_column), @Const(p_itype),
        use_mapped::Bool, megfac, tfrzT, n_megcomps::Int, n_mechcomps::Int, pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax
        T = eltype(out.vocflx_tot)
        # init outputs for this patch (unconditional, mirrors the host init loops)
        out.vocflx_tot[p] = zero(T)
        for imech in 1:n_mechcomps
            out.vocflx_patch[p, imech] = zero(T)
        end
        for imeg in 1:n_megcomps
            out.meg_flux_out[imeg, p] = zero(T)
            vmeg[p, imeg] = zero(T)
        end

        if mask_soilp[p] && p_itype[p] > 0
            g = p_gridcell[p]; c = p_column[p]; ivt = p_itype[p]
            # PAR [umol/m2/s] = W/m2 * 4.6
            par_sun    = (force.forc_solad_col[c,1] + force.fsun[p]    * force.forc_solai_grc[g,1]) * T(4.6)
            par24_sun  = (force.fsd24[p]            + force.fsun24[p]  * force.fsi24[p])             * T(4.6)
            par240_sun = (force.fsd240[p]           + force.fsun240[p] * force.fsi240[p])           * T(4.6)
            par_sha    = ((one(T) - force.fsun[p])    * force.forc_solai_grc[g,1]) * T(4.6)
            par24_sha  = ((one(T) - force.fsun24[p])  * force.fsi24[p])            * T(4.6)
            par240_sha = ((one(T) - force.fsun240[p]) * force.fsi240[p])           * T(4.6)

            gamma_l  = get_gamma_L(force.fsun240[p], force.elai[p])
            gamma_sm = get_gamma_SM(force.btran[p])

            for imeg in 1:n_megcomps
                if meg.is_isoprene[imeg] && use_mapped
                    epsilon = get_map_EF(ivt, g, out.efisop_grc)
                else
                    epsilon = meg.emis_factors[imeg, ivt]
                end
                class_num = meg.class_number[imeg]

                (gamma_p, cp, alpha) = get_gamma_P(par_sun, par24_sun, par240_sun,
                    par_sha, par24_sha, par240_sha, force.fsun[p], force.fsun240[p],
                    force.fsd240[p], force.fsi240[p], mf.LDF[class_num])
                (gamma_t, Eopt, topt) = get_gamma_T(force.t_veg240[p], force.t_veg24[p],
                    force.t_veg[p], mf.ct1[class_num], mf.ct2[class_num], mf.betaT[class_num],
                    mf.LDF[class_num], mf.Ceo[class_num], ivt; tfrz = tfrzT)
                gamma_a = get_gamma_A(ivt, force.elai240[p], force.elai[p], class_num, mf)
                if meg.is_isoprene[imeg]
                    co2_ppmv = T(1.0e6) * force.forc_pco2_grc[g] / force.forc_pbot_col[c]
                    gamma_c = get_gamma_C(force.cisun_z[p,1], force.cisha_z[p,1],
                                          force.forc_pbot_col[c], force.fsun[p], co2_ppmv)
                else
                    gamma_c = one(T)
                end

                gamma = gamma_l * gamma_sm * gamma_a * gamma_p * gamma_t * gamma_c

                if gamma >= zero(T) && gamma < T(100.0)
                    vmeg[p, imeg] = meg.coeff[imeg] * epsilon * gamma * megfac / meg.molec_weight[imeg]
                    out.meg_flux_out[imeg, p] = epsilon * gamma * megfac * T(1.0e-3)
                    if imeg == 1
                        out.gamma_out[p]  = gamma;  out.gammaP_out[p] = gamma_p; out.gammaT_out[p] = gamma_t
                        out.gammaA_out[p] = gamma_a; out.gammaS_out[p] = gamma_sm; out.gammaL_out[p] = gamma_l
                        out.gammaC_out[p] = gamma_c
                        out.paru_out[p]   = par_sun;  out.par24u_out[p] = par24_sun; out.par240u_out[p] = par240_sun
                        out.para_out[p]   = par_sha;  out.par24a_out[p] = par24_sha; out.par240a_out[p] = par240_sha
                        out.alpha_out[p]  = alpha;    out.cp_out[p]     = cp
                        out.topt_out[p]   = topt;     out.Eopt_out[p]   = Eopt
                    end
                end
            end

            # Sum mega-compound fluxes into mechanism compounds (local accum,
            # one write each — byte-identical to the host += and --check-bounds safe)
            tot = zero(T)
            for imech in 1:n_mechcomps
                acc = zero(T)
                for ii in 1:mech.n_megan[imech]
                    acc += vmeg[p, mech.megan_idx[imech, ii]]
                end
                out.vocflx_patch[p, imech] = acc
                tot += acc
            end
            out.vocflx_tot[p] = tot
        end
    end
end

function voc_emission!(
    # VOC data
    voc::VOCEmisData,
    # MEGAN compound descriptors
    meg_compounds::Vector{<:MEGANCompound},
    mech_comps::Vector{MEGANMechComp},
    mf::MEGANFactors,
    # Patch structure
    patch::PatchData,
    bounds_p::UnitRange{Int},
    mask_soilp::AbstractVector{Bool},
    # Atmospheric forcing (bare arrays, as in surface_radiation pattern)
    forc_solad_col::AbstractMatrix{<:Real},   # direct beam radiation (ncol, numrad) [W/m2]
    forc_solai_grc::AbstractMatrix{<:Real},   # diffuse radiation (ngrc, numrad) [W/m2]
    forc_pbot_col::AbstractVector{<:Real},    # atmospheric pressure (ncol) [Pa]
    forc_pco2_grc::AbstractVector{<:Real},    # partial pressure CO2 (ngrc) [Pa]
    # Time-averaged PAR forcing (patch level)
    fsd24_patch::AbstractVector{<:Real},      # direct beam 24hr avg [W/m2]
    fsd240_patch::AbstractVector{<:Real},     # direct beam 240hr avg [W/m2]
    fsi24_patch::AbstractVector{<:Real},      # diffuse 24hr avg [W/m2]
    fsi240_patch::AbstractVector{<:Real},     # diffuse 240hr avg [W/m2]
    # Canopy state
    fsun_patch::AbstractVector{<:Real},
    fsun24_patch::AbstractVector{<:Real},
    fsun240_patch::AbstractVector{<:Real},
    elai_patch::AbstractVector{<:Real},
    elai240_patch::AbstractVector{<:Real},
    # Photosynthesis
    cisun_z_patch::AbstractMatrix{<:Real},    # sunlit intracellular CO2 (np, nlevcan) [Pa]
    cisha_z_patch::AbstractMatrix{<:Real},    # shaded intracellular CO2 (np, nlevcan) [Pa]
    # Temperature
    t_veg_patch::AbstractVector{<:Real},
    t_veg24_patch::AbstractVector{<:Real},
    t_veg240_patch::AbstractVector{<:Real},
    # Energy flux
    btran_patch::AbstractVector{<:Real};
    # Keyword arguments for PFT indices
    use_mapped_emisfctrs::Bool = true,
    ndllf_evr_tmp_tree::Int = 1,
    ndllf_evr_brl_tree::Int = 2,
    ndllf_dcd_brl_tree::Int = 3,
    nbrdlf_evr_trp_tree::Int = 4,
    nbrdlf_evr_tmp_tree::Int = 5,
    nbrdlf_dcd_trp_tree::Int = 6,
    nbrdlf_dcd_tmp_tree::Int = 7,
    nbrdlf_dcd_brl_tree::Int = 8,
    nbrdlf_evr_shrub::Int = 9,
    nbrdlf_dcd_brl_shrub::Int = 11,
    nc3_arctic_grass::Int = 12,
    nc4_grass::Int = 14,
    nc3crop::Int = 15,
    noveg::Int = 0
)
    n_megcomps = length(meg_compounds)
    n_mechcomps = length(mech_comps)

    if n_mechcomps < 1
        return nothing
    end

    FT = eltype(voc.vocflx_tot_patch)
    # factor to convert MEGAN units [ug/m2/hr] to [g/m2/sec]
    megemis_units_factor = FT(1.0 / 3600.0 / 1.0e6)

    # --- Flatten the Vector{MEGANCompound}/Vector{MEGANMechComp} metadata (which
    #     carry String + nested-vector fields, not device-movable) into plain
    #     numeric arrays, then place them on the backend of the voc outputs.
    maxpft = 0
    for mc in meg_compounds; maxpft = max(maxpft, length(mc.emis_factors)); end
    maxnmg = 0
    for mc in mech_comps; maxnmg = max(maxnmg, mc.n_megan_comps); end
    is_iso_h = Bool[mc.name == "isoprene" for mc in meg_compounds]
    class_h  = Int[mc.class_number for mc in meg_compounds]
    coeff_h  = FT[mc.coeff for mc in meg_compounds]
    molw_h   = FT[mc.molec_weight for mc in meg_compounds]
    ef_h     = zeros(FT, n_megcomps, max(maxpft, 1))
    for (i, mc) in enumerate(meg_compounds), k in 1:length(mc.emis_factors)
        ef_h[i, k] = mc.emis_factors[k]
    end
    nmg_h  = Int[mc.n_megan_comps for mc in mech_comps]
    midx_h = zeros(Int, n_mechcomps, max(maxnmg, 1))
    for (i, mc) in enumerate(mech_comps), ii in 1:mc.n_megan_comps
        midx_h[i, ii] = mc.megan_indices[ii]
    end

    # place a host array onto the backend of the voc output arrays
    _md(h) = copyto!(similar(voc.vocflx_tot_patch, eltype(h), size(h)...), h)

    out = _VocOut(;
        vocflx_patch = voc.vocflx_patch, vocflx_tot = voc.vocflx_tot_patch,
        meg_flux_out = voc.meg_flux_out, efisop_grc = voc.efisop_grc,
        Eopt_out = voc.Eopt_out_patch, topt_out = voc.topt_out_patch,
        alpha_out = voc.alpha_out_patch, cp_out = voc.cp_out_patch,
        paru_out = voc.paru_out_patch, par24u_out = voc.par24u_out_patch, par240u_out = voc.par240u_out_patch,
        para_out = voc.para_out_patch, par24a_out = voc.par24a_out_patch, par240a_out = voc.par240a_out_patch,
        gamma_out = voc.gamma_out_patch, gammaL_out = voc.gammaL_out_patch, gammaT_out = voc.gammaT_out_patch,
        gammaP_out = voc.gammaP_out_patch, gammaA_out = voc.gammaA_out_patch, gammaS_out = voc.gammaS_out_patch,
        gammaC_out = voc.gammaC_out_patch)
    megb = _VocMeg(; is_isoprene = _md(is_iso_h), class_number = _md(class_h),
        coeff = _md(coeff_h), molec_weight = _md(molw_h), emis_factors = _md(ef_h))
    mechb = _VocMech(; n_megan = _md(nmg_h), megan_idx = _md(midx_h))
    mfb = _VocMF(; Anew = mf.Anew, Agro = mf.Agro, Amat = mf.Amat, Aold = mf.Aold,
        betaT = mf.betaT, ct1 = mf.ct1, ct2 = mf.ct2, Ceo = mf.Ceo, LDF = mf.LDF)
    force = _VocForce(;
        forc_solad_col = forc_solad_col, forc_solai_grc = forc_solai_grc,
        forc_pbot_col = forc_pbot_col, forc_pco2_grc = forc_pco2_grc,
        fsd24 = fsd24_patch, fsd240 = fsd240_patch, fsi24 = fsi24_patch, fsi240 = fsi240_patch,
        fsun = fsun_patch, fsun24 = fsun24_patch, fsun240 = fsun240_patch,
        elai = elai_patch, elai240 = elai240_patch,
        cisun_z = cisun_z_patch, cisha_z = cisha_z_patch,
        t_veg = t_veg_patch, t_veg24 = t_veg24_patch, t_veg240 = t_veg240_patch,
        btran = btran_patch)

    vmeg = fill!(similar(voc.vocflx_tot_patch, FT, length(voc.vocflx_tot_patch), n_megcomps), zero(FT))

    tfrzT = FT(TFRZ)   # convert the Float64 TFRZ global once on the host (Metal: no Float64)
    backend = _kernel_backend(voc.vocflx_tot_patch)
    _voc_emission_kernel!(backend)(out, megb, mechb, mfb, force, vmeg,
        mask_soilp, patch.gridcell, patch.column, patch.itype,
        use_mapped_emisfctrs, megemis_units_factor, tfrzT, n_megcomps, n_mechcomps,
        first(bounds_p), last(bounds_p); ndrange = length(voc.vocflx_tot_patch))
    KA.synchronize(backend)

    return nothing
end
