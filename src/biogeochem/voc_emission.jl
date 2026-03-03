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
Base.@kwdef mutable struct MEGANFactors
    n_classes::Int = 20
    # Leaf age fractions (indexed by class number, 1:n_classes)
    Agro ::Vector{Float64} = Float64[]  # growing leaves factor
    Amat ::Vector{Float64} = Float64[]  # mature leaves factor
    Anew ::Vector{Float64} = Float64[]  # new leaves factor
    Aold ::Vector{Float64} = Float64[]  # old leaves factor
    # Temperature parameters
    betaT::Vector{Float64} = Float64[]  # beta for T-dependent LIF
    ct1  ::Vector{Float64} = Float64[]  # ct1 coefficient for T activity factor
    ct2  ::Vector{Float64} = Float64[]  # ct2 coefficient for T activity factor
    Ceo  ::Vector{Float64} = Float64[]  # Ceo coefficient for Eopt
    # Light dependence fraction
    LDF  ::Vector{Float64} = Float64[]  # light-dependent fraction (0-1)
end

"""
    MEGANCompound

Descriptor for a single MEGAN mega-compound.
Ported from `shr_megan_megcomp_t` in `shr_megan_mod`.
"""
Base.@kwdef mutable struct MEGANCompound
    name          ::String          = ""        # compound name (e.g. "isoprene")
    index         ::Int             = 0         # index in mega-compound list
    class_number  ::Int             = 0         # MEGAN class number (1-20)
    molec_weight  ::Float64         = 0.0       # molecular weight [g/mol]
    coeff         ::Float64         = 1.0       # coefficient for mapping to mechanism compound
    emis_factors  ::Vector{Float64} = Float64[] # emission factors per PFT [ug m-2 h-1]
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
Base.@kwdef mutable struct VOCEmisData
    # Diagnostic coefficients (patch level)
    Eopt_out_patch    ::Vector{Float64} = Float64[]  # Eopt coefficient
    topt_out_patch    ::Vector{Float64} = Float64[]  # topt coefficient
    alpha_out_patch   ::Vector{Float64} = Float64[]  # alpha coefficient
    cp_out_patch      ::Vector{Float64} = Float64[]  # cp coefficient

    # PAR diagnostics (patch level)
    paru_out_patch    ::Vector{Float64} = Float64[]  # sunlit PAR [umol/m2/s]
    par24u_out_patch  ::Vector{Float64} = Float64[]  # sunlit PAR 24hr avg
    par240u_out_patch ::Vector{Float64} = Float64[]  # sunlit PAR 240hr avg
    para_out_patch    ::Vector{Float64} = Float64[]  # shade PAR [umol/m2/s]
    par24a_out_patch  ::Vector{Float64} = Float64[]  # shade PAR 24hr avg
    par240a_out_patch ::Vector{Float64} = Float64[]  # shade PAR 240hr avg

    # Gamma diagnostics (patch level)
    gamma_out_patch   ::Vector{Float64} = Float64[]  # total gamma
    gammaL_out_patch  ::Vector{Float64} = Float64[]  # gamma for LAI
    gammaT_out_patch  ::Vector{Float64} = Float64[]  # gamma for temperature
    gammaP_out_patch  ::Vector{Float64} = Float64[]  # gamma for PPFD
    gammaA_out_patch  ::Vector{Float64} = Float64[]  # gamma for leaf age
    gammaS_out_patch  ::Vector{Float64} = Float64[]  # gamma for soil moisture
    gammaC_out_patch  ::Vector{Float64} = Float64[]  # gamma for CO2

    # Fluxes (patch level)
    vocflx_tot_patch  ::Vector{Float64} = Float64[]  # total VOC flux [moles/m2/sec]
    vocflx_patch      ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # per-mechanism-compound flux [moles/m2/sec]

    # Gridcell isoprene emission factors (6 × ngrc)
    efisop_grc        ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # [ug m-2 h-1]

    # Per-mega-compound flux output (n_megcomps × np) [kg/m2/sec] for history
    meg_flux_out      ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
end

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
function vocemis_clean!(voc::VOCEmisData)
    voc.Eopt_out_patch    = Float64[]
    voc.topt_out_patch    = Float64[]
    voc.alpha_out_patch   = Float64[]
    voc.cp_out_patch      = Float64[]
    voc.paru_out_patch    = Float64[]
    voc.par24u_out_patch  = Float64[]
    voc.par240u_out_patch = Float64[]
    voc.para_out_patch    = Float64[]
    voc.par24a_out_patch  = Float64[]
    voc.par240a_out_patch = Float64[]
    voc.gamma_out_patch   = Float64[]
    voc.gammaL_out_patch  = Float64[]
    voc.gammaT_out_patch  = Float64[]
    voc.gammaP_out_patch  = Float64[]
    voc.gammaA_out_patch  = Float64[]
    voc.gammaS_out_patch  = Float64[]
    voc.gammaC_out_patch  = Float64[]
    voc.vocflx_tot_patch  = Float64[]
    voc.vocflx_patch      = Matrix{Float64}(undef, 0, 0)
    voc.efisop_grc        = Matrix{Float64}(undef, 0, 0)
    voc.meg_flux_out      = Matrix{Float64}(undef, 0, 0)
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
function megan_factors_init!(mf::MEGANFactors, n_classes::Int=20)
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
function get_map_EF(ivt::Int, g::Int, efisop_grc::Matrix{Float64};
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
    ef = 0.0

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
function get_gamma_P(par_sun::Float64, par24_sun::Float64, par240_sun::Float64,
                     par_sha::Float64, par24_sha::Float64, par240_sha::Float64,
                     fsun::Float64, fsun240::Float64,
                     forc_solad240::Float64, forc_solai240::Float64,
                     LDF_in::Float64)
    # Empirical coefficients
    ca1 = 0.004
    ca2 = 0.0005
    ca3 = 0.0468
    par0_sun = 200.0    # std conditions for past 24 hrs [umol/m2/s]
    par0_shade = 50.0   # std conditions for past 24 hrs [umol/m2/s]
    alpha_fix = 0.001   # empirical coefficient
    cp_fix = 1.21       # empirical coefficient

    cp = 0.0
    alpha = 0.0
    gamma_p_LDF = 0.0

    if fsun240 > 0.0 && fsun240 < 1.0 && forc_solad240 > 0.0 && forc_solai240 > 0.0
        # With alpha and cp calculated based on eq 6 and 7:
        # SUN:
        alpha = ca1 - ca2 * log(par240_sun)
        cp = ca3 * exp(ca2 * (par24_sun - par0_sun)) * par240_sun^0.6
        gamma_p_LDF = fsun * (cp * alpha * par_sun *
                      (1.0 + alpha * alpha * par_sun * par_sun)^(-0.5))
        # SHADE:
        alpha = ca1 - ca2 * log(par240_sha)
        cp = ca3 * exp(ca2 * (par_sha - par0_shade)) * par240_sha^0.6
        gamma_p_LDF = gamma_p_LDF + (1.0 - fsun) *
                      (cp * alpha * par_sha * (1.0 + alpha * alpha * par_sha * par_sha)^(-0.5))
    else
        # With fixed alpha and cp (from MEGAN User's Guide):
        alpha = alpha_fix
        cp = cp_fix
        # SUN: direct + diffuse
        gamma_p_LDF = fsun * (cp * alpha * par_sun *
                      (1.0 + alpha * alpha * par_sun * par_sun)^(-0.5))
        # SHADE: diffuse
        gamma_p_LDF = gamma_p_LDF + (1.0 - fsun) *
                      (cp * alpha * par_sha * (1.0 + alpha * alpha * par_sha * par_sha)^(-0.5))
    end

    # Total activity factor accounting for light-dependent fraction
    gamma_p = (1.0 - LDF_in) + LDF_in * gamma_p_LDF

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
function get_gamma_L(fsun240::Float64, elai::Float64)
    cce  = 0.30   # factor to set emissions to unity @ std (accumulated vars available)
    cce1 = 0.24   # same but for non-accumulated vars

    if fsun240 > 0.0 && fsun240 < 1.0e30
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
function get_gamma_SM(btran::Float64)
    a1 = -7.4463
    b1 = 3.2552
    btran_threshold = 0.2

    if btran >= 1.0
        return 1.0
    else
        return 1.0 / (1.0 + b1 * exp(a1 * (btran - btran_threshold)))
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
function get_gamma_T(t_veg240::Float64, t_veg24::Float64, t_veg::Float64,
                     ct1_in::Float64, ct2_in::Float64,
                     betaT_in::Float64, LDF_in::Float64, Ceo_in::Float64,
                     ivt::Int;
                     tfrz::Float64 = TFRZ,
                     nbrdlf_dcd_brl_shrub::Int = 11,
                     nc3_arctic_grass::Int = 12)
    # Empirical coefficients
    co1 = 313.0
    co2 = 0.6
    co4 = 0.05
    tstd0 = 297.0          # std temperature [K]
    topt_fix = 317.0       # std temperature [K]
    Eopt_fix = 2.26        # empirical coefficient
    ct3 = 0.00831          # empirical coefficient
    tstd = 303.15          # std temperature [K]
    bet = 0.09             # beta empirical coefficient [K-1]

    # Boreal-specific parameters
    std_act_energy_isopr = 95.0
    empirical_param_1 = 9.49
    empirical_param_2 = 0.53
    empirical_param_3 = 0.12
    empirical_param_4 = 7.9
    empirical_param_5 = 0.217
    bet_arc_c3_max = 300.0

    Eopt = 0.0
    topt = 0.0

    # Light dependent fraction (Guenther et al., 2006)
    if t_veg240 > 0.0 && t_veg240 < 1.0e30
        # topt and Eopt from eq 8 and 9:
        topt = co1 + co2 * (t_veg240 - tstd0)

        if ivt == nbrdlf_dcd_brl_shrub
            # boreal-deciduous-shrub (BEAR-oNS campaign willows)
            Eopt = empirical_param_4 * exp(empirical_param_5 * (t_veg24 - tfrz - 24.0))
        elseif ivt == nc3_arctic_grass
            # boreal-grass
            Eopt = exp(empirical_param_3 * (t_veg240 - tfrz - 15.0))
        else
            Eopt = Ceo_in * exp(co4 * (t_veg24 - tstd0)) * exp(co4 * (t_veg240 - tstd0))
        end
    else
        topt = topt_fix
        Eopt = Eopt_fix
    end

    x = ((1.0 / topt) - (1.0 / t_veg)) / ct3

    # For the boreal grass (BEAR-oNS campaign)
    if ivt == nc3_arctic_grass
        bet_arc_c3 = std_act_energy_isopr + empirical_param_1 *
                     exp(empirical_param_2 * (tfrz + 15.0 - t_veg240))
        bet_arc_c3 = min(bet_arc_c3, bet_arc_c3_max)
        gamma_t_LDF = Eopt * exp(bet_arc_c3 * ((1.0 / (tfrz + 30.0) - 1.0 / t_veg) / ct3))
    else
        gamma_t_LDF = Eopt * (ct2_in * exp(ct1_in * x) /
                      (ct2_in - ct1_in * (1.0 - exp(ct2_in * x))))
    end

    # Light independent fraction (of exp(beta T) form)
    gamma_t_LIF = exp(betaT_in * (t_veg - tstd))

    # Total activity factor
    gamma_t = (1.0 - LDF_in) * gamma_t_LIF + LDF_in * gamma_t_LDF

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
function get_gamma_A(ivt::Int, elai240::Float64, elai::Float64,
                     nclass::Int, mf::MEGANFactors;
                     ndllf_dcd_brl_tree::Int = 3,
                     nbrdlf_dcd_trp_tree::Int = 6)
    # non-evergreen test
    if ivt == ndllf_dcd_brl_tree || ivt >= nbrdlf_dcd_trp_tree
        if elai240 > 0.0 && elai240 < 1.0e30
            elai_prev = 2.0 * elai240 - elai  # accumulated average over last 10 days

            if elai_prev == elai
                fnew = 0.0
                fgro = 0.0
                fmat = 1.0
                fold = 0.0
            elseif elai_prev > elai
                fnew = 0.0
                fgro = 0.0
                fmat = 1.0 - (elai_prev - elai) / elai_prev
                fold = (elai_prev - elai) / elai_prev
            else  # elai_prev < elai
                fnew = 1.0 - (elai_prev / elai)
                fgro = 0.0
                fmat = elai_prev / elai
                fold = 0.0
            end

            return fnew * mf.Anew[nclass] + fgro * mf.Agro[nclass] +
                   fmat * mf.Amat[nclass] + fold * mf.Aold[nclass]
        else
            return 1.0
        end
    else
        return 1.0
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
function get_gamma_C(cisun::Float64, cisha::Float64,
                     forc_pbot::Float64, fsun::Float64, co2_ppmv::Float64)
    # Long-term exposure coefficients
    Ismax_ca   = 1.344
    h_ca       = 1.4614
    Cstar_ca   = 585.0
    CiCa_ratio = 0.7

    # LONG-TERM EXPOSURE (based on ambient CO2, Ca)
    gamma_ca = Ismax_ca - (Ismax_ca * (CiCa_ratio * co2_ppmv)^h_ca) /
               (Cstar_ca^h_ca + (CiCa_ratio * co2_ppmv)^h_ca)

    # SHORT-TERM EXPOSURE (based on intercellular CO2, Ci)
    # Determine coefficients based on long-term CO2 growth environment
    if co2_ppmv < 400.0
        Ismax = 1.072
        h     = 1.70
        Cstar = 1218.0
    elseif co2_ppmv > 400.0 && co2_ppmv < 600.0
        fint  = (co2_ppmv - 400.0) / 200.0
        Ismax = fint * 1.036 + (1.0 - fint) * 1.072
        h     = fint * 2.0125 + (1.0 - fint) * 1.70
        Cstar = fint * 1150.0 + (1.0 - fint) * 1218.0
    elseif co2_ppmv > 600.0 && co2_ppmv < 800.0
        fint  = (co2_ppmv - 600.0) / 200.0
        Ismax = fint * 1.046 + (1.0 - fint) * 1.036
        h     = fint * 1.5380 + (1.0 - fint) * 2.0125
        Cstar = fint * 2025.0 + (1.0 - fint) * 1150.0
    else  # co2_ppmv >= 800.0
        Ismax = 1.014
        h     = 2.861
        Cstar = 1525.0
    end

    # Intercellular CO2 concentrations
    # Fortran uses (x .eq. x) as NaN check; Julia: !isnan(x)
    if !isnan(cisun) && !isnan(cisha) && forc_pbot > 0.0 && fsun > 0.0
        ci = (fsun * cisun + (1.0 - fsun) * cisha) / forc_pbot * 1.0e6
        gamma_ci = Ismax - (Ismax * ci^h) / (Cstar^h + ci^h)
    elseif cisun > 0.0 && cisun < 1.0e30 && forc_pbot > 0.0 && fsun == 1.0
        ci = cisun / forc_pbot * 1.0e6
        gamma_ci = Ismax - (Ismax * ci^h) / (Cstar^h + ci^h)
    elseif cisha > 0.0 && cisha < 1.0e30 && forc_pbot > 0.0 && fsun == 0.0
        ci = cisha / forc_pbot * 1.0e6
        gamma_ci = Ismax - (Ismax * ci^h) / (Cstar^h + ci^h)
    else
        gamma_ci = 1.0
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
function voc_emission!(
    # VOC data
    voc::VOCEmisData,
    # MEGAN compound descriptors
    meg_compounds::Vector{MEGANCompound},
    mech_comps::Vector{MEGANMechComp},
    mf::MEGANFactors,
    # Patch structure
    patch::PatchData,
    bounds_p::UnitRange{Int},
    mask_soilp::BitVector,
    # Atmospheric forcing (bare arrays, as in surface_radiation pattern)
    forc_solad_col::Matrix{Float64},   # direct beam radiation (ncol, numrad) [W/m2]
    forc_solai_grc::Matrix{Float64},   # diffuse radiation (ngrc, numrad) [W/m2]
    forc_pbot_col::Vector{Float64},    # atmospheric pressure (ncol) [Pa]
    forc_pco2_grc::Vector{Float64},    # partial pressure CO2 (ngrc) [Pa]
    # Time-averaged PAR forcing (patch level)
    fsd24_patch::Vector{Float64},      # direct beam 24hr avg [W/m2]
    fsd240_patch::Vector{Float64},     # direct beam 240hr avg [W/m2]
    fsi24_patch::Vector{Float64},      # diffuse 24hr avg [W/m2]
    fsi240_patch::Vector{Float64},     # diffuse 240hr avg [W/m2]
    # Canopy state
    fsun_patch::Vector{Float64},
    fsun24_patch::Vector{Float64},
    fsun240_patch::Vector{Float64},
    elai_patch::Vector{Float64},
    elai240_patch::Vector{Float64},
    # Photosynthesis
    cisun_z_patch::Matrix{Float64},    # sunlit intracellular CO2 (np, nlevcan) [Pa]
    cisha_z_patch::Matrix{Float64},    # shaded intracellular CO2 (np, nlevcan) [Pa]
    # Temperature
    t_veg_patch::Vector{Float64},
    t_veg24_patch::Vector{Float64},
    t_veg240_patch::Vector{Float64},
    # Energy flux
    btran_patch::Vector{Float64};
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

    # factor to convert MEGAN units [ug/m2/hr] to [g/m2/sec]
    megemis_units_factor = 1.0 / 3600.0 / 1.0e6

    # Initialize output fluxes
    for p in bounds_p
        for imech in 1:n_mechcomps
            voc.vocflx_patch[p, imech] = 0.0
        end
        voc.vocflx_tot_patch[p] = 0.0
    end
    for imeg in 1:n_megcomps
        for p in bounds_p
            voc.meg_flux_out[imeg, p] = 0.0
        end
    end

    # Temporary per-mega-compound flux
    vocflx_meg = zeros(n_megcomps)

    # Begin loop over patches
    for p in bounds_p
        mask_soilp[p] || continue

        g = patch.gridcell[p]
        c = patch.column[p]

        # initialize EF and per-compound fluxes
        epsilon = 0.0
        fill!(vocflx_meg, 0.0)

        # calculate VOC emissions for non-bare ground patches
        if patch.itype[p] > 0

            # Calculate PAR: multiply W/m2 by 4.6 to get umol/m2/s
            # SUN:
            par_sun    = (forc_solad_col[c, 1] + fsun_patch[p]    * forc_solai_grc[g, 1]) * 4.6
            par24_sun  = (fsd24_patch[p]        + fsun24_patch[p]  * fsi24_patch[p])       * 4.6
            par240_sun = (fsd240_patch[p]        + fsun240_patch[p] * fsi240_patch[p])      * 4.6

            # SHADE:
            par_sha    = ((1.0 - fsun_patch[p])    * forc_solai_grc[g, 1]) * 4.6
            par24_sha  = ((1.0 - fsun24_patch[p])  * fsi24_patch[p])       * 4.6
            par240_sha = ((1.0 - fsun240_patch[p]) * fsi240_patch[p])      * 4.6

            # Activity factor for LAI (all species)
            gamma_l = get_gamma_L(fsun240_patch[p], elai_patch[p])

            # Activity factor for soil moisture
            gamma_sm = get_gamma_SM(btran_patch[p])

            # Loop through MEGAN compounds
            for imeg in 1:n_megcomps
                meg_cmp = meg_compounds[imeg]

                # Set emission factor
                if meg_cmp.name == "isoprene" && use_mapped_emisfctrs
                    epsilon = get_map_EF(patch.itype[p], g, voc.efisop_grc;
                                         ndllf_evr_tmp_tree=ndllf_evr_tmp_tree,
                                         ndllf_evr_brl_tree=ndllf_evr_brl_tree,
                                         ndllf_dcd_brl_tree=ndllf_dcd_brl_tree,
                                         nbrdlf_evr_trp_tree=nbrdlf_evr_trp_tree,
                                         nbrdlf_evr_tmp_tree=nbrdlf_evr_tmp_tree,
                                         nbrdlf_dcd_trp_tree=nbrdlf_dcd_trp_tree,
                                         nbrdlf_dcd_tmp_tree=nbrdlf_dcd_tmp_tree,
                                         nbrdlf_dcd_brl_tree=nbrdlf_dcd_brl_tree,
                                         nbrdlf_evr_shrub=nbrdlf_evr_shrub,
                                         nbrdlf_dcd_brl_shrub=nbrdlf_dcd_brl_shrub,
                                         nc3_arctic_grass=nc3_arctic_grass,
                                         nc4_grass=nc4_grass,
                                         nc3crop=nc3crop)
                else
                    epsilon = meg_cmp.emis_factors[patch.itype[p]]
                end

                class_num = meg_cmp.class_number

                # Activity factor for PPFD
                (gamma_p, cp, alpha) = get_gamma_P(
                    par_sun, par24_sun, par240_sun,
                    par_sha, par24_sha, par240_sha,
                    fsun_patch[p], fsun240_patch[p],
                    fsd240_patch[p], fsi240_patch[p],
                    mf.LDF[class_num])

                # Activity factor for temperature
                (gamma_t, Eopt, topt) = get_gamma_T(
                    t_veg240_patch[p], t_veg24_patch[p], t_veg_patch[p],
                    mf.ct1[class_num], mf.ct2[class_num],
                    mf.betaT[class_num], mf.LDF[class_num], mf.Ceo[class_num],
                    patch.itype[p];
                    tfrz=TFRZ,
                    nbrdlf_dcd_brl_shrub=nbrdlf_dcd_brl_shrub,
                    nc3_arctic_grass=nc3_arctic_grass)

                # Activity factor for leaf age
                gamma_a = get_gamma_A(
                    patch.itype[p], elai240_patch[p], elai_patch[p],
                    class_num, mf;
                    ndllf_dcd_brl_tree=ndllf_dcd_brl_tree,
                    nbrdlf_dcd_trp_tree=nbrdlf_dcd_trp_tree)

                # Activity factor for CO2 (only for isoprene)
                if meg_cmp.name == "isoprene"
                    co2_ppmv = 1.0e6 * forc_pco2_grc[g] / forc_pbot_col[c]
                    gamma_c = get_gamma_C(
                        cisun_z_patch[p, 1], cisha_z_patch[p, 1],
                        forc_pbot_col[c], fsun_patch[p], co2_ppmv)
                else
                    gamma_c = 1.0
                end

                # Total scaling factor
                gamma = gamma_l * gamma_sm * gamma_a * gamma_p * gamma_t * gamma_c

                if gamma >= 0.0 && gamma < 100.0
                    vocflx_meg[imeg] = meg_cmp.coeff * epsilon * gamma *
                                       megemis_units_factor / meg_cmp.molec_weight  # moles/m2/sec

                    # History output (not weighted by landfrac) [kg/m2/sec]
                    voc.meg_flux_out[imeg, p] += epsilon * gamma * megemis_units_factor * 1.0e-3

                    # Save diagnostics for first compound
                    if imeg == 1
                        voc.gamma_out_patch[p]  = gamma
                        voc.gammaP_out_patch[p] = gamma_p
                        voc.gammaT_out_patch[p] = gamma_t
                        voc.gammaA_out_patch[p] = gamma_a
                        voc.gammaS_out_patch[p] = gamma_sm
                        voc.gammaL_out_patch[p] = gamma_l
                        voc.gammaC_out_patch[p] = gamma_c

                        voc.paru_out_patch[p]    = par_sun
                        voc.par24u_out_patch[p]  = par24_sun
                        voc.par240u_out_patch[p] = par240_sun

                        voc.para_out_patch[p]    = par_sha
                        voc.par24a_out_patch[p]  = par24_sha
                        voc.par240a_out_patch[p] = par240_sha

                        voc.alpha_out_patch[p] = alpha
                        voc.cp_out_patch[p]    = cp

                        voc.topt_out_patch[p] = topt
                        voc.Eopt_out_patch[p] = Eopt
                    end
                end
            end  # imeg loop

            # Sum up megan compound fluxes for mechanism compounds
            for imech in 1:n_mechcomps
                for ii in 1:mech_comps[imech].n_megan_comps
                    meg_idx = mech_comps[imech].megan_indices[ii]
                    voc.vocflx_patch[p, imech] += vocflx_meg[meg_idx]
                end
                voc.vocflx_tot_patch[p] += voc.vocflx_patch[p, imech]  # moles/m2/sec
            end

        end  # patch.itype > 0
    end  # p loop

    return nothing
end
