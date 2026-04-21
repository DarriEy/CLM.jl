# ==========================================================================
# Ported from: src/soilbiogeochem/SoilBiogeochemNitrifDenitrifMod.F90
# Calculate nitrification and denitrification rates
#
# Public functions:
#   nitrif_denitrif_read_params!  — Read/set nitrif-denitrif parameters
#   nitrif_denitrif!              — Calculate nitrification and denitrification
# ==========================================================================

# ---------------------------------------------------------------------------
# NitrifDenitrifParams — module parameters
# Ported from params_type in SoilBiogeochemNitrifDenitrifMod.F90
# ---------------------------------------------------------------------------

"""
    NitrifDenitrifParams

Nitrification/denitrification parameters. Holds rate constants and
coefficients for the CENTURY nitrification scheme and Del Grosso
denitrification model.

Ported from `params_type` in `SoilBiogeochemNitrifDenitrifMod.F90`.
"""
Base.@kwdef mutable struct NitrifDenitrifParams
    k_nitr_max_perday                ::Float64 = 0.1       # maximum nitrification rate constant (1/day)
    surface_tension_water            ::Float64 = 0.073     # surface tension of water (J/m^2), Arah and Vinten 1995
    rij_kro_a                        ::Float64 = 1.5e-10   # Arah and Vinten 1995
    rij_kro_alpha                    ::Float64 = 1.26      # Arah and Vinten 1995
    rij_kro_beta                     ::Float64 = 0.6       # Arah and Vinten 1995
    rij_kro_gamma                    ::Float64 = 0.6       # Arah and Vinten 1995
    rij_kro_delta                    ::Float64 = 0.85      # Arah and Vinten 1995
    denitrif_respiration_coefficient ::Float64 = 0.1       # multiplier for het. resp. for max denitrif rates
    denitrif_respiration_exponent    ::Float64 = 1.3       # exponent for het. resp. for max denitrif rates
    denitrif_nitrateconc_coefficient ::Float64 = 0.1       # multiplier for nitrate conc. for max denitrif rates
    denitrif_nitrateconc_exponent    ::Float64 = 1.3       # exponent for nitrate conc. for max denitrif rates
    om_frac_sf                       ::Float64 = 1.0       # scale factor for organic matter fraction (unitless)
end

# ---------------------------------------------------------------------------
# nitrif_denitrif_read_params! — read/set parameters
# Ported from readParams in SoilBiogeochemNitrifDenitrifMod.F90
# ---------------------------------------------------------------------------

"""
    nitrif_denitrif_read_params!(params; kwargs...)

Set nitrification/denitrification parameters from keyword arguments.
Corresponds to `readParams` in the Fortran source (reads from params file).
"""
function nitrif_denitrif_read_params!(params::NitrifDenitrifParams;
                                       k_nitr_max_perday::Real,
                                       surface_tension_water::Real,
                                       rij_kro_a::Real,
                                       rij_kro_alpha::Real,
                                       rij_kro_beta::Real,
                                       rij_kro_gamma::Real,
                                       rij_kro_delta::Real,
                                       denitrif_respiration_coefficient::Real,
                                       denitrif_respiration_exponent::Real,
                                       denitrif_nitrateconc_coefficient::Real,
                                       denitrif_nitrateconc_exponent::Real,
                                       om_frac_sf::Real)
    params.k_nitr_max_perday                = k_nitr_max_perday
    params.surface_tension_water            = surface_tension_water
    params.rij_kro_a                        = rij_kro_a
    params.rij_kro_alpha                    = rij_kro_alpha
    params.rij_kro_beta                     = rij_kro_beta
    params.rij_kro_gamma                    = rij_kro_gamma
    params.rij_kro_delta                    = rij_kro_delta
    params.denitrif_respiration_coefficient = denitrif_respiration_coefficient
    params.denitrif_respiration_exponent    = denitrif_respiration_exponent
    params.denitrif_nitrateconc_coefficient = denitrif_nitrateconc_coefficient
    params.denitrif_nitrateconc_exponent    = denitrif_nitrateconc_exponent
    params.om_frac_sf                       = om_frac_sf
    return nothing
end

# ---------------------------------------------------------------------------
# nitrif_denitrif! — main calculation
# Ported from SoilBiogeochemNitrifDenitrif in
#   SoilBiogeochemNitrifDenitrifMod.F90
# ---------------------------------------------------------------------------

"""
    nitrif_denitrif!(nf, ns, cf, params, cn_params; kwargs...)

Calculate nitrification and denitrification rates.

Follows the CENTURY nitrification scheme (Parton et al., 2001, 1996) and
Del Grosso et al. (2000) denitrification model. Computes potential
nitrification (`pot_f_nit_vr`), potential denitrification (`pot_f_denit_vr`),
and the N2:N2O ratio from denitrification (`n2_n2o_ratio_denit_vr`).

Corresponds to `SoilBiogeochemNitrifDenitrif` in the Fortran source.
"""
function nitrif_denitrif!(
    nf::SoilBiogeochemNitrogenFluxData,
    ns::SoilBiogeochemNitrogenStateData,
    cf::SoilBiogeochemCarbonFluxData,
    params::NitrifDenitrifParams,
    cn_params::CNSharedParamsData;
    mask_bgc_soilc::BitVector,
    bounds::UnitRange{Int},
    nlevdecomp::Int,
    # Soil state arrays
    watsat::Matrix{<:Real},
    watfc::Matrix{<:Real},
    bd::Matrix{<:Real},
    bsw::Matrix{<:Real},
    cellorg::Matrix{<:Real},
    sucsat::Matrix{<:Real},
    soilpsi::Matrix{<:Real},
    # Water state arrays
    h2osoi_vol::Matrix{<:Real},
    h2osoi_liq::Matrix{<:Real},
    # Temperature
    t_soisno::Matrix{<:Real},
    # CH4 arrays (only needed when use_lch4=true)
    o2_decomp_depth_unsat::Matrix{<:Real} = zeros(0, 0),
    conc_o2_unsat::Matrix{<:Real} = zeros(0, 0),
    # Column geometry
    col_dz::Matrix{<:Real},
    # Control flags
    use_lch4::Bool = true,
    no_frozen_nitrif_denitrif::Bool = false)

    # Local parameter aliases
    surface_tension_water = params.surface_tension_water
    rij_kro_a             = params.rij_kro_a
    rij_kro_alpha         = params.rij_kro_alpha
    rij_kro_beta          = params.rij_kro_beta
    rij_kro_gamma         = params.rij_kro_gamma
    rij_kro_delta         = params.rij_kro_delta
    k_nitr_max_perday     = params.k_nitr_max_perday
    denit_resp_coef       = params.denitrif_respiration_coefficient
    denit_resp_exp        = params.denitrif_respiration_exponent
    denit_nitrate_coef    = params.denitrif_nitrateconc_coefficient
    denit_nitrate_exp     = params.denitrif_nitrateconc_exponent

    organic_max = cn_params.organic_max

    rho_w = 1.0e3  # density of water (kg/m3)

    for j in 1:nlevdecomp
        for c in bounds
            mask_bgc_soilc[c] || continue

            # pH set as placeholder (all soils same pH)
            pH_c = 6.5

            # ---- Calculate soil anoxia state ----

            # Gas diffusivity of soil at field capacity
            fc_air_frac = watsat[c, j] - watfc[c, j]
            fc_air_frac_as_frac_porosity = 1.0 - watfc[c, j] / watsat[c, j]

            if use_lch4

                # Calculate organic matter fraction
                if organic_max > 0.0
                    om_frac = min(params.om_frac_sf * cellorg[c, j] / organic_max, 1.0)
                else
                    om_frac = 1.0
                end

                # Diffusivity after Moldrup et al. (2003)
                # Eq. 8 in Riley et al. (2011, Biogeosciences)
                diffus_moldrup = fc_air_frac^2 * fc_air_frac_as_frac_porosity^(3.0 / bsw[c, j])

                # Diffusivity after Millington & Quirk (1961)
                # Eq. 9 in Riley et al. (2011, Biogeosciences)
                diffus_millingtonquirk = fc_air_frac^(10.0 / 3.0) / watsat[c, j]^2

                # Blended unitless diffusivity
                nf.diffus_col[c, j] =
                    om_frac * diffus_millingtonquirk +
                    (1.0 - om_frac) * diffus_moldrup

                # Calculate anoxic fraction using Rijtema and Kroess model
                # after Riley et al. (2000)
                r_min_val = 2.0 * surface_tension_water / (rho_w * GRAV * abs(soilpsi[c, j]))
                r_max = 2.0 * surface_tension_water / (rho_w * GRAV * 0.1)
                nf.r_psi_col[c, j] = sqrt(r_min_val * r_max)

                ratio_diffusivity_water_gas =
                    (D_CON_G[2, 1] + D_CON_G[2, 2] * t_soisno[c, j]) * 1.0e-4 /
                    ((D_CON_W[2, 1] + D_CON_W[2, 2] * t_soisno[c, j] + D_CON_W[2, 3] * t_soisno[c, j]^2) * 1.0e-9)

                if o2_decomp_depth_unsat[c, j] > 0.0
                    nf.anaerobic_frac_col[c, j] = exp(-rij_kro_a *
                        nf.r_psi_col[c, j]^(-rij_kro_alpha) *
                        o2_decomp_depth_unsat[c, j]^(-rij_kro_beta) *
                        conc_o2_unsat[c, j]^rij_kro_gamma *
                        (h2osoi_vol[c, j] + ratio_diffusivity_water_gas * watsat[c, j])^rij_kro_delta)
                else
                    nf.anaerobic_frac_col[c, j] = 0.0
                end

            else
                # NITRIF_DENITRIF requires Methane model to be active,
                # otherwise diffusivity will be zeroed out here.
                nf.anaerobic_frac_col[c, j] = 0.0
                nf.diffus_col[c, j] = 0.0
            end

            # ---- Nitrification ----
            # Follows CENTURY nitrification scheme (Parton et al., 2001, 1996)

            # Assume nitrification temp function equal to the HR scalar
            nf.k_nitr_t_vr_col[c, j] = min(cf.t_scalar_col[c, j], 1.0)

            # pH function from Parton et al. (2001, 1996)
            nf.k_nitr_ph_vr_col[c, j] = 0.56 + atan(RPI * 0.45 * (-5.0 + pH_c)) / RPI

            # Moisture function — same as limits heterotrophic respiration
            nf.k_nitr_h2o_vr_col[c, j] = cf.w_scalar_col[c, j]

            # Nitrification rate constant (converted from 1/day to 1/s)
            nf.k_nitr_vr_col[c, j] = k_nitr_max_perday / SECSPDAY *
                nf.k_nitr_t_vr_col[c, j] * nf.k_nitr_h2o_vr_col[c, j] * nf.k_nitr_ph_vr_col[c, j]

            # Potential nitrification flux (first-order decay of ammonium)
            nf.pot_f_nit_vr_col[c, j] = max(ns.smin_nh4_vr_col[c, j] * nf.k_nitr_vr_col[c, j], 0.0)

            # Limit to oxic fraction of soils
            nf.pot_f_nit_vr_col[c, j] = nf.pot_f_nit_vr_col[c, j] * (1.0 - nf.anaerobic_frac_col[c, j])

            # Limit to non-frozen soil layers
            if t_soisno[c, j] <= TFRZ && no_frozen_nitrif_denitrif
                nf.pot_f_nit_vr_col[c, j] = 0.0
            end

            # ---- Denitrification ----

            soil_hr_vr = cf.phr_vr_col[c, j]

            # CENTURY papers give denitrification in units of per gram soil;
            # need to convert from volumetric to mass-based units
            nf.soil_bulkdensity_col[c, j] = bd[c, j] + h2osoi_liq[c, j] / col_dz[c, j]

            g_per_m3__to__ug_per_gsoil = 1.0e3 / nf.soil_bulkdensity_col[c, j]

            g_per_m3_sec__to__ug_per_gsoil_day = g_per_m3__to__ug_per_gsoil * SECSPDAY

            nf.smin_no3_massdens_vr_col[c, j] = max(ns.smin_no3_vr_col[c, j], 0.0) * g_per_m3__to__ug_per_gsoil

            nf.soil_co2_prod_col[c, j] = soil_hr_vr * g_per_m3_sec__to__ug_per_gsoil_day

            # Maximum potential denitrification rates based on heterotrophic
            # respiration rates or nitrate concentrations (del Grosso et al., 2000)
            nf.fmax_denit_carbonsubstrate_vr_col[c, j] =
                (denit_resp_coef * (nf.soil_co2_prod_col[c, j]^denit_resp_exp)) /
                g_per_m3_sec__to__ug_per_gsoil_day

            nf.fmax_denit_nitrate_vr_col[c, j] =
                (denit_nitrate_coef * nf.smin_no3_massdens_vr_col[c, j]^denit_nitrate_exp) /
                g_per_m3_sec__to__ug_per_gsoil_day

            # Limiting denitrification rate
            nf.f_denit_base_vr_col[c, j] = max(min(nf.fmax_denit_carbonsubstrate_vr_col[c, j],
                                                     nf.fmax_denit_nitrate_vr_col[c, j]), 0.0)

            # Limit to non-frozen soil layers
            if t_soisno[c, j] <= TFRZ && no_frozen_nitrif_denitrif
                nf.f_denit_base_vr_col[c, j] = 0.0
            end

            # Limit to anoxic fraction of soils
            nf.pot_f_denit_vr_col[c, j] = nf.f_denit_base_vr_col[c, j] * nf.anaerobic_frac_col[c, j]

            # ---- N2:N2O ratio from denitrification ----
            # Following Del Grosso et al. (2000)

            # Diffusivity constant (figure 6b)
            # Note: diffus_col is still the unitless relative diffusivity here
            nf.ratio_k1_col[c, j] = max(1.7, 38.4 - 350.0 * nf.diffus_col[c, j])

            # Convert diffusivity to m2/s using temperature-dependent free-air
            # diffusion rate (using O2 coefficients)
            D0 = (D_CON_G[2, 1] + D_CON_G[2, 2] * t_soisno[c, j]) * 1.0e-4
            nf.diffus_col[c, j] = nf.diffus_col[c, j] * D0

            # Ratio function (figure 7c)
            if nf.soil_co2_prod_col[c, j] > 1.0e-9
                nf.ratio_no3_co2_col[c, j] = nf.smin_no3_massdens_vr_col[c, j] / nf.soil_co2_prod_col[c, j]
            else
                # Function saturates at large NO3/CO2 ratios
                nf.ratio_no3_co2_col[c, j] = 100.0
            end

            # Total water limitation function (Del Grosso et al., 2000, figure 7a)
            nf.wfps_vr_col[c, j] = max(min(h2osoi_vol[c, j] / watsat[c, j], 1.0), 0.0) * 100.0
            nf.fr_WFPS_col[c, j] = max(0.1, 0.015 * nf.wfps_vr_col[c, j] - 0.32)

            # Final N2:N2O ratio expression
            nf.n2_n2o_ratio_denit_vr_col[c, j] = max(0.16 * nf.ratio_k1_col[c, j],
                nf.ratio_k1_col[c, j] * exp(-0.8 * nf.ratio_no3_co2_col[c, j])) * nf.fr_WFPS_col[c, j]
        end
    end

    return nothing
end
