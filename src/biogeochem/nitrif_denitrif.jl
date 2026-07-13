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
# Per-column nitrification/denitrification kernel with an internal level loop.
# Each (c, j) output is fully independent (no scatter, no cross-level coupling),
# so one thread per column iterating j ascending is byte-identical to the host
# `for j; for c` nest. D_CON_G/D_CON_W table entries and physical consts
# (GRAV/RPI/SECSPDAY/TFRZ) are resolved to scalars on the host; config gates
# (use_lch4, no_frozen_nitrif_denitrif) are passed as Bool scalars.
### Device-view + scalar-bundle structs.
### Place these immediately above the @kernel definition in nitrif_denitrif.jl.

# 20 per-(column,level) output matrices, one kernel arg.
Base.@kwdef struct _NitrifOut{M}
    diffus_col::M
    r_psi_col::M
    anaerobic_frac_col::M
    k_nitr_t_vr_col::M
    k_nitr_ph_vr_col::M
    k_nitr_h2o_vr_col::M
    k_nitr_vr_col::M
    pot_f_nit_vr_col::M
    soil_bulkdensity_col::M
    smin_no3_massdens_vr_col::M
    soil_co2_prod_col::M
    fmax_denit_carbonsubstrate_vr_col::M
    fmax_denit_nitrate_vr_col::M
    f_denit_base_vr_col::M
    pot_f_denit_vr_col::M
    ratio_k1_col::M
    ratio_no3_co2_col::M
    wfps_vr_col::M
    fr_WFPS_col::M
    n2_n2o_ratio_denit_vr_col::M
end
Adapt.@adapt_structure _NitrifOut

# 17 per-(column,level) @Const input matrices, one kernel arg.
Base.@kwdef struct _NitrifIn{M}
    watsat::M
    watfc::M
    bd::M
    bsw::M
    cellorg::M
    soilpsi::M
    h2osoi_vol::M
    h2osoi_liq::M
    t_soisno::M
    o2_decomp_depth_unsat::M
    conc_o2_unsat::M
    col_dz::M
    smin_nh4_vr_col::M
    smin_no3_vr_col::M
    t_scalar_col::M
    w_scalar_col::M
    phr_vr_col::M
end
Adapt.@adapt_structure _NitrifIn

# All Float64 scalars + D_CON table entries + physical consts, at working precision T.
# isbits -> safe to materialize on Metal. No Adapt needed (no array fields).
Base.@kwdef struct _NitrifScalars{T}
    surface_tension_water::T
    rij_kro_a::T
    rij_kro_alpha::T
    rij_kro_beta::T
    rij_kro_gamma::T
    rij_kro_delta::T
    k_nitr_max_perday::T
    denit_resp_coef::T
    denit_resp_exp::T
    denit_nitrate_coef::T
    denit_nitrate_exp::T
    organic_max::T
    om_frac_sf::T
    rho_w::T
    dcong_21::T
    dcong_22::T
    dconw_21::T
    dconw_22::T
    dconw_23::T
    GRAV::T
    RPI::T
    SECSPDAY::T
    TFRZ::T
end

@kernel function _nitrifdenit_kernel!(
    out::_NitrifOut, in::_NitrifIn, @Const(mask_bgc_soilc),
    cmin::Int, cmax::Int, nlevdecomp::Int,
    use_lch4::Bool, no_frozen_nitrif_denitrif::Bool,
    p::_NitrifScalars)

    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask_bgc_soilc[c]
        # Working element type taken from an output array (matches Float64 on CPU,
        # Float32 on Metal). Every numeric literal below is converted through T so
        # no Float64 is materialized in the kernel.
        T = eltype(out.diffus_col)

        # Output array aliases (Fortran-named locals; body stays verbatim).
        diffus_col                         = out.diffus_col
        r_psi_col                          = out.r_psi_col
        anaerobic_frac_col                 = out.anaerobic_frac_col
        k_nitr_t_vr_col                    = out.k_nitr_t_vr_col
        k_nitr_ph_vr_col                   = out.k_nitr_ph_vr_col
        k_nitr_h2o_vr_col                  = out.k_nitr_h2o_vr_col
        k_nitr_vr_col                      = out.k_nitr_vr_col
        pot_f_nit_vr_col                   = out.pot_f_nit_vr_col
        soil_bulkdensity_col               = out.soil_bulkdensity_col
        smin_no3_massdens_vr_col           = out.smin_no3_massdens_vr_col
        soil_co2_prod_col                  = out.soil_co2_prod_col
        fmax_denit_carbonsubstrate_vr_col  = out.fmax_denit_carbonsubstrate_vr_col
        fmax_denit_nitrate_vr_col          = out.fmax_denit_nitrate_vr_col
        f_denit_base_vr_col                = out.f_denit_base_vr_col
        pot_f_denit_vr_col                 = out.pot_f_denit_vr_col
        ratio_k1_col                       = out.ratio_k1_col
        ratio_no3_co2_col                  = out.ratio_no3_co2_col
        wfps_vr_col                        = out.wfps_vr_col
        fr_WFPS_col                        = out.fr_WFPS_col
        n2_n2o_ratio_denit_vr_col          = out.n2_n2o_ratio_denit_vr_col

        # Input array aliases.
        watsat                = in.watsat
        watfc                 = in.watfc
        bd                    = in.bd
        bsw                   = in.bsw
        cellorg               = in.cellorg
        soilpsi               = in.soilpsi
        h2osoi_vol            = in.h2osoi_vol
        h2osoi_liq            = in.h2osoi_liq
        t_soisno              = in.t_soisno
        o2_decomp_depth_unsat = in.o2_decomp_depth_unsat
        conc_o2_unsat         = in.conc_o2_unsat
        col_dz                = in.col_dz
        smin_nh4_vr_col       = in.smin_nh4_vr_col
        smin_no3_vr_col       = in.smin_no3_vr_col
        t_scalar_col          = in.t_scalar_col
        w_scalar_col          = in.w_scalar_col
        phr_vr_col            = in.phr_vr_col

        # Scalar aliases (already at precision T).
        surface_tension_water = p.surface_tension_water
        rij_kro_a             = p.rij_kro_a
        rij_kro_alpha         = p.rij_kro_alpha
        rij_kro_beta          = p.rij_kro_beta
        rij_kro_gamma         = p.rij_kro_gamma
        rij_kro_delta         = p.rij_kro_delta
        k_nitr_max_perday     = p.k_nitr_max_perday
        denit_resp_coef       = p.denit_resp_coef
        denit_resp_exp        = p.denit_resp_exp
        denit_nitrate_coef    = p.denit_nitrate_coef
        denit_nitrate_exp     = p.denit_nitrate_exp
        organic_max           = p.organic_max
        om_frac_sf            = p.om_frac_sf
        rho_w                 = p.rho_w
        dcong_21              = p.dcong_21
        dcong_22              = p.dcong_22
        dconw_21              = p.dconw_21
        dconw_22              = p.dconw_22
        dconw_23              = p.dconw_23
        GRAV                  = p.GRAV
        RPI                   = p.RPI
        SECSPDAY              = p.SECSPDAY
        TFRZ                  = p.TFRZ

        for j in 1:nlevdecomp

            # pH set as placeholder (all soils same pH)
            pH_c = T(6.5)

            # ---- Calculate soil anoxia state ----

            # Gas diffusivity of soil at field capacity
            fc_air_frac = watsat[c, j] - watfc[c, j]
            fc_air_frac_as_frac_porosity = one(T) - watfc[c, j] / watsat[c, j]

            if use_lch4

                # Calculate organic matter fraction
                if organic_max > zero(T)
                    om_frac = min(om_frac_sf * cellorg[c, j] / organic_max, one(T))
                else
                    om_frac = one(T)
                end

                # Diffusivity after Moldrup et al. (2003)
                # Eq. 8 in Riley et al. (2011, Biogeosciences)
                diffus_moldrup = fc_air_frac^2 * fc_air_frac_as_frac_porosity^(T(3) / bsw[c, j])

                # Diffusivity after Millington & Quirk (1961)
                # Eq. 9 in Riley et al. (2011, Biogeosciences)
                diffus_millingtonquirk = fc_air_frac^(T(10) / T(3)) / watsat[c, j]^2

                # Blended unitless diffusivity
                diffus_col[c, j] =
                    om_frac * diffus_millingtonquirk +
                    (one(T) - om_frac) * diffus_moldrup

                # Calculate anoxic fraction using Rijtema and Kroess model
                # after Riley et al. (2000)
                r_min_val = T(2) * surface_tension_water / (rho_w * GRAV * abs(soilpsi[c, j]))
                r_max = T(2) * surface_tension_water / (rho_w * GRAV * T(0.1))
                r_psi_col[c, j] = sqrt(r_min_val * r_max)

                ratio_diffusivity_water_gas =
                    (dcong_21 + dcong_22 * t_soisno[c, j]) * T(1.0e-4) /
                    ((dconw_21 + dconw_22 * t_soisno[c, j] + dconw_23 * t_soisno[c, j]^2) * T(1.0e-9))

                if o2_decomp_depth_unsat[c, j] > zero(T)
                    anaerobic_frac_col[c, j] = exp(-rij_kro_a *
                        r_psi_col[c, j]^(-rij_kro_alpha) *
                        o2_decomp_depth_unsat[c, j]^(-rij_kro_beta) *
                        conc_o2_unsat[c, j]^rij_kro_gamma *
                        (h2osoi_vol[c, j] + ratio_diffusivity_water_gas * watsat[c, j])^rij_kro_delta)
                else
                    anaerobic_frac_col[c, j] = zero(T)
                end

            else
                # NITRIF_DENITRIF requires Methane model to be active,
                # otherwise diffusivity will be zeroed out here.
                anaerobic_frac_col[c, j] = zero(T)
                diffus_col[c, j] = zero(T)
            end

            # ---- Nitrification ----
            # Follows CENTURY nitrification scheme (Parton et al., 2001, 1996)

            # Assume nitrification temp function equal to the HR scalar
            k_nitr_t_vr_col[c, j] = min(t_scalar_col[c, j], one(T))

            # pH function from Parton et al. (2001, 1996)
            k_nitr_ph_vr_col[c, j] = T(0.56) + atan(RPI * T(0.45) * (T(-5) + pH_c)) / RPI

            # Moisture function — same as limits heterotrophic respiration
            k_nitr_h2o_vr_col[c, j] = w_scalar_col[c, j]

            # Nitrification rate constant (converted from 1/day to 1/s)
            k_nitr_vr_col[c, j] = k_nitr_max_perday / SECSPDAY *
                k_nitr_t_vr_col[c, j] * k_nitr_h2o_vr_col[c, j] * k_nitr_ph_vr_col[c, j]

            # Potential nitrification flux (first-order decay of ammonium)
            pot_f_nit_vr_col[c, j] = max(smin_nh4_vr_col[c, j] * k_nitr_vr_col[c, j], zero(T))

            # Limit to oxic fraction of soils
            pot_f_nit_vr_col[c, j] = pot_f_nit_vr_col[c, j] * (one(T) - anaerobic_frac_col[c, j])

            # Limit to non-frozen soil layers
            if t_soisno[c, j] <= TFRZ && no_frozen_nitrif_denitrif
                pot_f_nit_vr_col[c, j] = zero(T)
            end

            # ---- Denitrification ----

            soil_hr_vr = phr_vr_col[c, j]

            # CENTURY papers give denitrification in units of per gram soil;
            # need to convert from volumetric to mass-based units
            soil_bulkdensity_col[c, j] = bd[c, j] + h2osoi_liq[c, j] / col_dz[c, j]

            g_per_m3__to__ug_per_gsoil = T(1.0e3) / soil_bulkdensity_col[c, j]

            g_per_m3_sec__to__ug_per_gsoil_day = g_per_m3__to__ug_per_gsoil * SECSPDAY

            smin_no3_massdens_vr_col[c, j] = max(smin_no3_vr_col[c, j], zero(T)) * g_per_m3__to__ug_per_gsoil

            soil_co2_prod_col[c, j] = soil_hr_vr * g_per_m3_sec__to__ug_per_gsoil_day

            # Maximum potential denitrification rates based on heterotrophic
            # respiration rates or nitrate concentrations (del Grosso et al., 2000)
            fmax_denit_carbonsubstrate_vr_col[c, j] =
                (denit_resp_coef * (soil_co2_prod_col[c, j]^denit_resp_exp)) /
                g_per_m3_sec__to__ug_per_gsoil_day

            fmax_denit_nitrate_vr_col[c, j] =
                (denit_nitrate_coef * smin_no3_massdens_vr_col[c, j]^denit_nitrate_exp) /
                g_per_m3_sec__to__ug_per_gsoil_day

            # Limiting denitrification rate
            f_denit_base_vr_col[c, j] = max(min(fmax_denit_carbonsubstrate_vr_col[c, j],
                                                  fmax_denit_nitrate_vr_col[c, j]), zero(T))

            # Limit to non-frozen soil layers
            if t_soisno[c, j] <= TFRZ && no_frozen_nitrif_denitrif
                f_denit_base_vr_col[c, j] = zero(T)
            end

            # Limit to anoxic fraction of soils
            pot_f_denit_vr_col[c, j] = f_denit_base_vr_col[c, j] * anaerobic_frac_col[c, j]

            # ---- N2:N2O ratio from denitrification ----
            # Following Del Grosso et al. (2000)

            # Diffusivity constant (figure 6b)
            # Note: diffus_col is still the unitless relative diffusivity here
            ratio_k1_col[c, j] = max(T(1.7), T(38.4) - T(350) * diffus_col[c, j])

            # Convert diffusivity to m2/s using temperature-dependent free-air
            # diffusion rate (using O2 coefficients)
            D0 = (dcong_21 + dcong_22 * t_soisno[c, j]) * T(1.0e-4)
            diffus_col[c, j] = diffus_col[c, j] * D0

            # Ratio function (figure 7c)
            if soil_co2_prod_col[c, j] > T(1.0e-9)
                ratio_no3_co2_col[c, j] = smin_no3_massdens_vr_col[c, j] / soil_co2_prod_col[c, j]
            else
                # Function saturates at large NO3/CO2 ratios
                ratio_no3_co2_col[c, j] = T(100)
            end

            # Total water limitation function (Del Grosso et al., 2000, figure 7a)
            wfps_vr_col[c, j] = max(min(h2osoi_vol[c, j] / watsat[c, j], one(T)), zero(T)) * T(100)
            fr_WFPS_col[c, j] = max(T(0.1), T(0.015) * wfps_vr_col[c, j] - T(0.32))

            # Final N2:N2O ratio expression
            n2_n2o_ratio_denit_vr_col[c, j] = max(T(0.16) * ratio_k1_col[c, j],
                ratio_k1_col[c, j] * exp(-T(0.8) * ratio_no3_co2_col[c, j])) * fr_WFPS_col[c, j]
        end
    end
end

function nitrif_denitrif!(
    nf::SoilBiogeochemNitrogenFluxData,
    ns::SoilBiogeochemNitrogenStateData,
    cf::SoilBiogeochemCarbonFluxData,
    params::NitrifDenitrifParams,
    cn_params::CNSharedParamsData;
    mask_bgc_soilc::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevdecomp::Int,
    # Soil state arrays
    watsat::AbstractMatrix{<:Real},
    watfc::AbstractMatrix{<:Real},
    bd::AbstractMatrix{<:Real},
    bsw::AbstractMatrix{<:Real},
    cellorg::AbstractMatrix{<:Real},
    sucsat::AbstractMatrix{<:Real},
    soilpsi::AbstractMatrix{<:Real},
    # Water state arrays
    h2osoi_vol::AbstractMatrix{<:Real},
    h2osoi_liq::AbstractMatrix{<:Real},
    # Temperature
    t_soisno::AbstractMatrix{<:Real},
    # CH4 arrays (only needed when use_lch4=true)
    o2_decomp_depth_unsat::AbstractMatrix{<:Real} = zeros(0, 0),
    conc_o2_unsat::AbstractMatrix{<:Real} = zeros(0, 0),
    # Column geometry
    col_dz::AbstractMatrix{<:Real},
    # Control flags
    use_lch4::Bool = true,
    no_frozen_nitrif_denitrif::Bool = false)

    isempty(bounds) && return nothing

    # Working element type taken from an output array; all scalars are converted to
    # it so no Float64 reaches a Float32-only backend (Metal). On Float64 this is
    # byte-identical (T(x) === x).
    T = eltype(nf.diffus_col)

    organic_max = cn_params.organic_max
    rho_w = 1.0e3  # density of water (kg/m3)

    # Resolve constant-table entries (O2 coefficients) to scalars on the host.
    dcong_21 = D_CON_G[2, 1]; dcong_22 = D_CON_G[2, 2]
    dconw_21 = D_CON_W[2, 1]; dconw_22 = D_CON_W[2, 2]; dconw_23 = D_CON_W[2, 3]

    # Group output / input col matrices into device-view structs (one kernel arg each).
    out = _NitrifOut(;
        diffus_col = nf.diffus_col,
        r_psi_col = nf.r_psi_col,
        anaerobic_frac_col = nf.anaerobic_frac_col,
        k_nitr_t_vr_col = nf.k_nitr_t_vr_col,
        k_nitr_ph_vr_col = nf.k_nitr_ph_vr_col,
        k_nitr_h2o_vr_col = nf.k_nitr_h2o_vr_col,
        k_nitr_vr_col = nf.k_nitr_vr_col,
        pot_f_nit_vr_col = nf.pot_f_nit_vr_col,
        soil_bulkdensity_col = nf.soil_bulkdensity_col,
        smin_no3_massdens_vr_col = nf.smin_no3_massdens_vr_col,
        soil_co2_prod_col = nf.soil_co2_prod_col,
        fmax_denit_carbonsubstrate_vr_col = nf.fmax_denit_carbonsubstrate_vr_col,
        fmax_denit_nitrate_vr_col = nf.fmax_denit_nitrate_vr_col,
        f_denit_base_vr_col = nf.f_denit_base_vr_col,
        pot_f_denit_vr_col = nf.pot_f_denit_vr_col,
        ratio_k1_col = nf.ratio_k1_col,
        ratio_no3_co2_col = nf.ratio_no3_co2_col,
        wfps_vr_col = nf.wfps_vr_col,
        fr_WFPS_col = nf.fr_WFPS_col,
        n2_n2o_ratio_denit_vr_col = nf.n2_n2o_ratio_denit_vr_col)

    in = _NitrifIn(;
        watsat = watsat,
        watfc = watfc,
        bd = bd,
        bsw = bsw,
        cellorg = cellorg,
        soilpsi = soilpsi,
        h2osoi_vol = h2osoi_vol,
        h2osoi_liq = h2osoi_liq,
        t_soisno = t_soisno,
        o2_decomp_depth_unsat = o2_decomp_depth_unsat,
        conc_o2_unsat = conc_o2_unsat,
        col_dz = col_dz,
        smin_nh4_vr_col = ns.smin_nh4_vr_col,
        smin_no3_vr_col = ns.smin_no3_vr_col,
        t_scalar_col = cf.t_scalar_col,
        w_scalar_col = cf.w_scalar_col,
        phr_vr_col = cf.phr_vr_col)

    # All Float64 scalars + table entries + physical consts, converted to T.
    p = _NitrifScalars{T}(;
        surface_tension_water = T(params.surface_tension_water),
        rij_kro_a             = T(params.rij_kro_a),
        rij_kro_alpha         = T(params.rij_kro_alpha),
        rij_kro_beta          = T(params.rij_kro_beta),
        rij_kro_gamma         = T(params.rij_kro_gamma),
        rij_kro_delta         = T(params.rij_kro_delta),
        k_nitr_max_perday     = T(params.k_nitr_max_perday),
        denit_resp_coef       = T(params.denitrif_respiration_coefficient),
        denit_resp_exp        = T(params.denitrif_respiration_exponent),
        denit_nitrate_coef    = T(params.denitrif_nitrateconc_coefficient),
        denit_nitrate_exp     = T(params.denitrif_nitrateconc_exponent),
        organic_max           = T(organic_max),
        om_frac_sf            = T(params.om_frac_sf),
        rho_w                 = T(rho_w),
        dcong_21              = T(dcong_21),
        dcong_22              = T(dcong_22),
        dconw_21              = T(dconw_21),
        dconw_22              = T(dconw_22),
        dconw_23              = T(dconw_23),
        GRAV                  = T(GRAV),
        RPI                   = T(RPI),
        SECSPDAY              = T(SECSPDAY),
        TFRZ                  = T(TFRZ))

    # Struct-first kernel: manual backend launch + synchronize (the bundle args
    # carry no backend, so take it from an output array's first field).
    backend = _kernel_backend(out.diffus_col)
    _nitrifdenit_kernel!(backend)(
        out, in, mask_bgc_soilc,
        first(bounds), last(bounds), nlevdecomp,
        use_lch4, no_frozen_nitrif_denitrif, p;
        ndrange = last(bounds))
    KA.synchronize(backend)

    return nothing
end

# ---------------------------------------------------------------------------
# soilbiogeochem_n_state_update1! — mineral N (NH4/NO3) state update
# Ported from SoilBiogeochemNStateUpdate1Mod.F90 (the use_nitrif_denitrif=.true.,
# use_fun=.false. branch). Applies the per-step N fluxes — deposition/fixation,
# gross mineralization, immobilization, plant uptake, nitrification,
# denitrification, supplement — to smin_nh4_vr / smin_no3_vr and updates the
# diagnostic sminn_vr. (The decomp-N-pools sourcesink part of NStateUpdate1 is
# already handled in n_state_update1!.) Runs after n_state_update1! (Fortran
# CNDriverNoLeaching order, line 658).
#
# GPU: one KA kernel, one thread per column (the nlevdecomp loop runs in-thread,
# byte-identical to the host loop and race-free — each thread owns its column).
# The N-state outputs / flux inputs / N-profile inputs are grouped into three
# @adapt_structure device-view bundles to stay well under Metal's ~31-arg limit.
# ---------------------------------------------------------------------------

# Mutated NH4/NO3/total mineral-N state (per (col,level)).
Base.@kwdef struct _NSU1Out{M}
    smin_nh4_vr_col::M
    smin_no3_vr_col::M
    sminn_vr_col::M
end
Adapt.@adapt_structure _NSU1Out

# Per-step N fluxes: two per-column vectors (V) + ten per-(col,level) matrices (M).
Base.@kwdef struct _NSU1Flux{V,M}
    ndep_to_sminn_col::V
    nfix_to_sminn_col::V
    ffix_to_sminn_col::V
    gross_nmin_vr_col::M
    actual_immob_nh4_vr_col::M
    actual_immob_no3_vr_col::M
    sminn_to_plant_fun_nh4_vr_col::M
    sminn_to_plant_fun_no3_vr_col::M
    smin_nh4_to_plant_vr_col::M
    smin_no3_to_plant_vr_col::M
    f_nit_vr_col::M
    f_denit_vr_col::M
    supplement_to_sminn_vr_col::M
end
Adapt.@adapt_structure _NSU1Flux

# Deposition / fixation vertical profiles (per (col,level)).
Base.@kwdef struct _NSU1Prof{M}
    ndep_prof_col::M
    nfixation_prof_col::M
end
Adapt.@adapt_structure _NSU1Prof

@kernel function _nsu1_kernel!(out, flux, prof, @Const(mask_bgc_soilc),
        lo::Int, hi::Int, nlevdecomp::Int, dt, n2o, use_fun::Bool)
    c = @index(Global)
    @inbounds if lo <= c <= hi && mask_bgc_soilc[c]
        T = eltype(out.smin_nh4_vr_col)
        for j in 1:nlevdecomp
            # N deposition + fixation (all into NH4).
            # Fortran SoilBiogeochemNStateUpdate1Mod.F90:70-85 — under FUN the fixed
            # N goes straight to the plant, so the pool receives the FREE-LIVING
            # fixation (ffix) instead of the symbiotic one (nfix). Selecting the
            # wrong one here silently drops the entire FUN fixation input.
            _fix = use_fun ? flux.ffix_to_sminn_col[c] : flux.nfix_to_sminn_col[c]
            out.smin_nh4_vr_col[c, j] += flux.ndep_to_sminn_col[c] * dt * prof.ndep_prof_col[c, j]
            out.smin_nh4_vr_col[c, j] += _fix * dt * prof.nfixation_prof_col[c, j]
            # gross mineralization → NH4
            out.smin_nh4_vr_col[c, j] += flux.gross_nmin_vr_col[c, j] * dt
            # immobilization ← NH4 / NO3
            out.smin_nh4_vr_col[c, j] -= flux.actual_immob_nh4_vr_col[c, j] * dt
            out.smin_no3_vr_col[c, j] -= flux.actual_immob_no3_vr_col[c, j] * dt
            # plant uptake ← NH4 / NO3 (FUN: cost-based *actual* uptake)
            if use_fun
                out.smin_nh4_vr_col[c, j] -= flux.sminn_to_plant_fun_nh4_vr_col[c, j] * dt
                out.smin_no3_vr_col[c, j] -= flux.sminn_to_plant_fun_no3_vr_col[c, j] * dt
            else
                out.smin_nh4_vr_col[c, j] -= flux.smin_nh4_to_plant_vr_col[c, j] * dt
                out.smin_no3_vr_col[c, j] -= flux.smin_no3_to_plant_vr_col[c, j] * dt
            end
            # nitrification: NH4 → NO3 (minus N2O loss)
            out.smin_nh4_vr_col[c, j] -= flux.f_nit_vr_col[c, j] * dt
            out.smin_no3_vr_col[c, j] += flux.f_nit_vr_col[c, j] * dt * (one(T) - n2o)
            # denitrification ← NO3
            out.smin_no3_vr_col[c, j] -= flux.f_denit_vr_col[c, j] * dt
            # supplement (carbon-only N limitation relief) → NH4
            out.smin_nh4_vr_col[c, j] += flux.supplement_to_sminn_vr_col[c, j] * dt
            # diagnostic total
            out.sminn_vr_col[c, j] = out.smin_nh4_vr_col[c, j] + out.smin_no3_vr_col[c, j]
        end
    end
end

function soilbiogeochem_n_state_update1!(ns::SoilBiogeochemNitrogenStateData,
                                          nf::SoilBiogeochemNitrogenFluxData,
                                          st::SoilBiogeochemStateData;
                                          mask_bgc_soilc::AbstractVector{Bool},
                                          bounds_col::UnitRange{Int},
                                          nlevdecomp::Int,
                                          dt::Real,
                                          use_fun::Bool=false)
    isempty(bounds_col) && return nothing
    FT = eltype(ns.smin_nh4_vr_col)
    out = _NSU1Out(; smin_nh4_vr_col=ns.smin_nh4_vr_col,
                     smin_no3_vr_col=ns.smin_no3_vr_col,
                     sminn_vr_col=ns.sminn_vr_col)
    flux = _NSU1Flux(; ndep_to_sminn_col=nf.ndep_to_sminn_col,
                       nfix_to_sminn_col=nf.nfix_to_sminn_col,
                       ffix_to_sminn_col=nf.ffix_to_sminn_col,
                       gross_nmin_vr_col=nf.gross_nmin_vr_col,
                       actual_immob_nh4_vr_col=nf.actual_immob_nh4_vr_col,
                       actual_immob_no3_vr_col=nf.actual_immob_no3_vr_col,
                       sminn_to_plant_fun_nh4_vr_col=nf.sminn_to_plant_fun_nh4_vr_col,
                       sminn_to_plant_fun_no3_vr_col=nf.sminn_to_plant_fun_no3_vr_col,
                       smin_nh4_to_plant_vr_col=nf.smin_nh4_to_plant_vr_col,
                       smin_no3_to_plant_vr_col=nf.smin_no3_to_plant_vr_col,
                       f_nit_vr_col=nf.f_nit_vr_col,
                       f_denit_vr_col=nf.f_denit_vr_col,
                       supplement_to_sminn_vr_col=nf.supplement_to_sminn_vr_col)
    prof = _NSU1Prof(; ndep_prof_col=st.ndep_prof_col,
                       nfixation_prof_col=st.nfixation_prof_col)
    backend = _kernel_backend(out.smin_nh4_vr_col)
    _nsu1_kernel!(backend)(out, flux, prof, mask_bgc_soilc,
        first(bounds_col), last(bounds_col), nlevdecomp,
        FT(dt), FT(NITRIF_N2O_LOSS_FRAC), use_fun; ndrange = last(bounds_col))
    KA.synchronize(backend)
    return nothing
end
