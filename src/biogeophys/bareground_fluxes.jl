# ==========================================================================
# BareGroundFluxes — ported from BareGroundFluxesMod.F90
#
# Compute sensible and latent fluxes and their derivatives with respect
# to ground temperature using ground temperatures from previous time step.
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level parameters (read from parameter file in Fortran)
# ---------------------------------------------------------------------------

Base.@kwdef mutable struct BareGroundFluxesParamsData
    a_coef   ::Float64 = 0.0   # Drag coefficient under less dense canopy (unitless)
    a_exp    ::Float64 = 0.0   # Drag exponent under less dense canopy (unitless)
    wind_min ::Float64 = 0.0   # Minimum wind speed at atmospheric forcing height (m/s)
end

const bareground_fluxes_params = BareGroundFluxesParamsData()

"""
    bareground_fluxes_read_params!(; a_coef, a_exp, wind_min)

Set BareGroundFluxes module parameters (replaces Fortran readParams).
"""
function bareground_fluxes_read_params!(;
        a_coef   ::Real = bareground_fluxes_params.a_coef,
        a_exp    ::Real = bareground_fluxes_params.a_exp,
        wind_min ::Real = bareground_fluxes_params.wind_min)
    bareground_fluxes_params.a_coef   = a_coef
    bareground_fluxes_params.a_exp    = a_exp
    bareground_fluxes_params.wind_min = wind_min
    return nothing
end

# ---------------------------------------------------------------------------
# Stability-iteration inner kernel (per-patch stability/roughness update).
#
# Kernelized form of the inner `for fi in 1:fn` loop in Phase 3 of
# bareground_fluxes!. One thread per filtered patch; each patch writes only its
# own [p] indices (no reductions). The String z0param_method branch is resolved
# on the host into z0flag (1 = ZengWang2007, 2 = Meier2022) and branched on here.
# ---------------------------------------------------------------------------
@kernel function _bgf_stability_kernel!(
        # written arrays
        z0hg_patch, z0qg_patch, forc_hgt_u_patch, forc_hgt_t_patch, forc_hgt_q_patch,
        zeta_patch, um, obu, num_iter_patch,
        # read-only arrays
        @Const(filterp), @Const(column), @Const(gridcell),
        @Const(temp1), @Const(temp2), @Const(dth), @Const(dqh),
        @Const(z0mg_patch), @Const(displa_patch),
        @Const(forc_hgt_u_grc), @Const(forc_hgt_t_grc), @Const(forc_hgt_q_grc),
        @Const(forc_q_col), @Const(forc_th_col), @Const(zldis), @Const(ustar),
        @Const(ur), @Const(thv_col), @Const(beta_col), @Const(zii),
        # scalars
        iter::Int, z0flag::Int, a_coef, a_exp, zetamaxstable,
        nu_param, meier_param3, beta_param, vkc, grav)

    fi = @index(Global)
    @inbounds begin
        p = filterp[fi]
        c = column[p]
        g = gridcell[p]

        tstar = temp1[p] * dth[p]
        qstar = temp2[p] * dqh[p]

        if z0flag == 2
            # Meier2022 — after Yang et al. (2008)
            z0hg_patch[p] = meier_param3 * nu_param / ustar[p] *
                exp(-beta_param * ustar[p]^0.5 * smooth_abs(tstar)^0.25)
        else
            # ZengWang2007
            z0hg_patch[p] = z0mg_patch[p] /
                exp(a_coef * (ustar[p] * z0mg_patch[p] / nu_param)^a_exp)
        end

        z0qg_patch[p] = z0hg_patch[p]

        # Update the forcing heights for new roughness lengths
        forc_hgt_u_patch[p] = forc_hgt_u_grc[g] + z0mg_patch[p] + displa_patch[p]
        forc_hgt_t_patch[p] = forc_hgt_t_grc[g] + z0hg_patch[p] + displa_patch[p]
        forc_hgt_q_patch[p] = forc_hgt_q_grc[g] + z0qg_patch[p] + displa_patch[p]

        thvstar = tstar * (1.0 + 0.61 * forc_q_col[c]) + 0.61 * forc_th_col[c] * qstar
        zeta_patch[p] = zldis[p] * vkc * grav * thvstar / (ustar[p]^2 * thv_col[c])

        if zeta_patch[p] >= 0.0  # stable
            zeta_patch[p] = smooth_clamp(zeta_patch[p], 0.01, zetamaxstable)
            um[p] = smooth_max(ur[p], 0.1)
        else  # unstable
            zeta_patch[p] = smooth_clamp(zeta_patch[p], -100.0, -0.01)
            wc_arg = smooth_max(-grav * ustar[p] * thvstar * zii[c] / thv_col[c], 0.0)
            wc = beta_col[c] * wc_arg^0.333
            um[p] = sqrt(ur[p]^2 + wc^2)
        end
        obu[p] = zldis[p] / zeta_patch[p]

        num_iter_patch[p] = Float64(iter)
    end
end

"""
    bgf_stability_update!(frictionvel, canopystate, temperature, col_data,
                          patch_data, filterp, fn, temp1, temp2, dth, dqh,
                          zldis, ustar, ur, um, obu, forc_q_col, forc_th_col,
                          forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc,
                          iter, z0flag, params)

Launch the bare-ground stability inner-loop kernel over the `fn` filtered
patches. Backend-agnostic (CPU loop or GPU); one thread per filtered patch.
Replaces the inner `for fi in 1:fn` stability/roughness update in Phase 3.
"""
function bgf_stability_update!(frictionvel, canopystate, temperature, col_data,
        patch_data, filterp, fn::Int,
        temp1, temp2, dth, dqh, zldis, ustar, ur, um, obu,
        forc_q_col, forc_th_col, forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc,
        iter::Int, z0flag::Int, params)
    _launch!(_bgf_stability_kernel!, frictionvel.z0hg_patch,
        frictionvel.z0qg_patch, frictionvel.forc_hgt_u_patch,
        frictionvel.forc_hgt_t_patch, frictionvel.forc_hgt_q_patch,
        frictionvel.zeta_patch, um, obu, frictionvel.num_iter_patch,
        filterp, patch_data.column, patch_data.gridcell,
        temp1, temp2, dth, dqh,
        frictionvel.z0mg_patch, canopystate.displa_patch,
        forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc,
        forc_q_col, forc_th_col, zldis, ustar, ur,
        temperature.thv_col, temperature.beta_col, col_data.zii,
        iter, z0flag, params.a_coef, params.a_exp, frictionvel.zetamaxstable,
        NU_PARAM, MEIER_PARAM3, BETA_PARAM, VKC, GRAV;
        ndrange = fn)
    return nothing
end

# ---------------------------------------------------------------------------
# Main subroutine
# ---------------------------------------------------------------------------

"""
    bareground_fluxes!(...)

Compute sensible and latent fluxes and their derivatives with respect to
ground temperature using ground temperatures from previous time step.

Ported from: BareGroundFluxesMod.F90 :: BareGroundFluxes
"""
function bareground_fluxes!(
        # Data structures
        canopystate      ::CanopyStateData,
        energyflux       ::EnergyFluxData,
        frictionvel      ::FrictionVelocityData,
        temperature      ::TemperatureData,
        soilstate        ::SoilStateData,
        waterfluxbulk    ::WaterFluxBulkData,
        waterstatebulk   ::WaterStateBulkData,
        waterdiagbulk    ::WaterDiagnosticBulkData,
        photosyns        ::PhotosynthesisData,
        patch_data       ::PatchData,
        col_data         ::ColumnData,
        lun_data         ::LandunitData,
        mask_noexposedvegp ::BitVector,
        bounds_patch     ::UnitRange{Int},
        # Atmospheric forcing (column-level)
        forc_q_col       ::Vector{<:Real},
        forc_pbot_col    ::Vector{<:Real},
        forc_th_col      ::Vector{<:Real},
        forc_rho_col     ::Vector{<:Real},
        forc_t_col       ::Vector{<:Real},
        # Atmospheric forcing (gridcell-level)
        forc_u_grc       ::Vector{<:Real},
        forc_v_grc       ::Vector{<:Real},
        forc_hgt_t_grc   ::Vector{<:Real},
        forc_hgt_u_grc   ::Vector{<:Real},
        forc_hgt_q_grc   ::Vector{<:Real};
        # Feature flags
        use_lch4         ::Bool   = false,
        z0param_method   ::String = "ZengWang2007",
        # CH4 conductance output (optional)
        grnd_ch4_cond_patch ::Vector{<:Real} = Float64[])

    np   = length(bounds_patch)
    begp = first(bounds_patch)
    endp = last(bounds_patch)
    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno

    niters = 3  # maximum number of iterations for surface temperature

    params = bareground_fluxes_params

    # --- Build integer filter from mask ---
    filterp = zeros(Int, np)
    fn = 0
    for p in bounds_patch
        if mask_noexposedvegp[p]
            fn += 1
            filterp[fn] = p
        end
    end

    # --- Local work arrays (patch-indexed) ---
    FT = eltype(forc_t_col)
    zldis  = zeros(FT, endp)
    dth    = zeros(FT, endp)
    dqh    = zeros(FT, endp)
    obu    = zeros(FT, endp)
    ur     = zeros(FT, endp)
    um     = zeros(FT, endp)
    temp1  = zeros(FT, endp)
    temp12m = zeros(FT, endp)
    temp2  = zeros(FT, endp)
    temp22m = zeros(FT, endp)
    ustar  = zeros(FT, endp)
    fm     = zeros(FT, endp)

    # =========================================================================
    # Phase 1: Initialization — simple settings for no-exposed-veg patches
    # =========================================================================

    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]

        energyflux.btran_patch[p] = 0.0
        temperature.t_veg_patch[p] = forc_t_col[c]
        cf_bare = forc_pbot_col[c] / (RGAS * temperature.thm_patch[p]) * 1.0e06
        photosyns.rssun_patch[p] = 1.0 / 1.0e15 * cf_bare
        photosyns.rssha_patch[p] = 1.0 / 1.0e15 * cf_bare
        for j in 1:nlevgrnd
            soilstate.rootr_patch[p, j] = 0.0
            energyflux.rresis_patch[p, j] = 0.0
        end
    end

    # =========================================================================
    # Phase 2: Compute initial state and stability parameters
    # =========================================================================

    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]
        g = patch_data.gridcell[p]

        # Initialization variables
        canopystate.displa_patch[p] = 0.0
        frictionvel.z0mv_patch[p]   = 0.0
        frictionvel.z0hv_patch[p]   = 0.0
        frictionvel.z0qv_patch[p]   = 0.0
        energyflux.dlrad_patch[p]   = 0.0
        energyflux.ulrad_patch[p]   = 0.0
        energyflux.dhsdt_canopy_patch[p] = 0.0
        energyflux.eflx_sh_stem_patch[p] = 0.0

        ur[p]  = smooth_max(params.wind_min, sqrt(forc_u_grc[g]^2 + forc_v_grc[g]^2))
        dth[p] = temperature.thm_patch[p] - temperature.t_grnd_col[c]
        dqh[p] = forc_q_col[c] - waterdiagbulk.qg_col[c]
        dthv   = dth[p] * (1.0 + 0.61 * forc_q_col[c]) + 0.61 * forc_th_col[c] * dqh[p]
        zldis[p] = frictionvel.forc_hgt_u_patch[p]

        # Copy column roughness to local patch-level arrays
        frictionvel.z0mg_patch[p] = frictionvel.z0mg_col[c]
        frictionvel.z0hg_patch[p] = frictionvel.z0hg_col[c]
        frictionvel.z0qg_patch[p] = frictionvel.z0qg_col[c]

        # Initialize Monin-Obukhov length and wind speed
        (um_val, obu_val) = monin_obuk_ini(frictionvel.zetamaxstable,
            ur[p], temperature.thv_col[c], dthv, zldis[p], frictionvel.z0mg_patch[p])
        um[p]  = um_val
        obu[p] = obu_val

        # Initialize iteration counter
        frictionvel.num_iter_patch[p] = 0.0
    end

    # =========================================================================
    # Phase 3: Stability iteration
    # =========================================================================

    # Resolve the String z0param_method to an Int flag on the host (a kernel
    # cannot compare Strings): 1 = ZengWang2007 (default), 2 = Meier2022.
    z0flag = z0param_method == "Meier2022" ? 2 : 1

    for iter in 1:niters

        # Friction velocity calculation
        filt_arr  = filterp[1:fn]
        disp_vec  = canopystate.displa_patch
        z0m_vec   = frictionvel.z0mg_patch
        z0h_vec   = frictionvel.z0hg_patch
        z0q_vec   = frictionvel.z0qg_patch

        friction_velocity!(frictionvel, fn, filt_arr,
            disp_vec, z0m_vec, z0h_vec, z0q_vec,
            obu, iter, ur, um, ustar,
            temp1, temp2, temp12m, temp22m, fm)

        bgf_stability_update!(frictionvel, canopystate, temperature, col_data,
            patch_data, filterp, fn,
            temp1, temp2, dth, dqh, zldis, ustar, ur, um, obu,
            forc_q_col, forc_th_col, forc_hgt_u_grc, forc_hgt_t_grc, forc_hgt_q_grc,
            iter, z0flag, params)
    end  # end stability iteration

    # =========================================================================
    # Phase 4: Post-iteration — compute final fluxes and diagnostics
    # =========================================================================

    for fi in 1:fn
        p = filterp[fi]
        c = patch_data.column[p]
        g = patch_data.gridcell[p]
        l = patch_data.landunit[p]

        # Determine aerodynamic resistances
        ram  = 1.0 / (ustar[p]^2 / um[p])
        rah  = 1.0 / (temp1[p] * ustar[p])
        raw  = 1.0 / (temp2[p] * ustar[p])
        raih = forc_rho_col[c] * CPAIR / rah

        if use_lch4 && !isempty(grnd_ch4_cond_patch)
            grnd_ch4_cond_patch[p] = 1.0 / raw
        end

        # Soil evaporation resistance
        # Fortran index 1 = first soil layer = Julia index 1 + nlevsno
        www = (waterstatebulk.ws.h2osoi_liq_col[c, 1 + nlevsno] / DENH2O +
               waterstatebulk.ws.h2osoi_ice_col[c, 1 + nlevsno] / DENICE) /
              col_data.dz[c, 1 + nlevsno] / soilstate.watsat_col[c, 1]
        www = smooth_clamp(www, 0.0, 1.0)

        # Soil evaporation: beta or resistance method
        if dqh[p] > 0.0  # dew (beta not applied)
            raiw = forc_rho_col[c] / raw
        else
            if do_soilevap_beta()
                # Lee and Pielke 1992 beta
                raiw = soilstate.soilbeta_col[c] * forc_rho_col[c] / raw
            end
            if do_soil_resistance_sl14()
                # Swenson & Lawrence 2014 soil resistance
                raiw = forc_rho_col[c] / (raw + soilstate.soilresis_col[c])
            end
        end

        frictionvel.ram1_patch[p] = ram  # pass value to global variable
        frictionvel.ustar_patch[p] = ustar[p]

        # Derivative of fluxes with respect to ground temperature
        energyflux.cgrnds_patch[p] = raih
        energyflux.cgrndl_patch[p] = raiw * waterdiagbulk.dqgdT_col[c]
        energyflux.cgrnd_patch[p]  = energyflux.cgrnds_patch[p] +
            energyflux.htvp_col[c] * energyflux.cgrndl_patch[p]

        # Surface fluxes of momentum, sensible and latent heat
        energyflux.taux_patch[p]          = -forc_rho_col[c] * forc_u_grc[g] / ram
        energyflux.tauy_patch[p]          = -forc_rho_col[c] * forc_v_grc[g] / ram
        energyflux.eflx_sh_grnd_patch[p]  = -raih * dth[p]
        energyflux.eflx_sh_tot_patch[p]   = energyflux.eflx_sh_grnd_patch[p]

        # Compute sensible heat fluxes individually
        energyflux.eflx_sh_snow_patch[p]   = -raih * (temperature.thm_patch[p] -
            temperature.t_soisno_col[c, col_data.snl[c] + 1 + nlevsno])
        energyflux.eflx_sh_soil_patch[p]   = -raih * (temperature.thm_patch[p] -
            temperature.t_soisno_col[c, 1 + nlevsno])
        energyflux.eflx_sh_h2osfc_patch[p] = -raih * (temperature.thm_patch[p] -
            temperature.t_h2osfc_col[c])

        # Water fluxes from soil
        waterfluxbulk.wf.qflx_tran_veg_patch[p] = 0.0
        waterfluxbulk.wf.qflx_evap_veg_patch[p] = 0.0
        waterfluxbulk.wf.qflx_evap_soi_patch[p] = -raiw * dqh[p]
        waterfluxbulk.wf.qflx_evap_tot_patch[p] = waterfluxbulk.wf.qflx_evap_soi_patch[p]

        # Compute latent heat fluxes individually
        waterfluxbulk.qflx_ev_snow_patch[p]   = -raiw * (forc_q_col[c] - waterdiagbulk.qg_snow_col[c])
        waterfluxbulk.qflx_ev_soil_patch[p]   = -raiw * (forc_q_col[c] - waterdiagbulk.qg_soil_col[c])
        waterfluxbulk.qflx_ev_h2osfc_patch[p] = -raiw * (forc_q_col[c] - waterdiagbulk.qg_h2osfc_col[c])

        # 2 m height air temperature
        temperature.t_ref2m_patch[p] = temperature.thm_patch[p] +
            temp1[p] * dth[p] * (1.0 / temp12m[p] - 1.0 / temp1[p])

        # 2 m height specific humidity
        waterdiagbulk.q_ref2m_patch[p] = forc_q_col[c] +
            temp2[p] * dqh[p] * (1.0 / temp22m[p] - 1.0 / temp2[p])

        # 2 m height relative humidity
        (qsat_ref2m, e_ref2m, _, _) = qsat(temperature.t_ref2m_patch[p], forc_pbot_col[c])
        waterdiagbulk.rh_ref2m_patch[p] = smooth_min(100.0,
            waterdiagbulk.q_ref2m_patch[p] / qsat_ref2m * 100.0)

        if lun_data.itype[l] == ISTSOIL || lun_data.itype[l] == ISTCROP
            waterdiagbulk.rh_ref2m_r_patch[p] = waterdiagbulk.rh_ref2m_patch[p]
            temperature.t_ref2m_r_patch[p] = temperature.t_ref2m_patch[p]
        end

        frictionvel.kbm1_patch[p] = log(frictionvel.z0mg_patch[p] / frictionvel.z0hg_patch[p])

        # Copy local patch ground roughness back to column arrays for history output
        frictionvel.z0hg_col[c] = frictionvel.z0hg_patch[p]
        frictionvel.z0qg_col[c] = frictionvel.z0qg_patch[p]

        # Human heat stress indices (stub — HumanIndexMod not yet ported)
        # Would be called here in full implementation
    end

    return nothing
end
