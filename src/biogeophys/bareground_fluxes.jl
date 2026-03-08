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
        a_coef   ::Float64 = bareground_fluxes_params.a_coef,
        a_exp    ::Float64 = bareground_fluxes_params.a_exp,
        wind_min ::Float64 = bareground_fluxes_params.wind_min)
    bareground_fluxes_params.a_coef   = a_coef
    bareground_fluxes_params.a_exp    = a_exp
    bareground_fluxes_params.wind_min = wind_min
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
        forc_q_col       ::Vector{Float64},
        forc_pbot_col    ::Vector{Float64},
        forc_th_col      ::Vector{Float64},
        forc_rho_col     ::Vector{Float64},
        forc_t_col       ::Vector{Float64},
        # Atmospheric forcing (gridcell-level)
        forc_u_grc       ::Vector{Float64},
        forc_v_grc       ::Vector{Float64},
        forc_hgt_t_grc   ::Vector{Float64},
        forc_hgt_u_grc   ::Vector{Float64},
        forc_hgt_q_grc   ::Vector{Float64};
        # Feature flags
        use_lch4         ::Bool   = false,
        z0param_method   ::String = "ZengWang2007",
        # CH4 conductance output (optional)
        grnd_ch4_cond_patch ::Vector{Float64} = Float64[])

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
    zldis  = zeros(endp)
    dth    = zeros(endp)
    dqh    = zeros(endp)
    obu    = zeros(endp)
    ur     = zeros(endp)
    um     = zeros(endp)
    temp1  = zeros(endp)
    temp12m = zeros(endp)
    temp2  = zeros(endp)
    temp22m = zeros(endp)
    ustar  = zeros(endp)
    fm     = zeros(endp)

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

        ur[p]  = max(params.wind_min, sqrt(forc_u_grc[g]^2 + forc_v_grc[g]^2))
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

        for fi in 1:fn
            p = filterp[fi]
            c = patch_data.column[p]
            g = patch_data.gridcell[p]

            tstar = temp1[p] * dth[p]
            qstar = temp2[p] * dqh[p]

            if z0param_method == "ZengWang2007"
                frictionvel.z0hg_patch[p] = frictionvel.z0mg_patch[p] /
                    exp(params.a_coef * (ustar[p] * frictionvel.z0mg_patch[p] / NU_PARAM)^params.a_exp)
            elseif z0param_method == "Meier2022"
                # After Yang et al. (2008)
                frictionvel.z0hg_patch[p] = MEIER_PARAM3 * NU_PARAM / ustar[p] *
                    exp(-BETA_PARAM * ustar[p]^0.5 * abs(tstar)^0.25)
            end

            frictionvel.z0qg_patch[p] = frictionvel.z0hg_patch[p]

            # Update the forcing heights for new roughness lengths
            frictionvel.forc_hgt_u_patch[p] = forc_hgt_u_grc[g] + frictionvel.z0mg_patch[p] + canopystate.displa_patch[p]
            frictionvel.forc_hgt_t_patch[p] = forc_hgt_t_grc[g] + frictionvel.z0hg_patch[p] + canopystate.displa_patch[p]
            frictionvel.forc_hgt_q_patch[p] = forc_hgt_q_grc[g] + frictionvel.z0qg_patch[p] + canopystate.displa_patch[p]

            thvstar = tstar * (1.0 + 0.61 * forc_q_col[c]) + 0.61 * forc_th_col[c] * qstar
            frictionvel.zeta_patch[p] = zldis[p] * VKC * GRAV * thvstar / (ustar[p]^2 * temperature.thv_col[c])

            if frictionvel.zeta_patch[p] >= 0.0  # stable
                frictionvel.zeta_patch[p] = min(frictionvel.zetamaxstable, max(frictionvel.zeta_patch[p], 0.01))
                um[p] = max(ur[p], 0.1)
            else  # unstable
                frictionvel.zeta_patch[p] = max(-100.0, min(frictionvel.zeta_patch[p], -0.01))
                wc_arg = max(-GRAV * ustar[p] * thvstar * col_data.zii[c] / temperature.thv_col[c], 0.0)
                wc = temperature.beta_col[c] * wc_arg^0.333
                um[p] = sqrt(ur[p]^2 + wc^2)
            end
            obu[p] = zldis[p] / frictionvel.zeta_patch[p]

            frictionvel.num_iter_patch[p] = Float64(iter)
        end
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
        www = min(max(www, 0.0), 1.0)

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
        waterdiagbulk.rh_ref2m_patch[p] = min(100.0,
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
