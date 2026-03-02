# ==========================================================================
# UrbanFluxes — ported from UrbanFluxesMod.F90
#
# Calculate turbulent fluxes for urban landunit (roof, sunwall, shadewall,
# pervious and impervious road).
# ==========================================================================

# ---------------------------------------------------------------------------
# Module-level parameters (read from parameter file in Fortran)
# ---------------------------------------------------------------------------

Base.@kwdef mutable struct UrbanFluxesParamsData
    wind_min::Float64 = 0.0  # Minimum wind speed at atmospheric forcing height (m/s)
end

const urban_fluxes_params = UrbanFluxesParamsData()

"""
    urban_fluxes_read_params!(; wind_min)

Set UrbanFluxes module parameters (replaces Fortran readParams).
"""
function urban_fluxes_read_params!(; wind_min::Float64 = urban_fluxes_params.wind_min)
    urban_fluxes_params.wind_min = wind_min
    return nothing
end

# ---------------------------------------------------------------------------
# simple_wasteheatfromac — Calculate waste heat from AC (simple CLM4.5 method)
# ---------------------------------------------------------------------------

"""
    simple_wasteheatfromac(eflx_urban_ac, eflx_urban_heat)
        -> (eflx_wasteheat, eflx_heat_from_ac)

Calculate waste heat from air-conditioning with the simpler method (CLM4.5).

Ported from: UrbanFluxesMod.F90 :: simple_wasteheatfromac
"""
function simple_wasteheatfromac(eflx_urban_ac::Float64, eflx_urban_heat::Float64)
    # wasteheat from heating/cooling
    if urban_ctrl.urban_hac == URBAN_WASTEHEAT_ON
        eflx_wasteheat = AC_WASTEHEAT_FACTOR * eflx_urban_ac +
                          HT_WASTEHEAT_FACTOR * eflx_urban_heat
    else
        eflx_wasteheat = 0.0
    end

    # If air conditioning on, always replace heat removed with heat into canyon
    if urban_ctrl.urban_hac == URBAN_HAC_ON || urban_ctrl.urban_hac == URBAN_WASTEHEAT_ON
        eflx_heat_from_ac = abs(eflx_urban_ac)
    else
        eflx_heat_from_ac = 0.0
    end

    return (eflx_wasteheat, eflx_heat_from_ac)
end

# ---------------------------------------------------------------------------
# wasteheat! — Calculate wasteheat from urban heating/cooling
# ---------------------------------------------------------------------------

"""
    wasteheat!(energyflux, lun_data, num_urbanl, filter_urbanl,
               eflx_wasteheat_roof, eflx_wasteheat_sunwall,
               eflx_wasteheat_shadewall, eflx_heat_from_ac_roof,
               eflx_heat_from_ac_sunwall, eflx_heat_from_ac_shadewall)

Calculate the wasteheat flux from urban heating or air-conditioning.

Ported from: UrbanFluxesMod.F90 :: wasteheat
"""
function wasteheat!(
        energyflux       ::EnergyFluxData,
        lun_data         ::LandunitData,
        num_urbanl       ::Int,
        filter_urbanl    ::Vector{Int},
        eflx_wasteheat_roof      ::Vector{Float64},
        eflx_wasteheat_sunwall   ::Vector{Float64},
        eflx_wasteheat_shadewall ::Vector{Float64},
        eflx_heat_from_ac_roof      ::Vector{Float64},
        eflx_heat_from_ac_sunwall   ::Vector{Float64},
        eflx_heat_from_ac_shadewall ::Vector{Float64})

    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        g = lun_data.gridcell[l]

        if is_simple_build_temp()
            # Total waste heat and heat from AC is sum of heat for walls and roofs
            # accounting for different surface areas
            energyflux.eflx_wasteheat_lun[l] = lun_data.wtlunit_roof[l] * eflx_wasteheat_roof[l] +
                (1.0 - lun_data.wtlunit_roof[l]) * (lun_data.canyon_hwr[l] * (eflx_wasteheat_sunwall[l] +
                eflx_wasteheat_shadewall[l]))

        elseif is_prog_build_temp()
            # wasteheat from heating/cooling
            if urban_ctrl.urban_hac == URBAN_WASTEHEAT_ON
                energyflux.eflx_wasteheat_lun[l] = AC_WASTEHEAT_FACTOR * energyflux.eflx_urban_ac_lun[l] +
                    HT_WASTEHEAT_FACTOR * energyflux.eflx_urban_heat_lun[l]
            else
                energyflux.eflx_wasteheat_lun[l] = 0.0
            end
        end

        # Limit wasteheat to ensure no unrealistically strong positive feedbacks
        energyflux.eflx_wasteheat_lun[l] = min(energyflux.eflx_wasteheat_lun[l], WASTEHEAT_LIMIT)

        if is_simple_build_temp()
            energyflux.eflx_heat_from_ac_lun[l] = lun_data.wtlunit_roof[l] * eflx_heat_from_ac_roof[l] +
                (1.0 - lun_data.wtlunit_roof[l]) * (lun_data.canyon_hwr[l] * (eflx_heat_from_ac_sunwall[l] +
                eflx_heat_from_ac_shadewall[l]))

        elseif is_prog_build_temp()
            # If air conditioning on, always replace heat removed with heat into canyon
            if urban_ctrl.urban_hac == URBAN_HAC_ON || urban_ctrl.urban_hac == URBAN_WASTEHEAT_ON
                energyflux.eflx_heat_from_ac_lun[l] = abs(energyflux.eflx_urban_ac_lun[l])
            else
                energyflux.eflx_heat_from_ac_lun[l] = 0.0
            end
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# calc_simple_internal_building_temp! — Simple building temp (CLM4.5)
# ---------------------------------------------------------------------------

"""
    calc_simple_internal_building_temp!(temperature, col_data, lun_data,
        num_urbanc, filter_urbanc, num_urbanl, filter_urbanl)

Calculate the internal building temperature by simpler method (CLM4.5).

Ported from: UrbanFluxesMod.F90 :: calc_simple_internal_building_temp
"""
function calc_simple_internal_building_temp!(
        temperature   ::TemperatureData,
        col_data      ::ColumnData,
        lun_data      ::LandunitData,
        num_urbanc    ::Int,
        filter_urbanc ::Vector{Int},
        num_urbanl    ::Int,
        filter_urbanl ::Vector{Int})

    nlevurb = varpar.nlevurb

    nl = length(lun_data.gridcell)
    t_sunwall_innerl   = zeros(nl)
    t_shadewall_innerl = zeros(nl)
    t_roof_innerl      = zeros(nl)

    # Gather inner layer temperatures from urban columns
    for fc in 1:num_urbanc
        c = filter_urbanc[fc]
        l = col_data.landunit[c]

        if col_data.itype[c] == ICOL_ROOF
            t_roof_innerl[l] = temperature.t_soisno_col[c, nlevurb]
        elseif col_data.itype[c] == ICOL_SUNWALL
            t_sunwall_innerl[l] = temperature.t_soisno_col[c, nlevurb]
        elseif col_data.itype[c] == ICOL_SHADEWALL
            t_shadewall_innerl[l] = temperature.t_soisno_col[c, nlevurb]
        end
    end

    # Calculate internal building temperature
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]

        lngth_roof = (lun_data.ht_roof[l] / lun_data.canyon_hwr[l]) *
            lun_data.wtlunit_roof[l] / (1.0 - lun_data.wtlunit_roof[l])
        temperature.t_building_lun[l] = (lun_data.ht_roof[l] * (t_shadewall_innerl[l] + t_sunwall_innerl[l]) +
            lngth_roof * t_roof_innerl[l]) / (2.0 * lun_data.ht_roof[l] + lngth_roof)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Main subroutine: urban_fluxes!
# ---------------------------------------------------------------------------

"""
    urban_fluxes!(...)

Compute turbulent and momentum fluxes from urban canyon (roof, sunwall,
shadewall, pervious and impervious road).

Ported from: UrbanFluxesMod.F90 :: UrbanFluxes
"""
function urban_fluxes!(
        # Data structures
        energyflux       ::EnergyFluxData,
        frictionvel      ::FrictionVelocityData,
        temperature      ::TemperatureData,
        soilstate        ::SoilStateData,
        urbanparams      ::UrbanParamsData,
        waterfluxbulk    ::WaterFluxBulkData,
        waterstatebulk   ::WaterStateBulkData,
        waterdiagbulk    ::WaterDiagnosticBulkData,
        patch_data       ::PatchData,
        col_data         ::ColumnData,
        lun_data         ::LandunitData,
        # Filters
        num_nourbanl     ::Int,
        filter_nourbanl  ::Vector{Int},
        num_urbanl       ::Int,
        filter_urbanl    ::Vector{Int},
        num_urbanc       ::Int,
        filter_urbanc    ::Vector{Int},
        num_urbanp       ::Int,
        filter_urbanp    ::Vector{Int},
        # Bounds
        bounds_lun       ::UnitRange{Int},
        bounds_col       ::UnitRange{Int},
        bounds_patch     ::UnitRange{Int},
        # Atmospheric forcing (gridcell-level)
        forc_t_grc       ::Vector{Float64},
        forc_th_grc      ::Vector{Float64},
        forc_rho_grc     ::Vector{Float64},
        forc_q_grc       ::Vector{Float64},
        forc_pbot_grc    ::Vector{Float64},
        forc_u_grc       ::Vector{Float64},
        forc_v_grc       ::Vector{Float64};
        # Optional time step info
        dtime            ::Float64 = 3600.0,
        nstep            ::Int     = 1)

    params = urban_fluxes_params
    nlevgrnd = varpar.nlevgrnd

    begl = first(bounds_lun)
    endl = last(bounds_lun)
    begc = first(bounds_col)
    endc = last(bounds_col)
    begp = first(bounds_patch)
    endp = last(bounds_patch)

    lapse_rate = 0.0098   # Dry adiabatic lapse rate (K/m)
    niters = 3            # maximum number of iterations for surface temperature

    # --- Local work arrays (landunit-indexed) ---
    nl = endl
    nc = endc
    np_local = endp

    canyontop_wind      = zeros(nl)
    canyon_u_wind       = zeros(nl)
    canyon_wind         = zeros(nl)
    canyon_resistance   = zeros(nl)
    ur                  = zeros(nl)
    ustar_loc           = zeros(nl)
    ramu                = zeros(nl)
    rahu                = zeros(nl)
    rawu                = zeros(nl)
    temp1               = zeros(nl)
    temp12m             = zeros(nl)
    temp2               = zeros(nl)
    temp22m             = zeros(nl)
    thm_g               = zeros(nl)
    thv_g               = zeros(nl)
    dth                 = zeros(nl)
    dqh                 = zeros(nl)
    zldis_arr           = zeros(nl)
    zeta_lunit          = zeros(nl)
    um                  = zeros(nl)
    obu                 = zeros(nl)
    taf_numer           = zeros(nl)
    taf_denom           = zeros(nl)
    qaf_numer           = zeros(nl)
    qaf_denom           = zeros(nl)
    wtas                = zeros(nl)
    wtaq                = zeros(nl)
    wts_sum             = zeros(nl)
    wtq_sum             = zeros(nl)
    beta_arr            = zeros(nl)
    zii_arr             = zeros(nl)
    fm_arr              = zeros(nl)

    wtus_col            = zeros(nc)
    wtuq_col            = zeros(nc)

    wtus_roof           = zeros(nl)
    wtuq_roof           = zeros(nl)
    wtus_road_perv      = zeros(nl)
    wtuq_road_perv      = zeros(nl)
    wtus_road_imperv    = zeros(nl)
    wtuq_road_imperv    = zeros(nl)
    wtus_sunwall        = zeros(nl)
    wtuq_sunwall        = zeros(nl)
    wtus_shadewall      = zeros(nl)
    wtuq_shadewall      = zeros(nl)

    wtus_roof_unscl        = zeros(nl)
    wtuq_roof_unscl        = zeros(nl)
    wtus_road_perv_unscl   = zeros(nl)
    wtuq_road_perv_unscl   = zeros(nl)
    wtus_road_imperv_unscl = zeros(nl)
    wtuq_road_imperv_unscl = zeros(nl)
    wtus_sunwall_unscl     = zeros(nl)
    wtuq_sunwall_unscl     = zeros(nl)
    wtus_shadewall_unscl   = zeros(nl)
    wtuq_shadewall_unscl   = zeros(nl)

    eflx_sh_grnd_scale     = zeros(np_local)
    qflx_evap_soi_scale    = zeros(np_local)

    eflx_wasteheat_roof      = zeros(nl)
    eflx_wasteheat_sunwall   = zeros(nl)
    eflx_wasteheat_shadewall = zeros(nl)
    eflx_heat_from_ac_roof      = zeros(nl)
    eflx_heat_from_ac_sunwall   = zeros(nl)
    eflx_heat_from_ac_shadewall = zeros(nl)

    eflx_arr       = zeros(nl)
    qflx_arr       = zeros(nl)
    eflx_scale_arr = zeros(nl)
    qflx_scale_arr = zeros(nl)
    eflx_err       = zeros(nl)
    qflx_err       = zeros(nl)

    # =========================================================================
    # Set restart fields for non-urban landunits
    # =========================================================================
    for fl in 1:num_nourbanl
        l = filter_nourbanl[fl]
        temperature.taf_lun[l] = SPVAL
        waterdiagbulk.qaf_lun[l] = SPVAL
    end

    # =========================================================================
    # Set constants
    # =========================================================================
    for l in begl:endl
        beta_arr[l] = 1.0
        zii_arr[l]  = 1000.0
    end

    # =========================================================================
    # Compute canyontop wind using Masson (2000)
    # =========================================================================
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        g = lun_data.gridcell[l]

        # Error checks
        if lun_data.ht_roof[l] - lun_data.z_d_town[l] <= lun_data.z_0_town[l]
            error("aerodynamic parameter error in urban_fluxes!: ht_roof - z_d_town <= z_0_town " *
                  "ht_roof=$(lun_data.ht_roof[l]) z_d_town=$(lun_data.z_d_town[l]) z_0_town=$(lun_data.z_0_town[l])")
        end
        if frictionvel.forc_hgt_u_patch[lun_data.patchi[l]] - lun_data.z_d_town[l] <= lun_data.z_0_town[l]
            error("aerodynamic parameter error in urban_fluxes!: forc_hgt_u - z_d_town <= z_0_town " *
                  "forc_hgt_u=$(frictionvel.forc_hgt_u_patch[lun_data.patchi[l]]) z_d_town=$(lun_data.z_d_town[l]) z_0_town=$(lun_data.z_0_town[l])")
        end

        # Magnitude of atmospheric wind
        ur[l] = max(params.wind_min, sqrt(forc_u_grc[g]^2 + forc_v_grc[g]^2))

        # Canyon top wind
        canyontop_wind[l] = ur[l] *
            log((lun_data.ht_roof[l] - lun_data.z_d_town[l]) / lun_data.z_0_town[l]) /
            log((frictionvel.forc_hgt_u_patch[lun_data.patchi[l]] - lun_data.z_d_town[l]) / lun_data.z_0_town[l])

        # U component of canyon wind
        if lun_data.canyon_hwr[l] < 0.5  # isolated roughness flow
            canyon_u_wind[l] = canyontop_wind[l] * exp(-0.5 * lun_data.canyon_hwr[l] *
                (1.0 - (urbanparams.wind_hgt_canyon[l] / lun_data.ht_roof[l])))
        elseif lun_data.canyon_hwr[l] < 1.0  # wake interference flow
            canyon_u_wind[l] = canyontop_wind[l] * (1.0 + 2.0 * (2.0 / RPI - 1.0) *
                (lun_data.ht_roof[l] / (lun_data.ht_roof[l] / lun_data.canyon_hwr[l]) - 0.5)) *
                exp(-0.5 * lun_data.canyon_hwr[l] * (1.0 - (urbanparams.wind_hgt_canyon[l] / lun_data.ht_roof[l])))
        else  # skimming flow
            canyon_u_wind[l] = canyontop_wind[l] * (2.0 / RPI) *
                exp(-0.5 * lun_data.canyon_hwr[l] * (1.0 - (urbanparams.wind_hgt_canyon[l] / lun_data.ht_roof[l])))
        end
    end

    # =========================================================================
    # Compute fluxes — follows CLM approach for bare soils
    # =========================================================================
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        g = lun_data.gridcell[l]

        thm_g[l] = forc_t_grc[g] + lapse_rate * frictionvel.forc_hgt_t_patch[lun_data.patchi[l]]
        thv_g[l] = forc_th_grc[g] * (1.0 + 0.61 * forc_q_grc[g])
        dth[l]   = thm_g[l] - temperature.taf_lun[l]
        dqh[l]   = forc_q_grc[g] - waterdiagbulk.qaf_lun[l]
        dthv     = dth[l] * (1.0 + 0.61 * forc_q_grc[g]) + 0.61 * forc_th_grc[g] * dqh[l]
        zldis_arr[l] = frictionvel.forc_hgt_u_patch[lun_data.patchi[l]] - lun_data.z_d_town[l]

        # Initialize Monin-Obukhov length and wind speed
        (um_val, obu_val) = monin_obuk_ini(frictionvel.zetamaxstable,
            ur[l], thv_g[l], dthv, zldis_arr[l], lun_data.z_0_town[l])
        um[l]  = um_val
        obu[l] = obu_val
    end

    # =========================================================================
    # Stability iteration
    # =========================================================================
    for iter in 1:niters

        # Get friction velocity
        if num_urbanl > 0
            z_d_vec  = lun_data.z_d_town
            z_0_vec  = lun_data.z_0_town

            friction_velocity!(frictionvel, num_urbanl, filter_urbanl,
                z_d_vec, z_0_vec, z_0_vec, z_0_vec,
                obu, iter, ur, um, ustar_loc,
                temp1, temp2, temp12m, temp22m, fm_arr;
                landunit_index=true,
                lun_gridcell=lun_data.gridcell,
                lun_patchi=lun_data.patchi,
                lun_patchf=lun_data.patchf)
        end

        for fl in 1:num_urbanl
            l = filter_urbanl[fl]
            g = lun_data.gridcell[l]

            # Aerodynamic resistance from urban canopy air to atmosphere
            ramu[l] = 1.0 / (ustar_loc[l]^2 / um[l])
            rahu[l] = 1.0 / (temp1[l] * ustar_loc[l])
            rawu[l] = 1.0 / (temp2[l] * ustar_loc[l])

            # Canyon wind magnitude
            canyon_wind[l] = sqrt(canyon_u_wind[l]^2.0 + ustar_loc[l]^2.0)

            # Canyon resistance (Masson 2000)
            canyon_resistance[l] = CPAIR * forc_rho_grc[g] / (11.8 + 4.2 * canyon_wind[l])
        end

        # First term in taf/qaf equations
        for fl in 1:num_urbanl
            l = filter_urbanl[fl]
            g = lun_data.gridcell[l]

            taf_numer[l] = thm_g[l] / rahu[l]
            taf_denom[l] = 1.0 / rahu[l]
            qaf_numer[l] = forc_q_grc[g] / rawu[l]
            qaf_denom[l] = 1.0 / rawu[l]

            wtas[l] = 1.0 / rahu[l]
            wtaq[l] = 1.0 / rawu[l]
        end

        # Gather terms from urban columns
        for fc in 1:num_urbanc
            c = filter_urbanc[fc]
            l = col_data.landunit[c]

            if col_data.itype[c] == ICOL_ROOF
                # Scaled sensible heat conductance
                wtus_col[c] = lun_data.wtlunit_roof[l] / canyon_resistance[l]
                wtus_roof[l] = wtus_col[c]
                wtus_roof_unscl[l] = 1.0 / canyon_resistance[l]

                # Wetness fraction for roof
                if waterdiagbulk.snow_depth_col[c] > 0.0
                    fwet_roof = min(waterdiagbulk.snow_depth_col[c] / 0.05, 1.0)
                else
                    fwet_roof = (max(0.0, waterstatebulk.ws.h2osoi_liq_col[c, 1] +
                        waterstatebulk.ws.h2osoi_ice_col[c, 1]) / PONDMX_URBAN)^0.666666666666
                    fwet_roof = min(fwet_roof, 1.0)
                end
                if waterdiagbulk.qaf_lun[l] > waterdiagbulk.qg_col[c]
                    fwet_roof = 1.0
                end

                # Scaled latent heat conductance
                wtuq_col[c] = fwet_roof * (lun_data.wtlunit_roof[l] / canyon_resistance[l])
                wtuq_roof[l] = wtuq_col[c]
                wtuq_roof_unscl[l] = fwet_roof * (1.0 / canyon_resistance[l])

                if is_simple_build_temp()
                    (wh, hfac) = simple_wasteheatfromac(
                        energyflux.eflx_urban_ac_col[c], energyflux.eflx_urban_heat_col[c])
                    eflx_wasteheat_roof[l] = wh
                    eflx_heat_from_ac_roof[l] = hfac
                end

            elseif col_data.itype[c] == ICOL_ROAD_PERV
                # Scaled sensible heat conductance
                wtus_col[c] = lun_data.wtroad_perv[l] * (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtus_road_perv[l] = wtus_col[c]
                wtus_road_perv_unscl[l] = 1.0 / canyon_resistance[l]

                # Scaled latent heat conductance
                wtuq_col[c] = lun_data.wtroad_perv[l] * (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtuq_road_perv[l] = wtuq_col[c]
                wtuq_road_perv_unscl[l] = 1.0 / canyon_resistance[l]

            elseif col_data.itype[c] == ICOL_ROAD_IMPERV
                # Scaled sensible heat conductance
                wtus_col[c] = (1.0 - lun_data.wtroad_perv[l]) * (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtus_road_imperv[l] = wtus_col[c]
                wtus_road_imperv_unscl[l] = 1.0 / canyon_resistance[l]

                # Wetness fraction for impervious road
                if waterdiagbulk.snow_depth_col[c] > 0.0
                    fwet_road_imperv = min(waterdiagbulk.snow_depth_col[c] / 0.05, 1.0)
                else
                    fwet_road_imperv = (max(0.0, waterstatebulk.ws.h2osoi_liq_col[c, 1] +
                        waterstatebulk.ws.h2osoi_ice_col[c, 1]) / PONDMX_URBAN)^0.666666666666
                    fwet_road_imperv = min(fwet_road_imperv, 1.0)
                end
                if waterdiagbulk.qaf_lun[l] > waterdiagbulk.qg_col[c]
                    fwet_road_imperv = 1.0
                end

                # Scaled latent heat conductance
                wtuq_col[c] = fwet_road_imperv * (1.0 - lun_data.wtroad_perv[l]) *
                    (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtuq_road_imperv[l] = wtuq_col[c]
                wtuq_road_imperv_unscl[l] = fwet_road_imperv * (1.0 / canyon_resistance[l])

            elseif col_data.itype[c] == ICOL_SUNWALL
                # Scaled sensible heat conductance
                wtus_col[c] = lun_data.canyon_hwr[l] * (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtus_sunwall[l] = wtus_col[c]
                wtus_sunwall_unscl[l] = 1.0 / canyon_resistance[l]

                # Walls have zero latent heat conductance
                wtuq_col[c] = 0.0
                wtuq_sunwall[l] = 0.0
                wtuq_sunwall_unscl[l] = 0.0

                if is_simple_build_temp()
                    (wh, hfac) = simple_wasteheatfromac(
                        energyflux.eflx_urban_ac_col[c], energyflux.eflx_urban_heat_col[c])
                    eflx_wasteheat_sunwall[l] = wh
                    eflx_heat_from_ac_sunwall[l] = hfac
                end

            elseif col_data.itype[c] == ICOL_SHADEWALL
                # Scaled sensible heat conductance
                wtus_col[c] = lun_data.canyon_hwr[l] * (1.0 - lun_data.wtlunit_roof[l]) / canyon_resistance[l]
                wtus_shadewall[l] = wtus_col[c]
                wtus_shadewall_unscl[l] = 1.0 / canyon_resistance[l]

                # Walls have zero latent heat conductance
                wtuq_col[c] = 0.0
                wtuq_shadewall[l] = 0.0
                wtuq_shadewall_unscl[l] = 0.0

                if is_simple_build_temp()
                    (wh, hfac) = simple_wasteheatfromac(
                        energyflux.eflx_urban_ac_col[c], energyflux.eflx_urban_heat_col[c])
                    eflx_wasteheat_shadewall[l] = wh
                    eflx_heat_from_ac_shadewall[l] = hfac
                end

            else
                error("ERROR: ctype out of range in urban_fluxes! c=$c ctype=$(col_data.itype[c])")
            end

            taf_numer[l] = taf_numer[l] + temperature.t_grnd_col[c] * wtus_col[c]
            taf_denom[l] = taf_denom[l] + wtus_col[c]
            qaf_numer[l] = qaf_numer[l] + waterdiagbulk.qg_col[c] * wtuq_col[c]
            qaf_denom[l] = qaf_denom[l] + wtuq_col[c]
        end

        # Calculate new urban canopy air temperature and specific humidity
        wasteheat!(energyflux, lun_data, num_urbanl, filter_urbanl,
            eflx_wasteheat_roof, eflx_wasteheat_sunwall, eflx_wasteheat_shadewall,
            eflx_heat_from_ac_roof, eflx_heat_from_ac_sunwall, eflx_heat_from_ac_shadewall)

        for fl in 1:num_urbanl
            l = filter_urbanl[fl]
            g = lun_data.gridcell[l]

            # Calculate traffic heat flux (only from impervious road)
            energyflux.eflx_traffic_lun[l] = (1.0 - lun_data.wtlunit_roof[l]) *
                (1.0 - lun_data.wtroad_perv[l]) * urbanparams.eflx_traffic_factor[l]

            temperature.taf_lun[l] = taf_numer[l] / taf_denom[l]
            waterdiagbulk.qaf_lun[l] = qaf_numer[l] / qaf_denom[l]

            wts_sum[l] = wtas[l] + wtus_roof[l] + wtus_road_perv[l] +
                wtus_road_imperv[l] + wtus_sunwall[l] + wtus_shadewall[l]

            wtq_sum[l] = wtaq[l] + wtuq_roof[l] + wtuq_road_perv[l] +
                wtuq_road_imperv[l] + wtuq_sunwall[l] + wtuq_shadewall[l]
        end

        # Determine stability using new taf and qaf
        for fl in 1:num_urbanl
            l = filter_urbanl[fl]
            g = lun_data.gridcell[l]

            dth[l] = thm_g[l] - temperature.taf_lun[l]
            dqh[l] = forc_q_grc[g] - waterdiagbulk.qaf_lun[l]
            tstar = temp1[l] * dth[l]
            qstar = temp2[l] * dqh[l]
            thvstar = tstar * (1.0 + 0.61 * forc_q_grc[g]) + 0.61 * forc_th_grc[g] * qstar
            zeta_lunit[l] = zldis_arr[l] * VKC * GRAV * thvstar / (ustar_loc[l]^2 * thv_g[l])

            if zeta_lunit[l] >= 0.0  # stable
                zeta_lunit[l] = min(frictionvel.zetamaxstable, max(zeta_lunit[l], 0.01))
                um[l] = max(ur[l], 0.1)
            else  # unstable
                zeta_lunit[l] = max(-100.0, min(zeta_lunit[l], -0.01))
                wc = beta_arr[l] * (-GRAV * ustar_loc[l] * thvstar * zii_arr[l] / thv_g[l])^0.333
                um[l] = sqrt(ur[l]^2 + wc^2)
            end

            obu[l] = zldis_arr[l] / zeta_lunit[l]
        end
    end  # end stability iteration

    # =========================================================================
    # Determine fluxes from canyon surfaces
    # =========================================================================

    # Initialize scaled flux arrays
    for p in begp:endp
        eflx_sh_grnd_scale[p] = 0.0
        qflx_evap_soi_scale[p] = 0.0
    end

    for f in 1:num_urbanp
        p = filter_urbanp[f]
        c = patch_data.column[p]
        g = patch_data.gridcell[p]
        l = patch_data.landunit[p]

        frictionvel.ram1_patch[p] = ramu[l]
        frictionvel.zeta_patch[p] = zeta_lunit[l]

        # Upward and downward canopy longwave are zero
        energyflux.ulrad_patch[p] = 0.0
        energyflux.dlrad_patch[p] = 0.0

        # Derivative of sensible and latent heat fluxes wrt ground temperature
        if col_data.itype[c] == ICOL_ROOF
            energyflux.cgrnds_patch[p] = forc_rho_grc[g] * CPAIR *
                (wtas[l] + wtus_road_perv[l] + wtus_road_imperv[l] + wtus_sunwall[l] + wtus_shadewall[l]) *
                (wtus_roof_unscl[l] / wts_sum[l])
            energyflux.cgrndl_patch[p] = forc_rho_grc[g] *
                (wtaq[l] + wtuq_road_perv[l] + wtuq_road_imperv[l] + wtuq_sunwall[l] + wtuq_shadewall[l]) *
                (wtuq_roof_unscl[l] / wtq_sum[l]) * waterdiagbulk.dqgdT_col[c]
        elseif col_data.itype[c] == ICOL_ROAD_PERV
            energyflux.cgrnds_patch[p] = forc_rho_grc[g] * CPAIR *
                (wtas[l] + wtus_roof[l] + wtus_road_imperv[l] + wtus_sunwall[l] + wtus_shadewall[l]) *
                (wtus_road_perv_unscl[l] / wts_sum[l])
            energyflux.cgrndl_patch[p] = forc_rho_grc[g] *
                (wtaq[l] + wtuq_roof[l] + wtuq_road_imperv[l] + wtuq_sunwall[l] + wtuq_shadewall[l]) *
                (wtuq_road_perv_unscl[l] / wtq_sum[l]) * waterdiagbulk.dqgdT_col[c]
        elseif col_data.itype[c] == ICOL_ROAD_IMPERV
            energyflux.cgrnds_patch[p] = forc_rho_grc[g] * CPAIR *
                (wtas[l] + wtus_roof[l] + wtus_road_perv[l] + wtus_sunwall[l] + wtus_shadewall[l]) *
                (wtus_road_imperv_unscl[l] / wts_sum[l])
            energyflux.cgrndl_patch[p] = forc_rho_grc[g] *
                (wtaq[l] + wtuq_roof[l] + wtuq_road_perv[l] + wtuq_sunwall[l] + wtuq_shadewall[l]) *
                (wtuq_road_imperv_unscl[l] / wtq_sum[l]) * waterdiagbulk.dqgdT_col[c]
        elseif col_data.itype[c] == ICOL_SUNWALL
            energyflux.cgrnds_patch[p] = forc_rho_grc[g] * CPAIR *
                (wtas[l] + wtus_roof[l] + wtus_road_perv[l] + wtus_road_imperv[l] + wtus_shadewall[l]) *
                (wtus_sunwall_unscl[l] / wts_sum[l])
            energyflux.cgrndl_patch[p] = 0.0
        elseif col_data.itype[c] == ICOL_SHADEWALL
            energyflux.cgrnds_patch[p] = forc_rho_grc[g] * CPAIR *
                (wtas[l] + wtus_roof[l] + wtus_road_perv[l] + wtus_road_imperv[l] + wtus_sunwall[l]) *
                (wtus_shadewall_unscl[l] / wts_sum[l])
            energyflux.cgrndl_patch[p] = 0.0
        end
        energyflux.cgrnd_patch[p] = energyflux.cgrnds_patch[p] +
            energyflux.cgrndl_patch[p] * energyflux.htvp_col[c]

        # Surface fluxes of momentum
        energyflux.taux_patch[p] = -forc_rho_grc[g] * forc_u_grc[g] / ramu[l]
        energyflux.tauy_patch[p] = -forc_rho_grc[g] * forc_v_grc[g] / ramu[l]

        # Sensible heat flux from ground using new canopy air temperature
        dth[l] = temperature.taf_lun[l] - temperature.t_grnd_col[c]

        if col_data.itype[c] == ICOL_ROOF
            energyflux.eflx_sh_grnd_patch[p] = -forc_rho_grc[g] * CPAIR * wtus_roof_unscl[l] * dth[l]
        elseif col_data.itype[c] == ICOL_ROAD_PERV
            energyflux.eflx_sh_grnd_patch[p] = -forc_rho_grc[g] * CPAIR * wtus_road_perv_unscl[l] * dth[l]
        elseif col_data.itype[c] == ICOL_ROAD_IMPERV
            energyflux.eflx_sh_grnd_patch[p] = -forc_rho_grc[g] * CPAIR * wtus_road_imperv_unscl[l] * dth[l]
        elseif col_data.itype[c] == ICOL_SUNWALL
            energyflux.eflx_sh_grnd_patch[p] = -forc_rho_grc[g] * CPAIR * wtus_sunwall_unscl[l] * dth[l]
        elseif col_data.itype[c] == ICOL_SHADEWALL
            energyflux.eflx_sh_grnd_patch[p] = -forc_rho_grc[g] * CPAIR * wtus_shadewall_unscl[l] * dth[l]
        end
        energyflux.eflx_sh_snow_patch[p] = 0.0
        energyflux.eflx_sh_soil_patch[p] = 0.0
        energyflux.eflx_sh_h2osfc_patch[p] = 0.0

        energyflux.eflx_sh_tot_patch[p] = energyflux.eflx_sh_grnd_patch[p]
        energyflux.eflx_sh_tot_u_patch[p] = energyflux.eflx_sh_tot_patch[p]

        # Latent heat flux
        dqh[l] = waterdiagbulk.qaf_lun[l] - waterdiagbulk.qg_col[c]

        if col_data.itype[c] == ICOL_ROOF
            waterfluxbulk.wf.qflx_evap_soi_patch[p] = -forc_rho_grc[g] * wtuq_roof_unscl[l] * dqh[l]
        elseif col_data.itype[c] == ICOL_ROAD_PERV
            if dqh[l] > 0.0 || waterdiagbulk.frac_sno_col[c] > 0.0 || soilstate.soilalpha_u_col[c] <= 0.0
                waterfluxbulk.wf.qflx_evap_soi_patch[p] = -forc_rho_grc[g] * wtuq_road_perv_unscl[l] * dqh[l]
                waterfluxbulk.wf.qflx_tran_veg_patch[p] = 0.0
            else
                waterfluxbulk.wf.qflx_evap_soi_patch[p] = 0.0
                waterfluxbulk.wf.qflx_tran_veg_patch[p] = -forc_rho_grc[g] * wtuq_road_perv_unscl[l] * dqh[l]
            end
            waterfluxbulk.wf.qflx_evap_veg_patch[p] = waterfluxbulk.wf.qflx_tran_veg_patch[p]
        elseif col_data.itype[c] == ICOL_ROAD_IMPERV
            waterfluxbulk.wf.qflx_evap_soi_patch[p] = -forc_rho_grc[g] * wtuq_road_imperv_unscl[l] * dqh[l]
        elseif col_data.itype[c] == ICOL_SUNWALL
            waterfluxbulk.wf.qflx_evap_soi_patch[p] = 0.0
        elseif col_data.itype[c] == ICOL_SHADEWALL
            waterfluxbulk.wf.qflx_evap_soi_patch[p] = 0.0
        end

        # SCALED sensible and latent heat flux for error check
        eflx_sh_grnd_scale[p] = -forc_rho_grc[g] * CPAIR * wtus_col[c] * dth[l]
        qflx_evap_soi_scale[p] = -forc_rho_grc[g] * wtuq_col[c] * dqh[l]
    end

    # =========================================================================
    # Error checking: total fluxes should equal sum of scaled fluxes
    # =========================================================================
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        g = lun_data.gridcell[l]
        eflx_arr[l] = -(forc_rho_grc[g] * CPAIR / rahu[l]) * (thm_g[l] - temperature.taf_lun[l])
        qflx_arr[l] = -(forc_rho_grc[g] / rawu[l]) * (forc_q_grc[g] - waterdiagbulk.qaf_lun[l])
        eflx_scale_arr[l] = sum(eflx_sh_grnd_scale[lun_data.patchi[l]:lun_data.patchf[l]])
        qflx_scale_arr[l] = sum(qflx_evap_soi_scale[lun_data.patchi[l]:lun_data.patchf[l]])
        eflx_err[l] = eflx_scale_arr[l] - eflx_arr[l]
        qflx_err[l] = qflx_scale_arr[l] - qflx_arr[l]
    end

    # Check sensible heat flux error
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        if abs(eflx_err[l]) > 0.01
            @warn "Total sensible heat does not equal sum of scaled heat fluxes for urban columns" nstep l eflx_err=eflx_err[l]
            if abs(eflx_err[l]) > 0.01
                error("urban_fluxes! sensible heat flux error > 0.01 W/m**2: " *
                      "eflx_scale=$(eflx_scale_arr[l]) eflx=$(eflx_arr[l]) err=$(eflx_err[l])")
            end
        end
    end

    # Check water vapor flux error
    for fl in 1:num_urbanl
        l = filter_urbanl[fl]
        if abs(qflx_err[l]) > 4.0e-9
            @warn "Total water vapor flux does not equal sum of scaled fluxes for urban columns" nstep l qflx_err=qflx_err[l]
            if abs(qflx_err[l]) > 4.0e-9
                error("urban_fluxes! water vapor flux error > 4e-9 kg/m**2/s: " *
                      "qflx_scale=$(qflx_scale_arr[l]) qflx=$(qflx_arr[l]) err=$(qflx_err[l])")
            end
        end
    end

    # =========================================================================
    # Internal building temperature (simple method)
    # =========================================================================
    if is_simple_build_temp()
        calc_simple_internal_building_temp!(temperature, col_data, lun_data,
            num_urbanc, filter_urbanc, num_urbanl, filter_urbanl)
    end

    # =========================================================================
    # Roots for urban (only pervious road has roots)
    # =========================================================================
    for j in 1:nlevgrnd
        for f in 1:num_urbanp
            p = filter_urbanp[f]
            c = patch_data.column[p]
            if col_data.itype[c] == ICOL_ROAD_PERV
                soilstate.rootr_patch[p, j] = soilstate.rootr_road_perv_col[c, j]
            else
                soilstate.rootr_patch[p, j] = 0.0
            end
        end
    end

    # =========================================================================
    # 2-m temperature, humidity, and diagnostics
    # =========================================================================
    for f in 1:num_urbanp
        p = filter_urbanp[f]
        c = patch_data.column[p]
        g = patch_data.gridcell[p]
        l = patch_data.landunit[p]

        # Use urban canopy air temperature and specific humidity for 2-m
        temperature.t_ref2m_patch[p] = temperature.taf_lun[l]
        waterdiagbulk.q_ref2m_patch[p] = waterdiagbulk.qaf_lun[l]
        temperature.t_ref2m_u_patch[p] = temperature.taf_lun[l]

        # 2 m height relative humidity
        (qsat_ref2m, e_ref2m, _, _) = qsat(temperature.t_ref2m_patch[p], forc_pbot_grc[g])
        waterdiagbulk.rh_ref2m_patch[p] = min(100.0,
            waterdiagbulk.q_ref2m_patch[p] / qsat_ref2m * 100.0)
        waterdiagbulk.rh_ref2m_u_patch[p] = waterdiagbulk.rh_ref2m_patch[p]

        # Human heat stress indices — stub (HumanIndexMod not yet ported)

        # Variables needed by history tape
        temperature.t_veg_patch[p] = forc_t_grc[g]
    end

    return nothing
end
