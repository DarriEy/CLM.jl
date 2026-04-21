# ==========================================================================
# SoilFluxes — ported from SoilFluxesMod.F90
#
# Update surface fluxes based on the new ground temperature.
# ==========================================================================

"""
    soil_fluxes!(...)

Update surface fluxes based on the new ground temperature.

Corrects sensible heat, latent heat, and soil heat fluxes for the change
in ground temperature between the previous and current time steps.
Partitions evaporation into liquid and solid components, constrains snow
evaporation to available moisture, computes total fluxes, outgoing
longwave radiation, and the soil energy balance error.

Ported from: SoilFluxesMod.F90 :: SoilFluxes
"""
function soil_fluxes!(
        # Data structures
        energyflux       ::EnergyFluxData,
        temperature      ::TemperatureData,
        canopystate      ::CanopyStateData,
        waterstatebulk   ::WaterStateBulkData,
        waterdiagbulk    ::WaterDiagnosticBulkData,
        waterfluxbulk    ::WaterFluxBulkData,
        solarabs         ::SolarAbsorbedData,
        patch_data       ::PatchData,
        col_data         ::ColumnData,
        lun_data         ::LandunitData,
        # Masks
        mask_nolakec     ::BitVector,
        mask_nolakep     ::BitVector,
        mask_urbanp      ::BitVector,
        # Bounds
        bounds_col       ::UnitRange{Int},
        bounds_patch     ::UnitRange{Int},
        # Atmospheric forcing (column-level)
        forc_lwrad_col   ::Vector{<:Real},
        # Time step
        dtime            ::Real)

    nlevsno  = varpar.nlevsno
    nlevgrnd = varpar.nlevgrnd
    nlevurb  = varpar.nlevurb

    # Local work arrays
    endc = last(bounds_col)
    endp = last(bounds_patch)
    FT = eltype(temperature.t_grnd_col)
    tinc    = zeros(FT, endc)
    t_grnd0 = zeros(FT, endc)
    eflx_lwrad_del = zeros(FT, endp)

    # =========================================================================
    # Loop 1: Calculate temperature difference for flux corrections (column)
    # =========================================================================

    for c in bounds_col
        mask_nolakec[c] || continue

        if col_data.snl[c] < 0
            t_grnd0[c] = waterdiagbulk.frac_sno_eff_col[c] *
                    temperature.t_ssbef_col[c, col_data.snl[c] + 1 + nlevsno] +
                (1.0 - waterdiagbulk.frac_sno_eff_col[c] - waterdiagbulk.frac_h2osfc_col[c]) *
                    temperature.t_ssbef_col[c, 1 + nlevsno] +
                waterdiagbulk.frac_h2osfc_col[c] * temperature.t_h2osfc_bef_col[c]
        else
            t_grnd0[c] = (1.0 - waterdiagbulk.frac_h2osfc_col[c]) *
                    temperature.t_ssbef_col[c, 1 + nlevsno] +
                waterdiagbulk.frac_h2osfc_col[c] * temperature.t_h2osfc_bef_col[c]
        end

        tinc[c] = temperature.t_grnd_col[c] - t_grnd0[c]
    end

    # =========================================================================
    # Loop 2: Correct fluxes to present soil temperature (patch)
    # =========================================================================

    for p in bounds_patch
        mask_nolakep[p] || continue
        c = patch_data.column[p]

        energyflux.eflx_sh_grnd_patch[p] += tinc[c] * energyflux.cgrnds_patch[p]
        waterfluxbulk.wf.qflx_evap_soi_patch[p] += tinc[c] * energyflux.cgrndl_patch[p]

        # Set ev_soil, ev_h2osfc, ev_snow for urban landunits
        l = patch_data.landunit[p]
        if lun_data.urbpoi[l]
            waterfluxbulk.qflx_ev_soil_patch[p] = 0.0
            waterfluxbulk.qflx_ev_h2osfc_patch[p] = 0.0
            waterfluxbulk.qflx_ev_snow_patch[p] = waterfluxbulk.wf.qflx_evap_soi_patch[p]
        else
            waterfluxbulk.qflx_ev_snow_patch[p] += tinc[c] * energyflux.cgrndl_patch[p]
            waterfluxbulk.qflx_ev_soil_patch[p] += tinc[c] * energyflux.cgrndl_patch[p]
            waterfluxbulk.qflx_ev_h2osfc_patch[p] += tinc[c] * energyflux.cgrndl_patch[p]
        end
    end

    # =========================================================================
    # Loop 3: Partition evaporation into liquid and solid (patch)
    # =========================================================================

    for p in bounds_patch
        mask_nolakep[p] || continue
        c = patch_data.column[p]
        l = patch_data.landunit[p]
        j = col_data.snl[c] + 1       # Fortran-style top layer index
        j_jl = j + nlevsno            # Julia array index

        waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[p] = 0.0
        waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[p] = 0.0
        waterfluxbulk.wf.qflx_soliddew_to_top_layer_patch[p] = 0.0
        waterfluxbulk.wf.qflx_liqdew_to_top_layer_patch[p] = 0.0

        h2o_liq = waterstatebulk.ws.h2osoi_liq_col[c, j_jl]
        h2o_ice = waterstatebulk.ws.h2osoi_ice_col[c, j_jl]

        if !lun_data.urbpoi[l]
            if waterfluxbulk.qflx_ev_snow_patch[p] >= 0.0
                # Evaporation: partition between liquid evap and ice sublimation
                if (h2o_liq + h2o_ice) > 0.0
                    waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[p] = max(
                        waterfluxbulk.qflx_ev_snow_patch[p] * (h2o_liq / (h2o_liq + h2o_ice)), 0.0)
                else
                    waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[p] = 0.0
                end
                waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[p] =
                    waterfluxbulk.qflx_ev_snow_patch[p] -
                    waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[p]
            else
                # Dew/frost deposition
                if temperature.t_grnd_col[c] < TFRZ
                    waterfluxbulk.wf.qflx_soliddew_to_top_layer_patch[p] = abs(waterfluxbulk.qflx_ev_snow_patch[p])
                else
                    waterfluxbulk.wf.qflx_liqdew_to_top_layer_patch[p] = abs(waterfluxbulk.qflx_ev_snow_patch[p])
                end
            end

        else  # Urban columns
            if waterfluxbulk.wf.qflx_evap_soi_patch[p] >= 0.0
                if (h2o_liq + h2o_ice) > 0.0
                    waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[p] = max(
                        waterfluxbulk.wf.qflx_evap_soi_patch[p] * (h2o_liq / (h2o_liq + h2o_ice)), 0.0)
                else
                    waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[p] = 0.0
                end
                waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[p] =
                    waterfluxbulk.wf.qflx_evap_soi_patch[p] -
                    waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[p]
            else
                if temperature.t_grnd_col[c] < TFRZ
                    waterfluxbulk.wf.qflx_soliddew_to_top_layer_patch[p] = abs(waterfluxbulk.wf.qflx_evap_soi_patch[p])
                else
                    waterfluxbulk.wf.qflx_liqdew_to_top_layer_patch[p] = abs(waterfluxbulk.wf.qflx_evap_soi_patch[p])
                end
            end
        end
    end

    # =========================================================================
    # Loop 4: Constrain evaporation from snow to be <= available moisture
    # =========================================================================

    for p in bounds_patch
        mask_nolakep[p] || continue
        c = patch_data.column[p]
        j = col_data.snl[c] + 1
        j_jl = j + nlevsno

        h2o_liq = waterstatebulk.ws.h2osoi_liq_col[c, j_jl]
        h2o_ice = waterstatebulk.ws.h2osoi_ice_col[c, j_jl]

        # Snow layers; assumes for j < 1 that frac_sno_eff > 0
        if j < 1
            evaporation_limit = (h2o_ice + h2o_liq) / (waterdiagbulk.frac_sno_eff_col[c] * dtime)
            if waterfluxbulk.qflx_ev_snow_patch[p] > evaporation_limit
                evaporation_demand = waterfluxbulk.qflx_ev_snow_patch[p]
                waterfluxbulk.qflx_ev_snow_patch[p] = evaporation_limit
                waterfluxbulk.wf.qflx_evap_soi_patch[p] -=
                    waterdiagbulk.frac_sno_eff_col[c] * (evaporation_demand - evaporation_limit)
                waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[p] =
                    max(h2o_liq / (waterdiagbulk.frac_sno_eff_col[c] * dtime), 0.0)
                waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[p] =
                    max(h2o_ice / (waterdiagbulk.frac_sno_eff_col[c] * dtime), 0.0)
                # Conserve total energy flux
                energyflux.eflx_sh_grnd_patch[p] +=
                    waterdiagbulk.frac_sno_eff_col[c] * (evaporation_demand - evaporation_limit) * energyflux.htvp_col[c]
            end
        end

        # Top soil layer for urban columns (excluding pervious road)
        if lun_data.urbpoi[patch_data.landunit[p]] && (col_data.itype[c] != ICOL_ROAD_PERV) && (j == 1)
            evaporation_limit = (h2o_ice + h2o_liq) / dtime
            if waterfluxbulk.wf.qflx_evap_soi_patch[p] > evaporation_limit
                evaporation_demand = waterfluxbulk.wf.qflx_evap_soi_patch[p]
                waterfluxbulk.wf.qflx_evap_soi_patch[p] = evaporation_limit
                waterfluxbulk.qflx_ev_snow_patch[p] = waterfluxbulk.wf.qflx_evap_soi_patch[p]
                waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[p] = max(h2o_liq / dtime, 0.0)
                waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[p] = max(h2o_ice / dtime, 0.0)
                # Conserve total energy flux
                energyflux.eflx_sh_grnd_patch[p] +=
                    (evaporation_demand - evaporation_limit) * energyflux.htvp_col[c]
            end
        end

        # Limit only solid evaporation (sublimation) from top soil layer
        if j == 1 && waterdiagbulk.frac_h2osfc_col[c] < 1.0
            evaporation_limit = h2o_ice / (dtime * (1.0 - waterdiagbulk.frac_h2osfc_col[c]))
            if waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[p] >= evaporation_limit
                evaporation_demand = waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[p]
                waterfluxbulk.wf.qflx_solidevap_from_top_layer_patch[p] = evaporation_limit
                waterfluxbulk.wf.qflx_liqevap_from_top_layer_patch[p] +=
                    (evaporation_demand - evaporation_limit)
            end
        end
    end

    # =========================================================================
    # Loop 5: Ground heat flux, total fluxes, history variables (patch)
    # =========================================================================

    for p in bounds_patch
        mask_nolakep[p] || continue
        c = patch_data.column[p]
        l = patch_data.landunit[p]

        # Ground heat flux
        if !lun_data.urbpoi[l]
            lw_grnd = waterdiagbulk.frac_sno_eff_col[c] *
                    temperature.t_ssbef_col[c, col_data.snl[c] + 1 + nlevsno]^4 +
                (1.0 - waterdiagbulk.frac_sno_eff_col[c] - waterdiagbulk.frac_h2osfc_col[c]) *
                    temperature.t_ssbef_col[c, 1 + nlevsno]^4 +
                waterdiagbulk.frac_h2osfc_col[c] * temperature.t_h2osfc_bef_col[c]^4

            energyflux.eflx_soil_grnd_patch[p] =
                ((1.0 - waterdiagbulk.frac_sno_eff_col[c]) * solarabs.sabg_soil_patch[p] +
                 waterdiagbulk.frac_sno_eff_col[c] * solarabs.sabg_snow_patch[p]) +
                energyflux.dlrad_patch[p] +
                (1 - canopystate.frac_veg_nosno_patch[p]) * temperature.emg_col[c] * forc_lwrad_col[c] -
                temperature.emg_col[c] * SB * lw_grnd -
                temperature.emg_col[c] * SB * t_grnd0[c]^3 * (4.0 * tinc[c]) -
                (energyflux.eflx_sh_grnd_patch[p] + waterfluxbulk.wf.qflx_evap_soi_patch[p] * energyflux.htvp_col[c])

            if lun_data.itype[l] == ISTSOIL || lun_data.itype[l] == ISTCROP
                energyflux.eflx_soil_grnd_r_patch[p] = energyflux.eflx_soil_grnd_patch[p]
            end
        else
            # Urban columns
            eflx_lwrad_del[p] = 4.0 * temperature.emg_col[c] * SB * t_grnd0[c]^3 * tinc[c]

            # Include transpiration term because needed for pervious road
            # and wasteheat and traffic flux
            energyflux.eflx_soil_grnd_patch[p] = solarabs.sabg_patch[p] + energyflux.dlrad_patch[p] -
                energyflux.eflx_lwrad_net_patch[p] - eflx_lwrad_del[p] -
                (energyflux.eflx_sh_grnd_patch[p] +
                 waterfluxbulk.wf.qflx_evap_soi_patch[p] * energyflux.htvp_col[c] +
                 waterfluxbulk.wf.qflx_tran_veg_patch[p] * HVAP) +
                energyflux.eflx_wasteheat_patch[p] + energyflux.eflx_heat_from_ac_patch[p] +
                energyflux.eflx_traffic_patch[p] + energyflux.eflx_ventilation_patch[p]
            energyflux.eflx_soil_grnd_u_patch[p] = energyflux.eflx_soil_grnd_patch[p]
        end

        # Total fluxes (vegetation + ground)
        energyflux.eflx_sh_tot_patch[p] = energyflux.eflx_sh_veg_patch[p] + energyflux.eflx_sh_grnd_patch[p]
        if !lun_data.urbpoi[l]
            energyflux.eflx_sh_tot_patch[p] += energyflux.eflx_sh_stem_patch[p]
        end
        waterfluxbulk.wf.qflx_evap_tot_patch[p] =
            waterfluxbulk.wf.qflx_evap_veg_patch[p] + waterfluxbulk.wf.qflx_evap_soi_patch[p]

        energyflux.eflx_lh_tot_patch[p] =
            HVAP * waterfluxbulk.wf.qflx_evap_veg_patch[p] +
            energyflux.htvp_col[c] * waterfluxbulk.wf.qflx_evap_soi_patch[p]
        if lun_data.itype[l] == ISTSOIL || lun_data.itype[l] == ISTCROP
            energyflux.eflx_lh_tot_r_patch[p] = energyflux.eflx_lh_tot_patch[p]
            energyflux.eflx_sh_tot_r_patch[p] = energyflux.eflx_sh_tot_patch[p]
        elseif lun_data.urbpoi[l]
            energyflux.eflx_lh_tot_u_patch[p] = energyflux.eflx_lh_tot_patch[p]
            energyflux.eflx_sh_tot_u_patch[p] = energyflux.eflx_sh_tot_patch[p]
        end

        # Variables needed by history tape
        waterfluxbulk.wf.qflx_evap_can_patch[p] =
            waterfluxbulk.wf.qflx_evap_veg_patch[p] - waterfluxbulk.wf.qflx_tran_veg_patch[p]
        energyflux.eflx_lh_vege_patch[p] =
            (waterfluxbulk.wf.qflx_evap_veg_patch[p] - waterfluxbulk.wf.qflx_tran_veg_patch[p]) * HVAP
        energyflux.eflx_lh_vegt_patch[p] = waterfluxbulk.wf.qflx_tran_veg_patch[p] * HVAP
        energyflux.eflx_lh_grnd_patch[p] = waterfluxbulk.wf.qflx_evap_soi_patch[p] * energyflux.htvp_col[c]
    end

    # =========================================================================
    # Loop 6: Soil energy balance check (errsoi_patch)
    # =========================================================================

    for p in bounds_patch
        mask_nolakep[p] || continue
        c = patch_data.column[p]

        energyflux.errsoi_patch[p] = energyflux.eflx_soil_grnd_patch[p] -
            temperature.xmf_col[c] - temperature.xmf_h2osfc_col[c] -
            waterdiagbulk.frac_h2osfc_col[c] *
            (temperature.t_h2osfc_col[c] - temperature.t_h2osfc_bef_col[c]) *
            (temperature.c_h2osfc_col[c] / dtime)
        energyflux.errsoi_patch[p] += energyflux.eflx_h2osfc_to_snow_col[c]

        # For urban sunwall, shadewall, and roof columns, include building heat flux
        if col_data.itype[c] == ICOL_SUNWALL || col_data.itype[c] == ICOL_SHADEWALL || col_data.itype[c] == ICOL_ROOF
            energyflux.errsoi_patch[p] += energyflux.eflx_building_heat_errsoi_col[c]
        end
    end

    # Nested loop: errsoi contributions for non-urban columns (j from -nlevsno+1 to nlevgrnd)
    for j_f in (-nlevsno + 1):nlevgrnd
        j_jl = j_f + nlevsno  # Julia array index
        for p in bounds_patch
            mask_nolakep[p] || continue
            c = patch_data.column[p]

            if col_data.itype[c] != ICOL_SUNWALL && col_data.itype[c] != ICOL_SHADEWALL &&
               col_data.itype[c] != ICOL_ROOF
                # Area weight heat absorbed by snow layers
                if j_f >= col_data.snl[c] + 1 && j_f < 1
                    energyflux.errsoi_patch[p] -= waterdiagbulk.frac_sno_eff_col[c] *
                        (temperature.t_soisno_col[c, j_jl] - temperature.t_ssbef_col[c, j_jl]) /
                        temperature.fact_col[c, j_jl]
                end
                if j_f >= 1
                    energyflux.errsoi_patch[p] -=
                        (temperature.t_soisno_col[c, j_jl] - temperature.t_ssbef_col[c, j_jl]) /
                        temperature.fact_col[c, j_jl]
                end
            end
        end
    end

    # Nested loop: errsoi contributions for urban columns (j from -nlevsno+1 to nlevurb)
    for j_f in (-nlevsno + 1):nlevurb
        j_jl = j_f + nlevsno
        for p in bounds_patch
            mask_urbanp[p] || continue
            c = patch_data.column[p]

            if col_data.itype[c] == ICOL_SUNWALL || col_data.itype[c] == ICOL_SHADEWALL ||
               col_data.itype[c] == ICOL_ROOF
                # Area weight heat absorbed by snow layers
                if j_f >= col_data.snl[c] + 1 && j_f < 1
                    energyflux.errsoi_patch[p] -= waterdiagbulk.frac_sno_eff_col[c] *
                        (temperature.t_soisno_col[c, j_jl] - temperature.t_ssbef_col[c, j_jl]) /
                        temperature.fact_col[c, j_jl]
                end
                if j_f >= 1
                    energyflux.errsoi_patch[p] -=
                        (temperature.t_soisno_col[c, j_jl] - temperature.t_ssbef_col[c, j_jl]) /
                        temperature.fact_col[c, j_jl]
                end
            end
        end
    end

    # =========================================================================
    # Loop 7: Outgoing longwave radiation from vegetation + ground
    # =========================================================================

    for p in bounds_patch
        mask_nolakep[p] || continue
        c = patch_data.column[p]
        l = patch_data.landunit[p]

        if !lun_data.urbpoi[l]
            lw_grnd = waterdiagbulk.frac_sno_eff_col[c] *
                    temperature.t_ssbef_col[c, col_data.snl[c] + 1 + nlevsno]^4 +
                (1.0 - waterdiagbulk.frac_sno_eff_col[c] - waterdiagbulk.frac_h2osfc_col[c]) *
                    temperature.t_ssbef_col[c, 1 + nlevsno]^4 +
                waterdiagbulk.frac_h2osfc_col[c] * temperature.t_h2osfc_bef_col[c]^4

            energyflux.eflx_lwrad_out_patch[p] = energyflux.ulrad_patch[p] +
                (1 - canopystate.frac_veg_nosno_patch[p]) * (1.0 - temperature.emg_col[c]) * forc_lwrad_col[c] +
                (1 - canopystate.frac_veg_nosno_patch[p]) * temperature.emg_col[c] * SB * lw_grnd +
                4.0 * temperature.emg_col[c] * SB * t_grnd0[c]^3 * tinc[c]

            # Bare ground skin temperature
            if canopystate.frac_veg_nosno_patch[p] == 0
                temperature.t_skin_patch[p] = sqrt(sqrt(lw_grnd))
            end

            energyflux.eflx_lwrad_net_patch[p] = energyflux.eflx_lwrad_out_patch[p] - forc_lwrad_col[c]
            if lun_data.itype[l] == ISTSOIL || lun_data.itype[l] == ISTCROP
                energyflux.eflx_lwrad_net_r_patch[p] = energyflux.eflx_lwrad_out_patch[p] - forc_lwrad_col[c]
                energyflux.eflx_lwrad_out_r_patch[p] = energyflux.eflx_lwrad_out_patch[p]
            end
        else
            energyflux.eflx_lwrad_out_patch[p] += eflx_lwrad_del[p]
            energyflux.eflx_lwrad_net_patch[p] += eflx_lwrad_del[p]
            energyflux.eflx_lwrad_net_u_patch[p] += eflx_lwrad_del[p]
            energyflux.eflx_lwrad_out_u_patch[p] = energyflux.eflx_lwrad_out_patch[p]
        end
    end

    # =========================================================================
    # Patch-to-column averaging for errsoi (replaces Fortran p2c call)
    # =========================================================================

    for c in bounds_col
        mask_nolakec[c] || continue
        energyflux.errsoi_col[c] = 0.0
        for p in col_data.patchi[c]:col_data.patchf[c]
            energyflux.errsoi_col[c] += energyflux.errsoi_patch[p] * patch_data.wtcol[p]
        end
    end

    # =========================================================================
    # Assign column-level t_soisno(snl+1) to t_skin for urban patches
    # =========================================================================

    for p in bounds_patch
        mask_urbanp[p] || continue
        c = patch_data.column[p]
        temperature.t_skin_patch[p] = temperature.t_soisno_col[c, col_data.snl[c] + 1 + nlevsno]
    end

    return nothing
end
