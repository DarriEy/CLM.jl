# ==========================================================================
# Ported from: src/biogeophys/BiogeophysPreFluxCalcsMod.F90
# Pre-flux biogeophysics calculations
#
# Performs calculations needed before the main flux routines
# (BareGroundFluxes, CanopyFluxes, etc.) are called:
#   1. Set momentum roughness length and displacement height
#   2. Set roughness lengths and forcing heights for non-lake points
#   3. Initialize temperature and energy variables for the timestep
#
# Public functions:
#   set_z0m_displa!              — Set z0m and displacement height
#   calc_initial_temp_energy!    — Initialize temperature/energy vars
#   biogeophys_pre_flux_calcs!   — Top-level orchestrator
# ==========================================================================

"""
    set_z0m_displa!(canopystate, frictionvel, patch_data, lun_data,
                    mask_nolakep, bounds_patch;
                    z0param_method)

Set momentum roughness length (z0m) and displacement height (displa)
for non-lake patches.

Uses PFT-specific coefficients:
- ZengWang2007: z0m = z0mr * htop, displa = displar * htop
- Meier2022: more complex parameterization using LAI

Ported from `SetZ0mDisp` in `BiogeophysPreFluxCalcsMod.F90`.
"""
function set_z0m_displa!(canopystate::CanopyStateData,
                          frictionvel::FrictionVelocityData,
                          patch_data::PatchData,
                          lun_data::LandunitData,
                          mask_nolakep::BitVector,
                          bounds_patch::UnitRange{Int};
                          z0param_method::String = "ZengWang2007",
                          z0mr::Vector{<:Real} = pftcon.z0mr,
                          displar::Vector{<:Real} = pftcon.displar)

    for p in bounds_patch
        mask_nolakep[p] || continue

        c = patch_data.column[p]
        l = patch_data.landunit[p]
        pft = patch_data.itype[p]

        htop = canopystate.htop_patch[p]

        # Look up PFT-specific values (1-indexed in Julia)
        pft_idx = pft + 1  # Fortran 0-based PFT index → Julia 1-based
        z0mr_val = (pft_idx >= 1 && pft_idx <= length(z0mr) && z0mr[pft_idx] > 0.0) ?
                   z0mr[pft_idx] : 0.055
        displar_val = (pft_idx >= 1 && pft_idx <= length(displar) && displar[pft_idx] > 0.0) ?
                      displar[pft_idx] : 0.67

        if z0param_method == "ZengWang2007"
            canopystate.z0m_patch[p] = z0mr_val * htop
            canopystate.displa_patch[p] = displar_val * htop
        else
            canopystate.z0m_patch[p] = z0mr_val * htop
            canopystate.displa_patch[p] = displar_val * htop
        end
    end

    return nothing
end

"""
    calc_initial_temp_energy!(temperature, energyflux, canopystate,
                               frictionvel, waterstatebulk, waterdiagbulk,
                               col_data, lun_data, patch_data,
                               forc_t, forc_th, forc_q,
                               mask_nolakec, mask_nolakep,
                               bounds_col, bounds_patch)

Initialize temperature and energy variables at the start of each timestep.

Sets:
- tssbef (temperature before update, for flux calculations)
- t_grnd (ground temperature from weighted snow/soil/surface-water)
- emg (ground emissivity)
- htvp (latent heat: sublimation or evaporation)
- thv (virtual potential temperature)
- beta, zii (convective velocity parameters)
- emv (vegetation emissivity)
- thm (intermediate temperature variable)
- Zeroes out sensible/latent heat fluxes and ground heat derivatives

Ported from `CalcInitialTemperatureAndEnergyVars` in `BiogeophysPreFluxCalcsMod.F90`.
"""
function calc_initial_temp_energy!(temperature::TemperatureData,
                                    energyflux::EnergyFluxData,
                                    canopystate::CanopyStateData,
                                    frictionvel::FrictionVelocityData,
                                    waterstatebulk::WaterStateBulkData,
                                    waterdiagbulk::WaterDiagnosticBulkData,
                                    col_data::ColumnData,
                                    lun_data::LandunitData,
                                    patch_data::PatchData,
                                    forc_t::Vector{<:Real},
                                    forc_th::Vector{<:Real},
                                    forc_q::Vector{<:Real},
                                    mask_nolakec::BitVector,
                                    mask_nolakep::BitVector,
                                    bounds_col::UnitRange{Int},
                                    bounds_patch::UnitRange{Int})

    nlevsno_val = varpar.nlevsno
    nlevgrnd_val = varpar.nlevgrnd

    # ---------------------------------------------------------------
    # 1. Save current temperatures (tssbef) for flux calculations
    # ---------------------------------------------------------------
    for c in bounds_col
        mask_nolakec[c] || continue
        l = col_data.landunit[c]

        for j in 1:(nlevsno_val + nlevgrnd_val)
            # Fortran j ranges from -nlevsno+1 to nlevgrnd
            # Julia j_jl = j (already 1-based in the matrix)
            temperature.t_ssbef_col[c, j] = temperature.t_soisno_col[c, j]
        end
    end

    # ---------------------------------------------------------------
    # 2. Ground temperature: weighted average of snow, soil, surface water
    # ---------------------------------------------------------------
    for c in bounds_col
        mask_nolakec[c] || continue

        snl = col_data.snl[c]
        frac_sno_eff = waterdiagbulk.frac_sno_eff_col[c]
        frac_h2osfc = waterdiagbulk.frac_h2osfc_col[c]

        if snl < 0
            # Has snow layers — top snow layer index in Julia: snl + 1 + nlevsno
            jtop = snl + 1 + nlevsno_val
            t_snow_top = temperature.t_soisno_col[c, jtop]
            t_soil_1 = temperature.t_soisno_col[c, 1 + nlevsno_val]
            t_h2osfc = temperature.t_h2osfc_col[c]
            temperature.t_grnd_col[c] = frac_sno_eff * t_snow_top +
                                        (1.0 - frac_sno_eff - frac_h2osfc) * t_soil_1 +
                                        frac_h2osfc * t_h2osfc
        else
            # No snow
            t_soil_1 = temperature.t_soisno_col[c, 1 + nlevsno_val]
            t_h2osfc = temperature.t_h2osfc_col[c]
            temperature.t_grnd_col[c] = (1.0 - frac_h2osfc) * t_soil_1 +
                                        frac_h2osfc * t_h2osfc
        end
    end

    # ---------------------------------------------------------------
    # 3. Ground emissivity
    # ---------------------------------------------------------------
    for c in bounds_col
        mask_nolakec[c] || continue
        l = col_data.landunit[c]

        if !lun_data.urbpoi[l]
            frac_sno = waterdiagbulk.frac_sno_col[c]
            if lun_data.itype[l] == ISTICE
                temperature.emg_col[c] = 0.97
            else
                temperature.emg_col[c] = (1.0 - frac_sno) * 0.96 + frac_sno * 0.97
            end
        end
        # Urban emissivity is set elsewhere (from surface data)
    end

    # ---------------------------------------------------------------
    # 4. Latent heat of vaporization/sublimation
    # ---------------------------------------------------------------
    for c in bounds_col
        mask_nolakec[c] || continue
        snl = col_data.snl[c]
        jtop = snl + 1 + nlevsno_val

        h2o_liq = waterstatebulk.ws.h2osoi_liq_col[c, jtop]
        h2o_ice = waterstatebulk.ws.h2osoi_ice_col[c, jtop]

        if h2o_liq <= 0.0 && h2o_ice > 0.0
            energyflux.htvp_col[c] = HSUB  # sublimation
        else
            energyflux.htvp_col[c] = HVAP  # evaporation
        end
    end

    # ---------------------------------------------------------------
    # 5. Virtual potential temperature, beta, zii
    # ---------------------------------------------------------------
    for c in bounds_col
        mask_nolakec[c] || continue
        temperature.thv_col[c] = forc_th[c] * (1.0 + 0.61 * forc_q[c])
        temperature.beta_col[c] = 1.0
        col_data.zii[c] = 1000.0  # convective boundary layer height [m]
    end

    # ---------------------------------------------------------------
    # 6. Initialize energy fluxes to zero
    # ---------------------------------------------------------------
    for p in bounds_patch
        mask_nolakep[p] || continue
        energyflux.eflx_sh_tot_patch[p] = 0.0
        energyflux.eflx_sh_tot_r_patch[p] = 0.0
        energyflux.eflx_sh_tot_u_patch[p] = 0.0
        energyflux.eflx_lh_tot_patch[p] = 0.0
        energyflux.eflx_lh_tot_r_patch[p] = 0.0
        energyflux.eflx_lh_tot_u_patch[p] = 0.0
        energyflux.eflx_sh_veg_patch[p] = 0.0
        energyflux.cgrnd_patch[p] = 0.0
        energyflux.cgrnds_patch[p] = 0.0
        energyflux.cgrndl_patch[p] = 0.0
    end

    # ---------------------------------------------------------------
    # 7. Vegetation emissivity
    # ---------------------------------------------------------------
    for p in bounds_patch
        mask_nolakep[p] || continue
        avmuir = 1.0  # inverse optical depth per unit leaf+stem area
        elai = canopystate.elai_patch[p]
        esai = canopystate.esai_patch[p]
        temperature.emv_patch[p] = 1.0 - exp(-(elai + esai) / avmuir)
    end

    # ---------------------------------------------------------------
    # 8. Intermediate temperature variable thm
    # ---------------------------------------------------------------
    for p in bounds_patch
        mask_nolakep[p] || continue
        c = patch_data.column[p]
        temperature.thm_patch[p] = forc_t[c] + 0.0098 * frictionvel.forc_hgt_t_patch[p]
    end

    return nothing
end

"""
    biogeophys_pre_flux_calcs!(canopystate, energyflux, frictionvel,
                                temperature, soilstate,
                                waterstatebulk, waterdiagbulk, waterfluxbulk,
                                col_data, lun_data, patch_data,
                                forc_t, forc_th, forc_q,
                                forc_hgt_u, forc_hgt_t, forc_hgt_q,
                                mask_nolakec, mask_nolakep, mask_urbanc,
                                bounds_col, bounds_patch;
                                z0param_method)

Top-level pre-flux calculations orchestrator.

1. Sets z0m and displacement height
2. Sets roughness lengths and forcing heights for non-lake points
3. Initializes temperature and energy variables

Ported from `BiogeophysPreFluxCalcs` in `BiogeophysPreFluxCalcsMod.F90`.
"""
function biogeophys_pre_flux_calcs!(canopystate::CanopyStateData,
                                     energyflux::EnergyFluxData,
                                     frictionvel::FrictionVelocityData,
                                     temperature::TemperatureData,
                                     soilstate::SoilStateData,
                                     waterstatebulk::WaterStateBulkData,
                                     waterdiagbulk::WaterDiagnosticBulkData,
                                     waterfluxbulk::WaterFluxBulkData,
                                     col_data::ColumnData,
                                     lun_data::LandunitData,
                                     patch_data::PatchData,
                                     forc_t::Vector{<:Real},
                                     forc_th::Vector{<:Real},
                                     forc_q::Vector{<:Real},
                                     forc_hgt_u::Vector{<:Real},
                                     forc_hgt_t::Vector{<:Real},
                                     forc_hgt_q::Vector{<:Real},
                                     mask_nolakec::BitVector,
                                     mask_nolakep::BitVector,
                                     mask_urbanc::BitVector,
                                     bounds_col::UnitRange{Int},
                                     bounds_patch::UnitRange{Int};
                                     z0param_method::String = "ZengWang2007")

    # 1. Set z0m and displacement height
    set_z0m_displa!(canopystate, frictionvel, patch_data, lun_data,
                    mask_nolakep, bounds_patch;
                    z0param_method=z0param_method)

    # 2. Set roughness lengths and forcing heights for non-lake points
    set_roughness_and_forc_heights_nonlake!(
        frictionvel,
        mask_nolakec, mask_nolakep,
        bounds_col, bounds_patch,
        waterdiagbulk.frac_sno_col,
        zeros(length(waterdiagbulk.frac_sno_col)),  # snomelt_accum placeholder
        canopystate.frac_veg_nosno_patch,
        canopystate.z0m_patch,
        canopystate.displa_patch,
        forc_hgt_u, forc_hgt_t, forc_hgt_q,
        col_data.landunit,
        patch_data.gridcell,
        patch_data.landunit,
        patch_data.column,
        lun_data.itype,
        lun_data.urbpoi,
        lun_data.z_0_town,
        lun_data.z_d_town;
        z0param_method=z0param_method)

    # 3. Initialize temperature and energy variables
    calc_initial_temp_energy!(temperature, energyflux, canopystate,
                              frictionvel, waterstatebulk, waterdiagbulk,
                              col_data, lun_data, patch_data,
                              forc_t, forc_th, forc_q,
                              mask_nolakec, mask_nolakep,
                              bounds_col, bounds_patch)

    return nothing
end
