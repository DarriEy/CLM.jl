# ==========================================================================
# KernelAbstractions kernels for pre-flux per-column/per-patch setup.
# Each kernel iterates one independent element (column or patch). The masked
# Fortran filter loop becomes an in-kernel `if mask[i]` guard; globals/params
# are passed as @Const/scalar args. Backend taken from the (written) output
# array via the module-wide _launch! helper. See src/infrastructure/kernels.jl.
# ==========================================================================

# --------------------------------------------------------------------------
# z0m / displa for non-lake patches. PFT-specific coefficients are looked up
# from z0mr/displar with a fallback (matches the scalar code; both z0param
# branches were identical so the common expression is used).
# --------------------------------------------------------------------------
@kernel function _preflux_z0m_displa_kernel!(z0m_patch, displa_patch,
                                             @Const(mask), @Const(itype),
                                             @Const(htop_patch), @Const(z0mr),
                                             @Const(displar))
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(z0m_patch)
        htop = htop_patch[p]
        pft_idx = itype[p] + 1  # Fortran 0-based PFT index -> Julia 1-based
        z0mr_val = (pft_idx >= 1 && pft_idx <= length(z0mr) && z0mr[pft_idx] > zero(T)) ?
                   z0mr[pft_idx] : T(0.055)
        displar_val = (pft_idx >= 1 && pft_idx <= length(displar) && displar[pft_idx] > zero(T)) ?
                      displar[pft_idx] : T(0.67)
        z0m_patch[p] = z0mr_val * htop
        displa_patch[p] = displar_val * htop
    end
end

"""
    preflux_z0m_displa!(z0m_patch, displa_patch, mask, itype, htop_patch, z0mr, displar)

Set momentum roughness length and displacement height for non-lake patches,
one thread per patch. Backend-agnostic.
"""
preflux_z0m_displa!(z0m_patch, displa_patch, mask, itype, htop_patch, z0mr, displar) =
    _launch!(_preflux_z0m_displa_kernel!, z0m_patch, displa_patch, mask, itype,
             htop_patch, z0mr, displar)

# --------------------------------------------------------------------------
# Save current temperatures (tssbef) for flux calculations — 2D over
# (column, level), pure copy of t_soisno into t_ssbef on active columns.
# --------------------------------------------------------------------------
@kernel function _preflux_tssbef_kernel!(t_ssbef_col, @Const(mask), @Const(t_soisno_col))
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        t_ssbef_col[c, j] = t_soisno_col[c, j]
    end
end

"""
    preflux_tssbef!(t_ssbef_col, mask, t_soisno_col, nlev)

Copy t_soisno into t_ssbef over all (column, level) pairs for active columns.
2D backend-agnostic kernel.
"""
preflux_tssbef!(t_ssbef_col, mask, t_soisno_col, nlev::Int) =
    _launch!(_preflux_tssbef_kernel!, t_ssbef_col, mask, t_soisno_col;
             ndrange = (length(mask), nlev))

# --------------------------------------------------------------------------
# Ground temperature: weighted average of snow / soil / surface water.
# --------------------------------------------------------------------------
@kernel function _preflux_t_grnd_kernel!(t_grnd_col, t_h2osfc_bef_col, @Const(mask), @Const(snl),
                                         @Const(frac_sno_eff_col), @Const(frac_h2osfc_col),
                                         @Const(t_soisno_col), @Const(t_h2osfc_col),
                                         nlevsno::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(t_grnd_col)
        snl_c = snl[c]
        frac_sno_eff = frac_sno_eff_col[c]
        frac_h2osfc = frac_h2osfc_col[c]
        t_soil_1 = t_soisno_col[c, 1 + nlevsno]
        t_h2osfc = t_h2osfc_col[c]
        # Record surface-water temperature prior to the soil-temperature update
        # (Fortran BiogeophysPreFluxCalcsMod.F90:331 `t_h2osfc_bef(c) = t_h2osfc(c)`).
        # SoilFluxes reads t_h2osfc_bef in lw_grnd / t_grnd0 / the h2osfc heat-storage
        # term; without this per-step snapshot it stayed pinned at the cold-start value
        # (274 K), so on ponded forest columns (frac_h2osfc ~ 0.5) the emitted-longwave
        # assembly used a stale-cold surface-water temperature → a season-long ~35 W/m²
        # surface-energy-balance leak (FIRA/FSH/FGR inflated) once the canopy leafed out.
        t_h2osfc_bef_col[c] = t_h2osfc
        if snl_c < 0
            jtop = snl_c + 1 + nlevsno
            t_snow_top = t_soisno_col[c, jtop]
            t_grnd_col[c] = frac_sno_eff * t_snow_top +
                            (one(T) - frac_sno_eff - frac_h2osfc) * t_soil_1 +
                            frac_h2osfc * t_h2osfc
        else
            t_grnd_col[c] = (one(T) - frac_h2osfc) * t_soil_1 +
                            frac_h2osfc * t_h2osfc
        end
    end
end

"""
    preflux_t_grnd!(t_grnd_col, mask, snl, frac_sno_eff_col, frac_h2osfc_col,
                    t_soisno_col, t_h2osfc_col, nlevsno)

Ground temperature as snow/soil/surface-water weighted average, per column.
"""
preflux_t_grnd!(t_grnd_col, t_h2osfc_bef_col, mask, snl, frac_sno_eff_col, frac_h2osfc_col,
                t_soisno_col, t_h2osfc_col, nlevsno::Int) =
    _launch!(_preflux_t_grnd_kernel!, t_grnd_col, t_h2osfc_bef_col, mask, snl, frac_sno_eff_col,
             frac_h2osfc_col, t_soisno_col, t_h2osfc_col, nlevsno)

# --------------------------------------------------------------------------
# Ground emissivity (non-urban columns; urban set elsewhere).
# --------------------------------------------------------------------------
@kernel function _preflux_emg_kernel!(emg_col, @Const(mask), @Const(landunit),
                                      @Const(urbpoi), @Const(lun_itype),
                                      @Const(frac_sno_col), istice::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(emg_col)
        l = landunit[c]
        if !urbpoi[l]
            frac_sno = frac_sno_col[c]
            if lun_itype[l] == istice
                emg_col[c] = T(0.97)
            else
                emg_col[c] = (one(T) - frac_sno) * T(0.96) + frac_sno * T(0.97)
            end
        end
    end
end

"""
    preflux_emg!(emg_col, mask, landunit, urbpoi, lun_itype, frac_sno_col, istice)

Ground emissivity for non-urban columns, per column.
"""
preflux_emg!(emg_col, mask, landunit, urbpoi, lun_itype, frac_sno_col, istice::Int) =
    _launch!(_preflux_emg_kernel!, emg_col, mask, landunit, urbpoi, lun_itype,
             frac_sno_col, istice)

# --------------------------------------------------------------------------
# Latent heat: sublimation (HSUB) if top layer has ice but no liquid, else
# evaporation (HVAP).
# --------------------------------------------------------------------------
@kernel function _preflux_htvp_kernel!(htvp_col, @Const(mask), @Const(snl),
                                       @Const(h2osoi_liq_col), @Const(h2osoi_ice_col),
                                       nlevsno::Int, hsub, hvap)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(htvp_col)
        jtop = snl[c] + 1 + nlevsno
        h2o_liq = h2osoi_liq_col[c, jtop]
        h2o_ice = h2osoi_ice_col[c, jtop]
        htvp_col[c] = (h2o_liq <= zero(T) && h2o_ice > zero(T)) ? hsub : hvap
    end
end

"""
    preflux_htvp!(htvp_col, mask, snl, h2osoi_liq_col, h2osoi_ice_col, nlevsno, hsub, hvap)

Latent heat of vaporization/sublimation per column.
"""
function preflux_htvp!(htvp_col, mask, snl, h2osoi_liq_col, h2osoi_ice_col, nlevsno::Int, hsub, hvap)
    T = eltype(htvp_col)
    _launch!(_preflux_htvp_kernel!, htvp_col, mask, snl, h2osoi_liq_col,
             h2osoi_ice_col, nlevsno, T(hsub), T(hvap))
end

# --------------------------------------------------------------------------
# Virtual potential temperature, beta, convective BL height (zii).
# --------------------------------------------------------------------------
@kernel function _preflux_thv_kernel!(thv_col, beta_col, zii,
                                      @Const(mask), @Const(forc_th), @Const(forc_q))
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(thv_col)
        thv_col[c] = forc_th[c] * (one(T) + T(0.61) * forc_q[c])
        beta_col[c] = one(eltype(beta_col))
        zii[c] = T(1000.0)
    end
end

"""
    preflux_thv!(thv_col, beta_col, zii, mask, forc_th, forc_q)

Virtual potential temperature, beta, and convective boundary-layer height per column.
"""
preflux_thv!(thv_col, beta_col, zii, mask, forc_th, forc_q) =
    _launch!(_preflux_thv_kernel!, thv_col, beta_col, zii, mask, forc_th, forc_q)

# --------------------------------------------------------------------------
# Zero out per-patch sensible/latent heat fluxes and ground heat derivatives.
# --------------------------------------------------------------------------
@kernel function _preflux_zero_fluxes_kernel!(eflx_sh_tot, @Const(mask),
                                              eflx_sh_tot_r, eflx_sh_tot_u,
                                              eflx_lh_tot, eflx_lh_tot_r,
                                              eflx_lh_tot_u, eflx_sh_veg,
                                              cgrnd, cgrnds, cgrndl)
    p = @index(Global)
    @inbounds if mask[p]
        z = zero(eltype(eflx_sh_tot))
        eflx_sh_tot[p] = z
        eflx_sh_tot_r[p] = z
        eflx_sh_tot_u[p] = z
        eflx_lh_tot[p] = z
        eflx_lh_tot_r[p] = z
        eflx_lh_tot_u[p] = z
        eflx_sh_veg[p] = z
        cgrnd[p] = z
        cgrnds[p] = z
        cgrndl[p] = z
    end
end

"""
    preflux_zero_fluxes!(eflx_sh_tot, eflx_sh_tot_r, eflx_sh_tot_u, eflx_lh_tot,
                         eflx_lh_tot_r, eflx_lh_tot_u, eflx_sh_veg, cgrnd, cgrnds,
                         cgrndl, mask)

Zero per-patch flux/derivative arrays for active patches.
"""
function preflux_zero_fluxes!(eflx_sh_tot, eflx_sh_tot_r, eflx_sh_tot_u, eflx_lh_tot,
                              eflx_lh_tot_r, eflx_lh_tot_u, eflx_sh_veg, cgrnd, cgrnds,
                              cgrndl, mask)
    _launch!(_preflux_zero_fluxes_kernel!, eflx_sh_tot, mask, eflx_sh_tot_r,
             eflx_sh_tot_u, eflx_lh_tot, eflx_lh_tot_r, eflx_lh_tot_u, eflx_sh_veg,
             cgrnd, cgrnds, cgrndl)
end

# --------------------------------------------------------------------------
# Vegetation emissivity from exposed leaf + stem area.
# --------------------------------------------------------------------------
@kernel function _preflux_emv_kernel!(emv_patch, @Const(mask), @Const(elai_patch),
                                      @Const(esai_patch))
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(emv_patch)
        avmuir = one(T)  # inverse optical depth per unit leaf+stem area
        emv_patch[p] = one(T) - exp(-(elai_patch[p] + esai_patch[p]) / avmuir)
    end
end

"""
    preflux_emv!(emv_patch, mask, elai_patch, esai_patch)

Vegetation emissivity per patch.
"""
preflux_emv!(emv_patch, mask, elai_patch, esai_patch) =
    _launch!(_preflux_emv_kernel!, emv_patch, mask, elai_patch, esai_patch)

# --------------------------------------------------------------------------
# Intermediate temperature variable thm = forc_t(column) + lapse * forc_hgt_t.
# --------------------------------------------------------------------------
@kernel function _preflux_thm_kernel!(thm_patch, @Const(mask), @Const(column),
                                      @Const(forc_t), @Const(forc_hgt_t_patch))
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(thm_patch)
        thm_patch[p] = forc_t[column[p]] + T(0.0098) * forc_hgt_t_patch[p]
    end
end

"""
    preflux_thm!(thm_patch, mask, column, forc_t, forc_hgt_t_patch)

Intermediate temperature variable thm per patch.
"""
preflux_thm!(thm_patch, mask, column, forc_t, forc_hgt_t_patch) =
    _launch!(_preflux_thm_kernel!, thm_patch, mask, column, forc_t, forc_hgt_t_patch)

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
                          mask_nolakep::AbstractVector{Bool},
                          bounds_patch::UnitRange{Int};
                          z0param_method::String = "ZengWang2007",
                          z0mr::AbstractVector{<:Real} = pftcon.z0mr,
                          displar::AbstractVector{<:Real} = pftcon.displar)

    # z0mr/displar default to the host-resident pftcon globals; move them onto the
    # working backend so the kernel doesn't get host Arrays among device args (no-op on CPU).
    if !(canopystate.z0m_patch isa Array)
        T = eltype(canopystate.z0m_patch)
        z0mr = copyto!(similar(canopystate.z0m_patch, T, length(z0mr)), T.(z0mr))
        displar = copyto!(similar(canopystate.z0m_patch, T, length(displar)), T.(displar))
    end

    # Both z0param_method branches are identical, so the common expression is used.
    preflux_z0m_displa!(canopystate.z0m_patch, canopystate.displa_patch,
                        mask_nolakep, patch_data.itype, canopystate.htop_patch,
                        z0mr, displar)

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
                                    forc_t::AbstractVector{<:Real},
                                    forc_th::AbstractVector{<:Real},
                                    forc_q::AbstractVector{<:Real},
                                    mask_nolakec::AbstractVector{Bool},
                                    mask_nolakep::AbstractVector{Bool},
                                    bounds_col::UnitRange{Int},
                                    bounds_patch::UnitRange{Int})

    nlevsno_val = varpar.nlevsno
    nlevgrnd_val = varpar.nlevgrnd

    # ---------------------------------------------------------------
    # 1. Save current temperatures (tssbef) for flux calculations
    # ---------------------------------------------------------------
    preflux_tssbef!(temperature.t_ssbef_col, mask_nolakec, temperature.t_soisno_col,
                    nlevsno_val + nlevgrnd_val)

    # ---------------------------------------------------------------
    # 2. Ground temperature: weighted average of snow, soil, surface water
    # ---------------------------------------------------------------
    preflux_t_grnd!(temperature.t_grnd_col, temperature.t_h2osfc_bef_col, mask_nolakec,
                    col_data.snl, waterdiagbulk.frac_sno_eff_col, waterdiagbulk.frac_h2osfc_col,
                    temperature.t_soisno_col, temperature.t_h2osfc_col, nlevsno_val)

    # ---------------------------------------------------------------
    # 3. Ground emissivity
    # ---------------------------------------------------------------
    # Urban emissivity is set elsewhere (from surface data)
    preflux_emg!(temperature.emg_col, mask_nolakec, col_data.landunit,
                 lun_data.urbpoi, lun_data.itype, waterdiagbulk.frac_sno_col, ISTICE)

    # ---------------------------------------------------------------
    # 4. Latent heat of vaporization/sublimation
    # ---------------------------------------------------------------
    preflux_htvp!(energyflux.htvp_col, mask_nolakec, col_data.snl,
                  waterstatebulk.ws.h2osoi_liq_col, waterstatebulk.ws.h2osoi_ice_col,
                  nlevsno_val, HSUB, HVAP)

    # ---------------------------------------------------------------
    # 5. Virtual potential temperature, beta, zii
    # ---------------------------------------------------------------
    preflux_thv!(temperature.thv_col, temperature.beta_col, col_data.zii,
                 mask_nolakec, forc_th, forc_q)

    # ---------------------------------------------------------------
    # 6. Initialize energy fluxes to zero
    # ---------------------------------------------------------------
    preflux_zero_fluxes!(energyflux.eflx_sh_tot_patch, energyflux.eflx_sh_tot_r_patch,
                         energyflux.eflx_sh_tot_u_patch, energyflux.eflx_lh_tot_patch,
                         energyflux.eflx_lh_tot_r_patch, energyflux.eflx_lh_tot_u_patch,
                         energyflux.eflx_sh_veg_patch, energyflux.cgrnd_patch,
                         energyflux.cgrnds_patch, energyflux.cgrndl_patch, mask_nolakep)

    # ---------------------------------------------------------------
    # 7. Vegetation emissivity
    # ---------------------------------------------------------------
    preflux_emv!(temperature.emv_patch, mask_nolakep, canopystate.elai_patch,
                 canopystate.esai_patch)

    # ---------------------------------------------------------------
    # 8. Intermediate temperature variable thm
    # ---------------------------------------------------------------
    preflux_thm!(temperature.thm_patch, mask_nolakep, patch_data.column,
                 forc_t, frictionvel.forc_hgt_t_patch)

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
                                     forc_t::AbstractVector{<:Real},
                                     forc_th::AbstractVector{<:Real},
                                     forc_q::AbstractVector{<:Real},
                                     forc_hgt_u::AbstractVector{<:Real},
                                     forc_hgt_t::AbstractVector{<:Real},
                                     forc_hgt_q::AbstractVector{<:Real},
                                     mask_nolakec::AbstractVector{Bool},
                                     mask_nolakep::AbstractVector{Bool},
                                     mask_urbanc::AbstractVector{Bool},
                                     bounds_col::UnitRange{Int},
                                     bounds_patch::UnitRange{Int};
                                     z0param_method::String = "ZengWang2007")

    # 1. Set z0m and displacement height
    set_z0m_displa!(canopystate, frictionvel, patch_data, lun_data,
                    mask_nolakep, bounds_patch;
                    z0param_method=z0param_method)

    # 2. Set roughness lengths and forcing heights for non-lake points.
    # snomelt_accum is a placeholder (only read on the Meier2022 + snowmelt path);
    # allocate it device-resident (similar+fill!) so it matches the backend of the
    # state arrays and nothing host-side leaks onto the device.
    snomelt_accum = similar(waterdiagbulk.frac_sno_col)
    fill!(snomelt_accum, zero(eltype(snomelt_accum)))
    set_roughness_and_forc_heights_nonlake!(
        frictionvel,
        mask_nolakec, mask_nolakep,
        bounds_col, bounds_patch,
        waterdiagbulk.frac_sno_col,
        snomelt_accum,
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
