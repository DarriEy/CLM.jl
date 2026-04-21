# ==========================================================================
# Ported from: src/biogeochem/CNCStateUpdate2Mod.F90
# Carbon state variable update, mortality fluxes.
#
# When the matrix solution is being used (use_matrixcn and use_soil_matrixcn)
# only some state updates are done here; other state updates happen after the
# matrix is solved in VegMatrix and SoilMatrix.
#
# Public functions:
#   c_state_update2!  — Gap-phase mortality C state update
#   c_state_update2h! — Harvest mortality C state update
#   c_state_update2g! — Gross unrepresented landcover change mortality C state update
# ==========================================================================

# ---------------------------------------------------------------------------
# c_state_update2! — Gap-phase mortality
# ---------------------------------------------------------------------------

"""
    c_state_update2!(cs_veg, cf_veg, cs_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

On the radiation time step, update all prognostic carbon state variables
affected by gap-phase mortality fluxes.

Ported from `CStateUpdate2` in `CNCStateUpdate2Mod.F90`.

Note: The matrix-CN code paths (`use_soil_matrixcn`) for soil matrix
solutions are omitted since they require the full matrix infrastructure
(matrix_Cinput). The matrix-CN veg path (`use_matrixcn`) is a no-op in
the Fortran original (comment only), so it is omitted here as well.
"""
function c_state_update2!(cs_veg::CNVegCarbonStateData,
                           cf_veg::CNVegCarbonFluxData,
                           cs_soil::SoilBiogeochemCarbonStateData;
                           mask_soilc::BitVector,
                           mask_soilp::BitVector,
                           bounds_col::UnitRange{Int},
                           bounds_patch::UnitRange{Int},
                           nlevdecomp::Int,
                           i_litr_min::Int,
                           i_litr_max::Int,
                           i_cwd::Int,
                           use_matrixcn::Bool = false,
                           use_soil_matrixcn::Bool = false,
                           dt::Real)

    # --- Column loop: gap-phase mortality fluxes to soil decomposition pools ---
    for j in 1:nlevdecomp
        for c in bounds_col
            mask_soilc[c] || continue

            if !use_soil_matrixcn
                for i in i_litr_min:i_litr_max
                    cs_soil.decomp_cpools_vr_col[c, j, i] =
                        cs_soil.decomp_cpools_vr_col[c, j, i] +
                        cf_veg.gap_mortality_c_to_litr_c_col[c, j, i] * dt
                end
                # Currently i_cwd != i_litr_max + 1 if not fates and
                # i_cwd = 0 if fates, so not including in the i-loop
                cs_soil.decomp_cpools_vr_col[c, j, i_cwd] =
                    cs_soil.decomp_cpools_vr_col[c, j, i_cwd] +
                    cf_veg.gap_mortality_c_to_cwdc_col[c, j] * dt
            else
                # Matrix path: matrix_Cinput setup omitted (requires full matrix infrastructure)
            end
        end
    end

    # --- Patch loop: vegetation C state updates from gap-phase mortality ---
    for p in bounds_patch
        mask_soilp[p] || continue

        cs_veg.gresp_storage_patch[p] = cs_veg.gresp_storage_patch[p] -
            cf_veg.m_gresp_storage_to_litter_patch[p] * dt
        cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] -
            cf_veg.m_gresp_xfer_to_litter_patch[p] * dt

        if !use_matrixcn
            # displayed pools
            cs_veg.leafc_patch[p] = cs_veg.leafc_patch[p] -
                cf_veg.m_leafc_to_litter_patch[p] * dt
            cs_veg.frootc_patch[p] = cs_veg.frootc_patch[p] -
                cf_veg.m_frootc_to_litter_patch[p] * dt
            cs_veg.livestemc_patch[p] = cs_veg.livestemc_patch[p] -
                cf_veg.m_livestemc_to_litter_patch[p] * dt
            cs_veg.deadstemc_patch[p] = cs_veg.deadstemc_patch[p] -
                cf_veg.m_deadstemc_to_litter_patch[p] * dt
            cs_veg.livecrootc_patch[p] = cs_veg.livecrootc_patch[p] -
                cf_veg.m_livecrootc_to_litter_patch[p] * dt
            cs_veg.deadcrootc_patch[p] = cs_veg.deadcrootc_patch[p] -
                cf_veg.m_deadcrootc_to_litter_patch[p] * dt

            # storage pools
            cs_veg.leafc_storage_patch[p] = cs_veg.leafc_storage_patch[p] -
                cf_veg.m_leafc_storage_to_litter_patch[p] * dt
            cs_veg.frootc_storage_patch[p] = cs_veg.frootc_storage_patch[p] -
                cf_veg.m_frootc_storage_to_litter_patch[p] * dt
            cs_veg.livestemc_storage_patch[p] = cs_veg.livestemc_storage_patch[p] -
                cf_veg.m_livestemc_storage_to_litter_patch[p] * dt
            cs_veg.deadstemc_storage_patch[p] = cs_veg.deadstemc_storage_patch[p] -
                cf_veg.m_deadstemc_storage_to_litter_patch[p] * dt
            cs_veg.livecrootc_storage_patch[p] = cs_veg.livecrootc_storage_patch[p] -
                cf_veg.m_livecrootc_storage_to_litter_patch[p] * dt
            cs_veg.deadcrootc_storage_patch[p] = cs_veg.deadcrootc_storage_patch[p] -
                cf_veg.m_deadcrootc_storage_to_litter_patch[p] * dt

            # transfer pools
            cs_veg.leafc_xfer_patch[p] = cs_veg.leafc_xfer_patch[p] -
                cf_veg.m_leafc_xfer_to_litter_patch[p] * dt
            cs_veg.frootc_xfer_patch[p] = cs_veg.frootc_xfer_patch[p] -
                cf_veg.m_frootc_xfer_to_litter_patch[p] * dt
            cs_veg.livestemc_xfer_patch[p] = cs_veg.livestemc_xfer_patch[p] -
                cf_veg.m_livestemc_xfer_to_litter_patch[p] * dt
            cs_veg.deadstemc_xfer_patch[p] = cs_veg.deadstemc_xfer_patch[p] -
                cf_veg.m_deadstemc_xfer_to_litter_patch[p] * dt
            cs_veg.livecrootc_xfer_patch[p] = cs_veg.livecrootc_xfer_patch[p] -
                cf_veg.m_livecrootc_xfer_to_litter_patch[p] * dt
            cs_veg.deadcrootc_xfer_patch[p] = cs_veg.deadcrootc_xfer_patch[p] -
                cf_veg.m_deadcrootc_xfer_to_litter_patch[p] * dt
        else
            # Matrix version: handled in CNGapMortality
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# c_state_update2h! — Harvest mortality
# ---------------------------------------------------------------------------

"""
    c_state_update2h!(cs_veg, cf_veg, cs_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

Update all prognostic carbon state variables affected by harvest
mortality fluxes.

Ported from `CStateUpdate2h` in `CNCStateUpdate2Mod.F90`.

Note: The matrix-CN code paths are omitted (require full matrix infrastructure).
Wood to product pools state updates happen in CNProducts (not here).
"""
function c_state_update2h!(cs_veg::CNVegCarbonStateData,
                            cf_veg::CNVegCarbonFluxData,
                            cs_soil::SoilBiogeochemCarbonStateData;
                            mask_soilc::BitVector,
                            mask_soilp::BitVector,
                            bounds_col::UnitRange{Int},
                            bounds_patch::UnitRange{Int},
                            nlevdecomp::Int,
                            i_litr_min::Int,
                            i_litr_max::Int,
                            i_cwd::Int,
                            use_matrixcn::Bool = false,
                            use_soil_matrixcn::Bool = false,
                            dt::Real)

    # --- Column loop: harvest mortality fluxes to soil decomposition pools ---
    for j in 1:nlevdecomp
        for c in bounds_col
            mask_soilc[c] || continue

            if !use_soil_matrixcn
                for i in i_litr_min:i_litr_max
                    cs_soil.decomp_cpools_vr_col[c, j, i] =
                        cs_soil.decomp_cpools_vr_col[c, j, i] +
                        cf_veg.harvest_c_to_litr_c_col[c, j, i] * dt
                end
                # Currently i_cwd != i_litr_max + 1 if not fates and
                # i_cwd = 0 if fates, so not including in the i-loop
                cs_soil.decomp_cpools_vr_col[c, j, i_cwd] =
                    cs_soil.decomp_cpools_vr_col[c, j, i_cwd] +
                    cf_veg.harvest_c_to_cwdc_col[c, j] * dt
            else
                # Matrix path: matrix_Cinput setup omitted (requires full matrix infrastructure)
            end

            # wood to product pools - states updated in CNProducts
        end
    end

    # --- Patch loop: vegetation C state updates from harvest mortality ---
    for p in bounds_patch
        mask_soilp[p] || continue

        # xsmrpool
        cs_veg.xsmrpool_patch[p] = cs_veg.xsmrpool_patch[p] -
            cf_veg.hrv_xsmrpool_to_atm_patch[p] * dt

        # storage pools
        cs_veg.gresp_storage_patch[p] = cs_veg.gresp_storage_patch[p] -
            cf_veg.hrv_gresp_storage_to_litter_patch[p] * dt

        # transfer pools
        cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] -
            cf_veg.hrv_gresp_xfer_to_litter_patch[p] * dt

        if !use_matrixcn
            # displayed pools
            cs_veg.leafc_patch[p] = cs_veg.leafc_patch[p] -
                cf_veg.hrv_leafc_to_litter_patch[p] * dt
            cs_veg.frootc_patch[p] = cs_veg.frootc_patch[p] -
                cf_veg.hrv_frootc_to_litter_patch[p] * dt
            cs_veg.livestemc_patch[p] = cs_veg.livestemc_patch[p] -
                cf_veg.hrv_livestemc_to_litter_patch[p] * dt
            cs_veg.deadstemc_patch[p] = cs_veg.deadstemc_patch[p] -
                cf_veg.wood_harvestc_patch[p] * dt
            cs_veg.livecrootc_patch[p] = cs_veg.livecrootc_patch[p] -
                cf_veg.hrv_livecrootc_to_litter_patch[p] * dt
            cs_veg.deadcrootc_patch[p] = cs_veg.deadcrootc_patch[p] -
                cf_veg.hrv_deadcrootc_to_litter_patch[p] * dt

            # storage pools
            cs_veg.leafc_storage_patch[p] = cs_veg.leafc_storage_patch[p] -
                cf_veg.hrv_leafc_storage_to_litter_patch[p] * dt
            cs_veg.frootc_storage_patch[p] = cs_veg.frootc_storage_patch[p] -
                cf_veg.hrv_frootc_storage_to_litter_patch[p] * dt
            cs_veg.livestemc_storage_patch[p] = cs_veg.livestemc_storage_patch[p] -
                cf_veg.hrv_livestemc_storage_to_litter_patch[p] * dt
            cs_veg.deadstemc_storage_patch[p] = cs_veg.deadstemc_storage_patch[p] -
                cf_veg.hrv_deadstemc_storage_to_litter_patch[p] * dt
            cs_veg.livecrootc_storage_patch[p] = cs_veg.livecrootc_storage_patch[p] -
                cf_veg.hrv_livecrootc_storage_to_litter_patch[p] * dt
            cs_veg.deadcrootc_storage_patch[p] = cs_veg.deadcrootc_storage_patch[p] -
                cf_veg.hrv_deadcrootc_storage_to_litter_patch[p] * dt

            # transfer pools
            cs_veg.leafc_xfer_patch[p] = cs_veg.leafc_xfer_patch[p] -
                cf_veg.hrv_leafc_xfer_to_litter_patch[p] * dt
            cs_veg.frootc_xfer_patch[p] = cs_veg.frootc_xfer_patch[p] -
                cf_veg.hrv_frootc_xfer_to_litter_patch[p] * dt
            cs_veg.livestemc_xfer_patch[p] = cs_veg.livestemc_xfer_patch[p] -
                cf_veg.hrv_livestemc_xfer_to_litter_patch[p] * dt
            cs_veg.deadstemc_xfer_patch[p] = cs_veg.deadstemc_xfer_patch[p] -
                cf_veg.hrv_deadstemc_xfer_to_litter_patch[p] * dt
            cs_veg.livecrootc_xfer_patch[p] = cs_veg.livecrootc_xfer_patch[p] -
                cf_veg.hrv_livecrootc_xfer_to_litter_patch[p] * dt
            cs_veg.deadcrootc_xfer_patch[p] = cs_veg.deadcrootc_xfer_patch[p] -
                cf_veg.hrv_deadcrootc_xfer_to_litter_patch[p] * dt
        else
            # Matrix version: handled in CNHarvest
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# c_state_update2g! — Gross unrepresented landcover change mortality
# ---------------------------------------------------------------------------

"""
    c_state_update2g!(cs_veg, cf_veg, cs_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

Update all prognostic carbon state variables affected by gross
unrepresented landcover change mortality fluxes.

Ported from `CStateUpdate2g` in `CNCStateUpdate2Mod.F90`.

Note: The matrix-CN code paths are omitted (require full matrix infrastructure).
Wood to product pools state updates happen in CNProducts (not here).
"""
function c_state_update2g!(cs_veg::CNVegCarbonStateData,
                            cf_veg::CNVegCarbonFluxData,
                            cs_soil::SoilBiogeochemCarbonStateData;
                            mask_soilc::BitVector,
                            mask_soilp::BitVector,
                            bounds_col::UnitRange{Int},
                            bounds_patch::UnitRange{Int},
                            nlevdecomp::Int,
                            i_litr_min::Int,
                            i_litr_max::Int,
                            i_cwd::Int,
                            use_matrixcn::Bool = false,
                            use_soil_matrixcn::Bool = false,
                            dt::Real)

    # --- Column loop: gross unrepresented landcover change fluxes to soil ---
    for j in 1:nlevdecomp
        for c in bounds_col
            mask_soilc[c] || continue

            if !use_soil_matrixcn
                for i in i_litr_min:i_litr_max
                    cs_soil.decomp_cpools_vr_col[c, j, i] =
                        cs_soil.decomp_cpools_vr_col[c, j, i] +
                        cf_veg.gru_c_to_litr_c_col[c, j, i] * dt
                end
                cs_soil.decomp_cpools_vr_col[c, j, i_cwd] =
                    cs_soil.decomp_cpools_vr_col[c, j, i_cwd] +
                    cf_veg.gru_c_to_cwdc_col[c, j] * dt

                # wood to product pools - states updated in CNProducts
            else
                # Matrix path: matrix_Cinput setup omitted (requires full matrix infrastructure)
            end
        end
    end

    # --- Patch loop: vegetation C state updates from gross unrepresented landcover change ---
    for p in bounds_patch
        mask_soilp[p] || continue

        # xsmrpool
        cs_veg.xsmrpool_patch[p] = cs_veg.xsmrpool_patch[p] -
            cf_veg.gru_xsmrpool_to_atm_patch[p] * dt
        # gresp storage pool
        cs_veg.gresp_storage_patch[p] = cs_veg.gresp_storage_patch[p] -
            cf_veg.gru_gresp_storage_to_atm_patch[p] * dt
        # gresp transfer pool
        cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] -
            cf_veg.gru_gresp_xfer_to_atm_patch[p] * dt

        if !use_matrixcn
            # displayed pools
            cs_veg.leafc_patch[p] = cs_veg.leafc_patch[p] -
                cf_veg.gru_leafc_to_litter_patch[p] * dt
            cs_veg.frootc_patch[p] = cs_veg.frootc_patch[p] -
                cf_veg.gru_frootc_to_litter_patch[p] * dt
            cs_veg.livestemc_patch[p] = cs_veg.livestemc_patch[p] -
                cf_veg.gru_livestemc_to_atm_patch[p] * dt
            cs_veg.deadstemc_patch[p] = cs_veg.deadstemc_patch[p] -
                cf_veg.gru_deadstemc_to_atm_patch[p] * dt
            cs_veg.deadstemc_patch[p] = cs_veg.deadstemc_patch[p] -
                cf_veg.gru_wood_productc_gain_patch[p] * dt
            cs_veg.livecrootc_patch[p] = cs_veg.livecrootc_patch[p] -
                cf_veg.gru_livecrootc_to_litter_patch[p] * dt
            cs_veg.deadcrootc_patch[p] = cs_veg.deadcrootc_patch[p] -
                cf_veg.gru_deadcrootc_to_litter_patch[p] * dt

            # storage pools
            cs_veg.leafc_storage_patch[p] = cs_veg.leafc_storage_patch[p] -
                cf_veg.gru_leafc_storage_to_atm_patch[p] * dt
            cs_veg.frootc_storage_patch[p] = cs_veg.frootc_storage_patch[p] -
                cf_veg.gru_frootc_storage_to_atm_patch[p] * dt
            cs_veg.livestemc_storage_patch[p] = cs_veg.livestemc_storage_patch[p] -
                cf_veg.gru_livestemc_storage_to_atm_patch[p] * dt
            cs_veg.deadstemc_storage_patch[p] = cs_veg.deadstemc_storage_patch[p] -
                cf_veg.gru_deadstemc_storage_to_atm_patch[p] * dt
            cs_veg.livecrootc_storage_patch[p] = cs_veg.livecrootc_storage_patch[p] -
                cf_veg.gru_livecrootc_storage_to_atm_patch[p] * dt
            cs_veg.deadcrootc_storage_patch[p] = cs_veg.deadcrootc_storage_patch[p] -
                cf_veg.gru_deadcrootc_storage_to_atm_patch[p] * dt

            # transfer pools
            cs_veg.leafc_xfer_patch[p] = cs_veg.leafc_xfer_patch[p] -
                cf_veg.gru_leafc_xfer_to_atm_patch[p] * dt
            cs_veg.frootc_xfer_patch[p] = cs_veg.frootc_xfer_patch[p] -
                cf_veg.gru_frootc_xfer_to_atm_patch[p] * dt
            cs_veg.livestemc_xfer_patch[p] = cs_veg.livestemc_xfer_patch[p] -
                cf_veg.gru_livestemc_xfer_to_atm_patch[p] * dt
            cs_veg.deadstemc_xfer_patch[p] = cs_veg.deadstemc_xfer_patch[p] -
                cf_veg.gru_deadstemc_xfer_to_atm_patch[p] * dt
            cs_veg.livecrootc_xfer_patch[p] = cs_veg.livecrootc_xfer_patch[p] -
                cf_veg.gru_livecrootc_xfer_to_atm_patch[p] * dt
            cs_veg.deadcrootc_xfer_patch[p] = cs_veg.deadcrootc_xfer_patch[p] -
                cf_veg.gru_deadcrootc_xfer_to_atm_patch[p] * dt
        else
            # Matrix version: handled in dynGrossUnrepMod::CNGrossUnrep
        end
    end

    return nothing
end
