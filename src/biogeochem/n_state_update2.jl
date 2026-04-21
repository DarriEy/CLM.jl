# ==========================================================================
# Ported from: src/biogeochem/CNNStateUpdate2Mod.F90
# Nitrogen state variable update, mortality fluxes.
#
# When the matrix solution is being used (use_matrixcn and use_soil_matrixcn)
# only some state updates are done here, the other state updates happen
# after the matrix is solved in VegMatrix and SoilMatrix.
#
# Public functions:
#   n_state_update2!  — Gap-phase mortality N state update
#   n_state_update2h! — Harvest mortality N state update
#   n_state_update2g! — Gross unrepresented landcover change mortality N state update
# ==========================================================================

# ---------------------------------------------------------------------------
# n_state_update2! — Gap-phase mortality
# ---------------------------------------------------------------------------

"""
    n_state_update2!(ns_veg, nf_veg, ns_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

On the radiation time step, update all prognostic nitrogen state variables
affected by gap-phase mortality fluxes.

Ported from `NStateUpdate2` in `CNNStateUpdate2Mod.F90`.

Note: The matrix-CN code paths (`use_soil_matrixcn`) for soil matrix
solutions are omitted since they require the full matrix infrastructure
(matrix_Ninput). The matrix-CN veg path (`use_matrixcn`) is a no-op in
the Fortran original (comment only), so it is omitted here as well.
"""
function n_state_update2!(ns_veg::CNVegNitrogenStateData,
                           nf_veg::CNVegNitrogenFluxData,
                           ns_soil::SoilBiogeochemNitrogenStateData;
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
                    ns_soil.decomp_npools_vr_col[c, j, i] =
                        ns_soil.decomp_npools_vr_col[c, j, i] +
                        nf_veg.gap_mortality_n_to_litr_n_col[c, j, i] * dt
                end
                # Currently i_cwd != i_litr_max + 1 if not fates and
                # i_cwd = 0 if fates, so not including in the i-loop
                ns_soil.decomp_npools_vr_col[c, j, i_cwd] =
                    ns_soil.decomp_npools_vr_col[c, j, i_cwd] +
                    nf_veg.gap_mortality_n_to_cwdn_col[c, j] * dt
            else
                # Matrix path: matrix_Ninput setup omitted (requires full matrix infrastructure)
            end
        end
    end

    # --- Patch loop: vegetation N state updates from gap-phase mortality ---
    for p in bounds_patch
        mask_soilp[p] || continue

        if !use_matrixcn
            # displayed pools
            ns_veg.leafn_patch[p] = ns_veg.leafn_patch[p] -
                nf_veg.m_leafn_to_litter_patch[p] * dt
            ns_veg.frootn_patch[p] = ns_veg.frootn_patch[p] -
                nf_veg.m_frootn_to_litter_patch[p] * dt
            ns_veg.livestemn_patch[p] = ns_veg.livestemn_patch[p] -
                nf_veg.m_livestemn_to_litter_patch[p] * dt
            ns_veg.deadstemn_patch[p] = ns_veg.deadstemn_patch[p] -
                nf_veg.m_deadstemn_to_litter_patch[p] * dt
            ns_veg.livecrootn_patch[p] = ns_veg.livecrootn_patch[p] -
                nf_veg.m_livecrootn_to_litter_patch[p] * dt
            ns_veg.deadcrootn_patch[p] = ns_veg.deadcrootn_patch[p] -
                nf_veg.m_deadcrootn_to_litter_patch[p] * dt
            ns_veg.retransn_patch[p] = ns_veg.retransn_patch[p] -
                nf_veg.m_retransn_to_litter_patch[p] * dt

            # storage pools
            ns_veg.leafn_storage_patch[p] = ns_veg.leafn_storage_patch[p] -
                nf_veg.m_leafn_storage_to_litter_patch[p] * dt
            ns_veg.frootn_storage_patch[p] = ns_veg.frootn_storage_patch[p] -
                nf_veg.m_frootn_storage_to_litter_patch[p] * dt
            ns_veg.livestemn_storage_patch[p] = ns_veg.livestemn_storage_patch[p] -
                nf_veg.m_livestemn_storage_to_litter_patch[p] * dt
            ns_veg.deadstemn_storage_patch[p] = ns_veg.deadstemn_storage_patch[p] -
                nf_veg.m_deadstemn_storage_to_litter_patch[p] * dt
            ns_veg.livecrootn_storage_patch[p] = ns_veg.livecrootn_storage_patch[p] -
                nf_veg.m_livecrootn_storage_to_litter_patch[p] * dt
            ns_veg.deadcrootn_storage_patch[p] = ns_veg.deadcrootn_storage_patch[p] -
                nf_veg.m_deadcrootn_storage_to_litter_patch[p] * dt

            # transfer pools
            ns_veg.leafn_xfer_patch[p] = ns_veg.leafn_xfer_patch[p] -
                nf_veg.m_leafn_xfer_to_litter_patch[p] * dt
            ns_veg.frootn_xfer_patch[p] = ns_veg.frootn_xfer_patch[p] -
                nf_veg.m_frootn_xfer_to_litter_patch[p] * dt
            ns_veg.livestemn_xfer_patch[p] = ns_veg.livestemn_xfer_patch[p] -
                nf_veg.m_livestemn_xfer_to_litter_patch[p] * dt
            ns_veg.deadstemn_xfer_patch[p] = ns_veg.deadstemn_xfer_patch[p] -
                nf_veg.m_deadstemn_xfer_to_litter_patch[p] * dt
            ns_veg.livecrootn_xfer_patch[p] = ns_veg.livecrootn_xfer_patch[p] -
                nf_veg.m_livecrootn_xfer_to_litter_patch[p] * dt
            ns_veg.deadcrootn_xfer_patch[p] = ns_veg.deadcrootn_xfer_patch[p] -
                nf_veg.m_deadcrootn_xfer_to_litter_patch[p] * dt
        else
            # Matrix version: handled in CNGapMortality
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# n_state_update2h! — Harvest mortality
# ---------------------------------------------------------------------------

"""
    n_state_update2h!(ns_veg, nf_veg, ns_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

Update all prognostic nitrogen state variables affected by harvest
mortality fluxes.

Ported from `NStateUpdate2h` in `CNNStateUpdate2Mod.F90`.

Note: The matrix-CN code paths are omitted (require full matrix infrastructure).
"""
function n_state_update2h!(ns_veg::CNVegNitrogenStateData,
                            nf_veg::CNVegNitrogenFluxData,
                            ns_soil::SoilBiogeochemNitrogenStateData;
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
                    ns_soil.decomp_npools_vr_col[c, j, i] =
                        ns_soil.decomp_npools_vr_col[c, j, i] +
                        nf_veg.harvest_n_to_litr_n_col[c, j, i] * dt
                end
                # Currently i_cwd != i_litr_max + 1 if not fates and
                # i_cwd = 0 if fates, so not including in the i-loop
                ns_soil.decomp_npools_vr_col[c, j, i_cwd] =
                    ns_soil.decomp_npools_vr_col[c, j, i_cwd] +
                    nf_veg.harvest_n_to_cwdn_col[c, j] * dt
            else
                # Matrix path: matrix_Ninput setup omitted (requires full matrix infrastructure)
            end
        end
    end

    # --- Patch loop: vegetation N state updates from harvest mortality ---
    for p in bounds_patch
        mask_soilp[p] || continue

        if !use_matrixcn
            # displayed pools
            ns_veg.leafn_patch[p] = ns_veg.leafn_patch[p] -
                nf_veg.hrv_leafn_to_litter_patch[p] * dt
            ns_veg.frootn_patch[p] = ns_veg.frootn_patch[p] -
                nf_veg.hrv_frootn_to_litter_patch[p] * dt
            ns_veg.livestemn_patch[p] = ns_veg.livestemn_patch[p] -
                nf_veg.hrv_livestemn_to_litter_patch[p] * dt
            ns_veg.deadstemn_patch[p] = ns_veg.deadstemn_patch[p] -
                nf_veg.wood_harvestn_patch[p] * dt
            ns_veg.livecrootn_patch[p] = ns_veg.livecrootn_patch[p] -
                nf_veg.hrv_livecrootn_to_litter_patch[p] * dt
            ns_veg.deadcrootn_patch[p] = ns_veg.deadcrootn_patch[p] -
                nf_veg.hrv_deadcrootn_to_litter_patch[p] * dt
            ns_veg.retransn_patch[p] = ns_veg.retransn_patch[p] -
                nf_veg.hrv_retransn_to_litter_patch[p] * dt

            # storage pools
            ns_veg.leafn_storage_patch[p] = ns_veg.leafn_storage_patch[p] -
                nf_veg.hrv_leafn_storage_to_litter_patch[p] * dt
            ns_veg.frootn_storage_patch[p] = ns_veg.frootn_storage_patch[p] -
                nf_veg.hrv_frootn_storage_to_litter_patch[p] * dt
            ns_veg.livestemn_storage_patch[p] = ns_veg.livestemn_storage_patch[p] -
                nf_veg.hrv_livestemn_storage_to_litter_patch[p] * dt
            ns_veg.deadstemn_storage_patch[p] = ns_veg.deadstemn_storage_patch[p] -
                nf_veg.hrv_deadstemn_storage_to_litter_patch[p] * dt
            ns_veg.livecrootn_storage_patch[p] = ns_veg.livecrootn_storage_patch[p] -
                nf_veg.hrv_livecrootn_storage_to_litter_patch[p] * dt
            ns_veg.deadcrootn_storage_patch[p] = ns_veg.deadcrootn_storage_patch[p] -
                nf_veg.hrv_deadcrootn_storage_to_litter_patch[p] * dt

            # transfer pools
            ns_veg.leafn_xfer_patch[p] = ns_veg.leafn_xfer_patch[p] -
                nf_veg.hrv_leafn_xfer_to_litter_patch[p] * dt
            ns_veg.frootn_xfer_patch[p] = ns_veg.frootn_xfer_patch[p] -
                nf_veg.hrv_frootn_xfer_to_litter_patch[p] * dt
            ns_veg.livestemn_xfer_patch[p] = ns_veg.livestemn_xfer_patch[p] -
                nf_veg.hrv_livestemn_xfer_to_litter_patch[p] * dt
            ns_veg.deadstemn_xfer_patch[p] = ns_veg.deadstemn_xfer_patch[p] -
                nf_veg.hrv_deadstemn_xfer_to_litter_patch[p] * dt
            ns_veg.livecrootn_xfer_patch[p] = ns_veg.livecrootn_xfer_patch[p] -
                nf_veg.hrv_livecrootn_xfer_to_litter_patch[p] * dt
            ns_veg.deadcrootn_xfer_patch[p] = ns_veg.deadcrootn_xfer_patch[p] -
                nf_veg.hrv_deadcrootn_xfer_to_litter_patch[p] * dt
        else
            # Matrix version: handled in CNHarvest
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# n_state_update2g! — Gross unrepresented landcover change mortality
# ---------------------------------------------------------------------------

"""
    n_state_update2g!(ns_veg, nf_veg, ns_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

Update all prognostic nitrogen state variables affected by gross
unrepresented landcover change mortality fluxes.

Ported from `NStateUpdate2g` in `CNNStateUpdate2Mod.F90`.

Note: The matrix-CN code paths are omitted (require full matrix infrastructure).
"""
function n_state_update2g!(ns_veg::CNVegNitrogenStateData,
                            nf_veg::CNVegNitrogenFluxData,
                            ns_soil::SoilBiogeochemNitrogenStateData;
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
                    ns_soil.decomp_npools_vr_col[c, j, i] =
                        ns_soil.decomp_npools_vr_col[c, j, i] +
                        nf_veg.gru_n_to_litr_n_col[c, j, i] * dt
                end
                # Currently i_cwd != i_litr_max + 1 if not fates and
                # i_cwd = 0 if fates, so not including in the i-loop
                ns_soil.decomp_npools_vr_col[c, j, i_cwd] =
                    ns_soil.decomp_npools_vr_col[c, j, i_cwd] +
                    nf_veg.gru_n_to_cwdn_col[c, j] * dt
            else
                # Matrix path: matrix_Ninput setup omitted (requires full matrix infrastructure)
            end
        end
    end

    # --- Patch loop: vegetation N state updates from gross unrepresented landcover change ---
    for p in bounds_patch
        mask_soilp[p] || continue

        if !use_matrixcn
            # displayed pools
            ns_veg.leafn_patch[p] = ns_veg.leafn_patch[p] -
                nf_veg.gru_leafn_to_litter_patch[p] * dt
            ns_veg.frootn_patch[p] = ns_veg.frootn_patch[p] -
                nf_veg.gru_frootn_to_litter_patch[p] * dt
            ns_veg.livestemn_patch[p] = ns_veg.livestemn_patch[p] -
                nf_veg.gru_livestemn_to_atm_patch[p] * dt
            ns_veg.deadstemn_patch[p] = ns_veg.deadstemn_patch[p] -
                nf_veg.gru_deadstemn_to_atm_patch[p] * dt
            ns_veg.deadstemn_patch[p] = ns_veg.deadstemn_patch[p] -
                nf_veg.gru_wood_productn_gain_patch[p] * dt
            ns_veg.livecrootn_patch[p] = ns_veg.livecrootn_patch[p] -
                nf_veg.gru_livecrootn_to_litter_patch[p] * dt
            ns_veg.deadcrootn_patch[p] = ns_veg.deadcrootn_patch[p] -
                nf_veg.gru_deadcrootn_to_litter_patch[p] * dt
            ns_veg.retransn_patch[p] = ns_veg.retransn_patch[p] -
                nf_veg.gru_retransn_to_litter_patch[p] * dt

            # storage pools
            ns_veg.leafn_storage_patch[p] = ns_veg.leafn_storage_patch[p] -
                nf_veg.gru_leafn_storage_to_atm_patch[p] * dt
            ns_veg.frootn_storage_patch[p] = ns_veg.frootn_storage_patch[p] -
                nf_veg.gru_frootn_storage_to_atm_patch[p] * dt
            ns_veg.livestemn_storage_patch[p] = ns_veg.livestemn_storage_patch[p] -
                nf_veg.gru_livestemn_storage_to_atm_patch[p] * dt
            ns_veg.deadstemn_storage_patch[p] = ns_veg.deadstemn_storage_patch[p] -
                nf_veg.gru_deadstemn_storage_to_atm_patch[p] * dt
            ns_veg.livecrootn_storage_patch[p] = ns_veg.livecrootn_storage_patch[p] -
                nf_veg.gru_livecrootn_storage_to_atm_patch[p] * dt
            ns_veg.deadcrootn_storage_patch[p] = ns_veg.deadcrootn_storage_patch[p] -
                nf_veg.gru_deadcrootn_storage_to_atm_patch[p] * dt

            # transfer pools
            ns_veg.leafn_xfer_patch[p] = ns_veg.leafn_xfer_patch[p] -
                nf_veg.gru_leafn_xfer_to_atm_patch[p] * dt
            ns_veg.frootn_xfer_patch[p] = ns_veg.frootn_xfer_patch[p] -
                nf_veg.gru_frootn_xfer_to_atm_patch[p] * dt
            ns_veg.livestemn_xfer_patch[p] = ns_veg.livestemn_xfer_patch[p] -
                nf_veg.gru_livestemn_xfer_to_atm_patch[p] * dt
            ns_veg.deadstemn_xfer_patch[p] = ns_veg.deadstemn_xfer_patch[p] -
                nf_veg.gru_deadstemn_xfer_to_atm_patch[p] * dt
            ns_veg.livecrootn_xfer_patch[p] = ns_veg.livecrootn_xfer_patch[p] -
                nf_veg.gru_livecrootn_xfer_to_atm_patch[p] * dt
            ns_veg.deadcrootn_xfer_patch[p] = ns_veg.deadcrootn_xfer_patch[p] -
                nf_veg.gru_deadcrootn_xfer_to_atm_patch[p] * dt
        else
            # Matrix version: handled in dynGrossUnrepMod::CNGrossUnrep
        end
    end

    return nothing
end
