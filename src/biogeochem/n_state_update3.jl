# ==========================================================================
# Ported from: src/biogeochem/CNNStateUpdate3Mod.F90
# Nitrogen state variable update, mortality fluxes (fire).
# Also, sminn leaching flux.
#
# When the matrix solution is being used (use_matrixcn and use_soil_matrixcn)
# only some state updates are done here, the other state updates happen
# after the matrix is solved in VegMatrix and SoilMatrix.
#
# Public functions:
#   n_state_update3!          — Fire N state update
#   n_state_update_leaching!  — Sminn leaching N state update
# ==========================================================================

# ---------------------------------------------------------------------------
# n_state_update_leaching! — Sminn leaching
# ---------------------------------------------------------------------------

"""
    n_state_update_leaching!(ns_soil, nf_soil;
        mask_soilc, bounds_col, nlevdecomp,
        use_nitrif_denitrif, dt)

On the radiation time step, update all prognostic nitrogen state variables
affected by the Sminn leaching flux.

Ported from `NStateUpdateLeaching` in `CNNStateUpdate3Mod.F90`.

Note: This code was separated from gap mortality fluxes to make it
compatible with FATES.
"""
function n_state_update_leaching!(ns_soil::SoilBiogeochemNitrogenStateData,
                                   nf_soil::SoilBiogeochemNitrogenFluxData;
                                   mask_soilc::BitVector,
                                   bounds_col::UnitRange{Int},
                                   nlevdecomp::Int,
                                   use_nitrif_denitrif::Bool = false,
                                   dt::Real)

    for j in 1:nlevdecomp
        for c in bounds_col
            mask_soilc[c] || continue

            if !use_nitrif_denitrif
                # mineral N loss due to leaching
                ns_soil.sminn_vr_col[c, j] = ns_soil.sminn_vr_col[c, j] -
                    nf_soil.sminn_leached_vr_col[c, j] * dt
            else
                # mineral N loss due to leaching and runoff
                ns_soil.smin_no3_vr_col[c, j] = max(
                    ns_soil.smin_no3_vr_col[c, j] -
                    (nf_soil.smin_no3_leached_vr_col[c, j] + nf_soil.smin_no3_runoff_vr_col[c, j]) * dt,
                    0.0)

                ns_soil.sminn_vr_col[c, j] = ns_soil.smin_no3_vr_col[c, j] + ns_soil.smin_nh4_vr_col[c, j]
            end
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# n_state_update3! — Fire N state update
# ---------------------------------------------------------------------------

"""
    n_state_update3!(ns_veg, nf_veg, ns_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, ndecomp_pools, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

On the radiation time step, update all prognostic nitrogen state variables
affected by fire fluxes.

Ported from `NStateUpdate3` in `CNNStateUpdate3Mod.F90`.

Note: The matrix-CN code paths (`use_soil_matrixcn`) for soil matrix
solutions are omitted since they require the full matrix infrastructure
(matrix_Ninput). The matrix-CN veg path (`use_matrixcn`) is a no-op in
the Fortran original (comment only), so it is omitted here as well.
"""
function n_state_update3!(ns_veg::CNVegNitrogenStateData,
                           nf_veg::CNVegNitrogenFluxData,
                           ns_soil::SoilBiogeochemNitrogenStateData;
                           mask_soilc::BitVector,
                           mask_soilp::BitVector,
                           bounds_col::UnitRange{Int},
                           bounds_patch::UnitRange{Int},
                           nlevdecomp::Int,
                           ndecomp_pools::Int,
                           i_litr_min::Int,
                           i_litr_max::Int,
                           i_cwd::Int,
                           use_matrixcn::Bool = false,
                           use_soil_matrixcn::Bool = false,
                           dt::Real)

    # --- Column loop: fire mortality fluxes to soil decomposition pools ---
    for j in 1:nlevdecomp
        for c in bounds_col
            mask_soilc[c] || continue

            if !use_soil_matrixcn
                # patch-level wood to column-level CWD (uncombusted wood)
                ns_soil.decomp_npools_vr_col[c, j, i_cwd] =
                    ns_soil.decomp_npools_vr_col[c, j, i_cwd] +
                    nf_veg.fire_mortality_n_to_cwdn_col[c, j] * dt

                # patch-level wood to column-level litter (uncombusted wood)
                for k in i_litr_min:i_litr_max
                    ns_soil.decomp_npools_vr_col[c, j, k] =
                        ns_soil.decomp_npools_vr_col[c, j, k] +
                        nf_veg.m_n_to_litr_fire_col[c, j, k] * dt
                end
            else
                # Matrix path: matrix_Ninput setup omitted (requires full matrix infrastructure)
            end
        end
    end

    # --- Column loop: litter and CWD losses to fire ---
    if !use_soil_matrixcn
        for l in 1:ndecomp_pools
            for j in 1:nlevdecomp
                for c in bounds_col
                    mask_soilc[c] || continue
                    ns_soil.decomp_npools_vr_col[c, j, l] =
                        ns_soil.decomp_npools_vr_col[c, j, l] -
                        nf_veg.m_decomp_npools_to_fire_vr_col[c, j, l] * dt
                end
            end
        end
    end

    # --- Patch loop: vegetation N state updates from fire ---
    for p in bounds_patch
        mask_soilp[p] || continue

        if !use_matrixcn
            # from fire displayed pools
            ns_veg.leafn_patch[p] = ns_veg.leafn_patch[p] -
                nf_veg.m_leafn_to_fire_patch[p] * dt
            ns_veg.frootn_patch[p] = ns_veg.frootn_patch[p] -
                nf_veg.m_frootn_to_fire_patch[p] * dt
            ns_veg.livestemn_patch[p] = ns_veg.livestemn_patch[p] -
                nf_veg.m_livestemn_to_fire_patch[p] * dt
            ns_veg.deadstemn_patch[p] = ns_veg.deadstemn_patch[p] -
                nf_veg.m_deadstemn_to_fire_patch[p] * dt
            ns_veg.livecrootn_patch[p] = ns_veg.livecrootn_patch[p] -
                nf_veg.m_livecrootn_to_fire_patch[p] * dt
            ns_veg.deadcrootn_patch[p] = ns_veg.deadcrootn_patch[p] -
                nf_veg.m_deadcrootn_to_fire_patch[p] * dt

            ns_veg.leafn_patch[p] = ns_veg.leafn_patch[p] -
                nf_veg.m_leafn_to_litter_fire_patch[p] * dt
            ns_veg.frootn_patch[p] = ns_veg.frootn_patch[p] -
                nf_veg.m_frootn_to_litter_fire_patch[p] * dt
            ns_veg.livestemn_patch[p] = ns_veg.livestemn_patch[p] -
                nf_veg.m_livestemn_to_litter_fire_patch[p] * dt -
                nf_veg.m_livestemn_to_deadstemn_fire_patch[p] * dt
            ns_veg.deadstemn_patch[p] = ns_veg.deadstemn_patch[p] -
                nf_veg.m_deadstemn_to_litter_fire_patch[p] * dt +
                nf_veg.m_livestemn_to_deadstemn_fire_patch[p] * dt
            ns_veg.livecrootn_patch[p] = ns_veg.livecrootn_patch[p] -
                nf_veg.m_livecrootn_to_litter_fire_patch[p] * dt -
                nf_veg.m_livecrootn_to_deadcrootn_fire_patch[p] * dt
            ns_veg.deadcrootn_patch[p] = ns_veg.deadcrootn_patch[p] -
                nf_veg.m_deadcrootn_to_litter_fire_patch[p] * dt +
                nf_veg.m_livecrootn_to_deadcrootn_fire_patch[p] * dt

            # storage pools
            ns_veg.leafn_storage_patch[p] = ns_veg.leafn_storage_patch[p] -
                nf_veg.m_leafn_storage_to_fire_patch[p] * dt
            ns_veg.frootn_storage_patch[p] = ns_veg.frootn_storage_patch[p] -
                nf_veg.m_frootn_storage_to_fire_patch[p] * dt
            ns_veg.livestemn_storage_patch[p] = ns_veg.livestemn_storage_patch[p] -
                nf_veg.m_livestemn_storage_to_fire_patch[p] * dt
            ns_veg.deadstemn_storage_patch[p] = ns_veg.deadstemn_storage_patch[p] -
                nf_veg.m_deadstemn_storage_to_fire_patch[p] * dt
            ns_veg.livecrootn_storage_patch[p] = ns_veg.livecrootn_storage_patch[p] -
                nf_veg.m_livecrootn_storage_to_fire_patch[p] * dt
            ns_veg.deadcrootn_storage_patch[p] = ns_veg.deadcrootn_storage_patch[p] -
                nf_veg.m_deadcrootn_storage_to_fire_patch[p] * dt

            ns_veg.leafn_storage_patch[p] = ns_veg.leafn_storage_patch[p] -
                nf_veg.m_leafn_storage_to_litter_fire_patch[p] * dt
            ns_veg.frootn_storage_patch[p] = ns_veg.frootn_storage_patch[p] -
                nf_veg.m_frootn_storage_to_litter_fire_patch[p] * dt
            ns_veg.livestemn_storage_patch[p] = ns_veg.livestemn_storage_patch[p] -
                nf_veg.m_livestemn_storage_to_litter_fire_patch[p] * dt
            ns_veg.deadstemn_storage_patch[p] = ns_veg.deadstemn_storage_patch[p] -
                nf_veg.m_deadstemn_storage_to_litter_fire_patch[p] * dt
            ns_veg.livecrootn_storage_patch[p] = ns_veg.livecrootn_storage_patch[p] -
                nf_veg.m_livecrootn_storage_to_litter_fire_patch[p] * dt
            ns_veg.deadcrootn_storage_patch[p] = ns_veg.deadcrootn_storage_patch[p] -
                nf_veg.m_deadcrootn_storage_to_litter_fire_patch[p] * dt

            # transfer pools
            ns_veg.leafn_xfer_patch[p] = ns_veg.leafn_xfer_patch[p] -
                nf_veg.m_leafn_xfer_to_fire_patch[p] * dt
            ns_veg.frootn_xfer_patch[p] = ns_veg.frootn_xfer_patch[p] -
                nf_veg.m_frootn_xfer_to_fire_patch[p] * dt
            ns_veg.livestemn_xfer_patch[p] = ns_veg.livestemn_xfer_patch[p] -
                nf_veg.m_livestemn_xfer_to_fire_patch[p] * dt
            ns_veg.deadstemn_xfer_patch[p] = ns_veg.deadstemn_xfer_patch[p] -
                nf_veg.m_deadstemn_xfer_to_fire_patch[p] * dt
            ns_veg.livecrootn_xfer_patch[p] = ns_veg.livecrootn_xfer_patch[p] -
                nf_veg.m_livecrootn_xfer_to_fire_patch[p] * dt
            ns_veg.deadcrootn_xfer_patch[p] = ns_veg.deadcrootn_xfer_patch[p] -
                nf_veg.m_deadcrootn_xfer_to_fire_patch[p] * dt

            ns_veg.leafn_xfer_patch[p] = ns_veg.leafn_xfer_patch[p] -
                nf_veg.m_leafn_xfer_to_litter_fire_patch[p] * dt
            ns_veg.frootn_xfer_patch[p] = ns_veg.frootn_xfer_patch[p] -
                nf_veg.m_frootn_xfer_to_litter_fire_patch[p] * dt
            ns_veg.livestemn_xfer_patch[p] = ns_veg.livestemn_xfer_patch[p] -
                nf_veg.m_livestemn_xfer_to_litter_fire_patch[p] * dt
            ns_veg.deadstemn_xfer_patch[p] = ns_veg.deadstemn_xfer_patch[p] -
                nf_veg.m_deadstemn_xfer_to_litter_fire_patch[p] * dt
            ns_veg.livecrootn_xfer_patch[p] = ns_veg.livecrootn_xfer_patch[p] -
                nf_veg.m_livecrootn_xfer_to_litter_fire_patch[p] * dt
            ns_veg.deadcrootn_xfer_patch[p] = ns_veg.deadcrootn_xfer_patch[p] -
                nf_veg.m_deadcrootn_xfer_to_litter_fire_patch[p] * dt

            # retranslocated N pool
            ns_veg.retransn_patch[p] = ns_veg.retransn_patch[p] -
                nf_veg.m_retransn_to_fire_patch[p] * dt
            ns_veg.retransn_patch[p] = ns_veg.retransn_patch[p] -
                nf_veg.m_retransn_to_litter_fire_patch[p] * dt
        else
            # NOTE: The equivalent changes for matrix code are in CNFireBase and CNFireLi2014 codes
        end
    end

    return nothing
end
