# ==========================================================================
# Ported from: src/biogeochem/CNNStateUpdate1Mod.F90
# Nitrogen state variable update, non-mortality fluxes.
#
# When the matrix solution is being used (use_matrixcn and use_soil_matrixcn)
# only some state updates are done here, the other state updates happen
# after the matrix is solved in VegMatrix and SoilMatrix.
#
# Public functions:
#   n_state_update_dyn_patch! — Update N states from dyn_cnbal_patch fluxes
#   n_state_update1!         — Full prognostic N state update (non-mortality, non-fire)
# ==========================================================================

# ---------------------------------------------------------------------------
# n_state_update_dyn_patch! — Dynamic patch nitrogen state update
# ---------------------------------------------------------------------------

"""
    n_state_update_dyn_patch!(ns_veg, nf_veg, ns_soil;
        mask_soilc_with_inactive, bounds_col, bounds_grc,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd, dt)

Update nitrogen states based on fluxes from dyn_cnbal_patch.

Ported from `NStateUpdateDynPatch` in `CNNStateUpdate1Mod.F90`.
"""
function n_state_update_dyn_patch!(ns_veg::CNVegNitrogenStateData,
                                    nf_veg::CNVegNitrogenFluxData,
                                    ns_soil::SoilBiogeochemNitrogenStateData;
                                    mask_soilc_with_inactive::BitVector,
                                    bounds_col::UnitRange{Int},
                                    bounds_grc::UnitRange{Int},
                                    nlevdecomp::Int,
                                    i_litr_min::Int,
                                    i_litr_max::Int,
                                    i_cwd::Int,
                                    dt::Real)

    for j in 1:nlevdecomp
        for c in bounds_col
            mask_soilc_with_inactive[c] || continue
            for i in i_litr_min:i_litr_max
                ns_soil.decomp_npools_vr_col[c, j, i] =
                    ns_soil.decomp_npools_vr_col[c, j, i] +
                    nf_veg.dwt_frootn_to_litr_n_col[c, j, i] * dt
            end
            ns_soil.decomp_npools_vr_col[c, j, i_cwd] =
                ns_soil.decomp_npools_vr_col[c, j, i_cwd] +
                (nf_veg.dwt_livecrootn_to_cwdn_col[c, j] +
                 nf_veg.dwt_deadcrootn_to_cwdn_col[c, j]) * dt
        end
    end

    for g in bounds_grc
        ns_veg.seedn_grc[g] = ns_veg.seedn_grc[g] - nf_veg.dwt_seedn_to_leaf_grc[g] * dt
        ns_veg.seedn_grc[g] = ns_veg.seedn_grc[g] - nf_veg.dwt_seedn_to_deadstem_grc[g] * dt
    end

    return nothing
end

# ---------------------------------------------------------------------------
# n_state_update1! — Main prognostic N state update
# ---------------------------------------------------------------------------

"""
    n_state_update1!(ns_veg, nf_veg, nf_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        ivt, woody, col_is_fates,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        npcropmin, nrepr, repr_grain_min, repr_grain_max,
        repr_structure_min, repr_structure_max,
        use_matrixcn, use_soil_matrixcn, use_fun, dt)

On the radiation time step, update all prognostic nitrogen state
variables (except for gap-phase mortality and fire fluxes).

Ported from `NStateUpdate1` in `CNNStateUpdate1Mod.F90`.

Note: The FATES litter flux callback (`clm_fates%UpdateNLitterfluxes`) is
skipped in this port — FATES columns are handled by zeroing
`decomp_npools_sourcesink_col` via `col_is_fates`. The matrix-CN code paths
(`use_matrixcn`, `use_soil_matrixcn`) for soil/veg matrix solutions are
ported as-is (non-matrix path only; matrix input paths are no-ops since
they depend on the full matrix infrastructure).
"""
function n_state_update1!(ns_veg::CNVegNitrogenStateData,
                           nf_veg::CNVegNitrogenFluxData,
                           nf_soil::SoilBiogeochemNitrogenFluxData;
                           mask_soilc::BitVector,
                           mask_soilp::BitVector,
                           bounds_col::UnitRange{Int},
                           bounds_patch::UnitRange{Int},
                           ivt::Vector{Int},
                           woody::Vector{<:Real},
                           col_is_fates::Vector{Bool},
                           nlevdecomp::Int,
                           i_litr_min::Int,
                           i_litr_max::Int,
                           i_cwd::Int,
                           npcropmin::Int = NPCROPMIN,
                           nrepr::Int = NREPR,
                           repr_grain_min::Int = 1,
                           repr_grain_max::Int = 1,
                           repr_structure_min::Int = 1,
                           repr_structure_max::Int = 0,
                           use_matrixcn::Bool = false,
                           use_soil_matrixcn::Bool = false,
                           use_fun::Bool = false,
                           dt::Real)

    # --- Column loop: soil decomposition input fluxes ---
    for c in bounds_col
        mask_soilc[c] || continue

        if col_is_fates[c]
            # FATES columns: litter fluxes handled externally (skip)
            # In full CLM this calls clm_fates%UpdateNLitterfluxes
        else
            for j in 1:nlevdecomp
                if !use_soil_matrixcn
                    for i in i_litr_min:i_litr_max
                        nf_soil.decomp_npools_sourcesink_col[c, j, i] =
                            nf_veg.phenology_n_to_litr_n_col[c, j, i] * dt
                    end
                    # NOTE(wjs, 2017-01-02) This used to be set to a non-zero value, but the
                    # terms have been moved to NStateUpdateDynPatch. Explicitly zeroed here.
                    nf_soil.decomp_npools_sourcesink_col[c, j, i_cwd] = 0.0
                end
                # Note: matrix_Ninput path omitted (requires full matrix infrastructure)
            end
        end
    end

    # --- Patch loop: vegetation N state updates ---
    for p in bounds_patch
        mask_soilp[p] || continue
        ivt[p] >= 1 || continue  # skip bare ground (Fortran PFT index 0)

        # === Phenology: transfer growth fluxes ===
        if !use_matrixcn
            ns_veg.leafn_patch[p]       = ns_veg.leafn_patch[p]       + nf_veg.leafn_xfer_to_leafn_patch[p] * dt
            ns_veg.leafn_xfer_patch[p]  = ns_veg.leafn_xfer_patch[p]  - nf_veg.leafn_xfer_to_leafn_patch[p] * dt
            ns_veg.frootn_patch[p]      = ns_veg.frootn_patch[p]      + nf_veg.frootn_xfer_to_frootn_patch[p] * dt
            ns_veg.frootn_xfer_patch[p] = ns_veg.frootn_xfer_patch[p] - nf_veg.frootn_xfer_to_frootn_patch[p] * dt

            if woody[ivt[p]] == 1.0
                ns_veg.livestemn_patch[p]       = ns_veg.livestemn_patch[p]       + nf_veg.livestemn_xfer_to_livestemn_patch[p] * dt
                ns_veg.livestemn_xfer_patch[p]  = ns_veg.livestemn_xfer_patch[p]  - nf_veg.livestemn_xfer_to_livestemn_patch[p] * dt
                ns_veg.deadstemn_patch[p]       = ns_veg.deadstemn_patch[p]       + nf_veg.deadstemn_xfer_to_deadstemn_patch[p] * dt
                ns_veg.deadstemn_xfer_patch[p]  = ns_veg.deadstemn_xfer_patch[p]  - nf_veg.deadstemn_xfer_to_deadstemn_patch[p] * dt
                ns_veg.livecrootn_patch[p]      = ns_veg.livecrootn_patch[p]      + nf_veg.livecrootn_xfer_to_livecrootn_patch[p] * dt
                ns_veg.livecrootn_xfer_patch[p] = ns_veg.livecrootn_xfer_patch[p] - nf_veg.livecrootn_xfer_to_livecrootn_patch[p] * dt
                ns_veg.deadcrootn_patch[p]      = ns_veg.deadcrootn_patch[p]      + nf_veg.deadcrootn_xfer_to_deadcrootn_patch[p] * dt
                ns_veg.deadcrootn_xfer_patch[p] = ns_veg.deadcrootn_xfer_patch[p] - nf_veg.deadcrootn_xfer_to_deadcrootn_patch[p] * dt
            end

            if ivt[p] >= npcropmin  # skip 2 generic crops
                # lines here for consistency; the transfer terms are zero
                ns_veg.livestemn_patch[p]      = ns_veg.livestemn_patch[p]      + nf_veg.livestemn_xfer_to_livestemn_patch[p] * dt
                ns_veg.livestemn_xfer_patch[p] = ns_veg.livestemn_xfer_patch[p] - nf_veg.livestemn_xfer_to_livestemn_patch[p] * dt
                for k in 1:nrepr
                    ns_veg.reproductiven_patch[p, k] = ns_veg.reproductiven_patch[p, k] +
                        nf_veg.reproductiven_xfer_to_reproductiven_patch[p, k] * dt
                    ns_veg.reproductiven_xfer_patch[p, k] = ns_veg.reproductiven_xfer_patch[p, k] -
                        nf_veg.reproductiven_xfer_to_reproductiven_patch[p, k] * dt
                end
            end

            # phenology: litterfall and retranslocation fluxes
            ns_veg.leafn_patch[p]    = ns_veg.leafn_patch[p]    - nf_veg.leafn_to_litter_patch[p] * dt
            ns_veg.frootn_patch[p]   = ns_veg.frootn_patch[p]   - nf_veg.frootn_to_litter_patch[p] * dt
            ns_veg.leafn_patch[p]    = ns_veg.leafn_patch[p]    - nf_veg.leafn_to_retransn_patch[p] * dt
            ns_veg.retransn_patch[p] = ns_veg.retransn_patch[p] + nf_veg.leafn_to_retransn_patch[p] * dt
        end  # !use_matrixcn

        # live wood turnover and retranslocation fluxes
        if woody[ivt[p]] == 1.0
            if !use_matrixcn
                ns_veg.livestemn_patch[p]  = ns_veg.livestemn_patch[p]  - nf_veg.livestemn_to_deadstemn_patch[p] * dt
                ns_veg.deadstemn_patch[p]  = ns_veg.deadstemn_patch[p]  + nf_veg.livestemn_to_deadstemn_patch[p] * dt
                ns_veg.livestemn_patch[p]  = ns_veg.livestemn_patch[p]  - nf_veg.livestemn_to_retransn_patch[p] * dt
                ns_veg.retransn_patch[p]   = ns_veg.retransn_patch[p]   + nf_veg.livestemn_to_retransn_patch[p] * dt
                ns_veg.livecrootn_patch[p] = ns_veg.livecrootn_patch[p] - nf_veg.livecrootn_to_deadcrootn_patch[p] * dt
                ns_veg.deadcrootn_patch[p] = ns_veg.deadcrootn_patch[p] + nf_veg.livecrootn_to_deadcrootn_patch[p] * dt
                ns_veg.livecrootn_patch[p] = ns_veg.livecrootn_patch[p] - nf_veg.livecrootn_to_retransn_patch[p] * dt
                ns_veg.retransn_patch[p]   = ns_veg.retransn_patch[p]   + nf_veg.livecrootn_to_retransn_patch[p] * dt
                # WW change logic so livestem_retrans goes to npool (via free_retrans flux)
                if use_fun
                    nf_veg.free_retransn_to_npool_patch[p] = nf_veg.free_retransn_to_npool_patch[p] + nf_veg.livestemn_to_retransn_patch[p]
                    nf_veg.free_retransn_to_npool_patch[p] = nf_veg.free_retransn_to_npool_patch[p] + nf_veg.livecrootn_to_retransn_patch[p]
                end
            end  # !use_matrixcn
        end

        if ivt[p] >= npcropmin  # Beth adds retrans from froot
            if !use_matrixcn
                ns_veg.frootn_patch[p]    = ns_veg.frootn_patch[p]    - nf_veg.frootn_to_retransn_patch[p] * dt
                ns_veg.retransn_patch[p]  = ns_veg.retransn_patch[p]  + nf_veg.frootn_to_retransn_patch[p] * dt
                ns_veg.livestemn_patch[p] = ns_veg.livestemn_patch[p] - nf_veg.livestemn_to_litter_patch[p] * dt
                ns_veg.livestemn_patch[p] = ns_veg.livestemn_patch[p] -
                    (nf_veg.livestemn_to_biofueln_patch[p] + nf_veg.livestemn_to_removedresiduen_patch[p]) * dt
                ns_veg.leafn_patch[p]     = ns_veg.leafn_patch[p] -
                    (nf_veg.leafn_to_biofueln_patch[p] + nf_veg.leafn_to_removedresiduen_patch[p]) * dt
                ns_veg.livestemn_patch[p] = ns_veg.livestemn_patch[p] - nf_veg.livestemn_to_retransn_patch[p] * dt
                ns_veg.retransn_patch[p]  = ns_veg.retransn_patch[p]  + nf_veg.livestemn_to_retransn_patch[p] * dt
                for k in repr_grain_min:repr_grain_max
                    ns_veg.reproductiven_patch[p, k] = ns_veg.reproductiven_patch[p, k] -
                        (nf_veg.repr_grainn_to_food_patch[p, k] + nf_veg.repr_grainn_to_seed_patch[p, k]) * dt
                end
                for k in repr_structure_min:repr_structure_max
                    ns_veg.reproductiven_patch[p, k] = ns_veg.reproductiven_patch[p, k] -
                        (nf_veg.repr_structuren_to_cropprod_patch[p, k] + nf_veg.repr_structuren_to_litter_patch[p, k]) * dt
                end
            end  # !use_matrixcn
            ns_veg.cropseedn_deficit_patch[p] = ns_veg.cropseedn_deficit_patch[p] -
                nf_veg.crop_seedn_to_leaf_patch[p] * dt
            for k in repr_grain_min:repr_grain_max
                ns_veg.cropseedn_deficit_patch[p] = ns_veg.cropseedn_deficit_patch[p] +
                    nf_veg.repr_grainn_to_seed_patch[p, k] * dt
            end
        end

        # uptake from soil mineral N pool
        ns_veg.npool_patch[p] = ns_veg.npool_patch[p] + nf_veg.sminn_to_npool_patch[p] * dt

        # deployment from retranslocation pool
        ns_veg.npool_patch[p] = ns_veg.npool_patch[p] + nf_veg.retransn_to_npool_patch[p] * dt

        ns_veg.npool_patch[p] = ns_veg.npool_patch[p] + nf_veg.free_retransn_to_npool_patch[p] * dt

        # allocation fluxes
        ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_leafn_patch[p] * dt
        ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_leafn_storage_patch[p] * dt
        ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_frootn_patch[p] * dt
        ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_frootn_storage_patch[p] * dt

        if !use_matrixcn
            ns_veg.retransn_patch[p]       = ns_veg.retransn_patch[p]       - nf_veg.retransn_to_npool_patch[p] * dt
            ns_veg.retransn_patch[p]       = ns_veg.retransn_patch[p]       - nf_veg.free_retransn_to_npool_patch[p] * dt
            ns_veg.leafn_patch[p]          = ns_veg.leafn_patch[p]          + nf_veg.npool_to_leafn_patch[p] * dt
            ns_veg.leafn_storage_patch[p]  = ns_veg.leafn_storage_patch[p]  + nf_veg.npool_to_leafn_storage_patch[p] * dt
            ns_veg.frootn_patch[p]         = ns_veg.frootn_patch[p]         + nf_veg.npool_to_frootn_patch[p] * dt
            ns_veg.frootn_storage_patch[p] = ns_veg.frootn_storage_patch[p] + nf_veg.npool_to_frootn_storage_patch[p] * dt
        end

        if woody[ivt[p]] == 1.0
            ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_livestemn_patch[p] * dt
            ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_livestemn_storage_patch[p] * dt
            ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_deadstemn_patch[p] * dt
            ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_deadstemn_storage_patch[p] * dt
            ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_livecrootn_patch[p] * dt
            ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_livecrootn_storage_patch[p] * dt
            ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_deadcrootn_patch[p] * dt
            ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_deadcrootn_storage_patch[p] * dt

            if !use_matrixcn
                ns_veg.livestemn_patch[p]          = ns_veg.livestemn_patch[p]          + nf_veg.npool_to_livestemn_patch[p] * dt
                ns_veg.livestemn_storage_patch[p]  = ns_veg.livestemn_storage_patch[p]  + nf_veg.npool_to_livestemn_storage_patch[p] * dt
                ns_veg.deadstemn_patch[p]          = ns_veg.deadstemn_patch[p]          + nf_veg.npool_to_deadstemn_patch[p] * dt
                ns_veg.deadstemn_storage_patch[p]  = ns_veg.deadstemn_storage_patch[p]  + nf_veg.npool_to_deadstemn_storage_patch[p] * dt
                ns_veg.livecrootn_patch[p]         = ns_veg.livecrootn_patch[p]         + nf_veg.npool_to_livecrootn_patch[p] * dt
                ns_veg.livecrootn_storage_patch[p] = ns_veg.livecrootn_storage_patch[p] + nf_veg.npool_to_livecrootn_storage_patch[p] * dt
                ns_veg.deadcrootn_patch[p]         = ns_veg.deadcrootn_patch[p]         + nf_veg.npool_to_deadcrootn_patch[p] * dt
                ns_veg.deadcrootn_storage_patch[p] = ns_veg.deadcrootn_storage_patch[p] + nf_veg.npool_to_deadcrootn_storage_patch[p] * dt
            end
        end

        if ivt[p] >= npcropmin  # skip 2 generic crops
            ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_livestemn_patch[p] * dt
            ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_livestemn_storage_patch[p] * dt
            for k in 1:nrepr
                ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_reproductiven_patch[p, k] * dt
                ns_veg.npool_patch[p] = ns_veg.npool_patch[p] - nf_veg.npool_to_reproductiven_storage_patch[p, k] * dt
            end

            if !use_matrixcn
                ns_veg.livestemn_patch[p]         = ns_veg.livestemn_patch[p]         + nf_veg.npool_to_livestemn_patch[p] * dt
                ns_veg.livestemn_storage_patch[p] = ns_veg.livestemn_storage_patch[p] + nf_veg.npool_to_livestemn_storage_patch[p] * dt
                for k in 1:nrepr
                    ns_veg.reproductiven_patch[p, k] = ns_veg.reproductiven_patch[p, k] +
                        nf_veg.npool_to_reproductiven_patch[p, k] * dt
                    ns_veg.reproductiven_storage_patch[p, k] = ns_veg.reproductiven_storage_patch[p, k] +
                        nf_veg.npool_to_reproductiven_storage_patch[p, k] * dt
                end
            end
        end

        # move storage pools into transfer pools
        if !use_matrixcn
            ns_veg.leafn_storage_patch[p]  = ns_veg.leafn_storage_patch[p]  - nf_veg.leafn_storage_to_xfer_patch[p] * dt
            ns_veg.leafn_xfer_patch[p]     = ns_veg.leafn_xfer_patch[p]     + nf_veg.leafn_storage_to_xfer_patch[p] * dt
            ns_veg.frootn_storage_patch[p] = ns_veg.frootn_storage_patch[p] - nf_veg.frootn_storage_to_xfer_patch[p] * dt
            ns_veg.frootn_xfer_patch[p]    = ns_veg.frootn_xfer_patch[p]    + nf_veg.frootn_storage_to_xfer_patch[p] * dt

            if woody[ivt[p]] == 1.0
                ns_veg.livestemn_storage_patch[p]  = ns_veg.livestemn_storage_patch[p]  - nf_veg.livestemn_storage_to_xfer_patch[p] * dt
                ns_veg.livestemn_xfer_patch[p]     = ns_veg.livestemn_xfer_patch[p]     + nf_veg.livestemn_storage_to_xfer_patch[p] * dt
                ns_veg.deadstemn_storage_patch[p]  = ns_veg.deadstemn_storage_patch[p]  - nf_veg.deadstemn_storage_to_xfer_patch[p] * dt
                ns_veg.deadstemn_xfer_patch[p]     = ns_veg.deadstemn_xfer_patch[p]     + nf_veg.deadstemn_storage_to_xfer_patch[p] * dt
                ns_veg.livecrootn_storage_patch[p] = ns_veg.livecrootn_storage_patch[p] - nf_veg.livecrootn_storage_to_xfer_patch[p] * dt
                ns_veg.livecrootn_xfer_patch[p]    = ns_veg.livecrootn_xfer_patch[p]    + nf_veg.livecrootn_storage_to_xfer_patch[p] * dt
                ns_veg.deadcrootn_storage_patch[p] = ns_veg.deadcrootn_storage_patch[p] - nf_veg.deadcrootn_storage_to_xfer_patch[p] * dt
                ns_veg.deadcrootn_xfer_patch[p]    = ns_veg.deadcrootn_xfer_patch[p]    + nf_veg.deadcrootn_storage_to_xfer_patch[p] * dt
            end
        end  # !use_matrixcn

        if ivt[p] >= npcropmin  # skip 2 generic crops
            # lines here for consistency; the transfer terms are zero
            if !use_matrixcn
                ns_veg.livestemn_storage_patch[p] = ns_veg.livestemn_storage_patch[p] - nf_veg.livestemn_storage_to_xfer_patch[p] * dt
                ns_veg.livestemn_xfer_patch[p]    = ns_veg.livestemn_xfer_patch[p]    + nf_veg.livestemn_storage_to_xfer_patch[p] * dt
                for k in 1:nrepr
                    ns_veg.reproductiven_storage_patch[p, k] = ns_veg.reproductiven_storage_patch[p, k] -
                        nf_veg.reproductiven_storage_to_xfer_patch[p, k] * dt
                    ns_veg.reproductiven_xfer_patch[p, k] = ns_veg.reproductiven_xfer_patch[p, k] +
                        nf_veg.reproductiven_storage_to_xfer_patch[p, k] * dt
                end
            end  # !use_matrixcn
        end

    end  # end of patch loop

    return nothing
end
