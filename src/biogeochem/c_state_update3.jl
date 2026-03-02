# ==========================================================================
# Ported from: src/biogeochem/CNCStateUpdate3Mod.F90
# Carbon state variable update, fire fluxes.
#
# When the matrix solution is being used (use_matrixcn and use_soil_matrixcn)
# only some state updates are done here; other state updates happen after the
# matrix is solved in VegMatrix and SoilMatrix.
#
# Public functions:
#   c_state_update3! — Fire C state update
# ==========================================================================

# ---------------------------------------------------------------------------
# c_state_update3! — Fire carbon state update
# ---------------------------------------------------------------------------

"""
    c_state_update3!(cs_veg, cf_veg, cs_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, ndecomp_pools, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

On the radiation time step, update all prognostic carbon state variables
affected by fire fluxes.

Ported from `CStateUpdate3` in `CNCStateUpdate3Mod.F90`.

Note: The matrix-CN code paths (`use_soil_matrixcn`) for soil matrix
solutions are omitted since they require the full matrix infrastructure
(matrix_Cinput). The matrix-CN veg path (`use_matrixcn`) is a no-op in
the Fortran original (comment only), so it is omitted here as well.
"""
function c_state_update3!(cs_veg::CNVegCarbonStateData,
                           cf_veg::CNVegCarbonFluxData,
                           cs_soil::SoilBiogeochemCarbonStateData;
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
                           dt::Float64)

    # --- Column loop: fire mortality fluxes to soil decomposition pools ---
    for j in 1:nlevdecomp
        for c in bounds_col
            mask_soilc[c] || continue

            if !use_soil_matrixcn
                # patch-level wood to column-level CWD (uncombusted wood)
                cs_soil.decomp_cpools_vr_col[c, j, i_cwd] =
                    cs_soil.decomp_cpools_vr_col[c, j, i_cwd] +
                    cf_veg.fire_mortality_c_to_cwdc_col[c, j] * dt

                # patch-level wood to column-level litter (uncombusted wood)
                for i in i_litr_min:i_litr_max
                    cs_soil.decomp_cpools_vr_col[c, j, i] =
                        cs_soil.decomp_cpools_vr_col[c, j, i] +
                        cf_veg.m_c_to_litr_fire_col[c, j, i] * dt
                end
            else
                # Matrix path: matrix_Cinput setup omitted (requires full matrix infrastructure)
            end
        end
    end

    # --- Column loop: litter and CWD losses to fire ---
    if !use_soil_matrixcn
        for l in 1:ndecomp_pools
            for j in 1:nlevdecomp
                for c in bounds_col
                    mask_soilc[c] || continue
                    cs_soil.decomp_cpools_vr_col[c, j, l] =
                        cs_soil.decomp_cpools_vr_col[c, j, l] -
                        cf_veg.m_decomp_cpools_to_fire_vr_col[c, j, l] * dt
                end
            end
        end
    end

    # --- Patch loop: vegetation C state updates from fire ---
    for p in bounds_patch
        mask_soilp[p] || continue

        cs_veg.gresp_storage_patch[p] = cs_veg.gresp_storage_patch[p] -
            cf_veg.m_gresp_storage_to_fire_patch[p] * dt
        cs_veg.gresp_storage_patch[p] = cs_veg.gresp_storage_patch[p] -
            cf_veg.m_gresp_storage_to_litter_fire_patch[p] * dt
        cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] -
            cf_veg.m_gresp_xfer_to_fire_patch[p] * dt
        cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] -
            cf_veg.m_gresp_xfer_to_litter_fire_patch[p] * dt

        if !use_matrixcn
            # displayed pools
            cs_veg.leafc_patch[p] = cs_veg.leafc_patch[p] -
                cf_veg.m_leafc_to_fire_patch[p] * dt
            cs_veg.leafc_patch[p] = cs_veg.leafc_patch[p] -
                cf_veg.m_leafc_to_litter_fire_patch[p] * dt
            cs_veg.frootc_patch[p] = cs_veg.frootc_patch[p] -
                cf_veg.m_frootc_to_fire_patch[p] * dt
            cs_veg.frootc_patch[p] = cs_veg.frootc_patch[p] -
                cf_veg.m_frootc_to_litter_fire_patch[p] * dt
            cs_veg.livestemc_patch[p] = cs_veg.livestemc_patch[p] -
                cf_veg.m_livestemc_to_fire_patch[p] * dt
            cs_veg.livestemc_patch[p] = cs_veg.livestemc_patch[p] -
                cf_veg.m_livestemc_to_litter_fire_patch[p] * dt -
                cf_veg.m_livestemc_to_deadstemc_fire_patch[p] * dt
            cs_veg.deadstemc_patch[p] = cs_veg.deadstemc_patch[p] -
                cf_veg.m_deadstemc_to_fire_patch[p] * dt
            cs_veg.deadstemc_patch[p] = cs_veg.deadstemc_patch[p] -
                cf_veg.m_deadstemc_to_litter_fire_patch[p] * dt +
                cf_veg.m_livestemc_to_deadstemc_fire_patch[p] * dt
            cs_veg.livecrootc_patch[p] = cs_veg.livecrootc_patch[p] -
                cf_veg.m_livecrootc_to_fire_patch[p] * dt
            cs_veg.livecrootc_patch[p] = cs_veg.livecrootc_patch[p] -
                cf_veg.m_livecrootc_to_litter_fire_patch[p] * dt -
                cf_veg.m_livecrootc_to_deadcrootc_fire_patch[p] * dt
            cs_veg.deadcrootc_patch[p] = cs_veg.deadcrootc_patch[p] -
                cf_veg.m_deadcrootc_to_fire_patch[p] * dt
            cs_veg.deadcrootc_patch[p] = cs_veg.deadcrootc_patch[p] -
                cf_veg.m_deadcrootc_to_litter_fire_patch[p] * dt +
                cf_veg.m_livecrootc_to_deadcrootc_fire_patch[p] * dt

            # storage pools
            cs_veg.leafc_storage_patch[p] = cs_veg.leafc_storage_patch[p] -
                cf_veg.m_leafc_storage_to_fire_patch[p] * dt
            cs_veg.leafc_storage_patch[p] = cs_veg.leafc_storage_patch[p] -
                cf_veg.m_leafc_storage_to_litter_fire_patch[p] * dt
            cs_veg.frootc_storage_patch[p] = cs_veg.frootc_storage_patch[p] -
                cf_veg.m_frootc_storage_to_fire_patch[p] * dt
            cs_veg.frootc_storage_patch[p] = cs_veg.frootc_storage_patch[p] -
                cf_veg.m_frootc_storage_to_litter_fire_patch[p] * dt
            cs_veg.livestemc_storage_patch[p] = cs_veg.livestemc_storage_patch[p] -
                cf_veg.m_livestemc_storage_to_fire_patch[p] * dt
            cs_veg.livestemc_storage_patch[p] = cs_veg.livestemc_storage_patch[p] -
                cf_veg.m_livestemc_storage_to_litter_fire_patch[p] * dt
            cs_veg.deadstemc_storage_patch[p] = cs_veg.deadstemc_storage_patch[p] -
                cf_veg.m_deadstemc_storage_to_fire_patch[p] * dt
            cs_veg.deadstemc_storage_patch[p] = cs_veg.deadstemc_storage_patch[p] -
                cf_veg.m_deadstemc_storage_to_litter_fire_patch[p] * dt
            cs_veg.livecrootc_storage_patch[p] = cs_veg.livecrootc_storage_patch[p] -
                cf_veg.m_livecrootc_storage_to_fire_patch[p] * dt
            cs_veg.livecrootc_storage_patch[p] = cs_veg.livecrootc_storage_patch[p] -
                cf_veg.m_livecrootc_storage_to_litter_fire_patch[p] * dt
            cs_veg.deadcrootc_storage_patch[p] = cs_veg.deadcrootc_storage_patch[p] -
                cf_veg.m_deadcrootc_storage_to_fire_patch[p] * dt
            cs_veg.deadcrootc_storage_patch[p] = cs_veg.deadcrootc_storage_patch[p] -
                cf_veg.m_deadcrootc_storage_to_litter_fire_patch[p] * dt

            # transfer pools
            cs_veg.leafc_xfer_patch[p] = cs_veg.leafc_xfer_patch[p] -
                cf_veg.m_leafc_xfer_to_fire_patch[p] * dt
            cs_veg.leafc_xfer_patch[p] = cs_veg.leafc_xfer_patch[p] -
                cf_veg.m_leafc_xfer_to_litter_fire_patch[p] * dt
            cs_veg.frootc_xfer_patch[p] = cs_veg.frootc_xfer_patch[p] -
                cf_veg.m_frootc_xfer_to_fire_patch[p] * dt
            cs_veg.frootc_xfer_patch[p] = cs_veg.frootc_xfer_patch[p] -
                cf_veg.m_frootc_xfer_to_litter_fire_patch[p] * dt
            cs_veg.livestemc_xfer_patch[p] = cs_veg.livestemc_xfer_patch[p] -
                cf_veg.m_livestemc_xfer_to_fire_patch[p] * dt
            cs_veg.livestemc_xfer_patch[p] = cs_veg.livestemc_xfer_patch[p] -
                cf_veg.m_livestemc_xfer_to_litter_fire_patch[p] * dt
            cs_veg.deadstemc_xfer_patch[p] = cs_veg.deadstemc_xfer_patch[p] -
                cf_veg.m_deadstemc_xfer_to_fire_patch[p] * dt
            cs_veg.deadstemc_xfer_patch[p] = cs_veg.deadstemc_xfer_patch[p] -
                cf_veg.m_deadstemc_xfer_to_litter_fire_patch[p] * dt
            cs_veg.livecrootc_xfer_patch[p] = cs_veg.livecrootc_xfer_patch[p] -
                cf_veg.m_livecrootc_xfer_to_fire_patch[p] * dt
            cs_veg.livecrootc_xfer_patch[p] = cs_veg.livecrootc_xfer_patch[p] -
                cf_veg.m_livecrootc_xfer_to_litter_fire_patch[p] * dt
            cs_veg.deadcrootc_xfer_patch[p] = cs_veg.deadcrootc_xfer_patch[p] -
                cf_veg.m_deadcrootc_xfer_to_fire_patch[p] * dt
            cs_veg.deadcrootc_xfer_patch[p] = cs_veg.deadcrootc_xfer_patch[p] -
                cf_veg.m_deadcrootc_xfer_to_litter_fire_patch[p] * dt
        else
            # NOTE: The equivalent changes for matrix code are in CNFireBase and CNFireLi2014 codes
        end
    end

    return nothing
end
