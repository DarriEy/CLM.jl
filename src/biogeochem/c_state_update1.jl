# ==========================================================================
# Ported from: src/biogeochem/CNCStateUpdate1Mod.F90
# Carbon state variable update, non-mortality fluxes.
#
# Public functions:
#   c_state_update_dyn_patch! — Update C states from dyn_cnbal_patch fluxes
#   c_state_update0!         — Update cpool from photosynthesis on radiation timestep
#   c_state_update1!         — Full prognostic C state update (non-mortality, non-fire)
# ==========================================================================

# ---------------------------------------------------------------------------
# c_state_update_dyn_patch! — Dynamic patch carbon state update
# ---------------------------------------------------------------------------

"""
    c_state_update_dyn_patch!(cs_veg, cf_veg, cs_soil;
        mask_soilc_with_inactive, bounds_col, bounds_grc,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd, dt)

Update carbon states based on fluxes from dyn_cnbal_patch.
This routine is not called with FATES active.

Ported from `CStateUpdateDynPatch` in `CNCStateUpdate1Mod.F90`.
"""
function c_state_update_dyn_patch!(cs_veg::CNVegCarbonStateData,
                                    cf_veg::CNVegCarbonFluxData,
                                    cs_soil::SoilBiogeochemCarbonStateData;
                                    mask_soilc_with_inactive::BitVector,
                                    bounds_col::UnitRange{Int},
                                    bounds_grc::UnitRange{Int},
                                    nlevdecomp::Int,
                                    i_litr_min::Int,
                                    i_litr_max::Int,
                                    i_cwd::Int,
                                    dt::Float64)

    for j in 1:nlevdecomp
        for c in bounds_col
            mask_soilc_with_inactive[c] || continue
            for i in i_litr_min:i_litr_max
                cs_soil.decomp_cpools_vr_col[c, j, i] =
                    cs_soil.decomp_cpools_vr_col[c, j, i] +
                    cf_veg.dwt_frootc_to_litr_c_col[c, j, i] * dt
            end
            cs_soil.decomp_cpools_vr_col[c, j, i_cwd] =
                cs_soil.decomp_cpools_vr_col[c, j, i_cwd] +
                (cf_veg.dwt_livecrootc_to_cwdc_col[c, j] +
                 cf_veg.dwt_deadcrootc_to_cwdc_col[c, j]) * dt
        end
    end

    for g in bounds_grc
        cs_veg.seedc_grc[g] = cs_veg.seedc_grc[g] - cf_veg.dwt_seedc_to_leaf_grc[g] * dt
        cs_veg.seedc_grc[g] = cs_veg.seedc_grc[g] - cf_veg.dwt_seedc_to_deadstem_grc[g] * dt
    end

    return nothing
end

# ---------------------------------------------------------------------------
# c_state_update0! — Photosynthesis cpool update
# ---------------------------------------------------------------------------

"""
    c_state_update0!(cs_veg, cf_veg;
        mask_soilp, bounds_patch, dt)

On the radiation time step, update cpool carbon state from gross
photosynthesis fluxes.

Ported from `CStateUpdate0` in `CNCStateUpdate1Mod.F90`.
"""
function c_state_update0!(cs_veg::CNVegCarbonStateData,
                           cf_veg::CNVegCarbonFluxData;
                           mask_soilp::BitVector,
                           bounds_patch::UnitRange{Int},
                           dt::Float64)

    for p in bounds_patch
        mask_soilp[p] || continue
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] + cf_veg.psnsun_to_cpool_patch[p] * dt
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] + cf_veg.psnshade_to_cpool_patch[p] * dt
    end

    return nothing
end

# ---------------------------------------------------------------------------
# c_state_update1! — Main prognostic C state update
# ---------------------------------------------------------------------------

"""
    c_state_update1!(cs_veg, cf_veg, cf_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        patch_column, ivt, woody,
        cascade_donor_pool, cascade_receiver_pool,
        harvdate, col_is_fates,
        nlevdecomp, ndecomp_cascade_transitions,
        i_litr_min, i_litr_max, i_cwd,
        npcropmin, nrepr, repr_grain_min, repr_grain_max,
        repr_structure_min, repr_structure_max,
        use_matrixcn, use_soil_matrixcn,
        carbon_resp_opt, dribble_crophrv_xsmrpool_2atm, dt)

On the radiation time step, update all prognostic carbon state
variables (except for gap-phase mortality and fire fluxes).

Ported from `CStateUpdate1` in `CNCStateUpdate1Mod.F90`.

Note: The FATES litter flux callback (`clm_fates%UpdateCLitterfluxes`) is
skipped in this port — FATES columns are handled by zeroing
`decomp_cpools_sourcesink_col` via `col_is_fates`. The matrix-CN code paths
(`use_matrixcn`, `use_soil_matrixcn`) for soil/veg matrix solutions are
ported as-is, but the `matrix_update_phc` helper is replaced by a no-op
since it depends on the full matrix infrastructure.
"""
function c_state_update1!(cs_veg::CNVegCarbonStateData,
                           cf_veg::CNVegCarbonFluxData,
                           cf_soil::SoilBiogeochemCarbonFluxData;
                           mask_soilc::BitVector,
                           mask_soilp::BitVector,
                           bounds_col::UnitRange{Int},
                           bounds_patch::UnitRange{Int},
                           patch_column::Vector{Int},
                           ivt::Vector{Int},
                           woody::Vector{Float64},
                           cascade_donor_pool::Vector{Int},
                           cascade_receiver_pool::Vector{Int},
                           harvdate::Vector{Int},
                           col_is_fates::Vector{Bool},
                           nlevdecomp::Int,
                           ndecomp_cascade_transitions::Int,
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
                           carbon_resp_opt::Int = 0,
                           dribble_crophrv_xsmrpool_2atm::Bool = false,
                           dt::Float64)

    # Decay constant for 0.5-year product pool (1/s)
    kprod05 = 1.44e-7

    # --- Column loop: soil decomposition input fluxes ---
    for c in bounds_col
        mask_soilc[c] || continue

        if col_is_fates[c]
            # FATES columns: litter fluxes handled externally (skip)
            # In full CLM this calls clm_fates%UpdateCLitterfluxes
        else
            for j in 1:nlevdecomp
                if !use_soil_matrixcn
                    # phenology and dynamic land cover fluxes
                    for i in i_litr_min:i_litr_max
                        cf_soil.decomp_cpools_sourcesink_col[c, j, i] =
                            cf_veg.phenology_c_to_litr_c_col[c, j, i] * dt
                    end
                    # CWD: zeroed here (terms moved to CStateUpdateDynPatch)
                    cf_soil.decomp_cpools_sourcesink_col[c, j, i_cwd] = 0.0
                end
                # Note: matrix_Cinput path omitted (requires full matrix infrastructure)
            end
        end

        # Decomposition cascade HR and transfer fluxes
        for j in 1:nlevdecomp
            if !use_soil_matrixcn
                for k in 1:ndecomp_cascade_transitions
                    cf_soil.decomp_cpools_sourcesink_col[c, j, cascade_donor_pool[k]] =
                        cf_soil.decomp_cpools_sourcesink_col[c, j, cascade_donor_pool[k]] -
                        (cf_soil.decomp_cascade_hr_vr_col[c, j, k] +
                         cf_soil.decomp_cascade_ctransfer_vr_col[c, j, k]) * dt
                    if cascade_receiver_pool[k] != 0  # skip terminal transitions
                        cf_soil.decomp_cpools_sourcesink_col[c, j, cascade_receiver_pool[k]] =
                            cf_soil.decomp_cpools_sourcesink_col[c, j, cascade_receiver_pool[k]] +
                            cf_soil.decomp_cascade_ctransfer_vr_col[c, j, k] * dt
                    end
                end
            end
        end
    end

    # --- Patch loop: vegetation C state updates ---
    for p in bounds_patch
        mask_soilp[p] || continue
        c = patch_column[p]

        # === Phenology: transfer growth fluxes ===
        if !use_matrixcn
            cs_veg.leafc_patch[p]       = cs_veg.leafc_patch[p]       + cf_veg.leafc_xfer_to_leafc_patch[p] * dt
            cs_veg.leafc_xfer_patch[p]  = cs_veg.leafc_xfer_patch[p]  - cf_veg.leafc_xfer_to_leafc_patch[p] * dt
            cs_veg.frootc_patch[p]      = cs_veg.frootc_patch[p]      + cf_veg.frootc_xfer_to_frootc_patch[p] * dt
            cs_veg.frootc_xfer_patch[p] = cs_veg.frootc_xfer_patch[p] - cf_veg.frootc_xfer_to_frootc_patch[p] * dt

            if woody[ivt[p]] == 1.0
                cs_veg.livestemc_patch[p]       = cs_veg.livestemc_patch[p]       + cf_veg.livestemc_xfer_to_livestemc_patch[p] * dt
                cs_veg.livestemc_xfer_patch[p]  = cs_veg.livestemc_xfer_patch[p]  - cf_veg.livestemc_xfer_to_livestemc_patch[p] * dt
                cs_veg.deadstemc_patch[p]       = cs_veg.deadstemc_patch[p]       + cf_veg.deadstemc_xfer_to_deadstemc_patch[p] * dt
                cs_veg.deadstemc_xfer_patch[p]  = cs_veg.deadstemc_xfer_patch[p]  - cf_veg.deadstemc_xfer_to_deadstemc_patch[p] * dt
                cs_veg.livecrootc_patch[p]      = cs_veg.livecrootc_patch[p]      + cf_veg.livecrootc_xfer_to_livecrootc_patch[p] * dt
                cs_veg.livecrootc_xfer_patch[p] = cs_veg.livecrootc_xfer_patch[p] - cf_veg.livecrootc_xfer_to_livecrootc_patch[p] * dt
                cs_veg.deadcrootc_patch[p]      = cs_veg.deadcrootc_patch[p]      + cf_veg.deadcrootc_xfer_to_deadcrootc_patch[p] * dt
                cs_veg.deadcrootc_xfer_patch[p] = cs_veg.deadcrootc_xfer_patch[p] - cf_veg.deadcrootc_xfer_to_deadcrootc_patch[p] * dt
            end

            if ivt[p] >= npcropmin  # skip 2 generic crops
                # lines here for consistency; the transfer terms are zero
                cs_veg.livestemc_patch[p]      = cs_veg.livestemc_patch[p]      + cf_veg.livestemc_xfer_to_livestemc_patch[p] * dt
                cs_veg.livestemc_xfer_patch[p] = cs_veg.livestemc_xfer_patch[p] - cf_veg.livestemc_xfer_to_livestemc_patch[p] * dt
                for k in 1:nrepr
                    cs_veg.reproductivec_patch[p, k]      = cs_veg.reproductivec_patch[p, k] +
                        cf_veg.reproductivec_xfer_to_reproductivec_patch[p, k] * dt
                    cs_veg.reproductivec_xfer_patch[p, k] = cs_veg.reproductivec_xfer_patch[p, k] -
                        cf_veg.reproductivec_xfer_to_reproductivec_patch[p, k] * dt
                end
            end

            # phenology: litterfall fluxes
            cs_veg.leafc_patch[p]  = cs_veg.leafc_patch[p]  - cf_veg.leafc_to_litter_patch[p] * dt
            cs_veg.frootc_patch[p] = cs_veg.frootc_patch[p] - cf_veg.frootc_to_litter_patch[p] * dt

            # livewood turnover fluxes
            if woody[ivt[p]] == 1.0
                cs_veg.livestemc_patch[p]  = cs_veg.livestemc_patch[p]  - cf_veg.livestemc_to_deadstemc_patch[p] * dt
                cs_veg.deadstemc_patch[p]  = cs_veg.deadstemc_patch[p]  + cf_veg.livestemc_to_deadstemc_patch[p] * dt
                cs_veg.livecrootc_patch[p] = cs_veg.livecrootc_patch[p] - cf_veg.livecrootc_to_deadcrootc_patch[p] * dt
                cs_veg.deadcrootc_patch[p] = cs_veg.deadcrootc_patch[p] + cf_veg.livecrootc_to_deadcrootc_patch[p] * dt
            end

            if ivt[p] >= npcropmin  # skip 2 generic crops
                cs_veg.livestemc_patch[p] = cs_veg.livestemc_patch[p] - cf_veg.livestemc_to_litter_patch[p] * dt
                cs_veg.livestemc_patch[p] = cs_veg.livestemc_patch[p] -
                    (cf_veg.livestemc_to_biofuelc_patch[p] + cf_veg.livestemc_to_removedresiduec_patch[p]) * dt
                cs_veg.leafc_patch[p]     = cs_veg.leafc_patch[p] -
                    (cf_veg.leafc_to_biofuelc_patch[p] + cf_veg.leafc_to_removedresiduec_patch[p]) * dt
                cs_veg.cropseedc_deficit_patch[p] = cs_veg.cropseedc_deficit_patch[p] -
                    cf_veg.crop_seedc_to_leaf_patch[p] * dt
                for k in repr_grain_min:repr_grain_max
                    cs_veg.reproductivec_patch[p, k] = cs_veg.reproductivec_patch[p, k] -
                        (cf_veg.repr_grainc_to_food_patch[p, k] + cf_veg.repr_grainc_to_seed_patch[p, k]) * dt
                    cs_veg.cropseedc_deficit_patch[p] = cs_veg.cropseedc_deficit_patch[p] +
                        cf_veg.repr_grainc_to_seed_patch[p, k] * dt
                end
                for k in repr_structure_min:repr_structure_max
                    cs_veg.reproductivec_patch[p, k] = cs_veg.reproductivec_patch[p, k] -
                        (cf_veg.repr_structurec_to_cropprod_patch[p, k] + cf_veg.repr_structurec_to_litter_patch[p, k]) * dt
                end
            end
        else
            # Matrix version: only cropseedc_deficit updates
            if ivt[p] >= npcropmin
                cs_veg.cropseedc_deficit_patch[p] = cs_veg.cropseedc_deficit_patch[p] -
                    cf_veg.crop_seedc_to_leaf_patch[p] * dt
                for k in repr_grain_min:repr_grain_max
                    cs_veg.cropseedc_deficit_patch[p] = cs_veg.cropseedc_deficit_patch[p] +
                        cf_veg.repr_grainc_to_seed_patch[p, k] * dt
                end
            end
        end  # !use_matrixcn

        # === Maintenance respiration fluxes from cpool ===
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_xsmrpool_patch[p] * dt
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.leaf_curmr_patch[p] * dt
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.froot_curmr_patch[p] * dt

        if woody[ivt[p]] == 1.0
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.livestem_curmr_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.livecroot_curmr_patch[p] * dt
        end
        if ivt[p] >= npcropmin  # skip 2 generic crops
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.livestem_curmr_patch[p] * dt
            for k in 1:nrepr
                cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.reproductive_curmr_patch[p, k] * dt
            end
        end

        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_resp_patch[p] * dt

        # RF: carbon spent on uptake respiration
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.soilc_change_patch[p] * dt

        # === Maintenance respiration fluxes from xsmrpool ===
        cs_veg.xsmrpool_patch[p] = cs_veg.xsmrpool_patch[p] + cf_veg.cpool_to_xsmrpool_patch[p] * dt
        cs_veg.xsmrpool_patch[p] = cs_veg.xsmrpool_patch[p] - cf_veg.leaf_xsmr_patch[p] * dt
        cs_veg.xsmrpool_patch[p] = cs_veg.xsmrpool_patch[p] - cf_veg.froot_xsmr_patch[p] * dt
        if woody[ivt[p]] == 1.0
            cs_veg.xsmrpool_patch[p] = cs_veg.xsmrpool_patch[p] - cf_veg.livestem_xsmr_patch[p] * dt
            cs_veg.xsmrpool_patch[p] = cs_veg.xsmrpool_patch[p] - cf_veg.livecroot_xsmr_patch[p] * dt
        end

        # === Allocation fluxes ===
        if carbon_resp_opt == 1
            cf_veg.cpool_to_leafc_patch[p]          = cf_veg.cpool_to_leafc_patch[p]          - cf_veg.cpool_to_leafc_resp_patch[p]
            cf_veg.cpool_to_leafc_storage_patch[p]  = cf_veg.cpool_to_leafc_storage_patch[p]  - cf_veg.cpool_to_leafc_storage_resp_patch[p]
            cf_veg.cpool_to_frootc_patch[p]         = cf_veg.cpool_to_frootc_patch[p]         - cf_veg.cpool_to_frootc_resp_patch[p]
            cf_veg.cpool_to_frootc_storage_patch[p] = cf_veg.cpool_to_frootc_storage_patch[p] - cf_veg.cpool_to_frootc_storage_resp_patch[p]
        end

        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_leafc_patch[p] * dt
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_leafc_storage_patch[p] * dt
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_frootc_patch[p] * dt
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_frootc_storage_patch[p] * dt

        if !use_matrixcn
            cs_veg.leafc_patch[p]          = cs_veg.leafc_patch[p]          + cf_veg.cpool_to_leafc_patch[p] * dt
            cs_veg.leafc_storage_patch[p]  = cs_veg.leafc_storage_patch[p]  + cf_veg.cpool_to_leafc_storage_patch[p] * dt
            cs_veg.frootc_patch[p]         = cs_veg.frootc_patch[p]         + cf_veg.cpool_to_frootc_patch[p] * dt
            cs_veg.frootc_storage_patch[p] = cs_veg.frootc_storage_patch[p] + cf_veg.cpool_to_frootc_storage_patch[p] * dt
        end

        if woody[ivt[p]] == 1.0
            if carbon_resp_opt == 1
                cf_veg.cpool_to_livecrootc_patch[p]         = cf_veg.cpool_to_livecrootc_patch[p]         - cf_veg.cpool_to_livecrootc_resp_patch[p]
                cf_veg.cpool_to_livecrootc_storage_patch[p] = cf_veg.cpool_to_livecrootc_storage_patch[p] - cf_veg.cpool_to_livecrootc_storage_resp_patch[p]
                cf_veg.cpool_to_livestemc_patch[p]          = cf_veg.cpool_to_livestemc_patch[p]          - cf_veg.cpool_to_livestemc_resp_patch[p]
                cf_veg.cpool_to_livestemc_storage_patch[p]  = cf_veg.cpool_to_livestemc_storage_patch[p]  - cf_veg.cpool_to_livestemc_storage_resp_patch[p]
            end

            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_livestemc_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_livestemc_storage_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_deadstemc_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_deadstemc_storage_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_livecrootc_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_livecrootc_storage_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_deadcrootc_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_deadcrootc_storage_patch[p] * dt

            if !use_matrixcn
                cs_veg.livestemc_patch[p]          = cs_veg.livestemc_patch[p]          + cf_veg.cpool_to_livestemc_patch[p] * dt
                cs_veg.livestemc_storage_patch[p]  = cs_veg.livestemc_storage_patch[p]  + cf_veg.cpool_to_livestemc_storage_patch[p] * dt
                cs_veg.deadstemc_patch[p]          = cs_veg.deadstemc_patch[p]          + cf_veg.cpool_to_deadstemc_patch[p] * dt
                cs_veg.deadstemc_storage_patch[p]  = cs_veg.deadstemc_storage_patch[p]  + cf_veg.cpool_to_deadstemc_storage_patch[p] * dt
                cs_veg.livecrootc_patch[p]         = cs_veg.livecrootc_patch[p]         + cf_veg.cpool_to_livecrootc_patch[p] * dt
                cs_veg.livecrootc_storage_patch[p] = cs_veg.livecrootc_storage_patch[p] + cf_veg.cpool_to_livecrootc_storage_patch[p] * dt
                cs_veg.deadcrootc_patch[p]         = cs_veg.deadcrootc_patch[p]         + cf_veg.cpool_to_deadcrootc_patch[p] * dt
                cs_veg.deadcrootc_storage_patch[p] = cs_veg.deadcrootc_storage_patch[p] + cf_veg.cpool_to_deadcrootc_storage_patch[p] * dt
            end
        end

        if ivt[p] >= npcropmin  # skip 2 generic crops
            if carbon_resp_opt == 1
                cf_veg.cpool_to_livestemc_patch[p]         = cf_veg.cpool_to_livestemc_patch[p]         - cf_veg.cpool_to_livestemc_resp_patch[p]
                cf_veg.cpool_to_livestemc_storage_patch[p] = cf_veg.cpool_to_livestemc_storage_patch[p] - cf_veg.cpool_to_livestemc_storage_resp_patch[p]
            end

            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_livestemc_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_livestemc_storage_patch[p] * dt
            for k in 1:nrepr
                cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_reproductivec_patch[p, k] * dt
                cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_to_reproductivec_storage_patch[p, k] * dt
            end

            if !use_matrixcn
                cs_veg.livestemc_patch[p]         = cs_veg.livestemc_patch[p]         + cf_veg.cpool_to_livestemc_patch[p] * dt
                cs_veg.livestemc_storage_patch[p] = cs_veg.livestemc_storage_patch[p] + cf_veg.cpool_to_livestemc_storage_patch[p] * dt
                for k in 1:nrepr
                    cs_veg.reproductivec_patch[p, k]         = cs_veg.reproductivec_patch[p, k] +
                        cf_veg.cpool_to_reproductivec_patch[p, k] * dt
                    cs_veg.reproductivec_storage_patch[p, k] = cs_veg.reproductivec_storage_patch[p, k] +
                        cf_veg.cpool_to_reproductivec_storage_patch[p, k] * dt
                end
            end
        end

        # === Growth respiration fluxes for current growth ===
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_leaf_gr_patch[p] * dt
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_froot_gr_patch[p] * dt

        if woody[ivt[p]] == 1.0
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_livestem_gr_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_deadstem_gr_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_livecroot_gr_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_deadcroot_gr_patch[p] * dt
        end
        if ivt[p] >= npcropmin  # skip 2 generic crops
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_livestem_gr_patch[p] * dt
            for k in 1:nrepr
                cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_reproductive_gr_patch[p, k] * dt
            end
        end

        # === Growth respiration for transfer growth ===
        cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] - cf_veg.transfer_leaf_gr_patch[p] * dt
        cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] - cf_veg.transfer_froot_gr_patch[p] * dt
        if woody[ivt[p]] == 1.0
            cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] - cf_veg.transfer_livestem_gr_patch[p] * dt
            cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] - cf_veg.transfer_deadstem_gr_patch[p] * dt
            cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] - cf_veg.transfer_livecroot_gr_patch[p] * dt
            cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] - cf_veg.transfer_deadcroot_gr_patch[p] * dt
        end
        if ivt[p] >= npcropmin  # skip 2 generic crops
            cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] - cf_veg.transfer_livestem_gr_patch[p] * dt
            for k in 1:nrepr
                cs_veg.gresp_xfer_patch[p] = cs_veg.gresp_xfer_patch[p] - cf_veg.transfer_reproductive_gr_patch[p, k] * dt
            end
        end

        # === Growth respiration at time of storage ===
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_leaf_storage_gr_patch[p] * dt
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_froot_storage_gr_patch[p] * dt

        if woody[ivt[p]] == 1.0
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_livestem_storage_gr_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_deadstem_storage_gr_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_livecroot_storage_gr_patch[p] * dt
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_deadcroot_storage_gr_patch[p] * dt
        end
        if ivt[p] >= npcropmin  # skip 2 generic crops
            cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_livestem_storage_gr_patch[p] * dt
            for k in 1:nrepr
                cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] - cf_veg.cpool_reproductive_storage_gr_patch[p, k] * dt
            end
        end

        # === Growth respiration stored for release during transfer growth ===
        cs_veg.cpool_patch[p]         = cs_veg.cpool_patch[p]         - cf_veg.cpool_to_gresp_storage_patch[p] * dt
        cs_veg.gresp_storage_patch[p] = cs_veg.gresp_storage_patch[p] + cf_veg.cpool_to_gresp_storage_patch[p] * dt

        # === Move storage pools into transfer pools ===
        if !use_matrixcn
            cs_veg.leafc_storage_patch[p]  = cs_veg.leafc_storage_patch[p]  - cf_veg.leafc_storage_to_xfer_patch[p] * dt
            cs_veg.leafc_xfer_patch[p]     = cs_veg.leafc_xfer_patch[p]     + cf_veg.leafc_storage_to_xfer_patch[p] * dt
            cs_veg.frootc_storage_patch[p] = cs_veg.frootc_storage_patch[p] - cf_veg.frootc_storage_to_xfer_patch[p] * dt
            cs_veg.frootc_xfer_patch[p]    = cs_veg.frootc_xfer_patch[p]    + cf_veg.frootc_storage_to_xfer_patch[p] * dt
        end

        if woody[ivt[p]] == 1.0
            cs_veg.gresp_storage_patch[p] = cs_veg.gresp_storage_patch[p] - cf_veg.gresp_storage_to_xfer_patch[p] * dt
            cs_veg.gresp_xfer_patch[p]    = cs_veg.gresp_xfer_patch[p]    + cf_veg.gresp_storage_to_xfer_patch[p] * dt
            if !use_matrixcn
                cs_veg.livestemc_storage_patch[p]  = cs_veg.livestemc_storage_patch[p]  - cf_veg.livestemc_storage_to_xfer_patch[p] * dt
                cs_veg.livestemc_xfer_patch[p]     = cs_veg.livestemc_xfer_patch[p]     + cf_veg.livestemc_storage_to_xfer_patch[p] * dt
                cs_veg.deadstemc_storage_patch[p]  = cs_veg.deadstemc_storage_patch[p]  - cf_veg.deadstemc_storage_to_xfer_patch[p] * dt
                cs_veg.deadstemc_xfer_patch[p]     = cs_veg.deadstemc_xfer_patch[p]     + cf_veg.deadstemc_storage_to_xfer_patch[p] * dt
                cs_veg.livecrootc_storage_patch[p] = cs_veg.livecrootc_storage_patch[p] - cf_veg.livecrootc_storage_to_xfer_patch[p] * dt
                cs_veg.livecrootc_xfer_patch[p]    = cs_veg.livecrootc_xfer_patch[p]    + cf_veg.livecrootc_storage_to_xfer_patch[p] * dt
                cs_veg.deadcrootc_storage_patch[p] = cs_veg.deadcrootc_storage_patch[p] - cf_veg.deadcrootc_storage_to_xfer_patch[p] * dt
                cs_veg.deadcrootc_xfer_patch[p]    = cs_veg.deadcrootc_xfer_patch[p]    + cf_veg.deadcrootc_storage_to_xfer_patch[p] * dt
            end
        end

        if ivt[p] >= npcropmin  # skip 2 generic crops
            if !use_matrixcn
                # lines here for consistency; the transfer terms are zero
                cs_veg.livestemc_storage_patch[p] = cs_veg.livestemc_storage_patch[p] - cf_veg.livestemc_storage_to_xfer_patch[p] * dt
                cs_veg.livestemc_xfer_patch[p]    = cs_veg.livestemc_xfer_patch[p]    + cf_veg.livestemc_storage_to_xfer_patch[p] * dt
                for k in 1:nrepr
                    cs_veg.reproductivec_storage_patch[p, k] = cs_veg.reproductivec_storage_patch[p, k] -
                        cf_veg.reproductivec_storage_to_xfer_patch[p, k] * dt
                    cs_veg.reproductivec_xfer_patch[p, k]    = cs_veg.reproductivec_xfer_patch[p, k] +
                        cf_veg.reproductivec_storage_to_xfer_patch[p, k] * dt
                end
            end
        end

        # === Crop-specific xsmrpool and harvest logic ===
        if ivt[p] >= npcropmin  # skip 2 generic crops
            cs_veg.xsmrpool_patch[p] = cs_veg.xsmrpool_patch[p] - cf_veg.livestem_xsmr_patch[p] * dt
            for k in 1:nrepr
                cs_veg.xsmrpool_patch[p] = cs_veg.xsmrpool_patch[p] - cf_veg.reproductive_xsmr_patch[p, k] * dt
            end

            if harvdate[p] < 999  # beginning at harvest, send to atm
                if !dribble_crophrv_xsmrpool_2atm
                    cf_veg.xsmrpool_to_atm_patch[p] = cf_veg.xsmrpool_to_atm_patch[p] + cs_veg.xsmrpool_patch[p] / dt
                    cf_veg.xsmrpool_to_atm_patch[p] = cf_veg.xsmrpool_to_atm_patch[p] + cs_veg.cpool_patch[p] / dt
                    if !use_matrixcn
                        cf_veg.xsmrpool_to_atm_patch[p] = cf_veg.xsmrpool_to_atm_patch[p] + cs_veg.frootc_patch[p] / dt
                    end
                    # Note: matrix_update_phc path omitted (requires full matrix infrastructure)
                else
                    cs_veg.xsmrpool_loss_patch[p] = cs_veg.xsmrpool_loss_patch[p] +
                        cs_veg.xsmrpool_patch[p] +
                        cs_veg.cpool_patch[p]
                    if !use_matrixcn
                        cs_veg.xsmrpool_loss_patch[p] = cs_veg.xsmrpool_loss_patch[p] + cs_veg.frootc_patch[p]
                    end
                    # Note: matrix_update_phc path omitted (requires full matrix infrastructure)
                end
                if !use_matrixcn
                    cs_veg.frootc_patch[p] = 0.0
                end
                cs_veg.xsmrpool_patch[p] = 0.0
                cs_veg.cpool_patch[p]    = 0.0
            end

            # Slowly release xsmrpool to atmosphere
            if dribble_crophrv_xsmrpool_2atm
                cf_veg.xsmrpool_to_atm_patch[p] = cs_veg.xsmrpool_loss_patch[p] * kprod05
                cs_veg.xsmrpool_loss_patch[p]   = cs_veg.xsmrpool_loss_patch[p] - cf_veg.xsmrpool_to_atm_patch[p] * dt
            end
        end
    end  # end of patch loop

    return nothing
end
