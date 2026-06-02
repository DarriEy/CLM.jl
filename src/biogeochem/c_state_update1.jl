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
                                    dt::Real)

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
                           dt::Real)

    for p in bounds_patch
        mask_soilp[p] || continue
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] + cf_veg.psnsun_to_cpool_patch[p] * dt
        cs_veg.cpool_patch[p] = cs_veg.cpool_patch[p] + cf_veg.psnshade_to_cpool_patch[p] * dt
    end

    return nothing
end

# ---------------------------------------------------------------------------
# c_state_update1! — Main prognostic C state update
#
# GPU kernelization (Phase B pattern-setter for the C/N state-update cascade):
#   * Column loop  -> `_csu1_col_kernel!` (one thread per soil column; the
#     internal j/k decomposition-cascade RMW into decomp_cpools_sourcesink_col
#     runs sequentially in-thread, so it is byte-identical to the host loop and
#     race-free across threads — each thread owns its column).
#   * Patch loop   -> `_csu1_patch_kernel!` (one thread per patch; every write is
#     to the patch's own index, no patch->column scatter). The ~24+3 carbon-state
#     and ~73+13 carbon-flux arrays the body touches are grouped into two
#     immutable `@adapt_structure` device-view bundles (`_CSU1CS{V,M}` /
#     `_CSU1CF{V,M}`) so the launch stays well under Metal's ~31-arg limit; the
#     bundle field names mirror the state structs so the body reads verbatim.
# Float64 literals are eltype-converted (`zero(T)`, `T(...)`) so the kernels carry
# no Float64 on a Float32-only backend; on Float64 this is byte-identical.
# ---------------------------------------------------------------------------

# --- Device-view bundle: CNVeg carbon STATE arrays the patch loop reads/writes ---
# `V` = the 1D patch vectors, `M` = the 2D (patch x nrepr) reproductive matrices.
Base.@kwdef struct _CSU1CS{V,M}
    leafc_patch::V; leafc_xfer_patch::V; leafc_storage_patch::V
    frootc_patch::V; frootc_xfer_patch::V; frootc_storage_patch::V
    livestemc_patch::V; livestemc_xfer_patch::V; livestemc_storage_patch::V
    deadstemc_patch::V; deadstemc_xfer_patch::V; deadstemc_storage_patch::V
    livecrootc_patch::V; livecrootc_xfer_patch::V; livecrootc_storage_patch::V
    deadcrootc_patch::V; deadcrootc_xfer_patch::V; deadcrootc_storage_patch::V
    cpool_patch::V; xsmrpool_patch::V; xsmrpool_loss_patch::V
    gresp_storage_patch::V; gresp_xfer_patch::V; cropseedc_deficit_patch::V
    reproductivec_patch::M; reproductivec_storage_patch::M; reproductivec_xfer_patch::M
end
Adapt.@adapt_structure _CSU1CS

_csu1_cs(cs) = _CSU1CS(;
    leafc_patch = cs.leafc_patch, leafc_xfer_patch = cs.leafc_xfer_patch,
    leafc_storage_patch = cs.leafc_storage_patch,
    frootc_patch = cs.frootc_patch, frootc_xfer_patch = cs.frootc_xfer_patch,
    frootc_storage_patch = cs.frootc_storage_patch,
    livestemc_patch = cs.livestemc_patch, livestemc_xfer_patch = cs.livestemc_xfer_patch,
    livestemc_storage_patch = cs.livestemc_storage_patch,
    deadstemc_patch = cs.deadstemc_patch, deadstemc_xfer_patch = cs.deadstemc_xfer_patch,
    deadstemc_storage_patch = cs.deadstemc_storage_patch,
    livecrootc_patch = cs.livecrootc_patch, livecrootc_xfer_patch = cs.livecrootc_xfer_patch,
    livecrootc_storage_patch = cs.livecrootc_storage_patch,
    deadcrootc_patch = cs.deadcrootc_patch, deadcrootc_xfer_patch = cs.deadcrootc_xfer_patch,
    deadcrootc_storage_patch = cs.deadcrootc_storage_patch,
    cpool_patch = cs.cpool_patch, xsmrpool_patch = cs.xsmrpool_patch,
    xsmrpool_loss_patch = cs.xsmrpool_loss_patch,
    gresp_storage_patch = cs.gresp_storage_patch, gresp_xfer_patch = cs.gresp_xfer_patch,
    cropseedc_deficit_patch = cs.cropseedc_deficit_patch,
    reproductivec_patch = cs.reproductivec_patch,
    reproductivec_storage_patch = cs.reproductivec_storage_patch,
    reproductivec_xfer_patch = cs.reproductivec_xfer_patch)

# --- Device-view bundle: CNVeg carbon FLUX arrays the patch loop reads/writes ---
Base.@kwdef struct _CSU1CF{V,M}
    # 1D patch vectors
    leafc_xfer_to_leafc_patch::V; frootc_xfer_to_frootc_patch::V
    livestemc_xfer_to_livestemc_patch::V; deadstemc_xfer_to_deadstemc_patch::V
    livecrootc_xfer_to_livecrootc_patch::V; deadcrootc_xfer_to_deadcrootc_patch::V
    leafc_to_litter_patch::V; frootc_to_litter_patch::V
    livestemc_to_deadstemc_patch::V; livecrootc_to_deadcrootc_patch::V
    livestemc_to_litter_patch::V; livestemc_to_biofuelc_patch::V
    livestemc_to_removedresiduec_patch::V
    leafc_to_biofuelc_patch::V; leafc_to_removedresiduec_patch::V
    crop_seedc_to_leaf_patch::V
    cpool_to_xsmrpool_patch::V; leaf_curmr_patch::V; froot_curmr_patch::V
    livestem_curmr_patch::V; livecroot_curmr_patch::V
    cpool_to_resp_patch::V; soilc_change_patch::V
    leaf_xsmr_patch::V; froot_xsmr_patch::V; livestem_xsmr_patch::V; livecroot_xsmr_patch::V
    cpool_to_leafc_patch::V; cpool_to_leafc_resp_patch::V
    cpool_to_leafc_storage_patch::V; cpool_to_leafc_storage_resp_patch::V
    cpool_to_frootc_patch::V; cpool_to_frootc_resp_patch::V
    cpool_to_frootc_storage_patch::V; cpool_to_frootc_storage_resp_patch::V
    cpool_to_livecrootc_patch::V; cpool_to_livecrootc_resp_patch::V
    cpool_to_livecrootc_storage_patch::V; cpool_to_livecrootc_storage_resp_patch::V
    cpool_to_livestemc_patch::V; cpool_to_livestemc_resp_patch::V
    cpool_to_livestemc_storage_patch::V; cpool_to_livestemc_storage_resp_patch::V
    cpool_to_deadstemc_patch::V; cpool_to_deadstemc_storage_patch::V
    cpool_to_deadcrootc_patch::V; cpool_to_deadcrootc_storage_patch::V
    cpool_leaf_gr_patch::V; cpool_froot_gr_patch::V
    cpool_livestem_gr_patch::V; cpool_deadstem_gr_patch::V
    cpool_livecroot_gr_patch::V; cpool_deadcroot_gr_patch::V
    transfer_leaf_gr_patch::V; transfer_froot_gr_patch::V
    transfer_livestem_gr_patch::V; transfer_deadstem_gr_patch::V
    transfer_livecroot_gr_patch::V; transfer_deadcroot_gr_patch::V
    cpool_leaf_storage_gr_patch::V; cpool_froot_storage_gr_patch::V
    cpool_livestem_storage_gr_patch::V; cpool_deadstem_storage_gr_patch::V
    cpool_livecroot_storage_gr_patch::V; cpool_deadcroot_storage_gr_patch::V
    cpool_to_gresp_storage_patch::V
    leafc_storage_to_xfer_patch::V; frootc_storage_to_xfer_patch::V
    livestemc_storage_to_xfer_patch::V; deadstemc_storage_to_xfer_patch::V
    livecrootc_storage_to_xfer_patch::V; deadcrootc_storage_to_xfer_patch::V
    gresp_storage_to_xfer_patch::V
    xsmrpool_to_atm_patch::V
    # 2D (patch x nrepr) reproductive matrices
    reproductivec_xfer_to_reproductivec_patch::M
    repr_grainc_to_food_patch::M; repr_grainc_to_seed_patch::M
    repr_structurec_to_cropprod_patch::M; repr_structurec_to_litter_patch::M
    reproductive_curmr_patch::M; reproductive_xsmr_patch::M
    cpool_to_reproductivec_patch::M; cpool_to_reproductivec_storage_patch::M
    cpool_reproductive_gr_patch::M; cpool_reproductive_storage_gr_patch::M
    transfer_reproductive_gr_patch::M
    reproductivec_storage_to_xfer_patch::M
end
Adapt.@adapt_structure _CSU1CF

_csu1_cf(cf) = _CSU1CF(;
    leafc_xfer_to_leafc_patch = cf.leafc_xfer_to_leafc_patch,
    frootc_xfer_to_frootc_patch = cf.frootc_xfer_to_frootc_patch,
    livestemc_xfer_to_livestemc_patch = cf.livestemc_xfer_to_livestemc_patch,
    deadstemc_xfer_to_deadstemc_patch = cf.deadstemc_xfer_to_deadstemc_patch,
    livecrootc_xfer_to_livecrootc_patch = cf.livecrootc_xfer_to_livecrootc_patch,
    deadcrootc_xfer_to_deadcrootc_patch = cf.deadcrootc_xfer_to_deadcrootc_patch,
    leafc_to_litter_patch = cf.leafc_to_litter_patch,
    frootc_to_litter_patch = cf.frootc_to_litter_patch,
    livestemc_to_deadstemc_patch = cf.livestemc_to_deadstemc_patch,
    livecrootc_to_deadcrootc_patch = cf.livecrootc_to_deadcrootc_patch,
    livestemc_to_litter_patch = cf.livestemc_to_litter_patch,
    livestemc_to_biofuelc_patch = cf.livestemc_to_biofuelc_patch,
    livestemc_to_removedresiduec_patch = cf.livestemc_to_removedresiduec_patch,
    leafc_to_biofuelc_patch = cf.leafc_to_biofuelc_patch,
    leafc_to_removedresiduec_patch = cf.leafc_to_removedresiduec_patch,
    crop_seedc_to_leaf_patch = cf.crop_seedc_to_leaf_patch,
    cpool_to_xsmrpool_patch = cf.cpool_to_xsmrpool_patch,
    leaf_curmr_patch = cf.leaf_curmr_patch, froot_curmr_patch = cf.froot_curmr_patch,
    livestem_curmr_patch = cf.livestem_curmr_patch,
    livecroot_curmr_patch = cf.livecroot_curmr_patch,
    cpool_to_resp_patch = cf.cpool_to_resp_patch, soilc_change_patch = cf.soilc_change_patch,
    leaf_xsmr_patch = cf.leaf_xsmr_patch, froot_xsmr_patch = cf.froot_xsmr_patch,
    livestem_xsmr_patch = cf.livestem_xsmr_patch, livecroot_xsmr_patch = cf.livecroot_xsmr_patch,
    cpool_to_leafc_patch = cf.cpool_to_leafc_patch,
    cpool_to_leafc_resp_patch = cf.cpool_to_leafc_resp_patch,
    cpool_to_leafc_storage_patch = cf.cpool_to_leafc_storage_patch,
    cpool_to_leafc_storage_resp_patch = cf.cpool_to_leafc_storage_resp_patch,
    cpool_to_frootc_patch = cf.cpool_to_frootc_patch,
    cpool_to_frootc_resp_patch = cf.cpool_to_frootc_resp_patch,
    cpool_to_frootc_storage_patch = cf.cpool_to_frootc_storage_patch,
    cpool_to_frootc_storage_resp_patch = cf.cpool_to_frootc_storage_resp_patch,
    cpool_to_livecrootc_patch = cf.cpool_to_livecrootc_patch,
    cpool_to_livecrootc_resp_patch = cf.cpool_to_livecrootc_resp_patch,
    cpool_to_livecrootc_storage_patch = cf.cpool_to_livecrootc_storage_patch,
    cpool_to_livecrootc_storage_resp_patch = cf.cpool_to_livecrootc_storage_resp_patch,
    cpool_to_livestemc_patch = cf.cpool_to_livestemc_patch,
    cpool_to_livestemc_resp_patch = cf.cpool_to_livestemc_resp_patch,
    cpool_to_livestemc_storage_patch = cf.cpool_to_livestemc_storage_patch,
    cpool_to_livestemc_storage_resp_patch = cf.cpool_to_livestemc_storage_resp_patch,
    cpool_to_deadstemc_patch = cf.cpool_to_deadstemc_patch,
    cpool_to_deadstemc_storage_patch = cf.cpool_to_deadstemc_storage_patch,
    cpool_to_deadcrootc_patch = cf.cpool_to_deadcrootc_patch,
    cpool_to_deadcrootc_storage_patch = cf.cpool_to_deadcrootc_storage_patch,
    cpool_leaf_gr_patch = cf.cpool_leaf_gr_patch, cpool_froot_gr_patch = cf.cpool_froot_gr_patch,
    cpool_livestem_gr_patch = cf.cpool_livestem_gr_patch,
    cpool_deadstem_gr_patch = cf.cpool_deadstem_gr_patch,
    cpool_livecroot_gr_patch = cf.cpool_livecroot_gr_patch,
    cpool_deadcroot_gr_patch = cf.cpool_deadcroot_gr_patch,
    transfer_leaf_gr_patch = cf.transfer_leaf_gr_patch,
    transfer_froot_gr_patch = cf.transfer_froot_gr_patch,
    transfer_livestem_gr_patch = cf.transfer_livestem_gr_patch,
    transfer_deadstem_gr_patch = cf.transfer_deadstem_gr_patch,
    transfer_livecroot_gr_patch = cf.transfer_livecroot_gr_patch,
    transfer_deadcroot_gr_patch = cf.transfer_deadcroot_gr_patch,
    cpool_leaf_storage_gr_patch = cf.cpool_leaf_storage_gr_patch,
    cpool_froot_storage_gr_patch = cf.cpool_froot_storage_gr_patch,
    cpool_livestem_storage_gr_patch = cf.cpool_livestem_storage_gr_patch,
    cpool_deadstem_storage_gr_patch = cf.cpool_deadstem_storage_gr_patch,
    cpool_livecroot_storage_gr_patch = cf.cpool_livecroot_storage_gr_patch,
    cpool_deadcroot_storage_gr_patch = cf.cpool_deadcroot_storage_gr_patch,
    cpool_to_gresp_storage_patch = cf.cpool_to_gresp_storage_patch,
    leafc_storage_to_xfer_patch = cf.leafc_storage_to_xfer_patch,
    frootc_storage_to_xfer_patch = cf.frootc_storage_to_xfer_patch,
    livestemc_storage_to_xfer_patch = cf.livestemc_storage_to_xfer_patch,
    deadstemc_storage_to_xfer_patch = cf.deadstemc_storage_to_xfer_patch,
    livecrootc_storage_to_xfer_patch = cf.livecrootc_storage_to_xfer_patch,
    deadcrootc_storage_to_xfer_patch = cf.deadcrootc_storage_to_xfer_patch,
    gresp_storage_to_xfer_patch = cf.gresp_storage_to_xfer_patch,
    xsmrpool_to_atm_patch = cf.xsmrpool_to_atm_patch,
    reproductivec_xfer_to_reproductivec_patch = cf.reproductivec_xfer_to_reproductivec_patch,
    repr_grainc_to_food_patch = cf.repr_grainc_to_food_patch,
    repr_grainc_to_seed_patch = cf.repr_grainc_to_seed_patch,
    repr_structurec_to_cropprod_patch = cf.repr_structurec_to_cropprod_patch,
    repr_structurec_to_litter_patch = cf.repr_structurec_to_litter_patch,
    reproductive_curmr_patch = cf.reproductive_curmr_patch,
    reproductive_xsmr_patch = cf.reproductive_xsmr_patch,
    cpool_to_reproductivec_patch = cf.cpool_to_reproductivec_patch,
    cpool_to_reproductivec_storage_patch = cf.cpool_to_reproductivec_storage_patch,
    cpool_reproductive_gr_patch = cf.cpool_reproductive_gr_patch,
    cpool_reproductive_storage_gr_patch = cf.cpool_reproductive_storage_gr_patch,
    transfer_reproductive_gr_patch = cf.transfer_reproductive_gr_patch,
    reproductivec_storage_to_xfer_patch = cf.reproductivec_storage_to_xfer_patch)

# --- Kernel: column-level soil decomposition input fluxes (one thread per column) ---
@kernel function _csu1_col_kernel!(@Const(mask_soilc), @Const(col_is_fates),
        @Const(phenology_c_to_litr_c_col), decomp_cpools_sourcesink_col,
        @Const(decomp_cascade_hr_vr_col), @Const(decomp_cascade_ctransfer_vr_col),
        @Const(cascade_donor_pool), @Const(cascade_receiver_pool),
        nlevdecomp::Int, ndecomp_cascade_transitions::Int,
        i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        use_soil_matrixcn::Bool, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        T = typeof(dt)
        if !col_is_fates[c]
            for j in 1:nlevdecomp
                if !use_soil_matrixcn
                    # phenology and dynamic land cover fluxes
                    for i in i_litr_min:i_litr_max
                        decomp_cpools_sourcesink_col[c, j, i] =
                            phenology_c_to_litr_c_col[c, j, i] * dt
                    end
                    # CWD: zeroed here (terms moved to CStateUpdateDynPatch)
                    decomp_cpools_sourcesink_col[c, j, i_cwd] = zero(T)
                end
            end
        end

        # Decomposition cascade HR and transfer fluxes
        for j in 1:nlevdecomp
            if !use_soil_matrixcn
                for k in 1:ndecomp_cascade_transitions
                    decomp_cpools_sourcesink_col[c, j, cascade_donor_pool[k]] =
                        decomp_cpools_sourcesink_col[c, j, cascade_donor_pool[k]] -
                        (decomp_cascade_hr_vr_col[c, j, k] +
                         decomp_cascade_ctransfer_vr_col[c, j, k]) * dt
                    if cascade_receiver_pool[k] != 0  # skip terminal transitions
                        decomp_cpools_sourcesink_col[c, j, cascade_receiver_pool[k]] =
                            decomp_cpools_sourcesink_col[c, j, cascade_receiver_pool[k]] +
                            decomp_cascade_ctransfer_vr_col[c, j, k] * dt
                    end
                end
            end
        end
    end
end

# --- Kernel: patch-level vegetation C state updates (one thread per patch) ---
# `cs` / `cf` are the _CSU1CS / _CSU1CF device-view bundles; field names mirror the
# state structs so the body below is the host loop verbatim (cs_veg.->cs., cf_veg.->cf.).
@kernel function _csu1_patch_kernel!(@Const(mask_soilp), @Const(ivt), @Const(woody),
        @Const(harvdate), cs, cf,
        npcropmin::Int, nrepr::Int, repr_grain_min::Int, repr_grain_max::Int,
        repr_structure_min::Int, repr_structure_max::Int,
        use_matrixcn::Bool, carbon_resp_opt::Int,
        dribble_crophrv_xsmrpool_2atm::Bool, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p] && ivt[p] >= 1  # skip masked + bare ground (PFT index 0)
        T = typeof(dt)
        kprod05 = T(1.44e-7)            # decay constant for 0.5-year product pool (1/s)
        is_woody = woody[ivt[p]] == one(eltype(woody))

        # === Phenology: transfer growth fluxes ===
        if !use_matrixcn
            cs.leafc_patch[p]       = cs.leafc_patch[p]       + cf.leafc_xfer_to_leafc_patch[p] * dt
            cs.leafc_xfer_patch[p]  = cs.leafc_xfer_patch[p]  - cf.leafc_xfer_to_leafc_patch[p] * dt
            cs.frootc_patch[p]      = cs.frootc_patch[p]      + cf.frootc_xfer_to_frootc_patch[p] * dt
            cs.frootc_xfer_patch[p] = cs.frootc_xfer_patch[p] - cf.frootc_xfer_to_frootc_patch[p] * dt

            if is_woody
                cs.livestemc_patch[p]       = cs.livestemc_patch[p]       + cf.livestemc_xfer_to_livestemc_patch[p] * dt
                cs.livestemc_xfer_patch[p]  = cs.livestemc_xfer_patch[p]  - cf.livestemc_xfer_to_livestemc_patch[p] * dt
                cs.deadstemc_patch[p]       = cs.deadstemc_patch[p]       + cf.deadstemc_xfer_to_deadstemc_patch[p] * dt
                cs.deadstemc_xfer_patch[p]  = cs.deadstemc_xfer_patch[p]  - cf.deadstemc_xfer_to_deadstemc_patch[p] * dt
                cs.livecrootc_patch[p]      = cs.livecrootc_patch[p]      + cf.livecrootc_xfer_to_livecrootc_patch[p] * dt
                cs.livecrootc_xfer_patch[p] = cs.livecrootc_xfer_patch[p] - cf.livecrootc_xfer_to_livecrootc_patch[p] * dt
                cs.deadcrootc_patch[p]      = cs.deadcrootc_patch[p]      + cf.deadcrootc_xfer_to_deadcrootc_patch[p] * dt
                cs.deadcrootc_xfer_patch[p] = cs.deadcrootc_xfer_patch[p] - cf.deadcrootc_xfer_to_deadcrootc_patch[p] * dt
            end

            if ivt[p] >= npcropmin  # skip 2 generic crops
                # lines here for consistency; the transfer terms are zero
                cs.livestemc_patch[p]      = cs.livestemc_patch[p]      + cf.livestemc_xfer_to_livestemc_patch[p] * dt
                cs.livestemc_xfer_patch[p] = cs.livestemc_xfer_patch[p] - cf.livestemc_xfer_to_livestemc_patch[p] * dt
                for k in 1:nrepr
                    cs.reproductivec_patch[p, k]      = cs.reproductivec_patch[p, k] +
                        cf.reproductivec_xfer_to_reproductivec_patch[p, k] * dt
                    cs.reproductivec_xfer_patch[p, k] = cs.reproductivec_xfer_patch[p, k] -
                        cf.reproductivec_xfer_to_reproductivec_patch[p, k] * dt
                end
            end

            # phenology: litterfall fluxes
            cs.leafc_patch[p]  = cs.leafc_patch[p]  - cf.leafc_to_litter_patch[p] * dt
            cs.frootc_patch[p] = cs.frootc_patch[p] - cf.frootc_to_litter_patch[p] * dt

            # livewood turnover fluxes
            if is_woody
                cs.livestemc_patch[p]  = cs.livestemc_patch[p]  - cf.livestemc_to_deadstemc_patch[p] * dt
                cs.deadstemc_patch[p]  = cs.deadstemc_patch[p]  + cf.livestemc_to_deadstemc_patch[p] * dt
                cs.livecrootc_patch[p] = cs.livecrootc_patch[p] - cf.livecrootc_to_deadcrootc_patch[p] * dt
                cs.deadcrootc_patch[p] = cs.deadcrootc_patch[p] + cf.livecrootc_to_deadcrootc_patch[p] * dt
            end

            if ivt[p] >= npcropmin  # skip 2 generic crops
                cs.livestemc_patch[p] = cs.livestemc_patch[p] - cf.livestemc_to_litter_patch[p] * dt
                cs.livestemc_patch[p] = cs.livestemc_patch[p] -
                    (cf.livestemc_to_biofuelc_patch[p] + cf.livestemc_to_removedresiduec_patch[p]) * dt
                cs.leafc_patch[p]     = cs.leafc_patch[p] -
                    (cf.leafc_to_biofuelc_patch[p] + cf.leafc_to_removedresiduec_patch[p]) * dt
                cs.cropseedc_deficit_patch[p] = cs.cropseedc_deficit_patch[p] -
                    cf.crop_seedc_to_leaf_patch[p] * dt
                for k in repr_grain_min:repr_grain_max
                    cs.reproductivec_patch[p, k] = cs.reproductivec_patch[p, k] -
                        (cf.repr_grainc_to_food_patch[p, k] + cf.repr_grainc_to_seed_patch[p, k]) * dt
                    cs.cropseedc_deficit_patch[p] = cs.cropseedc_deficit_patch[p] +
                        cf.repr_grainc_to_seed_patch[p, k] * dt
                end
                for k in repr_structure_min:repr_structure_max
                    cs.reproductivec_patch[p, k] = cs.reproductivec_patch[p, k] -
                        (cf.repr_structurec_to_cropprod_patch[p, k] + cf.repr_structurec_to_litter_patch[p, k]) * dt
                end
            end
        else
            # Matrix version: only cropseedc_deficit updates
            if ivt[p] >= npcropmin
                cs.cropseedc_deficit_patch[p] = cs.cropseedc_deficit_patch[p] -
                    cf.crop_seedc_to_leaf_patch[p] * dt
                for k in repr_grain_min:repr_grain_max
                    cs.cropseedc_deficit_patch[p] = cs.cropseedc_deficit_patch[p] +
                        cf.repr_grainc_to_seed_patch[p, k] * dt
                end
            end
        end  # !use_matrixcn

        # === Maintenance respiration fluxes from cpool ===
        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_xsmrpool_patch[p] * dt
        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.leaf_curmr_patch[p] * dt
        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.froot_curmr_patch[p] * dt

        if is_woody
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.livestem_curmr_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.livecroot_curmr_patch[p] * dt
        end
        if ivt[p] >= npcropmin  # skip 2 generic crops
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.livestem_curmr_patch[p] * dt
            for k in 1:nrepr
                cs.cpool_patch[p] = cs.cpool_patch[p] - cf.reproductive_curmr_patch[p, k] * dt
            end
        end

        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_resp_patch[p] * dt

        # RF: carbon spent on uptake respiration
        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.soilc_change_patch[p] * dt

        # === Maintenance respiration fluxes from xsmrpool ===
        cs.xsmrpool_patch[p] = cs.xsmrpool_patch[p] + cf.cpool_to_xsmrpool_patch[p] * dt
        cs.xsmrpool_patch[p] = cs.xsmrpool_patch[p] - cf.leaf_xsmr_patch[p] * dt
        cs.xsmrpool_patch[p] = cs.xsmrpool_patch[p] - cf.froot_xsmr_patch[p] * dt
        if is_woody
            cs.xsmrpool_patch[p] = cs.xsmrpool_patch[p] - cf.livestem_xsmr_patch[p] * dt
            cs.xsmrpool_patch[p] = cs.xsmrpool_patch[p] - cf.livecroot_xsmr_patch[p] * dt
        end

        # === Allocation fluxes ===
        if carbon_resp_opt == 1
            cf.cpool_to_leafc_patch[p]          = cf.cpool_to_leafc_patch[p]          - cf.cpool_to_leafc_resp_patch[p]
            cf.cpool_to_leafc_storage_patch[p]  = cf.cpool_to_leafc_storage_patch[p]  - cf.cpool_to_leafc_storage_resp_patch[p]
            cf.cpool_to_frootc_patch[p]         = cf.cpool_to_frootc_patch[p]         - cf.cpool_to_frootc_resp_patch[p]
            cf.cpool_to_frootc_storage_patch[p] = cf.cpool_to_frootc_storage_patch[p] - cf.cpool_to_frootc_storage_resp_patch[p]
        end

        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_leafc_patch[p] * dt
        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_leafc_storage_patch[p] * dt
        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_frootc_patch[p] * dt
        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_frootc_storage_patch[p] * dt

        if !use_matrixcn
            cs.leafc_patch[p]          = cs.leafc_patch[p]          + cf.cpool_to_leafc_patch[p] * dt
            cs.leafc_storage_patch[p]  = cs.leafc_storage_patch[p]  + cf.cpool_to_leafc_storage_patch[p] * dt
            cs.frootc_patch[p]         = cs.frootc_patch[p]         + cf.cpool_to_frootc_patch[p] * dt
            cs.frootc_storage_patch[p] = cs.frootc_storage_patch[p] + cf.cpool_to_frootc_storage_patch[p] * dt
        end

        if is_woody
            if carbon_resp_opt == 1
                cf.cpool_to_livecrootc_patch[p]         = cf.cpool_to_livecrootc_patch[p]         - cf.cpool_to_livecrootc_resp_patch[p]
                cf.cpool_to_livecrootc_storage_patch[p] = cf.cpool_to_livecrootc_storage_patch[p] - cf.cpool_to_livecrootc_storage_resp_patch[p]
                cf.cpool_to_livestemc_patch[p]          = cf.cpool_to_livestemc_patch[p]          - cf.cpool_to_livestemc_resp_patch[p]
                cf.cpool_to_livestemc_storage_patch[p]  = cf.cpool_to_livestemc_storage_patch[p]  - cf.cpool_to_livestemc_storage_resp_patch[p]
            end

            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_livestemc_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_livestemc_storage_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_deadstemc_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_deadstemc_storage_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_livecrootc_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_livecrootc_storage_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_deadcrootc_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_deadcrootc_storage_patch[p] * dt

            if !use_matrixcn
                cs.livestemc_patch[p]          = cs.livestemc_patch[p]          + cf.cpool_to_livestemc_patch[p] * dt
                cs.livestemc_storage_patch[p]  = cs.livestemc_storage_patch[p]  + cf.cpool_to_livestemc_storage_patch[p] * dt
                cs.deadstemc_patch[p]          = cs.deadstemc_patch[p]          + cf.cpool_to_deadstemc_patch[p] * dt
                cs.deadstemc_storage_patch[p]  = cs.deadstemc_storage_patch[p]  + cf.cpool_to_deadstemc_storage_patch[p] * dt
                cs.livecrootc_patch[p]         = cs.livecrootc_patch[p]         + cf.cpool_to_livecrootc_patch[p] * dt
                cs.livecrootc_storage_patch[p] = cs.livecrootc_storage_patch[p] + cf.cpool_to_livecrootc_storage_patch[p] * dt
                cs.deadcrootc_patch[p]         = cs.deadcrootc_patch[p]         + cf.cpool_to_deadcrootc_patch[p] * dt
                cs.deadcrootc_storage_patch[p] = cs.deadcrootc_storage_patch[p] + cf.cpool_to_deadcrootc_storage_patch[p] * dt
            end
        end

        if ivt[p] >= npcropmin  # skip 2 generic crops
            if carbon_resp_opt == 1
                cf.cpool_to_livestemc_patch[p]         = cf.cpool_to_livestemc_patch[p]         - cf.cpool_to_livestemc_resp_patch[p]
                cf.cpool_to_livestemc_storage_patch[p] = cf.cpool_to_livestemc_storage_patch[p] - cf.cpool_to_livestemc_storage_resp_patch[p]
            end

            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_livestemc_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_livestemc_storage_patch[p] * dt
            for k in 1:nrepr
                cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_reproductivec_patch[p, k] * dt
                cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_to_reproductivec_storage_patch[p, k] * dt
            end

            if !use_matrixcn
                cs.livestemc_patch[p]         = cs.livestemc_patch[p]         + cf.cpool_to_livestemc_patch[p] * dt
                cs.livestemc_storage_patch[p] = cs.livestemc_storage_patch[p] + cf.cpool_to_livestemc_storage_patch[p] * dt
                for k in 1:nrepr
                    cs.reproductivec_patch[p, k]         = cs.reproductivec_patch[p, k] +
                        cf.cpool_to_reproductivec_patch[p, k] * dt
                    cs.reproductivec_storage_patch[p, k] = cs.reproductivec_storage_patch[p, k] +
                        cf.cpool_to_reproductivec_storage_patch[p, k] * dt
                end
            end
        end

        # === Growth respiration fluxes for current growth ===
        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_leaf_gr_patch[p] * dt
        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_froot_gr_patch[p] * dt

        if is_woody
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_livestem_gr_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_deadstem_gr_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_livecroot_gr_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_deadcroot_gr_patch[p] * dt
        end
        if ivt[p] >= npcropmin  # skip 2 generic crops
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_livestem_gr_patch[p] * dt
            for k in 1:nrepr
                cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_reproductive_gr_patch[p, k] * dt
            end
        end

        # === Growth respiration for transfer growth ===
        cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] - cf.transfer_leaf_gr_patch[p] * dt
        cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] - cf.transfer_froot_gr_patch[p] * dt
        if is_woody
            cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] - cf.transfer_livestem_gr_patch[p] * dt
            cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] - cf.transfer_deadstem_gr_patch[p] * dt
            cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] - cf.transfer_livecroot_gr_patch[p] * dt
            cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] - cf.transfer_deadcroot_gr_patch[p] * dt
        end
        if ivt[p] >= npcropmin  # skip 2 generic crops
            cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] - cf.transfer_livestem_gr_patch[p] * dt
            for k in 1:nrepr
                cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] - cf.transfer_reproductive_gr_patch[p, k] * dt
            end
        end

        # === Growth respiration at time of storage ===
        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_leaf_storage_gr_patch[p] * dt
        cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_froot_storage_gr_patch[p] * dt

        if is_woody
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_livestem_storage_gr_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_deadstem_storage_gr_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_livecroot_storage_gr_patch[p] * dt
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_deadcroot_storage_gr_patch[p] * dt
        end
        if ivt[p] >= npcropmin  # skip 2 generic crops
            cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_livestem_storage_gr_patch[p] * dt
            for k in 1:nrepr
                cs.cpool_patch[p] = cs.cpool_patch[p] - cf.cpool_reproductive_storage_gr_patch[p, k] * dt
            end
        end

        # === Growth respiration stored for release during transfer growth ===
        cs.cpool_patch[p]         = cs.cpool_patch[p]         - cf.cpool_to_gresp_storage_patch[p] * dt
        cs.gresp_storage_patch[p] = cs.gresp_storage_patch[p] + cf.cpool_to_gresp_storage_patch[p] * dt

        # === Move storage pools into transfer pools ===
        if !use_matrixcn
            cs.leafc_storage_patch[p]  = cs.leafc_storage_patch[p]  - cf.leafc_storage_to_xfer_patch[p] * dt
            cs.leafc_xfer_patch[p]     = cs.leafc_xfer_patch[p]     + cf.leafc_storage_to_xfer_patch[p] * dt
            cs.frootc_storage_patch[p] = cs.frootc_storage_patch[p] - cf.frootc_storage_to_xfer_patch[p] * dt
            cs.frootc_xfer_patch[p]    = cs.frootc_xfer_patch[p]    + cf.frootc_storage_to_xfer_patch[p] * dt
        end

        if is_woody
            cs.gresp_storage_patch[p] = cs.gresp_storage_patch[p] - cf.gresp_storage_to_xfer_patch[p] * dt
            cs.gresp_xfer_patch[p]    = cs.gresp_xfer_patch[p]    + cf.gresp_storage_to_xfer_patch[p] * dt
            if !use_matrixcn
                cs.livestemc_storage_patch[p]  = cs.livestemc_storage_patch[p]  - cf.livestemc_storage_to_xfer_patch[p] * dt
                cs.livestemc_xfer_patch[p]     = cs.livestemc_xfer_patch[p]     + cf.livestemc_storage_to_xfer_patch[p] * dt
                cs.deadstemc_storage_patch[p]  = cs.deadstemc_storage_patch[p]  - cf.deadstemc_storage_to_xfer_patch[p] * dt
                cs.deadstemc_xfer_patch[p]     = cs.deadstemc_xfer_patch[p]     + cf.deadstemc_storage_to_xfer_patch[p] * dt
                cs.livecrootc_storage_patch[p] = cs.livecrootc_storage_patch[p] - cf.livecrootc_storage_to_xfer_patch[p] * dt
                cs.livecrootc_xfer_patch[p]    = cs.livecrootc_xfer_patch[p]    + cf.livecrootc_storage_to_xfer_patch[p] * dt
                cs.deadcrootc_storage_patch[p] = cs.deadcrootc_storage_patch[p] - cf.deadcrootc_storage_to_xfer_patch[p] * dt
                cs.deadcrootc_xfer_patch[p]    = cs.deadcrootc_xfer_patch[p]    + cf.deadcrootc_storage_to_xfer_patch[p] * dt
            end
        end

        if ivt[p] >= npcropmin  # skip 2 generic crops
            if !use_matrixcn
                # lines here for consistency; the transfer terms are zero
                cs.livestemc_storage_patch[p] = cs.livestemc_storage_patch[p] - cf.livestemc_storage_to_xfer_patch[p] * dt
                cs.livestemc_xfer_patch[p]    = cs.livestemc_xfer_patch[p]    + cf.livestemc_storage_to_xfer_patch[p] * dt
                for k in 1:nrepr
                    cs.reproductivec_storage_patch[p, k] = cs.reproductivec_storage_patch[p, k] -
                        cf.reproductivec_storage_to_xfer_patch[p, k] * dt
                    cs.reproductivec_xfer_patch[p, k]    = cs.reproductivec_xfer_patch[p, k] +
                        cf.reproductivec_storage_to_xfer_patch[p, k] * dt
                end
            end
        end

        # === Crop-specific xsmrpool and harvest logic ===
        if ivt[p] >= npcropmin  # skip 2 generic crops
            cs.xsmrpool_patch[p] = cs.xsmrpool_patch[p] - cf.livestem_xsmr_patch[p] * dt
            for k in 1:nrepr
                cs.xsmrpool_patch[p] = cs.xsmrpool_patch[p] - cf.reproductive_xsmr_patch[p, k] * dt
            end

            if harvdate[p] < 999  # beginning at harvest, send to atm
                if !dribble_crophrv_xsmrpool_2atm
                    cf.xsmrpool_to_atm_patch[p] = cf.xsmrpool_to_atm_patch[p] + cs.xsmrpool_patch[p] / dt
                    cf.xsmrpool_to_atm_patch[p] = cf.xsmrpool_to_atm_patch[p] + cs.cpool_patch[p] / dt
                    if !use_matrixcn
                        cf.xsmrpool_to_atm_patch[p] = cf.xsmrpool_to_atm_patch[p] + cs.frootc_patch[p] / dt
                    end
                    # Note: matrix_update_phc path omitted (requires full matrix infrastructure)
                else
                    cs.xsmrpool_loss_patch[p] = cs.xsmrpool_loss_patch[p] +
                        cs.xsmrpool_patch[p] +
                        cs.cpool_patch[p]
                    if !use_matrixcn
                        cs.xsmrpool_loss_patch[p] = cs.xsmrpool_loss_patch[p] + cs.frootc_patch[p]
                    end
                    # Note: matrix_update_phc path omitted (requires full matrix infrastructure)
                end
                if !use_matrixcn
                    cs.frootc_patch[p] = zero(T)
                end
                cs.xsmrpool_patch[p] = zero(T)
                cs.cpool_patch[p]    = zero(T)
            end

            # Slowly release xsmrpool to atmosphere
            if dribble_crophrv_xsmrpool_2atm
                cf.xsmrpool_to_atm_patch[p] = cs.xsmrpool_loss_patch[p] * kprod05
                cs.xsmrpool_loss_patch[p]   = cs.xsmrpool_loss_patch[p] - cf.xsmrpool_to_atm_patch[p] * dt
            end
        end
    end
end

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

Ported from `CStateUpdate1` in `CNCStateUpdate1Mod.F90`. Runs as two
KernelAbstractions kernels (one per soil column, one per patch); see the
kernel comments above for the GPU mapping. The CPU path is byte-identical.

Note: The FATES litter flux callback (`clm_fates%UpdateCLitterfluxes`) is
skipped in this port — FATES columns are handled by zeroing
`decomp_cpools_sourcesink_col` via `col_is_fates`. The matrix-CN code paths
(`use_matrixcn`, `use_soil_matrixcn`) for soil/veg matrix solutions are
ported as-is, but the `matrix_update_phc` helper is replaced by a no-op
since it depends on the full matrix infrastructure. `patch_column` is accepted
for call-site compatibility but unused (the patch loop touches no column index).
"""
function c_state_update1!(cs_veg::CNVegCarbonStateData,
                           cf_veg::CNVegCarbonFluxData,
                           cf_soil::SoilBiogeochemCarbonFluxData;
                           mask_soilc::AbstractVector{Bool},
                           mask_soilp::AbstractVector{Bool},
                           bounds_col::UnitRange{Int},
                           bounds_patch::UnitRange{Int},
                           patch_column::AbstractVector{<:Integer},
                           ivt::AbstractVector{<:Integer},
                           woody::AbstractVector,
                           cascade_donor_pool::AbstractVector{<:Integer},
                           cascade_receiver_pool::AbstractVector{<:Integer},
                           harvdate::AbstractVector{<:Integer},
                           col_is_fates::AbstractVector{Bool},
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
                           dt::Real)

    # --- Column loop: soil decomposition input fluxes (one thread per column) ---
    _launch!(_csu1_col_kernel!, mask_soilc, col_is_fates,
        cf_veg.phenology_c_to_litr_c_col, cf_soil.decomp_cpools_sourcesink_col,
        cf_soil.decomp_cascade_hr_vr_col, cf_soil.decomp_cascade_ctransfer_vr_col,
        cascade_donor_pool, cascade_receiver_pool,
        nlevdecomp, ndecomp_cascade_transitions, i_litr_min, i_litr_max, i_cwd,
        use_soil_matrixcn, dt)

    # --- Patch loop: vegetation C state updates (one thread per patch) ---
    cs = _csu1_cs(cs_veg)
    cf = _csu1_cf(cf_veg)
    _launch!(_csu1_patch_kernel!, mask_soilp, ivt, woody, harvdate, cs, cf,
        npcropmin, nrepr, repr_grain_min, repr_grain_max,
        repr_structure_min, repr_structure_max,
        use_matrixcn, carbon_resp_opt, dribble_crophrv_xsmrpool_2atm, dt)

    return nothing
end
