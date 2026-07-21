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
#
# GPU kernelization (Phase B C/N state-update cascade, nitrogen analog of
# c_state_update1.jl):
#   * n_state_update_dyn_patch!: small per-column j/i RMW + per-gridcell loops ->
#     `_nsudp_col_kernel!` (one thread per column; the internal j/i decomp loop runs
#     sequentially in-thread, each thread owns its column) and `_nsudp_grc_kernel!`
#     (one thread per gridcell). Only ~6 array fields, so they are passed loose
#     (no device-view bundle).
#   * n_state_update1!: the column input-flux loop -> `_nsu1_col_kernel!`; the
#     field-heavy patch loop -> `_nsu1_patch_kernel!`, with the ~24 N-state and
#     ~48 N-flux arrays it touches grouped into two immutable @adapt_structure
#     device-view bundles (`_NSU1NS{V,M}` / `_NSU1NF{V,M}`) so the launch stays under
#     Metal's ~31-arg limit. Every write is to the loop's own index (no patch->column
#     scatter). Float64 literals are eltype-converted (`zero(T)`) so the kernels carry
#     no Float64 on a Float32-only backend; on Float64 this is byte-identical.
# ==========================================================================

# ---------------------------------------------------------------------------
# n_state_update_dyn_patch! — Dynamic patch nitrogen state update
# ---------------------------------------------------------------------------

# --- Kernel: column-level dyn-patch decomposition input (one thread per column) ---
@kernel function _nsudp_col_kernel!(@Const(mask_soilc_with_inactive),
        decomp_npools_vr_col, matrix_Ninput_col, @Const(dwt_frootn_to_litr_n_col),
        @Const(dwt_livecrootn_to_cwdn_col), @Const(dwt_deadcrootn_to_cwdn_col),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        use_soil_matrixcn::Bool, dt)
    c = @index(Global)
    @inbounds if mask_soilc_with_inactive[c]
        if use_soil_matrixcn
            # Matrix mode: accumulate dwt→litter/CWD N into the persistent B-input.
            for j in 1:nlevdecomp
                for i in i_litr_min:i_litr_max
                    matrix_Ninput_col[c, j + (i - 1) * nlevdecomp] =
                        matrix_Ninput_col[c, j + (i - 1) * nlevdecomp] +
                        dwt_frootn_to_litr_n_col[c, j, i] * dt
                end
                matrix_Ninput_col[c, j + (i_cwd - 1) * nlevdecomp] =
                    matrix_Ninput_col[c, j + (i_cwd - 1) * nlevdecomp] +
                    (dwt_livecrootn_to_cwdn_col[c, j] +
                     dwt_deadcrootn_to_cwdn_col[c, j]) * dt
            end
        else
            for j in 1:nlevdecomp
                for i in i_litr_min:i_litr_max
                    decomp_npools_vr_col[c, j, i] =
                        decomp_npools_vr_col[c, j, i] +
                        dwt_frootn_to_litr_n_col[c, j, i] * dt
                end
                decomp_npools_vr_col[c, j, i_cwd] =
                    decomp_npools_vr_col[c, j, i_cwd] +
                    (dwt_livecrootn_to_cwdn_col[c, j] +
                     dwt_deadcrootn_to_cwdn_col[c, j]) * dt
            end
        end
    end
end

# --- Kernel: gridcell-level seed-N update (one thread per gridcell) ---
@kernel function _nsudp_grc_kernel!(seedn_grc,
        @Const(dwt_seedn_to_leaf_grc), @Const(dwt_seedn_to_deadstem_grc), dt)
    g = @index(Global)
    @inbounds begin
        seedn_grc[g] = seedn_grc[g] - dwt_seedn_to_leaf_grc[g] * dt
        seedn_grc[g] = seedn_grc[g] - dwt_seedn_to_deadstem_grc[g] * dt
    end
end

"""
    n_state_update_dyn_patch!(ns_veg, nf_veg, ns_soil;
        mask_soilc_with_inactive, bounds_col, bounds_grc,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd, dt)

Update nitrogen states based on fluxes from dyn_cnbal_patch.

Ported from `NStateUpdateDynPatch` in `CNNStateUpdate1Mod.F90`. Runs as two
KernelAbstractions kernels (one per soil column, one per gridcell); the CPU
path is byte-identical.
"""
function n_state_update_dyn_patch!(ns_veg::CNVegNitrogenStateData,
                                    nf_veg::CNVegNitrogenFluxData,
                                    ns_soil::SoilBiogeochemNitrogenStateData;
                                    mask_soilc_with_inactive::AbstractVector{Bool},
                                    bounds_col::UnitRange{Int},
                                    bounds_grc::UnitRange{Int},
                                    nlevdecomp::Int,
                                    i_litr_min::Int,
                                    i_litr_max::Int,
                                    i_cwd::Int,
                                    use_soil_matrixcn::Bool = false,
                                    nf_soil::Union{SoilBiogeochemNitrogenFluxData, Nothing} = nothing,
                                    dt::Real)

    do_matrix = use_soil_matrixcn && nf_soil !== nothing &&
                size(nf_soil.matrix_Ninput_col, 2) > 0
    mNin = do_matrix ? nf_soil.matrix_Ninput_col :
           similar(ns_soil.decomp_npools_vr_col, size(ns_soil.decomp_npools_vr_col, 1), 1)
    _launch!(_nsudp_col_kernel!, mask_soilc_with_inactive,
        ns_soil.decomp_npools_vr_col, mNin, nf_veg.dwt_frootn_to_litr_n_col,
        nf_veg.dwt_livecrootn_to_cwdn_col, nf_veg.dwt_deadcrootn_to_cwdn_col,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd, do_matrix, dt)

    _launch!(_nsudp_grc_kernel!, ns_veg.seedn_grc,
        nf_veg.dwt_seedn_to_leaf_grc, nf_veg.dwt_seedn_to_deadstem_grc, dt)

    return nothing
end

# ---------------------------------------------------------------------------
# n_state_update1! — Main prognostic N state update
#
# GPU kernelization: column loop -> `_nsu1_col_kernel!` (one thread per soil
# column); patch loop -> `_nsu1_patch_kernel!` (one thread per patch, all writes
# to the patch's own index). The N-state and N-flux arrays the patch loop touches
# are grouped into the `_NSU1NS{V,M}` / `_NSU1NF{V,M}` device-view bundles whose
# field names mirror the state structs so the kernel body reads verbatim.
# ---------------------------------------------------------------------------

# --- Device-view bundle: CNVeg nitrogen STATE arrays the patch loop reads/writes ---
# `V` = the 1D patch vectors, `M` = the 2D (patch x nrepr) reproductive matrices.
Base.@kwdef struct _NSU1NS{V,M}
    leafn_patch::V; leafn_xfer_patch::V; leafn_storage_patch::V
    frootn_patch::V; frootn_xfer_patch::V; frootn_storage_patch::V
    livestemn_patch::V; livestemn_xfer_patch::V; livestemn_storage_patch::V
    deadstemn_patch::V; deadstemn_xfer_patch::V; deadstemn_storage_patch::V
    livecrootn_patch::V; livecrootn_xfer_patch::V; livecrootn_storage_patch::V
    deadcrootn_patch::V; deadcrootn_xfer_patch::V; deadcrootn_storage_patch::V
    retransn_patch::V; npool_patch::V; cropseedn_deficit_patch::V
    reproductiven_patch::M; reproductiven_storage_patch::M; reproductiven_xfer_patch::M
end
Adapt.@adapt_structure _NSU1NS

_nsu1_ns(ns) = _NSU1NS(;
    leafn_patch = ns.leafn_patch, leafn_xfer_patch = ns.leafn_xfer_patch,
    leafn_storage_patch = ns.leafn_storage_patch,
    frootn_patch = ns.frootn_patch, frootn_xfer_patch = ns.frootn_xfer_patch,
    frootn_storage_patch = ns.frootn_storage_patch,
    livestemn_patch = ns.livestemn_patch, livestemn_xfer_patch = ns.livestemn_xfer_patch,
    livestemn_storage_patch = ns.livestemn_storage_patch,
    deadstemn_patch = ns.deadstemn_patch, deadstemn_xfer_patch = ns.deadstemn_xfer_patch,
    deadstemn_storage_patch = ns.deadstemn_storage_patch,
    livecrootn_patch = ns.livecrootn_patch, livecrootn_xfer_patch = ns.livecrootn_xfer_patch,
    livecrootn_storage_patch = ns.livecrootn_storage_patch,
    deadcrootn_patch = ns.deadcrootn_patch, deadcrootn_xfer_patch = ns.deadcrootn_xfer_patch,
    deadcrootn_storage_patch = ns.deadcrootn_storage_patch,
    retransn_patch = ns.retransn_patch, npool_patch = ns.npool_patch,
    cropseedn_deficit_patch = ns.cropseedn_deficit_patch,
    reproductiven_patch = ns.reproductiven_patch,
    reproductiven_storage_patch = ns.reproductiven_storage_patch,
    reproductiven_xfer_patch = ns.reproductiven_xfer_patch)

# --- Device-view bundle: CNVeg nitrogen FLUX arrays the patch loop reads/writes ---
Base.@kwdef struct _NSU1NF{V,M}
    # 1D patch vectors
    leafn_xfer_to_leafn_patch::V; frootn_xfer_to_frootn_patch::V
    livestemn_xfer_to_livestemn_patch::V; deadstemn_xfer_to_deadstemn_patch::V
    livecrootn_xfer_to_livecrootn_patch::V; deadcrootn_xfer_to_deadcrootn_patch::V
    leafn_to_litter_patch::V; frootn_to_litter_patch::V; leafn_to_retransn_patch::V
    livestemn_to_deadstemn_patch::V; livestemn_to_retransn_patch::V
    livecrootn_to_deadcrootn_patch::V; livecrootn_to_retransn_patch::V
    free_retransn_to_npool_patch::V; frootn_to_retransn_patch::V
    livestemn_to_litter_patch::V; livestemn_to_biofueln_patch::V
    livestemn_to_removedresiduen_patch::V
    leafn_to_biofueln_patch::V; leafn_to_removedresiduen_patch::V
    crop_seedn_to_leaf_patch::V
    sminn_to_npool_patch::V; retransn_to_npool_patch::V
    npool_to_leafn_patch::V; npool_to_leafn_storage_patch::V
    npool_to_frootn_patch::V; npool_to_frootn_storage_patch::V
    npool_to_livestemn_patch::V; npool_to_livestemn_storage_patch::V
    npool_to_deadstemn_patch::V; npool_to_deadstemn_storage_patch::V
    npool_to_livecrootn_patch::V; npool_to_livecrootn_storage_patch::V
    npool_to_deadcrootn_patch::V; npool_to_deadcrootn_storage_patch::V
    leafn_storage_to_xfer_patch::V; frootn_storage_to_xfer_patch::V
    livestemn_storage_to_xfer_patch::V; deadstemn_storage_to_xfer_patch::V
    livecrootn_storage_to_xfer_patch::V; deadcrootn_storage_to_xfer_patch::V
    # 2D (patch x nrepr) reproductive matrices
    reproductiven_xfer_to_reproductiven_patch::M
    repr_grainn_to_food_patch::M; repr_grainn_to_seed_patch::M
    repr_structuren_to_cropprod_patch::M; repr_structuren_to_litter_patch::M
    npool_to_reproductiven_patch::M; npool_to_reproductiven_storage_patch::M
    reproductiven_storage_to_xfer_patch::M
end
Adapt.@adapt_structure _NSU1NF

_nsu1_nf(nf) = _NSU1NF(;
    leafn_xfer_to_leafn_patch = nf.leafn_xfer_to_leafn_patch,
    frootn_xfer_to_frootn_patch = nf.frootn_xfer_to_frootn_patch,
    livestemn_xfer_to_livestemn_patch = nf.livestemn_xfer_to_livestemn_patch,
    deadstemn_xfer_to_deadstemn_patch = nf.deadstemn_xfer_to_deadstemn_patch,
    livecrootn_xfer_to_livecrootn_patch = nf.livecrootn_xfer_to_livecrootn_patch,
    deadcrootn_xfer_to_deadcrootn_patch = nf.deadcrootn_xfer_to_deadcrootn_patch,
    leafn_to_litter_patch = nf.leafn_to_litter_patch,
    frootn_to_litter_patch = nf.frootn_to_litter_patch,
    leafn_to_retransn_patch = nf.leafn_to_retransn_patch,
    livestemn_to_deadstemn_patch = nf.livestemn_to_deadstemn_patch,
    livestemn_to_retransn_patch = nf.livestemn_to_retransn_patch,
    livecrootn_to_deadcrootn_patch = nf.livecrootn_to_deadcrootn_patch,
    livecrootn_to_retransn_patch = nf.livecrootn_to_retransn_patch,
    free_retransn_to_npool_patch = nf.free_retransn_to_npool_patch,
    frootn_to_retransn_patch = nf.frootn_to_retransn_patch,
    livestemn_to_litter_patch = nf.livestemn_to_litter_patch,
    livestemn_to_biofueln_patch = nf.livestemn_to_biofueln_patch,
    livestemn_to_removedresiduen_patch = nf.livestemn_to_removedresiduen_patch,
    leafn_to_biofueln_patch = nf.leafn_to_biofueln_patch,
    leafn_to_removedresiduen_patch = nf.leafn_to_removedresiduen_patch,
    crop_seedn_to_leaf_patch = nf.crop_seedn_to_leaf_patch,
    sminn_to_npool_patch = nf.sminn_to_npool_patch,
    retransn_to_npool_patch = nf.retransn_to_npool_patch,
    npool_to_leafn_patch = nf.npool_to_leafn_patch,
    npool_to_leafn_storage_patch = nf.npool_to_leafn_storage_patch,
    npool_to_frootn_patch = nf.npool_to_frootn_patch,
    npool_to_frootn_storage_patch = nf.npool_to_frootn_storage_patch,
    npool_to_livestemn_patch = nf.npool_to_livestemn_patch,
    npool_to_livestemn_storage_patch = nf.npool_to_livestemn_storage_patch,
    npool_to_deadstemn_patch = nf.npool_to_deadstemn_patch,
    npool_to_deadstemn_storage_patch = nf.npool_to_deadstemn_storage_patch,
    npool_to_livecrootn_patch = nf.npool_to_livecrootn_patch,
    npool_to_livecrootn_storage_patch = nf.npool_to_livecrootn_storage_patch,
    npool_to_deadcrootn_patch = nf.npool_to_deadcrootn_patch,
    npool_to_deadcrootn_storage_patch = nf.npool_to_deadcrootn_storage_patch,
    leafn_storage_to_xfer_patch = nf.leafn_storage_to_xfer_patch,
    frootn_storage_to_xfer_patch = nf.frootn_storage_to_xfer_patch,
    livestemn_storage_to_xfer_patch = nf.livestemn_storage_to_xfer_patch,
    deadstemn_storage_to_xfer_patch = nf.deadstemn_storage_to_xfer_patch,
    livecrootn_storage_to_xfer_patch = nf.livecrootn_storage_to_xfer_patch,
    deadcrootn_storage_to_xfer_patch = nf.deadcrootn_storage_to_xfer_patch,
    reproductiven_xfer_to_reproductiven_patch = nf.reproductiven_xfer_to_reproductiven_patch,
    repr_grainn_to_food_patch = nf.repr_grainn_to_food_patch,
    repr_grainn_to_seed_patch = nf.repr_grainn_to_seed_patch,
    repr_structuren_to_cropprod_patch = nf.repr_structuren_to_cropprod_patch,
    repr_structuren_to_litter_patch = nf.repr_structuren_to_litter_patch,
    npool_to_reproductiven_patch = nf.npool_to_reproductiven_patch,
    npool_to_reproductiven_storage_patch = nf.npool_to_reproductiven_storage_patch,
    reproductiven_storage_to_xfer_patch = nf.reproductiven_storage_to_xfer_patch)

# --- Kernel: column-level soil decomposition input fluxes (one thread per column) ---
@kernel function _nsu1_col_kernel!(@Const(mask_soilc), @Const(col_is_fates),
        @Const(phenology_n_to_litr_n_col), decomp_npools_sourcesink_col,
        @Const(decomp_cascade_ntransfer_vr_col), @Const(decomp_cascade_sminn_flux_vr_col),
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
                    for i in i_litr_min:i_litr_max
                        decomp_npools_sourcesink_col[c, j, i] =
                            phenology_n_to_litr_n_col[c, j, i] * dt
                    end
                    # NOTE(wjs, 2017-01-02) This used to be set to a non-zero value, but the
                    # terms have been moved to NStateUpdateDynPatch. Explicitly zeroed here.
                    decomp_npools_sourcesink_col[c, j, i_cwd] = zero(T)
                end
            end
        end

        # Decomposition-cascade N transfers into the organic decomp pools.
        # Ported from SoilBiogeochemNStateUpdate1Mod.F90:118-160 (the N analogue of
        # the carbon block in `_csu1_col_kernel!`). The MINERAL-N side of the cascade
        # (gross mineralization + immobilization) is applied separately in
        # `soilbiogeochem_n_state_update1!`. WITHOUT this block the net-mineralized N
        # was ADDED to the smin pool but never REMOVED from the organic pools, so a
        # column created N equal to the net mineralization each step (~3e-4 gN/m2/step
        # in the Bow use_cn summer; the whole errnb residual documented in
        # docs/CN_BALANCE_STATUS.md). Runs for all soil columns (fates handled apart),
        # matching the carbon kernel. Sequential k within each column thread => the
        # donor/receiver += into shared pool indices is race-free.
        if !use_soil_matrixcn
            for j in 1:nlevdecomp
                for k in 1:ndecomp_cascade_transitions
                    # N loss from the donor decomposing pool
                    decomp_npools_sourcesink_col[c, j, cascade_donor_pool[k]] =
                        decomp_npools_sourcesink_col[c, j, cascade_donor_pool[k]] -
                        decomp_cascade_ntransfer_vr_col[c, j, k] * dt
                    if cascade_receiver_pool[k] != 0  # skip terminal transitions
                        # N gain to the receiver pool: transferred N + cascade mineral-N flux
                        decomp_npools_sourcesink_col[c, j, cascade_receiver_pool[k]] =
                            decomp_npools_sourcesink_col[c, j, cascade_receiver_pool[k]] +
                            (decomp_cascade_ntransfer_vr_col[c, j, k] +
                             decomp_cascade_sminn_flux_vr_col[c, j, k]) * dt
                    else  # terminal transitions: debit donor by the cascade mineral-N flux
                        decomp_npools_sourcesink_col[c, j, cascade_donor_pool[k]] =
                            decomp_npools_sourcesink_col[c, j, cascade_donor_pool[k]] -
                            decomp_cascade_sminn_flux_vr_col[c, j, k] * dt
                    end
                end
            end
        end
    end
end

# --- Kernel: patch-level vegetation N state updates (one thread per patch) ---
# `ns` / `nf` are the _NSU1NS / _NSU1NF device-view bundles; field names mirror the
# state structs so the body below is the host loop verbatim (ns_veg.->ns., nf_veg.->nf.).
@kernel function _nsu1_patch_kernel!(@Const(mask_soilp), @Const(ivt), @Const(woody),
        ns, nf,
        npcropmin::Int, nrepr::Int, repr_grain_min::Int, repr_grain_max::Int,
        repr_structure_min::Int, repr_structure_max::Int,
        use_matrixcn::Bool, use_fun::Bool, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p] && ivt[p] >= 1  # skip masked + bare ground (PFT index 0)
        # OFF-BY-ONE (fixed) — the nitrogen mirror of the C bug; see the long note in
        # c_state_update1.jl. `ivt[p]` is the raw 0-based Fortran PFT index; pftcon
        # rows are 1-based, so the PFT is `ivt[p] + 1`. Indexing with the raw value read
        # the previous PFT's woody flag and made every TREE non-woody, skipping the
        # woody npool debits (livestem/livecroot/deadstem/deadcroot allocation) while
        # the demand/allocation code still counted them.
        is_woody = woody[ivt[p] + 1] == one(eltype(woody))

        # === Phenology: transfer growth fluxes ===
        if !use_matrixcn
            ns.leafn_patch[p]       = ns.leafn_patch[p]       + nf.leafn_xfer_to_leafn_patch[p] * dt
            ns.leafn_xfer_patch[p]  = ns.leafn_xfer_patch[p]  - nf.leafn_xfer_to_leafn_patch[p] * dt
            ns.frootn_patch[p]      = ns.frootn_patch[p]      + nf.frootn_xfer_to_frootn_patch[p] * dt
            ns.frootn_xfer_patch[p] = ns.frootn_xfer_patch[p] - nf.frootn_xfer_to_frootn_patch[p] * dt

            if is_woody
                ns.livestemn_patch[p]       = ns.livestemn_patch[p]       + nf.livestemn_xfer_to_livestemn_patch[p] * dt
                ns.livestemn_xfer_patch[p]  = ns.livestemn_xfer_patch[p]  - nf.livestemn_xfer_to_livestemn_patch[p] * dt
                ns.deadstemn_patch[p]       = ns.deadstemn_patch[p]       + nf.deadstemn_xfer_to_deadstemn_patch[p] * dt
                ns.deadstemn_xfer_patch[p]  = ns.deadstemn_xfer_patch[p]  - nf.deadstemn_xfer_to_deadstemn_patch[p] * dt
                ns.livecrootn_patch[p]      = ns.livecrootn_patch[p]      + nf.livecrootn_xfer_to_livecrootn_patch[p] * dt
                ns.livecrootn_xfer_patch[p] = ns.livecrootn_xfer_patch[p] - nf.livecrootn_xfer_to_livecrootn_patch[p] * dt
                ns.deadcrootn_patch[p]      = ns.deadcrootn_patch[p]      + nf.deadcrootn_xfer_to_deadcrootn_patch[p] * dt
                ns.deadcrootn_xfer_patch[p] = ns.deadcrootn_xfer_patch[p] - nf.deadcrootn_xfer_to_deadcrootn_patch[p] * dt
            end

            if ivt[p] >= npcropmin  # skip 2 generic crops
                # lines here for consistency; the transfer terms are zero
                ns.livestemn_patch[p]      = ns.livestemn_patch[p]      + nf.livestemn_xfer_to_livestemn_patch[p] * dt
                ns.livestemn_xfer_patch[p] = ns.livestemn_xfer_patch[p] - nf.livestemn_xfer_to_livestemn_patch[p] * dt
                for k in 1:nrepr
                    ns.reproductiven_patch[p, k] = ns.reproductiven_patch[p, k] +
                        nf.reproductiven_xfer_to_reproductiven_patch[p, k] * dt
                    ns.reproductiven_xfer_patch[p, k] = ns.reproductiven_xfer_patch[p, k] -
                        nf.reproductiven_xfer_to_reproductiven_patch[p, k] * dt
                end
            end

            # phenology: litterfall and retranslocation fluxes
            ns.leafn_patch[p]    = ns.leafn_patch[p]    - nf.leafn_to_litter_patch[p] * dt
            ns.frootn_patch[p]   = ns.frootn_patch[p]   - nf.frootn_to_litter_patch[p] * dt
            ns.leafn_patch[p]    = ns.leafn_patch[p]    - nf.leafn_to_retransn_patch[p] * dt
            ns.retransn_patch[p] = ns.retransn_patch[p] + nf.leafn_to_retransn_patch[p] * dt
        end  # !use_matrixcn

        # live wood turnover and retranslocation fluxes
        if is_woody
            if !use_matrixcn
                ns.livestemn_patch[p]  = ns.livestemn_patch[p]  - nf.livestemn_to_deadstemn_patch[p] * dt
                ns.deadstemn_patch[p]  = ns.deadstemn_patch[p]  + nf.livestemn_to_deadstemn_patch[p] * dt
                ns.livestemn_patch[p]  = ns.livestemn_patch[p]  - nf.livestemn_to_retransn_patch[p] * dt
                ns.retransn_patch[p]   = ns.retransn_patch[p]   + nf.livestemn_to_retransn_patch[p] * dt
                ns.livecrootn_patch[p] = ns.livecrootn_patch[p] - nf.livecrootn_to_deadcrootn_patch[p] * dt
                ns.deadcrootn_patch[p] = ns.deadcrootn_patch[p] + nf.livecrootn_to_deadcrootn_patch[p] * dt
                ns.livecrootn_patch[p] = ns.livecrootn_patch[p] - nf.livecrootn_to_retransn_patch[p] * dt
                ns.retransn_patch[p]   = ns.retransn_patch[p]   + nf.livecrootn_to_retransn_patch[p] * dt
                # WW change logic so livestem_retrans goes to npool (via free_retrans flux)
                if use_fun
                    nf.free_retransn_to_npool_patch[p] = nf.free_retransn_to_npool_patch[p] + nf.livestemn_to_retransn_patch[p]
                    nf.free_retransn_to_npool_patch[p] = nf.free_retransn_to_npool_patch[p] + nf.livecrootn_to_retransn_patch[p]
                end
            end  # !use_matrixcn
        end

        if ivt[p] >= npcropmin  # Beth adds retrans from froot
            if !use_matrixcn
                ns.frootn_patch[p]    = ns.frootn_patch[p]    - nf.frootn_to_retransn_patch[p] * dt
                ns.retransn_patch[p]  = ns.retransn_patch[p]  + nf.frootn_to_retransn_patch[p] * dt
                ns.livestemn_patch[p] = ns.livestemn_patch[p] - nf.livestemn_to_litter_patch[p] * dt
                ns.livestemn_patch[p] = ns.livestemn_patch[p] -
                    (nf.livestemn_to_biofueln_patch[p] + nf.livestemn_to_removedresiduen_patch[p]) * dt
                ns.leafn_patch[p]     = ns.leafn_patch[p] -
                    (nf.leafn_to_biofueln_patch[p] + nf.leafn_to_removedresiduen_patch[p]) * dt
                ns.livestemn_patch[p] = ns.livestemn_patch[p] - nf.livestemn_to_retransn_patch[p] * dt
                ns.retransn_patch[p]  = ns.retransn_patch[p]  + nf.livestemn_to_retransn_patch[p] * dt
                for k in repr_grain_min:repr_grain_max
                    ns.reproductiven_patch[p, k] = ns.reproductiven_patch[p, k] -
                        (nf.repr_grainn_to_food_patch[p, k] + nf.repr_grainn_to_seed_patch[p, k]) * dt
                end
                for k in repr_structure_min:repr_structure_max
                    ns.reproductiven_patch[p, k] = ns.reproductiven_patch[p, k] -
                        (nf.repr_structuren_to_cropprod_patch[p, k] + nf.repr_structuren_to_litter_patch[p, k]) * dt
                end
            end  # !use_matrixcn
            ns.cropseedn_deficit_patch[p] = ns.cropseedn_deficit_patch[p] -
                nf.crop_seedn_to_leaf_patch[p] * dt
            for k in repr_grain_min:repr_grain_max
                ns.cropseedn_deficit_patch[p] = ns.cropseedn_deficit_patch[p] +
                    nf.repr_grainn_to_seed_patch[p, k] * dt
            end
        end

        # uptake from soil mineral N pool
        ns.npool_patch[p] = ns.npool_patch[p] + nf.sminn_to_npool_patch[p] * dt

        # deployment from retranslocation pool
        ns.npool_patch[p] = ns.npool_patch[p] + nf.retransn_to_npool_patch[p] * dt

        ns.npool_patch[p] = ns.npool_patch[p] + nf.free_retransn_to_npool_patch[p] * dt

        # allocation fluxes
        ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_leafn_patch[p] * dt
        ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_leafn_storage_patch[p] * dt
        ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_frootn_patch[p] * dt
        ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_frootn_storage_patch[p] * dt

        if !use_matrixcn
            ns.retransn_patch[p]       = ns.retransn_patch[p]       - nf.retransn_to_npool_patch[p] * dt
            ns.retransn_patch[p]       = ns.retransn_patch[p]       - nf.free_retransn_to_npool_patch[p] * dt
            ns.leafn_patch[p]          = ns.leafn_patch[p]          + nf.npool_to_leafn_patch[p] * dt
            ns.leafn_storage_patch[p]  = ns.leafn_storage_patch[p]  + nf.npool_to_leafn_storage_patch[p] * dt
            ns.frootn_patch[p]         = ns.frootn_patch[p]         + nf.npool_to_frootn_patch[p] * dt
            ns.frootn_storage_patch[p] = ns.frootn_storage_patch[p] + nf.npool_to_frootn_storage_patch[p] * dt
        end

        if is_woody
            ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_livestemn_patch[p] * dt
            ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_livestemn_storage_patch[p] * dt
            ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_deadstemn_patch[p] * dt
            ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_deadstemn_storage_patch[p] * dt
            ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_livecrootn_patch[p] * dt
            ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_livecrootn_storage_patch[p] * dt
            ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_deadcrootn_patch[p] * dt
            ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_deadcrootn_storage_patch[p] * dt

            if !use_matrixcn
                ns.livestemn_patch[p]          = ns.livestemn_patch[p]          + nf.npool_to_livestemn_patch[p] * dt
                ns.livestemn_storage_patch[p]  = ns.livestemn_storage_patch[p]  + nf.npool_to_livestemn_storage_patch[p] * dt
                ns.deadstemn_patch[p]          = ns.deadstemn_patch[p]          + nf.npool_to_deadstemn_patch[p] * dt
                ns.deadstemn_storage_patch[p]  = ns.deadstemn_storage_patch[p]  + nf.npool_to_deadstemn_storage_patch[p] * dt
                ns.livecrootn_patch[p]         = ns.livecrootn_patch[p]         + nf.npool_to_livecrootn_patch[p] * dt
                ns.livecrootn_storage_patch[p] = ns.livecrootn_storage_patch[p] + nf.npool_to_livecrootn_storage_patch[p] * dt
                ns.deadcrootn_patch[p]         = ns.deadcrootn_patch[p]         + nf.npool_to_deadcrootn_patch[p] * dt
                ns.deadcrootn_storage_patch[p] = ns.deadcrootn_storage_patch[p] + nf.npool_to_deadcrootn_storage_patch[p] * dt
            end
        end

        if ivt[p] >= npcropmin  # skip 2 generic crops
            ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_livestemn_patch[p] * dt
            ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_livestemn_storage_patch[p] * dt
            for k in 1:nrepr
                ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_reproductiven_patch[p, k] * dt
                ns.npool_patch[p] = ns.npool_patch[p] - nf.npool_to_reproductiven_storage_patch[p, k] * dt
            end

            if !use_matrixcn
                ns.livestemn_patch[p]         = ns.livestemn_patch[p]         + nf.npool_to_livestemn_patch[p] * dt
                ns.livestemn_storage_patch[p] = ns.livestemn_storage_patch[p] + nf.npool_to_livestemn_storage_patch[p] * dt
                for k in 1:nrepr
                    ns.reproductiven_patch[p, k] = ns.reproductiven_patch[p, k] +
                        nf.npool_to_reproductiven_patch[p, k] * dt
                    ns.reproductiven_storage_patch[p, k] = ns.reproductiven_storage_patch[p, k] +
                        nf.npool_to_reproductiven_storage_patch[p, k] * dt
                end
            end
        end

        # move storage pools into transfer pools
        if !use_matrixcn
            ns.leafn_storage_patch[p]  = ns.leafn_storage_patch[p]  - nf.leafn_storage_to_xfer_patch[p] * dt
            ns.leafn_xfer_patch[p]     = ns.leafn_xfer_patch[p]     + nf.leafn_storage_to_xfer_patch[p] * dt
            ns.frootn_storage_patch[p] = ns.frootn_storage_patch[p] - nf.frootn_storage_to_xfer_patch[p] * dt
            ns.frootn_xfer_patch[p]    = ns.frootn_xfer_patch[p]    + nf.frootn_storage_to_xfer_patch[p] * dt

            if is_woody
                ns.livestemn_storage_patch[p]  = ns.livestemn_storage_patch[p]  - nf.livestemn_storage_to_xfer_patch[p] * dt
                ns.livestemn_xfer_patch[p]     = ns.livestemn_xfer_patch[p]     + nf.livestemn_storage_to_xfer_patch[p] * dt
                ns.deadstemn_storage_patch[p]  = ns.deadstemn_storage_patch[p]  - nf.deadstemn_storage_to_xfer_patch[p] * dt
                ns.deadstemn_xfer_patch[p]     = ns.deadstemn_xfer_patch[p]     + nf.deadstemn_storage_to_xfer_patch[p] * dt
                ns.livecrootn_storage_patch[p] = ns.livecrootn_storage_patch[p] - nf.livecrootn_storage_to_xfer_patch[p] * dt
                ns.livecrootn_xfer_patch[p]    = ns.livecrootn_xfer_patch[p]    + nf.livecrootn_storage_to_xfer_patch[p] * dt
                ns.deadcrootn_storage_patch[p] = ns.deadcrootn_storage_patch[p] - nf.deadcrootn_storage_to_xfer_patch[p] * dt
                ns.deadcrootn_xfer_patch[p]    = ns.deadcrootn_xfer_patch[p]    + nf.deadcrootn_storage_to_xfer_patch[p] * dt
            end
        end  # !use_matrixcn

        if ivt[p] >= npcropmin  # skip 2 generic crops
            # lines here for consistency; the transfer terms are zero
            if !use_matrixcn
                ns.livestemn_storage_patch[p] = ns.livestemn_storage_patch[p] - nf.livestemn_storage_to_xfer_patch[p] * dt
                ns.livestemn_xfer_patch[p]    = ns.livestemn_xfer_patch[p]    + nf.livestemn_storage_to_xfer_patch[p] * dt
                for k in 1:nrepr
                    ns.reproductiven_storage_patch[p, k] = ns.reproductiven_storage_patch[p, k] -
                        nf.reproductiven_storage_to_xfer_patch[p, k] * dt
                    ns.reproductiven_xfer_patch[p, k] = ns.reproductiven_xfer_patch[p, k] +
                        nf.reproductiven_storage_to_xfer_patch[p, k] * dt
                end
            end  # !use_matrixcn
        end
    end
end

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

Ported from `NStateUpdate1` in `CNNStateUpdate1Mod.F90`. Runs as two
KernelAbstractions kernels (one per soil column, one per patch); the CPU
path is byte-identical.

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
                           mask_soilc::AbstractVector{Bool},
                           mask_soilp::AbstractVector{Bool},
                           bounds_col::UnitRange{Int},
                           bounds_patch::UnitRange{Int},
                           ivt::AbstractVector{<:Integer},
                           woody::AbstractVector,
                           col_is_fates::AbstractVector{Bool},
                           nlevdecomp::Int,
                           i_litr_min::Int,
                           i_litr_max::Int,
                           i_cwd::Int,
                           cascade_donor_pool::AbstractVector{<:Integer} = Int[],
                           cascade_receiver_pool::AbstractVector{<:Integer} = Int[],
                           ndecomp_cascade_transitions::Int = 0,
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

    # --- Column loop: soil decomposition input fluxes (one thread per column) ---
    _FT = eltype(nf_soil.decomp_npools_sourcesink_col)
    _donor    = _to_backend_like(nf_soil.decomp_npools_sourcesink_col, _FT, cascade_donor_pool)
    _receiver = _to_backend_like(nf_soil.decomp_npools_sourcesink_col, _FT, cascade_receiver_pool)
    _launch!(_nsu1_col_kernel!, mask_soilc, col_is_fates,
        nf_veg.phenology_n_to_litr_n_col, nf_soil.decomp_npools_sourcesink_col,
        nf_soil.decomp_cascade_ntransfer_vr_col, nf_soil.decomp_cascade_sminn_flux_vr_col,
        _donor, _receiver,
        nlevdecomp, ndecomp_cascade_transitions, i_litr_min, i_litr_max, i_cwd,
        use_soil_matrixcn, _FT(dt))

    # --- Patch loop: vegetation N state updates (one thread per patch) ---
    ns = _nsu1_ns(ns_veg)
    nf = _nsu1_nf(nf_veg)
    _ivtp  = _to_backend_like(nf_soil.decomp_npools_sourcesink_col, _FT, ivt)
    _woodp = _to_backend_like(nf_soil.decomp_npools_sourcesink_col, _FT, woody)
    _launch!(_nsu1_patch_kernel!, mask_soilp, _ivtp, _woodp, ns, nf,
        npcropmin, nrepr, repr_grain_min, repr_grain_max,
        repr_structure_min, repr_structure_max,
        use_matrixcn, use_fun, _FT(dt))

    return nothing
end
