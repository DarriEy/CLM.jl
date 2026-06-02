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
#
# GPU kernelization (Phase B C/N state-update cascade):
#   * Each function's column loop -> a one-thread-per-column kernel; the internal
#     j decomposition-level loop and i litter-pool loop run sequentially in-thread
#     (each thread owns its column, so the RMW into decomp_npools_vr_col is
#     byte-identical to the host loop and race-free across threads).
#   * Each function's patch loop  -> a one-thread-per-patch kernel; every write is
#     to the patch's own index (no patch->column scatter). The ~20 N-state arrays +
#     ~19/21 N-flux arrays are grouped into immutable `@adapt_structure` device-view
#     bundles (`_NSU2NS` for the shared veg N state, `_NSU2NF`/`_NSU2HNF`/`_NSU2GNF`
#     for each function's veg N fluxes) so the launch stays under Metal's ~31-arg
#     limit; bundle field names mirror the state structs so the body reads verbatim.
#   Float64 literals are eltype-converted (none appear in the body besides `dt`,
#   which carries the working precision); on Float64 this is byte-identical.
# ==========================================================================

# --- Device-view bundle: CNVeg nitrogen STATE arrays the patch loops read/write.
# Shared across all three functions (identical 20 displayed/storage/transfer pools).
# `V` = the 1D patch vectors. ---
Base.@kwdef struct _NSU2NS{V}
    leafn_patch::V; frootn_patch::V; livestemn_patch::V; deadstemn_patch::V
    livecrootn_patch::V; deadcrootn_patch::V; retransn_patch::V
    leafn_storage_patch::V; frootn_storage_patch::V; livestemn_storage_patch::V
    deadstemn_storage_patch::V; livecrootn_storage_patch::V; deadcrootn_storage_patch::V
    leafn_xfer_patch::V; frootn_xfer_patch::V; livestemn_xfer_patch::V
    deadstemn_xfer_patch::V; livecrootn_xfer_patch::V; deadcrootn_xfer_patch::V
end
Adapt.@adapt_structure _NSU2NS

_nsu2_ns(ns) = _NSU2NS(;
    leafn_patch = ns.leafn_patch, frootn_patch = ns.frootn_patch,
    livestemn_patch = ns.livestemn_patch, deadstemn_patch = ns.deadstemn_patch,
    livecrootn_patch = ns.livecrootn_patch, deadcrootn_patch = ns.deadcrootn_patch,
    retransn_patch = ns.retransn_patch,
    leafn_storage_patch = ns.leafn_storage_patch,
    frootn_storage_patch = ns.frootn_storage_patch,
    livestemn_storage_patch = ns.livestemn_storage_patch,
    deadstemn_storage_patch = ns.deadstemn_storage_patch,
    livecrootn_storage_patch = ns.livecrootn_storage_patch,
    deadcrootn_storage_patch = ns.deadcrootn_storage_patch,
    leafn_xfer_patch = ns.leafn_xfer_patch,
    frootn_xfer_patch = ns.frootn_xfer_patch,
    livestemn_xfer_patch = ns.livestemn_xfer_patch,
    deadstemn_xfer_patch = ns.deadstemn_xfer_patch,
    livecrootn_xfer_patch = ns.livecrootn_xfer_patch,
    deadcrootn_xfer_patch = ns.deadcrootn_xfer_patch)

# ---------------------------------------------------------------------------
# n_state_update2! — Gap-phase mortality
# ---------------------------------------------------------------------------

# --- Device-view bundle: CNVeg nitrogen FLUX arrays the gap-mortality patch loop reads ---
Base.@kwdef struct _NSU2NF{V}
    m_leafn_to_litter_patch::V; m_frootn_to_litter_patch::V
    m_livestemn_to_litter_patch::V; m_deadstemn_to_litter_patch::V
    m_livecrootn_to_litter_patch::V; m_deadcrootn_to_litter_patch::V
    m_retransn_to_litter_patch::V
    m_leafn_storage_to_litter_patch::V; m_frootn_storage_to_litter_patch::V
    m_livestemn_storage_to_litter_patch::V; m_deadstemn_storage_to_litter_patch::V
    m_livecrootn_storage_to_litter_patch::V; m_deadcrootn_storage_to_litter_patch::V
    m_leafn_xfer_to_litter_patch::V; m_frootn_xfer_to_litter_patch::V
    m_livestemn_xfer_to_litter_patch::V; m_deadstemn_xfer_to_litter_patch::V
    m_livecrootn_xfer_to_litter_patch::V; m_deadcrootn_xfer_to_litter_patch::V
end
Adapt.@adapt_structure _NSU2NF

_nsu2_nf(nf) = _NSU2NF(;
    m_leafn_to_litter_patch = nf.m_leafn_to_litter_patch,
    m_frootn_to_litter_patch = nf.m_frootn_to_litter_patch,
    m_livestemn_to_litter_patch = nf.m_livestemn_to_litter_patch,
    m_deadstemn_to_litter_patch = nf.m_deadstemn_to_litter_patch,
    m_livecrootn_to_litter_patch = nf.m_livecrootn_to_litter_patch,
    m_deadcrootn_to_litter_patch = nf.m_deadcrootn_to_litter_patch,
    m_retransn_to_litter_patch = nf.m_retransn_to_litter_patch,
    m_leafn_storage_to_litter_patch = nf.m_leafn_storage_to_litter_patch,
    m_frootn_storage_to_litter_patch = nf.m_frootn_storage_to_litter_patch,
    m_livestemn_storage_to_litter_patch = nf.m_livestemn_storage_to_litter_patch,
    m_deadstemn_storage_to_litter_patch = nf.m_deadstemn_storage_to_litter_patch,
    m_livecrootn_storage_to_litter_patch = nf.m_livecrootn_storage_to_litter_patch,
    m_deadcrootn_storage_to_litter_patch = nf.m_deadcrootn_storage_to_litter_patch,
    m_leafn_xfer_to_litter_patch = nf.m_leafn_xfer_to_litter_patch,
    m_frootn_xfer_to_litter_patch = nf.m_frootn_xfer_to_litter_patch,
    m_livestemn_xfer_to_litter_patch = nf.m_livestemn_xfer_to_litter_patch,
    m_deadstemn_xfer_to_litter_patch = nf.m_deadstemn_xfer_to_litter_patch,
    m_livecrootn_xfer_to_litter_patch = nf.m_livecrootn_xfer_to_litter_patch,
    m_deadcrootn_xfer_to_litter_patch = nf.m_deadcrootn_xfer_to_litter_patch)

# --- Kernel: column-level gap-mortality fluxes to soil pools (one thread per column) ---
@kernel function _nsu2_col_kernel!(@Const(mask_soilc), decomp_npools_vr_col,
        @Const(gap_mortality_n_to_litr_n_col), @Const(gap_mortality_n_to_cwdn_col),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        use_soil_matrixcn::Bool, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        if !use_soil_matrixcn
            for j in 1:nlevdecomp
                for i in i_litr_min:i_litr_max
                    decomp_npools_vr_col[c, j, i] =
                        decomp_npools_vr_col[c, j, i] +
                        gap_mortality_n_to_litr_n_col[c, j, i] * dt
                end
                # Currently i_cwd != i_litr_max + 1 if not fates and
                # i_cwd = 0 if fates, so not including in the i-loop
                decomp_npools_vr_col[c, j, i_cwd] =
                    decomp_npools_vr_col[c, j, i_cwd] +
                    gap_mortality_n_to_cwdn_col[c, j] * dt
            end
        else
            # Matrix path: matrix_Ninput setup omitted (requires full matrix infrastructure)
        end
    end
end

# --- Kernel: patch-level veg N state updates from gap mortality (one thread per patch) ---
@kernel function _nsu2_patch_kernel!(@Const(mask_soilp), ns, nf,
        use_matrixcn::Bool, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        if !use_matrixcn
            # displayed pools
            ns.leafn_patch[p] = ns.leafn_patch[p] -
                nf.m_leafn_to_litter_patch[p] * dt
            ns.frootn_patch[p] = ns.frootn_patch[p] -
                nf.m_frootn_to_litter_patch[p] * dt
            ns.livestemn_patch[p] = ns.livestemn_patch[p] -
                nf.m_livestemn_to_litter_patch[p] * dt
            ns.deadstemn_patch[p] = ns.deadstemn_patch[p] -
                nf.m_deadstemn_to_litter_patch[p] * dt
            ns.livecrootn_patch[p] = ns.livecrootn_patch[p] -
                nf.m_livecrootn_to_litter_patch[p] * dt
            ns.deadcrootn_patch[p] = ns.deadcrootn_patch[p] -
                nf.m_deadcrootn_to_litter_patch[p] * dt
            ns.retransn_patch[p] = ns.retransn_patch[p] -
                nf.m_retransn_to_litter_patch[p] * dt

            # storage pools
            ns.leafn_storage_patch[p] = ns.leafn_storage_patch[p] -
                nf.m_leafn_storage_to_litter_patch[p] * dt
            ns.frootn_storage_patch[p] = ns.frootn_storage_patch[p] -
                nf.m_frootn_storage_to_litter_patch[p] * dt
            ns.livestemn_storage_patch[p] = ns.livestemn_storage_patch[p] -
                nf.m_livestemn_storage_to_litter_patch[p] * dt
            ns.deadstemn_storage_patch[p] = ns.deadstemn_storage_patch[p] -
                nf.m_deadstemn_storage_to_litter_patch[p] * dt
            ns.livecrootn_storage_patch[p] = ns.livecrootn_storage_patch[p] -
                nf.m_livecrootn_storage_to_litter_patch[p] * dt
            ns.deadcrootn_storage_patch[p] = ns.deadcrootn_storage_patch[p] -
                nf.m_deadcrootn_storage_to_litter_patch[p] * dt

            # transfer pools
            ns.leafn_xfer_patch[p] = ns.leafn_xfer_patch[p] -
                nf.m_leafn_xfer_to_litter_patch[p] * dt
            ns.frootn_xfer_patch[p] = ns.frootn_xfer_patch[p] -
                nf.m_frootn_xfer_to_litter_patch[p] * dt
            ns.livestemn_xfer_patch[p] = ns.livestemn_xfer_patch[p] -
                nf.m_livestemn_xfer_to_litter_patch[p] * dt
            ns.deadstemn_xfer_patch[p] = ns.deadstemn_xfer_patch[p] -
                nf.m_deadstemn_xfer_to_litter_patch[p] * dt
            ns.livecrootn_xfer_patch[p] = ns.livecrootn_xfer_patch[p] -
                nf.m_livecrootn_xfer_to_litter_patch[p] * dt
            ns.deadcrootn_xfer_patch[p] = ns.deadcrootn_xfer_patch[p] -
                nf.m_deadcrootn_xfer_to_litter_patch[p] * dt
        else
            # Matrix version: handled in CNGapMortality
        end
    end
end

"""
    n_state_update2!(ns_veg, nf_veg, ns_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

On the radiation time step, update all prognostic nitrogen state variables
affected by gap-phase mortality fluxes.

Ported from `NStateUpdate2` in `CNNStateUpdate2Mod.F90`. Runs as two
KernelAbstractions kernels (one per soil column, one per patch); the CPU path
is byte-identical.

Note: The matrix-CN code paths (`use_soil_matrixcn`) for soil matrix
solutions are omitted since they require the full matrix infrastructure
(matrix_Ninput). The matrix-CN veg path (`use_matrixcn`) is a no-op in
the Fortran original (comment only), so it is omitted here as well.
"""
function n_state_update2!(ns_veg::CNVegNitrogenStateData,
                           nf_veg::CNVegNitrogenFluxData,
                           ns_soil::SoilBiogeochemNitrogenStateData;
                           mask_soilc::AbstractVector{Bool},
                           mask_soilp::AbstractVector{Bool},
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
    _launch!(_nsu2_col_kernel!, mask_soilc, ns_soil.decomp_npools_vr_col,
        nf_veg.gap_mortality_n_to_litr_n_col, nf_veg.gap_mortality_n_to_cwdn_col,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd, use_soil_matrixcn, dt)

    # --- Patch loop: vegetation N state updates from gap-phase mortality ---
    ns = _nsu2_ns(ns_veg)
    nf = _nsu2_nf(nf_veg)
    _launch!(_nsu2_patch_kernel!, mask_soilp, ns, nf, use_matrixcn, dt)

    return nothing
end

# ---------------------------------------------------------------------------
# n_state_update2h! — Harvest mortality
# ---------------------------------------------------------------------------

# --- Device-view bundle: CNVeg nitrogen FLUX arrays the harvest patch loop reads ---
Base.@kwdef struct _NSU2HNF{V}
    hrv_leafn_to_litter_patch::V; hrv_frootn_to_litter_patch::V
    hrv_livestemn_to_litter_patch::V; wood_harvestn_patch::V
    hrv_livecrootn_to_litter_patch::V; hrv_deadcrootn_to_litter_patch::V
    hrv_retransn_to_litter_patch::V
    hrv_leafn_storage_to_litter_patch::V; hrv_frootn_storage_to_litter_patch::V
    hrv_livestemn_storage_to_litter_patch::V; hrv_deadstemn_storage_to_litter_patch::V
    hrv_livecrootn_storage_to_litter_patch::V; hrv_deadcrootn_storage_to_litter_patch::V
    hrv_leafn_xfer_to_litter_patch::V; hrv_frootn_xfer_to_litter_patch::V
    hrv_livestemn_xfer_to_litter_patch::V; hrv_deadstemn_xfer_to_litter_patch::V
    hrv_livecrootn_xfer_to_litter_patch::V; hrv_deadcrootn_xfer_to_litter_patch::V
end
Adapt.@adapt_structure _NSU2HNF

_nsu2h_nf(nf) = _NSU2HNF(;
    hrv_leafn_to_litter_patch = nf.hrv_leafn_to_litter_patch,
    hrv_frootn_to_litter_patch = nf.hrv_frootn_to_litter_patch,
    hrv_livestemn_to_litter_patch = nf.hrv_livestemn_to_litter_patch,
    wood_harvestn_patch = nf.wood_harvestn_patch,
    hrv_livecrootn_to_litter_patch = nf.hrv_livecrootn_to_litter_patch,
    hrv_deadcrootn_to_litter_patch = nf.hrv_deadcrootn_to_litter_patch,
    hrv_retransn_to_litter_patch = nf.hrv_retransn_to_litter_patch,
    hrv_leafn_storage_to_litter_patch = nf.hrv_leafn_storage_to_litter_patch,
    hrv_frootn_storage_to_litter_patch = nf.hrv_frootn_storage_to_litter_patch,
    hrv_livestemn_storage_to_litter_patch = nf.hrv_livestemn_storage_to_litter_patch,
    hrv_deadstemn_storage_to_litter_patch = nf.hrv_deadstemn_storage_to_litter_patch,
    hrv_livecrootn_storage_to_litter_patch = nf.hrv_livecrootn_storage_to_litter_patch,
    hrv_deadcrootn_storage_to_litter_patch = nf.hrv_deadcrootn_storage_to_litter_patch,
    hrv_leafn_xfer_to_litter_patch = nf.hrv_leafn_xfer_to_litter_patch,
    hrv_frootn_xfer_to_litter_patch = nf.hrv_frootn_xfer_to_litter_patch,
    hrv_livestemn_xfer_to_litter_patch = nf.hrv_livestemn_xfer_to_litter_patch,
    hrv_deadstemn_xfer_to_litter_patch = nf.hrv_deadstemn_xfer_to_litter_patch,
    hrv_livecrootn_xfer_to_litter_patch = nf.hrv_livecrootn_xfer_to_litter_patch,
    hrv_deadcrootn_xfer_to_litter_patch = nf.hrv_deadcrootn_xfer_to_litter_patch)

# --- Kernel: column-level harvest fluxes to soil pools (one thread per column) ---
@kernel function _nsu2h_col_kernel!(@Const(mask_soilc), decomp_npools_vr_col,
        @Const(harvest_n_to_litr_n_col), @Const(harvest_n_to_cwdn_col),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        use_soil_matrixcn::Bool, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        if !use_soil_matrixcn
            for j in 1:nlevdecomp
                for i in i_litr_min:i_litr_max
                    decomp_npools_vr_col[c, j, i] =
                        decomp_npools_vr_col[c, j, i] +
                        harvest_n_to_litr_n_col[c, j, i] * dt
                end
                # Currently i_cwd != i_litr_max + 1 if not fates and
                # i_cwd = 0 if fates, so not including in the i-loop
                decomp_npools_vr_col[c, j, i_cwd] =
                    decomp_npools_vr_col[c, j, i_cwd] +
                    harvest_n_to_cwdn_col[c, j] * dt
            end
        else
            # Matrix path: matrix_Ninput setup omitted (requires full matrix infrastructure)
        end
    end
end

# --- Kernel: patch-level veg N state updates from harvest mortality (one thread per patch) ---
@kernel function _nsu2h_patch_kernel!(@Const(mask_soilp), ns, nf,
        use_matrixcn::Bool, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        if !use_matrixcn
            # displayed pools
            ns.leafn_patch[p] = ns.leafn_patch[p] -
                nf.hrv_leafn_to_litter_patch[p] * dt
            ns.frootn_patch[p] = ns.frootn_patch[p] -
                nf.hrv_frootn_to_litter_patch[p] * dt
            ns.livestemn_patch[p] = ns.livestemn_patch[p] -
                nf.hrv_livestemn_to_litter_patch[p] * dt
            ns.deadstemn_patch[p] = ns.deadstemn_patch[p] -
                nf.wood_harvestn_patch[p] * dt
            ns.livecrootn_patch[p] = ns.livecrootn_patch[p] -
                nf.hrv_livecrootn_to_litter_patch[p] * dt
            ns.deadcrootn_patch[p] = ns.deadcrootn_patch[p] -
                nf.hrv_deadcrootn_to_litter_patch[p] * dt
            ns.retransn_patch[p] = ns.retransn_patch[p] -
                nf.hrv_retransn_to_litter_patch[p] * dt

            # storage pools
            ns.leafn_storage_patch[p] = ns.leafn_storage_patch[p] -
                nf.hrv_leafn_storage_to_litter_patch[p] * dt
            ns.frootn_storage_patch[p] = ns.frootn_storage_patch[p] -
                nf.hrv_frootn_storage_to_litter_patch[p] * dt
            ns.livestemn_storage_patch[p] = ns.livestemn_storage_patch[p] -
                nf.hrv_livestemn_storage_to_litter_patch[p] * dt
            ns.deadstemn_storage_patch[p] = ns.deadstemn_storage_patch[p] -
                nf.hrv_deadstemn_storage_to_litter_patch[p] * dt
            ns.livecrootn_storage_patch[p] = ns.livecrootn_storage_patch[p] -
                nf.hrv_livecrootn_storage_to_litter_patch[p] * dt
            ns.deadcrootn_storage_patch[p] = ns.deadcrootn_storage_patch[p] -
                nf.hrv_deadcrootn_storage_to_litter_patch[p] * dt

            # transfer pools
            ns.leafn_xfer_patch[p] = ns.leafn_xfer_patch[p] -
                nf.hrv_leafn_xfer_to_litter_patch[p] * dt
            ns.frootn_xfer_patch[p] = ns.frootn_xfer_patch[p] -
                nf.hrv_frootn_xfer_to_litter_patch[p] * dt
            ns.livestemn_xfer_patch[p] = ns.livestemn_xfer_patch[p] -
                nf.hrv_livestemn_xfer_to_litter_patch[p] * dt
            ns.deadstemn_xfer_patch[p] = ns.deadstemn_xfer_patch[p] -
                nf.hrv_deadstemn_xfer_to_litter_patch[p] * dt
            ns.livecrootn_xfer_patch[p] = ns.livecrootn_xfer_patch[p] -
                nf.hrv_livecrootn_xfer_to_litter_patch[p] * dt
            ns.deadcrootn_xfer_patch[p] = ns.deadcrootn_xfer_patch[p] -
                nf.hrv_deadcrootn_xfer_to_litter_patch[p] * dt
        else
            # Matrix version: handled in CNHarvest
        end
    end
end

"""
    n_state_update2h!(ns_veg, nf_veg, ns_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

Update all prognostic nitrogen state variables affected by harvest
mortality fluxes.

Ported from `NStateUpdate2h` in `CNNStateUpdate2Mod.F90`. Runs as two
KernelAbstractions kernels (one per soil column, one per patch); the CPU path
is byte-identical.

Note: The matrix-CN code paths are omitted (require full matrix infrastructure).
"""
function n_state_update2h!(ns_veg::CNVegNitrogenStateData,
                            nf_veg::CNVegNitrogenFluxData,
                            ns_soil::SoilBiogeochemNitrogenStateData;
                            mask_soilc::AbstractVector{Bool},
                            mask_soilp::AbstractVector{Bool},
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
    _launch!(_nsu2h_col_kernel!, mask_soilc, ns_soil.decomp_npools_vr_col,
        nf_veg.harvest_n_to_litr_n_col, nf_veg.harvest_n_to_cwdn_col,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd, use_soil_matrixcn, dt)

    # --- Patch loop: vegetation N state updates from harvest mortality ---
    ns = _nsu2_ns(ns_veg)
    nf = _nsu2h_nf(nf_veg)
    _launch!(_nsu2h_patch_kernel!, mask_soilp, ns, nf, use_matrixcn, dt)

    return nothing
end

# ---------------------------------------------------------------------------
# n_state_update2g! — Gross unrepresented landcover change mortality
# ---------------------------------------------------------------------------

# --- Device-view bundle: CNVeg nitrogen FLUX arrays the gross-unrep patch loop reads ---
Base.@kwdef struct _NSU2GNF{V}
    gru_leafn_to_litter_patch::V; gru_frootn_to_litter_patch::V
    gru_livestemn_to_atm_patch::V; gru_deadstemn_to_atm_patch::V
    gru_wood_productn_gain_patch::V
    gru_livecrootn_to_litter_patch::V; gru_deadcrootn_to_litter_patch::V
    gru_retransn_to_litter_patch::V
    gru_leafn_storage_to_atm_patch::V; gru_frootn_storage_to_atm_patch::V
    gru_livestemn_storage_to_atm_patch::V; gru_deadstemn_storage_to_atm_patch::V
    gru_livecrootn_storage_to_atm_patch::V; gru_deadcrootn_storage_to_atm_patch::V
    gru_leafn_xfer_to_atm_patch::V; gru_frootn_xfer_to_atm_patch::V
    gru_livestemn_xfer_to_atm_patch::V; gru_deadstemn_xfer_to_atm_patch::V
    gru_livecrootn_xfer_to_atm_patch::V; gru_deadcrootn_xfer_to_atm_patch::V
end
Adapt.@adapt_structure _NSU2GNF

_nsu2g_nf(nf) = _NSU2GNF(;
    gru_leafn_to_litter_patch = nf.gru_leafn_to_litter_patch,
    gru_frootn_to_litter_patch = nf.gru_frootn_to_litter_patch,
    gru_livestemn_to_atm_patch = nf.gru_livestemn_to_atm_patch,
    gru_deadstemn_to_atm_patch = nf.gru_deadstemn_to_atm_patch,
    gru_wood_productn_gain_patch = nf.gru_wood_productn_gain_patch,
    gru_livecrootn_to_litter_patch = nf.gru_livecrootn_to_litter_patch,
    gru_deadcrootn_to_litter_patch = nf.gru_deadcrootn_to_litter_patch,
    gru_retransn_to_litter_patch = nf.gru_retransn_to_litter_patch,
    gru_leafn_storage_to_atm_patch = nf.gru_leafn_storage_to_atm_patch,
    gru_frootn_storage_to_atm_patch = nf.gru_frootn_storage_to_atm_patch,
    gru_livestemn_storage_to_atm_patch = nf.gru_livestemn_storage_to_atm_patch,
    gru_deadstemn_storage_to_atm_patch = nf.gru_deadstemn_storage_to_atm_patch,
    gru_livecrootn_storage_to_atm_patch = nf.gru_livecrootn_storage_to_atm_patch,
    gru_deadcrootn_storage_to_atm_patch = nf.gru_deadcrootn_storage_to_atm_patch,
    gru_leafn_xfer_to_atm_patch = nf.gru_leafn_xfer_to_atm_patch,
    gru_frootn_xfer_to_atm_patch = nf.gru_frootn_xfer_to_atm_patch,
    gru_livestemn_xfer_to_atm_patch = nf.gru_livestemn_xfer_to_atm_patch,
    gru_deadstemn_xfer_to_atm_patch = nf.gru_deadstemn_xfer_to_atm_patch,
    gru_livecrootn_xfer_to_atm_patch = nf.gru_livecrootn_xfer_to_atm_patch,
    gru_deadcrootn_xfer_to_atm_patch = nf.gru_deadcrootn_xfer_to_atm_patch)

# --- Kernel: column-level gross-unrep fluxes to soil pools (one thread per column) ---
@kernel function _nsu2g_col_kernel!(@Const(mask_soilc), decomp_npools_vr_col,
        @Const(gru_n_to_litr_n_col), @Const(gru_n_to_cwdn_col),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        use_soil_matrixcn::Bool, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        if !use_soil_matrixcn
            for j in 1:nlevdecomp
                for i in i_litr_min:i_litr_max
                    decomp_npools_vr_col[c, j, i] =
                        decomp_npools_vr_col[c, j, i] +
                        gru_n_to_litr_n_col[c, j, i] * dt
                end
                # Currently i_cwd != i_litr_max + 1 if not fates and
                # i_cwd = 0 if fates, so not including in the i-loop
                decomp_npools_vr_col[c, j, i_cwd] =
                    decomp_npools_vr_col[c, j, i_cwd] +
                    gru_n_to_cwdn_col[c, j] * dt
            end
        else
            # Matrix path: matrix_Ninput setup omitted (requires full matrix infrastructure)
        end
    end
end

# --- Kernel: patch-level veg N state updates from gross-unrep LCC (one thread per patch) ---
@kernel function _nsu2g_patch_kernel!(@Const(mask_soilp), ns, nf,
        use_matrixcn::Bool, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        if !use_matrixcn
            # displayed pools
            ns.leafn_patch[p] = ns.leafn_patch[p] -
                nf.gru_leafn_to_litter_patch[p] * dt
            ns.frootn_patch[p] = ns.frootn_patch[p] -
                nf.gru_frootn_to_litter_patch[p] * dt
            ns.livestemn_patch[p] = ns.livestemn_patch[p] -
                nf.gru_livestemn_to_atm_patch[p] * dt
            ns.deadstemn_patch[p] = ns.deadstemn_patch[p] -
                nf.gru_deadstemn_to_atm_patch[p] * dt
            ns.deadstemn_patch[p] = ns.deadstemn_patch[p] -
                nf.gru_wood_productn_gain_patch[p] * dt
            ns.livecrootn_patch[p] = ns.livecrootn_patch[p] -
                nf.gru_livecrootn_to_litter_patch[p] * dt
            ns.deadcrootn_patch[p] = ns.deadcrootn_patch[p] -
                nf.gru_deadcrootn_to_litter_patch[p] * dt
            ns.retransn_patch[p] = ns.retransn_patch[p] -
                nf.gru_retransn_to_litter_patch[p] * dt

            # storage pools
            ns.leafn_storage_patch[p] = ns.leafn_storage_patch[p] -
                nf.gru_leafn_storage_to_atm_patch[p] * dt
            ns.frootn_storage_patch[p] = ns.frootn_storage_patch[p] -
                nf.gru_frootn_storage_to_atm_patch[p] * dt
            ns.livestemn_storage_patch[p] = ns.livestemn_storage_patch[p] -
                nf.gru_livestemn_storage_to_atm_patch[p] * dt
            ns.deadstemn_storage_patch[p] = ns.deadstemn_storage_patch[p] -
                nf.gru_deadstemn_storage_to_atm_patch[p] * dt
            ns.livecrootn_storage_patch[p] = ns.livecrootn_storage_patch[p] -
                nf.gru_livecrootn_storage_to_atm_patch[p] * dt
            ns.deadcrootn_storage_patch[p] = ns.deadcrootn_storage_patch[p] -
                nf.gru_deadcrootn_storage_to_atm_patch[p] * dt

            # transfer pools
            ns.leafn_xfer_patch[p] = ns.leafn_xfer_patch[p] -
                nf.gru_leafn_xfer_to_atm_patch[p] * dt
            ns.frootn_xfer_patch[p] = ns.frootn_xfer_patch[p] -
                nf.gru_frootn_xfer_to_atm_patch[p] * dt
            ns.livestemn_xfer_patch[p] = ns.livestemn_xfer_patch[p] -
                nf.gru_livestemn_xfer_to_atm_patch[p] * dt
            ns.deadstemn_xfer_patch[p] = ns.deadstemn_xfer_patch[p] -
                nf.gru_deadstemn_xfer_to_atm_patch[p] * dt
            ns.livecrootn_xfer_patch[p] = ns.livecrootn_xfer_patch[p] -
                nf.gru_livecrootn_xfer_to_atm_patch[p] * dt
            ns.deadcrootn_xfer_patch[p] = ns.deadcrootn_xfer_patch[p] -
                nf.gru_deadcrootn_xfer_to_atm_patch[p] * dt
        else
            # Matrix version: handled in dynGrossUnrepMod::CNGrossUnrep
        end
    end
end

"""
    n_state_update2g!(ns_veg, nf_veg, ns_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

Update all prognostic nitrogen state variables affected by gross
unrepresented landcover change mortality fluxes.

Ported from `NStateUpdate2g` in `CNNStateUpdate2Mod.F90`. Runs as two
KernelAbstractions kernels (one per soil column, one per patch); the CPU path
is byte-identical.

Note: The matrix-CN code paths are omitted (require full matrix infrastructure).
"""
function n_state_update2g!(ns_veg::CNVegNitrogenStateData,
                            nf_veg::CNVegNitrogenFluxData,
                            ns_soil::SoilBiogeochemNitrogenStateData;
                            mask_soilc::AbstractVector{Bool},
                            mask_soilp::AbstractVector{Bool},
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
    _launch!(_nsu2g_col_kernel!, mask_soilc, ns_soil.decomp_npools_vr_col,
        nf_veg.gru_n_to_litr_n_col, nf_veg.gru_n_to_cwdn_col,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd, use_soil_matrixcn, dt)

    # --- Patch loop: vegetation N state updates from gross unrepresented landcover change ---
    ns = _nsu2_ns(ns_veg)
    nf = _nsu2g_nf(nf_veg)
    _launch!(_nsu2g_patch_kernel!, mask_soilp, ns, nf, use_matrixcn, dt)

    return nothing
end
