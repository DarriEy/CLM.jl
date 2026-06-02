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
#
# GPU kernelization (Phase B C/N state-update cascade):
#   * n_state_update_leaching! — small per-column j loop -> `_nsu_leach_kernel!`
#     (one thread per soil column; the internal nlevdecomp j loop runs sequentially
#     in-thread, each thread owns its column). ~4 soil-N arrays passed loose.
#   * n_state_update3! — three kernels: two column kernels for the soil-pool fire
#     fluxes (`_nsu3_col_cwdlitr_kernel!` for the CWD+litter inputs, and
#     `_nsu3_col_firelosses_kernel!` for the per-pool fire losses; each is one
#     thread per soil column with internal sequential j/k/l loops), and one patch
#     kernel `_nsu3_patch_kernel!` (one thread per patch; every write is to the
#     patch's own index, no patch->column scatter). The ~38 vegetation N
#     state/flux arrays the patch loop touches are grouped into two immutable
#     `@adapt_structure` device-view bundles (`_NSU3NS{V}` / `_NSU3NF{V}`) so the
#     launch stays under Metal's ~31-arg limit; the bundle field names mirror the
#     state structs so the body reads verbatim.
# Float64 literals are eltype-converted (`zero(T)`) so the kernels carry no Float64
# on a Float32-only backend; on Float64 this is byte-identical.
# ==========================================================================

# ---------------------------------------------------------------------------
# n_state_update_leaching! — Sminn leaching
# ---------------------------------------------------------------------------

# --- Kernel: per-column sminn leaching (one thread per soil column) ---
@kernel function _nsu_leach_kernel!(@Const(mask_soilc),
        sminn_vr_col, smin_no3_vr_col, @Const(smin_nh4_vr_col),
        @Const(sminn_leached_vr_col),
        @Const(smin_no3_leached_vr_col), @Const(smin_no3_runoff_vr_col),
        nlevdecomp::Int, use_nitrif_denitrif::Bool, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        T = typeof(dt)
        for j in 1:nlevdecomp
            if !use_nitrif_denitrif
                # mineral N loss due to leaching
                sminn_vr_col[c, j] = sminn_vr_col[c, j] -
                    sminn_leached_vr_col[c, j] * dt
            else
                # mineral N loss due to leaching and runoff
                smin_no3_vr_col[c, j] = max(
                    smin_no3_vr_col[c, j] -
                    (smin_no3_leached_vr_col[c, j] + smin_no3_runoff_vr_col[c, j]) * dt,
                    zero(T))

                sminn_vr_col[c, j] = smin_no3_vr_col[c, j] + smin_nh4_vr_col[c, j]
            end
        end
    end
end

"""
    n_state_update_leaching!(ns_soil, nf_soil;
        mask_soilc, bounds_col, nlevdecomp,
        use_nitrif_denitrif, dt)

On the radiation time step, update all prognostic nitrogen state variables
affected by the Sminn leaching flux.

Ported from `NStateUpdateLeaching` in `CNNStateUpdate3Mod.F90`. Runs as one
KernelAbstractions kernel (one thread per soil column). The CPU path is
byte-identical.

Note: This code was separated from gap mortality fluxes to make it
compatible with FATES.
"""
function n_state_update_leaching!(ns_soil::SoilBiogeochemNitrogenStateData,
                                   nf_soil::SoilBiogeochemNitrogenFluxData;
                                   mask_soilc::AbstractVector{Bool},
                                   bounds_col::UnitRange{Int},
                                   nlevdecomp::Int,
                                   use_nitrif_denitrif::Bool = false,
                                   dt::Real)

    _launch!(_nsu_leach_kernel!, mask_soilc,
        ns_soil.sminn_vr_col, ns_soil.smin_no3_vr_col, ns_soil.smin_nh4_vr_col,
        nf_soil.sminn_leached_vr_col,
        nf_soil.smin_no3_leached_vr_col, nf_soil.smin_no3_runoff_vr_col,
        nlevdecomp, use_nitrif_denitrif, dt)

    return nothing
end

# ---------------------------------------------------------------------------
# n_state_update3! — Fire N state update
# ---------------------------------------------------------------------------

# --- Device-view bundle: CNVeg nitrogen STATE arrays the patch loop reads/writes ---
# `V` = the 1D patch vectors.
Base.@kwdef struct _NSU3NS{V}
    leafn_patch::V; frootn_patch::V; livestemn_patch::V
    deadstemn_patch::V; livecrootn_patch::V; deadcrootn_patch::V
    leafn_storage_patch::V; frootn_storage_patch::V; livestemn_storage_patch::V
    deadstemn_storage_patch::V; livecrootn_storage_patch::V; deadcrootn_storage_patch::V
    leafn_xfer_patch::V; frootn_xfer_patch::V; livestemn_xfer_patch::V
    deadstemn_xfer_patch::V; livecrootn_xfer_patch::V; deadcrootn_xfer_patch::V
    retransn_patch::V
end
Adapt.@adapt_structure _NSU3NS

_nsu3_ns(ns) = _NSU3NS(;
    leafn_patch = ns.leafn_patch, frootn_patch = ns.frootn_patch,
    livestemn_patch = ns.livestemn_patch, deadstemn_patch = ns.deadstemn_patch,
    livecrootn_patch = ns.livecrootn_patch, deadcrootn_patch = ns.deadcrootn_patch,
    leafn_storage_patch = ns.leafn_storage_patch,
    frootn_storage_patch = ns.frootn_storage_patch,
    livestemn_storage_patch = ns.livestemn_storage_patch,
    deadstemn_storage_patch = ns.deadstemn_storage_patch,
    livecrootn_storage_patch = ns.livecrootn_storage_patch,
    deadcrootn_storage_patch = ns.deadcrootn_storage_patch,
    leafn_xfer_patch = ns.leafn_xfer_patch, frootn_xfer_patch = ns.frootn_xfer_patch,
    livestemn_xfer_patch = ns.livestemn_xfer_patch,
    deadstemn_xfer_patch = ns.deadstemn_xfer_patch,
    livecrootn_xfer_patch = ns.livecrootn_xfer_patch,
    deadcrootn_xfer_patch = ns.deadcrootn_xfer_patch,
    retransn_patch = ns.retransn_patch)

# --- Device-view bundle: CNVeg nitrogen FLUX arrays the patch loop reads ---
Base.@kwdef struct _NSU3NF{V}
    # to combustion
    m_leafn_to_fire_patch::V; m_frootn_to_fire_patch::V
    m_livestemn_to_fire_patch::V; m_deadstemn_to_fire_patch::V
    m_livecrootn_to_fire_patch::V; m_deadcrootn_to_fire_patch::V
    m_leafn_storage_to_fire_patch::V; m_frootn_storage_to_fire_patch::V
    m_livestemn_storage_to_fire_patch::V; m_deadstemn_storage_to_fire_patch::V
    m_livecrootn_storage_to_fire_patch::V; m_deadcrootn_storage_to_fire_patch::V
    m_leafn_xfer_to_fire_patch::V; m_frootn_xfer_to_fire_patch::V
    m_livestemn_xfer_to_fire_patch::V; m_deadstemn_xfer_to_fire_patch::V
    m_livecrootn_xfer_to_fire_patch::V; m_deadcrootn_xfer_to_fire_patch::V
    m_retransn_to_fire_patch::V
    # to litter
    m_leafn_to_litter_fire_patch::V; m_frootn_to_litter_fire_patch::V
    m_livestemn_to_litter_fire_patch::V; m_deadstemn_to_litter_fire_patch::V
    m_livecrootn_to_litter_fire_patch::V; m_deadcrootn_to_litter_fire_patch::V
    m_leafn_storage_to_litter_fire_patch::V; m_frootn_storage_to_litter_fire_patch::V
    m_livestemn_storage_to_litter_fire_patch::V; m_deadstemn_storage_to_litter_fire_patch::V
    m_livecrootn_storage_to_litter_fire_patch::V; m_deadcrootn_storage_to_litter_fire_patch::V
    m_leafn_xfer_to_litter_fire_patch::V; m_frootn_xfer_to_litter_fire_patch::V
    m_livestemn_xfer_to_litter_fire_patch::V; m_deadstemn_xfer_to_litter_fire_patch::V
    m_livecrootn_xfer_to_litter_fire_patch::V; m_deadcrootn_xfer_to_litter_fire_patch::V
    m_retransn_to_litter_fire_patch::V
    # live-to-dead fire transfers
    m_livestemn_to_deadstemn_fire_patch::V; m_livecrootn_to_deadcrootn_fire_patch::V
end
Adapt.@adapt_structure _NSU3NF

_nsu3_nf(nf) = _NSU3NF(;
    m_leafn_to_fire_patch = nf.m_leafn_to_fire_patch,
    m_frootn_to_fire_patch = nf.m_frootn_to_fire_patch,
    m_livestemn_to_fire_patch = nf.m_livestemn_to_fire_patch,
    m_deadstemn_to_fire_patch = nf.m_deadstemn_to_fire_patch,
    m_livecrootn_to_fire_patch = nf.m_livecrootn_to_fire_patch,
    m_deadcrootn_to_fire_patch = nf.m_deadcrootn_to_fire_patch,
    m_leafn_storage_to_fire_patch = nf.m_leafn_storage_to_fire_patch,
    m_frootn_storage_to_fire_patch = nf.m_frootn_storage_to_fire_patch,
    m_livestemn_storage_to_fire_patch = nf.m_livestemn_storage_to_fire_patch,
    m_deadstemn_storage_to_fire_patch = nf.m_deadstemn_storage_to_fire_patch,
    m_livecrootn_storage_to_fire_patch = nf.m_livecrootn_storage_to_fire_patch,
    m_deadcrootn_storage_to_fire_patch = nf.m_deadcrootn_storage_to_fire_patch,
    m_leafn_xfer_to_fire_patch = nf.m_leafn_xfer_to_fire_patch,
    m_frootn_xfer_to_fire_patch = nf.m_frootn_xfer_to_fire_patch,
    m_livestemn_xfer_to_fire_patch = nf.m_livestemn_xfer_to_fire_patch,
    m_deadstemn_xfer_to_fire_patch = nf.m_deadstemn_xfer_to_fire_patch,
    m_livecrootn_xfer_to_fire_patch = nf.m_livecrootn_xfer_to_fire_patch,
    m_deadcrootn_xfer_to_fire_patch = nf.m_deadcrootn_xfer_to_fire_patch,
    m_retransn_to_fire_patch = nf.m_retransn_to_fire_patch,
    m_leafn_to_litter_fire_patch = nf.m_leafn_to_litter_fire_patch,
    m_frootn_to_litter_fire_patch = nf.m_frootn_to_litter_fire_patch,
    m_livestemn_to_litter_fire_patch = nf.m_livestemn_to_litter_fire_patch,
    m_deadstemn_to_litter_fire_patch = nf.m_deadstemn_to_litter_fire_patch,
    m_livecrootn_to_litter_fire_patch = nf.m_livecrootn_to_litter_fire_patch,
    m_deadcrootn_to_litter_fire_patch = nf.m_deadcrootn_to_litter_fire_patch,
    m_leafn_storage_to_litter_fire_patch = nf.m_leafn_storage_to_litter_fire_patch,
    m_frootn_storage_to_litter_fire_patch = nf.m_frootn_storage_to_litter_fire_patch,
    m_livestemn_storage_to_litter_fire_patch = nf.m_livestemn_storage_to_litter_fire_patch,
    m_deadstemn_storage_to_litter_fire_patch = nf.m_deadstemn_storage_to_litter_fire_patch,
    m_livecrootn_storage_to_litter_fire_patch = nf.m_livecrootn_storage_to_litter_fire_patch,
    m_deadcrootn_storage_to_litter_fire_patch = nf.m_deadcrootn_storage_to_litter_fire_patch,
    m_leafn_xfer_to_litter_fire_patch = nf.m_leafn_xfer_to_litter_fire_patch,
    m_frootn_xfer_to_litter_fire_patch = nf.m_frootn_xfer_to_litter_fire_patch,
    m_livestemn_xfer_to_litter_fire_patch = nf.m_livestemn_xfer_to_litter_fire_patch,
    m_deadstemn_xfer_to_litter_fire_patch = nf.m_deadstemn_xfer_to_litter_fire_patch,
    m_livecrootn_xfer_to_litter_fire_patch = nf.m_livecrootn_xfer_to_litter_fire_patch,
    m_deadcrootn_xfer_to_litter_fire_patch = nf.m_deadcrootn_xfer_to_litter_fire_patch,
    m_retransn_to_litter_fire_patch = nf.m_retransn_to_litter_fire_patch,
    m_livestemn_to_deadstemn_fire_patch = nf.m_livestemn_to_deadstemn_fire_patch,
    m_livecrootn_to_deadcrootn_fire_patch = nf.m_livecrootn_to_deadcrootn_fire_patch)

# --- Kernel: fire mortality wood -> column CWD + litter (one thread per column) ---
@kernel function _nsu3_col_cwdlitr_kernel!(@Const(mask_soilc),
        decomp_npools_vr_col,
        @Const(fire_mortality_n_to_cwdn_col), @Const(m_n_to_litr_fire_col),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        use_soil_matrixcn::Bool, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        if !use_soil_matrixcn
            for j in 1:nlevdecomp
                # patch-level wood to column-level CWD (uncombusted wood)
                decomp_npools_vr_col[c, j, i_cwd] =
                    decomp_npools_vr_col[c, j, i_cwd] +
                    fire_mortality_n_to_cwdn_col[c, j] * dt

                # patch-level wood to column-level litter (uncombusted wood)
                for k in i_litr_min:i_litr_max
                    decomp_npools_vr_col[c, j, k] =
                        decomp_npools_vr_col[c, j, k] +
                        m_n_to_litr_fire_col[c, j, k] * dt
                end
            end
        else
            # Matrix path: matrix_Ninput setup omitted (requires full matrix infrastructure)
        end
    end
end

# --- Kernel: litter and CWD losses to fire (one thread per column) ---
@kernel function _nsu3_col_firelosses_kernel!(@Const(mask_soilc),
        decomp_npools_vr_col, @Const(m_decomp_npools_to_fire_vr_col),
        nlevdecomp::Int, ndecomp_pools::Int, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        for l in 1:ndecomp_pools
            for j in 1:nlevdecomp
                decomp_npools_vr_col[c, j, l] =
                    decomp_npools_vr_col[c, j, l] -
                    m_decomp_npools_to_fire_vr_col[c, j, l] * dt
            end
        end
    end
end

# --- Kernel: patch-level vegetation N state updates from fire (one thread per patch) ---
# `ns` / `nf` are the _NSU3NS / _NSU3NF device-view bundles; field names mirror the
# state structs so the body below is the host loop verbatim (ns_veg.->ns., nf_veg.->nf.).
@kernel function _nsu3_patch_kernel!(@Const(mask_soilp), ns, nf,
        use_matrixcn::Bool, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        if !use_matrixcn
            # from fire displayed pools
            ns.leafn_patch[p] = ns.leafn_patch[p] -
                nf.m_leafn_to_fire_patch[p] * dt
            ns.frootn_patch[p] = ns.frootn_patch[p] -
                nf.m_frootn_to_fire_patch[p] * dt
            ns.livestemn_patch[p] = ns.livestemn_patch[p] -
                nf.m_livestemn_to_fire_patch[p] * dt
            ns.deadstemn_patch[p] = ns.deadstemn_patch[p] -
                nf.m_deadstemn_to_fire_patch[p] * dt
            ns.livecrootn_patch[p] = ns.livecrootn_patch[p] -
                nf.m_livecrootn_to_fire_patch[p] * dt
            ns.deadcrootn_patch[p] = ns.deadcrootn_patch[p] -
                nf.m_deadcrootn_to_fire_patch[p] * dt

            ns.leafn_patch[p] = ns.leafn_patch[p] -
                nf.m_leafn_to_litter_fire_patch[p] * dt
            ns.frootn_patch[p] = ns.frootn_patch[p] -
                nf.m_frootn_to_litter_fire_patch[p] * dt
            ns.livestemn_patch[p] = ns.livestemn_patch[p] -
                nf.m_livestemn_to_litter_fire_patch[p] * dt -
                nf.m_livestemn_to_deadstemn_fire_patch[p] * dt
            ns.deadstemn_patch[p] = ns.deadstemn_patch[p] -
                nf.m_deadstemn_to_litter_fire_patch[p] * dt +
                nf.m_livestemn_to_deadstemn_fire_patch[p] * dt
            ns.livecrootn_patch[p] = ns.livecrootn_patch[p] -
                nf.m_livecrootn_to_litter_fire_patch[p] * dt -
                nf.m_livecrootn_to_deadcrootn_fire_patch[p] * dt
            ns.deadcrootn_patch[p] = ns.deadcrootn_patch[p] -
                nf.m_deadcrootn_to_litter_fire_patch[p] * dt +
                nf.m_livecrootn_to_deadcrootn_fire_patch[p] * dt

            # storage pools
            ns.leafn_storage_patch[p] = ns.leafn_storage_patch[p] -
                nf.m_leafn_storage_to_fire_patch[p] * dt
            ns.frootn_storage_patch[p] = ns.frootn_storage_patch[p] -
                nf.m_frootn_storage_to_fire_patch[p] * dt
            ns.livestemn_storage_patch[p] = ns.livestemn_storage_patch[p] -
                nf.m_livestemn_storage_to_fire_patch[p] * dt
            ns.deadstemn_storage_patch[p] = ns.deadstemn_storage_patch[p] -
                nf.m_deadstemn_storage_to_fire_patch[p] * dt
            ns.livecrootn_storage_patch[p] = ns.livecrootn_storage_patch[p] -
                nf.m_livecrootn_storage_to_fire_patch[p] * dt
            ns.deadcrootn_storage_patch[p] = ns.deadcrootn_storage_patch[p] -
                nf.m_deadcrootn_storage_to_fire_patch[p] * dt

            ns.leafn_storage_patch[p] = ns.leafn_storage_patch[p] -
                nf.m_leafn_storage_to_litter_fire_patch[p] * dt
            ns.frootn_storage_patch[p] = ns.frootn_storage_patch[p] -
                nf.m_frootn_storage_to_litter_fire_patch[p] * dt
            ns.livestemn_storage_patch[p] = ns.livestemn_storage_patch[p] -
                nf.m_livestemn_storage_to_litter_fire_patch[p] * dt
            ns.deadstemn_storage_patch[p] = ns.deadstemn_storage_patch[p] -
                nf.m_deadstemn_storage_to_litter_fire_patch[p] * dt
            ns.livecrootn_storage_patch[p] = ns.livecrootn_storage_patch[p] -
                nf.m_livecrootn_storage_to_litter_fire_patch[p] * dt
            ns.deadcrootn_storage_patch[p] = ns.deadcrootn_storage_patch[p] -
                nf.m_deadcrootn_storage_to_litter_fire_patch[p] * dt

            # transfer pools
            ns.leafn_xfer_patch[p] = ns.leafn_xfer_patch[p] -
                nf.m_leafn_xfer_to_fire_patch[p] * dt
            ns.frootn_xfer_patch[p] = ns.frootn_xfer_patch[p] -
                nf.m_frootn_xfer_to_fire_patch[p] * dt
            ns.livestemn_xfer_patch[p] = ns.livestemn_xfer_patch[p] -
                nf.m_livestemn_xfer_to_fire_patch[p] * dt
            ns.deadstemn_xfer_patch[p] = ns.deadstemn_xfer_patch[p] -
                nf.m_deadstemn_xfer_to_fire_patch[p] * dt
            ns.livecrootn_xfer_patch[p] = ns.livecrootn_xfer_patch[p] -
                nf.m_livecrootn_xfer_to_fire_patch[p] * dt
            ns.deadcrootn_xfer_patch[p] = ns.deadcrootn_xfer_patch[p] -
                nf.m_deadcrootn_xfer_to_fire_patch[p] * dt

            ns.leafn_xfer_patch[p] = ns.leafn_xfer_patch[p] -
                nf.m_leafn_xfer_to_litter_fire_patch[p] * dt
            ns.frootn_xfer_patch[p] = ns.frootn_xfer_patch[p] -
                nf.m_frootn_xfer_to_litter_fire_patch[p] * dt
            ns.livestemn_xfer_patch[p] = ns.livestemn_xfer_patch[p] -
                nf.m_livestemn_xfer_to_litter_fire_patch[p] * dt
            ns.deadstemn_xfer_patch[p] = ns.deadstemn_xfer_patch[p] -
                nf.m_deadstemn_xfer_to_litter_fire_patch[p] * dt
            ns.livecrootn_xfer_patch[p] = ns.livecrootn_xfer_patch[p] -
                nf.m_livecrootn_xfer_to_litter_fire_patch[p] * dt
            ns.deadcrootn_xfer_patch[p] = ns.deadcrootn_xfer_patch[p] -
                nf.m_deadcrootn_xfer_to_litter_fire_patch[p] * dt

            # retranslocated N pool
            ns.retransn_patch[p] = ns.retransn_patch[p] -
                nf.m_retransn_to_fire_patch[p] * dt
            ns.retransn_patch[p] = ns.retransn_patch[p] -
                nf.m_retransn_to_litter_fire_patch[p] * dt
        else
            # NOTE: The equivalent changes for matrix code are in CNFireBase and CNFireLi2014 codes
        end
    end
end

"""
    n_state_update3!(ns_veg, nf_veg, ns_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, ndecomp_pools, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

On the radiation time step, update all prognostic nitrogen state variables
affected by fire fluxes.

Ported from `NStateUpdate3` in `CNNStateUpdate3Mod.F90`. Runs as three
KernelAbstractions kernels (two per soil column, one per patch); see the
kernel comments above for the GPU mapping. The CPU path is byte-identical.

Note: The matrix-CN code paths (`use_soil_matrixcn`) for soil matrix
solutions are omitted since they require the full matrix infrastructure
(matrix_Ninput). The matrix-CN veg path (`use_matrixcn`) is a no-op in
the Fortran original (comment only), so it is omitted here as well.
"""
function n_state_update3!(ns_veg::CNVegNitrogenStateData,
                           nf_veg::CNVegNitrogenFluxData,
                           ns_soil::SoilBiogeochemNitrogenStateData;
                           mask_soilc::AbstractVector{Bool},
                           mask_soilp::AbstractVector{Bool},
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
    _launch!(_nsu3_col_cwdlitr_kernel!, mask_soilc,
        ns_soil.decomp_npools_vr_col,
        nf_veg.fire_mortality_n_to_cwdn_col, nf_veg.m_n_to_litr_fire_col,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd, use_soil_matrixcn, dt)

    # --- Column loop: litter and CWD losses to fire ---
    if !use_soil_matrixcn
        _launch!(_nsu3_col_firelosses_kernel!, mask_soilc,
            ns_soil.decomp_npools_vr_col, nf_veg.m_decomp_npools_to_fire_vr_col,
            nlevdecomp, ndecomp_pools, dt)
    end

    # --- Patch loop: vegetation N state updates from fire (one thread per patch) ---
    ns = _nsu3_ns(ns_veg)
    nf = _nsu3_nf(nf_veg)
    _launch!(_nsu3_patch_kernel!, mask_soilp, ns, nf, use_matrixcn, dt)

    return nothing
end
