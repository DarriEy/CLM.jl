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
#
# GPU kernelization (Phase B C/N state-update cascade):
#   * Column loops -> `_csu3_col_kernel!` (one thread per soil column). The two
#     host column loops (fire mortality -> cwd/litter, then litter/CWD/all-pool
#     losses to fire) are fused per-column: for each column the additions over j
#     run first, then the per-pool losses over l/j run, preserving the host order
#     exactly. Internal j / l / i loops stay sequential in-thread, so the result is
#     byte-identical and race-free (each thread owns its column index).
#   * Patch loop -> `_csu3_patch_kernel!` (one thread per patch; every write is to
#     the patch's own index, no patch->column scatter). The ~20 carbon-state and
#     ~50 carbon-flux arrays the body touches are grouped into two immutable
#     `@adapt_structure` device-view bundles (`_CSU3CS{V}` / `_CSU3CF{V}`) so the
#     launch stays well under Metal's ~31-arg limit; bundle field names mirror the
#     state structs so the kernel body reads verbatim.
# Float64 literals are eltype-converted (none appear here beyond scalar arithmetic
# on `dt`, so the kernels carry no Float64 on a Float32-only backend); on Float64
# this is byte-identical.
# ==========================================================================

# --- Device-view bundle: CNVeg carbon STATE arrays the patch loop reads/writes ---
# `V` = the 1D patch vectors (all displayed/storage/transfer pools + gresp).
Base.@kwdef struct _CSU3CS{V}
    gresp_storage_patch::V; gresp_xfer_patch::V
    leafc_patch::V; frootc_patch::V
    livestemc_patch::V; deadstemc_patch::V
    livecrootc_patch::V; deadcrootc_patch::V
    leafc_storage_patch::V; frootc_storage_patch::V
    livestemc_storage_patch::V; deadstemc_storage_patch::V
    livecrootc_storage_patch::V; deadcrootc_storage_patch::V
    leafc_xfer_patch::V; frootc_xfer_patch::V
    livestemc_xfer_patch::V; deadstemc_xfer_patch::V
    livecrootc_xfer_patch::V; deadcrootc_xfer_patch::V
end
Adapt.@adapt_structure _CSU3CS

_csu3_cs(cs) = _CSU3CS(;
    gresp_storage_patch = cs.gresp_storage_patch, gresp_xfer_patch = cs.gresp_xfer_patch,
    leafc_patch = cs.leafc_patch, frootc_patch = cs.frootc_patch,
    livestemc_patch = cs.livestemc_patch, deadstemc_patch = cs.deadstemc_patch,
    livecrootc_patch = cs.livecrootc_patch, deadcrootc_patch = cs.deadcrootc_patch,
    leafc_storage_patch = cs.leafc_storage_patch, frootc_storage_patch = cs.frootc_storage_patch,
    livestemc_storage_patch = cs.livestemc_storage_patch,
    deadstemc_storage_patch = cs.deadstemc_storage_patch,
    livecrootc_storage_patch = cs.livecrootc_storage_patch,
    deadcrootc_storage_patch = cs.deadcrootc_storage_patch,
    leafc_xfer_patch = cs.leafc_xfer_patch, frootc_xfer_patch = cs.frootc_xfer_patch,
    livestemc_xfer_patch = cs.livestemc_xfer_patch, deadstemc_xfer_patch = cs.deadstemc_xfer_patch,
    livecrootc_xfer_patch = cs.livecrootc_xfer_patch, deadcrootc_xfer_patch = cs.deadcrootc_xfer_patch)

# --- Device-view bundle: CNVeg carbon FLUX arrays the patch loop reads ---
Base.@kwdef struct _CSU3CF{V}
    # gresp fire fluxes
    m_gresp_storage_to_fire_patch::V; m_gresp_storage_to_litter_fire_patch::V
    m_gresp_xfer_to_fire_patch::V; m_gresp_xfer_to_litter_fire_patch::V
    # displayed pools to fire / litter-fire
    m_leafc_to_fire_patch::V; m_leafc_to_litter_fire_patch::V
    m_frootc_to_fire_patch::V; m_frootc_to_litter_fire_patch::V
    m_livestemc_to_fire_patch::V; m_livestemc_to_litter_fire_patch::V
    m_livestemc_to_deadstemc_fire_patch::V
    m_deadstemc_to_fire_patch::V; m_deadstemc_to_litter_fire_patch::V
    m_livecrootc_to_fire_patch::V; m_livecrootc_to_litter_fire_patch::V
    m_livecrootc_to_deadcrootc_fire_patch::V
    m_deadcrootc_to_fire_patch::V; m_deadcrootc_to_litter_fire_patch::V
    # storage pools to fire / litter-fire
    m_leafc_storage_to_fire_patch::V; m_leafc_storage_to_litter_fire_patch::V
    m_frootc_storage_to_fire_patch::V; m_frootc_storage_to_litter_fire_patch::V
    m_livestemc_storage_to_fire_patch::V; m_livestemc_storage_to_litter_fire_patch::V
    m_deadstemc_storage_to_fire_patch::V; m_deadstemc_storage_to_litter_fire_patch::V
    m_livecrootc_storage_to_fire_patch::V; m_livecrootc_storage_to_litter_fire_patch::V
    m_deadcrootc_storage_to_fire_patch::V; m_deadcrootc_storage_to_litter_fire_patch::V
    # transfer pools to fire / litter-fire
    m_leafc_xfer_to_fire_patch::V; m_leafc_xfer_to_litter_fire_patch::V
    m_frootc_xfer_to_fire_patch::V; m_frootc_xfer_to_litter_fire_patch::V
    m_livestemc_xfer_to_fire_patch::V; m_livestemc_xfer_to_litter_fire_patch::V
    m_deadstemc_xfer_to_fire_patch::V; m_deadstemc_xfer_to_litter_fire_patch::V
    m_livecrootc_xfer_to_fire_patch::V; m_livecrootc_xfer_to_litter_fire_patch::V
    m_deadcrootc_xfer_to_fire_patch::V; m_deadcrootc_xfer_to_litter_fire_patch::V
end
Adapt.@adapt_structure _CSU3CF

_csu3_cf(cf) = _CSU3CF(;
    m_gresp_storage_to_fire_patch = cf.m_gresp_storage_to_fire_patch,
    m_gresp_storage_to_litter_fire_patch = cf.m_gresp_storage_to_litter_fire_patch,
    m_gresp_xfer_to_fire_patch = cf.m_gresp_xfer_to_fire_patch,
    m_gresp_xfer_to_litter_fire_patch = cf.m_gresp_xfer_to_litter_fire_patch,
    m_leafc_to_fire_patch = cf.m_leafc_to_fire_patch,
    m_leafc_to_litter_fire_patch = cf.m_leafc_to_litter_fire_patch,
    m_frootc_to_fire_patch = cf.m_frootc_to_fire_patch,
    m_frootc_to_litter_fire_patch = cf.m_frootc_to_litter_fire_patch,
    m_livestemc_to_fire_patch = cf.m_livestemc_to_fire_patch,
    m_livestemc_to_litter_fire_patch = cf.m_livestemc_to_litter_fire_patch,
    m_livestemc_to_deadstemc_fire_patch = cf.m_livestemc_to_deadstemc_fire_patch,
    m_deadstemc_to_fire_patch = cf.m_deadstemc_to_fire_patch,
    m_deadstemc_to_litter_fire_patch = cf.m_deadstemc_to_litter_fire_patch,
    m_livecrootc_to_fire_patch = cf.m_livecrootc_to_fire_patch,
    m_livecrootc_to_litter_fire_patch = cf.m_livecrootc_to_litter_fire_patch,
    m_livecrootc_to_deadcrootc_fire_patch = cf.m_livecrootc_to_deadcrootc_fire_patch,
    m_deadcrootc_to_fire_patch = cf.m_deadcrootc_to_fire_patch,
    m_deadcrootc_to_litter_fire_patch = cf.m_deadcrootc_to_litter_fire_patch,
    m_leafc_storage_to_fire_patch = cf.m_leafc_storage_to_fire_patch,
    m_leafc_storage_to_litter_fire_patch = cf.m_leafc_storage_to_litter_fire_patch,
    m_frootc_storage_to_fire_patch = cf.m_frootc_storage_to_fire_patch,
    m_frootc_storage_to_litter_fire_patch = cf.m_frootc_storage_to_litter_fire_patch,
    m_livestemc_storage_to_fire_patch = cf.m_livestemc_storage_to_fire_patch,
    m_livestemc_storage_to_litter_fire_patch = cf.m_livestemc_storage_to_litter_fire_patch,
    m_deadstemc_storage_to_fire_patch = cf.m_deadstemc_storage_to_fire_patch,
    m_deadstemc_storage_to_litter_fire_patch = cf.m_deadstemc_storage_to_litter_fire_patch,
    m_livecrootc_storage_to_fire_patch = cf.m_livecrootc_storage_to_fire_patch,
    m_livecrootc_storage_to_litter_fire_patch = cf.m_livecrootc_storage_to_litter_fire_patch,
    m_deadcrootc_storage_to_fire_patch = cf.m_deadcrootc_storage_to_fire_patch,
    m_deadcrootc_storage_to_litter_fire_patch = cf.m_deadcrootc_storage_to_litter_fire_patch,
    m_leafc_xfer_to_fire_patch = cf.m_leafc_xfer_to_fire_patch,
    m_leafc_xfer_to_litter_fire_patch = cf.m_leafc_xfer_to_litter_fire_patch,
    m_frootc_xfer_to_fire_patch = cf.m_frootc_xfer_to_fire_patch,
    m_frootc_xfer_to_litter_fire_patch = cf.m_frootc_xfer_to_litter_fire_patch,
    m_livestemc_xfer_to_fire_patch = cf.m_livestemc_xfer_to_fire_patch,
    m_livestemc_xfer_to_litter_fire_patch = cf.m_livestemc_xfer_to_litter_fire_patch,
    m_deadstemc_xfer_to_fire_patch = cf.m_deadstemc_xfer_to_fire_patch,
    m_deadstemc_xfer_to_litter_fire_patch = cf.m_deadstemc_xfer_to_litter_fire_patch,
    m_livecrootc_xfer_to_fire_patch = cf.m_livecrootc_xfer_to_fire_patch,
    m_livecrootc_xfer_to_litter_fire_patch = cf.m_livecrootc_xfer_to_litter_fire_patch,
    m_deadcrootc_xfer_to_fire_patch = cf.m_deadcrootc_xfer_to_fire_patch,
    m_deadcrootc_xfer_to_litter_fire_patch = cf.m_deadcrootc_xfer_to_litter_fire_patch)

# --- Kernel: column-level fire fluxes into soil decomposition pools (one thread per column) ---
# Fuses the two host column loops: per column, do the fire-mortality additions to
# cwd/litter over j first, then the per-pool fire losses over l/j. Order matches the
# host (the loss loop runs after the addition loop for every column).
@kernel function _csu3_col_kernel!(@Const(mask_soilc), decomp_cpools_vr_col,
        @Const(fire_mortality_c_to_cwdc_col), @Const(m_c_to_litr_fire_col),
        @Const(m_decomp_cpools_to_fire_vr_col),
        nlevdecomp::Int, ndecomp_pools::Int,
        i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        use_soil_matrixcn::Bool, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        # --- fire mortality fluxes to soil decomposition pools ---
        for j in 1:nlevdecomp
            if !use_soil_matrixcn
                # patch-level wood to column-level CWD (uncombusted wood)
                decomp_cpools_vr_col[c, j, i_cwd] =
                    decomp_cpools_vr_col[c, j, i_cwd] +
                    fire_mortality_c_to_cwdc_col[c, j] * dt

                # patch-level wood to column-level litter (uncombusted wood)
                for i in i_litr_min:i_litr_max
                    decomp_cpools_vr_col[c, j, i] =
                        decomp_cpools_vr_col[c, j, i] +
                        m_c_to_litr_fire_col[c, j, i] * dt
                end
            else
                # Matrix path: matrix_Cinput setup omitted (requires full matrix infrastructure)
            end
        end

        # --- litter and CWD losses to fire ---
        if !use_soil_matrixcn
            for l in 1:ndecomp_pools
                for j in 1:nlevdecomp
                    decomp_cpools_vr_col[c, j, l] =
                        decomp_cpools_vr_col[c, j, l] -
                        m_decomp_cpools_to_fire_vr_col[c, j, l] * dt
                end
            end
        end
    end
end

# --- Kernel: patch-level vegetation C state updates from fire (one thread per patch) ---
# `cs` / `cf` are the _CSU3CS / _CSU3CF device-view bundles; field names mirror the
# state structs so the body below is the host loop verbatim (cs_veg.->cs., cf_veg.->cf.).
@kernel function _csu3_patch_kernel!(@Const(mask_soilp), cs, cf,
        use_matrixcn::Bool, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        cs.gresp_storage_patch[p] = cs.gresp_storage_patch[p] -
            cf.m_gresp_storage_to_fire_patch[p] * dt
        cs.gresp_storage_patch[p] = cs.gresp_storage_patch[p] -
            cf.m_gresp_storage_to_litter_fire_patch[p] * dt
        cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] -
            cf.m_gresp_xfer_to_fire_patch[p] * dt
        cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] -
            cf.m_gresp_xfer_to_litter_fire_patch[p] * dt

        if !use_matrixcn
            # displayed pools
            cs.leafc_patch[p] = cs.leafc_patch[p] -
                cf.m_leafc_to_fire_patch[p] * dt
            cs.leafc_patch[p] = cs.leafc_patch[p] -
                cf.m_leafc_to_litter_fire_patch[p] * dt
            cs.frootc_patch[p] = cs.frootc_patch[p] -
                cf.m_frootc_to_fire_patch[p] * dt
            cs.frootc_patch[p] = cs.frootc_patch[p] -
                cf.m_frootc_to_litter_fire_patch[p] * dt
            cs.livestemc_patch[p] = cs.livestemc_patch[p] -
                cf.m_livestemc_to_fire_patch[p] * dt
            cs.livestemc_patch[p] = cs.livestemc_patch[p] -
                cf.m_livestemc_to_litter_fire_patch[p] * dt -
                cf.m_livestemc_to_deadstemc_fire_patch[p] * dt
            cs.deadstemc_patch[p] = cs.deadstemc_patch[p] -
                cf.m_deadstemc_to_fire_patch[p] * dt
            cs.deadstemc_patch[p] = cs.deadstemc_patch[p] -
                cf.m_deadstemc_to_litter_fire_patch[p] * dt +
                cf.m_livestemc_to_deadstemc_fire_patch[p] * dt
            cs.livecrootc_patch[p] = cs.livecrootc_patch[p] -
                cf.m_livecrootc_to_fire_patch[p] * dt
            cs.livecrootc_patch[p] = cs.livecrootc_patch[p] -
                cf.m_livecrootc_to_litter_fire_patch[p] * dt -
                cf.m_livecrootc_to_deadcrootc_fire_patch[p] * dt
            cs.deadcrootc_patch[p] = cs.deadcrootc_patch[p] -
                cf.m_deadcrootc_to_fire_patch[p] * dt
            cs.deadcrootc_patch[p] = cs.deadcrootc_patch[p] -
                cf.m_deadcrootc_to_litter_fire_patch[p] * dt +
                cf.m_livecrootc_to_deadcrootc_fire_patch[p] * dt

            # storage pools
            cs.leafc_storage_patch[p] = cs.leafc_storage_patch[p] -
                cf.m_leafc_storage_to_fire_patch[p] * dt
            cs.leafc_storage_patch[p] = cs.leafc_storage_patch[p] -
                cf.m_leafc_storage_to_litter_fire_patch[p] * dt
            cs.frootc_storage_patch[p] = cs.frootc_storage_patch[p] -
                cf.m_frootc_storage_to_fire_patch[p] * dt
            cs.frootc_storage_patch[p] = cs.frootc_storage_patch[p] -
                cf.m_frootc_storage_to_litter_fire_patch[p] * dt
            cs.livestemc_storage_patch[p] = cs.livestemc_storage_patch[p] -
                cf.m_livestemc_storage_to_fire_patch[p] * dt
            cs.livestemc_storage_patch[p] = cs.livestemc_storage_patch[p] -
                cf.m_livestemc_storage_to_litter_fire_patch[p] * dt
            cs.deadstemc_storage_patch[p] = cs.deadstemc_storage_patch[p] -
                cf.m_deadstemc_storage_to_fire_patch[p] * dt
            cs.deadstemc_storage_patch[p] = cs.deadstemc_storage_patch[p] -
                cf.m_deadstemc_storage_to_litter_fire_patch[p] * dt
            cs.livecrootc_storage_patch[p] = cs.livecrootc_storage_patch[p] -
                cf.m_livecrootc_storage_to_fire_patch[p] * dt
            cs.livecrootc_storage_patch[p] = cs.livecrootc_storage_patch[p] -
                cf.m_livecrootc_storage_to_litter_fire_patch[p] * dt
            cs.deadcrootc_storage_patch[p] = cs.deadcrootc_storage_patch[p] -
                cf.m_deadcrootc_storage_to_fire_patch[p] * dt
            cs.deadcrootc_storage_patch[p] = cs.deadcrootc_storage_patch[p] -
                cf.m_deadcrootc_storage_to_litter_fire_patch[p] * dt

            # transfer pools
            cs.leafc_xfer_patch[p] = cs.leafc_xfer_patch[p] -
                cf.m_leafc_xfer_to_fire_patch[p] * dt
            cs.leafc_xfer_patch[p] = cs.leafc_xfer_patch[p] -
                cf.m_leafc_xfer_to_litter_fire_patch[p] * dt
            cs.frootc_xfer_patch[p] = cs.frootc_xfer_patch[p] -
                cf.m_frootc_xfer_to_fire_patch[p] * dt
            cs.frootc_xfer_patch[p] = cs.frootc_xfer_patch[p] -
                cf.m_frootc_xfer_to_litter_fire_patch[p] * dt
            cs.livestemc_xfer_patch[p] = cs.livestemc_xfer_patch[p] -
                cf.m_livestemc_xfer_to_fire_patch[p] * dt
            cs.livestemc_xfer_patch[p] = cs.livestemc_xfer_patch[p] -
                cf.m_livestemc_xfer_to_litter_fire_patch[p] * dt
            cs.deadstemc_xfer_patch[p] = cs.deadstemc_xfer_patch[p] -
                cf.m_deadstemc_xfer_to_fire_patch[p] * dt
            cs.deadstemc_xfer_patch[p] = cs.deadstemc_xfer_patch[p] -
                cf.m_deadstemc_xfer_to_litter_fire_patch[p] * dt
            cs.livecrootc_xfer_patch[p] = cs.livecrootc_xfer_patch[p] -
                cf.m_livecrootc_xfer_to_fire_patch[p] * dt
            cs.livecrootc_xfer_patch[p] = cs.livecrootc_xfer_patch[p] -
                cf.m_livecrootc_xfer_to_litter_fire_patch[p] * dt
            cs.deadcrootc_xfer_patch[p] = cs.deadcrootc_xfer_patch[p] -
                cf.m_deadcrootc_xfer_to_fire_patch[p] * dt
            cs.deadcrootc_xfer_patch[p] = cs.deadcrootc_xfer_patch[p] -
                cf.m_deadcrootc_xfer_to_litter_fire_patch[p] * dt
        else
            # NOTE: The equivalent changes for matrix code are in CNFireBase and CNFireLi2014 codes
        end
    end
end

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

Ported from `CStateUpdate3` in `CNCStateUpdate3Mod.F90`. Runs as two
KernelAbstractions kernels (one per soil column, one per patch); see the
kernel comments above for the GPU mapping. The CPU path is byte-identical.

Note: The matrix-CN code paths (`use_soil_matrixcn`) for soil matrix
solutions are omitted since they require the full matrix infrastructure
(matrix_Cinput). The matrix-CN veg path (`use_matrixcn`) is a no-op in
the Fortran original (comment only), so it is omitted here as well.
"""
function c_state_update3!(cs_veg::CNVegCarbonStateData,
                           cf_veg::CNVegCarbonFluxData,
                           cs_soil::SoilBiogeochemCarbonStateData;
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

    # --- Column loop: fire fluxes into soil decomposition pools (one thread per column) ---
    _launch!(_csu3_col_kernel!, mask_soilc, cs_soil.decomp_cpools_vr_col,
        cf_veg.fire_mortality_c_to_cwdc_col, cf_veg.m_c_to_litr_fire_col,
        cf_veg.m_decomp_cpools_to_fire_vr_col,
        nlevdecomp, ndecomp_pools, i_litr_min, i_litr_max, i_cwd,
        use_soil_matrixcn, dt)

    # --- Patch loop: vegetation C state updates from fire (one thread per patch) ---
    cs = _csu3_cs(cs_veg)
    cf = _csu3_cf(cf_veg)
    _launch!(_csu3_patch_kernel!, mask_soilp, cs, cf, use_matrixcn, dt)

    return nothing
end
