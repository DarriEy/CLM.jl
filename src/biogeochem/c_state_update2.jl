# ==========================================================================
# Ported from: src/biogeochem/CNCStateUpdate2Mod.F90
# Carbon state variable update, mortality fluxes.
#
# When the matrix solution is being used (use_matrixcn and use_soil_matrixcn)
# only some state updates are done here; other state updates happen after the
# matrix is solved in VegMatrix and SoilMatrix.
#
# Public functions:
#   c_state_update2!  — Gap-phase mortality C state update
#   c_state_update2h! — Harvest mortality C state update
#   c_state_update2g! — Gross unrepresented landcover change mortality C state update
#
# GPU kernelization (Phase B; same pattern as c_state_update1.jl):
#   * Each column j/i decomposition loop -> a per-column kernel (one thread per
#     soil column; the internal j/i RMW accumulation into decomp_cpools_vr_col
#     runs sequentially in-thread, byte-identical and race-free — each thread
#     owns its column).
#   * Each per-patch veg-state loop -> a per-patch kernel (one thread per patch;
#     every write is to the patch's own index, no patch->column scatter). The
#     ~20 carbon-state + ~20 carbon-flux arrays each body touches are grouped
#     into immutable `@adapt_structure` device-view bundles (one cs + one cf per
#     function: `_CSU2CS/_CSU2CF`, `_CSU2HCS/_CSU2HCF`, `_CSU2GCS/_CSU2GCF`) so
#     the launch stays well under Metal's ~31-arg limit; the bundle field names
#     mirror the state structs so the body reads verbatim.
# Float64 literals are eltype-converted (`T(...)`) so the kernels carry no
# Float64 on a Float32-only backend; on Float64 this is byte-identical.
# ==========================================================================

# ---------------------------------------------------------------------------
# c_state_update2! — Gap-phase mortality
# ---------------------------------------------------------------------------

# --- Device-view bundle: CNVeg carbon STATE arrays the gap-mortality patch loop touches ---
Base.@kwdef struct _CSU2CS{V}
    leafc_patch::V; leafc_storage_patch::V; leafc_xfer_patch::V
    frootc_patch::V; frootc_storage_patch::V; frootc_xfer_patch::V
    livestemc_patch::V; livestemc_storage_patch::V; livestemc_xfer_patch::V
    deadstemc_patch::V; deadstemc_storage_patch::V; deadstemc_xfer_patch::V
    livecrootc_patch::V; livecrootc_storage_patch::V; livecrootc_xfer_patch::V
    deadcrootc_patch::V; deadcrootc_storage_patch::V; deadcrootc_xfer_patch::V
    gresp_storage_patch::V; gresp_xfer_patch::V
end
Adapt.@adapt_structure _CSU2CS

_csu2_cs(cs) = _CSU2CS(;
    leafc_patch = cs.leafc_patch, leafc_storage_patch = cs.leafc_storage_patch,
    leafc_xfer_patch = cs.leafc_xfer_patch,
    frootc_patch = cs.frootc_patch, frootc_storage_patch = cs.frootc_storage_patch,
    frootc_xfer_patch = cs.frootc_xfer_patch,
    livestemc_patch = cs.livestemc_patch, livestemc_storage_patch = cs.livestemc_storage_patch,
    livestemc_xfer_patch = cs.livestemc_xfer_patch,
    deadstemc_patch = cs.deadstemc_patch, deadstemc_storage_patch = cs.deadstemc_storage_patch,
    deadstemc_xfer_patch = cs.deadstemc_xfer_patch,
    livecrootc_patch = cs.livecrootc_patch, livecrootc_storage_patch = cs.livecrootc_storage_patch,
    livecrootc_xfer_patch = cs.livecrootc_xfer_patch,
    deadcrootc_patch = cs.deadcrootc_patch, deadcrootc_storage_patch = cs.deadcrootc_storage_patch,
    deadcrootc_xfer_patch = cs.deadcrootc_xfer_patch,
    gresp_storage_patch = cs.gresp_storage_patch, gresp_xfer_patch = cs.gresp_xfer_patch)

# --- Device-view bundle: CNVeg carbon FLUX arrays the gap-mortality patch loop reads ---
Base.@kwdef struct _CSU2CF{V}
    m_gresp_storage_to_litter_patch::V; m_gresp_xfer_to_litter_patch::V
    m_leafc_to_litter_patch::V; m_frootc_to_litter_patch::V
    m_livestemc_to_litter_patch::V; m_deadstemc_to_litter_patch::V
    m_livecrootc_to_litter_patch::V; m_deadcrootc_to_litter_patch::V
    m_leafc_storage_to_litter_patch::V; m_frootc_storage_to_litter_patch::V
    m_livestemc_storage_to_litter_patch::V; m_deadstemc_storage_to_litter_patch::V
    m_livecrootc_storage_to_litter_patch::V; m_deadcrootc_storage_to_litter_patch::V
    m_leafc_xfer_to_litter_patch::V; m_frootc_xfer_to_litter_patch::V
    m_livestemc_xfer_to_litter_patch::V; m_deadstemc_xfer_to_litter_patch::V
    m_livecrootc_xfer_to_litter_patch::V; m_deadcrootc_xfer_to_litter_patch::V
end
Adapt.@adapt_structure _CSU2CF

_csu2_cf(cf) = _CSU2CF(;
    m_gresp_storage_to_litter_patch = cf.m_gresp_storage_to_litter_patch,
    m_gresp_xfer_to_litter_patch = cf.m_gresp_xfer_to_litter_patch,
    m_leafc_to_litter_patch = cf.m_leafc_to_litter_patch,
    m_frootc_to_litter_patch = cf.m_frootc_to_litter_patch,
    m_livestemc_to_litter_patch = cf.m_livestemc_to_litter_patch,
    m_deadstemc_to_litter_patch = cf.m_deadstemc_to_litter_patch,
    m_livecrootc_to_litter_patch = cf.m_livecrootc_to_litter_patch,
    m_deadcrootc_to_litter_patch = cf.m_deadcrootc_to_litter_patch,
    m_leafc_storage_to_litter_patch = cf.m_leafc_storage_to_litter_patch,
    m_frootc_storage_to_litter_patch = cf.m_frootc_storage_to_litter_patch,
    m_livestemc_storage_to_litter_patch = cf.m_livestemc_storage_to_litter_patch,
    m_deadstemc_storage_to_litter_patch = cf.m_deadstemc_storage_to_litter_patch,
    m_livecrootc_storage_to_litter_patch = cf.m_livecrootc_storage_to_litter_patch,
    m_deadcrootc_storage_to_litter_patch = cf.m_deadcrootc_storage_to_litter_patch,
    m_leafc_xfer_to_litter_patch = cf.m_leafc_xfer_to_litter_patch,
    m_frootc_xfer_to_litter_patch = cf.m_frootc_xfer_to_litter_patch,
    m_livestemc_xfer_to_litter_patch = cf.m_livestemc_xfer_to_litter_patch,
    m_deadstemc_xfer_to_litter_patch = cf.m_deadstemc_xfer_to_litter_patch,
    m_livecrootc_xfer_to_litter_patch = cf.m_livecrootc_xfer_to_litter_patch,
    m_deadcrootc_xfer_to_litter_patch = cf.m_deadcrootc_xfer_to_litter_patch)

# --- Kernel: gap-mortality fluxes into soil decomposition pools (one thread per column) ---
@kernel function _csu2_col_kernel!(@Const(mask_soilc), decomp_cpools_vr_col,
        @Const(gap_mortality_c_to_litr_c_col), @Const(gap_mortality_c_to_cwdc_col),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        use_soil_matrixcn::Bool, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        if !use_soil_matrixcn
            for j in 1:nlevdecomp
                for i in i_litr_min:i_litr_max
                    decomp_cpools_vr_col[c, j, i] =
                        decomp_cpools_vr_col[c, j, i] +
                        gap_mortality_c_to_litr_c_col[c, j, i] * dt
                end
                # Currently i_cwd != i_litr_max + 1 if not fates and
                # i_cwd = 0 if fates, so not including in the i-loop
                decomp_cpools_vr_col[c, j, i_cwd] =
                    decomp_cpools_vr_col[c, j, i_cwd] +
                    gap_mortality_c_to_cwdc_col[c, j] * dt
            end
        else
            # Matrix path: matrix_Cinput setup omitted (requires full matrix infrastructure)
        end
    end
end

# --- Kernel: gap-mortality veg C state updates (one thread per patch) ---
@kernel function _csu2_patch_kernel!(@Const(mask_soilp), cs, cf, use_matrixcn::Bool, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        cs.gresp_storage_patch[p] = cs.gresp_storage_patch[p] -
            cf.m_gresp_storage_to_litter_patch[p] * dt
        cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] -
            cf.m_gresp_xfer_to_litter_patch[p] * dt

        if !use_matrixcn
            # displayed pools
            cs.leafc_patch[p] = cs.leafc_patch[p] -
                cf.m_leafc_to_litter_patch[p] * dt
            cs.frootc_patch[p] = cs.frootc_patch[p] -
                cf.m_frootc_to_litter_patch[p] * dt
            cs.livestemc_patch[p] = cs.livestemc_patch[p] -
                cf.m_livestemc_to_litter_patch[p] * dt
            cs.deadstemc_patch[p] = cs.deadstemc_patch[p] -
                cf.m_deadstemc_to_litter_patch[p] * dt
            cs.livecrootc_patch[p] = cs.livecrootc_patch[p] -
                cf.m_livecrootc_to_litter_patch[p] * dt
            cs.deadcrootc_patch[p] = cs.deadcrootc_patch[p] -
                cf.m_deadcrootc_to_litter_patch[p] * dt

            # storage pools
            cs.leafc_storage_patch[p] = cs.leafc_storage_patch[p] -
                cf.m_leafc_storage_to_litter_patch[p] * dt
            cs.frootc_storage_patch[p] = cs.frootc_storage_patch[p] -
                cf.m_frootc_storage_to_litter_patch[p] * dt
            cs.livestemc_storage_patch[p] = cs.livestemc_storage_patch[p] -
                cf.m_livestemc_storage_to_litter_patch[p] * dt
            cs.deadstemc_storage_patch[p] = cs.deadstemc_storage_patch[p] -
                cf.m_deadstemc_storage_to_litter_patch[p] * dt
            cs.livecrootc_storage_patch[p] = cs.livecrootc_storage_patch[p] -
                cf.m_livecrootc_storage_to_litter_patch[p] * dt
            cs.deadcrootc_storage_patch[p] = cs.deadcrootc_storage_patch[p] -
                cf.m_deadcrootc_storage_to_litter_patch[p] * dt

            # transfer pools
            cs.leafc_xfer_patch[p] = cs.leafc_xfer_patch[p] -
                cf.m_leafc_xfer_to_litter_patch[p] * dt
            cs.frootc_xfer_patch[p] = cs.frootc_xfer_patch[p] -
                cf.m_frootc_xfer_to_litter_patch[p] * dt
            cs.livestemc_xfer_patch[p] = cs.livestemc_xfer_patch[p] -
                cf.m_livestemc_xfer_to_litter_patch[p] * dt
            cs.deadstemc_xfer_patch[p] = cs.deadstemc_xfer_patch[p] -
                cf.m_deadstemc_xfer_to_litter_patch[p] * dt
            cs.livecrootc_xfer_patch[p] = cs.livecrootc_xfer_patch[p] -
                cf.m_livecrootc_xfer_to_litter_patch[p] * dt
            cs.deadcrootc_xfer_patch[p] = cs.deadcrootc_xfer_patch[p] -
                cf.m_deadcrootc_xfer_to_litter_patch[p] * dt
        else
            # Matrix version: handled in CNGapMortality
        end
    end
end

"""
    c_state_update2!(cs_veg, cf_veg, cs_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

On the radiation time step, update all prognostic carbon state variables
affected by gap-phase mortality fluxes.

Ported from `CStateUpdate2` in `CNCStateUpdate2Mod.F90`. Runs as two
KernelAbstractions kernels (one per soil column, one per patch); the CPU
path is byte-identical.

Note: The matrix-CN code paths (`use_soil_matrixcn`) for soil matrix
solutions are omitted since they require the full matrix infrastructure
(matrix_Cinput). The matrix-CN veg path (`use_matrixcn`) is a no-op in
the Fortran original (comment only), so it is omitted here as well.
"""
function c_state_update2!(cs_veg::CNVegCarbonStateData,
                           cf_veg::CNVegCarbonFluxData,
                           cs_soil::SoilBiogeochemCarbonStateData;
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
    _launch!(_csu2_col_kernel!, mask_soilc, cs_soil.decomp_cpools_vr_col,
        cf_veg.gap_mortality_c_to_litr_c_col, cf_veg.gap_mortality_c_to_cwdc_col,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd, use_soil_matrixcn, dt)

    # --- Patch loop: vegetation C state updates from gap-phase mortality ---
    cs = _csu2_cs(cs_veg)
    cf = _csu2_cf(cf_veg)
    _launch!(_csu2_patch_kernel!, mask_soilp, cs, cf, use_matrixcn, dt)

    return nothing
end

# ---------------------------------------------------------------------------
# c_state_update2h! — Harvest mortality
# ---------------------------------------------------------------------------

# --- Device-view bundle: CNVeg carbon STATE arrays the harvest patch loop touches ---
Base.@kwdef struct _CSU2HCS{V}
    xsmrpool_patch::V; gresp_storage_patch::V; gresp_xfer_patch::V
    leafc_patch::V; leafc_storage_patch::V; leafc_xfer_patch::V
    frootc_patch::V; frootc_storage_patch::V; frootc_xfer_patch::V
    livestemc_patch::V; livestemc_storage_patch::V; livestemc_xfer_patch::V
    deadstemc_patch::V; deadstemc_storage_patch::V; deadstemc_xfer_patch::V
    livecrootc_patch::V; livecrootc_storage_patch::V; livecrootc_xfer_patch::V
    deadcrootc_patch::V; deadcrootc_storage_patch::V; deadcrootc_xfer_patch::V
end
Adapt.@adapt_structure _CSU2HCS

_csu2h_cs(cs) = _CSU2HCS(;
    xsmrpool_patch = cs.xsmrpool_patch,
    gresp_storage_patch = cs.gresp_storage_patch, gresp_xfer_patch = cs.gresp_xfer_patch,
    leafc_patch = cs.leafc_patch, leafc_storage_patch = cs.leafc_storage_patch,
    leafc_xfer_patch = cs.leafc_xfer_patch,
    frootc_patch = cs.frootc_patch, frootc_storage_patch = cs.frootc_storage_patch,
    frootc_xfer_patch = cs.frootc_xfer_patch,
    livestemc_patch = cs.livestemc_patch, livestemc_storage_patch = cs.livestemc_storage_patch,
    livestemc_xfer_patch = cs.livestemc_xfer_patch,
    deadstemc_patch = cs.deadstemc_patch, deadstemc_storage_patch = cs.deadstemc_storage_patch,
    deadstemc_xfer_patch = cs.deadstemc_xfer_patch,
    livecrootc_patch = cs.livecrootc_patch, livecrootc_storage_patch = cs.livecrootc_storage_patch,
    livecrootc_xfer_patch = cs.livecrootc_xfer_patch,
    deadcrootc_patch = cs.deadcrootc_patch, deadcrootc_storage_patch = cs.deadcrootc_storage_patch,
    deadcrootc_xfer_patch = cs.deadcrootc_xfer_patch)

# --- Device-view bundle: CNVeg carbon FLUX arrays the harvest patch loop reads ---
Base.@kwdef struct _CSU2HCF{V}
    hrv_xsmrpool_to_atm_patch::V
    hrv_gresp_storage_to_litter_patch::V; hrv_gresp_xfer_to_litter_patch::V
    hrv_leafc_to_litter_patch::V; hrv_frootc_to_litter_patch::V
    hrv_livestemc_to_litter_patch::V; wood_harvestc_patch::V
    hrv_livecrootc_to_litter_patch::V; hrv_deadcrootc_to_litter_patch::V
    hrv_leafc_storage_to_litter_patch::V; hrv_frootc_storage_to_litter_patch::V
    hrv_livestemc_storage_to_litter_patch::V; hrv_deadstemc_storage_to_litter_patch::V
    hrv_livecrootc_storage_to_litter_patch::V; hrv_deadcrootc_storage_to_litter_patch::V
    hrv_leafc_xfer_to_litter_patch::V; hrv_frootc_xfer_to_litter_patch::V
    hrv_livestemc_xfer_to_litter_patch::V; hrv_deadstemc_xfer_to_litter_patch::V
    hrv_livecrootc_xfer_to_litter_patch::V; hrv_deadcrootc_xfer_to_litter_patch::V
end
Adapt.@adapt_structure _CSU2HCF

_csu2h_cf(cf) = _CSU2HCF(;
    hrv_xsmrpool_to_atm_patch = cf.hrv_xsmrpool_to_atm_patch,
    hrv_gresp_storage_to_litter_patch = cf.hrv_gresp_storage_to_litter_patch,
    hrv_gresp_xfer_to_litter_patch = cf.hrv_gresp_xfer_to_litter_patch,
    hrv_leafc_to_litter_patch = cf.hrv_leafc_to_litter_patch,
    hrv_frootc_to_litter_patch = cf.hrv_frootc_to_litter_patch,
    hrv_livestemc_to_litter_patch = cf.hrv_livestemc_to_litter_patch,
    wood_harvestc_patch = cf.wood_harvestc_patch,
    hrv_livecrootc_to_litter_patch = cf.hrv_livecrootc_to_litter_patch,
    hrv_deadcrootc_to_litter_patch = cf.hrv_deadcrootc_to_litter_patch,
    hrv_leafc_storage_to_litter_patch = cf.hrv_leafc_storage_to_litter_patch,
    hrv_frootc_storage_to_litter_patch = cf.hrv_frootc_storage_to_litter_patch,
    hrv_livestemc_storage_to_litter_patch = cf.hrv_livestemc_storage_to_litter_patch,
    hrv_deadstemc_storage_to_litter_patch = cf.hrv_deadstemc_storage_to_litter_patch,
    hrv_livecrootc_storage_to_litter_patch = cf.hrv_livecrootc_storage_to_litter_patch,
    hrv_deadcrootc_storage_to_litter_patch = cf.hrv_deadcrootc_storage_to_litter_patch,
    hrv_leafc_xfer_to_litter_patch = cf.hrv_leafc_xfer_to_litter_patch,
    hrv_frootc_xfer_to_litter_patch = cf.hrv_frootc_xfer_to_litter_patch,
    hrv_livestemc_xfer_to_litter_patch = cf.hrv_livestemc_xfer_to_litter_patch,
    hrv_deadstemc_xfer_to_litter_patch = cf.hrv_deadstemc_xfer_to_litter_patch,
    hrv_livecrootc_xfer_to_litter_patch = cf.hrv_livecrootc_xfer_to_litter_patch,
    hrv_deadcrootc_xfer_to_litter_patch = cf.hrv_deadcrootc_xfer_to_litter_patch)

# --- Kernel: harvest fluxes into soil decomposition pools (one thread per column) ---
@kernel function _csu2h_col_kernel!(@Const(mask_soilc), decomp_cpools_vr_col,
        @Const(harvest_c_to_litr_c_col), @Const(harvest_c_to_cwdc_col),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        use_soil_matrixcn::Bool, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        if !use_soil_matrixcn
            for j in 1:nlevdecomp
                for i in i_litr_min:i_litr_max
                    decomp_cpools_vr_col[c, j, i] =
                        decomp_cpools_vr_col[c, j, i] +
                        harvest_c_to_litr_c_col[c, j, i] * dt
                end
                # Currently i_cwd != i_litr_max + 1 if not fates and
                # i_cwd = 0 if fates, so not including in the i-loop
                decomp_cpools_vr_col[c, j, i_cwd] =
                    decomp_cpools_vr_col[c, j, i_cwd] +
                    harvest_c_to_cwdc_col[c, j] * dt
            end
            # wood to product pools - states updated in CNProducts
        else
            # Matrix path: matrix_Cinput setup omitted (requires full matrix infrastructure)
        end
    end
end

# --- Kernel: harvest veg C state updates (one thread per patch) ---
@kernel function _csu2h_patch_kernel!(@Const(mask_soilp), cs, cf, use_matrixcn::Bool, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        # xsmrpool
        cs.xsmrpool_patch[p] = cs.xsmrpool_patch[p] -
            cf.hrv_xsmrpool_to_atm_patch[p] * dt

        # storage pools
        cs.gresp_storage_patch[p] = cs.gresp_storage_patch[p] -
            cf.hrv_gresp_storage_to_litter_patch[p] * dt

        # transfer pools
        cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] -
            cf.hrv_gresp_xfer_to_litter_patch[p] * dt

        if !use_matrixcn
            # displayed pools
            cs.leafc_patch[p] = cs.leafc_patch[p] -
                cf.hrv_leafc_to_litter_patch[p] * dt
            cs.frootc_patch[p] = cs.frootc_patch[p] -
                cf.hrv_frootc_to_litter_patch[p] * dt
            cs.livestemc_patch[p] = cs.livestemc_patch[p] -
                cf.hrv_livestemc_to_litter_patch[p] * dt
            cs.deadstemc_patch[p] = cs.deadstemc_patch[p] -
                cf.wood_harvestc_patch[p] * dt
            cs.livecrootc_patch[p] = cs.livecrootc_patch[p] -
                cf.hrv_livecrootc_to_litter_patch[p] * dt
            cs.deadcrootc_patch[p] = cs.deadcrootc_patch[p] -
                cf.hrv_deadcrootc_to_litter_patch[p] * dt

            # storage pools
            cs.leafc_storage_patch[p] = cs.leafc_storage_patch[p] -
                cf.hrv_leafc_storage_to_litter_patch[p] * dt
            cs.frootc_storage_patch[p] = cs.frootc_storage_patch[p] -
                cf.hrv_frootc_storage_to_litter_patch[p] * dt
            cs.livestemc_storage_patch[p] = cs.livestemc_storage_patch[p] -
                cf.hrv_livestemc_storage_to_litter_patch[p] * dt
            cs.deadstemc_storage_patch[p] = cs.deadstemc_storage_patch[p] -
                cf.hrv_deadstemc_storage_to_litter_patch[p] * dt
            cs.livecrootc_storage_patch[p] = cs.livecrootc_storage_patch[p] -
                cf.hrv_livecrootc_storage_to_litter_patch[p] * dt
            cs.deadcrootc_storage_patch[p] = cs.deadcrootc_storage_patch[p] -
                cf.hrv_deadcrootc_storage_to_litter_patch[p] * dt

            # transfer pools
            cs.leafc_xfer_patch[p] = cs.leafc_xfer_patch[p] -
                cf.hrv_leafc_xfer_to_litter_patch[p] * dt
            cs.frootc_xfer_patch[p] = cs.frootc_xfer_patch[p] -
                cf.hrv_frootc_xfer_to_litter_patch[p] * dt
            cs.livestemc_xfer_patch[p] = cs.livestemc_xfer_patch[p] -
                cf.hrv_livestemc_xfer_to_litter_patch[p] * dt
            cs.deadstemc_xfer_patch[p] = cs.deadstemc_xfer_patch[p] -
                cf.hrv_deadstemc_xfer_to_litter_patch[p] * dt
            cs.livecrootc_xfer_patch[p] = cs.livecrootc_xfer_patch[p] -
                cf.hrv_livecrootc_xfer_to_litter_patch[p] * dt
            cs.deadcrootc_xfer_patch[p] = cs.deadcrootc_xfer_patch[p] -
                cf.hrv_deadcrootc_xfer_to_litter_patch[p] * dt
        else
            # Matrix version: handled in CNHarvest
        end
    end
end

"""
    c_state_update2h!(cs_veg, cf_veg, cs_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

Update all prognostic carbon state variables affected by harvest
mortality fluxes.

Ported from `CStateUpdate2h` in `CNCStateUpdate2Mod.F90`. Runs as two
KernelAbstractions kernels (one per soil column, one per patch); the CPU
path is byte-identical.

Note: The matrix-CN code paths are omitted (require full matrix infrastructure).
Wood to product pools state updates happen in CNProducts (not here).
"""
function c_state_update2h!(cs_veg::CNVegCarbonStateData,
                            cf_veg::CNVegCarbonFluxData,
                            cs_soil::SoilBiogeochemCarbonStateData;
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
    _launch!(_csu2h_col_kernel!, mask_soilc, cs_soil.decomp_cpools_vr_col,
        cf_veg.harvest_c_to_litr_c_col, cf_veg.harvest_c_to_cwdc_col,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd, use_soil_matrixcn, dt)

    # --- Patch loop: vegetation C state updates from harvest mortality ---
    cs = _csu2h_cs(cs_veg)
    cf = _csu2h_cf(cf_veg)
    _launch!(_csu2h_patch_kernel!, mask_soilp, cs, cf, use_matrixcn, dt)

    return nothing
end

# ---------------------------------------------------------------------------
# c_state_update2g! — Gross unrepresented landcover change mortality
# ---------------------------------------------------------------------------

# --- Device-view bundle: CNVeg carbon STATE arrays the gru patch loop touches ---
Base.@kwdef struct _CSU2GCS{V}
    xsmrpool_patch::V; gresp_storage_patch::V; gresp_xfer_patch::V
    leafc_patch::V; leafc_storage_patch::V; leafc_xfer_patch::V
    frootc_patch::V; frootc_storage_patch::V; frootc_xfer_patch::V
    livestemc_patch::V; livestemc_storage_patch::V; livestemc_xfer_patch::V
    deadstemc_patch::V; deadstemc_storage_patch::V; deadstemc_xfer_patch::V
    livecrootc_patch::V; livecrootc_storage_patch::V; livecrootc_xfer_patch::V
    deadcrootc_patch::V; deadcrootc_storage_patch::V; deadcrootc_xfer_patch::V
end
Adapt.@adapt_structure _CSU2GCS

_csu2g_cs(cs) = _CSU2GCS(;
    xsmrpool_patch = cs.xsmrpool_patch,
    gresp_storage_patch = cs.gresp_storage_patch, gresp_xfer_patch = cs.gresp_xfer_patch,
    leafc_patch = cs.leafc_patch, leafc_storage_patch = cs.leafc_storage_patch,
    leafc_xfer_patch = cs.leafc_xfer_patch,
    frootc_patch = cs.frootc_patch, frootc_storage_patch = cs.frootc_storage_patch,
    frootc_xfer_patch = cs.frootc_xfer_patch,
    livestemc_patch = cs.livestemc_patch, livestemc_storage_patch = cs.livestemc_storage_patch,
    livestemc_xfer_patch = cs.livestemc_xfer_patch,
    deadstemc_patch = cs.deadstemc_patch, deadstemc_storage_patch = cs.deadstemc_storage_patch,
    deadstemc_xfer_patch = cs.deadstemc_xfer_patch,
    livecrootc_patch = cs.livecrootc_patch, livecrootc_storage_patch = cs.livecrootc_storage_patch,
    livecrootc_xfer_patch = cs.livecrootc_xfer_patch,
    deadcrootc_patch = cs.deadcrootc_patch, deadcrootc_storage_patch = cs.deadcrootc_storage_patch,
    deadcrootc_xfer_patch = cs.deadcrootc_xfer_patch)

# --- Device-view bundle: CNVeg carbon FLUX arrays the gru patch loop reads ---
Base.@kwdef struct _CSU2GCF{V}
    gru_xsmrpool_to_atm_patch::V
    gru_gresp_storage_to_atm_patch::V; gru_gresp_xfer_to_atm_patch::V
    gru_leafc_to_litter_patch::V; gru_frootc_to_litter_patch::V
    gru_livestemc_to_atm_patch::V; gru_deadstemc_to_atm_patch::V
    gru_wood_productc_gain_patch::V
    gru_livecrootc_to_litter_patch::V; gru_deadcrootc_to_litter_patch::V
    gru_leafc_storage_to_atm_patch::V; gru_frootc_storage_to_atm_patch::V
    gru_livestemc_storage_to_atm_patch::V; gru_deadstemc_storage_to_atm_patch::V
    gru_livecrootc_storage_to_atm_patch::V; gru_deadcrootc_storage_to_atm_patch::V
    gru_leafc_xfer_to_atm_patch::V; gru_frootc_xfer_to_atm_patch::V
    gru_livestemc_xfer_to_atm_patch::V; gru_deadstemc_xfer_to_atm_patch::V
    gru_livecrootc_xfer_to_atm_patch::V; gru_deadcrootc_xfer_to_atm_patch::V
end
Adapt.@adapt_structure _CSU2GCF

_csu2g_cf(cf) = _CSU2GCF(;
    gru_xsmrpool_to_atm_patch = cf.gru_xsmrpool_to_atm_patch,
    gru_gresp_storage_to_atm_patch = cf.gru_gresp_storage_to_atm_patch,
    gru_gresp_xfer_to_atm_patch = cf.gru_gresp_xfer_to_atm_patch,
    gru_leafc_to_litter_patch = cf.gru_leafc_to_litter_patch,
    gru_frootc_to_litter_patch = cf.gru_frootc_to_litter_patch,
    gru_livestemc_to_atm_patch = cf.gru_livestemc_to_atm_patch,
    gru_deadstemc_to_atm_patch = cf.gru_deadstemc_to_atm_patch,
    gru_wood_productc_gain_patch = cf.gru_wood_productc_gain_patch,
    gru_livecrootc_to_litter_patch = cf.gru_livecrootc_to_litter_patch,
    gru_deadcrootc_to_litter_patch = cf.gru_deadcrootc_to_litter_patch,
    gru_leafc_storage_to_atm_patch = cf.gru_leafc_storage_to_atm_patch,
    gru_frootc_storage_to_atm_patch = cf.gru_frootc_storage_to_atm_patch,
    gru_livestemc_storage_to_atm_patch = cf.gru_livestemc_storage_to_atm_patch,
    gru_deadstemc_storage_to_atm_patch = cf.gru_deadstemc_storage_to_atm_patch,
    gru_livecrootc_storage_to_atm_patch = cf.gru_livecrootc_storage_to_atm_patch,
    gru_deadcrootc_storage_to_atm_patch = cf.gru_deadcrootc_storage_to_atm_patch,
    gru_leafc_xfer_to_atm_patch = cf.gru_leafc_xfer_to_atm_patch,
    gru_frootc_xfer_to_atm_patch = cf.gru_frootc_xfer_to_atm_patch,
    gru_livestemc_xfer_to_atm_patch = cf.gru_livestemc_xfer_to_atm_patch,
    gru_deadstemc_xfer_to_atm_patch = cf.gru_deadstemc_xfer_to_atm_patch,
    gru_livecrootc_xfer_to_atm_patch = cf.gru_livecrootc_xfer_to_atm_patch,
    gru_deadcrootc_xfer_to_atm_patch = cf.gru_deadcrootc_xfer_to_atm_patch)

# --- Kernel: gross unrepresented LCC fluxes into soil pools (one thread per column) ---
@kernel function _csu2g_col_kernel!(@Const(mask_soilc), decomp_cpools_vr_col,
        @Const(gru_c_to_litr_c_col), @Const(gru_c_to_cwdc_col),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, i_cwd::Int,
        use_soil_matrixcn::Bool, dt)
    c = @index(Global)
    @inbounds if mask_soilc[c]
        if !use_soil_matrixcn
            for j in 1:nlevdecomp
                for i in i_litr_min:i_litr_max
                    decomp_cpools_vr_col[c, j, i] =
                        decomp_cpools_vr_col[c, j, i] +
                        gru_c_to_litr_c_col[c, j, i] * dt
                end
                decomp_cpools_vr_col[c, j, i_cwd] =
                    decomp_cpools_vr_col[c, j, i_cwd] +
                    gru_c_to_cwdc_col[c, j] * dt

                # wood to product pools - states updated in CNProducts
            end
        else
            # Matrix path: matrix_Cinput setup omitted (requires full matrix infrastructure)
        end
    end
end

# --- Kernel: gross unrepresented LCC veg C state updates (one thread per patch) ---
@kernel function _csu2g_patch_kernel!(@Const(mask_soilp), cs, cf, use_matrixcn::Bool, dt)
    p = @index(Global)
    @inbounds if mask_soilp[p]
        # xsmrpool
        cs.xsmrpool_patch[p] = cs.xsmrpool_patch[p] -
            cf.gru_xsmrpool_to_atm_patch[p] * dt
        # gresp storage pool
        cs.gresp_storage_patch[p] = cs.gresp_storage_patch[p] -
            cf.gru_gresp_storage_to_atm_patch[p] * dt
        # gresp transfer pool
        cs.gresp_xfer_patch[p] = cs.gresp_xfer_patch[p] -
            cf.gru_gresp_xfer_to_atm_patch[p] * dt

        if !use_matrixcn
            # displayed pools
            cs.leafc_patch[p] = cs.leafc_patch[p] -
                cf.gru_leafc_to_litter_patch[p] * dt
            cs.frootc_patch[p] = cs.frootc_patch[p] -
                cf.gru_frootc_to_litter_patch[p] * dt
            cs.livestemc_patch[p] = cs.livestemc_patch[p] -
                cf.gru_livestemc_to_atm_patch[p] * dt
            cs.deadstemc_patch[p] = cs.deadstemc_patch[p] -
                cf.gru_deadstemc_to_atm_patch[p] * dt
            cs.deadstemc_patch[p] = cs.deadstemc_patch[p] -
                cf.gru_wood_productc_gain_patch[p] * dt
            cs.livecrootc_patch[p] = cs.livecrootc_patch[p] -
                cf.gru_livecrootc_to_litter_patch[p] * dt
            cs.deadcrootc_patch[p] = cs.deadcrootc_patch[p] -
                cf.gru_deadcrootc_to_litter_patch[p] * dt

            # storage pools
            cs.leafc_storage_patch[p] = cs.leafc_storage_patch[p] -
                cf.gru_leafc_storage_to_atm_patch[p] * dt
            cs.frootc_storage_patch[p] = cs.frootc_storage_patch[p] -
                cf.gru_frootc_storage_to_atm_patch[p] * dt
            cs.livestemc_storage_patch[p] = cs.livestemc_storage_patch[p] -
                cf.gru_livestemc_storage_to_atm_patch[p] * dt
            cs.deadstemc_storage_patch[p] = cs.deadstemc_storage_patch[p] -
                cf.gru_deadstemc_storage_to_atm_patch[p] * dt
            cs.livecrootc_storage_patch[p] = cs.livecrootc_storage_patch[p] -
                cf.gru_livecrootc_storage_to_atm_patch[p] * dt
            cs.deadcrootc_storage_patch[p] = cs.deadcrootc_storage_patch[p] -
                cf.gru_deadcrootc_storage_to_atm_patch[p] * dt

            # transfer pools
            cs.leafc_xfer_patch[p] = cs.leafc_xfer_patch[p] -
                cf.gru_leafc_xfer_to_atm_patch[p] * dt
            cs.frootc_xfer_patch[p] = cs.frootc_xfer_patch[p] -
                cf.gru_frootc_xfer_to_atm_patch[p] * dt
            cs.livestemc_xfer_patch[p] = cs.livestemc_xfer_patch[p] -
                cf.gru_livestemc_xfer_to_atm_patch[p] * dt
            cs.deadstemc_xfer_patch[p] = cs.deadstemc_xfer_patch[p] -
                cf.gru_deadstemc_xfer_to_atm_patch[p] * dt
            cs.livecrootc_xfer_patch[p] = cs.livecrootc_xfer_patch[p] -
                cf.gru_livecrootc_xfer_to_atm_patch[p] * dt
            cs.deadcrootc_xfer_patch[p] = cs.deadcrootc_xfer_patch[p] -
                cf.gru_deadcrootc_xfer_to_atm_patch[p] * dt
        else
            # Matrix version: handled in dynGrossUnrepMod::CNGrossUnrep
        end
    end
end

"""
    c_state_update2g!(cs_veg, cf_veg, cs_soil;
        mask_soilc, mask_soilp, bounds_col, bounds_patch,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd,
        use_matrixcn, use_soil_matrixcn, dt)

Update all prognostic carbon state variables affected by gross
unrepresented landcover change mortality fluxes.

Ported from `CStateUpdate2g` in `CNCStateUpdate2Mod.F90`. Runs as two
KernelAbstractions kernels (one per soil column, one per patch); the CPU
path is byte-identical.

Note: The matrix-CN code paths are omitted (require full matrix infrastructure).
Wood to product pools state updates happen in CNProducts (not here).
"""
function c_state_update2g!(cs_veg::CNVegCarbonStateData,
                            cf_veg::CNVegCarbonFluxData,
                            cs_soil::SoilBiogeochemCarbonStateData;
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
    _launch!(_csu2g_col_kernel!, mask_soilc, cs_soil.decomp_cpools_vr_col,
        cf_veg.gru_c_to_litr_c_col, cf_veg.gru_c_to_cwdc_col,
        nlevdecomp, i_litr_min, i_litr_max, i_cwd, use_soil_matrixcn, dt)

    # --- Patch loop: vegetation C state updates from gross unrepresented landcover change ---
    cs = _csu2g_cs(cs_veg)
    cf = _csu2g_cf(cf_veg)
    _launch!(_csu2g_patch_kernel!, mask_soilp, cs, cf, use_matrixcn, dt)

    return nothing
end
