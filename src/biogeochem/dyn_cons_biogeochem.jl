# ==========================================================================
# Ported from: src/biogeochem/dynConsBiogeochemMod.F90 (864 lines)
#
# Conserves CARBON and NITROGEN across land-cover-area changes — the
# biogeochem analogue of the already-ported dyn_cons_biogeophys.jl
# (water/energy). It handles `dwt` (delta-weight)-driven seeding of newly
# vegetated area, routing of C/N from shrinking patches to product pools +
# litter, and the conservation flux accounting that closes the C/N budget.
#
# Public API (preserves Fortran names):
#   dyn_cnbal_patch!  — patch-level C/N balance across weight change
#   dyn_cnbal_col!    — column-level C/N balance across column weight change
#
# Translation notes vs. Fortran:
# - The Fortran type-bound methods `DynamicPatchAdjustments` (on the cnveg
#   carbon/nitrogen STATE types) and `DynamicColumnAdjustments` (on the
#   soilbgc C/N STATE types) are NOT present as methods in the Julia port,
#   so their bodies are inlined here as private helpers
#   (`dynamic_patch_adjustments_carbon!` / `_nitrogen!` and the column loops
#   in `dyn_cnbal_col!`). They reuse the already-ported conservative
#   machinery: `update_patch_state!`,
#   `update_patch_state_partition_flux_by_type!` (dyn_patch_state_updater.jl)
#   and `update_column_state_no_special_handling!`
#   (dyn_column_state_updater.jl), plus `compute_seed_amounts!`
#   (veg_compute_seed.jl).
# - Fortran integer filters (filter_soilp_with_inactive) become BitVector
#   masks (mask_soilp_with_inactive), per CLAUDE.md conventions.
# - This port handles the BULK C12 carbon + nitrogen state only — the c13 /
#   c14 isotope tracers (use_c13 / use_c14) and the lake-methane column
#   adjustment (use_lch4) are NOT ported here, matching the rest of CLM.jl
#   (no isotopic / methane tracers wired through the dynamic-landcover path).
#   Crop reproductive pools are handled when `use_crop` is true.
# - The Fortran `dwt_dribbler` (annual smoothing of the patch weight delta)
#   is a pass-through with the default annual dynamics; `dwt_smoothed_patch`
#   is set equal to the raw `dwt` here (same simplification used in
#   dyn_cons_biogeophys.jl for the energy/water dribbler).
# - photosyns NewPatchInit (resetting per-patch photosynthesis history) is
#   not ported (no separate photosyns struct in CLM.jl); the canopy LAI and
#   cnveg-state re-initialization for newly-initiating patches IS done.
# ==========================================================================

"""
    DynConsBiogeochemState

Holds the constant seed amounts (gC/m2) used when seeding newly-vegetated
patch area, plus the land-model time step (s). Standalone — NOT part of
`CLMInstances` or any ForwardDiff-dual-copied struct.

Mirrors the Fortran module parameters `leafc_seed` / `deadstemc_seed` and the
`dt = get_step_size_real()` local.

- `leafc_seed`     : seed amount added to growing patches for leaf C (gC/m2)
- `deadstemc_seed` : seed amount added to growing patches for deadstem C (gC/m2)
"""
Base.@kwdef mutable struct DynConsBiogeochemState
    leafc_seed::Float64     = 1.0
    deadstemc_seed::Float64 = 0.1
end

# ---------------------------------------------------------------------------
# Private: dynamic_patch_adjustments_carbon!
# ---------------------------------------------------------------------------

"""
    dynamic_patch_adjustments_carbon!(cs, updater, bounds, patch, pftcon_data,
        mask_soilp, filterp;
        leafc_seed, deadstemc_seed,
        conv_cflux, wood_product_cflux, crop_product_cflux,
        dwt_frootc_to_litter, dwt_livecrootc_to_litter, dwt_deadcrootc_to_litter,
        dwt_leafc_seed, dwt_deadstemc_seed; use_crop=false, nrepr=NREPR)

Adjust patch-level vegetation CARBON state when patch areas change, computing
the associated conversion / wood-product / crop-product / fine-root-to-litter /
coarse-root-to-litter / seed fluxes. The flux accumulators are accumulated into
(as NEGATIVE quantities for the loss fluxes), matching the Fortran convention.

Inlined from `DynamicPatchAdjustments` in `CNVegCarbonStateType.F90`.
`filterp` is a vector of soil patch indices (including inactive points).
"""
function dynamic_patch_adjustments_carbon!(cs::CNVegCarbonStateData,
        updater::PatchStateUpdater, bounds::BoundsType,
        patch::PatchData, pftcon_data::PftconType,
        mask_soilp::BitVector, filterp::AbstractVector{<:Integer};
        leafc_seed::Real, deadstemc_seed::Real,
        conv_cflux::AbstractVector{<:Real},
        wood_product_cflux::AbstractVector{<:Real},
        crop_product_cflux::AbstractVector{<:Real},
        dwt_frootc_to_litter::AbstractVector{<:Real},
        dwt_livecrootc_to_litter::AbstractVector{<:Real},
        dwt_deadcrootc_to_litter::AbstractVector{<:Real},
        dwt_leafc_seed::AbstractVector{<:Real},
        dwt_deadstemc_seed::AbstractVector{<:Real},
        use_crop::Bool=false, nrepr::Int=NREPR)

    bp = bounds.begp:bounds.endp
    np = bounds.endp

    # Backend-aware scratch: allocate on cs.leafc_patch's backend so the whole
    # adjustment runs on-device when the CN state is device-resident (identity on
    # CPU → byte-identical). owz/pgrew return host BitVectors (the updater helpers
    # Array()-copy); collect + copy them to the working backend.
    _ref = cs.leafc_patch
    FT = eltype(_ref)
    _zeros(n) = fill!(similar(_ref, FT, n), zero(FT))
    _onbk_bool(a) = _ref isa Array ? a : copyto!(similar(_ref, Bool, length(a)), a)
    _onbk_f(a)    = _ref isa Array ? a : copyto!(similar(_ref, FT, length(a)), a)

    owz = old_weight_was_zero(updater, bounds)
    pgrew = patch_grew(updater, bounds)

    # ComputeSeedAmounts uses Bool compute/ignore flags (device-resident on GPU)
    compute_here = _onbk_bool(collect(Bool, pgrew))
    ignore_state = _onbk_bool(collect(Bool, owz))
    mask = _onbk_bool(collect(Bool, mask_soilp))

    seed_leafc          = _zeros(np)
    seed_leafc_storage  = _zeros(np)
    seed_leafc_xfer     = _zeros(np)
    seed_deadstemc      = _zeros(np)

    compute_seed_amounts!(mask, bp, patch, pftcon_data;
        species = CN_SPECIES_C12,
        leafc_seed = leafc_seed, deadstemc_seed = deadstemc_seed,
        leaf_patch = cs.leafc_patch,
        leaf_storage_patch = cs.leafc_storage_patch,
        leaf_xfer_patch = cs.leafc_xfer_patch,
        compute_here_patch = compute_here,
        ignore_current_state_patch = ignore_state,
        seed_leaf_patch = seed_leafc,
        seed_leaf_storage_patch = seed_leafc_storage,
        seed_leaf_xfer_patch = seed_leafc_xfer,
        seed_deadstem_patch = seed_deadstemc)

    # --- Leaf pools (seeded, conv flux) ---
    update_patch_state!(updater, bounds, patch, filterp, cs.leafc_patch;
        flux_out_grc_area = conv_cflux, seed = seed_leafc,
        seed_addition = dwt_leafc_seed)
    update_patch_state!(updater, bounds, patch, filterp, cs.leafc_storage_patch;
        flux_out_grc_area = conv_cflux, seed = seed_leafc_storage,
        seed_addition = dwt_leafc_seed)
    update_patch_state!(updater, bounds, patch, filterp, cs.leafc_xfer_patch;
        flux_out_grc_area = conv_cflux, seed = seed_leafc_xfer,
        seed_addition = dwt_leafc_seed)

    # --- Fine root (to litter, per COLUMN area) ---
    update_patch_state!(updater, bounds, patch, filterp, cs.frootc_patch;
        flux_out_col_area = dwt_frootc_to_litter)
    update_patch_state!(updater, bounds, patch, filterp, cs.frootc_storage_patch;
        flux_out_grc_area = conv_cflux)
    update_patch_state!(updater, bounds, patch, filterp, cs.frootc_xfer_patch;
        flux_out_grc_area = conv_cflux)

    # --- Live stem ---
    update_patch_state!(updater, bounds, patch, filterp, cs.livestemc_patch;
        flux_out_grc_area = conv_cflux)
    update_patch_state!(updater, bounds, patch, filterp, cs.livestemc_storage_patch;
        flux_out_grc_area = conv_cflux)
    update_patch_state!(updater, bounds, patch, filterp, cs.livestemc_xfer_patch;
        flux_out_grc_area = conv_cflux)

    # --- Dead stem: split by PFT type into conv (pconv) + wood product ---
    pconv = _onbk_f(build_flux_fraction_by_type(pftcon_data.pconv))
    update_patch_state_partition_flux_by_type!(updater, bounds, patch, filterp,
        pconv, cs.deadstemc_patch, conv_cflux, wood_product_cflux;
        seed = seed_deadstemc, seed_addition = dwt_deadstemc_seed)

    update_patch_state!(updater, bounds, patch, filterp, cs.deadstemc_storage_patch;
        flux_out_grc_area = conv_cflux)
    update_patch_state!(updater, bounds, patch, filterp, cs.deadstemc_xfer_patch;
        flux_out_grc_area = conv_cflux)

    # --- Live coarse root (to litter, per COLUMN area) ---
    update_patch_state!(updater, bounds, patch, filterp, cs.livecrootc_patch;
        flux_out_col_area = dwt_livecrootc_to_litter)
    update_patch_state!(updater, bounds, patch, filterp, cs.livecrootc_storage_patch;
        flux_out_grc_area = conv_cflux)
    update_patch_state!(updater, bounds, patch, filterp, cs.livecrootc_xfer_patch;
        flux_out_grc_area = conv_cflux)

    # --- Dead coarse root (to litter, per COLUMN area) ---
    update_patch_state!(updater, bounds, patch, filterp, cs.deadcrootc_patch;
        flux_out_col_area = dwt_deadcrootc_to_litter)
    update_patch_state!(updater, bounds, patch, filterp, cs.deadcrootc_storage_patch;
        flux_out_grc_area = conv_cflux)
    update_patch_state!(updater, bounds, patch, filterp, cs.deadcrootc_xfer_patch;
        flux_out_grc_area = conv_cflux)

    # --- Remaining C pools (all to conv flux) ---
    update_patch_state!(updater, bounds, patch, filterp, cs.gresp_storage_patch;
        flux_out_grc_area = conv_cflux)
    update_patch_state!(updater, bounds, patch, filterp, cs.gresp_xfer_patch;
        flux_out_grc_area = conv_cflux)
    update_patch_state!(updater, bounds, patch, filterp, cs.cpool_patch;
        flux_out_grc_area = conv_cflux)
    update_patch_state!(updater, bounds, patch, filterp, cs.xsmrpool_patch;
        flux_out_grc_area = conv_cflux)
    update_patch_state!(updater, bounds, patch, filterp, cs.ctrunc_patch;
        flux_out_grc_area = conv_cflux)

    if use_crop
        for k in 1:nrepr
            update_patch_state!(updater, bounds, patch, filterp,
                @view(cs.reproductivec_patch[:, k]);
                flux_out_grc_area = crop_product_cflux)
        end
        for k in 1:nrepr
            update_patch_state!(updater, bounds, patch, filterp,
                @view(cs.reproductivec_storage_patch[:, k]);
                flux_out_grc_area = conv_cflux)
        end
        for k in 1:nrepr
            update_patch_state!(updater, bounds, patch, filterp,
                @view(cs.reproductivec_xfer_patch[:, k]);
                flux_out_grc_area = conv_cflux)
        end
        # cropseedc_deficit is a negative pool: any unrepaid deficit is sucked
        # out of the atmosphere (so its flux_out goes into conv flux).
        update_patch_state!(updater, bounds, patch, filterp, cs.cropseedc_deficit_patch;
            flux_out_grc_area = conv_cflux)
        update_patch_state!(updater, bounds, patch, filterp, cs.xsmrpool_loss_patch;
            flux_out_grc_area = conv_cflux)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Private: dynamic_patch_adjustments_nitrogen!
# ---------------------------------------------------------------------------

"""
    dynamic_patch_adjustments_nitrogen!(ns, updater, bounds, patch, pftcon_data,
        mask_soilp, filterp;
        leafc_seed, deadstemc_seed,
        conv_nflux, wood_product_nflux, crop_product_nflux,
        dwt_frootn_to_litter, dwt_livecrootn_to_litter, dwt_deadcrootn_to_litter,
        dwt_leafn_seed, dwt_deadstemn_seed; use_crop=false, nrepr=NREPR)

Nitrogen analogue of `dynamic_patch_adjustments_carbon!`.
Inlined from `DynamicPatchAdjustments` in `CNVegNitrogenStateType.F90`.
"""
function dynamic_patch_adjustments_nitrogen!(ns::CNVegNitrogenStateData,
        updater::PatchStateUpdater, bounds::BoundsType,
        patch::PatchData, pftcon_data::PftconType,
        mask_soilp::BitVector, filterp::AbstractVector{<:Integer};
        leafc_seed::Real, deadstemc_seed::Real,
        conv_nflux::AbstractVector{<:Real},
        wood_product_nflux::AbstractVector{<:Real},
        crop_product_nflux::AbstractVector{<:Real},
        dwt_frootn_to_litter::AbstractVector{<:Real},
        dwt_livecrootn_to_litter::AbstractVector{<:Real},
        dwt_deadcrootn_to_litter::AbstractVector{<:Real},
        dwt_leafn_seed::AbstractVector{<:Real},
        dwt_deadstemn_seed::AbstractVector{<:Real},
        use_crop::Bool=false, nrepr::Int=NREPR)

    bp = bounds.begp:bounds.endp
    np = bounds.endp

    # Backend-aware scratch (identity on CPU; device-resident when ns is on GPU).
    _ref = ns.leafn_patch
    FT = eltype(_ref)
    _zeros(n) = fill!(similar(_ref, FT, n), zero(FT))
    _onbk_bool(a) = _ref isa Array ? a : copyto!(similar(_ref, Bool, length(a)), a)
    _onbk_f(a)    = _ref isa Array ? a : copyto!(similar(_ref, FT, length(a)), a)

    owz = old_weight_was_zero(updater, bounds)
    pgrew = patch_grew(updater, bounds)
    compute_here = _onbk_bool(collect(Bool, pgrew))
    ignore_state = _onbk_bool(collect(Bool, owz))
    mask = _onbk_bool(collect(Bool, mask_soilp))

    seed_leafn          = _zeros(np)
    seed_leafn_storage  = _zeros(np)
    seed_leafn_xfer     = _zeros(np)
    seed_deadstemn      = _zeros(np)

    compute_seed_amounts!(mask, bp, patch, pftcon_data;
        species = CN_SPECIES_N,
        leafc_seed = leafc_seed, deadstemc_seed = deadstemc_seed,
        leaf_patch = ns.leafn_patch,
        leaf_storage_patch = ns.leafn_storage_patch,
        leaf_xfer_patch = ns.leafn_xfer_patch,
        compute_here_patch = compute_here,
        ignore_current_state_patch = ignore_state,
        seed_leaf_patch = seed_leafn,
        seed_leaf_storage_patch = seed_leafn_storage,
        seed_leaf_xfer_patch = seed_leafn_xfer,
        seed_deadstem_patch = seed_deadstemn)

    # --- Leaf pools (seeded, conv flux) ---
    update_patch_state!(updater, bounds, patch, filterp, ns.leafn_patch;
        flux_out_grc_area = conv_nflux, seed = seed_leafn,
        seed_addition = dwt_leafn_seed)
    update_patch_state!(updater, bounds, patch, filterp, ns.leafn_storage_patch;
        flux_out_grc_area = conv_nflux, seed = seed_leafn_storage,
        seed_addition = dwt_leafn_seed)
    update_patch_state!(updater, bounds, patch, filterp, ns.leafn_xfer_patch;
        flux_out_grc_area = conv_nflux, seed = seed_leafn_xfer,
        seed_addition = dwt_leafn_seed)

    # --- Fine root (to litter, per COLUMN area) ---
    update_patch_state!(updater, bounds, patch, filterp, ns.frootn_patch;
        flux_out_col_area = dwt_frootn_to_litter)
    update_patch_state!(updater, bounds, patch, filterp, ns.frootn_storage_patch;
        flux_out_grc_area = conv_nflux)
    update_patch_state!(updater, bounds, patch, filterp, ns.frootn_xfer_patch;
        flux_out_grc_area = conv_nflux)

    # --- Live stem ---
    update_patch_state!(updater, bounds, patch, filterp, ns.livestemn_patch;
        flux_out_grc_area = conv_nflux)
    update_patch_state!(updater, bounds, patch, filterp, ns.livestemn_storage_patch;
        flux_out_grc_area = conv_nflux)
    update_patch_state!(updater, bounds, patch, filterp, ns.livestemn_xfer_patch;
        flux_out_grc_area = conv_nflux)

    # --- Dead stem: split by PFT type into conv (pconv) + wood product ---
    pconv = _onbk_f(build_flux_fraction_by_type(pftcon_data.pconv))
    update_patch_state_partition_flux_by_type!(updater, bounds, patch, filterp,
        pconv, ns.deadstemn_patch, conv_nflux, wood_product_nflux;
        seed = seed_deadstemn, seed_addition = dwt_deadstemn_seed)

    update_patch_state!(updater, bounds, patch, filterp, ns.deadstemn_storage_patch;
        flux_out_grc_area = conv_nflux)
    update_patch_state!(updater, bounds, patch, filterp, ns.deadstemn_xfer_patch;
        flux_out_grc_area = conv_nflux)

    # --- Live coarse root (to litter, per COLUMN area) ---
    update_patch_state!(updater, bounds, patch, filterp, ns.livecrootn_patch;
        flux_out_col_area = dwt_livecrootn_to_litter)
    update_patch_state!(updater, bounds, patch, filterp, ns.livecrootn_storage_patch;
        flux_out_grc_area = conv_nflux)
    update_patch_state!(updater, bounds, patch, filterp, ns.livecrootn_xfer_patch;
        flux_out_grc_area = conv_nflux)

    # --- Dead coarse root (to litter, per COLUMN area) ---
    update_patch_state!(updater, bounds, patch, filterp, ns.deadcrootn_patch;
        flux_out_col_area = dwt_deadcrootn_to_litter)
    update_patch_state!(updater, bounds, patch, filterp, ns.deadcrootn_storage_patch;
        flux_out_grc_area = conv_nflux)
    update_patch_state!(updater, bounds, patch, filterp, ns.deadcrootn_xfer_patch;
        flux_out_grc_area = conv_nflux)

    # --- Remaining N pools (all to conv flux) ---
    update_patch_state!(updater, bounds, patch, filterp, ns.retransn_patch;
        flux_out_grc_area = conv_nflux)
    update_patch_state!(updater, bounds, patch, filterp, ns.npool_patch;
        flux_out_grc_area = conv_nflux)
    update_patch_state!(updater, bounds, patch, filterp, ns.ntrunc_patch;
        flux_out_grc_area = conv_nflux)

    if use_crop
        for k in 1:nrepr
            update_patch_state!(updater, bounds, patch, filterp,
                @view(ns.reproductiven_patch[:, k]);
                flux_out_grc_area = crop_product_nflux)
        end
        for k in 1:nrepr
            update_patch_state!(updater, bounds, patch, filterp,
                @view(ns.reproductiven_storage_patch[:, k]);
                flux_out_grc_area = conv_nflux)
        end
        for k in 1:nrepr
            update_patch_state!(updater, bounds, patch, filterp,
                @view(ns.reproductiven_xfer_patch[:, k]);
                flux_out_grc_area = conv_nflux)
        end
        # cropseedn_deficit is a negative pool: deficit comes from atmosphere.
        update_patch_state!(updater, bounds, patch, filterp, ns.cropseedn_deficit_patch;
            flux_out_grc_area = conv_nflux)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Private helper: build the (MXPFT+1)-length, 1-based flux-fraction vector
# expected by update_patch_state_partition_flux_by_type! from a 1-based
# per-PFT pftcon vector (which is itself indexed by Fortran-PFT + 1).
# ---------------------------------------------------------------------------
function build_flux_fraction_by_type(pconv::AbstractVector{<:Real})
    n = MXPFT + 1
    out = zeros(Float64, n)
    for i in 1:min(n, length(pconv))
        out[i] = pconv[i]
    end
    return out
end

# ---------------------------------------------------------------------------
# Public: dyn_cnbal_patch!
# ---------------------------------------------------------------------------

"""
    dyn_cnbal_patch!(dynbal, bounds, mask_soilp_with_inactive, prior_pwtgcell,
        updater, patch, lun, col, pftcon_data,
        canopystate, cnveg_state, cnveg_carbonstate, cnveg_carbonflux,
        cnveg_nitrogenstate, cnveg_nitrogenflux, soilbiogeochem_state;
        dt, i_litr_min, i_litr_max, use_crop=false, nrepr=NREPR)

Modify patch-level C/N state and flux variables to maintain carbon and nitrogen
balance when patch weights change under transient land use.

For each soil/crop patch this:
1. Re-initializes phenology / allometry / flux-accumulator state for patches
   that are *initiating* (growing from zero area) and resets canopy LAI.
2. Calls the inlined carbon and nitrogen `DynamicPatchAdjustments`, which
   conservatively scale veg pools for growing/shrinking patches, seed new
   area, and accumulate the conversion / wood-product / crop-product /
   root-to-litter loss fluxes (as negative quantities).
3. Flips the loss-flux signs to positive and translates them into the
   per-time-step seeding, litter/CWD, product-pool, conversion, and slash
   flux fields on the carbon/nitrogen flux structs (per-patch and aggregated
   per-gridcell).

`dynbal` carries the `leafc_seed` / `deadstemc_seed` constants. `dt` is the
land-model time step (s). `mask_soilp_with_inactive` is a BitVector over
patches (true for soil/crop patches, including inactive). `prior_pwtgcell` is
the per-patch gridcell weight prior to the weight update.

Ported from `dyn_cnbal_patch` in `dynConsBiogeochemMod.F90`. Handles bulk C12
+ N only (no c13/c14 isotope tracers).
"""
function dyn_cnbal_patch!(dynbal::DynConsBiogeochemState,
        bounds::BoundsType,
        mask_soilp_with_inactive::BitVector,
        prior_pwtgcell::AbstractVector{<:Real},
        updater::PatchStateUpdater,
        patch::PatchData, lun::LandunitData, col::ColumnData,
        pftcon_data::PftconType,
        canopystate::CanopyStateData,
        cnveg_state::CNVegStateData,
        cnveg_carbonstate::CNVegCarbonStateData,
        cnveg_carbonflux::CNVegCarbonFluxData,
        cnveg_nitrogenstate::CNVegNitrogenStateData,
        cnveg_nitrogenflux::CNVegNitrogenFluxData,
        soilbiogeochem_state::SoilBiogeochemStateData;
        dt::Real,
        i_litr_min::Int, i_litr_max::Int,
        use_crop::Bool=false, nrepr::Int=NREPR)

    begp = bounds.begp
    endp = bounds.endp
    np = endp

    cs = cnveg_carbonstate
    cf = cnveg_carbonflux
    ns = cnveg_nitrogenstate
    nf = cnveg_nitrogenflux

    # patch-level local flux arrays (Fortran allocatables), all zeroed
    z() = zeros(Float64, np)
    dwt                       = z()
    dwt_leafc_seed            = z()
    dwt_leafn_seed            = z()
    dwt_deadstemc_seed        = z()
    dwt_deadstemn_seed        = z()
    dwt_frootc_to_litter      = z()
    dwt_livecrootc_to_litter  = z()
    dwt_deadcrootc_to_litter  = z()
    dwt_frootn_to_litter      = z()
    dwt_livecrootn_to_litter  = z()
    dwt_deadcrootn_to_litter  = z()
    conv_cflux                = z()
    wood_product_cflux        = z()
    crop_product_cflux        = z()
    conv_nflux                = z()
    wood_product_nflux        = z()
    crop_product_nflux        = z()

    initiating = patch_initiating(updater, bounds)

    # -----------------------------------------------------------------------
    # 1. Compute dwt and re-initialize newly-initiating patches
    # -----------------------------------------------------------------------
    for p in begp:endp
        c = patch.column[p]
        l = patch.landunit[p]
        lt = lun.itype[l]
        if lt == ISTSOIL || lt == ISTCROP
            dwt[p] = patch.wtgcell[p] - prior_pwtgcell[p]

            if initiating[p]
                canopystate.laisun_patch[p] = 0.0
                canopystate.laisha_patch[p] = 0.0

                cnveg_state.dormant_flag_patch[p]          = 1.0
                cnveg_state.days_active_patch[p]           = 0.0
                cnveg_state.onset_flag_patch[p]            = 0.0
                cnveg_state.onset_counter_patch[p]         = 0.0
                cnveg_state.onset_gddflag_patch[p]         = 0.0
                cnveg_state.onset_fdd_patch[p]             = 0.0
                cnveg_state.onset_gdd_patch[p]             = 0.0
                cnveg_state.onset_swi_patch[p]             = 0.0
                cnveg_state.offset_flag_patch[p]           = 0.0
                cnveg_state.offset_counter_patch[p]        = 0.0
                cnveg_state.offset_fdd_patch[p]            = 0.0
                cnveg_state.offset_swi_patch[p]            = 0.0
                cnveg_state.lgsf_patch[p]                  = 0.0
                cnveg_state.bglfr_patch[p]                 = 0.0
                cnveg_state.bgtr_patch[p]                  = 0.0
                cnveg_state.annavg_t2m_patch[p]            = cnveg_state.annavg_t2m_col[c]
                cnveg_state.tempavg_t2m_patch[p]           = 0.0
                cnveg_state.c_allometry_patch[p]           = 0.0
                cnveg_state.n_allometry_patch[p]           = 0.0
                cnveg_state.tempsum_potential_gpp_patch[p] = 0.0
                cnveg_state.annsum_potential_gpp_patch[p]  = 0.0
                cnveg_state.tempmax_retransn_patch[p]      = 0.0
                cnveg_state.annmax_retransn_patch[p]       = 0.0
                cnveg_state.downreg_patch[p]               = 0.0

                cf.xsmrpool_recover_patch[p]      = 0.0
                cf.plant_calloc_patch[p]          = 0.0
                cf.excess_cflux_patch[p]          = 0.0
                cf.prev_leafc_to_litter_patch[p]  = 0.0
                cf.prev_frootc_to_litter_patch[p] = 0.0
                cf.availc_patch[p]                = 0.0
                cf.gpp_before_downreg_patch[p]    = 0.0
                cf.tempsum_npp_patch[p]           = 0.0
                cf.annsum_npp_patch[p]            = 0.0

                nf.plant_ndemand_patch[p]  = 0.0
                nf.avail_retransn_patch[p] = 0.0
                nf.plant_nalloc_patch[p]   = 0.0
            end
        end
    end

    # Annually-smoothed (dribbled) change in weight: pass-through here.
    for p in begp:endp
        cnveg_state.dwt_smoothed_patch[p] = dwt[p]
    end

    # Soil patch filter (indices), including inactive points
    filterp = Int[p for p in begp:endp if mask_soilp_with_inactive[p]]

    # -----------------------------------------------------------------------
    # 2. Adjust patch C and N pools + accumulate loss fluxes
    # -----------------------------------------------------------------------
    dynamic_patch_adjustments_carbon!(cs, updater, bounds, patch, pftcon_data,
        mask_soilp_with_inactive, filterp;
        leafc_seed = dynbal.leafc_seed, deadstemc_seed = dynbal.deadstemc_seed,
        conv_cflux = conv_cflux,
        wood_product_cflux = wood_product_cflux,
        crop_product_cflux = crop_product_cflux,
        dwt_frootc_to_litter = dwt_frootc_to_litter,
        dwt_livecrootc_to_litter = dwt_livecrootc_to_litter,
        dwt_deadcrootc_to_litter = dwt_deadcrootc_to_litter,
        dwt_leafc_seed = dwt_leafc_seed,
        dwt_deadstemc_seed = dwt_deadstemc_seed,
        use_crop = use_crop, nrepr = nrepr)

    # These fluxes are computed as negative quantities, but are expected to be
    # positive, so flip the signs.
    for p in begp:endp
        dwt_frootc_to_litter[p]     = -dwt_frootc_to_litter[p]
        dwt_livecrootc_to_litter[p] = -dwt_livecrootc_to_litter[p]
        dwt_deadcrootc_to_litter[p] = -dwt_deadcrootc_to_litter[p]
    end

    dynamic_patch_adjustments_nitrogen!(ns, updater, bounds, patch, pftcon_data,
        mask_soilp_with_inactive, filterp;
        leafc_seed = dynbal.leafc_seed, deadstemc_seed = dynbal.deadstemc_seed,
        conv_nflux = conv_nflux,
        wood_product_nflux = wood_product_nflux,
        crop_product_nflux = crop_product_nflux,
        dwt_frootn_to_litter = dwt_frootn_to_litter,
        dwt_livecrootn_to_litter = dwt_livecrootn_to_litter,
        dwt_deadcrootn_to_litter = dwt_deadcrootn_to_litter,
        dwt_leafn_seed = dwt_leafn_seed,
        dwt_deadstemn_seed = dwt_deadstemn_seed,
        use_crop = use_crop, nrepr = nrepr)

    for p in begp:endp
        dwt_frootn_to_litter[p]     = -dwt_frootn_to_litter[p]
        dwt_livecrootn_to_litter[p] = -dwt_livecrootn_to_litter[p]
        dwt_deadcrootn_to_litter[p] = -dwt_deadcrootn_to_litter[p]
    end

    # -----------------------------------------------------------------------
    # 3. Column-/gridcell-level seeding fluxes
    # -----------------------------------------------------------------------
    for p in begp:endp
        g = patch.gridcell[p]

        cf.dwt_seedc_to_leaf_patch[p] = dwt_leafc_seed[p] / dt
        cf.dwt_seedc_to_leaf_grc[g]   = cf.dwt_seedc_to_leaf_grc[g] +
                                        cf.dwt_seedc_to_leaf_patch[p]

        cf.dwt_seedc_to_deadstem_patch[p] = dwt_deadstemc_seed[p] / dt
        cf.dwt_seedc_to_deadstem_grc[g]   = cf.dwt_seedc_to_deadstem_grc[g] +
                                            cf.dwt_seedc_to_deadstem_patch[p]

        nf.dwt_seedn_to_leaf_patch[p] = dwt_leafn_seed[p] / dt
        nf.dwt_seedn_to_leaf_grc[g]   = nf.dwt_seedn_to_leaf_grc[g] +
                                        nf.dwt_seedn_to_leaf_patch[p]

        nf.dwt_seedn_to_deadstem_patch[p] = dwt_deadstemn_seed[p] / dt
        nf.dwt_seedn_to_deadstem_grc[g]   = nf.dwt_seedn_to_deadstem_grc[g] +
                                            nf.dwt_seedn_to_deadstem_patch[p]
    end

    # -----------------------------------------------------------------------
    # 4. patch-to-column fluxes into litter and CWD pools
    # -----------------------------------------------------------------------
    nlevdecomp = size(soilbiogeochem_state.froot_prof_patch, 2)
    for j in 1:nlevdecomp
        for p in begp:endp
            c = patch.column[p]
            ji = patch.itype[p] + 1  # 1-based pftcon index

            # fine root litter C fluxes
            for i in i_litr_min:i_litr_max
                cf.dwt_frootc_to_litr_c_col[c, j, i] =
                    cf.dwt_frootc_to_litr_c_col[c, j, i] +
                    (dwt_frootc_to_litter[p] * pftcon_data.fr_f[ji, i]) / dt *
                    soilbiogeochem_state.froot_prof_patch[p, j]
            end

            # fine root litter N fluxes
            for i in i_litr_min:i_litr_max
                nf.dwt_frootn_to_litr_n_col[c, j, i] =
                    nf.dwt_frootn_to_litr_n_col[c, j, i] +
                    (dwt_frootn_to_litter[p] * pftcon_data.fr_f[ji, i]) / dt *
                    soilbiogeochem_state.froot_prof_patch[p, j]
            end

            # livecroot -> CWD
            cf.dwt_livecrootc_to_cwdc_col[c, j] =
                cf.dwt_livecrootc_to_cwdc_col[c, j] +
                dwt_livecrootc_to_litter[p] / dt *
                soilbiogeochem_state.croot_prof_patch[p, j]
            nf.dwt_livecrootn_to_cwdn_col[c, j] =
                nf.dwt_livecrootn_to_cwdn_col[c, j] +
                dwt_livecrootn_to_litter[p] / dt *
                soilbiogeochem_state.croot_prof_patch[p, j]

            # deadcroot -> CWD
            cf.dwt_deadcrootc_to_cwdc_col[c, j] =
                cf.dwt_deadcrootc_to_cwdc_col[c, j] +
                dwt_deadcrootc_to_litter[p] / dt *
                soilbiogeochem_state.croot_prof_patch[p, j]
            nf.dwt_deadcrootn_to_cwdn_col[c, j] =
                nf.dwt_deadcrootn_to_cwdn_col[c, j] +
                dwt_deadcrootn_to_litter[p] / dt *
                soilbiogeochem_state.croot_prof_patch[p, j]
        end
    end

    # -----------------------------------------------------------------------
    # 5. Store fluxes into product pools (positive values)
    # -----------------------------------------------------------------------
    for p in begp:endp
        cf.dwt_wood_productc_gain_patch[p] = -wood_product_cflux[p] / dt
        cf.dwt_crop_productc_gain_patch[p] = -crop_product_cflux[p] / dt
        nf.dwt_wood_productn_gain_patch[p] = -wood_product_nflux[p] / dt
        nf.dwt_crop_productn_gain_patch[p] = -crop_product_nflux[p] / dt
    end

    # -----------------------------------------------------------------------
    # 6. Column-level conversion fluxes
    # -----------------------------------------------------------------------
    for p in begp:endp
        g = patch.gridcell[p]

        cf.dwt_conv_cflux_patch[p] = -conv_cflux[p] / dt
        cf.dwt_conv_cflux_grc[g]   = cf.dwt_conv_cflux_grc[g] +
                                     cf.dwt_conv_cflux_patch[p]

        nf.dwt_conv_nflux_patch[p] = -conv_nflux[p] / dt
        nf.dwt_conv_nflux_grc[g]   = nf.dwt_conv_nflux_grc[g] +
                                     nf.dwt_conv_nflux_patch[p]
    end

    # -----------------------------------------------------------------------
    # 7. patch-to-gridcell slash fluxes into litter and CWD pools
    # -----------------------------------------------------------------------
    for p in begp:endp
        g = patch.gridcell[p]

        cf.dwt_slash_cflux_patch[p] =
            dwt_frootc_to_litter[p] / dt +
            dwt_livecrootc_to_litter[p] / dt +
            dwt_deadcrootc_to_litter[p] / dt
        cf.dwt_slash_cflux_grc[g] = cf.dwt_slash_cflux_grc[g] +
                                    cf.dwt_slash_cflux_patch[p]
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Public: dyn_cnbal_col!
# ---------------------------------------------------------------------------

"""
    dyn_cnbal_col!(bounds, clump_index, updater, col,
        soilbiogeochem_carbonstate, soilbiogeochem_nitrogenstate;
        use_nitrif_denitrif=false)

Modify column-level soil-biogeochem state variables to maintain carbon and
nitrogen balance when column weights change under transient land use.

Conservatively redistributes the vertically-resolved decomposition C/N pools
(`decomp_cpools_vr_col` / `decomp_npools_vr_col`), the C/N truncation sinks
(`ctrunc_vr_col` / `ntrunc_vr_col`), and the soil mineral N pool
(`sminn_vr_col`) using `update_column_state_no_special_handling!`, and tracks
the depth-integrated apparent state change in `dyn_cbal_adjustments_col` /
`dyn_nbal_adjustments_col` (weighted by `dzsoi_decomp`). When
`use_nitrif_denitrif` is true, the NO3 / NH4 pools are also redistributed and
tracked in `dyn_no3bal_adjustments_col` / `dyn_nh4bal_adjustments_col`.

Ported from `dyn_cnbal_col` in `dynConsBiogeochemMod.F90`. Handles bulk C12 +
N only (no c13/c14 isotopes, no methane/use_lch4 column adjustment).
"""
function dyn_cnbal_col!(bounds::BoundsType, clump_index::Int,
        updater::ColumnStateUpdater, col::ColumnData,
        soilbiogeochem_carbonstate::SoilBiogeochemCarbonStateData,
        soilbiogeochem_nitrogenstate::SoilBiogeochemNitrogenStateData;
        use_nitrif_denitrif::Bool=false)

    begc = bounds.begc
    endc = bounds.endc
    nc = endc

    cs = soilbiogeochem_carbonstate
    ns = soilbiogeochem_nitrogenstate

    dzs = dzsoi_decomp[]

    # Backend-aware adjustment scratch (device-resident when the soil BGC state
    # is on-device; the per-level depth-integration uses zero_col!/accumulate_
    # scaled_col! so no device array is scalar-indexed on the host).
    FT = eltype(cs.decomp_cpools_vr_col)
    adjustment = fill!(similar(cs.decomp_cpools_vr_col, FT, nc), zero(FT))

    # ----- Carbon -----
    zero_col!(cs.dyn_cbal_adjustments_col, begc, endc)

    ndecomp_pools_c = size(cs.decomp_cpools_vr_col, 3)
    nlevdecomp_c    = size(cs.decomp_cpools_vr_col, 2)
    for l in 1:ndecomp_pools_c
        for j in 1:nlevdecomp_c
            update_column_state_no_special_handling!(updater, bounds, clump_index,
                col, @view(cs.decomp_cpools_vr_col[:, j, l]); adjustment = adjustment)
            accumulate_scaled_col!(cs.dyn_cbal_adjustments_col, adjustment, dzs[j], begc, endc)
        end
    end

    nlevc_trunc = size(cs.ctrunc_vr_col, 2)
    for j in 1:nlevc_trunc
        update_column_state_no_special_handling!(updater, bounds, clump_index,
            col, @view(cs.ctrunc_vr_col[:, j]); adjustment = adjustment)
        accumulate_scaled_col!(cs.dyn_cbal_adjustments_col, adjustment, dzs[j], begc, endc)
    end

    # ----- Nitrogen -----
    zero_col!(ns.dyn_nbal_adjustments_col, begc, endc)

    ndecomp_pools_n = size(ns.decomp_npools_vr_col, 3)
    nlevdecomp_n    = size(ns.decomp_npools_vr_col, 2)
    for l in 1:ndecomp_pools_n
        for j in 1:nlevdecomp_n
            update_column_state_no_special_handling!(updater, bounds, clump_index,
                col, @view(ns.decomp_npools_vr_col[:, j, l]); adjustment = adjustment)
            accumulate_scaled_col!(ns.dyn_nbal_adjustments_col, adjustment, dzs[j], begc, endc)
        end
    end

    nlevn_trunc = size(ns.ntrunc_vr_col, 2)
    for j in 1:nlevn_trunc
        update_column_state_no_special_handling!(updater, bounds, clump_index,
            col, @view(ns.ntrunc_vr_col[:, j]); adjustment = adjustment)
        accumulate_scaled_col!(ns.dyn_nbal_adjustments_col, adjustment, dzs[j], begc, endc)

        update_column_state_no_special_handling!(updater, bounds, clump_index,
            col, @view(ns.sminn_vr_col[:, j]); adjustment = adjustment)
        accumulate_scaled_col!(ns.dyn_nbal_adjustments_col, adjustment, dzs[j], begc, endc)
    end

    if use_nitrif_denitrif
        nlevsmin = size(ns.smin_no3_vr_col, 2)
        for j in 1:nlevsmin
            # These pools aren't included in the overall N balance (totn), so
            # track them separately (reset each level, matching Fortran).
            zero_col!(ns.dyn_no3bal_adjustments_col, begc, endc)
            zero_col!(ns.dyn_nh4bal_adjustments_col, begc, endc)

            update_column_state_no_special_handling!(updater, bounds, clump_index,
                col, @view(ns.smin_no3_vr_col[:, j]); adjustment = adjustment)
            accumulate_scaled_col!(ns.dyn_no3bal_adjustments_col, adjustment, dzs[j], begc, endc)

            update_column_state_no_special_handling!(updater, bounds, clump_index,
                col, @view(ns.smin_nh4_vr_col[:, j]); adjustment = adjustment)
            accumulate_scaled_col!(ns.dyn_nh4bal_adjustments_col, adjustment, dzs[j], begc, endc)
        end
    end

    return nothing
end
