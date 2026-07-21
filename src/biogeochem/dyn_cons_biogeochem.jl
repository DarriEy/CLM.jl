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
        mask_soilp::AbstractVector{Bool}, filterp::AbstractVector{<:Integer};
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
        mask_soilp::AbstractVector{Bool}, filterp::AbstractVector{<:Integer};
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

# ===========================================================================
# GPU kernelization of the dyn_cnbal_patch! ORCHESTRATOR loops (the sub-functions
# dynamic_patch_adjustments_carbon!/_nitrogen! and dyn_cnbal_col! already run on
# host+Metal via the conservative helpers). One thread per patch; own-index writes
# are race-free, and the patch->gridcell / patch->column accumulations use the
# atomic `_scatter_add!` so concurrent patches on the same gridcell/column add
# correctly. Every kernel is byte-identical to the host loop on the KA CPU backend.
# Field bundles keep the big re-init launch under Metal's ~31-buffer limit.
# ===========================================================================

# --- Re-init bundle: cnveg_state fields reset for initiating patches (+ the
#     column-indexed annavg_t2m_col read). ---
struct _DCReinitVS{V}
    dormant_flag::V; days_active::V; onset_flag::V; onset_counter::V
    onset_gddflag::V; onset_fdd::V; onset_gdd::V; onset_swi::V
    offset_flag::V; offset_counter::V; offset_fdd::V; offset_swi::V
    lgsf::V; bglfr::V; bgtr::V; annavg_t2m::V; tempavg_t2m::V
    c_allometry::V; n_allometry::V; tempsum_potential_gpp::V; annsum_potential_gpp::V
    tempmax_retransn::V; annmax_retransn::V; downreg::V; annavg_t2m_col::V
end
Adapt.@adapt_structure _DCReinitVS

struct _DCReinitCF{V}
    xsmrpool_recover::V; plant_calloc::V; excess_cflux::V
    prev_leafc_to_litter::V; prev_frootc_to_litter::V; availc::V
    gpp_before_downreg::V; tempsum_npp::V; annsum_npp::V
end
Adapt.@adapt_structure _DCReinitCF

struct _DCReinitNF{V}
    plant_ndemand::V; avail_retransn::V; plant_nalloc::V
end
Adapt.@adapt_structure _DCReinitNF

# 1. dwt + re-initialize initiating patches (one thread per patch).
@kernel function _dcp_reinit_kernel!(dwt, dwt_smoothed, laisun, laisha, vs, cfb, nfb,
        @Const(col_of_p), @Const(lun_of_p), @Const(lun_itype),
        @Const(wtgcell), @Const(prior_pwtgcell), @Const(initiating),
        ISTSOIL_::Int, ISTCROP_::Int)
    p = @index(Global)
    T = eltype(dwt)
    @inbounds begin
        l = lun_of_p[p]
        lt = lun_itype[l]
        if lt == ISTSOIL_ || lt == ISTCROP_
            d = wtgcell[p] - prior_pwtgcell[p]
            dwt[p] = d
            dwt_smoothed[p] = d
            if initiating[p]
                c = col_of_p[p]
                z = zero(T); o = one(T)
                laisun[p] = z; laisha[p] = z
                vs.dormant_flag[p] = o; vs.days_active[p] = z
                vs.onset_flag[p] = z; vs.onset_counter[p] = z; vs.onset_gddflag[p] = z
                vs.onset_fdd[p] = z; vs.onset_gdd[p] = z; vs.onset_swi[p] = z
                vs.offset_flag[p] = z; vs.offset_counter[p] = z; vs.offset_fdd[p] = z
                vs.offset_swi[p] = z; vs.lgsf[p] = z; vs.bglfr[p] = z; vs.bgtr[p] = z
                vs.annavg_t2m[p] = vs.annavg_t2m_col[c]
                vs.tempavg_t2m[p] = z; vs.c_allometry[p] = z; vs.n_allometry[p] = z
                vs.tempsum_potential_gpp[p] = z; vs.annsum_potential_gpp[p] = z
                vs.tempmax_retransn[p] = z; vs.annmax_retransn[p] = z; vs.downreg[p] = z
                cfb.xsmrpool_recover[p] = z; cfb.plant_calloc[p] = z; cfb.excess_cflux[p] = z
                cfb.prev_leafc_to_litter[p] = z; cfb.prev_frootc_to_litter[p] = z
                cfb.availc[p] = z; cfb.gpp_before_downreg[p] = z
                cfb.tempsum_npp[p] = z; cfb.annsum_npp[p] = z
                nfb.plant_ndemand[p] = z; nfb.avail_retransn[p] = z; nfb.plant_nalloc[p] = z
            end
        else
            dwt_smoothed[p] = dwt[p]
        end
    end
end

# 2. Flip 3 loss-flux accumulators to positive (own index, all patches).
@kernel function _dcp_signflip3_kernel!(a, b, cc)
    p = @index(Global)
    @inbounds begin
        a[p] = -a[p]; b[p] = -b[p]; cc[p] = -cc[p]
    end
end

# 3. Seeding fluxes: per-patch value + patch->gridcell accumulation.
@kernel function _dcp_seed_kernel!(sc_leaf_p, sc_leaf_g, sc_ds_p, sc_ds_g,
        sn_leaf_p, sn_leaf_g, sn_ds_p, sn_ds_g,
        @Const(dwt_leafc_seed), @Const(dwt_deadstemc_seed),
        @Const(dwt_leafn_seed), @Const(dwt_deadstemn_seed), @Const(gridcell), dt)
    p = @index(Global)
    @inbounds begin
        g = gridcell[p]
        v1 = dwt_leafc_seed[p] / dt;     sc_leaf_p[p] = v1; _scatter_add!(sc_leaf_g, g, v1)
        v2 = dwt_deadstemc_seed[p] / dt; sc_ds_p[p]   = v2; _scatter_add!(sc_ds_g, g, v2)
        v3 = dwt_leafn_seed[p] / dt;     sn_leaf_p[p] = v3; _scatter_add!(sn_leaf_g, g, v3)
        v4 = dwt_deadstemn_seed[p] / dt; sn_ds_p[p]   = v4; _scatter_add!(sn_ds_g, g, v4)
    end
end

# 4. patch->column litter/CWD fluxes (one thread per patch; atomic into col pools).
@kernel function _dcp_p2c_kernel!(frootc_col, frootn_col, livecrootc_col, livecrootn_col,
        deadcrootc_col, deadcrootn_col,
        @Const(dwt_frootc_to_litter), @Const(dwt_frootn_to_litter),
        @Const(dwt_livecrootc_to_litter), @Const(dwt_livecrootn_to_litter),
        @Const(dwt_deadcrootc_to_litter), @Const(dwt_deadcrootn_to_litter),
        @Const(froot_prof), @Const(croot_prof), @Const(fr_f),
        @Const(col_of_p), @Const(itype_of_p),
        nlevdecomp::Int, i_litr_min::Int, i_litr_max::Int, dt)
    p = @index(Global)
    @inbounds begin
        c = col_of_p[p]
        ji = itype_of_p[p] + 1
        for j in 1:nlevdecomp
            fp = froot_prof[p, j]; cp = croot_prof[p, j]
            for i in i_litr_min:i_litr_max
                frf = fr_f[ji, i]
                _scatter_add!(frootc_col, c, j, i, (dwt_frootc_to_litter[p] * frf) / dt * fp)
                _scatter_add!(frootn_col, c, j, i, (dwt_frootn_to_litter[p] * frf) / dt * fp)
            end
            _scatter_add!(livecrootc_col, c, j, dwt_livecrootc_to_litter[p] / dt * cp)
            _scatter_add!(livecrootn_col, c, j, dwt_livecrootn_to_litter[p] / dt * cp)
            _scatter_add!(deadcrootc_col, c, j, dwt_deadcrootc_to_litter[p] / dt * cp)
            _scatter_add!(deadcrootn_col, c, j, dwt_deadcrootn_to_litter[p] / dt * cp)
        end
    end
end

# 5. Product-pool gain fluxes (own index, positive).
@kernel function _dcp_product_kernel!(wood_c_p, crop_c_p, wood_n_p, crop_n_p,
        @Const(wood_product_cflux), @Const(crop_product_cflux),
        @Const(wood_product_nflux), @Const(crop_product_nflux), dt)
    p = @index(Global)
    @inbounds begin
        wood_c_p[p] = -wood_product_cflux[p] / dt
        crop_c_p[p] = -crop_product_cflux[p] / dt
        wood_n_p[p] = -wood_product_nflux[p] / dt
        crop_n_p[p] = -crop_product_nflux[p] / dt
    end
end

# 6. Conversion fluxes: per-patch + patch->gridcell accumulation.
@kernel function _dcp_conv_kernel!(conv_c_p, conv_c_g, conv_n_p, conv_n_g,
        @Const(conv_cflux), @Const(conv_nflux), @Const(gridcell), dt)
    p = @index(Global)
    @inbounds begin
        g = gridcell[p]
        vc = -conv_cflux[p] / dt; conv_c_p[p] = vc; _scatter_add!(conv_c_g, g, vc)
        vn = -conv_nflux[p] / dt; conv_n_p[p] = vn; _scatter_add!(conv_n_g, g, vn)
    end
end

# 7. Slash fluxes (froot+livecroot+deadcroot to litter/CWD): per-patch + p->grc.
@kernel function _dcp_slash_kernel!(slash_c_p, slash_c_g,
        @Const(dwt_frootc_to_litter), @Const(dwt_livecrootc_to_litter),
        @Const(dwt_deadcrootc_to_litter), @Const(gridcell), dt)
    p = @index(Global)
    @inbounds begin
        g = gridcell[p]
        v = dwt_frootc_to_litter[p] / dt + dwt_livecrootc_to_litter[p] / dt +
            dwt_deadcrootc_to_litter[p] / dt
        slash_c_p[p] = v; _scatter_add!(slash_c_g, g, v)
    end
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
        mask_soilp_with_inactive::AbstractVector{Bool},
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

    # Backend-aware scratch (identity on CPU; device-resident when the CN state is
    # on GPU → the whole orchestrator runs on-device, byte-identical on KA CPU).
    _ref = cs.leafc_patch
    FT = eltype(_ref)
    _zeros(n) = fill!(similar(_ref, FT, n), zero(FT))
    _onbk_bool(a) = _ref isa Array ? a : copyto!(similar(_ref, Bool, length(a)), a)
    _onbk_int(a)  = _ref isa Array ? a : copyto!(similar(_ref, Int, length(a)), a)

    dwt                       = _zeros(np)
    dwt_leafc_seed            = _zeros(np)
    dwt_leafn_seed            = _zeros(np)
    dwt_deadstemc_seed        = _zeros(np)
    dwt_deadstemn_seed        = _zeros(np)
    dwt_frootc_to_litter      = _zeros(np)
    dwt_livecrootc_to_litter  = _zeros(np)
    dwt_deadcrootc_to_litter  = _zeros(np)
    dwt_frootn_to_litter      = _zeros(np)
    dwt_livecrootn_to_litter  = _zeros(np)
    dwt_deadcrootn_to_litter  = _zeros(np)
    conv_cflux                = _zeros(np)
    wood_product_cflux        = _zeros(np)
    crop_product_cflux        = _zeros(np)
    conv_nflux                = _zeros(np)
    wood_product_nflux        = _zeros(np)
    crop_product_nflux        = _zeros(np)

    dt_ft      = FT(dt)
    initiating = _onbk_bool(collect(Bool, patch_initiating(updater, bounds)))
    prior_bk   = _to_backend_like(_ref, FT, prior_pwtgcell)
    fr_f_bk    = _to_backend_like(_ref, FT, pftcon_data.fr_f)

    # -----------------------------------------------------------------------
    # 1. Compute dwt (+ dwt_smoothed pass-through) and re-init initiating patches
    # -----------------------------------------------------------------------
    vs = _DCReinitVS(cnveg_state.dormant_flag_patch, cnveg_state.days_active_patch,
        cnveg_state.onset_flag_patch, cnveg_state.onset_counter_patch,
        cnveg_state.onset_gddflag_patch, cnveg_state.onset_fdd_patch,
        cnveg_state.onset_gdd_patch, cnveg_state.onset_swi_patch,
        cnveg_state.offset_flag_patch, cnveg_state.offset_counter_patch,
        cnveg_state.offset_fdd_patch, cnveg_state.offset_swi_patch,
        cnveg_state.lgsf_patch, cnveg_state.bglfr_patch, cnveg_state.bgtr_patch,
        cnveg_state.annavg_t2m_patch, cnveg_state.tempavg_t2m_patch,
        cnveg_state.c_allometry_patch, cnveg_state.n_allometry_patch,
        cnveg_state.tempsum_potential_gpp_patch, cnveg_state.annsum_potential_gpp_patch,
        cnveg_state.tempmax_retransn_patch, cnveg_state.annmax_retransn_patch,
        cnveg_state.downreg_patch, cnveg_state.annavg_t2m_col)
    cfb = _DCReinitCF(cf.xsmrpool_recover_patch, cf.plant_calloc_patch,
        cf.excess_cflux_patch, cf.prev_leafc_to_litter_patch, cf.prev_frootc_to_litter_patch,
        cf.availc_patch, cf.gpp_before_downreg_patch, cf.tempsum_npp_patch, cf.annsum_npp_patch)
    nfb = _DCReinitNF(nf.plant_ndemand_patch, nf.avail_retransn_patch, nf.plant_nalloc_patch)
    _launch!(_dcp_reinit_kernel!, dwt, cnveg_state.dwt_smoothed_patch,
        canopystate.laisun_patch, canopystate.laisha_patch, vs, cfb, nfb,
        patch.column, patch.landunit, lun.itype, patch.wtgcell, prior_bk, initiating,
        Int(ISTSOIL), Int(ISTCROP))

    # Soil patch filter (indices), including inactive points — built on host then
    # moved to the state backend (control-flow array, O(np), done once).
    mask_host = _ref isa Array ? mask_soilp_with_inactive : Array(mask_soilp_with_inactive)
    filterp = _onbk_int(Int[p for p in begp:endp if mask_host[p]])

    # -----------------------------------------------------------------------
    # 2. Adjust patch C and N pools + accumulate loss fluxes (device-capable subfns)
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

    # Flip loss-flux signs to positive.
    _launch!(_dcp_signflip3_kernel!, dwt_frootc_to_litter, dwt_livecrootc_to_litter,
        dwt_deadcrootc_to_litter)

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

    _launch!(_dcp_signflip3_kernel!, dwt_frootn_to_litter, dwt_livecrootn_to_litter,
        dwt_deadcrootn_to_litter)

    # -----------------------------------------------------------------------
    # 3. Seeding fluxes (patch value + patch->gridcell accumulation)
    # -----------------------------------------------------------------------
    _launch!(_dcp_seed_kernel!, cf.dwt_seedc_to_leaf_patch, cf.dwt_seedc_to_leaf_grc,
        cf.dwt_seedc_to_deadstem_patch, cf.dwt_seedc_to_deadstem_grc,
        nf.dwt_seedn_to_leaf_patch, nf.dwt_seedn_to_leaf_grc,
        nf.dwt_seedn_to_deadstem_patch, nf.dwt_seedn_to_deadstem_grc,
        dwt_leafc_seed, dwt_deadstemc_seed, dwt_leafn_seed, dwt_deadstemn_seed,
        patch.gridcell, dt_ft)

    # -----------------------------------------------------------------------
    # 4. patch-to-column litter/CWD fluxes (atomic scatter into column pools)
    # -----------------------------------------------------------------------
    nlevdecomp = size(soilbiogeochem_state.froot_prof_patch, 2)
    _launch!(_dcp_p2c_kernel!, cf.dwt_frootc_to_litr_c_col, nf.dwt_frootn_to_litr_n_col,
        cf.dwt_livecrootc_to_cwdc_col, nf.dwt_livecrootn_to_cwdn_col,
        cf.dwt_deadcrootc_to_cwdc_col, nf.dwt_deadcrootn_to_cwdn_col,
        dwt_frootc_to_litter, dwt_frootn_to_litter,
        dwt_livecrootc_to_litter, dwt_livecrootn_to_litter,
        dwt_deadcrootc_to_litter, dwt_deadcrootn_to_litter,
        soilbiogeochem_state.froot_prof_patch, soilbiogeochem_state.croot_prof_patch,
        fr_f_bk, patch.column, patch.itype,
        nlevdecomp, i_litr_min, i_litr_max, dt_ft; ndrange = np)

    # -----------------------------------------------------------------------
    # 5. Product-pool gain fluxes (positive values)
    # -----------------------------------------------------------------------
    _launch!(_dcp_product_kernel!, cf.dwt_wood_productc_gain_patch,
        cf.dwt_crop_productc_gain_patch, nf.dwt_wood_productn_gain_patch,
        nf.dwt_crop_productn_gain_patch, wood_product_cflux, crop_product_cflux,
        wood_product_nflux, crop_product_nflux, dt_ft)

    # -----------------------------------------------------------------------
    # 6. Conversion fluxes (patch value + patch->gridcell accumulation)
    # -----------------------------------------------------------------------
    _launch!(_dcp_conv_kernel!, cf.dwt_conv_cflux_patch, cf.dwt_conv_cflux_grc,
        nf.dwt_conv_nflux_patch, nf.dwt_conv_nflux_grc,
        conv_cflux, conv_nflux, patch.gridcell, dt_ft)

    # -----------------------------------------------------------------------
    # 7. Slash fluxes (patch value + patch->gridcell accumulation)
    # -----------------------------------------------------------------------
    _launch!(_dcp_slash_kernel!, cf.dwt_slash_cflux_patch, cf.dwt_slash_cflux_grc,
        dwt_frootc_to_litter, dwt_livecrootc_to_litter, dwt_deadcrootc_to_litter,
        patch.gridcell, dt_ft)

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
    # Fortran DynamicColumnAdjustments loops the vertical index `do j = 1, nlevdecomp`
    # (SoilBiogeochemCarbon/NitrogenStateType.F90) — NOT over the full allocated
    # depth. The vr_col pools are allocated 1:nlevdecomp_full (25 on this grid) but
    # only 1:nlevdecomp (=length(dzsoi_decomp), 20) carry decomposition state, and
    # dzsoi_decomp itself is only nlevdecomp long. Iterating size(_,2) here overran
    # dzsoi_decomp on the real cold-start sizing (the wiring test happened to size
    # its arrays to nlevdecomp, so it never exposed this).
    nlevdecomp = length(dzs)

    # Backend-aware adjustment scratch (device-resident when the soil BGC state
    # is on-device; the per-level depth-integration uses zero_col!/accumulate_
    # scaled_col! so no device array is scalar-indexed on the host).
    FT = eltype(cs.decomp_cpools_vr_col)
    adjustment = fill!(similar(cs.decomp_cpools_vr_col, FT, nc), zero(FT))

    # ----- Carbon -----
    zero_col!(cs.dyn_cbal_adjustments_col, begc, endc)

    ndecomp_pools_c = size(cs.decomp_cpools_vr_col, 3)
    for l in 1:ndecomp_pools_c
        for j in 1:nlevdecomp
            update_column_state_no_special_handling!(updater, bounds, clump_index,
                col, @view(cs.decomp_cpools_vr_col[:, j, l]); adjustment = adjustment)
            accumulate_scaled_col!(cs.dyn_cbal_adjustments_col, adjustment, dzs[j], begc, endc)
        end
    end

    for j in 1:nlevdecomp
        update_column_state_no_special_handling!(updater, bounds, clump_index,
            col, @view(cs.ctrunc_vr_col[:, j]); adjustment = adjustment)
        accumulate_scaled_col!(cs.dyn_cbal_adjustments_col, adjustment, dzs[j], begc, endc)
    end

    # ----- Nitrogen -----
    zero_col!(ns.dyn_nbal_adjustments_col, begc, endc)

    ndecomp_pools_n = size(ns.decomp_npools_vr_col, 3)
    for l in 1:ndecomp_pools_n
        for j in 1:nlevdecomp
            update_column_state_no_special_handling!(updater, bounds, clump_index,
                col, @view(ns.decomp_npools_vr_col[:, j, l]); adjustment = adjustment)
            accumulate_scaled_col!(ns.dyn_nbal_adjustments_col, adjustment, dzs[j], begc, endc)
        end
    end

    for j in 1:nlevdecomp
        update_column_state_no_special_handling!(updater, bounds, clump_index,
            col, @view(ns.ntrunc_vr_col[:, j]); adjustment = adjustment)
        accumulate_scaled_col!(ns.dyn_nbal_adjustments_col, adjustment, dzs[j], begc, endc)

        update_column_state_no_special_handling!(updater, bounds, clump_index,
            col, @view(ns.sminn_vr_col[:, j]); adjustment = adjustment)
        accumulate_scaled_col!(ns.dyn_nbal_adjustments_col, adjustment, dzs[j], begc, endc)
    end

    if use_nitrif_denitrif
        for j in 1:nlevdecomp
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
