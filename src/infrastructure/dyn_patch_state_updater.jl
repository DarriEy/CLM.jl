# ==========================================================================
# Ported from:
#   src/dyn_subgrid/dynPatchStateUpdaterMod.F90
#
# Class for adjusting patch-level (aboveground) state variables conservatively
# when patch areas change during transient land use.
#
# In each time step, the object should be set up with:
#   - set_old_weights!  (before dyn subgrid weight updates)
#   - set_new_weights!  (after  dyn subgrid weight updates)
# Then each state variable can be updated with:
#   - update_patch_state!
#
# GPU: the shared conservative updaters (set_old/new_weights!, update_patch_state!,
# update_patch_state_partition_flux_by_type!) run through KernelAbstractions
# kernels on the _launch! path. Every write is own-index per-patch (no scatter /
# reduction), so the KA CPU backend is byte-identical to the original host loops
# and the device backend matches at working precision. The struct is
# Adapt-registered so it can move to the device with the CN state it updates.
# The three BitVector query helpers (old_weight_was_zero / patch_grew /
# patch_initiating) stay host — they feed the (host) dyn_cnbal orchestrator and
# BitVectors are CPU-pinned by convention.
# ==========================================================================

"""
    PatchStateUpdater

Holds the old/new patch weights on the gridcell plus derived per-patch
quantities used to conservatively adjust patch-level state when patch areas
change. Ported from `patch_state_updater_type` in dynPatchStateUpdaterMod.F90.

- `pwtgcell_old`         : old patch weights on the gridcell (one per patch)
- `pwtgcell_new`         : new patch weights on the gridcell (one per patch)
- `cwtgcell_old`         : old column weights on the gridcell (one per column)
- `dwt`                  : `pwtgcell_new - pwtgcell_old`
- `growing_old_fraction` : `pwtgcell_old / pwtgcell_new`; only valid for growing patches
- `growing_new_fraction` : `dwt / pwtgcell_new`; only valid for growing patches

The array field type `V` is loose (any `AbstractVector{<:Real}`) so the struct
can be adapted onto a GPU backend / working precision alongside the CN state.
"""
Base.@kwdef mutable struct PatchStateUpdater{V<:AbstractVector{<:Real}}
    pwtgcell_old         ::V = Float64[]  # old patch weights on the gridcell
    pwtgcell_new         ::V = Float64[]  # new patch weights on the gridcell
    cwtgcell_old         ::V = Float64[]  # old column weights on the gridcell
    dwt                  ::V = Float64[]  # (pwtgcell_new - pwtgcell_old)
    growing_old_fraction ::V = Float64[]  # (pwtgcell_old / pwtgcell_new); growing only
    growing_new_fraction ::V = Float64[]  # (dwt / pwtgcell_new); growing only
end
Adapt.@adapt_structure PatchStateUpdater

"""
    PatchStateUpdater(bounds)

Construct a `PatchStateUpdater` sized to the processor bounds, allocating one
entry per patch (`begp:endp`) and per column (`begc:endc`). Patch/column
quantities are initialized to NaN, matching the Fortran constructor.

Ported from `constructor` in dynPatchStateUpdaterMod.F90.
"""
function PatchStateUpdater(bounds::BoundsType)
    @assert bounds.level == BOUNDS_LEVEL_PROC "PatchStateUpdater: argument must be PROC-level bounds"
    return PatchStateUpdater(
        pwtgcell_old         = fill(NaN, bounds.endp),
        pwtgcell_new         = fill(NaN, bounds.endp),
        cwtgcell_old         = fill(NaN, bounds.endc),
        dwt                  = fill(NaN, bounds.endp),
        growing_old_fraction = fill(NaN, bounds.endp),
        growing_new_fraction = fill(NaN, bounds.endp),
    )
end

# --- set_old_weights! ------------------------------------------------------
@kernel function _psu_set_old_kernel!(pwt_old, cwt_old, @Const(p_wtgcell),
                                      @Const(c_wtgcell), @Const(column),
                                      begp::Int, endp::Int)
    p = @index(Global)
    @inbounds if begp <= p <= endp
        pwt_old[p] = p_wtgcell[p]
        c = column[p]
        # Multiple patches on the same column write the same value (idempotent).
        cwt_old[c] = c_wtgcell[c]
    end
end

"""
    set_old_weights!(updater, bounds, pch, col)

Snapshot subgrid weights before dyn subgrid updates: each patch's `wtgcell`
and each column's `wtgcell`, over the processor bounds.

Ported from `set_old_weights` in dynPatchStateUpdaterMod.F90.
"""
function set_old_weights!(updater::PatchStateUpdater, bounds::BoundsType,
                          pch::PatchData, col::ColumnData)
    _launch!(_psu_set_old_kernel!, updater.pwtgcell_old, updater.cwtgcell_old,
             pch.wtgcell, col.wtgcell, pch.column, bounds.begp, bounds.endp;
             ndrange = bounds.endp)
    return nothing
end

# --- set_new_weights! ------------------------------------------------------
@kernel function _psu_set_new_kernel!(pwt_new, dwt, gof, gnf,
                                      @Const(p_wtgcell), @Const(pwt_old),
                                      begp::Int, endp::Int)
    p = @index(Global)
    T = eltype(pwt_new)
    @inbounds if begp <= p <= endp
        pwt_new[p] = p_wtgcell[p]
        d = pwt_new[p] - pwt_old[p]
        dwt[p] = d
        if d > zero(T)
            gof[p] = pwt_old[p] / pwt_new[p]
            gnf[p] = d / pwt_new[p]
        else
            # Unused for non-growing patches; set to safe values (1, 0).
            gof[p] = one(T)
            gnf[p] = zero(T)
        end
    end
end

"""
    set_new_weights!(updater, bounds, pch)

Snapshot subgrid weights after dyn subgrid updates and compute the derived
per-patch quantities (`dwt`, and for growing patches `growing_old_fraction` /
`growing_new_fraction`). For non-growing patches the growing fractions are set
to safe (unused) values 1 and 0.

Ported from `set_new_weights` in dynPatchStateUpdaterMod.F90.
"""
function set_new_weights!(updater::PatchStateUpdater, bounds::BoundsType,
                          pch::PatchData)
    _launch!(_psu_set_new_kernel!, updater.pwtgcell_new, updater.dwt,
             updater.growing_old_fraction, updater.growing_new_fraction,
             pch.wtgcell, updater.pwtgcell_old, bounds.begp, bounds.endp;
             ndrange = bounds.endp)
    return nothing
end

# --- update_patch_state! ---------------------------------------------------
# Per-filter-patch kernel: every write is own-index (var[p], flux_*[p],
# seed_addition[p]) so no atomics/reductions are needed. Optional flux/seed
# outputs are resolved to Bool presence flags on the host; when an optional is
# absent, `var` is passed as a harmless placeholder that the false-guarded
# branch never touches.
@kernel function _update_patch_state_kernel!(var, flux_col, flux_grc, seed, seed_add,
        @Const(filter), @Const(column), @Const(dwt), @Const(gof), @Const(gnf),
        @Const(cwt_old),
        has_flux_col::Bool, has_flux_grc::Bool, has_seed::Bool, has_seed_add::Bool)
    i = @index(Global)
    T = eltype(var)
    @inbounds begin
        p = filter[i]
        c = column[p]
        if dwt[p] > zero(T)
            var[p] = var[p] * gof[p]
            if has_seed
                var[p] = var[p] + seed[p] * gnf[p]
                if has_seed_add
                    seed_add[p] = seed_add[p] + seed[p] * dwt[p]
                end
            end
        elseif dwt[p] < zero(T)
            if has_flux_grc
                flux_grc[p] = flux_grc[p] + var[p] * dwt[p]
            end
            if has_flux_col
                # No divide-by-0 risk: dwt < 0 implies cwtgcell_old > 0.
                flux_col[p] = flux_col[p] + var[p] * (dwt[p] / cwt_old[c])
            end
        end
    end
end

"""
    update_patch_state!(updater, bounds, pch, filterp_with_inactive, var;
                        flux_out_col_area=nothing, flux_out_grc_area=nothing,
                        seed=nothing, seed_addition=nothing)

Update a patch-level state variable `var` and compute associated fluxes based
on changing patch areas.

For growing patches, `var` is scaled by `growing_old_fraction` (preserving the
per-area amount over the larger area) and optionally augmented by `seed`
applied over the newly-grown fraction. For shrinking patches, the state that
leaves is accumulated (as a NEGATIVE quantity) into the optional flux outputs:

- `flux_out_col_area` : mass per unit area COLUMN, using the OLD column weight
  (appropriate when applying these fluxes to column-level state before the
  column-level state adjustment).
- `flux_out_grc_area` : mass per unit area GRIDCELL.

`filterp_with_inactive` is a list of patch indices (including inactive points,
so that patches that just became inactive are processed). Changes are only made
within this filter.

`seed` (optional) gives a seed amount added to growing patches' state over the
area into which they grow; it is ignored for constant/shrinking patches.
`seed_addition` (optional, requires `seed`) accumulates `seed(p) * dwt(p)`,
expressed as mass per unit area GRIDCELL.

Ported from `update_patch_state` in dynPatchStateUpdaterMod.F90.
"""
function update_patch_state!(updater::PatchStateUpdater, bounds::BoundsType,
                             pch::PatchData,
                             filterp_with_inactive::AbstractVector{<:Integer},
                             var::AbstractVector{<:Real};
                             flux_out_col_area::Union{Nothing,AbstractVector{<:Real}}=nothing,
                             flux_out_grc_area::Union{Nothing,AbstractVector{<:Real}}=nothing,
                             seed::Union{Nothing,AbstractVector{<:Real}}=nothing,
                             seed_addition::Union{Nothing,AbstractVector{<:Real}}=nothing)

    if seed_addition !== nothing && seed === nothing
        error("update_patch_state! ERROR: seed_addition can only be provided if seed is provided")
    end

    isempty(filterp_with_inactive) && return nothing

    has_flux_col = flux_out_col_area !== nothing
    has_flux_grc = flux_out_grc_area !== nothing
    has_seed     = seed !== nothing
    has_seed_add = seed_addition !== nothing

    # Placeholders for absent optionals (never indexed under a false flag).
    fcol = has_flux_col ? flux_out_col_area : var
    fgrc = has_flux_grc ? flux_out_grc_area : var
    sd   = has_seed     ? seed              : var
    sadd = has_seed_add ? seed_addition     : var

    _launch!(_update_patch_state_kernel!, var, fcol, fgrc, sd, sadd,
             filterp_with_inactive, pch.column, updater.dwt,
             updater.growing_old_fraction, updater.growing_new_fraction,
             updater.cwtgcell_old,
             has_flux_col, has_flux_grc, has_seed, has_seed_add;
             ndrange = length(filterp_with_inactive))
    return nothing
end

# --- update_patch_state_partition_flux_by_type! ----------------------------
@kernel function _psu_partition_kernel!(flux1_out, flux2_out, @Const(total_flux_out),
        @Const(filter), @Const(itype), @Const(frac_by_type))
    i = @index(Global)
    T = eltype(flux1_out)
    @inbounds begin
        p = filter[i]
        # Fortran itype is 0-based (0 = bare ground); shift +1 for 1-based Julia.
        f1 = frac_by_type[itype[p] + 1]
        flux1_out[p] = flux1_out[p] + total_flux_out[p] * f1
        flux2_out[p] = flux2_out[p] + total_flux_out[p] * (one(T) - f1)
    end
end

"""
    update_patch_state_partition_flux_by_type!(updater, bounds, pch,
        filterp_with_inactive, flux1_fraction_by_pft_type, var, flux1_out, flux2_out;
        seed=nothing, seed_addition=nothing)

Update a patch-level state variable and compute associated fluxes based on
changing patch areas, with the flux out of shrinking patches partitioned into
two fluxes (`flux1_out`, `flux2_out`) according to PFT type. Both flux outputs
are expressed as mass per unit area GRIDCELL (equivalent to the
`flux_out_grc_area` argument of `update_patch_state!`).

`flux1_fraction_by_pft_type` is indexed by PFT type. The Fortran array is
declared `(0:mxpft)` (0-based, bare ground = 0); here it is a 1-based vector of
length `MXPFT+1`, indexed as `flux1_fraction_by_pft_type[itype(p) + 1]`.

Ported from `update_patch_state_partition_flux_by_type` in
dynPatchStateUpdaterMod.F90.
"""
function update_patch_state_partition_flux_by_type!(updater::PatchStateUpdater,
        bounds::BoundsType, pch::PatchData,
        filterp_with_inactive::AbstractVector{<:Integer},
        flux1_fraction_by_pft_type::AbstractVector{<:Real},
        var::AbstractVector{<:Real},
        flux1_out::AbstractVector{<:Real},
        flux2_out::AbstractVector{<:Real};
        seed::Union{Nothing,AbstractVector{<:Real}}=nothing,
        seed_addition::Union{Nothing,AbstractVector{<:Real}}=nothing)

    @assert length(flux1_fraction_by_pft_type) == MXPFT + 1 "update_patch_state_partition_flux_by_type!: flux1_fraction_by_pft_type must have length MXPFT+1 (Fortran 0:mxpft)"

    # Device-resident scratch (lands on var's backend; zeroed like the original).
    total_flux_out = fill!(similar(var, eltype(var), bounds.endp), zero(eltype(var)))
    update_patch_state!(updater, bounds, pch, filterp_with_inactive, var;
                        flux_out_grc_area = total_flux_out,
                        seed = seed, seed_addition = seed_addition)

    isempty(filterp_with_inactive) && return nothing
    _launch!(_psu_partition_kernel!, flux1_out, flux2_out, total_flux_out,
             filterp_with_inactive, pch.itype, flux1_fraction_by_pft_type;
             ndrange = length(filterp_with_inactive))
    return nothing
end

"""
    old_weight_was_zero(updater, bounds) -> BitVector

Returns a patch-level `BitVector` (indexed `1:endp`) that is true wherever the
patch weight was zero prior to the weight updates.

Host-side (BitVector is CPU-pinned); consumed by the dyn_cnbal orchestrator.
Ported from `old_weight_was_zero` in dynPatchStateUpdaterMod.F90.
"""
function old_weight_was_zero(updater::PatchStateUpdater, bounds::BoundsType)
    result = falses(bounds.endp)
    pwtgcell_old = Array(updater.pwtgcell_old)
    for p in bounds.begp:bounds.endp
        result[p] = (pwtgcell_old[p] == 0.0)
    end
    return result
end

"""
    patch_grew(updater, bounds) -> BitVector

Returns a patch-level `BitVector` (indexed `1:endp`) that is true wherever the
patch grew in this time step (`dwt > 0`).

Ported from `patch_grew` in dynPatchStateUpdaterMod.F90.
"""
function patch_grew(updater::PatchStateUpdater, bounds::BoundsType)
    result = falses(bounds.endp)
    dwt = Array(updater.dwt)
    for p in bounds.begp:bounds.endp
        result[p] = (dwt[p] > 0.0)
    end
    return result
end

"""
    patch_initiating(updater, bounds) -> BitVector

Returns a patch-level `BitVector` (indexed `1:endp`) that is true wherever a
patch is initiating — i.e. growing from zero area to non-zero area — in this
time step.

Ported from `patch_initiating` in dynPatchStateUpdaterMod.F90.
"""
function patch_initiating(updater::PatchStateUpdater, bounds::BoundsType)
    result = falses(bounds.endp)
    pwtgcell_old = Array(updater.pwtgcell_old)
    pwtgcell_new = Array(updater.pwtgcell_new)
    for p in bounds.begp:bounds.endp
        result[p] = (pwtgcell_old[p] == 0.0 && pwtgcell_new[p] > 0.0)
    end
    return result
end
