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
# Standalone: no file I/O, not wired into the driver.
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
"""
Base.@kwdef mutable struct PatchStateUpdater
    pwtgcell_old         ::Vector{Float64} = Float64[]  # old patch weights on the gridcell
    pwtgcell_new         ::Vector{Float64} = Float64[]  # new patch weights on the gridcell
    cwtgcell_old         ::Vector{Float64} = Float64[]  # old column weights on the gridcell
    dwt                  ::Vector{Float64} = Float64[]  # (pwtgcell_new - pwtgcell_old)
    growing_old_fraction ::Vector{Float64} = Float64[]  # (pwtgcell_old / pwtgcell_new); growing only
    growing_new_fraction ::Vector{Float64} = Float64[]  # (dwt / pwtgcell_new); growing only
end

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

"""
    set_old_weights!(updater, bounds, pch, col)

Snapshot subgrid weights before dyn subgrid updates: each patch's `wtgcell`
and each column's `wtgcell`, over the processor bounds.

Ported from `set_old_weights` in dynPatchStateUpdaterMod.F90.
"""
function set_old_weights!(updater::PatchStateUpdater, bounds::BoundsType,
                          pch::PatchData, col::ColumnData)
    for p in bounds.begp:bounds.endp
        c = pch.column[p]
        updater.pwtgcell_old[p] = pch.wtgcell[p]
        updater.cwtgcell_old[c] = col.wtgcell[c]
    end
    return nothing
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
    for p in bounds.begp:bounds.endp
        updater.pwtgcell_new[p] = pch.wtgcell[p]
        updater.dwt[p] = updater.pwtgcell_new[p] - updater.pwtgcell_old[p]
        if updater.dwt[p] > 0.0
            updater.growing_old_fraction[p] = updater.pwtgcell_old[p] / updater.pwtgcell_new[p]
            updater.growing_new_fraction[p] = updater.dwt[p] / updater.pwtgcell_new[p]
        else
            # These values are unused in this case, but set them to something
            # reasonable for safety.
            updater.growing_old_fraction[p] = 1.0
            updater.growing_new_fraction[p] = 0.0
        end
    end
    return nothing
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

    for p in filterp_with_inactive
        c = pch.column[p]

        if updater.dwt[p] > 0.0
            var[p] = var[p] * updater.growing_old_fraction[p]
            if seed !== nothing
                var[p] = var[p] + seed[p] * updater.growing_new_fraction[p]
                if seed_addition !== nothing
                    seed_addition[p] = seed_addition[p] + seed[p] * updater.dwt[p]
                end
            end

        elseif updater.dwt[p] < 0.0
            if flux_out_grc_area !== nothing
                flux_out_grc_area[p] = flux_out_grc_area[p] + var[p] * updater.dwt[p]
            end
            if flux_out_col_area !== nothing
                # No need to check for divide by 0 here: if dwt < 0 then we must
                # have had cwtgcell_old > 0.
                flux_out_col_area[p] = flux_out_col_area[p] +
                    var[p] * (updater.dwt[p] / updater.cwtgcell_old[c])
            end
        end
    end

    return nothing
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

    total_flux_out = zeros(Float64, bounds.endp)
    update_patch_state!(updater, bounds, pch, filterp_with_inactive, var;
                        flux_out_grc_area = total_flux_out,
                        seed = seed, seed_addition = seed_addition)

    for p in filterp_with_inactive
        # Fortran patch%itype is 0-based (0 = bare ground); shift by +1 for the
        # 1-based Julia array.
        my_flux1_fraction = flux1_fraction_by_pft_type[pch.itype[p] + 1]
        flux1_out[p] = flux1_out[p] + total_flux_out[p] * my_flux1_fraction
        flux2_out[p] = flux2_out[p] + total_flux_out[p] * (1.0 - my_flux1_fraction)
    end

    return nothing
end

"""
    old_weight_was_zero(updater, bounds) -> BitVector

Returns a patch-level `BitVector` (indexed `1:endp`) that is true wherever the
patch weight was zero prior to the weight updates.

Ported from `old_weight_was_zero` in dynPatchStateUpdaterMod.F90.
"""
function old_weight_was_zero(updater::PatchStateUpdater, bounds::BoundsType)
    result = falses(bounds.endp)
    for p in bounds.begp:bounds.endp
        result[p] = (updater.pwtgcell_old[p] == 0.0)
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
    for p in bounds.begp:bounds.endp
        result[p] = (updater.dwt[p] > 0.0)
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
    for p in bounds.begp:bounds.endp
        result[p] = (updater.pwtgcell_old[p] == 0.0 && updater.pwtgcell_new[p] > 0.0)
    end
    return result
end
