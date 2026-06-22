# ==========================================================================
# Ported from:
#   src/dyn_subgrid/dynLandunitAreaMod.F90  (update_landunit_weights[_one_gcell])
#   src/dyn_subgrid/dynPriorWeightsMod.F90  (prior_weights_type + set_prior_weights)
#
# Subgrid WEIGHT bookkeeping after transient land-use areas change.
# Standalone: no file I/O, not wired into the driver.
# ==========================================================================

# --------------------------------------------------------------------------
# dynLandunitAreaMod
# --------------------------------------------------------------------------

# This parameter specifies the order in which landunit areas are decreased when
# the specified areas add to greater than 100%. Landunits not listed here can
# never be decreased unless the forcings say they should be decreased. In
# particular, note that istice doesn't appear here, so that the istice area
# always matches the areas specified by GLC. In general, the code will NOT be
# robust if more than one landunit is excluded from this list. Meaning: since
# istice is excluded from this list, all other landunits should appear here.
# Matches dynLandunitAreaMod.F90 decrease_order(7).
const DECREASE_ORDER = (ISTSOIL, ISTCROP, ISTURB_MD, ISTURB_HD, ISTURB_TBD,
                        ISTWET, ISTDLAK)

"""
    set_landunit_weight!(grc, lun, g, ltype, weight)

Set weight of landunit type `ltype` on gridcell `g`. No-op if that landunit
type does not exist on the gridcell. Mirror of `subgridWeightsMod`'s
`set_landunit_weight`, used by `update_landunit_weights!`.
"""
function set_landunit_weight!(grc::GridcellData, lun::LandunitData,
                              g::Int, ltype::Int, weight::Real)
    l = grc.landunit_indices[ltype, g]
    if l != ISPVAL
        lun.wtgcell[l] = weight
    end
    nothing
end

"""
    update_landunit_weights_one_gcell!(landunit_weights)

Update landunit weights for a single grid cell.

`landunit_weights` is a vector of length `MAX_LUNIT`, updated in-place. Element
`ltype` is the weight of landunit type `ltype` (e.g. `ISTSOIL`, `ISTCROP`, ...).
Landunits that do not exist in this grid cell should be given a weight of 0.

After execution, `sum(landunit_weights)` is 1 within a small roundoff-level
tolerance, achieved by growing or shrinking landunits as needed:
- if the sum is already ~1, do nothing;
- if the sum is < 1, increase natural vegetation (`ISTSOIL`) to absorb the deficit;
- if the sum is > 1, decrease areas in `DECREASE_ORDER` priority (never below 0).

Faithful to `update_landunit_weights_one_gcell` in dynLandunitAreaMod.F90.
"""
function update_landunit_weights_one_gcell!(landunit_weights::AbstractVector{<:Real})
    tol = 1.0e-14  # tolerance for making sure sum of landunit weights equals 1

    @assert length(landunit_weights) == MAX_LUNIT "update_landunit_weights_one_gcell!: landunit_weights must have length MAX_LUNIT"

    landunit_sum = sum(landunit_weights)

    if abs(landunit_sum - 1.0) <= tol
        # If landunits sum to ~ 100% already, we're done.

    elseif landunit_sum < 1.0
        # If landunits sum to < 100%, increase natural vegetation so sum is 100%.
        landunit_weights[ISTSOIL] = landunit_weights[ISTSOIL] + (1.0 - landunit_sum)

    else
        # If landunits sum to > 100%, decrease areas in priority order.
        decrease_index = 1
        excess = landunit_sum - 1.0
        while (excess > tol) && (decrease_index <= length(DECREASE_ORDER))
            # Decrease weight of the next landunit, but not below 0.
            cur_landunit = DECREASE_ORDER[decrease_index]
            landunit_weights[cur_landunit] = landunit_weights[cur_landunit] - excess
            if landunit_weights[cur_landunit] < 0.0
                landunit_weights[cur_landunit] = 0.0
            end

            # Update variables for next loop iteration.
            landunit_sum = sum(landunit_weights)
            excess = landunit_sum - 1.0
            decrease_index = decrease_index + 1
        end
    end

    # Confirm that landunit sum is now equal to 100%, within tolerance.
    landunit_sum = sum(landunit_weights)
    if abs(landunit_sum - 1.0) > tol
        error("update_landunit_weights_one_gcell! ERROR: After all landunit " *
              "adjustments, landunit weights still do not equal 100%; " *
              "landunit_sum = $landunit_sum, landunit_weights = $landunit_weights")
    end

    nothing
end

"""
    update_landunit_weights!(bounds, grc, lun)

Update landunit weights for all grid cells.

Assumes `lun.wtgcell` has been updated for all landunits whose areas are
specified by the dynamic subgrid code. Updates `lun.wtgcell` for all other
landunits (and may change some specified-area landunits, e.g. when glacier
area and crop area conflict) so that the landunit weights on each gridcell sum
to 1, applying the `DECREASE_ORDER` priority (natural veg / soil absorbs the
residual last).

Faithful to `update_landunit_weights` in dynLandunitAreaMod.F90.
"""
function update_landunit_weights!(bounds::BoundsType, grc::GridcellData,
                                  lun::LandunitData)
    landunit_weights = Vector{Float64}(undef, MAX_LUNIT)

    for g in bounds.begg:bounds.endg
        # Determine current landunit weights. Landunits that don't exist on
        # this grid cell get a weight of 0.
        for ltype in 1:MAX_LUNIT
            landunit_weights[ltype] = get_landunit_weight(grc, lun, g, ltype)
        end

        # Adjust weights so they sum to 100%.
        update_landunit_weights_one_gcell!(landunit_weights)

        # Put the new landunit weights back into lun.wtgcell.
        for ltype in 1:MAX_LUNIT
            set_landunit_weight!(grc, lun, g, ltype, landunit_weights[ltype])
        end
    end

    nothing
end


# --------------------------------------------------------------------------
# dynPriorWeightsMod
# --------------------------------------------------------------------------

"""
    PriorWeights

Holds prior subgrid weights (i.e. before the weight updates of this time step),
used later for state conservation. A passive data-holder: components should be
treated as read-only after `set_prior_weights!`.

Ported from `prior_weights_type` in dynPriorWeightsMod.F90.

- `pwtgcell` : prior patch weight on the gridcell (one entry per patch)
- `cactive`  : prior `col.active` flags (one entry per column)
"""
Base.@kwdef mutable struct PriorWeights
    pwtgcell ::Vector{Float64} = Float64[]  # prior patch weight on the gridcell
    cactive  ::BitVector       = falses(0)  # prior col.active flags
end

"""
    PriorWeights(bounds)

Construct a `PriorWeights` sized to the processor bounds, allocating one entry
per patch (`begp:endp`) and per column (`begc:endc`).

Ported from the `prior_weights_type` constructor in dynPriorWeightsMod.F90.
"""
function PriorWeights(bounds::BoundsType)
    @assert bounds.level == BOUNDS_LEVEL_PROC "PriorWeights: argument must be PROC-level bounds"
    return PriorWeights(
        pwtgcell = zeros(Float64, bounds.endp),
        cactive  = falses(bounds.endc),
    )
end

"""
    set_prior_weights!(prior, bounds, pch, col)

Snapshot current weights into `prior`: each patch's `wtgcell` and each column's
`active` flag, over the processor bounds.

Ported from `set_prior_weights` in dynPriorWeightsMod.F90.
"""
function set_prior_weights!(prior::PriorWeights, bounds::BoundsType,
                            pch::PatchData, col::ColumnData)
    for p in bounds.begp:bounds.endp
        prior.pwtgcell[p] = pch.wtgcell[p]
    end

    for c in bounds.begc:bounds.endc
        prior.cactive[c] = col.active[c]
    end

    nothing
end
