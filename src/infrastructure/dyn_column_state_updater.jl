# ==========================================================================
# Ported from:
#   src/dyn_subgrid/dynColumnStateUpdaterMod.F90
#
# Class for adjusting COLUMN-level state variables conservatively when transient
# subgrid column areas change (columns shrinking / growing within a grid cell).
#
# Usage each time step:
#   - call `set_old_weights!(this, bounds, grc, lun, col)`  (BEFORE weight updates)
#   - call `set_new_weights!(this, bounds, clump_index, col)` (AFTER weight updates)
#   - then update each state var with one of the `update_column_state_*` variants.
#
# All calls must be made from within a loop over clumps; `any_changes` is indexed
# by clump for thread-safety (matching the Fortran design).
#
# Five conservation-strategy variants are provided (see module-level docs in the
# Fortran for the full rationale):
#   1. update_column_state_no_special_handling!        — base conservative redistribution
#   2. update_column_state_fill_special_using_natveg!  — special lunits seed from natveg col
#   3. update_column_state_fill_using_fixed_values!     — per-landunit fixed seed values
#   4. update_column_state_fill_special_using_fixed_value! — all special lunits same fixed seed
#   5. update_column_state_with_optional_fractions!     — intermediate (fractional-area) layer
#
# Standalone: no file I/O, not wired into the driver.
# ==========================================================================

# For update_column_state_fill_using_fixed_values!, any landunit whose
# landunit_values(i) == FILLVAL_USE_EXISTING_VALUE will use the existing value in
# the state variable. (Matches the Fortran FILLVAL_USE_EXISTING_VALUE = spval.)
const FILLVAL_USE_EXISTING_VALUE = SPVAL

# conservation tolerance for the no-special-handling variant
const COLSTATE_CONSERVATION_TOLERANCE = 1.0e-12

"""
    ColumnStateUpdater

Holds the per-column old/new gridcell weights, the derived per-column area change
(`area_gained_col = cwtgcell_new - cwtgcell_old`), the natveg "template" column
for each column, and a per-clump `any_changes` flag.

Ported from `column_state_updater_type` in `dynColumnStateUpdaterMod.F90`.

- `cwtgcell_old`        : old column weights on the grid cell
- `cwtgcell_new`        : new column weights on the grid cell
- `area_gained_col`     : `cwtgcell_new - cwtgcell_old` from last `set_new_weights!`
- `natveg_template_col` : for each column, the first active natveg column in its grid cell
                          (active determined at `set_old_weights!` time → prior-step active)
- `any_changes`         : per-clump flag, true if any column changed area this step
"""
Base.@kwdef mutable struct ColumnStateUpdater
    cwtgcell_old        ::Vector{Float64} = Float64[]
    cwtgcell_new        ::Vector{Float64} = Float64[]
    area_gained_col     ::Vector{Float64} = Float64[]
    natveg_template_col ::Vector{Int}     = Int[]
    any_changes         ::Vector{Bool}    = Bool[]
end

"""
    ColumnStateUpdater(bounds, nclumps)

Construct a `ColumnStateUpdater` sized to the processor bounds (`begc:endc`),
with one `any_changes` slot per clump.

`cwtgcell_old`/`cwtgcell_new`/`area_gained_col` start as `NaN` (uninitialized,
matching the Fortran `nan` fill), `natveg_template_col` starts as
`TEMPLATE_NONE_FOUND`, and `any_changes` starts `false`.

Ported from the `constructor` in `dynColumnStateUpdaterMod.F90`.
"""
function ColumnStateUpdater(bounds::BoundsType, nclumps::Int)
    @assert bounds.level == BOUNDS_LEVEL_PROC "ColumnStateUpdater: argument must be PROC-level bounds"
    n = bounds.endc
    return ColumnStateUpdater(
        cwtgcell_old        = fill(NaN, n),
        cwtgcell_new        = fill(NaN, n),
        area_gained_col     = fill(NaN, n),
        natveg_template_col = fill(TEMPLATE_NONE_FOUND, n),
        any_changes         = fill(false, nclumps),
    )
end

"""
    set_old_weights!(this, bounds, grc, lun, col)

Set subgrid column weights BEFORE the dyn-subgrid weight updates, and recompute
the natveg template column for each column (using the current — i.e. prior-step
— `col.active` flags).

Ported from `set_old_weights` in `dynColumnStateUpdaterMod.F90`.
"""
function set_old_weights!(this::ColumnStateUpdater, bounds::BoundsType,
                          grc::GridcellData, lun::LandunitData, col::ColumnData)
    for c in bounds.begc:bounds.endc
        this.cwtgcell_old[c] = col.wtgcell[c]
    end

    template_col_from_natveg_array!(bounds, col.active, this.natveg_template_col,
                                    grc, lun, col)

    return nothing
end

"""
    set_new_weights!(this, bounds, clump_index, col)

Set subgrid column weights AFTER the dyn-subgrid weight updates, compute
`area_gained_col`, and record in `any_changes[clump_index]` whether any column
changed area this step.

Ported from `set_new_weights` in `dynColumnStateUpdaterMod.F90`.
"""
function set_new_weights!(this::ColumnStateUpdater, bounds::BoundsType,
                          clump_index::Int, col::ColumnData)
    this.any_changes[clump_index] = false
    for c in bounds.begc:bounds.endc
        this.cwtgcell_new[c]    = col.wtgcell[c]
        this.area_gained_col[c] = this.cwtgcell_new[c] - this.cwtgcell_old[c]
        if this.area_gained_col[c] != 0.0
            this.any_changes[clump_index] = true
        end
    end

    return nothing
end

# ========================================================================
# Public update methods
# ========================================================================

"""
    update_column_state_no_special_handling!(this, bounds, clump_index, col, var;
                                              fractional_area_old=nothing,
                                              fractional_area_new=nothing,
                                              adjustment=nothing)

Adjust the values of the column-level state variable `var` (indexed `begc:endc`,
updated in-place) due to changes in subgrid weights, with NO special handling of
any columns. Appropriate for state variables with valid values on all landunits;
shrinking-column mass is redistributed conservatively to growing columns.

Optional `fractional_area_old`/`fractional_area_new` give the fraction of each
column over which the state variable applies (both or neither). Optional
`adjustment` (indexed `begc:endc`) receives the apparent per-column state change.

Asserts that mass is conserved to within `COLSTATE_CONSERVATION_TOLERANCE`.

Ported from `update_column_state_no_special_handling` in `dynColumnStateUpdaterMod.F90`.
"""
function update_column_state_no_special_handling!(this::ColumnStateUpdater,
        bounds::BoundsType, clump_index::Int, col::ColumnData,
        var::AbstractVector{<:Real};
        fractional_area_old::Union{Nothing,AbstractVector{<:Real}}=nothing,
        fractional_area_new::Union{Nothing,AbstractVector{<:Real}}=nothing,
        adjustment::Union{Nothing,AbstractVector{<:Real}}=nothing)

    # Even if there's no work to be done, zero out adjustment (intent(out)).
    if adjustment !== nothing
        for c in bounds.begc:bounds.endc
            adjustment[c] = 0.0
        end
    end

    if this.any_changes[clump_index]
        nc = bounds.endc
        vals_input            = Vector{Float64}(undef, nc)
        vals_input_valid      = Vector{Bool}(undef, nc)
        has_prognostic_state  = Vector{Bool}(undef, nc)
        non_conserved_mass    = zeros(Float64, bounds.endg)

        for c in bounds.begc:bounds.endc
            vals_input[c]           = var[c]
            vals_input_valid[c]     = true
            has_prognostic_state[c] = true
        end

        update_column_state_with_optional_fractions!(this, bounds, col,
            vals_input, vals_input_valid, has_prognostic_state,
            var, non_conserved_mass;
            fractional_area_old=fractional_area_old,
            fractional_area_new=fractional_area_new,
            adjustment=adjustment)

        # No special handling here, so non_conserved_mass must be ~0 (allow roundoff).
        for g in bounds.begg:bounds.endg
            if abs(non_conserved_mass[g]) >= COLSTATE_CONSERVATION_TOLERANCE
                error("update_column_state_no_special_handling!: ERROR: failure to " *
                      "conserve mass when using no special handling; g=$g, " *
                      "non_conserved_mass=$(non_conserved_mass[g])")
            end
        end
    end

    return nothing
end

"""
    update_column_state_fill_special_using_natveg!(this, bounds, clump_index,
        grc, lun, col, var, non_conserved_mass_grc;
        fractional_area_old=nothing, fractional_area_new=nothing, adjustment=nothing)

Adjust `var` due to changing subgrid weights, where any shrinking column on a
SPECIAL landunit contributes state equal to the first natural-veg column on its
grid cell (its `natveg_template_col`). Special landunits do not have prognostic
state for this variable, so growing special columns throw away incoming mass and
shrinking special columns pull "virtual" mass out of thin air; both are tracked
in `non_conserved_mass_grc` (indexed `begg:endg`, accumulated into).

If a special column has no natveg template, its input value is invalid; this is
only an error if that column is actually shrinking (see `update_column_state!`).

Ported from `update_column_state_fill_special_using_natveg` in `dynColumnStateUpdaterMod.F90`.
"""
function update_column_state_fill_special_using_natveg!(this::ColumnStateUpdater,
        bounds::BoundsType, clump_index::Int,
        grc::GridcellData, lun::LandunitData, col::ColumnData,
        var::AbstractVector{<:Real}, non_conserved_mass_grc::AbstractVector{<:Real};
        fractional_area_old::Union{Nothing,AbstractVector{<:Real}}=nothing,
        fractional_area_new::Union{Nothing,AbstractVector{<:Real}}=nothing,
        adjustment::Union{Nothing,AbstractVector{<:Real}}=nothing)

    if adjustment !== nothing
        for c in bounds.begc:bounds.endc
            adjustment[c] = 0.0
        end
    end

    if this.any_changes[clump_index]
        nc = bounds.endc
        vals_input           = Vector{Float64}(undef, nc)
        vals_input_valid     = Vector{Bool}(undef, nc)
        has_prognostic_state = Vector{Bool}(undef, nc)

        for c in bounds.begc:bounds.endc
            l = col.landunit[c]
            if lun.ifspecial[l]
                has_prognostic_state[c] = false

                template_col = this.natveg_template_col[c]
                if template_col == TEMPLATE_NONE_FOUND
                    vals_input[c]       = SPVAL
                    vals_input_valid[c] = false
                else
                    vals_input[c]       = var[template_col]
                    vals_input_valid[c] = true
                end
            else
                has_prognostic_state[c] = true
                vals_input[c]           = var[c]
                vals_input_valid[c]     = true
            end
        end

        update_column_state_with_optional_fractions!(this, bounds, col,
            vals_input, vals_input_valid, has_prognostic_state,
            var, non_conserved_mass_grc;
            fractional_area_old=fractional_area_old,
            fractional_area_new=fractional_area_new,
            adjustment=adjustment)
    end

    return nothing
end

"""
    update_column_state_fill_using_fixed_values!(this, bounds, clump_index,
        lun, col, var, landunit_values, non_conserved_mass_grc;
        fractional_area_old=nothing, fractional_area_new=nothing, adjustment=nothing)

Adjust `var` due to changing subgrid weights, where any column on landunit type
`i` is assumed to have state equal to `landunit_values[i]` when shrinking. If
`landunit_values[i] == FILLVAL_USE_EXISTING_VALUE`, columns of that type keep
their existing value (and are prognostic — they can accept mass); otherwise the
landunit is specially treated (cannot accept mass; contributes/loses virtual
mass tracked in `non_conserved_mass_grc`, indexed `begg:endg`).

`landunit_values` must have length `MAX_LUNIT`.

Ported from `update_column_state_fill_using_fixed_values` in `dynColumnStateUpdaterMod.F90`.
"""
function update_column_state_fill_using_fixed_values!(this::ColumnStateUpdater,
        bounds::BoundsType, clump_index::Int,
        lun::LandunitData, col::ColumnData,
        var::AbstractVector{<:Real}, landunit_values::AbstractVector{<:Real},
        non_conserved_mass_grc::AbstractVector{<:Real};
        fractional_area_old::Union{Nothing,AbstractVector{<:Real}}=nothing,
        fractional_area_new::Union{Nothing,AbstractVector{<:Real}}=nothing,
        adjustment::Union{Nothing,AbstractVector{<:Real}}=nothing)

    @assert length(landunit_values) == MAX_LUNIT "update_column_state_fill_using_fixed_values!: must provide values for all landunits"

    if adjustment !== nothing
        for c in bounds.begc:bounds.endc
            adjustment[c] = 0.0
        end
    end

    if this.any_changes[clump_index]
        nc = bounds.endc
        vals_input           = Vector{Float64}(undef, nc)
        vals_input_valid     = Vector{Bool}(undef, nc)
        has_prognostic_state = Vector{Bool}(undef, nc)

        for c in bounds.begc:bounds.endc
            l = col.landunit[c]
            ltype = lun.itype[l]
            my_fillval = landunit_values[ltype]

            if my_fillval == FILLVAL_USE_EXISTING_VALUE
                vals_input[c]           = var[c]
                vals_input_valid[c]     = true
                has_prognostic_state[c] = true
            else
                vals_input[c]           = my_fillval
                vals_input_valid[c]     = true
                has_prognostic_state[c] = false
            end
        end

        update_column_state_with_optional_fractions!(this, bounds, col,
            vals_input, vals_input_valid, has_prognostic_state,
            var, non_conserved_mass_grc;
            fractional_area_old=fractional_area_old,
            fractional_area_new=fractional_area_new,
            adjustment=adjustment)
    end

    return nothing
end

"""
    update_column_state_fill_special_using_fixed_value!(this, bounds, clump_index,
        lun, col, var, special_value, non_conserved_mass_grc;
        fractional_area_old=nothing, fractional_area_new=nothing, adjustment=nothing)

Convenience wrapper to `update_column_state_fill_using_fixed_values!`: vegetated
(non-special) landunits keep their existing value (`FILLVAL_USE_EXISTING_VALUE`),
and ALL special landunits contribute the same fixed `special_value`.

Ported from `update_column_state_fill_special_using_fixed_value` in `dynColumnStateUpdaterMod.F90`.
"""
function update_column_state_fill_special_using_fixed_value!(this::ColumnStateUpdater,
        bounds::BoundsType, clump_index::Int,
        lun::LandunitData, col::ColumnData,
        var::AbstractVector{<:Real}, special_value::Real,
        non_conserved_mass_grc::AbstractVector{<:Real};
        fractional_area_old::Union{Nothing,AbstractVector{<:Real}}=nothing,
        fractional_area_new::Union{Nothing,AbstractVector{<:Real}}=nothing,
        adjustment::Union{Nothing,AbstractVector{<:Real}}=nothing)

    landunit_values = Vector{Float64}(undef, MAX_LUNIT)
    for ltype in 1:MAX_LUNIT
        if landunit_is_special(ltype)
            landunit_values[ltype] = special_value
        else
            landunit_values[ltype] = FILLVAL_USE_EXISTING_VALUE
        end
    end

    update_column_state_fill_using_fixed_values!(this, bounds, clump_index,
        lun, col, var, landunit_values, non_conserved_mass_grc;
        fractional_area_old=fractional_area_old,
        fractional_area_new=fractional_area_new,
        adjustment=adjustment)

    return nothing
end

# ========================================================================
# Private methods
# ========================================================================

"""
    update_column_state_with_optional_fractions!(this, bounds, col,
        vals_input, vals_input_valid, has_prognostic_state, var, non_conserved_mass;
        fractional_area_old=nothing, fractional_area_new=nothing, adjustment=nothing)

Intermediate routine between the public update variants and the core work routine
`update_column_state!`. Resolves the optional fractional areas (defaulting to 1
for every column when absent — both must be supplied or neither) then calls the
core routine.

Ported from `update_column_state_with_optional_fractions` in `dynColumnStateUpdaterMod.F90`.
"""
function update_column_state_with_optional_fractions!(this::ColumnStateUpdater,
        bounds::BoundsType, col::ColumnData,
        vals_input::AbstractVector{<:Real}, vals_input_valid::AbstractVector{Bool},
        has_prognostic_state::AbstractVector{Bool},
        var::AbstractVector{<:Real}, non_conserved_mass::AbstractVector{<:Real};
        fractional_area_old::Union{Nothing,AbstractVector{<:Real}}=nothing,
        fractional_area_new::Union{Nothing,AbstractVector{<:Real}}=nothing,
        adjustment::Union{Nothing,AbstractVector{<:Real}}=nothing)

    if (fractional_area_old !== nothing) && (fractional_area_new === nothing)
        error("update_column_state_with_optional_fractions!: ERROR: If fractional_area_old is provided, then fractional_area_new must be provided, too")
    end
    if (fractional_area_new !== nothing) && (fractional_area_old === nothing)
        error("update_column_state_with_optional_fractions!: ERROR: If fractional_area_new is provided, then fractional_area_old must be provided, too")
    end

    my_fractional_area_old = Vector{Float64}(undef, bounds.endc)
    my_fractional_area_new = Vector{Float64}(undef, bounds.endc)

    if fractional_area_old !== nothing
        for c in bounds.begc:bounds.endc
            my_fractional_area_old[c] = fractional_area_old[c]
        end
    else
        for c in bounds.begc:bounds.endc
            my_fractional_area_old[c] = 1.0
        end
    end

    if fractional_area_new !== nothing
        for c in bounds.begc:bounds.endc
            my_fractional_area_new[c] = fractional_area_new[c]
        end
    else
        for c in bounds.begc:bounds.endc
            my_fractional_area_new[c] = 1.0
        end
    end

    update_column_state!(this, bounds, col,
        vals_input, vals_input_valid, has_prognostic_state,
        my_fractional_area_old, my_fractional_area_new,
        var, non_conserved_mass; adjustment=adjustment)

    return nothing
end

"""
    update_column_state!(this, bounds, col,
        vals_input, vals_input_valid, has_prognostic_state,
        fractional_area_old, fractional_area_new, var, non_conserved_mass;
        adjustment=nothing)

Do the work of conservatively redistributing a column-level state variable when
subgrid column weights change. For each grid cell, the area-weighted mass lost by
shrinking columns is pooled and redistributed per unit area to the growing
columns. Columns without prognostic state (special handling) cannot accept mass;
their contribution/uptake is tracked in `non_conserved_mass` instead.

`var` is updated in-place; `non_conserved_mass` (indexed `begg:endg`) is
accumulated into; optional `adjustment` (indexed `begc:endc`) receives the
apparent per-column state change.

Ported from `update_column_state` in `dynColumnStateUpdaterMod.F90`.
"""
function update_column_state!(this::ColumnStateUpdater,
        bounds::BoundsType, col::ColumnData,
        vals_input::AbstractVector{<:Real}, vals_input_valid::AbstractVector{Bool},
        has_prognostic_state::AbstractVector{Bool},
        fractional_area_old::AbstractVector{<:Real},
        fractional_area_new::AbstractVector{<:Real},
        var::AbstractVector{<:Real}, non_conserved_mass::AbstractVector{<:Real};
        adjustment::Union{Nothing,AbstractVector{<:Real}}=nothing)

    # ------------------------------------------------------------------------
    # Error-checking on inputs
    # ------------------------------------------------------------------------
    # For the sake of conservation, vals_input must equal var wherever
    # has_prognostic_state is true (and vals_input is valid).
    for c in bounds.begc:bounds.endc
        if has_prognostic_state[c] && vals_input_valid[c]
            if vals_input[c] != var[c]
                error("update_column_state!: ERROR: where has_prognostic_state is true, " *
                      "vals_input must equal var; c=$c, vals_input=$(vals_input[c]), var=$(var[c])")
            end
        end
    end

    # ------------------------------------------------------------------------
    # Begin main work
    # ------------------------------------------------------------------------
    # Determine the total mass loss for each grid cell, plus the gross area loss
    # (which should match the gross area gain).
    total_loss_grc      = zeros(Float64, bounds.endg)
    total_area_lost_grc = zeros(Float64, bounds.endg)

    for c in bounds.begc:bounds.endc
        g = col.gridcell[c]
        if this.area_gained_col[c] < 0.0
            if !vals_input_valid[c]
                error("update_column_state!: ERROR: shrinking column without valid input value; c=$c")
            end
            area_lost = -1.0 * this.area_gained_col[c]
            total_area_lost_grc[g] += area_lost
            area_weighted_loss = area_lost * vals_input[c] * fractional_area_old[c]
            total_loss_grc[g] += area_weighted_loss

            if !has_prognostic_state[c]
                # This column doesn't model the state variable, so its vals_input is
                # a fictitious quantity. Track how much of it we added to the system.
                non_conserved_mass[g] -= area_weighted_loss
            end
        end
    end

    # Mass loss per unit area for each grid cell: lump all loss into a "loss" pool
    # per grid cell, to then distribute amongst growing columns.
    gain_per_unit_area_grc = zeros(Float64, bounds.endg)
    for g in bounds.begg:bounds.endg
        if total_area_lost_grc[g] > 0.0
            gain_per_unit_area_grc[g] = total_loss_grc[g] / total_area_lost_grc[g]
        else
            gain_per_unit_area_grc[g] = 0.0
        end
    end

    # Distribute gain to growing columns.
    for c in bounds.begc:bounds.endc
        g = col.gridcell[c]
        if this.area_gained_col[c] > 0.0
            mass_gained = this.area_gained_col[c] * gain_per_unit_area_grc[g]
            if has_prognostic_state[c]
                val_old = var[c]

                # Avoid divide-by-zero. fractional_area_new == 0 can only happen if
                # fractional_area_old(c) == 0 and the shrinking columns' fractional
                # areas were all 0 — in which case var is irrelevant for conservation.
                if fractional_area_new[c] != 0.0
                    var[c] = (this.cwtgcell_old[c] * var[c] * fractional_area_old[c] + mass_gained) /
                             (this.cwtgcell_new[c] * fractional_area_new[c])
                end

                if adjustment !== nothing
                    adjustment[c] = var[c] * fractional_area_new[c] -
                                    val_old * fractional_area_old[c]
                end
            else
                non_conserved_mass[g] += mass_gained
            end
        end
    end

    return nothing
end
