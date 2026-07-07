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
# GPU: the per-state-variable UPDATE PATH (the conservative redistribution called
# once per state variable) runs through KernelAbstractions kernels on the _launch!
# path. The grid-cell reductions (total mass/area lost, non_conserved_mass) use
# `_scatter_add!` (atomic on GPU; ascending-index on the KA CPU backend → the
# per-gridcell sums are byte-identical to the host loops). set_old_weights! /
# set_new_weights! stay HOST: they are per-step setup, and set_old_weights!'s
# natveg-template search is a per-gridcell sequential scan (reduction-natured,
# host-appropriate) and `any_changes` is a per-clump host gate flag. The struct's
# weight fields are Adapt-movable so the update kernels can run on-device.
#
# Five conservation-strategy variants are provided (see module-level docs in the
# Fortran for the full rationale):
#   1. update_column_state_no_special_handling!        — base conservative redistribution
#   2. update_column_state_fill_special_using_natveg!  — special lunits seed from natveg col
#   3. update_column_state_fill_using_fixed_values!     — per-landunit fixed seed values
#   4. update_column_state_fill_special_using_fixed_value! — all special lunits same fixed seed
#   5. update_column_state_with_optional_fractions!     — intermediate (fractional-area) layer
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

The float weight fields (`V`) and the Int `natveg_template_col` (`VI`) are loose
so the struct can be adapted onto a GPU backend for the update kernels;
`any_changes` stays a host `Vector{Bool}` (per-clump gate flag, read on the host).
"""
Base.@kwdef mutable struct ColumnStateUpdater{V<:AbstractVector{<:Real},
                                              VI<:AbstractVector{<:Integer}}
    cwtgcell_old        ::V           = Float64[]
    cwtgcell_new        ::V           = Float64[]
    area_gained_col     ::V           = Float64[]
    natveg_template_col ::VI          = Int[]
    any_changes         ::Vector{Bool} = Bool[]
end

# Custom adapt rule: move the weight/template arrays to the device but keep
# `any_changes` on the host (it is scalar-indexed by clump on the host).
function Adapt.adapt_structure(to, x::ColumnStateUpdater)
    return ColumnStateUpdater(
        cwtgcell_old        = Adapt.adapt(to, x.cwtgcell_old),
        cwtgcell_new        = Adapt.adapt(to, x.cwtgcell_new),
        area_gained_col     = Adapt.adapt(to, x.area_gained_col),
        natveg_template_col = Adapt.adapt(to, x.natveg_template_col),
        any_changes         = x.any_changes,
    )
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

Host-side per-step setup: the natveg-template search is a per-gridcell sequential
scan (reduction-natured, host-appropriate). Ported from `set_old_weights`.
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

Host-side per-step setup: `any_changes[clump_index]` is a host gate flag read on
the host by every update variant. Ported from `set_new_weights`.
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
# Kernels for the per-state-variable update path
# ========================================================================

# Zero a per-column output over [begc, endc] (e.g. adjustment, intent(out)).
@kernel function _csu_zero_col_kernel!(out, begc::Int, endc::Int)
    c = @index(Global)
    T = eltype(out)
    @inbounds if begc <= c <= endc
        out[c] = zero(T)
    end
end

# Accumulate acc[c] += src[c] * scale over [begc, endc] (used by dyn_cnbal_col!
# to depth-integrate the per-level apparent state change into the column total).
@kernel function _csu_accum_scaled_kernel!(acc, @Const(src), scale, begc::Int, endc::Int)
    c = @index(Global)
    @inbounds if begc <= c <= endc
        acc[c] = acc[c] + src[c] * scale
    end
end

# Device-capable helpers wrapping the two column bookkeeping kernels.
function zero_col!(out, begc::Int, endc::Int)
    _launch!(_csu_zero_col_kernel!, out, begc, endc; ndrange = endc)
    return nothing
end
function accumulate_scaled_col!(acc, src, scale::Real, begc::Int, endc::Int)
    _launch!(_csu_accum_scaled_kernel!, acc, src, eltype(acc)(scale), begc, endc;
             ndrange = endc)
    return nothing
end

# Variant 1 builder: every column is prognostic and seeds its own value.
@kernel function _csu_vals_nospecial_kernel!(vals_input, has_prog, @Const(var),
                                             begc::Int, endc::Int)
    c = @index(Global)
    @inbounds if begc <= c <= endc
        vals_input[c] = var[c]
        has_prog[c]   = true
    end
end

# Variant 2 builder: special landunits are non-prognostic and seed from the
# gridcell's natveg template column; natveg columns seed their own value.
@kernel function _csu_vals_natveg_kernel!(vals_input, has_prog,
        @Const(var), @Const(landunit), @Const(ifspecial), @Const(natveg_template),
        template_none::Int, spval, begc::Int, endc::Int)
    c = @index(Global)
    T = eltype(vals_input)
    @inbounds if begc <= c <= endc
        l = landunit[c]
        if ifspecial[l]
            has_prog[c] = false
            tc = natveg_template[c]
            if tc == template_none
                vals_input[c] = T(spval)   # invalid; only an error if this col shrinks
            else
                vals_input[c] = var[tc]
            end
        else
            has_prog[c]   = true
            vals_input[c] = var[c]
        end
    end
end

# Variant 3 builder: columns whose landunit fixed value == FILLVAL keep their own
# (prognostic) value; otherwise they seed the fixed value (non-prognostic).
@kernel function _csu_vals_fixed_kernel!(vals_input, has_prog,
        @Const(var), @Const(landunit), @Const(lun_itype), @Const(landunit_values),
        fillval, begc::Int, endc::Int)
    c = @index(Global)
    T = eltype(vals_input)
    @inbounds if begc <= c <= endc
        l = landunit[c]
        ltype = lun_itype[l]
        my_fillval = landunit_values[ltype]
        if my_fillval == T(fillval)
            vals_input[c] = var[c]
            has_prog[c]   = true
        else
            vals_input[c] = my_fillval
            has_prog[c]   = false
        end
    end
end

# Core K1: per-column reduction of shrinking-column mass/area into the gridcell
# (atomic scatter). Non-prognostic shrinkers pull virtual mass tracked in
# non_conserved_mass. The Fortran shrinking-without-valid-input error() is
# dropped (GPU can't throw); byte-identical on valid input.
@kernel function _csu_loss_kernel!(total_loss_grc, total_area_lost_grc, non_conserved_mass,
        @Const(area_gained_col), @Const(vals_input), @Const(frac_old),
        @Const(has_prog), @Const(gridcell), begc::Int, endc::Int)
    c = @index(Global)
    T = eltype(total_loss_grc)
    @inbounds if begc <= c <= endc
        if area_gained_col[c] < zero(T)
            g = gridcell[c]
            area_lost = -one(T) * area_gained_col[c]
            _scatter_add!(total_area_lost_grc, g, area_lost)
            area_weighted_loss = area_lost * vals_input[c] * frac_old[c]
            _scatter_add!(total_loss_grc, g, area_weighted_loss)
            if !has_prog[c]
                _scatter_add!(non_conserved_mass, g, -area_weighted_loss)
            end
        end
    end
end

# Core K2: per-gridcell mass-loss-per-unit-area (loss pooled over shrinkers).
@kernel function _csu_finalize_kernel!(gain_per_unit_area_grc, @Const(total_loss_grc),
        @Const(total_area_lost_grc), begg::Int, endg::Int)
    g = @index(Global)
    T = eltype(gain_per_unit_area_grc)
    @inbounds if begg <= g <= endg
        if total_area_lost_grc[g] > zero(T)
            gain_per_unit_area_grc[g] = total_loss_grc[g] / total_area_lost_grc[g]
        else
            gain_per_unit_area_grc[g] = zero(T)
        end
    end
end

# Core K3: distribute pooled gain to growing columns. Prognostic columns absorb
# mass into var[c] (own-index) + record adjustment[c]; non-prognostic growers
# throw incoming mass into non_conserved_mass (atomic scatter).
@kernel function _csu_distribute_kernel!(var, adjustment, non_conserved_mass,
        @Const(area_gained_col), @Const(gain_per_unit_area_grc), @Const(has_prog),
        @Const(cwt_old), @Const(cwt_new), @Const(frac_old), @Const(frac_new),
        @Const(gridcell), has_adjustment::Bool, begc::Int, endc::Int)
    c = @index(Global)
    T = eltype(var)
    @inbounds if begc <= c <= endc
        if area_gained_col[c] > zero(T)
            g = gridcell[c]
            mass_gained = area_gained_col[c] * gain_per_unit_area_grc[g]
            if has_prog[c]
                val_old = var[c]
                if frac_new[c] != zero(T)
                    var[c] = (cwt_old[c] * var[c] * frac_old[c] + mass_gained) /
                             (cwt_new[c] * frac_new[c])
                end
                if has_adjustment
                    adjustment[c] = var[c] * frac_new[c] - val_old * frac_old[c]
                end
            else
                _scatter_add!(non_conserved_mass, g, mass_gained)
            end
        end
    end
end

# ========================================================================
# Public update methods
# ========================================================================

"""
    update_column_state_no_special_handling!(this, bounds, clump_index, col, var;
        fractional_area_old=nothing, fractional_area_new=nothing, adjustment=nothing)

Adjust `var` (indexed `begc:endc`, in place) due to changing subgrid weights with
NO special handling; shrinking-column mass is redistributed conservatively to
growing columns. Asserts mass is conserved to `COLSTATE_CONSERVATION_TOLERANCE`.

Ported from `update_column_state_no_special_handling`.
"""
function update_column_state_no_special_handling!(this::ColumnStateUpdater,
        bounds::BoundsType, clump_index::Int, col::ColumnData,
        var::AbstractVector{<:Real};
        fractional_area_old::Union{Nothing,AbstractVector{<:Real}}=nothing,
        fractional_area_new::Union{Nothing,AbstractVector{<:Real}}=nothing,
        adjustment::Union{Nothing,AbstractVector{<:Real}}=nothing)

    if adjustment !== nothing
        _launch!(_csu_zero_col_kernel!, adjustment, bounds.begc, bounds.endc;
                 ndrange = bounds.endc)
    end

    if this.any_changes[clump_index]
        FT = eltype(var)
        vals_input           = similar(var, FT, bounds.endc)
        has_prognostic_state = similar(var, Bool, bounds.endc)
        non_conserved_mass   = fill!(similar(var, FT, bounds.endg), zero(FT))

        _launch!(_csu_vals_nospecial_kernel!, vals_input, has_prognostic_state, var,
                 bounds.begc, bounds.endc; ndrange = bounds.endc)

        update_column_state_with_optional_fractions!(this, bounds, col,
            vals_input, has_prognostic_state,
            var, non_conserved_mass;
            fractional_area_old=fractional_area_old,
            fractional_area_new=fractional_area_new,
            adjustment=adjustment)

        # No special handling here, so non_conserved_mass must be ~0 (host check).
        ncm = Array(non_conserved_mass)
        for g in bounds.begg:bounds.endg
            if abs(ncm[g]) >= COLSTATE_CONSERVATION_TOLERANCE
                error("update_column_state_no_special_handling!: ERROR: failure to " *
                      "conserve mass when using no special handling; g=$g, " *
                      "non_conserved_mass=$(ncm[g])")
            end
        end
    end

    return nothing
end

"""
    update_column_state_fill_special_using_natveg!(this, bounds, clump_index,
        grc, lun, col, var, non_conserved_mass_grc; kwargs...)

Adjust `var` where any shrinking column on a SPECIAL landunit contributes state
equal to the first natural-veg column on its grid cell (`natveg_template_col`).
Special-column virtual mass is tracked in `non_conserved_mass_grc`.

Ported from `update_column_state_fill_special_using_natveg`.
"""
function update_column_state_fill_special_using_natveg!(this::ColumnStateUpdater,
        bounds::BoundsType, clump_index::Int,
        grc::GridcellData, lun::LandunitData, col::ColumnData,
        var::AbstractVector{<:Real}, non_conserved_mass_grc::AbstractVector{<:Real};
        fractional_area_old::Union{Nothing,AbstractVector{<:Real}}=nothing,
        fractional_area_new::Union{Nothing,AbstractVector{<:Real}}=nothing,
        adjustment::Union{Nothing,AbstractVector{<:Real}}=nothing)

    if adjustment !== nothing
        _launch!(_csu_zero_col_kernel!, adjustment, bounds.begc, bounds.endc;
                 ndrange = bounds.endc)
    end

    if this.any_changes[clump_index]
        FT = eltype(var)
        vals_input           = similar(var, FT, bounds.endc)
        has_prognostic_state = similar(var, Bool, bounds.endc)

        _launch!(_csu_vals_natveg_kernel!, vals_input, has_prognostic_state,
                 var, col.landunit, lun.ifspecial, this.natveg_template_col,
                 TEMPLATE_NONE_FOUND, FT(SPVAL), bounds.begc, bounds.endc;
                 ndrange = bounds.endc)

        update_column_state_with_optional_fractions!(this, bounds, col,
            vals_input, has_prognostic_state,
            var, non_conserved_mass_grc;
            fractional_area_old=fractional_area_old,
            fractional_area_new=fractional_area_new,
            adjustment=adjustment)
    end

    return nothing
end

"""
    update_column_state_fill_using_fixed_values!(this, bounds, clump_index,
        lun, col, var, landunit_values, non_conserved_mass_grc; kwargs...)

Adjust `var` where any column on landunit type `i` has state `landunit_values[i]`
when shrinking; `landunit_values[i] == FILLVAL_USE_EXISTING_VALUE` keeps the
existing (prognostic) value. `landunit_values` must have length `MAX_LUNIT`.

Ported from `update_column_state_fill_using_fixed_values`.
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
        _launch!(_csu_zero_col_kernel!, adjustment, bounds.begc, bounds.endc;
                 ndrange = bounds.endc)
    end

    if this.any_changes[clump_index]
        FT = eltype(var)
        vals_input           = similar(var, FT, bounds.endc)
        has_prognostic_state = similar(var, Bool, bounds.endc)

        _launch!(_csu_vals_fixed_kernel!, vals_input, has_prognostic_state,
                 var, col.landunit, lun.itype, landunit_values,
                 FT(FILLVAL_USE_EXISTING_VALUE), bounds.begc, bounds.endc;
                 ndrange = bounds.endc)

        update_column_state_with_optional_fractions!(this, bounds, col,
            vals_input, has_prognostic_state,
            var, non_conserved_mass_grc;
            fractional_area_old=fractional_area_old,
            fractional_area_new=fractional_area_new,
            adjustment=adjustment)
    end

    return nothing
end

"""
    update_column_state_fill_special_using_fixed_value!(this, bounds, clump_index,
        lun, col, var, special_value, non_conserved_mass_grc; kwargs...)

Convenience wrapper to `update_column_state_fill_using_fixed_values!`: vegetated
(non-special) landunits keep their existing value, and ALL special landunits
contribute the same fixed `special_value`.

Ported from `update_column_state_fill_special_using_fixed_value`.
"""
function update_column_state_fill_special_using_fixed_value!(this::ColumnStateUpdater,
        bounds::BoundsType, clump_index::Int,
        lun::LandunitData, col::ColumnData,
        var::AbstractVector{<:Real}, special_value::Real,
        non_conserved_mass_grc::AbstractVector{<:Real};
        fractional_area_old::Union{Nothing,AbstractVector{<:Real}}=nothing,
        fractional_area_new::Union{Nothing,AbstractVector{<:Real}}=nothing,
        adjustment::Union{Nothing,AbstractVector{<:Real}}=nothing)

    # Host-built per-landunit seed table (MAX_LUNIT small); then reuse variant 3.
    landunit_values = Vector{Float64}(undef, MAX_LUNIT)
    for ltype in 1:MAX_LUNIT
        if landunit_is_special(ltype)
            landunit_values[ltype] = special_value
        else
            landunit_values[ltype] = FILLVAL_USE_EXISTING_VALUE
        end
    end
    # Move the seed table onto var's backend so the variant-3 kernel can read it.
    lv = var isa Array ? landunit_values :
         copyto!(similar(var, eltype(var), MAX_LUNIT), landunit_values)

    update_column_state_fill_using_fixed_values!(this, bounds, clump_index,
        lun, col, var, lv, non_conserved_mass_grc;
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
        vals_input, has_prognostic_state, var, non_conserved_mass; kwargs...)

Intermediate routine: resolve the optional fractional areas (defaulting to 1 for
every column when absent — both must be supplied or neither) then call the core
`update_column_state!`.

Ported from `update_column_state_with_optional_fractions`.
"""
function update_column_state_with_optional_fractions!(this::ColumnStateUpdater,
        bounds::BoundsType, col::ColumnData,
        vals_input::AbstractVector{<:Real},
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

    FT = eltype(var)
    # Default the fractional areas to 1 when absent (both or neither supplied).
    my_frac_old = fractional_area_old !== nothing ? fractional_area_old :
                  fill!(similar(var, FT, bounds.endc), one(FT))
    my_frac_new = fractional_area_new !== nothing ? fractional_area_new :
                  fill!(similar(var, FT, bounds.endc), one(FT))

    update_column_state!(this, bounds, col,
        vals_input, has_prognostic_state,
        my_frac_old, my_frac_new,
        var, non_conserved_mass; adjustment=adjustment)

    return nothing
end

"""
    update_column_state!(this, bounds, col, vals_input, has_prognostic_state,
        fractional_area_old, fractional_area_new, var, non_conserved_mass;
        adjustment=nothing)

Core work: for each grid cell, pool the area-weighted mass lost by shrinking
columns and redistribute it per unit area to growing columns. Non-prognostic
(special) columns cannot accept mass; their contribution/uptake is tracked in
`non_conserved_mass` instead.

The Fortran `has_prognostic_state`/`vals_input == var` consistency error() is a
debug assertion, dropped here (byte-identical on valid input); mass conservation
for the no-special variant is still checked on the host after the call.

Ported from `update_column_state`.
"""
function update_column_state!(this::ColumnStateUpdater,
        bounds::BoundsType, col::ColumnData,
        vals_input::AbstractVector{<:Real},
        has_prognostic_state::AbstractVector{Bool},
        fractional_area_old::AbstractVector{<:Real},
        fractional_area_new::AbstractVector{<:Real},
        var::AbstractVector{<:Real}, non_conserved_mass::AbstractVector{<:Real};
        adjustment::Union{Nothing,AbstractVector{<:Real}}=nothing)

    FT = eltype(var)
    total_loss_grc         = fill!(similar(var, FT, bounds.endg), zero(FT))
    total_area_lost_grc    = fill!(similar(var, FT, bounds.endg), zero(FT))
    gain_per_unit_area_grc = fill!(similar(var, FT, bounds.endg), zero(FT))

    _launch!(_csu_loss_kernel!, total_loss_grc, total_area_lost_grc, non_conserved_mass,
             this.area_gained_col, vals_input, fractional_area_old,
             has_prognostic_state, col.gridcell, bounds.begc, bounds.endc;
             ndrange = bounds.endc)

    _launch!(_csu_finalize_kernel!, gain_per_unit_area_grc, total_loss_grc,
             total_area_lost_grc, bounds.begg, bounds.endg; ndrange = bounds.endg)

    has_adj = adjustment !== nothing
    adj = has_adj ? adjustment : var
    _launch!(_csu_distribute_kernel!, var, adj, non_conserved_mass,
             this.area_gained_col, gain_per_unit_area_grc, has_prognostic_state,
             this.cwtgcell_old, this.cwtgcell_new, fractional_area_old, fractional_area_new,
             col.gridcell, has_adj, bounds.begc, bounds.endc; ndrange = bounds.endc)

    return nothing
end
