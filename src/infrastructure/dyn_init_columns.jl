# ==========================================================================
# Ported from:
#   src/dyn_subgrid/dynColumnTemplateMod.F90
#   src/dyn_subgrid/dynInitColumnsMod.F90
#
# Handle initialization of columns that just switched from inactive to active
# during transient land use (e.g. a crop landunit that transitions from 0 to
# >0 weight on a grid cell). For each newly-active column we SELECT an existing
# "template" (exemplar) column to copy state from, then COPY a defined subset
# of column-level state into the new column.
# ==========================================================================

# If no template column was found, this value is returned.
# (Matches the Fortran TEMPLATE_NONE_FOUND = ispval.)
const TEMPLATE_NONE_FOUND = ISPVAL

# --------------------------------------------------------------------------
# Functions that operate on a single column at a time
# --------------------------------------------------------------------------

"""
    template_col_from_landunit(bounds, c_target, landunit_type, cactive, grc, lun, col) -> c_template

Find a column to serve as a template for the state variables on the target
column `c_target`.

Looks for a landunit of the type given by `landunit_type` (e.g. `ISTSOIL`,
`ISTCROP`) in the same grid cell as `c_target`. Returns the index of the first
*active* column on that landunit (order within a landunit is arbitrary — given
by memory order). If there are no active columns on that landunit in this grid
cell, returns `TEMPLATE_NONE_FOUND`.

`cactive` is the column-level active flags — generally from the *prior* time
step, so that we don't identify a point that just became active this step.

Ported from `template_col_from_landunit` in `dynColumnTemplateMod.F90`.
"""
function template_col_from_landunit(bounds::BoundsType, c_target::Int,
                                    landunit_type::Int,
                                    cactive::AbstractVector{Bool},
                                    grc::GridcellData, lun::LandunitData,
                                    col::ColumnData)
    found = false
    g = col.gridcell[c_target]
    l = grc.landunit_indices[landunit_type, g]

    c = ISPVAL
    # If this landunit exists on this grid cell...
    if l != ISPVAL
        # Loop through columns on this landunit; stop as soon as we find an
        # active column: that will serve as the template.
        c = lun.coli[l]
        while !found && c <= lun.colf[l]
            if cactive[c]
                found = true
            else
                c += 1
            end
        end
    end

    if found
        c_template = c
    else
        c_template = TEMPLATE_NONE_FOUND
    end

    return c_template
end

# --------------------------------------------------------------------------
# Subroutines that operate on the whole column-level array at once
# --------------------------------------------------------------------------

"""
    template_col_from_natveg_array!(bounds, cactive, c_templates, grc, lun, col)

For each column, find a column to serve as a template for the state variables on
the target column by looking for the first active column on the natural-veg
(`ISTSOIL`) landunit in the same grid cell. Assigns `TEMPLATE_NONE_FOUND` where
there is no active natural-veg column in the grid cell.

Note: if there are multiple columns on the natural-veg landunit, then a given
natural-veg column may get a template column that differs from itself! The
caller is responsible for deciding when this template should be used.

Ported from `template_col_from_natveg_array` in `dynColumnTemplateMod.F90`.
"""
function template_col_from_natveg_array!(bounds::BoundsType,
                                         cactive::AbstractVector{Bool},
                                         c_templates::AbstractVector{<:Integer},
                                         grc::GridcellData, lun::LandunitData,
                                         col::ColumnData)
    for c in bounds.begc:bounds.endc
        c_templates[c] = template_col_from_landunit(bounds, c, ISTSOIL,
                                                     cactive, grc, lun, col)
    end

    nothing
end

# --------------------------------------------------------------------------
# Initialization of newly-active columns
# --------------------------------------------------------------------------

"""
    initialize_new_columns!(bounds, cactive_prior, temp, ws, grc, lun, col)

Do initialization for all columns that are newly-active in this time step.

For each column that is active now (`col.active`) but was not active in the
prior time step (`cactive_prior`), pick a template column via
[`initial_template_col_dispatcher`](@ref) and copy the template column's state
into it via [`copy_state!`](@ref). If no template column is found, the state
already in memory is kept (possibly arbitrary initial conditions) and a warning
is emitted.

`cactive_prior` is the column-level active flags from the prior time step (the
PriorWeights snapshot in CLM); here it is taken as a plain `BitVector` /
`Vector{Bool}` to stay decoupled.

Ported from `initialize_new_columns` in `dynInitColumnsMod.F90`.
"""
function initialize_new_columns!(bounds::BoundsType,
                                 cactive_prior::AbstractVector{Bool},
                                 temp::TemperatureData, ws::WaterStateData,
                                 grc::GridcellData, lun::LandunitData,
                                 col::ColumnData)
    subname = "initialize_new_columns!"

    for c in bounds.begc:bounds.endc
        # If this column is newly-active, then we need to initialize it.
        if col.active[c] && !cactive_prior[c]
            c_template = initial_template_col_dispatcher(bounds, c,
                                                         cactive_prior,
                                                         grc, lun, col)
            if c_template != TEMPLATE_NONE_FOUND
                copy_state!(c, c_template, temp, ws)
            else
                @warn "$subname WARNING: No template column found to " *
                      "initialize newly-active column -- keeping the state " *
                      "that was already in memory, possibly from arbitrary " *
                      "initialization (c=$c)"
            end
        end
    end

    nothing
end

"""
    initial_template_col_dispatcher(bounds, c_new, cactive_prior, grc, lun, col) -> c_template

Find a column to use as a template for the given column that has newly become
active; this dispatches to the appropriate routine based on the landunit type
of `c_new`. Returns `TEMPLATE_NONE_FOUND` if there is no column to use for
initialization.

`ISTICE`/`ISTDLAK`/`ISTWET`/urban landunits are not yet supported (these match
the Fortran `endrun` cases) and raise an error.

Ported from `initial_template_col_dispatcher` in `dynInitColumnsMod.F90`.
"""
function initial_template_col_dispatcher(bounds::BoundsType, c_new::Int,
                                         cactive_prior::AbstractVector{Bool},
                                         grc::GridcellData, lun::LandunitData,
                                         col::ColumnData)
    subname = "initial_template_col_dispatcher"

    l = col.landunit[c_new]
    ltype = lun.itype[l]
    if ltype == ISTSOIL
        c_template = initial_template_col_soil(c_new, col)
    elseif ltype == ISTCROP
        c_template = initial_template_col_crop(bounds, c_new, cactive_prior,
                                               grc, lun, col)
    elseif ltype == ISTICE
        error("$subname ERROR: Ability to initialize a newly-active glacier " *
              "mec column not yet implemented (expectation is that glacier " *
              "mec columns should be active from the start of the run " *
              "wherever they can grow) (c_new=$c_new)")
    elseif ltype == ISTDLAK
        error("$subname ERROR: Ability to initialize a newly-active lake " *
              "column not yet implemented (c_new=$c_new)")
    elseif ltype == ISTWET
        error("$subname ERROR: Ability to initialize a newly-active wetland " *
              "column not yet implemented (c_new=$c_new)")
    elseif ISTURB_MIN <= ltype <= ISTURB_MAX
        error("$subname ERROR: Ability to initialize a newly-active urban " *
              "column not yet implemented (c_new=$c_new)")
    else
        error("$subname ERROR: Unknown landunit type: $ltype (c_new=$c_new)")
    end

    return c_template
end

"""
    initial_template_col_soil(c_new, col) -> c_template

Find a column to use as a template for a vegetated (natural-veg / `ISTSOIL`)
column that has newly become active.

For now, the only vegetated columns that can newly become active are ones with
0 weight on the grid cell (i.e. virtual columns). For these, we simply keep the
state at its current value (likely arbitrary initial conditions) and return
`TEMPLATE_NONE_FOUND`. This assumption is checked here.

Ported from `initial_template_col_soil` in `dynInitColumnsMod.F90`.
"""
function initial_template_col_soil(c_new::Int, col::ColumnData)
    subname = "initial_template_col_soil"

    if col.wtgcell[c_new] > 0.0
        error("$subname ERROR: Expectation is that the only vegetated " *
              "columns that can newly become active are ones with 0 weight " *
              "on the grid cell (c_new=$c_new, wtgcell=$(col.wtgcell[c_new]))")
    end

    c_template = TEMPLATE_NONE_FOUND

    return c_template
end

"""
    initial_template_col_crop(bounds, c_new, cactive_prior, grc, lun, col) -> c_template

Find a column to use as a template for a crop (`ISTCROP`) column that has newly
become active.

First tries to find an active column on the vegetated (`ISTSOIL`) landunit; if
there is none, tries the first active column on the crop (`ISTCROP`) landunit;
if there is still none, returns `TEMPLATE_NONE_FOUND`.

Ported from `initial_template_col_crop` in `dynInitColumnsMod.F90`.
"""
function initial_template_col_crop(bounds::BoundsType, c_new::Int,
                                   cactive_prior::AbstractVector{Bool},
                                   grc::GridcellData, lun::LandunitData,
                                   col::ColumnData)
    c_template = template_col_from_landunit(bounds, c_new, ISTSOIL,
                                            cactive_prior, grc, lun, col)
    if c_template == TEMPLATE_NONE_FOUND
        c_template = template_col_from_landunit(bounds, c_new, ISTCROP,
                                                cactive_prior, grc, lun, col)
    end

    return c_template
end

"""
    copy_state!(c_new, c_template, temp, ws)

Copy a subset of column-level state variables from a template column
(`c_template`) to a newly-active column (`c_new`).

Only the *below-ground* portion of the multi-level variables is copied, not the
above-ground (snow) portion: it is challenging to initialize the snow pack in a
consistent state, so the new column's snow pack is left at cold-start
conditions. In the Julia SoA layout the snow layers occupy the leading
`nlevsno` columns of `t_soisno_col` / `h2osoi_liq_col` / `h2osoi_ice_col`, so
the below-ground portion (Fortran index `1:`) begins at Julia column
`nlevsno + 1`. `h2osoi_vol_col` carries no snow padding (`1:nlevmaxurbgrnd`),
so it is copied in full.

State copied (matching the Fortran `copy_state`):
- `temp.t_soisno_col`   (below-ground levels)
- `ws.h2osoi_liq_col`   (below-ground levels)
- `ws.h2osoi_ice_col`   (below-ground levels)
- `ws.h2osoi_vol_col`   (all levels)
- `ws.wa_col`           (scalar)

Ported from `copy_state` in `dynInitColumnsMod.F90`.
"""
function copy_state!(c_new::Int, c_template::Int,
                     temp::TemperatureData, ws::WaterStateData)
    nlevsno = varpar.nlevsno

    # Below-ground portion starts after the leading snow layers (Fortran 1:).
    soi_lb = nlevsno + 1
    soi_ub = size(temp.t_soisno_col, 2)

    @views temp.t_soisno_col[c_new, soi_lb:soi_ub] .=
        temp.t_soisno_col[c_template, soi_lb:soi_ub]

    # h2osoi_liq_col / h2osoi_ice_col share the snow-padded layout.
    liq_ub = size(ws.h2osoi_liq_col, 2)
    @views ws.h2osoi_liq_col[c_new, soi_lb:liq_ub] .=
        ws.h2osoi_liq_col[c_template, soi_lb:liq_ub]
    @views ws.h2osoi_ice_col[c_new, soi_lb:liq_ub] .=
        ws.h2osoi_ice_col[c_template, soi_lb:liq_ub]

    # h2osoi_vol_col has no snow padding (1:nlevmaxurbgrnd) -> copy in full.
    vol_ub = size(ws.h2osoi_vol_col, 2)
    @views ws.h2osoi_vol_col[c_new, 1:vol_ub] .=
        ws.h2osoi_vol_col[c_template, 1:vol_ub]

    ws.wa_col[c_new] = ws.wa_col[c_template]

    nothing
end
