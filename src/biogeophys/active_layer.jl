# ==========================================================================
# Ported from: src/biogeophys/ActiveLayerMod.F90
# Active layer depth calculation (permafrost diagnostics)
#
# Calculates the depth of the active layer (ALT = active layer thickness)
# and the depth to a frozen interface, using soil temperature profiles.
# Also tracks annual maxima and prior-year maxima for permafrost diagnostics.
# ==========================================================================

"""
    ActiveLayerData

Active layer state data structure. Holds current and annual maximum active
layer thickness (depth of thaw) at column level.

Ported from `active_layer_type` in `ActiveLayerMod.F90`.
"""
Base.@kwdef mutable struct ActiveLayerData
    # --- Public data members ---
    altmax_col               ::Vector{Float64} = Float64[]   # col maximum annual depth of thaw (m)
    altmax_lastyear_col      ::Vector{Float64} = Float64[]   # col prior year maximum annual depth of thaw (m)
    altmax_indx_col          ::Vector{Int}     = Int[]       # col maximum annual depth of thaw index
    altmax_lastyear_indx_col ::Vector{Int}     = Int[]       # col prior year maximum annual depth of thaw index

    # --- Private data members ---
    alt_col                  ::Vector{Float64} = Float64[]   # col current depth of thaw (m)
    alt_indx_col             ::Vector{Int}     = Int[]       # col current depth of thaw index
end

"""
    active_layer_init!(al::ActiveLayerData, ncols::Int)

Allocate and initialize all fields of an `ActiveLayerData` instance for
`ncols` columns. Real fields are initialized to `SPVAL`, integer fields
to `typemax(Int)`.

Ported from `active_layer_type%InitAllocate` in `ActiveLayerMod.F90`.
"""
function active_layer_init!(al::ActiveLayerData, ncols::Int)
    al.alt_col                  = fill(SPVAL, ncols)
    al.altmax_col               = fill(SPVAL, ncols)
    al.altmax_lastyear_col      = fill(SPVAL, ncols)
    al.alt_indx_col             = fill(typemax(Int), ncols)
    al.altmax_indx_col          = fill(typemax(Int), ncols)
    al.altmax_lastyear_indx_col = fill(typemax(Int), ncols)
    return nothing
end

"""
    active_layer_init_cold!(al::ActiveLayerData, col_data::ColumnData,
                            lun::LandunitData, bounds_col::UnitRange{Int})

Cold-start initialization for active layer data. Sets all active layer
fields to zero for soil and crop columns.

Ported from `active_layer_type%InitCold` in `ActiveLayerMod.F90`.
"""
function active_layer_init_cold!(al::ActiveLayerData,
                                  col_data::ColumnData,
                                  lun::LandunitData,
                                  bounds_col::UnitRange{Int})
    for c in bounds_col
        l = col_data.landunit[c]
        if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
            al.alt_col[c]                  = 0.0
            al.altmax_col[c]               = 0.0
            al.altmax_lastyear_col[c]      = 0.0
            al.alt_indx_col[c]             = 0
            al.altmax_indx_col[c]          = 0
            al.altmax_lastyear_indx_col[c] = 0
        end
    end
    return nothing
end

"""
    active_layer_clean!(al::ActiveLayerData)

Deallocate (reset to empty) all fields of an `ActiveLayerData` instance.

Ported from `active_layer_type%Clean`.
"""
function active_layer_clean!(al::ActiveLayerData)
    al.alt_col                  = Float64[]
    al.altmax_col               = Float64[]
    al.altmax_lastyear_col      = Float64[]
    al.alt_indx_col             = Int[]
    al.altmax_indx_col          = Int[]
    al.altmax_lastyear_indx_col = Int[]
    return nothing
end

"""
    alt_calc!(al::ActiveLayerData, mask_soil::BitVector,
              temperature::TemperatureData, col_data::ColumnData,
              grc::GridcellData;
              mon::Int, day::Int, sec::Int, dtime::Int)

Calculate active layer thickness (ALT) for soil columns.

Defines the active layer as the deepest thawed layer, starting from the
base of the soil profile and searching upward for the first layer above
freezing. Where the soil temperature crosses the freezing point between
two layers, linear interpolation is used to estimate the exact depth.

Also updates the annual maximum active layer thickness. The annual maximum
is reset on January 1 for Northern Hemisphere columns and July 1 for
Southern Hemisphere columns.

# Arguments
- `al`: ActiveLayerData instance (modified in place)
- `mask_soil`: BitVector mask for soil columns
- `temperature`: TemperatureData with soil/snow temperature profiles
- `col_data`: ColumnData with gridcell indices
- `grc`: GridcellData with latitude information
- `mon`: current month (1-12)
- `day`: current day of month (1-31)
- `sec`: seconds into current date
- `dtime`: time step length in seconds

Ported from `active_layer_type%alt_calc` in `ActiveLayerMod.F90`.
"""
function alt_calc!(al::ActiveLayerData,
                   mask_soil::BitVector,
                   temperature::TemperatureData,
                   col_data::ColumnData,
                   grc::GridcellData;
                   mon::Int,
                   day::Int,
                   sec::Int,
                   dtime::Int)

    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno

    # Offset for converting Fortran soil index j (1:nlevgrnd) to Julia
    # t_soisno_col index (j + joff). The Julia t_soisno_col matrix has
    # dimensions (ncols, nlevsno + nlevmaxurbgrnd), where the first
    # nlevsno columns correspond to snow layers.
    joff = nlevsno

    # Soil node depths (1:nlevgrnd), populated by varcon_init!()
    zsoi_vals = zsoi[]

    # Aliases for active layer arrays
    alt                  = al.alt_col
    altmax               = al.altmax_col
    altmax_lastyear      = al.altmax_lastyear_col
    alt_indx             = al.alt_indx_col
    altmax_indx          = al.altmax_indx_col
    altmax_lastyear_indx = al.altmax_lastyear_indx_col

    # Soil temperature (column, layer)
    t_soisno = temperature.t_soisno_col

    # ------------------------------------------------------------------
    # On a set annual timestep, update annual maxima.
    # January 1 for NH columns, July 1 for SH columns.
    # ------------------------------------------------------------------
    if mon == 1 && day == 1 && div(sec, dtime) == 1
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            g = col_data.gridcell[c]
            if grc.lat[g] > 0.0
                altmax_lastyear[c]      = altmax[c]
                altmax_lastyear_indx[c] = altmax_indx[c]
                altmax[c]               = 0.0
                altmax_indx[c]          = 0
            end
        end
    end

    if mon == 7 && day == 1 && div(sec, dtime) == 1
        for c in eachindex(mask_soil)
            mask_soil[c] || continue
            g = col_data.gridcell[c]
            if grc.lat[g] <= 0.0
                altmax_lastyear[c]      = altmax[c]
                altmax_lastyear_indx[c] = altmax_indx[c]
                altmax[c]               = 0.0
                altmax_indx[c]          = 0
            end
        end
    end

    # ------------------------------------------------------------------
    # Calculate ALT for each soil column at this timestep.
    # Start from the base of the soil and search upward for the first
    # thawed layer. If a thawed layer is found, linearly interpolate
    # between that layer and the one below to find the exact crossing
    # of the freezing point.
    # ------------------------------------------------------------------
    for c in eachindex(mask_soil)
        mask_soil[c] || continue

        if t_soisno[c, nlevgrnd + joff] > TFRZ
            # Deepest soil layer is above freezing — entire column thawed
            alt[c]      = zsoi_vals[nlevgrnd]
            alt_indx[c] = nlevgrnd
        else
            # Search upward from the second-deepest layer
            k_frz = 0
            found_thawlayer = false
            for j in (nlevgrnd - 1):-1:1
                if t_soisno[c, j + joff] > TFRZ && !found_thawlayer
                    k_frz = j
                    found_thawlayer = true
                end
            end

            if k_frz > 0
                # Linearly interpolate to find the depth where
                # temperature crosses TFRZ
                z1 = zsoi_vals[k_frz]
                z2 = zsoi_vals[k_frz + 1]
                t1 = t_soisno[c, k_frz + joff]
                t2 = t_soisno[c, k_frz + 1 + joff]
                alt[c]      = z1 + (t1 - TFRZ) * (z2 - z1) / (t1 - t2)
                alt_indx[c] = k_frz
            else
                # No thawed layer found — entire column frozen
                alt[c]      = 0.0
                alt_indx[c] = 0
            end
        end

        # Update annual maximum active layer thickness
        if alt[c] > altmax[c]
            altmax[c]      = alt[c]
            altmax_indx[c] = alt_indx[c]
        end
    end

    return nothing
end
