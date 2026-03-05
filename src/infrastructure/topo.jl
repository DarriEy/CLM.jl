# ==========================================================================
# Ported from: src/main/TopoMod.F90
# Topographic height handling for each column
# ==========================================================================

# --------------------------------------------------------------------------
# Main topo data type
# --------------------------------------------------------------------------

"""
    TopoData

Column-level topographic height data. Holds the surface elevation for each
column and a flag indicating whether the column needs atmospheric downscaling.

Ported from `topo_type` in `TopoMod.F90`.
"""
Base.@kwdef mutable struct TopoData{FT<:AbstractFloat}
    topo_col::Vector{FT} = Float64[]                  # surface elevation (m)
    needs_downscaling_col::Vector{Bool} = Bool[]           # whether a column needs to be downscaled
end

# --------------------------------------------------------------------------
# Init (mirrors Fortran Init = InitAllocate + InitHistory + InitCold)
# --------------------------------------------------------------------------

"""
    topo_init!(topo, ncols, col, lun;
               topo_glc_mec=nothing,
               use_hillslope=false,
               downscale_hillslope_meteorology=false)

Initialize topographic data: allocate, register history (no-op), and
cold-start initialize.

Ported from `Init` in `TopoMod.F90`.
"""
function topo_init!(topo::TopoData, ncols::Int,
                    col::ColumnData, lun::LandunitData;
                    topo_glc_mec::Union{Matrix{Float64}, Nothing} = nothing,
                    use_hillslope::Bool = false,
                    downscale_hillslope_meteorology::Bool = false)
    topo_init_allocate!(topo, ncols)
    topo_init_history!(topo)
    topo_init_cold!(topo, ncols, col, lun;
                    topo_glc_mec = topo_glc_mec,
                    use_hillslope = use_hillslope,
                    downscale_hillslope_meteorology = downscale_hillslope_meteorology)
    nothing
end

# --------------------------------------------------------------------------
# InitAllocate
# --------------------------------------------------------------------------

"""
    topo_init_allocate!(topo, ncols)

Allocate topo_col (initialized to NaN) and needs_downscaling_col (initialized
to false) for `ncols` columns.

Ported from `InitAllocate` in `TopoMod.F90`.
"""
function topo_init_allocate!(topo::TopoData{FT}, ncols::Int) where {FT}
    topo.topo_col = fill(FT(NaN), ncols)
    topo.needs_downscaling_col = fill(false, ncols)
    nothing
end

# --------------------------------------------------------------------------
# InitHistory (stub — Fortran history output not applicable in Julia)
# --------------------------------------------------------------------------

"""
    topo_init_history!(topo)

Stub for history field registration. No-op in Julia port.

Ported from `InitHistory` in `TopoMod.F90`.
"""
function topo_init_history!(topo::TopoData{FT}) where {FT}
    nothing
end

# --------------------------------------------------------------------------
# InitCold
# --------------------------------------------------------------------------

"""
    topo_init_cold!(topo, ncols, col, lun;
                    topo_glc_mec=nothing,
                    use_hillslope=false,
                    downscale_hillslope_meteorology=false)

Cold-start initialization of topographic heights. For ice landunits, uses
`topo_glc_mec` surface dataset values. For hillslope columns with downscaling
enabled, uses `col.hill_elev`. Otherwise initializes to 0.

Ported from `InitCold` in `TopoMod.F90`.
"""
function topo_init_cold!(topo::TopoData, ncols::Int,
                         col::ColumnData, lun::LandunitData;
                         topo_glc_mec::Union{Matrix{Float64}, Nothing} = nothing,
                         use_hillslope::Bool = false,
                         downscale_hillslope_meteorology::Bool = false)
    for c in 1:ncols
        l = col.landunit[c]
        g = col.gridcell[c]

        if lun.itype[l] == ISTICE
            # For ice landunits, initialize topo_col based on surface dataset;
            # this will get overwritten in the run loop by values sent from CISM
            ice_class = col_itype_to_ice_class(col.itype[c])
            if topo_glc_mec !== nothing
                topo.topo_col[c] = topo_glc_mec[g, ice_class]
            else
                topo.topo_col[c] = 0.0
            end
            topo.needs_downscaling_col[c] = true
        else
            # For other landunits, arbitrarily initialize topo_col to 0 m
            if col.is_hillslope_column[c] && downscale_hillslope_meteorology
                topo.topo_col[c] = col.hill_elev[c]
                topo.needs_downscaling_col[c] = true
            else
                topo.topo_col[c] = 0.0
                topo.needs_downscaling_col[c] = false
            end
        end
    end
    nothing
end

# --------------------------------------------------------------------------
# Restart (stub — Fortran restart I/O not applicable in Julia)
# --------------------------------------------------------------------------

"""
    topo_restart!(topo)

Stub for restart variable read/write. No-op in Julia port.

Ported from `Restart` in `TopoMod.F90`.
"""
function topo_restart!(topo::TopoData{FT}) where {FT}
    nothing
end

# --------------------------------------------------------------------------
# UpdateTopo
# --------------------------------------------------------------------------

"""
    topo_update!(topo, ncols, col, lun, bounds_l, atm_topo;
                 mask_ice=nothing,
                 ice_cols_need_downscaling!=nothing,
                 update_glc2lnd_topo!=nothing,
                 update_glc_classes!=nothing,
                 use_hillslope=false,
                 downscale_hillslope_meteorology=false)

Update topographic heights each time step.

Should be called after `update_glc2lnd_fracs` and before `downscale_forcings`.

External dependencies are passed as callback functions:
- `ice_cols_need_downscaling!`: `(mask_ice, needs_downscaling_col) -> nothing`
- `update_glc2lnd_topo!`: `(topo_col, needs_downscaling_col) -> nothing`
- `update_glc_classes!`: `(topo_col) -> nothing`

Ported from `UpdateTopo` in `TopoMod.F90`.
"""
function topo_update!(topo::TopoData, ncols::Int,
                      col::ColumnData, lun::LandunitData,
                      bounds_l::UnitRange{Int},
                      atm_topo::Vector{Float64};
                      mask_ice::Union{BitVector, Nothing} = nothing,
                      ice_cols_need_downscaling!::Union{Function, Nothing} = nothing,
                      update_glc2lnd_topo!::Union{Function, Nothing} = nothing,
                      update_glc_classes!::Union{Function, Nothing} = nothing,
                      use_hillslope::Bool = false,
                      downscale_hillslope_meteorology::Bool = false)

    # Reset needs_downscaling_col each time step
    topo.needs_downscaling_col[1:ncols] .= false

    # Let GLC behavior set which ice columns need downscaling
    if ice_cols_need_downscaling! !== nothing && mask_ice !== nothing
        ice_cols_need_downscaling!(mask_ice, topo.needs_downscaling_col)
    end

    # Update topo_col from GLC and set additional needs_downscaling_col elements
    if update_glc2lnd_topo! !== nothing
        update_glc2lnd_topo!(topo.topo_col, topo.needs_downscaling_col)
    end

    # Calculate area-weighted mean hillslope elevation on each landunit
    if use_hillslope
        nl = length(bounds_l)
        mean_hillslope_elevation = zeros(Float64, last(bounds_l))
        for l_idx in bounds_l
            mhe_norm = 0.0
            for c in lun.coli[l_idx]:lun.colf[l_idx]
                if col.is_hillslope_column[c]
                    mean_hillslope_elevation[l_idx] += col.hill_elev[c] * col.hill_area[c]
                    mhe_norm += col.hill_area[c]
                end
            end
            if mhe_norm > 0.0
                mean_hillslope_elevation[l_idx] /= mhe_norm
            end
        end
    end

    # Set topo_col for non-downscaled columns and apply hillslope adjustments
    for c in 1:ncols
        if !topo.needs_downscaling_col[c]
            # Set topo to atmosphere's topographic height
            g = col.gridcell[c]
            topo.topo_col[c] = atm_topo[g]
        end

        if col.is_hillslope_column[c] && downscale_hillslope_meteorology
            l = col.landunit[c]
            topo.topo_col[c] = topo.topo_col[c] +
                (col.hill_elev[c] - mean_hillslope_elevation[l])
            topo.needs_downscaling_col[c] = true
        end
    end

    # Update GLC classes based on new topo_col
    if update_glc_classes! !== nothing
        update_glc_classes!(topo.topo_col)
    end

    nothing
end

# --------------------------------------------------------------------------
# DownscaleFilterc
# --------------------------------------------------------------------------

"""
    topo_downscale_filter(topo, ncols) -> BitVector

Returns a column-level mask: which columns need downscaling.

The main reason for this filter (vs. downscaling all columns) is that for
downscaled fields that are normalized (like longwave radiation), adding a
downscaled column in a gridcell should not change answers for other columns.

Ported from `DownscaleFilterc` in `TopoMod.F90`.
"""
function topo_downscale_filter(topo::TopoData{FT}, ncols::Int) where {FT}
    return BitVector(topo.needs_downscaling_col[1:ncols])
end

# --------------------------------------------------------------------------
# Clean
# --------------------------------------------------------------------------

"""
    topo_clean!(topo)

Reset all array fields to empty. Mirrors Fortran deallocation.

Ported from `Clean` in `TopoMod.F90`.
"""
function topo_clean!(topo::TopoData{FT}) where {FT}
    topo.topo_col = FT[]
    topo.needs_downscaling_col = Bool[]
    nothing
end
