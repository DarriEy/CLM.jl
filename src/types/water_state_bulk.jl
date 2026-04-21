# ==========================================================================
# Ported from: src/biogeophys/WaterStateBulkType.F90
# Water state variables that apply only to bulk water.
# Extends the base WaterStateData with bulk-specific fields.
# ==========================================================================

"""
    WaterStateBulkData

Water state bulk data structure. Holds water state variables specific to bulk
water, in addition to all fields from the parent `WaterStateData`.

The Fortran `waterstatebulk_type` extends `waterstate_type` with two
column-level fields: `snow_persistence_col` and `int_snow_col`.

Ported from `waterstatebulk_type` in `WaterStateBulkType.F90`.
"""
Base.@kwdef mutable struct WaterStateBulkData{FT<:Real}
    # --- Parent water state fields (composition) ---
    ws::WaterStateData = WaterStateData()

    # --- Bulk-specific column-level 1D fields ---
    snow_persistence_col ::Vector{FT} = Float64[]   # col length of time that ground has had non-zero snow thickness (sec)
    int_snow_col         ::Vector{FT} = Float64[]   # col integrated snowfall (mm H2O)
end

"""
    waterstatebulk_init!(wsb, nc, np, nl, ng)

Allocate and initialize all fields of a `WaterStateBulkData` instance for
`nc` columns, `np` patches, `nl` landunits, and `ng` gridcells.
Calls parent `waterstate_init!` plus allocates bulk-specific fields.

Ported from `InitBulk` + `InitBulkAllocate` in `WaterStateBulkType.F90`.
"""
function waterstatebulk_init!(wsb::WaterStateBulkData{FT}, nc::Int, np::Int, nl::Int, ng::Int) where {FT}
    # Initialize parent fields
    waterstate_init!(wsb.ws, nc, np, nl, ng)

    # Bulk-specific fields (InitBulkAllocate)
    wsb.snow_persistence_col = fill(zero(FT), nc)
    wsb.int_snow_col         = fill(zero(FT), nc)

    return nothing
end

"""
    waterstatebulk_clean!(wsb)

Deallocate (reset to empty) all fields of a `WaterStateBulkData` instance.
"""
function waterstatebulk_clean!(wsb::WaterStateBulkData{FT}) where {FT}
    waterstate_clean!(wsb.ws)

    wsb.snow_persistence_col = FT[]
    wsb.int_snow_col         = FT[]

    return nothing
end

"""
    waterstatebulk_init_cold!(wsb, bounds_col; h2osno_input_col)

Initialize cold-start conditions for bulk-specific water state variables.

Sets `int_snow_col` to `h2osno_input_col` and `snow_persistence_col` to 0.

Ported from `InitBulkCold` in `WaterStateBulkType.F90`.
"""
function waterstatebulk_init_cold!(wsb::WaterStateBulkData,
                                     bounds_col::UnitRange{Int};
                                     h2osno_input_col::Vector{<:Real})
    for c in bounds_col
        wsb.int_snow_col[c]         = h2osno_input_col[c]
        wsb.snow_persistence_col[c] = 0.0
    end

    return nothing
end

# ==========================================================================
# Stubs for infrastructure-dependent subroutines
# ==========================================================================

"""
    waterstatebulk_init_history!(wsb, bounds_col)

Register water state bulk fields for history file output.

Ported from `InitBulkHistory` in `WaterStateBulkType.F90`.
Requires history infrastructure — stub until that module is ported.
"""
function waterstatebulk_init_history!(wsb::WaterStateBulkData,
                                       bounds_col::UnitRange{Int})
    return nothing
end

"""
    waterstatebulk_restart!(wsb, bounds_col; flag="read")

Read/write bulk water state from/to restart file.

Ported from `RestartBulk` in `WaterStateBulkType.F90`.
Requires NetCDF/restart infrastructure — stub until that module is ported.
"""
function waterstatebulk_restart!(wsb::WaterStateBulkData,
                                   bounds_col::UnitRange{Int};
                                   flag::String = "read")
    return nothing
end
