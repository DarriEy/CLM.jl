# ==========================================================================
# Ported from: src/biogeophys/WaterBalanceType.F90
# Water balance-related variables that apply to both bulk water and tracers.
# ==========================================================================

"""
    WaterBalanceData

Water balance data structure. Holds water balance variables at the column,
patch, and gridcell levels.

Ported from `waterbalance_type` in `WaterBalanceType.F90`.
"""
Base.@kwdef mutable struct WaterBalanceData
    # --- Column-level 1D fields ---
    h2osno_old_col                       ::Vector{Float64} = Float64[]  # col snow mass for previous time step (kg/m2)
    snow_sources_col                     ::Vector{Float64} = Float64[]  # col snow sources (mm H2O/s)
    snow_sinks_col                       ::Vector{Float64} = Float64[]  # col snow sinks (mm H2O/s)
    wa_reset_nonconservation_gain_col    ::Vector{Float64} = Float64[]  # col mass gained from resetting wa_col (mm) [negative = mass lost]
    begwb_col                            ::Vector{Float64} = Float64[]  # col water mass beginning of time step
    endwb_col                            ::Vector{Float64} = Float64[]  # col water mass end of time step
    errh2o_col                           ::Vector{Float64} = Float64[]  # col water conservation error (mm H2O)
    errh2osno_col                        ::Vector{Float64} = Float64[]  # col snow water conservation error (mm H2O)

    # --- Patch-level 1D fields ---
    errh2o_patch                         ::Vector{Float64} = Float64[]  # patch water conservation error (mm H2O)

    # --- Gridcell-level 1D fields ---
    liq1_grc                             ::Vector{Float64} = Float64[]  # grc initial gridcell total h2o liq content
    liq2_grc                             ::Vector{Float64} = Float64[]  # grc post land cover change total liq content
    ice1_grc                             ::Vector{Float64} = Float64[]  # grc initial gridcell total h2o ice content
    ice2_grc                             ::Vector{Float64} = Float64[]  # grc post land cover change total ice content
    begwb_grc                            ::Vector{Float64} = Float64[]  # grc water mass beginning of time step
    endwb_grc                            ::Vector{Float64} = Float64[]  # grc water mass end of time step
end

"""
    waterbalance_init!(wb, nc, np, ng)

Allocate and initialize all fields of a `WaterBalanceData` instance for
`nc` columns, `np` patches, and `ng` gridcells.

Ported from `Init` + `InitAllocate` + `InitCold` in `WaterBalanceType.F90`.
"""
function waterbalance_init!(wb::WaterBalanceData, nc::Int, np::Int, ng::Int)
    # --- Column 1D (InitAllocate) ---
    wb.h2osno_old_col                    = fill(NaN, nc)
    wb.snow_sources_col                  = fill(NaN, nc)
    wb.snow_sinks_col                    = fill(NaN, nc)
    wb.wa_reset_nonconservation_gain_col = fill(NaN, nc)
    wb.begwb_col                         = fill(NaN, nc)
    wb.endwb_col                         = fill(NaN, nc)
    wb.errh2o_col                        = fill(NaN, nc)
    wb.errh2osno_col                     = fill(NaN, nc)

    # --- Patch 1D (InitAllocate) ---
    wb.errh2o_patch                      = fill(NaN, np)

    # --- Gridcell 1D (InitAllocate) ---
    wb.liq1_grc                          = fill(NaN, ng)
    wb.liq2_grc                          = fill(NaN, ng)
    wb.ice1_grc                          = fill(NaN, ng)
    wb.ice2_grc                          = fill(NaN, ng)
    wb.begwb_grc                         = fill(NaN, ng)
    wb.endwb_grc                         = fill(NaN, ng)

    # --- InitCold ---
    wb.wa_reset_nonconservation_gain_col .= 0.0

    return nothing
end

"""
    waterbalance_clean!(wb)

Deallocate (reset to empty) all fields of a `WaterBalanceData` instance.
"""
function waterbalance_clean!(wb::WaterBalanceData)
    # Column 1D
    wb.h2osno_old_col                    = Float64[]
    wb.snow_sources_col                  = Float64[]
    wb.snow_sinks_col                    = Float64[]
    wb.wa_reset_nonconservation_gain_col = Float64[]
    wb.begwb_col                         = Float64[]
    wb.endwb_col                         = Float64[]
    wb.errh2o_col                        = Float64[]
    wb.errh2osno_col                     = Float64[]

    # Patch 1D
    wb.errh2o_patch                      = Float64[]

    # Gridcell 1D
    wb.liq1_grc                          = Float64[]
    wb.liq2_grc                          = Float64[]
    wb.ice1_grc                          = Float64[]
    wb.ice2_grc                          = Float64[]
    wb.begwb_grc                         = Float64[]
    wb.endwb_grc                         = Float64[]

    return nothing
end

"""
    waterbalance_init_cold!(wb, bounds_col)

Initialize cold-start conditions for water balance variables.

Sets `wa_reset_nonconservation_gain_col` to 0.0 for all columns in bounds.

Ported from `InitCold` in `WaterBalanceType.F90`.
"""
function waterbalance_init_cold!(wb::WaterBalanceData, bounds_col::UnitRange{Int})
    for c in bounds_col
        wb.wa_reset_nonconservation_gain_col[c] = 0.0
    end
    return nothing
end

# ==========================================================================
# Stubs for infrastructure-dependent subroutines
# ==========================================================================

"""
    waterbalance_init_history!(wb, bounds_col)

Register water balance fields for history file output.

Ported from `InitHistory` in `WaterBalanceType.F90`.
Requires history infrastructure — stub until that module is ported.
"""
function waterbalance_init_history!(wb::WaterBalanceData, bounds_col::UnitRange{Int})
    return nothing
end
