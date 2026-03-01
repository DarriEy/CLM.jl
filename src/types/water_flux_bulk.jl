# ==========================================================================
# Ported from: src/biogeophys/WaterFluxBulkType.F90
# Water flux variables that apply only to bulk water.
# Extends the base WaterFluxData with bulk-specific fields.
# ==========================================================================

"""
    WaterFluxBulkData

Water flux bulk data structure. Holds water flux variables specific to bulk
water, in addition to all fields from the parent `WaterFluxData`.

The Fortran `waterfluxbulk_type` extends `waterflux_type` with fields for
evaporation partitioning, advective fluxes, hydraulic redistribution,
infiltration/runoff partitioning, snow melt per layer, and ET accumulation.

All water fluxes are in units of mm/s unless otherwise noted.

Ported from `waterfluxbulk_type` in `WaterFluxBulkType.F90`.
"""
Base.@kwdef mutable struct WaterFluxBulkData
    # --- Parent water flux fields (composition) ---
    wf::WaterFluxData = WaterFluxData()

    # --- Bulk-specific patch-level 1D fields ---
    qflx_snowindunload_patch  ::Vector{Float64} = Float64[]  # patch canopy snow wind unloading (mm H2O/s)
    qflx_snotempunload_patch  ::Vector{Float64} = Float64[]  # patch canopy snow temp unloading (mm H2O/s)
    qflx_ev_snow_patch        ::Vector{Float64} = Float64[]  # patch evaporation heat flux from snow (mm H2O/s) [+ to atm]
    qflx_ev_soil_patch        ::Vector{Float64} = Float64[]  # patch evaporation heat flux from soil (mm H2O/s) [+ to atm]
    qflx_ev_h2osfc_patch      ::Vector{Float64} = Float64[]  # patch evaporation heat flux from h2osfc (mm H2O/s) [+ to atm]
    qflx_hydr_redist_patch    ::Vector{Float64} = Float64[]  # patch hydraulic redistribution (mm H2O/s)

    # --- Bulk-specific column-level 1D fields ---
    qflx_phs_neg_col          ::Vector{Float64} = Float64[]  # col sum of negative hydraulic redistribution fluxes (mm H2O/s) [+]
    qflx_ev_snow_col          ::Vector{Float64} = Float64[]  # col evaporation heat flux from snow (mm H2O/s) [+ to atm]
    qflx_ev_soil_col          ::Vector{Float64} = Float64[]  # col evaporation heat flux from soil (mm H2O/s) [+ to atm]
    qflx_ev_h2osfc_col        ::Vector{Float64} = Float64[]  # col evaporation heat flux from h2osfc (mm H2O/s) [+ to atm]
    qflx_sat_excess_surf_col  ::Vector{Float64} = Float64[]  # col surface runoff due to saturated surface (mm H2O/s)
    qflx_infl_excess_col      ::Vector{Float64} = Float64[]  # col infiltration excess runoff (mm H2O/s)
    qflx_infl_excess_surf_col ::Vector{Float64} = Float64[]  # col surface runoff due to infiltration excess (mm H2O/s)
    qflx_h2osfc_surf_col      ::Vector{Float64} = Float64[]  # col surface water runoff (mm H2O/s)
    qflx_in_soil_col          ::Vector{Float64} = Float64[]  # col surface input to soil (mm/s)
    qflx_in_soil_limited_col  ::Vector{Float64} = Float64[]  # col surface input to soil, limited by max infiltration rate (mm/s)
    qflx_h2osfc_drain_col     ::Vector{Float64} = Float64[]  # col bottom drainage from h2osfc (mm/s)
    qflx_top_soil_to_h2osfc_col ::Vector{Float64} = Float64[]  # col portion of qflx_top_soil going to h2osfc, minus evaporation (mm/s)
    qflx_in_h2osfc_col        ::Vector{Float64} = Float64[]  # col total surface input to h2osfc
    qflx_deficit_col          ::Vector{Float64} = Float64[]  # col water deficit to keep non-negative liquid water content (mm H2O)
    AnnET                     ::Vector{Float64} = Float64[]  # col annual average ET flux (mm H2O/s)

    # --- Bulk-specific column-level 2D fields ---
    qflx_adv_col              ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col advective flux across soil layer interfaces (mm H2O/s) [+ downward] (0:nlevsoi)
    qflx_rootsoi_col          ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col root and soil water exchange (mm H2O/s) [+ into root] (1:nlevsoi)
    qflx_snomelt_lyr_col      ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col snow melt in each layer (mm H2O/s) (-nlevsno+1:0)
    qflx_drain_vr_col         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col liquid water lost as drainage (m/time step) (1:nlevsoi)
end

"""
    waterfluxbulk_init!(wfb, nc, np, nl, ng)

Allocate and initialize all fields of a `WaterFluxBulkData` instance for
`nc` columns, `np` patches, `nl` landunits, and `ng` gridcells.
Calls parent `waterflux_init!` plus allocates bulk-specific fields.

Ported from `InitBulk` + `InitBulkAllocate` in `WaterFluxBulkType.F90`.
"""
function waterfluxbulk_init!(wfb::WaterFluxBulkData, nc::Int, np::Int, nl::Int, ng::Int)
    nlevsno = varpar.nlevsno
    nlevsoi = varpar.nlevsoi

    # Initialize parent fields
    waterflux_init!(wfb.wf, nc, np, nl, ng)

    # --- Patch 1D ---
    wfb.qflx_snowindunload_patch  = fill(NaN, np)
    wfb.qflx_snotempunload_patch  = fill(NaN, np)
    wfb.qflx_ev_snow_patch        = fill(NaN, np)
    wfb.qflx_ev_soil_patch        = fill(NaN, np)
    wfb.qflx_ev_h2osfc_patch      = fill(NaN, np)
    wfb.qflx_hydr_redist_patch    = fill(NaN, np)

    # --- Column 1D ---
    wfb.qflx_phs_neg_col          = fill(NaN, nc)
    wfb.qflx_ev_snow_col          = fill(NaN, nc)
    wfb.qflx_ev_soil_col          = fill(NaN, nc)
    wfb.qflx_ev_h2osfc_col        = fill(NaN, nc)
    wfb.qflx_sat_excess_surf_col  = fill(NaN, nc)
    wfb.qflx_infl_excess_col      = fill(NaN, nc)
    wfb.qflx_infl_excess_surf_col = fill(NaN, nc)
    wfb.qflx_h2osfc_surf_col      = fill(NaN, nc)
    wfb.qflx_in_soil_col          = fill(NaN, nc)
    wfb.qflx_in_soil_limited_col  = fill(NaN, nc)
    wfb.qflx_h2osfc_drain_col     = fill(NaN, nc)
    wfb.qflx_top_soil_to_h2osfc_col = fill(NaN, nc)
    wfb.qflx_in_h2osfc_col        = fill(NaN, nc)
    wfb.qflx_deficit_col          = fill(NaN, nc)
    wfb.AnnET                     = fill(NaN, nc)

    # --- Column 2D ---
    wfb.qflx_adv_col              = fill(NaN, nc, nlevsoi + 1)  # (0:nlevsoi) → nlevsoi+1 columns
    wfb.qflx_rootsoi_col          = fill(NaN, nc, nlevsoi)      # (1:nlevsoi)
    wfb.qflx_snomelt_lyr_col      = fill(NaN, nc, nlevsno)      # (-nlevsno+1:0)
    wfb.qflx_drain_vr_col         = fill(NaN, nc, nlevsoi)      # (1:nlevsoi)

    return nothing
end

"""
    waterfluxbulk_clean!(wfb)

Deallocate (reset to empty) all fields of a `WaterFluxBulkData` instance.
"""
function waterfluxbulk_clean!(wfb::WaterFluxBulkData)
    waterflux_clean!(wfb.wf)

    # Patch 1D
    wfb.qflx_snowindunload_patch  = Float64[]
    wfb.qflx_snotempunload_patch  = Float64[]
    wfb.qflx_ev_snow_patch        = Float64[]
    wfb.qflx_ev_soil_patch        = Float64[]
    wfb.qflx_ev_h2osfc_patch      = Float64[]
    wfb.qflx_hydr_redist_patch    = Float64[]

    # Column 1D
    wfb.qflx_phs_neg_col          = Float64[]
    wfb.qflx_ev_snow_col          = Float64[]
    wfb.qflx_ev_soil_col          = Float64[]
    wfb.qflx_ev_h2osfc_col        = Float64[]
    wfb.qflx_sat_excess_surf_col  = Float64[]
    wfb.qflx_infl_excess_col      = Float64[]
    wfb.qflx_infl_excess_surf_col = Float64[]
    wfb.qflx_h2osfc_surf_col      = Float64[]
    wfb.qflx_in_soil_col          = Float64[]
    wfb.qflx_in_soil_limited_col  = Float64[]
    wfb.qflx_h2osfc_drain_col     = Float64[]
    wfb.qflx_top_soil_to_h2osfc_col = Float64[]
    wfb.qflx_in_h2osfc_col        = Float64[]
    wfb.qflx_deficit_col          = Float64[]
    wfb.AnnET                     = Float64[]

    # Column 2D
    wfb.qflx_adv_col              = Matrix{Float64}(undef, 0, 0)
    wfb.qflx_rootsoi_col          = Matrix{Float64}(undef, 0, 0)
    wfb.qflx_snomelt_lyr_col      = Matrix{Float64}(undef, 0, 0)
    wfb.qflx_drain_vr_col         = Matrix{Float64}(undef, 0, 0)

    return nothing
end

"""
    waterfluxbulk_init_cold!(wfb, bounds_col, bounds_patch)

Initialize cold-start conditions for bulk-specific water flux variables.

Sets `qflx_snowindunload_patch`, `qflx_snotempunload_patch` to 0 for patches,
and `qflx_phs_neg_col`, `qflx_h2osfc_surf_col` to 0 for columns.

Ported from `InitBulkCold` in `WaterFluxBulkType.F90`.
"""
function waterfluxbulk_init_cold!(wfb::WaterFluxBulkData,
                                    bounds_col::UnitRange{Int},
                                    bounds_patch::UnitRange{Int})
    for p in bounds_patch
        wfb.qflx_snowindunload_patch[p] = 0.0
        wfb.qflx_snotempunload_patch[p] = 0.0
    end

    for c in bounds_col
        wfb.qflx_phs_neg_col[c]      = 0.0
        wfb.qflx_h2osfc_surf_col[c]  = 0.0
    end

    return nothing
end

# ==========================================================================
# Stubs for infrastructure-dependent subroutines
# ==========================================================================

"""
    waterfluxbulk_init_history!(wfb, bounds_col)

Register water flux bulk fields for history file output.

Ported from `InitBulkHistory` in `WaterFluxBulkType.F90`.
Requires history infrastructure — stub until that module is ported.
"""
function waterfluxbulk_init_history!(wfb::WaterFluxBulkData,
                                      bounds_col::UnitRange{Int})
    return nothing
end

"""
    waterfluxbulk_restart!(wfb, bounds_col; flag="read")

Read/write bulk water flux from/to restart file.

Ported from `RestartBulk` in `WaterFluxBulkType.F90`.
Requires NetCDF/restart infrastructure — stub until that module is ported.
"""
function waterfluxbulk_restart!(wfb::WaterFluxBulkData,
                                  bounds_col::UnitRange{Int};
                                  flag::String = "read")
    return nothing
end

"""
    waterfluxbulk_init_acc_buffer!(wfb, bounds_col)

Initialize accumulation buffer for AnnET (365-day running mean of total ET).

Ported from `InitAccBuffer` in `WaterFluxBulkType.F90`.
Requires accumulation infrastructure — stub until that module is ported.
"""
function waterfluxbulk_init_acc_buffer!(wfb::WaterFluxBulkData,
                                          bounds_col::UnitRange{Int})
    return nothing
end

"""
    waterfluxbulk_init_acc_vars!(wfb, bounds_col)

Initialize accumulation variables from restart.

Ported from `InitAccVars` in `WaterFluxBulkType.F90`.
Requires accumulation infrastructure — stub until that module is ported.
"""
function waterfluxbulk_init_acc_vars!(wfb::WaterFluxBulkData,
                                        bounds_col::UnitRange{Int})
    return nothing
end

"""
    waterfluxbulk_update_acc_vars!(wfb, bounds_col)

Update accumulation variables (AnnET from qflx_evap_tot_col).

Ported from `UpdateAccVars` in `WaterFluxBulkType.F90`.
Requires accumulation infrastructure — stub until that module is ported.
"""
function waterfluxbulk_update_acc_vars!(wfb::WaterFluxBulkData,
                                          bounds_col::UnitRange{Int})
    return nothing
end
