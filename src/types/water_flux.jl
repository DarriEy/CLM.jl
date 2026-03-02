# ==========================================================================
# Ported from: src/biogeophys/WaterFluxType.F90
# Water flux variables that apply to both bulk water and water tracers.
# ==========================================================================

"""
    WaterFluxData

Water flux data structure. Holds water flux variables at the column, patch,
landunit, and gridcell levels, including canopy throughfall, evaporation,
transpiration, runoff, snow melt, irrigation, and dynamic land cover fluxes.

All water fluxes are in units of mm/s unless otherwise noted.

Ported from `waterflux_type` in `WaterFluxType.F90`.
"""
Base.@kwdef mutable struct WaterFluxData
    # --- Patch-level 1D fields ---
    qflx_through_snow_patch             ::Vector{Float64} = Float64[]  # patch canopy throughfall of snow (mm H2O/s)
    qflx_through_liq_patch              ::Vector{Float64} = Float64[]  # patch canopy throughfall of liquid (rain+irrigation) (mm H2O/s)
    qflx_intercepted_snow_patch         ::Vector{Float64} = Float64[]  # patch canopy interception of snow (mm H2O/s)
    qflx_intercepted_liq_patch          ::Vector{Float64} = Float64[]  # patch canopy interception of liquid (rain+irrigation) (mm H2O/s)
    qflx_snocanfall_patch               ::Vector{Float64} = Float64[]  # patch rate of excess canopy snow falling off canopy (mm H2O/s)
    qflx_liqcanfall_patch               ::Vector{Float64} = Float64[]  # patch rate of excess canopy liquid falling off canopy (mm H2O/s)
    qflx_snow_unload_patch              ::Vector{Float64} = Float64[]  # patch rate of canopy snow unloading (mm H2O/s)
    qflx_solidevap_from_top_layer_patch ::Vector{Float64} = Float64[]  # patch rate of ice sublimation from top layer (mm H2O/s) [+]
    qflx_tran_veg_patch                 ::Vector{Float64} = Float64[]  # patch vegetation transpiration (mm H2O/s) (+ = to atm)
    qflx_liqdew_to_top_layer_patch      ::Vector{Float64} = Float64[]  # patch rate of liquid dew on top layer (mm H2O/s) [+]
    qflx_soliddew_to_top_layer_patch    ::Vector{Float64} = Float64[]  # patch rate of frost on top layer (mm H2O/s) [+]
    qflx_evap_veg_patch                 ::Vector{Float64} = Float64[]  # patch vegetation evaporation (mm H2O/s) (+ = to atm)
    qflx_evap_can_patch                 ::Vector{Float64} = Float64[]  # patch evaporation from leaves and stems (mm H2O/s) (+ = to atm)
    qflx_evap_soi_patch                 ::Vector{Float64} = Float64[]  # patch soil evaporation (mm H2O/s) (+ = to atm)
    qflx_evap_tot_patch                 ::Vector{Float64} = Float64[]  # patch total evaporation (mm H2O/s)
    qflx_liqevap_from_top_layer_patch   ::Vector{Float64} = Float64[]  # patch rate of liquid evaporation from top layer (mm H2O/s) [+]
    qflx_ev_snow_patch                  ::Vector{Float64} = Float64[]  # patch evaporation heat flux from snow (mm H2O/s)
    qflx_ev_soil_patch                  ::Vector{Float64} = Float64[]  # patch evaporation heat flux from soil (mm H2O/s)
    qflx_ev_h2osfc_patch                ::Vector{Float64} = Float64[]  # patch evaporation heat flux from surface water (mm H2O/s)
    qflx_irrig_drip_patch               ::Vector{Float64} = Float64[]  # patch drip irrigation (mm H2O/s)
    qflx_irrig_sprinkler_patch          ::Vector{Float64} = Float64[]  # patch sprinkler irrigation (mm H2O/s)

    # --- Column-level 1D fields ---
    qflx_liq_grnd_col                   ::Vector{Float64} = Float64[]  # col liquid on ground after interception (mm H2O/s) [+]
    qflx_snow_grnd_col                  ::Vector{Float64} = Float64[]  # col snow on ground after interception (mm H2O/s) [+]
    qflx_rain_plus_snomelt_col          ::Vector{Float64} = Float64[]  # col rain plus snow melt on soil (mm/s)
    qflx_solidevap_from_top_layer_col   ::Vector{Float64} = Float64[]  # col rate of ice sublimation from top layer (mm H2O/s) [+]
    qflx_snwcp_liq_col                  ::Vector{Float64} = Float64[]  # col excess liquid from snow capping (mm H2O/s)
    qflx_snwcp_ice_col                  ::Vector{Float64} = Float64[]  # col excess solid from snow capping (mm H2O/s)
    qflx_snwcp_discarded_liq_col        ::Vector{Float64} = Float64[]  # col discarded excess liquid from snow capping (mm H2O/s)
    qflx_snwcp_discarded_ice_col        ::Vector{Float64} = Float64[]  # col discarded excess solid from snow capping (mm H2O/s)
    qflx_glcice_col                     ::Vector{Float64} = Float64[]  # col net glacial ice flux (mm H2O/s)
    qflx_glcice_frz_col                 ::Vector{Float64} = Float64[]  # col ice growth (mm H2O/s) [+]
    qflx_glcice_melt_col                ::Vector{Float64} = Float64[]  # col ice melt (mm H2O/s) [+]
    qflx_glcice_dyn_water_flux_col      ::Vector{Float64} = Float64[]  # col dynamic water flux for glc routing (mm H2O/s)
    qflx_tran_veg_col                   ::Vector{Float64} = Float64[]  # col vegetation transpiration (mm H2O/s) (+ = to atm)
    qflx_evap_veg_col                   ::Vector{Float64} = Float64[]  # col vegetation evaporation (mm H2O/s) (+ = to atm)
    qflx_evap_can_col                   ::Vector{Float64} = Float64[]  # col evaporation from leaves and stems (mm H2O/s) (+ = to atm)
    qflx_evap_soi_col                   ::Vector{Float64} = Float64[]  # col soil evaporation (mm H2O/s) (+ = to atm)
    qflx_evap_tot_col                   ::Vector{Float64} = Float64[]  # col total evaporation (mm H2O/s)
    qflx_liqevap_from_top_layer_col     ::Vector{Float64} = Float64[]  # col rate of liquid evaporation from top layer (mm H2O/s) [+]
    qflx_liqdew_to_top_layer_col        ::Vector{Float64} = Float64[]  # col rate of liquid dew on top layer (mm H2O/s) [+]
    qflx_soliddew_to_top_layer_col      ::Vector{Float64} = Float64[]  # col rate of frost on top layer (mm H2O/s) [+]
    qflx_infl_col                       ::Vector{Float64} = Float64[]  # col infiltration (mm H2O/s)
    qflx_surf_col                       ::Vector{Float64} = Float64[]  # col total surface runoff (mm H2O/s)
    qflx_drain_col                      ::Vector{Float64} = Float64[]  # col sub-surface runoff (mm H2O/s)
    qflx_drain_perched_col              ::Vector{Float64} = Float64[]  # col perched water table drainage (mm H2O/s)
    qflx_latflow_in_col                 ::Vector{Float64} = Float64[]  # col hillslope lateral flow input (mm/s)
    qflx_latflow_out_col                ::Vector{Float64} = Float64[]  # col hillslope lateral flow output (mm/s)
    volumetric_discharge_col            ::Vector{Float64} = Float64[]  # col hillslope discharge (m3/s)
    qflx_top_soil_col                   ::Vector{Float64} = Float64[]  # col net water input into soil from top (mm/s)
    qflx_floodc_col                     ::Vector{Float64} = Float64[]  # col flood water flux
    qflx_sl_top_soil_col                ::Vector{Float64} = Float64[]  # col liquid+ice from above soil to top soil or qrgwl (mm H2O/s)
    qflx_snomelt_col                    ::Vector{Float64} = Float64[]  # col snow melt (mm H2O/s)
    qflx_qrgwl_col                      ::Vector{Float64} = Float64[]  # col surface runoff at glaciers/wetlands/lakes
    qflx_runoff_col                     ::Vector{Float64} = Float64[]  # col total runoff (mm H2O/s)
    qflx_runoff_r_col                   ::Vector{Float64} = Float64[]  # col rural total runoff (mm H2O/s)
    qflx_runoff_u_col                   ::Vector{Float64} = Float64[]  # col urban total runoff (mm H2O/s)
    qflx_rsub_sat_col                   ::Vector{Float64} = Float64[]  # col soil saturation excess (mm/s)
    qflx_snofrz_col                     ::Vector{Float64} = Float64[]  # col column-integrated snow freezing rate (kg m-2 s-1)
    qflx_snow_drain_col                 ::Vector{Float64} = Float64[]  # col drainage from snow pack
    qflx_ice_runoff_snwcp_col           ::Vector{Float64} = Float64[]  # col solid runoff from snow capping (mm H2O/s)
    qflx_ice_runoff_xs_col              ::Vector{Float64} = Float64[]  # col solid runoff from excess ice (mm H2O/s)
    qflx_h2osfc_to_ice_col              ::Vector{Float64} = Float64[]  # col conversion of h2osfc to ice
    qflx_snow_h2osfc_col                ::Vector{Float64} = Float64[]  # col snow falling on surface water
    qflx_too_small_h2osfc_to_soil_col   ::Vector{Float64} = Float64[]  # col h2osfc transferred to soil if below threshold (mm H2O/s)
    qflx_sfc_irrig_col                  ::Vector{Float64} = Float64[]  # col surface irrigation flux (mm H2O/s) [+]
    qflx_gw_uncon_irrig_col             ::Vector{Float64} = Float64[]  # col unconfined groundwater irrigation flux (mm H2O/s)
    qflx_gw_con_irrig_col               ::Vector{Float64} = Float64[]  # col confined groundwater irrigation flux (mm H2O/s)

    # --- Column-level 2D fields ---
    qflx_snofrz_lyr_col                 ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col snow freezing rate per layer (kg m-2 s-1) (-nlevsno+1:0)
    qflx_snomelt_lyr_col                ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col snow melt rate per layer (kg m-2 s-1) (-nlevsno+1:0)
    qflx_snow_percolation_col           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col liquid percolation from snow layer (mm H2O/s) (-nlevsno+1:0)
    qflx_gw_uncon_irrig_lyr_col         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col unconfined gw irrigation by layer (mm H2O/s) (1:nlevsoi)

    # --- Landunit-level 1D fields ---
    volumetric_streamflow_lun           ::Vector{Float64} = Float64[]  # lun stream discharge (m3/s)

    # --- Gridcell-level 1D fields ---
    qflx_liq_dynbal_grc                 ::Vector{Float64} = Float64[]  # grc liq dynamic land cover change runoff flux
    qflx_ice_dynbal_grc                 ::Vector{Float64} = Float64[]  # grc ice dynamic land cover change runoff flux
end

"""
    waterflux_init!(wf::WaterFluxData, nc::Int, np::Int, nl::Int, ng::Int)

Allocate and initialize all fields of a `WaterFluxData` instance for
`nc` columns, `np` patches, `nl` landunits, and `ng` gridcells.
Real fields are initialized to `NaN`.

Ported from `waterflux_type%InitAllocate` in `WaterFluxType.F90`.
"""
function waterflux_init!(wf::WaterFluxData, nc::Int, np::Int, nl::Int, ng::Int)
    nlevsno = varpar.nlevsno
    nlevsoi = varpar.nlevsoi

    # --- Patch 1D ---
    wf.qflx_through_snow_patch             = fill(NaN, np)
    wf.qflx_through_liq_patch              = fill(NaN, np)
    wf.qflx_intercepted_snow_patch         = fill(NaN, np)
    wf.qflx_intercepted_liq_patch          = fill(NaN, np)
    wf.qflx_snocanfall_patch               = fill(NaN, np)
    wf.qflx_liqcanfall_patch               = fill(NaN, np)
    wf.qflx_snow_unload_patch              = fill(NaN, np)
    wf.qflx_solidevap_from_top_layer_patch = fill(0.0, np)  # ival = 0.0 in Fortran
    wf.qflx_tran_veg_patch                 = fill(NaN, np)
    wf.qflx_liqdew_to_top_layer_patch      = fill(NaN, np)
    wf.qflx_soliddew_to_top_layer_patch    = fill(NaN, np)
    wf.qflx_evap_veg_patch                 = fill(NaN, np)
    wf.qflx_evap_can_patch                 = fill(NaN, np)
    wf.qflx_evap_soi_patch                 = fill(NaN, np)
    wf.qflx_evap_tot_patch                 = fill(NaN, np)
    wf.qflx_liqevap_from_top_layer_patch   = fill(NaN, np)
    wf.qflx_ev_snow_patch                  = fill(NaN, np)
    wf.qflx_ev_soil_patch                  = fill(NaN, np)
    wf.qflx_ev_h2osfc_patch                = fill(NaN, np)
    wf.qflx_irrig_drip_patch               = fill(NaN, np)
    wf.qflx_irrig_sprinkler_patch          = fill(NaN, np)

    # --- Column 1D ---
    wf.qflx_liq_grnd_col                   = fill(NaN, nc)
    wf.qflx_snow_grnd_col                  = fill(NaN, nc)
    wf.qflx_rain_plus_snomelt_col          = fill(NaN, nc)
    wf.qflx_solidevap_from_top_layer_col   = fill(0.0, nc)  # ival = 0.0 in Fortran
    wf.qflx_snwcp_liq_col                  = fill(NaN, nc)
    wf.qflx_snwcp_ice_col                  = fill(NaN, nc)
    wf.qflx_snwcp_discarded_liq_col        = fill(NaN, nc)
    wf.qflx_snwcp_discarded_ice_col        = fill(NaN, nc)
    wf.qflx_glcice_col                     = fill(NaN, nc)
    wf.qflx_glcice_frz_col                 = fill(NaN, nc)
    wf.qflx_glcice_melt_col                = fill(NaN, nc)
    wf.qflx_glcice_dyn_water_flux_col      = fill(NaN, nc)
    wf.qflx_tran_veg_col                   = fill(NaN, nc)
    wf.qflx_evap_veg_col                   = fill(NaN, nc)
    wf.qflx_evap_can_col                   = fill(NaN, nc)
    wf.qflx_evap_soi_col                   = fill(NaN, nc)
    wf.qflx_evap_tot_col                   = fill(NaN, nc)
    wf.qflx_liqevap_from_top_layer_col     = fill(NaN, nc)
    wf.qflx_liqdew_to_top_layer_col        = fill(NaN, nc)
    wf.qflx_soliddew_to_top_layer_col      = fill(NaN, nc)
    wf.qflx_infl_col                       = fill(NaN, nc)
    wf.qflx_surf_col                       = fill(NaN, nc)
    wf.qflx_drain_col                      = fill(NaN, nc)
    wf.qflx_drain_perched_col              = fill(NaN, nc)
    wf.qflx_latflow_in_col                 = fill(NaN, nc)
    wf.qflx_latflow_out_col                = fill(NaN, nc)
    wf.volumetric_discharge_col            = fill(NaN, nc)
    wf.qflx_top_soil_col                   = fill(NaN, nc)
    wf.qflx_floodc_col                     = fill(NaN, nc)
    wf.qflx_sl_top_soil_col                = fill(NaN, nc)
    wf.qflx_snomelt_col                    = fill(NaN, nc)
    wf.qflx_qrgwl_col                      = fill(NaN, nc)
    wf.qflx_runoff_col                     = fill(NaN, nc)
    wf.qflx_runoff_r_col                   = fill(NaN, nc)
    wf.qflx_runoff_u_col                   = fill(NaN, nc)
    wf.qflx_rsub_sat_col                   = fill(NaN, nc)
    wf.qflx_snofrz_col                     = fill(NaN, nc)
    wf.qflx_snow_drain_col                 = fill(NaN, nc)
    wf.qflx_ice_runoff_snwcp_col           = fill(NaN, nc)
    wf.qflx_ice_runoff_xs_col              = fill(NaN, nc)
    wf.qflx_h2osfc_to_ice_col             = fill(NaN, nc)
    wf.qflx_snow_h2osfc_col               = fill(NaN, nc)
    wf.qflx_too_small_h2osfc_to_soil_col  = fill(NaN, nc)
    wf.qflx_sfc_irrig_col                 = fill(NaN, nc)
    wf.qflx_gw_uncon_irrig_col            = fill(NaN, nc)
    wf.qflx_gw_con_irrig_col              = fill(NaN, nc)

    # --- Column 2D ---
    wf.qflx_snofrz_lyr_col               = fill(NaN, nc, nlevsno)       # (-nlevsno+1:0)
    wf.qflx_snomelt_lyr_col              = fill(NaN, nc, nlevsno)       # (-nlevsno+1:0)
    wf.qflx_snow_percolation_col          = fill(NaN, nc, nlevsno)       # (-nlevsno+1:0)
    wf.qflx_gw_uncon_irrig_lyr_col        = fill(NaN, nc, nlevsoi)       # (1:nlevsoi)

    # --- Landunit 1D ---
    wf.volumetric_streamflow_lun           = fill(NaN, nl)

    # --- Gridcell 1D ---
    wf.qflx_liq_dynbal_grc                = fill(NaN, ng)
    wf.qflx_ice_dynbal_grc                = fill(NaN, ng)

    return nothing
end

"""
    waterflux_clean!(wf::WaterFluxData)

Deallocate (reset to empty) all fields of a `WaterFluxData` instance.
"""
function waterflux_clean!(wf::WaterFluxData)
    # Patch 1D
    wf.qflx_through_snow_patch             = Float64[]
    wf.qflx_through_liq_patch              = Float64[]
    wf.qflx_intercepted_snow_patch         = Float64[]
    wf.qflx_intercepted_liq_patch          = Float64[]
    wf.qflx_snocanfall_patch               = Float64[]
    wf.qflx_liqcanfall_patch               = Float64[]
    wf.qflx_snow_unload_patch              = Float64[]
    wf.qflx_solidevap_from_top_layer_patch = Float64[]
    wf.qflx_tran_veg_patch                 = Float64[]
    wf.qflx_liqdew_to_top_layer_patch      = Float64[]
    wf.qflx_soliddew_to_top_layer_patch    = Float64[]
    wf.qflx_evap_veg_patch                 = Float64[]
    wf.qflx_evap_can_patch                 = Float64[]
    wf.qflx_evap_soi_patch                 = Float64[]
    wf.qflx_evap_tot_patch                 = Float64[]
    wf.qflx_liqevap_from_top_layer_patch   = Float64[]
    wf.qflx_ev_snow_patch                  = Float64[]
    wf.qflx_ev_soil_patch                  = Float64[]
    wf.qflx_ev_h2osfc_patch                = Float64[]
    wf.qflx_irrig_drip_patch               = Float64[]
    wf.qflx_irrig_sprinkler_patch          = Float64[]

    # Column 1D
    wf.qflx_liq_grnd_col                   = Float64[]
    wf.qflx_snow_grnd_col                  = Float64[]
    wf.qflx_rain_plus_snomelt_col          = Float64[]
    wf.qflx_solidevap_from_top_layer_col   = Float64[]
    wf.qflx_snwcp_liq_col                  = Float64[]
    wf.qflx_snwcp_ice_col                  = Float64[]
    wf.qflx_snwcp_discarded_liq_col        = Float64[]
    wf.qflx_snwcp_discarded_ice_col        = Float64[]
    wf.qflx_glcice_col                     = Float64[]
    wf.qflx_glcice_frz_col                 = Float64[]
    wf.qflx_glcice_melt_col                = Float64[]
    wf.qflx_glcice_dyn_water_flux_col      = Float64[]
    wf.qflx_tran_veg_col                   = Float64[]
    wf.qflx_evap_veg_col                   = Float64[]
    wf.qflx_evap_can_col                   = Float64[]
    wf.qflx_evap_soi_col                   = Float64[]
    wf.qflx_evap_tot_col                   = Float64[]
    wf.qflx_liqevap_from_top_layer_col     = Float64[]
    wf.qflx_liqdew_to_top_layer_col        = Float64[]
    wf.qflx_soliddew_to_top_layer_col      = Float64[]
    wf.qflx_infl_col                       = Float64[]
    wf.qflx_surf_col                       = Float64[]
    wf.qflx_drain_col                      = Float64[]
    wf.qflx_drain_perched_col              = Float64[]
    wf.qflx_latflow_in_col                 = Float64[]
    wf.qflx_latflow_out_col                = Float64[]
    wf.volumetric_discharge_col            = Float64[]
    wf.qflx_top_soil_col                   = Float64[]
    wf.qflx_floodc_col                     = Float64[]
    wf.qflx_sl_top_soil_col                = Float64[]
    wf.qflx_snomelt_col                    = Float64[]
    wf.qflx_qrgwl_col                      = Float64[]
    wf.qflx_runoff_col                     = Float64[]
    wf.qflx_runoff_r_col                   = Float64[]
    wf.qflx_runoff_u_col                   = Float64[]
    wf.qflx_rsub_sat_col                   = Float64[]
    wf.qflx_snofrz_col                     = Float64[]
    wf.qflx_snow_drain_col                 = Float64[]
    wf.qflx_ice_runoff_snwcp_col           = Float64[]
    wf.qflx_ice_runoff_xs_col              = Float64[]
    wf.qflx_h2osfc_to_ice_col             = Float64[]
    wf.qflx_snow_h2osfc_col               = Float64[]
    wf.qflx_too_small_h2osfc_to_soil_col  = Float64[]
    wf.qflx_sfc_irrig_col                 = Float64[]
    wf.qflx_gw_uncon_irrig_col            = Float64[]
    wf.qflx_gw_con_irrig_col              = Float64[]

    # Column 2D
    wf.qflx_snofrz_lyr_col               = Matrix{Float64}(undef, 0, 0)
    wf.qflx_snomelt_lyr_col              = Matrix{Float64}(undef, 0, 0)
    wf.qflx_snow_percolation_col          = Matrix{Float64}(undef, 0, 0)
    wf.qflx_gw_uncon_irrig_lyr_col        = Matrix{Float64}(undef, 0, 0)

    # Landunit 1D
    wf.volumetric_streamflow_lun           = Float64[]

    # Gridcell 1D
    wf.qflx_liq_dynbal_grc                = Float64[]
    wf.qflx_ice_dynbal_grc                = Float64[]

    return nothing
end

"""
    waterflux_init_cold!(wf::WaterFluxData,
                          bounds_col::UnitRange{Int},
                          bounds_patch::UnitRange{Int},
                          bounds_lun::UnitRange{Int};
                          landunit_col::Vector{Int},
                          lun_itype::Vector{Int},
                          use_hillslope_routing::Bool = false)

Initialize cold-start conditions for water flux variables.

Ported from `waterflux_type%InitCold` in `WaterFluxType.F90`.
"""
function waterflux_init_cold!(wf::WaterFluxData,
                               bounds_col::UnitRange{Int},
                               bounds_patch::UnitRange{Int},
                               bounds_lun::UnitRange{Int};
                               landunit_col::Vector{Int},
                               lun_itype::Vector{Int},
                               use_hillslope_routing::Bool = false)
    # Patch-level zeros
    for p in bounds_patch
        wf.qflx_snocanfall_patch[p]               = 0.0
        wf.qflx_liqcanfall_patch[p]               = 0.0
        wf.qflx_snow_unload_patch[p]              = 0.0
        wf.qflx_liqevap_from_top_layer_patch[p]   = 0.0
        wf.qflx_liqdew_to_top_layer_patch[p]      = 0.0
        wf.qflx_soliddew_to_top_layer_patch[p]    = 0.0
        wf.qflx_irrig_drip_patch[p]               = 0.0
        wf.qflx_irrig_sprinkler_patch[p]          = 0.0
        wf.qflx_tran_veg_patch[p]                 = 0.0
        wf.qflx_evap_veg_patch[p]                 = 0.0
    end

    # Column-level zeros
    for c in bounds_col
        wf.qflx_sfc_irrig_col[c]                  = 0.0
        wf.qflx_gw_uncon_irrig_col[c]             = 0.0
        wf.qflx_gw_con_irrig_col[c]               = 0.0
        wf.qflx_liqevap_from_top_layer_col[c]     = 0.0
        wf.qflx_liqdew_to_top_layer_col[c]        = 0.0
        wf.qflx_soliddew_to_top_layer_col[c]      = 0.0
        wf.qflx_snow_drain_col[c]                 = 0.0
        wf.qflx_ice_runoff_xs_col[c]              = 0.0
        wf.qflx_glcice_dyn_water_flux_col[c]      = 0.0
    end

    # Column 2D irrigation layers
    nlevsoi = varpar.nlevsoi
    for c in bounds_col
        for j in 1:nlevsoi
            wf.qflx_gw_uncon_irrig_lyr_col[c, j] = 0.0
        end
    end

    # Soil/crop columns: zero drain, surf, latflow, discharge
    for c in bounds_col
        l = landunit_col[c]
        if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
            wf.qflx_drain_col[c]            = 0.0
            wf.qflx_surf_col[c]             = 0.0
            wf.qflx_latflow_in_col[c]       = 0.0
            wf.qflx_latflow_out_col[c]      = 0.0
            wf.volumetric_discharge_col[c]   = 0.0
        end
    end

    # Hillslope routing: zero streamflow for soil/crop landunits
    if use_hillslope_routing
        for l in bounds_lun
            if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
                wf.volumetric_streamflow_lun[l] = 0.0
            end
        end
    end

    return nothing
end

# ==========================================================================
# Stubs for infrastructure-dependent subroutines
# ==========================================================================

"""
    waterflux_init_history!(wf, bounds_col)

Register water flux fields for history file output.

Ported from `waterflux_type%InitHistory` in `WaterFluxType.F90`.
Requires history infrastructure — stub until that module is ported.
"""
function waterflux_init_history!(wf::WaterFluxData,
                                  bounds_col::UnitRange{Int})
    return nothing
end

"""
    waterflux_restart!(wf, bounds_col; flag="read")

Read/write water flux from/to restart file.

Ported from `waterflux_type%Restart` in `WaterFluxType.F90`.
Requires NetCDF/restart infrastructure — stub until that module is ported.
"""
function waterflux_restart!(wf::WaterFluxData,
                              bounds_col::UnitRange{Int};
                              flag::String = "read")
    return nothing
end
