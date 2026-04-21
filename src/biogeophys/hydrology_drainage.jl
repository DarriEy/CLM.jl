# ==========================================================================
# Ported from: src/biogeophys/HydrologyDrainageMod.F90 (~244 lines)
# Calculates soil/snow hydrology with drainage (subsurface runoff).
#
# Public functions:
#   compute_wetland_ice_hydrology!   — Set fluxes for wetland/ice landunits
#   compute_urban_drainage_fluxes!   — Set drainage fluxes for urban non-pervious
#   compute_total_runoff!            — Compute total runoff
#   compute_water_mass_non_lake!     — Compute water mass (stub)
#   adjust_runoff_terms!             — Adjust runoff terms for glacier SMB (stub)
#   hydrology_drainage!              — Main orchestrator
# ==========================================================================

# =========================================================================
# compute_wetland_ice_hydrology!
# =========================================================================

"""
    compute_wetland_ice_hydrology!(
        qflx_drain, qflx_drain_perched, qflx_surf, qflx_infl,
        qflx_qrgwl, qflx_latflow_out, qflx_rsub_sat,
        qflx_evap_tot, qflx_snwcp_ice,
        qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq,
        qflx_ice_runoff_snwcp,
        forc_rain, forc_snow, qflx_floodg,
        begwb, endwb,
        col_landunit, col_gridcell, col_itype,
        lun_itype, lun_urbpoi,
        mask_nolake, bounds, dtime)

Set hydrological fluxes for wetland and land ice columns. For these
landunit types, drainage, infiltration, and surface runoff are zeroed,
and `qflx_qrgwl` absorbs the water balance residual.

Also handles urban columns that are not pervious roads: zeroes
perched drainage, sets rsub_sat to SPVAL, and zeroes infiltration.

Assigns `qflx_ice_runoff_snwcp` for all non-lake columns.

Ported from inline code in `HydrologyDrainage` in `HydrologyDrainageMod.F90`.
"""
function compute_wetland_ice_hydrology!(
    qflx_drain::Vector{<:Real},
    qflx_drain_perched::Vector{<:Real},
    qflx_surf::Vector{<:Real},
    qflx_infl::Vector{<:Real},
    qflx_qrgwl::Vector{<:Real},
    qflx_latflow_out::Vector{<:Real},
    qflx_rsub_sat::Vector{<:Real},
    qflx_evap_tot::Vector{<:Real},
    qflx_snwcp_ice::Vector{<:Real},
    qflx_snwcp_discarded_ice::Vector{<:Real},
    qflx_snwcp_discarded_liq::Vector{<:Real},
    qflx_ice_runoff_snwcp::Vector{<:Real},
    forc_rain::Vector{<:Real},
    forc_snow::Vector{<:Real},
    qflx_floodg::Vector{<:Real},
    begwb::Vector{<:Real},
    endwb::Vector{<:Real},
    col_landunit::Vector{Int},
    col_gridcell::Vector{Int},
    col_itype::Vector{Int},
    lun_itype::Vector{Int},
    lun_urbpoi::AbstractVector{Bool},
    mask_nolake::BitVector,
    bounds::UnitRange{Int},
    dtime::Real
)
    for c in bounds
        mask_nolake[c] || continue
        l = col_landunit[c]
        g = col_gridcell[c]

        if lun_itype[l] == ISTWET || lun_itype[l] == ISTICE
            qflx_latflow_out[c]   = 0.0
            qflx_drain[c]         = 0.0
            qflx_drain_perched[c] = 0.0
            qflx_surf[c]          = 0.0
            qflx_infl[c]          = 0.0
            qflx_qrgwl[c] = forc_rain[c] + forc_snow[c] + qflx_floodg[g] -
                             qflx_evap_tot[c] - qflx_snwcp_ice[c] -
                             qflx_snwcp_discarded_ice[c] - qflx_snwcp_discarded_liq[c] -
                             (endwb[c] - begwb[c]) / dtime
        elseif lun_urbpoi[l] && col_itype[c] != ICOL_ROAD_PERV
            qflx_drain_perched[c] = 0.0
            qflx_rsub_sat[c]      = SPVAL
            qflx_infl[c]          = 0.0
        end

        qflx_ice_runoff_snwcp[c] = qflx_snwcp_ice[c]
    end

    return nothing
end

# =========================================================================
# compute_total_runoff!
# =========================================================================

"""
    compute_total_runoff!(
        qflx_runoff, qflx_runoff_u, qflx_runoff_r,
        qflx_drain, qflx_surf, qflx_qrgwl, qflx_drain_perched,
        col_landunit, lun_itype, lun_urbpoi,
        mask_nolake, bounds)

Compute total runoff as the sum of drainage, surface runoff, glacier/wetland
runoff, and perched drainage. Assigns urban and rural runoff separately.

Ported from inline code in `HydrologyDrainage` in `HydrologyDrainageMod.F90`.
"""
function compute_total_runoff!(
    qflx_runoff::Vector{<:Real},
    qflx_runoff_u::Vector{<:Real},
    qflx_runoff_r::Vector{<:Real},
    qflx_drain::Vector{<:Real},
    qflx_surf::Vector{<:Real},
    qflx_qrgwl::Vector{<:Real},
    qflx_drain_perched::Vector{<:Real},
    col_landunit::Vector{Int},
    lun_itype::Vector{Int},
    lun_urbpoi::AbstractVector{Bool},
    mask_nolake::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask_nolake[c] || continue
        l = col_landunit[c]

        qflx_runoff[c] = qflx_drain[c] + qflx_surf[c] + qflx_qrgwl[c] + qflx_drain_perched[c]

        if lun_urbpoi[l]
            qflx_runoff_u[c] = qflx_runoff[c]
        elseif lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
            qflx_runoff_r[c] = qflx_runoff[c]
        end
    end

    return nothing
end

# =========================================================================
# compute_water_mass_non_lake! (stub)
# =========================================================================

"""
    compute_water_mass_non_lake!(endwb, waterstatebulk, waterdiagbulk,
        mask_nolake, bounds)

Compute total water mass for non-lake columns.

Ported from `ComputeWaterMassNonLake` in `TotalWaterAndHeatMod.F90`.
Stub: TotalWaterAndHeatMod is not yet ported. When ported, this function
should compute total water mass (liquid + ice + snow + surface water +
aquifer) for each non-lake column and store in endwb.
"""
function compute_water_mass_non_lake!(
    endwb::Vector{<:Real},
    waterstatebulk::WaterStateBulkData,
    waterdiagbulk::WaterDiagnosticBulkData,
    mask_nolake::BitVector,
    bounds::UnitRange{Int}
)
    # Stub: sets endwb to 0.0 for all non-lake columns.
    # When TotalWaterAndHeatMod is ported, replace with actual computation.
    for c in bounds
        mask_nolake[c] || continue
        endwb[c] = 0.0
    end

    return nothing
end

# =========================================================================
# adjust_runoff_terms! (stub)
# =========================================================================

"""
    adjust_runoff_terms!(qflx_qrgwl, qflx_ice_runoff_snwcp,
        waterfluxbulk, mask_do_smb, bounds)

Adjust runoff terms for glacier surface mass balance.

Ported from `glacier_smb_inst%AdjustRunoffTerms` in
`GlacierSurfaceMassBalanceMod.F90`.
Stub: GlacierSurfaceMassBalanceMod is not yet ported. When ported, this
function should adjust qflx_qrgwl and qflx_ice_runoff_snwcp based on
glacier SMB calculations.
"""
function adjust_runoff_terms!(
    qflx_qrgwl::Vector{<:Real},
    qflx_ice_runoff_snwcp::Vector{<:Real},
    waterfluxbulk::WaterFluxBulkData,
    mask_do_smb::BitVector,
    bounds::UnitRange{Int}
)
    # Stub: no-op until GlacierSurfaceMassBalanceMod is ported
    return nothing
end

# =========================================================================
# hydrology_drainage! — Main orchestrator
# =========================================================================

"""
    hydrology_drainage!(
        temperature, soilhydrology, soilstate,
        waterstatebulk, waterdiagbulk, waterbalancebulk, waterfluxbulk,
        col, lun,
        mask_nolake, mask_hydrology, mask_urban, mask_do_smb,
        forc_rain, forc_snow, qflx_floodg,
        bounds, dtime, nlevsno, nlevsoi, nlevgrnd, nlevurb;
        use_vichydro, use_aquifer_layer, use_hillslope_routing)

Main subroutine to calculate soil/snow hydrology with drainage
(subsurface runoff).

Calls drainage or lateral flow routines, updates volumetric soil water,
computes water balance, determines wetland/ice hydrology, adjusts
glacier SMB runoff terms, and computes total runoff.

Ported from `HydrologyDrainage` in `HydrologyDrainageMod.F90`.
"""
function hydrology_drainage!(
    temperature::TemperatureData,
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    waterstatebulk::WaterStateBulkData,
    waterdiagbulk::WaterDiagnosticBulkData,
    waterbalancebulk::WaterBalanceData,
    waterfluxbulk::WaterFluxBulkData,
    col::ColumnData,
    lun::LandunitData,
    mask_nolake::BitVector,
    mask_hydrology::BitVector,
    mask_urban::BitVector,
    mask_do_smb::BitVector,
    forc_rain::Vector{<:Real},
    forc_snow::Vector{<:Real},
    qflx_floodg::Vector{<:Real},
    bounds::UnitRange{Int},
    dtime::Real,
    nlevsno::Int,
    nlevsoi::Int,
    nlevgrnd::Int,
    nlevurb::Int;
    use_vichydro::Bool = false,
    use_aquifer_layer::Bool = true,
    use_hillslope_routing::Bool = false
)
    wf = waterfluxbulk.wf

    # --- VIC hydrology mapping ---
    if use_vichydro
        clm_vic_map!(
            soilhydrology,
            waterstatebulk.ws,
            col.dz, col.zi, col.z,
            mask_hydrology, bounds,
            nlevsoi, nlevgrnd,
            nlevsoi, nlevsoi + 1)  # nlayer = nlevsoi, nlayert = nlevsoi+1
    end

    # --- Drainage or lateral flow ---
    if use_aquifer_layer
        drainage!(
            temperature.t_soisno_col,
            soilhydrology, soilstate,
            waterstatebulk.ws, waterfluxbulk,
            col.dz, col.z, col.zi, col.snl,
            col.itype, col.landunit, col.topo_slope,
            lun.urbpoi,
            mask_hydrology, mask_urban, bounds,
            nlevsoi, dtime;
            nlevsno=nlevsno,
            use_vichydro=use_vichydro)
    else
        # Perched/subsurface lateral flow expects gridcell-indexed vectors.
        # In single-point mode, use safe defaults when no external routing data exist.
        FT = eltype(temperature.t_soisno_col)
        ng = isempty(bounds) ? 0 : maximum(col.gridcell[bounds])
        tdepth_grc = zeros(FT, ng)
        tdepthmax_grc = zeros(FT, ng)
        grc_area = ones(FT, ng)

        perched_lateral_flow!(
            soilhydrology, soilstate,
            waterstatebulk.ws, waterfluxbulk,
            col, lun,
            tdepth_grc, tdepthmax_grc,
            mask_hydrology, bounds,
            nlevsoi, dtime;
            use_hillslope_routing=use_hillslope_routing)

        subsurface_lateral_flow!(
            soilhydrology, soilstate,
            waterstatebulk.ws, waterfluxbulk,
            col, lun,
            tdepth_grc, tdepthmax_grc, grc_area,
            mask_hydrology, mask_urban, bounds,
            nlevsoi, dtime;
            use_hillslope_routing=use_hillslope_routing)

        # Hillslope routing (stub: not yet ported)
        # if use_hillslope_routing
        #     hillslope_stream_outflow!(...)
        #     hillslope_update_stream_water!(...)
        # end
    end

    # --- Update volumetric soil water ---
    update_h2osoi_vol!(
        waterstatebulk.ws.h2osoi_vol_col,
        waterstatebulk.ws.h2osoi_liq_col,
        waterstatebulk.ws.h2osoi_ice_col,
        col.dz, col.itype,
        mask_nolake, mask_urban, bounds,
        nlevgrnd, nlevurb, nlevsno)

    # --- Compute end-of-timestep water mass ---
    compute_water_mass_non_lake!(
        waterbalancebulk.endwb_col,
        waterstatebulk, waterdiagbulk,
        mask_nolake, bounds)

    # --- Wetland/ice hydrology and urban drainage fluxes ---
    compute_wetland_ice_hydrology!(
        wf.qflx_drain_col,
        wf.qflx_drain_perched_col,
        wf.qflx_surf_col,
        wf.qflx_infl_col,
        wf.qflx_qrgwl_col,
        wf.qflx_latflow_out_col,
        wf.qflx_rsub_sat_col,
        wf.qflx_evap_tot_col,
        wf.qflx_snwcp_ice_col,
        wf.qflx_snwcp_discarded_ice_col,
        wf.qflx_snwcp_discarded_liq_col,
        wf.qflx_ice_runoff_snwcp_col,
        forc_rain, forc_snow, qflx_floodg,
        waterbalancebulk.begwb_col,
        waterbalancebulk.endwb_col,
        col.landunit, col.gridcell, col.itype,
        lun.itype, lun.urbpoi,
        mask_nolake, bounds, dtime)

    # --- Adjust runoff terms for glacier SMB ---
    adjust_runoff_terms!(
        wf.qflx_qrgwl_col,
        wf.qflx_ice_runoff_snwcp_col,
        waterfluxbulk,
        mask_do_smb, bounds)

    # --- Compute total runoff ---
    compute_total_runoff!(
        wf.qflx_runoff_col,
        wf.qflx_runoff_u_col,
        wf.qflx_runoff_r_col,
        wf.qflx_drain_col,
        wf.qflx_surf_col,
        wf.qflx_qrgwl_col,
        wf.qflx_drain_perched_col,
        col.landunit,
        lun.itype, lun.urbpoi,
        mask_nolake, bounds)

    return nothing
end
