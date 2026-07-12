# ==========================================================================
# Ported from: src/biogeophys/HydrologyDrainageMod.F90 (~244 lines)
# Calculates soil/snow hydrology with drainage (subsurface runoff).
#
# Public functions:
#   compute_wetland_ice_hydrology!   — Set fluxes for wetland/ice landunits
#   compute_urban_drainage_fluxes!   — Set drainage fluxes for urban non-pervious
#   compute_total_runoff!            — Compute total runoff
#   hydrology_drainage!              — Main orchestrator
# (end-of-timestep water mass comes from compute_water_mass_non_lake! in
#  biogeophys/total_water_heat.jl, exactly as the Fortran calls
#  TotalWaterAndHeatMod::ComputeWaterMassNonLake)
# (adjust_runoff_terms! for glacier SMB lives in glacier_surface_mass_balance.jl)
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
# ---- compute_wetland_ice_hydrology! : per-column wetland/ice + urban fluxes ----
# Each column is fully independent (no loop-carried deps / cross-column coupling).
# One thread per column, mask-gated. dtime is passed as a working-eltype scalar
# and SPVAL is eltype-converted so no Float64 reaches a Float32-only backend
# (Metal); byte-identical on a Float64 CPU run. Bare 0.0 stores are left as-is.
@kernel function _hyddr_wetland_ice_kernel!(qflx_drain, qflx_drain_perched, qflx_surf,
        qflx_infl, qflx_qrgwl, qflx_latflow_out, qflx_rsub_sat,
        qflx_ice_runoff_snwcp, @Const(mask), @Const(col_landunit), @Const(col_gridcell),
        @Const(col_itype), @Const(lun_itype), @Const(lun_urbpoi),
        @Const(qflx_evap_tot), @Const(qflx_snwcp_ice), @Const(qflx_snwcp_discarded_ice),
        @Const(qflx_snwcp_discarded_liq), @Const(forc_rain), @Const(forc_snow),
        @Const(qflx_floodg), @Const(begwb), @Const(endwb),
        istwet::Int, istice::Int, icol_road_perv::Int, spval, dtime)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(qflx_drain)
        l = col_landunit[c]
        g = col_gridcell[c]

        if lun_itype[l] == istwet || lun_itype[l] == istice
            qflx_latflow_out[c]   = 0.0
            qflx_drain[c]         = 0.0
            qflx_drain_perched[c] = 0.0
            qflx_surf[c]          = 0.0
            qflx_infl[c]          = 0.0
            qflx_qrgwl[c] = forc_rain[c] + forc_snow[c] + qflx_floodg[g] -
                             qflx_evap_tot[c] - qflx_snwcp_ice[c] -
                             qflx_snwcp_discarded_ice[c] - qflx_snwcp_discarded_liq[c] -
                             (endwb[c] - begwb[c]) / T(dtime)
        elseif lun_urbpoi[l] && col_itype[c] != icol_road_perv
            qflx_drain_perched[c] = 0.0
            qflx_rsub_sat[c]      = T(spval)
            qflx_infl[c]          = 0.0
        end

        qflx_ice_runoff_snwcp[c] = qflx_snwcp_ice[c]
    end
end

function compute_wetland_ice_hydrology!(
    qflx_drain::AbstractVector{<:Real},
    qflx_drain_perched::AbstractVector{<:Real},
    qflx_surf::AbstractVector{<:Real},
    qflx_infl::AbstractVector{<:Real},
    qflx_qrgwl::AbstractVector{<:Real},
    qflx_latflow_out::AbstractVector{<:Real},
    qflx_rsub_sat::AbstractVector{<:Real},
    qflx_evap_tot::AbstractVector{<:Real},
    qflx_snwcp_ice::AbstractVector{<:Real},
    qflx_snwcp_discarded_ice::AbstractVector{<:Real},
    qflx_snwcp_discarded_liq::AbstractVector{<:Real},
    qflx_ice_runoff_snwcp::AbstractVector{<:Real},
    forc_rain::AbstractVector{<:Real},
    forc_snow::AbstractVector{<:Real},
    qflx_floodg::AbstractVector{<:Real},
    begwb::AbstractVector{<:Real},
    endwb::AbstractVector{<:Real},
    col_landunit::AbstractVector{<:Integer},
    col_gridcell::AbstractVector{<:Integer},
    col_itype::AbstractVector{<:Integer},
    lun_itype::AbstractVector{<:Integer},
    lun_urbpoi::AbstractVector{Bool},
    mask_nolake::AbstractVector{Bool},
    bounds::UnitRange{Int},
    dtime::Real
)
    T = eltype(qflx_drain)
    _launch!(_hyddr_wetland_ice_kernel!, qflx_drain, qflx_drain_perched, qflx_surf,
             qflx_infl, qflx_qrgwl, qflx_latflow_out, qflx_rsub_sat,
             qflx_ice_runoff_snwcp, mask_nolake, col_landunit, col_gridcell,
             col_itype, lun_itype, lun_urbpoi, qflx_evap_tot, qflx_snwcp_ice,
             qflx_snwcp_discarded_ice, qflx_snwcp_discarded_liq, forc_rain, forc_snow,
             qflx_floodg, begwb, endwb, ISTWET, ISTICE, ICOL_ROAD_PERV, T(SPVAL), T(dtime);
             ndrange = length(mask_nolake))

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
# ---- compute_total_runoff! : per-column total + urban/rural split ----
# Each column is fully independent. One thread per column, mask-gated. No Float64
# literals in arithmetic; byte-identical on a Float64 CPU run.
@kernel function _hyddr_total_runoff_kernel!(qflx_runoff, qflx_runoff_u, qflx_runoff_r,
        @Const(qflx_drain), @Const(qflx_surf), @Const(qflx_qrgwl),
        @Const(qflx_drain_perched), @Const(mask), @Const(col_landunit),
        @Const(lun_itype), @Const(lun_urbpoi), istsoil::Int, istcrop::Int)
    c = @index(Global)
    @inbounds if mask[c]
        l = col_landunit[c]
        qflx_runoff[c] = qflx_drain[c] + qflx_surf[c] + qflx_qrgwl[c] + qflx_drain_perched[c]

        if lun_urbpoi[l]
            qflx_runoff_u[c] = qflx_runoff[c]
        elseif lun_itype[l] == istsoil || lun_itype[l] == istcrop
            qflx_runoff_r[c] = qflx_runoff[c]
        end
    end
end

function compute_total_runoff!(
    qflx_runoff::AbstractVector{<:Real},
    qflx_runoff_u::AbstractVector{<:Real},
    qflx_runoff_r::AbstractVector{<:Real},
    qflx_drain::AbstractVector{<:Real},
    qflx_surf::AbstractVector{<:Real},
    qflx_qrgwl::AbstractVector{<:Real},
    qflx_drain_perched::AbstractVector{<:Real},
    col_landunit::AbstractVector{<:Integer},
    lun_itype::AbstractVector{<:Integer},
    lun_urbpoi::AbstractVector{Bool},
    mask_nolake::AbstractVector{Bool},
    bounds::UnitRange{Int}
)
    _launch!(_hyddr_total_runoff_kernel!, qflx_runoff, qflx_runoff_u, qflx_runoff_r,
             qflx_drain, qflx_surf, qflx_qrgwl, qflx_drain_perched, mask_nolake,
             col_landunit, lun_itype, lun_urbpoi, ISTSOIL, ISTCROP;
             ndrange = length(mask_nolake))

    return nothing
end

# =========================================================================
# adjust_runoff_terms! for glacier SMB is now ported in
# biogeophys/glacier_surface_mass_balance.jl (alongside handle_ice_melt! and
# compute_surface_mass_balance!). The orchestrator below calls it directly.
# =========================================================================

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
    mask_nolake::AbstractVector{Bool},
    mask_hydrology::AbstractVector{Bool},
    mask_urban::AbstractVector{Bool},
    mask_do_smb::AbstractVector{Bool},
    forc_rain::AbstractVector{<:Real},
    forc_snow::AbstractVector{<:Real},
    qflx_floodg::AbstractVector{<:Real},
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
        # gridcell ids: pull to host so the maximum() stays scalar-index-free on device.
        ng = isempty(bounds) ? 0 : maximum(Array(col.gridcell)[bounds])
        # tdepth/tdepthmax are only read on the (CPU-only) hillslope branch, so plain
        # host scratch is fine. grc_area is read inside the subsurface_lateral_flow!
        # device kernel, so it must live on the same backend as the state arrays.
        tdepth_grc = zeros(FT, ng)
        tdepthmax_grc = zeros(FT, ng)
        grc_area = fill!(similar(waterstatebulk.ws.h2osoi_liq_col, FT, ng), one(FT))

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
    # Fortran: call ComputeWaterMassNonLake(bounds, num_nolakec, filter_nolakec,
    #              waterstatebulk_inst, waterdiagnosticbulk_inst,
    #              subtract_dynbal_baselines = .false., water_mass = endwb(...))
    # The real implementation lives in biogeophys/total_water_heat.jl. Canopy water
    # (p2c of liqcan/snocan) is added by the driver's balance seam, consistent with
    # balance_check.jl's begwb/endwb path; endwb here is used by the wetland/ice
    # branch below, whose columns carry no canopy water.
    compute_water_mass_non_lake!(
        mask_nolake, col,
        waterstatebulk.ws, waterdiagbulk,
        false,
        waterbalancebulk.endwb_col)

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
    # glc_dyn_runoff_routing (dynamic ice-sheet coupling fraction, per gridcell)
    # is not yet ported (no glc2lnd coupling) → standalone-CLM case = 0 everywhere.
    FT_g = eltype(wf.qflx_qrgwl_col)
    ng_g = isempty(bounds) ? 0 : maximum(Array(col.gridcell)[bounds])
    glc_dyn_runoff_routing_grc =
        fill!(similar(wf.qflx_qrgwl_col, FT_g, ng_g), zero(FT_g))
    adjust_runoff_terms!(
        wf.qflx_qrgwl_col,
        wf.qflx_ice_runoff_snwcp_col,
        wf.qflx_glcice_frz_col,
        wf.qflx_glcice_melt_col,
        glc_dyn_runoff_routing_grc,
        col.gridcell,
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
