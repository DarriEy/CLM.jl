# ==========================================================================
# Ported from: src/biogeophys/BalanceCheckMod.F90
# Water and energy balance check.
# ==========================================================================

# --------------------------------------------------------------------------
# Module-level constants and state
# --------------------------------------------------------------------------

const BALANCE_CHECK_SKIP_SIZE = 3600.0  # Time steps to skip the balance check at startup (sec)

"""
    BalanceCheckData

Module-level state for balance checking. Holds the number of timesteps
to skip at startup.

Ported from module-level `skip_steps` in `BalanceCheckMod.F90`.
"""
Base.@kwdef mutable struct BalanceCheckData
    skip_steps::Int = -999
end

# --------------------------------------------------------------------------
# BalanceCheckInit
# --------------------------------------------------------------------------

"""
    balance_check_init!(bc::BalanceCheckData, dtime::Float64)

Initialize balance check. Computes the number of timesteps to skip
based on `dtime` (model timestep in seconds).

Ported from `BalanceCheckInit` in `BalanceCheckMod.F90`.
"""
function balance_check_init!(bc::BalanceCheckData, dtime::Real)
    bc.skip_steps = max(2, round(Int, BALANCE_CHECK_SKIP_SIZE / dtime)) + 1
    return nothing
end

# --------------------------------------------------------------------------
# BalanceCheckClean
# --------------------------------------------------------------------------

"""
    balance_check_clean!(bc::BalanceCheckData)

Clean up BalanceCheck.

Ported from `BalanceCheckClean` in `BalanceCheckMod.F90`.
"""
function balance_check_clean!(bc::BalanceCheckData)
    bc.skip_steps = -999
    return nothing
end

# --------------------------------------------------------------------------
# GetBalanceCheckSkipSteps
# --------------------------------------------------------------------------

"""
    get_balance_check_skip_steps(bc::BalanceCheckData) -> Int

Get the number of steps to skip for the balance check.

Ported from `GetBalanceCheckSkipSteps` in `BalanceCheckMod.F90`.
"""
function get_balance_check_skip_steps(bc::BalanceCheckData)
    if bc.skip_steps > 0
        return bc.skip_steps
    else
        error("GetBalanceCheckSkipSteps called before BalanceCheckInit")
    end
end

# --------------------------------------------------------------------------
# WaterGridcellBalance (wrapper over bulk + tracers)
# --------------------------------------------------------------------------

"""
    water_gridcell_balance!(water::WaterData,
        lakestate::LakeStateData,
        col_data::ColumnData, lun_data::LandunitData, grc_data::GridcellData,
        mask_nolake::BitVector, mask_lake::BitVector,
        bounds_c::UnitRange{Int}, bounds_l::UnitRange{Int}, bounds_g::UnitRange{Int},
        flag::String; kwargs...)

Grid cell-level water balance for bulk water and each water tracer.

Ported from `WaterGridcellBalance` in `BalanceCheckMod.F90`.

Since `WaterAtm2LndData` and `WaterLnd2AtmData` types are not yet ported,
the dribbler amounts (`qflx_liq_dynbal_left_to_dribble`,
`qflx_ice_dynbal_left_to_dribble`) must be pre-computed by the caller
for the appropriate flag (begwb or endwb).
"""
function water_gridcell_balance!(
    water::WaterData,
    lakestate::LakeStateData,
    col_data::ColumnData,
    lun_data::LandunitData,
    grc_data::GridcellData,
    mask_nolake::BitVector,
    mask_lake::BitVector,
    bounds_c::UnitRange{Int},
    bounds_l::UnitRange{Int},
    bounds_g::UnitRange{Int},
    flag::String;
    use_aquifer_layer::Bool=false,
    use_hillslope_routing::Bool=false,
    qflx_liq_dynbal_left_to_dribble::Vector{<:Real}=Float64[],
    qflx_ice_dynbal_left_to_dribble::Vector{<:Real}=Float64[]
)
    for i in water.bulk_and_tracers_beg:water.bulk_and_tracers_end
        bt = water.bulk_and_tracers[i]
        water_gridcell_balance_single!(
            bt.waterstate,
            bt.waterbalance,
            bt.waterflux,
            bt.waterdiagnostic,
            lakestate,
            col_data, lun_data, grc_data,
            mask_nolake, mask_lake,
            bounds_c, bounds_l, bounds_g,
            flag;
            use_aquifer_layer=use_aquifer_layer,
            use_hillslope_routing=use_hillslope_routing,
            qflx_liq_dynbal_left_to_dribble=qflx_liq_dynbal_left_to_dribble,
            qflx_ice_dynbal_left_to_dribble=qflx_ice_dynbal_left_to_dribble,
            wa_reset_nonconservation_gain_col=isnothing(bt.waterbalance) ? Float64[] :
                bt.waterbalance.wa_reset_nonconservation_gain_col
        )
    end
    return nothing
end

# --------------------------------------------------------------------------
# WaterGridcellBalanceSingle
# --------------------------------------------------------------------------

"""
    water_gridcell_balance_single!(waterstate, waterbalance, waterflux,
        lakestate, col_data, lun_data, grc_data,
        mask_nolake, mask_lake,
        bounds_c, bounds_l, bounds_g, flag; kwargs...)

Grid cell-level water balance for bulk or a single tracer
at beginning or end of time step as specified by `flag` ("begwb" or "endwb").

Ported from `WaterGridcellBalanceSingle` in `BalanceCheckMod.F90`.

Uses `c2g_unity!` for column-to-gridcell aggregation (equivalent to
Fortran `c2g` with `c2l_scale_type='urbanf', l2g_scale_type='unity'`).
"""
function water_gridcell_balance_single!(
    waterstate::Union{WaterStateData, WaterStateBulkData, Nothing},
    waterbalance::Union{WaterBalanceData, Nothing},
    waterflux::Union{WaterFluxData, WaterFluxBulkData, Nothing},
    waterdiagnostic::Union{WaterDiagnosticBulkData, Nothing},
    lakestate::LakeStateData,
    col_data::ColumnData,
    lun_data::LandunitData,
    grc_data::GridcellData,
    mask_nolake::BitVector,
    mask_lake::BitVector,
    bounds_c::UnitRange{Int},
    bounds_l::UnitRange{Int},
    bounds_g::UnitRange{Int},
    flag::String;
    use_aquifer_layer::Bool=false,
    use_hillslope_routing::Bool=false,
    qflx_liq_dynbal_left_to_dribble::Vector{<:Real}=Float64[],
    qflx_ice_dynbal_left_to_dribble::Vector{<:Real}=Float64[],
    wa_reset_nonconservation_gain_col::Vector{<:Real}=Float64[]
)
    isnothing(waterbalance) && return nothing
    isnothing(waterstate) && return nothing

    begwb_grc = waterbalance.begwb_grc
    endwb_grc = waterbalance.endwb_grc

    # Temporary arrays
    nc = length(bounds_c) > 0 ? last(bounds_c) : 0
    ng = length(bounds_g) > 0 ? last(bounds_g) : 0
    FT = eltype(begwb_grc)
    wb_col = zeros(FT, nc)
    wb_grc = zeros(FT, ng)

    # Compute water mass for non-lake columns
    ws_raw = waterstate isa WaterStateBulkData ? waterstate.ws : waterstate
    compute_water_mass_non_lake_bc!(wb_col, waterstate, waterdiagnostic, mask_nolake, bounds_c, col_data)

    # Compute water mass for lake columns
    compute_water_mass_lake_bc!(wb_col, waterstate, lakestate, mask_lake, bounds_c, col_data)

    # Column-to-gridcell aggregation
    c2g_unity!(wb_grc, wb_col, col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)

    # Add landunit-level state (stream water volume), convert from m3 to kg/m2
    if use_hillslope_routing
        ws = waterstate isa WaterStateBulkData ? waterstate.ws : waterstate
        stream_water_volume_lun = ws.stream_water_volume_lun
        for l in bounds_l
            g = lun_data.gridcell[l]
            if g in bounds_g
                wb_grc[g] += stream_water_volume_lun[l] * 1.0e3 / (grc_data.area[g] * 1.0e6)
            end
        end
    end

    # Subtract dynbal dribbler amounts
    if !isempty(qflx_liq_dynbal_left_to_dribble)
        for g in bounds_g
            wb_grc[g] -= qflx_liq_dynbal_left_to_dribble[g] +
                          qflx_ice_dynbal_left_to_dribble[g]
        end
    end

    # Map wb_grc to beginning/ending water balance according to flag
    if flag == "begwb"
        for g in bounds_g
            begwb_grc[g] = wb_grc[g]
        end
    elseif flag == "endwb"
        # endwb_grc requires wa_reset_nonconservation_gain adjustment
        wa_reset_grc = zeros(FT, ng)
        if use_aquifer_layer && !isempty(wa_reset_nonconservation_gain_col)
            c2g_unity!(wa_reset_grc, wa_reset_nonconservation_gain_col,
                       col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)
        end
        for g in bounds_g
            endwb_grc[g] = wb_grc[g] - wa_reset_grc[g]
        end
    else
        error("Unknown flag '$flag' passed to water_gridcell_balance_single!. " *
              "Expecting either 'begwb' or 'endwb'.")
    end

    return nothing
end

# --------------------------------------------------------------------------
# compute_water_mass_lake_bc!
# --------------------------------------------------------------------------

"""
    compute_water_mass_lake_bc!(water_mass, waterstate, lakestate,
        mask_lake, bounds, col_data)

Compute total water mass for lake columns.
Delegates to `compute_water_mass_lake!` in `total_water_heat.jl`.

Ported from `ComputeWaterMassLake` in `TotalWaterAndHeatMod.F90`.
"""
function compute_water_mass_lake_bc!(
    water_mass::Vector{<:Real},
    waterstate::Union{WaterStateData, WaterStateBulkData},
    lakestate::LakeStateData,
    mask_lake::BitVector,
    bounds::UnitRange{Int},
    col_data::ColumnData
)
    ws = waterstate isa WaterStateBulkData ? waterstate.ws : waterstate
    compute_water_mass_lake!(mask_lake, col_data, ws, lakestate, false, water_mass)
    return nothing
end

# --------------------------------------------------------------------------
# compute_water_mass_non_lake_bc!
# --------------------------------------------------------------------------

"""
    compute_water_mass_non_lake_bc!(water_mass, waterstate, waterdiagnostic,
        mask_nolake, bounds, col_data)

Compute total water mass for non-lake columns.
Delegates to `compute_water_mass_non_lake!` in `total_water_heat.jl`.

Ported from `ComputeWaterMassNonLake` in `TotalWaterAndHeatMod.F90`.
"""
function compute_water_mass_non_lake_bc!(
    water_mass::Vector{<:Real},
    waterstate::Union{WaterStateData, WaterStateBulkData},
    waterdiagnostic::Union{WaterDiagnosticBulkData, Nothing},
    mask_nolake::BitVector,
    bounds::UnitRange{Int},
    col_data::ColumnData
)
    ws = waterstate isa WaterStateBulkData ? waterstate.ws : waterstate
    if waterdiagnostic !== nothing
        compute_water_mass_non_lake!(mask_nolake, col_data, ws, waterdiagnostic, false, water_mass)
    else
        # Tracer path: no waterdiagnostic available, compute without plant stored water
        nc = length(water_mass)
        nlevsno = varpar.nlevsno
        for c in eachindex(mask_nolake)
            mask_nolake[c] || continue
            water_mass[c] = ws.h2osno_no_layers_col[c] + ws.h2osfc_col[c]
            for j in (col_data.snl[c] + 1):0
                jj = j + nlevsno
                water_mass[c] += ws.h2osoi_liq_col[c, jj] + ws.h2osoi_ice_col[c, jj]
            end
            if col_data.hydrologically_active[c]
                water_mass[c] += (ws.wa_col[c] - ws.aquifer_water_baseline)
            end
            for j in 1:varpar.nlevgrnd
                jj = j + nlevsno
                water_mass[c] += ws.h2osoi_liq_col[c, jj] + ws.h2osoi_ice_col[c, jj]
            end
        end
    end
    return nothing
end

# --------------------------------------------------------------------------
# BeginWaterColumnBalance (wrapper over bulk + tracers)
# --------------------------------------------------------------------------

"""
    begin_water_column_balance!(water::WaterData,
        soilhydrology::SoilHydrologyData,
        lakestate::LakeStateData,
        col_data::ColumnData, lun_data::LandunitData,
        mask_nolake::BitVector, mask_lake::BitVector,
        bounds_c::UnitRange{Int};
        use_aquifer_layer::Bool)

Initialize column-level water balance at beginning of time step, for bulk
water and each water tracer.

Ported from `BeginWaterColumnBalance` in `BalanceCheckMod.F90`.
"""
function begin_water_column_balance!(
    water::WaterData,
    soilhydrology::SoilHydrologyData,
    lakestate::LakeStateData,
    col_data::ColumnData,
    lun_data::LandunitData,
    mask_nolake::BitVector,
    mask_lake::BitVector,
    bounds_c::UnitRange{Int};
    use_aquifer_layer::Bool=false
)
    for i in water.bulk_and_tracers_beg:water.bulk_and_tracers_end
        bt = water.bulk_and_tracers[i]
        begin_water_column_balance_single!(
            bt.waterstate,
            bt.waterbalance,
            bt.waterdiagnostic,
            soilhydrology,
            lakestate,
            col_data,
            lun_data,
            mask_nolake, mask_lake,
            bounds_c;
            use_aquifer_layer=use_aquifer_layer
        )
    end
    return nothing
end

# --------------------------------------------------------------------------
# BeginWaterColumnBalanceSingle
# --------------------------------------------------------------------------

"""
    begin_water_column_balance_single!(waterstate, waterbalance,
        soilhydrology, lakestate, col_data, lun_data,
        mask_nolake, mask_lake, bounds_c;
        use_aquifer_layer)

Initialize column-level water balance at beginning of time step, for bulk
or a single tracer.

Ported from `BeginWaterColumnBalanceSingle` in `BalanceCheckMod.F90`.
"""
function begin_water_column_balance_single!(
    waterstate::Union{WaterStateData, WaterStateBulkData, Nothing},
    waterbalance::Union{WaterBalanceData, Nothing},
    waterdiagnostic::Union{WaterDiagnosticBulkData, Nothing},
    soilhydrology::SoilHydrologyData,
    lakestate::LakeStateData,
    col_data::ColumnData,
    lun_data::LandunitData,
    mask_nolake::BitVector,
    mask_lake::BitVector,
    bounds_c::UnitRange{Int};
    use_aquifer_layer::Bool=false
)
    isnothing(waterbalance) && return nothing
    isnothing(waterstate) && return nothing

    zi = col_data.zi
    zwt = soilhydrology.zwt_col
    nlevsoi_val = varpar.nlevsoi
    joff_zi = varpar.nlevsno + 1

    ws = waterstate isa WaterStateBulkData ? waterstate.ws : waterstate
    aquifer_water_baseline = ws.aquifer_water_baseline
    wa = ws.wa_col

    wa_reset_nonconservation_gain = waterbalance.wa_reset_nonconservation_gain_col
    begwb = waterbalance.begwb_col
    h2osno_old = waterbalance.h2osno_old_col

    # Reset aquifer water to baseline under certain conditions
    if use_aquifer_layer
        for c in bounds_c
            mask_nolake[c] || continue
            if col_data.active[c]
                if is_hydrologically_active(col_data.itype[c], lun_data.itype[col_data.landunit[c]])
                    # zi is indexed on combined snow+soil interfaces.
                    if zwt[c] <= zi[c, joff_zi + nlevsoi_val]
                        wa_reset_nonconservation_gain[c] = aquifer_water_baseline - wa[c]
                        wa[c] = aquifer_water_baseline
                    else
                        wa_reset_nonconservation_gain[c] = 0.0
                    end
                end
            end
        end
    end

    # Compute water mass for non-lake columns → begwb
    compute_water_mass_non_lake_bc!(begwb, waterstate, waterdiagnostic, mask_nolake, bounds_c, col_data)

    # Compute water mass for lake columns → begwb
    compute_water_mass_lake_bc!(begwb, waterstate, lakestate, mask_lake, bounds_c, col_data)

    # Calculate total h2osno for snow balance tracking
    snl_col = col_data.snl
    waterstate_calculate_total_h2osno!(ws, mask_nolake, bounds_c, snl_col, h2osno_old)
    waterstate_calculate_total_h2osno!(ws, mask_lake, bounds_c, snl_col, h2osno_old)

    return nothing
end

# --------------------------------------------------------------------------
# BalanceCheck — main water & energy balance check
# --------------------------------------------------------------------------

# Warning/error thresholds (from Fortran parameters)
const H2O_WARNING_THRESH    = 1.0e-9
const ENERGY_WARNING_THRESH = 1.0e-7
const BALANCE_ERROR_THRESH  = 1.0e-5

"""
    balance_check!(bc::BalanceCheckData, ...)

Main water and energy balance check subroutine.

Accumulates numerical truncation errors of the water and energy balance
calculation. The error for energy balance:
  error = abs(Net radiation - change of internal energy - Sensible heat - Latent heat)
The error for water balance:
  error = abs(precipitation - change of water storage - evaporation - runoff)

Ported from `BalanceCheck` in `BalanceCheckMod.F90`.

Since `Atm2LndData`, `WaterAtm2LndData`, and `WaterLnd2AtmData` types are
not yet ported, the forcing and flux fields are passed as keyword arguments.
The column-level `qflx_ice_runoff_col` (from `WaterLnd2AtmData`) must also
be passed as a keyword argument.
"""
function balance_check!(
    bc::BalanceCheckData,
    waterflux::Union{WaterFluxData, WaterFluxBulkData},
    waterstate::Union{WaterStateData, WaterStateBulkData},
    waterbalance::WaterBalanceData,
    waterdiagnosticbulk::WaterDiagnosticBulkData,
    energyflux::EnergyFluxData,
    solarabs::SolarAbsorbedData,
    canopystate::CanopyStateData,
    surfalb::SurfaceAlbedoData,
    col_data::ColumnData,
    lun_data::LandunitData,
    pat_data::PatchData,
    grc_data::GridcellData,
    mask_allc::BitVector,
    bounds_c::UnitRange{Int},
    bounds_p::UnitRange{Int},
    bounds_g::UnitRange{Int},
    nstep::Int,
    DAnstep::Int,
    dtime::Real;
    # --- Atmospheric forcing (from Atm2Lnd / WaterAtm2Lnd) ---
    forc_rain_col::Vector{<:Real}=Float64[],
    forc_snow_col::Vector{<:Real}=Float64[],
    forc_rain_grc::Vector{<:Real}=Float64[],
    forc_snow_grc::Vector{<:Real}=Float64[],
    forc_solad_col::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0),
    forc_solai_grc::Matrix{<:Real}=Matrix{Float64}(undef, 0, 0),
    forc_lwrad_col::Vector{<:Real}=Float64[],
    forc_flood_grc::Vector{<:Real}=Float64[],
    # --- Fluxes from WaterLnd2Atm (not yet ported as types) ---
    qflx_ice_runoff_col::Vector{<:Real}=Float64[],
    qflx_evap_tot_grc::Vector{<:Real}=Float64[],
    qflx_surf_grc::Vector{<:Real}=Float64[],
    qflx_qrgwl_grc::Vector{<:Real}=Float64[],
    qflx_drain_grc::Vector{<:Real}=Float64[],
    qflx_drain_perched_grc::Vector{<:Real}=Float64[],
    qflx_ice_runoff_grc::Vector{<:Real}=Float64[],
    qflx_sfc_irrig_grc::Vector{<:Real}=Float64[],
    qflx_streamflow_grc::Vector{<:Real}=Float64[],
    # --- Control flags ---
    use_fates_planthydro::Bool=false,
    use_soil_moisture_streams::Bool=false,
    use_hillslope_routing::Bool=false,
    for_testing_zero_dynbal_fluxes::Bool=false
)
    # Extract water flux fields
    wf = waterflux isa WaterFluxBulkData ? waterflux.wf : waterflux
    ws = waterstate isa WaterStateBulkData ? waterstate.ws : waterstate

    qflx_liq_grnd_col       = wf.qflx_liq_grnd_col
    qflx_snow_grnd_col      = wf.qflx_snow_grnd_col
    qflx_snwcp_liq          = wf.qflx_snwcp_liq_col
    qflx_snwcp_ice          = wf.qflx_snwcp_ice_col
    qflx_snwcp_discarded_liq_col = wf.qflx_snwcp_discarded_liq_col
    qflx_snwcp_discarded_ice_col = wf.qflx_snwcp_discarded_ice_col
    qflx_evap_tot_col       = wf.qflx_evap_tot_col
    qflx_soliddew_to_top_layer    = wf.qflx_soliddew_to_top_layer_col
    qflx_solidevap_from_top_layer = wf.qflx_solidevap_from_top_layer_col
    qflx_liqevap_from_top_layer   = wf.qflx_liqevap_from_top_layer_col
    qflx_liqdew_to_top_layer      = wf.qflx_liqdew_to_top_layer_col
    qflx_snow_h2osfc        = wf.qflx_snow_h2osfc_col
    qflx_h2osfc_to_ice      = wf.qflx_h2osfc_to_ice_col
    qflx_drain_perched_col  = wf.qflx_drain_perched_col
    qflx_flood_col          = wf.qflx_floodc_col
    qflx_snow_drain         = wf.qflx_snow_drain_col
    qflx_surf_col           = wf.qflx_surf_col
    qflx_qrgwl_col          = wf.qflx_qrgwl_col
    qflx_drain_col          = wf.qflx_drain_col
    qflx_sfc_irrig_col      = wf.qflx_sfc_irrig_col
    qflx_glcice_dyn_water_flux_col = wf.qflx_glcice_dyn_water_flux_col
    qflx_sl_top_soil        = wf.qflx_sl_top_soil_col

    # Water balance fields
    h2osno_old      = waterbalance.h2osno_old_col
    begwb_grc       = waterbalance.begwb_grc
    endwb_grc       = waterbalance.endwb_grc
    begwb_col       = waterbalance.begwb_col
    endwb_col       = waterbalance.endwb_col
    errh2o_col      = waterbalance.errh2o_col
    errh2osno       = waterbalance.errh2osno_col
    snow_sources    = waterbalance.snow_sources_col
    snow_sinks      = waterbalance.snow_sinks_col

    # Diagnostic bulk fields
    frac_sno_eff    = waterdiagnosticbulk.frac_sno_eff_col
    frac_sno        = waterdiagnosticbulk.frac_sno_col
    snow_depth      = waterdiagnosticbulk.snow_depth_col
    qflx_prec_grnd  = waterdiagnosticbulk.qflx_prec_grnd_col

    # Energy flux fields
    dhsdt_canopy    = energyflux.dhsdt_canopy_patch
    eflx_lwrad_out  = energyflux.eflx_lwrad_out_patch
    eflx_lwrad_net  = energyflux.eflx_lwrad_net_patch
    eflx_sh_tot     = energyflux.eflx_sh_tot_patch
    eflx_lh_tot     = energyflux.eflx_lh_tot_patch
    eflx_soil_grnd  = energyflux.eflx_soil_grnd_patch
    eflx_wasteheat_p = energyflux.eflx_wasteheat_patch
    eflx_ventilation_p = energyflux.eflx_ventilation_patch
    eflx_heat_from_ac_p = energyflux.eflx_heat_from_ac_patch
    eflx_traffic_p  = energyflux.eflx_traffic_patch
    errsoi_col      = energyflux.errsoi_col
    errsol          = energyflux.errsol_patch
    errseb          = energyflux.errseb_patch
    errlon          = energyflux.errlon_patch
    netrad          = energyflux.netrad_patch

    # Solar absorbed fields
    sabg_soil       = solarabs.sabg_soil_patch
    sabg_snow       = solarabs.sabg_snow_patch
    sabg_chk        = solarabs.sabg_chk_patch
    fsa             = solarabs.fsa_patch
    fsr             = solarabs.fsr_patch
    sabv            = solarabs.sabv_patch
    sabg            = solarabs.sabg_patch

    # Canopy state
    elai            = canopystate.elai_patch
    esai            = canopystate.esai_patch

    # Surface albedo
    fabd            = surfalb.fabd_patch
    fabi            = surfalb.fabi_patch
    albd            = surfalb.albd_patch
    albi            = surfalb.albi_patch
    ftdd            = surfalb.ftdd_patch
    ftid            = surfalb.ftid_patch
    ftii            = surfalb.ftii_patch

    skip_steps = bc.skip_steps

    # =====================================================================
    # Determine column level incoming snow and rain
    # Assume no incident precipitation on urban wall columns
    # =====================================================================

    nc = length(bounds_c) > 0 ? last(bounds_c) : 0
    FT_bc = eltype(errh2o_col)
    forc_rain_c = zeros(FT_bc, nc)
    forc_snow_c = zeros(FT_bc, nc)

    for c in bounds_c
        if col_data.itype[c] == ICOL_SUNWALL || col_data.itype[c] == ICOL_SHADEWALL
            forc_rain_c[c] = 0.0
            forc_snow_c[c] = 0.0
        else
            forc_rain_c[c] = forc_rain_col[c]
            forc_snow_c[c] = forc_snow_col[c]
        end
    end

    # =====================================================================
    # Water balance check at the column level
    # =====================================================================

    for c in bounds_c
        if col_data.active[c]
            errh2o_col[c] = endwb_col[c] - begwb_col[c] -
                (forc_rain_c[c] +
                 forc_snow_c[c] +
                 qflx_flood_col[c] +
                 qflx_sfc_irrig_col[c] +
                 qflx_glcice_dyn_water_flux_col[c] -
                 qflx_evap_tot_col[c] -
                 qflx_surf_col[c] -
                 qflx_qrgwl_col[c] -
                 qflx_drain_col[c] -
                 qflx_drain_perched_col[c] -
                 qflx_ice_runoff_col[c] -
                 qflx_snwcp_discarded_liq_col[c] -
                 qflx_snwcp_discarded_ice_col[c]) * dtime
        else
            errh2o_col[c] = 0.0
        end
    end

    errh2o_max_val = maximum(abs.(errh2o_col[bounds_c]))

    if errh2o_max_val > H2O_WARNING_THRESH
        indexc = bounds_c[argmax(abs.(errh2o_col[bounds_c]))]
        @warn "column-level water balance error" nstep indexc errh2o=errh2o_col[indexc]

        if errh2o_max_val > BALANCE_ERROR_THRESH && DAnstep > skip_steps
            if _is_ad_type(FT_bc)
                @warn "BalanceCheck: column water balance error exceeded threshold (AD mode, continuing)" maxlog=1
            else
                @error "Stopping: errh2o > $(BALANCE_ERROR_THRESH) mm" nstep indexc errh2o=errh2o_col[indexc] forc_rain=forc_rain_c[indexc]*dtime forc_snow=forc_snow_c[indexc]*dtime endwb=endwb_col[indexc] begwb=begwb_col[indexc] qflx_evap_tot=qflx_evap_tot_col[indexc]*dtime qflx_surf=qflx_surf_col[indexc]*dtime qflx_drain=qflx_drain_col[indexc]*dtime
                error("BalanceCheck: column water balance error exceeded threshold at c=$indexc")
            end
        end
    end

    # =====================================================================
    # Water balance check at the grid cell level
    # =====================================================================

    ng = length(bounds_g) > 0 ? last(bounds_g) : 0
    errh2o_grc = zeros(FT_bc, ng)
    qflx_glcice_dyn_water_flux_grc_arr = zeros(FT_bc, ng)
    qflx_snwcp_discarded_liq_grc_arr = zeros(FT_bc, ng)
    qflx_snwcp_discarded_ice_grc_arr = zeros(FT_bc, ng)

    c2g_unity!(qflx_glcice_dyn_water_flux_grc_arr, qflx_glcice_dyn_water_flux_col,
               col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)
    c2g_unity!(qflx_snwcp_discarded_liq_grc_arr, qflx_snwcp_discarded_liq_col,
               col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)
    c2g_unity!(qflx_snwcp_discarded_ice_grc_arr, qflx_snwcp_discarded_ice_col,
               col_data.gridcell, col_data.wtgcell, bounds_c, bounds_g)

    for g in bounds_g
        errh2o_grc[g] = endwb_grc[g] - begwb_grc[g] -
            (forc_rain_grc[g] +
             forc_snow_grc[g] +
             forc_flood_grc[g] +
             qflx_sfc_irrig_grc[g] +
             qflx_glcice_dyn_water_flux_grc_arr[g] -
             qflx_evap_tot_grc[g] -
             qflx_surf_grc[g] -
             qflx_qrgwl_grc[g] -
             qflx_drain_grc[g] -
             qflx_drain_perched_grc[g] -
             qflx_ice_runoff_grc[g] -
             qflx_snwcp_discarded_liq_grc_arr[g] -
             qflx_snwcp_discarded_ice_grc_arr[g]) * dtime
    end

    # Add landunit level flux (streamflow) for hillslope routing
    if use_hillslope_routing
        for g in bounds_g
            errh2o_grc[g] += qflx_streamflow_grc[g] * dtime
        end
    end

    errh2o_grc_max_val = maximum(abs.(errh2o_grc[bounds_g]))

    # BUG(rgk, 2021-04-13, ESCOMP/CTSM#1314) Temporarily bypassing gridcell-level check
    # with use_fates_planthydro
    if errh2o_grc_max_val > H2O_WARNING_THRESH && !use_fates_planthydro
        indexg = bounds_g[argmax(abs.(errh2o_grc[bounds_g]))]
        @warn "grid cell-level water balance error" nstep indexg errh2o_grc=errh2o_grc[indexg]

        if errh2o_grc_max_val > BALANCE_ERROR_THRESH && DAnstep > skip_steps &&
           !use_soil_moisture_streams && !for_testing_zero_dynbal_fluxes
            if _is_ad_type(FT_bc)
                @warn "BalanceCheck: gridcell water balance error exceeded threshold (AD mode, continuing)" maxlog=1
            else
                @error "Stopping: errh2o_grc > $(BALANCE_ERROR_THRESH) mm" nstep indexg errh2o_grc=errh2o_grc[indexg] forc_rain=forc_rain_grc[indexg]*dtime forc_snow=forc_snow_grc[indexg]*dtime endwb_grc=endwb_grc[indexg] begwb_grc=begwb_grc[indexg]
                error("BalanceCheck: gridcell water balance error exceeded threshold at g=$indexg")
            end
        end
    end

    # =====================================================================
    # Snow balance check at the column level
    # =====================================================================

    h2osno_total = zeros(FT_bc, nc)
    waterstate_calculate_total_h2osno!(ws, mask_allc, bounds_c, col_data.snl, h2osno_total)

    for c in bounds_c
        if col_data.active[c]
            l = col_data.landunit[c]

            if col_data.snl[c] < 0
                snow_sources[c] = qflx_prec_grnd[c] + qflx_soliddew_to_top_layer[c] +
                    qflx_liqdew_to_top_layer[c]
                snow_sinks[c] = qflx_solidevap_from_top_layer[c] +
                    qflx_liqevap_from_top_layer[c] +
                    qflx_snow_drain[c] + qflx_snwcp_ice[c] + qflx_snwcp_liq[c] +
                    qflx_snwcp_discarded_ice_col[c] + qflx_snwcp_discarded_liq_col[c] +
                    qflx_sl_top_soil[c]

                if lun_data.itype[l] == ISTDLAK
                    snow_sources[c] = qflx_snow_grnd_col[c] +
                        frac_sno_eff[c] * (qflx_liq_grnd_col[c] +
                        qflx_soliddew_to_top_layer[c] + qflx_liqdew_to_top_layer[c])
                    snow_sinks[c] = frac_sno_eff[c] * (qflx_solidevap_from_top_layer[c] +
                        qflx_liqevap_from_top_layer[c]) + qflx_snwcp_ice[c] + qflx_snwcp_liq[c] +
                        qflx_snwcp_discarded_ice_col[c] + qflx_snwcp_discarded_liq_col[c] +
                        qflx_snow_drain[c] + qflx_sl_top_soil[c]
                end

                if col_data.itype[c] == ICOL_ROAD_PERV || lun_data.itype[l] == ISTSOIL ||
                   lun_data.itype[l] == ISTCROP || lun_data.itype[l] == ISTWET ||
                   lun_data.itype[l] == ISTICE
                    snow_sources[c] = (qflx_snow_grnd_col[c] - qflx_snow_h2osfc[c]) +
                        frac_sno_eff[c] * (qflx_liq_grnd_col[c] +
                        qflx_soliddew_to_top_layer[c] + qflx_liqdew_to_top_layer[c]) +
                        qflx_h2osfc_to_ice[c]
                    snow_sinks[c] = frac_sno_eff[c] * (qflx_solidevap_from_top_layer[c] +
                        qflx_liqevap_from_top_layer[c]) + qflx_snwcp_ice[c] + qflx_snwcp_liq[c] +
                        qflx_snwcp_discarded_ice_col[c] + qflx_snwcp_discarded_liq_col[c] +
                        qflx_snow_drain[c] + qflx_sl_top_soil[c]
                end

                errh2osno[c] = (h2osno_total[c] - h2osno_old[c]) -
                    (snow_sources[c] - snow_sinks[c]) * dtime
            else
                snow_sources[c] = 0.0
                snow_sinks[c] = 0.0
                errh2osno[c] = 0.0
            end
        else
            errh2osno[c] = 0.0
        end
    end

    errh2osno_max_val = maximum(abs.(errh2osno[bounds_c]))

    if errh2osno_max_val > H2O_WARNING_THRESH
        indexc = bounds_c[argmax(abs.(errh2osno[bounds_c]))]
        l_idx = col_data.landunit[indexc]
        @warn "snow balance error" nstep indexc col_itype=col_data.itype[indexc] lun_itype=lun_data.itype[l_idx] errh2osno=errh2osno[indexc]

        if errh2osno_max_val > BALANCE_ERROR_THRESH && DAnstep > skip_steps
            if _is_ad_type(FT_bc)
                @warn "BalanceCheck: snow balance error exceeded threshold (AD mode, continuing)" maxlog=1
            else
                @error "Stopping: errh2osno > $(BALANCE_ERROR_THRESH) mm" nstep indexc errh2osno=errh2osno[indexc] snl=col_data.snl[indexc] snow_depth=snow_depth[indexc] h2osno=h2osno_total[indexc] h2osno_old=h2osno_old[indexc] snow_sources=snow_sources[indexc]*dtime snow_sinks=snow_sinks[indexc]*dtime
                error("BalanceCheck: snow balance error exceeded threshold at c=$indexc")
            end
        end
    end

    # =====================================================================
    # Energy balance checks
    # =====================================================================

    for p in bounds_p
        if pat_data.active[p]
            c = pat_data.column[p]
            l = pat_data.landunit[p]
            g = pat_data.gridcell[p]

            # Solar radiation energy balance
            if !lun_data.urbpoi[l]
                errsol[p] = fsa[p] + fsr[p] -
                    (forc_solad_col[c, 1] + forc_solad_col[c, 2] +
                     forc_solai_grc[g, 1] + forc_solai_grc[g, 2])
            else
                errsol[p] = SPVAL
            end

            # Longwave radiation energy balance
            if !lun_data.urbpoi[l]
                errlon[p] = eflx_lwrad_out[p] - eflx_lwrad_net[p] - forc_lwrad_col[c]
            else
                errlon[p] = SPVAL
            end

            # Surface energy balance
            if !lun_data.urbpoi[l]
                errseb[p] = sabv[p] + sabg_chk[p] + forc_lwrad_col[c] - eflx_lwrad_out[p] -
                    eflx_sh_tot[p] - eflx_lh_tot[p] - eflx_soil_grnd[p] - dhsdt_canopy[p]
            else
                errseb[p] = sabv[p] + sabg[p] -
                    eflx_lwrad_net[p] -
                    eflx_sh_tot[p] - eflx_lh_tot[p] - eflx_soil_grnd[p] +
                    eflx_wasteheat_p[p] + eflx_heat_from_ac_p[p] + eflx_traffic_p[p] +
                    eflx_ventilation_p[p]
            end

            # Net radiation
            netrad[p] = fsa[p] - eflx_lwrad_net[p]
        else
            errsol[p] = 0.0
            errlon[p] = 0.0
            errseb[p] = 0.0
        end
    end

    # Solar radiation energy balance check
    errsol_vals = [errsol[p] != SPVAL ? abs(errsol[p]) : 0.0 for p in bounds_p]
    errsol_max_val = isempty(errsol_vals) ? 0.0 : maximum(errsol_vals)

    if errsol_max_val > ENERGY_WARNING_THRESH && DAnstep > skip_steps
        indexp = bounds_p[argmax(errsol_vals)]
        @warn "solar radiation balance error (W/m2)" nstep errsol=errsol[indexp]

        if errsol_max_val > BALANCE_ERROR_THRESH
            if _is_ad_type(eltype(errsol))
                @warn "BalanceCheck: solar radiation balance error exceeded threshold (AD mode, continuing)" maxlog=1
            else
                @error "Stopping: errsol > $(BALANCE_ERROR_THRESH) W/m2" nstep indexp errsol=errsol[indexp] fsa=fsa[indexp] fsr=fsr[indexp]
                error("BalanceCheck: solar radiation balance error exceeded threshold at p=$indexp")
            end
        end
    end

    # Longwave radiation energy balance check
    errlon_vals = [errlon[p] != SPVAL ? abs(errlon[p]) : 0.0 for p in bounds_p]
    errlon_max_val = isempty(errlon_vals) ? 0.0 : maximum(errlon_vals)

    if errlon_max_val > ENERGY_WARNING_THRESH && DAnstep > skip_steps
        indexp = bounds_p[argmax(errlon_vals)]
        @warn "longwave energy balance error (W/m2)" nstep indexp errlon=errlon[indexp]

        if errlon_max_val > BALANCE_ERROR_THRESH
            if _is_ad_type(eltype(errlon))
                @warn "BalanceCheck: longwave energy balance error exceeded threshold (AD mode, continuing)" maxlog=1
            else
                @error "Stopping: errlon > $(BALANCE_ERROR_THRESH) W/m2" nstep indexp errlon=errlon[indexp]
                error("BalanceCheck: longwave energy balance error exceeded threshold at p=$indexp")
            end
        end
    end

    # Surface energy balance check
    errseb_max_val = isempty(bounds_p) ? 0.0 : maximum(abs.(errseb[bounds_p]))

    if errseb_max_val > ENERGY_WARNING_THRESH && DAnstep > skip_steps
        indexp = bounds_p[argmax(abs.(errseb[bounds_p]))]
        @warn "surface flux energy balance error (W/m2)" nstep errseb=errseb[indexp]

        if errseb_max_val > BALANCE_ERROR_THRESH
            if _is_ad_type(eltype(errseb))
                @warn "BalanceCheck: surface energy balance error exceeded threshold (AD mode, continuing)" maxlog=1
            else
                @error "Stopping: errseb > $(BALANCE_ERROR_THRESH) W/m2" nstep indexp errseb=errseb[indexp] sabv=sabv[indexp] sabg=sabg[indexp] eflx_lwrad_net=eflx_lwrad_net[indexp] eflx_sh_tot=eflx_sh_tot[indexp] eflx_lh_tot=eflx_lh_tot[indexp] eflx_soil_grnd=eflx_soil_grnd[indexp] dhsdt_canopy=dhsdt_canopy[indexp]
                error("BalanceCheck: surface energy balance error exceeded threshold at p=$indexp")
            end
        end
    end

    # Soil energy balance check
    errsoi_vals = [col_data.active[c] ? abs(errsoi_col[c]) : 0.0 for c in bounds_c]
    errsoi_col_max_val = isempty(errsoi_vals) ? 0.0 : maximum(errsoi_vals)

    if errsoi_col_max_val > 1.0e-5
        indexc = bounds_c[argmax(errsoi_vals)]
        @warn "soil balance error (W/m2)" nstep errsoi_col=errsoi_col[indexc]

        if errsoi_col_max_val > 1.0e-4 && DAnstep > skip_steps
            if _is_ad_type(eltype(errsoi_col))
                @warn "BalanceCheck: soil energy balance error exceeded threshold (AD mode, continuing)" maxlog=1
            else
                @error "Stopping: errsoi_col > 1.0e-4 W/m2" nstep indexc errsoi_col=errsoi_col[indexc]
                error("BalanceCheck: soil energy balance error exceeded threshold at c=$indexc")
            end
        end
    end

    return nothing
end
