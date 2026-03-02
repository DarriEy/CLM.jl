# ==========================================================================
# Ported from: src/biogeophys/SoilHydrologyMod.F90 (~2910 lines)
# Calculate soil hydrology
#
# Public functions:
#   set_soil_water_fractions!       — Set diagnostic water/ice fractions per layer
#   set_floodc!                     — Apply gridcell flood water flux to non-lake columns
#   set_qflx_inputs!               — Set flux of water into soil from top
#   route_infiltration_excess!      — Route infiltration excess runoff flux
#   infiltration!                   — Calculate total infiltration
#   total_surface_runoff!           — Calculate total surface runoff
#   update_urban_ponding!           — Update ponding on urban surfaces
#   water_table!                    — Calculate water table (pre-drainage)
#   drainage!                       — Calculate subsurface drainage
#   clm_vic_map!                    — Map CLM layers to VIC layers
#   perched_water_table!            — Calculate perched water table
#   perched_lateral_flow!           — Calculate lateral flow from perched saturated zone
#   theta_based_water_table!        — Calculate water table from soil moisture state
#   subsurface_lateral_flow!        — Calculate subsurface lateral flow
#   renew_condensation!             — Misc. condensation/sublimation corrections
#   calc_irrig_withdrawals!         — Calculate irrigation withdrawals from groundwater
#   withdraw_groundwater_irrigation! — Remove groundwater irrigation
# ==========================================================================

# ---- Module-level parameters (from params_type in Fortran) ----

"""
    SoilHydrologyParams

Parameters read from parameter file for soil hydrology.

Ported from `params_type` in `SoilHydrologyMod.F90`.
"""
Base.@kwdef mutable struct SoilHydrologyParams
    aq_sp_yield_min::Float64 = 0.01           # Minimum aquifer specific yield (unitless)
    n_baseflow::Float64 = 1.0                 # Drainage power law exponent (unitless)
    perched_baseflow_scalar::Float64 = 5.0e-4 # Scalar multiplier for perched base flow rate (kg/m2/s)
    e_ice::Float64 = 6.0                      # Soil ice impedance factor (unitless)
end

const soilhydrology_params = SoilHydrologyParams()

# ---- Module-level constants ----
const BASEFLOW_SCALAR = Ref(1.0e-2)                  # baseflow scalar (namelist adjustable)
const TOLERANCE_SOILHYDRO = 1.0e-12                   # tolerance for sublimation check

# Head gradient methods
const HEAD_GRADIENT_KINEMATIC = 0
const HEAD_GRADIENT_DARCY     = 1
# Transmissivity methods
const TRANSMISSIVITY_UNIFORM  = 0
const TRANSMISSIVITY_LAYERSUM = 1

# Module-level method settings (defaults)
const HEAD_GRADIENT_METHOD  = Ref(HEAD_GRADIENT_DARCY)
const TRANSMISSIVITY_METHOD = Ref(TRANSMISSIVITY_LAYERSUM)

"""
    SoilHydrologyConfig

Configuration for hillslope hydrology settings.
"""
Base.@kwdef mutable struct SoilHydrologyConfig
    head_gradient_method::Int  = HEAD_GRADIENT_DARCY
    transmissivity_method::Int = TRANSMISSIVITY_LAYERSUM
    baseflow_scalar::Float64   = 1.0e-2
end

"""
    init_soil_hydrology_config(; kwargs...) -> SoilHydrologyConfig

Create a soil hydrology configuration with defaults.
"""
function init_soil_hydrology_config(;
    head_gradient_method::Int = HEAD_GRADIENT_DARCY,
    transmissivity_method::Int = TRANSMISSIVITY_LAYERSUM,
    baseflow_scalar::Float64 = 1.0e-2
)
    cfg = SoilHydrologyConfig(
        head_gradient_method = head_gradient_method,
        transmissivity_method = transmissivity_method,
        baseflow_scalar = baseflow_scalar
    )
    # Update module-level refs
    HEAD_GRADIENT_METHOD[] = cfg.head_gradient_method
    TRANSMISSIVITY_METHOD[] = cfg.transmissivity_method
    BASEFLOW_SCALAR[] = cfg.baseflow_scalar
    return cfg
end

# =========================================================================
# SetSoilWaterFractions
# =========================================================================

"""
    set_soil_water_fractions!(soilhydrology, soilstate, waterstatebulk,
                               col_dz, mask_hydrology, bounds, nlevsoi)

Set diagnostic variables related to the fraction of water and ice in each layer.

Ported from `SetSoilWaterFractions` in `SoilHydrologyMod.F90`.
"""
function set_soil_water_fractions!(
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    waterstatebulk_ws::WaterStateData,
    col_dz::Matrix{Float64},
    mask_hydrology::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int
)
    watsat       = soilstate.watsat_col
    eff_porosity = soilstate.eff_porosity_col
    h2osoi_liq   = waterstatebulk_ws.h2osoi_liq_col
    h2osoi_ice   = waterstatebulk_ws.h2osoi_ice_col
    excess_ice   = waterstatebulk_ws.excess_ice_col
    icefrac      = soilhydrology.icefrac_col

    for j in 1:nlevsoi
        for c in bounds
            mask_hydrology[c] || continue

            dz_ext = col_dz[c, j] + excess_ice[c, j] / DENICE
            vol_ice_val = min(watsat[c, j], (h2osoi_ice[c, j] + excess_ice[c, j]) / (dz_ext * DENICE))
            eff_porosity[c, j] = max(0.01, watsat[c, j] - vol_ice_val)
            icefrac[c, j] = min(1.0, vol_ice_val / watsat[c, j])
        end
    end

    return nothing
end

# =========================================================================
# SetFloodc
# =========================================================================

"""
    set_floodc!(qflx_floodc, qflx_floodg, col_gridcell, col_itype,
                mask_nolake, bounds)

Apply gridcell flood water flux to non-lake columns.

Ported from `SetFloodc` in `SoilHydrologyMod.F90`.
"""
function set_floodc!(
    qflx_floodc::Vector{Float64},
    qflx_floodg::Vector{Float64},
    col_gridcell::Vector{Int},
    col_itype::Vector{Int},
    mask_nolake::BitVector,
    bounds::UnitRange{Int}
)
    for c in bounds
        mask_nolake[c] || continue
        g = col_gridcell[c]
        if col_itype[c] != ICOL_SUNWALL && col_itype[c] != ICOL_SHADEWALL
            qflx_floodc[c] = qflx_floodg[g]
        else
            qflx_floodc[c] = 0.0
        end
    end

    return nothing
end

# =========================================================================
# SetQflxInputs
# =========================================================================

"""
    set_qflx_inputs!(waterfluxbulk, waterdiagbulk, col_snl,
                      mask_hydrology, bounds)

Set various input fluxes of water.

Ported from `SetQflxInputs` in `SoilHydrologyMod.F90`.
"""
function set_qflx_inputs!(
    waterfluxbulk::WaterFluxBulkData,
    waterdiagbulk::WaterDiagnosticBulkData,
    col_snl::Vector{Int},
    mask_hydrology::BitVector,
    bounds::UnitRange{Int}
)
    wf = waterfluxbulk.wf

    qflx_top_soil           = wf.qflx_top_soil_col
    qflx_in_soil            = waterfluxbulk.qflx_in_soil_col
    qflx_top_soil_to_h2osfc = waterfluxbulk.qflx_top_soil_to_h2osfc_col
    qflx_rain_plus_snomelt  = wf.qflx_rain_plus_snomelt_col
    qflx_snow_h2osfc        = wf.qflx_snow_h2osfc_col
    qflx_floodc             = wf.qflx_floodc_col
    qflx_ev_soil            = waterfluxbulk.qflx_ev_soil_col
    qflx_liqevap_from_top_layer = wf.qflx_liqevap_from_top_layer_col
    qflx_ev_h2osfc          = waterfluxbulk.qflx_ev_h2osfc_col
    qflx_sat_excess_surf    = waterfluxbulk.qflx_sat_excess_surf_col

    frac_sno    = waterdiagbulk.frac_sno_eff_col
    frac_h2osfc = waterdiagbulk.frac_h2osfc_col

    for c in bounds
        mask_hydrology[c] || continue

        qflx_top_soil[c] = qflx_rain_plus_snomelt[c] + qflx_snow_h2osfc[c] + qflx_floodc[c]

        # Partition surface inputs between soil and h2osfc
        if col_snl[c] >= 0
            fsno = 0.0
            qflx_evap = qflx_liqevap_from_top_layer[c]
        else
            fsno = frac_sno[c]
            qflx_evap = qflx_ev_soil[c]
        end

        qflx_in_soil[c] = (1.0 - frac_h2osfc[c]) * (qflx_top_soil[c] - qflx_sat_excess_surf[c])
        qflx_top_soil_to_h2osfc[c] = frac_h2osfc[c] * (qflx_top_soil[c] - qflx_sat_excess_surf[c])

        # remove evaporation (snow treated in SnowHydrology)
        qflx_in_soil[c] = qflx_in_soil[c] - (1.0 - fsno - frac_h2osfc[c]) * qflx_evap
        qflx_top_soil_to_h2osfc[c] = qflx_top_soil_to_h2osfc[c] - frac_h2osfc[c] * qflx_ev_h2osfc[c]
    end

    return nothing
end

# =========================================================================
# RouteInfiltrationExcess
# =========================================================================

"""
    route_infiltration_excess!(waterfluxbulk, soilhydrology,
                                col_landunit, lun_itype,
                                mask_hydrology, bounds)

Route the infiltration excess runoff flux.

Ported from `RouteInfiltrationExcess` in `SoilHydrologyMod.F90`.
"""
function route_infiltration_excess!(
    waterfluxbulk::WaterFluxBulkData,
    soilhydrology::SoilHydrologyData,
    col_landunit::Vector{Int},
    lun_itype::Vector{Int},
    mask_hydrology::BitVector,
    bounds::UnitRange{Int}
)
    wf = waterfluxbulk.wf

    qflx_in_soil_limited    = waterfluxbulk.qflx_in_soil_limited_col
    qflx_in_h2osfc          = waterfluxbulk.qflx_in_h2osfc_col
    qflx_infl_excess_surf   = waterfluxbulk.qflx_infl_excess_surf_col
    qflx_in_soil            = waterfluxbulk.qflx_in_soil_col
    qflx_top_soil_to_h2osfc = waterfluxbulk.qflx_top_soil_to_h2osfc_col
    qflx_infl_excess        = waterfluxbulk.qflx_infl_excess_col

    h2osfcflag = soilhydrology.h2osfcflag

    for c in bounds
        mask_hydrology[c] || continue
        l = col_landunit[c]
        if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
            qflx_in_soil_limited[c] = qflx_in_soil[c] - qflx_infl_excess[c]
            if h2osfcflag != 0
                qflx_in_h2osfc[c] = qflx_top_soil_to_h2osfc[c] + qflx_infl_excess[c]
                qflx_infl_excess_surf[c] = 0.0
            else
                qflx_in_h2osfc[c] = qflx_top_soil_to_h2osfc[c]
                qflx_infl_excess_surf[c] = qflx_infl_excess[c]
            end
        else
            qflx_in_soil_limited[c] = qflx_in_soil[c]
            qflx_in_h2osfc[c] = 0.0
            qflx_infl_excess_surf[c] = 0.0
        end
    end

    return nothing
end

# =========================================================================
# Infiltration
# =========================================================================

"""
    infiltration!(waterfluxbulk, mask_hydrology, bounds)

Calculate total infiltration.

Ported from `Infiltration` in `SoilHydrologyMod.F90`.
"""
function infiltration!(
    waterfluxbulk::WaterFluxBulkData,
    mask_hydrology::BitVector,
    bounds::UnitRange{Int}
)
    wf = waterfluxbulk.wf

    qflx_infl            = wf.qflx_infl_col
    qflx_in_soil_limited = waterfluxbulk.qflx_in_soil_limited_col
    qflx_h2osfc_drain    = waterfluxbulk.qflx_h2osfc_drain_col

    for c in bounds
        mask_hydrology[c] || continue
        qflx_infl[c] = qflx_in_soil_limited[c] + qflx_h2osfc_drain[c]
    end

    return nothing
end

# =========================================================================
# TotalSurfaceRunoff
# =========================================================================

"""
    total_surface_runoff!(waterfluxbulk, soilhydrology, waterstatebulk_ws,
                           col_snl, col_itype, col_landunit, lun_urbpoi,
                           mask_hydrology, mask_urban, bounds, dtime)

Calculate total surface runoff.

Ported from `TotalSurfaceRunoff` in `SoilHydrologyMod.F90`.
"""
function total_surface_runoff!(
    waterfluxbulk::WaterFluxBulkData,
    soilhydrology::SoilHydrologyData,
    waterstatebulk_ws::WaterStateData,
    col_snl::Vector{Int},
    col_itype::Vector{Int},
    col_landunit::Vector{Int},
    lun_urbpoi::Vector{Bool},
    mask_hydrology::BitVector,
    mask_urban::BitVector,
    bounds::UnitRange{Int},
    dtime::Float64
)
    wf = waterfluxbulk.wf

    qflx_surf                = wf.qflx_surf_col
    qflx_infl_excess_surf    = waterfluxbulk.qflx_infl_excess_surf_col
    qflx_h2osfc_surf         = waterfluxbulk.qflx_h2osfc_surf_col
    qflx_rain_plus_snomelt   = wf.qflx_rain_plus_snomelt_col
    qflx_liqevap_from_top_layer = wf.qflx_liqevap_from_top_layer_col
    qflx_floodc              = wf.qflx_floodc_col
    qflx_sat_excess_surf     = waterfluxbulk.qflx_sat_excess_surf_col

    xs_urban    = soilhydrology.xs_urban_col
    h2osoi_liq  = waterstatebulk_ws.h2osoi_liq_col

    # Set qflx_surf for hydrologically-active columns
    for c in bounds
        mask_hydrology[c] || continue
        qflx_surf[c] = qflx_sat_excess_surf[c] + qflx_infl_excess_surf[c] + qflx_h2osfc_surf[c]
    end

    # Set qflx_surf for non-hydrologically-active urban columns
    for c in bounds
        mask_urban[c] || continue
        if col_itype[c] == ICOL_ROOF || col_itype[c] == ICOL_ROAD_IMPERV
            if col_snl[c] < 0
                qflx_surf[c] = max(0.0, qflx_rain_plus_snomelt[c])
            else
                xs_urban[c] = max(0.0,
                    h2osoi_liq[c, 1] / dtime + qflx_rain_plus_snomelt[c] -
                    qflx_liqevap_from_top_layer[c] - PONDMX_URBAN / dtime)
                qflx_surf[c] = xs_urban[c]
            end
            # send flood water flux to runoff for all urban columns
            qflx_surf[c] = qflx_surf[c] + qflx_floodc[c]
        elseif col_itype[c] == ICOL_SUNWALL || col_itype[c] == ICOL_SHADEWALL
            qflx_surf[c] = 0.0
        end
    end

    return nothing
end

# =========================================================================
# UpdateUrbanPonding
# =========================================================================

"""
    update_urban_ponding!(waterstatebulk_ws, soilhydrology, waterfluxbulk,
                           col_snl, col_itype, mask_urban, bounds, dtime)

Update the state variable representing ponding on urban surfaces.

Ported from `UpdateUrbanPonding` in `SoilHydrologyMod.F90`.
"""
function update_urban_ponding!(
    waterstatebulk_ws::WaterStateData,
    soilhydrology::SoilHydrologyData,
    waterfluxbulk::WaterFluxBulkData,
    col_snl::Vector{Int},
    col_itype::Vector{Int},
    mask_urban::BitVector,
    bounds::UnitRange{Int},
    dtime::Float64
)
    wf = waterfluxbulk.wf

    h2osoi_liq  = waterstatebulk_ws.h2osoi_liq_col
    xs_urban    = soilhydrology.xs_urban_col
    qflx_rain_plus_snomelt = wf.qflx_rain_plus_snomelt_col
    qflx_liqevap_from_top_layer = wf.qflx_liqevap_from_top_layer_col

    for c in bounds
        mask_urban[c] || continue

        if col_itype[c] == ICOL_ROOF || col_itype[c] == ICOL_ROAD_IMPERV
            if col_snl[c] >= 0
                if xs_urban[c] > 0.0
                    h2osoi_liq[c, 1] = PONDMX_URBAN
                else
                    h2osoi_liq[c, 1] = max(0.0, h2osoi_liq[c, 1] +
                        (qflx_rain_plus_snomelt[c] - qflx_liqevap_from_top_layer[c]) * dtime)
                end
            end
        end
    end

    return nothing
end

# =========================================================================
# WaterTable
# =========================================================================

"""
    water_table!(soilhydrology, soilstate, temperature, waterstatebulk_ws,
                  waterfluxbulk, col_dz, col_z, col_zi,
                  mask_hydrology, bounds, nlevsoi, dtime)

Calculate water table, considering aquifer recharge but no drainage.

Ported from `WaterTable` in `SoilHydrologyMod.F90`.
"""
function water_table!(
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    temperature_t_soisno::Matrix{Float64},
    waterstatebulk_ws::WaterStateData,
    waterfluxbulk::WaterFluxBulkData,
    col_dz::Matrix{Float64},
    col_z::Matrix{Float64},
    col_zi::Matrix{Float64},
    mask_hydrology::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Float64
)
    wf = waterfluxbulk.wf
    params = soilhydrology_params

    t_soisno   = temperature_t_soisno
    h2osfc     = waterstatebulk_ws.h2osfc_col
    h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col
    h2osoi_vol = waterstatebulk_ws.h2osoi_vol_col
    wa         = waterstatebulk_ws.wa_col

    bsw        = soilstate.bsw_col
    hksat      = soilstate.hksat_col
    sucsat     = soilstate.sucsat_col
    watsat     = soilstate.watsat_col
    eff_porosity = soilstate.eff_porosity_col

    zwt         = soilhydrology.zwt_col
    zwt_perched = soilhydrology.zwt_perched_col
    frost_table = soilhydrology.frost_table_col
    qcharge     = soilhydrology.qcharge_col

    qflx_drain         = wf.qflx_drain_col
    qflx_drain_perched = wf.qflx_drain_perched_col
    qflx_rsub_sat      = wf.qflx_rsub_sat_col

    nc = length(bounds)

    # Local arrays
    jwt = Vector{Int}(undef, last(bounds))

    # Convert layer thicknesses from m to mm (in-place temp)
    dzmm = Matrix{Float64}(undef, last(bounds), nlevsoi)
    for j in 1:nlevsoi
        for c in bounds
            mask_hydrology[c] || continue
            dzmm[c, j] = col_dz[c, j] * 1.0e3
        end
    end

    # Initialize drainage fluxes
    for c in bounds
        mask_hydrology[c] || continue
        qflx_drain[c]         = 0.0
        qflx_rsub_sat[c]      = 0.0
        qflx_drain_perched[c] = 0.0
    end

    # Find jwt: index of soil layer right above water table
    for c in bounds
        mask_hydrology[c] || continue
        jwt[c] = nlevsoi
        for j in 1:nlevsoi
            if zwt[c] <= col_zi[c, j]
                jwt[c] = j - 1
                break
            end
        end
    end

    # ========================== QCHARGE ====================================
    # Water table changes due to qcharge
    for c in bounds
        mask_hydrology[c] || continue

        # analytical expression for aquifer specific yield
        rous = watsat[c, nlevsoi] *
            (1.0 - (1.0 + 1.0e3 * zwt[c] / sucsat[c, nlevsoi])^(-1.0 / bsw[c, nlevsoi]))
        rous = max(rous, params.aq_sp_yield_min)

        if jwt[c] == nlevsoi
            # water table below soil column
            wa[c]  = wa[c] + qcharge[c] * dtime
            zwt[c] = zwt[c] - (qcharge[c] * dtime) / 1000.0 / rous
        else
            # water table within soil layers
            qcharge_tot = qcharge[c] * dtime
            if qcharge_tot > 0.0  # rising water table
                for j in (jwt[c]+1):-1:1
                    s_y = watsat[c, j] *
                        (1.0 - (1.0 + 1.0e3 * zwt[c] / sucsat[c, j])^(-1.0 / bsw[c, j]))
                    s_y = max(s_y, params.aq_sp_yield_min)

                    qcharge_layer = min(qcharge_tot, s_y * (zwt[c] - col_zi[c, j-1]) * 1.0e3)
                    qcharge_layer = max(qcharge_layer, 0.0)

                    if s_y > 0.0
                        zwt[c] = zwt[c] - qcharge_layer / s_y / 1000.0
                    end

                    qcharge_tot = qcharge_tot - qcharge_layer
                    if qcharge_tot <= 0.0
                        break
                    end
                end
            else  # deepening water table
                for j in (jwt[c]+1):nlevsoi
                    s_y = watsat[c, j] *
                        (1.0 - (1.0 + 1.0e3 * zwt[c] / sucsat[c, j])^(-1.0 / bsw[c, j]))
                    s_y = max(s_y, params.aq_sp_yield_min)

                    qcharge_layer = max(qcharge_tot, -(s_y * (col_zi[c, j] - zwt[c]) * 1.0e3))
                    qcharge_layer = min(qcharge_layer, 0.0)
                    qcharge_tot = qcharge_tot - qcharge_layer
                    if qcharge_tot >= 0.0
                        zwt[c] = zwt[c] - qcharge_layer / s_y / 1000.0
                        break
                    else
                        zwt[c] = col_zi[c, j]
                    end
                end
                if qcharge_tot > 0.0
                    zwt[c] = zwt[c] - qcharge_tot / 1000.0 / rous
                end
            end

            # recompute jwt
            jwt[c] = nlevsoi
            for j in 1:nlevsoi
                if zwt[c] <= col_zi[c, j]
                    jwt[c] = j - 1
                    break
                end
            end
        end
    end

    # ========================== BASEFLOW ====================================
    # perched water table code
    for c in bounds
        mask_hydrology[c] || continue

        # define frost table as first frozen layer with unfrozen layer above it
        if t_soisno[c, 1] > TFRZ
            k_frz = nlevsoi
        else
            k_frz = 1
        end

        for k in 2:nlevsoi
            if t_soisno[c, k-1] > TFRZ && t_soisno[c, k] <= TFRZ
                k_frz = k
                break
            end
        end

        frost_table[c] = col_z[c, k_frz]

        # initialize perched water table to frost table
        zwt_perched[c] = frost_table[c]

        # ======= water table above frost table ========
        if zwt[c] < frost_table[c] && t_soisno[c, k_frz] <= TFRZ
            # do nothing - handled in Drainage
        else
            # ======= water table below frost table ========
            sat_lev = 0.9

            k_perch = 1
            for k in k_frz:-1:1
                h2osoi_vol[c, k] = h2osoi_liq[c, k] / (col_dz[c, k] * DENH2O) +
                    h2osoi_ice[c, k] / (col_dz[c, k] * DENICE)

                if h2osoi_vol[c, k] / watsat[c, k] <= sat_lev
                    k_perch = k
                    break
                end
            end

            # if frost_table = nlevsoi, only compute perched water table if frozen
            if t_soisno[c, k_frz] > TFRZ
                k_perch = k_frz
            end

            # if perched water table exists
            if k_frz > k_perch
                s1 = (h2osoi_liq[c, k_perch] / (col_dz[c, k_perch] * DENH2O) +
                    h2osoi_ice[c, k_perch] / (col_dz[c, k_perch] * DENICE)) / watsat[c, k_perch]
                s2 = (h2osoi_liq[c, k_perch+1] / (col_dz[c, k_perch+1] * DENH2O) +
                    h2osoi_ice[c, k_perch+1] / (col_dz[c, k_perch+1] * DENICE)) / watsat[c, k_perch+1]

                m_val = (col_z[c, k_perch+1] - col_z[c, k_perch]) / (s2 - s1)
                b_val = col_z[c, k_perch+1] - m_val * s2
                zwt_perched[c] = max(0.0, m_val * sat_lev + b_val)
            end
        end
    end

    return nothing
end

# =========================================================================
# Drainage
# =========================================================================

"""
    drainage!(temperature_t_soisno, soilhydrology, soilstate,
              waterstatebulk_ws, waterfluxbulk, col_dz, col_z, col_zi,
              col_snl, col_itype, col_landunit, col_topo_slope,
              lun_urbpoi,
              mask_hydrology, mask_urban, bounds, nlevsoi, dtime;
              use_vichydro=false)

Calculate subsurface drainage.

Ported from `Drainage` in `SoilHydrologyMod.F90`.
"""
function drainage!(
    temperature_t_soisno::Matrix{Float64},
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    waterstatebulk_ws::WaterStateData,
    waterfluxbulk::WaterFluxBulkData,
    col_dz::Matrix{Float64},
    col_z::Matrix{Float64},
    col_zi::Matrix{Float64},
    col_snl::Vector{Int},
    col_itype::Vector{Int},
    col_landunit::Vector{Int},
    col_topo_slope::Vector{Float64},
    lun_urbpoi::Vector{Bool},
    mask_hydrology::BitVector,
    mask_urban::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Float64;
    use_vichydro::Bool = false
)
    wf = waterfluxbulk.wf
    params = soilhydrology_params

    t_soisno   = temperature_t_soisno
    h2osfc     = waterstatebulk_ws.h2osfc_col
    h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col
    wa         = waterstatebulk_ws.wa_col
    aquifer_water_baseline = waterstatebulk_ws.aquifer_water_baseline

    bsw          = soilstate.bsw_col
    hksat        = soilstate.hksat_col
    sucsat       = soilstate.sucsat_col
    watsat       = soilstate.watsat_col
    eff_porosity = soilstate.eff_porosity_col
    hk_l         = soilstate.hk_l_col

    icefrac      = soilhydrology.icefrac_col
    hkdepth      = soilhydrology.hkdepth_col
    frost_table  = soilhydrology.frost_table_col
    zwt          = soilhydrology.zwt_col
    zwt_perched  = soilhydrology.zwt_perched_col
    h2osfcflag   = soilhydrology.h2osfcflag

    qflx_snwcp_liq     = wf.qflx_snwcp_liq_col
    qflx_ice_runoff_xs = wf.qflx_ice_runoff_xs_col
    qflx_drain         = wf.qflx_drain_col
    qflx_qrgwl         = wf.qflx_qrgwl_col
    qflx_rsub_sat      = wf.qflx_rsub_sat_col
    qflx_drain_perched = wf.qflx_drain_perched_col

    # Local arrays
    nb = last(bounds)
    jwt      = Vector{Int}(undef, nb)
    rsub_top = Vector{Float64}(undef, nb)
    fff      = Vector{Float64}(undef, nb)
    xsi_arr  = Vector{Float64}(undef, nb)
    xs1_arr  = Vector{Float64}(undef, nb)
    xs_arr   = Vector{Float64}(undef, nb)

    dzmm = Matrix{Float64}(undef, nb, nlevsoi)

    # Convert layer thicknesses and compute icefrac
    for j in 1:nlevsoi
        for c in bounds
            mask_hydrology[c] || continue
            dzmm[c, j] = col_dz[c, j] * 1.0e3

            vol_ice_val = min(watsat[c, j], h2osoi_ice[c, j] / (col_dz[c, j] * DENICE))
            icefrac[c, j] = min(1.0, vol_ice_val / watsat[c, j])
        end
    end

    # Initial set
    for c in bounds
        mask_hydrology[c] || continue
        qflx_drain[c]    = 0.0
        qflx_rsub_sat[c] = 0.0
        rsub_top[c]       = 0.0
    end

    # jwt: layer index right above water table
    for c in bounds
        mask_hydrology[c] || continue
        jwt[c] = nlevsoi
        for j in 1:nlevsoi
            if zwt[c] <= col_zi[c, j]
                jwt[c] = j - 1
                break
            end
        end
    end

    rous = 0.2

    # ========================== BASEFLOW ====================================
    for c in bounds
        mask_hydrology[c] || continue

        q_perch_max = params.perched_baseflow_scalar * sin(col_topo_slope[c] * (RPI / 180.0))

        # define frost table
        if t_soisno[c, 1] > TFRZ
            k_frz = nlevsoi
        else
            k_frz = 1
        end

        for k in 2:nlevsoi
            if t_soisno[c, k-1] > TFRZ && t_soisno[c, k] <= TFRZ
                k_frz = k
                break
            end
        end

        frost_table[c] = col_z[c, k_frz]
        zwt_perched[c] = frost_table[c]
        qflx_drain_perched[c] = 0.0

        # ======= water table above frost table ========
        if zwt[c] < frost_table[c] && t_soisno[c, k_frz] <= TFRZ
            # compute drainage from perched saturated region
            wtsub = 0.0
            q_perch = 0.0
            for k in (jwt[c]+1):k_frz
                imped = 10.0^(-params.e_ice * (0.5 * (icefrac[c, k] + icefrac[c, min(nlevsoi, k+1)])))
                q_perch = q_perch + imped * hksat[c, k] * dzmm[c, k]
                wtsub = wtsub + dzmm[c, k]
            end
            if wtsub > 0.0
                q_perch = q_perch / wtsub
            end

            qflx_drain_perched[c] = q_perch_max * q_perch * (frost_table[c] - zwt[c])

            # remove drainage from perched saturated layers
            rsub_top_tot = -qflx_drain_perched[c] * dtime
            for k in (jwt[c]+1):k_frz
                rsub_top_layer = max(rsub_top_tot, -(h2osoi_liq[c, k] - WATMIN))
                rsub_top_layer = min(rsub_top_layer, 0.0)
                rsub_top_tot = rsub_top_tot - rsub_top_layer

                h2osoi_liq[c, k] = h2osoi_liq[c, k] + rsub_top_layer

                if rsub_top_tot >= 0.0
                    zwt[c] = zwt[c] - rsub_top_layer / eff_porosity[c, k] / 1000.0
                    break
                else
                    zwt[c] = col_zi[c, k]
                end
            end

            qflx_drain_perched[c] = qflx_drain_perched[c] + rsub_top_tot / dtime

            # recompute jwt
            jwt[c] = nlevsoi
            for j in 1:nlevsoi
                if zwt[c] <= col_zi[c, j]
                    jwt[c] = j - 1
                    break
                end
            end
        else
            # ======= water table below frost table ========
            sat_lev = 0.9

            k_perch = 1
            for k in k_frz:-1:1
                h2osoi_vol_val = h2osoi_liq[c, k] / (col_dz[c, k] * DENH2O) +
                    h2osoi_ice[c, k] / (col_dz[c, k] * DENICE)

                if h2osoi_vol_val / watsat[c, k] <= sat_lev
                    k_perch = k
                    break
                end
            end

            if t_soisno[c, k_frz] > TFRZ
                k_perch = k_frz
            end

            # if perched water table exists
            if k_frz > k_perch
                s1 = (h2osoi_liq[c, k_perch] / (col_dz[c, k_perch] * DENH2O) +
                    h2osoi_ice[c, k_perch] / (col_dz[c, k_perch] * DENICE)) / watsat[c, k_perch]
                s2 = (h2osoi_liq[c, k_perch+1] / (col_dz[c, k_perch+1] * DENH2O) +
                    h2osoi_ice[c, k_perch+1] / (col_dz[c, k_perch+1] * DENICE)) / watsat[c, k_perch+1]

                m_val = (col_z[c, k_perch+1] - col_z[c, k_perch]) / (s2 - s1)
                b_val = col_z[c, k_perch+1] - m_val * s2
                zwt_perched[c] = max(0.0, m_val * sat_lev + b_val)

                # compute drainage from perched saturated region
                wtsub = 0.0
                q_perch = 0.0
                for k in k_perch:k_frz
                    imped = 10.0^(-params.e_ice * (0.5 * (icefrac[c, k] + icefrac[c, min(nlevsoi, k+1)])))
                    q_perch = q_perch + imped * hksat[c, k] * dzmm[c, k]
                    wtsub = wtsub + dzmm[c, k]
                end
                if wtsub > 0.0
                    q_perch = q_perch / wtsub
                end

                qflx_drain_perched[c] = q_perch_max * q_perch * (frost_table[c] - zwt_perched[c])

                # remove drainage from perched saturated layers
                rsub_top_tot = -qflx_drain_perched[c] * dtime
                for k in (k_perch+1):k_frz
                    rsub_top_layer = max(rsub_top_tot, -(h2osoi_liq[c, k] - WATMIN))
                    rsub_top_layer = min(rsub_top_layer, 0.0)
                    rsub_top_tot = rsub_top_tot - rsub_top_layer

                    h2osoi_liq[c, k] = h2osoi_liq[c, k] + rsub_top_layer

                    if rsub_top_tot >= 0.0
                        zwt_perched[c] = zwt_perched[c] - rsub_top_layer / eff_porosity[c, k] / 1000.0
                        break
                    else
                        zwt_perched[c] = col_zi[c, k]
                    end
                end

                qflx_drain_perched[c] = qflx_drain_perched[c] + rsub_top_tot / dtime
            else
                qflx_drain_perched[c] = 0.0
            end

            # -- Topographic runoff --
            fff[c] = 1.0 / hkdepth[c]
            dzsum = 0.0
            icefracsum = 0.0
            for j in max(jwt[c], 1):nlevsoi
                dzsum = dzsum + dzmm[c, j]
                icefracsum = icefracsum + icefrac[c, j] * dzmm[c, j]
            end

            if !use_vichydro
                imped = 10.0^(-params.e_ice * (icefracsum / dzsum))
                rsub_top_max = 10.0 * sin((RPI / 180.0) * col_topo_slope[c])
            else
                # VIC hydrology path — not commonly used, stubbed
                imped = 10.0^(-params.e_ice * (icefracsum / dzsum))
                rsub_top_max = 10.0 * sin((RPI / 180.0) * col_topo_slope[c])
            end

            rsub_top[c] = imped * rsub_top_max * exp(-fff[c] * zwt[c])

            # analytical expression for aquifer specific yield
            rous_local = watsat[c, nlevsoi] *
                (1.0 - (1.0 + 1.0e3 * zwt[c] / sucsat[c, nlevsoi])^(-1.0 / bsw[c, nlevsoi]))
            rous_local = max(rous_local, params.aq_sp_yield_min)

            if jwt[c] == nlevsoi
                # water table below soil column
                wa[c] = wa[c] - rsub_top[c] * dtime
                zwt[c] = zwt[c] + (rsub_top[c] * dtime) / 1000.0 / rous_local
                h2osoi_liq[c, nlevsoi] = h2osoi_liq[c, nlevsoi] + max(0.0, wa[c] - aquifer_water_baseline)
                wa[c] = min(wa[c], aquifer_water_baseline)
            else
                # water table within soil layers
                rsub_top_tot = -rsub_top[c] * dtime

                if rsub_top_tot > 0.0
                    error("RSUB_TOP IS POSITIVE in Drainage!")
                else
                    for j in (jwt[c]+1):nlevsoi
                        s_y = watsat[c, j] *
                            (1.0 - (1.0 + 1.0e3 * zwt[c] / sucsat[c, j])^(-1.0 / bsw[c, j]))
                        s_y = max(s_y, params.aq_sp_yield_min)

                        rsub_top_layer = max(rsub_top_tot, -(s_y * (col_zi[c, j] - zwt[c]) * 1.0e3))
                        rsub_top_layer = min(rsub_top_layer, 0.0)
                        h2osoi_liq[c, j] = h2osoi_liq[c, j] + rsub_top_layer

                        rsub_top_tot = rsub_top_tot - rsub_top_layer

                        if rsub_top_tot >= 0.0
                            zwt[c] = zwt[c] - rsub_top_layer / s_y / 1000.0
                            break
                        else
                            zwt[c] = col_zi[c, j]
                        end
                    end

                    # remove residual
                    zwt[c] = zwt[c] - rsub_top_tot / 1000.0 / rous_local
                    wa[c] = wa[c] + rsub_top_tot
                end

                # recompute jwt
                jwt[c] = nlevsoi
                for j in 1:nlevsoi
                    if zwt[c] <= col_zi[c, j]
                        jwt[c] = j - 1
                        break
                    end
                end
            end

            zwt[c] = max(0.0, zwt[c])
            zwt[c] = min(80.0, zwt[c])
        end
    end

    # Excessive water above saturation: redistribute upward
    for j in nlevsoi:-1:2
        for c in bounds
            mask_hydrology[c] || continue
            xsi_arr[c]       = max(h2osoi_liq[c, j] - eff_porosity[c, j] * dzmm[c, j], 0.0)
            h2osoi_liq[c, j]   = min(eff_porosity[c, j] * dzmm[c, j], h2osoi_liq[c, j])
            h2osoi_liq[c, j-1] = h2osoi_liq[c, j-1] + xsi_arr[c]
        end
    end

    for c in bounds
        mask_hydrology[c] || continue

        # watmin addition to fix water balance errors
        xs1_arr[c] = max(max(h2osoi_liq[c, 1] - WATMIN, 0.0) -
            max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_ice[c, 1] - WATMIN), 0.0)
        h2osoi_liq[c, 1] = h2osoi_liq[c, 1] - xs1_arr[c]

        l = col_landunit[c]
        if lun_urbpoi[l]
            qflx_rsub_sat[c] = xs1_arr[c] / dtime
        else
            if h2osfcflag == 1
                h2osfc[c] = h2osfc[c] + xs1_arr[c]
                qflx_rsub_sat[c] = 0.0
            else
                qflx_rsub_sat[c] = xs1_arr[c] / dtime
            end
        end

        # ice check
        xs1_arr[c] = max(max(h2osoi_ice[c, 1], 0.0) -
            max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_liq[c, 1]), 0.0)
        h2osoi_ice[c, 1] = min(max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_liq[c, 1]),
            h2osoi_ice[c, 1])
        qflx_ice_runoff_xs[c] = xs1_arr[c] / dtime
    end

    # Limit h2osoi_liq >= watmin: pull from below
    for j in 1:(nlevsoi-1)
        for c in bounds
            mask_hydrology[c] || continue
            if h2osoi_liq[c, j] < WATMIN
                xs_arr[c] = WATMIN - h2osoi_liq[c, j]
                if j == jwt[c]
                    zwt[c] = zwt[c] + xs_arr[c] / eff_porosity[c, j] / 1000.0
                end
            else
                xs_arr[c] = 0.0
            end
            h2osoi_liq[c, j]   = h2osoi_liq[c, j]   + xs_arr[c]
            h2osoi_liq[c, j+1] = h2osoi_liq[c, j+1] - xs_arr[c]
        end
    end

    # Get water for bottom layer from layers above if possible
    j = nlevsoi
    for c in bounds
        mask_hydrology[c] || continue
        if h2osoi_liq[c, j] < WATMIN
            xs_arr[c] = WATMIN - h2osoi_liq[c, j]
            for i in (nlevsoi-1):-1:1
                available_h2osoi_liq = max(h2osoi_liq[c, i] - WATMIN - xs_arr[c], 0.0)
                if available_h2osoi_liq >= xs_arr[c]
                    h2osoi_liq[c, j] = h2osoi_liq[c, j] + xs_arr[c]
                    h2osoi_liq[c, i] = h2osoi_liq[c, i] - xs_arr[c]
                    xs_arr[c] = 0.0
                    break
                else
                    h2osoi_liq[c, j] = h2osoi_liq[c, j] + available_h2osoi_liq
                    h2osoi_liq[c, i] = h2osoi_liq[c, i] - available_h2osoi_liq
                    xs_arr[c] = xs_arr[c] - available_h2osoi_liq
                end
            end
        else
            xs_arr[c] = 0.0
        end
        h2osoi_liq[c, j] = h2osoi_liq[c, j] + xs_arr[c]
        rsub_top[c] = rsub_top[c] - xs_arr[c] / dtime
    end

    for c in bounds
        mask_hydrology[c] || continue
        qflx_drain[c] = qflx_rsub_sat[c] + rsub_top[c]
        qflx_qrgwl[c] = qflx_snwcp_liq[c]
    end

    # No drainage for urban columns (except pervious road)
    for c in bounds
        mask_urban[c] || continue
        if col_itype[c] != ICOL_ROAD_PERV
            qflx_drain[c] = 0.0
            qflx_qrgwl[c] = qflx_snwcp_liq[c]
        end
    end

    return nothing
end

# =========================================================================
# CLMVICMap
# =========================================================================

"""
    clm_vic_map!(soilhydrology, waterstatebulk_ws, col_dz, col_zi, col_z,
                  mask, bounds, nlevsoi, nlevgrnd, nlayer, nlayert)

Map CLM layers to VIC layers.

Ported from `CLMVICMap` in `SoilHydrologyMod.F90`.
"""
function clm_vic_map!(
    soilhydrology::SoilHydrologyData,
    waterstatebulk_ws::WaterStateData,
    col_dz::Matrix{Float64},
    col_zi::Matrix{Float64},
    col_z::Matrix{Float64},
    mask::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int,
    nlevgrnd::Int,
    nlayer::Int,
    nlayert::Int
)
    h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col
    h2osoi_vol = waterstatebulk_ws.h2osoi_vol_col

    depth         = soilhydrology.depth_col
    porosity      = soilhydrology.porosity_col
    max_moist     = soilhydrology.max_moist_col
    vic_clm_fract = soilhydrology.vic_clm_fract_col
    moist         = soilhydrology.moist_col
    ice           = soilhydrology.ice_col
    moist_vol     = soilhydrology.moist_vol_col
    top_moist     = soilhydrology.top_moist_col
    top_max_moist = soilhydrology.top_max_moist_col
    top_ice       = soilhydrology.top_ice_col
    top_moist_limited = soilhydrology.top_moist_limited_col

    for c in bounds
        mask[c] || continue

        for i in 1:nlayer
            ice0 = ice[c, i]
            moist0 = moist[c, i]
            ice[c, i] = 0.0
            moist[c, i] = 0.0
            for j in 1:nlevsoi
                ice[c, i] = ice[c, i] + h2osoi_ice[c, j] * vic_clm_fract[c, i, j]
                moist[c, i] = moist[c, i] + h2osoi_liq[c, j] * vic_clm_fract[c, i, j]
            end
            ice[c, i]       = min(moist0 + ice0, ice[c, i])
            ice[c, i]       = max(0.0, ice[c, i])
            moist[c, i]     = max(WATMIN, moist[c, i])
            moist[c, i]     = min(max_moist[c, i] - ice[c, i], moist[c, i])
            moist_vol[c, i] = moist[c, i] / (depth[c, i] * DENICE) + ice[c, i] / (depth[c, i] * DENH2O)
            moist_vol[c, i] = min(porosity[c, i], moist_vol[c, i])
            moist_vol[c, i] = max(0.01, moist_vol[c, i])
        end

        # hydrologic inactive layers
        for k in (nlayer+1):nlayert
            j_clm = nlevsoi + (k - nlayer)
            if j_clm <= nlevgrnd
                ice[c, k]       = h2osoi_ice[c, j_clm]
                moist[c, k]     = h2osoi_liq[c, j_clm]
                moist_vol[c, k] = h2osoi_vol[c, j_clm]
            end
        end
    end

    # Set values related to top VIC layers
    for c in bounds
        mask[c] || continue
        top_moist[c]     = 0.0
        top_ice[c]       = 0.0
        top_max_moist[c] = 0.0
    end

    for j in 1:(nlayer - 1)
        for c in bounds
            mask[c] || continue
            top_ice[c]       = top_ice[c] + ice[c, j]
            top_moist[c]     = top_moist[c] + moist[c, j] + ice[c, j]
            top_max_moist[c] = top_max_moist[c] + max_moist[c, j]
        end
    end

    for c in bounds
        mask[c] || continue
        top_moist_limited[c] = min(top_moist[c], top_max_moist[c])
    end

    return nothing
end

# =========================================================================
# PerchedWaterTable
# =========================================================================

"""
    perched_water_table!(soilhydrology, soilstate, temperature_t_soisno,
                          waterstatebulk_ws, col_dz, col_z, col_zi,
                          mask_hydrology, bounds, nlevsoi)

Calculate perched water table location.

Ported from `PerchedWaterTable` in `SoilHydrologyMod.F90`.
"""
function perched_water_table!(
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    temperature_t_soisno::Matrix{Float64},
    waterstatebulk_ws::WaterStateData,
    col_dz::Matrix{Float64},
    col_z::Matrix{Float64},
    col_zi::Matrix{Float64},
    mask_hydrology::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int
)
    t_soisno   = temperature_t_soisno
    h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col
    h2osoi_vol = waterstatebulk_ws.h2osoi_vol_col
    watsat     = soilstate.watsat_col
    zwt        = soilhydrology.zwt_col
    zwt_perched = soilhydrology.zwt_perched_col
    frost_table = soilhydrology.frost_table_col

    sat_lev = 0.9

    for c in bounds
        mask_hydrology[c] || continue

        # define frost table
        if t_soisno[c, 1] > TFRZ
            k_frz = nlevsoi
        else
            k_frz = 1
        end

        for k in 2:nlevsoi
            if t_soisno[c, k-1] > TFRZ && t_soisno[c, k] <= TFRZ
                k_frz = k
                break
            end
        end

        # frost table is top of frozen layer
        frost_table[c] = col_zi[c, k_frz - 1]
        zwt_perched[c] = frost_table[c]

        # water table above frost table: do nothing
        if zwt[c] < frost_table[c] && t_soisno[c, k_frz] <= TFRZ
            # do nothing
        elseif k_frz > 1
            # water table below frost table
            k_perch = 1
            for k in k_frz:-1:1
                h2osoi_vol[c, k] = h2osoi_liq[c, k] / (col_dz[c, k] * DENH2O) +
                    h2osoi_ice[c, k] / (col_dz[c, k] * DENICE)

                if h2osoi_vol[c, k] / watsat[c, k] <= sat_lev
                    k_perch = k
                    break
                end
            end

            if t_soisno[c, k_frz] > TFRZ
                k_perch = k_frz
            end

            if k_frz > k_perch
                s1 = (h2osoi_liq[c, k_perch] / (col_dz[c, k_perch] * DENH2O) +
                    h2osoi_ice[c, k_perch] / (col_dz[c, k_perch] * DENICE)) / watsat[c, k_perch]
                s2 = (h2osoi_liq[c, k_perch+1] / (col_dz[c, k_perch+1] * DENH2O) +
                    h2osoi_ice[c, k_perch+1] / (col_dz[c, k_perch+1] * DENICE)) / watsat[c, k_perch+1]

                if s1 > s2
                    zwt_perched[c] = col_zi[c, k_perch - 1]
                else
                    m_val = (col_z[c, k_perch+1] - col_z[c, k_perch]) / (s2 - s1)
                    b_val = col_z[c, k_perch+1] - m_val * s2
                    zwt_perched[c] = max(0.0, m_val * sat_lev + b_val)
                end
            end
        end
    end

    return nothing
end

# =========================================================================
# ThetaBasedWaterTable
# =========================================================================

"""
    theta_based_water_table!(soilhydrology, soilstate, waterstatebulk_ws,
                              col_dz, col_z, col_zi, col_nbedrock,
                              mask_hydrology, bounds, nlevsoi)

Calculate water table from soil moisture state.

Ported from `ThetaBasedWaterTable` in `SoilHydrologyMod.F90`.
"""
function theta_based_water_table!(
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    waterstatebulk_ws::WaterStateData,
    col_dz::Matrix{Float64},
    col_z::Matrix{Float64},
    col_zi::Matrix{Float64},
    col_nbedrock::Vector{Int},
    mask_hydrology::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int
)
    h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col
    h2osoi_vol = waterstatebulk_ws.h2osoi_vol_col
    watsat     = soilstate.watsat_col
    zwt        = soilhydrology.zwt_col

    sat_lev = 0.9

    for c in bounds
        mask_hydrology[c] || continue

        zwt[c] = col_zi[c, nlevsoi]

        k_zwt = col_nbedrock[c]
        sat_flag = 1
        for k in col_nbedrock[c]:-1:1
            h2osoi_vol[c, k] = h2osoi_liq[c, k] / (col_dz[c, k] * DENH2O) +
                h2osoi_ice[c, k] / (col_dz[c, k] * DENICE)

            if h2osoi_vol[c, k] / watsat[c, k] <= sat_lev
                k_zwt = k
                sat_flag = 0
                break
            end
        end
        if sat_flag == 1
            k_zwt = 1
        end

        if k_zwt == 1
            zwt[c] = col_zi[c, 1]
        elseif k_zwt < col_nbedrock[c]
            s1 = (h2osoi_liq[c, k_zwt] / (col_dz[c, k_zwt] * DENH2O) +
                h2osoi_ice[c, k_zwt] / (col_dz[c, k_zwt] * DENICE)) / watsat[c, k_zwt]
            s2 = (h2osoi_liq[c, k_zwt+1] / (col_dz[c, k_zwt+1] * DENH2O) +
                h2osoi_ice[c, k_zwt+1] / (col_dz[c, k_zwt+1] * DENICE)) / watsat[c, k_zwt+1]

            m_val = (col_z[c, k_zwt+1] - col_z[c, k_zwt]) / (s2 - s1)
            b_val = col_z[c, k_zwt+1] - m_val * s2
            zwt[c] = max(0.0, m_val * sat_lev + b_val)
        else
            zwt[c] = col_zi[c, col_nbedrock[c]]
        end
    end

    return nothing
end

# =========================================================================
# RenewCondensation
# =========================================================================

"""
    renew_condensation!(waterstatebulk_ws, waterdiagbulk, waterfluxbulk,
                         col_snl, col_itype, mask_hydrology, mask_urban,
                         bounds, dtime)

Renew ice and liquid mass due to condensation.

Ported from `RenewCondensation` in `SoilHydrologyMod.F90`.
"""
function renew_condensation!(
    waterstatebulk_ws::WaterStateData,
    waterdiagbulk::WaterDiagnosticBulkData,
    waterfluxbulk::WaterFluxBulkData,
    col_snl::Vector{Int},
    col_itype::Vector{Int},
    mask_hydrology::BitVector,
    mask_urban::BitVector,
    bounds::UnitRange{Int},
    dtime::Float64
)
    wf = waterfluxbulk.wf

    h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col
    frac_h2osfc = waterdiagbulk.frac_h2osfc_col
    qflx_liqdew_to_top_layer      = wf.qflx_liqdew_to_top_layer_col
    qflx_soliddew_to_top_layer    = wf.qflx_soliddew_to_top_layer_col
    qflx_solidevap_from_top_layer = wf.qflx_solidevap_from_top_layer_col

    for c in bounds
        mask_hydrology[c] || continue

        if col_snl[c] + 1 >= 1
            h2osoi_liq[c, 1] = h2osoi_liq[c, 1] +
                (1.0 - frac_h2osfc[c]) * qflx_liqdew_to_top_layer[c] * dtime
            h2osoi_ice[c, 1] = h2osoi_ice[c, 1] +
                (1.0 - frac_h2osfc[c]) * qflx_soliddew_to_top_layer[c] * dtime
            h2osoi_ice_before = h2osoi_ice[c, 1]
            h2osoi_ice[c, 1] = h2osoi_ice[c, 1] -
                (1.0 - frac_h2osfc[c]) * qflx_solidevap_from_top_layer[c] * dtime

            # truncate small values (simplified from Fortran call)
            if h2osoi_ice[c, 1] < 0.0 && abs(h2osoi_ice[c, 1]) < TOLERANCE_SOILHYDRO * abs(h2osoi_ice_before)
                h2osoi_ice[c, 1] = 0.0
            end

            if h2osoi_ice[c, 1] < 0.0
                error("In RenewCondensation, h2osoi_ice has gone significantly negative at c=$c")
            end
        end
    end

    # Urban columns
    for c in bounds
        mask_urban[c] || continue
        if col_itype[c] == ICOL_ROOF || col_itype[c] == ICOL_ROAD_IMPERV
            if col_snl[c] + 1 >= 1
                h2osoi_liq[c, 1] = h2osoi_liq[c, 1] +
                    qflx_liqdew_to_top_layer[c] * dtime
                h2osoi_ice[c, 1] = h2osoi_ice[c, 1] +
                    qflx_soliddew_to_top_layer[c] * dtime
                h2osoi_ice_before = h2osoi_ice[c, 1]
                h2osoi_ice[c, 1] = h2osoi_ice[c, 1] -
                    qflx_solidevap_from_top_layer[c] * dtime

                if h2osoi_ice[c, 1] < 0.0 && abs(h2osoi_ice[c, 1]) < TOLERANCE_SOILHYDRO * abs(h2osoi_ice_before)
                    h2osoi_ice[c, 1] = 0.0
                end

                if h2osoi_ice[c, 1] < 0.0
                    error("In RenewCondensation (urban), h2osoi_ice has gone significantly negative at c=$c")
                end
            end
        end
    end

    return nothing
end

# =========================================================================
# CalcIrrigWithdrawals
# =========================================================================

"""
    calc_irrig_withdrawals!(soilhydrology, soilstate,
                             qflx_gw_demand, qflx_gw_uncon_irrig_lyr, qflx_gw_con_irrig,
                             col_nbedrock, col_zi,
                             mask_soil, bounds, nlevsoi, dtime)

Calculate irrigation withdrawals from groundwater by layer.

Ported from `CalcIrrigWithdrawals` in `SoilHydrologyMod.F90`.
"""
function calc_irrig_withdrawals!(
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    qflx_gw_demand::Vector{Float64},
    qflx_gw_uncon_irrig_lyr::Matrix{Float64},
    qflx_gw_con_irrig::Vector{Float64},
    col_nbedrock::Vector{Int},
    col_zi::Matrix{Float64},
    mask_soil::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Float64
)
    params = soilhydrology_params
    bsw    = soilstate.bsw_col
    sucsat = soilstate.sucsat_col
    watsat = soilstate.watsat_col
    zwt    = soilhydrology.zwt_col

    # Initialize
    for j in 1:nlevsoi
        for c in bounds
            mask_soil[c] || continue
            qflx_gw_uncon_irrig_lyr[c, j] = 0.0
        end
    end

    for c in bounds
        mask_soil[c] || continue

        irrig_demand_remaining = qflx_gw_demand[c] * dtime

        if irrig_demand_remaining < 0.0
            error("negative groundwater irrigation demand at c=$c!")
        end

        jwt_c = nlevsoi
        for j in 1:nlevsoi
            if zwt[c] <= col_zi[c, j]
                jwt_c = j - 1
                break
            end
        end

        for j in (jwt_c+1):col_nbedrock[c]
            s_y = watsat[c, j] *
                (1.0 - (1.0 + 1.0e3 * zwt[c] / sucsat[c, j])^(-1.0 / bsw[c, j]))
            s_y = max(s_y, params.aq_sp_yield_min)

            if j == jwt_c + 1
                available_water_layer = max(0.0, s_y * (col_zi[c, j] - zwt[c]) * 1.0e3)
            else
                available_water_layer = max(0.0, s_y * (col_zi[c, j] - col_zi[c, j-1]) * 1.0e3)
            end

            irrig_layer = min(irrig_demand_remaining, available_water_layer)
            qflx_gw_uncon_irrig_lyr[c, j] = irrig_layer / dtime

            irrig_demand_remaining = irrig_demand_remaining - irrig_layer

            if irrig_demand_remaining <= 0.0
                break
            end
        end

        if irrig_demand_remaining > 0.0
            qflx_gw_con_irrig[c] = irrig_demand_remaining / dtime
        else
            qflx_gw_con_irrig[c] = 0.0
        end
    end

    return nothing
end

# =========================================================================
# WithdrawGroundwaterIrrigation
# =========================================================================

"""
    withdraw_groundwater_irrigation!(waterflux_wf, waterstate_ws,
                                      mask_soil, bounds, nlevsoi, dtime)

Remove groundwater irrigation from unconfined and confined aquifers.

Ported from `WithdrawGroundwaterIrrigation` in `SoilHydrologyMod.F90`.
"""
function withdraw_groundwater_irrigation!(
    waterflux_wf::WaterFluxData,
    waterstate_ws::WaterStateData,
    mask_soil::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Float64
)
    qflx_gw_uncon_irrig_lyr = waterflux_wf.qflx_gw_uncon_irrig_lyr_col
    qflx_gw_con_irrig       = waterflux_wf.qflx_gw_con_irrig_col
    wa                       = waterstate_ws.wa_col
    h2osoi_liq               = waterstate_ws.h2osoi_liq_col

    for j in 1:nlevsoi
        for c in bounds
            mask_soil[c] || continue
            h2osoi_liq[c, j] = h2osoi_liq[c, j] - qflx_gw_uncon_irrig_lyr[c, j] * dtime
        end
    end

    for c in bounds
        mask_soil[c] || continue
        wa[c] = wa[c] - qflx_gw_con_irrig[c] * dtime
    end

    return nothing
end
