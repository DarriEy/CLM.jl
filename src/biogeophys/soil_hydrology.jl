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
    baseflow_scalar::Real = 1.0e-2
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
# KernelAbstractions kernels for independent per-element loops
# (each (column) or (column,layer) iteration is fully independent — no
#  loop-carried dependencies, accumulation, or inter-column coupling).
# =========================================================================

# ---- set_soil_water_fractions! : per-(column,layer) diagnostic fractions ----
@kernel function _soilhyd_water_fractions_kernel!(eff_porosity, icefrac, @Const(mask),
                                                  @Const(watsat), @Const(col_dz),
                                                  @Const(h2osoi_ice), @Const(excess_ice),
                                                  joff::Int, denice)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        T = eltype(eff_porosity)
        dn = T(denice)
        dz_ext = col_dz[c, j + joff] + excess_ice[c, j] / dn
        vol_ice_val = smooth_min(watsat[c, j], (h2osoi_ice[c, j + joff] + excess_ice[c, j]) / (dz_ext * dn))
        eff_porosity[c, j] = smooth_max(T(0.01), watsat[c, j] - vol_ice_val)
        icefrac[c, j] = smooth_min(one(T), vol_ice_val / watsat[c, j])
    end
end

soilhyd_water_fractions!(eff_porosity, icefrac, mask, watsat, col_dz,
                         h2osoi_ice, excess_ice, joff::Int, nlevsoi::Int, denice) =
    _launch!(_soilhyd_water_fractions_kernel!, eff_porosity, icefrac, mask, watsat,
             col_dz, h2osoi_ice, excess_ice, joff, convert(eltype(eff_porosity), denice);
             ndrange = (length(mask), nlevsoi))

# ---- set_floodc! : apply gridcell flood flux to non-lake columns ----
@kernel function _soilhyd_floodc_kernel!(qflx_floodc, @Const(qflx_floodg),
                                         @Const(col_gridcell), @Const(col_itype),
                                         @Const(mask_nolake), icol_sunwall::Int,
                                         icol_shadewall::Int)
    c = @index(Global)
    @inbounds if mask_nolake[c]
        g = col_gridcell[c]
        if col_itype[c] != icol_sunwall && col_itype[c] != icol_shadewall
            qflx_floodc[c] = qflx_floodg[g]
        else
            qflx_floodc[c] = zero(eltype(qflx_floodc))
        end
    end
end

soilhyd_floodc!(qflx_floodc, qflx_floodg, col_gridcell, col_itype, mask_nolake,
                icol_sunwall::Int, icol_shadewall::Int) =
    _launch!(_soilhyd_floodc_kernel!, qflx_floodc, qflx_floodg, col_gridcell,
             col_itype, mask_nolake, icol_sunwall, icol_shadewall)

# ---- set_qflx_inputs! : partition surface inputs between soil and h2osfc ----
@kernel function _soilhyd_qflx_inputs_kernel!(qflx_top_soil, qflx_in_soil,
                                              qflx_top_soil_to_h2osfc, @Const(mask),
                                              @Const(col_snl), @Const(qflx_rain_plus_snomelt),
                                              @Const(qflx_snow_h2osfc), @Const(qflx_floodc),
                                              @Const(qflx_liqevap_from_top_layer),
                                              @Const(qflx_ev_soil), @Const(qflx_ev_h2osfc),
                                              @Const(qflx_sat_excess_surf), @Const(frac_sno),
                                              @Const(frac_h2osfc))
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(qflx_top_soil)
        qflx_top_soil[c] = qflx_rain_plus_snomelt[c] + qflx_snow_h2osfc[c] + qflx_floodc[c]

        if col_snl[c] >= 0
            fsno = zero(T)
            qflx_evap = qflx_liqevap_from_top_layer[c]
        else
            fsno = frac_sno[c]
            qflx_evap = qflx_ev_soil[c]
        end

        qflx_in_soil[c] = (one(T) - frac_h2osfc[c]) * (qflx_top_soil[c] - qflx_sat_excess_surf[c])
        qflx_top_soil_to_h2osfc[c] = frac_h2osfc[c] * (qflx_top_soil[c] - qflx_sat_excess_surf[c])

        qflx_in_soil[c] = qflx_in_soil[c] - (one(T) - fsno - frac_h2osfc[c]) * qflx_evap
        qflx_top_soil_to_h2osfc[c] = qflx_top_soil_to_h2osfc[c] - frac_h2osfc[c] * qflx_ev_h2osfc[c]
    end
end

soilhyd_qflx_inputs!(qflx_top_soil, qflx_in_soil, qflx_top_soil_to_h2osfc, mask,
                     col_snl, qflx_rain_plus_snomelt, qflx_snow_h2osfc, qflx_floodc,
                     qflx_liqevap_from_top_layer, qflx_ev_soil, qflx_ev_h2osfc,
                     qflx_sat_excess_surf, frac_sno, frac_h2osfc) =
    _launch!(_soilhyd_qflx_inputs_kernel!, qflx_top_soil, qflx_in_soil,
             qflx_top_soil_to_h2osfc, mask, col_snl, qflx_rain_plus_snomelt,
             qflx_snow_h2osfc, qflx_floodc, qflx_liqevap_from_top_layer,
             qflx_ev_soil, qflx_ev_h2osfc, qflx_sat_excess_surf, frac_sno, frac_h2osfc)

# ---- route_infiltration_excess! : route infiltration excess runoff ----
@kernel function _soilhyd_route_infl_excess_kernel!(qflx_in_soil_limited, qflx_in_h2osfc,
                                                    qflx_infl_excess_surf, @Const(mask),
                                                    @Const(col_landunit), @Const(lun_itype),
                                                    @Const(qflx_in_soil),
                                                    @Const(qflx_top_soil_to_h2osfc),
                                                    @Const(qflx_infl_excess),
                                                    h2osfcflag::Int, istsoil::Int, istcrop::Int)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(qflx_in_soil_limited)
        l = col_landunit[c]
        if lun_itype[l] == istsoil || lun_itype[l] == istcrop
            qflx_in_soil_limited[c] = qflx_in_soil[c] - qflx_infl_excess[c]
            if h2osfcflag != 0
                qflx_in_h2osfc[c] = qflx_top_soil_to_h2osfc[c] + qflx_infl_excess[c]
                qflx_infl_excess_surf[c] = zero(T)
            else
                qflx_in_h2osfc[c] = qflx_top_soil_to_h2osfc[c]
                qflx_infl_excess_surf[c] = qflx_infl_excess[c]
            end
        else
            qflx_in_soil_limited[c] = qflx_in_soil[c]
            qflx_in_h2osfc[c] = zero(T)
            qflx_infl_excess_surf[c] = zero(T)
        end
    end
end

soilhyd_route_infl_excess!(qflx_in_soil_limited, qflx_in_h2osfc, qflx_infl_excess_surf,
                           mask, col_landunit, lun_itype, qflx_in_soil,
                           qflx_top_soil_to_h2osfc, qflx_infl_excess,
                           h2osfcflag::Int, istsoil::Int, istcrop::Int) =
    _launch!(_soilhyd_route_infl_excess_kernel!, qflx_in_soil_limited, qflx_in_h2osfc,
             qflx_infl_excess_surf, mask, col_landunit, lun_itype, qflx_in_soil,
             qflx_top_soil_to_h2osfc, qflx_infl_excess, h2osfcflag, istsoil, istcrop)

# ---- infiltration! : total infiltration ----
@kernel function _soilhyd_infiltration_kernel!(qflx_infl, @Const(mask),
                                               @Const(qflx_in_soil_limited),
                                               @Const(qflx_h2osfc_drain))
    c = @index(Global)
    @inbounds if mask[c]
        qflx_infl[c] = qflx_in_soil_limited[c] + qflx_h2osfc_drain[c]
    end
end

soilhyd_infiltration!(qflx_infl, mask, qflx_in_soil_limited, qflx_h2osfc_drain) =
    _launch!(_soilhyd_infiltration_kernel!, qflx_infl, mask, qflx_in_soil_limited,
             qflx_h2osfc_drain)

# ---- total_surface_runoff! : qflx_surf for hydrologically-active columns ----
@kernel function _soilhyd_surf_runoff_kernel!(qflx_surf, @Const(mask),
                                              @Const(qflx_sat_excess_surf),
                                              @Const(qflx_infl_excess_surf),
                                              @Const(qflx_h2osfc_surf))
    c = @index(Global)
    @inbounds if mask[c]
        qflx_surf[c] = qflx_sat_excess_surf[c] + qflx_infl_excess_surf[c] + qflx_h2osfc_surf[c]
    end
end

soilhyd_surf_runoff!(qflx_surf, mask, qflx_sat_excess_surf, qflx_infl_excess_surf,
                     qflx_h2osfc_surf) =
    _launch!(_soilhyd_surf_runoff_kernel!, qflx_surf, mask, qflx_sat_excess_surf,
             qflx_infl_excess_surf, qflx_h2osfc_surf)

# ---- total_surface_runoff! : qflx_surf for non-hydrologic urban columns ----
@kernel function _soilhyd_surf_runoff_urban_kernel!(qflx_surf, xs_urban, @Const(mask_urban),
                                                    @Const(col_itype), @Const(col_snl),
                                                    @Const(qflx_rain_plus_snomelt),
                                                    @Const(h2osoi_liq),
                                                    @Const(qflx_liqevap_from_top_layer),
                                                    @Const(qflx_floodc),
                                                    icol_roof::Int, icol_road_imperv::Int,
                                                    icol_sunwall::Int, icol_shadewall::Int,
                                                    pondmx_urban, dtime)
    c = @index(Global)
    @inbounds if mask_urban[c]
        T = eltype(qflx_surf)
        dt = T(dtime)
        pmx = T(pondmx_urban)
        if col_itype[c] == icol_roof || col_itype[c] == icol_road_imperv
            if col_snl[c] < 0
                qflx_surf[c] = smooth_max(zero(T), qflx_rain_plus_snomelt[c])
            else
                xs_urban[c] = smooth_max(zero(T),
                    h2osoi_liq[c, 1] / dt + qflx_rain_plus_snomelt[c] -
                    qflx_liqevap_from_top_layer[c] - pmx / dt)
                qflx_surf[c] = xs_urban[c]
            end
            qflx_surf[c] = qflx_surf[c] + qflx_floodc[c]
        elseif col_itype[c] == icol_sunwall || col_itype[c] == icol_shadewall
            qflx_surf[c] = zero(T)
        end
    end
end

soilhyd_surf_runoff_urban!(qflx_surf, xs_urban, mask_urban, col_itype, col_snl,
                           qflx_rain_plus_snomelt, h2osoi_liq, qflx_liqevap_from_top_layer,
                           qflx_floodc, icol_roof::Int, icol_road_imperv::Int,
                           icol_sunwall::Int, icol_shadewall::Int, pondmx_urban, dtime) =
    _launch!(_soilhyd_surf_runoff_urban_kernel!, qflx_surf, xs_urban, mask_urban,
             col_itype, col_snl, qflx_rain_plus_snomelt, h2osoi_liq,
             qflx_liqevap_from_top_layer, qflx_floodc, icol_roof, icol_road_imperv,
             icol_sunwall, icol_shadewall, convert(eltype(qflx_surf), pondmx_urban),
             convert(eltype(qflx_surf), dtime))

# ---- update_urban_ponding! : ponding state on urban surfaces ----
@kernel function _soilhyd_urban_ponding_kernel!(h2osoi_liq, @Const(mask_urban),
                                                @Const(col_itype), @Const(col_snl),
                                                @Const(xs_urban),
                                                @Const(qflx_rain_plus_snomelt),
                                                @Const(qflx_liqevap_from_top_layer),
                                                icol_roof::Int, icol_road_imperv::Int,
                                                pondmx_urban, dtime)
    c = @index(Global)
    @inbounds if mask_urban[c]
        if col_itype[c] == icol_roof || col_itype[c] == icol_road_imperv
            if col_snl[c] >= 0
                if xs_urban[c] > 0.0
                    h2osoi_liq[c, 1] = pondmx_urban
                else
                    h2osoi_liq[c, 1] = smooth_max(0.0, h2osoi_liq[c, 1] +
                        (qflx_rain_plus_snomelt[c] - qflx_liqevap_from_top_layer[c]) * dtime)
                end
            end
        end
    end
end

soilhyd_urban_ponding!(h2osoi_liq, mask_urban, col_itype, col_snl, xs_urban,
                       qflx_rain_plus_snomelt, qflx_liqevap_from_top_layer,
                       icol_roof::Int, icol_road_imperv::Int, pondmx_urban, dtime) =
    _launch!(_soilhyd_urban_ponding_kernel!, h2osoi_liq, mask_urban, col_itype,
             col_snl, xs_urban, qflx_rain_plus_snomelt, qflx_liqevap_from_top_layer,
             icol_roof, icol_road_imperv, pondmx_urban, dtime)

# ---- withdraw_groundwater_irrigation! : per-(column,layer) liq withdrawal ----
@kernel function _soilhyd_withdraw_gw_lyr_kernel!(h2osoi_liq, @Const(mask),
                                                  @Const(qflx_gw_uncon_irrig_lyr), dtime)
    c, j = @index(Global, NTuple)
    @inbounds if mask[c]
        h2osoi_liq[c, j] = h2osoi_liq[c, j] - qflx_gw_uncon_irrig_lyr[c, j] * dtime
    end
end

soilhyd_withdraw_gw_lyr!(h2osoi_liq, mask, qflx_gw_uncon_irrig_lyr, nlevsoi::Int, dtime) =
    _launch!(_soilhyd_withdraw_gw_lyr_kernel!, h2osoi_liq, mask, qflx_gw_uncon_irrig_lyr,
             dtime; ndrange = (length(mask), nlevsoi))

# ---- withdraw_groundwater_irrigation! : per-column confined aquifer ----
@kernel function _soilhyd_withdraw_gw_con_kernel!(wa, @Const(mask),
                                                  @Const(qflx_gw_con_irrig), dtime)
    c = @index(Global)
    @inbounds if mask[c]
        wa[c] = wa[c] - qflx_gw_con_irrig[c] * dtime
    end
end

soilhyd_withdraw_gw_con!(wa, mask, qflx_gw_con_irrig, dtime) =
    _launch!(_soilhyd_withdraw_gw_con_kernel!, wa, mask, qflx_gw_con_irrig, dtime)

# =========================================================================
# SetSoilWaterFractions
# =========================================================================

"""
    set_soil_water_fractions!(soilhydrology, soilstate, waterstatebulk,
                               col_dz, mask_hydrology, bounds, nlevsoi, nlevsno)

Set diagnostic variables related to the fraction of water and ice in each layer.

Ported from `SetSoilWaterFractions` in `SoilHydrologyMod.F90`.
"""
function set_soil_water_fractions!(
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    waterstatebulk_ws::WaterStateData,
    col_dz::AbstractMatrix{<:Real},
    mask_hydrology::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsoi::Int,
    nlevsno::Int
)
    joff = nlevsno  # offset for snow+soil arrays (col_dz, h2osoi_ice/liq)

    watsat       = soilstate.watsat_col
    eff_porosity = soilstate.eff_porosity_col
    h2osoi_ice   = waterstatebulk_ws.h2osoi_ice_col
    excess_ice   = waterstatebulk_ws.excess_ice_col
    icefrac      = soilhydrology.icefrac_col

    soilhyd_water_fractions!(eff_porosity, icefrac, mask_hydrology, watsat, col_dz,
                             h2osoi_ice, excess_ice, joff, nlevsoi, DENICE)

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
    qflx_floodc::AbstractVector{<:Real},
    qflx_floodg::AbstractVector{<:Real},
    col_gridcell::AbstractVector{<:Integer},
    col_itype::AbstractVector{<:Integer},
    mask_nolake::AbstractVector{Bool},
    bounds::UnitRange{Int}
)
    soilhyd_floodc!(qflx_floodc, qflx_floodg, col_gridcell, col_itype, mask_nolake,
                    ICOL_SUNWALL, ICOL_SHADEWALL)

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
    col_snl::AbstractVector{<:Integer},
    mask_hydrology::AbstractVector{Bool},
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

    soilhyd_qflx_inputs!(qflx_top_soil, qflx_in_soil, qflx_top_soil_to_h2osfc,
                         mask_hydrology, col_snl, qflx_rain_plus_snomelt, qflx_snow_h2osfc,
                         qflx_floodc, qflx_liqevap_from_top_layer, qflx_ev_soil,
                         qflx_ev_h2osfc, qflx_sat_excess_surf, frac_sno, frac_h2osfc)

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
    col_landunit::AbstractVector{<:Integer},
    lun_itype::AbstractVector{<:Integer},
    mask_hydrology::AbstractVector{Bool},
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

    soilhyd_route_infl_excess!(qflx_in_soil_limited, qflx_in_h2osfc, qflx_infl_excess_surf,
                               mask_hydrology, col_landunit, lun_itype, qflx_in_soil,
                               qflx_top_soil_to_h2osfc, qflx_infl_excess,
                               h2osfcflag, ISTSOIL, ISTCROP)

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
    mask_hydrology::AbstractVector{Bool},
    bounds::UnitRange{Int}
)
    wf = waterfluxbulk.wf

    qflx_infl            = wf.qflx_infl_col
    qflx_in_soil_limited = waterfluxbulk.qflx_in_soil_limited_col
    qflx_h2osfc_drain    = waterfluxbulk.qflx_h2osfc_drain_col

    soilhyd_infiltration!(qflx_infl, mask_hydrology, qflx_in_soil_limited, qflx_h2osfc_drain)

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
    col_snl::AbstractVector{<:Integer},
    col_itype::AbstractVector{<:Integer},
    col_landunit::AbstractVector{<:Integer},
    lun_urbpoi::AbstractVector{Bool},
    mask_hydrology::AbstractVector{Bool},
    mask_urban::AbstractVector{Bool},
    bounds::UnitRange{Int},
    dtime::Real
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
    soilhyd_surf_runoff!(qflx_surf, mask_hydrology, qflx_sat_excess_surf,
                         qflx_infl_excess_surf, qflx_h2osfc_surf)

    # Set qflx_surf for non-hydrologically-active urban columns
    soilhyd_surf_runoff_urban!(qflx_surf, xs_urban, mask_urban, col_itype, col_snl,
                               qflx_rain_plus_snomelt, h2osoi_liq, qflx_liqevap_from_top_layer,
                               qflx_floodc, ICOL_ROOF, ICOL_ROAD_IMPERV, ICOL_SUNWALL,
                               ICOL_SHADEWALL, PONDMX_URBAN, dtime)

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
    dtime::Real
)
    wf = waterfluxbulk.wf

    h2osoi_liq  = waterstatebulk_ws.h2osoi_liq_col
    xs_urban    = soilhydrology.xs_urban_col
    qflx_rain_plus_snomelt = wf.qflx_rain_plus_snomelt_col
    qflx_liqevap_from_top_layer = wf.qflx_liqevap_from_top_layer_col

    soilhyd_urban_ponding!(h2osoi_liq, mask_urban, col_itype, col_snl, xs_urban,
                           qflx_rain_plus_snomelt, qflx_liqevap_from_top_layer,
                           ICOL_ROOF, ICOL_ROAD_IMPERV, PONDMX_URBAN, dtime)

    return nothing
end

# =========================================================================
# WaterTable
# =========================================================================

# --------------------------------------------------------------------------
# water_table! : ONE thread per (hydrology) column. The whole WaterTable
# subroutine is per-column with internal SEQUENTIAL nlevsoi level searches
# (jwt water-table search, qcharge rise/deepen cascade with loop-carried zwt,
# perched-table frost/saturation search) — exactly the loop-carried pattern
# from soil_temperature.jl / soil_water.jl. Each thread runs the full nested
# search in-thread and writes only its own column's outputs.
#
# The many soil-column arrays are grouped into two immutable device-view
# bundles (Adapt.@adapt_structure'd); field names mirror the original Julia
# locals so the kernel body reads verbatim. On the host path the same bundles
# are built from CPU arrays, so CPU stays byte-identical.
# --------------------------------------------------------------------------

# Soil-state / geometry inputs (read-only) the search loops touch.
Base.@kwdef struct _WTInDV{V,M}
    col_dz::M; col_z::M; col_zi::M
    t_soisno::M; watsat::M; sucsat::M; bsw::M
    qcharge::V
    h2osoi_liq::M; h2osoi_ice::M
end
Adapt.@adapt_structure _WTInDV

# State outputs (written and/or read-modified) the search loops touch.
Base.@kwdef struct _WTOutDV{V,M}
    zwt::V; wa::V
    zwt_perched::V; frost_table::V; h2osoi_vol::M
    qflx_drain::V; qflx_rsub_sat::V; qflx_drain_perched::V
end
Adapt.@adapt_structure _WTOutDV

@kernel function _water_table_kernel!(zwt_out, odv, @Const(idv), @Const(mask),
        nlevsoi::Int, joff::Int, joff_zi::Int, dtime, aq_sp_yield_min, tfrz, denh2o, denice)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(zwt_out)
        aqy   = T(aq_sp_yield_min)
        e3    = T(1.0e3)
        thou  = T(1000.0)
        zr    = zero(T)
        on    = one(T)
        TFRZl = T(tfrz)
        DH2O  = T(denh2o)
        DICE  = T(denice)

        # ---- init drainage fluxes ----
        odv.qflx_drain[c]         = zr
        odv.qflx_rsub_sat[c]      = zr
        odv.qflx_drain_perched[c] = zr

        # ---- find jwt: index of soil layer right above the water table ----
        jwt = nlevsoi
        for j in 1:nlevsoi
            if odv.zwt[c] <= idv.col_zi[c, joff_zi + j]
                jwt = j - 1
                break
            end
        end

        # ============================ QCHARGE ============================
        # analytical expression for aquifer specific yield
        rous = idv.watsat[c, nlevsoi] *
            (on - (on + e3 * odv.zwt[c] / idv.sucsat[c, nlevsoi])^(-on / idv.bsw[c, nlevsoi]))
        rous = smooth_max(rous, aqy)

        if jwt == nlevsoi
            # water table below soil column
            odv.wa[c]  = odv.wa[c] + idv.qcharge[c] * dtime
            odv.zwt[c] = odv.zwt[c] - (idv.qcharge[c] * dtime) / thou / rous
        else
            # water table within soil layers
            qcharge_tot = idv.qcharge[c] * dtime
            if qcharge_tot > zr  # rising water table
                for j in (jwt+1):-1:1
                    s_y = idv.watsat[c, j] *
                        (on - (on + e3 * odv.zwt[c] / idv.sucsat[c, j])^(-on / idv.bsw[c, j]))
                    s_y = smooth_max(s_y, aqy)

                    qcharge_layer = smooth_min(qcharge_tot,
                        s_y * (odv.zwt[c] - idv.col_zi[c, joff_zi + j - 1]) * e3)
                    qcharge_layer = smooth_max(qcharge_layer, zr)

                    if s_y > zr
                        odv.zwt[c] = odv.zwt[c] - qcharge_layer / s_y / thou
                    end

                    qcharge_tot = qcharge_tot - qcharge_layer
                    if qcharge_tot <= zr
                        break
                    end
                end
            else  # deepening water table
                for j in (jwt+1):nlevsoi
                    s_y = idv.watsat[c, j] *
                        (on - (on + e3 * odv.zwt[c] / idv.sucsat[c, j])^(-on / idv.bsw[c, j]))
                    s_y = smooth_max(s_y, aqy)

                    qcharge_layer = smooth_max(qcharge_tot,
                        -(s_y * (idv.col_zi[c, joff_zi + j] - odv.zwt[c]) * e3))
                    qcharge_layer = smooth_min(qcharge_layer, zr)
                    qcharge_tot = qcharge_tot - qcharge_layer
                    if qcharge_tot >= zr
                        odv.zwt[c] = odv.zwt[c] - qcharge_layer / s_y / thou
                        break
                    else
                        odv.zwt[c] = idv.col_zi[c, joff_zi + j]
                    end
                end
                if qcharge_tot > zr
                    odv.zwt[c] = odv.zwt[c] - qcharge_tot / thou / rous
                end
            end

            # recompute jwt
            jwt = nlevsoi
            for j in 1:nlevsoi
                if odv.zwt[c] <= idv.col_zi[c, joff_zi + j]
                    jwt = j - 1
                    break
                end
            end
        end

        # ============================ BASEFLOW ============================
        # define frost table as first frozen layer with unfrozen layer above it
        if idv.t_soisno[c, joff + 1] > TFRZl
            k_frz = nlevsoi
        else
            k_frz = 1
        end
        for k in 2:nlevsoi
            if idv.t_soisno[c, joff + k - 1] > TFRZl && idv.t_soisno[c, joff + k] <= TFRZl
                k_frz = k
                break
            end
        end

        odv.frost_table[c] = idv.col_z[c, joff + k_frz]
        odv.zwt_perched[c] = odv.frost_table[c]

        if odv.zwt[c] < odv.frost_table[c] && idv.t_soisno[c, joff + k_frz] <= TFRZl
            # do nothing - handled in Drainage
        else
            sat_lev = T(0.9)

            k_perch = 1
            for k in k_frz:-1:1
                odv.h2osoi_vol[c, k] = idv.h2osoi_liq[c, joff + k] / (idv.col_dz[c, joff + k] * DH2O) +
                    idv.h2osoi_ice[c, joff + k] / (idv.col_dz[c, joff + k] * DICE)
                if odv.h2osoi_vol[c, k] / idv.watsat[c, k] <= sat_lev
                    k_perch = k
                    break
                end
            end

            if idv.t_soisno[c, joff + k_frz] > TFRZl
                k_perch = k_frz
            end

            if k_frz > k_perch
                s1 = (idv.h2osoi_liq[c, joff + k_perch] / (idv.col_dz[c, joff + k_perch] * DH2O) +
                    idv.h2osoi_ice[c, joff + k_perch] / (idv.col_dz[c, joff + k_perch] * DICE)) / idv.watsat[c, k_perch]
                s2 = (idv.h2osoi_liq[c, joff + k_perch + 1] / (idv.col_dz[c, joff + k_perch + 1] * DH2O) +
                    idv.h2osoi_ice[c, joff + k_perch + 1] / (idv.col_dz[c, joff + k_perch + 1] * DICE)) / idv.watsat[c, k_perch + 1]

                m_val = (idv.col_z[c, joff + k_perch + 1] - idv.col_z[c, joff + k_perch]) / (s2 - s1)
                b_val = idv.col_z[c, joff + k_perch + 1] - m_val * s2
                odv.zwt_perched[c] = smooth_max(zr, m_val * sat_lev + b_val)
            end
        end
    end
end

function soilhyd_water_table!(odv, idv, mask, nlevsoi::Int, joff::Int, joff_zi::Int,
                              dtime, aq_sp_yield_min, tfrz, denh2o, denice)
    # Convert all scalar reals to the working precision of the device arrays so no
    # Float64 reaches the Float32-only Metal backend (byte-identical on Float64 CPU).
    T = eltype(odv.zwt)
    _launch!(_water_table_kernel!, odv.zwt, odv, idv, mask, nlevsoi, joff, joff_zi,
             T(dtime), T(aq_sp_yield_min), T(tfrz), T(denh2o), T(denice);
             ndrange = length(mask))
end

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
    temperature_t_soisno::AbstractMatrix{<:Real},
    waterstatebulk_ws::WaterStateData,
    waterfluxbulk::WaterFluxBulkData,
    col_dz::AbstractMatrix{<:Real},
    col_z::AbstractMatrix{<:Real},
    col_zi::AbstractMatrix{<:Real},
    mask_hydrology::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Real
)
    wf = waterfluxbulk.wf
    params = soilhydrology_params

    # Offsets for snow+soil arrays (col_dz, col_z, col_zi, h2osoi_liq/ice, t_soisno)
    nlevsno  = varpar.nlevsno
    joff     = nlevsno          # for z, dz, h2osoi_liq/ice, t_soisno
    joff_zi  = nlevsno + 1     # for zi (has extra element at top)

    # Device-view bundles (read-only inputs + read/write outputs). Each is built
    # from the live state arrays so all writes flow straight back; on the host
    # these are CPU arrays and the result is byte-identical.
    idv = _WTInDV(;
        col_dz     = col_dz,
        col_z      = col_z,
        col_zi     = col_zi,
        t_soisno   = temperature_t_soisno,
        watsat     = soilstate.watsat_col,
        sucsat     = soilstate.sucsat_col,
        bsw        = soilstate.bsw_col,
        qcharge    = soilhydrology.qcharge_col,
        h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col,
        h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col)

    odv = _WTOutDV(;
        zwt                = soilhydrology.zwt_col,
        wa                 = waterstatebulk_ws.wa_col,
        zwt_perched        = soilhydrology.zwt_perched_col,
        frost_table        = soilhydrology.frost_table_col,
        h2osoi_vol         = waterstatebulk_ws.h2osoi_vol_col,
        qflx_drain         = wf.qflx_drain_col,
        qflx_rsub_sat      = wf.qflx_rsub_sat_col,
        qflx_drain_perched = wf.qflx_drain_perched_col)

    soilhyd_water_table!(odv, idv, mask_hydrology, nlevsoi, joff, joff_zi,
                         dtime, params.aq_sp_yield_min, TFRZ, DENH2O, DENICE)

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
    temperature_t_soisno::Matrix{<:Real},
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    waterstatebulk_ws::WaterStateData,
    waterfluxbulk::WaterFluxBulkData,
    col_dz::Matrix{<:Real},
    col_z::Matrix{<:Real},
    col_zi::Matrix{<:Real},
    col_snl::Vector{Int},
    col_itype::Vector{Int},
    col_landunit::Vector{Int},
    col_topo_slope::Vector{<:Real},
    lun_urbpoi::Vector{Bool},
    mask_hydrology::BitVector,
    mask_urban::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Real;
    nlevsno::Int = varpar.nlevsno,
    use_vichydro::Bool = false
)
    wf = waterfluxbulk.wf
    params = soilhydrology_params

    # Offsets for snow+soil arrays (h2osoi_liq/ice, t_soisno, dz, z, zi)
    # Soil layer j in Fortran maps to index j + joff in Julia 1-based arrays
    joff = nlevsno
    joff_zi = nlevsno + 1

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
    FT_dr = eltype(col_dz)
    rsub_top = Vector{FT_dr}(undef, nb)
    fff      = Vector{FT_dr}(undef, nb)
    xsi_arr  = Vector{FT_dr}(undef, nb)
    xs1_arr  = Vector{FT_dr}(undef, nb)
    xs_arr   = Vector{FT_dr}(undef, nb)

    dzmm = Matrix{FT_dr}(undef, nb, nlevsoi)

    # Convert layer thicknesses and compute icefrac
    for j in 1:nlevsoi
        for c in bounds
            mask_hydrology[c] || continue
            dzmm[c, j] = col_dz[c, j + joff] * 1.0e3

            vol_ice_val = smooth_min(watsat[c, j], h2osoi_ice[c, j + joff] / (col_dz[c, j + joff] * DENICE))
            icefrac[c, j] = smooth_min(1.0, vol_ice_val / watsat[c, j])
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
            if zwt[c] <= col_zi[c, j + joff_zi]
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
        if t_soisno[c, 1 + joff] > TFRZ
            k_frz = nlevsoi
        else
            k_frz = 1
        end

        for k in 2:nlevsoi
            if t_soisno[c, k-1 + joff] > TFRZ && t_soisno[c, k + joff] <= TFRZ
                k_frz = k
                break
            end
        end

        frost_table[c] = col_z[c, k_frz + joff]
        zwt_perched[c] = frost_table[c]
        qflx_drain_perched[c] = 0.0

        # ======= water table above frost table ========
        if zwt[c] < frost_table[c] && t_soisno[c, k_frz + joff] <= TFRZ
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
                rsub_top_layer = smooth_max(rsub_top_tot, -(h2osoi_liq[c, k + joff] - WATMIN))
                rsub_top_layer = smooth_min(rsub_top_layer, 0.0)
                rsub_top_tot = rsub_top_tot - rsub_top_layer

                h2osoi_liq[c, k + joff] = h2osoi_liq[c, k + joff] + rsub_top_layer

                if rsub_top_tot >= 0.0
                    zwt[c] = zwt[c] - rsub_top_layer / eff_porosity[c, k] / 1000.0
                    break
                else
                    zwt[c] = col_zi[c, k + joff_zi]
                end
            end

            qflx_drain_perched[c] = qflx_drain_perched[c] + rsub_top_tot / dtime

            # recompute jwt
            jwt[c] = nlevsoi
            for j in 1:nlevsoi
                if zwt[c] <= col_zi[c, j + joff_zi]
                    jwt[c] = j - 1
                    break
                end
            end
        else
            # ======= water table below frost table ========
            sat_lev = 0.9

            k_perch = 1
            for k in k_frz:-1:1
                h2osoi_vol_val = h2osoi_liq[c, k + joff] / (col_dz[c, k + joff] * DENH2O) +
                    h2osoi_ice[c, k + joff] / (col_dz[c, k + joff] * DENICE)

                if h2osoi_vol_val / watsat[c, k] <= sat_lev
                    k_perch = k
                    break
                end
            end

            if t_soisno[c, k_frz + joff] > TFRZ
                k_perch = k_frz
            end

            # if perched water table exists
            if k_frz > k_perch
                s1 = (h2osoi_liq[c, k_perch + joff] / (col_dz[c, k_perch + joff] * DENH2O) +
                    h2osoi_ice[c, k_perch + joff] / (col_dz[c, k_perch + joff] * DENICE)) / watsat[c, k_perch]
                s2 = (h2osoi_liq[c, k_perch+1 + joff] / (col_dz[c, k_perch+1 + joff] * DENH2O) +
                    h2osoi_ice[c, k_perch+1 + joff] / (col_dz[c, k_perch+1 + joff] * DENICE)) / watsat[c, k_perch+1]

                m_val = (col_z[c, k_perch+1 + joff] - col_z[c, k_perch + joff]) / (s2 - s1)
                b_val = col_z[c, k_perch+1 + joff] - m_val * s2
                zwt_perched[c] = smooth_max(0.0, m_val * sat_lev + b_val)

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
                    rsub_top_layer = smooth_max(rsub_top_tot, -(h2osoi_liq[c, k + joff] - WATMIN))
                    rsub_top_layer = smooth_min(rsub_top_layer, 0.0)
                    rsub_top_tot = rsub_top_tot - rsub_top_layer

                    h2osoi_liq[c, k + joff] = h2osoi_liq[c, k + joff] + rsub_top_layer

                    if rsub_top_tot >= 0.0
                        zwt_perched[c] = zwt_perched[c] - rsub_top_layer / eff_porosity[c, k] / 1000.0
                        break
                    else
                        zwt_perched[c] = col_zi[c, k + joff_zi]
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
            _rous_base = 1.0 + 1.0e3 * zwt[c] / sucsat[c, nlevsoi]
            _rous_exp = -1.0 / bsw[c, nlevsoi]
            rous_local = watsat[c, nlevsoi] * (1.0 - _rous_base^_rous_exp)
            rous_local = smooth_max(rous_local, params.aq_sp_yield_min)

            if jwt[c] == nlevsoi
                # water table below soil column
                wa[c] = wa[c] - rsub_top[c] * dtime
                zwt[c] = zwt[c] + (rsub_top[c] * dtime) / 1000.0 / rous_local
                h2osoi_liq[c, nlevsoi + joff] = h2osoi_liq[c, nlevsoi + joff] + smooth_max(0.0, wa[c] - aquifer_water_baseline)
                wa[c] = smooth_min(wa[c], aquifer_water_baseline)
            else
                # water table within soil layers
                rsub_top_tot = -rsub_top[c] * dtime

                if rsub_top_tot > 0.0
                    if _is_ad_type(eltype(rsub_top))
                        @warn "RSUB_TOP IS POSITIVE in Drainage (AD mode, continuing)" maxlog=1
                    else
                        error("RSUB_TOP IS POSITIVE in Drainage!")
                    end
                else
                    for j in (jwt[c]+1):nlevsoi
                        s_y = watsat[c, j] *
                            (1.0 - (1.0 + 1.0e3 * zwt[c] / sucsat[c, j])^(-1.0 / bsw[c, j]))
                        s_y = smooth_max(s_y, params.aq_sp_yield_min)

                        rsub_top_layer = smooth_max(rsub_top_tot, -(s_y * (col_zi[c, j + joff_zi] - zwt[c]) * 1.0e3))
                        rsub_top_layer = smooth_min(rsub_top_layer, 0.0)
                        h2osoi_liq[c, j + joff] = h2osoi_liq[c, j + joff] + rsub_top_layer

                        rsub_top_tot = rsub_top_tot - rsub_top_layer

                        if rsub_top_tot >= 0.0
                            zwt[c] = zwt[c] - rsub_top_layer / s_y / 1000.0
                            break
                        else
                            zwt[c] = col_zi[c, j + joff_zi]
                        end
                    end

                    # remove residual
                    zwt[c] = zwt[c] - rsub_top_tot / 1000.0 / rous_local
                    wa[c] = wa[c] + rsub_top_tot
                end

                # recompute jwt
                jwt[c] = nlevsoi
                for j in 1:nlevsoi
                    if zwt[c] <= col_zi[c, j + joff_zi]
                        jwt[c] = j - 1
                        break
                    end
                end
            end

            zwt[c] = smooth_max(0.0, zwt[c])
            zwt[c] = smooth_min(80.0, zwt[c])
        end
    end

    # Excessive water above saturation: redistribute upward
    for j in nlevsoi:-1:2
        for c in bounds
            mask_hydrology[c] || continue
            xsi_arr[c]             = smooth_max(h2osoi_liq[c, j + joff] - eff_porosity[c, j] * dzmm[c, j], 0.0)
            h2osoi_liq[c, j + joff]   = smooth_min(eff_porosity[c, j] * dzmm[c, j], h2osoi_liq[c, j + joff])
            h2osoi_liq[c, j-1 + joff] = h2osoi_liq[c, j-1 + joff] + xsi_arr[c]
        end
    end

    for c in bounds
        mask_hydrology[c] || continue

        # watmin addition to fix water balance errors
        xs1_arr[c] = smooth_max(smooth_max(h2osoi_liq[c, 1 + joff] - WATMIN, 0.0) -
            smooth_max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_ice[c, 1 + joff] - WATMIN), 0.0)
        h2osoi_liq[c, 1 + joff] = h2osoi_liq[c, 1 + joff] - xs1_arr[c]

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
        xs1_arr[c] = smooth_max(smooth_max(h2osoi_ice[c, 1 + joff], 0.0) -
            smooth_max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_liq[c, 1 + joff]), 0.0)
        h2osoi_ice[c, 1 + joff] = smooth_min(smooth_max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_liq[c, 1 + joff]),
            h2osoi_ice[c, 1 + joff])
        qflx_ice_runoff_xs[c] = xs1_arr[c] / dtime
    end

    # Limit h2osoi_liq >= watmin: pull from below
    for j in 1:(nlevsoi-1)
        for c in bounds
            mask_hydrology[c] || continue
            if h2osoi_liq[c, j + joff] < WATMIN
                xs_arr[c] = WATMIN - h2osoi_liq[c, j + joff]
                if j == jwt[c]
                    zwt[c] = zwt[c] + xs_arr[c] / eff_porosity[c, j] / 1000.0
                end
            else
                xs_arr[c] = 0.0
            end
            h2osoi_liq[c, j + joff]   = h2osoi_liq[c, j + joff]   + xs_arr[c]
            h2osoi_liq[c, j+1 + joff] = h2osoi_liq[c, j+1 + joff] - xs_arr[c]
        end
    end

    # Get water for bottom layer from layers above if possible
    j = nlevsoi
    for c in bounds
        mask_hydrology[c] || continue
        if h2osoi_liq[c, j + joff] < WATMIN
            xs_arr[c] = WATMIN - h2osoi_liq[c, j + joff]
            for i in (nlevsoi-1):-1:1
                available_h2osoi_liq = smooth_max(h2osoi_liq[c, i + joff] - WATMIN - xs_arr[c], 0.0)
                if available_h2osoi_liq >= xs_arr[c]
                    h2osoi_liq[c, j + joff] = h2osoi_liq[c, j + joff] + xs_arr[c]
                    h2osoi_liq[c, i + joff] = h2osoi_liq[c, i + joff] - xs_arr[c]
                    xs_arr[c] = 0.0
                    break
                else
                    h2osoi_liq[c, j + joff] = h2osoi_liq[c, j + joff] + available_h2osoi_liq
                    h2osoi_liq[c, i + joff] = h2osoi_liq[c, i + joff] - available_h2osoi_liq
                    xs_arr[c] = xs_arr[c] - available_h2osoi_liq
                end
            end
        else
            xs_arr[c] = 0.0
        end
        h2osoi_liq[c, j + joff] = h2osoi_liq[c, j + joff] + xs_arr[c]
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
    col_dz::Matrix{<:Real},
    col_zi::Matrix{<:Real},
    col_z::Matrix{<:Real},
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
            ice[c, i]       = smooth_min(moist0 + ice0, ice[c, i])
            ice[c, i]       = smooth_max(0.0, ice[c, i])
            moist[c, i]     = smooth_max(WATMIN, moist[c, i])
            moist[c, i]     = smooth_min(max_moist[c, i] - ice[c, i], moist[c, i])
            moist_vol[c, i] = moist[c, i] / (depth[c, i] * DENICE) + ice[c, i] / (depth[c, i] * DENH2O)
            moist_vol[c, i] = smooth_min(porosity[c, i], moist_vol[c, i])
            moist_vol[c, i] = smooth_max(0.01, moist_vol[c, i])
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
        top_moist_limited[c] = smooth_min(top_moist[c], top_max_moist[c])
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
    temperature_t_soisno::Matrix{<:Real},
    waterstatebulk_ws::WaterStateData,
    col_dz::Matrix{<:Real},
    col_z::Matrix{<:Real},
    col_zi::Matrix{<:Real},
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
                    zwt_perched[c] = smooth_max(0.0, m_val * sat_lev + b_val)
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
    col_dz::Matrix{<:Real},
    col_z::Matrix{<:Real},
    col_zi::Matrix{<:Real},
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
            zwt[c] = smooth_max(0.0, m_val * sat_lev + b_val)
        else
            zwt[c] = col_zi[c, col_nbedrock[c]]
        end
    end

    return nothing
end

# =========================================================================
# RenewCondensation
# =========================================================================

# ---- renew_condensation! : per-column condensation/sublimation on top layer ----
# Hydrologically-active columns: weight dew/evap by (1 - frac_h2osfc).
@kernel function _soilhyd_renew_cond_kernel!(h2osoi_liq, h2osoi_ice, @Const(mask),
                                             @Const(col_snl), @Const(frac_h2osfc),
                                             @Const(qflx_liqdew_to_top_layer),
                                             @Const(qflx_soliddew_to_top_layer),
                                             @Const(qflx_solidevap_from_top_layer),
                                             nlevsno::Int, tol, dtime)
    c = @index(Global)
    @inbounds if mask[c]
        if col_snl[c] + 1 >= 1
            T = eltype(h2osoi_liq)
            dt = T(dtime)
            tl = T(tol)
            jj = col_snl[c] + 1 + nlevsno
            wgt = one(T) - frac_h2osfc[c]
            h2osoi_liq[c, jj] = h2osoi_liq[c, jj] + wgt * qflx_liqdew_to_top_layer[c] * dt
            h2osoi_ice[c, jj] = h2osoi_ice[c, jj] + wgt * qflx_soliddew_to_top_layer[c] * dt
            h2osoi_ice_before = h2osoi_ice[c, jj]
            h2osoi_ice[c, jj] = h2osoi_ice[c, jj] - wgt * qflx_solidevap_from_top_layer[c] * dt

            if h2osoi_ice[c, jj] < zero(T) && abs(h2osoi_ice[c, jj]) < tl * abs(h2osoi_ice_before)
                h2osoi_ice[c, jj] = zero(T)
            end
            if h2osoi_ice[c, jj] < zero(T)
                h2osoi_ice[c, jj] = zero(T)
            end
        end
    end
end

soilhyd_renew_cond!(h2osoi_liq, h2osoi_ice, mask, col_snl, frac_h2osfc,
                    qflx_liqdew_to_top_layer, qflx_soliddew_to_top_layer,
                    qflx_solidevap_from_top_layer, nlevsno::Int, tol, dtime) =
    _launch!(_soilhyd_renew_cond_kernel!, h2osoi_liq, h2osoi_ice, mask, col_snl,
             frac_h2osfc, qflx_liqdew_to_top_layer, qflx_soliddew_to_top_layer,
             qflx_solidevap_from_top_layer, nlevsno,
             convert(eltype(h2osoi_liq), tol), convert(eltype(h2osoi_liq), dtime);
             ndrange = length(mask))

# ---- renew_condensation! : per-column condensation on urban roof/road ----
# Unweighted (frac_h2osfc not applied) for impervious urban surfaces.
@kernel function _soilhyd_renew_cond_urban_kernel!(h2osoi_liq, h2osoi_ice, @Const(mask_urban),
                                                   @Const(col_itype), @Const(col_snl),
                                                   @Const(qflx_liqdew_to_top_layer),
                                                   @Const(qflx_soliddew_to_top_layer),
                                                   @Const(qflx_solidevap_from_top_layer),
                                                   icol_roof::Int, icol_road_imperv::Int,
                                                   nlevsno::Int, tol, dtime)
    c = @index(Global)
    @inbounds if mask_urban[c]
        if col_itype[c] == icol_roof || col_itype[c] == icol_road_imperv
            if col_snl[c] + 1 >= 1
                T = eltype(h2osoi_liq)
                dt = T(dtime)
                tl = T(tol)
                jj = col_snl[c] + 1 + nlevsno
                h2osoi_liq[c, jj] = h2osoi_liq[c, jj] + qflx_liqdew_to_top_layer[c] * dt
                h2osoi_ice[c, jj] = h2osoi_ice[c, jj] + qflx_soliddew_to_top_layer[c] * dt
                h2osoi_ice_before = h2osoi_ice[c, jj]
                h2osoi_ice[c, jj] = h2osoi_ice[c, jj] - qflx_solidevap_from_top_layer[c] * dt

                if h2osoi_ice[c, jj] < zero(T) && abs(h2osoi_ice[c, jj]) < tl * abs(h2osoi_ice_before)
                    h2osoi_ice[c, jj] = zero(T)
                end
                if h2osoi_ice[c, jj] < zero(T)
                    h2osoi_ice[c, jj] = zero(T)
                end
            end
        end
    end
end

soilhyd_renew_cond_urban!(h2osoi_liq, h2osoi_ice, mask_urban, col_itype, col_snl,
                          qflx_liqdew_to_top_layer, qflx_soliddew_to_top_layer,
                          qflx_solidevap_from_top_layer, icol_roof::Int,
                          icol_road_imperv::Int, nlevsno::Int, tol, dtime) =
    _launch!(_soilhyd_renew_cond_urban_kernel!, h2osoi_liq, h2osoi_ice, mask_urban,
             col_itype, col_snl, qflx_liqdew_to_top_layer, qflx_soliddew_to_top_layer,
             qflx_solidevap_from_top_layer, icol_roof, icol_road_imperv, nlevsno,
             convert(eltype(h2osoi_liq), tol), convert(eltype(h2osoi_liq), dtime);
             ndrange = length(mask_urban))

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
    col_snl::AbstractVector{<:Integer},
    col_itype::AbstractVector{<:Integer},
    mask_hydrology::AbstractVector{Bool},
    mask_urban::AbstractVector{Bool},
    bounds::UnitRange{Int},
    dtime::Real
)
    wf = waterfluxbulk.wf

    h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col
    frac_h2osfc = waterdiagbulk.frac_h2osfc_col
    qflx_liqdew_to_top_layer      = wf.qflx_liqdew_to_top_layer_col
    qflx_soliddew_to_top_layer    = wf.qflx_soliddew_to_top_layer_col
    qflx_solidevap_from_top_layer = wf.qflx_solidevap_from_top_layer_col

    # Fortran layer index: snl(c)+1 → Julia: snl(c)+1+nlevsno
    nlevsno = varpar.nlevsno

    soilhyd_renew_cond!(h2osoi_liq, h2osoi_ice, mask_hydrology, col_snl, frac_h2osfc,
                        qflx_liqdew_to_top_layer, qflx_soliddew_to_top_layer,
                        qflx_solidevap_from_top_layer, nlevsno, TOLERANCE_SOILHYDRO, dtime)

    # Urban columns
    soilhyd_renew_cond_urban!(h2osoi_liq, h2osoi_ice, mask_urban, col_itype, col_snl,
                              qflx_liqdew_to_top_layer, qflx_soliddew_to_top_layer,
                              qflx_solidevap_from_top_layer, ICOL_ROOF, ICOL_ROAD_IMPERV,
                              nlevsno, TOLERANCE_SOILHYDRO, dtime)

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
    qflx_gw_demand::Vector{<:Real},
    qflx_gw_uncon_irrig_lyr::Matrix{<:Real},
    qflx_gw_con_irrig::Vector{<:Real},
    col_nbedrock::Vector{Int},
    col_zi::Matrix{<:Real},
    mask_soil::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Real
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
            s_y = smooth_max(s_y, params.aq_sp_yield_min)

            if j == jwt_c + 1
                available_water_layer = smooth_max(0.0, s_y * (col_zi[c, j] - zwt[c]) * 1.0e3)
            else
                available_water_layer = smooth_max(0.0, s_y * (col_zi[c, j] - col_zi[c, j-1]) * 1.0e3)
            end

            irrig_layer = smooth_min(irrig_demand_remaining, available_water_layer)
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
    dtime::Real
)
    qflx_gw_uncon_irrig_lyr = waterflux_wf.qflx_gw_uncon_irrig_lyr_col
    qflx_gw_con_irrig       = waterflux_wf.qflx_gw_con_irrig_col
    wa                       = waterstate_ws.wa_col
    h2osoi_liq               = waterstate_ws.h2osoi_liq_col

    soilhyd_withdraw_gw_lyr!(h2osoi_liq, mask_soil, qflx_gw_uncon_irrig_lyr, nlevsoi, dtime)

    soilhyd_withdraw_gw_con!(wa, mask_soil, qflx_gw_con_irrig, dtime)

    return nothing
end

# =========================================================================
# PerchedLateralFlow
# =========================================================================

"""
    perched_lateral_flow!(soilhydrology, soilstate, waterstatebulk_ws,
                           waterfluxbulk, col_data, lun_data,
                           tdepth_grc, tdepthmax_grc,
                           mask_hydrology, bounds, nlevsoi, dtime;
                           use_hillslope_routing=false)

Calculate subsurface drainage from perched saturated zone.

Ported from `PerchedLateralFlow` in `SoilHydrologyMod.F90`.
"""
function perched_lateral_flow!(
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    waterstatebulk_ws::WaterStateData,
    waterfluxbulk::WaterFluxBulkData,
    col_data::ColumnData,
    lun_data::LandunitData,
    tdepth_grc::Vector{<:Real},
    tdepthmax_grc::Vector{<:Real},
    mask_hydrology::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Real;
    use_hillslope_routing::Bool = false
)
    wf = waterfluxbulk.wf
    params = soilhydrology_params
    nlevsno = varpar.nlevsno
    joff = nlevsno
    joff_zi = nlevsno + 1

    h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col
    stream_water_volume = waterstatebulk_ws.stream_water_volume_lun

    bsw    = soilstate.bsw_col
    hksat  = soilstate.hksat_col
    sucsat = soilstate.sucsat_col
    watsat = soilstate.watsat_col

    frost_table  = soilhydrology.frost_table_col
    zwt          = soilhydrology.zwt_col
    zwt_perched  = soilhydrology.zwt_perched_col

    qflx_drain_perched = wf.qflx_drain_perched_col

    col_dz          = col_data.dz
    col_z           = col_data.z
    col_zi          = col_data.zi
    col_nbedrock    = col_data.nbedrock
    col_landunit    = col_data.landunit
    col_gridcell    = col_data.gridcell
    col_is_hillslope = col_data.is_hillslope_column
    col_active      = col_data.active
    col_cold        = col_data.cold
    col_hill_slope  = col_data.hill_slope
    col_hill_elev   = col_data.hill_elev
    col_hill_distance = col_data.hill_distance
    col_hill_width  = col_data.hill_width
    col_hill_area   = col_data.hill_area
    col_topo_slope  = col_data.topo_slope

    lun_stream_channel_length = lun_data.stream_channel_length
    lun_stream_channel_width  = lun_data.stream_channel_width
    lun_stream_channel_depth  = lun_data.stream_channel_depth

    k_anisotropic = 1.0

    nb = last(bounds)

    # Local arrays
    k_frost_arr = Vector{Int}(undef, nb)
    k_perch_arr = Vector{Int}(undef, nb)
    FT_pl = eltype(col_data.dz)
    qflx_drain_perched_vol = Vector{FT_pl}(undef, nb)
    qflx_drain_perched_out = Vector{FT_pl}(undef, nb)

    # ---- Locate frost table and perched water table layers ----
    for c in bounds
        mask_hydrology[c] || continue
        k_frost_arr[c] = col_nbedrock[c]
        k_perch_arr[c] = col_nbedrock[c]

        for k in 1:col_nbedrock[c]
            # Fortran zi(c,0) = 0.0 (ground surface); handle k-1=0 case
            zi_km1 = k > 1 ? col_zi[c, k-1 + joff_zi] : 0.0
            if frost_table[c] >= zi_km1 && frost_table[c] < col_zi[c, k + joff_zi]
                k_frost_arr[c] = k
                break
            end
        end

        for k in 1:col_nbedrock[c]
            zi_km1 = k > 1 ? col_zi[c, k-1 + joff_zi] : 0.0
            if zwt_perched[c] >= zi_km1 && zwt_perched[c] < col_zi[c, k + joff_zi]
                k_perch_arr[c] = k
                break
            end
        end
    end

    # ---- Compute drainage from perched saturated region ----
    for c in bounds
        mask_hydrology[c] || continue
        l = col_landunit[c]
        g = col_gridcell[c]
        qflx_drain_perched[c]     = 0.0
        qflx_drain_perched_out[c] = 0.0
        qflx_drain_perched_vol[c] = 0.0

        if frost_table[c] > zwt_perched[c]
            # Hillslope columns
            if col_is_hillslope[c] && col_active[c]
                # Calculate head gradient
                if HEAD_GRADIENT_METHOD[] == HEAD_GRADIENT_KINEMATIC
                    head_gradient = col_hill_slope[c]
                elseif HEAD_GRADIENT_METHOD[] == HEAD_GRADIENT_DARCY
                    if col_cold[c] != ISPVAL
                        head_gradient = (col_hill_elev[c] - zwt_perched[c]) -
                            (col_hill_elev[col_cold[c]] - zwt_perched[col_cold[c]])
                        head_gradient = head_gradient /
                            (col_hill_distance[c] - col_hill_distance[col_cold[c]])
                    else
                        if use_hillslope_routing
                            stream_water_depth = stream_water_volume[l] /
                                lun_stream_channel_length[l] / lun_stream_channel_width[l]
                            stream_channel_depth = lun_stream_channel_depth[l]
                        else
                            stream_water_depth = tdepth_grc[g]
                            stream_channel_depth = tdepthmax_grc[g]
                        end

                        head_gradient = (col_hill_elev[c] - zwt_perched[c]) -
                            smooth_max(smooth_min(stream_water_depth - stream_channel_depth, 0.0),
                                col_hill_elev[c] - frost_table[c])
                        head_gradient = head_gradient / col_hill_distance[c]

                        if stream_water_depth <= 0.0
                            head_gradient = smooth_max(head_gradient, 0.0)
                        end
                    end
                else
                    error("head_gradient_method must be kinematic or darcy")
                end

                # Determine source and destination columns
                if head_gradient >= 0.0
                    c_src = c
                    c_dst = col_cold[c]
                else
                    c_src = col_cold[c]
                    c_dst = c
                end

                # Calculate transmissivity of source column
                transmis = 0.0

                if TRANSMISSIVITY_METHOD[] == TRANSMISSIVITY_LAYERSUM
                    if HEAD_GRADIENT_METHOD[] == HEAD_GRADIENT_KINEMATIC
                        if k_perch_arr[c_src] < k_frost_arr[c_src]
                            for k in k_perch_arr[c_src]:(k_frost_arr[c_src]-1)
                                if k == k_perch_arr[c_src]
                                    transmis += 1.0e-3 * hksat[c_src, k] *
                                        (col_zi[c_src, k + joff_zi] - zwt_perched[c_src])
                                else
                                    transmis += 1.0e-3 * hksat[c_src, k] * col_dz[c_src, k + joff]
                                end
                            end
                        end
                    elseif HEAD_GRADIENT_METHOD[] == HEAD_GRADIENT_DARCY
                        if c_src == ISPVAL
                            transmis = 1.0e-3 * hksat[c, k_perch_arr[c_dst]] * stream_water_depth
                        else
                            if k_perch_arr[c_src] < k_frost_arr[c_src]
                                for k in k_perch_arr[c_src]:(k_frost_arr[c_src]-1)
                                    if k == k_perch_arr[c_src]
                                        transmis += 1.0e-3 * hksat[c_src, k] *
                                            (col_zi[c_src, k + joff_zi] - zwt_perched[c_src])
                                    else
                                        if c_dst == ISPVAL
                                            if (col_hill_elev[c_src] - col_z[c_src, k + joff]) > (-stream_channel_depth)
                                                transmis += 1.0e-3 * hksat[c_src, k] * col_dz[c_src, k + joff]
                                            end
                                        else
                                            if (col_hill_elev[c_src] - col_z[c_src, k + joff]) >
                                                    (col_hill_elev[c_dst] - zwt_perched[c_dst])
                                                transmis += 1.0e-3 * hksat[c_src, k] * col_dz[c_src, k + joff]
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                elseif TRANSMISSIVITY_METHOD[] == TRANSMISSIVITY_UNIFORM
                    transmis = 1.0e-3 * hksat[c_src, k_perch_arr[c_src]] *
                        (col_zi[c_src, k_frost_arr[c_src] + joff_zi] - zwt_perched[c_src])
                end

                transmis = k_anisotropic * transmis

                qflx_drain_perched_vol[c] = transmis * col_hill_width[c] * head_gradient
                qflx_drain_perched_out[c] = 1.0e3 * qflx_drain_perched_vol[c] / col_hill_area[c]

            else
                # Non-hillslope columns
                q_perch_max = params.perched_baseflow_scalar *
                    sin(col_topo_slope[c] * (RPI / 180.0))

                wtsub = 0.0
                q_perch = 0.0
                for k in k_perch_arr[c]:(k_frost_arr[c]-1)
                    q_perch += hksat[c, k] * col_dz[c, k + joff]
                    wtsub += col_dz[c, k + joff]
                end
                if wtsub > 0.0
                    q_perch = q_perch / wtsub
                end

                qflx_drain_perched_out[c] = q_perch_max * q_perch *
                    (frost_table[c] - zwt_perched[c])
            end
        end
    end

    # ---- Compute net drainage from perched saturated region ----
    for c in bounds
        mask_hydrology[c] || continue
        qflx_drain_perched[c] += qflx_drain_perched_out[c]
        if col_is_hillslope[c] && col_active[c]
            if col_cold[c] != ISPVAL
                qflx_drain_perched[col_cold[c]] -= 1.0e3 *
                    qflx_drain_perched_vol[c] / col_hill_area[col_cold[c]]
            end
        end
    end

    # ---- Remove drainage from soil moisture storage ----
    for c in bounds
        mask_hydrology[c] || continue

        drainage_tot = qflx_drain_perched[c] * dtime
        for k in k_perch_arr[c]:(k_frost_arr[c]-1)
            s_y = watsat[c, k] *
                (1.0 - (1.0 + 1.0e3 * zwt_perched[c] / sucsat[c, k])^(-1.0 / bsw[c, k]))
            s_y = smooth_max(s_y, params.aq_sp_yield_min)

            if k == k_perch_arr[c]
                drainage_layer = smooth_min(drainage_tot, s_y * (col_zi[c, k + joff_zi] - zwt_perched[c]) * 1.0e3)
            else
                drainage_layer = smooth_min(drainage_tot, s_y * col_dz[c, k + joff] * 1.0e3)
            end

            drainage_layer = smooth_max(drainage_layer, 0.0)
            drainage_tot -= drainage_layer
            h2osoi_liq[c, k + joff] -= drainage_layer
        end

        qflx_drain_perched[c] -= drainage_tot / dtime
    end

    return nothing
end

# =========================================================================
# SubsurfaceLateralFlow
# =========================================================================

"""
    subsurface_lateral_flow!(soilhydrology, soilstate, waterstatebulk_ws,
                              waterfluxbulk, col_data, lun_data,
                              tdepth_grc, tdepthmax_grc, grc_area,
                              mask_hydrology, mask_urban, bounds,
                              nlevsoi, dtime;
                              use_hillslope_routing=false,
                              nhillslope=0)

Calculate subsurface lateral flow from saturated zone.

Ported from `SubsurfaceLateralFlow` in `SoilHydrologyMod.F90`.
"""
function subsurface_lateral_flow!(
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    waterstatebulk_ws::WaterStateData,
    waterfluxbulk::WaterFluxBulkData,
    col_data::ColumnData,
    lun_data::LandunitData,
    tdepth_grc::Vector{<:Real},
    tdepthmax_grc::Vector{<:Real},
    grc_area::Vector{<:Real},
    mask_hydrology::BitVector,
    mask_urban::BitVector,
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Real;
    use_hillslope_routing::Bool = false,
    nhillslope::Int = 0
)
    wf = waterfluxbulk.wf
    params = soilhydrology_params
    nlevsno = varpar.nlevsno
    joff = nlevsno
    joff_zi = nlevsno + 1

    h2osfc     = waterstatebulk_ws.h2osfc_col
    h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col
    h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col
    stream_water_volume = waterstatebulk_ws.stream_water_volume_lun

    bsw          = soilstate.bsw_col
    hksat        = soilstate.hksat_col
    sucsat       = soilstate.sucsat_col
    watsat       = soilstate.watsat_col
    eff_porosity = soilstate.eff_porosity_col
    hk_l         = soilstate.hk_l_col

    icefrac      = soilhydrology.icefrac_col
    frost_table  = soilhydrology.frost_table_col
    zwt          = soilhydrology.zwt_col

    qflx_latflow_out   = wf.qflx_latflow_out_col
    qflx_latflow_in    = wf.qflx_latflow_in_col
    volumetric_discharge = wf.volumetric_discharge_col
    qflx_snwcp_liq     = wf.qflx_snwcp_liq_col
    qflx_ice_runoff_xs = wf.qflx_ice_runoff_xs_col
    qflx_drain         = wf.qflx_drain_col
    qflx_qrgwl         = wf.qflx_qrgwl_col
    qflx_rsub_sat      = wf.qflx_rsub_sat_col

    col_dz           = col_data.dz
    col_z            = col_data.z
    col_zi           = col_data.zi
    col_itype        = col_data.itype
    col_nbedrock     = col_data.nbedrock
    col_landunit     = col_data.landunit
    col_gridcell     = col_data.gridcell
    col_is_hillslope = col_data.is_hillslope_column
    col_active       = col_data.active
    col_cold         = col_data.cold
    col_colu         = col_data.colu
    col_hill_slope   = col_data.hill_slope
    col_hill_elev    = col_data.hill_elev
    col_hill_distance = col_data.hill_distance
    col_hill_width   = col_data.hill_width
    col_hill_area    = col_data.hill_area
    col_topo_slope   = col_data.topo_slope
    col_wtgcell      = col_data.wtgcell

    lun_stream_channel_length = lun_data.stream_channel_length
    lun_stream_channel_width  = lun_data.stream_channel_width
    lun_stream_channel_depth  = lun_data.stream_channel_depth
    lun_stream_channel_number = lun_data.stream_channel_number
    lun_urbpoi = lun_data.urbpoi

    k_anisotropic = 1.0

    nb = last(bounds)

    # Local arrays
    FT_sl = eltype(col_data.dz)
    jwt      = Vector{Int}(undef, nb)
    drainage_arr = Vector{FT_sl}(undef, nb)
    xsi_arr  = Vector{FT_sl}(undef, nb)
    xs1_arr  = Vector{FT_sl}(undef, nb)
    xs_arr   = Vector{FT_sl}(undef, nb)
    ice_imped_col_arr = Vector{FT_sl}(undef, nb)
    ice_imped_arr = Matrix{FT_sl}(undef, nb, nlevsoi)
    qflx_latflow_out_vol = Vector{FT_sl}(undef, nb)
    qflx_net_latflow     = Vector{FT_sl}(undef, nb)

    dzmm = Matrix{FT_sl}(undef, nb, nlevsoi)

    # ---- Convert layer thicknesses and compute icefrac/impedance ----
    for j in 1:nlevsoi
        for c in bounds
            mask_hydrology[c] || continue
            dzmm[c, j] = col_dz[c, j + joff] * 1.0e3

            vol_ice_val = smooth_min(watsat[c, j], h2osoi_ice[c, j + joff] / (col_dz[c, j + joff] * DENICE))
            icefrac[c, j] = smooth_min(1.0, vol_ice_val / watsat[c, j])
            ice_imped_arr[c, j] = 10.0^(-params.e_ice * icefrac[c, j])
        end
    end

    # ---- Initial set ----
    for c in bounds
        mask_hydrology[c] || continue
        qflx_drain[c]       = 0.0
        qflx_rsub_sat[c]    = 0.0
        drainage_arr[c]      = 0.0
        qflx_latflow_in[c]  = 0.0
        qflx_latflow_out[c] = 0.0
        qflx_net_latflow[c] = 0.0
        volumetric_discharge[c] = 0.0
        qflx_latflow_out_vol[c] = 0.0
    end

    # ---- jwt: layer index right above water table ----
    for c in bounds
        mask_hydrology[c] || continue
        jwt[c] = nlevsoi
        for j in 1:nlevsoi
            if zwt[c] <= col_zi[c, j + joff_zi]
                jwt[c] = j - 1
                break
            end
        end
    end

    # ---- Calculate column-average ice impedance factor ----
    for c in bounds
        mask_hydrology[c] || continue
        dzsum = 0.0
        icefracsum = 0.0
        for j in max(jwt[c], 1):nlevsoi
            dzsum += dzmm[c, j]
            icefracsum += icefrac[c, j] * dzmm[c, j]
        end
        ice_imped_col_arr[c] = 10.0^(-params.e_ice * (icefracsum / dzsum))
    end

    # ---- Compute lateral flow ----
    for c in bounds
        mask_hydrology[c] || continue
        l = col_landunit[c]
        g = col_gridcell[c]

        # Hillslope columns
        if col_is_hillslope[c] && col_active[c]
            # Head gradient method
            if HEAD_GRADIENT_METHOD[] == HEAD_GRADIENT_KINEMATIC
                head_gradient = col_hill_slope[c]
            elseif HEAD_GRADIENT_METHOD[] == HEAD_GRADIENT_DARCY
                if col_cold[c] != ISPVAL
                    head_gradient = (col_hill_elev[c] - zwt[c]) -
                        (col_hill_elev[col_cold[c]] - zwt[col_cold[c]])
                    head_gradient /= (col_hill_distance[c] - col_hill_distance[col_cold[c]])
                else
                    if use_hillslope_routing
                        stream_water_depth = stream_water_volume[l] /
                            lun_stream_channel_length[l] / lun_stream_channel_width[l]
                        stream_channel_depth = lun_stream_channel_depth[l]
                    else
                        stream_water_depth = tdepth_grc[g]
                        stream_channel_depth = tdepthmax_grc[g]
                    end

                    head_gradient = (col_hill_elev[c] - zwt[c]) -
                        smooth_min(stream_water_depth - stream_channel_depth, 0.0)
                    head_gradient /= col_hill_distance[c]

                    if stream_water_depth <= 0.0
                        head_gradient = smooth_max(head_gradient, 0.0)
                    end

                    if head_gradient < 0.0
                        head_gradient -= 1.0 / k_anisotropic
                    end
                end
            else
                error("head_gradient_method must be kinematic or darcy")
            end

            # Cap maximum head_gradient
            head_gradient = smooth_min(smooth_max(head_gradient, -2.0), 2.0)

            # Determine source and destination columns
            if head_gradient >= 0.0
                c_src = c
                c_dst = col_cold[c]
            else
                c_src = col_cold[c]
                c_dst = c
            end

            # Calculate transmissivity of source column
            transmis = 0.0
            if c_src != ISPVAL
                if zwt[c_src] <= col_zi[c_src, col_nbedrock[c_src] + joff_zi]
                    if TRANSMISSIVITY_METHOD[] == TRANSMISSIVITY_LAYERSUM
                        for j in (jwt[c_src]+1):col_nbedrock[c_src]
                            if j == jwt[c_src] + 1
                                transmis += 1.0e-3 * ice_imped_arr[c_src, j] *
                                    hksat[c_src, j] * (col_zi[c_src, j + joff_zi] - zwt[c_src])
                            else
                                if c_dst == ISPVAL
                                    if (col_hill_elev[c_src] - col_z[c_src, j + joff]) > (-stream_channel_depth)
                                        transmis += 1.0e-3 * ice_imped_arr[c_src, j] *
                                            hksat[c_src, j] * col_dz[c_src, j + joff]
                                    end
                                else
                                    if (col_hill_elev[c_src] - col_z[c_src, j + joff]) >
                                            (col_hill_elev[c_dst] - zwt[c_dst])
                                        transmis += 1.0e-3 * ice_imped_arr[c_src, j] *
                                            hksat[c_src, j] * col_dz[c_src, j + joff]
                                    end
                                end
                            end
                        end
                    elseif TRANSMISSIVITY_METHOD[] == TRANSMISSIVITY_UNIFORM
                        transmis = 1.0e-3 * ice_imped_arr[c_src, jwt[c_src]+1] *
                            hksat[c_src, jwt[c_src]+1] *
                            (col_zi[c_src, col_nbedrock[c_src] + joff_zi] - zwt[c_src])
                    else
                        error("transmissivity_method must be LayerSum or Uniform")
                    end
                end
            else
                # transmissivity of losing stream
                transmis = 1.0e-3 * ice_imped_arr[c, jwt[c]+1] *
                    hksat[c, jwt[c]+1] * stream_water_depth
            end

            transmis = k_anisotropic * transmis

            qflx_latflow_out_vol[c] = transmis * col_hill_width[c] * head_gradient

            # Limit outflow by available stream channel water
            if use_hillslope_routing && qflx_latflow_out_vol[c] < 0.0
                available_stream_water = stream_water_volume[l] /
                    lun_stream_channel_number[l] / nhillslope
                if smooth_abs(qflx_latflow_out_vol[c]) * dtime > available_stream_water
                    qflx_latflow_out_vol[c] = -available_stream_water / dtime
                end
            end

            # volumetric_discharge from lowest column
            if col_cold[c] == ISPVAL
                volumetric_discharge[c] = qflx_latflow_out_vol[c] *
                    (grc_area[g] * 1.0e6 * col_wtgcell[c] / col_hill_area[c])
            end

            # convert volumetric flow to equivalent flux
            qflx_latflow_out[c] = 1.0e3 * qflx_latflow_out_vol[c] / col_hill_area[c]

            # hilltop column has no inflow
            if col_colu[c] == ISPVAL
                qflx_latflow_in[c] = 0.0
            end

            # current outflow is inflow to downhill column
            if col_cold[c] != ISPVAL
                qflx_latflow_in[col_cold[c]] += 1.0e3 *
                    qflx_latflow_out_vol[c] / col_hill_area[col_cold[c]]
            end

        else
            # Non-hillslope columns: power law baseflow
            if zwt[c] <= col_zi[c, col_nbedrock[c] + joff_zi]
                qflx_latflow_out[c] = ice_imped_col_arr[c] * BASEFLOW_SCALAR[] *
                    tan(RPI / 180.0 * col_topo_slope[c]) *
                    (col_zi[c, col_nbedrock[c] + joff_zi] - zwt[c])^params.n_baseflow
            end
            qflx_latflow_out_vol[c] = 1.0e-3 * qflx_latflow_out[c] *
                (grc_area[g] * 1.0e6 * col_wtgcell[c])
            volumetric_discharge[c] = qflx_latflow_out_vol[c]
        end
    end

    # ---- Topographic runoff: remove water via drainage ----
    for c in bounds
        mask_hydrology[c] || continue

        qflx_net_latflow[c] = qflx_latflow_out[c] - qflx_latflow_in[c]

        if zwt[c] <= col_zi[c, col_nbedrock[c] + joff_zi]
            drainage_arr[c] = qflx_net_latflow[c]
        else
            drainage_arr[c] = 0.0
        end

        drainage_tot = -drainage_arr[c] * dtime

        if drainage_tot > 0.0  # rising water table
            for j in (jwt[c]+1):-1:1
                if col_zi[c, j + joff_zi] < frost_table[c]
                    s_y = watsat[c, j] *
                        (1.0 - (1.0 + 1.0e3 * zwt[c] / sucsat[c, j])^(-1.0 / bsw[c, j]))
                    s_y = smooth_max(s_y, params.aq_sp_yield_min)

                    drainage_layer = smooth_min(drainage_tot, s_y * col_dz[c, j + joff] * 1.0e3)
                    drainage_layer = smooth_max(drainage_layer, 0.0)
                    h2osoi_liq[c, j + joff] += drainage_layer

                    drainage_tot -= drainage_layer

                    if drainage_tot <= 0.0
                        zwt[c] -= drainage_layer / s_y / 1000.0
                        break
                    else
                        # Fortran zi(c,0) = 0.0 (ground surface); handle j-1=0 case
                        zwt[c] = j > 1 ? col_zi[c, j - 1 + joff_zi] : 0.0
                    end
                end
            end

            # remove residual
            h2osfc[c] += drainage_tot

        else  # deepening water table
            for j in (jwt[c]+1):col_nbedrock[c]
                s_y = watsat[c, j] *
                    (1.0 - (1.0 + 1.0e3 * zwt[c] / sucsat[c, j])^(-1.0 / bsw[c, j]))
                s_y = smooth_max(s_y, params.aq_sp_yield_min)

                drainage_layer = smooth_max(drainage_tot, -(s_y * (col_zi[c, j + joff_zi] - zwt[c]) * 1.0e3))
                drainage_layer = smooth_min(drainage_layer, 0.0)
                h2osoi_liq[c, j + joff] += drainage_layer

                drainage_tot -= drainage_layer

                if drainage_tot >= 0.0
                    zwt[c] -= drainage_layer / s_y / 1000.0
                    break
                else
                    zwt[c] = col_zi[c, j + joff_zi]
                end
            end

            # remove residual
            drainage_arr[c] += drainage_tot / dtime
        end

        zwt[c] = smooth_max(0.0, zwt[c])
        zwt[c] = smooth_min(80.0, zwt[c])
    end

    # ---- Excessive water above saturation: redistribute upward ----
    for j in nlevsoi:-1:2
        for c in bounds
            mask_hydrology[c] || continue
            xsi_arr[c]        = smooth_max(h2osoi_liq[c, j + joff] - eff_porosity[c, j] * dzmm[c, j], 0.0)
            h2osoi_liq[c, j + joff]  = smooth_min(eff_porosity[c, j] * dzmm[c, j], h2osoi_liq[c, j + joff])
            h2osoi_liq[c, j - 1 + joff] += xsi_arr[c]
        end
    end

    for c in bounds
        mask_hydrology[c] || continue

        xs1_arr[c] = smooth_max(smooth_max(h2osoi_liq[c, 1 + joff] - WATMIN, 0.0) -
            smooth_max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_ice[c, 1 + joff] - WATMIN), 0.0)
        h2osoi_liq[c, 1 + joff] -= xs1_arr[c]

        l = col_landunit[c]
        if lun_urbpoi[l]
            qflx_rsub_sat[c] = xs1_arr[c] / dtime
        else
            h2osfc[c] += xs1_arr[c]
            qflx_rsub_sat[c] = 0.0
        end

        # ice check
        xs1_arr[c] = smooth_max(smooth_max(h2osoi_ice[c, 1 + joff], 0.0) -
            smooth_max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_liq[c, 1 + joff]), 0.0)
        h2osoi_ice[c, 1 + joff] = smooth_min(smooth_max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_liq[c, 1 + joff]),
            h2osoi_ice[c, 1 + joff])
        qflx_ice_runoff_xs[c] = xs1_arr[c] / dtime
    end

    # ---- Limit h2osoi_liq >= watmin ----
    for j in 1:(nlevsoi-1)
        for c in bounds
            mask_hydrology[c] || continue
            if h2osoi_liq[c, j + joff] < WATMIN
                xs_arr[c] = WATMIN - h2osoi_liq[c, j + joff]
                if j == jwt[c]
                    zwt[c] += xs_arr[c] / eff_porosity[c, j] / 1000.0
                end
            else
                xs_arr[c] = 0.0
            end
            h2osoi_liq[c, j + joff]   += xs_arr[c]
            h2osoi_liq[c, j + 1 + joff] -= xs_arr[c]
        end
    end

    # Bottom layer: get water from above
    j = nlevsoi
    for c in bounds
        mask_hydrology[c] || continue
        if h2osoi_liq[c, j + joff] < WATMIN
            xs_arr[c] = WATMIN - h2osoi_liq[c, j + joff]
            for i in (nlevsoi-1):-1:1
                available_h2osoi_liq = smooth_max(h2osoi_liq[c, i + joff] - WATMIN - xs_arr[c], 0.0)
                if available_h2osoi_liq >= xs_arr[c]
                    h2osoi_liq[c, j + joff] += xs_arr[c]
                    h2osoi_liq[c, i + joff] -= xs_arr[c]
                    xs_arr[c] = 0.0
                    break
                else
                    h2osoi_liq[c, j + joff] += available_h2osoi_liq
                    h2osoi_liq[c, i + joff] -= available_h2osoi_liq
                    xs_arr[c] -= available_h2osoi_liq
                end
            end
        else
            xs_arr[c] = 0.0
        end
        h2osoi_liq[c, j + joff] += xs_arr[c]
        qflx_rsub_sat[c] -= xs_arr[c] / dtime
    end

    for c in bounds
        mask_hydrology[c] || continue
        qflx_drain[c] = qflx_rsub_sat[c] + drainage_arr[c]
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
