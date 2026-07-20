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
# Baseflow scalar (namelist adjustable). CTSM keys this on
# lower_boundary_condition (namelist_defaults_ctsm.xml:195-197). For clm5_0 the
# chain soilwater_movement_method=1 (defaults:419) + use_bedrock=.true.
# (defaults:180) resolves lbc=2 (defaults:426), giving 0.001. `1.d-2` is the
# lbc=1 / clm4_5 value. SoilHydrologyMod.F90:71's `= 1.e-2_r8` is only the
# pre-namelist initialiser that build-namelist overwrites -- the same
# code-fallback-instead-of-namelist-default trap as #252/#259/#273.
# Corroborated by two CTSM-emitted lnd_in files (Bow, MerBleue), both
# `baseflow_scalar = 0.001d00`. NOTE: this Ref, the SoilHydrologyConfig field,
# the init_soil_hydrology_config kwarg and clm_run!'s kwarg must all agree --
# test_soil_hydrology_mod.jl calls bare init_soil_hydrology_config(), so a
# partial change makes the live global depend on test ordering.
const BASEFLOW_SCALAR = Ref(0.001)
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
    baseflow_scalar::Float64   = 0.001   # CTSM clm5_0 (lbc=2); see BASEFLOW_SCALAR
end

"""
    init_soil_hydrology_config(; kwargs...) -> SoilHydrologyConfig

Create a soil hydrology configuration with defaults.
"""
function init_soil_hydrology_config(;
    head_gradient_method::Int = HEAD_GRADIENT_DARCY,
    transmissivity_method::Int = TRANSMISSIVITY_LAYERSUM,
    baseflow_scalar::Real = 0.001   # CTSM clm5_0 (lbc=2); see BASEFLOW_SCALAR
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
        # Axis-scaled sharpness — see SOIL_HYDRAULIC_K (soil_water_movement.jl). These clamp a
        # volumetric water content [m3/m3] and a relative ice fraction [-] against CONSTANT bounds,
        # so the clamped branch carries no derivative anyway. At the generic k = 50 the width is
        # 0.0139 in those units: a fully-frozen layer gets eff_porosity = 0.0195 instead of the
        # 0.01 floor (+95%, and eff_porosity gates infiltration), and icefrac ≈ 0.965 instead of 1,
        # which makes the ice impedance 10^(-e_ice·icefrac) ~1.6x too high.
        kh = T(SOIL_HYDRAULIC_K)
        dz_ext = col_dz[c, j + joff] + excess_ice[c, j] / dn
        vol_ice_val = smooth_min(watsat[c, j], (h2osoi_ice[c, j + joff] + excess_ice[c, j]) / (dz_ext * dn); k = kh)
        eff_porosity[c, j] = smooth_max(T(0.01), watsat[c, j] - vol_ice_val; k = kh)
        icefrac[c, j] = smooth_min(one(T), vol_ice_val / watsat[c, j]; k = kh)
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
#
# `joff` (= nlevsno) converts Fortran's soil-layer index to CLM.jl's COMBINED
# snow+soil array index: Fortran `h2osoi_liq(c,1)` (top SOIL layer) is Julia
# `h2osoi_liq[c, joff+1]`. These two urban kernels used a bare `h2osoi_liq[c, 1]`
# — which is the DEEPEST SNOW slot, nlevsno(=12) rows above the soil — so the
# roof / impervious-road ponding store they read was always the (zero) snow slot
# instead of the real pond. See the water-balance write-up in update_urban_ponding!.
@kernel function _soilhyd_surf_runoff_urban_kernel!(qflx_surf, xs_urban, @Const(mask_urban),
                                                    @Const(col_itype), @Const(col_snl),
                                                    @Const(qflx_rain_plus_snomelt),
                                                    @Const(h2osoi_liq),
                                                    @Const(qflx_liqevap_from_top_layer),
                                                    @Const(qflx_floodc),
                                                    icol_roof::Int, icol_road_imperv::Int,
                                                    icol_sunwall::Int, icol_shadewall::Int,
                                                    joff::Int, pondmx_urban, dtime)
    c = @index(Global)
    @inbounds if mask_urban[c]
        T = eltype(qflx_surf)
        dt = T(dtime)
        pmx = T(pondmx_urban)
        if col_itype[c] == icol_roof || col_itype[c] == icol_road_imperv
            if col_snl[c] < 0
                qflx_surf[c] = max(zero(T), qflx_rain_plus_snomelt[c])
            else
                xs_urban[c] = max(zero(T),
                    h2osoi_liq[c, joff + 1] / dt + qflx_rain_plus_snomelt[c] -
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
                           icol_sunwall::Int, icol_shadewall::Int, pondmx_urban, dtime;
                           joff::Int = varpar.nlevsno) =
    _launch!(_soilhyd_surf_runoff_urban_kernel!, qflx_surf, xs_urban, mask_urban,
             col_itype, col_snl, qflx_rain_plus_snomelt, h2osoi_liq,
             qflx_liqevap_from_top_layer, qflx_floodc, icol_roof, icol_road_imperv,
             icol_sunwall, icol_shadewall, joff,
             convert(eltype(qflx_surf), pondmx_urban),
             convert(eltype(qflx_surf), dtime))

# ---- update_urban_ponding! : ponding state on urban surfaces ----
@kernel function _soilhyd_urban_ponding_kernel!(h2osoi_liq, @Const(mask_urban),
                                                @Const(col_itype), @Const(col_snl),
                                                @Const(xs_urban),
                                                @Const(qflx_rain_plus_snomelt),
                                                @Const(qflx_liqevap_from_top_layer),
                                                icol_roof::Int, icol_road_imperv::Int,
                                                joff::Int, pondmx_urban, dtime)
    c = @index(Global)
    @inbounds if mask_urban[c]
        T = eltype(h2osoi_liq)        # working precision (Float32 on Metal); no Float64 literals
        zr = zero(T)
        if col_itype[c] == icol_roof || col_itype[c] == icol_road_imperv
            if col_snl[c] >= 0
                if xs_urban[c] > zr
                    h2osoi_liq[c, joff + 1] = pondmx_urban
                else
                    h2osoi_liq[c, joff + 1] = max(zr, h2osoi_liq[c, joff + 1] +
                        (qflx_rain_plus_snomelt[c] - qflx_liqevap_from_top_layer[c]) * dtime)
                end
            end
        end
    end
end

# ndrange MUST be the COLUMN count, not `length(out)`. `_launch!` defaults to
# `ndrange = length(out)` and `out` here is the h2osoi_liq MATRIX (nc x nlevtot),
# so the default would run nc*nlevtot threads and index mask_urban/col_itype far
# past their ends — silently reading garbage under @inbounds, and throwing under
# CI's --check-bounds=yes. Latent only because this routine was never called.
soilhyd_urban_ponding!(h2osoi_liq, mask_urban, col_itype, col_snl, xs_urban,
                       qflx_rain_plus_snomelt, qflx_liqevap_from_top_layer,
                       icol_roof::Int, icol_road_imperv::Int, pondmx_urban, dtime;
                       joff::Int = varpar.nlevsno) =
    _launch!(_soilhyd_urban_ponding_kernel!, h2osoi_liq, mask_urban, col_itype,
             col_snl, xs_urban, qflx_rain_plus_snomelt, qflx_liqevap_from_top_layer,
             icol_roof, icol_road_imperv, joff,
             convert(eltype(h2osoi_liq), pondmx_urban),
             convert(eltype(h2osoi_liq), dtime);
             ndrange = length(mask_urban))

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
             convert(eltype(h2osoi_liq), dtime); ndrange = (length(mask), nlevsoi))

# ---- withdraw_groundwater_irrigation! : per-column confined aquifer ----
@kernel function _soilhyd_withdraw_gw_con_kernel!(wa, @Const(mask),
                                                  @Const(qflx_gw_con_irrig), dtime)
    c = @index(Global)
    @inbounds if mask[c]
        wa[c] = wa[c] - qflx_gw_con_irrig[c] * dtime
    end
end

soilhyd_withdraw_gw_con!(wa, mask, qflx_gw_con_irrig, dtime) =
    _launch!(_soilhyd_withdraw_gw_con_kernel!, wa, mask, qflx_gw_con_irrig,
             convert(eltype(wa), dtime))

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
    col_snl::AbstractVector{<:Integer},
    col_itype::AbstractVector{<:Integer},
    mask_urban::AbstractVector{Bool},
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
        nlevsoi::Int, joff::Int, joff_zi::Int, dtime, aq_sp_yield_min, tfrz, denh2o, denice,
        recompute_frost_table::Bool)
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
        rous = max(rous, aqy)

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
                    s_y = max(s_y, aqy)

                    qcharge_layer = min(qcharge_tot,
                        s_y * (odv.zwt[c] - idv.col_zi[c, joff_zi + j - 1]) * e3)
                    qcharge_layer = max(qcharge_layer, zr)

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
                    s_y = max(s_y, aqy)

                    qcharge_layer = max(qcharge_tot,
                        -(s_y * (idv.col_zi[c, joff_zi + j] - odv.zwt[c]) * e3))
                    qcharge_layer = min(qcharge_layer, zr)
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
        # Frost table + perched water table. In Fortran this lives in WaterTable
        # (HydrologyNoDrainageMod.F90:357), which is ONLY called when
        # use_aquifer_layer==true; the non-aquifer path uses ThetaBasedWaterTable,
        # which leaves frost_table/zwt_perched at their restart values (the perched
        # drainage is handled later by perched_lateral_flow!). Gate accordingly so we
        # do not spuriously recompute frost_table = col_z[k_frz] (node) every step.
        if recompute_frost_table
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
                    odv.zwt_perched[c] = max(zr, m_val * sat_lev + b_val)
                end
            end
        end
    end
end

function soilhyd_water_table!(odv, idv, mask, nlevsoi::Int, joff::Int, joff_zi::Int,
                              dtime, aq_sp_yield_min, tfrz, denh2o, denice;
                              recompute_frost_table::Bool=true)
    # Convert all scalar reals to the working precision of the device arrays so no
    # Float64 reaches the Float32-only Metal backend (byte-identical on Float64 CPU).
    T = eltype(odv.zwt)
    _launch!(_water_table_kernel!, odv.zwt, odv, idv, mask, nlevsoi, joff, joff_zi,
             T(dtime), T(aq_sp_yield_min), T(tfrz), T(denh2o), T(denice),
             recompute_frost_table;
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
    dtime::Real;
    recompute_frost_table::Bool=true
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
                         dtime, params.aq_sp_yield_min, TFRZ, DENH2O, DENICE;
                         recompute_frost_table=recompute_frost_table)

    return nothing
end

# =========================================================================
# Drainage
# =========================================================================

# --------------------------------------------------------------------------
# drainage! : ONE thread per column. The whole Drainage subroutine is per-column
# with internal SEQUENTIAL nlevsoi level sweeps (jwt search, perched/baseflow
# saturated-layer cascades with loop-carried zwt, the upward saturation-excess
# redistribution and the watmin pull-from-below sweeps). All sweeps live inside
# a column, so each thread runs the full subroutine in-thread and writes only
# its own column's outputs (no cross-column dependence). The per-column scratch
# vectors of the scalar version (jwt, rsub_top, fff, xs*_arr) become in-thread
# locals; dzmm is recomputed inline as col_dz*1e3.
#
# Read-only inputs and read/write outputs are grouped into two immutable
# device-view bundles (Adapt.@adapt_structure'd); field names mirror the
# original Julia locals so the kernel body reads verbatim. On the host path the
# same bundles are built from CPU arrays, so CPU stays byte-identical.
# --------------------------------------------------------------------------

# Soil-state / geometry / forcing inputs (read-only) the sweeps touch.
Base.@kwdef struct _DrInDV{VR,VI,VB,M}
    col_dz::M; col_z::M; col_zi::M; t_soisno::M
    watsat::M; sucsat::M; bsw::M; hksat::M; eff_porosity::M
    hkdepth::VR; col_topo_slope::VR; qflx_snwcp_liq::VR
    col_itype::VI; col_landunit::VI; lun_urbpoi::VB
end
Adapt.@adapt_structure _DrInDV

# State outputs (written and/or read-modified) the sweeps touch.
Base.@kwdef struct _DrOutDV{VR,M}
    zwt::VR; zwt_perched::VR; frost_table::VR; wa::VR; h2osfc::VR
    h2osoi_liq::M; h2osoi_ice::M; icefrac::M
    qflx_drain::VR; qflx_rsub_sat::VR; qflx_drain_perched::VR
    qflx_qrgwl::VR; qflx_ice_runoff_xs::VR
end
Adapt.@adapt_structure _DrOutDV

@kernel function _drainage_kernel!(zwt_out, odv, @Const(idv), @Const(mask_hydrology),
        @Const(mask_urban), nlevsoi::Int, joff::Int, joff_zi::Int, dtime,
        perched_baseflow_scalar, e_ice, aq_sp_yield_min, aquifer_water_baseline,
        h2osfcflag::Int, icol_road_perv::Int, tfrz, denh2o, denice, watmin, pondmx, rpi)
    c = @index(Global)
    @inbounds begin
        T     = eltype(zwt_out)
        zr    = zero(T)
        on    = one(T)
        e3    = T(1.0e3)
        thou  = T(1000.0)
        TFRZl = T(tfrz)
        DH2O  = T(denh2o)
        DICE  = T(denice)
        WMIN  = T(watmin)
        PMX   = T(pondmx)
        EICE  = T(e_ice)
        AQY   = T(aq_sp_yield_min)
        PBFS  = T(perched_baseflow_scalar)
        AQWB  = T(aquifer_water_baseline)
        deg2rad = T(rpi) / T(180.0)

        rsub_top = zr   # per-column carried scalar (init 0)

        if mask_hydrology[c]
            # ---- icefrac (per layer) ----
            for j in 1:nlevsoi
                vol_ice_val = min(idv.watsat[c, j],
                    odv.h2osoi_ice[c, j + joff] / (idv.col_dz[c, j + joff] * DICE))
                odv.icefrac[c, j] = min(on, vol_ice_val / idv.watsat[c, j])
            end

            # ---- init drainage fluxes ----
            odv.qflx_drain[c]    = zr
            odv.qflx_rsub_sat[c] = zr

            # ---- jwt: layer index right above water table ----
            jwt = nlevsoi
            for j in 1:nlevsoi
                if odv.zwt[c] <= idv.col_zi[c, j + joff_zi]
                    jwt = j - 1
                    break
                end
            end

            # ============================ BASEFLOW ============================
            q_perch_max = PBFS * sin(idv.col_topo_slope[c] * deg2rad)

            # define frost table
            if idv.t_soisno[c, 1 + joff] > TFRZl
                k_frz = nlevsoi
            else
                k_frz = 1
            end
            for k in 2:nlevsoi
                if idv.t_soisno[c, k-1 + joff] > TFRZl && idv.t_soisno[c, k + joff] <= TFRZl
                    k_frz = k
                    break
                end
            end

            odv.frost_table[c] = idv.col_z[c, k_frz + joff]
            odv.zwt_perched[c] = odv.frost_table[c]
            odv.qflx_drain_perched[c] = zr

            if odv.zwt[c] < odv.frost_table[c] && idv.t_soisno[c, k_frz + joff] <= TFRZl
                # ======= water table above frost table =======
                wtsub = zr
                q_perch = zr
                for k in (jwt+1):k_frz
                    imped = T(10.0)^(-EICE * (T(0.5) * (odv.icefrac[c, k] + odv.icefrac[c, min(nlevsoi, k+1)])))
                    q_perch = q_perch + imped * idv.hksat[c, k] * (idv.col_dz[c, k + joff] * e3)
                    wtsub = wtsub + (idv.col_dz[c, k + joff] * e3)
                end
                if wtsub > zr
                    q_perch = q_perch / wtsub
                end

                odv.qflx_drain_perched[c] = q_perch_max * q_perch * (odv.frost_table[c] - odv.zwt[c])

                rsub_top_tot = -odv.qflx_drain_perched[c] * dtime
                for k in (jwt+1):k_frz
                    rsub_top_layer = max(rsub_top_tot, -(odv.h2osoi_liq[c, k + joff] - WMIN))
                    rsub_top_layer = min(rsub_top_layer, zr)
                    rsub_top_tot = rsub_top_tot - rsub_top_layer

                    odv.h2osoi_liq[c, k + joff] = odv.h2osoi_liq[c, k + joff] + rsub_top_layer

                    if rsub_top_tot >= zr
                        odv.zwt[c] = odv.zwt[c] - rsub_top_layer / idv.eff_porosity[c, k] / thou
                        break
                    else
                        odv.zwt[c] = idv.col_zi[c, k + joff_zi]
                    end
                end

                odv.qflx_drain_perched[c] = odv.qflx_drain_perched[c] + rsub_top_tot / dtime

                # recompute jwt
                jwt = nlevsoi
                for j in 1:nlevsoi
                    if odv.zwt[c] <= idv.col_zi[c, j + joff_zi]
                        jwt = j - 1
                        break
                    end
                end
            else
                # ======= water table below frost table =======
                sat_lev = T(0.9)

                k_perch = 1
                for k in k_frz:-1:1
                    h2osoi_vol_val = odv.h2osoi_liq[c, k + joff] / (idv.col_dz[c, k + joff] * DH2O) +
                        odv.h2osoi_ice[c, k + joff] / (idv.col_dz[c, k + joff] * DICE)
                    if h2osoi_vol_val / idv.watsat[c, k] <= sat_lev
                        k_perch = k
                        break
                    end
                end

                if idv.t_soisno[c, k_frz + joff] > TFRZl
                    k_perch = k_frz
                end

                if k_frz > k_perch
                    s1 = (odv.h2osoi_liq[c, k_perch + joff] / (idv.col_dz[c, k_perch + joff] * DH2O) +
                        odv.h2osoi_ice[c, k_perch + joff] / (idv.col_dz[c, k_perch + joff] * DICE)) / idv.watsat[c, k_perch]
                    s2 = (odv.h2osoi_liq[c, k_perch+1 + joff] / (idv.col_dz[c, k_perch+1 + joff] * DH2O) +
                        odv.h2osoi_ice[c, k_perch+1 + joff] / (idv.col_dz[c, k_perch+1 + joff] * DICE)) / idv.watsat[c, k_perch+1]

                    m_val = (idv.col_z[c, k_perch+1 + joff] - idv.col_z[c, k_perch + joff]) / (s2 - s1)
                    b_val = idv.col_z[c, k_perch+1 + joff] - m_val * s2
                    odv.zwt_perched[c] = max(zr, m_val * sat_lev + b_val)

                    wtsub = zr
                    q_perch = zr
                    for k in k_perch:k_frz
                        imped = T(10.0)^(-EICE * (T(0.5) * (odv.icefrac[c, k] + odv.icefrac[c, min(nlevsoi, k+1)])))
                        q_perch = q_perch + imped * idv.hksat[c, k] * (idv.col_dz[c, k + joff] * e3)
                        wtsub = wtsub + (idv.col_dz[c, k + joff] * e3)
                    end
                    if wtsub > zr
                        q_perch = q_perch / wtsub
                    end

                    odv.qflx_drain_perched[c] = q_perch_max * q_perch * (odv.frost_table[c] - odv.zwt_perched[c])

                    rsub_top_tot = -odv.qflx_drain_perched[c] * dtime
                    for k in (k_perch+1):k_frz
                        rsub_top_layer = max(rsub_top_tot, -(odv.h2osoi_liq[c, k + joff] - WMIN))
                        rsub_top_layer = min(rsub_top_layer, zr)
                        rsub_top_tot = rsub_top_tot - rsub_top_layer

                        odv.h2osoi_liq[c, k + joff] = odv.h2osoi_liq[c, k + joff] + rsub_top_layer

                        if rsub_top_tot >= zr
                            odv.zwt_perched[c] = odv.zwt_perched[c] - rsub_top_layer / idv.eff_porosity[c, k] / thou
                            break
                        else
                            odv.zwt_perched[c] = idv.col_zi[c, k + joff_zi]
                        end
                    end

                    odv.qflx_drain_perched[c] = odv.qflx_drain_perched[c] + rsub_top_tot / dtime
                else
                    odv.qflx_drain_perched[c] = zr
                end

                # -- Topographic runoff --
                fff = on / idv.hkdepth[c]
                dzsum = zr
                icefracsum = zr
                for j in max(jwt, 1):nlevsoi
                    dzj = idv.col_dz[c, j + joff] * e3
                    dzsum = dzsum + dzj
                    icefracsum = icefracsum + odv.icefrac[c, j] * dzj
                end

                imped = T(10.0)^(-EICE * (icefracsum / dzsum))
                rsub_top_max = T(10.0) * sin(deg2rad * idv.col_topo_slope[c])

                rsub_top = imped * rsub_top_max * exp(-fff * odv.zwt[c])

                # analytical expression for aquifer specific yield
                _rous_base = on + e3 * odv.zwt[c] / idv.sucsat[c, nlevsoi]
                _rous_exp = -on / idv.bsw[c, nlevsoi]
                rous_local = idv.watsat[c, nlevsoi] * (on - _rous_base^_rous_exp)
                rous_local = max(rous_local, AQY)

                if jwt == nlevsoi
                    # water table below soil column
                    odv.wa[c] = odv.wa[c] - rsub_top * dtime
                    odv.zwt[c] = odv.zwt[c] + (rsub_top * dtime) / thou / rous_local
                    odv.h2osoi_liq[c, nlevsoi + joff] = odv.h2osoi_liq[c, nlevsoi + joff] +
                        max(zr, odv.wa[c] - AQWB)
                    odv.wa[c] = min(odv.wa[c], AQWB)
                else
                    # water table within soil layers
                    rsub_top_tot = -rsub_top * dtime

                    # rsub_top_tot <= 0 for valid inputs; positive branch is a no-op
                    # here (the scalar CPU version errors / warns — never reached for
                    # physically valid forcing, so this stays byte-identical).
                    if rsub_top_tot <= zr
                        for j in (jwt+1):nlevsoi
                            s_y = idv.watsat[c, j] *
                                (on - (on + e3 * odv.zwt[c] / idv.sucsat[c, j])^(-on / idv.bsw[c, j]))
                            s_y = max(s_y, AQY)

                            rsub_top_layer = max(rsub_top_tot,
                                -(s_y * (idv.col_zi[c, j + joff_zi] - odv.zwt[c]) * e3))
                            rsub_top_layer = min(rsub_top_layer, zr)
                            odv.h2osoi_liq[c, j + joff] = odv.h2osoi_liq[c, j + joff] + rsub_top_layer

                            rsub_top_tot = rsub_top_tot - rsub_top_layer

                            if rsub_top_tot >= zr
                                odv.zwt[c] = odv.zwt[c] - rsub_top_layer / s_y / thou
                                break
                            else
                                odv.zwt[c] = idv.col_zi[c, j + joff_zi]
                            end
                        end

                        odv.zwt[c] = odv.zwt[c] - rsub_top_tot / thou / rous_local
                        odv.wa[c] = odv.wa[c] + rsub_top_tot
                    end

                    # recompute jwt
                    jwt = nlevsoi
                    for j in 1:nlevsoi
                        if odv.zwt[c] <= idv.col_zi[c, j + joff_zi]
                            jwt = j - 1
                            break
                        end
                    end
                end

                odv.zwt[c] = max(zr, odv.zwt[c])
                odv.zwt[c] = min(T(80.0), odv.zwt[c])
            end

            # ---- Excessive water above saturation: redistribute upward ----
            for j in nlevsoi:-1:2
                effp_dz = idv.eff_porosity[c, j] * (idv.col_dz[c, j + joff] * e3)
                xsi = max(odv.h2osoi_liq[c, j + joff] - effp_dz, zr)
                odv.h2osoi_liq[c, j + joff]   = min(effp_dz, odv.h2osoi_liq[c, j + joff])
                odv.h2osoi_liq[c, j-1 + joff] = odv.h2osoi_liq[c, j-1 + joff] + xsi
            end

            # ---- watmin addition to fix water balance errors ----
            dzmm1 = idv.col_dz[c, 1 + joff] * e3
            xs1 = max(max(odv.h2osoi_liq[c, 1 + joff] - WMIN, zr) -
                max(zr, PMX + idv.watsat[c, 1] * dzmm1 - odv.h2osoi_ice[c, 1 + joff] - WMIN), zr)
            odv.h2osoi_liq[c, 1 + joff] = odv.h2osoi_liq[c, 1 + joff] - xs1

            l = idv.col_landunit[c]
            if idv.lun_urbpoi[l]
                odv.qflx_rsub_sat[c] = xs1 / dtime
            else
                if h2osfcflag == 1
                    odv.h2osfc[c] = odv.h2osfc[c] + xs1
                    odv.qflx_rsub_sat[c] = zr
                else
                    odv.qflx_rsub_sat[c] = xs1 / dtime
                end
            end

            # ---- ice check ----
            xs1 = max(max(odv.h2osoi_ice[c, 1 + joff], zr) -
                max(zr, PMX + idv.watsat[c, 1] * dzmm1 - odv.h2osoi_liq[c, 1 + joff]), zr)
            odv.h2osoi_ice[c, 1 + joff] = min(
                max(zr, PMX + idv.watsat[c, 1] * dzmm1 - odv.h2osoi_liq[c, 1 + joff]),
                odv.h2osoi_ice[c, 1 + joff])
            odv.qflx_ice_runoff_xs[c] = xs1 / dtime

            # ---- Limit h2osoi_liq >= watmin: pull from below ----
            for j in 1:(nlevsoi-1)
                if odv.h2osoi_liq[c, j + joff] < WMIN
                    xs = WMIN - odv.h2osoi_liq[c, j + joff]
                    if j == jwt
                        odv.zwt[c] = odv.zwt[c] + xs / idv.eff_porosity[c, j] / thou
                    end
                else
                    xs = zr
                end
                odv.h2osoi_liq[c, j + joff]   = odv.h2osoi_liq[c, j + joff]   + xs
                odv.h2osoi_liq[c, j+1 + joff] = odv.h2osoi_liq[c, j+1 + joff] - xs
            end

            # ---- Get water for bottom layer from layers above ----
            jb = nlevsoi
            if odv.h2osoi_liq[c, jb + joff] < WMIN
                xs = WMIN - odv.h2osoi_liq[c, jb + joff]
                for i in (nlevsoi-1):-1:1
                    available = max(odv.h2osoi_liq[c, i + joff] - WMIN - xs, zr)
                    if available >= xs
                        odv.h2osoi_liq[c, jb + joff] = odv.h2osoi_liq[c, jb + joff] + xs
                        odv.h2osoi_liq[c, i + joff]  = odv.h2osoi_liq[c, i + joff] - xs
                        xs = zr
                        break
                    else
                        odv.h2osoi_liq[c, jb + joff] = odv.h2osoi_liq[c, jb + joff] + available
                        odv.h2osoi_liq[c, i + joff]  = odv.h2osoi_liq[c, i + joff] - available
                        xs = xs - available
                    end
                end
            else
                xs = zr
            end
            odv.h2osoi_liq[c, jb + joff] = odv.h2osoi_liq[c, jb + joff] + xs
            rsub_top = rsub_top - xs / dtime

            odv.qflx_drain[c] = odv.qflx_rsub_sat[c] + rsub_top
            odv.qflx_qrgwl[c] = idv.qflx_snwcp_liq[c]
        end

        # ---- No drainage for urban columns (except pervious road) ----
        if mask_urban[c]
            if idv.col_itype[c] != icol_road_perv
                odv.qflx_drain[c] = zr
                odv.qflx_qrgwl[c] = idv.qflx_snwcp_liq[c]
            end
        end
    end
end

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
    temperature_t_soisno::AbstractMatrix{<:Real},
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    waterstatebulk_ws::WaterStateData,
    waterfluxbulk::WaterFluxBulkData,
    col_dz::AbstractMatrix{<:Real},
    col_z::AbstractMatrix{<:Real},
    col_zi::AbstractMatrix{<:Real},
    col_snl::AbstractVector{Int},
    col_itype::AbstractVector{Int},
    col_landunit::AbstractVector{Int},
    col_topo_slope::AbstractVector{<:Real},
    lun_urbpoi::AbstractVector{Bool},
    mask_hydrology::AbstractVector{Bool},
    mask_urban::AbstractVector{Bool},
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

    idv = _DrInDV(;
        col_dz         = col_dz,
        col_z          = col_z,
        col_zi         = col_zi,
        t_soisno       = temperature_t_soisno,
        watsat         = soilstate.watsat_col,
        sucsat         = soilstate.sucsat_col,
        bsw            = soilstate.bsw_col,
        hksat          = soilstate.hksat_col,
        eff_porosity   = soilstate.eff_porosity_col,
        hkdepth        = soilhydrology.hkdepth_col,
        col_topo_slope = col_topo_slope,
        qflx_snwcp_liq = wf.qflx_snwcp_liq_col,
        col_itype      = col_itype,
        col_landunit   = col_landunit,
        lun_urbpoi     = lun_urbpoi)

    odv = _DrOutDV(;
        zwt                = soilhydrology.zwt_col,
        zwt_perched        = soilhydrology.zwt_perched_col,
        frost_table        = soilhydrology.frost_table_col,
        wa                 = waterstatebulk_ws.wa_col,
        h2osfc             = waterstatebulk_ws.h2osfc_col,
        h2osoi_liq         = waterstatebulk_ws.h2osoi_liq_col,
        h2osoi_ice         = waterstatebulk_ws.h2osoi_ice_col,
        icefrac            = soilhydrology.icefrac_col,
        qflx_drain         = wf.qflx_drain_col,
        qflx_rsub_sat      = wf.qflx_rsub_sat_col,
        qflx_drain_perched = wf.qflx_drain_perched_col,
        qflx_qrgwl         = wf.qflx_qrgwl_col,
        qflx_ice_runoff_xs = wf.qflx_ice_runoff_xs_col)

    # Convert every scalar real to the working precision of the device arrays so
    # no Float64 literal reaches the Float32-only Metal backend (byte-identical on
    # a Float64 CPU run: T(x) === x there).
    T = eltype(odv.zwt)
    _launch!(_drainage_kernel!, odv.zwt, odv, idv, mask_hydrology, mask_urban,
             nlevsoi, joff, joff_zi, T(dtime),
             T(params.perched_baseflow_scalar), T(params.e_ice), T(params.aq_sp_yield_min),
             T(waterstatebulk_ws.aquifer_water_baseline), soilhydrology.h2osfcflag,
             ICOL_ROAD_PERV, T(TFRZ), T(DENH2O), T(DENICE), T(WATMIN), T(PONDMX), T(RPI);
             ndrange = length(mask_hydrology))

    return nothing
end

# =========================================================================
# CLMVICMap
# =========================================================================

# --------------------------------------------------------------------------
# clm_vic_map! : ONE thread per (hydrology) column. The CLM->VIC layer map is
# fully per-column: each thread sums its own CLM layers into VIC layers, then
# accumulates the per-column top-layer totals. No cross-column dependence, so
# the whole subroutine runs in-thread; writes touch only the thread's column.
#
# All read-only inputs + read/write outputs are grouped into two immutable
# device-view bundles (Adapt.@adapt_structure'd); field names mirror the
# original Julia locals so the kernel body reads verbatim. On the host path the
# same bundles are built from CPU arrays, so CPU stays byte-identical.
# --------------------------------------------------------------------------

# Read-only inputs the map touches.
Base.@kwdef struct _VicInDV{M,A3}
    h2osoi_liq::M; h2osoi_ice::M; h2osoi_vol::M
    depth::M; porosity::M; max_moist::M
    vic_clm_fract::A3
end
Adapt.@adapt_structure _VicInDV

# State outputs (written / read-modified) the map touches.
Base.@kwdef struct _VicOutDV{V,M}
    moist::M; ice::M; moist_vol::M
    top_moist::V; top_max_moist::V; top_ice::V; top_moist_limited::V
end
Adapt.@adapt_structure _VicOutDV

@kernel function _clm_vic_map_kernel!(moist_out, odv, @Const(idv), @Const(mask),
        nlevsoi::Int, nlevgrnd::Int, nlayer::Int, nlayert::Int, watmin, denh2o, denice)
    c = @index(Global)
    @inbounds if mask[c]
        T   = eltype(moist_out)
        zr  = zero(T)
        WMIN = T(watmin)
        DH2O = T(denh2o)
        DICE = T(denice)
        small = T(0.01)

        for i in 1:nlayer
            ice0   = odv.ice[c, i]
            moist0 = odv.moist[c, i]
            odv.ice[c, i]   = zr
            odv.moist[c, i] = zr
            for j in 1:nlevsoi
                odv.ice[c, i]   = odv.ice[c, i]   + idv.h2osoi_ice[c, j] * idv.vic_clm_fract[c, i, j]
                odv.moist[c, i] = odv.moist[c, i] + idv.h2osoi_liq[c, j] * idv.vic_clm_fract[c, i, j]
            end
            odv.ice[c, i]       = min(moist0 + ice0, odv.ice[c, i])
            odv.ice[c, i]       = max(zr, odv.ice[c, i])
            odv.moist[c, i]     = max(WMIN, odv.moist[c, i])
            odv.moist[c, i]     = min(idv.max_moist[c, i] - odv.ice[c, i], odv.moist[c, i])
            odv.moist_vol[c, i] = odv.moist[c, i] / (idv.depth[c, i] * DICE) +
                                  odv.ice[c, i]   / (idv.depth[c, i] * DH2O)
            odv.moist_vol[c, i] = min(idv.porosity[c, i], odv.moist_vol[c, i])
            odv.moist_vol[c, i] = max(small, odv.moist_vol[c, i])
        end

        # hydrologic inactive layers
        for k in (nlayer+1):nlayert
            j_clm = nlevsoi + (k - nlayer)
            if j_clm <= nlevgrnd
                odv.ice[c, k]       = idv.h2osoi_ice[c, j_clm]
                odv.moist[c, k]     = idv.h2osoi_liq[c, j_clm]
                odv.moist_vol[c, k] = idv.h2osoi_vol[c, j_clm]
            end
        end

        # top VIC layer totals
        odv.top_moist[c]     = zr
        odv.top_ice[c]       = zr
        odv.top_max_moist[c] = zr
        for j in 1:(nlayer - 1)
            odv.top_ice[c]       = odv.top_ice[c]   + odv.ice[c, j]
            odv.top_moist[c]     = odv.top_moist[c] + odv.moist[c, j] + odv.ice[c, j]
            odv.top_max_moist[c] = odv.top_max_moist[c] + idv.max_moist[c, j]
        end
        odv.top_moist_limited[c] = min(odv.top_moist[c], odv.top_max_moist[c])
    end
end

"""
    clm_vic_map!(soilhydrology, waterstatebulk_ws, col_dz, col_zi, col_z,
                  mask, bounds, nlevsoi, nlevgrnd, nlayer, nlayert)

Map CLM layers to VIC layers.

Ported from `CLMVICMap` in `SoilHydrologyMod.F90`.
"""
function clm_vic_map!(
    soilhydrology::SoilHydrologyData,
    waterstatebulk_ws::WaterStateData,
    col_dz::AbstractMatrix{<:Real},
    col_zi::AbstractMatrix{<:Real},
    col_z::AbstractMatrix{<:Real},
    mask::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsoi::Int,
    nlevgrnd::Int,
    nlayer::Int,
    nlayert::Int
)
    idv = _VicInDV(;
        h2osoi_liq    = waterstatebulk_ws.h2osoi_liq_col,
        h2osoi_ice    = waterstatebulk_ws.h2osoi_ice_col,
        h2osoi_vol    = waterstatebulk_ws.h2osoi_vol_col,
        depth         = soilhydrology.depth_col,
        porosity      = soilhydrology.porosity_col,
        max_moist     = soilhydrology.max_moist_col,
        vic_clm_fract = soilhydrology.vic_clm_fract_col)

    odv = _VicOutDV(;
        moist             = soilhydrology.moist_col,
        ice               = soilhydrology.ice_col,
        moist_vol         = soilhydrology.moist_vol_col,
        top_moist         = soilhydrology.top_moist_col,
        top_max_moist     = soilhydrology.top_max_moist_col,
        top_ice           = soilhydrology.top_ice_col,
        top_moist_limited = soilhydrology.top_moist_limited_col)

    # Convert scalar constants to the device working precision (no Float64 reaches
    # the Float32-only Metal backend; identity on a Float64 CPU run).
    T = eltype(odv.moist)
    _launch!(_clm_vic_map_kernel!, odv.moist, odv, idv, mask,
             nlevsoi, nlevgrnd, nlayer, nlayert, T(WATMIN), T(DENH2O), T(DENICE);
             ndrange = length(mask))

    return nothing
end

# =========================================================================
# PerchedWaterTable
# =========================================================================

# --------------------------------------------------------------------------
# perched_water_table! : ONE thread per (hydrology) column. The whole
# PerchedWaterTable subroutine is per-column with an internal SEQUENTIAL frost-
# table search (2→nlevsoi) carrying loop-carried k_frz state, then a bottom-up
# k_frz→1 saturation search carrying k_perch, then an interpolation branch.
# Mirrors water_table!'s device-view layout: read-only inputs and written
# outputs grouped into two immutable Adapt.@adapt_structure'd bundles whose
# field names mirror the host locals so the kernel body reads verbatim. On the
# host path the same bundles wrap CPU arrays, so CPU stays byte-identical.
#
# Indexing is direct [c,k] (no snow offset) — matching the host loop.
# --------------------------------------------------------------------------

# Read-only inputs (all matrices + the read-only zwt vector).
Base.@kwdef struct _PWTInDV{V,M}
    col_dz::M; col_z::M; col_zi::M
    t_soisno::M; watsat::M
    h2osoi_liq::M; h2osoi_ice::M
    zwt::V
end
Adapt.@adapt_structure _PWTInDV

# Written outputs (zwt_perched + frost_table vectors, h2osoi_vol matrix).
Base.@kwdef struct _PWTOutDV{V,M}
    zwt_perched::V; frost_table::V; h2osoi_vol::M
end
Adapt.@adapt_structure _PWTOutDV

@kernel function _perched_water_table_kernel!(zwt_perched_out, odv, @Const(idv), @Const(mask),
        nlevsoi::Int, joff::Int, joff_zi::Int, tfrz, denh2o, denice)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(zwt_perched_out)
        zr    = zero(T)
        TFRZl = T(tfrz)
        DH2O  = T(denh2o)
        DICE  = T(denice)
        sat_lev = T(0.9)

        # define frost table as first frozen layer with unfrozen layer above it
        if idv.t_soisno[c, 1 + joff] > TFRZl
            k_frz = nlevsoi
        else
            k_frz = 1
        end

        for k in 2:nlevsoi
            if idv.t_soisno[c, k-1 + joff] > TFRZl && idv.t_soisno[c, k + joff] <= TFRZl
                k_frz = k
                break
            end
        end

        # frost table is top of frozen layer
        odv.frost_table[c] = idv.col_zi[c, k_frz - 1 + joff_zi]
        odv.zwt_perched[c] = odv.frost_table[c]

        # water table above frost table: do nothing
        if idv.zwt[c] < odv.frost_table[c] && idv.t_soisno[c, k_frz + joff] <= TFRZl
            # do nothing
        elseif k_frz > 1
            # water table below frost table
            k_perch = 1
            for k in k_frz:-1:1
                odv.h2osoi_vol[c, k] = idv.h2osoi_liq[c, k + joff] / (idv.col_dz[c, k + joff] * DH2O) +
                    idv.h2osoi_ice[c, k + joff] / (idv.col_dz[c, k + joff] * DICE)

                if odv.h2osoi_vol[c, k] / idv.watsat[c, k] <= sat_lev
                    k_perch = k
                    break
                end
            end

            if idv.t_soisno[c, k_frz + joff] > TFRZl
                k_perch = k_frz
            end

            if k_frz > k_perch
                s1 = (idv.h2osoi_liq[c, k_perch + joff] / (idv.col_dz[c, k_perch + joff] * DH2O) +
                    idv.h2osoi_ice[c, k_perch + joff] / (idv.col_dz[c, k_perch + joff] * DICE)) / idv.watsat[c, k_perch]
                s2 = (idv.h2osoi_liq[c, k_perch+1 + joff] / (idv.col_dz[c, k_perch+1 + joff] * DH2O) +
                    idv.h2osoi_ice[c, k_perch+1 + joff] / (idv.col_dz[c, k_perch+1 + joff] * DICE)) / idv.watsat[c, k_perch+1]

                if s1 > s2
                    odv.zwt_perched[c] = idv.col_zi[c, k_perch - 1 + joff_zi]
                else
                    m_val = (idv.col_z[c, k_perch+1 + joff] - idv.col_z[c, k_perch + joff]) / (s2 - s1)
                    b_val = idv.col_z[c, k_perch+1 + joff] - m_val * s2
                    odv.zwt_perched[c] = max(zr, m_val * sat_lev + b_val)
                end
            end
        end
    end
end

function soilhyd_perched_water_table!(odv, idv, mask, nlevsoi::Int, joff::Int, joff_zi::Int, tfrz, denh2o, denice)
    T = eltype(odv.zwt_perched)
    _launch!(_perched_water_table_kernel!, odv.zwt_perched, odv, idv, mask, nlevsoi,
             joff, joff_zi, T(tfrz), T(denh2o), T(denice); ndrange = length(mask))
end

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
    temperature_t_soisno::AbstractMatrix{<:Real},
    waterstatebulk_ws::WaterStateData,
    col_dz::AbstractMatrix{<:Real},
    col_z::AbstractMatrix{<:Real},
    col_zi::AbstractMatrix{<:Real},
    mask_hydrology::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsoi::Int;
    joff::Int=varpar.nlevsno,
    joff_zi::Int=varpar.nlevsno + 1
)
    idv = _PWTInDV(;
        col_dz     = col_dz,
        col_z      = col_z,
        col_zi     = col_zi,
        t_soisno   = temperature_t_soisno,
        watsat     = soilstate.watsat_col,
        h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col,
        h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col,
        zwt        = soilhydrology.zwt_col)

    odv = _PWTOutDV(;
        zwt_perched = soilhydrology.zwt_perched_col,
        frost_table = soilhydrology.frost_table_col,
        h2osoi_vol  = waterstatebulk_ws.h2osoi_vol_col)

    soilhyd_perched_water_table!(odv, idv, mask_hydrology, nlevsoi, joff, joff_zi, TFRZ, DENH2O, DENICE)

    return nothing
end

# =========================================================================
# ThetaBasedWaterTable
# =========================================================================

# --------------------------------------------------------------------------
# theta_based_water_table! : ONE thread per (hydrology) column. The whole
# ThetaBasedWaterTable subroutine is per-column with an internal SEQUENTIAL
# nbedrock→1 saturation search carrying loop-carried k_zwt / sat_flag state,
# then a branch on k_zwt for the water-table depth. Mirrors water_table!'s
# device-view layout: read-only inputs and written outputs grouped into two
# immutable Adapt.@adapt_structure'd bundles whose field names mirror the
# host locals so the kernel body reads verbatim. On the host path the same
# bundles wrap CPU arrays, so CPU stays byte-identical.
#
# Indexing is direct [c,k] (no snow offset) — matching the host loop and the
# existing test geometry.
# --------------------------------------------------------------------------

# Read-only inputs the search loop touches (all matrices + the Int nbedrock vec).
Base.@kwdef struct _TBWTInDV{M,VI}
    col_dz::M; col_z::M; col_zi::M
    watsat::M; h2osoi_liq::M; h2osoi_ice::M
    col_nbedrock::VI
end
Adapt.@adapt_structure _TBWTInDV

# Written outputs (zwt vector + h2osoi_vol matrix).
Base.@kwdef struct _TBWTOutDV{V,M}
    zwt::V; h2osoi_vol::M
end
Adapt.@adapt_structure _TBWTOutDV

@kernel function _theta_based_water_table_kernel!(zwt_out, odv, @Const(idv), @Const(mask),
        nlevsoi::Int, joff::Int, joff_zi::Int, denh2o, denice)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(zwt_out)
        zr   = zero(T)
        DH2O = T(denh2o)
        DICE = T(denice)
        sat_lev = T(0.9)

        nbed = idv.col_nbedrock[c]

        odv.zwt[c] = idv.col_zi[c, nlevsoi + joff_zi]

        k_zwt = nbed
        sat_flag = 1
        for k in nbed:-1:1
            odv.h2osoi_vol[c, k] = idv.h2osoi_liq[c, k + joff] / (idv.col_dz[c, k + joff] * DH2O) +
                idv.h2osoi_ice[c, k + joff] / (idv.col_dz[c, k + joff] * DICE)

            if odv.h2osoi_vol[c, k] / idv.watsat[c, k] <= sat_lev
                k_zwt = k
                sat_flag = 0
                break
            end
        end
        if sat_flag == 1
            k_zwt = 1
        end

        if k_zwt == 1
            odv.zwt[c] = idv.col_zi[c, 1 + joff_zi]
        elseif k_zwt < nbed
            s1 = (idv.h2osoi_liq[c, k_zwt + joff] / (idv.col_dz[c, k_zwt + joff] * DH2O) +
                idv.h2osoi_ice[c, k_zwt + joff] / (idv.col_dz[c, k_zwt + joff] * DICE)) / idv.watsat[c, k_zwt]
            s2 = (idv.h2osoi_liq[c, k_zwt+1 + joff] / (idv.col_dz[c, k_zwt+1 + joff] * DH2O) +
                idv.h2osoi_ice[c, k_zwt+1 + joff] / (idv.col_dz[c, k_zwt+1 + joff] * DICE)) / idv.watsat[c, k_zwt+1]

            m_val = (idv.col_z[c, k_zwt+1 + joff] - idv.col_z[c, k_zwt + joff]) / (s2 - s1)
            b_val = idv.col_z[c, k_zwt+1 + joff] - m_val * s2
            odv.zwt[c] = max(zr, m_val * sat_lev + b_val)
        else
            odv.zwt[c] = idv.col_zi[c, nbed + joff_zi]
        end
    end
end

function soilhyd_theta_based_water_table!(odv, idv, mask, nlevsoi::Int, joff::Int, joff_zi::Int, denh2o, denice)
    T = eltype(odv.zwt)
    _launch!(_theta_based_water_table_kernel!, odv.zwt, odv, idv, mask, nlevsoi,
             joff, joff_zi, T(denh2o), T(denice); ndrange = length(mask))
end

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
    col_dz::AbstractMatrix{<:Real},
    col_z::AbstractMatrix{<:Real},
    col_zi::AbstractMatrix{<:Real},
    col_nbedrock::AbstractVector{<:Integer},
    mask_hydrology::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsoi::Int;
    joff::Int=varpar.nlevsno,
    joff_zi::Int=varpar.nlevsno + 1
)
    idv = _TBWTInDV(;
        col_dz       = col_dz,
        col_z        = col_z,
        col_zi       = col_zi,
        watsat       = soilstate.watsat_col,
        h2osoi_liq   = waterstatebulk_ws.h2osoi_liq_col,
        h2osoi_ice   = waterstatebulk_ws.h2osoi_ice_col,
        col_nbedrock = col_nbedrock)

    odv = _TBWTOutDV(;
        zwt        = soilhydrology.zwt_col,
        h2osoi_vol = waterstatebulk_ws.h2osoi_vol_col)

    soilhyd_theta_based_water_table!(odv, idv, mask_hydrology, nlevsoi, joff, joff_zi, DENH2O, DENICE)

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
# ---- calc_irrig_withdrawals! : per-column layered groundwater withdrawal ----
# Each column is independent: a per-layer init, a jwt search, then a SEQUENTIAL
# layered withdrawal cascade carrying loop-carried `irrig_demand_remaining`. One
# thread per column runs the whole subroutine in-thread and writes only its own
# column (the matrix init folds into the kernel's first sub-loop). Scalar args
# (dtime, aq_sp_yield_min) carry the working eltype and literals in arithmetic
# are eltype-converted so no Float64 reaches a Float32-only backend (Metal);
# byte-identical on a Float64 CPU run.
@kernel function _soilhyd_irrig_withdrawals_kernel!(qflx_gw_uncon_irrig_lyr,
        qflx_gw_con_irrig, @Const(mask), @Const(qflx_gw_demand), @Const(col_nbedrock),
        @Const(col_zi), @Const(zwt), @Const(watsat), @Const(sucsat), @Const(bsw),
        nlevsoi::Int, dtime, aq_sp_yield_min)
    c = @index(Global)
    @inbounds if mask[c]
        T   = eltype(qflx_gw_uncon_irrig_lyr)
        zr  = zero(T)
        on  = one(T)
        e3  = T(1.0e3)
        dt  = T(dtime)
        aqy = T(aq_sp_yield_min)

        # Initialize this column's layers
        for j in 1:nlevsoi
            qflx_gw_uncon_irrig_lyr[c, j] = zr
        end

        irrig_demand_remaining = qflx_gw_demand[c] * dt

        # (negative-demand error is a host-only pre-check; no error() in-kernel)

        jwt_c = nlevsoi
        for j in 1:nlevsoi
            if zwt[c] <= col_zi[c, j]
                jwt_c = j - 1
                break
            end
        end

        for j in (jwt_c+1):col_nbedrock[c]
            s_y = watsat[c, j] *
                (on - (on + e3 * zwt[c] / sucsat[c, j])^(-on / bsw[c, j]))
            s_y = max(s_y, aqy)

            if j == jwt_c + 1
                available_water_layer = max(zr, s_y * (col_zi[c, j] - zwt[c]) * e3)
            else
                available_water_layer = max(zr, s_y * (col_zi[c, j] - col_zi[c, j-1]) * e3)
            end

            irrig_layer = min(irrig_demand_remaining, available_water_layer)
            qflx_gw_uncon_irrig_lyr[c, j] = irrig_layer / dt

            irrig_demand_remaining = irrig_demand_remaining - irrig_layer

            if irrig_demand_remaining <= zr
                break
            end
        end

        if irrig_demand_remaining > zr
            qflx_gw_con_irrig[c] = irrig_demand_remaining / dt
        else
            qflx_gw_con_irrig[c] = zr
        end
    end
end

function calc_irrig_withdrawals!(
    soilhydrology::SoilHydrologyData,
    soilstate::SoilStateData,
    qflx_gw_demand::AbstractVector{<:Real},
    qflx_gw_uncon_irrig_lyr::AbstractMatrix{<:Real},
    qflx_gw_con_irrig::AbstractVector{<:Real},
    col_nbedrock::AbstractVector{<:Integer},
    col_zi::AbstractMatrix{<:Real},
    mask_soil::AbstractVector{Bool},
    bounds::UnitRange{Int},
    nlevsoi::Int,
    dtime::Real
)
    params = soilhydrology_params
    bsw    = soilstate.bsw_col
    sucsat = soilstate.sucsat_col
    watsat = soilstate.watsat_col
    zwt    = soilhydrology.zwt_col

    # Host-only pre-check for negative demand. The CPU path keeps the original
    # error() (the loop is over host Arrays); on a device backend this scalar
    # scan is skipped (valid forcing only), so the kernel stays scalar-index-free.
    if qflx_gw_demand isa Array
        for c in bounds
            mask_soil[c] || continue
            if qflx_gw_demand[c] * dtime < 0.0
                error("negative groundwater irrigation demand at c=$c!")
            end
        end
    end

    T = eltype(qflx_gw_uncon_irrig_lyr)
    _launch!(_soilhyd_irrig_withdrawals_kernel!, qflx_gw_uncon_irrig_lyr,
             qflx_gw_con_irrig, mask_soil, qflx_gw_demand, col_nbedrock, col_zi,
             zwt, watsat, sucsat, bsw, nlevsoi, T(dtime), T(params.aq_sp_yield_min);
             ndrange = length(mask_soil))

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
    mask_soil::AbstractVector{Bool},
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

# --------------------------------------------------------------------------
# perched_lateral_flow! : ONE thread per (hydrology) column for the
# NON-HILLSLOPE path. That path is fully per-column independent (no cross-
# column scatter — `cold[c] == ISPVAL` so the hillslope donor/receiver
# branch is never entered), so the whole subroutine maps to one thread per
# column running the in-thread frost/perch-table search, the q_perch
# drainage, and the per-layer storage removal. The hillslope path carries a
# downhill scatter (`qflx_drain_perched[cold[c]] -= …`) that needs exact
# ordering, so it stays on the sequential CPU branch below. Device arrays
# (Metal/CUDA/…) take the kernel; CPU arrays take the byte-identical scalar
# body verbatim.
#
# Soil-column arrays are grouped into two immutable device-view bundles
# (Adapt.@adapt_structure'd); field names mirror the original Julia locals so
# the kernel body reads like the scalar loops.
# --------------------------------------------------------------------------
Base.@kwdef struct _PLInDV{V,M,VI}
    col_dz::M; col_zi::M
    hksat::M; sucsat::M; watsat::M; bsw::M
    topo_slope::V; nbedrock::VI
    frost_table::V; zwt_perched::V
end
Adapt.@adapt_structure _PLInDV

Base.@kwdef struct _PLOutDV{V,M}
    qflx_drain_perched::V; h2osoi_liq::M
end
Adapt.@adapt_structure _PLOutDV

@kernel function _perched_lateral_flow_kernel!(qfdp_out, odv, @Const(idv), @Const(mask),
        nlevsoi::Int, joff::Int, joff_zi::Int, dtime, perched_baseflow_scalar,
        aq_sp_yield_min, rpi)
    c = @index(Global)
    @inbounds if mask[c]
        T   = eltype(qfdp_out)
        zr  = zero(T)
        on  = one(T)
        e3  = T(1.0e3)
        thou = T(1000.0)
        aqy = T(aq_sp_yield_min)
        nb  = idv.nbedrock[c]

        # ---- locate frost table and perched water table layers ----
        k_frost = nb
        for k in 1:nb
            zi_km1 = k > 1 ? idv.col_zi[c, k - 1 + joff_zi] : zr
            if idv.frost_table[c] >= zi_km1 && idv.frost_table[c] < idv.col_zi[c, k + joff_zi]
                k_frost = k
                break
            end
        end
        k_perch = nb
        for k in 1:nb
            zi_km1 = k > 1 ? idv.col_zi[c, k - 1 + joff_zi] : zr
            if idv.zwt_perched[c] >= zi_km1 && idv.zwt_perched[c] < idv.col_zi[c, k + joff_zi]
                k_perch = k
                break
            end
        end

        # ---- drainage from perched saturated region (non-hillslope) ----
        qfdp = zr
        out  = zr
        if idv.frost_table[c] > idv.zwt_perched[c]
            q_perch_max = T(perched_baseflow_scalar) * sin(idv.topo_slope[c] * (T(rpi) / T(180.0)))
            wtsub = zr
            q_perch = zr
            for k in k_perch:(k_frost - 1)
                q_perch += idv.hksat[c, k] * idv.col_dz[c, k + joff]
                wtsub += idv.col_dz[c, k + joff]
            end
            if wtsub > zr
                q_perch = q_perch / wtsub
            end
            out = q_perch_max * q_perch * (idv.frost_table[c] - idv.zwt_perched[c])
        end

        # ---- net drainage (no scatter on the non-hillslope path) ----
        qfdp += out

        # ---- remove drainage from soil moisture storage ----
        drainage_tot = qfdp * dtime
        for k in k_perch:(k_frost - 1)
            s_y = idv.watsat[c, k] *
                (on - (on + e3 * idv.zwt_perched[c] / idv.sucsat[c, k])^(-on / idv.bsw[c, k]))
            s_y = max(s_y, aqy)

            if k == k_perch
                drainage_layer = min(drainage_tot, s_y * (idv.col_zi[c, k + joff_zi] - idv.zwt_perched[c]) * e3)
            else
                drainage_layer = min(drainage_tot, s_y * idv.col_dz[c, k + joff] * e3)
            end
            drainage_layer = max(drainage_layer, zr)
            drainage_tot -= drainage_layer
            odv.h2osoi_liq[c, k + joff] -= drainage_layer
        end

        qfdp -= drainage_tot / dtime
        odv.qflx_drain_perched[c] = qfdp
    end
end

function soilhyd_perched_lateral_flow!(odv, idv, mask, nlevsoi::Int, joff::Int, joff_zi::Int,
                                       dtime, perched_baseflow_scalar, aq_sp_yield_min, rpi)
    T = eltype(odv.qflx_drain_perched)
    _launch!(_perched_lateral_flow_kernel!, odv.qflx_drain_perched, odv, idv, mask,
             nlevsoi, joff, joff_zi, T(dtime), T(perched_baseflow_scalar),
             T(aq_sp_yield_min), T(rpi); ndrange = length(mask))
end

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
    tdepth_grc::AbstractVector{<:Real},
    tdepthmax_grc::AbstractVector{<:Real},
    mask_hydrology::AbstractVector{Bool},
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

    # Device path (Metal/CUDA/…): the non-hillslope, per-column-independent
    # subroutine runs as one thread per column. The hillslope donor/receiver
    # scatter is not expressible as a race-free per-column kernel, so device
    # execution is restricted to the non-hillslope case (asserted below). The
    # CPU branch (KA.CPU backend) falls through to the byte-identical scalar
    # body that handles both paths.
    if !(KA.get_backend(soilhydrology.zwt_perched_col) isa KA.CPU)
        any(col_data.is_hillslope_column) &&
            error("perched_lateral_flow! device kernel supports non-hillslope columns only")
        idv = _PLInDV(;
            col_dz = col_data.dz, col_zi = col_data.zi,
            hksat = soilstate.hksat_col, sucsat = soilstate.sucsat_col,
            watsat = soilstate.watsat_col, bsw = soilstate.bsw_col,
            topo_slope = col_data.topo_slope, nbedrock = col_data.nbedrock,
            frost_table = soilhydrology.frost_table_col,
            zwt_perched = soilhydrology.zwt_perched_col)
        odv = _PLOutDV(;
            qflx_drain_perched = wf.qflx_drain_perched_col,
            h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col)
        soilhyd_perched_lateral_flow!(odv, idv, mask_hydrology, nlevsoi, joff, joff_zi,
            dtime, params.perched_baseflow_scalar, params.aq_sp_yield_min, RPI)
        return nothing
    end

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
                            max(min(stream_water_depth - stream_channel_depth, 0.0),
                                col_hill_elev[c] - frost_table[c])
                        head_gradient = head_gradient / col_hill_distance[c]

                        if stream_water_depth <= 0.0
                            head_gradient = max(head_gradient, 0.0)
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
            s_y = max(s_y, params.aq_sp_yield_min)

            if k == k_perch_arr[c]
                drainage_layer = min(drainage_tot, s_y * (col_zi[c, k + joff_zi] - zwt_perched[c]) * 1.0e3)
            else
                drainage_layer = min(drainage_tot, s_y * col_dz[c, k + joff] * 1.0e3)
            end

            drainage_layer = max(drainage_layer, 0.0)
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

# --------------------------------------------------------------------------
# subsurface_lateral_flow! : ONE thread per (hydrology) column for the
# NON-HILLSLOPE path. As with perched_lateral_flow!, that path is fully
# per-column independent (`cold[c] == ISPVAL` so the hillslope inflow scatter
# `qflx_latflow_in[cold[c]] += …` is never entered). Each thread runs the full
# per-column sequence in-thread: icefrac/impedance, jwt search, power-law
# baseflow, topographic-runoff water-table rise/deepen cascade (loop-carried
# zwt), the upward excess redistribution, the surface/ice overflow, and the
# watmin back-fill — writing only its own column. The hillslope downhill
# scatter needs exact ordering and stays on the sequential CPU branch.
#
# All soil-column arrays are grouped into Adapt'd device-view bundles whose
# field names mirror the scalar locals so the kernel body reads verbatim.
# --------------------------------------------------------------------------
Base.@kwdef struct _SLInDV{V,M,VI}
    col_dz::M; col_zi::M
    hksat::M; sucsat::M; watsat::M; bsw::M; eff_porosity::M
    topo_slope::V; wtgcell::V; nbedrock::VI
    itype::VI; landunit::VI; gridcell::VI
    grc_area::V; frost_table::V
    qflx_snwcp_liq::V
end
Adapt.@adapt_structure _SLInDV

Base.@kwdef struct _SLOutDV{V,M}
    zwt::V; h2osfc::V
    h2osoi_liq::M; h2osoi_ice::M; icefrac::M
    qflx_latflow_out::V; qflx_latflow_in::V; volumetric_discharge::V
    qflx_drain::V; qflx_qrgwl::V; qflx_rsub_sat::V; qflx_ice_runoff_xs::V
end
Adapt.@adapt_structure _SLOutDV

@kernel function _subsurface_lateral_flow_kernel!(zwt_out, odv, @Const(idv),
        @Const(mask), @Const(mask_urban), @Const(urbpoi),
        nlevsoi::Int, joff::Int, joff_zi::Int, dtime,
        e_ice, n_baseflow, aq_sp_yield_min, baseflow_scalar, rpi,
        denice, watmin, pondmx, icol_road_perv)
    c = @index(Global)
    @inbounds if mask[c]
        T    = eltype(zwt_out)
        zr   = zero(T)
        on   = one(T)
        e3   = T(1.0e3)
        em3  = T(1.0e-3)
        e6   = T(1.0e6)
        thou = T(1000.0)
        aqy  = T(aq_sp_yield_min)
        eice = T(e_ice)
        DICE = T(denice)
        WMIN = T(watmin)
        PMX  = T(pondmx)
        nb   = idv.nbedrock[c]
        l    = idv.landunit[c]
        g    = idv.gridcell[c]

        # ---- icefrac / impedance, column-average impedance, jwt ----
        zi_bedrock = idv.col_zi[c, nb + joff_zi]
        jwt = nlevsoi
        for j in 1:nlevsoi
            if odv.zwt[c] <= idv.col_zi[c, j + joff_zi]
                jwt = j - 1
                break
            end
        end

        dzsum = zr
        icefracsum = zr
        for j in max(jwt, 1):nlevsoi
            dzmm_j = idv.col_dz[c, j + joff] * e3
            vol_ice_val = min(idv.watsat[c, j], odv.h2osoi_ice[c, j + joff] / (idv.col_dz[c, j + joff] * DICE))
            ic = min(on, vol_ice_val / idv.watsat[c, j])
            odv.icefrac[c, j] = ic
            dzsum += dzmm_j
            icefracsum += ic * dzmm_j
        end
        ice_imped_col = T(10.0)^(-eice * (icefracsum / dzsum))

        # ---- init fluxes ----
        odv.qflx_drain[c]            = zr
        odv.qflx_rsub_sat[c]         = zr
        odv.qflx_latflow_in[c]       = zr
        odv.qflx_latflow_out[c]      = zr
        odv.volumetric_discharge[c]  = zr
        drainage = zr

        # ---- power-law baseflow (non-hillslope) ----
        if odv.zwt[c] <= zi_bedrock
            odv.qflx_latflow_out[c] = ice_imped_col * T(baseflow_scalar) *
                tan(T(rpi) / T(180.0) * idv.topo_slope[c]) *
                (zi_bedrock - odv.zwt[c])^T(n_baseflow)
        end
        qflx_latflow_out_vol = em3 * odv.qflx_latflow_out[c] *
            (idv.grc_area[g] * e6 * idv.wtgcell[c])
        odv.volumetric_discharge[c] = qflx_latflow_out_vol

        # ---- topographic runoff: remove water via drainage ----
        qflx_net_latflow = odv.qflx_latflow_out[c] - odv.qflx_latflow_in[c]
        drainage = odv.zwt[c] <= zi_bedrock ? qflx_net_latflow : zr

        drainage_tot = -drainage * dtime
        if drainage_tot > zr  # rising water table
            for j in (jwt + 1):-1:1
                if idv.col_zi[c, j + joff_zi] < idv.frost_table[c]
                    s_y = idv.watsat[c, j] *
                        (on - (on + e3 * odv.zwt[c] / idv.sucsat[c, j])^(-on / idv.bsw[c, j]))
                    s_y = max(s_y, aqy)
                    drainage_layer = min(drainage_tot, s_y * idv.col_dz[c, j + joff] * e3)
                    drainage_layer = max(drainage_layer, zr)
                    odv.h2osoi_liq[c, j + joff] += drainage_layer
                    drainage_tot -= drainage_layer
                    if drainage_tot <= zr
                        odv.zwt[c] -= drainage_layer / s_y / thou
                        break
                    else
                        odv.zwt[c] = j > 1 ? idv.col_zi[c, j - 1 + joff_zi] : zr
                    end
                end
            end
            odv.h2osfc[c] += drainage_tot
        else  # deepening water table
            for j in (jwt + 1):nb
                s_y = idv.watsat[c, j] *
                    (on - (on + e3 * odv.zwt[c] / idv.sucsat[c, j])^(-on / idv.bsw[c, j]))
                s_y = max(s_y, aqy)
                drainage_layer = max(drainage_tot, -(s_y * (idv.col_zi[c, j + joff_zi] - odv.zwt[c]) * e3))
                drainage_layer = min(drainage_layer, zr)
                odv.h2osoi_liq[c, j + joff] += drainage_layer
                drainage_tot -= drainage_layer
                if drainage_tot >= zr
                    odv.zwt[c] -= drainage_layer / s_y / thou
                    break
                else
                    odv.zwt[c] = idv.col_zi[c, j + joff_zi]
                end
            end
            drainage += drainage_tot / dtime
        end
        odv.zwt[c] = max(zr, odv.zwt[c])
        odv.zwt[c] = min(T(80.0), odv.zwt[c])

        # ---- excessive water above saturation: redistribute upward ----
        # HARD max/min, as a PAIR. Both clamp a layer water MASS [mm] against 0 / the layer's
        # capacity — constant-or-independent branches, so the derivative on the clamped side is
        # zero and smoothing recovers nothing. Two reasons this pair had to be hardened:
        #  1. the amount MOVED was wrong by log(2)/50 = 0.0139 mm per layer per step, and in a
        #     saturated column (where the ReLU sits exactly at its kink, liq - capacity == 0) that
        #     is ~0.2 mm/step pumped up into layer 1 and straight out as saturation-excess runoff;
        #  2. it is INVISIBLE to the water balance, because max(x-cap,0) and min(cap,x) carry
        #     equal and opposite errors that cancel in the sum — total mass is conserved exactly
        #     while the water is put in the wrong place. Keeping them as one hard pair preserves
        #     the identity min(cap,x) + max(x-cap,0) == x by construction.
        for j in nlevsoi:-1:2
            dzmm_j = idv.col_dz[c, j + joff] * e3
            xsi = max(odv.h2osoi_liq[c, j + joff] - idv.eff_porosity[c, j] * dzmm_j, zr)
            odv.h2osoi_liq[c, j + joff] = min(idv.eff_porosity[c, j] * dzmm_j, odv.h2osoi_liq[c, j + joff])
            odv.h2osoi_liq[c, j - 1 + joff] += xsi
        end

        # ---- top-layer surface / ice overflow ----
        dzmm_1 = idv.col_dz[c, 1 + joff] * e3
        xs1 = max(max(odv.h2osoi_liq[c, 1 + joff] - WMIN, zr) -
            max(zr, PMX + idv.watsat[c, 1] * dzmm_1 - odv.h2osoi_ice[c, 1 + joff] - WMIN), zr)
        odv.h2osoi_liq[c, 1 + joff] -= xs1
        if urbpoi[l]
            odv.qflx_rsub_sat[c] = xs1 / dtime
        else
            odv.h2osfc[c] += xs1
            odv.qflx_rsub_sat[c] = zr
        end

        xs1 = max(max(odv.h2osoi_ice[c, 1 + joff], zr) -
            max(zr, PMX + idv.watsat[c, 1] * dzmm_1 - odv.h2osoi_liq[c, 1 + joff]), zr)
        odv.h2osoi_ice[c, 1 + joff] = min(max(zr, PMX + idv.watsat[c, 1] * dzmm_1 - odv.h2osoi_liq[c, 1 + joff]),
            odv.h2osoi_ice[c, 1 + joff])
        odv.qflx_ice_runoff_xs[c] = xs1 / dtime

        # ---- limit h2osoi_liq >= watmin (top nlevsoi-1 layers cascade down) ----
        for j in 1:(nlevsoi - 1)
            if odv.h2osoi_liq[c, j + joff] < WMIN
                xs = WMIN - odv.h2osoi_liq[c, j + joff]
                if j == jwt
                    odv.zwt[c] += xs / idv.eff_porosity[c, j] / thou
                end
            else
                xs = zr
            end
            odv.h2osoi_liq[c, j + joff]     += xs
            odv.h2osoi_liq[c, j + 1 + joff] -= xs
        end

        # ---- bottom layer: pull water from above ----
        jb = nlevsoi
        if odv.h2osoi_liq[c, jb + joff] < WMIN
            xs = WMIN - odv.h2osoi_liq[c, jb + joff]
            for i in (nlevsoi - 1):-1:1
                avail = max(odv.h2osoi_liq[c, i + joff] - WMIN - xs, zr)
                if avail >= xs
                    odv.h2osoi_liq[c, jb + joff] += xs
                    odv.h2osoi_liq[c, i + joff]  -= xs
                    xs = zr
                    break
                else
                    odv.h2osoi_liq[c, jb + joff] += avail
                    odv.h2osoi_liq[c, i + joff]  -= avail
                    xs -= avail
                end
            end
        else
            xs = zr
        end
        odv.h2osoi_liq[c, jb + joff] += xs
        odv.qflx_rsub_sat[c] -= xs / dtime

        # ---- final fluxes ----
        odv.qflx_drain[c] = odv.qflx_rsub_sat[c] + drainage
        odv.qflx_qrgwl[c] = idv.qflx_snwcp_liq[c]

        # ---- urban: no drainage except pervious road ----
        if mask_urban[c] && idv.itype[c] != icol_road_perv
            odv.qflx_drain[c] = zr
            odv.qflx_qrgwl[c] = idv.qflx_snwcp_liq[c]
        end
    end
end

function soilhyd_subsurface_lateral_flow!(odv, idv, mask, mask_urban, urbpoi,
        nlevsoi::Int, joff::Int, joff_zi::Int, dtime,
        e_ice, n_baseflow, aq_sp_yield_min, baseflow_scalar, rpi,
        denice, watmin, pondmx, icol_road_perv)
    T = eltype(odv.zwt)
    _launch!(_subsurface_lateral_flow_kernel!, odv.zwt, odv, idv, mask, mask_urban, urbpoi,
             nlevsoi, joff, joff_zi, T(dtime), T(e_ice), T(n_baseflow), T(aq_sp_yield_min),
             T(baseflow_scalar), T(rpi), T(denice), T(watmin), T(pondmx), icol_road_perv;
             ndrange = length(mask))
end

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
    tdepth_grc::AbstractVector{<:Real},
    tdepthmax_grc::AbstractVector{<:Real},
    grc_area::AbstractVector{<:Real},
    mask_hydrology::AbstractVector{Bool},
    mask_urban::AbstractVector{Bool},
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

    # Device path (Metal/CUDA/…): one thread per column for the non-hillslope,
    # per-column-independent subroutine. The hillslope inflow scatter
    # (`qflx_latflow_in[cold[c]] += …`) needs exact ordering and is not a
    # race-free per-column kernel, so device execution is restricted to the
    # non-hillslope case. CPU arrays fall through to the byte-identical scalar
    # body that handles both paths.
    if !(KA.get_backend(soilhydrology.zwt_col) isa KA.CPU)
        any(col_data.is_hillslope_column) &&
            error("subsurface_lateral_flow! device kernel supports non-hillslope columns only")
        idv = _SLInDV(;
            col_dz = col_data.dz, col_zi = col_data.zi,
            hksat = soilstate.hksat_col, sucsat = soilstate.sucsat_col,
            watsat = soilstate.watsat_col, bsw = soilstate.bsw_col,
            eff_porosity = soilstate.eff_porosity_col,
            topo_slope = col_data.topo_slope, wtgcell = col_data.wtgcell,
            nbedrock = col_data.nbedrock, itype = col_data.itype,
            landunit = col_data.landunit, gridcell = col_data.gridcell,
            grc_area = grc_area, frost_table = soilhydrology.frost_table_col,
            qflx_snwcp_liq = wf.qflx_snwcp_liq_col)
        odv = _SLOutDV(;
            zwt = soilhydrology.zwt_col, h2osfc = waterstatebulk_ws.h2osfc_col,
            h2osoi_liq = waterstatebulk_ws.h2osoi_liq_col,
            h2osoi_ice = waterstatebulk_ws.h2osoi_ice_col,
            icefrac = soilhydrology.icefrac_col,
            qflx_latflow_out = wf.qflx_latflow_out_col,
            qflx_latflow_in = wf.qflx_latflow_in_col,
            volumetric_discharge = wf.volumetric_discharge_col,
            qflx_drain = wf.qflx_drain_col, qflx_qrgwl = wf.qflx_qrgwl_col,
            qflx_rsub_sat = wf.qflx_rsub_sat_col,
            qflx_ice_runoff_xs = wf.qflx_ice_runoff_xs_col)
        soilhyd_subsurface_lateral_flow!(odv, idv, mask_hydrology, mask_urban,
            lun_data.urbpoi, nlevsoi, joff, joff_zi, dtime,
            params.e_ice, params.n_baseflow, params.aq_sp_yield_min,
            BASEFLOW_SCALAR[], RPI, DENICE, WATMIN, PONDMX, ICOL_ROAD_PERV)
        return nothing
    end

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

            # Axis-scaled — see SOIL_HYDRAULIC_K. icefrac is the exponent of the ice impedance, so
            # a 0.0139 error in it is a factor 10^(e_ice·0.0139) ≈ 1.2-1.6x error in permeability.
            vol_ice_val = smooth_min(watsat[c, j], h2osoi_ice[c, j + joff] / (col_dz[c, j + joff] * DENICE);
                                     k = SOIL_HYDRAULIC_K)
            icefrac[c, j] = smooth_min(1.0, vol_ice_val / watsat[c, j]; k = SOIL_HYDRAULIC_K)
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
                        min(stream_water_depth - stream_channel_depth, 0.0)
                    head_gradient /= col_hill_distance[c]

                    if stream_water_depth <= 0.0
                        head_gradient = max(head_gradient, 0.0)
                    end

                    if head_gradient < 0.0
                        head_gradient -= 1.0 / k_anisotropic
                    end
                end
            else
                error("head_gradient_method must be kinematic or darcy")
            end

            # Cap maximum head_gradient
            head_gradient = min(max(head_gradient, -2.0), 2.0)

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
                    s_y = max(s_y, params.aq_sp_yield_min)

                    drainage_layer = min(drainage_tot, s_y * col_dz[c, j + joff] * 1.0e3)
                    drainage_layer = max(drainage_layer, 0.0)
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
                s_y = max(s_y, params.aq_sp_yield_min)

                drainage_layer = max(drainage_tot, -(s_y * (col_zi[c, j + joff_zi] - zwt[c]) * 1.0e3))
                drainage_layer = min(drainage_layer, 0.0)
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

        zwt[c] = max(0.0, zwt[c])
        zwt[c] = min(80.0, zwt[c])
    end

    # ---- Excessive water above saturation: redistribute upward ----
    for j in nlevsoi:-1:2
        for c in bounds
            mask_hydrology[c] || continue
            xsi_arr[c]        = max(h2osoi_liq[c, j + joff] - eff_porosity[c, j] * dzmm[c, j], 0.0)
            h2osoi_liq[c, j + joff]  = min(eff_porosity[c, j] * dzmm[c, j], h2osoi_liq[c, j + joff])
            h2osoi_liq[c, j - 1 + joff] += xsi_arr[c]
        end
    end

    for c in bounds
        mask_hydrology[c] || continue

        xs1_arr[c] = max(max(h2osoi_liq[c, 1 + joff] - WATMIN, 0.0) -
            max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_ice[c, 1 + joff] - WATMIN), 0.0)
        h2osoi_liq[c, 1 + joff] -= xs1_arr[c]

        l = col_landunit[c]
        if lun_urbpoi[l]
            qflx_rsub_sat[c] = xs1_arr[c] / dtime
        else
            h2osfc[c] += xs1_arr[c]
            qflx_rsub_sat[c] = 0.0
        end

        # ice check
        xs1_arr[c] = max(max(h2osoi_ice[c, 1 + joff], 0.0) -
            max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_liq[c, 1 + joff]), 0.0)
        h2osoi_ice[c, 1 + joff] = min(max(0.0, PONDMX + watsat[c, 1] * dzmm[c, 1] - h2osoi_liq[c, 1 + joff]),
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
                available_h2osoi_liq = max(h2osoi_liq[c, i + joff] - WATMIN - xs_arr[c], 0.0)
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
