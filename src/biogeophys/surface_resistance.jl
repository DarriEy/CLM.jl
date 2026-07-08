# ==========================================================================
# Ported from: src/biogeophys/SurfaceResistanceMod.F90 (446 lines)
# Surface resistance calculations for soil evaporation.
# Computes either a beta factor (Lee-Pielke 1992) or a soil resistance
# (Swenson & Lawrence 2014) depending on the selected method.
# ==========================================================================

# --- Method selection constants ---
const SOIL_RESIS_LEEPIELKE_1992 = 0
const SOIL_RESIS_SL_14 = 1

# --- Module-level state (equivalent to Fortran module variable) ---
Base.@kwdef mutable struct SurfaceResistanceControl
    soil_resis_method::Int = SOIL_RESIS_SL_14  # default method
end

const surface_resistance_ctrl = SurfaceResistanceControl()

# --- Parameters read from parameter file ---
Base.@kwdef mutable struct SurfaceResistanceParams
    d_max::Float64 = 15.0               # Dry surface layer parameter (mm)
    frac_sat_soil_dsl_init::Float64 = 0.8  # Fraction of saturated soil for moisture value at which DSL initiates (unitless)
end

const surface_resistance_params = SurfaceResistanceParams()

"""
    soil_resistance_read_nl!(; soil_resis_method=SOIL_RESIS_SL_14)

Set the soil resistance method. In Julia, namelist values are passed
directly instead of reading from a file.

Ported from `soil_resistance_readNL` in `SurfaceResistanceMod.F90`.
"""
function soil_resistance_read_nl!(; soil_resis_method::Int = SOIL_RESIS_SL_14)
    if soil_resis_method != SOIL_RESIS_LEEPIELKE_1992 && soil_resis_method != SOIL_RESIS_SL_14
        error("soil_resistance_read_nl!: invalid soil_resis_method = $soil_resis_method")
    end
    surface_resistance_ctrl.soil_resis_method = soil_resis_method
    return nothing
end

"""
    surface_resistance_read_params!(; d_max, frac_sat_soil_dsl_init)

Set the surface resistance parameters. In Julia, parameter values are passed
directly instead of reading from a NetCDF file.

Ported from `readParams` in `SurfaceResistanceMod.F90`.
"""
function surface_resistance_read_params!(; d_max::Real = 15.0,
                                           frac_sat_soil_dsl_init::Real = 0.8)
    surface_resistance_params.d_max = d_max
    surface_resistance_params.frac_sat_soil_dsl_init = frac_sat_soil_dsl_init
    return nothing
end

"""
    do_soilevap_beta() -> Bool

Return true if the moisture stress for soil evaporation is computed as a
beta factor (Lee-Pielke 1992 method).

Ported from `do_soilevap_beta` in `SurfaceResistanceMod.F90`.
"""
function do_soilevap_beta()
    return surface_resistance_ctrl.soil_resis_method == SOIL_RESIS_LEEPIELKE_1992
end

"""
    do_soil_resistance_sl14() -> Bool

Return true if the soil evaporative resistance is computed using a
dry surface layer (Swenson & Lawrence 2014 method).

Ported from `do_soil_resistance_sl14` in `SurfaceResistanceMod.F90`.
"""
function do_soil_resistance_sl14()
    return surface_resistance_ctrl.soil_resis_method == SOIL_RESIS_SL_14
end

"""
    calc_soilevap_resis!(col, lun, soilstate, waterstatebulk, waterdiagbulk,
                         temperature, mask_nolakec, bounds)

Compute the resistance factor for soil evaporation. Dispatches to either
the Lee-Pielke 1992 beta factor or the Swenson & Lawrence 2014 soil
resistance depending on `surface_resistance_ctrl.soil_resis_method`.

# Arguments
- `col::ColumnData`: column data (dz, landunit, itype)
- `lun::LandunitData`: landunit data (itype)
- `soilstate::SoilStateData`: soil state (output: soilbeta_col, dsl_col, soilresis_col;
   input: watsat_col, watfc_col, bsw_col, sucsat_col)
- `waterstatebulk::WaterStateBulkData`: water state (h2osoi_ice_col, h2osoi_liq_col via ws)
- `waterdiagbulk::WaterDiagnosticBulkData`: water diagnostics (frac_sno_col, frac_h2osfc_col)
- `temperature::TemperatureData`: temperature state (t_soisno_col)
- `mask_nolakec::AbstractVector{Bool}`: non-lake column mask
- `bounds::UnitRange{Int}`: column bounds

Ported from `calc_soilevap_resis` in `SurfaceResistanceMod.F90`.
"""
function calc_soilevap_resis!(col::ColumnData,
                               lun::LandunitData,
                               soilstate::SoilStateData,
                               waterstatebulk::WaterStateBulkData,
                               waterdiagbulk::WaterDiagnosticBulkData,
                               temperature::TemperatureData,
                               mask_nolakec::AbstractVector{Bool},
                               bounds::UnitRange{Int})
    if surface_resistance_ctrl.soil_resis_method == SOIL_RESIS_LEEPIELKE_1992
        calc_beta_leepielke1992!(col, lun, soilstate, waterstatebulk,
                                  waterdiagbulk, mask_nolakec, bounds)
    elseif surface_resistance_ctrl.soil_resis_method == SOIL_RESIS_SL_14
        calc_soil_resistance_sl14!(col, lun, soilstate, waterstatebulk,
                                    temperature, mask_nolakec, bounds)
    else
        error("calc_soilevap_resis!: a soilevap resis function must be specified!")
    end
    return nothing
end

# --------------------------------------------------------------------------
# Kernel: Lee-Pielke (1992) beta factor (per column, masked). Each column
# writes only soilbeta[c]; fully independent. Literals eltype-converted so no
# Float64 reaches a Float32-only backend (Metal); byte-identical on Float64.
# --------------------------------------------------------------------------
@kernel function _beta_leepielke_kernel!(soilbeta, @Const(mask), @Const(landunit),
                                         @Const(col_itype), @Const(lun_itype),
                                         @Const(dz), @Const(watsat), @Const(watfc),
                                         @Const(h2osoi_ice), @Const(h2osoi_liq),
                                         @Const(frac_sno), @Const(frac_h2osfc),
                                         joff::Int, denh2o, denice)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(soilbeta)
        l = landunit[c]
        if lun_itype[l] != ISTWET && lun_itype[l] != ISTICE
            if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
                wx  = (h2osoi_liq[c, 1 + joff] / denh2o + h2osoi_ice[c, 1 + joff] / denice) / dz[c, 1 + joff]
                if wx < watfc[c, 1]
                    fac_fc = min(one(T), wx / watfc[c, 1])
                    fac_fc = max(fac_fc, T(0.01))
                    soilbeta[c] = (one(T) - frac_sno[c] - frac_h2osfc[c]) *
                                  T(0.25) * (one(T) - cos(T(π) * fac_fc))^T(2.0) +
                                  frac_sno[c] + frac_h2osfc[c]
                else
                    soilbeta[c] = one(T)
                end
            elseif col_itype[c] == ICOL_ROAD_PERV
                soilbeta[c] = zero(T)
            elseif col_itype[c] == ICOL_SUNWALL || col_itype[c] == ICOL_SHADEWALL
                soilbeta[c] = zero(T)
            elseif col_itype[c] == ICOL_ROOF || col_itype[c] == ICOL_ROAD_IMPERV
                soilbeta[c] = zero(T)
            end
        else
            soilbeta[c] = one(T)
        end
    end
end

"""
    calc_beta_leepielke1992!(col, lun, soilstate, waterstatebulk, waterdiagbulk,
                              mask_nolakec, bounds)

Compute the Lee-Pielke (1992) beta factor to scale actual soil evaporation
from potential evaporation.

Output: `soilstate.soilbeta_col`

Ported from `calc_beta_leepielke1992` in `SurfaceResistanceMod.F90`.
"""
function calc_beta_leepielke1992!(col::ColumnData,
                                   lun::LandunitData,
                                   soilstate::SoilStateData,
                                   waterstatebulk::WaterStateBulkData,
                                   waterdiagbulk::WaterDiagnosticBulkData,
                                   mask_nolakec::AbstractVector{Bool},
                                   bounds::UnitRange{Int})
    nlevsno = varpar.nlevsno
    joff = nlevsno  # offset for combined snow+soil indexing

    soilbeta = soilstate.soilbeta_col
    T = eltype(soilbeta)
    _launch!(_beta_leepielke_kernel!, soilbeta, mask_nolakec, col.landunit,
             col.itype, lun.itype, col.dz, soilstate.watsat_col, soilstate.watfc_col,
             waterstatebulk.ws.h2osoi_ice_col, waterstatebulk.ws.h2osoi_liq_col,
             waterdiagbulk.frac_sno_col, waterdiagbulk.frac_h2osfc_col,
             joff, T(DENH2O), T(DENICE))

    return nothing
end

# --------------------------------------------------------------------------
# Kernel: Swenson & Lawrence (2014) dry-surface-layer resistance (per column,
# masked). Each column writes only dsl[c] and soilresis[c]; independent.
# d_max / frac_sat_soil_dsl_init are host-resolved module params, passed at the
# output eltype so no Float64 reaches a Float32-only backend.
# --------------------------------------------------------------------------
@kernel function _soilresist_sl14_kernel!(soilresis, @Const(mask), @Const(landunit),
                                          @Const(col_itype), @Const(lun_itype),
                                          dsl, @Const(dz), @Const(watsat),
                                          @Const(bsw), @Const(sucsat), @Const(t_soisno),
                                          @Const(h2osoi_ice), @Const(h2osoi_liq),
                                          joff::Int, d_max, frac_sat_soil_dsl_init,
                                          denh2o, denice)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(soilresis)
        l = landunit[c]
        if lun_itype[l] != ISTWET && lun_itype[l] != ISTICE
            if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
                vwc_liq = max(h2osoi_liq[c, 1 + joff], T(1.0e-6)) / (dz[c, 1 + joff] * denh2o)

                eff_por_top = max(T(0.01), watsat[c, 1] - min(watsat[c, 1],
                                  h2osoi_ice[c, 1 + joff] / (dz[c, 1 + joff] * denice)))

                aird = watsat[c, 1] * (sucsat[c, 1] / T(1.0e7))^(one(T) / bsw[c, 1])
                d0 = T(2.12e-5) * (t_soisno[c, 1 + joff] / T(273.15))^T(1.75)
                eps = watsat[c, 1] - aird
                dg = eps * d0 * (eps / watsat[c, 1])^(T(3.0) / max(T(3.0), bsw[c, 1]))

                dsl[c] = d_max * max(T(0.001), (frac_sat_soil_dsl_init * eff_por_top - vwc_liq)) /
                         max(T(0.001), (frac_sat_soil_dsl_init * watsat[c, 1] - aird))

                dsl[c] = max(dsl[c], zero(T))
                dsl[c] = min(dsl[c], T(200.0))

                soilresis[c] = dsl[c] / (dg * eps * T(1.0e3)) + T(20.0)
                soilresis[c] = min(T(1.0e6), soilresis[c])

            elseif col_itype[c] == ICOL_ROAD_PERV
                soilresis[c] = T(1.0e6)
            elseif col_itype[c] == ICOL_SUNWALL || col_itype[c] == ICOL_SHADEWALL
                soilresis[c] = T(1.0e6)
            elseif col_itype[c] == ICOL_ROOF || col_itype[c] == ICOL_ROAD_IMPERV
                soilresis[c] = T(1.0e6)
            end
        else
            soilresis[c] = zero(T)
        end
    end
end

"""
    calc_soil_resistance_sl14!(col, lun, soilstate, waterstatebulk,
                                temperature, mask_nolakec, bounds)

Compute soil evaporative resistance using the Swenson & Lawrence (2014)
dry surface layer method.

Output: `soilstate.dsl_col`, `soilstate.soilresis_col`

Ported from `calc_soil_resistance_sl14` in `SurfaceResistanceMod.F90`.
"""
function calc_soil_resistance_sl14!(col::ColumnData,
                                     lun::LandunitData,
                                     soilstate::SoilStateData,
                                     waterstatebulk::WaterStateBulkData,
                                     temperature::TemperatureData,
                                     mask_nolakec::AbstractVector{Bool},
                                     bounds::UnitRange{Int})
    nlevsno = varpar.nlevsno
    joff = nlevsno  # offset for combined snow+soil indexing

    soilresis = soilstate.soilresis_col
    T = eltype(soilresis)

    d_max = surface_resistance_params.d_max
    frac_sat_soil_dsl_init = surface_resistance_params.frac_sat_soil_dsl_init

    _launch!(_soilresist_sl14_kernel!, soilresis, mask_nolakec, col.landunit,
             col.itype, lun.itype, soilstate.dsl_col, col.dz, soilstate.watsat_col,
             soilstate.bsw_col, soilstate.sucsat_col, temperature.t_soisno_col,
             waterstatebulk.ws.h2osoi_ice_col, waterstatebulk.ws.h2osoi_liq_col,
             joff, T(d_max), T(frac_sat_soil_dsl_init), T(DENH2O), T(DENICE))

    return nothing
end
