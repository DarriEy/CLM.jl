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
function surface_resistance_read_params!(; d_max::Float64 = 15.0,
                                           frac_sat_soil_dsl_init::Float64 = 0.8)
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
- `mask_nolakec::BitVector`: non-lake column mask
- `bounds::UnitRange{Int}`: column bounds

Ported from `calc_soilevap_resis` in `SurfaceResistanceMod.F90`.
"""
function calc_soilevap_resis!(col::ColumnData,
                               lun::LandunitData,
                               soilstate::SoilStateData,
                               waterstatebulk::WaterStateBulkData,
                               waterdiagbulk::WaterDiagnosticBulkData,
                               temperature::TemperatureData,
                               mask_nolakec::BitVector,
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
                                   mask_nolakec::BitVector,
                                   bounds::UnitRange{Int})
    nlevsno = varpar.nlevsno
    joff = nlevsno  # offset for combined snow+soil indexing

    # Shorthand references (matching Fortran associate block)
    watsat      = soilstate.watsat_col
    watfc       = soilstate.watfc_col
    soilbeta    = soilstate.soilbeta_col
    h2osoi_ice  = waterstatebulk.ws.h2osoi_ice_col
    h2osoi_liq  = waterstatebulk.ws.h2osoi_liq_col
    frac_sno    = waterdiagbulk.frac_sno_col
    frac_h2osfc = waterdiagbulk.frac_h2osfc_col

    for c in bounds
        mask_nolakec[c] || continue
        l = col.landunit[c]

        if lun.itype[l] != ISTWET && lun.itype[l] != ISTICE
            if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                wx  = (h2osoi_liq[c, 1 + joff] / DENH2O + h2osoi_ice[c, 1 + joff] / DENICE) / col.dz[c, 1 + joff]
                fac = min(1.0, wx / watsat[c, 1])
                fac = max(fac, 0.01)
                # Lee and Pielke 1992 beta, added by K.Sakaguchi
                if wx < watfc[c, 1]  # water content of top layer < field capacity
                    fac_fc = min(1.0, wx / watfc[c, 1])  # eqn5.66 but divided by theta at F.C.
                    fac_fc = max(fac_fc, 0.01)
                    # modify soil beta by snow cover. soilbeta for snow surface is one
                    soilbeta[c] = (1.0 - frac_sno[c] - frac_h2osfc[c]) *
                                  0.25 * (1.0 - cos(π * fac_fc))^2.0 +
                                  frac_sno[c] + frac_h2osfc[c]
                else  # water content of top layer >= field capacity
                    soilbeta[c] = 1.0
                end
            elseif col.itype[c] == ICOL_ROAD_PERV
                soilbeta[c] = 0.0
            elseif col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL
                soilbeta[c] = 0.0
            elseif col.itype[c] == ICOL_ROOF || col.itype[c] == ICOL_ROAD_IMPERV
                soilbeta[c] = 0.0
            end
        else
            soilbeta[c] = 1.0
        end
    end

    return nothing
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
                                     mask_nolakec::BitVector,
                                     bounds::UnitRange{Int})
    nlevsno = varpar.nlevsno
    joff = nlevsno  # offset for combined snow+soil indexing

    # Shorthand references (matching Fortran associate block)
    dz         = col.dz
    watsat     = soilstate.watsat_col
    bsw        = soilstate.bsw_col
    sucsat     = soilstate.sucsat_col
    dsl        = soilstate.dsl_col
    soilresis  = soilstate.soilresis_col
    t_soisno   = temperature.t_soisno_col
    h2osoi_ice = waterstatebulk.ws.h2osoi_ice_col
    h2osoi_liq = waterstatebulk.ws.h2osoi_liq_col

    d_max = surface_resistance_params.d_max
    frac_sat_soil_dsl_init = surface_resistance_params.frac_sat_soil_dsl_init

    for c in bounds
        mask_nolakec[c] || continue
        l = col.landunit[c]

        if lun.itype[l] != ISTWET && lun.itype[l] != ISTICE
            if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                vwc_liq = max(h2osoi_liq[c, 1 + joff], 1.0e-6) / (dz[c, 1 + joff] * DENH2O)

                # eff_porosity not calculated til SoilHydrology
                eff_por_top = max(0.01, watsat[c, 1] - min(watsat[c, 1],
                                  h2osoi_ice[c, 1 + joff] / (dz[c, 1 + joff] * DENICE)))

                # calculate diffusivity and air free pore space
                aird = watsat[c, 1] * (sucsat[c, 1] / 1.0e7)^(1.0 / bsw[c, 1])
                d0 = 2.12e-5 * (t_soisno[c, 1 + joff] / 273.15)^1.75  # [Bitelli et al., JH, 08]
                eps = watsat[c, 1] - aird
                dg = eps * d0 * (eps / watsat[c, 1])^(3.0 / max(3.0, bsw[c, 1]))

                dsl[c] = d_max * max(0.001, (frac_sat_soil_dsl_init * eff_por_top - vwc_liq)) /
                         max(0.001, (frac_sat_soil_dsl_init * watsat[c, 1] - aird))

                dsl[c] = max(dsl[c], 0.0)
                dsl[c] = min(dsl[c], 200.0)

                soilresis[c] = dsl[c] / (dg * eps * 1.0e3) + 20.0
                soilresis[c] = min(1.0e6, soilresis[c])

            elseif col.itype[c] == ICOL_ROAD_PERV
                soilresis[c] = 1.0e6
            elseif col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL
                soilresis[c] = 1.0e6
            elseif col.itype[c] == ICOL_ROOF || col.itype[c] == ICOL_ROAD_IMPERV
                soilresis[c] = 1.0e6
            end
        else
            soilresis[c] = 0.0
        end
    end

    return nothing
end
