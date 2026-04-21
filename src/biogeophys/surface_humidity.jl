# ==========================================================================
# Ported from: src/biogeophys/SurfaceHumidityMod.F90 (241 lines)
# Calculate surface humidities and intermediate variables needed for
# humidity calculations
# ==========================================================================

"""
    calculate_surface_humidity!(col, lun, temperature, soilstate,
        waterstatebulk, waterdiagbulk, forc_pbot, forc_q,
        mask_nolakec, bounds)

Calculate surface humidities, as well as a few intermediate variables that are
needed in the humidity calculations.

# Arguments
- `col::ColumnData`: column data (snl, dz, landunit, itype)
- `lun::LandunitData`: landunit data (itype)
- `temperature::TemperatureData`: temperature state (t_h2osfc_col, t_soisno_col, t_grnd_col)
- `soilstate::SoilStateData`: soil state (input/output: soilalpha_col, soilalpha_u_col, rootr_road_perv_col; input: smpmin_col, sucsat_col, watsat_col, watdry_col, watopt_col, bsw_col, rootfr_road_perv_col)
- `waterstatebulk::WaterStateBulkData`: water state (h2osoi_ice_col, h2osoi_liq_col via ws)
- `waterdiagbulk::WaterDiagnosticBulkData`: water diagnostics (input: frac_h2osfc_col, frac_sno_eff_col; output: qg_snow_col, qg_soil_col, qg_col, qg_h2osfc_col, dqgdT_col)
- `forc_pbot::Vector{<:Real}`: atmospheric pressure, downscaled to column (Pa)
- `forc_q::Vector{<:Real}`: atmospheric specific humidity, downscaled to column (kg/kg)
- `mask_nolakec::BitVector`: non-lake column mask
- `bounds::UnitRange{Int}`: column bounds

Ported from `CalculateSurfaceHumidity` in `SurfaceHumidityMod.F90`.
"""
function calculate_surface_humidity!(col::ColumnData,
                                     lun::LandunitData,
                                     temperature::TemperatureData,
                                     soilstate::SoilStateData,
                                     waterstatebulk::WaterStateBulkData,
                                     waterdiagbulk::WaterDiagnosticBulkData,
                                     forc_pbot::Vector{<:Real},
                                     forc_q::Vector{<:Real},
                                     mask_nolakec::BitVector,
                                     bounds::UnitRange{Int})

    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno

    # Shorthand references (matching Fortran associate block)
    snl              = col.snl
    dz               = col.dz

    frac_h2osfc      = waterdiagbulk.frac_h2osfc_col
    frac_sno_eff     = waterdiagbulk.frac_sno_eff_col
    h2osoi_ice       = waterstatebulk.ws.h2osoi_ice_col
    h2osoi_liq       = waterstatebulk.ws.h2osoi_liq_col
    qg_snow          = waterdiagbulk.qg_snow_col
    qg_soil          = waterdiagbulk.qg_soil_col
    qg               = waterdiagbulk.qg_col
    qg_h2osfc        = waterdiagbulk.qg_h2osfc_col
    dqgdT            = waterdiagbulk.dqgdT_col

    smpmin           = soilstate.smpmin_col
    sucsat           = soilstate.sucsat_col
    watsat           = soilstate.watsat_col
    watdry           = soilstate.watdry_col
    watopt           = soilstate.watopt_col
    bsw              = soilstate.bsw_col
    rootfr_road_perv = soilstate.rootfr_road_perv_col
    rootr_road_perv  = soilstate.rootr_road_perv_col
    soilalpha        = soilstate.soilalpha_col
    soilalpha_u      = soilstate.soilalpha_u_col

    t_h2osfc         = temperature.t_h2osfc_col
    t_soisno         = temperature.t_soisno_col
    t_grnd           = temperature.t_grnd_col

    # Offset for converting Fortran soil-layer index to Julia combined snow+soil index
    # Fortran j ∈ {1,...,nlevgrnd} → Julia j + nlevsno
    joff = nlevsno

    for c in bounds
        mask_nolakec[c] || continue
        l = col.landunit[c]

        if col.itype[c] == ICOL_ROAD_PERV
            hr_road_perv = 0.0
        end

        # Saturated vapor pressure, specific humidity and their derivatives
        # at ground surface
        qred = 1.0
        hr = 0.0  # initialize; set properly in soil/crop branch

        if lun.itype[l] != ISTWET && lun.itype[l] != ISTICE

            if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                wx  = (h2osoi_liq[c, 1 + joff] / DENH2O + h2osoi_ice[c, 1 + joff] / DENICE) / dz[c, 1 + joff]
                fac = min(1.0, wx / watsat[c, 1])
                fac = max(fac, 0.01)
                psit = -sucsat[c, 1] * fac^(-bsw[c, 1])
                psit = max(smpmin[c], psit)
                # modify qred to account for h2osfc
                hr   = exp(psit / ROVERG / t_soisno[c, 1 + joff])
                qred = (1.0 - frac_sno_eff[c] - frac_h2osfc[c]) * hr +
                       frac_sno_eff[c] + frac_h2osfc[c]
                soilalpha[c] = qred

            elseif col.itype[c] == ICOL_ROAD_PERV
                # Pervious road depends on water in total soil column
                for j in 1:nlevgrnd
                    if t_soisno[c, j + joff] >= TFRZ
                        vol_ice = min(watsat[c, j], h2osoi_ice[c, j + joff] / (dz[c, j + joff] * DENICE))
                        eff_porosity = watsat[c, j] - vol_ice
                        vol_liq = min(eff_porosity, h2osoi_liq[c, j + joff] / (dz[c, j + joff] * DENH2O))
                        fac = min(max(vol_liq - watdry[c, j], 0.0) / (watopt[c, j] - watdry[c, j]), 1.0)
                    else
                        fac = 0.0
                    end
                    rootr_road_perv[c, j] = rootfr_road_perv[c, j] * fac
                    hr_road_perv = hr_road_perv + rootr_road_perv[c, j]
                end
                # Allows for sublimation of snow or dew on snow
                qred = (1.0 - frac_sno_eff[c]) * hr_road_perv + frac_sno_eff[c]

                # Normalize root resistances to get layer contribution to total ET
                if hr_road_perv > 0.0
                    for j in 1:nlevgrnd
                        rootr_road_perv[c, j] = rootr_road_perv[c, j] / hr_road_perv
                    end
                end
                soilalpha_u[c] = qred

            elseif col.itype[c] == ICOL_SUNWALL || col.itype[c] == ICOL_SHADEWALL
                qred = 0.0
                soilalpha_u[c] = SPVAL

            elseif col.itype[c] == ICOL_ROOF || col.itype[c] == ICOL_ROAD_IMPERV
                qred = 1.0
                soilalpha_u[c] = SPVAL
            end

        else
            soilalpha[c] = SPVAL
        end

        # compute humidities individually for snow, soil, h2osfc for vegetated landunits
        if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP

            (qsatg, _, qsatgdT_soil, _) = qsat(t_soisno[c, 1 + joff], forc_pbot[c])
            if qsatg > forc_q[c] && forc_q[c] > hr * qsatg
                qsatg = forc_q[c]
                qsatgdT_soil = 0.0
            end
            qg_soil[c] = hr * qsatg

            if snl[c] < 0
                (qsatg_snow, _, qsatgdT_snow, _) = qsat(t_soisno[c, snl[c] + 1 + joff], forc_pbot[c])
                qg_snow[c] = qsatg_snow
                dqgdT[c] = frac_sno_eff[c] * qsatgdT_snow +
                            (1.0 - frac_sno_eff[c] - frac_h2osfc[c]) * hr * qsatgdT_soil
            else
                # To be consistent with hs_top values in SoilTemp, set qg_snow to qg_soil
                # for snl = 0 case. This ensures hs_top_snow will equal hs_top_soil.
                qg_snow[c] = qg_soil[c]
                dqgdT[c] = (1.0 - frac_h2osfc[c]) * hr * qsatgdT_soil
            end

            if frac_h2osfc[c] > 0.0
                (qsatg_h2osfc, _, qsatgdT_h2osfc, _) = qsat(t_h2osfc[c], forc_pbot[c])
                qg_h2osfc[c] = qsatg_h2osfc
                dqgdT[c] = dqgdT[c] + frac_h2osfc[c] * qsatgdT_h2osfc
            else
                qg_h2osfc[c] = qg_soil[c]
            end

            qg[c] = frac_sno_eff[c] * qg_snow[c] +
                     (1.0 - frac_sno_eff[c] - frac_h2osfc[c]) * qg_soil[c] +
                     frac_h2osfc[c] * qg_h2osfc[c]

        else
            (qsatg, _, qsatgdT, _) = qsat(t_grnd[c], forc_pbot[c])
            qg[c] = qred * qsatg
            dqgdT[c] = qred * qsatgdT

            if qsatg > forc_q[c] && forc_q[c] > qred * qsatg
                qg[c] = forc_q[c]
                dqgdT[c] = 0.0
            end

            qg_snow[c] = qg[c]
            qg_soil[c] = qg[c]
            qg_h2osfc[c] = qg[c]
        end

    end  # end of columns loop

    return nothing
end
