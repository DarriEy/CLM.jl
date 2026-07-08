# ==========================================================================
# Ported from: src/biogeophys/SurfaceHumidityMod.F90 (241 lines)
# Calculate surface humidities and intermediate variables needed for
# humidity calculations
# ==========================================================================

# --------------------------------------------------------------------------
# Device-view bundles for calculate_surface_humidity! (per-column kernel with
# an internal sequential nlevgrnd loop for the pervious-road branch). Bundling
# keeps the kernel under Metal's ~31-arg cap. Field names alias the original
# Julia variable paths so the kernel body reads like the scalar loop; the host
# path builds the same bundles from CPU arrays (so CPU stays byte-identical).
# --------------------------------------------------------------------------
Base.@kwdef struct _SHSoilDV{M, V}
    smpmin::V
    sucsat::M
    watsat::M
    watdry::M
    watopt::M
    bsw::M
    rootfr_road_perv::M
    rootr_road_perv::M
    soilalpha::V
    soilalpha_u::V
end
Adapt.@adapt_structure _SHSoilDV

Base.@kwdef struct _SHWaterDV{M, V}
    frac_h2osfc::V
    frac_sno_eff::V
    h2osoi_ice::M
    h2osoi_liq::M
    qg_snow::V
    qg_soil::V
    qg_h2osfc::V
    dqgdT::V
end
Adapt.@adapt_structure _SHWaterDV

Base.@kwdef struct _SHTempDV{M, V}
    t_h2osfc::V
    t_soisno::M
    t_grnd::V
end
Adapt.@adapt_structure _SHTempDV

# --------------------------------------------------------------------------
# Kernel: surface humidity per column (masked). Each column writes only its
# own indices; the internal nlevgrnd loop (pervious road) is sequential within
# the thread (loop-carried hr_road_perv), then normalizes. Literals are
# eltype-converted so no Float64 reaches a Float32-only backend; the qsat()
# helper is already eltype-generic and device-callable.
# --------------------------------------------------------------------------
@kernel function _surface_humidity_kernel!(qg, @Const(mask), @Const(landunit),
                                           @Const(col_itype), @Const(lun_itype),
                                           @Const(snl), @Const(dz),
                                           ss::_SHSoilDV, wd::_SHWaterDV, td::_SHTempDV,
                                           @Const(forc_pbot), @Const(forc_q),
                                           nlevgrnd::Int, joff::Int,
                                           denh2o, denice, roverg, tfrz, spval)
    c = @index(Global)
    @inbounds if mask[c]
        T = eltype(qg)
        l = landunit[c]

        hr_road_perv = zero(T)
        qred = one(T)
        hr = zero(T)

        if lun_itype[l] != ISTWET && lun_itype[l] != ISTICE
            if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
                wx  = (wd.h2osoi_liq[c, 1 + joff] / denh2o + wd.h2osoi_ice[c, 1 + joff] / denice) / dz[c, 1 + joff]
                fac = min(one(T), wx / ss.watsat[c, 1])
                fac = max(fac, T(0.01))
                psit = -ss.sucsat[c, 1] * fac^(-ss.bsw[c, 1])
                psit = max(ss.smpmin[c], psit)
                hr   = exp(psit / roverg / td.t_soisno[c, 1 + joff])
                qred = (one(T) - wd.frac_sno_eff[c] - wd.frac_h2osfc[c]) * hr +
                       wd.frac_sno_eff[c] + wd.frac_h2osfc[c]
                ss.soilalpha[c] = qred

            elseif col_itype[c] == ICOL_ROAD_PERV
                for j in 1:nlevgrnd
                    if td.t_soisno[c, j + joff] >= tfrz
                        vol_ice = min(ss.watsat[c, j], wd.h2osoi_ice[c, j + joff] / (dz[c, j + joff] * denice))
                        eff_porosity = ss.watsat[c, j] - vol_ice
                        vol_liq = min(eff_porosity, wd.h2osoi_liq[c, j + joff] / (dz[c, j + joff] * denh2o))
                        fac = min(max(vol_liq - ss.watdry[c, j], zero(T)) / (ss.watopt[c, j] - ss.watdry[c, j]), one(T))
                    else
                        fac = zero(T)
                    end
                    ss.rootr_road_perv[c, j] = ss.rootfr_road_perv[c, j] * fac
                    hr_road_perv = hr_road_perv + ss.rootr_road_perv[c, j]
                end
                qred = (one(T) - wd.frac_sno_eff[c]) * hr_road_perv + wd.frac_sno_eff[c]

                if hr_road_perv > zero(T)
                    for j in 1:nlevgrnd
                        ss.rootr_road_perv[c, j] = ss.rootr_road_perv[c, j] / hr_road_perv
                    end
                end
                ss.soilalpha_u[c] = qred

            elseif col_itype[c] == ICOL_SUNWALL || col_itype[c] == ICOL_SHADEWALL
                qred = zero(T)
                ss.soilalpha_u[c] = spval

            elseif col_itype[c] == ICOL_ROOF || col_itype[c] == ICOL_ROAD_IMPERV
                qred = one(T)
                ss.soilalpha_u[c] = spval
            end
        else
            ss.soilalpha[c] = spval
        end

        if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
            (qsatg, _, qsatgdT_soil, _) = qsat(td.t_soisno[c, 1 + joff], forc_pbot[c])
            if qsatg > forc_q[c] && forc_q[c] > hr * qsatg
                qsatg = forc_q[c]
                qsatgdT_soil = zero(T)
            end
            wd.qg_soil[c] = hr * qsatg

            if snl[c] < 0
                (qsatg_snow, _, qsatgdT_snow, _) = qsat(td.t_soisno[c, snl[c] + 1 + joff], forc_pbot[c])
                wd.qg_snow[c] = qsatg_snow
                wd.dqgdT[c] = wd.frac_sno_eff[c] * qsatgdT_snow +
                            (one(T) - wd.frac_sno_eff[c] - wd.frac_h2osfc[c]) * hr * qsatgdT_soil
            else
                wd.qg_snow[c] = wd.qg_soil[c]
                wd.dqgdT[c] = (one(T) - wd.frac_h2osfc[c]) * hr * qsatgdT_soil
            end

            if wd.frac_h2osfc[c] > zero(T)
                (qsatg_h2osfc, _, qsatgdT_h2osfc, _) = qsat(td.t_h2osfc[c], forc_pbot[c])
                wd.qg_h2osfc[c] = qsatg_h2osfc
                wd.dqgdT[c] = wd.dqgdT[c] + wd.frac_h2osfc[c] * qsatgdT_h2osfc
            else
                wd.qg_h2osfc[c] = wd.qg_soil[c]
            end

            qg[c] = wd.frac_sno_eff[c] * wd.qg_snow[c] +
                     (one(T) - wd.frac_sno_eff[c] - wd.frac_h2osfc[c]) * wd.qg_soil[c] +
                     wd.frac_h2osfc[c] * wd.qg_h2osfc[c]

        else
            (qsatg, _, qsatgdT, _) = qsat(td.t_grnd[c], forc_pbot[c])
            qg[c] = qred * qsatg
            wd.dqgdT[c] = qred * qsatgdT

            if qsatg > forc_q[c] && forc_q[c] > qred * qsatg
                qg[c] = forc_q[c]
                wd.dqgdT[c] = zero(T)
            end

            wd.qg_snow[c] = qg[c]
            wd.qg_soil[c] = qg[c]
            wd.qg_h2osfc[c] = qg[c]
        end
    end
end

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
- `mask_nolakec::AbstractVector{Bool}`: non-lake column mask
- `bounds::UnitRange{Int}`: column bounds

Ported from `CalculateSurfaceHumidity` in `SurfaceHumidityMod.F90`.
"""
function calculate_surface_humidity!(col::ColumnData,
                                     lun::LandunitData,
                                     temperature::TemperatureData,
                                     soilstate::SoilStateData,
                                     waterstatebulk::WaterStateBulkData,
                                     waterdiagbulk::WaterDiagnosticBulkData,
                                     forc_pbot::AbstractVector{<:Real},
                                     forc_q::AbstractVector{<:Real},
                                     mask_nolakec::AbstractVector{Bool},
                                     bounds::UnitRange{Int})

    nlevgrnd = varpar.nlevgrnd
    nlevsno  = varpar.nlevsno
    # Offset for converting Fortran soil-layer index to Julia combined snow+soil index
    joff = nlevsno

    qg = waterdiagbulk.qg_col
    T = eltype(qg)

    ss = _SHSoilDV(; smpmin = soilstate.smpmin_col,
                     sucsat = soilstate.sucsat_col,
                     watsat = soilstate.watsat_col,
                     watdry = soilstate.watdry_col,
                     watopt = soilstate.watopt_col,
                     bsw    = soilstate.bsw_col,
                     rootfr_road_perv = soilstate.rootfr_road_perv_col,
                     rootr_road_perv  = soilstate.rootr_road_perv_col,
                     soilalpha   = soilstate.soilalpha_col,
                     soilalpha_u = soilstate.soilalpha_u_col)
    wd = _SHWaterDV(; frac_h2osfc  = waterdiagbulk.frac_h2osfc_col,
                      frac_sno_eff = waterdiagbulk.frac_sno_eff_col,
                      h2osoi_ice   = waterstatebulk.ws.h2osoi_ice_col,
                      h2osoi_liq   = waterstatebulk.ws.h2osoi_liq_col,
                      qg_snow      = waterdiagbulk.qg_snow_col,
                      qg_soil      = waterdiagbulk.qg_soil_col,
                      qg_h2osfc    = waterdiagbulk.qg_h2osfc_col,
                      dqgdT        = waterdiagbulk.dqgdT_col)
    td = _SHTempDV(; t_h2osfc = temperature.t_h2osfc_col,
                     t_soisno = temperature.t_soisno_col,
                     t_grnd   = temperature.t_grnd_col)

    _launch!(_surface_humidity_kernel!, qg, mask_nolakec, col.landunit,
             col.itype, lun.itype, col.snl, col.dz, ss, wd, td,
             forc_pbot, forc_q, nlevgrnd, joff,
             T(DENH2O), T(DENICE), T(ROVERG), T(TFRZ), T(SPVAL))

    return nothing
end
