# ==========================================================================
# Ported from: src/biogeophys/SurfaceWaterMod.F90
# Surface water (h2osfc) dynamics — submerged fraction and state updates
#
# Determines fraction of land surfaces submerged based on microtopography
# and calculates water fluxes into and out of the surface water store.
#
# Public functions:
#   update_frac_h2osfc!  — Main entry: determine submerged fraction
#   update_h2osfc!       — Calculate fluxes out of h2osfc and update state
# ==========================================================================

# Module-level parameters (read from params file in Fortran)
const SURFACE_WATER_PC = 0.5    # threshold probability for surface water connectivity
const SURFACE_WATER_MU = 1.0    # connectivity exponent
const MIN_H2OSFC = 1.0e-8       # minimum surface water for numerical stability [mm]

"""
    bulkdiag_frac_h2osfc!(mask_soilc, bounds_col, dtime,
                           micro_sigma, h2osno_total,
                           h2osfc, frac_sno, frac_sno_eff,
                           frac_h2osfc, frac_h2osfc_nosnow,
                           qflx_too_small_h2osfc_to_soil)

Core algorithm to determine submerged fraction using microtopography
and surface water storage via Newton-Raphson iteration.

Ported from `BulkDiag_FracH2oSfc` in `SurfaceWaterMod.F90`.
"""
function bulkdiag_frac_h2osfc!(mask_soilc::BitVector,
                                bounds_col::UnitRange{Int},
                                dtime::Float64,
                                micro_sigma::Vector{Float64},
                                h2osno_total::Vector{Float64},
                                h2osfc::Vector{Float64},
                                frac_sno::Vector{Float64},
                                frac_sno_eff::Vector{Float64},
                                frac_h2osfc::Vector{Float64},
                                frac_h2osfc_nosnow::Vector{Float64},
                                qflx_too_small_h2osfc_to_soil::Vector{Float64})

    for c in bounds_col
        mask_soilc[c] || continue

        if h2osfc[c] > MIN_H2OSFC
            # Convert microtopography sigma from m to mm
            sigma = 1.0e3 * micro_sigma[c]
            sigma = max(sigma, 1.0e-6)  # prevent division by zero

            # Newton-Raphson to find water depth d such that
            # integral of PDF gives h2osfc
            d = 0.0
            for _iter in 1:10
                # f(d) = 0.5*d*(1+erf(d/(sigma*sqrt(2)))) + sigma/sqrt(2pi)*exp(-d^2/(2*sigma^2)) - h2osfc
                arg = d / (sigma * sqrt(2.0))
                fd = 0.5 * d * (1.0 + erf(arg)) +
                     sigma / sqrt(2.0 * π) * exp(-d^2 / (2.0 * sigma^2)) -
                     h2osfc[c]
                dfdd = 0.5 * (1.0 + erf(arg))
                dfdd = max(dfdd, 1.0e-12)
                d = d - fd / dfdd
            end

            # Compute submerged fraction
            frac_h2osfc[c] = 0.5 * (1.0 + erf(d / (sigma * sqrt(2.0))))
            frac_h2osfc_nosnow[c] = frac_h2osfc[c]
            qflx_too_small_h2osfc_to_soil[c] = 0.0
        else
            # Surface water too small — transfer to soil
            frac_h2osfc[c] = 0.0
            frac_h2osfc_nosnow[c] = 0.0
            qflx_too_small_h2osfc_to_soil[c] = h2osfc[c] / dtime
        end

        # Adjust fractions so frac_sno + frac_h2osfc <= 1
        if frac_sno[c] > (1.0 - frac_h2osfc[c]) && h2osno_total[c] > 0.0
            if frac_h2osfc[c] > 0.01
                frac_h2osfc[c] = max(1.0 - frac_sno[c], 0.01)
                frac_sno[c] = 1.0 - frac_h2osfc[c]
            else
                frac_sno[c] = 1.0 - frac_h2osfc[c]
            end
            frac_sno_eff[c] = frac_sno[c]
        end
    end

    return nothing
end

"""
    update_state_too_small_h2osfc_to_soil!(mask_soilc, bounds_col, dtime,
                                            qflx_too_small_h2osfc_to_soil,
                                            h2osfc, h2osoi_liq)

Transfer surface water below threshold to top soil layer.

Ported from `UpdateState_TooSmallH2osfcToSoil` in `SurfaceWaterMod.F90`.
"""
function update_state_too_small_h2osfc_to_soil!(mask_soilc::BitVector,
                                                 bounds_col::UnitRange{Int},
                                                 dtime::Float64,
                                                 qflx_too_small_h2osfc_to_soil::Vector{Float64},
                                                 h2osfc::Vector{Float64},
                                                 h2osoi_liq::Matrix{Float64},
                                                 nlevsno::Int)

    for c in bounds_col
        mask_soilc[c] || continue
        flux = qflx_too_small_h2osfc_to_soil[c]
        if flux > 0.0
            h2osfc[c] = h2osfc[c] - flux * dtime
            h2osoi_liq[c, 1 + nlevsno] = h2osoi_liq[c, 1 + nlevsno] + flux * dtime
            # Truncate small values
            if abs(h2osfc[c]) < 1.0e-10
                h2osfc[c] = 0.0
            end
        end
    end

    return nothing
end

"""
    update_frac_h2osfc!(water, col_data, mask_soilc, bounds_col; dtime)

Main entry point to determine fraction of land surface submerged and
update surface water state.

1. Calculate total snow water
2. Compute submerged fraction via Newton-Raphson
3. Transfer small h2osfc amounts to soil

Ported from `UpdateFracH2oSfc` in `SurfaceWaterMod.F90`.
"""
function update_frac_h2osfc!(water::WaterData,
                              col_data::ColumnData,
                              mask_soilc::BitVector,
                              bounds_col::UnitRange{Int};
                              dtime::Float64 = 1800.0)

    nlevsno_val = varpar.nlevsno
    wsb = water.waterstatebulk_inst
    wdb = water.waterdiagnosticbulk_inst
    wfb = water.waterfluxbulk_inst

    # Calculate total snow water equivalent
    nc_size = length(wsb.ws.h2osfc_col)
    h2osno_total = zeros(nc_size)
    qflx_too_small = zeros(nc_size)  # local array (not stored in WaterFluxBulkData)
    for c in bounds_col
        mask_soilc[c] || continue
        h2osno_total[c] = wsb.ws.h2osno_no_layers_col[c]
    end

    # Compute submerged fractions
    bulkdiag_frac_h2osfc!(mask_soilc, bounds_col, dtime,
                          col_data.micro_sigma,
                          h2osno_total,
                          wsb.ws.h2osfc_col,
                          wdb.frac_sno_col,
                          wdb.frac_sno_eff_col,
                          wdb.frac_h2osfc_col,
                          wdb.frac_h2osfc_nosnow_col,
                          qflx_too_small)

    # Transfer small surface water to soil
    update_state_too_small_h2osfc_to_soil!(mask_soilc, bounds_col, dtime,
                                           qflx_too_small,
                                           wsb.ws.h2osfc_col,
                                           wsb.ws.h2osoi_liq_col,
                                           nlevsno_val)

    return nothing
end

"""
    qflx_h2osfc_surf!(mask_hydrologyc, bounds_col, dtime,
                       h2osfcflag, h2osfc, h2osfc_thresh,
                       frac_h2osfc_nosnow, topo_slope,
                       qflx_h2osfc_surf)

Compute surface runoff from h2osfc using connectivity-dependent
linear reservoir model.

Ported from `QflxH2osfcSurf` in `SurfaceWaterMod.F90`.
"""
function qflx_h2osfc_surf!(mask_hydrologyc::BitVector,
                            bounds_col::UnitRange{Int},
                            dtime::Float64,
                            h2osfcflag::Int,
                            h2osfc::Vector{Float64},
                            h2osfc_thresh::Vector{Float64},
                            frac_h2osfc_nosnow::Vector{Float64},
                            topo_slope::Vector{Float64},
                            qflx_h2osfc_surf::Vector{Float64})

    for c in bounds_col
        mask_hydrologyc[c] || continue

        if h2osfcflag != 1
            qflx_h2osfc_surf[c] = 0.0
            continue
        end

        # Fractional connectivity (power law)
        frac_nosnow = frac_h2osfc_nosnow[c]
        if frac_nosnow <= SURFACE_WATER_PC
            frac_infclust = 0.0
        else
            frac_infclust = (frac_nosnow - SURFACE_WATER_PC)^SURFACE_WATER_MU
        end

        # Compute runoff if above threshold
        if h2osfc[c] > h2osfc_thresh[c] && h2osfcflag != 0
            k_wet = 1.0e-4 * sin(π / 180.0 * topo_slope[c])
            k_wet = max(k_wet, 1.0e-7)
            qflx_h2osfc_surf[c] = k_wet * frac_infclust * (h2osfc[c] - h2osfc_thresh[c])
            # Limit to available water
            qflx_h2osfc_surf[c] = min(qflx_h2osfc_surf[c],
                                       (h2osfc[c] - h2osfc_thresh[c]) / dtime)
        else
            qflx_h2osfc_surf[c] = 0.0
        end

        # Cutoff small flows
        if qflx_h2osfc_surf[c] < 1.0e-8
            qflx_h2osfc_surf[c] = 0.0
        end
    end

    return nothing
end

"""
    qflx_h2osfc_drain!(mask_hydrologyc, bounds_col, dtime,
                        h2osfcflag, h2osfc, frac_h2osfc,
                        qinmax, qflx_h2osfc_drain)

Compute infiltration/drainage from h2osfc into soil.

Ported from `QflxH2osfcDrain` in `SurfaceWaterMod.F90`.
"""
function qflx_h2osfc_drain!(mask_hydrologyc::BitVector,
                              bounds_col::UnitRange{Int},
                              dtime::Float64,
                              h2osfcflag::Int,
                              h2osfc::Vector{Float64},
                              frac_h2osfc::Vector{Float64},
                              qinmax::Vector{Float64},
                              qflx_h2osfc_drain::Vector{Float64})

    for c in bounds_col
        mask_hydrologyc[c] || continue

        if h2osfc[c] < 0.0
            # Numerical error recovery
            qflx_h2osfc_drain[c] = h2osfc[c] / dtime
        else
            qflx_h2osfc_drain[c] = min(frac_h2osfc[c] * qinmax[c],
                                        h2osfc[c] / dtime)
            if h2osfcflag == 0
                qflx_h2osfc_drain[c] = max(0.0, h2osfc[c] / dtime)
            end
        end
    end

    return nothing
end

"""
    update_h2osfc!(col_data, soilhydrology, energyflux,
                    waterfluxbulk, waterstatebulk, waterdiagbulk,
                    mask_hydrologyc, bounds_col;
                    dtime, h2osfcflag)

Calculate fluxes out of h2osfc and update the h2osfc state variable.

Ported from `UpdateH2osfc` in `SurfaceWaterMod.F90`.
"""
function update_h2osfc!(col_data::ColumnData,
                         soilhydrology::SoilHydrologyData,
                         energyflux::EnergyFluxData,
                         waterfluxbulk::WaterFluxBulkData,
                         waterstatebulk::WaterStateBulkData,
                         waterdiagbulk::WaterDiagnosticBulkData,
                         mask_hydrologyc::BitVector,
                         bounds_col::UnitRange{Int};
                         dtime::Float64 = 1800.0,
                         h2osfcflag::Int = 1)

    nc = length(waterstatebulk.ws.h2osfc_col)
    qflx_h2osfc_surf_arr = zeros(nc)
    qflx_h2osfc_drain_arr = zeros(nc)
    h2osfc_thresh = zeros(nc)  # percolation threshold
    topo_slope = zeros(nc)
    qinmax = zeros(nc)

    # Initialize arrays from column data
    for c in bounds_col
        mask_hydrologyc[c] || continue
        topo_slope[c] = col_data.topo_slope[c]
        h2osfc_thresh[c] = soilhydrology.h2osfc_thresh_col[c]
    end

    # 1. Compute surface runoff
    qflx_h2osfc_surf!(mask_hydrologyc, bounds_col, dtime,
                       h2osfcflag, waterstatebulk.ws.h2osfc_col,
                       h2osfc_thresh, waterdiagbulk.frac_h2osfc_nosnow_col,
                       topo_slope, qflx_h2osfc_surf_arr)

    # 2. Partially update h2osfc
    for c in bounds_col
        mask_hydrologyc[c] || continue
        waterstatebulk.ws.h2osfc_col[c] += (waterfluxbulk.qflx_in_h2osfc_col[c] -
                                         qflx_h2osfc_surf_arr[c]) * dtime
        # Truncate small values
        if abs(waterstatebulk.ws.h2osfc_col[c]) < 1.0e-10
            waterstatebulk.ws.h2osfc_col[c] = 0.0
        end
    end

    # 3. Compute drainage
    qflx_h2osfc_drain!(mask_hydrologyc, bounds_col, dtime,
                        h2osfcflag, waterstatebulk.ws.h2osfc_col,
                        waterdiagbulk.frac_h2osfc_col,
                        qinmax, qflx_h2osfc_drain_arr)

    # 4. Final update
    for c in bounds_col
        mask_hydrologyc[c] || continue
        waterstatebulk.ws.h2osfc_col[c] -= qflx_h2osfc_drain_arr[c] * dtime
        if abs(waterstatebulk.ws.h2osfc_col[c]) < 1.0e-10
            waterstatebulk.ws.h2osfc_col[c] = 0.0
        end
        waterfluxbulk.qflx_h2osfc_surf_col[c] = qflx_h2osfc_surf_arr[c]
        waterfluxbulk.qflx_h2osfc_drain_col[c] = qflx_h2osfc_drain_arr[c]
    end

    return nothing
end
