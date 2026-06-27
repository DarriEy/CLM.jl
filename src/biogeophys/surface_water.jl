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
# GPU-safe erf: exact (SpecialFunctions) on Float64 / ForwardDiff.Dual (the generic
# fallback — keeps CPU byte-identical AND AD-differentiable); an Abramowitz-Stegun
# 7.1.26 polynomial on Float32 (max abs err ~1.5e-7) since Metal's erf pulls in Float64.
@inline _erf_gpu(x) = erf(x)
@inline function _erf_gpu(x::Float32)
    s = sign(x); a = abs(x)
    t = 1f0 / (1f0 + 0.3275911f0 * a)
    y = 1f0 - (((((1.061405429f0 * t - 1.453152027f0) * t) + 1.421413741f0) * t -
              0.284496736f0) * t + 0.254829592f0) * t * exp(-a * a)
    return s * y
end

# Per-column kernel: submerged-fraction Newton (erf/exp) + frac_sno/h2osfc cap.
# One thread per column, own-index writes only. Float64 consts/literals wrapped
# at the working precision; byte-identical on CPU (Float64).
@kernel function _surfwat_frac_h2osfc_kernel!(frac_h2osfc, frac_h2osfc_nosnow,
        qflx_too_small, frac_sno, frac_sno_eff, @Const(mask_soilc),
        @Const(micro_sigma), @Const(h2osno_total), @Const(h2osfc), dtime,
        cmin::Int, cmax::Int)
    c = @index(Global)
    T = eltype(frac_h2osfc)
    @inbounds if cmin <= c <= cmax && mask_soilc[c]
        if h2osfc[c] > T(MIN_H2OSFC)
            sigma = T(1.0e3) * micro_sigma[c]
            sigma = max(sigma, T(1.0e-6))
            d = zero(T)
            for _iter in 1:10
                arg = d / (sigma * sqrt(T(2.0)))
                fd = T(0.5) * d * (one(T) + _erf_gpu(arg)) +
                     sigma / sqrt(T(2.0) * T(π)) * exp(-d^2 / (T(2.0) * sigma^2)) - h2osfc[c]
                dfdd = max(T(0.5) * (one(T) + _erf_gpu(arg)), T(1.0e-12))
                d = d - fd / dfdd
            end
            fh = T(0.5) * (one(T) + _erf_gpu(d / (sigma * sqrt(T(2.0)))))
            frac_h2osfc[c] = fh
            frac_h2osfc_nosnow[c] = fh
            qflx_too_small[c] = zero(T)
        else
            frac_h2osfc[c] = zero(T)
            frac_h2osfc_nosnow[c] = zero(T)
            qflx_too_small[c] = h2osfc[c] / dtime
        end

        if frac_sno[c] > (one(T) - frac_h2osfc[c]) && h2osno_total[c] > zero(T)
            if frac_h2osfc[c] > T(0.01)
                frac_h2osfc[c] = max(one(T) - frac_sno[c], T(0.01))
                frac_sno[c] = one(T) - frac_h2osfc[c]
            else
                frac_sno[c] = one(T) - frac_h2osfc[c]
            end
            frac_sno_eff[c] = frac_sno[c]
        end
    end
end

function bulkdiag_frac_h2osfc!(mask_soilc::AbstractVector{Bool},
                                bounds_col::UnitRange{Int},
                                dtime::Real,
                                micro_sigma::AbstractVector{<:Real},
                                h2osno_total::AbstractVector{<:Real},
                                h2osfc::AbstractVector{<:Real},
                                frac_sno::AbstractVector{<:Real},
                                frac_sno_eff::AbstractVector{<:Real},
                                frac_h2osfc::AbstractVector{<:Real},
                                frac_h2osfc_nosnow::AbstractVector{<:Real},
                                qflx_too_small_h2osfc_to_soil::AbstractVector{<:Real})
    _launch!(_surfwat_frac_h2osfc_kernel!, frac_h2osfc, frac_h2osfc_nosnow,
             qflx_too_small_h2osfc_to_soil, frac_sno, frac_sno_eff, mask_soilc,
             micro_sigma, h2osno_total, h2osfc, eltype(frac_h2osfc)(dtime),
             first(bounds_col), last(bounds_col); ndrange = last(bounds_col))
    return nothing
end

"""
    update_state_too_small_h2osfc_to_soil!(mask_soilc, bounds_col, dtime,
                                            qflx_too_small_h2osfc_to_soil,
                                            h2osfc, h2osoi_liq)

Transfer surface water below threshold to top soil layer.

Ported from `UpdateState_TooSmallH2osfcToSoil` in `SurfaceWaterMod.F90`.
"""
# --------------------------------------------------------------------------
# Kernel: transfer below-threshold surface water to the top soil layer.
# One thread per column; writes h2osfc[c] and h2osoi_liq[c, 1+nlevsno].
# Fully independent across columns (no reduction/loop-carried dependence).
# --------------------------------------------------------------------------
@kernel function _surfwat_too_small_to_soil_kernel!(h2osfc, h2osoi_liq,
                                                    @Const(mask_soilc),
                                                    @Const(qflx_too_small_h2osfc_to_soil),
                                                    dtime, nlevsno::Int)
    c = @index(Global)
    T = eltype(h2osfc)
    @inbounds if mask_soilc[c]
        flux = qflx_too_small_h2osfc_to_soil[c]
        if flux > zero(T)
            hnew = h2osfc[c] - flux * dtime
            h2osoi_liq[c, 1 + nlevsno] = h2osoi_liq[c, 1 + nlevsno] + flux * dtime
            # Truncate small values
            if abs(hnew) < T(1.0e-10)
                hnew = zero(T)
            end
            h2osfc[c] = hnew
        end
    end
end

surfwat_too_small_to_soil!(h2osfc, h2osoi_liq, mask_soilc,
                           qflx_too_small_h2osfc_to_soil, dtime, nlevsno::Int) =
    _launch!(_surfwat_too_small_to_soil_kernel!, h2osfc, h2osoi_liq, mask_soilc,
             qflx_too_small_h2osfc_to_soil, dtime, nlevsno)

function update_state_too_small_h2osfc_to_soil!(mask_soilc::AbstractVector{Bool},
                                                 bounds_col::UnitRange{Int},
                                                 dtime::Real,
                                                 qflx_too_small_h2osfc_to_soil::AbstractVector{<:Real},
                                                 h2osfc::AbstractVector{<:Real},
                                                 h2osoi_liq::AbstractMatrix{<:Real},
                                                 nlevsno::Int)

    surfwat_too_small_to_soil!(h2osfc, h2osoi_liq, mask_soilc,
                               qflx_too_small_h2osfc_to_soil, dtime, nlevsno)

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
                              mask_soilc::AbstractVector{Bool},
                              bounds_col::UnitRange{Int};
                              dtime::Real = 1800.0)

    nlevsno_val = varpar.nlevsno
    wsb = water.waterstatebulk_inst
    wdb = water.waterdiagnosticbulk_inst
    wfb = water.waterfluxbulk_inst

    # Calculate total snow water equivalent.
    nc_size = length(wsb.ws.h2osfc_col)
    FT = eltype(wsb.ws.h2osfc_col)
    qflx_too_small = fill!(similar(wsb.ws.h2osfc_col, FT, nc_size), zero(FT))  # device-resident scratch
    # h2osno_total is only read (for masked cols) inside bulkdiag_frac_h2osfc!, where
    # it equals h2osno_no_layers_col — pass that directly (byte-identical, no scratch).

    # Compute submerged fractions
    bulkdiag_frac_h2osfc!(mask_soilc, bounds_col, dtime,
                          col_data.micro_sigma,
                          wsb.ws.h2osno_no_layers_col,
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
# --------------------------------------------------------------------------
# Kernel: surface runoff from h2osfc (connectivity-dependent linear reservoir).
# One thread per column; writes qflx_h2osfc_surf[c] only. Fully independent.
# `continue` early-outs are rewritten as per-element branches (single element).
# --------------------------------------------------------------------------
@kernel function _surfwat_qflx_surf_kernel!(qflx_h2osfc_surf, @Const(mask_hydrologyc),
                                            @Const(h2osfc), @Const(h2osfc_thresh),
                                            @Const(frac_h2osfc_nosnow), @Const(topo_slope),
                                            dtime, h2osfcflag::Int)
    c = @index(Global)
    T = eltype(qflx_h2osfc_surf)
    @inbounds if mask_hydrologyc[c]
        if h2osfcflag != 1
            qflx_h2osfc_surf[c] = zero(T)
        else
            # Fractional connectivity (power law)
            frac_nosnow = frac_h2osfc_nosnow[c]
            if frac_nosnow <= T(SURFACE_WATER_PC)
                frac_infclust = zero(T)
            else
                frac_infclust = (frac_nosnow - T(SURFACE_WATER_PC))^T(SURFACE_WATER_MU)
            end

            # Compute runoff if above threshold
            q = zero(T)
            if h2osfc[c] > h2osfc_thresh[c] && h2osfcflag != 0
                k_wet = T(1.0e-4) * sin(T(π / 180.0) * topo_slope[c])
                k_wet = max(k_wet, T(1.0e-7))
                q = k_wet * frac_infclust * (h2osfc[c] - h2osfc_thresh[c])
                # Limit to available water
                q = min(q, (h2osfc[c] - h2osfc_thresh[c]) / dtime)
            end

            # Cutoff small flows
            if q < T(1.0e-8)
                q = zero(T)
            end
            qflx_h2osfc_surf[c] = q
        end
    end
end

surfwat_qflx_surf!(qflx_h2osfc_surf, mask_hydrologyc, h2osfc, h2osfc_thresh,
                   frac_h2osfc_nosnow, topo_slope, dtime, h2osfcflag::Int) =
    _launch!(_surfwat_qflx_surf_kernel!, qflx_h2osfc_surf, mask_hydrologyc, h2osfc,
             h2osfc_thresh, frac_h2osfc_nosnow, topo_slope, dtime, h2osfcflag)

function qflx_h2osfc_surf!(mask_hydrologyc::AbstractVector{Bool},
                            bounds_col::UnitRange{Int},
                            dtime::Real,
                            h2osfcflag::Int,
                            h2osfc::AbstractVector{<:Real},
                            h2osfc_thresh::AbstractVector{<:Real},
                            frac_h2osfc_nosnow::AbstractVector{<:Real},
                            topo_slope::AbstractVector{<:Real},
                            qflx_h2osfc_surf::AbstractVector{<:Real})

    surfwat_qflx_surf!(qflx_h2osfc_surf, mask_hydrologyc, h2osfc, h2osfc_thresh,
                       frac_h2osfc_nosnow, topo_slope, dtime, h2osfcflag)

    return nothing
end

"""
    qflx_h2osfc_drain!(mask_hydrologyc, bounds_col, dtime,
                        h2osfcflag, h2osfc, frac_h2osfc,
                        qinmax, qflx_h2osfc_drain)

Compute infiltration/drainage from h2osfc into soil.

Ported from `QflxH2osfcDrain` in `SurfaceWaterMod.F90`.
"""
# --------------------------------------------------------------------------
# Kernel: infiltration/drainage from h2osfc into soil.
# One thread per column; writes qflx_h2osfc_drain[c] only. Fully independent.
# --------------------------------------------------------------------------
@kernel function _surfwat_qflx_drain_kernel!(qflx_h2osfc_drain, @Const(mask_hydrologyc),
                                             @Const(h2osfc), @Const(frac_h2osfc),
                                             @Const(qinmax), dtime, h2osfcflag::Int)
    c = @index(Global)
    T = eltype(qflx_h2osfc_drain)
    @inbounds if mask_hydrologyc[c]
        if h2osfc[c] < zero(T)
            # Numerical error recovery
            qflx_h2osfc_drain[c] = h2osfc[c] / dtime
        else
            q = min(frac_h2osfc[c] * qinmax[c], h2osfc[c] / dtime)
            if h2osfcflag == 0
                q = max(zero(T), h2osfc[c] / dtime)
            end
            qflx_h2osfc_drain[c] = q
        end
    end
end

surfwat_qflx_drain!(qflx_h2osfc_drain, mask_hydrologyc, h2osfc, frac_h2osfc,
                    qinmax, dtime, h2osfcflag::Int) =
    _launch!(_surfwat_qflx_drain_kernel!, qflx_h2osfc_drain, mask_hydrologyc, h2osfc,
             frac_h2osfc, qinmax, dtime, h2osfcflag)

function qflx_h2osfc_drain!(mask_hydrologyc::AbstractVector{Bool},
                              bounds_col::UnitRange{Int},
                              dtime::Real,
                              h2osfcflag::Int,
                              h2osfc::AbstractVector{<:Real},
                              frac_h2osfc::AbstractVector{<:Real},
                              qinmax::AbstractVector{<:Real},
                              qflx_h2osfc_drain::AbstractVector{<:Real})

    surfwat_qflx_drain!(qflx_h2osfc_drain, mask_hydrologyc, h2osfc, frac_h2osfc,
                        qinmax, dtime, h2osfcflag)

    return nothing
end

# --------------------------------------------------------------------------
# Kernel: step-2 partial h2osfc update (inflow minus surface runoff).
# One thread per column; writes h2osfc[c] only. Fully independent.
# --------------------------------------------------------------------------
@kernel function _surfwat_partial_update_kernel!(h2osfc, @Const(mask_hydrologyc),
                                                 @Const(qflx_in_h2osfc),
                                                 @Const(qflx_h2osfc_surf), dtime)
    c = @index(Global)
    T = eltype(h2osfc)
    @inbounds if mask_hydrologyc[c]
        h = h2osfc[c] + (qflx_in_h2osfc[c] - qflx_h2osfc_surf[c]) * dtime
        # Truncate small values
        if abs(h) < T(1.0e-10)
            h = zero(T)
        end
        h2osfc[c] = h
    end
end

surfwat_partial_update!(h2osfc, mask_hydrologyc, qflx_in_h2osfc,
                        qflx_h2osfc_surf, dtime) =
    _launch!(_surfwat_partial_update_kernel!, h2osfc, mask_hydrologyc,
             qflx_in_h2osfc, qflx_h2osfc_surf, dtime)

# --------------------------------------------------------------------------
# Kernel: step-4 final h2osfc update (subtract drainage) and store diagnostics.
# One thread per column; writes h2osfc[c], qflx_h2osfc_surf_col[c],
# qflx_h2osfc_drain_col[c]. Fully independent across columns.
# --------------------------------------------------------------------------
@kernel function _surfwat_final_update_kernel!(h2osfc, qflx_h2osfc_surf_col,
                                               qflx_h2osfc_drain_col,
                                               @Const(mask_hydrologyc),
                                               @Const(qflx_h2osfc_drain_arr),
                                               @Const(qflx_h2osfc_surf_arr), dtime)
    c = @index(Global)
    T = eltype(h2osfc)
    @inbounds if mask_hydrologyc[c]
        h = h2osfc[c] - qflx_h2osfc_drain_arr[c] * dtime
        if abs(h) < T(1.0e-10)
            h = zero(T)
        end
        h2osfc[c] = h
        qflx_h2osfc_surf_col[c] = qflx_h2osfc_surf_arr[c]
        qflx_h2osfc_drain_col[c] = qflx_h2osfc_drain_arr[c]
    end
end

surfwat_final_update!(h2osfc, qflx_h2osfc_surf_col, qflx_h2osfc_drain_col,
                      mask_hydrologyc, qflx_h2osfc_drain_arr,
                      qflx_h2osfc_surf_arr, dtime) =
    _launch!(_surfwat_final_update_kernel!, h2osfc, qflx_h2osfc_surf_col,
             qflx_h2osfc_drain_col, mask_hydrologyc, qflx_h2osfc_drain_arr,
             qflx_h2osfc_surf_arr, dtime)

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
                         mask_hydrologyc::AbstractVector{Bool},
                         bounds_col::UnitRange{Int};
                         dtime::Real = 1800.0,
                         h2osfcflag::Int = 1)

    nc = length(waterstatebulk.ws.h2osfc_col)
    # Avoid capturing a Type-valued `FT` in the closure: on Julia 1.12 that boxes
    # FT and emits a runtime `_typeof_captured_variable`/`has_free_typevars` check
    # that Enzyme reverse-AD cannot differentiate. `similar(arr, n)` keeps eltype +
    # backend (GPU-safe); `zero(eltype(arr))` is constant-folded (no Type capture).
    _h2sc(n) = fill!(similar(waterstatebulk.ws.h2osfc_col, n), zero(eltype(waterstatebulk.ws.h2osfc_col)))  # device-resident
    qflx_h2osfc_surf_arr = _h2sc(nc)
    qflx_h2osfc_drain_arr = _h2sc(nc)
    qinmax = _h2sc(nc)
    # topo_slope / h2osfc_thresh are read only for masked columns inside the kernels,
    # where they equal the column arrays — pass those directly (no host masked-copy).
    topo_slope    = col_data.topo_slope
    h2osfc_thresh = soilhydrology.h2osfc_thresh_col

    # 1. Compute surface runoff
    qflx_h2osfc_surf!(mask_hydrologyc, bounds_col, dtime,
                       h2osfcflag, waterstatebulk.ws.h2osfc_col,
                       h2osfc_thresh, waterdiagbulk.frac_h2osfc_nosnow_col,
                       topo_slope, qflx_h2osfc_surf_arr)

    # 2. Partially update h2osfc
    surfwat_partial_update!(waterstatebulk.ws.h2osfc_col, mask_hydrologyc,
                            waterfluxbulk.qflx_in_h2osfc_col, qflx_h2osfc_surf_arr,
                            dtime)

    # 3. Compute drainage
    qflx_h2osfc_drain!(mask_hydrologyc, bounds_col, dtime,
                        h2osfcflag, waterstatebulk.ws.h2osfc_col,
                        waterdiagbulk.frac_h2osfc_col,
                        qinmax, qflx_h2osfc_drain_arr)

    # 4. Final update
    surfwat_final_update!(waterstatebulk.ws.h2osfc_col,
                          waterfluxbulk.qflx_h2osfc_surf_col,
                          waterfluxbulk.qflx_h2osfc_drain_col,
                          mask_hydrologyc, qflx_h2osfc_drain_arr,
                          qflx_h2osfc_surf_arr, dtime)

    return nothing
end
