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

# Module-level parameters (read from the params file in Fortran — SurfaceWaterMod
# readParams reads `pc` and `mu`). These were previously hardcoded to 0.5/1.0, but
# clm5_params.nc ships pc=0.4, mu=0.13889. The wrong values raised the connectivity
# threshold (frac_h2osfc_nosnow>pc for any outflow) and linearized the connectivity
# (mu=1 vs the strongly concave 0.139), throttling surface-water runoff. Read them
# from the params file into Refs (dereferenced on the host and passed to the kernel
# as scalars so no Ref reaches a GPU kernel).
const SURFACE_WATER_PC = Ref(0.4)    # threshold probability for surface water connectivity (params `pc`)
const SURFACE_WATER_MU = Ref(0.13889)# connectivity exponent (params `mu`)
const MIN_H2OSFC = 1.0e-8            # minimum surface water for numerical stability [mm]

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

"""
    compute_h2osfc_thresh!(h2osfc_thresh, micro_sigma, mask, bounds_col;
                            h2osfcflag, pc)

Compute the "fill & spill" surface-water threshold h2osfc_thresh_col (mm) from
microtopography. Ported from `SoilHydrologyInitTimeConst` in
`SoilHydrologyInitTimeConstMod.F90` (lines 216-236): a 4-iteration Newton solve
for the depth `d` at which the submerged fraction equals `pc`, then integrating
the storage at that depth. h2osfc must exceed this threshold before surface-water
runoff (QflxH2osfcSurf) can occur. Cold-start left it at 0, so runoff triggered
with no threshold (part of the cold-site QOVER/QDRAI partition error).

This is a one-time time-constant init, so it runs on the host.
"""
function compute_h2osfc_thresh!(h2osfc_thresh::AbstractVector{<:Real},
                                 micro_sigma::AbstractVector{<:Real},
                                 mask::AbstractVector{Bool},
                                 bounds_col::UnitRange{Int};
                                 h2osfcflag::Int = 1,
                                 pc::Real = SURFACE_WATER_PC[])
    sqrt2   = sqrt(2.0)
    sqrt2pi = sqrt(2.0 * π)
    for c in bounds_col
        mask[c] || continue
        h2osfc_thresh[c] = 0.0
        if micro_sigma[c] > 1.0e-6 && h2osfcflag != 0
            d = 0.0
            for _ in 1:4
                fd   = 0.5 * (1.0 + erf(d / (micro_sigma[c] * sqrt2))) - pc
                dfdd = exp(-d^2 / (2.0 * micro_sigma[c]^2)) / (micro_sigma[c] * sqrt2pi)
                d    = d - fd / dfdd
            end
            thr = 0.5 * d * (1.0 + erf(d / (micro_sigma[c] * sqrt2))) +
                  micro_sigma[c] / sqrt2pi * exp(-d^2 / (2.0 * micro_sigma[c]^2))
            h2osfc_thresh[c] = 1.0e3 * thr   # m -> mm
        end
    end
    return nothing
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
            # NumericsMod.truncate_small_values uses the pre-update storage as
            # the scale: only residuals below 1e-13 * |baseline| are zeroed.
            if abs(hnew) < T(1.0e-13) * abs(h2osfc[c])
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

    # Fortran BulkDiag_FracH2oSfc gates its snow-cover clamp on h2osno_total > 0,
    # where h2osno_total is the TOTAL snow water (CalculateTotalH2osno: unresolved
    # h2osno_no_layers + all explicit snow-layer ice+liq). Passing only
    # h2osno_no_layers_col was wrong: once a layered snowpack forms it is ~0, so the
    # clamp silently skipped the entire melt season and frac_h2osfc stayed at the
    # unclamped micro_sigma value → cold-site FH2OSFC/QINFL over-prediction
    # (Krycklan/Abisko/Iceland). Compute the true total to match Fortran.
    h2osno_total = fill!(similar(wsb.ws.h2osfc_col, FT, nc_size), zero(FT))  # device-resident scratch
    waterstate_calculate_total_h2osno!(wsb.ws, mask_soilc, bounds_col,
                                       col_data.snl, h2osno_total)

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
# --------------------------------------------------------------------------
# Kernel: surface runoff from h2osfc (connectivity-dependent linear reservoir).
# One thread per column; writes qflx_h2osfc_surf[c] only. Fully independent.
# `continue` early-outs are rewritten as per-element branches (single element).
# --------------------------------------------------------------------------
@kernel function _surfwat_qflx_surf_kernel!(qflx_h2osfc_surf, @Const(mask_hydrologyc),
                                            @Const(h2osfc), @Const(h2osfc_thresh),
                                            @Const(frac_h2osfc_nosnow), @Const(topo_slope),
                                            dtime, h2osfcflag::Int, pc, mu)
    c = @index(Global)
    T = eltype(qflx_h2osfc_surf)
    @inbounds if mask_hydrologyc[c]
        if h2osfcflag != 1
            qflx_h2osfc_surf[c] = zero(T)
        else
            # Fractional connectivity (power law); pc/mu are params `pc`,`mu`.
            frac_nosnow = frac_h2osfc_nosnow[c]
            if frac_nosnow <= T(pc)
                frac_infclust = zero(T)
            else
                frac_infclust = (frac_nosnow - T(pc))^T(mu)
            end

            # Compute runoff if above threshold. Fortran SurfaceWaterMod:486 uses
            # k_wet = 1e-4*sin(topo_slope) with NO floor for non-hillslope columns
            # (the min-slope guard is only inside the is_hillslope branch), so no
            # max(k_wet, 1e-7) here.
            q = zero(T)
            if h2osfc[c] > h2osfc_thresh[c] && h2osfcflag != 0
                k_wet = T(1.0e-4) * sin(T(π / 180.0) * topo_slope[c])
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
             h2osfc_thresh, frac_h2osfc_nosnow, topo_slope, dtime, h2osfcflag,
             SURFACE_WATER_PC[], SURFACE_WATER_MU[])

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
        # Match NumericsMod.truncate_small_values with the pre-update h2osfc
        # as data_baseline (SurfaceWaterMod.F90).
        if abs(h) < T(1.0e-13) * abs(h2osfc[c])
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
        # Here the partially updated h2osfc is the Fortran baseline.
        if abs(h) < T(1.0e-13) * abs(h2osfc[c])
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
                         h2osfcflag::Int = 1,
                         qinmax::Union{AbstractVector{<:Real},Nothing} = nothing)

    nc = length(waterstatebulk.ws.h2osfc_col)
    # Avoid capturing a Type-valued `FT` in the closure: on Julia 1.12 that boxes
    # FT and emits a runtime `_typeof_captured_variable`/`has_free_typevars` check
    # that Enzyme reverse-AD cannot differentiate. `similar(arr, n)` keeps eltype +
    # backend (GPU-safe); `zero(eltype(arr))` is constant-folded (no Type capture).
    _h2sc(n) = fill!(similar(waterstatebulk.ws.h2osfc_col, n), zero(eltype(waterstatebulk.ws.h2osfc_col)))  # device-resident
    qflx_h2osfc_surf_arr = _h2sc(nc)
    qflx_h2osfc_drain_arr = _h2sc(nc)
    # Fortran UpdateH2osfc (SurfaceWaterMod.F90:371) reads the maximum infiltration
    # rate `qinmax` from infiltration_excess_runoff_inst%qinmax_col (populated by the
    # InfiltrationExcessRunoff call earlier this step). Passing a zero array here makes
    # qflx_h2osfc_drain = min(frac_h2osfc*0, …) = 0, so ponded surface water never
    # infiltrates: infiltration (QINFL) collapses and the water spills to surface
    # runoff (QOVER) instead — the cold/frozen-soil partition defect. Use the real
    # qinmax when provided; fall back to zeros only for legacy callers/tests.
    qinmax = qinmax === nothing ? _h2sc(nc) : qinmax
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
