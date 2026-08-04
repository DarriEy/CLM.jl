# ==========================================================================
# Ported from: src/main/atm2lndMod.F90
# Atmospheric forcing downscaling from gridcell to column
#
# Public functions:
#   downscale_forcings!    — Main downscaling routine
#   partition_precip!      — Rain/snow partitioning from temperature
#   downscale_longwave!    — Elevation-based LW downscaling
#
# GPU: every per-column loop is a KernelAbstractions kernel launched via
# `_launch!`, so the routine runs on whichever backend its arrays live on —
# KA.CPU (Float64, byte-identical to the original loop) on the host path, and
# the GPU backend (Metal Float32, …) when handed device arrays. One thread per
# column; `lo = bounds.begc` maps thread i → column c = lo + i - 1. Constants
# are pre-converted to the working eltype (no Float64 reaches Metal). Kernels
# that read gridcell fields index them with `col_gridcell[c]`; the CO2/O2
# rescale writes gridcell fields from columns and so assumes the parity-domain
# 1:1 column↔gridcell mapping (matches the original loop's semantics).
#
# Skipped: hillslope downscaling, water tracers
# ==========================================================================

# --------------------------------------------------------------------------
# Kernels
# --------------------------------------------------------------------------

# Step 1: copy the not-downscaled gridcell forcing to every column as the
# baseline (later steps overwrite the elevation-corrected fields in place).
@kernel function _downscale_baseline_kernel!(t_ds, th_ds, pbot_ds, rho_ds, lw_ds,
        solad_ds, solar_ds, @Const(col_gridcell), @Const(t_nd), @Const(th_nd),
        @Const(pbot_nd), @Const(rho_nd), @Const(lw_nd), @Const(solad_nd),
        @Const(solar_nd), lo::Int, numrad::Int)
    i = @index(Global)
    @inbounds begin
        c = lo + i - 1
        g = col_gridcell[c]
        t_ds[c]    = t_nd[g]
        th_ds[c]   = th_nd[g]
        pbot_ds[c] = pbot_nd[g]
        rho_ds[c]  = rho_nd[g]
        lw_ds[c]   = lw_nd[g]
        for b in 1:numrad
            solad_ds[c, b] = solad_nd[g, b]
        end
        solar_ds[c] = solar_nd[g]
    end
end

# Step 2: elevation-based downscaling of temperature, potential temperature,
# pressure (hydrostatic balance) and density (ratio method).
@kernel function _downscale_temp_kernel!(t_ds, th_ds, pbot_ds, rho_ds,
        @Const(col_gridcell), @Const(topo_col), @Const(topo_grc),
        @Const(t_nd), @Const(th_nd), @Const(pbot_nd), @Const(rho_nd), @Const(hgt_grc),
        lapse_rate, rair, grav, cpair, lo::Int)
    i = @index(Global)
    @inbounds begin
        T = eltype(t_ds)
        c = lo + i - 1
        g = col_gridcell[c]
        hsurf_g = topo_grc[g]
        hsurf_c = topo_col[c]
        if !(isnan(hsurf_c) || isnan(hsurf_g))
            tbot_g = t_nd[g]
            tbot_c = tbot_g - lapse_rate * (hsurf_c - hsurf_g)
            t_ds[c] = tbot_c

            zbot = hgt_grc[g]
            Hbot = rair * T(0.5) * (tbot_g + tbot_c) / grav
            if Hbot > zero(T)
                # Potential temperature correction
                th_ds[c] = th_nd[g] + (tbot_c - tbot_g) * exp((zbot / Hbot) * (rair / cpair))
                # Pressure via hydrostatic balance
                pbot_ds[c] = pbot_nd[g] * exp(-(hsurf_c - hsurf_g) / Hbot)
            end

            # Density correction via ratio method
            pbot_g = pbot_nd[g]
            pbot_c = pbot_ds[c]
            rho_est_g = pbot_g / (rair * tbot_g)
            rho_est_c = pbot_c / (rair * tbot_c)
            if rho_est_g > zero(T)
                rho_ds[c] = rho_nd[g] * (rho_est_c / rho_est_g)
            end
        end
    end
end

# Step 2b, pass 1: accumulate the gridcell mean of the downscaled column pbot.
#
# The weight is the column's area fraction. Two cases need care, and they are
# NOT the same case:
#
#   * NON-FINITE weight -> weight 1. `column_init!` allocates `wtgcell` as NaN
#     and CLM.jl's hand-built driver fixtures routinely leave it unset; falling
#     back to an unweighted mean keeps them working, and for their
#     single-column-per-gridcell domains the mean IS that column's pbot, so the
#     per-column rescale this replaced is reproduced exactly.
#
#   * ZERO weight -> weight 0, i.e. it contributes NOTHING. A zero-area column
#     is CTSM's "virtual" column: a glacier elevation-class column that exists
#     only so a surface mass balance can be computed at that class elevation.
#     Lumping those in (an earlier version of this kernel treated zero like
#     unset and gave them weight 1) drags the gridcell mean down to the mean
#     CLASS elevation of the whole ladder — several km up — and drove the
#     implied pco2 pressure ~10 kPa below the surface pressure in a gridcell
#     that is not glaciated at all.
#
# If every weight in a gridcell is zero, wsum is 0 and the rescale is skipped.
@kernel function _downscale_pbot_accum_kernel!(psum_grc, wsum_grc, @Const(col_gridcell),
        @Const(col_wtgcell), @Const(pbot_ds), cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax
        T = eltype(psum_grc)
        g = col_gridcell[c]
        w = col_wtgcell[c]
        w = isfinite(w) ? max(T(w), zero(T)) : one(T)
        _scatter_add!(psum_grc, g, T(pbot_ds[c]) * w)
        _scatter_add!(wsum_grc, g, w)
    end
end

# Step 2b, pass 2: rescale CO2/O2 partial pressures from the not-downscaled
# surface pbot to the gridcell-mean elevation-corrected pbot (see note at the
# call site). Indexed over GRIDCELLS — one rescale per gridcell.
#
# It MUST be one: this is a read-modify-write of a gridcell field, so a
# column-indexed version applied the ratio once per column and compounded it.
# With a glacier elevation-class subgrid (11 columns per gridcell) that drove
# forc_pco2 to ~pbot*(pbot_ds/pbot_nd)^11 — a ~50% CO2 error, invisible in the
# pco2/po2 molar-ratio invariant because both compound identically.
@kernel function _downscale_pco2_kernel!(pco2_grc, po2_grc, pc13o2_grc,
        @Const(pbot_nd), @Const(psum_grc), @Const(wsum_grc), gmin::Int, gmax::Int)
    g = @index(Global)
    @inbounds if gmin <= g <= gmax
        T = eltype(pco2_grc)
        pbot_nd_g = pbot_nd[g]
        wsum = wsum_grc[g]
        pbot_ds_g = wsum > zero(T) ? psum_grc[g] / wsum : zero(T)
        if pbot_nd_g > zero(T) && pbot_ds_g > zero(T)
            r = pbot_ds_g / pbot_nd_g
            pco2_grc[g] *= r
            po2_grc[g]  *= r
            # forc_pc13o2 is C13RATIO*forc_pco2 (the reader seeds it): rescale it by
            # the same pbot ratio so rc13_canair = pc13o2/(pco2-pc13o2) stays exactly
            # C13RATIO/(1-C13RATIO) at the elevation-corrected column pbot.
            pc13o2_grc[g] *= r
        end
    end
end

# Rain/snow partitioning from column temperature (ramp) or the gridcell split.
@kernel function _partition_precip_kernel!(rain_ds, snow_ds, @Const(col_gridcell),
        @Const(col_landunit), @Const(lun_itype), @Const(rain_nd), @Const(snow_nd),
        @Const(t_ds), repart::Bool, glc_snow_t, glc_slope, nonglc_snow_t, nonglc_slope,
        istice::Int, nrain::Int, nsnow::Int, lo::Int)
    i = @index(Global)
    @inbounds begin
        T = eltype(rain_ds)
        c = lo + i - 1
        g = col_gridcell[c]
        total_rain = g <= nrain ? rain_nd[g] : zero(T)
        total_snow = g <= nsnow ? snow_nd[g] : zero(T)
        total_precip = total_rain + total_snow
        t_col = t_ds[c]
        if repart && total_precip > zero(T)
            l = col_landunit[c]
            is_glc = lun_itype[l] == istice
            all_snow_t      = is_glc ? glc_snow_t : nonglc_snow_t
            frac_rain_slope = is_glc ? glc_slope  : nonglc_slope
            frac_rain = (t_col - all_snow_t) * frac_rain_slope
            frac_rain = clamp(frac_rain, zero(T), one(T))
            rain_ds[c] = total_precip * frac_rain
            snow_ds[c] = total_precip * (one(T) - frac_rain)
        else
            rain_ds[c] = total_rain
            snow_ds[c] = total_snow
        end
    end
end

# Specific humidity at column level from vapor pressure and downscaled pbot.
@kernel function _partition_humidity_kernel!(q_ds, @Const(col_gridcell),
        @Const(vp_grc), @Const(pbot_ds), lo::Int)
    i = @index(Global)
    @inbounds begin
        T = eltype(q_ds)
        c = lo + i - 1
        g = col_gridcell[c]
        vp = vp_grc[g]
        pbot = pbot_ds[c]
        q_ds[c] = T(0.622) * vp / max(pbot - T(0.378) * vp, one(T))
    end
end

# Elevation-based longwave downscaling with conservation limits.
@kernel function _downscale_lw_kernel!(lw_ds, @Const(col_gridcell), @Const(topo_col),
        @Const(topo_grc), @Const(lw_nd), lr_lw, limit, ntopo::Int, lo::Int)
    i = @index(Global)
    @inbounds begin
        T = eltype(lw_ds)
        c = lo + i - 1
        g = col_gridcell[c]
        hsurf_g = topo_grc[g]
        hsurf_c = c <= ntopo ? topo_col[c] : hsurf_g
        lwrad_g = lw_nd[g]
        lwrad_c = lwrad_g - lr_lw * (hsurf_c - hsurf_g)
        lwrad_min = lwrad_g * (one(T) - limit)
        lwrad_max = lwrad_g * (one(T) + limit)
        lwrad_c = clamp(lwrad_c, lwrad_min, lwrad_max)
        lw_ds[c] = max(lwrad_c, zero(T))
    end
end

# --------------------------------------------------------------------------
# Public routines
# --------------------------------------------------------------------------

"""
    downscale_forcings!(bounds, a2l, col, lun, topo)

Downscale atmospheric forcings from gridcell to column level based on
topographic differences. For a single-gridcell simulation, this mostly
copies gridcell values to columns with elevation-based corrections.

Ported from `downscale_forcings` in `atm2lndMod.F90`.
"""
function downscale_forcings!(bounds::BoundsType,
                             a2l::Atm2LndData,
                             col::ColumnData,
                             lun::LandunitData,
                             topo::TopoData)
    bc_col = bounds.begc:bounds.endc
    isempty(bc_col) && return nothing
    lo = bounds.begc
    n = length(bc_col)

    FT = eltype(a2l.forc_t_downscaled_col)
    lapse_rate = a2l.params.lapse_rate
    if isnan(lapse_rate)
        lapse_rate = 0.006  # default 6 K/km
    end

    # --- Step 1: Copy gridcell values to all columns as baseline ---
    _launch!(_downscale_baseline_kernel!, a2l.forc_t_downscaled_col,
        a2l.forc_th_downscaled_col, a2l.forc_pbot_downscaled_col,
        a2l.forc_rho_downscaled_col, a2l.forc_lwrad_downscaled_col,
        a2l.forc_solad_downscaled_col, a2l.forc_solar_downscaled_col,
        col.gridcell, a2l.forc_t_not_downscaled_grc, a2l.forc_th_not_downscaled_grc,
        a2l.forc_pbot_not_downscaled_grc, a2l.forc_rho_not_downscaled_grc,
        a2l.forc_lwrad_not_downscaled_grc, a2l.forc_solad_not_downscaled_grc,
        a2l.forc_solar_not_downscaled_grc, lo, Int(NUMRAD); ndrange = n)

    # --- Step 2: Elevation-based temperature/pressure/density downscaling ---
    _launch!(_downscale_temp_kernel!, a2l.forc_t_downscaled_col,
        a2l.forc_th_downscaled_col, a2l.forc_pbot_downscaled_col,
        a2l.forc_rho_downscaled_col, col.gridcell, topo.topo_col,
        a2l.forc_topo_grc, a2l.forc_t_not_downscaled_grc,
        a2l.forc_th_not_downscaled_grc, a2l.forc_pbot_not_downscaled_grc,
        a2l.forc_rho_not_downscaled_grc, a2l.forc_hgt_grc,
        FT(lapse_rate), FT(RAIR), FT(GRAV), FT(CPAIR), lo; ndrange = n)

    # --- Step 2b: CO2/O2 partial pressures → downscaled surface pressure ---
    # forc_pco2/forc_po2 are molar_fraction * surface_pbot. The forcing reader
    # seeds them from the NOT-downscaled grc pbot; rescale to the elevation-
    # corrected pbot that the photosynthesis solve uses for `cair` and cs
    # (Fortran builds forc_pco2 from the same surface pbot → cair = 367e-6*pbot).
    # Leaving cair on the not-downscaled pbot made it ~2.7% low at elevation-
    # corrected columns → stomata ~3% too open → ~2% excess transp (Krycklan
    # deep-soil drying → water table too deep → QDRAI -10%).
    #
    # NOTE this is a CLM.jl addition: CTSM's `downscale_forcings` does not touch
    # forc_pco2/forc_po2 at all (they have no column-level counterpart there).
    # Since the target IS a gridcell field, the elevation correction it can carry
    # is the gridcell mean of the downscaled column pressure — computed here in
    # two passes so it is applied exactly ONCE per gridcell. Doing it per column
    # compounded the ratio; see the kernel comment.
    bg = bounds.begg:bounds.endg
    if !isempty(bg) && length(a2l.pbot_downscaled_grc) >= last(bg)
        fill!(a2l.pbot_downscaled_grc, zero(FT))
        fill!(a2l.pbot_downscale_wsum_grc, zero(FT))
        _launch!(_downscale_pbot_accum_kernel!, a2l.pbot_downscaled_grc,
            a2l.pbot_downscale_wsum_grc, col.gridcell, col.wtgcell,
            a2l.forc_pbot_downscaled_col, lo, bounds.endc;
            ndrange = length(col.gridcell))
        _launch!(_downscale_pco2_kernel!, a2l.forc_pco2_grc, a2l.forc_po2_grc,
            a2l.forc_pc13o2_grc, a2l.forc_pbot_not_downscaled_grc,
            a2l.pbot_downscaled_grc, a2l.pbot_downscale_wsum_grc,
            first(bg), last(bg); ndrange = length(a2l.forc_pco2_grc))
    end

    # --- Step 3: Partition precipitation ---
    partition_precip!(bounds, a2l, col, lun)

    # --- Step 4: Longwave downscaling (if enabled) ---
    if a2l.params.glcmec_downscale_longwave
        downscale_longwave!(bounds, a2l, col, lun; topo_col=topo.topo_col)
    end

    return nothing
end

"""
    partition_precip!(bounds, a2l, col, lun)

Partition total precipitation into rain and snow based on column temperature.
Uses the repartitioning parameters if `repartition_rain_snow` is enabled;
otherwise keeps the gridcell rain/snow split unchanged.

Ported from `partition_precip` in `atm2lndMod.F90`.
"""
function partition_precip!(bounds::BoundsType,
                           a2l::Atm2LndData,
                           col::ColumnData,
                           lun::LandunitData)
    bc_col = bounds.begc:bounds.endc
    isempty(bc_col) && return nothing
    lo = bounds.begc
    n = length(bc_col)
    p = a2l.params
    FT = eltype(a2l.forc_rain_downscaled_col)

    # repartition_rain_snow disabled: the column keeps the gridcell rain/snow
    # split unchanged (Fortran atm2lndMod.F90:181-186 initializes
    # forc_rain_c=forc_rain_g, forc_snow_c=forc_snow_g, and partition_precip
    # ONLY modifies them when repartition_rain_snow is .true.). Preserve the
    # gridcell split here rather than imposing a hard TFRZ threshold (which
    # would discard the DATM linear-ramp phase). NOTE: the parity harness runs
    # with repartition_rain_snow=true, so the ramp branch is active there.
    _launch!(_partition_precip_kernel!, a2l.forc_rain_downscaled_col,
        a2l.forc_snow_downscaled_col, col.gridcell, col.landunit, lun.itype,
        a2l.forc_rain_not_downscaled_grc, a2l.forc_snow_not_downscaled_grc,
        a2l.forc_t_downscaled_col, p.repartition_rain_snow,
        FT(p.precip_repartition_glc_all_snow_t), FT(p.precip_repartition_glc_frac_rain_slope),
        FT(p.precip_repartition_nonglc_all_snow_t), FT(p.precip_repartition_nonglc_frac_rain_slope),
        Int(ISTICE), length(a2l.forc_rain_not_downscaled_grc),
        length(a2l.forc_snow_not_downscaled_grc), lo; ndrange = n)

    # Compute specific humidity at column level
    _launch!(_partition_humidity_kernel!, a2l.forc_q_downscaled_col,
        col.gridcell, a2l.forc_vp_grc, a2l.forc_pbot_downscaled_col, lo; ndrange = n)

    return nothing
end

"""
    downscale_longwave!(bounds, a2l, col, lun)

Downscale longwave radiation based on elevation using linear lapse rate
with conservation normalization.

Ported from `downscale_longwave` in `atm2lndMod.F90`.
"""
function downscale_longwave!(bounds::BoundsType,
                             a2l::Atm2LndData,
                             col::ColumnData,
                             lun::LandunitData;
                             topo_col::AbstractVector{<:Real}=Float64[])
    bc_col = bounds.begc:bounds.endc
    isempty(bc_col) && return nothing
    lo = bounds.begc
    n = length(bc_col)
    p = a2l.params
    FT = eltype(a2l.forc_lwrad_downscaled_col)

    _launch!(_downscale_lw_kernel!, a2l.forc_lwrad_downscaled_col, col.gridcell,
        topo_col, a2l.forc_topo_grc, a2l.forc_lwrad_not_downscaled_grc,
        FT(p.lapse_rate_longwave), FT(p.longwave_downscaling_limit),
        length(topo_col), lo; ndrange = n)

    return nothing
end
