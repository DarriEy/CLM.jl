# ==========================================================================
# Ported from: src/main/atm2lndMod.F90
# Atmospheric forcing downscaling from gridcell to column
#
# Public functions:
#   downscale_forcings!    — Main downscaling routine
#   partition_precip!      — Rain/snow partitioning from temperature
#   downscale_longwave!    — Elevation-based LW downscaling
#
# Skipped: hillslope downscaling, water tracers
# ==========================================================================

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
    bc_grc = bounds.begg:bounds.endg

    lapse_rate = a2l.params.lapse_rate
    if isnan(lapse_rate)
        lapse_rate = 0.006  # default 6 K/km
    end

    # --- Step 1: Copy gridcell values to all columns as baseline ---
    for c in bc_col
        g = col.gridcell[c]
        a2l.forc_t_downscaled_col[c]     = a2l.forc_t_not_downscaled_grc[g]
        a2l.forc_th_downscaled_col[c]    = a2l.forc_th_not_downscaled_grc[g]
        a2l.forc_pbot_downscaled_col[c]  = a2l.forc_pbot_not_downscaled_grc[g]
        a2l.forc_rho_downscaled_col[c]   = a2l.forc_rho_not_downscaled_grc[g]
        a2l.forc_lwrad_downscaled_col[c] = a2l.forc_lwrad_not_downscaled_grc[g]

        # Solar: copy direct beam
        for b in 1:NUMRAD
            a2l.forc_solad_downscaled_col[c, b] = a2l.forc_solad_not_downscaled_grc[g, b]
        end
        # Total solar
        a2l.forc_solar_downscaled_col[c] = a2l.forc_solar_not_downscaled_grc[g]
    end

    # --- Step 2: Elevation-based temperature downscaling ---
    for c in bc_col
        g = col.gridcell[c]
        hsurf_g = a2l.forc_topo_grc[g]
        hsurf_c = topo.topo_col[c]
        if isnan(hsurf_c) || isnan(hsurf_g)
            continue
        end

        tbot_g = a2l.forc_t_not_downscaled_grc[g]
        tbot_c = tbot_g - lapse_rate * (hsurf_c - hsurf_g)
        a2l.forc_t_downscaled_col[c] = tbot_c

        # Potential temperature correction
        zbot = a2l.forc_hgt_grc[g]
        Hbot = RAIR * 0.5 * (tbot_g + tbot_c) / GRAV
        if Hbot > 0.0
            a2l.forc_th_downscaled_col[c] = a2l.forc_th_not_downscaled_grc[g] +
                (tbot_c - tbot_g) * exp((zbot / Hbot) * (RAIR / CPAIR))
        end

        # Pressure via hydrostatic balance
        if Hbot > 0.0
            a2l.forc_pbot_downscaled_col[c] = a2l.forc_pbot_not_downscaled_grc[g] *
                exp(-(hsurf_c - hsurf_g) / Hbot)
        end

        # Density correction via ratio method
        pbot_g = a2l.forc_pbot_not_downscaled_grc[g]
        pbot_c = a2l.forc_pbot_downscaled_col[c]
        rho_est_g = pbot_g / (RAIR * tbot_g)
        rho_est_c = pbot_c / (RAIR * tbot_c)
        if rho_est_g > 0.0
            a2l.forc_rho_downscaled_col[c] = a2l.forc_rho_not_downscaled_grc[g] *
                (rho_est_c / rho_est_g)
        end
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
otherwise uses a simple freezing-point threshold.

Ported from `partition_precip` in `atm2lndMod.F90`.
"""
function partition_precip!(bounds::BoundsType,
                           a2l::Atm2LndData,
                           col::ColumnData,
                           lun::LandunitData)
    bc_col = bounds.begc:bounds.endc

    for c in bc_col
        g = col.gridcell[c]
        total_rain = length(a2l.forc_rain_not_downscaled_grc) >= g ?
            a2l.forc_rain_not_downscaled_grc[g] : 0.0
        total_snow = length(a2l.forc_snow_not_downscaled_grc) >= g ?
            a2l.forc_snow_not_downscaled_grc[g] : 0.0
        total_precip = total_rain + total_snow

        t_col = a2l.forc_t_downscaled_col[c]

        if a2l.params.repartition_rain_snow && total_precip > 0.0
            # Get ramp parameters based on landunit type
            l = col.landunit[c]
            is_glc = lun.itype[l] == ISTICE
            if is_glc
                all_snow_t = a2l.params.precip_repartition_glc_all_snow_t
                frac_rain_slope = a2l.params.precip_repartition_glc_frac_rain_slope
            else
                all_snow_t = a2l.params.precip_repartition_nonglc_all_snow_t
                frac_rain_slope = a2l.params.precip_repartition_nonglc_frac_rain_slope
            end

            frac_rain = (t_col - all_snow_t) * frac_rain_slope
            frac_rain = clamp(frac_rain, 0.0, 1.0)

            a2l.forc_rain_downscaled_col[c] = total_precip * frac_rain
            a2l.forc_snow_downscaled_col[c] = total_precip * (1.0 - frac_rain)
        else
            # Simple threshold: if T > TFRZ, it's rain; otherwise snow
            if t_col > TFRZ
                a2l.forc_rain_downscaled_col[c] = total_precip
                a2l.forc_snow_downscaled_col[c] = 0.0
            else
                a2l.forc_rain_downscaled_col[c] = 0.0
                a2l.forc_snow_downscaled_col[c] = total_precip
            end
        end
    end

    # Compute specific humidity at column level
    for c in bc_col
        g = col.gridcell[c]
        vp = a2l.forc_vp_grc[g]
        pbot = a2l.forc_pbot_downscaled_col[c]
        q = 0.622 * vp / max(pbot - 0.378 * vp, 1.0)
        a2l.forc_q_downscaled_col[c] = q
    end

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
                             topo_col::Vector{<:Real}=Float64[])
    bc_col = bounds.begc:bounds.endc
    bc_grc = bounds.begg:bounds.endg

    lr_lw = a2l.params.lapse_rate_longwave
    limit = a2l.params.longwave_downscaling_limit

    for c in bc_col
        g = col.gridcell[c]
        hsurf_g = a2l.forc_topo_grc[g]
        hsurf_c = length(topo_col) >= c ? topo_col[c] : hsurf_g
        lwrad_g = a2l.forc_lwrad_not_downscaled_grc[g]

        lwrad_c = lwrad_g - lr_lw * (hsurf_c - hsurf_g)

        # Bound within limit
        lwrad_min = lwrad_g * (1.0 - limit)
        lwrad_max = lwrad_g * (1.0 + limit)
        lwrad_c = clamp(lwrad_c, lwrad_min, lwrad_max)

        # Ensure non-negative
        a2l.forc_lwrad_downscaled_col[c] = max(lwrad_c, 0.0)
    end

    return nothing
end
