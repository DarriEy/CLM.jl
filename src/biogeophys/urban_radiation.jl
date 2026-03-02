# ==========================================================================
# Ported from: src/biogeophys/UrbanRadiationMod.F90
# Calculate solar and longwave radiation, and turbulent fluxes for urban landunit
#
# Public functions:
#   urban_radiation!  — Solar fluxes absorbed and reflected by roof and canyon
#   net_longwave!     — Net longwave radiation for road and walls in urban canyon
# ==========================================================================

# --- Module-level constants ---
const MPE_URBAN_RAD = 1.0e-06   # prevents overflow for division by zero
const SNOEM         = 0.97      # snow emissivity

"""
    net_longwave!(canyon_hwr, wtroad_perv, lwdown, em_roof, em_improad,
                  em_perroad, em_wall, t_roof, t_improad, t_perroad,
                  t_sunwall, t_shadewall,
                  lwnet_roof, lwnet_improad, lwnet_perroad,
                  lwnet_sunwall, lwnet_shadewall, lwnet_canyon,
                  lwup_roof, lwup_improad, lwup_perroad,
                  lwup_sunwall, lwup_shadewall, lwup_canyon,
                  urbanparams, mask_urbanl, bounds)

Net longwave radiation for road and both walls in urban canyon allowing for
multiple reflection. Also net longwave radiation for urban roof.

Ported from `net_longwave` in `UrbanRadiationMod.F90`.

# Arguments
- `canyon_hwr`    : ratio of building height to street width (landunit)
- `wtroad_perv`   : weight of pervious road wrt total road (landunit)
- `lwdown`        : atmospheric longwave radiation (W/m²) (landunit)
- `em_roof`       : roof emissivity with snow effects (landunit)
- `em_improad`    : impervious road emissivity with snow effects (landunit)
- `em_perroad`    : pervious road emissivity with snow effects (landunit)
- `em_wall`       : wall emissivity (landunit)
- `t_roof`        : roof temperature (K) (landunit)
- `t_improad`     : impervious road temperature (K) (landunit)
- `t_perroad`     : pervious road temperature (K) (landunit)
- `t_sunwall`     : sunlit wall temperature (K) (landunit)
- `t_shadewall`   : shaded wall temperature (K) (landunit)
- `lwnet_*`       : [output] net (outgoing-incoming) longwave radiation (W/m²)
- `lwup_*`        : [output] upward longwave radiation (W/m²)
- `urbanparams`   : UrbanParamsData (view factors)
- `mask_urbanl`   : BitVector mask for urban landunits
- `bounds`        : landunit index range
"""
function net_longwave!(canyon_hwr::Vector{Float64},
                        wtroad_perv::Vector{Float64},
                        lwdown::Vector{Float64},
                        em_roof::Vector{Float64},
                        em_improad::Vector{Float64},
                        em_perroad::Vector{Float64},
                        em_wall::Vector{Float64},
                        t_roof::Vector{Float64},
                        t_improad::Vector{Float64},
                        t_perroad::Vector{Float64},
                        t_sunwall::Vector{Float64},
                        t_shadewall::Vector{Float64},
                        lwnet_roof::Vector{Float64},
                        lwnet_improad::Vector{Float64},
                        lwnet_perroad::Vector{Float64},
                        lwnet_sunwall::Vector{Float64},
                        lwnet_shadewall::Vector{Float64},
                        lwnet_canyon::Vector{Float64},
                        lwup_roof::Vector{Float64},
                        lwup_improad::Vector{Float64},
                        lwup_perroad::Vector{Float64},
                        lwup_sunwall::Vector{Float64},
                        lwup_shadewall::Vector{Float64},
                        lwup_canyon::Vector{Float64},
                        urbanparams::UrbanParamsData,
                        mask_urbanl::BitVector,
                        bounds::UnitRange{Int})

    vf_sr = urbanparams.vf_sr
    vf_wr = urbanparams.vf_wr
    vf_sw = urbanparams.vf_sw
    vf_rw = urbanparams.vf_rw
    vf_ww = urbanparams.vf_ww

    n = 50  # number of iterations

    # Local working arrays (sized over bounds)
    nl = length(bounds)
    lwdown_road_v      = zeros(last(bounds))
    lwdown_sunwall_v   = zeros(last(bounds))
    lwdown_shadewall_v = zeros(last(bounds))
    lwtot_v            = zeros(last(bounds))
    wtroad_imperv_v    = zeros(last(bounds))

    improad_a_v           = zeros(last(bounds))
    improad_r_v           = zeros(last(bounds))
    improad_r_sky_v       = zeros(last(bounds))
    improad_r_sunwall_v   = zeros(last(bounds))
    improad_r_shadewall_v = zeros(last(bounds))
    improad_e_v           = zeros(last(bounds))
    improad_e_sky_v       = zeros(last(bounds))
    improad_e_sunwall_v   = zeros(last(bounds))
    improad_e_shadewall_v = zeros(last(bounds))

    perroad_a_v           = zeros(last(bounds))
    perroad_r_v           = zeros(last(bounds))
    perroad_r_sky_v       = zeros(last(bounds))
    perroad_r_sunwall_v   = zeros(last(bounds))
    perroad_r_shadewall_v = zeros(last(bounds))
    perroad_e_v           = zeros(last(bounds))
    perroad_e_sky_v       = zeros(last(bounds))
    perroad_e_sunwall_v   = zeros(last(bounds))
    perroad_e_shadewall_v = zeros(last(bounds))

    road_a_v              = zeros(last(bounds))
    road_r_v              = zeros(last(bounds))
    road_r_sky_v          = zeros(last(bounds))
    road_r_sunwall_v      = zeros(last(bounds))
    road_r_shadewall_v    = zeros(last(bounds))
    road_e_v              = zeros(last(bounds))
    road_e_sky_v          = zeros(last(bounds))
    road_e_sunwall_v      = zeros(last(bounds))
    road_e_shadewall_v    = zeros(last(bounds))

    sunwall_a_v              = zeros(last(bounds))
    sunwall_r_v              = zeros(last(bounds))
    sunwall_r_sky_v          = zeros(last(bounds))
    sunwall_r_road_v         = zeros(last(bounds))
    sunwall_r_shadewall_v    = zeros(last(bounds))
    sunwall_e_v              = zeros(last(bounds))
    sunwall_e_sky_v          = zeros(last(bounds))
    sunwall_e_road_v         = zeros(last(bounds))
    sunwall_e_shadewall_v    = zeros(last(bounds))

    shadewall_a_v            = zeros(last(bounds))
    shadewall_r_v            = zeros(last(bounds))
    shadewall_r_sky_v        = zeros(last(bounds))
    shadewall_r_road_v       = zeros(last(bounds))
    shadewall_r_sunwall_v    = zeros(last(bounds))
    shadewall_e_v            = zeros(last(bounds))
    shadewall_e_sky_v        = zeros(last(bounds))
    shadewall_e_road_v       = zeros(last(bounds))
    shadewall_e_sunwall_v    = zeros(last(bounds))

    # Calculate impervious road weight
    for l in bounds
        mask_urbanl[l] || continue
        wtroad_imperv_v[l] = 1.0 - wtroad_perv[l]
    end

    # Atmospheric longwave radiation incident on walls and road in urban canyon
    for l in bounds
        mask_urbanl[l] || continue

        lwdown_road_v[l]      = lwdown[l] * vf_sr[l]
        lwdown_sunwall_v[l]   = lwdown[l] * vf_sw[l]
        lwdown_shadewall_v[l] = lwdown[l] * vf_sw[l]

        # Conservation check
        err = lwdown[l] - (lwdown_road_v[l] + (lwdown_shadewall_v[l] + lwdown_sunwall_v[l]) * canyon_hwr[l])
        if abs(err) > 0.10
            error("urban incident atmospheric longwave radiation balance error: err=$err, l=$l")
        end
    end

    # Initial absorption, reflection, and emission
    for l in bounds
        mask_urbanl[l] || continue

        road_a_v[l] = 0.0
        road_r_v[l] = 0.0
        road_e_v[l] = 0.0

        improad_a_v[l]           =     em_improad[l]  * lwdown_road_v[l]
        improad_r_v[l]           = (1.0 - em_improad[l]) * lwdown_road_v[l]
        improad_r_sky_v[l]       = improad_r_v[l] * vf_sr[l]
        improad_r_sunwall_v[l]   = improad_r_v[l] * vf_wr[l]
        improad_r_shadewall_v[l] = improad_r_v[l] * vf_wr[l]
        improad_e_v[l]           = em_improad[l] * SB * (t_improad[l]^4)
        improad_e_sky_v[l]       = improad_e_v[l] * vf_sr[l]
        improad_e_sunwall_v[l]   = improad_e_v[l] * vf_wr[l]
        improad_e_shadewall_v[l] = improad_e_v[l] * vf_wr[l]
        road_a_v[l]              = road_a_v[l] + improad_a_v[l] * wtroad_imperv_v[l]
        road_r_v[l]              = road_r_v[l] + improad_r_v[l] * wtroad_imperv_v[l]
        road_e_v[l]              = road_e_v[l] + improad_e_v[l] * wtroad_imperv_v[l]

        perroad_a_v[l]           =     em_perroad[l]  * lwdown_road_v[l]
        perroad_r_v[l]           = (1.0 - em_perroad[l]) * lwdown_road_v[l]
        perroad_r_sky_v[l]       = perroad_r_v[l] * vf_sr[l]
        perroad_r_sunwall_v[l]   = perroad_r_v[l] * vf_wr[l]
        perroad_r_shadewall_v[l] = perroad_r_v[l] * vf_wr[l]
        perroad_e_v[l]           = em_perroad[l] * SB * (t_perroad[l]^4)
        perroad_e_sky_v[l]       = perroad_e_v[l] * vf_sr[l]
        perroad_e_sunwall_v[l]   = perroad_e_v[l] * vf_wr[l]
        perroad_e_shadewall_v[l] = perroad_e_v[l] * vf_wr[l]
        road_a_v[l]              = road_a_v[l] + perroad_a_v[l] * wtroad_perv[l]
        road_r_v[l]              = road_r_v[l] + perroad_r_v[l] * wtroad_perv[l]
        road_e_v[l]              = road_e_v[l] + perroad_e_v[l] * wtroad_perv[l]

        road_r_sky_v[l]          = road_r_v[l] * vf_sr[l]
        road_r_sunwall_v[l]      = road_r_v[l] * vf_wr[l]
        road_r_shadewall_v[l]    = road_r_v[l] * vf_wr[l]
        road_e_sky_v[l]          = road_e_v[l] * vf_sr[l]
        road_e_sunwall_v[l]      = road_e_v[l] * vf_wr[l]
        road_e_shadewall_v[l]    = road_e_v[l] * vf_wr[l]

        sunwall_a_v[l]           = em_wall[l] * lwdown_sunwall_v[l]
        sunwall_r_v[l]           = (1.0 - em_wall[l]) * lwdown_sunwall_v[l]
        sunwall_r_sky_v[l]       = sunwall_r_v[l] * vf_sw[l]
        sunwall_r_road_v[l]      = sunwall_r_v[l] * vf_rw[l]
        sunwall_r_shadewall_v[l] = sunwall_r_v[l] * vf_ww[l]
        sunwall_e_v[l]           = em_wall[l] * SB * (t_sunwall[l]^4)
        sunwall_e_sky_v[l]       = sunwall_e_v[l] * vf_sw[l]
        sunwall_e_road_v[l]      = sunwall_e_v[l] * vf_rw[l]
        sunwall_e_shadewall_v[l] = sunwall_e_v[l] * vf_ww[l]

        shadewall_a_v[l]         = em_wall[l] * lwdown_shadewall_v[l]
        shadewall_r_v[l]         = (1.0 - em_wall[l]) * lwdown_shadewall_v[l]
        shadewall_r_sky_v[l]     = shadewall_r_v[l] * vf_sw[l]
        shadewall_r_road_v[l]    = shadewall_r_v[l] * vf_rw[l]
        shadewall_r_sunwall_v[l] = shadewall_r_v[l] * vf_ww[l]
        shadewall_e_v[l]         = em_wall[l] * SB * (t_shadewall[l]^4)
        shadewall_e_sky_v[l]     = shadewall_e_v[l] * vf_sw[l]
        shadewall_e_road_v[l]    = shadewall_e_v[l] * vf_rw[l]
        shadewall_e_sunwall_v[l] = shadewall_e_v[l] * vf_ww[l]

        # Initialize sum of net and upward longwave radiation
        lwnet_improad[l]   = improad_e_v[l]   - improad_a_v[l]
        lwnet_perroad[l]   = perroad_e_v[l]   - perroad_a_v[l]
        lwnet_sunwall[l]   = sunwall_e_v[l]   - sunwall_a_v[l]
        lwnet_shadewall[l] = shadewall_e_v[l] - shadewall_a_v[l]

        lwup_improad[l]   = improad_r_sky_v[l]   + improad_e_sky_v[l]
        lwup_perroad[l]   = perroad_r_sky_v[l]   + perroad_e_sky_v[l]
        lwup_sunwall[l]   = sunwall_r_sky_v[l]   + sunwall_e_sky_v[l]
        lwup_shadewall[l] = shadewall_r_sky_v[l] + shadewall_e_sky_v[l]
    end

    # Multiple reflection iterations within canyon
    for l in bounds
        mask_urbanl[l] || continue

        converged = false
        for iter in 1:n
            # step (1): absorption and reflection
            lwtot_v[l] = (sunwall_r_road_v[l] + sunwall_e_road_v[l] +
                           shadewall_r_road_v[l] + shadewall_e_road_v[l]) * canyon_hwr[l]
            road_a_v[l]    = 0.0
            road_r_v[l]    = 0.0
            improad_r_v[l] = (1.0 - em_improad[l]) * lwtot_v[l]
            improad_a_v[l] =     em_improad[l]  * lwtot_v[l]
            road_a_v[l]    = road_a_v[l] + improad_a_v[l] * wtroad_imperv_v[l]
            road_r_v[l]    = road_r_v[l] + improad_r_v[l] * wtroad_imperv_v[l]
            perroad_r_v[l] = (1.0 - em_perroad[l]) * lwtot_v[l]
            perroad_a_v[l] =     em_perroad[l]  * lwtot_v[l]
            road_a_v[l]    = road_a_v[l] + perroad_a_v[l] * wtroad_perv[l]
            road_r_v[l]    = road_r_v[l] + perroad_r_v[l] * wtroad_perv[l]

            lwtot_v[l] = (road_r_sunwall_v[l] + road_e_sunwall_v[l]) / canyon_hwr[l] +
                           (shadewall_r_sunwall_v[l] + shadewall_e_sunwall_v[l])
            sunwall_a_v[l] =     em_wall[l]  * lwtot_v[l]
            sunwall_r_v[l] = (1.0 - em_wall[l]) * lwtot_v[l]

            lwtot_v[l] = (road_r_shadewall_v[l] + road_e_shadewall_v[l]) / canyon_hwr[l] +
                           (sunwall_r_shadewall_v[l] + sunwall_e_shadewall_v[l])
            shadewall_a_v[l] =     em_wall[l]  * lwtot_v[l]
            shadewall_r_v[l] = (1.0 - em_wall[l]) * lwtot_v[l]

            # Zero emission terms after first iteration
            sunwall_e_road_v[l]      = 0.0
            shadewall_e_road_v[l]    = 0.0
            road_e_sunwall_v[l]      = 0.0
            shadewall_e_sunwall_v[l] = 0.0
            road_e_shadewall_v[l]    = 0.0
            sunwall_e_shadewall_v[l] = 0.0

            # step (2): add net longwave for ith reflection to total
            lwnet_improad[l]   = lwnet_improad[l]   - improad_a_v[l]
            lwnet_perroad[l]   = lwnet_perroad[l]   - perroad_a_v[l]
            lwnet_sunwall[l]   = lwnet_sunwall[l]   - sunwall_a_v[l]
            lwnet_shadewall[l] = lwnet_shadewall[l] - shadewall_a_v[l]

            # step (3): distribute reflected radiation
            improad_r_sky_v[l]       = improad_r_v[l] * vf_sr[l]
            improad_r_sunwall_v[l]   = improad_r_v[l] * vf_wr[l]
            improad_r_shadewall_v[l] = improad_r_v[l] * vf_wr[l]

            perroad_r_sky_v[l]       = perroad_r_v[l] * vf_sr[l]
            perroad_r_sunwall_v[l]   = perroad_r_v[l] * vf_wr[l]
            perroad_r_shadewall_v[l] = perroad_r_v[l] * vf_wr[l]

            road_r_sky_v[l]          = road_r_v[l] * vf_sr[l]
            road_r_sunwall_v[l]      = road_r_v[l] * vf_wr[l]
            road_r_shadewall_v[l]    = road_r_v[l] * vf_wr[l]

            sunwall_r_sky_v[l]       = sunwall_r_v[l] * vf_sw[l]
            sunwall_r_road_v[l]      = sunwall_r_v[l] * vf_rw[l]
            sunwall_r_shadewall_v[l] = sunwall_r_v[l] * vf_ww[l]

            shadewall_r_sky_v[l]     = shadewall_r_v[l] * vf_sw[l]
            shadewall_r_road_v[l]    = shadewall_r_v[l] * vf_rw[l]
            shadewall_r_sunwall_v[l] = shadewall_r_v[l] * vf_ww[l]

            # step (4): add upward longwave radiation to sky
            lwup_improad[l]   = lwup_improad[l]   + improad_r_sky_v[l]
            lwup_perroad[l]   = lwup_perroad[l]   + perroad_r_sky_v[l]
            lwup_sunwall[l]   = lwup_sunwall[l]   + sunwall_r_sky_v[l]
            lwup_shadewall[l] = lwup_shadewall[l] + shadewall_r_sky_v[l]

            # step (5): convergence check
            crit = max(road_a_v[l], sunwall_a_v[l], shadewall_a_v[l])
            if crit < 0.001
                converged = true
                break
            end
        end

        if !converged
            error("urban net longwave radiation error: no convergence for landunit $l")
        end

        # Total net longwave radiation for canyon
        lwnet_canyon[l] = 0.0
        lwnet_canyon[l] = lwnet_canyon[l] + lwnet_improad[l] * wtroad_imperv_v[l]
        lwnet_canyon[l] = lwnet_canyon[l] + lwnet_perroad[l] * wtroad_perv[l]
        lwnet_canyon[l] = lwnet_canyon[l] + (lwnet_sunwall[l] + lwnet_shadewall[l]) * canyon_hwr[l]

        # Total emitted longwave for canyon
        lwup_canyon[l] = 0.0
        lwup_canyon[l] = lwup_canyon[l] + lwup_improad[l] * wtroad_imperv_v[l]
        lwup_canyon[l] = lwup_canyon[l] + lwup_perroad[l] * wtroad_perv[l]
        lwup_canyon[l] = lwup_canyon[l] + (lwup_sunwall[l] + lwup_shadewall[l]) * canyon_hwr[l]

        # Conservation check
        err = lwnet_canyon[l] - (lwup_canyon[l] - lwdown[l])
        if abs(err) > 0.10
            error("urban net longwave radiation balance error: err=$err for landunit $l")
        end
    end

    # Net longwave radiation for roof
    for l in bounds
        mask_urbanl[l] || continue
        lwup_roof[l] = em_roof[l] * SB * (t_roof[l]^4) + (1.0 - em_roof[l]) * lwdown[l]
        lwnet_roof[l] = lwup_roof[l] - lwdown[l]
    end

    return nothing
end

"""
    urban_radiation!(solarabs, energyflux, col, lun, pch, urbanparams, temperature,
                     waterdiag, forc_lwrad, forc_solad, forc_solai,
                     mask_nourbanl, mask_urbanl, mask_urbanc, mask_urbanp,
                     bounds_lun, bounds_patch)

Solar fluxes absorbed and reflected by roof and canyon (walls, road).
Also net and upward longwave fluxes.

Ported from `UrbanRadiation` in `UrbanRadiationMod.F90`.

# Arguments
- `solarabs`    : SolarAbsorbedData (input/output)
- `energyflux`  : EnergyFluxData (output: eflx_lwrad_out, eflx_lwrad_net, eflx_lwrad_net_u)
- `col`         : ColumnData (input: itype)
- `lun`         : LandunitData (input: coli, colf, canyon_hwr, wtroad_perv, gridcell)
- `pch`         : PatchData (input: column, landunit, gridcell)
- `urbanparams` : UrbanParamsData (input: em_roof, em_improad, em_perroad, em_wall, view factors)
- `temperature` : TemperatureData (input: t_grnd_col)
- `waterdiag`   : WaterDiagnosticBulkData (input: frac_sno_col)
- `forc_lwrad`  : downward longwave radiation per gridcell (Vector)
- `forc_solad`  : direct beam radiation per gridcell (Matrix: ngrc × numrad)
- `forc_solai`  : diffuse radiation per gridcell (Matrix: ngrc × numrad)
- `mask_nourbanl` : BitVector mask for non-urban landunits
- `mask_urbanl`   : BitVector mask for urban landunits
- `mask_urbanc`   : BitVector mask for urban columns (unused, kept for interface compat)
- `mask_urbanp`   : BitVector mask for urban patches
- `bounds_lun`    : landunit index range
- `bounds_patch`  : patch index range
"""
function urban_radiation!(solarabs::SolarAbsorbedData,
                            energyflux::EnergyFluxData,
                            col::ColumnData,
                            lun::LandunitData,
                            pch::PatchData,
                            urbanparams::UrbanParamsData,
                            temperature::TemperatureData,
                            waterdiag::WaterDiagnosticBulkData,
                            forc_lwrad::Vector{Float64},
                            forc_solad::Matrix{Float64},
                            forc_solai::Matrix{Float64},
                            mask_nourbanl::BitVector,
                            mask_urbanl::BitVector,
                            mask_urbanc::BitVector,
                            mask_urbanp::BitVector,
                            bounds_lun::UnitRange{Int},
                            bounds_patch::UnitRange{Int})

    nl = last(bounds_lun)

    # Local arrays (landunit-level)
    lwnet_roof      = zeros(nl)
    lwnet_improad   = zeros(nl)
    lwnet_perroad   = zeros(nl)
    lwnet_sunwall   = zeros(nl)
    lwnet_shadewall = zeros(nl)
    lwnet_canyon    = zeros(nl)
    lwup_roof       = zeros(nl)
    lwup_improad    = zeros(nl)
    lwup_perroad    = zeros(nl)
    lwup_sunwall    = zeros(nl)
    lwup_shadewall  = zeros(nl)
    lwup_canyon     = zeros(nl)
    t_roof_l        = zeros(nl)
    t_improad_l     = zeros(nl)
    t_perroad_l     = zeros(nl)
    t_sunwall_l     = zeros(nl)
    t_shadewall_l   = zeros(nl)
    lwdown_l        = zeros(nl)
    em_roof_s       = zeros(nl)
    em_improad_s    = zeros(nl)
    em_perroad_s    = zeros(nl)

    # Define fields that appear on the restart file for non-urban landunits
    for l in bounds_lun
        mask_nourbanl[l] || continue
        solarabs.sabs_roof_dir_lun[l, :]      .= SPVAL
        solarabs.sabs_roof_dif_lun[l, :]      .= SPVAL
        solarabs.sabs_sunwall_dir_lun[l, :]   .= SPVAL
        solarabs.sabs_sunwall_dif_lun[l, :]   .= SPVAL
        solarabs.sabs_shadewall_dir_lun[l, :] .= SPVAL
        solarabs.sabs_shadewall_dif_lun[l, :] .= SPVAL
        solarabs.sabs_improad_dir_lun[l, :]   .= SPVAL
        solarabs.sabs_improad_dif_lun[l, :]   .= SPVAL
        solarabs.sabs_perroad_dir_lun[l, :]   .= SPVAL
        solarabs.sabs_perroad_dif_lun[l, :]   .= SPVAL
    end

    # Set input forcing fields
    for l in bounds_lun
        mask_urbanl[l] || continue
        g = lun.gridcell[l]

        # Default temperatures
        t_roof_l[l]      = 19.0 + TFRZ
        t_sunwall_l[l]   = 19.0 + TFRZ
        t_shadewall_l[l] = 19.0 + TFRZ
        t_improad_l[l]   = 19.0 + TFRZ
        t_perroad_l[l]   = 19.0 + TFRZ

        # Initial emissivity assignment
        em_roof_s[l]    = urbanparams.em_roof[l]
        em_improad_s[l] = urbanparams.em_improad[l]
        em_perroad_s[l] = urbanparams.em_perroad[l]

        # Set urban temperatures and emissivity including snow effects
        for c in lun.coli[l]:lun.colf[l]
            if col.itype[c] == ICOL_ROOF
                t_roof_l[l] = temperature.t_grnd_col[c]
                em_roof_s[l] = urbanparams.em_roof[l] * (1.0 - waterdiag.frac_sno_col[c]) + SNOEM * waterdiag.frac_sno_col[c]
            elseif col.itype[c] == ICOL_ROAD_IMPERV
                t_improad_l[l] = temperature.t_grnd_col[c]
                em_improad_s[l] = urbanparams.em_improad[l] * (1.0 - waterdiag.frac_sno_col[c]) + SNOEM * waterdiag.frac_sno_col[c]
            elseif col.itype[c] == ICOL_ROAD_PERV
                t_perroad_l[l] = temperature.t_grnd_col[c]
                em_perroad_s[l] = urbanparams.em_perroad[l] * (1.0 - waterdiag.frac_sno_col[c]) + SNOEM * waterdiag.frac_sno_col[c]
            elseif col.itype[c] == ICOL_SUNWALL
                t_sunwall_l[l] = temperature.t_grnd_col[c]
            elseif col.itype[c] == ICOL_SHADEWALL
                t_shadewall_l[l] = temperature.t_grnd_col[c]
            end
        end
        lwdown_l[l] = forc_lwrad[g]
    end

    # Net longwave radiation for road and both walls in urban canyon
    any_urban = any(l -> mask_urbanl[l], bounds_lun)
    if any_urban
        net_longwave!(lun.canyon_hwr, lun.wtroad_perv,
                       lwdown_l, em_roof_s, em_improad_s, em_perroad_s,
                       urbanparams.em_wall,
                       t_roof_l, t_improad_l, t_perroad_l,
                       t_sunwall_l, t_shadewall_l,
                       lwnet_roof, lwnet_improad, lwnet_perroad,
                       lwnet_sunwall, lwnet_shadewall, lwnet_canyon,
                       lwup_roof, lwup_improad, lwup_perroad,
                       lwup_sunwall, lwup_shadewall, lwup_canyon,
                       urbanparams, mask_urbanl, bounds_lun)
    end

    # Determine variables needed for history output and communication with atm
    for p in bounds_patch
        mask_urbanp[p] || continue
        c = pch.column[p]
        l = pch.landunit[p]
        g = pch.gridcell[p]

        if col.itype[c] == ICOL_ROOF
            energyflux.eflx_lwrad_out_patch[p]   = lwup_roof[l]
            energyflux.eflx_lwrad_net_patch[p]   = lwnet_roof[l]
            energyflux.eflx_lwrad_net_u_patch[p] = lwnet_roof[l]
            solarabs.sabg_patch[p] = solarabs.sabs_roof_dir_lun[l, 1] * forc_solad[g, 1] +
                solarabs.sabs_roof_dif_lun[l, 1] * forc_solai[g, 1] +
                solarabs.sabs_roof_dir_lun[l, 2] * forc_solad[g, 2] +
                solarabs.sabs_roof_dif_lun[l, 2] * forc_solai[g, 2]

        elseif col.itype[c] == ICOL_SUNWALL
            energyflux.eflx_lwrad_out_patch[p]   = lwup_sunwall[l]
            energyflux.eflx_lwrad_net_patch[p]   = lwnet_sunwall[l]
            energyflux.eflx_lwrad_net_u_patch[p] = lwnet_sunwall[l]
            solarabs.sabg_patch[p] = solarabs.sabs_sunwall_dir_lun[l, 1] * forc_solad[g, 1] +
                solarabs.sabs_sunwall_dif_lun[l, 1] * forc_solai[g, 1] +
                solarabs.sabs_sunwall_dir_lun[l, 2] * forc_solad[g, 2] +
                solarabs.sabs_sunwall_dif_lun[l, 2] * forc_solai[g, 2]

        elseif col.itype[c] == ICOL_SHADEWALL
            energyflux.eflx_lwrad_out_patch[p]   = lwup_shadewall[l]
            energyflux.eflx_lwrad_net_patch[p]   = lwnet_shadewall[l]
            energyflux.eflx_lwrad_net_u_patch[p] = lwnet_shadewall[l]
            solarabs.sabg_patch[p] = solarabs.sabs_shadewall_dir_lun[l, 1] * forc_solad[g, 1] +
                solarabs.sabs_shadewall_dif_lun[l, 1] * forc_solai[g, 1] +
                solarabs.sabs_shadewall_dir_lun[l, 2] * forc_solad[g, 2] +
                solarabs.sabs_shadewall_dif_lun[l, 2] * forc_solai[g, 2]

        elseif col.itype[c] == ICOL_ROAD_PERV
            energyflux.eflx_lwrad_out_patch[p]   = lwup_perroad[l]
            energyflux.eflx_lwrad_net_patch[p]   = lwnet_perroad[l]
            energyflux.eflx_lwrad_net_u_patch[p] = lwnet_perroad[l]
            solarabs.sabg_patch[p] = solarabs.sabs_perroad_dir_lun[l, 1] * forc_solad[g, 1] +
                solarabs.sabs_perroad_dif_lun[l, 1] * forc_solai[g, 1] +
                solarabs.sabs_perroad_dir_lun[l, 2] * forc_solad[g, 2] +
                solarabs.sabs_perroad_dif_lun[l, 2] * forc_solai[g, 2]

        elseif col.itype[c] == ICOL_ROAD_IMPERV
            energyflux.eflx_lwrad_out_patch[p]   = lwup_improad[l]
            energyflux.eflx_lwrad_net_patch[p]   = lwnet_improad[l]
            energyflux.eflx_lwrad_net_u_patch[p] = lwnet_improad[l]
            solarabs.sabg_patch[p] = solarabs.sabs_improad_dir_lun[l, 1] * forc_solad[g, 1] +
                solarabs.sabs_improad_dif_lun[l, 1] * forc_solai[g, 1] +
                solarabs.sabs_improad_dir_lun[l, 2] * forc_solad[g, 2] +
                solarabs.sabs_improad_dif_lun[l, 2] * forc_solai[g, 2]
        end

        solarabs.sabv_patch[p]  = 0.0
        solarabs.fsa_patch[p]   = solarabs.sabv_patch[p] + solarabs.sabg_patch[p]
        solarabs.fsa_u_patch[p] = solarabs.fsa_patch[p]
    end

    return nothing
end
