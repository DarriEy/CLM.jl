# ==========================================================================
# Ported from: src/biogeophys/UrbanAlbedoMod.F90
# Calculate solar radiation and albedos for urban landunit
#
# Public functions:
#   urban_albedo!     — determine urban landunit component albedos
#
# Private functions:
#   snow_albedo!      — urban snow albedos
#   incident_direct!  — direct beam solar radiation incident on walls/road
#   incident_diffuse! — diffuse solar radiation incident on walls/road
#   net_solar!        — solar radiation absorbed by road/walls with multiple reflection
# ==========================================================================

# Snow albedo constants derived from Marshall (1989) assuming soot content
# of 1.5e-5 (three times what LSM uses globally).
# Snow age effects are ignored.
const SNAL0 = 0.66  # vis albedo of urban snow
const SNAL1 = 0.56  # nir albedo of urban snow

# --------------------------------------------------------------------------
# snow_albedo!
# --------------------------------------------------------------------------

"""
    snow_albedo!(mask_urbanc, col, coszen, ind, albsn_roof, albsn_improad,
                 albsn_perroad, h2osno_total)

Determine urban snow albedos.

# Arguments
- `mask_urbanc`: BitVector mask for urban columns
- `col`: ColumnData
- `coszen`: cosine solar zenith angle (landunit-indexed)
- `ind`: 0=direct beam, 1=diffuse radiation
- `albsn_roof`: output snow albedo for roof (landunit × numrad)
- `albsn_improad`: output snow albedo for impervious road (landunit × numrad)
- `albsn_perroad`: output snow albedo for pervious road (landunit × numrad)
- `h2osno_total`: total snow water per column (mm H2O)

Ported from `SnowAlbedo` in `UrbanAlbedoMod.F90`.
"""
function snow_albedo!(mask_urbanc::BitVector, col::ColumnData,
                      coszen::Vector{Float64}, ind::Int,
                      albsn_roof::Matrix{Float64},
                      albsn_improad::Matrix{Float64},
                      albsn_perroad::Matrix{Float64},
                      h2osno_total::Vector{Float64})

    for c in eachindex(mask_urbanc)
        mask_urbanc[c] || continue
        l = col.landunit[c]
        if coszen[l] > 0.0 && h2osno_total[c] > 0.0
            if col.itype[c] == ICOL_ROOF
                albsn_roof[l, 1] = SNAL0
                albsn_roof[l, 2] = SNAL1
            elseif col.itype[c] == ICOL_ROAD_IMPERV
                albsn_improad[l, 1] = SNAL0
                albsn_improad[l, 2] = SNAL1
            elseif col.itype[c] == ICOL_ROAD_PERV
                albsn_perroad[l, 1] = SNAL0
                albsn_perroad[l, 2] = SNAL1
            end
        else
            if col.itype[c] == ICOL_ROOF
                albsn_roof[l, 1] = 0.0
                albsn_roof[l, 2] = 0.0
            elseif col.itype[c] == ICOL_ROAD_IMPERV
                albsn_improad[l, 1] = 0.0
                albsn_improad[l, 2] = 0.0
            elseif col.itype[c] == ICOL_ROAD_PERV
                albsn_perroad[l, 1] = 0.0
                albsn_perroad[l, 2] = 0.0
            end
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# incident_direct!
# --------------------------------------------------------------------------

"""
    incident_direct!(mask_urbanl, canyon_hwr, coszen, zen, sdir,
                     sdir_road, sdir_sunwall, sdir_shadewall)

Direct beam solar radiation incident on walls and road in urban canyon.

Conservation check: Total incoming direct beam (sdir) =
    sdir_road + (sdir_shadewall + sdir_sunwall) * canyon_hwr

Source: Masson, V. (2000) A physically-based scheme for the urban energy
budget in atmospheric models. Boundary-Layer Meteorology 94:357-397.

Ported from `incident_direct` in `UrbanAlbedoMod.F90`.
"""
function incident_direct!(mask_urbanl::BitVector,
                           canyon_hwr::Vector{Float64},
                           coszen::Vector{Float64},
                           zen::Vector{Float64},
                           sdir::Matrix{Float64},
                           sdir_road::Matrix{Float64},
                           sdir_sunwall::Matrix{Float64},
                           sdir_shadewall::Matrix{Float64})

    nl = length(mask_urbanl)
    theta0 = zeros(nl)
    tanzen = zeros(nl)

    for l in eachindex(mask_urbanl)
        mask_urbanl[l] || continue
        if coszen[l] > 0.0
            theta0[l] = asin(min(1.0 / (canyon_hwr[l] * tan(max(zen[l], 0.000001))), 1.0))
            tanzen[l] = tan(zen[l])
        end
    end

    for ib in 1:NUMRAD
        for l in eachindex(mask_urbanl)
            mask_urbanl[l] || continue
            if coszen[l] > 0.0
                sdir_shadewall[l, ib] = 0.0

                # incident solar radiation on wall and road integrated over all
                # canyon orientations (0 <= theta <= pi/2)
                sdir_road[l, ib] = sdir[l, ib] *
                    (2.0 * theta0[l] / RPI - 2.0 / RPI * canyon_hwr[l] * tanzen[l] * (1.0 - cos(theta0[l])))
                sdir_sunwall[l, ib] = 2.0 * sdir[l, ib] * ((1.0 / canyon_hwr[l]) *
                    (0.5 - theta0[l] / RPI) + (1.0 / RPI) * tanzen[l] * (1.0 - cos(theta0[l])))

                # conservation check
                swall_projected = (sdir_shadewall[l, ib] + sdir_sunwall[l, ib]) * canyon_hwr[l]
                err1 = sdir[l, ib] - (sdir_road[l, ib] + swall_projected)
                if abs(err1) > 0.001
                    error("urban direct beam solar radiation balance error: $err1")
                end
            else
                sdir_road[l, ib] = 0.0
                sdir_sunwall[l, ib] = 0.0
                sdir_shadewall[l, ib] = 0.0
            end
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# incident_diffuse!
# --------------------------------------------------------------------------

"""
    incident_diffuse!(mask_urbanl, canyon_hwr, sdif, sdif_road,
                      sdif_sunwall, sdif_shadewall, urbanparams)

Diffuse solar radiation incident on walls and road in urban canyon.

Conservation check: Total incoming diffuse (sdif) =
    sdif_road + (sdif_shadewall + sdif_sunwall) * canyon_hwr

Ported from `incident_diffuse` in `UrbanAlbedoMod.F90`.
"""
function incident_diffuse!(mask_urbanl::BitVector,
                            canyon_hwr::Vector{Float64},
                            sdif::Matrix{Float64},
                            sdif_road::Matrix{Float64},
                            sdif_sunwall::Matrix{Float64},
                            sdif_shadewall::Matrix{Float64},
                            urbanparams::UrbanParamsData)

    vf_sr = urbanparams.vf_sr
    vf_sw = urbanparams.vf_sw

    for ib in 1:NUMRAD
        for l in eachindex(mask_urbanl)
            mask_urbanl[l] || continue
            sdif_road[l, ib]      = sdif[l, ib] * vf_sr[l]
            sdif_sunwall[l, ib]   = sdif[l, ib] * vf_sw[l]
            sdif_shadewall[l, ib] = sdif[l, ib] * vf_sw[l]

            swall_projected = (sdif_shadewall[l, ib] + sdif_sunwall[l, ib]) * canyon_hwr[l]
            err = sdif[l, ib] - (sdif_road[l, ib] + swall_projected)
            if abs(err) > 0.001
                error("urban diffuse solar radiation balance error: $err")
            end
        end
    end

    return nothing
end

# --------------------------------------------------------------------------
# net_solar!
# --------------------------------------------------------------------------

"""
    net_solar!(mask_urbanl, coszen, canyon_hwr, wtroad_perv,
               sdir, sdif,
               alb_improad_dir, alb_perroad_dir, alb_wall_dir, alb_roof_dir,
               alb_improad_dif, alb_perroad_dif, alb_wall_dif, alb_roof_dif,
               sdir_road, sdir_sunwall, sdir_shadewall,
               sdif_road, sdif_sunwall, sdif_shadewall,
               sref_improad_dir, sref_perroad_dir, sref_sunwall_dir,
               sref_shadewall_dir, sref_roof_dir,
               sref_improad_dif, sref_perroad_dif, sref_sunwall_dif,
               sref_shadewall_dif, sref_roof_dif,
               urbanparams, solarabs)

Solar radiation absorbed by road and both walls in urban canyon allowing
for multiple reflection.

Ported from `net_solar` in `UrbanAlbedoMod.F90`.
"""
function net_solar!(mask_urbanl::BitVector,
                    coszen::Vector{Float64},
                    canyon_hwr::Vector{Float64},
                    wtroad_perv::Vector{Float64},
                    sdir::Matrix{Float64},
                    sdif::Matrix{Float64},
                    alb_improad_dir::Matrix{Float64},
                    alb_perroad_dir::Matrix{Float64},
                    alb_wall_dir::Matrix{Float64},
                    alb_roof_dir::Matrix{Float64},
                    alb_improad_dif::Matrix{Float64},
                    alb_perroad_dif::Matrix{Float64},
                    alb_wall_dif::Matrix{Float64},
                    alb_roof_dif::Matrix{Float64},
                    sdir_road::Matrix{Float64},
                    sdir_sunwall::Matrix{Float64},
                    sdir_shadewall::Matrix{Float64},
                    sdif_road::Matrix{Float64},
                    sdif_sunwall::Matrix{Float64},
                    sdif_shadewall::Matrix{Float64},
                    sref_improad_dir::Matrix{Float64},
                    sref_perroad_dir::Matrix{Float64},
                    sref_sunwall_dir::Matrix{Float64},
                    sref_shadewall_dir::Matrix{Float64},
                    sref_roof_dir::Matrix{Float64},
                    sref_improad_dif::Matrix{Float64},
                    sref_perroad_dif::Matrix{Float64},
                    sref_sunwall_dif::Matrix{Float64},
                    sref_shadewall_dif::Matrix{Float64},
                    sref_roof_dif::Matrix{Float64},
                    urbanparams::UrbanParamsData,
                    solarabs::SolarAbsorbedData)

    vf_sr = urbanparams.vf_sr
    vf_wr = urbanparams.vf_wr
    vf_sw = urbanparams.vf_sw
    vf_rw = urbanparams.vf_rw
    vf_ww = urbanparams.vf_ww

    sabs_roof_dir      = solarabs.sabs_roof_dir_lun
    sabs_roof_dif      = solarabs.sabs_roof_dif_lun
    sabs_sunwall_dir   = solarabs.sabs_sunwall_dir_lun
    sabs_sunwall_dif   = solarabs.sabs_sunwall_dif_lun
    sabs_shadewall_dir = solarabs.sabs_shadewall_dir_lun
    sabs_shadewall_dif = solarabs.sabs_shadewall_dif_lun
    sabs_improad_dir   = solarabs.sabs_improad_dir_lun
    sabs_improad_dif   = solarabs.sabs_improad_dif_lun
    sabs_perroad_dir   = solarabs.sabs_perroad_dir_lun
    sabs_perroad_dif   = solarabs.sabs_perroad_dif_lun

    nl = length(mask_urbanl)
    n = 50              # number of iterations for multiple reflection
    errcrit = 0.00001   # error criteria for convergence

    # Per-landunit working arrays
    wtroad_imperv = zeros(nl)

    improad_a_dir = zeros(nl);  improad_r_dir = zeros(nl)
    improad_r_sky_dir = zeros(nl); improad_r_sunwall_dir = zeros(nl); improad_r_shadewall_dir = zeros(nl)
    improad_a_dif = zeros(nl);  improad_r_dif = zeros(nl)
    improad_r_sky_dif = zeros(nl); improad_r_sunwall_dif = zeros(nl); improad_r_shadewall_dif = zeros(nl)

    perroad_a_dir = zeros(nl);  perroad_r_dir = zeros(nl)
    perroad_r_sky_dir = zeros(nl); perroad_r_sunwall_dir = zeros(nl); perroad_r_shadewall_dir = zeros(nl)
    perroad_a_dif = zeros(nl);  perroad_r_dif = zeros(nl)
    perroad_r_sky_dif = zeros(nl); perroad_r_sunwall_dif = zeros(nl); perroad_r_shadewall_dif = zeros(nl)

    road_a_dir = zeros(nl); road_r_dir = zeros(nl)
    road_r_sky_dir = zeros(nl); road_r_sunwall_dir = zeros(nl); road_r_shadewall_dir = zeros(nl)
    road_a_dif = zeros(nl); road_r_dif = zeros(nl)
    road_r_sky_dif = zeros(nl); road_r_sunwall_dif = zeros(nl); road_r_shadewall_dif = zeros(nl)

    sunwall_a_dir = zeros(nl); sunwall_r_dir = zeros(nl)
    sunwall_r_sky_dir = zeros(nl); sunwall_r_road_dir = zeros(nl); sunwall_r_shadewall_dir = zeros(nl)
    sunwall_a_dif = zeros(nl); sunwall_r_dif = zeros(nl)
    sunwall_r_sky_dif = zeros(nl); sunwall_r_road_dif = zeros(nl); sunwall_r_shadewall_dif = zeros(nl)

    shadewall_a_dir = zeros(nl); shadewall_r_dir = zeros(nl)
    shadewall_r_sky_dir = zeros(nl); shadewall_r_road_dir = zeros(nl); shadewall_r_sunwall_dir = zeros(nl)
    shadewall_a_dif = zeros(nl); shadewall_r_dif = zeros(nl)
    shadewall_r_sky_dif = zeros(nl); shadewall_r_road_dif = zeros(nl); shadewall_r_sunwall_dif = zeros(nl)

    stot = zeros(nl)
    sref_canyon_dir = zeros(nl); sref_canyon_dif = zeros(nl)
    sabs_canyon_dir = zeros(nl); sabs_canyon_dif = zeros(nl)
    stot_dir = zeros(nl); stot_dif = zeros(nl)

    # Calculate impervious road weight
    for l in eachindex(mask_urbanl)
        mask_urbanl[l] || continue
        wtroad_imperv[l] = 1.0 - wtroad_perv[l]
    end

    for ib in 1:NUMRAD
        for l in eachindex(mask_urbanl)
            mask_urbanl[l] || continue
            if coszen[l] > 0.0

                # ---- Initial absorption and reflection ----
                # Direct beam
                road_a_dir[l]              = 0.0
                road_r_dir[l]              = 0.0
                improad_a_dir[l]           = (1.0 - alb_improad_dir[l, ib]) * sdir_road[l, ib]
                improad_r_dir[l]           =        alb_improad_dir[l, ib]  * sdir_road[l, ib]
                improad_r_sky_dir[l]       = improad_r_dir[l] * vf_sr[l]
                improad_r_sunwall_dir[l]   = improad_r_dir[l] * vf_wr[l]
                improad_r_shadewall_dir[l] = improad_r_dir[l] * vf_wr[l]
                road_a_dir[l]              = road_a_dir[l] + improad_a_dir[l] * wtroad_imperv[l]
                road_r_dir[l]              = road_r_dir[l] + improad_r_dir[l] * wtroad_imperv[l]

                perroad_a_dir[l]           = (1.0 - alb_perroad_dir[l, ib]) * sdir_road[l, ib]
                perroad_r_dir[l]           =        alb_perroad_dir[l, ib]  * sdir_road[l, ib]
                perroad_r_sky_dir[l]       = perroad_r_dir[l] * vf_sr[l]
                perroad_r_sunwall_dir[l]   = perroad_r_dir[l] * vf_wr[l]
                perroad_r_shadewall_dir[l] = perroad_r_dir[l] * vf_wr[l]
                road_a_dir[l]              = road_a_dir[l] + perroad_a_dir[l] * wtroad_perv[l]
                road_r_dir[l]              = road_r_dir[l] + perroad_r_dir[l] * wtroad_perv[l]

                road_r_sky_dir[l]          = road_r_dir[l] * vf_sr[l]
                road_r_sunwall_dir[l]      = road_r_dir[l] * vf_wr[l]
                road_r_shadewall_dir[l]    = road_r_dir[l] * vf_wr[l]

                sunwall_a_dir[l]           = (1.0 - alb_wall_dir[l, ib]) * sdir_sunwall[l, ib]
                sunwall_r_dir[l]           =        alb_wall_dir[l, ib]  * sdir_sunwall[l, ib]
                sunwall_r_sky_dir[l]       = sunwall_r_dir[l] * vf_sw[l]
                sunwall_r_road_dir[l]      = sunwall_r_dir[l] * vf_rw[l]
                sunwall_r_shadewall_dir[l] = sunwall_r_dir[l] * vf_ww[l]

                shadewall_a_dir[l]         = (1.0 - alb_wall_dir[l, ib]) * sdir_shadewall[l, ib]
                shadewall_r_dir[l]         =        alb_wall_dir[l, ib]  * sdir_shadewall[l, ib]
                shadewall_r_sky_dir[l]     = shadewall_r_dir[l] * vf_sw[l]
                shadewall_r_road_dir[l]    = shadewall_r_dir[l] * vf_rw[l]
                shadewall_r_sunwall_dir[l] = shadewall_r_dir[l] * vf_ww[l]

                # Diffuse
                road_a_dif[l]              = 0.0
                road_r_dif[l]              = 0.0
                improad_a_dif[l]           = (1.0 - alb_improad_dif[l, ib]) * sdif_road[l, ib]
                improad_r_dif[l]           =        alb_improad_dif[l, ib]  * sdif_road[l, ib]
                improad_r_sky_dif[l]       = improad_r_dif[l] * vf_sr[l]
                improad_r_sunwall_dif[l]   = improad_r_dif[l] * vf_wr[l]
                improad_r_shadewall_dif[l] = improad_r_dif[l] * vf_wr[l]
                road_a_dif[l]              = road_a_dif[l] + improad_a_dif[l] * wtroad_imperv[l]
                road_r_dif[l]              = road_r_dif[l] + improad_r_dif[l] * wtroad_imperv[l]

                perroad_a_dif[l]           = (1.0 - alb_perroad_dif[l, ib]) * sdif_road[l, ib]
                perroad_r_dif[l]           =        alb_perroad_dif[l, ib]  * sdif_road[l, ib]
                perroad_r_sky_dif[l]       = perroad_r_dif[l] * vf_sr[l]
                perroad_r_sunwall_dif[l]   = perroad_r_dif[l] * vf_wr[l]
                perroad_r_shadewall_dif[l] = perroad_r_dif[l] * vf_wr[l]
                road_a_dif[l]              = road_a_dif[l] + perroad_a_dif[l] * wtroad_perv[l]
                road_r_dif[l]              = road_r_dif[l] + perroad_r_dif[l] * wtroad_perv[l]

                road_r_sky_dif[l]          = road_r_dif[l] * vf_sr[l]
                road_r_sunwall_dif[l]      = road_r_dif[l] * vf_wr[l]
                road_r_shadewall_dif[l]    = road_r_dif[l] * vf_wr[l]

                sunwall_a_dif[l]           = (1.0 - alb_wall_dif[l, ib]) * sdif_sunwall[l, ib]
                sunwall_r_dif[l]           =        alb_wall_dif[l, ib]  * sdif_sunwall[l, ib]
                sunwall_r_sky_dif[l]       = sunwall_r_dif[l] * vf_sw[l]
                sunwall_r_road_dif[l]      = sunwall_r_dif[l] * vf_rw[l]
                sunwall_r_shadewall_dif[l] = sunwall_r_dif[l] * vf_ww[l]

                shadewall_a_dif[l]         = (1.0 - alb_wall_dif[l, ib]) * sdif_shadewall[l, ib]
                shadewall_r_dif[l]         =        alb_wall_dif[l, ib]  * sdif_shadewall[l, ib]
                shadewall_r_sky_dif[l]     = shadewall_r_dif[l] * vf_sw[l]
                shadewall_r_road_dif[l]    = shadewall_r_dif[l] * vf_rw[l]
                shadewall_r_sunwall_dif[l] = shadewall_r_dif[l] * vf_ww[l]

                # Initialize sums of absorption and reflection
                sabs_improad_dir[l, ib]   = improad_a_dir[l]
                sabs_perroad_dir[l, ib]   = perroad_a_dir[l]
                sabs_sunwall_dir[l, ib]   = sunwall_a_dir[l]
                sabs_shadewall_dir[l, ib] = shadewall_a_dir[l]

                sabs_improad_dif[l, ib]   = improad_a_dif[l]
                sabs_perroad_dif[l, ib]   = perroad_a_dif[l]
                sabs_sunwall_dif[l, ib]   = sunwall_a_dif[l]
                sabs_shadewall_dif[l, ib] = shadewall_a_dif[l]

                sref_improad_dir[l, ib]   = improad_r_sky_dir[l]
                sref_perroad_dir[l, ib]   = perroad_r_sky_dir[l]
                sref_sunwall_dir[l, ib]   = sunwall_r_sky_dir[l]
                sref_shadewall_dir[l, ib] = shadewall_r_sky_dir[l]

                sref_improad_dif[l, ib]   = improad_r_sky_dif[l]
                sref_perroad_dif[l, ib]   = perroad_r_sky_dif[l]
                sref_sunwall_dif[l, ib]   = sunwall_r_sky_dif[l]
                sref_shadewall_dif[l, ib] = shadewall_r_sky_dif[l]
            end
        end

        # ---- Multiple reflections ----
        for l in eachindex(mask_urbanl)
            mask_urbanl[l] || continue
            if coszen[l] > 0.0

                # Reflected direct beam
                iter_dir = 0
                for iter in 1:n
                    iter_dir = iter
                    # step (1)
                    stot[l] = (sunwall_r_road_dir[l] + shadewall_r_road_dir[l]) * canyon_hwr[l]

                    road_a_dir[l] = 0.0
                    road_r_dir[l] = 0.0
                    improad_a_dir[l] = (1.0 - alb_improad_dir[l, ib]) * stot[l]
                    improad_r_dir[l] =        alb_improad_dir[l, ib]  * stot[l]
                    road_a_dir[l]    = road_a_dir[l] + improad_a_dir[l] * wtroad_imperv[l]
                    road_r_dir[l]    = road_r_dir[l] + improad_r_dir[l] * wtroad_imperv[l]
                    perroad_a_dir[l] = (1.0 - alb_perroad_dir[l, ib]) * stot[l]
                    perroad_r_dir[l] =        alb_perroad_dir[l, ib]  * stot[l]
                    road_a_dir[l]    = road_a_dir[l] + perroad_a_dir[l] * wtroad_perv[l]
                    road_r_dir[l]    = road_r_dir[l] + perroad_r_dir[l] * wtroad_perv[l]

                    stot[l] = road_r_sunwall_dir[l] / canyon_hwr[l] + shadewall_r_sunwall_dir[l]
                    sunwall_a_dir[l] = (1.0 - alb_wall_dir[l, ib]) * stot[l]
                    sunwall_r_dir[l] =        alb_wall_dir[l, ib]  * stot[l]

                    stot[l] = road_r_shadewall_dir[l] / canyon_hwr[l] + sunwall_r_shadewall_dir[l]
                    shadewall_a_dir[l] = (1.0 - alb_wall_dir[l, ib]) * stot[l]
                    shadewall_r_dir[l] =        alb_wall_dir[l, ib]  * stot[l]

                    # step (2)
                    sabs_improad_dir[l, ib]   += improad_a_dir[l]
                    sabs_perroad_dir[l, ib]   += perroad_a_dir[l]
                    sabs_sunwall_dir[l, ib]   += sunwall_a_dir[l]
                    sabs_shadewall_dir[l, ib] += shadewall_a_dir[l]

                    # step (3)
                    improad_r_sky_dir[l]       = improad_r_dir[l] * vf_sr[l]
                    improad_r_sunwall_dir[l]   = improad_r_dir[l] * vf_wr[l]
                    improad_r_shadewall_dir[l] = improad_r_dir[l] * vf_wr[l]

                    perroad_r_sky_dir[l]       = perroad_r_dir[l] * vf_sr[l]
                    perroad_r_sunwall_dir[l]   = perroad_r_dir[l] * vf_wr[l]
                    perroad_r_shadewall_dir[l] = perroad_r_dir[l] * vf_wr[l]

                    road_r_sky_dir[l]          = road_r_dir[l] * vf_sr[l]
                    road_r_sunwall_dir[l]      = road_r_dir[l] * vf_wr[l]
                    road_r_shadewall_dir[l]    = road_r_dir[l] * vf_wr[l]

                    sunwall_r_sky_dir[l]       = sunwall_r_dir[l] * vf_sw[l]
                    sunwall_r_road_dir[l]      = sunwall_r_dir[l] * vf_rw[l]
                    sunwall_r_shadewall_dir[l] = sunwall_r_dir[l] * vf_ww[l]

                    shadewall_r_sky_dir[l]     = shadewall_r_dir[l] * vf_sw[l]
                    shadewall_r_road_dir[l]    = shadewall_r_dir[l] * vf_rw[l]
                    shadewall_r_sunwall_dir[l] = shadewall_r_dir[l] * vf_ww[l]

                    # step (4)
                    sref_improad_dir[l, ib]   += improad_r_sky_dir[l]
                    sref_perroad_dir[l, ib]   += perroad_r_sky_dir[l]
                    sref_sunwall_dir[l, ib]   += sunwall_r_sky_dir[l]
                    sref_shadewall_dir[l, ib] += shadewall_r_sky_dir[l]

                    # step (5)
                    crit = max(road_a_dir[l], sunwall_a_dir[l], shadewall_a_dir[l])
                    if crit < errcrit
                        break
                    end
                end
                if iter_dir >= n
                    error("urban net solar radiation error: no convergence, direct beam")
                end

                # Reflected diffuse
                iter_dif = 0
                for iter in 1:n
                    iter_dif = iter
                    # step (1)
                    stot[l] = (sunwall_r_road_dif[l] + shadewall_r_road_dif[l]) * canyon_hwr[l]
                    road_a_dif[l]    = 0.0
                    road_r_dif[l]    = 0.0
                    improad_a_dif[l] = (1.0 - alb_improad_dif[l, ib]) * stot[l]
                    improad_r_dif[l] =        alb_improad_dif[l, ib]  * stot[l]
                    road_a_dif[l]    = road_a_dif[l] + improad_a_dif[l] * wtroad_imperv[l]
                    road_r_dif[l]    = road_r_dif[l] + improad_r_dif[l] * wtroad_imperv[l]
                    perroad_a_dif[l] = (1.0 - alb_perroad_dif[l, ib]) * stot[l]
                    perroad_r_dif[l] =        alb_perroad_dif[l, ib]  * stot[l]
                    road_a_dif[l]    = road_a_dif[l] + perroad_a_dif[l] * wtroad_perv[l]
                    road_r_dif[l]    = road_r_dif[l] + perroad_r_dif[l] * wtroad_perv[l]

                    stot[l] = road_r_sunwall_dif[l] / canyon_hwr[l] + shadewall_r_sunwall_dif[l]
                    sunwall_a_dif[l] = (1.0 - alb_wall_dif[l, ib]) * stot[l]
                    sunwall_r_dif[l] =        alb_wall_dif[l, ib]  * stot[l]

                    stot[l] = road_r_shadewall_dif[l] / canyon_hwr[l] + sunwall_r_shadewall_dif[l]
                    shadewall_a_dif[l] = (1.0 - alb_wall_dif[l, ib]) * stot[l]
                    shadewall_r_dif[l] =        alb_wall_dif[l, ib]  * stot[l]

                    # step (2)
                    sabs_improad_dif[l, ib]   += improad_a_dif[l]
                    sabs_perroad_dif[l, ib]   += perroad_a_dif[l]
                    sabs_sunwall_dif[l, ib]   += sunwall_a_dif[l]
                    sabs_shadewall_dif[l, ib] += shadewall_a_dif[l]

                    # step (3)
                    improad_r_sky_dif[l]       = improad_r_dif[l] * vf_sr[l]
                    improad_r_sunwall_dif[l]   = improad_r_dif[l] * vf_wr[l]
                    improad_r_shadewall_dif[l] = improad_r_dif[l] * vf_wr[l]

                    perroad_r_sky_dif[l]       = perroad_r_dif[l] * vf_sr[l]
                    perroad_r_sunwall_dif[l]   = perroad_r_dif[l] * vf_wr[l]
                    perroad_r_shadewall_dif[l] = perroad_r_dif[l] * vf_wr[l]

                    road_r_sky_dif[l]          = road_r_dif[l] * vf_sr[l]
                    road_r_sunwall_dif[l]      = road_r_dif[l] * vf_wr[l]
                    road_r_shadewall_dif[l]    = road_r_dif[l] * vf_wr[l]

                    sunwall_r_sky_dif[l]       = sunwall_r_dif[l] * vf_sw[l]
                    sunwall_r_road_dif[l]      = sunwall_r_dif[l] * vf_rw[l]
                    sunwall_r_shadewall_dif[l] = sunwall_r_dif[l] * vf_ww[l]

                    shadewall_r_sky_dif[l]     = shadewall_r_dif[l] * vf_sw[l]
                    shadewall_r_road_dif[l]    = shadewall_r_dif[l] * vf_rw[l]
                    shadewall_r_sunwall_dif[l] = shadewall_r_dif[l] * vf_ww[l]

                    # step (4)
                    sref_improad_dif[l, ib]   += improad_r_sky_dif[l]
                    sref_perroad_dif[l, ib]   += perroad_r_sky_dif[l]
                    sref_sunwall_dif[l, ib]   += sunwall_r_sky_dif[l]
                    sref_shadewall_dif[l, ib] += shadewall_r_sky_dif[l]

                    # step (5)
                    crit = max(road_a_dif[l], sunwall_a_dif[l], shadewall_a_dif[l])
                    if crit < errcrit
                        break
                    end
                end
                if iter_dif >= n
                    error("urban net solar radiation error: no convergence, diffuse")
                end

                # ---- Total reflected and absorbed by canyon ----
                sref_canyon_dir[l] = 0.0
                sref_canyon_dif[l] = 0.0
                sref_canyon_dir[l] += sref_improad_dir[l, ib] * wtroad_imperv[l]
                sref_canyon_dif[l] += sref_improad_dif[l, ib] * wtroad_imperv[l]
                sref_canyon_dir[l] += sref_perroad_dir[l, ib] * wtroad_perv[l]
                sref_canyon_dif[l] += sref_perroad_dif[l, ib] * wtroad_perv[l]
                sref_canyon_dir[l] += (sref_sunwall_dir[l, ib] + sref_shadewall_dir[l, ib]) * canyon_hwr[l]
                sref_canyon_dif[l] += (sref_sunwall_dif[l, ib] + sref_shadewall_dif[l, ib]) * canyon_hwr[l]

                sabs_canyon_dir[l] = 0.0
                sabs_canyon_dif[l] = 0.0
                sabs_canyon_dir[l] += sabs_improad_dir[l, ib] * wtroad_imperv[l]
                sabs_canyon_dif[l] += sabs_improad_dif[l, ib] * wtroad_imperv[l]
                sabs_canyon_dir[l] += sabs_perroad_dir[l, ib] * wtroad_perv[l]
                sabs_canyon_dif[l] += sabs_perroad_dif[l, ib] * wtroad_perv[l]
                sabs_canyon_dir[l] += (sabs_sunwall_dir[l, ib] + sabs_shadewall_dir[l, ib]) * canyon_hwr[l]
                sabs_canyon_dif[l] += (sabs_sunwall_dif[l, ib] + sabs_shadewall_dif[l, ib]) * canyon_hwr[l]

                # Conservation check
                stot_dir[l] = sdir_road[l, ib] + (sdir_sunwall[l, ib] + sdir_shadewall[l, ib]) * canyon_hwr[l]
                stot_dif[l] = sdif_road[l, ib] + (sdif_sunwall[l, ib] + sdif_shadewall[l, ib]) * canyon_hwr[l]

                err = stot_dir[l] + stot_dif[l] -
                    (sabs_canyon_dir[l] + sabs_canyon_dif[l] + sref_canyon_dir[l] + sref_canyon_dif[l])
                if abs(err) > 0.001
                    error("urban net solar radiation balance error for ib=$ib err=$err")
                end
            end
        end

        # Reflected and absorbed solar radiation for roof
        for l in eachindex(mask_urbanl)
            mask_urbanl[l] || continue
            if coszen[l] > 0.0
                sref_roof_dir[l, ib] = alb_roof_dir[l, ib] * sdir[l, ib]
                sref_roof_dif[l, ib] = alb_roof_dif[l, ib] * sdif[l, ib]
                sabs_roof_dir[l, ib] = sdir[l, ib] - sref_roof_dir[l, ib]
                sabs_roof_dif[l, ib] = sdif[l, ib] - sref_roof_dif[l, ib]
            end
        end

    end  # ib loop

    return nothing
end

# --------------------------------------------------------------------------
# urban_albedo!
# --------------------------------------------------------------------------

"""
    urban_albedo!(mask_urbanl, mask_urbanc, mask_urbanp,
                  lun, col, pch,
                  waterstatebulk, waterdiagnosticbulk,
                  urbanparams, solarabs, surfalb)

Determine urban landunit component albedos.

This routine is called with "inactive_and_active" filters because the
variables computed here are needed over inactive points that might later
become active (due to landuse change).

Ported from `UrbanAlbedo` in `UrbanAlbedoMod.F90`.
"""
function urban_albedo!(mask_urbanl::BitVector,
                       mask_urbanc::BitVector,
                       mask_urbanp::BitVector,
                       lun::LandunitData,
                       col::ColumnData,
                       pch::PatchData,
                       waterstatebulk::WaterStateBulkData,
                       waterdiagnosticbulk::WaterDiagnosticBulkData,
                       urbanparams::UrbanParamsData,
                       solarabs::SolarAbsorbedData,
                       surfalb::SurfaceAlbedoData)

    nl = length(mask_urbanl)
    numrad = NUMRAD

    # Aliases for output arrays
    albgrd     = surfalb.albgrd_col
    albgri     = surfalb.albgri_col
    albd       = surfalb.albd_patch
    albi       = surfalb.albi_patch
    fabd       = surfalb.fabd_patch
    fabd_sun   = surfalb.fabd_sun_patch
    fabd_sha   = surfalb.fabd_sha_patch
    fabi       = surfalb.fabi_patch
    fabi_sun   = surfalb.fabi_sun_patch
    fabi_sha   = surfalb.fabi_sha_patch
    ftdd       = surfalb.ftdd_patch
    ftid       = surfalb.ftid_patch
    ftii       = surfalb.ftii_patch
    albd_hst   = surfalb.albd_hst_patch
    albi_hst   = surfalb.albi_hst_patch
    albgrd_hst = surfalb.albgrd_hst_col
    albgri_hst = surfalb.albgri_hst_col

    frac_sno   = waterdiagnosticbulk.frac_sno_col

    vf_sr      = urbanparams.vf_sr
    vf_sw      = urbanparams.vf_sw

    # ---- Cosine solar zenith angle ----
    coszen = zeros(nl)
    zen    = zeros(nl)
    for l in eachindex(mask_urbanl)
        mask_urbanl[l] || continue
        coszen[l] = surfalb.coszen_col[lun.coli[l]]
        zen[l]    = acos(coszen[l])
    end

    # ---- Initialize output ----
    for ib in 1:numrad
        for c in eachindex(mask_urbanc)
            mask_urbanc[c] || continue
            albgrd[c, ib] = 0.0
            albgri[c, ib] = 0.0
        end

        for p in eachindex(mask_urbanp)
            mask_urbanp[p] || continue
            l = pch.landunit[p]
            c = pch.column[p]
            if col.itype[c] == ICOL_SUNWALL
                albd[p, ib] = vf_sw[l]
                albi[p, ib] = vf_sw[l]
            elseif col.itype[c] == ICOL_SHADEWALL
                albd[p, ib] = vf_sw[l]
                albi[p, ib] = vf_sw[l]
            elseif col.itype[c] == ICOL_ROAD_PERV || col.itype[c] == ICOL_ROAD_IMPERV
                albd[p, ib] = vf_sr[l]
                albi[p, ib] = vf_sr[l]
            elseif col.itype[c] == ICOL_ROOF
                albd[p, ib] = 1.0
                albi[p, ib] = 1.0
            end
            fabd[p, ib]     = 0.0
            fabd_sun[p, ib] = 0.0
            fabd_sha[p, ib] = 0.0
            fabi[p, ib]     = 0.0
            fabi_sun[p, ib] = 0.0
            fabi_sha[p, ib] = 0.0
            if coszen[l] > 0.0
                ftdd[p, ib] = 1.0
            else
                ftdd[p, ib] = 0.0
            end
            ftid[p, ib] = 0.0
            if coszen[l] > 0.0
                ftii[p, ib] = 1.0
            else
                ftii[p, ib] = 0.0
            end
        end
    end

    # ---- Initialize solarabs and sref defaults ----
    sabs_roof_dir      = solarabs.sabs_roof_dir_lun
    sabs_roof_dif      = solarabs.sabs_roof_dif_lun
    sabs_sunwall_dir   = solarabs.sabs_sunwall_dir_lun
    sabs_sunwall_dif   = solarabs.sabs_sunwall_dif_lun
    sabs_shadewall_dir = solarabs.sabs_shadewall_dir_lun
    sabs_shadewall_dif = solarabs.sabs_shadewall_dif_lun
    sabs_improad_dir   = solarabs.sabs_improad_dir_lun
    sabs_improad_dif   = solarabs.sabs_improad_dif_lun
    sabs_perroad_dir   = solarabs.sabs_perroad_dir_lun
    sabs_perroad_dif   = solarabs.sabs_perroad_dif_lun

    # Local arrays for sref (landunit × numrad)
    sref_roof_dir      = ones(nl, numrad)
    sref_roof_dif      = ones(nl, numrad)
    sref_sunwall_dir   = zeros(nl, numrad)
    sref_sunwall_dif   = zeros(nl, numrad)
    sref_shadewall_dir = zeros(nl, numrad)
    sref_shadewall_dif = zeros(nl, numrad)
    sref_improad_dir   = zeros(nl, numrad)
    sref_improad_dif   = zeros(nl, numrad)
    sref_perroad_dir   = zeros(nl, numrad)
    sref_perroad_dif   = zeros(nl, numrad)

    for ib in 1:numrad
        for l in eachindex(mask_urbanl)
            mask_urbanl[l] || continue
            sabs_roof_dir[l, ib]      = 0.0
            sabs_roof_dif[l, ib]      = 0.0
            sabs_sunwall_dir[l, ib]   = 0.0
            sabs_sunwall_dif[l, ib]   = 0.0
            sabs_shadewall_dir[l, ib] = 0.0
            sabs_shadewall_dif[l, ib] = 0.0
            sabs_improad_dir[l, ib]   = 0.0
            sabs_improad_dif[l, ib]   = 0.0
            sabs_perroad_dir[l, ib]   = 0.0
            sabs_perroad_dif[l, ib]   = 0.0
            sref_roof_dir[l, ib]      = 1.0
            sref_roof_dif[l, ib]      = 1.0
            sref_sunwall_dir[l, ib]   = vf_sw[l]
            sref_sunwall_dif[l, ib]   = vf_sw[l]
            sref_shadewall_dir[l, ib] = vf_sw[l]
            sref_shadewall_dif[l, ib] = vf_sw[l]
            sref_improad_dir[l, ib]   = vf_sr[l]
            sref_improad_dif[l, ib]   = vf_sr[l]
            sref_perroad_dir[l, ib]   = vf_sr[l]
            sref_perroad_dif[l, ib]   = vf_sr[l]
        end
    end

    # ---- Check if any landunits have sun above horizon ----
    num_solar = 0
    for l in eachindex(mask_urbanl)
        mask_urbanl[l] || continue
        if coszen[l] > 0.0
            num_solar += 1
        end
    end

    if num_solar > 0
        # Set constants - solar fluxes are per unit incoming flux
        sdir = ones(nl, numrad)
        sdif = ones(nl, numrad)

        # Local arrays for incident radiation
        sdir_road      = zeros(nl, numrad)
        sdif_road      = zeros(nl, numrad)
        sdir_sunwall   = zeros(nl, numrad)
        sdif_sunwall   = zeros(nl, numrad)
        sdir_shadewall = zeros(nl, numrad)
        sdif_shadewall = zeros(nl, numrad)

        # Snow albedos
        albsnd_roof    = zeros(nl, numrad)
        albsni_roof    = zeros(nl, numrad)
        albsnd_improad = zeros(nl, numrad)
        albsni_improad = zeros(nl, numrad)
        albsnd_perroad = zeros(nl, numrad)
        albsni_perroad = zeros(nl, numrad)

        # Combined albedos with snow
        alb_roof_dir_s    = zeros(nl, numrad)
        alb_roof_dif_s    = zeros(nl, numrad)
        alb_improad_dir_s = zeros(nl, numrad)
        alb_perroad_dir_s = zeros(nl, numrad)
        alb_improad_dif_s = zeros(nl, numrad)
        alb_perroad_dif_s = zeros(nl, numrad)

        # Compute total h2osno for urban columns
        # (inline version of CalculateTotalH2osno)
        nlevsno = varpar.nlevsno
        nc = length(mask_urbanc)
        h2osno_total = zeros(nc)
        for c in eachindex(mask_urbanc)
            mask_urbanc[c] || continue
            h2osno_total[c] = waterstatebulk.ws.h2osno_no_layers_col[c]
            if col.snl[c] < 0
                for j in (col.snl[c] + nlevsno):(nlevsno)
                    h2osno_total[c] += waterstatebulk.ws.h2osoi_ice_col[c, j] +
                                       waterstatebulk.ws.h2osoi_liq_col[c, j]
                end
            end
        end

        # Incident direct beam radiation
        incident_direct!(mask_urbanl, lun.canyon_hwr, coszen, zen,
                         sdir, sdir_road, sdir_sunwall, sdir_shadewall)

        # Incident diffuse radiation
        incident_diffuse!(mask_urbanl, lun.canyon_hwr, sdif,
                          sdif_road, sdif_sunwall, sdif_shadewall,
                          urbanparams)

        # Snow albedos - direct beam (ind=0)
        snow_albedo!(mask_urbanc, col, coszen, 0,
                     albsnd_roof, albsnd_improad, albsnd_perroad,
                     h2osno_total)

        # Snow albedos - diffuse (ind=1)
        snow_albedo!(mask_urbanc, col, coszen, 1,
                     albsni_roof, albsni_improad, albsni_perroad,
                     h2osno_total)

        # Combine snow-free and snow albedos
        alb_roof_dir    = urbanparams.alb_roof_dir
        alb_roof_dif    = urbanparams.alb_roof_dif
        alb_improad_dir = urbanparams.alb_improad_dir
        alb_improad_dif = urbanparams.alb_improad_dif
        alb_perroad_dir = urbanparams.alb_perroad_dir
        alb_perroad_dif = urbanparams.alb_perroad_dif
        alb_wall_dir    = urbanparams.alb_wall_dir
        alb_wall_dif    = urbanparams.alb_wall_dif

        for ib in 1:numrad
            for c in eachindex(mask_urbanc)
                mask_urbanc[c] || continue
                l = col.landunit[c]
                if col.itype[c] == ICOL_ROOF
                    alb_roof_dir_s[l, ib] = alb_roof_dir[l, ib] * (1.0 - frac_sno[c]) +
                        albsnd_roof[l, ib] * frac_sno[c]
                    alb_roof_dif_s[l, ib] = alb_roof_dif[l, ib] * (1.0 - frac_sno[c]) +
                        albsni_roof[l, ib] * frac_sno[c]
                elseif col.itype[c] == ICOL_ROAD_IMPERV
                    alb_improad_dir_s[l, ib] = alb_improad_dir[l, ib] * (1.0 - frac_sno[c]) +
                        albsnd_improad[l, ib] * frac_sno[c]
                    alb_improad_dif_s[l, ib] = alb_improad_dif[l, ib] * (1.0 - frac_sno[c]) +
                        albsni_improad[l, ib] * frac_sno[c]
                elseif col.itype[c] == ICOL_ROAD_PERV
                    alb_perroad_dir_s[l, ib] = alb_perroad_dir[l, ib] * (1.0 - frac_sno[c]) +
                        albsnd_perroad[l, ib] * frac_sno[c]
                    alb_perroad_dif_s[l, ib] = alb_perroad_dif[l, ib] * (1.0 - frac_sno[c]) +
                        albsni_perroad[l, ib] * frac_sno[c]
                end
            end
        end

        # Net solar radiation with multiple reflections
        net_solar!(mask_urbanl, coszen, lun.canyon_hwr, lun.wtroad_perv,
                   sdir, sdif,
                   alb_improad_dir_s, alb_perroad_dir_s, alb_wall_dir, alb_roof_dir_s,
                   alb_improad_dif_s, alb_perroad_dif_s, alb_wall_dif, alb_roof_dif_s,
                   sdir_road, sdir_sunwall, sdir_shadewall,
                   sdif_road, sdif_sunwall, sdif_shadewall,
                   sref_improad_dir, sref_perroad_dir, sref_sunwall_dir,
                   sref_shadewall_dir, sref_roof_dir,
                   sref_improad_dif, sref_perroad_dif, sref_sunwall_dif,
                   sref_shadewall_dif, sref_roof_dif,
                   urbanparams, solarabs)

        # ---- Map urban output to surfalb_inst components ----
        for ib in 1:numrad
            for c in eachindex(mask_urbanc)
                mask_urbanc[c] || continue
                l = col.landunit[c]
                if col.itype[c] == ICOL_ROOF
                    albgrd[c, ib] = sref_roof_dir[l, ib]
                    albgri[c, ib] = sref_roof_dif[l, ib]
                elseif col.itype[c] == ICOL_SUNWALL
                    albgrd[c, ib] = sref_sunwall_dir[l, ib]
                    albgri[c, ib] = sref_sunwall_dif[l, ib]
                elseif col.itype[c] == ICOL_SHADEWALL
                    albgrd[c, ib] = sref_shadewall_dir[l, ib]
                    albgri[c, ib] = sref_shadewall_dif[l, ib]
                elseif col.itype[c] == ICOL_ROAD_PERV
                    albgrd[c, ib] = sref_perroad_dir[l, ib]
                    albgri[c, ib] = sref_perroad_dif[l, ib]
                elseif col.itype[c] == ICOL_ROAD_IMPERV
                    albgrd[c, ib] = sref_improad_dir[l, ib]
                    albgri[c, ib] = sref_improad_dif[l, ib]
                end
                if coszen[l] > 0.0
                    albgrd_hst[c, ib] = albgrd[c, ib]
                    albgri_hst[c, ib] = albgri[c, ib]
                end
            end
            for p in eachindex(mask_urbanp)
                mask_urbanp[p] || continue
                c = pch.column[p]
                l = pch.landunit[p]
                albd[p, ib] = albgrd[c, ib]
                albi[p, ib] = albgri[c, ib]
                if coszen[l] > 0.0
                    albd_hst[p, ib] = albd[p, ib]
                    albi_hst[p, ib] = albi[p, ib]
                end
            end
        end
    end

    return nothing
end
