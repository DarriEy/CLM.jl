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
# Kernels (KernelAbstractions) for fully-independent per-element loops.
# --------------------------------------------------------------------------

# Urban snow albedo, one thread per column. Each column writes only the
# matrix matching its itype (roof / improad / perroad), so the three output
# matrices never alias across columns -> independent iterations.
@kernel function _ualb_snow_albedo_kernel!(albsn_roof, @Const(mask_urbanc),
                                           @Const(landunit), @Const(itype),
                                           @Const(coszen), ind::Int,
                                           albsn_improad, albsn_perroad,
                                           @Const(h2osno_total),
                                           icol_roof::Int, icol_imp::Int,
                                           icol_per::Int, snal0, snal1)
    c = @index(Global)
    @inbounds if mask_urbanc[c]
        l = landunit[c]
        if coszen[l] > 0.0 && h2osno_total[c] > 0.0
            if itype[c] == icol_roof
                albsn_roof[l, 1] = snal0
                albsn_roof[l, 2] = snal1
            elseif itype[c] == icol_imp
                albsn_improad[l, 1] = snal0
                albsn_improad[l, 2] = snal1
            elseif itype[c] == icol_per
                albsn_perroad[l, 1] = snal0
                albsn_perroad[l, 2] = snal1
            end
        else
            if itype[c] == icol_roof
                albsn_roof[l, 1] = 0.0
                albsn_roof[l, 2] = 0.0
            elseif itype[c] == icol_imp
                albsn_improad[l, 1] = 0.0
                albsn_improad[l, 2] = 0.0
            elseif itype[c] == icol_per
                albsn_perroad[l, 1] = 0.0
                albsn_perroad[l, 2] = 0.0
            end
        end
    end
end

ualb_snow_albedo!(albsn_roof, mask_urbanc, landunit, itype, coszen, ind::Int,
                  albsn_improad, albsn_perroad, h2osno_total,
                  icol_roof::Int, icol_imp::Int, icol_per::Int, snal0, snal1) =
    _launch!(_ualb_snow_albedo_kernel!, albsn_roof, mask_urbanc, landunit, itype,
             coszen, ind, albsn_improad, albsn_perroad, h2osno_total,
             icol_roof, icol_imp, icol_per, snal0, snal1;
             ndrange = length(mask_urbanc))

# Cosine solar zenith angle + zenith angle, one thread per landunit.
@kernel function _ualb_coszen_kernel!(coszen, @Const(mask_urbanl),
                                      @Const(coszen_col), @Const(coli), zen)
    l = @index(Global)
    @inbounds if mask_urbanl[l]
        cz = coszen_col[coli[l]]
        coszen[l] = cz
        zen[l] = acos(cz)
    end
end

ualb_coszen!(coszen, mask_urbanl, coszen_col, coli, zen) =
    _launch!(_ualb_coszen_kernel!, coszen, mask_urbanl, coszen_col, coli, zen;
             ndrange = length(mask_urbanl))

# Zero ground albedo outputs over urban columns (2D: column x waveband).
@kernel function _ualb_init_albgr_kernel!(albgrd, @Const(mask_urbanc), albgri)
    c, ib = @index(Global, NTuple)
    @inbounds if mask_urbanc[c]
        albgrd[c, ib] = 0.0
        albgri[c, ib] = 0.0
    end
end

ualb_init_albgr!(albgrd, mask_urbanc, albgri, numrad::Int) =
    _launch!(_ualb_init_albgr_kernel!, albgrd, mask_urbanc, albgri;
             ndrange = (length(mask_urbanc), numrad))

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
function snow_albedo!(mask_urbanc::AbstractVector{Bool}, col::ColumnData,
                      coszen::Vector{<:Real}, ind::Int,
                      albsn_roof::Matrix{<:Real},
                      albsn_improad::Matrix{<:Real},
                      albsn_perroad::Matrix{<:Real},
                      h2osno_total::Vector{<:Real})

    ualb_snow_albedo!(albsn_roof, mask_urbanc, col.landunit, col.itype,
                      coszen, ind, albsn_improad, albsn_perroad, h2osno_total,
                      ICOL_ROOF, ICOL_ROAD_IMPERV, ICOL_ROAD_PERV, SNAL0, SNAL1)

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
# Direct-beam incidence, one thread per landunit; loops the NUMRAD wavebands
# in-thread. theta0/tanzen (host per-landunit temporaries) become thread-local
# scalars. The Fortran balance-check `error()` (an analytic self-consistency
# identity that never trips for physical input) is dropped — it cannot run on a
# GPU and the outputs are byte-identical without it.
@kernel function _incident_direct_kernel!(sdir_road, @Const(mask_urbanl),
        @Const(canyon_hwr), @Const(coszen), @Const(zen), @Const(sdir),
        sdir_sunwall, sdir_shadewall, numrad::Int)
    l = @index(Global)
    @inbounds if mask_urbanl[l]
        T = eltype(sdir_road)
        if coszen[l] > zero(T)
            rpi = T(RPI)
            theta0 = asin(min(one(T) / (canyon_hwr[l] * tan(max(zen[l], T(0.000001)))), one(T)))
            tanzen = tan(zen[l])
            for ib in 1:numrad
                sdir_shadewall[l, ib] = zero(T)
                sdir_road[l, ib] = sdir[l, ib] *
                    (T(2) * theta0 / rpi - T(2) / rpi * canyon_hwr[l] * tanzen * (one(T) - cos(theta0)))
                sdir_sunwall[l, ib] = T(2) * sdir[l, ib] * ((one(T) / canyon_hwr[l]) *
                    (T(0.5) - theta0 / rpi) + (one(T) / rpi) * tanzen * (one(T) - cos(theta0)))
            end
        else
            for ib in 1:numrad
                sdir_road[l, ib] = zero(T)
                sdir_sunwall[l, ib] = zero(T)
                sdir_shadewall[l, ib] = zero(T)
            end
        end
    end
end

function incident_direct!(mask_urbanl::AbstractVector{Bool},
                           canyon_hwr::AbstractVector{<:Real},
                           coszen::AbstractVector{<:Real},
                           zen::AbstractVector{<:Real},
                           sdir::AbstractMatrix{<:Real},
                           sdir_road::AbstractMatrix{<:Real},
                           sdir_sunwall::AbstractMatrix{<:Real},
                           sdir_shadewall::AbstractMatrix{<:Real})

    _launch!(_incident_direct_kernel!, sdir_road, mask_urbanl, canyon_hwr, coszen,
             zen, sdir, sdir_sunwall, sdir_shadewall, NUMRAD;
             ndrange = length(mask_urbanl))
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
# Diffuse incidence, one thread per landunit (loops NUMRAD in-thread). Analytic
# balance-check `error()` dropped (see incident_direct!); outputs byte-identical.
@kernel function _incident_diffuse_kernel!(sdif_road, @Const(mask_urbanl),
        @Const(sdif), @Const(vf_sr), @Const(vf_sw),
        sdif_sunwall, sdif_shadewall, numrad::Int)
    l = @index(Global)
    @inbounds if mask_urbanl[l]
        for ib in 1:numrad
            sdif_road[l, ib]      = sdif[l, ib] * vf_sr[l]
            sdif_sunwall[l, ib]   = sdif[l, ib] * vf_sw[l]
            sdif_shadewall[l, ib] = sdif[l, ib] * vf_sw[l]
        end
    end
end

function incident_diffuse!(mask_urbanl::AbstractVector{Bool},
                            canyon_hwr::AbstractVector{<:Real},
                            sdif::AbstractMatrix{<:Real},
                            sdif_road::AbstractMatrix{<:Real},
                            sdif_sunwall::AbstractMatrix{<:Real},
                            sdif_shadewall::AbstractMatrix{<:Real},
                            urbanparams::UrbanParamsData)

    _launch!(_incident_diffuse_kernel!, sdif_road, mask_urbanl, sdif,
             urbanparams.vf_sr, urbanparams.vf_sw, sdif_sunwall, sdif_shadewall,
             NUMRAD; ndrange = length(mask_urbanl))
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
# --------------------------------------------------------------------------
# Device-view bundles (30+ matrix args → 4 grouped structs to stay under
# Metal's ~31-buffer launch limit). Field names mirror the kernel locals.
# --------------------------------------------------------------------------
struct _NSIn{M}
    sdir::M; sdif::M
    alb_improad_dir::M; alb_perroad_dir::M; alb_wall_dir::M; alb_roof_dir::M
    alb_improad_dif::M; alb_perroad_dif::M; alb_wall_dif::M; alb_roof_dif::M
    sdir_road::M; sdir_sunwall::M; sdir_shadewall::M
    sdif_road::M; sdif_sunwall::M; sdif_shadewall::M
end
Adapt.@adapt_structure _NSIn

struct _NSSref{M}
    improad_dir::M; perroad_dir::M; sunwall_dir::M; shadewall_dir::M; roof_dir::M
    improad_dif::M; perroad_dif::M; sunwall_dif::M; shadewall_dif::M; roof_dif::M
end
Adapt.@adapt_structure _NSSref

struct _NSSabs{M}
    roof_dir::M; roof_dif::M; sunwall_dir::M; sunwall_dif::M
    shadewall_dir::M; shadewall_dif::M; improad_dir::M; improad_dif::M
    perroad_dir::M; perroad_dif::M
end
Adapt.@adapt_structure _NSSabs

struct _NSLun{V}
    coszen::V; canyon_hwr::V; wtroad_perv::V
    vf_sr::V; vf_wr::V; vf_sw::V; vf_rw::V; vf_ww::V
end
Adapt.@adapt_structure _NSLun

# net_solar! kernel — ONE thread per landunit, loops NUMRAD in-thread. The ~40
# per-landunit host working arrays collapse to thread-local scalars (each thread
# owns one landunit) → byte-identical arithmetic to the host loop. The convergence
# `break` (step 5) is preserved; only the non-convergence / balance `error()` calls
# are dropped (they cannot run on a GPU and never fire for physical input, so
# outputs are byte-identical). Literals carried at the working type T.
@kernel function _net_solar_kernel!(@Const(mask_urbanl), lun, inp, sref, sabs, numrad::Int)
    l = @index(Global)
    @inbounds if mask_urbanl[l]
        T = eltype(inp.sdir)
        coszen = lun.coszen
        if coszen[l] > zero(T)
            n = 50
            errcrit = T(0.00001)
            chwr = lun.canyon_hwr[l]
            wtperv = lun.wtroad_perv[l]
            wtimp = one(T) - wtperv
            vsr = lun.vf_sr[l]; vwr = lun.vf_wr[l]; vsw = lun.vf_sw[l]
            vrw = lun.vf_rw[l]; vww = lun.vf_ww[l]

            for ib in 1:numrad
                aim_dir = inp.alb_improad_dir[l, ib]; ape_dir = inp.alb_perroad_dir[l, ib]
                awl_dir = inp.alb_wall_dir[l, ib]
                aim_dif = inp.alb_improad_dif[l, ib]; ape_dif = inp.alb_perroad_dif[l, ib]
                awl_dif = inp.alb_wall_dif[l, ib]
                sr_road = inp.sdir_road[l, ib]; sr_sun = inp.sdir_sunwall[l, ib]; sr_shd = inp.sdir_shadewall[l, ib]
                sf_road = inp.sdif_road[l, ib]; sf_sun = inp.sdif_sunwall[l, ib]; sf_shd = inp.sdif_shadewall[l, ib]

                # ---- Initial absorption and reflection (direct) ----
                road_a_dir = zero(T); road_r_dir = zero(T)
                improad_a_dir = (one(T) - aim_dir) * sr_road
                improad_r_dir =            aim_dir  * sr_road
                improad_r_sky_dir       = improad_r_dir * vsr
                improad_r_sunwall_dir   = improad_r_dir * vwr
                improad_r_shadewall_dir = improad_r_dir * vwr
                road_a_dir += improad_a_dir * wtimp
                road_r_dir += improad_r_dir * wtimp

                perroad_a_dir = (one(T) - ape_dir) * sr_road
                perroad_r_dir =            ape_dir  * sr_road
                perroad_r_sky_dir       = perroad_r_dir * vsr
                perroad_r_sunwall_dir   = perroad_r_dir * vwr
                perroad_r_shadewall_dir = perroad_r_dir * vwr
                road_a_dir += perroad_a_dir * wtperv
                road_r_dir += perroad_r_dir * wtperv

                road_r_sky_dir       = road_r_dir * vsr
                road_r_sunwall_dir   = road_r_dir * vwr
                road_r_shadewall_dir = road_r_dir * vwr

                sunwall_a_dir = (one(T) - awl_dir) * sr_sun
                sunwall_r_dir =            awl_dir  * sr_sun
                sunwall_r_sky_dir       = sunwall_r_dir * vsw
                sunwall_r_road_dir      = sunwall_r_dir * vrw
                sunwall_r_shadewall_dir = sunwall_r_dir * vww

                shadewall_a_dir = (one(T) - awl_dir) * sr_shd
                shadewall_r_dir =            awl_dir  * sr_shd
                shadewall_r_sky_dir     = shadewall_r_dir * vsw
                shadewall_r_road_dir    = shadewall_r_dir * vrw
                shadewall_r_sunwall_dir = shadewall_r_dir * vww

                # ---- Initial absorption and reflection (diffuse) ----
                road_a_dif = zero(T); road_r_dif = zero(T)
                improad_a_dif = (one(T) - aim_dif) * sf_road
                improad_r_dif =            aim_dif  * sf_road
                improad_r_sky_dif       = improad_r_dif * vsr
                improad_r_sunwall_dif   = improad_r_dif * vwr
                improad_r_shadewall_dif = improad_r_dif * vwr
                road_a_dif += improad_a_dif * wtimp
                road_r_dif += improad_r_dif * wtimp

                perroad_a_dif = (one(T) - ape_dif) * sf_road
                perroad_r_dif =            ape_dif  * sf_road
                perroad_r_sky_dif       = perroad_r_dif * vsr
                perroad_r_sunwall_dif   = perroad_r_dif * vwr
                perroad_r_shadewall_dif = perroad_r_dif * vwr
                road_a_dif += perroad_a_dif * wtperv
                road_r_dif += perroad_r_dif * wtperv

                road_r_sky_dif       = road_r_dif * vsr
                road_r_sunwall_dif   = road_r_dif * vwr
                road_r_shadewall_dif = road_r_dif * vwr

                sunwall_a_dif = (one(T) - awl_dif) * sf_sun
                sunwall_r_dif =            awl_dif  * sf_sun
                sunwall_r_sky_dif       = sunwall_r_dif * vsw
                sunwall_r_road_dif      = sunwall_r_dif * vrw
                sunwall_r_shadewall_dif = sunwall_r_dif * vww

                shadewall_a_dif = (one(T) - awl_dif) * sf_shd
                shadewall_r_dif =            awl_dif  * sf_shd
                shadewall_r_sky_dif     = shadewall_r_dif * vsw
                shadewall_r_road_dif    = shadewall_r_dif * vrw
                shadewall_r_sunwall_dif = shadewall_r_dif * vww

                # Running absorption / reflection sums (local accumulators)
                sabs_improad_dir = improad_a_dir; sabs_perroad_dir = perroad_a_dir
                sabs_sunwall_dir = sunwall_a_dir; sabs_shadewall_dir = shadewall_a_dir
                sabs_improad_dif = improad_a_dif; sabs_perroad_dif = perroad_a_dif
                sabs_sunwall_dif = sunwall_a_dif; sabs_shadewall_dif = shadewall_a_dif
                sref_improad_dir = improad_r_sky_dir; sref_perroad_dir = perroad_r_sky_dir
                sref_sunwall_dir = sunwall_r_sky_dir; sref_shadewall_dir = shadewall_r_sky_dir
                sref_improad_dif = improad_r_sky_dif; sref_perroad_dif = perroad_r_sky_dif
                sref_sunwall_dif = sunwall_r_sky_dif; sref_shadewall_dif = shadewall_r_sky_dif

                # ---- Multiple reflections (direct) ----
                for iter in 1:n
                    stot = (sunwall_r_road_dir + shadewall_r_road_dir) * chwr
                    road_a_dir = zero(T); road_r_dir = zero(T)
                    improad_a_dir = (one(T) - aim_dir) * stot
                    improad_r_dir =            aim_dir  * stot
                    road_a_dir += improad_a_dir * wtimp
                    road_r_dir += improad_r_dir * wtimp
                    perroad_a_dir = (one(T) - ape_dir) * stot
                    perroad_r_dir =            ape_dir  * stot
                    road_a_dir += perroad_a_dir * wtperv
                    road_r_dir += perroad_r_dir * wtperv

                    stot = road_r_sunwall_dir / chwr + shadewall_r_sunwall_dir
                    sunwall_a_dir = (one(T) - awl_dir) * stot
                    sunwall_r_dir =            awl_dir  * stot

                    stot = road_r_shadewall_dir / chwr + sunwall_r_shadewall_dir
                    shadewall_a_dir = (one(T) - awl_dir) * stot
                    shadewall_r_dir =            awl_dir  * stot

                    sabs_improad_dir   += improad_a_dir
                    sabs_perroad_dir   += perroad_a_dir
                    sabs_sunwall_dir   += sunwall_a_dir
                    sabs_shadewall_dir += shadewall_a_dir

                    improad_r_sky_dir       = improad_r_dir * vsr
                    improad_r_sunwall_dir   = improad_r_dir * vwr
                    improad_r_shadewall_dir = improad_r_dir * vwr
                    perroad_r_sky_dir       = perroad_r_dir * vsr
                    perroad_r_sunwall_dir   = perroad_r_dir * vwr
                    perroad_r_shadewall_dir = perroad_r_dir * vwr
                    road_r_sky_dir          = road_r_dir * vsr
                    road_r_sunwall_dir      = road_r_dir * vwr
                    road_r_shadewall_dir    = road_r_dir * vwr
                    sunwall_r_sky_dir       = sunwall_r_dir * vsw
                    sunwall_r_road_dir      = sunwall_r_dir * vrw
                    sunwall_r_shadewall_dir = sunwall_r_dir * vww
                    shadewall_r_sky_dir     = shadewall_r_dir * vsw
                    shadewall_r_road_dir    = shadewall_r_dir * vrw
                    shadewall_r_sunwall_dir = shadewall_r_dir * vww

                    sref_improad_dir   += improad_r_sky_dir
                    sref_perroad_dir   += perroad_r_sky_dir
                    sref_sunwall_dir   += sunwall_r_sky_dir
                    sref_shadewall_dir += shadewall_r_sky_dir

                    crit = max(road_a_dir, sunwall_a_dir, shadewall_a_dir)
                    if crit < errcrit
                        break
                    end
                end

                # ---- Multiple reflections (diffuse) ----
                for iter in 1:n
                    stot = (sunwall_r_road_dif + shadewall_r_road_dif) * chwr
                    road_a_dif = zero(T); road_r_dif = zero(T)
                    improad_a_dif = (one(T) - aim_dif) * stot
                    improad_r_dif =            aim_dif  * stot
                    road_a_dif += improad_a_dif * wtimp
                    road_r_dif += improad_r_dif * wtimp
                    perroad_a_dif = (one(T) - ape_dif) * stot
                    perroad_r_dif =            ape_dif  * stot
                    road_a_dif += perroad_a_dif * wtperv
                    road_r_dif += perroad_r_dif * wtperv

                    stot = road_r_sunwall_dif / chwr + shadewall_r_sunwall_dif
                    sunwall_a_dif = (one(T) - awl_dif) * stot
                    sunwall_r_dif =            awl_dif  * stot

                    stot = road_r_shadewall_dif / chwr + sunwall_r_shadewall_dif
                    shadewall_a_dif = (one(T) - awl_dif) * stot
                    shadewall_r_dif =            awl_dif  * stot

                    sabs_improad_dif   += improad_a_dif
                    sabs_perroad_dif   += perroad_a_dif
                    sabs_sunwall_dif   += sunwall_a_dif
                    sabs_shadewall_dif += shadewall_a_dif

                    improad_r_sky_dif       = improad_r_dif * vsr
                    improad_r_sunwall_dif   = improad_r_dif * vwr
                    improad_r_shadewall_dif = improad_r_dif * vwr
                    perroad_r_sky_dif       = perroad_r_dif * vsr
                    perroad_r_sunwall_dif   = perroad_r_dif * vwr
                    perroad_r_shadewall_dif = perroad_r_dif * vwr
                    road_r_sky_dif          = road_r_dif * vsr
                    road_r_sunwall_dif      = road_r_dif * vwr
                    road_r_shadewall_dif    = road_r_dif * vwr
                    sunwall_r_sky_dif       = sunwall_r_dif * vsw
                    sunwall_r_road_dif      = sunwall_r_dif * vrw
                    sunwall_r_shadewall_dif = sunwall_r_dif * vww
                    shadewall_r_sky_dif     = shadewall_r_dif * vsw
                    shadewall_r_road_dif    = shadewall_r_dif * vrw
                    shadewall_r_sunwall_dif = shadewall_r_dif * vww

                    sref_improad_dif   += improad_r_sky_dif
                    sref_perroad_dif   += perroad_r_sky_dif
                    sref_sunwall_dif   += sunwall_r_sky_dif
                    sref_shadewall_dif += shadewall_r_sky_dif

                    crit = max(road_a_dif, sunwall_a_dif, shadewall_a_dif)
                    if crit < errcrit
                        break
                    end
                end

                # ---- Store canyon sums to the output arrays ----
                sabs.improad_dir[l, ib]   = sabs_improad_dir
                sabs.perroad_dir[l, ib]   = sabs_perroad_dir
                sabs.sunwall_dir[l, ib]   = sabs_sunwall_dir
                sabs.shadewall_dir[l, ib] = sabs_shadewall_dir
                sabs.improad_dif[l, ib]   = sabs_improad_dif
                sabs.perroad_dif[l, ib]   = sabs_perroad_dif
                sabs.sunwall_dif[l, ib]   = sabs_sunwall_dif
                sabs.shadewall_dif[l, ib] = sabs_shadewall_dif

                sref.improad_dir[l, ib]   = sref_improad_dir
                sref.perroad_dir[l, ib]   = sref_perroad_dir
                sref.sunwall_dir[l, ib]   = sref_sunwall_dir
                sref.shadewall_dir[l, ib] = sref_shadewall_dir
                sref.improad_dif[l, ib]   = sref_improad_dif
                sref.perroad_dif[l, ib]   = sref_perroad_dif
                sref.sunwall_dif[l, ib]   = sref_sunwall_dif
                sref.shadewall_dif[l, ib] = sref_shadewall_dif

                # Roof reflected + absorbed
                sref.roof_dir[l, ib] = inp.alb_roof_dir[l, ib] * inp.sdir[l, ib]
                sref.roof_dif[l, ib] = inp.alb_roof_dif[l, ib] * inp.sdif[l, ib]
                sabs.roof_dir[l, ib] = inp.sdir[l, ib] - sref.roof_dir[l, ib]
                sabs.roof_dif[l, ib] = inp.sdif[l, ib] - sref.roof_dif[l, ib]
            end  # ib
        end
    end
end

function net_solar!(mask_urbanl::AbstractVector{Bool},
                    coszen::AbstractVector{<:Real},
                    canyon_hwr::AbstractVector{<:Real},
                    wtroad_perv::AbstractVector{<:Real},
                    sdir::AbstractMatrix{<:Real},
                    sdif::AbstractMatrix{<:Real},
                    alb_improad_dir::AbstractMatrix{<:Real},
                    alb_perroad_dir::AbstractMatrix{<:Real},
                    alb_wall_dir::AbstractMatrix{<:Real},
                    alb_roof_dir::AbstractMatrix{<:Real},
                    alb_improad_dif::AbstractMatrix{<:Real},
                    alb_perroad_dif::AbstractMatrix{<:Real},
                    alb_wall_dif::AbstractMatrix{<:Real},
                    alb_roof_dif::AbstractMatrix{<:Real},
                    sdir_road::AbstractMatrix{<:Real},
                    sdir_sunwall::AbstractMatrix{<:Real},
                    sdir_shadewall::AbstractMatrix{<:Real},
                    sdif_road::AbstractMatrix{<:Real},
                    sdif_sunwall::AbstractMatrix{<:Real},
                    sdif_shadewall::AbstractMatrix{<:Real},
                    sref_improad_dir::AbstractMatrix{<:Real},
                    sref_perroad_dir::AbstractMatrix{<:Real},
                    sref_sunwall_dir::AbstractMatrix{<:Real},
                    sref_shadewall_dir::AbstractMatrix{<:Real},
                    sref_roof_dir::AbstractMatrix{<:Real},
                    sref_improad_dif::AbstractMatrix{<:Real},
                    sref_perroad_dif::AbstractMatrix{<:Real},
                    sref_sunwall_dif::AbstractMatrix{<:Real},
                    sref_shadewall_dif::AbstractMatrix{<:Real},
                    sref_roof_dif::AbstractMatrix{<:Real},
                    urbanparams::UrbanParamsData,
                    solarabs::SolarAbsorbedData)

    inp = _NSIn(sdir, sdif, alb_improad_dir, alb_perroad_dir, alb_wall_dir, alb_roof_dir,
                alb_improad_dif, alb_perroad_dif, alb_wall_dif, alb_roof_dif,
                sdir_road, sdir_sunwall, sdir_shadewall, sdif_road, sdif_sunwall, sdif_shadewall)
    sref = _NSSref(sref_improad_dir, sref_perroad_dir, sref_sunwall_dir, sref_shadewall_dir,
                   sref_roof_dir, sref_improad_dif, sref_perroad_dif, sref_sunwall_dif,
                   sref_shadewall_dif, sref_roof_dif)
    sabs = _NSSabs(solarabs.sabs_roof_dir_lun, solarabs.sabs_roof_dif_lun,
                   solarabs.sabs_sunwall_dir_lun, solarabs.sabs_sunwall_dif_lun,
                   solarabs.sabs_shadewall_dir_lun, solarabs.sabs_shadewall_dif_lun,
                   solarabs.sabs_improad_dir_lun, solarabs.sabs_improad_dif_lun,
                   solarabs.sabs_perroad_dir_lun, solarabs.sabs_perroad_dif_lun)
    lun = _NSLun(coszen, canyon_hwr, wtroad_perv, urbanparams.vf_sr, urbanparams.vf_wr,
                 urbanparams.vf_sw, urbanparams.vf_rw, urbanparams.vf_ww)

    _launch!(_net_solar_kernel!, mask_urbanl, lun, inp, sref, sabs, NUMRAD;
             ndrange = length(mask_urbanl))
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
function urban_albedo!(mask_urbanl::AbstractVector{Bool},
                       mask_urbanc::AbstractVector{Bool},
                       mask_urbanp::AbstractVector{Bool},
                       lun::LandunitData,
                       col::ColumnData,
                       pch::PatchData,
                       waterstatebulk::WaterStateBulkData,
                       waterdiagnosticbulk::WaterDiagnosticBulkData,
                       urbanparams::UrbanParamsData,
                       solarabs::SolarAbsorbedData,
                       surfalb::SurfaceAlbedoData)

    # No urban landunits/columns/patches → nothing to do. All work here is
    # urban-gated; the host loops would otherwise scalar-index a device mask even
    # with no urban present. Device-safe any() reduction; byte-identical (all-false loops).
    (any(mask_urbanl) || any(mask_urbanc) || any(mask_urbanp)) || return nothing

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
    FT = eltype(surfalb.coszen_col)
    coszen = zeros(FT, nl)
    zen    = zeros(FT, nl)
    ualb_coszen!(coszen, mask_urbanl, surfalb.coszen_col, lun.coli, zen)

    # ---- Initialize output ----
    ualb_init_albgr!(albgrd, mask_urbanc, albgri, numrad)
    for ib in 1:numrad
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
    sref_roof_dir      = ones(FT, nl, numrad)
    sref_roof_dif      = ones(FT, nl, numrad)
    sref_sunwall_dir   = zeros(FT, nl, numrad)
    sref_sunwall_dif   = zeros(FT, nl, numrad)
    sref_shadewall_dir = zeros(FT, nl, numrad)
    sref_shadewall_dif = zeros(FT, nl, numrad)
    sref_improad_dir   = zeros(FT, nl, numrad)
    sref_improad_dif   = zeros(FT, nl, numrad)
    sref_perroad_dir   = zeros(FT, nl, numrad)
    sref_perroad_dif   = zeros(FT, nl, numrad)

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
        sdir = ones(FT, nl, numrad)
        sdif = ones(FT, nl, numrad)

        # Local arrays for incident radiation
        sdir_road      = zeros(FT, nl, numrad)
        sdif_road      = zeros(FT, nl, numrad)
        sdir_sunwall   = zeros(FT, nl, numrad)
        sdif_sunwall   = zeros(FT, nl, numrad)
        sdir_shadewall = zeros(FT, nl, numrad)
        sdif_shadewall = zeros(FT, nl, numrad)

        # Snow albedos
        albsnd_roof    = zeros(FT, nl, numrad)
        albsni_roof    = zeros(FT, nl, numrad)
        albsnd_improad = zeros(FT, nl, numrad)
        albsni_improad = zeros(FT, nl, numrad)
        albsnd_perroad = zeros(FT, nl, numrad)
        albsni_perroad = zeros(FT, nl, numrad)

        # Combined albedos with snow
        alb_roof_dir_s    = zeros(FT, nl, numrad)
        alb_roof_dif_s    = zeros(FT, nl, numrad)
        alb_improad_dir_s = zeros(FT, nl, numrad)
        alb_perroad_dir_s = zeros(FT, nl, numrad)
        alb_improad_dif_s = zeros(FT, nl, numrad)
        alb_perroad_dif_s = zeros(FT, nl, numrad)

        # Compute total h2osno for urban columns
        # (inline version of CalculateTotalH2osno)
        nlevsno = varpar.nlevsno
        nc = length(mask_urbanc)
        h2osno_total = zeros(FT, nc)
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
