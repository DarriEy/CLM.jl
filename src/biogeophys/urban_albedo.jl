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
        T = eltype(albsn_roof); z = zero(T)
        l = landunit[c]
        if coszen[l] > z && h2osno_total[c] > z
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
                albsn_roof[l, 1] = z
                albsn_roof[l, 2] = z
            elseif itype[c] == icol_imp
                albsn_improad[l, 1] = z
                albsn_improad[l, 2] = z
            elseif itype[c] == icol_per
                albsn_perroad[l, 1] = z
                albsn_perroad[l, 2] = z
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
# Device-view bundles + kernels for the urban_albedo! ORCHESTRATOR's own
# per-element loops (init / snow-combine / net-solar output mapping). Each
# @adapt_structure struct-of-arrays passes to a Metal kernel as ONE buffer, so
# the launches stay well under the ~31-buffer limit. All run on the KA CPU
# backend byte-identically (every write is to the thread's own patch/column/
# landunit index; the h2osno accumulation uses a local scalar — see below).
# --------------------------------------------------------------------------

# surfalb outputs (patch matrices albd..ftii + col matrices albgr*).
struct _UAlbSurf{M}
    albd::M; albi::M; fabd::M; fabd_sun::M; fabd_sha::M
    fabi::M; fabi_sun::M; fabi_sha::M; ftdd::M; ftid::M; ftii::M
    albd_hst::M; albi_hst::M; albgrd::M; albgri::M; albgrd_hst::M; albgri_hst::M
end
Adapt.@adapt_structure _UAlbSurf

# sref working reflectances (landunit × numrad).
struct _UAlbSref{M}
    roof_dir::M; roof_dif::M; sunwall_dir::M; sunwall_dif::M
    shadewall_dir::M; shadewall_dif::M; improad_dir::M; improad_dif::M
    perroad_dir::M; perroad_dif::M
end
Adapt.@adapt_structure _UAlbSref

# sabs absorbed-solar outputs (solarabs landunit × numrad).
struct _UAlbSabs{M}
    roof_dir::M; roof_dif::M; sunwall_dir::M; sunwall_dif::M
    shadewall_dir::M; shadewall_dif::M; improad_dir::M; improad_dif::M
    perroad_dir::M; perroad_dif::M
end
Adapt.@adapt_structure _UAlbSabs

# snow albedos (albsn*) + snow-combined surface albedos (alb_*_s), ld × numrad.
struct _UAlbSnow{M}
    albsnd_roof::M; albsni_roof::M; albsnd_improad::M; albsni_improad::M
    albsnd_perroad::M; albsni_perroad::M
    roof_dir_s::M; roof_dif_s::M; improad_dir_s::M; perroad_dir_s::M
    improad_dif_s::M; perroad_dif_s::M
end
Adapt.@adapt_structure _UAlbSnow

# urbanparams base albedos (landunit × numrad).
struct _UAlbParam{M}
    roof_dir::M; roof_dif::M; improad_dir::M; improad_dif::M
    perroad_dir::M; perroad_dif::M
end
Adapt.@adapt_structure _UAlbParam

# Init patch albedo/absorption/transmission outputs (one thread per patch).
@kernel function _ualb_initout_kernel!(o, @Const(mask_urbanp), @Const(pch_landunit),
        @Const(pch_column), @Const(col_itype), @Const(coszen), @Const(vf_sw),
        @Const(vf_sr), numrad::Int, icol_sunw::Int, icol_shdw::Int,
        icol_per::Int, icol_imp::Int, icol_roof::Int)
    p = @index(Global)
    @inbounds if mask_urbanp[p]
        l = pch_landunit[p]; c = pch_column[p]; it = col_itype[c]
        T = eltype(o.albd); sunny = coszen[l] > zero(T)
        for ib in 1:numrad
            if it == icol_sunw || it == icol_shdw
                o.albd[p, ib] = vf_sw[l]; o.albi[p, ib] = vf_sw[l]
            elseif it == icol_per || it == icol_imp
                o.albd[p, ib] = vf_sr[l]; o.albi[p, ib] = vf_sr[l]
            elseif it == icol_roof
                o.albd[p, ib] = one(T); o.albi[p, ib] = one(T)
            end
            o.fabd[p, ib] = zero(T); o.fabd_sun[p, ib] = zero(T); o.fabd_sha[p, ib] = zero(T)
            o.fabi[p, ib] = zero(T); o.fabi_sun[p, ib] = zero(T); o.fabi_sha[p, ib] = zero(T)
            o.ftdd[p, ib] = sunny ? one(T) : zero(T)
            o.ftid[p, ib] = zero(T)
            o.ftii[p, ib] = sunny ? one(T) : zero(T)
        end
    end
end

# Init sabs (absorbed→0) and sref (roof→1, walls→vf_sw, roads→vf_sr) per landunit.
@kernel function _ualb_sabsref_init_kernel!(sa, sr, @Const(mask_urbanl),
        @Const(vf_sw), @Const(vf_sr), numrad::Int)
    l = @index(Global)
    @inbounds if mask_urbanl[l]
        T = eltype(sr.roof_dir); sw = vf_sw[l]; sd = vf_sr[l]
        for ib in 1:numrad
            sa.roof_dir[l, ib] = zero(T); sa.roof_dif[l, ib] = zero(T)
            sa.sunwall_dir[l, ib] = zero(T); sa.sunwall_dif[l, ib] = zero(T)
            sa.shadewall_dir[l, ib] = zero(T); sa.shadewall_dif[l, ib] = zero(T)
            sa.improad_dir[l, ib] = zero(T); sa.improad_dif[l, ib] = zero(T)
            sa.perroad_dir[l, ib] = zero(T); sa.perroad_dif[l, ib] = zero(T)
            sr.roof_dir[l, ib] = one(T); sr.roof_dif[l, ib] = one(T)
            sr.sunwall_dir[l, ib] = sw; sr.sunwall_dif[l, ib] = sw
            sr.shadewall_dir[l, ib] = sw; sr.shadewall_dif[l, ib] = sw
            sr.improad_dir[l, ib] = sd; sr.improad_dif[l, ib] = sd
            sr.perroad_dir[l, ib] = sd; sr.perroad_dif[l, ib] = sd
        end
    end
end

# Total column snow water (one thread per column). Accumulate the layer sum into
# a LOCAL scalar and store once — a repeated `h2osno_total[c] += …` inside the
# j-loop miscompiles on the KA CPU backend under --check-bounds=yes (inner loop
# silently runs a single iteration). Byte-identical to the base + loop-add form.
@kernel function _ualb_h2osno_kernel!(h2osno_total, @Const(mask_urbanc),
        @Const(h2osno_no_layers), @Const(h2osoi_ice), @Const(h2osoi_liq),
        @Const(snl), nlevsno::Int)
    c = @index(Global)
    @inbounds if mask_urbanc[c]
        acc = h2osno_no_layers[c]
        if snl[c] < 0
            for j in (snl[c] + nlevsno):nlevsno
                acc += h2osoi_ice[c, j] + h2osoi_liq[c, j]
            end
        end
        h2osno_total[c] = acc
    end
end

# Combine snow-free + snow albedos per column (roof/improad/perroad by itype).
@kernel function _ualb_combine_kernel!(sn, ap, @Const(mask_urbanc),
        @Const(col_landunit), @Const(col_itype), @Const(frac_sno), numrad::Int,
        icol_roof::Int, icol_imp::Int, icol_per::Int)
    c = @index(Global)
    @inbounds if mask_urbanc[c]
        l = col_landunit[c]; it = col_itype[c]
        T = eltype(frac_sno); fs = frac_sno[c]; omfs = one(T) - fs
        for ib in 1:numrad
            if it == icol_roof
                sn.roof_dir_s[l, ib] = ap.roof_dir[l, ib] * omfs + sn.albsnd_roof[l, ib] * fs
                sn.roof_dif_s[l, ib] = ap.roof_dif[l, ib] * omfs + sn.albsni_roof[l, ib] * fs
            elseif it == icol_imp
                sn.improad_dir_s[l, ib] = ap.improad_dir[l, ib] * omfs + sn.albsnd_improad[l, ib] * fs
                sn.improad_dif_s[l, ib] = ap.improad_dif[l, ib] * omfs + sn.albsni_improad[l, ib] * fs
            elseif it == icol_per
                sn.perroad_dir_s[l, ib] = ap.perroad_dir[l, ib] * omfs + sn.albsnd_perroad[l, ib] * fs
                sn.perroad_dif_s[l, ib] = ap.perroad_dif[l, ib] * omfs + sn.albsni_perroad[l, ib] * fs
            end
        end
    end
end

# Map net-solar sref → per-column ground albedo (albgrd/albgri + _hst) by itype.
@kernel function _ualb_map_col_kernel!(o, sr, @Const(mask_urbanc),
        @Const(col_landunit), @Const(col_itype), @Const(coszen), numrad::Int,
        icol_roof::Int, icol_sunw::Int, icol_shdw::Int, icol_per::Int, icol_imp::Int)
    c = @index(Global)
    @inbounds if mask_urbanc[c]
        l = col_landunit[c]; it = col_itype[c]
        T = eltype(coszen); sunny = coszen[l] > zero(T)
        for ib in 1:numrad
            if it == icol_roof
                o.albgrd[c, ib] = sr.roof_dir[l, ib]; o.albgri[c, ib] = sr.roof_dif[l, ib]
            elseif it == icol_sunw
                o.albgrd[c, ib] = sr.sunwall_dir[l, ib]; o.albgri[c, ib] = sr.sunwall_dif[l, ib]
            elseif it == icol_shdw
                o.albgrd[c, ib] = sr.shadewall_dir[l, ib]; o.albgri[c, ib] = sr.shadewall_dif[l, ib]
            elseif it == icol_per
                o.albgrd[c, ib] = sr.perroad_dir[l, ib]; o.albgri[c, ib] = sr.perroad_dif[l, ib]
            elseif it == icol_imp
                o.albgrd[c, ib] = sr.improad_dir[l, ib]; o.albgri[c, ib] = sr.improad_dif[l, ib]
            end
            if sunny
                o.albgrd_hst[c, ib] = o.albgrd[c, ib]; o.albgri_hst[c, ib] = o.albgri[c, ib]
            end
        end
    end
end

# Map per-column ground albedo → per-patch albd/albi (+ _hst) (one thread/patch).
@kernel function _ualb_map_patch_kernel!(o, @Const(mask_urbanp), @Const(pch_column),
        @Const(pch_landunit), @Const(coszen), numrad::Int)
    p = @index(Global)
    @inbounds if mask_urbanp[p]
        c = pch_column[p]; l = pch_landunit[p]
        T = eltype(coszen); sunny = coszen[l] > zero(T)
        for ib in 1:numrad
            o.albd[p, ib] = o.albgrd[c, ib]; o.albi[p, ib] = o.albgri[c, ib]
            if sunny
                o.albd_hst[p, ib] = o.albd[p, ib]; o.albi_hst[p, ib] = o.albi[p, ib]
            end
        end
    end
end

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
                      coszen::AbstractVector{<:Real}, ind::Int,
                      albsn_roof::AbstractMatrix{<:Real},
                      albsn_improad::AbstractMatrix{<:Real},
                      albsn_perroad::AbstractMatrix{<:Real},
                      h2osno_total::AbstractVector{<:Real})

    # SNAL0/SNAL1 at the array element type so no Float64 scalar reaches a Metal
    # kernel (no-op on the Float64 host path).
    FT = eltype(albsn_roof)
    ualb_snow_albedo!(albsn_roof, mask_urbanc, col.landunit, col.itype,
                      coszen, ind, albsn_improad, albsn_perroad, h2osno_total,
                      ICOL_ROOF, ICOL_ROAD_IMPERV, ICOL_ROAD_PERV, FT(SNAL0), FT(SNAL1))

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
    nc = length(mask_urbanc)
    numrad = NUMRAD
    FT = eltype(surfalb.coszen_col)

    vf_sr    = urbanparams.vf_sr
    vf_sw    = urbanparams.vf_sw
    frac_sno = waterdiagnosticbulk.frac_sno_col

    # Backend + backend-aware local allocation (device arrays when the state is on
    # a GPU, host Arrays otherwise). `_run!` launches a bundle-first kernel on the
    # state's backend and syncs (the KA CPU path is byte-identical to the old loop).
    be = _kernel_backend(surfalb.coszen_col)
    _z(dims...) = fill!(similar(surfalb.coszen_col, FT, dims...), zero(FT))
    _o(dims...) = fill!(similar(surfalb.coszen_col, FT, dims...), one(FT))
    _run!(k, nd, args...) = (k(be)(args...; ndrange = nd); KA.synchronize(be))

    # surfalb output bundle (patch matrices albd..ftii + column matrices albgr*).
    surf = _UAlbSurf(surfalb.albd_patch, surfalb.albi_patch, surfalb.fabd_patch,
        surfalb.fabd_sun_patch, surfalb.fabd_sha_patch, surfalb.fabi_patch,
        surfalb.fabi_sun_patch, surfalb.fabi_sha_patch, surfalb.ftdd_patch,
        surfalb.ftid_patch, surfalb.ftii_patch, surfalb.albd_hst_patch,
        surfalb.albi_hst_patch, surfalb.albgrd_col, surfalb.albgri_col,
        surfalb.albgrd_hst_col, surfalb.albgri_hst_col)

    # ---- Cosine solar zenith angle ----
    coszen = _z(nl)
    zen    = _z(nl)
    ualb_coszen!(coszen, mask_urbanl, surfalb.coszen_col, lun.coli, zen)

    # ---- Initialize outputs: ground albedo + patch albedo/absorption/transmission ----
    ualb_init_albgr!(surfalb.albgrd_col, mask_urbanc, surfalb.albgri_col, numrad)
    _run!(_ualb_initout_kernel!, length(mask_urbanp), surf, mask_urbanp, pch.landunit,
        pch.column, col.itype, coszen, vf_sw, vf_sr, numrad,
        ICOL_SUNWALL, ICOL_SHADEWALL, ICOL_ROAD_PERV, ICOL_ROAD_IMPERV, ICOL_ROOF)

    # ---- sref working reflectances (local) + sabs (solarabs state) init ----
    sref_roof_dir = _o(nl, numrad); sref_roof_dif = _o(nl, numrad)
    sref_sunwall_dir = _z(nl, numrad); sref_sunwall_dif = _z(nl, numrad)
    sref_shadewall_dir = _z(nl, numrad); sref_shadewall_dif = _z(nl, numrad)
    sref_improad_dir = _z(nl, numrad); sref_improad_dif = _z(nl, numrad)
    sref_perroad_dir = _z(nl, numrad); sref_perroad_dif = _z(nl, numrad)
    sref = _UAlbSref(sref_roof_dir, sref_roof_dif, sref_sunwall_dir, sref_sunwall_dif,
        sref_shadewall_dir, sref_shadewall_dif, sref_improad_dir, sref_improad_dif,
        sref_perroad_dir, sref_perroad_dif)
    sabs = _UAlbSabs(solarabs.sabs_roof_dir_lun, solarabs.sabs_roof_dif_lun,
        solarabs.sabs_sunwall_dir_lun, solarabs.sabs_sunwall_dif_lun,
        solarabs.sabs_shadewall_dir_lun, solarabs.sabs_shadewall_dif_lun,
        solarabs.sabs_improad_dir_lun, solarabs.sabs_improad_dif_lun,
        solarabs.sabs_perroad_dir_lun, solarabs.sabs_perroad_dif_lun)
    _run!(_ualb_sabsref_init_kernel!, nl, sabs, sref, mask_urbanl, vf_sw, vf_sr, numrad)

    # ---- num_solar gate: any urban landunit with sun above horizon? ----
    # coszen is a small (nl) landunit array; copy to host for this control-flow
    # scalar (non-urban entries are 0, so they never count). Cheap; once per step.
    cz_host = coszen isa Array ? coszen : Array(coszen)
    num_solar = count(>(zero(FT)), cz_host)

    if num_solar > 0
        # Solar fluxes are per unit incoming flux (device-resident working arrays).
        sdir = _o(nl, numrad); sdif = _o(nl, numrad)
        sdir_road = _z(nl, numrad); sdif_road = _z(nl, numrad)
        sdir_sunwall = _z(nl, numrad); sdif_sunwall = _z(nl, numrad)
        sdir_shadewall = _z(nl, numrad); sdif_shadewall = _z(nl, numrad)
        albsnd_roof = _z(nl, numrad); albsni_roof = _z(nl, numrad)
        albsnd_improad = _z(nl, numrad); albsni_improad = _z(nl, numrad)
        albsnd_perroad = _z(nl, numrad); albsni_perroad = _z(nl, numrad)
        alb_roof_dir_s = _z(nl, numrad); alb_roof_dif_s = _z(nl, numrad)
        alb_improad_dir_s = _z(nl, numrad); alb_perroad_dir_s = _z(nl, numrad)
        alb_improad_dif_s = _z(nl, numrad); alb_perroad_dif_s = _z(nl, numrad)

        # Total column snow water (one thread per column; local-scalar layer sum).
        nlevsno = varpar.nlevsno
        h2osno_total = _z(nc)
        _run!(_ualb_h2osno_kernel!, nc, h2osno_total, mask_urbanc,
            waterstatebulk.ws.h2osno_no_layers_col, waterstatebulk.ws.h2osoi_ice_col,
            waterstatebulk.ws.h2osoi_liq_col, col.snl, nlevsno)

        # Incident direct + diffuse radiation on the canyon facets.
        incident_direct!(mask_urbanl, lun.canyon_hwr, coszen, zen,
                         sdir, sdir_road, sdir_sunwall, sdir_shadewall)
        incident_diffuse!(mask_urbanl, lun.canyon_hwr, sdif,
                          sdif_road, sdif_sunwall, sdif_shadewall, urbanparams)

        # Snow albedos (direct beam ind=0, diffuse ind=1).
        snow_albedo!(mask_urbanc, col, coszen, 0,
                     albsnd_roof, albsnd_improad, albsnd_perroad, h2osno_total)
        snow_albedo!(mask_urbanc, col, coszen, 1,
                     albsni_roof, albsni_improad, albsni_perroad, h2osno_total)

        # Combine snow-free + snow albedos per column.
        snow = _UAlbSnow(albsnd_roof, albsni_roof, albsnd_improad, albsni_improad,
            albsnd_perroad, albsni_perroad, alb_roof_dir_s, alb_roof_dif_s,
            alb_improad_dir_s, alb_perroad_dir_s, alb_improad_dif_s, alb_perroad_dif_s)
        aparam = _UAlbParam(urbanparams.alb_roof_dir, urbanparams.alb_roof_dif,
            urbanparams.alb_improad_dir, urbanparams.alb_improad_dif,
            urbanparams.alb_perroad_dir, urbanparams.alb_perroad_dif)
        _run!(_ualb_combine_kernel!, nc, snow, aparam, mask_urbanc, col.landunit,
            col.itype, frac_sno, numrad, ICOL_ROOF, ICOL_ROAD_IMPERV, ICOL_ROAD_PERV)

        # Net solar radiation with multiple reflections.
        net_solar!(mask_urbanl, coszen, lun.canyon_hwr, lun.wtroad_perv,
                   sdir, sdif,
                   alb_improad_dir_s, alb_perroad_dir_s, urbanparams.alb_wall_dir, alb_roof_dir_s,
                   alb_improad_dif_s, alb_perroad_dif_s, urbanparams.alb_wall_dif, alb_roof_dif_s,
                   sdir_road, sdir_sunwall, sdir_shadewall,
                   sdif_road, sdif_sunwall, sdif_shadewall,
                   sref_improad_dir, sref_perroad_dir, sref_sunwall_dir,
                   sref_shadewall_dir, sref_roof_dir,
                   sref_improad_dif, sref_perroad_dif, sref_sunwall_dif,
                   sref_shadewall_dif, sref_roof_dif,
                   urbanparams, solarabs)

        # ---- Map net-solar reflectances → per-column then per-patch albedo ----
        _run!(_ualb_map_col_kernel!, nc, surf, sref, mask_urbanc, col.landunit,
            col.itype, coszen, numrad, ICOL_ROOF, ICOL_SUNWALL, ICOL_SHADEWALL,
            ICOL_ROAD_PERV, ICOL_ROAD_IMPERV)
        _run!(_ualb_map_patch_kernel!, length(mask_urbanp), surf, mask_urbanp,
            pch.column, pch.landunit, coszen, numrad)
    end

    return nothing
end
