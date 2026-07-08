# ==========================================================================
# Ported from: src/biogeophys/SurfaceRadiationMod.F90
# Calculate solar fluxes absorbed by vegetation and ground surface
#
# Public functions:
#   surfrad_init!          — allocate and initialize SurfaceRadiationData
#   surfrad_clean!         — deallocate SurfaceRadiationData
#   surfrad_init_history!  — register history fields (stub)
#   surfrad_init_cold!     — cold-start initialization (no-op)
#   is_near_local_noon     — check if current time is near local noon
#   canopy_sun_shade_fracs! — compute sun/shade fractions and PAR absorption
#   surface_radiation!     — solar fluxes absorbed by veg and ground surface
# ==========================================================================

# --- Module-level constants ---
const MPE_SURFRAD = 1.0e-06  # prevents overflow for division by zero

"""
    SurfaceRadiationData

Surface radiation diagnostic data structure. Holds aerosol forcing,
reflected/incident radiation diagnostics at the patch level.

Ported from `surfrad_type` in `SurfaceRadiationMod.F90`.
"""
Base.@kwdef mutable struct SurfaceRadiationData{FT<:Real,
                                    V<:AbstractVector{FT}}
    # --- Aerosol forcing (patch-level 1D) ---
    sfc_frc_aer_patch      ::V = Float64[]   # patch surface forcing of snow with all aerosols [W/m2]
    sfc_frc_bc_patch       ::V = Float64[]   # patch surface forcing of snow with BC [W/m2]
    sfc_frc_oc_patch       ::V = Float64[]   # patch surface forcing of snow with OC [W/m2]
    sfc_frc_dst_patch      ::V = Float64[]   # patch surface forcing of snow with dust [W/m2]
    sfc_frc_aer_sno_patch  ::V = Float64[]   # patch surface forcing of snow with all aerosols, snow-only avg [W/m2]
    sfc_frc_bc_sno_patch   ::V = Float64[]   # patch surface forcing of snow with BC, snow-only avg [W/m2]
    sfc_frc_oc_sno_patch   ::V = Float64[]   # patch surface forcing of snow with OC, snow-only avg [W/m2]
    sfc_frc_dst_sno_patch  ::V = Float64[]   # patch surface forcing of snow with dust, snow-only avg [W/m2]

    # --- Vegetation PAR at local noon (patch-level 1D) ---
    parveg_ln_patch        ::V = Float64[]   # patch absorbed par by vegetation at local noon [W/m**2]

    # --- Reflected solar from snow (patch-level 1D) ---
    fsr_sno_vd_patch       ::V = Float64[]   # patch reflected direct beam vis solar radiation from snow [W/m**2]
    fsr_sno_nd_patch       ::V = Float64[]   # patch reflected direct beam NIR solar radiation from snow [W/m**2]
    fsr_sno_vi_patch       ::V = Float64[]   # patch reflected diffuse vis solar radiation from snow [W/m**2]
    fsr_sno_ni_patch       ::V = Float64[]   # patch reflected diffuse NIR solar radiation from snow [W/m**2]

    # --- Reflected solar VIS (patch-level 1D) ---
    fsr_vis_d_patch        ::V = Float64[]   # patch reflected direct beam vis solar radiation [W/m**2]
    fsr_vis_i_patch        ::V = Float64[]   # patch reflected diffuse vis solar radiation [W/m**2]
    fsr_vis_d_ln_patch     ::V = Float64[]   # patch reflected direct beam vis solar radiation at local noon [W/m**2]

    # --- Snow-free reflected VIS diagnostics (patch-level 1D) ---
    fsrSF_vis_d_patch      ::V = Float64[]   # snow-free patch reflected direct beam vis solar radiation [W/m**2]
    fsrSF_vis_i_patch      ::V = Float64[]   # snow-free patch reflected diffuse vis solar radiation [W/m**2]
    fsrSF_vis_d_ln_patch   ::V = Float64[]   # snow-free patch reflected direct beam vis solar rad at local noon [W/m**2]

    # --- Snow radiative effect VIS (patch-level 1D) ---
    ssre_fsr_vis_d_patch   ::V = Float64[]   # snow radiative effect on direct vis reflected [W/m**2]
    ssre_fsr_vis_i_patch   ::V = Float64[]   # snow radiative effect on diffuse vis reflected [W/m**2]
    ssre_fsr_vis_d_ln_patch ::V = Float64[]  # snow radiative effect on direct vis reflected at local noon [W/m**2]

    # --- Incident solar on snow (patch-level 1D) ---
    fsds_sno_vd_patch      ::V = Float64[]   # patch incident visible, direct radiation on snow [W/m2]
    fsds_sno_nd_patch      ::V = Float64[]   # patch incident near-IR, direct radiation on snow [W/m2]
    fsds_sno_vi_patch      ::V = Float64[]   # patch incident visible, diffuse radiation on snow [W/m2]
    fsds_sno_ni_patch      ::V = Float64[]   # patch incident near-IR, diffuse radiation on snow [W/m2]

    # --- Incident solar VIS (patch-level 1D) ---
    fsds_vis_d_patch       ::V = Float64[]   # patch incident direct beam vis solar radiation [W/m**2]
    fsds_vis_i_patch       ::V = Float64[]   # patch incident diffuse vis solar radiation [W/m**2]
    fsds_vis_d_ln_patch    ::V = Float64[]   # patch incident direct beam vis solar rad at local noon [W/m**2]
    fsds_vis_i_ln_patch    ::V = Float64[]   # patch incident diffuse beam vis solar rad at local noon [W/m**2]
end

SurfaceRadiationData{FT}(; kwargs...) where {FT<:Real} =
    SurfaceRadiationData{FT, Vector{FT}}(; kwargs...)
Adapt.@adapt_structure SurfaceRadiationData


"""
    surfrad_init!(sr::SurfaceRadiationData, np::Int)

Allocate and initialize all fields of a `SurfaceRadiationData` instance for
`np` patches. Fields are initialized to `NaN`.

Ported from `surfrad_type%InitAllocate` in `SurfaceRadiationMod.F90`.
"""
function surfrad_init!(sr::SurfaceRadiationData{FT}, np::Int) where {FT}
    sr.sfc_frc_aer_patch       = fill(FT(NaN), np)
    sr.sfc_frc_bc_patch        = fill(FT(NaN), np)
    sr.sfc_frc_oc_patch        = fill(FT(NaN), np)
    sr.sfc_frc_dst_patch       = fill(FT(NaN), np)
    sr.sfc_frc_aer_sno_patch   = fill(FT(NaN), np)
    sr.sfc_frc_bc_sno_patch    = fill(FT(NaN), np)
    sr.sfc_frc_oc_sno_patch    = fill(FT(NaN), np)
    sr.sfc_frc_dst_sno_patch   = fill(FT(NaN), np)

    sr.parveg_ln_patch         = fill(FT(NaN), np)

    sr.fsr_vis_d_patch         = fill(FT(NaN), np)
    sr.fsr_vis_d_ln_patch      = fill(FT(NaN), np)
    sr.fsr_vis_i_patch         = fill(FT(NaN), np)
    sr.fsrSF_vis_d_patch       = fill(FT(NaN), np)
    sr.fsrSF_vis_d_ln_patch    = fill(FT(NaN), np)
    sr.fsrSF_vis_i_patch       = fill(FT(NaN), np)
    sr.ssre_fsr_vis_d_patch    = fill(FT(NaN), np)
    sr.ssre_fsr_vis_d_ln_patch = fill(FT(NaN), np)
    sr.ssre_fsr_vis_i_patch    = fill(FT(NaN), np)
    sr.fsr_sno_vd_patch        = fill(FT(NaN), np)
    sr.fsr_sno_nd_patch        = fill(FT(NaN), np)
    sr.fsr_sno_vi_patch        = fill(FT(NaN), np)
    sr.fsr_sno_ni_patch        = fill(FT(NaN), np)

    sr.fsds_vis_d_patch        = fill(FT(NaN), np)
    sr.fsds_vis_i_patch        = fill(FT(NaN), np)
    sr.fsds_vis_d_ln_patch     = fill(FT(NaN), np)
    sr.fsds_vis_i_ln_patch     = fill(FT(NaN), np)
    sr.fsds_sno_vd_patch       = fill(FT(NaN), np)
    sr.fsds_sno_nd_patch       = fill(FT(NaN), np)
    sr.fsds_sno_vi_patch       = fill(FT(NaN), np)
    sr.fsds_sno_ni_patch       = fill(FT(NaN), np)

    return nothing
end

"""
    surfrad_clean!(sr::SurfaceRadiationData)

Deallocate (reset to empty) all fields of a `SurfaceRadiationData` instance.
"""
function surfrad_clean!(sr::SurfaceRadiationData{FT}) where {FT}
    sr.sfc_frc_aer_patch       = FT[]
    sr.sfc_frc_bc_patch        = FT[]
    sr.sfc_frc_oc_patch        = FT[]
    sr.sfc_frc_dst_patch       = FT[]
    sr.sfc_frc_aer_sno_patch   = FT[]
    sr.sfc_frc_bc_sno_patch    = FT[]
    sr.sfc_frc_oc_sno_patch    = FT[]
    sr.sfc_frc_dst_sno_patch   = FT[]

    sr.parveg_ln_patch         = FT[]

    sr.fsr_vis_d_patch         = FT[]
    sr.fsr_vis_d_ln_patch      = FT[]
    sr.fsr_vis_i_patch         = FT[]
    sr.fsrSF_vis_d_patch       = FT[]
    sr.fsrSF_vis_d_ln_patch    = FT[]
    sr.fsrSF_vis_i_patch       = FT[]
    sr.ssre_fsr_vis_d_patch    = FT[]
    sr.ssre_fsr_vis_d_ln_patch = FT[]
    sr.ssre_fsr_vis_i_patch    = FT[]
    sr.fsr_sno_vd_patch        = FT[]
    sr.fsr_sno_nd_patch        = FT[]
    sr.fsr_sno_vi_patch        = FT[]
    sr.fsr_sno_ni_patch        = FT[]

    sr.fsds_vis_d_patch        = FT[]
    sr.fsds_vis_i_patch        = FT[]
    sr.fsds_vis_d_ln_patch     = FT[]
    sr.fsds_vis_i_ln_patch     = FT[]
    sr.fsds_sno_vd_patch       = FT[]
    sr.fsds_sno_nd_patch       = FT[]
    sr.fsds_sno_vi_patch       = FT[]
    sr.fsds_sno_ni_patch       = FT[]

    return nothing
end

"""
    surfrad_init_history!(sr::SurfaceRadiationData, bounds_patch::UnitRange{Int};
                          use_snicar_frc::Bool=false, use_SSRE::Bool=false)

Register surface radiation fields for history output. Sets fields to SPVAL.

Ported from `surfrad_type%InitHistory` in `SurfaceRadiationMod.F90`.
Requires history infrastructure — stub until histFileMod is ported.
"""
function surfrad_init_history!(sr::SurfaceRadiationData, bounds_patch::UnitRange{Int};
                                use_snicar_frc::Bool = false,
                                use_SSRE::Bool = false)
    for p in bounds_patch
        if use_snicar_frc
            sr.sfc_frc_aer_patch[p]     = SPVAL
            sr.sfc_frc_aer_sno_patch[p] = SPVAL
            sr.sfc_frc_bc_patch[p]      = SPVAL
            sr.sfc_frc_bc_sno_patch[p]  = SPVAL
            sr.sfc_frc_oc_patch[p]      = SPVAL
            sr.sfc_frc_oc_sno_patch[p]  = SPVAL
            sr.sfc_frc_dst_patch[p]     = SPVAL
            sr.sfc_frc_dst_sno_patch[p] = SPVAL
        end
        sr.fsds_vis_d_patch[p]     = SPVAL
        sr.fsds_vis_i_patch[p]     = SPVAL
        sr.fsr_vis_d_patch[p]      = SPVAL
        sr.fsr_vis_i_patch[p]      = SPVAL
        if use_SSRE
            sr.fsrSF_vis_d_patch[p]       = SPVAL
            sr.fsrSF_vis_i_patch[p]       = SPVAL
            sr.ssre_fsr_vis_d_patch[p]    = SPVAL
            sr.ssre_fsr_vis_i_patch[p]    = SPVAL
        end
        sr.fsds_vis_d_ln_patch[p]  = SPVAL
        sr.fsds_vis_i_ln_patch[p]  = SPVAL
        sr.parveg_ln_patch[p]      = SPVAL
        sr.fsr_vis_d_ln_patch[p]   = SPVAL
        if use_SSRE
            sr.fsrSF_vis_d_ln_patch[p]    = SPVAL
            sr.ssre_fsr_vis_d_ln_patch[p] = SPVAL
        end
        sr.fsds_sno_vd_patch[p]    = SPVAL
        sr.fsds_sno_nd_patch[p]    = SPVAL
        sr.fsds_sno_vi_patch[p]    = SPVAL
        sr.fsds_sno_ni_patch[p]    = SPVAL
        sr.fsr_sno_vd_patch[p]     = SPVAL
        sr.fsr_sno_nd_patch[p]     = SPVAL
        sr.fsr_sno_vi_patch[p]     = SPVAL
        sr.fsr_sno_ni_patch[p]     = SPVAL
    end
    return nothing
end

"""
    surfrad_init_cold!(sr::SurfaceRadiationData, bounds_patch::UnitRange{Int})

Initialize cold-start conditions for surface radiation (currently a no-op).

Ported from `surfrad_type%InitCold` in `SurfaceRadiationMod.F90`.
"""
function surfrad_init_cold!(sr::SurfaceRadiationData{FT}, bounds_patch::UnitRange{Int}) where {FT}
    # nothing for now (matches Fortran)
    return nothing
end

# ==========================================================================
# is_near_local_noon helper
# ==========================================================================

"""
    is_near_local_noon(londeg::Float64; deltasec::Int=1800,
                       solar_noon_secs::Float64=43200.0,
                       current_tod::Float64=43200.0)

Check if the current time-of-day is near local solar noon for a given longitude.

In the full model this is provided by `clm_time_manager`. Here we provide a
standalone version that can be called with pre-computed time-of-day.

# Arguments
- `londeg`       : longitude in degrees
- `deltasec`     : half-window around local noon (seconds)
- `solar_noon_secs` : seconds since midnight at solar noon for reference (default 43200 = 12h UTC)
- `current_tod`  : current time-of-day in seconds since midnight (UTC)

Returns `true` if within `deltasec` of local solar noon.
"""
function is_near_local_noon(londeg::Real; deltasec::Int = 1800,
                             current_tod::Real = 43200.0)
    # Local solar noon occurs when UTC time = 12h - longitude/15h
    # londeg > 0 is east, so local noon is earlier in UTC
    local_noon = 43200.0 - londeg * 240.0  # 240 = 3600/15
    # Wrap to [0, 86400)
    local_noon = mod(local_noon, 86400.0)
    diff = abs(current_tod - local_noon)
    # Handle day wrap
    if diff > 43200.0
        diff = 86400.0 - diff
    end
    return diff <= deltasec
end

# ==========================================================================
# Device helper: near-local-noon test that runs inside a kernel.
# Mirrors is_near_local_noon but eltype-generic (no Float64 literals on a
# Float32 backend) and branch-free except the day-wrap. `deltasec` is passed
# in working precision; on Float64 this is byte-identical to the host version.
# ==========================================================================
@inline function _near_local_noon(londeg::T, deltasec::T, current_tod::T) where {T}
    local_noon = T(43200) - londeg * T(240)
    local_noon = mod(local_noon, T(86400))
    d = abs(current_tod - local_noon)
    if d > T(43200)
        d = T(86400) - d
    end
    return d <= deltasec
end

# ==========================================================================
# CanopySunShadeFracs
# ==========================================================================

# --------------------------------------------------------------------------
# Device-view bundle for canopy_sun_shade_fracs!: the per-(patch, layer) and
# per-patch arrays the loop touches, aliased so the kernel body reads verbatim.
# --------------------------------------------------------------------------
Base.@kwdef struct _SSFDV{M,V}
    parsun_z_patch::M; parsha_z_patch::M
    laisun_z_patch::M; laisha_z_patch::M
    laisun_patch::V; laisha_patch::V; fsun_patch::V; elai_patch::V
    tlai_z_patch::M; fsun_z_patch::M
    fabd_sun_z_patch::M; fabd_sha_z_patch::M; fabi_sun_z_patch::M; fabi_sha_z_patch::M
end
Adapt.@adapt_structure _SSFDV

_ssf_dv(surfalb, canopystate, solarabs) = _SSFDV(;
    parsun_z_patch   = solarabs.parsun_z_patch,
    parsha_z_patch   = solarabs.parsha_z_patch,
    laisun_z_patch   = canopystate.laisun_z_patch,
    laisha_z_patch   = canopystate.laisha_z_patch,
    laisun_patch     = canopystate.laisun_patch,
    laisha_patch     = canopystate.laisha_patch,
    fsun_patch       = canopystate.fsun_patch,
    elai_patch       = canopystate.elai_patch,
    tlai_z_patch     = surfalb.tlai_z_patch,
    fsun_z_patch     = surfalb.fsun_z_patch,
    fabd_sun_z_patch = surfalb.fabd_sun_z_patch,
    fabd_sha_z_patch = surfalb.fabd_sha_z_patch,
    fabi_sun_z_patch = surfalb.fabi_sun_z_patch,
    fabi_sha_z_patch = surfalb.fabi_sha_z_patch)

# One thread per patch; the internal nrad loop is loop-carried (laisun/laisha
# sums), so it stays sequential inside the thread — fully independent per patch.
@kernel function _canopy_sun_shade_kernel!(@Const(_out), @Const(mask),
        @Const(nrad_patch), @Const(column), @Const(gridcell),
        @Const(forc_solad_col), @Const(forc_solai), dv, ipar::Int)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(dv.laisun_patch)
        nrad_p = nrad_patch[p]

        for iv in 1:nrad_p
            dv.parsun_z_patch[p, iv] = zero(T)
            dv.parsha_z_patch[p, iv] = zero(T)
            dv.laisun_z_patch[p, iv] = zero(T)
            dv.laisha_z_patch[p, iv] = zero(T)
        end

        dv.laisun_patch[p] = zero(T)
        dv.laisha_patch[p] = zero(T)
        for iv in 1:nrad_p
            dv.laisun_z_patch[p, iv] = dv.tlai_z_patch[p, iv] * dv.fsun_z_patch[p, iv]
            dv.laisha_z_patch[p, iv] = dv.tlai_z_patch[p, iv] * (one(T) - dv.fsun_z_patch[p, iv])
            dv.laisun_patch[p] += dv.laisun_z_patch[p, iv]
            dv.laisha_patch[p] += dv.laisha_z_patch[p, iv]
        end
        if dv.elai_patch[p] > zero(T)
            dv.fsun_patch[p] = dv.laisun_patch[p] / dv.elai_patch[p]
        else
            dv.fsun_patch[p] = zero(T)
        end

        g = gridcell[p]
        c = column[p]
        for iv in 1:nrad_p
            dv.parsun_z_patch[p, iv] = forc_solad_col[c, ipar] * dv.fabd_sun_z_patch[p, iv] +
                                       forc_solai[g, ipar] * dv.fabi_sun_z_patch[p, iv]
            dv.parsha_z_patch[p, iv] = forc_solad_col[c, ipar] * dv.fabd_sha_z_patch[p, iv] +
                                       forc_solai[g, ipar] * dv.fabi_sha_z_patch[p, iv]
        end
    end
end

"""
    canopy_sun_shade_fracs!(surfalb::SurfaceAlbedoData,
                            canopystate::CanopyStateData,
                            solarabs::SolarAbsorbedData,
                            forc_solad_col::Matrix{<:Real},
                            forc_solai::Matrix{<:Real},
                            pch::PatchData,
                            mask_nourbanp::AbstractVector{Bool},
                            bounds::UnitRange{Int})

Calculate sun/shade fractions and absorbed PAR for each canopy layer.

Computes:
1. Absorbed PAR for sunlit leaves in canopy layer
2. Absorbed PAR for shaded leaves in canopy layer
3. Sunlit leaf area
4. Shaded leaf area
5. Sunlit leaf area for canopy layer
6. Shaded leaf area for canopy layer
7. Sunlit fraction of canopy

Ported from `CanopySunShadeFracs` in `SurfaceRadiationMod.F90`.

# Arguments
- `surfalb`        : surface albedo data (input: nrad, tlai_z, fsun_z, fabd/fabi_sun/sha_z)
- `canopystate`    : canopy state (output: laisun, laisha, laisun_z, laisha_z, fsun; input: elai)
- `solarabs`       : solar absorbed data (output: parsun_z, parsha_z)
- `forc_solad_col` : direct beam radiation per column (ncol, numrad) [W/m**2]
- `forc_solai`     : diffuse radiation per gridcell (ngrc, numrad) [W/m**2]
- `pch`            : patch data (input: column, gridcell mappings)
- `mask_nourbanp`  : mask for non-urban patches
- `bounds`         : patch index range
"""
function canopy_sun_shade_fracs!(surfalb::SurfaceAlbedoData,
                                  canopystate::CanopyStateData,
                                  solarabs::SolarAbsorbedData,
                                  forc_solad_col::AbstractMatrix{<:Real},
                                  forc_solai::AbstractMatrix{<:Real},
                                  pch::PatchData,
                                  mask_nourbanp::AbstractVector{Bool},
                                  bounds::UnitRange{Int})
    ipar = 1  # band index for PAR

    dv = _ssf_dv(surfalb, canopystate, solarabs)
    _launch!(_canopy_sun_shade_kernel!, solarabs.parsun_z_patch, mask_nourbanp,
             surfalb.nrad_patch, pch.column, pch.gridcell,
             forc_solad_col, forc_solai, dv, ipar;
             ndrange = length(mask_nourbanp))

    return nothing
end

# --------------------------------------------------------------------------
# Zero-out canopy sunlit fraction for urban patches (masked, one thread per
# patch). Fully independent per patch — no accumulation or cross-patch reads.
# --------------------------------------------------------------------------
@kernel function _surfrad_zero_fsun_urban_kernel!(fsun_patch, @Const(mask_urbanp))
    p = @index(Global)
    @inbounds if mask_urbanp[p]
        fsun_patch[p] = zero(eltype(fsun_patch))
    end
end

"""
    surfrad_zero_fsun_urban!(fsun_patch, mask_urbanp)

Set the canopy sunlit fraction to zero for urban patches. Backend-agnostic
(CPU loop or GPU); one thread per patch.
"""
surfrad_zero_fsun_urban!(fsun_patch, mask_urbanp) =
    _launch!(_surfrad_zero_fsun_urban_kernel!, fsun_patch, mask_urbanp)

# ==========================================================================
# SurfaceRadiation — device-view bundles + per-patch kernels
#
# GPU note: every per-patch pass of SurfaceRadiation becomes a KernelAbstractions
# kernel (one thread per patch); the nband accumulation and the snow-layer pass
# keep their internal sequential loops (loop-carried sums) inside each thread, so
# there is no cross-patch dependence. The many state arrays each pass touches are
# grouped into immutable device-view bundles (Adapt.@adapt_structure'd) whose
# field names mirror the original Julia paths, so the kernel bodies read verbatim
# and CPU stays byte-identical. Scalar constants and config flags are passed as
# kernel args (ints / working-precision reals), keeping every launch under
# Metal's ~31-arg cap.
# ==========================================================================

# SolarAbsorbedData outputs the passes write/read.
Base.@kwdef struct _SRSolDV{V,M}
    sabg_soil_patch::V; sabg_snow_patch::V; sabg_patch::V; sabv_patch::V
    fsa_patch::V; fsa_r_patch::V; sabg_pen_patch::V; sub_surf_abs_SW_patch::V
    sabg_lyr_patch::M
    fsr_patch::V; fsrSF_patch::V; ssre_fsr_patch::V
    fsds_nir_d_patch::V; fsds_nir_i_patch::V; fsds_nir_d_ln_patch::V
    fsr_nir_d_patch::V; fsr_nir_i_patch::V; fsr_nir_d_ln_patch::V
    fsrSF_nir_d_patch::V; fsrSF_nir_i_patch::V; fsrSF_nir_d_ln_patch::V
    ssre_fsr_nir_d_patch::V; ssre_fsr_nir_i_patch::V
end
Adapt.@adapt_structure _SRSolDV

_sr_sol_dv(sa) = _SRSolDV(;
    sabg_soil_patch = sa.sabg_soil_patch, sabg_snow_patch = sa.sabg_snow_patch,
    sabg_patch = sa.sabg_patch, sabv_patch = sa.sabv_patch,
    fsa_patch = sa.fsa_patch, fsa_r_patch = sa.fsa_r_patch,
    sabg_pen_patch = sa.sabg_pen_patch, sub_surf_abs_SW_patch = sa.sub_surf_abs_SW_patch,
    sabg_lyr_patch = sa.sabg_lyr_patch,
    fsr_patch = sa.fsr_patch, fsrSF_patch = sa.fsrSF_patch, ssre_fsr_patch = sa.ssre_fsr_patch,
    fsds_nir_d_patch = sa.fsds_nir_d_patch, fsds_nir_i_patch = sa.fsds_nir_i_patch,
    fsds_nir_d_ln_patch = sa.fsds_nir_d_ln_patch,
    fsr_nir_d_patch = sa.fsr_nir_d_patch, fsr_nir_i_patch = sa.fsr_nir_i_patch,
    fsr_nir_d_ln_patch = sa.fsr_nir_d_ln_patch,
    fsrSF_nir_d_patch = sa.fsrSF_nir_d_patch, fsrSF_nir_i_patch = sa.fsrSF_nir_i_patch,
    fsrSF_nir_d_ln_patch = sa.fsrSF_nir_d_ln_patch,
    ssre_fsr_nir_d_patch = sa.ssre_fsr_nir_d_patch, ssre_fsr_nir_i_patch = sa.ssre_fsr_nir_i_patch)

# SurfaceRadiationData diagnostic outputs.
Base.@kwdef struct _SRRadDV{V}
    sfc_frc_bc_patch::V; sfc_frc_oc_patch::V; sfc_frc_dst_patch::V; sfc_frc_aer_patch::V
    sfc_frc_bc_sno_patch::V; sfc_frc_oc_sno_patch::V
    sfc_frc_dst_sno_patch::V; sfc_frc_aer_sno_patch::V
    fsds_vis_d_patch::V; fsds_vis_i_patch::V; fsds_vis_d_ln_patch::V; fsds_vis_i_ln_patch::V
    fsr_vis_d_patch::V; fsr_vis_i_patch::V; fsr_vis_d_ln_patch::V
    fsrSF_vis_d_patch::V; fsrSF_vis_i_patch::V; fsrSF_vis_d_ln_patch::V
    ssre_fsr_vis_d_patch::V; ssre_fsr_vis_i_patch::V
    parveg_ln_patch::V
    fsds_sno_vd_patch::V; fsds_sno_nd_patch::V; fsds_sno_vi_patch::V; fsds_sno_ni_patch::V
    fsr_sno_vd_patch::V; fsr_sno_nd_patch::V; fsr_sno_vi_patch::V; fsr_sno_ni_patch::V
end
Adapt.@adapt_structure _SRRadDV

_sr_rad_dv(sr) = _SRRadDV(;
    sfc_frc_bc_patch = sr.sfc_frc_bc_patch, sfc_frc_oc_patch = sr.sfc_frc_oc_patch,
    sfc_frc_dst_patch = sr.sfc_frc_dst_patch, sfc_frc_aer_patch = sr.sfc_frc_aer_patch,
    sfc_frc_bc_sno_patch = sr.sfc_frc_bc_sno_patch, sfc_frc_oc_sno_patch = sr.sfc_frc_oc_sno_patch,
    sfc_frc_dst_sno_patch = sr.sfc_frc_dst_sno_patch, sfc_frc_aer_sno_patch = sr.sfc_frc_aer_sno_patch,
    fsds_vis_d_patch = sr.fsds_vis_d_patch, fsds_vis_i_patch = sr.fsds_vis_i_patch,
    fsds_vis_d_ln_patch = sr.fsds_vis_d_ln_patch, fsds_vis_i_ln_patch = sr.fsds_vis_i_ln_patch,
    fsr_vis_d_patch = sr.fsr_vis_d_patch, fsr_vis_i_patch = sr.fsr_vis_i_patch,
    fsr_vis_d_ln_patch = sr.fsr_vis_d_ln_patch,
    fsrSF_vis_d_patch = sr.fsrSF_vis_d_patch, fsrSF_vis_i_patch = sr.fsrSF_vis_i_patch,
    fsrSF_vis_d_ln_patch = sr.fsrSF_vis_d_ln_patch,
    ssre_fsr_vis_d_patch = sr.ssre_fsr_vis_d_patch, ssre_fsr_vis_i_patch = sr.ssre_fsr_vis_i_patch,
    parveg_ln_patch = sr.parveg_ln_patch,
    fsds_sno_vd_patch = sr.fsds_sno_vd_patch, fsds_sno_nd_patch = sr.fsds_sno_nd_patch,
    fsds_sno_vi_patch = sr.fsds_sno_vi_patch, fsds_sno_ni_patch = sr.fsds_sno_ni_patch,
    fsr_sno_vd_patch = sr.fsr_sno_vd_patch, fsr_sno_nd_patch = sr.fsr_sno_nd_patch,
    fsr_sno_vi_patch = sr.fsr_sno_vi_patch, fsr_sno_ni_patch = sr.fsr_sno_ni_patch)

# SurfaceAlbedoData inputs (col- and patch-level albedos / canopy transfer / SNICAR).
Base.@kwdef struct _SRAlbDV{M}
    albgrd_col::M; albgri_col::M; albsod_col::M; albsoi_col::M
    albsnd_hst_col::M; albsni_hst_col::M
    albgrd_bc_col::M; albgri_bc_col::M; albgrd_oc_col::M; albgri_oc_col::M
    albgrd_dst_col::M; albgri_dst_col::M; albgrd_pur_col::M; albgri_pur_col::M
    fabd_patch::M; fabi_patch::M; ftdd_patch::M; ftid_patch::M; ftii_patch::M
    albd_patch::M; albi_patch::M; albdSF_patch::M; albiSF_patch::M
    flx_absdv_col::M; flx_absdn_col::M; flx_absiv_col::M; flx_absin_col::M
end
Adapt.@adapt_structure _SRAlbDV

_sr_alb_dv(sa) = _SRAlbDV(;
    albgrd_col = sa.albgrd_col, albgri_col = sa.albgri_col,
    albsod_col = sa.albsod_col, albsoi_col = sa.albsoi_col,
    albsnd_hst_col = sa.albsnd_hst_col, albsni_hst_col = sa.albsni_hst_col,
    albgrd_bc_col = sa.albgrd_bc_col, albgri_bc_col = sa.albgri_bc_col,
    albgrd_oc_col = sa.albgrd_oc_col, albgri_oc_col = sa.albgri_oc_col,
    albgrd_dst_col = sa.albgrd_dst_col, albgri_dst_col = sa.albgri_dst_col,
    albgrd_pur_col = sa.albgrd_pur_col, albgri_pur_col = sa.albgri_pur_col,
    fabd_patch = sa.fabd_patch, fabi_patch = sa.fabi_patch,
    ftdd_patch = sa.ftdd_patch, ftid_patch = sa.ftid_patch, ftii_patch = sa.ftii_patch,
    albd_patch = sa.albd_patch, albi_patch = sa.albi_patch,
    albdSF_patch = sa.albdSF_patch, albiSF_patch = sa.albiSF_patch,
    flx_absdv_col = sa.flx_absdv_col, flx_absdn_col = sa.flx_absdn_col,
    flx_absiv_col = sa.flx_absiv_col, flx_absin_col = sa.flx_absin_col)

# Per-patch scratch arrays carried between passes (device-resident).
Base.@kwdef struct _SRTmpDV{M,V}
    trd::M; tri::M; cad::M; cai::M
    parveg::V; sabg_pur::V; sabg_bc::V; sabg_oc::V; sabg_dst::V
end
Adapt.@adapt_structure _SRTmpDV

# --------------------------------------------------------------------------
# Pass 0: zero the per-patch flux accumulators (non-urban, masked).
# --------------------------------------------------------------------------
@kernel function _surfrad_init_kernel!(@Const(_out), @Const(mask), @Const(landunit),
        @Const(lun_itype), sol, tmp, nlevsno::Int, istsoil::Int, istcrop::Int)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(sol.sabg_patch)
        l = landunit[p]
        sol.sabg_soil_patch[p] = zero(T)
        sol.sabg_snow_patch[p] = zero(T)
        sol.sabg_patch[p]      = zero(T)
        sol.sabv_patch[p]      = zero(T)
        sol.fsa_patch[p]       = zero(T)
        if lun_itype[l] == istsoil || lun_itype[l] == istcrop
            sol.fsa_r_patch[p] = zero(T)
        end
        for j in axes(sol.sabg_lyr_patch, 2)
            sol.sabg_lyr_patch[p, j] = zero(T)
        end
        tmp.parveg[p]   = zero(T)
        tmp.sabg_pur[p] = zero(T)
        tmp.sabg_bc[p]  = zero(T)
        tmp.sabg_oc[p]  = zero(T)
        tmp.sabg_dst[p] = zero(T)
    end
end

# --------------------------------------------------------------------------
# Pass 1: nband accumulation (per patch; internal sequential ib loop). Mirrors
# the Fortran do-ib loop exactly — canopy/ground absorption sums, transmitted
# fluxes, snow/soil/total ground absorption, SNICAR no-aerosol sums.
# --------------------------------------------------------------------------
@kernel function _surfrad_nband_kernel!(@Const(_out), @Const(mask), @Const(column),
        @Const(landunit), @Const(gridcell), @Const(lun_itype), @Const(snl),
        @Const(forc_solad_col), @Const(forc_solai),
        alb, sol, tmp, nband::Int, istsoil::Int, istcrop::Int, istdlak::Int,
        use_subgrid_fluxes::Bool, use_snicar_frc::Bool)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(sol.sabg_patch)
        c = column[p]
        l = landunit[p]
        g = gridcell[p]
        is_rural = (lun_itype[l] == istsoil || lun_itype[l] == istcrop)

        for ib in 1:nband
            cad = forc_solad_col[c, ib] * alb.fabd_patch[p, ib]
            cai = forc_solai[g, ib] * alb.fabi_patch[p, ib]
            tmp.cad[p, ib] = cad
            tmp.cai[p, ib] = cai
            sol.sabv_patch[p] += cad + cai
            sol.fsa_patch[p]  += cad + cai
            if ib == 1
                tmp.parveg[p] = cad + cai
            end
            if is_rural
                sol.fsa_r_patch[p] += cad + cai
            end

            trd = forc_solad_col[c, ib] * alb.ftdd_patch[p, ib]
            tri = forc_solad_col[c, ib] * alb.ftid_patch[p, ib] +
                  forc_solai[g, ib] * alb.ftii_patch[p, ib]
            tmp.trd[p, ib] = trd
            tmp.tri[p, ib] = tri

            absrad = trd * (one(T) - alb.albsod_col[c, ib]) +
                     tri * (one(T) - alb.albsoi_col[c, ib])
            sol.sabg_soil_patch[p] += absrad

            absrad = trd * (one(T) - alb.albsnd_hst_col[c, ib]) +
                     tri * (one(T) - alb.albsni_hst_col[c, ib])
            sol.sabg_snow_patch[p] += absrad

            absrad = trd * (one(T) - alb.albgrd_col[c, ib]) +
                     tri * (one(T) - alb.albgri_col[c, ib])
            sol.sabg_patch[p] += absrad
            sol.fsa_patch[p]  += absrad
            if is_rural
                sol.fsa_r_patch[p] += absrad
            end

            if snl[c] == 0
                sol.sabg_snow_patch[p] = sol.sabg_patch[p]
                sol.sabg_soil_patch[p] = sol.sabg_patch[p]
            end
            if !use_subgrid_fluxes || lun_itype[l] == istdlak
                sol.sabg_snow_patch[p] = sol.sabg_patch[p]
                sol.sabg_soil_patch[p] = sol.sabg_patch[p]
            end

            if use_snicar_frc
                absrad_bc = trd * (one(T) - alb.albgrd_bc_col[c, ib]) +
                            tri * (one(T) - alb.albgri_bc_col[c, ib])
                tmp.sabg_bc[p] += absrad_bc
                absrad_oc = trd * (one(T) - alb.albgrd_oc_col[c, ib]) +
                            tri * (one(T) - alb.albgri_oc_col[c, ib])
                tmp.sabg_oc[p] += absrad_oc
                absrad_dst = trd * (one(T) - alb.albgrd_dst_col[c, ib]) +
                             tri * (one(T) - alb.albgri_dst_col[c, ib])
                tmp.sabg_dst[p] += absrad_dst
                absrad_pur = trd * (one(T) - alb.albgrd_pur_col[c, ib]) +
                             tri * (one(T) - alb.albgri_pur_col[c, ib])
                tmp.sabg_pur[p] += absrad_pur
            end
        end
    end
end

# --------------------------------------------------------------------------
# Pass 2: absorbed flux per snow/soil layer + SNICAR forcings (per patch;
# internal sequential layer loop). Mirrors the Fortran snow-layer block.
# --------------------------------------------------------------------------
@kernel function _surfrad_snowlyr_kernel!(@Const(_out), @Const(mask), @Const(column),
        @Const(landunit), @Const(lun_itype), @Const(snl),
        @Const(snow_depth_col), @Const(frac_sno_col),
        alb, sol, srd, tmp, nlevsno::Int, istsoil::Int, istcrop::Int, istdlak::Int,
        use_subgrid_fluxes::Bool, use_snicar_frc::Bool, do_sno_oc::Bool, spval)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(sol.sabg_patch)
        c = column[p]
        l = landunit[p]
        sabg_snl_sum = zero(T)
        sol.sub_surf_abs_SW_patch[p] = zero(T)

        if snl[c] == 0
            for j in axes(sol.sabg_lyr_patch, 2)
                sol.sabg_lyr_patch[p, j] = zero(T)
            end
            sol.sabg_lyr_patch[p, nlevsno + 1] = sol.sabg_patch[p]
            sabg_snl_sum = sol.sabg_lyr_patch[p, nlevsno + 1]
        else
            for i in (-nlevsno + 1):1
                j = i + nlevsno
                sol.sabg_lyr_patch[p, j] =
                    alb.flx_absdv_col[c, j] * tmp.trd[p, 1] +
                    alb.flx_absdn_col[c, j] * tmp.trd[p, 2] +
                    alb.flx_absiv_col[c, j] * tmp.tri[p, 1] +
                    alb.flx_absin_col[c, j] * tmp.tri[p, 2]
                if i >= snl[c] + 1
                    sabg_snl_sum += sol.sabg_lyr_patch[p, j]
                end
                if i > snl[c] + 1
                    sol.sub_surf_abs_SW_patch[p] += sol.sabg_lyr_patch[p, j]
                end
            end

            if sabg_snl_sum != zero(T)
                sol.sub_surf_abs_SW_patch[p] /= sabg_snl_sum
            else
                sol.sub_surf_abs_SW_patch[p] = zero(T)
            end

            if abs(sabg_snl_sum - sol.sabg_snow_patch[p]) > T(0.00001)
                if snl[c] == 0
                    for jj in 1:nlevsno
                        sol.sabg_lyr_patch[p, jj] = zero(T)
                    end
                    sol.sabg_lyr_patch[p, nlevsno + 1] = sol.sabg_patch[p]
                elseif snl[c] == -1
                    for jj in 1:(nlevsno - 1)
                        sol.sabg_lyr_patch[p, jj] = zero(T)
                    end
                    sol.sabg_lyr_patch[p, nlevsno]     = sol.sabg_snow_patch[p] * T(0.6)
                    sol.sabg_lyr_patch[p, nlevsno + 1] = sol.sabg_snow_patch[p] * T(0.4)
                else
                    for jj in axes(sol.sabg_lyr_patch, 2)
                        sol.sabg_lyr_patch[p, jj] = zero(T)
                    end
                    j1 = snl[c] + 1 + nlevsno
                    j2 = snl[c] + 2 + nlevsno
                    sol.sabg_lyr_patch[p, j1] = sol.sabg_snow_patch[p] * T(0.75)
                    sol.sabg_lyr_patch[p, j2] = sol.sabg_snow_patch[p] * T(0.25)
                end
            end

            if !use_subgrid_fluxes || lun_itype[l] == istdlak
                if snow_depth_col[c] < T(0.10)
                    if snl[c] == 0
                        for jj in 1:nlevsno
                            sol.sabg_lyr_patch[p, jj] = zero(T)
                        end
                        sol.sabg_lyr_patch[p, nlevsno + 1] = sol.sabg_patch[p]
                    elseif snl[c] == -1
                        for jj in 1:(nlevsno - 1)
                            sol.sabg_lyr_patch[p, jj] = zero(T)
                        end
                        sol.sabg_lyr_patch[p, nlevsno]     = sol.sabg_patch[p]
                        sol.sabg_lyr_patch[p, nlevsno + 1] = zero(T)
                    else
                        for jj in axes(sol.sabg_lyr_patch, 2)
                            sol.sabg_lyr_patch[p, jj] = zero(T)
                        end
                        j1 = snl[c] + 1 + nlevsno
                        j2 = snl[c] + 2 + nlevsno
                        sol.sabg_lyr_patch[p, j1] = sol.sabg_patch[p] * T(0.75)
                        sol.sabg_lyr_patch[p, j2] = sol.sabg_patch[p] * T(0.25)
                    end
                end
            end
        end

        if lun_itype[l] == istsoil || lun_itype[l] == istcrop
            j_top = snl[c] + 1 + nlevsno
            sol.sabg_pen_patch[p] = sol.sabg_patch[p] - sol.sabg_lyr_patch[p, j_top]
        end

        if use_snicar_frc
            srd.sfc_frc_bc_patch[p] = sol.sabg_patch[p] - tmp.sabg_bc[p]
            if do_sno_oc
                srd.sfc_frc_oc_patch[p] = sol.sabg_patch[p] - tmp.sabg_oc[p]
            else
                srd.sfc_frc_oc_patch[p] = zero(T)
            end
            srd.sfc_frc_dst_patch[p] = sol.sabg_patch[p] - tmp.sabg_dst[p]
            srd.sfc_frc_aer_patch[p] = sol.sabg_patch[p] - tmp.sabg_pur[p]

            if frac_sno_col[c] > zero(T)
                srd.sfc_frc_bc_sno_patch[p]  = srd.sfc_frc_bc_patch[p] / frac_sno_col[c]
                srd.sfc_frc_oc_sno_patch[p]  = srd.sfc_frc_oc_patch[p] / frac_sno_col[c]
                srd.sfc_frc_dst_sno_patch[p] = srd.sfc_frc_dst_patch[p] / frac_sno_col[c]
                srd.sfc_frc_aer_sno_patch[p] = srd.sfc_frc_aer_patch[p] / frac_sno_col[c]
            else
                srd.sfc_frc_bc_sno_patch[p]  = T(spval)
                srd.sfc_frc_oc_sno_patch[p]  = T(spval)
                srd.sfc_frc_dst_sno_patch[p] = T(spval)
                srd.sfc_frc_aer_sno_patch[p] = T(spval)
            end
        end
    end
end

# --------------------------------------------------------------------------
# Pass 3: radiation diagnostics (per patch, non-urban). NDVI/reflected solar,
# incident/reflected VIS+NIR, local-noon diagnostics, snow diagnostics.
# --------------------------------------------------------------------------
@kernel function _surfrad_diag_kernel!(@Const(_out), @Const(mask), @Const(column),
        @Const(gridcell), @Const(snl), @Const(londeg),
        @Const(forc_solad_col), @Const(forc_solai),
        alb, sol, srd, tmp, deltasec, current_tod,
        use_SSRE::Bool, spval)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(sol.sabg_patch)
        g = gridcell[p]
        c = column[p]

        rvis = alb.albd_patch[p, 1] * forc_solad_col[c, 1] +
               alb.albi_patch[p, 1] * forc_solai[g, 1]
        rnir = alb.albd_patch[p, 2] * forc_solad_col[c, 2] +
               alb.albi_patch[p, 2] * forc_solai[g, 2]
        sol.fsr_patch[p] = rvis + rnir

        if use_SSRE
            rvisSF = alb.albdSF_patch[p, 1] * forc_solad_col[c, 1] +
                     alb.albiSF_patch[p, 1] * forc_solai[g, 1]
            rnirSF = alb.albdSF_patch[p, 2] * forc_solad_col[c, 2] +
                     alb.albiSF_patch[p, 2] * forc_solai[g, 2]
            sol.fsrSF_patch[p] = rvisSF + rnirSF
            sol.ssre_fsr_patch[p] = sol.fsr_patch[p] - sol.fsrSF_patch[p]
        end

        srd.fsds_vis_d_patch[p] = forc_solad_col[c, 1]
        sol.fsds_nir_d_patch[p] = forc_solad_col[c, 2]
        srd.fsds_vis_i_patch[p] = forc_solai[g, 1]
        sol.fsds_nir_i_patch[p] = forc_solai[g, 2]

        srd.fsr_vis_d_patch[p]  = alb.albd_patch[p, 1] * forc_solad_col[c, 1]
        sol.fsr_nir_d_patch[p]  = alb.albd_patch[p, 2] * forc_solad_col[c, 2]
        srd.fsr_vis_i_patch[p]  = alb.albi_patch[p, 1] * forc_solai[g, 1]
        sol.fsr_nir_i_patch[p]  = alb.albi_patch[p, 2] * forc_solai[g, 2]

        if use_SSRE
            srd.fsrSF_vis_d_patch[p]  = alb.albdSF_patch[p, 1] * forc_solad_col[c, 1]
            sol.fsrSF_nir_d_patch[p]  = alb.albdSF_patch[p, 2] * forc_solad_col[c, 2]
            srd.fsrSF_vis_i_patch[p]  = alb.albiSF_patch[p, 1] * forc_solai[g, 1]
            sol.fsrSF_nir_i_patch[p]  = alb.albiSF_patch[p, 2] * forc_solai[g, 2]

            srd.ssre_fsr_vis_d_patch[p]  = srd.fsrSF_vis_d_patch[p] - srd.fsr_vis_d_patch[p]
            sol.ssre_fsr_nir_d_patch[p]  = sol.fsrSF_nir_d_patch[p] - sol.fsr_nir_d_patch[p]
            srd.ssre_fsr_vis_i_patch[p]  = srd.fsrSF_vis_i_patch[p] - srd.fsr_vis_i_patch[p]
            sol.ssre_fsr_nir_i_patch[p]  = sol.fsrSF_nir_i_patch[p] - sol.fsr_nir_i_patch[p]
        end

        near_noon = _near_local_noon(londeg[g], deltasec, current_tod)
        if near_noon
            srd.fsds_vis_d_ln_patch[p] = forc_solad_col[c, 1]
            sol.fsds_nir_d_ln_patch[p] = forc_solad_col[c, 2]
            srd.fsr_vis_d_ln_patch[p]  = alb.albd_patch[p, 1] * forc_solad_col[c, 1]
            sol.fsr_nir_d_ln_patch[p]  = alb.albd_patch[p, 2] * forc_solad_col[c, 2]
            srd.fsds_vis_i_ln_patch[p] = forc_solai[g, 1]
            srd.parveg_ln_patch[p]     = tmp.parveg[p]
        else
            srd.fsds_vis_d_ln_patch[p] = T(spval)
            sol.fsds_nir_d_ln_patch[p] = T(spval)
            srd.fsr_vis_d_ln_patch[p]  = T(spval)
            sol.fsr_nir_d_ln_patch[p]  = T(spval)
            srd.fsds_vis_i_ln_patch[p] = T(spval)
            srd.parveg_ln_patch[p]     = T(spval)
        end

        if use_SSRE
            if near_noon
                srd.fsrSF_vis_d_ln_patch[p] = alb.albdSF_patch[p, 1] * forc_solad_col[c, 1]
                sol.fsrSF_nir_d_ln_patch[p] = alb.albdSF_patch[p, 2] * forc_solad_col[c, 2]
            else
                srd.fsrSF_vis_d_ln_patch[p] = T(spval)
                sol.fsrSF_nir_d_ln_patch[p] = T(spval)
            end
        end

        if snl[c] < 0
            srd.fsds_sno_vd_patch[p] = forc_solad_col[c, 1]
            srd.fsds_sno_nd_patch[p] = forc_solad_col[c, 2]
            srd.fsds_sno_vi_patch[p] = forc_solai[g, 1]
            srd.fsds_sno_ni_patch[p] = forc_solai[g, 2]
            srd.fsr_sno_vd_patch[p] = srd.fsds_vis_d_patch[p] * alb.albsnd_hst_col[c, 1]
            srd.fsr_sno_nd_patch[p] = sol.fsds_nir_d_patch[p] * alb.albsnd_hst_col[c, 2]
            srd.fsr_sno_vi_patch[p] = srd.fsds_vis_i_patch[p] * alb.albsni_hst_col[c, 1]
            srd.fsr_sno_ni_patch[p] = sol.fsds_nir_i_patch[p] * alb.albsni_hst_col[c, 2]
        else
            srd.fsds_sno_vd_patch[p] = T(spval)
            srd.fsds_sno_nd_patch[p] = T(spval)
            srd.fsds_sno_vi_patch[p] = T(spval)
            srd.fsds_sno_ni_patch[p] = T(spval)
            srd.fsr_sno_vd_patch[p] = T(spval)
            srd.fsr_sno_nd_patch[p] = T(spval)
            srd.fsr_sno_vi_patch[p] = T(spval)
            srd.fsr_sno_ni_patch[p] = T(spval)
        end
    end
end

# --------------------------------------------------------------------------
# Pass 4: urban patches (per patch, urban mask). Incident/reflected solar +
# local-noon diagnostics; fsr is summed from VIS+NIR reflected.
# --------------------------------------------------------------------------
@kernel function _surfrad_urban_kernel!(@Const(_out), @Const(mask), @Const(column),
        @Const(gridcell), @Const(londeg),
        @Const(forc_solad_col), @Const(forc_solai),
        alb, sol, srd, deltasec, current_tod, spval)
    p = @index(Global)
    @inbounds if mask[p]
        T = eltype(sol.fsr_patch)
        g = gridcell[p]
        c = column[p]

        srd.fsds_vis_d_patch[p] = forc_solad_col[c, 1]
        sol.fsds_nir_d_patch[p] = forc_solad_col[c, 2]
        srd.fsds_vis_i_patch[p] = forc_solai[g, 1]
        sol.fsds_nir_i_patch[p] = forc_solai[g, 2]

        near_noon = _near_local_noon(londeg[g], deltasec, current_tod)
        if near_noon
            srd.fsds_vis_d_ln_patch[p] = forc_solad_col[c, 1]
            sol.fsds_nir_d_ln_patch[p] = forc_solad_col[c, 2]
            srd.fsds_vis_i_ln_patch[p] = forc_solai[g, 1]
            srd.parveg_ln_patch[p]     = zero(T)
        else
            srd.fsds_vis_d_ln_patch[p] = T(spval)
            sol.fsds_nir_d_ln_patch[p] = T(spval)
            srd.fsds_vis_i_ln_patch[p] = T(spval)
            srd.parveg_ln_patch[p]     = T(spval)
        end

        srd.fsr_vis_d_patch[p]  = alb.albd_patch[p, 1] * forc_solad_col[c, 1]
        sol.fsr_nir_d_patch[p]  = alb.albd_patch[p, 2] * forc_solad_col[c, 2]
        srd.fsr_vis_i_patch[p]  = alb.albi_patch[p, 1] * forc_solai[g, 1]
        sol.fsr_nir_i_patch[p]  = alb.albi_patch[p, 2] * forc_solai[g, 2]

        if near_noon
            srd.fsr_vis_d_ln_patch[p] = srd.fsr_vis_d_patch[p]
            sol.fsr_nir_d_ln_patch[p] = sol.fsr_nir_d_patch[p]
        else
            srd.fsr_vis_d_ln_patch[p] = T(spval)
            sol.fsr_nir_d_ln_patch[p] = T(spval)
        end

        sol.fsr_patch[p] = srd.fsr_vis_d_patch[p] + sol.fsr_nir_d_patch[p] +
                           srd.fsr_vis_i_patch[p] + sol.fsr_nir_i_patch[p]
    end
end

# ==========================================================================
# SurfaceRadiation
# ==========================================================================

"""
    surface_radiation!(surfalb::SurfaceAlbedoData,
                       canopystate::CanopyStateData,
                       solarabs::SolarAbsorbedData,
                       surfrad::SurfaceRadiationData,
                       waterdiag::WaterDiagnosticBulkData,
                       col::ColumnData,
                       lun::LandunitData,
                       grc::GridcellData,
                       pch::PatchData,
                       forc_solad_col::Matrix{<:Real},
                       forc_solai::Matrix{<:Real},
                       mask_nourbanp::AbstractVector{Bool},
                       mask_urbanp::AbstractVector{Bool},
                       bounds::UnitRange{Int};
                       dtime::Float64 = 3600.0,
                       current_tod::Float64 = 43200.0,
                       use_subgrid_fluxes::Bool = true,
                       use_snicar_frc::Bool = false,
                       use_SSRE::Bool = false,
                       do_sno_oc::Bool = false)

Calculate solar fluxes absorbed by vegetation and ground surface.

Ported from `SurfaceRadiation` in `SurfaceRadiationMod.F90`.

# Arguments
- `surfalb`        : surface albedo data (input)
- `canopystate`    : canopy state (input: tlai, elai, esai; output: fsun)
- `solarabs`       : solar absorbed data (output: fsa, fsr, sabv, sabg, sabg_lyr, etc.)
- `surfrad`        : surface radiation diagnostics (output)
- `waterdiag`      : water diagnostic bulk data (input: snow_depth, frac_sno, frac_sno_eff)
- `col`            : column data (input: snl)
- `lun`            : landunit data (input: itype)
- `grc`            : gridcell data (input: londeg)
- `pch`            : patch data (input: column, gridcell, landunit mappings)
- `forc_solad_col` : direct beam radiation per column (ncol, numrad) [W/m**2]
- `forc_solai`     : diffuse radiation per gridcell (ngrc, numrad) [W/m**2]
- `mask_nourbanp`  : mask for non-urban patches
- `mask_urbanp`    : mask for urban patches
- `bounds`         : patch index range
- `dtime`          : time step size in seconds
- `current_tod`    : current time-of-day in seconds since midnight (UTC)
- `use_subgrid_fluxes` : whether to use subgrid fluxes
- `use_snicar_frc` : whether to use SNICAR aerosol forcing
- `use_SSRE`       : whether to compute snow-surface radiative effect
- `do_sno_oc`      : whether to compute OC snow forcing
"""
function surface_radiation!(surfalb::SurfaceAlbedoData,
                              canopystate::CanopyStateData,
                              solarabs::SolarAbsorbedData,
                              surfrad::SurfaceRadiationData,
                              waterdiag::WaterDiagnosticBulkData,
                              col::ColumnData,
                              lun::LandunitData,
                              grc::GridcellData,
                              pch::PatchData,
                              forc_solad_col::AbstractMatrix{<:Real},
                              forc_solai::AbstractMatrix{<:Real},
                              mask_nourbanp::AbstractVector{Bool},
                              mask_urbanp::AbstractVector{Bool},
                              bounds::UnitRange{Int};
                              dtime::Real = 3600.0,
                              current_tod::Real = 43200.0,
                              use_subgrid_fluxes::Bool = true,
                              use_snicar_frc::Bool = false,
                              use_SSRE::Bool = false,
                              do_sno_oc::Bool = false)
    nband = NUMRAD
    nlevsno = varpar.nlevsno

    # Working precision / backend prototype taken from a state array (so device
    # arrays drive scratch allocation); literal config flags resolved on host.
    FT = eltype(solarabs.sabg_patch)
    endp = last(bounds)
    proto = solarabs.sabg_patch

    # Per-patch scratch carried across passes — device-resident (similar()+fill!,
    # NOT zeros()), so a Float32/Metal state runs the whole function on-device.
    mkV() = fill!(similar(proto, endp), zero(eltype(proto)))
    mkM() = fill!(similar(proto, endp, nband), zero(eltype(proto)))
    tmp = _SRTmpDV(; trd = mkM(), tri = mkM(), cad = mkM(), cai = mkM(),
                    parveg = mkV(), sabg_pur = mkV(), sabg_bc = mkV(),
                    sabg_oc = mkV(), sabg_dst = mkV())

    # Device-view bundles (shared array refs → writes flow back into the structs).
    sol = _sr_sol_dv(solarabs)
    srd = _sr_rad_dv(surfrad)
    alb = _sr_alb_dv(surfalb)

    # Scalar constants/flags at working precision (no Float64 on a Float32 backend).
    spval     = convert(FT, SPVAL)
    current_t = convert(FT, current_tod)
    deltasec  = convert(FT, div(round(Int, dtime), 2))

    nd = length(mask_nourbanp)

    # --- Pass 0: initialize fluxes (non-urban) ---
    _launch!(_surfrad_init_kernel!, solarabs.sabg_patch, mask_nourbanp,
             pch.landunit, lun.itype, sol, tmp, nlevsno, ISTSOIL, ISTCROP;
             ndrange = nd)

    # Zero-out fsun for urban patches
    surfrad_zero_fsun_urban!(canopystate.fsun_patch, mask_urbanp)

    # --- Pass 1: nband accumulation (non-urban) ---
    _launch!(_surfrad_nband_kernel!, solarabs.sabg_patch, mask_nourbanp,
             pch.column, pch.landunit, pch.gridcell, lun.itype, col.snl,
             forc_solad_col, forc_solai, alb, sol, tmp, nband,
             ISTSOIL, ISTCROP, ISTDLAK, use_subgrid_fluxes, use_snicar_frc;
             ndrange = nd)

    # --- Pass 2: absorbed flux per snow/soil layer + SNICAR forcings ---
    _launch!(_surfrad_snowlyr_kernel!, solarabs.sabg_patch, mask_nourbanp,
             pch.column, pch.landunit, lun.itype, col.snl,
             waterdiag.snow_depth_col, waterdiag.frac_sno_col,
             alb, sol, srd, tmp, nlevsno, ISTSOIL, ISTCROP, ISTDLAK,
             use_subgrid_fluxes, use_snicar_frc, do_sno_oc, spval;
             ndrange = nd)

    # --- Pass 3: radiation diagnostics (non-urban) ---
    _launch!(_surfrad_diag_kernel!, solarabs.sabg_patch, mask_nourbanp,
             pch.column, pch.gridcell, col.snl, grc.londeg,
             forc_solad_col, forc_solai, alb, sol, srd, tmp,
             deltasec, current_t, use_SSRE, spval; ndrange = nd)

    # --- Pass 4: urban patches ---
    _launch!(_surfrad_urban_kernel!, solarabs.fsr_patch, mask_urbanp,
             pch.column, pch.gridcell, grc.londeg,
             forc_solad_col, forc_solai, alb, sol, srd,
             deltasec, current_t, spval; ndrange = length(mask_urbanp))

    return nothing
end
