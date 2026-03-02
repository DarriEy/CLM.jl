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
Base.@kwdef mutable struct SurfaceRadiationData
    # --- Aerosol forcing (patch-level 1D) ---
    sfc_frc_aer_patch      ::Vector{Float64} = Float64[]   # patch surface forcing of snow with all aerosols [W/m2]
    sfc_frc_bc_patch       ::Vector{Float64} = Float64[]   # patch surface forcing of snow with BC [W/m2]
    sfc_frc_oc_patch       ::Vector{Float64} = Float64[]   # patch surface forcing of snow with OC [W/m2]
    sfc_frc_dst_patch      ::Vector{Float64} = Float64[]   # patch surface forcing of snow with dust [W/m2]
    sfc_frc_aer_sno_patch  ::Vector{Float64} = Float64[]   # patch surface forcing of snow with all aerosols, snow-only avg [W/m2]
    sfc_frc_bc_sno_patch   ::Vector{Float64} = Float64[]   # patch surface forcing of snow with BC, snow-only avg [W/m2]
    sfc_frc_oc_sno_patch   ::Vector{Float64} = Float64[]   # patch surface forcing of snow with OC, snow-only avg [W/m2]
    sfc_frc_dst_sno_patch  ::Vector{Float64} = Float64[]   # patch surface forcing of snow with dust, snow-only avg [W/m2]

    # --- Vegetation PAR at local noon (patch-level 1D) ---
    parveg_ln_patch        ::Vector{Float64} = Float64[]   # patch absorbed par by vegetation at local noon [W/m**2]

    # --- Reflected solar from snow (patch-level 1D) ---
    fsr_sno_vd_patch       ::Vector{Float64} = Float64[]   # patch reflected direct beam vis solar radiation from snow [W/m**2]
    fsr_sno_nd_patch       ::Vector{Float64} = Float64[]   # patch reflected direct beam NIR solar radiation from snow [W/m**2]
    fsr_sno_vi_patch       ::Vector{Float64} = Float64[]   # patch reflected diffuse vis solar radiation from snow [W/m**2]
    fsr_sno_ni_patch       ::Vector{Float64} = Float64[]   # patch reflected diffuse NIR solar radiation from snow [W/m**2]

    # --- Reflected solar VIS (patch-level 1D) ---
    fsr_vis_d_patch        ::Vector{Float64} = Float64[]   # patch reflected direct beam vis solar radiation [W/m**2]
    fsr_vis_i_patch        ::Vector{Float64} = Float64[]   # patch reflected diffuse vis solar radiation [W/m**2]
    fsr_vis_d_ln_patch     ::Vector{Float64} = Float64[]   # patch reflected direct beam vis solar radiation at local noon [W/m**2]

    # --- Snow-free reflected VIS diagnostics (patch-level 1D) ---
    fsrSF_vis_d_patch      ::Vector{Float64} = Float64[]   # snow-free patch reflected direct beam vis solar radiation [W/m**2]
    fsrSF_vis_i_patch      ::Vector{Float64} = Float64[]   # snow-free patch reflected diffuse vis solar radiation [W/m**2]
    fsrSF_vis_d_ln_patch   ::Vector{Float64} = Float64[]   # snow-free patch reflected direct beam vis solar rad at local noon [W/m**2]

    # --- Snow radiative effect VIS (patch-level 1D) ---
    ssre_fsr_vis_d_patch   ::Vector{Float64} = Float64[]   # snow radiative effect on direct vis reflected [W/m**2]
    ssre_fsr_vis_i_patch   ::Vector{Float64} = Float64[]   # snow radiative effect on diffuse vis reflected [W/m**2]
    ssre_fsr_vis_d_ln_patch ::Vector{Float64} = Float64[]  # snow radiative effect on direct vis reflected at local noon [W/m**2]

    # --- Incident solar on snow (patch-level 1D) ---
    fsds_sno_vd_patch      ::Vector{Float64} = Float64[]   # patch incident visible, direct radiation on snow [W/m2]
    fsds_sno_nd_patch      ::Vector{Float64} = Float64[]   # patch incident near-IR, direct radiation on snow [W/m2]
    fsds_sno_vi_patch      ::Vector{Float64} = Float64[]   # patch incident visible, diffuse radiation on snow [W/m2]
    fsds_sno_ni_patch      ::Vector{Float64} = Float64[]   # patch incident near-IR, diffuse radiation on snow [W/m2]

    # --- Incident solar VIS (patch-level 1D) ---
    fsds_vis_d_patch       ::Vector{Float64} = Float64[]   # patch incident direct beam vis solar radiation [W/m**2]
    fsds_vis_i_patch       ::Vector{Float64} = Float64[]   # patch incident diffuse vis solar radiation [W/m**2]
    fsds_vis_d_ln_patch    ::Vector{Float64} = Float64[]   # patch incident direct beam vis solar rad at local noon [W/m**2]
    fsds_vis_i_ln_patch    ::Vector{Float64} = Float64[]   # patch incident diffuse beam vis solar rad at local noon [W/m**2]
end

"""
    surfrad_init!(sr::SurfaceRadiationData, np::Int)

Allocate and initialize all fields of a `SurfaceRadiationData` instance for
`np` patches. Fields are initialized to `NaN`.

Ported from `surfrad_type%InitAllocate` in `SurfaceRadiationMod.F90`.
"""
function surfrad_init!(sr::SurfaceRadiationData, np::Int)
    sr.sfc_frc_aer_patch       = fill(NaN, np)
    sr.sfc_frc_bc_patch        = fill(NaN, np)
    sr.sfc_frc_oc_patch        = fill(NaN, np)
    sr.sfc_frc_dst_patch       = fill(NaN, np)
    sr.sfc_frc_aer_sno_patch   = fill(NaN, np)
    sr.sfc_frc_bc_sno_patch    = fill(NaN, np)
    sr.sfc_frc_oc_sno_patch    = fill(NaN, np)
    sr.sfc_frc_dst_sno_patch   = fill(NaN, np)

    sr.parveg_ln_patch         = fill(NaN, np)

    sr.fsr_vis_d_patch         = fill(NaN, np)
    sr.fsr_vis_d_ln_patch      = fill(NaN, np)
    sr.fsr_vis_i_patch         = fill(NaN, np)
    sr.fsrSF_vis_d_patch       = fill(NaN, np)
    sr.fsrSF_vis_d_ln_patch    = fill(NaN, np)
    sr.fsrSF_vis_i_patch       = fill(NaN, np)
    sr.ssre_fsr_vis_d_patch    = fill(NaN, np)
    sr.ssre_fsr_vis_d_ln_patch = fill(NaN, np)
    sr.ssre_fsr_vis_i_patch    = fill(NaN, np)
    sr.fsr_sno_vd_patch        = fill(NaN, np)
    sr.fsr_sno_nd_patch        = fill(NaN, np)
    sr.fsr_sno_vi_patch        = fill(NaN, np)
    sr.fsr_sno_ni_patch        = fill(NaN, np)

    sr.fsds_vis_d_patch        = fill(NaN, np)
    sr.fsds_vis_i_patch        = fill(NaN, np)
    sr.fsds_vis_d_ln_patch     = fill(NaN, np)
    sr.fsds_vis_i_ln_patch     = fill(NaN, np)
    sr.fsds_sno_vd_patch       = fill(NaN, np)
    sr.fsds_sno_nd_patch       = fill(NaN, np)
    sr.fsds_sno_vi_patch       = fill(NaN, np)
    sr.fsds_sno_ni_patch       = fill(NaN, np)

    return nothing
end

"""
    surfrad_clean!(sr::SurfaceRadiationData)

Deallocate (reset to empty) all fields of a `SurfaceRadiationData` instance.
"""
function surfrad_clean!(sr::SurfaceRadiationData)
    sr.sfc_frc_aer_patch       = Float64[]
    sr.sfc_frc_bc_patch        = Float64[]
    sr.sfc_frc_oc_patch        = Float64[]
    sr.sfc_frc_dst_patch       = Float64[]
    sr.sfc_frc_aer_sno_patch   = Float64[]
    sr.sfc_frc_bc_sno_patch    = Float64[]
    sr.sfc_frc_oc_sno_patch    = Float64[]
    sr.sfc_frc_dst_sno_patch   = Float64[]

    sr.parveg_ln_patch         = Float64[]

    sr.fsr_vis_d_patch         = Float64[]
    sr.fsr_vis_d_ln_patch      = Float64[]
    sr.fsr_vis_i_patch         = Float64[]
    sr.fsrSF_vis_d_patch       = Float64[]
    sr.fsrSF_vis_d_ln_patch    = Float64[]
    sr.fsrSF_vis_i_patch       = Float64[]
    sr.ssre_fsr_vis_d_patch    = Float64[]
    sr.ssre_fsr_vis_d_ln_patch = Float64[]
    sr.ssre_fsr_vis_i_patch    = Float64[]
    sr.fsr_sno_vd_patch        = Float64[]
    sr.fsr_sno_nd_patch        = Float64[]
    sr.fsr_sno_vi_patch        = Float64[]
    sr.fsr_sno_ni_patch        = Float64[]

    sr.fsds_vis_d_patch        = Float64[]
    sr.fsds_vis_i_patch        = Float64[]
    sr.fsds_vis_d_ln_patch     = Float64[]
    sr.fsds_vis_i_ln_patch     = Float64[]
    sr.fsds_sno_vd_patch       = Float64[]
    sr.fsds_sno_nd_patch       = Float64[]
    sr.fsds_sno_vi_patch       = Float64[]
    sr.fsds_sno_ni_patch       = Float64[]

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
function surfrad_init_cold!(sr::SurfaceRadiationData, bounds_patch::UnitRange{Int})
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
function is_near_local_noon(londeg::Float64; deltasec::Int = 1800,
                             current_tod::Float64 = 43200.0)
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
# CanopySunShadeFracs
# ==========================================================================

"""
    canopy_sun_shade_fracs!(surfalb::SurfaceAlbedoData,
                            canopystate::CanopyStateData,
                            solarabs::SolarAbsorbedData,
                            forc_solad_col::Matrix{Float64},
                            forc_solai::Matrix{Float64},
                            pch::PatchData,
                            mask_nourbanp::BitVector,
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
                                  forc_solad_col::Matrix{Float64},
                                  forc_solai::Matrix{Float64},
                                  pch::PatchData,
                                  mask_nourbanp::BitVector,
                                  bounds::UnitRange{Int})
    ipar = 1  # band index for PAR

    for p in bounds
        mask_nourbanp[p] || continue

        nrad_p = surfalb.nrad_patch[p]

        # Initialize layer values
        for iv in 1:nrad_p
            solarabs.parsun_z_patch[p, iv] = 0.0
            solarabs.parsha_z_patch[p, iv] = 0.0
            canopystate.laisun_z_patch[p, iv] = 0.0
            canopystate.laisha_z_patch[p, iv] = 0.0
        end

        # Compute laisun_z and laisha_z for each layer.
        # Derive canopy laisun, laisha, and fsun from layer sums.
        canopystate.laisun_patch[p] = 0.0
        canopystate.laisha_patch[p] = 0.0
        for iv in 1:nrad_p
            canopystate.laisun_z_patch[p, iv] = surfalb.tlai_z_patch[p, iv] * surfalb.fsun_z_patch[p, iv]
            canopystate.laisha_z_patch[p, iv] = surfalb.tlai_z_patch[p, iv] * (1.0 - surfalb.fsun_z_patch[p, iv])
            canopystate.laisun_patch[p] += canopystate.laisun_z_patch[p, iv]
            canopystate.laisha_patch[p] += canopystate.laisha_z_patch[p, iv]
        end
        if canopystate.elai_patch[p] > 0.0
            canopystate.fsun_patch[p] = canopystate.laisun_patch[p] / canopystate.elai_patch[p]
        else
            canopystate.fsun_patch[p] = 0.0
        end

        # Absorbed PAR profile through canopy
        g = pch.gridcell[p]
        c = pch.column[p]

        for iv in 1:nrad_p
            solarabs.parsun_z_patch[p, iv] = forc_solad_col[c, ipar] * surfalb.fabd_sun_z_patch[p, iv] +
                                              forc_solai[g, ipar] * surfalb.fabi_sun_z_patch[p, iv]
            solarabs.parsha_z_patch[p, iv] = forc_solad_col[c, ipar] * surfalb.fabd_sha_z_patch[p, iv] +
                                              forc_solai[g, ipar] * surfalb.fabi_sha_z_patch[p, iv]
        end
    end

    return nothing
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
                       forc_solad_col::Matrix{Float64},
                       forc_solai::Matrix{Float64},
                       mask_nourbanp::BitVector,
                       mask_urbanp::BitVector,
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
                              forc_solad_col::Matrix{Float64},
                              forc_solai::Matrix{Float64},
                              mask_nourbanp::BitVector,
                              mask_urbanp::BitVector,
                              bounds::UnitRange{Int};
                              dtime::Float64 = 3600.0,
                              current_tod::Float64 = 43200.0,
                              use_subgrid_fluxes::Bool = true,
                              use_snicar_frc::Bool = false,
                              use_SSRE::Bool = false,
                              do_sno_oc::Bool = false)
    nband = NUMRAD
    nlevsno = varpar.nlevsno

    # Number of patches
    np = length(bounds)

    # Temporary arrays (patch-level)
    trd = zeros(last(bounds), nband)
    tri = zeros(last(bounds), nband)
    cad = zeros(last(bounds), nband)
    cai = zeros(last(bounds), nband)
    parveg = zeros(last(bounds))
    sabg_pur = zeros(last(bounds))
    sabg_bc  = zeros(last(bounds))
    sabg_oc  = zeros(last(bounds))
    sabg_dst = zeros(last(bounds))

    # --- Initialize fluxes ---
    for p in bounds
        mask_nourbanp[p] || continue
        l = pch.landunit[p]

        solarabs.sabg_soil_patch[p] = 0.0
        solarabs.sabg_snow_patch[p] = 0.0
        solarabs.sabg_patch[p]      = 0.0
        solarabs.sabv_patch[p]      = 0.0
        solarabs.fsa_patch[p]       = 0.0
        if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
            solarabs.fsa_r_patch[p] = 0.0
        end
        solarabs.sabg_lyr_patch[p, :] .= 0.0
        sabg_pur[p] = 0.0
        sabg_bc[p]  = 0.0
        sabg_oc[p]  = 0.0
        sabg_dst[p] = 0.0
    end

    # Zero-out fsun for urban patches
    for p in bounds
        mask_urbanp[p] || continue
        canopystate.fsun_patch[p] = 0.0
    end

    # --- Loop over nband wavebands ---
    for ib in 1:nband
        for p in bounds
            mask_nourbanp[p] || continue
            c = pch.column[p]
            l = pch.landunit[p]
            g = pch.gridcell[p]

            # Absorbed by canopy
            cad[p, ib] = forc_solad_col[c, ib] * surfalb.fabd_patch[p, ib]
            cai[p, ib] = forc_solai[g, ib] * surfalb.fabi_patch[p, ib]
            solarabs.sabv_patch[p] += cad[p, ib] + cai[p, ib]
            solarabs.fsa_patch[p]  += cad[p, ib] + cai[p, ib]
            if ib == 1
                parveg[p] = cad[p, ib] + cai[p, ib]
            end
            if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                solarabs.fsa_r_patch[p] += cad[p, ib] + cai[p, ib]
            end

            # Transmitted = solar fluxes incident on ground
            trd[p, ib] = forc_solad_col[c, ib] * surfalb.ftdd_patch[p, ib]
            tri[p, ib] = forc_solad_col[c, ib] * surfalb.ftid_patch[p, ib] +
                          forc_solai[g, ib] * surfalb.ftii_patch[p, ib]

            # Solar radiation absorbed by ground surface
            # Calculate absorbed solar by soil/snow separately
            absrad = trd[p, ib] * (1.0 - surfalb.albsod_col[c, ib]) +
                     tri[p, ib] * (1.0 - surfalb.albsoi_col[c, ib])
            solarabs.sabg_soil_patch[p] += absrad

            absrad = trd[p, ib] * (1.0 - surfalb.albsnd_hst_col[c, ib]) +
                     tri[p, ib] * (1.0 - surfalb.albsni_hst_col[c, ib])
            solarabs.sabg_snow_patch[p] += absrad

            absrad = trd[p, ib] * (1.0 - surfalb.albgrd_col[c, ib]) +
                     tri[p, ib] * (1.0 - surfalb.albgri_col[c, ib])
            solarabs.sabg_patch[p] += absrad
            solarabs.fsa_patch[p]  += absrad
            if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
                solarabs.fsa_r_patch[p] += absrad
            end

            if col.snl[c] == 0
                solarabs.sabg_snow_patch[p] = solarabs.sabg_patch[p]
                solarabs.sabg_soil_patch[p] = solarabs.sabg_patch[p]
            end

            # if no subgrid fluxes, make sure to set both components equal to weighted average
            if !use_subgrid_fluxes || lun.itype[l] == ISTDLAK
                solarabs.sabg_snow_patch[p] = solarabs.sabg_patch[p]
                solarabs.sabg_soil_patch[p] = solarabs.sabg_patch[p]
            end

            if use_snicar_frc
                # Solar radiation absorbed by ground surface without BC
                absrad_bc = trd[p, ib] * (1.0 - surfalb.albgrd_bc_col[c, ib]) +
                            tri[p, ib] * (1.0 - surfalb.albgri_bc_col[c, ib])
                sabg_bc[p] += absrad_bc

                # Solar radiation absorbed by ground surface without OC
                absrad_oc = trd[p, ib] * (1.0 - surfalb.albgrd_oc_col[c, ib]) +
                            tri[p, ib] * (1.0 - surfalb.albgri_oc_col[c, ib])
                sabg_oc[p] += absrad_oc

                # Solar radiation absorbed by ground surface without dust
                absrad_dst = trd[p, ib] * (1.0 - surfalb.albgrd_dst_col[c, ib]) +
                             tri[p, ib] * (1.0 - surfalb.albgri_dst_col[c, ib])
                sabg_dst[p] += absrad_dst

                # Solar radiation absorbed by ground surface without any aerosols
                absrad_pur = trd[p, ib] * (1.0 - surfalb.albgrd_pur_col[c, ib]) +
                             tri[p, ib] * (1.0 - surfalb.albgri_pur_col[c, ib])
                sabg_pur[p] += absrad_pur
            end
        end # end patch loop
    end # end nband loop

    # --- Compute absorbed flux in each snow layer and top soil layer ---
    # Based on flux factors computed in the radiative transfer portion of SNICAR.
    # Fortran indices: sabg_lyr(p, -nlevsno+1 : 1)
    # Julia indices:   sabg_lyr_patch[p, 1 : nlevsno+1]
    # Mapping: fortran_i → julia_j = fortran_i + nlevsno
    # So fortran -nlevsno+1 → 1, fortran 0 → nlevsno, fortran 1 → nlevsno+1

    for p in bounds
        mask_nourbanp[p] || continue
        c = pch.column[p]
        l = pch.landunit[p]
        sabg_snl_sum = 0.0

        solarabs.sub_surf_abs_SW_patch[p] = 0.0

        # CASE1: No snow layers: all energy is absorbed in top soil layer
        if col.snl[c] == 0
            solarabs.sabg_lyr_patch[p, :] .= 0.0
            solarabs.sabg_lyr_patch[p, nlevsno + 1] = solarabs.sabg_patch[p]
            sabg_snl_sum = solarabs.sabg_lyr_patch[p, nlevsno + 1]

        # CASE 2: Snow layers present
        else
            for i in (-nlevsno + 1):1
                j = i + nlevsno  # Julia 1-based index
                solarabs.sabg_lyr_patch[p, j] =
                    surfalb.flx_absdv_col[c, j] * trd[p, 1] +
                    surfalb.flx_absdn_col[c, j] * trd[p, 2] +
                    surfalb.flx_absiv_col[c, j] * tri[p, 1] +
                    surfalb.flx_absin_col[c, j] * tri[p, 2]

                # summed radiation in active snow layers
                if i >= col.snl[c] + 1
                    sabg_snl_sum += solarabs.sabg_lyr_patch[p, j]
                end
                # accumulate subsurface flux as diagnostic
                if i > col.snl[c] + 1
                    solarabs.sub_surf_abs_SW_patch[p] += solarabs.sabg_lyr_patch[p, j]
                end
            end

            # Divide absorbed by total, to get fraction absorbed in subsurface
            if sabg_snl_sum != 0.0
                solarabs.sub_surf_abs_SW_patch[p] /= sabg_snl_sum
            else
                solarabs.sub_surf_abs_SW_patch[p] = 0.0
            end

            # Error handling: redistribute absorbed energy when snow layers changed
            if abs(sabg_snl_sum - solarabs.sabg_snow_patch[p]) > 0.00001
                if col.snl[c] == 0
                    # indices 1:nlevsno → fortran -nlevsno+1:0
                    for jj in 1:nlevsno
                        solarabs.sabg_lyr_patch[p, jj] = 0.0
                    end
                    solarabs.sabg_lyr_patch[p, nlevsno + 1] = solarabs.sabg_patch[p]
                elseif col.snl[c] == -1
                    for jj in 1:(nlevsno - 1)
                        solarabs.sabg_lyr_patch[p, jj] = 0.0
                    end
                    # fortran 0 → nlevsno, fortran 1 → nlevsno+1
                    solarabs.sabg_lyr_patch[p, nlevsno]     = solarabs.sabg_snow_patch[p] * 0.6
                    solarabs.sabg_lyr_patch[p, nlevsno + 1] = solarabs.sabg_snow_patch[p] * 0.4
                else
                    solarabs.sabg_lyr_patch[p, :] .= 0.0
                    # snl(c)+1 in fortran → snl(c)+1+nlevsno in Julia
                    j1 = col.snl[c] + 1 + nlevsno
                    j2 = col.snl[c] + 2 + nlevsno
                    solarabs.sabg_lyr_patch[p, j1] = solarabs.sabg_snow_patch[p] * 0.75
                    solarabs.sabg_lyr_patch[p, j2] = solarabs.sabg_snow_patch[p] * 0.25
                end
            end

            # If shallow snow depth, all solar radiation absorbed in top or top two snow layers
            if !use_subgrid_fluxes || lun.itype[l] == ISTDLAK
                if waterdiag.snow_depth_col[c] < 0.10
                    if col.snl[c] == 0
                        for jj in 1:nlevsno
                            solarabs.sabg_lyr_patch[p, jj] = 0.0
                        end
                        solarabs.sabg_lyr_patch[p, nlevsno + 1] = solarabs.sabg_patch[p]
                    elseif col.snl[c] == -1
                        for jj in 1:(nlevsno - 1)
                            solarabs.sabg_lyr_patch[p, jj] = 0.0
                        end
                        solarabs.sabg_lyr_patch[p, nlevsno]     = solarabs.sabg_patch[p]
                        solarabs.sabg_lyr_patch[p, nlevsno + 1] = 0.0
                    else
                        solarabs.sabg_lyr_patch[p, :] .= 0.0
                        j1 = col.snl[c] + 1 + nlevsno
                        j2 = col.snl[c] + 2 + nlevsno
                        solarabs.sabg_lyr_patch[p, j1] = solarabs.sabg_patch[p] * 0.75
                        solarabs.sabg_lyr_patch[p, j2] = solarabs.sabg_patch[p] * 0.25
                    end
                end
            end
        end

        # Diagnostic: shortwave penetrating ground (e.g. top layer)
        if lun.itype[l] == ISTSOIL || lun.itype[l] == ISTCROP
            # snl(c)+1 in fortran → snl(c)+1+nlevsno in Julia
            j_top = col.snl[c] + 1 + nlevsno
            solarabs.sabg_pen_patch[p] = solarabs.sabg_patch[p] - solarabs.sabg_lyr_patch[p, j_top]
        end

        if use_snicar_frc
            # BC aerosol forcing (patch-level):
            surfrad.sfc_frc_bc_patch[p] = solarabs.sabg_patch[p] - sabg_bc[p]

            # OC aerosol forcing (patch-level):
            if do_sno_oc
                surfrad.sfc_frc_oc_patch[p] = solarabs.sabg_patch[p] - sabg_oc[p]
            else
                surfrad.sfc_frc_oc_patch[p] = 0.0
            end

            # dust aerosol forcing (patch-level):
            surfrad.sfc_frc_dst_patch[p] = solarabs.sabg_patch[p] - sabg_dst[p]

            # all-aerosol forcing (patch-level):
            surfrad.sfc_frc_aer_patch[p] = solarabs.sabg_patch[p] - sabg_pur[p]

            # forcings averaged only over snow:
            if waterdiag.frac_sno_col[c] > 0.0
                surfrad.sfc_frc_bc_sno_patch[p]  = surfrad.sfc_frc_bc_patch[p] / waterdiag.frac_sno_col[c]
                surfrad.sfc_frc_oc_sno_patch[p]  = surfrad.sfc_frc_oc_patch[p] / waterdiag.frac_sno_col[c]
                surfrad.sfc_frc_dst_sno_patch[p] = surfrad.sfc_frc_dst_patch[p] / waterdiag.frac_sno_col[c]
                surfrad.sfc_frc_aer_sno_patch[p] = surfrad.sfc_frc_aer_patch[p] / waterdiag.frac_sno_col[c]
            else
                surfrad.sfc_frc_bc_sno_patch[p]  = SPVAL
                surfrad.sfc_frc_oc_sno_patch[p]  = SPVAL
                surfrad.sfc_frc_dst_sno_patch[p] = SPVAL
                surfrad.sfc_frc_aer_sno_patch[p] = SPVAL
            end
        end
    end

    # --- Radiation diagnostics ---
    for p in bounds
        mask_nourbanp[p] || continue
        g = pch.gridcell[p]
        c = pch.column[p]

        # NDVI and reflected solar radiation
        rvis = surfalb.albd_patch[p, 1] * forc_solad_col[c, 1] +
               surfalb.albi_patch[p, 1] * forc_solai[g, 1]
        rnir = surfalb.albd_patch[p, 2] * forc_solad_col[c, 2] +
               surfalb.albi_patch[p, 2] * forc_solai[g, 2]
        solarabs.fsr_patch[p] = rvis + rnir

        if use_SSRE
            rvisSF = surfalb.albdSF_patch[p, 1] * forc_solad_col[c, 1] +
                     surfalb.albiSF_patch[p, 1] * forc_solai[g, 1]
            rnirSF = surfalb.albdSF_patch[p, 2] * forc_solad_col[c, 2] +
                     surfalb.albiSF_patch[p, 2] * forc_solai[g, 2]
            solarabs.fsrSF_patch[p] = rvisSF + rnirSF
            solarabs.ssre_fsr_patch[p] = solarabs.fsr_patch[p] - solarabs.fsrSF_patch[p]
        end

        surfrad.fsds_vis_d_patch[p] = forc_solad_col[c, 1]
        solarabs.fsds_nir_d_patch[p] = forc_solad_col[c, 2]
        surfrad.fsds_vis_i_patch[p] = forc_solai[g, 1]
        solarabs.fsds_nir_i_patch[p] = forc_solai[g, 2]

        surfrad.fsr_vis_d_patch[p]  = surfalb.albd_patch[p, 1] * forc_solad_col[c, 1]
        solarabs.fsr_nir_d_patch[p] = surfalb.albd_patch[p, 2] * forc_solad_col[c, 2]
        surfrad.fsr_vis_i_patch[p]  = surfalb.albi_patch[p, 1] * forc_solai[g, 1]
        solarabs.fsr_nir_i_patch[p] = surfalb.albi_patch[p, 2] * forc_solai[g, 2]

        if use_SSRE
            surfrad.fsrSF_vis_d_patch[p]  = surfalb.albdSF_patch[p, 1] * forc_solad_col[c, 1]
            solarabs.fsrSF_nir_d_patch[p] = surfalb.albdSF_patch[p, 2] * forc_solad_col[c, 2]
            surfrad.fsrSF_vis_i_patch[p]  = surfalb.albiSF_patch[p, 1] * forc_solai[g, 1]
            solarabs.fsrSF_nir_i_patch[p] = surfalb.albiSF_patch[p, 2] * forc_solai[g, 2]

            surfrad.ssre_fsr_vis_d_patch[p]  = surfrad.fsrSF_vis_d_patch[p] - surfrad.fsr_vis_d_patch[p]
            solarabs.ssre_fsr_nir_d_patch[p] = solarabs.fsrSF_nir_d_patch[p] - solarabs.fsr_nir_d_patch[p]
            surfrad.ssre_fsr_vis_i_patch[p]  = surfrad.fsrSF_vis_i_patch[p] - surfrad.fsr_vis_i_patch[p]
            solarabs.ssre_fsr_nir_i_patch[p] = solarabs.fsrSF_nir_i_patch[p] - solarabs.fsr_nir_i_patch[p]
        end

        deltasec = div(round(Int, dtime), 2)
        if is_near_local_noon(grc.londeg[g]; deltasec=deltasec, current_tod=current_tod)
            surfrad.fsds_vis_d_ln_patch[p] = forc_solad_col[c, 1]
            solarabs.fsds_nir_d_ln_patch[p] = forc_solad_col[c, 2]
            surfrad.fsr_vis_d_ln_patch[p] = surfalb.albd_patch[p, 1] * forc_solad_col[c, 1]
            solarabs.fsr_nir_d_ln_patch[p] = surfalb.albd_patch[p, 2] * forc_solad_col[c, 2]
            surfrad.fsds_vis_i_ln_patch[p] = forc_solai[g, 1]
            surfrad.parveg_ln_patch[p] = parveg[p]
        else
            surfrad.fsds_vis_d_ln_patch[p] = SPVAL
            solarabs.fsds_nir_d_ln_patch[p] = SPVAL
            surfrad.fsr_vis_d_ln_patch[p] = SPVAL
            solarabs.fsr_nir_d_ln_patch[p] = SPVAL
            surfrad.fsds_vis_i_ln_patch[p] = SPVAL
            surfrad.parveg_ln_patch[p] = SPVAL
        end

        if use_SSRE
            if is_near_local_noon(grc.londeg[g]; deltasec=deltasec, current_tod=current_tod)
                surfrad.fsrSF_vis_d_ln_patch[p] = surfalb.albdSF_patch[p, 1] * forc_solad_col[c, 1]
                solarabs.fsrSF_nir_d_ln_patch[p] = surfalb.albdSF_patch[p, 2] * forc_solad_col[c, 2]
            else
                surfrad.fsrSF_vis_d_ln_patch[p] = SPVAL
                solarabs.fsrSF_nir_d_ln_patch[p] = SPVAL
            end
        end

        # Diagnostic variables for history files (OPTIONAL)
        if col.snl[c] < 0
            surfrad.fsds_sno_vd_patch[p] = forc_solad_col[c, 1]
            surfrad.fsds_sno_nd_patch[p] = forc_solad_col[c, 2]
            surfrad.fsds_sno_vi_patch[p] = forc_solai[g, 1]
            surfrad.fsds_sno_ni_patch[p] = forc_solai[g, 2]

            surfrad.fsr_sno_vd_patch[p] = surfrad.fsds_vis_d_patch[p] * surfalb.albsnd_hst_col[c, 1]
            surfrad.fsr_sno_nd_patch[p] = solarabs.fsds_nir_d_patch[p] * surfalb.albsnd_hst_col[c, 2]
            surfrad.fsr_sno_vi_patch[p] = surfrad.fsds_vis_i_patch[p] * surfalb.albsni_hst_col[c, 1]
            surfrad.fsr_sno_ni_patch[p] = solarabs.fsds_nir_i_patch[p] * surfalb.albsni_hst_col[c, 2]
        else
            surfrad.fsds_sno_vd_patch[p] = SPVAL
            surfrad.fsds_sno_nd_patch[p] = SPVAL
            surfrad.fsds_sno_vi_patch[p] = SPVAL
            surfrad.fsds_sno_ni_patch[p] = SPVAL

            surfrad.fsr_sno_vd_patch[p] = SPVAL
            surfrad.fsr_sno_nd_patch[p] = SPVAL
            surfrad.fsr_sno_vi_patch[p] = SPVAL
            surfrad.fsr_sno_ni_patch[p] = SPVAL
        end
    end

    # --- Urban patches ---
    for p in bounds
        mask_urbanp[p] || continue
        g = pch.gridcell[p]
        c = pch.column[p]

        # Solar incident
        surfrad.fsds_vis_d_patch[p] = forc_solad_col[c, 1]
        solarabs.fsds_nir_d_patch[p] = forc_solad_col[c, 2]
        surfrad.fsds_vis_i_patch[p] = forc_solai[g, 1]
        solarabs.fsds_nir_i_patch[p] = forc_solai[g, 2]

        # Determine local noon incident solar
        deltasec = div(round(Int, dtime), 2)
        if is_near_local_noon(grc.londeg[g]; deltasec=deltasec, current_tod=current_tod)
            surfrad.fsds_vis_d_ln_patch[p] = forc_solad_col[c, 1]
            solarabs.fsds_nir_d_ln_patch[p] = forc_solad_col[c, 2]
            surfrad.fsds_vis_i_ln_patch[p] = forc_solai[g, 1]
            surfrad.parveg_ln_patch[p] = 0.0
        else
            surfrad.fsds_vis_d_ln_patch[p] = SPVAL
            solarabs.fsds_nir_d_ln_patch[p] = SPVAL
            surfrad.fsds_vis_i_ln_patch[p] = SPVAL
            surfrad.parveg_ln_patch[p] = SPVAL
        end

        # Solar reflected
        surfrad.fsr_vis_d_patch[p]  = surfalb.albd_patch[p, 1] * forc_solad_col[c, 1]
        solarabs.fsr_nir_d_patch[p] = surfalb.albd_patch[p, 2] * forc_solad_col[c, 2]
        surfrad.fsr_vis_i_patch[p]  = surfalb.albi_patch[p, 1] * forc_solai[g, 1]
        solarabs.fsr_nir_i_patch[p] = surfalb.albi_patch[p, 2] * forc_solai[g, 2]

        # Determine local noon reflected solar
        if is_near_local_noon(grc.londeg[g]; deltasec=deltasec, current_tod=current_tod)
            surfrad.fsr_vis_d_ln_patch[p] = surfrad.fsr_vis_d_patch[p]
            solarabs.fsr_nir_d_ln_patch[p] = solarabs.fsr_nir_d_patch[p]
        else
            surfrad.fsr_vis_d_ln_patch[p] = SPVAL
            solarabs.fsr_nir_d_ln_patch[p] = SPVAL
        end

        solarabs.fsr_patch[p] = surfrad.fsr_vis_d_patch[p] + solarabs.fsr_nir_d_patch[p] +
                                 surfrad.fsr_vis_i_patch[p] + solarabs.fsr_nir_i_patch[p]
    end

    return nothing
end
