# ==========================================================================
# Ported from: src/biogeophys/SurfaceAlbedoType.F90
# Surface albedo state data type allocation and initialization
# ==========================================================================

"""
    SurfaceAlbedoData

Surface albedo state data structure. Holds all albedo variables at gridcell,
column, and patch levels, including ground/snow/soil albedos, canopy radiative
transfer variables, and SNICAR aerosol forcing diagnostics.

Ported from `surfalb_type` in `SurfaceAlbedoType.F90`.
"""
Base.@kwdef mutable struct SurfaceAlbedoData
    # --- Gridcell-level (1D) ---
    azsun_grc                ::Vector{Float64} = Float64[]   # gridcell azimuth angle of sun
    coszen_grc               ::Vector{Float64} = Float64[]   # gridcell cosine of solar zenith angle
    coszen_col               ::Vector{Float64} = Float64[]   # col cosine of solar zenith angle

    # --- Patch-level albedos (2D: npatch × numrad) ---
    albd_patch               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch surface albedo (direct)
    albi_patch               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch surface albedo (diffuse)
    albdSF_patch             ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch snow-free surface albedo (direct)
    albiSF_patch             ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch snow-free surface albedo (diffuse)

    # --- Column-level ground/snow albedos (2D: ncol × numrad) ---
    albgrd_pur_col           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col pure snow ground direct albedo
    albgri_pur_col           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col pure snow ground diffuse albedo
    albgrd_bc_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground direct albedo without BC
    albgri_bc_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground diffuse albedo without BC
    albgrd_oc_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground direct albedo without OC
    albgri_oc_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground diffuse albedo without OC
    albgrd_dst_col           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground direct albedo without dust
    albgri_dst_col           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground diffuse albedo without dust
    albgrd_col               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground albedo (direct)
    albgri_col               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground albedo (diffuse)
    albsod_col               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col soil albedo: direct
    albsoi_col               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col soil albedo: diffuse
    albsnd_hst_col           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col snow albedo, direct, for history files
    albsni_hst_col           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col snow albedo, diffuse, for history files

    # --- SNICAR history output variables (2D: ncol/npatch × numrad) ---
    albd_hst_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch surface albedo (direct) for history
    albi_hst_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch surface albedo (diffuse) for history
    albgrd_pur_hst_col       ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col pure snow ground direct albedo for history
    albgri_pur_hst_col       ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col pure snow ground diffuse albedo for history
    albgrd_bc_hst_col        ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground direct albedo without BC for history
    albgri_bc_hst_col        ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground diffuse albedo without BC for history
    albgrd_oc_hst_col        ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground direct albedo without OC for history
    albgri_oc_hst_col        ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground diffuse albedo without OC for history
    albgrd_dst_hst_col       ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground direct albedo without dust for history
    albgri_dst_hst_col       ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground diffuse albedo without dust for history
    albgrd_hst_col           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground albedo (direct) for history
    albgri_hst_col           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col ground albedo (diffuse) for history
    albsnd_hst2_col          ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col snow albedo, direct, for history (snicar)
    albsni_hst2_col          ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col snow albedo, diffuse, for history (snicar)

    # --- Canopy radiative transfer (patch-level, 2D: npatch × numrad) ---
    ftdd_patch               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch down direct flux below canopy per unit direct flx
    ftid_patch               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch down diffuse flux below canopy per unit direct flx
    ftii_patch               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch down diffuse flux below canopy per unit diffuse flx
    fabd_patch               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch flux absorbed by canopy per unit direct flux
    fabd_sun_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch flux absorbed by sunlit canopy per unit direct flux
    fabd_sha_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch flux absorbed by shaded canopy per unit direct flux
    fabi_patch               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch flux absorbed by canopy per unit diffuse flux
    fabi_sun_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch flux absorbed by sunlit canopy per unit diffuse flux
    fabi_sha_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch flux absorbed by shaded canopy per unit diffuse flux

    # --- Canopy layer PAR absorption (patch-level, 2D: npatch × nlevcan) ---
    fabd_sun_z_patch         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch absorbed sunlit leaf direct PAR per canopy layer
    fabd_sha_z_patch         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch absorbed shaded leaf direct PAR per canopy layer
    fabi_sun_z_patch         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch absorbed sunlit leaf diffuse PAR per canopy layer
    fabi_sha_z_patch         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch absorbed shaded leaf diffuse PAR per canopy layer

    # --- Snow layer flux absorption (col-level, 2D: ncol × (nlevsno+1)) ---
    flx_absdv_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col absorbed flux per unit incident direct flux: VIS
    flx_absdn_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col absorbed flux per unit incident direct flux: NIR
    flx_absiv_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col absorbed flux per unit incident diffuse flux: VIS
    flx_absin_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col absorbed flux per unit incident diffuse flux: NIR

    # --- Canopy structure (patch-level) ---
    fsun_z_patch             ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch sunlit fraction of canopy layer (npatch × nlevcan)
    tlai_z_patch             ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch tlai increment for canopy layer (npatch × nlevcan)
    tsai_z_patch             ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch tsai increment for canopy layer (npatch × nlevcan)
    ncan_patch               ::Vector{Int}     = Int[]       # patch number of canopy layers
    nrad_patch               ::Vector{Int}     = Int[]       # patch number of canopy layers above snow for radiative transfer

    # --- Canopy scaling coefficients (patch-level, 1D) ---
    vcmaxcintsun_patch       ::Vector{Float64} = Float64[]   # patch leaf to canopy scaling coefficient, sunlit leaf vcmax
    vcmaxcintsha_patch       ::Vector{Float64} = Float64[]   # patch leaf to canopy scaling coefficient, shaded leaf vcmax
end

"""
    surfalb_init!(sa::SurfaceAlbedoData, np::Int, nc::Int, ng::Int)

Allocate and initialize all fields of a `SurfaceAlbedoData` instance for
`np` patches, `nc` columns, and `ng` gridcells.
Real fields are initialized to `NaN`, `SPVAL`, or `0.0` as in Fortran.

Ported from `surfalb_type%InitAllocate` in `SurfaceAlbedoType.F90`.
"""
function surfalb_init!(sa::SurfaceAlbedoData, np::Int, nc::Int, ng::Int)
    numrad  = NUMRAD
    nlevcan = NLEVCAN
    nlevsno = varpar.nlevsno

    # Fortran dimension -nlevsno+1:1 → nlevsno+1 levels in Julia
    nlev_sno1 = nlevsno + 1

    # --- Gridcell-level ---
    sa.azsun_grc             = fill(NaN, ng)
    sa.coszen_grc            = fill(NaN, ng)
    sa.coszen_col            = fill(NaN, nc)

    # --- Column-level ground/snow albedos (ncol × numrad) ---
    sa.albgrd_col            = fill(NaN, nc, numrad)
    sa.albgri_col            = fill(NaN, nc, numrad)
    sa.albsnd_hst_col        = fill(SPVAL, nc, numrad)
    sa.albsni_hst_col        = fill(SPVAL, nc, numrad)
    sa.albsod_col            = fill(SPVAL, nc, numrad)
    sa.albsoi_col            = fill(SPVAL, nc, numrad)
    sa.albgrd_pur_col        = fill(NaN, nc, numrad)
    sa.albgri_pur_col        = fill(NaN, nc, numrad)
    sa.albgrd_bc_col         = fill(NaN, nc, numrad)
    sa.albgri_bc_col         = fill(NaN, nc, numrad)
    sa.albgrd_oc_col         = fill(NaN, nc, numrad)
    sa.albgri_oc_col         = fill(NaN, nc, numrad)
    sa.albgrd_dst_col        = fill(NaN, nc, numrad)
    sa.albgri_dst_col        = fill(NaN, nc, numrad)

    # --- Patch-level albedos (npatch × numrad) ---
    sa.albd_patch            = fill(NaN, np, numrad)
    sa.albi_patch            = fill(NaN, np, numrad)
    sa.albdSF_patch          = fill(NaN, np, numrad)
    sa.albiSF_patch          = fill(NaN, np, numrad)

    # --- Canopy radiative transfer (npatch × numrad) ---
    sa.ftdd_patch            = fill(NaN, np, numrad)
    sa.ftid_patch            = fill(NaN, np, numrad)
    sa.ftii_patch            = fill(NaN, np, numrad)
    sa.fabd_patch            = fill(NaN, np, numrad)
    sa.fabd_sun_patch        = fill(NaN, np, numrad)
    sa.fabd_sha_patch        = fill(NaN, np, numrad)
    sa.fabi_patch            = fill(NaN, np, numrad)
    sa.fabi_sun_patch        = fill(NaN, np, numrad)
    sa.fabi_sha_patch        = fill(NaN, np, numrad)

    # --- Canopy layer PAR absorption (npatch × nlevcan) ---
    sa.fabd_sun_z_patch      = fill(0.0, np, nlevcan)
    sa.fabd_sha_z_patch      = fill(0.0, np, nlevcan)
    sa.fabi_sun_z_patch      = fill(0.0, np, nlevcan)
    sa.fabi_sha_z_patch      = fill(0.0, np, nlevcan)

    # --- Snow layer flux absorption (ncol × (nlevsno+1)) ---
    sa.flx_absdv_col         = fill(SPVAL, nc, nlev_sno1)
    sa.flx_absdn_col         = fill(SPVAL, nc, nlev_sno1)
    sa.flx_absiv_col         = fill(SPVAL, nc, nlev_sno1)
    sa.flx_absin_col         = fill(SPVAL, nc, nlev_sno1)

    # --- Canopy structure (npatch × nlevcan) ---
    sa.fsun_z_patch          = fill(0.0, np, nlevcan)
    sa.tlai_z_patch          = fill(0.0, np, nlevcan)
    sa.tsai_z_patch          = fill(0.0, np, nlevcan)
    sa.ncan_patch            = fill(0, np)
    sa.nrad_patch            = fill(0, np)
    sa.vcmaxcintsun_patch    = fill(NaN, np)
    sa.vcmaxcintsha_patch    = fill(NaN, np)

    # --- SNICAR history output variables ---
    sa.albgrd_hst_col        = fill(SPVAL, nc, numrad)
    sa.albgri_hst_col        = fill(SPVAL, nc, numrad)
    sa.albsnd_hst2_col       = fill(SPVAL, nc, numrad)
    sa.albsni_hst2_col       = fill(SPVAL, nc, numrad)
    sa.albgrd_pur_hst_col    = fill(SPVAL, nc, numrad)
    sa.albgri_pur_hst_col    = fill(SPVAL, nc, numrad)
    sa.albgrd_bc_hst_col     = fill(SPVAL, nc, numrad)
    sa.albgri_bc_hst_col     = fill(SPVAL, nc, numrad)
    sa.albgrd_oc_hst_col     = fill(SPVAL, nc, numrad)
    sa.albgri_oc_hst_col     = fill(SPVAL, nc, numrad)
    sa.albgrd_dst_hst_col    = fill(SPVAL, nc, numrad)
    sa.albgri_dst_hst_col    = fill(SPVAL, nc, numrad)
    sa.albd_hst_patch        = fill(SPVAL, np, numrad)
    sa.albi_hst_patch        = fill(SPVAL, np, numrad)

    return nothing
end

"""
    surfalb_clean!(sa::SurfaceAlbedoData)

Deallocate (reset to empty) all fields of a `SurfaceAlbedoData` instance.

Ported from deallocation in `SurfaceAlbedoType.F90`.
"""
function surfalb_clean!(sa::SurfaceAlbedoData)
    # Gridcell-level
    sa.azsun_grc             = Float64[]
    sa.coszen_grc            = Float64[]
    sa.coszen_col            = Float64[]

    # Column-level (2D)
    sa.albgrd_col            = Matrix{Float64}(undef, 0, 0)
    sa.albgri_col            = Matrix{Float64}(undef, 0, 0)
    sa.albsnd_hst_col        = Matrix{Float64}(undef, 0, 0)
    sa.albsni_hst_col        = Matrix{Float64}(undef, 0, 0)
    sa.albsod_col            = Matrix{Float64}(undef, 0, 0)
    sa.albsoi_col            = Matrix{Float64}(undef, 0, 0)
    sa.albgrd_pur_col        = Matrix{Float64}(undef, 0, 0)
    sa.albgri_pur_col        = Matrix{Float64}(undef, 0, 0)
    sa.albgrd_bc_col         = Matrix{Float64}(undef, 0, 0)
    sa.albgri_bc_col         = Matrix{Float64}(undef, 0, 0)
    sa.albgrd_oc_col         = Matrix{Float64}(undef, 0, 0)
    sa.albgri_oc_col         = Matrix{Float64}(undef, 0, 0)
    sa.albgrd_dst_col        = Matrix{Float64}(undef, 0, 0)
    sa.albgri_dst_col        = Matrix{Float64}(undef, 0, 0)

    # Patch-level (2D)
    sa.albd_patch            = Matrix{Float64}(undef, 0, 0)
    sa.albi_patch            = Matrix{Float64}(undef, 0, 0)
    sa.albdSF_patch          = Matrix{Float64}(undef, 0, 0)
    sa.albiSF_patch          = Matrix{Float64}(undef, 0, 0)

    # Canopy radiative transfer
    sa.ftdd_patch            = Matrix{Float64}(undef, 0, 0)
    sa.ftid_patch            = Matrix{Float64}(undef, 0, 0)
    sa.ftii_patch            = Matrix{Float64}(undef, 0, 0)
    sa.fabd_patch            = Matrix{Float64}(undef, 0, 0)
    sa.fabd_sun_patch        = Matrix{Float64}(undef, 0, 0)
    sa.fabd_sha_patch        = Matrix{Float64}(undef, 0, 0)
    sa.fabi_patch            = Matrix{Float64}(undef, 0, 0)
    sa.fabi_sun_patch        = Matrix{Float64}(undef, 0, 0)
    sa.fabi_sha_patch        = Matrix{Float64}(undef, 0, 0)

    # Canopy layer PAR
    sa.fabd_sun_z_patch      = Matrix{Float64}(undef, 0, 0)
    sa.fabd_sha_z_patch      = Matrix{Float64}(undef, 0, 0)
    sa.fabi_sun_z_patch      = Matrix{Float64}(undef, 0, 0)
    sa.fabi_sha_z_patch      = Matrix{Float64}(undef, 0, 0)

    # Snow layer flux absorption
    sa.flx_absdv_col         = Matrix{Float64}(undef, 0, 0)
    sa.flx_absdn_col         = Matrix{Float64}(undef, 0, 0)
    sa.flx_absiv_col         = Matrix{Float64}(undef, 0, 0)
    sa.flx_absin_col         = Matrix{Float64}(undef, 0, 0)

    # Canopy structure
    sa.fsun_z_patch          = Matrix{Float64}(undef, 0, 0)
    sa.tlai_z_patch          = Matrix{Float64}(undef, 0, 0)
    sa.tsai_z_patch          = Matrix{Float64}(undef, 0, 0)
    sa.ncan_patch            = Int[]
    sa.nrad_patch            = Int[]
    sa.vcmaxcintsun_patch    = Float64[]
    sa.vcmaxcintsha_patch    = Float64[]

    # SNICAR history output
    sa.albgrd_hst_col        = Matrix{Float64}(undef, 0, 0)
    sa.albgri_hst_col        = Matrix{Float64}(undef, 0, 0)
    sa.albsnd_hst2_col       = Matrix{Float64}(undef, 0, 0)
    sa.albsni_hst2_col       = Matrix{Float64}(undef, 0, 0)
    sa.albgrd_pur_hst_col    = Matrix{Float64}(undef, 0, 0)
    sa.albgri_pur_hst_col    = Matrix{Float64}(undef, 0, 0)
    sa.albgrd_bc_hst_col     = Matrix{Float64}(undef, 0, 0)
    sa.albgri_bc_hst_col     = Matrix{Float64}(undef, 0, 0)
    sa.albgrd_oc_hst_col     = Matrix{Float64}(undef, 0, 0)
    sa.albgri_oc_hst_col     = Matrix{Float64}(undef, 0, 0)
    sa.albgrd_dst_hst_col    = Matrix{Float64}(undef, 0, 0)
    sa.albgri_dst_hst_col    = Matrix{Float64}(undef, 0, 0)
    sa.albd_hst_patch        = Matrix{Float64}(undef, 0, 0)
    sa.albi_hst_patch        = Matrix{Float64}(undef, 0, 0)

    return nothing
end

"""
    surfalb_init_cold!(sa::SurfaceAlbedoData, bounds_col::UnitRange{Int},
                       bounds_patch::UnitRange{Int};
                       use_SSRE::Bool=false)

Initialize cold-start conditions for surface albedo state variables.
Sets reasonable default albedo values for ground, snow, soil, and canopy.

Ported from `surfalb_type%InitCold` in `SurfaceAlbedoType.F90`.
"""
function surfalb_init_cold!(sa::SurfaceAlbedoData,
                            bounds_col::UnitRange{Int},
                            bounds_patch::UnitRange{Int};
                            use_SSRE::Bool = false)
    numrad = NUMRAD

    # Column-level ground/snow/soil albedos
    for c in bounds_col
        for ib in 1:numrad
            sa.albgrd_col[c, ib]     = 0.2
            sa.albgri_col[c, ib]     = 0.2
            sa.albsod_col[c, ib]     = 0.2
            sa.albsoi_col[c, ib]     = 0.2
            sa.albsnd_hst_col[c, ib] = 0.6
            sa.albsni_hst_col[c, ib] = 0.6
            sa.albgrd_pur_col[c, ib] = 0.2
            sa.albgri_pur_col[c, ib] = 0.2
            sa.albgrd_bc_col[c, ib]  = 0.2
            sa.albgri_bc_col[c, ib]  = 0.2
            sa.albgrd_oc_col[c, ib]  = 0.2
            sa.albgri_oc_col[c, ib]  = 0.2
            sa.albgrd_dst_col[c, ib] = 0.2
            sa.albgri_dst_col[c, ib] = 0.2
        end
    end

    # Patch-level albedos and canopy fluxes
    for p in bounds_patch
        for ib in 1:numrad
            sa.albd_patch[p, ib]     = 0.2
            sa.albi_patch[p, ib]     = 0.2

            sa.fabi_patch[p, ib]     = 0.0
            sa.fabd_patch[p, ib]     = 0.0
            sa.fabi_sun_patch[p, ib] = 0.0
            sa.fabd_sun_patch[p, ib] = 0.0
            sa.fabd_sha_patch[p, ib] = 0.0
            sa.fabi_sha_patch[p, ib] = 0.0
            sa.ftdd_patch[p, ib]     = 1.0
            sa.ftid_patch[p, ib]     = 0.0
            sa.ftii_patch[p, ib]     = 1.0
        end
    end

    # Snow-free albedos (SSRE diagnostic)
    if use_SSRE
        for p in bounds_patch
            for ib in 1:numrad
                sa.albdSF_patch[p, ib] = 0.2
                sa.albiSF_patch[p, ib] = 0.2
            end
        end
    end

    return nothing
end

# ==========================================================================
# The following subroutines depend on infrastructure modules that are not yet
# ported (history, restart/IO). They are provided as stubs that document the
# Fortran interface and can be filled in when those modules become available.
# ==========================================================================

"""
    surfalb_init_history!(sa::SurfaceAlbedoData, bounds_col::UnitRange{Int},
                          bounds_patch::UnitRange{Int}, bounds_grc::UnitRange{Int})

Register surface albedo fields for history file output.

Ported from `surfalb_type%InitHistory` in `SurfaceAlbedoType.F90`.
Requires history infrastructure (histFileMod) -- stub until that module is ported.
"""
function surfalb_init_history!(sa::SurfaceAlbedoData,
                               bounds_col::UnitRange{Int},
                               bounds_patch::UnitRange{Int},
                               bounds_grc::UnitRange{Int})
    # Set fields to SPVAL as Fortran does in InitHistory
    for g in bounds_grc
        sa.azsun_grc[g]  = SPVAL
        sa.coszen_grc[g] = SPVAL
    end
    for c in bounds_col
        sa.coszen_col[c] = SPVAL
        for ib in 1:NUMRAD
            sa.albgrd_col[c, ib] = SPVAL
            sa.albgri_col[c, ib] = SPVAL
        end
    end
    for p in bounds_patch
        for ib in 1:NUMRAD
            sa.albd_patch[p, ib] = SPVAL
            sa.albi_patch[p, ib] = SPVAL
        end
    end
    # Stub: history field registration will be added when histFileMod is ported.
    # Fields that would be registered:
    #   AZSUN, COSZEN_GRC, COSZEN, ALBGRD, ALBGRI, ALBDSF, ALBISF, ALBD, ALBI,
    #   ALBD_HIST, ALBI_HIST, ALBGRD_HIST, ALBGRI_HIST, ALBGRD_PUR_HIST, ALBGRI_PUR_HIST,
    #   ALBGRD_BC_HIST, ALBGRI_BC_HIST, ALBGRD_OC_HIST, ALBGRI_OC_HIST,
    #   ALBGRD_DST_HIST, ALBGRI_DST_HIST, ALBSND_HIST, ALBSNI_HIST
    return nothing
end

"""
    surfalb_restart!(sa::SurfaceAlbedoData, bounds_col::UnitRange{Int},
                     bounds_patch::UnitRange{Int})

Read/write surface albedo state from/to restart file.

Ported from `surfalb_type%Restart` in `SurfaceAlbedoType.F90`.
Requires NetCDF/restart infrastructure -- stub until that module is ported.
"""
function surfalb_restart!(sa::SurfaceAlbedoData,
                          bounds_col::UnitRange{Int},
                          bounds_patch::UnitRange{Int})
    # Stub: restart variable I/O will be added when restUtilMod/ncdio_pio is ported.
    # Variables that would be read/written:
    #   coszen_grc, coszen, albd, albi, albdSF, albiSF, albgrd, albgri,
    #   albsod, albsoi, albsnd_hst, albsni_hst, tlai_z, tsai_z, ncan, nrad,
    #   fsun_z, vcmaxcintsun, vcmaxcintsha, albgrd_bc, albgri_bc, albgrd_pur,
    #   albgri_pur, albgrd_oc, albgri_oc, albgrd_dst, albgri_dst,
    #   albd_hist, albi_hist, albgrd_hist, albgri_hist, albsnd_hst2, albsni_hst2,
    #   albgrd_bc_hist, albgri_bc_hist, albgrd_pur_hist, albgri_pur_hist,
    #   albgrd_oc_hist, albgri_oc_hist, albgrd_dst_hist, albgri_dst_hist,
    #   fabd, fabi, fabd_sun, fabd_sha, fabi_sun, fabi_sha,
    #   fabd_sun_z, fabd_sha_z, fabi_sun_z, fabi_sha_z,
    #   ftdd, ftid, ftii, flx_absdv, flx_absdn, flx_absiv, flx_absin
    return nothing
end
