# ==========================================================================
# Ported from: src/biogeophys/SolarAbsorbedType.F90
# Solar absorption state data type allocation and initialization
# ==========================================================================

"""
    SolarAbsorbedData

Solar absorption state data structure. Holds all solar absorbed/reflected
variables at patch and landunit levels, including vegetation absorption,
ground absorption, urban surface absorption, and NIR diagnostics.

Ported from `solarabs_type` in `SolarAbsorbedType.F90`.
"""
Base.@kwdef mutable struct SolarAbsorbedData{FT<:AbstractFloat}
    # --- Solar reflected (patch-level 1D) ---
    fsr_patch                ::Vector{FT} = Float64[]   # patch solar radiation reflected (W/m**2)
    fsrSF_patch              ::Vector{FT} = Float64[]   # diagnostic snow-free patch solar radiation reflected (W/m**2)
    ssre_fsr_patch           ::Vector{FT} = Float64[]   # snow radiative effect on patch solar radiation reflected (W/m**2)

    # --- Solar absorbed (patch-level 1D) ---
    fsa_patch                ::Vector{FT} = Float64[]   # patch solar radiation absorbed (total) (W/m**2)
    fsa_u_patch              ::Vector{FT} = Float64[]   # patch urban solar radiation absorbed (total) (W/m**2)
    fsa_r_patch              ::Vector{FT} = Float64[]   # patch rural solar radiation absorbed (total) (W/m**2)

    # --- Canopy PAR absorption (patch-level 2D: npatch × nlevcan) ---
    parsun_z_patch           ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # patch absorbed PAR for sunlit leaves in canopy layer (W/m**2)
    parsha_z_patch           ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # patch absorbed PAR for shaded leaves in canopy layer (W/m**2)
    par240d_z_patch          ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # 10-day running mean of daytime patch absorbed PAR (W/m**2)
    par240x_z_patch          ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # 10-day running mean of maximum patch absorbed PAR (W/m**2)
    par24d_z_patch           ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # daily accumulated absorbed PAR (J/m**2)
    par24x_z_patch           ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # daily max of patch absorbed PAR (W/m**2)

    # --- Ground absorption (patch-level 1D) ---
    sabg_soil_patch          ::Vector{FT} = Float64[]   # patch solar radiation absorbed by soil (W/m**2)
    sabg_snow_patch          ::Vector{FT} = Float64[]   # patch solar radiation absorbed by snow (W/m**2)
    sabg_patch               ::Vector{FT} = Float64[]   # patch solar radiation absorbed by ground (W/m**2)
    sabg_chk_patch           ::Vector{FT} = Float64[]   # patch fsno weighted sum (W/m**2)
    sabg_pen_patch           ::Vector{FT} = Float64[]   # patch shortwave radiation penetrating top soisno layer (W/m2)
    sub_surf_abs_SW_patch    ::Vector{FT} = Float64[]   # patch fraction of solar radiation absorbed below first snow layer

    # --- Ground absorption per layer (patch-level 2D: npatch × (nlevsno+1)) ---
    sabg_lyr_patch           ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # patch absorbed radiation in each snow layer and top soil layer (W/m2)

    # --- Vegetation absorption (patch-level 1D) ---
    sabv_patch               ::Vector{FT} = Float64[]   # patch solar radiation absorbed by vegetation (W/m**2)

    # --- Urban surface absorption (landunit-level 2D: nlun × numrad) ---
    sabs_roof_dir_lun        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun direct solar absorbed by roof
    sabs_roof_dif_lun        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun diffuse solar absorbed by roof
    sabs_sunwall_dir_lun     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun direct solar absorbed by sunwall
    sabs_sunwall_dif_lun     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun diffuse solar absorbed by sunwall
    sabs_shadewall_dir_lun   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun direct solar absorbed by shadewall
    sabs_shadewall_dif_lun   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun diffuse solar absorbed by shadewall
    sabs_improad_dir_lun     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun direct solar absorbed by impervious road
    sabs_improad_dif_lun     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun diffuse solar absorbed by impervious road
    sabs_perroad_dir_lun     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun direct solar absorbed by pervious road
    sabs_perroad_dif_lun     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # lun diffuse solar absorbed by pervious road

    # --- NIR diagnostics (patch-level 1D) ---
    fsds_nir_d_patch         ::Vector{FT} = Float64[]   # patch incident direct beam nir solar radiation (W/m**2)
    fsds_nir_i_patch         ::Vector{FT} = Float64[]   # patch incident diffuse nir solar radiation (W/m**2)
    fsds_nir_d_ln_patch      ::Vector{FT} = Float64[]   # patch incident direct beam nir solar radiation at local noon (W/m**2)
    fsr_nir_d_patch          ::Vector{FT} = Float64[]   # patch reflected direct beam nir solar radiation (W/m**2)
    fsr_nir_i_patch          ::Vector{FT} = Float64[]   # patch reflected diffuse nir solar radiation (W/m**2)
    fsr_nir_d_ln_patch       ::Vector{FT} = Float64[]   # patch reflected direct beam nir solar radiation at local noon (W/m**2)

    # --- Snow-free NIR diagnostics (patch-level 1D) ---
    fsrSF_nir_d_patch        ::Vector{FT} = Float64[]   # snow-free patch reflected direct beam nir solar radiation (W/m**2)
    fsrSF_nir_i_patch        ::Vector{FT} = Float64[]   # snow-free patch reflected diffuse nir solar radiation (W/m**2)
    fsrSF_nir_d_ln_patch     ::Vector{FT} = Float64[]   # snow-free patch reflected direct beam nir solar radiation at local noon (W/m**2)
    ssre_fsr_nir_d_patch     ::Vector{FT} = Float64[]   # snow radiative effect on direct beam nir reflected (W/m**2)
    ssre_fsr_nir_i_patch     ::Vector{FT} = Float64[]   # snow radiative effect on diffuse nir reflected (W/m**2)
    ssre_fsr_nir_d_ln_patch  ::Vector{FT} = Float64[]   # snow radiative effect on direct beam nir reflected at local noon (W/m**2)
end

"""
    solarabs_init!(sa::SolarAbsorbedData, np::Int, nl::Int;
                   use_luna::Bool=false)

Allocate and initialize all fields of a `SolarAbsorbedData` instance for
`np` patches and `nl` landunits.
Real fields are initialized to `NaN` or `SPVAL` as in Fortran InitAllocate.

Ported from `solarabs_type%InitAllocate` in `SolarAbsorbedType.F90`.
"""
function solarabs_init!(sa::SolarAbsorbedData, np::Int, nl::Int;
                        use_luna::Bool = false)
    numrad  = NUMRAD
    nlevcan = NLEVCAN
    nlevsno = varpar.nlevsno

    # Fortran dimension -nlevsno+1:1 → nlevsno+1 levels in Julia
    nlev_sno1 = nlevsno + 1

    # --- Patch-level 1D (NaN init) ---
    sa.fsa_patch              = fill(NaN, np)
    sa.fsa_u_patch            = fill(NaN, np)
    sa.fsa_r_patch            = fill(NaN, np)
    sa.sabv_patch             = fill(NaN, np)
    sa.sabg_patch             = fill(NaN, np)
    sa.sabg_pen_patch         = fill(NaN, np)
    sa.sabg_soil_patch        = fill(NaN, np)
    sa.sabg_snow_patch        = fill(NaN, np)
    sa.sabg_chk_patch         = fill(NaN, np)
    sa.sub_surf_abs_SW_patch  = fill(NaN, np)
    sa.fsr_patch              = fill(NaN, np)
    sa.fsrSF_patch            = fill(NaN, np)
    sa.ssre_fsr_patch         = fill(NaN, np)

    # --- Patch-level 2D: canopy PAR (NaN init) ---
    sa.parsun_z_patch         = fill(NaN, np, nlevcan)
    sa.parsha_z_patch         = fill(NaN, np, nlevcan)

    # --- LUNA-dependent PAR arrays (SPVAL init) ---
    if use_luna
        sa.par240d_z_patch    = fill(SPVAL, np, nlevcan)
        sa.par240x_z_patch    = fill(SPVAL, np, nlevcan)
        sa.par24d_z_patch     = fill(SPVAL, np, nlevcan)
        sa.par24x_z_patch     = fill(SPVAL, np, nlevcan)
    else
        sa.par240d_z_patch    = Matrix{Float64}(undef, 0, 0)
        sa.par240x_z_patch    = Matrix{Float64}(undef, 0, 0)
        sa.par24d_z_patch     = Matrix{Float64}(undef, 0, 0)
        sa.par24x_z_patch     = Matrix{Float64}(undef, 0, 0)
    end

    # --- Patch-level 2D: snow layer absorption (NaN init) ---
    sa.sabg_lyr_patch         = fill(NaN, np, nlev_sno1)

    # --- Landunit-level 2D: urban surface absorption (NaN init) ---
    sa.sabs_roof_dir_lun      = fill(NaN, nl, numrad)
    sa.sabs_roof_dif_lun      = fill(NaN, nl, numrad)
    sa.sabs_sunwall_dir_lun   = fill(NaN, nl, numrad)
    sa.sabs_sunwall_dif_lun   = fill(NaN, nl, numrad)
    sa.sabs_shadewall_dir_lun = fill(NaN, nl, numrad)
    sa.sabs_shadewall_dif_lun = fill(NaN, nl, numrad)
    sa.sabs_improad_dir_lun   = fill(NaN, nl, numrad)
    sa.sabs_improad_dif_lun   = fill(NaN, nl, numrad)
    sa.sabs_perroad_dir_lun   = fill(NaN, nl, numrad)
    sa.sabs_perroad_dif_lun   = fill(NaN, nl, numrad)

    # --- NIR diagnostics (NaN init) ---
    sa.fsr_nir_d_patch        = fill(NaN, np)
    sa.fsr_nir_i_patch        = fill(NaN, np)
    sa.fsr_nir_d_ln_patch     = fill(NaN, np)
    sa.fsds_nir_d_patch       = fill(NaN, np)
    sa.fsds_nir_i_patch       = fill(NaN, np)
    sa.fsds_nir_d_ln_patch    = fill(NaN, np)

    # --- Snow-free NIR diagnostics (NaN init) ---
    sa.fsrSF_nir_d_patch      = fill(NaN, np)
    sa.fsrSF_nir_i_patch      = fill(NaN, np)
    sa.fsrSF_nir_d_ln_patch   = fill(NaN, np)
    sa.ssre_fsr_nir_d_patch   = fill(NaN, np)
    sa.ssre_fsr_nir_i_patch   = fill(NaN, np)
    sa.ssre_fsr_nir_d_ln_patch = fill(NaN, np)

    return nothing
end

"""
    solarabs_clean!(sa::SolarAbsorbedData)

Deallocate (reset to empty) all fields of a `SolarAbsorbedData` instance.

Ported from deallocation in `SolarAbsorbedType.F90`.
"""
function solarabs_clean!(sa::SolarAbsorbedData{FT}) where {FT}
    # Patch-level 1D
    sa.fsr_patch              = FT[]
    sa.fsrSF_patch            = FT[]
    sa.ssre_fsr_patch         = FT[]
    sa.fsa_patch              = FT[]
    sa.fsa_u_patch            = FT[]
    sa.fsa_r_patch            = FT[]
    sa.sabg_soil_patch        = FT[]
    sa.sabg_snow_patch        = FT[]
    sa.sabg_patch             = FT[]
    sa.sabg_chk_patch         = FT[]
    sa.sabg_pen_patch         = FT[]
    sa.sub_surf_abs_SW_patch  = FT[]
    sa.sabv_patch             = FT[]

    # Patch-level 2D
    sa.parsun_z_patch         = Matrix{FT}(undef, 0, 0)
    sa.parsha_z_patch         = Matrix{FT}(undef, 0, 0)
    sa.par240d_z_patch        = Matrix{FT}(undef, 0, 0)
    sa.par240x_z_patch        = Matrix{FT}(undef, 0, 0)
    sa.par24d_z_patch         = Matrix{FT}(undef, 0, 0)
    sa.par24x_z_patch         = Matrix{FT}(undef, 0, 0)
    sa.sabg_lyr_patch         = Matrix{FT}(undef, 0, 0)

    # Landunit-level 2D
    sa.sabs_roof_dir_lun      = Matrix{FT}(undef, 0, 0)
    sa.sabs_roof_dif_lun      = Matrix{FT}(undef, 0, 0)
    sa.sabs_sunwall_dir_lun   = Matrix{FT}(undef, 0, 0)
    sa.sabs_sunwall_dif_lun   = Matrix{FT}(undef, 0, 0)
    sa.sabs_shadewall_dir_lun = Matrix{FT}(undef, 0, 0)
    sa.sabs_shadewall_dif_lun = Matrix{FT}(undef, 0, 0)
    sa.sabs_improad_dir_lun   = Matrix{FT}(undef, 0, 0)
    sa.sabs_improad_dif_lun   = Matrix{FT}(undef, 0, 0)
    sa.sabs_perroad_dir_lun   = Matrix{FT}(undef, 0, 0)
    sa.sabs_perroad_dif_lun   = Matrix{FT}(undef, 0, 0)

    # NIR diagnostics
    sa.fsds_nir_d_patch       = FT[]
    sa.fsds_nir_i_patch       = FT[]
    sa.fsds_nir_d_ln_patch    = FT[]
    sa.fsr_nir_d_patch        = FT[]
    sa.fsr_nir_i_patch        = FT[]
    sa.fsr_nir_d_ln_patch     = FT[]
    sa.fsrSF_nir_d_patch      = FT[]
    sa.fsrSF_nir_i_patch      = FT[]
    sa.fsrSF_nir_d_ln_patch   = FT[]
    sa.ssre_fsr_nir_d_patch   = FT[]
    sa.ssre_fsr_nir_i_patch   = FT[]
    sa.ssre_fsr_nir_d_ln_patch = FT[]

    return nothing
end

"""
    solarabs_init_cold!(sa::SolarAbsorbedData, bounds_lun::UnitRange{Int})

Initialize cold-start conditions for solar absorption state variables.
Sets urban surface absorption fields to zero.

Ported from `solarabs_type%InitCold` in `SolarAbsorbedType.F90`.
"""
function solarabs_init_cold!(sa::SolarAbsorbedData{FT}, bounds_lun::UnitRange{Int}) where {FT}
    for l in bounds_lun
        for ib in 1:NUMRAD
            sa.sabs_roof_dir_lun[l, ib]      = 0.0
            sa.sabs_roof_dif_lun[l, ib]      = 0.0
            sa.sabs_sunwall_dir_lun[l, ib]   = 0.0
            sa.sabs_sunwall_dif_lun[l, ib]   = 0.0
            sa.sabs_shadewall_dir_lun[l, ib] = 0.0
            sa.sabs_shadewall_dif_lun[l, ib] = 0.0
            sa.sabs_improad_dir_lun[l, ib]   = 0.0
            sa.sabs_improad_dif_lun[l, ib]   = 0.0
            sa.sabs_perroad_dir_lun[l, ib]   = 0.0
            sa.sabs_perroad_dif_lun[l, ib]   = 0.0
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
    solarabs_init_history!(sa::SolarAbsorbedData, bounds_patch::UnitRange{Int};
                           use_SSRE::Bool=false, use_luna::Bool=false)

Register solar absorption fields for history file output.
Sets relevant fields to SPVAL as Fortran does in InitHistory.

Ported from `solarabs_type%InitHistory` in `SolarAbsorbedType.F90`.
Requires history infrastructure (histFileMod) -- stub until that module is ported.
"""
function solarabs_init_history!(sa::SolarAbsorbedData,
                                 bounds_patch::UnitRange{Int};
                                 use_SSRE::Bool = false,
                                 use_luna::Bool = false)
    # Set fields to SPVAL as Fortran does in InitHistory
    for p in bounds_patch
        sa.fsa_patch[p]              = SPVAL
        sa.fsa_r_patch[p]           = SPVAL
        sa.fsa_u_patch[p]           = SPVAL
        sa.fsr_patch[p]             = SPVAL
        sa.sabv_patch[p]            = SPVAL
        sa.sabg_patch[p]            = SPVAL
        sa.sabg_pen_patch[p]        = SPVAL
        sa.fsds_nir_d_patch[p]      = SPVAL
        sa.fsds_nir_i_patch[p]      = SPVAL
        sa.fsds_nir_d_ln_patch[p]   = SPVAL
        sa.fsr_nir_d_patch[p]       = SPVAL
        sa.fsr_nir_i_patch[p]       = SPVAL
        sa.fsr_nir_d_ln_patch[p]    = SPVAL
        sa.sub_surf_abs_SW_patch[p] = SPVAL
    end

    # Set sabg_lyr_patch snow layers to SPVAL (Fortran sets -nlevsno+1:0)
    nlevsno = varpar.nlevsno
    for p in bounds_patch
        for j in 1:nlevsno
            sa.sabg_lyr_patch[p, j] = SPVAL
        end
    end

    if use_SSRE
        for p in bounds_patch
            sa.fsrSF_patch[p]            = SPVAL
            sa.ssre_fsr_patch[p]         = SPVAL
            sa.fsrSF_nir_d_patch[p]      = SPVAL
            sa.fsrSF_nir_i_patch[p]      = SPVAL
            sa.fsrSF_nir_d_ln_patch[p]   = SPVAL
            sa.ssre_fsr_nir_d_patch[p]   = SPVAL
            sa.ssre_fsr_nir_i_patch[p]   = SPVAL
            sa.ssre_fsr_nir_d_ln_patch[p] = SPVAL
        end
    end

    # Stub: history field registration will be added when histFileMod is ported.
    # Fields that would be registered:
    #   FSA, FSA_ICE, FSA_R, FSA_U, FSR, SWup, FSR_ICE, SNO_ABS, SNO_ABS_ICE,
    #   SABV, SABG, SABG_PEN, FSDSND, FSDSNI, FSDSNDLN, FSRND, FSRNI, FSRNDLN,
    #   FSRSF, SSRE_FSR, FSRSFND, FSRSFNI, FSRSFNDLN, SSRE_FSRND, SSRE_FSRNI,
    #   SSRE_FSRNDLN, SNOINTABS, PAR240DZ, PAR240XZ
    return nothing
end

"""
    solarabs_restart!(sa::SolarAbsorbedData, bounds_lun::UnitRange{Int},
                      bounds_patch::UnitRange{Int})

Read/write solar absorption state from/to restart file.

Ported from `solarabs_type%Restart` in `SolarAbsorbedType.F90`.
Requires NetCDF/restart infrastructure -- stub until that module is ported.
"""
function solarabs_restart!(sa::SolarAbsorbedData,
                            bounds_lun::UnitRange{Int},
                            bounds_patch::UnitRange{Int})
    # Stub: restart variable I/O will be added when restUtilMod/ncdio_pio is ported.
    # Variables that would be read/written:
    #   sabs_roof_dir, sabs_roof_dif, sabs_sunwall_dir, sabs_sunwall_dif,
    #   sabs_shadewall_dir, sabs_shadewall_dif, sabs_improad_dir, sabs_improad_dif,
    #   sabs_perroad_dir, sabs_perroad_dif,
    #   par240d, par24d, par240x, par24x, parsun, parsha (LUNA only)
    return nothing
end
