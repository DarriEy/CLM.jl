# ==========================================================================
# Ported from: src/biogeophys/SoilStateType.F90
# Soil state data type allocation and initialization
# ==========================================================================

"""
    SoilStateData

Soil state data structure. Holds all soil state variables at patch and
column levels, including sand/clay fractions, hydraulic properties,
thermal conductivity/heat capacity, van Genuchten parameters, and
root distributions.

Ported from `soilstate_type` in `SoilStateType.F90`.
"""
Base.@kwdef mutable struct SoilStateData{FT<:Real}
    # --- Sand / clay / organic matter ---
    sandfrac_patch           ::Vector{FT} = Float64[]   # patch sand fraction
    clayfrac_patch           ::Vector{FT} = Float64[]   # patch clay fraction
    mss_frc_cly_vld_col      ::Vector{FT} = Float64[]   # col mass fraction clay limited to 0.20
    cellorg_col              ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col organic matter (ncols, nlevsoi)
    cellsand_col             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col sand value (ncols, nlevsoi)
    cellclay_col             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col clay value (ncols, nlevsoi)
    bd_col                   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col bulk density of dry soil [kg/m^3] (ncols, nlevgrnd)

    # --- Hydraulic properties ---
    hksat_col                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col hydraulic conductivity at saturation (mm H2O/s) (ncols, nlevgrnd)
    hksat_min_col            ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col mineral hydraulic conductivity at saturation (mm/s) (ncols, nlevgrnd)
    hk_l_col                 ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col hydraulic conductivity (mm/s) (ncols, nlevgrnd)
    smp_l_col                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col soil matric potential (mm) (ncols, nlevgrnd)
    smpmin_col               ::Vector{FT} = Float64[]   # col restriction for min of soil potential (mm)
    bsw_col                  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col Clapp and Hornberger "b" (ncols, nlevgrnd)
    watsat_col               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col volumetric soil water at saturation (porosity) (ncols, nlevmaxurbgrnd)
    watdry_col               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col btran parameter for btran=0 (ncols, nlevgrnd)
    watopt_col               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col btran parameter for btran=1 (ncols, nlevgrnd)
    watfc_col                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col volumetric soil water at field capacity (ncols, nlevgrnd)
    sucsat_col               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col minimum soil suction (mm) (ncols, nlevgrnd)
    dsl_col                  ::Vector{FT} = Float64[]   # col dry surface layer thickness (mm)
    soilresis_col            ::Vector{FT} = Float64[]   # col soil evaporative resistance S&L14 (s/m)
    soilbeta_col             ::Vector{FT} = Float64[]   # col factor that reduces ground evaporation L&P1992 (-)
    soilalpha_col            ::Vector{FT} = Float64[]   # col factor that reduces ground saturated specific humidity (-)
    soilalpha_u_col          ::Vector{FT} = Float64[]   # col urban factor that reduces ground saturated specific humidity (-)
    soilpsi_col              ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col soil water potential in each soil layer (MPa) (ncols, nlevgrnd)
    wtfact_col               ::Vector{FT} = Float64[]   # col maximum saturated fraction for a gridcell
    porosity_col             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col soil porosity (1-bulk_density/soil_density) (VIC) (ncols, nlayer)
    eff_porosity_col         ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col effective porosity = porosity - vol_ice (ncols, nlevgrnd)
    gwc_thr_col              ::Vector{FT} = Float64[]   # col threshold soil moisture based on clay content

    # --- Van Genuchten parameters ---
    msw_col                  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col vanGenuchtenClapp "m" (ncols, nlevgrnd)
    nsw_col                  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col vanGenuchtenClapp "n" (ncols, nlevgrnd)
    alphasw_col              ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col vanGenuchtenClapp "alpha" (ncols, nlevgrnd)
    watres_col               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col residual soil water content (ncols, nlevgrnd)

    # --- Thermal conductivity / heat capacity ---
    thk_col                  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col thermal conductivity [W/m-K] (ncols, nlevsno+nlevmaxurbgrnd)
    tkmg_col                 ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col thermal conductivity, soil minerals [W/m-K] (ncols, nlevgrnd)
    tkdry_col                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col thermal conductivity, dry soil (W/m/K) (ncols, nlevgrnd)
    tksatu_col               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col thermal conductivity, saturated soil [W/m-K] (ncols, nlevgrnd)
    csol_col                 ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col heat capacity, soil solids (J/m^3/K) (ncols, nlevgrnd)

    # --- Roots ---
    rootr_patch              ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # patch effective fraction of roots in each soil layer (SMS) (npatches, nlevgrnd)
    rootr_col                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col effective fraction of roots in each soil layer (SMS) (ncols, nlevgrnd)
    rootfr_col               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col fraction of roots in each soil layer (ncols, nlevgrnd)
    rootfr_patch             ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # patch fraction of roots for water (npatches, nlevgrnd)
    crootfr_patch            ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # patch fraction of roots for carbon (npatches, nlevgrnd)
    root_depth_patch         ::Vector{FT} = Float64[]   # patch root depth
    rootr_road_perv_col      ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col effective root fraction in urban pervious road (ncols, nlevgrnd)
    rootfr_road_perv_col     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # col root fraction in urban pervious road (ncols, nlevgrnd)
    k_soil_root_patch        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # patch soil-root interface conductance [mm/s] (npatches, nlevsoi)
    root_conductance_patch   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # patch root conductance [mm/s] (npatches, nlevsoi)
    soil_conductance_patch   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # patch soil conductance [mm/s] (npatches, nlevsoi)
end

"""
    soilstate_init!(ss::SoilStateData, np::Int, nc::Int)

Allocate and initialize all fields of a `SoilStateData` instance for
`np` patches and `nc` columns.
Real fields are initialized to `NaN` or `SPVAL` matching the Fortran source.

Ported from `soilstate_type%InitAllocate` in `SoilStateType.F90`.
"""
function soilstate_init!(ss::SoilStateData{FT}, np::Int, nc::Int) where {FT}
    nlevsoi        = varpar.nlevsoi
    nlevgrnd       = varpar.nlevgrnd
    nlevsno        = varpar.nlevsno
    nlevmaxurbgrnd = varpar.nlevmaxurbgrnd

    # Vertical dimension size for thk_col (Fortran offset -nlevsno+1:nlevmaxurbgrnd)
    nlev_soisno = nlevsno + nlevmaxurbgrnd

    # --- Sand / clay / organic matter ---
    ss.sandfrac_patch        = fill(FT(NaN), np)
    ss.clayfrac_patch        = fill(FT(NaN), np)
    ss.mss_frc_cly_vld_col   = fill(FT(NaN), nc)
    ss.cellorg_col           = fill(FT(NaN), nc, nlevsoi)
    ss.cellsand_col          = fill(FT(NaN), nc, nlevsoi)
    ss.cellclay_col          = fill(FT(NaN), nc, nlevsoi)
    ss.bd_col                = fill(FT(NaN), nc, nlevgrnd)

    # --- Hydraulic properties ---
    ss.hksat_col             = fill(FT(SPVAL), nc, nlevgrnd)
    ss.hksat_min_col         = fill(FT(SPVAL), nc, nlevgrnd)
    ss.hk_l_col              = fill(FT(NaN), nc, nlevgrnd)
    ss.smp_l_col             = fill(FT(NaN), nc, nlevgrnd)
    ss.smpmin_col            = fill(FT(NaN), nc)
    ss.bsw_col               = fill(FT(NaN), nc, nlevgrnd)
    ss.watsat_col            = fill(FT(NaN), nc, nlevmaxurbgrnd)
    ss.watdry_col            = fill(FT(SPVAL), nc, nlevgrnd)
    ss.watopt_col            = fill(FT(SPVAL), nc, nlevgrnd)
    ss.watfc_col             = fill(FT(NaN), nc, nlevgrnd)
    ss.sucsat_col            = fill(FT(SPVAL), nc, nlevgrnd)
    ss.dsl_col               = fill(FT(SPVAL), nc)
    ss.soilresis_col         = fill(FT(SPVAL), nc)
    ss.soilbeta_col          = fill(FT(NaN), nc)
    ss.soilalpha_col         = fill(FT(NaN), nc)
    ss.soilalpha_u_col       = fill(FT(NaN), nc)
    ss.soilpsi_col           = fill(FT(NaN), nc, nlevgrnd)
    ss.wtfact_col            = fill(FT(NaN), nc)
    ss.porosity_col          = fill(FT(SPVAL), nc, NLAYER)
    ss.eff_porosity_col      = fill(FT(SPVAL), nc, nlevgrnd)
    ss.gwc_thr_col           = fill(FT(NaN), nc)

    # --- Van Genuchten parameters ---
    ss.msw_col               = fill(FT(NaN), nc, nlevgrnd)
    ss.nsw_col               = fill(FT(NaN), nc, nlevgrnd)
    ss.alphasw_col           = fill(FT(NaN), nc, nlevgrnd)
    ss.watres_col            = fill(FT(NaN), nc, nlevgrnd)

    # --- Thermal conductivity / heat capacity ---
    ss.thk_col               = fill(FT(NaN), nc, nlev_soisno)
    ss.tkmg_col              = fill(FT(NaN), nc, nlevgrnd)
    ss.tkdry_col             = fill(FT(NaN), nc, nlevgrnd)
    ss.tksatu_col            = fill(FT(NaN), nc, nlevgrnd)
    ss.csol_col              = fill(FT(NaN), nc, nlevgrnd)

    # --- Roots ---
    ss.rootr_patch           = fill(FT(NaN), np, nlevgrnd)
    ss.root_depth_patch      = fill(FT(NaN), np)
    ss.rootr_col             = fill(FT(NaN), nc, nlevgrnd)
    ss.rootr_road_perv_col   = fill(FT(NaN), nc, nlevgrnd)
    ss.rootfr_patch          = fill(FT(NaN), np, nlevgrnd)
    ss.crootfr_patch         = fill(FT(NaN), np, nlevgrnd)
    ss.rootfr_col            = fill(FT(NaN), nc, nlevgrnd)
    ss.rootfr_road_perv_col  = fill(FT(NaN), nc, nlevgrnd)
    ss.k_soil_root_patch     = fill(FT(NaN), np, nlevsoi)
    ss.root_conductance_patch = fill(FT(NaN), np, nlevsoi)
    ss.soil_conductance_patch = fill(FT(NaN), np, nlevsoi)

    return nothing
end

"""
    soilstate_clean!(ss::SoilStateData)

Deallocate (reset to empty) all fields of a `SoilStateData` instance.
"""
function soilstate_clean!(ss::SoilStateData{FT}) where {FT}
    # Patch-level 1D
    ss.sandfrac_patch        = FT[]
    ss.clayfrac_patch        = FT[]
    ss.root_depth_patch      = FT[]

    # Column-level 1D
    ss.mss_frc_cly_vld_col   = FT[]
    ss.smpmin_col            = FT[]
    ss.dsl_col               = FT[]
    ss.soilresis_col         = FT[]
    ss.soilbeta_col          = FT[]
    ss.soilalpha_col         = FT[]
    ss.soilalpha_u_col       = FT[]
    ss.wtfact_col            = FT[]
    ss.gwc_thr_col           = FT[]

    # Column-level 2D
    ss.cellorg_col           = Matrix{FT}(undef, 0, 0)
    ss.cellsand_col          = Matrix{FT}(undef, 0, 0)
    ss.cellclay_col          = Matrix{FT}(undef, 0, 0)
    ss.bd_col                = Matrix{FT}(undef, 0, 0)
    ss.hksat_col             = Matrix{FT}(undef, 0, 0)
    ss.hksat_min_col         = Matrix{FT}(undef, 0, 0)
    ss.hk_l_col              = Matrix{FT}(undef, 0, 0)
    ss.smp_l_col             = Matrix{FT}(undef, 0, 0)
    ss.bsw_col               = Matrix{FT}(undef, 0, 0)
    ss.watsat_col            = Matrix{FT}(undef, 0, 0)
    ss.watdry_col            = Matrix{FT}(undef, 0, 0)
    ss.watopt_col            = Matrix{FT}(undef, 0, 0)
    ss.watfc_col             = Matrix{FT}(undef, 0, 0)
    ss.sucsat_col            = Matrix{FT}(undef, 0, 0)
    ss.soilpsi_col           = Matrix{FT}(undef, 0, 0)
    ss.porosity_col          = Matrix{FT}(undef, 0, 0)
    ss.eff_porosity_col      = Matrix{FT}(undef, 0, 0)
    ss.msw_col               = Matrix{FT}(undef, 0, 0)
    ss.nsw_col               = Matrix{FT}(undef, 0, 0)
    ss.alphasw_col           = Matrix{FT}(undef, 0, 0)
    ss.watres_col            = Matrix{FT}(undef, 0, 0)
    ss.thk_col               = Matrix{FT}(undef, 0, 0)
    ss.tkmg_col              = Matrix{FT}(undef, 0, 0)
    ss.tkdry_col             = Matrix{FT}(undef, 0, 0)
    ss.tksatu_col            = Matrix{FT}(undef, 0, 0)
    ss.csol_col              = Matrix{FT}(undef, 0, 0)
    ss.rootr_col             = Matrix{FT}(undef, 0, 0)
    ss.rootfr_col            = Matrix{FT}(undef, 0, 0)
    ss.rootr_road_perv_col   = Matrix{FT}(undef, 0, 0)
    ss.rootfr_road_perv_col  = Matrix{FT}(undef, 0, 0)

    # Patch-level 2D
    ss.rootr_patch           = Matrix{FT}(undef, 0, 0)
    ss.rootfr_patch          = Matrix{FT}(undef, 0, 0)
    ss.crootfr_patch         = Matrix{FT}(undef, 0, 0)
    ss.k_soil_root_patch     = Matrix{FT}(undef, 0, 0)
    ss.root_conductance_patch = Matrix{FT}(undef, 0, 0)
    ss.soil_conductance_patch = Matrix{FT}(undef, 0, 0)

    return nothing
end

"""
    soilstate_init_cold!(ss::SoilStateData, bounds_col::UnitRange{Int})

Initialize cold-start conditions for soil state variables.
Sets soil matric potential to -1000 mm and hydraulic conductivity to 0.

Ported from `soilstate_type%InitCold` in `SoilStateType.F90`.
"""
function soilstate_init_cold!(ss::SoilStateData{FT}, bounds_col::UnitRange{Int}) where {FT}
    nlevgrnd = varpar.nlevgrnd

    for c in bounds_col
        for j in 1:nlevgrnd
            ss.smp_l_col[c, j] = -1000.0
            ss.hk_l_col[c, j]  = 0.0
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
    soilstate_init_history!(ss, bounds_col, bounds_patch)

Register soil state fields for history file output.

Ported from `soilstate_type%InitHistory` in `SoilStateType.F90`.
Requires history infrastructure (histFileMod) — stub until that module is ported.
"""
function soilstate_init_history!(ss::SoilStateData,
                                  bounds_col::UnitRange{Int},
                                  bounds_patch::UnitRange{Int})
    # Stub: history field registration will be added when histFileMod is ported.
    # All fields that would be registered:
    #   SMP, KROOT, KSOIL, bsw, ROOTR, ROOTR_COLUMN, SOILPSI,
    #   SNO_TK, SNO_TK_ICE, HK, SoilAlpha, SoilAlpha_U,
    #   watsat, EFF_POROSITY, watfc, SOILRESIS, DSL
    return nothing
end

"""
    soilstate_restart!(ss, bounds_col, bounds_patch; flag="read")

Read/write soil state from/to restart file.

Ported from `soilstate_type%Restart` in `SoilStateType.F90`.
Requires NetCDF/restart infrastructure — stub until that module is ported.
"""
function soilstate_restart!(ss::SoilStateData,
                             bounds_col::UnitRange{Int},
                             bounds_patch::UnitRange{Int};
                             flag::String = "read")
    # Stub: restart variable I/O will be added when restUtilMod/ncdio_pio is ported.
    # Variables that would be read/written:
    #   DSL, SOILRESIS, SMP, HK
    # On read, rootfr_patch and crootfr_patch are initialized via init_vegrootfr
    return nothing
end
