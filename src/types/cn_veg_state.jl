# ==========================================================================
# Ported from: src/biogeochem/CNVegStateType.F90
# Vegetation CN state data type allocation and initialization
# ==========================================================================

"""
    CNVegStateData

Vegetation carbon-nitrogen state data structure. Holds prognostic crop model
variables, fire-related variables, phenology flags, allocation coefficients,
and annual accumulator variables at patch and column levels.

Ported from `cnveg_state_type` in `CNVegStateType.F90`.
"""
Base.@kwdef mutable struct CNVegStateData{FT<:Real}
    # --- Patch-level integer fields ---
    burndate_patch              ::Vector{Int}     = Int[]       # patch crop burn date
    peaklai_patch               ::Vector{Int}     = Int[]       # patch 1: max allowed lai; 0: not at max
    idop_patch                  ::Vector{Int}     = Int[]       # patch date of planting (day of year)
    iyop_patch                  ::Vector{Int}     = Int[]       # patch year of planting

    # --- Dynamic weight smoothing ---
    # dwt_dribbler_patch is an annual_flux_dribbler_type in Fortran;
    # stub: will be added when AnnualFluxDribbler is ported
    dwt_smoothed_patch          ::Vector{FT} = Float64[]   # patch change in patch weight (-1 to 1) on gridcell

    # --- Prognostic crop model (patch-level) ---
    hdidx_patch                 ::Vector{FT} = Float64[]   # patch cold hardening index
    cumvd_patch                 ::Vector{FT} = Float64[]   # patch cumulative vernalization dependence
    gddmaturity_patch           ::Vector{FT} = Float64[]   # patch growing degree days needed to harvest (ddays)
    gddmaturity_thisyr          ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # at-harvest GDD this year (patch x mxharvests)
    huileaf_patch               ::Vector{FT} = Float64[]   # patch heat unit index needed from planting to leaf emergence
    huigrain_patch              ::Vector{FT} = Float64[]   # patch heat unit index needed to reach vegetative maturity
    aleafi_patch                ::Vector{FT} = Float64[]   # patch saved leaf allocation coefficient from phase 2
    astemi_patch                ::Vector{FT} = Float64[]   # patch saved stem allocation coefficient from phase 2

    # --- Allocation coefficients (patch-level) ---
    aleaf_patch                 ::Vector{FT} = Float64[]   # patch leaf allocation coefficient
    astem_patch                 ::Vector{FT} = Float64[]   # patch stem allocation coefficient
    aroot_patch                 ::Vector{FT} = Float64[]   # patch root allocation coefficient
    arepr_patch                 ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # patch reproductive allocation coefficient(s) (patch x nrepr)

    # --- AgSys nitrogen allocation coefficients (patch-level) ---
    aleaf_n_patch               ::Vector{FT} = Float64[]   # patch leaf allocation coefficient for N
    astem_n_patch               ::Vector{FT} = Float64[]   # patch stem allocation coefficient for N
    aroot_n_patch               ::Vector{FT} = Float64[]   # patch root allocation coefficient for N
    arepr_n_patch               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # patch reproductive allocation coefficient(s) for N (patch x nrepr)

    # --- Crop height/LAI (patch-level) ---
    htmx_patch                  ::Vector{FT} = Float64[]   # patch max hgt attained by a crop during yr (m)

    # --- Column-level GDP/population fire factors ---
    lgdp_col                    ::Vector{FT} = Float64[]   # col gdp limitation factor for fire occurrence (0-1)
    lgdp1_col                   ::Vector{FT} = Float64[]   # col gdp limitation factor for fire spreading (0-1)
    lpop_col                    ::Vector{FT} = Float64[]   # col pop limitation factor for fire spreading (0-1)

    # --- Temperature averages ---
    tempavg_t2m_patch           ::Vector{FT} = Float64[]   # patch temporary average 2m air temperature (K)
    annavg_t2m_patch            ::Vector{FT} = Float64[]   # patch annual average 2m air temperature (K)
    annavg_t2m_col              ::Vector{FT} = Float64[]   # col annual average of 2m air temperature (K)
    annsum_counter_col          ::Vector{FT} = Float64[]   # col seconds since last annual accumulator turnover

    # --- Fire (column-level) ---
    nfire_col                   ::Vector{FT} = Float64[]   # col fire counts (count/km2/sec)
    fsr_col                     ::Vector{FT} = Float64[]   # col fire spread rate (m/s)
    fd_col                      ::Vector{FT} = Float64[]   # col fire duration (hr)
    lfc_col                     ::Vector{FT} = Float64[]   # col conversion area fraction of BET and BDT that haven't burned before (/timestep)
    lfc2_col                    ::Vector{FT} = Float64[]   # col conversion area fraction of BET and BDT that burned (/sec)
    dtrotr_col                  ::Vector{FT} = Float64[]   # col annual decreased fraction coverage of BET on gridcell (0-1)
    trotr1_col                  ::Vector{FT} = Float64[]   # col patch weight of BET on column (0-1)
    trotr2_col                  ::Vector{FT} = Float64[]   # col patch weight of BDT on column (0-1)
    cropf_col                   ::Vector{FT} = Float64[]   # col crop fraction in veg column (0-1)
    baf_crop_col                ::Vector{FT} = Float64[]   # col baf for cropland (/sec)
    baf_peatf_col               ::Vector{FT} = Float64[]   # col baf for peatland (/sec)
    fbac_col                    ::Vector{FT} = Float64[]   # col total burned area out of conversion (/sec)
    fbac1_col                   ::Vector{FT} = Float64[]   # col burned area out of conversion region due to land use fire (/sec)
    wtlf_col                    ::Vector{FT} = Float64[]   # col fractional coverage of non-crop patches (0-1)
    lfwt_col                    ::Vector{FT} = Float64[]   # col fractional coverage of non-crop and non-bare-soil patches (0-1)
    farea_burned_col            ::Vector{FT} = Float64[]   # col fractional area burned (/sec)

    # --- Phenology flags (patch-level) ---
    dormant_flag_patch          ::Vector{FT} = Float64[]   # patch dormancy flag
    days_active_patch           ::Vector{FT} = Float64[]   # patch number of days since last dormancy
    onset_flag_patch            ::Vector{FT} = Float64[]   # patch onset flag
    onset_counter_patch         ::Vector{FT} = Float64[]   # patch onset days counter
    onset_gddflag_patch         ::Vector{FT} = Float64[]   # patch onset flag for growing degree day sum
    onset_fdd_patch             ::Vector{FT} = Float64[]   # patch onset freezing degree days counter
    onset_gdd_patch             ::Vector{FT} = Float64[]   # patch onset growing degree days
    onset_swi_patch             ::Vector{FT} = Float64[]   # patch onset soil water index
    offset_flag_patch           ::Vector{FT} = Float64[]   # patch offset flag
    offset_counter_patch        ::Vector{FT} = Float64[]   # patch offset days counter
    offset_fdd_patch            ::Vector{FT} = Float64[]   # patch offset freezing degree days counter
    offset_swi_patch            ::Vector{FT} = Float64[]   # patch offset soil water index
    grain_flag_patch            ::Vector{FT} = Float64[]   # patch 1: grain fill stage; 0: not
    lgsf_patch                  ::Vector{FT} = Float64[]   # patch long growing season factor [0-1]
    bglfr_patch                 ::Vector{FT} = Float64[]   # patch background litterfall rate (1/s)
    bgtr_patch                  ::Vector{FT} = Float64[]   # patch background transfer growth rate (1/s)

    # --- Allometry and GPP accumulators (patch-level) ---
    c_allometry_patch           ::Vector{FT} = Float64[]   # patch C allocation index (DIM)
    n_allometry_patch           ::Vector{FT} = Float64[]   # patch N allocation index (DIM)
    tempsum_potential_gpp_patch ::Vector{FT} = Float64[]   # patch temporary annual sum of potential GPP
    annsum_potential_gpp_patch  ::Vector{FT} = Float64[]   # patch annual sum of potential GPP
    tempmax_retransn_patch      ::Vector{FT} = Float64[]   # patch temporary annual max of retranslocated N pool (gN/m2)
    annmax_retransn_patch       ::Vector{FT} = Float64[]   # patch annual max of retranslocated N pool (gN/m2)
    downreg_patch               ::Vector{FT} = Float64[]   # patch fractional reduction in GPP due to N limitation (DIM)
    leafcn_offset_patch         ::Vector{FT} = Float64[]   # patch leaf C:N used by FUN
    plantCN_patch               ::Vector{FT} = Float64[]   # patch plant C:N used by FUN
end

"""
    cnveg_state_init!(vs::CNVegStateData, np::Int, nc::Int;
                      use_crop_agsys::Bool=false)

Allocate and initialize all fields of a `CNVegStateData` instance for
`np` patches and `nc` columns. Matches the Fortran `InitAllocate`.

Arguments:
- `vs`: CNVegStateData instance
- `np`: number of patches
- `nc`: number of columns
- `use_crop_agsys`: whether to allocate AgSys nitrogen allocation fields
"""
function cnveg_state_init!(vs::CNVegStateData, np::Int, nc::Int;
                           use_crop_agsys::Bool=false)
    mxharvests = MXHARVESTS
    nrepr = NREPR

    # --- Patch-level integer fields ---
    vs.burndate_patch      = fill(ISPVAL, np)
    vs.peaklai_patch       = fill(0, np)
    vs.idop_patch          = fill(typemax(Int), np)
    vs.iyop_patch          = fill(ISPVAL, np)

    # --- Dynamic weight smoothing ---
    vs.dwt_smoothed_patch  = fill(NaN, np)

    # --- Prognostic crop model ---
    vs.hdidx_patch         = fill(NaN, np)
    vs.cumvd_patch         = fill(NaN, np)
    vs.gddmaturity_patch   = fill(SPVAL, np)
    vs.gddmaturity_thisyr  = fill(SPVAL, np, mxharvests)
    vs.huileaf_patch       = fill(NaN, np)
    vs.huigrain_patch      = fill(0.0, np)
    vs.aleafi_patch        = fill(NaN, np)
    vs.astemi_patch        = fill(NaN, np)

    # --- Allocation coefficients ---
    vs.aleaf_patch         = fill(NaN, np)
    vs.astem_patch         = fill(NaN, np)
    vs.aroot_patch         = fill(NaN, np)
    vs.arepr_patch         = fill(NaN, np, nrepr)

    # --- AgSys nitrogen allocation coefficients ---
    if use_crop_agsys
        vs.aleaf_n_patch   = fill(NaN, np)
        vs.astem_n_patch   = fill(NaN, np)
        vs.aroot_n_patch   = fill(NaN, np)
        vs.arepr_n_patch   = fill(NaN, np, nrepr)
    else
        vs.aleaf_n_patch   = Float64[]
        vs.astem_n_patch   = Float64[]
        vs.aroot_n_patch   = Float64[]
        vs.arepr_n_patch   = Matrix{Float64}(undef, 0, 0)
    end

    # --- Crop height ---
    vs.htmx_patch          = fill(0.0, np)

    # --- Column-level GDP/population factors ---
    vs.lgdp_col            = Vector{Float64}(undef, nc)
    vs.lgdp1_col           = Vector{Float64}(undef, nc)
    vs.lpop_col            = Vector{Float64}(undef, nc)

    # --- Temperature averages ---
    vs.tempavg_t2m_patch   = fill(NaN, np)
    vs.annavg_t2m_patch    = fill(NaN, np)
    vs.annavg_t2m_col      = fill(NaN, nc)
    vs.annsum_counter_col  = fill(NaN, nc)

    # --- Fire (column-level) ---
    vs.nfire_col           = fill(SPVAL, nc)
    vs.fsr_col             = fill(NaN, nc)
    vs.fd_col              = fill(NaN, nc)
    vs.lfc_col             = fill(SPVAL, nc)
    vs.lfc2_col            = fill(0.0, nc)
    vs.dtrotr_col          = fill(0.0, nc)
    vs.trotr1_col          = fill(0.0, nc)
    vs.trotr2_col          = fill(0.0, nc)
    vs.cropf_col           = fill(NaN, nc)
    vs.baf_crop_col        = fill(NaN, nc)
    vs.baf_peatf_col       = fill(NaN, nc)
    vs.fbac_col            = fill(NaN, nc)
    vs.fbac1_col           = fill(NaN, nc)
    vs.wtlf_col            = fill(NaN, nc)
    vs.lfwt_col            = fill(NaN, nc)
    vs.farea_burned_col    = fill(NaN, nc)

    # --- Phenology flags ---
    vs.dormant_flag_patch          = fill(NaN, np)
    vs.days_active_patch           = fill(NaN, np)
    vs.onset_flag_patch            = fill(NaN, np)
    vs.onset_counter_patch         = fill(NaN, np)
    vs.onset_gddflag_patch         = fill(NaN, np)
    vs.onset_fdd_patch             = fill(NaN, np)
    vs.onset_gdd_patch             = fill(NaN, np)
    vs.onset_swi_patch             = fill(NaN, np)
    vs.offset_flag_patch           = fill(NaN, np)
    vs.offset_counter_patch        = fill(NaN, np)
    vs.offset_fdd_patch            = fill(NaN, np)
    vs.offset_swi_patch            = fill(NaN, np)
    vs.grain_flag_patch            = fill(NaN, np)
    vs.lgsf_patch                  = fill(NaN, np)
    vs.bglfr_patch                 = fill(NaN, np)
    vs.bgtr_patch                  = fill(NaN, np)

    # --- Allometry and GPP accumulators ---
    vs.c_allometry_patch           = fill(NaN, np)
    vs.n_allometry_patch           = fill(NaN, np)
    vs.tempsum_potential_gpp_patch = fill(NaN, np)
    vs.annsum_potential_gpp_patch  = fill(NaN, np)
    vs.tempmax_retransn_patch      = fill(NaN, np)
    vs.annmax_retransn_patch       = fill(NaN, np)
    vs.downreg_patch               = fill(NaN, np)
    vs.leafcn_offset_patch         = fill(NaN, np)
    vs.plantCN_patch               = fill(NaN, np)

    return nothing
end

"""
    cnveg_state_clean!(vs::CNVegStateData)

Deallocate (reset to empty) all fields of a `CNVegStateData` instance.
"""
function cnveg_state_clean!(vs::CNVegStateData{FT}) where {FT}
    # Integer vectors
    vs.burndate_patch      = Int[]
    vs.peaklai_patch       = Int[]
    vs.idop_patch          = Int[]
    vs.iyop_patch          = Int[]

    # Float64 vectors - patch level
    vs.dwt_smoothed_patch  = FT[]
    vs.hdidx_patch         = FT[]
    vs.cumvd_patch         = FT[]
    vs.gddmaturity_patch   = FT[]
    vs.huileaf_patch       = FT[]
    vs.huigrain_patch      = FT[]
    vs.aleafi_patch        = FT[]
    vs.astemi_patch        = FT[]
    vs.aleaf_patch         = FT[]
    vs.astem_patch         = FT[]
    vs.aroot_patch         = FT[]
    vs.aleaf_n_patch       = FT[]
    vs.astem_n_patch       = FT[]
    vs.aroot_n_patch       = FT[]
    vs.htmx_patch          = FT[]
    vs.tempavg_t2m_patch   = FT[]
    vs.annavg_t2m_patch    = FT[]
    vs.dormant_flag_patch  = FT[]
    vs.days_active_patch   = FT[]
    vs.onset_flag_patch    = FT[]
    vs.onset_counter_patch = FT[]
    vs.onset_gddflag_patch = FT[]
    vs.onset_fdd_patch     = FT[]
    vs.onset_gdd_patch     = FT[]
    vs.onset_swi_patch     = FT[]
    vs.offset_flag_patch   = FT[]
    vs.offset_counter_patch = FT[]
    vs.offset_fdd_patch    = FT[]
    vs.offset_swi_patch    = FT[]
    vs.grain_flag_patch    = FT[]
    vs.lgsf_patch          = FT[]
    vs.bglfr_patch         = FT[]
    vs.bgtr_patch          = FT[]
    vs.c_allometry_patch   = FT[]
    vs.n_allometry_patch   = FT[]
    vs.tempsum_potential_gpp_patch = FT[]
    vs.annsum_potential_gpp_patch  = FT[]
    vs.tempmax_retransn_patch      = FT[]
    vs.annmax_retransn_patch       = FT[]
    vs.downreg_patch       = FT[]
    vs.leafcn_offset_patch = FT[]
    vs.plantCN_patch       = FT[]

    # Float64 vectors - column level
    vs.lgdp_col            = FT[]
    vs.lgdp1_col           = FT[]
    vs.lpop_col            = FT[]
    vs.annavg_t2m_col      = FT[]
    vs.annsum_counter_col  = FT[]
    vs.nfire_col           = FT[]
    vs.fsr_col             = FT[]
    vs.fd_col              = FT[]
    vs.lfc_col             = FT[]
    vs.lfc2_col            = FT[]
    vs.dtrotr_col          = FT[]
    vs.trotr1_col          = FT[]
    vs.trotr2_col          = FT[]
    vs.cropf_col           = FT[]
    vs.baf_crop_col        = FT[]
    vs.baf_peatf_col       = FT[]
    vs.fbac_col            = FT[]
    vs.fbac1_col           = FT[]
    vs.wtlf_col            = FT[]
    vs.lfwt_col            = FT[]
    vs.farea_burned_col    = FT[]

    # Matrices
    vs.gddmaturity_thisyr  = Matrix{FT}(undef, 0, 0)
    vs.arepr_patch         = Matrix{FT}(undef, 0, 0)
    vs.arepr_n_patch       = Matrix{FT}(undef, 0, 0)

    return nothing
end

"""
    cnveg_state_init_cold!(vs::CNVegStateData, bounds_col::UnitRange{Int},
                           bounds_patch::UnitRange{Int};
                           col_landunit::Union{Vector{Int},Nothing}=nothing,
                           patch_landunit::Union{Vector{Int},Nothing}=nothing,
                           lun_ifspecial::Union{Vector{Bool},Nothing}=nothing,
                           lun_itype::Union{Vector{Int},Nothing}=nothing)

Initialize cold-start conditions for CN vegetation state variables.
Sets fire, phenology, and accumulator variables to their cold-start values.

Ported from `cnveg_state_type%InitCold` in `CNVegStateType.F90`.
"""
function cnveg_state_init_cold!(vs::CNVegStateData,
                                bounds_col::UnitRange{Int},
                                bounds_patch::UnitRange{Int};
                                col_landunit::Union{Vector{Int},Nothing}=nothing,
                                patch_landunit::Union{Vector{Int},Nothing}=nothing,
                                lun_ifspecial::Union{Vector{Bool},Nothing}=nothing,
                                lun_itype::Union{Vector{Int},Nothing}=nothing)

    have_landunit_info = (col_landunit !== nothing && patch_landunit !== nothing &&
                          lun_ifspecial !== nothing && lun_itype !== nothing)

    # --- Column-level cold start ---
    for c in bounds_col
        if have_landunit_info
            l = col_landunit[c]

            # Special landunits get spval
            if lun_ifspecial[l]
                vs.annsum_counter_col[c] = SPVAL
                vs.annavg_t2m_col[c]     = SPVAL
                vs.nfire_col[c]          = SPVAL
                vs.baf_crop_col[c]       = SPVAL
                vs.baf_peatf_col[c]      = SPVAL
                vs.fbac_col[c]           = SPVAL
                vs.fbac1_col[c]          = SPVAL
                vs.farea_burned_col[c]   = SPVAL
            end

            # Soil/crop landunits
            if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
                vs.annsum_counter_col[c] = 0.0
                vs.annavg_t2m_col[c]     = 280.0
                vs.baf_crop_col[c]       = 0.0
                vs.baf_peatf_col[c]      = 0.0
                vs.fbac_col[c]           = 0.0
                vs.fbac1_col[c]          = 0.0
                vs.farea_burned_col[c]   = 0.0
                vs.nfire_col[c]          = 0.0
            end
        else
            # Default: assume soil/crop
            vs.annsum_counter_col[c] = 0.0
            vs.annavg_t2m_col[c]     = 280.0
            vs.baf_crop_col[c]       = 0.0
            vs.baf_peatf_col[c]      = 0.0
            vs.fbac_col[c]           = 0.0
            vs.fbac1_col[c]          = 0.0
            vs.farea_burned_col[c]   = 0.0
            vs.nfire_col[c]          = 0.0
        end
    end

    # --- Patch-level cold start ---
    for p in bounds_patch
        if have_landunit_info
            l = patch_landunit[p]

            # Special landunits get spval
            if lun_ifspecial[l]
                vs.annavg_t2m_patch[p]          = SPVAL
                vs.tempavg_t2m_patch[p]         = SPVAL
                vs.dormant_flag_patch[p]        = SPVAL
                vs.days_active_patch[p]         = SPVAL
                vs.onset_flag_patch[p]          = SPVAL
                vs.onset_counter_patch[p]       = SPVAL
                vs.onset_gddflag_patch[p]       = SPVAL
                vs.onset_fdd_patch[p]           = SPVAL
                vs.onset_gdd_patch[p]           = SPVAL
                vs.onset_swi_patch[p]           = SPVAL
                vs.offset_flag_patch[p]         = SPVAL
                vs.offset_counter_patch[p]      = SPVAL
                vs.offset_fdd_patch[p]          = SPVAL
                vs.offset_swi_patch[p]          = SPVAL
                vs.grain_flag_patch[p]          = SPVAL
                vs.lgsf_patch[p]                = SPVAL
                vs.bglfr_patch[p]               = SPVAL
                vs.bgtr_patch[p]                = SPVAL
                vs.c_allometry_patch[p]         = SPVAL
                vs.n_allometry_patch[p]         = SPVAL
                vs.tempsum_potential_gpp_patch[p] = SPVAL
                vs.annsum_potential_gpp_patch[p]  = SPVAL
                vs.tempmax_retransn_patch[p]    = SPVAL
                vs.annmax_retransn_patch[p]     = SPVAL
                vs.downreg_patch[p]             = SPVAL
                vs.leafcn_offset_patch[p]       = SPVAL
                vs.plantCN_patch[p]             = SPVAL
            end

            # Soil/crop landunits
            if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
                # phenology variables
                vs.dormant_flag_patch[p]   = 1.0
                vs.days_active_patch[p]    = 0.0
                vs.onset_flag_patch[p]     = 0.0
                vs.onset_counter_patch[p]  = 0.0
                vs.onset_gddflag_patch[p]  = 0.0
                vs.onset_fdd_patch[p]      = 0.0
                vs.onset_gdd_patch[p]      = 0.0
                vs.onset_swi_patch[p]      = 0.0
                vs.offset_flag_patch[p]    = 0.0
                vs.offset_counter_patch[p] = 0.0
                vs.offset_fdd_patch[p]     = 0.0
                vs.offset_swi_patch[p]     = 0.0
                vs.lgsf_patch[p]           = 0.0
                vs.bglfr_patch[p]          = 0.0
                vs.bgtr_patch[p]           = 0.0
                vs.annavg_t2m_patch[p]     = 280.0
                vs.tempavg_t2m_patch[p]    = 0.0
                vs.grain_flag_patch[p]     = 0.0

                # non-phenology variables
                vs.c_allometry_patch[p]           = 0.0
                vs.n_allometry_patch[p]           = 0.0
                vs.tempsum_potential_gpp_patch[p] = 0.0
                vs.annsum_potential_gpp_patch[p]  = 0.0
                vs.tempmax_retransn_patch[p]      = 0.0
                vs.annmax_retransn_patch[p]       = 0.0
                vs.downreg_patch[p]               = 0.0
                vs.leafcn_offset_patch[p]         = SPVAL
                vs.plantCN_patch[p]               = SPVAL
            end
        else
            # Default: assume soil/crop
            vs.dormant_flag_patch[p]   = 1.0
            vs.days_active_patch[p]    = 0.0
            vs.onset_flag_patch[p]     = 0.0
            vs.onset_counter_patch[p]  = 0.0
            vs.onset_gddflag_patch[p]  = 0.0
            vs.onset_fdd_patch[p]      = 0.0
            vs.onset_gdd_patch[p]      = 0.0
            vs.onset_swi_patch[p]      = 0.0
            vs.offset_flag_patch[p]    = 0.0
            vs.offset_counter_patch[p] = 0.0
            vs.offset_fdd_patch[p]     = 0.0
            vs.offset_swi_patch[p]     = 0.0
            vs.lgsf_patch[p]           = 0.0
            vs.bglfr_patch[p]          = 0.0
            vs.bgtr_patch[p]           = 0.0
            vs.annavg_t2m_patch[p]     = 280.0
            vs.tempavg_t2m_patch[p]    = 0.0
            vs.grain_flag_patch[p]     = 0.0
            vs.c_allometry_patch[p]           = 0.0
            vs.n_allometry_patch[p]           = 0.0
            vs.tempsum_potential_gpp_patch[p] = 0.0
            vs.annsum_potential_gpp_patch[p]  = 0.0
            vs.tempmax_retransn_patch[p]      = 0.0
            vs.annmax_retransn_patch[p]       = 0.0
            vs.downreg_patch[p]               = 0.0
            vs.leafcn_offset_patch[p]         = SPVAL
            vs.plantCN_patch[p]               = SPVAL
        end
    end

    # Fire variable: lfc2 initialized to zero for all columns
    for c in bounds_col
        vs.lfc2_col[c] = 0.0
    end

    return nothing
end

# ==========================================================================
# The following subroutines depend on infrastructure modules that are not yet
# ported (history, restart/IO). They are provided as stubs that document the
# Fortran interface and can be filled in when those modules become available.
# ==========================================================================

"""
    cnveg_state_init_history!(vs, bounds_patch, bounds_col)

Register CN vegetation state fields for history file output.

Ported from `cnveg_state_type%InitHistory` in `CNVegStateType.F90`.
Requires history infrastructure (histFileMod) - stub until that module is ported.
"""
function cnveg_state_init_history!(vs::CNVegStateData,
                                   bounds_patch::UnitRange{Int},
                                   bounds_col::UnitRange{Int})
    # Stub: history field registration will be added when histFileMod is ported.
    # Fields that would be registered:
    #   GDDHARV, GDDHARV_PERHARV, LFC2, ANNSUM_COUNTER, CANNAVG_T2M,
    #   NFIRE, FAREA_BURNED, BAF_CROP, BAF_PEATF, ANNAVG_T2M, TEMPAVG_T2M,
    #   DORMANT_FLAG, DAYS_ACTIVE, ONSET_FLAG, ONSET_COUNTER, ONSET_GDDFLAG,
    #   ONSET_FDD, ONSET_GDD, ONSET_SWI, OFFSET_FLAG, OFFSET_COUNTER,
    #   OFFSET_FDD, OFFSET_SWI, LGSF, BGLFR, BGTR, C_ALLOMETRY, N_ALLOMETRY,
    #   TEMPSUM_POTENTIAL_GPP, ANNSUM_POTENTIAL_GPP, TEMPMAX_RETRANSN,
    #   ANNMAX_RETRANSN, DOWNREG, LEAFCN_OFFSET, PLANTCN
    return nothing
end

"""
    cnveg_state_restart!(vs, bounds_patch, bounds_col; flag="read")

Read/write CN vegetation state from/to restart file.

Ported from `cnveg_state_type%Restart` in `CNVegStateType.F90`.
Requires NetCDF/restart infrastructure - stub until that module is ported.
"""
function cnveg_state_restart!(vs::CNVegStateData,
                               bounds_patch::UnitRange{Int},
                               bounds_col::UnitRange{Int};
                               flag::String="read")
    # Stub: restart variable I/O will be added when restUtilMod/ncdio_pio is ported.
    # Variables that would be read/written:
    #   dormant_flag, days_active, onset_flag, onset_counter, onset_gddflag,
    #   onset_fdd, onset_gdd, onset_swi, offset_flag, offset_counter,
    #   offset_fdd, offset_swi, lgsf, bglfr, bgtr, annavg_t2m, tempavg_t2m,
    #   c_allometry, n_allometry, tempsum_potential_gpp, annsum_potential_gpp,
    #   tempmax_retransn, annmax_retransn, downreg, leafcn_offset, plantCN,
    #   annsum_counter, burndate, lfc, cannavg_t2m,
    #   htmx, peaklai, idop, aleaf, aleafi, astem, astemi, hdidx, cumvd,
    #   gddmaturity, huileaf, huigrain, grain_flag, gddmaturity_thisyr
    return nothing
end
