# ==========================================================================
# Ported from: src/biogeophys/SoilHydrologyType.F90
# Soil hydrology data type allocation and initialization
# ==========================================================================

"""
    SoilHydrologyData

Soil hydrology data structure. Holds all soil hydrology state variables
at the column level, including water table depths, ice fractions, aquifer
recharge, and VIC (Variable Infiltration Capacity) model parameters.

Ported from `soilhydrology_type` in `SoilHydrologyType.F90`.
"""
Base.@kwdef mutable struct SoilHydrologyData
    # --- Control flag ---
    h2osfcflag::Int = 1                # true => surface water is active (namelist)

    # --- NON-VIC state ---
    num_substeps_col     ::Vector{Float64} = Float64[]   # col adaptive timestep counter
    frost_table_col      ::Vector{Float64} = Float64[]   # col frost table depth (m)
    zwt_col              ::Vector{Float64} = Float64[]   # col water table depth (m)
    zwts_col             ::Vector{Float64} = Float64[]   # col water table depth, shallower of two water depths (m)
    zwt_perched_col      ::Vector{Float64} = Float64[]   # col perched water table depth (m)
    qcharge_col          ::Vector{Float64} = Float64[]   # col aquifer recharge rate (mm/s)
    icefrac_col          ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col fraction of ice (ncols, nlevgrnd)
    h2osfc_thresh_col    ::Vector{Float64} = Float64[]   # col level at which h2osfc "percolates" (time constant)
    xs_urban_col         ::Vector{Float64} = Float64[]   # col excess soil water above urban ponding limit

    # --- VIC parameters ---
    hkdepth_col          ::Vector{Float64} = Float64[]   # col VIC decay factor (m) (time constant)
    b_infil_col          ::Vector{Float64} = Float64[]   # col VIC b infiltration parameter (time constant)
    ds_col               ::Vector{Float64} = Float64[]   # col VIC fraction of Dsmax where non-linear baseflow begins (time constant)
    dsmax_col            ::Vector{Float64} = Float64[]   # col VIC max velocity of baseflow (mm/day) (time constant)
    Wsvic_col            ::Vector{Float64} = Float64[]   # col VIC fraction of max soil moisture where non-linear baseflow occurs (time constant)
    porosity_col         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col VIC porosity (1-bulk_density/soil_density) (ncols, nlayer)
    vic_clm_fract_col    ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)  # col VIC fraction of VIC layers in CLM layers (ncols, nlayer, nlevsoi)
    depth_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col VIC layer depth of upper layer (ncols, nlayert)
    c_param_col          ::Vector{Float64} = Float64[]   # col VIC baseflow exponent (Qb)
    expt_col             ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col VIC pore-size distribution related parameter (Q12) (ncols, nlayer)
    ksat_col             ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col VIC saturated hydrologic conductivity (ncols, nlayer)
    phi_s_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col VIC soil moisture diffusion parameter (ncols, nlayer)
    moist_col            ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col VIC soil moisture (kg/m2) for VIC soil layers (ncols, nlayert)
    moist_vol_col        ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col VIC volumetric soil moisture for VIC soil layers (ncols, nlayert)
    max_moist_col        ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col VIC max layer moist + ice (mm) (ncols, nlayer)
    top_moist_col        ::Vector{Float64} = Float64[]   # col VIC soil moisture in top layers
    top_max_moist_col    ::Vector{Float64} = Float64[]   # col VIC maximum soil moisture in top layers
    top_ice_col          ::Vector{Float64} = Float64[]   # col VIC ice len in top layers
    top_moist_limited_col::Vector{Float64} = Float64[]   # col VIC soil moisture in top layers, limited to no greater than top_max_moist_col
    ice_col              ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # col VIC soil ice (kg/m2) for VIC soil layers (ncols, nlayert)
end

"""
    soilhydrology_init!(sh::SoilHydrologyData, nc::Int)

Allocate and initialize all fields of a `SoilHydrologyData` instance for
`nc` columns. Real fields are initialized to `NaN` matching the Fortran
source (`InitAllocate` sets all to `nan`).

Ported from `soilhydrology_type%InitAllocate` in `SoilHydrologyType.F90`.
"""
function soilhydrology_init!(sh::SoilHydrologyData, nc::Int)
    nlevgrnd = varpar.nlevgrnd
    nlevsoi  = varpar.nlevsoi
    nlayer   = NLAYER
    nlayert  = varpar.nlayert

    # --- NON-VIC 1D ---
    sh.num_substeps_col      = fill(NaN, nc)
    sh.frost_table_col       = fill(NaN, nc)
    sh.zwt_col               = fill(NaN, nc)
    sh.zwt_perched_col       = fill(NaN, nc)
    sh.zwts_col              = fill(NaN, nc)
    sh.qcharge_col           = fill(NaN, nc)
    sh.h2osfc_thresh_col     = fill(NaN, nc)
    sh.xs_urban_col          = fill(NaN, nc)

    # --- NON-VIC 2D ---
    sh.icefrac_col           = fill(NaN, nc, nlevgrnd)

    # --- VIC 1D ---
    sh.hkdepth_col           = fill(NaN, nc)
    sh.b_infil_col           = fill(NaN, nc)
    sh.ds_col                = fill(NaN, nc)
    sh.dsmax_col             = fill(NaN, nc)
    sh.Wsvic_col             = fill(NaN, nc)
    sh.c_param_col           = fill(NaN, nc)
    sh.top_moist_col         = fill(NaN, nc)
    sh.top_max_moist_col     = fill(NaN, nc)
    sh.top_ice_col           = fill(NaN, nc)
    sh.top_moist_limited_col = fill(NaN, nc)

    # --- VIC 2D ---
    sh.porosity_col          = fill(NaN, nc, nlayer)
    sh.depth_col             = fill(NaN, nc, nlayert)
    sh.expt_col              = fill(NaN, nc, nlayer)
    sh.ksat_col              = fill(NaN, nc, nlayer)
    sh.phi_s_col             = fill(NaN, nc, nlayer)
    sh.moist_col             = fill(NaN, nc, nlayert)
    sh.moist_vol_col         = fill(NaN, nc, nlayert)
    sh.max_moist_col         = fill(NaN, nc, nlayer)
    sh.ice_col               = fill(NaN, nc, nlayert)

    # --- VIC 3D ---
    sh.vic_clm_fract_col     = fill(NaN, nc, nlayer, nlevsoi)

    return nothing
end

"""
    soilhydrology_clean!(sh::SoilHydrologyData)

Deallocate (reset to empty) all fields of a `SoilHydrologyData` instance.
"""
function soilhydrology_clean!(sh::SoilHydrologyData)
    # 1D vectors
    sh.num_substeps_col      = Float64[]
    sh.frost_table_col       = Float64[]
    sh.zwt_col               = Float64[]
    sh.zwts_col              = Float64[]
    sh.zwt_perched_col       = Float64[]
    sh.qcharge_col           = Float64[]
    sh.h2osfc_thresh_col     = Float64[]
    sh.xs_urban_col          = Float64[]
    sh.hkdepth_col           = Float64[]
    sh.b_infil_col           = Float64[]
    sh.ds_col                = Float64[]
    sh.dsmax_col             = Float64[]
    sh.Wsvic_col             = Float64[]
    sh.c_param_col           = Float64[]
    sh.top_moist_col         = Float64[]
    sh.top_max_moist_col     = Float64[]
    sh.top_ice_col           = Float64[]
    sh.top_moist_limited_col = Float64[]

    # 2D matrices
    sh.icefrac_col           = Matrix{Float64}(undef, 0, 0)
    sh.porosity_col          = Matrix{Float64}(undef, 0, 0)
    sh.depth_col             = Matrix{Float64}(undef, 0, 0)
    sh.expt_col              = Matrix{Float64}(undef, 0, 0)
    sh.ksat_col              = Matrix{Float64}(undef, 0, 0)
    sh.phi_s_col             = Matrix{Float64}(undef, 0, 0)
    sh.moist_col             = Matrix{Float64}(undef, 0, 0)
    sh.moist_vol_col         = Matrix{Float64}(undef, 0, 0)
    sh.max_moist_col         = Matrix{Float64}(undef, 0, 0)
    sh.ice_col               = Matrix{Float64}(undef, 0, 0)

    # 3D array
    sh.vic_clm_fract_col     = Array{Float64}(undef, 0, 0, 0)

    return nothing
end

"""
    soilhydrology_init_cold!(sh::SoilHydrologyData, bounds_col::UnitRange{Int};
                             use_aquifer_layer::Bool=true,
                             zi_col::Union{Matrix{Float64},Nothing}=nothing,
                             wa_col::Union{Vector{Float64},Nothing}=nothing,
                             nlevsoi::Int=0,
                             landunit_col::Union{Vector{Int},Nothing}=nothing,
                             lakpoi::Union{BitVector,Nothing}=nothing,
                             urbpoi::Union{BitVector,Nothing}=nothing,
                             itype_col::Union{Vector{Int},Nothing}=nothing,
                             nbedrock_col::Union{Vector{Int},Nothing}=nothing)

Initialize cold-start conditions for soil hydrology variables.
Sets zwt, frost_table, and zwt_perched based on column types.

Ported from `soilhydrology_type%InitCold` in `SoilHydrologyType.F90`.
"""
function soilhydrology_init_cold!(sh::SoilHydrologyData, bounds_col::UnitRange{Int};
                                   use_aquifer_layer::Bool = true,
                                   zi_col::Union{Matrix{Float64},Nothing} = nothing,
                                   wa_col::Union{Vector{Float64},Nothing} = nothing,
                                   nlevsoi::Int = varpar.nlevsoi,
                                   landunit_col::Union{Vector{Int},Nothing} = nothing,
                                   lakpoi::Union{BitVector,Nothing} = nothing,
                                   urbpoi::Union{BitVector,Nothing} = nothing,
                                   itype_col::Union{Vector{Int},Nothing} = nothing,
                                   nbedrock_col::Union{Vector{Int},Nothing} = nothing)
    # Initialize num_substeps_col to SPVAL (needed for accum field averaging)
    for c in bounds_col
        sh.num_substeps_col[c] = SPVAL
    end

    # Initialize zwt_col to zero
    for c in bounds_col
        sh.zwt_col[c] = 0.0
    end

    # Full cold-start initialization requires column/landunit info.
    # When called without column metadata (basic init), just set defaults.
    if zi_col === nothing || landunit_col === nothing
        for c in bounds_col
            sh.zwt_perched_col[c] = 0.0
            sh.frost_table_col[c] = 0.0
        end
        return nothing
    end

    for c in bounds_col
        l = landunit_col[c]
        if !(lakpoi !== nothing && lakpoi[l])  # not lake
            if urbpoi !== nothing && urbpoi[l]
                # Urban column
                if itype_col !== nothing && itype_col[c] == ICOL_ROAD_PERV
                    if use_aquifer_layer
                        sh.zwt_col[c] = (25.0 + zi_col[c, nlevsoi]) - wa_col[c] / 0.2 / 1000.0
                    else
                        sh.zwt_col[c] = zi_col[c, nbedrock_col[c]]
                    end
                else
                    sh.zwt_col[c] = SPVAL
                end
                sh.zwt_perched_col[c] = SPVAL
                sh.frost_table_col[c] = SPVAL
            else
                # Non-urban, non-lake
                if use_aquifer_layer
                    sh.zwt_col[c] = (25.0 + zi_col[c, nlevsoi]) - wa_col[c] / 0.2 / 1000.0
                else
                    sh.zwt_col[c] = zi_col[c, nbedrock_col[c]]
                end
                sh.zwt_perched_col[c] = zi_col[c, nlevsoi]
                sh.frost_table_col[c] = zi_col[c, nlevsoi]
            end
        end
    end

    return nothing
end

# ==========================================================================
# The following subroutines depend on infrastructure modules that are not yet
# ported (history, restart/IO, namelist reading). They are provided as stubs
# that document the Fortran interface and can be filled in when those modules
# become available.
# ==========================================================================

"""
    soilhydrology_init_history!(sh, bounds_col; use_aquifer_layer=true)

Register soil hydrology fields for history file output.

Ported from `soilhydrology_type%InitHistory` in `SoilHydrologyType.F90`.
Requires history infrastructure (histFileMod) — stub until that module is ported.
"""
function soilhydrology_init_history!(sh::SoilHydrologyData,
                                      bounds_col::UnitRange{Int};
                                      use_aquifer_layer::Bool = true)
    # Stub: history field registration will be added when histFileMod is ported.
    # Fields that would be registered:
    #   QCHARGE (if use_aquifer_layer), NSUBSTEPS, FROST_TABLE, ZWT, ZWT_PERCH
    return nothing
end

"""
    soilhydrology_restart!(sh, bounds_col; flag="read")

Read/write soil hydrology state from/to restart file.

Ported from `soilhydrology_type%Restart` in `SoilHydrologyType.F90`.
Requires NetCDF/restart infrastructure — stub until that module is ported.
"""
function soilhydrology_restart!(sh::SoilHydrologyData,
                                 bounds_col::UnitRange{Int};
                                 flag::String = "read")
    # Stub: restart variable I/O will be added when restUtilMod/ncdio_pio is ported.
    # Variables that would be read/written:
    #   FROST_TABLE, ZWT, ZWT_PERCH
    return nothing
end

"""
    soilhydrology_read_nl!(sh; h2osfcflag=1)

Read namelist parameters for SoilHydrology.

Ported from `soilhydrology_type%ReadNL` in `SoilHydrologyType.F90`.
In Julia, namelist reading is replaced by keyword arguments.
"""
function soilhydrology_read_nl!(sh::SoilHydrologyData; h2osfcflag::Int = 1)
    sh.h2osfcflag = h2osfcflag
    return nothing
end
