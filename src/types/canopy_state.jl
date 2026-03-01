# ==========================================================================
# Ported from: src/biogeophys/CanopyStateType.F90
# Canopy state data type allocation and initialization
# ==========================================================================

"""
    CanopyStateData

Canopy state data structure. Holds all canopy state variables at the patch
level, including LAI, SAI, canopy heights, roughness lengths, sunlit/shaded
fractions, vegetation water potentials, and accumulated quantities.

Ported from `canopystate_type` in `CanopyStateType.F90`.
"""
Base.@kwdef mutable struct CanopyStateData
    # --- Integer patch-level fields ---
    frac_veg_nosno_patch     ::Vector{Int}     = Int[]       # patch fraction of vegetation not covered by snow (0 OR 1) [-]
    frac_veg_nosno_alb_patch ::Vector{Int}     = Int[]       # patch fraction of vegetation not covered by snow (0 OR 1) [-]

    # --- LAI/SAI (patch-level, 1D) ---
    tlai_patch               ::Vector{Float64} = Float64[]   # patch canopy one-sided leaf area index, no burying by snow
    tsai_patch               ::Vector{Float64} = Float64[]   # patch canopy one-sided stem area index, no burying by snow
    elai_patch               ::Vector{Float64} = Float64[]   # patch canopy one-sided leaf area index with burying by snow
    esai_patch               ::Vector{Float64} = Float64[]   # patch canopy one-sided stem area index with burying by snow

    # --- SP mode history fields ---
    tlai_hist_patch          ::Vector{Float64} = Float64[]   # patch canopy one-sided leaf area index, for SP mode
    tsai_hist_patch          ::Vector{Float64} = Float64[]   # patch canopy one-sided stem area index, for SP mode
    htop_hist_patch          ::Vector{Float64} = Float64[]   # patch canopy height, for SP mode

    # --- Sunlit/shaded LAI (patch-level) ---
    elai240_patch            ::Vector{Float64} = Float64[]   # patch canopy one-sided LAI with burying by snow avg over 10days
    laisun_patch             ::Vector{Float64} = Float64[]   # patch sunlit projected leaf area index
    laisha_patch             ::Vector{Float64} = Float64[]   # patch shaded projected leaf area index
    mlaidiff_patch           ::Vector{Float64} = Float64[]   # patch difference between lai month one and month two

    # --- Sunlit/shaded LAI by canopy layer (patch-level, 2D) ---
    laisun_z_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch sunlit leaf area for canopy layer (np, nlevcan)
    laisha_z_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch shaded leaf area for canopy layer (np, nlevcan)

    # --- Monthly LAI (patch-level, 2D) ---
    annlai_patch             ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch 12 months of monthly lai (np, 12)

    # --- Biomass (patch-level) ---
    stem_biomass_patch       ::Vector{Float64} = Float64[]   # Aboveground stem biomass (kg/m**2)
    leaf_biomass_patch       ::Vector{Float64} = Float64[]   # Aboveground leaf biomass (kg/m**2)

    # --- Canopy geometry (patch-level) ---
    htop_patch               ::Vector{Float64} = Float64[]   # patch canopy top (m)
    hbot_patch               ::Vector{Float64} = Float64[]   # patch canopy bottom (m)
    z0m_patch                ::Vector{Float64} = Float64[]   # patch momentum roughness length (m)
    displa_patch             ::Vector{Float64} = Float64[]   # patch displacement height (m)

    # --- Sunlit fraction (patch-level) ---
    fsun_patch               ::Vector{Float64} = Float64[]   # patch sunlit fraction of canopy
    fsun24_patch             ::Vector{Float64} = Float64[]   # patch 24hr average of sunlit fraction of canopy
    fsun240_patch            ::Vector{Float64} = Float64[]   # patch 240hr average of sunlit fraction of canopy

    # --- Leaf properties (patch-level) ---
    dleaf_patch              ::Vector{Float64} = Float64[]   # patch characteristic leaf width (diameter) [m]
    rscanopy_patch           ::Vector{Float64} = Float64[]   # patch canopy stomatal resistance (s/m) (ED specific)

    # --- Vegetation water potential (patch-level, 2D) ---
    vegwp_patch              ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch vegetation water matric potential (mm) (np, nvegwcs)
    vegwp_ln_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch vegetation water matric potential at local noon (mm) (np, nvegwcs)
    vegwp_pd_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # patch predawn vegetation water matric potential (mm) (np, nvegwcs)

    # --- Namelist parameter ---
    leaf_mr_vcm              ::Float64 = SPVAL               # Scalar constant of leaf respiration with Vcmax
end

"""
    canopystate_init!(cs::CanopyStateData, np::Int)

Allocate and initialize all fields of a `CanopyStateData` instance for
`np` patches. Integer fields are initialized to `typemax(Int)` or `0`,
real fields to `NaN`, matching the Fortran `InitAllocate`.

Ported from `canopystate_type%InitAllocate` in `CanopyStateType.F90`.
"""
function canopystate_init!(cs::CanopyStateData, np::Int)
    nlevcan = NLEVCAN
    nvegwcs = NVEGWCS

    # --- Integer patch-level ---
    cs.frac_veg_nosno_patch     = fill(typemax(Int), np)
    cs.frac_veg_nosno_alb_patch = fill(0, np)

    # --- LAI/SAI 1D ---
    cs.tlai_hist_patch          = fill(NaN, np)
    cs.tsai_hist_patch          = fill(NaN, np)
    cs.htop_hist_patch          = fill(NaN, np)
    cs.tlai_patch               = fill(NaN, np)
    cs.tsai_patch               = fill(NaN, np)
    cs.elai_patch               = fill(NaN, np)
    cs.elai240_patch            = fill(NaN, np)
    cs.esai_patch               = fill(NaN, np)
    cs.laisun_patch             = fill(NaN, np)
    cs.laisha_patch             = fill(NaN, np)
    cs.mlaidiff_patch           = fill(NaN, np)

    # --- Sunlit/shaded LAI by canopy layer 2D ---
    cs.laisun_z_patch           = fill(NaN, np, nlevcan)
    cs.laisha_z_patch           = fill(NaN, np, nlevcan)

    # --- Monthly LAI 2D (Fortran: 12 × np, Julia: np × 12) ---
    cs.annlai_patch             = fill(NaN, np, 12)

    # --- Biomass ---
    cs.stem_biomass_patch       = fill(NaN, np)
    cs.leaf_biomass_patch       = fill(NaN, np)

    # --- Canopy geometry ---
    cs.htop_patch               = fill(NaN, np)
    cs.hbot_patch               = fill(NaN, np)
    cs.z0m_patch                = fill(NaN, np)
    cs.displa_patch             = fill(NaN, np)

    # --- Sunlit fraction ---
    cs.fsun_patch               = fill(NaN, np)
    cs.fsun24_patch             = fill(NaN, np)
    cs.fsun240_patch            = fill(NaN, np)

    # --- Leaf properties ---
    cs.dleaf_patch              = fill(NaN, np)
    cs.rscanopy_patch           = fill(NaN, np)

    # --- Vegetation water potential 2D ---
    cs.vegwp_patch              = fill(NaN, np, nvegwcs)
    cs.vegwp_ln_patch           = fill(NaN, np, nvegwcs)
    cs.vegwp_pd_patch           = fill(NaN, np, nvegwcs)

    return nothing
end

"""
    canopystate_clean!(cs::CanopyStateData)

Deallocate (reset to empty) all fields of a `CanopyStateData` instance.
"""
function canopystate_clean!(cs::CanopyStateData)
    # Integer vectors
    cs.frac_veg_nosno_patch     = Int[]
    cs.frac_veg_nosno_alb_patch = Int[]

    # Float64 vectors
    cs.tlai_patch               = Float64[]
    cs.tsai_patch               = Float64[]
    cs.elai_patch               = Float64[]
    cs.esai_patch               = Float64[]
    cs.tlai_hist_patch          = Float64[]
    cs.tsai_hist_patch          = Float64[]
    cs.htop_hist_patch          = Float64[]
    cs.elai240_patch            = Float64[]
    cs.laisun_patch             = Float64[]
    cs.laisha_patch             = Float64[]
    cs.mlaidiff_patch           = Float64[]
    cs.stem_biomass_patch       = Float64[]
    cs.leaf_biomass_patch       = Float64[]
    cs.htop_patch               = Float64[]
    cs.hbot_patch               = Float64[]
    cs.z0m_patch                = Float64[]
    cs.displa_patch             = Float64[]
    cs.fsun_patch               = Float64[]
    cs.fsun24_patch             = Float64[]
    cs.fsun240_patch            = Float64[]
    cs.dleaf_patch              = Float64[]
    cs.rscanopy_patch           = Float64[]

    # Matrices
    cs.laisun_z_patch           = Matrix{Float64}(undef, 0, 0)
    cs.laisha_z_patch           = Matrix{Float64}(undef, 0, 0)
    cs.annlai_patch             = Matrix{Float64}(undef, 0, 0)
    cs.vegwp_patch              = Matrix{Float64}(undef, 0, 0)
    cs.vegwp_ln_patch           = Matrix{Float64}(undef, 0, 0)
    cs.vegwp_pd_patch           = Matrix{Float64}(undef, 0, 0)

    return nothing
end

"""
    canopystate_init_cold!(cs::CanopyStateData, bounds_patch::UnitRange{Int};
                            landunit_patch::Union{Vector{Int},Nothing}=nothing,
                            lun_itype::Union{Vector{Int},Nothing}=nothing)

Initialize cold-start conditions for canopy state variables.
Sets LAI, SAI, canopy heights, biomass, vegetation water potential,
and sunlit/shaded fractions to their cold-start values.

Arguments:
- `cs`: CanopyStateData instance
- `bounds_patch`: UnitRange{Int} for patches (1:npatches)
- `landunit_patch`: landunit index for each patch (optional)
- `lun_itype`: landunit type for each landunit (optional)

Ported from `canopystate_type%InitCold` in `CanopyStateType.F90`.
"""
function canopystate_init_cold!(cs::CanopyStateData, bounds_patch::UnitRange{Int};
                                 landunit_patch::Union{Vector{Int},Nothing} = nothing,
                                 lun_itype::Union{Vector{Int},Nothing} = nothing)
    for p in bounds_patch
        cs.tlai_patch[p]         = 0.0
        cs.tsai_patch[p]         = 0.0
        cs.elai_patch[p]         = 0.0
        cs.esai_patch[p]         = 0.0
        cs.stem_biomass_patch[p] = 0.0
        cs.leaf_biomass_patch[p] = 0.0
        cs.htop_patch[p]         = 0.0
        cs.hbot_patch[p]         = 0.0
        cs.vegwp_patch[p, :]    .= -2.5e4

        # Sunlit/shaded LAI conditional on landunit type
        if landunit_patch !== nothing && lun_itype !== nothing
            l = landunit_patch[p]
            if lun_itype[l] == ISTSOIL || lun_itype[l] == ISTCROP
                cs.laisun_patch[p] = 0.0
                cs.laisha_patch[p] = 0.0
            end
        else
            # If no landunit info provided, initialize all
            cs.laisun_patch[p] = 0.0
            cs.laisha_patch[p] = 0.0
        end

        cs.tlai_hist_patch[p] = 0.0
        cs.tsai_hist_patch[p] = 0.0
        cs.htop_hist_patch[p] = 0.0

        # fsun must be initialized to spval for accumulation field averaging
        cs.fsun_patch[p] = SPVAL
    end

    return nothing
end

"""
    canopystate_read_nml!(cs::CanopyStateData; leaf_mr_vcm::Float64=0.015)

Set namelist parameters for canopy state.

Ported from `canopystate_type%ReadNML` in `CanopyStateType.F90`.
In Julia, namelist values are passed directly instead of reading from a file.
"""
function canopystate_read_nml!(cs::CanopyStateData; leaf_mr_vcm::Float64 = 0.015)
    cs.leaf_mr_vcm = leaf_mr_vcm
    return nothing
end

"""
    canopystate_set_nml_for_testing!(cs::CanopyStateData)

Set namelist parameters for unit-testing (leaf_mr_vcm = 0.015).

Ported from `canopystate_type%SetNMLForTesting` in `CanopyStateType.F90`.
"""
function canopystate_set_nml_for_testing!(cs::CanopyStateData)
    cs.leaf_mr_vcm = 0.015
    return nothing
end

# ==========================================================================
# The following subroutines depend on infrastructure modules that are not yet
# ported (history, restart/IO, accumulator). They are provided as stubs that
# document the Fortran interface and can be filled in when those modules
# become available.
# ==========================================================================

"""
    canopystate_init_history!(cs, bounds_patch)

Register canopy state fields for history file output.

Ported from `canopystate_type%InitHistory` in `CanopyStateType.F90`.
Requires history infrastructure (histFileMod) — stub until that module is ported.
"""
function canopystate_init_history!(cs::CanopyStateData,
                                    bounds_patch::UnitRange{Int})
    # Stub: history field registration will be added when histFileMod is ported.
    # Fields that would be registered:
    #   ELAI, ESAI, LAISUN, LAISHA, AGSB, AGLB, FSUN, HBOT, DISPLA, HTOP,
    #   TLAI, TSAI, Z0MV_DENSE, FSUN24, FSUN240, LAI240, RSCANOPY,
    #   VEGWP, VEGWPLN, VEGWPPD
    return nothing
end

"""
    canopystate_restart!(cs, bounds_patch; flag="read")

Read/write canopy state from/to restart file.

Ported from `canopystate_type%Restart` in `CanopyStateType.F90`.
Requires NetCDF/restart infrastructure — stub until that module is ported.
"""
function canopystate_restart!(cs::CanopyStateData,
                               bounds_patch::UnitRange{Int};
                               flag::String = "read")
    # Stub: restart variable I/O will be added when restUtilMod/ncdio_pio is ported.
    # Variables that would be read/written:
    #   FRAC_VEG_NOSNO_ALB, tlai, tsai, elai, esai, stem_biomass,
    #   leaf_biomass, htop, hbot, mlaidiff, fsun, vegwp, VEGWPLN, VEGWPPD
    return nothing
end

"""
    canopystate_init_acc_buffer!(cs, bounds_patch)

Initialize accumulation buffer for canopy state accumulated fields.

Ported from `canopystate_type%InitAccBuffer` in `CanopyStateType.F90`.
Requires accumulation infrastructure (accumulMod) — stub until that module is ported.
"""
function canopystate_init_acc_buffer!(cs::CanopyStateData,
                                       bounds_patch::UnitRange{Int})
    # Stub: accumulation field definitions will be added when accumulMod is ported.
    # Fields that would be initialized:
    #   FSUN24 (runmean, -1 day), FSUN240 (runmean, -10 days),
    #   LAI240 (runmean, -10 days)

    # Initialize values to SPVAL as in Fortran
    for p in bounds_patch
        cs.fsun24_patch[p]  = SPVAL
        cs.fsun240_patch[p] = SPVAL
        cs.elai240_patch[p] = SPVAL
    end

    return nothing
end

"""
    canopystate_init_acc_vars!(cs, bounds_patch)

Initialize accumulated variables from the accumulation buffer.
Called for both initial and restart runs.

Ported from `canopystate_type%InitAccVars` in `CanopyStateType.F90`.
Requires accumulation infrastructure (accumulMod) — stub until that module is ported.
"""
function canopystate_init_acc_vars!(cs::CanopyStateData,
                                     bounds_patch::UnitRange{Int})
    # Stub: When accumulMod is ported, this will:
    # 1. Extract FSUN24, FSUN240, LAI240 from accumulation buffer
    # 2. Copy extracted values into fsun24_patch, fsun240_patch, elai240_patch
    return nothing
end

"""
    canopystate_update_acc_vars!(cs, bounds_patch)

Update accumulated canopy state variables each timestep.
Handles fsun24, fsun240, and elai240 running means.

Ported from `canopystate_type%UpdateAccVars` in `CanopyStateType.F90`.
Core logic is ported; accumulator calls are stubs until accumulMod is ported.
"""
function canopystate_update_acc_vars!(cs::CanopyStateData,
                                       bounds_patch::UnitRange{Int})
    # When accumulMod is available, the logic is:
    #
    # 1. Copy fsun_patch into buffer, then:
    #    update_accum_field('FSUN24',  buffer, nstep)
    #    extract_accum_field('FSUN24', fsun24_patch, nstep)
    #    update_accum_field('FSUN240', buffer, nstep)
    #    extract_accum_field('FSUN240', fsun240_patch, nstep)
    #
    # 2. Copy elai_patch into buffer, then:
    #    update_accum_field('LAI240',  buffer, nstep)
    #    extract_accum_field('LAI240', elai240_patch, nstep)

    return nothing
end
