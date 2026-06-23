# FatesPatchMod.jl
# Julia port of FATES src/fates/biogeochem/FatesPatchMod.F90
#
# The PATCH type — a FATES patch is a collection of cohorts that share a common
# disturbance history (and therefore age). A patch owns:
#   * the height-ordered cohort linked list (`tallest`/`shortest` head/tail
#     pointers into `fates_cohort_type` nodes connected by their `taller`/`shorter`
#     pointers),
#   * the older/younger patch linked-list pointers that order patches by age on a
#     site,
#   * the canopy-layer / PFT / leaf-layer leaf-area and radiation profile arrays
#     (tlai/elai/tsai_profile, canopy_area_profile, fabd/fabi/parsun/parsha, the
#     sun-fraction and PAR profiles),
#   * the per-element `litter` array (one `litter_type` per PARTEH element),
#   * the SPITFIRE `fuel` object,
#   * the two-stream radiation object `twostr`,
#   * the per-patch running means (24-hr veg temperature, leaf-photosynthesis
#     acclimation means, and the TRS seedling-layer means).
#
# Translation notes (per project conventions):
#   * `fates_r8` -> Float64; patch integers -> Int. NaN initializations preserved
#     as `NaN`; `fates_unset_int`/`fates_unset_r8` flushes follow Fortran exactly.
#   * `type(fates_cohort_type), pointer :: tallest/shortest` ->
#     `Union{fates_cohort_type,Nothing}` (the cohort linked-list head/tail).
#   * `type(fates_patch_type), pointer :: older/younger` ->
#     self-referential `Union{fates_patch_type,Nothing}` fields.
#   * `type(litter_type), pointer :: litter(:)` -> `Vector{litter_type}`.
#   * `type(fuel_type), pointer :: fuel` -> `Union{fuel_type,Nothing}` (allocated
#     in Create).
#   * `type(twostream_type) :: twostr` -> a `twostream_type` value (not a pointer).
#   * The polymorphic running-mean pointers (`class(rmean_type), pointer`) ->
#     `Union{rmean_type,Nothing}`; the per-PFT array pointers
#     (`class(rmean_arr_type), pointer :: (:)`) -> `Vector{rmean_arr_type}`.
#   * Allocatable profile arrays -> SoA Julia arrays sized at allocate time
#     (ReAllocateDynamics! / NanValues! allocate them; "not allocated" maps to a
#     zero-sized array, so the Fortran `allocated()` guard maps to `length(...) > 0`
#     / `!isempty(...)`).
#   * Type-bound procedures -> bang-functions dispatching on `fates_patch_type`,
#     preserving the Fortran names (Init!/ReAllocateDynamics!/NanValues!/
#     ZeroValues!/NanDynamics!/ZeroDynamics!/InitRunningMeans!/InitLitter!/Create!/
#     UpdateTreeGrassArea!/UpdateLiveGrass!/FreeMemory!/Dump/CheckVars).
#
# Deps: FatesConstantsMod (fates_unset_int/_r8, primaryland/secondaryland,
#   n_landuse_cats, TRS_regeneration, itrue/ifalse, nocomp_bareground,
#   N_DBH_BINS, N_DIST_TYPES, t_water_freeze_k_1atm), EDParamsMod
#   (nlevleaf, nclmax, maxpft), FatesRadiationMemMod (num_swb,
#   num_rad_stream_types), FatesInterfaceTypesMod (numpft), FatesCohortMod
#   (fates_cohort_type, FreeMemory), FatesRunningMeanMod (rmean_type,
#   rmean_arr_type, rmean_def_type, InitRMean!, moving_ema_window, fixed_window),
#   FatesLitterMod (litter_type, InitAllocate!/ZeroFlux!/InitConditions!/
#   DeallocateLitt!), FatesFuelMod (fuel_type, init_fuel!), TwoStreamMLPEMod
#   (twostream_type, DeallocTwoStream!), PRTGenericMod (num_elements,
#   element_list, carbon12_element, struct_organ/leaf_organ/sapw_organ, GetState),
#   PRTParametersMod (prt_params), FatesGlobals (fates_log).
# Standalone — NOT added to CLMInstances or any dual-copied struct.

# ---------------------------------------------------------------------------
# Patch-level running-mean definitions
# ---------------------------------------------------------------------------
# In Fortran these `rmean_def_type` instances are module globals in
# FatesRunningMeanMod, configured once at proc init from parameter-file
# timescales. That configuration step is not yet ported; we provide
# self-contained default definitions here so a patch's running means can be
# constructed and initialized (`InitRunningMeans!`). They can be re-pointed to
# the real definitions once the timescale-parameter wiring is ported.
const _patch_fixed_24hr   = rmean_def_type(mem_period = 86400.0, up_period = 1800.0, n_mem = 48, method = fixed_window)
const _patch_ema_lpa      = rmean_def_type(mem_period = 2.592e6, up_period = 86400.0, n_mem = 30, method = moving_ema_window)
const _patch_ema_longterm = rmean_def_type(mem_period = 3.1536e7, up_period = 86400.0, n_mem = 365, method = moving_ema_window)
const _patch_ema_sdlng_emerg_h2o = rmean_def_type(mem_period = 1.296e6, up_period = 86400.0, n_mem = 15, method = moving_ema_window)
const _patch_ema_sdlng_mort_par  = rmean_def_type(mem_period = 1.296e6, up_period = 86400.0, n_mem = 15, method = moving_ema_window)
const _patch_ema_sdlng2sap_par   = rmean_def_type(mem_period = 6.3072e6, up_period = 86400.0, n_mem = 73, method = moving_ema_window)
const _patch_ema_sdlng_mdd       = rmean_def_type(mem_period = 1.296e6, up_period = 86400.0, n_mem = 15, method = moving_ema_window)

# ===========================================================================
# FATES PATCH TYPE
# ===========================================================================

"""
    fates_patch_type

The per-patch state of a FATES grid cell. Field names are preserved from the
Fortran `fates_patch_type` for traceability. Scalar reals default to `NaN` and
integers to `fates_unset_int`, matching the Fortran NanValues initialization; the
`tallest`/`shortest`/`older`/`younger` pointers default to `nothing`, the
allocatable profile arrays default to zero-sized (i.e. "not allocated"), and the
`litter`/`fuel`/running-mean container objects default to empty/`nothing`.

Construct an empty patch and then call [`Create!`](@ref) (or [`Init!`](@ref) for a
bare allocation) to populate it. The `twostr` two-stream radiation object is held
by value and is default-constructed empty.
"""
Base.@kwdef mutable struct fates_patch_type

    # POINTERS --------------------------------------------------------------
    tallest::Union{fates_cohort_type,Nothing}  = nothing  # pointer to patch's tallest cohort
    shortest::Union{fates_cohort_type,Nothing} = nothing  # pointer to patch's shortest cohort
    older::Union{fates_patch_type,Nothing}     = nothing  # pointer to next older patch
    younger::Union{fates_patch_type,Nothing}   = nothing  # pointer to next younger patch

    # INDICES ---------------------------------------------------------------
    patchno::Int          = fates_unset_int  # unique number given to each new patch
    nocomp_pft_label::Int = fates_unset_int  # nocomp PFT ID for the patch (0 = bareground)

    # PATCH INFO ------------------------------------------------------------
    age::Float64                          = NaN              # average patch age [years]
    age_class::Int                        = fates_unset_int  # age class (history binning)
    area::Float64                         = NaN              # patch area [m2]
    countcohorts::Int                     = fates_unset_int  # number of cohorts in patch
    ncl_p::Int                            = fates_unset_int  # number of occupied canopy layers
    land_use_label::Int                   = fates_unset_int  # land use classification label
    age_since_anthro_disturbance::Float64 = NaN              # secondary-forest age since last anthro disturbance [years]
    changed_landuse_this_ts::Bool         = false           # flag: patch just underwent land use change

    # RUNNING MEANS ---------------------------------------------------------
    tveg24::Union{rmean_type,Nothing}               = nothing  # 24-hr mean veg temperature [K]
    tveg_lpa::Union{rmean_type,Nothing}             = nothing  # leaf-photosynthesis-acclimation mean veg T [K]
    tveg_longterm::Union{rmean_type,Nothing}        = nothing  # long-term mean veg T (T_home) [K]
    seedling_layer_par24::Union{rmean_type,Nothing} = nothing  # 24-hr mean PAR at seedling layer [W/m2]
    sdlng_emerg_smp::Vector{rmean_arr_type}         = rmean_arr_type[]  # per-PFT seedling-emergence soil matric potential mean
    sdlng_mort_par::Union{rmean_type,Nothing}       = nothing  # mean PAR at seedling layer (mortality timescale)
    sdlng_mdd::Vector{rmean_arr_type}               = rmean_arr_type[]  # per-PFT seedling moisture-deficit-days mean
    sdlng2sap_par::Union{rmean_type,Nothing}        = nothing  # mean PAR at seedling layer (sdlng->sap timescale)

    # LEAF ORGANIZATION -----------------------------------------------------
    pft_agb_profile::Matrix{Float64}  = fill(NaN, maxpft, N_DBH_BINS)  # binned AGB for patch fusion [kgC/m2]
    canopy_layer_tlai::Vector{Float64} = fill(NaN, nclmax)             # total LAI of each canopy layer [m2/m2]
    total_canopy_area::Float64 = NaN  # area covered by vegetation [m2]
    total_tree_area::Float64   = NaN  # area covered by woody vegetation [m2]
    total_grass_area::Float64  = NaN  # area covered by non-woody vegetation [m2]
    zstar::Float64             = NaN  # height of smallest canopy tree (strict-PPA mode) [m]

    # exposed/total leaf+stem area in each canopy layer, pft, and leaf layer
    elai_profile::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # exposed leaf area
    esai_profile::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # exposed stem area
    tlai_profile::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # total leaf area (incl. under snow)
    tsai_profile::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # total stem area (incl. under snow)
    canopy_area_profile::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # fraction of crown area per canopy area

    canopy_mask::Matrix{Int} = fill(fates_unset_int, nclmax, maxpft)  # is any of this pft in this canopy layer?
    nrad::Matrix{Int}        = fill(fates_unset_int, nclmax, maxpft)  # number of exposed veg layers
    nleaf::Matrix{Int}       = fill(fates_unset_int, nclmax, maxpft)  # number of total leaf layers
    c_stomata::Float64       = NaN  # mean stomatal conductance of all leaves [umol/m2/s]
    c_lblayer::Float64       = NaN  # mean boundary-layer conductance of all leaves [umol/m2/s]

    nrmlzd_parprof_pft_dir_z::Array{Float64,4} = Array{Float64,4}(undef, 0, 0, 0, 0)
    nrmlzd_parprof_pft_dif_z::Array{Float64,4} = Array{Float64,4}(undef, 0, 0, 0, 0)

    # RADIATION -------------------------------------------------------------
    rad_error::Vector{Float64}   = fill(NaN, num_swb)  # radiation consv error by band [W/m2]
    fcansno::Float64             = NaN                 # fraction of canopy covered in snow [0-1]
    solar_zenith_flag::Bool      = false               # daylight flag (based on zenith angle)
    solar_zenith_angle::Float64  = NaN                 # solar zenith angle [radians]
    gnd_alb_dif::Vector{Float64} = fill(NaN, num_swb)  # ground albedo for diffuse rad [0-1]
    gnd_alb_dir::Vector{Float64} = fill(NaN, num_swb)  # ground albedo for direct rad [0-1]

    # organized by canopy layer, pft, and leaf layer
    fabd_sun_z::Array{Float64,3}  = Array{Float64,3}(undef, 0, 0, 0)  # sun fraction of direct light absorbed
    fabd_sha_z::Array{Float64,3}  = Array{Float64,3}(undef, 0, 0, 0)  # shade fraction of direct light absorbed
    fabi_sun_z::Array{Float64,3}  = Array{Float64,3}(undef, 0, 0, 0)  # sun fraction of indirect light absorbed
    fabi_sha_z::Array{Float64,3}  = Array{Float64,3}(undef, 0, 0, 0)  # shade fraction of indirect light absorbed
    ed_parsun_z::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # PAR absorbed in the sun [W/m2]
    ed_parsha_z::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # PAR absorbed in the shade [W/m2]
    f_sun::Array{Float64,3}       = Array{Float64,3}(undef, 0, 0, 0)  # fraction of leaves in the sun [0-1]
    ed_laisun_z::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # sunlit LAI by layer
    ed_laisha_z::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # shaded LAI by layer

    # radiation profiles for comparison against observations
    parprof_pft_dir_z::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # direct-beam PAR profile [W/m2]
    parprof_pft_dif_z::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # diffuse PAR profile [W/m2]

    tr_soil_dir::Vector{Float64}     = Float64[]  # frac incoming direct -> soil as direct [0-1] (by numSWB)
    tr_soil_dif::Vector{Float64}     = Float64[]  # frac incoming diffuse -> soil as diffuse [0-1]
    tr_soil_dir_dif::Vector{Float64} = Float64[]  # frac incoming direct -> soil as diffuse [0-1]
    fab::Vector{Float64}             = Float64[]  # frac incoming total radiation absorbed by canopy
    fabd::Vector{Float64}            = Float64[]  # frac incoming direct radiation absorbed by canopy
    fabi::Vector{Float64}            = Float64[]  # frac incoming diffuse radiation absorbed by canopy
    sabs_dir::Vector{Float64}        = Float64[]  # frac incoming direct radiation absorbed (sabs)
    sabs_dif::Vector{Float64}        = Float64[]  # frac incoming diffuse radiation absorbed (sabs)

    # Two-stream data structures
    twostr::twostream_type = twostream_type()  # holds all two-stream data and procedures

    # ROOTS -----------------------------------------------------------------
    btran_ft::Vector{Float64}       = fill(NaN, maxpft)  # btran per PFT
    bstress_sal_ft::Vector{Float64} = fill(NaN, maxpft)  # salinity bstress per PFT

    # EXTERNAL SEED RAIN ----------------------------------------------------
    nitr_repro_stoich::Vector{Float64} = fill(NaN, maxpft)  # N:C ratio of a new recruit
    phos_repro_stoich::Vector{Float64} = fill(NaN, maxpft)  # P:C ratio of a new recruit

    # DISTURBANCE -----------------------------------------------------------
    disturbance_rates::Vector{Float64}       = fill(NaN, N_DIST_TYPES)   # disturbance rate [0-1/day] by type
    landuse_transition_rates::Vector{Float64} = fill(NaN, n_landuse_cats)  # land use transition rate
    fract_ldist_not_harvested::Float64        = NaN  # frac logged canopy-tree area not harvested [0-1]

    # LITTER AND COARSE WOODY DEBRIS ----------------------------------------
    litter::Vector{litter_type}            = litter_type[]  # litter (one per element)
    fuel::Union{fuel_type,Nothing}         = nothing        # SPITFIRE fuel class
    fragmentation_scaler::Vector{Float64}  = Float64[]      # litter fragmentation scaler by soil layer [0-1]

    # FUELS AND FIRE --------------------------------------------------------
    livegrass::Float64 = NaN  # total aboveground grass biomass in patch [kgC/m2]
    ros_front::Float64 = NaN  # rate of forward spread of fire [m/min]
    ros_back::Float64  = NaN  # rate of backward spread of fire [m/min]
    tau_l::Float64     = NaN  # duration of lethal heating [min]
    fi::Float64        = NaN  # average fire intensity of flaming front [kW/m]
    fire::Int          = fates_unset_int  # is there a fire? [1=yes, 0=no]
    fd::Float64        = NaN  # fire duration [min]
    frac_burnt::Float64 = NaN # fraction of patch burnt by fire
    scorch_ht::Vector{Float64} = fill(NaN, maxpft)  # scorch height [m]
    tfc_ros::Float64   = NaN  # total intensity-relevant fuel consumed - no trunks [kgC/m2/day]
end

# ===========================================================================
# Patch routines (Fortran type-bound procedures)
# ===========================================================================

"""
    Init!(this::fates_patch_type, num_swb, num_levsoil)

Initialize a new patch: allocate the fixed-size radiation/fragmentation arrays,
then NaN everything and zero what must be zeroed. Mirrors the Fortran `Init`.
"""
function Init!(this::fates_patch_type, num_swb::Integer, num_levsoil::Integer)
    # allocate arrays
    this.tr_soil_dir          = Vector{Float64}(undef, num_swb)
    this.tr_soil_dif          = Vector{Float64}(undef, num_swb)
    this.tr_soil_dir_dif      = Vector{Float64}(undef, num_swb)
    this.fab                  = Vector{Float64}(undef, num_swb)
    this.fabd                 = Vector{Float64}(undef, num_swb)
    this.fabi                 = Vector{Float64}(undef, num_swb)
    this.sabs_dir             = Vector{Float64}(undef, num_swb)
    this.sabs_dif             = Vector{Float64}(undef, num_swb)
    this.fragmentation_scaler = Vector{Float64}(undef, num_levsoil)

    # initialize all values to nan
    NanValues!(this)

    # zero values that should be zeroed
    ZeroValues!(this)

    return nothing
end

"""
    ReAllocateDynamics!(this::fates_patch_type)

Allocate / re-allocate the potentially large dynamic patch arrays (leaf-area and
radiation profiles), sized by the current number of canopy layers (`ncl_p`) and
the maximum number of leaf layers (`maxval(nleaf)`). Re-allocates only when the
canopy-layer count changed or the leaf-layer count grew (or shrank a good bit),
to avoid churn. Mirrors the Fortran `ReAllocateDynamics`.
"""
function ReAllocateDynamics!(this::fates_patch_type)
    ncan = this.ncl_p
    nveg = maximum(this.nleaf)

    # Assume we will need to allocate, unless the arrays already exist at the
    # right size.
    re_allocate = true

    # If the large patch arrays are not new, decide whether to re-allocate them.
    if !isempty(this.elai_profile)
        prev_ncan = size(this.tlai_profile, 1)
        prev_nveg = size(this.tlai_profile, 3)

        # Re-allocate if the number of canopy layers changed, the number of
        # vegetation layers grew, or it shrank a good bit. Otherwise keep them.
        if prev_ncan != ncan || nveg > prev_nveg || nveg < prev_nveg - 2
            # (Julia is GC'd; just drop the references by re-allocating below.)
        else
            re_allocate = false
        end
    end

    if re_allocate
        # Add a little buffer to nveg so it doesn't need reallocating as much.
        nveg = nveg + 1

        this.tlai_profile        = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.tsai_profile        = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.elai_profile        = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.esai_profile        = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.canopy_area_profile = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.f_sun               = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.fabd_sun_z          = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.fabd_sha_z          = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.fabi_sun_z          = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.fabi_sha_z          = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.nrmlzd_parprof_pft_dir_z = Array{Float64,4}(undef, num_rad_stream_types, ncan, numpft[], nveg)
        this.nrmlzd_parprof_pft_dif_z = Array{Float64,4}(undef, num_rad_stream_types, ncan, numpft[], nveg)
        this.ed_parsun_z         = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.ed_parsha_z         = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.ed_laisun_z         = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.ed_laisha_z         = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.parprof_pft_dir_z   = Array{Float64,3}(undef, ncan, numpft[], nveg)
        this.parprof_pft_dif_z   = Array{Float64,3}(undef, ncan, numpft[], nveg)
    end

    return nothing
end

"""
    NanDynamics!(this::fates_patch_type)

Set the dynamic leaf-area/radiation profile arrays to NaN. Mirrors the Fortran
`NanDynamics`.
"""
function NanDynamics!(this::fates_patch_type)
    fill!(this.elai_profile, NaN)
    fill!(this.esai_profile, NaN)
    fill!(this.tlai_profile, NaN)
    fill!(this.tsai_profile, NaN)
    fill!(this.canopy_area_profile, NaN)
    fill!(this.nrmlzd_parprof_pft_dir_z, NaN)
    fill!(this.nrmlzd_parprof_pft_dif_z, NaN)

    fill!(this.fabd_sun_z, NaN)
    fill!(this.fabd_sha_z, NaN)
    fill!(this.fabi_sun_z, NaN)
    fill!(this.fabi_sha_z, NaN)
    fill!(this.ed_laisun_z, NaN)
    fill!(this.ed_laisha_z, NaN)
    fill!(this.ed_parsun_z, NaN)
    fill!(this.ed_parsha_z, NaN)
    fill!(this.f_sun, NaN)
    fill!(this.parprof_pft_dir_z, NaN)
    fill!(this.parprof_pft_dif_z, NaN)

    return nothing
end

"""
    NanValues!(this::fates_patch_type)

Set all (fixed-size) values in the patch to NaN / unset so they are not used
before being defined. Pointers are nulled (`nothing`). Mirrors the Fortran
`NanValues`.
"""
function NanValues!(this::fates_patch_type)
    # set pointers to null
    this.tallest  = nothing
    this.shortest = nothing
    this.older    = nothing
    this.younger  = nothing

    # INDICES
    this.patchno          = fates_unset_int
    this.nocomp_pft_label = fates_unset_int

    # PATCH INFO
    this.age                          = NaN
    this.age_class                    = fates_unset_int
    this.area                         = NaN
    this.countcohorts                 = fates_unset_int
    this.ncl_p                        = fates_unset_int
    this.land_use_label               = fates_unset_int
    this.age_since_anthro_disturbance = NaN

    # LEAF ORGANIZATION
    fill!(this.pft_agb_profile, NaN)
    fill!(this.canopy_layer_tlai, NaN)
    this.total_canopy_area = NaN
    this.total_tree_area   = NaN
    this.total_grass_area  = NaN
    this.zstar             = NaN

    fill!(this.canopy_mask, fates_unset_int)
    fill!(this.nrad, fates_unset_int)
    fill!(this.nleaf, fates_unset_int)
    this.c_stomata = NaN
    this.c_lblayer = NaN

    # RADIATION
    fill!(this.rad_error, NaN)
    this.fcansno            = NaN
    this.solar_zenith_flag  = false
    this.solar_zenith_angle = NaN
    fill!(this.gnd_alb_dif, NaN)
    fill!(this.gnd_alb_dir, NaN)

    fill!(this.tr_soil_dir, NaN)
    fill!(this.tr_soil_dif, NaN)
    fill!(this.tr_soil_dir_dif, NaN)
    fill!(this.fab, NaN)
    fill!(this.fabd, NaN)
    fill!(this.fabi, NaN)
    fill!(this.sabs_dir, NaN)
    fill!(this.sabs_dif, NaN)

    # ROOTS
    fill!(this.btran_ft, NaN)
    fill!(this.bstress_sal_ft, NaN)

    # EXTERNAL SEED RAIN
    fill!(this.nitr_repro_stoich, NaN)
    fill!(this.phos_repro_stoich, NaN)

    # DISTURBANCE
    fill!(this.disturbance_rates, NaN)
    this.fract_ldist_not_harvested = NaN

    # LAND USE
    fill!(this.landuse_transition_rates, NaN)

    # LITTER AND COARSE WOODY DEBRIS
    fill!(this.fragmentation_scaler, NaN)

    # FUELS AND FIRE
    this.livegrass  = NaN
    this.ros_front  = NaN
    this.ros_back   = NaN
    this.tau_l      = NaN
    this.fi         = NaN
    this.fire       = fates_unset_int
    this.fd         = NaN
    fill!(this.scorch_ht, NaN)
    this.tfc_ros    = NaN
    this.frac_burnt = NaN

    return nothing
end

"""
    ZeroDynamics!(this::fates_patch_type)

Zero the dynamic leaf-area / radiation profile arrays (values that get
incremented). Mirrors the Fortran `ZeroDynamics`.
"""
function ZeroDynamics!(this::fates_patch_type)
    fill!(this.f_sun, 0.0)
    fill!(this.fabd_sun_z, 0.0)
    fill!(this.fabi_sun_z, 0.0)
    fill!(this.fabd_sha_z, 0.0)
    fill!(this.fabi_sha_z, 0.0)
    fill!(this.nrmlzd_parprof_pft_dir_z, 0.0)
    fill!(this.nrmlzd_parprof_pft_dif_z, 0.0)

    fill!(this.elai_profile, 0.0)
    fill!(this.esai_profile, 0.0)
    fill!(this.tlai_profile, 0.0)
    fill!(this.tsai_profile, 0.0)
    fill!(this.canopy_area_profile, 0.0)

    fill!(this.ed_laisun_z, 0.0)
    fill!(this.ed_laisha_z, 0.0)
    fill!(this.ed_parsun_z, 0.0)
    fill!(this.ed_parsha_z, 0.0)
    fill!(this.parprof_pft_dir_z, 0.0)
    fill!(this.parprof_pft_dif_z, 0.0)

    return nothing
end

"""
    ZeroValues!(this::fates_patch_type)

Zero the specific (incremented) patch variables so other uninitialized variables
can still be caught by their NaNs. Mirrors the Fortran `ZeroValues`.
"""
function ZeroValues!(this::fates_patch_type)
    # LEAF ORGANIZATION
    fill!(this.canopy_layer_tlai, 0.0)
    this.total_tree_area  = 0.0
    this.total_grass_area = 0.0
    this.zstar            = 0.0

    this.c_stomata = 0.0
    this.c_lblayer = 0.0

    # RADIATION
    fill!(this.rad_error, 0.0)
    fill!(this.tr_soil_dir_dif, 0.0)
    fill!(this.fab, 0.0)
    fill!(this.fabi, 0.0)
    fill!(this.fabd, 0.0)
    fill!(this.sabs_dir, 0.0)
    fill!(this.sabs_dif, 0.0)

    # ROOTS
    fill!(this.btran_ft, 0.0)

    # DISTURBANCE
    fill!(this.disturbance_rates, 0.0)
    this.fract_ldist_not_harvested = 0.0

    # LAND USE
    fill!(this.landuse_transition_rates, 0.0)

    # LITTER AND COARSE WOODY DEBRIS
    fill!(this.fragmentation_scaler, 0.0)

    # FIRE
    this.livegrass  = 0.0
    this.ros_front  = 0.0
    this.ros_back   = 0.0
    this.tau_l      = 0.0
    this.fi         = 0.0
    this.fd         = 0.0
    fill!(this.scorch_ht, 0.0)
    this.tfc_ros    = 0.0
    this.frac_burnt = 0.0

    return nothing
end

"""
    InitRunningMeans!(this::fates_patch_type, current_tod, regeneration_model, numpft)

Set initial values for the patch running means (24-hr veg temperature and the
leaf-photosynthesis-acclimation means; under the TRS regeneration model also the
seedling-layer means). Mirrors the Fortran `InitRunningMeans`. The running-mean
*definitions* used here are the module-local defaults
([`_patch_fixed_24hr`](@ref) etc.).
"""
function InitRunningMeans!(this::fates_patch_type, current_tod::Integer,
                           regeneration_model::Integer, numpft::Integer)
    # PARAMETERS
    temp_init_veg     = 15.0 + t_water_freeze_k_1atm  # default veg temp [K]
    init_seedling_par = 5.0       # arbitrary seedling-layer init [MJ m-2 d-1]
    init_seedling_smp = -26652.0  # arbitrary smp init [mm]

    this.tveg24        = rmean_type()
    this.tveg_lpa      = rmean_type()
    this.tveg_longterm = rmean_type()

    # set initial values for running means
    InitRMean!(this.tveg24, _patch_fixed_24hr; init_value = temp_init_veg,
               init_offset = Float64(current_tod))
    InitRMean!(this.tveg_lpa, _patch_ema_lpa; init_value = temp_init_veg)
    InitRMean!(this.tveg_longterm, _patch_ema_longterm; init_value = temp_init_veg)

    if regeneration_model == TRS_regeneration
        this.seedling_layer_par24 = rmean_type()
        this.sdlng_mdd            = [rmean_arr_type() for _ in 1:numpft]
        this.sdlng_emerg_smp      = [rmean_arr_type() for _ in 1:numpft]
        this.sdlng_mort_par       = rmean_type()
        this.sdlng2sap_par        = rmean_type()

        InitRMean!(this.seedling_layer_par24, _patch_fixed_24hr;
                   init_value = init_seedling_par, init_offset = Float64(current_tod))
        InitRMean!(this.sdlng_mort_par, _patch_ema_sdlng_mort_par;
                   init_value = temp_init_veg)
        InitRMean!(this.sdlng2sap_par, _patch_ema_sdlng2sap_par;
                   init_value = init_seedling_par)

        for pft in 1:numpft
            this.sdlng_mdd[pft].p       = rmean_type()
            this.sdlng_emerg_smp[pft].p = rmean_type()

            InitRMean!(this.sdlng_mdd[pft].p, _patch_ema_sdlng_mdd; init_value = 0.0)
            InitRMean!(this.sdlng_emerg_smp[pft].p, _patch_ema_sdlng_emerg_h2o;
                       init_value = init_seedling_smp)
        end
    end

    return nothing
end

"""
    InitLitter!(this::fates_patch_type, num_pft, num_levsoil)

Set initial values for the per-element litter pools: allocate one `litter_type`
per PARTEH element, allocate its pool arrays, zero the fluxes, and set the
prognostic pools to the unset sentinel. Mirrors the Fortran `InitLitter`.
"""
function InitLitter!(this::fates_patch_type, num_pft::Integer, num_levsoil::Integer)
    this.litter = [litter_type() for _ in 1:num_elements[]]

    for el in 1:num_elements[]
        InitAllocate!(this.litter[el], num_pft, num_levsoil, element_list[el])
        ZeroFlux!(this.litter[el])
        InitConditions!(this.litter[el], fates_unset_r8, fates_unset_r8,
                        fates_unset_r8, fates_unset_r8, fates_unset_r8, fates_unset_r8)
    end

    return nothing
end

"""
    Create!(this, age, area, land_use_label, nocomp_pft, num_swb, num_pft,
            num_levsoil, current_tod, regeneration_model)

Create a new patch with input and default values: initialize (NaN then zero),
initialize the running means, the per-element litter, and the fuel object, mark
the two-stream object as un-initialized (so the radiation module will allocate
it), then assign the known patch attributes (age/area/land-use/nocomp label,
transmitted-radiation defaults, one canopy layer). Mirrors the Fortran `Create`.
"""
function Create!(this::fates_patch_type, age::Real, area::Real,
                 land_use_label::Integer, nocomp_pft::Integer, num_swb::Integer,
                 num_pft::Integer, num_levsoil::Integer, current_tod::Integer,
                 regeneration_model::Integer)
    # initialize patch (sets all values to nan, then some values to zero)
    Init!(this, num_swb, num_levsoil)

    # initialize running means for patch
    InitRunningMeans!(this, current_tod, regeneration_model, num_pft)

    # initialize litter
    InitLitter!(this, num_pft, num_levsoil)

    # initialize fuel
    this.fuel = fuel_type()
    init_fuel!(this.fuel)

    # Leave the two-stream scattering elements un-allocated; the radiation module
    # checks this and will then initialize and allocate (Fortran nulls
    # twostr%scelg here; our default twostream_type starts with an empty scelg).
    this.twostr = twostream_type()

    # assign known patch attributes
    this.age       = age
    this.age_class = 1
    this.area      = area

    # assign anthropogenic disturbance category and label
    this.land_use_label = land_use_label
    if land_use_label == secondaryland
        this.age_since_anthro_disturbance = age
    else
        this.age_since_anthro_disturbance = fates_unset_r8
    end
    this.nocomp_pft_label = nocomp_pft

    fill!(this.tr_soil_dir, 1.0)
    fill!(this.tr_soil_dif, 1.0)
    this.ncl_p = 1

    this.changed_landuse_this_ts = false

    return nothing
end

"""
    UpdateTreeGrassArea!(this::fates_patch_type)

Calculate and update the total tree area and grass area (summed over crown areas
of woody / non-woody cohorts, capped at the patch area). Mirrors the Fortran
`UpdateTreeGrassArea`. Skips bareground patches.
"""
function UpdateTreeGrassArea!(this::fates_patch_type)
    if this.nocomp_pft_label != nocomp_bareground
        tree_area  = 0.0
        grass_area = 0.0

        currentCohort = this.tallest
        while currentCohort !== nothing
            if prt_params.woody[currentCohort.pft] == itrue
                tree_area += currentCohort.c_area
            else
                grass_area += currentCohort.c_area
            end
            currentCohort = currentCohort.shorter
        end

        this.total_tree_area  = min(tree_area, this.area)
        this.total_grass_area = min(grass_area, this.area)
    end

    return nothing
end

"""
    UpdateLiveGrass!(this::fates_patch_type)

Calculate the sum of live grass biomass [kgC/m2] on the patch: for each non-woody
cohort, sum the leaf + sapwood + structural carbon, scale by number density and
divide by patch area. Mirrors the Fortran `UpdateLiveGrass`.
"""
function UpdateLiveGrass!(this::fates_patch_type)
    live_grass = 0.0
    currentCohort = this.tallest
    while currentCohort !== nothing
        # for grasses sum all aboveground tissues
        if prt_params.woody[currentCohort.pft] == ifalse
            live_grass += (GetState(currentCohort.prt, leaf_organ, carbon12_element) +
                           GetState(currentCohort.prt, sapw_organ, carbon12_element) +
                           GetState(currentCohort.prt, struct_organ, carbon12_element)) *
                          currentCohort.n / this.area
        end
        currentCohort = currentCohort.shorter
    end

    this.livegrass = live_grass

    return nothing
end

"""
    FreeMemory!(this::fates_patch_type, regeneration_model, numpft)

Release the dynamic memory and objects held within the patch (cohorts, two-stream
scattering elements, litter, fuel, profile arrays, running means). Julia is
garbage-collected, so we delegate to the contained objects' deallocators and then
drop the references (reset arrays to empty / pointers to `nothing`). Mirrors the
Fortran `FreeMemory` (which does NOT deallocate the patch structure itself).
"""
function FreeMemory!(this::fates_patch_type, regeneration_model::Integer,
                     numpft::Integer)
    # first free the cohorts (walk shortest -> taller, freeing each)
    ccohort = this.shortest
    while ccohort !== nothing
        ncohort = ccohort.taller
        FreeMemory(ccohort)
        ccohort = ncohort
    end
    this.tallest  = nothing
    this.shortest = nothing

    # Deallocate radiation scattering elements (if allocated)
    if !isempty(this.twostr.scelg)
        DeallocTwoStream!(this.twostr)
    end

    # deallocate all litter objects
    for el in 1:num_elements[]
        DeallocateLitt!(this.litter[el])
    end
    this.litter = litter_type[]

    # deallocate fuel
    this.fuel = nothing

    # deallocate the allocatable arrays
    this.tr_soil_dir          = Float64[]
    this.tr_soil_dif          = Float64[]
    this.tr_soil_dir_dif      = Float64[]
    this.fab                  = Float64[]
    this.fabd                 = Float64[]
    this.fabi                 = Float64[]
    this.sabs_dir             = Float64[]
    this.sabs_dif             = Float64[]
    this.fragmentation_scaler = Float64[]

    # The dynamic profile arrays may not have been allocated (patches can be
    # spawned and destroyed before EDCanopyStructureMod sizes them); guard.
    if !isempty(this.elai_profile)
        this.tlai_profile        = Array{Float64,3}(undef, 0, 0, 0)
        this.tsai_profile        = Array{Float64,3}(undef, 0, 0, 0)
        this.elai_profile        = Array{Float64,3}(undef, 0, 0, 0)
        this.esai_profile        = Array{Float64,3}(undef, 0, 0, 0)
        this.f_sun               = Array{Float64,3}(undef, 0, 0, 0)
        this.fabd_sun_z          = Array{Float64,3}(undef, 0, 0, 0)
        this.fabd_sha_z          = Array{Float64,3}(undef, 0, 0, 0)
        this.fabi_sun_z          = Array{Float64,3}(undef, 0, 0, 0)
        this.fabi_sha_z          = Array{Float64,3}(undef, 0, 0, 0)
        this.nrmlzd_parprof_pft_dir_z = Array{Float64,4}(undef, 0, 0, 0, 0)
        this.nrmlzd_parprof_pft_dif_z = Array{Float64,4}(undef, 0, 0, 0, 0)
        this.ed_parsun_z         = Array{Float64,3}(undef, 0, 0, 0)
        this.ed_parsha_z         = Array{Float64,3}(undef, 0, 0, 0)
        this.ed_laisun_z         = Array{Float64,3}(undef, 0, 0, 0)
        this.ed_laisha_z         = Array{Float64,3}(undef, 0, 0, 0)
        this.parprof_pft_dir_z   = Array{Float64,3}(undef, 0, 0, 0)
        this.parprof_pft_dif_z   = Array{Float64,3}(undef, 0, 0, 0)
        this.canopy_area_profile = Array{Float64,3}(undef, 0, 0, 0)
    end

    # deallocate running means
    this.tveg24        = nothing
    this.tveg_lpa      = nothing
    this.tveg_longterm = nothing

    if regeneration_model == TRS_regeneration
        this.seedling_layer_par24 = nothing
        this.sdlng_mort_par       = nothing
        this.sdlng2sap_par        = nothing
        this.sdlng_mdd            = rmean_arr_type[]
        this.sdlng_emerg_smp      = rmean_arr_type[]
    end

    return nothing
end

"""
    Dump(this::fates_patch_type)

Print out the (scalar) attributes of a patch and the per-element litter pool
sums, for diagnostics (omitting the large arrays). Mirrors the Fortran `Dump`.
"""
function Dump(this::fates_patch_type)
    log = fates_log()
    println(log, "----------------------------------------")
    println(log, " Dumping Patch Information              ")
    println(log, " (omitting arrays)                      ")
    println(log, "----------------------------------------")
    println(log, "pa%patchno            = ", this.patchno)
    println(log, "pa%age                = ", this.age)
    println(log, "pa%age_class          = ", this.age_class)
    println(log, "pa%area               = ", this.area)
    println(log, "pa%countcohorts       = ", this.countcohorts)
    println(log, "pa%ncl_p              = ", this.ncl_p)
    println(log, "pa%total_canopy_area  = ", this.total_canopy_area)
    println(log, "pa%total_tree_area    = ", this.total_tree_area)
    println(log, "pa%total_grass_area   = ", this.total_grass_area)
    println(log, "pa%zstar              = ", this.zstar)
    println(log, "pa%solar_zenith_flag  = ", this.solar_zenith_flag)
    println(log, "pa%solar_zenith_angle = ", this.solar_zenith_angle)
    println(log, "pa%gnd_alb_dif        = ", this.gnd_alb_dif)
    println(log, "pa%gnd_alb_dir        = ", this.gnd_alb_dir)
    println(log, "pa%c_stomata          = ", this.c_stomata)
    println(log, "pa%c_lblayer          = ", this.c_lblayer)
    println(log, "pa%disturbance_rates  = ", this.disturbance_rates)
    println(log, "pa%land_use_label     = ", this.land_use_label)
    println(log, "----------------------------------------")

    for el in 1:num_elements[]
        println(log, "element id: ", element_list[el])
        println(log, "seed mass: ", sum(this.litter[el].seed))
        println(log, "seed germ mass: ", sum(this.litter[el].seed_germ))
        println(log, "leaf fines(pft): ", sum(this.litter[el].leaf_fines))
        println(log, "root fines(pft,sl): ", sum(this.litter[el].root_fines))
        println(log, "ag_cwd(c): ", sum(this.litter[el].ag_cwd))
        println(log, "bg_cwd(c,sl): ", sum(this.litter[el].bg_cwd))
    end

    return nothing
end

"""
    CheckVars(this::fates_patch_type, var_aliases::AbstractString) -> Int

Perform numerical checks on patch (and contained-cohort) variables of interest.
`var_aliases` is a colon-separated registry string ('VAR1:VAR2:...'). Returns 0
if all checks pass, nonzero (and dumps the offending patch/cohort) otherwise.
Mirrors the Fortran `CheckVars` (`return_code` becomes the return value).
"""
function CheckVars(this::fates_patch_type, var_aliases::AbstractString)
    return_code = 0

    if occursin("co_n", var_aliases)
        currentCohort = this.shortest
        while currentCohort !== nothing
            return_code = _check_var_real(currentCohort.n, "cohort%n")
            if return_code != 0
                Dump(this)
                Dump(currentCohort)
                return return_code
            end
            currentCohort = currentCohort.taller
        end
    end

    if occursin("co_dbh", var_aliases)
        currentCohort = this.shortest
        while currentCohort !== nothing
            return_code = _check_var_real(currentCohort.dbh, "cohort%dbh")
            if return_code != 0
                Dump(this)
                Dump(currentCohort)
                return return_code
            end
            currentCohort = currentCohort.taller
        end
    end

    if occursin("pa_area", var_aliases)
        return_code = _check_var_real(this.area, "patch%area")
        if return_code != 0
            Dump(this)
            return return_code
        end
    end

    return return_code
end

# Local numerical check mirroring FatesUtilsMod::check_var_real:
#   0  = fine
#   1  = NaN detected
#   10 = overflow (Inf)
#   100 = underflow (denormal)
function _check_var_real(rval::Real, varname::AbstractString)
    if isnan(rval)
        @warn "check_var_real: NaN detected in $varname"
        return 1
    elseif isinf(rval)
        @warn "check_var_real: overflow detected in $varname"
        return 10
    elseif rval != 0.0 && abs(rval) < floatmin(Float64)
        @warn "check_var_real: underflow detected in $varname"
        return 100
    end
    return 0
end
