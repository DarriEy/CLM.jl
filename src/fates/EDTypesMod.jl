# EDTypesMod.jl
# Julia port of FATES src/fates/main/EDTypesMod.F90
#
# The TOP of the FATES type system. Defines `ed_site_type` — the FATES
# site/column-level structure that owns the patch linked list (oldest_patch /
# youngest_patch) and ALL site-level state: site-level plant hydraulics object,
# fire-weather object, the disturbance / mortality / recruitment accumulators,
# demotion/promotion rates, the site-level mass-balance + flux diagnostics,
# phenology state, running means, and the per-size/pft history accumulation
# arrays. Also defines the smaller site-scoped helper types
# (ed_resources_management_type, elem_diag_type, site_ifluxbal_type,
# site_fluxdiags_type, site_massbal_type) and the module-level FATES
# constants/parameters that live here.
#
# Translation notes (per project conventions):
#   * `fates_r8` -> Float64; integers -> Int. Allocatable arrays -> SoA Julia
#     arrays defaulting to zero-sized ("not allocated"); fixed-dimension Fortran
#     arrays (e.g. `dstatus(maxpft)`, `vegtemp_memory(num_vegtemp_mem)`,
#     `wood_product_harvest(maxpft)`) -> sized Julia arrays.
#   * `type(fates_patch_type), pointer :: oldest_patch/youngest_patch` ->
#     `Union{fates_patch_type,Nothing}` (the patch age-ordered linked-list heads).
#   * `type(ed_site_hydr_type), pointer :: si_hydr` ->
#     `Union{ed_site_hydr_type,Nothing}` (allocatable site hydraulics object).
#   * `class(fire_weather), pointer :: fireWeather` ->
#     `Union{fire_weather,Nothing}` (the abstract fire-weather object).
#   * The pointer-array members (mass_balance, iflux_balance, flux_diags%elem)
#     -> `Vector{...}` (allocated per-element); default zero-sized.
#   * Type-bound procedures -> bang-functions dispatching on the owning type,
#     preserving the Fortran names (ZeroFluxDiags!, ZeroMassBalState!,
#     ZeroMassBalFlux!, CalculateTreeGrassAreaSite!, dump_site,
#     get_current_landuse_statevector, get_secondary_young_fraction).
#
# Constant-collision note: the following constants used by EDTypesMod are
# ALREADY merged elsewhere and are referenced here, NOT redefined:
#   - FatesConstantsMod: itrue, ifalse, nocomp_bareground, secondaryland,
#     secondary_age_threshold, nearzero, days_per_year, fates_unset_r8,
#     n_landuse_cats, N_DIST_TYPES, N_DBH_BINS, nocomp_bareground_land
#   - EDParamsMod: nclmax, nlevleaf, maxpft
#   - FatesLitterMod: ncwd
#   - PRTGenericMod: num_elements (a Ref{Int})
#   - FatesInterfaceTypesMod: numpft, nlevsclass, hlm_parteh_mode (Ref{Int})
# Everything below (init_recruit_trim, n_rad_stream_types, idirect/idiffuse,
# do_fates_salinity, init_spread_*, area/area_inv, numWaterMem, numlevsoil_max,
# num_vegtemp_mem, the phen_*stat flags, patch/cohort fusion + termination
# thresholds, homogenize_seed_pfts) is NEW and declared by EDTypesMod itself.
# Standalone — NOT added to CLMInstances or any dual-copied struct.

# ---------------------------------------------------------------------------
# Module-level FATES constants/parameters owned by EDTypesMod
# ---------------------------------------------------------------------------

# This is the initial trimming value that new recruits start with
const init_recruit_trim = 0.8

# Radiation parameters (one option only, so they live here for now)
const n_rad_stream_types = 2  # number of radiation streams (direct/diffuse)
const idirect  = 1            # array index for direct radiation
const idiffuse = 2            # array index for diffuse radiation

# Flag to turn on/off salinity effects on the effective "btran" stress function
const do_fates_salinity = false

# Community-level amount of spread expected in nearly-bare-ground / inventory
# starting modes. 1 => crown areas at maximum spread; 0 => least spread.
const init_spread_near_bare_ground = 1.0
const init_spread_inventory        = 0.0

# MODEL PARAMETERS
const area     = 10000.0  # Notional area of simulated forest [m2]
const area_inv = 1.0e-4   # Inverse of the notional area (faster math)

const numWaterMem    = 10  # water memory saved as site level var
const numlevsoil_max = 30  # scratch space for static arrays (max # soil layers)

# BIOLOGY / BIOGEOCHEMISTRY
const num_vegtemp_mem = 10  # window over which we track temp for cold senescence (days)

# Phenology status flag definitions (cold type is cstat, dry type is dstat)
const phen_cstat_nevercold = 0  # has not experienced a cold period, leaves dropped, non-cold region
const phen_cstat_iscold    = 1  # in a cold-state where leaves should have fallen
const phen_cstat_notcold   = 2  # in a warm-state where leaves are allowed to flush

const phen_dstat_timeoff  = 0  # leaves off due to time exceedance (drought phenology)
const phen_dstat_moistoff = 1  # leaves off due to moisture avail  (drought phenology)
const phen_dstat_moiston  = 2  # leaves on due to moisture avail    (drought phenology)
const phen_dstat_timeon   = 3  # leaves on due to time exceedance   (drought phenology)
const phen_dstat_pshed    = 4  # leaves partially abscissing        (drought phenology)

# PATCH FUSION
const force_patchfuse_min_biomass = 0.005  # min biomass [kg/m2] below which to force-fuse patches
const patch_fusion_tolerance_relaxation_increment = 1.1  # increment for patch fusion threshold
const max_age_of_second_oldest_patch = 200.0  # age [yr] above which to combine all patches

# COHORT FUSION
const HEIGHTMAX     = 30.0  # max dbh value used in hgt profile comparison
const N_HEIGHT_BINS = 60    # no. of height bins used to distribute LAI

# COHORT TERMINATION
const min_npm2       = 1.0e-7   # minimum cohort number density per m2 before termination
const min_patch_area = 0.01     # smallest allowable patch area before termination
const min_patch_area_forced = 0.0001  # patch termination will not fuse youngest patch below this
const min_nppatch    = min_npm2 * min_patch_area  # minimum number of cohorts per patch
const min_n_safemath = 1.0e-12  # immediately remove super small number densities to prevent FPEs

# special mode to cause PFTs to create seed mass of all currently-existing PFTs
const homogenize_seed_pfts = false

# ===========================================================================
# Resources management type (YX)
# ===========================================================================

"""
    ed_resources_management_type

Site-level resources-management accounting (KgC/site). Field names preserved
from the Fortran `ed_resources_management_type`.
"""
Base.@kwdef mutable struct ed_resources_management_type
    trunk_product_site::Float64  = 0.0  # actual trunk product at site level [KgC/site]
    harvest_debt::Float64        = 0.0  # kgC per site not successfully harvested
    harvest_debt_sec::Float64    = 0.0  # kgC per site from secondary patches not harvested
    # debug variables
    delta_litter_stock::Float64  = 0.0  # [kgC/site = kgC/ha]
    delta_biomass_stock::Float64 = 0.0  # [kgC/site]
    delta_individual::Float64    = 0.0
end

# ===========================================================================
# Element flux diagnostics type
# ===========================================================================

"""
    elem_diag_type

Per-element (C/N/P) diagnostics of fluxes. Acts as an intermediary to write
fluxes to history after number densities of plants have changed, and to rebuild
the history flux diagnostics on restart. Litter fluxes are the total from living
plant turnover, non-disturbance mortality mass transfer, and disturbance
mortality mass transfer [kg/ha/day]. Field names preserved from Fortran.
"""
Base.@kwdef mutable struct elem_diag_type
    cwd_ag_input::Vector{Float64} = zeros(ncwd)  # [kg/ha/day]
    cwd_bg_input::Vector{Float64} = zeros(ncwd)  # [kg/ha/day]
    surf_fine_litter_input::Vector{Float64} = Float64[]  # allocatable
    root_litter_input::Vector{Float64}      = Float64[]  # allocatable

    tot_seed_turnover::Float64 = 0.0  # decay of living seed bank to fragmented litter [kg/m2/day]
    exported_harvest::Float64  = 0.0  # mass of harvested veg exported, not sent to litter [kg/m2/day]
    burned_liveveg::Float64    = 0.0  # mass burned from living plants [kg/m2/day]

    # Integrated Error Terms ( Int. Flux - State )
    err_liveveg::Float64 = 0.0  # error from comparing [state-integrated flux] in live veg [kg/m2]
    err_litter::Float64  = 0.0  # net change in litter [kg/m2]
end

# ===========================================================================
# Integrated mass-flux balance type
# ===========================================================================

"""
    site_ifluxbal_type

Tracks integrated daily flux vs instantaneous state for live vegetation and
litter, so that state == initial condition + integrated fluxes over the run.
Field names preserved from Fortran.
"""
Base.@kwdef mutable struct site_ifluxbal_type
    state_liveveg::Float64 = 0.0  # assessed instantaneously [kg/m2]
    iflux_liveveg::Float64 = 0.0  # integrated daily         [kg/m2]
    state_litter::Float64  = 0.0  # assessed instantaneously [kg/m2]
    iflux_litter::Float64  = 0.0  # integrated daily         [kg/m2]
end

# ===========================================================================
# Site flux diagnostics type
# ===========================================================================

"""
    site_fluxdiags_type

Site flux diagnostics used to write history output based on fluxes calculated
before mortality and cohort/patch restructuring but written after patch
structure. `elem` is the per-element array of [`elem_diag_type`]. The remaining
scalars are uniform over all elements; the `*_scpf` arrays are size x pft
delineated nutrient flux diagnostics (allocated only when complex diagnostics +
the species of interest are requested). Field names preserved from Fortran.
"""
Base.@kwdef mutable struct site_fluxdiags_type
    elem::Vector{elem_diag_type} = elem_diag_type[]  # type(elem_diag_type),pointer :: elem(:)

    npp::Float64         = 0.0  # kg m-2 day-1

    # Nutrient Flux Diagnostics
    resp_excess::Float64 = 0.0  # plant C respired due to C overflow [kg m-2 s-1]
    nh4_uptake::Float64  = 0.0  # plant nh4 uptake [kg m-2 s-1]
    no3_uptake::Float64  = 0.0  # plant no3 uptake [kg m-2 s-1]
    sym_nfix::Float64    = 0.0  # plant N uptake via symbiotic fixation [kg m-2 s-1]
    n_efflux::Float64    = 0.0  # efflux of unusable N to soil labile pool [kg m-2 s-1]
    p_uptake::Float64    = 0.0  # po4 uptake [kg m-2 s-1]
    p_efflux::Float64    = 0.0  # efflux of unusable P to soil labile pool [kg m-2 s-1]

    # Size by PFT delineated nutrient flux diagnostics (allocatable)
    nh4_uptake_scpf::Vector{Float64} = Float64[]
    no3_uptake_scpf::Vector{Float64} = Float64[]
    sym_nfix_scpf::Vector{Float64}   = Float64[]
    n_efflux_scpf::Vector{Float64}   = Float64[]
    p_uptake_scpf::Vector{Float64}   = Float64[]
    p_efflux_scpf::Vector{Float64}   = Float64[]
end

# ===========================================================================
# Site mass-balance type
# ===========================================================================

"""
    site_massbal_type

Per-element accounting type used to ensure mass is not lost or created. One unit
of "site" is nominally equivalent to 1 hectare. These mass checks are for
INCREMENTAL checks during the dynamics step. Field names preserved from Fortran.
"""
Base.@kwdef mutable struct site_massbal_type
    old_stock::Float64 = 0.0  # remember biomass stock from last time [Kg/site]
    err_fates::Float64 = 0.0  # total mass balance error for FATES processes [kg/site]

    # Group 3: Components of the total site level mass fluxes
    gpp_acc::Float64         = 0.0  # accumulated gross primary productivity [kg/site/day]
    aresp_acc::Float64       = 0.0  # accumulated autotrophic respiration [kg/site/day]
    net_root_uptake::Float64 = 0.0  # net uptake of C or nutrients through roots [kg/site/day]
    seed_in::Float64         = 0.0  # total mass of external seed rain into site [kg/site/day]
    seed_out::Float64        = 0.0  # total mass of seeds exported outside site [kg/site/day]
    frag_out::Float64        = 0.0  # litter and CWD fragmentation flux [kg/site/day]
    wood_product_harvest::Vector{Float64}       = zeros(maxpft)  # exported wood product from harvest [kg/site/day]
    wood_product_landusechange::Vector{Float64} = zeros(maxpft)  # exported wood product from land-use change [kg/site/day]
    burn_flux_to_atm::Float64 = 0.0  # total mass burned and exported to atmosphere [kg/site/day]
    flux_generic_in::Float64  = 0.0  # prescribed/artificial input fluxes + init [kg/site/day]
    flux_generic_out::Float64 = 0.0  # prescribed/artificial output fluxes [kg/site/day]
    patch_resize_err::Float64 = 0.0  # mass gained/lost due to re-sizing patches [kg/site/day]
end

# ===========================================================================
# SITE TYPE
# ===========================================================================

"""
    ed_site_type

The FATES site (nominally 1 hectare / a column). The TOP of the FATES type
system. Owns the age-ordered patch linked list via `oldest_patch` /
`youngest_patch` (`fates_patch_type` nodes connected by their `older`/`younger`
pointers), the site-level plant hydraulics object `si_hydr`, the fire-weather
object `fireWeather`, the per-element mass-balance / integrated-flux-balance
arrays, the flux diagnostics, the full phenology state, the fire driver state,
soil layering, and the termination/recruitment/demotion/promotion/disturbance
accumulators (many on size x pft history arrays). Field names preserved from the
Fortran `ed_site_type` for traceability.

Construct an empty site, then attach a patch linked list and (optionally) call
the zero/init helpers ([`ZeroMassBalState!`](@ref), [`ZeroMassBalFlux!`](@ref),
[`ZeroFluxDiags!`](@ref)). The `*_scpf`/size-x-pft and soil-layer arrays default
to zero-sized ("not allocated"); fixed-dimension arrays
(`dstatus`/`vegtemp_memory`/`recruitment_rate`/`liqvol_memory`/`smp_memory`) are
sized at construction.
"""
Base.@kwdef mutable struct ed_site_type
    # POINTERS — the patch age-ordered linked list
    oldest_patch::Union{fates_patch_type,Nothing}   = nothing  # oldest patch at the site
    youngest_patch::Union{fates_patch_type,Nothing} = nothing  # youngest patch at the site

    # Resource management
    resources_management::ed_resources_management_type = ed_resources_management_type()

    # Global index of this site in the history output file
    h_gid::Int = fates_unset_int

    # INDICES
    lat::Float64 = NaN  # latitude:  degrees
    lon::Float64 = NaN  # longitude: degrees

    # Fixed Biogeography mode inputs
    area_PFT::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # area to PFTs by land use class [ha/ha non-bareground]
    area_bareground::Float64  = NaN  # area to bare ground in nocomp (HLM PFT 0) [ha/ha]
    use_this_pft::Vector{Int} = Int[]  # is area_PFT > 0 ? (1=yes, 0=no)

    area_by_age::Vector{Float64} = Float64[]  # total area of patches in each age bin [m2]

    # Nutrient relevant
    rec_l2fr::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # running mean l2fr for new recruits, pft x canopy_layer
    ema_npp::Float64 = NaN  # exponential moving average of NPP [gC/m2/year]

    # Two-stream scratch arrays
    omega_2str::Matrix{Float64}   = Matrix{Float64}(undef, 0, 0)  # matrix inverted in two-stream solve
    taulambda_2str::Vector{Float64} = Float64[]  # two-stream linear-system coefficients (tau / lambda)
    ipiv_2str::Vector{Int}        = Int[]        # pivot indices for the lapack 2str solver

    # SP mode target PFT level variables
    sp_tlai::Vector{Float64} = Float64[]  # target TLAI per FATES pft
    sp_tsai::Vector{Float64} = Float64[]  # target TSAI per FATES pft
    sp_htop::Vector{Float64} = Float64[]  # target HTOP per FATES pft

    # Instantaneous Mass Balance (one per element)
    mass_balance::Vector{site_massbal_type} = site_massbal_type[]

    # Integrated Mass Balance checks (live vegetation + litter)
    iflux_balance::Vector{site_ifluxbal_type} = site_ifluxbal_type[]

    # Flux diagnostics
    flux_diags::site_fluxdiags_type = site_fluxdiags_type()

    # PHENOLOGY
    grow_deg_days::Float64 = NaN  # phenology growing degree days
    snow_depth::Float64    = NaN  # site-level snow depth (used for ELAI/TLAI calcs)

    cstatus::Int = fates_unset_int  # leaves on/off for cold decid (see phen_cstat_*)
    dstatus::Vector{Int} = fill(fates_unset_int, maxpft)  # leaves on/off for drought decid (see phen_dstat_*)
    nchilldays::Int = fates_unset_int  # num chilling days (botta gdd threshold calc)
    ncolddays::Int  = fates_unset_int  # num cold days (must exceed threshold to drop leaves)
    vegtemp_memory::Vector{Float64} = zeros(num_vegtemp_mem)  # last 10 days temperature for senescence [degC]
    cleafondate::Int   = fates_unset_int  # model date (day int) of leaf on  (cold)
    cleafoffdate::Int  = fates_unset_int  # model date (day int) of leaf off (cold)
    cndaysleafon::Int  = fates_unset_int  # days since leaf on period started  (cold)
    cndaysleafoff::Int = fates_unset_int  # days since leaf off period started (cold)
    dleafondate::Vector{Int}   = fill(fates_unset_int, maxpft)  # model date of leaf on  (drought)
    dleafoffdate::Vector{Int}  = fill(fates_unset_int, maxpft)  # model date of leaf off (drought)
    dndaysleafon::Vector{Int}  = fill(fates_unset_int, maxpft)  # days since leaf on  (drought)
    dndaysleafoff::Vector{Int} = fill(fates_unset_int, maxpft)  # days since leaf off (drought)
    elong_factor::Vector{Float64} = fill(NaN, maxpft)  # elongation factor (ED2-like phenology) [0..1]
    phen_model_date::Int = fates_unset_int  # current model date (day int), continuous over restarts

    liqvol_memory::Matrix{Float64} = zeros(numWaterMem, maxpft)  # last 10 days soil liquid water vol (drought)
    smp_memory::Matrix{Float64}    = zeros(numWaterMem, maxpft)  # last 10 days soil matric potential (drought)

    # FIRE
    wind::Float64 = NaN  # daily wind [m/min] for Spitfire units
    fdi::Float64  = NaN  # daily probability an ignition event starts a fire
    NF::Float64   = NaN  # daily ignitions in km2
    NF_successful::Float64 = NaN  # daily ignitions in km2 that actually lead to fire
    fireWeather::Union{fire_weather,Nothing} = nothing  # fire weather object

    # PLANT HYDRAULICS
    si_hydr::Union{ed_site_hydr_type,Nothing} = nothing

    # Soil Layering
    nlevsoil::Int = fates_unset_int  # number of soil layers in this site
    zi_soil::Vector{Float64} = Float64[]  # interface level below a "z" level [m] (zero index for surface)
    dz_soil::Vector{Float64} = Float64[]  # layer thickness [m]
    z_soil::Vector{Float64}  = Float64[]  # layer depth [m]
    rootfrac_scr::Vector{Float64} = Float64[]  # scratch space for root fractions (NOT thread-safe)

    # DIAGNOSTICS — TERMINATION, RECRUITMENT, DEMOTION, DISTURBANCE
    term_crownarea_canopy::Float64 = NaN  # crownarea from termination mortality, per canopy level
    term_crownarea_ustory::Float64 = NaN  # crownarea from termination mortality, per canopy level
    imort_crownarea::Float64       = NaN  # crownarea killed by impact mortality per year [m2 day]
    fmort_crownarea_canopy::Float64 = NaN # crownarea of canopy indivs killed by fire per year [m2/sec]
    fmort_crownarea_ustory::Float64 = NaN # crownarea of understory indivs killed by fire per year [m2/sec]

    term_nindivs_canopy::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # canopy indivs terminated, type x size x pft
    term_nindivs_ustory::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # understory indivs terminated, type x size x pft

    term_carbonflux_canopy::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # C flux live->dead, termination, type x pft [kgC/ha/day]
    term_carbonflux_ustory::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # C flux live->dead, termination, type x pft [kgC/ha/day]
    imort_carbonflux::Vector{Float64}        = Float64[]  # biomass killed by impact mortality, by pft [kgC/m2/sec]
    fmort_carbonflux_canopy::Vector{Float64} = Float64[]  # biomass of canopy indivs killed by fire, by pft [gC/m2/sec]
    fmort_carbonflux_ustory::Vector{Float64} = Float64[]  # biomass of understory indivs killed by fire, by pft [gC/m2/sec]

    term_abg_flux::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)  # aboveground biomass lost to termination, size x pft
    imort_abg_flux::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # aboveground biomass lost to impact mortality, size x pft [kgC/m2/sec]
    fmort_abg_flux::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # aboveground biomass lost to fire mortality, size x pft

    demotion_carbonflux::Float64  = NaN  # biomass of demoted individuals canopy->understory [kgC/ha/day]
    promotion_carbonflux::Float64 = NaN  # biomass of promoted individuals understory->canopy [kgC/ha/day]
    recruitment_rate::Vector{Float64} = zeros(maxpft)  # individuals recruited into new cohorts
    demotion_rate::Vector{Float64}  = Float64[]  # rate demoted canopy->understory per FATES timestep
    promotion_rate::Vector{Float64} = Float64[]  # rate promoted understory->canopy per FATES timestep
    imort_rate::Matrix{Float64}     = Matrix{Float64}(undef, 0, 0)  # rate killed by impact mortality per year, size x pft

    fmort_rate_canopy::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # rate canopy indivs killed by fire per year, size x pft
    fmort_rate_ustory::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # rate understory indivs killed by fire per year, size x pft
    fmort_rate_cambial::Matrix{Float64} = Matrix{Float64}(undef, 0, 0) # rate killed by fire from cambial damage per year, size x pft
    fmort_rate_crown::Matrix{Float64}   = Matrix{Float64}(undef, 0, 0) # rate killed by fire from crown damage per year, size x pft

    imort_rate_damage::Array{Float64,3}        = Array{Float64,3}(undef, 0, 0, 0)  # indivs per damage class dying from impact mortality
    term_nindivs_canopy_damage::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # indivs per damage class, termination - canopy
    term_nindivs_ustory_damage::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # indivs per damage class, termination - ustory
    fmort_rate_canopy_damage::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # indivs per damage class dying from fire - canopy
    fmort_rate_ustory_damage::Array{Float64,3} = Array{Float64,3}(undef, 0, 0, 0)  # indivs per damage class dying from fire - ustory
    fmort_cflux_canopy_damage::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # cflux per damage class dying from fire - canopy
    fmort_cflux_ustory_damage::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # cflux per damage class dying from fire - ustory
    imort_cflux_damage::Matrix{Float64}        = Matrix{Float64}(undef, 0, 0)  # C flux from impact mortality by damage class [kgC/m2/sec]
    term_cflux_canopy_damage::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)  # C flux from termination mortality by damage class
    term_cflux_ustory_damage::Matrix{Float64}  = Matrix{Float64}(undef, 0, 0)  # C flux from termination mortality by damage class

    growthflux_fusion::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # rate indivs moving into a size class via fusion/day, size x pft

    crownarea_canopy_damage::Float64 = NaN  # crown area of canopy damaged annually
    crownarea_ustory_damage::Float64 = NaN  # crown area of understory damaged annually

    # Canopy Spread
    spread::Float64 = NaN  # dynamic canopy allometric term [unitless]

    # Seed dispersal
    seed_out::Vector{Float64} = Float64[]  # amount of seed leaving the site [kg/site/day]
    seed_in::Vector{Float64}  = Float64[]  # amount of seed dispersed into the site [kg/site/day]

    # Disturbance rates (actual and "potential")
    disturbance_rates::Array{Float64,3} = zeros(N_DIST_TYPES, n_landuse_cats, n_landuse_cats)  # actual disturbance rates [m2/m2/day]
    primary_land_patchfusion_error::Float64 = NaN  # error in total area of primary patches from fusion [m2/m2/day]
    landuse_transition_matrix::Matrix{Float64} = zeros(n_landuse_cats, n_landuse_cats)  # land use transition matrix [m2/m2/year]

    min_allowed_landuse_fraction::Float64 = NaN  # min land-use type below which patches would be too small [m2/m2]
    landuse_vector_gt_min::Vector{Bool}   = Bool[]  # is the land use state vector > the minimum?
    transition_landuse_from_off_to_on::Bool = false  # restart flag: triggers land use init procedure
end

# ===========================================================================
# Type-bound / module procedures (site init / zero / diagnostics)
# ===========================================================================

"""
    ZeroFluxDiags!(this::site_fluxdiags_type)

Zero the site flux diagnostics (all per-element litter inputs and the uniform
nutrient flux scalars/arrays). `gpp_prev_scpf` is intentionally NOT zeroed in
Fortran (assigned at the end of the daily history write), so it is omitted here.
"""
function ZeroFluxDiags!(this::site_fluxdiags_type)
    for el in 1:num_elements[]
        e = this.elem[el]
        fill!(e.cwd_ag_input, 0.0)
        fill!(e.cwd_bg_input, 0.0)
        isempty(e.surf_fine_litter_input) || fill!(e.surf_fine_litter_input, 0.0)
        isempty(e.root_litter_input) || fill!(e.root_litter_input, 0.0)
        e.burned_liveveg    = 0.0
        e.tot_seed_turnover = 0.0
        e.exported_harvest  = 0.0
        e.err_liveveg       = 0.0
        e.err_litter        = 0.0
    end

    this.npp         = 0.0
    this.resp_excess = 0.0
    this.nh4_uptake  = 0.0
    this.no3_uptake  = 0.0
    this.sym_nfix    = 0.0
    this.n_efflux    = 0.0
    this.p_uptake    = 0.0
    this.p_efflux    = 0.0

    isempty(this.nh4_uptake_scpf) || fill!(this.nh4_uptake_scpf, 0.0)
    isempty(this.no3_uptake_scpf) || fill!(this.no3_uptake_scpf, 0.0)
    isempty(this.sym_nfix_scpf)   || fill!(this.sym_nfix_scpf, 0.0)
    isempty(this.n_efflux_scpf)   || fill!(this.n_efflux_scpf, 0.0)
    isempty(this.p_uptake_scpf)   || fill!(this.p_uptake_scpf, 0.0)
    isempty(this.p_efflux_scpf)   || fill!(this.p_efflux_scpf, 0.0)
    return nothing
end

"""
    ZeroMassBalState!(this::site_massbal_type)

Zero the per-element mass-balance state (the remembered stock and accumulated
error).
"""
function ZeroMassBalState!(this::site_massbal_type)
    this.old_stock = 0.0
    this.err_fates = 0.0
    return nothing
end

"""
    ZeroMassBalFlux!(this::site_massbal_type)

Zero the per-element mass-balance flux components.
"""
function ZeroMassBalFlux!(this::site_massbal_type)
    this.gpp_acc         = 0.0
    this.aresp_acc       = 0.0
    this.net_root_uptake = 0.0
    this.seed_in         = 0.0
    this.seed_out        = 0.0
    this.frag_out        = 0.0
    fill!(this.wood_product_harvest, 0.0)
    fill!(this.wood_product_landusechange, 0.0)
    this.burn_flux_to_atm = 0.0
    this.flux_generic_in  = 0.0
    this.flux_generic_out = 0.0
    this.patch_resize_err = 0.0
    return nothing
end

"""
    CalculateTreeGrassAreaSite!(csite::ed_site_type)
        -> (tree_fraction, grass_fraction, bare_fraction)

Calculate total tree, grass, and bare fractions for a site by traversing the
patch linked list (oldest -> youngest) and accumulating each patch's
tree/grass area (normalized by [`area`](@ref)). Bareground patches
(`nocomp_pft_label == nocomp_bareground`) are skipped. Mirrors the Fortran
subroutine, which calls `UpdateTreeGrassArea` on each non-bareground patch.
"""
function CalculateTreeGrassAreaSite!(csite::ed_site_type)
    tree_fraction = 0.0
    grass_fraction = 0.0

    currentPatch = csite.oldest_patch
    while currentPatch !== nothing
        if currentPatch.nocomp_pft_label != nocomp_bareground
            UpdateTreeGrassArea!(currentPatch)
            tree_fraction  += currentPatch.total_tree_area / area
            grass_fraction += currentPatch.total_grass_area / area
        end
        currentPatch = currentPatch.younger
    end

    # if cover > 1.0, grasses are under the trees
    grass_fraction = min(grass_fraction, 1.0 - tree_fraction)
    bare_fraction = 1.0 - tree_fraction - grass_fraction

    return tree_fraction, grass_fraction, bare_fraction
end

"""
    dump_site(csite::ed_site_type)

Print the site coordinates (latitude / longitude) to the FATES log. Mirrors the
Fortran `dump_site` diagnostic.
"""
function dump_site(csite::ed_site_type)
    log = fates_log()
    println(log, "----------------------------------------")
    println(log, " Site Coordinates                       ")
    println(log, "----------------------------------------")
    println(log, "latitude                    = ", csite.lat)
    println(log, "longitude                   = ", csite.lon)
    println(log, "----------------------------------------")
    return nothing
end

"""
    get_current_landuse_statevector(this::ed_site_type) -> Vector{Float64}

Calculate how much of a site is each land use category (length
`n_landuse_cats`). Does not include bare ground when nocomp + fixed biogeography
is on (so will not sum to one in that case; otherwise it does). Traverses the
patch list oldest -> youngest.
"""
function get_current_landuse_statevector(this::ed_site_type)
    current_state_vector = zeros(n_landuse_cats)

    currentPatch = this.oldest_patch
    while currentPatch !== nothing
        if currentPatch.land_use_label > nocomp_bareground_land
            current_state_vector[currentPatch.land_use_label] +=
                currentPatch.area / area
        end
        currentPatch = currentPatch.younger
    end

    return current_state_vector
end

"""
    get_secondary_young_fraction(this::ed_site_type) -> Float64

Calculate how much of the secondary area is "young" (below
[`secondary_age_threshold`](@ref)). Returns -1 if there is no secondary patch
area at all. Traverses the patch list oldest -> youngest.
"""
function get_secondary_young_fraction(this::ed_site_type)
    secondary_young_area = 0.0
    secondary_old_area = 0.0

    currentPatch = this.oldest_patch
    while currentPatch !== nothing
        if currentPatch.land_use_label == secondaryland
            if currentPatch.age >= secondary_age_threshold
                secondary_old_area += currentPatch.area
            else
                secondary_young_area += currentPatch.area
            end
        end
        currentPatch = currentPatch.younger
    end

    if (secondary_young_area + secondary_old_area) > nearzero
        secondary_young_fraction = secondary_young_area /
            (secondary_young_area + secondary_old_area)
    else
        secondary_young_fraction = -1.0
    end

    return secondary_young_fraction
end
