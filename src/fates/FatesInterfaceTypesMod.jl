# FatesInterfaceTypesMod.jl
# Julia port of FATES src/fates/main/FatesInterfaceTypesMod.F90
#
# The type definitions for the host land model (HLM) <-> FATES boundary:
#   - bc_in_type     boundary conditions IN  (host -> FATES)
#   - bc_out_type    boundary conditions OUT (FATES -> host)
#   - bc_pconst_type parameter constants (set once, never modified)
# plus the module-level control/dimension parameters carried across the
# interface (the Fortran `public` module variables) and the history-output
# dimension map vectors.
#
# Pure type definitions + their allocate/init. fates_r8 -> Float64. Fortran
# allocatable/pointer arrays -> Julia arrays (SoA, undef-0 defaults). The
# Fortran allocate routines become @kwdef defaults + an allocate/init
# bang-function sized from the relevant dimensions.
#
# Standalone — nothing here is added to CLMInstances or any dual-copied struct.
# Deps: FatesConstantsMod (fates_r8, itrue, ifalse, fates_unset_int).

# ===========================================================================
# Parameters dictated by the Host Land Model.
# THESE ARE NOT DYNAMIC. Set once during initialization. The Fortran module
# carries these as bare public module variables; here we hold them as Refs so
# they are mutable globals with module scope.
# ===========================================================================

const hlm_numSWb        = Ref{Int}(fates_unset_int)  # number of shortwave broad-bands (typ. 2, VIS/NIR)
const hlm_ivis          = Ref{Int}(fates_unset_int)  # array index of the visible band
const hlm_inir          = Ref{Int}(fates_unset_int)  # array index of the NIR band
const hlm_maxlevsoil    = Ref{Int}(fates_unset_int)  # max number of soil layers
const hlm_is_restart    = Ref{Int}(fates_unset_int)  # 1=restart, 0=not
const hlm_name          = Ref{String}("")            # HLM name, drives IO filtering (CLM/ALM/ATS)
const hlm_decomp        = Ref{String}("")            # decomposition scheme: CENTURY/MIMICS/CTC
const hlm_nu_com        = Ref{String}("")            # nutrient competition scheme: RD/ECA/NONE
const hlm_nitrogen_spec   = Ref{Int}(fates_unset_int)  # 0:none 1:nh4 2:nh4+no3
const hlm_phosphorus_spec = Ref{Int}(fates_unset_int)  # 0:none 1:p on
const hlm_stepsize      = Ref{Float64}(fates_unset_r8)  # HLM step-size (s)
const hlm_hio_ignore_val = Ref{Float64}(fates_unset_r8) # flush value ignored in history averages
const hlm_masterproc    = Ref{Int}(fates_unset_int)  # 1=master proc, 0=not
const hlm_ipedof        = Ref{Int}(fates_unset_int)  # HLM pedotransfer index (plant hydraulics)
const hlm_parteh_mode   = Ref{Int}(fates_unset_int)  # which PARTEH allocation hypothesis
const hlm_seeddisp_cadence = Ref{Int}(fates_unset_int) # 0:none 1:daily 2:monthly 3:yearly
const hlm_use_ch4       = Ref{Int}(fates_unset_int)  # methane model active?
const hlm_use_vertsoilc = Ref{Int}(fates_unset_int)  # vertically discretized soil carbon?
const hlm_spitfire_mode = Ref{Int}(fates_unset_int)  # SPITFIRE mode
const hlm_use_lu_harvest = Ref{Int}(fates_unset_int) # use HLM land-use harvest data?
const hlm_num_lu_harvest_cats = Ref{Int}(fates_unset_int) # number of HLM harvest categories
const hlm_use_luh           = Ref{Int}(fates_unset_int)  # use luh2 drivers?
const hlm_use_potentialveg  = Ref{Int}(fates_unset_int)  # potential vegetation only?
const hlm_num_luh2_states      = Ref{Int}(fates_unset_int) # number of LUH2 state types
const hlm_num_luh2_transitions = Ref{Int}(fates_unset_int) # number of LUH2 transition types
const hlm_sf_nofire_def               = Ref{Int}(fates_unset_int) # no-fire case definition
const hlm_sf_scalar_lightning_def     = Ref{Int}(fates_unset_int) # scalar-lightning case definition
const hlm_sf_successful_ignitions_def = Ref{Int}(fates_unset_int) # successful-ignition dataset case
const hlm_sf_anthro_ignitions_def     = Ref{Int}(fates_unset_int) # anthropogenic-ignition dataset case
const hlm_use_logging      = Ref{Int}(fates_unset_int)  # use logging module?
const hlm_use_planthydro   = Ref{Int}(fates_unset_int)  # use plant hydraulics?
const hlm_use_cohort_age_tracking = Ref{Int}(fates_unset_int) # cohort age tracking?
const hlm_use_tree_damage  = Ref{Int}(fates_unset_int)  # tree damage module?
const hlm_use_ed_st3       = Ref{Int}(fates_unset_int)  # static stand structure mode?
const hlm_use_ed_prescribed_phys = Ref{Int}(fates_unset_int) # prescribed physiology?
const hlm_use_inventory_init     = Ref{Int}(fates_unset_int) # init from inventory file?
const hlm_inventory_ctrl_file    = Ref{String}("")  # full path to inventory control file
const hlm_use_fixed_biogeog = Ref{Int}(fates_unset_int) # fixed biogeography mode?
const hlm_use_nocomp        = Ref{Int}(fates_unset_int) # no competition mode?
const hlm_use_sp            = Ref{Int}(fates_unset_int) # satellite phenology (LAI) mode?
const hlm_hist_level_dynam  = Ref{Int}(fates_unset_int) # history fields on dynamics (daily) step
const hlm_hist_level_hifrq  = Ref{Int}(fates_unset_int) # history fields on high-frequency step

# ===========================================================================
# Parameters dictated by FATES, required knowledge for the HLMs.
# Mostly used for dimensioning HLM array spaces.
# ===========================================================================

const fates_maxElementsPerPatch = Ref{Int}(fates_unset_int) # largest array size per patch (restart)
const fates_maxElementsPerSite  = Ref{Int}(fates_unset_int) # max individual items per grid cell
const fates_maxPatchesPerSite   = Ref{Int}(fates_unset_int) # patches FATES wants the HLM to allocate
const max_comp_per_site         = Ref{Int}(fates_unset_int) # max nutrient-acquisition competitors per site
const fates_dispersal_kernel_mode = Ref{Int}(fates_unset_int) # seed-dispersal kernel type

# Seed dispersal kernel modes (parameters)
const fates_dispersal_kernel_exponential = 1  # exponential dispersal kernel
const fates_dispersal_kernel_exppower    = 2  # exponential power (ExP) dispersal kernel
const fates_dispersal_kernel_logsech     = 3  # logistic-sech (LogS) dispersal kernel

# Seed dispersal cadences (parameters)
const fates_dispersal_cadence_none    = 0  # no dispersal (use seed rain only)
const fates_dispersal_cadence_daily   = 1  # disperse seeds daily
const fates_dispersal_cadence_monthly = 2  # disperse seeds monthly
const fates_dispersal_cadence_yearly  = 3  # disperse seeds yearly

# ===========================================================================
# History-output dimension mapping vectors.
# CLM/ALM have limited multi-dimensional history support, so FATES multiplexes
# multiple dimensions into one. These hold the dimension definitions and the
# maps from multiplexed indices back to their component dimensions.
# Fortran `allocatable` module arrays -> Julia Refs to Vectors (undef-0 default,
# allocated by allocate_bcin/allocate via the host once dimensions are known).
# r8 maps -> Vector{Float64}; integer maps -> Vector{Int}.
# ===========================================================================

const fates_hdim_levcoage        = Ref(Float64[])  # cohort age class lower bound dimension
const fates_hdim_pfmap_levcapf   = Ref(Int[])      # pfts into cohort age class x pft
const fates_hdim_camap_levcapf   = Ref(Int[])      # cohort age class into cohort age x pft
const fates_hdim_levsclass       = Ref(Float64[])  # plant size class lower bound dimension
const fates_hdim_pfmap_levscpf   = Ref(Int[])      # pfts into size-class x pft
const fates_hdim_scmap_levscpf   = Ref(Int[])      # size-class into size-class x pft
const fates_hdim_levage          = Ref(Float64[])  # patch age lower bound dimension
const fates_hdim_levheight       = Ref(Float64[])  # height lower bound dimension
const fates_hdim_levpft          = Ref(Int[])      # plant pft dimension
const fates_hdim_levlanduse      = Ref(Int[])      # land use label dimension
const fates_hdim_levfuel         = Ref(Int[])      # fire fuel size class dimension
const fates_hdim_levcwdsc        = Ref(Int[])      # cwd class dimension
const fates_hdim_levcan          = Ref(Int[])      # canopy-layer dimension
const fates_hdim_levleaf         = Ref(Float64[])  # leaf-layer dimension, integrated VAI [m2/m2]
const fates_hdim_levelem         = Ref(Int[])      # element dimension
const fates_hdim_canmap_levcnlf  = Ref(Int[])      # canopy-layer map into can-layer x leaf-layer
const fates_hdim_lfmap_levcnlf   = Ref(Int[])      # leaf-layer map into can-layer x leaf-layer
const fates_hdim_canmap_levcnlfpf = Ref(Int[])     # can-layer map into can-layer x pft x leaf-layer
const fates_hdim_lfmap_levcnlfpf  = Ref(Int[])     # leaf-layer map into can-layer x pft x leaf-layer
const fates_hdim_pftmap_levcnlfpf = Ref(Int[])     # pft map into can-layer x pft x leaf-layer
const fates_hdim_scmap_levscag   = Ref(Int[])      # size-class into size-class x patch age
const fates_hdim_agmap_levscag   = Ref(Int[])      # patch-age into size-class x patch age
const fates_hdim_scmap_levscagpft = Ref(Int[])     # size-class into size-class x patch age x pft
const fates_hdim_agmap_levscagpft = Ref(Int[])     # patch-age into size-class x patch age x pft
const fates_hdim_pftmap_levscagpft = Ref(Int[])    # pft into size-class x patch age x pft
const fates_hdim_agmap_levagepft  = Ref(Int[])     # patch-age into patch age x pft
const fates_hdim_pftmap_levagepft = Ref(Int[])     # pft into patch age x pft
const fates_hdim_agmap_levagefuel = Ref(Int[])     # patch-age into patch age x fsc
const fates_hdim_fscmap_levagefuel = Ref(Int[])    # fuel size-class into patch age x fsc
const fates_hdim_elmap_levelpft  = Ref(Int[])      # elements in element x pft
const fates_hdim_elmap_levelcwd  = Ref(Int[])      # elements in element x cwd
const fates_hdim_elmap_levelage  = Ref(Int[])      # elements in element x age
const fates_hdim_pftmap_levelpft = Ref(Int[])      # pfts in element x pft
const fates_hdim_cwdmap_levelcwd = Ref(Int[])      # cwds in element x cwd
const fates_hdim_agemap_levelage = Ref(Int[])      # ages in element x age
const fates_hdim_pftmap_levcdpf  = Ref(Int[])      # pfts into size x crowndamage x pft
const fates_hdim_cdmap_levcdpf   = Ref(Int[])      # crowndamage into size x crowndamage x pft
const fates_hdim_scmap_levcdpf   = Ref(Int[])      # size into size x crowndamage x pft
const fates_hdim_cdmap_levcdsc   = Ref(Int[])      # crowndamage into size x crowndamage
const fates_hdim_scmap_levcdsc   = Ref(Int[])      # size into size x crowndamage
const fates_hdim_levdamage       = Ref(Float64[])  # plant damage class lower bound dimension

# ===========================================================================
# Scalar timing variables. All sites on a machine are synchronous; the HLM
# controls time.
# ===========================================================================

const hlm_current_year   = Ref{Int}(fates_unset_int)   # current year
const hlm_current_month  = Ref{Int}(fates_unset_int)   # month of year
const hlm_current_day    = Ref{Int}(fates_unset_int)   # day of month
const hlm_current_tod    = Ref{Int}(fates_unset_int)   # time of day (seconds past 0Z)
const hlm_current_date   = Ref{Int}(fates_unset_int)   # current date
const hlm_reference_date = Ref{Int}(fates_unset_int)   # YYYYMMDD reference date
const hlm_model_day      = Ref{Float64}(fates_unset_r8) # elapsed days between current and ref date
const hlm_day_of_year    = Ref{Int}(fates_unset_int)   # integer day of the year
const hlm_days_per_year  = Ref{Int}(fates_unset_int)   # days per year (HLM may include leap)
const hlm_freq_day       = Ref{Float64}(fates_unset_r8) # fraction of year for a daily step

# ===========================================================================
# Constant parameters dictated by the fates parameter file.
# ===========================================================================

const numpft      = Ref{Int}(fates_unset_int)  # total number of PFTs in the simulation
const nlevsclass  = Ref{Int}(fates_unset_int)  # number of cohort size class bins (history)
const nlevage     = Ref{Int}(fates_unset_int)  # number of patch age bins (history)
const nlevheight  = Ref{Int}(fates_unset_int)  # number of height bins (history)
const nlevcoage   = Ref{Int}(fates_unset_int)  # number of cohort age bins (history)
const nleafage    = Ref{Int}(fates_unset_int)  # number of leaf age classes
const nlevdamage  = Ref{Int}(fates_unset_int)  # number of damage classes

# ===========================================================================
# Structured boundary conditions (SITE/PATCH scale).
# Naming conventions:  _si site (scalar), _pa patch, _rb radiation band,
#                      _sl soil layer, _sisl site x soil layer.
# ===========================================================================

"""
    bc_in_type

Boundary conditions passed IN from the host land model to FATES. SoA layout;
allocatable arrays default to undef-0 and are sized by [`allocate_bcin!`](@ref).
Field names match the Fortran `bc_in_type` exactly.
"""
Base.@kwdef mutable struct bc_in_type
    # The actual number of FATES' ED patches
    npatches::Int = fates_unset_int

    # Soil layer structure
    nlevsoil::Int   = fates_unset_int  # number of soil layers in this column
    nlevdecomp::Int = fates_unset_int  # number of biogeochemically active soil layers
    zi_sisl::Vector{Float64} = Float64[]  # interface level below a "z" level (m); incl. zero index for surface
    dz_sisl::Vector{Float64} = Float64[]  # layer thickness (m)
    z_sisl::Vector{Float64}  = Float64[]  # layer depth (m)

    # Decomposition layer structure
    dz_decomp_sisl::Vector{Float64} = Float64[]  # decomp layer thickness (m); == dz_sisl unless single layer
    decomp_id::Vector{Int} = Int[]               # soil layer -> decomposition layer index map

    # Decomposition fractions
    w_scalar_sisl::Vector{Float64} = Float64[]  # moisture limitation on decomposition (fraction)
    t_scalar_sisl::Vector{Float64} = Float64[]  # temperature limitation on decomposition (fraction)

    # Fire model
    lightning24::Vector{Float64}   = Float64[]  # 24-hour lightning or ignitions [#/km2/day]
    pop_density::Vector{Float64}   = Float64[]  # population density [#/km2]
    precip24_pa::Vector{Float64}   = Float64[]  # avg precipitation over last 24h [mm/s]
    relhumid24_pa::Vector{Float64} = Float64[]  # avg relative humidity over last 24h [-]
    wind24_pa::Vector{Float64}     = Float64[]  # patch 24-hour running mean wind (m/s)

    # Radiation variables for sun/shade fractions (patch, radiation-band)
    solad_parb::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # downwelling direct beam radiation [W/m2]
    solai_parb::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # downwelling diffuse radiation [W/m2]

    # Nutrient input fluxes (integrated over the day) — Fortran pointers
    plant_nh4_uptake_flux::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # ammonium uptake per competitor [gN/m2/day]
    plant_no3_uptake_flux::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # nitrate uptake per competitor [gN/m2/day]
    plant_p_uptake_flux::Matrix{Float64}   = Matrix{Float64}(undef, 0, 0)  # phosphorus input per competitor [gP/m2/day]

    # Photosynthesis variables
    filter_photo_pa::Vector{Int} = Int[]  # patch photosynthesis filter flag (1/2/3)
    forc_pbot::Float64 = fates_unset_r8    # atmospheric pressure (Pa)
    dayl_factor_pa::Vector{Float64} = Float64[]  # daylength scaling factor (0-1)
    esat_tv_pa::Vector{Float64} = Float64[]      # saturation vapor pressure at t_veg (Pa)
    eair_pa::Vector{Float64}    = Float64[]      # vapor pressure of canopy air (Pa)
    oair_pa::Vector{Float64}    = Float64[]      # atmospheric O2 partial pressure (Pa)
    cair_pa::Vector{Float64}    = Float64[]      # atmospheric CO2 partial pressure (Pa)
    rb_pa::Vector{Float64}      = Float64[]      # boundary layer resistance (s/m)
    t_veg_pa::Vector{Float64}   = Float64[]      # vegetation temperature (K)
    tgcm_pa::Vector{Float64}    = Float64[]      # air temperature at agcm reference height (K)
    t_soisno_sl::Vector{Float64} = Float64[]     # soil temperature (K)

    # Canopy radiation boundaries
    filter_vegzen_pa::Vector{Bool} = Bool[]  # filter for patches with positive zenith angle (daylight)
    coszen_pa::Vector{Float64} = Float64[]   # cosine of the zenith angle (0-1), by patch
    fcansno_pa::Vector{Float64} = Float64[]  # fraction of canopy covered in snow
    albgr_dir_rb::Vector{Float64} = Float64[]  # ground albedo, direct, by site broadband (0-1)
    albgr_dif_rb::Vector{Float64} = Float64[]  # ground albedo, diffuse, by site broadband (0-1)

    # LitterFlux boundaries
    max_rooting_depth_index_col::Int = fates_unset_int  # deepest soil level where roots may be

    # BGC accounting
    tot_het_resp::Float64 = fates_unset_r8  # total heterotrophic respiration (gC/m2/s)
    tot_somc::Float64     = fates_unset_r8  # total soil organic matter carbon (gC/m2)
    tot_litc::Float64     = fates_unset_r8  # total litter carbon tracked in the HLM (gC/m2)

    # Canopy structure
    snow_depth_si::Float64   = fates_unset_r8  # depth of snow in snowy areas of site (m)
    frac_sno_eff_si::Float64 = fates_unset_r8  # fraction of ground covered by snow (0-1)

    # Hydrology variables for BTRAN
    smp_sl::Vector{Float64}         = Float64[]  # soil suction potential of layers, negative [mm]
    salinity_sl::Vector{Float64}    = Float64[]  # soil salinity of layers [ppt]
    eff_porosity_sl::Vector{Float64} = Float64[] # effective porosity = porosity - vol_ic [-]
    watsat_sl::Vector{Float64}      = Float64[]  # volumetric soil water at saturation (porosity)
    tempk_sl::Vector{Float64}       = Float64[]  # temperature of soil layers [K]
    h2o_liqvol_sl::Vector{Float64}  = Float64[]  # liquid volume in soil layer (m3/m3)
    filter_btran::Bool = false                   # site level filter for uptake response functions

    # Plant-Hydro (allocated on rhizosphere levels)
    qflx_transp_pa::Vector{Float64} = Float64[]  # transpiration flux from HLM canopy solver [mm H2O/s, + into root]
    swrad_net_pa::Vector{Float64}   = Float64[]  # net absorbed shortwave radiation (W/m2)
    lwrad_net_pa::Vector{Float64}   = Float64[]  # net absorbed longwave radiation (W/m2)
    watsat_sisl::Vector{Float64}    = Float64[]  # volumetric soil water at saturation (porosity)
    watres_sisl::Vector{Float64}    = Float64[]  # volumetric residual soil water
    sucsat_sisl::Vector{Float64}    = Float64[]  # minimum soil suction (mm)
    bsw_sisl::Vector{Float64}       = Float64[]  # Clapp and Hornberger "b"
    hksat_sisl::Vector{Float64}     = Float64[]  # hydraulic conductivity at saturation (mm H2O/s)
    h2o_liq_sisl::Vector{Float64}   = Float64[]  # liquid water mass in each layer (kg/m2)
    smpmin_si::Float64 = fates_unset_r8          # restriction for min of soil potential (mm)

    # Land use
    hlm_harvest_rates::Vector{Float64}      = Float64[]  # annual harvest rate per category for a site
    hlm_harvest_catnames::Vector{String}    = String[]   # names of hlm_harvest d1
    hlm_luh_states::Vector{Float64}         = Float64[]  # LUH2 states
    hlm_luh_state_names::Vector{String}     = String[]   # names of LUH2 states
    hlm_luh_transitions::Vector{Float64}    = Float64[]  # LUH2 transitions
    hlm_luh_transition_names::Vector{String} = String[]  # names of LUH2 transitions
    hlm_harvest_units::Int = fates_unset_int             # units of harvest rates (area vs carbon)
    site_area::Float64 = fates_unset_r8                  # actual area of current site [m2]

    # Fixed biogeography mode
    pft_areafrac::Vector{Float64} = Float64[]            # fractional area occupied by each PFT
    pft_areafrac_lu::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # fractional area per PFT per land use
    baregroundfrac::Float64 = fates_unset_r8             # fractional area held as bare-ground

    # Satellite Phenology (SP) input variables (each patch has one PFT)
    hlm_sp_tlai::Vector{Float64} = Float64[]  # interpolated daily total LAI per patch/pft
    hlm_sp_tsai::Vector{Float64} = Float64[]  # interpolated daily total SAI per patch/pft
    hlm_sp_htop::Vector{Float64} = Float64[]  # interpolated daily canopy height per patch/pft
end


"""
    bc_out_type

Boundary conditions passed OUT from FATES to the host land model. SoA layout;
allocatable arrays default to undef-0 and are sized by [`allocate_bcout!`](@ref).
Field names match the Fortran `bc_out_type` exactly.
"""
Base.@kwdef mutable struct bc_out_type
    # Sun/shade canopy
    fsun_pa::Vector{Float64}   = Float64[]  # sunlit fraction of the canopy for this patch [0-1]
    laisun_pa::Vector{Float64} = Float64[]  # sunlit canopy LAI
    laisha_pa::Vector{Float64} = Float64[]  # shaded canopy LAI

    active_suction_sl::Vector{Bool} = Bool[]  # whether a soil layer can have water uptake by plants
    rootr_pasl::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # effective fraction of roots per soil layer
    btran_pa::Vector{Float64}  = Float64[]  # integrated transpiration wetness factor (0-1), diagnostic
    rssun_pa::Vector{Float64}  = Float64[]  # sunlit canopy resistance [s/m]
    rssha_pa::Vector{Float64}  = Float64[]  # shaded canopy resistance [s/m]

    # Canopy radiation boundaries (patch, radiation-band)
    albd_parb::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # surface albedo (direct)
    albi_parb::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # surface albedo (diffuse)
    fabd_parb::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # flux absorbed by canopy per unit direct flux
    fabi_parb::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # flux absorbed by canopy per unit diffuse flux
    ftdd_parb::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # down direct flux below canopy per unit direct
    ftid_parb::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # down diffuse flux below canopy per unit direct
    ftii_parb::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # down diffuse flux below canopy per unit diffuse

    # Mass fluxes to BGC from fragmentation of litter into decomposing pools [g/m3/s]
    litt_flux_cel_c_si::Vector{Float64} = Float64[]  # cellulose carbon litter
    litt_flux_lig_c_si::Vector{Float64} = Float64[]  # lignin carbon litter
    litt_flux_lab_c_si::Vector{Float64} = Float64[]  # labile carbon litter
    litt_flux_cel_n_si::Vector{Float64} = Float64[]  # cellulose nitrogen litter
    litt_flux_lig_n_si::Vector{Float64} = Float64[]  # lignin nitrogen litter
    litt_flux_lab_n_si::Vector{Float64} = Float64[]  # labile nitrogen litter
    litt_flux_cel_p_si::Vector{Float64} = Float64[]  # cellulose phosphorus litter
    litt_flux_lig_p_si::Vector{Float64} = Float64[]  # lignin phosphorus litter
    litt_flux_lab_p_si::Vector{Float64} = Float64[]  # labile phosphorus litter

    # MIMICS boundary conditions
    litt_flux_ligc_per_n::Float64 = fates_unset_r8  # lignin carbon per total nitrogen in frag flux [g/g]

    # Nutrient competition boundary conditions (Fortran pointers)
    num_plant_comps::Int = fates_unset_int  # number of unique competitors
    source_nh4::Vector{Float64} = Float64[]  # FATES source of ammonium to mineralized N pool [gN/m3]
    source_p::Vector{Float64}   = Float64[]  # FATES source of phosphorus to mineralized P pool [gP/m3]
    veg_rootc::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # fine-root carbon per competitor [gC/m3] (maxcohort x nlevdecomp)
    decompmicc::Vector{Float64} = Float64[]  # microbial decomposer biomass [gC/m3] (numpft x nlevdecomp_full)
    ft_index::Vector{Int} = Int[]            # functional type index of each competitor (maxcohort)
    cn_scalar::Vector{Float64} = Float64[]   # C:N scaling factor for root N uptake kinetics
    cp_scalar::Vector{Float64} = Float64[]   # C:P scaling factor for root P uptake kinetics

    # CH4 boundary conditions (Fortran pointers)
    annavg_agnpp_pa::Vector{Float64} = Float64[]    # annual average patch npp above ground (gC/m2/s)
    annavg_bgnpp_pa::Vector{Float64} = Float64[]    # annual average patch npp below ground (gC/m2/s)
    annsum_npp_pa::Vector{Float64}   = Float64[]    # annual sum patch npp (gC/m2/yr)
    frootc_pa::Vector{Float64}       = Float64[]    # carbon in fine roots (gC/m2)
    root_resp::Vector{Float64}       = Float64[]    # root respiration (gC/m2/s)
    rootfr_pa::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # rooting fraction with depth
    woody_frac_aere_pa::Vector{Float64} = Float64[]  # woody plant fraction (by crown area) for aerenchyma porosity
    ema_npp::Float64 = fates_unset_r8                # site-level NPP smoothed over time

    # Canopy structure
    elai_pa::Vector{Float64} = Float64[]  # exposed leaf area index
    esai_pa::Vector{Float64} = Float64[]  # exposed stem area index
    tlai_pa::Vector{Float64} = Float64[]  # total leaf area index
    tsai_pa::Vector{Float64} = Float64[]  # total stem area index
    htop_pa::Vector{Float64} = Float64[]  # top of the canopy [m]
    hbot_pa::Vector{Float64} = Float64[]  # bottom of canopy [m]
    z0m_pa::Vector{Float64}    = Float64[]  # roughness length [m]
    displa_pa::Vector{Float64} = Float64[]  # displacement height [m]
    dleaf_pa::Vector{Float64}  = Float64[]  # leaf characteristic dimension/width/diameter [m]
    canopy_fraction_pa::Vector{Float64} = Float64[]  # area fraction of each patch in the site [0-1]
    frac_veg_nosno_alb_pa::Vector{Float64} = Float64[]  # binary: any vegetation in patch exposed [0,1]
    nocomp_pft_label_pa::Vector{Int} = Int[]  # in nocomp/SP mode, the PFT identity of each patch

    # FATES hydraulics
    plant_stored_h2o_si::Float64 = fates_unset_r8  # stored water in live+dead vegetation (kg/m2 H2O)
    qflx_soil2root_sisl::Vector{Float64} = Float64[]  # water flux from soil into root [mm H2O/s, + into root]
    qflx_ro_sisl::Vector{Float64} = Float64[]         # runoff from root->soil super-saturation [mm H2O/s]

    # FATES LULCC
    hrv_deadstemc_to_prod10c::Float64  = fates_unset_r8  # harvested C flux to 10-yr wood product pool [gC/m2/s]
    hrv_deadstemc_to_prod100c::Float64 = fates_unset_r8  # harvested C flux to 100-yr wood product pool [gC/m2/s]
    gpp_site::Float64 = fates_unset_r8  # site level GPP, for NBP diagnosis in HLM [gC/m2/s]
    ar_site::Float64  = fates_unset_r8  # site level autotrophic resp, for NBP diagnosis in HLM [gC/m2/s]
end


"""
    bc_pconst_type

Parameter constants for the host<->FATES boundary. Specified once and never
modified. The ECA-hypothesis nutrient-competition vectors are dimensioned by
PFT. Field names match the Fortran `bc_pconst_type` exactly.
"""
Base.@kwdef mutable struct bc_pconst_type
    max_plant_comps::Int = fates_unset_int

    vmax_nh4::Vector{Float64} = Float64[]
    vmax_no3::Vector{Float64} = Float64[]
    vmax_p::Vector{Float64}   = Float64[]

    eca_km_nh4::Vector{Float64} = Float64[]
    eca_km_no3::Vector{Float64} = Float64[]
    eca_km_p::Vector{Float64}   = Float64[]

    eca_km_ptase::Vector{Float64}     = Float64[]
    eca_vmax_ptase::Vector{Float64}   = Float64[]
    eca_alpha_ptase::Vector{Float64}  = Float64[]
    eca_lambda_ptase::Vector{Float64} = Float64[]
    eca_plant_escalar::Float64 = fates_unset_r8

    j_uptake::Vector{Int} = Int[]  # decomposition-layer -> uptake-layer map
end


# ===========================================================================
# Allocate / init routines.
# The Fortran allocate routines for bc_in/bc_out live in the *interface* module
# (FatesInterfaceMod), but the sizing is fully determined by the dimensions
# carried here (npatches, nlevsoil, nlevdecomp, number of radiation bands,
# competitors, etc.). We provide allocators sized from those dimensions so the
# types are usable standalone. Allocated arrays are zero-initialized (Fortran
# allocate leaves them undefined; FATES zeroes them in the subsequent set
# routines). Field names and shapes follow the Fortran declarations.
# ===========================================================================

"""
    allocate_bcin!(bc_in; npatches, nlevsoil, nlevdecomp, num_rad_bands=2,
                   max_comp=0, num_harvest_cats=0, num_luh2_states=0,
                   num_luh2_transitions=0, num_pft=0, num_lu=0, num_sp=0)

Allocate the host->FATES boundary-condition arrays of `bc_in` from a dimension
set. Patch arrays get length `npatches`, soil-layer arrays `nlevsoil`,
site x soil-layer arrays `nlevsoil`, and radiation arrays (npatches x
num_rad_bands). Scalar dimensions are stored on the struct. Returns `bc_in`.
"""
function allocate_bcin!(bc_in::bc_in_type; npatches::Integer, nlevsoil::Integer,
                        nlevdecomp::Integer, num_rad_bands::Integer=2,
                        max_comp::Integer=0, num_harvest_cats::Integer=0,
                        num_luh2_states::Integer=0, num_luh2_transitions::Integer=0,
                        num_pft::Integer=0, num_lu::Integer=0, num_sp::Integer=npatches)
    bc_in.npatches   = npatches
    bc_in.nlevsoil   = nlevsoil
    bc_in.nlevdecomp = nlevdecomp

    # Soil layer structure. zi_sisl carries a zero index for the surface -> length nlevsoil+1.
    bc_in.zi_sisl = zeros(Float64, nlevsoil + 1)
    bc_in.dz_sisl = zeros(Float64, nlevsoil)
    bc_in.z_sisl  = zeros(Float64, nlevsoil)

    # Decomposition layer structure
    bc_in.dz_decomp_sisl = zeros(Float64, nlevdecomp)
    bc_in.decomp_id      = zeros(Int, nlevsoil)
    bc_in.w_scalar_sisl  = zeros(Float64, nlevdecomp)
    bc_in.t_scalar_sisl  = zeros(Float64, nlevdecomp)

    # Fire model (patch + site)
    bc_in.lightning24   = zeros(Float64, npatches)
    bc_in.pop_density   = zeros(Float64, npatches)
    bc_in.precip24_pa   = zeros(Float64, npatches)
    bc_in.relhumid24_pa = zeros(Float64, npatches)
    bc_in.wind24_pa     = zeros(Float64, npatches)

    # Radiation (patch x band)
    bc_in.solad_parb = zeros(Float64, npatches, num_rad_bands)
    bc_in.solai_parb = zeros(Float64, npatches, num_rad_bands)

    # Nutrient input fluxes (competitor x 1; see Fortran note)
    bc_in.plant_nh4_uptake_flux = zeros(Float64, max_comp, 1)
    bc_in.plant_no3_uptake_flux = zeros(Float64, max_comp, 1)
    bc_in.plant_p_uptake_flux   = zeros(Float64, max_comp, 1)

    # Photosynthesis (patch)
    bc_in.filter_photo_pa = zeros(Int, npatches)
    bc_in.dayl_factor_pa  = zeros(Float64, npatches)
    bc_in.esat_tv_pa      = zeros(Float64, npatches)
    bc_in.eair_pa         = zeros(Float64, npatches)
    bc_in.oair_pa         = zeros(Float64, npatches)
    bc_in.cair_pa         = zeros(Float64, npatches)
    bc_in.rb_pa           = zeros(Float64, npatches)
    bc_in.t_veg_pa        = zeros(Float64, npatches)
    bc_in.tgcm_pa         = zeros(Float64, npatches)
    bc_in.t_soisno_sl     = zeros(Float64, nlevsoil)

    # Canopy radiation
    bc_in.filter_vegzen_pa = falses(npatches)
    bc_in.coszen_pa    = zeros(Float64, npatches)
    bc_in.fcansno_pa   = zeros(Float64, npatches)
    bc_in.albgr_dir_rb = zeros(Float64, num_rad_bands)
    bc_in.albgr_dif_rb = zeros(Float64, num_rad_bands)

    # Hydrology for BTRAN (soil layer)
    bc_in.smp_sl          = zeros(Float64, nlevsoil)
    bc_in.salinity_sl     = zeros(Float64, nlevsoil)
    bc_in.eff_porosity_sl = zeros(Float64, nlevsoil)
    bc_in.watsat_sl       = zeros(Float64, nlevsoil)
    bc_in.tempk_sl        = zeros(Float64, nlevsoil)
    bc_in.h2o_liqvol_sl   = zeros(Float64, nlevsoil)

    # Plant-Hydro (patch + site x soil layer)
    bc_in.qflx_transp_pa = zeros(Float64, npatches)
    bc_in.swrad_net_pa   = zeros(Float64, npatches)
    bc_in.lwrad_net_pa   = zeros(Float64, npatches)
    bc_in.watsat_sisl    = zeros(Float64, nlevsoil)
    bc_in.watres_sisl    = zeros(Float64, nlevsoil)
    bc_in.sucsat_sisl    = zeros(Float64, nlevsoil)
    bc_in.bsw_sisl       = zeros(Float64, nlevsoil)
    bc_in.hksat_sisl     = zeros(Float64, nlevsoil)
    bc_in.h2o_liq_sisl   = zeros(Float64, nlevsoil)

    # Land use
    bc_in.hlm_harvest_rates        = zeros(Float64, num_harvest_cats)
    bc_in.hlm_harvest_catnames     = fill("", num_harvest_cats)
    bc_in.hlm_luh_states           = zeros(Float64, num_luh2_states)
    bc_in.hlm_luh_state_names      = fill("", num_luh2_states)
    bc_in.hlm_luh_transitions      = zeros(Float64, num_luh2_transitions)
    bc_in.hlm_luh_transition_names = fill("", num_luh2_transitions)

    # Fixed biogeography
    bc_in.pft_areafrac    = zeros(Float64, num_pft)
    bc_in.pft_areafrac_lu = zeros(Float64, num_pft, num_lu)

    # Satellite Phenology
    bc_in.hlm_sp_tlai = zeros(Float64, num_sp)
    bc_in.hlm_sp_tsai = zeros(Float64, num_sp)
    bc_in.hlm_sp_htop = zeros(Float64, num_sp)

    return bc_in
end

"""
    allocate_bcout!(bc_out; npatches, nlevsoil, nlevdecomp, num_rad_bands=2,
                    max_comp=0, num_pft=0)

Allocate the FATES->host boundary-condition arrays of `bc_out` from a dimension
set. Patch arrays get length `npatches`; site-level litter/source arrays
`nlevdecomp`; radiation arrays (npatches x num_rad_bands); rootr/rootfr
(npatches x nlevsoil); veg_rootc (max_comp x nlevdecomp). Returns `bc_out`.
"""
function allocate_bcout!(bc_out::bc_out_type; npatches::Integer, nlevsoil::Integer,
                         nlevdecomp::Integer, num_rad_bands::Integer=2,
                         max_comp::Integer=0, num_pft::Integer=0)
    # Sun/shade canopy (patch)
    bc_out.fsun_pa   = zeros(Float64, npatches)
    bc_out.laisun_pa = zeros(Float64, npatches)
    bc_out.laisha_pa = zeros(Float64, npatches)

    bc_out.active_suction_sl = falses(nlevsoil)
    bc_out.rootr_pasl = zeros(Float64, npatches, nlevsoil)
    bc_out.btran_pa   = zeros(Float64, npatches)
    bc_out.rssun_pa   = zeros(Float64, npatches)
    bc_out.rssha_pa   = zeros(Float64, npatches)

    # Canopy radiation (patch x band)
    bc_out.albd_parb = zeros(Float64, npatches, num_rad_bands)
    bc_out.albi_parb = zeros(Float64, npatches, num_rad_bands)
    bc_out.fabd_parb = zeros(Float64, npatches, num_rad_bands)
    bc_out.fabi_parb = zeros(Float64, npatches, num_rad_bands)
    bc_out.ftdd_parb = zeros(Float64, npatches, num_rad_bands)
    bc_out.ftid_parb = zeros(Float64, npatches, num_rad_bands)
    bc_out.ftii_parb = zeros(Float64, npatches, num_rad_bands)

    # Litter mass fluxes (site x decomp layer)
    bc_out.litt_flux_cel_c_si = zeros(Float64, nlevdecomp)
    bc_out.litt_flux_lig_c_si = zeros(Float64, nlevdecomp)
    bc_out.litt_flux_lab_c_si = zeros(Float64, nlevdecomp)
    bc_out.litt_flux_cel_n_si = zeros(Float64, nlevdecomp)
    bc_out.litt_flux_lig_n_si = zeros(Float64, nlevdecomp)
    bc_out.litt_flux_lab_n_si = zeros(Float64, nlevdecomp)
    bc_out.litt_flux_cel_p_si = zeros(Float64, nlevdecomp)
    bc_out.litt_flux_lig_p_si = zeros(Float64, nlevdecomp)
    bc_out.litt_flux_lab_p_si = zeros(Float64, nlevdecomp)

    # Nutrient competition
    bc_out.source_nh4 = zeros(Float64, nlevdecomp)
    bc_out.source_p   = zeros(Float64, nlevdecomp)
    bc_out.veg_rootc  = zeros(Float64, max_comp, nlevdecomp)
    bc_out.decompmicc = zeros(Float64, nlevdecomp)
    bc_out.ft_index   = zeros(Int, max_comp)
    bc_out.cn_scalar  = zeros(Float64, max_comp)
    bc_out.cp_scalar  = zeros(Float64, max_comp)

    # CH4 boundary conditions (patch + patch x soil)
    bc_out.annavg_agnpp_pa    = zeros(Float64, npatches)
    bc_out.annavg_bgnpp_pa    = zeros(Float64, npatches)
    bc_out.annsum_npp_pa      = zeros(Float64, npatches)
    bc_out.frootc_pa          = zeros(Float64, npatches)
    bc_out.root_resp          = zeros(Float64, npatches)
    bc_out.rootfr_pa          = zeros(Float64, npatches, nlevsoil)
    bc_out.woody_frac_aere_pa = zeros(Float64, npatches)

    # Canopy structure (patch)
    bc_out.elai_pa = zeros(Float64, npatches)
    bc_out.esai_pa = zeros(Float64, npatches)
    bc_out.tlai_pa = zeros(Float64, npatches)
    bc_out.tsai_pa = zeros(Float64, npatches)
    bc_out.htop_pa = zeros(Float64, npatches)
    bc_out.hbot_pa = zeros(Float64, npatches)
    bc_out.z0m_pa    = zeros(Float64, npatches)
    bc_out.displa_pa = zeros(Float64, npatches)
    bc_out.dleaf_pa  = zeros(Float64, npatches)
    bc_out.canopy_fraction_pa    = zeros(Float64, npatches)
    bc_out.frac_veg_nosno_alb_pa = zeros(Float64, npatches)
    bc_out.nocomp_pft_label_pa   = zeros(Int, npatches)

    # FATES hydraulics (site x soil layer)
    bc_out.qflx_soil2root_sisl = zeros(Float64, nlevsoil)
    bc_out.qflx_ro_sisl        = zeros(Float64, nlevsoil)

    return bc_out
end

"""
    allocate_bcpconst!(bc_pconst; num_pft, max_plant_comps=0, nlevdecomp=0)

Allocate the PFT-dimensioned ECA nutrient-competition parameter vectors of
`bc_pconst`, and the decomposition->uptake layer map `j_uptake` (nlevdecomp).
Returns `bc_pconst`.
"""
function allocate_bcpconst!(bc_pconst::bc_pconst_type; num_pft::Integer,
                            max_plant_comps::Integer=0, nlevdecomp::Integer=0)
    bc_pconst.max_plant_comps = max_plant_comps

    bc_pconst.vmax_nh4 = zeros(Float64, num_pft)
    bc_pconst.vmax_no3 = zeros(Float64, num_pft)
    bc_pconst.vmax_p   = zeros(Float64, num_pft)

    bc_pconst.eca_km_nh4 = zeros(Float64, num_pft)
    bc_pconst.eca_km_no3 = zeros(Float64, num_pft)
    bc_pconst.eca_km_p   = zeros(Float64, num_pft)

    bc_pconst.eca_km_ptase     = zeros(Float64, num_pft)
    bc_pconst.eca_vmax_ptase   = zeros(Float64, num_pft)
    bc_pconst.eca_alpha_ptase  = zeros(Float64, num_pft)
    bc_pconst.eca_lambda_ptase = zeros(Float64, num_pft)

    bc_pconst.j_uptake = zeros(Int, nlevdecomp)

    return bc_pconst
end
