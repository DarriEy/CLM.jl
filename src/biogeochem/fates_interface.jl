# ============================================================================
# FATES Interface — CLM-FATES boundary condition types and integration stubs
#
# Ported from:
#   src/utils/clmfates_interfaceMod.F90  (~3,896 lines)
#   src/fates/main/FatesInterfaceTypesMod.F90  (~821 lines)
#
# This file defines all CLM<->FATES data exchange types and provides
# documented stub functions at every integration point. The stubs are
# no-ops that make it easy to incrementally port FATES dynamics later.
#
# FATES (Functionally Assembled Terrestrial Ecosystem Simulator) is a
# cohort-based vegetation demography model that replaces CLM's big-leaf
# approach with explicit tracking of plant cohorts within patches of
# different disturbance ages. Key FATES subsystems include:
#   - ED (Ecosystem Demography): growth, mortality, recruitment, disturbance
#   - SPITFIRE: fire behavior and effects
#   - Plant hydraulics: explicit plant water transport
#   - Canopy radiation: multi-layer canopy radiative transfer
#   - Carbon/nutrient allocation (PARTEH)
#
# The Fortran FATES source tree (~10,000+ lines) is organized as:
#   fates/main/         — core types, interface, parameters
#   fates/biogeophys/   — radiation, hydraulics, photosynthesis wrappers
#   fates/biogeochem/   — carbon/nutrient cycling, litter, soil BGC coupling
#   fates/fire/         — SPITFIRE fire model
#   fates/radiation/    — two-stream canopy RT
#   fates/parteh/       — plant allocation and reactive transport
#
# Architecture notes:
#   - In CLM-Fortran, each site (column) has its own bc_in/bc_out pair
#   - FATES patches are sub-column entities with variable number per site
#   - The HLM (Host Land Model) controls time, soil physics, and atmosphere
#   - FATES controls vegetation structure, demographics, and fire
# ============================================================================

# ============================================================================
#                     FATES GLOBAL CONFIGURATION
# ============================================================================

"""
    FATESConfig

Global configuration parameters for FATES, set once during initialization.
Corresponds to the `hlm_*` global variables in `FatesInterfaceTypesMod.F90`.
"""
Base.@kwdef mutable struct FATESConfig
    # --- HLM parameters passed to FATES ---
    hlm_numSWb::Int = 2                    # Number of shortwave broadbands (VIS/NIR)
    hlm_ivis::Int = 1                      # Array index for visible band
    hlm_inir::Int = 2                      # Array index for NIR band
    hlm_maxlevsoil::Int = 20               # Max number of soil layers
    hlm_stepsize::Float64 = 1800.0         # HLM timestep [s]
    hlm_name::String = "CLM"               # Host land model name
    hlm_decomp::String = "CENTURY"         # Soil decomposition scheme (CENTURY/MIMICS/CTC/NONE)
    hlm_nu_com::String = "RD"              # Nutrient competition scheme (RD/ECA/NONE)

    # --- FATES feature flags (0=false, 1=true to match Fortran convention) ---
    hlm_use_planthydro::Bool = false       # Plant hydraulics (experimental)
    hlm_use_ed_st3::Bool = false           # Static stand structure mode
    hlm_use_ed_prescribed_phys::Bool = false # Prescribed physiology
    hlm_use_inventory_init::Bool = false   # Inventory-based initialization
    hlm_use_fixed_biogeog::Bool = false    # Fixed biogeography mode
    hlm_use_nocomp::Bool = false           # No-competition mode
    hlm_use_sp::Bool = false               # Satellite phenology mode
    hlm_use_logging::Bool = false          # Logging module
    hlm_use_lu_harvest::Bool = false       # Land-use harvest
    hlm_use_luh::Bool = false              # LUH2 land-use drivers
    hlm_use_potentialveg::Bool = false     # Potential vegetation only
    hlm_use_tree_damage::Bool = false      # Tree damage module
    hlm_use_cohort_age_tracking::Bool = false  # Cohort age tracking
    hlm_use_ch4::Bool = false              # Methane model coupling
    hlm_use_vertsoilc::Bool = true         # Vertically discretized soil C

    # --- Fire configuration ---
    hlm_spitfire_mode::Int = 0             # SPITFIRE mode (0=no fire)

    # --- PARTEH allocation mode ---
    hlm_parteh_mode::Int = 1               # 1=carbon-only, 2=CNP

    # --- Seed dispersal ---
    hlm_seeddisp_cadence::Int = 0          # 0=none, 1=daily, 2=monthly, 3=yearly

    # --- Nitrogen/phosphorus ---
    hlm_nitrogen_spec::Int = 0             # 0=none, 1=NH4 only, 2=NH4+NO3
    hlm_phosphorus_spec::Int = 0           # 0=none, 1=P is on

    # --- History output levels ---
    hlm_hist_level_dynam::Int = 1          # Dynamics-step history level
    hlm_hist_level_hifrq::Int = 1          # High-frequency history level

    # --- FATES-determined dimensions ---
    numpft::Int = 0                        # Number of FATES PFTs
    nlevsclass::Int = 0                    # Number of size classes
    nlevage::Int = 0                       # Number of patch age bins
    nlevheight::Int = 0                    # Number of height bins
    nlevcoage::Int = 0                     # Number of cohort age bins
    nleafage::Int = 0                      # Number of leaf age classes
    nlevdamage::Int = 0                    # Number of damage classes
    fates_maxPatchesPerSite::Int = 0       # Max patches per site
    fates_maxElementsPerPatch::Int = 0     # Max elements per patch (for restart sizing)
    fates_maxElementsPerSite::Int = 0      # Max elements per site
end

# ============================================================================
#                     BOUNDARY CONDITION TYPES
# ============================================================================

"""
    FATESBoundaryCondIn{FT}

Input boundary conditions passed from CLM to FATES for a single site (column).
Corresponds to `bc_in_type` in `FatesInterfaceTypesMod.F90`.

Fields are organized by subsystem. Patch-dimensioned arrays use index `[ipatch]`,
soil-layer arrays use `[isoil]`, and radiation-band arrays use `[iband]`.
"""
Base.@kwdef mutable struct FATESBoundaryCondIn{FT<:Real}

    # --- Patch count ---
    npatches::Int = 0                              # Actual number of FATES patches at this site

    # --- Soil layer structure ---
    nlevsoil::Int = 0                              # Number of soil layers
    nlevdecomp::Int = 0                            # Number of biogeochemically active soil layers
    zi_sisl::Vector{FT} = Float64[]                # [m] Interface level below z-level (0:nlevsoil)
    dz_sisl::Vector{FT} = Float64[]                # [m] Layer thickness (1:nlevsoil)
    z_sisl::Vector{FT} = Float64[]                 # [m] Layer depth / midpoint (1:nlevsoil)
    dz_decomp_sisl::Vector{FT} = Float64[]         # [m] Decomposition layer thickness
    decomp_id::Vector{Int} = Int[]                 # Decomp layer index each soil layer maps to

    # --- Decomposition scalars ---
    w_scalar_sisl::Vector{FT} = Float64[]          # [-] Moisture limitation on decomposition (1:nlevsoil)
    t_scalar_sisl::Vector{FT} = Float64[]          # [-] Temperature limitation on decomposition (1:nlevsoil)

    # --- Fire inputs ---
    lightning24::Vector{FT} = Float64[]            # [#/km2/day] 24-hour lightning rate (1:npatches)
    pop_density::Vector{FT} = Float64[]            # [#/km2] Population density (1:npatches)
    precip24_pa::Vector{FT} = Float64[]            # [mm/s] 24-hour mean precipitation (1:npatches)
    relhumid24_pa::Vector{FT} = Float64[]          # [-] 24-hour mean relative humidity (1:npatches)
    wind24_pa::Vector{FT} = Float64[]              # [m/s] 24-hour mean wind speed (1:npatches)

    # --- Radiation inputs ---
    solad_parb::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # [W/m2] Direct beam downwelling (npatches x nbands)
    solai_parb::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # [W/m2] Diffuse downwelling (npatches x nbands)
    filter_vegzen_pa::Vector{Bool} = Bool[]        # Daylight filter for patches
    coszen_pa::Vector{FT} = Float64[]              # [-] Cosine of solar zenith angle (1:npatches)
    fcansno_pa::Vector{FT} = Float64[]             # [-] Fraction of canopy covered by snow (1:npatches)
    albgr_dir_rb::Vector{FT} = Float64[]           # [-] Ground albedo, direct (1:nbands)
    albgr_dif_rb::Vector{FT} = Float64[]           # [-] Ground albedo, diffuse (1:nbands)

    # --- Photosynthesis inputs ---
    filter_photo_pa::Vector{Int} = Int[]           # Photosynthesis filter flag (1=not called, 2=active, 3=done)
    forc_pbot::FT = 0.0                            # [Pa] Atmospheric pressure
    dayl_factor_pa::Vector{FT} = Float64[]         # [-] Daylength scaling factor (1:npatches)
    esat_tv_pa::Vector{FT} = Float64[]             # [Pa] Saturation vapor pressure at t_veg (1:npatches)
    eair_pa::Vector{FT} = Float64[]                # [Pa] Vapor pressure of canopy air (1:npatches)
    oair_pa::Vector{FT} = Float64[]                # [Pa] Atmospheric O2 partial pressure (1:npatches)
    cair_pa::Vector{FT} = Float64[]                # [Pa] Atmospheric CO2 partial pressure (1:npatches)
    rb_pa::Vector{FT} = Float64[]                  # [s/m] Boundary layer resistance (1:npatches)
    t_veg_pa::Vector{FT} = Float64[]               # [K] Vegetation temperature (1:npatches)
    tgcm_pa::Vector{FT} = Float64[]                # [K] Air temp at agcm reference height (1:npatches)

    # --- Soil temperature and moisture ---
    t_soisno_sl::Vector{FT} = Float64[]            # [K] Soil temperature (1:nlevsoil)

    # --- BTRAN / hydrology inputs ---
    smp_sl::Vector{FT} = Float64[]                 # [mm] Soil matric potential, negative (1:nlevsoil)
    salinity_sl::Vector{FT} = Float64[]            # [ppt] Soil salinity (1:nlevsoil)
    eff_porosity_sl::Vector{FT} = Float64[]        # [-] Effective porosity (1:nlevsoil)
    watsat_sl::Vector{FT} = Float64[]              # [-] Volumetric soil water at saturation (1:nlevsoil)
    tempk_sl::Vector{FT} = Float64[]               # [K] Soil temperature for BTRAN (1:nlevsoil)
    h2o_liqvol_sl::Vector{FT} = Float64[]          # [m3/m3] Liquid volume in soil layer (1:nlevsoil)
    filter_btran::Bool = false                     # Site-level filter for BTRAN calculations

    # --- Plant hydraulics inputs ---
    qflx_transp_pa::Vector{FT} = Float64[]         # [mm H2O/s] Transpiration flux (1:npatches)
    swrad_net_pa::Vector{FT} = Float64[]           # [W/m2] Net absorbed shortwave (1:npatches)
    lwrad_net_pa::Vector{FT} = Float64[]           # [W/m2] Net absorbed longwave (1:npatches)
    watsat_sisl::Vector{FT} = Float64[]            # [-] Porosity for hydraulics (1:nlevsoil)
    watres_sisl::Vector{FT} = Float64[]            # [-] Residual soil water (1:nlevsoil)
    sucsat_sisl::Vector{FT} = Float64[]            # [mm] Minimum soil suction (1:nlevsoil)
    bsw_sisl::Vector{FT} = Float64[]              # [-] Clapp-Hornberger "b" (1:nlevsoil)
    hksat_sisl::Vector{FT} = Float64[]            # [mm/s] Hydraulic conductivity at saturation (1:nlevsoil)
    h2o_liq_sisl::Vector{FT} = Float64[]          # [kg/m2] Liquid water mass per layer (1:nlevsoil)
    smpmin_si::FT = 0.0                           # [mm] Minimum soil potential

    # --- Canopy structure / snow ---
    snow_depth_si::FT = 0.0                        # [m] Snow depth in snowy areas
    frac_sno_eff_si::FT = 0.0                      # [-] Fraction of ground covered by snow

    # --- Litter / BGC accounting ---
    max_rooting_depth_index_col::Int = 0           # Deepest soil level index with roots
    tot_het_resp::FT = 0.0                         # [gC/m2/s] Total heterotrophic respiration
    tot_somc::FT = 0.0                             # [gC/m2] Total soil organic matter carbon
    tot_litc::FT = 0.0                             # [gC/m2] Total litter carbon

    # --- Nutrient uptake (accumulated over sub-timesteps) ---
    plant_nh4_uptake_flux::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # [gN/m2/day] NH4 uptake per competitor
    plant_no3_uptake_flux::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # [gN/m2/day] NO3 uptake per competitor
    plant_p_uptake_flux::Matrix{FT} = Matrix{Float64}(undef, 0, 0)    # [gP/m2/day] P uptake per competitor

    # --- Land use / harvest ---
    hlm_harvest_rates::Vector{FT} = Float64[]      # Annual harvest rate per category
    hlm_harvest_catnames::Vector{String} = String[] # Names of harvest categories
    hlm_harvest_units::Int = 0                     # Harvest rate units (area fraction vs carbon)
    site_area::FT = 0.0                            # [m2] Site area for carbon-based harvest
    hlm_luh_states::Vector{FT} = Float64[]         # LUH2 land-use state fractions
    hlm_luh_state_names::Vector{String} = String[]
    hlm_luh_transitions::Vector{FT} = Float64[]    # LUH2 transition rates
    hlm_luh_transition_names::Vector{String} = String[]

    # --- Fixed biogeography ---
    pft_areafrac::Vector{FT} = Float64[]           # [-] Fractional area per PFT
    baregroundfrac::FT = 0.0                       # [-] Bare ground fraction

    # --- Satellite phenology (SP mode) inputs ---
    hlm_sp_tlai::Vector{FT} = Float64[]            # SP total LAI per PFT
    hlm_sp_tsai::Vector{FT} = Float64[]            # SP total SAI per PFT
    hlm_sp_htop::Vector{FT} = Float64[]            # SP canopy height per PFT
end


"""
    FATESBoundaryCondOut{FT}

Output boundary conditions returned from FATES to CLM for a single site (column).
Corresponds to `bc_out_type` in `FatesInterfaceTypesMod.F90`.
"""
Base.@kwdef mutable struct FATESBoundaryCondOut{FT<:Real}

    # --- Sun/shade fractions ---
    fsun_pa::Vector{FT} = Float64[]                # [-] Sunlit fraction of canopy (1:npatches)
    laisun_pa::Vector{FT} = Float64[]              # [m2/m2] Sunlit canopy LAI (1:npatches)
    laisha_pa::Vector{FT} = Float64[]              # [m2/m2] Shaded canopy LAI (1:npatches)

    # --- BTRAN / root water uptake ---
    active_suction_sl::Vector{Bool} = Bool[]       # Soil layers with active water uptake
    rootr_pasl::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # [-] Root fraction by patch x soil layer
    btran_pa::Vector{FT} = Float64[]               # [-] Transpiration wetness factor (1:npatches)

    # --- Stomatal resistance ---
    rssun_pa::Vector{FT} = Float64[]               # [s/m] Sunlit canopy resistance (1:npatches)
    rssha_pa::Vector{FT} = Float64[]               # [s/m] Shaded canopy resistance (1:npatches)

    # --- Canopy radiation outputs ---
    albd_parb::Matrix{FT} = Matrix{Float64}(undef, 0, 0)   # [-] Surface albedo, direct (npatches x nbands)
    albi_parb::Matrix{FT} = Matrix{Float64}(undef, 0, 0)   # [-] Surface albedo, diffuse (npatches x nbands)
    fabd_parb::Matrix{FT} = Matrix{Float64}(undef, 0, 0)   # [-] Canopy-absorbed fraction, direct
    fabi_parb::Matrix{FT} = Matrix{Float64}(undef, 0, 0)   # [-] Canopy-absorbed fraction, diffuse
    ftdd_parb::Matrix{FT} = Matrix{Float64}(undef, 0, 0)   # [-] Below-canopy direct per unit direct
    ftid_parb::Matrix{FT} = Matrix{Float64}(undef, 0, 0)   # [-] Below-canopy diffuse per unit direct
    ftii_parb::Matrix{FT} = Matrix{Float64}(undef, 0, 0)   # [-] Below-canopy diffuse per unit diffuse

    # --- Litter fluxes to soil BGC ---
    litt_flux_cel_c_si::Vector{FT} = Float64[]     # [gC/m3/s] Cellulose carbon litter
    litt_flux_lig_c_si::Vector{FT} = Float64[]     # [gC/m3/s] Lignin carbon litter
    litt_flux_lab_c_si::Vector{FT} = Float64[]     # [gC/m3/s] Labile carbon litter
    litt_flux_cel_n_si::Vector{FT} = Float64[]     # [gN/m3/s] Cellulose nitrogen litter
    litt_flux_lig_n_si::Vector{FT} = Float64[]     # [gN/m3/s] Lignin nitrogen litter
    litt_flux_lab_n_si::Vector{FT} = Float64[]     # [gN/m3/s] Labile nitrogen litter
    litt_flux_cel_p_si::Vector{FT} = Float64[]     # [gP/m3/s] Cellulose phosphorus litter
    litt_flux_lig_p_si::Vector{FT} = Float64[]     # [gP/m3/s] Lignin phosphorus litter
    litt_flux_lab_p_si::Vector{FT} = Float64[]     # [gP/m3/s] Labile phosphorus litter
    litt_flux_ligc_per_n::FT = 0.0                 # [g/g] Lignin-C per total-N in fragmentation flux

    # --- Nutrient competition outputs ---
    num_plant_comps::Int = 0                       # Number of unique competitors
    source_nh4::Vector{FT} = Float64[]             # [gN/m3] FATES NH4 source to mineralization pool
    source_p::Vector{FT} = Float64[]               # [gP/m3] FATES P source to mineralization pool
    veg_rootc::Matrix{FT} = Matrix{Float64}(undef, 0, 0)  # [gC/m3] Fine-root C per competitor
    ft_index::Vector{Int} = Int[]                  # PFT index per competitor
    cn_scalar::Vector{FT} = Float64[]              # C:N scaling factor for root N uptake
    cp_scalar::Vector{FT} = Float64[]              # C:P scaling factor for root P uptake

    # --- CH4 boundary conditions ---
    annavg_agnpp_pa::Vector{FT} = Float64[]        # [gC/m2/s] Annual avg above-ground NPP
    annavg_bgnpp_pa::Vector{FT} = Float64[]        # [gC/m2/s] Annual avg below-ground NPP
    annsum_npp_pa::Vector{FT} = Float64[]          # [gC/m2/yr] Annual sum NPP
    frootc_pa::Vector{FT} = Float64[]              # [gC/m2] Fine root carbon
    root_resp::Vector{FT} = Float64[]              # [gC/m2/s] Root respiration
    rootfr_pa::Matrix{FT} = Matrix{Float64}(undef, 0, 0)   # [-] Rooting fraction with depth
    woody_frac_aere_pa::Vector{FT} = Float64[]     # [-] Woody plant fraction for aerenchyma
    ema_npp::FT = 0.0                              # [gC/m2/s] Site-level smoothed NPP

    # --- Canopy structure outputs ---
    elai_pa::Vector{FT} = Float64[]                # [m2/m2] Exposed leaf area index (1:npatches)
    esai_pa::Vector{FT} = Float64[]                # [m2/m2] Exposed stem area index (1:npatches)
    tlai_pa::Vector{FT} = Float64[]                # [m2/m2] Total leaf area index (1:npatches)
    tsai_pa::Vector{FT} = Float64[]                # [m2/m2] Total stem area index (1:npatches)
    htop_pa::Vector{FT} = Float64[]                # [m] Top of canopy (1:npatches)
    hbot_pa::Vector{FT} = Float64[]                # [m] Bottom of canopy (1:npatches)
    z0m_pa::Vector{FT} = Float64[]                 # [m] Roughness length (1:npatches)
    displa_pa::Vector{FT} = Float64[]              # [m] Displacement height (1:npatches)
    dleaf_pa::Vector{FT} = Float64[]               # [m] Leaf characteristic dimension (1:npatches)
    canopy_fraction_pa::Vector{FT} = Float64[]     # [-] Area fraction per patch (1:npatches)
    frac_veg_nosno_alb_pa::Vector{FT} = Float64[]  # [0/1] Vegetation present flag (1:npatches)
    nocomp_pft_label_pa::Vector{Int} = Int[]       # PFT label per patch (nocomp/SP mode)

    # --- Plant hydraulics outputs ---
    plant_stored_h2o_si::FT = 0.0                  # [kg/m2] Stored water in vegetation
    qflx_soil2root_sisl::Vector{FT} = Float64[]    # [mm/s] Soil-to-root water flux per layer
    qflx_ro_sisl::Vector{FT} = Float64[]           # [mm/s] Runoff from root-soil super-saturation

    # --- LULCC wood products ---
    hrv_deadstemc_to_prod10c::FT = 0.0             # [gC/m2/s] Harvest C to 10-yr product pool
    hrv_deadstemc_to_prod100c::FT = 0.0            # [gC/m2/s] Harvest C to 100-yr product pool
    gpp_site::FT = 0.0                             # [gC/m2/s] Site-level GPP for NBP diagnosis
    ar_site::FT = 0.0                              # [gC/m2/s] Site-level autotrophic respiration
end


"""
    FATESParameterConstants{FT}

Constant parameters set once during initialization for nutrient competition.
Corresponds to `bc_pconst_type` in `FatesInterfaceTypesMod.F90`.
"""
Base.@kwdef mutable struct FATESParameterConstants{FT<:Real}
    max_plant_comps::Int = 0
    vmax_nh4::Vector{FT} = Float64[]               # Max NH4 uptake rate per PFT
    vmax_no3::Vector{FT} = Float64[]               # Max NO3 uptake rate per PFT
    vmax_p::Vector{FT} = Float64[]                 # Max P uptake rate per PFT
    eca_km_nh4::Vector{FT} = Float64[]             # ECA half-sat for NH4 per PFT
    eca_km_no3::Vector{FT} = Float64[]             # ECA half-sat for NO3 per PFT
    eca_km_p::Vector{FT} = Float64[]               # ECA half-sat for P per PFT
    eca_km_ptase::Vector{FT} = Float64[]           # ECA half-sat for phosphatase per PFT
    eca_vmax_ptase::Vector{FT} = Float64[]         # ECA Vmax phosphatase per PFT
    eca_alpha_ptase::Vector{FT} = Float64[]        # ECA alpha phosphatase per PFT
    eca_lambda_ptase::Vector{FT} = Float64[]       # ECA lambda phosphatase per PFT
    eca_plant_escalar::FT = 0.0                    # ECA plant efficiency scalar
    j_uptake::Vector{Int} = Int[]                  # Mapping: decomp layers -> uptake layers
end


# ============================================================================
#                     FATES SITE AND STATE TYPES
# ============================================================================

"""
    FATESSiteMap

Mapping between FATES sites and CLM columns/gridcells.
Corresponds to `f2hmap_type` in `clmfates_interfaceMod.F90`.
"""
Base.@kwdef mutable struct FATESSiteMap
    fcolumn::Vector{Int} = Int[]                   # Column index for each FATES site
    hsites::Vector{Int} = Int[]                    # Site index for each HLM column (sparse, 0=non-site)
end


"""
    FATESState

Simplified container for FATES ecosystem state at a single site.
In the full Fortran FATES, this is the `ed_site_type` with linked-list
patch/cohort structures. Here we provide a minimal placeholder that
tracks the key dimensions and state needed for the CLM interface.

When FATES dynamics are incrementally ported, this struct will be
expanded to hold the full patch/cohort hierarchy.
"""
Base.@kwdef mutable struct FATESState{FT<:Real}
    # --- Site dimensions ---
    nsites::Int = 0                                # Number of sites in this clump

    # --- Per-site boundary conditions ---
    bc_in::Vector{FATESBoundaryCondIn{FT}} = FATESBoundaryCondIn{FT}[]
    bc_out::Vector{FATESBoundaryCondOut{FT}} = FATESBoundaryCondOut{FT}[]
    bc_pconst::Vector{FATESParameterConstants{FT}} = FATESParameterConstants{FT}[]

    # --- Site mapping ---
    sitemap::FATESSiteMap = FATESSiteMap()

    # --- Placeholder for patch/cohort state ---
    # When porting FATES dynamics, add:
    #   sites::Vector{FATESSiteData}   (contains patch linked-list)
    #   Each patch contains cohort linked-list with:
    #     dbh, height, n (density), pft, canopy_layer, carbon pools, etc.
end


"""
    HLMFATESInterface

Top-level FATES interface object held on CLMInstances.
Corresponds to `hlm_fates_interface_type` in `clmfates_interfaceMod.F90`.
"""
Base.@kwdef mutable struct HLMFATESInterface{FT<:Real}
    config::FATESConfig = FATESConfig()
    state::FATESState{FT} = FATESState{FT}()
    initialized::Bool = false
end


# ============================================================================
#                     INITIALIZATION FUNCTIONS
# ============================================================================

"""
    fates_init!(iface::HLMFATESInterface, nsites, nlevsoil, npatches_max, nbands;
                config_kwargs...)

Initialize the FATES interface for a given number of sites.

**What this needs to implement (Fortran `init` subroutine, ~260 lines):**
- Call `FatesInterfaceInit()` to set up FATES logging and verbosity
- Set all FATES control parameters via `set_fates_ctrlparms()`
- Call `SetFatesGlobalElements1/2()` to determine dimension sizes
- Call `allocate_bcin/bcout/bcpconst()` for each site
- Initialize fire data method based on `hlm_spitfire_mode`
- Set up site-to-column mapping (`f2hmap`)
- Initialize FATES history output
- Call `init_soil_depths()` to set vertical grid in bc_in
- Call `init_coldstart()` or `restart()` depending on simulation type
"""
function fates_init!(iface::HLMFATESInterface{FT}, nsites::Int, nlevsoil::Int,
                     npatches_max::Int, nbands::Int;
                     config_kwargs...) where {FT}
    # Apply any config overrides
    for (k, v) in config_kwargs
        if hasproperty(iface.config, k)
            setproperty!(iface.config, k, v)
        end
    end

    # Allocate boundary conditions for each site
    state = iface.state
    state.nsites = nsites
    state.bc_in = [FATESBoundaryCondIn{FT}() for _ in 1:nsites]
    state.bc_out = [FATESBoundaryCondOut{FT}() for _ in 1:nsites]
    state.bc_pconst = [FATESParameterConstants{FT}() for _ in 1:nsites]

    for s in 1:nsites
        _allocate_bc_in!(state.bc_in[s], nlevsoil, npatches_max, nbands)
        _allocate_bc_out!(state.bc_out[s], nlevsoil, npatches_max, nbands)
    end

    iface.initialized = true
    return nothing
end


"""
    _allocate_bc_in!(bc, nlevsoil, npatches, nbands)

Allocate arrays in a `FATESBoundaryCondIn` struct.
"""
function _allocate_bc_in!(bc::FATESBoundaryCondIn{FT}, nlevsoil::Int,
                          npatches::Int, nbands::Int) where {FT}
    bc.nlevsoil = nlevsoil
    bc.nlevdecomp = nlevsoil  # default: every soil layer is biogeochemically active

    # Soil structure
    bc.zi_sisl = zeros(FT, nlevsoil + 1)  # 0:nlevsoil
    bc.dz_sisl = zeros(FT, nlevsoil)
    bc.z_sisl = zeros(FT, nlevsoil)
    bc.dz_decomp_sisl = zeros(FT, nlevsoil)
    bc.decomp_id = collect(1:nlevsoil)

    # Decomposition scalars
    bc.w_scalar_sisl = zeros(FT, nlevsoil)
    bc.t_scalar_sisl = zeros(FT, nlevsoil)

    # Fire
    bc.lightning24 = zeros(FT, npatches)
    bc.pop_density = zeros(FT, npatches)
    bc.precip24_pa = zeros(FT, npatches)
    bc.relhumid24_pa = zeros(FT, npatches)
    bc.wind24_pa = zeros(FT, npatches)

    # Radiation
    bc.solad_parb = zeros(FT, npatches, nbands)
    bc.solai_parb = zeros(FT, npatches, nbands)
    bc.filter_vegzen_pa = falses(npatches)
    bc.coszen_pa = zeros(FT, npatches)
    bc.fcansno_pa = zeros(FT, npatches)
    bc.albgr_dir_rb = zeros(FT, nbands)
    bc.albgr_dif_rb = zeros(FT, nbands)

    # Photosynthesis
    bc.filter_photo_pa = ones(Int, npatches)
    bc.dayl_factor_pa = zeros(FT, npatches)
    bc.esat_tv_pa = zeros(FT, npatches)
    bc.eair_pa = zeros(FT, npatches)
    bc.oair_pa = zeros(FT, npatches)
    bc.cair_pa = zeros(FT, npatches)
    bc.rb_pa = zeros(FT, npatches)
    bc.t_veg_pa = zeros(FT, npatches)
    bc.tgcm_pa = zeros(FT, npatches)

    # Soil temperature / moisture
    bc.t_soisno_sl = zeros(FT, nlevsoil)
    bc.smp_sl = zeros(FT, nlevsoil)
    bc.salinity_sl = zeros(FT, nlevsoil)
    bc.eff_porosity_sl = zeros(FT, nlevsoil)
    bc.watsat_sl = zeros(FT, nlevsoil)
    bc.tempk_sl = zeros(FT, nlevsoil)
    bc.h2o_liqvol_sl = zeros(FT, nlevsoil)

    # Plant hydraulics
    bc.qflx_transp_pa = zeros(FT, npatches)
    bc.swrad_net_pa = zeros(FT, npatches)
    bc.lwrad_net_pa = zeros(FT, npatches)
    bc.watsat_sisl = zeros(FT, nlevsoil)
    bc.watres_sisl = zeros(FT, nlevsoil)
    bc.sucsat_sisl = zeros(FT, nlevsoil)
    bc.bsw_sisl = zeros(FT, nlevsoil)
    bc.hksat_sisl = zeros(FT, nlevsoil)
    bc.h2o_liq_sisl = zeros(FT, nlevsoil)

    # SP mode
    bc.hlm_sp_tlai = zeros(FT, npatches)
    bc.hlm_sp_tsai = zeros(FT, npatches)
    bc.hlm_sp_htop = zeros(FT, npatches)

    return nothing
end


"""
    _allocate_bc_out!(bc, nlevsoil, npatches, nbands)

Allocate arrays in a `FATESBoundaryCondOut` struct.
"""
function _allocate_bc_out!(bc::FATESBoundaryCondOut{FT}, nlevsoil::Int,
                           npatches::Int, nbands::Int) where {FT}
    # Sun/shade
    bc.fsun_pa = zeros(FT, npatches)
    bc.laisun_pa = zeros(FT, npatches)
    bc.laisha_pa = zeros(FT, npatches)

    # BTRAN
    bc.active_suction_sl = falses(nlevsoil)
    bc.rootr_pasl = zeros(FT, npatches, nlevsoil)
    bc.btran_pa = zeros(FT, npatches)

    # Stomatal resistance
    bc.rssun_pa = fill(FT(2.0e4), npatches)  # rsmax0 = 2e4 s/m
    bc.rssha_pa = fill(FT(2.0e4), npatches)

    # Canopy radiation
    bc.albd_parb = zeros(FT, npatches, nbands)
    bc.albi_parb = zeros(FT, npatches, nbands)
    bc.fabd_parb = zeros(FT, npatches, nbands)
    bc.fabi_parb = zeros(FT, npatches, nbands)
    bc.ftdd_parb = zeros(FT, npatches, nbands)
    bc.ftid_parb = zeros(FT, npatches, nbands)
    bc.ftii_parb = zeros(FT, npatches, nbands)

    # Litter fluxes
    bc.litt_flux_cel_c_si = zeros(FT, nlevsoil)
    bc.litt_flux_lig_c_si = zeros(FT, nlevsoil)
    bc.litt_flux_lab_c_si = zeros(FT, nlevsoil)
    bc.litt_flux_cel_n_si = zeros(FT, nlevsoil)
    bc.litt_flux_lig_n_si = zeros(FT, nlevsoil)
    bc.litt_flux_lab_n_si = zeros(FT, nlevsoil)
    bc.litt_flux_cel_p_si = zeros(FT, nlevsoil)
    bc.litt_flux_lig_p_si = zeros(FT, nlevsoil)
    bc.litt_flux_lab_p_si = zeros(FT, nlevsoil)

    # Canopy structure
    bc.elai_pa = zeros(FT, npatches)
    bc.esai_pa = zeros(FT, npatches)
    bc.tlai_pa = zeros(FT, npatches)
    bc.tsai_pa = zeros(FT, npatches)
    bc.htop_pa = zeros(FT, npatches)
    bc.hbot_pa = zeros(FT, npatches)
    bc.z0m_pa = zeros(FT, npatches)
    bc.displa_pa = zeros(FT, npatches)
    bc.dleaf_pa = zeros(FT, npatches)
    bc.canopy_fraction_pa = zeros(FT, npatches)
    bc.frac_veg_nosno_alb_pa = zeros(FT, npatches)
    bc.nocomp_pft_label_pa = zeros(Int, npatches)

    # CH4
    bc.annavg_agnpp_pa = zeros(FT, npatches)
    bc.annavg_bgnpp_pa = zeros(FT, npatches)
    bc.annsum_npp_pa = zeros(FT, npatches)
    bc.frootc_pa = zeros(FT, npatches)
    bc.root_resp = zeros(FT, npatches)
    bc.rootfr_pa = zeros(FT, npatches, nlevsoil)
    bc.woody_frac_aere_pa = zeros(FT, npatches)

    # Plant hydraulics
    bc.qflx_soil2root_sisl = zeros(FT, nlevsoil)
    bc.qflx_ro_sisl = zeros(FT, nlevsoil)

    return nothing
end


# ============================================================================
#               CLM <-> FATES DATA TRANSFER FUNCTIONS
# ============================================================================

"""
    clm_to_fates_soil!(bc_in, col, ss, nlevsoil, c)

Transfer soil structure from CLM column to FATES boundary conditions.

**What this needs to implement (Fortran `init_soil_depths`, ~45 lines):**
- Copy `col.zi`, `col.dz`, `col.z` into `bc_in.zi_sisl/dz_sisl/z_sisl`
- Set `bc_in.dz_decomp_sisl` from `dzsoi_decomp`
- Set `bc_in.decomp_id` mapping (identity when using full vertical resolution)
"""
function clm_to_fates_soil!(bc_in::FATESBoundaryCondIn, col, ss, nlevsoil::Int, c::Int)
    # STUB: No-op. When implemented, copy vertical grid from CLM column c
    # into bc_in soil layer arrays.
    return nothing
end


"""
    clm_to_fates_dynamics!(bc_in, inst, c, npatches, nlevsoil, config;
                           year, month, day, tod)

Transfer all CLM state needed for FATES daily dynamics into bc_in.

**What this needs to implement (Fortran `dynamics_driv` Part I, ~230 lines):**
- Set time information via `GetAndSetTime()`
- Copy decomposition scalars: `w_scalar_sisl`, `t_scalar_sisl`
- Copy soil moisture: `h2o_liqvol_sl`, `tempk_sl`
- Compute `max_rooting_depth_index_col` from active layer
- Compute soil matric potential `smp_sl` via SWRC for active suction layers
- Copy fire inputs: `precip24_pa`, `relhumid24_pa`, `wind24_pa`, `lightning24`, `pop_density`
- Copy SP inputs: `hlm_sp_tlai/tsai/htop` (if use_fates_sp)
- Copy plant hydraulics inputs: `hksat/watsat/watres/sucsat/bsw/h2o_liq` (if use_planthydro)
- Copy harvest and land-use data (if harvest/LUH modes active)
"""
function clm_to_fates_dynamics!(bc_in::FATESBoundaryCondIn, inst, c::Int,
                                npatches::Int, nlevsoil::Int, config;
                                year::Int=0, month::Int=0, day::Int=0, tod::Int=0)
    # STUB: No-op.
    return nothing
end


"""
    clm_to_fates_radiation!(bc_in, alb, a2l, coszen, fcansno, c, npatches, filter_vegsol)

Transfer radiation inputs from CLM to FATES for canopy radiation calculations.

**What this needs to implement (Fortran `wrap_canopy_radiation` input section, ~30 lines):**
- Set `bc_in.filter_vegzen_pa` based on filter membership
- Copy `coszen_pa`, `fcansno_pa` from patch arrays
- Copy `albgr_dir_rb`, `albgr_dif_rb` from column ground albedo
"""
function clm_to_fates_radiation!(bc_in::FATESBoundaryCondIn, alb, a2l,
                                 coszen, fcansno, c::Int, npatches::Int,
                                 filter_vegsol)
    # STUB: No-op.
    return nothing
end


"""
    clm_to_fates_sunfrac!(bc_in, a2l, c, g, npatches)

Transfer downwelling radiation from CLM to FATES for sun/shade fraction calculation.

**What this needs to implement (Fortran `wrap_sunfrac` input section, ~8 lines):**
- Copy `forc_solad_not_downscaled_grc` -> `bc_in.solad_parb`
- Copy `forc_solai_grc` -> `bc_in.solai_parb`
"""
function clm_to_fates_sunfrac!(bc_in::FATESBoundaryCondIn, a2l, c::Int, g::Int,
                                npatches::Int)
    # STUB: No-op.
    return nothing
end


"""
    clm_to_fates_btran!(bc_in, ss, wdb, temp, c, nlevsoil, filter_exposed)

Transfer soil state for BTRAN calculation from CLM to FATES.

**What this needs to implement (Fortran `wrap_btran` input section, ~25 lines):**
- Check if column is in exposed vegetation filter -> set `bc_in.filter_btran`
- Copy `tempk_sl`, `h2o_liqvol_sl`, `eff_porosity_sl`, `watsat_sl`
- After FATES returns active_suction_sl, compute `smp_sl` via SWRC
"""
function clm_to_fates_btran!(bc_in::FATESBoundaryCondIn, ss, wdb, temp,
                              c::Int, nlevsoil::Int, filter_exposed)
    # STUB: No-op.
    return nothing
end


"""
    clm_to_fates_photosynthesis!(bc_in, temp, a2l, esat_tv, eair, oair, cair, rb,
                                  dayl_factor, c, npatches, filterp)

Transfer photosynthesis inputs from CLM to FATES.

**What this needs to implement (Fortran `wrap_photosynthesis` input section, ~35 lines):**
- Copy soil temperatures: `t_soisno_sl`
- Copy atmospheric pressure: `forc_pbot`
- For each patch in filter: copy `esat_tv_pa`, `eair_pa`, `oair_pa`, `cair_pa`,
  `rb_pa`, `t_veg_pa`, `tgcm_pa`, `dayl_factor_pa`
- Set `filter_photo_pa` to 2 for active patches
"""
function clm_to_fates_photosynthesis!(bc_in::FATESBoundaryCondIn, temp, a2l,
                                       esat_tv, eair, oair, cair, rb,
                                       dayl_factor, c::Int, npatches::Int, filterp)
    # STUB: No-op.
    return nothing
end


"""
    clm_to_fates_hydraulics!(bc_in, ss, wsb, wdb, wfb, solarabs, ef,
                              c, nlevsoil, npatches, filterp)

Transfer hydraulics inputs from CLM to FATES plant hydraulics.

**What this needs to implement (Fortran `wrap_hydraulics_drive` input section, ~40 lines):**
- Copy soil properties: `smpmin_si`, `watsat/watres/sucsat/bsw/hksat/h2o_liq/eff_porosity`
- Copy patch-level: `swrad_net_pa`, `lwrad_net_pa`, `qflx_transp_pa`
"""
function clm_to_fates_hydraulics!(bc_in::FATESBoundaryCondIn, ss, wsb, wdb, wfb,
                                   solarabs, ef, c::Int, nlevsoil::Int,
                                   npatches::Int, filterp)
    # STUB: No-op.
    return nothing
end


"""
    clm_to_fates_snow!(bc_in, wdb, c)

Transfer snow state from CLM to FATES for canopy structure updates.

**What this needs to implement (Fortran `wrap_update_hlmfates_dyn` input section, ~5 lines):**
- Copy `snow_depth_col(c)` -> `bc_in.snow_depth_si`
- Copy `frac_sno_eff_col(c)` -> `bc_in.frac_sno_eff_si`
"""
function clm_to_fates_snow!(bc_in::FATESBoundaryCondIn, wdb, c::Int)
    # STUB: No-op.
    return nothing
end


"""
    clm_to_fates_hifreq_hist!(bc_in, soilbgc_cf, soilbgc_cs, c)

Transfer high-frequency history inputs (het resp, soil/litter C) from CLM to FATES.

**What this needs to implement (Fortran `wrap_update_hifrq_hist` input section, ~10 lines):**
- Copy `hr_col(c)` -> `bc_in.tot_het_resp`
- Copy `totsomc_col(c)` -> `bc_in.tot_somc`
- Copy `totlitc_col(c)` -> `bc_in.tot_litc`
"""
function clm_to_fates_hifreq_hist!(bc_in::FATESBoundaryCondIn, soilbgc_cf, soilbgc_cs, c::Int)
    # STUB: No-op.
    return nothing
end


# ============================================================================
#               FATES -> CLM DATA TRANSFER FUNCTIONS
# ============================================================================

"""
    fates_to_clm_sunfrac!(bc_out, cs, c, npatches)

Transfer sun/shade fractions from FATES back to CLM.

**What this needs to implement (Fortran `wrap_sunfrac` output section, ~10 lines):**
- Copy `bc_out.fsun_pa` -> `canopystate.fsun_patch`
- Copy `bc_out.laisun_pa` -> `canopystate.laisun_patch`
- Copy `bc_out.laisha_pa` -> `canopystate.laisha_patch`
"""
function fates_to_clm_sunfrac!(bc_out::FATESBoundaryCondOut, cs, c::Int, npatches::Int)
    # STUB: No-op.
    return nothing
end


"""
    fates_to_clm_btran!(bc_out, ef, ss, c, npatches, nlevsoil)

Transfer BTRAN and root distribution from FATES back to CLM.

**What this needs to implement (Fortran `wrap_btran` output section, ~20 lines):**
- Copy `bc_out.btran_pa(ifp)` -> `energyflux.btran_patch(p)`
- Copy `bc_out.rootr_pasl(ifp,j)` -> `soilstate.rootr_patch(p,j)`
- Set `rresis_patch` to -999.9 (not correctly calculated by FATES)
"""
function fates_to_clm_btran!(bc_out::FATESBoundaryCondOut, ef, ss,
                              c::Int, npatches::Int, nlevsoil::Int)
    # STUB: No-op.
    return nothing
end


"""
    fates_to_clm_photosynthesis!(bc_out, photosyns, c, npatches, filterp)

Transfer stomatal resistance from FATES back to CLM.

**What this needs to implement (Fortran `wrap_photosynthesis` output section, ~15 lines):**
- Copy `bc_out.rssun_pa(ifp)` -> `photosyns.rssun_patch(p)`
- Copy `bc_out.rssha_pa(ifp)` -> `photosyns.rssha_patch(p)`
- Set `psnsun_patch`, `psnsha_patch` to SPVAL (not used by CLM)
- Verify filter_photo_pa transitions from 2 -> 3
"""
function fates_to_clm_photosynthesis!(bc_out::FATESBoundaryCondOut, photosyns,
                                       c::Int, npatches::Int, filterp)
    # STUB: No-op.
    return nothing
end


"""
    fates_to_clm_radiation!(bc_out, alb, c, npatches, filter_vegsol)

Transfer canopy radiation from FATES back to CLM.

**What this needs to implement (Fortran `wrap_canopy_radiation` output section, ~15 lines):**
- Copy `albd/albi/fabd/fabi/ftdd/ftid/ftii` from bc_out to surfalb patch arrays
"""
function fates_to_clm_radiation!(bc_out::FATESBoundaryCondOut, alb,
                                  c::Int, npatches::Int, filter_vegsol)
    # STUB: No-op.
    return nothing
end


"""
    fates_to_clm_canopy_structure!(bc_out, cs, wdb, pch, c, npatches; use_fates_sp=false)

Transfer canopy structure from FATES to CLM patch arrays.

**What this needs to implement (Fortran `wrap_update_hlmfates_dyn` output section, ~140 lines):**
- Update `patch.is_veg`, `patch.is_bareground`, `patch.wt_ed`
- Copy `elai/esai/tlai/tsai/htop/hbot/z0m/displa/dleaf` to CLM patch arrays
- Copy `frac_veg_nosno_alb` to canopy state
- Copy `canopy_fraction_pa` for area weighting
- Handle SP-mode tlai_hist/tsai_hist/htop_hist separately
- Verify area fractions sum to 1.0 within tolerance
"""
function fates_to_clm_canopy_structure!(bc_out::FATESBoundaryCondOut, cs, wdb, pch,
                                         c::Int, npatches::Int; use_fates_sp::Bool=false)
    # STUB: No-op.
    return nothing
end


"""
    fates_to_clm_z0m_displa!(bc_out, z0m_patch, displa_patch, c, npatches)

Transfer roughness length and displacement height from FATES to CLM.

**What this needs to implement (Fortran `TransferZ0mDisp`, ~38 lines):**
- Zero z0m/displa for non-bareground patches
- Copy `bc_out.z0m_pa(ifp)` -> `z0m_patch(p)`
- Copy `bc_out.displa_pa(ifp)` -> `displa_patch(p)`
"""
function fates_to_clm_z0m_displa!(bc_out::FATESBoundaryCondOut, z0m_patch, displa_patch,
                                   c::Int, npatches::Int)
    # STUB: No-op.
    return nothing
end


"""
    fates_to_clm_litter_c!(bc_out, soilbgc_cf, c, nlevsoil)

Transfer carbon litter fluxes from FATES to CLM soil BGC.

**What this needs to implement (Fortran `UpdateCLitterFluxes`, ~74 lines):**
- Map FATES litter fractions (cellulose/lignin/labile) to CLM decomposition pools
- Copy `litt_flux_cel_c_si`, `litt_flux_lig_c_si`, `litt_flux_lab_c_si` ->
  `soilbgc_cf.decomp_cpools_sourcesink_col`
- Copy `litt_flux_ligc_per_n` -> `soilbgc_cf.litr_lig_c_to_n_col`
"""
function fates_to_clm_litter_c!(bc_out::FATESBoundaryCondOut, soilbgc_cf,
                                 c::Int, nlevsoil::Int)
    # STUB: No-op.
    return nothing
end


"""
    fates_to_clm_litter_n!(bc_out, soilbgc_nf, c, nlevsoil)

Transfer nitrogen litter fluxes from FATES to CLM soil BGC.

**What this needs to implement (Fortran `UpdateNLitterFluxes`, ~69 lines):**
- Map FATES N litter fractions to CLM decomposition pools
- Copy `litt_flux_cel_n_si`, `litt_flux_lig_n_si`, `litt_flux_lab_n_si`
"""
function fates_to_clm_litter_n!(bc_out::FATESBoundaryCondOut, soilbgc_nf,
                                 c::Int, nlevsoil::Int)
    # STUB: No-op.
    return nothing
end


"""
    fates_to_clm_hydraulics!(bc_out, wfb, wdb, c, nlevsoil)

Transfer plant hydraulics outputs from FATES to CLM water fluxes.

**What this needs to implement (Fortran `ComputeRootSoilFlux`, ~58 lines + `wrap_hydraulics_drive` output, ~15 lines):**
- Copy `bc_out.qflx_soil2root_sisl` -> `wfb.qflx_rootsoi_col`
- Copy `bc_out.plant_stored_h2o_si` -> `wdb.total_plant_stored_h2o_col`
"""
function fates_to_clm_hydraulics!(bc_out::FATESBoundaryCondOut, wfb, wdb,
                                   c::Int, nlevsoil::Int)
    # STUB: No-op.
    return nothing
end


"""
    fates_to_clm_wood_products!(bc_out, c_products, n_products, c, g)

Transfer harvested wood product fluxes from FATES to CLM product pools.

**What this needs to implement (Fortran `wrap_WoodProducts`, ~50 lines):**
- Accumulate `bc_out.hrv_deadstemc_to_prod10c/100c` into gridcell product pools
"""
function fates_to_clm_wood_products!(bc_out::FATESBoundaryCondOut, c_products,
                                      n_products, c::Int, g::Int)
    # STUB: No-op.
    return nothing
end


# ============================================================================
#                     CORE FATES DYNAMICS STUBS
# ============================================================================

"""
    fates_dynamics!(iface::HLMFATESInterface, inst, config, dtime;
                    is_beg_curr_day, year, month, day, tod)

Main FATES daily dynamics driver. Called once per day from `clm_drv!`.

**What this needs to implement (Fortran `dynamics_driv` Parts I-IV, ~290 lines):**

Part I - Prepare boundary conditions:
- Call `clm_to_fates_dynamics!()` for each site
- Unpack nutrient acquisition BCs: `UnPackNutrientAquisitionBCs()`
- Distribute seeds from neighbors: `WrapUpdateFatesSeedInOut()`

Part II - Call FATES ecosystem dynamics:
- `ed_ecosystem_dynamics()` — the main FATES daily step (~5000+ lines internally):
  - Phenology (leaf on/off, carbon allocation to leaves)
  - Seed germination and recruitment of new cohorts
  - Growth (diameter increment, height allometry, carbon allocation)
  - Mortality (background, carbon starvation, hydraulic failure, fire, logging)
  - Disturbance (fire, treefall, logging) creating new patches
  - Patch fusion (merge similar patches to control complexity)
  - Cohort fusion (merge similar cohorts)
  - Cohort promotion/demotion between canopy layers
- `ed_update_site()` — post-dynamics site-level updates

Part III - Update HLM diagnostics:
- Call `fates_to_clm_canopy_structure!()` to update CLM patch arrays
- Update litter C:N ratio for soil BGC
- Update plant stored water (if hydraulics on)

Part IV - Update FATES history:
- `fates_hist.update_history_dyn()` — accumulate FATES-specific history output
"""
function fates_dynamics!(iface::HLMFATESInterface, inst, config, dtime::Real;
                         is_beg_curr_day::Bool=false, year::Int=0, month::Int=0,
                         day::Int=0, tod::Int=0)
    # STUB: No-op. Only runs on the daily timestep (is_beg_curr_day).
    if !is_beg_curr_day
        return nothing
    end
    # When implemented:
    # 1. clm_to_fates_dynamics! for each site
    # 2. ed_ecosystem_dynamics for each site
    # 3. ed_update_site for each site
    # 4. fates_to_clm_canopy_structure! for each site
    return nothing
end


"""
    fates_update_running_means!(iface::HLMFATESInterface, temp, dtime)

Update FATES running mean temperature for growth calculations.

**What this needs to implement (Fortran `WrapUpdateFatesRmean`, ~18 lines):**
- For each site/patch: update `bc_in.t_veg_pa` from `temp.t_veg_patch`
- Call `UpdateFatesRMeansTStep()` which updates internal FATES running means
  used for cold-hardiness, growth decisions, and phenology triggers
"""
function fates_update_running_means!(iface::HLMFATESInterface, temp, dtime::Real)
    # STUB: No-op.
    return nothing
end


"""
    fates_update_hifreq_hist!(iface::HLMFATESInterface, soilbgc_cf, soilbgc_cs, dtime)

Update FATES high-frequency (sub-daily) history diagnostics.

**What this needs to implement (Fortran `wrap_update_hifrq_hist`, ~50 lines):**
- Transfer het_resp, totsomc, totlitc from CLM to bc_in
- Call `fates_hist.update_history_hifrq()` for sub-daily FATES diagnostics
"""
function fates_update_hifreq_hist!(iface::HLMFATESInterface, soilbgc_cf, soilbgc_cs,
                                    dtime::Real)
    # STUB: No-op.
    return nothing
end


"""
    fates_wrap_sunfrac!(iface::HLMFATESInterface, a2l, cs)

Compute sun/shade fractions for FATES canopy structure.

**What this needs to implement (Fortran `wrap_sunfrac`, ~88 lines):**
1. Transfer CLM radiation inputs to bc_in (solad/solai per patch)
2. Call `FatesSunShadeFracs()` — computes multi-layer canopy sunlit fractions
   using FATES two-stream radiation solution
3. Transfer `fsun_pa`, `laisun_pa`, `laisha_pa` back to CLM canopy state

This replaces the standard CLM `canopy_sun_shade_fracs!()` when use_fates=true.
"""
function fates_wrap_sunfrac!(iface::HLMFATESInterface, a2l, cs)
    # STUB: No-op.
    return nothing
end


"""
    fates_wrap_btran!(iface::HLMFATESInterface, ss, wdb, temp, ef, filter_exposed)

Compute BTRAN (transpiration wetness factor) using FATES root distribution.

**What this needs to implement (Fortran `wrap_btran`, ~167 lines):**
1. Transfer soil state to bc_in (temp, moisture, porosity)
2. Call `get_active_suction_layers()` — FATES determines which layers can uptake water
3. Compute soil suction potential via CLM SWRC for active layers
4. Call `btran_ed()` — FATES BTRAN using PFT-specific root profiles and
   cohort-weighted root water uptake fractions
5. Transfer btran_pa, rootr_pasl back to CLM energy flux / soil state

This replaces the standard CLM `calc_root_moist_stress!()` when use_fates=true.
"""
function fates_wrap_btran!(iface::HLMFATESInterface, ss, wdb, temp, ef, filter_exposed)
    # STUB: No-op.
    return nothing
end


"""
    fates_prep_canopyfluxes!(iface::HLMFATESInterface, photosyns)

Prepare FATES patches for photosynthesis iteration loop.

**What this needs to implement (Fortran `prep_canopyfluxes`, ~35 lines):**
- Reset `filter_photo_pa` to 1 (not yet called) for all patches
- Zero `qflx_transp_pa` if plant hydraulics is on
"""
function fates_prep_canopyfluxes!(iface::HLMFATESInterface, photosyns)
    # STUB: No-op.
    return nothing
end


"""
    fates_wrap_photosynthesis!(iface::HLMFATESInterface, bounds, filterp,
                               esat_tv, eair, oair, cair, rb, dayl_factor,
                               a2l, temp, cs, photosyns, dtime)

Compute photosynthesis for FATES cohorts within patches.

**What this needs to implement (Fortran `wrap_photosynthesis`, ~126 lines):**
1. Transfer atmospheric state to bc_in (pbot, esat_tv, eair, oair, cair, rb, t_veg, tgcm)
2. Set filter_photo_pa = 2 for patches in the active filter
3. Call `FatesPlantRespPhotosynthDrive()` — computes photosynthesis for each
   cohort within each patch using Farquhar model with cohort-specific leaf traits
4. Transfer rssun_pa, rssha_pa back to CLM photosyns
5. Set filter_photo_pa = 3 to mark completion

This replaces the standard CLM `photosynthesis!()` call when use_fates=true.
The stomatal conductance is computed at the cohort level and then aggregated
to a patch-level sunlit/shaded average for the CLM canopy flux iteration.
"""
function fates_wrap_photosynthesis!(iface::HLMFATESInterface, bounds, filterp,
                                    esat_tv, eair, oair, cair, rb, dayl_factor,
                                    a2l, temp, cs, photosyns, dtime::Real)
    # STUB: No-op.
    return nothing
end


"""
    fates_wrap_accumulatefluxes!(iface::HLMFATESInterface, filterp, dtime)

Accumulate FATES carbon/water fluxes over the high-frequency timestep.

**What this needs to implement (Fortran `wrap_accumulatefluxes`, ~37 lines):**
- Verify filter_photo_pa == 3 for all patches in filter
- Call `AccumulateFluxes_ED()` — accumulates GPP, respiration, and transpiration
  from the high-frequency (half-hourly) timestep into daily totals used by
  the daily dynamics step
"""
function fates_wrap_accumulatefluxes!(iface::HLMFATESInterface, filterp, dtime::Real)
    # STUB: No-op.
    return nothing
end


"""
    fates_wrap_canopy_radiation!(iface::HLMFATESInterface, alb, coszen, fcansno,
                                 filter_vegsol)

Compute canopy radiation transfer using FATES two-stream model.

**What this needs to implement (Fortran `wrap_canopy_radiation`, ~90 lines):**
1. Transfer radiation inputs: coszen, fcansno, ground albedos, filter
2. Call `FatesNormalizedCanopyRadiation()` — multi-layer two-stream canopy
   radiative transfer that computes albedo, absorbed, and transmitted
   radiation fractions accounting for crown area, LAI profile, and
   leaf optical properties at each canopy layer
3. Transfer albd/albi/fabd/fabi/ftdd/ftid/ftii back to CLM surface albedo

This replaces the standard CLM `SNICAR + TwoStream` when use_fates=true.
"""
function fates_wrap_canopy_radiation!(iface::HLMFATESInterface, alb, coszen, fcansno,
                                      filter_vegsol)
    # STUB: No-op.
    return nothing
end


"""
    fates_wrap_hydraulics_drive!(iface::HLMFATESInterface, ss, wsb, wdb, wfb,
                                 solarabs, ef, filterp, dtime)

Drive FATES plant hydraulics model.

**What this needs to implement (Fortran `wrap_hydraulics_drive`, ~108 lines):**
1. Transfer soil and atmospheric inputs to bc_in
2. Transfer transpiration flux from canopy fluxes
3. Call `hydraulics_drive()` — computes plant water potential, xylem and
   rhizosphere water transport, and soil-root water exchange for each
   cohort based on its size, root profile, and sapwood area
4. Transfer plant_stored_h2o and qflx_soil2root back to CLM
5. Update hydraulics history diagnostics
"""
function fates_wrap_hydraulics_drive!(iface::HLMFATESInterface, ss, wsb, wdb, wfb,
                                      solarabs, ef, filterp, dtime::Real)
    # STUB: No-op. Only runs when use_fates_planthydro=true.
    return nothing
end


"""
    fates_wrap_seed_dispersal!(iface::HLMFATESInterface; is_restart::Bool=false)

Global seed dispersal across grid cells.

**What this needs to implement (Fortran `WrapGlobalSeedDispersal`, ~89 lines):**
- Only executes at beginning of each day
- Collects outgoing seeds from all sites
- Performs MPI_Allgatherv to distribute globally (in serial Julia: direct copy)
- Computes incoming seeds for each site based on neighbor kernel weights
- Calls `WrapUpdateFatesSeedInOut()` to distribute to local sites

**Note:** In a single-site Julia simulation, this simplifies to a local
seed rain with no inter-gridcell transport.
"""
function fates_wrap_seed_dispersal!(iface::HLMFATESInterface; is_restart::Bool=false)
    # STUB: No-op.
    return nothing
end


"""
    fates_update_acc_vars!(iface::HLMFATESInterface, bounds)

Update FATES accumulation variables (fire weather, etc.).

**What this needs to implement (Fortran `UpdateAccVars`, ~10 lines):**
- Call `fates_fire_data_method.UpdateAccVars()` to update fire weather
  accumulation buffers (lightning, precipitation, humidity, wind)
"""
function fates_update_acc_vars!(iface::HLMFATESInterface, bounds)
    # STUB: No-op.
    return nothing
end


"""
    fates_coldstart!(iface::HLMFATESInterface, wsb, wdb, cs, ss, col, config)

Cold-start FATES initialization.

**What this needs to implement (Fortran `init_coldstart`, ~180 lines):**
- Call `zero_site()` for each site — zero all site-level variables
- Call `init_site_vars()` — initialize site properties (soil, climate)
- Call `set_site_properties()` — set up site characteristics
- Call `init_patches()` — create initial patch structure with:
  - Single near-bare-ground patch (for cold start)
  - Or patch structure from inventory file (for inventory init)
- Initialize plant hydraulics via `HydrSiteColdStart()` if active
- Call `fates_to_clm_canopy_structure!()` to initialize CLM arrays
"""
function fates_coldstart!(iface::HLMFATESInterface, wsb, wdb, cs, ss, col, config)
    # STUB: No-op.
    return nothing
end


"""
    fates_restart!(iface::HLMFATESInterface, ncid, flag)

Read or write FATES restart data.

**What this needs to implement (Fortran `restart`, ~640 lines):**
- Define all restart variables via `restartvar()` calls
- On read: restore full patch/cohort hierarchy from flat restart arrays
- On write: serialize patch/cohort hierarchy to flat restart arrays
- Call `RestartHydrStates()` if plant hydraulics is active
- Call `fates_to_clm_canopy_structure!()` after restart read to
  re-initialize CLM patch arrays from FATES state
"""
function fates_restart!(iface::HLMFATESInterface, ncid, flag::Symbol)
    # STUB: No-op. flag is :read or :write.
    return nothing
end


"""
    fates_transfer_z0m_displa!(iface::HLMFATESInterface, z0m_patch, displa_patch)

Transfer roughness length and displacement height each timestep.

**What this needs to implement (Fortran `TransferZ0mDisp`, ~38 lines):**
- Called every timestep because CLM resets these values
- Copy bc_out.z0m_pa and bc_out.displa_pa to CLM patch arrays
"""
function fates_transfer_z0m_displa!(iface::HLMFATESInterface, z0m_patch, displa_patch)
    # STUB: No-op.
    return nothing
end


"""
    fates_update_litter_fluxes!(iface::HLMFATESInterface, soilbgc_cf, soilbgc_nf)

Transfer FATES litter fluxes to CLM soil BGC model.

**What this needs to implement (Fortran `UpdateCLitterFluxes` + `UpdateNLitterFluxes`, ~143 lines):**
- Map FATES 3-pool litter (cellulose/lignin/labile) to CLM decomposition pool structure
- Update `decomp_cpools_sourcesink_col` for carbon
- Update `decomp_npools_sourcesink_col` for nitrogen
- Apply to all soil layers based on vertical litter distribution
"""
function fates_update_litter_fluxes!(iface::HLMFATESInterface, soilbgc_cf, soilbgc_nf)
    # STUB: No-op.
    return nothing
end


"""
    fates_wrap_wood_products!(iface::HLMFATESInterface, c_products, n_products, bounds)

Transfer FATES harvest wood product fluxes to CLM product pools.

**What this needs to implement (Fortran `wrap_WoodProducts`, ~50 lines):**
- For each site: accumulate hrv_deadstemc_to_prod10c/100c into gridcell product pools
"""
function fates_wrap_wood_products!(iface::HLMFATESInterface, c_products, n_products, bounds)
    # STUB: No-op.
    return nothing
end


"""
    fates_compute_root_soil_flux!(iface::HLMFATESInterface, ss, wfb, filter_soilc)

Transfer FATES root-soil water flux to CLM soil hydrology.

**What this needs to implement (Fortran `ComputeRootSoilFlux`, ~58 lines):**
- Only active when use_fates_planthydro=true
- Copy bc_out.qflx_soil2root_sisl -> wfb.qflx_rootsoi_col for each site
"""
function fates_compute_root_soil_flux!(iface::HLMFATESInterface, ss, wfb, filter_soilc)
    # STUB: No-op.
    return nothing
end


"""
    fates_sp_phenology!(iface::HLMFATESInterface, cs, wdb, pch, filter)

FATES Satellite Phenology (SP) mode phenology update.

**What this needs to implement:**
- When use_fates_sp=true and doalb=true, apply prescribed LAI/SAI/height
  from the HLM surface dataset to FATES patches
- This is analogous to the standard CLM `satellite_phenology!()` but routes
  through FATES patch structure so that the SP values are properly distributed
  to FATES cohorts
"""
function fates_sp_phenology!(iface::HLMFATESInterface, cs, wdb, pch, filter)
    # STUB: No-op.
    return nothing
end
