# ==========================================================================
# Ported from: src/main/clm_varctl.F90 (585 lines)
# Control variables and flags for CLM configuration
# ==========================================================================

# Undefined sentinel values
const IUNDEF = -9999999
const RUNDEF = -9999999.0

# Filename length
const FNAME_LEN = 256

# Startup type constants
const NSRSTARTUP  = 0
const NSRCONTINUE = 1
const NSRBRANCH   = 2

"""
Configuration struct holding all CLM control variables.
Set once at initialization, constant during simulation.
Mirrors Fortran module-level variables in clm_varctl.F90.
"""
Base.@kwdef mutable struct VarCtl
    # --- Case identification ---
    caseid::String = ""
    ctitle::String = ""
    nsrest::Int = IUNDEF
    is_cold_start::Bool = false
    brnch_retain_casename::Bool = false
    version::String = ""
    hostname::String = ""
    username::String = ""
    source::String = "Community Terrestrial Systems Model"
    conventions::String = "CF-1.0"
    compname::String = "clm2"

    # --- I/O control ---
    iulog::Int = 6  # stdout
    outnc_large_files::Bool = true

    # --- Input file paths ---
    finidat::String = ""
    fsurdat::String = ""
    hillslope_file::String = ""
    paramfile::String = ""
    nrevsn::String = ""
    fsnowoptics::String = ""
    fsnowaging::String = ""
    fatmlndfrc::String = ""
    finidat_interp_source::String = ""
    finidat_interp_dest::String = ""

    # --- Atmospheric forcing ---
    ndep_from_cpl::Bool = false
    co2_type::String = "constant"
    co2_ppmv::Float64 = 355.0

    # --- Hydrology control ---
    bound_h2osoi::Bool = true
    irrigate::Bool = false
    crop_fsat_equals_zero::Bool = false
    crop_residue_removal_frac::Float64 = 0.0

    # --- Crop/vegetation control ---
    use_crop::Bool = false
    create_crop_landunit::Bool = false

    # --- Subgrid control ---
    all_active::Bool = false
    convert_ocean_to_land::Bool = false
    collapse_urban::Bool = false
    run_zero_weight_urban::Bool = false
    n_dom_landunits::Int = -1
    n_dom_pfts::Int = -1
    use_subgrid_fluxes::Bool = true

    # --- Landunit size thresholds ---
    toosmall_soil::Float64 = -1.0
    toosmall_crop::Float64 = -1.0
    toosmall_glacier::Float64 = -1.0
    toosmall_lake::Float64 = -1.0
    toosmall_wetland::Float64 = -1.0
    toosmall_urban::Float64 = -1.0

    # --- Hillslope control ---
    nhillslope::Int = 0
    max_columns_hillslope::Int = 1
    use_hillslope::Bool = false
    downscale_hillslope_meteorology::Bool = false
    use_hillslope_routing::Bool = false
    hillslope_fsat_equals_zero::Bool = false

    # --- Biogeochemistry control ---
    spinup_state::Int = 0
    anoxia::Bool = true
    override_bgc_restart_mismatch_dump::Bool = false
    nfix_timeconst::Float64 = -1.2345

    # --- Carbon/nitrogen model flags ---
    use_cn::Bool = false
    use_cndv::Bool = false
    use_lch4::Bool = true
    use_nitrif_denitrif::Bool = true
    use_grainproduct::Bool = false
    use_fertilizer::Bool = false
    use_flexibleCN::Bool = false
    MM_Nuptake_opt::Bool = false
    CNratio_floating::Bool = false
    lnc_opt::Bool = false
    reduce_dayl_factor::Bool = false
    vcmax_opt::Int = 0
    CN_evergreen_phenology_opt::Int = 0
    carbon_resp_opt::Int = 0

    # --- Isotope flags ---
    use_c13::Bool = false
    use_c14::Bool = false

    # --- Spinup matrix method ---
    spinup_matrixcn::Bool = false
    hist_wrt_matrixcn_diag::Bool = false
    nyr_forcing::Int = 10
    nyr_SASU::Int = 1
    iloop_avg::Int = -999

    # --- Ozone ---
    o3_veg_stress_method::String = "unset"
    o3_ppbv::Float64 = 100.0

    # --- Snow optics (SNICAR) ---
    snicar_numrad_snw::Int = 5
    snicar_solarspec::String = "mid_latitude_winter"
    snicar_dust_optics::String = "sahara"
    snicar_use_aerosol::Bool = true
    snicar_snw_shape::String = "hexagonal_plate"
    snicar_snobc_intmix::Bool = false
    snicar_snodst_intmix::Bool = false
    do_sno_oc::Bool = false
    use_snicar_frc::Bool = false
    use_SSRE::Bool = false

    # --- Snow cover and thermal ---
    snow_cover_fraction_method::String = ""
    snow_thermal_cond_method::String = "Jordan1991"
    use_z0m_snowmelt::Bool = false
    z0param_method::String = ""
    h2osno_max::Float64 = -999.0

    # --- LUNA ---
    use_luna::Bool = false

    # --- FATES ---
    use_fates::Bool = false
    use_fates_sp::Bool = false
    use_fates_bgc::Bool = false
    use_fates_planthydro::Bool = false
    use_fates_tree_damage::Bool = false
    use_fates_cohort_age_tracking::Bool = false
    use_fates_ed_st3::Bool = false
    use_fates_ed_prescribed_phys::Bool = false
    use_fates_inventory_init::Bool = false
    use_fates_fixed_biogeog::Bool = false
    use_fates_nocomp::Bool = false
    use_fates_luh::Bool = false
    use_fates_lupft::Bool = false
    use_fates_potentialveg::Bool = false
    fates_seeddisp_cadence::Int = 0
    fates_parteh_mode::Int = 0
    fates_spitfire_mode::Int = 0
    fates_harvest_mode::String = ""
    fates_history_dimlevel::Int = 0
    fluh_timeseries::String = ""
    flandusepftdat::String = ""
    fates_inventory_ctrl_filename::String = ""
    fates_paramfile::String = ""

    # --- Streams ---
    use_soil_moisture_streams::Bool = false
    use_lai_streams::Bool = false
    use_cropcal_streams::Bool = false
    use_cropcal_rx_swindows::Bool = false
    use_cropcal_rx_cultivar_gdds::Bool = false
    adapt_cropcal_rx_cultivar_gdds::Bool = false
    flush_gdd20::Bool = false

    # --- Miscellaneous ---
    use_biomass_heat_storage::Bool = false
    use_bedrock::Bool = false
    use_excess_ice::Bool = false
    use_hydrstress::Bool = false
    use_extralakelayers::Bool = false
    use_vichydro::Bool = false
    use_nguardrail::Bool = false
    use_noio::Bool = false
    use_vancouver::Bool = false
    use_mexicocity::Bool = false
    hist_wrtch4diag::Bool = false
    hist_fields_list_file::Bool = false

    # --- Soil layer structure ---
    soil_layerstruct_predefined::String = "UNSET"
    soil_layerstruct_userdefined::Vector{Float64} = fill(RUNDEF, 99)
    soil_layerstruct_userdefined_nlevsoi::Int = IUNDEF

    # --- Glacier dynamics ---
    glc_do_dynglacier::Bool = false
    glc_snow_persistence_max_days::Int = 7300

    # --- Single column mode ---
    single_column::Bool = false
    scmlat::Float64 = RUNDEF
    scmlon::Float64 = RUNDEF

    # --- Instance info (multi-instance) ---
    inst_index::Int = 1
    inst_name::String = ""
    inst_suffix::String = ""

    # --- Decomposition ---
    nsegspc::Int = 20

    # --- Restart pointer ---
    rpntdir::String = "."
    rpntfil::String = "rpointer.lnd"

    # --- Testing flags ---
    for_testing_run_ncdiopio_tests::Bool = false
    for_testing_use_second_grain_pool::Bool = false
    for_testing_use_repr_structure_pool::Bool = false
    for_testing_no_crop_seed_replenishment::Bool = false
    for_testing_allow_interp_non_ciso_to_ciso::Bool = false
end

# Global configuration instance
const varctl = VarCtl()
