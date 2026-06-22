# EDParamsMod.jl
# Julia port of FATES src/fates/main/EDParamsMod.F90
#
# ED-model scalar (and a few small vector) parameters, plus their
# registration with / retrieval from the FATES parameters interface
# (fates_parameters_type) and init/report API.
#
# In Fortran these are module-level `protected, save` variables. Here we collect
# the *settable* ones into a single mutable holder (`EDParams`, a const Ref to an
# `ed_params_type` instance) so they can be re-initialized at runtime, mirroring
# the FatesGlobals Ref-holder pattern. True compile-time `parameter`s (e.g.
# soil_tfrz_thresh, nclmax, nlevleaf, maxpft) stay as `const`.
#
# The parameter-file READ stays stubbed (no NetCDF), exactly like the foundation
# FatesParametersInterface: FatesRegisterParams! registers the names + shapes,
# FatesReceiveParams! retrieves values from an already-populated
# fates_parameters_type, and the parameter names are preserved verbatim.
#
# Deps: FatesConstantsMod (fates_r8, nearzero, fates_unset_r8,
#       cstarvation_model_lin, n_landuse_cats), FatesGlobals (fates_log,
#       fates_endrun), FatesParametersInterface (fates_parameters_type,
#       param_string_length + RegisterParameter!/RetrieveParameter*).

# ---------------------------------------------------------------------------
# Compile-time parameters (Fortran `parameter`, public) -> const
# ---------------------------------------------------------------------------

# Soil temperature threshold below which hydraulic failure mortality is off
# (non-hydro only) in degrees C
const soil_tfrz_thresh = -2.0

# Maximum number of canopy layers (used only for scratch arrays)
const nclmax = 2

# Number of leaf+stem layers in each canopy layer (VAI/LAI+SAI bins)
const nlevleaf = 30

# Maximum number of PFTs allowed
const maxpft = 16

# ---------------------------------------------------------------------------
# Parameter NAMES (Fortran character `parameter`s) -> const Strings.
# Preserve the FATES parameter-file names verbatim.
# ---------------------------------------------------------------------------
const fates_name_active_crown_fire          = "fates_fire_active_crown_fire"
const fates_name_cg_strikes                  = "fates_fire_cg_strikes"
const name_dev_arbitrary                     = "fates_dev_arbitrary"

const ED_name_sdlng_emerg_h2o_timescale      = "fates_trs_seedling_emerg_h2o_timescale"
const ED_name_sdlng_mort_par_timescale       = "fates_trs_seedling_mort_par_timescale"
const ED_name_sdlng_mdd_timescale            = "fates_trs_seedling_mdd_timescale"
const ED_name_sdlng2sap_par_timescale        = "fates_trs_seedling2sap_par_timescale"
const ED_name_photo_temp_acclim_timescale    = "fates_leaf_photo_temp_acclim_timescale"
const ED_name_photo_temp_acclim_thome_time   = "fates_leaf_photo_temp_acclim_thome_time"
const name_photo_tempsens_model              = "fates_leaf_photo_tempsens_model"
const name_maintresp_model                   = "fates_maintresp_leaf_model"
const name_radiation_model                   = "fates_rad_model"
const ED_name_hydr_htftype_node              = "fates_hydro_htftype_node"
const ED_name_mort_disturb_frac              = "fates_mort_disturb_frac"
const ED_name_mort_cstarvation_model         = "fates_mort_cstarvation_model"
const ED_name_comp_excln                     = "fates_comp_excln"
const ED_name_vai_top_bin_width              = "fates_vai_top_bin_width"
const ED_name_vai_width_increase_factor      = "fates_vai_width_increase_factor"
const ED_name_nignitions                     = "fates_fire_nignitions"
const ED_name_understorey_death              = "fates_mort_understorey_death"
const ED_name_cwd_fcel                       = "fates_frag_cwd_fcel"
const ED_name_cwd_flig                       = "fates_frag_cwd_flig"
const fates_name_maintresp_nonleaf_baserate  = "fates_maintresp_nonleaf_baserate"
const ED_name_phen_a                         = "fates_phen_gddthresh_a"
const ED_name_phen_b                         = "fates_phen_gddthresh_b"
const ED_name_phen_c                         = "fates_phen_gddthresh_c"
const ED_name_phen_chiltemp                  = "fates_phen_chilltemp"
const ED_name_phen_mindayson                 = "fates_phen_mindayson"
const ED_name_phen_ncolddayslim              = "fates_phen_ncolddayslim"
const ED_name_phen_coldtemp                  = "fates_phen_coldtemp"
const ED_name_cohort_size_fusion_tol         = "fates_cohort_size_fusion_tol"
const ED_name_cohort_age_fusion_tol          = "fates_cohort_age_fusion_tol"
const ED_name_patch_fusion_tol               = "fates_patch_fusion_tol"
const ED_name_canopy_closure_thresh          = "fates_canopy_closure_thresh"
const ED_name_stomatal_model                 = "fates_leaf_stomatal_model"
const ED_name_dayl_switch                    = "fates_daylength_factor_switch"
const ED_name_regeneration_model             = "fates_regeneration_model"

const name_theta_cj_c3                       = "fates_leaf_theta_cj_c3"
const name_theta_cj_c4                       = "fates_leaf_theta_cj_c4"

const fates_name_q10_mr                      = "fates_q10_mr"
const fates_name_q10_froz                    = "fates_q10_froz"

# non-scalar parameter names
const ED_name_history_sizeclass_bin_edges    = "fates_history_sizeclass_bin_edges"
const ED_name_history_ageclass_bin_edges     = "fates_history_ageclass_bin_edges"
const ED_name_history_height_bin_edges       = "fates_history_height_bin_edges"
const ED_name_history_coageclass_bin_edges   = "fates_history_coageclass_bin_edges"
const ED_name_history_damage_bin_edges       = "fates_history_damage_bin_edges"
const ED_name_crop_lu_pft_vector             = "fates_landuse_crop_lu_pft_vector"
const ED_name_maxpatches_by_landuse          = "fates_maxpatches_by_landuse"
const ED_name_max_nocomp_pfts_by_landuse     = "fates_max_nocomp_pfts_by_landuse"

# Hydraulics control parameter names
const hydr_name_kmax_rsurf1                  = "fates_hydro_kmax_rsurf1"
const hydr_name_kmax_rsurf2                  = "fates_hydro_kmax_rsurf2"
const hydr_name_psi0                         = "fates_hydro_psi0"
const hydr_name_psicap                       = "fates_hydro_psicap"
const hydr_name_solver                       = "fates_hydro_solver"

# Soil BGC parameter names
const bgc_name_soil_salinity                 = "fates_soil_salinity"

const stomatal_assim_name                    = "fates_leaf_stomatal_assim_model"
const damage_name_event_code                 = "fates_damage_event_code"
const damage_name_canopy_layer_code          = "fates_damage_canopy_layer_code"
const maxcohort_name                         = "fates_maxcohort"

# Logging control parameter names
const logging_name_dbhmin                    = "fates_landuse_logging_dbhmin"
const logging_name_dbhmax                    = "fates_landuse_logging_dbhmax"
const logging_name_collateral_frac           = "fates_landuse_logging_collateral_frac"
const logging_name_coll_under_frac           = "fates_landuse_logging_coll_under_frac"
const logging_name_direct_frac               = "fates_landuse_logging_direct_frac"
const logging_name_mechanical_frac           = "fates_landuse_logging_mechanical_frac"
const logging_name_event_code                = "fates_landuse_logging_event_code"
const logging_name_dbhmax_infra              = "fates_landuse_logging_dbhmax_infra"
const logging_name_export_frac               = "fates_landuse_logging_export_frac"

const eca_name_plant_escalar                 = "fates_cnp_eca_plant_escalar"

# ---------------------------------------------------------------------------
# ed_params_type — the settable ED-model parameters.
#
# Holds every Fortran module variable that is set inside FatesParamsInit /
# FatesReceiveParams. Reals default to `fates_unset_r8` and integer switches to
# the Fortran `-9` sentinel (matching FatesParamsInit). The const `EDParams`
# Ref below is the live, module-global instance the API mutates.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct ed_params_type
    # Radiative-transfer / VAI binning
    vai_top_bin_width::Float64                 = fates_unset_r8
    vai_width_increase_factor::Float64         = fates_unset_r8

    # Photosynthesis / respiration acclimation
    photo_temp_acclim_timescale::Float64       = fates_unset_r8
    photo_temp_acclim_thome_time::Float64      = fates_unset_r8
    photo_tempsens_model::Int                  = -9
    maintresp_leaf_model::Int                  = -9
    radiation_model::Int                       = -9

    # Seedling (tree recruitment scheme) timescales
    sdlng_emerg_h2o_timescale::Float64         = fates_unset_r8
    sdlng_mort_par_timescale::Float64          = fates_unset_r8
    sdlng_mdd_timescale::Float64               = fates_unset_r8
    sdlng2sap_par_timescale::Float64           = fates_unset_r8

    # Mortality / disturbance
    mort_cstarvation_model::Int                = -9
    fates_mortality_disturbance_fraction::Float64 = fates_unset_r8

    # ED_val_* tuned parameters
    ED_val_comp_excln::Float64                 = fates_unset_r8
    ED_val_vai_top_bin_width::Float64          = fates_unset_r8
    ED_val_vai_width_increase_factor::Float64  = fates_unset_r8
    ED_val_nignitions::Float64                 = fates_unset_r8
    ED_val_understorey_death::Float64          = fates_unset_r8
    ED_val_cwd_fcel::Float64                   = fates_unset_r8
    ED_val_cwd_flig::Float64                   = fates_unset_r8
    maintresp_nonleaf_baserate::Float64        = fates_unset_r8
    ED_val_phen_a::Float64                     = fates_unset_r8
    ED_val_phen_b::Float64                     = fates_unset_r8
    ED_val_phen_c::Float64                     = fates_unset_r8
    ED_val_phen_chiltemp::Float64              = fates_unset_r8
    ED_val_phen_mindayson::Float64             = fates_unset_r8
    ED_val_phen_ncolddayslim::Float64          = fates_unset_r8
    ED_val_phen_coldtemp::Float64              = fates_unset_r8
    ED_val_cohort_size_fusion_tol::Float64     = fates_unset_r8
    ED_val_cohort_age_fusion_tol::Float64      = fates_unset_r8
    ED_val_patch_fusion_tol::Float64           = fates_unset_r8
    ED_val_canopy_closure_thresh::Float64      = fates_unset_r8

    # Stomatal / day-length / regeneration switches
    stomatal_model::Int                        = -9
    dayl_switch::Int                           = -9
    regeneration_model::Int                    = -9
    stomatal_assim_model::Int                  = -9

    # Cohort packing
    max_cohort_per_patch::Int                  = -9

    # Fire
    active_crown_fire::Bool                    = false
    cg_strikes::Float64                        = fates_unset_r8

    # Co-limitation curvature
    theta_cj_c3::Float64                       = fates_unset_r8
    theta_cj_c4::Float64                       = fates_unset_r8

    # Nutrient interaction modes (set by the HLM, not the param file)
    n_uptake_mode::Int                         = fates_unset_int
    p_uptake_mode::Int                         = fates_unset_int

    # Respiration Q10s
    q10_mr::Float64                            = fates_unset_r8
    q10_froz::Float64                          = fates_unset_r8

    # Developer free parameter
    dev_arbitrary::Float64                     = fates_unset_r8

    # VAI bin arrays (sized nlevleaf)
    dinc_vai::Vector{Float64}                  = fill(fates_unset_r8, nlevleaf)
    dlower_vai::Vector{Float64}                = fill(fates_unset_r8, nlevleaf)

    # Hydraulics
    hydr_kmax_rsurf1::Float64                  = fates_unset_r8
    hydr_kmax_rsurf2::Float64                  = fates_unset_r8
    hydr_psi0::Float64                         = fates_unset_r8
    hydr_psicap::Float64                       = fates_unset_r8
    hydr_solver::Int                           = -9

    # Soil BGC
    bgc_soil_salinity::Float64                 = fates_unset_r8

    # Damage
    damage_event_code::Int                     = -9
    damage_canopy_layer_code::Int              = -9

    # Patch/landuse allocation (sized n_landuse_cats)
    maxpatches_by_landuse::Vector{Int}         = zeros(Int, n_landuse_cats)
    max_nocomp_pfts_by_landuse::Vector{Int}    = zeros(Int, n_landuse_cats)
    maxpatch_total::Int                        = 0
    crop_lu_pft_vector::Vector{Int}            = zeros(Int, n_landuse_cats)

    # Logging
    logging_dbhmin::Float64                    = fates_unset_r8
    logging_dbhmax::Float64                    = fates_unset_r8
    logging_collateral_frac::Float64           = fates_unset_r8
    logging_coll_under_frac::Float64           = fates_unset_r8
    logging_direct_frac::Float64               = fates_unset_r8
    logging_mechanical_frac::Float64           = fates_unset_r8
    logging_event_code::Float64                = fates_unset_r8
    logging_dbhmax_infra::Float64              = fates_unset_r8
    logging_export_frac::Float64               = fates_unset_r8

    # ECA
    eca_plant_escalar::Float64                 = fates_unset_r8

    # Allocatable history-bin edges + hydro node types (sized from the param file
    # at read; empty until then).
    ED_val_history_sizeclass_bin_edges::Vector{Float64}  = Float64[]
    ED_val_history_ageclass_bin_edges::Vector{Float64}   = Float64[]
    ED_val_history_height_bin_edges::Vector{Float64}     = Float64[]
    ED_val_history_coageclass_bin_edges::Vector{Float64} = Float64[]
    ED_val_history_damage_bin_edges::Vector{Float64}     = Float64[]
    hydr_htftype_node::Vector{Int}                       = Int[]
end

# The live module-global parameter instance the API mutates (Fortran module
# `save` state). A const Ref so the held object can be swapped/reset.
const EDParams = Ref{ed_params_type}(ed_params_type())

# Convenience accessor for the live parameter object.
ed_params() = EDParams[]

# =====================================================================================

"""
    FatesParamsInit!()

Reset all settable ED parameters to their unset sentinels (reals ->
`fates_unset_r8`, integer switches -> -9), mirroring the Fortran
FatesParamsInit which sets everything to NaN/-9. Implemented by installing a
freshly-defaulted `ed_params_type`.
"""
function FatesParamsInit!()
    EDParams[] = ed_params_type()
    return nothing
end

# =====================================================================================

"""
    FatesRegisterParams!(fates_params::fates_parameters_type)

Register the parameters FATES wants the host to provide. Calls FatesParamsInit!
first (matching Fortran), then registers each scalar and 1-D parameter name with
its dimension shape on the given `fates_parameters_type`.
"""
function FatesRegisterParams!(fates_params::fates_parameters_type)

    dim_names_scalar       = [dimension_name_scalar]
    dim_names_sizeclass    = [dimension_name_history_size_bins]
    dim_names_ageclass     = [dimension_name_history_age_bins]
    dim_names_height       = [dimension_name_history_height_bins]
    dim_names_coageclass   = [dimension_name_history_coage_bins]
    dim_names_hydro_organs = [dimension_name_hydr_organs]
    dim_names_damageclass  = [dimension_name_history_damage_bins]
    dim_names_landuse      = [dimension_name_landuse]

    FatesParamsInit!()

    # --- scalar parameters ---
    scalar_names = (
        ED_name_photo_temp_acclim_timescale,
        ED_name_sdlng_emerg_h2o_timescale,
        ED_name_sdlng_mort_par_timescale,
        ED_name_sdlng_mdd_timescale,
        ED_name_sdlng2sap_par_timescale,
        ED_name_photo_temp_acclim_thome_time,
        name_photo_tempsens_model,
        name_radiation_model,
        name_maintresp_model,
        name_theta_cj_c3,
        name_theta_cj_c4,
        ED_name_mort_disturb_frac,
        ED_name_mort_cstarvation_model,
        ED_name_comp_excln,
        ED_name_vai_top_bin_width,
        ED_name_vai_width_increase_factor,
        ED_name_nignitions,
        ED_name_understorey_death,
        ED_name_cwd_fcel,
        ED_name_cwd_flig,
        fates_name_maintresp_nonleaf_baserate,
        ED_name_phen_a,
        ED_name_phen_b,
        ED_name_phen_c,
        ED_name_phen_chiltemp,
        ED_name_phen_mindayson,
        ED_name_phen_ncolddayslim,
        ED_name_phen_coldtemp,
        ED_name_cohort_size_fusion_tol,
        ED_name_cohort_age_fusion_tol,
        ED_name_patch_fusion_tol,
        ED_name_canopy_closure_thresh,
        ED_name_stomatal_model,
        ED_name_dayl_switch,
        ED_name_regeneration_model,
        stomatal_assim_name,
        maxcohort_name,
        hydr_name_solver,
        hydr_name_kmax_rsurf1,
        hydr_name_kmax_rsurf2,
        hydr_name_psi0,
        hydr_name_psicap,
        bgc_name_soil_salinity,
        logging_name_dbhmin,
        logging_name_dbhmax,
        logging_name_collateral_frac,
        logging_name_coll_under_frac,
        logging_name_direct_frac,
        logging_name_mechanical_frac,
        logging_name_event_code,
        logging_name_dbhmax_infra,
        logging_name_export_frac,
        eca_name_plant_escalar,
        fates_name_q10_mr,
        fates_name_q10_froz,
        name_dev_arbitrary,
        damage_name_event_code,
        damage_name_canopy_layer_code,
        fates_name_active_crown_fire,
        fates_name_cg_strikes,
    )
    for nm in scalar_names
        RegisterParameter!(fates_params, nm, dimension_shape_scalar, dim_names_scalar)
    end

    # --- non-scalar (1-D) parameters ---
    RegisterParameter!(fates_params, ED_name_hydr_htftype_node, dimension_shape_1d, dim_names_hydro_organs)
    RegisterParameter!(fates_params, ED_name_history_sizeclass_bin_edges, dimension_shape_1d, dim_names_sizeclass)
    RegisterParameter!(fates_params, ED_name_history_ageclass_bin_edges, dimension_shape_1d, dim_names_ageclass)
    RegisterParameter!(fates_params, ED_name_history_height_bin_edges, dimension_shape_1d, dim_names_height)
    RegisterParameter!(fates_params, ED_name_history_coageclass_bin_edges, dimension_shape_1d, dim_names_coageclass)
    RegisterParameter!(fates_params, ED_name_history_damage_bin_edges, dimension_shape_1d, dim_names_damageclass)
    RegisterParameter!(fates_params, ED_name_crop_lu_pft_vector, dimension_shape_1d, dim_names_landuse)
    RegisterParameter!(fates_params, ED_name_maxpatches_by_landuse, dimension_shape_1d, dim_names_landuse)
    RegisterParameter!(fates_params, ED_name_max_nocomp_pfts_by_landuse, dimension_shape_1d, dim_names_landuse)

    return nothing
end

# =====================================================================================

"""
    FatesReceiveParams!(fates_params::fates_parameters_type)

Retrieve the previously-registered parameter values out of an already-populated
`fates_parameters_type` into the live `EDParams` instance. Integer-coded
parameters are read as reals and rounded (Fortran `nint`); `active_crown_fire`
is the boolean test `abs(val-1) < nearzero`. The allocatable history-bin edges
and hydro node types are copied/rounded from the stored vectors.
"""
function FatesReceiveParams!(fates_params::fates_parameters_type)
    p = EDParams[]

    p.photo_temp_acclim_timescale = RetrieveParameterScalar(fates_params, ED_name_photo_temp_acclim_timescale)
    p.sdlng_emerg_h2o_timescale   = RetrieveParameterScalar(fates_params, ED_name_sdlng_emerg_h2o_timescale)
    p.sdlng_mort_par_timescale    = RetrieveParameterScalar(fates_params, ED_name_sdlng_mort_par_timescale)
    p.sdlng_mdd_timescale         = RetrieveParameterScalar(fates_params, ED_name_sdlng_mdd_timescale)
    p.sdlng2sap_par_timescale     = RetrieveParameterScalar(fates_params, ED_name_sdlng2sap_par_timescale)
    p.photo_temp_acclim_thome_time = RetrieveParameterScalar(fates_params, ED_name_photo_temp_acclim_thome_time)

    p.photo_tempsens_model = round(Int, RetrieveParameterScalar(fates_params, name_photo_tempsens_model))
    p.radiation_model      = round(Int, RetrieveParameterScalar(fates_params, name_radiation_model))
    p.maintresp_leaf_model = round(Int, RetrieveParameterScalar(fates_params, name_maintresp_model))

    p.fates_mortality_disturbance_fraction = RetrieveParameterScalar(fates_params, ED_name_mort_disturb_frac)
    p.mort_cstarvation_model = round(Int, RetrieveParameterScalar(fates_params, ED_name_mort_cstarvation_model))

    p.ED_val_comp_excln                = RetrieveParameterScalar(fates_params, ED_name_comp_excln)
    p.ED_val_vai_top_bin_width         = RetrieveParameterScalar(fates_params, ED_name_vai_top_bin_width)
    p.ED_val_vai_width_increase_factor = RetrieveParameterScalar(fates_params, ED_name_vai_width_increase_factor)
    p.ED_val_nignitions                = RetrieveParameterScalar(fates_params, ED_name_nignitions)
    p.ED_val_understorey_death         = RetrieveParameterScalar(fates_params, ED_name_understorey_death)
    p.ED_val_cwd_fcel                  = RetrieveParameterScalar(fates_params, ED_name_cwd_fcel)
    p.ED_val_cwd_flig                  = RetrieveParameterScalar(fates_params, ED_name_cwd_flig)
    p.maintresp_nonleaf_baserate       = RetrieveParameterScalar(fates_params, fates_name_maintresp_nonleaf_baserate)
    p.ED_val_phen_a                    = RetrieveParameterScalar(fates_params, ED_name_phen_a)
    p.ED_val_phen_b                    = RetrieveParameterScalar(fates_params, ED_name_phen_b)
    p.ED_val_phen_c                    = RetrieveParameterScalar(fates_params, ED_name_phen_c)
    p.ED_val_phen_chiltemp             = RetrieveParameterScalar(fates_params, ED_name_phen_chiltemp)
    p.ED_val_phen_mindayson            = RetrieveParameterScalar(fates_params, ED_name_phen_mindayson)
    p.ED_val_phen_ncolddayslim         = RetrieveParameterScalar(fates_params, ED_name_phen_ncolddayslim)
    p.ED_val_phen_coldtemp             = RetrieveParameterScalar(fates_params, ED_name_phen_coldtemp)
    p.ED_val_cohort_size_fusion_tol    = RetrieveParameterScalar(fates_params, ED_name_cohort_size_fusion_tol)
    p.ED_val_cohort_age_fusion_tol     = RetrieveParameterScalar(fates_params, ED_name_cohort_age_fusion_tol)
    p.ED_val_patch_fusion_tol          = RetrieveParameterScalar(fates_params, ED_name_patch_fusion_tol)
    p.ED_val_canopy_closure_thresh     = RetrieveParameterScalar(fates_params, ED_name_canopy_closure_thresh)

    p.stomatal_model       = round(Int, RetrieveParameterScalar(fates_params, ED_name_stomatal_model))
    p.dayl_switch          = round(Int, RetrieveParameterScalar(fates_params, ED_name_dayl_switch))
    p.regeneration_model   = round(Int, RetrieveParameterScalar(fates_params, ED_name_regeneration_model))
    p.stomatal_assim_model = round(Int, RetrieveParameterScalar(fates_params, stomatal_assim_name))
    p.max_cohort_per_patch = round(Int, RetrieveParameterScalar(fates_params, maxcohort_name))

    p.hydr_kmax_rsurf1 = RetrieveParameterScalar(fates_params, hydr_name_kmax_rsurf1)
    p.hydr_kmax_rsurf2 = RetrieveParameterScalar(fates_params, hydr_name_kmax_rsurf2)
    p.hydr_psi0        = RetrieveParameterScalar(fates_params, hydr_name_psi0)
    p.hydr_psicap      = RetrieveParameterScalar(fates_params, hydr_name_psicap)
    p.hydr_solver      = round(Int, RetrieveParameterScalar(fates_params, hydr_name_solver))

    p.bgc_soil_salinity = RetrieveParameterScalar(fates_params, bgc_name_soil_salinity)

    p.logging_dbhmin          = RetrieveParameterScalar(fates_params, logging_name_dbhmin)
    p.logging_dbhmax          = RetrieveParameterScalar(fates_params, logging_name_dbhmax)
    p.logging_collateral_frac = RetrieveParameterScalar(fates_params, logging_name_collateral_frac)
    p.logging_coll_under_frac = RetrieveParameterScalar(fates_params, logging_name_coll_under_frac)
    p.logging_direct_frac     = RetrieveParameterScalar(fates_params, logging_name_direct_frac)
    p.logging_mechanical_frac = RetrieveParameterScalar(fates_params, logging_name_mechanical_frac)
    p.logging_event_code      = RetrieveParameterScalar(fates_params, logging_name_event_code)
    p.logging_dbhmax_infra    = RetrieveParameterScalar(fates_params, logging_name_dbhmax_infra)
    p.logging_export_frac     = RetrieveParameterScalar(fates_params, logging_name_export_frac)

    p.eca_plant_escalar = RetrieveParameterScalar(fates_params, eca_name_plant_escalar)

    p.theta_cj_c3 = RetrieveParameterScalar(fates_params, name_theta_cj_c3)
    p.theta_cj_c4 = RetrieveParameterScalar(fates_params, name_theta_cj_c4)

    p.q10_mr   = RetrieveParameterScalar(fates_params, fates_name_q10_mr)
    p.q10_froz = RetrieveParameterScalar(fates_params, fates_name_q10_froz)

    p.dev_arbitrary = RetrieveParameterScalar(fates_params, name_dev_arbitrary)

    tmpreal = RetrieveParameterScalar(fates_params, fates_name_active_crown_fire)
    p.active_crown_fire = (abs(tmpreal - 1.0) < nearzero)

    p.cg_strikes = RetrieveParameterScalar(fates_params, fates_name_cg_strikes)

    p.damage_event_code        = round(Int, RetrieveParameterScalar(fates_params, damage_name_event_code))
    p.damage_canopy_layer_code = round(Int, RetrieveParameterScalar(fates_params, damage_name_canopy_layer_code))

    # parameters that are arrays of size defined within the params file
    p.ED_val_history_sizeclass_bin_edges  = RetrieveParameter1DAllocate(fates_params, ED_name_history_sizeclass_bin_edges)
    p.ED_val_history_ageclass_bin_edges   = RetrieveParameter1DAllocate(fates_params, ED_name_history_ageclass_bin_edges)
    p.ED_val_history_height_bin_edges     = RetrieveParameter1DAllocate(fates_params, ED_name_history_height_bin_edges)
    p.ED_val_history_coageclass_bin_edges = RetrieveParameter1DAllocate(fates_params, ED_name_history_coageclass_bin_edges)
    p.ED_val_history_damage_bin_edges     = RetrieveParameter1DAllocate(fates_params, ED_name_history_damage_bin_edges)

    tmp_vector_by_landuse1 = RetrieveParameter1DAllocate(fates_params, ED_name_crop_lu_pft_vector)
    p.crop_lu_pft_vector = round.(Int, tmp_vector_by_landuse1)

    tmp_vector_by_landuse2 = RetrieveParameter1DAllocate(fates_params, ED_name_maxpatches_by_landuse)
    p.maxpatches_by_landuse = round.(Int, tmp_vector_by_landuse2)
    p.maxpatch_total = sum(p.maxpatches_by_landuse)

    tmp_vector_by_landuse3 = RetrieveParameter1DAllocate(fates_params, ED_name_max_nocomp_pfts_by_landuse)
    p.max_nocomp_pfts_by_landuse = round.(Int, tmp_vector_by_landuse3)

    hydr_htftype_real = RetrieveParameter1DAllocate(fates_params, ED_name_hydr_htftype_node)
    p.hydr_htftype_node = round.(Int, hydr_htftype_real)

    return nothing
end

# =====================================================================================

"""
    FatesReportParams(is_master::Bool; debug_report::Bool=false)

Write the FATES scalar parameters to the log. As in Fortran, gated by both
`is_master` and a `debug_report` flag (default false), so this is a no-op unless
explicitly enabled.
"""
function FatesReportParams(is_master::Bool; debug_report::Bool=false)
    (debug_report && is_master) || return nothing
    p = EDParams[]
    println("-----------  FATES Scalar Parameters -----------------")
    println("vai_top_bin_width = ", p.vai_top_bin_width)
    println("vai_width_increase_factor = ", p.vai_width_increase_factor)
    println("photo_temp_acclim_timescale = ", p.photo_temp_acclim_timescale)
    println("sdlng_emerg_h2o_timescale = ", p.sdlng_emerg_h2o_timescale)
    println("sdlng_mort_par_timescale = ", p.sdlng_mort_par_timescale)
    println("sdlng_mdd_timescale = ", p.sdlng_mdd_timescale)
    println("sdlng2sap_par_timescale = ", p.sdlng2sap_par_timescale)
    println("photo_temp_acclim_thome_time (years) = ", p.photo_temp_acclim_thome_time)
    println("hydr_htftype_node = ", p.hydr_htftype_node)
    println("fates_mortality_disturbance_fraction = ", p.fates_mortality_disturbance_fraction)
    println("ED_val_comp_excln = ", p.ED_val_comp_excln)
    println("ED_val_nignitions = ", p.ED_val_nignitions)
    println("ED_val_understorey_death = ", p.ED_val_understorey_death)
    println("ED_val_cwd_fcel = ", p.ED_val_cwd_fcel)
    println("ED_val_cwd_flig = ", p.ED_val_cwd_flig)
    println("maintresp_nonleaf_baserate = ", p.maintresp_nonleaf_baserate)
    println("ED_val_phen_a = ", p.ED_val_phen_a)
    println("ED_val_phen_b = ", p.ED_val_phen_b)
    println("ED_val_phen_c = ", p.ED_val_phen_c)
    println("ED_val_phen_chiltemp = ", p.ED_val_phen_chiltemp)
    println("ED_val_phen_mindayson = ", p.ED_val_phen_mindayson)
    println("ED_val_phen_ncolddayslim = ", p.ED_val_phen_ncolddayslim)
    println("ED_val_phen_coldtemp = ", p.ED_val_phen_coldtemp)
    println("ED_val_cohort_size_fusion_tol = ", p.ED_val_cohort_size_fusion_tol)
    println("ED_val_cohort_age_fusion_tol = ", p.ED_val_cohort_age_fusion_tol)
    println("ED_val_patch_fusion_tol = ", p.ED_val_patch_fusion_tol)
    println("ED_val_canopy_closure_thresh = ", p.ED_val_canopy_closure_thresh)
    println("regeneration_model = ", p.regeneration_model)
    println("dayl_switch = ", p.dayl_switch)
    println("stomatal_model = ", p.stomatal_model)
    println("stomatal_assim_model = ", p.stomatal_assim_model)
    println("hydro_kmax_rsurf1 = ", p.hydr_kmax_rsurf1)
    println("hydro_kmax_rsurf2 = ", p.hydr_kmax_rsurf2)
    println("hydro_psi0 = ", p.hydr_psi0)
    println("hydro_psicap = ", p.hydr_psicap)
    println("hydro_solver = ", p.hydr_solver)
    println("bgc_soil_salinity = ", p.bgc_soil_salinity)
    println("logging_dbhmin = ", p.logging_dbhmin)
    println("logging_dbhmax = ", p.logging_dbhmax)
    println("logging_collateral_frac = ", p.logging_collateral_frac)
    println("logging_coll_under_frac = ", p.logging_coll_under_frac)
    println("logging_direct_frac = ", p.logging_direct_frac)
    println("logging_mechanical_frac = ", p.logging_mechanical_frac)
    println("logging_event_code = ", p.logging_event_code)
    println("logging_dbhmax_infra = ", p.logging_dbhmax_infra)
    println("eca_plant_escalar = ", p.eca_plant_escalar)
    println("q10_mr = ", p.q10_mr)
    println("q10_froz = ", p.q10_froz)
    println("cg_strikes = ", p.cg_strikes)
    println("active_crown_fire = ", p.active_crown_fire)
    println("damage_event_code = ", p.damage_event_code)
    println("damage_canopy_layer_code = ", p.damage_canopy_layer_code)
    println("------------------------------------------------------")
    return nothing
end
