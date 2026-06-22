# EDPftvarcon.jl
# Julia port of FATES src/fates/main/EDPftvarcon.F90
#
# The FATES PFT-indexed parameter container. Holds the large set of per-PFT
# trait / allometry / allocation / mortality / phenology / fire / hydraulics
# parameter arrays, their registration with the FATES parameter reader, their
# retrieval, the parameter-consistency (`FatesCheckParams`) checks, the optional
# debug report, and the `GetDecompyFrac` litter-decomposability lookup.
#
# This is FATES's OWN PFT param table (EDPftvarcon_inst) — NOT the host CLM
# `pftcon`. It is standalone: nothing here is added to CLMInstances or any
# dual-copied struct.
#
# Translation notes:
#   * fates_r8 -> Float64; per-PFT arrays -> Vector{Float64}; PFT x {leafage,
#     hydr-organ, swb} arrays -> Matrix{Float64} (SoA, PFT is dim 1).
#   * The Fortran `type(EDPftvarcon_type),public :: EDPftvarcon_inst` module
#     singleton -> a const Ref{EDPftvarcon_type} (`EDPftvarcon_inst`), mirroring
#     the EDParamsMod `EDParams` Ref pattern. `edpftvarcon_inst()` accesses it.
#   * The parameter-file READ stays stubbed: Register! registers names+shapes on
#     the FatesParametersInterface registry; Receive! retrieves values from an
#     already-populated fates_parameters_type. Parameter name strings ("fates_*")
#     are preserved verbatim.
#   * The Fortran type-bound procedures (Init/Register/Receive + private
#     Register_PFT/Receive_PFT/..._numrad/..._hydr_organs/..._leafage) become
#     plain functions dispatching on EDPftvarcon_type.
#   * endrun(...) -> Julia error(...).
#   * Module-level EDParamsMod scalars that are protected/save in Fortran are, in
#     this port, fields of the live `ed_params()` instance (radiation_model,
#     dayl_switch, regeneration_model, logging_*_frac, max_nocomp_pfts_by_landuse,
#     maxpatches_by_landuse) — read via ed_params() here.
#
# Deps: FatesRadiationMemMod (num_swb,ivis,inir,norman_solver,twostr_solver),
#       FatesConstantsMod (nearzero,itrue,ifalse,fates_check_param_set,
#         t_water_freeze_k_1atm,lmr_r_1,lmr_r_2,n_landuse_cats,
#         default_regeneration,TRS_regeneration,TRS_no_seedling_dyn),
#       PRTParametersMod (prt_params), FatesGlobals (fates_log,fates_endrun),
#       FatesLitterMod (ilabile,icellulose,ilignin),
#       PRTGenericMod (leaf_organ,fnrt_organ,prt_cnp_flex_allom_hyp,
#         prt_carbon_allom_hyp),
#       FatesInterfaceTypesMod (hlm_nitrogen_spec,hlm_phosphorus_spec,
#         hlm_parteh_mode,hlm_nu_com,hlm_use_fixed_biogeog,hlm_use_sp,
#         hlm_use_inventory_init,hlm_use_nocomp),
#       EDParamsMod (ed_params + regeneration_model field),
#       FatesParametersInterface (fates_parameters_type + Register/Retrieve API).

# Lower-bound constants (Fortran `parameter, public`). The Julia arrays are
# always 1-based; these are kept for traceability and used as the numrad/hydr
# bounds origin.
const lower_bound_pft = 1
const lower_bound_general = 1

# ---------------------------------------------------------------------------
# EDPftvarcon_type — the per-PFT FATES parameter table.
#
# Every Fortran `real(r8), allocatable :: x(:)` -> `x::Vector{Float64}` (empty
# until Receive!); every `x(:,:)` -> `x::Matrix{Float64}`. Field names mirror the
# Fortran members verbatim.
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct EDPftvarcon_type
    # Core trait / turbulence / stomatal
    freezetol::Vector{Float64}          = Float64[]   # minimum temperature tolerance
    hgt_min::Vector{Float64}            = Float64[]   # sapling height (m)
    dleaf::Vector{Float64}              = Float64[]   # leaf characteristic dimension (m)
    z0mr::Vector{Float64}               = Float64[]   # roughness-length / height ratio (-)
    displar::Vector{Float64}            = Float64[]   # displacement / canopy-top-height ratio
    bark_scaler::Vector{Float64}        = Float64[]   # dbh->bark-thickness scaler (fire)
    crown_kill::Vector{Float64}         = Float64[]   # fire-death scaler
    initd::Vector{Float64}              = Float64[]   # initial seedling density
    seed_suppl::Vector{Float64}         = Float64[]   # external seed supplement
    bb_slope::Vector{Float64}           = Float64[]   # Ball-Berry slope
    medlyn_slope::Vector{Float64}       = Float64[]   # Medlyn slope (KPa^0.5)
    stomatal_intercept::Vector{Float64} = Float64[]   # stomatal-conductance intercept

    # Litter fractions
    lf_flab::Vector{Float64} = Float64[]  # Leaf litter labile fraction
    lf_fcel::Vector{Float64} = Float64[]  # Leaf litter cellulose fraction
    lf_flig::Vector{Float64} = Float64[]  # Leaf litter lignin fraction
    fr_flab::Vector{Float64} = Float64[]  # Fine-root litter labile fraction
    fr_fcel::Vector{Float64} = Float64[]  # Fine-root litter cellulose fraction
    fr_flig::Vector{Float64} = Float64[]  # Fine-root litter lignin fraction

    xl::Vector{Float64}             = Float64[]  # Leaf-stem orientation index
    clumping_index::Vector{Float64} = Float64[]  # self-occlusion clumping factor
    c3psn::Vector{Float64}          = Float64[]  # photosynthetic pathway (C4=0, C3=1)
    smpso::Vector{Float64}          = Float64[]  # soil water potential at full stomatal open (mm)
    smpsc::Vector{Float64}          = Float64[]  # soil water potential at full stomatal close (mm)

    # Maintenance-respiration reduction / base rates
    maintresp_reduction_curvature::Vector{Float64}     = Float64[]
    maintresp_reduction_intercept::Vector{Float64}     = Float64[]
    maintresp_reduction_upthresh::Vector{Float64}      = Float64[]
    maintresp_leaf_atkin2017_baserate::Vector{Float64} = Float64[]
    maintresp_leaf_ryan1991_baserate::Vector{Float64}  = Float64[]
    maintresp_leaf_vert_scaler_coeff1::Vector{Float64} = Float64[]
    maintresp_leaf_vert_scaler_coeff2::Vector{Float64} = Float64[]

    # Mortality
    bmort::Vector{Float64}                   = Float64[]
    mort_ip_size_senescence::Vector{Float64} = Float64[]
    mort_r_size_senescence::Vector{Float64}  = Float64[]
    mort_ip_age_senescence::Vector{Float64}  = Float64[]
    mort_r_age_senescence::Vector{Float64}   = Float64[]
    mort_scalar_coldstress::Vector{Float64}  = Float64[]
    mort_scalar_cstarvation::Vector{Float64} = Float64[]
    mort_scalar_hydrfailure::Vector{Float64} = Float64[]
    mort_upthresh_cstarvation::Vector{Float64} = Float64[]
    hf_sm_threshold::Vector{Float64}         = Float64[]
    hf_flc_threshold::Vector{Float64}        = Float64[]

    # Photosynthesis temperature response
    vcmaxha::Vector{Float64} = Float64[]
    jmaxha::Vector{Float64}  = Float64[]
    vcmaxhd::Vector{Float64} = Float64[]
    jmaxhd::Vector{Float64}  = Float64[]
    vcmaxse::Vector{Float64} = Float64[]
    jmaxse::Vector{Float64}  = Float64[]

    # Seed / regeneration
    germination_rate::Vector{Float64}         = Float64[]
    seed_decay_rate::Vector{Float64}          = Float64[]
    seed_dispersal_pdf_scale::Vector{Float64} = Float64[]
    seed_dispersal_pdf_shape::Vector{Float64} = Float64[]
    seed_dispersal_max_dist::Vector{Float64}  = Float64[]
    seed_dispersal_fraction::Vector{Float64}  = Float64[]
    repro_frac_seed::Vector{Float64}          = Float64[]
    a_emerg::Vector{Float64}                  = Float64[]
    b_emerg::Vector{Float64}                  = Float64[]
    par_crit_germ::Vector{Float64}            = Float64[]
    seedling_psi_emerg::Vector{Float64}       = Float64[]
    seedling_psi_crit::Vector{Float64}        = Float64[]
    seedling_light_rec_a::Vector{Float64}     = Float64[]
    seedling_light_rec_b::Vector{Float64}     = Float64[]
    seedling_mdd_crit::Vector{Float64}        = Float64[]
    seedling_h2o_mort_a::Vector{Float64}      = Float64[]
    seedling_h2o_mort_b::Vector{Float64}      = Float64[]
    seedling_h2o_mort_c::Vector{Float64}      = Float64[]
    seedling_root_depth::Vector{Float64}      = Float64[]
    seedling_light_mort_a::Vector{Float64}    = Float64[]
    seedling_light_mort_b::Vector{Float64}    = Float64[]
    background_seedling_mort::Vector{Float64} = Float64[]

    # Trimming / radiation (PFT x num_swb)
    trim_limit::Vector{Float64} = Float64[]
    trim_inc::Vector{Float64}   = Float64[]
    rhol::Matrix{Float64}       = Matrix{Float64}(undef, 0, 0)  # Leaf reflectance (vis,nir)
    rhos::Matrix{Float64}       = Matrix{Float64}(undef, 0, 0)  # Stem reflectance (vis,nir)
    taul::Matrix{Float64}       = Matrix{Float64}(undef, 0, 0)  # Leaf transmittance (vis,nir)
    taus::Matrix{Float64}       = Matrix{Float64}(undef, 0, 0)  # Stem transmittance (vis,nir)

    # Fire
    fire_alpha_SH::Vector{Float64} = Float64[]   # SPITFIRE scorch-height alpha

    # Non-PARTEH allometry
    allom_frbstor_repro::Vector{Float64} = Float64[]  # bstore fraction to reproduction at mortality

    # Prescribed physiology mode
    prescribed_npp_canopy::Vector{Float64}           = Float64[]
    prescribed_npp_understory::Vector{Float64}       = Float64[]
    prescribed_mortality_canopy::Vector{Float64}     = Float64[]
    prescribed_mortality_understory::Vector{Float64} = Float64[]
    prescribed_recruitment::Vector{Float64}          = Float64[]

    # Damage
    damage_frac::Vector{Float64}            = Float64[]
    damage_mort_p1::Vector{Float64}         = Float64[]
    damage_mort_p2::Vector{Float64}         = Float64[]
    damage_recovery_scalar::Vector{Float64} = Float64[]

    # Nutrient acquisition (ECA & RD)
    decompmicc::Vector{Float64}       = Float64[]
    vmax_nh4::Vector{Float64}         = Float64[]
    vmax_no3::Vector{Float64}         = Float64[]
    vmax_p::Vector{Float64}           = Float64[]
    eca_km_nh4::Vector{Float64}       = Float64[]
    eca_km_no3::Vector{Float64}       = Float64[]
    eca_km_p::Vector{Float64}         = Float64[]
    eca_km_ptase::Vector{Float64}     = Float64[]
    eca_vmax_ptase::Vector{Float64}   = Float64[]
    eca_alpha_ptase::Vector{Float64}  = Float64[]
    eca_lambda_ptase::Vector{Float64} = Float64[]

    # Phenology
    phenflush_fraction::Vector{Float64}       = Float64[]
    phen_cold_size_threshold::Vector{Float64} = Float64[]

    # Prescribed N/P uptake
    prescribed_nuptake::Vector{Float64} = Float64[]
    prescribed_puptake::Vector{Float64} = Float64[]

    # Developer free PFT parameter
    dev_arbitrary_pft::Vector{Float64} = Float64[]

    # PFT x leaf-age
    vcmax25top::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # vcmax25, canopy top, by leaf age

    # Plant hydraulics — PFT dimension
    hydr_p_taper::Vector{Float64}    = Float64[]
    hydr_rs2::Vector{Float64}        = Float64[]
    hydr_srl::Vector{Float64}        = Float64[]
    hydr_rfrac_stem::Vector{Float64} = Float64[]
    hydr_avuln_gs::Vector{Float64}   = Float64[]
    hydr_p50_gs::Vector{Float64}     = Float64[]
    hydr_k_lwp::Vector{Float64}      = Float64[]

    # Plant hydraulics — PFT x organ (1=leaf,2=stem,3=transp root,4=abs root)
    hydr_vg_alpha_node::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    hydr_vg_m_node::Matrix{Float64}     = Matrix{Float64}(undef, 0, 0)
    hydr_vg_n_node::Matrix{Float64}     = Matrix{Float64}(undef, 0, 0)
    hydr_avuln_node::Matrix{Float64}    = Matrix{Float64}(undef, 0, 0)
    hydr_p50_node::Matrix{Float64}      = Matrix{Float64}(undef, 0, 0)
    hydr_epsil_node::Matrix{Float64}    = Matrix{Float64}(undef, 0, 0)
    hydr_pitlp_node::Matrix{Float64}    = Matrix{Float64}(undef, 0, 0)
    hydr_fcap_node::Matrix{Float64}     = Matrix{Float64}(undef, 0, 0)
    hydr_pinot_node::Matrix{Float64}    = Matrix{Float64}(undef, 0, 0)
    hydr_kmax_node::Matrix{Float64}     = Matrix{Float64}(undef, 0, 0)
    hydr_resid_node::Matrix{Float64}    = Matrix{Float64}(undef, 0, 0)
    hydr_thetas_node::Matrix{Float64}   = Matrix{Float64}(undef, 0, 0)

    # HLM-PFT -> FATES-PFT area-fraction map (PFT x HLM-PFT)
    hlm_pft_map::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # Land-use / land-use-change PFT parameters
    harvest_pprod10::Vector{Float64}          = Float64[]
    landusechange_frac_burned::Vector{Float64}   = Float64[]
    landusechange_frac_exported::Vector{Float64} = Float64[]
    landusechange_pprod10::Vector{Float64}       = Float64[]
end

# The live module-global instance (Fortran `EDPftvarcon_inst`). A const Ref so
# the held object can be swapped/reset (e.g. across tests / re-init).
const EDPftvarcon_inst = Ref{EDPftvarcon_type}(EDPftvarcon_type())

# Convenience accessor for the live PFT-parameter object.
edpftvarcon_inst() = EDPftvarcon_inst[]

# =====================================================================================

"""
    EDpftconInit!(this::EDPftvarcon_type)

Init hook (Fortran `EDpftconInit`). The Fortran body is empty; kept for API
parity.
"""
function EDpftconInit!(this::EDPftvarcon_type)
    return nothing
end

# =====================================================================================
# Registration. Registers PFT-dimensioned, numrad, hydr-organ, and leaf-age
# parameters with the FATES parameter registry. Names are verbatim.
# =====================================================================================

"""
    Register!(this::EDPftvarcon_type, fates_params::fates_parameters_type)

Register every EDPftvarcon parameter (1-D PFT, numrad-as-1-D, PFT x organ 2-D,
PFT x leaf-age 2-D) with `fates_params`.
"""
function Register!(this::EDPftvarcon_type, fates_params::fates_parameters_type)
    Register_PFT!(this, fates_params)
    Register_PFT_numrad!(this, fates_params)
    Register_PFT_hydr_organs!(this, fates_params)
    Register_PFT_leafage!(this, fates_params)
    return nothing
end

"""
    Receive!(this::EDPftvarcon_type, fates_params::fates_parameters_type)

Retrieve all EDPftvarcon parameter values out of an already-populated
`fates_params` into `this`.
"""
function Receive!(this::EDPftvarcon_type, fates_params::fates_parameters_type)
    Receive_PFT!(this, fates_params)
    Receive_PFT_numrad!(this, fates_params)
    Receive_PFT_hydr_organs!(this, fates_params)
    Receive_PFT_leafage!(this, fates_params)
    return nothing
end

# -------------------------------------------------------------------------------------

# The complete, ordered list of 1-D PFT-dimensioned parameter name strings. The
# order matches the Fortran Register_PFT registration order (significant only for
# parameter-registry indexing, but kept faithful).
const _edpft_pft_param_names = (
    "fates_mort_freezetol",
    "fates_recruit_height_min",
    "fates_fire_bark_scaler",
    "fates_fire_crown_kill",
    "fates_recruit_init_density",
    "fates_recruit_seed_supplement",
    "fates_leaf_stomatal_slope_ballberry",
    "fates_leaf_stomatal_slope_medlyn",
    "fates_leaf_stomatal_intercept",
    "fates_frag_leaf_flab",
    "fates_frag_leaf_fcel",
    "fates_frag_leaf_flig",
    "fates_frag_fnrt_flab",
    "fates_frag_fnrt_fcel",
    "fates_frag_fnrt_flig",
    "fates_rad_leaf_xl",
    "fates_rad_leaf_clumping_index",
    "fates_leaf_c3psn",
    "fates_nonhydro_smpso",
    "fates_nonhydro_smpsc",
    "fates_maintresp_reduction_curvature",
    "fates_maintresp_reduction_intercept",
    "fates_maintresp_reduction_upthresh",
    "fates_maintresp_leaf_atkin2017_baserate",
    "fates_maintresp_leaf_ryan1991_baserate",
    "fates_maintresp_leaf_vert_scaler_coeff1",
    "fates_maintresp_leaf_vert_scaler_coeff2",
    "fates_prescribed_npp_canopy",
    "fates_prescribed_npp_understory",
    "fates_mort_prescribed_canopy",
    "fates_mort_prescribed_understory",
    "fates_recruit_prescribed_rate",
    "fates_damage_frac",
    "fates_damage_mort_p1",
    "fates_damage_mort_p2",
    "fates_damage_recovery_scalar",
    "fates_fire_alpha_SH",
    "fates_allom_frbstor_repro",
    "fates_hydro_p_taper",
    "fates_hydro_rs2",
    "fates_hydro_srl",
    "fates_hydro_rfrac_stem",
    "fates_hydro_avuln_gs",
    "fates_hydro_p50_gs",
    "fates_hydro_k_lwp",
    "fates_mort_bmort",
    "fates_mort_r_size_senescence",
    "fates_mort_ip_size_senescence",
    "fates_mort_r_age_senescence",
    "fates_mort_ip_age_senescence",
    "fates_mort_scalar_coldstress",
    "fates_mort_scalar_cstarvation",
    "fates_mort_scalar_hydrfailure",
    "fates_mort_upthresh_cstarvation",
    "fates_mort_hf_sm_threshold",
    "fates_mort_hf_flc_threshold",
    "fates_leaf_vcmaxha",
    "fates_leaf_jmaxha",
    "fates_leaf_vcmaxhd",
    "fates_leaf_jmaxhd",
    "fates_leaf_vcmaxse",
    "fates_leaf_jmaxse",
    "fates_recruit_seed_germination_rate",
    "fates_trs_repro_frac_seed",
    "fates_trs_seedling_a_emerg",
    "fates_trs_seedling_b_emerg",
    "fates_trs_seedling_par_crit_germ",
    "fates_trs_seedling_psi_emerg",
    "fates_trs_seedling_psi_crit",
    "fates_trs_seedling_light_rec_a",
    "fates_trs_seedling_light_rec_b",
    "fates_trs_seedling_mdd_crit",
    "fates_trs_seedling_h2o_mort_a",
    "fates_trs_seedling_h2o_mort_b",
    "fates_trs_seedling_h2o_mort_c",
    "fates_trs_seedling_root_depth",
    "fates_trs_seedling_light_mort_a",
    "fates_trs_seedling_light_mort_b",
    "fates_trs_seedling_background_mort",
    "fates_frag_seed_decay_rate",
    "fates_seed_dispersal_pdf_scale",
    "fates_seed_dispersal_pdf_shape",
    "fates_seed_dispersal_max_dist",
    "fates_seed_dispersal_fraction",
    "fates_trim_limit",
    "fates_trim_inc",
    "fates_turb_leaf_diameter",
    "fates_turb_z0mr",
    "fates_turb_displar",
    "fates_phen_flush_fraction",
    "fates_phen_cold_size_threshold",
    "fates_cnp_eca_decompmicc",
    "fates_cnp_eca_km_nh4",
    "fates_cnp_vmax_nh4",
    "fates_cnp_eca_km_no3",
    "fates_cnp_vmax_no3",
    "fates_cnp_eca_km_p",
    "fates_cnp_vmax_p",
    "fates_cnp_eca_km_ptase",
    "fates_cnp_eca_vmax_ptase",
    "fates_cnp_eca_alpha_ptase",
    "fates_cnp_eca_lambda_ptase",
    "fates_cnp_prescribed_nuptake",
    "fates_cnp_prescribed_puptake",
    "fates_landuse_harvest_pprod10",
    "fates_landuse_luc_frac_burned",
    "fates_landuse_luc_frac_exported",
    "fates_landuse_luc_pprod10",
    "fates_dev_arbitrary_pft",
)

# numrad parameters (stored on the param file as separate 1-D arrays, assembled
# into 2-D rhol/rhos/taul/taus in Receive).
const _edpft_numrad_param_names = (
    "fates_rad_leaf_rhovis",
    "fates_rad_leaf_rhonir",
    "fates_rad_stem_rhovis",
    "fates_rad_stem_rhonir",
    "fates_rad_leaf_tauvis",
    "fates_rad_leaf_taunir",
    "fates_rad_stem_tauvis",
    "fates_rad_stem_taunir",
)

# PFT x organ (hydraulics) 2-D parameter names.
const _edpft_hydr_organ_param_names = (
    "fates_hydro_vg_alpha_node",
    "fates_hydro_vg_m_node",
    "fates_hydro_vg_n_node",
    "fates_hydro_avuln_node",
    "fates_hydro_p50_node",
    "fates_hydro_thetas_node",
    "fates_hydro_epsil_node",
    "fates_hydro_pitlp_node",
    "fates_hydro_resid_node",
    "fates_hydro_fcap_node",
    "fates_hydro_pinot_node",
    "fates_hydro_kmax_node",
)

function Register_PFT!(this::EDPftvarcon_type, fates_params::fates_parameters_type)
    dim_names = [dimension_name_pft]
    dim_lower_bound = [lower_bound_pft]

    for nm in _edpft_pft_param_names
        RegisterParameter!(fates_params, nm, dimension_shape_1d, dim_names;
                           lower_bounds=dim_lower_bound)
    end

    # 2-D HLM-PFT map (FATES PFTno x HLM PFTno)
    pftmap_dim_names = [dimension_name_pft, dimension_name_hlm_pftno]
    RegisterParameter!(fates_params, "fates_hlm_pft_map", dimension_shape_2d,
                       pftmap_dim_names; lower_bounds=dim_lower_bound)
    return nothing
end

function Receive_PFT!(this::EDPftvarcon_type, fates_params::fates_parameters_type)
    R(nm) = RetrieveParameter1DAllocate(fates_params, nm)

    this.freezetol            = R("fates_mort_freezetol")
    this.hgt_min              = R("fates_recruit_height_min")
    this.bark_scaler          = R("fates_fire_bark_scaler")
    this.crown_kill           = R("fates_fire_crown_kill")
    this.initd                = R("fates_recruit_init_density")
    this.seed_suppl           = R("fates_recruit_seed_supplement")
    this.bb_slope             = R("fates_leaf_stomatal_slope_ballberry")
    this.medlyn_slope         = R("fates_leaf_stomatal_slope_medlyn")
    this.stomatal_intercept   = R("fates_leaf_stomatal_intercept")
    this.lf_flab              = R("fates_frag_leaf_flab")
    this.lf_fcel              = R("fates_frag_leaf_fcel")
    this.lf_flig              = R("fates_frag_leaf_flig")
    this.fr_flab              = R("fates_frag_fnrt_flab")
    this.fr_fcel              = R("fates_frag_fnrt_fcel")
    this.fr_flig              = R("fates_frag_fnrt_flig")
    this.xl                   = R("fates_rad_leaf_xl")
    this.clumping_index       = R("fates_rad_leaf_clumping_index")
    this.c3psn                = R("fates_leaf_c3psn")
    this.smpso                = R("fates_nonhydro_smpso")
    this.smpsc                = R("fates_nonhydro_smpsc")
    this.maintresp_reduction_curvature     = R("fates_maintresp_reduction_curvature")
    this.maintresp_reduction_intercept     = R("fates_maintresp_reduction_intercept")
    this.maintresp_reduction_upthresh      = R("fates_maintresp_reduction_upthresh")
    this.maintresp_leaf_atkin2017_baserate = R("fates_maintresp_leaf_atkin2017_baserate")
    this.maintresp_leaf_ryan1991_baserate  = R("fates_maintresp_leaf_ryan1991_baserate")
    this.maintresp_leaf_vert_scaler_coeff1 = R("fates_maintresp_leaf_vert_scaler_coeff1")
    this.maintresp_leaf_vert_scaler_coeff2 = R("fates_maintresp_leaf_vert_scaler_coeff2")
    this.prescribed_npp_canopy           = R("fates_prescribed_npp_canopy")
    this.prescribed_npp_understory       = R("fates_prescribed_npp_understory")
    this.prescribed_mortality_canopy     = R("fates_mort_prescribed_canopy")
    this.prescribed_mortality_understory = R("fates_mort_prescribed_understory")
    this.prescribed_recruitment          = R("fates_recruit_prescribed_rate")
    this.damage_frac            = R("fates_damage_frac")
    this.damage_mort_p1         = R("fates_damage_mort_p1")
    this.damage_mort_p2         = R("fates_damage_mort_p2")
    this.damage_recovery_scalar = R("fates_damage_recovery_scalar")
    this.fire_alpha_SH          = R("fates_fire_alpha_SH")
    this.allom_frbstor_repro    = R("fates_allom_frbstor_repro")
    this.hydr_p_taper           = R("fates_hydro_p_taper")
    this.hydr_rs2               = R("fates_hydro_rs2")
    this.hydr_srl               = R("fates_hydro_srl")
    this.hydr_rfrac_stem        = R("fates_hydro_rfrac_stem")
    this.hydr_avuln_gs          = R("fates_hydro_avuln_gs")
    this.hydr_p50_gs            = R("fates_hydro_p50_gs")
    this.hydr_k_lwp             = R("fates_hydro_k_lwp")
    this.bmort                  = R("fates_mort_bmort")
    this.mort_scalar_coldstress = R("fates_mort_scalar_coldstress")
    this.mort_scalar_cstarvation = R("fates_mort_scalar_cstarvation")
    this.mort_scalar_hydrfailure = R("fates_mort_scalar_hydrfailure")
    this.mort_upthresh_cstarvation = R("fates_mort_upthresh_cstarvation")
    this.mort_ip_size_senescence = R("fates_mort_ip_size_senescence")
    this.mort_r_size_senescence  = R("fates_mort_r_size_senescence")
    this.mort_ip_age_senescence  = R("fates_mort_ip_age_senescence")
    this.mort_r_age_senescence   = R("fates_mort_r_age_senescence")
    this.hf_sm_threshold        = R("fates_mort_hf_sm_threshold")
    this.hf_flc_threshold       = R("fates_mort_hf_flc_threshold")
    this.vcmaxha                = R("fates_leaf_vcmaxha")
    this.jmaxha                 = R("fates_leaf_jmaxha")
    this.vcmaxhd                = R("fates_leaf_vcmaxhd")
    this.jmaxhd                 = R("fates_leaf_jmaxhd")
    this.vcmaxse                = R("fates_leaf_vcmaxse")
    this.jmaxse                 = R("fates_leaf_jmaxse")
    this.germination_rate       = R("fates_recruit_seed_germination_rate")
    this.repro_frac_seed        = R("fates_trs_repro_frac_seed")
    this.a_emerg                = R("fates_trs_seedling_a_emerg")
    this.b_emerg                = R("fates_trs_seedling_b_emerg")
    this.par_crit_germ          = R("fates_trs_seedling_par_crit_germ")
    this.seedling_psi_emerg     = R("fates_trs_seedling_psi_emerg")
    this.seedling_psi_crit      = R("fates_trs_seedling_psi_crit")
    this.seedling_light_rec_a   = R("fates_trs_seedling_light_rec_a")
    this.seedling_light_rec_b   = R("fates_trs_seedling_light_rec_b")
    this.seedling_mdd_crit      = R("fates_trs_seedling_mdd_crit")
    this.seedling_h2o_mort_a    = R("fates_trs_seedling_h2o_mort_a")
    this.seedling_h2o_mort_b    = R("fates_trs_seedling_h2o_mort_b")
    this.seedling_h2o_mort_c    = R("fates_trs_seedling_h2o_mort_c")
    this.seedling_root_depth    = R("fates_trs_seedling_root_depth")
    this.seedling_light_mort_a  = R("fates_trs_seedling_light_mort_a")
    this.seedling_light_mort_b  = R("fates_trs_seedling_light_mort_b")
    this.background_seedling_mort = R("fates_trs_seedling_background_mort")
    this.seed_decay_rate        = R("fates_frag_seed_decay_rate")
    this.seed_dispersal_pdf_scale = R("fates_seed_dispersal_pdf_scale")
    this.seed_dispersal_pdf_shape = R("fates_seed_dispersal_pdf_shape")
    this.seed_dispersal_max_dist  = R("fates_seed_dispersal_max_dist")
    this.seed_dispersal_fraction  = R("fates_seed_dispersal_fraction")
    this.trim_limit             = R("fates_trim_limit")
    this.trim_inc               = R("fates_trim_inc")
    this.dleaf                  = R("fates_turb_leaf_diameter")
    this.z0mr                   = R("fates_turb_z0mr")
    this.displar                = R("fates_turb_displar")
    this.phenflush_fraction     = R("fates_phen_flush_fraction")
    this.phen_cold_size_threshold = R("fates_phen_cold_size_threshold")
    this.prescribed_nuptake     = R("fates_cnp_prescribed_nuptake")
    this.prescribed_puptake     = R("fates_cnp_prescribed_puptake")
    this.dev_arbitrary_pft      = R("fates_dev_arbitrary_pft")
    this.decompmicc             = R("fates_cnp_eca_decompmicc")
    this.eca_km_nh4             = R("fates_cnp_eca_km_nh4")
    this.vmax_nh4               = R("fates_cnp_vmax_nh4")
    this.eca_km_no3             = R("fates_cnp_eca_km_no3")
    this.vmax_no3               = R("fates_cnp_vmax_no3")
    this.eca_km_p               = R("fates_cnp_eca_km_p")
    this.vmax_p                 = R("fates_cnp_vmax_p")
    this.eca_km_ptase           = R("fates_cnp_eca_km_ptase")
    this.eca_vmax_ptase         = R("fates_cnp_eca_vmax_ptase")
    this.eca_alpha_ptase        = R("fates_cnp_eca_alpha_ptase")
    this.eca_lambda_ptase       = R("fates_cnp_eca_lambda_ptase")
    this.harvest_pprod10            = R("fates_landuse_harvest_pprod10")
    this.landusechange_frac_burned   = R("fates_landuse_luc_frac_burned")
    this.landusechange_frac_exported = R("fates_landuse_luc_frac_exported")
    this.landusechange_pprod10       = R("fates_landuse_luc_pprod10")

    # 2-D HLM-PFT map
    this.hlm_pft_map = RetrieveParameter2DAllocate(fates_params, "fates_hlm_pft_map")
    return nothing
end

# -------------------------------------------------------------------------------------
# numrad. These are stored on the param file as separate 1-D arrays (vis/nir),
# but FATES wants them as PFT x num_swb matrices with ivis=1, inir=2.

function Register_PFT_numrad!(this::EDPftvarcon_type, fates_params::fates_parameters_type)
    dim_names = [dimension_name_pft]
    for nm in _edpft_numrad_param_names
        RegisterParameter!(fates_params, nm, dimension_shape_1d, dim_names)
    end
    return nothing
end

function Receive_PFT_numrad!(this::EDPftvarcon_type, fates_params::fates_parameters_type)
    # All numrad params share the same PFT extent; size from a representative one.
    npft = size(FindParameterData(fates_params, "fates_rad_leaf_rhovis"), 1)

    this.rhol = Matrix{Float64}(undef, npft, num_swb)
    this.rhos = Matrix{Float64}(undef, npft, num_swb)
    this.taul = Matrix{Float64}(undef, npft, num_swb)
    this.taus = Matrix{Float64}(undef, npft, num_swb)

    this.rhol[:, ivis] = RetrieveParameter1DAllocate(fates_params, "fates_rad_leaf_rhovis")
    this.rhol[:, inir] = RetrieveParameter1DAllocate(fates_params, "fates_rad_leaf_rhonir")
    this.rhos[:, ivis] = RetrieveParameter1DAllocate(fates_params, "fates_rad_stem_rhovis")
    this.rhos[:, inir] = RetrieveParameter1DAllocate(fates_params, "fates_rad_stem_rhonir")
    this.taul[:, ivis] = RetrieveParameter1DAllocate(fates_params, "fates_rad_leaf_tauvis")
    this.taul[:, inir] = RetrieveParameter1DAllocate(fates_params, "fates_rad_leaf_taunir")
    this.taus[:, ivis] = RetrieveParameter1DAllocate(fates_params, "fates_rad_stem_tauvis")
    this.taus[:, inir] = RetrieveParameter1DAllocate(fates_params, "fates_rad_stem_taunir")
    return nothing
end

# Helper: return the stored data matrix for `name` (used to size numrad).
function FindParameterData(fates_params::fates_parameters_type, name::AbstractString)
    i = FindIndex(fates_params, name)
    return fates_params.parameters[i].data
end

# -------------------------------------------------------------------------------------
# leaf-age (PFT x leaf-age, 2-D)

function Register_PFT_leafage!(this::EDPftvarcon_type, fates_params::fates_parameters_type)
    dim_names = [dimension_name_pft, dimension_name_leaf_age]
    dim_lower_bound = [lower_bound_pft, lower_bound_general]
    RegisterParameter!(fates_params, "fates_leaf_vcmax25top", dimension_shape_2d,
                       dim_names; lower_bounds=dim_lower_bound)
    return nothing
end

function Receive_PFT_leafage!(this::EDPftvarcon_type, fates_params::fates_parameters_type)
    this.vcmax25top = RetrieveParameter2DAllocate(fates_params, "fates_leaf_vcmax25top")
    return nothing
end

# -------------------------------------------------------------------------------------
# hydraulic organs (PFT x organ, 2-D)

function Register_PFT_hydr_organs!(this::EDPftvarcon_type, fates_params::fates_parameters_type)
    dim_names = [dimension_name_pft, dimension_name_hydr_organs]
    dim_lower_bound = [lower_bound_pft, lower_bound_general]
    for nm in _edpft_hydr_organ_param_names
        RegisterParameter!(fates_params, nm, dimension_shape_2d, dim_names;
                           lower_bounds=dim_lower_bound)
    end
    return nothing
end

function Receive_PFT_hydr_organs!(this::EDPftvarcon_type, fates_params::fates_parameters_type)
    A(nm) = RetrieveParameter2DAllocate(fates_params, nm)
    this.hydr_vg_alpha_node = A("fates_hydro_vg_alpha_node")
    this.hydr_vg_m_node     = A("fates_hydro_vg_m_node")
    this.hydr_vg_n_node     = A("fates_hydro_vg_n_node")
    this.hydr_avuln_node    = A("fates_hydro_avuln_node")
    this.hydr_p50_node      = A("fates_hydro_p50_node")
    this.hydr_thetas_node   = A("fates_hydro_thetas_node")
    this.hydr_epsil_node    = A("fates_hydro_epsil_node")
    this.hydr_pitlp_node    = A("fates_hydro_pitlp_node")
    this.hydr_resid_node    = A("fates_hydro_resid_node")
    this.hydr_fcap_node     = A("fates_hydro_fcap_node")
    this.hydr_pinot_node    = A("fates_hydro_pinot_node")
    this.hydr_kmax_node     = A("fates_hydro_kmax_node")
    return nothing
end

# =====================================================================================

"""
    FatesReportPFTParams(is_master::Bool; debug_report::Bool=false)

Write the FATES PFT parameters to the log. As in Fortran, gated by both
`is_master` and a `debug_report` flag (default false), so this is a no-op unless
explicitly enabled. Operates on the live `EDPftvarcon_inst`.
"""
function FatesReportPFTParams(is_master::Bool; debug_report::Bool=false)
    (debug_report && is_master) || return nothing
    p = EDPftvarcon_inst[]
    npft = length(p.initd)
    if npft > 100
        println("you are trying to report pft parameters during initialization")
        println("but you have so many that it is over-running the format spec")
        fates_endrun("FatesReportPFTParams: too many PFTs for report format")
    end
    println("-----------  FATES PFT Parameters -----------------")
    println("freezetol = ", p.freezetol)
    println("hgt_min = ", p.hgt_min)
    println("dleaf = ", p.dleaf)
    println("z0mr = ", p.z0mr)
    println("displar = ", p.displar)
    println("bark_scaler = ", p.bark_scaler)
    println("crown_kill = ", p.crown_kill)
    println("initd = ", p.initd)
    println("c3psn = ", p.c3psn)
    println("vcmax25top = ", p.vcmax25top)
    println("phen_flush_fraction = ", p.phenflush_fraction)
    println("fire_alpha_SH = ", p.fire_alpha_SH)
    println("hlm_pft_map = ", p.hlm_pft_map)
    println("-------------------------------------------------")
    return nothing
end

# =====================================================================================

"""
    FatesCheckParams(is_master::Bool)

Cross-check user-supplied FATES PFT (and a few EDParamsMod) parameters for
logical consistency, raising a Julia `error` (Fortran `endrun`) on any failure.
Mirrors the Fortran `FatesCheckParams`. No-op (returns) unless `is_master`.
Operates on the live `EDPftvarcon_inst` and the live `ed_params()` / `prt_params`.
"""
function FatesCheckParams(is_master::Bool)
    p   = EDPftvarcon_inst[]
    edp = ed_params()

    npft = length(p.freezetol)

    is_master || return nothing

    # Radiation model must be Norman or two-stream.
    if !(edp.radiation_model in (norman_solver, twostr_solver))
        println("The only available canopy radiation models are Norman and Two-stream.")
        println("You specified fates_rad_model = ", edp.radiation_model)
        fates_endrun("FatesCheckParams: invalid fates_rad_model")
    end

    # Regeneration model must be a known type.
    if !(edp.regeneration_model in (default_regeneration, TRS_regeneration, TRS_no_seedling_dyn))
        println("The regeneration model must be set to a known model type (1, 2, or 3).")
        println("You specified fates_regeneration_model = ", edp.regeneration_model)
        fates_endrun("FatesCheckParams: invalid fates_regeneration_model")
    end

    # Day-length switch must be 0 or 1.
    if !(edp.dayl_switch in (itrue, ifalse))
        println("fates_daylength_factor_switch must be 0 or 1.")
        println("You specified fates_daylength_factor_switch = ", edp.dayl_switch)
        fates_endrun("FatesCheckParams: invalid fates_daylength_factor_switch")
    end

    # PARTEH-mode-dependent nutrient checks.
    parteh_mode = hlm_parteh_mode[]
    if parteh_mode == prt_cnp_flex_allom_hyp
        nu_com = strip(hlm_nu_com[])
        if !(nu_com == "RD" || nu_com == "ECA")
            println("FATES PARTEH CNP must have a valid BGC model: RD or ECA.")
            println("nu_comp: ", nu_com)
            fates_endrun("FatesCheckParams: CNP requires RD or ECA")
        end

        if hlm_nitrogen_spec[] > 0
            if nu_com == "ECA"
                if any(<(0.0), p.eca_km_nh4)
                    println("ECA with nitrogen is on; bad ECA km value(s) for nh4: ", p.eca_km_nh4)
                    fates_endrun("FatesCheckParams: bad eca_km_nh4")
                end
                if hlm_nitrogen_spec[] == 2
                    if any(<(0.0), p.eca_km_no3)
                        println("ECA with nit/denitr is on; bad ECA km value(s) for no3: ", p.eca_km_no3)
                        fates_endrun("FatesCheckParams: bad eca_km_no3")
                    end
                end
            end
        end

        # If any PFT has prescribed N uptake, all must.
        if any(x -> x < -nearzero, p.prescribed_nuptake) ||
           any(x -> x > 10.0, p.prescribed_nuptake)
            println("Negative (or >10) values for prescribed_nuptake are not allowed.")
            fates_endrun("FatesCheckParams: invalid prescribed_nuptake")
        elseif any(x -> abs(x) > nearzero, p.prescribed_nuptake)
            if !all(x -> abs(x) > nearzero, p.prescribed_nuptake)
                println("If any PFTs have prescribed N uptake, then they all must.")
                println("prescribed_nuptake: ", p.prescribed_nuptake)
                fates_endrun("FatesCheckParams: mixed prescribed_nuptake")
            end
        end

        # Simple phosphatase model: lambda and alpha must be 0.
        if hlm_phosphorus_spec[] > 0 && strip(hlm_nu_com[]) == "ECA"
            if any(x -> abs(x) > nearzero, p.eca_lambda_ptase)
                println("Critical Values for phosphatase in ECA not enabled; set fates_eca_lambda_ptase = 0.")
                fates_endrun("FatesCheckParams: eca_lambda_ptase must be 0")
            end
            if any(x -> abs(x) > nearzero, p.eca_alpha_ptase)
                println("No preferential P uptake from phosphatase; set fates_eca_alpha_ptase = 0.")
                fates_endrun("FatesCheckParams: eca_alpha_ptase must be 0")
            end
        end

    elseif parteh_mode == prt_carbon_allom_hyp
        # No additional checks needed.
    else
        println("FATES PARTEH supports only allometric carbon (1) and CNP (2).")
        println("fates_parteh_mode must be 1 or 2; got ", parteh_mode)
        fates_endrun("FatesCheckParams: invalid fates_parteh_mode")
    end

    # Logging fractions must sum to <= 1.
    if (edp.logging_mechanical_frac + edp.logging_collateral_frac + edp.logging_direct_frac) > 1.0
        println("logging_mechanical_frac + logging_collateral_frac + logging_direct_frac must be <= 1.")
        fates_endrun("FatesCheckParams: logging fractions sum > 1")
    end

    # Prescribed P uptake checks.
    if any(x -> x < -nearzero, p.prescribed_puptake) ||
       any(x -> x > 10.0, p.prescribed_puptake)
        println("Negative (or >10) values for prescribed_puptake are not allowed.")
        fates_endrun("FatesCheckParams: invalid prescribed_puptake")
    elseif any(x -> abs(x) > nearzero, p.prescribed_puptake)
        if !all(x -> abs(x) > nearzero, p.prescribed_puptake)
            println("If any PFTs have prescribed P uptake, then they all must.")
            println("prescribed_puptake: ", p.prescribed_puptake)
            fates_endrun("FatesCheckParams: mixed prescribed_puptake")
        end
    end

    # Per-PFT checks.
    for ipft in 1:npft
        # xl in [-0.4, 0.6]
        if p.xl[ipft] < -0.4 || p.xl[ipft] > 0.6
            println("fates_rad_leaf_xl for pft ", ipft, " is outside [-0.4, 0.6]")
            fates_endrun("FatesCheckParams: xl out of range")
        end

        # Seed-dispersal: if one of the four is set, all must be.
        if p.seed_dispersal_pdf_scale[ipft] < fates_check_param_set &&
           (p.seed_dispersal_max_dist[ipft]  > fates_check_param_set ||
            p.seed_dispersal_pdf_shape[ipft] > fates_check_param_set ||
            p.seed_dispersal_fraction[ipft]  > fates_check_param_set)
            println("Seed dispersal on (pdf_scale set); provide all seed_dispersal parameters.")
            fates_endrun("FatesCheckParams: incomplete seed_dispersal params (pdf_scale)")
        end
        if p.seed_dispersal_pdf_shape[ipft] < fates_check_param_set &&
           (p.seed_dispersal_max_dist[ipft]  > fates_check_param_set ||
            p.seed_dispersal_pdf_scale[ipft] > fates_check_param_set ||
            p.seed_dispersal_fraction[ipft]  > fates_check_param_set)
            println("Seed dispersal on (pdf_shape set); provide all seed_dispersal parameters.")
            fates_endrun("FatesCheckParams: incomplete seed_dispersal params (pdf_shape)")
        end
        if p.seed_dispersal_max_dist[ipft] < fates_check_param_set &&
           (p.seed_dispersal_pdf_shape[ipft] > fates_check_param_set ||
            p.seed_dispersal_pdf_scale[ipft] > fates_check_param_set ||
            p.seed_dispersal_fraction[ipft]  > fates_check_param_set)
            println("Seed dispersal on (max_dist set); provide all seed_dispersal parameters.")
            fates_endrun("FatesCheckParams: incomplete seed_dispersal params (max_dist)")
        end
        if p.seed_dispersal_fraction[ipft] < fates_check_param_set &&
           (p.seed_dispersal_pdf_shape[ipft] > fates_check_param_set ||
            p.seed_dispersal_pdf_scale[ipft] > fates_check_param_set ||
            p.seed_dispersal_max_dist[ipft]  > fates_check_param_set)
            println("Seed dispersal on (fraction set); provide all seed_dispersal parameters.")
            fates_endrun("FatesCheckParams: incomplete seed_dispersal params (fraction)")
        end

        # Seed-dispersal fraction (when set) must be within [0, 1].
        if p.seed_dispersal_fraction[ipft] < fates_check_param_set &&
           (p.seed_dispersal_fraction[ipft] > 1.0 || p.seed_dispersal_fraction[ipft] < 0.0)
            println("Seed dispersal fraction must be between 0 and 1.")
            fates_endrun("FatesCheckParams: seed_dispersal_fraction out of [0,1]")
        end

        # Age-dependent mortality: if ip set, r must be set.
        if p.mort_ip_age_senescence[ipft] < fates_check_param_set &&
           p.mort_r_age_senescence[ipft] > fates_check_param_set
            println("Age-dependent mortality is on; please also set mort_r_age_senescence.")
            fates_endrun("FatesCheckParams: missing mort_r_age_senescence")
        end

        # Size-dependent mortality: if ip set, r must be set.
        if p.mort_ip_size_senescence[ipft] < fates_check_param_set &&
           p.mort_r_size_senescence[ipft] > fates_check_param_set
            println("Size-dependent mortality is on; please also set mort_r_size_senescence.")
            fates_endrun("FatesCheckParams: missing mort_r_size_senescence")
        end

        # Size-dependent mortality parameters must be non-negative.
        if p.mort_ip_size_senescence[ipft] < 0.0 || p.mort_r_size_senescence[ipft] < 0.0
            println("mort_ip_size_senescence / mort_r_size_senescence cannot be negative.")
            fates_endrun("FatesCheckParams: negative size-senescence params")
        end

        # Deciduous plants must flush 0 < phenflush_fraction <= 1 of storage.
        if prt_params.evergreen[ipft] == ifalse
            if p.phenflush_fraction[ipft] < nearzero || p.phenflush_fraction[ipft] > 1
                println("Deciduous plants must flush some storage carbon on bud-burst.")
                println(" PFT#: ", ipft, " phenflush_fraction: ", p.phenflush_fraction[ipft])
                fates_endrun("FatesCheckParams: invalid phenflush_fraction for deciduous PFT")
            end
        end

        # Freezing tolerance within reasonable bounds.
        if p.freezetol[ipft] > 60.0 || p.freezetol[ipft] < -273.1
            println("Freezing tolerance set to a strange value (degrees C).")
            println(" PFT#: ", ipft, " freezetol: ", p.freezetol[ipft])
            fates_endrun("FatesCheckParams: freezetol out of range")
        end

        # In a cold-start run initial density cannot be zero.
        if hlm_use_inventory_init[] == ifalse && abs(p.initd[ipft]) < nearzero
            println("In a cold start run initial density cannot be zero (PFT ", ipft, ").")
            fates_endrun("FatesCheckParams: zero initd in cold start")
        end

        # Grass sapwood allometry (smode 2) must not be a woody plant.
        if prt_params.allom_smode[ipft] == 2 && prt_params.woody[ipft] == itrue
            println("Allometry mode 2 is for grasses only; pft ", ipft, " is woody.")
            fates_endrun("FatesCheckParams: woody PFT with grass allom_smode")
        end

        # Fraction of storage to reproduction must be in [0, 1].
        if p.allom_frbstor_repro[ipft] < 0.0 || p.allom_frbstor_repro[ipft] > 1.0
            println("allom_frbstor_repro must be between 0 and 1 (PFT ", ipft, ").")
            fates_endrun("FatesCheckParams: allom_frbstor_repro out of [0,1]")
        end

        # Photosynthetic pathway must be C3 (1) or C4 (0).
        if p.c3psn[ipft] < 0.0 || p.c3psn[ipft] > 1.0
            println("c3psn must be 0 (C4) or 1 (C3); PFT ", ipft, " = ", p.c3psn[ipft])
            fates_endrun("FatesCheckParams: c3psn out of [0,1]")
        end

        # Fixed-biogeography: each HLM PFT's distribution into FATES PFTs sums to 1.
        if hlm_use_fixed_biogeog[] == itrue
            for hlm_pft in 1:size(p.hlm_pft_map, 2)
                sumarea = sum(@view p.hlm_pft_map[1:npft, hlm_pft])
                if abs(sumarea - 1.0) > nearzero
                    println("HLM PFT ", hlm_pft, " distribution into FATES PFTs does not sum to 1.")
                    println("Error is: ", sumarea - 1.0)
                    fates_endrun("FatesCheckParams: hlm_pft_map column does not sum to 1")
                end
            end
        end
    end

    # Nocomp: max nocomp PFTs per land use must be <= max patches per land use,
    # unless this is an SP run.
    if hlm_use_nocomp[] == itrue && hlm_use_sp[] == ifalse
        for i_lu in 1:n_landuse_cats
            if edp.max_nocomp_pfts_by_landuse[i_lu] > edp.maxpatches_by_landuse[i_lu]
                println("max nocomp PFTs must be <= number of patches, per land use type.")
                println("land use index: ", i_lu,
                        " max_nocomp_pfts: ", edp.max_nocomp_pfts_by_landuse[i_lu],
                        " maxpatches: ", edp.maxpatches_by_landuse[i_lu])
                fates_endrun("FatesCheckParams: nocomp PFTs exceed patches")
            end
        end
    end

    # Informational: temperature at which Atkin-2017 Rdark goes negative per PFT.
    # (Fortran writes to the log; no abort. Reproduced as a silent computation so
    #  the array accesses still validate.)
    for ipft in 1:npft
        r_0 = p.maintresp_leaf_atkin2017_baserate[ipft]
        lnc_top = prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[leaf_organ]]
        neg_lmr_temp = (-1.0 * (r_0 + lmr_r_1 * lnc_top)) / lmr_r_2
        # neg_lmr_temp is purely diagnostic; suppress unused warning.
        neg_lmr_temp === neg_lmr_temp
    end

    return nothing
end

# =====================================================================================

"""
    GetDecompyFrac(pft::Integer, organ_id::Integer, dcmpy::Integer) -> Float64

Match the correct litter-decomposability pool (labile/cellulose/lignin) with the
PFT parameter data for leaf or fine-root organs. Errors on an unknown pool or
organ index. Reads the live `EDPftvarcon_inst`.
"""
function GetDecompyFrac(pft::Integer, organ_id::Integer, dcmpy::Integer)
    p = EDPftvarcon_inst[]
    if organ_id == leaf_organ
        if dcmpy == ilabile
            return p.lf_flab[pft]
        elseif dcmpy == icellulose
            return p.lf_fcel[pft]
        elseif dcmpy == ilignin
            return p.lf_flig[pft]
        else
            fates_endrun("GetDecompyFrac: Unknown decompositibility pool index: $(dcmpy)")
        end
    elseif organ_id == fnrt_organ
        if dcmpy == ilabile
            return p.fr_flab[pft]
        elseif dcmpy == icellulose
            return p.fr_fcel[pft]
        elseif dcmpy == ilignin
            return p.fr_flig[pft]
        else
            fates_endrun("GetDecompyFrac: Unknown decompositibility pool index: $(dcmpy)")
        end
    else
        fates_endrun("GetDecompyFrac: Unknown parteh organ index: $(organ_id)")
    end
end
