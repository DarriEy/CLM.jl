# fates_driver_init.jl
# FATES live-driver wiring — W1 (init globals/alloc) + W2 (cold-start chain).
#
# This is the host-side bootstrap that mirrors the Fortran sequence
#   CLMFatesGlobals1  -> set_fates_ctrlparms(...) + SetFatesGlobalElements1
#   CLMFatesGlobals2  -> SetFatesGlobalElements2
#   clm_fates%Init    -> allocate bc_in/bc_out/bc_pconst + zero_bcs
#   init_coldstart    -> init_site_vars! / zero_site! / set_site_properties! /
#                        init_patches! (-> init_cohorts!) / ed_update_site
# (src/utils/clmfates_interfaceMod.F90).
#
# It builds a finite, cold-started, carbon-only single FATES site and attaches it
# to a CLMInstances via the `fates::Union{fates_interface_type,Nothing}` field.
# That field is excluded from the AD dual-copy (_make_dual_instances) and from the
# GPU @adapt — FATES columns are CPU-only / non-differentiable by design.
#
# Minimal config (the simplest branch that cold-starts with NO external data):
#   * use_fates = true, carbon-only PARTEH (prt_carbon_allom_hyp)
#   * planthydro OFF by default (use_planthydro=true runs InitHydroGlobals /
#     InitHydrSites / HydrSiteColdStart), SPITFIRE/fire OFF (no_fire), inventory-init OFF
#   * SP OFF, nocomp OFF, fixed-biogeog OFF, LUH OFF, harvest/logging OFF
#   * tree-damage OFF, cohort-age-tracking OFF, ch4 OFF, vertsoilc OFF (=> nlevdecomp 1)
# This is the default NBG (natural-bare-ground) cold-start: init_patches! builds a
# single primaryland patch holding the whole site area, init_cohorts! seeds one
# small cohort per PFT.  It needs none of the SP-LAI / inventory / biogeog inputs.
#
# SUPERSEDED: the PFT / allometry / ED parameter tables are now read from the REAL
# FATES default parameter file (data/fates/fates_params_default.cdl) via
# read_fates_params! (fates_params_reader.jl), which clm_fates_init! calls below.
# The synthetic in-memory table _fates_spike_setup_pft! (a 2-PFT woody-evergreen
# carbon table, the EDInitMod unit-test table) is kept below for reference but is
# NO LONGER CALLED by the live cold start.

"""
    _fates_spike_setup_pft!(npft) -> npft

SUPERSEDED by [`read_fates_params!`](@ref) (real FATES param-file reader). Kept as
a reference for the minimal synthetic table only — NOT called by `clm_fates_init!`.

Populate the in-memory FATES parameter tables (prt_params allometry, param_derived,
EDPftvarcon, ed_params bins/tols, sf_params CWD frac) for a minimal `npft`-PFT
woody-evergreen carbon-only cold start, and set the size/age/coage/damage level
counts.  Mirrors the synthetic table the EDInitMod test builds.
"""
function _fates_spike_setup_pft!(npft::Int)
    p = prt_params
    allocate_prt_params!(p, npft, num_organ_types, 1)  # 1 leaf age class

    p.c2b          .= 2.0
    p.wood_density .= 0.6
    p.slatop       .= 0.012
    p.slamax       .= 0.020
    p.allom_agb_frac .= 0.6
    p.allom_dbh_maxheight .= 90.0
    p.allom_la_per_sa_int .= 0.8e3
    p.allom_la_per_sa_slp .= 0.0
    p.allom_sai_scaler    .= 0.1
    p.allom_l2fr          .= 1.0
    p.cushion             .= 1.0
    p.allom_blca_expnt_diff      .= 0.0
    p.allom_d2ca_coefficient_min .= 0.3
    p.allom_d2ca_coefficient_max .= 0.6
    p.allom_h2cd1 .= 0.5
    p.allom_h2cd2 .= 1.0

    p.allom_dmode .= 1
    p.woody       .= 1
    p.allom_cmode .= 1
    p.allom_smode .= 1
    p.allom_fmode .= 1
    p.allom_stmode .= 1
    p.allom_hmode .= 1
    p.allom_amode .= 1
    p.allom_lmode .= 1

    p.allom_d2h1 .= 0.64
    p.allom_d2h2 .= 0.37
    p.allom_d2h3 .= -0.034

    p.allom_agb1 .= 0.06896; p.allom_agb2 .= 0.572; p.allom_agb3 .= 1.94; p.allom_agb4 .= 0.931
    p.allom_d2bl1 .= 0.07; p.allom_d2bl2 .= 1.3; p.allom_d2bl3 .= 0.55

    p.fnrt_prof_mode .= 1
    p.fnrt_prof_a    .= 0.976
    p.fnrt_prof_b    .= 0.0

    # Evergreen: the "fully flushed" cohort branch.
    p.stress_decid .= 0
    p.season_decid .= 0
    p.evergreen    .= itrue
    p.phen_fnrt_drop_fraction .= 0.0
    p.phen_stem_drop_fraction .= 0.0
    p.leaf_stor_priority .= 0.8
    p.leaf_long          .= 1.5
    p.leaf_long_ustory   .= 1.5

    p.seed_alloc         .= 0.1
    p.seed_alloc_mature  .= 0.0
    p.dbh_repro_threshold .= 1000.0
    p.repro_alloc_a      .= 0.0
    p.repro_alloc_b      .= 0.0

    # PRT N stoichiometry (PFT x organ) + the organ->param-column reverse lookup.
    # Read by the maintenance-respiration path in FatesPlantRespPhotosynthDrive
    # (live_stem_n/live_croot_n/fnrt_n) and by EDInitMod cold-start mass init.
    # Identity organ_param_id (organ index == nitr_stoich_p1 column) + nominal
    # broadleaf N:C fractions per organ (leaf/fnrt/sapw/store/repro/struct).
    p.organ_param_id .= collect(1:num_organ_types)
    p.nitr_stoich_p1[:, leaf_organ]   .= 0.040
    p.nitr_stoich_p1[:, fnrt_organ]   .= 0.030
    p.nitr_stoich_p1[:, sapw_organ]   .= 0.0047
    p.nitr_stoich_p1[:, store_organ]  .= 0.040
    p.nitr_stoich_p1[:, repro_organ]  .= 0.024
    p.nitr_stoich_p1[:, struct_organ] .= 0.0047

    pd = param_derived_type()
    pd.branch_frac = fill(0.25, npft)
    pd.jmax25top = fill(85.0, npft, 1)
    pd.tpu25top  = fill(8.0, npft, 1)
    pd.kp25top   = fill(0.6, npft, 1)
    ParamDerived[] = pd

    ev = EDPftvarcon_type()
    ev.vcmax25top = fill(50.0, npft, 1)
    ev.damage_recovery_scalar = fill(0.0, npft)
    ev.lf_flab = fill(0.25, npft); ev.lf_fcel = fill(0.50, npft); ev.lf_flig = fill(0.25, npft)
    ev.fr_flab = fill(0.25, npft); ev.fr_fcel = fill(0.50, npft); ev.fr_flig = fill(0.25, npft)
    ev.allom_frbstor_repro = fill(0.0, npft)
    ev.landusechange_frac_burned   = fill(0.0, npft)
    ev.landusechange_frac_exported = fill(0.0, npft)
    # ECA nutrient-competition vectors — unused in carbon-only, but set_bcpconst
    # copies them, so allocate as zeros to keep the (inert) copy in-bounds.
    ev.vmax_nh4 = zeros(npft); ev.vmax_no3 = zeros(npft); ev.vmax_p = zeros(npft)
    ev.eca_km_nh4 = zeros(npft); ev.eca_km_no3 = zeros(npft); ev.eca_km_p = zeros(npft)
    ev.eca_km_ptase = zeros(npft); ev.eca_vmax_ptase = zeros(npft)
    ev.eca_alpha_ptase = zeros(npft); ev.eca_lambda_ptase = zeros(npft)
    # Cold-start init knobs read by init_cohorts!:
    ev.initd   = fill(0.2, npft)   # positive => initial recruit density
    ev.hgt_min = fill(1.25, npft)  # sapling height (m)
    ev.hlm_pft_map = [Float64(i == j) for i in 1:npft, j in 1:npft]  # identity HLM->FATES map

    # ---- Radiation optical params (read by PatchNormanRadiation) -------------
    # Leaf/stem reflectance + transmittance (npft, num_swb), orientation index xl,
    # and self-occlusion clumping. Representative broadleaf values.
    ev.rhol = fill(0.10, npft, num_swb)   # leaf reflectance
    ev.rhos = fill(0.16, npft, num_swb)   # stem reflectance
    ev.taul = fill(0.05, npft, num_swb)   # leaf transmittance
    ev.taus = fill(0.001, npft, num_swb)  # stem transmittance
    ev.xl   = fill(0.01, npft)            # leaf/stem orientation index
    ev.clumping_index = fill(0.85, npft)  # canopy clumping

    # ---- BTRAN params (read by btran_ed!) ------------------------------------
    ev.smpso = fill(-66000.0, npft)   # soil potential at full stomatal open  [mm]
    ev.smpsc = fill(-255000.0, npft)  # soil potential at full stomatal close [mm]

    # ---- Photosynthesis params (read by FatesPlantRespPhotosynthDrive) -------
    ev.c3psn              = fill(1.0, npft)        # C3 pathway
    ev.bb_slope           = fill(9.0, npft)        # Ball-Berry slope
    ev.medlyn_slope       = fill(4.1, npft)        # Medlyn slope
    ev.stomatal_intercept = fill(10000.0, npft)    # C3 Ball-Berry intercept
    ev.maintresp_leaf_ryan1991_baserate  = fill(2.525e-6, npft)
    ev.maintresp_leaf_atkin2017_baserate = fill(1.756, npft)
    ev.maintresp_leaf_vert_scaler_coeff1 = fill(1.5, npft)
    ev.maintresp_leaf_vert_scaler_coeff2 = fill(3.0, npft)
    ev.maintresp_reduction_curvature  = fill(0.01, npft)
    ev.maintresp_reduction_intercept  = fill(1.0, npft)
    ev.vcmaxha = fill(65330.0, npft);  ev.jmaxha = fill(43540.0, npft)
    ev.vcmaxhd = fill(149250.0, npft); ev.jmaxhd = fill(152040.0, npft)
    ev.vcmaxse = fill(485.0, npft);    ev.jmaxse = fill(495.0, npft)
    ev.hydr_k_lwp = fill(0.0, npft)    # plant hydraulics off

    # ---- Mortality params (read by mortality_rates in the daily demographic
    # step, W5) -------------------------------------------------------------
    # Senescence OFF: setting mort_ip_*_senescence >= fates_check_param_set drives
    # the (smort/asmort) senescence branches to 0.  Background / carbon-starvation /
    # cold-stress / hydraulic-failure scalars set to representative small rates so
    # the demographic step has finite, physical mortality.
    ev.bmort                   = fill(0.014, npft)   # background mortality (/yr)
    ev.mort_ip_size_senescence = fill(fates_check_param_set, npft)  # senescence off
    ev.mort_r_size_senescence  = fill(0.0, npft)
    ev.mort_ip_age_senescence  = fill(fates_check_param_set, npft)  # senescence off
    ev.mort_r_age_senescence   = fill(0.0, npft)
    ev.hf_sm_threshold         = fill(1.0e-6, npft)  # btran hydraulic-failure thresh
    ev.hf_flc_threshold        = fill(0.5, npft)
    ev.mort_scalar_hydrfailure = fill(0.6, npft)
    ev.mort_scalar_cstarvation = fill(0.6, npft)
    ev.mort_upthresh_cstarvation = fill(0.1, npft)
    ev.mort_scalar_coldstress  = fill(3.0, npft)
    ev.freezetol               = fill(-80.0, npft)   # very cold-tolerant => no frmort

    # ---- Seed / recruitment params (read by SeedUpdate / recruitment in the
    # daily demographic step, W5) ------------------------------------------
    ev.seed_suppl              = fill(0.0, npft)   # no external seed supplement
    ev.seed_decay_rate         = fill(0.51, npft)  # seed pool -> litter decay (/yr)
    ev.germination_rate        = fill(0.5, npft)   # seed -> seedling germination
    ev.repro_frac_seed         = fill(0.5, npft)   # TRS only (unused: default regen)
    ev.seed_dispersal_fraction = fill(0.0, npft)   # seed dispersal off (cadence none)
    # Phenology / trimming knobs read on the integrate / recruitment path.
    ev.phenflush_fraction      = fill(0.5, npft)
    ev.phen_cold_size_threshold = fill(0.0, npft)
    ev.trim_inc                = fill(0.03, npft)
    ev.trim_limit              = fill(0.3, npft)

    EDPftvarcon_inst[] = ev

    edp = ed_params_type()
    edp.regeneration_model = default_regeneration
    # Norman (big-leaf) radiation solver — the simplest path that needs no
    # two-stream `twostr.band` allocation. Required for the W3 radiation driver
    # hooks (FatesNormalizedCanopyRadiation / FatesSunShadeFracs). Default param
    # is -9 ("unset"); without this the SunShadeFracs two-stream else-branch is
    # entered with an empty band vector.
    edp.radiation_model = norman_solver
    # Photosynthesis / maintenance-respiration model selectors + scalars read by
    # FatesPlantRespPhotosynthDrive (Ball-Berry stomata, net-assim, no acclim).
    edp.stomatal_model        = ballberry_model
    edp.stomatal_assim_model  = net_assim_model
    edp.photo_tempsens_model  = photosynth_acclim_model_none
    edp.dayl_switch           = 0
    edp.maintresp_leaf_model  = lmrmodel_ryan_1991
    edp.maintresp_nonleaf_baserate = 2.525e-6
    edp.q10_mr                = 1.5
    edp.theta_cj_c3           = 0.999
    edp.theta_cj_c4           = 0.999
    edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
    edp.ED_val_history_sizeclass_bin_edges  = [0.0, 5.0, 10.0, 20.0, 50.0]
    edp.ED_val_history_coageclass_bin_edges = [0.0, 3.0, 6.0, 9.0, 12.0]
    edp.ED_val_history_ageclass_bin_edges   = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    edp.ED_val_cohort_size_fusion_tol = 1.0e-6
    edp.ED_val_cohort_age_fusion_tol  = 1.0e-6
    edp.max_cohort_per_patch          = 100
    edp.ED_val_patch_fusion_tol       = 1.0e-6
    edp.fates_mortality_disturbance_fraction = 1.0
    edp.ED_val_understorey_death      = 0.55983
    edp.logging_coll_under_frac       = 0.55983
    # logging_event_code is read (Int(trunc(...))) at the TOP of IsItLoggingTime
    # — BEFORE the use_logging early-return — so it must be finite even with
    # logging OFF.  1.0 == "logging off" event code.  The daily demographic step
    # (W5) reaches IsItLoggingTime, unlike the W3/W4 per-step hooks.
    edp.logging_event_code            = 1.0
    edp.mort_cstarvation_model        = cstarvation_model_lin  # linear C-starvation
    edp.maxpatches_by_landuse         = fill(10, n_landuse_cats)
    edp.max_nocomp_pfts_by_landuse    = fill(npft, n_landuse_cats)
    edp.crop_lu_pft_vector            = fill(1, n_landuse_cats)
    edp.maxpatch_total                = sum(edp.maxpatches_by_landuse)
    nlv = nlevleaf
    edp.dinc_vai   = fill(1.0, nlv)
    edp.dlower_vai = collect(1.0:1.0:Float64(nlv))
    EDParams[] = edp

    nlevsclass[] = length(edp.ED_val_history_sizeclass_bin_edges)
    nlevcoage[]  = length(edp.ED_val_history_coageclass_bin_edges)
    nlevage[]    = length(edp.ED_val_history_ageclass_bin_edges)
    nlevdamage[] = length(edp.ED_val_history_damage_bin_edges)

    sfp = sf_params_type()
    sfp.SF_val_CWD_frac = [0.045, 0.075, 0.21, 0.67]
    SFParams[] = sfp

    return npft
end

"""
    _fates_spike_set_ctrlparms!(; numpft_in, nlevsoil, hist_dimlevel=1,
                                  spitfire_mode=0, parteh_mode=prt_carbon_allom_hyp,
                                  nitrogen_spec=0, phosphorus_spec=0,
                                  use_tree_damage=false, use_planthydro=false)

Set the FATES `hlm_*` control Refs for the carbon-only cold start (the defaults),
mirroring the `set_fates_ctrlparms` tag sequence in `CLMFatesGlobals1`/`CLMFatesGlobals2`
(see clmfates_interfaceMod.F90).  Brackets with `flush_to_unset` + `check_allset`.

Mode kwargs (all gated, defaulting to the carbon-only / off values so existing callers
are unchanged):
  * `spitfire_mode` — SPITFIRE fire mode (0 = no_fire; 1 = scalar_lightning;
    2 = successful_ignitions; 3 = anthro_ignitions). The `sf_*_def` enumeration Refs
    are set so `spitfire_mode > sf_nofire_def` (=0) enables the fire chain.
  * `parteh_mode` — PARTEH allocation hypothesis (`prt_carbon_allom_hyp` = carbon-only;
    `prt_cnp_flex_allom_hyp` = CNP). For CNP, pass `nitrogen_spec`/`phosphorus_spec`
    (the nutrient-acquisition specs; nonzero for an active N/P cycle).
  * `use_tree_damage` — enable the tree-damage disturbance class state.
  * `use_planthydro` — enable plant hydraulics (Tier B BC solve).
  * `use_ed_st3` — static stand structure (no growth / mortality / recruitment /
    disturbance; only the fast physiology runs).
  * `use_ed_prescribed_phys` — prescribed NPP / mortality / recruitment from the
    `fates_prescribed_*` parameters (incompatible with ST3).
  * `use_cohort_age_tracking` — track cohort age (coage) and use it in fusion.
  * `seeddisp_cadence` — inter-gridcell seed dispersal cadence (0 off, 1 daily,
    2 monthly, 3 yearly).
  * `use_inventory_init` / `inventory_ctrl_file` — initialize the sites from
    pss/css inventory files listed in the control file.
"""
function _fates_spike_set_ctrlparms!(; numpft_in::Int, nlevsoil::Int, hist_dimlevel::Int = 1,
                                       spitfire_mode::Int = 0,
                                       parteh_mode::Int = prt_carbon_allom_hyp,
                                       nitrogen_spec::Int = 0, phosphorus_spec::Int = 0,
                                       use_tree_damage::Bool = false,
                                       use_planthydro::Bool = false,
                                       use_fixed_biogeog::Bool = false,
                                       use_nocomp::Bool = false,
                                       use_sp::Bool = false,
                                       use_vertsoilc::Bool = false,
                                       decomp_method::AbstractString = "CTC",
                                       harvest_mode::AbstractString = fates_harvest_no_logging,
                                       use_luh::Bool = false,
                                       use_potentialveg::Bool = false,
                                       use_ed_st3::Bool = false,
                                       use_ed_prescribed_phys::Bool = false,
                                       use_cohort_age_tracking::Bool = false,
                                       seeddisp_cadence::Int = 0,
                                       use_inventory_init::Bool = false,
                                       inventory_ctrl_file::AbstractString = "none",
                                       do_check::Bool = true)
    set_fates_ctrlparms("flush_to_unset")

    # Biogeography / competition modes.  Fixed-biogeog is OFF by default (the NBG
    # cold start seeds every PFT); when the caller enables it the per-site
    # `use_this_pft` filter (fed from `bc_in.pft_areafrac`) restricts cold-start
    # seeding + runtime recruitment + seed rain to the climatically-present PFTs.
    set_fates_ctrlparms("use_fixed_biogeog"; ival = use_fixed_biogeog ? itrue : ifalse)
    # nocomp / SP mirror CLMFatesGlobals1 (pass_nocomp / pass_sp); FATES' check_allset
    # enforces SP => nocomp .and. fixed_biogeog.
    set_fates_ctrlparms("use_nocomp";        ival = use_nocomp ? itrue : ifalse)
    set_fates_ctrlparms("use_sp";            ival = use_sp ? itrue : ifalse)
    set_fates_ctrlparms("masterproc";        ival = itrue)

    # Radiation bands / soil dims / HLM identity.
    set_fates_ctrlparms("num_sw_bbands"; ival = num_swb)
    set_fates_ctrlparms("vis_sw_index";  ival = ivis)
    set_fates_ctrlparms("nir_sw_index";  ival = inir)
    set_fates_ctrlparms("num_lev_soil";  ival = nlevsoil)
    set_fates_ctrlparms("hlm_name";      cval = "CLM")
    set_fates_ctrlparms("hio_ignore_val"; rval = -9.9e30)
    set_fates_ctrlparms("soilwater_ipedof"; ival = 0)

    # PARTEH mode (carbon-only by default; prt_cnp_flex_allom_hyp for CNP) +
    # seed-dispersal cadence + history levels.
    set_fates_ctrlparms("parteh_mode";       ival = parteh_mode)
    set_fates_ctrlparms("seeddisp_cadence";  ival = seeddisp_cadence)
    set_fates_ctrlparms("hist_hifrq_dimlevel"; ival = hist_dimlevel)
    set_fates_ctrlparms("hist_dynam_dimlevel"; ival = hist_dimlevel)

    # Nutrient competition / decomposition.  Carbon-only => N/P spec off, CTC decomp.
    # For CNP, nitrogen_spec/phosphorus_spec select the (nonzero) acquisition specs.
    set_fates_ctrlparms("nu_com";        cval = "RD")
    # Host soil-decomposition method (clmfates_interfaceMod.F90:431-437 passes
    # MIMICS / CENTURY / NONE); the legacy no-soil-BGC configuration keeps "CTC".
    set_fates_ctrlparms("decomp_method"; cval = decomp_method)
    set_fates_ctrlparms("use_tree_damage"; ival = use_tree_damage ? itrue : ifalse)
    set_fates_ctrlparms("nitrogen_spec";   ival = nitrogen_spec)
    set_fates_ctrlparms("phosphorus_spec"; ival = phosphorus_spec)

    # SPITFIRE / fire.  spitfire_mode selects the mode (0 = no_fire); the four
    # definition Refs carry the HLM's mode enumeration, which CTSM takes from
    # FATESFireFactoryMod.F90:24-30 (clmfates_interfaceMod.F90:464-467):
    #   no_fire=0, scalar_lightning=1, lightning_from_data=2,
    #   successful_ignitions=3, anthro_ignitions=4, anthro_suppression=5.
    # (The port used to pass successful_ignitions=2 / anthro_ignitions=3, which made
    # mode 2 behave as successful ignitions (FDI forced to 1) and mode 3 as
    # anthropogenic ignitions.) spitfire_mode > sf_nofire_def enables the fire chain.
    set_fates_ctrlparms("spitfire_mode"; ival = spitfire_mode)
    set_fates_ctrlparms("sf_nofire_def"; ival = 0)
    set_fates_ctrlparms("sf_scalar_lightning_def";    ival = 1)
    set_fates_ctrlparms("sf_successful_ignitions_def"; ival = 3)
    set_fates_ctrlparms("sf_anthro_ignitions_def";    ival = 4)

    set_fates_ctrlparms("is_restart"; ival = ifalse)
    set_fates_ctrlparms("use_ch4";        ival = ifalse)
    # CTSM always passes use_vertsoilc=1 (clmfates_interfaceMod.F90:486-488) with
    # nlevdecomp = nlevsoi under FATES soil BGC; the legacy one-layer path passes 0.
    set_fates_ctrlparms("use_vertsoilc";  ival = use_vertsoilc ? itrue : ifalse)
    # clmfates_interfaceMod.F90:490-516 (pass_ed_st3 / pass_ed_prescribed_phys /
    # pass_cohort_age_tracking).
    set_fates_ctrlparms("use_ed_st3";     ival = use_ed_st3 ? itrue : ifalse)
    set_fates_ctrlparms("use_ed_prescribed_phys"; ival = use_ed_prescribed_phys ? itrue : ifalse)
    set_fates_ctrlparms("use_planthydro"; ival = use_planthydro ? itrue : ifalse)
    set_fates_ctrlparms("use_cohort_age_tracking"; ival = use_cohort_age_tracking ? itrue : ifalse)

    # Land-use harvest / logging / LUH (CLMFatesGlobals1, clmfates_interfaceMod.F90
    # 518-561). All off for the default `no_harvest` / no-LUH configuration.
    pass_logging, pass_lu_harvest, pass_num_lu_harvest_cats =
        fates_harvest_ctrl(harvest_mode)
    set_fates_ctrlparms("use_lu_harvest";     ival = pass_lu_harvest)
    set_fates_ctrlparms("num_lu_harvest_cats"; ival = pass_num_lu_harvest_cats)
    set_fates_ctrlparms("use_logging";        ival = pass_logging)
    set_fates_ctrlparms("use_luh2";           ival = use_luh ? itrue : ifalse)
    set_fates_ctrlparms("num_luh2_states";      ival = use_luh ? num_landuse_state_vars : 0)
    set_fates_ctrlparms("num_luh2_transitions"; ival = use_luh ? num_landuse_transition_vars : 0)
    set_fates_ctrlparms("use_fates_potentialveg"; ival = use_potentialveg ? itrue : ifalse)

    # Inventory init (clmfates_interfaceMod.F90:563-570). OFF by default, so the
    # cold-start NBG path runs. The control file must be "set" to a non-"unset" string
    # to pass check_allset even when it's unused.
    set_fates_ctrlparms("use_inventory_init"; ival = use_inventory_init ? itrue : ifalse)
    set_fates_ctrlparms("inventory_ctrl_file"; cval = inventory_ctrl_file)

    # The `check_allset` validation reads `edpftvarcon_inst().prescribed_*uptake`
    # for the CNP (parteh mode 2) "forced uptake" gate. Those PFT arrays are only
    # populated by `read_fates_params!`, which runs AFTER this in `clm_fates_init!`.
    # For carbon-only the gate is skipped (parteh_mode != 2), so an inline check is
    # fine; for CNP the caller passes `do_check=false` and runs `check_allset`
    # itself once params are loaded + prescribed uptake is set (mirroring the Fortran
    # `CLMFatesGlobals2` ordering, where check happens after the param read).
    do_check && set_fates_ctrlparms("check_allset")
    return nothing
end

# Cumulative NOLEAP days before each month (CTSM's default calendar).
const _FATES_NOLEAP_CUMDAYS = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
_fates_noleap_doy0(mon::Integer, day::Integer) = _FATES_NOLEAP_CUMDAYS[mon] + min(day, mon == 2 ? 28 : 31) - 1
_fates_yyyymmdd_split(d::Integer) = (d ÷ 10000, (d ÷ 100) % 100, d % 100)

"""
    fates_get_and_set_time!(year, mon, day, tod)

Set the FATES clock module-globals from the host date. Mirrors the CTSM host
routine `GetAndSetTime` (clmfates_interfaceMod.F90:3730-3786), which runs at FATES
init and at the top of every daily `dynamics_driv` (line 1068):

* `day_of_year = floor(timemgr_datediff(jan01_curr_year, now))` — 0-based;
* `model_day   = timemgr_datediff(reference_date, now)` [days, fractional];
* `days_per_year = get_curr_days_per_year()`, `freq_day = 1/days_per_year`,

all in the HOST calendar, [`FATES_HLM_CALENDAR`](@ref): `:noleap` (CTSM's default,
`calendar = NO_LEAP`) or `:gregorian` (`calendar = GREGORIAN` in the case's
nuopc.runconfig CLOCK_attributes; the fates_parity_1pt case runs GREGORIAN).

`hlm_reference_date[]` (YYYYMMDD) must be set first (`clm_fates_init!` sets it to
the run start date, CTSM's `ref_ymd` default). A `year <= 0` means "the caller does
not track the year": the year is then carried forward from the previous call and
advanced when the day of year wraps (`clm_drv!`'s `year` kwarg defaults to 0).
"""
function fates_get_and_set_time!(year::Integer, mon::Integer, day::Integer, tod::Integer)
    if year <= 0
        year = hlm_current_year[]
        if hlm_day_of_year[] != fates_unset_int &&
           _fates_noleap_doy0(mon, day) < _fates_noleap_doy0(hlm_current_month[], hlm_current_day[])
            year += 1
        end
    end
    ry, rm, rd = _fates_yyyymmdd_split(hlm_reference_date[])
    if FATES_HLM_CALENDAR[] === :gregorian
        now = Date(year, mon, day)
        doy = Dates.dayofyear(now) - 1
        dpy = Dates.daysinyear(year)
        model_day = Float64(Dates.value(now - Date(ry, rm, rd))) + tod / sec_per_day
    else
        doy = _fates_noleap_doy0(mon, day)
        dpy = 365
        model_day = 365.0 * (year - ry) + (doy - _fates_noleap_doy0(rm, rd)) + tod / sec_per_day
    end
    SetFatesTime(year, mon, day, tod, year * 10000 + mon * 100 + day,
                 hlm_reference_date[], model_day, doy, dpy, 1.0 / dpy)
    return nothing
end

"""
    FATES_HLM_CALENDAR :: Ref{Symbol}

The host model calendar the FATES clock follows: `:noleap` (default; CTSM
`NO_LEAP`) or `:gregorian`. Set it before `clm_initialize!` / `clm_fates_init!`
to match the case being reproduced.
"""
const FATES_HLM_CALENDAR = Ref(:noleap)

"""
    fates_check_niche_modes(; use_fates_tree_damage=false, use_fates_ed_st3=false,
                              use_fates_ed_prescribed_phys=false,
                              use_fates_cohort_age_tracking=false,
                              fates_seeddisp_cadence=0, use_fates_inventory_init=false,
                              fates_inventory_ctrl_filename="", fates_parteh_mode=1)

Host-side validity checks on the FATES research-mode namelist switches, mirroring
CTSM. FATES's own `check_allset` repeats the FATES-side rules after the parameter
read (FatesInterfaceMod.F90:1550-1590, 1644-1655); they are checked here too so a
bad combination fails before any state is built:

* `fates_seeddisp_cadence` ∈ {0,1,2,3} (namelist_definition_ctsm.xml valid_values);
* `use_fates_inventory_init` needs `fates_inventory_ctrl_filename` set and pointing
  at an existing file (CLMBuildNamelist.pm `setup_logic_fates`);
* ST3 and prescribed physiology cannot both be on (FatesInterfaceMod.F90:1570);
* inventory init cannot be combined with cohort age tracking (:1576);
* tree damage is not compatible with CNP allocation, parteh mode 2 (:1648).
"""
function fates_check_niche_modes(; use_fates_tree_damage::Bool = false,
                                   use_fates_ed_st3::Bool = false,
                                   use_fates_ed_prescribed_phys::Bool = false,
                                   use_fates_cohort_age_tracking::Bool = false,
                                   fates_seeddisp_cadence::Integer = 0,
                                   use_fates_inventory_init::Bool = false,
                                   fates_inventory_ctrl_filename::AbstractString = "",
                                   fates_parteh_mode::Integer = prt_carbon_allom_hyp)
    fates_seeddisp_cadence in 0:3 ||
        error("fates_seeddisp_cadence = $fates_seeddisp_cadence: valid values are " *
              "0 (off), 1 (daily), 2 (monthly), 3 (yearly)")
    if use_fates_inventory_init
        isempty(strip(fates_inventory_ctrl_filename)) &&
            error("fates_inventory_ctrl_filename is required when use_fates_inventory_init is set")
        isfile(fates_inventory_ctrl_filename) ||
            error("$fates_inventory_ctrl_filename does NOT point to a valid filename")
        use_fates_cohort_age_tracking &&
            error("Fates inventory init cannot be used with age dependent mortality: " *
                  "set use_fates_cohort_age_tracking to false or turn off inventory init")
    end
    (use_fates_ed_st3 && use_fates_ed_prescribed_phys) &&
        error("FATES ST3 and prescribed physiology cannot both be turned on")
    (use_fates_tree_damage && fates_parteh_mode == prt_cnp_flex_allom_hyp) &&
        error("FATES tree damage (use_fates_tree_damage) is not (yet) compatible with " *
              "CNP allocation (fates_parteh_mode = 2)")
    return nothing
end

"""
    clm_fates_init!(inst::CLMInstances; nsites=1, numpft_in=nothing, nlevsoil,
                    nlevdecomp, params_path=FATES_PARAMS_DEFAULT_CDL,
                    current_year=1, current_month=1, current_day=1)
        -> fates_interface_type

Bootstrap a finite, cold-started, carbon-only FATES interface with `nsites` sites
and attach it to `inst.fates`.  Mirrors the Fortran host sequence
`CLMFatesGlobals1 -> CLMFatesGlobals2 -> clm_fates%Init -> init_coldstart`.

Carbon-only minimal config by default (planthydro / SPITFIRE / inventory / LUH OFF; SP /
nocomp / fixed-biogeog OFF) — the NBG cold start that needs no external data. The
reduced-complexity modes are opt-in: `use_nocomp`, `use_sp` (requires `use_nocomp` and
`hlm_pft_areafrac`) and the host fixed biogeography `hlm_pft_areafrac` (per site, the
surfdata `wt_nat_patch` over HLM PFTs 0:surfpft_ub, mapped by `fates_hlm_pft_map`); SP
also needs the host satellite phenology to have run (see `clm_initialize!`).

Parameters are read from the REAL FATES default parameter file
(`data/fates/fates_params_default.cdl`, via [`read_fates_params!`](@ref)) — this
is the physically-meaningful replacement for the old synthetic 2-PFT table. The
PFT count (`numpft[]`) comes from the file (14 for the default file); pass
`numpft_in` only to override/validate it (a mismatch is an error). `params_path`
selects an alternative parameter file.

This sets the FATES module-global `hlm_*` control Refs and the parameter globals
as a side effect. The default test suite never calls this, so the globals are
left in the minimal-cold-start state on return; callers that mix FATES with other
FATES tests should save/restore the globals themselves.
"""
function clm_fates_init!(inst::CLMInstances; nsites::Int = 1,
                         numpft_in::Union{Int,Nothing} = nothing,
                         nlevsoil::Int, nlevdecomp::Int,
                         params_path::AbstractString = FATES_PARAMS_DEFAULT_CDL,
                         current_year::Int = 1, current_month::Int = 1,
                         current_day::Int = 1, hist_dimlevel::Int = 1,
                         spitfire_mode::Int = 0,
                         parteh_mode::Int = prt_carbon_allom_hyp,
                         nitrogen_spec::Int = 0, phosphorus_spec::Int = 0,
                         use_tree_damage::Bool = false,
                         use_planthydro::Bool = false,
                         fates_pft_areafrac::Union{Nothing,AbstractVector} = nothing,
                         fates_biogeog_screen::Symbol = :none,
                         use_sp::Bool = false,
                         use_nocomp::Bool = false,
                         hlm_pft_areafrac::Union{Nothing,AbstractVector} = nothing,
                         surf_numpft::Int = 0, surf_numcft::Int = 0,
                         use_vertsoilc::Bool = false,
                         decomp_method::AbstractString = "CTC",
                         harvest_mode::AbstractString = fates_harvest_no_logging,
                         use_luh::Bool = false,
                         use_potentialveg::Bool = false,
                         fluh_timeseries::AbstractString = "",
                         flandusepftdat::AbstractString = "",
                         use_ed_st3::Bool = false,
                         use_ed_prescribed_phys::Bool = false,
                         use_cohort_age_tracking::Bool = false,
                         seeddisp_cadence::Int = 0,
                         use_inventory_init::Bool = false,
                         inventory_ctrl_file::AbstractString = "none",
                         dtime::Real = 1800.0)

    # ---- Fixed-biogeography (climate-appropriate PFT screening) config ----------
    # Decide up front whether to run the cold start in fixed-biogeography mode. When
    # active, only a chosen subset of PFTs is seeded at the site (and recruited /
    # seed-rained thereafter), via the ported `use_this_pft` machinery. The subset
    # is either supplied explicitly (`fates_pft_areafrac`, a per-FATES-PFT area
    # vector) or derived by a named climate screen (`fates_biogeog_screen`):
    #   * :none               — OFF (default); every PFT seeded (byte-identical path)
    #   * :drop_cold_deciduous — present = PFTs with season_decid==false, i.e. a
    #                            never-cold (tropical/warm) site cannot host the
    #                            cold-deciduous PFTs whose forced leaf-off drives the
    #                            multi-year boom→bust die-back. Yields evergreen +
    #                            stress(drought)-deciduous PFTs.
    # Both feed `bc_in.pft_areafrac`; the switch must be known BEFORE the ctrlparms
    # set (it sizes the bc allocation) — its CONTENT is filled after the param read.
    if !(fates_biogeog_screen in (:none, :drop_cold_deciduous))
        error("clm_fates_init!: unknown fates_biogeog_screen=$(fates_biogeog_screen) " *
              "(expected :none or :drop_cold_deciduous)")
    end
    # Host (CTSM) fixed biogeography: `hlm_pft_areafrac[s]` is site s's surfdata
    # wt_nat_patch over HLM PFTs 0:surfpft_ub (element k+1 <-> HLM PFT k), the
    # clmfates_interfaceMod `bc_in%pft_areafrac` source. It maps onto FATES PFTs
    # through the parameter file's fates_hlm_pft_map (NOT the identity crosswalk of
    # the legacy FATES-PFT-space `fates_pft_areafrac` / screen path below).
    host_biogeog = hlm_pft_areafrac !== nothing
    host_biogeog && (fates_pft_areafrac !== nothing || fates_biogeog_screen != :none) &&
        error("clm_fates_init!: hlm_pft_areafrac cannot be combined with " *
              "fates_pft_areafrac / fates_biogeog_screen")
    use_sp && !(use_nocomp && host_biogeog) &&
        error("clm_fates_init!: use_sp requires use_nocomp and the host fixed " *
              "biogeography (hlm_pft_areafrac)")
    use_fixed_biogeog = host_biogeog || (fates_pft_areafrac !== nothing) ||
                        (fates_biogeog_screen != :none)
    # Land use (CTSM build-namelist / controlMod checks on the FATES land-use options).
    use_potentialveg && !use_luh &&
        error("clm_fates_init!: use_fates_luh must be true when use_fates_potentialveg is true")
    use_potentialveg && harvest_mode != fates_harvest_no_logging &&
        error("clm_fates_init!: fates_harvest_mode must be no_harvest when " *
              "use_fates_potentialveg is true")
    harvest_mode in (fates_harvest_luh_area, fates_harvest_luh_mass) && !use_luh &&
        error("clm_fates_init!: use_fates_luh is required to be true when " *
              "fates_harvest_mode is luhdata_mass or luhdata_area")
    use_luh && use_sp &&
        error("clm_fates_init!: use_fates_luh can NOT be true when use_fates_sp is true")
    use_luh && use_fixed_biogeog && !host_biogeog &&
        error("clm_fates_init!: use_fates_luh with fixed biogeography needs the host " *
              "biogeography (hlm_pft_areafrac + flandusepftdat)")

    # ---- W1.1: FATES globals + parameters (CLMFatesGlobals1) ----
    FatesInterfaceInit(6, false)

    # Control parameters (set_fates_ctrlparms tag sequence). MUST precede the param
    # read AND the element registry: it sets hlm_parteh_mode (carbon-only by default,
    # or CNP), which read_fates_params! branches on (PRTDerivedParams! organ map + the
    # PARTEH/PFT consistency checks) and which InitPARTEHGlobals dispatches on below.
    # nlevsoil is the only numeric arg it needs; numpft is resolved from the file.
    is_cnp = parteh_mode == prt_cnp_flex_allom_hyp
    _fates_spike_set_ctrlparms!(; numpft_in = 0, nlevsoil = nlevsoil,
                                  hist_dimlevel = hist_dimlevel,
                                  spitfire_mode = spitfire_mode,
                                  parteh_mode = parteh_mode,
                                  nitrogen_spec = nitrogen_spec,
                                  phosphorus_spec = phosphorus_spec,
                                  use_tree_damage = use_tree_damage,
                                  use_planthydro = use_planthydro,
                                  use_fixed_biogeog = use_fixed_biogeog,
                                  use_nocomp = use_nocomp,
                                  use_sp = use_sp,
                                  use_vertsoilc = use_vertsoilc,
                                  decomp_method = decomp_method,
                                  harvest_mode = harvest_mode,
                                  use_luh = use_luh,
                                  use_potentialveg = use_potentialveg,
                                  use_ed_st3 = use_ed_st3,
                                  use_ed_prescribed_phys = use_ed_prescribed_phys,
                                  use_cohort_age_tracking = use_cohort_age_tracking,
                                  seeddisp_cadence = seeddisp_cadence,
                                  use_inventory_init = use_inventory_init,
                                  inventory_ctrl_file = inventory_ctrl_file,
                                  # check_allset reads PFT parameters (the CNP forced-
                                  # uptake gate, the age-senescence/coage gate), so it
                                  # runs after the parameter read, as in CLMFatesGlobals2.
                                  do_check = false)

    # Read the REAL FATES default parameter file (replaces _fates_spike_setup_pft!).
    # This populates prt_params / EDParams / EDPftvarcon_inst / SFParams /
    # ParamDerived and the numpft[] / nlev*class[] / nleafage[] dim Refs. The N/P
    # stoichiometry (nitr_stoich_p1 / phos_stoich_p1) and the fire / damage params
    # are all read here regardless of mode (the CDL carries real values for them);
    # the modes only gate their RUNTIME use.
    read_fates_params!(; path = params_path)

    # ---- CNP "forced uptake" configuration (CLM only supports prescribed uptake) ----
    # The ctsm5.3 host never fills bc_in%plant_nh4/no3/p_uptake_flux (no FATES branch
    # in SoilBiogeochemCompetition, plant_ndemand_col = 0 for FATES columns,
    # CNDriverMod.F90:466-468), so check_allset endruns a CLM CNP run unless some
    # fates_cnp_prescribed_n/puptake > 0 (FatesInterfaceMod.F90:1728-1736). CTSM's
    # FatesColdPRT2 testmod sets both to 1.0 for all PFTs in the parameter file.
    # A parameter file that already prescribes uptake is used as-is (the uptake modes
    # and competitor scaling then follow from it in SetFatesGlobalElements2); the
    # shipped default file has both at 0, for which the FatesColdPRT2 values (1.0) are
    # applied here instead of aborting.
    if is_cnp
        evcon = edpftvarcon_inst()
        n = numpft[]
        if sum(abs, evcon.prescribed_nuptake) < nearzero &&
           sum(abs, evcon.prescribed_puptake) < nearzero
            evcon.prescribed_nuptake = fill(1.0, n)
            evcon.prescribed_puptake = fill(1.0, n)
        end
        # FatesInterfaceMod.F90:901-911 (SetFatesGlobalElements2): the uptake modes.
        edp_cnp = ed_params()
        edp_cnp.n_uptake_mode = any(x -> abs(x) > nearzero, evcon.prescribed_nuptake) ?
                                prescribed_n_uptake : coupled_n_uptake
        edp_cnp.p_uptake_mode = any(x -> abs(x) > nearzero, evcon.prescribed_puptake) ?
                                prescribed_p_uptake : coupled_p_uptake
    end
    # The deferred validation (params + any CNP prescribed uptake are in place).
    set_fates_ctrlparms("check_allset")
    numpft_resolved = numpft[]
    if numpft_in !== nothing && numpft_in != numpft_resolved
        error("clm_fates_init!: numpft_in=$numpft_in does not match the parameter " *
              "file's fates_pft dimension ($numpft_resolved). Omit numpft_in to " *
              "use the file value.")
    end

    # ---- Build the fixed-biogeog per-FATES-PFT area vector (post param read) -----
    # `season_decid` / `numpft` are only known once the param file is read.  The
    # vector is expressed directly in FATES-PFT space; to feed it through the ported
    # `_set_site_area_pft!` (which maps HLM-PFT areas onto FATES PFTs via
    # `hlm_pft_map`) unchanged, we set `hlm_pft_map` to the identity so
    # `area_PFT[ft] == pft_areafrac[ft]`.  This is the self-consistent "host already
    # supplies areas in FATES-PFT space" configuration and only runs on this gated
    # path (default path never touches the param-file crosswalk).
    fates_areafrac_vec = Float64[]
    if use_fixed_biogeog && !host_biogeog
        fates_areafrac_vec = zeros(Float64, numpft_resolved)
        if fates_pft_areafrac !== nothing
            length(fates_pft_areafrac) == numpft_resolved ||
                error("clm_fates_init!: fates_pft_areafrac has length " *
                      "$(length(fates_pft_areafrac)); expected numpft=$numpft_resolved.")
            any(<(0.0), fates_pft_areafrac) &&
                error("clm_fates_init!: fates_pft_areafrac has negative entries.")
            fates_areafrac_vec .= fates_pft_areafrac
        elseif fates_biogeog_screen == :drop_cold_deciduous
            # A never-cold site cannot host cold(season)-deciduous PFTs.  Give every
            # non-cold-deciduous PFT equal presence; normalization happens downstream.
            for ft in 1:numpft_resolved
                fates_areafrac_vec[ft] = (prt_params.season_decid[ft] == ifalse) ? 1.0 : 0.0
            end
        end
        sum(fates_areafrac_vec) > nearzero ||
            error("clm_fates_init!: fixed-biogeog PFT area vector is all zero " *
                  "(screen=$(fates_biogeog_screen)) — no PFT would be seeded.")
        # Identity HLM<->FATES crosswalk so pft_areafrac maps 1:1 onto area_PFT.
        edpftvarcon_inst().hlm_pft_map =
            [Float64(i == j) for i in 1:numpft_resolved, j in 1:numpft_resolved]
    end

    # PARTEH element registry + global hypothesis state. InitPARTEHGlobals dispatches
    # on hlm_parteh_mode (set above): carbon-only => num_elements=1 / carbon12 only /
    # InitPRTGlobalAllometricCarbon!; CNP => num_elements=3 / C,N,P / element_pos /
    # InitPRTGlobalAllometricCNP. The daily demographic step (W5) indexes
    # site.mass_balance by element_pos, so this MUST run before the cold start.
    InitPARTEHGlobals()
    nleafage[] = size(prt_params.leaf_long, 2)

    # SetFatesGlobalElements2 (FatesInterfaceMod.F90:398-400): under plant hydraulics,
    # build the plant water transfer functions (WRF/WKF per organ x PFT + stomata)
    # from the parameter file. Not needed (and not built) with hydraulics off.
    use_planthydro && InitHydroGlobals!(numpft[])

    # Running-mean time-averaging globals. Fortran FatesInterfaceMod:1041-1047 defines
    # every rmean_def with `up_period = hlm_stepsize` -- the HOST's actual timestep --
    # and derives `n_mem = nint(mem_period/up_period)`. Both `hlm_stepsize` and
    # `InitTimeAveragingGlobals()` were ported and then never called (the "ported then
    # never wired" class), so every patch fell back to the hardcoded `_patch_*`
    # definitions in FatesPatchMod.jl, which assume a 1800 s step. On the reference
    # FATES case (dtime = 3600 s) that makes the "24-hour" fixed window 48 samples
    # long = 48 HOURS: `tveg24` advances only every second day, `vegtemp_memory` fills
    # at half rate, and its unfilled 0.0 slots count as cold days (0 < phen_coldtemp),
    # firing a spurious cold leaf-off on all cold-deciduous PFTs in mid-July (D3).
    hlm_stepsize[] = Float64(dtime)
    InitTimeAveragingGlobals()

    # SetFatesGlobalElements1 sizes maxpatch_total / fates_maxPatchesPerSite from the
    # ed_params.maxpatches_by_landuse table we set above (non-SP / non-nocomp branch:
    # fates_maxPatchesPerSite = maxpatch_total + 1).
    p = ed_params()
    p.maxpatch_total = sum(p.maxpatches_by_landuse)
    if use_sp
        # SP: the primary-land patches hold one PFT each; the HLM needs enough patches
        # to carry the surfdata LAI of every surface PFT/CFT (surf_numpft incl. bare).
        p.maxpatches_by_landuse[primaryland] = numpft_resolved
        p.maxpatches_by_landuse[secondaryland:n_landuse_cats] .= 0
        p.maxpatch_total = numpft_resolved
        fates_maxPatchesPerSite[] = max(surf_numpft + surf_numcft, p.maxpatch_total + 1)
    else
        if use_nocomp
            p.maxpatches_by_landuse[primaryland] =
                max(p.maxpatches_by_landuse[primaryland], numpft_resolved)
            p.maxpatch_total = sum(p.maxpatches_by_landuse)
        end
        fates_maxPatchesPerSite[] = p.maxpatch_total + 1
    end

    # SetFatesGlobalElements2 derived dims the cold-start chain + bc allocators read.
    # The full SetFatesGlobalElements2 (InitHydroGlobals / fates_history_maps / VAI
    # bins) is not needed for the cold start, so set just these directly.
    # Restart-vector sizing (SetFatesGlobalElements2, non-SP branch); only the
    # restart interface reads these.
    fates_maxElementsPerPatch[] = max(num_swb, p.max_cohort_per_patch,
                                      ndcmpy * nlevsoil, ncwd * nlevsoil)
    fates_maxElementsPerSite[] = max(fates_maxPatchesPerSite[] * fates_maxElementsPerPatch[],
                                     numWaterMem * numpft[], num_vegtemp_mem, num_elements[],
                                     fates_maxPatchesPerSite[] * num_swb,
                                     nlevsclass[] * numpft[] * n_term_mort_types)
    # Nutrient competitors (FatesInterfaceMod.F90:913-926): one / trivial scaling
    # except CNP with a coupled (non-prescribed) N or P uptake mode.
    if is_cnp && (ed_params().n_uptake_mode == coupled_n_uptake ||
                  ed_params().p_uptake_mode == coupled_p_uptake)
        max_comp_per_site[]     = fates_maxElementsPerSite[]
        fates_np_comp_scaling[] = coupled_np_comp_scaling
    else
        max_comp_per_site[]     = 1
        fates_np_comp_scaling[] = trivial_np_comp_scaling
    end

    # ---- W1.2: FATES clock (SetFatesTime) ----
    # Fortran `GetAndSetTime` (clmfates_interfaceMod.F90:3730-3786), called here at
    # init and again at the top of every daily `dynamics_driv`. The reference date is
    # the run start (CTSM's ref_ymd defaults to start_ymd), so model_day starts at 0,
    # and day_of_year is the 0-based NOLEAP day (`floor(timemgr_datediff(jan01, now))`).
    # This used to hardcode day_of_year = 1 and reference_date = 10101 and was never
    # advanced (see `fates_get_and_set_time!`).
    hlm_reference_date[] = current_year * 10000 + current_month * 100 + current_day
    fates_get_and_set_time!(current_year, current_month, current_day, 0)

    # ---- W1.3: build the interface container + per-site boundary conditions ----
    fates = fates_interface_type()
    fates.nsites = nsites
    fates.sites  = [ed_site_type()  for _ in 1:nsites]
    fates.bc_in  = [bc_in_type()    for _ in 1:nsites]
    fates.bc_out = [bc_out_type()   for _ in 1:nsites]
    fates.bc_pconst = bc_pconst_type()

    # Parameter-constant block (ECA vectors + j_uptake map).
    allocate_bcpconst(fates.bc_pconst, nlevdecomp)
    set_bcpconst(fates.bc_pconst, nlevdecomp)

    # Surface-PFT span of the per-site pft_areafrac / hlm_sp_* arrays. The host path
    # carries HLM PFTs 1:surfpft_ub (bare ground, HLM PFT 0, goes to baregroundfrac);
    # the legacy FATES-PFT-space path keeps 1:numpft.
    surfpft_ub_bc = host_biogeog ? length(hlm_pft_areafrac[1]) - 1 : numpft_resolved
    host_biogeog && length(hlm_pft_areafrac) != nsites &&
        error("clm_fates_init!: hlm_pft_areafrac has $(length(hlm_pft_areafrac)) " *
              "site vectors; expected nsites=$nsites")
    # Host LUH2 data (clm_initializeMod: dynFatesLandUseInit before clm_fates%Init;
    # GetLandusePFTData inside Init when fixed biogeography + LUH). One gridcell per
    # FATES site's column.
    fcol_lu = (use_luh || hlm_use_lu_harvest[] == itrue) ?
              fates_site_columns(inst.column, inst.landunit, nsites) : Int[]
    site_g = [fcol_lu[s] == 0 ? 0 : inst.column.gridcell[fcol_lu[s]] for s in eachindex(fcol_lu)]
    if use_luh
        ng = length(inst.gridcell.area)
        fates.luh = dyn_fates_landuse_init(ng, fluh_timeseries; harvest_mode = harvest_mode,
                                           use_potentialveg = use_potentialveg,
                                           current_year = current_year)
        use_fixed_biogeog && get_landuse_pft_data!(fates.luh, flandusepftdat, ng)
    end
    for s in 1:nsites
        allocate_bcin(fates.bc_in[s], nlevsoil, nlevdecomp,
                      hlm_num_lu_harvest_cats[],   # num_lu_harvest_cats
                      hlm_num_luh2_states[],       # num_luh2_states
                      hlm_num_luh2_transitions[],  # num_luh2_transitions
                      1, surfpft_ub_bc)  # surfpft_lb, surfpft_ub
        allocate_bcout(fates.bc_out[s], nlevsoil, nlevdecomp)
        zero_bcs(fates, s)
        set_bcs(fates.bc_in[s])

        # Fill the minimal bc_in fields the cold-start chain reads.
        bc = fates.bc_in[s]
        bc.nlevsoil   = nlevsoil
        bc.nlevdecomp = nlevdecomp
        # Soil layering (init_site_vars! copies these into the site; PrepCH4/Nutrient
        # read dz_decomp / decomp_id during ed_update_site).
        bc.zi_sisl  .= collect(0.0:0.1:0.1 * nlevsoil)
        bc.dz_sisl  .= fill(0.1, nlevsoil)
        bc.z_sisl   .= collect(0.05:0.1:0.1 * nlevsoil)
        bc.dz_decomp_sisl .= fill(0.1, nlevdecomp)
        # non-vertsoilc => one decomp layer; every soil layer maps to decomp layer 1.
        fill!(bc.decomp_id, 1)
        bc.max_rooting_depth_index_col = min(3, nlevsoil)
        # (hlm_harvest_rates / catnames are sized by allocate_bcin: 0 unless an LUH2
        # harvest mode is on.)
        bc.hlm_harvest_units    = fates_unset_int
        bc.site_area            = area

        # Pass the grid-cell derived attributes to the site (clmfates_interfaceMod.F90
        # :880-884: `sites(s)%lat = grc%latdeg(g)`, `sites(s)%lon = grc%londeg(g)` with
        # g = col%gridcell(fcolumn(s))). This was never ported, so `site.lat` stayed NaN
        # and FATES phenology (`currentSite%lat > 0`, EDPhysiologyMod.F90:1033, :1087) took the
        # SOUTHERN-hemisphere branch everywhere (gddstart = 181, ncdstart = 120). The
        # other assignments in that Fortran block -- `h_gid` (set below) and the
        # fixed-biogeog `pft_areafrac` (set just below) -- are already wired. Site s
        # maps onto column s here (one site per column, `clm_initialize!`).
        let colg = inst.column.gridcell, grc = inst.gridcell
            if s <= length(colg)
                g = colg[s]
                latd = grc.latdeg[g]; lond = grc.londeg[g]
                # latdeg/londeg can be unfilled on a hand-built surfdata; derive them
                # from the radian lat/lon as the harnesses do.
                isfinite(latd) || (latd = rad2deg(grc.lat[g]))
                isfinite(lond) || (lond = rad2deg(grc.lon[g]))
                fates.sites[s].lat = latd
                fates.sites[s].lon = lond
            end
        end

        # Fixed-biogeog PFT presence (climate-appropriate seeding).  `pft_areafrac`
        # is allocated (length numpft) by `allocate_bcin` only when fixed-biogeog is
        # on; `baregroundfrac` must be finite (the no-LUH `_set_site_area_pft!` branch
        # has no NaN guard).  A fully-vegetated FATES column => bareground 0.
        if host_biogeog && use_luh
            # clmfates_interfaceMod init, LUH branch: the static landuse x pft map
            # (natural land-use types) and bare ground from flandusepftdat.
            g = site_g[s] == 0 ? 1 : site_g[s]
            for k in 1:num_landuse_pft_vars, ft in 1:size(bc.pft_areafrac_lu, 1)
                bc.pft_areafrac_lu[ft, k] = fates.luh.landuse_pft_map[g, ft, k]
            end
            bc.baregroundfrac = fates.luh.landuse_bareground[g]
        elseif host_biogeog
            # clmfates_interfaceMod init: pft_areafrac(ft) = wt_nat_patch(g, ft) over
            # surfpft_lb:surfpft_ub, endrun unless the areas (bare ground included) sum
            # to one. The ported _set_site_area_pft! reads bare ground from
            # baregroundfrac (Fortran: pft_areafrac(0)).
            wts = hlm_pft_areafrac[s]
            abs(sum(wts) - 1.0) > 1.0e-13 &&   # clm_varcon sum_to_1_tol
                error("clm_fates_init!: pft_area error in interfc, site $s: " *
                      "sum(wt_nat_patch) - 1 = $(sum(wts) - 1.0)")
            for k in 1:surfpft_ub_bc
                bc.pft_areafrac[k] = wts[k + 1]
            end
            bc.baregroundfrac = wts[1]
        elseif use_fixed_biogeog
            bc.pft_areafrac   .= fates_areafrac_vec
            bc.baregroundfrac  = 0.0
        end
    end

    # Fortran `init` follows the per-site bc allocation with init_soil_depths, BEFORE
    # init_coldstart, so the cold-start site copies the HLM's real soil layering. The
    # legacy one-decomp-layer path keeps its synthetic 0.1 m layering unchanged.
    use_vertsoilc && fates_init_soil_depths!(fates, inst)

    # clm_fates%init (clmfates_interfaceMod.F90:913-915): after init_soil_depths,
    # allocate each site's hydraulics object and its rhizosphere layering.
    # The rhizosphere layers are built from the host soil layering, so the
    # legacy synthetic 0.1 m layering (no init_soil_depths) is refused here.
    if use_planthydro
        use_vertsoilc || error("clm_fates_init!: use_planthydro needs the host soil " *
                               "layering (init_soil_depths), i.e. use_vertsoilc=true " *
                               "(the CTSM FATES soil-BGC configuration).")
        InitHydrSites(fates.sites, fates.bc_in)
    end

    # ---- W2: cold-start chain (init_coldstart) ----
    sites = fates.sites
    bc_in = fates.bc_in
    for s in 1:nsites
        init_site_vars!(sites[s], bc_in[s], fates.bc_out[s])
        zero_site!(sites[s])
    end
    set_site_properties!(nsites, sites, bc_in)
    # SP: hand FATES the HLM satellite-phenology LAI/SAI/height before init_patches
    # (init_coldstart; the host ran SatellitePhenology just before).
    if use_sp
        fcol_sp = fates_site_columns(inst.column, inst.landunit, nsites)
        for s in 1:nsites
            fcol_sp[s] == 0 && continue
            fates_pack_bcin_sp!(inst; s = s, c = fcol_sp[s])
        end
    end
    # init_coldstart (clmfates_interfaceMod.F90:2119-2157): the site rhizosphere
    # shells must exist before init_patches! creates the cohorts. Pack the host soil
    # (incl. the ice-reduced effective porosity) and cold-start the site hydraulics.
    if use_planthydro
        hcol = fates_site_columns(inst.column, inst.landunit, nsites)
        for s in 1:nsites
            hcol[s] == 0 && continue
            fates_pack_bcin_hydro_soil!(bc_in[s], inst; c=hcol[s], nlevsoil=nlevsoil,
                                        coldstart=true)
        end
        HydrSiteColdStart(sites, bc_in)
    end
    # init_coldstart: hand FATES the LUH2 states / transitions (and LUH harvest) of
    # the site's gridcell before init_patches.
    if use_luh
        # (Sites beyond the host's soil columns -- nsites is the column count -- never
        # run the daily step; they take gridcell 1's data so their cold start is valid.)
        for s in 1:nsites
            fates_pack_bcin_landuse!(bc_in[s], fates.luh, site_g[s] == 0 ? 1 : site_g[s])
        end
    end
    if use_inventory_init
        # Inventory cold start (EDInitMod init_patches -> initialize_sites_by_inventory).
        # The Fortran matches the inventory list against sites(s)%lat/lon, set from the
        # host gridcell (clmfates_interfaceMod.F90:884-885); the port's site lat/lon
        # are not set on the default path, so the host coordinates are passed in.
        grc = inst.gridcell
        site_coords = map(1:nsites) do s
            g = s <= length(inst.column.gridcell) ? inst.column.gridcell[s] : 1
            isempty(grc.lat) ? (NaN, NaN) : (rad2deg(grc.lat[g]), rad2deg(grc.lon[g]))
        end
        init_patches!(nsites, sites, bc_in;
                      sitelist_lines = read_lines_file(inventory_ctrl_file),
                      site_coords = site_coords)
    else
        init_patches!(nsites, sites, bc_in)
    end
    for s in 1:nsites
        # CLM has not computed last year's maximum rooting depth yet (it will at the
        # start of the timestep loop), so the litter-flux boundary starts trivial.
        bc_in[s].max_rooting_depth_index_col = bc_in[s].nlevsoil
        ed_update_site(sites[s], bc_in[s], fates.bc_out[s], false)
    end
    # Fortran `init_coldstart` follows the ed_update_site loop with
    #     call this%wrap_update_hlmfates_dyn(nc, bounds_clump, ..., .false.)
    # which runs canopy_summarization + update_hlm_dynamics. Without it the cold-start
    # cohorts have NO leaf-area profile: `nv` (number of leaf layers) stays 0 and the
    # patch `total_canopy_area` stays NaN, so the FATES canopy has zero leaf layers on
    # the very first timestep — photosynthesis and leaf dark respiration are identically
    # zero for that step, and the HLM canopy diagnostics are NaN. (Caught by the
    # time-stepped Fortran-FATES parity harness: Fortran has nv=2 / tcanarea=6260 at
    # nstep=0 where the port had nv=0 / NaN.) The HLM unpack is NOT done here: the
    # cold-start caller owns the canopystate/filter wiring.
    canopy_summarization!(nsites, sites, bc_in)

    # ---- W6: history-output interface (registry + buffers) ----
    # Each FATES site `s` writes into history-output slot `h_gid = s` (the Fortran
    # `f2hmap` site<->io-column map). For our single-clump / nsites layout the io
    # column space is just 1:nsites, so h_gid = s.
    for s in 1:nsites
        sites[s].h_gid = s
    end
    fates.hist = _build_fates_history_interface(nsites)

    # ---- Seed dispersal (clm_fates%init, clmfates_interfaceMod.F90:756-766) ----
    seeddisp_cadence != fates_dispersal_cadence_none &&
        fates_seed_dispersal_init!(fates, inst)

    inst.fates = fates

    # ---- Cold-start canopy transfer into the HLM (wrap_update_hlmfates_dyn) ----
    # Fortran's init_coldstart calls wrap_update_hlmfates_dyn (clmfates_interfaceMod.F90:2202)
    # so that the canopy state the HLM reads on the FIRST timestep — before any daily
    # demographic step has run — is FATES's, not the cold-start NaN. Without this the
    # first albedo/radiation pass on a FATES column classifies its patches off NaN
    # elai/esai. Flush the whole patch range (bare ground -> exactly zero), then fill
    # the vegetated patches from bc_out, then rebuild the HLM patch weights.
    col = inst.column
    if !isempty(col.is_fates)
        site_col0 = zeros(Int, nsites)
        s0 = 0
        for c in 1:length(col.is_fates)
            col.is_fates[c] || continue
            s0 += 1
            s0 <= nsites || break
            site_col0[s0] = c
        end
        canopy_summarization!(nsites, sites, bc_in)
        update_hlm_dynamics!(nsites, sites, site_col0, fates.bc_out)
        for s in 1:nsites
            c = site_col0[s]
            c == 0 && continue
            fates_flush_hlm_canopy!(inst; c=c)
            for (ifp, p) in fates_veg_patches(sites[s], c, col)
                fates_unpack_bcout_canopy_structure!(inst; s=s, c=c, p=p, ifp=ifp)
            end
            fates_set_filters!(inst; c=c, s=s)
        end
    end

    # wrap_update_hlmfates_dyn (clmfates_interfaceMod.F90:1552-1558): the cold-start
    # plant water enters the host column water from the first step's begwb on. Mapped
    # through the natural-veg columns (f2hmap%fcolumn), which does not depend on the
    # caller having flagged col.is_fates yet.
    if use_planthydro
        hcol = fates_site_columns(inst.column, inst.landunit, nsites)
        for s in 1:nsites
            hcol[s] == 0 && continue
            inst.water.waterdiagnosticbulk_inst.total_plant_stored_h2o_col[hcol[s]] =
                fates.bc_out[s].plant_stored_h2o_si
        end
    end

    return fates
end

"""
    fates_site_columns(col, lun, nsites) -> Vector{Int}

HLM column hosting each FATES site: the `s`-th natural-vegetation (`ISTSOIL`) column,
the Fortran `f2hmap%fcolumn` for a single clump (one FATES site per natural-veg column).
Sites without a host column map to 0.
"""
function fates_site_columns(col, lun, nsites::Integer)
    fc = zeros(Int, nsites)
    s = 0
    for c in 1:length(col.landunit)
        lun.itype[col.landunit[c]] == ISTSOIL || continue
        s += 1
        s <= nsites || break
        fc[s] = c
    end
    return fc
end

"""
    fates_init_soil_depths!(fates, inst)

Port of `init_soil_depths` (clmfates_interfaceMod.F90:3397-3442): copy the host
column's soil layering into each FATES site's `bc_in` — interfaces `zi_sisl(0:nlevsoil)`,
thicknesses `dz_sisl`, node depths `z_sisl`, the decomposition-layer thicknesses
`dz_decomp_sisl = dzsoi_decomp(1:nlevdecomp)` and the identity soil→decomp map
`decomp_id(j) = j`. Fortran aborts if a decomposition layer does not coincide with its
soil layer (FATES assumes they are equal when vertical decomposition is on).
"""
function fates_init_soil_depths!(fates, inst::CLMInstances)
    col = inst.column
    fcol = fates_site_columns(col, inst.landunit, fates.nsites)
    joff = varpar.nlevsno
    dzd = dzsoi_decomp[]
    for s in 1:fates.nsites
        c = fcol[s]
        c == 0 && continue
        bc = fates.bc_in[s]
        nlevsoil = bc.nlevsoil
        nlevdecomp = bc.nlevdecomp
        for j in 0:nlevsoil
            bc.zi_sisl[j + 1] = col.zi[c, joff + 1 + j]
        end
        for j in 1:nlevsoil
            bc.dz_sisl[j] = col.dz[c, joff + j]
            bc.z_sisl[j]  = col.z[c, joff + j]
        end
        for j in 1:nlevdecomp
            bc.dz_decomp_sisl[j] = dzd[j]
        end
        for j in 1:nlevsoil
            bc.decomp_id[j] = j
            if abs(bc.dz_decomp_sisl[j] - bc.dz_sisl[j]) > 1.0e-10
                error("fates_init_soil_depths!: with vertical soil decomposition FATES " *
                      "assumes decomposition depths equal soil depths (layer $j: " *
                      "dz_decomp=$(bc.dz_decomp_sisl[j]), dz=$(bc.dz_sisl[j]))")
            end
        end
    end
    return fates
end

"""
    _build_fates_history_interface(nsites) -> fates_history_interface_type

Instantiate and fully build the FATES history-output interface for `nsites`
sites: dimension bookkeeping (`Init!`), per-thread bounds, the dim-kind maps
(`assemble_history_output_types!`), and the two-pass variable registry +
buffer allocation (`initialize_history_vars!`). Mirrors the Fortran host
sequence (`fates_hist%Init` -> `SetThreadBoundsEach` -> `assemble_*` ->
`initialize_history_vars`).

The non-column dimension sizes are taken from the FATES class-count Refs
(`nlevsclass[]`, `nlevage[]`, ...) populated by the param read / ctrlparms set;
only the `column` (site) dimension drives the site-scalar buffers that the
W1 fills (dyn1 core + hifrq1) actually write. The level-2 (scpf/scag/...) dims
size the not-yet-filled dyn2/hifrq2 buffers.
"""
function _build_fates_history_interface(nsites::Integer)
    # Resolve class-count dims with safe fallbacks (these Refs are set by the
    # param read; guard against an unset Ref so the registry still allocates).
    _dim(r::Ref{Int}, fallback::Int) = (r[] > 0 ? r[] : fallback)
    npft   = _dim(numpft, 1)
    nscls  = _dim(nlevsclass, 1)
    nage   = _dim(nlevage, 1)
    nheight= _dim(nlevheight, 1)
    ncoage = _dim(nlevcoage, 1)
    ndamage= _dim(nlevdamage, 1)
    nelem  = max(num_elements[], 1)
    ncwd_n = ncwd                  # coarse woody debris size classes (=4)
    nfuel  = num_fuel_classes      # fuel size classes (=6)
    ncan   = nclmax                # number of canopy layers
    nlu    = n_landuse_cats

    fb = fates_bounds_type(
        column_begin = 1,                     column_end = nsites,
        soil_begin = 1,                       soil_end = max(nlevsoil_hist(), 1),
        sizepft_class_begin = 1,              sizepft_class_end = nscls * npft,
        size_class_begin = 1,                 size_class_end = nscls,
        coage_class_begin = 1,                coage_class_end = ncoage,
        coagepf_class_begin = 1,              coagepf_class_end = ncoage * npft,
        pft_class_begin = 1,                  pft_class_end = npft,
        age_class_begin = 1,                  age_class_end = nage,
        fuel_begin = 1,                       fuel_end = nfuel,
        cwdsc_begin = 1,                      cwdsc_end = ncwd_n,
        can_begin = 1,                        can_end = ncan,
        cnlf_begin = 1,                       cnlf_end = ncan * nlevleaf,
        cnlfpft_begin = 1,                    cnlfpft_end = ncan * nlevleaf * npft,
        cdpf_begin = 1,                       cdpf_end = ndamage * nscls * npft,
        cdsc_begin = 1,                       cdsc_end = ndamage * nscls,
        cdam_begin = 1,                       cdam_end = ndamage,
        sizeage_class_begin = 1,              sizeage_class_end = nscls * nage,
        sizeagepft_class_begin = 1,           sizeagepft_class_end = nscls * nage * npft,
        agepft_class_begin = 1,               agepft_class_end = nage * npft,
        height_begin = 1,                     height_end = nheight,
        elem_begin = 1,                       elem_end = nelem,
        elpft_begin = 1,                      elpft_end = nelem * npft,
        elcwd_begin = 1,                      elcwd_end = nelem * ncwd_n,
        elage_begin = 1,                      elage_end = nelem * nage,
        agefuel_begin = 1,                    agefuel_end = nage * nfuel,
        clscpf_begin = 1,                     clscpf_end = ncan * nscls * npft,
        landuse_begin = 1,                    landuse_end = nlu,
        lulu_begin = 1,                       lulu_end = nlu * nlu,
        lupft_begin = 1,                      lupft_end = nlu * npft,
    )

    hist = fates_history_interface_type()
    Init!(hist, 1, fb)                       # dim bookkeeping (1 thread)
    SetThreadBoundsEach!(hist, 1, fb)        # per-thread (clump) bounds
    assemble_history_output_types!(hist)     # dim-kind maps + dim indices
    initialize_history_vars!(hist)           # count + allocate + Init buffers
    return hist
end

# Soil-layer dim for the history buffers. The site-scalar level-1 fills
# (dyn1 core + hifrq1) don't index soil layers, so a size of 1 is sufficient
# for them; the soil-dimensioned dyn2 vars are not yet filled.
nlevsoil_hist() = 1
