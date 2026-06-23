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
#   * planthydro OFF, SPITFIRE/fire OFF (no_fire), inventory-init OFF
#   * SP OFF, nocomp OFF, fixed-biogeog OFF, LUH OFF, harvest/logging OFF
#   * tree-damage OFF, cohort-age-tracking OFF, ch4 OFF, vertsoilc OFF (=> nlevdecomp 1)
# This is the default NBG (natural-bare-ground) cold-start: init_patches! builds a
# single primaryland patch holding the whole site area, init_cohorts! seeds one
# small cohort per PFT.  It needs none of the SP-LAI / inventory / biogeog inputs.
#
# Since there is no FATES NetCDF parameter file wired into this port, the PFT /
# allometry / ED parameter tables are populated synthetically in-memory (a 2-PFT
# woody evergreen carbon table — the same table the EDInitMod unit test uses).

"""
    _fates_spike_setup_pft!(npft) -> npft

Populate the in-memory FATES parameter tables (prt_params allometry, param_derived,
EDPftvarcon, ed_params bins/tols, sf_params CWD frac) for a minimal `npft`-PFT
woody-evergreen carbon-only cold start, and set the size/age/coage/damage level
counts.  Mirrors the synthetic table the EDInitMod test builds — used in lieu of a
FATES parameter NetCDF file (not wired into this port).
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
    EDPftvarcon_inst[] = ev

    edp = ed_params_type()
    edp.regeneration_model = default_regeneration
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
    _fates_spike_set_ctrlparms!(; numpft_in, nlevsoil)

Set the FATES `hlm_*` control Refs for the minimal carbon-only cold start, mirroring
the `set_fates_ctrlparms` tag sequence in `CLMFatesGlobals1`/`CLMFatesGlobals2`
(see clmfates_interfaceMod.F90).  Brackets with `flush_to_unset` + `check_allset`.
"""
function _fates_spike_set_ctrlparms!(; numpft_in::Int, nlevsoil::Int)
    set_fates_ctrlparms("flush_to_unset")

    # Biogeography / competition modes — all OFF (default NBG cold start).
    set_fates_ctrlparms("use_fixed_biogeog"; ival = ifalse)
    set_fates_ctrlparms("use_nocomp";        ival = ifalse)
    set_fates_ctrlparms("use_sp";            ival = ifalse)
    set_fates_ctrlparms("masterproc";        ival = itrue)

    # Radiation bands / soil dims / HLM identity.
    set_fates_ctrlparms("num_sw_bbands"; ival = num_swb)
    set_fates_ctrlparms("vis_sw_index";  ival = ivis)
    set_fates_ctrlparms("nir_sw_index";  ival = inir)
    set_fates_ctrlparms("num_lev_soil";  ival = nlevsoil)
    set_fates_ctrlparms("hlm_name";      cval = "CLM")
    set_fates_ctrlparms("hio_ignore_val"; rval = -9.9e30)
    set_fates_ctrlparms("soilwater_ipedof"; ival = 0)

    # PARTEH carbon-only + seed-dispersal cadence + history levels.
    set_fates_ctrlparms("parteh_mode";       ival = prt_carbon_allom_hyp)
    set_fates_ctrlparms("seeddisp_cadence";  ival = 0)
    set_fates_ctrlparms("hist_hifrq_dimlevel"; ival = 1)
    set_fates_ctrlparms("hist_dynam_dimlevel"; ival = 1)

    # Nutrient competition / decomposition.  Carbon-only => N/P spec off, CTC decomp.
    set_fates_ctrlparms("nu_com";        cval = "RD")
    set_fates_ctrlparms("decomp_method"; cval = "CTC")
    set_fates_ctrlparms("use_tree_damage"; ival = ifalse)
    set_fates_ctrlparms("nitrogen_spec";   ival = 0)
    set_fates_ctrlparms("phosphorus_spec"; ival = 0)

    # SPITFIRE / fire — OFF (no_fire mode).  These four definition Refs select the
    # mode enumerations; 0 selects the "off" branch of each.
    set_fates_ctrlparms("spitfire_mode"; ival = 0)
    set_fates_ctrlparms("sf_nofire_def"; ival = 0)
    set_fates_ctrlparms("sf_scalar_lightning_def";    ival = 1)
    set_fates_ctrlparms("sf_successful_ignitions_def"; ival = 2)
    set_fates_ctrlparms("sf_anthro_ignitions_def";    ival = 3)

    set_fates_ctrlparms("is_restart"; ival = ifalse)
    set_fates_ctrlparms("use_ch4";        ival = ifalse)
    set_fates_ctrlparms("use_vertsoilc";  ival = ifalse)
    set_fates_ctrlparms("use_ed_st3";     ival = ifalse)
    set_fates_ctrlparms("use_ed_prescribed_phys"; ival = ifalse)
    set_fates_ctrlparms("use_planthydro"; ival = ifalse)
    set_fates_ctrlparms("use_cohort_age_tracking"; ival = ifalse)

    # Land-use harvest / logging / LUH — OFF.
    set_fates_ctrlparms("use_lu_harvest";     ival = ifalse)
    set_fates_ctrlparms("num_lu_harvest_cats"; ival = 0)
    set_fates_ctrlparms("use_logging";        ival = ifalse)
    set_fates_ctrlparms("use_luh2";           ival = ifalse)
    set_fates_ctrlparms("num_luh2_states";      ival = 0)
    set_fates_ctrlparms("num_luh2_transitions"; ival = 0)
    set_fates_ctrlparms("use_fates_potentialveg"; ival = ifalse)

    # Inventory-init OFF (so cold-start NBG path runs).  The control file must be
    # "set" to a non-"unset" string to pass check_allset even though it's unused.
    set_fates_ctrlparms("use_inventory_init"; ival = ifalse)
    set_fates_ctrlparms("inventory_ctrl_file"; cval = "none")

    set_fates_ctrlparms("check_allset")
    return nothing
end

"""
    clm_fates_init!(inst::CLMInstances; nsites=1, numpft_in, nlevsoil, nlevdecomp,
                    current_year=1, current_month=1, current_day=1)
        -> fates_interface_type

Bootstrap a finite, cold-started, carbon-only FATES interface with `nsites` sites
and attach it to `inst.fates`.  Mirrors the Fortran host sequence
`CLMFatesGlobals1 -> CLMFatesGlobals2 -> clm_fates%Init -> init_coldstart`.

Carbon-only minimal config (planthydro / SPITFIRE / inventory all OFF; SP / nocomp /
fixed-biogeog / LUH all OFF) — the default NBG cold start that needs no external
SP-LAI / inventory / biogeography data.

This sets the FATES module-global `hlm_*` control Refs and the in-memory parameter
tables as a side effect (there is no FATES parameter NetCDF wired into this port).
The default test suite never calls this, so the globals are left in the minimal-
cold-start state on return; callers that mix FATES with other FATES tests should
save/restore the globals themselves.
"""
function clm_fates_init!(inst::CLMInstances; nsites::Int = 1, numpft_in::Int,
                         nlevsoil::Int, nlevdecomp::Int,
                         current_year::Int = 1, current_month::Int = 1,
                         current_day::Int = 1)

    # ---- W1.1: FATES globals + parameters (CLMFatesGlobals1) ----
    FatesInterfaceInit(6, false)

    # Populate the in-memory parameter tables (stand-in for the FATES param file).
    _fates_spike_setup_pft!(numpft_in)

    # Carbon-only element registry.
    num_elements[] = 1
    empty!(element_list)
    push!(element_list, carbon12_element)

    # Control parameters (set_fates_ctrlparms tag sequence).
    _fates_spike_set_ctrlparms!(; numpft_in = numpft_in, nlevsoil = nlevsoil)

    # PARTEH carbon hypothesis global state + numpft / nleafage / max-patch dims.
    InitPRTGlobalAllometricCarbon!()
    numpft[]   = numpft_in
    nleafage[] = size(prt_params.leaf_long, 2)

    # SetFatesGlobalElements1 sizes maxpatch_total / fates_maxPatchesPerSite from the
    # ed_params.maxpatches_by_landuse table we set above (non-SP / non-nocomp branch:
    # fates_maxPatchesPerSite = maxpatch_total + 1).
    p = ed_params()
    p.maxpatch_total = sum(p.maxpatches_by_landuse)
    fates_maxPatchesPerSite[] = p.maxpatch_total + 1

    # SetFatesGlobalElements2 derived dims the cold-start chain + bc allocators read.
    # Carbon-only => one nutrient competitor / trivial competition scaling; the full
    # SetFatesGlobalElements2 (InitHydroGlobals / fates_history_maps / VAI bins) is
    # not needed for a finite carbon-only cold start, so set just these directly.
    max_comp_per_site[]     = 1
    fates_np_comp_scaling[] = trivial_np_comp_scaling
    fates_maxElementsPerSite[] = max(fates_maxPatchesPerSite[] * num_swb,
                                     nlevsclass[] * numpft[] * n_term_mort_types)

    # ---- W1.2: FATES clock (SetFatesTime) ----
    SetFatesTime(current_year, current_month, current_day,
                 0,                       # current_tod (seconds into day)
                 current_year * 10000 + current_month * 100 + current_day,  # current_date
                 10101,                   # reference_date
                 0.0,                     # model_day
                 1,                       # day_of_year
                 365,                     # days_per_year
                 1.0 / 365.0)             # freq_day

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

    for s in 1:nsites
        allocate_bcin(fates.bc_in[s], nlevsoil, nlevdecomp,
                      0,   # num_lu_harvest_cats
                      0,   # num_luh2_states
                      0,   # num_luh2_transitions
                      1, numpft_in)  # surfpft_lb, surfpft_ub
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
        bc.hlm_harvest_rates    = Float64[]
        bc.hlm_harvest_catnames = String[]
        bc.hlm_harvest_units    = fates_unset_int
        bc.site_area            = area
    end

    # ---- W2: cold-start chain (init_coldstart) ----
    sites = fates.sites
    bc_in = fates.bc_in
    for s in 1:nsites
        init_site_vars!(sites[s], bc_in[s], fates.bc_out[s])
        zero_site!(sites[s])
    end
    set_site_properties!(nsites, sites, bc_in)
    init_patches!(nsites, sites, bc_in)
    for s in 1:nsites
        ed_update_site(sites[s], bc_in[s], fates.bc_out[s], false)
    end

    inst.fates = fates
    return fates
end
