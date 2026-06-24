# FatesHistoryInterfaceMod.jl
# Julia port of FATES src/fates/main/FatesHistoryInterfaceMod.F90 (9332 lines).
#
# The `fates_history_interface_type` owns:
#   * `hvars`     — the registry Vector{fates_history_variable_type} (allocated to
#                   exactly the number of variables active for this HLM build),
#   * `dim_kinds` — fates_history_num_dim_kinds variable-kind descriptors,
#   * `dim_bounds`— fates_history_num_dimensions dimension descriptors,
#   * the per-dimension integer index handles (column_index_, levsoil_index_, ...),
#   * the per-variable `ih_*` index handles (stored in `ih`, a Dict mapping the
#     handle name -> hvars slot, so update routines can address a var by name).
#
# What is ported in FULL and faithfully:
#   * The type + all 30 dimension index handles + all 475 ih_* var handles.
#   * Init (dim_bounds setup), SetThreadBoundsEach, assemble_history_output_types,
#     set_dim_indices, init_dim_kinds_maps, the index get/set accessors.
#   * set_history_var, flush_hvars/flush_all_hvars/zero_site_hvars.
#   * define_history_vars: the COMPLETE registry of all 472 active history
#     variables (the 1 commented-out FATES_L2FR_CLSZPF is intentionally omitted,
#     matching Fortran) is reproduced 1:1 via the data-driven `_HISTORY_VAR_REGISTRY`
#     table — names/units/long/use_default/avgflag/vtype/hlms/upfreq/handle/
#     flush_to_zero and the exact HLM gating (hist-level / planthydro / tree-damage /
#     nocomp / element presence / not-st3-not-sp) are preserved, in source order.
#
# What is ported as a representative update path:
#   * update_history_dyn! orchestrator + update_history_dyn1! covering the
#     patch/cohort COUNTING and the simple site-level state copies that map onto
#     already-ported site/patch/cohort fields. The deep per-cohort PRT / element /
#     hydraulics fills (which require FATES state not yet wired here) are marked
#     `# TODO Batch 18-followup:` in update_history_dyn2!/hifrq/hydraulics, which
#     are present as orchestrator stubs.
#
# Deps: FatesHistoryVariableType, FatesIODimensionsMod, FatesIOVariableKindMod,
#       FatesInterfaceTypesMod (hlm_* flags, numpft, ...), EDTypesMod (ed_site_type,
#       AREA, AREA_INV), FatesPatchMod, FatesCohortMod, FatesConstantsMod,
#       PRTGenericMod (element_list, *_element), FatesGlobals.
#
# Upstream-Fortran quirks preserved:
#   * fates_history_num_dimensions == fates_history_num_dim_kinds == 50 (static),
#     but only ~30 of each slot is used (Init/init_dim_kinds_maps fill the first
#     ~30; the rest stay default). We allocate exactly 50 of each to match.
#   * FATES_RECRUITMENT_CFLUX_PF is registered with a literal `upfreq=1`
#     (== group_dyna_simple); we keep upfreq=group_dyna_simple. (banked)
#   * 3 ih_* handles are declared but never registered (ih_l2fr_clscpf,
#     ih_supsub_scpf, ih_crownarea_clll) — they stay 0 in `ih`, matching Fortran.

const fates_history_num_dimensions = 50
const fates_history_num_dim_kinds = 50

# ---------------------------------------------------------------------------
# The 30 dimension-name constants, in the exact order the Fortran Init pushes
# them onto dim_count. The Julia IODimensions port names them `dimname_*`.
# (dim_index, dimname, fates_bounds field-pair accessor) tuples.
# ---------------------------------------------------------------------------
const _DIM_SETUP = Tuple{Symbol,String,Symbol,Symbol}[
    (:column,     dimname_column,     :column_begin,          :column_end),
    (:levsoil,    dimname_levsoil,    :soil_begin,            :soil_end),
    (:levscpf,    dimname_levscpf,    :sizepft_class_begin,   :sizepft_class_end),
    (:levscls,    dimname_levscls,    :size_class_begin,      :size_class_end),
    (:levcacls,   dimname_levcacls,   :coage_class_begin,     :coage_class_end),
    (:levcapf,    dimname_levcapf,    :coagepf_class_begin,   :coagepf_class_end),
    (:levpft,     dimname_levpft,     :pft_class_begin,       :pft_class_end),
    (:levage,     dimname_levage,     :age_class_begin,       :age_class_end),
    (:levfuel,    dimname_levfuel,    :fuel_begin,            :fuel_end),
    (:levcwdsc,   dimname_levcwdsc,   :cwdsc_begin,           :cwdsc_end),
    (:levcan,     dimname_levcan,     :can_begin,             :can_end),
    (:levcnlf,    dimname_levcnlf,    :cnlf_begin,            :cnlf_end),
    (:levcnlfpft, dimname_levcnlfpft, :cnlfpft_begin,         :cnlfpft_end),
    (:levcdpf,    dimname_levcdpf,    :cdpf_begin,            :cdpf_end),
    (:levcdsc,    dimname_levcdsc,    :cdsc_begin,            :cdsc_end),
    (:levcdam,    dimname_levcdam,    :cdam_begin,            :cdam_end),
    (:levscag,    dimname_levscag,    :sizeage_class_begin,   :sizeage_class_end),
    (:levscagpft, dimname_levscagpft, :sizeagepft_class_begin,:sizeagepft_class_end),
    (:levagepft,  dimname_levagepft,  :agepft_class_begin,    :agepft_class_end),
    (:levheight,  dimname_levheight,  :height_begin,          :height_end),
    (:levelem,    dimname_levelem,    :elem_begin,            :elem_end),
    (:levelpft,   dimname_levelpft,   :elpft_begin,           :elpft_end),
    (:levelcwd,   dimname_levelcwd,   :elcwd_begin,           :elcwd_end),
    (:levelage,   dimname_levelage,   :elage_begin,           :elage_end),
    (:levagefuel, dimname_levagefuel, :agefuel_begin,         :agefuel_end),
    (:levclscpf,  dimname_levclscpf,  :clscpf_begin,          :clscpf_end),
    (:levlanduse, dimname_levlanduse, :landuse_begin,         :landuse_end),
    (:levlulu,    dimname_levlulu,    :lulu_begin,            :lulu_end),
    (:levlupft,   dimname_levlupft,   :lupft_begin,           :lupft_end),
]

# ---------------------------------------------------------------------------
# init_dim_kinds_maps registration order: (slot, variable-kind name, ndims).
# Matches the Fortran init_dim_kinds_maps exactly (29 kinds).
# ---------------------------------------------------------------------------
const _DIM_KIND_SETUP = Tuple{String,Int}[
    (site_r8, 1), (site_soil_r8, 2), (site_size_pft_r8, 2), (site_size_r8, 2),
    (site_coage_pft_r8, 2), (site_coage_r8, 2), (site_pft_r8, 2), (site_age_r8, 2),
    (site_fuel_r8, 2), (site_cwdsc_r8, 2), (site_can_r8, 2), (site_cnlf_r8, 2),
    (site_cnlfpft_r8, 2), (site_cdpf_r8, 2), (site_cdsc_r8, 2), (site_cdam_r8, 2),
    (site_scag_r8, 2), (site_scagpft_r8, 2), (site_agepft_r8, 2), (site_height_r8, 2),
    (site_elem_r8, 2), (site_elpft_r8, 2), (site_elcwd_r8, 2), (site_elage_r8, 2),
    (site_agefuel_r8, 2), (site_clscpf_r8, 2), (site_landuse_r8, 2), (site_lulu_r8, 2),
    (site_lupft_r8, 2),
]

# ---------------------------------------------------------------------------
# assemble_history_output_types dim-index map: (kind, idim, dim-handle-symbol).
# ---------------------------------------------------------------------------
const _ASSEMBLE_MAP = Tuple{String,Int,Symbol}[
    (site_r8, 1, :column),
    (site_soil_r8, 1, :column), (site_soil_r8, 2, :levsoil),
    (site_size_pft_r8, 1, :column), (site_size_pft_r8, 2, :levscpf),
    (site_size_r8, 1, :column), (site_size_r8, 2, :levscls),
    (site_coage_r8, 1, :column), (site_coage_r8, 2, :levcacls),
    (site_coage_pft_r8, 1, :column), (site_coage_pft_r8, 2, :levcapf),
    (site_pft_r8, 1, :column), (site_pft_r8, 2, :levpft),
    (site_age_r8, 1, :column), (site_age_r8, 2, :levage),
    (site_fuel_r8, 1, :column), (site_fuel_r8, 2, :levfuel),
    (site_cwdsc_r8, 1, :column), (site_cwdsc_r8, 2, :levcwdsc),
    (site_can_r8, 1, :column), (site_can_r8, 2, :levcan),
    (site_cnlf_r8, 1, :column), (site_cnlf_r8, 2, :levcnlf),
    (site_cnlfpft_r8, 1, :column), (site_cnlfpft_r8, 2, :levcnlfpft),
    (site_cdpf_r8, 1, :column), (site_cdpf_r8, 2, :levcdpf),
    (site_cdsc_r8, 1, :column), (site_cdsc_r8, 2, :levcdsc),
    (site_cdam_r8, 1, :column), (site_cdam_r8, 2, :levcdam),
    (site_scag_r8, 1, :column), (site_scag_r8, 2, :levscag),
    (site_scagpft_r8, 1, :column), (site_scagpft_r8, 2, :levscagpft),
    (site_agepft_r8, 1, :column), (site_agepft_r8, 2, :levagepft),
    (site_height_r8, 1, :column), (site_height_r8, 2, :levheight),
    (site_elem_r8, 1, :column), (site_elem_r8, 2, :levelem),
    (site_elpft_r8, 1, :column), (site_elpft_r8, 2, :levelpft),
    (site_elcwd_r8, 1, :column), (site_elcwd_r8, 2, :levelcwd),
    (site_elage_r8, 1, :column), (site_elage_r8, 2, :levelage),
    (site_agefuel_r8, 1, :column), (site_agefuel_r8, 2, :levagefuel),
    (site_clscpf_r8, 1, :column), (site_clscpf_r8, 2, :levclscpf),
    (site_landuse_r8, 1, :column), (site_landuse_r8, 2, :levlanduse),
    (site_lulu_r8, 1, :column), (site_lulu_r8, 2, :levlulu),
    (site_lupft_r8, 1, :column), (site_lupft_r8, 2, :levlupft),
]


# All 475 ih_* index-handle fields, declared in module scope in the Fortran.
# Stored on the interface type so each registered var records its hvars slot.
const _IH_HANDLE_NAMES = Symbol[
    :ih_storec_si, :ih_storectfrac_si, :ih_storectfrac_canopy_scpf, :ih_storectfrac_ustory_scpf, :ih_leafc_si, :ih_sapwc_si,
    :ih_fnrtc_si, :ih_fnrtc_sl, :ih_reproc_si, :ih_totvegc_si, :ih_storen_si, :ih_leafn_si,
    :ih_sapwn_si, :ih_fnrtn_si, :ih_repron_si, :ih_totvegn_si, :ih_storentfrac_si, :ih_totvegn_scpf,
    :ih_leafn_scpf, :ih_fnrtn_scpf, :ih_storen_scpf, :ih_sapwn_scpf, :ih_repron_scpf, :ih_storentfrac_canopy_scpf,
    :ih_storentfrac_understory_scpf, :ih_storep_si, :ih_leafp_si, :ih_sapwp_si, :ih_fnrtp_si, :ih_reprop_si,
    :ih_totvegp_si, :ih_storeptfrac_si, :ih_totvegp_scpf, :ih_leafp_scpf, :ih_fnrtp_scpf, :ih_reprop_scpf,
    :ih_storep_scpf, :ih_sapwp_scpf, :ih_storeptfrac_canopy_scpf, :ih_storeptfrac_understory_scpf, :ih_l2fr_si, :ih_l2fr_clscpf,
    :ih_recl2fr_canopy_pf, :ih_recl2fr_ustory_pf, :ih_nh4uptake_scpf, :ih_no3uptake_scpf, :ih_puptake_scpf, :ih_nh4uptake_si,
    :ih_no3uptake_si, :ih_puptake_si, :ih_nefflux_si, :ih_pefflux_si, :ih_nefflux_scpf, :ih_pefflux_scpf,
    :ih_nfix_si, :ih_nfix_scpf, :ih_ndemand_si, :ih_ndemand_scpf, :ih_pdemand_si, :ih_pdemand_scpf,
    :ih_trimming_si, :ih_area_plant_si, :ih_area_trees_si, :ih_litter_in_elem, :ih_litter_out_elem, :ih_seed_bank_elem,
    :ih_fates_fraction_si, :ih_litter_in_si, :ih_litter_out_si, :ih_seed_bank_si, :ih_seeds_in_si, :ih_seeds_in_local_si,
    :ih_ungerm_seed_bank_si, :ih_seedling_pool_si, :ih_ba_weighted_height_si, :ih_ca_weighted_height_si, :ih_seeds_in_local_elem, :ih_seeds_in_extern_elem,
    :ih_seed_decay_elem, :ih_seed_germ_elem, :ih_fines_ag_elem, :ih_fines_bg_elem, :ih_cwd_ag_elem, :ih_cwd_bg_elem,
    :ih_cwd_elcwd, :ih_burn_flux_elem, :ih_bstor_canopy_si_scpf, :ih_bstor_understory_si_scpf, :ih_bleaf_canopy_si_scpf, :ih_bleaf_understory_si_scpf,
    :ih_lai_canopy_si_scpf, :ih_lai_understory_si_scpf, :ih_crownarea_canopy_si_scpf, :ih_crownarea_understory_si_scpf, :ih_totvegc_scpf, :ih_leafc_scpf,
    :ih_fnrtc_scpf, :ih_storec_scpf, :ih_sapwc_scpf, :ih_reproc_scpf, :ih_bdead_si, :ih_balive_si,
    :ih_agb_si, :ih_npp_si, :ih_gpp_si, :ih_aresp_si, :ih_maint_resp_si, :ih_growth_resp_si,
    :ih_excess_resp_si, :ih_ar_canopy_si, :ih_gpp_canopy_si, :ih_ar_understory_si, :ih_gpp_understory_si, :ih_canopy_biomass_si,
    :ih_understory_biomass_si, :ih_maint_resp_unreduced_si, :ih_npp_secondary_si, :ih_gpp_secondary_si, :ih_aresp_secondary_si, :ih_maint_resp_secondary_si,
    :ih_growth_resp_secondary_si, :ih_primaryland_fusion_error_si, :ih_area_si_landuse, :ih_disturbance_rate_si_lulu, :ih_transition_matrix_si_lulu, :ih_fire_disturbance_rate_si,
    :ih_logging_disturbance_rate_si, :ih_fall_disturbance_rate_si, :ih_harvest_debt_si, :ih_harvest_debt_sec_si, :ih_harvest_woodprod_carbonflux_si, :ih_luchange_woodprod_carbonflux_si,
    :ih_nplant_si_scag, :ih_nplant_canopy_si_scag, :ih_nplant_understory_si_scag, :ih_ddbh_canopy_si_scag, :ih_ddbh_understory_si_scag, :ih_mortality_canopy_si_scag,
    :ih_mortality_understory_si_scag, :ih_nplant_si_scagpft, :ih_biomass_si_agepft, :ih_npp_si_agepft, :ih_scorch_height_si_agepft, :ih_tveg24_si,
    :ih_tlongterm_si, :ih_tgrowth_si, :ih_tveg_si, :ih_nep_si, :ih_hr_si, :ih_c_stomata_si,
    :ih_c_lblayer_si, :ih_vis_rad_err_si, :ih_nir_rad_err_si, :ih_fire_c_to_atm_si, :ih_interr_liveveg_elem, :ih_interr_litter_elem,
    :ih_cbal_err_fates_si, :ih_err_fates_elem, :ih_npatches_si, :ih_npatches_sec_si, :ih_ncohorts_si, :ih_ncohorts_sec_si,
    :ih_demotion_carbonflux_si, :ih_promotion_carbonflux_si, :ih_canopy_mortality_carbonflux_si, :ih_understory_mortality_carbonflux_si, :ih_canopy_mortality_crownarea_si, :ih_understory_mortality_crownarea_si,
    :ih_canopy_spread_si, :ih_npp_leaf_si, :ih_npp_seed_si, :ih_npp_stem_si, :ih_npp_froot_si, :ih_npp_croot_si,
    :ih_npp_stor_si, :ih_leaf_mr_si, :ih_froot_mr_si, :ih_livestem_mr_si, :ih_livecroot_mr_si, :ih_fraction_secondary_forest_si,
    :ih_biomass_secondary_forest_si, :ih_woodproduct_si, :ih_h2oveg_si, :ih_h2oveg_dead_si, :ih_h2oveg_recruit_si, :ih_h2oveg_growturn_err_si,
    :ih_h2oveg_hydro_err_si, :ih_lai_si, :ih_elai_si, :ih_site_cstatus_si, :ih_gdd_si, :ih_site_nchilldays_si,
    :ih_site_ncolddays_si, :ih_cleafoff_si, :ih_cleafon_si, :ih_nesterov_fire_danger_si, :ih_fire_nignitions_si, :ih_fire_fdi_si,
    :ih_fire_intensity_area_product_si, :ih_spitfire_ros_si, :ih_effect_wspeed_si, :ih_tfc_ros_si, :ih_fire_intensity_si, :ih_fire_area_si,
    :ih_fire_fuel_bulkd_si, :ih_fire_fuel_eff_moist_si, :ih_fire_fuel_sav_si, :ih_fire_fuel_mef_si, :ih_sum_fuel_si, :ih_fragmentation_scaler_sl,
    :ih_nplant_si_scpf, :ih_gpp_si_scpf, :ih_npp_totl_si_scpf, :ih_npp_leaf_si_scpf, :ih_npp_seed_si_scpf, :ih_npp_fnrt_si_scpf,
    :ih_npp_bgsw_si_scpf, :ih_npp_bgdw_si_scpf, :ih_npp_agsw_si_scpf, :ih_npp_agdw_si_scpf, :ih_npp_stor_si_scpf, :ih_mortality_canopy_si_scpf,
    :ih_mortality_understory_si_scpf, :ih_m3_mortality_canopy_si_scpf, :ih_m3_mortality_understory_si_scpf, :ih_nplant_canopy_si_scpf, :ih_nplant_understory_si_scpf, :ih_ddbh_canopy_si_scpf,
    :ih_ddbh_understory_si_scpf, :ih_gpp_canopy_si_scpf, :ih_gpp_understory_si_scpf, :ih_ar_canopy_si_scpf, :ih_ar_understory_si_scpf, :ih_ddbh_si_scpf,
    :ih_growthflux_si_scpf, :ih_growthflux_fusion_si_scpf, :ih_ba_si_scpf, :ih_agb_si_scpf, :ih_m1_si_scpf, :ih_m2_si_scpf,
    :ih_m3_si_scpf, :ih_m4_si_scpf, :ih_m5_si_scpf, :ih_m6_si_scpf, :ih_m7_si_scpf, :ih_m8_si_scpf,
    :ih_m9_si_scpf, :ih_m10_si_scpf, :ih_m11_si_scpf, :ih_crownfiremort_si_scpf, :ih_cambialfiremort_si_scpf, :ih_abg_mortality_cflux_si_scpf,
    :ih_abg_productivity_cflux_si_scpf, :ih_m10_si_capf, :ih_nplant_si_capf, :ih_ar_si_scpf, :ih_ar_grow_si_scpf, :ih_ar_maint_si_scpf,
    :ih_ar_darkm_si_scpf, :ih_ar_agsapm_si_scpf, :ih_ar_crootm_si_scpf, :ih_ar_frootm_si_scpf, :ih_c13disc_si_scpf, :ih_ba_si_scls,
    :ih_nplant_si_scls, :ih_nplant_canopy_si_scls, :ih_nplant_understory_si_scls, :ih_lai_canopy_si_scls, :ih_lai_understory_si_scls, :ih_sai_canopy_si_scls,
    :ih_sai_understory_si_scls, :ih_mortality_canopy_si_scls, :ih_mortality_understory_si_scls, :ih_m3_mortality_canopy_si_scls, :ih_m3_mortality_understory_si_scls, :ih_demotion_rate_si_scls,
    :ih_promotion_rate_si_scls, :ih_trimming_canopy_si_scls, :ih_trimming_understory_si_scls, :ih_crown_area_canopy_si_scls, :ih_crown_area_understory_si_scls, :ih_ddbh_canopy_si_scls,
    :ih_ddbh_understory_si_scls, :ih_agb_si_scls, :ih_biomass_si_scls, :ih_mortality_canopy_secondary_si_scls, :ih_m1_si_scls, :ih_m2_si_scls,
    :ih_m3_si_scls, :ih_m4_si_scls, :ih_m5_si_scls, :ih_m6_si_scls, :ih_m7_si_scls, :ih_m8_si_scls,
    :ih_m9_si_scls, :ih_m10_si_scls, :ih_m1_sec_si_scls, :ih_m2_sec_si_scls, :ih_m3_sec_si_scls, :ih_m7_sec_si_scls,
    :ih_m8_sec_si_scls, :ih_m9_sec_si_scls, :ih_m10_sec_si_scls, :ih_m10_si_cacls, :ih_nplant_si_cacls, :ih_rdark_canopy_si_scls,
    :ih_livestem_mr_canopy_si_scls, :ih_livecroot_mr_canopy_si_scls, :ih_froot_mr_canopy_si_scls, :ih_resp_g_canopy_si_scls, :ih_resp_m_canopy_si_scls, :ih_leaf_md_canopy_si_scls,
    :ih_root_md_canopy_si_scls, :ih_carbon_balance_canopy_si_scls, :ih_bstore_md_canopy_si_scls, :ih_bdead_md_canopy_si_scls, :ih_bsw_md_canopy_si_scls, :ih_seed_prod_canopy_si_scls,
    :ih_npp_leaf_canopy_si_scls, :ih_npp_fnrt_canopy_si_scls, :ih_npp_sapw_canopy_si_scls, :ih_npp_dead_canopy_si_scls, :ih_npp_seed_canopy_si_scls, :ih_npp_stor_canopy_si_scls,
    :ih_rdark_understory_si_scls, :ih_livestem_mr_understory_si_scls, :ih_livecroot_mr_understory_si_scls, :ih_froot_mr_understory_si_scls, :ih_resp_g_understory_si_scls, :ih_resp_m_understory_si_scls,
    :ih_leaf_md_understory_si_scls, :ih_root_md_understory_si_scls, :ih_carbon_balance_understory_si_scls, :ih_bsw_md_understory_si_scls, :ih_bdead_md_understory_si_scls, :ih_bstore_md_understory_si_scls,
    :ih_seed_prod_understory_si_scls, :ih_npp_leaf_understory_si_scls, :ih_npp_fnrt_understory_si_scls, :ih_npp_sapw_understory_si_scls, :ih_npp_dead_understory_si_scls, :ih_npp_seed_understory_si_scls,
    :ih_npp_stor_understory_si_scls, :ih_yesterdaycanopylevel_canopy_si_scls, :ih_yesterdaycanopylevel_understory_si_scls, :ih_biomass_si_pft, :ih_biomass_sec_si_pft, :ih_leafbiomass_si_pft,
    :ih_storebiomass_si_pft, :ih_nindivs_si_pft, :ih_nindivs_sec_si_pft, :ih_recruitment_si_pft, :ih_recruitment_cflux_si_pft, :ih_mortality_si_pft,
    :ih_mortality_carbonflux_si_pft, :ih_hydraulicmortality_carbonflux_si_pft, :ih_cstarvmortality_carbonflux_si_pft, :ih_firemortality_carbonflux_si_pft, :ih_cstarvmortality_continuous_carbonflux_si_pft, :ih_crownarea_si_pft,
    :ih_canopycrownarea_si_pft, :ih_crownarea_si_cnlf, :ih_gpp_si_pft, :ih_gpp_sec_si_pft, :ih_npp_si_pft, :ih_npp_sec_si_pft,
    :ih_site_dstatus_si_pft, :ih_dleafoff_si_pft, :ih_dleafon_si_pft, :ih_meanliqvol_si_pft, :ih_meansmp_si_pft, :ih_elong_factor_si_pft,
    :ih_nocomp_pftpatchfraction_si_pft, :ih_nocomp_pftnpatches_si_pft, :ih_nocomp_pftburnedarea_si_pft, :ih_seeds_out_gc_si_pft, :ih_seeds_in_gc_si_pft, :ih_area_si_age,
    :ih_lai_si_age, :ih_canopy_area_si_age, :ih_gpp_si_age, :ih_npp_si_age, :ih_ncl_si_age, :ih_npatches_si_age,
    :ih_zstar_si_age, :ih_biomass_si_age, :ih_c_stomata_si_age, :ih_c_lblayer_si_age, :ih_agesince_anthrodist_si_age, :ih_secondarylands_area_si_age,
    :ih_area_burnt_si_age, :ih_fire_intensity_si_age, :ih_fire_sum_fuel_si_age, :ih_lai_secondary_si, :ih_canopy_height_dist_si_height, :ih_leaf_height_dist_si_height,
    :ih_errh2o_scpf, :ih_tran_scpf, :ih_sapflow_scpf, :ih_sapflow_si, :ih_iterh1_scpf, :ih_iterh2_scpf,
    :ih_supsub_scpf, :ih_ath_scpf, :ih_tth_scpf, :ih_sth_scpf, :ih_lth_scpf, :ih_awp_scpf,
    :ih_twp_scpf, :ih_swp_scpf, :ih_lwp_scpf, :ih_aflc_scpf, :ih_tflc_scpf, :ih_sflc_scpf,
    :ih_lflc_scpf, :ih_btran_scpf, :ih_rootwgt_soilvwc_si, :ih_rootwgt_soilvwcsat_si, :ih_rootwgt_soilmatpot_si, :ih_soilmatpot_sl,
    :ih_soilvwc_sl, :ih_soilvwcsat_sl, :ih_rootuptake_si, :ih_rootuptake_sl, :ih_rootuptake0_scpf, :ih_rootuptake10_scpf,
    :ih_rootuptake50_scpf, :ih_rootuptake100_scpf, :ih_litter_moisture_si_fuel, :ih_burnt_frac_litter_si_fuel, :ih_fuel_amount_si_fuel, :ih_cwd_ag_si_cwdsc,
    :ih_cwd_bg_si_cwdsc, :ih_cwd_ag_in_si_cwdsc, :ih_cwd_bg_in_si_cwdsc, :ih_cwd_ag_out_si_cwdsc, :ih_cwd_bg_out_si_cwdsc, :ih_parsun_z_si_cnlf,
    :ih_parsha_z_si_cnlf, :ih_laisun_z_si_cnlf, :ih_laisha_z_si_cnlf, :ih_ts_net_uptake_si_cnlf, :ih_crownarea_clll, :ih_parprof_dir_si_cnlf,
    :ih_parprof_dif_si_cnlf, :ih_parsun_z_si_cnlfpft, :ih_parsha_z_si_cnlfpft, :ih_laisun_clllpf, :ih_laisha_clllpf, :ih_parprof_dir_si_cnlfpft,
    :ih_parprof_dif_si_cnlfpft, :ih_crownfrac_clllpf, :ih_nplant_si_cdpf, :ih_nplant_canopy_si_cdpf, :ih_nplant_understory_si_cdpf, :ih_mortality_si_cdpf,
    :ih_mortality_canopy_si_cdpf, :ih_mortality_understory_si_cdpf, :ih_m3_si_cdpf, :ih_m11_si_cdpf, :ih_m3_mortality_canopy_si_cdpf, :ih_m3_mortality_understory_si_cdpf,
    :ih_m11_mortality_canopy_si_cdpf, :ih_m11_mortality_understory_si_cdpf, :ih_ddbh_si_cdpf, :ih_ddbh_canopy_si_cdpf, :ih_ddbh_understory_si_cdpf, :ih_crownarea_canopy_damage_si,
    :ih_crownarea_ustory_damage_si, :ih_parsun_si_can, :ih_parsha_si_can, :ih_laisun_si_can, :ih_laisha_si_can, :ih_crownarea_cl,
    :ih_fuel_amount_age_fuel,
]

# Auto-derived registration table from FATES define_history_vars (faithful 1:1, 472 vars).
# Tuple fields: (gates::Vector{Symbol}, vname, units, long, use_default, avgflag,
#                vtype, hlms, upfreq, handle_field::Symbol, flush_to_zero::Bool)
# `gates` are ANDed; each is checked by `_gate_passes` against the live hlm flags.
const _HISTORY_VAR_REGISTRY = Tuple{Vector{Symbol},String,String,String,String,String,String,String,Int,Symbol,Bool}[
    (Symbol[:dynam0], "FATES_NPATCHES", "", "total number of patches per site", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_npatches_si, false),
    (Symbol[:dynam0], "FATES_NCOHORTS", "", "total number of cohorts per site", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_ncohorts_si, false),
    (Symbol[:dynam0], "FATES_NPATCHES_SECONDARY", "", "total number of patches per site", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_npatches_sec_si, false),
    (Symbol[:dynam0], "FATES_NCOHORTS_SECONDARY", "", "total number of cohorts per site", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_ncohorts_sec_si, false),
    (Symbol[:dynam0], "FATES_TRIMMING", "1", "degree to which canopy expansion is limited by leaf economics (0-1)", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_trimming_si, false),
    (Symbol[:dynam0], "FATES_AREA_PLANTS", "m2 m-2", "area occupied by all plants per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_area_plant_si, false),
    (Symbol[:dynam0], "FATES_AREA_TREES", "m2 m-2", "area occupied by woody plants per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_area_trees_si, false),
    (Symbol[:dynam0], "FATES_FRACTION", "m2 m-2", "total gridcell fraction which FATES is running over", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fates_fraction_si, true),
    (Symbol[:dynam0], "FATES_BA_WEIGHTED_HEIGHT", "m", "basal area-weighted mean height of woody plants", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_ba_weighted_height_si, false),
    (Symbol[:dynam0], "FATES_CA_WEIGHTED_HEIGHT", "m", "crown area-weighted mean height of canopy plants", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_ca_weighted_height_si, false),
    (Symbol[:dynam0], "FATES_COLD_STATUS", "", "'site-level cold status", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_site_cstatus_si, false),
    (Symbol[:dynam0], "FATES_GDD", "degree_Celsius", "site-level growing degree days", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_gdd_si, false),
    (Symbol[:dynam0], "FATES_NCHILLDAYS", "days", "site-level number of chill days", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_site_nchilldays_si, false),
    (Symbol[:dynam0], "FATES_NCOLDDAYS", "days", "site-level number of cold days", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_site_ncolddays_si, false),
    (Symbol[:dynam0], "FATES_DAYSINCE_COLDLEAFOFF", "days", "site-level days elapsed since cold leaf drop", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_cleafoff_si, false),
    (Symbol[:dynam0], "FATES_DAYSINCE_COLDLEAFON", "days", "site-level days elapsed since cold leaf flush", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_cleafon_si, false),
    (Symbol[:dynam0], "FATES_CANOPY_SPREAD", "", "scaling factor (0-1) between tree basal area and canopy area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_canopy_spread_si, false),
    (Symbol[:dynam0], "FATES_LAI", "m2 m-2", "total leaf area index per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_lai_si, false),
    (Symbol[:dynam0], "FATES_ELAI", "m2 m-2", "exposed (non snow-occluded) leaf area index per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_elai_si, false),
    (Symbol[:dynam0], "FATES_LAI_SECONDARY", "m2 m-2", "'leaf area index per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_lai_secondary_si, false),
    (Symbol[:dynam0], "FATES_SECONDARY_FOREST_FRACTION", "m2 m-2", "secondary forest fraction", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fraction_secondary_forest_si, false),
    (Symbol[:dynam0], "FATES_WOOD_PRODUCT", "kg m-2", "total wood product from logging in kg carbon per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_woodproduct_si, false),
    (Symbol[:dynam0], "FATES_SECONDARY_FOREST_VEGC", "kg m-2", "biomass on secondary lands in kg carbon per m2 land area (mult by FATES_SECONDARY_FOREST_FRACTION to get per secondary forest area)", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_biomass_secondary_forest_si, false),
    (Symbol[:dynam0], "FATES_NESTEROV_INDEX", "", "nesterov fire danger index", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_nesterov_fire_danger_si, false),
    (Symbol[:dynam0], "FATES_IGNITIONS", "m-2 s-1", "number of successful fire ignitions per m2 land area per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fire_nignitions_si, false),
    (Symbol[:dynam0], "FATES_FDI", "1", "Fire Danger Index (probability that an ignition will lead to a fire)", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fire_fdi_si, false),
    (Symbol[:dynam0], "FATES_ROS", "m s-1", "fire rate of spread in meters per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_spitfire_ros_si, false),
    (Symbol[:dynam0], "FATES_EFFECT_WSPEED", "m s-1", "effective wind speed for fire spread in meters per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_effect_wspeed_si, false),
    (Symbol[:dynam0], "FATES_FUELCONSUMED", "kg m-2", "total fuel consumed in kg carbon per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_tfc_ros_si, false),
    (Symbol[:dynam0], "FATES_FIRE_INTENSITY", "J m-1 s-1", "spitfire surface fireline intensity in J per m per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fire_intensity_si, false),
    (Symbol[:dynam0], "FATES_FIRE_INTENSITY_BURNFRAC", "J m-1 s-1", "product of surface fire intensity and burned area fraction -- divide by FATES_BURNFRAC to get area-weighted mean intensity", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fire_intensity_area_product_si, false),
    (Symbol[:dynam0], "FATES_BURNFRAC", "s-1", "burned area fraction per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fire_area_si, false),
    (Symbol[:dynam0], "FATES_FUEL_MEF", "m3 m-3", "fuel moisture of extinction (volumetric)", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fire_fuel_mef_si, false),
    (Symbol[:dynam0], "FATES_FUEL_BULKD", "kg m-3", "fuel bulk density in kg per m3", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fire_fuel_bulkd_si, false),
    (Symbol[:dynam0], "FATES_FUEL_EFF_MOIST", "m3 m-3", "spitfire fuel moisture (volumetric)", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fire_fuel_eff_moist_si, false),
    (Symbol[:dynam0], "FATES_FUEL_SAV", "m-1", "spitfire fuel surface area to volume ratio", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fire_fuel_sav_si, false),
    (Symbol[:dynam0], "FATES_FUEL_AMOUNT", "kg m-2", "total ground fuel related to FATES_ROS (omits 1000hr fuels) in kg C per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_sum_fuel_si, false),
    (Symbol[:dynam0], "FATES_LITTER_IN", "kg m-2 s-1", "litter flux in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_litter_in_si, false),
    (Symbol[:dynam0], "FATES_LITTER_OUT", "kg m-2 s-1", "litter flux out in kg carbon (exudation, fragmentation, seed decay)", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_litter_out_si, false),
    (Symbol[:dynam0], "FATES_SEED_BANK", "kg m-2", "total seed mass of all PFTs in kg carbon per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_seed_bank_si, false),
    (Symbol[:dynam0], "FATES_UNGERM_SEED_BANK", "kg m-2", "ungerminated seed mass of all PFTs in kg carbon per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_ungerm_seed_bank_si, false),
    (Symbol[:dynam0], "FATES_SEEDLING_POOL", "kg m-2", "total seedling (ie germinated seeds) mass of all PFTs in kg carbon per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_seedling_pool_si, false),
    (Symbol[:dynam0], "FATES_SEEDS_IN", "kg m-2 s-1", "seed production rate in kg carbon per m2 second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_seeds_in_si, false),
    (Symbol[:dynam0], "FATES_SEEDS_IN_LOCAL", "kg m-2 s-1", "local seed production rate in kg carbon per m2 second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_seeds_in_local_si, false),
    (Symbol[:dynam0], "FATES_STOREC", "kg m-2", "total biomass in live plant storage in kg carbon per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_storec_si, false),
    (Symbol[:dynam0], "FATES_STOREC_TF", "kg kg-1", "Storage C fraction of target", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_storectfrac_si, false),
    (Symbol[:dynam0], "FATES_VEGC", "kg m-2", "total biomass in live plants in kg carbon per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_totvegc_si, false),
    (Symbol[:dynam0], "FATES_SAPWOODC", "kg m-2", "total biomass in live plant sapwood in kg carbon per m2", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_sapwc_si, false),
    (Symbol[:dynam0], "FATES_LEAFC", "kg m-2", "total biomass in live plant leaves in kg carbon per m2", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_leafc_si, false),
    (Symbol[:dynam0], "FATES_FROOTC", "kg m-2", "total biomass in live plant fine roots in kg carbon per m2", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fnrtc_si, false),
    (Symbol[:dynam0], "FATES_REPROC", "kg m-2", "total biomass in live plant reproductive tissues in kg carbon per m2", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_reproc_si, false),
    (Symbol[:dynam0], "FATES_L2FR", "kg kg-1", "The leaf to fineroot biomass multiplier for target allometry", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_l2fr_si, false),
    (Symbol[:dynam0, :has_n], "FATES_NH4UPTAKE", "kg m-2 s-1", "ammonium uptake rate by plants in kg NH4 per m2 per second", "active", "A", site_r8, "CLM:ALM", group_nflx_simple, :ih_nh4uptake_si, false),
    (Symbol[:dynam0, :has_n], "FATES_NO3UPTAKE", "kg m-2 s-1", "nitrate uptake rate by plants in kg NO3 per m2 per second", "active", "A", site_r8, "CLM:ALM", group_nflx_simple, :ih_no3uptake_si, false),
    (Symbol[:dynam0, :has_n], "FATES_NEFFLUX", "kg m-2 s-1", "nitrogen effluxed from plant in kg N per m2 per second (unused)", "active", "A", site_r8, "CLM:ALM", group_nflx_simple, :ih_nefflux_si, false),
    (Symbol[:dynam0, :has_n], "FATES_NDEMAND", "kg m-2 s-1", "plant nitrogen need (algorithm dependent) in kg N per m2 per second", "active", "A", site_r8, "CLM:ALM", group_nflx_simple, :ih_ndemand_si, false),
    (Symbol[:dynam0, :has_n], "FATES_NFIX_SYM", "kg m-2 s-1", "symbiotic dinitrogen fixation in kg N per m2 per second", "active", "A", site_r8, "CLM:ALM", group_nflx_simple, :ih_nfix_si, false),
    (Symbol[:dynam0, :has_n], "FATES_STOREN", "kg m-2", "total nitrogen in live plant storage", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_storen_si, false),
    (Symbol[:dynam0, :has_n], "FATES_STOREN_TF", "1", "storage N fraction of target", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_storentfrac_si, false),
    (Symbol[:dynam0, :has_n], "FATES_VEGN", "kg m-2", "total nitrogen in live plants", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_totvegn_si, false),
    (Symbol[:dynam0, :has_n], "FATES_SAPWOODN", "kg m-2", "total nitrogen in live plant sapwood", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_sapwn_si, false),
    (Symbol[:dynam0, :has_n], "FATES_LEAFN", "kg m-2", "total nitrogen in live plant leaves", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_leafn_si, false),
    (Symbol[:dynam0, :has_n], "FATES_FROOTN", "kg m-2", "total nitrogen in live plant fine-roots", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fnrtn_si, false),
    (Symbol[:dynam0, :has_n], "FATES_REPRON", "kg m-2", "total nitrogen in live plant reproductive tissues", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_repron_si, false),
    (Symbol[:dynam0, :has_p], "FATES_STOREP", "kg m-2", "total phosphorus in live plant storage", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_storep_si, false),
    (Symbol[:dynam0, :has_p], "FATES_STOREP_TF", "1", "storage P fraction of target", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_storeptfrac_si, false),
    (Symbol[:dynam0, :has_p], "FATES_VEGP", "kg m-2", "total phosphorus in live plants", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_totvegp_si, false),
    (Symbol[:dynam0, :has_p], "FATES_SAPWOODP", "kg m-2", "Total phosphorus in live plant sapwood", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_sapwp_si, false),
    (Symbol[:dynam0, :has_p], "FATES_LEAFP", "kg m-2", "total phosphorus in live plant leaves", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_leafp_si, false),
    (Symbol[:dynam0, :has_p], "FATES_FROOTP", "kg m-2", "total phosphorus in live plant fine roots", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fnrtp_si, false),
    (Symbol[:dynam0, :has_p], "FATES_REPROP", "kg m-2", "total phosphorus in live plant reproductive tissues", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_reprop_si, false),
    (Symbol[:dynam0, :has_p], "FATES_PUPTAKE", "kg m-2 s-1", "mineralized phosphorus uptake rate of plants in kg P per m2 per second", "active", "A", site_r8, "CLM:ALM", group_nflx_simple, :ih_puptake_si, false),
    (Symbol[:dynam0, :has_p], "FATES_PEFFLUX", "kg m-2 s-1", "phosphorus effluxed from plant in kg P per m2 per second (unused)", "active", "A", site_r8, "CLM:ALM", group_nflx_simple, :ih_pefflux_si, false),
    (Symbol[:dynam0, :has_p], "FATES_PDEMAND", "kg m-2 s-1", "plant phosphorus need (algorithm dependent) in kg P per m2 per second", "active", "A", site_r8, "CLM:ALM", group_nflx_simple, :ih_pdemand_si, false),
    (Symbol[:dynam0], "FATES_STRUCTC", "kg m-2", "structural biomass in kg carbon per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_bdead_si, false),
    (Symbol[:dynam0], "FATES_NONSTRUCTC", "kg m-2", "non-structural biomass (sapwood + leaf + fineroot) in kg carbon per m2", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_balive_si, false),
    (Symbol[:dynam0], "FATES_VEGC_ABOVEGROUND", "kg m-2", "aboveground biomass in kg carbon per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_agb_si, false),
    (Symbol[:dynam0], "FATES_CANOPY_VEGC", "kg m-2", "biomass of canopy plants in kg carbon per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_canopy_biomass_si, false),
    (Symbol[:dynam0], "FATES_USTORY_VEGC", "kg m-2", "biomass of understory plants in kg carbon per m2 land area", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_understory_biomass_si, false),
    (Symbol[:dynam0], "FATES_PRIMARY_PATCHFUSION_ERR", "m2 m-2 yr-1", "error in total primary lands associated with patch fusion", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_primaryland_fusion_error_si, false),
    (Symbol[:dynam0], "FATES_DISTURBANCE_RATE_FIRE", "m2 m-2 yr-1", "disturbance rate from fire", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fire_disturbance_rate_si, false),
    (Symbol[:dynam0], "FATES_DISTURBANCE_RATE_LOGGING", "m2 m-2 yr-1", "disturbance rate from logging", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_logging_disturbance_rate_si, false),
    (Symbol[:dynam0], "FATES_DISTURBANCE_RATE_TREEFALL", "m2 m-2 yr-1", "disturbance rate from treefall", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fall_disturbance_rate_si, false),
    (Symbol[:dynam0], "FATES_HARVEST_WOODPROD_C_FLUX", "kg m-2 yr-1", "harvest-associated wood product carbon flux in kg C per m2 per year", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_harvest_woodprod_carbonflux_si, false),
    (Symbol[:dynam0], "FATES_LUCHANGE_WOODPROD_C_FLUX", "kg m-2 yr-1", "land-use-change-associated wood product carbon flux in kg C per m2 per year", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_luchange_woodprod_carbonflux_si, false),
    (Symbol[:dynam0], "FATES_TVEG24", "degree_Celsius", "fates 24-hr running mean vegetation temperature by site", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_tveg24_si, false),
    (Symbol[:dynam0], "FATES_TLONGTERM", "degree_Celsius", "fates 30-year running mean vegetation temperature by site", "inactive", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_tlongterm_si, false),
    (Symbol[:dynam0], "FATES_TGROWTH", "degree_Celsius", "fates long-term running mean vegetation temperature by site", "inactive", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_tgrowth_si, false),
    (Symbol[:dynam0], "FATES_HARVEST_DEBT", "kg C", "Accumulated carbon failed to be harvested", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_harvest_debt_si, false),
    (Symbol[:dynam0], "FATES_HARVEST_DEBT_SEC", "kg C", "Accumulated carbon failed to be harvested from secondary patches", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_harvest_debt_sec_si, false),
    (Symbol[:dynam0], "FATES_EXCESS_RESP", "kg m-2 s-1", "respiration of un-allocatable carbon gain", "active", "A", site_r8, "CLM:ALM", group_nflx_simple, :ih_excess_resp_si, false),
    (Symbol[:dynam0], "FATES_DEMOTION_CARBONFLUX", "kg m-2 s-1", "demotion-associated biomass carbon flux from canopy to understory in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_demotion_carbonflux_si, false),
    (Symbol[:dynam0], "FATES_PROMOTION_CARBONFLUX", "kg m-2 s-1", "promotion-associated biomass carbon flux from understory to canopy in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_promotion_carbonflux_si, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CFLUX_CANOPY", "kg m-2 s-1", "flux of biomass carbon from live to dead pools from mortality of canopy plants in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_canopy_mortality_carbonflux_si, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CFLUX_USTORY", "kg m-2 s-1", "flux of biomass carbon from live to dead pools from mortality of understory plants in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_understory_mortality_carbonflux_si, false),
    (Symbol[:dynam0], "MORTALITY_CROWNAREA_CANOPY", "m2/ha/year", "Crown area of canopy trees that died", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_canopy_mortality_crownarea_si, false),
    (Symbol[:dynam0], "MORTALITY_CROWNAREA_UNDERSTORY", "m2/ha/year", "Crown aera of understory trees that died", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_understory_mortality_crownarea_si, false),
    (Symbol[:dynam0], "FATES_FIRE_CLOSS", "kg m-2 s-1", "carbon loss to atmosphere from fire in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_fire_c_to_atm_si, false),
    (Symbol[:dynam0], "FATES_CBALANCE_ERROR", "kg s-1", "total carbon error in kg carbon per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_cbal_err_fates_si, false),
    (Symbol[:dynam0], "FATES_LEAF_ALLOC", "kg m-2 s-1", "allocation to leaves in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_npp_leaf_si, false),
    (Symbol[:dynam0], "FATES_SEED_ALLOC", "kg m-2 s-1", "allocation to seeds in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_npp_seed_si, false),
    (Symbol[:dynam0], "FATES_STEM_ALLOC", "kg m-2 s-1", "allocation to stem in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_npp_stem_si, false),
    (Symbol[:dynam0], "FATES_FROOT_ALLOC", "kg m-2 s-1", "allocation to fine roots in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_npp_froot_si, false),
    (Symbol[:dynam0], "FATES_CROOT_ALLOC", "kg m-2 s-1", "allocation to coarse roots in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_npp_croot_si, false),
    (Symbol[:dynam0], "FATES_STORE_ALLOC", "kg m-2 s-1", "allocation to storage tissues in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_npp_stor_si, false),
    (Symbol[:dynam0, :planthydro], "FATES_VEGH2O_DEAD", "kg m-2", "cumulative water stored in dead biomass due to mortality", "inactive", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_h2oveg_dead_si, false),
    (Symbol[:dynam0, :planthydro], "FATES_VEGH2O_RECRUIT", "kg m-2", "amount of water in new recruits", "inactive", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_h2oveg_recruit_si, false),
    (Symbol[:dynam0, :planthydro], "FATES_VEGH2O_GROWTURN_ERR", "kg m-2", "cumulative net borrowed (+) or lost (-) from water storage due to combined growth turnover", "inactive", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_h2oveg_growturn_err_si, false),
    (Symbol[:dynam0, :treedamage], "FATES_CROWNAREA_CANOPY_CD", "m2 m-2 yr-1", "crownarea lost to damage each year", "inactive", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_crownarea_canopy_damage_si, false),
    (Symbol[:dynam0, :treedamage], "FATES_CROWNAREA_USTORY_CD", "m2 m-2 yr-1", "crownarea lost to damage each year", "inactive", "A", site_r8, "CLM:ALM", group_dyna_simple, :ih_crownarea_ustory_damage_si, false),
    (Symbol[:dynam0, :dynam1], "FATES_PATCHAREA_LU", "m2 m-2", "patch area by land use type", "active", "A", site_landuse_r8, "CLM:ALM", group_dyna_complx, :ih_area_si_landuse, false),
    (Symbol[:dynam0, :dynam1], "FATES_TRANSITION_MATRIX_LULU", "m2 m-2 yr-1", "land use transition matrix", "active", "A", site_lulu_r8, "CLM:ALM", group_dyna_complx, :ih_transition_matrix_si_lulu, false),
    (Symbol[:dynam0, :dynam1], "FATES_DISTURBANCE_RATE_MATRIX_LULU", "m2 m-2 yr-1", "disturbance rates by land use type x land use type matrix", "active", "A", site_lulu_r8, "CLM:ALM", group_dyna_complx, :ih_disturbance_rate_si_lulu, false),
    (Symbol[:dynam0, :dynam1], "FATES_VEGC_PF", "kg m-2", "total PFT-level biomass in kg of carbon per land area", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_biomass_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_VEGC_SE_PF", "kg m-2", "'total PFT-level biomass in kg of carbon per land area", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_biomass_sec_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_RECRUITMENT_CFLUX_PF", "kg m-2 yr-1", "total PFT-level biomass of new recruits in kg of carbon per land area", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_simple, :ih_recruitment_cflux_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_LEAFC_PF", "kg m-2", "total PFT-level leaf biomass in kg carbon per m2 land area", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_leafbiomass_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_STOREC_PF", "kg m-2", "total PFT-level stored biomass in kg carbon per m2 land area", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_storebiomass_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_CROWNAREA_PF", "m2 m-2", "total PFT-level crown area per m2 land area", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_crownarea_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_CANOPYCROWNAREA_PF", "m2 m-2", "total PFT-level canopy-layer crown area per m2 land area", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_canopycrownarea_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_GPP_PF", "kg m-2 s-1", "total PFT-level GPP in kg carbon per m2 land area per second", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_gpp_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_NPP_PF", "kg m-2 s-1", "total PFT-level NPP in kg carbon per m2 land area per second", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_GPP_SE_PF", "kg m-2 s-1", "'total PFT-level GPP in kg carbon per m2 land area per second", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_gpp_sec_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_NPP_SE_PF", "kg m-2 s-1", "'total PFT-level NPP in kg carbon per m2 land area per second", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_sec_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_NPLANT_PF", "m-2", "total PFT-level number of individuals per m2 land area", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_nindivs_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_NPLANT_SEC_PF", "m-2", "'total PFT-level number of individuals per m2 land area", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_nindivs_sec_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_RECRUITMENT_PF", "m-2 yr-1", "PFT-level recruitment rate in number of individuals per m2 land area per year", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_recruitment_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_SEEDS_IN_GRIDCELL_PF", "kg", "Site-level seed mass input from neighboring gridcells per pft", "inactive", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_seeds_in_gc_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_SEEDS_OUT_GRIDCELL_PF", "kg", "Site-level seed mass output to neighboring gridcells per pft", "inactive", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_seeds_out_gc_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_MORTALITY_PF", "m-2 yr-1", "PFT-level mortality rate in number of individuals per m2 land area per year", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_DROUGHT_STATUS_PF", "", "'PFT-level drought status", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_site_dstatus_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_DAYSINCE_DROUGHTLEAFOFF_PF", "days", "PFT-level days elapsed since drought leaf drop", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_dleafoff_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_DAYSINCE_DROUGHTLEAFON_PF", "days", "PFT-level days elapsed since drought leaf flush", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_dleafon_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_MEANLIQVOL_DROUGHTPHEN_PF", "m3 m-3", "PFT-level mean liquid water volume for drought phenolgy", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_meanliqvol_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_MEANSMP_DROUGHTPHEN_PF", "Pa", "PFT-level mean soil matric potential for drought phenology", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_meansmp_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_ELONG_FACTOR_PF", "1", "PFT-level mean elongation factor (partial flushing/abscission)", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_elong_factor_si_pft, false),
    (Symbol[:dynam0, :dynam1, :nocomp], "FATES_NOCOMP_NPATCHES_PF", "", "number of patches per PFT (nocomp-mode-only)", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_nocomp_pftnpatches_si_pft, false),
    (Symbol[:dynam0, :dynam1, :nocomp], "FATES_NOCOMP_PATCHAREA_PF", "m2 m-2", "total patch area allowed per PFT (nocomp-mode-only)", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_nocomp_pftpatchfraction_si_pft, false),
    (Symbol[:dynam0, :dynam1, :nocomp], "FATES_NOCOMP_BURNEDAREA_PF", "s-1", "total burned area of PFT-labeled patch area (nocomp-mode-only)", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_nocomp_pftburnedarea_si_pft, false),
    (Symbol[:dynam0, :dynam1], "FATES_PATCHAREA_AP", "m2 m-2", "patch area by age bin per m2 land area", "active", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_area_si_age, false),
    (Symbol[:dynam0, :dynam1], "FATES_LAI_AP", "m2 m-2", "total leaf area index by age bin per m2 land area", "active", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_lai_si_age, false),
    (Symbol[:dynam0, :dynam1], "FATES_CANOPYAREA_AP", "m2 m-2", "canopy area by age bin per m2 land area", "active", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_canopy_area_si_age, false),
    (Symbol[:dynam0, :dynam1], "FATES_NCL_AP", "", "number of canopy levels by age bin", "inactive", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_ncl_si_age, false),
    (Symbol[:dynam0, :dynam1], "FATES_NPATCH_AP", "", "number of patches by age bin", "inactive", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_npatches_si_age, false),
    (Symbol[:dynam0], "FATES_ZSTAR_AP", "m", "product of zstar and patch area by age bin (divide by FATES_PATCHAREA_AP to get mean zstar)", "trim(tempstring)", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_zstar_si_age, false),
    (Symbol[:dynam0], "FATES_CANOPYAREA_HT", "m2 m-2", "canopy area height distribution", "active", "A", site_height_r8, "CLM:ALM", group_dyna_complx, :ih_canopy_height_dist_si_height, false),
    (Symbol[:dynam0], "FATES_LEAFAREA_HT", "m2 m-2", "leaf area height distribution", "active", "A", site_height_r8, "CLM:ALM", group_dyna_complx, :ih_leaf_height_dist_si_height, false),
    (Symbol[:dynam0], "FATES_VEGC_AP", "kg m-2", "total biomass within a given patch age bin in kg carbon per m2 land area", "inactive", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_biomass_si_age, false),
    (Symbol[:dynam0], "FATES_SECONDAREA_ANTHRODIST_AP", "m2 m-2", "secondary forest patch area age distribution since anthropgenic disturbance", "inactive", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_agesince_anthrodist_si_age, false),
    (Symbol[:dynam0], "FATES_SECONDAREA_DIST_AP", "m2 m-2", "secondary forest patch area age distribution since any kind of disturbance", "inactive", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_secondarylands_area_si_age, false),
    (Symbol[:dynam0], "FATES_FRAGMENTATION_SCALER_SL", "", "factor (0-1) by which litter/cwd fragmentation proceeds relative to max rate by soil layer", "active", "A", site_soil_r8, "CLM:ALM", group_dyna_complx, :ih_fragmentation_scaler_sl, false),
    (Symbol[:dynam0], "FATES_FUEL_MOISTURE_FC", "m3 m-3", "spitfire fuel class-level fuel moisture (volumetric)", "active", "A", site_fuel_r8, "CLM:ALM", group_dyna_complx, :ih_litter_moisture_si_fuel, false),
    (Symbol[:dynam0], "FATES_FUEL_AMOUNT_FC", "kg m-2", "spitfire fuel-class level fuel amount in kg carbon per m2 land area", "active", "A", site_fuel_r8, "CLM:ALM", group_dyna_complx, :ih_fuel_amount_si_fuel, false),
    (Symbol[:dynam0], "FATES_FUEL_AMOUNT_APFC", "kg m-2", "spitfire fuel quantity in each age x fuel class in kg carbon per m2 land area", "inactive", "A", site_agefuel_r8, "CLM:ALM", group_dyna_complx, :ih_fuel_amount_age_fuel, false),
    (Symbol[:dynam0], "FATES_BURNFRAC_AP", "s-1", "spitfire fraction area burnt (per second) by patch age", "active", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_area_burnt_si_age, false),
    (Symbol[:dynam0], "FATES_FIRE_INTENSITY_BURNFRAC_AP", "J m-1 s-1", "'product of fire intensity and burned fraction", "active", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_fire_intensity_si_age, false),
    (Symbol[:dynam0], "FATES_FUEL_AMOUNT_AP", "kg m-2", "spitfire ground fuel (kg carbon per m2) related to FATES_ROS (omits 1000hr fuels) within each patch age bin (divide by FATES_PATCHAREA_AP to get fuel per unit area of that-age patch)", "active", "A", site_age_r8, "CLM:ALM", group_dyna_complx, :ih_fire_sum_fuel_si_age, false),
    (Symbol[:dynam0], "FATES_FUEL_BURNT_BURNFRAC_FC", "1", "product of fraction (0-1) of fuel burnt and burnt fraction (divide by FATES_BURNFRAC to get burned-area-weighted mean fraction fuel burnt)", "active", "A", site_fuel_r8, "CLM:ALM", group_dyna_complx, :ih_burnt_frac_litter_si_fuel, false),
    (Symbol[:dynam0], "FATES_LITTER_IN_EL", "kg m-2 s-1", "litter flux in in kg element per m2 per second", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_litter_in_elem, false),
    (Symbol[:dynam0], "FATES_LITTER_OUT_EL", "kg m-2 s-1", "litter flux out (exudation, fragmentation and seed decay) in kg element", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_litter_out_elem, false),
    (Symbol[:dynam0], "FATES_SEED_BANK_EL", "kg m-2", "element-level total seed mass of all PFTs in kg element per m2", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_seed_bank_elem, false),
    (Symbol[:dynam0], "FATES_SEEDS_IN_LOCAL_EL", "kg m-2 s-1", "'within-site", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_seeds_in_local_elem, false),
    (Symbol[:dynam0], "FATES_SEEDS_IN_EXTERN_EL", "kg m-2 s-1", "external seed influx rate in kg element per m2 per second", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_seeds_in_extern_elem, false),
    (Symbol[:dynam0], "FATES_SEED_GERM_EL", "kg m-2", "element-level total germinated seed mass of all PFTs in kg element per m2", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_seed_germ_elem, false),
    (Symbol[:dynam0], "FATES_SEED_DECAY_EL", "kg m-2 s-1", "seed mass decay (germinated and un-germinated) in kg element per m2 per second", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_seed_decay_elem, false),
    (Symbol[:dynam0], "FATES_STOREC_TF_USTORY_SZPF", "kg kg-1", "'Storage C fraction of target by size x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_storectfrac_ustory_scpf, false),
    (Symbol[:dynam0], "FATES_STOREC_TF_CANOPY_SZPF", "kg kg-1", "'Storage C fraction of target by size x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_storectfrac_canopy_scpf, false),
    (Symbol[:dynam0], "FATES_FROOTC_SL", "kg m-3", "Total carbon in live plant fine-roots over depth", "active", "A", site_soil_r8, "CLM:ALM", group_dyna_complx, :ih_fnrtc_sl, false),
    (Symbol[:dynam0], "FATES_L2FR_CANOPY_REC_PF", "kg kg-1", "The leaf to fineroot biomass multiplier for recruits (canopy)", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_recl2fr_canopy_pf, false),
    (Symbol[:dynam0], "FATES_L2FR_USTORY_REC_PF", "kg kg-1", "The leaf to fineroot biomass multiplier for recruits (understory)", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_recl2fr_ustory_pf, false),
    (Symbol[:dynam0, :has_n], "FATES_NH4UPTAKE_SZPF", "kg m-2 s-1", "ammonium uptake rate by plants by size-class x pft in kg NH4 per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_nflx_complx, :ih_nh4uptake_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_NO3UPTAKE_SZPF", "kg m-2 s-1", "nitrate uptake rate by plants by size-class x pft in kg NO3 per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_nflx_complx, :ih_no3uptake_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_NEFFLUX_SZPF", "kg m-2 s-1", "'nitrogen efflux", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_nflx_complx, :ih_nefflux_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_NDEMAND_SZPF", "kg m-2 s-1", "'plant N need (algorithm dependent)", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_nflx_complx, :ih_ndemand_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_NFIX_SYM_SZPF", "kg m-2 s-1", "'symbiotic dinitrogen fixation", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_nflx_complx, :ih_nfix_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_VEGN_SZPF", "kg m-2", "total (live) vegetation nitrogen mass by size-class x pft in kg N per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_totvegn_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_LEAFN_SZPF", "kg m-2", "leaf nitrogen mass by size-class x pft in kg N per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_leafn_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_FROOTN_SZPF", "kg m-2", "fine-root nitrogen mass by size-class x pft in kg N per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_fnrtn_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_SAPWOODN_SZPF", "kg m-2", "sapwood nitrogen mass by size-class x pft in kg N per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_sapwn_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_STOREN_SZPF", "kg m-2", "storage nitrogen mass by size-class x pft in kg N per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_storen_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_STOREN_TF_CANOPY_SZPF", "1", "'storage nitrogen fraction (0-1) of target", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_storentfrac_canopy_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_STOREN_TF_USTORY_SZPF", "1", "'storage nitrogen fraction (0-1) of target", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_storentfrac_understory_scpf, false),
    (Symbol[:dynam0, :has_n], "FATES_REPRON_SZPF", "kg m-2", "reproductive nitrogen mass (on plant) by size-class x pft in kg N per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_repron_scpf, false),
    (Symbol[:dynam0, :has_p], "FATES_VEGP_SZPF", "kg m-2", "total (live) vegetation phosphorus mass by size-class x pft in kg P per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_totvegp_scpf, false),
    (Symbol[:dynam0, :has_p], "FATES_LEAFP_SZPF", "kg m-2", "leaf phosphorus mass by size-class x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_leafp_scpf, false),
    (Symbol[:dynam0, :has_p], "FATES_FROOTP_SZPF", "kg m-2", "fine-root phosphorus mass by size-class x pft in kg P per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_fnrtp_scpf, false),
    (Symbol[:dynam0, :has_p], "FATES_SAPWOODP_SZPF", "kg m-2", "sapwood phosphorus mass by size-class x pft in kg P per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_sapwp_scpf, false),
    (Symbol[:dynam0, :has_p], "FATES_STOREP_SZPF", "kg m-2", "storage phosphorus mass by size-class x pft in kg P per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_storep_scpf, false),
    (Symbol[:dynam0, :has_p], "FATES_STOREP_TF_CANOPY_SZPF", "1", "'storage phosphorus fraction (0-1) of target", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_storeptfrac_canopy_scpf, false),
    (Symbol[:dynam0, :has_p], "FATES_STOREP_TF_USTORY_SZPF", "1", "'storage phosphorus fraction (0-1) of target", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_storeptfrac_understory_scpf, false),
    (Symbol[:dynam0, :has_p], "FATES_REPROP_SZPF", "kg m-2", "reproductive phosphorus mass (on plant) by size-class x pft in kg P per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_reprop_scpf, false),
    (Symbol[:dynam0, :has_p], "FATES_PUPTAKE_SZPF", "kg m-2 s-1", "'phosphorus uptake rate by plants", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_nflx_complx, :ih_puptake_scpf, false),
    (Symbol[:dynam0, :has_p], "FATES_PEFFLUX_SZPF", "kg m-2 s-1", "'phosphorus efflux", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_nflx_complx, :ih_pefflux_scpf, false),
    (Symbol[:dynam0, :has_p], "FATES_PDEMAND_SZPF", "kg m-2 s-1", "'plant P need (algorithm dependent)", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_nflx_complx, :ih_pdemand_scpf, false),
    (Symbol[:dynam0], "FATES_CROWNAREA_CLLL", "m2 m-2", "area fraction of the total ground occupied by each canopy-leaf layer", "inactive", "A", site_cnlf_r8, "CLM:ALM", group_dyna_complx, :ih_crownarea_si_cnlf, false),
    (Symbol[:dynam0], "FATES_CROWNAREA_CL", "m2 m-2", "area fraction of the canopy footprint occupied by each canopy-leaf layer", "active", "A", site_can_r8, "CLM:ALM", group_dyna_complx, :ih_crownarea_cl, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CFLUX_PF", "kg m-2 s-1", "PFT-level flux of biomass carbon from live to dead pool from mortality", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_carbonflux_si_pft, false),
    (Symbol[:dynam0], "FATES_MORTALITY_FIRE_CFLUX_PF", "kg m-2 s-1", "PFT-level flux of biomass carbon from live to dead pool from fire mortality", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_firemortality_carbonflux_si_pft, false),
    (Symbol[:dynam0], "FATES_MORTALITY_HYDRO_CFLUX_PF", "kg m-2 s-1", "PFT-level flux of biomass carbon from live to dead pool from hydraulic failure mortality", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_hydraulicmortality_carbonflux_si_pft, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CSTARV_CFLUX_PF", "kg m-2 s-1", "PFT-level flux of biomass carbon from live to dead pool from carbon starvation mortality (both continuous and termination)", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_cstarvmortality_carbonflux_si_pft, false),
    (Symbol[:dynam0], "FATES_MORT_CSTARV_CONT_CFLUX_PF", "kg m-2 s-1", "PFT-level flux of biomass carbon from live to dead pool from carbon starvation mortality (Continuous-only, without termination)", "active", "A", site_pft_r8, "CLM:ALM", group_dyna_complx, :ih_cstarvmortality_continuous_carbonflux_si_pft, false),
    (Symbol[:dynam0], "FATES_ABOVEGROUND_MORT_SZPF", "kg m-2 s-1", "Aboveground flux of carbon from AGB to necromass due to mortality", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_abg_mortality_cflux_si_scpf, false),
    (Symbol[:dynam0], "FATES_ABOVEGROUND_PROD_SZPF", "kg m-2 s-1", "Aboveground carbon productivity", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_abg_productivity_cflux_si_scpf, false),
    (Symbol[:dynam0], "FATES_NPLANT_SZAP", "m-2", "number of plants per m2 in each size x age class", "inactive", "A", site_scag_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_si_scag, false),
    (Symbol[:dynam0], "FATES_NPLANT_CANOPY_SZAP", "m-2", "number of plants per m2 in canopy in each size x age class", "inactive", "A", site_scag_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_canopy_si_scag, false),
    (Symbol[:dynam0], "FATES_NPLANT_USTORY_SZAP", "m-2", "number of plants per m2 in understory in each size x age class", "inactive", "A", site_scag_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_understory_si_scag, false),
    (Symbol[:dynam0], "FATES_DDBH_CANOPY_SZAP", "m m-2 yr-1", "growth rate of canopy plants in meters DBH per m2 per year in canopy in each size x age class", "inactive", "A", site_scag_r8, "CLM:ALM", group_dyna_complx, :ih_ddbh_canopy_si_scag, false),
    (Symbol[:dynam0], "FATES_DDBH_USTORY_SZAP", "m m-2 yr-1", "growth rate of understory plants in meters DBH per m2 per year in each size x age class", "inactive", "A", site_scag_r8, "CLM:ALM", group_dyna_complx, :ih_ddbh_understory_si_scag, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CANOPY_SZAP", "m-2 yr-1", "mortality rate of canopy plants in number of plants per m2 per year in each size x age class", "inactive", "A", site_scag_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_canopy_si_scag, false),
    (Symbol[:dynam0], "FATES_MORTALITY_USTORY_SZAP", "m-2 yr-1", "mortality rate of understory plants in number of plants per m2 per year in each size x age class", "inactive", "A", site_scag_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_understory_si_scag, false),
    (Symbol[:dynam0], "FATES_NPLANT_SZAPPF", "m-2", "number of plants per m2 in each size x age x pft class", "inactive", "A", site_scagpft_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_si_scagpft, false),
    (Symbol[:dynam0], "FATES_NPP_APPF", "kg m-2 s-1", "NPP per PFT in each age bin in kg carbon per m2 per second", "inactive", "A", site_agepft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_si_agepft, false),
    (Symbol[:dynam0], "FATES_VEGC_APPF", "kg m-2", "biomass per PFT in each age bin in kg carbon per m2", "inactive", "A", site_agepft_r8, "CLM:ALM", group_dyna_complx, :ih_biomass_si_agepft, false),
    (Symbol[:dynam0], "FATES_SCORCH_HEIGHT_APPF", "m", "SPITFIRE flame Scorch Height (calculated per PFT in each patch age bin)", "inactive", "A", site_agepft_r8, "CLM:ALM", group_dyna_complx, :ih_scorch_height_si_agepft, false),
    (Symbol[:dynam0], "FATES_GPP_SZPF", "kg m-2 s-1", "gross primary production by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_gpp_si_scpf, false),
    (Symbol[:dynam0], "FATES_GPP_CANOPY_SZPF", "kg m-2 s-1", "gross primary production of canopy plants by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_gpp_canopy_si_scpf, false),
    (Symbol[:dynam0], "FATES_AUTORESP_CANOPY_SZPF", "kg m-2 s-1", "autotrophic respiration of canopy plants by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_ar_canopy_si_scpf, false),
    (Symbol[:dynam0], "FATES_GPP_USTORY_SZPF", "kg m-2 s-1", "gross primary production of understory plants by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_gpp_understory_si_scpf, false),
    (Symbol[:dynam0], "FATES_AUTORESP_USTORY_SZPF", "kg m-2 s-1", "autotrophic respiration of understory plants by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_ar_understory_si_scpf, false),
    (Symbol[:dynam0], "FATES_NPP_SZPF", "kg m-2 s-1", "total net primary production by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_totl_si_scpf, false),
    (Symbol[:dynam0], "FATES_LEAF_ALLOC_SZPF", "kg m-2 s-1", "allocation to leaves by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_leaf_si_scpf, false),
    (Symbol[:dynam0], "FATES_SEED_ALLOC_SZPF", "kg m-2 s-1", "allocation to seeds by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_seed_si_scpf, false),
    (Symbol[:dynam0], "FATES_FROOT_ALLOC_SZPF", "kg m-2 s-1", "allocation to fine roots by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_fnrt_si_scpf, false),
    (Symbol[:dynam0], "FATES_BGSAPWOOD_ALLOC_SZPF", "kg m-2 s-1", "allocation to below-ground sapwood by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_bgsw_si_scpf, false),
    (Symbol[:dynam0], "FATES_BGSTRUCT_ALLOC_SZPF", "kg m-2 s-1", "allocation to below-ground structural (deadwood) by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_bgdw_si_scpf, false),
    (Symbol[:dynam0], "FATES_AGSAPWOOD_ALLOC_SZPF", "kg m-2 s-1", "allocation to above-ground sapwood by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_agsw_si_scpf, false),
    (Symbol[:dynam0], "FATES_AGSTRUCT_ALLOC_SZPF", "kg m-2 s-1", "allocation to above-ground structural (deadwood) by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_agdw_si_scpf, false),
    (Symbol[:dynam0], "FATES_STORE_ALLOC_SZPF", "kg m-2 s-1", "allocation to storage C by pft/size in kg carbon per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_npp_stor_si_scpf, false),
    (Symbol[:dynam0], "FATES_DDBH_SZPF", "m m-2 yr-1", "diameter growth increment by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_ddbh_si_scpf, false),
    (Symbol[:dynam0], "FATES_GROWTHFLUX_SZPF", "m-2 yr-1", "flux of individuals into a given size class bin via growth and recruitment", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_growthflux_si_scpf, false),
    (Symbol[:dynam0], "FATES_GROWTHFLUX_FUSION_SZPF", "m-2 yr-1", "flux of individuals into a given size class bin via fusion", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_growthflux_fusion_si_scpf, false),
    (Symbol[:dynam0], "FATES_DDBH_CANOPY_SZPF", "m m-2 yr-1", "diameter growth increment by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_ddbh_canopy_si_scpf, false),
    (Symbol[:dynam0], "FATES_DDBH_USTORY_SZPF", "m m-2 yr-1", "diameter growth increment by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_ddbh_understory_si_scpf, false),
    (Symbol[:dynam0], "FATES_BASALAREA_SZPF", "m2 m-2", "basal area by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_ba_si_scpf, false),
    (Symbol[:dynam0], "FATES_VEGC_ABOVEGROUND_SZPF", "kg m-2", "aboveground biomass by pft/size in kg carbon per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_agb_si_scpf, false),
    (Symbol[:dynam0], "FATES_NPLANT_SZPF", "m-2", "stem number density by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_si_scpf, false),
    (Symbol[:dynam0], "FATES_NPLANT_ACPF", "m-2", "stem number density by pft and age class", "inactive", "A", site_coage_pft_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_si_capf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_BACKGROUND_SZPF", "m-2 yr-1", "background mortality by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m1_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_HYDRAULIC_SZPF", "m-2 yr-1", "hydraulic mortality by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m2_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CSTARV_SZPF", "m-2 yr-1", "carbon starvation mortality by pft/size in number of plants per m2 per year (both continous and termination)", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m3_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_IMPACT_SZPF", "m-2 yr-1", "impact mortality by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m4_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_FIRE_SZPF", "m-2 yr-1", "fire mortality by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m5_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CROWNSCORCH_SZPF", "m-2 yr-1", "fire mortality from crown scorch by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_crownfiremort_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CAMBIALBURN_SZPF", "m-2 yr-1", "fire mortality from cambial burn by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_cambialfiremort_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_TERMINATION_SZPF", "m-2 yr-1", "termination mortality (excluding C-starvation) by pft/size in number pf plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m6_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_LOGGING_SZPF", "m-2 yr-1", "logging mortality by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m7_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_FREEZING_SZPF", "m-2 yr-1", "freezing mortality by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m8_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_SENESCENCE_SZPF", "m-2 yr-1", "senescence mortality by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m9_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_AGESCEN_SZPF", "m-2 yr-1", "age senescence mortality by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m10_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_AGESCEN_ACPF", "m-2 yr-1", "age senescence mortality by pft/cohort age in number of plants per m2 per year", "inactive", "A", site_coage_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m10_si_capf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CANOPY_SZPF", "m-2 yr-1", "total mortality of canopy plants by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_canopy_si_scpf, false),
    (Symbol[:dynam0], "FATES_M3_MORTALITY_CANOPY_SZPF", "m-2 yr-1", "C starvation mortality of canopy plants by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m3_mortality_canopy_si_scpf, false),
    (Symbol[:dynam0], "FATES_M3_MORTALITY_USTORY_SZPF", "m-2 yr-1", "C starvation mortality of understory plants by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m3_mortality_understory_si_scpf, false),
    (Symbol[:dynam0], "FATES_C13DISC_SZPF", "per mil", "C13 discrimination by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_c13disc_si_scpf, false),
    (Symbol[:dynam0], "FATES_STOREC_CANOPY_SZPF", "kg m-2", "biomass in storage pools of canopy plants by pft/size in kg carbon per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_bstor_canopy_si_scpf, false),
    (Symbol[:dynam0], "FATES_LEAFC_CANOPY_SZPF", "kg m-2", "biomass in leaves of canopy plants by pft/size in kg carbon per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_bleaf_canopy_si_scpf, false),
    (Symbol[:dynam0], "FATES_LAI_CANOPY_SZPF", "m2 m-2", "Leaf area index (LAI) of canopy plants by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_lai_canopy_si_scpf, false),
    (Symbol[:dynam0], "FATES_CROWNAREA_CANOPY_SZPF", "m2 m-2", "Total crown area of canopy plants by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_crownarea_canopy_si_scpf, false),
    (Symbol[:dynam0], "FATES_CROWNAREA_USTORY_SZPF", "m2 m-2", "Total crown area of understory plants by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_crownarea_understory_si_scpf, false),
    (Symbol[:dynam0], "FATES_NPLANT_CANOPY_SZPF", "m-2", "number of canopy plants by size/pft per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_canopy_si_scpf, false),
    (Symbol[:dynam0], "FATES_MORTALITY_USTORY_SZPF", "m-2 yr-1", "total mortality of understory plants by pft/size in number of plants per m2 per year", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_understory_si_scpf, false),
    (Symbol[:dynam0], "FATES_STOREC_USTORY_SZPF", "kg m-2", "biomass in storage pools of understory plants by pft/size in kg carbon per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_bstor_understory_si_scpf, false),
    (Symbol[:dynam0], "FATES_LEAFC_USTORY_SZPF", "kg m-2", "biomass in leaves of understory plants by pft/size in kg carbon per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_bleaf_understory_si_scpf, false),
    (Symbol[:dynam0], "FATES_LAI_USTORY_SZPF", "m2 m-2", "Leaf area index (LAI) of understory plants by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_lai_understory_si_scpf, false),
    (Symbol[:dynam0], "FATES_NPLANT_USTORY_SZPF", "m-2", "density of understory plants by pft/size in number of plants per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_understory_si_scpf, false),
    (Symbol[:dynam0], "FATES_CWD_ABOVEGROUND_DC", "kg m-2", "debris class-level aboveground coarse woody debris stocks in kg carbon per m2", "inactive", "A", site_cwdsc_r8, "CLM:ALM", group_dyna_complx, :ih_cwd_ag_si_cwdsc, false),
    (Symbol[:dynam0], "FATES_CWD_BELOWGROUND_DC", "kg m-2", "debris class-level belowground coarse woody debris stocks in kg carbon per m2", "inactive", "A", site_cwdsc_r8, "CLM:ALM", group_dyna_complx, :ih_cwd_bg_si_cwdsc, false),
    (Symbol[:dynam0], "FATES_CWD_ABOVEGROUND_IN_DC", "kg m-2 s-1", "debris class-level aboveground coarse woody debris input in kg carbon per m2 per second", "inactive", "A", site_cwdsc_r8, "CLM:ALM", group_dyna_complx, :ih_cwd_ag_in_si_cwdsc, false),
    (Symbol[:dynam0], "FATES_CWD_BELOWGROUND_IN_DC", "kg m-2 s-1", "debris class-level belowground coarse woody debris input in kg carbon per m2 per second", "inactive", "A", site_cwdsc_r8, "CLM:ALM", group_dyna_complx, :ih_cwd_bg_in_si_cwdsc, false),
    (Symbol[:dynam0], "FATES_CWD_ABOVEGROUND_OUT_DC", "kg m-2 s-1", "debris class-level aboveground coarse woody debris output in kg carbon per m2 per second", "inactive", "A", site_cwdsc_r8, "CLM:ALM", group_dyna_complx, :ih_cwd_ag_out_si_cwdsc, false),
    (Symbol[:dynam0], "FATES_CWD_BELOWGROUND_OUT_DC", "kg m-2 s-1", "debris class-level belowground coarse woody debris output in kg carbon per m2 per second", "inactive", "A", site_cwdsc_r8, "CLM:ALM", group_dyna_complx, :ih_cwd_bg_out_si_cwdsc, false),
    (Symbol[:dynam0], "FATES_DDBH_CANOPY_SZ", "m m-2 yr-1", "diameter growth increment by size of canopy plants", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_ddbh_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_DDBH_USTORY_SZ", "m m-2 yr-1", "diameter growth increment by size of understory plants", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_ddbh_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_YESTCANLEV_CANOPY_SZ", "m-2", "yesterdays canopy level for canopy plants by size class in number of plants per m2", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_yesterdaycanopylevel_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_YESTCANLEV_USTORY_SZ", "m-2", "yesterdays canopy level for understory plants by size class in number of plants per m2", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_yesterdaycanopylevel_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_BASALAREA_SZ", "m2 m-2", "basal area by size class", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_ba_si_scls, false),
    (Symbol[:dynam0], "FATES_VEGC_ABOVEGROUND_SZ", "kg m-2", "aboveground biomass by size class in kg carbon per m2", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_agb_si_scls, false),
    (Symbol[:dynam0], "FATES_VEGC_SZ", "kg m-2", "total biomass by size class in kg carbon per m2", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_biomass_si_scls, false),
    (Symbol[:dynam0], "FATES_DEMOTION_RATE_SZ", "m-2 yr-1", "demotion rate from canopy to understory by size class in number of plants per m2 per year", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_demotion_rate_si_scls, false),
    (Symbol[:dynam0], "FATES_PROMOTION_RATE_SZ", "m-2 yr-1", "promotion rate from understory to canopy by size class", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_promotion_rate_si_scls, false),
    (Symbol[:dynam0], "FATES_NPLANT_CANOPY_SZ", "m-2", "number of canopy plants per m2 by size class", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_LAI_CANOPY_SZ", "m2 m-2", "leaf area index (LAI) of canopy plants by size class", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_lai_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_SAI_CANOPY_SZ", "m2 m-2", "stem area index (SAI) of canopy plants by size class", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_sai_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CANOPY_SZ", "m-2 yr-1", "total mortality of canopy trees by size class in number of plants per m2", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CANOPY_SE_SZ", "m-2 yr-1", "'total mortality of canopy trees by size class in number of plants per m2", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_canopy_secondary_si_scls, false),
    (Symbol[:dynam0], "FATES_NPLANT_USTORY_SZ", "m-2", "number of understory plants per m2 by size class", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_M3_MORTALITY_CANOPY_SZ", "m-2 yr-1", "C starvation mortality of canopy plants by size", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m3_mortality_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_M3_MORTALITY_USTORY_SZ", "m-2 yr-1", "C starvation mortality of understory plants by size", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m3_mortality_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_LAI_USTORY_SZ", "m2 m-2", "leaf area index (LAI) of understory plants by size class", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_lai_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_SAI_USTORY_SZ", "m2 m-2", "stem area index (SAI) of understory plants by size class", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_sai_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_NPLANT_SZ", "m-2", "number of plants per m2 by size class", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_si_scls, false),
    (Symbol[:dynam0], "FATES_NPLANT_AC", "m-2", "number of plants per m2 by cohort age class", "active", "A", site_coage_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_si_cacls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_BACKGROUND_SZ", "m-2 yr-1", "background mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m1_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_HYDRAULIC_SZ", "m-2 yr-1", "hydraulic mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m2_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CSTARV_SZ", "m-2 yr-1", "carbon starvation mortality by size in number of plants per m2 per year (both continous and termination)", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m3_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_BACKGROUND_SE_SZ", "m-2 yr-1", "'background mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m1_sec_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_HYDRAULIC_SE_SZ", "m-2 yr-1", "'hydraulic mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m2_sec_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_CSTARV_SE_SZ", "m-2 yr-1", "'carbon starvation mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m3_sec_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_IMPACT_SZ", "m-2 yr-1", "impact mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m4_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_FIRE_SZ", "m-2 yr-1", "fire mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m5_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_TERMINATION_SZ", "m-2 yr-1", "termination mortality (excluding C-starvation) by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m6_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_LOGGING_SZ", "m-2 yr-1", "logging mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m7_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_FREEZING_SZ", "m-2 yr-1", "freezing mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m8_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_SENESCENCE_SZ", "m-2 yr-1", "senescence mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m9_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_AGESCEN_SZ", "m-2 yr-1", "age senescence mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m10_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_AGESCEN_AC", "m-2 yr-1", "age senescence mortality by cohort age in number of plants per m2 per year", "active", "A", site_coage_r8, "CLM:ALM", group_dyna_complx, :ih_m10_si_cacls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_LOGGING_SE_SZ", "m-2 yr-1", "'logging mortality by size in number of plants per m2 per event", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m7_sec_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_FREEZING_SE_SZ", "m-2 event-1", "'freezing mortality by size in number of plants per m2 per event", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m8_sec_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_SENESCENCE_SE_SZ", "m-2 yr-1", "'senescence mortality by size in number of plants per m2 per event", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m9_sec_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_AGESCEN_SE_SZ", "m-2 yr-1", "'age senescence mortality by size in number of plants per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_m10_sec_si_scls, false),
    (Symbol[:dynam0], "FATES_NPP_CANOPY_SZ", "kg m-2 s-1", "NPP of canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_carbon_balance_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_NPP_USTORY_SZ", "kg m-2 s-1", "NPP of understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_carbon_balance_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_MORTALITY_USTORY_SZ", "m-2 yr-1", "total mortality of understory trees by size class in individuals per m2 per year", "active", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_TRIMMING_CANOPY_SZ", "m-2", "'trimming term of canopy plants weighted by plant density", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_trimming_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_TRIMMING_USTORY_SZ", "m-2", "'trimming term of understory plants weighted by plant density", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_trimming_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_CROWNAREA_CANOPY_SZ", "m2 m-2", "total crown area of canopy plants by size class", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_crown_area_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_CROWNAREA_USTORY_SZ", "m2 m-2", "total crown area of understory plants by size class", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_crown_area_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_LEAFCTURN_CANOPY_SZ", "kg m-2 s-1", "leaf turnover (non-mortal) for canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_leaf_md_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_FROOTCTURN_CANOPY_SZ", "kg m-2 s-1", "fine root turnover (non-mortal) for canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_root_md_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_STORECTURN_CANOPY_SZ", "kg m-2 s-1", "storage turnover (non-mortal) for canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_bstore_md_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_STRUCTCTURN_CANOPY_SZ", "kg m-2 s-1", "structural C turnover (non-mortal) for canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_bdead_md_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_SAPWOODCTURN_CANOPY_SZ", "kg m-2 s-1", "sapwood turnover (non-mortal) for canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_bsw_md_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_SEED_PROD_CANOPY_SZ", "kg m-2 s-1", "seed production of canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_seed_prod_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_LEAF_ALLOC_CANOPY_SZ", "kg m-2 s-1", "allocation to leaves for canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_leaf_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_FROOT_ALLOC_CANOPY_SZ", "kg m-2 s-1", "allocation to fine root C for canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_fnrt_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_SAPWOOD_ALLOC_CANOPY_SZ", "kg m-2 s-1", "allocation to sapwood C for canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_sapw_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_STRUCT_ALLOC_CANOPY_SZ", "kg m-2 s-1", "allocation to structural C for canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_dead_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_SEED_ALLOC_CANOPY_SZ", "kg m-2 s-1", "allocation to reproductive C for canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_seed_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_STORE_ALLOC_CANOPY_SZ", "kg m-2 s-1", "allocation to storage C for canopy plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_stor_canopy_si_scls, false),
    (Symbol[:dynam0], "FATES_LEAFCTURN_USTORY_SZ", "kg m-2 s-1", "leaf turnover (non-mortal) for understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_leaf_md_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_FROOTCTURN_USTORY_SZ", "kg m-2 s-1", "fine root turnover (non-mortal) for understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_root_md_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_STORECTURN_USTORY_SZ", "kg m-2 s-1", "storage C turnover (non-mortal) for understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_bstore_md_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_STRUCTCTURN_USTORY_SZ", "kg m-2 s-1", "structural C turnover (non-mortal) for understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_bdead_md_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_SAPWOODCTURN_USTORY_SZ", "kg m-2 s-1", "sapwood C turnover (non-mortal) for understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_bsw_md_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_SEED_PROD_USTORY_SZ", "kg m-2 s-1", "seed production of understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_seed_prod_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_LEAF_ALLOC_USTORY_SZ", "kg m-2 s-1", "allocation to leaves for understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_leaf_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_FROOT_ALLOC_USTORY_SZ", "kg m-2 s-1", "allocation to fine roots for understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_fnrt_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_SAPWOOD_ALLOC_USTORY_SZ", "kg m-2 s-1", "allocation to sapwood C for understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_sapw_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_STRUCT_ALLOC_USTORY_SZ", "kg m-2 s-1", "allocation to structural C for understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_dead_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_SEED_ALLOC_USTORY_SZ", "kg m-2 s-1", "allocation to reproductive C for understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_seed_understory_si_scls, false),
    (Symbol[:dynam0], "FATES_STORE_ALLOC_USTORY_SZ", "kg m-2 s-1", "allocation to storage C for understory plants by size class in kg carbon per m2 per second", "inactive", "A", site_size_r8, "CLM:ALM", group_dyna_complx, :ih_npp_stor_understory_si_scls, false),
    (Symbol[:dynam0, :treedamage], "FATES_NPLANT_CDPF", "m-2", "N. plants per damage x size x pft class", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_NPLANT_CANOPY_CDPF", "m-2", "N. plants per damage x size x pft class", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_canopy_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_NPLANT_USTORY_CDPF", "m-2", "N. plants in the understory per damage x size x pft class", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_nplant_understory_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_M3_CDPF", "m-2 yr-1", "carbon starvation mortality by damaage/pft/size", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_m3_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_M11_SZPF", "m-2 yr-1", "damage mortality by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_m11_si_scpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_M11_CDPF", "m-2 yr-1", "damage mortality by damaage/pft/size", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_m11_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_MORTALITY_CDPF", "m-2 yr-1", "mortality by damage class by size by pft", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_M3_MORTALITY_CANOPY_CDPF", "m-2 yr-1", "C starvation mortality of canopy plants by damage/pft/size", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_m3_mortality_canopy_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_M3_MORTALITY_USTORY_CDPF", "m-2 yr-1", "C starvation mortality of understory plants by pft/size", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_m3_mortality_understory_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_M11_MORTALITY_CANOPY_CDPF", "m-2 yr-1", "damage mortality of canopy plants by damage/pft/size", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_m11_mortality_canopy_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_M11_MORTALITY_USTORY_CDPF", "m-2 yr-1", "damage mortality of understory plants by pft/size", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_m11_mortality_understory_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_MORTALITY_CANOPY_CDPF", "m-2 yr-1", "mortality of canopy plants by damage/pft/size", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_canopy_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_MORTALITY_USTORY_CDPF", "m-2 yr-1", "mortality of understory plants by pft/size", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_mortality_understory_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_DDBH_CDPF", "m m-2 yr-1", "ddbh annual increment growth by damage x size pft", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_ddbh_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_DDBH_CANOPY_CDPF", "m m-2 yr-1", "ddbh annual canopy increment growth by damage x size pft", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_ddbh_canopy_si_cdpf, false),
    (Symbol[:dynam0, :treedamage], "FATES_DDBH_USTORY_CDPF", "m m-2 yr-1", "ddbh annual understory increment growth by damage x size pft", "inactive", "A", site_cdpf_r8, "CLM:ALM", group_dyna_complx, :ih_ddbh_understory_si_cdpf, false),
    (Symbol[:dynam0], "FATES_FIRE_FLUX_EL", "kg m-2 s-1", "loss to atmosphere from fire by element in kg element per m2 per s", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_burn_flux_elem, false),
    (Symbol[:dynam0, :not_st3_not_sp], "FATES_ERROR_EL", "kg s-1", "total mass-balance error in kg per second by element", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_err_fates_elem, false),
    (Symbol[:dynam0, :not_st3_not_sp], "FATES_INTERR_LIVEVEG_EL", "kg m-2", "Bias error between integrated flux and (minus) state in live vegetation ", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_interr_liveveg_elem, false),
    (Symbol[:dynam0, :not_st3_not_sp], "FATES_INTERR_LITTER_EL", "kg m-2", "Bias error between integrated flux and (minus) state in litter ", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_interr_litter_elem, false),
    (Symbol[:dynam0], "FATES_LITTER_AG_FINE_EL", "kg m-2", "mass of aboveground litter in fines (leaves, nonviable seed) by element", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_fines_ag_elem, false),
    (Symbol[:dynam0], "FATES_LITTER_BG_FINE_EL", "kg m-2", "mass of belowground litter in fines (fineroots) by element", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_fines_bg_elem, false),
    (Symbol[:dynam0], "FATES_LITTER_BG_CWD_EL", "kg m-2", "mass of belowground litter in coarse woody debris (coarse roots) by element", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_cwd_bg_elem, false),
    (Symbol[:dynam0], "FATES_LITTER_AG_CWD_EL", "kg m-2", "mass of aboveground litter in coarse woody debris (trunks/branches/twigs) by element", "active", "A", site_elem_r8, "CLM:ALM", group_dyna_complx, :ih_cwd_ag_elem, false),
    (Symbol[:dynam0], "FATES_LITTER_CWD_ELDC", "kg m-2", "total mass of litter in coarse woody debris by element and coarse woody debris size", "active", "A", site_elcwd_r8, "CLM:ALM", group_dyna_complx, :ih_cwd_elcwd, false),
    (Symbol[:dynam0], "FATES_VEGC_SZPF", "kg m-2", "total vegetation biomass in live plants by size-class x pft in kg carbon per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_totvegc_scpf, false),
    (Symbol[:dynam0], "FATES_LEAFC_SZPF", "kg m-2", "leaf carbon mass by size-class x pft in kg carbon per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_leafc_scpf, false),
    (Symbol[:dynam0], "FATES_FROOTC_SZPF", "kg m-2", "fine-root carbon mass by size-class x pft in kg carbon per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_fnrtc_scpf, false),
    (Symbol[:dynam0], "FATES_SAPWOODC_SZPF", "kg m-2", "sapwood carbon mass by size-class x pft in kg carbon per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_sapwc_scpf, false),
    (Symbol[:dynam0], "FATES_STOREC_SZPF", "kg m-2", "storage carbon mass by size-class x pft in kg carbon per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_storec_scpf, false),
    (Symbol[:dynam0], "FATES_REPROC_SZPF", "kg m-2", "reproductive carbon mass (on plant) by size-class x pft in kg carbon per m2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_dyna_complx, :ih_reproc_scpf, false),
    (Symbol[:hifrq0], "FATES_STOMATAL_COND", "mol m-2 s-1", "mean stomatal conductance", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_c_stomata_si, false),
    (Symbol[:hifrq0], "FATES_LBLAYER_COND", "mol m-2 s-1", "mean leaf boundary layer conductance", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_c_lblayer_si, false),
    (Symbol[:hifrq0], "FATES_TVEG", "degree_Celsius", "fates instantaneous mean vegetation temperature by site", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_tveg_si, false),
    (Symbol[:hifrq0], "FATES_VIS_RAD_ERROR", "-", "mean two-stream solver error for VIS", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_vis_rad_err_si, false),
    (Symbol[:hifrq0], "FATES_NIR_RAD_ERROR", "-", "mean two-stream solver error for NIR", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_nir_rad_err_si, false),
    (Symbol[:hifrq0], "FATES_NPP", "kg m-2 s-1", "net primary production in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_npp_si, false),
    (Symbol[:hifrq0], "FATES_NPP_SECONDARY", "kg m-2 s-1", "'net primary production in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_npp_secondary_si, false),
    (Symbol[:hifrq0], "FATES_GPP", "kg m-2 s-1", "gross primary production in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_gpp_si, false),
    (Symbol[:hifrq0], "FATES_GPP_SECONDARY", "kg m-2 s-1", "'gross primary production in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_gpp_secondary_si, false),
    (Symbol[:hifrq0], "FATES_AUTORESP", "kg m-2 s-1", "autotrophic respiration in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_aresp_si, false),
    (Symbol[:hifrq0], "FATES_AUTORESP_SECONDARY", "kg m-2 s-1", "'autotrophic respiration in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_aresp_secondary_si, false),
    (Symbol[:hifrq0], "FATES_GROWTH_RESP", "kg m-2 s-1", "growth respiration in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_growth_resp_si, false),
    (Symbol[:hifrq0], "FATES_GROWTH_RESP_SECONDARY", "kg m-2 s-1", "'growth respiration in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_growth_resp_secondary_si, false),
    (Symbol[:hifrq0], "FATES_MAINT_RESP", "kg m-2 s-1", "maintenance respiration in kg carbon per m2 land area per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_maint_resp_si, false),
    (Symbol[:hifrq0], "FATES_MAINT_RESP_UNREDUCED", "kg m-2 s-1", "diagnostic maintenance respiration if the low-carbon-storage reduction is ignored", "unactive", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_maint_resp_unreduced_si, false),
    (Symbol[:hifrq0], "FATES_MAINT_RESP_SECONDARY", "kg m-2 s-1", "'maintenance respiration in kg carbon per m2 land area per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_maint_resp_secondary_si, false),
    (Symbol[:hifrq0], "FATES_GPP_CANOPY", "kg m-2 s-1", "gross primary production of canopy plants in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_gpp_canopy_si, false),
    (Symbol[:hifrq0], "FATES_AUTORESP_CANOPY", "kg m-2 s-1", "autotrophic respiration of canopy plants in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_ar_canopy_si, false),
    (Symbol[:hifrq0], "FATES_GPP_USTORY", "kg m-2 s-1", "gross primary production of understory plants in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_gpp_understory_si, false),
    (Symbol[:hifrq0], "FATES_AUTORESP_USTORY", "kg m-2 s-1", "autotrophic respiration of understory plants in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_ar_understory_si, false),
    (Symbol[:hifrq0], "FATES_LEAFMAINTAR", "kg m-2 s-1", "leaf maintenance autotrophic respiration in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_leaf_mr_si, false),
    (Symbol[:hifrq0], "FATES_FROOTMAINTAR", "kg m-2 s-1", "fine root maintenance autotrophic respiration in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_froot_mr_si, false),
    (Symbol[:hifrq0], "FATES_CROOTMAINTAR", "kg m-2 s-1", "live coarse root maintenance autotrophic respiration in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_livecroot_mr_si, false),
    (Symbol[:hifrq0], "FATES_LSTEMMAINTAR", "kg m-2 s-1", "live stem maintenance autotrophic respiration in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_livestem_mr_si, false),
    (Symbol[:hifrq0], "FATES_NEP", "kg m-2 s-1", "net ecosystem production in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_nep_si, false),
    (Symbol[:hifrq0], "FATES_HET_RESP", "kg m-2 s-1", "heterotrophic respiration in kg carbon per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hifr_simple, :ih_hr_si, false),
    (Symbol[:hifrq0, :planthydro], "FATES_SAPFLOW", "kg m-2 s-1", "areal sap flow rate in kg per m2 per second", "active", "A", site_r8, "CLM:ALM", group_hydr_simple, :ih_sapflow_si, false),
    (Symbol[:hifrq0, :planthydro], "FATES_ROOTWGT_SOILVWC", "m3 m-3", "'soil volumetric water content", "active", "A", site_r8, "CLM:ALM", group_hydr_simple, :ih_rootwgt_soilvwc_si, false),
    (Symbol[:hifrq0, :planthydro], "FATES_ROOTWGT_SOILVWCSAT", "m3 m-3", "'soil saturated volumetric water content", "active", "A", site_r8, "CLM:ALM", group_hydr_simple, :ih_rootwgt_soilvwcsat_si, false),
    (Symbol[:hifrq0, :planthydro], "FATES_ROOTWGT_SOILMATPOT", "Pa", "'soil matric potential", "active", "A", site_r8, "CLM:ALM", group_hydr_simple, :ih_rootwgt_soilmatpot_si, false),
    (Symbol[:hifrq0, :planthydro], "FATES_ROOTUPTAKE", "kg m-2 s-1", "root water uptake rate", "active", "A", site_r8, "CLM:ALM", group_hydr_simple, :ih_rootuptake_si, false),
    (Symbol[:hifrq0, :planthydro], "FATES_VEGH2O", "kg m-2", "water stored inside vegetation tissues (leaf, stem, roots)", "inactive", "A", site_r8, "CLM:ALM", group_hydr_simple, :ih_h2oveg_si, false),
    (Symbol[:hifrq0, :planthydro], "FATES_VEGH2O_HYDRO_ERR", "kg m-2", "cumulative net borrowed (+) from plant_stored_h2o due to plant hydrodynamics", "inactive", "A", site_r8, "CLM:ALM", group_hydr_simple, :ih_h2oveg_hydro_err_si, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_NPP_AP", "kg m-2 s-1", "net primary productivity by age bin in kg carbon per m2 per second", "inactive", "A", site_age_r8, "CLM:ALM", group_hifr_complx, :ih_npp_si_age, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_GPP_AP", "kg m-2 s-1", "gross primary productivity by age bin in kg carbon per m2 per second", "inactive", "A", site_age_r8, "CLM:ALM", group_hifr_complx, :ih_gpp_si_age, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_RDARK_USTORY_SZ", "kg m-2 s-1", "dark respiration for understory plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_rdark_understory_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_LSTEMMAINTAR_USTORY_SZ", "kg m-2 s-1", "live stem maintenance autotrophic respiration for understory plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_livestem_mr_understory_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_CROOTMAINTAR_USTORY_SZ", "kg m-2 s-1", "live coarse root maintenance autotrophic respiration for understory plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_livecroot_mr_understory_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_FROOTMAINTAR_USTORY_SZ", "kg m-2 s-1", "fine root maintenance autotrophic respiration for understory plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_froot_mr_understory_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_GROWAR_USTORY_SZ", "kg m-2 s-1", "growth autotrophic respiration of understory plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_resp_g_understory_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_MAINTAR_USTORY_SZ", "kg m-2 s-1", "maintenance autotrophic respiration of understory plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_resp_m_understory_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_RDARK_CANOPY_SZ", "kg m-2 s-1", "dark respiration for canopy plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_rdark_canopy_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_CROOTMAINTAR_CANOPY_SZ", "kg m-2 s-1", "live coarse root maintenance autotrophic respiration for canopy plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_livecroot_mr_canopy_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_FROOTMAINTAR_CANOPY_SZ", "kg m-2 s-1", "live coarse root maintenance autotrophic respiration for canopy plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_froot_mr_canopy_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_GROWAR_CANOPY_SZ", "kg m-2 s-1", "growth autotrophic respiration of canopy plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_resp_g_canopy_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_MAINTAR_CANOPY_SZ", "kg m-2 s-1", "maintenance autotrophic respiration of canopy plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_resp_m_canopy_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_LSTEMMAINTAR_CANOPY_SZ", "kg m-2 s-1", "live stem maintenance autotrophic respiration for canopy plants in kg carbon per m2 per second by size", "inactive", "A", site_size_r8, "CLM:ALM", group_hifr_complx, :ih_livestem_mr_canopy_si_scls, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_AUTORESP_SZPF", "kg m-2 s-1", "total autotrophic respiration in kg carbon per m2 per second by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hifr_complx, :ih_ar_si_scpf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_GROWAR_SZPF", "kg m-2 s-1", "growth autotrophic respiration in kg carbon per m2 per second by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hifr_complx, :ih_ar_grow_si_scpf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_MAINTAR_SZPF", "kg m-2 s-1", "maintenance autotrophic respiration in kg carbon per m2 per second by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hifr_complx, :ih_ar_maint_si_scpf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_RDARK_SZPF", "kg m-2 s-1", "dark portion of maintenance autotrophic respiration in kg carbon per m2 per second by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hifr_complx, :ih_ar_darkm_si_scpf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_AGSAPMAINTAR_SZPF", "kg m-2 s-1", "above-ground sapwood maintenance autotrophic respiration in kg carbon per m2 per second by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hifr_complx, :ih_ar_agsapm_si_scpf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_BGSAPMAINTAR_SZPF", "kg m-2 s-1", "below-ground sapwood maintenance autotrophic respiration in kg carbon per m2 per second by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hifr_complx, :ih_ar_crootm_si_scpf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_FROOTMAINTAR_SZPF", "kg m-2 s-1", "fine root maintenance autotrophic respiration in kg carbon per m2 per second by pft/size", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hifr_complx, :ih_ar_frootm_si_scpf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_PARSUN_CLLL", "W m-2", "PAR absorbed in the sun by each canopy and leaf layer", "inactive", "A", site_cnlf_r8, "CLM:ALM", group_hifr_complx, :ih_parsun_z_si_cnlf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_PARSHA_CLLL", "W m-2", "PAR absorbed in the shade by each canopy and leaf layer", "inactive", "A", site_cnlf_r8, "CLM:ALM", group_hifr_complx, :ih_parsha_z_si_cnlf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_PARSUN_CLLLPF", "W m-2", "'PAR absorbed in the sun by each canopy", "inactive", "A", site_cnlfpft_r8, "CLM:ALM", group_hifr_complx, :ih_parsun_z_si_cnlfpft, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_PARSHA_CLLLPF", "W m-2", "'PAR absorbed in the shade by each canopy", "inactive", "A", site_cnlfpft_r8, "CLM:ALM", group_hifr_complx, :ih_parsha_z_si_cnlfpft, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_PARSUN_CL", "W m-2", "PAR absorbed by sunlit leaves in each canopy layer", "inactive", "A", site_can_r8, "CLM:ALM", group_hifr_complx, :ih_parsun_si_can, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_PARSHA_CL", "W m-2", "PAR absorbed by shaded leaves in each canopy layer", "inactive", "A", site_can_r8, "CLM:ALM", group_hifr_complx, :ih_parsha_si_can, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_LAISUN_CLLL", "m2 m-2", "LAI in the sun by each canopy and leaf layer", "inactive", "A", site_cnlf_r8, "CLM:ALM", group_hifr_complx, :ih_laisun_z_si_cnlf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_LAISHA_CLLL", "m2 m-2", "LAI in the shade by each canopy and leaf layer", "inactive", "A", site_cnlf_r8, "CLM:ALM", group_hifr_complx, :ih_laisha_z_si_cnlf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_LAISUN_CLLLPF", "m2 m-2", "'Sunlit leaf area by each canopy", "inactive", "A", site_cnlfpft_r8, "CLM:ALM", group_hifr_complx, :ih_laisun_clllpf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_LAISHA_CLLLPF", "m2 m-2", "'Shaded leaf area by each canopy", "inactive", "A", site_cnlfpft_r8, "CLM:ALM", group_hifr_complx, :ih_laisha_clllpf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_PARPROF_DIR_CLLLPF", "W m-2", "'radiative profile of direct PAR through each canopy", "inactive", "A", site_cnlfpft_r8, "CLM:ALM", group_hifr_complx, :ih_parprof_dir_si_cnlfpft, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_PARPROF_DIF_CLLLPF", "W m-2", "'radiative profile of diffuse PAR through each canopy", "inactive", "A", site_cnlfpft_r8, "CLM:ALM", group_hifr_complx, :ih_parprof_dif_si_cnlfpft, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_LAISUN_CL", "m2 m-2", "LAI of sunlit leaves by canopy layer", "inactive", "A", site_can_r8, "CLM:ALM", group_hifr_complx, :ih_laisun_si_can, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_LAISHA_CL", "m2 m-2", "LAI of shaded leaves by canopy layer", "inactive", "A", site_can_r8, "CLM:ALM", group_hifr_complx, :ih_laisha_si_can, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_PARPROF_DIR_CLLL", "W m-2", "radiative profile of direct PAR through each canopy and leaf layer (averaged across PFTs)", "inactive", "A", site_cnlf_r8, "CLM:ALM", group_hifr_complx, :ih_parprof_dir_si_cnlf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_PARPROF_DIF_CLLL", "W m-2", "radiative profile of diffuse PAR through each canopy and leaf layer (averaged across PFTs)", "inactive", "A", site_cnlf_r8, "CLM:ALM", group_hifr_complx, :ih_parprof_dif_si_cnlf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_NET_C_UPTAKE_CLLL", "kg m-2 s-1", "net carbon uptake in kg carbon per m2 per second by each canopy and leaf layer per unit ground area (i.e. divide by CROWNAREA_CLLL to make per leaf area)", "inactive", "A", site_cnlf_r8, "CLM:ALM", group_hifr_complx, :ih_ts_net_uptake_si_cnlf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_CROWNFRAC_CLLLPF", "m2 m-2", "area fraction of the canopy footprint occupied by each canopy-leaf-pft layer", "inactive", "A", site_cnlfpft_r8, "CLM:ALM", group_hifr_complx, :ih_crownfrac_clllpf, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_LBLAYER_COND_AP", "mol m-2 s-1", "mean leaf boundary layer conductance - by patch age", "inactive", "A", site_age_r8, "CLM:ALM", group_hifr_complx, :ih_c_lblayer_si_age, false),
    (Symbol[:hifrq0, :hifrq1], "FATES_STOMATAL_COND_AP", "mol m-2 s-1", "mean stomatal conductance - by patch age", "inactive", "A", site_age_r8, "CLM:ALM", group_hifr_complx, :ih_c_stomata_si_age, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_ERRH2O_SZPF", "kg s-1", "mean individual water balance error in kg per individual per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_errh2o_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_TRAN_SZPF", "kg s-1", "mean individual transpiration rate in kg per individual per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_tran_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_SAPFLOW_SZPF", "kg m-2 s-1", "areal sap flow rate dimensioned by size x pft in kg per m2 per second", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_sapflow_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_ITERH1_SZPF", "count indiv-1 step-1", "water balance error iteration diagnostic 1", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_iterh1_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_ITERH2_SZPF", "count indiv-1 step-1", "water balance error iteration diagnostic 2", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_iterh2_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_ABSROOT_H2O_SZPF", "m3 m-3", "absorbing volumetric root water content by size class x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_ath_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_TRANSROOT_H2O_SZPF", "m3 m-3", "transporting volumetric root water content by size class x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_tth_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_STEM_H2O_SZPF", "m3 m-3", "stem volumetric water content by size class x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_sth_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_LEAF_H2O_SZPF", "m3 m-3", "leaf volumetric water content by size class x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_lth_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_ABSROOT_H2OPOT_SZPF", "Pa", "absorbing root water potential by size class x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_awp_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_TRANSROOT_H2OPOT_SZPF", "Pa", "transporting root water potential by size class x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_twp_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_STEM_H2OPOT_SZPF", "Pa", "stem water potential by size class x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_swp_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_LEAF_H2OPOT_SZPF", "Pa", "leaf water potential by size class x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_lwp_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_ABSROOT_CONDFRAC_SZPF", "1", "absorbing root fraction (0-1) of condutivity by size class x pft", "active", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_aflc_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_TRANSROOT_CONDFRAC_SZPF", "1", "transporting root fraction (0-1) of condutivity by size class x pft", "active", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_tflc_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_STEM_CONDFRAC_SZPF", "1", "stem water fraction (0-1) of condutivity by size class x pft", "active", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_sflc_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_LEAF_CONDFRAC_SZPF", "1", "leaf water fraction (0-1) of condutivity by size class x pft", "active", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_lflc_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_BTRAN_SZPF", "1", "mean individual level BTRAN by size class x pft", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_btran_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_SOILMATPOT_SL", "Pa", "soil water matric potenial by soil layer", "inactive", "A", site_soil_r8, "CLM:ALM", group_hydr_complx, :ih_soilmatpot_sl, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_SOILVWC_SL", "m3 m-3", "soil volumetric water content by soil layer", "inactive", "A", site_soil_r8, "CLM:ALM", group_hydr_complx, :ih_soilvwc_sl, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_SOILVWCSAT_SL", "m3 m-3", "soil saturated volumetric water content by soil layer", "inactive", "A", site_soil_r8, "CLM:ALM", group_hydr_complx, :ih_soilvwcsat_sl, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_ROOTUPTAKE_SL", "kg m-2 s-1", "root water uptake rate by soil layer", "inactive", "A", site_soil_r8, "CLM:ALM", group_hydr_complx, :ih_rootuptake_sl, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_ROOTUPTAKE0_SZPF", "kg m-2 m-1 s-1", "'root water uptake from 0 to to 10 cm depth", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_rootuptake0_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_ROOTUPTAKE10_SZPF", "kg m-2 m-1 s-1", "'root water uptake from 10 to to 50 cm depth", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_rootuptake10_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_ROOTUPTAKE50_SZPF", "kg m-2 m-1 s-1", "'root water uptake from 50 to to 100 cm depth", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_rootuptake50_scpf, false),
    (Symbol[:hifrq0, :hifrq1, :planthydro], "FATES_ROOTUPTAKE100_SZPF", "kg m-2 m-1 s-1", "'root water uptake below 100 cm depth", "inactive", "A", site_size_pft_r8, "CLM:ALM", group_hydr_complx, :ih_rootuptake100_scpf, false),
]


# ---------------------------------------------------------------------------
# fates_history_interface_type
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct fates_history_interface_type
    hvars::Vector{fates_history_variable_type} = fates_history_variable_type[]
    num_history_vars_::Int = 0
    dim_kinds::Vector{fates_io_variable_kind_type} =
        [fates_io_variable_kind_type() for _ in 1:fates_history_num_dim_kinds]
    dim_bounds::Vector{fates_io_dimension_type} =
        [fates_io_dimension_type() for _ in 1:fates_history_num_dimensions]
    # Per-dimension index handles (1-based slot into dim_bounds/dim_kinds).
    # Keyed by the dimension symbol used in _DIM_SETUP / _ASSEMBLE_MAP.
    dim_index::Dict{Symbol,Int} = Dict{Symbol,Int}()
    # Per-variable history-output index handles (ih_* -> hvars slot, 0 if unused).
    ih::Dict{Symbol,Int} = Dict{Symbol,Int}(h => 0 for h in _IH_HANDLE_NAMES)
end

# =====================================================================================
# Dimension index get/set accessors (Fortran column_index()/set_column_index(), ...).
# Generic Dict-backed forms; the Fortran had one named pair per dimension.
# =====================================================================================

set_dim_index!(this::fates_history_interface_type, name::Symbol, index::Integer) =
    (this.dim_index[name] = Int(index); nothing)

dim_index_of(this::fates_history_interface_type, name::Symbol) = this.dim_index[name]

# Named accessors mirroring the Fortran type-bound getters.
column_index(this::fates_history_interface_type) = this.dim_index[:column]
levsoil_index(this::fates_history_interface_type) = this.dim_index[:levsoil]
levscpf_index(this::fates_history_interface_type) = this.dim_index[:levscpf]
levscls_index(this::fates_history_interface_type) = this.dim_index[:levscls]
levcacls_index(this::fates_history_interface_type) = this.dim_index[:levcacls]
levcapf_index(this::fates_history_interface_type) = this.dim_index[:levcapf]
levpft_index(this::fates_history_interface_type) = this.dim_index[:levpft]
levage_index(this::fates_history_interface_type) = this.dim_index[:levage]
levfuel_index(this::fates_history_interface_type) = this.dim_index[:levfuel]
levcwdsc_index(this::fates_history_interface_type) = this.dim_index[:levcwdsc]
levcan_index(this::fates_history_interface_type) = this.dim_index[:levcan]
levcnlf_index(this::fates_history_interface_type) = this.dim_index[:levcnlf]
levcnlfpft_index(this::fates_history_interface_type) = this.dim_index[:levcnlfpft]
levcdpf_index(this::fates_history_interface_type) = this.dim_index[:levcdpf]
levcdsc_index(this::fates_history_interface_type) = this.dim_index[:levcdsc]
levcdam_index(this::fates_history_interface_type) = this.dim_index[:levcdam]
levscag_index(this::fates_history_interface_type) = this.dim_index[:levscag]
levscagpft_index(this::fates_history_interface_type) = this.dim_index[:levscagpft]
levagepft_index(this::fates_history_interface_type) = this.dim_index[:levagepft]
levheight_index(this::fates_history_interface_type) = this.dim_index[:levheight]
levelem_index(this::fates_history_interface_type) = this.dim_index[:levelem]
levelpft_index(this::fates_history_interface_type) = this.dim_index[:levelpft]
levelcwd_index(this::fates_history_interface_type) = this.dim_index[:levelcwd]
levelage_index(this::fates_history_interface_type) = this.dim_index[:levelage]
levagefuel_index(this::fates_history_interface_type) = this.dim_index[:levagefuel]
levclscpf_index(this::fates_history_interface_type) = this.dim_index[:levclscpf]
levlanduse_index(this::fates_history_interface_type) = this.dim_index[:levlanduse]
levlulu_index(this::fates_history_interface_type) = this.dim_index[:levlulu]
levlupft_index(this::fates_history_interface_type) = this.dim_index[:levlupft]

num_history_vars(this::fates_history_interface_type) = this.num_history_vars_

# =====================================================================================

"""
    Init!(this, num_threads, fates_bounds)

Set up the `dim_bounds` descriptors (one per FATES IO dimension) and record each
dimension's 1-based slot index. Mirrors the Fortran `Init` (which increments
`dim_count` per dimension and calls `set_<dim>_index` + `dim_bounds(i)%Init`).
"""
function Init!(this::fates_history_interface_type, num_threads::Integer,
               fates_bounds::fates_bounds_type)
    dim_count = 0
    for (sym, dname, lo_field, hi_field) in _DIM_SETUP
        dim_count += 1
        set_dim_index!(this, sym, dim_count)
        lo = getfield(fates_bounds, lo_field)
        hi = getfield(fates_bounds, hi_field)
        Init!(this.dim_bounds[dim_count], dname, num_threads, lo, hi)
    end
    return nothing
end

# =====================================================================================

"""
    SetThreadBoundsEach!(this, thread_index, thread_bounds)

Set the per-thread (clump) bounds of every dimension descriptor. Mirrors the
Fortran `SetThreadBoundsEach`.
"""
function SetThreadBoundsEach!(this::fates_history_interface_type, thread_index::Integer,
                              thread_bounds::fates_bounds_type)
    for (sym, _dname, lo_field, hi_field) in _DIM_SETUP
        idx = this.dim_index[sym]
        lo = getfield(thread_bounds, lo_field)
        hi = getfield(thread_bounds, hi_field)
        SetThreadBounds!(this.dim_bounds[idx], thread_index, lo, hi)
    end
    return nothing
end

# =====================================================================================

"""
    init_dim_kinds_maps!(this)

Initialize the `dim_kinds` registry (29 variable kinds, in Fortran order).
Mirrors `init_dim_kinds_maps`.
"""
function init_dim_kinds_maps!(this::fates_history_interface_type)
    index = 0
    for (name, ndims) in _DIM_KIND_SETUP
        index += 1
        Init!(this.dim_kinds[index], name, ndims)
    end
    return nothing
end

# =====================================================================================

"""
    set_dim_indices!(this, dk_name, idim, dim_index)

Associate dimension `idim` of variable-kind `dk_name` with dimension descriptor
`dim_index`, and set the kind's `dimsize[idim]` from the descriptor bounds.
Mirrors `set_dim_indices`.
"""
function set_dim_indices!(this::fates_history_interface_type, dk_name::AbstractString,
                          idim::Integer, dim_index::Integer)
    ityp = iotype_index(strip(dk_name), fates_history_num_dim_kinds, this.dim_kinds)
    dk = this.dim_kinds[ityp]
    if dk.ndims < idim
        fates_endrun("set_dim_indices!: dim index $idim does not exist for type $dk_name (ndims=$(dk.ndims))")
    end
    if idim == 1
        dk.dim1_index = dim_index
    elseif idim == 2
        dk.dim2_index = dim_index
    end
    dk.dimsize[idim] = this.dim_bounds[dim_index].upper_bound -
                       this.dim_bounds[dim_index].lower_bound + 1
    return nothing
end

# =====================================================================================

"""
    assemble_history_output_types!(this)

Build the dim-kinds registry then wire every kind's dimension indices. Mirrors
`assemble_history_output_types`.
"""
function assemble_history_output_types!(this::fates_history_interface_type)
    init_dim_kinds_maps!(this)
    for (kind, idim, dimsym) in _ASSEMBLE_MAP
        set_dim_indices!(this, kind, idim, this.dim_index[dimsym])
    end
    return nothing
end

# =====================================================================================
# HLM gate evaluation for the registry. Each registry row carries a Vector of
# gate symbols (ANDed). _gate_passes maps a gate symbol to the live hlm flag.
# =====================================================================================

function _gate_passes(g::Symbol)
    if g === :dynam0
        return hlm_hist_level_dynam[] > 0
    elseif g === :dynam1
        return hlm_hist_level_dynam[] > 1
    elseif g === :hifrq0
        return hlm_hist_level_hifrq[] > 0
    elseif g === :hifrq1
        return hlm_hist_level_hifrq[] > 1
    elseif g === :planthydro
        return hlm_use_planthydro[] == itrue
    elseif g === :treedamage
        return hlm_use_tree_damage[] == itrue
    elseif g === :nocomp
        return hlm_use_nocomp[] == itrue
    elseif g === :has_n
        return any(==(nitrogen_element), element_list)
    elseif g === :has_p
        return any(==(phosphorus_element), element_list)
    elseif g === :not_st3_not_sp
        return (hlm_use_ed_st3[] == ifalse) && (hlm_use_sp[] == ifalse)
    else
        error("FatesHistoryInterfaceMod: unknown gate symbol $g")
    end
end

_all_gates_pass(gates::AbstractVector{Symbol}) = all(_gate_passes, gates)

# =====================================================================================

# hlm_name list membership check (Fortran FatesUtilsMod.check_hlm_list).
# `hlms` is a colon-delimited list (e.g. "CLM:ALM"); we register the var if the
# active HLM name appears in it. The HLM name lives in FatesInterfaceTypesMod.
function _check_hlm_list(hlms::AbstractString, hlm::AbstractString)
    isempty(strip(hlm)) && return true   # standalone/unset HLM => register all
    return occursin(strip(hlm), hlms)
end

"""
    set_history_var!(this, vname, units, long, use_default, avgflag, vtype, hlms,
                     upfreq, ivar, initialize, handle; flush_to_zero=false) -> (ivar, index)

Register one history variable. Increments `ivar` and, in `initialize` mode,
`Init!`s its `hvars` slot. Returns the updated counter and the var's index
(0 when the var is not active for this HLM). Mirrors `set_history_var`.
"""
function set_history_var!(this::fates_history_interface_type,
                          vname, units, long, use_default, avgflag, vtype, hlms,
                          upfreq::Integer, ivar::Integer, initialize::Bool,
                          handle::Symbol; flush_to_zero::Bool=false)
    flushval = flush_to_zero ? 0.0 : hlm_hio_ignore_val[]
    write_var = _check_hlm_list(hlms, _hlm_name())
    if write_var
        ivar += 1
        index = ivar
        if initialize
            Init!(this.hvars[ivar], vname, units, long, use_default, vtype, avgflag,
                  flushval, upfreq, fates_history_num_dim_kinds, this.dim_kinds, this.dim_bounds)
        end
    else
        index = 0
    end
    this.ih[handle] = index
    return (ivar, index)
end

# The active HLM name. FATES standalone leaves this unset; default to "" so the
# full registry is built (every var lists CLM:ALM). `hlm_name` is the module-level
# Ref{String} from FatesInterfaceTypesMod.
_hlm_name() = String(hlm_name[])

# =====================================================================================

"""
    define_history_vars!(this, initialize_variables)

Walk the full `_HISTORY_VAR_REGISTRY` (472 vars, source order), registering each
whose HLM gates pass. In count mode (`initialize_variables=false`) it only counts
(`num_history_vars_`); in initialize mode it `Init!`s each `hvars` slot. Mirrors
`define_history_vars`.
"""
function define_history_vars!(this::fates_history_interface_type, initialize_variables::Bool)
    ivar = 0
    for (gates, vname, units, long, use_default, avgflag, vtype, hlms, upfreq, handle, ftz) in _HISTORY_VAR_REGISTRY
        _all_gates_pass(gates) || continue
        (ivar, _idx) = set_history_var!(this, vname, units, long, use_default, avgflag,
                                        vtype, hlms, upfreq, ivar, initialize_variables, handle;
                                        flush_to_zero=ftz)
    end
    this.num_history_vars_ = ivar   # must be last (matches Fortran)
    return nothing
end

# =====================================================================================

"""
    initialize_history_vars!(this)

Two-pass build: count active vars, allocate `hvars`, then initialize them.
Mirrors `initialize_history_vars`.
"""
function initialize_history_vars!(this::fates_history_interface_type)
    define_history_vars!(this, false)                       # count
    this.hvars = [fates_history_variable_type() for _ in 1:this.num_history_vars_]
    define_history_vars!(this, true)                        # initialize
    return nothing
end

# =====================================================================================
# Flush / zero
# =====================================================================================

"""
    flush_hvars!(this, nc, upfreq_in)

Flush every variable in update-group `upfreq_in` back to its flush value.
Mirrors `flush_hvars`.
"""
function flush_hvars!(this::fates_history_interface_type, nc::Integer, upfreq_in::Integer)
    for ivar in eachindex(this.hvars)
        if this.hvars[ivar].upfreq == upfreq_in
            HFlush!(this.hvars[ivar], nc, this.dim_bounds, this.dim_kinds)
        end
    end
    return nothing
end

"""
    flush_all_hvars!(this, nc)

Flush all active update-groups, gated by the live history-level / planthydro
flags. Mirrors `flush_all_hvars`.
"""
function flush_all_hvars!(this::fates_history_interface_type, nc::Integer)
    if hlm_hist_level_hifrq[] > 0
        flush_hvars!(this, nc, group_hifr_simple)
        hlm_use_planthydro[] == itrue && flush_hvars!(this, nc, group_hydr_simple)
        if hlm_hist_level_hifrq[] > 1
            flush_hvars!(this, nc, group_hifr_complx)
            hlm_use_planthydro[] == itrue && flush_hvars!(this, nc, group_hydr_complx)
        end
    end
    if hlm_hist_level_dynam[] > 0
        flush_hvars!(this, nc, group_dyna_simple)
        flush_hvars!(this, nc, group_nflx_simple)
        if hlm_hist_level_dynam[] > 1
            flush_hvars!(this, nc, group_dyna_complx)
            flush_hvars!(this, nc, group_nflx_complx)
        end
    end
    return nothing
end

"""
    zero_site_hvars!(this, currentSite, upfreq_in)

Zero (not flush) the `currentSite` slot of every variable in update-group
`upfreq_in`, prior to filling. Mirrors `zero_site_hvars`. The site is addressed
by `h_gid` (its history output array index).
"""
function zero_site_hvars!(this::fates_history_interface_type, currentSite::ed_site_type,
                          upfreq_in::Integer)
    io_si = currentSite.h_gid
    for ivar in eachindex(this.hvars)
        hv = this.hvars[ivar]
        hv.upfreq == upfreq_in || continue
        dk = this.dim_kinds[hv.dim_kinds_index]
        if strip(dk.name) == site_int
            fates_endrun("zero_site_hvars!: add zeroing provision for SI_INT")
        end
        if dk.ndims == 1
            hv.r81d[io_si] = 0.0
        elseif dk.ndims == 2
            @views hv.r82d[io_si, :] .= 0.0
        elseif dk.ndims == 3
            @views hv.r83d[io_si, :, :] .= 0.0
        end
    end
    return nothing
end

# =====================================================================================
# update_history_* — fill the output buffers from site/patch/cohort state.
# =====================================================================================

"""
    update_history_dyn!(this, nc, nsites, sites, bc_in)

Top-level dynamics (daily) history update. Aborts in ST3 mode; otherwise calls
the level-1 (and, if hist-level>1, level-2) updates. Mirrors `update_history_dyn`.
"""
function update_history_dyn!(this::fates_history_interface_type, nc::Integer,
                             nsites::Integer, sites::AbstractVector,
                             bc_in::AbstractVector)
    hlm_use_ed_st3[] == itrue && return nothing
    if hlm_hist_level_dynam[] > 0
        update_history_dyn1!(this, nc, nsites, sites, bc_in)
        if hlm_hist_level_dynam[] > 1
            update_history_dyn2!(this, nc, nsites, sites, bc_in)
        end
    end
    return nothing
end

"""
    update_history_dyn1!(this, nc, nsites, sites, bc_in)

Level-1 (daily) dynamics history update. Faithful port of the Fortran
`update_history_dyn1` (F90 2353–3024) covering the **W1 core**:

  * counts: NPATCHES/_SECONDARY, NCOHORTS/_SECONDARY
  * patch-area aggregates: TRIMMING, AREA_PLANTS, AREA_TREES, LAI, ELAI,
    secondary-forest fraction/biomass/LAI, the running-mean veg temps
    (TVEG24/TLONGTERM/TGROWTH)
  * site scalars: COLD_STATUS, GDD, NCHILLDAYS, NCOLDDAYS, DAYSINCE_COLDLEAFOFF/_ON,
    CANOPY_SPREAD, FATES_FRACTION, WOODPRODUCT, HARVEST_DEBT/_SEC,
    PRIMARYLAND_FUSION_ERROR, CBAL_ERR_FATES
  * disturbance-rate diagnostics (fire/logging/fall) + harvest/luchange wood-product C flux
  * mortality C-flux + crown-area (canopy/understory), demotion/promotion C flux
  * litter: LITTER_IN (from flux_diags) / LITTER_OUT (fragmentation),
    seed bank / ungerminated / seedling pool / seeds_in (+local)
  * per-cohort PRT biomass pools by element (storec/leafc/fnrtc/reproc/sapwc/totvegc
    + N/P analogues + storetfrac), AGB, BDEAD, BALIVE, canopy/understory biomass, L2FR
  * organ-partitioned NPP fluxes (leaf/seed/stem/froot/croot/stor)
  * basal-area / crown-area weighted height

# TODO Batch 18-followup-W1b: the **fire** group (SPITFIRE_ROS, TFC_ROS,
# FIRE_INTENSITY/_AREA_PRODUCT, FIRE_AREA, FIRE_FUEL_*, SUM_FUEL, NESTEROV,
# NIGNITIONS, FDI, EFFECT_WSPEED) and the **planthydro** (H2OVEG_*) +
# **treedamage** (CROWNAREA_*_DAMAGE) groups are deferred: the patch `fuel`
# object is `nothing` and the fire-state / si_hydr fields are unset (NaN) in the
# carbon-only cold-start path, and the hydro/damage modules are gated off. These
# vars stay at their flush value.
"""
function update_history_dyn1!(this::fates_history_interface_type, nc::Integer,
                              nsites::Integer, sites::AbstractVector,
                              bc_in::AbstractVector)
    # -- counts / patch-area aggregates --
    h_npatches      = this.hvars[this.ih[:ih_npatches_si]].r81d
    h_npatches_sec  = this.hvars[this.ih[:ih_npatches_sec_si]].r81d
    h_ncohorts      = this.hvars[this.ih[:ih_ncohorts_si]].r81d
    h_ncohorts_sec  = this.hvars[this.ih[:ih_ncohorts_sec_si]].r81d
    h_trimming      = this.hvars[this.ih[:ih_trimming_si]].r81d
    h_area_plant    = this.hvars[this.ih[:ih_area_plant_si]].r81d
    h_area_trees    = this.hvars[this.ih[:ih_area_trees_si]].r81d
    h_lai           = this.hvars[this.ih[:ih_lai_si]].r81d
    h_elai          = this.hvars[this.ih[:ih_elai_si]].r81d
    h_tveg24        = this.hvars[this.ih[:ih_tveg24_si]].r81d
    h_tlongterm     = this.hvars[this.ih[:ih_tlongterm_si]].r81d
    h_tgrowth       = this.hvars[this.ih[:ih_tgrowth_si]].r81d
    h_frac_sec      = this.hvars[this.ih[:ih_fraction_secondary_forest_si]].r81d
    h_lai_sec       = this.hvars[this.ih[:ih_lai_secondary_si]].r81d
    h_bio_sec       = this.hvars[this.ih[:ih_biomass_secondary_forest_si]].r81d
    # -- site scalars --
    h_cstatus       = this.hvars[this.ih[:ih_site_cstatus_si]].r81d
    h_gdd           = this.hvars[this.ih[:ih_gdd_si]].r81d
    h_nchill        = this.hvars[this.ih[:ih_site_nchilldays_si]].r81d
    h_ncold         = this.hvars[this.ih[:ih_site_ncolddays_si]].r81d
    h_cleafoff      = this.hvars[this.ih[:ih_cleafoff_si]].r81d
    h_cleafon       = this.hvars[this.ih[:ih_cleafon_si]].r81d
    h_fates_frac    = this.hvars[this.ih[:ih_fates_fraction_si]].r81d
    h_spread        = this.hvars[this.ih[:ih_canopy_spread_si]].r81d
    h_woodprod      = this.hvars[this.ih[:ih_woodproduct_si]].r81d
    h_hdebt         = this.hvars[this.ih[:ih_harvest_debt_si]].r81d
    h_hdebt_sec     = this.hvars[this.ih[:ih_harvest_debt_sec_si]].r81d
    h_pl_fus_err    = this.hvars[this.ih[:ih_primaryland_fusion_error_si]].r81d
    h_cbal_err      = this.hvars[this.ih[:ih_cbal_err_fates_si]].r81d
    h_fire_dist     = this.hvars[this.ih[:ih_fire_disturbance_rate_si]].r81d
    h_log_dist      = this.hvars[this.ih[:ih_logging_disturbance_rate_si]].r81d
    h_fall_dist     = this.hvars[this.ih[:ih_fall_disturbance_rate_si]].r81d
    h_hv_wp_cflux   = this.hvars[this.ih[:ih_harvest_woodprod_carbonflux_si]].r81d
    h_lu_wp_cflux   = this.hvars[this.ih[:ih_luchange_woodprod_carbonflux_si]].r81d
    h_fire_c_to_atm = this.hvars[this.ih[:ih_fire_c_to_atm_si]].r81d
    # -- mortality / demotion / promotion fluxes --
    h_can_mort_cflux = this.hvars[this.ih[:ih_canopy_mortality_carbonflux_si]].r81d
    h_ust_mort_cflux = this.hvars[this.ih[:ih_understory_mortality_carbonflux_si]].r81d
    h_can_mort_ca    = this.hvars[this.ih[:ih_canopy_mortality_crownarea_si]].r81d
    h_ust_mort_ca    = this.hvars[this.ih[:ih_understory_mortality_crownarea_si]].r81d
    h_demote_cflux   = this.hvars[this.ih[:ih_demotion_carbonflux_si]].r81d
    h_promote_cflux  = this.hvars[this.ih[:ih_promotion_carbonflux_si]].r81d
    # -- biomass aggregates --
    h_bdead         = this.hvars[this.ih[:ih_bdead_si]].r81d
    h_balive        = this.hvars[this.ih[:ih_balive_si]].r81d
    h_agb           = this.hvars[this.ih[:ih_agb_si]].r81d
    h_can_bio       = this.hvars[this.ih[:ih_canopy_biomass_si]].r81d
    h_ust_bio       = this.hvars[this.ih[:ih_understory_biomass_si]].r81d
    h_baw_height    = this.hvars[this.ih[:ih_ba_weighted_height_si]].r81d
    h_caw_height    = this.hvars[this.ih[:ih_ca_weighted_height_si]].r81d
    # -- per-cohort C pools --
    h_storec        = this.hvars[this.ih[:ih_storec_si]].r81d
    h_leafc         = this.hvars[this.ih[:ih_leafc_si]].r81d
    h_fnrtc         = this.hvars[this.ih[:ih_fnrtc_si]].r81d
    h_reproc        = this.hvars[this.ih[:ih_reproc_si]].r81d
    h_sapwc         = this.hvars[this.ih[:ih_sapwc_si]].r81d
    h_totvegc       = this.hvars[this.ih[:ih_totvegc_si]].r81d
    h_storectfrac   = this.hvars[this.ih[:ih_storectfrac_si]].r81d
    h_l2fr          = this.hvars[this.ih[:ih_l2fr_si]].r81d
    # -- organ-partitioned NPP --
    h_npp_leaf      = this.hvars[this.ih[:ih_npp_leaf_si]].r81d
    h_npp_seed      = this.hvars[this.ih[:ih_npp_seed_si]].r81d
    h_npp_stem      = this.hvars[this.ih[:ih_npp_stem_si]].r81d
    h_npp_froot     = this.hvars[this.ih[:ih_npp_froot_si]].r81d
    h_npp_croot     = this.hvars[this.ih[:ih_npp_croot_si]].r81d
    h_npp_stor      = this.hvars[this.ih[:ih_npp_stor_si]].r81d
    # -- litter / seed pools --
    h_litter_in     = this.hvars[this.ih[:ih_litter_in_si]].r81d
    h_litter_out    = this.hvars[this.ih[:ih_litter_out_si]].r81d
    h_seed_bank     = this.hvars[this.ih[:ih_seed_bank_si]].r81d
    h_ungerm_seed   = this.hvars[this.ih[:ih_ungerm_seed_bank_si]].r81d
    h_seedling      = this.hvars[this.ih[:ih_seedling_pool_si]].r81d
    h_seeds_in      = this.hvars[this.ih[:ih_seeds_in_si]].r81d
    h_seeds_in_loc  = this.hvars[this.ih[:ih_seeds_in_local_si]].r81d

    cpos = element_pos[carbon12_element]   # carbon's slot in mass_balance / litter / flux_diags

    for s in 1:nsites
        site = sites[s]
        io_si = site.h_gid

        site_ba = 0.0
        site_ca = 0.0

        # zero this site's dynamics-group slots before filling
        zero_site_hvars!(this, site, group_dyna_simple)

        # FATES fraction (1 on fates columns; gridcell-average gives the fraction)
        h_fates_frac[io_si] = 1.0

        # Site mass-balance / flux-diagnostic arrays are only allocated once the
        # element registry + cold-start chain has run (the live use_fates path).
        # The port's empty-array convention marks "unallocated"; guard the deep
        # reads so the registry/aggregation unit path (synthetic bare site) is
        # still exercisable. In the live run these are always populated.
        has_massbal  = cpos >= 1 && length(site.mass_balance) >= cpos
        has_fluxdiag = cpos >= 1 && length(site.flux_diags.elem) >= cpos

        # carbon-model error [kgC/day -> kgC/s]
        if has_massbal
            h_cbal_err[io_si] = site.mass_balance[cpos].err_fates / sec_per_day
            # carbon lost to atmosphere from burning [kgC/site/day -> kgC/m2/s]
            h_fire_c_to_atm[io_si] = site.mass_balance[cpos].burn_flux_to_atm *
                                     ha_per_m2 * days_per_sec
        end

        # site scalar state (Fortran casts integer counters to real)
        h_spread[io_si]  = site.spread
        h_cstatus[io_si] = Float64(site.cstatus)
        h_nchill[io_si]  = Float64(site.nchilldays)
        h_ncold[io_si]   = Float64(site.ncolddays)
        h_gdd[io_si]     = site.grow_deg_days
        h_cleafoff[io_si] = Float64(site.phen_model_date - site.cleafoffdate)
        h_cleafon[io_si]  = Float64(site.phen_model_date - site.cleafondate)

        # total wood product accumulation [kgC/site -> kgC/m2]
        h_woodprod[io_si] = site.resources_management.trunk_product_site * AREA_INV

        h_hdebt[io_si]     = site.resources_management.harvest_debt
        h_hdebt_sec[io_si] = site.resources_management.harvest_debt_sec

        # primary-land patch-fusion error [m2/m2/day -> m2/m2/yr] (finite after zero_site!)
        if isfinite(site.primary_land_patchfusion_error)
            h_pl_fus_err[io_si] = site.primary_land_patchfusion_error * days_per_year
        end

        # site-level disturbance rates [m2/m2/day -> m2/m2/yr] (always allocated)
        h_fire_dist[io_si] = sum(@view site.disturbance_rates[dtype_ifire, :, :]) * days_per_year
        h_log_dist[io_si]  = sum(@view site.disturbance_rates[dtype_ilog,  :, :]) * days_per_year
        h_fall_dist[io_si] = sum(@view site.disturbance_rates[dtype_ifall, :, :]) * days_per_year

        # harvest / land-use-change wood-product C flux [kgC/m2]
        if has_massbal && length(site.mass_balance[cpos].wood_product_harvest) >= numpft[]
            h_hv_wp_cflux[io_si] = AREA_INV * sum(@view site.mass_balance[cpos].wood_product_harvest[1:numpft[]])
            h_lu_wp_cflux[io_si] = AREA_INV * sum(@view site.mass_balance[cpos].wood_product_landusechange[1:numpft[]])
        end

        # mortality C flux + crown-area (site mort arrays allocated by zero_site!)
        if !isempty(site.fmort_carbonflux_canopy)
            # fire (gC->kgC) + impact + termination (kgC/ha/day -> kgC/m2/s)
            h_can_mort_cflux[io_si] += sum(site.fmort_carbonflux_canopy) / g_per_kg
            h_ust_mort_cflux[io_si] += sum(site.fmort_carbonflux_ustory) / g_per_kg
            h_ust_mort_cflux[io_si] += sum(site.imort_carbonflux)
            h_demote_cflux[io_si]   = site.demotion_carbonflux * ha_per_m2 * days_per_sec
            h_promote_cflux[io_si]  = site.promotion_carbonflux * ha_per_m2 * days_per_sec
            h_can_mort_cflux[io_si] += sum(site.term_carbonflux_canopy) * days_per_sec * ha_per_m2
            h_ust_mort_cflux[io_si] += sum(site.term_carbonflux_ustory) * days_per_sec * ha_per_m2

            h_can_mort_ca[io_si] += site.fmort_crownarea_canopy + site.term_crownarea_canopy * days_per_year
            h_ust_mort_ca[io_si] += site.fmort_crownarea_ustory + site.term_crownarea_ustory * days_per_year +
                                    site.imort_crownarea
        end

        # litter-input flux (from carbon flux diagnostics) [kgC/m2/s]
        if has_fluxdiag
            edc = site.flux_diags.elem[cpos]
            h_litter_in[io_si] = (sum(edc.cwd_ag_input) + sum(edc.cwd_bg_input) +
                                  sum(edc.surf_fine_litter_input) + sum(edc.root_litter_input)) *
                                 AREA_INV * days_per_sec
        end

        # patch / cohort aggregation
        cpatch = site.oldest_patch
        while cpatch !== nothing
            h_npatches[io_si] += 1.0
            cpatch.land_use_label == secondaryland && (h_npatches_sec[io_si] += 1.0)

            # total_canopy_area is left NaN by the carbon-only cold-start path on
            # a patch whose canopy structure hasn't been summarized; treat an
            # unset value as zero canopy (the port's "NaN == unset" convention),
            # so the area/LAI aggregates stay finite.
            tcanp = isfinite(cpatch.total_canopy_area) ? cpatch.total_canopy_area : 0.0
            ttree = isfinite(cpatch.total_tree_area)   ? cpatch.total_tree_area   : 0.0
            # LAI / ELAI = canopy-area-weighted profile integral [m2/m2]
            if !isempty(cpatch.canopy_area_profile)
                h_lai[io_si]  += sum(cpatch.canopy_area_profile .* cpatch.tlai_profile) * tcanp * AREA_INV
                h_elai[io_si] += sum(cpatch.canopy_area_profile .* cpatch.elai_profile) * tcanp * AREA_INV
            end

            # running-mean veg temperatures, patch-area weighted [degC]
            aw = cpatch.area * AREA_INV
            if cpatch.tveg24 !== nothing
                h_tveg24[io_si]    += (GetMean(cpatch.tveg24)        - t_water_freeze_k_1atm) * aw
            end
            if cpatch.tveg_longterm !== nothing
                h_tlongterm[io_si] += (GetMean(cpatch.tveg_longterm) - t_water_freeze_k_1atm) * aw
            end
            if cpatch.tveg_lpa !== nothing
                h_tgrowth[io_si]   += (GetMean(cpatch.tveg_lpa)      - t_water_freeze_k_1atm) * aw
            end

            # secondary-forest diagnostics
            if cpatch.land_use_label == secondaryland
                h_frac_sec[io_si] += cpatch.area * AREA_INV
                if !isempty(cpatch.tlai_profile)
                    h_lai_sec[io_si] += sum(cpatch.tlai_profile) * tcanp
                end
            end

            # canopy trimming (area-weighted, tallest cohort's trim)
            if cpatch.tallest !== nothing
                h_trimming[io_si] += cpatch.tallest.canopy_trim * cpatch.area * AREA_INV
            end

            # plant / tree area occupied [m2/m2]
            h_area_plant[io_si] += min(tcanp, cpatch.area) * AREA_INV
            h_area_trees[io_si] += min(ttree, cpatch.area) * AREA_INV

            # ---- litter / seed pools (carbon) ----
            if has_fluxdiag && length(cpatch.litter) >= cpos
                litt = cpatch.litter[cpos]
                area_frac = cpatch.area * AREA_INV
                h_litter_out[io_si] += (sum(litt.leaf_fines_frag) + sum(litt.root_fines_frag) +
                                        sum(litt.ag_cwd_frag) + sum(litt.bg_cwd_frag) +
                                        sum(litt.seed_decay) + sum(litt.seed_germ_decay)) *
                                       area_frac * days_per_sec
                h_seed_bank[io_si]   += (sum(litt.seed) + sum(litt.seed_germ)) * area_frac
                h_ungerm_seed[io_si] += sum(litt.seed) * area_frac
                h_seedling[io_si]    += sum(litt.seed_germ) * area_frac
                h_seeds_in[io_si]    += (sum(litt.seed_in_local) + sum(litt.seed_in_extern)) *
                                        area_frac * days_per_sec
                h_seeds_in_loc[io_si] += sum(litt.seed_in_local) * area_frac * days_per_sec
            end

            # ---- cohort loop (shortest -> taller, matching Fortran) ----
            ccohort = cpatch.shortest
            while ccohort !== nothing
                ft = ccohort.pft
                n_perm2 = ccohort.n * AREA_INV

                h_ncohorts[io_si] += 1.0
                cpatch.land_use_label == secondaryland && (h_ncohorts_sec[io_si] += 1.0)

                # Cohort PRT pools are only allocated on cohorts built by the
                # cold-start chain (InitPRTObject). Skip the biomass/NPP/mortality
                # fills for a bare unit-test cohort (prt === nothing); counts above
                # are still accumulated.
                if ccohort.prt === nothing
                    ccohort = ccohort.taller
                    continue
                end

                # biomass pools by element
                local total_m = 0.0
                local store_m = 0.0
                for el in 1:num_elements[]
                    eid    = element_list[el]
                    sapw_m   = GetState(ccohort.prt, sapw_organ,   eid)
                    struct_m = GetState(ccohort.prt, struct_organ, eid)
                    leaf_m   = GetState(ccohort.prt, leaf_organ,   eid)
                    fnrt_m   = GetState(ccohort.prt, fnrt_organ,   eid)
                    store_e  = GetState(ccohort.prt, store_organ,  eid)
                    repro_m  = GetState(ccohort.prt, repro_organ,  eid)
                    alive_m  = leaf_m + fnrt_m + sapw_m
                    tot_m    = alive_m + store_e + struct_m

                    if eid == carbon12_element
                        total_m = tot_m
                        store_m = store_e
                        h_storec[io_si]  += ccohort.n * store_e / m2_per_ha
                        h_leafc[io_si]   += ccohort.n * leaf_m  / m2_per_ha
                        h_fnrtc[io_si]   += ccohort.n * fnrt_m  / m2_per_ha
                        h_reproc[io_si]  += ccohort.n * repro_m / m2_per_ha
                        h_sapwc[io_si]   += ccohort.n * sapw_m  / m2_per_ha
                        h_totvegc[io_si] += ccohort.n * tot_m   / m2_per_ha

                        store_max, _ = bstore_allom(ccohort.dbh, ccohort.pft,
                                                    ccohort.crowndamage, ccohort.canopy_trim)
                        h_storectfrac[io_si] += ccohort.n * store_max / m2_per_ha

                        h_bdead[io_si]  += n_perm2 * struct_m
                        h_balive[io_si] += n_perm2 * alive_m
                        h_agb[io_si]    += n_perm2 *
                            (leaf_m + (sapw_m + struct_m + store_e) * prt_params.allom_agb_frac[ft])

                        if cpatch.land_use_label == secondaryland
                            h_bio_sec[io_si] += tot_m * ccohort.n * AREA_INV
                        end

                        if hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
                            h_l2fr[io_si] += ccohort.l2fr * ccohort.n * fnrt_m / m2_per_ha
                        else
                            h_l2fr[io_si] += prt_params.allom_l2fr[ft] * ccohort.n * fnrt_m / m2_per_ha
                        end
                    elseif eid == nitrogen_element
                        store_max = GetNutrientTarget(ccohort.prt, eid, store_organ, stoich_growth_min)
                        this.hvars[this.ih[:ih_storen_si]].r81d[io_si]      += ccohort.n * store_e / m2_per_ha
                        this.hvars[this.ih[:ih_storentfrac_si]].r81d[io_si] += ccohort.n * store_max / m2_per_ha
                        this.hvars[this.ih[:ih_leafn_si]].r81d[io_si]       += ccohort.n * leaf_m / m2_per_ha
                        this.hvars[this.ih[:ih_fnrtn_si]].r81d[io_si]       += ccohort.n * fnrt_m / m2_per_ha
                        this.hvars[this.ih[:ih_repron_si]].r81d[io_si]      += ccohort.n * repro_m / m2_per_ha
                        this.hvars[this.ih[:ih_sapwn_si]].r81d[io_si]       += ccohort.n * sapw_m / m2_per_ha
                        this.hvars[this.ih[:ih_totvegn_si]].r81d[io_si]     += ccohort.n * tot_m / m2_per_ha
                    elseif eid == phosphorus_element
                        store_max = GetNutrientTarget(ccohort.prt, eid, store_organ, stoich_growth_min)
                        this.hvars[this.ih[:ih_storep_si]].r81d[io_si]      += ccohort.n * store_e / m2_per_ha
                        this.hvars[this.ih[:ih_storeptfrac_si]].r81d[io_si] += ccohort.n * store_max / m2_per_ha
                        this.hvars[this.ih[:ih_leafp_si]].r81d[io_si]       += ccohort.n * leaf_m / m2_per_ha
                        this.hvars[this.ih[:ih_fnrtp_si]].r81d[io_si]       += ccohort.n * fnrt_m / m2_per_ha
                        this.hvars[this.ih[:ih_reprop_si]].r81d[io_si]      += ccohort.n * repro_m / m2_per_ha
                        this.hvars[this.ih[:ih_sapwp_si]].r81d[io_si]       += ccohort.n * sapw_m / m2_per_ha
                        this.hvars[this.ih[:ih_totvegp_si]].r81d[io_si]     += ccohort.n * tot_m / m2_per_ha
                    end
                end

                # FLUXES — only meaningful after a cohort has lived a day
                if !ccohort.isnew
                    # organ net-alloc [kgC/day] * [day/yr] = [kgC/yr]
                    sapw_na   = GetNetAlloc(ccohort.prt, sapw_organ,   carbon12_element) * days_per_year
                    store_na  = GetNetAlloc(ccohort.prt, store_organ,  carbon12_element) * days_per_year
                    leaf_na   = GetNetAlloc(ccohort.prt, leaf_organ,   carbon12_element) * days_per_year
                    fnrt_na   = GetNetAlloc(ccohort.prt, fnrt_organ,   carbon12_element) * days_per_year
                    struct_na = GetNetAlloc(ccohort.prt, struct_organ, carbon12_element) * days_per_year
                    repro_na  = GetNetAlloc(ccohort.prt, repro_organ,  carbon12_element) * days_per_year

                    agbf = prt_params.allom_agb_frac[ft]
                    # [kgC/yr] -> [kgC/m2/s]
                    h_npp_leaf[io_si]  += leaf_na  * n_perm2 / days_per_year / sec_per_day
                    h_npp_seed[io_si]  += repro_na * n_perm2 / days_per_year / sec_per_day
                    h_npp_stem[io_si]  += (sapw_na + struct_na) * n_perm2 * agbf / days_per_year / sec_per_day
                    h_npp_froot[io_si] += fnrt_na  * n_perm2 / days_per_year / sec_per_day
                    h_npp_croot[io_si] += (sapw_na + struct_na) * n_perm2 * (1.0 - agbf) / days_per_year / sec_per_day
                    h_npp_stor[io_si]  += store_na * n_perm2 / days_per_year / sec_per_day

                    # basal-area / crown-area weighted height
                    if prt_params.woody[ft] == itrue
                        cohort_ba = 0.25 * pi_const * ((ccohort.dbh / 100.0)^2) * ccohort.n
                        h_baw_height[io_si] += ccohort.height * cohort_ba
                        site_ba += cohort_ba
                    end
                    h_caw_height[io_si] += ccohort.height * ccohort.c_area / m2_per_ha
                    site_ca += ccohort.c_area / m2_per_ha

                    # canopy- / understory-separated biomass + mortality fluxes
                    sum_mort = ccohort.bmort + ccohort.hmort + ccohort.cmort +
                               ccohort.frmort + ccohort.smort + ccohort.asmort + ccohort.dgmort
                    sum_lmort = ccohort.lmort_direct + ccohort.lmort_collateral + ccohort.lmort_infra
                    if ccohort.canopy_layer == 1
                        h_can_bio[io_si] += n_perm2 * total_m
                        h_can_mort_cflux[io_si] +=
                            sum_mort * total_m * ccohort.n * days_per_sec * years_per_day * ha_per_m2 +
                            sum_lmort * total_m * ccohort.n * ha_per_m2
                        h_can_mort_ca[io_si] +=
                            sum_mort * ccohort.c_area +
                            sum_lmort * ccohort.c_area * sec_per_day * days_per_year
                    else
                        h_ust_bio[io_si] += n_perm2 * total_m
                        h_ust_mort_cflux[io_si] +=
                            sum_mort * total_m * ccohort.n * days_per_sec * years_per_day * ha_per_m2 +
                            sum_lmort * total_m * ccohort.n * ha_per_m2
                        h_ust_mort_ca[io_si] +=
                            sum_mort * ccohort.c_area +
                            sum_lmort * ccohort.c_area * sec_per_day * days_per_year
                    end
                end

                ccohort = ccohort.taller
            end

            cpatch = cpatch.younger
        end

        # ---- normalizations ----
        if h_frac_sec[io_si] > nearzero
            h_lai_sec[io_si] /= (h_frac_sec[io_si] * area)
        end
        if site_ca > nearzero
            h_caw_height[io_si] /= site_ca
        end
        if site_ba > nearzero
            h_baw_height[io_si] /= site_ba
        end

        # storage-fraction = stored / target (per element)
        for el in 1:num_elements[]
            eid = element_list[el]
            if eid == carbon12_element && h_storectfrac[io_si] > nearzero
                h_storectfrac[io_si] = h_storec[io_si] / h_storectfrac[io_si]
            elseif eid == nitrogen_element
                hn = this.hvars[this.ih[:ih_storentfrac_si]].r81d
                hn[io_si] > nearzero && (hn[io_si] = this.hvars[this.ih[:ih_storen_si]].r81d[io_si] / hn[io_si])
            elseif eid == phosphorus_element
                hp = this.hvars[this.ih[:ih_storeptfrac_si]].r81d
                hp[io_si] > nearzero && (hp[io_si] = this.hvars[this.ih[:ih_storep_si]].r81d[io_si] / hp[io_si])
            end
        end

        if h_fnrtc[io_si] > nearzero
            h_l2fr[io_si] /= h_fnrtc[io_si]
        else
            h_l2fr[io_si] = hlm_hio_ignore_val[]
        end

        # zero the site-level termination/damage accumulators (consumed)
        fill!(site.term_carbonflux_canopy, 0.0)
        fill!(site.term_carbonflux_ustory, 0.0)
        site.crownarea_canopy_damage = 0.0
        site.crownarea_ustory_damage = 0.0
    end
    return nothing
end

"""
    update_history_dyn2!(this, nc, nsites, sites, bc_in)

Level-2 (complex-dimension) dynamics update — the per-scpf/scag/cnlf/element
disaggregated fills.

# TODO Batch 18-followup: not yet ported. Requires per-cohort PRT pools and the
# size/age/element class-index helpers wired against live FATES state. The
# variable registry for these (site_size_pft_r8, site_scag_r8, site_elem_r8, ...)
# IS fully registered; only the buffer-fill logic is deferred.
"""
function update_history_dyn2!(this::fates_history_interface_type, nc::Integer,
                              nsites::Integer, sites::AbstractVector,
                              bc_in::AbstractVector)
    # TODO Batch 18-followup: port the ~1850-line per-class aggregation.
    return nothing
end

"""
    update_history_nutrflux!(this, csite)

Nutrient-flux (N/P uptake, efflux, demand, fixation) dynamics update.

# TODO Batch 18-followup: not yet ported (requires PRT CNP flux diagnostics).
"""
function update_history_nutrflux!(this::fates_history_interface_type, csite)
    # TODO Batch 18-followup
    return nothing
end

"""
    update_history_hifrq!(this, nc, nsites, sites, bc_in, bc_out, dt_tstep)

Top-level high-frequency (sub-daily) history update orchestrator. Mirrors
`update_history_hifrq`.

# TODO Batch 18-followup: hifrq1/hifrq2 buffer fills (GPP/AR/radiation/temperature
# per cnlf/can/scpf) not yet ported; the variable registry is complete.
"""
function update_history_hifrq!(this::fates_history_interface_type, nc::Integer,
                               nsites::Integer, sites::AbstractVector,
                               bc_in::AbstractVector, bc_out::AbstractVector, dt_tstep::Real)
    if hlm_hist_level_hifrq[] > 0
        update_history_hifrq1!(this, nc, nsites, sites, bc_in, bc_out, dt_tstep)
        if hlm_hist_level_hifrq[] > 1
            update_history_hifrq2!(this, nc, nsites, sites, bc_in, bc_out, dt_tstep)
        end
    end
    return nothing
end

"""
    update_history_hifrq1!(this, nc, nsites, sites, bc_in, bc_out, dt_tstep)

Level-1 high-frequency (per-timestep) history update. Faithful port of the
Fortran `update_history_hifrq1` (F90 4913–5170): site-level sub-daily carbon
fluxes from the per-step cohort state.

Filled:
  * NEP/HR (from `bc_in.tot_het_resp`), GPP/NPP/ARESP (+ _SECONDARY),
    GROWTH_RESP / MAINT_RESP (+ _SECONDARY) / MAINT_RESP_UNREDUCED
  * organ MR: LEAF_MR (rdark), FROOT_MR, LIVECROOT_MR, LIVESTEM_MR
  * canopy/understory GPP & AR
  * C_STOMATA / C_LBLAYER (patch conductances, canopy-area weighted)
  * TVEG (instantaneous veg temp from `bc_in.t_veg_pa`)
  * VIS/NIR radiation conservation error (age-area weighted across patches)

Where a per-step source is unavailable in the carbon-only path it is handled
exactly as Fortran would with a zero/ignore value: `tot_het_resp` defaults to
zero (no use_fates_bgc heterotrophic respiration), so NEP == GPP-AR and HR == 0.
"""
function update_history_hifrq1!(this::fates_history_interface_type, nc, nsites, sites, bc_in, bc_out, dt_tstep)
    h_gpp        = this.hvars[this.ih[:ih_gpp_si]].r81d
    h_gpp_sec    = this.hvars[this.ih[:ih_gpp_secondary_si]].r81d
    h_npp        = this.hvars[this.ih[:ih_npp_si]].r81d
    h_npp_sec    = this.hvars[this.ih[:ih_npp_secondary_si]].r81d
    h_aresp      = this.hvars[this.ih[:ih_aresp_si]].r81d
    h_aresp_sec  = this.hvars[this.ih[:ih_aresp_secondary_si]].r81d
    h_maint      = this.hvars[this.ih[:ih_maint_resp_si]].r81d
    h_maint_sec  = this.hvars[this.ih[:ih_maint_resp_secondary_si]].r81d
    h_growth     = this.hvars[this.ih[:ih_growth_resp_si]].r81d
    h_growth_sec = this.hvars[this.ih[:ih_growth_resp_secondary_si]].r81d
    h_cstom      = this.hvars[this.ih[:ih_c_stomata_si]].r81d
    h_clbl       = this.hvars[this.ih[:ih_c_lblayer_si]].r81d
    h_visrad_err = this.hvars[this.ih[:ih_vis_rad_err_si]].r81d
    h_nirrad_err = this.hvars[this.ih[:ih_nir_rad_err_si]].r81d
    h_nep        = this.hvars[this.ih[:ih_nep_si]].r81d
    h_hr         = this.hvars[this.ih[:ih_hr_si]].r81d
    h_gpp_can    = this.hvars[this.ih[:ih_gpp_canopy_si]].r81d
    h_ar_can     = this.hvars[this.ih[:ih_ar_canopy_si]].r81d
    h_gpp_ust    = this.hvars[this.ih[:ih_gpp_understory_si]].r81d
    h_ar_ust     = this.hvars[this.ih[:ih_ar_understory_si]].r81d
    h_leaf_mr    = this.hvars[this.ih[:ih_leaf_mr_si]].r81d
    h_froot_mr   = this.hvars[this.ih[:ih_froot_mr_si]].r81d
    h_livecr_mr  = this.hvars[this.ih[:ih_livecroot_mr_si]].r81d
    h_livest_mr  = this.hvars[this.ih[:ih_livestem_mr_si]].r81d
    h_maint_unr  = this.hvars[this.ih[:ih_maint_resp_unreduced_si]].r81d
    h_tveg       = this.hvars[this.ih[:ih_tveg_si]].r81d

    flush_hvars!(this, nc, group_hifr_simple)
    dt_tstep_inv = 1.0 / dt_tstep

    for s in 1:nsites
        site = sites[s]
        bcin = bc_in[s]
        io_si = site.h_gid

        zero_site_hvars!(this, site, group_hifr_simple)

        thr = isfinite(bcin.tot_het_resp) ? bcin.tot_het_resp : 0.0
        h_nep[io_si] = -thr * kg_per_g
        h_hr[io_si]  =  thr * kg_per_g

        # ---- radiation-error diagnostics (age-area weighted over solved patches) ----
        # Only patches whose VIS radiation solver ran carry a finite rad_error.
        sum_area_rad = 0.0
        cpatch = site.oldest_patch
        while cpatch !== nothing
            rv = cpatch.rad_error[ivis]
            if isfinite(rv) && abs(rv) > nearzero
                sum_area_rad += cpatch.total_canopy_area
            end
            cpatch = cpatch.younger
        end

        if sum_area_rad < nearzero
            h_visrad_err[io_si] = hlm_hio_ignore_val[]
            h_nirrad_err[io_si] = hlm_hio_ignore_val[]
        else
            h_visrad_err[io_si] = 0.0
            h_nirrad_err[io_si] = 0.0
            cpatch = site.oldest_patch
            while cpatch !== nothing
                rv = cpatch.rad_error[ivis]
                if isfinite(rv) && abs(rv) > nearzero
                    w = cpatch.total_canopy_area / sum_area_rad
                    h_visrad_err[io_si] += cpatch.rad_error[ivis] * w
                    h_nirrad_err[io_si] += cpatch.rad_error[inir] * w
                end
                cpatch = cpatch.younger
            end
        end

        # ---- vegetated-area normalizer ----
        if hlm_use_nocomp[] == itrue && hlm_use_fixed_biogeog[] == itrue
            site_area_veg = area - site.area_bareground * area
        else
            site_area_veg = 0.0
            cpatch = site.oldest_patch
            while cpatch !== nothing
                if isfinite(cpatch.total_canopy_area)
                    site_area_veg += cpatch.total_canopy_area
                end
                cpatch = cpatch.younger
            end
        end

        if site_area_veg < nearzero
            h_cstom[io_si] = hlm_hio_ignore_val[]
            h_clbl[io_si]  = hlm_hio_ignore_val[]
            h_tveg[io_si]  = hlm_hio_ignore_val[]
            continue
        end

        site_area_veg_inv = 1.0 / site_area_veg

        cpatch = site.oldest_patch
        while cpatch !== nothing
            tcanp = isfinite(cpatch.total_canopy_area) ? cpatch.total_canopy_area : 0.0
            if isfinite(cpatch.c_stomata)
                h_cstom[io_si] += cpatch.c_stomata * tcanp * mol_per_umol * site_area_veg_inv
            end
            if isfinite(cpatch.c_lblayer)
                h_clbl[io_si]  += cpatch.c_lblayer * tcanp * mol_per_umol * site_area_veg_inv
            end

            # instantaneous veg temperature for vegetated patches only
            if cpatch.patchno != 0 && cpatch.patchno >= 1 &&
               length(bcin.t_veg_pa) >= cpatch.patchno
                h_tveg[io_si] += (bcin.t_veg_pa[cpatch.patchno] - t_water_freeze_k_1atm) *
                                 tcanp * site_area_veg_inv
            end

            is_sec = cpatch.land_use_label == secondaryland

            ccohort = cpatch.shortest
            while ccohort !== nothing
                if ccohort.isnew
                    ccohort = ccohort.taller
                    continue
                end
                n_perm2 = ccohort.n * AREA_INV

                # [kg/plant/timestep] -> [kg/m2/s]
                h_npp[io_si]   += ccohort.npp_tstep   * n_perm2 * dt_tstep_inv
                h_nep[io_si]   += ccohort.npp_tstep   * n_perm2 * dt_tstep_inv
                h_gpp[io_si]   += ccohort.gpp_tstep   * n_perm2 * dt_tstep_inv
                h_aresp[io_si] += ccohort.resp_tstep  * n_perm2 * dt_tstep_inv
                h_growth[io_si] += ccohort.resp_g_tstep * n_perm2 * dt_tstep_inv
                h_maint[io_si] += ccohort.resp_m      * n_perm2 * dt_tstep_inv
                h_maint_unr[io_si] += ccohort.resp_m_unreduced * n_perm2 * dt_tstep_inv

                if is_sec
                    h_npp_sec[io_si]   += ccohort.npp_tstep   * n_perm2 * dt_tstep_inv
                    h_gpp_sec[io_si]   += ccohort.gpp_tstep   * n_perm2 * dt_tstep_inv
                    h_aresp_sec[io_si] += ccohort.resp_tstep  * n_perm2 * dt_tstep_inv
                    h_growth_sec[io_si] += ccohort.resp_g_tstep * n_perm2 * dt_tstep_inv
                    h_maint_sec[io_si] += ccohort.resp_m      * n_perm2 * dt_tstep_inv
                end

                # organ maintenance respiration [kg/plant/s] -> [kg/m2/s]
                h_leaf_mr[io_si]   += ccohort.rdark        * n_perm2
                h_froot_mr[io_si]  += ccohort.froot_mr     * n_perm2
                h_livecr_mr[io_si] += ccohort.livecroot_mr * n_perm2
                h_livest_mr[io_si] += ccohort.livestem_mr  * n_perm2

                if ccohort.canopy_layer == 1
                    h_gpp_can[io_si] += ccohort.gpp_tstep  * n_perm2 * dt_tstep_inv
                    h_ar_can[io_si]  += ccohort.resp_tstep * n_perm2 * dt_tstep_inv
                else
                    h_gpp_ust[io_si] += ccohort.gpp_tstep  * n_perm2 * dt_tstep_inv
                    h_ar_ust[io_si]  += ccohort.resp_tstep * n_perm2 * dt_tstep_inv
                end

                ccohort = ccohort.taller
            end
            cpatch = cpatch.younger
        end
    end
    return nothing
end

function update_history_hifrq2!(this::fates_history_interface_type, nc, nsites, sites, bc_in, bc_out, dt_tstep)
    # TODO Batch 18-followup
    return nothing
end

"""
    update_history_hydraulics!(this, nc, nsites, sites, bc_in, dt_tstep)

Plant-hydraulics high-frequency history update.

# TODO Batch 18-followup: not yet ported (requires the hydraulics state si_hydr/
# co_hydr). The registry's hydro vars (group_hydr_*) are fully registered.
"""
function update_history_hydraulics!(this::fates_history_interface_type, nc, nsites, sites, bc_in, dt_tstep)
    # TODO Batch 18-followup
    return nothing
end
