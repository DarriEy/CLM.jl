# ==========================================================================
# Ported from: src/biogeochem/CNVegNitrogenFluxType.F90
# Vegetation nitrogen flux data type allocation and initialization
# ==========================================================================

"""
    CNVegNitrogenFluxData

Vegetation nitrogen flux data structure. Holds nitrogen flux variables at patch,
column, and gridcell levels.

Ported from `cnveg_nitrogenflux_type` in `CNVegNitrogenFluxType.F90`.
"""
Base.@kwdef mutable struct CNVegNitrogenFluxData
    # --- Gap mortality fluxes (gN/m2/s) patch-level ---
    m_leafn_to_litter_patch                   ::Vector{Float64} = Float64[]
    m_frootn_to_litter_patch                  ::Vector{Float64} = Float64[]
    m_leafn_storage_to_litter_patch           ::Vector{Float64} = Float64[]
    m_frootn_storage_to_litter_patch          ::Vector{Float64} = Float64[]
    m_livestemn_storage_to_litter_patch       ::Vector{Float64} = Float64[]
    m_deadstemn_storage_to_litter_patch       ::Vector{Float64} = Float64[]
    m_livecrootn_storage_to_litter_patch      ::Vector{Float64} = Float64[]
    m_deadcrootn_storage_to_litter_patch      ::Vector{Float64} = Float64[]
    m_leafn_xfer_to_litter_patch              ::Vector{Float64} = Float64[]
    m_frootn_xfer_to_litter_patch             ::Vector{Float64} = Float64[]
    m_livestemn_xfer_to_litter_patch          ::Vector{Float64} = Float64[]
    m_deadstemn_xfer_to_litter_patch          ::Vector{Float64} = Float64[]
    m_livecrootn_xfer_to_litter_patch         ::Vector{Float64} = Float64[]
    m_deadcrootn_xfer_to_litter_patch         ::Vector{Float64} = Float64[]
    m_livestemn_to_litter_patch               ::Vector{Float64} = Float64[]
    m_deadstemn_to_litter_patch               ::Vector{Float64} = Float64[]
    m_livecrootn_to_litter_patch              ::Vector{Float64} = Float64[]
    m_deadcrootn_to_litter_patch              ::Vector{Float64} = Float64[]
    m_retransn_to_litter_patch                ::Vector{Float64} = Float64[]

    # --- Harvest fluxes (gN/m2/s) patch-level ---
    hrv_leafn_to_litter_patch                 ::Vector{Float64} = Float64[]
    hrv_frootn_to_litter_patch                ::Vector{Float64} = Float64[]
    hrv_leafn_storage_to_litter_patch         ::Vector{Float64} = Float64[]
    hrv_frootn_storage_to_litter_patch        ::Vector{Float64} = Float64[]
    hrv_livestemn_storage_to_litter_patch     ::Vector{Float64} = Float64[]
    hrv_deadstemn_storage_to_litter_patch     ::Vector{Float64} = Float64[]
    hrv_livecrootn_storage_to_litter_patch    ::Vector{Float64} = Float64[]
    hrv_deadcrootn_storage_to_litter_patch    ::Vector{Float64} = Float64[]
    hrv_leafn_xfer_to_litter_patch            ::Vector{Float64} = Float64[]
    hrv_frootn_xfer_to_litter_patch           ::Vector{Float64} = Float64[]
    hrv_livestemn_xfer_to_litter_patch        ::Vector{Float64} = Float64[]
    hrv_deadstemn_xfer_to_litter_patch        ::Vector{Float64} = Float64[]
    hrv_livecrootn_xfer_to_litter_patch       ::Vector{Float64} = Float64[]
    hrv_deadcrootn_xfer_to_litter_patch       ::Vector{Float64} = Float64[]
    hrv_livestemn_to_litter_patch             ::Vector{Float64} = Float64[]
    hrv_livecrootn_to_litter_patch            ::Vector{Float64} = Float64[]
    hrv_deadcrootn_to_litter_patch            ::Vector{Float64} = Float64[]
    hrv_retransn_to_litter_patch              ::Vector{Float64} = Float64[]
    crop_harvestn_to_cropprodn_patch          ::Vector{Float64} = Float64[]
    crop_harvestn_to_cropprodn_col            ::Vector{Float64} = Float64[]

    # --- Fire N fluxes (gN/m2/s) patch-level ---
    m_leafn_to_fire_patch                     ::Vector{Float64} = Float64[]
    m_leafn_storage_to_fire_patch             ::Vector{Float64} = Float64[]
    m_leafn_xfer_to_fire_patch                ::Vector{Float64} = Float64[]
    m_livestemn_to_fire_patch                 ::Vector{Float64} = Float64[]
    m_livestemn_storage_to_fire_patch         ::Vector{Float64} = Float64[]
    m_livestemn_xfer_to_fire_patch            ::Vector{Float64} = Float64[]
    m_deadstemn_to_fire_patch                 ::Vector{Float64} = Float64[]
    m_deadstemn_storage_to_fire_patch         ::Vector{Float64} = Float64[]
    m_deadstemn_xfer_to_fire_patch            ::Vector{Float64} = Float64[]
    m_frootn_to_fire_patch                    ::Vector{Float64} = Float64[]
    m_frootn_storage_to_fire_patch            ::Vector{Float64} = Float64[]
    m_frootn_xfer_to_fire_patch               ::Vector{Float64} = Float64[]
    m_livecrootn_to_fire_patch                ::Vector{Float64} = Float64[]
    m_livecrootn_storage_to_fire_patch        ::Vector{Float64} = Float64[]
    m_livecrootn_xfer_to_fire_patch           ::Vector{Float64} = Float64[]
    m_deadcrootn_to_fire_patch                ::Vector{Float64} = Float64[]
    m_deadcrootn_storage_to_fire_patch        ::Vector{Float64} = Float64[]
    m_deadcrootn_xfer_to_fire_patch           ::Vector{Float64} = Float64[]
    m_retransn_to_fire_patch                  ::Vector{Float64} = Float64[]

    # --- Fire-to-litter fluxes (gN/m2/s) patch-level ---
    m_leafn_to_litter_fire_patch              ::Vector{Float64} = Float64[]
    m_leafn_storage_to_litter_fire_patch      ::Vector{Float64} = Float64[]
    m_leafn_xfer_to_litter_fire_patch         ::Vector{Float64} = Float64[]
    m_livestemn_to_litter_fire_patch          ::Vector{Float64} = Float64[]
    m_livestemn_storage_to_litter_fire_patch  ::Vector{Float64} = Float64[]
    m_livestemn_xfer_to_litter_fire_patch     ::Vector{Float64} = Float64[]
    m_livestemn_to_deadstemn_fire_patch       ::Vector{Float64} = Float64[]
    m_deadstemn_to_litter_fire_patch          ::Vector{Float64} = Float64[]
    m_deadstemn_storage_to_litter_fire_patch  ::Vector{Float64} = Float64[]
    m_deadstemn_xfer_to_litter_fire_patch     ::Vector{Float64} = Float64[]
    m_frootn_to_litter_fire_patch             ::Vector{Float64} = Float64[]
    m_frootn_storage_to_litter_fire_patch     ::Vector{Float64} = Float64[]
    m_frootn_xfer_to_litter_fire_patch        ::Vector{Float64} = Float64[]
    m_livecrootn_to_litter_fire_patch         ::Vector{Float64} = Float64[]
    m_livecrootn_storage_to_litter_fire_patch ::Vector{Float64} = Float64[]
    m_livecrootn_xfer_to_litter_fire_patch    ::Vector{Float64} = Float64[]
    m_livecrootn_to_deadcrootn_fire_patch     ::Vector{Float64} = Float64[]
    m_deadcrootn_to_litter_fire_patch         ::Vector{Float64} = Float64[]
    m_deadcrootn_storage_to_litter_fire_patch ::Vector{Float64} = Float64[]
    m_deadcrootn_xfer_to_litter_fire_patch    ::Vector{Float64} = Float64[]
    m_retransn_to_litter_fire_patch           ::Vector{Float64} = Float64[]

    # --- Fire summary (gN/m2/s) ---
    fire_nloss_patch                          ::Vector{Float64} = Float64[]
    fire_nloss_col                            ::Vector{Float64} = Float64[]
    fire_nloss_p2c_col                        ::Vector{Float64} = Float64[]
    fire_mortality_n_to_cwdn_col              ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # --- Column-level decomposition fire fluxes ---
    m_decomp_npools_to_fire_vr_col            ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    m_decomp_npools_to_fire_col               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    m_n_to_litr_fire_col                      ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)

    # --- Phenology fluxes from transfer pool (gN/m2/s) ---
    reproductiven_xfer_to_reproductiven_patch ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    leafn_xfer_to_leafn_patch                 ::Vector{Float64} = Float64[]
    frootn_xfer_to_frootn_patch               ::Vector{Float64} = Float64[]
    livestemn_xfer_to_livestemn_patch         ::Vector{Float64} = Float64[]
    deadstemn_xfer_to_deadstemn_patch         ::Vector{Float64} = Float64[]
    livecrootn_xfer_to_livecrootn_patch       ::Vector{Float64} = Float64[]
    deadcrootn_xfer_to_deadcrootn_patch       ::Vector{Float64} = Float64[]

    # --- Litterfall fluxes (gN/m2/s) ---
    livestemn_to_litter_patch                 ::Vector{Float64} = Float64[]
    repr_grainn_to_food_patch                 ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    repr_grainn_to_food_perharv_patch         ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    repr_grainn_to_food_thisyr_patch          ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    repr_structuren_to_cropprod_patch         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    repr_structuren_to_litter_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    leafn_to_biofueln_patch                   ::Vector{Float64} = Float64[]
    livestemn_to_biofueln_patch               ::Vector{Float64} = Float64[]
    leafn_to_removedresiduen_patch            ::Vector{Float64} = Float64[]
    livestemn_to_removedresiduen_patch        ::Vector{Float64} = Float64[]
    repr_grainn_to_seed_patch                 ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    repr_grainn_to_seed_perharv_patch         ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    repr_grainn_to_seed_thisyr_patch          ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    leafn_to_litter_patch                     ::Vector{Float64} = Float64[]
    leafn_to_retransn_patch                   ::Vector{Float64} = Float64[]
    frootn_to_retransn_patch                  ::Vector{Float64} = Float64[]
    frootn_to_litter_patch                    ::Vector{Float64} = Float64[]

    # --- Allocation fluxes (gN/m2/s) ---
    retransn_to_npool_patch                   ::Vector{Float64} = Float64[]
    free_retransn_to_npool_patch              ::Vector{Float64} = Float64[]
    sminn_to_npool_patch                      ::Vector{Float64} = Float64[]
    npool_to_reproductiven_patch              ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    npool_to_reproductiven_storage_patch      ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    npool_to_leafn_patch                      ::Vector{Float64} = Float64[]
    npool_to_leafn_storage_patch              ::Vector{Float64} = Float64[]
    npool_to_frootn_patch                     ::Vector{Float64} = Float64[]
    npool_to_frootn_storage_patch             ::Vector{Float64} = Float64[]
    npool_to_livestemn_patch                  ::Vector{Float64} = Float64[]
    npool_to_livestemn_storage_patch          ::Vector{Float64} = Float64[]
    npool_to_deadstemn_patch                  ::Vector{Float64} = Float64[]
    npool_to_deadstemn_storage_patch          ::Vector{Float64} = Float64[]
    npool_to_livecrootn_patch                 ::Vector{Float64} = Float64[]
    npool_to_livecrootn_storage_patch         ::Vector{Float64} = Float64[]
    npool_to_deadcrootn_patch                 ::Vector{Float64} = Float64[]
    npool_to_deadcrootn_storage_patch         ::Vector{Float64} = Float64[]

    # --- Storage to transfer turnover (gN/m2/s) ---
    reproductiven_storage_to_xfer_patch       ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    leafn_storage_to_xfer_patch               ::Vector{Float64} = Float64[]
    frootn_storage_to_xfer_patch              ::Vector{Float64} = Float64[]
    livestemn_storage_to_xfer_patch           ::Vector{Float64} = Float64[]
    deadstemn_storage_to_xfer_patch           ::Vector{Float64} = Float64[]
    livecrootn_storage_to_xfer_patch          ::Vector{Float64} = Float64[]
    deadcrootn_storage_to_xfer_patch          ::Vector{Float64} = Float64[]
    fert_patch                                ::Vector{Float64} = Float64[]
    fert_counter_patch                        ::Vector{Float64} = Float64[]
    soyfixn_patch                             ::Vector{Float64} = Float64[]

    # --- Livewood to deadwood turnover (gN/m2/s) ---
    livestemn_to_deadstemn_patch              ::Vector{Float64} = Float64[]
    livestemn_to_retransn_patch               ::Vector{Float64} = Float64[]
    livecrootn_to_deadcrootn_patch            ::Vector{Float64} = Float64[]
    livecrootn_to_retransn_patch              ::Vector{Float64} = Float64[]

    # --- Summary/diagnostic ---
    ndeploy_patch                             ::Vector{Float64} = Float64[]
    wood_harvestn_patch                       ::Vector{Float64} = Float64[]
    wood_harvestn_col                         ::Vector{Float64} = Float64[]

    # --- Column-level phenology/mortality/harvest fluxes (gN/m3/s) ---
    phenology_n_to_litr_n_col                 ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    gap_mortality_n_to_litr_n_col             ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    gap_mortality_n_to_cwdn_col               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    harvest_n_to_litr_n_col                   ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    harvest_n_to_cwdn_col                     ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # --- Dynamic landcover fluxes ---
    dwt_seedn_to_leaf_patch                   ::Vector{Float64} = Float64[]
    dwt_seedn_to_leaf_grc                     ::Vector{Float64} = Float64[]
    dwt_seedn_to_deadstem_patch               ::Vector{Float64} = Float64[]
    dwt_seedn_to_deadstem_grc                 ::Vector{Float64} = Float64[]
    dwt_conv_nflux_patch                      ::Vector{Float64} = Float64[]
    dwt_conv_nflux_grc                        ::Vector{Float64} = Float64[]
    dwt_wood_productn_gain_patch              ::Vector{Float64} = Float64[]
    dwt_crop_productn_gain_patch              ::Vector{Float64} = Float64[]
    dwt_frootn_to_litr_n_col                  ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    dwt_livecrootn_to_cwdn_col                ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    dwt_deadcrootn_to_cwdn_col                ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # --- Gross unrepresented landcover fluxes ---
    gru_leafn_to_litter_patch                 ::Vector{Float64} = Float64[]
    gru_leafn_storage_to_atm_patch            ::Vector{Float64} = Float64[]
    gru_leafn_xfer_to_atm_patch               ::Vector{Float64} = Float64[]
    gru_frootn_to_litter_patch                ::Vector{Float64} = Float64[]
    gru_frootn_storage_to_atm_patch           ::Vector{Float64} = Float64[]
    gru_frootn_xfer_to_atm_patch              ::Vector{Float64} = Float64[]
    gru_livestemn_to_atm_patch                ::Vector{Float64} = Float64[]
    gru_livestemn_storage_to_atm_patch        ::Vector{Float64} = Float64[]
    gru_livestemn_xfer_to_atm_patch           ::Vector{Float64} = Float64[]
    gru_deadstemn_to_atm_patch                ::Vector{Float64} = Float64[]
    gru_deadstemn_storage_to_atm_patch        ::Vector{Float64} = Float64[]
    gru_deadstemn_xfer_to_atm_patch           ::Vector{Float64} = Float64[]
    gru_livecrootn_to_litter_patch            ::Vector{Float64} = Float64[]
    gru_livecrootn_storage_to_atm_patch       ::Vector{Float64} = Float64[]
    gru_livecrootn_xfer_to_atm_patch          ::Vector{Float64} = Float64[]
    gru_deadcrootn_to_litter_patch            ::Vector{Float64} = Float64[]
    gru_deadcrootn_storage_to_atm_patch       ::Vector{Float64} = Float64[]
    gru_deadcrootn_xfer_to_atm_patch          ::Vector{Float64} = Float64[]
    gru_retransn_to_litter_patch              ::Vector{Float64} = Float64[]
    gru_conv_nflux_patch                      ::Vector{Float64} = Float64[]
    gru_conv_nflux_col                        ::Vector{Float64} = Float64[]
    gru_conv_nflux_grc                        ::Vector{Float64} = Float64[]
    gru_wood_productn_gain_patch              ::Vector{Float64} = Float64[]
    gru_wood_productn_gain_col                ::Vector{Float64} = Float64[]
    gru_wood_productn_gain_grc                ::Vector{Float64} = Float64[]
    gru_n_to_litr_n_col                       ::Array{Float64,3} = Array{Float64}(undef, 0, 0, 0)
    gru_n_to_cwdn_col                         ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

    # --- Crop fluxes ---
    crop_seedn_to_leaf_patch                  ::Vector{Float64} = Float64[]

    # --- Misc / FUN fluxes (gN/m2/s) ---
    plant_ndemand_patch                       ::Vector{Float64} = Float64[]
    avail_retransn_patch                      ::Vector{Float64} = Float64[]
    plant_nalloc_patch                        ::Vector{Float64} = Float64[]
    plant_ndemand_retrans_patch               ::Vector{Float64} = Float64[]
    plant_ndemand_season_patch                ::Vector{Float64} = Float64[]
    plant_ndemand_stress_patch                ::Vector{Float64} = Float64[]
    Nactive_patch                             ::Vector{Float64} = Float64[]
    Nnonmyc_patch                             ::Vector{Float64} = Float64[]
    Nam_patch                                 ::Vector{Float64} = Float64[]
    Necm_patch                                ::Vector{Float64} = Float64[]
    Nactive_no3_patch                         ::Vector{Float64} = Float64[]
    Nactive_nh4_patch                         ::Vector{Float64} = Float64[]
    Nnonmyc_no3_patch                         ::Vector{Float64} = Float64[]
    Nnonmyc_nh4_patch                         ::Vector{Float64} = Float64[]
    Nam_no3_patch                             ::Vector{Float64} = Float64[]
    Nam_nh4_patch                             ::Vector{Float64} = Float64[]
    Necm_no3_patch                            ::Vector{Float64} = Float64[]
    Necm_nh4_patch                            ::Vector{Float64} = Float64[]
    Nfix_patch                                ::Vector{Float64} = Float64[]
    Npassive_patch                            ::Vector{Float64} = Float64[]
    Nretrans_patch                            ::Vector{Float64} = Float64[]
    Nretrans_org_patch                        ::Vector{Float64} = Float64[]
    Nretrans_season_patch                     ::Vector{Float64} = Float64[]
    Nretrans_stress_patch                     ::Vector{Float64} = Float64[]
    Nuptake_patch                             ::Vector{Float64} = Float64[]
    sminn_to_plant_fun_patch                  ::Vector{Float64} = Float64[]
    sminn_to_plant_fun_vr_patch               ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    sminn_to_plant_fun_no3_vr_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    sminn_to_plant_fun_nh4_vr_patch           ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    cost_nfix_patch                           ::Vector{Float64} = Float64[]
    cost_nactive_patch                        ::Vector{Float64} = Float64[]
    cost_nretrans_patch                       ::Vector{Float64} = Float64[]
    nuptake_npp_fraction_patch                ::Vector{Float64} = Float64[]

    # --- Matrix CN arrays ---
    matrix_nalloc_patch                       ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    matrix_Ninput_patch                       ::Vector{Float64} = Float64[]
    matrix_nphtransfer_patch                  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    matrix_nphturnover_patch                  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    matrix_nphtransfer_doner_patch            ::Vector{Int} = Int[]
    matrix_nphtransfer_receiver_patch         ::Vector{Int} = Int[]
    matrix_ngmtransfer_patch                  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    matrix_ngmturnover_patch                  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    matrix_ngmtransfer_doner_patch            ::Vector{Int} = Int[]
    matrix_ngmtransfer_receiver_patch         ::Vector{Int} = Int[]
    matrix_nfitransfer_patch                  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    matrix_nfiturnover_patch                  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    matrix_nfitransfer_doner_patch            ::Vector{Int} = Int[]
    matrix_nfitransfer_receiver_patch         ::Vector{Int} = Int[]

    # --- Matrix index scalars (phenology) ---
    ileaf_to_iretransn_ph::Int = 0
    ileafst_to_ileafxf_ph::Int = 0
    ileafxf_to_ileaf_ph::Int = 0
    ifroot_to_iretransn_ph::Int = 0
    ifrootst_to_ifrootxf_ph::Int = 0
    ifrootxf_to_ifroot_ph::Int = 0
    ilivestem_to_ideadstem_ph::Int = 0
    ilivestem_to_iretransn_ph::Int = 0
    ilivestemst_to_ilivestemxf_ph::Int = 0
    ilivestemxf_to_ilivestem_ph::Int = 0
    ideadstemst_to_ideadstemxf_ph::Int = 0
    ideadstemxf_to_ideadstem_ph::Int = 0
    ilivecroot_to_ideadcroot_ph::Int = 0
    ilivecroot_to_iretransn_ph::Int = 0
    ilivecrootst_to_ilivecrootxf_ph::Int = 0
    ilivecrootxf_to_ilivecroot_ph::Int = 0
    ideadcrootst_to_ideadcrootxf_ph::Int = 0
    ideadcrootxf_to_ideadcroot_ph::Int = 0
    iretransn_to_ileaf_ph::Int = 0
    iretransn_to_ileafst_ph::Int = 0
    iretransn_to_ifroot_ph::Int = 0
    iretransn_to_ifrootst_ph::Int = 0
    iretransn_to_ilivestem_ph::Int = 0
    iretransn_to_ilivestemst_ph::Int = 0
    iretransn_to_ideadstem_ph::Int = 0
    iretransn_to_ideadstemst_ph::Int = 0
    iretransn_to_ilivecroot_ph::Int = 0
    iretransn_to_ilivecrootst_ph::Int = 0
    iretransn_to_ideadcroot_ph::Int = 0
    iretransn_to_ideadcrootst_ph::Int = 0
    iretransn_to_igrain_ph::Int = 0
    iretransn_to_igrainst_ph::Int = 0
    ileaf_to_iout_ph::Int = 0
    ifroot_to_iout_ph::Int = 0
    ilivestem_to_iout_ph::Int = 0
    ileaf_to_iretransn_ph_idx::Int = 0  # separate from ileaf_to_iretransn_ph above
    ifroot_to_iretransn_ph_idx::Int = 0
    ilivestem_to_iretransn_ph_idx::Int = 0
    ilivecroot_to_iretransn_ph_idx::Int = 0
    igrain_to_iout_ph::Int = 0
    iretransn_to_iout_ph::Int = 0

    # --- Matrix index scalars (gap mortality) ---
    ileaf_to_iout_gm::Int = 0
    ileafst_to_iout_gm::Int = 0
    ileafxf_to_iout_gm::Int = 0
    ifroot_to_iout_gm::Int = 0
    ifrootst_to_iout_gm::Int = 0
    ifrootxf_to_iout_gm::Int = 0
    ilivestem_to_iout_gm::Int = 0
    ilivestemst_to_iout_gm::Int = 0
    ilivestemxf_to_iout_gm::Int = 0
    ideadstem_to_iout_gm::Int = 0
    ideadstemst_to_iout_gm::Int = 0
    ideadstemxf_to_iout_gm::Int = 0
    ilivecroot_to_iout_gm::Int = 0
    ilivecrootst_to_iout_gm::Int = 0
    ilivecrootxf_to_iout_gm::Int = 0
    ideadcroot_to_iout_gm::Int = 0
    ideadcrootst_to_iout_gm::Int = 0
    ideadcrootxf_to_iout_gm::Int = 0
    iretransn_to_iout_gm::Int = 0

    # --- Matrix index scalars (fire) ---
    ilivestem_to_ideadstem_fi::Int = 0
    ilivecroot_to_ideadcroot_fi::Int = 0
    ileaf_to_iout_fi::Int = 0
    ileafst_to_iout_fi::Int = 0
    ileafxf_to_iout_fi::Int = 0
    ifroot_to_iout_fi::Int = 0
    ifrootst_to_iout_fi::Int = 0
    ifrootxf_to_iout_fi::Int = 0
    ilivestem_to_iout_fi::Int = 0
    ilivestemst_to_iout_fi::Int = 0
    ilivestemxf_to_iout_fi::Int = 0
    ideadstem_to_iout_fi::Int = 0
    ideadstemst_to_iout_fi::Int = 0
    ideadstemxf_to_iout_fi::Int = 0
    ilivecroot_to_iout_fi::Int = 0
    ilivecrootst_to_iout_fi::Int = 0
    ilivecrootxf_to_iout_fi::Int = 0
    ideadcroot_to_iout_fi::Int = 0
    ideadcrootst_to_iout_fi::Int = 0
    ideadcrootxf_to_iout_fi::Int = 0
    iretransn_to_iout_fi::Int = 0
end

# ==========================================================================
# Init function — allocates all arrays
# Ported from: InitAllocate subroutine
# ==========================================================================

"""
    cnveg_nitrogen_flux_init!(nf, np, nc, ng; ...)

Allocate and initialize all nitrogen flux arrays.
"""
function cnveg_nitrogen_flux_init!(nf::CNVegNitrogenFluxData, np::Int, nc::Int, ng::Int;
                                   use_matrixcn::Bool=false,
                                   nrepr::Int=NREPR,
                                   nlevdecomp_full::Int=1,
                                   ndecomp_pools::Int=1,
                                   mxharvests::Int=MXHARVESTS,
                                   nvegnpool::Int=NVEGPOOL_NATVEG,
                                   i_litr_max::Int=I_LITR1)
    nanvec(n) = fill(NaN, n)
    nanmat(r, c) = fill(NaN, r, c)
    nan3d(a, b, c) = fill(NaN, a, b, c)

    # --- Gap mortality ---
    nf.m_leafn_to_litter_patch              = nanvec(np)
    nf.m_frootn_to_litter_patch             = nanvec(np)
    nf.m_leafn_storage_to_litter_patch      = nanvec(np)
    nf.m_frootn_storage_to_litter_patch     = nanvec(np)
    nf.m_livestemn_storage_to_litter_patch  = nanvec(np)
    nf.m_deadstemn_storage_to_litter_patch  = nanvec(np)
    nf.m_livecrootn_storage_to_litter_patch = nanvec(np)
    nf.m_deadcrootn_storage_to_litter_patch = nanvec(np)
    nf.m_leafn_xfer_to_litter_patch         = nanvec(np)
    nf.m_frootn_xfer_to_litter_patch        = nanvec(np)
    nf.m_livestemn_xfer_to_litter_patch     = nanvec(np)
    nf.m_deadstemn_xfer_to_litter_patch     = nanvec(np)
    nf.m_livecrootn_xfer_to_litter_patch    = nanvec(np)
    nf.m_deadcrootn_xfer_to_litter_patch    = nanvec(np)
    nf.m_livestemn_to_litter_patch          = nanvec(np)
    nf.m_deadstemn_to_litter_patch          = nanvec(np)
    nf.m_livecrootn_to_litter_patch         = nanvec(np)
    nf.m_deadcrootn_to_litter_patch         = nanvec(np)
    nf.m_retransn_to_litter_patch           = nanvec(np)

    # --- Harvest ---
    nf.hrv_leafn_to_litter_patch            = nanvec(np)
    nf.hrv_frootn_to_litter_patch           = nanvec(np)
    nf.hrv_leafn_storage_to_litter_patch    = nanvec(np)
    nf.hrv_frootn_storage_to_litter_patch   = nanvec(np)
    nf.hrv_livestemn_storage_to_litter_patch = nanvec(np)
    nf.hrv_deadstemn_storage_to_litter_patch = nanvec(np)
    nf.hrv_livecrootn_storage_to_litter_patch = nanvec(np)
    nf.hrv_deadcrootn_storage_to_litter_patch = nanvec(np)
    nf.hrv_leafn_xfer_to_litter_patch       = nanvec(np)
    nf.hrv_frootn_xfer_to_litter_patch      = nanvec(np)
    nf.hrv_livestemn_xfer_to_litter_patch   = nanvec(np)
    nf.hrv_deadstemn_xfer_to_litter_patch   = nanvec(np)
    nf.hrv_livecrootn_xfer_to_litter_patch  = nanvec(np)
    nf.hrv_deadcrootn_xfer_to_litter_patch  = nanvec(np)
    nf.hrv_livestemn_to_litter_patch        = nanvec(np)
    nf.hrv_livecrootn_to_litter_patch       = nanvec(np)
    nf.hrv_deadcrootn_to_litter_patch       = nanvec(np)
    nf.hrv_retransn_to_litter_patch         = nanvec(np)

    # --- Fire ---
    nf.m_leafn_to_fire_patch                = nanvec(np)
    nf.m_leafn_storage_to_fire_patch        = nanvec(np)
    nf.m_leafn_xfer_to_fire_patch           = nanvec(np)
    nf.m_livestemn_to_fire_patch            = nanvec(np)
    nf.m_livestemn_storage_to_fire_patch    = nanvec(np)
    nf.m_livestemn_xfer_to_fire_patch       = nanvec(np)
    nf.m_deadstemn_to_fire_patch            = nanvec(np)
    nf.m_deadstemn_storage_to_fire_patch    = nanvec(np)
    nf.m_deadstemn_xfer_to_fire_patch       = nanvec(np)
    nf.m_frootn_to_fire_patch               = nanvec(np)
    nf.m_frootn_storage_to_fire_patch       = nanvec(np)
    nf.m_frootn_xfer_to_fire_patch          = nanvec(np)
    nf.m_livecrootn_to_fire_patch           = nanvec(np)
    nf.m_livecrootn_storage_to_fire_patch   = nanvec(np)
    nf.m_livecrootn_xfer_to_fire_patch      = nanvec(np)
    nf.m_deadcrootn_to_fire_patch           = nanvec(np)
    nf.m_deadcrootn_storage_to_fire_patch   = nanvec(np)
    nf.m_deadcrootn_xfer_to_fire_patch      = nanvec(np)
    nf.m_retransn_to_fire_patch             = nanvec(np)

    # --- Fire-to-litter ---
    nf.m_leafn_to_litter_fire_patch              = nanvec(np)
    nf.m_leafn_storage_to_litter_fire_patch      = nanvec(np)
    nf.m_leafn_xfer_to_litter_fire_patch         = nanvec(np)
    nf.m_livestemn_to_litter_fire_patch          = nanvec(np)
    nf.m_livestemn_storage_to_litter_fire_patch  = nanvec(np)
    nf.m_livestemn_xfer_to_litter_fire_patch     = nanvec(np)
    nf.m_livestemn_to_deadstemn_fire_patch       = nanvec(np)
    nf.m_deadstemn_to_litter_fire_patch          = nanvec(np)
    nf.m_deadstemn_storage_to_litter_fire_patch  = nanvec(np)
    nf.m_deadstemn_xfer_to_litter_fire_patch     = nanvec(np)
    nf.m_frootn_to_litter_fire_patch             = nanvec(np)
    nf.m_frootn_storage_to_litter_fire_patch     = nanvec(np)
    nf.m_frootn_xfer_to_litter_fire_patch        = nanvec(np)
    nf.m_livecrootn_to_litter_fire_patch         = nanvec(np)
    nf.m_livecrootn_storage_to_litter_fire_patch = nanvec(np)
    nf.m_livecrootn_xfer_to_litter_fire_patch    = nanvec(np)
    nf.m_livecrootn_to_deadcrootn_fire_patch     = nanvec(np)
    nf.m_deadcrootn_to_litter_fire_patch         = nanvec(np)
    nf.m_deadcrootn_storage_to_litter_fire_patch = nanvec(np)
    nf.m_deadcrootn_xfer_to_litter_fire_patch    = nanvec(np)
    nf.m_retransn_to_litter_fire_patch           = nanvec(np)

    # --- Phenology ---
    nf.reproductiven_xfer_to_reproductiven_patch = nanmat(np, nrepr)
    nf.leafn_xfer_to_leafn_patch            = nanvec(np)
    nf.frootn_xfer_to_frootn_patch          = nanvec(np)
    nf.livestemn_xfer_to_livestemn_patch    = nanvec(np)
    nf.deadstemn_xfer_to_deadstemn_patch    = nanvec(np)
    nf.livecrootn_xfer_to_livecrootn_patch  = nanvec(np)
    nf.deadcrootn_xfer_to_deadcrootn_patch  = nanvec(np)

    # --- Litterfall / crop ---
    nf.leafn_to_litter_patch                = nanvec(np)
    nf.leafn_to_retransn_patch              = nanvec(np)
    nf.frootn_to_retransn_patch             = nanvec(np)
    nf.frootn_to_litter_patch               = nanvec(np)
    nf.retransn_to_npool_patch              = nanvec(np)
    nf.free_retransn_to_npool_patch         = nanvec(np)
    nf.sminn_to_npool_patch                 = nanvec(np)
    nf.npool_to_leafn_patch                 = nanvec(np)
    nf.npool_to_leafn_storage_patch         = nanvec(np)
    nf.npool_to_frootn_patch                = nanvec(np)
    nf.npool_to_frootn_storage_patch        = nanvec(np)
    nf.npool_to_livestemn_patch             = nanvec(np)
    nf.npool_to_livestemn_storage_patch     = nanvec(np)
    nf.npool_to_deadstemn_patch             = nanvec(np)
    nf.npool_to_deadstemn_storage_patch     = nanvec(np)
    nf.npool_to_livecrootn_patch            = nanvec(np)
    nf.npool_to_livecrootn_storage_patch    = nanvec(np)
    nf.npool_to_deadcrootn_patch            = nanvec(np)
    nf.npool_to_deadcrootn_storage_patch    = nanvec(np)
    nf.npool_to_reproductiven_patch         = nanmat(np, nrepr)
    nf.npool_to_reproductiven_storage_patch = nanmat(np, nrepr)
    nf.reproductiven_storage_to_xfer_patch  = nanmat(np, nrepr)
    nf.leafn_storage_to_xfer_patch          = nanvec(np)
    nf.frootn_storage_to_xfer_patch         = nanvec(np)
    nf.livestemn_storage_to_xfer_patch      = nanvec(np)
    nf.deadstemn_storage_to_xfer_patch      = nanvec(np)
    nf.livecrootn_storage_to_xfer_patch     = nanvec(np)
    nf.deadcrootn_storage_to_xfer_patch     = nanvec(np)
    nf.livestemn_to_deadstemn_patch         = nanvec(np)
    nf.livestemn_to_retransn_patch          = nanvec(np)
    nf.livecrootn_to_deadcrootn_patch       = nanvec(np)
    nf.livecrootn_to_retransn_patch         = nanvec(np)
    nf.ndeploy_patch                        = nanvec(np)
    nf.wood_harvestn_patch                  = nanvec(np)
    nf.fire_nloss_patch                     = nanvec(np)
    nf.livestemn_to_litter_patch            = nanvec(np)
    nf.repr_grainn_to_food_patch            = nanmat(np, nrepr)
    nf.repr_grainn_to_food_perharv_patch    = nan3d(np, mxharvests, nrepr)
    nf.repr_grainn_to_food_thisyr_patch     = nanmat(np, nrepr)
    nf.repr_structuren_to_cropprod_patch    = nanmat(np, nrepr)
    nf.repr_structuren_to_litter_patch      = nanmat(np, nrepr)
    nf.leafn_to_biofueln_patch              = nanvec(np)
    nf.livestemn_to_biofueln_patch          = nanvec(np)
    nf.leafn_to_removedresiduen_patch       = nanvec(np)
    nf.livestemn_to_removedresiduen_patch   = nanvec(np)
    nf.repr_grainn_to_seed_patch            = nanmat(np, nrepr)
    nf.repr_grainn_to_seed_perharv_patch    = nan3d(np, mxharvests, nrepr)
    nf.repr_grainn_to_seed_thisyr_patch     = nanmat(np, nrepr)
    nf.fert_patch                           = nanvec(np)
    nf.fert_counter_patch                   = nanvec(np)
    nf.soyfixn_patch                        = nanvec(np)

    nf.crop_harvestn_to_cropprodn_patch     = nanvec(np)
    nf.crop_harvestn_to_cropprodn_col       = nanvec(nc)

    # --- Fire summary ---
    nf.fire_nloss_col                       = nanvec(nc)
    nf.fire_nloss_p2c_col                   = nanvec(nc)

    # --- Column-level decomp fire ---
    nf.m_n_to_litr_fire_col                 = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    nf.m_decomp_npools_to_fire_vr_col       = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    nf.m_decomp_npools_to_fire_col          = nanmat(nc, ndecomp_pools)

    # --- Column-level phenology/mortality/harvest ---
    nf.phenology_n_to_litr_n_col            = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    nf.gap_mortality_n_to_litr_n_col        = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    nf.gap_mortality_n_to_cwdn_col          = nanmat(nc, nlevdecomp_full)
    nf.fire_mortality_n_to_cwdn_col         = nanmat(nc, nlevdecomp_full)
    nf.harvest_n_to_litr_n_col              = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    nf.harvest_n_to_cwdn_col                = nanmat(nc, nlevdecomp_full)

    # --- Dynamic landcover ---
    nf.dwt_seedn_to_leaf_patch              = nanvec(np)
    nf.dwt_seedn_to_leaf_grc                = nanvec(ng)
    nf.dwt_seedn_to_deadstem_patch          = nanvec(np)
    nf.dwt_seedn_to_deadstem_grc            = nanvec(ng)
    nf.dwt_conv_nflux_patch                 = nanvec(np)
    nf.dwt_conv_nflux_grc                   = nanvec(ng)
    nf.dwt_wood_productn_gain_patch         = nanvec(np)
    nf.dwt_crop_productn_gain_patch         = nanvec(np)
    nf.wood_harvestn_col                    = nanvec(nc)
    nf.dwt_frootn_to_litr_n_col             = nan3d(nc, nlevdecomp_full, i_litr_max)
    nf.dwt_livecrootn_to_cwdn_col           = nanmat(nc, nlevdecomp_full)
    nf.dwt_deadcrootn_to_cwdn_col           = nanmat(nc, nlevdecomp_full)

    # --- GRU ---
    nf.gru_leafn_to_litter_patch            = nanvec(np)
    nf.gru_leafn_storage_to_atm_patch       = nanvec(np)
    nf.gru_leafn_xfer_to_atm_patch          = nanvec(np)
    nf.gru_frootn_to_litter_patch           = nanvec(np)
    nf.gru_frootn_storage_to_atm_patch      = nanvec(np)
    nf.gru_frootn_xfer_to_atm_patch         = nanvec(np)
    nf.gru_livestemn_to_atm_patch           = nanvec(np)
    nf.gru_livestemn_storage_to_atm_patch   = nanvec(np)
    nf.gru_livestemn_xfer_to_atm_patch      = nanvec(np)
    nf.gru_deadstemn_to_atm_patch           = nanvec(np)
    nf.gru_deadstemn_storage_to_atm_patch   = nanvec(np)
    nf.gru_deadstemn_xfer_to_atm_patch      = nanvec(np)
    nf.gru_livecrootn_storage_to_atm_patch  = nanvec(np)
    nf.gru_livecrootn_xfer_to_atm_patch     = nanvec(np)
    nf.gru_livecrootn_to_litter_patch       = nanvec(np)
    nf.gru_deadcrootn_storage_to_atm_patch  = nanvec(np)
    nf.gru_deadcrootn_xfer_to_atm_patch     = nanvec(np)
    nf.gru_deadcrootn_to_litter_patch       = nanvec(np)
    nf.gru_retransn_to_litter_patch         = nanvec(np)
    nf.gru_conv_nflux_patch                 = nanvec(np)
    nf.gru_conv_nflux_col                   = nanvec(nc)
    nf.gru_conv_nflux_grc                   = nanvec(ng)
    nf.gru_wood_productn_gain_patch         = nanvec(np)
    nf.gru_wood_productn_gain_col           = nanvec(nc)
    nf.gru_wood_productn_gain_grc           = nanvec(ng)
    nf.gru_n_to_litr_n_col                  = nan3d(nc, nlevdecomp_full, i_litr_max)
    nf.gru_n_to_cwdn_col                    = nanmat(nc, nlevdecomp_full)

    # --- Crop ---
    nf.crop_seedn_to_leaf_patch             = nanvec(np)

    # --- Misc / FUN ---
    nf.plant_ndemand_patch                  = nanvec(np)
    nf.avail_retransn_patch                 = nanvec(np)
    nf.plant_nalloc_patch                   = nanvec(np)
    nf.plant_ndemand_retrans_patch          = nanvec(np)
    nf.plant_ndemand_season_patch           = nanvec(np)
    nf.plant_ndemand_stress_patch           = nanvec(np)
    nf.Nactive_patch                        = nanvec(np)
    nf.Nnonmyc_patch                        = nanvec(np)
    nf.Nam_patch                            = nanvec(np)
    nf.Necm_patch                           = nanvec(np)
    nf.Nactive_no3_patch                    = nanvec(np)
    nf.Nactive_nh4_patch                    = nanvec(np)
    nf.Nnonmyc_no3_patch                    = nanvec(np)
    nf.Nnonmyc_nh4_patch                    = nanvec(np)
    nf.Nam_no3_patch                        = nanvec(np)
    nf.Nam_nh4_patch                        = nanvec(np)
    nf.Necm_no3_patch                       = nanvec(np)
    nf.Necm_nh4_patch                       = nanvec(np)
    nf.Npassive_patch                       = nanvec(np)
    nf.Nfix_patch                           = nanvec(np)
    nf.Nretrans_patch                       = nanvec(np)
    nf.Nretrans_org_patch                   = nanvec(np)
    nf.Nretrans_season_patch                = nanvec(np)
    nf.Nretrans_stress_patch                = nanvec(np)
    nf.Nuptake_patch                        = nanvec(np)
    nf.sminn_to_plant_fun_patch             = nanvec(np)
    nf.sminn_to_plant_fun_vr_patch          = nanmat(np, nlevdecomp_full)
    nf.sminn_to_plant_fun_no3_vr_patch      = nanmat(np, nlevdecomp_full)
    nf.sminn_to_plant_fun_nh4_vr_patch      = nanmat(np, nlevdecomp_full)
    nf.cost_nfix_patch                      = nanvec(np)
    nf.cost_nactive_patch                   = nanvec(np)
    nf.cost_nretrans_patch                  = nanvec(np)
    nf.nuptake_npp_fraction_patch           = nanvec(np)

    # --- Matrix CN (optional) ---
    if use_matrixcn
        nf.matrix_Ninput_patch              = nanvec(np)
        nf.matrix_nalloc_patch              = nanmat(np, nvegnpool)
        nf.matrix_nphturnover_patch         = nanmat(np, nvegnpool)
        nf.matrix_ngmturnover_patch         = nanmat(np, nvegnpool)
        nf.matrix_nfiturnover_patch         = nanmat(np, nvegnpool)
    end

    nothing
end

# ==========================================================================
# SetValues — bulk-set all flux variables for given masks
# Ported from: SetValues subroutine
# ==========================================================================

"""
    cnveg_nitrogen_flux_set_values!(nf, mask_patch, value_patch, mask_col, value_column; ...)

Set all nitrogen flux variables to given values for masked patches/columns.
"""
function cnveg_nitrogen_flux_set_values!(nf::CNVegNitrogenFluxData,
                                          mask_patch::BitVector, value_patch::Float64,
                                          mask_col::BitVector, value_column::Float64;
                                          use_matrixcn::Bool=false,
                                          use_crop::Bool=false,
                                          nrepr::Int=NREPR,
                                          nlevdecomp_full::Int=1,
                                          ndecomp_pools::Int=1,
                                          i_litr_max::Int=I_LITR1,
                                          nvegnpool::Int=NVEGPOOL_NATVEG)
    # Patch-level 1D fields
    for i in eachindex(mask_patch)
        mask_patch[i] || continue

        nf.m_leafn_to_litter_patch[i]              = value_patch
        nf.m_frootn_to_litter_patch[i]             = value_patch
        nf.m_leafn_storage_to_litter_patch[i]      = value_patch
        nf.m_frootn_storage_to_litter_patch[i]     = value_patch
        nf.m_livestemn_storage_to_litter_patch[i]  = value_patch
        nf.m_deadstemn_storage_to_litter_patch[i]  = value_patch
        nf.m_livecrootn_storage_to_litter_patch[i] = value_patch
        nf.m_deadcrootn_storage_to_litter_patch[i] = value_patch
        nf.m_leafn_xfer_to_litter_patch[i]         = value_patch
        nf.m_frootn_xfer_to_litter_patch[i]        = value_patch
        nf.m_livestemn_xfer_to_litter_patch[i]     = value_patch
        nf.m_deadstemn_xfer_to_litter_patch[i]     = value_patch
        nf.m_livecrootn_xfer_to_litter_patch[i]    = value_patch
        nf.m_deadcrootn_xfer_to_litter_patch[i]    = value_patch
        nf.m_livestemn_to_litter_patch[i]          = value_patch
        nf.m_deadstemn_to_litter_patch[i]          = value_patch
        nf.m_livecrootn_to_litter_patch[i]         = value_patch
        nf.m_deadcrootn_to_litter_patch[i]         = value_patch
        nf.m_retransn_to_litter_patch[i]           = value_patch

        nf.hrv_leafn_to_litter_patch[i]            = value_patch
        nf.hrv_frootn_to_litter_patch[i]           = value_patch
        nf.hrv_leafn_storage_to_litter_patch[i]    = value_patch
        nf.hrv_frootn_storage_to_litter_patch[i]   = value_patch
        nf.hrv_livestemn_storage_to_litter_patch[i] = value_patch
        nf.hrv_deadstemn_storage_to_litter_patch[i] = value_patch
        nf.hrv_livecrootn_storage_to_litter_patch[i] = value_patch
        nf.hrv_deadcrootn_storage_to_litter_patch[i] = value_patch
        nf.hrv_leafn_xfer_to_litter_patch[i]       = value_patch
        nf.hrv_frootn_xfer_to_litter_patch[i]      = value_patch
        nf.hrv_livestemn_xfer_to_litter_patch[i]   = value_patch
        nf.hrv_deadstemn_xfer_to_litter_patch[i]   = value_patch
        nf.hrv_livecrootn_xfer_to_litter_patch[i]  = value_patch
        nf.hrv_deadcrootn_xfer_to_litter_patch[i]  = value_patch
        nf.hrv_livestemn_to_litter_patch[i]        = value_patch
        nf.hrv_livecrootn_to_litter_patch[i]       = value_patch
        nf.hrv_deadcrootn_to_litter_patch[i]       = value_patch
        nf.hrv_retransn_to_litter_patch[i]         = value_patch

        nf.gru_leafn_to_litter_patch[i]            = value_patch
        nf.gru_leafn_storage_to_atm_patch[i]       = value_patch
        nf.gru_leafn_xfer_to_atm_patch[i]          = value_patch
        nf.gru_frootn_to_litter_patch[i]           = value_patch
        nf.gru_frootn_storage_to_atm_patch[i]      = value_patch
        nf.gru_frootn_xfer_to_atm_patch[i]         = value_patch
        nf.gru_livestemn_to_atm_patch[i]           = value_patch
        nf.gru_livestemn_storage_to_atm_patch[i]   = value_patch
        nf.gru_livestemn_xfer_to_atm_patch[i]      = value_patch
        nf.gru_deadstemn_to_atm_patch[i]           = value_patch
        nf.gru_deadstemn_storage_to_atm_patch[i]   = value_patch
        nf.gru_deadstemn_xfer_to_atm_patch[i]      = value_patch
        nf.gru_livecrootn_to_litter_patch[i]       = value_patch
        nf.gru_livecrootn_storage_to_atm_patch[i]  = value_patch
        nf.gru_livecrootn_xfer_to_atm_patch[i]     = value_patch
        nf.gru_deadcrootn_storage_to_atm_patch[i]  = value_patch
        nf.gru_deadcrootn_xfer_to_atm_patch[i]     = value_patch
        nf.gru_deadcrootn_to_litter_patch[i]       = value_patch
        nf.gru_retransn_to_litter_patch[i]         = value_patch
        nf.gru_conv_nflux_patch[i]                 = value_patch
        nf.gru_wood_productn_gain_patch[i]         = value_patch

        nf.m_leafn_to_fire_patch[i]                = value_patch
        nf.m_leafn_storage_to_fire_patch[i]        = value_patch
        nf.m_leafn_xfer_to_fire_patch[i]           = value_patch
        nf.m_livestemn_to_fire_patch[i]            = value_patch
        nf.m_livestemn_storage_to_fire_patch[i]    = value_patch
        nf.m_livestemn_xfer_to_fire_patch[i]       = value_patch
        nf.m_deadstemn_to_fire_patch[i]            = value_patch
        nf.m_deadstemn_storage_to_fire_patch[i]    = value_patch
        nf.m_deadstemn_xfer_to_fire_patch[i]       = value_patch
        nf.m_frootn_to_fire_patch[i]               = value_patch
        nf.m_frootn_storage_to_fire_patch[i]       = value_patch
        nf.m_frootn_xfer_to_fire_patch[i]          = value_patch
        nf.m_livecrootn_to_fire_patch[i]           = value_patch
        nf.m_livecrootn_storage_to_fire_patch[i]   = value_patch
        nf.m_livecrootn_xfer_to_fire_patch[i]      = value_patch
        nf.m_deadcrootn_to_fire_patch[i]           = value_patch
        nf.m_deadcrootn_storage_to_fire_patch[i]   = value_patch
        nf.m_deadcrootn_xfer_to_fire_patch[i]      = value_patch
        nf.m_retransn_to_fire_patch[i]             = value_patch

        nf.m_leafn_to_litter_fire_patch[i]              = value_patch
        nf.m_leafn_storage_to_litter_fire_patch[i]      = value_patch
        nf.m_leafn_xfer_to_litter_fire_patch[i]         = value_patch
        nf.m_livestemn_to_litter_fire_patch[i]          = value_patch
        nf.m_livestemn_storage_to_litter_fire_patch[i]  = value_patch
        nf.m_livestemn_xfer_to_litter_fire_patch[i]     = value_patch
        nf.m_livestemn_to_deadstemn_fire_patch[i]       = value_patch
        nf.m_deadstemn_to_litter_fire_patch[i]          = value_patch
        nf.m_deadstemn_storage_to_litter_fire_patch[i]  = value_patch
        nf.m_deadstemn_xfer_to_litter_fire_patch[i]     = value_patch
        nf.m_frootn_to_litter_fire_patch[i]             = value_patch
        nf.m_frootn_storage_to_litter_fire_patch[i]     = value_patch
        nf.m_frootn_xfer_to_litter_fire_patch[i]        = value_patch
        nf.m_livecrootn_to_litter_fire_patch[i]         = value_patch
        nf.m_livecrootn_storage_to_litter_fire_patch[i] = value_patch
        nf.m_livecrootn_xfer_to_litter_fire_patch[i]    = value_patch
        nf.m_livecrootn_to_deadcrootn_fire_patch[i]     = value_patch
        nf.m_deadcrootn_to_litter_fire_patch[i]         = value_patch
        nf.m_deadcrootn_storage_to_litter_fire_patch[i] = value_patch
        nf.m_deadcrootn_xfer_to_litter_fire_patch[i]    = value_patch
        nf.m_retransn_to_litter_fire_patch[i]           = value_patch

        nf.leafn_xfer_to_leafn_patch[i]            = value_patch
        nf.frootn_xfer_to_frootn_patch[i]          = value_patch
        nf.livestemn_xfer_to_livestemn_patch[i]    = value_patch
        nf.deadstemn_xfer_to_deadstemn_patch[i]    = value_patch
        nf.livecrootn_xfer_to_livecrootn_patch[i]  = value_patch
        nf.deadcrootn_xfer_to_deadcrootn_patch[i]  = value_patch
        nf.leafn_to_litter_patch[i]                = value_patch
        nf.leafn_to_retransn_patch[i]              = value_patch
        nf.frootn_to_litter_patch[i]               = value_patch
        nf.retransn_to_npool_patch[i]              = value_patch
        nf.free_retransn_to_npool_patch[i]         = value_patch
        nf.sminn_to_npool_patch[i]                 = value_patch
        nf.npool_to_leafn_patch[i]                 = value_patch
        nf.npool_to_leafn_storage_patch[i]         = value_patch
        nf.npool_to_frootn_patch[i]                = value_patch
        nf.npool_to_frootn_storage_patch[i]        = value_patch
        nf.npool_to_livestemn_patch[i]             = value_patch
        nf.npool_to_livestemn_storage_patch[i]     = value_patch
        nf.npool_to_deadstemn_patch[i]             = value_patch
        nf.npool_to_deadstemn_storage_patch[i]     = value_patch
        nf.npool_to_livecrootn_patch[i]            = value_patch
        nf.npool_to_livecrootn_storage_patch[i]    = value_patch
        nf.npool_to_deadcrootn_patch[i]            = value_patch
        nf.npool_to_deadcrootn_storage_patch[i]    = value_patch
        nf.leafn_storage_to_xfer_patch[i]          = value_patch
        nf.frootn_storage_to_xfer_patch[i]         = value_patch
        nf.livestemn_storage_to_xfer_patch[i]      = value_patch
        nf.deadstemn_storage_to_xfer_patch[i]      = value_patch
        nf.livecrootn_storage_to_xfer_patch[i]     = value_patch
        nf.deadcrootn_storage_to_xfer_patch[i]     = value_patch
        nf.livestemn_to_deadstemn_patch[i]         = value_patch
        nf.livestemn_to_retransn_patch[i]          = value_patch
        nf.livecrootn_to_deadcrootn_patch[i]       = value_patch
        nf.livecrootn_to_retransn_patch[i]         = value_patch
        nf.ndeploy_patch[i]                        = value_patch
        nf.wood_harvestn_patch[i]                  = value_patch
        nf.fire_nloss_patch[i]                     = value_patch
        nf.crop_seedn_to_leaf_patch[i]             = value_patch
        nf.crop_harvestn_to_cropprodn_patch[i]     = value_patch
    end

    # Crop-specific patch-level
    if use_crop
        for i in eachindex(mask_patch)
            mask_patch[i] || continue
            nf.livestemn_to_litter_patch[i]         = value_patch
            nf.leafn_to_biofueln_patch[i]           = value_patch
            nf.livestemn_to_biofueln_patch[i]       = value_patch
            nf.leafn_to_removedresiduen_patch[i]    = value_patch
            nf.livestemn_to_removedresiduen_patch[i] = value_patch
            nf.soyfixn_patch[i]                     = value_patch
            nf.frootn_to_retransn_patch[i]          = value_patch
        end

        for k in 1:nrepr
            for i in eachindex(mask_patch)
                mask_patch[i] || continue
                nf.reproductiven_xfer_to_reproductiven_patch[i,k] = value_patch
                nf.npool_to_reproductiven_patch[i,k]         = value_patch
                nf.npool_to_reproductiven_storage_patch[i,k] = value_patch
                nf.reproductiven_storage_to_xfer_patch[i,k]  = value_patch
                nf.repr_grainn_to_food_patch[i,k]            = value_patch
                nf.repr_grainn_to_seed_patch[i,k]            = value_patch
                nf.repr_structuren_to_cropprod_patch[i,k]    = value_patch
                nf.repr_structuren_to_litter_patch[i,k]      = value_patch
            end
        end
    end

    # Column-level 3D decomposition fields
    for j in 1:nlevdecomp_full
        for i in eachindex(mask_col)
            mask_col[i] || continue
            for k in 1:i_litr_max
                nf.phenology_n_to_litr_n_col[i,j,k]     = value_column
                nf.gap_mortality_n_to_litr_n_col[i,j,k]  = value_column
                nf.harvest_n_to_litr_n_col[i,j,k]        = value_column
                nf.m_n_to_litr_fire_col[i,j,k]           = value_column
                nf.gru_n_to_litr_n_col[i,j,k]            = value_column
            end
            nf.gap_mortality_n_to_cwdn_col[i,j]          = value_column
            nf.fire_mortality_n_to_cwdn_col[i,j]         = value_column
            nf.harvest_n_to_cwdn_col[i,j]                = value_column
            nf.gru_n_to_cwdn_col[i,j]                    = value_column
        end
    end

    # Column-level 1D
    for i in eachindex(mask_col)
        mask_col[i] || continue
        nf.crop_harvestn_to_cropprodn_col[i]    = value_column
        nf.fire_nloss_col[i]                    = value_column
        nf.fire_nloss_p2c_col[i]                = value_column
        nf.wood_harvestn_col[i]                 = value_column
        nf.gru_conv_nflux_col[i]                = value_column
        nf.gru_wood_productn_gain_col[i]        = value_column
    end

    # Column-level decomp pools
    for k in 1:ndecomp_pools
        for i in eachindex(mask_col)
            mask_col[i] || continue
            nf.m_decomp_npools_to_fire_col[i,k] = value_column
        end
    end

    # Decomp fire vr
    for k in 1:ndecomp_pools
        for j in 1:nlevdecomp_full
            for i in eachindex(mask_col)
                mask_col[i] || continue
                nf.m_decomp_npools_to_fire_vr_col[i,j,k] = value_column
            end
        end
    end

    # Matrix (optional)
    if use_matrixcn
        for j in 1:nvegnpool
            for i in eachindex(mask_patch)
                mask_patch[i] || continue
                nf.matrix_nalloc_patch[i,j]       = value_patch
                nf.matrix_nphturnover_patch[i,j]  = value_patch
                nf.matrix_ngmturnover_patch[i,j]  = value_patch
                nf.matrix_nfiturnover_patch[i,j]  = value_patch
            end
        end
    end

    nothing
end

# ==========================================================================
# ZeroDWT — zero dynamic landuse fluxes
# Ported from: ZeroDwt subroutine
# ==========================================================================

function cnveg_nitrogen_flux_zero_dwt!(nf::CNVegNitrogenFluxData,
                                        bounds_grc::UnitRange{Int},
                                        bounds_col::UnitRange{Int};
                                        nlevdecomp_full::Int=1,
                                        i_litr_max::Int=I_LITR1)
    for g in bounds_grc
        nf.dwt_seedn_to_leaf_grc[g]     = 0.0
        nf.dwt_seedn_to_deadstem_grc[g] = 0.0
        nf.dwt_conv_nflux_grc[g]        = 0.0
    end

    for j in 1:nlevdecomp_full
        for c in bounds_col
            for k in 1:i_litr_max
                nf.dwt_frootn_to_litr_n_col[c,j,k] = 0.0
            end
            nf.dwt_livecrootn_to_cwdn_col[c,j] = 0.0
            nf.dwt_deadcrootn_to_cwdn_col[c,j] = 0.0
        end
    end

    nothing
end

# ==========================================================================
# ZeroGRU — zero gross unrepresented landcover change fluxes
# Ported from: ZeroGru subroutine
# ==========================================================================

function cnveg_nitrogen_flux_zero_gru!(nf::CNVegNitrogenFluxData,
                                        bounds_grc::UnitRange{Int})
    for g in bounds_grc
        nf.gru_conv_nflux_grc[g]         = 0.0
        nf.gru_wood_productn_gain_grc[g]  = 0.0
    end

    nothing
end

# ==========================================================================
# Summary — compute diagnostic summary fluxes
# Ported from: Summary_nitrogenflux subroutine
# ==========================================================================

function cnveg_nitrogen_flux_summary!(nf::CNVegNitrogenFluxData,
                                       mask_patch::BitVector,
                                       bounds_patch::UnitRange{Int})
    for p in bounds_patch
        mask_patch[p] || continue

        # total N deployment
        nf.ndeploy_patch[p] =
            nf.sminn_to_npool_patch[p] +
            nf.retransn_to_npool_patch[p] +
            nf.free_retransn_to_npool_patch[p]

        # total patch-level fire N losses
        nf.fire_nloss_patch[p] =
            nf.m_leafn_to_fire_patch[p]              +
            nf.m_leafn_storage_to_fire_patch[p]      +
            nf.m_leafn_xfer_to_fire_patch[p]         +
            nf.m_frootn_to_fire_patch[p]             +
            nf.m_frootn_storage_to_fire_patch[p]     +
            nf.m_frootn_xfer_to_fire_patch[p]        +
            nf.m_livestemn_to_fire_patch[p]          +
            nf.m_livestemn_storage_to_fire_patch[p]  +
            nf.m_livestemn_xfer_to_fire_patch[p]     +
            nf.m_deadstemn_to_fire_patch[p]          +
            nf.m_deadstemn_storage_to_fire_patch[p]  +
            nf.m_deadstemn_xfer_to_fire_patch[p]     +
            nf.m_livecrootn_to_fire_patch[p]         +
            nf.m_livecrootn_storage_to_fire_patch[p] +
            nf.m_livecrootn_xfer_to_fire_patch[p]    +
            nf.m_deadcrootn_to_fire_patch[p]         +
            nf.m_deadcrootn_storage_to_fire_patch[p] +
            nf.m_deadcrootn_xfer_to_fire_patch[p]    +
            nf.m_retransn_to_fire_patch[p]

        # GRU conversion flux
        nf.gru_conv_nflux_patch[p] =
            nf.gru_livestemn_to_atm_patch[p]         +
            nf.gru_deadstemn_to_atm_patch[p]         +
            nf.gru_leafn_storage_to_atm_patch[p]     +
            nf.gru_frootn_storage_to_atm_patch[p]    +
            nf.gru_livestemn_storage_to_atm_patch[p] +
            nf.gru_deadstemn_storage_to_atm_patch[p] +
            nf.gru_livecrootn_storage_to_atm_patch[p] +
            nf.gru_deadcrootn_storage_to_atm_patch[p] +
            nf.gru_leafn_xfer_to_atm_patch[p]        +
            nf.gru_frootn_xfer_to_atm_patch[p]       +
            nf.gru_livestemn_xfer_to_atm_patch[p]    +
            nf.gru_deadstemn_xfer_to_atm_patch[p]    +
            nf.gru_livecrootn_xfer_to_atm_patch[p]   +
            nf.gru_deadcrootn_xfer_to_atm_patch[p]
    end

    nothing
end

# ==========================================================================
# InitCold — cold-start initialization
# Ported from: InitCold subroutine (simplified — no landunit filtering)
# ==========================================================================

function cnveg_nitrogen_flux_init_cold!(nf::CNVegNitrogenFluxData,
                                         bounds_patch::UnitRange{Int},
                                         bounds_col::UnitRange{Int};
                                         nlevdecomp_full::Int=1,
                                         ndecomp_pools::Int=1,
                                         i_litr_max::Int=I_LITR1)
    for p in bounds_patch
        nf.fert_counter_patch[p] = 0.0
        nf.fert_patch[p]        = 0.0
        nf.soyfixn_patch[p]     = 0.0
    end

    # Zero DWT/GRU column-level arrays
    for j in 1:nlevdecomp_full
        for c in bounds_col
            for k in 1:i_litr_max
                nf.dwt_frootn_to_litr_n_col[c,j,k] = 0.0
                nf.gru_n_to_litr_n_col[c,j,k]      = 0.0
            end
            nf.dwt_livecrootn_to_cwdn_col[c,j] = 0.0
            nf.dwt_deadcrootn_to_cwdn_col[c,j] = 0.0
            nf.gru_n_to_cwdn_col[c,j]          = 0.0
        end
    end

    nothing
end

# ==========================================================================
# Stubs for InitHistory and Restart (no-ops in Julia port)
# ==========================================================================

cnveg_nitrogen_flux_init_history!(nf::CNVegNitrogenFluxData) = nothing
cnveg_nitrogen_flux_restart!(nf::CNVegNitrogenFluxData) = nothing
