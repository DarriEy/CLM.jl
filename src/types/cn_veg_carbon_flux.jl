# ==========================================================================
# Ported from: src/biogeochem/CNVegCarbonFluxType.F90
# Vegetation carbon flux data type allocation and initialization
# ==========================================================================

"""
    CNVegCarbonFluxData

Vegetation carbon flux data structure. Holds carbon flux variables at patch,
column, and gridcell levels.

Ported from `cnveg_carbonflux_type` in `CNVegCarbonFluxType.F90`.
"""
Base.@kwdef mutable struct CNVegCarbonFluxData{FT<:Real,
                                   V<:AbstractVector{FT},
                                   M<:AbstractMatrix{FT},
                                   A3<:AbstractArray{FT,3},
                                   VI<:AbstractVector{<:Integer}}
    dribble_crophrv_xsmrpool_2atm::Bool = false

    # --- Gap mortality fluxes (gC/m2/s) patch-level ---
    m_leafc_to_litter_patch                   ::V = Float64[]
    m_leafc_storage_to_litter_patch           ::V = Float64[]
    m_leafc_xfer_to_litter_patch              ::V = Float64[]
    m_frootc_to_litter_patch                  ::V = Float64[]
    m_frootc_storage_to_litter_patch          ::V = Float64[]
    m_frootc_xfer_to_litter_patch             ::V = Float64[]
    m_livestemc_to_litter_patch               ::V = Float64[]
    m_livestemc_storage_to_litter_patch       ::V = Float64[]
    m_livestemc_xfer_to_litter_patch          ::V = Float64[]
    m_deadstemc_to_litter_patch               ::V = Float64[]
    m_deadstemc_storage_to_litter_patch       ::V = Float64[]
    m_deadstemc_xfer_to_litter_patch          ::V = Float64[]
    m_livecrootc_to_litter_patch              ::V = Float64[]
    m_livecrootc_storage_to_litter_patch      ::V = Float64[]
    m_livecrootc_xfer_to_litter_patch         ::V = Float64[]
    m_deadcrootc_to_litter_patch              ::V = Float64[]
    m_deadcrootc_storage_to_litter_patch      ::V = Float64[]
    m_deadcrootc_xfer_to_litter_patch         ::V = Float64[]
    m_gresp_storage_to_litter_patch           ::V = Float64[]
    m_gresp_xfer_to_litter_patch              ::V = Float64[]

    # --- Harvest mortality fluxes (gC/m2/s) patch-level ---
    hrv_leafc_to_litter_patch                 ::V = Float64[]
    hrv_leafc_storage_to_litter_patch         ::V = Float64[]
    hrv_leafc_xfer_to_litter_patch            ::V = Float64[]
    hrv_frootc_to_litter_patch                ::V = Float64[]
    hrv_frootc_storage_to_litter_patch        ::V = Float64[]
    hrv_frootc_xfer_to_litter_patch           ::V = Float64[]
    hrv_livestemc_to_litter_patch             ::V = Float64[]
    hrv_livestemc_storage_to_litter_patch     ::V = Float64[]
    hrv_livestemc_xfer_to_litter_patch        ::V = Float64[]
    hrv_deadstemc_storage_to_litter_patch     ::V = Float64[]
    hrv_deadstemc_xfer_to_litter_patch        ::V = Float64[]
    hrv_livecrootc_to_litter_patch            ::V = Float64[]
    hrv_livecrootc_storage_to_litter_patch    ::V = Float64[]
    hrv_livecrootc_xfer_to_litter_patch       ::V = Float64[]
    hrv_deadcrootc_to_litter_patch            ::V = Float64[]
    hrv_deadcrootc_storage_to_litter_patch    ::V = Float64[]
    hrv_deadcrootc_xfer_to_litter_patch       ::V = Float64[]
    hrv_gresp_storage_to_litter_patch         ::V = Float64[]
    hrv_gresp_xfer_to_litter_patch            ::V = Float64[]
    hrv_xsmrpool_to_atm_patch                 ::V = Float64[]

    # --- Fire fluxes (gC/m2/s) patch-level ---
    m_leafc_to_fire_patch                     ::V = Float64[]
    m_leafc_storage_to_fire_patch             ::V = Float64[]
    m_leafc_xfer_to_fire_patch                ::V = Float64[]
    m_livestemc_to_fire_patch                 ::V = Float64[]
    m_livestemc_storage_to_fire_patch         ::V = Float64[]
    m_livestemc_xfer_to_fire_patch            ::V = Float64[]
    m_deadstemc_to_fire_patch                 ::V = Float64[]
    m_deadstemc_storage_to_fire_patch         ::V = Float64[]
    m_deadstemc_xfer_to_fire_patch            ::V = Float64[]
    m_frootc_to_fire_patch                    ::V = Float64[]
    m_frootc_storage_to_fire_patch            ::V = Float64[]
    m_frootc_xfer_to_fire_patch               ::V = Float64[]
    m_livecrootc_to_fire_patch                ::V = Float64[]
    m_livecrootc_storage_to_fire_patch        ::V = Float64[]
    m_livecrootc_xfer_to_fire_patch           ::V = Float64[]
    m_deadcrootc_to_fire_patch                ::V = Float64[]
    m_deadcrootc_storage_to_fire_patch        ::V = Float64[]
    m_deadcrootc_xfer_to_fire_patch           ::V = Float64[]
    m_gresp_storage_to_fire_patch             ::V = Float64[]
    m_gresp_xfer_to_fire_patch                ::V = Float64[]

    # --- Fire-to-litter fluxes (gC/m2/s) patch-level ---
    m_leafc_to_litter_fire_patch              ::V = Float64[]
    m_leafc_storage_to_litter_fire_patch      ::V = Float64[]
    m_leafc_xfer_to_litter_fire_patch         ::V = Float64[]
    m_livestemc_to_litter_fire_patch          ::V = Float64[]
    m_livestemc_storage_to_litter_fire_patch  ::V = Float64[]
    m_livestemc_xfer_to_litter_fire_patch     ::V = Float64[]
    m_livestemc_to_deadstemc_fire_patch       ::V = Float64[]
    m_deadstemc_to_litter_fire_patch          ::V = Float64[]
    m_deadstemc_storage_to_litter_fire_patch  ::V = Float64[]
    m_deadstemc_xfer_to_litter_fire_patch     ::V = Float64[]
    m_frootc_to_litter_fire_patch             ::V = Float64[]
    m_frootc_storage_to_litter_fire_patch     ::V = Float64[]
    m_frootc_xfer_to_litter_fire_patch        ::V = Float64[]
    m_livecrootc_to_litter_fire_patch         ::V = Float64[]
    m_livecrootc_storage_to_litter_fire_patch ::V = Float64[]
    m_livecrootc_xfer_to_litter_fire_patch    ::V = Float64[]
    m_livecrootc_to_deadcrootc_fire_patch     ::V = Float64[]
    m_deadcrootc_to_litter_fire_patch         ::V = Float64[]
    m_deadcrootc_storage_to_litter_fire_patch ::V = Float64[]
    m_deadcrootc_xfer_to_litter_fire_patch    ::V = Float64[]
    m_gresp_storage_to_litter_fire_patch      ::V = Float64[]
    m_gresp_xfer_to_litter_fire_patch         ::V = Float64[]

    # --- Phenology fluxes (gC/m2/s) ---
    reproductivec_xfer_to_reproductivec_patch ::M = Matrix{Float64}(undef, 0, 0)
    leafc_xfer_to_leafc_patch                 ::V = Float64[]
    frootc_xfer_to_frootc_patch               ::V = Float64[]
    livestemc_xfer_to_livestemc_patch         ::V = Float64[]
    deadstemc_xfer_to_deadstemc_patch         ::V = Float64[]
    livecrootc_xfer_to_livecrootc_patch       ::V = Float64[]
    deadcrootc_xfer_to_deadcrootc_patch       ::V = Float64[]

    # --- Leaf/froot litterfall and crop fluxes (gC/m2/s) ---
    leafc_to_litter_patch                     ::V = Float64[]
    leafc_to_litter_fun_patch                 ::V = Float64[]
    frootc_to_litter_patch                    ::V = Float64[]
    livestemc_to_litter_patch                 ::V = Float64[]
    repr_grainc_to_food_patch                 ::M = Matrix{Float64}(undef, 0, 0)
    repr_grainc_to_food_perharv_patch         ::A3 = Array{Float64}(undef, 0, 0, 0)
    repr_grainc_to_food_thisyr_patch          ::M = Matrix{Float64}(undef, 0, 0)
    repr_structurec_to_cropprod_patch         ::M = Matrix{Float64}(undef, 0, 0)
    repr_structurec_to_litter_patch           ::M = Matrix{Float64}(undef, 0, 0)
    leafc_to_biofuelc_patch                   ::V = Float64[]
    livestemc_to_biofuelc_patch               ::V = Float64[]
    leafc_to_removedresiduec_patch            ::V = Float64[]
    livestemc_to_removedresiduec_patch        ::V = Float64[]
    repr_grainc_to_seed_patch                 ::M = Matrix{Float64}(undef, 0, 0)
    repr_grainc_to_seed_perharv_patch         ::A3 = Array{Float64}(undef, 0, 0, 0)
    repr_grainc_to_seed_thisyr_patch          ::M = Matrix{Float64}(undef, 0, 0)

    # --- Maintenance respiration (gC/m2/s) ---
    cpool_to_resp_patch                       ::V = Float64[]
    cpool_to_leafc_resp_patch                 ::V = Float64[]
    cpool_to_leafc_storage_resp_patch         ::V = Float64[]
    cpool_to_frootc_resp_patch                ::V = Float64[]
    cpool_to_frootc_storage_resp_patch        ::V = Float64[]
    cpool_to_livecrootc_resp_patch            ::V = Float64[]
    cpool_to_livecrootc_storage_resp_patch    ::V = Float64[]
    cpool_to_livestemc_resp_patch             ::V = Float64[]
    cpool_to_livestemc_storage_resp_patch     ::V = Float64[]
    leaf_mr_patch                             ::V = Float64[]
    froot_mr_patch                            ::V = Float64[]
    livestem_mr_patch                         ::V = Float64[]
    livecroot_mr_patch                        ::V = Float64[]
    reproductive_mr_patch                     ::M = Matrix{Float64}(undef, 0, 0)
    leaf_curmr_patch                          ::V = Float64[]
    froot_curmr_patch                         ::V = Float64[]
    livestem_curmr_patch                      ::V = Float64[]
    livecroot_curmr_patch                     ::V = Float64[]
    reproductive_curmr_patch                  ::M = Matrix{Float64}(undef, 0, 0)
    leaf_xsmr_patch                           ::V = Float64[]
    froot_xsmr_patch                          ::V = Float64[]
    livestem_xsmr_patch                       ::V = Float64[]
    livecroot_xsmr_patch                      ::V = Float64[]
    reproductive_xsmr_patch                   ::M = Matrix{Float64}(undef, 0, 0)

    # --- Photosynthesis (gC/m2/s) ---
    psnsun_to_cpool_patch                     ::V = Float64[]
    psnshade_to_cpool_patch                   ::V = Float64[]

    # --- Allocation from cpool (gC/m2/s) ---
    cpool_to_xsmrpool_patch                   ::V = Float64[]
    cpool_to_reproductivec_patch              ::M = Matrix{Float64}(undef, 0, 0)
    cpool_to_reproductivec_storage_patch      ::M = Matrix{Float64}(undef, 0, 0)
    cpool_to_leafc_patch                      ::V = Float64[]
    cpool_to_leafc_storage_patch              ::V = Float64[]
    cpool_to_frootc_patch                     ::V = Float64[]
    cpool_to_frootc_storage_patch             ::V = Float64[]
    cpool_to_livestemc_patch                  ::V = Float64[]
    cpool_to_livestemc_storage_patch          ::V = Float64[]
    cpool_to_deadstemc_patch                  ::V = Float64[]
    cpool_to_deadstemc_storage_patch          ::V = Float64[]
    cpool_to_livecrootc_patch                 ::V = Float64[]
    cpool_to_livecrootc_storage_patch         ::V = Float64[]
    cpool_to_deadcrootc_patch                 ::V = Float64[]
    cpool_to_deadcrootc_storage_patch         ::V = Float64[]
    cpool_to_gresp_storage_patch              ::V = Float64[]

    # --- Growth respiration (gC/m2/s) ---
    xsmrpool_to_atm_patch                     ::V = Float64[]
    xsmrpool_to_atm_col                       ::V = Float64[]
    xsmrpool_to_atm_grc                       ::V = Float64[]
    cpool_leaf_gr_patch                       ::V = Float64[]
    cpool_leaf_storage_gr_patch               ::V = Float64[]
    transfer_leaf_gr_patch                    ::V = Float64[]
    cpool_froot_gr_patch                      ::V = Float64[]
    cpool_froot_storage_gr_patch              ::V = Float64[]
    transfer_froot_gr_patch                   ::V = Float64[]
    cpool_livestem_gr_patch                   ::V = Float64[]
    cpool_livestem_storage_gr_patch           ::V = Float64[]
    transfer_livestem_gr_patch                ::V = Float64[]
    cpool_deadstem_gr_patch                   ::V = Float64[]
    cpool_deadstem_storage_gr_patch           ::V = Float64[]
    transfer_deadstem_gr_patch                ::V = Float64[]
    cpool_livecroot_gr_patch                  ::V = Float64[]
    cpool_livecroot_storage_gr_patch          ::V = Float64[]
    transfer_livecroot_gr_patch               ::V = Float64[]
    cpool_deadcroot_gr_patch                  ::V = Float64[]
    cpool_deadcroot_storage_gr_patch          ::V = Float64[]
    transfer_deadcroot_gr_patch               ::V = Float64[]
    cpool_reproductive_gr_patch               ::M = Matrix{Float64}(undef, 0, 0)
    cpool_reproductive_storage_gr_patch       ::M = Matrix{Float64}(undef, 0, 0)
    transfer_reproductive_gr_patch            ::M = Matrix{Float64}(undef, 0, 0)

    # --- Storage to transfer turnover (gC/m2/s) ---
    reproductivec_storage_to_xfer_patch       ::M = Matrix{Float64}(undef, 0, 0)
    leafc_storage_to_xfer_patch               ::V = Float64[]
    frootc_storage_to_xfer_patch              ::V = Float64[]
    livestemc_storage_to_xfer_patch           ::V = Float64[]
    deadstemc_storage_to_xfer_patch           ::V = Float64[]
    livecrootc_storage_to_xfer_patch          ::V = Float64[]
    deadcrootc_storage_to_xfer_patch          ::V = Float64[]
    gresp_storage_to_xfer_patch               ::V = Float64[]

    # --- Livewood to deadwood turnover (gC/m2/s) ---
    livestemc_to_deadstemc_patch              ::V = Float64[]
    livecrootc_to_deadcrootc_patch            ::V = Float64[]

    # --- Column-level decomposition fluxes (gC/m3/s) ---
    phenology_c_to_litr_c_col                 ::A3 = Array{Float64}(undef, 0, 0, 0)
    gap_mortality_c_to_litr_c_col             ::A3 = Array{Float64}(undef, 0, 0, 0)
    gap_mortality_c_to_cwdc_col               ::M = Matrix{Float64}(undef, 0, 0)
    fire_mortality_c_to_cwdc_col              ::M = Matrix{Float64}(undef, 0, 0)
    harvest_c_to_litr_c_col                   ::A3 = Array{Float64}(undef, 0, 0, 0)
    harvest_c_to_cwdc_col                     ::M = Matrix{Float64}(undef, 0, 0)
    crop_harvestc_to_cropprodc_patch          ::V = Float64[]
    crop_harvestc_to_cropprodc_col            ::V = Float64[]
    m_decomp_cpools_to_fire_vr_col            ::A3 = Array{Float64}(undef, 0, 0, 0)
    m_decomp_cpools_to_fire_col               ::M = Matrix{Float64}(undef, 0, 0)
    m_c_to_litr_fire_col                      ::A3 = Array{Float64}(undef, 0, 0, 0)

    # --- Dynamic landcover fluxes ---
    dwt_seedc_to_leaf_patch                   ::V = Float64[]
    dwt_seedc_to_leaf_grc                     ::V = Float64[]
    dwt_seedc_to_deadstem_patch               ::V = Float64[]
    dwt_seedc_to_deadstem_grc                 ::V = Float64[]
    dwt_conv_cflux_patch                      ::V = Float64[]
    dwt_conv_cflux_grc                        ::V = Float64[]
    dwt_conv_cflux_dribbled_grc               ::V = Float64[]
    dwt_wood_productc_gain_patch              ::V = Float64[]
    dwt_crop_productc_gain_patch              ::V = Float64[]
    dwt_slash_cflux_patch                     ::V = Float64[]
    dwt_slash_cflux_grc                       ::V = Float64[]
    dwt_frootc_to_litr_c_col                  ::A3 = Array{Float64}(undef, 0, 0, 0)
    dwt_livecrootc_to_cwdc_col                ::M = Matrix{Float64}(undef, 0, 0)
    dwt_deadcrootc_to_cwdc_col                ::M = Matrix{Float64}(undef, 0, 0)

    # --- Gross unrepresented landcover fluxes ---
    gru_leafc_to_litter_patch                 ::V = Float64[]
    gru_leafc_storage_to_atm_patch            ::V = Float64[]
    gru_leafc_xfer_to_atm_patch               ::V = Float64[]
    gru_frootc_to_litter_patch                ::V = Float64[]
    gru_frootc_storage_to_atm_patch           ::V = Float64[]
    gru_frootc_xfer_to_atm_patch              ::V = Float64[]
    gru_livestemc_to_atm_patch                ::V = Float64[]
    gru_livestemc_storage_to_atm_patch        ::V = Float64[]
    gru_livestemc_xfer_to_atm_patch           ::V = Float64[]
    gru_deadstemc_to_atm_patch                ::V = Float64[]
    gru_deadstemc_storage_to_atm_patch        ::V = Float64[]
    gru_deadstemc_xfer_to_atm_patch           ::V = Float64[]
    gru_livecrootc_to_litter_patch            ::V = Float64[]
    gru_livecrootc_storage_to_atm_patch       ::V = Float64[]
    gru_livecrootc_xfer_to_atm_patch          ::V = Float64[]
    gru_deadcrootc_to_litter_patch            ::V = Float64[]
    gru_deadcrootc_storage_to_atm_patch       ::V = Float64[]
    gru_deadcrootc_xfer_to_atm_patch          ::V = Float64[]
    gru_gresp_storage_to_atm_patch            ::V = Float64[]
    gru_gresp_xfer_to_atm_patch               ::V = Float64[]
    gru_xsmrpool_to_atm_patch                 ::V = Float64[]
    gru_conv_cflux_patch                      ::V = Float64[]
    gru_conv_cflux_col                        ::V = Float64[]
    gru_conv_cflux_grc                        ::V = Float64[]
    gru_conv_cflux_dribbled_grc               ::V = Float64[]
    gru_wood_productc_gain_patch              ::V = Float64[]
    gru_wood_productc_gain_col                ::V = Float64[]
    gru_slash_cflux_patch                     ::V = Float64[]
    gru_c_to_litr_c_col                       ::A3 = Array{Float64}(undef, 0, 0, 0)
    gru_c_to_cwdc_col                         ::M = Matrix{Float64}(undef, 0, 0)

    # --- Crop fluxes ---
    crop_seedc_to_leaf_patch                  ::V = Float64[]

    # --- Summary/diagnostic fluxes (gC/m2/s) ---
    gpp_before_downreg_patch                  ::V = Float64[]
    current_gr_patch                          ::V = Float64[]
    transfer_gr_patch                         ::V = Float64[]
    storage_gr_patch                          ::V = Float64[]
    plant_calloc_patch                        ::V = Float64[]
    excess_cflux_patch                        ::V = Float64[]
    prev_leafc_to_litter_patch                ::V = Float64[]
    prev_frootc_to_litter_patch               ::V = Float64[]
    availc_patch                              ::V = Float64[]
    xsmrpool_recover_patch                    ::V = Float64[]
    xsmrpool_c13ratio_patch                   ::V = Float64[]
    cwdc_loss_col                             ::V = Float64[]
    litterc_loss_col                          ::V = Float64[]
    frootc_alloc_patch                        ::V = Float64[]
    frootc_loss_patch                         ::V = Float64[]
    leafc_alloc_patch                         ::V = Float64[]
    leafc_loss_patch                          ::V = Float64[]
    woodc_alloc_patch                         ::V = Float64[]
    woodc_loss_patch                          ::V = Float64[]
    gpp_patch                                 ::V = Float64[]
    gpp_col                                   ::V = Float64[]
    rr_patch                                  ::V = Float64[]
    rr_col                                    ::V = Float64[]
    mr_patch                                  ::V = Float64[]
    gr_patch                                  ::V = Float64[]
    ar_patch                                  ::V = Float64[]
    ar_col                                    ::V = Float64[]
    npp_patch                                 ::V = Float64[]
    npp_col                                   ::V = Float64[]
    agnpp_patch                               ::V = Float64[]
    bgnpp_patch                               ::V = Float64[]
    litfall_patch                             ::V = Float64[]
    wood_harvestc_patch                       ::V = Float64[]
    wood_harvestc_col                         ::V = Float64[]
    slash_harvestc_patch                      ::V = Float64[]
    cinputs_patch                             ::V = Float64[]
    coutputs_patch                            ::V = Float64[]
    sr_col                                    ::V = Float64[]
    er_col                                    ::V = Float64[]
    litfire_col                               ::V = Float64[]
    somfire_col                               ::V = Float64[]
    totfire_col                               ::V = Float64[]
    hrv_xsmrpool_to_atm_col                  ::V = Float64[]
    hrv_xsmrpool_to_atm_grc                  ::V = Float64[]

    # --- Fire summary ---
    fire_closs_patch                          ::V = Float64[]
    fire_closs_p2c_col                        ::V = Float64[]
    fire_closs_col                            ::V = Float64[]

    # --- Annual sums ---
    tempsum_litfall_patch                     ::V = Float64[]
    annsum_litfall_patch                      ::V = Float64[]
    tempsum_npp_patch                         ::V = Float64[]
    annsum_npp_patch                          ::V = Float64[]
    annsum_npp_col                            ::V = Float64[]
    lag_npp_col                               ::V = Float64[]

    # --- Summary C fluxes ---
    nep_col                                   ::V = Float64[]
    nbp_grc                                   ::V = Float64[]
    nee_grc                                   ::V = Float64[]
    landuseflux_grc                           ::V = Float64[]

    # --- FUN fluxes (gC/m2/s) patch-level ---
    npp_Nactive_patch                         ::V = Float64[]
    npp_burnedoff_patch                       ::V = Float64[]
    npp_Nnonmyc_patch                         ::V = Float64[]
    npp_Nam_patch                             ::V = Float64[]
    npp_Necm_patch                            ::V = Float64[]
    npp_Nactive_no3_patch                     ::V = Float64[]
    npp_Nactive_nh4_patch                     ::V = Float64[]
    npp_Nnonmyc_no3_patch                     ::V = Float64[]
    npp_Nnonmyc_nh4_patch                     ::V = Float64[]
    npp_Nam_no3_patch                         ::V = Float64[]
    npp_Nam_nh4_patch                         ::V = Float64[]
    npp_Necm_no3_patch                        ::V = Float64[]
    npp_Necm_nh4_patch                        ::V = Float64[]
    npp_Nfix_patch                            ::V = Float64[]
    npp_Nretrans_patch                        ::V = Float64[]
    npp_Nuptake_patch                         ::V = Float64[]
    npp_growth_patch                          ::V = Float64[]
    leafc_change_patch                        ::V = Float64[]
    soilc_change_patch                        ::V = Float64[]

    # --- Matrix CN arrays ---
    matrix_Cinput_patch                       ::V = Float64[]
    matrix_C13input_patch                     ::V = Float64[]
    matrix_C14input_patch                     ::V = Float64[]
    matrix_alloc_patch                        ::M = Matrix{Float64}(undef, 0, 0)
    matrix_phtransfer_patch                   ::M = Matrix{Float64}(undef, 0, 0)
    matrix_phturnover_patch                   ::M = Matrix{Float64}(undef, 0, 0)
    matrix_phtransfer_doner_patch             ::VI = Int[]
    matrix_phtransfer_receiver_patch          ::VI = Int[]
    matrix_gmtransfer_patch                   ::M = Matrix{Float64}(undef, 0, 0)
    matrix_gmturnover_patch                   ::M = Matrix{Float64}(undef, 0, 0)
    matrix_gmtransfer_doner_patch             ::VI = Int[]
    matrix_gmtransfer_receiver_patch          ::VI = Int[]
    matrix_fitransfer_patch                   ::M = Matrix{Float64}(undef, 0, 0)
    matrix_fiturnover_patch                   ::M = Matrix{Float64}(undef, 0, 0)
    matrix_fitransfer_doner_patch             ::VI = Int[]
    matrix_fitransfer_receiver_patch          ::VI = Int[]

    # --- Matrix index scalars ---
    ileafst_to_ileafxf_ph::Int = 0
    ileafxf_to_ileaf_ph::Int = 0
    ifrootst_to_ifrootxf_ph::Int = 0
    ifrootxf_to_ifroot_ph::Int = 0
    ilivestemst_to_ilivestemxf_ph::Int = 0
    ilivestemxf_to_ilivestem_ph::Int = 0
    ideadstemst_to_ideadstemxf_ph::Int = 0
    ideadstemxf_to_ideadstem_ph::Int = 0
    ilivecrootst_to_ilivecrootxf_ph::Int = 0
    ilivecrootxf_to_ilivecroot_ph::Int = 0
    ideadcrootst_to_ideadcrootxf_ph::Int = 0
    ideadcrootxf_to_ideadcroot_ph::Int = 0
    ilivestem_to_ideadstem_ph::Int = 0
    ilivecroot_to_ideadcroot_ph::Int = 0
    ileaf_to_iout_ph::Int = 0
    ifroot_to_iout_ph::Int = 0
    ilivestem_to_iout_ph::Int = 0
    igrain_to_iout_ph::Int = 0
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
end

CNVegCarbonFluxData{FT}(; kwargs...) where {FT<:Real} =
    CNVegCarbonFluxData{FT, Vector{FT}, Matrix{FT}, Array{FT,3}, Vector{Int}}(; kwargs...)
Adapt.@adapt_structure CNVegCarbonFluxData


# ==========================================================================
# Init function — allocates all arrays
# Ported from: InitAllocate subroutine
# ==========================================================================

"""
    cnveg_carbon_flux_init!(cf, np, nc, ng; use_matrixcn=false, nrepr=NREPR,
                            nlevdecomp_full=1, ndecomp_pools=1, mxharvests=MXHARVESTS)

Allocate and initialize all carbon flux arrays.
"""
function cnveg_carbon_flux_init!(cf::CNVegCarbonFluxData{FT}, np::Int, nc::Int, ng::Int;
                                  use_matrixcn::Bool=false,
                                  nrepr::Int=NREPR,
                                  nlevdecomp_full::Int=1,
                                  ndecomp_pools::Int=1,
                                  mxharvests::Int=MXHARVESTS,
                                  nvegcpool::Int=NVEGPOOL_NATVEG) where {FT}
    # Helper for NaN-filled arrays
    nanvec(n) = fill(FT(NaN), n)
    nanmat(r, c) = fill(FT(NaN), r, c)
    nan3d(a, b, c) = fill(FT(NaN), a, b, c)

    # --- Gap mortality ---
    cf.m_leafc_to_litter_patch              = nanvec(np)
    cf.m_leafc_storage_to_litter_patch      = nanvec(np)
    cf.m_leafc_xfer_to_litter_patch         = nanvec(np)
    cf.m_frootc_to_litter_patch             = nanvec(np)
    cf.m_frootc_storage_to_litter_patch     = nanvec(np)
    cf.m_frootc_xfer_to_litter_patch        = nanvec(np)
    cf.m_livestemc_to_litter_patch          = nanvec(np)
    cf.m_livestemc_storage_to_litter_patch  = nanvec(np)
    cf.m_livestemc_xfer_to_litter_patch     = nanvec(np)
    cf.m_deadstemc_to_litter_patch          = nanvec(np)
    cf.m_deadstemc_storage_to_litter_patch  = nanvec(np)
    cf.m_deadstemc_xfer_to_litter_patch     = nanvec(np)
    cf.m_livecrootc_to_litter_patch         = nanvec(np)
    cf.m_livecrootc_storage_to_litter_patch = nanvec(np)
    cf.m_livecrootc_xfer_to_litter_patch    = nanvec(np)
    cf.m_deadcrootc_to_litter_patch         = nanvec(np)
    cf.m_deadcrootc_storage_to_litter_patch = nanvec(np)
    cf.m_deadcrootc_xfer_to_litter_patch    = nanvec(np)
    cf.m_gresp_storage_to_litter_patch      = nanvec(np)
    cf.m_gresp_xfer_to_litter_patch         = nanvec(np)

    # --- Harvest mortality ---
    cf.hrv_leafc_to_litter_patch            = nanvec(np)
    cf.hrv_leafc_storage_to_litter_patch    = nanvec(np)
    cf.hrv_leafc_xfer_to_litter_patch       = nanvec(np)
    cf.hrv_frootc_to_litter_patch           = nanvec(np)
    cf.hrv_frootc_storage_to_litter_patch   = nanvec(np)
    cf.hrv_frootc_xfer_to_litter_patch      = nanvec(np)
    cf.hrv_livestemc_to_litter_patch        = nanvec(np)
    cf.hrv_livestemc_storage_to_litter_patch = nanvec(np)
    cf.hrv_livestemc_xfer_to_litter_patch   = nanvec(np)
    cf.hrv_deadstemc_storage_to_litter_patch = nanvec(np)
    cf.hrv_deadstemc_xfer_to_litter_patch   = nanvec(np)
    cf.hrv_livecrootc_to_litter_patch       = nanvec(np)
    cf.hrv_livecrootc_storage_to_litter_patch = nanvec(np)
    cf.hrv_livecrootc_xfer_to_litter_patch  = nanvec(np)
    cf.hrv_deadcrootc_to_litter_patch       = nanvec(np)
    cf.hrv_deadcrootc_storage_to_litter_patch = nanvec(np)
    cf.hrv_deadcrootc_xfer_to_litter_patch  = nanvec(np)
    cf.hrv_gresp_storage_to_litter_patch    = nanvec(np)
    cf.hrv_gresp_xfer_to_litter_patch       = nanvec(np)
    cf.hrv_xsmrpool_to_atm_patch            = zeros(FT, np)

    # --- Fire ---
    cf.m_leafc_to_fire_patch                = nanvec(np)
    cf.m_leafc_storage_to_fire_patch        = nanvec(np)
    cf.m_leafc_xfer_to_fire_patch           = nanvec(np)
    cf.m_livestemc_to_fire_patch            = nanvec(np)
    cf.m_livestemc_storage_to_fire_patch    = nanvec(np)
    cf.m_livestemc_xfer_to_fire_patch       = nanvec(np)
    cf.m_deadstemc_to_fire_patch            = nanvec(np)
    cf.m_deadstemc_storage_to_fire_patch    = nanvec(np)
    cf.m_deadstemc_xfer_to_fire_patch       = nanvec(np)
    cf.m_frootc_to_fire_patch               = nanvec(np)
    cf.m_frootc_storage_to_fire_patch       = nanvec(np)
    cf.m_frootc_xfer_to_fire_patch          = nanvec(np)
    cf.m_livecrootc_to_fire_patch           = nanvec(np)
    cf.m_livecrootc_storage_to_fire_patch   = nanvec(np)
    cf.m_livecrootc_xfer_to_fire_patch      = nanvec(np)
    cf.m_deadcrootc_to_fire_patch           = nanvec(np)
    cf.m_deadcrootc_storage_to_fire_patch   = nanvec(np)
    cf.m_deadcrootc_xfer_to_fire_patch      = nanvec(np)
    cf.m_gresp_storage_to_fire_patch        = nanvec(np)
    cf.m_gresp_xfer_to_fire_patch           = nanvec(np)

    # --- Fire-to-litter ---
    cf.m_leafc_to_litter_fire_patch              = nanvec(np)
    cf.m_leafc_storage_to_litter_fire_patch      = nanvec(np)
    cf.m_leafc_xfer_to_litter_fire_patch         = nanvec(np)
    cf.m_livestemc_to_litter_fire_patch          = nanvec(np)
    cf.m_livestemc_storage_to_litter_fire_patch  = nanvec(np)
    cf.m_livestemc_xfer_to_litter_fire_patch     = nanvec(np)
    cf.m_livestemc_to_deadstemc_fire_patch       = nanvec(np)
    cf.m_deadstemc_to_litter_fire_patch          = nanvec(np)
    cf.m_deadstemc_storage_to_litter_fire_patch  = nanvec(np)
    cf.m_deadstemc_xfer_to_litter_fire_patch     = nanvec(np)
    cf.m_frootc_to_litter_fire_patch             = nanvec(np)
    cf.m_frootc_storage_to_litter_fire_patch     = nanvec(np)
    cf.m_frootc_xfer_to_litter_fire_patch        = nanvec(np)
    cf.m_livecrootc_to_litter_fire_patch         = nanvec(np)
    cf.m_livecrootc_storage_to_litter_fire_patch = nanvec(np)
    cf.m_livecrootc_xfer_to_litter_fire_patch    = nanvec(np)
    cf.m_livecrootc_to_deadcrootc_fire_patch     = nanvec(np)
    cf.m_deadcrootc_to_litter_fire_patch         = nanvec(np)
    cf.m_deadcrootc_storage_to_litter_fire_patch = nanvec(np)
    cf.m_deadcrootc_xfer_to_litter_fire_patch    = nanvec(np)
    cf.m_gresp_storage_to_litter_fire_patch      = nanvec(np)
    cf.m_gresp_xfer_to_litter_fire_patch         = nanvec(np)

    # --- Phenology ---
    cf.reproductivec_xfer_to_reproductivec_patch = nanmat(np, nrepr)
    cf.leafc_xfer_to_leafc_patch            = nanvec(np)
    cf.frootc_xfer_to_frootc_patch          = nanvec(np)
    cf.livestemc_xfer_to_livestemc_patch    = nanvec(np)
    cf.deadstemc_xfer_to_deadstemc_patch    = nanvec(np)
    cf.livecrootc_xfer_to_livecrootc_patch  = nanvec(np)
    cf.deadcrootc_xfer_to_deadcrootc_patch  = nanvec(np)

    # --- Litterfall / crop ---
    cf.leafc_to_litter_patch                = nanvec(np)
    cf.leafc_to_litter_fun_patch            = nanvec(np)
    cf.frootc_to_litter_patch               = nanvec(np)
    cf.livestemc_to_litter_patch            = nanvec(np)
    cf.repr_grainc_to_food_patch            = nanmat(np, nrepr)
    cf.repr_grainc_to_food_perharv_patch    = nan3d(np, mxharvests, nrepr)
    cf.repr_grainc_to_food_thisyr_patch     = nanmat(np, nrepr)
    cf.repr_structurec_to_cropprod_patch    = nanmat(np, nrepr)
    cf.repr_structurec_to_litter_patch      = nanmat(np, nrepr)
    cf.leafc_to_biofuelc_patch              = nanvec(np)
    cf.livestemc_to_biofuelc_patch          = nanvec(np)
    cf.leafc_to_removedresiduec_patch       = nanvec(np)
    cf.livestemc_to_removedresiduec_patch   = nanvec(np)
    cf.repr_grainc_to_seed_patch            = nanmat(np, nrepr)
    cf.repr_grainc_to_seed_perharv_patch    = nan3d(np, mxharvests, nrepr)
    cf.repr_grainc_to_seed_thisyr_patch     = nanmat(np, nrepr)

    # --- Maintenance respiration ---
    cf.cpool_to_resp_patch                  = nanvec(np)
    cf.cpool_to_leafc_resp_patch            = nanvec(np)
    cf.cpool_to_leafc_storage_resp_patch    = nanvec(np)
    cf.cpool_to_frootc_resp_patch           = nanvec(np)
    cf.cpool_to_frootc_storage_resp_patch   = nanvec(np)
    cf.cpool_to_livecrootc_resp_patch       = nanvec(np)
    cf.cpool_to_livecrootc_storage_resp_patch = nanvec(np)
    cf.cpool_to_livestemc_resp_patch        = nanvec(np)
    cf.cpool_to_livestemc_storage_resp_patch = nanvec(np)
    cf.leaf_mr_patch                        = nanvec(np)
    cf.froot_mr_patch                       = nanvec(np)
    cf.livestem_mr_patch                    = nanvec(np)
    cf.livecroot_mr_patch                   = nanvec(np)
    cf.reproductive_mr_patch                = nanmat(np, nrepr)
    cf.leaf_curmr_patch                     = nanvec(np)
    cf.froot_curmr_patch                    = nanvec(np)
    cf.livestem_curmr_patch                 = nanvec(np)
    cf.livecroot_curmr_patch                = nanvec(np)
    cf.reproductive_curmr_patch             = nanmat(np, nrepr)
    cf.leaf_xsmr_patch                      = nanvec(np)
    cf.froot_xsmr_patch                     = nanvec(np)
    cf.livestem_xsmr_patch                  = nanvec(np)
    cf.livecroot_xsmr_patch                 = nanvec(np)
    cf.reproductive_xsmr_patch              = nanmat(np, nrepr)

    # --- Photosynthesis ---
    cf.psnsun_to_cpool_patch                = nanvec(np)
    cf.psnshade_to_cpool_patch              = nanvec(np)

    # --- Allocation ---
    cf.cpool_to_xsmrpool_patch              = nanvec(np)
    cf.cpool_to_reproductivec_patch         = nanmat(np, nrepr)
    cf.cpool_to_reproductivec_storage_patch = nanmat(np, nrepr)
    cf.cpool_to_leafc_patch                 = nanvec(np)
    cf.cpool_to_leafc_storage_patch         = nanvec(np)
    cf.cpool_to_frootc_patch                = nanvec(np)
    cf.cpool_to_frootc_storage_patch        = nanvec(np)
    cf.cpool_to_livestemc_patch             = nanvec(np)
    cf.cpool_to_livestemc_storage_patch     = nanvec(np)
    cf.cpool_to_deadstemc_patch             = nanvec(np)
    cf.cpool_to_deadstemc_storage_patch     = nanvec(np)
    cf.cpool_to_livecrootc_patch            = nanvec(np)
    cf.cpool_to_livecrootc_storage_patch    = nanvec(np)
    cf.cpool_to_deadcrootc_patch            = nanvec(np)
    cf.cpool_to_deadcrootc_storage_patch    = nanvec(np)
    cf.cpool_to_gresp_storage_patch         = nanvec(np)

    # --- Growth respiration ---
    cf.xsmrpool_to_atm_patch               = zeros(FT, np)
    cf.xsmrpool_to_atm_col                  = zeros(FT, nc)
    cf.xsmrpool_to_atm_grc                  = zeros(FT, ng)
    cf.cpool_leaf_gr_patch                  = nanvec(np)
    cf.cpool_leaf_storage_gr_patch          = nanvec(np)
    cf.transfer_leaf_gr_patch               = nanvec(np)
    cf.cpool_froot_gr_patch                 = nanvec(np)
    cf.cpool_froot_storage_gr_patch         = nanvec(np)
    cf.transfer_froot_gr_patch              = nanvec(np)
    cf.cpool_livestem_gr_patch              = nanvec(np)
    cf.cpool_livestem_storage_gr_patch      = nanvec(np)
    cf.transfer_livestem_gr_patch           = nanvec(np)
    cf.cpool_deadstem_gr_patch              = nanvec(np)
    cf.cpool_deadstem_storage_gr_patch      = nanvec(np)
    cf.transfer_deadstem_gr_patch           = nanvec(np)
    cf.cpool_livecroot_gr_patch             = nanvec(np)
    cf.cpool_livecroot_storage_gr_patch     = nanvec(np)
    cf.transfer_livecroot_gr_patch          = nanvec(np)
    cf.cpool_deadcroot_gr_patch             = nanvec(np)
    cf.cpool_deadcroot_storage_gr_patch     = nanvec(np)
    cf.transfer_deadcroot_gr_patch          = nanvec(np)
    cf.cpool_reproductive_gr_patch          = nanmat(np, nrepr)
    cf.cpool_reproductive_storage_gr_patch  = nanmat(np, nrepr)
    cf.transfer_reproductive_gr_patch       = nanmat(np, nrepr)

    # --- Storage to transfer ---
    cf.reproductivec_storage_to_xfer_patch  = nanmat(np, nrepr)
    cf.leafc_storage_to_xfer_patch          = nanvec(np)
    cf.frootc_storage_to_xfer_patch         = nanvec(np)
    cf.livestemc_storage_to_xfer_patch      = nanvec(np)
    cf.deadstemc_storage_to_xfer_patch      = nanvec(np)
    cf.livecrootc_storage_to_xfer_patch     = nanvec(np)
    cf.deadcrootc_storage_to_xfer_patch     = nanvec(np)
    cf.gresp_storage_to_xfer_patch          = nanvec(np)

    # --- Livewood to deadwood ---
    cf.livestemc_to_deadstemc_patch         = nanvec(np)
    cf.livecrootc_to_deadcrootc_patch       = nanvec(np)

    # --- Column-level decomposition ---
    cf.phenology_c_to_litr_c_col            = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    cf.gap_mortality_c_to_litr_c_col        = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    cf.gap_mortality_c_to_cwdc_col          = nanmat(nc, nlevdecomp_full)
    cf.fire_mortality_c_to_cwdc_col         = nanmat(nc, nlevdecomp_full)
    cf.harvest_c_to_litr_c_col              = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    cf.harvest_c_to_cwdc_col                = nanmat(nc, nlevdecomp_full)
    cf.crop_harvestc_to_cropprodc_patch     = nanvec(np)
    cf.crop_harvestc_to_cropprodc_col       = nanvec(nc)
    cf.m_decomp_cpools_to_fire_vr_col       = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    cf.m_decomp_cpools_to_fire_col          = nanmat(nc, ndecomp_pools)
    cf.m_c_to_litr_fire_col                 = nan3d(nc, nlevdecomp_full, ndecomp_pools)

    # --- Dynamic landcover ---
    cf.dwt_seedc_to_leaf_patch              = nanvec(np)
    cf.dwt_seedc_to_leaf_grc                = nanvec(ng)
    cf.dwt_seedc_to_deadstem_patch          = nanvec(np)
    cf.dwt_seedc_to_deadstem_grc            = nanvec(ng)
    cf.dwt_conv_cflux_patch                 = nanvec(np)
    cf.dwt_conv_cflux_grc                   = nanvec(ng)
    cf.dwt_conv_cflux_dribbled_grc          = nanvec(ng)
    cf.dwt_wood_productc_gain_patch         = nanvec(np)
    cf.dwt_crop_productc_gain_patch         = nanvec(np)
    cf.dwt_slash_cflux_patch                = nanvec(np)
    cf.dwt_slash_cflux_grc                  = nanvec(ng)
    cf.dwt_frootc_to_litr_c_col             = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    cf.dwt_livecrootc_to_cwdc_col           = nanmat(nc, nlevdecomp_full)
    cf.dwt_deadcrootc_to_cwdc_col           = nanmat(nc, nlevdecomp_full)

    # --- GRU fluxes ---
    cf.gru_leafc_to_litter_patch            = nanvec(np)
    cf.gru_leafc_storage_to_atm_patch       = nanvec(np)
    cf.gru_leafc_xfer_to_atm_patch          = nanvec(np)
    cf.gru_frootc_to_litter_patch           = nanvec(np)
    cf.gru_frootc_storage_to_atm_patch      = nanvec(np)
    cf.gru_frootc_xfer_to_atm_patch         = nanvec(np)
    cf.gru_livestemc_to_atm_patch           = nanvec(np)
    cf.gru_livestemc_storage_to_atm_patch   = nanvec(np)
    cf.gru_livestemc_xfer_to_atm_patch      = nanvec(np)
    cf.gru_deadstemc_to_atm_patch           = nanvec(np)
    cf.gru_deadstemc_storage_to_atm_patch   = nanvec(np)
    cf.gru_deadstemc_xfer_to_atm_patch      = nanvec(np)
    cf.gru_livecrootc_to_litter_patch       = nanvec(np)
    cf.gru_livecrootc_storage_to_atm_patch  = nanvec(np)
    cf.gru_livecrootc_xfer_to_atm_patch     = nanvec(np)
    cf.gru_deadcrootc_to_litter_patch       = nanvec(np)
    cf.gru_deadcrootc_storage_to_atm_patch  = nanvec(np)
    cf.gru_deadcrootc_xfer_to_atm_patch     = nanvec(np)
    cf.gru_gresp_storage_to_atm_patch       = nanvec(np)
    cf.gru_gresp_xfer_to_atm_patch          = nanvec(np)
    cf.gru_xsmrpool_to_atm_patch            = nanvec(np)
    cf.gru_conv_cflux_patch                 = nanvec(np)
    cf.gru_conv_cflux_col                   = nanvec(nc)
    cf.gru_conv_cflux_grc                   = nanvec(ng)
    cf.gru_conv_cflux_dribbled_grc          = nanvec(ng)
    cf.gru_wood_productc_gain_patch         = nanvec(np)
    cf.gru_wood_productc_gain_col           = nanvec(nc)
    cf.gru_slash_cflux_patch                = nanvec(np)
    cf.gru_c_to_litr_c_col                  = nan3d(nc, nlevdecomp_full, ndecomp_pools)
    cf.gru_c_to_cwdc_col                    = nanmat(nc, nlevdecomp_full)

    # --- Crop ---
    cf.crop_seedc_to_leaf_patch             = nanvec(np)

    # --- Summary/diagnostic ---
    cf.gpp_before_downreg_patch             = nanvec(np)
    cf.current_gr_patch                     = nanvec(np)
    cf.transfer_gr_patch                    = nanvec(np)
    cf.storage_gr_patch                     = nanvec(np)
    cf.plant_calloc_patch                   = nanvec(np)
    cf.excess_cflux_patch                   = nanvec(np)
    cf.prev_leafc_to_litter_patch           = nanvec(np)
    cf.prev_frootc_to_litter_patch          = nanvec(np)
    cf.availc_patch                         = nanvec(np)
    cf.xsmrpool_recover_patch               = nanvec(np)
    cf.xsmrpool_c13ratio_patch              = nanvec(np)
    cf.cwdc_loss_col                        = nanvec(nc)
    cf.litterc_loss_col                     = nanvec(nc)
    cf.frootc_alloc_patch                   = nanvec(np)
    cf.frootc_loss_patch                    = nanvec(np)
    cf.leafc_alloc_patch                    = nanvec(np)
    cf.leafc_loss_patch                     = nanvec(np)
    cf.woodc_alloc_patch                    = nanvec(np)
    cf.woodc_loss_patch                     = nanvec(np)
    cf.rr_patch                             = nanvec(np)
    cf.mr_patch                             = nanvec(np)
    cf.gr_patch                             = nanvec(np)
    cf.ar_patch                             = nanvec(np)
    cf.npp_patch                            = nanvec(np)
    cf.agnpp_patch                          = nanvec(np)
    cf.bgnpp_patch                          = nanvec(np)
    cf.litfall_patch                        = nanvec(np)
    cf.wood_harvestc_patch                  = nanvec(np)
    cf.slash_harvestc_patch                 = nanvec(np)
    cf.cinputs_patch                        = nanvec(np)
    cf.coutputs_patch                       = nanvec(np)
    cf.gpp_patch                            = nanvec(np)
    cf.fire_closs_patch                     = nanvec(np)
    cf.sr_col                               = nanvec(nc)
    cf.er_col                               = nanvec(nc)
    cf.litfire_col                          = nanvec(nc)
    cf.somfire_col                          = nanvec(nc)
    cf.totfire_col                          = nanvec(nc)
    cf.rr_col                               = nanvec(nc)
    cf.ar_col                               = nanvec(nc)
    cf.gpp_col                              = nanvec(nc)
    cf.npp_col                              = nanvec(nc)
    cf.fire_closs_p2c_col                   = nanvec(nc)
    cf.fire_closs_col                       = nanvec(nc)
    cf.wood_harvestc_col                    = nanvec(nc)
    cf.hrv_xsmrpool_to_atm_col             = zeros(FT, nc)
    cf.hrv_xsmrpool_to_atm_grc             = zeros(FT, ng)
    cf.tempsum_npp_patch                    = nanvec(np)
    cf.annsum_npp_patch                     = nanvec(np)
    cf.tempsum_litfall_patch                = nanvec(np)
    cf.annsum_litfall_patch                 = nanvec(np)
    cf.annsum_npp_col                       = nanvec(nc)
    cf.lag_npp_col                          = fill(SPVAL, nc)
    cf.nep_col                              = nanvec(nc)
    cf.nbp_grc                              = nanvec(ng)
    cf.nee_grc                              = nanvec(ng)
    cf.landuseflux_grc                      = nanvec(ng)

    # --- FUN ---
    cf.npp_Nactive_patch                    = nanvec(np)
    cf.npp_burnedoff_patch                  = nanvec(np)
    cf.npp_Nnonmyc_patch                    = nanvec(np)
    cf.npp_Nam_patch                        = nanvec(np)
    cf.npp_Necm_patch                       = nanvec(np)
    cf.npp_Nactive_no3_patch                = nanvec(np)
    cf.npp_Nactive_nh4_patch                = nanvec(np)
    cf.npp_Nnonmyc_no3_patch                = nanvec(np)
    cf.npp_Nnonmyc_nh4_patch                = nanvec(np)
    cf.npp_Nam_no3_patch                    = nanvec(np)
    cf.npp_Nam_nh4_patch                    = nanvec(np)
    cf.npp_Necm_no3_patch                   = nanvec(np)
    cf.npp_Necm_nh4_patch                   = nanvec(np)
    cf.npp_Nfix_patch                       = nanvec(np)
    cf.npp_Nretrans_patch                   = nanvec(np)
    cf.npp_Nuptake_patch                    = nanvec(np)
    cf.npp_growth_patch                     = nanvec(np)
    cf.leafc_change_patch                   = nanvec(np)
    cf.soilc_change_patch                   = nanvec(np)

    # --- Matrix CN (optional) ---
    if use_matrixcn
        cf.matrix_Cinput_patch              = nanvec(np)
        cf.matrix_C13input_patch            = nanvec(np)
        cf.matrix_C14input_patch            = nanvec(np)
        cf.matrix_alloc_patch               = nanmat(np, nvegcpool)
        cf.matrix_phturnover_patch          = nanmat(np, nvegcpool)
        cf.matrix_gmturnover_patch          = nanmat(np, nvegcpool)
        cf.matrix_fiturnover_patch          = nanmat(np, nvegcpool)
    end

    nothing
end

# ==========================================================================
# SetValues — bulk-set all flux variables for given masks
# Ported from: SetValues subroutine
# ==========================================================================

"""
    cnveg_carbon_flux_set_values!(cf, mask_patch, value_patch, mask_col, value_column;
                                  use_matrixcn=false, use_crop=false, nrepr=NREPR, ...)

Set all carbon flux variables to given values for masked patches/columns.
"""
function cnveg_carbon_flux_set_values!(cf::CNVegCarbonFluxData,
                                       mask_patch::BitVector, value_patch::Real,
                                       mask_col::BitVector, value_column::Real;
                                       use_matrixcn::Bool=false,
                                       use_crop::Bool=false,
                                       nrepr::Int=NREPR,
                                       nlevdecomp_full::Int=1,
                                       ndecomp_pools::Int=1,
                                       nvegcpool::Int=NVEGPOOL_NATVEG) where {FT}
    # Patch-level 1D fields
    for i in eachindex(mask_patch)
        mask_patch[i] || continue

        cf.m_leafc_to_litter_patch[i]              = value_patch
        cf.m_frootc_to_litter_patch[i]             = value_patch
        cf.m_leafc_storage_to_litter_patch[i]      = value_patch
        cf.m_frootc_storage_to_litter_patch[i]     = value_patch
        cf.m_livestemc_storage_to_litter_patch[i]  = value_patch
        cf.m_deadstemc_storage_to_litter_patch[i]  = value_patch
        cf.m_livecrootc_storage_to_litter_patch[i] = value_patch
        cf.m_deadcrootc_storage_to_litter_patch[i] = value_patch
        cf.m_leafc_xfer_to_litter_patch[i]         = value_patch
        cf.m_frootc_xfer_to_litter_patch[i]        = value_patch
        cf.m_livestemc_xfer_to_litter_patch[i]     = value_patch
        cf.m_deadstemc_xfer_to_litter_patch[i]     = value_patch
        cf.m_livecrootc_xfer_to_litter_patch[i]    = value_patch
        cf.m_deadcrootc_xfer_to_litter_patch[i]    = value_patch
        cf.m_livestemc_to_litter_patch[i]          = value_patch
        cf.m_deadstemc_to_litter_patch[i]          = value_patch
        cf.m_livecrootc_to_litter_patch[i]         = value_patch
        cf.m_deadcrootc_to_litter_patch[i]         = value_patch
        cf.m_gresp_storage_to_litter_patch[i]      = value_patch
        cf.m_gresp_xfer_to_litter_patch[i]         = value_patch

        cf.hrv_leafc_to_litter_patch[i]            = value_patch
        cf.hrv_leafc_storage_to_litter_patch[i]    = value_patch
        cf.hrv_leafc_xfer_to_litter_patch[i]       = value_patch
        cf.hrv_frootc_to_litter_patch[i]           = value_patch
        cf.hrv_frootc_storage_to_litter_patch[i]   = value_patch
        cf.hrv_frootc_xfer_to_litter_patch[i]      = value_patch
        cf.hrv_livestemc_to_litter_patch[i]        = value_patch
        cf.hrv_livestemc_storage_to_litter_patch[i] = value_patch
        cf.hrv_livestemc_xfer_to_litter_patch[i]   = value_patch
        cf.hrv_deadstemc_storage_to_litter_patch[i] = value_patch
        cf.hrv_deadstemc_xfer_to_litter_patch[i]   = value_patch
        cf.hrv_livecrootc_to_litter_patch[i]       = value_patch
        cf.hrv_livecrootc_storage_to_litter_patch[i] = value_patch
        cf.hrv_livecrootc_xfer_to_litter_patch[i]  = value_patch
        cf.hrv_deadcrootc_to_litter_patch[i]       = value_patch
        cf.hrv_deadcrootc_storage_to_litter_patch[i] = value_patch
        cf.hrv_deadcrootc_xfer_to_litter_patch[i]  = value_patch
        cf.hrv_gresp_storage_to_litter_patch[i]    = value_patch
        cf.hrv_gresp_xfer_to_litter_patch[i]       = value_patch
        cf.hrv_xsmrpool_to_atm_patch[i]            = value_patch

        cf.gru_leafc_to_litter_patch[i]            = value_patch
        cf.gru_leafc_storage_to_atm_patch[i]       = value_patch
        cf.gru_leafc_xfer_to_atm_patch[i]          = value_patch
        cf.gru_frootc_to_litter_patch[i]           = value_patch
        cf.gru_frootc_storage_to_atm_patch[i]      = value_patch
        cf.gru_frootc_xfer_to_atm_patch[i]         = value_patch
        cf.gru_livestemc_to_atm_patch[i]           = value_patch
        cf.gru_livestemc_storage_to_atm_patch[i]   = value_patch
        cf.gru_livestemc_xfer_to_atm_patch[i]      = value_patch
        cf.gru_deadstemc_to_atm_patch[i]           = value_patch
        cf.gru_deadstemc_storage_to_atm_patch[i]   = value_patch
        cf.gru_deadstemc_xfer_to_atm_patch[i]      = value_patch
        cf.gru_livecrootc_to_litter_patch[i]       = value_patch
        cf.gru_livecrootc_storage_to_atm_patch[i]  = value_patch
        cf.gru_livecrootc_xfer_to_atm_patch[i]     = value_patch
        cf.gru_deadcrootc_to_litter_patch[i]       = value_patch
        cf.gru_deadcrootc_storage_to_atm_patch[i]  = value_patch
        cf.gru_deadcrootc_xfer_to_atm_patch[i]     = value_patch
        cf.gru_gresp_storage_to_atm_patch[i]       = value_patch
        cf.gru_gresp_xfer_to_atm_patch[i]          = value_patch
        cf.gru_xsmrpool_to_atm_patch[i]            = value_patch
        cf.gru_conv_cflux_patch[i]                 = value_patch
        cf.gru_wood_productc_gain_patch[i]         = value_patch
        cf.gru_slash_cflux_patch[i]                = value_patch

        cf.m_leafc_to_fire_patch[i]                = value_patch
        cf.m_leafc_storage_to_fire_patch[i]        = value_patch
        cf.m_leafc_xfer_to_fire_patch[i]           = value_patch
        cf.m_livestemc_to_fire_patch[i]            = value_patch
        cf.m_livestemc_storage_to_fire_patch[i]    = value_patch
        cf.m_livestemc_xfer_to_fire_patch[i]       = value_patch
        cf.m_deadstemc_to_fire_patch[i]            = value_patch
        cf.m_deadstemc_storage_to_fire_patch[i]    = value_patch
        cf.m_deadstemc_xfer_to_fire_patch[i]       = value_patch
        cf.m_frootc_to_fire_patch[i]               = value_patch
        cf.m_frootc_storage_to_fire_patch[i]       = value_patch
        cf.m_frootc_xfer_to_fire_patch[i]          = value_patch
        cf.m_livecrootc_to_fire_patch[i]           = value_patch
        cf.m_livecrootc_storage_to_fire_patch[i]   = value_patch
        cf.m_livecrootc_xfer_to_fire_patch[i]      = value_patch
        cf.m_deadcrootc_to_fire_patch[i]           = value_patch
        cf.m_deadcrootc_storage_to_fire_patch[i]   = value_patch
        cf.m_deadcrootc_xfer_to_fire_patch[i]      = value_patch
        cf.m_gresp_storage_to_fire_patch[i]        = value_patch
        cf.m_gresp_xfer_to_fire_patch[i]           = value_patch

        cf.m_leafc_to_litter_fire_patch[i]              = value_patch
        cf.m_leafc_storage_to_litter_fire_patch[i]      = value_patch
        cf.m_leafc_xfer_to_litter_fire_patch[i]         = value_patch
        cf.m_livestemc_to_litter_fire_patch[i]          = value_patch
        cf.m_livestemc_storage_to_litter_fire_patch[i]  = value_patch
        cf.m_livestemc_xfer_to_litter_fire_patch[i]     = value_patch
        cf.m_livestemc_to_deadstemc_fire_patch[i]       = value_patch
        cf.m_deadstemc_to_litter_fire_patch[i]          = value_patch
        cf.m_deadstemc_storage_to_litter_fire_patch[i]  = value_patch
        cf.m_deadstemc_xfer_to_litter_fire_patch[i]     = value_patch
        cf.m_frootc_to_litter_fire_patch[i]             = value_patch
        cf.m_frootc_storage_to_litter_fire_patch[i]     = value_patch
        cf.m_frootc_xfer_to_litter_fire_patch[i]        = value_patch
        cf.m_livecrootc_to_litter_fire_patch[i]         = value_patch
        cf.m_livecrootc_storage_to_litter_fire_patch[i] = value_patch
        cf.m_livecrootc_xfer_to_litter_fire_patch[i]    = value_patch
        cf.m_livecrootc_to_deadcrootc_fire_patch[i]     = value_patch
        cf.m_deadcrootc_to_litter_fire_patch[i]         = value_patch
        cf.m_deadcrootc_storage_to_litter_fire_patch[i] = value_patch
        cf.m_deadcrootc_xfer_to_litter_fire_patch[i]    = value_patch
        cf.m_gresp_storage_to_litter_fire_patch[i]      = value_patch
        cf.m_gresp_xfer_to_litter_fire_patch[i]         = value_patch

        cf.leafc_xfer_to_leafc_patch[i]            = value_patch
        cf.frootc_xfer_to_frootc_patch[i]          = value_patch
        cf.livestemc_xfer_to_livestemc_patch[i]    = value_patch
        cf.deadstemc_xfer_to_deadstemc_patch[i]    = value_patch
        cf.livecrootc_xfer_to_livecrootc_patch[i]  = value_patch
        cf.deadcrootc_xfer_to_deadcrootc_patch[i]  = value_patch
        cf.leafc_to_litter_patch[i]                = value_patch
        cf.frootc_to_litter_patch[i]               = value_patch

        cf.cpool_to_resp_patch[i]                  = value_patch
        cf.cpool_to_leafc_resp_patch[i]            = value_patch
        cf.cpool_to_leafc_storage_resp_patch[i]    = value_patch
        cf.cpool_to_frootc_resp_patch[i]           = value_patch
        cf.cpool_to_frootc_storage_resp_patch[i]   = value_patch
        cf.cpool_to_livecrootc_resp_patch[i]       = value_patch
        cf.cpool_to_livecrootc_storage_resp_patch[i] = value_patch
        cf.cpool_to_livestemc_resp_patch[i]        = value_patch
        cf.cpool_to_livestemc_storage_resp_patch[i] = value_patch
        cf.leaf_mr_patch[i]                        = value_patch
        cf.froot_mr_patch[i]                       = value_patch
        cf.livestem_mr_patch[i]                    = value_patch
        cf.livecroot_mr_patch[i]                   = value_patch
        cf.leaf_curmr_patch[i]                     = value_patch
        cf.froot_curmr_patch[i]                    = value_patch
        cf.livestem_curmr_patch[i]                 = value_patch
        cf.livecroot_curmr_patch[i]                = value_patch
        cf.leaf_xsmr_patch[i]                      = value_patch
        cf.froot_xsmr_patch[i]                     = value_patch
        cf.livestem_xsmr_patch[i]                  = value_patch
        cf.livecroot_xsmr_patch[i]                 = value_patch
        cf.psnsun_to_cpool_patch[i]                = value_patch
        cf.psnshade_to_cpool_patch[i]              = value_patch
        cf.cpool_to_xsmrpool_patch[i]              = value_patch
        cf.cpool_to_leafc_patch[i]                 = value_patch
        cf.cpool_to_leafc_storage_patch[i]         = value_patch
        cf.cpool_to_frootc_patch[i]                = value_patch
        cf.cpool_to_frootc_storage_patch[i]        = value_patch
        cf.cpool_to_livestemc_patch[i]             = value_patch
        cf.cpool_to_livestemc_storage_patch[i]     = value_patch
        cf.cpool_to_deadstemc_patch[i]             = value_patch
        cf.cpool_to_deadstemc_storage_patch[i]     = value_patch
        cf.cpool_to_livecrootc_patch[i]            = value_patch
        cf.cpool_to_livecrootc_storage_patch[i]    = value_patch
        cf.cpool_to_deadcrootc_patch[i]            = value_patch
        cf.cpool_to_deadcrootc_storage_patch[i]    = value_patch
        cf.cpool_to_gresp_storage_patch[i]         = value_patch

        cf.cpool_leaf_gr_patch[i]                  = value_patch
        cf.cpool_leaf_storage_gr_patch[i]          = value_patch
        cf.transfer_leaf_gr_patch[i]               = value_patch
        cf.cpool_froot_gr_patch[i]                 = value_patch
        cf.cpool_froot_storage_gr_patch[i]         = value_patch
        cf.transfer_froot_gr_patch[i]              = value_patch
        cf.cpool_livestem_gr_patch[i]              = value_patch
        cf.cpool_livestem_storage_gr_patch[i]      = value_patch
        cf.transfer_livestem_gr_patch[i]           = value_patch
        cf.cpool_deadstem_gr_patch[i]              = value_patch
        cf.cpool_deadstem_storage_gr_patch[i]      = value_patch
        cf.transfer_deadstem_gr_patch[i]           = value_patch
        cf.cpool_livecroot_gr_patch[i]             = value_patch
        cf.cpool_livecroot_storage_gr_patch[i]     = value_patch
        cf.transfer_livecroot_gr_patch[i]          = value_patch
        cf.cpool_deadcroot_gr_patch[i]             = value_patch
        cf.cpool_deadcroot_storage_gr_patch[i]     = value_patch
        cf.transfer_deadcroot_gr_patch[i]          = value_patch

        cf.leafc_storage_to_xfer_patch[i]          = value_patch
        cf.frootc_storage_to_xfer_patch[i]         = value_patch
        cf.livestemc_storage_to_xfer_patch[i]      = value_patch
        cf.deadstemc_storage_to_xfer_patch[i]      = value_patch
        cf.livecrootc_storage_to_xfer_patch[i]     = value_patch
        cf.deadcrootc_storage_to_xfer_patch[i]     = value_patch
        cf.gresp_storage_to_xfer_patch[i]          = value_patch
        cf.livestemc_to_deadstemc_patch[i]         = value_patch
        cf.livecrootc_to_deadcrootc_patch[i]       = value_patch

        cf.current_gr_patch[i]                     = value_patch
        cf.transfer_gr_patch[i]                    = value_patch
        cf.storage_gr_patch[i]                     = value_patch
        cf.frootc_alloc_patch[i]                   = value_patch
        cf.frootc_loss_patch[i]                    = value_patch
        cf.leafc_alloc_patch[i]                    = value_patch
        cf.leafc_loss_patch[i]                     = value_patch
        cf.woodc_alloc_patch[i]                    = value_patch
        cf.woodc_loss_patch[i]                     = value_patch
        cf.crop_seedc_to_leaf_patch[i]             = value_patch
        cf.crop_harvestc_to_cropprodc_patch[i]     = value_patch

        cf.gpp_patch[i]                            = value_patch
        cf.mr_patch[i]                             = value_patch
        cf.gr_patch[i]                             = value_patch
        cf.ar_patch[i]                             = value_patch
        cf.rr_patch[i]                             = value_patch
        cf.npp_patch[i]                            = value_patch
        cf.agnpp_patch[i]                          = value_patch
        cf.bgnpp_patch[i]                          = value_patch
        cf.litfall_patch[i]                        = value_patch
        cf.wood_harvestc_patch[i]                  = value_patch
        cf.slash_harvestc_patch[i]                 = value_patch
        cf.cinputs_patch[i]                        = value_patch
        cf.coutputs_patch[i]                       = value_patch
        cf.fire_closs_patch[i]                     = value_patch

        cf.npp_Nactive_patch[i]                    = value_patch
        cf.npp_burnedoff_patch[i]                  = value_patch
        cf.npp_Nnonmyc_patch[i]                    = value_patch
        cf.npp_Nam_patch[i]                        = value_patch
        cf.npp_Necm_patch[i]                       = value_patch
        cf.npp_Nactive_no3_patch[i]                = value_patch
        cf.npp_Nactive_nh4_patch[i]                = value_patch
        cf.npp_Nnonmyc_no3_patch[i]                = value_patch
        cf.npp_Nnonmyc_nh4_patch[i]                = value_patch
        cf.npp_Nam_no3_patch[i]                    = value_patch
        cf.npp_Nam_nh4_patch[i]                    = value_patch
        cf.npp_Necm_no3_patch[i]                   = value_patch
        cf.npp_Necm_nh4_patch[i]                   = value_patch
        cf.npp_Nfix_patch[i]                       = value_patch
        cf.npp_Nretrans_patch[i]                   = value_patch
        cf.npp_Nuptake_patch[i]                    = value_patch
        cf.npp_growth_patch[i]                     = value_patch
        cf.leafc_change_patch[i]                   = value_patch
        cf.soilc_change_patch[i]                   = value_patch

        if use_matrixcn
            cf.matrix_Cinput_patch[i]              = value_patch
            cf.matrix_C13input_patch[i]            = value_patch
            cf.matrix_C14input_patch[i]            = value_patch
        end
    end

    # 2D patch fields (reproductive)
    for k in 1:nrepr
        for i in eachindex(mask_patch)
            mask_patch[i] || continue
            cf.reproductive_mr_patch[i,k]          = value_patch
            cf.reproductive_curmr_patch[i,k]       = value_patch
            cf.reproductive_xsmr_patch[i,k]        = value_patch
        end
    end

    if use_crop
        for i in eachindex(mask_patch)
            mask_patch[i] || continue
            cf.xsmrpool_to_atm_patch[i]            = value_patch
            cf.livestemc_to_litter_patch[i]        = value_patch
            cf.leafc_to_biofuelc_patch[i]          = value_patch
            cf.livestemc_to_biofuelc_patch[i]      = value_patch
            cf.leafc_to_removedresiduec_patch[i]   = value_patch
            cf.livestemc_to_removedresiduec_patch[i] = value_patch
        end

        for k in 1:nrepr
            for i in eachindex(mask_patch)
                mask_patch[i] || continue
                cf.reproductivec_xfer_to_reproductivec_patch[i,k] = value_patch
                cf.cpool_to_reproductivec_patch[i,k]              = value_patch
                cf.cpool_to_reproductivec_storage_patch[i,k]      = value_patch
                cf.cpool_reproductive_gr_patch[i,k]               = value_patch
                cf.cpool_reproductive_storage_gr_patch[i,k]       = value_patch
                cf.transfer_reproductive_gr_patch[i,k]            = value_patch
                cf.reproductivec_storage_to_xfer_patch[i,k]       = value_patch
                cf.repr_grainc_to_food_patch[i,k]                 = value_patch
                cf.repr_grainc_to_seed_patch[i,k]                 = value_patch
                cf.repr_structurec_to_cropprod_patch[i,k]         = value_patch
                cf.repr_structurec_to_litter_patch[i,k]           = value_patch
            end
        end
    end

    # Column-level fields
    for j in 1:nlevdecomp_full
        for i in eachindex(mask_col)
            mask_col[i] || continue
            cf.gap_mortality_c_to_cwdc_col[i,j]    = value_column
            cf.fire_mortality_c_to_cwdc_col[i,j]   = value_column
            cf.harvest_c_to_cwdc_col[i,j]          = value_column
            cf.gru_c_to_cwdc_col[i,j]              = value_column
            for k in 1:ndecomp_pools
                cf.phenology_c_to_litr_c_col[i,j,k]     = value_column
                cf.gap_mortality_c_to_litr_c_col[i,j,k]  = value_column
                cf.harvest_c_to_litr_c_col[i,j,k]        = value_column
                cf.m_c_to_litr_fire_col[i,j,k]           = value_column
                cf.gru_c_to_litr_c_col[i,j,k]            = value_column
            end
        end
    end

    for k in 1:ndecomp_pools
        for j in 1:nlevdecomp_full
            for i in eachindex(mask_col)
                mask_col[i] || continue
                cf.m_decomp_cpools_to_fire_vr_col[i,j,k] = value_column
            end
        end
    end

    for k in 1:ndecomp_pools
        for i in eachindex(mask_col)
            mask_col[i] || continue
            cf.m_decomp_cpools_to_fire_col[i,k]   = value_column
        end
    end

    for i in eachindex(mask_col)
        mask_col[i] || continue
        cf.crop_harvestc_to_cropprodc_col[i]   = value_column
        cf.cwdc_loss_col[i]                    = value_column
        cf.litterc_loss_col[i]                 = value_column
        cf.sr_col[i]                           = value_column
        cf.er_col[i]                           = value_column
        cf.litfire_col[i]                      = value_column
        cf.somfire_col[i]                      = value_column
        cf.totfire_col[i]                      = value_column
        cf.rr_col[i]                           = value_column
        cf.ar_col[i]                           = value_column
        cf.gpp_col[i]                          = value_column
        cf.npp_col[i]                          = value_column
        cf.fire_closs_col[i]                   = value_column
        cf.wood_harvestc_col[i]                = value_column
        cf.hrv_xsmrpool_to_atm_col[i]         = value_column
        cf.gru_conv_cflux_col[i]               = value_column
        cf.gru_wood_productc_gain_col[i]       = value_column
        cf.nep_col[i]                          = value_column
        if use_crop
            cf.xsmrpool_to_atm_col[i]          = value_column
        end
    end

    nothing
end

# ==========================================================================
# ZeroDwt — zero dynamic landcover fluxes at start of timestep
# Ported from: ZeroDwt subroutine
# ==========================================================================

function cnveg_carbon_flux_zero_dwt!(cf::CNVegCarbonFluxData,
                                      bounds_grc::UnitRange{Int},
                                      bounds_col::UnitRange{Int};
                                      nlevdecomp_full::Int=1,
                                      i_litr_max::Int=I_LITR1)
    for g in bounds_grc
        cf.dwt_seedc_to_leaf_grc[g]     = 0.0
        cf.dwt_seedc_to_deadstem_grc[g] = 0.0
        cf.dwt_conv_cflux_grc[g]        = 0.0
        cf.dwt_slash_cflux_grc[g]       = 0.0
    end

    # litr index runs i_litr_min:i_litr_max (i_litr_min == 1), matching Fortran
    # CNVegCarbonFluxType::ZeroDwt — NOT 1:ndecomp_pools (dwt_frootc_to_litr_c_col's
    # 3rd dim is the litter sub-pool count, not the full decomp-pool count).
    for j in 1:nlevdecomp_full
        for c in bounds_col
            for k in 1:i_litr_max
                cf.dwt_frootc_to_litr_c_col[c,j,k] = 0.0
            end
            cf.dwt_livecrootc_to_cwdc_col[c,j] = 0.0
            cf.dwt_deadcrootc_to_cwdc_col[c,j] = 0.0
        end
    end

    nothing
end

# ==========================================================================
# ZeroGru — zero GRU fluxes at start of timestep
# Ported from: ZeroGru subroutine
# ==========================================================================

function cnveg_carbon_flux_zero_gru!(cf::CNVegCarbonFluxData,
                                      bounds_grc::UnitRange{Int})
    for g in bounds_grc
        cf.gru_conv_cflux_grc[g] = 0.0
    end
    nothing
end

# ==========================================================================
# InitCold — cold-start initialization
# Ported from: InitCold subroutine
# ==========================================================================

function cnveg_carbon_flux_init_cold!(cf::CNVegCarbonFluxData,
                                      bounds_patch::UnitRange{Int},
                                      bounds_col::UnitRange{Int};
                                      use_matrixcn::Bool=false,
                                      nlevdecomp_full::Int=1,
                                      ndecomp_pools::Int=1)
    for p in bounds_patch
        cf.gpp_before_downreg_patch[p] = 0.0
        cf.gpp_patch[p]                = 0.0
        cf.availc_patch[p]             = 0.0
        cf.xsmrpool_recover_patch[p]   = 0.0
        cf.excess_cflux_patch[p]       = 0.0
        cf.prev_leafc_to_litter_patch[p] = 0.0
        cf.leafc_to_litter_fun_patch[p]  = 0.0
        cf.prev_frootc_to_litter_patch[p] = 0.0
        cf.plant_calloc_patch[p]       = 0.0
        cf.tempsum_npp_patch[p]        = 0.0
        cf.annsum_npp_patch[p]         = 0.0
        cf.tempsum_litfall_patch[p]    = 0.0
        cf.annsum_litfall_patch[p]     = 0.0
        if use_matrixcn
            cf.matrix_Cinput_patch[p]   = 0.0
            cf.matrix_C13input_patch[p] = 0.0
            cf.matrix_C14input_patch[p] = 0.0
        end
    end

    for c in bounds_col
        cf.annsum_npp_col[c] = 0.0
        for j in 1:nlevdecomp_full
            for k in 1:ndecomp_pools
                cf.dwt_frootc_to_litr_c_col[c,j,k] = 0.0
                cf.gru_c_to_litr_c_col[c,j,k]      = 0.0
            end
            cf.dwt_livecrootc_to_cwdc_col[c,j] = 0.0
            cf.dwt_deadcrootc_to_cwdc_col[c,j] = 0.0
            cf.gru_c_to_cwdc_col[c,j]          = 0.0
        end
    end

    nothing
end

# ==========================================================================
# Summary — compute patch-level summary carbon flux variables
# Ported from: Summary_carbonflux subroutine
# ==========================================================================

"""
    cnveg_carbon_flux_summary!(cf, mask_patch, bounds_patch;
                                use_crop=false, use_fun=false,
                                carbon_resp_opt=0,
                                patch_itype=nothing, npcropmin=0,
                                nrepr=NREPR)

Compute patch-level carbon flux summary variables (MR, GR, AR, GPP, NPP, etc.).
Column-level p2c aggregation is stubbed pending subgridAveMod port.

Ported from `cnveg_carbonflux_type%Summary_carbonflux`.
"""
# --- GPU kernelization (one thread per patch). Every total is a LOCAL scalar written
#     ONCE; the crop reproductive `for k` loops accumulate into the mr/current_gr/
#     transfer_gr/storage_gr LOCALS (safe — the KA trap is only array-element `+=` in a
#     loop — and byte-identical, exact per-k order preserved). All touched flux fields
#     grouped into one @adapt_structure {V,M} bundle (one Metal buffer arg).
struct _CVCFSumView{V,M}
    # maintenance / growth respiration inputs
    leaf_mr::V; froot_mr::V; livestem_mr::V; livecroot_mr::V; cpool_to_resp::V
    cpool_leaf_gr::V; cpool_froot_gr::V; cpool_livestem_gr::V; cpool_deadstem_gr::V
    cpool_livecroot_gr::V; cpool_deadcroot_gr::V
    transfer_leaf_gr::V; transfer_froot_gr::V; transfer_livestem_gr::V
    transfer_deadstem_gr::V; transfer_livecroot_gr::V; transfer_deadcroot_gr::V
    cpool_leaf_storage_gr::V; cpool_froot_storage_gr::V; cpool_livestem_storage_gr::V
    cpool_deadstem_storage_gr::V; cpool_livecroot_storage_gr::V; cpool_deadcroot_storage_gr::V
    soilc_change::V; xsmrpool_to_atm::V; psnsun_to_cpool::V; psnshade_to_cpool::V
    # allocation inputs
    cpool_to_leafc::V; cpool_to_leafc_storage::V; leafc_xfer_to_leafc::V
    cpool_to_livestemc::V; cpool_to_livestemc_storage::V; livestemc_xfer_to_livestemc::V
    cpool_to_deadstemc::V; cpool_to_deadstemc_storage::V; deadstemc_xfer_to_deadstemc::V
    cpool_to_frootc::V; cpool_to_frootc_storage::V; frootc_xfer_to_frootc::V
    cpool_to_livecrootc::V; cpool_to_livecrootc_storage::V; livecrootc_xfer_to_livecrootc::V
    cpool_to_deadcrootc::V; cpool_to_deadcrootc_storage::V; deadcrootc_xfer_to_deadcrootc::V
    # litterfall inputs
    leafc_to_litter::V; frootc_to_litter::V
    m_leafc_to_litter::V; m_frootc_to_litter::V; m_leafc_storage_to_litter::V
    m_frootc_storage_to_litter::V; m_leafc_xfer_to_litter::V; m_frootc_xfer_to_litter::V
    m_livestemc_to_litter::V; m_livestemc_storage_to_litter::V; m_livestemc_xfer_to_litter::V
    m_deadstemc_to_litter::V; m_deadstemc_storage_to_litter::V; m_deadstemc_xfer_to_litter::V
    m_livecrootc_to_litter::V; m_livecrootc_storage_to_litter::V; m_livecrootc_xfer_to_litter::V
    m_deadcrootc_to_litter::V; m_deadcrootc_storage_to_litter::V; m_deadcrootc_xfer_to_litter::V
    m_gresp_storage_to_litter::V; m_gresp_xfer_to_litter::V
    # fire inputs
    m_leafc_to_fire::V; m_frootc_to_fire::V; m_leafc_storage_to_fire::V
    m_frootc_storage_to_fire::V; m_leafc_xfer_to_fire::V; m_frootc_xfer_to_fire::V
    m_livestemc_to_fire::V; m_livestemc_storage_to_fire::V; m_livestemc_xfer_to_fire::V
    m_deadstemc_to_fire::V; m_deadstemc_storage_to_fire::V; m_deadstemc_xfer_to_fire::V
    m_livecrootc_to_fire::V; m_livecrootc_storage_to_fire::V; m_livecrootc_xfer_to_fire::V
    m_deadcrootc_to_fire::V; m_deadcrootc_storage_to_fire::V; m_deadcrootc_xfer_to_fire::V
    m_gresp_storage_to_fire::V; m_gresp_xfer_to_fire::V
    # reproductive (crop) matrices
    reproductive_mr::M; cpool_reproductive_gr::M; transfer_reproductive_gr::M
    cpool_reproductive_storage_gr::M
    # outputs
    mr::V; current_gr::V; transfer_gr::V; storage_gr::V; gr::V; ar::V; gpp::V; npp::V
    rr::V; agnpp::V; bgnpp::V; litfall::V; fire_closs::V
    frootc_alloc::V; frootc_loss::V; leafc_alloc::V; leafc_loss::V; woodc_alloc::V; woodc_loss::V
end
Adapt.@adapt_structure _CVCFSumView

@kernel function _cvcf_summary_kernel!(@Const(mask), b, @Const(patch_itype),
        lo::Int, hi::Int, do_crop::Bool, use_fun::Bool, carbon_resp_opt::Int,
        npcropmin::Int, nrepr::Int)
    p = @index(Global)
    @inbounds if lo <= p <= hi && mask[p]
        crop_p = do_crop && patch_itype[p] >= npcropmin

        # --- Maintenance respiration (MR) ---
        mr = b.leaf_mr[p] + b.froot_mr[p] + b.livestem_mr[p] + b.livecroot_mr[p]
        if carbon_resp_opt == 1
            mr += b.cpool_to_resp[p]
        end
        if crop_p
            for k in 1:nrepr
                mr += b.reproductive_mr[p, k]
            end
        end

        # --- Current / transfer / storage growth respiration ---
        current_gr = b.cpool_leaf_gr[p] + b.cpool_froot_gr[p] + b.cpool_livestem_gr[p] +
                     b.cpool_deadstem_gr[p] + b.cpool_livecroot_gr[p] + b.cpool_deadcroot_gr[p]
        transfer_gr = b.transfer_leaf_gr[p] + b.transfer_froot_gr[p] + b.transfer_livestem_gr[p] +
                      b.transfer_deadstem_gr[p] + b.transfer_livecroot_gr[p] + b.transfer_deadcroot_gr[p]
        storage_gr = b.cpool_leaf_storage_gr[p] + b.cpool_froot_storage_gr[p] +
                     b.cpool_livestem_storage_gr[p] + b.cpool_deadstem_storage_gr[p] +
                     b.cpool_livecroot_storage_gr[p] + b.cpool_deadcroot_storage_gr[p]
        if crop_p
            for k in 1:nrepr
                current_gr  += b.cpool_reproductive_gr[p, k]
                transfer_gr += b.transfer_reproductive_gr[p, k]
                storage_gr  += b.cpool_reproductive_storage_gr[p, k]
            end
        end

        gr = current_gr + transfer_gr + storage_gr        # total growth respiration
        ar = mr + gr                                        # autotrophic respiration
        if use_fun
            ar += b.soilc_change[p]
        end
        if crop_p
            ar += b.xsmrpool_to_atm[p]
        end

        gpp = b.psnsun_to_cpool[p] + b.psnshade_to_cpool[p]
        npp = gpp - ar

        rr = b.froot_mr[p] + b.cpool_froot_gr[p] + b.cpool_livecroot_gr[p] +
             b.cpool_deadcroot_gr[p] + b.transfer_froot_gr[p] + b.transfer_livecroot_gr[p] +
             b.transfer_deadcroot_gr[p] + b.cpool_froot_storage_gr[p] +
             b.cpool_livecroot_storage_gr[p] + b.cpool_deadcroot_storage_gr[p]

        agnpp = b.cpool_to_leafc[p] + b.cpool_to_leafc_storage[p] + b.leafc_xfer_to_leafc[p] +
                b.cpool_to_livestemc[p] + b.cpool_to_livestemc_storage[p] + b.livestemc_xfer_to_livestemc[p] +
                b.cpool_to_deadstemc[p] + b.cpool_to_deadstemc_storage[p] + b.deadstemc_xfer_to_deadstemc[p]

        bgnpp = b.cpool_to_frootc[p] + b.cpool_to_frootc_storage[p] + b.frootc_xfer_to_frootc[p] +
                b.cpool_to_livecrootc[p] + b.cpool_to_livecrootc_storage[p] + b.livecrootc_xfer_to_livecrootc[p] +
                b.cpool_to_deadcrootc[p] + b.cpool_to_deadcrootc_storage[p] + b.deadcrootc_xfer_to_deadcrootc[p]

        litfall = b.leafc_to_litter[p] + b.frootc_to_litter[p] +
                  b.m_leafc_to_litter[p] + b.m_frootc_to_litter[p] + b.m_leafc_storage_to_litter[p] +
                  b.m_frootc_storage_to_litter[p] + b.m_leafc_xfer_to_litter[p] + b.m_frootc_xfer_to_litter[p] +
                  b.m_livestemc_to_litter[p] + b.m_livestemc_storage_to_litter[p] + b.m_livestemc_xfer_to_litter[p] +
                  b.m_deadstemc_to_litter[p] + b.m_deadstemc_storage_to_litter[p] + b.m_deadstemc_xfer_to_litter[p] +
                  b.m_livecrootc_to_litter[p] + b.m_livecrootc_storage_to_litter[p] + b.m_livecrootc_xfer_to_litter[p] +
                  b.m_deadcrootc_to_litter[p] + b.m_deadcrootc_storage_to_litter[p] + b.m_deadcrootc_xfer_to_litter[p] +
                  b.m_gresp_storage_to_litter[p] + b.m_gresp_xfer_to_litter[p]

        fire_closs = b.m_leafc_to_fire[p] + b.m_frootc_to_fire[p] + b.m_leafc_storage_to_fire[p] +
                     b.m_frootc_storage_to_fire[p] + b.m_leafc_xfer_to_fire[p] + b.m_frootc_xfer_to_fire[p] +
                     b.m_livestemc_to_fire[p] + b.m_livestemc_storage_to_fire[p] + b.m_livestemc_xfer_to_fire[p] +
                     b.m_deadstemc_to_fire[p] + b.m_deadstemc_storage_to_fire[p] + b.m_deadstemc_xfer_to_fire[p] +
                     b.m_livecrootc_to_fire[p] + b.m_livecrootc_storage_to_fire[p] + b.m_livecrootc_xfer_to_fire[p] +
                     b.m_deadcrootc_to_fire[p] + b.m_deadcrootc_storage_to_fire[p] + b.m_deadcrootc_xfer_to_fire[p] +
                     b.m_gresp_storage_to_fire[p] + b.m_gresp_xfer_to_fire[p]

        b.mr[p] = mr; b.current_gr[p] = current_gr; b.transfer_gr[p] = transfer_gr
        b.storage_gr[p] = storage_gr; b.gr[p] = gr; b.ar[p] = ar; b.gpp[p] = gpp
        b.npp[p] = npp; b.rr[p] = rr; b.agnpp[p] = agnpp; b.bgnpp[p] = bgnpp
        b.litfall[p] = litfall; b.fire_closs[p] = fire_closs

        b.frootc_alloc[p] = b.cpool_to_frootc[p] + b.cpool_to_frootc_storage[p]
        b.frootc_loss[p]  = b.leafc_to_litter[p] + b.frootc_to_litter[p]
        b.leafc_alloc[p]  = b.cpool_to_leafc[p] + b.cpool_to_leafc_storage[p]
        b.leafc_loss[p]   = b.leafc_to_litter[p]
        b.woodc_alloc[p]  = b.cpool_to_livestemc[p] + b.cpool_to_deadstemc[p] +
                            b.cpool_to_livecrootc[p] + b.cpool_to_deadcrootc[p] +
                            b.cpool_to_livestemc_storage[p] + b.cpool_to_deadstemc_storage[p] +
                            b.cpool_to_livecrootc_storage[p] + b.cpool_to_deadcrootc_storage[p]
        b.woodc_loss[p]   = b.m_livestemc_to_litter[p] + b.m_deadstemc_to_litter[p] +
                            b.m_livecrootc_to_litter[p] + b.m_deadcrootc_to_litter[p]
    end
end

function cnveg_carbon_flux_summary!(cf::CNVegCarbonFluxData,
                                     mask_patch::AbstractVector{Bool},
                                     bounds_patch::UnitRange{Int};
                                     use_crop::Bool=false,
                                     use_fun::Bool=false,
                                     carbon_resp_opt::Int=0,
                                     patch_itype::Union{AbstractVector{<:Integer},Nothing}=nothing,
                                     npcropmin::Int=0,
                                     nrepr::Int=NREPR)
    isempty(bounds_patch) && return nothing
    b = _CVCFSumView(
        cf.leaf_mr_patch, cf.froot_mr_patch, cf.livestem_mr_patch, cf.livecroot_mr_patch,
        cf.cpool_to_resp_patch, cf.cpool_leaf_gr_patch, cf.cpool_froot_gr_patch,
        cf.cpool_livestem_gr_patch, cf.cpool_deadstem_gr_patch, cf.cpool_livecroot_gr_patch,
        cf.cpool_deadcroot_gr_patch, cf.transfer_leaf_gr_patch, cf.transfer_froot_gr_patch,
        cf.transfer_livestem_gr_patch, cf.transfer_deadstem_gr_patch, cf.transfer_livecroot_gr_patch,
        cf.transfer_deadcroot_gr_patch, cf.cpool_leaf_storage_gr_patch, cf.cpool_froot_storage_gr_patch,
        cf.cpool_livestem_storage_gr_patch, cf.cpool_deadstem_storage_gr_patch,
        cf.cpool_livecroot_storage_gr_patch, cf.cpool_deadcroot_storage_gr_patch,
        cf.soilc_change_patch, cf.xsmrpool_to_atm_patch, cf.psnsun_to_cpool_patch, cf.psnshade_to_cpool_patch,
        cf.cpool_to_leafc_patch, cf.cpool_to_leafc_storage_patch, cf.leafc_xfer_to_leafc_patch,
        cf.cpool_to_livestemc_patch, cf.cpool_to_livestemc_storage_patch, cf.livestemc_xfer_to_livestemc_patch,
        cf.cpool_to_deadstemc_patch, cf.cpool_to_deadstemc_storage_patch, cf.deadstemc_xfer_to_deadstemc_patch,
        cf.cpool_to_frootc_patch, cf.cpool_to_frootc_storage_patch, cf.frootc_xfer_to_frootc_patch,
        cf.cpool_to_livecrootc_patch, cf.cpool_to_livecrootc_storage_patch, cf.livecrootc_xfer_to_livecrootc_patch,
        cf.cpool_to_deadcrootc_patch, cf.cpool_to_deadcrootc_storage_patch, cf.deadcrootc_xfer_to_deadcrootc_patch,
        cf.leafc_to_litter_patch, cf.frootc_to_litter_patch,
        cf.m_leafc_to_litter_patch, cf.m_frootc_to_litter_patch, cf.m_leafc_storage_to_litter_patch,
        cf.m_frootc_storage_to_litter_patch, cf.m_leafc_xfer_to_litter_patch, cf.m_frootc_xfer_to_litter_patch,
        cf.m_livestemc_to_litter_patch, cf.m_livestemc_storage_to_litter_patch, cf.m_livestemc_xfer_to_litter_patch,
        cf.m_deadstemc_to_litter_patch, cf.m_deadstemc_storage_to_litter_patch, cf.m_deadstemc_xfer_to_litter_patch,
        cf.m_livecrootc_to_litter_patch, cf.m_livecrootc_storage_to_litter_patch, cf.m_livecrootc_xfer_to_litter_patch,
        cf.m_deadcrootc_to_litter_patch, cf.m_deadcrootc_storage_to_litter_patch, cf.m_deadcrootc_xfer_to_litter_patch,
        cf.m_gresp_storage_to_litter_patch, cf.m_gresp_xfer_to_litter_patch,
        cf.m_leafc_to_fire_patch, cf.m_frootc_to_fire_patch, cf.m_leafc_storage_to_fire_patch,
        cf.m_frootc_storage_to_fire_patch, cf.m_leafc_xfer_to_fire_patch, cf.m_frootc_xfer_to_fire_patch,
        cf.m_livestemc_to_fire_patch, cf.m_livestemc_storage_to_fire_patch, cf.m_livestemc_xfer_to_fire_patch,
        cf.m_deadstemc_to_fire_patch, cf.m_deadstemc_storage_to_fire_patch, cf.m_deadstemc_xfer_to_fire_patch,
        cf.m_livecrootc_to_fire_patch, cf.m_livecrootc_storage_to_fire_patch, cf.m_livecrootc_xfer_to_fire_patch,
        cf.m_deadcrootc_to_fire_patch, cf.m_deadcrootc_storage_to_fire_patch, cf.m_deadcrootc_xfer_to_fire_patch,
        cf.m_gresp_storage_to_fire_patch, cf.m_gresp_xfer_to_fire_patch,
        cf.reproductive_mr_patch, cf.cpool_reproductive_gr_patch, cf.transfer_reproductive_gr_patch,
        cf.cpool_reproductive_storage_gr_patch,
        cf.mr_patch, cf.current_gr_patch, cf.transfer_gr_patch, cf.storage_gr_patch, cf.gr_patch,
        cf.ar_patch, cf.gpp_patch, cf.npp_patch, cf.rr_patch, cf.agnpp_patch, cf.bgnpp_patch,
        cf.litfall_patch, cf.fire_closs_patch, cf.frootc_alloc_patch, cf.frootc_loss_patch,
        cf.leafc_alloc_patch, cf.leafc_loss_patch, cf.woodc_alloc_patch, cf.woodc_loss_patch)
    T = eltype(cf.mr_patch)
    do_crop = use_crop && patch_itype !== nothing
    pit = do_crop ? _to_backend_like(cf.mr_patch, T, patch_itype) :
                    _to_backend_like(cf.mr_patch, T, Int[])
    _launch!(_cvcf_summary_kernel!, mask_patch, b, pit,
        first(bounds_patch), last(bounds_patch), do_crop, use_fun, carbon_resp_opt, npcropmin, nrepr)

    # The COLUMN/GRIDCELL half of Fortran's Summary_carbonflux lives in
    # `cnveg_carbon_flux_summary_col!` below (it needs the soil-BGC HR, which is
    # not available here). CNDriverSummarizeFluxes calls both, in that order.
    return nothing
end

# ==========================================================================
# cnveg_carbon_flux_summary_col! — the COLUMN + GRIDCELL half of
# `Summary_carbonflux` (CNVegCarbonFluxType.F90:5322-5505).
#
# This half was a one-line stub ("Column-level p2c aggregation is stubbed
# pending subgridAveMod port") for the entire life of the port — even though
# subgridAveMod *is* ported (p2c_1d_filter! / c2g). The consequence: EVERY
# column and gridcell carbon flux was dead. `gpp_col`, `ar_col`, `rr_col`,
# `er_col`, `nep_col`, `fire_closs_col`, `hrv_xsmrpool_to_atm_col`,
# `gru_conv_cflux_col`, `nbp_grc`, `nee_grc` and `landuseflux_grc` were
# allocated to NaN, zeroed once per step by SetValues, and never written
# again. That is exactly why the carbon half of the CN balance check could
# not run: it reads gpp_col and er_col, which were structurally 0.0.
#
# Verified against the live model before the fix: ar_patch = 1.5e-7 but
# ar_col = 0.0 — the p2c was not merely wrong, it was absent.
# ==========================================================================

"""
    cnveg_carbon_flux_summary_col!(cf, mask_bgc_soilc, bounds_col, bounds_grc,
        col, pch; soilbgc_hr_col, ...)

Column- and gridcell-level half of the CNVeg carbon-flux summary: the patch→column
aggregations (`gpp_col`, `ar_col`, `rr_col`, `npp_col`, `fire_closs_p2c_col`,
`hrv_xsmrpool_to_atm_col`, `gru_conv_cflux_col`), the vertically-integrated column
fire losses, the derived column fluxes (`sr_col`, `er_col`, `nep_col`,
`fire_closs_col`), and the gridcell aggregates (`nee_grc`, `landuseflux_grc`,
`nbp_grc`).

Must be called AFTER `soil_bgc_carbon_flux_summary!` (it consumes `hr_col`,
`cwdhr_col`, `lithr_col`) and after the patch-level `cnveg_carbon_flux_summary!`.

Ported from the column/gridcell section of `Summary_carbonflux` in
`CNVegCarbonFluxType.F90`.

!!! note "Annual flux dribblers are not ported"
    Fortran routes `hrv_xsmrpool_to_atm`, `dwt_conv_cflux` and `gru_conv_cflux`
    through `annual_flux_dribbler_gridcell` before they enter NEE/land-use flux,
    which spreads each year's pulse evenly over the following year. That subsystem
    is not ported (see `dyn_cons_biogeochem.jl`), so here the **dribbled flux is
    the instantaneous flux** and the dribbler's `amount_left` is identically zero.

    This is a *timing* difference in the NBP/NEE **diagnostic**, not a conservation
    hole: whatever leaves the column is counted in NBP in the same step it leaves,
    so the gridcell budget still closes exactly (the balance check's `amount_left`
    terms are zero on both the begin and end side, consistently). In every
    configuration validated here (no land-use change, no wood harvest, no crop)
    all three source fluxes are identically zero, so this is also bit-identical to
    Fortran. `c_balance_check!` asserts that (see `dribbler_flux_tol`).
"""
function cnveg_carbon_flux_summary_col!(
        cf::CNVegCarbonFluxData,
        mask_bgc_soilc::AbstractVector{Bool},
        bounds_col::UnitRange{Int},
        bounds_grc::UnitRange{Int},
        col::ColumnData,
        pch::PatchData;
        soilbgc_hr_col::AbstractVector{<:Real},
        soilbgc_cwdhr_col::AbstractVector{<:Real},
        soilbgc_lithr_col::AbstractVector{<:Real},
        soilbgc_decomp_cascade_ctransfer_col::AbstractMatrix{<:Real},
        product_closs_grc::AbstractVector{<:Real},
        nlevdecomp::Int,
        ndecomp_pools::Int,
        ndecomp_cascade_transitions::Int,
        dzsoi_decomp_vals::AbstractVector{<:Real},
        cascade_donor_pool::AbstractVector{Int},
        is_cwd::AbstractVector{Bool},
        is_litter::AbstractVector{Bool},
        use_crop::Bool=false,
        dribble_crophrv_xsmrpool_2atm::Bool=false)

    isempty(bounds_col) && return nothing

    # ---- p2c: patch → column (Fortran CNVegCarbonFluxType.F90:5328-5366) ----
    p2c_1d_filter!(cf.hrv_xsmrpool_to_atm_col, cf.hrv_xsmrpool_to_atm_patch,
                   mask_bgc_soilc, col, pch)
    if use_crop && dribble_crophrv_xsmrpool_2atm
        p2c_1d_filter!(cf.xsmrpool_to_atm_col, cf.xsmrpool_to_atm_patch,
                       mask_bgc_soilc, col, pch)
        c2g_unity!(cf.xsmrpool_to_atm_grc, cf.xsmrpool_to_atm_col,
                   col.gridcell, col.wtgcell, bounds_col, bounds_grc)
    end
    p2c_1d_filter!(cf.gru_conv_cflux_col, cf.gru_conv_cflux_patch, mask_bgc_soilc, col, pch)
    p2c_1d_filter!(cf.fire_closs_p2c_col, cf.fire_closs_patch,   mask_bgc_soilc, col, pch)
    p2c_1d_filter!(cf.npp_col,            cf.npp_patch,          mask_bgc_soilc, col, pch)
    p2c_1d_filter!(cf.rr_col,             cf.rr_patch,           mask_bgc_soilc, col, pch)
    p2c_1d_filter!(cf.ar_col,             cf.ar_patch,           mask_bgc_soilc, col, pch)
    p2c_1d_filter!(cf.gpp_col,            cf.gpp_patch,          mask_bgc_soilc, col, pch)

    # ---- vertically integrate column-level carbon fire losses ----
    if !isempty(cf.m_decomp_cpools_to_fire_vr_col)
        for l in 1:ndecomp_pools, j in 1:nlevdecomp, c in bounds_col
            mask_bgc_soilc[c] || continue
            cf.m_decomp_cpools_to_fire_col[c, l] =
                cf.m_decomp_cpools_to_fire_col[c, l] +
                cf.m_decomp_cpools_to_fire_vr_col[c, j, l] * dzsoi_decomp_vals[j]
        end
    end

    # ---- derived column fluxes ----
    for c in bounds_col
        mask_bgc_soilc[c] || continue

        cf.litfire_col[c] = 0.0
        cf.somfire_col[c] = 0.0
        cf.totfire_col[c] = cf.litfire_col[c] + cf.somfire_col[c]

        # carbon losses to fire, including patch losses
        fc_c = cf.fire_closs_p2c_col[c]
        if !isempty(cf.m_decomp_cpools_to_fire_col)
            for l in 1:ndecomp_pools
                fc_c += cf.m_decomp_cpools_to_fire_col[c, l]
            end
        end
        cf.fire_closs_col[c] = fc_c

        # total soil respiration, heterotrophic + root respiration (SR)
        cf.sr_col[c] = cf.rr_col[c] + soilbgc_hr_col[c]
        # total ecosystem respiration, autotrophic + heterotrophic (ER)
        cf.er_col[c] = cf.ar_col[c] + soilbgc_hr_col[c]
        # net ecosystem production (NEP) — excludes fire, landcover change, products
        cf.nep_col[c] = cf.gpp_col[c] - cf.er_col[c]
    end

    # ---- c2g: column → gridcell (c2l_scale_type='unity', l2g_scale_type='unity') ----
    if !isempty(bounds_grc)
        FT = eltype(cf.nep_col)
        ng = length(cf.nbp_grc)
        nep_grc        = fill!(similar(cf.nbp_grc, FT, ng), zero(FT))
        fire_closs_grc = fill!(similar(cf.nbp_grc, FT, ng), zero(FT))

        c2g_unity!(nep_grc, cf.nep_col, col.gridcell, col.wtgcell, bounds_col, bounds_grc)
        c2g_unity!(fire_closs_grc, cf.fire_closs_col, col.gridcell, col.wtgcell,
                   bounds_col, bounds_grc)
        c2g_unity!(cf.hrv_xsmrpool_to_atm_grc, cf.hrv_xsmrpool_to_atm_col,
                   col.gridcell, col.wtgcell, bounds_col, bounds_grc)
        c2g_unity!(cf.gru_conv_cflux_grc, cf.gru_conv_cflux_col,
                   col.gridcell, col.wtgcell, bounds_col, bounds_grc)

        # Annual flux dribblers are NOT ported: dribbled == instantaneous. See the
        # docstring note — this is a diagnostic-timing difference, not a budget hole,
        # and all three fluxes are identically zero in every validated config.
        for g in bounds_grc
            cf.dwt_conv_cflux_dribbled_grc[g] = cf.dwt_conv_cflux_grc[g]
            cf.gru_conv_cflux_dribbled_grc[g] = cf.gru_conv_cflux_grc[g]
        end

        for g in bounds_grc
            # net ecosystem exchange, incl. fire + hrv_xsmrpool; positive = source
            cf.nee_grc[g] = -nep_grc[g] + fire_closs_grc[g] +
                            cf.hrv_xsmrpool_to_atm_grc[g]

            cf.landuseflux_grc[g] = cf.dwt_conv_cflux_dribbled_grc[g] +
                                    cf.gru_conv_cflux_dribbled_grc[g] +
                                    product_closs_grc[g]

            # net biome production; positive = sink
            cf.nbp_grc[g] = -cf.nee_grc[g] - cf.landuseflux_grc[g]
            if dribble_crophrv_xsmrpool_2atm
                cf.nbp_grc[g] -= cf.xsmrpool_to_atm_grc[g]
            end
        end
    end

    # ---- coarse woody debris + litter C loss (diagnostics) ----
    for c in bounds_col
        mask_bgc_soilc[c] || continue
        cf.cwdc_loss_col[c]    = soilbgc_cwdhr_col[c]
        cf.litterc_loss_col[c] = soilbgc_lithr_col[c]
    end
    if !isempty(cf.m_decomp_cpools_to_fire_col)
        for l in 1:ndecomp_pools
            is_cwd[l] || continue
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                cf.cwdc_loss_col[c] += cf.m_decomp_cpools_to_fire_col[c, l]
            end
        end
        for l in 1:ndecomp_pools
            is_litter[l] || continue
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                cf.litterc_loss_col[c] += cf.m_decomp_cpools_to_fire_col[c, l]
            end
        end
    end
    for k in 1:ndecomp_cascade_transitions
        donor = cascade_donor_pool[k]
        (donor >= 1 && donor <= ndecomp_pools) || continue
        if is_cwd[donor]
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                cf.cwdc_loss_col[c] += soilbgc_decomp_cascade_ctransfer_col[c, k]
            end
        end
        if is_litter[donor]
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                cf.litterc_loss_col[c] += soilbgc_decomp_cascade_ctransfer_col[c, k]
            end
        end
    end

    return nothing
end

# ==========================================================================
# Stub functions — history, restart, dynamic patch adjustments
# ==========================================================================

function cnveg_carbon_flux_init_history!(cf::CNVegCarbonFluxData{FT}, args...) where {FT}
    nothing
end

function cnveg_carbon_flux_restart!(cf::CNVegCarbonFluxData{FT}, args...) where {FT}
    nothing
end
