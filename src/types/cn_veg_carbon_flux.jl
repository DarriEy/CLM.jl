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
Base.@kwdef mutable struct CNVegCarbonFluxData{FT<:AbstractFloat}
    dribble_crophrv_xsmrpool_2atm::Bool = false

    # --- Gap mortality fluxes (gC/m2/s) patch-level ---
    m_leafc_to_litter_patch                   ::Vector{FT} = Float64[]
    m_leafc_storage_to_litter_patch           ::Vector{FT} = Float64[]
    m_leafc_xfer_to_litter_patch              ::Vector{FT} = Float64[]
    m_frootc_to_litter_patch                  ::Vector{FT} = Float64[]
    m_frootc_storage_to_litter_patch          ::Vector{FT} = Float64[]
    m_frootc_xfer_to_litter_patch             ::Vector{FT} = Float64[]
    m_livestemc_to_litter_patch               ::Vector{FT} = Float64[]
    m_livestemc_storage_to_litter_patch       ::Vector{FT} = Float64[]
    m_livestemc_xfer_to_litter_patch          ::Vector{FT} = Float64[]
    m_deadstemc_to_litter_patch               ::Vector{FT} = Float64[]
    m_deadstemc_storage_to_litter_patch       ::Vector{FT} = Float64[]
    m_deadstemc_xfer_to_litter_patch          ::Vector{FT} = Float64[]
    m_livecrootc_to_litter_patch              ::Vector{FT} = Float64[]
    m_livecrootc_storage_to_litter_patch      ::Vector{FT} = Float64[]
    m_livecrootc_xfer_to_litter_patch         ::Vector{FT} = Float64[]
    m_deadcrootc_to_litter_patch              ::Vector{FT} = Float64[]
    m_deadcrootc_storage_to_litter_patch      ::Vector{FT} = Float64[]
    m_deadcrootc_xfer_to_litter_patch         ::Vector{FT} = Float64[]
    m_gresp_storage_to_litter_patch           ::Vector{FT} = Float64[]
    m_gresp_xfer_to_litter_patch              ::Vector{FT} = Float64[]

    # --- Harvest mortality fluxes (gC/m2/s) patch-level ---
    hrv_leafc_to_litter_patch                 ::Vector{FT} = Float64[]
    hrv_leafc_storage_to_litter_patch         ::Vector{FT} = Float64[]
    hrv_leafc_xfer_to_litter_patch            ::Vector{FT} = Float64[]
    hrv_frootc_to_litter_patch                ::Vector{FT} = Float64[]
    hrv_frootc_storage_to_litter_patch        ::Vector{FT} = Float64[]
    hrv_frootc_xfer_to_litter_patch           ::Vector{FT} = Float64[]
    hrv_livestemc_to_litter_patch             ::Vector{FT} = Float64[]
    hrv_livestemc_storage_to_litter_patch     ::Vector{FT} = Float64[]
    hrv_livestemc_xfer_to_litter_patch        ::Vector{FT} = Float64[]
    hrv_deadstemc_storage_to_litter_patch     ::Vector{FT} = Float64[]
    hrv_deadstemc_xfer_to_litter_patch        ::Vector{FT} = Float64[]
    hrv_livecrootc_to_litter_patch            ::Vector{FT} = Float64[]
    hrv_livecrootc_storage_to_litter_patch    ::Vector{FT} = Float64[]
    hrv_livecrootc_xfer_to_litter_patch       ::Vector{FT} = Float64[]
    hrv_deadcrootc_to_litter_patch            ::Vector{FT} = Float64[]
    hrv_deadcrootc_storage_to_litter_patch    ::Vector{FT} = Float64[]
    hrv_deadcrootc_xfer_to_litter_patch       ::Vector{FT} = Float64[]
    hrv_gresp_storage_to_litter_patch         ::Vector{FT} = Float64[]
    hrv_gresp_xfer_to_litter_patch            ::Vector{FT} = Float64[]
    hrv_xsmrpool_to_atm_patch                 ::Vector{FT} = Float64[]

    # --- Fire fluxes (gC/m2/s) patch-level ---
    m_leafc_to_fire_patch                     ::Vector{FT} = Float64[]
    m_leafc_storage_to_fire_patch             ::Vector{FT} = Float64[]
    m_leafc_xfer_to_fire_patch                ::Vector{FT} = Float64[]
    m_livestemc_to_fire_patch                 ::Vector{FT} = Float64[]
    m_livestemc_storage_to_fire_patch         ::Vector{FT} = Float64[]
    m_livestemc_xfer_to_fire_patch            ::Vector{FT} = Float64[]
    m_deadstemc_to_fire_patch                 ::Vector{FT} = Float64[]
    m_deadstemc_storage_to_fire_patch         ::Vector{FT} = Float64[]
    m_deadstemc_xfer_to_fire_patch            ::Vector{FT} = Float64[]
    m_frootc_to_fire_patch                    ::Vector{FT} = Float64[]
    m_frootc_storage_to_fire_patch            ::Vector{FT} = Float64[]
    m_frootc_xfer_to_fire_patch               ::Vector{FT} = Float64[]
    m_livecrootc_to_fire_patch                ::Vector{FT} = Float64[]
    m_livecrootc_storage_to_fire_patch        ::Vector{FT} = Float64[]
    m_livecrootc_xfer_to_fire_patch           ::Vector{FT} = Float64[]
    m_deadcrootc_to_fire_patch                ::Vector{FT} = Float64[]
    m_deadcrootc_storage_to_fire_patch        ::Vector{FT} = Float64[]
    m_deadcrootc_xfer_to_fire_patch           ::Vector{FT} = Float64[]
    m_gresp_storage_to_fire_patch             ::Vector{FT} = Float64[]
    m_gresp_xfer_to_fire_patch                ::Vector{FT} = Float64[]

    # --- Fire-to-litter fluxes (gC/m2/s) patch-level ---
    m_leafc_to_litter_fire_patch              ::Vector{FT} = Float64[]
    m_leafc_storage_to_litter_fire_patch      ::Vector{FT} = Float64[]
    m_leafc_xfer_to_litter_fire_patch         ::Vector{FT} = Float64[]
    m_livestemc_to_litter_fire_patch          ::Vector{FT} = Float64[]
    m_livestemc_storage_to_litter_fire_patch  ::Vector{FT} = Float64[]
    m_livestemc_xfer_to_litter_fire_patch     ::Vector{FT} = Float64[]
    m_livestemc_to_deadstemc_fire_patch       ::Vector{FT} = Float64[]
    m_deadstemc_to_litter_fire_patch          ::Vector{FT} = Float64[]
    m_deadstemc_storage_to_litter_fire_patch  ::Vector{FT} = Float64[]
    m_deadstemc_xfer_to_litter_fire_patch     ::Vector{FT} = Float64[]
    m_frootc_to_litter_fire_patch             ::Vector{FT} = Float64[]
    m_frootc_storage_to_litter_fire_patch     ::Vector{FT} = Float64[]
    m_frootc_xfer_to_litter_fire_patch        ::Vector{FT} = Float64[]
    m_livecrootc_to_litter_fire_patch         ::Vector{FT} = Float64[]
    m_livecrootc_storage_to_litter_fire_patch ::Vector{FT} = Float64[]
    m_livecrootc_xfer_to_litter_fire_patch    ::Vector{FT} = Float64[]
    m_livecrootc_to_deadcrootc_fire_patch     ::Vector{FT} = Float64[]
    m_deadcrootc_to_litter_fire_patch         ::Vector{FT} = Float64[]
    m_deadcrootc_storage_to_litter_fire_patch ::Vector{FT} = Float64[]
    m_deadcrootc_xfer_to_litter_fire_patch    ::Vector{FT} = Float64[]
    m_gresp_storage_to_litter_fire_patch      ::Vector{FT} = Float64[]
    m_gresp_xfer_to_litter_fire_patch         ::Vector{FT} = Float64[]

    # --- Phenology fluxes (gC/m2/s) ---
    reproductivec_xfer_to_reproductivec_patch ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    leafc_xfer_to_leafc_patch                 ::Vector{FT} = Float64[]
    frootc_xfer_to_frootc_patch               ::Vector{FT} = Float64[]
    livestemc_xfer_to_livestemc_patch         ::Vector{FT} = Float64[]
    deadstemc_xfer_to_deadstemc_patch         ::Vector{FT} = Float64[]
    livecrootc_xfer_to_livecrootc_patch       ::Vector{FT} = Float64[]
    deadcrootc_xfer_to_deadcrootc_patch       ::Vector{FT} = Float64[]

    # --- Leaf/froot litterfall and crop fluxes (gC/m2/s) ---
    leafc_to_litter_patch                     ::Vector{FT} = Float64[]
    leafc_to_litter_fun_patch                 ::Vector{FT} = Float64[]
    frootc_to_litter_patch                    ::Vector{FT} = Float64[]
    livestemc_to_litter_patch                 ::Vector{FT} = Float64[]
    repr_grainc_to_food_patch                 ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    repr_grainc_to_food_perharv_patch         ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    repr_grainc_to_food_thisyr_patch          ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    repr_structurec_to_cropprod_patch         ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    repr_structurec_to_litter_patch           ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    leafc_to_biofuelc_patch                   ::Vector{FT} = Float64[]
    livestemc_to_biofuelc_patch               ::Vector{FT} = Float64[]
    leafc_to_removedresiduec_patch            ::Vector{FT} = Float64[]
    livestemc_to_removedresiduec_patch        ::Vector{FT} = Float64[]
    repr_grainc_to_seed_patch                 ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    repr_grainc_to_seed_perharv_patch         ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    repr_grainc_to_seed_thisyr_patch          ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)

    # --- Maintenance respiration (gC/m2/s) ---
    cpool_to_resp_patch                       ::Vector{FT} = Float64[]
    cpool_to_leafc_resp_patch                 ::Vector{FT} = Float64[]
    cpool_to_leafc_storage_resp_patch         ::Vector{FT} = Float64[]
    cpool_to_frootc_resp_patch                ::Vector{FT} = Float64[]
    cpool_to_frootc_storage_resp_patch        ::Vector{FT} = Float64[]
    cpool_to_livecrootc_resp_patch            ::Vector{FT} = Float64[]
    cpool_to_livecrootc_storage_resp_patch    ::Vector{FT} = Float64[]
    cpool_to_livestemc_resp_patch             ::Vector{FT} = Float64[]
    cpool_to_livestemc_storage_resp_patch     ::Vector{FT} = Float64[]
    leaf_mr_patch                             ::Vector{FT} = Float64[]
    froot_mr_patch                            ::Vector{FT} = Float64[]
    livestem_mr_patch                         ::Vector{FT} = Float64[]
    livecroot_mr_patch                        ::Vector{FT} = Float64[]
    reproductive_mr_patch                     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    leaf_curmr_patch                          ::Vector{FT} = Float64[]
    froot_curmr_patch                         ::Vector{FT} = Float64[]
    livestem_curmr_patch                      ::Vector{FT} = Float64[]
    livecroot_curmr_patch                     ::Vector{FT} = Float64[]
    reproductive_curmr_patch                  ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    leaf_xsmr_patch                           ::Vector{FT} = Float64[]
    froot_xsmr_patch                          ::Vector{FT} = Float64[]
    livestem_xsmr_patch                       ::Vector{FT} = Float64[]
    livecroot_xsmr_patch                      ::Vector{FT} = Float64[]
    reproductive_xsmr_patch                   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)

    # --- Photosynthesis (gC/m2/s) ---
    psnsun_to_cpool_patch                     ::Vector{FT} = Float64[]
    psnshade_to_cpool_patch                   ::Vector{FT} = Float64[]

    # --- Allocation from cpool (gC/m2/s) ---
    cpool_to_xsmrpool_patch                   ::Vector{FT} = Float64[]
    cpool_to_reproductivec_patch              ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    cpool_to_reproductivec_storage_patch      ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    cpool_to_leafc_patch                      ::Vector{FT} = Float64[]
    cpool_to_leafc_storage_patch              ::Vector{FT} = Float64[]
    cpool_to_frootc_patch                     ::Vector{FT} = Float64[]
    cpool_to_frootc_storage_patch             ::Vector{FT} = Float64[]
    cpool_to_livestemc_patch                  ::Vector{FT} = Float64[]
    cpool_to_livestemc_storage_patch          ::Vector{FT} = Float64[]
    cpool_to_deadstemc_patch                  ::Vector{FT} = Float64[]
    cpool_to_deadstemc_storage_patch          ::Vector{FT} = Float64[]
    cpool_to_livecrootc_patch                 ::Vector{FT} = Float64[]
    cpool_to_livecrootc_storage_patch         ::Vector{FT} = Float64[]
    cpool_to_deadcrootc_patch                 ::Vector{FT} = Float64[]
    cpool_to_deadcrootc_storage_patch         ::Vector{FT} = Float64[]
    cpool_to_gresp_storage_patch              ::Vector{FT} = Float64[]

    # --- Growth respiration (gC/m2/s) ---
    xsmrpool_to_atm_patch                     ::Vector{FT} = Float64[]
    xsmrpool_to_atm_col                       ::Vector{FT} = Float64[]
    xsmrpool_to_atm_grc                       ::Vector{FT} = Float64[]
    cpool_leaf_gr_patch                       ::Vector{FT} = Float64[]
    cpool_leaf_storage_gr_patch               ::Vector{FT} = Float64[]
    transfer_leaf_gr_patch                    ::Vector{FT} = Float64[]
    cpool_froot_gr_patch                      ::Vector{FT} = Float64[]
    cpool_froot_storage_gr_patch              ::Vector{FT} = Float64[]
    transfer_froot_gr_patch                   ::Vector{FT} = Float64[]
    cpool_livestem_gr_patch                   ::Vector{FT} = Float64[]
    cpool_livestem_storage_gr_patch           ::Vector{FT} = Float64[]
    transfer_livestem_gr_patch                ::Vector{FT} = Float64[]
    cpool_deadstem_gr_patch                   ::Vector{FT} = Float64[]
    cpool_deadstem_storage_gr_patch           ::Vector{FT} = Float64[]
    transfer_deadstem_gr_patch                ::Vector{FT} = Float64[]
    cpool_livecroot_gr_patch                  ::Vector{FT} = Float64[]
    cpool_livecroot_storage_gr_patch          ::Vector{FT} = Float64[]
    transfer_livecroot_gr_patch               ::Vector{FT} = Float64[]
    cpool_deadcroot_gr_patch                  ::Vector{FT} = Float64[]
    cpool_deadcroot_storage_gr_patch          ::Vector{FT} = Float64[]
    transfer_deadcroot_gr_patch               ::Vector{FT} = Float64[]
    cpool_reproductive_gr_patch               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    cpool_reproductive_storage_gr_patch       ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    transfer_reproductive_gr_patch            ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)

    # --- Storage to transfer turnover (gC/m2/s) ---
    reproductivec_storage_to_xfer_patch       ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    leafc_storage_to_xfer_patch               ::Vector{FT} = Float64[]
    frootc_storage_to_xfer_patch              ::Vector{FT} = Float64[]
    livestemc_storage_to_xfer_patch           ::Vector{FT} = Float64[]
    deadstemc_storage_to_xfer_patch           ::Vector{FT} = Float64[]
    livecrootc_storage_to_xfer_patch          ::Vector{FT} = Float64[]
    deadcrootc_storage_to_xfer_patch          ::Vector{FT} = Float64[]
    gresp_storage_to_xfer_patch               ::Vector{FT} = Float64[]

    # --- Livewood to deadwood turnover (gC/m2/s) ---
    livestemc_to_deadstemc_patch              ::Vector{FT} = Float64[]
    livecrootc_to_deadcrootc_patch            ::Vector{FT} = Float64[]

    # --- Column-level decomposition fluxes (gC/m3/s) ---
    phenology_c_to_litr_c_col                 ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    gap_mortality_c_to_litr_c_col             ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    gap_mortality_c_to_cwdc_col               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    fire_mortality_c_to_cwdc_col              ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    harvest_c_to_litr_c_col                   ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    harvest_c_to_cwdc_col                     ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    crop_harvestc_to_cropprodc_patch          ::Vector{FT} = Float64[]
    crop_harvestc_to_cropprodc_col            ::Vector{FT} = Float64[]
    m_decomp_cpools_to_fire_vr_col            ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    m_decomp_cpools_to_fire_col               ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    m_c_to_litr_fire_col                      ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)

    # --- Dynamic landcover fluxes ---
    dwt_seedc_to_leaf_patch                   ::Vector{FT} = Float64[]
    dwt_seedc_to_leaf_grc                     ::Vector{FT} = Float64[]
    dwt_seedc_to_deadstem_patch               ::Vector{FT} = Float64[]
    dwt_seedc_to_deadstem_grc                 ::Vector{FT} = Float64[]
    dwt_conv_cflux_patch                      ::Vector{FT} = Float64[]
    dwt_conv_cflux_grc                        ::Vector{FT} = Float64[]
    dwt_conv_cflux_dribbled_grc               ::Vector{FT} = Float64[]
    dwt_wood_productc_gain_patch              ::Vector{FT} = Float64[]
    dwt_crop_productc_gain_patch              ::Vector{FT} = Float64[]
    dwt_slash_cflux_patch                     ::Vector{FT} = Float64[]
    dwt_slash_cflux_grc                       ::Vector{FT} = Float64[]
    dwt_frootc_to_litr_c_col                  ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    dwt_livecrootc_to_cwdc_col                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    dwt_deadcrootc_to_cwdc_col                ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)

    # --- Gross unrepresented landcover fluxes ---
    gru_leafc_to_litter_patch                 ::Vector{FT} = Float64[]
    gru_leafc_storage_to_atm_patch            ::Vector{FT} = Float64[]
    gru_leafc_xfer_to_atm_patch               ::Vector{FT} = Float64[]
    gru_frootc_to_litter_patch                ::Vector{FT} = Float64[]
    gru_frootc_storage_to_atm_patch           ::Vector{FT} = Float64[]
    gru_frootc_xfer_to_atm_patch              ::Vector{FT} = Float64[]
    gru_livestemc_to_atm_patch                ::Vector{FT} = Float64[]
    gru_livestemc_storage_to_atm_patch        ::Vector{FT} = Float64[]
    gru_livestemc_xfer_to_atm_patch           ::Vector{FT} = Float64[]
    gru_deadstemc_to_atm_patch                ::Vector{FT} = Float64[]
    gru_deadstemc_storage_to_atm_patch        ::Vector{FT} = Float64[]
    gru_deadstemc_xfer_to_atm_patch           ::Vector{FT} = Float64[]
    gru_livecrootc_to_litter_patch            ::Vector{FT} = Float64[]
    gru_livecrootc_storage_to_atm_patch       ::Vector{FT} = Float64[]
    gru_livecrootc_xfer_to_atm_patch          ::Vector{FT} = Float64[]
    gru_deadcrootc_to_litter_patch            ::Vector{FT} = Float64[]
    gru_deadcrootc_storage_to_atm_patch       ::Vector{FT} = Float64[]
    gru_deadcrootc_xfer_to_atm_patch          ::Vector{FT} = Float64[]
    gru_gresp_storage_to_atm_patch            ::Vector{FT} = Float64[]
    gru_gresp_xfer_to_atm_patch               ::Vector{FT} = Float64[]
    gru_xsmrpool_to_atm_patch                 ::Vector{FT} = Float64[]
    gru_conv_cflux_patch                      ::Vector{FT} = Float64[]
    gru_conv_cflux_col                        ::Vector{FT} = Float64[]
    gru_conv_cflux_grc                        ::Vector{FT} = Float64[]
    gru_conv_cflux_dribbled_grc               ::Vector{FT} = Float64[]
    gru_wood_productc_gain_patch              ::Vector{FT} = Float64[]
    gru_wood_productc_gain_col                ::Vector{FT} = Float64[]
    gru_slash_cflux_patch                     ::Vector{FT} = Float64[]
    gru_c_to_litr_c_col                       ::Array{FT,3} = Array{Float64}(undef, 0, 0, 0)
    gru_c_to_cwdc_col                         ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)

    # --- Crop fluxes ---
    crop_seedc_to_leaf_patch                  ::Vector{FT} = Float64[]

    # --- Summary/diagnostic fluxes (gC/m2/s) ---
    gpp_before_downreg_patch                  ::Vector{FT} = Float64[]
    current_gr_patch                          ::Vector{FT} = Float64[]
    transfer_gr_patch                         ::Vector{FT} = Float64[]
    storage_gr_patch                          ::Vector{FT} = Float64[]
    plant_calloc_patch                        ::Vector{FT} = Float64[]
    excess_cflux_patch                        ::Vector{FT} = Float64[]
    prev_leafc_to_litter_patch                ::Vector{FT} = Float64[]
    prev_frootc_to_litter_patch               ::Vector{FT} = Float64[]
    availc_patch                              ::Vector{FT} = Float64[]
    xsmrpool_recover_patch                    ::Vector{FT} = Float64[]
    xsmrpool_c13ratio_patch                   ::Vector{FT} = Float64[]
    cwdc_loss_col                             ::Vector{FT} = Float64[]
    litterc_loss_col                          ::Vector{FT} = Float64[]
    frootc_alloc_patch                        ::Vector{FT} = Float64[]
    frootc_loss_patch                         ::Vector{FT} = Float64[]
    leafc_alloc_patch                         ::Vector{FT} = Float64[]
    leafc_loss_patch                          ::Vector{FT} = Float64[]
    woodc_alloc_patch                         ::Vector{FT} = Float64[]
    woodc_loss_patch                          ::Vector{FT} = Float64[]
    gpp_patch                                 ::Vector{FT} = Float64[]
    gpp_col                                   ::Vector{FT} = Float64[]
    rr_patch                                  ::Vector{FT} = Float64[]
    rr_col                                    ::Vector{FT} = Float64[]
    mr_patch                                  ::Vector{FT} = Float64[]
    gr_patch                                  ::Vector{FT} = Float64[]
    ar_patch                                  ::Vector{FT} = Float64[]
    ar_col                                    ::Vector{FT} = Float64[]
    npp_patch                                 ::Vector{FT} = Float64[]
    npp_col                                   ::Vector{FT} = Float64[]
    agnpp_patch                               ::Vector{FT} = Float64[]
    bgnpp_patch                               ::Vector{FT} = Float64[]
    litfall_patch                             ::Vector{FT} = Float64[]
    wood_harvestc_patch                       ::Vector{FT} = Float64[]
    wood_harvestc_col                         ::Vector{FT} = Float64[]
    slash_harvestc_patch                      ::Vector{FT} = Float64[]
    cinputs_patch                             ::Vector{FT} = Float64[]
    coutputs_patch                            ::Vector{FT} = Float64[]
    sr_col                                    ::Vector{FT} = Float64[]
    er_col                                    ::Vector{FT} = Float64[]
    litfire_col                               ::Vector{FT} = Float64[]
    somfire_col                               ::Vector{FT} = Float64[]
    totfire_col                               ::Vector{FT} = Float64[]
    hrv_xsmrpool_to_atm_col                  ::Vector{FT} = Float64[]

    # --- Fire summary ---
    fire_closs_patch                          ::Vector{FT} = Float64[]
    fire_closs_p2c_col                        ::Vector{FT} = Float64[]
    fire_closs_col                            ::Vector{FT} = Float64[]

    # --- Annual sums ---
    tempsum_litfall_patch                     ::Vector{FT} = Float64[]
    annsum_litfall_patch                      ::Vector{FT} = Float64[]
    tempsum_npp_patch                         ::Vector{FT} = Float64[]
    annsum_npp_patch                          ::Vector{FT} = Float64[]
    annsum_npp_col                            ::Vector{FT} = Float64[]
    lag_npp_col                               ::Vector{FT} = Float64[]

    # --- Summary C fluxes ---
    nep_col                                   ::Vector{FT} = Float64[]
    nbp_grc                                   ::Vector{FT} = Float64[]
    nee_grc                                   ::Vector{FT} = Float64[]
    landuseflux_grc                           ::Vector{FT} = Float64[]

    # --- FUN fluxes (gC/m2/s) patch-level ---
    npp_Nactive_patch                         ::Vector{FT} = Float64[]
    npp_burnedoff_patch                       ::Vector{FT} = Float64[]
    npp_Nnonmyc_patch                         ::Vector{FT} = Float64[]
    npp_Nam_patch                             ::Vector{FT} = Float64[]
    npp_Necm_patch                            ::Vector{FT} = Float64[]
    npp_Nactive_no3_patch                     ::Vector{FT} = Float64[]
    npp_Nactive_nh4_patch                     ::Vector{FT} = Float64[]
    npp_Nnonmyc_no3_patch                     ::Vector{FT} = Float64[]
    npp_Nnonmyc_nh4_patch                     ::Vector{FT} = Float64[]
    npp_Nam_no3_patch                         ::Vector{FT} = Float64[]
    npp_Nam_nh4_patch                         ::Vector{FT} = Float64[]
    npp_Necm_no3_patch                        ::Vector{FT} = Float64[]
    npp_Necm_nh4_patch                        ::Vector{FT} = Float64[]
    npp_Nfix_patch                            ::Vector{FT} = Float64[]
    npp_Nretrans_patch                        ::Vector{FT} = Float64[]
    npp_Nuptake_patch                         ::Vector{FT} = Float64[]
    npp_growth_patch                          ::Vector{FT} = Float64[]
    leafc_change_patch                        ::Vector{FT} = Float64[]
    soilc_change_patch                        ::Vector{FT} = Float64[]

    # --- Matrix CN arrays ---
    matrix_Cinput_patch                       ::Vector{FT} = Float64[]
    matrix_C13input_patch                     ::Vector{FT} = Float64[]
    matrix_C14input_patch                     ::Vector{FT} = Float64[]
    matrix_alloc_patch                        ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    matrix_phtransfer_patch                   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    matrix_phturnover_patch                   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    matrix_phtransfer_doner_patch             ::Vector{Int} = Int[]
    matrix_phtransfer_receiver_patch          ::Vector{Int} = Int[]
    matrix_gmtransfer_patch                   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    matrix_gmturnover_patch                   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    matrix_gmtransfer_doner_patch             ::Vector{Int} = Int[]
    matrix_gmtransfer_receiver_patch          ::Vector{Int} = Int[]
    matrix_fitransfer_patch                   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    matrix_fiturnover_patch                   ::Matrix{FT} = Matrix{Float64}(undef, 0, 0)
    matrix_fitransfer_doner_patch             ::Vector{Int} = Int[]
    matrix_fitransfer_receiver_patch          ::Vector{Int} = Int[]

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

# ==========================================================================
# Init function — allocates all arrays
# Ported from: InitAllocate subroutine
# ==========================================================================

"""
    cnveg_carbon_flux_init!(cf, np, nc, ng; use_matrixcn=false, nrepr=NREPR,
                            nlevdecomp_full=1, ndecomp_pools=1, mxharvests=MXHARVESTS)

Allocate and initialize all carbon flux arrays.
"""
function cnveg_carbon_flux_init!(cf::CNVegCarbonFluxData, np::Int, nc::Int, ng::Int;
                                  use_matrixcn::Bool=false,
                                  nrepr::Int=NREPR,
                                  nlevdecomp_full::Int=1,
                                  ndecomp_pools::Int=1,
                                  mxharvests::Int=MXHARVESTS,
                                  nvegcpool::Int=NVEGPOOL_NATVEG)
    # Helper for NaN-filled arrays
    nanvec(n) = fill(NaN, n)
    nanmat(r, c) = fill(NaN, r, c)
    nan3d(a, b, c) = fill(NaN, a, b, c)

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
    cf.hrv_xsmrpool_to_atm_patch            = zeros(np)

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
    cf.xsmrpool_to_atm_patch               = zeros(np)
    cf.xsmrpool_to_atm_col                  = zeros(nc)
    cf.xsmrpool_to_atm_grc                  = zeros(ng)
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
    cf.hrv_xsmrpool_to_atm_col             = zeros(nc)
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
                                       mask_patch::BitVector, value_patch::Float64,
                                       mask_col::BitVector, value_column::Float64;
                                       use_matrixcn::Bool=false,
                                       use_crop::Bool=false,
                                       nrepr::Int=NREPR,
                                       nlevdecomp_full::Int=1,
                                       ndecomp_pools::Int=1,
                                       nvegcpool::Int=NVEGPOOL_NATVEG)
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
                                      ndecomp_pools::Int=1)
    for g in bounds_grc
        cf.dwt_seedc_to_leaf_grc[g]     = 0.0
        cf.dwt_seedc_to_deadstem_grc[g] = 0.0
        cf.dwt_conv_cflux_grc[g]        = 0.0
        cf.dwt_slash_cflux_grc[g]       = 0.0
    end

    for j in 1:nlevdecomp_full
        for c in bounds_col
            for k in 1:ndecomp_pools
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
function cnveg_carbon_flux_summary!(cf::CNVegCarbonFluxData,
                                     mask_patch::BitVector,
                                     bounds_patch::UnitRange{Int};
                                     use_crop::Bool=false,
                                     use_fun::Bool=false,
                                     carbon_resp_opt::Int=0,
                                     patch_itype::Union{Vector{Int},Nothing}=nothing,
                                     npcropmin::Int=0,
                                     nrepr::Int=NREPR)
    for p in bounds_patch
        mask_patch[p] || continue

        # --- Maintenance respiration (MR) ---
        cf.mr_patch[p] = cf.leaf_mr_patch[p] +
                         cf.froot_mr_patch[p] +
                         cf.livestem_mr_patch[p] +
                         cf.livecroot_mr_patch[p]

        if carbon_resp_opt == 1
            cf.mr_patch[p] += cf.cpool_to_resp_patch[p]
        end

        if use_crop && patch_itype !== nothing && patch_itype[p] >= npcropmin
            for k in 1:nrepr
                cf.mr_patch[p] += cf.reproductive_mr_patch[p, k]
            end
        end

        # --- Current growth respiration ---
        cf.current_gr_patch[p] = cf.cpool_leaf_gr_patch[p] +
                                  cf.cpool_froot_gr_patch[p] +
                                  cf.cpool_livestem_gr_patch[p] +
                                  cf.cpool_deadstem_gr_patch[p] +
                                  cf.cpool_livecroot_gr_patch[p] +
                                  cf.cpool_deadcroot_gr_patch[p]

        # --- Transfer growth respiration ---
        cf.transfer_gr_patch[p] = cf.transfer_leaf_gr_patch[p] +
                                   cf.transfer_froot_gr_patch[p] +
                                   cf.transfer_livestem_gr_patch[p] +
                                   cf.transfer_deadstem_gr_patch[p] +
                                   cf.transfer_livecroot_gr_patch[p] +
                                   cf.transfer_deadcroot_gr_patch[p]

        # --- Storage growth respiration ---
        cf.storage_gr_patch[p] = cf.cpool_leaf_storage_gr_patch[p] +
                                  cf.cpool_froot_storage_gr_patch[p] +
                                  cf.cpool_livestem_storage_gr_patch[p] +
                                  cf.cpool_deadstem_storage_gr_patch[p] +
                                  cf.cpool_livecroot_storage_gr_patch[p] +
                                  cf.cpool_deadcroot_storage_gr_patch[p]

        if use_crop && patch_itype !== nothing && patch_itype[p] >= npcropmin
            for k in 1:nrepr
                cf.current_gr_patch[p]  += cf.cpool_reproductive_gr_patch[p, k]
                cf.transfer_gr_patch[p] += cf.transfer_reproductive_gr_patch[p, k]
                cf.storage_gr_patch[p]  += cf.cpool_reproductive_storage_gr_patch[p, k]
            end
        end

        # --- Total growth respiration (GR) ---
        cf.gr_patch[p] = cf.current_gr_patch[p] +
                          cf.transfer_gr_patch[p] +
                          cf.storage_gr_patch[p]

        # --- Autotrophic respiration (AR) ---
        cf.ar_patch[p] = cf.mr_patch[p] + cf.gr_patch[p]

        if use_fun
            cf.ar_patch[p] += cf.soilc_change_patch[p]
        end

        if use_crop && patch_itype !== nothing && patch_itype[p] >= npcropmin
            cf.ar_patch[p] += cf.xsmrpool_to_atm_patch[p]
        end

        # --- Gross primary production (GPP) ---
        cf.gpp_patch[p] = cf.psnsun_to_cpool_patch[p] +
                           cf.psnshade_to_cpool_patch[p]

        # --- Net primary production (NPP) ---
        cf.npp_patch[p] = cf.gpp_patch[p] - cf.ar_patch[p]

        # --- Root respiration (RR) ---
        cf.rr_patch[p] = cf.froot_mr_patch[p] +
                          cf.cpool_froot_gr_patch[p] +
                          cf.cpool_livecroot_gr_patch[p] +
                          cf.cpool_deadcroot_gr_patch[p] +
                          cf.transfer_froot_gr_patch[p] +
                          cf.transfer_livecroot_gr_patch[p] +
                          cf.transfer_deadcroot_gr_patch[p] +
                          cf.cpool_froot_storage_gr_patch[p] +
                          cf.cpool_livecroot_storage_gr_patch[p] +
                          cf.cpool_deadcroot_storage_gr_patch[p]

        # --- Aboveground NPP (AGNPP) ---
        cf.agnpp_patch[p] = cf.cpool_to_leafc_patch[p] +
                             cf.cpool_to_leafc_storage_patch[p] +
                             cf.leafc_xfer_to_leafc_patch[p] +
                             cf.cpool_to_livestemc_patch[p] +
                             cf.cpool_to_livestemc_storage_patch[p] +
                             cf.livestemc_xfer_to_livestemc_patch[p] +
                             cf.cpool_to_deadstemc_patch[p] +
                             cf.cpool_to_deadstemc_storage_patch[p] +
                             cf.deadstemc_xfer_to_deadstemc_patch[p]

        # --- Belowground NPP (BGNPP) ---
        cf.bgnpp_patch[p] = cf.cpool_to_frootc_patch[p] +
                              cf.cpool_to_frootc_storage_patch[p] +
                              cf.frootc_xfer_to_frootc_patch[p] +
                              cf.cpool_to_livecrootc_patch[p] +
                              cf.cpool_to_livecrootc_storage_patch[p] +
                              cf.livecrootc_xfer_to_livecrootc_patch[p] +
                              cf.cpool_to_deadcrootc_patch[p] +
                              cf.cpool_to_deadcrootc_storage_patch[p] +
                              cf.deadcrootc_xfer_to_deadcrootc_patch[p]

        # --- Litterfall (LITFALL) ---
        cf.litfall_patch[p] = cf.leafc_to_litter_patch[p] +
                               cf.frootc_to_litter_patch[p] +
                               cf.m_leafc_to_litter_patch[p] +
                               cf.m_frootc_to_litter_patch[p] +
                               cf.m_leafc_storage_to_litter_patch[p] +
                               cf.m_frootc_storage_to_litter_patch[p] +
                               cf.m_leafc_xfer_to_litter_patch[p] +
                               cf.m_frootc_xfer_to_litter_patch[p] +
                               cf.m_livestemc_to_litter_patch[p] +
                               cf.m_livestemc_storage_to_litter_patch[p] +
                               cf.m_livestemc_xfer_to_litter_patch[p] +
                               cf.m_deadstemc_to_litter_patch[p] +
                               cf.m_deadstemc_storage_to_litter_patch[p] +
                               cf.m_deadstemc_xfer_to_litter_patch[p] +
                               cf.m_livecrootc_to_litter_patch[p] +
                               cf.m_livecrootc_storage_to_litter_patch[p] +
                               cf.m_livecrootc_xfer_to_litter_patch[p] +
                               cf.m_deadcrootc_to_litter_patch[p] +
                               cf.m_deadcrootc_storage_to_litter_patch[p] +
                               cf.m_deadcrootc_xfer_to_litter_patch[p] +
                               cf.m_gresp_storage_to_litter_patch[p] +
                               cf.m_gresp_xfer_to_litter_patch[p]

        # --- Fire carbon loss ---
        cf.fire_closs_patch[p] = cf.m_leafc_to_fire_patch[p] +
                                  cf.m_frootc_to_fire_patch[p] +
                                  cf.m_leafc_storage_to_fire_patch[p] +
                                  cf.m_frootc_storage_to_fire_patch[p] +
                                  cf.m_leafc_xfer_to_fire_patch[p] +
                                  cf.m_frootc_xfer_to_fire_patch[p] +
                                  cf.m_livestemc_to_fire_patch[p] +
                                  cf.m_livestemc_storage_to_fire_patch[p] +
                                  cf.m_livestemc_xfer_to_fire_patch[p] +
                                  cf.m_deadstemc_to_fire_patch[p] +
                                  cf.m_deadstemc_storage_to_fire_patch[p] +
                                  cf.m_deadstemc_xfer_to_fire_patch[p] +
                                  cf.m_livecrootc_to_fire_patch[p] +
                                  cf.m_livecrootc_storage_to_fire_patch[p] +
                                  cf.m_livecrootc_xfer_to_fire_patch[p] +
                                  cf.m_deadcrootc_to_fire_patch[p] +
                                  cf.m_deadcrootc_storage_to_fire_patch[p] +
                                  cf.m_deadcrootc_xfer_to_fire_patch[p] +
                                  cf.m_gresp_storage_to_fire_patch[p] +
                                  cf.m_gresp_xfer_to_fire_patch[p]

        # --- Allocation/loss summaries ---
        cf.frootc_alloc_patch[p] = cf.cpool_to_frootc_patch[p] +
                                    cf.cpool_to_frootc_storage_patch[p]

        cf.frootc_loss_patch[p] = cf.leafc_to_litter_patch[p] +
                                   cf.frootc_to_litter_patch[p]

        cf.leafc_alloc_patch[p] = cf.cpool_to_leafc_patch[p] +
                                   cf.cpool_to_leafc_storage_patch[p]

        cf.leafc_loss_patch[p] = cf.leafc_to_litter_patch[p]

        cf.woodc_alloc_patch[p] = cf.cpool_to_livestemc_patch[p] +
                                   cf.cpool_to_deadstemc_patch[p] +
                                   cf.cpool_to_livecrootc_patch[p] +
                                   cf.cpool_to_deadcrootc_patch[p] +
                                   cf.cpool_to_livestemc_storage_patch[p] +
                                   cf.cpool_to_deadstemc_storage_patch[p] +
                                   cf.cpool_to_livecrootc_storage_patch[p] +
                                   cf.cpool_to_deadcrootc_storage_patch[p]

        cf.woodc_loss_patch[p] = cf.m_livestemc_to_litter_patch[p] +
                                  cf.m_deadstemc_to_litter_patch[p] +
                                  cf.m_livecrootc_to_litter_patch[p] +
                                  cf.m_deadcrootc_to_litter_patch[p]
    end

    # Column-level p2c aggregation is stubbed pending subgridAveMod port
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
