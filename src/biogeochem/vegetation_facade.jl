# ==========================================================================
# Ported from: src/biogeochem/CNVegetationFacade.F90
# Vegetation facade: unified interface to the CN Vegetation subsystem.
#
# A "facade" provides a higher-level interface that makes the subsystem
# easier to use. This module coordinates calls to vegetation state,
# carbon/nitrogen state and flux, products, balance checking, fire,
# dynamic vegetation, and the CN driver.
#
# Public types:
#   CNVegetationConfig — control flags for the vegetation facade
#   CNVegetationData   — aggregated state for the CN vegetation subsystem
#
# Public functions:
#   cn_vegetation_init!                    — Initialize the vegetation facade
#   cn_vegetation_init2!                   — Post-subgrid-weight initialization
#   cn_vegetation_init_each_timestep!      — Zero fluxes at start of timestep
#   cn_vegetation_init_column_balance!     — Set starting point for col balance
#   cn_vegetation_init_gridcell_balance!   — Set starting point for gc balance
#   cn_vegetation_ecosystem_pre_drainage!  — Main CN science before drainage
#   cn_vegetation_ecosystem_post_drainage! — Main CN science after drainage
#   cn_vegetation_balance_check!           — C & N balance check
#   cn_vegetation_end_of_timestep!         — End-of-timestep veg dynamics
#   get_net_carbon_exchange_grc            — Gridcell net C exchange
#   get_leafn_patch                        — Patch leaf nitrogen
#   get_downreg_patch                      — Patch N downregulation
#   get_root_respiration_patch             — Patch root respiration
#   get_annsum_npp_patch                   — Patch annual sum NPP
#   get_agnpp_patch                        — Patch aboveground NPP
#   get_bgnpp_patch                        — Patch belowground NPP
#   get_froot_carbon_patch                 — Patch fine root carbon
#   get_croot_carbon_patch                 — Patch coarse root carbon
#   get_totvegc_col                        — Column total veg carbon
# ==========================================================================

# ---------------------------------------------------------------------------
# CNVegetationConfig — control flags for the vegetation facade
# ---------------------------------------------------------------------------

"""
    CNVegetationConfig

Configuration flags for the CN Vegetation facade.
Combines control variables from `clm_varctl`, namelist settings, and
module-level parameters from `CNVegetationFacade.F90`.
"""
Base.@kwdef mutable struct CNVegetationConfig
    use_cn::Bool = false
    use_c13::Bool = false
    use_c14::Bool = false
    use_cndv::Bool = false
    use_fates_bgc::Bool = false
    use_fun::Bool = false
    use_crop::Bool = false
    use_crop_agsys::Bool = false
    use_nitrif_denitrif::Bool = false
    use_matrixcn::Bool = false
    use_soil_matrixcn::Bool = false
    reseed_dead_plants::Bool = false
    dribble_crophrv_xsmrpool_2atm::Bool = false
    skip_steps::Int = 0
end

# ---------------------------------------------------------------------------
# CNVegetationData — aggregated state for the CN vegetation subsystem
# ---------------------------------------------------------------------------

"""
    CNVegetationData

Aggregated state for the CN Vegetation subsystem. Holds all vegetation
carbon/nitrogen state/flux types, products, balance checking, and
configuration.

Ported from `cn_vegetation_type` in `CNVegetationFacade.F90`.
"""
Base.@kwdef mutable struct CNVegetationData
    # Configuration
    config::CNVegetationConfig = CNVegetationConfig()

    # CN driver config (passed to cn_driver functions)
    driver_config::CNDriverConfig = CNDriverConfig()

    # Vegetation state
    cnveg_state_inst::CNVegStateData = CNVegStateData()

    # Carbon state (C12, C13, C14)
    cnveg_carbonstate_inst::CNVegCarbonStateData = CNVegCarbonStateData()
    c13_cnveg_carbonstate_inst::CNVegCarbonStateData = CNVegCarbonStateData()
    c14_cnveg_carbonstate_inst::CNVegCarbonStateData = CNVegCarbonStateData()

    # Carbon flux (C12, C13, C14)
    cnveg_carbonflux_inst::CNVegCarbonFluxData = CNVegCarbonFluxData()
    c13_cnveg_carbonflux_inst::CNVegCarbonFluxData = CNVegCarbonFluxData()
    c14_cnveg_carbonflux_inst::CNVegCarbonFluxData = CNVegCarbonFluxData()

    # Nitrogen state and flux
    cnveg_nitrogenstate_inst::CNVegNitrogenStateData = CNVegNitrogenStateData()
    cnveg_nitrogenflux_inst::CNVegNitrogenFluxData = CNVegNitrogenFluxData()
end

# ---------------------------------------------------------------------------
# cn_vegetation_init! — Initialize the vegetation facade
# Ported from cn_vegetation_type%Init in CNVegetationFacade.F90
# ---------------------------------------------------------------------------

"""
    cn_vegetation_init!(veg, np, nc, ng; kwargs...)

Initialize a `CNVegetationData` instance with given dimensions.

Allocates and initializes all sub-type instances (vegetation state,
carbon state/flux, nitrogen state/flux). Mirrors the Fortran `Init`
subroutine which calls Init on each sub-type conditionally based on
`use_cn`, `use_c13`, `use_c14`, etc.

Arguments:
- `veg`: CNVegetationData instance
- `np`: number of patches
- `nc`: number of columns
- `ng`: number of gridcells
- `nlevdecomp`: number of decomposition levels
- `ndecomp_pools`: number of decomposition pools
- `ndecomp_cascade_transitions`: number of cascade transitions
- `nrepr`: number of reproductive pools
"""
function cn_vegetation_init!(veg::CNVegetationData, np::Int, nc::Int, ng::Int;
                              nlevdecomp::Int = 1,
                              ndecomp_pools::Int = 7,
                              ndecomp_cascade_transitions::Int = 5,
                              nrepr::Int = 1,
                              i_litr_max::Int = 3)
    cfg = veg.config

    # Always initialize vegetation state
    cnveg_state_init!(veg.cnveg_state_inst, np, nc;
                       use_crop_agsys=cfg.use_crop_agsys)

    if cfg.use_cn || cfg.use_fates_bgc
        # Carbon state (C12)
        cnveg_carbon_state_init!(veg.cnveg_carbonstate_inst, np, nc, ng;
                                  nrepr=nrepr)

        # C13 carbon state
        if cfg.use_c13
            cnveg_carbon_state_init!(veg.c13_cnveg_carbonstate_inst, np, nc, ng;
                                      nrepr=nrepr)
        end

        # C14 carbon state
        if cfg.use_c14
            cnveg_carbon_state_init!(veg.c14_cnveg_carbonstate_inst, np, nc, ng;
                                      nrepr=nrepr)
        end

        # Carbon flux (C12)
        cnveg_carbon_flux_init!(veg.cnveg_carbonflux_inst, np, nc, ng;
                                 nrepr=nrepr,
                                 nlevdecomp_full=nlevdecomp,
                                 ndecomp_pools=ndecomp_pools)

        # C13 carbon flux
        if cfg.use_c13
            cnveg_carbon_flux_init!(veg.c13_cnveg_carbonflux_inst, np, nc, ng;
                                     nrepr=nrepr,
                                     nlevdecomp_full=nlevdecomp,
                                     ndecomp_pools=ndecomp_pools)
        end

        # C14 carbon flux
        if cfg.use_c14
            cnveg_carbon_flux_init!(veg.c14_cnveg_carbonflux_inst, np, nc, ng;
                                     nrepr=nrepr,
                                     nlevdecomp_full=nlevdecomp,
                                     ndecomp_pools=ndecomp_pools)
        end

        # Nitrogen state
        cnveg_nitrogen_state_init!(veg.cnveg_nitrogenstate_inst, np, nc, ng;
                                    nrepr=nrepr)

        # Nitrogen flux
        cnveg_nitrogen_flux_init!(veg.cnveg_nitrogenflux_inst, np, nc, ng;
                                   nrepr=nrepr,
                                   nlevdecomp_full=nlevdecomp,
                                   ndecomp_pools=ndecomp_pools,
                                   i_litr_max=i_litr_max)
    end

    # Synchronize driver config with facade config
    _sync_driver_config!(veg)

    return nothing
end

# ---------------------------------------------------------------------------
# cn_vegetation_init2! — Post-subgrid-weight initialization
# Ported from cn_vegetation_type%Init2 in CNVegetationFacade.F90
# ---------------------------------------------------------------------------

"""
    cn_vegetation_init2!(veg; bounds)

Do initialization needed after subgrid weights are determined.
Calls `cn_driver_init!` and (when use_cndv) dynamic vegetation init.

Should only be called if `use_cn` is true.

Ported from `Init2` in `CNVegetationFacade.F90`.
"""
function cn_vegetation_init2!(veg::CNVegetationData;
                                bounds::UnitRange{Int} = 1:0)
    cn_driver_init!(veg.driver_config; bounds=bounds)

    # dynCNDV_init — not yet ported
    # if veg.config.use_cndv
    #     dynCNDV_init(bounds, veg.dgvs_inst)
    # end

    return nothing
end

# ---------------------------------------------------------------------------
# cn_vegetation_init_each_timestep! — Zero fluxes at start of timestep
# Ported from cn_vegetation_type%InitEachTimeStep in CNVegetationFacade.F90
# ---------------------------------------------------------------------------

"""
    cn_vegetation_init_each_timestep!(veg; mask_soilc, mask_soilp,
                                      bounds_col, bounds_patch)

Zero DWT and GRU flux fields at the start of each time step.

Should only be called if `use_cn` is true.

Ported from `InitEachTimeStep` in `CNVegetationFacade.F90`.

Note: In the Fortran, this calls ZeroDWT and ZeroGRU on the carbon flux,
nitrogen flux, carbon state, and nitrogen state types. Those methods are
not yet ported; this function documents the call sequence as a placeholder.
"""
function cn_vegetation_init_each_timestep!(veg::CNVegetationData;
                                            mask_soilc::BitVector,
                                            mask_soilp::BitVector,
                                            bounds_col::UnitRange{Int},
                                            bounds_patch::UnitRange{Int})
    cfg = veg.config

    # ZeroDWT on carbon fluxes
    # cnveg_carbonflux_inst%ZeroDWT(bounds) — not yet ported
    # c13_cnveg_carbonflux_inst%ZeroDWT(bounds) — if use_c13
    # c14_cnveg_carbonflux_inst%ZeroDWT(bounds) — if use_c14
    # cnveg_nitrogenflux_inst%ZeroDWT(bounds) — not yet ported
    # cnveg_carbonstate_inst%ZeroDWT(bounds) — not yet ported
    # cnveg_nitrogenstate_inst%ZeroDWT(bounds) — not yet ported

    # ZeroGRU on carbon fluxes
    # cnveg_carbonflux_inst%ZeroGRU(bounds) — not yet ported
    # c13_cnveg_carbonflux_inst%ZeroGRU(bounds) — if use_c13
    # c14_cnveg_carbonflux_inst%ZeroGRU(bounds) — if use_c14
    # cnveg_nitrogenflux_inst%ZeroGRU(bounds) — not yet ported

    return nothing
end

# ---------------------------------------------------------------------------
# cn_vegetation_init_column_balance! — Set starting point for col balance
# Ported from cn_vegetation_type%InitColumnBalance in CNVegetationFacade.F90
# ---------------------------------------------------------------------------

"""
    cn_vegetation_init_column_balance!(veg; kwargs...)

Set the starting point for column-level balance checks.
Calls `cn_driver_summarize_states!` and then BeginCNColumnBalance.

Should be called after DynamicAreaConservation.

Ported from `InitColumnBalance` in `CNVegetationFacade.F90`.
"""
function cn_vegetation_init_column_balance!(veg::CNVegetationData;
        mask_bgc_soilc::BitVector,
        mask_bgc_vegp::BitVector,
        mask_allc::BitVector = mask_bgc_soilc,
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData)

    cn_driver_summarize_states!(veg.driver_config;
        mask_bgc_soilc=mask_bgc_soilc,
        mask_bgc_vegp=mask_bgc_vegp,
        mask_allc=mask_allc,
        bounds_col=bounds_col,
        bounds_patch=bounds_patch,
        cnveg_cs=veg.cnveg_carbonstate_inst,
        cnveg_ns=veg.cnveg_nitrogenstate_inst,
        soilbgc_cs=soilbgc_cs,
        soilbgc_ns=soilbgc_ns)

    # cn_balance_inst%BeginCNColumnBalance — not yet ported
    # Would set the starting column totals for balance checking

    return nothing
end

# ---------------------------------------------------------------------------
# cn_vegetation_init_gridcell_balance! — Set starting point for gc balance
# Ported from cn_vegetation_type%InitGridcellBalance in CNVegetationFacade.F90
# ---------------------------------------------------------------------------

"""
    cn_vegetation_init_gridcell_balance!(veg; kwargs...)

Set the starting point for gridcell-level balance checks.
Calls `cn_driver_summarize_states!`, column-to-gridcell aggregation,
and BeginCNGridcellBalance.

Should be called before DynamicAreaConservation.

Ported from `InitGridcellBalance` in `CNVegetationFacade.F90`.
"""
function cn_vegetation_init_gridcell_balance!(veg::CNVegetationData;
        mask_bgc_soilc::BitVector,
        mask_bgc_vegp::BitVector,
        mask_allc::BitVector = mask_bgc_soilc,
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData)

    cn_driver_summarize_states!(veg.driver_config;
        mask_bgc_soilc=mask_bgc_soilc,
        mask_bgc_vegp=mask_bgc_vegp,
        mask_allc=mask_allc,
        bounds_col=bounds_col,
        bounds_patch=bounds_patch,
        cnveg_cs=veg.cnveg_carbonstate_inst,
        cnveg_ns=veg.cnveg_nitrogenstate_inst,
        soilbgc_cs=soilbgc_cs,
        soilbgc_ns=soilbgc_ns)

    # c2g (column to gridcell) for total C and N — not yet ported
    # cn_balance_inst%BeginCNGridcellBalance — not yet ported

    return nothing
end

# ---------------------------------------------------------------------------
# cn_vegetation_ecosystem_pre_drainage! — Main CN science before drainage
# Ported from cn_vegetation_type%EcosystemDynamicsPreDrainage
# ---------------------------------------------------------------------------

"""
    cn_vegetation_ecosystem_pre_drainage!(veg; kwargs...)

Execute the main biogeochemistry science that needs to be done before
hydrology-drainage. Orchestrates calls to:
- Crop year increment (not yet ported)
- CNDriverNoLeaching (already ported)
- Fire carbon emissions update (not yet ported)
- Annual update (not yet ported)

Can be called for either `use_cn` or `use_fates_bgc`.

Ported from `EcosystemDynamicsPreDrainage` in `CNVegetationFacade.F90`.
"""
function cn_vegetation_ecosystem_pre_drainage!(veg::CNVegetationData;
        mask_bgc_soilc::BitVector,
        mask_bgc_vegp::BitVector,
        mask_pcropp::BitVector = falses(0),
        mask_soilnopcropp::BitVector = falses(0),
        mask_exposedvegp::BitVector = falses(0),
        mask_noexposedvegp::BitVector = falses(0),
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        nlevdecomp::Int,
        nlevdecomp_full::Int = nlevdecomp,
        ndecomp_pools::Int,
        ndecomp_cascade_transitions::Int,
        i_litr_min::Int,
        i_litr_max::Int,
        i_cwd::Int,
        npcropmin::Int = 17,
        nrepr::Int = 1,
        patch_column::Vector{Int},
        ivt::Vector{Int},
        woody::Vector{<:Real},
        harvdate::Vector{Int},
        col_is_fates::Vector{Bool},
        cascade_donor_pool::Vector{Int},
        cascade_receiver_pool::Vector{Int},
        dt::Real,
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_cf::SoilBiogeochemCarbonFluxData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData,
        soilbgc_nf::SoilBiogeochemNitrogenFluxData,
        soilbgc_state::SoilBiogeochemStateData,
        # Decomposition infrastructure (optional, for full BGC)
        cascade_con::Union{DecompCascadeConData, Nothing} = nothing,
        decomp_bgc_state::Union{DecompBGCState, Nothing} = nothing,
        decomp_bgc_params::Union{DecompBGCParams, Nothing} = nothing,
        cn_shared_params::Union{CNSharedParamsData, Nothing} = nothing,
        decomp_params::Union{DecompParams, Nothing} = nothing,
        competition_state::Union{SoilBGCCompetitionState, Nothing} = nothing,
        competition_params::Union{SoilBGCCompetitionParams, Nothing} = nothing,
        litter_params::Union{LitterVertTranspParams, Nothing} = nothing,
        t_soisno::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        soilpsi::Union{Matrix{<:Real}, Nothing} = nothing,
        col::Union{ColumnData, Nothing} = nothing,
        grc::Union{GridcellData, Nothing} = nothing,
        active_layer::Union{ActiveLayerData, Nothing} = nothing,
        dzsoi_decomp::Union{Vector{<:Real}, Nothing} = nothing,
        zsoi_vals::Union{Vector{<:Real}, Nothing} = nothing,
        zisoi_vals::Union{Vector{<:Real}, Nothing} = nothing,
        mask_actfirec::BitVector = falses(length(bounds_col)),
        mask_actfirep::BitVector = falses(length(bounds_patch)))

    # crop_inst%CropIncrementYear — not yet ported

    # CNDriverNoLeaching — already ported
    cn_driver_no_leaching!(veg.driver_config;
        mask_bgc_soilc=mask_bgc_soilc,
        mask_bgc_vegp=mask_bgc_vegp,
        mask_pcropp=mask_pcropp,
        mask_soilnopcropp=mask_soilnopcropp,
        mask_exposedvegp=mask_exposedvegp,
        mask_noexposedvegp=mask_noexposedvegp,
        bounds_col=bounds_col,
        bounds_patch=bounds_patch,
        nlevdecomp=nlevdecomp,
        nlevdecomp_full=nlevdecomp_full,
        ndecomp_pools=ndecomp_pools,
        ndecomp_cascade_transitions=ndecomp_cascade_transitions,
        i_litr_min=i_litr_min,
        i_litr_max=i_litr_max,
        i_cwd=i_cwd,
        npcropmin=npcropmin,
        nrepr=nrepr,
        patch_column=patch_column,
        ivt=ivt,
        woody=woody,
        harvdate=harvdate,
        col_is_fates=col_is_fates,
        cascade_donor_pool=cascade_donor_pool,
        cascade_receiver_pool=cascade_receiver_pool,
        dt=dt,
        cnveg_cs=veg.cnveg_carbonstate_inst,
        cnveg_cf=veg.cnveg_carbonflux_inst,
        cnveg_ns=veg.cnveg_nitrogenstate_inst,
        cnveg_nf=veg.cnveg_nitrogenflux_inst,
        soilbgc_cs=soilbgc_cs,
        soilbgc_cf=soilbgc_cf,
        soilbgc_ns=soilbgc_ns,
        soilbgc_nf=soilbgc_nf,
        soilbgc_state=soilbgc_state,
        cascade_con=cascade_con,
        decomp_bgc_state=decomp_bgc_state,
        decomp_bgc_params=decomp_bgc_params,
        cn_shared_params=cn_shared_params,
        decomp_params=decomp_params,
        competition_state=competition_state,
        competition_params=competition_params,
        litter_params=litter_params,
        t_soisno=t_soisno,
        soilpsi=soilpsi,
        col=col,
        grc=grc,
        active_layer=active_layer,
        dzsoi_decomp=dzsoi_decomp,
        zsoi_vals=zsoi_vals,
        zisoi_vals=zisoi_vals,
        mask_actfirec=mask_actfirec,
        mask_actfirep=mask_actfirep)

    # CNFireEmisUpdate — not yet ported
    # CNAnnualUpdate — not yet ported

    return nothing
end

# ---------------------------------------------------------------------------
# cn_vegetation_ecosystem_post_drainage! — Main CN science after drainage
# Ported from cn_vegetation_type%EcosystemDynamicsPostDrainage
# ---------------------------------------------------------------------------

"""
    cn_vegetation_ecosystem_post_drainage!(veg; kwargs...)

Execute the main biogeochemistry science after hydrology-drainage.
Orchestrates calls to:
- CNDriverLeaching (already ported)
- CNPrecisionControl (not yet ported)
- SoilBiogeochemPrecisionControl (not yet ported)
- CNDriverSummarizeStates (already ported)
- CNDriverSummarizeFluxes (already ported)
- CNVegStructUpdate (not yet ported)

Should only be called if `use_cn` or `use_fates_bgc` is true.

Ported from `EcosystemDynamicsPostDrainage` in `CNVegetationFacade.F90`.
"""
function cn_vegetation_ecosystem_post_drainage!(veg::CNVegetationData;
        mask_bgc_soilc::BitVector,
        mask_bgc_vegp::BitVector,
        mask_allc::BitVector = mask_bgc_soilc,
        mask_actfirec::BitVector = falses(0),
        mask_actfirep::BitVector = falses(0),
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        nlevdecomp::Int,
        ndecomp_pools::Int = 7,
        ndecomp_cascade_transitions::Int = 7,
        i_litr_min::Int = 1,
        i_litr_max::Int = 3,
        i_cwd::Int = 4,
        dt::Real,
        doalb::Bool = false,
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_cf::SoilBiogeochemCarbonFluxData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData,
        soilbgc_nf::SoilBiogeochemNitrogenFluxData)

    # CNDriverLeaching — already ported
    cn_driver_leaching!(veg.driver_config;
        mask_bgc_soilc=mask_bgc_soilc,
        mask_bgc_vegp=mask_bgc_vegp,
        mask_actfirec=mask_actfirec,
        mask_actfirep=mask_actfirep,
        bounds_col=bounds_col,
        bounds_patch=bounds_patch,
        nlevdecomp=nlevdecomp,
        ndecomp_pools=ndecomp_pools,
        ndecomp_cascade_transitions=ndecomp_cascade_transitions,
        i_litr_min=i_litr_min,
        i_litr_max=i_litr_max,
        i_cwd=i_cwd,
        dt=dt,
        soilbgc_ns=soilbgc_ns,
        soilbgc_nf=soilbgc_nf,
        cnveg_ns=veg.cnveg_nitrogenstate_inst,
        cnveg_nf=veg.cnveg_nitrogenflux_inst)

    # CNPrecisionControl — not yet ported
    # SoilBiogeochemPrecisionControl — not yet ported

    # CNDriverSummarizeStates
    cn_driver_summarize_states!(veg.driver_config;
        mask_bgc_soilc=mask_bgc_soilc,
        mask_bgc_vegp=mask_bgc_vegp,
        mask_allc=mask_allc,
        bounds_col=bounds_col,
        bounds_patch=bounds_patch,
        cnveg_cs=veg.cnveg_carbonstate_inst,
        cnveg_ns=veg.cnveg_nitrogenstate_inst,
        soilbgc_cs=soilbgc_cs,
        soilbgc_ns=soilbgc_ns)

    # CNDriverSummarizeFluxes
    cn_driver_summarize_fluxes!(veg.driver_config;
        mask_bgc_soilc=mask_bgc_soilc,
        mask_bgc_vegp=mask_bgc_vegp,
        bounds_col=bounds_col,
        bounds_patch=bounds_patch,
        cnveg_cf=veg.cnveg_carbonflux_inst,
        cnveg_nf=veg.cnveg_nitrogenflux_inst,
        soilbgc_cf=soilbgc_cf,
        soilbgc_nf=soilbgc_nf)

    # CNVegStructUpdate — not yet ported
    # if num_bgc_vegp > 0 && doalb
    #     CNVegStructUpdate(...)
    # end

    return nothing
end

# ---------------------------------------------------------------------------
# cn_vegetation_balance_check! — C & N balance check
# Ported from cn_vegetation_type%BalanceCheck in CNVegetationFacade.F90
# ---------------------------------------------------------------------------

"""
    cn_vegetation_balance_check!(veg; kwargs...)

Check the carbon and nitrogen balance. Skips check during the first
`skip_steps` timesteps after startup.

Should only be called if `use_cn` or `use_fates_bgc` is true.

Ported from `BalanceCheck` in `CNVegetationFacade.F90`.
"""
function cn_vegetation_balance_check!(veg::CNVegetationData;
        mask_bgc_soilc::BitVector,
        bounds_col::UnitRange{Int},
        nstep_since_startup::Int = 1,
        soilbgc_cf::SoilBiogeochemCarbonFluxData,
        soilbgc_nf::SoilBiogeochemNitrogenFluxData,
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData)

    if nstep_since_startup <= veg.config.skip_steps
        # Skip balance check for first timesteps after startup
        return nothing
    end

    # cn_balance_inst%CBalanceCheck — not yet ported
    # cn_balance_inst%NBalanceCheck — not yet ported

    return nothing
end

# ---------------------------------------------------------------------------
# cn_vegetation_end_of_timestep! — End-of-timestep veg dynamics
# Ported from cn_vegetation_type%EndOfTimeStepVegDynamics
# ---------------------------------------------------------------------------

"""
    cn_vegetation_end_of_timestep!(veg; kwargs...)

Do vegetation dynamics that should be done at the end of each time step.
When `use_cndv` is true and it is the end of the current year (and not
the first step), calls the dynamic vegetation driver.

Should only be called if `use_cn` is true.

Ported from `EndOfTimeStepVegDynamics` in `CNVegetationFacade.F90`.
"""
function cn_vegetation_end_of_timestep!(veg::CNVegetationData;
        bounds_patch::UnitRange{Int} = 1:0,
        is_end_curr_year::Bool = false,
        is_first_step::Bool = false)

    if veg.config.use_cndv
        if is_end_curr_year && !is_first_step
            # CNDVDriver — not yet ported
            # dynCNDV driver would be called here
        end
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Getter functions — provide access to internal state arrays
# Ported from cn_vegetation_type getters in CNVegetationFacade.F90
# ---------------------------------------------------------------------------

"""
    get_net_carbon_exchange_grc(veg, bounds_grc)

Get gridcell-level net carbon exchange (positive = source).
Returns -nbp_grc when use_cn, otherwise zeros.

Ported from `get_net_carbon_exchange_grc` in `CNVegetationFacade.F90`.
"""
function get_net_carbon_exchange_grc(veg::CNVegetationData,
                                      bounds_grc::UnitRange{Int})
    if veg.config.use_cn
        return [-veg.cnveg_carbonflux_inst.nbp_grc[g] for g in bounds_grc]
    else
        FT = eltype(veg.cnveg_carbonflux_inst.nbp_grc)
        return zeros(FT, length(bounds_grc))
    end
end

"""
    get_leafn_patch(veg, bounds_patch)

Get patch-level leaf nitrogen array (gN/m2).
Returns leafn_patch when use_cn, otherwise NaN.

Ported from `get_leafn_patch` in `CNVegetationFacade.F90`.
"""
function get_leafn_patch(veg::CNVegetationData,
                           bounds_patch::UnitRange{Int})
    if veg.config.use_cn
        return veg.cnveg_nitrogenstate_inst.leafn_patch[bounds_patch]
    else
        return fill(NaN, length(bounds_patch))
    end
end

"""
    get_downreg_patch(veg, bounds_patch)

Get patch-level fractional reduction in GPP due to N limitation.
Returns downreg_patch when use_cn, otherwise NaN.

Ported from `get_downreg_patch` in `CNVegetationFacade.F90`.
"""
function get_downreg_patch(veg::CNVegetationData,
                             bounds_patch::UnitRange{Int})
    if veg.config.use_cn
        return veg.cnveg_state_inst.downreg_patch[bounds_patch]
    else
        return fill(NaN, length(bounds_patch))
    end
end

"""
    get_root_respiration_patch(veg, bounds_patch)

Get patch-level root respiration (gC/m2/s).
Returns rr_patch when use_cn, otherwise NaN.

Ported from `get_root_respiration_patch` in `CNVegetationFacade.F90`.
"""
function get_root_respiration_patch(veg::CNVegetationData,
                                      bounds_patch::UnitRange{Int})
    if veg.config.use_cn
        return veg.cnveg_carbonflux_inst.rr_patch[bounds_patch]
    else
        return fill(NaN, length(bounds_patch))
    end
end

"""
    get_annsum_npp_patch(veg, bounds_patch)

Get patch-level annual sum NPP (gC/m2/yr).
Returns annsum_npp_patch when use_cn, otherwise NaN.

Ported from `get_annsum_npp_patch` in `CNVegetationFacade.F90`.
"""
function get_annsum_npp_patch(veg::CNVegetationData,
                                bounds_patch::UnitRange{Int})
    if veg.config.use_cn
        return veg.cnveg_carbonflux_inst.annsum_npp_patch[bounds_patch]
    else
        return fill(NaN, length(bounds_patch))
    end
end

"""
    get_agnpp_patch(veg, bounds_patch)

Get patch-level aboveground NPP (gC/m2/s).
Returns agnpp_patch when use_cn, otherwise NaN.

Ported from `get_agnpp_patch` in `CNVegetationFacade.F90`.
"""
function get_agnpp_patch(veg::CNVegetationData,
                           bounds_patch::UnitRange{Int})
    if veg.config.use_cn
        return veg.cnveg_carbonflux_inst.agnpp_patch[bounds_patch]
    else
        return fill(NaN, length(bounds_patch))
    end
end

"""
    get_bgnpp_patch(veg, bounds_patch)

Get patch-level belowground NPP (gC/m2/s).
Returns bgnpp_patch when use_cn, otherwise NaN.

Ported from `get_bgnpp_patch` in `CNVegetationFacade.F90`.
"""
function get_bgnpp_patch(veg::CNVegetationData,
                           bounds_patch::UnitRange{Int})
    if veg.config.use_cn
        return veg.cnveg_carbonflux_inst.bgnpp_patch[bounds_patch]
    else
        return fill(NaN, length(bounds_patch))
    end
end

"""
    get_froot_carbon_patch(veg, bounds_patch; tlai, slatop, froot_leaf, ivt)

Get patch-level fine root carbon (gC/m2).
When use_cn, returns frootc_patch. Otherwise estimates from LAI and PFT
traits: froot_c = (LAI / slatop) * froot_leaf.

Ported from `get_froot_carbon_patch` in `CNVegetationFacade.F90`.
"""
function get_froot_carbon_patch(veg::CNVegetationData,
                                  bounds_patch::UnitRange{Int};
                                  tlai::Vector{<:Real} = Float64[],
                                  slatop::Vector{<:Real} = Float64[],
                                  froot_leaf::Vector{<:Real} = Float64[],
                                  ivt::Vector{Int} = Int[])
    if veg.config.use_cn
        return veg.cnveg_carbonstate_inst.frootc_patch[bounds_patch]
    else
        FT = eltype(tlai)
        result = zeros(FT, length(bounds_patch))
        for (i, p) in enumerate(bounds_patch)
            pft = ivt[p]
            if slatop[pft] > 0.0
                result[i] = tlai[p] / slatop[pft] * froot_leaf[pft]
            end
        end
        return result
    end
end

"""
    get_croot_carbon_patch(veg, bounds_patch; tlai, slatop, stem_leaf, croot_stem, ivt)

Get patch-level live coarse root carbon (gC/m2).
When use_cn, returns livecrootc_patch. Otherwise estimates from LAI and
PFT traits: croot_c = (LAI / slatop) * stem_leaf * croot_stem.

Ported from `get_croot_carbon_patch` in `CNVegetationFacade.F90`.
"""
function get_croot_carbon_patch(veg::CNVegetationData,
                                  bounds_patch::UnitRange{Int};
                                  tlai::Vector{<:Real} = Float64[],
                                  slatop::Vector{<:Real} = Float64[],
                                  stem_leaf::Vector{<:Real} = Float64[],
                                  croot_stem::Vector{<:Real} = Float64[],
                                  ivt::Vector{Int} = Int[])
    if veg.config.use_cn
        return veg.cnveg_carbonstate_inst.livecrootc_patch[bounds_patch]
    else
        FT = eltype(tlai)
        result = zeros(FT, length(bounds_patch))
        for (i, p) in enumerate(bounds_patch)
            pft = ivt[p]
            if slatop[pft] > 0.0
                result[i] = tlai[p] / slatop[pft] * stem_leaf[pft] * croot_stem[pft]
            end
        end
        return result
    end
end

"""
    get_totvegc_col(veg, bounds_col)

Get column-level total vegetation carbon (gC/m2).
Returns totvegc_col when use_cn, otherwise NaN.

Ported from `get_totvegc_col` in `CNVegetationFacade.F90`.
"""
function get_totvegc_col(veg::CNVegetationData,
                           bounds_col::UnitRange{Int})
    if veg.config.use_cn
        return veg.cnveg_carbonstate_inst.totvegc_col[bounds_col]
    else
        return fill(NaN, length(bounds_col))
    end
end

# ---------------------------------------------------------------------------
# Internal helper: synchronize driver config from facade config
# ---------------------------------------------------------------------------

"""
    _sync_driver_config!(veg)

Copy configuration flags from the facade config to the CN driver config.
"""
function _sync_driver_config!(veg::CNVegetationData)
    fc = veg.config
    dc = veg.driver_config

    dc.use_cn = fc.use_cn
    dc.use_c13 = fc.use_c13
    dc.use_c14 = fc.use_c14
    dc.use_fun = fc.use_fun
    dc.use_crop = fc.use_crop
    dc.use_crop_agsys = fc.use_crop_agsys
    dc.use_nitrif_denitrif = fc.use_nitrif_denitrif
    dc.use_fates_bgc = fc.use_fates_bgc
    dc.use_matrixcn = fc.use_matrixcn
    dc.use_soil_matrixcn = fc.use_soil_matrixcn
    dc.dribble_crophrv_xsmrpool_2atm = fc.dribble_crophrv_xsmrpool_2atm

    return nothing
end
