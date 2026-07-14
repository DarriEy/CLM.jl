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
#
# STATUS NOTE (read before believing any "not yet ported" comment below — the
# ones that were stale have been corrected in place): this facade is LIVE — the
# driver calls cn_vegetation_init_each_timestep!/_init_*_balance!/
# _ecosystem_pre_drainage!/_post_drainage!/_balance_check!/_end_of_timestep!
# from clm_drv!. But several Fortran calls it makes are NOT re-issued from here
# even though the Julia routine EXISTS. Where a routine is ported, the comment
# below now names the file — and, when the live driver invokes it on another
# path (e.g. CNVegStructUpdate, CNPrecisionControl, CNDVDriver), says so. Read
# "ported, not wired here" as exactly that: the code exists, this facade does
# not call it.
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
    use_flexiblecn::Bool = false
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
# Parametric on working precision FT. FT is a *construction-precision tag*: the keyword
# constructor builds the state children at FT, but the fields stay LOOSE (un-pinned
# UnionAll) types so the AD path can swap in `Dual`-typed sub-instances on a Float64
# container. NOT @kwdef: a single type parameter can't both build children at FT via
# `CNVegetationData{FT}()` and keep @kwdef's synthesised constructors (same signature),
# so we use explicit constructors instead.
mutable struct CNVegetationData{FT<:Real}
    # Configuration
    config::CNVegetationConfig

    # CN driver config (passed to cn_driver functions)
    driver_config::CNDriverConfig

    # Vegetation state (loose types: see note above)
    cnveg_state_inst::CNVegStateData

    # Carbon state (C12, C13, C14)
    cnveg_carbonstate_inst::CNVegCarbonStateData
    c13_cnveg_carbonstate_inst::CNVegCarbonStateData
    c14_cnveg_carbonstate_inst::CNVegCarbonStateData

    # Carbon flux (C12, C13, C14)
    cnveg_carbonflux_inst::CNVegCarbonFluxData
    c13_cnveg_carbonflux_inst::CNVegCarbonFluxData
    c14_cnveg_carbonflux_inst::CNVegCarbonFluxData

    # Nitrogen state and flux
    cnveg_nitrogenstate_inst::CNVegNitrogenStateData
    cnveg_nitrogenflux_inst::CNVegNitrogenFluxData

    # C/N mass-conservation check (Fortran cn_balance_type). Holds the
    # beginning/end-of-step column + gridcell C and N masses.
    cn_balance_inst::CNBalanceData

    # Wood/crop product pools (Fortran c_products_inst / n_products_inst). The
    # balance check needs them as a C/N sink; they are also the destination of
    # harvest/gross-unrepresented-disturbance fluxes.
    c_products_inst::CNProductsData
    n_products_inst::CNProductsData
end

# Keyword constructor (field-name keyword API replacing @kwdef's). State children
# default to working precision FT; config sub-structs are overridable as before.
function CNVegetationData{FT}(;
        config::CNVegetationConfig = CNVegetationConfig(),
        driver_config::CNDriverConfig = CNDriverConfig(),
        cnveg_state_inst::CNVegStateData            = CNVegStateData{FT}(),
        cnveg_carbonstate_inst::CNVegCarbonStateData     = CNVegCarbonStateData{FT}(),
        c13_cnveg_carbonstate_inst::CNVegCarbonStateData = CNVegCarbonStateData{FT}(),
        c14_cnveg_carbonstate_inst::CNVegCarbonStateData = CNVegCarbonStateData{FT}(),
        cnveg_carbonflux_inst::CNVegCarbonFluxData       = CNVegCarbonFluxData{FT}(),
        c13_cnveg_carbonflux_inst::CNVegCarbonFluxData   = CNVegCarbonFluxData{FT}(),
        c14_cnveg_carbonflux_inst::CNVegCarbonFluxData   = CNVegCarbonFluxData{FT}(),
        cnveg_nitrogenstate_inst::CNVegNitrogenStateData = CNVegNitrogenStateData{FT}(),
        cnveg_nitrogenflux_inst::CNVegNitrogenFluxData   = CNVegNitrogenFluxData{FT}(),
        cn_balance_inst::CNBalanceData                   = CNBalanceData{FT, Vector{FT}}(),
        c_products_inst::CNProductsData                  = CNProductsData{FT}(),
        n_products_inst::CNProductsData                  = CNProductsData{FT}()) where {FT<:Real}
    CNVegetationData{FT}(config, driver_config, cnveg_state_inst,
        cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst,
        cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,
        cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,
        cn_balance_inst, c_products_inst, n_products_inst)
end

# Default precision is Float64; legacy callers that don't specify FT get a Float64 tree.
CNVegetationData(; kwargs...) = CNVegetationData{Float64}(; kwargs...)

# Device-movable: adapt each child (all CNVeg* state/flux sub-structs already carry
# Adapt.@adapt_structure; the config/driver_config structs hold no arrays and pass
# through as identity). Without this the CN-veg state tree stays on the host when the
# enclosing CLMInstances is moved to a GPU, and the CN kernels receive host arrays.
# Written by hand rather than via @adapt_structure because the facade's FT type
# parameter isn't the type of any field, so the macro's UnionAll positional
# reconstruction can't infer it; we recover FT from an adapted child instead.
function Adapt.adapt_structure(to, x::CNVegetationData)
    cs    = Adapt.adapt(to, x.cnveg_carbonstate_inst)
    c13cs = Adapt.adapt(to, x.c13_cnveg_carbonstate_inst)
    c14cs = Adapt.adapt(to, x.c14_cnveg_carbonstate_inst)
    cf    = Adapt.adapt(to, x.cnveg_carbonflux_inst)
    c13cf = Adapt.adapt(to, x.c13_cnveg_carbonflux_inst)
    c14cf = Adapt.adapt(to, x.c14_cnveg_carbonflux_inst)
    st    = Adapt.adapt(to, x.cnveg_state_inst)
    ns    = Adapt.adapt(to, x.cnveg_nitrogenstate_inst)
    nf    = Adapt.adapt(to, x.cnveg_nitrogenflux_inst)
    bal   = Adapt.adapt(to, x.cn_balance_inst)
    cprod = Adapt.adapt(to, x.c_products_inst)
    nprod = Adapt.adapt(to, x.n_products_inst)
    FT    = typeof(cs).parameters[1]
    return CNVegetationData{FT}(x.config, x.driver_config, st,
        cs, c13cs, c14cs, cf, c13cf, c14cf, ns, nf,
        bal, cprod, nprod)
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

        # C/N mass-balance check state + the wood/crop product pools it needs as a
        # C/N sink. Fortran allocates all three in CNVegetationFacade::Init.
        cn_balance_init!(veg.cn_balance_inst, nc, ng)
        cn_products_init!(veg.c_products_inst, ng)
        cn_products_init!(veg.n_products_inst, ng)
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

    # dynCNDV_init — PORTED as dyn_cndv_init! (biogeochem/cndv.jl) and called from
    # the live init path (driver/clm_initialize.jl:197), not from here.
    # if veg.config.use_cndv
    #     dyn_cndv_init!(veg.dgvs_inst, patch, bounds)
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
nitrogen flux, carbon state, and nitrogen state types. Those methods ARE ported
(`cnveg_carbon_flux_zero_dwt!` / `_zero_gru!`, `cnveg_nitrogen_flux_zero_dwt!` /
`_zero_gru!`, `cnveg_carbon_state_zero_dwt!`, `cnveg_nitrogen_state_zero_dwt!`,
in `src/types/cn_veg_*.jl`) — they are simply not invoked from this facade.
Ported, not wired here.
"""
function cn_vegetation_init_each_timestep!(veg::CNVegetationData;
                                            mask_soilc::AbstractVector{Bool},
                                            mask_soilp::AbstractVector{Bool},
                                            bounds_col::UnitRange{Int},
                                            bounds_patch::UnitRange{Int})
    cfg = veg.config

    # ZeroDWT / ZeroGRU — all four ported (see the docstring), none called here.
    # cnveg_carbonflux_inst%ZeroDWT   → cnveg_carbon_flux_zero_dwt!   (ported, not wired here)
    # c13_/c14_ carbonflux ZeroDWT    → same routine on the isotope insts (if use_c13/use_c14)
    # cnveg_nitrogenflux_inst%ZeroDWT → cnveg_nitrogen_flux_zero_dwt! (ported, not wired here)
    # cnveg_carbonstate_inst%ZeroDWT  → cnveg_carbon_state_zero_dwt!  (ported, not wired here)
    # cnveg_nitrogenstate_inst%ZeroDWT→ cnveg_nitrogen_state_zero_dwt!(ported, not wired here)
    # cnveg_carbonflux_inst%ZeroGRU   → cnveg_carbon_flux_zero_gru!   (ported, not wired here)
    # cnveg_nitrogenflux_inst%ZeroGRU → cnveg_nitrogen_flux_zero_gru! (ported, not wired here)

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
        mask_bgc_soilc::AbstractVector{Bool},
        mask_bgc_vegp::AbstractVector{Bool},
        mask_allc::AbstractVector{Bool} = mask_bgc_soilc,
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData,
        # Soil-BGC state-summary inputs; needed so totc_col/totn_col (the balance
        # check's BEGIN mass) are recomputed here, as Fortran does.
        cascade_con::Union{DecompCascadeConData, Nothing} = nothing,
        col::Union{ColumnData, Nothing} = nothing,
        patch::Union{PatchData, Nothing} = nothing,
        nlevdecomp::Int = 0,
        ndecomp_pools::Int = 0,
        dzsoi_decomp::Union{AbstractVector{<:Real}, Nothing} = nothing,
        zisoi_vals::Union{AbstractVector{<:Real}, Nothing} = nothing)

    cn_driver_summarize_states!(veg.driver_config;
        mask_bgc_soilc=mask_bgc_soilc,
        mask_bgc_vegp=mask_bgc_vegp,
        mask_allc=mask_allc,
        bounds_col=bounds_col,
        bounds_patch=bounds_patch,
        cnveg_cs=veg.cnveg_carbonstate_inst,
        cnveg_ns=veg.cnveg_nitrogenstate_inst,
        soilbgc_cs=soilbgc_cs,
        soilbgc_ns=soilbgc_ns,
        cascade_con=cascade_con, col=col, patch=patch,
        nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
        dzsoi_decomp=dzsoi_decomp, zisoi_vals=zisoi_vals)

    # cn_balance_inst%BeginCNColumnBalance — WIRED. Seeds begcb_col/begnb_col from
    # the totc_col/totn_col just recomputed by cn_driver_summarize_states! above.
    # Without this the balance check has no "before" mass and cannot run at all —
    # which is why the check that would have caught the dead N cycle was itself off.
    begin_cn_column_balance!(veg.cn_balance_inst, soilbgc_cs, soilbgc_ns,
                             mask_bgc_soilc, bounds_col)

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
        mask_bgc_soilc::AbstractVector{Bool},
        mask_bgc_vegp::AbstractVector{Bool},
        mask_allc::AbstractVector{Bool} = mask_bgc_soilc,
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData,
        # Needed for the c2g aggregation below; omitted => gridcell balance not seeded.
        col::Union{ColumnData, Nothing} = nothing,
        bounds_grc::UnitRange{Int} = 1:0,
        cascade_con::Union{DecompCascadeConData, Nothing} = nothing,
        patch::Union{PatchData, Nothing} = nothing,
        nlevdecomp::Int = 0,
        ndecomp_pools::Int = 0,
        dzsoi_decomp::Union{AbstractVector{<:Real}, Nothing} = nothing,
        zisoi_vals::Union{AbstractVector{<:Real}, Nothing} = nothing)

    cn_driver_summarize_states!(veg.driver_config;
        mask_bgc_soilc=mask_bgc_soilc,
        mask_bgc_vegp=mask_bgc_vegp,
        mask_allc=mask_allc,
        bounds_col=bounds_col,
        bounds_patch=bounds_patch,
        cnveg_cs=veg.cnveg_carbonstate_inst,
        cnveg_ns=veg.cnveg_nitrogenstate_inst,
        soilbgc_cs=soilbgc_cs,
        soilbgc_ns=soilbgc_ns,
        cascade_con=cascade_con, col=col, patch=patch,
        nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
        dzsoi_decomp=dzsoi_decomp, zisoi_vals=zisoi_vals)

    # c2g (column -> gridcell) for total C and N, then BeginCNGridcellBalance.
    # WIRED. Fortran InitGridcellBalance does exactly this pair; n_balance_check!
    # recomputes the c2g for the END-of-step mass, so the BEGIN side must be
    # seeded the same way or begnb_grc stays NaN and the gridcell check is dead.
    if col !== nothing && !isempty(bounds_grc)
        c2g_unity!(soilbgc_cs.totc_grc, soilbgc_cs.totc_col,
                   col.gridcell, col.wtgcell, bounds_col, bounds_grc)
        c2g_unity!(soilbgc_ns.totn_grc, soilbgc_ns.totn_col,
                   col.gridcell, col.wtgcell, bounds_col, bounds_grc)
        begin_cn_gridcell_balance!(veg.cn_balance_inst, soilbgc_cs, soilbgc_ns,
                                   veg.c_products_inst, veg.n_products_inst,
                                   bounds_grc)
    end

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
- Crop year increment (`crop_increment_year!`, types/crop.jl — ported, not wired here)
- CNDriverNoLeaching (already ported)
- Fire carbon/trace-gas EMISSIONS update (CNFireEmissionsMod — genuinely not
  ported; note the fire model itself IS ported: biogeochem/fire_li2014.jl …
  fire_li2024.jl, dispatched via fire_factory.jl)
- Annual update (`cn_annual_update!`, biogeochem/cn_annual_update.jl — ported,
  not wired here)

Can be called for either `use_cn` or `use_fates_bgc`.

Ported from `EcosystemDynamicsPreDrainage` in `CNVegetationFacade.F90`.
"""
function cn_vegetation_ecosystem_pre_drainage!(veg::CNVegetationData;
        mask_bgc_soilc::AbstractVector{Bool},
        mask_bgc_vegp::AbstractVector{Bool},
        mask_pcropp::AbstractVector{Bool} = falses(0),
        mask_soilnopcropp::AbstractVector{Bool} = falses(0),
        mask_exposedvegp::AbstractVector{Bool} = falses(0),
        mask_noexposedvegp::AbstractVector{Bool} = falses(0),
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
        patch_column::AbstractVector{<:Integer},
        ivt::AbstractVector{<:Integer},
        woody::AbstractVector{<:Real},
        harvdate::AbstractVector{<:Integer},
        col_is_fates::AbstractVector{Bool},
        cascade_donor_pool::AbstractVector{<:Integer},
        cascade_receiver_pool::AbstractVector{<:Integer},
        dt::Real,
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_cf::SoilBiogeochemCarbonFluxData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData,
        soilbgc_nf::SoilBiogeochemNitrogenFluxData,
        soilbgc_state::SoilBiogeochemStateData,
        # C13/C14 soil-BGC state/flux (the CIsoFlux* cascade sinks; supplied by the
        # driver only when use_c13/use_c14 + the state is allocated — else nothing).
        c13_soilbgc_cs::Union{SoilBiogeochemCarbonStateData, Nothing} = nothing,
        c13_soilbgc_cf::Union{SoilBiogeochemCarbonFluxData, Nothing} = nothing,
        c14_soilbgc_cs::Union{SoilBiogeochemCarbonStateData, Nothing} = nothing,
        c14_soilbgc_cf::Union{SoilBiogeochemCarbonFluxData, Nothing} = nothing,
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
        # Vegetation-flux inputs (threaded to cn_driver_no_leaching! for the veg-CN
        # flux chain: maintenance respiration, gpp/allocation, …).
        patch::Union{PatchData, Nothing} = nothing,
        pftcon_main::Union{Any, Nothing} = nothing,
        crop::Union{CropData, Nothing} = nothing,
        photosyns::Union{PhotosynthesisData, Nothing} = nothing,
        canopystate::Union{CanopyStateData, Nothing} = nothing,
        soilstate::Union{SoilStateData, Nothing} = nothing,
        temperature::Union{TemperatureData, Nothing} = nothing,
        water_diag::Union{WaterDiagnosticBulkData, Nothing} = nothing,
        waterstate::Union{WaterStateData, Nothing} = nothing,
        gridcell::Union{GridcellData, Nothing} = nothing,
        is_first_step::Bool = false,
        # Accumulated forcing inputs (Fortran wateratm2lndbulk_type). These feed the
        # Li fire schemes (fuel moisture / ignition) and CNPhenology's rain-triggered
        # stress-deciduous onset. They used to be absent here entirely, so cn_driver's
        # kwarg defaults (empty Float64[]) reached the fire modules — which then read
        # off the end of a zero-length array under @inbounds.
        prec10_patch::AbstractVector{<:Real} = Float64[],
        prec30_patch::AbstractVector{<:Real} = Float64[],
        prec60_patch::AbstractVector{<:Real} = Float64[],
        rh30_patch::AbstractVector{<:Real} = Float64[],
        forc_rh_grc::AbstractVector{<:Real} = Float64[],
        forc_wind_grc::AbstractVector{<:Real} = Float64[],
        # Mineral-N inputs (deposition + fixation) — see cn_driver_no_leaching!.
        forc_ndep::AbstractVector{<:Real} = Float64[],
        AnnET::AbstractVector{<:Real} = Float64[],
        nfix_timeconst::Real = 0.0,
        dayspyr::Real = 365.0,
        ndyn_params::Union{NDynamicsParams, Nothing} = nothing,
        h2osoi_vol::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        h2osoi_liq::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        mask_actfirec::AbstractVector{Bool} = falses(length(bounds_col)),
        mask_actfirep::AbstractVector{Bool} = falses(length(bounds_patch)))

    # crop_inst%CropIncrementYear — ported as crop_increment_year! (types/crop.jl);
    # not called from here.

    # C13/C14 cascade activation: run the CIsoFlux* cascade only when the parallel
    # veg state (facade) AND soil state (driver-supplied) are both stood up.
    # cn_driver_no_leaching! additionally gates on driver_config.use_c13/use_c14, so
    # passing the bundles is a no-op unless the tracer is on. Pass nothing otherwise
    # so an unallocated instance is never indexed.
    _c13_on = !isempty(veg.c13_cnveg_carbonstate_inst.cpool_patch) && c13_soilbgc_cs !== nothing
    _c14_on = !isempty(veg.c14_cnveg_carbonstate_inst.cpool_patch) && c14_soilbgc_cs !== nothing

    # CNDriverNoLeaching — already ported
    cn_driver_no_leaching!(veg.driver_config;
        c13_cnveg_cs = _c13_on ? veg.c13_cnveg_carbonstate_inst : nothing,
        c13_cnveg_cf = _c13_on ? veg.c13_cnveg_carbonflux_inst : nothing,
        c13_soilbgc_cs = _c13_on ? c13_soilbgc_cs : nothing,
        c13_soilbgc_cf = _c13_on ? c13_soilbgc_cf : nothing,
        c14_cnveg_cs = _c14_on ? veg.c14_cnveg_carbonstate_inst : nothing,
        c14_cnveg_cf = _c14_on ? veg.c14_cnveg_carbonflux_inst : nothing,
        c14_soilbgc_cs = _c14_on ? c14_soilbgc_cs : nothing,
        c14_soilbgc_cf = _c14_on ? c14_soilbgc_cf : nothing,
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
        cnveg_state=veg.cnveg_state_inst,
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
        patch=patch,
        pftcon_main=pftcon_main,
        crop=crop,
        photosyns=photosyns,
        canopystate=canopystate,
        soilstate=soilstate,
        temperature=temperature,
        water_diag=water_diag,
        waterstate=waterstate,
        gridcell=gridcell,
        is_first_step=is_first_step,
        prec10_patch=prec10_patch,
        prec30_patch=prec30_patch,
        prec60_patch=prec60_patch,
        rh30_patch=rh30_patch,
        forc_rh_grc=forc_rh_grc,
        forc_wind_grc=forc_wind_grc,
        forc_ndep=forc_ndep,
        AnnET=AnnET,
        nfix_timeconst=nfix_timeconst,
        dayspyr=dayspyr,
        ndyn_params=ndyn_params,
        h2osoi_vol=h2osoi_vol,
        h2osoi_liq=h2osoi_liq,
        mask_actfirec=mask_actfirec,
        mask_actfirep=mask_actfirep)

    # CNFireEmisUpdate — genuinely NOT ported (CNFireEmissionsMod: fire trace-gas
    # emissions). The fire model itself is ported (fire_li2014.jl … fire_li2024.jl).

    # CNAnnualUpdate — WIRED. Fortran calls it as the LAST thing in
    # EcosystemDynamicsPreDrainage (CNVegetationFacade.F90:1049). It rolls the
    # tempsum_* running sums into the annsum_* annual means at end-of-year
    # (annsum_npp, annsum_potential_gpp, annmax_retransn, annsum_litfall).
    #
    # This was dead, and it is the second half of the fixation chain:
    # `annsum_npp_col` is what CNNFixation reads on the non-lagged branch, and
    # `annsum_potential_gpp`/`annmax_retransn` scale the plant N demand. Without
    # it the annual means never advance from their restart values.
    if col !== nothing && patch !== nothing
        cn_annual_update!(mask_bgc_soilc, mask_bgc_vegp, bounds_col, bounds_patch,
                          col, patch, veg.cnveg_state_inst, veg.cnveg_carbonflux_inst;
                          dt = dt, days_per_year = dayspyr)
    end

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
- CNPrecisionControl (`cn_precision_control!`, biogeochem/cn_precision_control.jl
  — ported, and run on the live path from `cn_driver_no_leaching!`)
- SoilBiogeochemPrecisionControl (`soil_bgc_precision_control!`,
  biogeochem/decomp_precision_control.jl — ported, not wired here)
- CNDriverSummarizeStates (already ported)
- CNDriverSummarizeFluxes (already ported)
- CNVegStructUpdate (`cn_veg_struct_update!` — ported, and run on the live path
  from `clm_drv!`, clm_driver.jl:2078)

Should only be called if `use_cn` or `use_fates_bgc` is true.

Ported from `EcosystemDynamicsPostDrainage` in `CNVegetationFacade.F90`.
"""
function cn_vegetation_ecosystem_post_drainage!(veg::CNVegetationData;
        mask_bgc_soilc::AbstractVector{Bool},
        mask_bgc_vegp::AbstractVector{Bool},
        mask_allc::AbstractVector{Bool} = mask_bgc_soilc,
        mask_actfirec::AbstractVector{Bool} = falses(0),
        mask_actfirep::AbstractVector{Bool} = falses(0),
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
        soilbgc_nf::SoilBiogeochemNitrogenFluxData,
        patch_itype::Union{AbstractVector{<:Integer},Nothing}=nothing,
        # N-leaching inputs (see cn_driver_leaching!); omitted => no leaching loss.
        nleach_params::Union{NLeachingParams, Nothing} = nothing,
        h2osoi_liq::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        qflx_drain::Union{AbstractVector{<:Real}, Nothing} = nothing,
        qflx_surf::Union{AbstractVector{<:Real}, Nothing} = nothing,
        col_dz::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        zisoi::Union{AbstractVector{<:Real}, Nothing} = nothing,
        nlevsoi::Int = 0,
        # Column NPP / lagged NPP update (feeds next step's CNNFixation).
        col::Union{ColumnData, Nothing} = nothing,
        patch::Union{PatchData, Nothing} = nothing,
        nfix_timeconst::Real = 0.0,
        # CNDriverSummarizeFluxes: the soil-BGC flux summaries + the column/gridcell
        # halves of the CNVeg summaries need the decomp cascade, the subgrid maps and
        # the product pools. Omitted => those summaries are skipped.
        bounds_grc::UnitRange{Int} = 1:0,
        decomp = nothing,
        dzsoi_decomp_vals::Union{AbstractVector{<:Real}, Nothing} = nothing)

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
        cnveg_nf=veg.cnveg_nitrogenflux_inst,
        nleach_params=nleach_params,
        h2osoi_liq=h2osoi_liq,
        qflx_drain=qflx_drain,
        qflx_surf=qflx_surf,
        col_dz=col_dz,
        zisoi=zisoi,
        nlevsoi=nlevsoi)

    # CNPrecisionControl — ported (cn_precision_control!); the live path runs it
    # inside cn_driver_no_leaching!, not here.
    # SoilBiogeochemPrecisionControl — ported (soil_bgc_precision_control!,
    # biogeochem/decomp_precision_control.jl); not called from here.

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
        soilbgc_ns=soilbgc_ns,
        patch_itype=patch_itype)

    # CNDriverSummarizeFluxes
    cn_driver_summarize_fluxes!(veg.driver_config;
        mask_bgc_soilc=mask_bgc_soilc,
        mask_bgc_vegp=mask_bgc_vegp,
        bounds_col=bounds_col,
        bounds_patch=bounds_patch,
        cnveg_cf=veg.cnveg_carbonflux_inst,
        cnveg_nf=veg.cnveg_nitrogenflux_inst,
        soilbgc_cf=soilbgc_cf,
        soilbgc_nf=soilbgc_nf,
        patch_itype=patch_itype,
        col=col,
        patch=patch,
        bounds_grc=bounds_grc,
        decomp=decomp,
        dzsoi_decomp_vals=dzsoi_decomp_vals,
        nlevdecomp=nlevdecomp,
        ndecomp_pools=ndecomp_pools,
        ndecomp_cascade_transitions=ndecomp_cascade_transitions,
        c_products=veg.c_products_inst,
        dt=dt,
        nfix_timeconst=nfix_timeconst)

    # CNVegStructUpdate — ported (cn_veg_struct_update!); the live path runs it
    # from clm_drv! (clm_driver.jl:2078), not here.
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
        mask_bgc_soilc::AbstractVector{Bool},
        bounds_col::UnitRange{Int},
        nstep_since_startup::Int = 1,
        soilbgc_cf::SoilBiogeochemCarbonFluxData,
        soilbgc_nf::SoilBiogeochemNitrogenFluxData,
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData,
        # Subgrid maps + dt, needed by the checks. Omitted => check stays off
        # (a caller that cannot supply them gets the previous no-op behaviour).
        col::Union{ColumnData, Nothing} = nothing,
        grc::Union{GridcellData, Nothing} = nothing,
        bounds_grc::UnitRange{Int} = 1:0,
        dt::Real = 0.0,
        # Opt-in: running the check is what exposes conservation bugs, but it
        # `error()`s on failure, so callers gate it explicitly.
        run_check::Bool = false)

    if nstep_since_startup <= veg.config.skip_steps
        # Skip balance check for first timesteps after startup
        return nothing
    end

    # cn_balance_inst%CBalanceCheck / NBalanceCheck — WIRED. These were dead: the
    # very check that would have caught nitrogen entering and leaving the
    # ecosystem nowhere was itself never called.
    (run_check && col !== nothing && grc !== nothing && dt > 0) || return nothing

    c_balance_check!(veg.cn_balance_inst, soilbgc_cf, soilbgc_cs,
                     veg.cnveg_carbonflux_inst, veg.c_products_inst,
                     col, grc, mask_bgc_soilc, bounds_col, bounds_grc, dt)

    n_balance_check!(veg.cn_balance_inst, soilbgc_nf, soilbgc_ns,
                     veg.cnveg_nitrogenflux_inst, veg.n_products_inst,
                     col, grc, mask_bgc_soilc, bounds_col, bounds_grc, dt;
                     use_nitrif_denitrif = veg.driver_config.use_nitrif_denitrif,
                     use_crop            = veg.driver_config.use_crop,
                     use_fun             = veg.driver_config.use_fun)

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
            # CNDVDriver — ported as cndv_driver! (biogeochem/cndv.jl); the live path
            # runs it from clm_drv! (clm_driver.jl:2626, gated on use_cndv), not here.
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
                                  tlai::AbstractVector{<:Real} = Float64[],
                                  slatop::AbstractVector{<:Real} = Float64[],
                                  froot_leaf::AbstractVector{<:Real} = Float64[],
                                  ivt::AbstractVector{<:Integer} = Int[])
    if veg.config.use_cn
        return veg.cnveg_carbonstate_inst.frootc_patch[bounds_patch]
    else
        FT = eltype(tlai)
        # Small SP-mode estimate: compute on host (inputs gathered if on device), then
        # place the result on the inputs' backend so it composes with device downstream.
        # No-op on the host path (Array-of-Array is identity) → byte-identical.
        tlaiH = tlai isa Array ? tlai : Array(tlai)
        slatopH = slatop isa Array ? slatop : Array(slatop)
        frootH = froot_leaf isa Array ? froot_leaf : Array(froot_leaf)
        ivtH = ivt isa Array ? ivt : Array(ivt)
        result = zeros(FT, length(bounds_patch))
        for (i, p) in enumerate(bounds_patch)
            pft = ivtH[p]
            if slatopH[pft] > 0.0
                result[i] = tlaiH[p] / slatopH[pft] * frootH[pft]
            end
        end
        return _to_backend_like(tlai, FT, result)
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
                                  tlai::AbstractVector{<:Real} = Float64[],
                                  slatop::AbstractVector{<:Real} = Float64[],
                                  stem_leaf::AbstractVector{<:Real} = Float64[],
                                  croot_stem::AbstractVector{<:Real} = Float64[],
                                  ivt::AbstractVector{<:Integer} = Int[])
    if veg.config.use_cn
        return veg.cnveg_carbonstate_inst.livecrootc_patch[bounds_patch]
    else
        FT = eltype(tlai)
        tlaiH = tlai isa Array ? tlai : Array(tlai)
        slatopH = slatop isa Array ? slatop : Array(slatop)
        stemH = stem_leaf isa Array ? stem_leaf : Array(stem_leaf)
        crootH = croot_stem isa Array ? croot_stem : Array(croot_stem)
        ivtH = ivt isa Array ? ivt : Array(ivt)
        result = zeros(FT, length(bounds_patch))
        for (i, p) in enumerate(bounds_patch)
            pft = ivtH[p]
            if slatopH[pft] > 0.0
                result[i] = tlaiH[p] / slatopH[pft] * stemH[pft] * crootH[pft]
            end
        end
        return _to_backend_like(tlai, FT, result)
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
    dc.use_flexiblecn = fc.use_flexiblecn
    dc.use_crop = fc.use_crop
    dc.use_crop_agsys = fc.use_crop_agsys
    dc.use_nitrif_denitrif = fc.use_nitrif_denitrif
    dc.use_fates_bgc = fc.use_fates_bgc
    dc.use_matrixcn = fc.use_matrixcn
    dc.use_soil_matrixcn = fc.use_soil_matrixcn
    dc.dribble_crophrv_xsmrpool_2atm = fc.dribble_crophrv_xsmrpool_2atm

    return nothing
end
