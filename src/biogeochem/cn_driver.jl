# ==========================================================================
# Ported from: src/biogeochem/CNDriverMod.F90
# CN driver: ecosystem dynamics orchestrator for phenology, vegetation,
# decomposition, nitrogen cycling, fire, and state updates.
#
# Public functions:
#   cn_driver_init!              — Initialization
#   cn_driver_no_leaching!       — Main CN driver (before N leaching)
#   cn_driver_leaching!          — N leaching and state updates
#   cn_driver_summarize_states!  — Summarize state variables
#   cn_driver_summarize_fluxes!  — Summarize flux variables
# ==========================================================================

# ---------------------------------------------------------------------------
# CNDriverConfig — control flags for the CN driver
# Holds all configuration flags referenced throughout the driver.
# ---------------------------------------------------------------------------

"""
    CNDriverConfig

Configuration flags for the CN driver module. Aggregates control variables
from `clm_varctl`, `CNSharedParamsMod`, and `SoilBiogeochemDecompCascadeConType`
that determine which code paths are active in the driver.

Ported from module-level `use` statements in `CNDriverMod.F90`.
"""
Base.@kwdef mutable struct CNDriverConfig
    use_c13::Bool = false
    use_c14::Bool = false
    use_cn::Bool = false
    use_fun::Bool = false
    use_flexiblecn::Bool = false     # flexible leaf C:N (FUN cost adjustment)
    carbon_resp_opt::Int = 0         # 1 => FlexibleCN high-C:N tissue C turnover
    use_crop::Bool = false
    use_crop_agsys::Bool = false
    use_nitrif_denitrif::Bool = false
    use_nguardrail::Bool = false
    use_fates::Bool = false
    use_fates_bgc::Bool = false
    use_matrixcn::Bool = false
    use_soil_matrixcn::Bool = false
    decomp_method::Int = 1           # 1=century_decomp, 2=mimics_decomp
    dribble_crophrv_xsmrpool_2atm::Bool = false
    do_harvest::Bool = false
    do_grossunrep::Bool = false
    # CN fire method (mirrors the `cnfire_method` namelist string mapped to a
    # symbol by `cnfire_method_symbol`). Default `:nofire` keeps the CN driver
    # byte-identical to the historical no-fire path; the live Li-family fire
    # (cnfire_area! → cnfire_fluxes_dispatch! → fire fluxes feeding the fire
    # state update) runs only when this is set to :li2014/:li2016/:li2021/:li2024
    # AND the fire-input bundle is supplied to cn_driver_no_leaching!.
    cnfire_method::Symbol = :nofire
end

# Decomposition method constants (from SoilBiogeochemDecompCascadeConType)
const CENTURY_DECOMP = 1
const MIMICS_DECOMP  = 2

# ---------------------------------------------------------------------------
# cn_driver_init! — Initialization
# Ported from CNDriverInit in CNDriverMod.F90
# ---------------------------------------------------------------------------

"""
    cn_driver_init!(config; kwargs...)

Initialize the CN ecosystem dynamics.

Fortran's `CNDriverInit` calls `SoilBiogeochemCompetitionInit`, `CNPhenologyInit`
and `FireInit`. All three ARE ported — they are just initialized elsewhere in
CLM.jl, so this hook has nothing left to do:
- `SoilBiogeochemCompetitionInit` → `soil_bgc_competition_init!`
  (`biogeochem/decomp_competition.jl`), called from `clm_initialize.jl:366`.
- `CNPhenologyInit` → `cn_phenology_init!` (`biogeochem/phenology.jl`).
- `FireInit` → the Li-family fire modules (`biogeochem/fire_base.jl`,
  `fire_li2014.jl` … `fire_li2024.jl`, dispatched by `fire_factory.jl`). Their
  params/state (`CNFireParams`, `CNFireBaseData`, …) are constructed by the
  caller and passed into `cn_driver_no_leaching!`, so there is no separate
  Fortran-style init hook.

This function is therefore a no-op retained for call-structure traceability.

Ported from `CNDriverInit` in `CNDriverMod.F90`.
"""
function cn_driver_init!(config::CNDriverConfig;
                          bounds::UnitRange{Int} = 1:0)
    # No-op: see the docstring — competition/phenology/fire are all initialized
    # on their own paths (clm_initialize.jl / caller-supplied fire bundles).
    if config.use_cn
        # (CNPhenologyInit / FireInit equivalents run outside this hook.)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# cn_driver_no_leaching! — Main CN driver (before N leaching)
# Ported from CNDriverNoLeaching in CNDriverMod.F90
# ---------------------------------------------------------------------------

"""
    cn_driver_no_leaching!(config; kwargs...)

The core CN code. Calculates fluxes for maintenance respiration,
decomposition, allocation, phenology, and growth respiration.
These routines happen on the radiation time step so that canopy structure
stays synchronized with albedo calculations.

The function orchestrates calls to all CN subroutines in the correct order.
Nutrient competition, fire, vertical transport and the isotope fluxes are all
ported now and wired here (fire behind `config.cnfire_method`, default
`:nofire`; the isotope/matrix paths behind their own `use_*` flags), so the
"documented as placeholder" comments this docstring used to carry are gone.
Anything still marked as not wired below names the specific missing piece.

Ported from `CNDriverNoLeaching` in `CNDriverMod.F90`.
"""
function cn_driver_no_leaching!(
        config::CNDriverConfig;
        # Masks (replace Fortran filter arrays)
        mask_bgc_soilc::AbstractVector{Bool},
        mask_bgc_vegp::AbstractVector{Bool},
        mask_pcropp::AbstractVector{Bool} = falses(0),
        mask_soilnopcropp::AbstractVector{Bool} = falses(0),
        mask_exposedvegp::AbstractVector{Bool} = falses(0),
        mask_noexposedvegp::AbstractVector{Bool} = falses(0),
        # Bounds
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        # Decomposition dimensions
        nlevdecomp::Int,
        nlevdecomp_full::Int = nlevdecomp,
        ndecomp_pools::Int,
        ndecomp_cascade_transitions::Int,
        i_litr_min::Int,
        i_litr_max::Int,
        i_cwd::Int,
        # Vegetation parameters
        npcropmin::Int = 17,
        nrepr::Int = 1,
        # Patch/column metadata (from PatchData, PftCon, CropData, etc.)
        patch_column::AbstractVector{<:Integer},
        ivt::AbstractVector{<:Integer},
        woody::AbstractVector{<:Real},
        harvdate::AbstractVector{<:Integer},
        col_is_fates::AbstractVector{Bool},
        cascade_donor_pool::AbstractVector{<:Integer},
        cascade_receiver_pool::AbstractVector{<:Integer},
        # Time step
        dt::Real,
        # Carbon/nitrogen state and flux data structures
        cnveg_cs::CNVegCarbonStateData,
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_ns::CNVegNitrogenStateData,
        cnveg_nf::CNVegNitrogenFluxData,
        cnveg_state::Union{CNVegStateData, Nothing} = nothing,
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_cf::SoilBiogeochemCarbonFluxData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData,
        soilbgc_nf::SoilBiogeochemNitrogenFluxData,
        soilbgc_state::SoilBiogeochemStateData,
        # Carbon-isotope (C13/C14) tracer instances. Optional: only used when the
        # matching use_c13/use_c14 flag is set AND the instance bundle is supplied;
        # the CIsoFlux1/2/2h/2g/3 cascades + isotope litter transport are then run.
        # Default (both nothing) → no isotope work → byte-identical to the bulk path.
        c13_cnveg_cs::Union{CNVegCarbonStateData, Nothing} = nothing,
        c13_cnveg_cf::Union{CNVegCarbonFluxData, Nothing} = nothing,
        c13_soilbgc_cs::Union{SoilBiogeochemCarbonStateData, Nothing} = nothing,
        c13_soilbgc_cf::Union{SoilBiogeochemCarbonFluxData, Nothing} = nothing,
        c14_cnveg_cs::Union{CNVegCarbonStateData, Nothing} = nothing,
        c14_cnveg_cf::Union{CNVegCarbonFluxData, Nothing} = nothing,
        c14_soilbgc_cs::Union{SoilBiogeochemCarbonStateData, Nothing} = nothing,
        c14_soilbgc_cf::Union{SoilBiogeochemCarbonFluxData, Nothing} = nothing,
        # Decomposition cascade and params (for rate constants + potential + decomp)
        cascade_con::Union{DecompCascadeConData, Nothing} = nothing,
        decomp_bgc_state::Union{DecompBGCState, Nothing} = nothing,
        decomp_bgc_params::Union{DecompBGCParams, Nothing} = nothing,
        cn_shared_params::Union{CNSharedParamsData, Nothing} = nothing,
        decomp_params::Union{DecompParams, Nothing} = nothing,
        competition_state::Union{SoilBGCCompetitionState, Nothing} = nothing,
        competition_params::Union{SoilBGCCompetitionParams, Nothing} = nothing,
        litter_params::Union{LitterVertTranspParams, Nothing} = nothing,
        # Physical state for decomposition
        t_soisno::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        soilpsi::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        col::Union{ColumnData, Nothing} = nothing,
        grc::Union{GridcellData, Nothing} = nothing,
        active_layer::Union{ActiveLayerData, Nothing} = nothing,
        # Vertical grid info
        dzsoi_decomp::Union{AbstractVector{<:Real}, Nothing} = nothing,
        zsoi_vals::Union{AbstractVector{<:Real}, Nothing} = nothing,
        zisoi_vals::Union{AbstractVector{<:Real}, Nothing} = nothing,
        # Vegetation-flux inputs (for the veg-CN flux chain: mresp, gpp/alloc, …).
        # Optional so SP / decomp-only callers are unaffected; the veg-flux routines
        # are wired only when these are supplied.
        patch::Union{PatchData, Nothing} = nothing,
        pftcon_main::Union{Any, Nothing} = nothing,
        # Transient wood-harvest + gross-unrepresented-landcover mortality (CNHarvest /
        # CNGrossUnrep). Wired only when the corresponding state (annual rates, read by
        # the dyn_subgrid driver) + pft params are supplied and config.do_harvest /
        # do_grossunrep are on.
        dynharvest_state = nothing,
        dyngrossunrep_state = nothing,
        pftcon_grossunrep = nothing,
        # Caller-owned persistent soil-matrix workspace (CNSoilMatrixState). When
        # supplied it is reused across steps so the sparse index memoization survives;
        # a fresh one is built each step otherwise (correct but slower).
        soil_matrix_state = nothing,
        # Caller-owned persistent veg-matrix structure workspaces (CNVegMatrixSolveState).
        # Supplied for a GPU matrix solve (built once on the host, reused device-side); a
        # fresh rebuild each step otherwise. See cn_veg_matrix_solve_c!/_n!'s `state`.
        veg_c_solve_state = nothing,
        veg_n_solve_state = nothing,
        crop::Union{CropData, Nothing} = nothing,
        photosyns::Union{PhotosynthesisData, Nothing} = nothing,
        canopystate::Union{CanopyStateData, Nothing} = nothing,
        soilstate::Union{SoilStateData, Nothing} = nothing,
        temperature::Union{TemperatureData, Nothing} = nothing,
        water_diag::Union{WaterDiagnosticBulkData, Nothing} = nothing,
        waterstate::Union{WaterStateData, Nothing} = nothing,   # for FUN (h2osoi_liq_col)
        gridcell::Union{GridcellData, Nothing} = nothing,
        is_first_step::Bool = false,
        # ---- Mineral-N INPUTS (CNNDynamicsMod) ----------------------------
        # forc_ndep: gridcell atmospheric N deposition (gN/m2/s), from the ndep
        # stream via atm2lnd.forc_ndep_grc. Empty => n_deposition! is skipped and
        # ndep_to_sminn stays at whatever it was (0), i.e. the historical
        # behaviour. AnnET: column annual-mean ET (mm/s), only used by the FUN
        # free-living fixation. nfix_timeconst: clm_varctl (10 d when
        # use_nitrif_denitrif, else 0) — selects the lagged-NPP vs annual-mean-NPP
        # branch of CNNFixation.
        forc_ndep::AbstractVector{<:Real} = Float64[],
        AnnET::AbstractVector{<:Real} = Float64[],
        nfix_timeconst::Real = 0.0,
        dayspyr::Real = 365.0,
        ndyn_params::Union{NDynamicsParams, Nothing} = nothing,
        # Soil-water state (soil-layer slices) for nitrif/denitrif anoxia.
        h2osoi_vol::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        h2osoi_liq::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        # ---- CN fire inputs (Li2014-family) -------------------------------
        # The live fire path (cnfire_area! burned-area + cnfire_fluxes_dispatch!
        # C/N fluxes) runs only when config.cnfire_method !== :nofire AND this
        # full bundle is supplied; otherwise the driver is byte-identical to the
        # historical no-fire placeholder. All defaulted to `nothing` so SP /
        # decomp-only / no-fire callers are unaffected.
        fire_data::Union{CNFireBaseData, Nothing} = nothing,
        fire_li2014::Union{CNFireLi2014Data, Nothing} = nothing,
        cnfire_const::Union{CNFireConstData, Nothing} = nothing,
        cnfire_params::Union{CNFireParams, Nothing} = nothing,
        pftcon_fire::Union{PftConFireBase, Nothing} = nothing,
        pftcon_fire_li2014::Union{PftConFireLi2014, Nothing} = nothing,
        dgvs_fire::Union{DgvsFireData, Nothing} = nothing,
        # Exposed/non-exposed-veg patch masks for the burned-area calc (default to
        # the BGC veg mask / its complement when unset).
        mask_exposedvegp_fire::Union{AbstractVector{Bool}, Nothing} = nothing,
        mask_noexposedvegp_fire::Union{AbstractVector{Bool}, Nothing} = nothing,
        # Fire forcing / running-mean / water inputs (gridcell/column/patch).
        forc_rh_grc::AbstractVector{<:Real} = Float64[],
        forc_wind_grc::AbstractVector{<:Real} = Float64[],
        forc_t_fire_col::AbstractVector{<:Real} = Float64[],
        forc_rain_fire_col::AbstractVector{<:Real} = Float64[],
        forc_snow_fire_col::AbstractVector{<:Real} = Float64[],
        prec60_patch::AbstractVector{<:Real} = Float64[],
        prec10_patch::AbstractVector{<:Real} = Float64[],
        prec30_patch::AbstractVector{<:Real} = Float64[],
        rh30_patch::AbstractVector{<:Real} = Float64[],
        fsat_fire_col::AbstractVector{<:Real} = Float64[],
        wf_fire_col::AbstractVector{<:Real} = Float64[],
        wf2_fire_col::AbstractVector{<:Real} = Float64[],
        # Fire date/step scalars.
        fire_kmo::Int = 1, fire_kda::Int = 1, fire_mcsec::Int = 0, fire_nstep::Int = 1,
        nlevgrnd_fire::Int = 10,
        transient_landcover::Bool = false,
        # Output: fire masks (populated by fire routines)
        mask_actfirec::AbstractVector{Bool} = falses(length(bounds_col)),
        mask_actfirep::AbstractVector{Bool} = falses(length(bounds_patch)))

    num_bgc_vegp = count(mask_bgc_vegp)

    # --------------------------------------------------
    # Carbon-isotope (C13/C14) cascade gating
    # --------------------------------------------------
    # The CIsoFlux* cascades scale each bulk C flux into its isotopic counterpart
    # using the donor-pool ratio (iso_state/tot_state). They run only when the flag
    # is set AND the full isotope-instance bundle (veg + soil C state/flux) is
    # supplied. Build a list of (isotope-tag, instances) to drive at each hook.
    _iso_active = Tuple{String, CNVegCarbonStateData, CNVegCarbonFluxData,
                        SoilBiogeochemCarbonStateData, SoilBiogeochemCarbonFluxData}[]
    if config.use_c13 && c13_cnveg_cs !== nothing && c13_cnveg_cf !== nothing &&
       c13_soilbgc_cs !== nothing && c13_soilbgc_cf !== nothing
        push!(_iso_active, ("c13", c13_cnveg_cs, c13_cnveg_cf, c13_soilbgc_cs, c13_soilbgc_cf))
    end
    if config.use_c14 && c14_cnveg_cs !== nothing && c14_cnveg_cf !== nothing &&
       c14_soilbgc_cs !== nothing && c14_soilbgc_cf !== nothing
        push!(_iso_active, ("c14", c14_cnveg_cs, c14_cnveg_cf, c14_soilbgc_cs, c14_soilbgc_cf))
    end
    _have_iso = !isempty(_iso_active)
    # Patch/pft metadata used by the patch→column scatters inside the cascades.
    _iso_pwtcol = (patch !== nothing) ? patch.wtcol : Float64[]
    _iso_lf_f   = (pftcon_main !== nothing) ? Float64.(pftcon_main.lf_f) : Matrix{Float64}(undef, 0, 0)
    _iso_fr_f   = (pftcon_main !== nothing) ? Float64.(pftcon_main.fr_f) : Matrix{Float64}(undef, 0, 0)
    _iso_pcol   = collect(Int, patch_column)
    _iso_pivt   = collect(Int, ivt)

    # Convert the timestep to the working precision of the state once, so every
    # sub-module kernel receives an FT scalar (Metal rejects Float64 args). On the
    # CPU eltype is Float64 so dt is unchanged (byte-identical).
    dt = eltype(cnveg_cs.cpool_patch)(dt)

    # --------------------------------------------------
    # Zero the column-level C and N fluxes
    # --------------------------------------------------
    # In Fortran: soilbiogeochem_carbonflux_inst%SetValues(...)
    _zero_soilbgc_cflux!(soilbgc_cf; mask=mask_bgc_soilc, bounds=bounds_col)

    if num_bgc_vegp > 0
        _zero_cnveg_cflux!(cnveg_cf; mask=mask_bgc_vegp, bounds=bounds_patch,
                           mask_col=mask_bgc_soilc, bounds_col=bounds_col)
        _zero_cnveg_nflux!(cnveg_nf; mask=mask_bgc_vegp, bounds=bounds_patch,
                           mask_col=mask_bgc_soilc, bounds_col=bounds_col)
    end

    _zero_soilbgc_nflux!(soilbgc_nf; mask=mask_bgc_soilc, bounds=bounds_col)

    # --------------------------------------------------
    # Nitrogen Deposition, Fixation and Respiration
    # --------------------------------------------------
    # WIRED. Fortran CNDriverNoLeaching lines 301-329: CNNDeposition, then EITHER
    # CNFreeLivingFixation (use_fun) OR CNNFixation, then the crop pair
    # (CNNFert/CNSoyfix). These are the ONLY inputs of N to the ecosystem; before
    # this they were dead ports (definitions with no call site), so ndep_to_sminn
    # and nfix_to_sminn were identically zero and nitrogen entered nowhere.
    #
    # Each call is skipped when its input is absent, so a caller that does not
    # supply the forcing gets exactly the previous (zero-input) behaviour.

    # CNNDeposition — atmospheric N deposition -> soil mineral N.
    if !isempty(forc_ndep) && col !== nothing
        n_deposition!(soilbgc_nf;
            forc_ndep    = forc_ndep,
            col_gridcell = col.gridcell,
            bounds       = bounds_col)
    end

    # CNFreeLivingFixation (FUN) / CNNFixation (non-FUN) — mutually exclusive in
    # Fortran, and they write DIFFERENT fields (ffix_to_sminn vs nfix_to_sminn),
    # which SoilBiogeochemNStateUpdate1 then selects between on use_fun.
    _ndp = ndyn_params === nothing ? NDynamicsParams() : ndyn_params
    if config.use_fun
        # AnnET is a 365-day running-mean ET ACCUMULATOR. Where the accumulator is
        # not yet maintained it sits at its NaN allocation, and feeding that in would
        # propagate NaN straight into the mineral-N pool. Require a finite value
        # rather than silently poisoning the N budget: on the parity/injection path
        # AnnET comes from the Fortran restart (AnnET_VALUE), and in a free run it
        # comes from the accumulator once that is live.
        _annet_ok = !isempty(AnnET) &&
                    all(c -> !mask_bgc_soilc[c] || isfinite(AnnET[c]), bounds_col)
        if _annet_ok
            n_free_living_fixation!(soilbgc_nf, _ndp;
                mask_soilc = mask_bgc_soilc,
                bounds     = bounds_col,
                AnnET      = AnnET,
                dayspyr    = dayspyr)
        end
    else
        # CNNFixation reads the LAGGED NPP when nfix_timeconst ∈ (0,500) (the CLM5
        # default, 10 d) and the ANNUAL-MEAN NPP otherwise. Fortran guards the
        # former with its own SPVAL check; the latter it does not guard, because
        # Fortran's annsum_npp is always initialised (to 0 in InitCold). Julia
        # allocates it as NaN, so an un-summarised annsum_npp would push NaN
        # straight into the mineral-N pool. Require a finite input for whichever
        # branch will actually be taken — same rule as AnnET above: never poison
        # the N budget with an accumulator that has not been computed yet.
        _npp_src = (nfix_timeconst > 0.0 && nfix_timeconst < 500.0) ?
                   cnveg_cf.lag_npp_col : cnveg_cf.annsum_npp_col
        _npp_ok = !isempty(_npp_src) &&
                  all(c -> !mask_bgc_soilc[c] || isfinite(_npp_src[c]), bounds_col)
        if _npp_ok
            n_fixation!(soilbgc_nf, cnveg_cf;
                mask_soilc     = mask_bgc_soilc,
                bounds         = bounds_col,
                col_is_fates   = col_is_fates,
                dayspyr        = dayspyr,
                nfix_timeconst = nfix_timeconst,
                use_fun        = config.use_fun)
        end
    end

    # CNNFert / CNSoyfix — crop-only (Fortran gates on use_crop). NOT wired: Bow
    # has no crop CFT, so there is no Fortran reference here to validate against,
    # and #218's rule stands — do not wire blind. Both are ported and unit-tested
    # (n_fert!, n_soyfix! in n_dynamics.jl); they need a crop-active Fortran
    # reference run before they go live. See docs/N_CYCLE_PARITY.md.

    # Maintenance respiration (CNMResp) — WIRED (veg-flux chain). Runs when the
    # veg-flux inputs (photosyns/canopystate/soilstate/temperature/patch/pftcon)
    # are supplied; produces cnveg_cf.{leaf,froot,livestem,livecroot,reproductive}_mr.
    if num_bgc_vegp > 0 && cn_shared_params !== nothing && pftcon_main !== nothing &&
       patch !== nothing && canopystate !== nothing && soilstate !== nothing &&
       temperature !== nothing && photosyns !== nothing
        _mrp = MaintRespParams(); maint_resp_read_params!(_mrp; br_root = cn_shared_params.br_root)
        _pftmr = PftConMaintResp{Float64}(woody = Float64.(pftcon_main.woody))
        cn_mresp!(mask_bgc_soilc, mask_bgc_vegp, bounds_col, bounds_patch,
                  _mrp, cn_shared_params, _pftmr,
                  patch, canopystate, soilstate, temperature, photosyns,
                  cnveg_cf, cnveg_ns;
                  npcropmin = npcropmin, nrepr = nrepr)
    end

    # GPP, maintenance-resp totals, and available C (calc_gpp_mr_availc!) — WIRED.
    # Produces cnveg_cf.{availc, gpp_before_downreg, psnsun_to_cpool, …} consumed by
    # plant N demand + allocation. Builds PftConAllocation as a subset of the main
    # pftcon (all fields present); reused below by calc_allometry!.
    _pfta = (pftcon_main === nothing) ? nothing : PftConAllocation{Float64}(
        woody = Float64.(pftcon_main.woody), froot_leaf = Float64.(pftcon_main.froot_leaf),
        croot_stem = Float64.(pftcon_main.croot_stem), stem_leaf = Float64.(pftcon_main.stem_leaf),
        flivewd = Float64.(pftcon_main.flivewd), leafcn = Float64.(pftcon_main.leafcn),
        frootcn = Float64.(pftcon_main.frootcn), livewdcn = Float64.(pftcon_main.livewdcn),
        deadwdcn = Float64.(pftcon_main.deadwdcn), graincn = Float64.(pftcon_main.graincn),
        grperc = Float64.(pftcon_main.grperc))
    if num_bgc_vegp > 0 && cn_shared_params !== nothing && _pfta !== nothing &&
       patch !== nothing && crop !== nothing && photosyns !== nothing && canopystate !== nothing
        _ap = AllocationParams()
        calc_gpp_mr_availc!(mask_bgc_vegp, bounds_patch, _ap, _pfta, cn_shared_params,
                  patch, crop, photosyns, canopystate, cnveg_cs, cnveg_cf;
                  use_c13 = config.use_c13, use_c14 = config.use_c14,
                  npcropmin = npcropmin, nrepr = nrepr)
    end

    # Crop allocation fractions (calc_crop_allocation_fractions!) — crops only;
    # Bow has no crop patches (mask_pcropp empty) so this is a no-op here. (Wire
    # later when crop infra is needed.)

    # C/N allometry (calc_allometry!) — WIRED. Produces cnveg_state.{c_allometry,
    # n_allometry} used by plant N demand + allocation.
    if num_bgc_vegp > 0 && cn_shared_params !== nothing && _pfta !== nothing &&
       patch !== nothing && cnveg_state !== nothing
        calc_allometry!(mask_bgc_vegp, bounds_patch, _pfta, cn_shared_params,
                  patch, cnveg_cf, cnveg_state;
                  npcropmin = npcropmin, nrepr = nrepr)
    end

    # Plant nitrogen demand (calc_plant_nutrient_demand!) — WIRED. Non-crop call
    # over the veg-patch filter (Bow has no prognostic crops). Produces
    # cnveg_nf.plant_ndemand_patch, then aggregated patch→column into
    # soilbgc_state.plant_ndemand_col so the (already-wired) soil_bgc_competition!
    # below computes a real fpg_col. PftConNutrientCompetition is a subset of the
    # main pftcon (all fields present).
    _pftnc = (pftcon_main === nothing) ? nothing : PftConNutrientCompetition{Float64}(
        woody = Float64.(pftcon_main.woody), froot_leaf = Float64.(pftcon_main.froot_leaf),
        croot_stem = Float64.(pftcon_main.croot_stem), stem_leaf = Float64.(pftcon_main.stem_leaf),
        flivewd = Float64.(pftcon_main.flivewd), leafcn = Float64.(pftcon_main.leafcn),
        frootcn = Float64.(pftcon_main.frootcn), livewdcn = Float64.(pftcon_main.livewdcn),
        deadwdcn = Float64.(pftcon_main.deadwdcn), fcur = Float64.(pftcon_main.fcur),
        graincn = Float64.(pftcon_main.graincn), grperc = Float64.(pftcon_main.grperc),
        grpnow = Float64.(pftcon_main.grpnow), fleafcn = Float64.(pftcon_main.fleafcn),
        ffrootcn = Float64.(pftcon_main.ffrootcn), fstemcn = Float64.(pftcon_main.fstemcn),
        astemf = Float64.(pftcon_main.astemf), season_decid = Float64.(pftcon_main.season_decid),
        stress_decid = Float64.(pftcon_main.stress_decid), evergreen = Float64.(pftcon_main.evergreen))
    if num_bgc_vegp > 0 && cn_shared_params !== nothing && _pftnc !== nothing &&
       patch !== nothing && crop !== nothing && cnveg_state !== nothing
        # FlexibleCN dispatches to the parallel N-demand method (Michaelis-Menten
        # uptake kinetics + flexible leaf C:N) when use_flexiblecn is set AND the
        # extra soilbgc/canopystate inputs it needs are available; otherwise the
        # CLM4.5 default method runs (default path is byte-identical).
        if config.use_flexiblecn && canopystate !== nothing && dzsoi_decomp !== nothing
            calc_plant_nutrient_demand_flexiblecn!(mask_bgc_vegp, bounds_patch, false,
                _pftnc, cn_shared_params, patch, crop, canopystate, cnveg_state,
                cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf, soilbgc_ns, soilbgc_cf;
                dt = dt, dzsoi_decomp = dzsoi_decomp, nlevdecomp = nlevdecomp,
                npcropmin = npcropmin, nrepr = nrepr)
        else
            calc_plant_nutrient_demand!(mask_bgc_vegp, bounds_patch, false,
                _pftnc, cn_shared_params, patch, crop, cnveg_state,
                cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
                dt = dt, npcropmin = npcropmin, nrepr = nrepr)
        end
        # patch→column "unity" aggregation of plant N demand for competition
        _Tnd = eltype(soilbgc_state.plant_ndemand_col)
        for c in bounds_col
            mask_bgc_soilc[c] && (soilbgc_state.plant_ndemand_col[c] = zero(_Tnd))
        end
        for p in bounds_patch
            mask_bgc_vegp[p] || continue
            w = patch.wtcol[p]; isfinite(w) || continue
            cc = patch_column[p]
            soilbgc_state.plant_ndemand_col[cc] += cnveg_nf.plant_ndemand_patch[p] * w
        end
    end

    # --------------------------------------------------
    # Soil Biogeochemistry
    # --------------------------------------------------
    # Decomposition rate constants + potential + competition + actual decomp
    _has_decomp = (cascade_con !== nothing && decomp_bgc_state !== nothing &&
                   decomp_bgc_params !== nothing && cn_shared_params !== nothing &&
                   t_soisno !== nothing && soilpsi !== nothing &&
                   dzsoi_decomp !== nothing && zsoi_vals !== nothing)

    nc = length(bounds_col)

    # Working arrays for the decomposition chain — device-resident (similar() off a
    # soil-BGC state field) + working precision FT, so the decomp kernels receive
    # device FT scratch rather than host Float64 (which would force the CPU backend).
    _FT = eltype(soilbgc_cs.decomp_cpools_vr_col)
    _z3(d3) = fill!(similar(soilbgc_cs.decomp_cpools_vr_col, _FT, nc, nlevdecomp, d3), zero(_FT))
    cn_decomp_pools       = _z3(ndecomp_pools)
    p_decomp_cpool_loss   = _z3(ndecomp_cascade_transitions)
    p_decomp_cn_gain      = _z3(ndecomp_pools)
    pmnf_decomp_cascade   = _z3(ndecomp_cascade_transitions)
    p_decomp_npool_to_din = _z3(ndecomp_cascade_transitions)

    # Soil matrix-CN turnover diagonal Ksoil (= decomp_k_col·dt, fpi-corrected inside
    # soil_biogeochem_decomp!). Allocated only in matrix mode; handed to cn_soil_matrix!.
    Ksoil_mat = config.use_soil_matrixcn ?
        fill!(similar(soilbgc_cs.decomp_cpools_vr_col, _FT, nc, ndecomp_pools * nlevdecomp), zero(_FT)) :
        nothing

    if _has_decomp
        # 1. Decomposition rate constants (T/moisture-dependent)
        decomp_rate_constants_bgc!(soilbgc_cf,
            decomp_bgc_state, decomp_bgc_params, cn_shared_params, cascade_con;
            mask_bgc_soilc=mask_bgc_soilc,
            bounds=bounds_col,
            nlevdecomp=nlevdecomp,
            t_soisno=t_soisno,
            soilpsi=soilpsi,
            days_per_year=365.0,
            dt=dt,
            zsoi_vals=zsoi_vals,
            # Materialize the soil-layer slice as a CONTIGUOUS array (not a @view)
            # so it is a valid device kernel arg; the empty fallback is device-resident
            # (similar()) rather than a host Matrix{Float64} (the single-level decomp
            # path, nlevdecomp==1, indexes col_dz on the device).
            col_dz=(col !== nothing ? col.dz[:, (varpar.nlevsno+1):end] :
                    similar(soilbgc_cs.decomp_cpools_vr_col, _FT, length(bounds_col), 0)))

        # 2. Potential decomposition and mineral N flux
        soil_bgc_potential!(soilbgc_cf, soilbgc_cs, soilbgc_nf, soilbgc_ns,
            soilbgc_state, cascade_con;
            mask_bgc_soilc=mask_bgc_soilc,
            bounds=bounds_col,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            ndecomp_cascade_transitions=ndecomp_cascade_transitions,
            cn_decomp_pools=cn_decomp_pools,
            p_decomp_cpool_loss=p_decomp_cpool_loss,
            p_decomp_cn_gain=p_decomp_cn_gain,
            pmnf_decomp_cascade=pmnf_decomp_cascade,
            p_decomp_npool_to_din=p_decomp_npool_to_din)

        # 2b. Nitrification / denitrification (SoilBiogeochemNitrifDenitrif) — WIRED.
        # Produces f_nit/f_denit/pot_f_nit (+ NO3/NH4 immob/uptake demand) consumed
        # by the competition's nitrif split and the smin nh4/no3 state update. Runs
        # BEFORE competition, per CNDriverNoLeaching (Fortran order, line 375).
        if config.use_nitrif_denitrif && soilstate !== nothing && cn_shared_params !== nothing &&
           soilpsi !== nothing && t_soisno !== nothing && h2osoi_vol !== nothing &&
           h2osoi_liq !== nothing && col !== nothing
            _ndp = NitrifDenitrifParams()
            nitrif_denitrif_read_params!(_ndp;
                k_nitr_max_perday=0.1, surface_tension_water=0.073,
                rij_kro_a=1.5e-10, rij_kro_alpha=1.26, rij_kro_beta=0.6,
                rij_kro_gamma=0.6, rij_kro_delta=0.85,
                denitrif_respiration_coefficient=0.1, denitrif_respiration_exponent=1.3,
                denitrif_nitrateconc_coefficient=1.15, denitrif_nitrateconc_exponent=0.57,
                om_frac_sf=1.0)
            nitrif_denitrif!(soilbgc_nf, soilbgc_ns, soilbgc_cf, _ndp, cn_shared_params;
                mask_bgc_soilc=mask_bgc_soilc, bounds=bounds_col, nlevdecomp=nlevdecomp,
                watsat=soilstate.watsat_col, watfc=soilstate.watfc_col,
                bd=soilstate.bd_col, bsw=soilstate.bsw_col, cellorg=soilstate.cellorg_col,
                sucsat=soilstate.sucsat_col, soilpsi=Array(soilpsi),
                h2osoi_vol=Array(h2osoi_vol), h2osoi_liq=Array(h2osoi_liq),
                t_soisno=Array(t_soisno),   # _NitrifIn groups inputs as Matrix; materialize views
                col_dz=Array(col.dz[:, (varpar.nlevsno+1):end]), use_lch4=false)
        end

        # 2c. FUN (Fixation & Uptake of Nitrogen) hook. When use_fun, the plant's
        # actual soil-N uptake is a carbon-cost optimization, not availability-
        # weighted. cnfun! runs INSIDE soil_bgc_competition! (Fortran calls CNFUN
        # inside SoilBiogeochemCompetition): after the per-layer loop sets the
        # *offered* N, the hook computes the cost-based actual uptake
        # (sminn_to_plant_fun_{nh4,no3}_vr_patch) and p2c's it to column.
        _fun_hook = nothing
        if config.use_fun && cnveg_state !== nothing && patch !== nothing &&
           pftcon_main !== nothing && soilstate !== nothing && temperature !== nothing &&
           canopystate !== nothing && h2osoi_liq !== nothing
            _pftcon_fun = pftcon_fun_from(pftcon_main)
            _fun_params = FUNParams()
            # Cold-start safety: leafcn_offset is a FUN prognostic (leaf C:N target).
            # When started from a Fortran restart it is injected; on a cold start it
            # is NaN/SPVAL — seed it with the base leaf C:N (Fortran InitCold default)
            # so FUN's litterfall-N / retranslocation stay finite.
            for p in bounds_patch
                mask_bgc_vegp[p] || continue
                lo = cnveg_state.leafcn_offset_patch[p]
                # Seed when unset (cold start → NaN), SPVAL, or non-positive (a
                # zeroed leaf C:N divides to NaN inside FUN). leaf C:N is always > 0.
                if !isfinite(lo) || lo <= 0.0 || lo > 1.0e30
                    _vt = patch.itype[p] + 1   # local PFT index — NOT the ivt kwarg vector
                    cnveg_state.leafcn_offset_patch[p] = _pftcon_fun.leafcn[_vt]
                end
            end
            # cnfun! reads only waterstate.h2osoi_liq_col, indexed by soil layer
            # 1:nlevdecomp. The driver's h2osoi_liq is already the soil slice
            # (snow padding stripped), so wrap it in a minimal WaterStateData.
            _ws_fun = WaterStateData()
            _ws_fun.h2osoi_liq_col = h2osoi_liq
            _fun_hook = function ()
                cnfun!(mask_bgc_vegp, mask_bgc_soilc, bounds_patch, bounds_col,
                       _fun_params, _pftcon_fun, patch, _ws_fun, temperature,
                       soilstate, cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
                       soilbgc_nf, soilbgc_cf, canopystate, soilbgc_ns;
                       dt=dt, nlevdecomp=nlevdecomp, dzsoi_decomp_vals=dzsoi_decomp,
                       use_flexiblecn=config.use_flexiblecn,
                       use_matrixcn=config.use_matrixcn, npcropmin=npcropmin)
                # p2c (unity): patch FUN uptake → column, per layer.
                _fun_p2c!(soilbgc_nf.sminn_to_plant_fun_no3_vr_col,
                          cnveg_nf.sminn_to_plant_fun_no3_vr_patch,
                          patch, mask_bgc_soilc, bounds_col, bounds_patch, nlevdecomp)
                _fun_p2c!(soilbgc_nf.sminn_to_plant_fun_nh4_vr_col,
                          cnveg_nf.sminn_to_plant_fun_nh4_vr_patch,
                          patch, mask_bgc_soilc, bounds_col, bounds_patch, nlevdecomp)
            end
        end

        # 3. Mineral N competition (plant vs decomposers)
        if competition_state !== nothing && competition_params !== nothing
            soil_bgc_competition!(soilbgc_state, soilbgc_nf, soilbgc_cf, soilbgc_ns,
                competition_state, competition_params;
                mask_bgc_soilc=mask_bgc_soilc,
                bounds=bounds_col,
                nlevdecomp=nlevdecomp,
                ndecomp_cascade_transitions=ndecomp_cascade_transitions,
                dzsoi_decomp=dzsoi_decomp,
                pmnf_decomp_cascade=pmnf_decomp_cascade,
                p_decomp_cn_gain=p_decomp_cn_gain,
                cascade_receiver_pool=cascade_receiver_pool,
                use_nitrif_denitrif=config.use_nitrif_denitrif,
                use_fun=config.use_fun,
                fun_hook=_fun_hook)
        end

        # 3b. Plant C/N allocation (calc_plant_nutrient_competition!) — WIRED.
        # Consumes the fpg_col just produced by soil_bgc_competition! and the
        # availc/allometry from the veg-flux chain above; produces cnveg_cf.cpool_to_*c
        # (consumed by c_state_update1! → makes leafc finite) + cnveg_nf alloc fluxes.
        if num_bgc_vegp > 0 && _pftnc !== nothing && cn_shared_params !== nothing &&
           patch !== nothing && crop !== nothing && cnveg_state !== nothing
            # FlexibleCN allocation: plant_calloc = availc (no GPP downregulation),
            # npool_to_* drawn from the existing npool by fractional demand, plus
            # carbon_resp_opt high-C:N tissue turnover. Default path byte-identical.
            if config.use_flexiblecn && canopystate !== nothing
                calc_plant_nutrient_competition_flexiblecn!(mask_bgc_vegp, bounds_patch,
                    _pftnc, cn_shared_params, patch, crop, canopystate, cnveg_state,
                    cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
                    fpg_col = soilbgc_state.fpg_col, dt = dt,
                    use_c13 = config.use_c13, use_c14 = config.use_c14,
                    carbon_resp_opt = config.carbon_resp_opt,
                    use_crop_agsys = config.use_crop_agsys,
                    npcropmin = npcropmin, nrepr = nrepr)
            else
                calc_plant_nutrient_competition!(mask_bgc_vegp, bounds_patch,
                    _pftnc, cn_shared_params, patch, crop, cnveg_state,
                    cnveg_cs, cnveg_cf, cnveg_nf;
                    fpg_col = soilbgc_state.fpg_col,
                    cnveg_ns = cnveg_ns, dt = dt,
                    use_c13 = config.use_c13, use_c14 = config.use_c14,
                    npcropmin = npcropmin, nrepr = nrepr)
            end
        end

        # Soil matrix-CN: seed the turnover diagonal Ksoil = decomp_k_col·dt (before
        # decomp!, which then applies the fpi_vr immobilization correction in place).
        if config.use_soil_matrixcn
            _dk = soilbgc_cf.decomp_k_col
            for i in 1:ndecomp_pools, j in 1:nlevdecomp, c in bounds_col
                if mask_bgc_soilc[c]
                    Ksoil_mat[c, j + (i - 1) * nlevdecomp] = _dk[c, j, i] * dt
                end
            end
        end

        # 4. Actual decomposition
        if decomp_params !== nothing
            soil_biogeochem_decomp!(soilbgc_cf, soilbgc_cs, soilbgc_nf, soilbgc_ns,
                soilbgc_state, cascade_con, decomp_params;
                mask_bgc_soilc=mask_bgc_soilc,
                bounds=bounds_col,
                nlevdecomp=nlevdecomp,
                ndecomp_pools=ndecomp_pools,
                ndecomp_cascade_transitions=ndecomp_cascade_transitions,
                cn_decomp_pools=cn_decomp_pools,
                p_decomp_cpool_loss=p_decomp_cpool_loss,
                pmnf_decomp_cascade=pmnf_decomp_cascade,
                p_decomp_npool_to_din=p_decomp_npool_to_din,
                dzsoi_decomp=dzsoi_decomp,
                use_nitrif_denitrif=config.use_nitrif_denitrif,
                use_soil_matrixcn=config.use_soil_matrixcn, Ksoil_DM=Ksoil_mat)
        end
    end

    # --------------------------------------------------
    # Plant/heterotroph competition for mineral N
    # --------------------------------------------------
    if num_bgc_vegp > 0
        # Phenology phase 1 (if use_fun) — already ported: cn_phenology!(...; phase=1)
        # CNFUNInit — already ported: cn_fun_init!(...)
        # Allocation — already ported:
        #   calc_gpp_mr_availc!(...), calc_crop_allocation_fractions!(...), calc_allometry!(...)
    end
    # NOTE (updated — the old "not yet wired" note here was stale): plant nutrient
    # demand, soil-BGC competition and plant CN allocation ARE ported
    # (biogeochem/nutrient_competition.jl, decomp_competition.jl) and ARE wired in
    # this driver — calc_plant_nutrient_demand! and calc_plant_nutrient_competition!
    # (→ calc_plant_cn_alloc!) both run earlier in this function, whenever the
    # caller supplies the pftcon/crop/shared-params/competition-state bundle.
    # Decomposition runs N-limited on that path.

    # --------------------------------------------------
    # Actual decomposition
    # --------------------------------------------------
    # SoilBiogeochemDecomp — already ported: soil_biogeochem_decomp!(...)

    # --------------------------------------------------
    # Phenology (cn_phenology!) phase 1 then phase 2 — WIRED. For use_fun=false
    # both phases run here (after decomposition, before growth respiration), per
    # CNDriverNoLeaching. Phase 1 updates onset/offset/dormancy state; phase 2
    # emits growth-from-storage, litterfall, live-wood turnover and distributes
    # patch litter to the column via leaf/froot vertical profiles.
    # --------------------------------------------------
    if num_bgc_vegp > 0 && pftcon_main !== nothing && cn_shared_params !== nothing &&
       patch !== nothing && crop !== nothing && cnveg_state !== nothing &&
       temperature !== nothing && water_diag !== nothing && canopystate !== nothing &&
       soilstate !== nothing && gridcell !== nothing
        # Use the module-level params instance, which readParams_CNPhenology! fills
        # from the parameter file. Constructing a fresh PhenologyParams() here meant
        # the run silently used Fortran's *unit-test* defaults (e.g. ndays_off=30)
        # instead of the file values (ndays_off=15) — a 2x-too-long leaf-offset ramp.
        _phparams = cn_phenology_params
        _phstate  = PhenologyState()
        # Resolve the phenology-trigger soil layer from phenology_soil_depth (0.08 m),
        # mirroring Fortran find_soil_layer_containing_depth: the first layer whose
        # lower interface >= depth (CLM5 layers -> layer 3, node ~0.09 m). The
        # cn_phenology_init! default stub returns 1; using the real layer matches
        # Fortran's onset timing (layer 1 is too shallow/responsive -> early onset).
        # zisoi[] is 1-based with zisoi[1]=0 (top) and zisoi[i+1] = lower interface of
        # soil layer i; empty before the vertical grid is built -> fall back to 1.
        _find_phen_layer = function (depth)
            zi = zisoi[]
            isempty(zi) && return 1
            for i in 1:varpar.nlevgrnd
                (i + 1) <= length(zi) && depth <= zi[i + 1] && return i
            end
            return varpar.nlevgrnd
        end
        cn_phenology_init!(_phstate, _phparams, dt; find_soil_layer_fn = _find_phen_layer)
        _pftph = PftConPhenology{Float64}(
            evergreen = Float64.(pftcon_main.evergreen), season_decid = Float64.(pftcon_main.season_decid),
            season_decid_temperate = Float64.(pftcon_main.season_decid_temperate),
            stress_decid = Float64.(pftcon_main.stress_decid), woody = Float64.(pftcon_main.woody),
            leaf_long = Float64.(pftcon_main.leaf_long), leafcn = Float64.(pftcon_main.leafcn),
            frootcn = Float64.(pftcon_main.frootcn), lflitcn = Float64.(pftcon_main.lflitcn),
            livewdcn = Float64.(pftcon_main.livewdcn), deadwdcn = Float64.(pftcon_main.deadwdcn),
            ndays_on = Float64.(pftcon_main.ndays_on), crit_onset_gdd_sf = Float64.(pftcon_main.crit_onset_gdd_sf),
            lf_f = Float64.(pftcon_main.lf_f), fr_f = Float64.(pftcon_main.fr_f),
            biofuel_harvfrac = Float64.(pftcon_main.biofuel_harvfrac),
            repr_structure_harvfrac = Float64.(pftcon_main.repr_structure_harvfrac),
            minplanttemp = Float64.(pftcon_main.minplanttemp), planttemp = Float64.(pftcon_main.planttemp),
            gddmin = Float64.(pftcon_main.gddmin), lfemerg = Float64.(pftcon_main.lfemerg),
            grnfill = Float64.(pftcon_main.grnfill), hybgdd = Float64.(pftcon_main.hybgdd),
            mxmat = Int.(pftcon_main.mxmat), manunitro = Float64.(pftcon_main.manunitro),
            is_pft_known_to_model = Bool.(pftcon_main.is_pft_known_to_model))
        _lprof = soilbgc_state.leaf_prof_patch
        _fprof = soilbgc_state.froot_prof_patch
        for _phase in (1, 2)
            cn_phenology!(_phstate, _phparams, _pftph,
                mask_bgc_vegp, mask_pcropp, mask_bgc_soilc,
                temperature, water_diag, canopystate, soilstate,
                cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
                crop, patch, gridcell, cn_shared_params,
                _lprof, _fprof, _phase;
                varctl = varctl, is_first_step = is_first_step, avg_dayspyr = 365.0,
                prec10_patch = prec10_patch,
                # REQUIRED: without these, cn_litter_to_column! defaults to
                # nlevdecomp=1 and destroys ~86% of all phenological leaf/fine-root
                # litterfall C and N (see phenology.jl).
                nlevdecomp = nlevdecomp, i_litr_min = i_litr_min,
                i_litr_max = i_litr_max, npcropmin = npcropmin)
        end
    end

    # --------------------------------------------------
    # Growth respiration (cn_gresp!) — WIRED. Reads cpool_to_*c (allocation) +
    # *_xfer_to_*c (phenology); writes the *_gr growth-respiration fluxes consumed
    # by c_state_update1!. PftConGrowthResp is a subset of the main pftcon.
    # --------------------------------------------------
    if num_bgc_vegp > 0 && pftcon_main !== nothing && patch !== nothing
        _pftgr = PftConGrowthResp{Float64}(
            woody = Float64.(pftcon_main.woody), grperc = Float64.(pftcon_main.grperc),
            grpnow = Float64.(pftcon_main.grpnow))
        cn_gresp!(mask_bgc_vegp, bounds_patch, _pftgr, patch, cnveg_cf;
                  npcropmin = npcropmin, nrepr = nrepr)
    end

    # --------------------------------------------------
    # CStateUpdate0
    # --------------------------------------------------
    c_state_update0!(cnveg_cs, cnveg_cf;
        mask_soilp=mask_bgc_vegp,
        bounds_patch=bounds_patch,
        dt=dt)

    # --------------------------------------------------
    # CIsoFlux1 — non-mortality carbon-isotope fluxes (WIRED).
    # Scales every veg + decomposition C flux into its C13/C14 counterpart by the
    # donor-pool ratio, then gathers patch-level litterfall to columns.
    # --------------------------------------------------
    if _have_iso && num_bgc_vegp > 0
        _lprof = soilbgc_state.leaf_prof_patch
        _fprof = soilbgc_state.froot_prof_patch
        for (_tag, _ics, _icf, _iscs, _iscf) in _iso_active
            c_iso_flux1!(soilbgc_state, soilbgc_cf, soilbgc_cs, cnveg_cf, cnveg_cs,
                _iscf, _iscs, _icf, _ics,
                mask_bgc_soilc, mask_bgc_vegp, bounds_col, bounds_patch,
                collect(Int, cascade_donor_pool), nlevdecomp,
                ndecomp_cascade_transitions, _tag;
                use_crop=config.use_crop, nrepr=nrepr, npcropmin=npcropmin,
                patch_column=_iso_pcol, patch_itype=_iso_pivt, patch_wtcol=_iso_pwtcol,
                lf_f=_iso_lf_f, fr_f=_iso_fr_f, leaf_prof=_lprof, froot_prof=_fprof)
        end
    end

    # --------------------------------------------------
    # CStateUpdate1
    # --------------------------------------------------
    c_state_update1!(cnveg_cs, cnveg_cf, soilbgc_cf;
        mask_soilc=mask_bgc_soilc,
        mask_soilp=mask_bgc_vegp,
        bounds_col=bounds_col,
        bounds_patch=bounds_patch,
        patch_column=patch_column,
        ivt=ivt,
        woody=woody,
        cascade_donor_pool=cascade_donor_pool,
        cascade_receiver_pool=cascade_receiver_pool,
        harvdate=harvdate,
        col_is_fates=col_is_fates,
        nlevdecomp=nlevdecomp,
        ndecomp_cascade_transitions=ndecomp_cascade_transitions,
        i_litr_min=i_litr_min,
        i_litr_max=i_litr_max,
        i_cwd=i_cwd,
        npcropmin=npcropmin,
        nrepr=nrepr,
        use_matrixcn=config.use_matrixcn,
        use_soil_matrixcn=config.use_soil_matrixcn,
        dribble_crophrv_xsmrpool_2atm=config.dribble_crophrv_xsmrpool_2atm,
        dt=dt)

    # --------------------------------------------------
    # NStateUpdate1
    # --------------------------------------------------
    n_state_update1!(cnveg_ns, cnveg_nf, soilbgc_nf;
        mask_soilc=mask_bgc_soilc,
        mask_soilp=mask_bgc_vegp,
        bounds_col=bounds_col,
        bounds_patch=bounds_patch,
        ivt=ivt,
        woody=woody,
        col_is_fates=col_is_fates,
        nlevdecomp=nlevdecomp,
        i_litr_min=i_litr_min,
        i_litr_max=i_litr_max,
        i_cwd=i_cwd,
        npcropmin=npcropmin,
        nrepr=nrepr,
        use_matrixcn=config.use_matrixcn,
        use_soil_matrixcn=config.use_soil_matrixcn,
        use_fun=config.use_fun,
        dt=dt)

    # SoilBiogeochemNStateUpdate1 — WIRED. Applies the mineral-N fluxes
    # (dep/fix, gross mineralization, immobilization, plant uptake, nitrification,
    # denitrification, supplement) to smin_nh4_vr/smin_no3_vr (the nitrif path).
    # Fortran CNDriverNoLeaching line 658, after NStateUpdate1.
    if config.use_nitrif_denitrif && _has_decomp
        soilbiogeochem_n_state_update1!(soilbgc_ns, soilbgc_nf, soilbgc_state;
            mask_bgc_soilc=mask_bgc_soilc, bounds_col=bounds_col,
            nlevdecomp=nlevdecomp, dt=dt, use_fun=config.use_fun)
    end

    # CNPrecisionControl — WIRED
    # Called with integer filter built from mask (cn_precision_control! uses filter array)
    if num_bgc_vegp > 0
        filter_bgc_vegp = findall(mask_bgc_vegp)
        cn_precision_control!(cnveg_cs, cnveg_ns, filter_bgc_vegp, ivt;
            use_c13=config.use_c13, use_c14=config.use_c14,
            use_crop=config.use_crop, use_nguardrail=config.use_nguardrail,
            use_matrixcn=config.use_matrixcn)
    end
    # SoilBiogeochemLittVertTransp — litter vertical transport
    # Skip when nlevdecomp <= 1: vertical transport is a no-op with a single level
    if _has_decomp && nlevdecomp > 1 && litter_params !== nothing && col !== nothing &&
       grc !== nothing && active_layer !== nothing && zisoi_vals !== nothing
        litter_vert_transp!(soilbgc_cs, soilbgc_cf, soilbgc_ns, soilbgc_nf,
            soilbgc_state, col, grc, cascade_con, active_layer, litter_params;
            mask_bgc_soilc=mask_bgc_soilc,
            bounds=bounds_col,
            dtime=dt,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            zsoi_vals=zsoi_vals,
            dzsoi_decomp_vals=dzsoi_decomp,
            zisoi_vals=zisoi_vals,
            use_soil_matrixcn=config.use_soil_matrixcn,
            use_c13=(config.use_c13 && c13_soilbgc_cs !== nothing && c13_soilbgc_cf !== nothing),
            use_c14=(config.use_c14 && c14_soilbgc_cs !== nothing && c14_soilbgc_cf !== nothing),
            c13_cs=c13_soilbgc_cs, c13_cf=c13_soilbgc_cf,
            c14_cs=c14_soilbgc_cs, c14_cf=c14_soilbgc_cf)
    end

    # --------------------------------------------------
    # Gap mortality and Update2
    # --------------------------------------------------
    if num_bgc_vegp > 0
        # CNGapMortality — WIRED. Background gap mortality (r_mort=0.02/yr CLM5
        # default; DGVS path off so dgvs is unused/zero) → m_*_to_litter patch
        # fluxes, then gathered patch→column via the vertical profiles. Feeds
        # c_state_update2!/n_state_update2! below.
        if pftcon_main !== nothing && patch !== nothing && canopystate !== nothing
            _np = length(bounds_patch)
            _gmp = GapMortalityParams{Float64}(k_mort = 0.3,
                r_mort = fill(0.02, length(pftcon_main.woody)))
            _pftgm = PftConGapMort{Float64}(
                woody = Float64.(pftcon_main.woody), leafcn = Float64.(pftcon_main.leafcn),
                livewdcn = Float64.(pftcon_main.livewdcn),
                lf_f = Float64.(pftcon_main.lf_f), fr_f = Float64.(pftcon_main.fr_f))
            _dgvs = DgvsGapMortData{Float64}(greffic_patch = zeros(_np),
                heatstress_patch = zeros(_np), nind_patch = zeros(_np))
            cn_gap_mortality!(mask_bgc_vegp, bounds_patch, _gmp, _pftgm, _dgvs, patch,
                canopystate, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
                dt = dt, days_per_year = 365.0, use_cndv = false,
                use_matrixcn = config.use_matrixcn, spinup_state = 0, npcropmin = npcropmin)
            cn_gap_patch_to_column!(mask_bgc_vegp, bounds_patch, _pftgm, patch,
                cnveg_cf, cnveg_nf,
                soilbgc_state.leaf_prof_patch, soilbgc_state.froot_prof_patch,
                soilbgc_state.croot_prof_patch, soilbgc_state.stem_prof_patch;
                nlevdecomp = nlevdecomp, i_litr_min = i_litr_min,
                i_litr_max = i_litr_max, i_met_lit = i_litr_min)
        end
        # CIsoFlux2 — gap-mortality carbon-isotope fluxes (WIRED).
        if _have_iso
            _mvegp = BitVector(mask_bgc_vegp)
            for (_tag, _ics, _icf, _iscs, _iscf) in _iso_active
                c_iso_flux2!(soilbgc_state, cnveg_cf, cnveg_cs, _icf, _ics,
                    _mvegp, bounds_patch, nlevdecomp, _tag;
                    patch_column=_iso_pcol, patch_itype=_iso_pivt,
                    patch_wtcol=_iso_pwtcol, lf_f=_iso_lf_f, fr_f=_iso_fr_f)
            end
        end

        c_state_update2!(cnveg_cs, cnveg_cf, soilbgc_cs;
            mask_soilc=mask_bgc_soilc,
            mask_soilp=mask_bgc_vegp,
            bounds_col=bounds_col,
            bounds_patch=bounds_patch,
            nlevdecomp=nlevdecomp,
            i_litr_min=i_litr_min,
            i_litr_max=i_litr_max,
            i_cwd=i_cwd,
            use_matrixcn=config.use_matrixcn,
            use_soil_matrixcn=config.use_soil_matrixcn,
            dt=dt)

        n_state_update2!(cnveg_ns, cnveg_nf, soilbgc_ns;
            mask_soilc=mask_bgc_soilc,
            mask_soilp=mask_bgc_vegp,
            bounds_col=bounds_col,
            bounds_patch=bounds_patch,
            nlevdecomp=nlevdecomp,
            i_litr_min=i_litr_min,
            i_litr_max=i_litr_max,
            i_cwd=i_cwd,
            use_matrixcn=config.use_matrixcn,
            use_soil_matrixcn=config.use_soil_matrixcn,
            dt=dt)

        # Harvest (Update2h) — CNHarvest WIRED. Compute the patch-level wood-harvest
        # mortality fluxes then gather patch→column into harvest_c/n_to_litr/cwd, which
        # c_state_update2h! (or, in matrix mode, the soil-matrix B-input) consumes.
        if config.do_harvest && dynharvest_state !== nothing &&
           patch !== nothing && pftcon_main !== nothing
            cn_harvest!(dynharvest_state, mask_bgc_vegp, patch, pftcon_main,
                soilbgc_state, cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf;
                dt=dt, nlevdecomp=nlevdecomp, i_litr_min=i_litr_min,
                i_litr_max=i_litr_max, i_met_lit=i_litr_min)
            cn_harvest_pft_to_column!(mask_bgc_vegp, patch, pftcon_main, soilbgc_state,
                cnveg_cf, cnveg_nf; nlevdecomp=nlevdecomp, i_litr_min=i_litr_min,
                i_litr_max=i_litr_max, i_met_lit=i_litr_min)
        end
        # CIsoFlux2h — harvest-mortality carbon-isotope fluxes (WIRED).
        if _have_iso
            _mvegp = BitVector(mask_bgc_vegp)
            for (_tag, _ics, _icf, _iscs, _iscf) in _iso_active
                c_iso_flux2h!(soilbgc_state, cnveg_cf, cnveg_cs, _icf, _ics,
                    _mvegp, bounds_patch, nlevdecomp, _tag;
                    patch_column=_iso_pcol, patch_itype=_iso_pivt,
                    patch_wtcol=_iso_pwtcol, lf_f=_iso_lf_f, fr_f=_iso_fr_f)
            end
        end

        c_state_update2h!(cnveg_cs, cnveg_cf, soilbgc_cs;
            mask_soilc=mask_bgc_soilc,
            mask_soilp=mask_bgc_vegp,
            bounds_col=bounds_col,
            bounds_patch=bounds_patch,
            nlevdecomp=nlevdecomp,
            i_litr_min=i_litr_min,
            i_litr_max=i_litr_max,
            i_cwd=i_cwd,
            use_matrixcn=config.use_matrixcn,
            use_soil_matrixcn=config.use_soil_matrixcn,
            dt=dt)

        n_state_update2h!(cnveg_ns, cnveg_nf, soilbgc_ns;
            mask_soilc=mask_bgc_soilc,
            mask_soilp=mask_bgc_vegp,
            bounds_col=bounds_col,
            bounds_patch=bounds_patch,
            nlevdecomp=nlevdecomp,
            i_litr_min=i_litr_min,
            i_litr_max=i_litr_max,
            i_cwd=i_cwd,
            use_matrixcn=config.use_matrixcn,
            use_soil_matrixcn=config.use_soil_matrixcn,
            dt=dt)

        # Gross unrepresented landcover change (Update2g) — CNGrossUnrep WIRED.
        if config.do_grossunrep && dyngrossunrep_state !== nothing &&
           patch !== nothing && pftcon_grossunrep !== nothing
            cn_gross_unrep!(mask_bgc_vegp, dyngrossunrep_state, pftcon_grossunrep, patch,
                cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf, soilbgc_state;
                dt=dt, nlevdecomp=nlevdecomp, i_litr_min=i_litr_min,
                i_litr_max=i_litr_max, i_met_lit=i_litr_min)
            cn_gross_unrep_pft_to_column!(mask_bgc_vegp, pftcon_grossunrep, patch,
                cnveg_cf, cnveg_nf, soilbgc_state; nlevdecomp=nlevdecomp,
                i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_met_lit=i_litr_min)
        end
        # CIsoFlux2g — gross-unrepresented-LCC carbon-isotope fluxes (WIRED).
        if _have_iso
            _mvegp = BitVector(mask_bgc_vegp)
            for (_tag, _ics, _icf, _iscs, _iscf) in _iso_active
                c_iso_flux2g!(soilbgc_state, cnveg_cf, cnveg_cs, _icf, _ics,
                    _mvegp, bounds_patch, nlevdecomp, _tag;
                    patch_column=_iso_pcol, patch_itype=_iso_pivt,
                    patch_wtcol=_iso_pwtcol, lf_f=_iso_lf_f, fr_f=_iso_fr_f)
            end
        end

        c_state_update2g!(cnveg_cs, cnveg_cf, soilbgc_cs;
            mask_soilc=mask_bgc_soilc,
            mask_soilp=mask_bgc_vegp,
            bounds_col=bounds_col,
            bounds_patch=bounds_patch,
            nlevdecomp=nlevdecomp,
            i_litr_min=i_litr_min,
            i_litr_max=i_litr_max,
            i_cwd=i_cwd,
            use_matrixcn=config.use_matrixcn,
            use_soil_matrixcn=config.use_soil_matrixcn,
            dt=dt)

        n_state_update2g!(cnveg_ns, cnveg_nf, soilbgc_ns;
            mask_soilc=mask_bgc_soilc,
            mask_soilp=mask_bgc_vegp,
            bounds_col=bounds_col,
            bounds_patch=bounds_patch,
            nlevdecomp=nlevdecomp,
            i_litr_min=i_litr_min,
            i_litr_max=i_litr_max,
            i_cwd=i_cwd,
            use_matrixcn=config.use_matrixcn,
            use_soil_matrixcn=config.use_soil_matrixcn,
            dt=dt)
    end  # num_bgc_vegp > 0

    # CNPrecisionControl (post gap-mortality) — WIRED
    if num_bgc_vegp > 0
        filter_bgc_vegp2 = findall(mask_bgc_vegp)
        cn_precision_control!(cnveg_cs, cnveg_ns, filter_bgc_vegp2, ivt;
            use_c13=config.use_c13, use_c14=config.use_c14,
            use_crop=config.use_crop, use_nguardrail=config.use_nguardrail,
            use_matrixcn=config.use_matrixcn)
    end
    # NOTE: the wood-product pools ARE ported (types/cn_products.jl +
    # biogeochem/cn_products_mod.jl, initialized as CLMInstances.cn_products), but
    # cn_products_update! is not CALLED from this driver — it needs the harvest
    # fluxes threaded in. Ported, not wired.

    # --------------------------------------------------
    # Fire and Update3
    # --------------------------------------------------
    if num_bgc_vegp > 0
        # --------------------------------------------------
        # Fire: burned area (cnfire_area!) then C/N fluxes (cnfire_fluxes_dispatch!)
        # WIRED — runs only when config.cnfire_method !== :nofire AND the full fire
        # bundle is supplied. The fire fluxes (m_*_to_fire patch fluxes +
        # m_decomp_*_to_fire_vr / m_c_to_litr_fire / fire_mortality_c_to_cwdc column
        # fluxes) are then consumed by c_state_update3! below — exactly the Fortran
        # CNDriverNoLeaching order (CNFireArea → CNFireFluxes → CStateUpdate3).
        # The column soil temperature the old placeholder said it needed is
        # temperature.t_soi17cm_col (top-17cm soil T), which the driver already
        # receives via the `temperature` arg.
        _fire_active = config.cnfire_method !== :nofire &&
            fire_li2014 !== nothing && pftcon_fire_li2014 !== nothing &&
            fire_data !== nothing && cnfire_const !== nothing &&
            cnfire_params !== nothing && pftcon_fire !== nothing &&
            patch !== nothing && col !== nothing && grc !== nothing &&
            soilstate !== nothing && h2osoi_vol !== nothing &&
            cnveg_state !== nothing && cascade_con !== nothing &&
            temperature !== nothing && dgvs_fire !== nothing
        if _fire_active
            _mexp  = mask_exposedvegp_fire === nothing ? mask_bgc_vegp : mask_exposedvegp_fire
            _mnoex = mask_noexposedvegp_fire === nothing ?
                     [mask_bgc_vegp[p] ? false : false for p in eachindex(mask_bgc_vegp)] :
                     mask_noexposedvegp_fire
            # totlitc_col / totsomc_col fuel inputs: sum the column litter/SOM pools
            # over (layer, pool) per the decomp cascade is_litter / !is_litter split.
            _is_litr = cascade_con.is_litter
            _Tc = eltype(soilbgc_cs.decomp_cpools_vr_col)
            _totlitc = fill!(similar(soilbgc_cs.decomp_cpools_vr_col, _Tc, nc), zero(_Tc))
            _totsomc = fill!(similar(soilbgc_cs.decomp_cpools_vr_col, _Tc, nc), zero(_Tc))
            for c in bounds_col
                mask_bgc_soilc[c] || continue
                lit = zero(_Tc); som = zero(_Tc)
                for k in 1:ndecomp_pools, j in 1:nlevdecomp
                    v = soilbgc_cs.decomp_cpools_vr_col[c, j, k] *
                        (dzsoi_decomp === nothing ? one(_Tc) : _Tc(dzsoi_decomp[j]))
                    if k <= length(_is_litr) && _is_litr[k]
                        lit += v
                    else
                        som += v
                    end
                end
                _totlitc[c] = lit; _totsomc[c] = som
            end

            # 1. Burned area — factory dispatch on config.cnfire_method.
            cnfire_area!(config.cnfire_method,
                fire_li2014, pftcon_fire_li2014, fire_data, cnfire_const,
                cnfire_params, pftcon_fire,
                mask_bgc_soilc, mask_bgc_vegp, _mexp, _mnoex,
                bounds_col, bounds_patch,
                patch, col, grc, soilstate, h2osoi_vol,
                cnveg_state, cnveg_cs, cascade_con,
                _totlitc, soilbgc_cs.decomp_cpools_vr_col, temperature.t_soi17cm_col;
                forc_rh_grc = forc_rh_grc, forc_wind_grc = forc_wind_grc,
                forc_t_col = forc_t_fire_col, forc_rain_col = forc_rain_fire_col,
                forc_snow_col = forc_snow_fire_col,
                prec60_patch = prec60_patch, prec10_patch = prec10_patch,
                prec30_patch = prec30_patch, rh30_patch = rh30_patch,
                fsat_col = fsat_fire_col, wf_col = wf_fire_col, wf2_col = wf2_fire_col,
                dt = dt, dayspyr = 365.0,
                kmo = fire_kmo, kda = fire_kda, mcsec = fire_mcsec, nstep = fire_nstep,
                nlevgrnd = nlevgrnd_fire, nlevdecomp = nlevdecomp,
                ndecomp_pools = ndecomp_pools, transient_landcover = transient_landcover)

            # 2. Fire C/N fluxes — factory dispatch; produces the m_*_to_fire patch
            #    fluxes + column decomp/litter fire fluxes consumed by c_state_update3!.
            #    Returns the active-fire column/patch masks (threaded to the leaching
            #    driver's fire N state update via mask_actfirec/mask_actfirep).
            _maskc, _maskp = cnfire_fluxes_dispatch!(config.cnfire_method,
                mask_bgc_soilc, mask_bgc_vegp, bounds_col, bounds_patch,
                cnfire_const, pftcon_fire, patch, col, grc, dgvs_fire,
                cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
                soilbgc_cf, cascade_con,
                soilbgc_state.leaf_prof_patch, soilbgc_state.froot_prof_patch,
                soilbgc_state.croot_prof_patch, soilbgc_state.stem_prof_patch,
                _totsomc, soilbgc_cs.decomp_cpools_vr_col,
                soilbgc_ns.decomp_npools_vr_col, soilbgc_cf.somc_fire_col;
                dt = dt, dayspyr = 365.0, nlevdecomp = nlevdecomp,
                ndecomp_pools = ndecomp_pools, i_met_lit = i_litr_min,
                i_litr_max = i_litr_max, transient_landcover = transient_landcover,
                use_matrixcn = config.use_matrixcn,
                kmo = fire_kmo, kda = fire_kda, mcsec = fire_mcsec)
            # Copy the returned active-fire masks into the caller-provided outputs.
            for c in eachindex(mask_actfirec); mask_actfirec[c] = _maskc[c]; end
            for p in eachindex(mask_actfirep); mask_actfirep[p] = _maskp[p]; end
        end
        # CIsoFlux3 — fire-mortality carbon-isotope fluxes (WIRED). Scales the fire C
        # fluxes into isotopic counterparts (zero here until fire fluxes are wired,
        # but the cascade stays complete and mass/ratio-consistent).
        if _have_iso
            _lprof = soilbgc_state.leaf_prof_patch
            _fprof = soilbgc_state.froot_prof_patch
            _cprof = soilbgc_state.croot_prof_patch
            _sprof = soilbgc_state.stem_prof_patch
            for (_tag, _ics, _icf, _iscs, _iscf) in _iso_active
                c_iso_flux3!(soilbgc_state, soilbgc_cs, cnveg_cf, cnveg_cs,
                    _icf, _ics, _iscs,
                    mask_bgc_vegp, bounds_patch, nlevdecomp, ndecomp_pools,
                    i_litr_min, i_litr_max, _tag;
                    patch_column=_iso_pcol, patch_itype=_iso_pivt,
                    patch_wtcol=_iso_pwtcol, lf_f=_iso_lf_f, fr_f=_iso_fr_f,
                    leaf_prof=_lprof, froot_prof=_fprof,
                    stem_prof=_sprof, croot_prof=_cprof)
            end
        end

        c_state_update3!(cnveg_cs, cnveg_cf, soilbgc_cs;
            mask_soilc=mask_bgc_soilc,
            mask_soilp=mask_bgc_vegp,
            bounds_col=bounds_col,
            bounds_patch=bounds_patch,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            i_litr_min=i_litr_min,
            i_litr_max=i_litr_max,
            i_cwd=i_cwd,
            use_matrixcn=config.use_matrixcn,
            use_soil_matrixcn=config.use_soil_matrixcn,
            dt=dt)

        # Matrix-CN veg-C solve — WIRED (gated on use_matrixcn, phase 1: natveg).
        # In matrix mode the c_state_update1/2/3 blocks above skipped every veg-pool
        # increment (only cpool, the input source, decremented), so the 18 veg C
        # pools are still at start-of-step and every allocation/phenology/gap/fire
        # transfer has been computed into its flux field. Assemble those into the
        # transfer matrix + B-input and advance all pools in one solve (== the
        # sequential update at matrixcheck=false; validated in test_cn_veg_matrix_wiring).
        if config.use_matrixcn
            # Crop adds the 3 grain pools (nvegcpool 18→21) + shifts retransn to pool 22.
            _uc = config.use_crop
            _mcounts = veg_matrix_transfer_counts(_uc)
            _ncpool = _uc ? NVEGPOOL_NATVEG + NVEGPOOL_CROP : NVEGPOOL_NATVEG   # 21 / 18
            _nnpool = _ncpool + 1                                              # 22 / 19 (incl. retransn)
            # Device backend: when the state arrays live on a GPU, derive the ref prototype
            # + working precision for the matrix solve (nothing/Float64 keeps the host path).
            _mref = cnveg_cs.leafc_patch isa Array ? nothing : cnveg_cs.leafc_patch
            _mFT = eltype(cnveg_cs.leafc_patch)
            if cnveg_cf.ileafst_to_ileafxf_ph == 0
                cn_veg_matrix_c_topology!(cnveg_cf; use_crop=_uc, nvegcpool=_ncpool)
            end
            cn_veg_matrix_alloc_c!(cnveg_cf, mask_bgc_vegp, bounds_patch; use_crop=_uc)
            cn_veg_matrix_accumulate_ph_c!(cnveg_cf, cnveg_cs, mask_bgc_vegp, bounds_patch; dt=dt, use_crop=_uc)
            cn_veg_matrix_accumulate_gm_c!(cnveg_cf, cnveg_cs, mask_bgc_vegp, bounds_patch; dt=dt)
            cn_veg_matrix_accumulate_fi_c!(cnveg_cf, cnveg_cs, mask_bgc_vegp, bounds_patch; dt=dt)
            cn_veg_matrix_solve_c!(cnveg_cs, cnveg_cf;
                mask_soilp=mask_bgc_vegp, bounds_patch=bounds_patch,
                ivt=ivt, woody=woody, npcropmin=npcropmin, nvegcpool=_ncpool,
                counts=_mcounts, dt=dt, num_actfirep=count(mask_actfirep),
                ref=_mref, FT=_mFT, state=veg_c_solve_state)

            # Matrix-CN veg-N solve (phase 2) — same as C, plus the retranslocation
            # pool (nvegnpool = nvegcpool+1 incl. retransn). The N fire fluxes are computed
            # by the fire block above; n_state_update1/2/3 all skip veg-N pool increments
            # in matrix mode (leaving pools at start-of-step), so this one solve subsumes
            # them (== sequential; validated in test_cn_veg_matrix_wiring_n[/_crop_n]).
            if cnveg_nf.ileaf_to_iretransn_ph == 0
                cn_veg_matrix_n_topology!(cnveg_nf; use_crop=_uc, nvegnpool=_nnpool)
            end
            cn_veg_matrix_alloc_n!(cnveg_nf, mask_bgc_vegp, bounds_patch; use_crop=_uc)
            cn_veg_matrix_accumulate_phn!(cnveg_nf, cnveg_ns, mask_bgc_vegp, bounds_patch; dt=dt, use_crop=_uc, nvegnpool=_nnpool)
            cn_veg_matrix_accumulate_gmn!(cnveg_nf, cnveg_ns, mask_bgc_vegp, bounds_patch; dt=dt, nvegnpool=_nnpool)
            cn_veg_matrix_accumulate_fin!(cnveg_nf, cnveg_ns, mask_bgc_vegp, bounds_patch; dt=dt, nvegnpool=_nnpool)
            cn_veg_matrix_solve_n!(cnveg_ns, cnveg_nf;
                mask_soilp=mask_bgc_vegp, bounds_patch=bounds_patch,
                ivt=ivt, npcropmin=npcropmin, nvegnpool=_nnpool,
                counts=_mcounts, dt=dt, num_actfirep=count(mask_actfirep),
                ref=_mref, FT=_mFT, state=veg_n_solve_state)
        end

        # C14Decay — ported (c14_decay!, biogeochem/carbon_isotopes.jl) and called
        # from clm_drv! (clm_driver.jl:2161) under use_c14, not from here.
    end  # num_bgc_vegp > 0

    # Soil matrix-CN solve (phase 3) — WIRED (gated on use_soil_matrixcn). Runs after
    # all veg→soil litter fluxes are computed (phenology/gap/fire → litr/cwd). In
    # matrix mode the sequential decomp cascade + vertical transport + litter-input
    # state updates all skipped (leaving decomp_cpools_vr at start-of-step); this one
    # solve advances the pools by (A·K + V)·X + B, with K=Ksoil (decomp_k·dt, fpi-
    # corrected in decomp!), V=tri_ma_vr (built by litter_vert_transp!), B=the litter
    # inputs. Validated: cn_soil_matrix == sequential cascade (test_cn_soil_matrix) +
    # tri_ma_vr tracks sequential transport (test_lvt_tri_ma_vr).
    if _has_decomp && config.use_soil_matrixcn && decomp_params !== nothing
        cn_soil_matrix_advance!(cascade_con, soilbgc_cs, soilbgc_ns, soilbgc_cf,
            cnveg_cf, cnveg_nf; Ksoil=Ksoil_mat,
            mask_soilc=mask_bgc_soilc, bounds_col=bounds_col,
            nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
            ndecomp_cascade_transitions=ndecomp_cascade_transitions,
            i_litr_min=i_litr_min, i_litr_max=i_litr_max, i_cwd=i_cwd, dt=dt,
            num_actfirec=count(mask_actfirec),
            # Harvest + gross-unrep litter/CWD inputs enter the B-input when transient
            # landcover is active; their state updates (c/n_state_update2h/2g) skip the
            # direct decomp-pool addition under use_soil_matrixcn. dwt enters via the
            # persistent soilbgc_nf/cf.matrix_Cinput/Ninput_col (accumulated in dyn_subgrid).
            transient_landcover=transient_landcover, soilbgc_nf=soilbgc_nf,
            soil_matrix_state=soil_matrix_state,
            ref=(soilbgc_cs.decomp_cpools_vr_col isa Array ? nothing : soilbgc_cs.decomp_cpools_vr_col),
            FT=eltype(soilbgc_cs.decomp_cpools_vr_col))
    end

    # CNPrecisionControl (post fire) — WIRED
    if num_bgc_vegp > 0
        filter_bgc_vegp3 = findall(mask_bgc_vegp)
        cn_precision_control!(cnveg_cs, cnveg_ns, filter_bgc_vegp3, ivt;
            use_c13=config.use_c13, use_c14=config.use_c14,
            use_crop=config.use_crop, use_nguardrail=config.use_nguardrail,
            use_matrixcn=config.use_matrixcn)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# cn_driver_leaching! — N leaching and state updates
# Ported from CNDriverLeaching in CNDriverMod.F90
# ---------------------------------------------------------------------------

"""
    cn_driver_leaching!(config; kwargs...)

Update the nitrogen leaching rate as a function of soluble mineral N and
total soil water outflow. Also update nitrogen state variables.

Ported from `CNDriverLeaching` in `CNDriverMod.F90`.
"""
function cn_driver_leaching!(
        config::CNDriverConfig;
        # Masks
        mask_bgc_soilc::AbstractVector{Bool},
        mask_bgc_vegp::AbstractVector{Bool},
        mask_actfirec::AbstractVector{Bool} = falses(0),
        mask_actfirep::AbstractVector{Bool} = falses(0),
        # Bounds
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        # Decomposition dimensions
        nlevdecomp::Int,
        ndecomp_pools::Int = 7,
        ndecomp_cascade_transitions::Int = 7,
        i_litr_min::Int = 1,
        i_litr_max::Int = 3,
        i_cwd::Int = 4,
        # Time step
        dt::Real,
        # State/flux data
        soilbgc_ns::SoilBiogeochemNitrogenStateData,
        soilbgc_nf::SoilBiogeochemNitrogenFluxData,
        cnveg_ns::CNVegNitrogenStateData,
        cnveg_nf::CNVegNitrogenFluxData,
        # ---- N leaching inputs (SoilBiogeochemNLeachingMod) ----------------
        # Supplied by the facade from the water state/flux + column geometry.
        # All default to `nothing`/empty so a caller that does not provide them
        # gets exactly the previous behaviour (no leaching loss).
        nleach_params::Union{NLeachingParams, Nothing} = nothing,
        h2osoi_liq::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        qflx_drain::Union{AbstractVector{<:Real}, Nothing} = nothing,
        qflx_surf::Union{AbstractVector{<:Real}, Nothing} = nothing,
        col_dz::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        zisoi::Union{AbstractVector{<:Real}, Nothing} = nothing,
        nlevsoi::Int = 0)

    num_bgc_vegp = count(mask_bgc_vegp)

    # SoilBiogeochemNLeaching — WIRED. Fortran CNDriverLeaching line 1070: the N
    # LOSS term (mineral N leached with drainage + lost with surface runoff),
    # called immediately before NStateUpdateLeaching, which subtracts it from the
    # pool. Previously a dead port: the loss was never computed, so mineral N had
    # an input-free, loss-free budget.
    if h2osoi_liq !== nothing && qflx_drain !== nothing && qflx_surf !== nothing &&
       col_dz !== nothing && zisoi !== nothing && nlevsoi > 0
        n_leaching!(soilbgc_nf, soilbgc_ns,
            nleach_params === nothing ? NLeachingParams() : nleach_params;
            mask_bgc_soilc      = mask_bgc_soilc,
            bounds              = bounds_col,
            nlevdecomp          = nlevdecomp,
            nlevsoi             = nlevsoi,
            dt                  = dt,
            h2osoi_liq          = h2osoi_liq,
            qflx_drain          = qflx_drain,
            qflx_surf           = qflx_surf,
            col_dz              = col_dz,
            zisoi               = zisoi,
            use_nitrif_denitrif = config.use_nitrif_denitrif)
    end

    # NStateUpdateLeaching
    n_state_update_leaching!(soilbgc_ns, soilbgc_nf;
        mask_soilc=mask_bgc_soilc,
        bounds_col=bounds_col,
        nlevdecomp=nlevdecomp,
        use_nitrif_denitrif=config.use_nitrif_denitrif,
        dt=dt)

    # NStateUpdate3 (fire N)
    if num_bgc_vegp > 0
        n_state_update3!(cnveg_ns, cnveg_nf, soilbgc_ns;
            mask_soilc=mask_bgc_soilc,
            mask_soilp=mask_bgc_vegp,
            bounds_col=bounds_col,
            bounds_patch=bounds_patch,
            nlevdecomp=nlevdecomp,
            ndecomp_pools=ndecomp_pools,
            i_litr_min=i_litr_min,
            i_litr_max=i_litr_max,
            i_cwd=i_cwd,
            use_matrixcn=config.use_matrixcn,
            use_soil_matrixcn=config.use_soil_matrixcn,
            dt=dt)
    end

    # Matrix solutions (CNVegMatrix / CNSoilMatrix) — PORTED AND WIRED, higher up in
    # this function: the veg-C/veg-N solves (cn_veg_matrix_solve_c!/_n!) run under
    # `config.use_matrixcn` and the soil-C/N solve (cn_soil_matrix_advance!) under
    # `config.use_soil_matrixcn`. Both default off → byte-identical to the
    # sequential path. (The old "not yet ported" note here was stale.)

    return nothing
end

# ---------------------------------------------------------------------------
# cn_driver_summarize_states! — Summarize state variables
# Ported from CNDriverSummarizeStates in CNDriverMod.F90
# ---------------------------------------------------------------------------

"""
    cn_driver_summarize_states!(config; kwargs...)

Call to all CN and SoilBiogeochem summary routines for state variables.

Ported from `CNDriverSummarizeStates` in `CNDriverMod.F90`.

Note: the CNVeg carbon/nitrogen state Summary methods ARE ported and ARE called
here (`cnveg_carbon_state_summary!` / `cnveg_nitrogen_state_summary!`, in
`types/cn_veg_carbon_state.jl` / `cn_veg_nitrogen_state.jl`). The SoilBiogeochem
state summaries are ported too (`soil_bgc_carbon_state_summary!` /
`soil_bgc_nitrogen_state_summary!`) but are NOT called from here yet — see the
inline comment below for what they still need threaded in.
"""
function cn_driver_summarize_states!(
        config::CNDriverConfig;
        mask_bgc_soilc::AbstractVector{Bool},
        mask_bgc_vegp::AbstractVector{Bool},
        mask_allc::AbstractVector{Bool} = mask_bgc_soilc,
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        cnveg_cs::CNVegCarbonStateData,
        cnveg_ns::CNVegNitrogenStateData,
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData,
        patch_itype::Union{AbstractVector{<:Integer},Nothing}=nothing,
        # Soil-BGC state-summary inputs (decomp infrastructure). Supplied => the
        # soil summaries run here too, as in Fortran. Omitted => skipped (previous
        # behaviour), so callers without the cascade are unaffected.
        cascade_con::Union{DecompCascadeConData, Nothing} = nothing,
        col::Union{ColumnData, Nothing} = nothing,
        patch::Union{PatchData, Nothing} = nothing,
        nlevdecomp::Int = 0,
        ndecomp_pools::Int = 0,
        dzsoi_decomp::Union{AbstractVector{<:Real}, Nothing} = nothing,
        zisoi_vals::Union{AbstractVector{<:Real}, Nothing} = nothing)

    if count(mask_bgc_vegp) > 0
        # cnveg_carbonstate/nitrogenstate Summary — WIRED. These set the patch-level
        # veg pool diagnostics: dispvegc/storvegc/totvegc_patch and the N equivalents.
        # Pure diagnostics (no downstream physics reads them); without this they stay
        # at restart-init NaN/0 in history output. (npcropmin/nrepr only matter for the
        # use_crop reproductive pools; CNDriverConfig has no such field -> standard
        # CLM5 literals, moot when use_crop=false.)
        cnveg_carbon_state_summary!(cnveg_cs, mask_bgc_vegp, bounds_patch;
            use_crop=config.use_crop, patch_itype=patch_itype, npcropmin=17, nrepr=NREPR)
        cnveg_nitrogen_state_summary!(cnveg_ns, mask_bgc_vegp, bounds_patch;
            use_crop=config.use_crop, patch_itype=patch_itype, npcropmin=17, nrepr=NREPR)
    end
    # soilbiogeochem carbon/nitrogen STATE summaries — WIRED. Fortran's
    # CNDriverSummarizeStates calls these, and it matters: they produce totc_col /
    # totn_col, which are the BEGIN and END masses of the CN mass-conservation
    # check. Omitting them here (the previous behaviour) left the check with no
    # "before" mass at all, so BeginCNColumnBalance could only ever seed zero — one
    # of the reasons the CN balance check was unusable.
    if cascade_con !== nothing && col !== nothing && patch !== nothing &&
       dzsoi_decomp !== nothing && zisoi_vals !== nothing &&
       nlevdecomp > 0 && ndecomp_pools > 0
        _nc  = length(mask_allc)
        _tvc = fill!(similar(soilbgc_cs.totc_col, _nc), 0)
        _tvn = fill!(similar(soilbgc_ns.totn_col, _nc), 0)
        # TWO DIFFERENT veg aggregates (SoilBiogeochemCarbonStateType.F90:1627-1644):
        #   totecosysc_col (TOTECOSYSC) <- totvegc_col : p2c of totvegc_patch, EXCLUDES cpool
        #   totc_col       (TOTCOLC)    <- totc_p2c_col: p2c of totc_patch,    INCLUDES cpool
        # The CN balance check integrates totc_col, so it needs the cpool-INCLUSIVE
        # aggregate. Passing only totvegc_col left totc_p2c_col at its zeros() default,
        # so totc_col carried no vegetation carbon at all — see the long note at the
        # matching block in clm_driver.jl. This is the BEGIN-of-step summary (it seeds
        # begcb_col), so getting it wrong here while the END-of-step summary was right
        # produced a constant ~totvegc offset (391.7 gC/m2 at Bow) in every step's errcb.
        _tp2c = fill!(similar(soilbgc_cs.totc_col, _nc), 0)
        _tn2c = fill!(similar(soilbgc_ns.totn_col, _nc), 0)
        p2c_1d_filter!(_tvc,  cnveg_cs.totvegc_patch, mask_allc, col, patch)
        p2c_1d_filter!(_tvn,  cnveg_ns.totvegn_patch, mask_allc, col, patch)
        p2c_1d_filter!(_tp2c, cnveg_cs.totc_patch,    mask_allc, col, patch)
        p2c_1d_filter!(_tn2c, cnveg_ns.totn_patch,    mask_allc, col, patch)
        if !isempty(cnveg_cs.totc_p2c_col); copyto!(cnveg_cs.totc_p2c_col, _tp2c); end
        if !isempty(cnveg_ns.totn_p2c_col); copyto!(cnveg_ns.totn_p2c_col, _tn2c); end
        soil_bgc_carbon_state_summary!(soilbgc_cs, mask_allc, bounds_col;
            nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
            dzsoi_decomp_vals=dzsoi_decomp, zisoi_vals=zisoi_vals,
            is_litter=cascade_con.is_litter, is_soil=cascade_con.is_soil,
            is_cwd=cascade_con.is_cwd, totvegc_col=_tvc, totc_p2c_col=_tp2c)
        soil_bgc_nitrogen_state_summary!(soilbgc_ns, mask_allc, bounds_col;
            nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
            dzsoi_decomp_vals=dzsoi_decomp, zisoi_vals=zisoi_vals,
            is_litter=cascade_con.is_litter, is_soil=cascade_con.is_soil,
            is_cwd=cascade_con.is_cwd, totvegn_col=_tvn, totn_p2c_col=_tn2c)

        # Special (non-soil-BGC) columns carry no C/N budget: Fortran zeroes their state
        # (SoilBiogeochemCarbonStateType.F90:686). Leaving them NaN makes the gridcell
        # c2g compute NaN*0.0 == NaN, which silently disables the gridcell balance check
        # (NaN > cerror is false). This is the BEGIN-of-step summary that seeds
        # begcb_grc/begnb_grc, so it needs the same treatment as the END-of-step one.
        for c in bounds_col
            mask_bgc_soilc[c] && continue
            soilbgc_cs.totc_col[c] = 0
            soilbgc_cs.totecosysc_col[c] = 0
            soilbgc_ns.totn_col[c] = 0
            soilbgc_ns.totecosysn_col[c] = 0
        end
    end
    # c13/c14 isotope variants — config-gated off in the default path (not ported).

    return nothing
end

# ---------------------------------------------------------------------------
# cn_driver_summarize_fluxes! — Summarize flux variables
# Ported from CNDriverSummarizeFluxes in CNDriverMod.F90
# ---------------------------------------------------------------------------

"""
    cn_driver_summarize_fluxes!(config; kwargs...)

Call to all CN and SoilBiogeochem summary routines for flux variables.

Ported from `CNDriverSummarizeFluxes` in `CNDriverMod.F90`.

Calls, in Fortran's order:
1. `soil_bgc_carbon_flux_summary!`   — vertically integrates HR / decomp cascade,
   producing `hr_col`, `somhr_col`, `lithr_col`, `cwdhr_col`, `som_c_leached_col`.
2. `soil_bgc_nitrogen_flux_summary!` — produces `denit_col`, `f_n2o_nit_col`,
   `smin_no3_leached_col`, `smin_no3_runoff_col`, `som_n_leached_col`.
3. `cnveg_carbon_flux_summary!` (patch) then `cnveg_carbon_flux_summary_col!`
   (column+gridcell — consumes step 1's `hr_col`).
4. `cnveg_nitrogen_flux_summary!` (patch) then `cnveg_nitrogen_flux_summary_col!`.

Steps 1, 2, and the column halves of 3 and 4 were **ported but never called** —
every column/gridcell C and N flux (`gpp_col`, `er_col`, `hr_col`, `f_n2o_nit_col`,
`nbp_grc`, ...) was structurally dead, which is precisely why the carbon half of
the CN balance check could not run. See `cnveg_carbon_flux_summary_col!`.

The soil-BGC summaries need the decomp cascade metadata (`decomp`, `dzsoi_decomp_vals`);
the column/gridcell halves need the subgrid maps (`col`, `patch`, `lun`, `bounds_grc`)
and the product pools (`c_products`). When those are omitted the corresponding
summary is skipped — that keeps the many unit-test call sites (which pass only the
flux structs) working, and the live driver always supplies them.
"""
function cn_driver_summarize_fluxes!(
        config::CNDriverConfig;
        mask_bgc_soilc::AbstractVector{Bool},
        mask_bgc_vegp::AbstractVector{Bool},
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_nf::CNVegNitrogenFluxData,
        soilbgc_cf::SoilBiogeochemCarbonFluxData,
        soilbgc_nf::SoilBiogeochemNitrogenFluxData,
        patch_itype::Union{AbstractVector{<:Integer},Nothing}=nothing,
        # Subgrid maps — needed for every p2c/c2g. Omitted => column halves skipped.
        col::Union{ColumnData, Nothing} = nothing,
        patch::Union{PatchData, Nothing} = nothing,
        bounds_grc::UnitRange{Int} = 1:0,
        # Decomp cascade metadata — needed by the soil-BGC summaries.
        decomp = nothing,
        dzsoi_decomp_vals::Union{AbstractVector{<:Real}, Nothing} = nothing,
        nlevdecomp::Int = 0,
        ndecomp_pools::Int = 0,
        ndecomp_cascade_transitions::Int = 0,
        # Wood/crop product pools (product_closs_grc feeds landuseflux_grc → nbp_grc).
        c_products::Union{CNProductsData, Nothing} = nothing,
        dt::Real = 0.0,
        nfix_timeconst::Real = 0.0)

    num_bgc_vegp = count(mask_bgc_vegp)
    have_decomp = decomp !== nothing && dzsoi_decomp_vals !== nothing &&
                  nlevdecomp > 0 && ndecomp_pools > 0 && ndecomp_cascade_transitions > 0
    have_subgrid = col !== nothing && patch !== nothing && !isempty(bounds_grc)
    dzv = have_decomp ? collect(Float64.(dzsoi_decomp_vals)) : Float64[]

    # ---- 1. soilbiogeochem_carbonflux_inst%Summary ----
    # (c13/c14 variants: config-gated off in the default path.)
    # `is_microbe` is a MIMICS-only pool flag with no field on DecompCascadeData
    # (decomp_method=1/century is the default); falses() reproduces Fortran there.
    if have_decomp
        soil_bgc_carbon_flux_summary!(soilbgc_cf, BitVector(mask_bgc_soilc), bounds_col;
            ndecomp_cascade_transitions=ndecomp_cascade_transitions,
            nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
            dzsoi_decomp_vals=dzv,
            cascade_donor_pool=decomp.cascade_donor_pool,
            is_soil=decomp.is_soil, is_litter=decomp.is_litter,
            is_cwd=decomp.is_cwd, is_microbe=falses(ndecomp_pools))

        # ---- 2. soilbiogeochem_nitrogenflux_inst%Summary ----
        soil_bgc_nitrogen_flux_summary!(soilbgc_nf, BitVector(mask_bgc_soilc), bounds_col;
            ndecomp_cascade_transitions=ndecomp_cascade_transitions,
            nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
            dzsoi_decomp_vals=dzv,
            use_nitrif_denitrif=config.use_nitrif_denitrif)
    end

    if num_bgc_vegp > 0
        # ---- 3. cnveg_carbonflux_inst%Summary — patch half ----
        # npcropmin/nrepr only matter for the use_crop branch (crop reproductive MR);
        # CNDriverConfig has no such field, so use the standard CLM5 values — harmless
        # when use_crop=false (the default path), correct for crops (npcropmin=17).
        cnveg_carbon_flux_summary!(cnveg_cf, mask_bgc_vegp, bounds_patch;
            use_crop=config.use_crop, use_fun=config.use_fun,
            patch_itype=patch_itype, npcropmin=17, nrepr=NREPR)

        # ---- 3b. column + gridcell half (consumes hr_col from step 1) ----
        if have_decomp && have_subgrid && c_products !== nothing
            cnveg_carbon_flux_summary_col!(cnveg_cf, mask_bgc_soilc, bounds_col,
                bounds_grc, col, patch;
                soilbgc_hr_col=soilbgc_cf.hr_col,
                soilbgc_cwdhr_col=soilbgc_cf.cwdhr_col,
                soilbgc_lithr_col=soilbgc_cf.lithr_col,
                soilbgc_decomp_cascade_ctransfer_col=soilbgc_cf.decomp_cascade_ctransfer_col,
                product_closs_grc=c_products.product_loss_grc,
                nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                ndecomp_cascade_transitions=ndecomp_cascade_transitions,
                dzsoi_decomp_vals=dzv,
                cascade_donor_pool=decomp.cascade_donor_pool,
                is_cwd=decomp.is_cwd, is_litter=decomp.is_litter,
                use_crop=config.use_crop,
                dribble_crophrv_xsmrpool_2atm=config.dribble_crophrv_xsmrpool_2atm)
        end

        # ---- 4. cnveg_nitrogenflux_inst%Summary — patch then column ----
        cnveg_nitrogen_flux_summary!(cnveg_nf, BitVector(mask_bgc_vegp), bounds_patch)
        if have_decomp && have_subgrid
            cnveg_nitrogen_flux_summary_col!(cnveg_nf, mask_bgc_soilc, bounds_col,
                bounds_grc, col, patch;
                nlevdecomp=nlevdecomp, ndecomp_pools=ndecomp_pools,
                dzsoi_decomp_vals=dzv)
        end
    end

    # The exponentially-relaxed lagged NPP. Fortran does this inside
    # CNVegCarbonFluxType::Summary right after the npp_patch p2c
    # (CNVegCarbonFluxType.F90:5368-5385); `npp_col` itself is now p2c'd by
    # cnveg_carbon_flux_summary_col! above, so this only runs the relaxation.
    #
    # `lag_npp_col` is what CNNFixation reads when nfix_timeconst ∈ (0,500) (the
    # CLM5 default, 10 d).
    if col !== nothing && patch !== nothing && !isempty(cnveg_cf.npp_col)
        if !(have_decomp && have_subgrid)
            # Column half did not run (unit-test call shape): still supply npp_col.
            p2c_1d_filter!(cnveg_cf.npp_col, cnveg_cf.npp_patch, mask_bgc_soilc, col, patch)
        end
        cn_lag_npp_update!(cnveg_cf, mask_bgc_soilc;
                           dt = dt, nfix_timeconst = nfix_timeconst)
    end

    return nothing
end

# ---------------------------------------------------------------------------
# Internal helpers: zero flux arrays
# These replace the Fortran SetValues calls that zero flux arrays
# at the start of each timestep.
# ---------------------------------------------------------------------------

"""
    _zero_soilbgc_cflux!(cf; mask, bounds)

Zero soil biogeochem carbon flux fields over masked columns.
Replaces `soilbiogeochem_carbonflux_inst%SetValues(...)` in the Fortran.
"""
# One thread per column; inner loops over (layer, pool) dims taken from size() so
# an unallocated (0-sized) field simply zeros nothing — matches the host length>0 guard.
@kernel function _zero_soilbgc_cf_kernel!(hr_vr, ctransfer_vr, phr_vr, @Const(mask),
        nj_hr::Int, nk_hr::Int, nj_ct::Int, nk_ct::Int, nj_phr::Int, cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        Th = eltype(hr_vr); Tc = eltype(ctransfer_vr); Tp = eltype(phr_vr)
        for k in 1:nk_hr, j in 1:nj_hr; hr_vr[c, j, k] = zero(Th); end
        for k in 1:nk_ct, j in 1:nj_ct; ctransfer_vr[c, j, k] = zero(Tc); end
        for j in 1:nj_phr; phr_vr[c, j] = zero(Tp); end
    end
end

function _zero_soilbgc_cflux!(cf::SoilBiogeochemCarbonFluxData;
                                mask::AbstractVector{Bool},
                                bounds::UnitRange{Int})
    isempty(bounds) && return nothing
    _zero_cnveg_flux_arrays!(cf)   # full per-step reset (ann* accumulators kept)
    return nothing
end

"""
    _zero_cnveg_cflux!(cf; mask, bounds, mask_col, bounds_col)

Zero CN veg carbon flux fields over masked patches and columns.
Replaces `cnveg_carbonflux_inst%SetValues(...)` in the Fortran.
"""
@kernel function _zero_cnveg_cf_kernel!(psnsun, psnshade, @Const(mask), pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        T = eltype(psnsun)
        psnsun[p] = zero(T); psnshade[p] = zero(T)
    end
end

# Zero every per-step flux array field (matches Fortran CNVegCarbonFluxType%SetValues,
# which resets all fluxes each step before the producer routines fill them). Annual
# accumulators (ann*) are preserved. fill! is device-safe; the field iteration is a
# one-shot reset per step. This replaces the earlier 2-field stub — the full driver
# needs ALL fluxes consumed by the c_state_update* cascade reset to 0, otherwise the
# unwired producers leave NaN-init fluxes that poison the state updates.
function _zero_cnveg_flux_arrays!(cf)
    for name in fieldnames(typeof(cf))
        startswith(String(name), "ann") && continue   # keep annual accumulators
        # leafc_to_litter_fun persists step-to-step (Fortran SetValues does NOT
        # zero it): phenology phase-2 sets it, and FUN — which runs BEFORE phase-2
        # — reads the previous step's value for retranslocation accounting.
        name === :leafc_to_litter_fun_patch && continue
        # prev_{leafc,frootc}_to_litter carry the seasonal/stress-deciduous OFFSET
        # LITTERFALL RAMP MEMORY across steps. Fortran's CNVegCarbonFluxType::SetValues
        # does NOT zero them (they are set only in InitCold — CNVegCarbonFluxType.F90:
        # 3986/3988 — and are restart vars, ibid. 4291/4296). They are reset to 0 in
        # exactly two places, both in CNPhenologyMod: at the offset TRIGGER and at the
        # ramp-END cleanup — which phenology.jl already does.
        #
        # Zeroing them here collapsed CNOffsetLitterfall's ramp
        #     leafc_to_litter = prev + t1*(leafc - prev*offset_counter)
        # to just t1*leafc — only ~2-7% of the correct flux mid-ramp. Leaves were then
        # not shed gradually over ndays_off; the pool instead sat ~full for the whole
        # offset period and was dumped in a single timestep by the final-step branch.
        # Invisible in a summer window (offset_flag == 0 all summer); it is the entire
        # story in the autumn leaf-offset window.
        (name === :prev_leafc_to_litter_patch || name === :prev_frootc_to_litter_patch) && continue
        # lag_npp_col is the EXPONENTIALLY-RELAXED NPP memory that CNNFixation reads
        # (nfix_timeconst ∈ (0,500) — the CLM5 default is 10 days). Fortran's
        # CNVegCarbonFluxType::SetValues does NOT zero it: it is spval-init in
        # InitAllocate (CNVegCarbonFluxType.F90:1149), it is a RESTART variable
        # (ibid. 4314), and it is written only by the Summary relaxation
        #     lag_npp = lag_npp*exp(-dt/tau) + npp*(1-exp(-dt/tau))     (ibid. 5375-5381)
        # Blanket-zeroing it here destroyed that memory every single step, collapsing
        # the relaxation to `npp*(1-exp(-dt/tau))` — i.e. ~0.4% of NPP at the 10-day
        # default, with no history at all. The lag is a whole-year smoother; a
        # step-local fraction is a different quantity entirely.
        name === :lag_npp_col && continue
        # Persistent soil-matrix B-input: accumulated across drivers (the dyn_subgrid
        # dwt inputs are added BEFORE this per-step reset) and zeroed only after the
        # soil-matrix solve — so it must survive the step-start reset.
        (name === :matrix_Cinput_col || name === :matrix_Ninput_col) && continue
        arr = getfield(cf, name)
        (arr isa AbstractArray && !isempty(arr) && eltype(arr) <: Real) || continue
        fill!(arr, zero(eltype(arr)))
    end
    return nothing
end

function _zero_cnveg_cflux!(cf::CNVegCarbonFluxData;
                              mask::AbstractVector{Bool},
                              bounds::UnitRange{Int},
                              mask_col::AbstractVector{Bool},
                              bounds_col::UnitRange{Int})
    isempty(bounds) && return nothing
    _zero_cnveg_flux_arrays!(cf)
    return nothing
end

"""
    _zero_cnveg_nflux!(nf; mask, bounds, mask_col, bounds_col)

Zero CN veg nitrogen flux fields over masked patches and columns.
Replaces `cnveg_nitrogenflux_inst%SetValues(...)` in the Fortran.
"""
@kernel function _zero_cnveg_nf_kernel!(plant_ndemand, @Const(mask), pmin::Int, pmax::Int)
    p = @index(Global)
    @inbounds if pmin <= p <= pmax && mask[p]
        plant_ndemand[p] = zero(eltype(plant_ndemand))
    end
end

function _zero_cnveg_nflux!(nf::CNVegNitrogenFluxData;
                              mask::AbstractVector{Bool},
                              bounds::UnitRange{Int},
                              mask_col::AbstractVector{Bool},
                              bounds_col::UnitRange{Int})
    isempty(bounds) && return nothing
    _zero_cnveg_flux_arrays!(nf)   # full per-step reset (ann* accumulators kept)
    return nothing
end

"""
    _zero_soilbgc_nflux!(nf; mask, bounds)

Zero soil biogeochem nitrogen flux fields over masked columns.
Replaces `soilbiogeochem_nitrogenflux_inst%SetValues(...)` in the Fortran.
"""
@kernel function _zero_soilbgc_nf_kernel!(ndep, nfix, fert, soyfixn, ffix, @Const(mask),
        cmin::Int, cmax::Int)
    c = @index(Global)
    @inbounds if cmin <= c <= cmax && mask[c]
        T = eltype(ndep)
        ndep[c] = zero(T); nfix[c] = zero(T); fert[c] = zero(T)
        soyfixn[c] = zero(T); ffix[c] = zero(T)
    end
end

function _zero_soilbgc_nflux!(nf::SoilBiogeochemNitrogenFluxData;
                                mask::AbstractVector{Bool},
                                bounds::UnitRange{Int})
    isempty(bounds) && return nothing
    _zero_cnveg_flux_arrays!(nf)   # full per-step reset incl. leached/runoff (ann* kept)
    return nothing
end
