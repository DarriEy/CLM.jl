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

Calls:
- `SoilBiogeochemCompetitionInit` (not yet ported — placeholder)
- `CNPhenologyInit` (via `cn_phenology_init!`)
- `FireInit` (not yet ported — placeholder)

Ported from `CNDriverInit` in `CNDriverMod.F90`.
"""
function cn_driver_init!(config::CNDriverConfig;
                          bounds::UnitRange{Int} = 1:0)
    # SoilBiogeochemCompetitionInit — not yet ported
    # Placeholder: would initialize competition parameters

    if config.use_cn
        # CNPhenologyInit — already ported in phenology.jl
        # In full CLM, called as: cn_phenology_init!(...)
        # Fire initialization — not yet ported
        # Placeholder: cnfire_method%FireInit(bounds, NLFilename)
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
Already-ported state update functions are called directly; subroutines
that require not-yet-ported infrastructure (nutrient competition, fire,
vertical transport, isotope fluxes, etc.) are documented as placeholders.

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
        crop::Union{CropData, Nothing} = nothing,
        photosyns::Union{PhotosynthesisData, Nothing} = nothing,
        canopystate::Union{CanopyStateData, Nothing} = nothing,
        soilstate::Union{SoilStateData, Nothing} = nothing,
        temperature::Union{TemperatureData, Nothing} = nothing,
        water_diag::Union{WaterDiagnosticBulkData, Nothing} = nothing,
        waterstate::Union{WaterStateData, Nothing} = nothing,   # for FUN (h2osoi_liq_col)
        gridcell::Union{GridcellData, Nothing} = nothing,
        is_first_step::Bool = false,
        # Soil-water state (soil-layer slices) for nitrif/denitrif anoxia.
        h2osoi_vol::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        h2osoi_liq::Union{AbstractMatrix{<:Real}, Nothing} = nothing,
        # Output: fire masks (populated by fire routines)
        mask_actfirec::AbstractVector{Bool} = falses(length(bounds_col)),
        mask_actfirep::AbstractVector{Bool} = falses(length(bounds_patch)))

    num_bgc_vegp = count(mask_bgc_vegp)

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
    # N deposition — already ported: n_deposition!(...)
    # N fixation — already ported: n_fixation!(...) or n_free_living_fixation!(...)
    # Crop N fertilization — already ported: n_fert!(...), n_soyfix!(...)

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
        stress_decid = Float64.(pftcon_main.stress_decid))
    if num_bgc_vegp > 0 && cn_shared_params !== nothing && _pftnc !== nothing &&
       patch !== nothing && crop !== nothing && cnveg_state !== nothing
        calc_plant_nutrient_demand!(mask_bgc_vegp, bounds_patch, false,
            _pftnc, cn_shared_params, patch, crop, cnveg_state,
            cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf;
            dt = dt, npcropmin = npcropmin, nrepr = nrepr)
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
            calc_plant_nutrient_competition!(mask_bgc_vegp, bounds_patch,
                _pftnc, cn_shared_params, patch, crop, cnveg_state,
                cnveg_cs, cnveg_cf, cnveg_nf;
                fpg_col = soilbgc_state.fpg_col,
                cnveg_ns = cnveg_ns, dt = dt,
                use_c13 = config.use_c13, use_c14 = config.use_c14,
                npcropmin = npcropmin, nrepr = nrepr)
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
                use_nitrif_denitrif=config.use_nitrif_denitrif)
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
    # NOTE: Plant nutrient demand, soil BGC competition, and plant CN allocation
    # are not yet wired. These require full nutrient competition infrastructure
    # (PatchData, CropData, CNSharedParamsData, SoilBGCCompetitionState/Params).
    # Decomposition proceeds without N limitation in the current implementation.

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
        _phparams = PhenologyParams()
        _phstate  = PhenologyState()
        cn_phenology_init!(_phstate, _phparams, dt)
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
                varctl = varctl, is_first_step = is_first_step, avg_dayspyr = 365.0)
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
    # CIsoFlux1 — not yet ported (carbon isotope flux calculations)
    # --------------------------------------------------

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
            zisoi_vals=zisoi_vals)
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
        # CIsoFlux2 — not yet ported

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

        # Harvest (Update2h)
        # CNHarvest — not yet ported
        # CIsoFlux2h — not yet ported

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

        # Gross unrepresented landcover change (Update2g)
        # CNGrossUnrep — not yet ported
        # CIsoFlux2g — not yet ported

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
    # NOTE: Wood product pools (cn_products) not yet wired — requires harvest fluxes

    # --------------------------------------------------
    # Fire and Update3
    # --------------------------------------------------
    if num_bgc_vegp > 0
        # NOTE: Fire area/fluxes (cnfire_area_li2014!, cnfire_fluxes_li2014!) are
        # implemented but not wired here — requires column temperature data in driver args.
        # C isotope flux phase 3 (CIsoFlux3) not yet ported.

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

        # C14Decay — not yet ported
    end  # num_bgc_vegp > 0

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
        mask_bgc_soilc::BitVector,
        mask_bgc_vegp::BitVector,
        mask_actfirec::BitVector = falses(0),
        mask_actfirep::BitVector = falses(0),
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
        cnveg_nf::CNVegNitrogenFluxData)

    num_bgc_vegp = count(mask_bgc_vegp)

    # N leaching — already ported in n_leaching.jl
    # Called externally: n_leaching!(soilbgc_nf, soilbgc_ns, params; ...)

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

    # Matrix solutions — not yet ported (CNVegMatrix, CNSoilMatrix)

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

Note: The actual Summary methods on the state types are not yet ported.
This function documents the call sequence.
"""
function cn_driver_summarize_states!(
        config::CNDriverConfig;
        mask_bgc_soilc::BitVector,
        mask_bgc_vegp::BitVector,
        mask_allc::BitVector = mask_bgc_soilc,
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        cnveg_cs::CNVegCarbonStateData,
        cnveg_ns::CNVegNitrogenStateData,
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData)

    # cnveg_carbonstate_inst%Summary — not yet ported
    # c13/c14 variants as well
    # soilbiogeochem_carbonstate_inst%summary — not yet ported
    # c13/c14 variants as well
    # cnveg_nitrogenstate_inst%Summary — not yet ported
    # soilbiogeochem_nitrogenstate_inst%summary — not yet ported

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

Note: The actual Summary methods on the flux types are not yet ported.
This function documents the call sequence.
"""
function cn_driver_summarize_fluxes!(
        config::CNDriverConfig;
        mask_bgc_soilc::BitVector,
        mask_bgc_vegp::BitVector,
        bounds_col::UnitRange{Int},
        bounds_patch::UnitRange{Int},
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_nf::CNVegNitrogenFluxData,
        soilbgc_cf::SoilBiogeochemCarbonFluxData,
        soilbgc_nf::SoilBiogeochemNitrogenFluxData)

    num_bgc_vegp = count(mask_bgc_vegp)

    # soilbiogeochem_carbonflux_inst%Summary — not yet ported
    # c13/c14 variants as well
    # soilbiogeochem_nitrogenflux_inst%Summary — not yet ported
    if num_bgc_vegp > 0
        # cnveg_carbonflux_inst%Summary — not yet ported
        # c13/c14 variants as well
        # cnveg_nitrogenflux_inst%Summary — not yet ported
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
