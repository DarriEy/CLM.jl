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
        mask_bgc_soilc::BitVector,
        mask_bgc_vegp::BitVector,
        mask_pcropp::BitVector = falses(0),
        mask_soilnopcropp::BitVector = falses(0),
        mask_exposedvegp::BitVector = falses(0),
        mask_noexposedvegp::BitVector = falses(0),
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
        patch_column::Vector{Int},
        ivt::Vector{Int},
        woody::Vector{Float64},
        harvdate::Vector{Int},
        col_is_fates::Vector{Bool},
        cascade_donor_pool::Vector{Int},
        cascade_receiver_pool::Vector{Int},
        # Time step
        dt::Float64,
        # Carbon/nitrogen state and flux data structures
        cnveg_cs::CNVegCarbonStateData,
        cnveg_cf::CNVegCarbonFluxData,
        cnveg_ns::CNVegNitrogenStateData,
        cnveg_nf::CNVegNitrogenFluxData,
        soilbgc_cs::SoilBiogeochemCarbonStateData,
        soilbgc_cf::SoilBiogeochemCarbonFluxData,
        soilbgc_ns::SoilBiogeochemNitrogenStateData,
        soilbgc_nf::SoilBiogeochemNitrogenFluxData,
        soilbgc_state::SoilBiogeochemStateData,
        # Output: fire masks (populated by fire routines)
        mask_actfirec::BitVector = falses(length(bounds_col)),
        mask_actfirep::BitVector = falses(length(bounds_patch)))

    num_bgc_vegp = count(mask_bgc_vegp)

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
    # Maintenance respiration — already ported: cn_mresp!(...)
    # (These are called externally before/after this driver)

    # --------------------------------------------------
    # Soil Biogeochemistry
    # --------------------------------------------------
    # Decomposition rate constants — already ported:
    #   decomp_rate_constants_bgc!(...) or decomp_rates_mimics!(...)
    # SoilBiogeochemPotential — not yet ported
    # SoilBiogeochemVerticalProfile — not yet ported
    # SoilBiogeochemNitrifDenitrif — already ported: nitrif_denitrif!(...)

    # --------------------------------------------------
    # Plant/heterotroph competition for mineral N
    # --------------------------------------------------
    if num_bgc_vegp > 0
        # Phenology phase 1 (if use_fun) — already ported: cn_phenology!(...; phase=1)
        # CNFUNInit — already ported: cn_fun_init!(...)
        # Allocation — already ported:
        #   calc_gpp_mr_availc!(...), calc_crop_allocation_fractions!(...), calc_allometry!(...)
    end
    # calc_plant_nutrient_demand — not yet ported
    # SoilBiogeochemCompetition — not yet ported
    # calc_plant_nutrient_competition — not yet ported

    # --------------------------------------------------
    # Actual decomposition
    # --------------------------------------------------
    # SoilBiogeochemDecomp — already ported: soil_biogeochem_decomp!(...)

    # --------------------------------------------------
    # Phenology phase 2
    # --------------------------------------------------
    if num_bgc_vegp > 0
        # CNPhenology phase 1 (if !use_fun) and phase 2 — already ported
    end

    # --------------------------------------------------
    # Growth respiration — already ported: cn_gresp!(...)
    # --------------------------------------------------

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

    # CNPrecisionControl — not yet ported
    # SoilBiogeochemNStateUpdate1 — not yet ported
    # SoilBiogeochemLittVertTransp — not yet ported

    # --------------------------------------------------
    # Gap mortality and Update2
    # --------------------------------------------------
    if num_bgc_vegp > 0
        # CNGapMortality — already ported: cn_gap_mortality!(...), cn_gap_patch_to_column!(...)
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

    # CNPrecisionControl — not yet ported
    # Wood products — not yet ported

    # --------------------------------------------------
    # Fire and Update3
    # --------------------------------------------------
    if num_bgc_vegp > 0
        # CNFireArea, CNFireFluxes — not yet ported
        # CIsoFlux3 — not yet ported

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

    # CNPrecisionControl — not yet ported

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
        dt::Float64,
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
function _zero_soilbgc_cflux!(cf::SoilBiogeochemCarbonFluxData;
                                mask::BitVector,
                                bounds::UnitRange{Int})
    for c in bounds
        mask[c] || continue
        if length(cf.decomp_cascade_hr_vr_col) > 0
            cf.decomp_cascade_hr_vr_col[c, :, :] .= 0.0
        end
        if length(cf.decomp_cascade_ctransfer_vr_col) > 0
            cf.decomp_cascade_ctransfer_vr_col[c, :, :] .= 0.0
        end
        if length(cf.phr_vr_col) > 0
            cf.phr_vr_col[c, :] .= 0.0
        end
    end
    return nothing
end

"""
    _zero_cnveg_cflux!(cf; mask, bounds, mask_col, bounds_col)

Zero CN veg carbon flux fields over masked patches and columns.
Replaces `cnveg_carbonflux_inst%SetValues(...)` in the Fortran.
"""
function _zero_cnveg_cflux!(cf::CNVegCarbonFluxData;
                              mask::BitVector,
                              bounds::UnitRange{Int},
                              mask_col::BitVector,
                              bounds_col::UnitRange{Int})
    for p in bounds
        mask[p] || continue
        cf.psnsun_to_cpool_patch[p] = 0.0
        cf.psnshade_to_cpool_patch[p] = 0.0
    end
    return nothing
end

"""
    _zero_cnveg_nflux!(nf; mask, bounds, mask_col, bounds_col)

Zero CN veg nitrogen flux fields over masked patches and columns.
Replaces `cnveg_nitrogenflux_inst%SetValues(...)` in the Fortran.
"""
function _zero_cnveg_nflux!(nf::CNVegNitrogenFluxData;
                              mask::BitVector,
                              bounds::UnitRange{Int},
                              mask_col::BitVector,
                              bounds_col::UnitRange{Int})
    for p in bounds
        mask[p] || continue
        nf.plant_ndemand_patch[p] = 0.0
    end
    return nothing
end

"""
    _zero_soilbgc_nflux!(nf; mask, bounds)

Zero soil biogeochem nitrogen flux fields over masked columns.
Replaces `soilbiogeochem_nitrogenflux_inst%SetValues(...)` in the Fortran.
"""
function _zero_soilbgc_nflux!(nf::SoilBiogeochemNitrogenFluxData;
                                mask::BitVector,
                                bounds::UnitRange{Int})
    for c in bounds
        mask[c] || continue
        nf.ndep_to_sminn_col[c] = 0.0
        nf.nfix_to_sminn_col[c] = 0.0
        nf.fert_to_sminn_col[c] = 0.0
        nf.soyfixn_to_sminn_col[c] = 0.0
        nf.ffix_to_sminn_col[c] = 0.0
    end
    return nothing
end
