# ==========================================================================
# Ported from: src/biogeochem/CNPhenologyMod.F90
# Phenology routines for coupled carbon-nitrogen code (CN).
# Handles evergreen, seasonal deciduous, stress deciduous, and crop phenology.
#
# Public functions:
#   cn_phenology_set_params!        — Set parameter defaults (for unit tests)
#   cn_phenology_set_nml!           — Set namelist settings (for unit tests)
#   cn_phenology_init!              — Initialization
#   cn_phenology!                   — Main driver (two-phase)
#   crop_phase!                     — Get current crop phase
#   days_past_planting              — Days since planting
#   seasonal_decid_onset            — Seasonal deciduous onset test
#   seasonal_critical_daylength     — Critical daylength for offset
#   get_swindow                     — Get next sowing window
#   was_sown_in_this_window         — Check if crop sown in current window
#
# Private functions:
#   cn_phenology_climate!           — Climate averaging
#   cn_evergreen_phenology!         — Evergreen phenology
#   cn_season_decid_phenology!      — Seasonal deciduous phenology
#   cn_stress_decid_phenology!      — Stress deciduous phenology
#   crop_phenology!                 — Crop phenology
#   crop_phenology_init!            — Crop phenology init
#   plant_crop!                     — Initialize crop at planting
#   vernalization!                  — Vernalization for winter cereal
#   cn_onset_growth!                — Transfer → display during onset
#   cn_offset_litterfall!           — Display → litter during offset
#   cn_background_litterfall!       — Background litterfall
#   cn_livewood_turnover!           — Live wood → dead wood turnover
#   cn_crop_harvest_to_product_pools! — Crop harvest to product pools
#   cn_litter_to_column!            — Patch litter → column level
# ==========================================================================

# ---------------------------------------------------------------------------
# Parameters (replaces Fortran params_type)
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct PhenologyParams
    crit_dayl             ::Float64 = 39200.0    # critical daylength for senescence (s)
    crit_dayl_at_high_lat ::Float64 = 54000.0    # critical daylength at high latitudes (s)
    crit_dayl_lat_slope   ::Float64 = 720.0      # slope of critical daylength with latitude (s/deg)
    ndays_off             ::Float64 = 30.0       # number of days to complete leaf offset
    fstor2tran            ::Float64 = 0.5        # fraction of storage to move to transfer
    crit_onset_fdd        ::Float64 = 15.0       # critical freezing days for onset GDD trigger
    crit_onset_swi        ::Float64 = 15.0       # critical soil water index for onset
    soilpsi_on            ::Float64 = -0.6       # soil water potential for onset (MPa)
    crit_offset_fdd       ::Float64 = 15.0       # critical freezing days for offset trigger
    crit_offset_swi       ::Float64 = 15.0       # critical soil water index for offset
    soilpsi_off           ::Float64 = -0.8       # soil water potential for offset (MPa)
    lwtop                 ::Float64 = 0.7        # live wood turnover proportion (annual)
    phenology_soil_depth  ::Float64 = 0.08       # soil depth for phenology triggers (m)
    snow5d_thresh_for_onset::Float64 = 0.2       # 5-day snow depth threshold for onset (m)
end

# ---------------------------------------------------------------------------
# Module-level state (replaces Fortran module variables)
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct PhenologyState
    dt                    ::Float64 = 0.0        # time step (s)
    fracday               ::Float64 = 0.0        # fraction of day per timestep
    crit_dayl             ::Float64 = 0.0        # critical daylength for offset (s)
    ndays_off             ::Float64 = 0.0        # number of days to complete offset
    fstor2tran            ::Float64 = 0.0        # fraction of storage → transfer
    crit_onset_fdd        ::Float64 = 0.0        # critical freezing days for onset
    crit_onset_swi        ::Float64 = 0.0        # critical soil water index for onset
    soilpsi_on            ::Float64 = 0.0        # water potential for onset (MPa)
    crit_offset_fdd       ::Float64 = 0.0        # critical freezing days for offset
    crit_offset_swi       ::Float64 = 0.0        # critical soil water index for offset
    soilpsi_off           ::Float64 = 0.0        # water potential for offset (MPa)
    lwtop                 ::Float64 = 0.0        # live wood turnover rate (per second)
    phenology_soil_layer  ::Int     = 1          # soil layer index for phenology triggers
    # Crop constants
    p1d                   ::Float64 = 0.004      # photoperiod factor for vernalization
    p1v                   ::Float64 = 0.003      # vernalization factor constant
    hti                   ::Float64 = 1.0        # cold hardening index threshold
    tbase                 ::Float64 = 0.0        # base temperature for vernalization
    # Crop state arrays
    inhemi                ::Vector{Int}     = Int[]          # hemisphere per patch (1=NH, 2=SH)
    minplantjday          ::Matrix{Int}     = Matrix{Int}(undef, 0, 0)  # min planting jday(pft, hemi)
    maxplantjday          ::Matrix{Int}     = Matrix{Int}(undef, 0, 0)  # max planting jday(pft, hemi)
    jdayyrstart           ::Vector{Int}     = [1, 182]       # julian day start of year per hemisphere
end

# ---------------------------------------------------------------------------
# PFT constants needed by phenology (from pftcon)
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct PftConPhenology
    evergreen             ::Vector{Float64} = Float64[]   # binary evergreen flag
    season_decid          ::Vector{Float64} = Float64[]   # binary seasonal deciduous flag
    season_decid_temperate::Vector{Float64} = Float64[]   # binary temperate seasonal deciduous flag
    stress_decid          ::Vector{Float64} = Float64[]   # binary stress deciduous flag
    woody                 ::Vector{Float64} = Float64[]   # binary woody flag
    leaf_long             ::Vector{Float64} = Float64[]   # leaf longevity (yrs)
    leafcn                ::Vector{Float64} = Float64[]   # leaf C:N (gC/gN)
    frootcn               ::Vector{Float64} = Float64[]   # fine root C:N (gC/gN)
    lflitcn               ::Vector{Float64} = Float64[]   # leaf litter C:N (gC/gN)
    livewdcn              ::Vector{Float64} = Float64[]   # live wood C:N (gC/gN)
    deadwdcn              ::Vector{Float64} = Float64[]   # dead wood C:N (gC/gN)
    ndays_on              ::Vector{Float64} = Float64[]   # days to complete onset
    crit_onset_gdd_sf     ::Vector{Float64} = Float64[]   # scale factor for crit_onset_gdd
    lf_f                  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # leaf litter fractions (pft, litr)
    fr_f                  ::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)  # fine root litter fractions (pft, litr)
    biofuel_harvfrac      ::Vector{Float64} = Float64[]   # biofuel harvest fraction
    repr_structure_harvfrac::Matrix{Float64} = Matrix{Float64}(undef, 0, 0) # repr structure harvest frac
    # Crop-specific PFT parameters
    minplanttemp          ::Vector{Float64} = Float64[]
    planttemp             ::Vector{Float64} = Float64[]
    gddmin                ::Vector{Float64} = Float64[]
    lfemerg               ::Vector{Float64} = Float64[]
    grnfill               ::Vector{Float64} = Float64[]
    hybgdd                ::Vector{Float64} = Float64[]
    mxmat                 ::Vector{Int}     = Int[]
    manunitro             ::Vector{Float64} = Float64[]
    is_pft_known_to_model ::Vector{Bool}    = Bool[]
    mnNHplantdate         ::Vector{Float64} = Float64[]
    mxNHplantdate         ::Vector{Float64} = Float64[]
    mnSHplantdate         ::Vector{Float64} = Float64[]
    mxSHplantdate         ::Vector{Float64} = Float64[]
end

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
const NOT_Planted   = 999
const NOT_Harvested = 999
const inNH = 1
const inSH = 2

# Critical daylight method constants
const critical_daylight_constant           = 1
const critical_daylight_depends_on_lat     = 2
const critical_daylight_depends_on_veg     = 3
const critical_daylight_depends_on_latnveg = 4

# Critical offset high latitude (degrees)
const critical_offset_high_lat = 65.0

# Harvest reason constants
const HARVEST_REASON_MATURE          = 1.0
const HARVEST_REASON_MAXSEASLENGTH   = 2.0
const HARVEST_REASON_SOWNBADDEC31    = 3.0
const HARVEST_REASON_SOWTODAY        = 4.0
const HARVEST_REASON_SOWTOMORROW     = 5.0
const HARVEST_REASON_IDOPTOMORROW    = 6.0
const HARVEST_REASON_VERNFREEZEKILL  = 7.0

# Namelist-controlled module variables
const _initial_seed_at_planting = Ref(3.0)
const _onset_thresh_depends_on_veg = Ref(false)
const _critical_daylight_method = Ref(critical_daylight_constant)
const _generate_crop_gdds = Ref(false)
const _use_mxmat = Ref(true)
const _min_gddmaturity = Ref(1.0)
const _min_gdd20_baseline = Ref(0.0)

# Minimum critical daylength for onset (s)
const _min_critical_daylength_onset = 39300.0 / 2.0

# ==========================================================================
# cn_phenology_set_params!  — Set parameter defaults for unit testing
# ==========================================================================
function cn_phenology_set_params!(params::PhenologyParams)
    params.crit_dayl             = 39200.0
    params.crit_dayl_at_high_lat = 54000.0
    params.crit_dayl_lat_slope   = 720.0
    params.ndays_off             = 30.0
    params.fstor2tran            = 0.5
    params.crit_onset_fdd        = 15.0
    params.crit_onset_swi        = 15.0
    params.soilpsi_on            = -0.6
    params.crit_offset_fdd       = 15.0
    params.crit_offset_swi       = 15.0
    params.soilpsi_off           = -0.8
    params.lwtop                 = 0.7
    params.phenology_soil_depth  = 0.08
    params.snow5d_thresh_for_onset = 0.2
    return nothing
end

# ==========================================================================
# cn_phenology_set_nml!  — Set namelist items for unit testing
# ==========================================================================
function cn_phenology_set_nml!(;
        onset_thresh_depends_on_veg::Bool = false,
        critical_daylight_method_in::Int  = critical_daylight_constant)
    _onset_thresh_depends_on_veg[] = onset_thresh_depends_on_veg
    if critical_daylight_method_in < critical_daylight_constant ||
       critical_daylight_method_in > critical_daylight_depends_on_latnveg
        error("ERROR critical_daylight_method out of range")
    end
    _critical_daylight_method[] = critical_daylight_method_in
    return nothing
end

# ==========================================================================
# cn_phenology_init!  — Initialization (after time manager & pftcon ready)
# ==========================================================================
function cn_phenology_init!(pstate::PhenologyState, params::PhenologyParams,
                            dt_in::Float64;
                            use_crop::Bool=false,
                            pftcon::PftConPhenology=PftConPhenology(),
                            patch_data::PatchData=PatchData(),
                            gridcell::GridcellData=GridcellData(),
                            npcropmin::Int=17, npcropmax::Int=78, maxveg::Int=78,
                            find_soil_layer_fn::Function = (depth) -> 1)

    pstate.dt      = dt_in
    pstate.fracday = dt_in / SECSPDAY

    # Set constants from parameters
    pstate.crit_dayl     = params.crit_dayl
    pstate.ndays_off     = params.ndays_off
    pstate.fstor2tran    = params.fstor2tran

    pstate.phenology_soil_layer = find_soil_layer_fn(params.phenology_soil_depth)

    pstate.crit_onset_fdd  = params.crit_onset_fdd
    pstate.crit_onset_swi  = params.crit_onset_swi
    pstate.soilpsi_on      = params.soilpsi_on

    pstate.crit_offset_fdd = params.crit_offset_fdd
    pstate.crit_offset_swi = params.crit_offset_swi
    pstate.soilpsi_off     = params.soilpsi_off

    # Live wood turnover: annual fraction → per second
    pstate.lwtop = params.lwtop / 31536000.0

    # Crop-specific initialization
    if use_crop
        crop_phenology_init!(pstate, pftcon, patch_data, gridcell,
                             npcropmin, npcropmax, maxveg)
    end

    # Error checking for daylength methods
    meth = _critical_daylight_method[]
    if meth == critical_daylight_depends_on_lat ||
       meth == critical_daylight_depends_on_veg ||
       meth == critical_daylight_depends_on_latnveg
        if params.crit_dayl_at_high_lat < params.crit_dayl
            error("ERROR crit_dayl_at_high_lat should be higher than crit_dayl")
        end
        if params.crit_dayl_at_high_lat >= SECSPDAY
            error("ERROR crit_dayl_at_high_lat >= seconds in a day")
        end
        if params.crit_dayl >= SECSPDAY
            error("ERROR crit_dayl >= seconds in a day")
        end
    end
    if meth == critical_daylight_depends_on_lat ||
       meth == critical_daylight_depends_on_latnveg
        if params.crit_dayl_lat_slope <= 0.0
            error("ERROR crit_dayl_lat_slope must be > 0")
        end
    end

    return nothing
end

# ==========================================================================
# cn_phenology!  — Main driver (two-phase)
# ==========================================================================
function cn_phenology!(pstate::PhenologyState, params::PhenologyParams,
                       pftcon::PftConPhenology,
                       mask_soilp::BitVector,   # soil patches
                       mask_pcropp::BitVector,   # prognostic crop patches
                       mask_soilc::BitVector,    # soil columns
                       temperature::TemperatureData,
                       water_diag::WaterDiagnosticBulkData,
                       canopy_state::CanopyStateData,
                       soil_state::SoilStateData,
                       cnveg_state::CNVegStateData,
                       cnveg_cs::CNVegCarbonStateData,
                       cnveg_cf::CNVegCarbonFluxData,
                       cnveg_ns::CNVegNitrogenStateData,
                       cnveg_nf::CNVegNitrogenFluxData,
                       crop::CropData,
                       patch_data::PatchData,
                       gridcell::GridcellData,
                       cn_params::CNSharedParamsData,
                       leaf_prof_patch::Matrix{Float64},
                       froot_prof_patch::Matrix{Float64},
                       phase::Int;
                       varctl::VarCtl=VarCtl(),
                       is_first_step::Bool=false,
                       avg_dayspyr::Float64=365.0)

    if phase == 1
        cn_phenology_climate!(pstate, mask_soilp, temperature, cnveg_state, crop,
                              patch_data, pftcon)

        cn_evergreen_phenology!(pstate, mask_soilp, pftcon, cnveg_state,
                                cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                                patch_data, avg_dayspyr)

        cn_season_decid_phenology!(pstate, params, mask_soilp, pftcon,
                                   temperature, water_diag, cnveg_state,
                                   cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                                   patch_data, gridcell; use_cndv=varctl.use_cndv)

        cn_stress_decid_phenology!(pstate, mask_soilp, pftcon,
                                   soil_state, temperature,
                                   cnveg_state, cnveg_cs, cnveg_ns,
                                   cnveg_cf, cnveg_nf,
                                   patch_data, gridcell, cn_params;
                                   avg_dayspyr=avg_dayspyr)

        if any(mask_pcropp) && !is_first_step
            crop_phenology!(pstate, params, mask_pcropp, pftcon,
                            water_diag, temperature, crop, canopy_state,
                            cnveg_state, cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                            patch_data, gridcell;
                            varctl=varctl, avg_dayspyr=avg_dayspyr)
        end

    elseif phase == 2
        cn_onset_growth!(pstate, mask_soilp, pftcon,
                         cnveg_state, cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                         patch_data)

        cn_offset_litterfall!(pstate, mask_soilp, pftcon,
                              cnveg_state, cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                              crop, patch_data; use_fun=cn_params.use_fun,
                              CNratio_floating=varctl.CNratio_floating,
                              for_testing_no_crop_seed_replenishment=varctl.for_testing_no_crop_seed_replenishment)

        cn_background_litterfall!(pstate, mask_soilp, pftcon,
                                  cnveg_state, cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                                  patch_data; use_fun=cn_params.use_fun,
                                  CNratio_floating=varctl.CNratio_floating)

        cn_livewood_turnover!(pstate, mask_soilp, pftcon,
                              cnveg_cs, cnveg_ns, cnveg_cf, cnveg_nf,
                              patch_data; CNratio_floating=varctl.CNratio_floating)

        cn_crop_harvest_to_product_pools!(mask_soilp, mask_soilc,
                                          cnveg_cf, cnveg_nf, patch_data;
                                          use_crop=varctl.use_crop)

        cn_litter_to_column!(mask_soilp, pftcon,
                             cnveg_state, cnveg_cf, cnveg_nf,
                             patch_data, leaf_prof_patch, froot_prof_patch;
                             use_grainproduct=false)
    else
        error("bad phase: $phase")
    end

    return nothing
end

# ==========================================================================
# cn_phenology_climate!  — Climate averaging for phenology triggers
# ==========================================================================
function cn_phenology_climate!(pstate::PhenologyState,
                               mask_soilp::BitVector,
                               temperature::TemperatureData,
                               cnveg_state::CNVegStateData,
                               crop::CropData,
                               patch_data::PatchData,
                               pftcon::PftConPhenology)
    dt      = pstate.dt
    fracday = pstate.fracday

    for p in eachindex(mask_soilp)
        mask_soilp[p] || continue

        ivt = patch_data.itype[p]

        # Update tempavg_t2m accumulator
        # This is a 10-day running mean accumulator update
        if temperature.t_ref2m_patch[p] != SPVAL
            cnveg_state.tempavg_t2m_patch[p] = cnveg_state.tempavg_t2m_patch[p] +
                fracday * (temperature.t_ref2m_patch[p] - cnveg_state.tempavg_t2m_patch[p])
        end
    end

    return nothing
end

# ==========================================================================
# cn_evergreen_phenology!
# ==========================================================================
function cn_evergreen_phenology!(pstate::PhenologyState,
                                 mask_soilp::BitVector,
                                 pftcon::PftConPhenology,
                                 cnveg_state::CNVegStateData,
                                 cnveg_cs::CNVegCarbonStateData,
                                 cnveg_ns::CNVegNitrogenStateData,
                                 cnveg_cf::CNVegCarbonFluxData,
                                 cnveg_nf::CNVegNitrogenFluxData,
                                 patch_data::PatchData,
                                 avg_dayspyr::Float64)
    dt         = pstate.dt
    fstor2tran = pstate.fstor2tran

    for p in eachindex(mask_soilp)
        mask_soilp[p] || continue

        ivt = patch_data.itype[p]
        if pftcon.evergreen[ivt] == 1.0
            # set background litterfall rate
            cnveg_state.bglfr_patch[p] = 1.0 / (pftcon.leaf_long[ivt] * avg_dayspyr * SECSPDAY)
            cnveg_state.bgtr_patch[p]  = 0.0
            cnveg_state.lgsf_patch[p]  = 0.0

            # move storage pools to transfer pools
            cnveg_cf.leafc_storage_to_xfer_patch[p]  = fstor2tran * cnveg_cs.leafc_storage_patch[p] / dt
            cnveg_cf.frootc_storage_to_xfer_patch[p] = fstor2tran * cnveg_cs.frootc_storage_patch[p] / dt

            if pftcon.woody[ivt] == 1.0
                cnveg_cf.livestemc_storage_to_xfer_patch[p]  = fstor2tran * cnveg_cs.livestemc_storage_patch[p] / dt
                cnveg_cf.deadstemc_storage_to_xfer_patch[p]  = fstor2tran * cnveg_cs.deadstemc_storage_patch[p] / dt
                cnveg_cf.livecrootc_storage_to_xfer_patch[p] = fstor2tran * cnveg_cs.livecrootc_storage_patch[p] / dt
                cnveg_cf.deadcrootc_storage_to_xfer_patch[p] = fstor2tran * cnveg_cs.deadcrootc_storage_patch[p] / dt
                cnveg_cf.gresp_storage_to_xfer_patch[p]      = fstor2tran * cnveg_cs.gresp_storage_patch[p] / dt
            end

            cnveg_nf.leafn_storage_to_xfer_patch[p]  = fstor2tran * cnveg_ns.leafn_storage_patch[p] / dt
            cnveg_nf.frootn_storage_to_xfer_patch[p] = fstor2tran * cnveg_ns.frootn_storage_patch[p] / dt
            if pftcon.woody[ivt] == 1.0
                cnveg_nf.livestemn_storage_to_xfer_patch[p]  = fstor2tran * cnveg_ns.livestemn_storage_patch[p] / dt
                cnveg_nf.deadstemn_storage_to_xfer_patch[p]  = fstor2tran * cnveg_ns.deadstemn_storage_patch[p] / dt
                cnveg_nf.livecrootn_storage_to_xfer_patch[p] = fstor2tran * cnveg_ns.livecrootn_storage_patch[p] / dt
                cnveg_nf.deadcrootn_storage_to_xfer_patch[p] = fstor2tran * cnveg_ns.deadcrootn_storage_patch[p] / dt
            end
        end
    end

    return nothing
end

# ==========================================================================
# cn_season_decid_phenology!
# ==========================================================================
function cn_season_decid_phenology!(pstate::PhenologyState,
                                    params::PhenologyParams,
                                    mask_soilp::BitVector,
                                    pftcon::PftConPhenology,
                                    temperature::TemperatureData,
                                    water_diag::WaterDiagnosticBulkData,
                                    cnveg_state::CNVegStateData,
                                    cnveg_cs::CNVegCarbonStateData,
                                    cnveg_ns::CNVegNitrogenStateData,
                                    cnveg_cf::CNVegCarbonFluxData,
                                    cnveg_nf::CNVegNitrogenFluxData,
                                    patch_data::PatchData,
                                    gridcell::GridcellData;
                                    use_cndv::Bool=false)
    dt         = pstate.dt
    fracday    = pstate.fracday
    crit_dayl  = pstate.crit_dayl
    ndays_off  = pstate.ndays_off
    fstor2tran = pstate.fstor2tran
    soil_layer = pstate.phenology_soil_layer

    for p in eachindex(mask_soilp)
        mask_soilp[p] || continue

        ivt = patch_data.itype[p]
        c   = patch_data.column[p]
        g   = patch_data.gridcell[p]

        if pftcon.season_decid[ivt] != 1.0
            continue
        end

        # set background rates to 0
        cnveg_state.bglfr_patch[p] = 0.0
        cnveg_state.bgtr_patch[p]  = 0.0
        cnveg_state.lgsf_patch[p]  = 0.0

        # onset GDD sum
        crit_onset_gdd = pftcon.crit_onset_gdd_sf[ivt] *
            exp(4.8 + 0.13 * (cnveg_state.annavg_t2m_patch[p] - TFRZ))

        # winter→summer flag
        ws_flag = gridcell.dayl[g] >= gridcell.prev_dayl[g] ? 1.0 : 0.0

        # update offset_counter
        if cnveg_state.offset_flag_patch[p] == 1.0
            cnveg_state.offset_counter_patch[p] -= dt

            if cnveg_state.offset_counter_patch[p] < dt / 2.0
                cnveg_state.offset_flag_patch[p]    = 0.0
                cnveg_state.offset_counter_patch[p]  = 0.0
                cnveg_state.dormant_flag_patch[p]    = 1.0
                cnveg_state.days_active_patch[p]     = 0.0

                cnveg_cf.prev_leafc_to_litter_patch[p]  = 0.0
                cnveg_cf.prev_frootc_to_litter_patch[p] = 0.0
            end
        end

        # update onset_counter
        if cnveg_state.onset_flag_patch[p] == 1.0
            cnveg_state.onset_counter_patch[p] -= dt

            if cnveg_state.onset_counter_patch[p] < dt / 2.0
                cnveg_state.onset_flag_patch[p]    = 0.0
                cnveg_state.onset_counter_patch[p] = 0.0
                # zero transfer growth rates
                cnveg_cf.leafc_xfer_to_leafc_patch[p]   = 0.0
                cnveg_cf.frootc_xfer_to_frootc_patch[p] = 0.0
                cnveg_nf.leafn_xfer_to_leafn_patch[p]   = 0.0
                cnveg_nf.frootn_xfer_to_frootn_patch[p] = 0.0
                if pftcon.woody[ivt] == 1.0
                    cnveg_cf.livestemc_xfer_to_livestemc_patch[p]   = 0.0
                    cnveg_cf.deadstemc_xfer_to_deadstemc_patch[p]   = 0.0
                    cnveg_cf.livecrootc_xfer_to_livecrootc_patch[p] = 0.0
                    cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch[p] = 0.0
                    cnveg_nf.livestemn_xfer_to_livestemn_patch[p]   = 0.0
                    cnveg_nf.deadstemn_xfer_to_deadstemn_patch[p]   = 0.0
                    cnveg_nf.livecrootn_xfer_to_livecrootn_patch[p] = 0.0
                    cnveg_nf.deadcrootn_xfer_to_deadcrootn_patch[p] = 0.0
                end
                # zero transfer pools
                cnveg_cs.leafc_xfer_patch[p]  = 0.0
                cnveg_ns.leafn_xfer_patch[p]  = 0.0
                cnveg_cs.frootc_xfer_patch[p] = 0.0
                cnveg_ns.frootn_xfer_patch[p] = 0.0
                if pftcon.woody[ivt] == 1.0
                    cnveg_cs.livestemc_xfer_patch[p]  = 0.0
                    cnveg_ns.livestemn_xfer_patch[p]  = 0.0
                    cnveg_cs.deadstemc_xfer_patch[p]  = 0.0
                    cnveg_ns.deadstemn_xfer_patch[p]  = 0.0
                    cnveg_cs.livecrootc_xfer_patch[p] = 0.0
                    cnveg_ns.livecrootn_xfer_patch[p] = 0.0
                    cnveg_cs.deadcrootc_xfer_patch[p] = 0.0
                    cnveg_ns.deadcrootn_xfer_patch[p] = 0.0
                end
            end
        end

        # test dormant → growth
        if cnveg_state.dormant_flag_patch[p] == 1.0
            soilt = temperature.t_soisno_col[c, soil_layer]
            snow_5day = water_diag.snow_5day_col[c]
            soila10 = temperature.soila10_col[c]
            t_a5min = temperature.t_a5min_patch[p]

            do_onset = seasonal_decid_onset(
                cnveg_state.onset_gdd_patch[p],
                cnveg_state.onset_gddflag_patch[p],
                soilt, soila10, t_a5min, gridcell.dayl[g],
                snow_5day, ws_flag, crit_onset_gdd,
                pftcon.season_decid_temperate[ivt],
                fracday, params.snow5d_thresh_for_onset)

            # update in/out args
            cnveg_state.onset_gdd_patch[p]     = do_onset.onset_gdd
            cnveg_state.onset_gddflag_patch[p]  = do_onset.onset_gddflag

            if do_onset.result
                cnveg_state.onset_flag_patch[p]     = 1.0
                cnveg_state.dormant_flag_patch[p]    = 0.0
                cnveg_state.onset_gddflag_patch[p]   = 0.0
                cnveg_state.onset_gdd_patch[p]       = 0.0
                cnveg_state.onset_counter_patch[p]   = pftcon.ndays_on[ivt] * SECSPDAY

                # storage → transfer (non-matrix)
                cnveg_cf.leafc_storage_to_xfer_patch[p]  = fstor2tran * cnveg_cs.leafc_storage_patch[p] / dt
                cnveg_cf.frootc_storage_to_xfer_patch[p] = fstor2tran * cnveg_cs.frootc_storage_patch[p] / dt
                if pftcon.woody[ivt] == 1.0
                    cnveg_cf.livestemc_storage_to_xfer_patch[p]  = fstor2tran * cnveg_cs.livestemc_storage_patch[p] / dt
                    cnveg_cf.deadstemc_storage_to_xfer_patch[p]  = fstor2tran * cnveg_cs.deadstemc_storage_patch[p] / dt
                    cnveg_cf.livecrootc_storage_to_xfer_patch[p] = fstor2tran * cnveg_cs.livecrootc_storage_patch[p] / dt
                    cnveg_cf.deadcrootc_storage_to_xfer_patch[p] = fstor2tran * cnveg_cs.deadcrootc_storage_patch[p] / dt
                    cnveg_cf.gresp_storage_to_xfer_patch[p]      = fstor2tran * cnveg_cs.gresp_storage_patch[p] / dt
                end
                cnveg_nf.leafn_storage_to_xfer_patch[p]  = fstor2tran * cnveg_ns.leafn_storage_patch[p] / dt
                cnveg_nf.frootn_storage_to_xfer_patch[p] = fstor2tran * cnveg_ns.frootn_storage_patch[p] / dt
                if pftcon.woody[ivt] == 1.0
                    cnveg_nf.livestemn_storage_to_xfer_patch[p]  = fstor2tran * cnveg_ns.livestemn_storage_patch[p] / dt
                    cnveg_nf.deadstemn_storage_to_xfer_patch[p]  = fstor2tran * cnveg_ns.deadstemn_storage_patch[p] / dt
                    cnveg_nf.livecrootn_storage_to_xfer_patch[p] = fstor2tran * cnveg_ns.livecrootn_storage_patch[p] / dt
                    cnveg_nf.deadcrootn_storage_to_xfer_patch[p] = fstor2tran * cnveg_ns.deadcrootn_storage_patch[p] / dt
                end
            end

        # test growth → offset
        elseif cnveg_state.offset_flag_patch[p] == 0.0
            if use_cndv
                cnveg_state.days_active_patch[p] += fracday
            end

            crit_daylat = seasonal_critical_daylength(
                g, p, pstate.crit_dayl, params, pftcon, gridcell, patch_data)

            if ws_flag == 0.0 && gridcell.dayl[g] < crit_daylat
                cnveg_state.offset_flag_patch[p]    = 1.0
                cnveg_state.offset_counter_patch[p]  = ndays_off * SECSPDAY
                cnveg_cf.prev_leafc_to_litter_patch[p]  = 0.0
                cnveg_cf.prev_frootc_to_litter_patch[p] = 0.0
            end
        end

    end # patch loop

    return nothing
end

# ==========================================================================
# seasonal_critical_daylength — Critical daylength for seasonal deciduous offset
# ==========================================================================
function seasonal_critical_daylength(g::Int, p::Int,
                                     crit_dayl_val::Float64,
                                     params::PhenologyParams,
                                     pftcon::PftConPhenology,
                                     gridcell::GridcellData,
                                     patch_data::PatchData)
    method = _critical_daylight_method[]
    ivt = patch_data.itype[p]

    if method == critical_daylight_depends_on_latnveg
        if pftcon.season_decid_temperate[ivt] == 1.0
            return crit_dayl_val
        else
            cd = params.crit_dayl_at_high_lat -
                 params.crit_dayl_lat_slope * (critical_offset_high_lat - abs(gridcell.latdeg[g]))
            return max(cd, crit_dayl_val)
        end
    elseif method == critical_daylight_depends_on_veg
        if pftcon.season_decid_temperate[ivt] == 1.0
            return crit_dayl_val
        else
            return params.crit_dayl_at_high_lat
        end
    elseif method == critical_daylight_depends_on_lat
        cd = params.crit_dayl_at_high_lat -
             params.crit_dayl_lat_slope * (critical_offset_high_lat - abs(gridcell.latdeg[g]))
        return max(cd, crit_dayl_val)
    elseif method == critical_daylight_constant
        return crit_dayl_val
    else
        error("ERROR: critical_daylight_method not implemented")
    end
end

# ==========================================================================
# seasonal_decid_onset — Determine if seasonal deciduous onset should happen
# ==========================================================================
function seasonal_decid_onset(onset_gdd::Float64, onset_gddflag::Float64,
                              soilt::Float64, soila10::Float64, t_a5min::Float64,
                              dayl::Float64, snow_5day::Float64,
                              ws_flag::Float64, crit_onset_gdd::Float64,
                              season_decid_temperate::Float64,
                              fracday::Float64, snow5d_thresh::Float64)
    og      = onset_gdd
    ogf     = onset_gddflag
    result  = false

    # switch on GDD sum at winter solstice
    if ogf == 0.0 && ws_flag == 1.0
        ogf = 1.0
        og  = 0.0
    end

    # reset if past summer solstice without reaching threshold
    if ogf == 1.0 && ws_flag == 0.0
        ogf = 0.0
        og  = 0.0
    end

    # accumulate GDD
    if ogf == 1.0 && soilt > TFRZ
        og += (soilt - TFRZ) * fracday
    end

    if _onset_thresh_depends_on_veg[]
        if og > crit_onset_gdd && season_decid_temperate == 1.0
            result = true
        elseif season_decid_temperate == 0.0 && ogf == 1.0 &&
               soila10 > TFRZ && t_a5min > TFRZ && ws_flag == 1.0 &&
               dayl > _min_critical_daylength_onset &&
               snow_5day < snow5d_thresh
            result = true
        end
    else
        if og > crit_onset_gdd
            result = true
        end
    end

    return (result=result, onset_gdd=og, onset_gddflag=ogf)
end

# ==========================================================================
# cn_stress_decid_phenology!
# ==========================================================================
function cn_stress_decid_phenology!(pstate::PhenologyState,
                                    mask_soilp::BitVector,
                                    pftcon::PftConPhenology,
                                    soil_state::SoilStateData,
                                    temperature::TemperatureData,
                                    cnveg_state::CNVegStateData,
                                    cnveg_cs::CNVegCarbonStateData,
                                    cnveg_ns::CNVegNitrogenStateData,
                                    cnveg_cf::CNVegCarbonFluxData,
                                    cnveg_nf::CNVegNitrogenFluxData,
                                    patch_data::PatchData,
                                    gridcell::GridcellData,
                                    cn_params::CNSharedParamsData;
                                    avg_dayspyr::Float64=365.0)
    dt              = pstate.dt
    fracday         = pstate.fracday
    ndays_off_val   = pstate.ndays_off
    fstor2tran      = pstate.fstor2tran
    soil_layer      = pstate.phenology_soil_layer
    crit_onset_fdd  = pstate.crit_onset_fdd
    crit_onset_swi  = pstate.crit_onset_swi
    soilpsi_on      = pstate.soilpsi_on
    crit_offset_fdd = pstate.crit_offset_fdd
    crit_offset_swi = pstate.crit_offset_swi
    soilpsi_off     = pstate.soilpsi_off
    secspqtrday     = SECSPDAY / 4.0
    rain_threshold  = 20.0

    for p in eachindex(mask_soilp)
        mask_soilp[p] || continue

        ivt = patch_data.itype[p]
        c   = patch_data.column[p]
        g   = patch_data.gridcell[p]

        if pftcon.stress_decid[ivt] != 1.0
            continue
        end

        soilt = temperature.t_soisno_col[c, soil_layer]
        psi   = soil_state.soilpsi_col[c, soil_layer]

        crit_onset_gdd = pftcon.crit_onset_gdd_sf[ivt] *
            exp(4.8 + 0.13 * (cnveg_state.annavg_t2m_patch[p] - TFRZ))

        # offset counter
        if cnveg_state.offset_flag_patch[p] == 1.0
            cnveg_state.offset_counter_patch[p] -= dt
            if cnveg_state.offset_counter_patch[p] < dt / 2.0
                cnveg_state.offset_flag_patch[p]    = 0.0
                cnveg_state.offset_counter_patch[p]  = 0.0
                cnveg_state.dormant_flag_patch[p]    = 1.0
                cnveg_state.days_active_patch[p]     = 0.0
                cnveg_cf.prev_leafc_to_litter_patch[p]  = 0.0
                cnveg_cf.prev_frootc_to_litter_patch[p] = 0.0
            end
        end

        # onset counter
        if cnveg_state.onset_flag_patch[p] == 1.0
            cnveg_state.onset_counter_patch[p] -= dt
            if cnveg_state.onset_counter_patch[p] < dt / 2.0
                cnveg_state.onset_flag_patch[p]    = 0.0
                cnveg_state.onset_counter_patch[p] = 0.0
                cnveg_cf.leafc_xfer_to_leafc_patch[p]   = 0.0
                cnveg_cf.frootc_xfer_to_frootc_patch[p] = 0.0
                cnveg_nf.leafn_xfer_to_leafn_patch[p]   = 0.0
                cnveg_nf.frootn_xfer_to_frootn_patch[p] = 0.0
                if pftcon.woody[ivt] == 1.0
                    cnveg_cf.livestemc_xfer_to_livestemc_patch[p]   = 0.0
                    cnveg_cf.deadstemc_xfer_to_deadstemc_patch[p]   = 0.0
                    cnveg_cf.livecrootc_xfer_to_livecrootc_patch[p] = 0.0
                    cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch[p] = 0.0
                    cnveg_nf.livestemn_xfer_to_livestemn_patch[p]   = 0.0
                    cnveg_nf.deadstemn_xfer_to_deadstemn_patch[p]   = 0.0
                    cnveg_nf.livecrootn_xfer_to_livecrootn_patch[p] = 0.0
                    cnveg_nf.deadcrootn_xfer_to_deadcrootn_patch[p] = 0.0
                end
                cnveg_cs.leafc_xfer_patch[p]  = 0.0
                cnveg_ns.leafn_xfer_patch[p]  = 0.0
                cnveg_cs.frootc_xfer_patch[p] = 0.0
                cnveg_ns.frootn_xfer_patch[p] = 0.0
                if pftcon.woody[ivt] == 1.0
                    cnveg_cs.livestemc_xfer_patch[p]  = 0.0
                    cnveg_ns.livestemn_xfer_patch[p]  = 0.0
                    cnveg_cs.deadstemc_xfer_patch[p]  = 0.0
                    cnveg_ns.deadstemn_xfer_patch[p]  = 0.0
                    cnveg_cs.livecrootc_xfer_patch[p] = 0.0
                    cnveg_ns.livecrootn_xfer_patch[p] = 0.0
                    cnveg_cs.deadcrootc_xfer_patch[p] = 0.0
                    cnveg_ns.deadcrootn_xfer_patch[p] = 0.0
                end
            end
        end

        # test dormant → growth
        if cnveg_state.dormant_flag_patch[p] == 1.0
            # freezing degree days
            if cnveg_state.onset_gddflag_patch[p] == 0.0 && soilt < TFRZ
                cnveg_state.onset_fdd_patch[p] += fracday
            end
            if cnveg_state.onset_fdd_patch[p] > crit_onset_fdd
                cnveg_state.onset_gddflag_patch[p] = 1.0
                cnveg_state.onset_fdd_patch[p]     = 0.0
                cnveg_state.onset_swi_patch[p]     = 0.0
            end
            # GDD accumulation
            if cnveg_state.onset_gddflag_patch[p] == 1.0 && soilt > TFRZ
                cnveg_state.onset_gdd_patch[p] += (soilt - TFRZ) * fracday
            end
            # soil water index
            additional_onset_condition = true
            if cn_params.constrain_stress_deciduous_onset
                if hasfield(typeof(cnveg_state), :prec10_patch) &&
                   length(cnveg_state.prec10_patch) >= p
                    # prec10 check
                else
                    # skip (no prec10 available yet)
                end
            end
            if psi >= soilpsi_on
                cnveg_state.onset_swi_patch[p] += fracday
            end
            if cnveg_state.onset_swi_patch[p] > crit_onset_swi && additional_onset_condition
                cnveg_state.onset_flag_patch[p] = 1.0
                if cnveg_state.onset_gddflag_patch[p] == 1.0 &&
                   cnveg_state.onset_gdd_patch[p] < crit_onset_gdd
                    cnveg_state.onset_flag_patch[p] = 0.0
                end
            end
            # only allow onset if dayl > 6hrs
            if cnveg_state.onset_flag_patch[p] == 1.0 && gridcell.dayl[g] <= secspqtrday
                cnveg_state.onset_flag_patch[p] = 0.0
            end
            # if onset triggered, reset and do storage→transfer
            if cnveg_state.onset_flag_patch[p] == 1.0
                cnveg_state.dormant_flag_patch[p]    = 0.0
                cnveg_state.days_active_patch[p]     = 0.0
                cnveg_state.onset_gddflag_patch[p]   = 0.0
                cnveg_state.onset_fdd_patch[p]       = 0.0
                cnveg_state.onset_gdd_patch[p]       = 0.0
                cnveg_state.onset_swi_patch[p]       = 0.0
                cnveg_state.onset_counter_patch[p]   = pftcon.ndays_on[ivt] * SECSPDAY

                cnveg_cf.leafc_storage_to_xfer_patch[p]  = fstor2tran * cnveg_cs.leafc_storage_patch[p] / dt
                cnveg_cf.frootc_storage_to_xfer_patch[p] = fstor2tran * cnveg_cs.frootc_storage_patch[p] / dt
                if pftcon.woody[ivt] == 1.0
                    cnveg_cf.livestemc_storage_to_xfer_patch[p]  = fstor2tran * cnveg_cs.livestemc_storage_patch[p] / dt
                    cnveg_cf.deadstemc_storage_to_xfer_patch[p]  = fstor2tran * cnveg_cs.deadstemc_storage_patch[p] / dt
                    cnveg_cf.livecrootc_storage_to_xfer_patch[p] = fstor2tran * cnveg_cs.livecrootc_storage_patch[p] / dt
                    cnveg_cf.deadcrootc_storage_to_xfer_patch[p] = fstor2tran * cnveg_cs.deadcrootc_storage_patch[p] / dt
                    cnveg_cf.gresp_storage_to_xfer_patch[p]      = fstor2tran * cnveg_cs.gresp_storage_patch[p] / dt
                end
                cnveg_nf.leafn_storage_to_xfer_patch[p]  = fstor2tran * cnveg_ns.leafn_storage_patch[p] / dt
                cnveg_nf.frootn_storage_to_xfer_patch[p] = fstor2tran * cnveg_ns.frootn_storage_patch[p] / dt
                if pftcon.woody[ivt] == 1.0
                    cnveg_nf.livestemn_storage_to_xfer_patch[p]  = fstor2tran * cnveg_ns.livestemn_storage_patch[p] / dt
                    cnveg_nf.deadstemn_storage_to_xfer_patch[p]  = fstor2tran * cnveg_ns.deadstemn_storage_patch[p] / dt
                    cnveg_nf.livecrootn_storage_to_xfer_patch[p] = fstor2tran * cnveg_ns.livecrootn_storage_patch[p] / dt
                    cnveg_nf.deadcrootn_storage_to_xfer_patch[p] = fstor2tran * cnveg_ns.deadcrootn_storage_patch[p] / dt
                end
            end

        # test growth → offset
        elseif cnveg_state.offset_flag_patch[p] == 0.0
            if psi <= soilpsi_off
                cnveg_state.offset_swi_patch[p] += fracday
                if cnveg_state.offset_swi_patch[p] >= crit_offset_swi &&
                   cnveg_state.onset_flag_patch[p] == 0.0
                    cnveg_state.offset_flag_patch[p] = 1.0
                end
            elseif psi >= soilpsi_on
                cnveg_state.offset_swi_patch[p] -= fracday
                cnveg_state.offset_swi_patch[p] = max(cnveg_state.offset_swi_patch[p], 0.0)
            end
            # freezing day accumulation
            if cnveg_state.offset_fdd_patch[p] > 0.0 && soilt > TFRZ
                cnveg_state.offset_fdd_patch[p] -= fracday
                cnveg_state.offset_fdd_patch[p] = max(0.0, cnveg_state.offset_fdd_patch[p])
            end
            if soilt <= TFRZ
                cnveg_state.offset_fdd_patch[p] += fracday
                if cnveg_state.offset_fdd_patch[p] > crit_offset_fdd &&
                   cnveg_state.onset_flag_patch[p] == 0.0
                    cnveg_state.offset_flag_patch[p] = 1.0
                end
            end
            # force offset if dayl < 6 hrs
            if gridcell.dayl[g] <= secspqtrday
                cnveg_state.offset_flag_patch[p] = 1.0
            end
            # set offset params
            if cnveg_state.offset_flag_patch[p] == 1.0
                cnveg_state.offset_fdd_patch[p]     = 0.0
                cnveg_state.offset_swi_patch[p]     = 0.0
                cnveg_state.offset_counter_patch[p]  = ndays_off_val * SECSPDAY
                cnveg_cf.prev_leafc_to_litter_patch[p]  = 0.0
                cnveg_cf.prev_frootc_to_litter_patch[p] = 0.0
            end
        end

        # days active tracking
        if cnveg_state.dormant_flag_patch[p] == 0.0
            cnveg_state.days_active_patch[p] += fracday
        end

        # long growing season factor
        cnveg_state.lgsf_patch[p] = max(min(
            3.0 * (cnveg_state.days_active_patch[p] - pftcon.leaf_long[ivt] * avg_dayspyr) / avg_dayspyr,
            1.0), 0.0)

        # background litterfall rate
        if cnveg_state.offset_flag_patch[p] == 1.0
            cnveg_state.bglfr_patch[p] = 0.0
        else
            cnveg_state.bglfr_patch[p] = (1.0 / (pftcon.leaf_long[ivt] * avg_dayspyr * SECSPDAY)) *
                cnveg_state.lgsf_patch[p]
        end

        # background transfer rate
        if cnveg_state.onset_flag_patch[p] == 1.0
            cnveg_state.bgtr_patch[p] = 0.0
        else
            bgtr_val = (1.0 / (avg_dayspyr * SECSPDAY)) * cnveg_state.lgsf_patch[p]
            cnveg_state.bgtr_patch[p] = bgtr_val

            # non-matrix storage → transfer
            cnveg_cf.leafc_storage_to_xfer_patch[p]  = max(0.0, cnveg_cs.leafc_storage_patch[p] - cnveg_cs.leafc_patch[p]) * bgtr_val
            cnveg_cf.frootc_storage_to_xfer_patch[p] = max(0.0, cnveg_cs.frootc_storage_patch[p] - cnveg_cs.frootc_patch[p]) * bgtr_val
            if pftcon.woody[ivt] == 1.0
                cnveg_cf.livestemc_storage_to_xfer_patch[p]  = cnveg_cs.livestemc_storage_patch[p] * bgtr_val
                cnveg_cf.deadstemc_storage_to_xfer_patch[p]  = cnveg_cs.deadstemc_storage_patch[p] * bgtr_val
                cnveg_cf.livecrootc_storage_to_xfer_patch[p] = cnveg_cs.livecrootc_storage_patch[p] * bgtr_val
                cnveg_cf.deadcrootc_storage_to_xfer_patch[p] = cnveg_cs.deadcrootc_storage_patch[p] * bgtr_val
                cnveg_cf.gresp_storage_to_xfer_patch[p]      = cnveg_cs.gresp_storage_patch[p] * bgtr_val
            end
            cnveg_nf.leafn_storage_to_xfer_patch[p]  = cnveg_ns.leafn_storage_patch[p] * bgtr_val
            cnveg_nf.frootn_storage_to_xfer_patch[p] = cnveg_ns.frootn_storage_patch[p] * bgtr_val
            if pftcon.woody[ivt] == 1.0
                cnveg_nf.livestemn_storage_to_xfer_patch[p]  = cnveg_ns.livestemn_storage_patch[p] * bgtr_val
                cnveg_nf.deadstemn_storage_to_xfer_patch[p]  = cnveg_ns.deadstemn_storage_patch[p] * bgtr_val
                cnveg_nf.livecrootn_storage_to_xfer_patch[p] = cnveg_ns.livecrootn_storage_patch[p] * bgtr_val
                cnveg_nf.deadcrootn_storage_to_xfer_patch[p] = cnveg_ns.deadcrootn_storage_patch[p] * bgtr_val
            end
        end

    end # patch loop

    return nothing
end

# ==========================================================================
# get_swindow — Determine next sowing window
# ==========================================================================
function get_swindow(jday::Int, rx_starts::AbstractVector{Int},
                     rx_ends::AbstractVector{Int},
                     param_start::Int, param_end::Int)
    # No prescribed sowing windows
    if maximum(rx_starts) < 1
        return (w=1, start_w=param_start, end_w=param_end)
    end
    # Today is after latest sowing window end
    if jday > maximum(rx_ends)
        return (w=1, start_w=rx_starts[1], end_w=rx_ends[1])
    end
    # Find first window whose end >= today
    for w in eachindex(rx_starts)
        if min(rx_starts[w], rx_ends[w]) < 1
            break
        end
        if jday <= rx_ends[w]
            return (w=w, start_w=rx_starts[w], end_w=rx_ends[w])
        end
    end
    error("get_swindow(): No sowing window found")
end

# ==========================================================================
# was_sown_in_this_window
# ==========================================================================
function was_sown_in_this_window(sowing_window_startdate::Int,
                                 sowing_window_enddate::Int,
                                 jday::Int, idop::Int,
                                 sown_in_this_window::Bool)
    result = sown_in_this_window

    # Check if in sowing window
    is_in_sw = _is_doy_in_interval(sowing_window_startdate, sowing_window_enddate, jday)
    if !is_in_sw
        return false
    end
    # Check if planting date is in the window
    idop_in_sw = _is_doy_in_interval(sowing_window_startdate, sowing_window_enddate, idop)
    if is_in_sw && !idop_in_sw
        return false
    end
    # Check for same window vs different occurrence
    if sowing_window_startdate < sowing_window_enddate && idop > jday
        result = false
    elseif sowing_window_startdate > sowing_window_enddate
        if jday <= sowing_window_enddate && idop <= sowing_window_enddate && idop > jday
            result = false
        elseif jday >= sowing_window_startdate && (idop > jday || idop <= sowing_window_enddate)
            result = false
        end
    end
    return result
end

# Helper: is day-of-year in interval (wrapping around year boundary)
function _is_doy_in_interval(start_doy::Int, end_doy::Int, doy::Int)
    if start_doy <= end_doy
        return doy >= start_doy && doy <= end_doy
    else
        return doy >= start_doy || doy <= end_doy
    end
end

# ==========================================================================
# crop_phenology_init! — Crop-specific initialization
# ==========================================================================
function crop_phenology_init!(pstate::PhenologyState,
                              pftcon::PftConPhenology,
                              patch_data::PatchData,
                              gridcell::GridcellData,
                              npcropmin::Int, npcropmax::Int, maxveg::Int)
    np = length(patch_data.itype)
    pstate.inhemi = zeros(Int, np)

    pstate.minplantjday = fill(typemax(Int), maxveg + 1, inSH)  # 0:maxveg mapped to 1:maxveg+1
    pstate.maxplantjday = fill(typemax(Int), maxveg + 1, inSH)

    pstate.jdayyrstart = [1, 182]

    # Convert planting dates
    for n in npcropmin:npcropmax
        if n <= length(pftcon.is_pft_known_to_model) && pftcon.is_pft_known_to_model[n]
            pstate.minplantjday[n, inNH] = round(Int, pftcon.mnNHplantdate[n])
            pstate.maxplantjday[n, inNH] = round(Int, pftcon.mxNHplantdate[n])
            pstate.minplantjday[n, inSH] = round(Int, pftcon.mnSHplantdate[n])
            pstate.maxplantjday[n, inSH] = round(Int, pftcon.mxSHplantdate[n])
        end
    end

    # Determine hemisphere for each patch
    for p in 1:np
        g = patch_data.gridcell[p]
        if gridcell.latdeg[g] > 0.0
            pstate.inhemi[p] = inNH
        else
            pstate.inhemi[p] = inSH
        end
    end

    # Vernalization constants
    pstate.p1d   = 0.004
    pstate.p1v   = 0.003
    pstate.hti   = 1.0
    pstate.tbase = 0.0

    return nothing
end

# ==========================================================================
# crop_phenology! — Crop lifecycle management (simplified)
# ==========================================================================
function crop_phenology!(pstate::PhenologyState, params::PhenologyParams,
                         mask_pcropp::BitVector,
                         pftcon::PftConPhenology,
                         water_diag::WaterDiagnosticBulkData,
                         temperature::TemperatureData,
                         crop::CropData,
                         canopy_state::CanopyStateData,
                         cnveg_state::CNVegStateData,
                         cnveg_cs::CNVegCarbonStateData,
                         cnveg_ns::CNVegNitrogenStateData,
                         cnveg_cf::CNVegCarbonFluxData,
                         cnveg_nf::CNVegNitrogenFluxData,
                         patch_data::PatchData,
                         gridcell::GridcellData;
                         varctl::VarCtl=VarCtl(),
                         avg_dayspyr::Float64=365.0,
                         jday::Int=1, kyr::Int=1, dayspyr::Float64=365.0,
                         use_fertilizer::Bool=false)
    dt = pstate.dt

    for p in eachindex(mask_pcropp)
        mask_pcropp[p] || continue

        ivt = patch_data.itype[p]
        c   = patch_data.column[p]
        g   = patch_data.gridcell[p]
        h   = pstate.inhemi[p]

        # reset background rates
        cnveg_state.bglfr_patch[p] = 0.0
        cnveg_state.bgtr_patch[p]  = 0.0
        cnveg_state.lgsf_patch[p]  = 0.0

        cnveg_state.onset_flag_patch[p]  = 0.0
        cnveg_state.offset_flag_patch[p] = 0.0

        if crop.croplive_patch[p]
            crop.cphase_patch[p] = cphase_planted
            cnveg_state.onset_counter_patch[p] -= dt

            # grain fill phase: background litterfall
            if crop.hui_patch[p] >= cnveg_state.huigrain_patch[p]
                crop.cphase_patch[p] = cphase_grainfill
                cnveg_state.bglfr_patch[p] = 1.0 / (pftcon.leaf_long[ivt] * avg_dayspyr * SECSPDAY)
            end
        end

    end # patch loop

    return nothing
end

# ==========================================================================
# crop_phase! — Get current crop phase
# ==========================================================================
function crop_phase!(mask_pcropp::BitVector, crop::CropData,
                     cnveg_state::CNVegStateData,
                     crop_phase_out::Vector{Float64})
    for p in eachindex(mask_pcropp)
        mask_pcropp[p] || continue

        if crop.croplive_patch[p]
            crop_phase_out[p] = cphase_planted
            if crop.gddtsoi_patch[p] >= cnveg_state.huileaf_patch[p] &&
               crop.hui_patch[p] < cnveg_state.huigrain_patch[p]
                crop_phase_out[p] = cphase_leafemerge
            elseif crop.hui_patch[p] >= cnveg_state.huigrain_patch[p]
                crop_phase_out[p] = cphase_grainfill
            end
        end
    end
    return nothing
end

# ==========================================================================
# plant_crop! — Initialize crop at planting
# ==========================================================================
function plant_crop!(p::Int, leafcn_in::Float64, jday::Int, kyr::Int,
                     crop::CropData, cnveg_state::CNVegStateData,
                     cnveg_cs::CNVegCarbonStateData,
                     cnveg_ns::CNVegNitrogenStateData,
                     cnveg_cf::CNVegCarbonFluxData,
                     cnveg_nf::CNVegNitrogenFluxData,
                     dt::Float64)
    crop.croplive_patch[p]      = true
    crop.sown_in_this_window[p] = true
    cnveg_state.idop_patch[p]   = jday
    cnveg_state.iyop_patch[p]   = kyr
    crop.harvdate_patch[p]      = NOT_Harvested
    crop.sowing_count[p]       += 1

    seed = _initial_seed_at_planting[]
    cnveg_cs.leafc_xfer_patch[p] = seed
    cnveg_ns.leafn_xfer_patch[p] = seed / leafcn_in
    cnveg_cf.crop_seedc_to_leaf_patch[p] = seed / dt
    cnveg_nf.crop_seedn_to_leaf_patch[p] = (seed / leafcn_in) / dt

    # Initialize allocation coefficients
    cnveg_state.aleaf_patch[p]  = 1.0
    cnveg_state.aleafi_patch[p] = 1.0
    cnveg_state.astem_patch[p]  = 0.0
    cnveg_state.astemi_patch[p] = 0.0
    cnveg_state.aroot_patch[p]  = 0.0

    return nothing
end

# ==========================================================================
# days_past_planting — Calculate days since planting
# ==========================================================================
function days_past_planting(idop::Int, jday::Int; dayspyr::Int=365)
    if jday >= idop
        return jday - idop
    else
        return jday - idop + dayspyr
    end
end

# ==========================================================================
# vernalization! — Vernalization for winter cereal
# ==========================================================================
function vernalization!(p::Int, pstate::PhenologyState,
                        canopy_state::CanopyStateData,
                        temperature::TemperatureData,
                        water_diag::WaterDiagnosticBulkData,
                        cnveg_state::CNVegStateData,
                        crop::CropData,
                        patch_data::PatchData)
    c = patch_data.column[p]
    p1v   = pstate.p1v
    hti   = pstate.hti
    tbase = pstate.tbase

    t_ref2m     = temperature.t_ref2m_patch[p]
    t_ref2m_min = temperature.t_ref2m_min_patch[p]
    t_ref2m_max = temperature.t_ref2m_max_patch[p]
    snow_depth  = water_diag.snow_depth_col[c]

    force_harvest = false

    # crown temperature
    if t_ref2m < TFRZ
        tcrown = 2.0 + (t_ref2m - TFRZ) * (0.4 + 0.0018 *
            (min(snow_depth * 100.0, 15.0) - 15.0)^2)
    else
        tcrown = t_ref2m - TFRZ
    end

    # vernalization factor
    if t_ref2m_max > TFRZ
        if t_ref2m_min <= TFRZ + 15.0
            vd1 = 1.4 - 0.0778 * tcrown
            vd2 = 0.5 + 13.44 / ((t_ref2m_max - t_ref2m_min + 3.0)^2) * tcrown
            vd  = max(0.0, min(1.0, vd1, vd2))
            cnveg_state.cumvd_patch[p] += vd
        end
        if cnveg_state.cumvd_patch[p] < 10.0 && t_ref2m_max > TFRZ + 30.0
            cnveg_state.cumvd_patch[p] -= 0.5 * (t_ref2m_max - TFRZ - 30.0)
        end
        cnveg_state.cumvd_patch[p] = max(0.0, cnveg_state.cumvd_patch[p])

        crop.vf_patch[p] = 1.0 - p1v * (50.0 - cnveg_state.cumvd_patch[p])
        crop.vf_patch[p] = max(0.0, min(crop.vf_patch[p], 1.0))
    end

    # cold hardening
    hdidx = cnveg_state.hdidx_patch[p]
    if t_ref2m_min <= TFRZ - 3.0 || hdidx != 0.0
        if hdidx >= hti
            hdidx += 0.083
            hdidx = min(hdidx, hti * 2.0)
        end
        if t_ref2m_max >= tbase + TFRZ + 10.0
            hdidx -= 0.02 * (t_ref2m_max - tbase - TFRZ - 10.0)
            if hdidx > hti
                hdidx -= 0.02 * (t_ref2m_max - tbase - TFRZ - 10.0)
            end
            hdidx = max(0.0, hdidx)
        end
    elseif tcrown >= tbase - 1.0
        if tcrown <= tbase + 8.0
            hdidx += 0.1 - (tcrown - tbase + 3.5)^2 / 506.0
            if hdidx >= hti && tcrown <= tbase + 0.0
                hdidx += 0.083
                hdidx = min(hdidx, hti * 2.0)
            end
        end
        if t_ref2m_max >= tbase + TFRZ + 10.0
            hdidx -= 0.02 * (t_ref2m_max - tbase - TFRZ - 10.0)
            if hdidx > hti
                hdidx -= 0.02 * (t_ref2m_max - tbase - TFRZ - 10.0)
            end
            hdidx = max(0.0, hdidx)
        end
    end
    cnveg_state.hdidx_patch[p] = hdidx

    # killing temperature check
    if t_ref2m_min <= TFRZ - 6.0
        tkil = (tbase - 6.0) - 6.0 * hdidx
        if tkil >= tcrown
            if (0.95 - 0.02 * (tcrown - tkil)^2) < 0.02 && canopy_state.tlai_patch[p] > 0.0
                force_harvest = true
            end
        end
    end

    return force_harvest
end

# ==========================================================================
# cn_onset_growth! — Transfer → display during onset
# ==========================================================================
function cn_onset_growth!(pstate::PhenologyState,
                          mask_soilp::BitVector,
                          pftcon::PftConPhenology,
                          cnveg_state::CNVegStateData,
                          cnveg_cs::CNVegCarbonStateData,
                          cnveg_ns::CNVegNitrogenStateData,
                          cnveg_cf::CNVegCarbonFluxData,
                          cnveg_nf::CNVegNitrogenFluxData,
                          patch_data::PatchData)
    dt = pstate.dt

    for p in eachindex(mask_soilp)
        mask_soilp[p] || continue

        ivt = patch_data.itype[p]

        # onset period transfer
        if cnveg_state.onset_flag_patch[p] == 1.0
            if abs(cnveg_state.onset_counter_patch[p] - dt) <= dt / 2.0
                t1 = 1.0 / dt
            else
                t1 = 2.0 / cnveg_state.onset_counter_patch[p]
            end

            cnveg_cf.leafc_xfer_to_leafc_patch[p]   = t1 * cnveg_cs.leafc_xfer_patch[p]
            cnveg_cf.frootc_xfer_to_frootc_patch[p] = t1 * cnveg_cs.frootc_xfer_patch[p]
            cnveg_nf.leafn_xfer_to_leafn_patch[p]   = t1 * cnveg_ns.leafn_xfer_patch[p]
            cnveg_nf.frootn_xfer_to_frootn_patch[p] = t1 * cnveg_ns.frootn_xfer_patch[p]

            if pftcon.woody[ivt] == 1.0
                cnveg_cf.livestemc_xfer_to_livestemc_patch[p]   = t1 * cnveg_cs.livestemc_xfer_patch[p]
                cnveg_cf.deadstemc_xfer_to_deadstemc_patch[p]   = t1 * cnveg_cs.deadstemc_xfer_patch[p]
                cnveg_cf.livecrootc_xfer_to_livecrootc_patch[p] = t1 * cnveg_cs.livecrootc_xfer_patch[p]
                cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch[p] = t1 * cnveg_cs.deadcrootc_xfer_patch[p]
                cnveg_nf.livestemn_xfer_to_livestemn_patch[p]   = t1 * cnveg_ns.livestemn_xfer_patch[p]
                cnveg_nf.deadstemn_xfer_to_deadstemn_patch[p]   = t1 * cnveg_ns.deadstemn_xfer_patch[p]
                cnveg_nf.livecrootn_xfer_to_livecrootn_patch[p] = t1 * cnveg_ns.livecrootn_xfer_patch[p]
                cnveg_nf.deadcrootn_xfer_to_deadcrootn_patch[p] = t1 * cnveg_ns.deadcrootn_xfer_patch[p]
            end
        end

        # background transfer growth
        if cnveg_state.bgtr_patch[p] > 0.0
            cnveg_cf.leafc_xfer_to_leafc_patch[p]   = cnveg_cs.leafc_xfer_patch[p] / dt
            cnveg_cf.frootc_xfer_to_frootc_patch[p] = cnveg_cs.frootc_xfer_patch[p] / dt
            cnveg_nf.leafn_xfer_to_leafn_patch[p]   = cnveg_ns.leafn_xfer_patch[p] / dt
            cnveg_nf.frootn_xfer_to_frootn_patch[p] = cnveg_ns.frootn_xfer_patch[p] / dt

            if pftcon.woody[ivt] == 1.0
                cnveg_cf.livestemc_xfer_to_livestemc_patch[p]   = cnveg_cs.livestemc_xfer_patch[p] / dt
                cnveg_cf.deadstemc_xfer_to_deadstemc_patch[p]   = cnveg_cs.deadstemc_xfer_patch[p] / dt
                cnveg_cf.livecrootc_xfer_to_livecrootc_patch[p] = cnveg_cs.livecrootc_xfer_patch[p] / dt
                cnveg_cf.deadcrootc_xfer_to_deadcrootc_patch[p] = cnveg_cs.deadcrootc_xfer_patch[p] / dt
                cnveg_nf.livestemn_xfer_to_livestemn_patch[p]   = cnveg_ns.livestemn_xfer_patch[p] / dt
                cnveg_nf.deadstemn_xfer_to_deadstemn_patch[p]   = cnveg_ns.deadstemn_xfer_patch[p] / dt
                cnveg_nf.livecrootn_xfer_to_livecrootn_patch[p] = cnveg_ns.livecrootn_xfer_patch[p] / dt
                cnveg_nf.deadcrootn_xfer_to_deadcrootn_patch[p] = cnveg_ns.deadcrootn_xfer_patch[p] / dt
            end
        end

    end # patch loop

    return nothing
end

# ==========================================================================
# cn_offset_litterfall! — Display → litter during offset
# ==========================================================================
function cn_offset_litterfall!(pstate::PhenologyState,
                               mask_soilp::BitVector,
                               pftcon::PftConPhenology,
                               cnveg_state::CNVegStateData,
                               cnveg_cs::CNVegCarbonStateData,
                               cnveg_ns::CNVegNitrogenStateData,
                               cnveg_cf::CNVegCarbonFluxData,
                               cnveg_nf::CNVegNitrogenFluxData,
                               crop::CropData,
                               patch_data::PatchData;
                               use_fun::Bool=false,
                               CNratio_floating::Bool=false,
                               for_testing_no_crop_seed_replenishment::Bool=false)
    dt = pstate.dt

    for p in eachindex(mask_soilp)
        mask_soilp[p] || continue

        ivt = patch_data.itype[p]

        if cnveg_state.offset_flag_patch[p] != 1.0
            continue
        end

        if abs(cnveg_state.offset_counter_patch[p] - dt) <= dt / 2.0
            # last timestep of offset: dump all remaining
            t1 = 1.0 / dt
            cnveg_cf.frootc_to_litter_patch[p] = t1 * cnveg_cs.frootc_patch[p]
            cnveg_cf.leafc_to_litter_patch[p]  = t1 * cnveg_cs.leafc_patch[p]
        else
            t1 = dt * 2.0 / (cnveg_state.offset_counter_patch[p]^2)
            cnveg_cf.leafc_to_litter_patch[p]  = cnveg_cf.prev_leafc_to_litter_patch[p] +
                t1 * (cnveg_cs.leafc_patch[p] - cnveg_cf.prev_leafc_to_litter_patch[p] *
                      cnveg_state.offset_counter_patch[p])
            cnveg_cf.frootc_to_litter_patch[p] = cnveg_cf.prev_frootc_to_litter_patch[p] +
                t1 * (cnveg_cs.frootc_patch[p] - cnveg_cf.prev_frootc_to_litter_patch[p] *
                      cnveg_state.offset_counter_patch[p])
        end

        # Nitrogen litterfall
        if CNratio_floating
            fr_leafn_to_litter = 0.5
            if cnveg_cs.leafc_patch[p] == 0.0
                ntovr_leaf = 0.0
            else
                ntovr_leaf = cnveg_cf.leafc_to_litter_patch[p] *
                    (cnveg_ns.leafn_patch[p] / cnveg_cs.leafc_patch[p])
            end
            cnveg_nf.leafn_to_litter_patch[p]   = fr_leafn_to_litter * ntovr_leaf
            cnveg_nf.leafn_to_retransn_patch[p] = ntovr_leaf - cnveg_nf.leafn_to_litter_patch[p]
            if cnveg_cs.frootc_patch[p] == 0.0
                cnveg_nf.frootn_to_litter_patch[p] = 0.0
            else
                cnveg_nf.frootn_to_litter_patch[p] = cnveg_cf.frootc_to_litter_patch[p] *
                    (cnveg_ns.frootn_patch[p] / cnveg_cs.frootc_patch[p])
            end
        else
            cnveg_nf.leafn_to_litter_patch[p]   = cnveg_cf.leafc_to_litter_patch[p] / pftcon.lflitcn[ivt]
            cnveg_nf.leafn_to_retransn_patch[p] = (cnveg_cf.leafc_to_litter_patch[p] / pftcon.leafcn[ivt]) -
                cnveg_nf.leafn_to_litter_patch[p]
            cnveg_nf.frootn_to_litter_patch[p]  = cnveg_cf.frootc_to_litter_patch[p] / pftcon.frootcn[ivt]
        end

        # save current litterfall fluxes
        cnveg_cf.prev_leafc_to_litter_patch[p]  = cnveg_cf.leafc_to_litter_patch[p]
        cnveg_cf.prev_frootc_to_litter_patch[p] = cnveg_cf.frootc_to_litter_patch[p]

    end # patch loop

    return nothing
end

# ==========================================================================
# cn_background_litterfall! — Background litterfall
# ==========================================================================
function cn_background_litterfall!(pstate::PhenologyState,
                                   mask_soilp::BitVector,
                                   pftcon::PftConPhenology,
                                   cnveg_state::CNVegStateData,
                                   cnveg_cs::CNVegCarbonStateData,
                                   cnveg_ns::CNVegNitrogenStateData,
                                   cnveg_cf::CNVegCarbonFluxData,
                                   cnveg_nf::CNVegNitrogenFluxData,
                                   patch_data::PatchData;
                                   use_fun::Bool=false,
                                   CNratio_floating::Bool=false)
    for p in eachindex(mask_soilp)
        mask_soilp[p] || continue

        ivt = patch_data.itype[p]

        if cnveg_state.bglfr_patch[p] > 0.0
            cnveg_cf.leafc_to_litter_patch[p]  = cnveg_state.bglfr_patch[p] * cnveg_cs.leafc_patch[p]
            cnveg_cf.frootc_to_litter_patch[p] = cnveg_state.bglfr_patch[p] * cnveg_cs.frootc_patch[p]

            if CNratio_floating
                fr_leafn_to_litter = 0.5
                if cnveg_cs.leafc_patch[p] == 0.0
                    ntovr_leaf = 0.0
                else
                    ntovr_leaf = cnveg_cf.leafc_to_litter_patch[p] *
                        (cnveg_ns.leafn_patch[p] / cnveg_cs.leafc_patch[p])
                end
                cnveg_nf.leafn_to_litter_patch[p]   = fr_leafn_to_litter * ntovr_leaf
                cnveg_nf.leafn_to_retransn_patch[p] = ntovr_leaf - cnveg_nf.leafn_to_litter_patch[p]
                if cnveg_cs.frootc_patch[p] == 0.0
                    cnveg_nf.frootn_to_litter_patch[p] = 0.0
                else
                    cnveg_nf.frootn_to_litter_patch[p] = cnveg_cf.frootc_to_litter_patch[p] *
                        (cnveg_ns.frootn_patch[p] / cnveg_cs.frootc_patch[p])
                end
            else
                cnveg_nf.leafn_to_litter_patch[p]   = cnveg_cf.leafc_to_litter_patch[p] / pftcon.lflitcn[ivt]
                cnveg_nf.leafn_to_retransn_patch[p] = (cnveg_cf.leafc_to_litter_patch[p] / pftcon.leafcn[ivt]) -
                    cnveg_nf.leafn_to_litter_patch[p]
                cnveg_nf.frootn_to_litter_patch[p]  = cnveg_cf.frootc_to_litter_patch[p] / pftcon.frootcn[ivt]
            end
        end

    end # patch loop

    return nothing
end

# ==========================================================================
# cn_livewood_turnover! — Live wood → dead wood turnover
# ==========================================================================
function cn_livewood_turnover!(pstate::PhenologyState,
                               mask_soilp::BitVector,
                               pftcon::PftConPhenology,
                               cnveg_cs::CNVegCarbonStateData,
                               cnveg_ns::CNVegNitrogenStateData,
                               cnveg_cf::CNVegCarbonFluxData,
                               cnveg_nf::CNVegNitrogenFluxData,
                               patch_data::PatchData;
                               CNratio_floating::Bool=false)
    lwtop = pstate.lwtop

    for p in eachindex(mask_soilp)
        mask_soilp[p] || continue

        ivt = patch_data.itype[p]

        if pftcon.woody[ivt] <= 0.0
            continue
        end

        # live stem → dead stem
        ctovr = cnveg_cs.livestemc_patch[p] * lwtop
        ntovr = ctovr / pftcon.livewdcn[ivt]
        cnveg_cf.livestemc_to_deadstemc_patch[p] = ctovr
        cnveg_nf.livestemn_to_deadstemn_patch[p] = ctovr / pftcon.deadwdcn[ivt]

        if CNratio_floating
            if cnveg_cs.livestemc_patch[p] == 0.0
                ntovr = 0.0
                cnveg_nf.livestemn_to_deadstemn_patch[p] = 0.0
            else
                ntovr = ctovr * (cnveg_ns.livestemn_patch[p] / cnveg_cs.livestemc_patch[p])
                cnveg_nf.livestemn_to_deadstemn_patch[p] = ctovr / pftcon.deadwdcn[ivt]
            end
        end
        cnveg_nf.livestemn_to_retransn_patch[p] = ntovr - cnveg_nf.livestemn_to_deadstemn_patch[p]

        # live coarse root → dead coarse root
        ctovr = cnveg_cs.livecrootc_patch[p] * lwtop
        ntovr = ctovr / pftcon.livewdcn[ivt]
        cnveg_cf.livecrootc_to_deadcrootc_patch[p] = ctovr
        cnveg_nf.livecrootn_to_deadcrootn_patch[p] = ctovr / pftcon.deadwdcn[ivt]

        if CNratio_floating
            if cnveg_cs.livecrootc_patch[p] == 0.0
                ntovr = 0.0
                cnveg_nf.livecrootn_to_deadcrootn_patch[p] = 0.0
            else
                ntovr = ctovr * (cnveg_ns.livecrootn_patch[p] / cnveg_cs.livecrootc_patch[p])
                cnveg_nf.livecrootn_to_deadcrootn_patch[p] = ctovr / pftcon.deadwdcn[ivt]
            end
        end
        cnveg_nf.livecrootn_to_retransn_patch[p] = ntovr - cnveg_nf.livecrootn_to_deadcrootn_patch[p]

    end # patch loop

    return nothing
end

# ==========================================================================
# cn_crop_harvest_to_product_pools! — Crop harvest to product pools
# ==========================================================================
function cn_crop_harvest_to_product_pools!(mask_soilp::BitVector,
                                           mask_soilc::BitVector,
                                           cnveg_cf::CNVegCarbonFluxData,
                                           cnveg_nf::CNVegNitrogenFluxData,
                                           patch_data::PatchData;
                                           use_crop::Bool=false)
    if !use_crop
        return nothing
    end

    for p in eachindex(mask_soilp)
        mask_soilp[p] || continue

        cnveg_cf.crop_harvestc_to_cropprodc_patch[p] =
            cnveg_cf.leafc_to_biofuelc_patch[p] +
            cnveg_cf.livestemc_to_biofuelc_patch[p] +
            cnveg_cf.leafc_to_removedresiduec_patch[p] +
            cnveg_cf.livestemc_to_removedresiduec_patch[p]

        cnveg_nf.crop_harvestn_to_cropprodn_patch[p] =
            cnveg_nf.leafn_to_biofueln_patch[p] +
            cnveg_nf.livestemn_to_biofueln_patch[p] +
            cnveg_nf.leafn_to_removedresiduen_patch[p] +
            cnveg_nf.livestemn_to_removedresiduen_patch[p]
    end

    return nothing
end

# ==========================================================================
# cn_litter_to_column! — Aggregate patch litter to column level
# ==========================================================================
function cn_litter_to_column!(mask_soilp::BitVector,
                              pftcon::PftConPhenology,
                              cnveg_state::CNVegStateData,
                              cnveg_cf::CNVegCarbonFluxData,
                              cnveg_nf::CNVegNitrogenFluxData,
                              patch_data::PatchData,
                              leaf_prof::Matrix{Float64},
                              froot_prof::Matrix{Float64};
                              nlevdecomp::Int=1,
                              i_litr_min::Int=1, i_litr_max::Int=3,
                              npcropmin::Int=17,
                              use_grainproduct::Bool=false)
    for j in 1:nlevdecomp
        for p in eachindex(mask_soilp)
            mask_soilp[p] || continue

            ivt = patch_data.itype[p]
            c   = patch_data.column[p]
            wt  = patch_data.wtcol[p]

            for i in i_litr_min:i_litr_max
                # leaf litter C and N
                if size(cnveg_cf.phenology_c_to_litr_c_col, 1) >= c
                    cnveg_cf.phenology_c_to_litr_c_col[c, j, i] +=
                        cnveg_cf.leafc_to_litter_patch[p] * pftcon.lf_f[ivt, i] * wt * leaf_prof[p, j]

                    cnveg_nf.phenology_n_to_litr_n_col[c, j, i] +=
                        cnveg_nf.leafn_to_litter_patch[p] * pftcon.lf_f[ivt, i] * wt * leaf_prof[p, j]

                    # fine root litter C and N
                    cnveg_cf.phenology_c_to_litr_c_col[c, j, i] +=
                        cnveg_cf.frootc_to_litter_patch[p] * pftcon.fr_f[ivt, i] * wt * froot_prof[p, j]

                    cnveg_nf.phenology_n_to_litr_n_col[c, j, i] +=
                        cnveg_nf.frootn_to_litter_patch[p] * pftcon.fr_f[ivt, i] * wt * froot_prof[p, j]
                end
            end

            # crop stem litter uses leaf litter fractions
            if ivt >= npcropmin
                for i in i_litr_min:i_litr_max
                    if size(cnveg_cf.phenology_c_to_litr_c_col, 1) >= c
                        cnveg_cf.phenology_c_to_litr_c_col[c, j, i] +=
                            cnveg_cf.livestemc_to_litter_patch[p] * pftcon.lf_f[ivt, i] * wt * leaf_prof[p, j]
                        cnveg_nf.phenology_n_to_litr_n_col[c, j, i] +=
                            cnveg_nf.livestemn_to_litter_patch[p] * pftcon.lf_f[ivt, i] * wt * leaf_prof[p, j]
                    end
                end
            end

        end # patch loop
    end # decomp level loop

    return nothing
end
