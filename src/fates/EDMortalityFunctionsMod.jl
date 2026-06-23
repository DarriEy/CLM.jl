# EDMortalityFunctionsMod.jl
# Julia port of FATES src/fates/biogeochem/EDMortalityFunctionsMod.F90
#
# Functions that control per-cohort mortality. This module computes the
# individual (non-disturbance) mortality rates that act on a cohort's number
# density:
#   * carbon-starvation mortality (cmort)  -- storage-pool driven
#   * hydraulic-failure mortality (hmort)  -- btran / plant-hydro driven
#   * background mortality (bmort)         -- constant PFT rate
#   * freezing / cold-stress mortality (frmort)
#   * size-dependent senescence (smort)
#   * age-dependent senescence (asmort)
#   * damage-dependent mortality (dgmort)
# plus the `mortality_rates` aggregator and `Mortality_Derivative`, which turns
# the rates (plus understory logging) into d(n)/dt for the cohort-number ODE, and
# the small `ExemptTreefallDist` predicate.
#
# Original authors: created 10/30/09 Rosie Fisher; refactored 02/20/18 Ryan Knox.
#
# Deps (already ported):
#   FatesConstantsMod  : fates_r8->Float64, itrue, ifalse, nearzero,
#                        cstarvation_model_lin, cstarvation_model_exp,
#                        ihard_stress_decid, isemi_stress_decid, leaves_off,
#                        t_water_freeze_k_1atm (tfrz), fates_check_param_set
#   FatesGlobals       : fates_endrun, fates_log
#   EDPftvarcon        : edpftvarcon_inst (mort_* / hf_* / freezetol / bmort /
#                        prescribed_mortality_* parameter vectors)
#   FatesCohortMod     : fates_cohort_type
#   EDTypesMod         : ed_site_type
#   EDParamsMod        : maxpft, mort_cstarvation_model (via ed_params()),
#                        soil_tfrz_thresh, fates_mortality_disturbance_fraction
#                        (via ed_params()), ed_params()
#   FatesAllometryMod  : bleaf, storage_fraction_of_target
#   FatesInterfaceTypesMod : bc_in_type, hlm_use_ed_prescribed_phys, hlm_freq_day,
#                        hlm_use_planthydro, hlm_use_tree_damage
#   EDLoggingMortalityMod : LoggingMortality_frac
#   PRTGenericMod      : carbon12_element, store_organ, GetState
#   PRTParametersMod   : prt_params (stress_decid, season_decid, woody)
#   DamageMainMod      : GetDamageMortality
#
# NOT-yet-ported / gated-off dependency (stubbed below):
#   * Plant-hydraulics hmort path reads cohort_in%co_hydr%ftc_* . The cohort
#     hydraulics arrays exist (ed_cohort_hydr_type) but the ftc_* field names are
#     part of the FatesPlantHydraulicsMod state that is only populated when
#     hlm_use_planthydro == itrue (off in the default + test path). The hydro
#     branch is ported faithfully but guarded; see the TODO in mortality_rates.
#
# Standalone -- NOT added to CLMInstances or any dual-copied struct.

# ===========================================================================
# mortality_rates
# ===========================================================================
"""
    mortality_rates(cohort_in, bc_in, btran_ft, mean_temp)
        -> (cmort, hmort, bmort, frmort, smort, asmort, dgmort)

Calculate the per-cohort mortality rates [fraction per year] from carbon storage,
hydraulic cavitation, background, freezing, and size/age-dependent senescence (and
damage). Returns the seven rates as a tuple (the Fortran `intent(out)` arguments).
Mirrors the Fortran `mortality_rates`, preserving variable names and the
prescribed-physiology / plant-hydro / carbon-starvation-model branches.

`btran_ft` is the per-PFT soil-moisture transpiration stress factor (length
`maxpft`); `mean_temp` is the daily-mean temperature [K].
"""
function mortality_rates(cohort_in::fates_cohort_type, bc_in::bc_in_type,
                         btran_ft::AbstractVector{<:Real}, mean_temp::Real)

    tfrz = t_water_freeze_k_1atm
    p    = edpftvarcon_inst()
    edp  = ed_params()
    ipft = cohort_in.pft

    # 5deg buffer for freezing mortality.
    frost_mort_buffer = 5.0
    # Developer test which may help debug carbon imbalances (off, as in Fortran).
    test_zero_mortality = false

    # Outputs (defaults match the Fortran intent(out) flow).
    bmort  = 0.0
    cmort  = 0.0
    hmort  = 0.0
    frmort = 0.0
    smort  = 0.0
    asmort = 0.0
    dgmort = 0.0

    # Check if the PFT is deciduous and leaves are completely abscised. If so we
    # prevent hydraulic-failure mortality, as leafless plants cannot die of it.
    # Both drought-deciduous and cold-deciduous are considered.
    is_decid_dormant =
        ((prt_params.stress_decid[ipft] == ihard_stress_decid ||   # drought deciduous
          prt_params.stress_decid[ipft] == isemi_stress_decid ||   # semi-deciduous
          prt_params.season_decid[ipft] == itrue)             &&   # cold deciduous
         (cohort_in.status_coh == leaves_off))                     # fully abscised

    # Size Dependent Senescence: rate (r) and inflection point (ip) define the
    # increase in mortality rate with dbh.
    mort_r_size_senescence  = p.mort_r_size_senescence[ipft]
    mort_ip_size_senescence = p.mort_ip_size_senescence[ipft]

    if mort_ip_size_senescence < fates_check_param_set
        smort = 1.0 / (1.0 + exp(-1.0 * mort_r_size_senescence *
                                 (cohort_in.dbh - mort_ip_size_senescence)))
    else
        smort = 0.0
    end

    # Age Dependent Senescence.
    mort_r_age_senescence  = p.mort_r_age_senescence[ipft]
    mort_ip_age_senescence = p.mort_ip_age_senescence[ipft]

    if mort_ip_age_senescence < fates_check_param_set
        asmort = 1.0 / (1.0 + exp(-1.0 * mort_r_age_senescence *
                                  (cohort_in.coage - mort_ip_age_senescence)))
    else
        asmort = 0.0
    end

    # Damage dependent mortality.
    if hlm_use_tree_damage[] == itrue
        dgmort = GetDamageMortality(cohort_in.crowndamage, cohort_in.pft)
    else
        dgmort = 0.0
    end

    if hlm_use_ed_prescribed_phys[] == ifalse

        # 'Background' mortality (constant PFT rate here).
        bmort = p.bmort[ipft]

        # Proxy for hydraulic-failure-induced mortality.
        hf_sm_threshold  = p.hf_sm_threshold[ipft]
        hf_flc_threshold = p.hf_flc_threshold[ipft]

        if hlm_use_planthydro[] == itrue
            # TODO Batch NN: the plant-hydraulics flc path reads the cohort's
            # ftc_ag / ftc_troot / ftc_aroot fractional-loss-of-conductivity
            # arrays, which are populated by FatesPlantHydraulicsMod (only when
            # hlm_use_planthydro == itrue, off in the default + test path). The
            # branch is ported faithfully; it requires those co_hydr fields to be
            # present on the cohort hydraulics object.
            #   note the flc is set as the fraction of max conductivity in hydro
            min_fmc_ag = minimum(cohort_in.co_hydr.ftc_ag)
            min_fmc_tr = cohort_in.co_hydr.ftc_troot
            min_fmc_ar = minimum(cohort_in.co_hydr.ftc_aroot)
            min_fmc = min(min_fmc_ag, min_fmc_tr)
            min_fmc = min(min_fmc, min_fmc_ar)
            flc = 1.0 - min_fmc
            if flc >= hf_flc_threshold && hf_flc_threshold < 1.0
                hmort = (flc - hf_flc_threshold) / (1.0 - hf_flc_threshold) *
                        p.mort_scalar_hydrfailure[ipft]
            else
                hmort = 0.0
            end
        else
            # When FATES-Hydro is off, hydraulic-failure mortality occurs only when
            # btran falls below a threshold and plants have leaves.
            if (!is_decid_dormant) &&
               (btran_ft[ipft] <= hf_sm_threshold) &&
               ((minimum(bc_in.t_soisno_sl) - tfrz) > soil_tfrz_thresh)
                hmort = p.mort_scalar_hydrfailure[ipft]
            else
                hmort = 0.0
            end
        end

        # Carbon Starvation induced mortality.
        if cohort_in.dbh > 0.0
            # Ratio between storage biomass and (fully-flushed, damage/trim-
            # accounted) target leaf biomass, used to define carbon starvation.
            target_leaf_c, _ = bleaf(cohort_in.dbh, cohort_in.pft,
                                     cohort_in.crowndamage, cohort_in.canopy_trim, 1.0)
            store_c = GetState(cohort_in.prt, store_organ, carbon12_element)
            frac = storage_fraction_of_target(target_leaf_c, store_c)

            # Select the carbon-starvation mortality model (linear or exponential).
            if edp.mort_cstarvation_model == cstarvation_model_lin
                # Linear model. Zero when frac >= mort_upthresh_cstarvation, rising
                # to the max (mort_scalar_cstarvation) when frac = 0.
                cmort = p.mort_scalar_cstarvation[ipft] *
                    max(0.0, (p.mort_upthresh_cstarvation[ipft] - frac) /
                             p.mort_upthresh_cstarvation[ipft])

            elseif edp.mort_cstarvation_model == cstarvation_model_exp
                # Exponential model. Max (mort_scalar_cstarvation) at frac=0;
                # mort_upthresh_cstarvation controls the e-folding decay.
                cmort = p.mort_scalar_cstarvation[ipft] *
                    exp(-frac / p.mort_upthresh_cstarvation[ipft])

            else
                fates_endrun("Invalid carbon starvation model (" *
                             string(edp.mort_cstarvation_model) * ").")
            end

            # Make sure the mortality is set to zero when tiny.
            if cmort <= nearzero
                cmort = 0.0
            end
        else
            fates_endrun("dbh problem in mortality_rates: dbh=$(cohort_in.dbh) " *
                         "pft=$(cohort_in.pft) n=$(cohort_in.n) " *
                         "canopy_layer=$(cohort_in.canopy_layer)")
        end

        # ------------------------------------------------------------------------
        # Mortality due to cold and freezing stress (frmort), based on ED2 and
        # Albani et al. (2006), Glob. Change Biol., 12, 2370-2390.
        # ------------------------------------------------------------------------
        temp_in_C = mean_temp - tfrz

        temp_dep_fraction = max(0.0, min(1.0,
            1.0 - (temp_in_C - p.freezetol[ipft]) / frost_mort_buffer))
        frmort = p.mort_scalar_coldstress[ipft] * temp_dep_fraction

    else  # hlm_use_ed_prescribed_phys is true

        if cohort_in.canopy_layer == 1
            bmort = p.prescribed_mortality_canopy[ipft]
        else
            bmort = p.prescribed_mortality_understory[ipft]
        end
        cmort  = 0.0
        hmort  = 0.0
        frmort = 0.0
    end

    if test_zero_mortality
        cmort  = 0.0
        hmort  = 0.0
        frmort = 0.0
        bmort  = 0.0
        smort  = 0.0
        asmort = 0.0
        dgmort = 0.0
    end

    return cmort, hmort, bmort, frmort, smort, asmort, dgmort
end

# ===========================================================================
# Mortality_Derivative
# ===========================================================================
"""
    Mortality_Derivative(currentSite, currentCohort, bc_in, btran_ft, mean_temp,
                         land_use_label, age_since_anthro_disturbance,
                         frac_site_primary, frac_site_secondary,
                         harvestable_forest_c, harvest_tag)

Calculate the change in number density per unit time (`currentCohort.dndt`) from
the contributing (non-disturbance) mortality rates plus understory logging. These
are NOT disturbance-inducing rates (canopy logging disturbance is handled
elsewhere). Updates `currentCohort.dndt` and the cohort's logging-mortality
fields (`lmort_*`/`l_degrad`) in place, and writes `harvest_tag` in place. Mirrors
the Fortran `Mortality_Derivative`.

`harvest_tag` is the per-harvest-category status vector (modified in place).
"""
function Mortality_Derivative(currentSite::ed_site_type,
                              currentCohort::fates_cohort_type, bc_in::bc_in_type,
                              btran_ft::AbstractVector{<:Real}, mean_temp::Real,
                              land_use_label::Integer,
                              age_since_anthro_disturbance::Real,
                              frac_site_primary::Real, frac_site_secondary::Real,
                              harvestable_forest_c::AbstractVector{<:Real},
                              harvest_tag::AbstractVector{<:Integer})

    ipft = currentCohort.pft

    # Mortality for trees in the understorey. If trees are in the canopy, then
    # their death is 'disturbance' (handled elsewhere).
    cmort, hmort, bmort, frmort, smort, asmort, dgmort =
        mortality_rates(currentCohort, bc_in, btran_ft, mean_temp)
    currentCohort.cmort  = cmort
    currentCohort.hmort  = hmort
    currentCohort.bmort  = bmort
    currentCohort.frmort = frmort
    currentCohort.smort  = smort
    currentCohort.asmort = asmort
    currentCohort.dgmort = dgmort

    lmort_direct, lmort_collateral, lmort_infra, l_degrad, htag =
        LoggingMortality_frac(currentSite, bc_in, ipft, currentCohort.dbh,
            currentCohort.canopy_layer, bc_in.hlm_harvest_rates,
            bc_in.hlm_harvest_catnames, bc_in.hlm_harvest_units,
            land_use_label, age_since_anthro_disturbance,
            frac_site_primary, frac_site_secondary, harvestable_forest_c)

    currentCohort.lmort_direct     = lmort_direct
    currentCohort.lmort_collateral = lmort_collateral
    currentCohort.lmort_infra      = lmort_infra
    currentCohort.l_degrad         = l_degrad
    # Mirror the Fortran intent(out) harvest_tag write-back into the caller's array.
    @views harvest_tag[1:length(htag)] .= htag

    if currentCohort.canopy_layer > 1
        # Include understory logging mortality rates not associated with
        # disturbance.
        dndt_logging = (currentCohort.lmort_direct +
                        currentCohort.lmort_collateral +
                        currentCohort.lmort_infra) / hlm_freq_day[]

        currentCohort.dndt = -1.0 *
            (cmort + hmort + bmort + frmort + smort + asmort + dgmort + dndt_logging) *
            currentCohort.n
    else
        # Mortality from logging in the canopy is ONLY disturbance generating; don't
        # update number densities via non-disturbance death here. For plants whose
        # death is not considered disturbance (i.e. grasses), include all of their
        # mortality here.
        currentCohort.dndt =
            -(cmort + hmort + bmort + frmort + smort + asmort + dgmort) *
            currentCohort.n
        if !ExemptTreefallDist(currentCohort)
            currentCohort.dndt =
                (1.0 - ed_params().fates_mortality_disturbance_fraction) *
                currentCohort.dndt
        end
    end

    return nothing
end

# ===========================================================================
# ExemptTreefallDist
# ===========================================================================
"""
    ExemptTreefallDist(ccohort) -> Bool

Determine whether to consider some fraction of a cohort's crown area as disturbed
patch area when individuals of that cohort die. Current logic only exempts
non-woody plants (grasses): if the cohort's PFT is non-woody, all of its mortality
is treated as non-disturbance-generating (`true`). Mirrors the Fortran
`ExemptTreefallDist`.
"""
function ExemptTreefallDist(ccohort::fates_cohort_type)
    return prt_params.woody[ccohort.pft] == ifalse
end
