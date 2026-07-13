# EDMainMod.jl
# Julia port of FATES src/fates/main/EDMainMod.F90 (Batch 16).
#
# The DAILY FATES ECOSYSTEM-DYNAMICS DRIVER. `ed_ecosystem_dynamics` orchestrates
# the demographic ("dynamics") timestep for one site:
#
#   zero mass-balance / flux diagnostics
#     -> IsItLoggingTime / IsItDamageTime!
#     -> ZeroAllocationRates / ZeroLitterFluxes
#     -> TotalBalanceCheck(0)
#     -> phenology (or satellite_phenology in SP mode)
#     -> [non-ST3, non-SP]:  fire_model! -> disturbance_rates
#                            -> ed_integrate_state_variables  (growth/allocation
#                                                              integration to daily)
#        [else]:             bypass_dynamics
#     -> recruitment (per patch) -> TotalBalanceCheck(1)
#     -> sort/terminate/fuse/terminate cohorts (per patch)
#     -> TotalBalanceCheck(2)
#     -> [patch dynamics on]: spawn_patches -> TotalBalanceCheck(3)
#                             -> fuse_patches -> (rhiz hydro) -> TotalBalanceCheck(4)
#                             -> terminate_patches
#     -> TotalBalanceCheck(5)
#
# `ed_integrate_state_variables` advances every cohort one daily step: mortality
# derivative, maintenance turnover (PRTMaintTurnover!), leaf ageing (AgeLeaves!),
# the three-phase daily PARTEH allocation (DailyPRT! phases 1/2/3 with damage
# recovery between 2 and 3), efflux -> litter, NPP/GPP/AR bookkeeping, height /
# dbh / cohort-age updates, then SeedUpdate, the litter-flux generation pass, and
# finally the cohort number-density update (n += dndt*dt).
#
# `ed_update_site` recomputes site-level diagnostics after the demographic step:
# canopy_spread! -> TotalBalanceCheck(6) -> canopy_structure! ->
# TotalBalanceCheck(final) -> SetRecruitL2FR -> per-patch terminate_cohorts +
# count_cohorts + area_by_age -> CheckIntegratedMassPools ->
# PrepNutrientAquisitionBCs / PrepCH4BCs -> (second-to-last day) trim_canopy.
#
# `TotalBalanceCheck` audits carbon/N/P mass conservation: change-in-stock
# (SiteMassStock) vs net (flux_in - flux_out). `bypass_dynamics` is the
# no-dynamics fast path (ST3 / SP) that just sets trivial fluxes/rates/flags.
#
# Translation notes (per project conventions):
#   * fates_r8 -> Float64; the Fortran pointer linked-list traversals
#     (`currentPatch => currentPatch%younger`, `currentCohort => currentCohort%taller`)
#     -> `while p !== nothing ... p = p.younger`.
#   * Module flag globals (hlm_use_sp / hlm_use_ed_st3 / hlm_use_planthydro /
#     hlm_use_cohort_age_tracking / hlm_use_ed_prescribed_phys / hlm_use_tree_damage /
#     hlm_use_nocomp / hlm_masterproc / hlm_parteh_mode / hlm_freq_day /
#     hlm_days_per_year / hlm_day_of_year / numpft) are Ref{Int}/Ref{Float64}
#     -> dereferenced with `[]`.
#   * `currentSite%mass_balance(el)%ZeroMassBalFlux()` -> `ZeroMassBalFlux!(site.mass_balance[el])`
#     (type-bound bang functions). Likewise ZeroFluxDiags!, GetMean(tveg24).
#   * `element_pos(carbon12_element)` is a reverse-lookup ARRAY in the Julia port
#     (`element_pos[carbon12_element]`), not a function.
#   * Fortran out-args become Julia returns: EvaluateAndCorrectDBH -> (delta_dbh,
#     delta_height); DamageRecovery -> (newly_recovered, rcohort).
#   * Bang names of the reused ported routines: fire_model! / IsItDamageTime! /
#     PRTMaintTurnover! / AgeLeaves! / DailyPRT! / canopy_structure! / canopy_spread! /
#     UpdateSizeDepPlantHydProps!.
#   * AREA/AREA_INV are the lowercase EDTypesMod consts `area` / `area_inv`.
#
# PLANT HYDRAULICS (gated, inert by default — hlm_use_planthydro defaults to ifalse):
#   * UpdateSizeDepRhizHydProps / UpdateSizeDepPlantHydStates /
#     AccumulateMortalityWaterStorage (Tier A, FatesPlantHydraulicsMod) are now wired,
#     guarded behind `hlm_use_planthydro[] == itrue`.
#
# STUBS (gated, inert by default):
#   * `fates_hist%update_history_nutrflux` (FatesHistoryInterfaceMod) is NOT ported;
#     it only runs under the prt_cnp_flex_allom_hyp (CNP) mode, which is naturally
#     inert in the default carbon-only path. Guarded + stubbed.
#
# Sibling-module calls (ported in PARALLEL in this same batch, resolve at parent
# integration time — called here by their EXACT Fortran name, NOT redefined):
#   none directly invoked by EDMainMod (the driver reaches photosynthesis /
#   radiation / inventory only indirectly through the host driver, not from these
#   four routines).
#
# Upstream-Fortran quirks preserved (see inline comments):
#   * TotalBalanceCheck error tolerance is `error_frac > 10e-6` (i.e. 1e-5, NOT
#     1e-6 — the Fortran literal is `10e-6_r8`) and triggers on NaN (`error /= error`).
#   * error_frac is only computed when `change_in_stock > 0` (else 0).
#   * Only the `final_check_id == -1` call updates `old_stock` / `err_fates`.
#   * TotalBalanceCheck is entirely skipped in SP mode.
#   * ed_integrate_state_variables ORDERING is load-bearing: maintenance turnover +
#     leaf ageing happen BEFORE EvaluateAndCorrectDBH and the 3 DailyPRT phases;
#     DamageRecovery is invoked BETWEEN phase 2 and phase 3; the n += dndt*dt
#     update is LAST (after CWD/seed input, which assume pre-mortality n).
#   * The recruitment `bc_out` arg is commented out upstream (kept as a comment).
#   * The `FluxIntoLitterPools` call is "unnecessary for CLM coupling" per RGK but
#     retained verbatim.
#   * `currentSite%transition_landuse_from_off_to_on` is unset (false) right after
#     ed_integrate_state_variables on the dynamics path.
# Standalone — NOT added to CLMInstances or any dual-copied struct.

# Module-local parameters (Fortran module parameters).
const _edmain_debug          = false   # Fortran `debug` (off)
const _edmain_final_check_id = -1       # Fortran `final_check_id`

# ===========================================================================
# ed_ecosystem_dynamics — the daily demographic driver
# ===========================================================================

"""
    ed_ecosystem_dynamics(currentSite::ed_site_type, bc_in, bc_out)

Core of the FATES demographic model: orchestrates one daily dynamics step for a
site (phenology -> fire/disturbance -> growth/allocation integration -> mortality
-> recruitment -> cohort sort/terminate/fuse -> patch spawn/fuse/terminate),
interleaved with `TotalBalanceCheck` mass-conservation audits. Bypassed (via
`bypass_dynamics`) in ST3 mode; patch dynamics are turned off in SP/ST3 modes.
Mirrors the Fortran `ed_ecosystem_dynamics`.
"""
function ed_ecosystem_dynamics(currentSite::ed_site_type, bc_in, bc_out)

    # Zero the per-element mass-balance fluxes and the flux diagnostics.
    # (Consider moving towards the end — some are integrated over the short step.)
    for el in 1:num_elements[]
        ZeroMassBalFlux!(currentSite.mass_balance[el])
    end
    ZeroFluxDiags!(currentSite.flux_diags)

    # Identify if logging / damage should occur this step (global events).
    # (IsItLoggingTime has no bang in the Julia port; IsItDamageTime! does.)
    IsItLoggingTime(hlm_masterproc[], currentSite)
    IsItDamageTime!(hlm_masterproc[])

    # ----------------------------------------------------------------------
    # Fire, growth, biogeochemistry.
    # ----------------------------------------------------------------------

    # Zero turnover rates / growth diagnostics and litter-pool fluxes.
    ZeroAllocationRates(currentSite)
    ZeroLitterFluxes(currentSite)

    # Zero the mass balance.
    TotalBalanceCheck(currentSite, 0)
    fates_parity_hook(currentSite, bc_in, 10)   # parity: dyn_in (daily-step input state)

    # Phenology. Not allowed while in ST3 mode (litter fluxes of flushing /
    # turning over leaves are not plugged in for non-dynamics runs).
    if hlm_use_ed_st3[] == ifalse
        if hlm_use_sp[] == ifalse
            phenology(currentSite, bc_in)
        else
            satellite_phenology(currentSite, bc_in)
        end
    end
    fates_parity_hook(currentSite, bc_in, 11)   # parity: after phenology

    if hlm_use_ed_st3[] == ifalse && hlm_use_sp[] == ifalse   # bypass if ST3 / SP

        # Skip the fire model if the site is solely a single bareground patch
        # (the youngest patch's patchno is 0). With multiple patches, the
        # bareground patch is avoided inside the fire model itself.
        if currentSite.youngest_patch.patchno != 0
            fire_model!(currentSite, bc_in)
        end

        # Disturbance + mortality from previous-timestep vegetation.
        disturbance_rates(currentSite, bc_in)
        fates_parity_hook(currentSite, bc_in, 12)   # parity: after fire + disturbance_rates

        # Integrate state variables from annual rates to the daily timestep.
        ed_integrate_state_variables(currentSite, bc_in, bc_out)
        fates_parity_hook(currentSite, bc_in, 13)   # parity: after growth/allocation/PRT

        # The transition flag is no longer needed once we have integrated.
        if currentSite.transition_landuse_from_off_to_on
            currentSite.transition_landuse_from_off_to_on = false
        end

    else
        # ed_integrate_state_variables is where the new-cohort flag is cleared.
        # If we are not entering it, we must clear the flag here (mark cohorts
        # as non-recruits) via bypass_dynamics.
        bypass_dynamics(currentSite)
    end

    # ----------------------------------------------------------------------
    # Reproduction, recruitment, and cohort dynamics.
    # ----------------------------------------------------------------------

    if hlm_use_ed_st3[] == ifalse && hlm_use_sp[] == ifalse
        currentPatch = currentSite.oldest_patch
        while currentPatch !== nothing
            # adds a small cohort of each PFT
            recruitment(currentSite, currentPatch, bc_in)
            # YL: call recruitment(currentSite, currentPatch, bc_in, bc_out)
            currentPatch = currentPatch.younger
        end

        TotalBalanceCheck(currentSite, 1)
        fates_parity_hook(currentSite, bc_in, 14)   # parity: after recruitment

        currentPatch = currentSite.oldest_patch
        while currentPatch !== nothing
            sort_cohorts(currentPatch)                                  # right order
            terminate_cohorts(currentSite, currentPatch, 1, 10, bc_in)  # too-few
            fuse_cohorts(currentSite, currentPatch, bc_in)              # similar
            terminate_cohorts(currentSite, currentPatch, 2, 10, bc_in)  # other reasons
            currentPatch = currentPatch.younger
        end
    end

    TotalBalanceCheck(currentSite, 2)
    fates_parity_hook(currentSite, bc_in, 15)   # parity: after cohort sort/terminate/fuse

    # ----------------------------------------------------------------------
    # Patch dynamics: fusion, new-patch spawning, termination.
    # ----------------------------------------------------------------------

    # Turn off patch dynamics if SP or ST3 modes are in use.
    do_patch_dynamics = itrue
    if hlm_use_ed_st3[] == itrue || hlm_use_sp[] == itrue
        do_patch_dynamics = ifalse
    end

    if do_patch_dynamics == itrue
        # make new patches from disturbed land
        spawn_patches(currentSite, bc_in)

        TotalBalanceCheck(currentSite, 3)
        fates_parity_hook(currentSite, bc_in, 16)   # parity: after spawn_patches

        # fuse on the spawned patches.
        fuse_patches(currentSite, bc_in)

        # If using FATES hydraulics, update the rhizosphere geometry based on the
        # new cohort-patch structure.
        if (hlm_use_planthydro[] == itrue) && do_growthrecruiteffects
            UpdateSizeDepRhizHydProps(currentSite, bc_in)
            # !! UpdateSizeDepRhizHydStates(currentSite, bc_in) (RGK 12-2021)
        end

        # SP has changes in leaf carbon but we don't expect them to balance.
        TotalBalanceCheck(currentSite, 4)

        # kill patches that are too small
        terminate_patches(currentSite, bc_in)
    end

    # Final instantaneous mass balance check.
    TotalBalanceCheck(currentSite, 5)
    fates_parity_hook(currentSite, bc_in, 17)   # parity: end of ed_ecosystem_dynamics

    return nothing
end

# ===========================================================================
# ed_integrate_state_variables — daily growth / allocation integration
# ===========================================================================

"""
    ed_integrate_state_variables(currentSite::ed_site_type, bc_in, bc_out)

Advance every cohort by one daily timestep: compute the mortality derivative,
maintenance turnover, leaf ageing, the three-phase daily PARTEH allocation (with
damage recovery between phases 2 and 3), efflux into litter, NPP/GPP/AR mass
bookkeeping, then update height / dbh / cohort-age. Across the site, also updates
recruit L2FR / stoich, runs SeedUpdate, generates the litter fluxes, and finally
updates each cohort's number density (`n += dndt*hlm_freq_day`). Mirrors the
Fortran `ed_integrate_state_variables`.
"""
function ed_integrate_state_variables(currentSite::ed_site_type, bc_in, bc_out)

    current_fates_landuse_state_vector = get_current_landuse_statevector(currentSite)

    # Clear site GPP and AR passing to the HLM.
    bc_out.gpp_site = 0.0
    bc_out.ar_site  = 0.0

    # Patch-level biomass required for C-based harvest.
    harvestable_forest_c = get_harvestable_carbon(currentSite, bc_in.site_area,
                                                  bc_in.hlm_harvest_catnames)
    harvest_tag = fill(2, hlm_num_lu_harvest_cats[])

    # Pointer to this site's carbon12 mass balance.
    site_cmass = currentSite.mass_balance[element_pos[carbon12_element]]

    # Update the total stoichiometry assessment for a new recruit (prior to the
    # growth sequence where reproductive tissues are allocated).
    UpdateRecruitStoich(currentSite)

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing

        currentPatch.age += hlm_freq_day[]
        if currentPatch.age < 0.0
            fates_endrun("negative patch age? $(currentPatch.age) " *
                         "patchno=$(currentPatch.patchno) area=$(currentPatch.area)")
        end

        # Age increment for secondary forest patches.
        if currentPatch.land_use_label != primaryland
            currentPatch.age_since_anthro_disturbance += hlm_freq_day[]
        end

        # Has the patch moved to the next age class?
        currentPatch.age_class = get_age_class_index(currentPatch.age)

        # Within this loop we may create new cohorts (copies at reduced damage
        # class) and want to bypass some of the things here for them.
        newly_recovered = false

        currentCohort = currentPatch.shortest
        while currentCohort !== nothing

            ft = currentCohort.pft

            if !newly_recovered

                # Mortality derivatives.
                mean_temp = GetMean(currentPatch.tveg24)
                Mortality_Derivative(currentSite, currentCohort, bc_in,
                    currentPatch.btran_ft, mean_temp,
                    currentPatch.land_use_label,
                    currentPatch.age_since_anthro_disturbance,
                    current_fates_landuse_state_vector[primaryland],
                    current_fates_landuse_state_vector[secondaryland],
                    harvestable_forest_c, harvest_tag)

                # ----------------------------------------------------------
                # Identify net carbon gain for this dynamics interval.
                # ----------------------------------------------------------
                if hlm_use_ed_prescribed_phys[] == itrue
                    edpft = edpftvarcon_inst()
                    if currentCohort.canopy_layer == 1
                        currentCohort.npp_acc = edpft.prescribed_npp_canopy[ft] *
                            currentCohort.c_area / currentCohort.n / hlm_days_per_year[]
                    else
                        currentCohort.npp_acc = edpft.prescribed_npp_understory[ft] *
                            currentCohort.c_area / currentCohort.n / hlm_days_per_year[]
                    end
                    # No explicit respiration for prescribed phys; pass mass
                    # balance by saying respiration is zero.
                    currentCohort.gpp_acc  = currentCohort.npp_acc
                    currentCohort.resp_acc = 0.0
                end

                # ----------------------------------------------------------
                # Save NPP/GPP/R into the "hold" diagnostics (kgC/indiv/year),
                # which persist for I/O after the _acc vars are zeroed.
                # ----------------------------------------------------------
                currentCohort.npp_acc_hold  = currentCohort.npp_acc  * Float64(hlm_days_per_year[])
                currentCohort.gpp_acc_hold  = currentCohort.gpp_acc  * Float64(hlm_days_per_year[])
                currentCohort.resp_acc_hold = currentCohort.resp_acc * Float64(hlm_days_per_year[])

                # Pass gpp/ar to the HLM.
                bc_out.gpp_site += currentCohort.gpp_acc_hold * area_inv *
                    currentCohort.n / hlm_days_per_year[] / sec_per_day
                bc_out.ar_site  += currentCohort.resp_acc_hold * area_inv *
                    currentCohort.n / hlm_days_per_year[] / sec_per_day

                # Maintenance turnover (PARTEH).
                if _edmain_debug
                    CheckMassConservation(currentCohort.prt, ft, 3)
                end
                if any(currentSite.dstatus[ft] .== (phen_dstat_moiston, phen_dstat_timeon))
                    is_drought = false
                else
                    is_drought = true
                end

                PRTMaintTurnover!(currentCohort.prt, ft, currentCohort.canopy_layer, is_drought)

                # Advance leaves in age (move a portion of leaf mass into the
                # next age bin; not the oldest-bin -> litter movement).
                AgeLeaves!(currentCohort.prt, ft, currentCohort.canopy_layer, sec_per_day)

                # Combine N sources into daily_n_gain (used by allocation).
                currentCohort.daily_n_gain = currentCohort.daily_nh4_uptake +
                    currentCohort.daily_no3_uptake + currentCohort.sym_nfix_daily

                currentCohort.resp_excess = 0.0

            end  # if !newly_recovered

            # Correct dbh if it is below what is allometrically consistent with
            # the structural biomass. (delta_dbh/delta_height returned but unused.)
            _, _ = EvaluateAndCorrectDBH(currentCohort)

            # Save pre-PRT values (also for the newly recovered cohort).
            height_old = currentCohort.height
            dbh_old    = currentCohort.dbh

            # ----------------------------------------------------------------
            # Growth and allocation (PARTEH), in three phases.
            #   phase 1: prioritized once-per-day allocation (no stature growth)
            #   phase 2: repeatable non-stature allocation (recovering plants
            #            update targets here)
            #   phase 3: stature growth from any left-over resources
            # ----------------------------------------------------------------
            # Refresh the PARTEH boundary-condition Refs from the live cohort fields
            # before allocating. The port backs each BC by a Ref captured at cohort
            # registration (Julia stand-in for the Fortran aliasing pointer); nothing
            # else refreshes them, so without this DailyPRT reads a STALE carbon input
            # (netdc/npp_acc ≡ cold-start 0) → allocates zero → biomass never grows
            # (dbh frozen) and the fixed carbon (gpp_acc) has no sink → mass-balance error.
            sync_cohort_to_prt_bcs!(currentCohort)

            if !newly_recovered
                DailyPRT!(currentCohort.prt, 1)
            end

            DailyPRT!(currentCohort.prt, 2)

            if (!newly_recovered) && (hlm_use_tree_damage[] == itrue)
                # The recovered cohort (with larger targets) is created in
                # DamageRecovery and inserted into the next (taller) position.
                newly_recovered, _ = DamageRecovery(currentSite, currentPatch, currentCohort)
            else
                newly_recovered = false
            end

            DailyPRT!(currentCohort.prt, 3)

            # Copy the updated inout/out Refs (grown dbh, spent netdc, effluxes, CNP
            # state) BACK into the cohort fields — DailyPRT wrote the Refs; the
            # height/allometry/diagnostics/mass-balance below read the fields.
            sync_prt_bcs_to_cohort!(currentCohort)

            # Update mass-balance tracking for the daily nutrient uptake flux,
            # then zero the daily uptakes (they have been used).
            EffluxIntoLitterPools(currentSite, currentPatch, currentCohort, bc_in)

            if element_pos[nitrogen_element] > 0
                currentSite.mass_balance[element_pos[nitrogen_element]].net_root_uptake +=
                    (currentCohort.daily_n_gain - currentCohort.daily_n_efflux) * currentCohort.n
            end
            if element_pos[phosphorus_element] > 0
                currentSite.mass_balance[element_pos[phosphorus_element]].net_root_uptake +=
                    (currentCohort.daily_p_gain - currentCohort.daily_p_efflux) * currentCohort.n
            end

            # Mass balance for C efflux (if any).
            currentSite.mass_balance[element_pos[carbon12_element]].net_root_uptake -=
                currentCohort.daily_c_efflux * currentCohort.n

            # Save NPP diagnostic for flux accounting [kg/m2/day].
            currentSite.flux_diags.npp +=
                (currentCohort.npp_acc_hold / hlm_days_per_year[] - currentCohort.resp_excess) *
                currentCohort.n * area_inv

            # Add the input fluxes to the mass-balance accounting.
            site_cmass.gpp_acc   += currentCohort.gpp_acc * currentCohort.n
            site_cmass.aresp_acc += (currentCohort.resp_acc + currentCohort.resp_excess) *
                currentCohort.n

            CheckMassConservation(currentCohort.prt, ft, 5)

            # Update the leaf biophysical rates from the proportion of leaf mass
            # in the different age classes (won't change again after growth/turnover).
            UpdateCohortBioPhysRates(currentCohort)

            # This cohort has grown; it is no longer "new".
            currentCohort.isnew = false

            # Update plant height (if it has grown).
            # `h_allom` returns (height, dh/ddbh) — the Fortran is
            #     call h_allom(currentCohort%dbh, ft, currentCohort%height)
            # i.e. the HEIGHT is the output. The destructuring here used to be
            # reversed (`_, height = ...`), which assigned the allometric DERIVATIVE
            # dh/ddbh to cohort height on EVERY daily growth step — height then drifted
            # away from its own allometry (e.g. a PFT-1 sapling at dbh=0.764 cm has
            # h=1.30 m and dh/ddbh=1.36, so height jumped 1.30 -> 1.36 m on day one with
            # ddbhdt == 0) and kept climbing, corrupting canopy layering, crown area and
            # the tallest->shortest cohort ordering. Every other h_allom call site in the
            # port already uses `h, _ = h_allom(...)`. Caught by the time-stepped
            # Fortran-FATES parity harness (scripts/fates_fortran_parity.jl).
            currentCohort.height, _ = h_allom(currentCohort.dbh, ft)

            currentCohort.dhdt   = (currentCohort.height - height_old) / hlm_freq_day[]
            currentCohort.ddbhdt = (currentCohort.dbh - dbh_old) / hlm_freq_day[]

            # Carbon assimilate has been spent; safe to zero.
            currentCohort.npp_acc  = 0.0
            currentCohort.gpp_acc  = 0.0
            currentCohort.resp_acc = 0.0

            # Update tree hydraulic geometry.
            if (hlm_use_planthydro[] == itrue) && do_growthrecruiteffects
                UpdateSizeDepPlantHydProps!(currentSite, currentCohort, bc_in)
                UpdateSizeDepPlantHydStates(currentSite, currentCohort)
            end

            # Age-dependent mortality mode: update cohort age + age classes.
            if hlm_use_cohort_age_tracking[] == itrue
                currentCohort.coage += hlm_freq_day[]
                if currentCohort.coage < 0.0
                    fates_endrun("negative cohort age? $(currentCohort.coage)")
                end
                currentCohort.coage_class, currentCohort.coage_by_pft_class =
                    coagetype_class_index(currentCohort.coage, currentCohort.pft)
            end

            currentCohort = currentCohort.taller
        end  # cohort loop

        currentPatch = currentPatch.younger
    end  # patch loop

    # Record the L2FRs of near-recruit-size plants (per pft / canopy layer) to
    # set the L2FRs of newly recruited plants.
    UpdateRecruitL2FR(currentSite)

    # Update Nutrient history diagnostics (if any).
    if hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
        # TODO Batch NN: fates_hist%update_history_nutrflux not ported
        # (FatesHistoryInterfaceMod). Inert in the default carbon-only mode.
        # update_history_nutrflux(fates_hist, currentSite)
    end

    # When plants die, the water goes with them (affects the water balance).
    if hlm_use_planthydro[] == itrue
        currentPatch = currentSite.youngest_patch
        while currentPatch !== nothing
            currentCohort = currentPatch.shortest
            while currentCohort !== nothing
                AccumulateMortalityWaterStorage(currentSite, currentCohort,
                    -1.0 * currentCohort.dndt * hlm_freq_day[])
                currentCohort = currentCohort.taller
            end
            currentPatch = currentPatch.older
        end
    end

    # With growth + mortality rates calculated, determine the seed-rain fluxes.
    # (Potentially cross-patch mixing, so it is calculated as a group.)
    SeedUpdate(currentSite)

    # Calculate all other litter fluxes.
    currentPatch = currentSite.youngest_patch
    while currentPatch !== nothing
        GenerateDamageAndLitterFluxes(currentSite, currentPatch, bc_in)
        PreDisturbanceLitterFluxes(currentSite, currentPatch, bc_in)
        PreDisturbanceIntegrateLitter(currentPatch)
        currentPatch = currentPatch.older
    end

    # RGK: unnecessary for CLM coupling (could move to ELM's UpdateLitterFluxes),
    # but retained verbatim.
    FluxIntoLitterPools(currentSite, bc_in, bc_out)

    # Update cohort number. Must happen AFTER the CWD/seed input calculations as
    # they assume the pre-mortality currentCohort.n.
    currentPatch = currentSite.youngest_patch
    while currentPatch !== nothing
        currentCohort = currentPatch.shortest
        while currentCohort !== nothing
            currentCohort.n = max(0.0, currentCohort.n + currentCohort.dndt * hlm_freq_day[])
            currentCohort.sym_nfix_daily = 0.0
            currentCohort = currentCohort.taller
        end
        currentPatch = currentPatch.older
    end

    return nothing
end

# ===========================================================================
# ed_update_site — recompute site-level diagnostics after the demographic step
# ===========================================================================

"""
    ed_update_site(currentSite::ed_site_type, bc_in, bc_out, is_restarting::Bool)

Consolidate the ED growth process: update canopy spread + structure (assign
cohorts to canopy layers), audit balance, set recruit L2FRs, then per patch
terminate spurious cohorts, count cohorts (for the photosynthesis loop) and
accumulate `area_by_age`. Then check the integrated mass pools, prepare the
nutrient-acquisition / CH4 boundary conditions to the HLM, and (on the
second-to-last day of the year) trim the canopy. Skips the dynamics-only steps
when `is_restarting`. Mirrors the Fortran `ed_update_site`.
"""
function ed_update_site(currentSite::ed_site_type, bc_in, bc_out, is_restarting::Bool)

    if hlm_use_sp[] == ifalse && !is_restarting
        canopy_spread!(currentSite)
    end

    TotalBalanceCheck(currentSite, 6)

    if hlm_use_sp[] == ifalse && !is_restarting
        canopy_structure!(currentSite, bc_in)
    end

    TotalBalanceCheck(currentSite, _edmain_final_check_id)

    # Update recruit L2FRs based on the new canopy position.
    SetRecruitL2FR(currentSite)

    fill!(currentSite.area_by_age, 0.0)

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing

        if !is_restarting
            terminate_cohorts(currentSite, currentPatch, 1, 11, bc_in)
            terminate_cohorts(currentSite, currentPatch, 2, 11, bc_in)
        end

        # Used in the photosynthesis loop.
        count_cohorts(currentPatch)

        # Update the per-age-class area total.
        currentSite.area_by_age[currentPatch.age_class] += currentPatch.area

        currentPatch = currentPatch.younger
    end

    # Check that the time-integrated fluxes match the state. Don't call when
    # restarting (it would double count the flux).
    if !is_restarting
        CheckIntegratedMassPools(currentSite)
    end

    # The HLMs need to know about nutrient demand / root mass + affinities.
    PrepNutrientAquisitionBCs(currentSite, bc_in, bc_out)

    # The HLM methane module needs rooting mass / distributions / respiration / NPP.
    PrepCH4BCs(currentSite, bc_in, bc_out)

    # FIX(RF): this should be monthly, not annual. On the second-to-last day of
    # the year, perform trimming.
    if hlm_day_of_year[] == hlm_days_per_year[] - 1 && !is_restarting
        if hlm_use_sp[] == ifalse
            trim_canopy(currentSite)
        end
    end

    is_restarting || fates_parity_hook(currentSite, bc_in, 18)   # parity: end of daily step

    return nothing
end

# ===========================================================================
# TotalBalanceCheck — carbon / N / P mass-conservation audit
# ===========================================================================

"""
    TotalBalanceCheck(currentSite::ed_site_type, call_index::Integer)

Compare the mass flux in/out of FATES against the change in total stocks (states)
for every element. Fluxes in are seed-in + net root uptake + GPP + generic-in +
patch-resize-err; fluxes out are wood-product harvest/LUC + burn + seed-out +
generic-out + fragmentation + autotrophic respiration. Errors beyond the
tolerance (`error_frac > 10e-6`, i.e. 1e-5) or any NaN abort the run with a full
per-patch / per-cohort dump. The `call_index == -1` (final) call updates the
remembered `old_stock` / `err_fates`. Entirely skipped in SP mode. Mirrors the
Fortran `TotalBalanceCheck`.

Upstream quirk: the Fortran tolerance literal is `10e-6_r8` (= 1e-5), NOT 1e-6;
preserved verbatim. `error_frac` is only computed when `change_in_stock > 0`.
"""
function TotalBalanceCheck(currentSite::ed_site_type, call_index::Integer)

    hlm_use_sp[] == ifalse || return nothing   # skip entirely in SP mode

    change_in_stock = 0.0

    for el in 1:num_elements[]

        site_mass = currentSite.mass_balance[el]

        total_stock, biomass_stock, litter_stock, seed_stock =
            SiteMassStock(currentSite, el)

        change_in_stock = total_stock - site_mass.old_stock

        flux_in = site_mass.seed_in +
                  site_mass.net_root_uptake +
                  site_mass.gpp_acc +
                  site_mass.flux_generic_in +
                  site_mass.patch_resize_err

        flux_out = sum(site_mass.wood_product_harvest) +
                   sum(site_mass.wood_product_landusechange) +
                   site_mass.burn_flux_to_atm +
                   site_mass.seed_out +
                   site_mass.flux_generic_out +
                   site_mass.frag_out +
                   site_mass.aresp_acc

        net_flux = flux_in - flux_out
        error    = abs(net_flux - change_in_stock)

        if change_in_stock > 0.0
            error_frac = error / abs(total_stock)
        else
            error_frac = 0.0
        end

        # Upstream tolerance literal: 10e-6_r8 (= 1e-5), NOT 1e-6. NaN also fails.
        if error_frac > 10e-6 || isnan(error)
            @info "mass balance error detected"
            @info "element type (see PRTGenericMod): $(element_list[el])"
            @info "error fraction relative to biomass stock: $(error_frac)"
            @info "absolute error (flux in - change): $(net_flux - change_in_stock)"
            @info "call index: $(call_index)"
            @info "net: $(net_flux)"
            @info "dstock: $(change_in_stock)"
            @info "seed_in: $(site_mass.seed_in)"
            @info "net_root_uptake: $(site_mass.net_root_uptake)"
            @info "gpp_acc: $(site_mass.gpp_acc)"
            @info "flux_generic_in: $(site_mass.flux_generic_in)"
            @info "wood_product_harvest: $(site_mass.wood_product_harvest)"
            @info "wood_product_landusechange: $(site_mass.wood_product_landusechange)"
            @info "error from patch resizing: $(site_mass.patch_resize_err)"
            @info "burn_flux_to_atm: $(site_mass.burn_flux_to_atm)"
            @info "seed_out: $(site_mass.seed_out)"
            @info "flux_generic_out: $(site_mass.flux_generic_out)"
            @info "frag_out: $(site_mass.frag_out)"
            @info "aresp_acc: $(site_mass.aresp_acc)"
            @info "error=net_flux-dstock: $(error)"
            @info "biomass: $(biomass_stock)"
            @info "litter: $(litter_stock)"
            @info "seeds: $(seed_stock)"
            @info "total stock: $(total_stock)"
            @info "previous total: $(site_mass.old_stock)"
            @info "lat lon: $(currentSite.lat) $(currentSite.lon)"

            # Per-patch / per-cohort dump (upstream gates this on the first day of
            # simulation, but that guard is commented out in the Fortran).
            currentPatch = currentSite.oldest_patch
            while currentPatch !== nothing
                litt = currentPatch.litter[el]
                @info "---------------------------------------"
                @info "patch area: $(currentPatch.area)"
                @info "AG CWD: $(sum(litt.ag_cwd))"
                @info "BG CWD (sum): $(sum(litt.bg_cwd))"
                @info "leaf litter: $(sum(litt.leaf_fines))"
                @info "root litter (sum): $(sum(litt.root_fines))"
                @info "land_use_label: $(currentPatch.land_use_label)"

                currentCohort = currentPatch.tallest
                while currentCohort !== nothing
                    @info "pft: $(currentCohort.pft) dbh: $(currentCohort.dbh)"
                    leaf_m   = GetState(currentCohort.prt, leaf_organ,   element_list[el])
                    struct_m = GetState(currentCohort.prt, struct_organ, element_list[el])
                    store_m  = GetState(currentCohort.prt, store_organ,  element_list[el])
                    fnrt_m   = GetState(currentCohort.prt, fnrt_organ,   element_list[el])
                    repro_m  = GetState(currentCohort.prt, repro_organ,  element_list[el])
                    sapw_m   = GetState(currentCohort.prt, sapw_organ,   element_list[el])
                    @info "leaf: $(leaf_m) structure: $(struct_m) store: $(store_m)"
                    @info "fineroot: $(fnrt_m) repro: $(repro_m) sapwood: $(sapw_m)"
                    @info "num plant: $(currentCohort.n)"
                    currentCohort = currentCohort.shorter
                end

                currentPatch = currentPatch.younger
            end
            fates_endrun("TotalBalanceCheck: mass balance error (call $(call_index)) " *
                         "on date $(hlm_current_year[]) $(hlm_current_month[]) $(hlm_current_day[])")
        end

        # Last check of the sequence: update the total error check + final stock.
        if call_index == _edmain_final_check_id
            site_mass.old_stock = total_stock
            site_mass.err_fates = net_flux - change_in_stock
        end

    end  # element loop

    return nothing
end

# ===========================================================================
# bypass_dynamics — the no-dynamics fast path (ST3 / SP)
# ===========================================================================

"""
    bypass_dynamics(currentSite::ed_site_type)

When dynamics are bypassed (ST3 mode), set the various cohort fluxes, rates and
flags to trivial values: clear the new-cohort flag, roll npp/gpp/resp into the
"_hold" diagnostics and zero the accumulators, and zero all the mortality rates
and the size derivatives. WARNING (upstream): turning off dynamics is
experimental and this trivial-value setting may not be complete. Mirrors the
Fortran `bypass_dynamics`.
"""
function bypass_dynamics(currentSite::ed_site_type)
    currentPatch = currentSite.youngest_patch
    while currentPatch !== nothing
        currentCohort = currentPatch.shortest
        while currentCohort !== nothing

            currentCohort.isnew = false

            currentCohort.npp_acc_hold  = currentCohort.npp_acc  * Float64(hlm_days_per_year[])
            currentCohort.gpp_acc_hold  = currentCohort.gpp_acc  * Float64(hlm_days_per_year[])
            currentCohort.resp_acc_hold = currentCohort.resp_acc * Float64(hlm_days_per_year[])

            currentCohort.npp_acc  = 0.0
            currentCohort.gpp_acc  = 0.0
            currentCohort.resp_acc = 0.0

            # The "net_art" terms are zeroed at the beginning of the daily step;
            # if DailyPRT / maintenance / phenology aren't called they stay zero.

            currentCohort.bmort  = 0.0
            currentCohort.hmort  = 0.0
            currentCohort.cmort  = 0.0
            currentCohort.frmort = 0.0
            currentCohort.smort  = 0.0
            currentCohort.asmort = 0.0
            currentCohort.dgmort = 0.0

            currentCohort.dndt   = 0.0
            currentCohort.dhdt   = 0.0
            currentCohort.ddbhdt = 0.0

            # Nutrient fluxes should already be zero (no uptake in ST3 mode).

            currentCohort = currentCohort.taller
        end
        currentPatch = currentPatch.older
    end
    return nothing
end
