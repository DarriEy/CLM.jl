# EDPatchDynamicsMod.jl
# Julia port of FATES src/fates/biogeochem/EDPatchDynamicsMod.F90
#
# The PATCH DISTURBANCE ENGINE of FATES — controls formation, creation, fusing
# and termination of patch-level processes on the age-ordered patch linked list
# (`oldest_patch`/`youngest_patch` connected by each patch's `older`/`younger`
# pointers). It:
#   * disturbance_rates        — per-patch fire/treefall-mortality/logging/land-use
#                                disturbance rates (per daily timestep).
#   * spawn_patches            — resolve disturbance: create new disturbed patches,
#                                move surviving/killed cohorts + litter into them.
#   * split_patch              — split a patch into two patches identical but for
#                                area (cohorts split by number density).
#   * check_patch_area         — assert patch areas sum to the site area, fixing
#                                tiny precision drift onto the largest patch.
#   * set_patchno              — number patches oldest->youngest (bareground=0).
#   * TransLitterNewPatch      — move existing donor litter (partial burn) into a
#                                new patch, conserving mass.
#   * fire_litter_fluxes       — fire-killed plant biomass + partial burn -> litter
#                                / atmosphere.
#   * mortality_litter_fluxes  — treefall (gap + impact) killed biomass -> CWD.
#   * landusechange_litter_fluxes — felled-on-LUC biomass -> CWD / wood products /
#                                atmosphere (only when the clearing matrix says so).
#   * fuse_patches             — fuse patches with similar size profiles (per land
#                                use category and nocomp PFT), area-conserving.
#   * fuse_2_patches           — fuse one donor patch into a recipient and relink.
#   * terminate_patches        — fuse away patches that are too small.
#   * DistributeSeeds          — spread a seed mass over all patches' seed pools.
#   * patch_pft_size_profile   — binned AGB profile of a patch (for fusion).
#   * countPatches             — census of patches across sites.
#   * GetPseudoPatchAge        — sort key encoding land-use/pft label + age.
#  plus the helpers InsertPatch and CopyPatchMeansTimers used by spawn/split.
#
# Translation notes (per project conventions):
#   * `fates_r8` -> Float64; cohort/patch/site/litter structs are REUSED VERBATIM
#     from FatesCohortMod / FatesPatchMod / EDTypesMod / FatesLitterMod — this
#     module does NOT redefine any of them.
#   * Fortran `allocate(newPatch)`/`allocate(nc)` + pointer splice -> construct an
#     empty `fates_patch_type()` / `fates_cohort_type()` and walk/relink the
#     `Union{...,Nothing}` linked-list pointers. `deallocate(...)` -> Julia is
#     GC'd; the node is already unlinked + FreeMemory()'d, so we just drop it.
#   * Cohort-list head/tail pass-by-pointer updates use `Ref`s for
#     insert_cohort's storebigcohort/storesmallcohort, exactly as the ported
#     EDCohortDynamicsMod callers do.
#   * Module flag globals (hlm_use_planthydro / hlm_use_nocomp / hlm_use_luh /
#     hlm_use_fixed_biogeog / hlm_freq_day / hlm_current_tod / numpft / num_swb)
#     are `Ref`s -> dereferenced with `[]`. `regeneration_model` reads from
#     `ed_params().regeneration_model`. `prt_params.woody[...]` etc. are direct.
#   * Reused ported helpers (NOT reimplemented): create_cohort, terminate_cohorts,
#     fuse_cohorts, sort_cohorts, insert_cohort, Copy/ZeroValues/FreeMemory
#     (EDCohortDynamicsMod / FatesCohortMod, Batch 12); mortality_rates /
#     ExemptTreefallDist (EDMortalityFunctionsMod, Batch 13); logging_litter_fluxes!
#     / logging_time / get_harvest_rate_* / get_harvestable_carbon /
#     get_harvest_debt! / LoggingMortality_frac (EDLoggingMortalityMod, Batch 12);
#     carea_allom / set_root_fraction (FatesAllometryMod); get_age_class_index
#     (FatesSizeAgeTypeIndicesMod); GetState / PRTBurnLosses! (PRT*); adjust_SF_CWD_frac
#     / FuseLitter! / InitConditions! (FatesLitterMod); PatchMassStock
#     (ChecksBalancesMod); fuel fuse!/frac_burnt, fuel_classes accessors
#     (FatesFuelMod / FatesFuelClassesMod); the running-mean FuseRMean!/CopyFromDonor!.
#   * GetLanduseChangeRules / GetLanduseTransitionRates! / GetInitLanduseTransitionRates!
#     / GetLUHStatedata / GetInitLanduseHarvestRate (FatesLandUseChangeMod) are
#     ported and wired; all the land-use / LUH paths are inert when hlm_use_luh
#     defaults off.
#
# Stubs (clearly flagged `# TODO Batch NN:` inline):
#   * hlm_use_planthydro==itrue paths call AccumulateMortalityWaterStorage /
#     InitHydrCohort — only partially ported in FatesPlantHydraulicsMod; guarded,
#     inert by default (hlm_use_planthydro defaults to ifalse). Same convention as
#     Batches 12-13.
#
# Upstream-Fortran quirks preserved (see inline comments):
#   * existing_litt_localization=1 / treefall=burn=0 / landusechange=1 (litter
#     localization constants).
#   * fire_litter_fluxes early-returns if there was no fire in the donor patch.
#   * mortality_litter_fluxes' `DistributeSeeds` call is commented out upstream
#     ("BREAKING MASS CONSERVATION RIGHT NOW") — kept commented, seed mass goes
#     straight to the seed litter pool instead.
#   * fuse_2_patches: "rp%area = rp%area + dp%area  !THIS MUST COME AT THE END!".
#   * landusechange_litter_fluxes does nothing unless clearing_matrix_element.
#   * GetPseudoPatchAge negates the (label-weighted) categorical term so lower
#     labels sort "younger" (higher patchno).
# Standalone — NOT added to CLMInstances or any dual-copied struct.

# ---------------------------------------------------------------------------
# Litter localization constants (Fortran module parameters)
# ---------------------------------------------------------------------------
# When creating new patches, how much of the litter from the old patch (and from
# plants killed in the disturbance) is sent to the NEW patch vs. retained by the
# donor. localization=1 => all to the new patch; 0 => dispensed randomly across
# the combined area.
const existing_litt_localization = 1.0
const treefall_localization      = 0.0
const burn_localization          = 0.0
const landusechange_localization = 1.0

# ===========================================================================
# disturbance_rates
# ===========================================================================

"""
    disturbance_rates(site_in, bc_in)

Calculate the fire / treefall-mortality / logging / land-use disturbance rates
for each patch (all per daily timestep). First re-derives every cohort's
mortality + logging-mortality fractions (via the ported `mortality_rates` /
`LoggingMortality_frac`), then accumulates the per-patch `disturbance_rates`
(by `dtype_*`) and `landuse_transition_rates`, reducing them proportionally if
they would exceed the patch area in one day. Mirrors the Fortran
`disturbance_rates`.
"""
function disturbance_rates(site_in::ed_site_type, bc_in::bc_in_type)
    max_daily_disturbance_rate = 0.999

    # Fraction of the site that is each land use category.
    current_fates_landuse_state_vector = get_current_landuse_statevector(site_in)

    # Fraction of secondary land that is young secondary land.
    secondary_young_fraction = get_secondary_young_fraction(site_in)

    # Error-check the transition_landuse_from_off_to_on flag.
    if site_in.transition_landuse_from_off_to_on
        if sum(@view current_fates_landuse_state_vector[secondaryland:cropland]) > nearzero
            fates_endrun("transition_landuse_from_off_to_on true but site is not entirely primaryland")
        end
    end

    # Available biomass for harvest for all patches.
    harvestable_forest_c = get_harvestable_carbon(site_in, bc_in.site_area,
                                                  bc_in.hlm_harvest_catnames)

    # --- Per-cohort mortality + logging-mortality fractions -----------------
    currentPatch = site_in.oldest_patch
    while currentPatch !== nothing
        currentCohort = currentPatch.shortest
        while currentCohort !== nothing
            mean_temp = GetMean(currentPatch.tveg24)
            cmort, hmort, bmort, frmort, smort, asmort, dgmort =
                mortality_rates(currentCohort, bc_in, currentPatch.btran_ft, mean_temp)
            currentCohort.dmort = cmort + hmort + bmort + frmort + smort + asmort + dgmort
            _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
                site_in.spread, currentCohort.pft, currentCohort.crowndamage)

            # Diagnostic mortality rates.
            currentCohort.cmort  = cmort
            currentCohort.bmort  = bmort
            currentCohort.hmort  = hmort
            currentCohort.frmort = frmort
            currentCohort.smort  = smort
            currentCohort.asmort = asmort
            currentCohort.dgmort = dgmort

            lmort_direct, lmort_collateral, lmort_infra, l_degrad, _ =
                LoggingMortality_frac(site_in, bc_in, currentCohort.pft,
                    currentCohort.dbh, currentCohort.canopy_layer,
                    bc_in.hlm_harvest_rates, bc_in.hlm_harvest_catnames,
                    bc_in.hlm_harvest_units, currentPatch.land_use_label,
                    currentPatch.age_since_anthro_disturbance,
                    current_fates_landuse_state_vector[primaryland],
                    current_fates_landuse_state_vector[secondaryland],
                    harvestable_forest_c)

            currentCohort.lmort_direct     = lmort_direct
            currentCohort.lmort_collateral = lmort_collateral
            currentCohort.lmort_infra      = lmort_infra
            currentCohort.l_degrad         = l_degrad

            currentCohort = currentCohort.taller
        end
        currentPatch = currentPatch.younger
    end

    get_harvest_debt!(site_in, bc_in, fill(2, hlm_num_lu_harvest_cats[]))

    if hlm_use_luh[] == itrue
        if !site_in.transition_landuse_from_off_to_on
            GetLanduseTransitionRates!(bc_in, site_in.min_allowed_landuse_fraction,
                site_in.landuse_transition_matrix, site_in.landuse_vector_gt_min)
        else
            GetInitLanduseTransitionRates!(bc_in, site_in.min_allowed_landuse_fraction,
                site_in.landuse_transition_matrix, site_in.landuse_vector_gt_min)
        end
    else
        fill!(site_in.landuse_transition_matrix, 0.0)
    end

    # --- Recalculate total canopy area prior to resolving disturbance -------
    currentPatch = site_in.oldest_patch
    while currentPatch !== nothing
        currentPatch.total_canopy_area = 0.0
        currentCohort = currentPatch.shortest
        while currentCohort !== nothing
            if currentCohort.canopy_layer == 1
                currentPatch.total_canopy_area += currentCohort.c_area
            end
            currentCohort = currentCohort.taller
        end
        currentPatch = currentPatch.younger
    end

    # Info needed to decide whether to apply land-use change.
    site_secondaryland_first_exceeding_min = false
    if hlm_use_luh[] == itrue
        state_vector = GetLUHStatedata(bc_in)
        site_secondaryland_first_exceeding_min =
            (state_vector[secondaryland] > site_in.min_allowed_landuse_fraction) &&
            (!site_in.landuse_vector_gt_min[secondaryland])
    else
        state_vector = current_fates_landuse_state_vector
    end

    currentPatch = site_in.oldest_patch
    while currentPatch !== nothing
        currentPatch.disturbance_rates[dtype_ifall] = 0.0
        currentPatch.disturbance_rates[dtype_ilog]  = 0.0
        currentPatch.disturbance_rates[dtype_ifire] = 0.0

        dist_rate_ldist_notharvested = 0.0

        # Convert the transition matrix to area per unit area of that land-use
        # type per time. Avoid the divide-by-zero for off/bareground.
        if hlm_use_luh[] == itrue && currentPatch.land_use_label > nocomp_bareground_land
            for k in 1:n_landuse_cats
                currentPatch.landuse_transition_rates[k] = min(1.0,
                    site_in.landuse_transition_matrix[currentPatch.land_use_label, k] *
                    (1.0 - site_in.area_bareground) /
                    current_fates_landuse_state_vector[currentPatch.land_use_label])
            end
        else
            fill!(currentPatch.landuse_transition_rates, 0.0)
        end

        currentCohort = currentPatch.shortest
        while currentCohort !== nothing
            if currentCohort.canopy_layer == 1
                # Treefall disturbance rate (trees only, not grasses).
                if !ExemptTreefallDist(currentCohort)
                    currentPatch.disturbance_rates[dtype_ifall] +=
                        fates_mortality_disturbance_fraction() *
                        min(1.0, currentCohort.dmort) * hlm_freq_day[] *
                        currentCohort.c_area / currentPatch.area
                end

                # Logging disturbance rate.
                currentPatch.disturbance_rates[dtype_ilog] +=
                    min(1.0, currentCohort.lmort_direct + currentCohort.lmort_collateral +
                        currentCohort.lmort_infra + currentCohort.l_degrad) *
                    currentCohort.c_area / currentPatch.area

                # Non-harvested part of the logging disturbance rate.
                dist_rate_ldist_notharvested += currentCohort.l_degrad *
                    currentCohort.c_area / currentPatch.area
            end
            currentCohort = currentCohort.taller
        end

        # For non-closed canopies subject to logging, add the interstitial
        # ground-area transfer to new secondary lands.
        if (logging_time[] || site_in.transition_landuse_from_off_to_on ||
            site_secondaryland_first_exceeding_min) &&
           (currentPatch.area - currentPatch.total_canopy_area) > fates_tiny

            if !site_in.transition_landuse_from_off_to_on
                if bc_in.hlm_harvest_units == hlm_harvest_carbon
                    harvest_tag = fill(2, hlm_num_lu_harvest_cats[])
                    harvest_rate = get_harvest_rate_carbon(currentPatch.land_use_label,
                        bc_in.hlm_harvest_catnames, bc_in.hlm_harvest_rates,
                        currentPatch.age_since_anthro_disturbance, harvestable_forest_c,
                        harvest_tag)
                else
                    harvest_rate = get_harvest_rate_area(currentPatch.land_use_label,
                        bc_in.hlm_harvest_catnames, bc_in.hlm_harvest_rates,
                        current_fates_landuse_state_vector[primaryland],
                        current_fates_landuse_state_vector[secondaryland],
                        secondary_young_fraction, currentPatch.age_since_anthro_disturbance)
                end

                if state_vector[secondaryland] <= site_in.min_allowed_landuse_fraction
                    harvest_rate = 0.0
                elseif currentPatch.land_use_label == primaryland &&
                       !site_in.landuse_vector_gt_min[secondaryland]
                    harvest_rate = state_vector[secondaryland] / sum(state_vector)
                else
                    harvest_rate = 0.0
                end
            else
                harvest_rate = GetInitLanduseHarvestRate(bc_in,
                    site_in.min_allowed_landuse_fraction, site_in.landuse_vector_gt_min)
            end

            currentPatch.disturbance_rates[dtype_ilog] +=
                (currentPatch.area - currentPatch.total_canopy_area) * harvest_rate /
                currentPatch.area

            dist_rate_ldist_notharvested +=
                (currentPatch.area - currentPatch.total_canopy_area) * harvest_rate /
                currentPatch.area
        end

        # nocomp mode: prevent producing too-small patches.
        if hlm_use_nocomp[] == itrue &&
           currentPatch.disturbance_rates[dtype_ilog] * currentPatch.area < min_patch_area_forced
            currentPatch.disturbance_rates[dtype_ilog] = 0.0
        end

        # Fraction of the logging disturbance rate that is non-harvested.
        if currentPatch.disturbance_rates[dtype_ilog] > nearzero
            currentPatch.fract_ldist_not_harvested =
                dist_rate_ldist_notharvested / currentPatch.disturbance_rates[dtype_ilog]
        end

        # Fire disturbance rate.
        currentPatch.disturbance_rates[dtype_ifire] = currentPatch.frac_burnt

        # Reduce all rates proportionally if they would exceed the patch area today.
        if (sum(currentPatch.disturbance_rates) +
            sum(@view currentPatch.landuse_transition_rates[1:n_landuse_cats])) >
           max_daily_disturbance_rate
            tempsum = sum(currentPatch.disturbance_rates) +
                sum(@view currentPatch.landuse_transition_rates[1:n_landuse_cats])
            for i_dist in 1:N_DIST_TYPES
                currentPatch.disturbance_rates[i_dist] =
                    max_daily_disturbance_rate * currentPatch.disturbance_rates[i_dist] / tempsum
            end
            for i_dist in 1:n_landuse_cats
                currentPatch.landuse_transition_rates[i_dist] =
                    max_daily_disturbance_rate * currentPatch.landuse_transition_rates[i_dist] / tempsum
            end
        end

        currentPatch = currentPatch.younger
    end

    # Set the secondaryland-exceeds-minimum flag.
    if hlm_use_luh[] == itrue &&
       state_vector[secondaryland] > site_in.min_allowed_landuse_fraction &&
       !site_in.landuse_vector_gt_min[secondaryland]
        site_in.landuse_vector_gt_min[secondaryland] = true
    end

    return nothing
end

# `fates_mortality_disturbance_fraction` is an EDParams scalar (named like the
# Fortran module parameter). Wrap the field read so the call-sites read like
# the Fortran.
fates_mortality_disturbance_fraction() = ed_params().fates_mortality_disturbance_fraction

# ===========================================================================
# CopyPatchMeansTimers
# ===========================================================================

"""
    CopyPatchMeansTimers(dp, rp)

Copy the running means / timers from a donor patch `dp` into a recipient patch
`rp` (the recipient inherits all the donor's running-mean state). Mirrors the
Fortran `CopyPatchMeansTimers`.
"""
function CopyPatchMeansTimers(dp::fates_patch_type, rp::fates_patch_type)
    CopyFromDonor!(rp.tveg24, dp.tveg24)
    CopyFromDonor!(rp.tveg_lpa, dp.tveg_lpa)
    CopyFromDonor!(rp.tveg_longterm, dp.tveg_longterm)

    if ed_params().regeneration_model == TRS_regeneration
        CopyFromDonor!(rp.seedling_layer_par24, dp.seedling_layer_par24)
        CopyFromDonor!(rp.sdlng_mort_par, dp.sdlng_mort_par)
        CopyFromDonor!(rp.sdlng2sap_par, dp.sdlng2sap_par)
        for ipft in 1:numpft[]
            CopyFromDonor!(rp.sdlng_emerg_smp[ipft].p, dp.sdlng_emerg_smp[ipft].p)
            CopyFromDonor!(rp.sdlng_mdd[ipft].p, dp.sdlng_mdd[ipft].p)
        end
    end
    return nothing
end

# ===========================================================================
# GetPseudoPatchAge
# ===========================================================================

"""
    GetPseudoPatchAge(currentPatch) -> Float64

Pseudo-age used to sort patches by (land-use label, nocomp PFT label, age). The
categorical labels are weighted by large numbers and NEGATED so that lower
labels sort "younger" (higher patchno); the continuous age is added normally.
Mirrors the Fortran `GetPseudoPatchAge`.
"""
function GetPseudoPatchAge(currentPatch::fates_patch_type)
    max_actual_age         = 1.0e4   # hard to imagine a patch older than 10,000 years
    max_actual_age_squared = 1.0e8
    return -1.0 * (Float64(currentPatch.land_use_label) * max_actual_age_squared +
                   Float64(currentPatch.nocomp_pft_label) * max_actual_age) +
           currentPatch.age
end

# ===========================================================================
# InsertPatch
# ===========================================================================

"""
    InsertPatch(currentSite, newPatch)

Insert `newPatch` into the site's age-ordered patch linked list at the position
given by its [`GetPseudoPatchAge`](@ref) (land-use label, then nocomp PFT label,
then continuous age). Mirrors the Fortran `InsertPatch`.
"""
function InsertPatch(currentSite::ed_site_type, newPatch::fates_patch_type)
    patch_inserted = false

    if GetPseudoPatchAge(newPatch) <= GetPseudoPatchAge(currentSite.youngest_patch)
        # Insert at the head (youngest).
        newPatch.older   = currentSite.youngest_patch
        newPatch.younger = nothing
        currentSite.youngest_patch.younger = newPatch
        currentSite.youngest_patch = newPatch
        patch_inserted = true
    elseif GetPseudoPatchAge(newPatch) >= GetPseudoPatchAge(currentSite.oldest_patch)
        # Insert at the tail (oldest).
        newPatch.younger = currentSite.oldest_patch
        newPatch.older   = nothing
        currentSite.oldest_patch.older = newPatch
        currentSite.oldest_patch = newPatch
        patch_inserted = true
    else
        # Somewhere within the list: find the first patch with an older pseudo-age.
        currentPatch = currentSite.youngest_patch
        while currentPatch !== nothing && !patch_inserted
            if GetPseudoPatchAge(newPatch) < GetPseudoPatchAge(currentPatch)
                newPatch.older   = currentPatch
                newPatch.younger = currentPatch.younger
                currentPatch.younger.older = newPatch
                currentPatch.younger = newPatch
                patch_inserted = true
            end
            currentPatch = currentPatch.older
        end
    end

    if !patch_inserted
        fates_endrun("something has gone wrong in the patch insertion (no place found)")
    end
    return nothing
end

# ===========================================================================
# set_patchno
# ===========================================================================

"""
    set_patchno(currentSite)

Number the patches oldest -> youngest (1, 2, ...). Under fixed-biogeography +
nocomp, the bareground patch gets `patchno = 0` and is excluded from the
veg-patch numbering. Mirrors the Fortran `set_patchno`.
"""
function set_patchno(currentSite::ed_site_type)
    patchno = 1
    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        currentPatch.patchno = patchno
        patchno += 1
        currentPatch = currentPatch.younger
    end

    if hlm_use_fixed_biogeog[] == itrue && hlm_use_nocomp[] == itrue
        patchno = 1
        currentPatch = currentSite.oldest_patch
        while currentPatch !== nothing
            if currentPatch.nocomp_pft_label == nocomp_bareground
                currentPatch.patchno = 0
            else
                currentPatch.patchno = patchno
                patchno += 1
            end
            currentPatch = currentPatch.younger
        end
    end
    return nothing
end

# ===========================================================================
# check_patch_area
# ===========================================================================

"""
    check_patch_area(currentSite)

Verify the patch areas sum to the site area ([`area`](@ref)). A tiny precision
drift is folded onto the largest patch (with the corresponding seed+litter mass
gain tallied into `patch_resize_err`); a drift beyond `area_error_fail` is fatal.
Mirrors the Fortran `check_patch_area`.
"""
function check_patch_area(currentSite::ed_site_type)
    area_error_fail = 1.0e-6

    areatot = 0.0
    largest_area = 0.0
    largestPatch = nothing
    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        areatot += currentPatch.area
        if currentPatch.area > largest_area
            largestPatch = currentPatch
            largest_area = currentPatch.area
        end
        currentPatch = currentPatch.younger
    end

    if abs(areatot - area) > nearzero
        if abs(areatot - area) > area_error_fail
            fates_endrun("Patch areas do not sum to 10000 within tolerance (total $(areatot), error $(areatot - area))")
        end

        for el in 1:num_elements[]
            # Total mass on the largest patch for its current area [kg].
            _, seed_stock, litter_stock = PatchMassStock(largestPatch, el)

            # Scale the total mass by the added area.
            mass_gain = (seed_stock + litter_stock) * (area - areatot) / largestPatch.area
            currentSite.mass_balance[el].patch_resize_err += mass_gain
        end

        largestPatch.area += (area - areatot)
    end
    return nothing
end

# ===========================================================================
# patch_pft_size_profile
# ===========================================================================

"""
    patch_pft_size_profile(currentPatch)

Build the patch's binned (pft x dbh) aboveground structural-biomass density
profile `pft_agb_profile` used by patch fusion. Mirrors the Fortran
`patch_pft_size_profile`.
"""
function patch_pft_size_profile(currentPatch::fates_patch_type)
    gigantictrees = 1.0e8

    fill!(currentPatch.pft_agb_profile, 0.0)

    mind = Vector{Float64}(undef, N_DBH_BINS)
    maxd = Vector{Float64}(undef, N_DBH_BINS)
    for j in 1:N_DBH_BINS
        if j == N_DBH_BINS
            mind[j] = patchfusion_dbhbin_loweredges[j]
            maxd[j] = gigantictrees
        else
            mind[j] = patchfusion_dbhbin_loweredges[j]
            maxd[j] = patchfusion_dbhbin_loweredges[j+1]
        end
    end

    currentCohort = currentPatch.shortest
    while currentCohort !== nothing
        for j in 1:N_DBH_BINS
            if currentCohort.dbh > mind[j] && currentCohort.dbh <= maxd[j]
                currentPatch.pft_agb_profile[currentCohort.pft, j] +=
                    GetState(currentCohort.prt, struct_organ, carbon12_element) *
                    currentCohort.n / currentPatch.area
            end
        end
        currentCohort = currentCohort.taller
    end
    return nothing
end

# ===========================================================================
# countPatches
# ===========================================================================

"""
    countPatches(nsites, sites) -> Int

Count the total number of patches across all `sites`. Mirrors the Fortran
`countPatches`.
"""
function countPatches(nsites::Integer, sites::AbstractVector{ed_site_type})
    totNumPatches = 0
    for s in 1:nsites
        currentPatch = sites[s].oldest_patch
        while currentPatch !== nothing
            totNumPatches += 1
            currentPatch = currentPatch.younger
        end
    end
    return totNumPatches
end

# ===========================================================================
# DistributeSeeds
# ===========================================================================

"""
    DistributeSeeds(currentSite, seed_mass, el, pft)

Distribute a seed mass [kg] over the seed pools of every patch on the site
(area-normalized by the site area). With `homogenize_seed_pfts`, the mass is
split evenly across all PFTs. Mirrors the Fortran `DistributeSeeds`.
"""
function DistributeSeeds(currentSite::ed_site_type, seed_mass::Real, el::Integer,
                         pft::Integer)
    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        litt = currentPatch.litter[el]
        if homogenize_seed_pfts
            litt.seed .+= seed_mass / (area * Float64(numpft[]))
        else
            litt.seed[pft] += seed_mass / area
        end
        currentPatch = currentPatch.younger
    end
    return nothing
end

# ===========================================================================
# TransLitterNewPatch
# ===========================================================================

"""
    TransLitterNewPatch(currentSite, currentPatch, newPatch, patch_site_areadis)

Transfer the existing litter (and fragmentation/seed-decay flux diagnostics)
from a donor patch into a newly created patch, conserving mass. Partially burned
litter goes to the atmosphere (`fuel.frac_burnt` must be zeroed by the caller
when the disturbance is not fire). The remaining donatable mass is split between
the new patch and the donor by `existing_litt_localization`. Mirrors the Fortran
`TransLitterNewPatch`.
"""
function TransLitterNewPatch(currentSite::ed_site_type, currentPatch::fates_patch_type,
                             newPatch::fates_patch_type, patch_site_areadis::Real)
    nlevsoil = currentSite.nlevsoil

    for el in 1:num_elements[]
        site_mass = currentSite.mass_balance[el]
        curr_litt = currentPatch.litter[el]
        new_litt  = newPatch.litter[el]

        # Fragmentation litter flux rates (diagnostic).
        for c in 1:ncwd
            new_litt.ag_cwd_frag[c] += curr_litt.ag_cwd_frag[c] * patch_site_areadis / newPatch.area
            for sl in 1:nlevsoil
                new_litt.bg_cwd_frag[c, sl] += curr_litt.bg_cwd_frag[c, sl] * patch_site_areadis / newPatch.area
            end
        end

        for dcmpy in 1:ndcmpy
            new_litt.leaf_fines_frag[dcmpy] += curr_litt.leaf_fines_frag[dcmpy] * patch_site_areadis / newPatch.area
            for sl in 1:nlevsoil
                new_litt.root_fines_frag[dcmpy, sl] += curr_litt.root_fines_frag[dcmpy, sl] * patch_site_areadis / newPatch.area
            end
        end

        for pft in 1:numpft[]
            new_litt.seed_decay[pft]      += curr_litt.seed_decay[pft]      * patch_site_areadis / newPatch.area
            new_litt.seed_germ_decay[pft] += curr_litt.seed_germ_decay[pft] * patch_site_areadis / newPatch.area
        end

        # Distribute the existing litter: some burns, some goes to the new patch,
        # some is retained by the donor.
        remainder_area = currentPatch.area - patch_site_areadis

        retain_frac = (1.0 - existing_litt_localization) *
            remainder_area / (newPatch.area + remainder_area)
        # donate_frac = 1.0 - retain_frac   # (computed implicitly below)

        if remainder_area > rsnbl_math_prec
            retain_m2 = retain_frac / remainder_area
            donate_m2 = (1.0 - retain_frac) / newPatch.area
        else
            retain_m2 = 0.0
            donate_m2 = 1.0 / newPatch.area
        end

        for c in 1:ncwd
            # Above ground CWD (a fraction may burn).
            donatable_mass = curr_litt.ag_cwd[c] * patch_site_areadis *
                (1.0 - currentPatch.fuel.frac_burnt[c])
            burned_mass = curr_litt.ag_cwd[c] * patch_site_areadis *
                currentPatch.fuel.frac_burnt[c]

            new_litt.ag_cwd[c]  += donatable_mass * donate_m2
            curr_litt.ag_cwd[c] += donatable_mass * retain_m2
            site_mass.burn_flux_to_atm += burned_mass

            # Below ground CWD (none burns).
            for sl in 1:nlevsoil
                donatable_mass = curr_litt.bg_cwd[c, sl] * patch_site_areadis
                new_litt.bg_cwd[c, sl]  += donatable_mass * donate_m2
                curr_litt.bg_cwd[c, sl] += donatable_mass * retain_m2
            end
        end

        for dcmpy in 1:ndcmpy
            # Leaf fines (a fraction may burn).
            donatable_mass = curr_litt.leaf_fines[dcmpy] * patch_site_areadis *
                (1.0 - currentPatch.fuel.frac_burnt[dead_leaves(fuel_classes)])
            burned_mass = curr_litt.leaf_fines[dcmpy] * patch_site_areadis *
                currentPatch.fuel.frac_burnt[dead_leaves(fuel_classes)]

            new_litt.leaf_fines[dcmpy]  += donatable_mass * donate_m2
            curr_litt.leaf_fines[dcmpy] += donatable_mass * retain_m2
            site_mass.burn_flux_to_atm += burned_mass

            # Root fines (none burns).
            for sl in 1:nlevsoil
                donatable_mass = curr_litt.root_fines[dcmpy, sl] * patch_site_areadis
                new_litt.root_fines[dcmpy, sl]  += donatable_mass * donate_m2
                curr_litt.root_fines[dcmpy, sl] += donatable_mass * retain_m2
            end
        end

        for pft in 1:numpft[]
            # Seeds (currently never burned).
            donatable_mass = curr_litt.seed[pft] * patch_site_areadis
            new_litt.seed[pft]  += donatable_mass * donate_m2
            curr_litt.seed[pft] += donatable_mass * retain_m2

            donatable_mass = curr_litt.seed_germ[pft] * patch_site_areadis
            new_litt.seed_germ[pft]  += donatable_mass * donate_m2
            curr_litt.seed_germ[pft] += donatable_mass * retain_m2
        end
    end
    return nothing
end

# ===========================================================================
# fire_litter_fluxes
# ===========================================================================

"""
    fire_litter_fluxes(currentSite, currentPatch, newPatch, patch_site_areadis, bc_in)

Move the biomass of fire-killed plants (and the partially burned crown of
survivors) from the donor patch into the new patch's litter / CWD pools, with the
burned fraction sent to the atmosphere. Early-returns if there was no fire in the
donor patch. Mirrors the Fortran `fire_litter_fluxes`. (NB: the donor cohort
number densities are scaled down for area AFTER this routine, by the caller.)
"""
function fire_litter_fluxes(currentSite::ed_site_type, currentPatch::fates_patch_type,
                            newPatch::fates_patch_type, patch_site_areadis::Real,
                            bc_in::bc_in_type)
    # Only do this if there was a fire in this actual patch.
    currentPatch.fire == ifalse && return nothing

    # If plant hydraulics are on, account for water leaving the plant-soil mass
    # balance through the fire-killed trees.
    if hlm_use_planthydro[] == itrue
        currentCohort = currentPatch.shortest
        while currentCohort !== nothing
            num_dead_trees = currentCohort.fire_mort * currentCohort.n *
                patch_site_areadis / currentPatch.area
            AccumulateMortalityWaterStorage(currentSite, currentCohort, num_dead_trees)
            currentCohort = currentCohort.taller
        end
    end

    nlevsoil = currentSite.nlevsoil
    SF_val_CWD_frac = sf_params().SF_val_CWD_frac
    SF_val_CWD_frac_adj = zeros(Float64, ncwd)

    remainder_area = currentPatch.area - patch_site_areadis
    retain_frac = (1.0 - burn_localization) *
        remainder_area / (newPatch.area + remainder_area)

    if remainder_area > rsnbl_math_prec
        retain_m2 = retain_frac / remainder_area
        donate_m2 = (1.0 - retain_frac) / newPatch.area
    else
        retain_m2 = 0.0
        donate_m2 = 1.0 / newPatch.area
    end

    for el in 1:num_elements[]
        element_id   = element_list[el]
        site_mass    = currentSite.mass_balance[el]
        elflux_diags = currentSite.flux_diags.elem[el]
        curr_litt    = currentPatch.litter[el]
        new_litt     = newPatch.litter[el]

        currentCohort = currentPatch.shortest
        while currentCohort !== nothing
            pft = currentCohort.pft

            fnrt_m  = GetState(currentCohort.prt, fnrt_organ, element_id)
            store_m = GetState(currentCohort.prt, store_organ, element_id)
            repro_m = GetState(currentCohort.prt, repro_organ, element_id)

            if prt_params.woody[currentCohort.pft] == itrue
                leaf_m   = GetState(currentCohort.prt, leaf_organ, element_id)
                sapw_m   = GetState(currentCohort.prt, sapw_organ, element_id)
                struct_m = GetState(currentCohort.prt, struct_organ, element_id)
            else
                leaf_m   = GetState(currentCohort.prt, leaf_organ, element_id) +
                           GetState(currentCohort.prt, sapw_organ, element_id) +
                           GetState(currentCohort.prt, struct_organ, element_id)
                sapw_m   = 0.0
                struct_m = 0.0
            end

            # Absolute number of dead trees being transferred with the donated area.
            num_dead_trees = currentCohort.fire_mort * currentCohort.n *
                patch_site_areadis / currentPatch.area

            # Dead trees -> leaf litter (and the burned part -> atmosphere).
            donatable_mass = num_dead_trees * (leaf_m + repro_m) *
                (1.0 - currentCohort.fraction_crown_burned)
            burned_mass = num_dead_trees * (leaf_m + repro_m) *
                currentCohort.fraction_crown_burned

            for dcmpy in 1:ndcmpy
                dcmpy_frac = GetDecompyFrac(pft, leaf_organ, dcmpy)
                new_litt.leaf_fines[dcmpy]  += donatable_mass * donate_m2 * dcmpy_frac
                curr_litt.leaf_fines[dcmpy] += donatable_mass * retain_m2 * dcmpy_frac
            end
            site_mass.burn_flux_to_atm += burned_mass

            set_root_fraction(currentSite.rootfrac_scr, pft, currentSite.zi_soil;
                              max_nlevroot = bc_in.max_rooting_depth_index_col)

            # Dead trees -> root litter (no root burn flux).
            for dcmpy in 1:ndcmpy
                dcmpy_frac = GetDecompyFrac(pft, fnrt_organ, dcmpy)
                for sl in 1:nlevsoil
                    donatable_mass = num_dead_trees * (fnrt_m + store_m) * currentSite.rootfrac_scr[sl]
                    new_litt.root_fines[dcmpy, sl]  += donatable_mass * donate_m2 * dcmpy_frac
                    curr_litt.root_fines[dcmpy, sl] += donatable_mass * retain_m2 * dcmpy_frac
                end
            end

            elflux_diags.surf_fine_litter_input[pft] +=
                num_dead_trees * (leaf_m + repro_m) * (1.0 - currentCohort.fraction_crown_burned)
            elflux_diags.root_litter_input[pft] += (fnrt_m + store_m) * num_dead_trees

            # Coarse root biomass per tree -> below ground CWD.
            bcroot = (sapw_m + struct_m) * (1.0 - prt_params.allom_agb_frac[pft])
            adjust_SF_CWD_frac(currentCohort.dbh, ncwd, SF_val_CWD_frac, SF_val_CWD_frac_adj)

            for c in 1:ncwd
                for sl in 1:nlevsoil
                    donatable_mass = num_dead_trees * SF_val_CWD_frac_adj[c] *
                        bcroot * currentSite.rootfrac_scr[sl]
                    new_litt.bg_cwd[c, sl]  += donatable_mass * donate_m2
                    curr_litt.bg_cwd[c, sl] += donatable_mass * retain_m2
                    elflux_diags.cwd_bg_input[c] += donatable_mass
                end
            end

            # Stem biomass per tree -> above ground CWD (twig/sbranch may burn).
            bstem = (sapw_m + struct_m) * prt_params.allom_agb_frac[pft]
            for c in 1:ncwd
                donatable_mass = num_dead_trees * SF_val_CWD_frac_adj[c] * bstem
                if c == 1 || c == 2
                    donatable_mass *= (1.0 - currentCohort.fraction_crown_burned)
                    burned_mass = num_dead_trees * SF_val_CWD_frac_adj[c] * bstem *
                        currentCohort.fraction_crown_burned
                    site_mass.burn_flux_to_atm += burned_mass
                end
                new_litt.ag_cwd[c]  += donatable_mass * donate_m2
                curr_litt.ag_cwd[c] += donatable_mass * retain_m2
                elflux_diags.cwd_ag_input[c] += donatable_mass
            end

            currentCohort = currentCohort.taller
        end
    end
    return nothing
end

# ===========================================================================
# mortality_litter_fluxes
# ===========================================================================

"""
    mortality_litter_fluxes(currentSite, currentPatch, newPatch, patch_site_areadis, bc_in)

Move the biomass of plants killed by disturbance-associated treefall mortality
(canopy gap formers + understory impact kills) into the new patch's CWD / fine
litter pools, splitting between new and donor patch by `treefall_localization`.
Mirrors the Fortran `mortality_litter_fluxes`.

Upstream quirk preserved: the `DistributeSeeds(...)` call for the
storage->reproduction shunt is commented out in the Fortran ("BREAKING MASS
CONSERVATION RIGHT NOW"); the seed mass goes straight into the seed litter pool.
"""
function mortality_litter_fluxes(currentSite::ed_site_type, currentPatch::fates_patch_type,
                                 newPatch::fates_patch_type, patch_site_areadis::Real,
                                 bc_in::bc_in_type)
    nlevsoil = currentSite.nlevsoil
    SF_val_CWD_frac = sf_params().SF_val_CWD_frac
    SF_val_CWD_frac_adj = zeros(Float64, ncwd)
    ev = EDPftvarcon_inst[]

    remainder_area = currentPatch.area - patch_site_areadis
    retain_frac = (1.0 - treefall_localization) *
        remainder_area / (newPatch.area + remainder_area)

    if remainder_area > rsnbl_math_prec
        retain_m2 = retain_frac / remainder_area
        donate_m2 = (1.0 - retain_frac) / newPatch.area
    else
        retain_m2 = 0.0
        donate_m2 = 1.0 / newPatch.area
    end

    for el in 1:num_elements[]
        element_id   = element_list[el]
        site_mass    = currentSite.mass_balance[el]
        elflux_diags = currentSite.flux_diags.elem[el]
        curr_litt    = currentPatch.litter[el]
        new_litt     = newPatch.litter[el]

        currentCohort = currentPatch.shortest
        while currentCohort !== nothing
            pft = currentCohort.pft

            fnrt_m  = GetState(currentCohort.prt, fnrt_organ, element_id)
            store_m = GetState(currentCohort.prt, store_organ, element_id)
            repro_m = GetState(currentCohort.prt, repro_organ, element_id)

            if prt_params.woody[currentCohort.pft] == itrue
                leaf_m   = GetState(currentCohort.prt, leaf_organ, element_id)
                sapw_m   = GetState(currentCohort.prt, sapw_organ, element_id)
                struct_m = GetState(currentCohort.prt, struct_organ, element_id)
            else
                leaf_m   = GetState(currentCohort.prt, leaf_organ, element_id) +
                           GetState(currentCohort.prt, sapw_organ, element_id) +
                           GetState(currentCohort.prt, struct_organ, element_id)
                sapw_m   = 0.0
                struct_m = 0.0
            end

            if currentCohort.canopy_layer == 1
                # Upper-canopy trees: dead from disturbance-generating mortality.
                num_dead = currentCohort.n * min(1.0,
                    currentCohort.dmort * hlm_freq_day[] * fates_mortality_disturbance_fraction())
            elseif prt_params.woody[pft] == itrue
                # Understory trees: dead from survivorship + disturbed area.
                num_dead = ed_params().ED_val_understorey_death * currentCohort.n *
                    (patch_site_areadis / currentPatch.area)
            else
                # Understory grasses are not killed by tree-fall disturbance.
                num_dead = 0.0
            end

            # Update water balance by removing dead plant water, but only once
            # (use the carbon element id).
            if element_id == carbon12_element && hlm_use_planthydro[] == itrue
                AccumulateMortalityWaterStorage(currentSite, currentCohort, num_dead)
            end

            # Leaves (+ seeds) of dying trees -> leaf litter.
            for dcmpy in 1:ndcmpy
                dcmpy_frac = GetDecompyFrac(pft, leaf_organ, dcmpy)
                new_litt.leaf_fines[dcmpy]  += num_dead * (leaf_m + repro_m) * donate_m2 * dcmpy_frac
                curr_litt.leaf_fines[dcmpy] += num_dead * (leaf_m + repro_m) * retain_m2 * dcmpy_frac
            end

            ag_wood = num_dead * (struct_m + sapw_m) * prt_params.allom_agb_frac[pft]
            bg_wood = num_dead * (struct_m + sapw_m) * (1.0 - prt_params.allom_agb_frac[pft])

            set_root_fraction(currentSite.rootfrac_scr, pft, currentSite.zi_soil;
                              max_nlevroot = bc_in.max_rooting_depth_index_col)
            adjust_SF_CWD_frac(currentCohort.dbh, ncwd, SF_val_CWD_frac, SF_val_CWD_frac_adj)

            for c in 1:ncwd
                new_litt.ag_cwd[c]  += ag_wood * SF_val_CWD_frac_adj[c] * donate_m2
                curr_litt.ag_cwd[c] += ag_wood * SF_val_CWD_frac_adj[c] * retain_m2
                for sl in 1:nlevsoil
                    new_litt.bg_cwd[c, sl]  += bg_wood * currentSite.rootfrac_scr[sl] *
                        SF_val_CWD_frac_adj[c] * donate_m2
                    curr_litt.bg_cwd[c, sl] += bg_wood * currentSite.rootfrac_scr[sl] *
                        SF_val_CWD_frac_adj[c] * retain_m2
                end
            end

            # Fine roots (+ the non-reproductive storage) -> below ground litter.
            for dcmpy in 1:ndcmpy
                dcmpy_frac = GetDecompyFrac(pft, fnrt_organ, dcmpy)
                for sl in 1:nlevsoil
                    new_litt.root_fines[dcmpy, sl] += num_dead * currentSite.rootfrac_scr[sl] *
                        (fnrt_m + store_m * (1.0 - ev.allom_frbstor_repro[pft])) * donate_m2 * dcmpy_frac
                    curr_litt.root_fines[dcmpy, sl] += num_dead * currentSite.rootfrac_scr[sl] *
                        (fnrt_m + store_m * (1.0 - ev.allom_frbstor_repro[pft])) * retain_m2 * dcmpy_frac
                end
            end

            # Storage shunted to reproduction on death -> seed pool.
            seed_mass = num_dead * store_m * ev.allom_frbstor_repro[pft]
            # Upstream quirk: DistributeSeeds is commented out (mass conservation).
            # DistributeSeeds(currentSite, seed_mass, el, pft)
            new_litt.seed[pft]  += seed_mass * donate_m2
            curr_litt.seed[pft] += seed_mass * retain_m2

            for c in 1:ncwd
                elflux_diags.cwd_ag_input[c] += SF_val_CWD_frac_adj[c] * ag_wood
                elflux_diags.cwd_bg_input[c] += SF_val_CWD_frac_adj[c] * bg_wood
            end
            elflux_diags.surf_fine_litter_input[pft] += num_dead * (leaf_m + repro_m)
            elflux_diags.root_litter_input[pft] +=
                num_dead * (fnrt_m + store_m * (1.0 - ev.allom_frbstor_repro[pft]))

            currentCohort = currentCohort.taller
        end
    end
    return nothing
end

# ===========================================================================
# landusechange_litter_fluxes
# ===========================================================================

"""
    landusechange_litter_fluxes(currentSite, currentPatch, newPatch,
                                patch_site_areadis, bc_in, clearing_matrix_element)

Move felled-plant biomass during a land-use change into the new patch's CWD /
fine litter pools, with a fraction exported as wood products and a fraction
burned. Does nothing unless `clearing_matrix_element` is true (i.e. the LU
transition clears vegetation). Mirrors the Fortran `landusechange_litter_fluxes`.
"""
function landusechange_litter_fluxes(currentSite::ed_site_type, currentPatch::fates_patch_type,
                                     newPatch::fates_patch_type, patch_site_areadis::Real,
                                     bc_in::bc_in_type, clearing_matrix_element::Bool)
    clearing_matrix_element || return nothing

    nlevsoil = currentSite.nlevsoil
    SF_val_CWD_frac = sf_params().SF_val_CWD_frac
    ev = EDPftvarcon_inst[]

    # If plant hydraulics are on, account for water leaving the plant-soil mass
    # balance through the cleared (dead) trees.
    if hlm_use_planthydro[] == itrue
        currentCohort = currentPatch.shortest
        while currentCohort !== nothing
            num_dead_trees = currentCohort.n * patch_site_areadis / currentPatch.area
            AccumulateMortalityWaterStorage(currentSite, currentCohort, num_dead_trees)
            currentCohort = currentCohort.taller
        end
    end

    remainder_area = currentPatch.area - patch_site_areadis
    retain_frac = (1.0 - landusechange_localization) *
        remainder_area / (newPatch.area + remainder_area)

    if remainder_area > rsnbl_math_prec
        retain_m2 = retain_frac / remainder_area
        donate_m2 = (1.0 - retain_frac) / newPatch.area
    else
        retain_m2 = 0.0
        donate_m2 = 1.0 / newPatch.area
    end

    for el in 1:num_elements[]
        trunk_product_site = 0.0

        element_id   = element_list[el]
        site_mass    = currentSite.mass_balance[el]
        elflux_diags = currentSite.flux_diags.elem[el]
        curr_litt    = currentPatch.litter[el]
        new_litt     = newPatch.litter[el]

        currentCohort = currentPatch.shortest
        while currentCohort !== nothing
            pft = currentCohort.pft

            fnrt_m  = GetState(currentCohort.prt, fnrt_organ, element_id)
            store_m = GetState(currentCohort.prt, store_organ, element_id)
            repro_m = GetState(currentCohort.prt, repro_organ, element_id)

            if prt_params.woody[currentCohort.pft] == itrue
                leaf_m   = GetState(currentCohort.prt, leaf_organ, element_id)
                sapw_m   = GetState(currentCohort.prt, sapw_organ, element_id)
                struct_m = GetState(currentCohort.prt, struct_organ, element_id)
            else
                leaf_m   = GetState(currentCohort.prt, leaf_organ, element_id) +
                           GetState(currentCohort.prt, sapw_organ, element_id) +
                           GetState(currentCohort.prt, struct_organ, element_id)
                sapw_m   = 0.0
                struct_m = 0.0
            end

            num_dead_trees = currentCohort.n * patch_site_areadis / currentPatch.area

            # Leaves (a fraction may burn).
            donatable_mass = num_dead_trees * (leaf_m + repro_m) *
                (1.0 - ev.landusechange_frac_burned[pft])
            burned_mass = num_dead_trees * (leaf_m + repro_m) * ev.landusechange_frac_burned[pft]

            for dcmpy in 1:ndcmpy
                dcmpy_frac = GetDecompyFrac(pft, leaf_organ, dcmpy)
                new_litt.leaf_fines[dcmpy]  += donatable_mass * donate_m2 * dcmpy_frac
                curr_litt.leaf_fines[dcmpy] += donatable_mass * retain_m2 * dcmpy_frac
            end
            site_mass.burn_flux_to_atm += burned_mass

            set_root_fraction(currentSite.rootfrac_scr, pft, currentSite.zi_soil;
                              max_nlevroot = bc_in.max_rooting_depth_index_col)

            # Roots (no burn flux).
            for dcmpy in 1:ndcmpy
                dcmpy_frac = GetDecompyFrac(pft, fnrt_organ, dcmpy)
                for sl in 1:nlevsoil
                    donatable_mass = num_dead_trees * (fnrt_m + store_m) * currentSite.rootfrac_scr[sl]
                    new_litt.root_fines[dcmpy, sl]  += donatable_mass * donate_m2 * dcmpy_frac
                    curr_litt.root_fines[dcmpy, sl] += donatable_mass * retain_m2 * dcmpy_frac
                end
            end

            elflux_diags.surf_fine_litter_input[pft] +=
                num_dead_trees * (leaf_m + repro_m) * (1.0 - ev.landusechange_frac_burned[pft])
            elflux_diags.root_litter_input[pft] += (fnrt_m + store_m) * num_dead_trees

            bcroot = (sapw_m + struct_m) * (1.0 - prt_params.allom_agb_frac[pft])
            for c in 1:ncwd
                for sl in 1:nlevsoil
                    donatable_mass = num_dead_trees * SF_val_CWD_frac[c] *
                        bcroot * currentSite.rootfrac_scr[sl]
                    new_litt.bg_cwd[c, sl]  += donatable_mass * donate_m2
                    curr_litt.bg_cwd[c, sl] += donatable_mass * retain_m2
                    elflux_diags.cwd_bg_input[c] += donatable_mass
                end
            end

            bstem = (sapw_m + struct_m) * prt_params.allom_agb_frac[pft]
            for c in 1:ncwd
                donatable_mass = num_dead_trees * SF_val_CWD_frac[c] * bstem
                if c == 1 || c == 2   # these pools can burn
                    donatable_mass *= (1.0 - ev.landusechange_frac_burned[pft])
                    burned_mass = num_dead_trees * SF_val_CWD_frac[c] * bstem *
                        ev.landusechange_frac_burned[pft]
                    site_mass.burn_flux_to_atm += burned_mass
                else  # can also be exported as timber products
                    donatable_mass *= (1.0 - ev.landusechange_frac_exported[pft]) *
                        (1.0 - ev.landusechange_frac_burned[pft])
                    burned_mass = num_dead_trees * SF_val_CWD_frac[c] * bstem *
                        (1.0 - ev.landusechange_frac_exported[pft]) * ev.landusechange_frac_burned[pft]
                    woodproduct_mass = num_dead_trees * SF_val_CWD_frac[c] * bstem *
                        ev.landusechange_frac_exported[pft]
                    site_mass.burn_flux_to_atm += burned_mass
                    trunk_product_site += woodproduct_mass
                    elflux_diags.exported_harvest += woodproduct_mass * area_inv
                    site_mass.wood_product_landusechange[pft] += woodproduct_mass
                end
                new_litt.ag_cwd[c]  += donatable_mass * donate_m2
                curr_litt.ag_cwd[c] += donatable_mass * retain_m2
                elflux_diags.cwd_ag_input[c] += donatable_mass
            end

            currentCohort = currentCohort.taller
        end

        # Carbon exported from the site as wood product.
        if element_id == carbon12_element
            currentSite.resources_management.trunk_product_site += trunk_product_site
        end
    end
    return nothing
end

# ===========================================================================
# split_patch
# ===========================================================================

"""
    split_patch(currentSite, currentPatch, new_patch, fraction_to_keep, area_to_remove=nothing)

Split a patch into two identical patches differing only in area: `currentPatch`
keeps `fraction_to_keep` of its area (and number density), `new_patch` receives
the rest. Cohorts are split across the two by number density (`nc.n = n * (1 -
fraction_to_keep)`). Litter is moved by [`TransLitterNewPatch`](@ref). Mirrors the
Fortran `split_patch`.
"""
function split_patch(currentSite::ed_site_type, currentPatch::fates_patch_type,
                     new_patch::fates_patch_type, fraction_to_keep::Real,
                     area_to_remove=nothing)
    temp_area = area_to_remove === nothing ?
        currentPatch.area - (currentPatch.area * fraction_to_keep) :
        area_to_remove

    Create!(new_patch, 0.0, temp_area, currentPatch.land_use_label,
            currentPatch.nocomp_pft_label, num_swb, numpft[], currentSite.nlevsoil,
            hlm_current_tod[], ed_params().regeneration_model)

    for el in 1:num_elements[]
        InitConditions!(new_patch.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
    new_patch.tallest  = nothing
    new_patch.shortest = nothing

    CopyPatchMeansTimers(currentPatch, new_patch)
    TransLitterNewPatch(currentSite, currentPatch, new_patch, temp_area)

    fill!(currentPatch.fuel.frac_burnt, 0.0)

    currentCohort = currentPatch.shortest
    while currentCohort !== nothing
        nc = fates_cohort_type()
        if hlm_use_planthydro[] == itrue
            InitHydrCohort(currentSite, nc)
        end
        nc.prt = InitPRTObject!()
        InitPRTBoundaryConditions(nc)
        ZeroValues(nc)
        Copy(currentCohort, nc)

        # Number of members in the new patch (the rest stays in the donor).
        nc.n = currentCohort.n * (1.0 - fraction_to_keep)
        currentCohort.n = currentCohort.n * fraction_to_keep

        _insert_new_cohort!(new_patch, nc)

        currentCohort = currentCohort.taller
    end

    sort_cohorts(currentPatch)

    # Update area of donor patch.
    currentPatch.area = currentPatch.area - temp_area
    return nothing
end

# Splice a fully-formed cohort `nc` into `patch`'s height-ordered list, carrying
# the head/tail through `insert_cohort` (the storebigcohort/storesmallcohort
# Ref dance the Fortran does with pass-by-pointer). Mirrors the inline block
# repeated in spawn_patches / split_patch / fuse_2_patches.
function _insert_new_cohort!(patch::fates_patch_type, nc::fates_cohort_type)
    storebigcohort   = Ref{Union{fates_cohort_type,Nothing}}(patch.tallest)
    storesmallcohort = Ref{Union{fates_cohort_type,Nothing}}(patch.shortest)

    if patch.tallest !== nothing
        tnull = 0
    else
        tnull = 1
        patch.tallest = nc
        nc.taller = nothing
    end
    if patch.shortest !== nothing
        snull = 0
    else
        snull = 1
        patch.shortest = nc
        nc.shorter = nothing
    end

    insert_cohort(patch, nc, patch.tallest, patch.shortest, tnull, snull,
                  storebigcohort, storesmallcohort)

    patch.tallest  = storebigcohort[]
    patch.shortest = storesmallcohort[]
    return nothing
end

# ===========================================================================
# spawn_patches
# ===========================================================================

"""
    spawn_patches(currentSite, bc_in)

Resolve the per-patch disturbance computed by [`disturbance_rates`](@ref): for
each (nocomp-PFT, disturbance type, donor land-use label, receiver land-use
label) combination, create a new patch absorbing all the disturbed area from
matching donor patches, transfer in the donor litter and the killed/disturbed
plant litter, copy surviving cohorts (with disturbance-type-specific
survivorship) into the new patch, splice it into the list, then enforce patch
area + numbering. Mirrors the Fortran `spawn_patches`.

The nocomp+LUH patch-relabelling block at the end is ported but inert unless both
flags are on (they default off).
"""
function spawn_patches(currentSite::ed_site_type, bc_in::bc_in_type)
    if hlm_use_nocomp[] == itrue
        min_nocomp_pft = 0
        max_nocomp_pft = numpft[]
    else
        min_nocomp_pft = fates_unset_int
        max_nocomp_pft = fates_unset_int
    end

    # Zero the diagnostic disturbance rate fields.
    fill!(currentSite.disturbance_rates, 0.0)

    # Rules for vegetation clearing during land use change.
    clearing_matrix = GetLanduseChangeRules()

    for i_nocomp_pft in min_nocomp_pft:max_nocomp_pft
        for i_disturbance_type in 1:N_DIST_TYPES
            for i_donorpatch_landuse_type in 1:n_landuse_cats

                # Receiver land-use label loop bounds.
                if i_disturbance_type == dtype_ifire || i_disturbance_type == dtype_ifall
                    start_receiver_lulabel = i_donorpatch_landuse_type
                    end_receiver_lulabel   = i_donorpatch_landuse_type
                elseif i_disturbance_type == dtype_ilog
                    start_receiver_lulabel = secondaryland
                    end_receiver_lulabel   = secondaryland
                elseif i_disturbance_type == dtype_ilandusechange
                    start_receiver_lulabel = 1
                    end_receiver_lulabel   = n_landuse_cats
                else
                    fates_endrun("unknown disturbance mode? $(i_disturbance_type)")
                end

                for i_landusechange_receiverpatchlabel in start_receiver_lulabel:end_receiver_lulabel

                    # --- Total disturbed area meeting these criteria ---------
                    site_areadis = 0.0
                    currentPatch = currentSite.youngest_patch
                    while currentPatch !== nothing
                        if (hlm_use_nocomp[] == ifalse || currentPatch.nocomp_pft_label == i_nocomp_pft) &&
                           currentPatch.land_use_label == i_donorpatch_landuse_type

                            if i_disturbance_type != dtype_ilandusechange
                                disturbance_rate = currentPatch.disturbance_rates[i_disturbance_type]
                            else
                                disturbance_rate = currentPatch.landuse_transition_rates[i_landusechange_receiverpatchlabel]
                            end

                            if disturbance_rate > (1.0 + rsnbl_math_prec)
                                fates_endrun("patch disturbance rate > 1 ? $(disturbance_rate)")
                            elseif disturbance_rate > 1.0
                                disturbance_rate = 1.0
                            end

                            if (currentPatch.area * disturbance_rate) > nearzero
                                site_areadis += currentPatch.area * disturbance_rate
                                currentSite.disturbance_rates[i_disturbance_type,
                                    i_donorpatch_landuse_type, i_landusechange_receiverpatchlabel] +=
                                    currentPatch.area * disturbance_rate * area_inv
                            end
                        end
                        currentPatch = currentPatch.older
                    end

                    # It is possible no disturbance area was generated.
                    newPatch = nothing
                    if site_areadis > nearzero
                        newPatch = fates_patch_type()
                        Create!(newPatch, 0.0, site_areadis, i_landusechange_receiverpatchlabel,
                                i_nocomp_pft, num_swb, numpft[], currentSite.nlevsoil,
                                hlm_current_tod[], ed_params().regeneration_model)
                        for el in 1:num_elements[]
                            InitConditions!(newPatch.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                        end
                        newPatch.tallest  = nothing
                        newPatch.shortest = nothing
                    end

                    # --- Loop the donor patches, transfer litter + cohorts ---
                    currentPatch = currentSite.oldest_patch
                    while currentPatch !== nothing
                        if (hlm_use_nocomp[] == ifalse || currentPatch.nocomp_pft_label == i_nocomp_pft) &&
                           currentPatch.land_use_label == i_donorpatch_landuse_type

                            if i_disturbance_type != dtype_ilandusechange
                                disturbance_rate = currentPatch.disturbance_rates[i_disturbance_type]
                            else
                                disturbance_rate = currentPatch.landuse_transition_rates[i_landusechange_receiverpatchlabel]
                            end

                            patch_site_areadis = currentPatch.area * disturbance_rate

                            if patch_site_areadis > nearzero
                                newPatch === nothing && fates_endrun("Patch spawning pointed to an un-allocated patch")

                                # Average in age-since-anthro-disturbance for
                                # non-primary donors under non-anthro disturbance.
                                if currentPatch.land_use_label > primaryland &&
                                   i_disturbance_type < dtype_ilog
                                    newPatch.age_since_anthro_disturbance +=
                                        currentPatch.age_since_anthro_disturbance *
                                        (patch_site_areadis / site_areadis)
                                end

                                # Zero the donor burn fraction unless this is fire.
                                if i_disturbance_type != dtype_ifire
                                    fill!(currentPatch.fuel.frac_burnt, 0.0)
                                end

                                CopyPatchMeansTimers(currentPatch, newPatch)
                                TransLitterNewPatch(currentSite, currentPatch, newPatch, patch_site_areadis)

                                # Killed/disturbed plant litter fluxes.
                                if i_disturbance_type == dtype_ilog
                                    logging_litter_fluxes!(currentSite, currentPatch, newPatch,
                                        patch_site_areadis, bc_in)
                                    if i_donorpatch_landuse_type == primaryland
                                        newPatch.changed_landuse_this_ts = true
                                    end
                                elseif i_disturbance_type == dtype_ifire
                                    fire_litter_fluxes(currentSite, currentPatch, newPatch,
                                        patch_site_areadis, bc_in)
                                elseif i_disturbance_type == dtype_ifall
                                    mortality_litter_fluxes(currentSite, currentPatch, newPatch,
                                        patch_site_areadis, bc_in)
                                elseif i_disturbance_type == dtype_ilandusechange
                                    landusechange_litter_fluxes(currentSite, currentPatch, newPatch,
                                        patch_site_areadis, bc_in,
                                        clearing_matrix[i_donorpatch_landuse_type, i_landusechange_receiverpatchlabel])
                                    newPatch.changed_landuse_this_ts = true
                                else
                                    fates_endrun("unknown disturbance mode? $(i_disturbance_type)")
                                end

                                _spawn_move_cohorts!(currentSite, currentPatch, newPatch,
                                    patch_site_areadis, i_disturbance_type,
                                    i_donorpatch_landuse_type, i_landusechange_receiverpatchlabel,
                                    clearing_matrix, bc_in)

                                sort_cohorts(currentPatch)

                                # Update donor patch area + rescale unresolved rates.
                                oldarea = currentPatch.area
                                currentPatch.area = currentPatch.area - patch_site_areadis

                                if i_disturbance_type < N_DIST_TYPES
                                    for i_dist2 in (i_disturbance_type+1):(N_DIST_TYPES-1)
                                        currentPatch.disturbance_rates[i_dist2] *= oldarea / currentPatch.area
                                    end
                                    for i_dist2 in 1:n_landuse_cats
                                        currentPatch.landuse_transition_rates[i_dist2] *= oldarea / currentPatch.area
                                    end
                                else
                                    for i_dist2 in (i_landusechange_receiverpatchlabel+1):n_landuse_cats
                                        currentPatch.landuse_transition_rates[i_dist2] *= oldarea / currentPatch.area
                                    end
                                end

                                # Cull / fuse / re-sort the donor patch.
                                terminate_cohorts(currentSite, currentPatch, 1, 16, bc_in)
                                fuse_cohorts(currentSite, currentPatch, bc_in)
                                terminate_cohorts(currentSite, currentPatch, 2, 16, bc_in)
                                sort_cohorts(currentPatch)
                            end
                        end
                        currentPatch = currentPatch.younger
                    end

                    # --- Insert the new patch into the linked list ----------
                    if site_areadis > nearzero
                        InsertPatch(currentSite, newPatch)
                        terminate_cohorts(currentSite, newPatch, 1, 17, bc_in)
                        fuse_cohorts(currentSite, newPatch, bc_in)
                        terminate_cohorts(currentSite, newPatch, 2, 17, bc_in)
                        sort_cohorts(newPatch)
                    end

                    check_patch_area(currentSite)
                    set_patchno(currentSite)
                end
            end
        end
    end

    # --- nocomp + LUH: remap nocomp PFT identities of newly disturbed patches
    if hlm_use_nocomp[] == itrue && hlm_use_luh[] == itrue
        _spawn_remap_nocomp_luh!(currentSite)
    else
        currentPatch = currentSite.oldest_patch
        while currentPatch !== nothing
            currentPatch.changed_landuse_this_ts = false
            currentPatch = currentPatch.younger
        end
    end

    # Zero disturbance-rate trackers on all patches.
    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        fill!(currentPatch.disturbance_rates, 0.0)
        currentPatch.fract_ldist_not_harvested = 0.0
        currentPatch = currentPatch.younger
    end

    return nothing
end

# Move all of the donor patch's cohorts into the new patch, applying the
# disturbance-type-specific survivorship to both the new (nc) and the donor
# (currentCohort) copies. Factored out of spawn_patches for readability; the
# logic is a faithful transcription of the Fortran cohortloop.
function _spawn_move_cohorts!(currentSite::ed_site_type, currentPatch::fates_patch_type,
                              newPatch::fates_patch_type, patch_site_areadis::Real,
                              i_disturbance_type::Integer, i_donorpatch_landuse_type::Integer,
                              i_landusechange_receiverpatchlabel::Integer,
                              clearing_matrix::AbstractMatrix{Bool}, bc_in)
    ED_val_understorey_death = ed_params().ED_val_understorey_death
    logging_coll_under_frac  = ed_params().logging_coll_under_frac

    currentCohort = currentPatch.shortest
    while currentCohort !== nothing
        nc = fates_cohort_type()
        if hlm_use_planthydro[] == itrue
            InitHydrCohort(currentSite, nc)
        end
        nc.prt = InitPRTObject!()
        InitPRTBoundaryConditions(nc)
        ZeroValues(nc)
        Copy(currentCohort, nc)

        # New patch probably doesn't have a closed canopy; canopy_structure sorts it.
        nc.canopy_layer = 1
        nc.canopy_layer_yesterday = 1.0

        sapw_c   = GetState(currentCohort.prt, sapw_organ, carbon12_element)
        struct_c = GetState(currentCohort.prt, struct_organ, carbon12_element)
        leaf_c   = GetState(currentCohort.prt, leaf_organ, carbon12_element)
        fnrt_c   = GetState(currentCohort.prt, fnrt_organ, carbon12_element)
        store_c  = GetState(currentCohort.prt, store_organ, carbon12_element)
        total_c  = sapw_c + struct_c + leaf_c + fnrt_c + store_c

        if i_disturbance_type == dtype_ifall
            _spawn_surv_treefall!(currentSite, currentPatch, nc, currentCohort,
                patch_site_areadis, total_c, sapw_c, struct_c, store_c, leaf_c,
                ED_val_understorey_death)
        elseif i_disturbance_type == dtype_ifire
            _spawn_surv_fire!(currentSite, currentPatch, nc, currentCohort,
                patch_site_areadis, total_c, sapw_c, struct_c, store_c, leaf_c)
        elseif i_disturbance_type == dtype_ilog
            _spawn_surv_logging!(currentSite, currentPatch, nc, currentCohort,
                patch_site_areadis, total_c, sapw_c, struct_c, store_c, leaf_c,
                logging_coll_under_frac)
        elseif i_disturbance_type == dtype_ilandusechange
            nc.n = currentCohort.n * patch_site_areadis / currentPatch.area
            currentCohort.n = currentCohort.n * (1.0 - patch_site_areadis / currentPatch.area)
            if clearing_matrix[i_donorpatch_landuse_type, i_landusechange_receiverpatchlabel]
                nc.n = 0.0  # kill everything
            end
        else
            fates_endrun("unknown disturbance mode? $(i_disturbance_type)")
        end

        # If any plants in the new cohort survived, insert it; else drop it.
        if nc.n > 0.0
            _insert_new_cohort!(newPatch, nc)
        else
            FreeMemory(nc)
        end

        currentCohort = currentCohort.taller
    end
    return nothing
end

# Treefall survivorship (dtype_ifall).
function _spawn_surv_treefall!(currentSite, currentPatch, nc, currentCohort,
                               patch_site_areadis, total_c, sapw_c, struct_c, store_c,
                               leaf_c, ED_val_understorey_death)
    if currentCohort.canopy_layer == 1
        # Donor keeps fewer trees (the part of the patch where no trees fell).
        currentCohort.n = currentCohort.n * (1.0 -
            fates_mortality_disturbance_fraction() *
            min(1.0, currentCohort.dmort * hlm_freq_day[]))
        nc.n = 0.0  # kill the gap-forming trees
        # Mortality diagnostics -> NaN (the cohort should disappear).
        nc.cmort = NaN; nc.hmort = NaN; nc.bmort = NaN; nc.frmort = NaN
        nc.smort = NaN; nc.asmort = NaN; nc.dgmort = NaN
        nc.lmort_direct = NaN; nc.lmort_collateral = NaN; nc.lmort_infra = NaN
        nc.l_degrad = NaN
    elseif prt_params.woody[currentCohort.pft] == itrue
        # Understory woody plants: step 1 reduce for area change.
        nc.n = currentCohort.n * patch_site_areadis / currentPatch.area

        currentSite.imort_rate[currentCohort.size_class, currentCohort.pft] +=
            nc.n * ED_val_understorey_death / hlm_freq_day[]
        currentSite.imort_carbonflux[currentCohort.pft] +=
            (nc.n * ED_val_understorey_death / hlm_freq_day[]) *
            total_c * g_per_kg * days_per_sec * years_per_day * ha_per_m2
        currentSite.imort_abg_flux[currentCohort.size_class, currentCohort.pft] +=
            (nc.n * ED_val_understorey_death / hlm_freq_day[]) *
            ((sapw_c + struct_c + store_c) * prt_params.allom_agb_frac[currentCohort.pft] + leaf_c) *
            g_per_kg * days_per_sec * years_per_day * ha_per_m2

        # Step 2 apply understory death survivorship.
        nc.n = nc.n * (1.0 - ED_val_understorey_death)
        _carry_mort_diags!(nc, currentCohort)

        # Donor loses only the area-change part (not mortality).
        currentCohort.n = currentCohort.n * (1.0 - patch_site_areadis / currentPatch.area)
    else
        # Grass is not killed by mortality disturbance; just split it.
        nc.n = currentCohort.n * patch_site_areadis / currentPatch.area
        currentCohort.n = currentCohort.n * (1.0 - patch_site_areadis / currentPatch.area)
        _carry_mort_diags!(nc, currentCohort)
    end
    return nothing
end

# Fire survivorship (dtype_ifire).
function _spawn_surv_fire!(currentSite, currentPatch, nc, currentCohort,
                           patch_site_areadis, total_c, sapw_c, struct_c, store_c, leaf_c)
    # Members in the new patch before fire survivorship.
    nc.n = currentCohort.n * patch_site_areadis / currentPatch.area
    # Donor loses individuals from area shrinking.
    currentCohort.n = currentCohort.n * (1.0 - patch_site_areadis / currentPatch.area)

    levcan = currentCohort.canopy_layer
    if levcan == ican_upper
        currentSite.fmort_rate_canopy[currentCohort.size_class, currentCohort.pft] +=
            nc.n * currentCohort.fire_mort / hlm_freq_day[]
        currentSite.fmort_carbonflux_canopy[currentCohort.pft] +=
            (nc.n * currentCohort.fire_mort) * total_c * g_per_kg * days_per_sec * ha_per_m2
    else
        currentSite.fmort_rate_ustory[currentCohort.size_class, currentCohort.pft] +=
            nc.n * currentCohort.fire_mort / hlm_freq_day[]
        currentSite.fmort_carbonflux_ustory[currentCohort.pft] +=
            (nc.n * currentCohort.fire_mort) * total_c * g_per_kg * days_per_sec * ha_per_m2
    end

    currentSite.fmort_abg_flux[currentCohort.size_class, currentCohort.pft] +=
        (nc.n * currentCohort.fire_mort) *
        ((sapw_c + struct_c + store_c) * prt_params.allom_agb_frac[currentCohort.pft] + leaf_c) *
        g_per_kg * days_per_sec * ha_per_m2

    currentSite.fmort_rate_cambial[currentCohort.size_class, currentCohort.pft] +=
        nc.n * currentCohort.cambial_mort / hlm_freq_day[]
    currentSite.fmort_rate_crown[currentCohort.size_class, currentCohort.pft] +=
        nc.n * currentCohort.crownfire_mort / hlm_freq_day[]

    # Fire kills some of the individuals in the new patch.
    nc.n = nc.n * (1.0 - currentCohort.fire_mort)
    _carry_mort_diags!(nc, currentCohort)

    # Burn off some leaf mass from living plants -> atmosphere.
    if prt_params.woody[currentCohort.pft] == itrue
        leaf_burn_frac = currentCohort.fraction_crown_burned
    else
        leaf_burn_frac = currentPatch.fuel.frac_burnt[live_grass(fuel_classes)]
    end

    if leaf_burn_frac < 0.0 || leaf_burn_frac > 1.0 ||
       currentCohort.fire_mort < 0.0 || currentCohort.fire_mort > 1.0
        fates_endrun("unexpected fire fractions")
    end

    for el in 1:num_elements[]
        if prt_params.woody[currentCohort.pft] == itrue
            leaf_m = GetState(nc.prt, leaf_organ, element_list[el])
        else
            leaf_m = GetState(nc.prt, leaf_organ, element_list[el]) +
                     GetState(nc.prt, sapw_organ, element_list[el]) +
                     GetState(nc.prt, struct_organ, element_list[el])
        end
        currentSite.mass_balance[el].burn_flux_to_atm += leaf_burn_frac * leaf_m * nc.n
        currentSite.flux_diags.elem[el].burned_liveveg += leaf_burn_frac * leaf_m * nc.n * area_inv
    end

    # Remove the burned mass from the plant.
    if prt_params.woody[currentCohort.pft] == itrue
        PRTBurnLosses!(nc.prt, leaf_organ, leaf_burn_frac)
    else
        PRTBurnLosses!(nc.prt, leaf_organ, leaf_burn_frac)
        PRTBurnLosses!(nc.prt, sapw_organ, leaf_burn_frac)
        PRTBurnLosses!(nc.prt, struct_organ, leaf_burn_frac)
    end

    currentCohort.fraction_crown_burned = 0.0
    nc.fraction_crown_burned = 0.0
    return nothing
end

# Logging survivorship (dtype_ilog).
function _spawn_surv_logging!(currentSite, currentPatch, nc, currentCohort,
                              patch_site_areadis, total_c, sapw_c, struct_c, store_c,
                              leaf_c, logging_coll_under_frac)
    if currentCohort.canopy_layer == 1
        # Upper-canopy: survivorship of degraded (non-harvested) trees.
        nc.n = currentCohort.n * currentCohort.l_degrad
        # Reduce the donor patch per the logging rate.
        currentCohort.n = currentCohort.n * (1.0 - min(1.0,
            currentCohort.lmort_direct + currentCohort.lmort_collateral +
            currentCohort.lmort_infra + currentCohort.l_degrad))
        _carry_mort_diags_nolog!(nc, currentCohort)
        # These weren't logged -> zero the logging mortality rates.
        nc.lmort_direct = 0.0; nc.lmort_collateral = 0.0; nc.lmort_infra = 0.0
    elseif prt_params.woody[currentCohort.pft] == itrue
        # Understory woody plants: step 1 reduce for area change.
        nc.n = currentCohort.n * patch_site_areadis / currentPatch.area

        currentSite.imort_rate[currentCohort.size_class, currentCohort.pft] +=
            nc.n * currentPatch.fract_ldist_not_harvested * logging_coll_under_frac / hlm_freq_day[]
        currentSite.imort_carbonflux[currentCohort.pft] +=
            (nc.n * currentPatch.fract_ldist_not_harvested * logging_coll_under_frac / hlm_freq_day[]) *
            total_c * g_per_kg * days_per_sec * years_per_day * ha_per_m2
        currentSite.imort_abg_flux[currentCohort.size_class, currentCohort.pft] +=
            (nc.n * currentPatch.fract_ldist_not_harvested * logging_coll_under_frac / hlm_freq_day[]) *
            ((sapw_c + struct_c + store_c) * prt_params.allom_agb_frac[currentCohort.pft] + leaf_c) *
            days_per_sec * years_per_day * ha_per_m2

        # Step 2 apply understory survivorship.
        nc.n = nc.n * (1.0 -
            (1.0 - currentPatch.fract_ldist_not_harvested) * logging_coll_under_frac)
        # Step 3 reduce donor count for area change.
        currentCohort.n = currentCohort.n * (1.0 - patch_site_areadis / currentPatch.area)
        _carry_mort_diags!(nc, currentCohort)
    else
        # Grass is not killed; just split it.
        nc.n = currentCohort.n * patch_site_areadis / currentPatch.area
        currentCohort.n = currentCohort.n * (1.0 - patch_site_areadis / currentPatch.area)
        _carry_mort_diags!(nc, currentCohort)
    end
    return nothing
end

# Carry over the donor cohort's mortality + logging diagnostics to the new cohort.
function _carry_mort_diags!(nc::fates_cohort_type, currentCohort::fates_cohort_type)
    nc.cmort  = currentCohort.cmort;  nc.hmort  = currentCohort.hmort
    nc.bmort  = currentCohort.bmort;  nc.frmort = currentCohort.frmort
    nc.smort  = currentCohort.smort;  nc.asmort = currentCohort.asmort
    nc.dgmort = currentCohort.dgmort; nc.dmort  = currentCohort.dmort
    nc.lmort_direct     = currentCohort.lmort_direct
    nc.lmort_collateral = currentCohort.lmort_collateral
    nc.lmort_infra      = currentCohort.lmort_infra
    return nothing
end

# Carry over mortality diagnostics WITHOUT the logging rates (logging upper canopy).
function _carry_mort_diags_nolog!(nc::fates_cohort_type, currentCohort::fates_cohort_type)
    nc.cmort  = currentCohort.cmort;  nc.hmort  = currentCohort.hmort
    nc.bmort  = currentCohort.bmort;  nc.frmort = currentCohort.frmort
    nc.smort  = currentCohort.smort;  nc.asmort = currentCohort.asmort
    nc.dgmort = currentCohort.dgmort; nc.dmort  = currentCohort.dmort
    return nothing
end

# nocomp + LUH patch-relabelling block at the tail of spawn_patches. Inert unless
# both flags are on (defaults off). Faithful transcription of the Fortran
# nocomp_and_luh_if block (buffer-patch carving + split_patch + fuse_2_patches).
function _spawn_remap_nocomp_luh!(currentSite::ed_site_type)
    for i_land_use_label in n_landuse_cats:-1:1
        nocomp_pft_area_vector        = zeros(numpft[])
        nocomp_pft_area_vector_filled = zeros(numpft[])

        copyPatch = nothing
        currentPatch = currentSite.oldest_patch
        while currentPatch !== nothing
            if currentPatch.changed_landuse_this_ts && currentPatch.land_use_label == i_land_use_label
                nocomp_pft_area_vector[currentPatch.nocomp_pft_label] += currentPatch.area
                copyPatch = currentPatch
            end
            currentPatch = currentPatch.younger
        end

        # How many PFTs on each land use type.
        n_pfts_by_landuse = 0
        which_pft_allowed = fates_unset_int
        for i_pft in 1:numpft[]
            if currentSite.area_PFT[i_pft, i_land_use_label] > nearzero
                n_pfts_by_landuse += 1
                which_pft_allowed = i_pft
            end
        end
        if n_pfts_by_landuse != 1
            which_pft_allowed = fates_unset_int
        end

        if sum(nocomp_pft_area_vector) > nearzero
            if n_pfts_by_landuse != 1
                _spawn_remap_buffer!(currentSite, i_land_use_label, copyPatch,
                    nocomp_pft_area_vector, nocomp_pft_area_vector_filled)
            else
                # Only one PFT allowed: relabel the changed patches.
                currentPatch = currentSite.oldest_patch
                while currentPatch !== nothing
                    if currentPatch.changed_landuse_this_ts && currentPatch.land_use_label == i_land_use_label
                        currentPatch.nocomp_pft_label = which_pft_allowed
                        currentPatch.changed_landuse_this_ts = false
                    end
                    currentPatch = currentPatch.younger
                end
            end
        end
        check_patch_area(currentSite)
    end
    return nothing
end

# The buffer-patch carving sub-block (n_pfts_by_landuse != 1). Kept as a faithful
# transcription; exercised only under nocomp+LUH (off by default).
function _spawn_remap_buffer!(currentSite, i_land_use_label, copyPatch,
                              nocomp_pft_area_vector, nocomp_pft_area_vector_filled)
    buffer_patch = fates_patch_type()
    Create!(buffer_patch, 0.0, 0.0, i_land_use_label, 0, num_swb, numpft[],
            currentSite.nlevsoil, hlm_current_tod[], ed_params().regeneration_model)
    for el in 1:num_elements[]
        InitConditions!(buffer_patch.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
    buffer_patch.tallest  = nothing
    buffer_patch.shortest = nothing
    copyPatch !== nothing && CopyPatchMeansTimers(copyPatch, buffer_patch)

    buffer_patch_in_linked_list = false
    buffer_patch_used = false

    currentPatch = currentSite.oldest_patch
    while currentPatch !== nothing
        if currentPatch.changed_landuse_this_ts && currentPatch.land_use_label == i_land_use_label
            area_to_keep = currentSite.area_PFT[currentPatch.nocomp_pft_label, i_land_use_label] *
                sum(nocomp_pft_area_vector) -
                nocomp_pft_area_vector_filled[currentPatch.nocomp_pft_label]
            newp_area = currentPatch.area - area_to_keep
            fraction_to_keep = area_to_keep / currentPatch.area

            if fraction_to_keep <= nearzero || area_to_keep < rsnbl_math_prec
                currentPatch.nocomp_pft_label = 0
                previousPatch = currentPatch.older !== nothing ? currentPatch.older : currentPatch
                fuse_2_patches(currentSite, currentPatch, buffer_patch)
                currentPatch = previousPatch
                buffer_patch_used = true
            elseif area_to_keep >= rsnbl_math_prec && newp_area >= rsnbl_math_prec
                temp_patch = fates_patch_type()
                split_patch(currentSite, currentPatch, temp_patch, fraction_to_keep, newp_area)
                temp_patch.nocomp_pft_label = 0
                fuse_2_patches(currentSite, temp_patch, buffer_patch)
                nocomp_pft_area_vector_filled[currentPatch.nocomp_pft_label] += currentPatch.area
                currentPatch.changed_landuse_this_ts = false
                buffer_patch_used = true
            else
                nocomp_pft_area_vector_filled[currentPatch.nocomp_pft_label] += currentPatch.area
                currentPatch.changed_landuse_this_ts = false
            end
        end
        currentPatch = currentPatch.younger
    end

    if buffer_patch_used
        newp_area_vector = (currentSite.area_PFT[:, i_land_use_label] .* sum(nocomp_pft_area_vector)) .-
            nocomp_pft_area_vector_filled
        newp_area_buffer_frac = newp_area_vector ./ buffer_patch.area
        max_val = maximum(newp_area_buffer_frac)

        if abs(sum(newp_area_buffer_frac) - max_val) <= nearzero
            i_pft = 1
            while !buffer_patch_in_linked_list
                if abs(newp_area_buffer_frac[i_pft] - max_val) <= nearzero
                    buffer_patch.nocomp_pft_label = i_pft
                    nocomp_pft_area_vector_filled[i_pft] += buffer_patch.area
                    InsertPatch(currentSite, buffer_patch)
                    buffer_patch_in_linked_list = true
                end
                i_pft += 1
            end
        end

        for i_pft in 1:numpft[]
            if currentSite.area_PFT[i_pft, i_land_use_label] > nearzero && !buffer_patch_in_linked_list
                nocomp_pft_area_vector_alt = copy(nocomp_pft_area_vector)
                nocomp_pft_area_vector_alt[i_pft] = 0.0
                newp_area = (currentSite.area_PFT[i_pft, i_land_use_label] * nocomp_pft_area_vector[i_pft]) -
                    nocomp_pft_area_vector_filled[i_pft]
                newp_area += sum(currentSite.area_PFT[i_pft, i_land_use_label] .* nocomp_pft_area_vector_alt)

                area_to_keep = buffer_patch.area - newp_area
                fraction_to_keep = area_to_keep / buffer_patch.area

                if newp_area > rsnbl_math_prec * 0.01
                    if area_to_keep > rsnbl_math_prec
                        temp_patch = fates_patch_type()
                        split_patch(currentSite, buffer_patch, temp_patch, fraction_to_keep, newp_area)
                        temp_patch.nocomp_pft_label = i_pft
                        nocomp_pft_area_vector_filled[i_pft] += temp_patch.area
                        InsertPatch(currentSite, temp_patch)
                    else
                        buffer_patch.nocomp_pft_label = i_pft
                        nocomp_pft_area_vector_filled[i_pft] += buffer_patch.area
                        InsertPatch(currentSite, buffer_patch)
                        buffer_patch_in_linked_list = true
                    end
                end
            end
        end

        if !buffer_patch_in_linked_list && buffer_patch.area >= rsnbl_math_prec
            fates_endrun("Buffer patch still has area and it wasn't put into the linked list")
        end
    end
    # (buffer_patch with no area is simply dropped; Julia is GC'd.)
    return nothing
end

# ===========================================================================
# fuse_2_patches
# ===========================================================================

"""
    fuse_2_patches(csite, dp, rp)

Fuse the donor patch `dp` into the recipient patch `rp` (area-weighted averages
of all the patch-level means + running means + litter + fuel, then splice all of
`dp`'s cohorts into `rp`'s height-ordered list), then unlink `dp` from the
site's age-ordered list (relinking its older/younger neighbors). Mirrors the
Fortran `fuse_2_patches`.

Upstream quirk preserved: `rp.area += dp.area` MUST come last (after the
area-weighted averages, whose weights use the un-summed areas).
"""
function fuse_2_patches(csite::ed_site_type, dp::fates_patch_type, rp::fates_patch_type)
    inv_sum_area = 1.0 / (dp.area + rp.area)

    rp.age = (dp.age * dp.area + rp.age * rp.area) * inv_sum_area
    rp.age_since_anthro_disturbance =
        (dp.age_since_anthro_disturbance * dp.area +
         rp.age_since_anthro_disturbance * rp.area) * inv_sum_area
    rp.age_class = get_age_class_index(rp.age)

    for el in 1:num_elements[]
        FuseLitter!(rp.litter[el], rp.area, dp.area, dp.litter[el])
    end
    fuse!(rp.fuel, rp.area, dp.area, dp.fuel)

    if rp.land_use_label != dp.land_use_label
        fates_endrun("trying to fuse patches with different land_use_label values")
    end
    if hlm_use_nocomp[] == itrue && rp.nocomp_pft_label != dp.nocomp_pft_label
        fates_endrun("trying to fuse patches with different nocomp_pft_label values")
    end

    # Weighted mean of the running means.
    FuseRMean!(rp.tveg24, dp.tveg24, rp.area * inv_sum_area)
    FuseRMean!(rp.tveg_lpa, dp.tveg_lpa, rp.area * inv_sum_area)
    if ed_params().regeneration_model == TRS_regeneration
        FuseRMean!(rp.seedling_layer_par24, dp.seedling_layer_par24, rp.area * inv_sum_area)
        FuseRMean!(rp.sdlng_mort_par, dp.sdlng_mort_par, rp.area * inv_sum_area)
        FuseRMean!(rp.sdlng2sap_par, dp.sdlng2sap_par, rp.area * inv_sum_area)
        for pft in 1:numpft[]
            FuseRMean!(rp.sdlng_emerg_smp[pft].p, dp.sdlng_emerg_smp[pft].p, rp.area * inv_sum_area)
            FuseRMean!(rp.sdlng_mdd[pft].p, dp.sdlng_mdd[pft].p, rp.area * inv_sum_area)
        end
    end
    FuseRMean!(rp.tveg_longterm, dp.tveg_longterm, rp.area * inv_sum_area)

    rp.livegrass  = (dp.livegrass * dp.area + rp.livegrass * rp.area) * inv_sum_area
    rp.ros_front  = (dp.ros_front * dp.area + rp.ros_front * rp.area) * inv_sum_area
    rp.tau_l      = (dp.tau_l * dp.area + rp.tau_l * rp.area) * inv_sum_area
    rp.tfc_ros    = (dp.tfc_ros * dp.area + rp.tfc_ros * rp.area) * inv_sum_area
    rp.fi         = (dp.fi * dp.area + rp.fi * rp.area) * inv_sum_area
    rp.fd         = (dp.fd * dp.area + rp.fd * rp.area) * inv_sum_area
    rp.ros_back   = (dp.ros_back * dp.area + rp.ros_back * rp.area) * inv_sum_area
    rp.scorch_ht  .= (dp.scorch_ht .* dp.area .+ rp.scorch_ht .* rp.area) .* inv_sum_area
    rp.frac_burnt = (dp.frac_burnt * dp.area + rp.frac_burnt * rp.area) * inv_sum_area
    rp.btran_ft   .= (dp.btran_ft .* dp.area .+ rp.btran_ft .* rp.area) .* inv_sum_area
    rp.zstar      = (dp.zstar * dp.area + rp.zstar * rp.area) * inv_sum_area
    rp.c_stomata  = (dp.c_stomata * dp.area + rp.c_stomata * rp.area) * inv_sum_area
    rp.c_lblayer  = (dp.c_lblayer * dp.area + rp.c_lblayer * rp.area) * inv_sum_area
    rp.rad_error[1] = (dp.rad_error[1] * dp.area + rp.rad_error[1] * rp.area) * inv_sum_area
    rp.rad_error[2] = (dp.rad_error[2] * dp.area + rp.rad_error[2] * rp.area) * inv_sum_area

    rp.area = rp.area + dp.area  # THIS MUST COME AT THE END!

    # Insert donor cohorts into the recipient patch.
    if dp.shortest !== nothing
        currentCohort = dp.shortest
        nextc = currentCohort.taller

        while dp.shortest !== nothing
            _insert_new_cohort!(rp, currentCohort)

            currentCohort = nextc
            dp.shortest = currentCohort
            if currentCohort !== nothing
                nextc = currentCohort.taller
            end
        end
    end

    patch_pft_size_profile(rp)

    # Relink dp's neighbors before dropping it.
    olderp   = dp.older
    youngerp = dp.younger

    FreeMemory!(dp, ed_params().regeneration_model, numpft[])

    if youngerp !== nothing || olderp !== nothing
        if youngerp !== nothing
            youngerp.older = olderp
        else
            csite.youngest_patch = olderp
            olderp.younger = nothing
        end
        if olderp !== nothing
            olderp.younger = youngerp
        else
            csite.oldest_patch = youngerp
            youngerp.older = nothing
        end
    end
    return nothing
end

# ===========================================================================
# fuse_patches
# ===========================================================================

"""
    fuse_patches(csite, bc_in)

Fuse patches whose (pft x dbh) aboveground-biomass profiles are sufficiently
similar (within the patch-fusion tolerance), per land-use category and (under
nocomp) per nocomp PFT label, relaxing the tolerance until the patch count is at
or below `maxpatches_by_landuse`. Two old patches (both above
`max_age_of_second_oldest_patch`), or two with tiny biomass
(`force_patchfuse_min_biomass`), are force-fused. Mirrors the Fortran
`fuse_patches`.
"""
function fuse_patches(csite::ed_site_type, bc_in::bc_in_type)
    currentSite = csite
    profiletol = ed_params().ED_val_patch_fusion_tol

    primary_land_fraction_beforefusion = 0.0
    primary_land_fraction_afterfusion  = 0.0

    nopatches = zeros(Int, n_landuse_cats)
    num_bareground_patches = 0

    currentPatch = currentSite.youngest_patch
    while currentPatch !== nothing
        if currentPatch.land_use_label > nocomp_bareground_land
            nopatches[currentPatch.land_use_label] += 1
            if currentPatch.land_use_label == primaryland
                primary_land_fraction_beforefusion += currentPatch.area * area_inv
            end
        else
            num_bareground_patches += 1
        end
        currentPatch = currentPatch.older
    end

    if num_bareground_patches > 1
        fates_endrun("somehow there is more than one bare ground patch")
    end

    pftlabelmin = 0
    pftlabelmax = hlm_use_nocomp[] == itrue ? numpft[] : 0

    for i_lulabel in 1:n_landuse_cats
        iterate = 1
        while iterate == 1
            for i_pftlabel in pftlabelmin:pftlabelmax
                # Recompute the biomass profile of each patch.
                currentPatch = currentSite.youngest_patch
                while currentPatch !== nothing
                    patch_pft_size_profile(currentPatch)
                    currentPatch = currentPatch.older
                end

                currentPatch = currentSite.youngest_patch
                while currentPatch !== nothing
                    tpp = currentSite.youngest_patch
                    while tpp !== nothing
                        if tpp.land_use_label == i_lulabel && currentPatch.land_use_label == i_lulabel &&
                           (hlm_use_nocomp[] == ifalse ||
                            (tpp.nocomp_pft_label == i_pftlabel && currentPatch.nocomp_pft_label == i_pftlabel))

                            fuse_flag = 1
                            if currentPatch.patchno != tpp.patchno  # not the same patch
                                if tpp.age <= max_age_of_second_oldest_patch ||
                                   currentPatch.age <= max_age_of_second_oldest_patch

                                    if sum(currentPatch.pft_agb_profile) > force_patchfuse_min_biomass ||
                                       sum(tpp.pft_agb_profile) > force_patchfuse_min_biomass

                                        for ft in 1:numpft[]
                                            for z in 1:N_DBH_BINS
                                                if currentPatch.pft_agb_profile[ft, z] > 0.0 ||
                                                   tpp.pft_agb_profile[ft, z] > 0.0
                                                    norm = abs(currentPatch.pft_agb_profile[ft, z] -
                                                        tpp.pft_agb_profile[ft, z]) /
                                                        (0.5 * (currentPatch.pft_agb_profile[ft, z] +
                                                                tpp.pft_agb_profile[ft, z]))
                                                    if norm > profiletol
                                                        fuse_flag = 0  # keep apart
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end

                                if fuse_flag == 1
                                    tmpptr = currentPatch.older
                                    fuse_2_patches(csite, currentPatch, tpp)
                                    fuse_cohorts(csite, tpp, bc_in)
                                    sort_cohorts(tpp)
                                    currentPatch = tmpptr
                                    # Reset tolerance after any fusion in this loop.
                                    profiletol = ed_params().ED_val_patch_fusion_tol
                                end
                            end
                        end

                        tpp = tpp.older
                    end

                    currentPatch = currentPatch === nothing ? nothing : currentPatch.older
                end
            end

            # Count patches in this land-use category.
            nopatches[i_lulabel] = 0
            currentPatch = currentSite.youngest_patch
            while currentPatch !== nothing
                if currentPatch.land_use_label == i_lulabel
                    nopatches[i_lulabel] += 1
                end
                currentPatch = currentPatch.older
            end

            if nopatches[i_lulabel] > maxpatches_by_landuse()[i_lulabel]
                iterate = 1
                profiletol *= patch_fusion_tolerance_relaxation_increment
                if profiletol > 100.0
                    fates_endrun("profile tolerance is too big, this shouldn't happen")
                end
            else
                iterate = 0
            end
        end
    end

    currentPatch = currentSite.youngest_patch
    while currentPatch !== nothing
        if currentPatch.land_use_label == primaryland
            primary_land_fraction_afterfusion += currentPatch.area * area_inv
        end
        currentPatch = currentPatch.older
    end

    currentSite.primary_land_patchfusion_error =
        primary_land_fraction_afterfusion - primary_land_fraction_beforefusion
    return nothing
end

# `maxpatches_by_landuse` is an EDParams field (Vector{Int}). Wrap the read so
# the call-sites read like the Fortran module parameter.
maxpatches_by_landuse() = ed_params().maxpatches_by_landuse

# ===========================================================================
# terminate_patches
# ===========================================================================

"""
    terminate_patches(currentSite, bc_in)

Fuse away patches whose area is below [`min_patch_area`](@ref) into a neighbor
(preferring the older patch of the same land-use label; the youngest patch is
spared unless it is below `min_patch_area_forced`). Mirrors the Fortran
`terminate_patches`.

NOTE: the Fortran's max-cycles fallback (a stuck tiny secondary/primary-only
patch that cannot fuse, triggering a land-use-type removal sweep) is reached only
in degenerate LUH configurations; we keep the cycle-count machinery but, like the
Fortran's "you can test your luck by disabling the endrun", raise on the genuine
no-progress case rather than running the elaborate LUH removal sweep (which needs
GetLUHStatedata to be meaningfully populated). Inert when LUH is off.
"""
function terminate_patches(currentSite::ed_site_type, bc_in::bc_in_type)
    max_cycles = 10
    count_cycles = 0

    currentPatch = currentSite.youngest_patch
    while currentPatch !== nothing
        if currentPatch.area <= min_patch_area
            if hlm_use_nocomp[] == itrue
                gotfused = false
                patchpointer = currentSite.youngest_patch
                while patchpointer !== nothing
                    if patchpointer !== currentPatch &&
                       patchpointer.nocomp_pft_label == currentPatch.nocomp_pft_label &&
                       patchpointer.land_use_label == currentPatch.land_use_label &&
                       !gotfused
                        fuse_2_patches(currentSite, patchpointer, currentPatch)
                        gotfused = true
                    else
                        patchpointer = patchpointer.older
                    end
                end
            else
                # Spare the youngest patch unless excessively small.
                if currentPatch !== currentSite.youngest_patch ||
                   currentPatch.area <= min_patch_area_forced

                    gotfused = false
                    if currentPatch.older !== nothing
                        olderPatch = currentPatch.older
                        if currentPatch.land_use_label == olderPatch.land_use_label
                            fuse_2_patches(currentSite, olderPatch, currentPatch)
                            gotfused = true
                        elseif count_cycles > 0
                            currentPatch.land_use_label = olderPatch.land_use_label
                            currentPatch.age_since_anthro_disturbance =
                                olderPatch.age_since_anthro_disturbance
                            fuse_2_patches(currentSite, olderPatch, currentPatch)
                            gotfused = true
                        end
                    end

                    if !gotfused && currentPatch.younger !== nothing
                        youngerPatch = currentPatch.younger
                        if currentPatch.land_use_label == youngerPatch.land_use_label
                            fuse_2_patches(currentSite, youngerPatch, currentPatch)
                            gotfused = true
                        elseif count_cycles > 0
                            currentPatch.land_use_label = youngerPatch.land_use_label
                            currentPatch.age_since_anthro_disturbance =
                                youngerPatch.age_since_anthro_disturbance
                            fuse_2_patches(currentSite, youngerPatch, currentPatch)
                            gotfused = true
                        end
                    end
                end
            end
        end

        # Don't move forward until enough area has merged into this patch.
        if currentPatch.area > min_patch_area_forced
            currentPatch = currentPatch.older
            count_cycles = 0
        else
            count_cycles += 1
        end

        if count_cycles > max_cycles
            # The Fortran here runs an elaborate LUH land-use-type-removal sweep
            # (only meaningful with populated LUH driver data) or else endruns.
            # We endrun on the genuine no-progress case.
            fates_endrun("FATES is having difficulties fusing very small patches.")
        end
    end

    check_patch_area(currentSite)
    return nothing
end
