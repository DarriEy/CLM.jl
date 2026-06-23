# EDLoggingMortalityMod.jl
# Julia port of FATES src/fates/biogeochem/EDLoggingMortalityMod.F90
#
# Logging / harvest disturbance mortality. Purpose:
#   1. Create logging mortalities (cohort level):
#      (a) direct logging mortality, (b) collateral mortality,
#      (c) infrastructure (mechanical) mortality.
#   2. Move the logged trunk fluxes from live into the (exported) product pool.
#   3. Move logging-associated mortality fluxes from live to CWD / litter.
#   4. Keep the carbon balance (resources-management + mass-balance accounting).
#
# Original authors: Yi Xu & M. Huang (09/2017).
#
# Deps (already ported):
#   FatesConstantsMod  : fates_r8->Float64, rsnbl_math_prec, primaryland,
#                        secondaryland, secondary_age_threshold, n_landuse_cats,
#                        fates_tiny, itrue, ifalse, months_per_year, days_per_sec,
#                        years_per_day, g_per_kg, fates_check_param_set
#   FatesCohortMod     : fates_cohort_type
#   FatesPatchMod      : fates_patch_type
#   EDTypesMod         : ed_site_type, site_massbal_type, site_fluxdiags_type,
#                        elem_diag_type, ed_resources_management_type, area, area_inv,
#                        get_secondary_young_fraction
#   FatesLitterMod     : ncwd, ndcmpy, litter_type, adjust_SF_CWD_frac
#   EDPftvarcon        : edpftvarcon_inst (harvest_pprod10, landusechange_pprod10),
#                        GetDecompyFrac
#   EDParamsMod        : ed_params() (logging_* parameter fields)
#   PRTParametersMod   : prt_params (woody, allom_agb_frac)
#   PRTGenericMod      : num_elements, element_list, element_pos, carbon12_element,
#                        sapw_organ/struct_organ/leaf_organ/fnrt_organ/store_organ/
#                        repro_organ, GetState
#   SFParamsMod        : sf_params() (SF_val_CWD_frac)
#   FatesAllometryMod  : set_root_fraction, carea_allom
#   FatesInterfaceTypesMod : bc_in_type, bc_out_type, hlm_* flags + dates, numpft
#   FatesLandUseChangeMod : GetInitLanduseHarvestRate, GetLUHStatedata
#
# NOT-yet-ported dependency (stubbed below, only reached when plant hydraulics
# is enabled — off in the default/test path):
#   FatesPlantHydraulicsMod : AccumulateMortalityWaterStorage
#
# Standalone — NOT added to CLMInstances or any dual-copied struct.

# ---------------------------------------------------------------------------
# Module-level state / parameters (Fortran module-save variables)
# ---------------------------------------------------------------------------

# If true, logging should be performed during the current time-step. Mirrors the
# Fortran `logging_time` module-save flag (a Ref so the set-by-IsItLoggingTime,
# read-elsewhere pattern is preserved). Other modules read it as logging_time[].
const logging_time = Ref{Bool}(false)

# Harvest litter localization: how much of the litter from a falling tree lands
# within the newly generated patch vs the original patch. 0 => no preference
# (mass distributed equally by area); 1 => completely local to the new patch.
const harvest_litter_localization = 0.0

# Constants used in hlm_harvest_units comparisons. The Fortran imports these from
# FatesConstantsMod; they are interface unit-type identifiers.
const hlm_harvest_area_fraction = 1  # area-fraction-based harvest [m2/m2]
const hlm_harvest_carbon        = 2  # carbon/biomass-based harvest

# TODO Batch NN: AccumulateMortalityWaterStorage lives in FatesPlantHydraulicsMod,
# which is partially ported but does not yet expose this routine. It is only
# called when hlm_use_planthydro == itrue (off in the default + test path). A
# minimal no-op stub keeps logging_litter_fluxes self-consistent until the real
# routine lands.
function AccumulateMortalityWaterStorage(currentSite, currentCohort, delta_n::Real)
    # TODO Batch NN: real implementation in FatesPlantHydraulicsMod.
    return nothing
end

# ===========================================================================
# IsItLoggingTime
# ===========================================================================
"""
    IsItLoggingTime(is_master, currentSite)

Determine whether the current dynamics step should enact the logging module, by
comparing the current model time to the logging event code (`ed_params()
.logging_event_code`). Sets the module flag [`logging_time`](@ref) and resets the
site's per-event resources-management diagnostics. Mirrors the Fortran
`IsItLoggingTime`.

Event-code semantics (icode = int(logging_event_code)):
  * 1                : logging off
  * 2                : logging on the first model day
  * 3                : logging every day
  * 4                : logging on the first day of each month
  * -1..-365         : logging once a year on |icode| day-of-year
  * > 10000          : a specific event YYYYMMDD
"""
function IsItLoggingTime(is_master::Integer, currentSite::ed_site_type)
    logging_time[] = false
    icode = Int(trunc(ed_params().logging_event_code))

    # This is true for either hlm harvest or fates logging.
    if hlm_use_logging[] == ifalse
        return nothing
    end

    if icode == 1
        # Logging is turned off
        logging_time[] = false
    elseif icode == 2
        # Logging event on the first step
        if hlm_model_day[] == 1
            logging_time[] = true
        end
    elseif icode == 3
        # Logging event every day
        logging_time[] = true
    elseif icode == 4
        # Logging event once a month
        if hlm_current_day[] == 1
            logging_time[] = true
        end
    elseif icode < 0 && icode > -366
        # Logging event every year on a specific day of year
        if hlm_day_of_year[] == abs(icode)
            logging_time[] = true
        end
    elseif icode > 10000
        # Specific Event: YYYYMMDD
        log_date  = icode - Int(100 * floor(icode / 100))
        log_year  = Int(floor(icode / 10000))
        log_month = Int(floor(icode / 100)) - log_year * 100

        if hlm_current_day[] == log_date &&
           hlm_current_month[] == log_month &&
           hlm_current_year[] == log_year
            logging_time[] = true
        end
    else
        # Bad logging event flag
        fates_endrun("An invalid logging code was specified in fates_params. " *
                     "Check EDLoggingMortalityMod.jl:IsItLoggingTime for a " *
                     "breakdown of the valid codes and change " *
                     "fates_logging_event_code in the file accordingly.")
    end

    # Initialize some site-level diagnostics that are calculated for each event.
    currentSite.resources_management.delta_litter_stock  = 0.0
    currentSite.resources_management.delta_biomass_stock = 0.0
    currentSite.resources_management.delta_individual    = 0.0

    return nothing
end

# ===========================================================================
# LoggingMortality_frac
# ===========================================================================
"""
    LoggingMortality_frac(currentSite, bc_in, pft_i, dbh, canopy_layer,
                          hlm_harvest_rates, hlm_harvest_catnames, hlm_harvest_units,
                          patch_land_use_label, secondary_age, frac_site_primary,
                          frac_site_secondary, harvestable_forest_c)
        -> (lmort_direct, lmort_collateral, lmort_infra, l_degrad, harvest_tag)

Compute the per-cohort logging mortality fractions for the given PFT / dbh /
canopy layer under the active harvest scenario. Returns:
  * `lmort_direct`     : direct (harvestable) mortality fraction
  * `lmort_collateral` : collateral-damage mortality fraction
  * `lmort_infra`      : infrastructure (mechanical) mortality fraction
  * `l_degrad`         : fraction of trees not killed but moved to degraded
                         (newly-anthro-disturbed secondary) forest
  * `harvest_tag`      : per-harvest-category status tag (0 success, 1 not enough
                         carbon, 2 not applicable), length `hlm_num_lu_harvest_cats`

Out-arguments of the Fortran subroutine are returned as a tuple (the cohort-level
`lmort_*`/`l_degrad` fields and `harvest_tag` are written by the caller). Mirrors
the Fortran `LoggingMortality_frac`, preserving the unintuitive dbhmax logic.
"""
function LoggingMortality_frac(currentSite::ed_site_type, bc_in::bc_in_type,
                               pft_i::Integer, dbh::Real, canopy_layer::Integer,
                               hlm_harvest_rates::AbstractVector{<:Real},
                               hlm_harvest_catnames::AbstractVector{<:AbstractString},
                               hlm_harvest_units::Integer,
                               patch_land_use_label::Integer, secondary_age::Real,
                               frac_site_primary::Real, frac_site_secondary::Real,
                               harvestable_forest_c::AbstractVector{<:Real})
    edp = ed_params()

    # Outputs (defaults match Fortran intent(out) flow)
    lmort_direct     = 0.0
    lmort_collateral = 0.0
    lmort_infra      = 0.0
    l_degrad         = 0.0
    harvest_tag      = fill(2, hlm_num_lu_harvest_cats[])

    cur_harvest_tag = 2
    harvest_rate    = 0.0

    # The transition_landuse_from_off_to_on special case handles the first
    # timestep after leaving potential-vegetation mode: all prior historical
    # land-use (including harvest) is applied on that first day.
    if !currentSite.transition_landuse_from_off_to_on

        # Check if the secondaryland exceeds the minimum if in landuse mode.
        site_secondaryland_first_exceeding_min = false
        if hlm_use_luh[] == itrue
            state_vector = GetLUHStatedata(bc_in)
            site_secondaryland_first_exceeding_min =
                (state_vector[secondaryland] > currentSite.min_allowed_landuse_fraction) &&
                (!currentSite.landuse_vector_gt_min[secondaryland])
        end

        if site_secondaryland_first_exceeding_min
            # Special logic when the intended secondary area first exceeds the
            # too-small-patch minimum.
            state_vector = GetLUHStatedata(bc_in)
            if patch_land_use_label == primaryland
                harvest_rate = state_vector[secondaryland] / state_vector[primaryland]
            else
                harvest_rate = 0.0
            end
            # For area-based harvest, harvest_tag shall always be 2 (not applicable).
            fill!(harvest_tag, 2)
            cur_harvest_tag = 2

        elseif logging_time[]

            # Pass logging rates to cohort level
            if hlm_use_lu_harvest[] == ifalse
                # 0 = use fates logging parameters directly: harvest the whole
                # cohort area.
                harvest_rate = 1.0
                # UPSTREAM-FORTRAN QUIRK: on this path the Fortran never assigns
                # cur_harvest_tag (it is an uninitialized local), yet the direct-
                # logging block below gates on `cur_harvest_tag == 0`. The fates-
                # parameter path is intended to directly log, so we make the
                # "successful harvest" tag (0) explicit here.
                cur_harvest_tag = 0

            elseif hlm_use_lu_harvest[] == itrue &&
                   hlm_harvest_units == hlm_harvest_area_fraction
                # 1 = use area fraction from hlm.
                secondary_young_fraction = get_secondary_young_fraction(currentSite)
                harvest_rate = get_harvest_rate_area(patch_land_use_label,
                    hlm_harvest_catnames, hlm_harvest_rates, frac_site_primary,
                    frac_site_secondary, secondary_young_fraction, secondary_age)
                # For area-based harvest, harvest_tag shall always be 2.
                fill!(harvest_tag, 2)
                cur_harvest_tag = 2

            elseif hlm_use_lu_harvest[] == itrue &&
                   hlm_harvest_units == hlm_harvest_carbon
                # 2 = use carbon from hlm.
                harvest_rate, cur_harvest_tag = get_harvest_rate_carbon(
                    patch_land_use_label, hlm_harvest_catnames, hlm_harvest_rates,
                    secondary_age, harvestable_forest_c, harvest_tag)
            end

        else
            harvest_rate = 0.0
            fill!(harvest_tag, 2)
            cur_harvest_tag = 2
        end

        # transfer of area to secondary land is based on overall area affected,
        # not just logged crown area; l_degrad accounts for the affected area
        # between logged crowns.
        if prt_params.woody[pft_i] == itrue  # only set logging rates for trees
            if cur_harvest_tag == 0
                # direct logging rates, based on dbh min and max criteria.
                # NOTE: the dbhmax comparison logic is deliberately unintuitive in
                # the Fortran (the .and. .not. after the first conditional means
                # the dbh:dbhmax comparison is inverted) — preserved verbatim so
                # the dbhmax check can be turned off entirely via fates_check_param_set.
                if dbh >= edp.logging_dbhmin && !(
                       (edp.logging_dbhmax < fates_check_param_set) &&
                       (dbh >= edp.logging_dbhmax))
                    lmort_direct = harvest_rate * edp.logging_direct_frac
                else
                    lmort_direct = 0.0
                end
            else
                lmort_direct = 0.0
            end

            # infrastructure (roads, skid trails, etc.) mortality rates
            if dbh >= edp.logging_dbhmax_infra
                lmort_infra = 0.0
            else
                lmort_infra = harvest_rate * edp.logging_mechanical_frac
            end

            # Collateral damage to smaller plants below the direct-logging size
            # threshold is applied via "understory_death" in the disturbance
            # algorithm; here only the canopy layer takes collateral damage.
            if canopy_layer == 1
                lmort_collateral = harvest_rate * edp.logging_collateral_frac
            else
                lmort_collateral = 0.0
            end

        else  # non-woody plants still killed by infrastructure
            lmort_direct     = 0.0
            lmort_collateral = 0.0
            lmort_infra      = harvest_rate * edp.logging_mechanical_frac
        end

        # the area occupied by canopy plants that aren't killed is still disturbed
        # at the harvest rate.
        if canopy_layer == 1
            l_degrad = harvest_rate - (lmort_direct + lmort_infra + lmort_collateral)
        else
            l_degrad = 0.0
        end

    else
        # First timestep after leaving potential-vegetation mode.
        harvest_rate = GetInitLanduseHarvestRate(bc_in,
            currentSite.min_allowed_landuse_fraction, currentSite.landuse_vector_gt_min)
        lmort_direct     = 0.0
        lmort_collateral = 0.0
        lmort_infra      = 0.0
        l_degrad         = 0.0
        if prt_params.woody[pft_i] == itrue
            lmort_direct = harvest_rate
        elseif canopy_layer == 1
            l_degrad = harvest_rate
        end
    end

    return lmort_direct, lmort_collateral, lmort_infra, l_degrad, harvest_tag
end

# ===========================================================================
# get_harvest_rate_area
# ===========================================================================
"""
    get_harvest_rate_area(patch_land_use_label, hlm_harvest_catnames,
                          hlm_harvest_rates, frac_site_primary, frac_site_secondary,
                          secondary_young_fraction, secondary_age) -> harvest_rate

Get the area-based harvest rate from the host boundary-condition info, given the
patch land-use history. Assumes `logging_time == true`. Forest categories only
(non-forest harvest has a geographical mismatch with LUH2). Normalizes by the
site primary / secondary fraction (and the young/old split for secondary), caps
at 1, and applies the per-day/-month divisor for periodic event codes. Mirrors
the Fortran `get_harvest_rate_area`.
"""
function get_harvest_rate_area(patch_land_use_label::Integer,
                               hlm_harvest_catnames::AbstractVector{<:AbstractString},
                               hlm_harvest_rates::AbstractVector{<:Real},
                               frac_site_primary::Real, frac_site_secondary::Real,
                               secondary_young_fraction::Real, secondary_age::Real)
    # Loop over harvest categories to determine the annual hlm harvest rate.
    harvest_rate = 0.0
    for h_index in 1:hlm_num_lu_harvest_cats[]
        if patch_land_use_label == primaryland
            if hlm_harvest_catnames[h_index] == "HARVEST_VH1" ||
               hlm_harvest_catnames[h_index] == "HARVEST_VH2"
                harvest_rate += hlm_harvest_rates[h_index]
            end
        elseif patch_land_use_label == secondaryland &&
               secondary_age >= secondary_age_threshold
            if hlm_harvest_catnames[h_index] == "HARVEST_SH1"
                harvest_rate += hlm_harvest_rates[h_index]
            end
        elseif patch_land_use_label == secondaryland &&
               secondary_age < secondary_age_threshold
            if hlm_harvest_catnames[h_index] == "HARVEST_SH2" ||
               hlm_harvest_catnames[h_index] == "HARVEST_SH3"
                harvest_rate += hlm_harvest_rates[h_index]
            end
        end
    end

    # Normalize by site-level primary or secondary forest fraction (harvest_rate
    # is specified as a fraction of the gridcell), capping so we do not harvest
    # more primary/secondary area than there is. For secondary, also normalize by
    # the young/old fraction.
    if patch_land_use_label == primaryland
        if frac_site_primary > fates_tiny
            harvest_rate = min(harvest_rate / frac_site_primary, 1.0)
        else
            harvest_rate = 0.0
        end
    elseif patch_land_use_label == secondaryland
        # frac_site_secondary returns -1 if no secondary area, hence the .gt. -0.5.
        if frac_site_secondary > fates_tiny && frac_site_secondary > -0.5
            if secondary_age < secondary_age_threshold
                harvest_rate = min(harvest_rate /
                    (frac_site_secondary * secondary_young_fraction), 1.0)
            else
                harvest_rate = min(harvest_rate /
                    (frac_site_secondary * (1.0 - secondary_young_fraction)), 1.0)
            end
        else
            harvest_rate = 0.0
        end
    else
        harvest_rate = 0.0
    end

    # Apply today's harvest rate. Whether to harvest today has already been
    # determined by IsItLoggingTime. For icode == 2, < 0, > 10000 apply the
    # annual rate one time (no calc). Bad event flag is caught in IsItLoggingTime.
    icode = Int(trunc(ed_params().logging_event_code))
    if icode == 1
        harvest_rate = 0.0
    elseif icode == 3
        # Logging event every day
        harvest_rate = harvest_rate / hlm_days_per_year[]
    elseif icode == 4
        # logging event once a month
        if hlm_current_day[] == 1
            harvest_rate = harvest_rate / months_per_year
        end
    end

    return harvest_rate
end

# ===========================================================================
# get_harvestable_carbon
# ===========================================================================
"""
    get_harvestable_carbon(csite, site_area, hlm_harvest_catnames)
        -> harvestable_forest_c

Total carbon available for harvest for the three FATES harvest categories
(primary, secondary mature, secondary young), aggregating all woody cohorts that
meet the dbhmin/dbhmax criteria. Called outside the patch loop; the output drives
the carbon-based harvest rate. Returns a vector of length
`hlm_num_lu_harvest_cats` [kgC/site]. Mirrors the Fortran `get_harvestable_carbon`.
"""
function get_harvestable_carbon(csite::ed_site_type, site_area::Real,
                                hlm_harvest_catnames::AbstractVector{<:AbstractString})
    edp = ed_params()
    SF_val_CWD_frac = sf_params().SF_val_CWD_frac

    harvestable_forest_c = zeros(Float64, hlm_num_lu_harvest_cats[])

    currentPatch = csite.oldest_patch
    while currentPatch !== nothing
        harvestable_patch_c = 0.0
        currentCohort = currentPatch.tallest

        while currentCohort !== nothing
            pft = currentCohort.pft

            if Int(prt_params.woody[pft]) == 1  # only set logging rates for trees
                sapw_m   = GetState(currentCohort.prt, sapw_organ, carbon12_element)
                struct_m = GetState(currentCohort.prt, struct_organ, carbon12_element)
                # unit: [kgC] = [kgC/plant] * [plant/ha] * [ha/10k m2] * [m2 area]
                harvestable_cohort_c = edp.logging_direct_frac * (sapw_m + struct_m) *
                    prt_params.allom_agb_frac[currentCohort.pft] *
                    SF_val_CWD_frac[ncwd] * edp.logging_export_frac *
                    currentCohort.n * area_inv * site_area

                # No harvest for trees without canopy.
                if currentCohort.canopy_layer >= 1
                    if currentCohort.dbh >= edp.logging_dbhmin && !(
                           (edp.logging_dbhmax < fates_check_param_set) &&
                           (currentCohort.dbh >= edp.logging_dbhmax))
                        harvestable_patch_c += harvestable_cohort_c
                    end
                end
            end
            currentCohort = currentCohort.shorter
        end

        # Judge which category the current patch belongs to. Since forest/non-forest
        # is not separated, all carbon belongs to the forest categories.
        for h_index in 1:hlm_num_lu_harvest_cats[]
            if currentPatch.land_use_label == primaryland
                if hlm_harvest_catnames[h_index] == "HARVEST_VH1"
                    harvestable_forest_c[h_index] += harvestable_patch_c
                end
            elseif currentPatch.land_use_label == secondaryland &&
                   currentPatch.age_since_anthro_disturbance >= secondary_age_threshold
                if hlm_harvest_catnames[h_index] == "HARVEST_SH1"
                    harvestable_forest_c[h_index] += harvestable_patch_c
                end
            elseif currentPatch.land_use_label == secondaryland &&
                   currentPatch.age_since_anthro_disturbance < secondary_age_threshold
                if hlm_harvest_catnames[h_index] == "HARVEST_SH2"
                    harvestable_forest_c[h_index] += harvestable_patch_c
                end
            end
        end
        currentPatch = currentPatch.younger
    end

    return harvestable_forest_c
end

# ===========================================================================
# get_harvest_rate_carbon
# ===========================================================================
"""
    get_harvest_rate_carbon(patch_land_use_label, hlm_harvest_catnames,
                            hlm_harvest_rates, secondary_age, harvestable_forest_c,
                            harvest_tag) -> (harvest_rate, cur_harvest_tag)

Get the carbon-based harvest rate, converting a carbon demand into an area
fraction using the available harvestable forest carbon. Assumes
`logging_time == true`. Updates `harvest_tag` in place (0 success / 1 not enough
carbon / 2 not applicable) and returns the area-fraction `harvest_rate` plus the
minimum tag over categories. Mirrors the Fortran `get_harvest_rate_carbon`.
"""
function get_harvest_rate_carbon(patch_land_use_label::Integer,
                                 hlm_harvest_catnames::AbstractVector{<:AbstractString},
                                 hlm_harvest_rates::AbstractVector{<:Real},
                                 secondary_age::Real,
                                 harvestable_forest_c::AbstractVector{<:Real},
                                 harvest_tag::AbstractVector{<:Integer})
    harvest_rate        = 0.0
    harvest_rate_c      = 0.0
    harvest_rate_supply = 0.0
    fill!(harvest_tag, 2)

    # Get the carbon-demand harvest rate from HLM (three FATES logging types:
    # primary, secondary mature, secondary young).
    for h_index in 1:hlm_num_lu_harvest_cats[]
        if patch_land_use_label == primaryland
            if hlm_harvest_catnames[h_index] == "HARVEST_VH1" ||
               hlm_harvest_catnames[h_index] == "HARVEST_VH2"
                harvest_rate_c += hlm_harvest_rates[h_index]
            end
        elseif patch_land_use_label == secondaryland &&
               secondary_age >= secondary_age_threshold
            if hlm_harvest_catnames[h_index] == "HARVEST_SH1"
                harvest_rate_c += hlm_harvest_rates[h_index]
            end
        elseif patch_land_use_label == secondaryland &&
               secondary_age < secondary_age_threshold
            if hlm_harvest_catnames[h_index] == "HARVEST_SH2" ||
               hlm_harvest_catnames[h_index] == "HARVEST_SH3"
                harvest_rate_c += hlm_harvest_rates[h_index]
            end
        end
    end

    # Determine harvest status (successful or not). Only three categories used.
    for h_index in 1:hlm_num_lu_harvest_cats[]
        if patch_land_use_label == primaryland
            if hlm_harvest_catnames[h_index] == "HARVEST_VH1"
                if harvestable_forest_c[h_index] >= harvest_rate_c
                    harvest_rate_supply += harvestable_forest_c[h_index]
                    harvest_tag[h_index] = 0
                else
                    harvest_tag[h_index] = 1
                end
            end
        elseif patch_land_use_label == secondaryland &&
               secondary_age >= secondary_age_threshold
            if hlm_harvest_catnames[h_index] == "HARVEST_SH1"
                if harvestable_forest_c[h_index] >= harvest_rate_c
                    harvest_rate_supply += harvestable_forest_c[h_index]
                    harvest_tag[h_index] = 0
                else
                    harvest_tag[h_index] = 1
                end
            end
        elseif patch_land_use_label == secondaryland &&
               secondary_age < secondary_age_threshold
            if hlm_harvest_catnames[h_index] == "HARVEST_SH2"
                if harvestable_forest_c[h_index] >= harvest_rate_c
                    harvest_rate_supply += harvestable_forest_c[h_index]
                    harvest_tag[h_index] = 0
                else
                    harvest_tag[h_index] = 1
                end
            end
        end
    end

    # If any harvest category available, the cur_harvest_tag triggers the event.
    cur_harvest_tag = minimum(harvest_tag)

    # Transfer carbon-based harvest rate to area-based harvest rate.
    if harvest_rate_supply > rsnbl_math_prec && harvest_rate_supply > harvest_rate_c
        harvest_rate = harvest_rate_c / harvest_rate_supply
    else
        # If we forced harvest_rate to 1 when we don't have enough C, we'd produce
        # a primary patch with no area (cannot be terminated under nocomp). Keep 0.
        harvest_rate = 0.0
    end

    # Prevent the generation of tiny secondary patches.
    if harvest_rate < 1e-8
        harvest_rate = 0.0
    end

    # Apply today's harvest rate (no site-fraction normalization for carbon mode).
    icode = Int(trunc(ed_params().logging_event_code))
    if icode == 1
        harvest_rate = 0.0
    elseif icode == 3
        harvest_rate = harvest_rate / hlm_days_per_year[]
    elseif icode == 4
        if hlm_current_day[] == 1
            harvest_rate = harvest_rate / months_per_year
        end
    end

    return harvest_rate, cur_harvest_tag
end

# ===========================================================================
# logging_litter_fluxes
# ===========================================================================
"""
    logging_litter_fluxes!(currentSite, currentPatch, newPatch, patch_site_areadis,
                           bc_in)

Route the biomass of logging-killed trees into CWD / fine-litter pools (split
between the donor `currentPatch` and the newly-disturbed `newPatch` by area), and
export the directly-logged above-ground trunk (times `logging_export_frac`) to
the wood-product pool. Generates the fluxes used in carbon-balance checking and
the resources-management diagnostics. Only called when logging disturbance is the
dominant disturbance. Mirrors the Fortran `logging_litter_fluxes`.

The litter losses are almost exactly the natural tree-fall case; the differences
are the mortality rates governing the fluxes and the exported trunk product.
"""
function logging_litter_fluxes!(currentSite::ed_site_type,
                                currentPatch::fates_patch_type,
                                newPatch::fates_patch_type,
                                patch_site_areadis::Real, bc_in::bc_in_type)
    edp = ed_params()
    SF_val_CWD_frac = sf_params().SF_val_CWD_frac
    SF_val_CWD_frac_adj = zeros(Float64, ncwd)

    nlevsoil = currentSite.nlevsoil

    # When sending litter fluxes to the old patch we divide the mass sent there by
    # the area it will have remaining after it donates area.
    remainder_area = currentPatch.area - patch_site_areadis

    # Fraction of litter retained in the donor patch vs donated to the new patch.
    retain_frac = (1.0 - harvest_litter_localization) *
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
        cur_litt     = currentPatch.litter[el]  # litter pool of the "current" patch
        new_litt     = newPatch.litter[el]      # litter pool of the "new" patch

        # Zero some site-level accumulator diagnostics.
        trunk_product_site  = 0.0
        delta_litter_stock  = 0.0
        delta_biomass_stock = 0.0
        delta_individual    = 0.0

        # -----------------------------------------------------------------------
        # Part 1: Send parts of dying plants to the litter pool.
        # -----------------------------------------------------------------------
        currentCohort = currentPatch.shortest
        while currentCohort !== nothing
            pft = currentCohort.pft

            sapw_m   = GetState(currentCohort.prt, sapw_organ, element_id)
            struct_m = GetState(currentCohort.prt, struct_organ, element_id)
            leaf_m   = GetState(currentCohort.prt, leaf_organ, element_id)
            fnrt_m   = GetState(currentCohort.prt, fnrt_organ, element_id)
            store_m  = GetState(currentCohort.prt, store_organ, element_id)
            repro_m  = GetState(currentCohort.prt, repro_organ, element_id)

            if currentCohort.canopy_layer == 1
                direct_dead   = currentCohort.n * currentCohort.lmort_direct
                indirect_dead = currentCohort.n *
                    (currentCohort.lmort_collateral + currentCohort.lmort_infra)
            else
                # This routine is only called during disturbance. No direct dead
                # can occur in the understory here; indirect are impacts.
                if prt_params.woody[pft] == itrue
                    direct_dead   = 0.0
                    indirect_dead = edp.logging_coll_under_frac *
                        (1.0 - currentPatch.fract_ldist_not_harvested) * currentCohort.n *
                        (patch_site_areadis / currentPatch.area)  # kgC/site/day
                else
                    # Grass experiences no logging-disturbance mortality.
                    direct_dead   = 0.0
                    indirect_dead = 0.0
                end
            end

            if element_id == carbon12_element && hlm_use_planthydro[] == itrue
                AccumulateMortalityWaterStorage(currentSite, currentCohort,
                                                direct_dead + indirect_dead)
            end

            # Distribute woody litter (non-bole components) between current and new
            # patches. Get the per-layer root fraction for the below-ground split.
            set_root_fraction(currentSite.rootfrac_scr, pft, currentSite.zi_soil;
                              max_nlevroot = bc_in.max_rooting_depth_index_col)

            ag_wood = (direct_dead + indirect_dead) * (struct_m + sapw_m) *
                prt_params.allom_agb_frac[currentCohort.pft]
            bg_wood = (direct_dead + indirect_dead) * (struct_m + sapw_m) *
                (1.0 - prt_params.allom_agb_frac[currentCohort.pft])

            # Adjust how wood is partitioned between CWD classes based on dbh.
            adjust_SF_CWD_frac(currentCohort.dbh, ncwd, SF_val_CWD_frac, SF_val_CWD_frac_adj)

            for c in 1:(ncwd - 1)
                new_litt.ag_cwd[c] += ag_wood * SF_val_CWD_frac_adj[c] * donate_m2
                cur_litt.ag_cwd[c] += ag_wood * SF_val_CWD_frac_adj[c] * retain_m2

                for ilyr in 1:nlevsoil
                    new_litt.bg_cwd[c, ilyr] += bg_wood *
                        currentSite.rootfrac_scr[ilyr] * SF_val_CWD_frac_adj[c] * donate_m2
                    cur_litt.bg_cwd[c, ilyr] += bg_wood *
                        currentSite.rootfrac_scr[ilyr] * SF_val_CWD_frac_adj[c] * retain_m2
                end

                # Diagnostics on fluxes into the AG and BG CWD pools.
                elflux_diags.cwd_ag_input[c] += SF_val_CWD_frac_adj[c] * ag_wood
                elflux_diags.cwd_bg_input[c] += SF_val_CWD_frac_adj[c] * bg_wood

                if element_id == carbon12_element
                    delta_litter_stock += (ag_wood + bg_wood) * SF_val_CWD_frac_adj[c]
                end
            end

            # Trunk wood of infrastructure and collateral-damage mortality.
            ag_wood = indirect_dead * (struct_m + sapw_m) *
                prt_params.allom_agb_frac[currentCohort.pft]
            bg_wood = indirect_dead * (struct_m + sapw_m) *
                (1.0 - prt_params.allom_agb_frac[currentCohort.pft])

            new_litt.ag_cwd[ncwd] += ag_wood * SF_val_CWD_frac_adj[ncwd] * donate_m2
            cur_litt.ag_cwd[ncwd] += ag_wood * SF_val_CWD_frac_adj[ncwd] * retain_m2

            for ilyr in 1:nlevsoil
                new_litt.bg_cwd[ncwd, ilyr] += bg_wood *
                    currentSite.rootfrac_scr[ilyr] * SF_val_CWD_frac_adj[ncwd] * donate_m2
                cur_litt.bg_cwd[ncwd, ilyr] += bg_wood *
                    currentSite.rootfrac_scr[ilyr] * SF_val_CWD_frac_adj[ncwd] * retain_m2
            end

            elflux_diags.cwd_ag_input[ncwd] += SF_val_CWD_frac_adj[ncwd] * ag_wood
            elflux_diags.cwd_bg_input[ncwd] += SF_val_CWD_frac_adj[ncwd] * bg_wood

            if element_id == carbon12_element
                delta_litter_stock += (ag_wood + bg_wood) * SF_val_CWD_frac_adj[ncwd]
            end

            # Below-ground trunk flux for directly-logged trees (c = ncwd).
            bg_wood = direct_dead * (struct_m + sapw_m) * SF_val_CWD_frac_adj[ncwd] *
                (1.0 - prt_params.allom_agb_frac[currentCohort.pft])

            for ilyr in 1:nlevsoil
                new_litt.bg_cwd[ncwd, ilyr] += bg_wood *
                    currentSite.rootfrac_scr[ilyr] * donate_m2
                cur_litt.bg_cwd[ncwd, ilyr] += bg_wood *
                    currentSite.rootfrac_scr[ilyr] * retain_m2
            end

            elflux_diags.cwd_bg_input[ncwd] += bg_wood

            # Harvest (export) flux for the above-ground boles. A fraction
            # (export_frac) of the boles from direct logging are exported off-site;
            # the remainder (1-export_frac) is added to the litter pools.
            ag_wood = direct_dead * (struct_m + sapw_m) *
                prt_params.allom_agb_frac[currentCohort.pft] * SF_val_CWD_frac_adj[ncwd]

            trunk_product_site += ag_wood * edp.logging_export_frac

            # For the total mass balance [kg/site/day].
            site_mass.wood_product_harvest[pft] += ag_wood * edp.logging_export_frac

            new_litt.ag_cwd[ncwd] += ag_wood * (1.0 - edp.logging_export_frac) * donate_m2
            cur_litt.ag_cwd[ncwd] += ag_wood * (1.0 - edp.logging_export_frac) * retain_m2

            # Fluxes of leaf, root, and storage carbon into litter pools (none
            # exported).
            leaf_litter = (direct_dead + indirect_dead) * (leaf_m + repro_m)
            root_litter = (direct_dead + indirect_dead) * (fnrt_m + store_m)

            for dcmpy in 1:ndcmpy
                dcmpy_frac = GetDecompyFrac(pft, leaf_organ, dcmpy)
                new_litt.leaf_fines[dcmpy] += leaf_litter * donate_m2 * dcmpy_frac
                cur_litt.leaf_fines[dcmpy] += leaf_litter * retain_m2 * dcmpy_frac

                dcmpy_frac = GetDecompyFrac(pft, fnrt_organ, dcmpy)
                for ilyr in 1:nlevsoil
                    new_litt.root_fines[dcmpy, ilyr] += root_litter *
                        currentSite.rootfrac_scr[ilyr] * dcmpy_frac * donate_m2
                    cur_litt.root_fines[dcmpy, ilyr] += root_litter *
                        currentSite.rootfrac_scr[ilyr] * dcmpy_frac * retain_m2
                end
            end

            # Track as diagnostic fluxes.
            elflux_diags.surf_fine_litter_input[pft] += leaf_litter
            elflux_diags.root_litter_input[pft]      += root_litter

            # Logging-specific diagnostics (litter stock also has CWD terms above).
            if element_id == carbon12_element
                delta_litter_stock  += leaf_litter + root_litter
                delta_biomass_stock += leaf_litter + root_litter +
                    (direct_dead + indirect_dead) * (struct_m + sapw_m)
                delta_individual    += direct_dead + indirect_dead
            end

            currentCohort = currentCohort.taller
        end

        # Amount of trunk mass exported off site [kg/m2].
        elflux_diags.exported_harvest += trunk_product_site * area_inv

        # Update the carbon exported from the site through logging operations.
        if element_id == carbon12_element
            currentSite.resources_management.trunk_product_site  += trunk_product_site
            currentSite.resources_management.delta_litter_stock  += delta_litter_stock
            currentSite.resources_management.delta_biomass_stock += delta_biomass_stock
            currentSite.resources_management.delta_individual    += delta_individual
        end
    end

    # Recompute crown area of the new patch's cohorts (rgk 06-2019: not sure why
    # this is called here, but it can't hurt).
    currentCohort = newPatch.shortest
    while currentCohort !== nothing
        _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
            currentSite.spread, currentCohort.pft, currentCohort.crowndamage)
        currentCohort = currentCohort.taller
    end

    return nothing
end

# ===========================================================================
# UpdateHarvestC
# ===========================================================================
"""
    UpdateHarvestC!(currentSite, bc_out)

Update the harvested-C flux passed back to the HLM (`hrv_deadstemc_to_prod10c` /
`hrv_deadstemc_to_prod100c`), splitting the per-PFT harvest and land-use-change
wood-product pools into 10-yr and 100-yr product pools by the PFT
`harvest_pprod10` / `landusechange_pprod10` fractions. Converts kgC m-2 day-1 to
gC m-2 s-1. Mirrors the Fortran `UpdateHarvestC`. Added by Shijie Shu.
"""
function UpdateHarvestC!(currentSite::ed_site_type, bc_out::bc_out_type)
    edpf = edpftvarcon_inst()

    # Flush the older value before update.
    bc_out.hrv_deadstemc_to_prod10c  = 0.0
    bc_out.hrv_deadstemc_to_prod100c = 0.0

    # Unit transfer factor (from kgC m-2 day-1 to gC m-2 s-1).
    unit_trans_factor = g_per_kg * days_per_sec

    el = element_pos[carbon12_element]
    site_mass = currentSite.mass_balance[el]

    # harvest-associated wood product pools
    for i_pft in 1:numpft[]
        bc_out.hrv_deadstemc_to_prod10c += site_mass.wood_product_harvest[i_pft] *
            area_inv * edpf.harvest_pprod10[i_pft] * unit_trans_factor
        bc_out.hrv_deadstemc_to_prod100c += site_mass.wood_product_harvest[i_pft] *
            area_inv * (1.0 - edpf.harvest_pprod10[i_pft]) * unit_trans_factor
    end

    # land-use-change-associated wood product pools
    for i_pft in 1:numpft[]
        bc_out.hrv_deadstemc_to_prod10c += site_mass.wood_product_landusechange[i_pft] *
            area_inv * edpf.landusechange_pprod10[i_pft] * unit_trans_factor
        bc_out.hrv_deadstemc_to_prod100c += site_mass.wood_product_landusechange[i_pft] *
            area_inv * (1.0 - edpf.landusechange_pprod10[i_pft]) * unit_trans_factor
    end

    return nothing
end

# ===========================================================================
# get_harvest_debt
# ===========================================================================
"""
    get_harvest_debt!(site_in, bc_in, harvest_tag)

Accumulate the site-level harvest debt (kgC/site not successfully harvested
because the carbon available was below the forcing-data harvest rate) for
primary, secondary-mature, and secondary-young land, using the per-category
`harvest_tag` (1 => not enough carbon). Only acts when `logging_time == true`.
Mirrors the Fortran `get_harvest_debt`.
"""
function get_harvest_debt!(site_in::ed_site_type, bc_in::bc_in_type,
                           harvest_tag::AbstractVector{<:Integer})
    if !logging_time[]
        return nothing
    end

    harvest_debt_pri        = 0.0
    harvest_debt_sec_mature = 0.0
    harvest_debt_sec_young  = 0.0

    # First get the harvest rate demand for all three categories.
    for h_index in 1:hlm_num_lu_harvest_cats[]
        if bc_in.hlm_harvest_catnames[h_index] == "HARVEST_VH1" ||
           bc_in.hlm_harvest_catnames[h_index] == "HARVEST_VH2"
            harvest_debt_pri += bc_in.hlm_harvest_rates[h_index]
        elseif bc_in.hlm_harvest_catnames[h_index] == "HARVEST_SH1"
            harvest_debt_sec_mature += bc_in.hlm_harvest_rates[h_index]
        elseif bc_in.hlm_harvest_catnames[h_index] == "HARVEST_SH2" ||
               bc_in.hlm_harvest_catnames[h_index] == "HARVEST_SH3"
            harvest_debt_sec_young += bc_in.hlm_harvest_rates[h_index]
        end
    end

    # Next accumulate the harvest debt through the harvest tag.
    for h_index in 1:hlm_num_lu_harvest_cats[]
        if harvest_tag[h_index] == 1
            if bc_in.hlm_harvest_catnames[h_index] == "HARVEST_VH1"
                site_in.resources_management.harvest_debt += harvest_debt_pri
            elseif bc_in.hlm_harvest_catnames[h_index] == "HARVEST_SH1"
                site_in.resources_management.harvest_debt += harvest_debt_sec_mature
                site_in.resources_management.harvest_debt_sec += harvest_debt_sec_mature
            elseif bc_in.hlm_harvest_catnames[h_index] == "HARVEST_SH2"
                site_in.resources_management.harvest_debt += harvest_debt_sec_young
                site_in.resources_management.harvest_debt_sec += harvest_debt_sec_young
            end
        end
    end

    return nothing
end
