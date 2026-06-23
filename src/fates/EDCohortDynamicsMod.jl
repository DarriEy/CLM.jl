# EDCohortDynamicsMod.jl
# Julia port of FATES src/fates/biogeochem/EDCohortDynamicsMod.F90
#
# The demographic ENGINE of FATES — manages the cohort lifecycle on the
# self-referential `taller`/`shorter` cohort linked list attached to each patch:
#
#   * create_cohort        — allocate, Create(), and splice a new cohort into the
#                            height-ordered list (recruits / cold-start / inventory
#                            / restart).
#   * InitPRTObject!       — allocate the PARTEH allocation object for the chosen
#                            hypothesis (carbon-only or flexible-CNP).
#   * terminate_cohorts    — sweep the list and remove cohorts that are too small,
#                            outside the canopy, or terminally carbon-depleted
#                            (level 1 = numerically dangerous; level 2 = the rest).
#   * terminate_cohort     — remove one cohort: tally site-level termination
#                            fluxes, send its mass to litter, relink the list, free.
#   * SendCohortToLitter    — transfer a select number of plants' organ masses into
#                            the patch litter / CWD pools (NOT turnover; does NOT
#                            change per-plant pools or number density).
#   * fuse_cohorts          — merge cohorts that are within size (and age, and
#                            damage, and canopy-layer, and pft, and is-new) fusion
#                            tolerances; conserves carbon (and crown area + number,
#                            re-solving dbh) while N-weighting every other scalar.
#                            Iterates with a growing tolerance until the patch is
#                            under `max_cohort_per_patch`.
#   * sort_cohorts          — rebuild the height-ordered list by re-inserting every
#                            cohort (DO NOT CHANGE — order matters for fusion).
#   * insert_cohort         — splice one cohort into the height-ordered list,
#                            updating the patch tallest/shortest head/tail.
#   * count_cohorts         — count the list both directions, set countcohorts.
#   * EvaluateAndCorrectDBH — bump dbh/height up if structural (or, grasses, leaf)
#                            carbon exceeds the allometric target (integration fix).
#   * DamageRecovery        — split a damaged cohort, creating a recovered (lower
#                            damage class) copy if it has the resources.
#
# Translation notes (per project conventions):
#   * `fates_r8` -> Float64; cohort/patch/site/litter structs are reused VERBATIM
#     from FatesCohortMod / FatesPatchMod / EDTypesMod / FatesLitterMod — this
#     module does NOT redefine any of them.
#   * Fortran `allocate(newCohort)` + pointer splice -> construct an empty
#     `fates_cohort_type()` and walk/relink the `Union{...,Nothing}` list pointers.
#   * `deallocate(currentCohort)` after `terminate_cohort` -> Julia is GC'd; the
#     cohort is already unlinked + FreeMemory()'d, so we just drop it.
#   * `insert_cohort` carries `storebigcohort`/`storesmallcohort` as `Ref`s so the
#     pass-by-pointer head/tail updates survive (Fortran intent(inout) pointer
#     args). The optional args are always supplied by our callers.
#   * The Fortran `select case(cohort_fusion_conservation_method)` is a module
#     parameter fixed to `conserve_crownarea_and_number_not_dbh`; preserved.
#   * Plant-hydraulics paths (hlm_use_planthydro==itrue) call the ported
#     FatesPlantHydraulicsMod helpers where they exist; the few cohort-fusion /
#     recruit-constraint / mortality-water helpers that are not yet ported are
#     guarded and flagged `# TODO` (they are inert when plant hydraulics is off,
#     which is the default and what the test exercises).
#   * Module flag globals (hlm_use_planthydro / hlm_use_cohort_age_tracking /
#     hlm_use_sp / hlm_parteh_mode) are `Ref{Int}` -> dereferenced with `[]`.
#
# Upstream-Fortran quirks preserved (see inline comments):
#   * `year_net_uptake == 999` sentinel handling during fusion (min, not weighted).
#   * The c13disc_acc GPP-weighted mean guards a 0/0 by zeroing.
#   * `daily_n_demand`/`daily_p_demand` carry the -9 sentinel through copy/fuse.
#   * fuse_cohorts compares EVERY pair (the inner `nextc` loop restarts from the
#     tallest each time) — O(n^2) by design; do not "optimize" the traversal.
#   * sort_cohorts header literally says "DO NOT CHANGE THIS IT WILL BREAK".
#
# Deps: FatesCohortMod (fates_cohort_type + Create/Copy/Init/ZeroValues/FreeMemory/
#   UpdateCohortBioPhysRates/InitPRTBoundaryConditions), FatesPatchMod
#   (fates_patch_type), EDTypesMod (ed_site_type, min_npm2/min_nppatch/
#   min_n_safemath, site_fluxdiags_type/elem_diag_type), FatesLitterMod
#   (litter_type, ncwd, ndcmpy, adjust_SF_CWD_frac), SFParamsMod (sf_params().SF_val_CWD_frac),
#   EDPftvarcon (edpftvarcon_inst, GetDecompyFrac), PRTParametersMod (prt_params),
#   EDParamsMod (nclmax, nlevleaf, EDParams[].max_cohort_per_patch,
#   ED_val_cohort_size/age_fusion_tol), PRTGenericMod (num_elements, element_list,
#   organ/element ids, GetState, WeightedFusePRTVartypes!, StorageNutrientTarget,
#   the carbon/CNP hypothesis ids + InitPRTVartype! + callom/cnp_allom_prt_vartypes),
#   FatesAllometryMod (carea_allom, h_allom, ForceDBH, bsap/bagw/bbgw/bdead/bleaf/
#   bfineroot/bstore_allom, set_root_fraction), FatesSizeAgeTypeIndicesMod
#   (sizetype_class_index, coagetype_class_index), DamageMainMod (undamaged_class),
#   FatesConstantsMod (itrue/ifalse, nearzero, calloc_abs_error, ican_upper,
#   i_term_mort_type_*, leaves_off/leaves_shedding, ihard/isemi_stress_decid),
#   FatesPlantHydraulicsMod (the hydro helpers, guarded), FatesGlobals (fates_endrun).
# Standalone — NOT added to CLMInstances or any dual-copied struct.

# ---------------------------------------------------------------------------
# Module-level fusion-conservation method selector (Fortran module parameters)
# ---------------------------------------------------------------------------
const conserve_crownarea_and_number_not_dbh = 1
const conserve_dbh_and_number_not_crownarea = 2

# Fixed choice, matching the Fortran module parameter.
const cohort_fusion_conservation_method = conserve_crownarea_and_number_not_dbh

# ===========================================================================
# create_cohort
# ===========================================================================

"""
    create_cohort(currentSite, patchptr, pft, nn, height, coage, dbh, prt,
                  elongf_leaf, elongf_fnrt, elongf_stem, status, recruitstatus,
                  ctrim, carea, clayer, crowndamage, spread, bc_in)

Create a new cohort and insert it at the correct height-ordered position in the
patch's cohort linked list. Called when initializing cohorts at cold start, when
recruiting during dynamics, during an inventory read, and during restart. The
PARTEH `prt` object must already be allocated and initialized. Mirrors the
Fortran `create_cohort`. `bc_in` is only used for the plant-hydraulics recruit
constraint.
"""
function create_cohort(currentSite::ed_site_type, patchptr::fates_patch_type,
                       pft::Integer, nn::Real, height::Real, coage::Real,
                       dbh::Real, prt::AbstractPRTVartypes, elongf_leaf::Real,
                       elongf_fnrt::Real, elongf_stem::Real, status::Integer,
                       recruitstatus::Integer, ctrim::Real, carea::Real,
                       clayer::Integer, crowndamage::Integer, spread::Real,
                       bc_in=nothing)

    # create new cohort
    newCohort = fates_cohort_type()
    Create(newCohort, prt, pft, nn, height, coage, dbh, status, ctrim, carea,
           clayer, crowndamage, spread, patchptr.canopy_layer_tlai,
           elongf_leaf, elongf_fnrt, elongf_stem)

    # Put cohort at the right place in the linked list.
    # storebigcohort / storesmallcohort carry the (possibly nil) head/tail.
    storebigcohort   = Ref{Union{fates_cohort_type,Nothing}}(patchptr.tallest)
    storesmallcohort = Ref{Union{fates_cohort_type,Nothing}}(patchptr.shortest)

    if patchptr.tallest !== nothing
        tnull = 0
    else
        tnull = 1
        patchptr.tallest = newCohort
    end

    if patchptr.shortest !== nothing
        snull = 0
    else
        snull = 1
        patchptr.shortest = newCohort
    end

    # Plant hydraulics initialization for the new cohort.
    if hlm_use_planthydro[] == itrue
        _create_cohort_inithydro!(currentSite, patchptr, newCohort, recruitstatus, bc_in)
    end

    insert_cohort(patchptr, newCohort, patchptr.tallest, patchptr.shortest,
                  tnull, snull, storebigcohort, storesmallcohort)

    patchptr.tallest  = storebigcohort[]
    patchptr.shortest = storesmallcohort[]

    return nothing
end

# Plant-hydraulics new-cohort initialization (mirrors the Fortran create_cohort
# hydro block). Only entered when hlm_use_planthydro==itrue.
function _create_cohort_inithydro!(currentSite::ed_site_type,
                                   patchptr::fates_patch_type,
                                   newCohort::fates_cohort_type,
                                   recruitstatus::Integer, bc_in)
    # Allocate the cohort hydro object + per-layer arrays.
    InitHydrCohort(currentSite, newCohort)
    newCohort.co_hydr.errh2o = 0.0

    # Node heights, then volumes/lengths/kmax (UpdateSizeDepPlantHydProps! bundles
    # SavePreviousCompartmentVolumes + UpdatePlantHydrNodes + UpdatePlantHydrLenVol
    # + UpdatePlantKmax). Re-snapshot the just-computed volumes as the "init".
    UpdatePlantHydrNodes!(newCohort, newCohort.pft, newCohort.height, currentSite.si_hydr)
    UpdateSizeDepPlantHydProps!(currentSite, newCohort, bc_in)
    SavePreviousCompartmentVolumes!(newCohort.co_hydr)

    # Starter suctions/water contents from the soil state.
    InitPlantHydStates!(currentSite, newCohort)

    if recruitstatus == 1
        newCohort.co_hydr.is_newly_recruited = true
        # Constrain the recruit number density to the water the rhizosphere can supply.
        rmean_temp = GetMean(patchptr.tveg24)
        ConstrainRecruitNumber(currentSite, newCohort, patchptr, bc_in, rmean_temp)
    end
    return nothing
end

# ===========================================================================
# InitPRTObject!
# ===========================================================================

"""
    InitPRTObject!() -> AbstractPRTVartypes

Allocate and initialize the PARTEH allocation object for the configured
hypothesis (`hlm_parteh_mode`): the carbon-only `callom_prt_vartypes` or the
flexible-CNP `cnp_allom_prt_vartypes`. The internal mappings are set up via
`InitPRTVartype!`; initial/boundary conditions are set by the caller. Mirrors the
Fortran `InitPRTObject` (which returns through a pointer argument).
"""
function InitPRTObject!()
    mode = hlm_parteh_mode[]
    if mode == prt_carbon_allom_hyp
        prt = callom_prt_vartypes()
    elseif mode == prt_cnp_flex_allom_hyp
        prt = cnp_allom_prt_vartypes()
    else
        fates_endrun("You specified an unknown PRT module. Aborting.")
    end

    InitPRTVartype!(prt)

    return prt
end

# ===========================================================================
# terminate_cohorts
# ===========================================================================

"""
    terminate_cohorts(currentSite, currentPatch, level, call_index, bc_in)

Sweep the patch cohort list (shortest -> taller) and remove cohorts that are too
small or terminally depleted. `level == 1` is called BEFORE fusion and removes
only the numerically dangerous (sub-`min_n_safemath`) cohorts; `level == 2`
removes the rest (low number density, outside the canopy, depleted live or total
biomass) but never recruits. Mirrors the Fortran `terminate_cohorts`.
"""
function terminate_cohorts(currentSite::ed_site_type, currentPatch::fates_patch_type,
                           level::Integer, call_index::Integer, bc_in=nothing)
    currentCohort = currentPatch.shortest
    while currentCohort !== nothing
        terminate        = ifalse
        termination_type = 0
        tallerCohort     = currentCohort.taller

        leaf_c   = GetState(currentCohort.prt, leaf_organ, carbon12_element)
        store_c  = GetState(currentCohort.prt, store_organ, carbon12_element)
        sapw_c   = GetState(currentCohort.prt, sapw_organ, carbon12_element)
        fnrt_c   = GetState(currentCohort.prt, fnrt_organ, carbon12_element)
        struct_c = GetState(currentCohort.prt, struct_organ, carbon12_element)

        # Check if number density is so low it breaks math (level 1).
        if currentCohort.n < min_n_safemath && level == 1
            terminate        = itrue
            termination_type = i_term_mort_type_numdens
        end

        # The rest are only allowed if we are not dealing with a recruit (level 2).
        if !currentCohort.isnew && level == 2

            # Not enough n or dbh.
            if (currentCohort.n / currentPatch.area <= min_npm2 ||
                currentCohort.n <= min_nppatch ||
                (currentCohort.dbh < 0.00001 && store_c < 0.0))
                terminate        = itrue
                termination_type = i_term_mort_type_numdens
            end

            # Outside the maximum canopy layer.
            if currentCohort.canopy_layer > nclmax
                terminate        = itrue
                termination_type = i_term_mort_type_canlev
            end

            # Live biomass pools are terminally depleted.
            if (sapw_c + leaf_c + fnrt_c) < 1e-10 || store_c < 1e-10
                terminate        = itrue
                termination_type = i_term_mort_type_cstarv
            end

            # Total cohort biomass is negative.
            if (struct_c + sapw_c + leaf_c + fnrt_c + store_c) < 0.0
                terminate        = itrue
                termination_type = i_term_mort_type_cstarv
            end
        end

        if terminate == itrue
            terminate_cohort(currentSite, currentPatch, currentCohort, bc_in,
                             termination_type)
            # FreeMemory + unlink already done in terminate_cohort; Julia GCs it.
        end

        currentCohort = tallerCohort
    end

    return nothing
end

# ===========================================================================
# terminate_cohort
# ===========================================================================

"""
    terminate_cohort(currentSite, currentPatch, currentCohort, bc_in, termination_type)

Terminate one cohort: tally the site-level termination number-density / carbon
fluxes for the appropriate canopy level, send the cohort's biomass to litter,
unlink it from the patch list (relinking neighbors and updating tallest/shortest),
and free its memory. Mirrors the Fortran `terminate_cohort`.
"""
function terminate_cohort(currentSite::ed_site_type, currentPatch::fates_patch_type,
                          currentCohort::fates_cohort_type, bc_in,
                          termination_type::Integer)
    # termination_type should never be 0.
    if termination_type == 0
        fates_endrun("termination_type=0")
    end

    leaf_c   = GetState(currentCohort.prt, leaf_organ, carbon12_element)
    store_c  = GetState(currentCohort.prt, store_organ, carbon12_element)
    sapw_c   = GetState(currentCohort.prt, sapw_organ, carbon12_element)
    fnrt_c   = GetState(currentCohort.prt, fnrt_organ, carbon12_element)
    struct_c = GetState(currentCohort.prt, struct_organ, carbon12_element)
    repro_c  = GetState(currentCohort.prt, repro_organ, carbon12_element)

    levcan = currentCohort.canopy_layer

    if hlm_use_planthydro[] == itrue
        AccumulateMortalityWaterStorage(currentSite, currentCohort, currentCohort.n)
    end

    sc   = currentCohort.size_class
    ipft = currentCohort.pft

    # Update the site-level termination tallies for the appropriate canopy layer.
    if levcan == ican_upper
        currentSite.term_nindivs_canopy[termination_type, sc, ipft] +=
            currentCohort.n
        currentSite.term_carbonflux_canopy[termination_type, ipft] +=
            currentCohort.n * (struct_c + sapw_c + leaf_c + fnrt_c + store_c + repro_c)
    else
        currentSite.term_nindivs_ustory[termination_type, sc, ipft] +=
            currentCohort.n
        currentSite.term_carbonflux_ustory[termination_type, ipft] +=
            currentCohort.n * (struct_c + sapw_c + leaf_c + fnrt_c + store_c + repro_c)
    end

    currentSite.term_abg_flux[sc, ipft] +=
        currentCohort.n * ((struct_c + sapw_c + store_c) *
        prt_params.allom_agb_frac[ipft] + leaf_c)

    # Put the litter from the terminated cohort straight into the fragmenting pools.
    if currentCohort.n > 0.0
        SendCohortToLitter(currentSite, currentPatch, currentCohort,
                           currentCohort.n, bc_in)
    end

    # Set pointers and remove the current cohort from the list.
    shorterCohort = currentCohort.shorter
    tallerCohort  = currentCohort.taller

    if tallerCohort === nothing
        currentPatch.tallest = shorterCohort
        shorterCohort !== nothing && (shorterCohort.taller = nothing)
    else
        tallerCohort.shorter = shorterCohort
    end

    if shorterCohort === nothing
        currentPatch.shortest = tallerCohort
        tallerCohort !== nothing && (tallerCohort.shorter = nothing)
    else
        shorterCohort.taller = tallerCohort
    end

    FreeMemory(currentCohort)

    return nothing
end

# ===========================================================================
# SendCohortToLitter
# ===========================================================================

"""
    SendCohortToLitter(csite, cpatch, ccohort, nplant, bc_in)

Transfer the existing organ masses (all pools, all elements) of `nplant` plants of
a cohort into the patch litter / CWD pools and accumulate the site flux diagnostics.

IMPORTANT (from the Fortran): this is NOT turnover and NOT a partial transfer; it
does NOT touch per-plant pools or the cohort number density; it is not used for
disturbance. Mirrors the Fortran `SendCohortToLitter`. `bc_in` supplies the
maximum rooting depth index for the root-fraction profile.
"""
function SendCohortToLitter(csite::ed_site_type, cpatch::fates_patch_type,
                            ccohort::fates_cohort_type, nplant::Real, bc_in=nothing)
    pft        = ccohort.pft
    plant_dens = nplant / cpatch.area

    # root-fraction profile into csite.rootfrac_scr. The Fortran passes
    # bc_in%max_rooting_depth_index_col positionally; the Julia port takes it as
    # the `max_nlevroot` keyword (defaulting to the full profile when absent).
    max_rooting_depth_index_col = bc_in === nothing ? csite.nlevsoil :
        bc_in.max_rooting_depth_index_col
    set_root_fraction(csite.rootfrac_scr, pft, csite.zi_soil;
                      max_nlevroot = max_rooting_depth_index_col)

    SF_val_CWD_frac_adj = zeros(Float64, ncwd)

    for el in 1:num_elements[]
        store_m = GetState(ccohort.prt, store_organ, element_list[el])
        fnrt_m  = GetState(ccohort.prt, fnrt_organ, element_list[el])
        repro_m = GetState(ccohort.prt, repro_organ, element_list[el])
        if prt_params.woody[ccohort.pft] == itrue
            leaf_m   = GetState(ccohort.prt, leaf_organ, element_list[el])
            sapw_m   = GetState(ccohort.prt, sapw_organ, element_list[el])
            struct_m = GetState(ccohort.prt, struct_organ, element_list[el])
        else
            # Non-woody: lump leaf+sapwood+structural all into "leaf" fines.
            leaf_m   = GetState(ccohort.prt, leaf_organ, element_list[el]) +
                       GetState(ccohort.prt, sapw_organ, element_list[el]) +
                       GetState(ccohort.prt, struct_organ, element_list[el])
            sapw_m   = 0.0
            struct_m = 0.0
        end

        litt         = cpatch.litter[el]
        elflux_diags = csite.flux_diags.elem[el]

        # Adjust how wood is partitioned between the cwd classes based on dbh.
        adjust_SF_CWD_frac(ccohort.dbh, ncwd, sf_params().SF_val_CWD_frac,
                           SF_val_CWD_frac_adj)

        for c in 1:ncwd
            # above ground CWD
            litt.ag_cwd[c] += plant_dens * (struct_m + sapw_m) *
                SF_val_CWD_frac_adj[c] * prt_params.allom_agb_frac[pft]

            # below ground CWD
            for sl in 1:csite.nlevsoil
                litt.bg_cwd[c, sl] += plant_dens * (struct_m + sapw_m) *
                    SF_val_CWD_frac_adj[c] *
                    (1.0 - prt_params.allom_agb_frac[pft]) * csite.rootfrac_scr[sl]
            end

            # above ground flux diagnostic
            elflux_diags.cwd_ag_input[c] += (struct_m + sapw_m) *
                SF_val_CWD_frac_adj[c] * prt_params.allom_agb_frac[pft] * nplant

            # below ground flux diagnostic
            elflux_diags.cwd_bg_input[c] += (struct_m + sapw_m) *
                SF_val_CWD_frac_adj[c] *
                (1.0 - prt_params.allom_agb_frac[pft]) * nplant
        end

        for dcmpy in 1:ndcmpy
            dcmpy_frac = GetDecompyFrac(pft, leaf_organ, dcmpy)
            litt.leaf_fines[dcmpy] += plant_dens * (leaf_m + repro_m) * dcmpy_frac

            dcmpy_frac = GetDecompyFrac(pft, fnrt_organ, dcmpy)
            for sl in 1:csite.nlevsoil
                litt.root_fines[dcmpy, sl] += plant_dens * (fnrt_m + store_m) *
                    csite.rootfrac_scr[sl] * dcmpy_frac
            end
        end

        elflux_diags.surf_fine_litter_input[pft] += (leaf_m + repro_m) * nplant
        elflux_diags.root_litter_input[pft]      += (fnrt_m + store_m) * nplant
    end

    return nothing
end

# ===========================================================================
# fuse_cohorts
# ===========================================================================

"""
    fuse_cohorts(currentSite, currentPatch, bc_in)

Join similar cohorts to reduce the total cohort count. Two cohorts fuse if they
are within the (dynamic) size and age fusion tolerances AND share pft, crown
damage class, canopy layer, and is-new status. On fusion all mass pools are
weighted-fused, crown area + number are conserved (and dbh re-solved from the
crown-area allometry, or weight-averaged + ForceDBH if above the cap), and every
remaining scalar is number-density-weighted. The fusion tolerance grows by 10%
each pass until the patch is under `max_cohort_per_patch`. Mirrors the Fortran
`fuse_cohorts`; sorts the list afterward if any fusion occurred.
"""
function fuse_cohorts(currentSite::ed_site_type, currentPatch::fates_patch_type,
                      bc_in=nothing)
    dynamic_size_fusion_tolerance = EDParams[].ED_val_cohort_size_fusion_tol
    dynamic_age_fusion_tolerance  = EDParams[].ED_val_cohort_age_fusion_tol

    iterate           = 1
    fusion_took_place = 0

    if currentPatch.shortest !== nothing
        while iterate == 1
            currentCohort = currentPatch.tallest

            # Loop until current cohort IS the shortest (the last has already been
            # compared by then).
            while currentCohort !== currentPatch.shortest && currentCohort !== nothing
                nextc = currentPatch.tallest

                while nextc !== nothing
                    nextnextc = nextc.shorter
                    diff = abs((currentCohort.dbh - nextc.dbh) /
                               (0.5 * (currentCohort.dbh + nextc.dbh)))

                    if diff < dynamic_size_fusion_tolerance
                        # Age tolerance (0 if equal, to avoid divide-by-zero).
                        if abs(currentCohort.coage - nextc.coage) < nearzero
                            coage_diff = 0.0
                        else
                            coage_diff = abs((currentCohort.coage - nextc.coage) /
                                (0.5 * (currentCohort.coage + nextc.coage)))
                        end

                        if coage_diff <= dynamic_age_fusion_tolerance &&
                           currentCohort !== nextc &&             # don't fuse with self
                           currentCohort.pft == nextc.pft &&
                           currentCohort.crowndamage == nextc.crowndamage &&
                           currentCohort.canopy_layer == nextc.canopy_layer &&
                           currentCohort.isnew == nextc.isnew

                            newn = currentCohort.n + nextc.n
                            fusion_took_place = 1

                            _fuse_pair!(currentSite, currentPatch, currentCohort,
                                        nextc, newn, bc_in)

                            currentCohort.n = newn

                            # Unlink nextc from the list.
                            shorterCohort = nextc.shorter
                            tallerCohort  = nextc.taller

                            if tallerCohort === nothing
                                currentPatch.tallest = shorterCohort
                                shorterCohort !== nothing && (shorterCohort.taller = nothing)
                            else
                                tallerCohort.shorter = shorterCohort
                            end

                            if shorterCohort === nothing
                                currentPatch.shortest = tallerCohort
                                tallerCohort !== nothing && (tallerCohort.shorter = nothing)
                            else
                                shorterCohort.taller = tallerCohort
                            end

                            if hlm_use_planthydro[] == itrue
                                # update hydraulics quantities that are functions
                                # of height + biomasses (helper ported).
                                UpdateSizeDepPlantHydProps!(currentSite, currentCohort, bc_in)
                            end

                            FreeMemory(nextc)
                        end
                    end

                    nextc = nextnextc
                end

                # Advance current cohort downward (unless it is already shortest).
                if currentCohort.shorter !== nothing
                    currentCohort = currentCohort.shorter
                else
                    break
                end
            end

            # Count the cohorts.
            nocohorts = 0
            cc = currentPatch.tallest
            while cc !== nothing
                nocohorts += 1
                cc = cc.shorter
            end

            if hlm_use_cohort_age_tracking[] == itrue
                if nocohorts > EDParams[].max_cohort_per_patch
                    iterate = 1
                    dynamic_size_fusion_tolerance *= 1.1
                    dynamic_age_fusion_tolerance  *= 1.1
                else
                    iterate = 0
                end
            else
                if nocohorts > EDParams[].max_cohort_per_patch
                    iterate = 1
                    dynamic_size_fusion_tolerance *= 1.1
                else
                    iterate = 0
                end
            end

            if dynamic_size_fusion_tolerance > 100.0
                fates_endrun("exceeded reasonable expectation of cohort fusion.")
            end
        end
    end

    if fusion_took_place == 1
        sort_cohorts(currentPatch)
    end

    return nothing
end

# Fuse the mass + scalar state of `nextc` (donor) into `currentCohort` (recipient).
# Carries out the carbon/crown-area-conserving dbh solve and all N-weighted means.
# `newn` is the (already-computed) combined number density.
function _fuse_pair!(currentSite::ed_site_type, currentPatch::fates_patch_type,
                     currentCohort::fates_cohort_type, nextc::fates_cohort_type,
                     newn::Real, bc_in)
    # new cohort age is the N-weighted mean of the two
    currentCohort.coage =
        (currentCohort.coage * (currentCohort.n / (currentCohort.n + nextc.n))) +
        (nextc.coage * (nextc.n / (currentCohort.n + nextc.n)))

    if hlm_use_cohort_age_tracking[] == itrue
        currentCohort.coage_class, currentCohort.coage_by_pft_class =
            coagetype_class_index(currentCohort.coage, currentCohort.pft)
    end

    # Fuse all mass pools (recipient weight = currentCohort.n / newn).
    WeightedFusePRTVartypes!(currentCohort.prt, nextc.prt, currentCohort.n / newn)

    # Leaf biophysical rates (leaf-mass weighting handled internally).
    UpdateCohortBioPhysRates(currentCohort)

    currentCohort.l2fr = (currentCohort.n * currentCohort.l2fr +
                          nextc.n * nextc.l2fr) / newn
    currentCohort.canopy_trim = (currentCohort.n * currentCohort.canopy_trim +
                                 nextc.n * nextc.canopy_trim) / newn

    # c13disc_acc: GPP-weighted mean, guarding the 0/0 case (quirk preserved).
    if (currentCohort.n * currentCohort.gpp_acc + nextc.n * nextc.gpp_acc) == 0.0
        currentCohort.c13disc_acc = 0.0
    else
        currentCohort.c13disc_acc =
            (currentCohort.n * currentCohort.gpp_acc * currentCohort.c13disc_acc +
             nextc.n * nextc.gpp_acc * nextc.c13disc_acc) /
            (currentCohort.n * currentCohort.gpp_acc + nextc.n * nextc.gpp_acc)
    end

    if cohort_fusion_conservation_method == conserve_crownarea_and_number_not_dbh
        # Conserve total crown area; solve dbh that conserves crown area + the
        # dbh->crown-area allometry. If above a capped allometry, fall back to the
        # weighted-mean dbh and ForceDBH.
        _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
            currentSite.spread, currentCohort.pft, currentCohort.crowndamage)
        _, nextc_carea = carea_allom(nextc.dbh, nextc.n, currentSite.spread,
            nextc.pft, nextc.crowndamage)
        currentCohort.c_area = currentCohort.c_area + nextc_carea

        dbh, _ = carea_allom(currentCohort.dbh, newn, currentSite.spread,
            currentCohort.pft, currentCohort.crowndamage, currentCohort.c_area;
            inverse = true)

        if abs(dbh - fates_unset_r8) < nearzero
            currentCohort.dbh = (currentCohort.n * currentCohort.dbh +
                                 nextc.n * nextc.dbh) / newn
            if prt_params.woody[currentCohort.pft] == itrue
                currentCohort.dbh, currentCohort.height = ForceDBH(currentCohort.pft,
                    currentCohort.crowndamage, currentCohort.canopy_trim,
                    currentCohort.efleaf_coh, currentCohort.efstem_coh,
                    currentCohort.dbh;
                    bdead = GetState(currentCohort.prt, struct_organ, carbon12_element))
            end
            _, currentCohort.c_area = carea_allom(currentCohort.dbh, newn,
                currentSite.spread, currentCohort.pft, currentCohort.crowndamage)
        else
            currentCohort.dbh = dbh
        end

        currentCohort.height, _ = h_allom(currentCohort.dbh, currentCohort.pft)

    elseif cohort_fusion_conservation_method == conserve_dbh_and_number_not_crownarea
        currentCohort.dbh = (currentCohort.n * currentCohort.dbh +
                             nextc.n * nextc.dbh) / newn
        currentCohort.height, _ = h_allom(currentCohort.dbh, currentCohort.pft)
        if prt_params.woody[currentCohort.pft] == itrue
            currentCohort.dbh, currentCohort.height = ForceDBH(currentCohort.pft,
                currentCohort.crowndamage, currentCohort.canopy_trim,
                currentCohort.efleaf_coh, currentCohort.efstem_coh,
                currentCohort.dbh;
                bdead = GetState(currentCohort.prt, struct_organ, carbon12_element))
        end
        _, currentCohort.c_area = carea_allom(currentCohort.dbh, newn,
            currentSite.spread, currentCohort.pft, currentCohort.crowndamage)
    else
        fates_endrun("FATES: Invalid choice for cohort_fusion_conservation_method")
    end

    currentCohort.size_class, currentCohort.size_by_pft_class =
        sizetype_class_index(currentCohort.dbh, currentCohort.pft)

    if hlm_use_planthydro[] == itrue
        # Conserve the fused cohort's plant water and recompute psi/ftc/btran.
        FuseCohortHydraulics(currentSite, currentCohort, nextc, bc_in, newn)
    end

    # recent canopy history
    currentCohort.canopy_layer_yesterday =
        (currentCohort.n * currentCohort.canopy_layer_yesterday +
         nextc.n * nextc.canopy_layer_yesterday) / newn

    # Size-class growth-flux diagnostic via fusion.
    if currentCohort.size_class_lasttimestep != nextc.size_class_lasttimestep
        if currentCohort.size_class_lasttimestep > nextc.size_class_lasttimestep
            largersc  = currentCohort.size_class_lasttimestep
            smallersc = nextc.size_class_lasttimestep
            larger_n  = currentCohort.n
            smaller_n = nextc.n
        else
            largersc  = nextc.size_class_lasttimestep
            smallersc = currentCohort.size_class_lasttimestep
            larger_n  = nextc.n
            smaller_n = currentCohort.n
        end

        # positive-growth case
        for sc_i in (smallersc + 1):currentCohort.size_class
            currentSite.growthflux_fusion[sc_i, currentCohort.pft] += smaller_n
        end
        # negative-growth case
        for sc_i in (currentCohort.size_class + 1):largersc
            currentSite.growthflux_fusion[sc_i, currentCohort.pft] -= larger_n
        end
        currentCohort.size_class_lasttimestep = currentCohort.size_class
    end

    # Flux + biophysics variables only meaningful for non-new cohorts.
    if !currentCohort.isnew
        currentCohort.seed_prod     = _wmean(currentCohort.n, currentCohort.seed_prod, nextc.n, nextc.seed_prod, newn)
        currentCohort.gpp_acc       = _wmean(currentCohort.n, currentCohort.gpp_acc, nextc.n, nextc.gpp_acc, newn)
        currentCohort.npp_acc       = _wmean(currentCohort.n, currentCohort.npp_acc, nextc.n, nextc.npp_acc, newn)
        currentCohort.resp_acc      = _wmean(currentCohort.n, currentCohort.resp_acc, nextc.n, nextc.resp_acc, newn)
        currentCohort.resp_acc_hold = _wmean(currentCohort.n, currentCohort.resp_acc_hold, nextc.n, nextc.resp_acc_hold, newn)
        currentCohort.npp_acc_hold  = _wmean(currentCohort.n, currentCohort.npp_acc_hold, nextc.n, nextc.npp_acc_hold, newn)
        currentCohort.gpp_acc_hold  = _wmean(currentCohort.n, currentCohort.gpp_acc_hold, nextc.n, nextc.gpp_acc_hold, newn)
        currentCohort.resp_excess   = _wmean(currentCohort.n, currentCohort.resp_excess, nextc.n, nextc.resp_excess, newn)
        currentCohort.dmort         = _wmean(currentCohort.n, currentCohort.dmort, nextc.n, nextc.dmort, newn)
        currentCohort.fire_mort     = _wmean(currentCohort.n, currentCohort.fire_mort, nextc.n, nextc.fire_mort, newn)

        # mortality diagnostics
        currentCohort.cmort  = _wmean(currentCohort.n, currentCohort.cmort, nextc.n, nextc.cmort, newn)
        currentCohort.hmort  = _wmean(currentCohort.n, currentCohort.hmort, nextc.n, nextc.hmort, newn)
        currentCohort.bmort  = _wmean(currentCohort.n, currentCohort.bmort, nextc.n, nextc.bmort, newn)
        currentCohort.smort  = _wmean(currentCohort.n, currentCohort.smort, nextc.n, nextc.smort, newn)
        currentCohort.asmort = _wmean(currentCohort.n, currentCohort.asmort, nextc.n, nextc.asmort, newn)
        currentCohort.frmort = _wmean(currentCohort.n, currentCohort.frmort, nextc.n, nextc.frmort, newn)

        # Nutrients (flexible-CNP only).
        if hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
            if nextc.n > currentCohort.n
                currentCohort.cnp_limiter = nextc.cnp_limiter
            end
            currentCohort.cx_int    = _wmean(currentCohort.n, currentCohort.cx_int, nextc.n, nextc.cx_int, newn)
            currentCohort.ema_dcxdt = _wmean(currentCohort.n, currentCohort.ema_dcxdt, nextc.n, nextc.ema_dcxdt, newn)
            currentCohort.cx0       = _wmean(currentCohort.n, currentCohort.cx0, nextc.n, nextc.cx0, newn)
            currentCohort.daily_nh4_uptake = _wmean(currentCohort.n, currentCohort.daily_nh4_uptake, nextc.n, nextc.daily_nh4_uptake, newn)
            currentCohort.daily_no3_uptake = _wmean(currentCohort.n, currentCohort.daily_no3_uptake, nextc.n, nextc.daily_no3_uptake, newn)
            currentCohort.sym_nfix_daily   = _wmean(currentCohort.n, currentCohort.sym_nfix_daily, nextc.n, nextc.sym_nfix_daily, newn)
            currentCohort.daily_n_gain     = _wmean(currentCohort.n, currentCohort.daily_n_gain, nextc.n, nextc.daily_n_gain, newn)
            currentCohort.daily_p_gain     = _wmean(currentCohort.n, currentCohort.daily_p_gain, nextc.n, nextc.daily_p_gain, newn)
            currentCohort.daily_p_demand   = _wmean(currentCohort.n, currentCohort.daily_p_demand, nextc.n, nextc.daily_p_demand, newn)
            currentCohort.daily_n_demand   = _wmean(currentCohort.n, currentCohort.daily_n_demand, nextc.n, nextc.daily_n_demand, newn)
            currentCohort.daily_c_efflux   = _wmean(currentCohort.n, currentCohort.daily_c_efflux, nextc.n, nextc.daily_c_efflux, newn)
            currentCohort.daily_n_efflux   = _wmean(currentCohort.n, currentCohort.daily_n_efflux, nextc.n, nextc.daily_n_efflux, newn)
            currentCohort.daily_p_efflux   = _wmean(currentCohort.n, currentCohort.daily_p_efflux, nextc.n, nextc.daily_p_efflux, newn)
        end

        # logging mortality
        currentCohort.lmort_direct     = _wmean(currentCohort.n, currentCohort.lmort_direct, nextc.n, nextc.lmort_direct, newn)
        currentCohort.lmort_collateral = _wmean(currentCohort.n, currentCohort.lmort_collateral, nextc.n, nextc.lmort_collateral, newn)
        currentCohort.lmort_infra      = _wmean(currentCohort.n, currentCohort.lmort_infra, nextc.n, nextc.lmort_infra, newn)
        currentCohort.l_degrad         = _wmean(currentCohort.n, currentCohort.l_degrad, nextc.n, nextc.l_degrad, newn)

        # biomass + dbh tendencies
        currentCohort.ddbhdt = _wmean(currentCohort.n, currentCohort.ddbhdt, nextc.n, nextc.ddbhdt, newn)

        # year_net_uptake: the 999 sentinel takes the MIN, not the weighted mean.
        for i in 1:nlevleaf
            if currentCohort.year_net_uptake[i] == 999.0 || nextc.year_net_uptake[i] == 999.0
                currentCohort.year_net_uptake[i] =
                    min(nextc.year_net_uptake[i], currentCohort.year_net_uptake[i])
            else
                currentCohort.year_net_uptake[i] =
                    (currentCohort.n * currentCohort.year_net_uptake[i] +
                     nextc.n * nextc.year_net_uptake[i]) / newn
            end
        end
    end

    return nothing
end

# N-weighted mean helper: (n1*a + n2*b)/newn.
@inline _wmean(n1::Real, a::Real, n2::Real, b::Real, newn::Real) = (n1 * a + n2 * b) / newn

# ===========================================================================
# sort_cohorts
# ===========================================================================

"""
    sort_cohorts(patchptr)

Re-sort the patch's cohorts into descending-height order by walking the existing
list and re-inserting every cohort into a freshly built list. Mirrors the Fortran
`sort_cohorts` (whose comment reads "DO NOT CHANGE THIS IT WILL BREAK").
"""
function sort_cohorts(patchptr::fates_patch_type)
    storebigcohort   = Ref{Union{fates_cohort_type,Nothing}}(nothing)
    storesmallcohort = Ref{Union{fates_cohort_type,Nothing}}(nothing)
    current_c = patchptr.tallest

    while current_c !== nothing
        next_c    = current_c.shorter
        tallestc  = storebigcohort[]
        shortestc = storesmallcohort[]

        if tallestc !== nothing
            tnull = 0
        else
            tnull = 1
            tallestc = current_c
        end

        if shortestc !== nothing
            snull = 0
        else
            snull = 1
            shortestc = current_c
        end

        insert_cohort(patchptr, current_c, tallestc, shortestc, tnull, snull,
                      storebigcohort, storesmallcohort)

        patchptr.tallest  = storebigcohort[]
        patchptr.shortest = storesmallcohort[]
        current_c = next_c
    end

    return nothing
end

# ===========================================================================
# insert_cohort
# ===========================================================================

"""
    insert_cohort(currentPatch, pcc, ptall, pshort, tnull, snull,
                  storebigcohort, storesmallcohort)

Insert cohort `pcc` into the height-ordered cohort linked list, finding the
position just below the next-taller cohort and updating the neighbors' links and
the patch tallest/shortest head/tail (returned via the `storebigcohort`/
`storesmallcohort` `Ref`s). `tnull`/`snull` flag that the supplied tallest/shortest
were actually null. Mirrors the Fortran `insert_cohort`.
"""
function insert_cohort(currentPatch::fates_patch_type, pcc::fates_cohort_type,
                       ptall::Union{fates_cohort_type,Nothing},
                       pshort::Union{fates_cohort_type,Nothing},
                       tnull::Integer, snull::Integer,
                       storebigcohort::Union{Ref,Nothing} = nothing,
                       storesmallcohort::Union{Ref,Nothing} = nothing)
    ptallest  = ptall
    pshortest = pshort

    if tnull == 1
        ptallest = nothing
    end
    if snull == 1
        pshortest = nothing
    end

    icohort = pcc
    tsp = icohort.height

    # Starting with the shortest tree, find the tree just taller than this one.
    current = pshortest
    exitloop = 0
    if current !== nothing
        while current !== nothing && exitloop == 0
            if current.height < tsp
                current = current.taller
            else
                exitloop = 1
            end
        end
    end

    if current !== nothing
        tallptr     = current
        tallptrnull = 0
    else
        tallptr     = nothing
        tallptrnull = 1
    end

    # new cohort is tallest
    if tallptr === nothing
        shortptr = ptallest          # new shorter = old tallest
        ptallest = icohort
        storebigcohort !== nothing && (storebigcohort[] = icohort)
        currentPatch.tallest = icohort
    else
        # next shorter = the next shorter to the cohort just taller
        shortptr = tallptr.shorter
        tallptr.shorter = icohort
    end

    # new cohort is shortest
    if shortptr === nothing
        pshortest = icohort
        storesmallcohort !== nothing && (storesmallcohort[] = icohort)
        currentPatch.shortest = icohort
    else
        shortptr.taller = icohort
    end

    # assign taller and shorter links for the new cohort
    icohort.taller = tallptr
    if tallptrnull == 1
        icohort.taller = nothing
    end
    icohort.shorter = shortptr

    return nothing
end

# ===========================================================================
# count_cohorts
# ===========================================================================

"""
    count_cohorts(currentPatch)

Count the cohorts on the patch (walking shortest -> taller) and store the result
in `currentPatch.countcohorts`. Mirrors the Fortran `count_cohorts` (the reverse
walk is a symmetry check; Fortran only errors under its `debug` flag — we keep the
check active since the linked list should always be symmetrical).
"""
function count_cohorts(currentPatch::fates_patch_type)
    currentCohort = currentPatch.shortest
    currentPatch.countcohorts = 0
    while currentCohort !== nothing
        currentPatch.countcohorts += 1
        currentCohort = currentCohort.taller
    end

    # Reverse-walk symmetry check.
    backcount = 0
    currentCohort = currentPatch.tallest
    while currentCohort !== nothing
        backcount += 1
        currentCohort = currentCohort.shorter
    end

    if backcount != currentPatch.countcohorts
        fates_endrun("problem with linked list, not symmetrical")
    end

    return nothing
end

# ===========================================================================
# EvaluateAndCorrectDBH
# ===========================================================================

"""
    EvaluateAndCorrectDBH(currentCohort) -> (delta_dbh, delta_height)

If a cohort's diameter is smaller than what is allometrically consistent with its
structural carbon (woody) or leaf carbon (grass), increase dbh (and height) to
match, returning the change. Mirrors the Fortran `EvaluateAndCorrectDBH` (whose
two `intent(out)` args are returned as a tuple).
"""
function EvaluateAndCorrectDBH(currentCohort::fates_cohort_type)
    dbh          = currentCohort.dbh
    ipft         = currentCohort.pft
    icrowndamage = currentCohort.crowndamage
    canopy_trim  = currentCohort.canopy_trim
    elongf_leaf  = currentCohort.efleaf_coh
    elongf_stem  = currentCohort.efstem_coh

    delta_dbh    = 0.0
    delta_height = 0.0

    if prt_params.woody[currentCohort.pft] == itrue
        struct_c = GetState(currentCohort.prt, struct_organ, carbon12_element)

        _, target_sapw_c, _ = bsap_allom(dbh, ipft, icrowndamage, canopy_trim, elongf_stem)
        target_agw_c, _     = bagw_allom(dbh, ipft, icrowndamage, elongf_stem)
        target_bgw_c, _     = bbgw_allom(dbh, ipft, elongf_stem)
        target_struct_c, _  = bdead_allom(target_agw_c, target_bgw_c, target_sapw_c, ipft)

        if (struct_c - target_struct_c) > calloc_abs_error
            dbh, height_out = ForceDBH(ipft, icrowndamage, canopy_trim,
                                       elongf_leaf, elongf_stem, dbh; bdead = struct_c)
            delta_dbh    = dbh - currentCohort.dbh
            delta_height = height_out - currentCohort.height
            currentCohort.dbh    = dbh
            currentCohort.height = height_out
        end
    else
        leaf_c = GetState(currentCohort.prt, leaf_organ, carbon12_element)
        target_leaf_c, _ = bleaf(dbh, ipft, icrowndamage, canopy_trim, elongf_leaf)

        if (leaf_c - target_leaf_c) > calloc_abs_error
            dbh, height_out = ForceDBH(ipft, icrowndamage, canopy_trim,
                                       elongf_leaf, elongf_stem, dbh; bl = leaf_c)
            delta_dbh    = dbh - currentCohort.dbh
            delta_height = height_out - currentCohort.height
            currentCohort.dbh    = dbh
            currentCohort.height = height_out
        end
    end

    return delta_dbh, delta_height
end

# ===========================================================================
# DamageRecovery
# ===========================================================================

"""
    DamageRecovery(csite, cpatch, ccohort) -> (newly_recovered, rcohort)

If a damaged cohort has enough excess resources, split off a fraction of its
plants into a NEW cohort with a lower (improved) damage class, inserting that
recovered cohort just taller than the donor in the list. Returns whether a new
cohort was created (and the new cohort, or `nothing`). Mirrors the Fortran
`DamageRecovery` (whose `newly_recovered` `intent(out)` is returned).
"""
function DamageRecovery(csite::ed_site_type, cpatch::fates_patch_type,
                        ccohort::fates_cohort_type)
    dbh         = ccohort.dbh
    ipft        = ccohort.pft
    canopy_trim = ccohort.canopy_trim
    elongf_leaf = ccohort.efleaf_coh
    elongf_fnrt = ccohort.effnrt_coh
    elongf_stem = ccohort.efstem_coh

    ev = edpftvarcon_inst()

    # Undamaged, or recovery disabled -> nothing to do.
    if ccohort.crowndamage == undamaged_class ||
       ev.damage_recovery_scalar[ipft] < nearzero
        return false, nothing
    end

    # Drought-deciduous + dormant -> cannot allocate to recovery; wait.
    is_hydecid_dormant =
        (prt_params.stress_decid[ipft] == ihard_stress_decid ||
         prt_params.stress_decid[ipft] == isemi_stress_decid) &&
        (ccohort.status_coh == leaves_off || ccohort.status_coh == leaves_shedding)
    is_sedecid_dormant =
        (prt_params.season_decid[ipft] == itrue) &&
        (ccohort.status_coh == leaves_off || ccohort.status_coh == leaves_shedding)

    if is_hydecid_dormant
        return false, nothing
    end

    # Target pools at the next (less) damage class.
    _, target_sapw_c, _ = bsap_allom(dbh, ipft, ccohort.crowndamage - 1, canopy_trim, elongf_stem)
    target_agw_c, _     = bagw_allom(dbh, ipft, ccohort.crowndamage - 1, elongf_stem)
    target_bgw_c, _     = bbgw_allom(dbh, ipft, elongf_stem)
    target_struct_c, _  = bdead_allom(target_agw_c, target_bgw_c, target_sapw_c, ipft)
    target_fnrt_c, _    = bfineroot(dbh, ipft, canopy_trim, ccohort.l2fr, elongf_fnrt)
    target_store_c, _   = bstore_allom(dbh, ipft, ccohort.crowndamage - 1, canopy_trim)
    target_leaf_c, _    = bleaf(dbh, ipft, ccohort.crowndamage - 1, canopy_trim, elongf_leaf)

    # Cold-deciduous + dormant: do not recover leaves (back-compat).
    if is_sedecid_dormant
        target_leaf_c = 0.0
    end

    nplant_recover = 1.0e10

    for el in 1:num_elements[]
        leaf_m   = GetState(ccohort.prt, leaf_organ, element_list[el])
        store_m  = GetState(ccohort.prt, store_organ, element_list[el])
        sapw_m   = GetState(ccohort.prt, sapw_organ, element_list[el])
        fnrt_m   = GetState(ccohort.prt, fnrt_organ, element_list[el])
        struct_m = GetState(ccohort.prt, struct_organ, element_list[el])
        repro_m  = GetState(ccohort.prt, repro_organ, element_list[el])

        elem = element_list[el]
        if elem == carbon12_element
            target_store_m  = target_store_c
            target_leaf_m   = target_leaf_c
            target_fnrt_m   = target_fnrt_c
            target_struct_m = target_struct_c
            target_sapw_m   = target_sapw_c
            available_m     = ccohort.npp_acc
        elseif elem == nitrogen_element
            target_struct_m = target_struct_c *
                prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[struct_organ]]
            target_leaf_m = target_leaf_c *
                prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[leaf_organ]]
            target_fnrt_m = target_fnrt_c *
                prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[fnrt_organ]]
            target_sapw_m = target_sapw_c *
                prt_params.nitr_stoich_p1[ipft, prt_params.organ_param_id[sapw_organ]]
            target_store_m = StorageNutrientTarget(ipft, elem,
                target_leaf_m, target_fnrt_m, target_sapw_m, target_struct_m)
            # uptake goes straight to storage: swap store -> available, zero store.
            available_m = store_m
            store_m     = 0.0
        else  # phosphorus_element
            target_struct_m = target_struct_c *
                prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[struct_organ]]
            target_leaf_m = target_leaf_c *
                prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[leaf_organ]]
            target_fnrt_m = target_fnrt_c *
                prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[fnrt_organ]]
            target_sapw_m = target_sapw_c *
                prt_params.phos_stoich_p1[ipft, prt_params.organ_param_id[sapw_organ]]
            target_store_m = StorageNutrientTarget(ipft, elem,
                target_leaf_m, target_fnrt_m, target_sapw_m, target_struct_m)
            available_m = store_m
            store_m     = 0.0
        end

        # mass at this damage class, and at the next-better class
        mass_d = leaf_m + store_m + sapw_m + fnrt_m + struct_m + repro_m
        mass_dminus1 = max(leaf_m, target_leaf_m) + max(fnrt_m, target_fnrt_m) +
                       max(store_m, target_store_m) + max(sapw_m, target_sapw_m) +
                       max(struct_m, target_struct_m)

        recovery_demand = mass_dminus1 - mass_d
        max_recover_nplant = available_m * ccohort.n / recovery_demand
        nplant_recover = min(nplant_recover,
            min(ccohort.n, max(0.0, max_recover_nplant * ev.damage_recovery_scalar[ipft])))
    end

    if nplant_recover < nearzero
        return false, nothing
    end

    # Build the recovered cohort as a copy of the donor.
    rcohort = fates_cohort_type()
    if hlm_use_planthydro[] == itrue
        InitHydrCohort(csite, rcohort)
    end
    rcohort.prt = InitPRTObject!()
    InitPRTBoundaryConditions(rcohort)
    Copy(ccohort, rcohort)

    rcohort.n           = nplant_recover
    rcohort.crowndamage = ccohort.crowndamage - 1

    # Adjust crown area (not per-individual).
    _, rcohort.c_area = carea_allom(dbh, rcohort.n, csite.spread, ipft, rcohort.crowndamage)

    # Update the un-recovered donor cohort.
    ccohort.n      = ccohort.n - rcohort.n
    ccohort.c_area = ccohort.c_area * ccohort.n / (ccohort.n + rcohort.n)

    # Insert the recovered cohort just taller than the donor.
    rcohort.shorter = ccohort
    if ccohort.taller !== nothing
        rcohort.taller         = ccohort.taller
        ccohort.taller.shorter = rcohort
    else
        cpatch.tallest = rcohort
        rcohort.taller = nothing
    end
    ccohort.taller = rcohort

    return true, rcohort
end
