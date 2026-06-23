# EDCanopyStructureMod.jl
# Julia port of FATES src/fates/biogeochem/EDCanopyStructureMod.F90 (Batch 15).
#
# The CANOPY STRUCTURE ENGINE. This module arranges cohorts into discrete canopy
# layers (1 = canopy/overstory, 2 = understory, ...) so that the total crown area
# in each layer fits the patch area, following the "Perfect Plasticity
# Approximation" (Purves et al. 2009; Fisher et al. 2010). The top organizer
# `canopy_structure!` iterates a demotion phase (move cohorts down from over-full
# layers via `DemoteFromLayer!`) and a promotion phase (move cohorts up into
# under-full layers via `PromoteIntoLayer!`) until each layer's crown area is
# within tolerance of the patch area. `canopy_spread!` updates the site-level
# crown-area scaling factor based on canopy closure. `canopy_summarization!` and
# `leaf_area_profile!` build the per-patch, canopy-layer x PFT x leaf-layer
# LAI/SAI/TAI profiles. `update_hlm_dynamics!` packs the FATES canopy state back to
# the host land model (LAI/SAI, canopy fraction, roughness, displacement, etc.).
# The remaining routines are area/LAI helpers.
#
# Translation notes (per project conventions):
#   * fates_r8 -> Float64; canopy-layer / leaf-layer indices -> Int.
#   * Fortran linked-list cohort traversal (`currentCohort => currentCohort%shorter`)
#     -> `while currentCohort !== nothing ... currentCohort = currentCohort.shorter`.
#   * Demotion/promotion cohort SPLITTING allocates a copy cohort, attaches a fresh
#     PARTEH object (InitPRTObject! + Copy + InitPRTBoundaryConditions!), and splices
#     it into the `taller`/`shorter` doubly-linked list exactly as the Fortran does.
#   * `carea_allom` returns `(dbh, c_area)` — we take only the c_area
#     (`_, c.c_area = carea_allom(...)`), mirroring FatesCohortMod's usage.
#   * Module-flag globals (hlm_use_planthydro/hlm_use_sp/hlm_use_cohort_age_tracking/
#     numpft) are Ref{Int} in FatesInterfaceTypesMod -> dereferenced with `[]`.
#   * ED params: nclmax/nlevleaf are `const`; ED_val_comp_excln/
#     ED_val_canopy_closure_thresh/dinc_vai/dlower_vai/radiation_model live on the
#     mutable `ed_params()` instance. AREA / N_HEIGHT_BINS / min_patch_area are
#     EDTypesMod consts (lowercase `area` etc.).
#   * Patch type-bound procedures are bang-functions: ReAllocateDynamics! /
#     NanDynamics! / ZeroDynamics!.
#   * The `preserve_b4b=true` Fortran flag (kept to preserve bit-for-bit baseline
#     test results vs. the two-stream code) is preserved here as a const; the b4b
#     leaf-area-profile path is the one that runs by default.
#
# STUBS (gated, inert by default):
#   * hlm_use_planthydro paths in DemoteFromLayer!/PromoteIntoLayer!/
#     update_hlm_dynamics! call InitHydrCohort/UpdateH2OVeg/RecruitWaterStorage,
#     which are not ported in FatesPlantHydraulicsMod. Guarded behind
#     `hlm_use_planthydro[] == itrue` (defaults to ifalse) — same convention as
#     Batches 12-14. See `# TODO Batch NN:` markers.
#
# Reused (NOT reimplemented) helpers: carea_allom/tree_lai/tree_sai/CrownDepth/
# VegAreaLayer (FatesAllometryMod); terminate_cohorts/terminate_cohort/
# fuse_cohorts/InitPRTObject! (EDCohortDynamicsMod); Copy/InitPRTBoundaryConditions!
# (FatesCohortMod); set_patchno (EDPatchDynamicsMod); sizetype_class_index/
# coagetype_class_index (FatesSizeAgeTypeIndicesMod); UpdateHarvestC!
# (EDLoggingMortalityMod); FatesConstructRadElements! (FatesTwoStreamUtilsMod);
# GetState + organ/element ids (PRTGenericMod); prt_params (PRTParametersMod);
# edpftvarcon_inst (EDPftvarcon).
# Standalone — NOT added to CLMInstances or any dual-copied struct.

# ---------------------------------------------------------------------------
# Module-local precision / tolerance parameters (Fortran module parameters)
# ---------------------------------------------------------------------------
const _ecs_area_target_precision    = 1.0e-11   # area conservation target
const _ecs_area_check_precision      = 1.0e-7    # absolute area-check tolerance
const _ecs_area_check_rel_precision  = 1.0e-4    # relative area-check tolerance
const _ecs_similar_height_tol        = 1.0e-3    # trees within 1mm treated as the same height
const _ecs_preserve_b4b              = true      # keep b4b leaf-area-profile path
const _ecs_canopy_debug              = false     # Fortran `debug` parameter
const _ecs_max_patch_iterations      = 10        # outer area-balance iteration cap

# ===========================================================================
# canopy_structure! — top-level organizer
# ===========================================================================

"""
    canopy_structure!(currentSite::ed_site_type, bc_in)

Allocate the `canopy_layer` attribute to each cohort so each canopy layer's total
crown area fits within the patch area. For every patch (oldest -> youngest) this
iteratively demotes cohorts from over-full layers (`DemoteFromLayer!`) and
promotes cohorts up into under-full layers (`PromoteIntoLayer!`), interleaving
cohort termination and fusion, until each layer's area is within tolerance of the
patch area (or `_ecs_max_patch_iterations` is exceeded). Stores the resulting
number of canopy layers in `patch.ncl_p`, and (under strict-PPA, when
`ED_val_comp_excln < 0`) the z* canopy-closure height. Mirrors the Fortran
`canopy_structure`.
"""
function canopy_structure!(currentSite::ed_site_type, bc_in)

    ED_val_comp_excln = ed_params().ED_val_comp_excln

    currentPatch = currentSite.oldest_patch

    # zero site-level demotion / promotion tracking info
    fill!(currentSite.demotion_rate, 0.0)
    fill!(currentSite.promotion_rate, 0.0)
    currentSite.demotion_carbonflux  = 0.0
    currentSite.promotion_carbonflux = 0.0

    # Section 1: Check total canopy area.
    while currentPatch !== nothing  # Patch loop

        # canopy layer has a special bounds check
        currentCohort = currentPatch.tallest
        while currentCohort !== nothing
            if currentCohort.canopy_layer < 1 || currentCohort.canopy_layer > nclmax + 1
                fates_endrun("BOGUS CANOPY LAYER: $(currentCohort.canopy_layer) " *
                             "(lat=$(currentSite.lat) lon=$(currentSite.lon))")
            end
            currentCohort = currentCohort.shorter
        end

        # Does any layer have excess area in it? Keep going until it does not...
        patch_area_counter = 0
        area_not_balanced  = true
        z = 1
        arealayer = zeros(Float64, nclmax + 2)

        while area_not_balanced

            # -----------------------------------------------------------------
            # Demotion Phase: identify over-full upper layers, demote downward.
            # -----------------------------------------------------------------

            # Its possible some cohort numbers are very low before we even enter
            # this scheme. Terminate them.
            terminate_cohorts(currentSite, currentPatch, 1, 12, bc_in)

            # Calculate how many layers we have in this canopy.
            z = NumPotentialCanopyLayers(currentPatch, currentSite.spread; include_substory=false)

            for i_lyr in 1:z   # Loop around the currently occupied canopy layers.
                DemoteFromLayer!(currentSite, currentPatch, i_lyr, bc_in)
            end

            # After demotions we may again have very sparse cohorts; remove them.
            terminate_cohorts(currentSite, currentPatch, 1, 13, bc_in)
            fuse_cohorts(currentSite, currentPatch, bc_in)
            terminate_cohorts(currentSite, currentPatch, 2, 13, bc_in)

            # -----------------------------------------------------------------
            # Promotion Phase: identify under-full upper layers, promote upward.
            # -----------------------------------------------------------------

            # Re-calculate number of layers without the false substory.
            z = NumPotentialCanopyLayers(currentPatch, currentSite.spread; include_substory=false)

            # We only promote if we have at least two layers.
            if z > 1
                for i_lyr in 1:(z - 1)
                    PromoteIntoLayer!(currentSite, currentPatch, i_lyr)
                end
                terminate_cohorts(currentSite, currentPatch, 1, 14, bc_in)
                fuse_cohorts(currentSite, currentPatch, bc_in)
                terminate_cohorts(currentSite, currentPatch, 2, 14, bc_in)
            end

            # -----------------------------------------------------------------
            # Check on layer area (if differences are not small, keep going).
            # -----------------------------------------------------------------
            z = NumPotentialCanopyLayers(currentPatch, currentSite.spread; include_substory=false)
            area_not_balanced = false
            for i_lyr in 1:z
                arealayer[i_lyr] = CanopyLayerArea(currentPatch, currentSite.spread, i_lyr)
                if ((arealayer[i_lyr] - currentPatch.area) / currentPatch.area > _ecs_area_check_rel_precision) ||
                   ((arealayer[i_lyr] - currentPatch.area) > _ecs_area_check_precision)
                    area_not_balanced = true
                end
            end

            # Gracefully exit if too many iterations have gone by.
            patch_area_counter += 1
            if patch_area_counter > _ecs_max_patch_iterations && area_not_balanced
                fates_endrun("PATCH AREA CHECK NOT CLOSING in canopy_structure! " *
                             "(patch area=$(currentPatch.area), spread=$(currentSite.spread), " *
                             "lat=$(currentSite.lat), lon=$(currentSite.lon))")
            end

        end  # while area_not_balanced

        # Save number of canopy layers to the patch structure.
        if z > nclmax
            fates_endrun("Termination should have ensured number of canopy layers " *
                         "was not larger than nclmax. Predicted: $(z), nclmax: $(nclmax). " *
                         "Consider increasing nclmax.")
        else
            currentPatch.ncl_p = z
        end

        # -------------------------------------------------------------------
        # Strict-PPA z* : height of the smallest tree in the canopy.
        # Loop top->bottom; find the shortest level-1 cohort whose shorter
        # neighbor is in level 2 and set zstar to that cohort's taller height.
        # -------------------------------------------------------------------
        if ED_val_comp_excln < 0.0
            currentPatch.zstar = 0.0
            currentCohort = currentPatch.tallest
            while currentCohort !== nothing
                if currentCohort.canopy_layer == 2
                    if currentCohort.taller !== nothing
                        if currentCohort.taller.canopy_layer == 1
                            currentPatch.zstar = currentCohort.taller.height
                        end
                    end
                end
                currentCohort = currentCohort.shorter
            end
        end

        currentPatch = currentPatch.younger
    end  # patch loop

    return nothing
end

# ===========================================================================
# DemoteFromLayer! — push cohorts down from an over-full layer
# ===========================================================================

"""
    DemoteFromLayer!(currentSite, currentPatch, i_lyr::Integer, bc_in)

If canopy layer `i_lyr` is over-occupied (total crown area exceeds patch area),
work out which cohorts to demote and move that excess crown area into layer
`i_lyr+1`. Demotion weight is either stochastic (inverse-height to the
`ED_val_comp_excln` power) when `ED_val_comp_excln >= 0`, or rank-ordered
deterministic (shortest-first, treating same-height cohorts as a group) when it is
negative. A cohort fully demoted just changes layer index; a partially demoted
cohort is split (the demoted part stays the original; an identical copy remains in
the upper layer). Mirrors the Fortran `DemoteFromLayer`.
"""
function DemoteFromLayer!(currentSite::ed_site_type, currentPatch::fates_patch_type,
                          i_lyr::Integer, bc_in)

    ED_val_comp_excln = ed_params().ED_val_comp_excln

    # First, determine how much total canopy area we have in this layer.
    arealayer = CanopyLayerArea(currentPatch, currentSite.spread, i_lyr)

    demote_area = arealayer - currentPatch.area

    if demote_area > _ecs_area_target_precision

        # This layer is over-occupied. Work out which cohorts to demote, going
        # from shortest to tallest for ranked demotion.

        sumweights = 0.0
        currentCohort = currentPatch.shortest
        while currentCohort !== nothing
            _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
                                                  currentSite.spread, currentCohort.pft,
                                                  currentCohort.crowndamage)

            if _ecs_canopy_debug && currentCohort.c_area < 0.0
                fates_endrun("negative c_area stage 1d in DemoteFromLayer!")
            end

            if currentCohort.canopy_layer == i_lyr

                if ED_val_comp_excln >= 0.0
                    # Stochastic: weight by inverse height to a constant power.
                    currentCohort.excl_weight = 1.0 / (currentCohort.height^ED_val_comp_excln)
                    sumweights += currentCohort.excl_weight
                else
                    # Rank-ordered deterministic. Same-height cohorts in the same
                    # layer are demoted as a single group.
                    total_crownarea_of_tied_cohorts = currentCohort.c_area
                    tied_size_with_neighbors = false
                    nextc = currentCohort.taller
                    while nextc !== nothing
                        if abs(nextc.height - currentCohort.height) < _ecs_similar_height_tol
                            if nextc.canopy_layer == currentCohort.canopy_layer
                                tied_size_with_neighbors = true
                                total_crownarea_of_tied_cohorts += nextc.c_area
                            end
                        else
                            break
                        end
                        nextc = nextc.taller
                    end

                    if tied_size_with_neighbors
                        currentCohort.excl_weight =
                            max(0.0, min(currentCohort.c_area,
                                (currentCohort.c_area / total_crownarea_of_tied_cohorts) *
                                (demote_area - sumweights)))
                        sumequal = currentCohort.excl_weight

                        nextc = currentCohort.taller
                        while nextc !== nothing
                            if abs(nextc.height - currentCohort.height) < _ecs_similar_height_tol
                                if nextc.canopy_layer == currentCohort.canopy_layer
                                    nextc.excl_weight =
                                        max(0.0, min(nextc.c_area,
                                            (nextc.c_area / total_crownarea_of_tied_cohorts) *
                                            (demote_area - sumweights)))
                                    sumequal += nextc.excl_weight
                                end
                            else
                                break
                            end
                            nextc = nextc.taller
                        end

                        # advance current pointer to the last similar cohort
                        if nextc !== nothing
                            currentCohort = nextc.shorter
                        else
                            currentCohort = currentPatch.tallest
                        end
                        sumweights += sumequal
                    else
                        currentCohort.excl_weight =
                            max(min(currentCohort.c_area, demote_area - sumweights), 0.0)
                        sumweights += currentCohort.excl_weight
                    end
                end
            end
            currentCohort = currentCohort.taller
        end

        # Stochastic case: normalize and pre-scale the demotion areas, capping
        # at each cohort's available area. See the tech note on promotion/demotion.
        if ED_val_comp_excln >= 0.0

            scale_factor_min = 1.0e10
            scale_factor     = 0.0
            currentCohort = currentPatch.tallest
            while currentCohort !== nothing
                if currentCohort.canopy_layer == i_lyr
                    currentCohort.excl_weight = currentCohort.excl_weight / sumweights
                    if 1.0 / currentCohort.excl_weight < scale_factor_min
                        scale_factor_min = 1.0 / currentCohort.excl_weight
                    end
                    scale_factor += currentCohort.excl_weight * currentCohort.c_area
                end
                currentCohort = currentCohort.shorter
            end

            # Factor to multiply demotion probabilities so the sum equals the
            # total amount to demote.
            scale_factor = demote_area / scale_factor

            if scale_factor <= scale_factor_min
                # Trivial case: all demotion fractions are less than 1.
                currentCohort = currentPatch.tallest
                while currentCohort !== nothing
                    if currentCohort.canopy_layer == i_lyr
                        currentCohort.excl_weight =
                            currentCohort.c_area * currentCohort.excl_weight * scale_factor
                    end
                    currentCohort = currentCohort.shorter
                end
            else
                # Non-trivial: at least one cohort's demotion would exceed its area.
                area_res         = 0.0
                scale_factor_res = 0.0
                currentCohort = currentPatch.tallest
                while currentCohort !== nothing
                    if currentCohort.canopy_layer == i_lyr
                        area_res += currentCohort.c_area * currentCohort.excl_weight * scale_factor_min
                        scale_factor_res += currentCohort.c_area *
                            (1.0 - (currentCohort.excl_weight * scale_factor_min))
                    end
                    currentCohort = currentCohort.shorter
                end

                area_res = demote_area - area_res
                scale_factor_res = area_res / scale_factor_res

                currentCohort = currentPatch.tallest
                while currentCohort !== nothing
                    if currentCohort.canopy_layer == i_lyr
                        currentCohort.excl_weight = currentCohort.c_area *
                            (currentCohort.excl_weight * scale_factor_min +
                             (1.0 - (currentCohort.excl_weight * scale_factor_min)) * scale_factor_res)
                    end
                    currentCohort = currentCohort.shorter
                end
            end
        end

        # Check the demotions meet the demand.
        sumweights = 0.0
        currentCohort = currentPatch.tallest
        while currentCohort !== nothing
            if currentCohort.canopy_layer == i_lyr
                sumweights += currentCohort.excl_weight
            end
            currentCohort = currentCohort.shorter
        end

        if abs(sumweights - demote_area) > _ecs_area_check_precision
            fates_endrun("demotions dont add up: sum=$(sumweights) needed=$(demote_area)")
        end

        # Weights have been calculated. Now move them to the lower layer.
        currentCohort = currentPatch.tallest
        while currentCohort !== nothing
            nextc = currentCohort.shorter

            if currentCohort.canopy_layer == i_lyr
                cc_loss  = currentCohort.excl_weight
                leaf_c   = GetState(currentCohort.prt, leaf_organ, carbon12_element)
                store_c  = GetState(currentCohort.prt, store_organ, carbon12_element)
                fnrt_c   = GetState(currentCohort.prt, fnrt_organ, carbon12_element)
                sapw_c   = GetState(currentCohort.prt, sapw_organ, carbon12_element)
                struct_c = GetState(currentCohort.prt, struct_organ, carbon12_element)

                if (cc_loss - currentCohort.c_area) > -nearzero &&
                   (cc_loss - currentCohort.c_area) < _ecs_area_target_precision

                    # Whole cohort demoted: just change its layer index.
                    currentCohort.canopy_layer = i_lyr + 1

                    currentSite.demotion_rate[currentCohort.size_class] += currentCohort.n
                    currentSite.demotion_carbonflux +=
                        (leaf_c + store_c + fnrt_c + sapw_c + struct_c) * currentCohort.n

                elseif (cc_loss < currentCohort.c_area) && (cc_loss > _ecs_area_target_precision)

                    # Partial demotion: split the cohort. The copy stays in the
                    # upper story; the original is demoted to the understory.
                    copyc = fates_cohort_type()
                    copyc.prt = InitPRTObject!()

                    if hlm_use_planthydro[] == itrue
                        InitHydrCohort(currentSite, copyc)
                    end

                    Copy(currentCohort, copyc)
                    InitPRTBoundaryConditions(copyc)

                    newarea = currentCohort.c_area - cc_loss
                    copyc.n = currentCohort.n * newarea / currentCohort.c_area
                    currentCohort.n = currentCohort.n - copyc.n

                    copyc.canopy_layer = i_lyr               # taller cohort is the copy
                    currentCohort.canopy_layer = i_lyr + 1   # demote the original

                    currentSite.demotion_rate[currentCohort.size_class] += currentCohort.n
                    currentSite.demotion_carbonflux +=
                        (leaf_c + store_c + fnrt_c + sapw_c + struct_c) * currentCohort.n

                    _, copyc.c_area = carea_allom(copyc.dbh, copyc.n, currentSite.spread,
                                                  copyc.pft, copyc.crowndamage)
                    _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
                                                          currentSite.spread, currentCohort.pft,
                                                          currentCohort.crowndamage)

                    # Insert the copy into the linked list (above currentCohort).
                    copyc.shorter = currentCohort
                    if currentCohort.taller !== nothing
                        copyc.taller = currentCohort.taller
                        currentCohort.taller.shorter = copyc
                    else
                        currentPatch.tallest = copyc
                        copyc.taller = nothing
                    end
                    currentCohort.taller = copyc

                elseif cc_loss > currentCohort.c_area
                    fates_endrun("more area than the cohort has is being demoted: " *
                                 "loss=$(cc_loss) existing=$(currentCohort.c_area)")
                end

                # Kill cohorts demoted into disallowed (too-deep) layers.
                if currentCohort.canopy_layer > nclmax
                    terminate_cohort(currentSite, currentPatch, currentCohort, bc_in,
                                     i_term_mort_type_canlev)
                    # Julia is GC'd; the cohort is unlinked inside terminate_cohort.
                else
                    _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
                                                          currentSite.spread, currentCohort.pft,
                                                          currentCohort.crowndamage)
                end
            end

            # We don't use the typical "point to shorter" here because the
            # current cohort may have been terminated.
            currentCohort = nextc
        end

        # Update the layer-area calculation for the current layer.
        arealayer = CanopyLayerArea(currentPatch, currentSite.spread, i_lyr)

        if (abs(arealayer - currentPatch.area) / arealayer > _ecs_area_check_rel_precision) ||
           (abs(arealayer - currentPatch.area) > _ecs_area_check_precision)
            fates_endrun("demotion did not trim area within tolerance: " *
                         "arealayer=$(arealayer) patch.area=$(currentPatch.area) ilayer=$(i_lyr)")
        end
    end

    return nothing
end

# ===========================================================================
# PromoteIntoLayer! — pull cohorts up into an under-full layer
# ===========================================================================

"""
    PromoteIntoLayer!(currentSite, currentPatch, i_lyr::Integer)

If canopy layer `i_lyr` is under-full (its crown area is less than the patch area),
promote crown area from the layer below (`i_lyr+1`). If the whole lower layer fits,
promote all of it; otherwise weight which cohorts to promote (stochastic by height
to the `ED_val_comp_excln` power, or rank-ordered deterministic treating tied-height
cohorts as a group), capping each cohort's promotion at its available area. A fully
promoted cohort changes layer index; a partially promoted one is split. Mirrors the
Fortran `PromoteIntoLayer`.
"""
function PromoteIntoLayer!(currentSite::ed_site_type, currentPatch::fates_patch_type,
                           i_lyr::Integer)

    ED_val_comp_excln = ed_params().ED_val_comp_excln

    arealayer_current = CanopyLayerArea(currentPatch, currentSite.spread, i_lyr)
    arealayer_below   = CanopyLayerArea(currentPatch, currentSite.spread, i_lyr + 1)

    # how much do we need to gain?
    promote_area = currentPatch.area - arealayer_current

    if promote_area > _ecs_area_target_precision

        if arealayer_below <= promote_area
            # Promote ALL cohorts from the layer below (it's smaller than the gain).
            currentCohort = currentPatch.tallest
            while currentCohort !== nothing
                if currentCohort.canopy_layer == i_lyr + 1
                    leaf_c   = GetState(currentCohort.prt, leaf_organ, carbon12_element)
                    store_c  = GetState(currentCohort.prt, store_organ, carbon12_element)
                    fnrt_c   = GetState(currentCohort.prt, fnrt_organ, carbon12_element)
                    sapw_c   = GetState(currentCohort.prt, sapw_organ, carbon12_element)
                    struct_c = GetState(currentCohort.prt, struct_organ, carbon12_element)

                    currentCohort.canopy_layer = i_lyr
                    _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
                                                          currentSite.spread, currentCohort.pft,
                                                          currentCohort.crowndamage)
                    currentSite.promotion_rate[currentCohort.size_class] += currentCohort.n
                    currentSite.promotion_carbonflux +=
                        (leaf_c + fnrt_c + store_c + sapw_c + struct_c) * currentCohort.n
                end
                currentCohort = currentCohort.shorter
            end

        else
            # Non-trivial: lower layer can accommodate more than necessary.
            # Figure out promotion weights (opposite of the demotion weighting).
            sumweights = 0.0
            currentCohort = currentPatch.tallest
            while currentCohort !== nothing
                _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
                                                      currentSite.spread, currentCohort.pft,
                                                      currentCohort.crowndamage)
                if currentCohort.canopy_layer == i_lyr + 1

                    if ED_val_comp_excln >= 0.0
                        currentCohort.prom_weight = currentCohort.height^ED_val_comp_excln
                        sumweights += currentCohort.prom_weight
                    else
                        # Rank-ordered deterministic; tied-height group handling.
                        total_crownarea_of_tied_cohorts = currentCohort.c_area
                        tied_size_with_neighbors = false
                        nextc = currentCohort.shorter
                        while nextc !== nothing
                            if abs(nextc.height - currentCohort.height) < _ecs_similar_height_tol
                                if nextc.canopy_layer == currentCohort.canopy_layer
                                    tied_size_with_neighbors = true
                                    total_crownarea_of_tied_cohorts += nextc.c_area
                                end
                            else
                                break
                            end
                            nextc = nextc.shorter
                        end

                        if tied_size_with_neighbors
                            currentCohort.prom_weight =
                                max(0.0, min(currentCohort.c_area,
                                    (currentCohort.c_area / total_crownarea_of_tied_cohorts) *
                                    (promote_area - sumweights)))
                            sumequal = currentCohort.prom_weight

                            nextc = currentCohort.shorter
                            while nextc !== nothing
                                if abs(nextc.height - currentCohort.height) < _ecs_similar_height_tol
                                    if nextc.canopy_layer == currentCohort.canopy_layer
                                        nextc.prom_weight =
                                            max(0.0, min(nextc.c_area,
                                                (nextc.c_area / total_crownarea_of_tied_cohorts) *
                                                (promote_area - sumweights)))
                                        sumequal += nextc.prom_weight
                                    end
                                else
                                    break
                                end
                                nextc = nextc.shorter
                            end

                            if nextc !== nothing
                                currentCohort = nextc.taller
                            else
                                currentCohort = currentPatch.shortest
                            end
                            sumweights += sumequal
                        else
                            currentCohort.prom_weight =
                                max(min(currentCohort.c_area, promote_area - sumweights), 0.0)
                            sumweights += currentCohort.prom_weight
                        end
                    end
                end
                currentCohort = currentCohort.shorter
            end

            # Stochastic case: normalize + pre-scale, capping at available area.
            if ED_val_comp_excln >= 0.0

                scale_factor_min = 1.0e10
                scale_factor     = 0.0
                currentCohort = currentPatch.tallest
                while currentCohort !== nothing
                    if currentCohort.canopy_layer == (i_lyr + 1)
                        currentCohort.prom_weight = currentCohort.prom_weight / sumweights
                        if 1.0 / currentCohort.prom_weight < scale_factor_min
                            scale_factor_min = 1.0 / currentCohort.prom_weight
                        end
                        scale_factor += currentCohort.prom_weight * currentCohort.c_area
                    end
                    currentCohort = currentCohort.shorter
                end

                scale_factor = promote_area / scale_factor

                if scale_factor <= scale_factor_min
                    currentCohort = currentPatch.tallest
                    while currentCohort !== nothing
                        if currentCohort.canopy_layer == (i_lyr + 1)
                            currentCohort.prom_weight =
                                currentCohort.c_area * currentCohort.prom_weight * scale_factor
                        end
                        currentCohort = currentCohort.shorter
                    end
                else
                    area_res         = 0.0
                    scale_factor_res = 0.0
                    currentCohort = currentPatch.tallest
                    while currentCohort !== nothing
                        if currentCohort.canopy_layer == (i_lyr + 1)
                            area_res += currentCohort.c_area * currentCohort.prom_weight * scale_factor_min
                            scale_factor_res += currentCohort.c_area *
                                (1.0 - (currentCohort.prom_weight * scale_factor_min))
                        end
                        currentCohort = currentCohort.shorter
                    end

                    area_res = promote_area - area_res
                    scale_factor_res = area_res / scale_factor_res

                    currentCohort = currentPatch.tallest
                    while currentCohort !== nothing
                        if currentCohort.canopy_layer == (i_lyr + 1)
                            currentCohort.prom_weight = currentCohort.c_area *
                                (currentCohort.prom_weight * scale_factor_min +
                                 (1.0 - (currentCohort.prom_weight * scale_factor_min)) * scale_factor_res)
                        end
                        currentCohort = currentCohort.shorter
                    end
                end
            end

            # Check the promotions meet the demand (debug-gated in Fortran).
            if _ecs_canopy_debug
                sumweights = 0.0
                currentCohort = currentPatch.tallest
                while currentCohort !== nothing
                    if currentCohort.canopy_layer == (i_lyr + 1)
                        sumweights += currentCohort.prom_weight
                    end
                    currentCohort = currentCohort.shorter
                end
                if abs(sumweights - promote_area) > _ecs_area_check_precision
                    fates_endrun("promotions dont add up: sum=$(sumweights) needed=$(promote_area)")
                end
            end

            currentCohort = currentPatch.tallest
            while currentCohort !== nothing
                if currentCohort.canopy_layer == i_lyr + 1

                    cc_gain  = currentCohort.prom_weight
                    leaf_c   = GetState(currentCohort.prt, leaf_organ, carbon12_element)
                    store_c  = GetState(currentCohort.prt, store_organ, carbon12_element)
                    fnrt_c   = GetState(currentCohort.prt, fnrt_organ, carbon12_element)
                    sapw_c   = GetState(currentCohort.prt, sapw_organ, carbon12_element)
                    struct_c = GetState(currentCohort.prt, struct_organ, carbon12_element)

                    if (cc_gain - currentCohort.c_area) > -nearzero &&
                       (cc_gain - currentCohort.c_area) < _ecs_area_target_precision

                        currentCohort.canopy_layer = i_lyr
                        currentSite.promotion_rate[currentCohort.size_class] += currentCohort.n
                        currentSite.promotion_carbonflux +=
                            (leaf_c + fnrt_c + store_c + sapw_c + struct_c) * currentCohort.n

                    elseif (cc_gain < currentCohort.c_area) && (cc_gain > _ecs_area_target_precision)

                        copyc = fates_cohort_type()
                        copyc.prt = InitPRTObject!()

                        if hlm_use_planthydro[] == itrue
                            InitHydrCohort(currentSite, copyc)
                        end

                        Copy(currentCohort, copyc)
                        InitPRTBoundaryConditions(copyc)

                        # new area of the existing (remaining) cohort
                        newarea = currentCohort.c_area - cc_gain

                        _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
                                                              currentSite.spread, currentCohort.pft,
                                                              currentCohort.crowndamage)

                        # number of individuals in the promoted copy
                        copyc.n = currentCohort.n * cc_gain / currentCohort.c_area
                        # number remaining in the understory
                        currentCohort.n = currentCohort.n - copyc.n

                        currentCohort.canopy_layer = i_lyr + 1   # remaining stays in understory
                        copyc.canopy_layer = i_lyr               # promote the copy

                        currentSite.promotion_rate[copyc.size_class] += copyc.n
                        currentSite.promotion_carbonflux +=
                            (leaf_c + fnrt_c + store_c + sapw_c + struct_c) * copyc.n

                        _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
                                                              currentSite.spread, currentCohort.pft,
                                                              currentCohort.crowndamage)
                        _, copyc.c_area = carea_allom(copyc.dbh, copyc.n, currentSite.spread,
                                                      copyc.pft, copyc.crowndamage)

                        # Insert the copy into the linked list (above currentCohort).
                        copyc.shorter = currentCohort
                        if currentCohort.taller !== nothing
                            copyc.taller = currentCohort.taller
                            currentCohort.taller.shorter = copyc
                        else
                            currentPatch.tallest = copyc
                            copyc.taller = nothing
                        end
                        currentCohort.taller = copyc

                    elseif cc_gain > currentCohort.c_area
                        fates_endrun("more area than the cohort has is being promoted: " *
                                     "gain=$(cc_gain) existing=$(currentCohort.c_area)")
                    end
                end
                currentCohort = currentCohort.shorter
            end

            arealayer_current = CanopyLayerArea(currentPatch, currentSite.spread, i_lyr)

            if (abs(arealayer_current - currentPatch.area) / arealayer_current > _ecs_area_check_rel_precision) ||
               (abs(arealayer_current - currentPatch.area) > _ecs_area_check_precision)
                fates_endrun("promotion did not bring area within tolerance: " *
                             "arealayer=$(arealayer_current) patch.area=$(currentPatch.area)")
            end
        end
    end

    return nothing
end

# ===========================================================================
# canopy_spread! — crown-area scaling factor
# ===========================================================================

"""
    canopy_spread!(currentSite::ed_site_type)

Update the site-level `spread` factor controlling the spatial spread of tree
canopies, based on canopy closure. The site-level top-layer woody crown area is
summed over all patches; if it exceeds `ED_val_canopy_closure_thresh` of the
notional area, `spread` is decremented (squashing canopies taller/thinner),
otherwise incremented, then clamped to [0,1]. Mirrors the Fortran `canopy_spread`.
"""
function canopy_spread!(currentSite::ed_site_type)
    ED_val_canopy_closure_thresh = ed_params().ED_val_canopy_closure_thresh

    inc = 0.05

    currentPatch = currentSite.oldest_patch
    sitelevel_canopyarea = 0.0
    while currentPatch !== nothing
        currentCohort = currentPatch.tallest
        while currentCohort !== nothing
            _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
                                                  currentSite.spread, currentCohort.pft,
                                                  currentCohort.crowndamage)
            if (prt_params.woody[currentCohort.pft] == itrue) &&
               (currentCohort.canopy_layer == 1)
                sitelevel_canopyarea += currentCohort.c_area
            end
            currentCohort = currentCohort.shorter
        end
        currentPatch = currentPatch.younger
    end

    # If the canopy area is approaching closure, squash the tree canopies and
    # make them taller and thinner.
    if sitelevel_canopyarea / area > ED_val_canopy_closure_thresh
        currentSite.spread = currentSite.spread - inc
    else
        currentSite.spread = currentSite.spread + inc
    end

    # Keep within [0,1].
    currentSite.spread = max(min(currentSite.spread, 1.0), 0.0)

    return nothing
end

# ===========================================================================
# canopy_summarization! — per-patch totals + leaf-area profile
# ===========================================================================

"""
    canopy_summarization!(nsites::Integer, sites::AbstractVector, bc_in::AbstractVector)

For each site: set the patch numbering (`set_patchno`), then for each patch
accumulate the total canopy area and total tree (woody) area over its top-layer
cohorts, updating each cohort's size/age class indices and crown area, before
building the per-patch leaf-area profile via [`leaf_area_profile!`](@ref) and
(when the two-stream radiation solver is active) the radiation elements. Mirrors
the Fortran `canopy_summarization`.
"""
function canopy_summarization!(nsites::Integer, sites::AbstractVector, bc_in::AbstractVector)

    radiation_model = ed_params().radiation_model

    for s in 1:nsites
        # Set patch indices (oldest = 1).
        set_patchno(sites[s])

        currentPatch = sites[s].oldest_patch
        while currentPatch !== nothing
            # zero cohort-summed variables
            currentPatch.total_canopy_area = 0.0
            currentPatch.total_tree_area   = 0.0

            currentCohort = currentPatch.shortest
            while currentCohort !== nothing
                ft       = currentCohort.pft
                leaf_c   = GetState(currentCohort.prt, leaf_organ, carbon12_element)
                sapw_c   = GetState(currentCohort.prt, sapw_organ, carbon12_element)
                struct_c = GetState(currentCohort.prt, struct_organ, carbon12_element)
                fnrt_c   = GetState(currentCohort.prt, fnrt_organ, carbon12_element)
                store_c  = GetState(currentCohort.prt, store_organ, carbon12_element)

                # Update size/SCPF classification.
                currentCohort.size_class, currentCohort.size_by_pft_class =
                    sizetype_class_index(currentCohort.dbh, currentCohort.pft)

                if hlm_use_cohort_age_tracking[] == itrue
                    currentCohort.coage_class, currentCohort.coage_by_pft_class =
                        coagetype_class_index(currentCohort.coage, currentCohort.pft)
                end

                if hlm_use_sp[] == ifalse
                    _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n,
                                                          sites[s].spread, currentCohort.pft,
                                                          currentCohort.crowndamage)
                end

                if currentCohort.canopy_layer == 1
                    currentPatch.total_canopy_area += currentCohort.c_area
                    if prt_params.woody[ft] == itrue
                        currentPatch.total_tree_area += currentCohort.c_area
                    end
                end

                # Checks for SP and NOCOMP modes.
                if currentPatch.nocomp_pft_label == nocomp_bareground
                    fates_endrun("cohorts in barepatch (canopy_summarization!)")
                end

                if hlm_use_sp[] == itrue
                    if currentPatch.tallest !== nothing && currentPatch.tallest.shorter !== nothing
                        fates_endrun("more than one cohort in SP mode (canopy_summarization!)")
                    end
                    if currentPatch.total_canopy_area - currentPatch.area > area_error_1
                        fates_endrun("too much canopy in summary (SP mode)")
                    end
                end

                # Erroneous-zero checks.
                if currentCohort.dbh <= 0.0 || currentCohort.n == 0.0
                    fates_endrun("FATES: dbh or n is zero in canopy_summarization: " *
                                 "dbh=$(currentCohort.dbh) n=$(currentCohort.n)")
                end
                if currentCohort.pft == 0 || currentCohort.canopy_trim <= 0.0
                    fates_endrun("FATES: PFT or trim is zero in canopy_summarization: " *
                                 "pft=$(currentCohort.pft) trim=$(currentCohort.canopy_trim)")
                end
                if (sapw_c + leaf_c + fnrt_c) <= 0.0
                    fates_endrun("FATES: alive biomass is zero in canopy_summarization")
                end

                currentCohort = currentCohort.taller
            end

            if currentPatch.total_canopy_area > currentPatch.area
                if currentPatch.total_canopy_area - currentPatch.area > 0.001
                    fates_endrun("FATES: canopy area bigger than area: " *
                                 "$(currentPatch.total_canopy_area) $(currentPatch.area)")
                end
                currentPatch.total_canopy_area = currentPatch.area
            end

            currentPatch = currentPatch.younger
        end

        leaf_area_profile!(sites[s])

        if radiation_model == twostr_solver
            FatesConstructRadElements!(sites[s], bc_in[s].fcansno_pa, bc_in[s].coszen_pa)
        end
    end

    return nothing
end

# ===========================================================================
# UpdateFatesAvgSnowDepth! — host snow depth -> FATES occlusion depth
# ===========================================================================

"""
    UpdateFatesAvgSnowDepth!(sites::AbstractVector, bc_in::AbstractVector)

Update the per-site snow depth used to occlude vegetation, as the host snow depth
weighted by the snow areal coverage fraction. Mirrors the Fortran
`UpdateFatesAvgSnowDepth`.
"""
function UpdateFatesAvgSnowDepth!(sites::AbstractVector, bc_in::AbstractVector)
    for s in 1:length(sites)
        sites[s].snow_depth = bc_in[s].snow_depth_si * bc_in[s].frac_sno_eff_si
    end
    return nothing
end

# ===========================================================================
# leaf_area_profile! — vertical/horizontal leaf+stem area distribution
# ===========================================================================

"""
    leaf_area_profile!(currentSite::ed_site_type)

Distribute leaf and stem area in vertical (leaf-layer) and horizontal
(canopy-layer x PFT) space for every patch in the site. Updates each patch's
`tlai_profile`/`elai_profile`/`tsai_profile`/`esai_profile` (leaf/stem area per
unit of the PFT's own radiative footprint), `canopy_area_profile` (fractional crown
area of each leaf layer), `nleaf`/`nrad` (vegetation-layer counts) and
`canopy_mask`. Cohort `treelai`/`treesai`/`NV` are refreshed via
[`UpdatePatchLAI!`](@ref). Uses the `preserve_b4b` integration path by default
(weighted leaf area by VAI-bin remainder), with snow occlusion of layers below the
snow surface. Mirrors the Fortran `leaf_area_profile`.
"""
function leaf_area_profile!(currentSite::ed_site_type)

    dinc_vai   = ed_params().dinc_vai
    dlower_vai = ed_params().dlower_vai

    cpatch = currentSite.oldest_patch
    while cpatch !== nothing

        fill!(cpatch.nleaf, 0)
        fill!(cpatch.canopy_layer_tlai, 0.0)

        # Updates cohort treelai/treesai/NV and the %nleaf array.
        UpdatePatchLAI!(cpatch)

        # Allocate / resize the large dynamic profile arrays, then NaN + zero them.
        ReAllocateDynamics!(cpatch)
        NanDynamics!(cpatch)
        ZeroDynamics!(cpatch)

        fill!(cpatch.canopy_mask, 0)

        # NRAD currently equals NLEAF (no snow filtering applied to it).
        cpatch.nrad .= cpatch.nleaf

        # It is remotely possible (deserts) to have no canopy area at all.
        if cpatch.total_canopy_area > nearzero

            currentCohort = cpatch.shortest
            while currentCohort !== nothing
                ft = currentCohort.pft
                cl = currentCohort.canopy_layer

                if _ecs_preserve_b4b
                    # preserve_b4b: round-off-stable baseline path (RGK 12-27-23).
                    lai = currentCohort.treelai * currentCohort.c_area / cpatch.total_canopy_area
                    sai = currentCohort.treesai * currentCohort.c_area / cpatch.total_canopy_area
                    if (currentCohort.treelai + currentCohort.treesai) > nearzero
                        # See fates issue #899: fleaf computed from layer lai/sai.
                        fleaf = lai / (lai + sai)
                    else
                        fleaf = 0.0
                    end

                    cpatch.nrad[cl, ft] = cpatch.nleaf[cl, ft]

                    if cpatch.nrad[cl, ft] > nlevleaf
                        fates_endrun("Number of radiative leaf layers larger than max: " *
                                     "cl=$(cl) ft=$(ft) nlevleaf=$(nlevleaf) nrad=$(cpatch.nrad[cl, ft])")
                    end

                    crown_depth = CrownDepth(currentCohort.height, currentCohort.pft)

                    for iv in 1:currentCohort.nv
                        layer_top_height = currentCohort.height -
                            (Float64(iv - 1) / currentCohort.nv * crown_depth)
                        layer_bottom_height = currentCohort.height -
                            (Float64(iv) / currentCohort.nv * crown_depth)

                        fraction_exposed = 1.0
                        if currentSite.snow_depth > layer_top_height
                            fraction_exposed = 0.0
                        end
                        if currentSite.snow_depth < layer_bottom_height
                            fraction_exposed = 1.0
                        end
                        if currentSite.snow_depth >= layer_bottom_height &&
                           currentSite.snow_depth <= layer_top_height
                            fraction_exposed = 1.0 - max(0.0, min(1.0,
                                (currentSite.snow_depth - layer_bottom_height) /
                                (layer_top_height - layer_bottom_height)))
                        end

                        if iv == currentCohort.nv
                            remainder = (currentCohort.treelai + currentCohort.treesai) -
                                        (dlower_vai[iv] - dinc_vai[iv])
                            if remainder > dinc_vai[iv]
                                fates_endrun("ED: issue with remainder: " *
                                             "treelai=$(currentCohort.treelai) treesai=$(currentCohort.treesai) " *
                                             "dinc=$(dinc_vai[iv]) NV=$(currentCohort.nv) remainder=$(remainder)")
                            end
                        else
                            remainder = dinc_vai[iv]
                        end

                        frac = currentCohort.c_area / cpatch.total_canopy_area

                        cpatch.tlai_profile[cl, ft, iv] += remainder * fleaf * frac
                        cpatch.elai_profile[cl, ft, iv] += remainder * fleaf * frac * fraction_exposed
                        cpatch.tsai_profile[cl, ft, iv] += remainder * (1.0 - fleaf) * frac
                        cpatch.esai_profile[cl, ft, iv] += remainder * (1.0 - fleaf) * frac * fraction_exposed
                        cpatch.canopy_area_profile[cl, ft, iv] += frac
                    end

                else  # non-b4b path (VegAreaLayer integration)
                    for iv in 1:currentCohort.nv
                        _, _, elai_layer, esai_layer, tlai_layer, tsai_layer =
                            VegAreaLayer(currentCohort.treelai, currentCohort.treesai,
                                         currentCohort.height, iv, currentCohort.nv,
                                         currentCohort.pft, currentSite.snow_depth)

                        frac = currentCohort.c_area / cpatch.total_canopy_area
                        cpatch.tlai_profile[cl, ft, iv] += tlai_layer * frac
                        cpatch.elai_profile[cl, ft, iv] += elai_layer * frac
                        cpatch.tsai_profile[cl, ft, iv] += tsai_layer * frac
                        cpatch.esai_profile[cl, ft, iv] += esai_layer * frac
                        cpatch.canopy_area_profile[cl, ft, iv] += frac
                    end
                end

                currentCohort = currentCohort.taller
            end

            # If there is an upper-story, the top canopy layer should have
            # exactly 1.0 in its top leaf layer.
            if (cpatch.ncl_p > 1) && (sum(@view cpatch.canopy_area_profile[1, :, 1]) < 0.9999)
                fates_endrun("FATES: canopy_area_profile was less than 1 at the canopy top: " *
                             "$(sum(@view cpatch.canopy_area_profile[1, :, 1]))")
            end

            # Normalize the area profiles from "per vegetated canopy area" into
            # "per area of the PFT's own radiative column". Also bounds-check.
            for cl in 1:cpatch.ncl_p
                # (debug-only > 1.0 check elided; preserved as Fortran behaviour)

                for ft in 1:numpft[]
                    for iv in 1:cpatch.nleaf[cl, ft]
                        if cpatch.canopy_area_profile[cl, ft, iv] > nearzero
                            cap = cpatch.canopy_area_profile[cl, ft, iv]
                            cpatch.tlai_profile[cl, ft, iv] /= cap
                            cpatch.tsai_profile[cl, ft, iv] /= cap
                            cpatch.elai_profile[cl, ft, iv] /= cap
                            cpatch.esai_profile[cl, ft, iv] /= cap
                        end
                    end
                end
            end

            # Set the radiation scattering-element mask.
            fill!(cpatch.canopy_mask, 0)
            if _ecs_preserve_b4b
                for cl in 1:cpatch.ncl_p
                    for ft in 1:numpft[]
                        for iv in 1:cpatch.nrad[cl, ft]
                            if cpatch.canopy_area_profile[cl, ft, iv] > 0.0
                                cpatch.canopy_mask[cl, ft] = 1
                            end
                        end
                    end
                end
            else
                for cl in 1:cpatch.ncl_p
                    for ft in 1:numpft[]
                        if cpatch.canopy_area_profile[cl, ft, 1] > 0.0
                            cpatch.canopy_mask[cl, ft] = 1
                        end
                    end
                end
            end
        end

        cpatch = cpatch.younger
    end

    return nothing
end

# ===========================================================================
# update_hlm_dynamics! — pack FATES canopy state to the host land model
# ===========================================================================

"""
    update_hlm_dynamics!(nsites::Integer, sites::AbstractVector, fcolumn::AbstractVector, bc_out::AbstractVector)

Package the vegetation-coverage output boundary conditions to the host land model:
per (vegetated) patch top-of-canopy / bottom-of-canopy height, crown-area-weighted
roughness length (`z0m`) and displacement height (`displa`), leaf/stem-area-weighted
leaf characteristic dimension (`dleaf`), patch canopy fraction, exposed/total LAI &
SAI via [`calc_areaindex`](@ref), the snow-free-vegetation flag, and the nocomp PFT
label. Applies a patch-area renormalization if the total drifts from 1.0. The
plant-hydraulics water-conservation/recruitment accounting is stubbed (off by
default). Mirrors the Fortran `update_hlm_dynamics`.
"""
function update_hlm_dynamics!(nsites::Integer, sites::AbstractVector,
                              fcolumn::AbstractVector, bc_out::AbstractVector)

    edpft = edpftvarcon_inst()

    for s in 1:nsites
        ifp = 0
        total_patch_area  = 0.0
        total_canopy_area = 0.0
        fill!(bc_out[s].canopy_fraction_pa, 0.0)
        fill!(bc_out[s].dleaf_pa, 0.0)
        fill!(bc_out[s].z0m_pa, 0.0)
        fill!(bc_out[s].displa_pa, 0.0)

        currentPatch = sites[s].oldest_patch
        while currentPatch !== nothing

            if currentPatch.nocomp_pft_label != nocomp_bareground
                # ignore the bare-ground-PFT patch entirely for these BC outs
                ifp += 1

                if currentPatch.total_canopy_area - currentPatch.area > 0.000001
                    currentPatch.total_canopy_area = currentPatch.area
                end

                if currentPatch.tallest !== nothing
                    bc_out[s].htop_pa[ifp] = currentPatch.tallest.height
                else
                    # FIX(RF,040113) — possible minimum-veg-height parameter?
                    bc_out[s].htop_pa[ifp] = 0.1
                end

                bc_out[s].hbot_pa[ifp] = max(0.0, min(0.2, bc_out[s].htop_pa[ifp] - 1.0))

                if currentPatch.total_canopy_area > nearzero
                    # Crown-area-weighted z0m / displa over top-layer cohorts.
                    currentCohort = currentPatch.shortest
                    while currentCohort !== nothing
                        if currentCohort.canopy_layer == 1
                            weight = min(1.0, currentCohort.c_area / currentPatch.total_canopy_area)
                            bc_out[s].z0m_pa[ifp] += edpft.z0mr[currentCohort.pft] *
                                                     currentCohort.height * weight
                            bc_out[s].displa_pa[ifp] += edpft.displar[currentCohort.pft] *
                                                        currentCohort.height * weight
                        end
                        currentCohort = currentCohort.taller
                    end

                    # dleaf weighted by total LAI+SAI in the patch.
                    total_patch_leaf_stem_area = 0.0
                    currentCohort = currentPatch.shortest
                    while currentCohort !== nothing
                        total_patch_leaf_stem_area +=
                            (currentCohort.treelai + currentCohort.treesai) * currentCohort.c_area
                        currentCohort = currentCohort.taller
                    end

                    if total_patch_leaf_stem_area > nearzero
                        currentCohort = currentPatch.shortest
                        while currentCohort !== nothing
                            weight = (currentCohort.treelai + currentCohort.treesai) *
                                     currentCohort.c_area / total_patch_leaf_stem_area
                            bc_out[s].dleaf_pa[ifp] += edpft.dleaf[currentCohort.pft] * weight
                            currentCohort = currentCohort.taller
                        end
                    else
                        bc_out[s].dleaf_pa[ifp] = edpft.dleaf[1]
                    end
                else
                    # No canopy: dummy (first PFT) aerodynamic properties.
                    bc_out[s].z0m_pa[ifp]    = edpft.z0mr[1] * bc_out[s].htop_pa[ifp]
                    bc_out[s].displa_pa[ifp] = edpft.displar[1] * bc_out[s].htop_pa[ifp]
                    bc_out[s].dleaf_pa[ifp]  = edpft.dleaf[1]
                end

                # Grass assumed to be located under tree canopies.
                if currentPatch.area > 0.0
                    bc_out[s].canopy_fraction_pa[ifp] =
                        min(1.0, currentPatch.total_canopy_area / currentPatch.area) *
                        (currentPatch.area / area)
                else
                    bc_out[s].canopy_fraction_pa[ifp] = 0.0
                end

                bare_frac_area = (1.0 - min(1.0, currentPatch.total_canopy_area / currentPatch.area)) *
                                 (currentPatch.area / area)

                total_patch_area  += bc_out[s].canopy_fraction_pa[ifp] + bare_frac_area
                total_canopy_area += bc_out[s].canopy_fraction_pa[ifp]

                bc_out[s].nocomp_pft_label_pa[ifp] = currentPatch.nocomp_pft_label

                # Area indices for output (profiles assumed already updated).
                bc_out[s].elai_pa[ifp] = calc_areaindex(currentPatch, "elai")
                bc_out[s].tlai_pa[ifp] = calc_areaindex(currentPatch, "tlai")
                bc_out[s].esai_pa[ifp] = calc_areaindex(currentPatch, "esai")
                bc_out[s].tsai_pa[ifp] = calc_areaindex(currentPatch, "tsai")

                # Snow-free vegetation flag (host uses it to gate photosynthesis).
                if (bc_out[s].elai_pa[ifp] + bc_out[s].esai_pa[ifp]) > 0.0
                    bc_out[s].frac_veg_nosno_alb_pa[ifp] = 1.0
                else
                    bc_out[s].frac_veg_nosno_alb_pa[ifp] = 0.0
                end

            else  # nocomp or SP bareground patch
                total_patch_area += currentPatch.area / area
            end
            currentPatch = currentPatch.younger
        end

        # Patch / canopy area corrections.
        if abs(total_patch_area - 1.0) > rsnbl_math_prec
            if abs(total_patch_area - 1.0) > 1.0e-8
                fates_endrun("total area is wrong in update_hlm_dynamics!: $(total_patch_area)")
            end

            currentPatch = sites[s].oldest_patch
            ifp = 0
            while currentPatch !== nothing
                if currentPatch.nocomp_pft_label != nocomp_bareground
                    ifp += 1
                    bc_out[s].canopy_fraction_pa[ifp] =
                        bc_out[s].canopy_fraction_pa[ifp] / total_patch_area
                end
                currentPatch = currentPatch.younger
            end
        end

        if hlm_use_planthydro[] == itrue
            # TODO Batch NN: UpdateH2OVeg not yet ported (FatesPlantHydraulicsMod).
        end

        # Pass FATES harvested C to bc_out.
        UpdateHarvestC!(sites[s], bc_out[s])
    end

    # Diagnose the site-level recruit water storage pool (no mass moved here).
    if hlm_use_planthydro[] == itrue
        RecruitWaterStorage(nsites, sites, bc_out)
    end

    return nothing
end

# ===========================================================================
# calc_areaindex — patch exposed/total LAI/SAI
# ===========================================================================

"""
    calc_areaindex(cpatch::fates_patch_type, ai_type::AbstractString) -> ai

Compute a patch-level area index (`"elai"`/`"tlai"`/`"esai"`/`"tsai"`) by
integrating the per-(canopy-layer x PFT x leaf-layer) area profile weighted by the
crown-area fraction. The result is m2 of leaf/stem per m2 of ground. A legacy
`ai_min = 0.1` floor is applied (kept for b4b; flagged for removal upstream).
Mirrors the Fortran `calc_areaindex`.
"""
function calc_areaindex(cpatch::fates_patch_type, ai_type::AbstractString)
    # TODO: ai_min is a long-ago testing artifact, kept for b4b (see Fortran note).
    ai_min = 0.1

    ai = 0.0
    if ai_type == "elai"
        prof = cpatch.elai_profile
    elseif ai_type == "tlai"
        prof = cpatch.tlai_profile
    elseif ai_type == "esai"
        prof = cpatch.esai_profile
    elseif ai_type == "tsai"
        prof = cpatch.tsai_profile
    else
        fates_endrun("Unsupported area index sent to calc_areaindex: $(ai_type)")
    end

    for cl in 1:cpatch.ncl_p
        for ft in 1:numpft[]
            nr = cpatch.nrad[cl, ft]
            for iv in 1:nr
                ai += cpatch.canopy_area_profile[cl, ft, iv] * prof[cl, ft, iv]
            end
        end
    end

    return max(ai_min, ai)
end

# ===========================================================================
# CanopyLayerArea — total crown area of one canopy layer
# ===========================================================================

"""
    CanopyLayerArea(currentPatch::fates_patch_type, site_spread::Real, layer_index::Integer) -> layer_area

Total crown-area footprint (m2, same units as `patch.area`) of all cohorts in
canopy layer `layer_index`. Recomputes each cohort's `c_area` from allometry as a
side effect. Mirrors the Fortran `CanopyLayerArea` (returned via an out-argument
there; returned directly here).
"""
function CanopyLayerArea(currentPatch::fates_patch_type, site_spread::Real,
                         layer_index::Integer)
    layer_area = 0.0
    currentCohort = currentPatch.tallest
    while currentCohort !== nothing
        _, currentCohort.c_area = carea_allom(currentCohort.dbh, currentCohort.n, site_spread,
                                              currentCohort.pft, currentCohort.crowndamage)
        if currentCohort.canopy_layer == layer_index
            layer_area += currentCohort.c_area
        end
        currentCohort = currentCohort.shorter
    end
    return layer_area
end

# ===========================================================================
# UpdatePatchLAI! — per-patch canopy-layer LAI accumulation
# ===========================================================================

"""
    UpdatePatchLAI!(currentPatch::fates_patch_type)

Walk the patch by canopy layer (top down) and, for each cohort, refresh its LAI/SAI
via [`UpdateCohortLAI!`](@ref), update the patch `nleaf(cl,ft)` vegetation-layer
count, and accumulate the layer LAI (`canopy_layer_tlai`). Because some understory
cohorts can be taller than top-layer cohorts, this must iterate canopy layer first
so above-layer LAI is known before a cohort's tree_lai is computed. Mirrors the
Fortran `UpdatePatchLAI`.
"""
function UpdatePatchLAI!(currentPatch::fates_patch_type)
    for cl in 1:nclmax
        currentCohort = currentPatch.tallest
        while currentCohort !== nothing
            if currentCohort.canopy_layer == cl
                ft = currentCohort.pft
                UpdateCohortLAI!(currentCohort, currentPatch.canopy_layer_tlai,
                                 currentPatch.total_canopy_area)

                currentPatch.nleaf[cl, ft] = max(currentPatch.nleaf[cl, ft], currentCohort.nv)

                currentPatch.canopy_layer_tlai[cl] += currentCohort.treelai *
                    currentCohort.c_area / currentPatch.total_canopy_area
            end
            currentCohort = currentCohort.shorter
        end
    end
    return nothing
end

# ===========================================================================
# UpdateCohortLAI! — per-cohort treelai/treesai/NV
# ===========================================================================

"""
    UpdateCohortLAI!(currentCohort::fates_cohort_type, canopy_layer_tlai::AbstractVector, total_canopy_area::Real)

Update a single cohort's `treelai` (per `tree_lai`), `treesai` (per `tree_sai`,
unless SP mode), and number of leaf layers `nv`. `canopy_layer_tlai` is the
per-layer total LAI of the layers above (used to attenuate LAI). Mirrors the
Fortran `UpdateCohortLAI`.
"""
function UpdateCohortLAI!(currentCohort::fates_cohort_type, canopy_layer_tlai::AbstractVector,
                          total_canopy_area::Real)
    dlower_vai = ed_params().dlower_vai

    leaf_c = GetState(currentCohort.prt, leaf_organ, carbon12_element)

    # tree_lai has an internal check on the canopy location.
    currentCohort.treelai = tree_lai(leaf_c, currentCohort.pft, currentCohort.c_area,
                                     currentCohort.n, currentCohort.canopy_layer,
                                     canopy_layer_tlai, currentCohort.vcmax25top)

    if hlm_use_sp[] == ifalse
        currentCohort.treesai = tree_sai(currentCohort.pft, currentCohort.dbh,
                                         currentCohort.crowndamage, currentCohort.canopy_trim,
                                         currentCohort.efstem_coh, currentCohort.c_area,
                                         currentCohort.n, currentCohort.canopy_layer,
                                         canopy_layer_tlai, currentCohort.treelai,
                                         currentCohort.vcmax25top, 4)
    end

    # Number of actual vegetation layers in this cohort's crown.
    currentCohort.nv =
        count(x -> (currentCohort.treelai + currentCohort.treesai) > x, dlower_vai) + 1

    return nothing
end

# ===========================================================================
# NumPotentialCanopyLayers — count canopy layers in a patch
# ===========================================================================

"""
    NumPotentialCanopyLayers(currentPatch::fates_patch_type, site_spread::Real; include_substory::Bool=false) -> z

Number of canopy layers in a patch, taken as the maximum `canopy_layer` over its
cohorts. If `include_substory`, also account for a temporary sub-understory layer
needed to hold demotions when the bottom layer's crown area already exceeds the
patch area (`z` is incremented in that case). Mirrors the Fortran
`NumPotentialCanopyLayers`.
"""
function NumPotentialCanopyLayers(currentPatch::fates_patch_type, site_spread::Real;
                                  include_substory::Bool=false)
    z = 1
    currentCohort = currentPatch.tallest
    while currentCohort !== nothing
        z = max(z, currentCohort.canopy_layer)
        currentCohort = currentCohort.shorter
    end

    if include_substory
        arealayer = 0.0
        currentCohort = currentPatch.tallest
        while currentCohort !== nothing
            if currentCohort.canopy_layer == z
                _, c_area = carea_allom(currentCohort.dbh, currentCohort.n, site_spread,
                                        currentCohort.pft, currentCohort.crowndamage)
                arealayer += c_area
            end
            currentCohort = currentCohort.shorter
        end

        # Does the bottom layer have more than a full canopy? Make another layer.
        if arealayer > currentPatch.area
            z = z + 1
        end
    end

    return z
end
