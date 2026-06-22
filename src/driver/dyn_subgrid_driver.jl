# ==========================================================================
# Ported from: src/dyn_subgrid/dynSubgridDriverMod.F90 (400 lines)
#
# Top-level orchestrator for the CLM transient-land-use subsystem. Three public
# routines (preserving the Fortran names):
#
#   - dynSubgrid_init!                  — initialize transient land cover: call the
#       per-dataset init routines for whichever transient aspects are enabled, build
#       the prior-weights + state-updater + conservation objects, do the first
#       interp (cold-start / init_interp initial weights), and wrap up the weight
#       changes.
#   - dynSubgrid_driver!                — per-timestep top-level driver: snapshot
#       prior weights + heat/water baselines, call each enabled interp reader to
#       update weights for the year, reconcile landunit areas + higher-order weights
#       + active flags, run the patch/column state updaters + biogeophys/biogeochem
#       conservation, and initialize newly-active columns.
#   - dynSubgrid_wrapup_weight_changes! — reconcile derived variables after subgrid
#       weights change (recompute higher-order weights, set_active / re-make filters).
#
# All the per-dataset readers, state updaters, conservation routines, prior-weights
# and init-columns machinery ALREADY EXIST on the base (see src/CLM.jl) and are
# reused verbatim here — this file only orchestrates them.
#
# DELIBERATE OMISSIONS / GATING (matching the conventions in the sibling dyn
# modules and the Fortran's optional paths; clearly flagged):
#   - FATES (use_fates / use_fates_luh, dyn_ED, dynFatesLandUseInterp) — DEFERRED;
#     gated off / no-op'd. The base has no FATES port.
#   - CNDV (bgc_vegetation_inst%UpdateSubgridWeights) — DEFERRED.
#   - The Fortran calls `reweight_wrapup` (reweightMod) which rebuilds filters and
#     sets active flags. There is no `reweight_wrapup` port on the base; the
#     available analogue is `set_active!` (subgrid_weights.jl), which re-derives the
#     landunit/column/patch active flags — i.e. the "filters are re-made" step. We
#     call that.
#   - set_subgrid_diagnostic_fields (history diagnostics) — no port; skipped.
#   - The biogeophys conservation (dyn_hwcontent_init!/dyn_hwcontent_final!) and
#     biogeochem conservation (dyn_cnbal_*) are OPTIONAL on this driver: they are
#     run only if the corresponding state objects are supplied (so the driver can
#     be exercised standalone with just the weight machinery). When the full state
#     bundle is provided they run, matching the Fortran's unconditional water/energy
#     accounting + use_cn-gated C/N accounting.
#
# As required, the orchestration state is bundled into a single standalone
# `DynSubgridState` @kwdef mutable struct. It is NOT added to CLMInstances or any
# ForwardDiff-dual-copied struct.
# ==========================================================================

"""
    DynSubgridState

Bundles all of the per-run orchestration objects for the transient-land-use
subsystem (the module-private singletons in `dynSubgridDriverMod.F90` plus the
per-dataset reader objects that the Fortran keeps in the individual file modules):

  - `ctl`                  : `DynSubgridControl` (which transient aspects are on)
  - `prior_weights`        : `PriorWeights`        (pre-update weight snapshot)
  - `patch_state_updater`  : `PatchStateUpdater`   (conservative patch-state update)
  - `column_state_updater` : `ColumnStateUpdater`  (conservative column-state update)
  - `dynbal`               : `DynConsBiogeophysState` (heat/water conservation)
  - per-dataset reader state objects (only those that are enabled are non-`nothing`):
      `dynpft_state` / `dyncrop_state` / `dynlake_state` / `dynurban_state` /
      `dynharvest_state` / `dyngrossunrep_state`

Standalone holder — NOT part of `CLMInstances` or any dual-copied struct.
"""
Base.@kwdef mutable struct DynSubgridState
    ctl::DynSubgridControl
    prior_weights::PriorWeights
    patch_state_updater::PatchStateUpdater
    column_state_updater::ColumnStateUpdater
    dynbal::DynConsBiogeophysState

    # Per-dataset readers (nothing when the corresponding transient aspect is off).
    dynpft_state       = nothing   # ::Union{DynpftState,Nothing}
    dyncrop_state      = nothing   # ::Union{DyncropState,Nothing}
    dynlake_state      = nothing   # ::Union{DynLakeFile,Nothing}
    dynurban_state     = nothing   # ::Union{DynUrbanFile,Nothing}
    dynharvest_state   = nothing   # ::Union{DynHarvestState,Nothing}
    dyngrossunrep_state = nothing  # ::Union{DynGrossUnrepState,Nothing}
end

# --------------------------------------------------------------------------
# dynSubgrid_init!
# --------------------------------------------------------------------------
"""
    dynSubgrid_init!(bounds_proc, ctl, grc, lun, col, pch;
                     nclumps=1, current_year,
                     natpft_size=0, cft_size=0,
                     wt_nat_patch=nothing, check_dynpft_consistency=true,
                     use_crop=false, crop_inst=nothing,
                     ng=bounds_proc.endg) -> DynSubgridState

Initialize objects needed for prescribed transient PFTs and/or dynamic landunits
(Fortran `dynSubgrid_init`).

Builds the `PriorWeights`, `PatchStateUpdater`, `ColumnStateUpdater` and
`DynConsBiogeophysState` objects (all sized to the PROC-level `bounds_proc`), then
calls the per-dataset init routines for whichever transient aspects are enabled in
`ctl`, sets the initial subgrid weights from file (the first `*_interp` call — for
cold start / init_interp), and finally wraps up the weight changes once.

`current_year` is the model year for the first interp. `natpft_size` / `cft_size`
are the file dimension sizes (required when the corresponding transient aspect is
on). All transient aspects default off (per `ctl`).

Note (matching Fortran): `dynHarvest_interp` is NOT called here — harvest data
aren't needed until the run loop and have nothing to do with subgrid weights.
"""
function dynSubgrid_init!(bounds_proc::BoundsType, ctl::DynSubgridControl,
                          grc::GridcellData, lun::LandunitData,
                          col::ColumnData, pch::PatchData;
                          nclumps::Int = 1,
                          current_year::Int,
                          natpft_size::Int = 0,
                          cft_size::Int = 0,
                          wt_nat_patch = nothing,
                          check_dynpft_consistency::Bool = true,
                          use_crop::Bool = false,
                          crop_inst = nothing,
                          collapse_crops::Bool = false,
                          glc_behavior = nothing,
                          ng::Int = bounds_proc.endg)

    @assert bounds_proc.level == BOUNDS_LEVEL_PROC "dynSubgrid_init!: argument must be PROC-level bounds"

    flanduse = get_flanduse_timeseries(ctl)

    # Module-private singleton objects (dynSubgridDriverMod private types).
    prior_weights        = PriorWeights(bounds_proc)
    patch_state_updater  = PatchStateUpdater(bounds_proc)
    column_state_updater = ColumnStateUpdater(bounds_proc, nclumps)
    dynbal               = dyn_cons_biogeophys_state_init(ng)

    state = DynSubgridState(
        ctl = ctl,
        prior_weights = prior_weights,
        patch_state_updater = patch_state_updater,
        column_state_updater = column_state_updater,
        dynbal = dynbal,
    )

    # ------------------------------------------------------------------
    # Initialize the per-dataset readers (Fortran *_init calls).
    # ------------------------------------------------------------------
    if get_do_transient_pfts(ctl)
        state.dynpft_state = dynpft_init(flanduse;
            ngridcells = ng, natpft_size = natpft_size, current_year = current_year,
            wt_nat_patch = wt_nat_patch,
            check_dynpft_consistency = check_dynpft_consistency)
    end

    if get_do_transient_crops(ctl)
        state.dyncrop_state = dyncrop_init(flanduse;
            ngridcells = ng, cft_size = cft_size, current_year = current_year)
    end

    # Harvest. (Currently the harvest data live on the flanduse_timeseries file.)
    if get_do_harvest(ctl)
        state.dynharvest_state = dynHarvest_init(bounds_proc.begg, bounds_proc.endg,
            flanduse; current_year = current_year)
    end

    # Gross unrepresented landuse change.
    if get_do_grossunrep(ctl)
        state.dyngrossunrep_state = dynGrossUnrep_init(bounds_proc.begg,
            bounds_proc.endg, natpft_size, flanduse; current_year = current_year)
    end

    # Prescribed transient lakes.
    if get_do_transient_lakes(ctl)
        state.dynlake_state = dynlake_init(flanduse, bounds_proc.begg,
            bounds_proc.endg; current_year = current_year)
    end

    # Prescribed transient urban.
    if get_do_transient_urban(ctl)
        state.dynurban_state = dynurban_init(flanduse, bounds_proc.begg,
            bounds_proc.endg; current_year = current_year)
    end

    # ------------------------------------------------------------------
    # Set initial subgrid weights for aspects read from file (cold start /
    # use_init_interp). These will be overwritten on a restart run.
    # ------------------------------------------------------------------
    if get_do_transient_pfts(ctl)
        dynpft_interp!(state.dynpft_state, bounds_proc, pch, lun; year = current_year)
    end

    if get_do_transient_crops(ctl)
        dyncrop_interp!(state.dyncrop_state, bounds_proc, grc, lun, col, pch;
            year = current_year, use_crop = use_crop, crop_inst = crop_inst,
            collapse_crops = collapse_crops)
    end

    if get_do_transient_lakes(ctl)
        dynlake_interp(state.dynlake_state, grc, lun, bounds_proc.begg,
            bounds_proc.endg, current_year)
    end

    if get_do_transient_urban(ctl)
        dynurban_interp(state.dynurban_state, grc, lun, bounds_proc.begg,
            bounds_proc.endg, current_year)
    end

    # (No dynHarvest_interp here — matching Fortran.)

    # Wrap up the weight changes. In Fortran this loops over clumps; here it runs
    # once over the proc bounds (single-clump dev path). Safe to always run.
    dynSubgrid_wrapup_weight_changes!(bounds_proc, grc, lun, col, pch;
                                      glc_behavior = glc_behavior)

    return state
end

# --------------------------------------------------------------------------
# dynSubgrid_driver!
# --------------------------------------------------------------------------
"""
    dynSubgrid_driver!(state, bounds_proc, grc, lun, col, pch;
                       year, clump_index=1,
                       use_crop=false, crop_inst=nothing, collapse_crops=false,
                       glc_behavior=nothing,
                       temp=nothing, ws=nothing,
                       cons_bgp_args=nothing, mask_nolakec=nothing,
                       mask_lakec=nothing)

Per-timestep top-level driver for transient land cover (Fortran
`dynSubgrid_driver`). Should be called once per time step from outside any clump
loop.

Sequence (mirroring the Fortran):
  1. (pre-change) snapshot heat/water baselines (`dyn_hwcontent_init!`) if the
     biogeophys state is supplied; snapshot prior weights and old weights into the
     patch/column state updaters.
  2. (I/O) call each enabled `*_interp` reader to update subgrid weights for `year`.
  3. (post-change, no I/O) wrap up weight changes (landunit reconcile + higher-order
     weights + active flags), snapshot new weights into the updaters, initialize
     newly-active columns (if `temp`/`ws` supplied), run the heat/water dynbal
     conservation final pass (if supplied).

`clump_index` selects the `any_changes` slot in the column state updater. FATES and
CNDV paths are deferred (no-op). The C/N (`use_cn`) biogeochem conservation is not
wired here (it requires the full CN vegetation facade); supply the biogeophys state
to exercise the water/energy conservation.

`cons_bgp` (a NamedTuple or struct with the heat/water state objects) enables the
biogeophys conservation; see [`dynSubgrid_run_conservation_init!`] /
[`dynSubgrid_run_conservation_final!`].
"""
function dynSubgrid_driver!(state::DynSubgridState, bounds_proc::BoundsType,
                            grc::GridcellData, lun::LandunitData,
                            col::ColumnData, pch::PatchData;
                            year::Int,
                            clump_index::Int = 1,
                            use_crop::Bool = false,
                            crop_inst = nothing,
                            collapse_crops::Bool = false,
                            glc_behavior = nothing,
                            temp::Union{TemperatureData,Nothing} = nothing,
                            ws::Union{WaterStateData,Nothing} = nothing,
                            cons_bgp = nothing,
                            mask_nolakec = nothing,
                            mask_lakec = nothing)

    @assert bounds_proc.level == BOUNDS_LEVEL_PROC "dynSubgrid_driver!: argument must be PROC-level bounds"

    ctl = state.ctl

    # ======================================================================
    # Initialization, prior to land cover change.
    # ======================================================================
    if cons_bgp !== nothing
        dyn_hwcontent_init!(state.dynbal, bounds_proc,
            mask_nolakec, mask_lakec, col, lun,
            cons_bgp.urbanparams, cons_bgp.soilstate,
            cons_bgp.waterstatebulk, cons_bgp.waterdiagnosticbulk,
            cons_bgp.temperature, cons_bgp.lakestate)
    end

    set_prior_weights!(state.prior_weights, bounds_proc, pch, col)
    set_old_weights!(state.patch_state_updater, bounds_proc, pch, col)
    set_old_weights!(state.column_state_updater, bounds_proc, grc, lun, col)

    # ======================================================================
    # Land cover change that requires I/O (outside threaded region in Fortran).
    # ======================================================================
    if get_do_transient_pfts(ctl)
        dynpft_interp!(state.dynpft_state, bounds_proc, pch, lun; year = year)
    end

    if get_do_transient_crops(ctl)
        dyncrop_interp!(state.dyncrop_state, bounds_proc, grc, lun, col, pch;
            year = year, use_crop = use_crop, crop_inst = crop_inst,
            collapse_crops = collapse_crops)
    end

    if get_do_harvest(ctl)
        # use_fates path deferred; CLM harvest interp always runs when enabled.
        dynHarvest_interp!(state.dynharvest_state, bounds_proc.begg,
            bounds_proc.endg, year)
    end

    if get_do_grossunrep(ctl)
        dynGrossUnrep_interp!(state.dyngrossunrep_state, year)
    end

    if get_do_transient_lakes(ctl)
        dynlake_interp(state.dynlake_state, grc, lun, bounds_proc.begg,
            bounds_proc.endg, year)
    end

    if get_do_transient_urban(ctl)
        dynurban_interp(state.dynurban_state, grc, lun, bounds_proc.begg,
            bounds_proc.endg, year)
    end

    # use_fates_luh path (dynFatesLandUseInterp) — DEFERRED (no FATES port).

    # ======================================================================
    # Land cover change that does not require I/O.
    # ======================================================================

    # CNDV bgc_vegetation_inst%UpdateSubgridWeights — DEFERRED.
    # use_fates dyn_ED — DEFERRED.
    # glc2lnd_inst%update_glc2lnd_fracs — no glc2lnd transient path on the base
    #   for this standalone driver; skipped.

    # Wrap up after land cover change.
    dynSubgrid_wrapup_weight_changes!(bounds_proc, grc, lun, col, pch;
                                      glc_behavior = glc_behavior)

    set_new_weights!(state.patch_state_updater, bounds_proc, pch)
    set_new_weights!(state.column_state_updater, bounds_proc, clump_index, col)

    # set_subgrid_diagnostic_fields — no history-diagnostics port; skipped.

    # Initialize newly-active columns (requires the temperature + water state).
    if temp !== nothing && ws !== nothing
        initialize_new_columns!(bounds_proc,
            state.prior_weights.cactive[bounds_proc.begc:bounds_proc.endc],
            temp, ws, grc, lun, col)
    end

    # Biogeophysics (water & energy) conservation, after land cover change.
    if cons_bgp !== nothing
        dyn_hwcontent_final!(state.dynbal, bounds_proc,
            mask_nolakec, mask_lakec, col, lun,
            cons_bgp.urbanparams, cons_bgp.soilstate,
            cons_bgp.waterstatebulk, cons_bgp.waterdiagnosticbulk,
            cons_bgp.temperature, cons_bgp.lakestate, ctl)
    end

    # use_cn DynamicAreaConservation (biogeochem C/N) — requires the full CN
    # vegetation facade; not wired into this standalone driver. The pieces
    # (dyn_cnbal_patch! / dyn_cnbal_col!) exist on the base and can be threaded in
    # when the CN facade is available.

    return nothing
end

# --------------------------------------------------------------------------
# dynSubgrid_wrapup_weight_changes!
# --------------------------------------------------------------------------
"""
    dynSubgrid_wrapup_weight_changes!(bounds, grc, lun, col, pch; glc_behavior=nothing)

Reconcile various variables after subgrid weights change (Fortran
`dynSubgrid_wrapup_weight_changes`):

  1. `update_landunit_weights!`     — reconcile landunit areas so they sum to 1.
  2. `compute_higher_order_weights!`— recompute col%wtgcell, patch%wtlunit/wtgcell.
  3. `set_active!`                   — re-derive active flags (the Fortran's
     `reweight_wrapup`, which re-makes filters; the available analogue on the base).

`glc_behavior` is accepted for signature compatibility with the Fortran but unused
(the base `set_active!` has no glc_behavior dependency).
"""
function dynSubgrid_wrapup_weight_changes!(bounds::BoundsType,
                                           grc::GridcellData, lun::LandunitData,
                                           col::ColumnData, pch::PatchData;
                                           glc_behavior = nothing)
    update_landunit_weights!(bounds, grc, lun)

    compute_higher_order_weights!(bounds, col, lun, pch)

    # Here: filters are re-made (Fortran reweight_wrapup). The available analogue
    # is set_active!, which re-derives the landunit/column/patch active flags.
    set_active!(bounds, lun, col, pch)

    return nothing
end
