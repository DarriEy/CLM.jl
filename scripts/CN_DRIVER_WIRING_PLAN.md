# CN veg-flux driver wiring plan (Track 2, #8)

Goal: wire the vegetation-CN flux routines into `cn_driver_no_leaching!`
(src/biogeochem/cn_driver.jl) so a `use_cn` `clm_drv!` step computes the C/N
fluxes that the already-wired state-update cascade consumes. Validate per-step vs
the spun-up Fortran BGC reference (`scripts/fortran_parity_cn_step.jl`,
dumps in `clm_bgc_spinup/bgc_ref_dumps_2202/`, nstep 1753153).

## Current state (2026-06-15)
- WIRED in cn_driver_no_leaching!: flux zeroing, decomposition chain (cond. on
  inst infra), c_state_update0/1/2/2h/2g/3, n_state_update1/2/2h/2g,
  cn_precision_control, litter_vert_transp.
- NOT WIRED (comments): cn_mresp!, calc_gpp_mr_availc!, calc_allometry!,
  calc_crop_allocation_fractions!, calc_plant_nitrogen_demand!, plant N
  competition→fpg_col, calc_plant_cn_alloc!, cn_phenology! (phase 1+2),
  cn_gresp!, cn_gap_mortality!+cn_gap_patch_to_column!.
- The decomp infra IS on inst + built in clm_initialize! when use_cn=true
  (inst.decomp_cascade/decomp_bgc_state/decomp_bgc_params/cn_shared_params/
  competition_state/params/litter_params/decomp_params) and passed conditionally
  at the clm_driver.jl:1352 call site (gated on _decomp_initialized).
- The veg-CN PARAM structs + PFTCON VARIANTS are NOT on inst and NOT read in
  production — only built ad-hoc in tests. THIS IS THE FIRST BUILD TASK.

## Fortran canonical order (CNDriverNoLeaching.F90) — the wiring order
1. cn_mresp!                         (maint resp)            [params: MaintRespParams, PftConMaintResp{woody}; reads photosyns lmrsun/sha, canopystate laisun/sha, soilstate crootfr, temperature, cnveg_ns; writes cnveg_cf *_mr]
2. cn_phenology! phase=1 (if use_fun) — our use_fun=false → at step 10
3. calc_gpp_mr_availc!               (gpp/mr/availc)         [AllocationParams, PftConAllocation, crop, photosyns psnsun/sha, canopystate laisun/sha; writes cnveg_cf availc/gpp/psn*_to_cpool]
4. calc_crop_allocation_fractions!   (crops only; Bow mask_pcropp empty → no-op)
5. calc_allometry!                   (c/n allometry)         [PftConAllocation; writes cnveg_state c_allometry/n_allometry]
6. calc_plant_nitrogen_demand! ×2    (non-crop, crop)        [PftConNutrientCompetition, crop; writes cnveg_nf plant_ndemand + retransn fluxes]
7. plant N competition → fpg_col     (SoilBiogeochemCompetition plant side) — PRODUCES fpg_col (frac of plant N demand met). The decomposer-side soil_bgc_competition! is already wired; need the PLANT fpg_col output. CHECK: does soil_bgc_competition! already compute fpg/fpg_col, or is a separate routine needed?
8. calc_plant_cn_alloc!              (C/N allocation)        [PftConNutrientCompetition, fpg_col; writes cnveg_cf cpool_to_*c (consumed by c_state_update1!), cnveg_nf sminn_to_npool]
9. soil_biogeochem_decomp!           (already wired)
10. cn_phenology! phase=1 (if !use_fun) — OUR case             [PhenologyState, PhenologyParams, PftConPhenology, temperature, water_diag, canopy_state, soil_state, gridcell, leaf_prof/froot_prof]
11. cn_phenology! phase=2 (all)
12. cn_gresp!                        (growth resp)           [PftConGrowthResp; reads cpool_to_* from step 8; writes cnveg_cf *_gr]
13. c_state_update0/1, n_state_update1 (WIRED)
14. cn_gap_mortality! + cn_gap_patch_to_column! (mortality)  [GapMortalityParams, PftConGapMort, DgvsGapMortData, leaf/froot/croot/stem_prof; writes m_*_to_litter]

## Build tasks (in order)
### A. Param + pftcon infrastructure (the prerequisite)
Add to CNVegetationData (vegetation_facade.jl) + populate in a new
`cn_vegetation_read_params!(veg, paramfile, pftcon)` called from clm_initialize!
(when use_cn): MaintRespParams, AllocationParams, GrowthRespParams,
GapMortalityParams, PhenologyParams (via their *_read_params! readers — already
exist), and the pftcon variants PftConMaintResp/Allocation/NutrientCompetition/
GrowthResp/GapMort/Phenology (each is a SUBSET of the main pftcon → write a
builder that copies the relevant fields from inst's main pftcon; the main pftcon
already loads woody/leafcn/slatop/froot_leaf/grperc/etc.). Also allocate
PhenologyState + DgvsGapMortData (or DGVSData) on the veg struct + cold-init.
Verify field availability: each pftcon-variant field must exist on the main
pftcon (most do; flag any missing).

### B. Thread inputs through the facade → driver
Add to cn_vegetation_ecosystem_pre_drainage! + cn_driver_no_leaching! signatures
and the clm_driver.jl:1352 call site: photosyns (inst.photosyns), canopystate
(inst.canopystate), soilstate (inst.soilstate), temperature (inst.temperature),
water_diag (inst.water.waterdiagnosticbulk_inst), gridcell (inst.gridcell), and
the leaf_prof/froot_prof/croot_prof/stem_prof matrices (from where? — check
phenology/the cnveg_state for the vertical profiles; in Fortran they come from
the canopy/decomp vertical-profile calc, already wired as soil_bgc_vertical_profile!).

### C. Wire routines in cn_driver_no_leaching! in the order above
Each guarded on `num_bgc_vegp > 0` and the relevant params !== nothing. Reuse the
NaN×0 guards already added (allocation.jl:185, maint_resp.jl:116).

### D. fpg_col (the hard part)
Determine how fpg_col (plant N uptake fraction) is produced. In Fortran,
SoilBiogeochemCompetition computes fpg (plant) and fpi (decomposer). Check if the
Julia soil_bgc_competition! already outputs fpg_col on soilbgc_state/competition_state;
if not, port that piece. calc_plant_cn_alloc! needs it.

## Validation per increment
After each routine wired: clean-rebuild cache
(`rm -rf ~/.julia/compiled/v1.12/CLM`; NOTE Google-Drive FUSE mtime is unreliable
→ ALWAYS nuke the cache to force recompile), then
`julia +1.12 --project=. scripts/cn_nan_probe.jl` → the corresponding flux/state
field should go from NaN→finite. Final: `scripts/fortran_parity_cn_step.jl` →
per-pool max|rel| vs the after_hydrologydrainage dump.

## PROGRESS 2026-06-15
- DONE: threaded veg-flux inputs (patch/pftcon_main/crop/photosyns/canopystate/
  soilstate/temperature) through clm_driver.jl:1352 call site → facade
  (cn_vegetation_ecosystem_pre_drainage!) → cn_driver_no_leaching! (all optional
  =nothing so SP/decomp-only callers unaffected).
- WIRED + VERIFIED (scripts/cn_nan_probe.jl, clean-cache rebuild):
  1. cn_mresp! → leaf_mr nan 4/4→1/4 (col-1 finite [0,0,1.6e-13]; the 1 NaN is
     col-2 inactive patch4). Builds MaintRespParams() (reader defaults br=2.525e-6)
     + PftConMaintResp{Float64}(woody=pftcon_main.woody).
  2. calc_gpp_mr_availc! → availc nan 4/4→1/4 (col-1 finite [0,0,0]). Builds
     PftConAllocation{Float64} from main pftcon (11 fields all present) + AllocationParams().
- fpg_col UNBLOCK: it's `soilbgc_state.fpg_col`, written by soil_bgc_competition!
  (ALREADY wired in the decomp block) as sminn_to_plant/plant_ndemand. BUT
  **ordering**: Fortran runs calc_plant_nitrogen_demand! (→plant_ndemand) BEFORE
  SoilBiogeochemCompetition (→fpg), then calc_plant_cn_alloc! AFTER. The current
  Julia driver runs soil_bgc_competition! EARLY (inside the decomp block, lines
  ~256) before any plant_ndemand. **To wire allocation correctly the driver must
  be REORDERED**: mresp→gpp→allometry→ndemand→[move soil_bgc_competition! here]→
  calc_plant_cn_alloc!→decomp(rest). This is the next structural step.
- REMAINING (next session): thread cnveg_state (veg.cnveg_state_inst) to the
  driver; wire calc_allometry! (have _pfta) → c_allometry/n_allometry; build
  PftConNutrientCompetition (subset of main pftcon) + wire calc_plant_nitrogen_demand!
  ×2 → plant_ndemand; REORDER soil_bgc_competition! after ndemand; wire
  calc_plant_cn_alloc! (fpg_col) → cpool_to_*c (then c_state_update1! makes leafc
  finite — THE milestone); then cn_phenology! ph1/ph2 (PhenologyState/Params/profiles),
  cn_gresp! (PftConGrowthResp), cn_gap_mortality!+patch_to_column (DgvsGapMortData/
  GapMortalityParams/profiles). Then parity-tune vs the dump.

## MILESTONE 2026-06-15b — veg-CN chain wired, use_cn step at NEAR-PARITY (~1e-4)
Wired into cn_driver_no_leaching! (in Fortran order), each builds its param/pftcon-
subset in-driver from the main pftcon (all fields present):
- cn_mresp! (after flux zeroing)
- calc_gpp_mr_availc!  → availc
- calc_allometry!      → c_allometry/n_allometry
- calc_plant_nutrient_demand! (non-crop) → plant_ndemand_patch + a patch→column
  "unity" aggregation (wtcol-weighted manual loop) into soilbgc_state.plant_ndemand_col
- calc_plant_nutrient_competition! (= calc_plant_cn_alloc!) INSERTED in the decomp
  block AFTER soil_bgc_competition! (which now sees real plant_ndemand → real fpg_col)
  and BEFORE soil_biogeochem_decomp! → produces cpool_to_*c.
Also threaded cnveg_state (veg.cnveg_state_inst) to the driver.
**KEY FIX — flux zeroing was a STUB**: _zero_cnveg_cflux!/_nflux! only zeroed 2/1
fields → all other per-step fluxes stayed NaN-init → c_state_update1! made leafc NaN.
Replaced with `_zero_cnveg_flux_arrays!` = full per-step reset (fill! every Real array
field EXCEPT ann* accumulators), matching Fortran CNVeg*FluxType%SetValues. Device-safe.
**RESULT**: use_cn step runs; ALL col-1 CN pools FINITE (leafc=[0,7.49,0.205] etc., =injected
since winter-night step barely changes pools); soil decomp pools finite (the 140/280 NaN
= col-2 inactive). scripts/fortran_parity_cn_step.jl: per-pool max|rel| vs the
after_hydrologydrainage dump = **1e-4 (xsmrpool) … 1e-7, several exact**. From NaN→near-parity.
No regressions (test_allocation 54/54, test_growth_resp 44/44).
REMAINING (parity tail): wire cn_phenology! ph1/ph2 + cn_gresp! (winter-night fluxes ~0 so
small now, needed for growing-season); use param-FILE values instead of reader defaults
(br etc.) to tighten ~1e-4 diffs; fix 1 deep-soil-N NaN level in col-1 sminn_vr (level >
nlevdecomp not updated); cn_gap_mortality!; then full parity-tune + a growing-season window.

## Scope honesty
This is the scoping report's multi-week tail. Each of A–D is substantial; D
(plant N competition fpg_col) + getting allocation to parity is the deepest.
Increment + verify finiteness at each step; don't expect one-pass parity.
