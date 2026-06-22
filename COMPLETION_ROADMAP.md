# CLM.jl — Roadmap to 100%

Companion to `PRD_CLM_JULIA_PORT.md`. Derived from the module-level audit of the
CTSM Fortran tree (`installs/clm/src`, excl. FATES) against `src/`:
**261 core modules, 172 ported (~66% by file), 89 not** — but the unported set is
mostly roads-not-taken (alternative methods), infrastructure/IO, the transient
land-use subsystem, and FATES. The *single-point CLM5 process model is ~complete
and Fortran-validated*; this roadmap covers everything between that and "100%".

## What "100%" means — pick the target tier
- **Tier A — Complete single-point CLM5 physics** (process fidelity for one column/gridcell). Closest to "done" for science use.
- **Tier B — Usable standalone model** (Tier A + NetCDF history output + restart I/O).
- **Tier C — Transient capability** (Tier B + dynamic/transient land use). *Dependency-ordered batch plan (whole `dyn_subgrid` tree, ~2,900 Julia lines, was fully unported — only base geometry `init_subgrid`/`subgrid_weights`/`subgrid_ave` existed):*
  - **C-Batch 1 — Foundation ✅ done (PR #53):** `dyn_subgrid_control.jl` (DynSubgridControl flags + namelist consistency checks, 48 tests) + `dyn_file_io.jl` (DynTimeInfo year↔file-index mapping + DynFile NetCDF wrapper + DynVarTimeUninterp step-wise reader, 68 tests).
  - **C-Batch 2 — Data readers + weight updates ✅ done (PR #57):** `dyn_pft_crop_file.jl` (dynpft/dyncrop PCT_NAT_PFT/PCT_CROP/PCT_CFT/FERTNITRO readers, 21 tests), `dyn_lake_urban_file.jl` (dynlake/dynurban PCT_LAKE + 2D PCT_URBAN density-class readers, 30 tests), `dyn_landunit_area.jl` (dynLandunitArea weight normalization w/ DECREASE_ORDER priority + ISTICE preserved, dynPriorWeights snapshot, canonical `set_landunit_weight!`, 38 tests), `dyn_init_columns.jl` (dynColumnTemplate exemplar selection + dynInitColumns new-column state copy w/ snow-layer offset, 25 tests).
  - **C-Batch 3 — State conservation ✅ done (PR #59):** `dyn_patch_state_updater.jl` (dynPatchStateUpdater — conservative patch-state redistribution on weight change + flux-partition-by-PFT variant + grew/initiating/was-zero queries, 44 tests), `dyn_column_state_updater.jl` (dynColumnStateUpdater — all 5 fill variants: no-special/fill-natveg/fill-fixed-values/fill-special-fixed/optional-fractions, reuses `template_col_from_natveg_array!`, 33 tests), `dyn_cons_biogeophys.jl` (dynConsBiogeophys — water/heat baselines + before/after content integrals → qflx/eflx dynbal fluxes, reuses `total_water_heat.jl` integrators + `c2g`, 26 tests).
  - **C-Batch 4a — CN dynamics ✅ done (PR #60):** `dyn_harvest.jl` (dynHarvest — 5 HARVEST_* readers + CNHarvest pool-by-pool C/N mortality → product/litter + patch→col, 38 tests), `dyn_gross_unrep.jl` (dynGrossUnrep — UNREPRESENTED_PFT_LULCC reader + CNGrossUnrep disturbance C/N → litter/atm, reuses existing `gru_*` flux fields, 45 tests), `dyn_cons_biogeochem.jl` (dynConsBiogeochem — `dyn_cnbal_patch!`/`dyn_cnbal_col!` C/N conservation across cover change, reuses patch/column state updaters + `compute_seed_amounts!`, 19 tests).
  - **C-Batch 4b — Orchestration ✅ done (PR #61):** `src/driver/dyn_subgrid_driver.jl` (dynSubgridDriver — `dynSubgrid_init!`/`dynSubgrid_driver!`/`dynSubgrid_wrapup_weight_changes!`, bundling all reader/updater/conservation state into a standalone `DynSubgridState`; init → per-year interp → prior-weights snapshot → landunit-area reconcile → patch/column state updaters → new-column init → optional biogeophys conservation; FATES/CNDV paths gated off, 23 tests).
  - **→ TIER C COMPLETE:** the whole `dyn_subgrid` transient-land-use subsystem is ported (control/time/file readers → pft/crop/lake/urban data → weight normalization + prior snapshot + new-column init → patch/column state conservation + biogeophys/biogeochem dynbal → harvest/gross-unrep CN → top-level driver). ~2,900 Julia lines across PRs #53/#57/#59/#60/#61; suite 16,650 → 17,108. The driver is independently callable; **not yet wired into the real `clm_drv!` timestep loop** (a single follow-up: thread `DynSubgridState` through the annual step + supply the full CN-veg facade for the `dyn_cnbal_*` path).
  - **Deferred (later tiers):** FATES land-use (dynFATESLandUseChange/dynED) + CNDV (dynCNDV) + full MPI broadcast.
- **Tier D — Project goals met** (reverse-AD complete, multi-GPU, AD-smoothing) — the CLAUDE.md mandate.
- **Tier E — Full CTSM** (Tier C/D + FATES + alternative method options + MPI). The literal 100%.
- **Tier F — FATES** (Functionally Assembled Terrestrial Ecosystem Simulator) — ~77k lines of Fortran / 68 files under `installs/clm/src/fates/` (main/biogeochem/biogeophys/parteh/radiation/fire). Greenfield → `src/fates/`. Dependency-ordered wave plan (each wave merges before the next; modules within a wave are mutually independent → parallel agents):
  - **F-Batch 0 — Foundation ✅ done (PR #62):** `src/fates/` — FatesConstantsMod, FatesGlobals, FatesIntegratorsMod (Euler+RKF45), FatesUtilsMod, FatesRunningMeanMod, FatesIODimensionsMod, FatesIOVariableKindMod, FatesParametersInterface, FatesSynchronizedParamsMod (9-module bedrock closure; `fates_r8`→Float64; logging→@warn/error; I/O stubbed; 98 tests).
  - **F-Batch 1 — Core leaf modules ✅ done (PR #63):** EDParamsMod (ed_params Ref-holder + registration, 72t), FatesHydroWTFMod (van Genuchten/Campbell/smooth-CCH/TFS water-retention+conductance fns w/ analytic derivs, 158t — ⚠️ banked: upstream Fortran VG `dftcdpsi_from_psi` cross-term sign disagrees w/ FD ~2-10%, preserved verbatim), FatesInterfaceTypesMod (bc_in/bc_out/bc_pconst host-boundary types, 99t), FatesLitterMod+FatesRadiationMemMod (litter pools + rad band indices, 84t), PRTParametersMod+SFFireWeatherMod (PARTEH param arrays + fire-weather base, 172t), TwoStreamMLPEMod (multi-layer 2-stream RT solver, LU via LinearAlgebra, energy-conserving on unit-incident contract, 38t). +623 tests; suite 17,206→17,829.
  - **F-Batch 2 — ✅ done (PR #64):** PRTGenericMod (generic PARTEH element/organ framework + abstract hypothesis interface, 77t), FatesHydraulicsMemMod (ed_site/cohort_hydr state, `Vector{WRFType/WKFType}` for the soil WTFs, 104t), FatesSizeAgeTypeIndicesMod (13 size/age/class flat-index fns, 404t), FatesFuelClassesMod+SFNesterovMod (fuel-class enum + Nesterov fire-weather index, 32t). +617 tests; suite 17,829→18,446. **F-Batch 3 — ✅ done (PR #65):** EDPftvarcon (FATES PFT param table, Ref-held, 41t), PRTLossFluxesMod (PARTEH turnover/retransloc/repro/burn/damage loss fluxes, 47t), SFParamsMod (SPITFIRE params + CWD-sum consistency, 87t), FatesFuelMod (fuel loading/moisture/bulk-density, 54t). +229 tests; suite 18,446→18,675. **F-Batch 4 (2):** FatesDispersalMod, FatesParameterDerivedMod. **F-Batch 5 (1):** DamageMainMod.
  - **F-Batch 6:** FatesAllometryMod (3.5k, allometry engine). **F-Batch 7 (2):** PRTAllometricCNPMod, PRTAllometricCarbonMod. **F-Batch 8:** FatesCohortMod. **F-Batch 9:** FatesPatchMod. **F-Batch 10:** EDTypesMod → **TYPE SYSTEM COMPLETE** checkpoint.
  - **F-Batch 11 (8, ~11.5k):** ChecksBalances, EDAccumulateFluxes, FatesLandUseChange, FatesPlantHydraulicsMod (6.4k), FatesSoilBGCFlux, FatesTwoStreamUtils, PRTParamsFATES, SFMainMod. **F-Batch 12 (3):** EDBtran, EDCohortDynamics, EDLoggingMortality. **F-Batch 13 (3):** EDMortalityFunctions, EDPhysiologyMod (3.4k), FatesBstress → physiology checkpoint.
  - **F-Batch 14:** EDPatchDynamicsMod (3.9k). **F-Batch 15:** EDCanopyStructureMod → structure checkpoint. **F-Batch 16 (4):** EDMainMod (driver), FatesInventoryInit, FatesNormanRad, FatesPlantRespPhotosynth. **F-Batch 17 (2):** EDInitMod, FatesRadiationDrive → init checkpoint.
  - **F-Batch 18 — host coupling/I/O (5, ~16k, PORT LAST/deferrable):** FatesHistoryVariableType, FatesHistoryInterfaceMod (9.3k), FatesRestartVariableType, FatesRestartInterfaceMod, FatesInterfaceMod (the CLM↔FATES coupling point).

Recommended definition of "done" for this project: **Tier A + Tier D** (a complete, differentiable, GPU single-point CLM5). Tiers B/C/E are large and only needed for specific use cases.

---

## Phase A — Complete single-point CLM5 physics  *(highest value/effort)*
Close the genuine process gaps + verify depth. Target: Tier A.

> **Phase A status (this round):** A1 ✅, A6 ✅ (was a false gap — already ported, added tests), A7 ✅, A8 ✅ (snow already ported; dust thresholds ported), A2/A9 ✅. Remaining real single-point gaps shrank to **three small items** (see A10). See `AUDIT_SUBROUTINE.md`.

**A1. Prognostic urban building temperature** — ✅ DONE. Ported `UrbBuildTempOleson2015Mod` → `urb_build_temp_oleson2015.jl` (5×5 energy-balance solve), wired into the `soil_temperature!` prog branch. 20 tests (energy-balance closure + ForwardDiff check); urban smoke finite in SIMPLE + PROG.

**A10. Remaining small single-point gaps — ✅ ALL DONE.**
- `species_from_string` ✅ (PR #40) — C12/C13/C14/N parser.
- `UpdateAccVars_CropGDDs` ✅ (PR #40) — crop-GDD accumulation (accumulMod was already ported); now also **wired into `temperature_update_acc_vars!` + GDD020/820/1020 runmeans** (this round).
- **MEGAN/VOC** ✅ (this round, AD-safe) — `megan_factors_get`/`gen_hashkey` (byte-exact vs Fortran) + `voc_emission!` wired behind a `use_voc` flag; the MEGAN descriptors live on a `MEGANConfig` side-struct on `CLMDriverConfig` (NOT the dual-copied `CLMInstances` — that placement was the AD regression that deferred it the first time). Still needs a MEGAN namelist parser to activate by default (no-op until then).

**→ With A10 closed, the single-point CLM5 physics port (Tier A) is effectively complete.** Remaining work is the big tiers: B (I/O), C (transient dyn_subgrid), D (reverse-AD finish + CUDA/AMD + AD-smoothing), E (alt-method options), F (FATES).

**A2. Subroutine-level audit** — ✅ DONE (`AUDIT_SUBROUTINE.md`): of 138 ported physics modules, **102 full / 23 partial / 13 stub** (most stubs are pure type-defs). The real partials it surfaced become the tasks below (A6–A9):

**A6. Soil lateral flow + irrigation routing** (`SoilHydrologyMod` 9/16 — the biggest gap) — port `PerchedLateralFlow`, `SubsurfaceLateralFlow`, `RouteInfiltrationExcess`, `SetFloodc`, `WithdrawGroundwaterIrrigation`. Core water-table/drainage is ported; lateral redistribution + irrigation withdrawal aren't. *DoD:* lateral flow + groundwater irrigation active, water balance still closes.

**A7. Glacier surface mass balance** (`GlacierSurfaceMassBalanceMod` 1/3) — port `HandleIceMelt`, `ComputeSurfaceMassBalance`. *DoD:* glacier SMB on the istice path, gated-tested.

**A8. Snow capping + dust thresholds** — `SnowHydrologyMod::SnowCapping` (excess-snow→ice/runoff) and `SoilStateInitTimeConstMod` dust soil-moisture thresholds (`ThresholdSoilMoist{Zender2003,Kok2014}`). *DoD:* capping conserves mass; dust emission threshold active.

**A9. Verify suspected false-negatives + h2osfc updates** — confirm `SurfaceWaterMod` h2osfc-fraction updates aren't silently dropped (may be folded into the hydrology driver); confirm the audit's likely false-negatives (`SaturatedExcessRunoff` kernel, `UrbanInput`=ported this session, MEGAN/species/GDD inlined utilities). *DoD:* each either confirmed present or filed as a task.

*(Minor, low-priority: parameter-reader routines — `readParams` in CNAllocation/CNFUN/NLeaching/InfiltrationExcess, `SnowOptics_init` — the port loads these via its own `read_params.jl`; alternative methods — `soilwater_zengdecker2009`, VIC init — are Phase E.)*

**A3. Close documented physics residuals** —
  - the ~64 mm/yr snowmelt water-balance leak (dead in-model `endwb_col` check; `scripts/longhorizon_conservation.jl`);
  - LUNA via injection → finish the real Rubisco-N optimum (the `luna-wiring-status` residual).
  *DoD:* water balance closes < a few mm/yr; LUNA vcmax matches Fortran without injection.

**A4. Small physics sub-features** — `TillageMod` (crop management), C-isotope atm forcing reader (`CNCIsoAtmTimeSeriesReadMod`; c-iso flux physics already ported), optionally `HumanIndexMod` (urban heat-stress diagnostic). *DoD:* each gated-tested.

**A5. Water isotope tracers** — `WaterTracer*` / `WaterInfoTracerType` plumbing (bulk water is ported; isotopic-tracer variants are not). *DoD:* a tracer (e.g. H2_18O) advects through the water cycle. *(Skip unless isotope science is needed.)*

---

## Phase B — Usable standalone model (I/O)  *(large, mechanical)*
Make it run + write output + restart without the parity harness's injection.

> **Tier A activations COMPLETE:** crop-GDD threaded into `clm_drv_core!` (AccumManager on `CLMDriverConfig`, AD-safe) + MEGAN namelist parser (`megan_exp_parse`/`megan_config_from_nl`) — both use_crop/use_voc paths now functional.

**B0. Standalone run harness** — ✅ done (`src/driver/run_clm.jl`): `run_clm!(...)` orchestrates init(±`finidat` restart) → driver loop → periodic h0 write → end restart, no parity harness. 15-test gated round-trip (run→write h0+restart→finidat-continue exact). Pure orchestration, AD-safe.
**B1. History output** — ✅ done (`src/infrastructure/history_io.jl`): add/accumulate/write (A/I/X/M) + multi-tape `HistoryTapeSet` (h0/h1/h2, per-tape nhtfrq) + `fincl`/`fexcl` w/ `name:flag` + master field list, **+ per-subgrid-level dims (PATCH/COLUMN/LANDUNIT/GRIDCELL — mixed-level tapes no longer DimensionMismatch) + p2c/c2l/l2g gridcell remap (`hist_dov2xy`) + multi-level 2D fields (TSOI etc.)** (112 tests). DEFERRED (minor): time_bounds, mfilt rollover, ndens.
**B2. Restart I/O** — ✅ done (`src/infrastructure/restart_io.jl`): `write_restart!`/`read_restart!` round-trip of prognostic biogeophys + (use_cn) C/N pools **+ CN flux/accumulators, crop/phenology counters, C13/C14 isotopes** (flag-gated), CLM-faithful names (93 tests). DEFERRED: full header metadata.
**B3. Init/domain plumbing** — ✅ `init_interp` done (`src/infrastructure/init_interp.jl`): cross-grid finidat — subgrid type-match + vertical level interpolation; reduces to `read_restart!` on exact-grid match (650 tests). *(Single-process mask-based init via `clm_initialize!` already works; the Fortran MPI-decomp modules `decompInitMod`/`domainMod`/`subgridMod` are N/A for this port.)*

**→ Tier B (usable standalone model) effectively COMPLETE:** standalone run harness, history (multi-tape + subgrid + 2D), restart (full prognostic state), and finidat/init_interp all work. Remaining = the minor history metadata bits + multi-gridcell spatial regridding (not needed for the single/few-column eval domains).

---

## Phase C — Transient capability (dynamic land use)
Port the `dyn_subgrid/` subsystem (19 modules, 0 ported): transient PFT/crop/urban/lake weights, harvest, land-use change, conservation on weight shifts. *DoD:* a transient (e.g. 1850–2015 land-use) run conserves C/N/water across weight changes. *(Needed only for transient/historical experiments.)*

---

## Phase D — Project goals (GPU + AD)  *(the CLAUDE.md mandate)*
**D1. Reverse-AD (Enzyme)** — finish the BGC phase chain (in progress: `wire-cn-state-summaries`, `driver_reverse.jl`). Extend the FD-validated reverse through the full `clm_drv!`. *DoD:* end-to-end reverse-AD gradient of a driver output w.r.t. parameters, FD-validated.
**D2. Multi-GPU backends** — CUDA + AMD (Metal done). Mostly backend-swap + Float64-literal/scalar-arg hardening already done for Metal. *DoD:* the biogeophys+BGC driver runs on CUDA at parity.
**D3. AD smoothing (PRD Phase 3)** — smooth remaining discontinuities (snow merge/split, phase change) for clean gradients. *DoD:* gradients finite + smooth across the documented discontinuities.

---

## Phase E — Alternative method options *(as-needed, not gaps)*
CTSM ships multiple implementations; the port chose one each. Port others only if a use case needs them: matrix CN solver (`CNVegMatrixMod` 3851 + `CNSoilMatrixMod` + `SparseMatrixMultiplyMod`), other fire models (`CNFireLi2016/2021/2024`), `DustEmisLeung2023`, `NutrientCompetitionFlexibleCNMod`, ozone base/factory. *DoD:* per-option parity + a config switch.

---

## Phase F — FATES  *(a project unto itself, ~52K lines)*
The demographic vegetation model (`src/fates/`). Explicitly deferred in the PRD as Phase 2. Separate PRD + sustained effort. *DoD:* FATES vegetation dynamics coupled through the host interface.

---

## Continuous (across all phases)
- **Eval-domain validation** (in progress): finish the Symfluence builds, run `run_clm_streamflow.jl all`, score the 7 wired gauges (Bow/Iceland/Baltimore/Massa/Tagus/Krycklan/Abisko).
- **Fortran parity expansion**: extend per-step parity beyond Bow to the multisite/eval domains.
- Keep the full test suite green (currently 15,700+).

---

## Suggested sequencing
1. **Phase A** first — it's the cheapest path to a defensible "100% of single-point CLM5," and A2 (subroutine audit) tells us how much of the 172 is *really* done.
2. **Phase D** in parallel — it's the project's stated differentiator (GPU + AD) and is already underway.
3. **Phase B** when standalone output/restart is needed (e.g. for longer production runs).
4. **Phase C / E / F** on demand, by science need.

**Start now:** A2 (subroutine-level audit) — it converts "~66% by file" into a precise, prioritized task list and may reveal Phase A is smaller (or larger) than it looks.
