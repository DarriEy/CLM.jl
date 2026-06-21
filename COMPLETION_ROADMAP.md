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
- **Tier C — Transient capability** (Tier B + dynamic/transient land use).
- **Tier D — Project goals met** (reverse-AD complete, multi-GPU, AD-smoothing) — the CLAUDE.md mandate.
- **Tier E — Full CTSM** (Tier C/D + FATES + alternative method options + MPI). The literal 100%.

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

**B1. History output** — port `histFileMod.F90` (~6,100 lines): field registration, time-averaging, NetCDF write. Largest single module. *DoD:* a run writes a CTSM-comparable `h0` history file.
**B2. Restart I/O** — `restFileMod.F90` (~1,034) + `ncdio_*` + `subgridRestMod`. *DoD:* write a restart, read it back, bit-for-bit continue.
**B3. Init/domain plumbing** — `decompInitMod`, `domainMod`, `subgridMod`, `init_interp` (11 files) as needed for standalone init. *DoD:* cold-start + finidat without the harness.

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
