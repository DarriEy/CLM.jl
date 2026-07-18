# Fortran-dependent validation backlog

**What this is.** The list of CLM.jl subsystems whose Fortran parity is *not yet
established* because it requires **generating new Fortran (CTSM) ground truth** — a new
site, a non-soil landunit, or a config-flag rerun — that does not exist on disk today.

**When to do it.** The CTSM source is re-clonable and buildable on any box
(`git clone -b ctsm5.3.012 https://github.com/ESCOMP/CTSM.git`, commit
`ab466d6f9789ca3df2c72bda46cf7afed2d04102`), so this is **anytime** work — it is *not*
gated by any particular machine. (Contrast: Metal-GPU validation, which is gated by the
Apple-Silicon dev box and must happen before it is returned.)

**What is deliberately NOT here** (doable anytime, needs no new Fortran run):
- Fire `NFIRE`/`FAREA_BURNED` — blocked only on porting accumulator-restart I/O; scores
  against the existing `bgc_ref_firech4` dump. (In progress.)
- Diffing the **existing-but-unread** mid-CN-driver dumps `after_competition` /
  `after_ecosysdyn_predrain` (ground truth already on disk, no harness reads it).
- Albedo bands (`albd/albi/fabd/fabi`), LUNA-recomputed `vcmx25_z/jmx25_z` (harness
  currently *injects* it), aerosol/dust/ozone — all present in existing dumps, just not
  diffed. Pure harness work.

**Run-verified harness status (2026-07-17):** see
[`docs/PARITY_COVERAGE_2026-07.md`](PARITY_COVERAGE_2026-07.md) for the refreshed
scorecard — every `scripts/fortran_parity_*.jl` added since the 2026-06-17 audit was
executed against the surviving local dumps. Headlines: active-layer / soilresis /
albedo (no-injection) are now **VALIDATED**; aerosol / daylength are green-by-skip
(need a snowy / `dayl` dump); **LUNA's EOD vcmax update DIVERGES ~5.3–5.6% high** (a
documented Rubisco-N optimum residual, jmax matches exactly); lake was validated
historically but its dump domain is gone; isotopes still need a c13/c14 spinup.

Cross-references: `docs/HARNESS_COVERAGE_AUDIT.md` (per-module coverage, tables A/B/C),
`docs/CH4_FIRE_PARITY.md` (§7 "not validated"), `docs/N_CYCLE_PARITY.md`,
`docs/PARITY_STATUS.md`. Ground-truth recipe: the `restFile_write_dump` /
`bgcdumpMod.F90` SourceMods instrumentation under `scripts/validation/`, run an
instrumented `cesm.exe` **plainly** (never under lldb — macOS xzone `brk` guards) with
`BGCDUMP_NSTEP_LO/HI` + `PDUMP_NSTEP_LO/HI` set to the diff window.

---

## A. Needs a new SITE (surfdata + forcing)

| # | Item | Julia module(s) | What a Fortran run must exercise | Candidate site (migrated) |
|---|------|-----------------|----------------------------------|---------------------------|
| A1 | **Methane as a SOURCE** (transport / aerenchyma / ebullition / `finundated>0`) — the regime CH4 actually matters for | `methane.jl` | `use_lch4=.true.` with `finundated>0`; diff `CH4_EBUL_*`, `CH4_AERE_*`, `CH4_SURF_*`, `CONC_CH4/O2_*`, `TOTCOLCH4`. **LARGELY UNBLOCKED 2026-07-18 (round 3, docs/CH4_FIRE_PARITY.md §11).** The unblock was NOT a wetter site but a wetter finundation METHOD: #238's `finundation_mtd=TWS_inversion` (Prigent satellite regression, decoupled from the column's own pond) gives `finundated≈0.045>0` at **Bow's own gridcell** — and Bow has a CONVERGED BGC restart (no cold-CN confound). Reference `clm_bgc_spinup/bow_ref_ch4tws`, harness `scripts/fortran_parity_ch4_tws.jl`. The sat source regime now agrees to the single-step-oracle floor (~1e-3–1e-4): `CH4_AERE_SAT` 0.61→**4.7e-3**, `CH4_SURF_AERE_SAT`→**1.1e-3**, `CONC_O2_SAT`→**2.5e-4**, whole UNSAT column EXACT. Two causes found: (1) a harness datm misalignment (the #233 lesson repeats — align to `clmforc.2002` −1h); (2) a REAL mis-port FIXED — grass/crop patches were denied the full aerenchyma porosity (`poros_tiller` applied `nongrassporosratio` unconditionally). REMAINING: `CH4_PROD_SAT` 7.9% (warm-CN somhr-not-dumped confound, not a kernel bug), and **ebullition is still physically zero at dry Bow** (below the bubble threshold — needs a wet, high-production column WITH a converged BGC restart, which no migrated site provides). MerBleue does NOT inundate under TWS_inversion either (its cell needs `TWS>2242`). **EBULLITION precise blocker CONFIRMED 2026-07-18 (§13):** the `_meth_ebul_kernel!` port is line-by-line-verified (`bubble_f=0.57`, `vgc_max=vgc_min=0.15`, threshold `vgc>0.0855`); the Bow-TWS reference itself gives `CH4_EBUL_TOTAL_SAT≡0` (peak `CONC_CH4_SAT≈0.029 mol/m³` ≪ threshold) → nothing to diff. Missing asset = a wetland/high-WT site WITH a CN/CH4 spinup restart (Bow-TWS has inundation+restart but low production; MerBleue has production potential but no inundation and no BGC restart). | `domain_Bow_at_Banff_lumped` (**inundates ~4.5% under TWS_inversion**); MerBleue does not |
| A2 | **Fire crop/peat/tropical/land-use branches** | `fire_li2016.jl`, `fire_base.jl` | `baf_crop` (crop CFT), `baf_peatf` (peatf>0), `dtrotr`/`TROTR` (tropical tree), `lfc` (land-use transition) — all `≡0` at Bow. **PEAT ROW DONE 2026-07-18 (§12): `baf_peatf` BIT-EXACT + found/fixed a real `wf2` double-count bug.** The unblock was a synthetic `peatf=0.5` in a copy of Bow's surfdata (peatf feeds ONLY CNFire → surgical, warm-CN), NOT a peatland site — ALL 21 domain surfdatas have peatf=0. Reference `clm_bgc_spinup/bow_ref_firepeat`, harness `scripts/fortran_parity_firepeat.jl`. **Bug: CTSM `HydrologyNoDrainageMod` accumulates the 0.17m `wf2` sums ON TOP of the 0.05m `wf` sums (no reset → top layers double-counted); Julia's `compute_wf!` omitted it → `baf_peatf` 7.7% high at the `max(wf2,0)` kink. Fixed via `compute_wf2!` (default byte-identical — `wf2_col` is fire-only).** STILL BLOCKED: `dtrotr` (tropical **deforestation** fire) needs `transient_landcover` (a `flanduse_timeseries`, absent from ALL domains → `dwt≡0`) AND a BGC restart at a tropical site (Aripuana/Leticia have `trotr>0.6` but are cold-CN); `baf_crop` shares A3's blocker — `_firea2016_crop_fire_kernel!` is reachable + gated on `itype>nc4_grass` but `≡0` with no crop CFT (docs/CROP_PARITY.md). | Mead (crop), ~~MerBleue (peat)~~ **peatf=0 everywhere → synthetic peatf at Bow**, Aripuana/Leticia (tropical, need transient LU) |
| A3 | **Crops** (`crop.jl` + crop phenology/allocation) + **crop N** (`n_fert!`, `n_soyfix!` — currently UNWIRED, correctly refused blind-wiring per #218) | `crop.jl`, crop phen/alloc, `n_dynamics.jl` | `use_crop=.true.` + surfdata with a crop CFT + crop spinup; diff crop C/N pools, GDD, harvest, fertilizer. **BLOCKED 2026-07-18 (docs/CROP_PARITY.md): the migrated Mead assets are an SP run (`use_crop/use_cn/irrigate=.false.`, surfdata `cft=2`, snow-only pdumps); build is `-bgc sp`. `use_crop=.true.` needs a crop-resolved (`cft`=crop CFT) surfdata — none on disk, no `inputdata`/raw crop datasets to build one, no global crop surfdata to subset. Julia side VERIFIED: crop accumulators (HUI/GDD) now live (SPVAL bug fixed #218), but `crop_phenology!` is a simplified cphase/bglfr setter — `plant_crop!`/`vernalization!`/harvest not driven; needs the reference to wire+validate.** | `domain_Cropland_Mead_USA` (SP-only — no crop CFT) |
| A4 | **Irrigation** | `irrigation.jl` | crop-bearing site + `irrigate=.true.`; at Bow `n_irrig_steps_left==0`, `irrig_rate` NaN (no crop CFT). ~~Driver call sites were empty stubs — verify wired.~~ **VERIFIED WIRED 2026-07-18 (docs/CROP_PARITY.md): `set_irrig_method!`/`calc_irrigation_needed!`/`calc_irrigation_fluxes!` are real, gated calls — the "empty stub" concern is resolved.** Behaviour still un-validatable: shares A3's crop-CFT-surfdata blocker. | `domain_Cropland_Mead_USA` (SP-only — no crop CFT) |

## B. Needs a non-soil LANDUNIT run

| # | Item | Julia module(s) | Fortran run | Notes |
|---|------|-----------------|-------------|-------|
| B1 | **Lake** | `lake_*.jl` (con/fluxes/hydrology/temperature) | lake-landunit column | **REFERENCE REGENERATED + harness re-run (2026-07-18).** The removed `clm_lake_run` h0 was reconstructed by `scripts/validation/gen_lake_ref.sh` — a pure namelist reconfig of the existing cesm.exe (NO SourceMods/rebuild): cold-start SP on the CTSM `test_inputs/lake/surfdata_lake100.nc` (PCT_LAKE=100), 48 hourly steps from 2003-01-01, based on the 2003 SP `clm_parity_run` staging that survives on the Drive migration copy. `fortran_parity_lake.jl` now runs green (rc=0): **TLAKE/LAKEICE thermodynamics track (~1%)**; the residual is the lake **SURFACE TURBULENT-FLUX solve** (FSH/EFLX_LH ~10%, TSA ~7–11%) — the same open residual documented historically (the FSA ~13× spikes are relative-diff artifacts on near-zero dawn/dusk solar). So: reference restored, harness green-runs, thermodynamics validated, surface-flux residual still open. See `docs/PARITY_COVERAGE_2026-07.md`. |
| B2 | **Urban** | `urban_*.jl` (albedo/fluxes/radiation) | urban-landunit column | Metal-validated; `urban_radiation` is called but no-ops in Bow |
| B3 | **Glacier (glcmec)** | glacier path | glcmec-landunit column | only surfdata/init coverage today |
| B4 | **Hillslope** | `hillslope_hydrology.jl` | hillslope routing config | unit-tested only |

## C. Needs a config-flag rerun

| # | Item | Julia module(s) | Fortran run | Notes |
|---|------|-----------------|-------------|-------|
| C1 | **Isotopes** (C13/C14) | `carbon_isotopes.jl`, `c_iso_flux.jl` | `use_c13/use_c14=.true.` rerun | **REFERENCE GENERATED + subsystem un-blocked (2026-07-18).** The Fortran `use_c13/c14=.true.` **spinup now runs** (`clm_bgc_spinup/bgc_ref_iso`, nstep window 4700-4725): the "restart-read crash" was CTSM's #2119 hard-abort ("cannot init c13 from a non-c13 restart"), bypassed by a case SourceMods patch (`scripts/validation/CNVegCarbonStateType.ciso_reseed_bypass.diff`) so the reseed-from-bulk logic runs — iso pools cold-start = bulk×ratio (leaf δ13C ≈ −28‰). Dumps carry 26 c13 + 26 c14 veg pools. `fortran_parity_isotopes.jl` rewritten to inject + diff them and **found 3 real port bugs, all FIXED**: (1) `use_c13/c14` never reached the veg facade config → the veg isotope carbon state was size-0 and the whole veg isotope path **silently no-op'd** (`clm_initialize.jl`); (2)+(3) the CIsoFlux litter/gap/harvest/gru p2c scatters used the `varpar` `i_litr_min`/`i_met_lit` **sentinel −9** → BoundsError once the state was live (`c_iso_flux.jl` + `cn_driver.jl` forward the runtime values). **NOW FIXED (2026-07-18, next PR):** the per-isotope `c_state_update{0,1,2,2h,2g,3}!` cascade is wired (each `CStateUpdate` runs 3× — bulk/c13/c14), so isotope fluxes are APPLIED to the veg pools; and the isotope psn→cpool input was repaired (`calc_gpp_mr_availc!` now receives the c13/c14 carbonflux instances — it was dropping to the bulk array, leaving `cpool_13` uninitialized/NaN once update0 ran). Result vs `bgc_ref_iso` n4712: post-step pool parity ~9e-6, `cpool_13/14` EXACT, and **δ13C discrimination parity <1e-4 ‰** (the definitive isotope-physics test, decoupled from bulk-C). The one-step increment |rel| hits the single-step bulk-C oracle floor on near-static veg pools (leafc net Δ ~1e-7 < the inherited summer bulk-C residual). C1 is now **VALIDATED (isotope physics) / bulk-C oracle floor (increment)**. See `docs/PARITY_COVERAGE_2026-07.md`. |
| C1a | **LUNA vcmax EOD update** (Rubisco-N optimum) | `luna.jl` `update_photosynthesis_capacity!` / `nitrogen_allocation!` | Fortran `Allocation` **internal** instrumentation (emit `Ncb`/`vcmx25_opt`, not just post `vcmx25_z`) | **KNOWN DIVERGE (2026-07-17):** `fortran_parity_luna_update.jl` shows Julia's post-EOD `vcmx25_z` runs **+5.2–5.6% high** (jmx25_z matches exactly). Constants byte-identical to `LunaMod.F90`; the residual is in the joint N-allocation optimizer's carboxylation pool `Ncb`. Not tuned. See `docs/PARITY_COVERAGE_2026-07.md` §"LUNA vcmax residual". |
| C2 | **CNDV** (dynamic vegetation) | `cndv.jl` | `use_cndv` rerun | Metal-validated only |
| C3 | **VOC / MEGAN** | `voc_emission.jl` | non-empty MEGAN compound list | not called in Bow (empty compound list) |
| C4 | **MIMICS decomposition** | `decomp_mimics.jl` | `decomp_method='MIMICS'` rerun | Bow uses CENTURY/BGC |

## D. FATES deeper validation

| # | Item | What a Fortran run must exercise |
|---|------|----------------------------------|
| D1 | Cold-start (#213, 27/27) + time-stepped 5-day dynamics (#217/#223, phase-attributed) are DONE. **Extended to 20 days (2026-07-18):** re-confirmed the committed 5-day oracle reproduces bit-for-bit on current `main` (the gated #197/#227 changes do not regress it), then generated a **fresh 480-step Fortran reference** with the surviving instrumented `cesm.exe` (no rebuild). It surfaced a divergence invisible at 5 days — by day 19 the port's `maxdbh` is **22 % low** and site carbon **5.5 % low** — **root-caused and attributed to the HOST**: the CLM top soil layer over-dries (`h2ovol1` 0.11–0.14 vs Fortran's 0.32–0.34, persistent from day 1 under bit-identical forcing) → FATES `btran_ft` correctly collapses to 0 on the mid-July dry days → PARTEH stature growth stalls. No FATES port bug; closing it is a CLM soil-hydrology task. See `scripts/validation/fates_fortran_parity/README.md` §"Extended-window validation". **Still open:** a **multi-year** FATES run (needed for the #197 screen's distinctive no-boom-bust value; the screen's cold start is already a subset of the Fortran-validated 14-PFT cold start), and a Fortran oracle for the **#227** boreal moisture BC (no clean short-horizon oracle — it is a harness-only prescription compensating a shared cold-start artifact, itself now Fortran-confirmed as a real over-drying bias). |

---

## Notes on prioritization

- **A1 (methane-as-source at MerBleue)** is the single highest-value scientific gap — it
  validates the half of the CH4 model that is the point of the CH4 model.
- **A2/A3/A4 (crop + fire crop/peat)** can share Fortran runs (a crop-bearing and a
  peatland run exercise several rows at once).
- **B (landunits)** are robustness/coverage rather than a specific science claim; batch
  them if a multi-landunit case is convenient.
- Several rows have an **existing `fortran_parity_*.jl` harness** (lake, isotopes) — the
  first action for those is to *run the existing harness* and see whether a dump already
  exists, before generating new ground truth.
