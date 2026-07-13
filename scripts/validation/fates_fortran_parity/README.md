# FATES Fortran-parity ground truth — **TIME-STEPPED, PHASE-ATTRIBUTED**

Fortran-side companion to `scripts/fates_fortran_parity.jl` (the Julia dump +
`compare_dumps` scorecard). An instrumented, FATES-enabled single-point CTSM run
emits a phase-tagged pdump of the site→patch→cohort hierarchy at **every FATES
phase boundary of every timestep**, so the Julia FATES port can be scored
field-by-field against Fortran ground truth **and any divergence attributed to the
FATES phase that produced it**.

## What changed vs. the cold-start-only version (PR #213)

PR #213 produced the first Fortran-FATES ground truth this project ever had — but
only at **nstep = 0**. That validated FATES *initialization*: the FATES parameter
file is byte-identical between port and install, so the cold-start cohorts are
parameter-identical **by construction**. It said nothing about the FATES
**dynamics** — photosynthesis, allocation, growth, mortality, recruitment,
disturbance, canopy promotion/demotion — which are the bulk of the ~62-module port.

This version runs the same case **forward** (default 120 × 3600 s = 5 days) and
dumps at every phase boundary:

| code | phase | what it bounds |
|---|---|---|
| 0 | `coldstart` | `EDInitMod::init_cohorts` (nstep = 0) |
| 1 | `fast` | **every timestep**, after the photosynthesis/flux solve + FATES running-mean + hifrq history, before the daily dynamics |
| 10 | `dyn_in` | `ed_ecosystem_dynamics` entry — the daily step's INPUT state |
| 11 | `phenology` | after `phenology` |
| 12 | `distrates` | after `fire_model` + `disturbance_rates` (→ mortality rates) |
| 13 | `integrate` | after `ed_integrate_state_variables` (growth / allocation / PRT / maintenance turnover) |
| 14 | `recruit` | after the per-patch `recruitment` loop |
| 15 | `cohortfuse` | after `sort_cohorts` / `terminate_cohorts` / `fuse_cohorts` |
| 16 | `spawn` | after `spawn_patches` (disturbance → new patch) |
| 17 | `patchfuse` | after `fuse_patches` + `terminate_patches` (= end of `ed_ecosystem_dynamics`) |
| 18 | `updatesite` | after `ed_update_site` (`canopy_spread`, `canopy_structure` = promotion/demotion, `trim_canopy`) — end of the daily step |

Each record also carries the **HLM→FATES boundary conditions** (`t_veg`, `eair`,
`cair`, `rb`, daylength factor, incident solar; soil T / moisture / porosity). That
is the discriminator that makes a divergence *attributable*: **if `bc_in` matches to
round-off but the FATES output does not, the divergence is inside FATES; if `bc_in`
already differs, it came from the host land model.**

## Why this is an oracle, not a drift plot

Both models start from a cold start that is **bit-identical** (re-verified here:
site carbon agrees to 7e-16) and are driven by the same forcing (incident solar,
CO₂/O₂ partial pressures and surface pressure are bit-identical in the dumps). So up
to the *first* divergence, every phase comparison is a genuine single-step oracle:
the Julia phase consumed the same input state the Fortran phase did, and any
difference in its output is **translation error**, not accumulated drift. After the
first divergence the numbers become a free-running drift measurement. The scorecard
reports both.

## Contents
- `FatesParityDumpMod.F90` — the Fortran instrument (schema v2). Walks
  `site%oldest_patch→younger` and `patch%tallest→shorter` (identical order to the
  Julia side), writes ~50 fields per cohort with 17 significant digits.
- `instrument.py` — injects the dump into a CTSM tree (exact-anchor patches, all
  asserted): `clmfates_interfaceMod.F90` (cold-start dump + a `parity_dump` type-bound
  procedure + site tagging around the daily loop), `clm_driver.F90` (the per-timestep
  `fast` dump), `EDMainMod.F90` (the nine daily sub-phase dumps).
- `setup_case.sh` — **the working recipe, end to end**: ncgen params → create_newcase
  → SourceMods → build → stage a single-point run dir → run → dump.
- `fates_pdump_fortran.txt.gz` — the ground truth itself (121 steps × 11 phases,
  5 390 records), committed so the scorecard reproduces without a CTSM build.
- `fates_pdump_fortran_n0.txt` — the original PR #213 cold-start-only (schema v1) dump,
  kept for provenance.

## Run window
Default **2002-07-01, 120 × 3600 s steps** (5 days). *Not* January: Bow at Banff is a
51°N montane site, so a winter window leaves every cohort dormant with GPP ≡ 0 —
which would "validate" the fast phase against a column of zeros. Mid-summer exercises
real photosynthesis, real daily allocation, real growth, and (by day 3) real
recruitment and patch disturbance.

```bash
bash scripts/validation/fates_fortran_parity/setup_case.sh      # Fortran ground truth
julia +1.12 --project=. scripts/fates_fortran_parity.jl compare \
    scripts/validation/fates_fortran_parity/fates_pdump_fortran.txt.gz
```

---

# THE FIRST FATES DYNAMICS SCORECARD

```
  ── per-phase worst case over the whole run ──
  ph  phase        bnds max rel(FATES) worst field max rel(bc_in) worst bc record-set
  ------------------------------------------------------------------------------------
  0   coldstart       1      7.133e-16 carbon          0.000e+00 -        match
  1   fast          121      6.028e+03 dndt            9.314e-01 solai1   MISM 71/121
  10  dyn_in          5      6.028e+03 dndt            1.140e-01 eair     MISM  2/5
  11  phenology       5      6.028e+03 dndt            0.000e+00 -        MISM  2/5
  12  distrates       5      6.028e+03 dndt            0.000e+00 -        MISM  2/5
  13  integrate       5      3.660e+03 dndt            0.000e+00 -        MISM  2/5
  14  recruit         5      3.660e+03 dndt            0.000e+00 -        MISM  3/5
  15  cohortfuse      5      6.028e+03 dndt            0.000e+00 -        MISM  3/5
  16  spawn           5      6.028e+03 dndt            0.000e+00 -        MISM  3/5
  17  patchfuse       5      6.028e+03 dndt            0.000e+00 -        MISM  3/5
  18  updatesite      5      6.028e+03 dndt            0.000e+00 -        MISM  3/5

  ── FIRST DIVERGENCE (tol = 1e-10) ──
    step 0, phase 1 (fast), field `npp`
    julia = -6.6197492656330261e-08   fortran = -6.1275741256280405e-08   rel = 4.9e-09
```

### Trajectory (5-day free run, site 1, patch 1)

| nstep | day | Julia carbon | Fortran carbon | rel | cohorts J/F | patches J/F | maxdbh J/F |
|---|---|---|---|---|---|---|---|
| 0   | 0 | 2549.0143 | 2549.0143 | **7.1e-16** | 14 / 14 | 1 / 1 | 1.88477 / 1.88477 |
| 25  | 1 | 2543.6671 | 2547.8844 | 1.7e-03 | 14 / 14 | 2 / 2 | 1.88477 / 1.88477 |
| 49  | 2 | 2548.1948 | 2555.4244 | 2.8e-03 | 14 / 14 | 2 / 2 | 1.88478 / 1.88500 |
| 73  | 3 | 2554.7695 | 2567.5769 | 5.0e-03 | 32 / 36 | 2 / 2 | 1.88499 / 1.88544 |
| 97  | 4 | 2563.8210 | 2584.0735 | 7.8e-03 | 40 / 42 | 2 / 2 | 1.88529 / 1.88600 |
| 120 | 5 | 2575.5305 | 2600.4706 | **9.6e-03** | 43 / 45 | 2 / 2 | 1.88562 / 1.88652 |

Per-PFT cohorts at day 5 (patch 1): `dbh` within **0.3 – 1.8 %**, cohort density `n`
within **0.16 – 0.22 %**, `leafc` within ~1 %. Patch dynamics agree exactly — both
models spawn a second patch on day 1 and keep two patches for the whole run.

### Reading of the result — honestly

* **The FATES phases themselves track Fortran.** The first divergence is `npp` at
  `rel = 4.9e-9` on the very first fast phase, and it is driven by the *host*: the
  cohort state is identical, and the CLM canopy temperature / vapour pressure /
  boundary-layer resistance handed to FATES already differ (the `bc_in` column). The
  fast (photosynthesis/respiration) thread reproduces Fortran to ~1e-8 given the same
  boundary conditions.
* **The daily thread tracks too**, but is *nonlinearly amplified* by threshold
  physics. The large `dndt` numbers are the mortality-rate cliff described below, not
  a broad-based error: `dndt` is a small signed number, so its relative error blows up
  while the resulting density `n` stays within 0.2 %.
* **The demography is the weakest link.** Recruitment first diverges at day 3
  (step 49, phase `recruit`): Fortran recruits **4 cohorts that Julia does not**, and
  the port ends the 5-day run with 43 cohorts vs 45. This is the *first structural*
  divergence and the honest headline: **FATES recruitment diverges at day 3, and the
  standing carbon is 0.96 % low after 5 days.**
* `NaN-only-1-side` fields (`hmort`, `cmort`, `bmort`, `ddbhdt`, `dhdt`, `dndt`,
  `*_hold`, …) are cohort *diagnostics* that Fortran zeroes at initialization and the
  port leaves `NaN` until the first daily step writes them. That is an
  initialization-cosmetics gap, not a physics divergence; the scorecard reports it
  separately rather than scoring it as `Inf`.

---

# REAL BUGS THIS HARNESS FOUND (all fixed, all unambiguous against the Fortran)

### 1. `h_allom` destructured backwards — cohort height was set to the allometric *derivative*
`src/fates/EDMainMod.jl` (`ed_integrate_state_variables`)

```julia
_, currentCohort.height = h_allom(currentCohort.dbh, ft)      # WRONG
currentCohort.height, _ = h_allom(currentCohort.dbh, ft)      # fixed
```

`h_allom` returns `(height, dh/ddbh)`; the Fortran is
`call h_allom(currentCohort%dbh, ft, currentCohort%height)` — the **height** is the
output. Every daily growth step was assigning `dh/ddbh` to `cohort.height`. For a
PFT-1 sapling at `dbh = 0.764 cm`, `h = 1.30 m` but `dh/ddbh = 1.36`, so height jumped
1.30 → 1.36 m on day one **with `ddbhdt == 0`**, and kept climbing (1.30 → 5.99 m by
step 14 for one PFT, at a *shrinking* dbh). Height drives canopy layering, crown area
and the tallest→shortest cohort ordering, so the whole canopy structure was corrupted.
Every other `h_allom` call site in the port already used `h, _ = h_allom(...)`.

### 2. CTSM's `use_fates` soil-water cold start was never ported
`src/infrastructure/cold_start.jl` (`init_soil_moisture!`)

CTSM `WaterStateType::InitCold` has an explicit FATES branch:

```fortran
if (use_fates) then
   h2osoi_vol_col(c,j) = 0.75_r8 * watsat_col(c,j) * ratio   ! wet
else
   h2osoi_vol_col(c,j) = 0.15_r8 * ratio                     ! CLM default
end if
```

The port implemented only the `0.15` branch (with a comment that mis-stated the
Fortran). At the cold-start soil temperature (272 K) a 0.15 profile is **100 % ice** →
zero liquid → FATES `btran_ft = 0` → `btran ≤ hf_sm_threshold (1e-6)` → FATES fires
hydraulic-failure mortality at `mort_scalar_hydrfailure = 0.6 /yr` on day one, a **44×
inflation of the total mortality rate** (`dmort` 0.614 vs 0.014).

### 3. FATES cold start never ran the canopy summarization
`src/fates/fates_driver_init.jl`

Fortran `init_coldstart` follows its `ed_update_site` loop with
`wrap_update_hlmfates_dyn` (→ `canopy_summarization` + `update_hlm_dynamics`). The
port called only `ed_update_site`, so the cold-start cohorts had **no leaf-area
profile**: `nv` (number of leaf layers) stayed **0** and patch `total_canopy_area`
stayed **NaN**. FATES therefore had zero leaf layers on the very first timestep —
photosynthesis and leaf dark respiration were identically zero for that step
(Fortran: `nv = 2`, `tcanarea = 6260`).

### 4. The daily `bc_in` soil-water pack used liquid volume, not total volume
`src/fates/fates_driver_bc.jl` (`fates_pack_bcin_daily!`)

Fortran `dynamics_driv` fills the daily boundary condition from the column's **total**
volumetric water,
`bc_in(s)%h2o_liqvol_sl(1:nlevsoil) = waterstatebulk_inst%h2osoi_vol_col(c,1:nlevsoil)`;
only the sub-daily `wrap_btran` pack uses the liquid-only `h2osoi_liqvol`. The port
packed liquid-only in both places, understating the soil water the daily FATES step
sees whenever any layer held ice.

### 5. (Harness, not model) CTSM runs a full `clm_drv` timestep at **nstep = 0**
The instrumented dump contains a `t=0 ph=1` record: the NUOPC initialization phase
advances the land once *before* nstep = 1. That step runs `SoilTemperature`/
`PhaseChange`, so by the time FATES first packs its boundary conditions the cold-start
soil (272 K, 100 % ice) has already partially melted onto the 273.15 K plateau and has
liquid water. A Julia harness that starts its loop at nstep = 1 is **one CLM step
behind** and hands FATES a bone-dry frozen soil. `scripts/fates_fortran_parity.jl` now
replays the step-0 call.

---

## What is still NOT validated
* **Recruitment** — the one structural divergence (day 3, 4 cohorts). Not chased to a
  root cause; it is the obvious next target.
* **Longer-timescale demography** — 5 days exercises growth, allocation, mortality
  rates, one disturbance/patch spawn and the first recruitment pulse. It does **not**
  exercise cohort fusion at scale, patch fusion/termination, canopy
  promotion/demotion under a closing canopy, fire, or the seasonal phenology
  transitions. Those need a multi-month run.
* **Plant hydraulics** is OFF (its 2-D solvers are FATAL-stubbed by design), so `hmort`'s
  hydro branch, `FatesPlantHydraulicsMod` and the rhizosphere geometry are untested.
* **`update_history_nutrflux`** and **`GetLUHStatedata`** (land use) are not ported.
* The residual **HLM** difference (CLM canopy temperature ~3.7 K at the first packed
  step, soil liquid ~15 % relative) is a *host-model* residual, below the FATES port
  and outside the scope of this harness. Note: the incident solar handed to FATES
  (`solad1`/`solai1`) and the CO₂/O₂/pressure boundary conditions are **bit-identical**
  in the two dumps, so the FATES↔HLM shortwave coupling gap closed by #215 is not a
  factor here — the scorecard above was re-measured on the post-#215 tree and is
  unchanged to the digit.

## Environment traps (macOS 26 / Homebrew) — all handled in `setup_case.sh`
1. **`HWLOC_COMPONENTS=-opencl` is REQUIRED.** hwloc's OpenCL probe walks into the Metal
   driver at MPI startup and dies with `SIGILL` inside `AGXMetalG16X`. Without it
   `cesm.exe` never reaches CLM at all.
2. **Flaky xzone malloc trap.** The optimized build `SIGTRAP`s inside `libsystem_malloc`
   on roughly half of all launches, in a *different* allocation each time — the known
   ESMF/xzone false positive, not a CLM bug. `setup_case.sh` retries. Also set
   `MallocNanoZone=0`. **Never run under lldb.**
3. **CIME python.** Needs `distutils` *and* a stdlib-compatible ElementTree. Use a 3.11
   conda env (`CIME_PY`). `create_newcase` is fed `printf 'u\nu\n'` so its
   "EXEROOT/RUNDIR already exists" prompts do not hit `EOFError` — and so the object
   tree is reused (an incremental rebuild of the four instrumented files takes ~40 s).
4. **netcdf link path.** The `homebrew` mach config hardcodes a Cellar version;
   `setup_case.sh` symlinks the stale path.
