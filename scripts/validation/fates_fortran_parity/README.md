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

# THE FATES DYNAMICS SCORECARD

**Current (post-#223 — the recruitment divergence root-caused and fixed).**

```
  ── per-phase worst case over the whole run ──
  ph  phase        bnds max rel(FATES) worst field max rel(bc_in) worst bc record-set
  ------------------------------------------------------------------------------------
  0   coldstart       1      7.133e-16 carbon          0.000e+00 -        match
  1   fast          121      2.069e+02 dndt            9.314e-01 solai1   MISM 47/121
  10  dyn_in          5      2.069e+02 dndt            9.032e-02 h2ovol1  MISM  1/5
  11  phenology       5      2.069e+02 dndt            0.000e+00 -        MISM  1/5
  12  distrates       5      2.069e+02 dndt            0.000e+00 -        MISM  1/5
  13  integrate       5      2.069e+02 dndt            0.000e+00 -        MISM  1/5
  14  recruit         5      2.069e+02 dndt            0.000e+00 -        MISM  1/5
  15  cohortfuse      5      2.069e+02 dndt            0.000e+00 -        MISM  2/5
  16  spawn           5      2.069e+02 dndt            0.000e+00 -        MISM  2/5
  17  patchfuse       5      2.069e+02 dndt            0.000e+00 -        MISM  2/5
  18  updatesite      5      2.069e+02 dndt            0.000e+00 -        MISM  2/5

  ── FIRST DIVERGENCE (tol = 1e-10) ──
    step 0, phase 1 (fast), field `npp`
    julia = -6.6197492656330261e-08   fortran = -6.1275741256280405e-08   rel = 4.9e-09
```

### Trajectory (5-day free run, site 1)

| nstep | day | Julia carbon | Fortran carbon | rel | cohorts J/F | patches J/F | maxdbh J/F |
|---|---|---|---|---|---|---|---|
| 0   | 0 | 2549.0143 | 2549.0143 | **7.1e-16** | 14 / 14 | 1 / 1 | 1.88477 / 1.88477 |
| 25  | 1 | 2548.0850 | 2547.8844 | 7.9e-05 | 14 / 14 | 2 / 2 | 1.88477 / 1.88477 |
| 49  | 2 | 2560.9846 | 2555.4244 | 2.2e-03 | 14 / 14 | 2 / 2 | 1.88512 / 1.88500 |
| 73  | 3 | 2572.9907 | 2567.5769 | 2.1e-03 | 36 / 36 | 2 / 2 | 1.88557 / 1.88544 |
| 97  | 4 | 2588.0635 | 2584.0735 | 1.5e-03 | 36 / 42 | 2 / 2 | 1.88612 / 1.88600 |
| 120 | 5 | 2605.8573 | 2600.4706 | **2.1e-03** | 44 / 45 | 2 / 2 | 1.88670 / 1.88652 |

(Was, before #223: day-5 carbon **9.6e-03 low**, 43/45 cohorts, and 4 cohorts
Fortran recruited that the port did not.)

**Recruitment now agrees structurally.** Every daily `recruit` phase produces the
**same 22 cohorts** (11 PFTs × 2 patches) in both models — the day-2 "4 missing
cohorts" are gone. Per-PFT persistent cohorts at day 5 (patch 1): `dbh` within
**0.01 – 0.93 %** (was 0.3 – 1.8 %), cohort density `n` within **0.5 %** (0.00 % for
11 of 14 PFTs), `leafc` within **0.06 – 1.8 %**. Patch dynamics agree exactly.

### Reading of the result — honestly

* **The FATES phases track Fortran.** The first divergence is `npp` at `rel = 4.9e-9`
  on the very first fast phase, and it is driven by the *host*: the cohort state is
  identical and the CLM canopy temperature / vapour pressure / boundary-layer
  resistance handed to FATES already differ (the `bc_in` column). The fast
  (photosynthesis/respiration) thread reproduces Fortran to ~1e-8 given the same
  boundary conditions.
* **The residual is now the HOST, and the dump proves it.** With the root profile
  fixed (see bug 6) FATES's entire BTRAN comes from soil layer 1, so the FATES water
  stress is a near-direct readout of one host number — and that number differs: the
  Julia CLM top layer leaves the 273.15 K melt plateau ~4 h early and holds
  **`h2ovol1` 0.31 vs 0.35** (day 1: 0.25 vs 0.17). That is now the worst `bc_in`
  field in the daily thread (`dyn_in`, 9.0e-2), and it is what remains of the
  daily-carbon difference: the port's daily GPP is now ~25 % **high** where it used to
  be ~45 % **low**. The FATES side of that boundary is exonerated; the CLM
  soil-temperature / phase-change residual owns it.
* **The daily thread is nonlinearly amplified by threshold physics.** The large
  `dndt` numbers are the mortality-rate cliff, not a broad-based error: `dndt` is a
  small signed number, so its relative error blows up while the resulting density `n`
  stays within 0.5 %.
* **The first STRUCTURAL divergence moved from `recruit` to `cohortfuse` (day 3).**
  Both models recruit 58 cohorts; Fortran's `sort/terminate/fuse` leaves 43 and the
  port's leaves 38. That difference lives *entirely in patch 2* — the 0.03 %-of-area
  patch spawned by disturbance, whose cohort densities are O(1e-7) and sit right on
  the `terminate_cohorts` threshold. Patch 1 (99.97 % of the area) is structurally
  **identical** through the whole run (36 → 25 cohorts, same PFT multiset). The
  fusion criterion itself was diffed against the Fortran and matches; the divergence
  is the threshold sensitivity of a tiny patch, downstream of the host GPP bias.
* `NaN-only-1-side` fields (`hmort`, `cmort`, `bmort`, `ddbhdt`, `dhdt`, `dndt`,
  `*_hold`, …) are cohort *diagnostics* that Fortran zeroes at initialization and the
  port leaves `NaN` until the first daily step writes them. That is an
  initialization-cosmetics gap, not a physics divergence; the scorecard reports it
  separately rather than scoring it as `Inf`.

### Was recruitment stochastic? No.

`EDPhysiologyMod::recruitment` is fully deterministic — `cohort_n` is
`min(mass_avail/mass_demand)` over the elements, gated on `cohort_n > min_n_safemath`.
There is no random draw anywhere in the seed → germination → recruitment chain (the
FATES stochasticity that does exist, in the fire/lightning stream, is not on this
path). Bit-parity of recruitment is therefore a legitimate target, not an unreachable
one, and any disagreement is a real bug.

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

### 6. THE RECRUITMENT DIVERGENCE: FATES's root profile was never told how deep the roots go
`src/fates/fates_driver_bc.jl` (`fates_pack_bcin_daily!`), `src/driver/clm_driver.jl`

**This is the root cause of the day-3 recruitment divergence.** It is NOT in the
recruitment code — `recruitment`, `SeedGermination`, `SeedDecay` and `SeedUpdate` are
all faithful and all live. The **seed pool entering recruitment already differed**:
`cohort_n = mass_avail/mass_demand` with `mass_avail = patch%area * litter%seed_germ(ft)`,
and the port's `seed_germ` was ~2× low for *every* PFT. Walking that upstream through
the dump: seed pool ← `seed_prod` (2–27× low on **day 1**) ← daily `npp_acc`
(**35–48 % low**) ← instantaneous `gpp` (37–63 % low at every daylight step, while
`solad1`/`solai1` were bit-identical) ← **`btran_ft`**, roughly *half* of Fortran's and
plateauing at odd values (0.43, 0.62) where Fortran reached exactly 1.0. The 4 "missing
cohorts" were just the 4 PFT×patch slots whose starved seed pool had not yet crossed
`min_n_safemath` — a symptom, not the disease.

`bc_in%max_rooting_depth_index_col` is the culprit. Fortran sets it **every day**, in
the daily pack:

```fortran
! clmfates_interfaceMod.F90 :: dynamics_driv
this%fates(nc)%bc_in(s)%max_rooting_depth_index_col = &
     min(nlevsoil, active_layer_inst%altmax_lastyear_indx_col(c))
```

The port set it **once**, at init, to `nlevsoil`, and never updated it — the classic
"ported then never wired" bug. It is the `max_nlevroot` argument to
`set_root_fraction`, which normalizes the root profile over layers `1..max_nlevroot`.
And CLM **cold-starts `altmax_lastyear_indx_col = 0`** (`ActiveLayerMod::InitCold`), so
Fortran passes `max_nlevroot = 0` — and `set_root_fraction`'s normalization correction
(`corr_id = maxloc(root_fraction)` on an all-zero profile → index 1) drops the entire
**1.0 of root fraction into soil layer 1**. Fortran's FATES therefore draws water only
from the top (thawed, wet, July) layer and gets `btran_ft = 1.0`. The port spread its
roots over the whole column, most of which is still frozen at this montane site, so
every deep layer failed `check_layer_water` and its root fraction was silently deleted
from BTRAN. Hence the ~50 % water-stress deficit, hence the halved GPP, hence the
starved seed bank, hence the missing recruits.

### 7. FATES's soil suction was computed against total porosity, not *effective* porosity
`src/fates/fates_driver_bc.jl` (`fates_pack_bcin_btran!`, `fates_pack_bcin_daily!`)

Same chain, same sign, independent bug. Fortran (`wrap_btran` and `dynamics_driv`,
both) computes the suction it hands FATES as

```fortran
s_node = max(h2osoi_liqvol(c,j)/eff_porosity(c,j), 0.01_r8)   ! EFFECTIVE porosity
call soil_water_retention_curve%soil_suction(c,j,s_node, soilstate_inst, smp_node)
this%fates(nc)%bc_in(s)%smp_sl(j) = smp_node
```

The port divided by **`watsat`** (total porosity) and additionally clamped `s_node` to
`≤ 1`. In a partially-frozen layer `eff_porosity` is the *unfrozen* pore space, so a
layer whose remaining pores are full of liquid is at near-zero suction and is **not**
water-stressed; dividing by the total porosity instead reports that same layer as bone
dry, clamps its suction at `smpsc`, and zeroes its BTRAN contribution. Fixed to the
Fortran formula, on the Fortran's layer set: the suction is now computed **only on
`active_suction_sl` layers**, which required wiring two more dead pieces —
`bc_in%filter_btran` (never set anywhere in the port; Fortran sets it from the
exposed-veg column filter and fills non-exposed columns with `-999`) and
`get_active_suction_layers` (ported, never called).

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
* ~~**Recruitment**~~ — **DONE (#223).** Root-caused (bugs 6 + 7 above), fixed, and
  re-scored: the recruit phase is now structurally identical (same 22 cohorts every
  day, 11 PFTs × 2 patches), and the day-5 carbon drift fell 9.6e-3 → 2.1e-3.
  Recruitment is deterministic in the Fortran, so this is a real oracle, not a
  distributional one.
* **`terminate_cohorts` / `fuse_cohorts` on a marginal patch** — the new first
  structural divergence (day 3, phase `cohortfuse`: 58 → 43 Fortran vs 58 → 38 port),
  confined to the tiny disturbance-spawned patch 2 whose cohort densities are O(1e-7),
  right on the termination threshold. The fusion *criterion* was diffed against the
  Fortran and matches line-for-line, and patch 1 (99.97 % of the area) is structurally
  identical — but the termination thresholds on a near-empty patch are not yet
  independently validated, and this sits downstream of the host GPP bias, so it cannot
  be cleanly attributed until that bias is closed.
* **The host CLM top-soil-layer phase change** is now the binding residual on the whole
  FATES daily thread (see "Reading of the result"). It is a *CLM* problem, not a FATES
  one, but it is what stands between this scorecard and 1e-8 agreement — because at
  cold start FATES's entire root profile sits in layer 1, so layer 1's liquid water is
  the FATES water stress.
* **Longer-timescale demography** — 5 days exercises growth, allocation, mortality
  rates, one disturbance/patch spawn and the first recruitment pulses. It does **not**
  exercise cohort fusion at scale, patch fusion/termination, canopy
  promotion/demotion under a closing canopy, fire, or the seasonal phenology
  transitions. Those need a multi-month run — and note that a multi-month run crosses
  the first year boundary, which is when `altmax_lastyear_indx_col` finally becomes
  non-zero and the FATES root profile stops being a single layer. That transition is
  now wired (bug 6) but is itself unvalidated.
* **Plant hydraulics** — *superseded 2026-10-29: the 1-D path is now validated, see
  §"PLANT HYDRAULICS" at the end.* (Originally: OFF; its 2-D solvers are FATAL-stubbed by design), so `hmort`'s
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

---

## Appendix: the attribution chain that found bugs 6 + 7

The phase-tagged dump makes this mechanical. Each arrow is a dump field that
localizes the previous one:

```
  4 cohorts Fortran recruits that the port does not     (recruit, day 2)
    └─ NOT the recruitment decision: `recruitment`, `SeedGermination`, `SeedDecay`,
       `SeedUpdate` all diffed clean against the Fortran, and all are live (call
       sites verified, not just definitions).
    └─ mass_avail = patch%area * seed_germ(ft)   →  the SEED POOL was ~2x low
         └─ seed_prod  2–27x low on DAY ONE            (COHORT.seed_prod)
             └─ npp_acc 35–48 % low                    (COHORT.npp_acc @ dyn_in)
                 └─ gpp   37–63 % low, EVERY daylight step   (COHORT.gpp @ fast)
                     └─ solad1/solai1 BIT-IDENTICAL     (PATCH bc_in)  → not the light
                     └─ btran_ft ~HALF, plateaus at 0.43/0.62 vs Fortran's 1.0
                         └─ btran_ed! itself diffed clean → the BC's are wrong
                             ├─ bug 6: max_rooting_depth_index_col never updated
                             └─ bug 7: suction vs watsat instead of eff_porosity
```

The single most useful diagnostic was `btran_ft` **plateauing below 1.0**: BTRAN is a
root-fraction-weighted sum of per-layer resistances, each capped at 1, so a plateau at
0.62 says "38 % of my root profile is in layers that are contributing nothing" — which
is a *root profile* statement, not a *soil moisture* one. That is what pointed at
`set_root_fraction`'s `max_nlevroot` rather than at the hydrology.

---

# EXTENDED-WINDOW (20-day) VALIDATION + D1 STATUS (2026-07-18, backlog D1)

The scorecard above is the committed **5-day** (120-step) oracle. This section records
the **D1 "deeper FATES validation"** pass: (1) re-confirming that oracle still holds on
current `main` after the gated FATES config work landed (#197 fixed-biogeog screen, #227
boreal cold-start moisture BC), and (2) a **fresh 20-day (480-step) Fortran run** that
extends the demographic window and surfaces a divergence the 5-day window is blind to.

## 1. The committed 5-day oracle still holds on current `main`

Re-ran `fates_fortran_parity.jl compare fates_pdump_fortran.txt.gz` on `main` (post
#197/#227). The scorecard reproduces the committed one **to the digit**: cold-start
`carbon` 7.133e-16 (match); first divergence `npp` rel **4.9218e-09** at step 0 phase 1
(host-driven, identical bytes); per-phase `dndt` 2.069e+02 / `fast` MISM 47/121 — every
line unchanged. So the two gated changes did **not** regress the validated dynamics: #227
is **scripts-only** (a prescribed growing-season soil-moisture BC in
`fates_multisite_validation.jl`, gated on `SiteCfg.soil_h2o`; not on this Bow path), and
#197's `src/` edits are gated behind `fates_biogeog_screen` (default `:none` → the
all-14-PFT NBG cold start is byte-identical, which this harness runs).

## 2. Fresh 20-day Fortran run — a HOST soil-moisture divergence, amplified

Re-ran the **surviving instrumented FATES `cesm.exe`** (`fates_parity_1pt`, no rebuild —
the 4-file SourceMods build from #217/#223 is intact) with `stop_n = 480` (20 days) and
scored the port against it at `FATES_PARITY_STEPS=480`. Cold start still 7.133e-16; the
`fast`/`npp` first divergence and the day-3 `cohortfuse` marginal-patch-2 structural
divergence are unchanged. But the **site-level trajectory diverges over the longer
window** in a way 5 days cannot show:

| day | Julia carbon | Fortran carbon | relC | Julia maxdbh | Fortran maxdbh | relDBH |
|---|---|---|---|---|---|---|
| 3  | 2588.06 | 2584.07 | 1.5e-3 | 1.8861 | 1.8860 | 6.7e-5 |
| 12 | ~2660   | ~2660   | ~1e-3  | 1.8900 | 1.8906 | 3.4e-4 |
| 19 | 2730.31 | 2888.45 | **5.5e-2** | 1.9754 | 2.5199 | **2.2e-1** |

`gpp`/`npp` stay tight the whole run (worst rel ~4e-4). The break is **stature growth**:
around **day 13** Fortran's largest cohort — a **PFT-13 (stress-deciduous, grass-form:
`structc=0`, height decoupled from dbh) cohort in patch 1 at full density n≈2000** —
enters vigorous dbh growth (1.89 → 2.52 by day 19). The port's identical cohort **stalls**
its dbh/leafc/storec around days 11–13, then only partly recovers.

### Attribution — it is the HOST (CLM top-soil-layer moisture), not the FATES port

The phase-tagged `bc_in` makes this mechanical, and it is the **same class** as the 5-day
residual noted above ("the host CLM top-soil-layer phase change is the binding residual"),
now amplified:

```
  PFT-13 patch-1 dbh STALLS days 11–13   (COHORT.dbh @ updatesite)
    └─ gpp_acc COLLAPSES: day 12 rel 0.92, day 13 = 0.000   (COHORT.gpp_acc @ dyn_in)
        └─ leaf status STAYS 2 (leaf-on)  → NOT a phenology leaf-off
        └─ btran_ft → 0.0 on days 11–13   (PATCH.btran_ft @ dyn_in; Fortran stays 1.0)
            └─ h2ovol1 (top soil layer handed to FATES) J 0.11–0.14 vs F 0.32–0.34
                — a 3× dry-down, PERSISTENT from day 1 (day1 J 0.302 vs F 0.350)
                └─ incident solar BIT-ALIGNED: daily-integrated solad1 rel = 0.00
                   EVERY day (forcing-misalignment ruled out — the #233 lesson)
```

At cold start FATES's entire root profile sits in soil layer 1 (README bug 6), so layer
1's liquid water **is** the FATES water stress. The port's CLM top layer over-dries under
bit-identical forcing (soil-evaporation / phase-change / drainage residual, present from
day 1); on the mid-July dry days 11–13 it drops below wilting, `btran_ft` correctly
collapses to 0, GPP → 0, and PARTEH's phase-3 stature-growth gate (`carbon_balance > 0`)
is never entered, so dbh freezes. **This is a CLM host divergence: the `bc_in` soil
moisture FATES is handed already differs 3×, and `btran_ft = 0` is the correct FATES
response.** No FATES-port mis-port was found — closing it is a CLM soil-hydrology task
(out of FATES scope; there is a plausible GPP↔transpiration↔soil-moisture feedback seeded
by the documented early ~25 %-high GPP, but the primary driver is the top-layer moisture
bias). Fix attempts belong to the host hydrology, not `src/fates/`.

> **CORRECTION (D2, 2026-07-18): it was NOT a host soil-hydrology mis-port — it was a
> `use_bedrock` HARNESS↔REFERENCE namelist mismatch.** See section 4 below. The
> soil-hydrology port (`drainage!`, `clamp_zwt_to_bedrock!`, `init_bedrock!`, the `zwt`
> cold start) is faithful; with the harness config matched to the Fortran case the
> over-drying, the `btran_ft=0` collapse, and the day-13 stature-growth stall all vanish.

**Reproduce** (surviving build, no CTSM rebuild needed):
```bash
# Fortran 20-day reference (edit stop_n=480 in the run dir's nuopc.runconfig):
cd <fates_parity_1pt run dir> && sed -i '' 's/stop_n = 120/stop_n = 480/' nuopc.runconfig
MallocNanoZone=0 HWLOC_COMPONENTS=-opencl <EXEROOT>/cesm.exe   # retry on xzone trap
# Julia side, scored against it:
FATES_PARITY_STEPS=480 julia +1.12 --project=. scripts/fates_fortran_parity.jl \
    compare <run dir>/fates_pdump_fortran.txt
```
(The 480-step dump is ~8 MB gzipped and is **not** committed; the 120-step
`fates_pdump_fortran.txt.gz` remains the committed anchor.)

## 3. #227 (boreal cold-start) and #197 (fixed-biogeog screen) vs Fortran

Both are **validated internally only**; the committed ground truth (unscreened Bow, 51°N)
does **not** cover either, and neither has a clean short-horizon Fortran oracle:

* **#227 — boreal cold-start moisture prescription.** A **harness-only** BC (no `src`
  change): it re-asserts growing-season root-zone liquid water for a *thawed*-but-drained
  cold-start boreal column, compensating for a freeze/thaw/drain **short-spin artifact**
  that a Fortran run cold-started the same way would *share* (a spun-up boreal soil holds
  field-capacity liquid; a cold start does not). Because the fix deliberately **overrides**
  soil initialization, there is no independent Fortran oracle for it short of a spun-up
  boreal FATES restart, which does not exist. The 20-day Bow result here is the same
  phenomenon from the other side: the *unprimed* Bow column also over-dries its top layer
  vs Fortran — i.e. #227 is treating a real, now-Fortran-confirmed CLM cold-start
  soil-moisture bias, not inventing water.

* **#197 — fixed-biogeog PFT screen (`:drop_cold_deciduous`).** Gated cold-start seeding
  that drops 5 of the 14 default PFTs via the ported `use_this_pft` machinery — the **same
  mechanism** FATES's `use_fixed_biogeog` uses. Its cold start is therefore a **subset** of
  the already-Fortran-validated all-14-PFT cold start (27/27 at nstep 0), and the demographic
  engine it feeds is the one validated bit-for-bit here. Its *distinctive* value — a
  bounded, no-boom-bust **multi-year** stand — only manifests over a multi-year run, which
  is beyond tractable Fortran-side wall time on this box. A dedicated screened Fortran run
  (`use_fates_fixed_biogeog=.true.` with a matched `pft_areafrac`) is the remaining item,
  and carries a real config-crosswalk risk (the #233/#240/#243 "harness input, not port"
  trap) that must be controlled before any divergence is trusted.

---

# D2 ROOT CAUSE (2026-07-18): the top-layer over-drying was a `use_bedrock` harness mismatch

The D1 section (above) isolated the over-drying to the *host* CLM top-soil moisture and,
correctly, exonerated `src/fates/` — but it left the CLM cause open ("a CLM soil-hydrology
task"). D2 closes it, and the answer is **not a soil-hydrology mis-port**: it is the
`#233/#240/#248` trap the D1 note itself flagged — a **namelist mismatch between the Julia
harness and its Fortran reference case**.

## The measurement that found it

A new **per-column water-budget dump** was added on both sides (`SoilBudgetDumpMod.F90` on
the Fortran side, injected by `instrument.py` at the same "fast" point as the FATES dump —
which is *after* `HydrologyDrainage`; `scripts/soil_budget_probe.jl` on the Julia side,
reusing the harness build/step loop). It dumps, per soil column per step, the per-layer
`h2osoi`, the water table `zwt`, the SL14 resistance, and **every flux that sets the top
layer**: infiltration, surface runoff, **subsurface drainage**, soil evaporation,
transpiration and per-layer root extraction.

Diffing the 20-day budget layer-by-layer, the diverging process was unambiguous:

```
day |  hvol1 (top layer)  |   qflx_drain (mm/day)   |  zwt (m)
    | Julia(bug) Fortran  |  Julia(bug)   Fortran   | Julia  Fortran
  3 |   0.236    0.348    |    15.48       0.00      |  2.28   8.60
  8 |   0.170    0.328    |    17.21       0.00      |  2.28   8.60
 13 |   0.122    0.330    |     9.40       0.00      |  2.28   8.60
```

**Fortran drains ZERO every day; the port drains 3–17 mm/day** — ~150 mm over the window,
which bleeds the whole column (all three top layers dry together), and the thin top layer
follows fastest. (Fortran's evaporation and transpiration are actually *higher* — it is the
wetter model — so those are exonerated: the sole cause is the spurious drainage.)

## The chain

`qflx_drain` in `SoilHydrologyMod::Drainage` is `rsub_top = imped·rsub_top_max·exp(-zwt/hkdepth)`
— exponentially sensitive to the water table depth. Fortran's `zwt = 8.60 m` makes
`exp(-8.60/0.4) ≈ 5e-10` → no drainage. The port's `zwt = 2.28 m` makes `exp(-2.28/0.4) ≈ 3e-3`
→ ~10 mm/day. Why 2.28 vs 8.60?

* The Fortran case's `lnd_in` (CLM `build-namelist` for `I2000Clm50FatesRs`) sets
  **`use_bedrock = .false.`** → `initVerticalMod` forces `nbedrock = nlevsoi` (`zi = 8.6 m`),
  so the per-step `zwt→bedrock` clamp keeps the water table at 8.6 m.
* `clm_initialize!` **defaults `use_bedrock = true`**, and `scripts/fates_fortran_parity.jl::build()`
  left it unset — so the port read the Bow surfdata's shallow `zbedrock ≈ 2.28 m`
  (`nbedrock = 12`), and `clamp_zwt_to_bedrock!` (faithful to CTSM, called every step) dragged
  `zwt` from 8.6 m down to 2.28 m — igniting the drainage.

## The fix (harness config, not `src/`)

`build()` now passes **`use_bedrock=false`** to match the reference `lnd_in`. That single
line makes `zwt` hold 8.6 m, `qflx_drain ≡ 0` like Fortran, and the top-layer `h2ovol1`
track the Fortran ground truth to **~0.01–0.02** (was a 3× dry-down). Because the fix is a
harness namelist and touches **no `src/` soil-hydrology path**, the shared Bow water balance,
the 69/69 multi-biome scorecard and `longhorizon_conservation` are structurally unaffected;
the 5-day FATES scorecard is unchanged-to-improved (cold start still `7.133e-16`; the
worst `dyn_in` `bc_in` field moved off `h2ovol1`).

**This does not invent water** — it is the same phenomenon as `#227` seen from the other
side: the port's soil-hydrology is faithful, and once its cold-start config matches the
reference's, the Bow column retains water exactly as Fortran does.

**Reproduce:**
```bash
# the divergence (harness default before the fix):
USE_BEDROCK=true  julia +1.12 --project=. scripts/soil_budget_probe.jl 480
# matched to the Fortran case (after the fix):
USE_BEDROCK=false julia +1.12 --project=. scripts/soil_budget_probe.jl 480
# Fortran side (surviving fates_parity_1pt build; SoilBudgetDumpMod now injected):
#   bash scripts/validation/fates_fortran_parity/setup_case.sh   (SKIP_BUILD=1 to reuse)
#   -> writes soil_budget_fortran.txt in the run dir
```

---

# D3 (2026-07-19): full flag audit — harness vs. reference `lnd_in`, flag by flag

Several CLM driver defaults moved in the days after #247 was measured (#252
`use_bedrock`, #259 `use_aquifer_layer`, #225 `h2osfcflag`, #265
`create_crop_landunit`, #267–#269 `use_luna`, #273 `use_hydrstress`). Because
`scripts/fates_fortran_parity.jl::build()` **inherits** most of them rather than
setting them, every one of those moves silently re-specifies this harness. A
namelist mismatch of exactly this kind has been the actual cause 8 times in this
repo (#233, #240, #243, #248, #251/#252 …), so the audit is run **before** any
number is trusted.

Reference = `<run dir>/lnd_in` from the surviving `fates_parity_1pt` case
(`I2000Clm50FatesRs`, CLM `build-namelist`). Julia = the flags actually resolved by
`clm_initialize!(use_fates=true, use_bedrock=false)` as `build()` calls it, read back
out of `varctl` / `SoilHydrologyData` after init.

| flag | Fortran `lnd_in` | Julia (resolved) | verdict |
|---|---|---|---|
| `use_fates` | `.true.` | `true` | match |
| `use_cn` | `.false.` | `false` | match |
| `use_crop` | `.false.` | `false` | match |
| `create_crop_landunit` | `.false.` | `false` | match — `clm_initialize.jl:159` sets `!use_fates`, so #265 lands on the correct side here |
| `use_bedrock` | `.false.` | `false` | match — set explicitly by `build()` (#251/#252); with #252's conditional default (`use_fates` ⇒ `.false.`) the explicit pass is now belt-and-braces |
| `use_aquifer_layer` | derived `.false.` (`lower_boundary_condition = 2` = `bc_zero_flux`; `SoilWaterMovementMod.F90:221-236`) | `false` (#259) | match |
| `h2osfcflag` | not in `lnd_in` ⇒ CTSM default `1` (`SoilHydrologyType.F90:339`) | `1` (#225) | match |
| `use_luna` | `.false.` | `false` | match — CTSM `endrun`s on LUNA+FATES (`controlMod.F90:505`), mirrored at `control.jl:104`; #267's `nothing` default resolves to `false` under FATES |
| `use_hydrstress` | `.false.` | `false` | match — also mutually exclusive with FATES (`control.jl:102`) |
| `use_fun` | `.false.` | n/a (no `varctl` field; FATES path never sets `bgc_vegetation.config.use_fun`) | match by construction |
| `use_nitrif_denitrif` | `.true.` | `true` | match (inert at `use_cn=false`) |
| `use_subgrid_fluxes` | `.true.` | `true` | match |
| `irrigate` | `.false.` | `false` | match |
| `use_cndv` / `use_c13` / `use_c14` | `.false.` | `false` | match |
| all 13 `use_fates_*` sub-flags | `.false.` (incl. `fixed_biogeog`, `nocomp`, `sp`, `planthydro`, `cohort_age_tracking`, `tree_damage`, `ed_st3`, `ed_prescribed_phys`, `inventory_init`, `luh`, `lupft`, `potentialveg`) | all `false` | match — in particular `use_fates_fixed_biogeog=.false.` confirms the reference is the **unscreened** 14-PFT case #197 is *not* covered by |
| `use_lch4` | `.false.` | `varctl.use_lch4 = true` **(≠)**, but `CLMDriverConfig.use_lch4 = false` | **benign, latent trap** — see below |

**Verdict: no active mismatch.** The one discrepancy is cosmetic *on this path* but is
a real latent trap worth recording: `clm_initialize!`'s kwarg default is
`use_lch4=false` and it is what drives allocation and `CLMDriverConfig`, but the
kwarg is **never written back to the `varctl` global**, which keeps its own
independent default of `true` (`varctl.jl:110`). Nothing on the FATES path reads
`varctl.use_lch4` (only `lake_temperature.jl:1158` and `control.jl`'s validator do,
and this site has no lake column), so the 20-day comparison is unaffected — but any
future harness that reads the global instead of the config would silently run a
different model. Same class as #252/#259: a struct default that is CTSM's *code
fallback* rather than its *namelist* value.

---

# D3 (2026-07-19): the 20-day comparison RE-RUN on the corrected footing

#247 measured the 20-day divergence with the harness running `use_bedrock=true`
against a reference that runs `.false.` — the mismatch #251/#252 then found and fixed.
This section re-measures the same 480-step window with the harness now matched
(flag-audited above, no active mismatch), against the **same** surviving Fortran
reference dump (`<run dir>/fates_pdump_fortran.txt`, `stop_n = 480`).

```bash
FATES_PARITY_STEPS=480 julia +1.12 --project=. scripts/fates_fortran_parity.jl \
    compare <run dir>/fates_pdump_fortran.txt
```

## Result: the over-drying is GONE and the growth stall went with it

Site state at phase 18 (`updatesite`, end of each daily step), one row per day:

| day | C julia | C fortran | relC | dbh julia | dbh fortran | relDBH | h2ovol1 J/F | ncoh J/F |
|---|---|---|---|---|---|---|---|---|
| 1  | 2548.10 | 2547.88 | +8.3e-05 | 1.8848 | 1.8848 | 0.0 | 0.3995 / 0.3979 | 14/14 |
| 5  | 2605.83 | 2600.47 | +2.1e-03 | 1.8867 | 1.8865 | +9.5e-05 | 0.3530 / 0.3538 | 44/45 |
| 10 | 2691.97 | 2683.57 | +3.1e-03 | 1.8893 | 1.8891 | +1.4e-04 | 0.3610 / 0.3612 | 50/50 |
| 11 | 2558.50 | 2703.25 | **−5.4e-02** | 1.8899 | 1.8897 | +1.4e-04 | 0.3476 / 0.3475 | **43/50** |
| 15 | 2576.87 | 2782.33 | −7.4e-02 | 2.0507 | 2.0474 | +1.6e-03 | 0.3615 / 0.3576 | 45/53 |
| 20 | 2602.02 | 2888.45 | −9.9e-02 | 2.4888 | 2.5199 | **−1.2e-02** | 0.3314 / 0.3198 | 45/57 |

Against #247's numbers (which stopped at day 19):

| quantity | #247 (bedrock mismatch) | D3 (corrected) | verdict |
|---|---|---|---|
| `h2ovol1`, days 1–20 | J **0.11–0.14** vs F 0.32–0.34 (3× dry) | J 0.324–0.400 vs F 0.320–0.422, tracking to **~0.003–0.01** | **CLOSED** |
| `maxdbh` day 19 | 1.9754 vs 2.5199 → **−22 %** | 2.3870 vs 2.4127 → **−1.1 %** | **CLOSED (20× better)** — the port now *enters* the day-13 vigorous-growth phase (1.89 → 2.49) instead of stalling |
| `btran_ft` days 11–13 | J collapses to 0.0, F stays 1.0 | no collapse | **CLOSED** |
| site carbon day 19 | −5.5 % | **−9.4 %** | **SURVIVES, and is a different mechanism** |

So D2's fix is confirmed end-to-end: the host over-drying, the `btran_ft` collapse and
the stature-growth stall were all one `use_bedrock` artefact, and with it gone the port
tracks Fortran's dbh trajectory to ~1 %. The carbon gap, however, did **not** close —
it got *worse*, which is the tell that it was never the same bug.

## The surviving carbon gap, attributed: a **dead running-mean wiring** in the port

Phase attribution (the #217 template) localised it in one step. Site carbon by phase
on the day it breaks (day 11 = step 241) versus the day before:

| phase | day 10 J / F | day 11 J / F |
|---|---|---|
| 10 `dyn_in` | 2673.26 / 2664.96 | 2691.97 / 2683.57 |
| 11 `phenology` | 2673.26 / 2664.96 | **2537.27** / 2683.57 |
| 13 `integrate` | 2692.05 / 2683.66 | 2558.59 / 2703.34 |
| 15 `cohortfuse` (ncoh) | 50 / 50 | **43** / 50 |

The entire loss appears **inside `phenology`**, where Fortran is a no-op. Cohort dumps
show why: on 11 July the port flips `status` 2 → 1 and zeroes `leafc` on every
**cold-deciduous** PFT (3, 6, 11) — a mid-summer leaf-off Fortran never performs.

### Root cause: `InitTimeAveragingGlobals()` was ported and never called

Instrumenting the site phenology state (`scripts/` probe, `vegtemp_memory` / `ncolddays`
/ `cstatus` per day) showed `tveg24` advancing **only every second day** — each daily
`phenology` call saw the same 24-hour mean twice, so `vegtemp_memory` filled at half
rate and its still-**uninitialised 0.0 slots counted as cold days** (`0.0 <
phen_coldtemp = 7.5`). On day 11 — the first day the `model_day > num_vegtemp_mem = 10`
gate opens — `ncolddays = 6 > ncolddayslim = 5` and the cold leaf-off fired:

```
        pre-fix           post-fix
day   tveg24  ncold     tveg24  ncold
  1   15.000    10      15.000    10
  3    3.993     8       5.152     9
  5    7.451     8       7.900     9
  7    6.692     8       7.700     8
  9   12.359     8      12.383     6
 11   11.887     6 -> LEAF OFF    13.388  4  (no leaf-off)
```
(the pre-fix `tveg24` column repeats each value for two consecutive days — the tell)

Fortran `FatesInterfaceMod.F90:1041-1047` defines every running-mean with
**`up_period = hlm_stepsize`**, the host's actual timestep, and derives
`n_mem = nint(mem_period/up_period)`. The port has all of this — `hlm_stepsize` and
`InitTimeAveragingGlobals()` in `FatesInterfaceMod.jl` — and **neither was ever
called or consumed anywhere in `src/`**. Patches instead used the `_patch_*`
fallbacks in `FatesPatchMod.jl`, which hardcode `up_period = 1800 s, n_mem = 48`.
On the reference case's `dtime = 3600 s` that is a **48-HOUR "24-hour" window**.

This is the `dead-initcold-systemic` class: ported, correct, never wired. The
source comment even advertised it ("not yet ported… can be re-pointed to the real
definitions"), and the harness ran for the port's whole life on the stub.

**Fix** (`fates(rmean)`): `clm_fates_init!` takes `dtime`, sets `hlm_stepsize[]` and
calls `InitTimeAveragingGlobals()`; `InitRunningMeans!` selects the defined global
via `_rmdef`, falling back to the `_patch_*` constant only when FATES init never ran
(unit-test fixtures), so no fixture changes behaviour.

### Result — the 20-day carbon gap closes

| day | relC before | relC **after** | relDBH before | relDBH **after** | ncoh J/F before | **after** |
|---|---|---|---|---|---|---|
| 11 | −5.4e-02 | **+3.0e-03** | +1.4e-04 | +1.4e-04 | 43/50 | **50/50** |
| 15 | −7.4e-02 | **+3.0e-03** | +1.6e-03 | +1.1e-02 | 45/53 | **53/53** |
| 20 | −9.9e-02 | **+2.4e-03** | −1.2e-02 | +7.1e-03 | 45/57 | **57/57** |

Site carbon tracks Fortran to **0.25 %** across the whole 20-day window (was −9.9 %),
the demography is **cohort-count-exact against Fortran every day**, and the
record-set mismatch rate drops from **311/481 to 96/481** at phase `fast` and
**12/20 to 4/20** at the daily phases. Cold start is unchanged at `7.133e-16` and the
first divergence is still the same host-driven `npp` `4.9218e-09` at step 0.

**What still survives:** `maxdbh` runs **+0.7 – 1.1 % high** from day 14 (it was
−1.2 % low before; the sign flipped, so this is a genuine residual and not the same
defect), and the `hmort`-family `NaN-only-one-side` diagnostics are unchanged and
pre-existing. Both are smaller than the two bugs already removed and are the next
D-tier item, not a blocker.

---

# D4 (2026-07-20): re-evaluating the #197 screen and the #227 moisture BC after D3

D3 (#274) closed the 20-day Bow window and, in doing so, cast doubt on two gated
features that had never been validated against Fortran:

* **#197** — `fates_biogeog_screen = :drop_cold_deciduous`, introduced to suppress a
  multi-year **boom-bust die-back** at tropical Aripuanã.
* **#227** — a harness-only growing-season root-zone **moisture prescription**,
  introduced to stop boreal Krycklan collapsing to zero carbon.

D3 suspected both of being *compensating workarounds for defects since fixed*. This
section re-measures both on current `main`.

## 1. #197 — D3's suspicion is WRONG: the running-mean bug was inert at `dtime = 1800`

D3's hypothesis was that the die-back the screen suppresses was "partly the
running-mean bug" — the 48-hour `tveg24` window firing a spurious cold leaf-off on
exactly the cold-deciduous PFTs the screen drops. **That cannot be true, and the
reason is structural, not statistical.**

The bug's magnitude is a function of the host timestep. `_patch_fixed_24hr` hardcoded
`up_period = 1800 s, n_mem = 48`; `InitTimeAveragingGlobals()` derives
`n_mem = nint(86400/dt)` at `up_period = dt`. Comparing the two definitions directly:

| `dtime` | global `fixed_24hr` | `_patch_fixed_24hr` fallback | fallback window in REAL time | identical? |
|---|---|---|---|---|
| **1800 s** | `up=1800, n_mem=48` | `up=1800, n_mem=48` | **24 h** | **YES — bit-identical** |
| 3600 s | `up=3600, n_mem=24` | `up=1800, n_mem=48` | **48 h** | no |

**Every number in #197 was produced by `scripts/fates_longhorizon.jl`, which runs at
`dtime = 1800`** — the one timestep at which the pre-#274 fallback was *exactly
correct*. The cold-leaf-off decision path reads only `tveg24` (`EDPhysiologyMod.jl:541`
→ `vegtemp_memory` → `ncolddays`), and `tveg24` is defined by `fixed_24hr`. So the
entire phenology branch D3 fixed is **byte-identical before and after #274 at
`dtime = 1800`**. The 48-hour window was a property of the *Fortran reference case's*
`dtime = 3600`, not of #197's evidence base.

Corollary: the die-back mechanism #197 actually describes is a different branch
entirely — the **400-day `nevercold` forced leaf reset** at `EDPhysiologyMod.jl:637`
(`cstatus == phen_cstat_notcold && cndaysleafoff > 400`). That is a pure day counter.
It cannot fire before model day 400, it is untouched by #274, and it is the reason the
artifact is a *multi-year* one.

**The screen is therefore NOT compensating for the defect #274 fixed.** It may still be
compensating for something else, but D3's specific reasoning for retiring it does not
hold and should not be acted on.

## 2. #197 — but its evidence no longer reproduces AT ALL, in either arm

Re-running #197's own recipe on current `main` (Aripuanã, `dtime = 1800`, 1460 days):

| run | `fates_biogeog_screen` | result |
|---|---|---|
| baseline | `:none` | **dies at day 74** — `errlon = -412.97 W/m2` at `p=7`, FATAL |
| screened | `:drop_cold_deciduous` | **dies at day 98** — `errlon = -405.41 W/m2` at `p=7`, FATAL |

With balance hard-errors degraded to warnings (`FATES_SOFT_BALANCE=1`, added here,
default off) the run gets one day further and then hits a hard FATES `endrun`:

```
✗ ERROR at step 3621 (day 75): ENDRUN: imaginary roots detected in
                               QuadraticRootsSridharachary   (FatesUtilsMod.jl:201)
```

So the energy-balance blow-up corrupts the canopy state and FATES aborts inside
photosynthesis. **Neither arm reaches model day 100, let alone the day-400 `nevercold`
timer or the 4-year horizon #197's claim rests on.** The cold-start carbon still
reproduces #197's number exactly (`2549`), so the setup is right — the model no longer
survives the horizon.

**Verdict on #197: NEITHER retire NOR validate it yet.** D3's stated ground for
retirement is disproved (§1), and the measurement that would settle the question is
blocked by §4. It is gated and default-off; leave it, with this section as the record
of *why* it is unresolved rather than unexamined.

## 3. #227 — this one IS compensating for something already fixed

**RETIRED (#278, 2026-07-20).** D4's recommendation was carried out, but only after
re-measuring it on a **controlled A/B** that D4 did not run — and that A/B both
strengthens the verdict and **corrects one of D4's two arguments for it.**

### 3a. D4's comparison was confounded; the clean A/B is same-site ± prime

D4 compared `krycklan_baseline` (prime OFF) against `krycklan_screened` (prime ON).
Those two sites **differ in the PFT set as well as the prime** — baseline is all-14-PFT,
screened is the boreal set `{2,3,6}` — so the comparison cannot isolate the prescription.
`FATES_SOIL_H2O=off` was added to disable the prime on a site that carries one, giving a
true one-variable A/B. All rows below are **the same commit, the same site, the same PFT
set, 365 days**, differing only in the prime:

| `krycklan_screened` | cold | final | min | peak | 2nd-half band | maxdbh | births/deaths | checks |
|---|---|---|---|---|---|---|---|---|
| prime **ON** (#227 default) | `749.6` | `1986` | `511.4` | `2569` | `[1194, 2569]` | `2.225 cm` | 44 / 13 | **10/10 PASS** |
| prime **OFF** (`FATES_SOIL_H2O=off`) | `749.6` | `1737` | `512.6` | `2341` | `[1169, 2341]` | `2.135 cm` | 40 / 11 | **10/10 PASS** |

Both arms ran **364/365 days**, carbon conserved 364/364.

### 3b. The collapse is gone — confirmed independently, and on a *completing* run

**#227's premise no longer reproduces.** Its stated justification was a monotone
`749.6 → 0` with `btran = 0` every growing-season day, on exactly this configuration.
With the prime off, that configuration now runs a full year and passes every check,
carbon `749.6 → 1737`. This is stronger evidence than D4 had: D4's OFF arm died at day
235 on the §4 `errlon`, so it could only argue "not a monotone collapse from day 1",
whereas the screened OFF arm **completes the year**.

D4's `krycklan_baseline` measurement also reproduces exactly — cold `2549` → `4040`,
min `2008`, no-collapse and no-boom-bust both PASS, dying at **day 235** on the same
fatal `errlon = -351 W/m2` at `p=7` (§4).

### 3c. Correction to D4: the cold start did NOT move

D4's second argument was that "the cold start itself is different — `749.6` then vs
`2549` now", inferring #227 had been measured on a materially different configuration.
**That inference is wrong, and it is an artifact of the confound in 3a.**
`krycklan_screened` still cold-starts at **exactly `749.6`**, identical to #227's report;
the `2549` is simply the all-14-PFT `krycklan_baseline` cold start. The two numbers are
two different PFT sets, not two different eras. There is no evidence here of cold-start
drift, and the `use_bedrock` (#251/#252) attribution built on it is **unsupported** —
it may still be the true cause of the collapse's disappearance, but this measurement
does not show it.

### 3d. It was NOT inert — it was actively distorting a now-correct model

The question that decides *retire* vs *leave it*: with the underlying defect fixed, does
enabling the prescription still change anything? **Yes.** Same site, same PFT set:

| quantity | prime ON | prime OFF | effect of the prime |
|---|---|---|---|
| final carbon | `1986` | `1737` | **+14.3 %** |
| peak carbon | `2569` | `2341` | +9.7 % |
| final maxdbh | `2.225 cm` | `2.135 cm` | +4.2 % |
| elai max | `1.102` | `1.084` | +1.7 % |
| births / deaths | 44 / 13 | 40 / 11 | +10 % / +18 % |

So it is **dead weight in the worst way**: not harmless dead code, but a live water
injection inflating boreal carbon by ~14 % on a model that no longer needs correcting.
Note that **carbon conservation held 364/364 in both arms** — the balance check is
completely blind to this, exactly the `conservation-is-not-accuracy` class. A green
balance was never evidence for #227 and is not evidence against it.

### 3e. What was done

- `SiteCfg.soil_h2o` is now `nothing` for **every** site — `krycklan_screened`'s `0.75`
  removed. No site is primed by default.
- `prime_fates_soil_moisture!` is **demoted, not deleted** (docs still reference it, and
  it has a legitimate diagnostic use: prime the soil, and if a symptom survives, water is
  not the lever). It is opt-in via `FATES_SOIL_H2O=<frac>` — `0.75` restores #227's
  setting — and prints a multi-line `!! DIAGNOSTIC ... NOT a validated config` banner
  whenever it activates, quoting the +14 % distortion.
- **Default path confirmed unaffected**: re-running `krycklan_screened` post-retirement
  with no env override reproduces the explicit-`off` arm on every reported statistic
  (`749.6 → 1737`, band `[1169, 2341]`, maxdbh `2.135`, 40/11, 10/10 PASS).
- **No regression on unprimed sites**: `hubbardbrook_screened` `1784 → 2895`, band
  `[2049, 3159]`, 364/364 balance, 10/10 PASS.
- No test fixtures encoded the boreal numbers (`grep` over `test/` finds no Krycklan
  FATES expectation), so re-baselining is harness + docs only.

### 3f. Re-baselined boreal band

The boreal reference for `krycklan_screened` (365 d, un-primed) is now:

```
cold 749.6 → final 1737   min 512.6   peak 2341   2nd-half band [1169, 2341] (hi/lo 2.00)
ncoh [3, 33]   npatch max 2   births 40 / deaths 11   elai [0.100, 1.084]   maxdbh 2.135 cm
364/365 days, carbon conserved 364/364, 10/10 checks PASS
```

The configured `[200, 6000]` biome band still holds and is left unchanged. The Fortran
oracle D3 asked for is still not the first need — the artifact it was meant to measure
has stopped reproducing.

**Scope limit:** this is a 365-day result. Multi-year FATES remains broken on `main`
(§4; root-caused in PR #280 — filter masks built once from surfdata and never rebuilt
after FATES daily dynamics), so #227's 1460-day claims are not re-measurable and are not
re-measured here. Pre-fix multi-year boreal numbers are additionally suspect because
`fates_longhorizon.jl` used to freeze on stale Dec-31 forcing after year one.

## 4. The actual blocker: a fatal longwave balance error kills EVERY multi-year FATES run

Not site-specific, not config-specific, and it is what stopped both investigations above:

| site | config | dies at | `errlon` | patch |
|---|---|---|---|---|
| Aripuanã (tropical) | screen `:none` | day 74 | `-412.97 W/m2` | 7 |
| Aripuanã (tropical) | screen `:drop_cold_deciduous` | day 98 | `-405.41 W/m2` | 7 |
| Krycklan (boreal) | baseline, no moisture prime | day 235 | `-351.36 W/m2` | 7 |
| Bow (the Fortran-validated site) | `FATES_WARMSOIL=0` | day 147 | `-294.76 W/m2` | 4 |

Three sites, four configurations, errors of 300-400 W/m² on a longwave flux. The check
itself is original (`ff6543c`, the first balance-check port), so what changed is the
model, not the assertion.

**It is a REGRESSION, and it is bracketed.** The identical Aripuanã baseline run at
#197's own merge commit (`712cb5c`) completes cleanly:

```
ran 149/150 days …  carbon: cold=2549 final=3822 …  daily balance held: 149/149
★ FATES stable + demographically active over 150 days      (8/8 PASS, no errlon)
```

versus current `main`, which dies at **day 74**. Bisected over `712cb5c..37db702`
(90-day Aripuanã baseline, 6 points):

| commit | | 90-day run |
|---|---|---|
| `d7fb711` | #251 (docs only) | **89/90 PASS** |
| `f2d805c` | **#252 — `use_bedrock` conditional default** | **dies day 74, `errlon` `p=7`** |
| `6405ca1` | #225 | dies day 74 |
| `bac5908` | #255 | dies day 74 |

**The culprit is #252**, and the failure appears at *exactly* the same step (3601) from
that commit onward. #252 made `use_bedrock`'s default conditional (`use_fates ⇒ false`);
it had been `true`. `fates_fortran_parity.jl` passes the flag **explicitly**, which is
why D3's flag audit came back clean — but `fates_longhorizon.jl` and
`fates_multisite_validation.jl` **inherit** it, so #252 silently flipped both harnesses
from `use_bedrock = true` to `false`.

### Pinning the flag: the mechanism holds at #252 — and exposes a SECOND defect

Re-running with `use_bedrock` pinned (`FATES_USE_BEDROCK`, added here) separates the two:

| | `use_bedrock` inherited (`false`) | `use_bedrock = true` (pre-#252 value) |
|---|---|---|
| at `f2d805c` (#252) | dies day 74 | **89/90 PASS** — `ncoh` max 69, `npatch` 2 |
| current `main` | dies day 74 | **still dies day 75** — `ncoh` max 363, `npatch` 5 |

So:

1. **At #252 the flip IS the mechanism** — restoring `use_bedrock = true` there fully
   rescues the run. The defect it exposes is *older* than #252; #252's only `src/`
   changes are the kwarg becoming `Union{Bool,Nothing}` in `clm_initialize.jl` /
   `clm_run.jl`, no physics.
2. **On current `main` that rescue no longer works**, so a **second, independent cause**
   entered in `f2d805c..37db702`. The demography confirms it is a distinct fault: with
   the flag pinned identically, `main` reaches `ncoh = 363 / npatch = 5` where #252
   reaches `69 / 2` (and #197-era holds `122 / 2` out at day 149) — recruitment has run
   away by roughly 5×.

**#252 is correct and must not be reverted** (it matches CTSM
`namelist_defaults_ctsm.xml:178-181`). There are now **two** defects to fix, and the
first has been latent since long before either.

Likewise, #227's collapse disappearing is a **measured fact but an unattributed one**.
It is tempting to credit the same flip; given that the same inference just failed on
`main`, it is recorded here as correlation with #252's window, not as mechanism. Worth noting the demography also diverged sharply over that
range — at day ~144 the #197-era run holds `ncoh = 118, npatch = 2`, while current
`main` had already reached `ncoh = 446, npatch = 6` by day 74. Note the daily FATES carbon balance holds *right up to the
failure* in every case (`235/235`, `147/147`, `75/75` days) — a textbook instance of
the `conservation-is-not-accuracy` class: the carbon books balance perfectly while the
energy books are off by 400 W/m².

**This is now the top D-tier item**, ahead of both #197 and #227: D1's long-standing
"still open: a multi-year FATES run" is not merely unstarted, it is *impossible* on
`main` today. Neither gated feature can be validated until it is fixed.

## Reproducing

```bash
SD="…/My Drive/data/SYMFLUENCE_data"   # NOT the compHydro root — see below

# #197 both arms (Aripuana, dtime=1800)
SYMFLUENCE_DATA=$SD FATES_NDAYS=1460                  julia +1.12 --project=. scripts/fates_longhorizon.jl
SYMFLUENCE_DATA=$SD FATES_NDAYS=1460 FATES_BIOGEOG=1  julia +1.12 --project=. scripts/fates_longhorizon.jl
# add FATES_SOFT_BALANCE=1 to run past the fatal errlon (diagnostic only —
# absolute carbon is NOT trustworthy with a known-open energy balance)

# #227 — the CONTROLLED A/B (#278): same site, same PFT set, prime the only variable.
# (The krycklan_baseline-vs-krycklan_screened pair D4 used is confounded — it also
#  changes the PFT set. See §3a.)
SYMFLUENCE_DATA=$SD SITE=krycklan_screened FATES_NDAYS=365 \
    FATES_SOIL_H2O=0.75 julia +1.12 --project=. scripts/fates_multisite_validation.jl  # prime ON (#227)
SYMFLUENCE_DATA=$SD SITE=krycklan_screened FATES_NDAYS=365 \
    julia +1.12 --project=. scripts/fates_multisite_validation.jl                      # prime OFF (default now)

# D4's original (confounded) pair, kept because §3b/§4 quote its numbers
SYMFLUENCE_DATA=$SD SITE=krycklan_baseline FATES_NDAYS=365 julia +1.12 --project=. scripts/fates_multisite_validation.jl
```

**Harness trap found while doing this:** `fates_multisite_validation.jl` **hardcoded**
the data root, so it could only ever see `/Users/…/compHydro/SYMFLUENCE_data` — which
holds none of the three FATES domains. They live under the Drive root. This is why #227
reported "Tropical Aripuana data is absent in this environment" and skipped the tropical
arm: the domain existed all along. Now honours `SYMFLUENCE_DATA` like every other
harness (`test/testdata.jl::symfluence_data_root`).

---

# SOIL-BGC COUPLING (2026-10-28): FATES litter -> CLM decomposition, validated

## What was missing

CTSM runs every non-SP FATES case with `use_fates_bgc = .not. use_fates_sp`
(`controlMod.F90:474-479`) and `soil_decomp_method = 'CENTURYKoven2013'`. The archived
`fates_parity_1pt` case (and therefore every reference dump in this directory) has the
coupling on. The port did not. `CLMDriverConfig(use_fates=true)` defaulted
`use_fates_bgc=false`. FATES computed `bc_out%litt_flux_{lab,cel,lig}_c_si`, but nothing
outside `src/fates/` read them. No decomposition cascade was built for FATES columns, and
soil HR (and so `FATES_HET_RESP`) was identically 0. The scorecards above covered FATES
*vegetation* only.

## What is connected now (Fortran reference -> port)

| Fortran (ctsm5.3.012) | port |
|---|---|
| `controlMod` `use_fates_bgc = .not. use_fates_sp`; suplnitro must be `'ALL'` for `fates_parteh_mode=1` | `CLMDriverConfig` / `clm_initialize!(use_fates=true)` default `use_fates_bgc`; `soil_bgc_competition_init!(suplnitro=SUPLN_ALL)`, and `carbon_only` is now passed to `soil_bgc_competition!` (it never was, so fpi = 0 blocked every litter transition) |
| `init_decomp_cascade_constants`: 6 pools / 8 transitions under FATES-BGC (no CWD) | `CLMDriverConfig.ndecomp_pools/transitions`; `init_decomp_cascade_bgc!(use_fates=true)`; the unused CWD slot of the CN-sized arrays is held at 0 |
| `clm_varpar`: `nlevdecomp = nlevsoi`; `use_vertsoilc=1`, `decomp_method='CENTURY'` passed to FATES | `clm_fates_init!(use_vertsoilc=true, decomp_method="CENTURY")`, `bc_in.nlevdecomp = nlevsoi` |
| `init_soil_depths` (real `zi/dz/z`, `dz_decomp_sisl`, `decomp_id(j)=j`) | `fates_init_soil_depths!` (replaces the synthetic 0.1 m layering under FATES-BGC) |
| `CStateUpdate1` -> `UpdateCLitterFluxes` -> `FluxIntoLitterPools` -> `decomp_cpools_sourcesink` (+ `fates_litter_flux`) | `fates_update_c_litter_fluxes!`, called from `cn_driver_no_leaching!` just before `c_state_update1!` |
| `NStateUpdate1` -> `UpdateNLitterFluxes` (zeroes N source/sink; N transfer commented out upstream) | `fates_update_n_litter_fluxes!` |
| `wrap_update_hifrq_hist` (after EcosystemDynamicsPostDrainage: `tot_het_resp/tot_somc/tot_litc`) | `fates_wrap_update_hifrq_hist!`; the hifrq fill moved to the same point of the step |
| `wrap_update_hlmfates_dyn` `litr_lig_c_to_n_col <- litt_flux_ligc_per_n` | `fates_update_litr_lig_c_to_n!` in the daily step |
| `CNDriverSummarizeStates` soil C summary at end of step, also for FATES-BGC | the driver's end-of-step soil summary is now gated `use_cn || use_fates_bgc` (FATES columns carry no CN veg C). Without it the reported pools lagged by one step |

The dynamics-driv `w_scalar/t_scalar` pack (already ported, gated on `use_fates_bgc`) now
receives live decomposition scalars, so FATES litter fragmentation runs on the host's
soil temperature and moisture, as in Fortran. Default non-FATES paths are unchanged.

## Fortran reference and result

`stage_soilbgc_run.sh` re-runs the archived case with the **uninstrumented**
`installs/clm/bin/cesm.exe` (FATES is always compiled into CTSM). The setup matches
`setup_case.sh` (Bow single point, 2002-07-01, 3600 s, `use_lch4` off), with an hourly
h0 history of HR, pools, scalars and FATES litter fields. Two independent launches
gave byte-identical output. The scored series are committed as
`fates_soilbgc_fortran.csv` (240 steps = 10 days):

```bash
julia +1.12 --project=. scripts/validation/fates_fortran_parity/fates_soilbgc_parity.jl
```

| field | end (day 10) Fortran / Julia | per-step rel. error from day 2 (max / mean) |
|---|---|---|
| HR | 3.8352e-06 / 3.8400e-06 gC/m2/s | 5.2e-2 / 1.9e-2 |
| LITTERC_HR | 1.5679e-09 / 1.5749e-09 | 8.6e-2 / 1.6e-2 |
| TOTLITC | 1.6561e-02 / 1.6583e-02 gC/m2 | 2.2e-2 / 6.4e-3 |
| LIT_MET/CEL/LIG_C | within 2.7e-2 | |
| TOTSOMC | 910.5446 / 910.5238 gC/m2 | 2.3e-5 / 1.2e-5 |
| FATES_LITTER_OUT (FATES litter flux into the soil) | 4.4043e-11 / 4.3978e-11 kg/m2/s | 2.7e-2 / 7.1e-3 |
| FATES_HET_RESP | = HR x 1e-3, both sides | |

Daily HR totals agree within 3.5 % (day 2), tightening to 0.3 % by day 10. Litter enters
on the same step in both models (nstep 26, the step after the first daily fragmentation).
**Day 1 is the one open residual.** During the cold-start thaw (nsteps 2-8) Julia's HR is
0-85 % low, because its top-layer `W_SCALAR` on the 273.15 K melt plateau is 0.02-0.14
against Fortran's 0.09-0.23. That is the host phase-change residual already documented
above, not the coupling. `TSOI_1..3` and `T_SCALAR_1..3` track to < 1 % mean.

### Effect on the FATES dynamics scorecard (120-step oracle, same tree +/- this change)

The FATES vegetation now sees the same soil layering, decomposition layers and
fragmentation scalars as the Fortran reference, and the oracle moves toward Fortran:

| | before | after |
|---|---|---|
| site carbon day 1 / day 5 (Fortran 2547.884 / 2600.471) | 2548.095 / 2605.611 (8.3e-5 / 2.0e-3) | 2547.931 / 2605.458 (1.8e-5 / 1.9e-3) |
| worst field, phases 13-14 (`dndt`) | 2.069e+02 | 2.475e+00 |
| cohorts day 4 (Fortran 42) | 43 | 43 |

### Still not validated
* FATES N/P litter coupling: commented out upstream in ctsm5.3 too, so N litter input is 0
  on both sides and `suplnitro='ALL'` supplements immobilization.
* MIMICS under FATES (`clm_initialize!` refuses `decomp_method=2` with FATES-BGC).
* FATES-BGC restarts, and windows longer than 10 days for the soil pools.

---

# FATES CNP (`fates_parteh_mode = 2`, 2026-10-28)

## What the ctsm5.3.012 host actually couples

CNP PARTEH is FATES-internal; the CTSM host connects almost nothing to it, and the port
now mirrors that rather than inventing a coupling:

| Fortran (ctsm5.3.012) | port |
|---|---|
| `bc_in%plant_nh4/no3/p_uptake_flux` are allocated but **never filled** by any host file (no FATES branch in `SoilBiogeochemCompetitionMod`; `plant_ndemand_col = 0` for FATES columns, `CNDriverMod.F90:466-468`) | unchanged: no host uptake. Soil mineral N is not drawn down by FATES plants on either side |
| `check_allset` endruns CLM CNP unless some `fates_cnp_prescribed_n/puptake > 0` (`FatesInterfaceMod.F90:1728-1736`); testmod `FatesColdPRT2` sets both to 1.0 | `clm_fates_init!`: a parameter file that prescribes uptake is used as-is; the shipped default (both 0) gets the FatesColdPRT2 values instead of aborting |
| `SetFatesGlobalElements2` uptake modes / `max_comp_per_site` / `fates_np_comp_scaling` (`:901-926`) | derived from the parameter values (previously hard-coded prescribed/trivial) |
| `UpdateNLitterFluxes`: N litter transfer commented out; `fates_litter_flux = 0` | unchanged (`fates_update_n_litter_fluxes!`, #336) |
| `wrap_WoodProducts` N products: commented out | nothing to port |
| `nitrogen_spec` = 1 with nitrif/denitrif else 2, `phosphorus_spec` = 0 (`clmfates_interfaceMod.F90:453-460`) | `clm_initialize!(fates_parteh_mode=2)` -> `_fates_parteh_kwargs` |
| `FastPRTAllometricCNP` | a no-op stub in Fortran too (`PRTAllometricCNPMod.F90:2005-2015`) |
| `EDMainMod.F90:707-709` `fates_hist%update_history_nutrflux` | now called (history passed down as `ed_ecosystem_dynamics(...; fates_hist)`). Its never-run port had a level-1 `UndefVarError` (scpf handles only bound at level 2), fixed |
| `update_history_dyn2` N/P scpf pools and storage fractions (`FatesHistoryInterfaceMod.F90:3338-3341, 3615-3650, 4627-4865`) | filled (`FATES_{VEG,LEAF,FROOT,SAPWOOD,STORE,REPRO}{N,P}_SZPF`, `FATES_STORE{N,P}_TF_{CANOPY,USTORY}_SZPF`) |

Carbon-only runs never reach any of this (all gated on `hlm_parteh_mode == 2` or
`fates_parteh_mode == 2`).

## Fortran reference and result

`stage_cnp_run.sh` (FatesColdPRT2 settings, 720 hourly steps = 30 days;
`PARTEH=1` stages the same case carbon-only as a control). The scored series are
committed as `fates_cnp_fortran.csv`:

```bash
julia +1.12 --project=. scripts/validation/fates_fortran_parity/fates_cnp_parity.jl
FATES_CNP_PARTEH=1 julia ... fates_cnp_parity.jl <carbon-only h0.nc> 170   # control
```

Days 1-10 (steps 1-240), per-step relative error max / mean:

| field | max | mean |
|---|---|---|
| FATES_VEGN / VEGP, STOREN / STOREP, LEAFN / LEAFP, FROOTN / FROOTP | 4.3e-3 | 2.9e-3 |
| FATES_SAPWOODN / P | 9.4e-3 | 6.4e-3 |
| FATES_STOREN_TF / STOREP_TF | 3.2e-6 | 1.5e-6 |
| FATES_NH4UPTAKE, NO3UPTAKE, NDEMAND, PDEMAND | 4.0e-3 | 2.4e-3 |
| FATES_PUPTAKE (net P allocated incl. storage draw-down, a difference of pools) | 1.7e-1 | 4.7e-2 |
| FATES_VEGC / LEAFC / STOREC / FROOTC | 2.9e-3 / 4.0e-3 / 4.4e-3 / 4.0e-3 | 2.0e-3 / 2.8e-3 / 3.0e-3 / 2.8e-3 |
| FATES_L2FR | 1.0e-6 | 5.4e-7 |
| FATES_LAI | 5.4e-3 | 3.7e-3 |
| FATES_NCOHORTS | 50 / 50 at day 10 | |
| SMINN (mean), TOTSOMC | 6.4e-2, 2.2e-5 | |

That is the same level as the carbon-only drift on this case, and the CNP-specific
increment matches: at day 4 (CNP minus carbon-only) VEGC is +2.38e-7 kg/m2 in both
models, LEAFC +4.0e-8 vs +3.9e-8, STOREC +5.4e-8 vs +5.6e-8, LAI +9.9e-7 vs +1.0e-6.

**Day 11 onward is a threshold bifurcation, not a CNP divergence.** Fortran CNP drops the
leaves of the five cold-deciduous PFTs on day 11 (COLD_STATUS 2 -> 1, LAI 0.66 -> 0.23).
Fortran carbon-only does not, and neither does the port in either mode. The cause is
`EDPhysiologyMod.F90:1150-1153`: leaf-off fires once `model_day > 10` if more than
`phen_ncolddayslim` (5) of the last 10 daily TVEG24 values are below `phen_coldtemp`
(7.5 C). Day 6 is the knife edge. Fortran carbon-only has 7.514 C, Fortran CNP has
7.369 C, and the port has 7.700 C in both modes. (Update: the port is now 7.512 C. The CNP/carbon-only split is a Fortran
round-off flip of the trace-snow `z0mg` branch; see "Cold canopy, canopy snow and the
winter warm bias".) The two Fortran runs agree to about
1e-7 until nstep 139, then jump together in GPP *and soil HR* (a host-side event that
amplifies the round-off-level CNP perturbation). The port is warmer than both by the
known ~0.2-0.8 C host vegetation-temperature residual, so it stays on the warm side.
Over 30 days that gives LAI 0.86 (port) vs 0.28 (Fortran CNP) vs 0.86 (Fortran
carbon-only).

Not scored: `SMIN_NH4` / `SMIN_NO3`. The port's soil-N summary never fills the column
NH4/NO3 totals (`use_nitrif_denitrif` is not passed to `soil_bgc_nitrogen_state_summary!`,
which is a host diagnostic gap outside CNP). The vertically resolved pools and `SMINN` are
fine.

---

# REDUCED-COMPLEXITY MODES (2026-10-28): SP, nocomp, fixed biogeography

`clm_initialize!(; use_fates=true, use_fates_sp, use_fates_nocomp,
use_fates_fixed_biogeog)` now switches on the three CTSM `&clm_inparm` FATES modes
(previously pinned off in `_fates_spike_set_ctrlparms!`). Pass the same `use_fates_sp` to
`CLMDriverConfig`; `clm_drv!` errors if the two disagree.

## What is connected (Fortran -> port)

| CTSM / FATES | port |
|---|---|
| controlMod: `use_fates_bgc = .not. use_fates_sp`; FATES `check_allset`: SP needs nocomp + fixed biogeog | `clm_initialize!` validation; `use_fates_bgc` default |
| `CLMFatesGlobals1` `use_nocomp` / `use_sp` | `_fates_spike_set_ctrlparms!` |
| `SetFatesGlobalElements1` patch sizing (nocomp: primary land >= numpft; SP: `max(surf_numpft+surf_numcft, numpft+1)`) | `clm_fates_init!`, `fates_max_patches_per_site` (HLM slot reservation) |
| `clm_varpar_init` FATES branch: `surfpft_ub = surf_numpft+surf_numcft-1`, `maxveg = maxsoil_patches-1` | `clm_initialize!` (mode paths only) |
| clmfates_interfaceMod init: `pft_areafrac = wt_nat_patch(g, 0:surfpft_ub)`, sum-to-1 check; `fates_hlm_pft_map` from the param file | `hlm_pft_areafrac` -> `bc_in.pft_areafrac` / `baregroundfrac` |
| initGridCellsMod: SP builds every natpft slot with `ptype = m` | `set_landunit_veg_compete!` SP branch |
| clm_initializeMod: `interpMonthlyVeg` + `SatellitePhenology(inactive_and_active%soilp)` before `init_coldstart`; clm_driver.F90:1052 the same every `doalb` step | `_fates_sp_host_phenology!`, driver SP block |
| `hlm_sp_tlai/tsai/htop(ft) = tlai/tsai/htop(patchi+ft)` (htop floor 0.01) in `init_coldstart` and `dynamics_driv` | `fates_pack_bcin_sp!` |
| `wrap_update_hlmfates_dyn` SP: FATES totals -> `tlai/tsai/htop_hist_patch` | `fates_flush_hlm_canopy!`, `fates_unpack_bcout_canopy_structure!` |
| history `FATES_NOCOMP_PATCHAREA_PF` / `_NPATCHES_PF` / `_BURNEDAREA_PF` | `update_history_dyn2!` |

Two defects surfaced on the way and are fixed here: (1) the non-methane
`decomp_rate_constants_bgc!` call never passed `use_fates`, so under FATES soil BGC the
CWD->litter kernel stored at cascade index -1 inside `@inbounds`, corrupting the heap (it
showed up as intermittent GC/inference segfaults; `--check-bounds=yes` finds it); (2) the
dyn2 history handle fetch indexed `hvars[0]` for variables the registry skips in SP mode.
The FATES clock fix of PR #341 (`fates_get_and_set_time!`) is cherry-picked; the numbers
below include it.

## Fortran reference and result

`stage_fates_modes_run.sh` (`MODE=sp|nocomp|fixedbg`) runs the archived case with the
mode's namelist overrides (SP also needs `soil_decomp_method='None'`) and an hourly h0;
the scored series are committed as `fates_<mode>_fortran.csv`. Score with
`fates_modes_parity.jl <mode> [ref] [nsteps]`. 240 steps (10 days from 2002-07-01),
relative error of the window mean (|J-F|/|F|):

| field | SP | nocomp + fixed | fixed only |
|---|---|---|---|
| patches / cohorts | 5/4 every step, = F | 6/8 vs F 6/7 at day 10 | 2/11 every step, = F |
| patch area by PFT | 5e-9 - 2e-8 | 5e-9 - 2e-8 | (1 patch) |
| NPLANT by PFT | 2e-9 - 2e-8 | <= 4e-7 | <= 4e-6 |
| FATES_VEGC | 7e-9 | 3.1e-3 | 3.4e-3 |
| FATES_LAI | 2e-8 | 4.5e-3 | 4.5e-3 |
| TLAI | 1.2e-3 | 8.2e-3 | 5.0e-3 |
| FATES_GPP | +9.0 % | +3.2 % | +3.7 % |
| FSH | +2.0 % | +3.2 % | +6.0 % |
| EFLX_LH_TOT | +2.9 % | -0.8 % | -0.5 % |

* **SP** reproduces the FATES state to round-off: with prescribed LAI there is no
  allocation feedback, so it isolates the host. The GPP and flux gaps are the same host
  residual the default mode shows (top-layer thaw on day 1, see the soil-BGC section).
* **nocomp**: every patch has the right PFT and area; the one-cohort difference at day 10
  is a recruitment threshold crossed one day apart under the ~3 % GPP offset.
* **fixed only**: identical census for all 10 days; only the 4 surfdata PFTs that survive
  `min_nocomp_pftfrac_perlanduse` are seeded, as in Fortran.

### Still not validated
SP/nocomp restarts, LUH-driven fixed biogeography (`use_fates_luh`), dry deposition with
FATES-SP (CTSM disables it too), and runs longer than 10 days.

## The host bias behind the SP GPP gap (2026-10-29): four HLM<->FATES coupling defects

The table above was taken before this fix. SP mode was the test bed: prescribed LAI, FATES
state matched to ~1e-8, 10 days. A second SP Fortran run with `hist_dov2xy = .false.`
(per-patch hourly h0: `TV SABV FSA FSR FSRVD/ND/VI/NI FCEV FCTR FSH_V USTAR DISPLA SNOCAN
LIQCAN QINTR QDRIP ...`) was compared against the port patch by patch. Each diverging
quantity in the chain, in the order found:

1. **Radiation: visible reflectance 0.08 vs 0.20** (FSRVD/FSDSVD). NIR matched, so the
   leaf/stem optics were fine. `fates_pack_bcin_radiation!` hard-coded
   `bc_in%fcansno_pa = 0`, but Fortran `wrap_canopy_radiation`
   (clmfates_interfaceMod.F90:2798) passes the host `fcansno(p)` (CanopyHydrologyMod
   `BulkDiag_FracWet`), and `FatesNormanRadMod` (F90:207-212) blends `rho/tau` toward
   the snow optics with it. Early July 2002 at Bow has snow on the canopy, so the port's
   canopy absorbed about 7 % too much (SABV) and GPP was 9 % high. Fix: pass
   `waterdiagnosticbulk.fcansno_patch[p]`.
2. **Turbulence: DISPLA 1.10 vs 2.62 m.** `set_z0m_displa!` overwrote the FATES
   `z0m/displa` of FATES patches every step with the host `z0mr/displar * htop(itype)`.
   Fortran calls `clm_fates%TransferZ0mDisp` and then sets only `.not. patch%is_fates`
   patches (BiogeophysPreFluxCalcsMod.F90:162-173; `patch%is_fates = patchi+1:patchf`,
   clmfates_interfaceMod.F90:934). Fix: `fates_transfer_z0m_disp!` plus a
   `patch_is_fates` mask passed from `clm_drv!`. Per-patch DISPLA and USTAR now match.
   (The gridcell `USTAR` still differs because Fortran leaves the bare patch's `fv` at
   spval, so its average skips that patch.)
3. **Canopy water: LIQCAN 0.04 vs 0.42 mm** on the first rain step, with QINTR/RAIN at
   23 % vs 100 %. `use_clm5_fpi` was set only inside `clm_run!`. Every harness that
   drives `clm_initialize!` + `clm_drv!` directly (all the FATES parity scripts) ran the
   CLM4 interception form `0.25*(1-exp(-0.5*lai))` (CanopyHydrologyMod.F90:536-538).
   That is CTSM's code fallback; the clm5 namelist default is `.true.`. Fix:
   `clm_initialize!(; use_clm5_fpi = true)`. This was behind FCEV -45 % and FSH_V +15 %
   on the tree patches.
4. **Cold start: the HLM ran on the surfdata patch layout for the first step(s).**
   `col.is_fates` was only set by the scripts after `clm_initialize!` returned, so the
   cold-start `wrap_update_hlmfates_dyn` transfer (FATES patch weights, elai/esai,
   z0m/displa) never ran. Fortran `clmfates_interfaceMod::init` (F90:814) flags every
   istsoil column before `init_coldstart`, and `clm_initialize!` now does the same.

5. **Daily-step re-weighting happened mid-step.** `fates_set_filters!` wrote the new
   FATES patch weights into `wtcol` during the daily step. Fortran
   `wrap_update_hlmfates_dyn` writes only `patch%wt_ed`; `wtcol = wt_ed` is applied by
   `dynSubgrid_driver -> dyn_ED` (dynEDMod.F90:30-38) at the top of the next step,
   before `BeginWaterColumnBalance`. Re-weighting mid-step moved per-patch canopy water
   and plant storage between `begwb` and `endwb`. Once (3) put real water on the
   canopy, that showed up as a 1.08e-4 mm errh2o at the step-25 daily step in
   `test_fates_planthydro_coupling.jl`. `fates_apply_deferred_activation!` now performs
   `dyn_ED`. The SP scores are unchanged.

6. **Gridcell weights were never recomputed after `dyn_ED`.** CTSM follows `dyn_ED`
   (dynSubgridDriverMod.F90:313) with `dynSubgrid_wrapup_weight_changes` (:333), whose
   `compute_higher_order_weights` (:389; subgridWeightsMod.F90:216-242) re-derives
   `patch%wtlunit`, `patch%wtgcell` and `col%wtgcell` from the new `wtcol` before
   `reweight_wrapup` (:395). The port ran `dyn_ED` (`fates_apply_deferred_activation!`)
   but not the recompute, so the FATES patches kept their **surfdata** `wtgcell`. An
   inactive slot kept a nonzero gridcell weight (nocomp p6: 0.408), and the bare and
   vegetated patches kept weights unrelated to their FATES area (nocomp bare patch:
   `wtgcell` 0.350 vs `wtcol` 0.880). The port's history writer
   (`_aggregate_to_gridcell`) averaged with those weights and without the
   `patch%active` test of `p2g_1d` (subgridAveMod.F90:531-538). Fix:
   `compute_higher_order_weights!` right after `fates_apply_deferred_activation!` in
   `clm_drv!`. The history writer and `hist_dov2xy` now skip inactive units, zero
   weights and spval, as p2g/c2g do, and `HistFieldDef` carries an `l2g_scale_type`.
   Under FATES-SP, TLAI/TSAI point at `tlai/tsai_hist_patch` (CanopyStateType.F90:232-241).
   Other `wtgcell`/`wtlunit` readers (lnd2atm, drydep, dust) also see the corrected
   weights. The prognostic state should not move, because column averages already used
   `wtcol`. The column fields below are identical before and after to the printed digits.

   Port history-writer gridcell value vs the Fortran h0 (`hist_dov2xy` default), 240
   hourly steps. The figure is the relative error of the window mean (`HIST_*` in
   `fates_modes_parity.jl`).

   | field | SP before | SP after | nocomp before | nocomp after |
   |---|---|---|---|---|
   | BTRAN (instantaneous, l2g `veg`) | 0.82 | 6.6e-4 | 1.51 | 6.8e-5 |
   | TV | 4.5e-3 (1.3 K) | 4.1e-6 | 5.4e-3 | 3.4e-5 |
   | FSH | 0.53 | 5.5e-4 | 0.33 | 1.3e-2 |
   | EFLX_LH_TOT | 0.67 | 9.8e-4 | 0.16 | 6.8e-3 |
   | FSA | 0.52 | 2.5e-6 | 0.37 | 2.0e-6 |
   | FIRE | 0.53 | 2.6e-5 | 0.37 | 8.1e-4 |
   | SABV | 0.75 | 4.2e-5 | 1.22 | 2.6e-4 |
   | ELAI | 0.62 | 1.6e-8 | 0.94 | 2.0e-4 |

   Column-level fields (QRUNOFF, H2OSNO, TSOI_10CM, ...) were already right and are
   unchanged. `test_fates_higher_order_weights.jl` pins the weight invariants and the
   SP history values.

**Default (non-FATES) path.** Items 1, 2, 4 and 5 are FATES-only. Item 3 is a real
default-behaviour change for any caller of `clm_initialize!` + `clm_drv!` that did not
set `use_clm5_fpi` itself. It is not a change for `clm_run!` (and so not for
`parity_run_domain.jl`) or for `fortran_parity_common.jl`; both already set it. Bow
2003 and Stillwater 2003 annual parity are bit-identical before and after. Of the
direct-`clm_drv!` harnesses, Stillwater (single-step and 24-step; snow season) and
the lake harness are identical. Aripuana changes, toward Fortran at the rain steps:

| `fortran_parity_aripuana.jl` step | 3 | 4 | 5 | 6 | 18 | 20 | 22 |
|---|---|---|---|---|---|---|---|
| max rel, before | 1.3e-2 | 3.0e-3 | 7.7e-3 | 5.5e-3 | 5.5e-3 | 7.9e-3 | 1.0e-2 |
| max rel, after | 5.7e-7 | 4.7e-7 | 5.5e-7 | 7.1e-7 | 4.8e-3 | 7.0e-3 | 4.9e-3 |

In `fortran_parity_aripuana_coldstart.jl`, the soil-water parity (`Wtot`) gets worse
at steps 1-5 (step 2: 4.2e-4 -> 7.5e-4) and better at steps 18-24 (step 22: 1.2e-3 ->
6.1e-4). This is not a balance leak; the column still conserves. It comes from the
known cold-start artifact: Fortran's veg filters are empty on the first steps, while
the port intercepts from step 1, and with CLM5 interception the port now holds more
of the early rain on the canopy.

### After (240 steps, relative error of the window mean; before in parentheses)

| field | SP | nocomp + fixed | fixed only |
|---|---|---|---|
| FATES_GPP | +0.69 % (+9.0 %) | +0.56 % (+3.2 %) | +0.99 % (+3.7 %) |
| FSH | +0.15 % (+2.0 %) | +2.0 % (+3.2 %) | +2.5 % (+6.0 %) |
| EFLX_LH_TOT | +0.75 % (+2.9 %) | -0.79 % (-0.8 %) | -0.92 % (-0.5 %) |
| TLAI | 1.7e-8 (1.2e-3) | 4.6e-3 (8.2e-3) | 1.8e-3 (5.0e-3) |
| FATES_VEGC | 7e-9 | 7.4e-4 (3.1e-3) | 1.0e-3 (3.4e-3) |
| cohorts at day 10 | = F | 7 = F (was 8 vs 7) | = F |

SP canopy temperature (TV, 10-day mean): +0.55 K before, +0.02 K after. Per patch, every
patch is within 2e-5 relative.

**Residual:** SP `SNOCAN` about 8 % low, and the first two SP steps (Fortran
`SABV = 0`). Both are closed in the next section.

To reproduce, stage `MODE=sp` with `stage_fates_modes_run.sh`, adding
`hist_dov2xy = .false.` and the extra `hist_fincl1` fields to the staged `lnd_in`.

## Cold canopy, canopy snow and the winter warm bias (2026-10-29)

Three host residuals were left after the section above. In default FATES (carbon only),
day-6 `FATES_TVEG24` was 7.399 C in the port and 7.514 C in Fortran. In SP, `SNOCAN` was
about 8 % low on the tree patches. In the 412-day SPITFIRE run, the port's leaf was 0.5 to
3 K warm from October to May. All three were localized with per-patch hourly h0
(`hist_dov2xy = .false.`) and, where h0 was not enough, a relinked `cesm.exe` with
`write` statements (`build_patched_cesm.jl`: CanopyFluxesMod per iteration,
FrictionVelocityMod `frac_sno`/`z0mg`, SnowCoverFraction SL2012 inputs). There were four
causes: two in `src/` and two in the harness.

1. **Canopy water phase bookkeeping** (`canopy_fluxes.jl`, default path). The port
   re-split the evaporation loss in proportion to the existing `liqcan`/`snocan`.
   CanopyFluxesMod.F90:1602-1620 takes it from `liqcan` first when `t_veg > tfrz`, and
   from `snocan` first otherwise, with the remainder from the other phase. So a warm leaf
   under a snow-laden canopy lost snow that Fortran keeps. It is now
   `canopy_water_evap_update`: the Fortran pair written branch-free and blended by the
   phase step, which is exact for Float64 and smooth (and still mass-conserving) for AD.
   It also adds the `truncate_small_values(snocan, 1e-10)` that was never ported
   (F90:1624-1629). It also passes the FATES `tgcm = thm_patch` (clmfates_interfaceMod.F90:2573,
   2613; the port passed `forc_t`), which is negligible here.
2. **Ground snow roughness after melt-out** (`friction_velocity.jl`, default path).
   `z0mg = zsno` whenever `frac_sno > 0` (FrictionVelocityMod.F90:598). After melt-out
   CTSM's `frac_sno` decays through trace values that come from canopy-snow unloading
   (Fortran, Bow 2002-07-06: 1.8e-5, 9.9e-7, 2.2e-12, 1.5e-17, 1.5e-27 ...). A 1e-11
   floor added in 5cf89ff (for a Savanna subnormal that no longer occurs: zero
   `0 < frac_sno <= 1e-11` events in the 2017 Donga year) switched the port to `zlnd`
   one step early. At nstep 137 that raised `z0mv` from 0.0308 to 0.0437 m, `ustar` from
   0.0644 to 0.0668 and lowered `rah1` from 364 to 338. The leaf ran 0.74 K warm and FSH
   18 W/m2 low that step, and 0.9 K cold for the rest of the afternoon. That was the
   whole day-6 `TVEG24` gap. The floor is removed.
3. **NUOPC albedo timing in the harness.** `fates_modes/cnp/landuse/spitfire_parity.jl`
   passed `doalb = true` from nstep 0. CTSM's cap first solves albedo at the end of
   nstep 2 (`nuopc_doalb`, already used by the other harnesses), so Fortran has `SABV = 0`
   on nsteps 0-2. With `true`, the port's canopy was 5.5 K warm on nstep 1, which put day 1
   `TVEG24` 0.35 K high.
4. **SNICAR files missing from the harness.** No FATES harness passed `fsnowoptics` or
   `fsnowaging`, so the SNICAR Mie tables stayed at zero. With full snow cover the bare
   patch absorbed 12-34 W/m2 more SW (SABG) from October on, and TG/TV ran 1-3 K warm.
   `fates_fortran_parity.jl` now defines `FSNOWOPTICS` / `FSNOWAGING` (the files in the
   Fortran `lnd_in`), and every FATES harness passes them.

Results (Fortran ctsm5.3.012; before -> after):

| quantity | before | after |
|---|---|---|
| carbon-only `FATES_TVEG24` day 1 / day 6 (F 2.071 / 7.514) | 2.425 / 7.399 | 2.071 / 7.512 |
| carbon-only `FATES_TVEG24` days 8-12, max abs diff | | 0.005 K |
| SP per patch, 10-day mean rel.: SNOCAN trees / grass | -9.0 % / -6.7 % | +0.003 % / +0.03 % |
| SP per patch: FCEV trees, FSH_V grass, SABV trees | +3.4 %, +9.3 %, +0.94 % | +0.02 %, -0.24 %, -0.001 % |
| SP `fates_modes_parity.jl`, rel(mean): FATES_GPP / FSH / EFLX_LH_TOT | +0.69 % / +0.15 % / +0.75 % | 5e-5 / 5e-4 / 1e-3 |
| SPITFIRE 412 d, TVEG24 J-F by 30-day block, Oct-Jan | +0.4 to +1.4 K | +0.04 to +0.17 K |
| SPITFIRE 412 d, same, Feb-May | +1.3 to +3.2 K | +0.3 to +1.1 K |

The SP before/after numbers include the harness `doalb` fix.

**The CNP leaf-off bifurcation is Fortran round-off on the same trace branch.** Fortran
CNP and Fortran carbon-only differ only at round-off level. In Fortran CNP the trace
`frac_sno` collapses to exactly 0 at nstep 139, so `z0mg` switches to `zlnd` there. In
carbon-only it holds until nstep 148 (instrumented runs of both cases). That one switch
is the 7.369 vs 7.514 C day-6 split, and so the day-11 leaf-off in Fortran CNP. The port
now follows the carbon-only path in both modes (7.512 C), so it does not reproduce
Fortran CNP's leaf-off.

**Still open.**
* Day-7 `TVEG24` is 12.404 C in the port vs 12.266 C in Fortran carbon-only (12.149 in
  Fortran CNP). It is the same trace branch again: the trace `frac_sno` comes from
  unloading of about 1e-16 mm of snow on the tiny FATES patch 3. At nstep 138 that
  residue is a near-total cancellation (Fortran 6.6e-18, port 1.9e-18), so the port's
  `frac_sno` reaches 0 at nstep 145 and Fortran's at nstep 148. Both models do the same
  arithmetic; the result is ill-conditioned.
* Melt season (Feb-May) is still +0.3 to +1.1 K warm, with single days up to 5 K. In
  March the tree patch absorbs about 7 W/m2 more SABV and holds more canopy snow.
  Not localized further. **Closed below ("Melt season").**
* 412 days into SPITFIRE: `FATES_VEGC` +7 %, live-grass fuel +35 %.

Reproduce: stage the carbon-only control (`fates_cnp_ctrl_c1` namelists) with
`hist_dov2xy = .false.` and the per-patch fields above. Score the port with
`FATES_CNP_PARTEH=1 fates_cnp_parity.jl <h0> 240`.

## Melt season: snow burial, per-patch `tveg24`, case snow namelist (2026-10-29)

Localized with a per-patch hourly Fortran tape over the SPITFIRE case (fates_fire_1pt
namelists, 335 days, `hist_fincl2` hourly with `hist_dov2xy = .false.`: TV, SABV, SABG,
ELAI/TLAI, SNOCAN, unloading, SNOWDP, H2OSNO, FSNO, ALBD/ALBGRD, ...) against the same
quantities dumped from the port each step. In March 2003 both FATES patches had the same
TLAI, HTOP, SWE and snow cover in both models, but the port's exposed LAI was far higher
(tree patch `elai` 0.41-0.55 vs 0.24-0.35; the small patch 0.14-0.21 vs the 0.1 floor),
and its snow was 6-15 % shallower. Three causes:

1. **The snowpack never reached FATES** (`src/`, FATES path). Fortran
   `wrap_update_hlmfates_dyn` (clmfates_interfaceMod.F90:1494-1504) opens with
   `bc_in%snow_depth_si = snow_depth(c)`, `bc_in%frac_sno_eff_si = frac_sno_eff(c)` and
   `UpdateFatesAvgSnowDepth` (EDCanopyStructureMod.F90:1478, `site%snow_depth = depth *
   fsno`). The port ported `UpdateFatesAvgSnowDepth!` but never called it or packed the
   two fields, so `site%snow_depth` stayed 0: the leaf-area profile buried nothing,
   `elai = tlai` all winter, and the canopy absorbed 2-15x the Fortran SABV
   (small patch, March: 11 vs 0.7 W/m2). `fates_update_snow_burial!` now runs where
   Fortran calls `wrap_update_hlmfates_dyn`: after the daily dynamics, at the cold
   start, and at a restart (pack only; the restarted `snow_depth` is kept).
2. **`tveg24` per patch** (`src/`, FATES path). `WrapUpdateFatesRmean`
   (clmfates_interfaceMod.F90:3227-3245) packs every patch's `t_veg_patch` into
   `t_veg_pa(ifp)` before `UpdateFatesRMeansTStep`, which rolls each patch's
   `tveg24`/`tveg_lpa`/`tveg_longterm` from its own value (FatesInterfaceMod.F90:2119-2121).
   The port drove every patch from the first finite `t_veg_pa`, the value the
   photosynthesis iteration had left for patch 1, not each patch's final leaf
   temperature. After fix 1 this was the remaining +0.17..0.34 K Feb-May bias.
3. **Snow namelist in the FATES harnesses** (harness). Every Fortran FATES case lnd_in
   (fates_parity_1pt, fates_fire_1pt, fates_cnp_*, fates_bgc_1pt_v2, fates_sp/nocomp/
   fixedbg/hydro/lu_*) sets `snow_overburden_compaction_method = 'Vionnet2012'` and
   `wind_dependent_snow_density = .true.`. The port defaults are Anderson1976 and no wind
   term (the smoother AD path); `parity_run_domain.jl` already set the Fortran values,
   the FATES harnesses did not. `fates_case_snow_namelist!` (fates_fortran_parity.jl)
   is now called by every FATES build (build, build_cnp, build_mode, build_hydro, the
   land-use build). Snow depth, March: 0.851/0.851 m (was 0.76/0.89).

SPITFIRE mode-1 case, `FATES_TVEG24` J - F by 30-day block (Fortran: the committed
`fates_spitfire_fortran.csv`):

| days | before (main) | + case snow namelist | + snow burial | + per-patch tveg24 |
|---|---|---|---|---|
| 181-210 | +0.17 (max 0.43) | +0.12 | +0.08 | **+0.00 (max 0.04)** |
| 211-240 | +0.33 (0.73) | +0.28 | +0.17 | **-0.00 (0.01)** |
| 241-270 | +0.79 (1.74) | +0.76 | +0.26 | **-0.01 (0.13)** |
| 271-300 | +1.15 (2.92) | +1.07 | +0.34 | **+0.00 (0.00)** |
| 301-330 | +0.86 (5.27) | +0.77 | +0.29 | **+0.00 (0.06)** |
| 331-360 | | | +0.09 | +0.05 (0.55) |
| `FATES_VEGC` J/F, day 330 / 412 | 1.048 / | 1.047 | 1.0002 | **1.0000 / 0.9996** |

The per-patch hourly means (20-day blocks) then agree to the last digit shown for TV,
TG, SABV, SABG, FSA, SNOCAN, SNOWDP, H2OSNO and FSNO on both FATES patches.

**Effect on the 4-year demography** (docs/FATES_MULTIYEAR_ENERGY.md, last follow-up):
the year-2..4 excess was this canopy. Recruitment J/F goes from 1.16/1.57/1.53/2.16 to
1.00/0.90/0.99/1.01, NPP from 1.00/1.30/1.42/1.40 to 0.99/0.98/0.99/0.98, the seed bank
from 1.8-2.4x to 0.94-1.02x, and the PFT-14 summer collapse (day 721: 0.74 vs 9.2 m-2)
to 8.3 vs 9.2 m-2.

**SPITFIRE fuel** (412 days, mode 1, `fates_spitfire_parity.jl`): the fuel-load residual of
the section "Fuel-load residual localized" is gone. Days 380-412 live grass (FC6)
0.0490 vs 0.0488 (was 0.0709 vs 0.0488), fuel amount mean rel. error 1.1e-3, peak fire
intensity 24.6 vs 24.3 kW/m (was 32.4), no burn days in either model.

Reproduce (scratch tooling, not committed): stage the fates_fire_1pt namelists with an
hourly `hist_fincl2` / `hist_dov2xy = .false.` tape and run the stock `cesm.exe`; dump
the port from `fates_spitfire_parity.jl`'s `build()` with the same per-patch fields.

---


---

# PLANT HYDRAULICS (2026-10-29): `use_fates_planthydro`, 1-D Taylor solver

## What was missing
The solver, the water transfer functions and the FATES-side lifecycle were ported, but
the host never ran them: no InitHydroGlobals / InitHydrSites / HydrSiteColdStart at init,
`hydraulics_drive` ran without the host's `total_plant_stored_h2o_col` update, the FATES
root uptake never reached CLM's soil (CLM kept removing `rootr x qflx_tran_veg`), and
`UpdateH2OVeg` / `RestartHydrStates` were not called / not ported.

## What is connected now (Fortran reference -> port)
| Fortran | port |
|---|---|
| `SetFatesGlobalElements2` -> `InitHydroGlobals` | `clm_fates_init!(use_planthydro=true)` |
| `clm_fates%init` -> `InitHydrSites` (after `init_soil_depths`) | `clm_fates_init!` |
| `init_coldstart` soil pack + `HydrSiteColdStart` (before `init_patches`) | `fates_pack_bcin_hydro_soil!(coldstart=true)` + `HydrSiteColdStart` |
| `prep_canopyfluxes`: `qflx_transp_pa = 0` | `clm_drv!` canopy-flux setup |
| `wrap_hydraulics_drive` (CanopyFluxesMod.F90:1635) | `fates_hydraulics_step!` (pack, drive, `total_plant_stored_h2o_col`, hydraulics history) |
| `ComputeRootSoilFlux` (HydrologyNoDrainageMod.F90:343) | `fates_compute_root_soil_flux!` between the root sink and `soil_water!` |
| `dynamics_driv` soil hydraulic pack (:1184) | `fates_daily_dynamics_step!` |
| `update_hlm_dynamics` -> `UpdateH2OVeg` (EDCanopyStructureMod.F90:2073) | `update_hlm_dynamics!` |
| `wrap_update_hlmfates_dyn` plant water (:1552) | daily step + cold start |
| `RestartHydrStates` (:312-528) | `RestartHydrStates` (shares `_hydr_set_soil_wtfs!` with the cold start) |
| BalanceCheck gridcell bypass under planthydro (CTSM#1314) | `balance_check!(use_fates_planthydro=...)` |

Two Fortran behaviours the port reproduces deliberately:
* `call UpdateH2OVeg(sites(s),bc_out(s),bc_out(s)%plant_stored_h2o_si,1)` passes the
  "previous" water by reference, aliased to the value UpdateH2OVeg zeroes and rewrites,
  so its conservation check compares the new value with itself and never fires (it could
  not hold at the cold start anyway). The port re-diagnoses without the check.
* Grasses (`allom_smode = 2`): ctsm5.3.012's FATES `bsap_allom` never assigns its
  intent(out) `sapw_area`, and nothing guards it (no `woody` checks in the hydraulics, no
  conductance clamp, no `-finit-*` flags in the `-O` build), so `UpdatePlantKmax` /
  `UpdatePlantHydrLenVol` read uninitialized stack memory. Measured by relinking the
  stock `cesm.exe` from its own objects: an unmodified relink reproduces the stock h0
  bit-for-bit; adding one `write` of `a_sapwood` shows 7.46 m2 in UpdatePlantKmax and
  ~1e-314 (sometimes ~0.3) in UpdatePlantHydrLenVol, and that write alone changes the run
  (H2OSOI by up to 0.024). The stock planthydro run is therefore not reproducible and
  cannot be a parity target. Upstream FATES fixed it with `SapwoodAreaGrass` (basal
  cross-section pi (d/2)^2); the port carries that fix (FatesAllometryMod.jl).

## The reference is a PATCHED Fortran
Because the stock ctsm5.3.012 result depends on uninitialized memory, the planthydro
reference is produced by `cesm.exe` relinked with one SourceMod,
`SourceMods_planthydro/FatesAllometryMod.F90`: the ctsm5.3.012 file plus upstream
FATES's `SapwoodAreaGrass` copied verbatim from NGEET/fates `main` and called from the
`case(2)` branch of `bsap_allom` (`diff` against the stock file shows only those lines).
Woody PFTs and every non-planthydro path are unaffected (nothing else reads the grass
`sapw_area`).

```bash
julia scripts/validation/fates_fortran_parity/build_patched_cesm.jl \
    $SYMFLUENCE_DATA/runs/fates_hydro_patched/exeA \
    scripts/validation/fates_fortran_parity/SourceMods_planthydro/FatesAllometryMod.F90
EXE=$SYMFLUENCE_DATA/runs/fates_hydro_patched/exeA/cesm.exe STEPS=240 \
    bash scripts/validation/fates_fortran_parity/stage_planthydro_run.sh
julia --project=. scripts/validation/fates_fortran_parity/fates_planthydro_parity.jl \
    export $RUN/fates_hydro_1pt.clm2.h0.2002-07-01-00000.nc
```
`build_patched_cesm.jl` recompiles the given sources with the exact commands in the stock
case's build log and relinks a private executable; the stock build tree is only read.
Determinism (`compare_h0_bitwise.jl`, all 591 h0 variables): two runs of the patched
executable are bit-identical, and so is a third built with an extra `write` of
`a_sapwood` in FatesPlantHydraulicsMod (which prints the grass area as pi (d/2)^2, e.g.
1.644e-5 m2 at d = 0.458 cm).

## Fortran reference and result (patched reference)
`fates_planthydro_fortran.csv` is the 240-step hourly extract of the patched run. The
table gives the median per-step relative error, with ratio = sum Julia / sum Fortran over
the window in parentheses. Values are after the solver fix below; where the pre-fix value
differs it follows in brackets.

| field | day 1 (n 2-24) | days 2-3 (n 25-72) | days 3-10 (n 73-240) |
|---|---|---|---|
| leaf psi, PFT 1 | 1.1e-7 | 2.9e-5 [1.0e-3] | 0.37 |
| leaf psi, woody mean | 1.9e-7 | 2.0e-3 [1.2e-2] | 0.86 [0.69] |
| leaf psi, grass PFT 12 / 14 | 9e-7 / 1.3e-6 | 4.2e-3 / 1.8e-2 [6.2e-2 / 0.16] | 0.48 / 0.09 |
| leaf / stem water content | 4e-8 / 2e-8 | 5.4e-4 / 3.2e-4 [2.5e-3 / 1.5e-3] | 0.09 / 0.06 |
| cohort hydraulic btran | 4.9e-7 | 2.3e-3 [3.2e-2] | 0.57 [0.73] |
| FATES_VEGH2O | 4.6e-8 | 1.2e-3 | 3.4e-2 |
| QROOTSINK layer 1 (ratio) | 2.6e-3 (1.00) | 2.7e-3 (0.99) | 0.18 (1.05) [(0.97)] |
| SOILLIQ layer 1 / 4 | 5e-5 / 1.8e-4 | 3.6e-4 / 4.7e-3 | 6e-4 / 5e-4 |
| EFLX_LH_TOT (ratio) | 5e-4 (1.01) | 1.3e-3 (1.01) | 5.5e-3 (1.00) |
| FATES_GPP (ratio) | 4.8e-4 (1.00) | 7.1e-3 (1.02) [2.8e-2 (0.88)] | 0.39 (**0.93**) [0.62 (0.31)] |
| QVEGT (ratio) | 2.4e-3 (1.01) [(0.89)] | 2.3e-2 (0.98) [0.11 (0.71)] | 0.60 (**0.85**) [0.77 (0.15)] |
| BTRAN (host, ratio) | 4.6e-8 (1.00) [0.60 (1.60)] | 9.7e-6 (1.00) [0.61 (1.66)] | 0.51 (0.97) [0.82 (1.54)] |

Over the whole 10 days GPP is 0.94x Fortran (was 0.33x) and transpiration 0.85x (was
0.16x). These numbers were measured on the integrated tree, which includes #350's
cold-start albedo / thaw fixes, so day 1 is tighter than the numbers first published with
#345.

### Root cause of the post-day-3 collapse: Fortran argument aliasing
`ImTaylorSolve1D` builds each path's Taylor terms with
`call GetImTaylorKAB(..., dftc_dtheta_node(i_up), dftc_dtheta_node(i_dn), ...)`.
* Those two dummies are `intent(inout)` (FatesPlantHydraulicsMod.F90:4117).
* The upstream-K branch zeroes the downstream-side derivative in place (:4144-4152).
* On exit the routine restores `ftc_up` / `ftc_dn` (:4181-4182), but not the derivatives.

So in Fortran, when path `j` has `h_up <= h_dn`, node `j+1` enters path `j+1` (where it
is the downstream node) with `dftc_dtheta = 0`, which changes path `j+1`'s A term. The
port passed values, so path `j+1` saw the full derivative. `_imtaylor_kab_path!` now
reproduces the aliasing.

How it was found:
1. `hydraulics_bc` was instrumented with a per-call, per-cohort dump
   (`instrument_planthydro.py`). It records the
   shells after `FillDrainRhizShells`, the host liquid, `qflx_transp_pa`, the cohort
   pre/post water contents, psi, btran, iteration counts, geometry and kmax. The relinked
   `cesm.exe` reproduces the reference h0 bit for bit.
2. The port was then run as a teacher-forced replay (`fates_planthydro_replay.jl`). At every call the shells, cohort
   water and patch transpiration were overwritten with the Fortran pre-solve values, and
   the post-solve change in water content was compared.

Results of the replay:
* **Before the fix** the solve differs at the very first call (8.6e-5 relative in
  d(theta)), reaches 0.29 by call 3, and reaches up to 2.6e4 by call 23.
* **After the fix** it matches to 1e-7 to 1e-10 for calls 1-22.
* **From call 23** the remaining differences follow the host's transpiration input. They
  are <= 1e-3, and larger only where d(theta) is ~0. `qflx_transp_pa` is 0.1-0.3 % off
  (a host canopy residual). `g_sb_laweight` carries a uniform 1.1e-3 offset, which cancels
  in the cohort partition.

Why one derivative term mattered so much:
* The Fortran solve takes a single 1-hour substep on every call (`iterh1 = iterh2 = 1`).
* Its water-balance test passes by construction, because the tridiagonal flux form
  conserves mass. So the substep is never refined, even when a leaf's water content jumps
  0.24 -> 0.61 in one step (PFT 3, nstep 73).
* The grass leaves have tiny volumes and oscillate hour to hour, even at night.
* Fortran's "erratic" per-PFT btran bursts (0 -> 0.74) are this solver's own behaviour,
  and the port reproduces them once the A terms match.

After day 3 those bursts carry most of the Fortran GPP and transpiration. That is why the
copy semantics showed up as a 0.3x / 0.16x shortfall.

### Still open
* **Post-day-3 decorrelation.** After day 3 the per-step plant states decorrelate (the
  median errors in the table), although the window sums agree. The solve is
  ill-conditioned enough that the 0.1-0.3 % host transpiration residual shifts the timing
  of individual bursts, so day-by-day GPP ratios scatter between 0.5 and 1.9.
* **Host `BTRAN`: closed, it was the harness.** The bracketed BTRAN values above are the
  old harness, which averaged `btran_patch` over the FATES `is_veg` patches only. The
  history `BTRAN` is `l2g_scale_type='veg'` (EnergyFluxType.F90:663), which selects the
  soil/crop landunits (subgridAveMod.F90:1077-1079), not `patch%is_veg`. p2c skips only
  inactive, zero-weight or spval patches (subgridAveMod.F90:203-204). So the FATES
  bare-ground patch counts. Its btran is 0, not spval: BareGroundFluxes zeroes it on the
  no-exposed-veg filter (BareGroundFluxesMod.F90:284), and `wrap_btran` writes only
  `patchi+1..patchi+npatch` (clmfates_interfaceMod.F90:2509-2519). Host BTRAN is
  therefore the vegetated btran times the canopy fraction (0.625 here, = FATES_AREA_PLANTS),
  which is the 1/0.625 = 1.60 the old harness showed. The port's host state was already
  right (bare patch active, wtcol 0.375, btran 0; vegetated patch = `bc_out%btran_pa`,
  the balive-weighted `BTranForHLMDiagnosticsFromCohortHydr`). With the harness fixed,
  day 1 matches to float32 precision. The days 3-10 error is the plant-state
  decorrelation above.
* **2-D solvers.** `MatSolve2D` and `PicardSolve2D` are still not ported.
* **Near-zero fluxes.** The site is extremely water-stressed (cold-start frozen soil; leaf
  psi ~ -15 MPa). Transpiration is ~1e-8 mm/s, so the QVEGT and BTRAN relative errors are
  on near-zero values.
* **Water balance.** With hydraulics on (plant water in the budget), the column balance
  |errh2o| is < 1e-12 mm on every step n >= 1 of the 240-step run. Fortran's daily
  plant-water conservation check is the inert one described above.

---

# LAND USE (2026-10-29): logging, LUH2 transitions and harvest, lupft, potential veg

`clm_initialize!` takes the CTSM `&clm_inparm` land-use options: `fates_harvest_mode`
(`no_harvest` | `event_code` | `luhdata_area` | `luhdata_mass`), `use_fates_luh`,
`use_fates_lupft`, `use_fates_potentialveg`, `fluh_timeseries`, `flandusepftdat`, plus
`fates_paramfile` (a FATES `.cdl`). All default off; the default path is unchanged.

## What is connected (Fortran -> port)

| CTSM / FATES | port |
|---|---|
| CLMFatesGlobals1 `fates_harvest_mode` -> `use_logging` / `use_lu_harvest` / `num_lu_harvest_cats`; `use_luh2`, `num_luh2_states/transitions`, `use_fates_potentialveg` | `fates_harvest_ctrl`, `_fates_spike_set_ctrlparms!` |
| `dynFATESLandUseChangeMod` (`dynFatesLandUseInit` / `dynFatesLandUseInterp`: 12 states, 108 transitions, 5 harvest vars, uninterpolated, END_OF_TIMESTEP year, zero before the series) | `src/fates/fates_driver_landuse.jl` (`DynFatesLandUse` in `inst.fates.luh`, on the ported `dyn_file_io.jl`) |
| `GetLandusePFTData` (`frac_primr/secnd/pastr/range`, `frac_brgnd`) -> `bc_in%pft_areafrac_lu` / `baregroundfrac` under fixed biogeog + LUH | `get_landuse_pft_data!`, `clm_fates_init!` |
| `init_coldstart` / `dynamics_driv` pack of `hlm_luh_*` and `hlm_harvest_*` | `fates_pack_bcin_landuse!` (cold start and every daily step) |
| build-namelist: `use_fates_lupft` => luh + nocomp + fixed biogeog (+ flandusepftdat); LUH harvest modes and potentialveg => luh; potentialveg reads no file | `clm_initialize!` |

`use_fates_lupft` has no runtime switch in CTSM either (no `hlm_use_lupft`): it only sets
those namelist defaults. `fates_harvest_mode = 'landuse_timeseries'` (harvest from the CLM
`flanduse_timeseries` via dynHarvestMod) is not connected and errors.

## Data

No LUH2 files were staged. The files `namelist_defaults_ctsm.xml` names for FATES are on
the CESM inputdata server at 4x5: the 1850-2015 LUH2 timeseries (339 MB) and the landuse x
PFT map (1.9 MB). Both were fetched and cut to the 4x5 cell nearest the site (50N 245E) by
`make_luh_1pt.jl`. The 1x1 results (226 kB) are committed in `data/fates/luh_1pt/` and used
by both models. The 2002 states of that cell are 35 % primary forest, 28 % secondary
forest, 18 % C3 annual crop and smaller pasture / rangeland / urban / other-crop shares, so the cold
start builds five land-use patches. Annual transitions are 1e-8 to 7e-4, and the cell has
`primf_harv` = 3.1e-3, `primf_bioh` = 9.2e6 kgC.

## Fortran reference and result

`stage_fates_landuse_run.sh` (`MODE=logging|luh|luh_area|luh_mass|lupft|potveg`) runs
the archived case. The FATES parameter file is the default with the logging event moved
to 2002-07-04 (`fates_landuse_logging_event_code = 20020704`; the default -30 is day 30)
and `fates_landuse_logging_dbhmin = 0.1` cm, so that the cold-start saplings qualify. The
same `.cdl` edit is applied on the Julia side. The scored series are committed as
`fates_lu_<mode>_fortran.csv`; score them with `fates_landuse_parity.jl <mode>`. The land-use
diagnostics are recomputed from the port's site state with the Fortran history formulas.
Relative error of the 10-day window mean:

| field | logging | luh (= luh_area = luh_mass) | lupft | potveg |
|---|---|---|---|---|
| patches / secondary patches | = F every step | 7 / 2, = F | 13 / 4, = F | 2 / 0, = F |
| patch area by land use | 1.5e-3 - 1.9e-3 | 1e-9 - 4e-9 | 2e-9 - 2e-8 | (primary only) |
| secondary forest fraction | 1.9e-3 | 1e-9 | 2e-9 | - |
| logging disturbance rate | 1.9e-3 | (no event) | - | - |
| LU-transition matrix / LUC disturbance | - | 1.4e-7 / 7e-9 | 1.4e-7 / 3e-8 | - |
| LUC wood product flux / stock | - | 1.2e-3 / 9e-4 | 1.6e-3 / 1.2e-3 | - |
| FATES_VEGC | 1.9e-3 | 1.7e-3 | 1.2e-3 | 2.0e-3 |
| FATES_LITTER_IN | 2.1e-3 | 4.7e-2 | 1.0e-3 | 5.9e-2 |
| FATES_GPP | +3.9 % | +3.2 % | +3.1 % | +3.9 % |

* **Logging.** The 2002-07-04 event moves 63.4 % of the site to secondary land (F
  0.63439, J 0.63559). Its litter pulse matches to 1e-3. One fix was needed: in event-code
  mode ctsm5.3.012's FATES gates direct logging on `cur_harvest_tag == 0` and never sets
  the tag (an uninitialized local). The port had chosen 0, which turns direct logging on;
  the reference `cesm.exe` behaves as "not applicable" (2), so there is no direct
  mortality. With 0 the port lost 11 % more vegetation C and had 2.4x the logging-day
  litter; with 2 both agree to 2e-3. Later FATES sets a bypass tag here and does log
  directly. Harvested wood is 0 in both models: the saplings carry no trunk-size CWD.
* **LUH2 harvest is a no-op in ctsm5.3.012.** The host passes the LUH2 names
  (`primf_harv`, ..., `primf_bioh`, ...) as `hlm_harvest_catnames`, but FATES'
  `get_harvest_rate_area` / `get_harvest_rate_carbon` match `HARVEST_VH1` ... `HARVEST_SH3`.
  Every rate comes out 0, and `luhdata_area` / `luhdata_mass` give the same result as `luh`
  in Fortran, and in the port, which keeps the names verbatim.
* **lupft** builds 13 nocomp patches across the five land-use types from the landuse x
  PFT map, with areas matching to ~1e-8.
* **potveg.** Two Fortran runs of the identical potveg case differed: the first one
  diverged from day 2 (litter -13 %, treefall -29 % by day 3), and a rerun matches the
  no-LUH case. That is nondeterminism in the reference build, not in the port, and
  the rerun is committed. Against the rerun, potveg matches the default mode, with one
  primary patch and no transitions.

The GPP / FSH offsets, and the 5 % day-2 litter gap in the non-logging modes, are the host
residual that the default mode also shows (see the soil-BGC section). Not validated: restarts with land use, potential-veg -> LUH switch-on
(`transition_landuse_from_off_to_on`), carbon-mode harvest debt with nonzero rates (they
cannot occur in this release), and runs longer than 10 days.

---

# SPITFIRE FIRE (2026-10-28): host fire weather + mode 1 vs Fortran

## What changed on the host
* **24-hour fire-weather means.** SPITFIRE reads `bc_in%precip24_pa / relhumid24_pa /
  wind24_pa` (and the FATES `tveg24` running mean). The port used to fill them from a
  single-step stand-in in `clm_driver.jl` (this step's rain, a Tetens RH, `|u,v|`), the
  same value for every patch. They now come from the host accumulators, as in Fortran:
  `PREC24` / `RH24` (Wateratm2lndBulkType) and `WIND24` (atm2lndType) are accumulMod
  `runmean`s with a 1-day period and init values 0 / 100 / 0. Only active patches update,
  each with its own step counter (`fates24_nsteps_patch`), so a FATES patch that switches
  on mid-run starts a fresh mean. `dynamics_driv` copies `prec24/rh24/wind24_patch(ifp +
  patchi(c))` into `bc_in` for every vegetated FATES patch (`fates_pack_bcin_daily!`).
  The fields are allocated only under `use_fates`.
* **Modes 2-5 host (FATESFireDataMod).** `clm_initialize!(fates_spitfire_mode>1, flnfm,
  fhdm)` reads the lightning / population-density streams (`FATESFireData`),
  `FireInterp` runs every step, `lnfm24` is a 1-day runmean, and the daily pack sets
  `lightning24 = lnfm24*24` (all modes > 1) and `pop_density = forc_hdm` (modes >= 3).
* `clm_initialize!(fates_spitfire_mode=...)` now reaches `clm_fates_init!`.
* The FATES clock fix from #341 (`fates_get_and_set_time!`, cherry-picked) is required:
  without it FATES sits on 2 January and the fire-season comparison is meaningless.

## Fortran reference
`stage_spitfire_run.sh` runs the archived `fates_parity_1pt` case (Bow at Banff,
2002-07-01, 3600 s, use_fates_bgc) with `fates_spitfire_mode=1` for 730 days and a
daily-mean h0 of the fire fields. Variants: `CDL=<cdl>` swaps the FATES parameter file;
`SPITFIRE_MODE=2..5` adds the synthetic uniform lightning/popdens streams (0.000457
strikes/km2/hr). Committed extracts: `fates_spitfire_fortran.csv` (mode 1, default
params), `fates_spitfire_thr5_fortran.csv` (mode 1, `fates_fire_threshold` 50 -> 5 kW/m,
450 days) and `fates_spitfire_m2thr5_fortran.csv` (mode 2, same threshold).
`fates_spitfire_parity.jl` steps the port over the same window and scores it.

With default parameters the Fortran site never burns: fire intensity peaks at
24.3 kW/m (2003-08-04..16, days 400-412, the 2003 dry spell), under the 50 kW/m
`fates_fire_threshold`. The threshold-5 variant burns on days 402-412 and exercises
area burnt, crown scorch, cambial kill, fire mortality and the fire carbon flux.

Note on forcing alignment: DATM indexes the standard-calendar `clmforc.2004.nc` by
elapsed time under the NO_LEAP model calendar, so from 2004-03-01 the model reads the
previous calendar day. The parity script reproduces this (`datm_time`); aligning by
date puts every day after 2004-02-28 one day off (TBOT ~2 K).

## Result (mode 1, default parameters, 730 days)
Host drivers (daily means): WIND, PRECIP, TBOT within 6e-8, RH within 1.1e-3.

Fire weather and fuel moisture, which depend only on weather and `tveg24`, track:

| window | NI mean F / J | FDI mean F / J | fuel moisture by class (mean rel) |
|---|---|---|---|
| days 1-40 | 275.9 / 280.1 (1.5 %) | 0.0876 / 0.0888 | 2e-5 .. 4e-3 |
| days 380-420 (2003 fire season) | 2673 / 2672 (0.04 %) | 0.5415 / 0.5401 | |

The Nesterov index rises, peaks (5611 / 5627 on day 412) and resets (day 413) on the
same days. The day-1-10 NI gap comes from `TVEG24` (2.07 / 2.84 C on day 2), the
cold-start canopy temperature residual, not from the fire code.

Fuel **amount** and **composition** do not match, and the gap grows with the vegetation
gap. FATES_VEGC and FATES_NPATCHES separate after about 5 weeks (Julia makes 9 patches by
day 90, Fortran waits until day ~400; Fortran drops 0.02 kgC/m2 of live carbon into leaf
litter on day 39, Julia does not). This happens identically with fire off (Fortran fire-on
VEGC equals fire-off VEGC to every digit until the first burn), so it is the FATES
demography residual tracked in `docs/FATES_MULTIYEAR_ENERGY.md`, not SPITFIRE. During the
2003 fire season Julia carries 27 % more ground fuel (live grass +25 %), a 7 % lower
effective wind (fewer trees), and therefore about 1.7x the ROS and 2.2x the intensity.
It crosses the 50 kW/m threshold on days 409-412 while Fortran (24 kW/m) does not:

| | Fortran | Julia |
|---|---|---|
| days with FI > 0 | 13 (400-412) | 13 (400-412) |
| peak FI | 24.3 kW/m (day 411) | 46.5 kW/m (day 411) |
| burn days (FI > 50 kW/m) | 0 | 4 (409-412), total burned fraction 1.2e-4 |

SFMainMod (fire weather, ROS, fuel consumption, area burnt, scorch, crown/cambial damage,
post-fire mortality), SFNesterovMod and SFFireWeatherMod were re-read line by line
against the Fortran; no discrepancy was found. FatesFuelMod was not re-read, but its
per-class fuel moisture output matches to < 0.5 %. The difference in the burn decision
comes from the fuel load.

## Burning chain (mode 1, `fates_fire_threshold` = 5 kW/m, 450 days)
With the lowered threshold both models burn during the 2003 dry spell:

| | Fortran | Julia |
|---|---|---|
| burn days | 11 (402-412) | 12 (401-412) |
| total burned fraction (days 400-413) | 1.29e-4 | 3.42e-4 (2.65x) |
| fire C loss to atmosphere | 4.6e-6 kgC/m2 | 1.8e-5 kgC/m2 (3.8x) |
| fuel consumed (sum of daily) | 1.089 | 1.376 (1.26x) |
| fire mortality, FATES_MORTALITY_FIRE_SZ summed | 0.0129 | 0.0495 m-2 yr-1 |

The ratios follow from the fuel load. Julia carries 27 % more ground fuel, which gives
1.6-2x the ROS. Burned area scales as (ROS x duration)^2, hence 2.65x. Higher intensity
raises scorch height, which drives the crown and cambial kill. The fire chain switches on
and off on the same days in both models.
(`FATES_MORTALITY_FIRE_SZ` is a level-2 history field. The level-1 harness never reaches
the dyn2 reset of the site `fmort_rate_*` accumulators, so the script reports the daily
increment.)

## Mode 2 (lightning from data), threshold 5 kW/m, 450 days
Run with the synthetic uniform lightning stream (0.000457 strikes/km2/hr, so
`lightning24` = 0.01097/km2/day).

**Bug found and fixed.** The port passed the SPITFIRE mode enumeration shifted by one:
`sf_successful_ignitions_def=2`, `sf_anthro_ignitions_def=3`. CTSM uses
FATESFireFactoryMod: 2 = lightning_from_data, 3 = successful_ignitions,
4 = anthro_ignitions. As a result, `fates_spitfire_mode=2` ran as "successful
ignitions" (FDI forced to 1, cloud-to-ground fraction 1), and mode 3 ran as
anthropogenic ignitions. With the fix, FDI matches Fortran (day 405: 0.8292 / 0.8304),
and the mode-2 to mode-1 burned-area ratio is the same in both models:

| | Fortran | Julia |
|---|---|---|
| burn days | 11 (402-412) | 12 (401-412) |
| total burned fraction | 3.43e-5 | 9.10e-5 |
| burned fraction, mode 2 / mode 1 | 0.267 | 0.266 |

The absolute burned area again carries the 2.65x fuel-driven factor. The ratio isolates
the ignition path, which is the host lightning stream -> `lnfm24` -> `lightning24` ->
`NF`, and it agrees to 0.4 %.

Not compared against Fortran: modes 3-5 (population density, anthropogenic ignitions).
The host path is ported and unit-tested, and the enumeration fix makes the mode numbers
correct. The real `clmforc.Li_2012_climo1995-2011.T62.lnfm_Total_c140423.nc` climatology
is not staged on this box; only the synthetic uniform stream is.

## Fuel-load residual localized (2026-10-29): litter and live grass, not fuel

The script now also scores the aboveground litter the fuel is built from,
`FATES_LITTER_AG_FINE_EL` and `FATES_LITTER_AG_CWD_EL` (extract refreshed; the old
columns are byte-identical). FatesFuelMod and SFMainMod `UpdateFuelCharacteristics`
(SFMainMod.F90:179-197) were re-read against the Fortran: class map (FatesFuelClassesMod.F90:14-19,
trunks = 4 excluded), `UpdateLiveGrass` (FatesPatchMod.F90:784-812), the
`(1 - SF_val_miner_total)` reduction in `rate_of_spread` that history reads
(FatesHistoryInterfaceMod.F90:2691), patch area weighting, fuse, spawn and the call
order in `ed_ecosystem_dynamics` all match. In both models each fuel class is its
litter pool x 0.945, so the fuel gap is a litter / live-grass gap. A unit test pins
this mapping (`test_fates_sfmain.jl`).

Integration state (`integration/wave-2026-10`), mode 1, default parameters, days 380-412:

| | Fortran | Julia | J/F |
|---|---|---|---|
| fuel (`FATES_FUEL_AMOUNT`) | 0.1625 | 0.1829 | 1.13 |
| twigs + small branches (FC1+FC2) | 0.0468 | 0.0405 | 0.87 |
| dead leaves (FC5) | 0.0670 | 0.0716 | 1.07 |
| live grass (FC6) | 0.0488 | 0.0709 | 1.45 |
| day 410 ag fine / ag CWD litter | 0.0968 / 0.0577 | 0.1025 / 0.0516 | |
| peak FI, burn days (> 50 kW/m) | 24.3 kW/m, 0 | 32.4 kW/m, 2 (410-411) | |

The first divergence is the cold-deciduous leaf-off date. Fortran drops the PFT 3/6
leaves on day 39; Julia drops them on day 11. EDPhysiologyMod `phenology` drops when
more than `phen_ncolddayslim` = 5 of the last 10 daily `tveg24` means are below
`phen_coldtemp` = 7.5 C. On day 11 the memory holds days 2-11, and day 7 decides it:
`FATES_TVEG24` = 7.514 C in Fortran and 7.399 C in Julia, a 0.115 K canopy-temperature
difference. So Julia counts 6 cold days and Fortran 5. `tveg24` is the 24-h mean of the
host `t_veg_patch` (clmfates_interfaceMod.F90:2612/3239 -> FatesInterfaceMod.F90:2119),
so this is a host canopy-temperature residual, not a FATES or fuel defect.

Diagnostic run: `PHEN_TMEM_FROM_FORTRAN=1` replaces each day's phenology memory entry
with the Fortran `FATES_TVEG24`. Julia then drops its leaves on day 39 as Fortran does,
and the fuel tracks Fortran for ten months (day 40 0.0232 / 0.0232, day 240 0.0760 /
0.0766, day 300 0.0914 / 0.0925). The second divergence is spring 2003 live grass. In
days 241-330 (Feb-May) Julia's `TVEG24` runs +3.1 K warm (up to 9 K on single days),
over a +0.5 to +1.4 K warm bias from October onward. Grass flushes and grows earlier
(FC6 2.3x on days 301-330, 3x on days 331-360), and the run lands on the PR #347
picture: fuel x1.23, peak FI 45.8 kW/m, 4 burn days. The +27 % reported before is this
second residual. On the current integration state the early leaf-off partly offsets it
(less woody litter, fewer leaves dropped).

TVEG24 by 30-day block (J - F): days 1-90 +0.00 K, days 91-210 +0.5 to +0.85 K, days 211-240
+1.4 K, days 241-330 +2.7 to +3.2 K, days 331-412 -0.3 to +0.6 K. The remaining fuel residual
is the host winter/spring canopy temperature over snow. Fix it there; SPITFIRE needs no
change.

---

# RESEARCH MODES (2026-10-29): tree damage, ST3, prescribed physiology, cohort age, seed dispersal, inventory init

These six FATES modes were ported but never run. `clm_fates_init!` forced them off.
They are now switchable from `clm_initialize!` with CTSM's namelist names
(`use_fates_tree_damage`, `use_fates_ed_st3`, `use_fates_ed_prescribed_phys`,
`use_fates_cohort_age_tracking`, `fates_seeddisp_cadence`, `use_fates_inventory_init` +
`fates_inventory_ctrl_filename`, and `fates_paramfile` for a parameter variant).
`fates_check_niche_modes` applies CTSM's host-side rules before any state is built:
cadence must be in 0..3, the inventory control file must be set and exist, ST3 excludes
prescribed physiology, inventory excludes cohort age, and tree damage excludes CNP. FATES's
own `check_allset` re-checks the FATES rules. It now runs after the parameter read, as in
`CLMFatesGlobals2`. With every mode off, the 120-step parity dump is byte-identical to the
tree before this change.

## Recipe

```bash
# Fortran: stage a mode from the fates_bgc_1pt_v2 run dir (lnd_in override + daily
# history) and run cesm.exe in it (~3 s for 14 days)
python3 stage_niche_mode.py runs/fates_bgc_1pt_v2 runs/fates_niche_<mode> 14 \
    use_fates_ed_prescribed_phys=.true. [paramfile=<.nc>] [mesh=<.nc>] [fields=A,B]
# parameter variants (seed dispersal, damage events): edit the CDL, ncgen it for Fortran
python3 make_param_variant.py data/fates/fates_params_default.cdl out.cdl name=value ...
# Julia: same case, one CSV row per day
FATES_NICHE_DAYS=14 julia --project=. scripts/fates_niche_modes.jl <mode> out.csv
# score: Julia vs Fortran, and the mode effect (mode - base) in each model
julia --project=. fates_niche_compare.jl jl_mode.csv f_mode.h0.nc jl_base.csv f_base.h0.nc
```

The base runs differ by the known host residual (vegc 0.3 %, FSH 11 % daily mean), so each
mode is scored on its own state and on its effect relative to the matching base run.

## Results (14 days from 2002-07-01, Bow at Banff, all 14 PFTs, FATES soil BGC)

| mode | Julia runs | Fortran parity | fixes |
|---|---|---|---|
| prescribed physiology | yes | ncoh identical every day (14 to 95); vegc/leafc/storec/structc/seed bank within 5e-8 (the float32 history precision) | none needed |
| inventory init | yes | day 0: 9 cohorts, 2 patches, every pool within 5e-8. 14 days: ncoh identical every day (57/57), npatch 9/9, vegc 1.6e-4, structc 1.3e-4 | NaN inventory coage blocked all fusion (Julia 335 cohorts vs 57); host coordinates passed to the proximity check |
| seed dispersal (cadence 1) | yes | seed-bank effect −4.80e-4 vs −4.64e-4 kgC/m² (3.4 %); per-PFT daily `seed_out` and its daily reset track Fortran | host wiring was missing (`WrapUpdateFatesSeedInOut`, outgoing buffer, `WrapGlobalSeedDispersal`) |
| cohort age tracking | yes | ncoh effect +98 vs +95 by day 9; vegc effect < 2e-7 in both | none needed |
| ST3 | yes | Fortran writes no FATES history under ST3 (`update_history_dyn` returns). Julia holds the stand bit-static (14 cohorts, cold-start carbon); host effect at day 14: FSH +2.33 vs +2.18 W/m², HR +3.13e-8 vs +2.97e-8 gC/m²/s; the TV effect is under 0.05 K and differs in sign | none needed |
| tree damage | yes | default `fates_damage_event_code = 1` means no events: both models are identical to their base runs on every field. With events (code 2 or 3) CTSM hangs at the first daily step (>13 min CPU at nstep 1), so there is no Fortran reference | `branch_loss_frac` was scoped inside the element loop, so the first damage event threw `UndefVarError`; Julia now runs with events (code 2 moves 0.8 % of plants to damage class 2, mass conserved) |

Notes:

* **Seed dispersal is single-point here.** One gridcell has no neighbours
  (`DetermineGridCellNeighbors` loops over pairs), so the incoming seed is always 0 and
  the mode reduces to exporting `fates_seed_dispersal_fraction` of the seed rain. The
  inter-gridcell kernel needs a multi-gridcell run, which this case cannot provide. The
  default parameter file leaves the four `fates_seed_dispersal_*` parameters unset, so the
  runs use `make_param_variant.py` with fraction 0.2.
* **Inventory init cannot run unmodified in Fortran at this site.** FATES converts negative
  inventory longitudes to [0, 360] but compares them with the host longitude as-is, so any
  site on a [-180, 180] mesh (Bow is at -116.01) aborts at FatesInventoryInitMod.F90:247.
  The Fortran reference uses a copy of the site mesh with longitude +360 (`mesh=`). That
  copy is not neutral for the host (FSH up to 45 %), which is why inventory is scored
  against `base_lon360` and mainly on the day-0 stand. The port wraps the host longitude
  into [0, 360) instead of aborting. The synthetic pss/css pair (2 patches, 6 cohorts of
  PFTs 2/6/12) is in `runs/fates_niche_inputs/`.
* **Cohort age:** after day 10 the Fortran run takes a knife-edge cold-deciduous leaf-off
  (PFTs 3/6/11 lose ~25 % vegc in one day) that its base run does not take. The unrelated
  `base_lon360` run takes the identical event, so it is a threshold flip, not an age-tracking
  effect. The mode is scored on days 0-10.
* **Tree damage hang:** reproduced with `fates_damage_event_code` 2 and 3 on this case. It
  was not traced further (debugger runs are unreliable on this machine; see setup_case.sh).

## Restart round trip

`write_fates_restart` packs the FATES state (`set_restart_vectors!`) to NetCDF on the
`cohort`/`column` dimensions. `read_fates_restart!` follows the Fortran `restart(flag='read')`
sequence. It rebuilds fresh sites with `create_patchcohort_structure!`, which now calls
`init_site_vars!` + `zero_site!` (the open TODO, FatesRestartInterfaceMod.F90:2856) and
`NanValues` + `ZeroValues` on each new cohort. Then it runs `get_restart_vectors!`,
`ed_update_site(is_restarting=true)`, the `wrap_update_hlmfates_dyn` hand-off,
`update_3dpatch_radiation!`, the dynamics history and the restart seed-dispersal call. Plant
hydraulics restarts (`RestartHydrStates`) are not wired and refuse.

Test: base case, FATES restart written after day 7 and the sites rebuilt from the file,
compared with the continuous run. **All 97 daily diagnostics (FATES state and host
FSH/TV/HR) are bit-identical on all 15 rows.** The first attempt failed the day-8 mass
balance with NaN litter fragmentation. The cause was that restarted cohorts skipped
`ZeroValues`, leaving the fields it sets (fire mortality, the daily N/P demand/uptake
sentinels, ...) as NaN; that is now fixed. The host CLM state is carried in memory, not through a host restart file.

## Found, not fixed (outside these modes)

`site.lat/lon` are never set in the port (NaN), while Fortran sets them from the host
gridcell (clmfates_interfaceMod.F90:884-885). Phenology tests `site.lat > 0`, so at Bow
(51°N) the port takes the southern-hemisphere branch: `gddstart = 181` resets GDD on 1 July
and `ncdstart = 120`. Fixing it changes the default output, so it is left for the default-path
owners.
