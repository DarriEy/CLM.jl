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

Unlike #197, #227's premise can be tested directly, because "prescription OFF" is
already a configured site (`krycklan_baseline`, `soil_h2o = nothing`).

| run | #227's report (then) | current `main` (now) |
|---|---|---|
| `krycklan_baseline` (prescription **OFF**) | cold `749.6` → **`0`**, monotone collapse; `btran = 0` every growing-season day | cold **`2549`** → `4040`; min `2008`; **no collapse, no boom-bust** (both checks PASS at day 235) |
| `krycklan_screened` (prescription **ON**) | carbon `749.6` → `1990`, 10/10 pass | 10/10 pass over the full 365 days |

Two things changed, and both matter:

1. **The collapse is gone without the prescription.** Carbon *rises* 2549 → 4040 and is
   still rising when the run stops. #227's entire justification was a monotone decline
   to zero. It does not happen any more.
2. **The cold start itself is different** — `749.6` then vs `2549` now. #227 was
   measured on a materially different configuration. That is consistent with the driver
   defaults that moved since (`use_bedrock` #251/#252 above all, plus #259, #225, #265,
   #267-#269, #273): #227's root cause was a bone-dry root zone, and `use_bedrock`
   truncates the soil column, which is exactly the mechanism #251/#252 corrected for the
   Bow window. `fates_multisite_validation.jl` **inherits** those defaults rather than
   setting them, so it was silently re-specified by each of those merges.

**Verdict on #227: it is a compensating workaround for a defect that has since been
fixed elsewhere.** It is harness-only (no `src/` change) and gated on `SiteCfg.soil_h2o`,
so it distorts nothing by sitting there — but it should be **retired, or demoted to an
explicitly-labelled diagnostic**, and the boreal band re-baselined from the unprimed run.
The Fortran oracle D3 asked for (a boreal single-point reference cold-started identically
with the prescription off) is no longer the first thing needed: the artifact it was meant
to measure has already stopped reproducing.

Caveat: the baseline run stopped at day 235 of 365 on the same fatal `errlon` (§4), so
this is 235 days of evidence, not a full year. It is nonetheless decisive against a
*monotone collapse from day 1*, which is what #227 described.

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
model, not the assertion. Note the daily FATES carbon balance holds *right up to the
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

# #227 both arms (Krycklan)
SYMFLUENCE_DATA=$SD SITE=krycklan_baseline FATES_NDAYS=365 julia +1.12 --project=. scripts/fates_multisite_validation.jl
SYMFLUENCE_DATA=$SD SITE=krycklan_screened FATES_NDAYS=365 julia +1.12 --project=. scripts/fates_multisite_validation.jl
```

**Harness trap found while doing this:** `fates_multisite_validation.jl` **hardcoded**
the data root, so it could only ever see `/Users/…/compHydro/SYMFLUENCE_data` — which
holds none of the three FATES domains. They live under the Drive root. This is why #227
reported "Tropical Aripuana data is absent in this environment" and skipped the tropical
arm: the domain existed all along. Now honours `SYMFLUENCE_DATA` like every other
harness (`test/testdata.jl::symfluence_data_root`).
