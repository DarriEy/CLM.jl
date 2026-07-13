# Biogeochemistry (CN/BGC) parity scorecard — Bow at Banff

**Status: partially validated.** This document says exactly what is and is not validated.
It deliberately does not generalise beyond what was measured.

Everything below is a Julia-vs-Fortran diff against per-timestep dumps from the *same*
instrumented CTSM `cesm.exe`, at the Bow-at-Banff single column, from a converged BGC
spinup (model year 2202). Nothing was tuned to match.

---

## 1. What the previous CN validation could not see

Before this work the entire CN/BGC evidence base was **one contiguous SUMMER window**
(28 steps, 2202-07-15/16), single-step-injected. Bow's PFTs are:

| patch | PFT | phenology |
|---|---|---|
| 1 | 0 — bare | none |
| 2 | 1 — needleleaf evergreen temperate tree | **evergreen** |
| 3 | 12 — C3 arctic grass | **season_decid = 1** |

In a summer window `offset_flag == 0` on every patch, for every step. So the summer suite
**cannot execute a single line of the leaf-offset / senescence code path** — not the
trigger, not the litterfall ramp, not the ramp-end cleanup. That is precisely where PR
#206's `ndays_off` finding (30 d vs the params file's 15 d) actually bites.

This work generated the **autumn / leaf-offset window** that had never existed, and it
immediately exposed three real bugs (§4).

---

## 2. Fortran reference windows (all newly regenerated or reused)

| window | nsteps | model dates | what it exercises |
|---|---|---|---|
| winter | 1753153 | 2202-01-01 | dormant CN |
| summer | 1757845–1757872 (28) | 2202-07-15/16 | peak growing season |
| **autumn (NEW)** | **1759897–1760376 (480)** | **2202-10-09 → 10-28** | **leaf offset: trigger, 15-day ramp, ramp end, dormancy** |
| **cold start (NEW)** | 0 | 2202-01-01 | `use_cn` InitCold |

The autumn and cold-start references required **restoring the `restFile_write_dump`
instrumentation** to CTSM (see §6) — the writer that produced the original BGC dumps had
been replaced by a 55-variable snow/PHS-only `pdumpMod` that carries **zero CN fields**.

### The autumn window is the right window

Bow's grass is seasonal-deciduous, so leaf offset triggers on `dayl < crit_dayl` past the
solstice. At 51.36°N, daylength crosses `crit_dayl` (39300 s) on **doy 284 = Oct 11**.
The Fortran reference fires offset at **nstep 1759931 = 2202-10-10 10:00**, with
`offset_counter = 15.0 days` exactly, and the ramp ends at nstep 1760291 (10-25 10:00) —
**direct confirmation that `ndays_off = 15`, not 30**. `leafc` falls 35.12 → 0.19 gC/m².

> Note: `crit_dayl` 39200 → 39300 (also fixed in #206) is a **no-op at Bow** — daylength
> falls 233 s/day there, so both thresholds are crossed on the same timestep. `ndays_off`
> 30 → 15 is the change that actually matters, and it is confirmed against Fortran.

---

## 3. Scorecard

### 3a. Single-step parity (inject Fortran `before_step`, run 1 step, diff `after_hydrologydrainage`)

> Methodology note (banked trap): CN probes **must** diff against
> `after_hydrologydrainage`, not `before_step` — `before_step` holds the *previous*
> step's values for the allocation-pipeline fields. `after_hydrologydrainage` is after
> `EcosystemDynamicsPreDrainage` (CN pools updated) but before `PostDrainage` (leaching).

| window | global CN max\|rel\| | worst fields |
|---|---|---|
| winter (1 step) | **5.90e-05** | xsmrpool, smin_nh4_vr |
| summer (1 step) | **2.18e-03** | cpool, smin_nh4_vr / sminn_vr |
| **autumn (10 probes across the ramp)** | **2.16e-03** | frootc_xfer, leafc_xfer, cpool |

Autumn is now **the same order as summer**, and dominated by the *same* fields (the
temporary C buffer pools), i.e. no new physics error is left in the offset path.

**Autumn per-probe (51 CN fields each):**

| probe | max\|rel\| | leafc rel | offset_flag | offset_counter |
|---|---|---|---|---|
| pre-trigger (10-09 03:00) | 8.58e-04 | 1.47e-09 | ok | ok |
| step before trigger (10-10 09:00) | 8.56e-04 | 1.47e-09 | ok | ok |
| **TRIGGER (10-10 10:00)** | 8.56e-04 | 1.47e-09 | **ok** | **ok** |
| ramp d0.6 (10-11) | 8.54e-04 | 1.47e-09 | ok | ok |
| ramp d3.6 (10-14) | 8.99e-04 | 1.47e-09 | ok | ok |
| ramp d7.6 (10-18) | 1.66e-03 | 1.46e-09 | ok | ok |
| ramp d12.6 (10-23) | 2.16e-03 | 1.38e-09 | ok | ok |
| step before ramp end (10-25 09:00) | 1.34e-03 | 2.44e-10 | ok | ok |
| **RAMP END (10-25 10:00)** | 8.23e-04 | 2.44e-10 | **ok** | **ok** |
| post-ramp dormant (10-28) | 8.16e-04 | 2.44e-10 | ok | ok |

The **offset litterfall flux itself matches to `ratio J/F = 1.0000`** at every step of the
ramp (`scripts/probe_offset_litterfall.jl`).

### 3b. Multi-step drift (inject ONCE, free-run, compare trajectory)

**Autumn, 480 steps (20 days), free-running straight through senescence:**

| step | sminn_vr | leafc | soil1c | leafc grass J/F | offset_flag J/F |
|---|---|---|---|---|---|
| 1 | 2.0e-05 | 4.3e-05 | 2.6e-06 | 35.123 / 35.123 | 0/0 |
| 49 | 8.2e-04 | 2.1e-03 | 1.3e-04 | 35.054 / 35.054 | 1/1 |
| 145 | 2.3e-03 | 6.1e-03 | 3.8e-04 | 31.744 / 31.744 | 1/1 |
| 289 | 4.5e-03 | 1.2e-02 | 7.6e-04 | 17.438 / 17.438 | 1/1 |
| 385 | 6.2e-03 | 1.6e-02 | 1.0e-03 | 1.717 / 1.717 | 1/1 |
| 480 | 8.0e-03 | 2.0e-02 | 1.3e-03 | 0.192 / 0.192 | 0/0 |

Drift is **bounded and near-linear** — the compounding of the known per-step residual, not
a runaway. Julia reproduces the whole 15-day senescence trajectory. (Before the §4 fixes
this same run reached a leafc relative error of **182** and never triggered offset at all.)

Summer, 28 steps: sminn_vr 0.66 %, leafc 1.2e-03, soil1c 7.5e-05 — unchanged from the
previously banked numbers (no regression).

### 3c. `use_cn` cold start (never previously possible)

Against a genuine Fortran `start_type=startup` CN cold start:

- **Finiteness: PASS.** Every CN state array on the soil-BGC column is fully finite. (The
  deep-lake column holds NaN where Fortran holds `spval` — the same "inactive" convention;
  the CN filters exclude it.) This is the #212 regression guard.
- **Values: BIT-EXACT (worst rel = 0.0)** across all three patches and every C and N pool —
  after the §4.4 fixes. Verified: bare = all zero; evergreen tree `leafc = frootc = 100`,
  storages 0, `deadstemc = 0.1`, `leafn = 100/58 = 1.7241`; deciduous grass `leafc = 0`,
  storages 100, `leafn_storage = 100/20.7 = 4.8309`.

---

## 4. Real bugs found by the autumn window (all fixed, all verified against Fortran)

### 4.1 The offset litterfall ramp memory was wiped every step — **physics bug**

`_zero_cnveg_flux_arrays!` (`cn_driver.jl`) blanket-`fill!(0)`s every array field of the CN
veg flux struct, with a hand-maintained exception list. `prev_leafc_to_litter` /
`prev_frootc_to_litter` were **not** on it, so they were zeroed at the start of every CN
step. Fortran's `CNVegCarbonFluxType::SetValues` never touches them — they are set only in
`InitCold` and are restart variables (`CNVegCarbonFluxType.F90:3986/3988`, `4291/4296`),
reset in exactly two places (offset trigger, ramp end), both in `CNPhenologyMod`.

They carry the ramp memory:

```
t1 = dt*2/offset_counter²
leafc_to_litter = prev_leafc_to_litter + t1*(leafc − prev_leafc_to_litter*offset_counter)
```

With `prev ≡ 0` this collapses to `t1*leafc` — measured at **1.7 %–7 % of the correct flux**
mid-ramp. Leaves were therefore not shed gradually over `ndays_off`; the pool sat ~full for
the whole offset period. **Completely invisible in a summer window.**

### 4.2 `read_fortran_restart!` dropped `prev_{leafc,frootc}_to_litter`

Genuine Fortran restart fields, never mapped by the Julia reader — so any restart/injection
resumed a mid-ramp senescence with the memory term at zero. Same bug class as #206/#212.
(Both this and 4.1 had to be fixed; either alone leaves the ramp broken.)

### 4.3 Final-step leaf dump applied to non-crop PFTs

Julia's `CNOffsetLitterfall` unconditionally set `leafc_to_litter = leafc/dt` on the final
offset step. Fortran sets it there **only inside `if (ivt(p) >= npcropmin)` — crops only**;
for natural PFTs it leaves the flux at 0, so a small leaf residual carries into dormancy.
Julia shed that residual a step early (0.19 gC/m², ~0.5 % of peak `leafc`, once a year).

### 4.4 The `use_cn` cold start ignored PFT type — **physics bug**

`cnveg_carbon_state_init_cold!` seeded **every** patch as deciduous ("Default: deciduous
non-crop"), never branching on `evergreen`/`noveg`/crop. Consequences:

- **evergreen PFTs cold-started with ZERO displayed leaf and root carbon** (all of it parked
  in storage) — an evergreen forest began with no canopy;
- bare-soil patches were given storage C they should not have;
- woody PFTs missed Fortran's `deadstemc = 0.1*ratio` seed;
- `cnveg_nitrogen_state_init_cold!` seeded **all N pools to zero**, so the cold vegetation
  had carbon but **no nitrogen** (Fortran derives N from the cold C pools via
  `leafcn`/`frootcn`/`deadwdcn`);
- `initial_vegC` was hardcoded to 20 (the CTSM *default*) with no way to set it. It is a
  **namelist** variable (`&cnvegcarbonstate`); Bow's `lnd_in` sets `initial_vegc = 100`.

All fixed against `CNVegCarbonStateType.F90` / `CNVegNitrogenStateType.F90`. Cold start is
now bit-exact vs Fortran (§3c).

### 4.5 Harness bugs found and fixed (not model bugs)

- `run_one_parity_step!` built the orbital/daylength state from the calday **before**
  `advance_timestep!`. CTSM's `get_curr_calday()` returns the calday at the **end** of the
  timestep, which is what `UpdateDaylength` uses — and what the production driver
  (`clm_run.jl`) already did correctly. The harness was one step behind, shifting the offset
  trigger by one timestep.

  Careful with the companion term: `nextsw_cday` (which drives the albedo) is the *next*
  radiation time, i.e. the END of the step. Pre-advance that was written
  `calday_start + dtime/SECSPDAY` — **numerically the same as the post-advance calday**. So
  once calday is taken post-advance, `nextsw_cday = calday` with **no further `+dtime`**,
  or the albedo is computed a whole step late (this regressed
  `test_fortran_parity_freewins.jl`'s `surface_albedo` from 2.2e-03 to 5.3e-03 before it was
  caught). With both terms consistent, daylength and albedo are now each correct.
- `fortran_parity_drift.jl` called `init_daylength!` **every** step. `update_daylength!`
  (inside `clm_drv!`, when `!is_first_step`) shifts `prev_dayl = dayl` and recomputes `dayl`
  from the same declination — so re-seeding first made `prev_dayl == dayl`, pinning
  `ws_flag = 1` and making the offset test (`ws_flag == 0 && dayl < crit_dayl`)
  **unable to fire at all** in a free run. Production never calls it in the loop; the drift
  harness now matches.

---

## 5. Open finding — NOT fixed, reported only

**Julia has no nitrogen deposition.** `atm2lnd.forc_ndep_grc` is initialised to `0.0` and
**never assigned anywhere in the codebase**, so `ndep_to_sminn ≡ 0`. CTSM reads an N-deposition
stream (`&ndepdyn_nml`). Per-step this is small (below the current ~2e-3 CN residual), which
is why it has hidden here, but over a spinup it is a systematic missing N input.

Not fixed: wiring an ndep stream reader is a new input path, not a clearly-correct one-line
fix. Flagged for follow-up.

*Consequence for the autumn/cold-start references:* the CMIP6 `fndep_*` file had been deleted
from the local inputdata tree and could not be re-fetched (the server's certificate does not
cover the inputdata hostname; disabling TLS verification was declined). The Fortran reference
runs therefore used a **synthetic zero-N-deposition stream**. Since Julia's ndep is identically
zero, this makes the comparison *cleaner* — both sides see the same (zero) deposition — but it
means the autumn window's absolute mineral-N is not directly comparable to the older summer
reference run, which had real ndep. Within-window Julia-vs-Fortran parity is unaffected.

---

## 6. Reproducing

**Fortran reference** (instrumentation restored in
`cases/symfluence_build/SourceMods/src.clm/{restFileMod,clm_driver}.F90`; the dump window is
env-gated, so a plain run emits nothing):

```bash
cd /Users/.../SYMFLUENCE_data/clm_bgc_spinup
source <case>/.env_mach_specific.sh
export DYLD_INSERT_LIBRARIES=/opt/homebrew/lib/libmimalloc.dylib   # macOS 26 xzone SIGTRAP
env PDUMP_NSTEP_LO=1759897 PDUMP_NSTEP_HI=1760376 <case>/bld/cesm.exe   # run PLAINLY, never lldb
```

**Julia harnesses:**

```bash
julia +1.12 --project=. scripts/fortran_parity_cn_autumn.jl      # autumn scorecard
julia +1.12 --project=. scripts/probe_offset_litterfall.jl       # the ramp, flux by flux
julia +1.12 --project=. scripts/fortran_parity_drift.jl autumn   # 480-step free-run drift
julia +1.12 --project=. scripts/fortran_parity_cn_coldstart.jl   # cold start
julia +1.12 --project=. scripts/fortran_parity_cn_summer.jl      # regression
```

---

## 7. What is still NOT validated

Stated plainly, so this scorecard is not over-read:

- **Only Bow.** One column, 3 patches, one PFT of each relevant phenology class. Nothing
  here generalises to other biomes, and no other site has a BGC reference.
- **Only the CN core.** **Fire, methane (`use_lch4`), isotopes (C13/C14), CNDV, MIMICS, and
  VOC/MEGAN have still never been diffed against Fortran at all.**
- **No crop BGC** (`use_crop = false` throughout). The crop branches of the offset litterfall
  and the cold start are ported but unvalidated.
- **No multi-year trajectory.** The longest free-run diff is 480 steps (20 days). Annual and
  spinup-scale conservation are not covered here.
- **N deposition is missing in Julia** (§5), so any long CN integration has a systematic N
  deficit relative to CTSM.
