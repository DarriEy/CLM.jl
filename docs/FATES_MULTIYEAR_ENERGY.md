# FATES multi-year energy failure — root-cause log

Follow-on to #277 (`docs/FORTRAN_VALIDATION_BACKLOG.md` D1, PR body §4).

**Status: Defect 1 CLOSED (#280). Defect 2 root-caused against a 4-year Fortran reference and fixed (the FATES clock was never advanced). Two residuals remain open; see the last section.** This file is the running record; it is written
incrementally so the reproduction survives a lost session.

## The defect, as handed over

No multi-year FATES run survives `main`. Four configs across three sites all die
with a fatal longwave imbalance of 300–400 W/m²:

| site | config | dies | `errlon` |
|---|---|---|---|
| Aripuanã | screen `:none` | day 74 | −413 W/m² |
| Aripuanã | screen `:drop_cold_deciduous` | day 98 | −405 W/m² |
| Krycklan | baseline, unprimed | day 235 | −351 W/m² |
| Bow | `FATES_WARMSOIL=0` | day 147 | −295 W/m² |

Bisected (6 points, `712cb5c..37db702`) to `f2d805c` = **#252**, which made
`use_bedrock`'s default conditional (`use_fates ⇒ .false.`, matching CTSM
`namelist_defaults_ctsm.xml:178-181`). **#252 is correct and is not to be
reverted** — the correct default *exposed* defects the wrong one masked.

Pinning `use_bedrock=true` separates two independent defects:

| | inherited (`false`) | pinned `true` |
|---|---|---|
| at `f2d805c` (#252) | dies day 74 | 89/90 PASS, `ncoh` 69 / `npatch` 2 |
| current `main` | dies day 74 | still dies day 75, `ncoh` 363 / `npatch` 5 |

* **Defect 1** — latent long before #252, unmasked by it.
* **Defect 2** — entered in `f2d805c..37db702`; recruitment additionally runs
  away ~5× (`main` 363 cohorts vs #252's 69 under an identical flag).

The daily FATES carbon balance holds 235/235, 147/147, 75/75 right up to every
failure. Per `conservation-is-not-accuracy`, a green balance is **no evidence**
here.

## Why it went unnoticed

`scripts/fates_fortran_parity.jl` sets `use_bedrock` **explicitly** and was
immune. The two multi-year harnesses — `scripts/fates_longhorizon.jl` and
`scripts/fates_multisite_validation.jl` — **inherit** the default, so #252
silently re-specified both. The inherited-default trap, one harness over from
where #274's flag audit was pointed.

## Log

(appended as the investigation proceeds)

### 2026-07-20 — validation matrix opened (rebased on `main` @ b2fef5b, post-#281/#282)

Reproduction environment pinned:

```
SD="…/My Drive/data/SYMFLUENCE_data"   # NOT the compHydro root — holds all 3 FATES domains
```

Four configs, both arms (`before` = `main` b2fef5b, `after` = this branch), 1460 d:

| key | harness | env |
|---|---|---|
| `arip_none` | `fates_longhorizon.jl` | (defaults) |
| `arip_screened` | `fates_longhorizon.jl` | `FATES_BIOGEOG=1` |
| `bow_nowarm` | `fates_longhorizon.jl` | `FATES_WARMSOIL=0` + Bow era5 fsurdat/param/forcing |
| `krycklan_baseline` | `fates_multisite_validation.jl` | `SITE=krycklan_baseline` |

**First result:** `arip_none` on this branch reaches **day 120 clean, 8/8 PASS**, against a
`main` death at **day 74**. Per `conservation-is-not-accuracy` the accompanying
`119/119` daily carbon balance is NOT cited as evidence — the survival past the
`errlon` step is.

### `before` arm — all four configs reproduce on current `main` (b2fef5b)

1460-day horizon requested; each run dies at the `errlon` assertion and stops.

| config | dies (day) | step | patch | `errlon` [W/m²] | `ncoh` range | `npatch` max |
|---|---|---|---|---|---|---|
| `arip_none` | **74** | 3601 | 7 | −412.97 | [14, 363] | 5 |
| `arip_screened` | **98** | 4753 | 7 | −405.41 | [9, 343] | 5 |
| `bow_nowarm` | **148** | 42241 | 7 | −280.20 | [14, 470] | 5 |
| `krycklan_baseline` | **235** | 11329 | 7 | −351.36 | [14, 331] | 5 |

Three of the four match #277's table to the day. **Two deltas worth recording**, both
from `main` having moved (#278/#281/#282) since #277 measured:

* Bow is now **day 148, patch 7, −280.20** where #277 recorded day 147, patch **4**,
  −294.76. The failing patch index moved, so the Bow row is the same *defect* but not
  the same *step* — do not treat #277's Bow numbers as a current baseline.
* `krycklan_baseline` is unchanged at day 235 / −351.36, confirming #282's retirement of
  the #227 prime did not move this config (it was already the un-primed arm).

`npatch` max is **5 in every one of the four**, and the failing patch index is **7 in
every one of the four** — independent of site, PFT set, forcing, and day of death. That
is the signature the root cause predicts: the failure fires when the FATES patch
population crosses the reserved-and-active HLM slot count, not at any site-specific
physical threshold.

### `after` arm — `arip_none`: 1459/1460 days, 8/8 PASS

```
ran 1459/1460 days   carbon cold=2549 final=1.728e+04   ncoh [14, 893]   npatch max=9
births 5196 / deaths 4795   8/8 PASS, no errlon
```

against `main`'s death at day 74. **Defect 1 is closed for this config.**

**The fix is a no-op until the defect fires.** Day 60 is *identical* on both arms —
`ncoh 152`, `npatch 2`, `carbon 3396.3`, `maxdbh 5.0122`, `GPP 34678.025` — so #280
changes nothing about the model until the FATES patch population crosses the
reserved-and-active slot count. That is the strongest single piece of evidence that it
is a filter/activation fix and not a physics change.

Note the daily carbon balance held **1459/1459**. Per `conservation-is-not-accuracy`
that is NOT offered as evidence of anything; see the next section for what it is blind to.

### But the surviving run is not a physical one — `maxdbh` reaches 9.35 m

The trajectory the fix makes reachable (day / `ncoh` / `npatch` / carbon / `maxdbh` cm):

```
  60  152 2   3396.3    5.01     900  349 4   3180.7  342.67
 120  878 9   4748.4   22.76    1020  287 4   8585.4  506.55
 180  881 9   4115     88.55    1140  294 4  12127    655.33
 240  490 5   4460    219.01    1260  310 4  14679    774.57
 360  460 5   6624.3  418.35    1380  363 5  15771    866.38
 600  840 9  11332    367.32    1440  408 5  17006    935.19
```

`maxdbh` goes **1.88 cm → 935 cm** in 1459 days ≈ **234 cm/yr**. A fast tropical pioneer
does 1–5 cm/yr. This is **50–200× too fast**, and it ends at a 9.35 m diameter tree.

**#280 does not cause it.** The rate is already wrong inside the byte-identical shared
prefix: 1.8848 → 5.0122 cm over the first 60 days is **19 cm/yr**, on the step range
where both arms agree to every printed digit. #280 only makes the later trajectory
*reachable*; the growth defect was always there, hidden behind a day-74 death.

Site carbon rises **monotonically** 2549 → 17006 (`min` = the cold start, `max` = the
final day) with one die-back at ~day 900 (11988 → 3180.7). The harness's
`[PASS] carbon bounded (no blow-up/collapse)` passes on that trajectory, as does
`carbon conserved EVERY day` (1459/1459) — **neither check can see a growth rate two
orders of magnitude too fast**, because carbon that is mis-allocated to stem is still
conserved. Textbook `conservation-is-not-accuracy`.

### `after` arm — `arip_screened`: 1459/1460 days, 8/8 PASS

```
ran 1459/1460 days   carbon cold=2006 final=3.190e+04   ncoh [9, 848]   npatch max=9
births 5209 / deaths 4434   8/8 PASS, no errlon
```

against `main`'s death at day 98. Same picture as `arip_none`: the run completes, and
the same unphysical growth rides along — `maxdbh` **1247 cm at day 1320** (12.5 m),
carbon again **monotone** 2006 → 31900 with `min` = the cold start.

### `after` arm — `bow_nowarm`: 1459/1460 days, 8/8 PASS

```
ran 1459/1460 days   carbon cold=2549 final=1991  min=1711  max=3937
ncoh [14, 889]   npatch max=9   births 2625 / deaths 2211   8/8 PASS, no errlon
```

against `main`'s death at day 148.

**This config is the informative one for the growth question.** Bow with
`FATES_WARMSOIL=0` is the cold, GPP-limited site, and it is the only config whose
carbon does *not* run away: it oscillates in **[1711, 3937]** and ends *below* its cold
start, which is what a bounded quasi-steady demography looks like. `maxdbh` reaches
159 cm by day 1020 — still too fast, but an order of magnitude short of Aripuanã's
1247 cm over the same window.

So the runaway scales with productivity, not with the filter fix: the same code produces
a bounded trajectory where GPP is small and an explosive one where GPP is large. That is
consistent with too much assimilate reaching structural/diameter growth, and it is
independent of #280.

### `after` arm — `krycklan_baseline`: 1459/1460 days, 10/10 PASS

```
ran 1459/1460 days   carbon cold=2549 final=1505  min=1443  max=4140  peak=4140
ncoh [14, 900]   npatch max=9   births 3100 / deaths 2926   10/10 PASS, no errlon
```

against `main`'s death at day 235. Like Bow, carbon stays bounded.

## Result: defect 1 is CLOSED on all four configs

| config | `before` (`main` b2fef5b) | `after` (#280) |
|---|---|---|
| `arip_none` | dies **day 74**, `errlon` −412.97 | **1459/1460**, 8/8 PASS |
| `arip_screened` | dies **day 98**, `errlon` −405.41 | **1459/1460**, 8/8 PASS |
| `bow_nowarm` | dies **day 148**, `errlon` −280.20 | **1459/1460**, 8/8 PASS |
| `krycklan_baseline` | dies **day 235**, `errlon` −351.36 | **1459/1460**, 10/10 PASS |

Three sites, four configs, 1460-day horizon, **no `errlon` anywhere**. The
"no multi-year FATES run survives `main`" defect is closed.

(The 1459 rather than 1460 is the harness's end-of-run day boundary, identical in
every arm — not a failure.)

## Defect 2 is NOT closed, and #280 cannot close it

#277 separated two defects by pinning `use_bedrock=true`. Re-measured here on current
`main` with that pin, 1460 d requested:

```
dies day 75, step 3649, p=7, errlon = -411.22
ncoh range [14, 363]   npatch max 5   births 368 / deaths 19
```

That reproduces #277's number exactly: **`ncoh` 363 under the pin, against the
#252-era's 69 — the ~5× recruitment runaway, alive on current `main`.**

**#280 provably cannot touch it**, for a reason stronger than any measurement: the two
arms are *byte-identical* over the whole window in which the runaway establishes itself.
Day 60 of `arip_none` prints the same `ncoh 152`, `npatch 2`, `carbon 3396.3`,
`maxdbh 5.0122`, `GPP 34678.025` on both arms, and the runaway is already at
`ncoh 363` by day 74-75 — entirely inside that shared prefix. #280 defers a patch
activation and rebuilds a filter; it changes nothing before the patch count crosses
the slot threshold, so it cannot change how many cohorts recruit before then.

**Defect 2 remains open and unattributed.** #280 closes defect 1 only.

## Defect 2 — measured against Fortran, root-caused, fixed (2026-10-28)

Everything above compares Julia with Julia (`main` vs the #252 era). This section
adds the missing oracle: a **4-year Fortran CTSM-FATES run of the `fates_parity_1pt`
case** (Bow at Banff, start 2002-07-01, dtime 3600 s, all 14 PFTs, `use_bedrock =
.false.`; `cesm.exe` ctsm5.3.012, uninstrumented, daily h0 history), and the
matching Julia run (`scripts/fates_parity_multiyear.jl`, which reuses the
flag-audited `build()` of `scripts/fates_fortran_parity.jl`, steps a NOLEAP calendar
and rolls the forcing file each 1 January).

### Reproduce

```bash
# Fortran (≈5 min for 4 years): stage a copy of the fates_parity_1pt run dir, then
python3 scripts/validation/fates_fortran_parity/stage_multiyear.py <rundir> 4
cd <rundir> && MallocNanoZone=0 HWLOC_COMPONENTS=-opencl $SYMFLUENCE_DATA/installs/clm/bin/cesm.exe
OUT=<rundir>/fortran_daily.csv julia --project=. \
    scripts/validation/fates_fortran_parity/fates_multiyear_extract.jl <rundir>/*.clm2.h0.*.nc
# (the result is committed: scripts/validation/fates_fortran_parity/fates_multiyear_fortran_daily.csv.gz)

# Julia (≈40 min for 4 years) and the scorecard:
FATES_MY_DAYS=1460 FATES_MY_OUT=/path/julia.csv julia --project=. scripts/fates_parity_multiyear.jl
julia --project=. scripts/fates_multiyear_compare.jl /path/julia.csv
```

The Fortran history has no per-cohort dbh, so "max dbh" is compared as the largest
occupied `FATES_NPLANT_SZPF` class (edges 0, 5, 10, 15, 20, 30, …, 100 cm).

### What Fortran says the trajectory should look like

| model yr | recruitment [/m² yr] | NPP [gC/m² d] | max ncoh | woody max dbh class | grass max dbh class |
|---|---|---|---|---|---|
| 1 | 6.37 | 0.45 | 177 | 0–5 cm | 10 cm |
| 2 | 256 | 1.58 | 880 | 0–5 cm | 30 cm |
| 3 | 396 | 1.49 | 672 | 0–5 cm | **≥ 100 cm** |
| 4 | 263 | 1.74 | 282 | 0–5 cm | **≥ 100 cm** |

Two things in that table change the reading of this defect:

* **Grass "dbh" in the hundreds of cm is FATES behaviour, not a port bug.** The C4
  grass (PFT 14) enters the ≥ 100 cm class around model day 750 in *Fortran*. For
  grasses `structc = 0` and the leaf/fine-root targets saturate at
  `allom_dbh_maxheight` (20–30 cm), so once past it nearly the whole carbon
  increment goes to the dbh proxy. The woody PFTs, by contrast, never leave the
  0–5 cm class in four years. The Aripuanã "maxdbh 1.88 → 935 cm" quoted above is
  the grass proxy (the pft-13/14 cohort is the largest-dbh cohort from day 1), and
  is of the same kind as this Fortran result. It does not on its own show a defect.
* **Hundreds of cohorts and ~1 recruit /m²/day are normal from year 2 on** (Fortran
  ncoh 880, recruitment 1.2 /m²/d in summer 2003), so the 363 vs 69 cohort counts in
  the Julia-vs-Julia comparison above are not by themselves evidence of a runaway.

### Before the fix, Julia fell far *below* Fortran in year 2

| window | recruitment J / F | NPP [gC/m² d] J / F | grass max dbh J / F class |
|---|---|---|---|
| year 1 (days 1–365) | 5.1 / 6.37 (0.80×) | 0.55 / 0.45 | 9.0 cm / 10 |
| days 366–669 | **2.39 / 204 (0.012×)** | **0.17 / 1.56** | 8.9 cm / 30 |

Year 1 tracks Fortran, with per-PFT recruitment matching to ~5 % through day 120.
In spring 2003 Julia recruits almost nothing and is barely productive.

### Root cause: the FATES clock was never advanced

CTSM's FATES host calls `GetAndSetTime` (clmfates_interfaceMod.F90:3730-3786) at
init (:1774) and **at the top of every daily `dynamics_driv` (:1068)**. It sets
`hlm_day_of_year = floor(timemgr_datediff(jan01, now))` (0-based), `hlm_model_day`,
and the current date. The port called `SetFatesTime` **once**, from
`clm_fates_init!`, with `day_of_year = 1` and `reference_date = 10101`
hardcoded, and never called it again. For the whole run FATES therefore believed it
was 2 January:

* `EDPhysiologyMod.F90:1073`: `if (hlm_day_of_year == gddstart) grow_deg_days = 0`.
  With NH `gddstart = 1`, the growing-degree-day sum was **reset every day**, so after
  the first cold leaf-off the cold-deciduous leaf-on test
  (`grow_deg_days > gdd_threshold`, :1130) could never pass. The `ncdstart`
  chill-day reset (:1045) and the `hlm_day_of_year > 180` GDD guard (:1091) were
  wrong as well.
* `EDMainMod.F90:853`: the annual `trim_canopy` (`hlm_day_of_year ==
  hlm_days_per_year-1`) **never fired**.
* `EDInitMod.F90:449-454`: the cold-start leaf-on/off dates (`cleafon -
  hlm_day_of_year`, …) were offset by 180 days for a July start.

The defect is older than the `f2d805c..37db702` window. It does not explain the
Julia-vs-Julia 5× cohort count; that comparison has no Fortran baseline and is not
re-attributed here.

**Fix.** `fates_get_and_set_time!(year, mon, day, tod)` (`fates_driver_init.jl`)
mirrors `GetAndSetTime`. It uses a NOLEAP 0-based day of year, sets `model_day`
relative to the run-start reference date, and carries the year forward when the host
does not pass one. `clm_fates_init!` calls it with the start date, and
`fates_daily_dynamics_step!` calls it every day from the driver's
`(year, mon, day, secs)`. Callers that drive the daily step directly without a date
(unit-test fixtures) keep the old behaviour. Non-FATES paths are untouched: the only
driver edit is a kwarg inside the `config.use_fates` daily block.

The committed 5-day parity window
(`scripts/fates_fortran_parity.jl compare fates_pdump_fortran.txt.gz`) produces a
Julia dump that is **byte-identical** before and after the fix: cold start
7.133e-16, first divergence `npp` 4.9218e-09 at step 0. Over those five days nothing
that is dumped depends on the calendar.

### After the fix

| model yr | recruitment J / F (ratio) | NPP J / F | max ncoh J / F | woody max dbh J | grass max dbh J [class] / F class |
|---|---|---|---|---|---|
| 1 | 6.44 / 6.37 (1.01) | 0.84 / 0.45 | 859 / 177 | 1.93 cm | 20 [20] / 10 |
| 2 | 367 / 256 (1.43) | 2.32 / 1.58 | 873 / 880 | 2.81 cm | 70 [60] / 30 |
| 3 | 638 / 396 (1.61) | 2.52 / 1.49 | 856 / 672 | 3.88 cm | 241 [100] / 100 |
| 4 | 568 / 263 (2.16) | 2.56 / 1.74 | 821 / 282 | 4.97 cm | 366 [100] / 100 |

On days 366–669, year-2 recruitment went from **0.012× to 1.36×** of Fortran.
Woody growth stays in Fortran's 0–5 cm class (≤ 1.1 cm/yr). The grass proxy reaches
the ≥ 100 cm class in the same model year as Fortran.

### What is still open (residuals, not fixed here)

1. **Julia NPP runs 1.5–1.9× above Fortran from year 1**, and recruitment drifts with
   it. The first divergence is the **2002-08-09 cold leaf-off**. Fortran fires it on
   model day 40 (6 of the last 10 days with `tveg24 < 7.5 °C`). Julia's `tveg24`
   during that cold spell is **0.2–0.6 K warmer** (days 32–39: 7.53/6.12/6.73/4.05/
   6.71/8.77/7.95/7.58 °C against Fortran's 6.91/5.76/6.11/3.67/6.31/8.91/7.72/7.40),
   so it counts only 4 cold days and the cold-deciduous PFTs 3, 6, 9 and 11 keep
   their leaves into autumn. The FATES phenology code matches the Fortran line for
   line. The cause is a host canopy-temperature residual that lands on a hard
   threshold.
2. **Patch count.** From model day 62 Julia keeps one new disturbance patch per day
   (npatch 2 → 9 by day 79), while Fortran stays at 2 patches until year 2. The new
   patch's structural-biomass density is 0.00515 kg/m², just above
   `force_patchfuse_min_biomass = 0.005`, so it is no longer force-fused. Its profile
   then differs from the parent's by more than `patch_fusion_tol`, so it survives.
   `fuse_patches`, `patch_pft_size_profile` and `terminate_patches` were diffed
   against the Fortran and match. This looks like a threshold crossing driven by the
   ~1 % canopy-layering differences already documented for the 20-day window, but it
   cannot be confirmed without a per-patch Fortran dump of this window (the
   instrumented `fates_parity_1pt` build is not present on this box). It is the main
   source of the 859 vs 177 year-1 cohort count.
3. Aripuanã has no Fortran reference, and it was not re-run here.

### Regression guard

* `test/test_fates_clock.jl` (always runs) checks the clock semantics: day of year,
  `model_day`, year wrap, the `trim_canopy` gate, and that the NH GDD reset fires
  once a year rather than daily.
* `test/test_fates_multiyear_regression.jl` always self-checks the committed Fortran
  reference. When given a harness CSV (`FATES_MULTIYEAR_CSV=…`) or asked to run one
  (`FATES_MULTIYEAR_RUN=1`), it scores each model year:
  - recruitment within 3× of Fortran (the pre-fix year 2 was 0.012×);
  - woody max-dbh increment ≤ 5 cm/yr;
  - woody size class equal to Fortran's;
  - grass ≥ 60 cm when Fortran reaches ≥ 100 cm.

  Without an input CSV the scoring block reports `@test_skip`; it never passes
  vacuously.

### Follow-up: `site.lat` / `site.lon` were never set (2026-10-29)

`clmfates_interfaceMod.F90:883-884` sets `sites(s)%lat = grc%latdeg(g)` and
`sites(s)%lon = grc%londeg(g)` while the FATES sites are being set up. The port never
did this, so `site.lat` stayed NaN. `currentSite%lat > 0` (`EDPhysiologyMod.F90:1033,
:1087`) was therefore false, and the phenology at Bow (51.36°N) used the
**southern-hemisphere** branch: `gddstart = 181`, so GDD was reset on 1 July, and
`ncdstart = 120`. `clm_fates_init!` now copies the gridcell lat/lon onto each site.
The other assignments in that Fortran block were checked:

* `h_gid` was already set.
* The fixed-biogeography `pft_areafrac` was already set.
* The bc allocation and zeroing were already done.

Rerun of the 4-year Bow comparison with the clock fix and the latitude fix:

| model yr | recruitment J / F (ratio) | NPP J / F (ratio) | max ncoh J / F | woody max dbh J | grass max dbh J / F class |
|---|---|---|---|---|---|
| 1 | 6.40 / 6.37 (1.01) | 0.99 / 0.45 (2.2) | 891 / 177 | 1.93 cm | 22 / 10 |
| 2 | 126 / 256 (0.49) | 1.88 / 1.58 (1.19) | 898 / 880 | 2.79 cm | 91 / 30 |
| 3 | 159 / 396 (0.40) | 1.88 / 1.49 (1.26) | 881 / 672 | 3.95 cm | 232 / 100 |
| 4 | 180 / 263 (0.69) | 2.06 / 1.74 (1.19) | 863 / 282 | 5.13 cm | 267 / 100 |

Compared with the run that had only the clock fix:

* **Years 2–4 NPP moved toward Fortran**: the excess fell from 1.5–1.9× to 1.19–1.26×.
* **Recruitment flipped from 1.4–2.2× high to 0.4–0.7× low.**
* Year 1 is unchanged up to day 120. The missed day-40 cold leaf-off (residual 1)
  and the patch-count divergence (residual 2) are both still present.

The regression test's woody-size-class check now allows one class of slack: Julia
reaches 5.13 cm in year 4, against Fortran's < 5 cm.

### Follow-up: host calendar, and the #349 host-coupling fixes (2026-10-29)

**Calendar.** The `fates_parity_1pt` case runs `calendar = GREGORIAN`. The evidence:
`nuopc.runconfig` CLOCK_attributes; `lnd.log` reports "Calendar type: GREGORIAN";
the h0 file has 1462 daily records for 2002-07-01..2006-07-01.

The earlier harness stepped NOLEAP while reading the forcing at real dates, so
every Julia day from 2004-03-01 on was one day off the reference. Two changes fix it:

* `FATES_HLM_CALENDAR` (`:noleap` default, `:gregorian`) now drives
  `fates_get_and_set_time!`, so `days_per_year` is 366 in 2004, as CTSM's time
  manager has it.
* The harness takes `FATES_MY_CALENDAR` (default `gregorian`). Under `noleap` it
  maps forcing records by elapsed NO_LEAP time, which is what DATM does (#347).

**#349.** Commit 163a4ab was cherry-picked. It fixes four host-coupling defects:

* FATES `fcansno` was passed as 0;
* FATES `z0m`/`displa` were overwritten every step;
* `use_clm5_fpi` was not set in `clm_initialize!`;
* `is_fates` was flagged after the cold start.

`use_clm5_fpi = true` now also applies to non-FATES callers of `clm_initialize!`,
which is #349's change and intended.

Julia / Fortran ratios, years 1–4 (Gregorian, 1461 days):

| run | recruitment | NPP | year-1 max ncoh J/F |
|---|---|---|---|
| (a) leap fix only | 1.01 / 0.49 / 0.40 / 0.70 | 2.21 / 1.18 / 1.27 / 1.18 | 891 / 177 |
| (b) leap fix + #349 | 1.17 / 1.42 / 1.66 / 2.20 | 1.03 / 1.34 / 1.59 / 1.40 | **195 / 177** |

Run (a) is the same as the NOLEAP run to within 1 %, so the leap-day misalignment
was not material.

With #349 applied:

* **Year-1 NPP matches Fortran** (0.460 vs 0.448 gC/m²/d).
* **Residual 2 (patch count) is gone.** Julia now holds 2 patches until spring 2003
  (9 by day 381), as Fortran does.
* **The cold leaf-off still does not match.** Julia now fires it on **model day 11**;
  Fortran fires it on day 40. Over days 2–7 Julia's `tveg24` is 2.33 / 4.46 / 6.78 /
  7.44 / 5.09 / **7.39 °C**, against Fortran's daily TVEG of 2.61 / 4.68 / 6.78 /
  7.38 / 5.01 / **7.90**. That is 6 days below 7.5 °C against Fortran's 5, which trips
  `ncolddays > 5` once the 10-day memory fills. The canopy is now within a few tenths
  of a kelvin of Fortran, but it sits on a hard threshold, so the leaf-off date
  remains sensitive.

From year 2 on, recruitment runs 1.4–2.2× high and NPP 1.3–1.6× high. The regression
test passes on run (b).

### Follow-up: last year's active layer never reached FATES (2026-10-29)

**Where years 2+ first diverge.** Per PFT, the year-2 excess is in the hard
drought-deciduous grasses (PFTs 13 and 14). They switch status, and with it
germination (`SeedGermination`, EDPhysiologyMod.F90:2445-2450), at the wrong time:

| event (PFT 14 recruitment on/off) | Fortran day | Julia day, before | Julia day, after |
|---|---|---|---|
| autumn 2002 off | 123 | 121 | 121 |
| spring 2003 on | 254 | 265 | 265 |
| **autumn 2003 off** | **555** (6 Jan 2004) | **497** (9 Nov 2003) | **550** |
| **spring 2004 on** | **650** | **600** (`prolonged_off`, day 335 of the window) | **651** |
| autumn 2005 off | 1259 | 1252 | 1254 |

The drought status reads a 10-day soil-matric-potential memory. That memory uses
layers `2:nlevroot` with `nlevroot = max(2, min(nlevsoil,
max_rooting_depth_index_col))` (EDPhysiologyMod.F90:1207-1245), and frozen layers
count as `smp_lwr_bound = -1e6 mm`. Julia went "moisture off" as soon as the top
soil layers froze.

**Root cause.** `max_rooting_depth_index_col` is
`min(nlevsoil, altmax_lastyear_indx_col)` (clmfates_interfaceMod.F90:1131-1132).
`altmax_lastyear_indx_col` is rolled from `altmax_indx_col` on the first step of
1 January (`mon == 1 .and. day == 1 .and. sec/dtime == 1`,
ActiveLayerMod.F90:114). `clm_drv!` called `alt_calc!` with a constant
`sec = 0`, so the roll never fired. The index kept its cold-start 0
(ActiveLayerMod.F90:305) for the whole run, and every FATES root profile
(drought phenology, btran, root litter) stayed in the top 2 layers in every year.
In year 1 this is also what Fortran does, which is why year 1 matched. At Bow the
roll gives index 11 from 1 January 2003.

**Fix.** Under `use_fates`, `clm_drv!` now passes the real `secs` to `alt_calc!`.
Non-FATES paths keep `sec = 0` and stay byte-identical. CN's
`soil_bgc_vertical_profile!` has the same missed roll; that is left for a separate
change.

Also fixed in this branch: the integration merge had left two copies of
`fates_get_and_set_time!`. The later, NOLEAP-only copy overrode the Gregorian one
and broke precompilation. This was dropped upstream in 99b7b51, and the branch is
rebased on it.

**Result** (Gregorian, 1461 d; "before" = this tree with `altmax_lastyear_indx_col`
forced to 0 every step, which reproduces run (b) above):

| model yr | recruitment J/F before → after | NPP J/F before → after | max ncoh J / F after |
|---|---|---|---|
| 1 | 1.17 → 1.16 | 1.03 → 1.00 | 197 / 177 |
| 2 | 1.43 → 1.57 | 1.35 → 1.30 | 874 / 880 |
| 3 | 1.74 → 1.53 | 1.59 → 1.42 | 852 / 672 |
| 4 | 2.00 → 2.16 | 1.40 → 1.40 | 822 / 282 |

The phenology timing now tracks Fortran to within 1-5 days in years 2 and 4. Through
winter 2003-04 the cohort count and grass `nplant` follow Fortran closely
(day 541: ncoh 643 / 623; before: 481 / 623). **The yearly recruitment ratio does not
close, so the remaining excess has a different cause:**

* **Seed bank.** Recruits per unit seed bank are now slightly *below* Fortran in
  years 2-4 (10.8 vs 12.2 /yr per kg C m⁻²). The excess comes from a seed bank
  (`FATES_SEED_BANK` = seed + seed_germ) that is 1.8× Fortran by autumn 2003 and
  2.2-2.4× by year 4. Year-2 NPP is 1.30× Fortran, and the excess starts in spring
  2003: grass leaf carbon is +14 % by late May. Seed production from mature grass
  cohorts (`seed_alloc_mature` 0.2-0.25) is the likely amplifier. Attributing it
  needs per-PFT NPP and repro output from both sides.
* **C4 grass population.** Julia's PFT-14 `nplant` collapses each summer
  (day 721: 0.15 vs 4.5 /m²) while its NPP stays high, giving fewer, larger plants
  (grass max dbh 71 cm vs Fortran's 30 cm class in year 2). This is a grass
  mortality/competition residual.
* **Year 1.** The seed bank is 0.72× Fortran from day 11, the day of Julia's early
  cold leaf-off (Fortran: day 40; the host `tveg24` threshold residual above). That
  removes seed input from the cold-deciduous woody PFTs.
* The year-3 off/on (Fortran: lifespan `timeoff` on day 1016, re-flush on day 1047;
  Julia: moisture-off on day 930, `timeon` on day 1031) still differs, because
  winter 2004-05 soil liquid water in the root zone differs.

Regression guard: `test/test_fates_rooting_depth_roll.jl` (isolated subprocess).
(A) checks the `alt_calc!` roll semantics. (B) is data-gated: a live FATES
`clm_drv!` step at 1 January 01:00 rolls the index, and the daily bc pack hands it
to FATES as `max_rooting_depth_index_col`.

### Follow-up: the melt-season canopy closes the year-2+ excess (2026-10-29)

The remaining year-2..4 excess (recruitment 1.5-2.2x, NPP 1.3-1.4x, seed bank
1.8-2.4x, PFT-14 `nplant` collapse) was the host spring canopy. Three causes, found on
a per-patch hourly Fortran tape of the SPITFIRE case (README "Melt season: snow
burial, per-patch tveg24, case snow namelist"):

1. the host snowpack never reached FATES (`wrap_update_hlmfates_dyn`,
   clmfates_interfaceMod.F90:1494-1504 -> `UpdateFatesAvgSnowDepth`), so the canopy
   was never buried, `elai = tlai` all winter and the leaf absorbed 2-15x the SW;
2. every patch's `tveg24` was driven from patch 1's photosynthesis-iteration
   `t_veg_pa` instead of its own `t_veg_patch` (`WrapUpdateFatesRmean`,
   clmfates_interfaceMod.F90:3227-3245);
3. the FATES harnesses ran the port's Anderson1976 / no-wind snow defaults, while every
   Fortran FATES lnd_in has Vionnet2012 + `wind_dependent_snow_density`.

No FATES demography code changed. The multiyear harness now also writes per-PFT NPP,
leaf C, storage C, reproductive C and seed (`npp_*`, `leafc_*`, `storec_*`, `repro_*`,
`seed_*`), comparable with `FATES_NPP_PF` / `FATES_LEAFC_PF` / `FATES_STOREC_PF`
of a Fortran rerun with those fields added (`stage_multiyear.py` field list +
`FATES_LEAFC_PF, FATES_STOREC_PF, ...`).

4-year Gregorian Bow run, J / F (before = #357, `julia_fix.csv`):

| model yr | recruitment before -> after | NPP before -> after | max ncoh J / F after |
|---|---|---|---|
| 1 | 1.16 -> **1.00** | 1.00 -> **0.99** | 178 / 177 |
| 2 | 1.57 -> **0.90** | 1.30 -> **0.98** | 879 / 880 |
| 3 | 1.53 -> **0.99** | 1.42 -> **0.99** | 757 / 672 |
| 4 | 2.16 -> **1.01** | 1.40 -> **0.98** | 282 / 282 |

* Seed bank J/F by year: 1.00 / 0.94 / 0.99 / 1.02 (before 0.72 in year 1, 1.8-2.4x later).
* Grass NPP per PFT (13, 14), J/F by year: 1.006/1.006, 1.000/1.005, 0.999/0.999,
  0.998/0.972. Grass leaf C days 240-330 (spring 2003): 0.993 (PFT 13), 0.980 (PFT 14),
  before +14 %.
* PFT-14 `nplant` [m-2], J / F: day 700 32.2 / 35.3, day 721 8.30 / 9.21 (before 0.74),
  day 760 4.09 / 4.55 (before 0.12), day 1100 4.05 / 3.93. A summer drop is Fortran
  behaviour too (in Fortran it is almost all `FATES_MORTALITY_TERMINATION_SZPF`,
  days 700-760); the port's deeper collapse went away with the canopy fixes, with no
  change to the mortality code.
* Year 1 now matches: the cold-deciduous leaf-off fires on model day 39 (Fortran 40,
  residual 1 above; before: day 11).

Still open: year-2 recruitment 0.90 and PFT-12 (arctic C3 grass, a trace PFT,
1e-12 kgC m-2 s-1 NPP) 1.2x in year 4.
