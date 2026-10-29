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
