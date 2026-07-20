# FATES multi-year energy failure — root-cause log

Follow-on to #277 (`docs/FORTRAN_VALIDATION_BACKLOG.md` D1, PR body §4).

**Status: IN PROGRESS.** This file is the running record; it is written
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
