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
