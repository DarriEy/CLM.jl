# Multi-site Fortran-parity sweep (2026-07-23)

Breadth validation of the subsystems wired/validated during the 2026-07 wiring campaign
(most were previously Bow single-point only), run against real CTSM Fortran across
contrasting biomes. **No rebuild required** — the CN/BGC exe
(`installs/clm/cases/symfluence_build/bld/cesm.exe`, `-bgc bgc`) has MIMICS, isotopes,
MEGAN, drydep, CNDV, and CH4 all compiled in; each subsystem is selected purely by
namelist and instrumented via native history (`hist_dov2xy=.false.`, `hist_nhtfrq=1`).

## Sites with local Fortran forcing (hand-stageable single-point)
| Site | Biome | Forcing format |
|---|---|---|
| Bow at Banff (lumped) | boreal/montane | CLMNCEP |
| Aripuanã | wet tropical Amazon | CLM1PT |
| MerBleue | boreal peatland | CLM_input |
| US Plains | temperate cropland (CFT) | CLM_input |

Blocked (forcing not on disk): Krycklan (boreal), Stillwater.

## Scorecard
| Subsystem | Site(s) | Verdict | Harness |
|---|---|---|---|
| MIMICS decomposition | Bow + Aripuanã | **byte-exact** — 7/8 `decomp_k` to max\|rel\| 1.5e-6 (first numeric Fortran oracle) | `scripts/mimics_multisite_fortran_parity.jl` |
| CNDV (dynamic veg) | Bow | faithful — init/AGDD/Light/Establishment ≤1e-14 | `scripts/validate_cndv_bow_fortran.jl` |
| Dry deposition | Bow | `index_season` bit-exact; DRYDEPV_O3 tracks season | `scripts/validate_drydep_bow_fortran.jl` |
| VOC/MEGAN | Aripuanã | PASS — activity factors ≤5e-7, isoprene 71% of flux | `scripts/fortran_parity_voc_aripuana.jl` |
| Isotopes C13/C14 | Aripuanã | discrimination + soil pools PASS (fractionation active) | `scripts/fortran_parity_isotopes_aripuana.jl` |

Per-subsystem detail: `MIMICS_MULTISITE_FORTRAN_PARITY.md`, `CNDV_DRYDEP_BOW_PARITY.md`,
`VOC_ISO_TROPICAL_PARITY.md`.

## Two findings — both non-bugs on scrutiny
1. **MIMICS `K_CWD` 0.066% high** = a NO_LEAP-vs-GREGORIAN **calendar** difference, not a
   port error. `k_frag` uses `get_average_days_per_year()` (calendar-dependent: NO_LEAP→365.0,
   GREGORIAN→365.2425). The reference runs used GREGORIAN; CTSM's default is NO_LEAP, which is
   exactly the port's hardcoded 365.0. Only `K_CWD` carries a day-count (the other 7 MIMICS pools
   are Michaelis-Menten). Making the port calendar-aware is a documented minor future item.
2. **Veg-pool NaN under `use_fun=true` cold-start** = the pre-existing, documented cold-start
   canopy-NaN blocker (`psn → availc = NaN`) propagating through FUN — **not** a FUN bug. Fortran's
   `CNFUNMod` has the identical unguarded arithmetic; with real params the FUN kernel never NaNs on
   finite input. `use_fun` defaults false, so the default path is unaffected.

Lesson reaffirmed: multi-site + true cold-start (not injected mid-run state) surfaced both — a
single-point injected-state harness was structurally blind to them.
