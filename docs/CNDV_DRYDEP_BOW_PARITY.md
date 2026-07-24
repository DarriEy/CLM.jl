# CNDV + dry-deposition Fortran parity at Bow-at-Banff (2026-07-23)

Fortran ground truth: `symfluence_build/bld/cesm.exe` (`-bgc bgc`, has CNDV + drydep
compiled). Bow single-point (1 col / 3 phys PFTs: bare itype0, needleleaf tree
itype1, C3 grass itype12). Runs staged in `/private/tmp/claude-501/{drydep_run,cndv_run}`
(shared build/case left PRISTINE — copied rundirs, no rebuild, no `--clean`).

## Part B — dry deposition (`DryDepVelocity.F90` / port `dry_dep_velocity.jl`, #298/#300)

Run: branch from Bow 2006-01-01 restart, full year 2006, daily native history
(`hist_dov2xy=.false.`), `drv_flds_in &drydep_inparm drydep_list='O3','SO2','NO2','CO','H2O2'`
+ `dep_data_c20221208.nc`. `DRYDEPV_<spec>` (cm/s) are history fields → no rebuild.
`SUCCESSFUL TERMINATION`. Harness `scripts/validate_drydep_bow_fortran.jl`.

`index_season` is a LOCAL Fortran var (not output) but a DETERMINISTIC function of the
published daily state (elai, snow_depth, annlai min/max, mlaidiff, landunit) —
`DryDepVelocity.F90:384-418`. Reconstructed each day two ways: (a) the port
`_wesely_index_season` (#300 code under test), (b) an INDEPENDENT inline transcription
of the Fortran rule — both fed Fortran's own daily ELAI/SNOW_DEPTH + surfdata annlai.

| field | scope | max abs | max rel | verdict |
|-------|-------|---------|---------|---------|
| index_season (port fn vs independent Fortran-rule oracle) | tree + grass, 365 d each | 0 (0/365 days) | 0 | PASS — bit-faithful transcription |
| index_season seasonal variation (boreal) | tree {1,4}, grass {1,2,4} | — | — | PASS — 3 distinct seasons, not stuck |
| DRYDEPV_O3 (Fortran ground truth, non-vacuous) | tree | 0.0116–0.485 cm/s | — | PASS — nonzero, tracks season |
| DRYDEPV_SO2 (Fortran ground truth) | tree/grass | means 0.15 / 0.08 cm/s | — | PASS — nonzero |

Season transitions are physically correct and independently corroborated by the Fortran
DRYDEPV: winter (season 4, snow_depth>0) DRYDEPV_O3 ≈ 0.09 cm/s → snow-free lush summer
(season 1, elai>0.5·maxlai) ≈ 0.24 cm/s → snow returns → season 4. The evergreen tree
(minlai 2.5) never hits dormant season 3; the grass (minlai 0) shows the autumn shoulder
(season 2, mlaidiff>0) for 3 days. All faithful to `DryDepVelocity.F90:398-416`.

NOTE on velocity magnitude parity: a per-day numeric DRYDEPV diff would need Fortran's
instantaneous `ram1`/`fv` (not standard history fields), so it is not teacher-forced here.
The port-governed part of the velocity — `index_season` → the Wesely `rc` tables — is
bit-faithful, and the Fortran DRYDEPV values are the ground-truth range/seasonality check.

## Part A — CNDV (`CNDVDriverMod`/`CNDVType.F90` / port `cndv.jl`)

CNDV updates vegetation ANNUALLY (year-end). Run: cold start (`finidat=''`),
`use_cn=.true. use_cndv=.true.`, `suplnitro='ALL'`, `soil_decomp_method='CENTURYKoven2013'`,
`fire_method='nofire'`, ndep = ZERO stream (nn-regridded, N supplemental so value inert),
`min_critical_dayl_method='Constant'`; 1 model year (2002) at Bow. `SUCCESSFUL TERMINATION`;
wrote `clm2.r.2003-01-01` with the post-year-1 CNDV state. Harness
`scripts/validate_cndv_bow_fortran.jl` (+ kernel oracles in `validate_cndv_fortran_parity.jl`).

| field | port vs Fortran | max abs | max rel | verdict |
|-------|-----------------|---------|---------|---------|
| cold-start init: present/crownarea/nind/agdd20 | dgvs_init_cold! vs CNDVType InitCold | 0 | 0 | PASS — all 0 |
| cold-start init: tmomin20 (=TFRZ−5) | 268.15 vs SHR_CONST_TKFRZ−5 | 0 | 0 | PASS |
| AGDD accumulator formula | port kernel vs CNDVType.F90:514 | 0 (identical) | 0 | PASS |
| AGDD year-1 trajectory (0→506.10 K) | Fortran history vs restart AGDD_VALUE | ~0.004 K | 2e-3 | PASS — non-vacuous |
| Light / Establishment / prec365 kernels | port vs independent oracles | ≤1e-14 | ≤1e-10 | PASS (20/20) |

Fortran year-1 cold-start outcome (the honest horizon): from bare ground the first
annual CNDV update establishes ONLY the C3 grass — restart `present[itype12]=1`,
`nind=1`, `fpcgrid=1.67e-4`; bare `fpcgrid=0.99983`; all trees still absent. AGDD reaches
506.10 K. Non-trivial multi-PFT vegetation needs a multi-year spin-up (many annual
updates); a 1-year anchor cannot exercise the developed state, and no such dump was
fabricated. The establishment CODE PATH that produced grass-only (grass→nind=1, fpc
recompute, bare fills remainder) is the exact path value-checked to 1e-10 in the oracle
harness.

## Reproduce
- drydep: `julia +1.12 --project=. scripts/validate_drydep_bow_fortran.jl`  (needs the 2006 h0)
- cndv:   `julia +1.12 --project=. scripts/validate_cndv_bow_fortran.jl`     (needs the 2002 h0)
- cndv kernels: `julia +1.12 --project=. scripts/validate_cndv_fortran_parity.jl`  (self-contained, 20/20)
