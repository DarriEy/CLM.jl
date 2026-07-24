# MIMICS soil-decomposition cascade — first real numeric Fortran parity (2 biomes)

**Status (2026-07-23):** the MIMICS (Wieder 2015, microbial-explicit) decomposition
kinetics in the CLM.jl port are **byte-exact vs a stock CTSM Fortran run at two
contrasting biomes** — all 8 `decomp_k` coefficients match to `max|rel| = 1.5e-6`
(the single-precision history floor). This is the FIRST numeric Fortran oracle for
MIMICS; prior validation (#294) was only an in-session scalar re-implementation.

## Oracle (no source changes, no rebuild)
A stock `cesm.exe` (`-bgc bgc`, MIMICS compiled in) was run per site with
`soil_decomp_method='MIMICSWieder2015'`, `use_cn=.true.`, cold start. The MIMICS
kinetics are captured **straight from native history** — CTSM registers them
(`SoilBiogeochemCarbonFluxType.F90:305`) as `default='inactive'`, so they are
requested via `hist_fincl1`:
- `K_LIT_MET K_LIT_STR K_SOM_AVL K_SOM_CHEM K_SOM_PHYS K_MIC_COP K_MIC_OLI K_CWD`
  — the 8 `decomp_k` coefficients (units 1/s, per column × levdcmp);
- `W_SCALAR`, `O_SCALAR` (moisture / anoxia scalars);
- `TSOI` (soil temperature) + `HR`/`SOILC_HR`.

Decomp-time input state (the 8 vertically-resolved C pools, `DZSOI`, `TSOI`) comes
from the existing `restFile_write_dump` `pdump_before_step_n<step>.nc` snapshots
(PDUMP_NSTEP_LO/HI env window); `ligninNratioAvg` from the restart. **`before_step`
(pre-cascade) pools are essential** — the microbial pools change ~1.7e-3 per step,
so `after_hydrologydrainage` pools inflate the residual to ~1e-3.

## Port comparison (`scripts/mimics_multisite_fortran_parity.jl`)
The port's `decomp_rates_mimics!` is fed the EXACT Fortran inputs (MIMICS params
read from `clm50_params.c260305.nc`: densdep=1.2, vint=6.6, vslope=0.0756,
tau_mod pinned 1.5). `soilpsi` is inverted from the Fortran `W_SCALAR` so the port's
internal moisture scalar equals Fortran's, isolating the MM/Arrhenius/tau/desorption
kinetics. `O_SCALAR=1` (anoxia off, `use_lch4=.false.`).

## Results — max|rel| over the thawed levels (worst of 3 steps/site)

| field       | Bow (boreal summer) | Aripuana (wet tropical) | verdict |
|-------------|:-------------------:|:-----------------------:|:-------:|
| K_LIT_MET   | 8.6e-7 | 9.2e-7 | IN-TOL |
| K_LIT_STR   | 7.8e-7 | 9.2e-7 | IN-TOL |
| K_SOM_AVL   | 8.4e-7 | 8.1e-7 | IN-TOL |
| K_SOM_CHEM  | 8.0e-7 | 8.5e-7 | IN-TOL |
| K_SOM_PHYS  | 1.4e-6 | 1.5e-6 | IN-TOL |
| K_MIC_COP   | 6.3e-8 | 6.2e-8 | IN-TOL |
| K_MIC_OLI   | 6.9e-8 | 7.2e-8 | IN-TOL |
| K_CWD       | 6.9e-8 | 5.4e-8 | IN-TOL |
| W_SCALAR    | 4.9e-16 (identity) | 3.5e-15 (identity) | IN-TOL |
| O_SCALAR    | 0 (exact) | 0 (exact) | IN-TOL |

Regimes (genuinely contrasting): Bow TSOI 279–305 K, W_SCALAR max **0.46**;
Aripuana TSOI 272–304 K, W_SCALAR max **0.74** (warmer + wetter). Non-vacuous:
8/8 pools nonzero, MM microbial kinetics fully exercised. Residual = Float32 history
precision; a formula bug would be O(1), not 1e-6 tracking the history LSB.

Microbial C:N (`cn_col`) has no direct history field, but is a function of `fmet`
(from `ligninNratioAvg`) — the same `fmet` drives `K_MIC_COP/OLI`, which match to
6e-8, so `cn_col` (port: MIC_COP=4.75, MIC_OLI=7.92) is validated by construction.

## One real finding — `days_per_year` in the CWD fragmentation constant
`K_CWD = k_frag·W·O`, `k_frag = 1/(SECSPDAY·days_per_year·tau_cwd)`. Fortran uses
`get_average_days_per_year()` = **365.2425** (`SoilBiogeochemDecompCascadeMIMICSMod.F90:898,905`).
The LIVE port hardcodes **`days_per_year=365.0`** (`src/biogeochem/cn_driver.jl:651,665,1165`).
→ the live port's `K_CWD` is **0.066 % too high** (the only pool affected). With
`days_per_year=365.2425` the harness matches K_CWD to 6.9e-8 (shown). Physically
negligible (CWD fragmentation), but a real fidelity gap: the live call site should
pass `get_average_days_per_year()`-equivalent, not a hardcoded 365.0.

## Reference run recipe (both sites)
1. Clone a single-point CN-plumbing rundir (Bow `dds_run_1/final_evaluation`);
   flip `use_cn=.true.`, `soil_decomp_method='MIMICSWieder2015'`,
   `use_nitrif_denitrif=.true.`, `paramfile=clm50_params.c260305.nc`, cold start.
2. Required CN namelist that SP configs lack: `&cnphenology min_critical_dayl_method='Constant'`,
   `&cnfire_inparm fire_method='nofire'` (empty popd/light streams → malformed-URL abort),
   `&ndepdyn_nml` pointing at a real (ZERO) ndep file + mesh.
3. `hist_dov2xy=.false. hist_nhtfrq=1`, `hist_fincl1` = the K_*/W_SCALAR/O_SCALAR/TSOI list.
4. Run `./cesm.exe` directly with `MallocNanoZone=0` (fixes the macOS-26 xzone SIGTRAP),
   `HWLOC_COMPONENTS=-opencl`, `PDUMP_NSTEP_LO/HI` around the target steps.
5. Bow (boreal) is FROZEN in Jan (W_SCALAR≡0, vacuous) — run cold-start to 2002-07-16
   (~4720 steps) for thawed summer soil. Aripuana (tropical) is warm/wet from the
   start — 48 steps suffice.

Rundirs: `/private/tmp/claude-501/mimics_bow_run`, `/private/tmp/claude-501/mimics_ari_run`.
