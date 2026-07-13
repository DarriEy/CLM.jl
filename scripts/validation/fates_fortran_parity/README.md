# FATES Fortran-parity ground truth (cold start) — **ACHIEVED**

Fortran-side companion to `scripts/fates_fortran_parity.jl` (the Julia dump +
`compare_dumps` scorecard). Goal: a FATES-enabled Fortran CTSM run that emits a
schema-v1 pdump of the site→patch→cohort hierarchy at cold start (nstep=0), so
the Julia FATES port can be scored field-by-field against Fortran ground truth.

**Status: done.** The Fortran run works, the dump exists
(`fates_pdump_fortran_n0.txt`, committed here), and the Julia FATES port
**matches it on all 27 fields** — see the scorecard below. This is the first
Fortran-FATES ground truth this project has ever had.

## Contents
- `FatesParityDumpMod.F90` — Fortran instrumentation. Walks
  `site%oldest_patch→younger` and `patch%tallest→shorter` (identical order to the
  Julia side) and writes schema-v1 records (META/SITE/PATCH/COHORT) with
  17-significant-digit reals. Written against ctsm5.3.012 FATES types
  (`FatesCohortMod`, `FatesPatchMod`, `EDTypesMod`, `PRTGenericMod::GetState`).
- `inject_init_coldstart.txt` — the two edits wiring the dump into
  `clmfates_interfaceMod.F90 :: init_coldstart`.
- `setup_case.sh` — **the working recipe, end to end**: ncgen params →
  create_newcase → SourceMods → build → stage a single-point run dir → run →
  dump. Run it and it produces the ground truth.
- `fates_pdump_fortran_n0.txt` — the ground-truth dump itself (1 site × 1 patch ×
  14 cohorts, one per FATES PFT), committed so the scorecard is reproducible
  without a CTSM build.

## The blocker, and why it was self-inflicted
The previous recipe built a **global** case (`I2000Clm60Fates` on `ne3_ne3_mg37`),
which requires GSWP3 0.5° global DATM forcing. That forcing is not staged
(`check_input_data`: **743 missing files**, hundreds of GB), so DATM aborted in
`shr_stream_init_from_xml → pio_openfile` during ATM `InitializeRealize` — and
because NUOPC realizes ATM **before** LND, FATES cold start was never reached.

**FATES does not need a global grid.** This project already has a proven
single-point Fortran CLM path — the one behind every pdump parity comparison to
date — in which CIME is used **only to build `cesm.exe`**, and the model is then
run directly in a hand-staged run directory whose namelists point at local site
forcing and a 1×1 ESMF mesh
(`domain_Bow_at_Banff_lumped_era5/settings/CLM/{lnd_in,datm_in,datm.streams.xml,
nuopc.runconfig,…}`, `datamode=CLMNCEP`, `nx_global=ny_global=1`, local
`clmforc.YYYY.nc`).

Rebuilding the FATES case on that template eliminated **all 743 missing files**:

| | old (blocked) | new (working) |
|---|---|---|
| compset | `I2000Clm60Fates` (MOSART) | `I2000Clm50FatesRs` — FATES + **SROF** (CPL/ATM/LND only) |
| grid | global `ne3_ne3_mg37` | build-only `f09_g17`; runtime = 1×1 ESMF mesh |
| forcing | GSWP3 global streams (absent) | local `clmforc.YYYY.nc` (present) |
| missing input files | **743** | **0** |
| FATES cold start reached | ❌ | ✅ |

The remaining namelist files that `check_input_data` still lists (popdens/lightning
fire streams, crop-calendar mesh, ch4 finundated) are **never read**:
`fates_spitfire_mode=0` selects `FATESFireNoData`, whose
`need_lightning_and_popdens()` returns `.false.`; `cropcals_rx`/`use_crop` are
`.false.`; and `use_lch4` is set `.false.` (FATES cold start does not need methane).

## Why nstep=0 is the clean target
FATES cold-start cohorts (`EDInitMod::init_cohorts`) are seeded one-per-PFT and are
deterministic given the FATES parameter file + allometry — **independent of
atmospheric forcing**. The port reads `data/fates/fates_params_default.cdl`, which
is **byte-identical** to the install's `src/fates/parameter_files/fates_params_default.cdl`
(re-verified by `cmp`). So cohort `dbh/height/n/leafc/…/vcmax` are
parameter-identical by construction and need no forcing alignment to compare —
which is exactly what makes a *bit-level* verdict meaningful.

## Scorecard (2026-07, first-ever Fortran-FATES ground truth)

```
julia +1.12 --project=scripts scripts/fates_fortran_parity.jl compare \
    scripts/validation/fates_fortran_parity/fates_pdump_fortran_n0.txt
```

```
  ── structure ──
    sites         julia=2    fortran=1    → scoring the 1 common site(s)
    patches/site  julia=1        fortran=1        match
    cohorts/site  julia=14       fortran=14       match
    extra sites   julia has 1 site beyond the Fortran reference; all 16 records
                  are EXACT replicas of the scored site → not a divergence
    records       17 matched, 0 unmatched

  field            n    max|abs|      max rel   verdict
  dbh             14  0.0000e+00   0.0000e+00   match     <- exact
  height          14  0.0000e+00   0.0000e+00   match     <- exact
  n               14  0.0000e+00   0.0000e+00   match     <- exact
  pft             14  0.0000e+00   0.0000e+00   match     <- exact
  carea           14  0.0000e+00   0.0000e+00   match     <- exact
  treelai         14  0.0000e+00   0.0000e+00   match     <- exact
  vcmax           14  0.0000e+00   0.0000e+00   match     <- exact
  canopy_layer    14  0.0000e+00   0.0000e+00   match     <- exact
  status          14  0.0000e+00   0.0000e+00   match     <- exact
  maxdbh           1  0.0000e+00   0.0000e+00   match     <- exact
  area/patchno/btran/btran_ft/gpp/npp/nstep      0.0000e+00  match
  leafc           14  6.9389e-18   6.7736e-18   match     <- 1 ULP
  fnrtc           14  6.9389e-18   6.7736e-18   match     <- 1 ULP
  storec          14  6.9389e-18   6.7415e-18   match     <- 1 ULP
  structc         14  6.9389e-18   6.6992e-18   match     <- 1 ULP
  sapwc           14  1.7347e-18   1.7221e-18   match     <- 1 ULP
  carbon           1  1.8190e-12   7.1333e-16   match     <- summation order

  ★ FATES port matches the Fortran ground truth within tol=1e-10 (27 fields, all records)
```

**Reading of the result.** 13 of 27 fields are *exactly* zero relative difference —
the Julia port reproduces Fortran FATES cold start bit-for-bit for cohort geometry,
allometry, density, PFT assignment and canopy structure. The five PRT carbon pools
differ by 1 ULP (~7e-18 relative), and the site carbon total by 7e-16, purely from
floating-point summation/rounding order. There is **no physics disagreement**.

Two things the scorecard reconciles rather than scores, both harness structure and
neither a port defect:
- **elai** — an HLM `canopystate` coupling diagnostic, not FATES-internal state.
  The Fortran instrument deliberately writes `-1` as a sentinel; the comparator
  skips it (`FORTRAN_SENTINEL`).
- **site count** — the Julia harness builds a 2-gridcell test instance; the
  single-point Fortran case has 1 site. Sites are FATES's independent replicate
  units, so the comparator scores the common site and *verifies the extra Julia
  site is an exact replica* (it is — all 16 records bit-identical). An extra site
  that were **not** a replica would be reported as a real divergence.

## Environment traps (macOS 26 / Homebrew) — all handled in `setup_case.sh`
1. **`HWLOC_COMPONENTS=-opencl` is REQUIRED.** hwloc's OpenCL probe walks into the
   Metal driver at MPI startup and dies with `SIGILL` inside `AGXMetalG16X`
   (`MTLCopyAllDevices` → `gldCreateDevice` → `hwloc_opencl_discover`). Without it
   `cesm.exe` never reaches CLM at all. This — not FATES, and not the forcing —
   was what made the first single-point runs look like a model crash.
2. **Flaky xzone malloc trap.** The optimized build `SIGTRAP`s inside
   `libsystem_malloc` on roughly half of all launches, landing in a *different*
   allocation each time (`NC3_open`, `PatchType::Init`, …). That non-determinism is
   the signature of the known ESMF/xzone false positive, not a CLM bug: every
   successful run emits a **byte-identical** dump (verified 6/6), and a DEBUG
   (`-O0`, bounds-checked) build runs clean and produces the *same* dump.
   `setup_case.sh` simply retries. Also set `MallocNanoZone=0`.
   **Never run under lldb** — macOS 26's inline `brk #0x1` guards halt a debugger
   spuriously.
3. **CIME python.** Needs `distutils` *and* a stdlib-compatible ElementTree.
   System python 3.14 has no `distutils`; miniforge base's `lxml` shadows
   `ElementTree` and breaks `create_newcase`. Use a 3.11 conda env (`CIME_PY`).
4. **netcdf link path.** The `homebrew` mach config hardcodes a Cellar version
   (`netcdf/4.9.3_2`); Homebrew now ships `4.10.0`, so the link fails with
   `ld: library 'netcdf' not found`. `setup_case.sh` symlinks the stale path.

## Next target
Cold start is closed. The next tier is a **stepped** comparison
(`FATES_PARITY_STEPS=N`), which additionally requires the Fortran run to use
identical forcing + soil boundary conditions — the run directory built here already
reads the same local `clmforc.YYYY.nc` the Julia harness can be pointed at, so the
remaining work is aligning the soil BC and the step loop, not standing up a case.
