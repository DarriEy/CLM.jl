# FATES Fortran-parity ground truth (cold start)

Fortran-side companion to `scripts/fates_fortran_parity.jl` (the Julia dump +
`compare_dumps` scorecard, PR #195). Goal: stand up a FATES-enabled Fortran
CTSM run that emits a schema-v1 pdump of the site→patch→cohort hierarchy at
cold start (nstep=0), so the Julia FATES port can be scored field-by-field
against Fortran ground truth.

## Contents
- `FatesParityDumpMod.F90` — Fortran instrumentation module. Walks
  `site%oldest_patch→younger` and `patch%tallest→shorter` (identical order to
  the Julia side) and writes the schema-v1 records (META/SITE/PATCH/COHORT) with
  17-significant-digit reals. Authored against the verified ctsm5.3.012 FATES
  type/API names (`FatesCohortMod`, `FatesPatchMod`, `EDTypesMod`,
  `PRTGenericMod::GetState`).
- `inject_init_coldstart.txt` — the two edits that wire the dump into
  `clmfates_interfaceMod.F90 :: init_coldstart` (fires after the cohort
  hierarchy + `wrap_update_hlmfates_dyn` are fully populated).
- `setup_case.sh` — runnable recipe: `ncgen` the FATES params, `create_newcase`,
  `case.setup`, namelist config, and SourceMods install. Every step here is
  VERIFIED to work against the local install.

## Why nstep=0 is the clean first target
FATES cold-start cohorts (`EDInitMod::init_cohorts`) are seeded one-per-PFT and
are deterministic given the FATES parameter file + allometry — independent of
atmospheric forcing. The port reads `data/fates/fates_params_default.cdl`, which
is **byte-identical** to the install's
`src/fates/parameter_files/fates_params_default.cdl` (verified). So cohort
`dbh/height/n/leafc/fnrtc/sapwc/storec/structc/vcmax/…` need no forcing
alignment to compare.

## Status reached (2026-07)
| step | result |
|------|--------|
| FATES param `.nc` from shared `.cdl` (`ncgen`) | ✅ generated (byte-identical cdl confirmed) |
| `create_newcase` `I2000Clm60Fates` `ne3_ne3_mg37` | ✅ case created |
| `case.setup` | ✅ |
| `user_nl_clm` fsurdat + fates_paramfile → staged ne3np4 data | ✅ |
| `preview_namelists` (all components) | ✅ generate cleanly |
| Julia cold-start dump (`fates_fortran_parity.jl`) | ✅ 2 sites × 1 patch × 14 cohorts |
| SourceMods instrumentation installed + injected | ✅ (`FatesParityDumpMod` + `init_coldstart` call) |
| `case.build` → `cesm.exe` | ✅ **MODEL BUILD FINISHED SUCCESSFULLY** — the CLM/FATES lib **compiles the SourceMod**; `nm cesm.exe` shows `_fatesparitydumpmod_MOD_fates_parity_dump`. (One env fix needed: the build config hardcodes a stale `netcdf/4.9.3_2` path — Homebrew now ships `netcdf/4.10.0`; see setup_case.sh.) |
| direct `cesm.exe` run (NTASKS=1, mpiuni) | ⚠️ launches + initializes, then **aborts in DATM before FATES init** (see below) |
| Fortran dump / scorecard | ❌ **blocked — no DATM forcing** |

## BLOCKER — DATM (GSWP3) atmospheric forcing is not staged
`./check_input_data` reports **743 missing input files**. The dominant, blocking
set is the GSWP3 v1 0.5° global DATM forcing that `DATM%GSWP3v1` requires:

```
Model datm missing file mesh  = .../ptclm-data/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/
                                 clmforc.GSWP3.c2011.0.5x0.5.TPQWL.SCRIP.210520_ESMFmesh.nc
Model datm missing file fileN = .../Solar/clmforc.GSWP3.c2011.0.5x0.5.Solr.YYYY-MM.nc   (× ~170)
                              +  .../Precip/...Prec.YYYY-MM.nc                           (× ~170)
                              +  .../TPHWL/...TPQWL.YYYY-MM.nc                            (× ~170)
```
- The forcing root `/Users/darri.eythorsson/projects/ptclm-data/` **does not
  exist**; GSWP3 global 0.5° forcing for 1991–2004 is hundreds of GB and is not
  downloadable/appropriate in this environment.
- Staged CLM boundary data (1.8 GB, 87 files) covers only: the `ne3np4`
  surfdata / ESMF mesh / domain, `clm50_params`, snicar, urban — **no DATM
  streams**.
- Secondary (droppable): MOSART routing `frivinp` + `r05` rof mesh — avoid by
  using an `SROF` compset (`I2000Clm60FatesRs`).

Even a 1-step cold-start run needs the DATM component to initialize its stream
mesh + first data slice, so the run cannot start without this forcing.

### Run-proven (not just predicted)
The instrumented `cesm.exe` was built and launched directly (NTASKS=1, ESMF
`mpiuni`; per the macOS-26 note it is run plainly, not under lldb). It initializes
the driver and then **aborts inside the DATM component**, before the LND/FATES
component is realized — so `init_coldstart` (and the dump) is never reached:

```
check_netcdf2  <-  __piolib_mod::pio_openfile
               <-  __dshr_stream_mod::shr_stream_getcalendar
               <-  __dshr_stream_mod::shr_stream_init_from_xml
               <-  __atm_comp_nuopc::InitializeRealize      (DATM realize)
  => MPI_ABORT (opening a nonexistent GSWP3 stream file)
```

The NUOPC driver realizes ATM (DATM) **before** LND, so there is no init-ordering
trick to emit the FATES dump ahead of the DATM stream open. Staging DATM forcing
is therefore an unavoidable prerequisite — this is the single hard wall.

## To finish once forcing is staged
1. Stage GSWP3 (or point `DATM` at any staged single-point forcing via a `1PT`
   compset + matching single-point surfdata).
2. `bash setup_case.sh` (regenerates the instrumented case).
3. `./check_input_data --download && ./case.build && ./case.submit`
   — run `cesm.exe` **directly, not under lldb** (macOS-26 xzone `brk #0x1`
   guards halt lldb but run fine plainly).
4. Score:
   ```
   julia +1.12 --project=scripts scripts/fates_fortran_parity.jl compare \
       <RUNDIR>/fates_pdump_fortran_n0.txt
   ```
