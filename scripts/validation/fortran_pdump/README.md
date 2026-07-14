# Fortran CTSM pdump instrumentation — for generating T1 parity reference dumps

The validation harness's **T1 `:parity` oracle** (`scripts/validation/validate.jl` →
`scripts/fortran_parity_common.jl`: `run_one_parity_step!` / `compare_inst_to_dump`)
diffs the Julia port's single-step state against Fortran CTSM reference snapshots
`pdump_<boundary>_n<nstep>.nc`. This directory holds the **reconstructed Fortran
instrumentation** that writes those snapshots, so new-config refs (e.g. `use_cn`)
can be generated — broadening parity beyond the existing Bow SP anchor.

## What's here

- **`pdumpMod.F90`** — a CTSM SourceMods module exposing `pdump_write(bounds, label)`.
  It writes the **16 fields the oracle actually reads** (`T_SOISNO`, `H2OSOI_LIQ/ICE`,
  `T_GRND`, `WA`, `H2OSFC`, `ZWT`, `ZWT_PERCH`, `INT_SNOW`, `SNOW_DEPTH`, `frac_sno`,
  `T_VEG`, `elai`, `tlai`, `LIQCAN`, `SNOCAN`) + `nstep`/date metadata, with the same
  dim names (`column`, `pft`, `levtot`) and the `levtot = j + nlevsno` layout as the
  legacy full-restart pdumps. Gated to an nstep window (default 13455–13470, two named
  constants). Modeled on the in-tree `dump_cv_sidecar`. Works for ANY config (SP/CN/…)
  since all 16 fields exist everywhere. It is a *minimal oracle-compatible* writer, not
  a byte-for-byte 319-var restart replica.
- **`clm_driver.hooks.diff`** — the 3 lines to add to the driver: the `use pdumpMod`
  statement + the `before_step` and `after_hydrologydrainage` call hooks (the two
  boundaries the oracle consumes). Apply to a copy of stock `src/main/clm_driver.F90`
  placed in `SourceMods/src.clm/`.

## Status (2026-06-27)

Reconstructed, validated against the dump layout + the real type defs/inst visibility,
and **compile-verified**: dropped into the instrumented Bow case's `SourceMods/src.clm/`,
`./case.build` finished successfully (`pdumpmod.mod` produced, exe relinked). The final
gated-run byte-verification did **not** complete in that session — the CTSM install tree
(`installs/clm`) disappeared from the filesystem mid-run (an environment/drive issue, not
an instrumentation fault). Re-run the recipe below once the CTSM case is restored.

## Recipe (run on the box with the instrumented CTSM build)

```bash
CASE=.../installs/clm/cases/symfluence_build      # the instrumented single-point Bow case
# 1. Drop in the instrumentation
cp pdumpMod.F90 "$CASE/SourceMods/src.clm/"
cp "$CASE/.../src/main/clm_driver.F90" "$CASE/SourceMods/src.clm/"   # stock driver
#    then apply clm_driver.hooks.diff to "$CASE/SourceMods/src.clm/clm_driver.F90"

# 2. Build — ⚠ a new SourceMods module needs the source/dependency lists regenerated,
#    else clm_driver compiles before pdumpmod.mod exists ("Cannot open module file
#    'pdumpmod.mod'"). Remove the stale lists first, then build:
OBJ="$CASE/bld/.../clm/obj"
rm -f "$OBJ/Srcfiles" "$OBJ/Depends"
cd "$CASE" && ./case.build          # ~10–15 min incremental; relinks the exe

# 3. Run into an ISOLATED RUNDIR (never clobber the validated dumps in clm_parity_run/).
#    Copy the working Bow namelists + the branch restart, shorten the run to cross the
#    dump window, repoint inputdata if relocated, then run the exe (mpiuni = no mpirun):
RUN=.../clm_parity_run/pdump_verify_run ; mkdir -p "$RUN" ; cd "$RUN"
cp .../clm_parity_run/{lnd_in,datm_in,drv_in,drv_flds_in,nuopc.runconfig,nuopc.runseq,\
   datm.streams.xml,fd.yaml,rpointer.*} .
cp .../clm_parity_run/Bow_at_Banff_lumped.{clm2,datm,cpl}.r.2003-07-16-21600.nc .
sed -i 's/stop_n = 4710/stop_n = 20/' nuopc.runconfig    # branch nstep ~13456 → window
# if inputdata was relocated, repoint it (urbantv/snicar live under installs/cesm-inputdata):
sed -i 's|/Users/.../projects/cesm-inputdata|/Users/.../installs/cesm-inputdata|g' lnd_in datm_in datm.streams.xml
source "$CASE/.env_mach_specific.sh" ; export NETCDF_PATH=/usr/local
"$CASE/bld/cesm.exe"                  # writes pdump_*_n134{55..70}.nc into $RUN

# 4. Byte-verify the new dump matches the existing reference, then wire a parity row:
#    ncdump the new pdump_after_hydrologydrainage_n13461.nc vs the existing one (16 vars).
#    For a use_cn ref, set use_cn=.true. in user_nl_clm and branch from the BGC-spun-up
#    restart (clm_bgc_spinup/.../clm2.r.2202-01-01, spinup_state=0), not the SP restart.
```

No Julia-repo changes are needed to *consume* new dumps — `oracle_parity` /
`compare_inst_to_dump` already read exactly the 16 fields `pdumpMod` emits; just add a
`vcfg("...-parity"; oracles=[:parity], ...)` row pointing at the new nstep.

---

## Flux references: the pdumps CANNOT score fluxes

`restFile_write_dump` writes a **restart-format** snapshot, so it carries only
*prognostic state*. Fluxes (`qflx_surf`, `qflx_infl`, `fsat`, `frac_h2osfc`, …)
are simply not in it, and no dump boundary can ever produce them.

For flux parity use CTSM's **history** mechanism at per-step instantaneous
resolution instead — no `restFileMod` change needed:

```
hist_nhtfrq          = -8760, 1     ! tape 2 = every step
hist_mfilt           = 20, 400
hist_avgflag_pertape = 'A', 'I'     ! tape 2 = INSTANTANEOUS, not a time mean
hist_fincl2          = 'QOVER', 'QINFL', 'FSAT', ...
```

For fluxes with no stock `hist_addfld` registration, add one in a SourceMods copy
of the owning type. `WaterFluxBulkType.runoff_hist.diff` in this directory is the
worked example: it registers the eight surface-runoff-chain internals
(`QSATEXCS`, `QINFLEXC`, `QINFLEXCS`, `QINSOIL`, `QINSOILL`, `QTS2SFC`,
`QINH2OSFC`, `QSFCDRAIN`) that let `qflx_surf` be bisected term by term.
Consumed by `scripts/fortran_parity_qflx_surf.jl`; see `docs/SURFACE_RUNOFF_PARITY.md`.

### Three traps that cost real time here

1. **A CTSM branch run needs the `.rh0` restart-history file**, not just `.r`.
   Without it the exe dies with a bare `SIGTRAP` immediately after the
   decomposition printout and **no CLM error message** — it looks like a
   compiler/allocator crash, but it is a missing input.
2. **`case.build` requires Python ≤ 3.11** (CIME imports `distutils`). On 3.12+
   it dies inside `standard_script_setup`; if you pipe its output through `tail`
   you will see a misleading exit 0 and then silently run a **stale exe**.
3. **Do not add `hist_addfld1d` to the base `WaterFluxType.F90`.** That type is
   instantiated per water tracer and a registration there traps during init. Put
   bulk-only fields in `WaterFluxBulkType.F90`.
