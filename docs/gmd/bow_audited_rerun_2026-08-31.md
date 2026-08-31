# Audited Bow annual rerun and unchanged 207-cell gate — 2026-08-31

Status: **rerun complete; gate reproduces 204/207 exactly; the Bow BTRAN cell is
unchanged and its disposition now rests solely on independent scientific review.**

## Run

`qualification_2026-08-25.md` item 4 required regenerating the Bow annual trajectory
with the audited restart fields and namelist controls, then rerunning the unchanged
207-cell gate. Executed at commit `46ef5cb`:

```sh
DOMAIN=Bow SYMFLUENCE_DATA=/private/tmp/SYMFLUENCE_data \
CLM_OUTPUT_DIR=<results> julia --project=. scripts/parity_run_domain.jl
```

- `scripts/parity_run_domain.jl` carries the audited domain registry (the
  reference-audit fixes through `8595fa0`/`b9e558c`: CTSM fresh-snow density control,
  prior snow drainage, hourly 3600 s protocol, calibrated `baseflow`/`int_snow`).
- The **full 8760-step year completes** (213.5 s wall on Apple M5) with conservation
  checks active — consistent with C02.

## Gate

`scripts/parity_gate.py` (unchanged), `CLM_PARITY_DOMAINS=Bow,Aripuana,Stillwater`,
fresh Bow output plus the existing Aripuanã/Stillwater qualification outputs:

- **204/207 strict cells; Aripuanã 69/69; Bow 67/69; Stillwater 68/69.**
- The two documented snow-covered-area `SNOW_DEPTH` exceptions recur identically.
- The single undocumented failure is the same **Bow BTRAN +4.2% annual / 0.21 nRMSE**.

## Consequence

The audited-controls rerun neither created nor removed any cell: the BTRAN residual is
robust to the restart/namelist audit, exactly as the three-cell compiler-sensitivity
matrix showed it robust to optimization level. Everything mechanical around C01 is now
exhausted; per the standing instruction ("do not add a BTRAN exception without
independent scientific and numerical review"), C01's remaining blocker is that review
and its prespecified disposition. The `open_release_requirements` list is updated
accordingly.
