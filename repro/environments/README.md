# Reproduction environments

The development package deliberately does not commit a root `Manifest.toml`; CI resolves
the declared Julia compatibility range. The paper artifact will additionally contain a
locked CPU environment and its checksum so the exact submitted dependency graph remains
reconstructable.

Capture a machine-readable record with:

```sh
julia --project=. scripts/gmd/record_environment.jl > environment.json
```

Generated records contain host paths and a hostname. Review them before publication and
redact only privacy-sensitive identifiers, never software versions or scientific
provenance. A dirty Git tree is recorded explicitly and is not acceptable for a definitive
experiment.

## Bow BTRAN Linux sensitivity cell

The cell runner intentionally does not create a CIME machine definition or acquire model
inputs. Those are controlled, separately audited artifacts. On a Linux host with an
already-configured case and isolated run directory:

```sh
export CTSM_CASE=/archive/cases/bow-release
export CTSM_ROOT=/archive/src/ctsm5.3.012
export CTSM_RUN_DIR=/scratch/bow-release-run
export CTSM_LAUNCH=/scratch/bow-release-build/cesm.exe
export SYMFLUENCE_DATA=/archive/SYMFLUENCE_data
export BTRAN_CELL_ID=linux-gnu-release
scripts/gmd/run_btran_sensitivity_cell.sh
```

Run the no-optimization cell from a distinct case whose recorded `Macros.make` uses
`FFLAGS_NOOPT=-O0`; set `BTRAN_CELL_ID=linux-gnu-o0`. Never mutate and reuse the release
case in place, because that makes the executable/configuration association ambiguous.

The runner requires Linux, explicit staged paths, and a clean CLM.jl worktree. It builds
the supplied case, executes CTSM, runs the step-2901 Julia comparator, and writes
`result.json` plus complete logs and an environment record. Review the generated record
before moving it into the experiment archive.

After all cells finish, aggregate them without ranking or filtering:

```sh
python3 scripts/gmd/aggregate_btran_sensitivity.py \
  /archive/results/linux-gnu-release/result.json \
  /archive/results/linux-gnu-o0/result.json \
  /archive/results/macos-gnu-release/result.json \
  --output /archive/results/matrix.json
```

The aggregate records every cell's error and the full min/max range. It deliberately
cannot declare acceptance or select a best-matching compiler configuration.
