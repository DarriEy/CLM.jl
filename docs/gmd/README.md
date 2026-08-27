# GMD submission workspace

This directory is the source of truth for preparing the first CLM.jl submission to
*Geoscientific Model Development* (GMD). Historical manuscript drafts and figures are
not inputs to the submission: every reported number, table, and figure must be regenerated
from the frozen release candidate and experiment protocol.

## Proposed paper

Working title:

> **CLM.jl v0.1.0: verification-gated translation of CLM5 biogeophysics into
> differentiable, accelerator-capable Julia**

Primary claim boundary: the configured single-site/single-column CLM5 biogeophysics path,
validated against independently generated Fortran reference results, plus automatic-
differentiation and accelerator demonstrations using that validated path.

Biogeochemistry, methane, fire, FATES, every CTSM configuration, and production/global
readiness are outside the primary validation claim unless promoted by a frozen experiment
that satisfies the same evidence requirements.

## Submission gates

- [x] Clean clone outside a cloud-synced directory.
- [x] Submission branch created from current public `main`.
- [x] Full macOS CPU suite passes with all staged scientific data required.
- [x] Full Linux CPU/data suite passes in the locked environment (complete run plus audited
  crop-CFT focused closure; definitive single-command rerun remains required).
- [x] Exact upstream CTSM source tag/commit identified.
- [ ] Exact upstream CTSM build and run configuration reproduced and archived.
- [ ] Experiment protocol reviewed and frozen.
- [x] Three-site annual qualification runs complete without disabled conservation checks.
- [ ] Three-site strict scientific-parity gate passes (currently 204/207 cells).
- [ ] One machine-readable results index drives all paper claims.
- [ ] Clean-room reproduction succeeds on Linux.
- [ ] Code, inputs, outputs, and analysis have licence/redistribution clearance.
- [ ] Exact software and experiment releases archived with persistent identifiers.
- [ ] Human domain and numerical reviews completed and recorded.

The current qualification checkpoint is indexed in `repro/results/index.json` and checked
by `scripts/gmd/verify_results_index.py`. The machine-readable-results item remains open
until the frozen definitive campaign adds command, environment, input, checksum, raw-output,
and persistent-archive links for every manuscript value.

The committed-data redistribution ledger is `repro/manifests/redistribution.csv`. Its
2026-08-27 audit clears the verbatim FATES parameter file after restoring the separate
FATES licence, but holds six derived NetCDF artifacts pending source-data permission and
privacy-clean regeneration. The release checklist therefore remains open.

## Rules

1. Never copy a headline value from README files, handoff notes, or an old manuscript.
2. Every paper value must resolve through the results index to a command, commit,
   configuration, inputs, checksums, environment, and raw output.
3. Thresholds and exclusions are fixed before the definitive campaign.
4. Failed and unresolved cases remain visible; no conservation check is disabled to make a
   run complete.
5. “Bit-identical”, “parity within tolerance”, and “scientifically acceptable agreement”
   are distinct claims and must not be interchanged.
6. GitHub is the development location, not the persistent paper archive.

## Directory plan

```text
docs/gmd/
  README.md                 submission scope and gates
  provenance.md             source, repository, environment, and licensing record
  redistribution_audit_2026-08-27.md  committed-data clearance audit
  experiment_protocol.md    frozen scientific protocol
  claims.csv                claim-to-evidence ledger
repro/
  environments/             locked paper environments and hardware records
  configs/                  complete model/experiment configurations
  manifests/                checksums and provenance manifests
  results/                  machine-readable summaries (not large raw output)
scripts/gmd/
  prepare.jl
  simulate.jl
  analyse.jl
  figures.jl
  verify.jl
```

Large inputs and outputs will live in a DOI-backed experiment archive rather than Git.

Audit status: the upstream CTSM tag and commit are now verified, but configuration and
reference-generation provenance remain open; see `reference_audit.md`.

The Linux sensitivity-cell execution contract is implemented in
`scripts/gmd/run_btran_sensitivity_cell.sh`; its JSON result schema is produced by
`scripts/gmd/summarize_btran_oracle.py`. Case creation and input acquisition remain
deliberately outside that runner until the CIME machine definition and redistribution
inventory are frozen.
`scripts/gmd/aggregate_btran_sensitivity.py` combines the prescribed cells while retaining
every result and intentionally has no mechanism for declaring a tolerance pass.
The complete 2026-08-27 macOS/Linux release/Linux `-O0` matrix is recorded in
`repro/configs/bow_btran_oracle.toml` and `provenance.md`. It excludes compiler optimization
as the source of the stable `1.02723e-5` grass residual; independent review and threshold
prespecification remain open.

The corrected strict external-data suite passes its macOS CPU execution gate (27,779 pass,
0 fail, 3 explicit hardware skips); see `strict_suite_2026-08-25.md`. Derivative acceptance,
initialization scope, lake/FATES diagnostics, and synthetic multi-parameter recovery remain
scientific qualification failures and cannot support manuscript claims yet.

The locked Linux ARM64 execution gate also covers 27,779 passing assertions with the same
three hardware skips; see `strict_suite_linux_2026-08-27.md`. Its complete run exposed one
crop-CFT test that bypassed strict-data routing. The corrected group passed 27/27 and its
negative missing-data check fails as required. A definitive one-command rerun remains part
of the frozen paper campaign.

The AD false-positive has now been converted into a real gate: five forcing regimes pass
125/125, while winter fails 9 finite-reference assertions. See
`ad_qualification_2026-08-27.md`. Seasonal-cycle differentiability remains outside the
submission claim unless that trajectory is repaired and independently requalified.

The synthetic calibration false-positive is also closed: scalar `csoilc` recovery passes,
but the joint Medlyn/Vcmax experiment now fails four identifiability/recovery assertions
because its starting objective and gradient are exactly zero at the wrong parameters. See
`calibration_qualification_2026-08-27.md`. Joint parameter-recovery claims remain excluded.
