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
- [ ] Full CPU suite passes from a clean environment.
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

The strict external-data suite currently fails (27,766 pass, 1 fail, 3 broken); see
`strict_suite_2026-08-25.md`. In particular, strict-mode coverage, derivative acceptance,
and synthetic parameter-recovery criteria must be repaired before they support claims.
