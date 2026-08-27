# Licence and redistribution audit — 2026-08-27

Status: **repository code attribution improved; committed-data release gate FAIL with six
artifacts on hold**.

## Scope and method

This audit covers every tracked file under `data/`, `test/reference_data/`, and
`test_inputs/`, plus the repository `LICENSE` and `NOTICE`. SHA-256 identity, embedded
NetCDF provenance, repository generators, and the staged CTSM 5.3.012 licences were
checked. It is a provenance audit, not legal advice and not permission inferred from a
file being publicly obtainable.

The authoritative file-level ledger is `repro/manifests/redistribution.csv`; run
`python3 scripts/gmd/verify_redistribution_manifest.py` after any fixture change.

## Findings

- Original repository contributions state MIT terms in `LICENSE`; CTSM-derived portions
  retain the exact CTSM 5.3.012 BSD notice in `NOTICE`.
- `data/fates/fates_params_default.cdl` is byte-identical to the file at the frozen CTSM
  revision. FATES has separate BSD terms that were absent from this distribution. Those
  terms are now reproduced at `LICENSES/FATES.txt` and referenced by `NOTICE`; D01 is
  cleared provided that attribution remains packaged.
- Both committed hillslope files declare that real coordinates were copied from external
  surface data. Neither contains licence evidence, and both expose a developer-specific
  absolute path in a global attribute. D02-D03 remain on hold.
- The Bow Fortran reference is generated model output, but its metadata identifies a
  SYMFLUENCE DDS workflow without recording licences for forcing, parameters, restart, or
  surface inputs. D04 remains on hold; generated status alone is insufficient evidence.
- The glacier and two lake fixtures were mechanically derived from Bow surface data. Their
  global `source = SYMFLUENCE` attribute is not a licence. D05-D07 remain on hold.

No staged scientific archive under `SYMFLUENCE_DATA` is approved for redistribution by
this audit. The definitive archive must extend the ledger to every forcing, restart,
surface, parameter, observation, reference, and generated-output artifact.

## Release disposition

Do not mark the GMD redistribution checklist complete. Before a public artifact release,
either obtain and record source-specific permission for D02-D07 and regenerate privacy-
clean files, or exclude those files and provide lawful acquisition/generation recipes that
still permit reproduction. Human review must confirm whether input terms attach to the
derived NetCDF products.
