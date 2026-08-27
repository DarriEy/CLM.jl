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
  generators first byte-copy S01 and change only two land-unit variables; every other
  variable was verified identical. Their global `source = SYMFLUENCE` attribute is not a
  licence. D05-D07 remain on hold.

`repro/manifests/redistribution_sources.json` now records privacy-safe logical source IDs,
source hashes, copied fields, and exact transformation evidence. S01 is the Bow surface
file declared by the local campaign configuration to incorporate Copernicus DEM and MODIS
land-cover inputs. S02 is the Aripuana surface file; S03 is the exact CTSM input-data
surface file. D04's metadata and surviving final-evaluation directory identify an RDRS,
WSC, Copernicus, MODIS, CTSM, and SYMFLUENCE chain, but no exact extraction recipe for the
committed 365-value subset was found.

The three fixture generators now require an explicit source or `SYMFLUENCE_DATA`, record a
source checksum, and avoid writing absolute source paths. This prevents future privacy
leakage but does not retroactively clear or alter the held artifacts.

Primary-provider terms were reviewed separately in `upstream_terms_2026-08-27.md` and
`repro/manifests/upstream_terms.json`. RDRS, WSC, NASA MODIS, and Copernicus DEM have usable
reuse terms subject to attribution and exact-product confirmation. The exact OpenGeoHub
MODIS mosaic record displays no licence, and no file-level redistribution statement was
found for the exact CTSM input surface file. Those unresolved terms keep every held
derivation on hold.

No staged scientific archive under `SYMFLUENCE_DATA` is approved for redistribution by
this audit. The definitive archive must extend the ledger to every forcing, restart,
surface, parameter, observation, reference, and generated-output artifact.

## Release disposition

Do not mark the GMD redistribution checklist complete. Before a public artifact release,
either obtain and record source-specific permission for D02-D07 and regenerate privacy-
clean files, or exclude those files and provide lawful acquisition/generation recipes that
still permit reproduction. Human review must confirm whether input terms attach to the
derived NetCDF products.
