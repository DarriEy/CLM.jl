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
WSC, Copernicus, MODIS, CTSM, and SYMFLUENCE chain. The exact source history file and
extraction have now been recovered: every one of the committed subset's 28 arrays is
bit-identical to `Bow_at_Banff_lumped.clm2.h0.2008-12-30-00000.nc` at the recorded S04
checksum. `scripts/gmd/extract_bow_reference.jl` reproduces the subset without aggregation
or recalculation. Regeneration matches dimensions, raw array values, variable attributes,
and global attributes; the NetCDF container checksum is not expected to match across
serialization libraries. This closes content provenance, not redistribution permission.

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

## Addendum — 2026-09-02 synthesis pass

D02, D03, D05, D06, and D07 are **CLEARED by content replacement**, not by permission:

- The three surfdata fixtures (D05-D07) were regenerated from scratch by
  `scripts/gen_synthetic_surfdata.jl` — every value is a constant, a standard
  pedotransfer formula on a synthetic loam texture, or an analytic seasonal cycle, at a
  generic synthetic point (45 N, 262.5 E). No S01 bytes remain. Their consumers pass
  unchanged (lake water balance 4/4, eddy ks 14/14, 2-m diagnostics 48/48, glacier
  robustness 3/3, init-cold wiring 36/36 + NaN gate 3/3, run_clm harness 15/15), and the
  subgrid structure they exercise is identical to the replaced fixtures.
- The two hillslope fixtures (D02-D03) had exactly two upstream-derived arrays each
  (LONGXY/LATIXY) plus a leaked absolute base path in a global attribute.
  `scripts/gmd/synthesize_hillslope_grid.jl` replaced the coordinates with analytic
  values (generic (-10, 300) for the single-point file; a Fibonacci sphere lattice of
  the same 488-cell count for the unstructured file) and removed the leaked attribute.
  The catena geometry was always synthetic. Their consumers pass unchanged (hillslope
  hydrology 45/45, end-to-end 44/44, full-driver E2E 21/21).
- Real-site/real-grid variants for Fortran-side parity remain locally regenerable via
  `scripts/gen_lake_surfdata.jl`, `scripts/gen_glacier_surfdata.jl`, and
  `scripts/make_hillslope_surfdata.jl` from `SYMFLUENCE_DATA`; those outputs are not
  tracked.
- **D04 keeps HOLD** with a recorded disposition: exclude from public archives and ship
  the lawful regeneration path (`scripts/gmd/extract_bow_reference.jl` plus the banked
  CTSM no-rebuild recipe), stating the restriction in Code and data availability with a
  reviewer-access route — the GMD-policy-compliant form of a restriction outside the
  authors' control.

`verify_redistribution_manifest.py` and `verify_source_release_policy.py` (both
unmodified) pass: 6 cleared / 1 held, 6 includable artifacts, release gate still
fail-closed on D04 as designed. The remaining human item from the original audit is
unchanged: confirm the exact upstream products used for the big experiment archive
match the reviewed provider terms before the experiment deposit is registered.
