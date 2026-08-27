# Data acquisition for the public source artifact

The public source artifact excludes every data file whose redistribution status is `HOLD`
in `repro/manifests/redistribution.csv`. Exclusion prevents redistribution; it does not make
the scientific campaign independently reproducible. Exact checksums and dependency chains
remain visible in `redistribution_sources.json`.

## Excluded artifacts

- D02: regenerate the Aripuana hillslope file with
  `scripts/make_hillslope_surfdata.jl` only after acquiring S02 under confirmed terms.
- D03: acquire exact S03 through the standard CTSM input-data workflow, verify its SHA-256,
  then run `scripts/make_hillslope_surfdata.jl`. Do not republish S03 based only on public
  download availability.
- D04: acquire the exact S04 history file under the complete upstream terms, verify its
  SHA-256, then run `scripts/gmd/extract_bow_reference.jl`.
- D05-D07: regenerate with `scripts/gen_glacier_surfdata.jl` and
  `scripts/gen_lake_surfdata.jl` only after S01's OpenGeoHub mosaic dependency is replaced
  by a confirmed direct-NASA path or receives explicit licence clarification.

Strict tests that require excluded files must fail as missing in a source-only archive.
They are enabled only after a separately archived, checksum-verified scientific-data bundle
has passed the complete redistribution review. Do not substitute approximately similar
inputs and report the resulting run as reproduction of the paper campaign.
