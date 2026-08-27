# Artifact manifests

The definitive campaign will place SHA-256 manifests here for configuration, input,
reference, result, analysis, and figure artifacts. Large artifacts belong in the linked
persistent experiment archive, not in Git.

Each manifest row must include an artifact ID, role, relative archive path, byte count,
SHA-256 digest, licence/redistribution status, source citation, and producing run ID where
applicable.

`redistribution.csv` is the current fail-closed ledger for data already committed to this
repository. It is narrower than the definitive archive manifest and deliberately retains
unresolved artifacts as `HOLD`. Verify its coverage and checksums with:

```sh
python3 scripts/gmd/verify_redistribution_manifest.py
```

`redistribution_sources.json` connects held artifacts to privacy-safe logical source IDs
and exact hashes. `upstream_terms.json` records the primary-provider terms review and
keeps unresolved source terms machine visible.
