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
