# Machine-readable GMD results

Definitive small summaries are promoted here after review. Raw cell directories produced by
`scripts/gmd/run_btran_sensitivity_cell.sh` are ignored because they contain logs and host
metadata that must first be checked for private paths. The complete reviewed bundle belongs
in the DOI-backed experiment archive; only its index and checksums are committed.

`index.json` is the current qualification checkpoint that drives the claim ledger. It is
not the definitive campaign index. Validate its evidence paths, cross-file metrics, and
failed-claim guards with:

```sh
python3 scripts/gmd/verify_results_index.py
```
