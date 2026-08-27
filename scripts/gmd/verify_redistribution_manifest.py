#!/usr/bin/env python3
"""Verify the committed-data redistribution ledger and its fail-closed status."""

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "repro/manifests/redistribution.csv"
SOURCES = ROOT / "repro/manifests/redistribution_sources.json"
TERMS = ROOT / "repro/manifests/upstream_terms.json"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    with MANIFEST.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    source_trace = json.loads(SOURCES.read_text())
    terms = json.loads(TERMS.read_text())
    if terms["release_decision"]["all_held_artifacts_cleared"]:
        raise AssertionError("upstream terms must not claim full clearance")
    unresolved = set(terms["release_decision"]["unresolved_terms"])
    if not unresolved:
        raise AssertionError("unresolved upstream terms unexpectedly empty")
    expected = set()
    for prefix in ("data", "test/reference_data", "test_inputs"):
        expected.update(path.relative_to(ROOT).as_posix() for path in (ROOT / prefix).rglob("*")
                        if path.is_file())
    indexed = {row["path"] for row in rows}
    if indexed != expected:
        raise AssertionError(f"manifest coverage differs: missing={expected-indexed}, extra={indexed-expected}")
    for row in rows:
        path = ROOT / row["path"]
        if sha256(path) != row["sha256"]:
            raise AssertionError(f"checksum mismatch: {row['path']}")
        if row["redistribution_status"] not in {"CLEARED", "HOLD", "EXCLUDE"}:
            raise AssertionError(f"invalid status: {row['artifact_id']}")
        if row["redistribution_status"] == "CLEARED":
            evidence = ROOT / row["license_evidence"]
            if not evidence.is_file():
                raise AssertionError(f"missing licence evidence: {row['artifact_id']}")
        elif row["redistribution_status"] == "HOLD":
            evidence = ROOT / row["license_evidence"]
            if not evidence.is_file():
                raise AssertionError(f"missing hold evidence: {row['artifact_id']}")
            derivation = source_trace["derivations"].get(row["artifact_id"])
            if derivation is None:
                raise AssertionError(f"missing source derivation: {row['artifact_id']}")
            for source_id in derivation["source_ids"]:
                source = source_trace["sources"][source_id]
                if source["license_evidence_complete"]:
                    raise AssertionError(f"held source unexpectedly claims clearance: {source_id}")
                if not unresolved.intersection(source["terms_ids"]):
                    raise AssertionError(f"held source has no unresolved upstream term: {source_id}")
    statuses = {row["redistribution_status"] for row in rows}
    if "HOLD" not in statuses:
        raise AssertionError("release must not appear cleared while unresolved data remain")
    print(f"verified {len(rows)} committed data artifacts; "
          f"{sum(row['redistribution_status'] == 'CLEARED' for row in rows)} cleared, "
          f"{sum(row['redistribution_status'] == 'HOLD' for row in rows)} on hold")


if __name__ == "__main__":
    main()
