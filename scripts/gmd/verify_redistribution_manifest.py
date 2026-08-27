#!/usr/bin/env python3
"""Verify the committed-data redistribution ledger and its fail-closed status."""

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "repro/manifests/redistribution.csv"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    with MANIFEST.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
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
    statuses = {row["redistribution_status"] for row in rows}
    if "HOLD" not in statuses:
        raise AssertionError("release must not appear cleared while unresolved data remain")
    print(f"verified {len(rows)} committed data artifacts; "
          f"{sum(row['redistribution_status'] == 'CLEARED' for row in rows)} cleared, "
          f"{sum(row['redistribution_status'] == 'HOLD' for row in rows)} on hold")


if __name__ == "__main__":
    main()
