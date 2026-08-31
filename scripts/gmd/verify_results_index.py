#!/usr/bin/env python3
"""Check that the GMD qualification index is internally honest and resolvable."""

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INDEX_PATH = ROOT / "repro/results/index.json"
CLAIMS_PATH = ROOT / "docs/gmd/claims.csv"


def require_path(relative):
    path = ROOT / relative
    if not path.is_file():
        raise AssertionError(f"indexed evidence does not exist: {relative}")


def main():
    index = json.loads(INDEX_PATH.read_text())
    if index["schema_version"] != 1:
        raise AssertionError("unsupported results-index schema")

    gates = index["gates"]
    for gate in gates.values():
        for key in ("evidence", "matrix"):
            if key in gate:
                require_path(gate[key])

    matrix = json.loads((ROOT / gates["bow_btran_process_oracle"]["matrix"]).read_text())
    if matrix["cell_count"] != gates["bow_btran_process_oracle"]["cells"]:
        raise AssertionError("BTRAN cell count differs between matrix and index")
    if matrix["acceptance_threshold_frozen"]:
        raise AssertionError("BTRAN matrix unexpectedly declares a frozen threshold")

    with CLAIMS_PATH.open(newline="") as stream:
        claims = {row["claim_id"]: row for row in csv.DictReader(stream)}

    # Invariants that hold until their specific evidence exists — update these
    # TOGETHER with the recorded evidence, never to tidy a red result away.
    if claims["C01"]["status"] != "FAIL-QUALIFICATION":
        raise AssertionError("C01 stays failed until the audited Bow annual rerun "
                             "passes the unchanged 207-cell gate")
    if claims["C08"]["status"] != "PARTIAL-QUALIFICATION":
        raise AssertionError("reproducibility must remain partial before the DOI campaign")
    if gates["strict_scientific_parity"]["claim_supported"]:
        raise AssertionError("204/207 parity cannot support the scoped parity claim")

    # Data-driven consistency: a gate may support a claim only from a passing
    # status, and every PASS-QUALIFICATION claim row must resolve to an
    # existing evidence file.
    for name, gate in gates.items():
        supported = gate.get("claim_supported") or \
            gate.get("seasonal_cycle_claim_supported") or \
            gate.get("joint_recovery_claim_supported")
        status = gate.get("status", "")
        if supported and not (status.startswith("pass") or status == "evidence-recorded"):
            raise AssertionError(f"gate {name} supports a claim from status {status!r}")
    for claim_id, row in claims.items():
        if row["status"].startswith("PASS-QUALIFICATION"):
            require_path(row["evidence_artifact"])

    print(f"verified {len(gates)} gates and {len(claims)} claims")


if __name__ == "__main__":
    main()
