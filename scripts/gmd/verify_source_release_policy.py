#!/usr/bin/env python3
"""Verify that the source-release policy excludes every held data artifact."""

import csv
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LEDGER = ROOT / "repro/manifests/redistribution.csv"
POLICY = ROOT / "repro/manifests/source_release_policy.json"


def main():
    with LEDGER.open(newline="") as stream:
        rows = {row["artifact_id"]: row for row in csv.DictReader(stream)}
    policy = json.loads(POLICY.read_text())
    actions = policy["artifact_actions"]
    if set(actions) != set(rows):
        raise AssertionError("release policy and redistribution ledger IDs differ")
    for artifact_id, row in rows.items():
        action = actions[artifact_id]["action"]
        if row["redistribution_status"] == "HOLD" and action != "EXCLUDE":
            raise AssertionError(f"held artifact is not excluded: {artifact_id}")
        if action == "INCLUDE" and row["redistribution_status"] != "CLEARED":
            raise AssertionError(f"uncleared artifact is included: {artifact_id}")
        if not actions[artifact_id]["reason"]:
            raise AssertionError(f"release action lacks reason: {artifact_id}")
    for relative in policy["required_archive_files"]:
        if not (ROOT / relative).is_file():
            raise AssertionError(f"required archive file missing: {relative}")
    tracked = set(subprocess.check_output(
        ["git", "-C", str(ROOT), "ls-files"], text=True
    ).splitlines())
    path_exclusions = policy["repository_path_exclusions"]
    for relative, reason in path_exclusions.items():
        if relative not in tracked:
            raise AssertionError(f"repository-path exclusion is not tracked: {relative}")
        if not reason:
            raise AssertionError(f"repository-path exclusion lacks reason: {relative}")
    held_paths = {
        row["path"] for row in rows.values()
        if row["redistribution_status"] == "HOLD"
    }
    overlap = held_paths & set(path_exclusions)
    if overlap:
        raise AssertionError(f"held data must be classified by artifact ID: {sorted(overlap)}")
    if set(policy["required_archive_files"]) & set(path_exclusions):
        raise AssertionError("required archive files cannot be repository-path exclusions")
    if policy["release_gate_complete"]:
        raise AssertionError("candidate policy must not claim complete release clearance")
    included = sum(item["action"] == "INCLUDE" for item in actions.values())
    excluded = sum(item["action"] == "EXCLUDE" for item in actions.values())
    print(f"verified source-release policy: {included} data artifact included, "
          f"{excluded} excluded; {len(path_exclusions)} development paths excluded")


if __name__ == "__main__":
    main()
