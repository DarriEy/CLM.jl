#!/usr/bin/env python3
"""Summarize a Metal fleet run (fleet.jsonl) into a results-index record.

Usage:
  python3 scripts/gmd/summarize_metal_fleet.py <results_dir> [--out <json>]

Reads <results_dir>/fleet.jsonl (written by run_metal_fleet.sh) and emits a
single JSON record with per-harness rows and environment provenance. It never
parses log text for pass/fail — the exit status is the verdict (structural
classification). Known environment-gated harnesses are annotated, not
excluded: every harness appears in the output.
"""
import argparse
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

# Harnesses that cannot pass on a Metal-only box for hardware reasons.
# They stay in the rows; this map only annotates the reason.
HARDWARE_GATED = {
    "gpu_validate_cuda": "requires NVIDIA CUDA hardware",
    "gpu_validate_amdgpu": "requires AMD ROCm hardware",
    "gpu_ad_reverse_driver_validate": "CUDA-targeted reverse-AD driver harness"
                                      " (errors 'CUDA is not functional' on a Metal box)",
}


def git(args, cwd):
    return subprocess.run(["git", *args], cwd=cwd, capture_output=True,
                          text=True, check=True).stdout.strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rdir = Path(args.results_dir)
    jsonl = rdir / "fleet.jsonl"
    if not jsonl.exists():
        sys.exit(f"missing {jsonl}")

    repo = Path(__file__).resolve().parents[2]
    rows = [json.loads(line) for line in jsonl.read_text().splitlines() if line.strip()]
    for r in rows:
        r["hardware_gated"] = HARDWARE_GATED.get(r["harness"])

    n_pass = sum(1 for r in rows if r["exit_status"] == 0)
    n_gated = sum(1 for r in rows if r["exit_status"] != 0 and r["hardware_gated"])
    failures = [r for r in rows if r["exit_status"] != 0 and not r["hardware_gated"]]

    dirty = bool(git(["status", "--porcelain"], repo))
    record = {
        "schema_version": 1,
        "kind": "metal-device-fleet",
        "recorded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_commit": git(["rev-parse", "HEAD"], repo),
        "tree_dirty": dirty,
        "platform": {
            "os": f"{platform.system()} {platform.release()}",
            "machine": platform.machine(),
        },
        "harness_project": "scripts (see scripts/Manifest.toml at source_commit)",
        "totals": {
            "harnesses": len(rows),
            "exit_zero": n_pass,
            "hardware_gated_nonzero": n_gated,
            "unexplained_nonzero": len(failures),
        },
        "unexplained_failures": [r["harness"] for r in failures],
        "rows": rows,
    }

    out = Path(args.out) if args.out else rdir / "summary.json"
    out.write_text(json.dumps(record, indent=2) + "\n")
    print(f"{out}: {n_pass}/{len(rows)} exit-zero, "
          f"{n_gated} hardware-gated, {len(failures)} unexplained")
    for r in failures:
        print(f"  UNEXPLAINED: {r['harness']} exit={r['exit_status']} "
              f"timed_out={r['timed_out']}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
