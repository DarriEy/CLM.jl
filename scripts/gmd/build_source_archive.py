#!/usr/bin/env python3
"""Build a source archive that cannot contain held data artifacts."""

import argparse
import csv
import hashlib
import io
import json
import subprocess
import tarfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
POLICY = ROOT / "repro/manifests/source_release_policy.json"
LEDGER = ROOT / "repro/manifests/redistribution.csv"


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_policy():
    policy = json.loads(POLICY.read_text())
    with LEDGER.open(newline="") as stream:
        rows = {row["artifact_id"]: row for row in csv.DictReader(stream)}
    subprocess.run(["python3", str(ROOT / "scripts/gmd/verify_redistribution_manifest.py")],
                   cwd=ROOT, check=True)
    subprocess.run(["python3", str(ROOT / "scripts/gmd/verify_source_release_policy.py")],
                   cwd=ROOT, check=True)
    excluded = {rows[key]["path"] for key, value in policy["artifact_actions"].items()
                if value["action"] == "EXCLUDE"}
    return policy, rows, excluded


def build(output):
    policy, rows, excluded = load_policy()
    if output == ROOT or ROOT in output.parents:
        raise RuntimeError("source archive output must be outside the checkout")
    if output.exists():
        raise RuntimeError(f"refusing to overwrite archive: {output}")
    status = git("status", "--porcelain=v1", "--untracked-files=all")
    if status:
        raise RuntimeError("source worktree is dirty")
    commit = git("rev-parse", "HEAD")
    files = [line for line in git("ls-files").splitlines() if line not in excluded]
    missing_required = set(policy["required_archive_files"]) - set(files)
    if missing_required:
        raise RuntimeError(f"required files absent from archive: {sorted(missing_required)}")
    prefix = f"CLM.jl-{commit[:12]}"
    manifest = {
        "schema_version": 1,
        "source_commit": commit,
        "included_file_count": len(files),
        "included_files_sha256": {relative: sha256(ROOT / relative) for relative in files},
        "excluded_artifacts": [
            {"artifact_id": key, "path": rows[key]["path"],
             "reason": value["reason"]}
            for key, value in policy["artifact_actions"].items()
            if value["action"] == "EXCLUDE"
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(output, "w:gz") as archive:
        for relative in files:
            archive.add(ROOT / relative, arcname=f"{prefix}/{relative}", recursive=False)
        payload = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
        info = tarfile.TarInfo(f"{prefix}/SOURCE_ARCHIVE_MANIFEST.json")
        info.size = len(payload)
        info.mode = 0o644
        archive.addfile(info, io.BytesIO(payload))
    print(json.dumps({"archive": str(output), "sha256": sha256(output), **manifest}, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-policy", action="store_true")
    args = parser.parse_args()
    if args.check_policy:
        load_policy()
        print("source archive policy is internally valid")
        return
    if args.output is None:
        parser.error("--output is required unless --check-policy is used")
    build(args.output.resolve())


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"error: {error}")
