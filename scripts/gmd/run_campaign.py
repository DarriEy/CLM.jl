#!/usr/bin/env python3
"""Run a frozen GMD campaign and emit checksummed, machine-readable records."""

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SPEC = ROOT / "repro/configs/campaign_candidate.json"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()


def validate_spec(spec):
    required = {"schema_version", "campaign_id", "protocol_frozen", "source_commit",
                "required_platform", "required_environment", "fixed_environment", "steps"}
    missing = sorted(required - spec.keys())
    if missing:
        raise ValueError(f"missing campaign fields: {', '.join(missing)}")
    if spec["schema_version"] != 1:
        raise ValueError("unsupported campaign schema")
    if not isinstance(spec["protocol_frozen"], bool):
        raise ValueError("protocol_frozen must be boolean")
    if not isinstance(spec["required_environment"], list):
        raise ValueError("required_environment must be a list")
    if not all(isinstance(name, str) and name for name in spec["required_environment"]):
        raise ValueError("required_environment entries must be nonempty strings")
    if not isinstance(spec["fixed_environment"], dict):
        raise ValueError("fixed_environment must be an object")
    if not all(isinstance(key, str) and isinstance(value, str)
               for key, value in spec["fixed_environment"].items()):
        raise ValueError("fixed_environment keys and values must be strings")
    if not spec["steps"]:
        raise ValueError("campaign has no steps")
    seen = set()
    for step in spec["steps"]:
        if set(step) != {"id", "argv"} or not step["id"] or not step["argv"]:
            raise ValueError("each step must contain only nonempty id and argv fields")
        if step["id"] in seen:
            raise ValueError(f"duplicate step id: {step['id']}")
        seen.add(step["id"])
        if not all(isinstance(arg, str) and arg for arg in step["argv"]):
            raise ValueError(f"step argv must be nonempty strings: {step['id']}")


def assert_ready(spec):
    blockers = []
    if not spec["protocol_frozen"]:
        blockers.append("protocol_frozen is false")
    expected = spec["source_commit"]
    if not re.fullmatch(r"[0-9a-f]{40}", expected):
        blockers.append("source_commit is not a full frozen Git SHA")
    if platform.system().lower() != spec["required_platform"].lower():
        blockers.append(f"platform is {platform.system()}, expected {spec['required_platform']}")
    status = git("status", "--porcelain=v1", "--untracked-files=all")
    if status:
        blockers.append("source worktree is dirty")
    head = git("rev-parse", "HEAD")
    if re.fullmatch(r"[0-9a-f]{40}", expected) and head != expected:
        # A committed spec cannot contain its own commit hash, so the original
        # HEAD == source_commit rule was unsatisfiable for any frozen spec.
        # The frozen source_commit is the PARENT of the freeze commit — the
        # release-candidate code state — and the freeze commit itself may only
        # flip the campaign configuration under repro/configs/. Anything else
        # between source_commit and HEAD still fails the preflight.
        parent = git("rev-parse", "HEAD^")
        if parent != expected:
            blockers.append(f"HEAD {head} and its parent {parent} both differ "
                            f"from frozen source_commit {expected}")
        else:
            changed = git("diff", "--name-only", f"{expected}..HEAD").splitlines()
            outside = [f for f in changed if not f.startswith("repro/configs/")]
            if outside:
                blockers.append("freeze commit touches files outside repro/configs/: "
                                + ", ".join(outside))
    for name in spec["required_environment"]:
        value = os.environ.get(name, "")
        if not value:
            blockers.append(f"required environment variable {name} is unset")
        elif not Path(value).exists():
            blockers.append(f"{name} path does not exist: {value}")
    if blockers:
        raise RuntimeError("campaign preflight failed:\n- " + "\n- ".join(blockers))
    return head


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def run(spec_path, output):
    spec = json.loads(spec_path.read_text())
    validate_spec(spec)
    try:
        committed_spec = spec_path.relative_to(ROOT)
    except ValueError as error:
        raise RuntimeError("campaign specification must be inside the source checkout") from error
    commit = assert_ready(spec)
    if output == ROOT or ROOT in output.parents:
        raise RuntimeError("campaign output must be outside the source checkout")
    if output.exists() and any(output.iterdir()):
        raise RuntimeError(f"output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env.update({str(key): str(value) for key, value in spec["fixed_environment"].items()})
    provenance = {
        "schema_version": 1,
        "campaign_id": spec["campaign_id"],
        "started_utc": utc_now(),
        "source_commit": commit,
        "spec": str(committed_spec),
        "spec_sha256": sha256(spec_path),
        "project_sha256": sha256(ROOT / "Project.toml"),
        "manifest_sha256": sha256(ROOT / "Manifest.toml") if (ROOT / "Manifest.toml").is_file() else None,
        "required_environment": {name: env[name] for name in spec["required_environment"]},
        "steps": [],
    }
    write_json(output / "campaign.json", provenance)

    environment_path = output / "environment.json"
    environment_error_path = output / "environment.stderr.log"
    with environment_path.open("wb") as stream, environment_error_path.open("wb") as errors:
        environment_result = subprocess.run(
            ["julia", f"--project={ROOT}", str(ROOT / "scripts/gmd/record_environment.jl")],
            cwd=ROOT, env=env, stdout=stream, stderr=errors, check=False)
    provenance["environment_exit_code"] = environment_result.returncode
    provenance["environment_sha256"] = sha256(environment_path)
    provenance["environment_stderr_sha256"] = sha256(environment_error_path)
    if environment_result.returncode != 0:
        provenance["finished_utc"] = utc_now()
        provenance["status"] = "fail"
        write_json(output / "campaign.json", provenance)
        return 1

    campaign_ok = True
    for step in spec["steps"]:
        log_path = output / f"{step['id']}.log"
        started = utc_now()
        clock = time.monotonic()
        with log_path.open("wb") as log:
            result = subprocess.run(step["argv"], cwd=ROOT, env=env, stdout=log,
                                    stderr=subprocess.STDOUT, check=False)
        record = {
            "id": step["id"], "argv": step["argv"], "started_utc": started,
            "elapsed_seconds": time.monotonic() - clock, "exit_code": result.returncode,
            "log": log_path.name, "log_sha256": sha256(log_path),
        }
        provenance["steps"].append(record)
        campaign_ok = campaign_ok and result.returncode == 0
        write_json(output / "campaign.json", provenance)
        if result.returncode != 0:
            break

    provenance["finished_utc"] = utc_now()
    provenance["status"] = "pass" if campaign_ok and len(provenance["steps"]) == len(spec["steps"]) else "fail"
    write_json(output / "campaign.json", provenance)
    return 0 if provenance["status"] == "pass" else 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-spec", action="store_true")
    args = parser.parse_args()
    spec_path = args.spec.resolve()
    spec = json.loads(spec_path.read_text())
    validate_spec(spec)
    if args.check_spec:
        print(f"valid campaign candidate: {spec['campaign_id']} ({len(spec['steps'])} steps)")
        return 0
    if args.output is None:
        parser.error("--output is required unless --check-spec is used")
    return run(spec_path, args.output.resolve())


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(2)
