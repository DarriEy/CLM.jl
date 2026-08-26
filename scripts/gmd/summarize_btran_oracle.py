#!/usr/bin/env python3
"""Create one machine-readable Bow BTRAN sensitivity-cell result."""

import argparse
import hashlib
import json
import platform
import re
import subprocess
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command(*args):
    try:
        return subprocess.check_output(args, text=True, stderr=subprocess.STDOUT).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def git_record(root):
    status = command("git", "-C", str(root), "status", "--porcelain=v1",
                     "--untracked-files=all")
    return {
        "commit": command("git", "-C", str(root), "rev-parse", "HEAD"),
        "describe": command("git", "-C", str(root), "describe", "--always", "--dirty", "--tags"),
        "clean": status == "",
        "remote": command("git", "-C", str(root), "remote", "get-url", "origin"),
    }


def parse_comparison(path):
    rows = []
    lines = path.read_text().splitlines()
    start = next((i for i, line in enumerate(lines)
                  if line.strip() == "patch,Julia,Fortran,abs_diff,rel_diff"), None)
    if start is None:
        raise ValueError("comparison log has no BTRAN result table")
    for line in lines[start + 1:]:
        if not re.match(r"^\d+,", line):
            if rows:
                break
            continue
        patch, julia, fortran, absolute, relative = line.split(",")
        rows.append({
            "patch": int(patch),
            "julia": float(julia),
            "fortran": float(fortran),
            "absolute_error": float(absolute),
            "relative_error": float(relative),
        })
    if not rows:
        raise ValueError("comparison log has an empty BTRAN result table")
    initial_match = re.search(
        r"injected initial state global max\|rel\| = ([0-9.eE+-]+)", "\n".join(lines))
    return rows, (float(initial_match.group(1)) if initial_match else None)


def read_trace(path):
    lines = [line.split() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError(f"empty canopy trace: {path}")
    header = lines[0]
    return header, [[float(value) for value in row] for row in lines[1:]]


def trace_metrics(julia_path, fortran_path):
    jhead, jrows = read_trace(julia_path)
    fhead, frows = read_trace(fortran_path)
    if jhead != fhead:
        raise ValueError("Julia and Fortran canopy trace schemas differ")
    index = {name: i for i, name in enumerate(jhead)}
    required = ("itlef", "p", "efe")
    missing = [name for name in required if name not in index]
    if missing:
        raise ValueError(f"canopy trace is missing fields: {missing}")

    def by_key(rows):
        return {(int(row[index["itlef"]]), int(row[index["p"]])): row for row in rows}

    jmap, fmap = by_key(jrows), by_key(frows)
    if set(jmap) != set(fmap):
        raise ValueError("Julia and Fortran iteration keys differ")
    patches = sorted({patch for _, patch in jmap})
    return {
        "records": {"julia": len(jrows), "fortran": len(frows)},
        "iterations_by_patch": {
            str(patch): max(iteration for iteration, p in jmap if p == patch)
            for patch in patches
        },
        "first_iteration_efe_absolute_error": {
            str(patch): abs(jmap[(1, patch)][index["efe"]] -
                            fmap[(1, patch)][index["efe"]])
            for patch in patches
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-id", required=True)
    parser.add_argument("--comparison-log", type=Path, required=True)
    parser.add_argument("--julia-trace", type=Path, required=True)
    parser.add_argument("--fortran-trace", type=Path, required=True)
    parser.add_argument("--ctsm-root", type=Path, required=True)
    parser.add_argument("--ctsm-executable", type=Path, required=True)
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for path in (args.comparison_log, args.julia_trace, args.fortran_trace):
        if not path.is_file():
            parser.error(f"missing input: {path}")
    if not args.ctsm_executable.is_file():
        parser.error(f"missing CTSM executable: {args.ctsm_executable}")

    case_artifacts = {}
    for name in ("CaseStatus", "Macros.make", "env_build.xml", "env_case.xml",
                 "env_mach_pes.xml", "env_run.xml"):
        path = args.case_dir / name
        if path.is_file():
            case_artifacts[name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}

    btran, initial_error = parse_comparison(args.comparison_log)
    result = {
        "schema_version": 1,
        "experiment": "bow-btran-n2901-numerical-sensitivity",
        "cell_id": args.cell_id,
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "gfortran": command("gfortran", "--version"),
            "julia": command("julia", "--version"),
        },
        "ctsm": {
            "git": git_record(args.ctsm_root),
            "executable": {
                "path": str(args.ctsm_executable),
                "bytes": args.ctsm_executable.stat().st_size,
                "sha256": sha256(args.ctsm_executable),
            },
            "case_artifacts": case_artifacts,
        },
        "initial_state_global_max_relative_error": initial_error,
        "btran": btran,
        "canopy_trace": trace_metrics(args.julia_trace, args.fortran_trace),
        "artifacts": {
            str(path): {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in (args.comparison_log, args.julia_trace, args.fortran_trace)
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
