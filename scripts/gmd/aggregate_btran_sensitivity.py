#!/usr/bin/env python3
"""Aggregate all prescribed Bow BTRAN sensitivity cells without selecting a winner."""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cells = [json.loads(path.read_text()) for path in args.results]
    ids = [cell.get("cell_id") for cell in cells]
    if len(set(ids)) != len(ids):
        parser.error("cell IDs must be unique")
    for path, cell in zip(args.results, cells):
        if cell.get("schema_version") != 1 or cell.get("experiment") != \
                "bow-btran-n2901-numerical-sensitivity":
            parser.error(f"incompatible result schema: {path}")

    patches = sorted({row["patch"] for cell in cells for row in cell["btran"]})
    summary = {}
    for patch in patches:
        values = [row["absolute_error"] for cell in cells for row in cell["btran"]
                  if row["patch"] == patch]
        if len(values) != len(cells):
            parser.error(f"patch {patch} is missing from one or more cells")
        summary[str(patch)] = {
            "absolute_error_min": min(values),
            "absolute_error_max": max(values),
            "all_absolute_errors": dict(zip(ids, values)),
        }

    output = {
        "schema_version": 1,
        "experiment": "bow-btran-n2901-numerical-sensitivity-matrix",
        "cell_ids": ids,
        "cell_count": len(cells),
        "btran_by_patch": summary,
        "acceptance_threshold_frozen": False,
        "note": "All supplied cells are retained; this file does not choose a best configuration.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
