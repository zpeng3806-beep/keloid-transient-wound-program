#!/usr/bin/env python3
"""Prepare deterministic MAGMA inputs without significance filtering."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sumstats", type=Path, required=True)
    parser.add_argument("--pval-output", type=Path, required=True)
    parser.add_argument("--qc-output", type=Path, required=True)
    args = parser.parse_args()

    args.pval_output.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    total = valid = missing_or_invalid_rsid = invalid_p = duplicate_rsid = 0

    with gzip.open(args.sumstats, "rt", newline="") as src, args.pval_output.open(
        "w", newline=""
    ) as dst:
        reader = csv.DictReader(src, delimiter="\t")
        required = {"rsid", "p_value"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise SystemExit(f"missing required fields: {sorted(required)}")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "P"])
        for row in reader:
            total += 1
            rsid = (row.get("rsid") or "").strip()
            if not (rsid.startswith("rs") and rsid[2:].isdigit()):
                missing_or_invalid_rsid += 1
                continue
            if rsid in seen:
                duplicate_rsid += 1
                continue
            seen.add(rsid)
            p_text = (row.get("p_value") or "").strip()
            try:
                p_value = float(p_text)
            except ValueError:
                invalid_p += 1
                continue
            if not 0.0 <= p_value <= 1.0:
                invalid_p += 1
                continue
            writer.writerow([rsid, p_text])
            valid += 1

    with args.qc_output.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(
            [
                ["input_variant_n", total],
                ["valid_unique_rsid_p_n", valid],
                ["missing_or_invalid_rsid_n", missing_or_invalid_rsid],
                ["invalid_p_n", invalid_p],
                ["duplicate_rsid_first_retained_n", duplicate_rsid],
                ["significance_filter_applied", "NO"],
            ]
        )


if __name__ == "__main__":
    main()
