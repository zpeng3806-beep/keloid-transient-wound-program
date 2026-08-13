#!/usr/bin/env python3
"""Freeze the two prelocked primary competitive MAGMA results."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_gsa(path: Path) -> dict[str, str]:
    with path.open() as src:
        header = None
        for line in src:
            if line.startswith("VARIABLE"):
                header = line.split()
                continue
            if header and line.strip() and not line.startswith("#"):
                values = line.split()
                return dict(zip(header, values))
    raise ValueError(f"no result row in {path}")


def classify(eur: dict[str, str], afr: dict[str, str]) -> str:
    betas = [float(eur["BETA"]), float(afr["BETA"])]
    ps = [float(eur["P"]), float(afr["P"])]
    if all(beta > 0 for beta in betas):
        if all(p < 0.025 for p in ps):
            return "VERY_STRONG_CROSS_ANCESTRY_SUPPORT"
        if any(p < 0.025 for p in ps):
            return "STRONG_PROGRAM_LEVEL_GENETIC_SUPPORT"
        return "MODERATE_PROGRAM_LEVEL_SUPPORT"
    if betas[0] * betas[1] < 0:
        return "CONTRADICTORY_PROGRAM_LEVEL_EVIDENCE"
    if all(beta <= 0 for beta in betas):
        return "NO_PROGRAM_LEVEL_SUPPORT"
    return "WEAK_PROGRAM_LEVEL_SUPPORT"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eur", type=Path, required=True)
    parser.add_argument("--afr", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = {"EUR": read_gsa(args.eur), "AFR": read_gsa(args.afr)}
    classification = classify(rows["EUR"], rows["AFR"])
    accession = {"EUR": "GCST90652488", "AFR": "GCST90652489"}
    with args.output.open("w", newline="") as out:
        fields = ["ancestry", "accession", "CORE100_testable_n", "beta", "SE", "P", "direction", "threshold", "classification"]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for ancestry in ("EUR", "AFR"):
            row = rows[ancestry]
            beta = float(row["BETA"])
            writer.writerow(
                {
                    "ancestry": ancestry,
                    "accession": accession[ancestry],
                    "CORE100_testable_n": row["NGENES"],
                    "beta": row["BETA"],
                    "SE": row["SE"],
                    "P": row["P"],
                    "direction": "POSITIVE" if beta > 0 else "NEGATIVE" if beta < 0 else "ZERO",
                    "threshold": "P<0.025",
                    "classification": classification,
                }
            )


if __name__ == "__main__":
    main()
