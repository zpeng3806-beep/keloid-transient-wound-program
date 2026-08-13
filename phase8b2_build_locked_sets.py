#!/usr/bin/env python3
"""Map locked signatures by exact symbol to NCBI37.3 Entrez IDs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_gene_locations(path: Path) -> tuple[dict[str, str], set[str]]:
    symbol_to_ids: dict[str, list[str]] = {}
    all_ids: set[str] = set()
    with path.open() as src:
        for line in src:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            gene_id, symbol = fields[0], fields[5]
            all_ids.add(gene_id)
            symbol_to_ids.setdefault(symbol, []).append(gene_id)
    unique = {symbol: ids[0] for symbol, ids in symbol_to_ids.items() if len(ids) == 1}
    return unique, all_ids


def read_signature(path: Path) -> list[str]:
    with path.open(newline="") as src:
        rows = list(csv.DictReader(src, delimiter="\t"))
    return [row["gene"].strip() for row in rows]


def read_tested_gene_ids(path: Path) -> set[str]:
    tested: set[str] = set()
    with path.open() as src:
        header = src.readline().split()
        gene_index = header.index("GENE") if "GENE" in header else 0
        for line in src:
            fields = line.split()
            if len(fields) > gene_index:
                tested.add(fields[gene_index])
    return tested


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-loc", type=Path, required=True)
    parser.add_argument("--core100", type=Path, required=True)
    parser.add_argument("--core25", type=Path, required=True)
    parser.add_argument("--core50", type=Path, required=True)
    parser.add_argument("--full", type=Path, required=True)
    parser.add_argument("--eur-genes", type=Path)
    parser.add_argument("--afr-genes", type=Path)
    parser.add_argument("--set-output", type=Path, required=True)
    parser.add_argument("--primary-set-output", type=Path)
    parser.add_argument("--sensitivity-set-output", type=Path)
    parser.add_argument("--mapping-output", type=Path, required=True)
    args = parser.parse_args()

    symbol_map, _ = read_gene_locations(args.gene_loc)
    signatures = {
        "CORE100": read_signature(args.core100),
        "CORE25": read_signature(args.core25),
        "CORE50": read_signature(args.core50),
        "FULL": read_signature(args.full),
    }
    args.set_output.parent.mkdir(parents=True, exist_ok=True)
    with args.set_output.open("w") as out:
        for name, symbols in signatures.items():
            mapped = [symbol_map[symbol] for symbol in symbols if symbol in symbol_map]
            out.write(" ".join([name, *mapped]) + "\n")
    if args.primary_set_output:
        mapped = [symbol_map[symbol] for symbol in signatures["CORE100"] if symbol in symbol_map]
        with args.primary_set_output.open("w") as out:
            out.write(" ".join(["CORE100", *mapped]) + "\n")
    if args.sensitivity_set_output:
        with args.sensitivity_set_output.open("w") as out:
            for name in ("CORE25", "CORE50", "FULL"):
                mapped = [symbol_map[symbol] for symbol in signatures[name] if symbol in symbol_map]
                out.write(" ".join([name, *mapped]) + "\n")

    eur_tested = read_tested_gene_ids(args.eur_genes) if args.eur_genes else set()
    afr_tested = read_tested_gene_ids(args.afr_genes) if args.afr_genes else set()
    with args.mapping_output.open("w", newline="") as out:
        fields = [
            "locked_gene",
            "mapped_gene_id",
            "MAGMA_gene_id",
            "EUR_testable",
            "AFR_testable",
            "mapping_status",
            "notes",
        ]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for symbol in signatures["CORE100"]:
            gene_id = symbol_map.get(symbol, "")
            mapped = bool(gene_id)
            writer.writerow(
                {
                    "locked_gene": symbol,
                    "mapped_gene_id": gene_id,
                    "MAGMA_gene_id": gene_id,
                    "EUR_testable": "YES" if gene_id in eur_tested else "NO",
                    "AFR_testable": "YES" if gene_id in afr_tested else "NO",
                    "mapping_status": "EXACT_SYMBOL_MAPPED" if mapped else "UNMAPPED_EXACT_SYMBOL_ONLY",
                    "notes": "NCBI37.3 exact symbol; no aliases or substitutions" if mapped else "No exact symbol in locked protein-coding NCBI37.3 resource",
                }
            )


if __name__ == "__main__":
    main()
