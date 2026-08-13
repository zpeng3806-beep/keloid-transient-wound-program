#!/usr/bin/env python3
"""Mechanically freeze Gate A3 signatures from discovery-only Gate A2 outputs."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "06_results" / "gateA"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


ranking = read_tsv(RESULTS / "A2_MODULE_CORE_RANKING.tsv")
robustness = read_tsv(RESULTS / "A2_GENE_ROBUSTNESS.tsv")
metrics = read_tsv(RESULTS / "A2_MODULE_METRICS.tsv")

metric_by_module = {row["module"]: row for row in metrics}
primary_metric = metric_by_module.get("D7_M3")
assert primary_metric is not None, "D7_M3 does not exist"
assert primary_metric["QUALITY"] == "HIGH_QUALITY", "D7_M3 is not HIGH_QUALITY"
assert float(primary_metric["ACTIVATION_EFFECT"]) > 0, "activation effect is not positive"
assert float(primary_metric["ATTENUATION_EFFECT"]) > 0, "attenuation effect is not positive"
assert int(primary_metric["activation_sign_preserved"]) == 3, "activation LODO failed"
assert int(primary_metric["attenuation_sign_preserved"]) == 3, "attenuation LODO failed"
assert primary_metric["SUBTYPE_SENSITIVITY"] == "BOTH_MAJOR_SUBTYPES", "subtype preservation failed"

robust_by_gene = {row["gene"]: row for row in robustness}
primary_ranked = sorted(
    (row for row in ranking if row["module"] == "D7_M3"),
    key=lambda row: int(row["MODULE_CORE_RANK"]),
)
assert len(primary_ranked) >= 100, "fewer than 100 D7_M3 ranked genes"
assert [int(row["MODULE_CORE_RANK"]) for row in primary_ranked] == list(
    range(1, len(primary_ranked) + 1)
), "D7_M3 ranks are not complete and unique"

signature_fields = [
    "rank",
    "gene",
    "module",
    "activation",
    "attenuation",
    "activation_concordance",
    "attenuation_concordance",
    "LODO_retention",
    "rank_stability",
    "technical_flag",
]


def signature_rows(n: int | None) -> list[dict[str, object]]:
    selected = primary_ranked if n is None else primary_ranked[:n]
    output = []
    for ranked in selected:
        gene = ranked["gene"]
        robust = robust_by_gene.get(gene)
        assert robust is not None, f"missing robustness row for {gene}"
        output.append(
            {
                "rank": int(ranked["MODULE_CORE_RANK"]),
                "gene": gene,
                "module": "D7_M3",
                "activation": robust["activation"],
                "attenuation": robust["attenuation"],
                "activation_concordance": robust["activation_concordance"],
                "attenuation_concordance": robust["attenuation_concordance"],
                "LODO_retention": robust["number_of_LODO_runs_retaining_candidate"],
                "rank_stability": robust["rank_stability"],
                "technical_flag": robust["TECHNICAL_GENE_FLAG"],
            }
        )
    return output


outputs = {
    "TWPS_PRIMARY_D7_M3_CORE100.tsv": signature_rows(100),
    "TWPS_SENSITIVITY_D7_M3_CORE50.tsv": signature_rows(50),
    "TWPS_SENSITIVITY_D7_M3_CORE25.tsv": signature_rows(25),
    "TWPS_SENSITIVITY_D7_M3_FULL.tsv": signature_rows(None),
}
for filename, rows in outputs.items():
    write_tsv(RESULTS / filename, rows, signature_fields)

secondary_order = ["D7_M1", "D7_M2", "D1_M2", "D1_M1"]
secondary_rows: list[dict[str, object]] = []
for module in secondary_order:
    metric = metric_by_module[module]
    module_rows = sorted(
        (row for row in ranking if row["module"] == module),
        key=lambda row: int(row["MODULE_CORE_RANK"]),
    )
    assert len(module_rows) == int(metric["N_GENES"]), f"membership mismatch for {module}"
    for row in module_rows:
        secondary_rows.append(
            {
                "module": module,
                "gene": row["gene"],
                "module_rank": int(row["MODULE_CORE_RANK"]),
                "quality": metric["QUALITY"],
                "peak_time": metric["PEAK_TIME"],
            }
        )

write_tsv(
    RESULTS / "TWPS_SECONDARY_MODULES_LOCKED.tsv",
    secondary_rows,
    ["module", "gene", "module_rank", "quality", "peak_time"],
)

assert len(outputs["TWPS_PRIMARY_D7_M3_CORE100.tsv"]) == 100
assert len(outputs["TWPS_SENSITIVITY_D7_M3_CORE50.tsv"]) == 50
assert len(outputs["TWPS_SENSITIVITY_D7_M3_CORE25.tsv"]) == 25
assert len(outputs["TWPS_SENSITIVITY_D7_M3_FULL.tsv"]) == int(primary_metric["N_GENES"])
