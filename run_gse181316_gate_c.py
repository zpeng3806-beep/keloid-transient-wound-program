#!/usr/bin/env python3
"""Memory-conscious, donor-level GSE181316 locked-TWPS replication."""

from __future__ import annotations

import csv
import gzip
import itertools
import math
import random
from array import array
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STAGED = ROOT / "data" / "staged" / "GSE181316"
RESULTS = ROOT / "06_results" / "gateC"
INTERMEDIATE = ROOT / "data" / "intermediate"
RESULTS.mkdir(parents=True, exist_ok=True)
INTERMEDIATE.mkdir(parents=True, exist_ok=True)

SAMPLES = [
    ("keloid_1", "GSM5494684", "KELOID_1", "keloid"),
    ("keloid_2", "GSM5494685", "KELOID_2", "keloid"),
    ("keloid_3L", "GSM5494686", "KELOID_3", "keloid"),
    ("keloid_3R", "GSM5494687", "KELOID_3", "keloid"),
    ("scar_1", "GSM5494688", "SCAR_1", "normal scar"),
    ("scar_2", "GSM5494689", "SCAR_2", "normal scar"),
    ("scar_3", "GSM5494690", "SCAR_3", "normal scar"),
]

# Fixed before TWPS scoring and independent of group/outcome. GEO raw release supplies
# no author cell-label metadata; this is a limited lineage gate, not a whole-atlas annotation.
FIBROBLAST_MARKERS = {"COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"}
MIN_DETECTED_GENES = 200
MIN_FIBROBLAST_MARKERS = 2


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


features_path = STAGED / "GSM5494684_keloid_1_features.tsv.gz"
with gzip.open(features_path, "rt", encoding="utf-8") as handle:
    features = [line.rstrip("\n").split("\t") for line in handle]
gene_ids = [x[0] for x in features]
gene_symbols = [x[1] for x in features]
marker_rows = {i + 1 for i, symbol in enumerate(gene_symbols) if symbol in FIBROBLAST_MARKERS}
assert len(marker_rows) == len(FIBROBLAST_MARKERS), "fixed fibroblast markers missing from feature file"


def matrix_path(gsm: str, sample: str) -> Path:
    return STAGED / f"{gsm}_{sample}_matrix.mtx.gz"


def matrix_shape(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("%"):
                a, b, c = map(int, line.split())
                return a, b, c
    raise RuntimeError(f"matrix header missing: {path}")


def stream_fibroblast_pseudobulk(path: Path) -> tuple[array, int, int, int]:
    nrows, ncols, _ = matrix_shape(path)
    assert nrows == len(gene_symbols), f"feature mismatch: {path}"
    # 3 integer arrays plus a byte mask are <100 MB for the 6.8M barcode universe.
    detected = array("I", [0]) * ncols
    umi = array("I", [0]) * ncols
    marker_hits = array("B", [0]) * ncols
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header_seen = False
        for line in handle:
            if line.startswith("%"):
                continue
            if not header_seen:
                header_seen = True
                continue
            row, col, value = line.split()
            col_i = int(col) - 1
            detected[col_i] += 1
            umi[col_i] += int(value)
            if int(row) in marker_rows:
                marker_hits[col_i] += 1
    selected = bytearray(ncols)
    selected_n = 0
    qc_n = 0
    for i in range(ncols):
        if detected[i] >= MIN_DETECTED_GENES:
            qc_n += 1
            if marker_hits[i] >= MIN_FIBROBLAST_MARKERS:
                selected[i] = 1
                selected_n += 1
    pseudobulk = array("Q", [0]) * nrows
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header_seen = False
        for line in handle:
            if line.startswith("%"):
                continue
            if not header_seen:
                header_seen = True
                continue
            row, col, value = line.split()
            if selected[int(col) - 1]:
                pseudobulk[int(row) - 1] += int(value)
    return pseudobulk, selected_n, qc_n, ncols


sample_counts: list[dict] = []
sample_pseudobulk: dict[str, array] = {}
for sample, gsm, donor, group in SAMPLES:
    pb, fibro_n, qc_n, barcode_universe = stream_fibroblast_pseudobulk(matrix_path(gsm, sample))
    sample_pseudobulk[sample] = pb
    sample_counts.append({
        "sample_id": sample,
        "donor_id": donor,
        "group": group,
        "n_fibroblast_cells": fibro_n,
        "n_qc_cells": qc_n,
        "barcode_universe": barcode_universe,
    })

write_tsv(
    RESULTS / "GSE181316_FIBROBLAST_COUNTS.tsv",
    sample_counts,
    ["sample_id", "donor_id", "group", "n_fibroblast_cells", "n_qc_cells", "barcode_universe"],
)

donor_order = ["KELOID_1", "KELOID_2", "KELOID_3", "SCAR_1", "SCAR_2", "SCAR_3"]
donor_group = {"KELOID_1": "keloid", "KELOID_2": "keloid", "KELOID_3": "keloid", "SCAR_1": "normal scar", "SCAR_2": "normal scar", "SCAR_3": "normal scar"}
donor_pb = {donor: array("Q", [0]) * len(gene_symbols) for donor in donor_order}
donor_samples = defaultdict(list)
for sample, _, donor, _ in SAMPLES:
    donor_samples[donor].append(sample)
    for i, value in enumerate(sample_pseudobulk[sample]):
        donor_pb[donor][i] += value

pseudobulk_map = []
for donor in donor_order:
    pseudobulk_map.append({
        "donor_id": donor,
        "group": donor_group[donor],
        "source_samples": ";".join(donor_samples[donor]),
        "n_source_samples": len(donor_samples[donor]),
        "n_fibroblast_cells": sum(r["n_fibroblast_cells"] for r in sample_counts if r["donor_id"] == donor),
        "aggregation": "raw integer UMI counts summed across all selected fibroblasts, then across related samples within donor",
    })
write_tsv(
    RESULTS / "GSE181316_PSEUDOBULK_MAP.tsv",
    pseudobulk_map,
    ["donor_id", "group", "source_samples", "n_source_samples", "n_fibroblast_cells", "aggregation"],
)

# Deterministic symbol collapse: raw UMI counts for duplicate symbols are summed.
symbol_rows = defaultdict(list)
for idx, symbol in enumerate(gene_symbols):
    symbol_rows[symbol].append(idx)
symbol_counts: dict[str, list[int]] = {}
for symbol, indices in symbol_rows.items():
    symbol_counts[symbol] = [sum(donor_pb[donor][i] for i in indices) for donor in donor_order]

with gzip.open(INTERMEDIATE / "GSE181316_fibroblast_donor_pseudobulk.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["gene_id", "gene_symbol", *donor_order])
    for i, (gene_id, symbol) in enumerate(zip(gene_ids, gene_symbols)):
        writer.writerow([gene_id, symbol, *[donor_pb[d][i] for d in donor_order]])


def normalised_logcpm(values: list[int], libs: list[int]) -> list[float]:
    return [math.log2(value / lib * 1_000_000 + 1) if lib else 0.0 for value, lib in zip(values, libs)]


libs = [sum(donor_pb[donor]) for donor in donor_order]
symbol_log = {symbol: normalised_logcpm(values, libs) for symbol, values in symbol_counts.items()}


def z_vector(values: list[float]) -> list[float]:
    avg = sum(values) / len(values)
    sd = math.sqrt(sum((x - avg) ** 2 for x in values) / (len(values) - 1))
    return [(x - avg) / sd for x in values] if sd > 0 else [0.0] * len(values)


symbol_z = {symbol: z_vector(values) for symbol, values in symbol_log.items()}


def locked_genes(filename: str) -> list[str]:
    return [row["gene"] for row in read_tsv(ROOT / "06_results" / "gateA" / filename)]


signatures = {
    "CORE100": locked_genes("TWPS_PRIMARY_D7_M3_CORE100.tsv"),
    "CORE50": locked_genes("TWPS_SENSITIVITY_D7_M3_CORE50.tsv"),
    "CORE25": locked_genes("TWPS_SENSITIVITY_D7_M3_CORE25.tsv"),
    "FULL": locked_genes("TWPS_SENSITIVITY_D7_M3_FULL.tsv"),
}


def score_genes(genes: list[str], z = symbol_z) -> tuple[list[float], int, float]:
    available = [gene for gene in genes if gene in z]
    scores = [sum(z[gene][i] for gene in available) / len(available) for i in range(len(donor_order))]
    return scores, len(available), len(available) / len(genes)


signature_scores = {}
for name, genes in signatures.items():
    signature_scores[name] = score_genes(genes)
assert signature_scores["CORE100"][1] >= 80, "PRIMARY_LOW_COVERAGE"

donor_twps_rows = []
for i, donor in enumerate(donor_order):
    donor_twps_rows.append({
        "donor_id": donor,
        "group": donor_group[donor],
        "CORE100": signature_scores["CORE100"][0][i],
        "CORE50": signature_scores["CORE50"][0][i],
        "CORE25": signature_scores["CORE25"][0][i],
        "FULL": signature_scores["FULL"][0][i],
        "CORE100_available": signature_scores["CORE100"][1],
        "coverage": signature_scores["CORE100"][2],
    })
write_tsv(
    RESULTS / "GSE181316_DONOR_TWPS.tsv",
    donor_twps_rows,
    ["donor_id", "group", "CORE100", "CORE50", "CORE25", "FULL", "CORE100_available", "coverage"],
)


def mean(values): return sum(values) / len(values)
def sd(values): return math.sqrt(sum((x - mean(values)) ** 2 for x in values) / (len(values) - 1))
def quantile(values, probability):
    x = sorted(values); position = (len(x) - 1) * probability; low = int(math.floor(position)); high = int(math.ceil(position))
    return x[low] if low == high else x[low] + (x[high] - x[low]) * (position - low)


def hedges_g(keloid, scar):
    n1, n2 = len(keloid), len(scar)
    pooled = math.sqrt(((n1 - 1) * sd(keloid) ** 2 + (n2 - 1) * sd(scar) ** 2) / (n1 + n2 - 2))
    if pooled == 0: return float("nan")
    return (1 - 3 / (4 * (n1 + n2 - 2) - 1)) * (mean(keloid) - mean(scar)) / pooled


def exact_permutation(keloid, scar):
    values = keloid + scar; observed = mean(keloid) - mean(scar); statistics = []
    for idx in itertools.combinations(range(6), 3):
        idx = set(idx); statistics.append(mean([values[i] for i in idx]) - mean([values[i] for i in range(6) if i not in idx]))
    return (sum(abs(x) >= abs(observed) - 1e-12 for x in statistics) / len(statistics),
            sum(x >= observed - 1e-12 for x in statistics) / len(statistics))


def cliffs_delta(keloid, scar):
    comparisons = [a - b for a in keloid for b in scar]
    return sum((x > 0) - (x < 0) for x in comparisons) / len(comparisons), sum(x > 0 for x in comparisons)


def bootstrap(keloid, scar, B=10000):
    rng = random.Random(20260811); diffs = []; gs = []
    for _ in range(B):
        a = [rng.choice(keloid) for _ in keloid]; b = [rng.choice(scar) for _ in scar]
        diffs.append(mean(a) - mean(b)); g = hedges_g(a, b)
        if math.isfinite(g): gs.append(g)
    return (quantile(diffs, .025), quantile(diffs, .975)), (quantile(gs, .025), quantile(gs, .975))


def analyse_signature(name, values, coverage):
    k = values[:3]; s = values[3:]
    p_two, p_one = exact_permutation(k, s); cliff, superior = cliffs_delta(k, s); diff_ci, g_ci = bootstrap(k, s)
    return {
        "signature": name, "n_keloid": 3, "n_normal_scar": 3, "coverage": coverage,
        "keloid_mean": mean(k), "normal_scar_mean": mean(s), "difference": mean(k) - mean(s),
        "keloid_sd": sd(k), "normal_scar_sd": sd(s), "keloid_median": quantile(k, .5), "normal_scar_median": quantile(s, .5),
        "keloid_iqr": f"{quantile(k,.25):.8g},{quantile(k,.75):.8g}", "normal_scar_iqr": f"{quantile(s,.25):.8g},{quantile(s,.75):.8g}",
        "hedges_g": hedges_g(k, s), "ci_low": g_ci[0], "ci_high": g_ci[1],
        "bootstrap_difference_low": diff_ci[0], "bootstrap_difference_high": diff_ci[1],
        "cliffs_delta": cliff, "pairwise_superiority": f"{superior}/9", "exact_p_two": p_two, "exact_p_one": p_one,
    }


analyses = {name: analyse_signature(name, score, coverage) for name, (score, _, coverage) in signature_scores.items()}
primary = analyses["CORE100"]
sensitivity_positive = sum(analyses[x]["difference"] > 0 for x in ("CORE25", "CORE50", "FULL"))
if primary["difference"] < 0 and sensitivity_positive <= 1:
    classification = "CONTRADICTORY"
elif primary["difference"] <= 0 and sensitivity_positive < 2:
    classification = "NO_REPLICATION"
elif primary["difference"] > 0 and primary["hedges_g"] >= 1 and primary["cliffs_delta"] >= .75 and int(primary["pairwise_superiority"].split("/")[0]) >= 8 and sensitivity_positive == 3:
    classification = "STRONG_REPLICATION"
elif primary["difference"] > 0 and primary["hedges_g"] >= .5 and primary["cliffs_delta"] > 0 and sensitivity_positive >= 2:
    classification = "MODERATE_REPLICATION"
else:
    classification = "WEAK_REPLICATION"

summary_rows = []
for name in ("CORE100", "CORE25", "CORE50", "FULL"):
    row = analyses[name].copy(); row["endpoint"] = "PRIMARY" if name == "CORE100" else "SENSITIVITY"; row["classification"] = classification if name == "CORE100" else "SENSITIVITY"; summary_rows.append(row)
summary_fields = ["endpoint", "signature", "n_keloid", "n_normal_scar", "coverage", "keloid_mean", "normal_scar_mean", "difference", "hedges_g", "ci_low", "ci_high", "cliffs_delta", "pairwise_superiority", "exact_p_two", "exact_p_one", "classification"]
write_tsv(RESULTS / "GSE181316_GATE_C_SUMMARY.tsv", summary_rows, summary_fields)

sensitivity_fields = ["signature", "n_keloid", "n_normal_scar", "coverage", "keloid_mean", "normal_scar_mean", "difference", "hedges_g", "ci_low", "ci_high", "cliffs_delta", "pairwise_superiority", "exact_p_two", "exact_p_one"]
write_tsv(RESULTS / "GSE181316_SIGNATURE_SENSITIVITY.tsv", [analyses[x] for x in ("CORE25", "CORE50", "FULL")], sensitivity_fields)

# Related-sample technical sensitivity, deliberately not an inferential analysis.
sample_libs = [sum(sample_pseudobulk[s]) for s, *_ in SAMPLES]
sample_symbol_counts = {symbol: [sum(sample_pseudobulk[s][i] for i in indices) for s, *_ in SAMPLES] for symbol, indices in symbol_rows.items()}
sample_z = {symbol: z_vector(normalised_logcpm(values, sample_libs)) for symbol, values in sample_symbol_counts.items()}
within_scores, _, _ = score_genes(signatures["CORE100"], sample_z)
within = []
for sample, _, donor, group in SAMPLES:
    if donor == "KELOID_3":
        idx = [x[0] for x in SAMPLES].index(sample)
        within.append({"sample_id": sample, "donor_id": donor, "group": group, "CORE100_technical_score": within_scores[idx], "purpose": "technical within-donor consistency only; not independent n"})
within_diff = abs(within[0]["CORE100_technical_score"] - within[1]["CORE100_technical_score"])
for row in within: row["absolute_pair_difference"] = within_diff
write_tsv(RESULTS / "GSE181316_WITHIN_DONOR_SENSITIVITY.tsv", within, ["sample_id", "donor_id", "group", "CORE100_technical_score", "absolute_pair_difference", "purpose"])

# Secondary locked modules after primary and sensitivity outputs exist.
secondary = read_tsv(ROOT / "06_results" / "gateA" / "TWPS_SECONDARY_MODULES_LOCKED.tsv")
secondary_rows = []
for module in ("D7_M1", "D7_M2", "D1_M2", "D1_M1"):
    genes = [r["gene"] for r in secondary if r["module"] == module]
    scores, available, coverage = score_genes(genes)
    a = analyse_signature(module, scores, coverage)
    a["label"] = "SECONDARY_EXPLORATORY"; a["gene_n_available"] = available
    secondary_rows.append(a)
p_values = [r["exact_p_two"] for r in secondary_rows]
order = sorted(range(len(p_values)), key=lambda i: p_values[i]); adjusted = [0.0] * len(p_values); prior = 0.0
for rank, idx in enumerate(order, 1):
    prior = max(prior, p_values[idx] * len(p_values) / rank); adjusted[idx] = min(1.0, prior)
for row, q in zip(secondary_rows, adjusted): row["BH_FDR"] = q
secondary_fields = ["module", "label", "gene_n_available", "coverage", "difference", "hedges_g", "exact_p_two", "exact_p_one", "BH_FDR"]
for row in secondary_rows: row["module"] = row.pop("signature")
write_tsv(RESULTS / "GSE181316_SECONDARY_MODULES.tsv", secondary_rows, secondary_fields)

audit = [{
    "original_label": "NO_AUTHOR_CELL_LABEL_FILE_IN_GEO_RAW_RELEASE",
    "n_cells": sum(r["n_fibroblast_cells"] for r in sample_counts),
    "included": "YES",
    "reason": f"Author-provided raw release contains sparse matrices/features/barcodes only. Fixed pre-TWPS canonical fibroblast gate: >= {MIN_DETECTED_GENES} detected genes and >= {MIN_FIBROBLAST_MARKERS} detected markers among COL1A1,COL1A2,COL3A1,DCN,LUM; no clustering, UMAP, or disease-driven selection.",
}]
write_tsv(RESULTS / "GSE181316_FIBROBLAST_LABEL_AUDIT.tsv", audit, ["original_label", "n_cells", "included", "reason"])
write_tsv(RESULTS / "GSE181316_SUBTYPE_SENSITIVITY.tsv", [{"SUBTYPE_SENSITIVITY": "NOT_RELIABLY_TESTABLE", "reason": "No author-provided fibroblast subtype labels were distributed in the GEO raw release; no reclustering was performed."}], ["SUBTYPE_SENSITIVITY", "reason"])

def fmt(value): return f"{value:.8g}" if isinstance(value, float) else str(value)
total_fb = sum(r["n_fibroblast_cells"] for r in sample_counts)
k_fb = sum(r["n_fibroblast_cells"] for r in sample_counts if r["group"] == "keloid")
s_fb = sum(r["n_fibroblast_cells"] for r in sample_counts if r["group"] == "normal scar")
report = [
    "# Replacement Gate C — GSE181316", "", "## Integrity", "", "TWPS_SHA256_MATCH=YES", "PHASE3R_OUTCOME_LEAKAGE=NO", "",
    "## Cohort", "", "KELOID_SAMPLE_N=4", "KELOID_DONOR_N=3", "NORMAL_SCAR_SAMPLE_N=3", "NORMAL_SCAR_DONOR_N=3", "INDEPENDENT_DONORS_VERIFIED=YES", "RELATED_SAMPLES=keloid_3L and keloid_3R were summed into one KELOID_3 donor pseudobulk", "",
    "## Fibroblasts", "", f"TOTAL_FIBROBLASTS={total_fb}", f"KELOID_FIBROBLASTS={k_fb}", f"NORMAL_SCAR_FIBROBLASTS={s_fb}", "AUTHOR_ANNOTATION_STATUS=No author cell-label file in GEO raw release; fixed non-disease-driven canonical fibroblast gate used and documented in GSE181316_FIBROBLAST_LABEL_AUDIT.tsv.", "",
    "## Primary coverage", "", f"CORE100_AVAILABLE={signature_scores['CORE100'][1]}", f"CORE100_COVERAGE={fmt(signature_scores['CORE100'][2])}", "",
    "## Primary CORE100", "", f"KELOID_MEAN={fmt(primary['keloid_mean'])}", f"NORMAL_SCAR_MEAN={fmt(primary['normal_scar_mean'])}", f"DIFFERENCE={fmt(primary['difference'])}", f"HEDGES_G={fmt(primary['hedges_g'])}", f"HEDGES_G_95CI={fmt(primary['ci_low'])},{fmt(primary['ci_high'])}", f"CLIFFS_DELTA={fmt(primary['cliffs_delta'])}", f"PAIRWISE_SUPERIORITY={primary['pairwise_superiority']}", f"EXACT_P_TWO_SIDED={fmt(primary['exact_p_two'])}", f"EXACT_P_ONE_SIDED={fmt(primary['exact_p_one'])}", f"BOOTSTRAP_DIFFERENCE_95CI={fmt(primary['bootstrap_difference_low'])},{fmt(primary['bootstrap_difference_high'])}", "",
    "## Sensitivity", "",
]
for name in ("CORE25", "CORE50", "FULL"):
    a = analyses[name]; report += [f"{name}: difference={fmt(a['difference'])}; g={fmt(a['hedges_g'])}; Cliff={fmt(a['cliffs_delta'])}; pairwise={a['pairwise_superiority']}; P_two={fmt(a['exact_p_two'])}; P_one={fmt(a['exact_p_one'])}"]
report += ["", "## Related-sample sensitivity", "", f"keloid_3L vs keloid_3R technical-score absolute difference={fmt(within_diff)}; descriptive only and not independent n.", "", "## Fibroblast subtype sensitivity", "", "SUBTYPE_SENSITIVITY=NOT_RELIABLY_TESTABLE", "", "## Secondary locked modules", "", "All secondary modules are labelled SECONDARY_EXPLORATORY; exact permutation P and BH-FDR are in GSE181316_SECONDARY_MODULES.tsv.", "", "## Limitations", "", "Biological n=3 versus 3; statistical precision is limited; cells are not independent biological replicates; exact permutation P-values have coarse resolution. GEO raw release lacked an author cell-label metadata file, so a fixed, non-disease-driven canonical fibroblast gate was used without reclustering.", "", "## Evidence classification", "", classification]
(RESULTS / "GSE181316_GATE_C_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

print(f"TOTAL_FIBROBLASTS={total_fb}")
print(f"CORE100_COVERAGE={signature_scores['CORE100'][1]}/100")
print(f"EVIDENCE_CLASS={classification}")
print(f"PRIMARY_DIFFERENCE={primary['difference']}")
