#!/usr/bin/env python3
"""Mechanical discovery-only GO BP/Reactome enrichment for fixed Gate A2 modules."""
from pathlib import Path
import csv
import json
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "06_results" / "gateA"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


robust = {r["gene"]: r for r in read_tsv(RES / "A2_GENE_ROBUSTNESS.tsv")}
membership = read_tsv(RES / "A2_MODULE_MEMBERSHIP.tsv")
metrics = read_tsv(RES / "A2_MODULE_METRICS.tsv")
lodo = read_tsv(RES / "A2_MODULE_LODO.tsv")
ranking = read_tsv(RES / "A2_MODULE_CORE_RANKING.tsv")

module_genes = {}
for row in membership:
    module = row["module"]
    if module == "NOT_CLUSTERED":
        continue
    symbol = robust[row["gene"]].get("gene_symbol") or row["gene"]
    module_genes.setdefault(module, []).append(symbol)


def post_json(url, fields):
    boundary = "----TWPSGateA2Boundary"
    parts = []
    for name, value in fields.items():
        parts.extend([
            f"--{boundary}\r\n",
            f'Content-Disposition: form-data; name="{name}"\r\n\r\n',
            str(value),
            "\r\n",
        ])
    parts.append(f"--{boundary}--\r\n")
    data = "".join(parts).encode()
    req = urllib.request.Request(url, data=data, headers={"Content-Type": f"multipart/form-data; boundary={boundary}", "User-Agent": "TWPS-GateA2/1.0"})
    with urllib.request.urlopen(req, timeout=60) as response:
        return json.load(response)


libraries = ["GO_Biological_Process_2025", "Reactome_Pathways_2024"]
enrichment = []
for module, genes in sorted(module_genes.items()):
    payload = post_json("https://maayanlab.cloud/Enrichr/addList", {
        "list": "\n".join(sorted(set(genes))),
        "description": f"{module} discovery-only fixed module",
    })
    user_list_id = payload["userListId"]
    for library in libraries:
        url = "https://maayanlab.cloud/Enrichr/enrich?" + urllib.parse.urlencode({
            "userListId": user_list_id,
            "backgroundType": library,
        })
        with urllib.request.urlopen(url, timeout=60) as response:
            result = json.load(response).get(library, [])
        kept = [x for x in result if float(x[6]) < 0.05][:20]
        if not kept:
            kept = result[:5]
        for x in kept:
            overlap_genes = x[5] if isinstance(x[5], list) else str(x[5]).split(";")
            enrichment.append({
                "module": module,
                "database": "GO Biological Process" if library.startswith("GO_") else "Reactome",
                "term": x[1],
                "overlap": f"{len(overlap_genes)}/{len(set(genes))}",
                "effect_or_odds_ratio_if_available": f"combined_score={float(x[4]):.6g}",
                "pvalue": f"{float(x[2]):.8g}",
                "FDR": f"{float(x[6]):.8g}",
                "overlap_genes": ";".join(overlap_genes),
            })

with (RES / "A2_MODULE_ENRICHMENT.tsv").open("w", newline="") as handle:
    fields = ["module", "database", "term", "overlap", "effect_or_odds_ratio_if_available", "pvalue", "FDR", "overlap_genes"]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(enrichment)

enrich_by_module = {}
for row in enrichment:
    if float(row["FDR"]) < 0.05:
        enrich_by_module.setdefault(row["module"], []).append(row)

lodo_by_module = {}
for row in lodo:
    lodo_by_module.setdefault(row["module"], []).append(row)

metric_by_module = {r["module"]: r for r in metrics}
quality_order = {"HIGH_QUALITY": 3, "MODERATE_QUALITY": 2, "LOW_QUALITY": 1, "FAIL": 0}
ordered = sorted(metrics, key=lambda r: (
    -quality_order[r["QUALITY"]],
    -sum(float(x["trajectory_correlation"]) for x in lodo_by_module[r["module"]]) / 3,
    -(float(r["ACTIVATION_EFFECT"]) + float(r["ATTENUATION_EFFECT"])),
    float(r["RETURN_TOWARD_SKIN_RATIO"]),
    r["module"],
))

candidate_rows = []
for i, row in enumerate(ordered, 1):
    module = row["module"]
    lr = lodo_by_module[module]
    top_terms = [x["term"] for x in enrich_by_module.get(module, [])[:3]]
    donor_concordance = f"activation {row['ACTIVATION_CONCORDANCE']}/3; attenuation {row['ATTENUATION_CONCORDANCE']}/3"
    lodo_stability = "mean trajectory r={:.3f}; activation signs {}/3; attenuation signs {}/3".format(
        sum(float(x["trajectory_correlation"]) for x in lr) / len(lr),
        sum(x["activation_sign_preserved"] == "TRUE" for x in lr),
        sum(x["attenuation_sign_preserved"] == "TRUE" for x in lr),
    )
    n = int(row["N_GENES"])
    candidate_rows.append({
        "candidate_id": f"A2_C{i}",
        "module": module,
        "full_gene_n": n,
        "core25_available": "YES" if n >= 25 else "NO",
        "core50_available": "YES" if n >= 50 else "NO",
        "core100_available": "YES" if n >= 100 else "NO",
        "quality": row["QUALITY"],
        "activation": row["ACTIVATION_EFFECT"],
        "attenuation": row["ATTENUATION_EFFECT"],
        "donor_concordance": donor_concordance,
        "LODO_stability": lodo_stability,
        "subtype_robustness": row["SUBTYPE_SENSITIVITY"],
        "top_enrichment_terms": ";".join(top_terms) if top_terms else "NO_FDR_SIGNIFICANT_TERM",
    })

with (RES / "A2_TWPS_CANDIDATES_NOT_LOCKED.tsv").open("w", newline="") as handle:
    handle.write("# STATUS=NOT_LOCKED\n")
    handle.write("# VALIDATION_DATA_VIEWED=NO\n")
    fields = list(candidate_rows[0])
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(candidate_rows)
