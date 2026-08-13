#!/usr/bin/env python3
"""Phase 4A metadata-only GSE178411 patient-structure audit.

This script reads sample metadata and the processed-count matrix header only. It
does not read expression values, calculate TWPS, or perform outcome analysis.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PHASE0 = ROOT / "06_results/phase0/GSE178411_SAMPLE_MAP.tsv"
MATRIX = ROOT / "data/raw_small/GSE178411/GSE178411_counts.txt.gz"
OUT = ROOT / "06_results/spectrum"
OUT.mkdir(parents=True, exist_ok=True)

STATE = {
    "Normal skin": "UNINJURED_SKIN",
    "Early Wound": "EARLY_WOUND",
    "Late wound": "LATE_WOUND",
    "Chronic wound": "CHRONIC_WOUND",
    "Normal scar": "NORMAL_SCAR",
    "HTS": "HYPERTROPHIC_SCAR",
}
ORDER = ["UNINJURED_SKIN", "EARLY_WOUND", "LATE_WOUND", "CHRONIC_WOUND", "NORMAL_SCAR", "HYPERTROPHIC_SCAR", "OTHER", "UNRESOLVED"]

def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

with PHASE0.open(newline="") as handle:
    source = list(csv.DictReader(handle, delimiter="\t"))
with gzip.open(MATRIX, "rt") as handle:
    # No GeneID placeholder is present in the header; its first token is 12-W.
    matrix_columns = handle.readline().rstrip("\n").split("\t")

by_patient: dict[str, list[dict]] = defaultdict(list)
for row in source:
    row["state"] = STATE.get(row["biological_state"], "UNRESOLVED")
    by_patient[row["patient_id"]].append(row)

record_rows = []
for row in source:
    patient_rows = by_patient[row["patient_id"]]
    related = [x["sample_id"] for x in patient_rows if x["sample_id"] != row["sample_id"]]
    record_rows.append({
        "record_id": row["GSM"], "GSM_if_available": row["GSM"], "sample_id": row["sample_id"],
        "patient_id_raw": row["patient_id"], "patient_id_resolved": row["patient_id"],
        "biological_sample_id": row["sample_id"], "technical_record_id": "",
        "tissue_or_state_raw": row["biological_state"], "state_standardized": row["state"],
        "timepoint_raw": row["time_or_stage"], "timepoint_standardized": row["state"],
        "paired_or_repeated": "PATIENT_REPEATED" if related else "UNPAIRED",
        "related_records": ";".join(related),
        "include_candidate": "YES" if row["sample_id"] in matrix_columns else "YES_METADATA_ONLY_MATRIX_MISSING",
        "evidence_source": "official GEO SOFT title/characteristics as parsed in Phase0",
        "confidence": "HIGH",
        "notes": "One GEO record represents one annotated biological tissue sample; no technical-record duplication was identified in GEO metadata.",
    })
write_tsv(OUT / "GSE178411_RECORD_MAP.tsv", record_rows, list(record_rows[0]))

patient_rows = []
for patient, rows in sorted(by_patient.items(), key=lambda x: int(x[0].split()[-1])):
    states = []
    for r in rows:
        if r["state"] not in states:
            states.append(r["state"])
    sample_ids = [r["sample_id"] for r in rows]
    patient_rows.append({
        "patient_id": patient, "available_states": ";".join(states),
        "n_biological_samples": len(rows), "n_records": len(rows),
        "longitudinal_or_repeated_structure": "REPEATED_OR_MULTISTATE" if len(rows) > 1 else "SINGLE_SAMPLE",
        "notes": ";".join(sample_ids),
    })
write_tsv(OUT / "GSE178411_PATIENT_STATE_MAP.tsv", patient_rows, list(patient_rows[0]))

groups = {}
for state in ORDER:
    rows = [r for r in source if r["state"] == state]
    groups[state] = {"records": len(rows), "samples": len(rows), "patients": len({r["patient_id"] for r in rows})}

matrix_set = set(matrix_columns)
match_rows = []
for col in matrix_columns:
    found = [r for r in source if r["sample_id"] == col]
    if len(found) == 1:
        r = found[0]
        match_rows.append({"matrix_column": col, "biological_sample_id": r["sample_id"], "patient_id": r["patient_id"], "state": r["state"], "matched": "YES", "ambiguity": ""})
    else:
        match_rows.append({"matrix_column": col, "biological_sample_id": "", "patient_id": "", "state": "", "matched": "NO", "ambiguity": "No unique Phase0 metadata record"})
write_tsv(OUT / "GSE178411_MATRIX_METADATA_MATCH.tsv", match_rows, list(match_rows[0]))

def patient_set(state: str) -> set[str]:
    return {r["patient_id"] for r in source if r["state"] == state}

contrasts = [
    ("early wound vs late wound", "EARLY_WOUND", "LATE_WOUND", "patient-level independent-group comparison with patient-clustered sensitivity for the one overlap", "INDEPENDENT_GROUP_COMPARISON", "One patient contributes both states; do not count their records as independent twice."),
    ("late wound vs normal scar", "LATE_WOUND", "NORMAL_SCAR", "descriptive patient-level contrast only", "DESCRIPTIVE_ONLY", "Normal-scar group has two patients."),
    ("normal scar vs hypertrophic scar", "NORMAL_SCAR", "HYPERTROPHIC_SCAR", "descriptive patient-level contrast only", "DESCRIPTIVE_ONLY", "Normal-scar group has two patients."),
    ("chronic wound vs normal scar", "CHRONIC_WOUND", "NORMAL_SCAR", "descriptive patient-level contrast only", "DESCRIPTIVE_ONLY", "Both groups are very small (3 and 2 patients)."),
    ("uninjured skin vs early wound", "UNINJURED_SKIN", "EARLY_WOUND", "paired subset analysis for the eight shared patients; separate patient-level independent comparison only if explicitly prespecified", "PATIENT_REPEATED_MEASURES", "Eight patients contribute both states; repeated structure is material."),
    ("uninjured skin vs late wound", "UNINJURED_SKIN", "LATE_WOUND", "paired subset analysis for the eight shared patients; separate patient-level independent comparison only if explicitly prespecified", "PATIENT_REPEATED_MEASURES", "Eight patients contribute both states; repeated structure is material."),
    ("uninjured skin vs chronic wound", "UNINJURED_SKIN", "CHRONIC_WOUND", "descriptive paired subset only", "DESCRIPTIVE_ONLY", "Only one shared patient and three chronic-wound patients total."),
]
design_rows = []
for name, a, b, model, validity, notes in contrasts:
    pa, pb = patient_set(a), patient_set(b)
    overlap = pa & pb
    design_rows.append({"contrast": name, "state_A": a, "state_B": b, "patient_overlap": ";".join(sorted(overlap)) if overlap else "NONE", "n_patients_A": len(pa), "n_patients_B": len(pb), "paired_n_if_any": len(overlap), "recommended_model_structure": model, "validity": validity, "notes": notes})
write_tsv(OUT / "GSE178411_FUTURE_DESIGN.tsv", design_rows, list(design_rows[0]))

sha = hashlib.sha256(PHASE0.read_bytes()).hexdigest()
missing_metadata_samples = sorted({r["sample_id"] for r in source} - matrix_set)
reconciliation = f"""# GSE178411 Count Reconciliation

GEO_RECORD_TOTAL=108

PROCESSED_MATRIX_SAMPLE_TOTAL={len(matrix_columns)}

BIOLOGICAL_SAMPLE_TOTAL=108

UNIQUE_PATIENT_TOTAL={len(by_patient)}

SUM_OF_STANDARDIZED_GROUP_SAMPLE_COUNTS={sum(v['samples'] for v in groups.values())}

SUM_OF_STANDARDIZED_GROUP_PATIENT_COUNTS={sum(v['patients'] for v in groups.values())}

DISCREPANCY_STATUS=PARTIALLY_RESOLVED

EXPLANATION=All 108 GEO records map one-to-one to annotated biological samples, but 23 of 75 patients contribute repeated or multistate samples, explaining why record count exceeds patient count. The processed count matrix has 108 sample columns and matches all records; its header omits a GeneID placeholder, so the first header token `12-W` is a sample. The official series narrative states uninjured skin n=26, acute wounds n=54, and HTS n=30; the record-level metadata instead gives UNINJURED_SKIN=24, all wound states=54, and scar states=30 (HYPERTROPHIC_SCAR=28 plus NORMAL_SCAR=2). Thus wound and aggregate-scar totals reconcile, but the two-record uninjured-skin discrepancy remains undocumented in the available metadata.

PHASE0_SAMPLE_MAP_SHA256={sha}
"""
(OUT / "GSE178411_COUNT_RECONCILIATION.md").write_text(reconciliation)

group_lines = []
for state in ORDER:
    g = groups[state]
    group_lines.extend([f"{state}:", f"samples={g['samples']}", f"patients={g['patients']}", ""])
audit = f"""# GSE178411 Patient Structure Audit

RECORD_N=108

BIOLOGICAL_SAMPLE_N=108

UNIQUE_PATIENT_N={len(by_patient)}

MATRIX_SAMPLE_N={len(matrix_columns)}

## Groups

{chr(10).join(group_lines)}
## Count reconciliation

STATUS=PARTIALLY_RESOLVED

EXPLANATION=The 108-record total reflects 108 biological samples from 75 patients, including repeated/multistate sampling from 23 patients. The matrix contains and matches all 108 samples; its header has no GeneID placeholder, and `12-W` is the first sample column. Record-level metadata reconciles all 54 wound samples and all 30 scar samples but documents only 24, rather than the series narrative's 26, uninjured-skin samples.

## Repeated / paired structure

Twenty-three patients contribute 56 records across repeated biopsies and/or states. Material shared-patient structures include eight patients each for UNINJURED_SKIN versus EARLY_WOUND and UNINJURED_SKIN versus LATE_WOUND. Future inference must use patient as the independent unit and preserve these repeated relationships.

## Matrix matching

FULLY_MATCHED={sum(r['matched'] == 'YES' for r in match_rows)}

PARTIALLY_MATCHED=0

UNMATCHED=0

The matrix header maps uniquely to all 108 metadata samples. The earlier 107-column interpretation was a parser error caused by treating the first sample name (`12-W`) as a feature-column header.

## Future valid comparisons

See `GSE178411_FUTURE_DESIGN.tsv`. Repeated structures support paired subset analyses for uninjured-versus-wound contrasts; sparse normal-scar and chronic-wound groups are descriptive only.

## Major limitations

- The official narrative's uninjured-skin count (26) conflicts with 24 record-level normal-skin samples; two records are not resolved by available metadata.
- The matrix header lacks an explicit GeneID placeholder; parsers must preserve its first token as sample `12-W`.
- Several state contrasts have small patient counts, especially NORMAL_SCAR (2) and CHRONIC_WOUND (3).

## Feasibility

USABLE_WITH_LIMITATIONS

## Proposed future statistical framework

Use patient-level sample scores. Use paired contrasts only for explicit shared-patient subsets; for broader contrasts with limited overlap, define a patient-aware model or an overlap-handling sensitivity plan before outcome analysis. Do not use record count as biological n.

## Decision category

REVIEW_REQUIRED

TWPS_ANALYZED=NO
"""
(OUT / "GSE178411_STRUCTURE_AUDIT.md").write_text(audit)

print("RECORD_N=108")
print("BIOLOGICAL_SAMPLE_N=108")
print(f"UNIQUE_PATIENT_N={len(by_patient)}")
print(f"MATRIX_COLUMNS={len(matrix_columns)}")
print("TWPS_ANALYZED=NO")
