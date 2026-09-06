"""Create a compact, reproducible metadata QA report for the five datasets."""
from __future__ import annotations

import csv
import gzip
import re
from collections import Counter
from pathlib import Path

BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)


def read_tsv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else Path.open
    kwargs = {"encoding": "utf-8", "errors": "replace"}
    with opener(path, "rt", **kwargs) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def gse47460_audit():
    rows = read_tsv(BASE / "GSE47460" / "sample_meta_final.tsv")
    keep = [r for r in rows if r["disease state"] in {"Control", "Chronic Obstructive Lung Disease"}]
    fields = ["smoker?", "%predicted fev1 (pre-bd)", "%predicted fvc (pre-bd)", "%predicted dlco", "age", "Sex"]
    missing = {f: sum(not r[f].strip() for r in keep) for f in fields}
    disease = Counter("COPD" if r["disease state"] == "Chronic Obstructive Lung Disease" else "Control" for r in keep)
    platform = Counter(r["gpl"] for r in keep)
    smoking = Counter()
    for r in keep:
        raw = r["smoker?"].strip()
        if not raw:
            smoking["missing"] += 1
        elif "Never" in raw:
            smoking["Never"] += 1
        elif "Current" in raw:
            smoking["Current"] += 1
        else:
            smoking["Ever"] += 1
    return {
        "dataset": "GSE47460",
        "source": "GSE47460/sample_meta_final.tsv",
        "all_rows": len(rows),
        "analytic_rows": len(keep),
        "disease_counts": dict(disease),
        "platform_counts": dict(platform),
        "smoking_mapping": dict(smoking),
        "missing_counts": missing,
        "field_definition": {
            "fev1": "%predicted fev1 (pre-bd), numeric, COPD-only association models",
            "smoking": "Never if label contains Never; Current if label contains Current; other non-empty labels map to Ever",
            "age": "numeric conversion of age",
            "sex": "factor from Sex; source values 1-Male and 2-Female",
        },
        "missing_handling": "genes observed in at least 80% of analytic samples retained; remaining gene-level missing values are row-mean imputed; covariate models use complete cases",
    }


def gse192483_audit():
    samples = []
    for path in sorted((BASE / "GSE192483" / "extracted").glob("*.h5")):
        m = re.search(r"_(SP\d+)([HL])\.filtered_feature_bc_matrix$", path.stem)
        if not m:
            raise ValueError(f"Unexpected sample filename: {path.name}")
        samples.append((m.group(1), m.group(2), path.name))
    by_patient = {}
    for patient, lesion, name in samples:
        by_patient.setdefault(patient, {})[lesion] = name
    paired = sorted(p for p, x in by_patient.items() if {"H", "L"} <= set(x))
    return {
        "dataset": "GSE192483",
        "source": "GSE192483/extracted/*.h5",
        "sample_count": len(samples),
        "patient_count": len(by_patient),
        "paired_patient_count": len(paired),
        "paired_patients": paired,
        "single_sample_patients": sorted(p for p, x in by_patient.items() if len(x) == 1),
        "sample_layout": [{"patient": p, "H": x.get("H", ""), "L": x.get("L", "")} for p, x in sorted(by_patient.items())],
        "analytic_unit": "patient-level H versus L paired comparison; cells are descriptive only",
        "current_qc": "n_genes_by_counts >= 200 and pct_counts_mt < 20",
    }


def gse136831_audit():
    rows = read_tsv(BASE / "GSE136831" / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz")
    disease = Counter(r["Disease_Identity"] for r in rows)
    subjects = {r["Subject_Identity"]: r["Disease_Identity"] for r in rows}
    types = Counter(r["Manuscript_Identity"] for r in rows if r["Manuscript_Identity"])
    pb = read_tsv(BASE / "GSE136831" / "pseudobulk_counts" / "pseudobulk_summary.tsv")
    low_replication = [
        {k: r[k] for k in ("celltype", "n_pseudobulk", "n_COPD", "n_Control")}
        for r in pb if int(r["n_COPD"]) < 3 or int(r["n_Control"]) < 3
    ]
    return {
        "dataset": "GSE136831",
        "source": "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz",
        "cell_count": len(rows),
        "disease_counts": dict(disease),
        "subject_count": len(subjects),
        "subjects_by_disease": dict(Counter(subjects.values())),
        "annotation_field": "Manuscript_Identity",
        "subject_field": "Subject_Identity",
        "disease_field": "Disease_Identity",
        "celltype_count": len(types),
        "celltype_cell_counts_top": dict(types.most_common(10)),
        "pseudobulk_replication_rule": "subject x cell type aggregation; Welch t test requires at least 3 subjects per disease group",
        "low_replication_celltypes": low_replication,
        "cell_level_results_note": "cell-level tests are descriptive diagnostics; formal disease inference uses subject-level pseudobulk outputs",
    }


records = [gse47460_audit(), gse192483_audit(), gse136831_audit()]
lines = [
    "# Five-dataset metadata QA",
    "",
    "> Generated from local metadata and analysis inputs on 2026-09-05. This report records evidence and analytic rules; it does not add unavailable source-paper claims.",
    "",
]
for rec in records:
    lines += [f"## {rec['dataset']}", ""]
    for key, value in rec.items():
        if key in {"dataset", "source"}:
            continue
        if isinstance(value, dict):
            value = "; ".join(f"{k}={v}" for k, v in value.items())
        lines.append(f"- **{key}**: {value}")
    lines += [f"- **source_file**: `{rec['source']}`", ""]
(OUT / "dataset_metadata_QA_2026-09-05.md").write_text("\n".join(lines), encoding="utf-8")
with (OUT / "dataset_metadata_QA_2026-09-05.tsv").open("w", encoding="utf-8", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["dataset", "metric", "value"])
    for rec in records:
        for key, value in rec.items():
            if key not in {"dataset", "source"}:
                writer.writerow([rec["dataset"], key, repr(value)])
print("Wrote metadata QA outputs.")
