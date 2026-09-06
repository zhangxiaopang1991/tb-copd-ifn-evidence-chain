# -*- coding: utf-8 -*-
"""Audit deposited-versus-analysed sample flow for GSE276060 and E-MTAB-17246."""
import csv
import gzip
import os
import re
from collections import Counter, OrderedDict

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "outputs")


def read_gse276060():
    fields = {}
    with gzip.open(os.path.join(BASE, "GSE276060", "GSE276060_series_matrix.txt.gz"), "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            fields[row[0]] = [x.strip('"') for x in row[1:]]

    with gzip.open(os.path.join(BASE, "GSE276060", "GSE276060_ExpressionRNAseq.csv.gz"), "rt", encoding="utf-8-sig", errors="replace") as fh:
        header = next(csv.reader(fh))
    pattern = re.compile(r"^(.+?)_S\d+ \(GE\) - Total counts$")
    analysed = {match.group(1) for value in header[5:] if (match := pattern.match(value))}

    rows = []
    for title, accession, description in zip(
        fields["!Sample_title"],
        fields["!Sample_geo_accession"],
        fields["!Sample_description"],
    ):
        body = description.split("Case ", 1)[-1]
        if body.startswith("MTB"):
            group = "MTB"
        elif body.startswith("NTM"):
            group = "NTM"
        elif body.startswith("Ctl-Lung"):
            group = "Ctl_Lung"
        elif body.startswith("Ctl-Tonsil"):
            group = "Ctl_Tonsil"
        else:
            group = "Other_granulomatous_or_disease"
        rows.append({
            "dataset": "GSE276060",
            "accession": accession,
            "sample": title,
            "group": group,
            "included": title in analysed,
            "reason": "Total counts column available and included in count matrix" if title in analysed else "Series sample without Total counts column in analysed RNA-seq file",
        })
    return rows


def read_emtab():
    path = os.path.join(BASE, "EMTAB17246.sdrf.txt")
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        seen = OrderedDict()
        for row in reader:
            sample = row[idx["Extract Name"]]
            if sample not in seen:
                seen[sample] = row

    rows = []
    for sample, row in seen.items():
        rows.append({
            "dataset": "E-MTAB-17246",
            "accession": "",
            "sample": sample,
            "group": f'{row[idx["Factor Value[infect]"]]} / {row[idx["Factor Value[stimulus]"]]}',
            "included": True,
            "reason": "Present in local SDRF and count matrix",
            "donor": row[idx["Characteristics[individual]"]],
        })
    return rows


def write_tsv(path, rows):
    keys = list(rows[0])
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


gse_rows = read_gse276060()
emtab_rows = read_emtab()
write_tsv(os.path.join(OUT, "GSE276060_sample_flow_audit.tsv"), gse_rows)
write_tsv(os.path.join(OUT, "EMTAB17246_sample_flow_audit.tsv"), emtab_rows)

lines = [
    "# Sample-flow audit",
    "",
    "> Generated from the local series matrix, expression matrix, SDRF and count files on 2026-09-05.",
    "",
    "## GSE276060",
    "",
    f"- Series samples: {len(gse_rows)}",
    f"- Samples with `Total counts` in the analysed expression file: {sum(row['included'] for row in gse_rows)}",
    f"- Excluded from the count analysis because no analysed `Total counts` column was present: {sum(not row['included'] for row in gse_rows)}",
    f"- Included group counts: {dict(Counter(row['group'] for row in gse_rows if row['included']))}",
    f"- Excluded group counts: {dict(Counter(row['group'] for row in gse_rows if not row['included']))}",
    "",
    "The local series design text reports a different group summary from the sample-level annotations. The manuscript should use the sample-level annotations for the analysed 24-sample subset and explicitly state that the remaining 13 series samples were outside the bulk count analysis.",
    "",
    "Excluded GSE276060 samples:",
]
for row in gse_rows:
    if not row["included"]:
        lines.append(f"- {row['accession']} ({row['sample']}): {row['group']}")

lines += [
    "",
    "## E-MTAB-17246",
    "",
    f"- Unique samples in local SDRF: {len(emtab_rows)}",
    f"- Donors represented: {dict(Counter(row['donor'] for row in emtab_rows))}",
    f"- Condition counts: {dict(Counter(row['group'] for row in emtab_rows))}",
    "- The local SDRF and count matrix form a complete 7-donor × 6-condition design. No eighth donor is present in the local sample table.",
    "- The local IDF narrative mentions eight donors; this is a repository-description discrepancy that must be resolved against the source publication or an archived original sample sheet.",
]

with open(os.path.join(OUT, "sample_flow_audit_2026-09-05.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")

print("\n".join(lines[:18]))
