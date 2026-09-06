# -*- coding: utf-8 -*-
"""Compare cell-type proportions between COPD and Control at subject level."""
import os
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

BASE = os.path.dirname(os.path.abspath(__file__))
meta = pd.read_csv(os.path.join(BASE, "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"), sep="\t")
meta = meta[meta["Disease_Identity"].isin(["COPD", "Control"])].copy()
meta["celltype"] = meta["Manuscript_Identity"].astype(str)
meta = meta[meta["celltype"] != "Multiplet"].copy()
subject_info = meta[["Subject_Identity", "Disease_Identity"]].drop_duplicates()
all_types = sorted(meta["celltype"].unique())
grid = (
    subject_info.assign(_key=1)
    .merge(pd.DataFrame({"celltype": all_types, "_key": 1}), on="_key")
    .drop(columns="_key")
)
observed = (
    meta.groupby(["Subject_Identity", "Disease_Identity", "celltype"])
    .size()
    .rename("n_cells")
    .reset_index()
)
counts = grid.merge(
    observed,
    on=["Subject_Identity", "Disease_Identity", "celltype"],
    how="left",
)
counts["n_cells"] = counts["n_cells"].fillna(0).astype(int)
totals = (
    counts.groupby(["Subject_Identity", "Disease_Identity"])["n_cells"]
    .sum()
    .rename("subject_total")
    .reset_index()
)
counts = counts.merge(totals, on=["Subject_Identity", "Disease_Identity"])
counts["proportion"] = counts["n_cells"] / counts["subject_total"]
rows = []
for ct in all_types:
    s = counts[counts["celltype"] == ct]
    copd = s.loc[s["Disease_Identity"] == "COPD", "proportion"]
    control = s.loc[s["Disease_Identity"] == "Control", "proportion"]
    if len(copd) < 3 or len(control) < 3:
        continue
    _, p = mannwhitneyu(copd, control, alternative="two-sided")
    rows.append({"celltype": ct, "n_COPD": len(copd), "n_Control": len(control),
                 "mean_prop_COPD": copd.mean(), "mean_prop_Control": control.mean(),
                 "mean_diff_COPD_minus_Control": copd.mean() - control.mean(), "pval": p})
res = pd.DataFrame(rows)
res["padj"] = multipletests(res["pval"], method="fdr_bh")[1]
counts.to_csv(os.path.join(BASE, "results", "GSE136831_subject_celltype_composition.csv"), index=False)
res.sort_values("padj").to_csv(os.path.join(BASE, "results", "GSE136831_celltype_composition_test.csv"), index=False)
print(res.sort_values("padj").head(20).to_string(index=False))
print(f"tested_celltypes={len(res)} fdr_lt_005={(res.padj < 0.05).sum()}")
