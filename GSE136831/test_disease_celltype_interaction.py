# -*- coding: utf-8 -*-
"""Test disease-by-cell-type interaction using subject-level IFN summaries."""
import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

BASE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(BASE, "results")
INPUT = os.path.join(RESULTS, "GSE136831_pseudobulk_subject_celltype.csv")
SCORES = ["PTLD_IFN", "IFNa", "IFNg"]

pb = pd.read_csv(INPUT)
pb = pb[pb["Disease"].isin(["COPD", "Control"])].copy()
pb["Disease"] = pd.Categorical(pb["Disease"], categories=["Control", "COPD"])

# Restrict the omnibus interaction to cell types with at least three subjects
# in each disease group, matching the project's formal replication rule.
counts = (
    pb.groupby(["Manuscript_Identity", "Disease"], observed=False)["Subject_Identity"]
    .nunique()
    .unstack(fill_value=0)
)
valid_types = counts.index[
    (counts.get("COPD", 0) >= 3) & (counts.get("Control", 0) >= 3)
].tolist()
pb = pb[pb["Manuscript_Identity"].isin(valid_types)].copy()
pb["celltype"] = pd.Categorical(pb["Manuscript_Identity"],
                                categories=sorted(valid_types))

summary_rows = []
detail_rows = []
for score in SCORES:
    dat = pb[["Subject_Identity", "Disease", "celltype", score]].dropna().copy()
    formula = f"{score} ~ C(Disease) * C(celltype)"
    fit = smf.ols(formula, data=dat).fit(
        cov_type="cluster", cov_kwds={"groups": dat["Subject_Identity"]}
    )
    interaction_names = [
        name for name in fit.params.index
        if "C(Disease)[T.COPD]:C(celltype)" in name
    ]
    if not interaction_names:
        raise RuntimeError(f"No interaction terms found for {score}")
    r = np.zeros((len(interaction_names), len(fit.params)))
    for i, name in enumerate(interaction_names):
        r[i, list(fit.params.index).index(name)] = 1.0
    wald = fit.wald_test(r, use_f=False)
    pval = float(np.asarray(wald.pvalue).reshape(-1)[0])
    summary_rows.append({
        "score": score,
        "n_subject_celltype": len(dat),
        "n_subjects": dat["Subject_Identity"].nunique(),
        "n_celltypes": dat["celltype"].nunique(),
        "n_interaction_terms": len(interaction_names),
        "interaction_chi2": float(np.asarray(wald.statistic).reshape(-1)[0]),
        "interaction_pval": pval,
        "model_r2": fit.rsquared,
    })

    for name in interaction_names:
        celltype = name.split("C(celltype)[T.", 1)[1].rstrip("]")
        detail_rows.append({
            "score": score,
            "celltype": celltype,
            "interaction_coefficient": fit.params[name],
            "se_cluster_robust": fit.bse[name],
            "pval_wald_term": fit.pvalues[name],
            "ci95_low": fit.conf_int().loc[name, 0],
            "ci95_high": fit.conf_int().loc[name, 1],
        })

summary = pd.DataFrame(summary_rows)
summary["padj_BH_across_scores"] = sm.stats.multipletests(
    summary["interaction_pval"], method="fdr_bh"
)[1]
summary["interaction_fdr_significant"] = (
    summary["padj_BH_across_scores"] < 0.05
)
detail = pd.DataFrame(detail_rows)

summary.to_csv(
    os.path.join(RESULTS, "GSE136831_disease_celltype_interaction_summary.csv"),
    index=False,
)
detail.to_csv(
    os.path.join(RESULTS, "GSE136831_disease_celltype_interaction_terms.csv"),
    index=False,
)
print(summary.to_string(index=False))
print(f"valid_celltypes={len(valid_types)}")
