# -*- coding: utf-8 -*-
"""Aggregate GSE136831 raw counts to subject x cell type for edgeR."""
import gzip
import os
import re
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "pseudobulk_counts")
os.makedirs(OUT, exist_ok=True)

genes = pd.read_csv(os.path.join(BASE, "GSE136831_AllCells.GeneIDs.txt.gz"), sep="\t")
gene_symbols = genes.iloc[:, 1].astype(str).to_numpy()
barcodes = pd.read_csv(os.path.join(BASE, "GSE136831_AllCells.cellBarcodes.txt.gz"), header=None)[0].astype(str)
meta = pd.read_csv(os.path.join(BASE, "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"), sep="\t")
meta = meta.set_index("CellBarcode_Identity").loc[barcodes]

keep = meta["Disease_Identity"].isin(["COPD", "Control"])
meta = meta.loc[keep].copy()

with gzip.open(os.path.join(BASE, "GSE136831_RawCounts_Sparse.mtx.gz"), "rb") as f:
    mtx = scipy.io.mmread(f).tocsr()
if mtx.shape == (len(gene_symbols), len(barcodes)):
    X = mtx[:, keep.to_numpy()].T.tocsr()
elif mtx.shape == (len(barcodes), len(gene_symbols)):
    X = mtx[keep.to_numpy(), :].tocsr()
else:
    raise ValueError(f"Unexpected matrix shape: {mtx.shape}")

genes_out = pd.DataFrame({"gene_id": genes.iloc[:, 0].astype(str), "symbol": gene_symbols})
genes_out.to_csv(os.path.join(OUT, "genes.tsv"), sep="\t", index=False)

meta = meta.reset_index().rename(columns={"index": "CellBarcode_Identity"})
meta.to_csv(os.path.join(OUT, "cell_metadata_copd_control.tsv"), sep="\t", index=False)

group_keys = meta["Subject_Identity"].astype(str) + "|" + meta["Manuscript_Identity"].astype(str)
meta["group_key"] = group_keys.to_numpy()
groups = meta[["group_key", "Subject_Identity", "Manuscript_Identity", "Disease_Identity"]].drop_duplicates().sort_values(["Manuscript_Identity", "Subject_Identity"])

summary = []
for celltype in sorted(meta["Manuscript_Identity"].dropna().unique()):
    rows = np.flatnonzero(meta["Manuscript_Identity"].to_numpy() == celltype)
    sub = meta.iloc[rows].copy()
    sub_groups = sub[["group_key", "Subject_Identity", "Manuscript_Identity", "Disease_Identity"]].drop_duplicates().sort_values("Subject_Identity")
    group_index = {key: i for i, key in enumerate(sub_groups["group_key"])}
    codes = sub["group_key"].map(group_index).to_numpy()
    indicator = sp.csr_matrix((np.ones(len(rows), dtype=np.int8), (codes, np.arange(len(rows)))),
                              shape=(len(sub_groups), len(rows)))
    agg = (indicator @ X[rows]).tocsr()
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", celltype)
    path = os.path.join(OUT, f"counts_{safe}.mtx.gz")
    with gzip.open(path, "wb") as f:
        scipy.io.mmwrite(f, agg.T)
    sub_groups = sub_groups.reset_index(drop=True)
    sub_groups.insert(0, "sample_id", [f"{safe}__{i+1:03d}" for i in range(len(sub_groups))])
    sub_groups["n_cells"] = np.asarray(indicator.sum(axis=1)).ravel().astype(int)
    sub_groups.to_csv(os.path.join(OUT, f"metadata_{safe}.tsv"), sep="\t", index=False)
    summary.append({
        "celltype": celltype,
        "safe_name": safe,
        "n_pseudobulk": len(sub_groups),
        "n_COPD": int((sub_groups["Disease_Identity"] == "COPD").sum()),
        "n_Control": int((sub_groups["Disease_Identity"] == "Control").sum()),
        "total_cells": int(len(rows)),
    })
    print(f"{celltype}: {len(sub_groups)} pseudobulk samples ({(sub_groups['Disease_Identity']=='COPD').sum()} COPD, {(sub_groups['Disease_Identity']=='Control').sum()} Control)")

pd.DataFrame(summary).to_csv(os.path.join(OUT, "pseudobulk_summary.tsv"), sep="\t", index=False)
print(f"Wrote pseudobulk counts to {OUT}")
