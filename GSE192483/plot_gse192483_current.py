# -*- coding: utf-8 -*-
"""Plot current GSE192483 directional summaries from the rerun CSV files."""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(BASE, "results")

cell = pd.read_csv(os.path.join(RESULTS, "GSE192483_IFN_score_H_vs_L_celllevel.csv"))
paired = pd.read_csv(os.path.join(RESULTS, "GSE192483_IFN_score_H_vs_L_pseudobulk_paired.csv"))

sns.set_theme(style="whitegrid", font_scale=0.85)

# Show the largest directional cell-level contrasts, without significance labels.
top = cell.assign(abs_diff=cell["mean_diff_HmL"].abs()).sort_values("abs_diff", ascending=False).head(18)
top = top.sort_values("mean_diff_HmL")
fig, ax = plt.subplots(figsize=(8.5, 6.5))
colors = ["#3B82F6" if x < 0 else "#D95F02" for x in top["mean_diff_HmL"]]
ax.barh([f"{c} | {s}" for c, s in zip(top["celltype"], top["score"])], top["mean_diff_HmL"], color=colors)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlabel("Mean score difference (H - L)")
ax.set_title("GSE192483 current rerun: directional cell-level contrasts")
fig.tight_layout()
fig.savefig(os.path.join(RESULTS, "GSE192483_H_vs_L_IFN_bar.png"), dpi=300)
plt.close(fig)

# Heatmap uses donor-aware mean differences; it is descriptive because no row is FDR-significant.
heat = paired.pivot(index="celltype", columns="score", values="mean_diff_HmL").fillna(0)
heat = heat.loc[heat.abs().max(axis=1).sort_values(ascending=False).index]
fig, ax = plt.subplots(figsize=(6.5, max(5.0, 0.32 * len(heat))))
sns.heatmap(heat, center=0, cmap="RdBu_r", annot=True, fmt="+.2f", linewidths=0.3, ax=ax, cbar_kws={"label": "Mean paired difference (H - L)"})
ax.set_title("GSE192483 current rerun: donor-aware directional summary")
ax.set_xlabel("IFN module")
ax.set_ylabel("Cell type")
fig.tight_layout()
fig.savefig(os.path.join(RESULTS, "GSE192483_IFN_heatmap_H_vs_L.png"), dpi=300)
plt.close(fig)

print("wrote current GSE192483 plots")
