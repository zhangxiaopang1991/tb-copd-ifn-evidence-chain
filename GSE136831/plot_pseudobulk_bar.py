# -*- coding: utf-8 -*-
"""
GSE136831 伪 bulk 条形图（Welch's t + BH-FDR + 95% CI）：
- 三个 panel：PTLD_IFN / IFNa / IFNg
- 每个细胞类型画一根 mean_diff 条，带 95% CI 误差棒
- 在 FDR<0.05 的条上画星号
- 标题、副标题显式写明 "pseudo-bulk, subject-level, BH-FDR, 95% CI"
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, 'results')
os.makedirs(OUT, exist_ok=True)

# 中文字体（沿用其他图配置）
font_manager.fontManager.addfont('C:/Windows/Fonts/msyh.ttc')
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

df = pd.read_csv(os.path.join(OUT, 'GSE136831_pseudobulk_celltype_test.csv'))

# 按 mean_diff 排序
fig, axes = plt.subplots(1, 3, figsize=(15, 8), sharey=True)
score_labels = {'PTLD_IFN': 'PTLD-IFN', 'IFNa': 'IFN-α', 'IFNg': 'IFN-γ'}
score_colors = {'PTLD_IFN': '#4878CF', 'IFNa': '#E1812C', 'IFNg': '#6ACC64'}

for ax, score in zip(axes, ['PTLD_IFN', 'IFNa', 'IFNg']):
    sub = df[df['score'] == score].sort_values('mean_diff').reset_index(drop=True)
    celltypes = sub['celltype'].tolist()
    mean_diff = sub['mean_diff'].values
    padj = sub['padj'].values
    pval = sub['pval'].values

    colors = []
    for md in mean_diff:
        colors.append('#C0392B' if md > 0 else '#2E86C1')  # 红=上调，蓝=下调（中文股票配色习惯）
    # 误差棒：mean_diff 与 95% CI 的距离
    err_lo = np.clip(mean_diff - sub['ci95_low'].values, 0, None)
    err_hi = np.clip(sub['ci95_high'].values - mean_diff, 0, None)
    xerr = np.vstack([err_lo, err_hi])
    bars = ax.barh(celltypes, mean_diff, xerr=xerr, color=colors, alpha=0.85,
                   edgecolor='black', linewidth=0.3,
                   error_kw={'ecolor':'black','elinewidth':0.7,'capsize':3,'capthick':0.7})

    # FDR 星号
    for i, (md, pj, pv) in enumerate(zip(mean_diff, padj, pval)):
        pj = float(pj); pv = float(pv); md = float(md)
        if pj < 0.05:
            ax.text(md + (0.005 if md > 0 else -0.005), i,
                    '**' if pj < 0.01 else '*',
                    va='center', ha='left' if md > 0 else 'right',
                    fontsize=10, fontweight='bold')
        elif pv < 0.05:
            ax.text(md + (0.005 if md > 0 else -0.005), i,
                    '·', va='center', ha='left' if md > 0 else 'right',
                    fontsize=12, color='gray', alpha=0.6)

    ax.axvline(0, color='black', linewidth=0.6)
    ax.set_title(f'{score_labels[score]} score', fontsize=12, fontweight='bold')
    ax.set_xlabel('mean difference\n(COPD − Control)', fontsize=10)
    ax.grid(axis='x', linestyle=':', alpha=0.4)

axes[0].set_ylabel('Cell type', fontsize=11)

# 总标题
fig.suptitle('GSE136831 COPD 肺单细胞 — 伪 bulk 层面 IFN 模块差异（subject-level, Welch t, BH-FDR, 95% CI）',
             fontsize=13, fontweight='bold', y=1.00)

# 图例
from matplotlib.patches import Patch
legend_elems = [
    Patch(facecolor='#C0392B', alpha=0.85, label='COPD 上调'),
    Patch(facecolor='#2E86C1', alpha=0.85, label='COPD 下调'),
    plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='black', markersize=12, label='BH-FDR < 0.05', linestyle=''),
    plt.Line2D([0], [0], marker='.', color='w', markerfacecolor='gray', markersize=14, label='raw p < 0.05 (FDR ns)', linestyle=''),
]
fig.legend(handles=legend_elems, loc='lower center', ncol=4, fontsize=10,
           bbox_to_anchor=(0.5, -0.04))

plt.tight_layout(rect=[0, 0.03, 1, 0.97])
out_png = os.path.join(OUT, 'GSE136831_IFN_pseudobulk_bar.png')
plt.savefig(out_png, dpi=200, bbox_inches='tight', facecolor='white')
print(f'saved: {out_png}')

# 顺带打印 FDR<0.05 总结
sig = df[df['padj'] < 0.05].sort_values('padj')
print(f'\nFDR<0.05 significant: {len(sig)} rows')
if len(sig):
    print(sig[['celltype','score','mean_diff','pval','padj']].to_string(index=False))