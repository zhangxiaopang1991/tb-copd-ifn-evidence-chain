# -*- coding: utf-8 -*-
"""
GSE136831 COPD 肺单细胞：IFN 模块评分 × 细胞类型 × 疾病状态。
目标：回答「COPD 稳定期 IFN 下调发生在哪些细胞群」，并与 GSE47460 bulk 结论互证。
"""
import gzip, os, sys
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import scanpy as sc

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, 'results')
os.makedirs(OUT, exist_ok=True)

# ---------- 1. 读数据 ----------
print('[1] reading data ...')
genes = pd.read_csv(os.path.join(BASE, 'GSE136831_AllCells.GeneIDs.txt.gz'),
                    sep='\t', header=0)
symbols = genes.iloc[:, 1].astype(str).values   # HGNC symbol
barcodes = pd.read_csv(os.path.join(BASE, 'GSE136831_AllCells.cellBarcodes.txt.gz'),
                       header=None)[0].astype(str).values
meta = pd.read_csv(os.path.join(BASE, 'GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz'),
                   sep='\t')

print(f'  genes: {len(symbols)}, barcodes: {len(barcodes)}, meta rows: {len(meta)}')

# 读 mtx（基因 × 细胞）
with gzip.open(os.path.join(BASE, 'GSE136831_RawCounts_Sparse.mtx.gz'), 'rb') as f:
    mtx = scipy.io.mmread(f)
print(f'  raw mtx shape (genes x cells?): {mtx.shape}')

# 转置为 细胞 × 基因
if mtx.shape[0] == len(symbols) and mtx.shape[1] == len(barcodes):
    X = mtx.T.tocsr()
elif mtx.shape[1] == len(symbols) and mtx.shape[0] == len(barcodes):
    X = mtx.tocsr()
else:
    raise ValueError(f'mtx shape {mtx.shape} 与 genes/barcodes 不匹配')

adata = sc.AnnData(X=X)
adata.obs_names = barcodes
adata.var_names = pd.Index(symbols)
adata.var_names_make_unique()

# 对齐 metadata
meta = meta.set_index('CellBarcode_Identity')
adata.obs = meta.loc[adata.obs_names]
print(f'[1] anndata: {adata.n_obs} cells x {adata.n_vars} genes')

# ---------- 2. 筛选 COPD + Control ----------
print('[2] subset COPD + Control ...')
adata = adata[adata.obs['Disease_Identity'].isin(['COPD', 'Control'])].copy()
adata.obs['Disease'] = adata.obs['Disease_Identity'].astype(str)
print(f'  cells: {adata.n_obs} (COPD {sum(adata.obs.Disease=="COPD")}, Control {sum(adata.obs.Disease=="Control")})')

# ---------- 3. normalize + log ----------
print('[3] normalize + log1p ...')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ---------- 4. IFN 模块评分 ----------
print('[4] scoring IFN gene sets ...')
gs_ptld = ['IFI6','OAS1','OAS2','ISG15','MX1','MX2','IFIT1','IFIT2','IFIT3','STAT1',
           'STAT2','IRF7','IRF9','GBP1','GBP2','GBP5','OASL','IFI44','IFI27','IFI35',
           'IFITM1','IFITM3','RSAD2','USP18','XAF1','BST2','EPSTI1','DDX58','IFIH1','EIF2AK2']
gs_ifna = ['ISG15','MX1','OAS1','OAS2','IFIT1','IFIT2','IFIT3','IFI6','IFI27','IFI44',
           'IFI44L','IFITM1','IFITM2','IFITM3','RSAD2','USP18','XAF1','BST2','EPSTI1',
           'EIF2AK2','GBP1','GBP2','DDX58','IFIH1','IRF7','IRF9','STAT1','STAT2',
           'SAMD9','SAMD9L','PARP9','DTX3L','TRIM22','TRIM25','TNFSF10','ISG20','OASL']
gs_ifng = ['STAT1','STAT2','IRF1','IRF7','IRF9','GBP1','GBP2','GBP5','CXCL9','CXCL10',
           'CXCL11','HLA-A','HLA-B','HLA-C','TAP1','TAP2','B2M','PSMB8','PSMB9','PSME1',
           'PSME2','NLRC5','CIITA','IFITM1','IFITM3','OAS1','OAS2','OASL','MX1','MX2',
           'ISG15','ISG20','IFI35','SOCS1','SOCS3','JAK1','JAK2','CCL2','CCL5','CCL7']

for name, gs in [('PTLD_IFN', gs_ptld), ('IFNa', gs_ifna), ('IFNg', gs_ifng)]:
    gs = [g for g in gs if g in adata.var_names]
    print(f'  {name}: {len(gs)} genes matched')
    sc.tl.score_genes(adata, gene_list=gs, score_name=f'{name}_score', ctrl_size=100)

# ---------- 5. 细胞类型 × 疾病 比较 ----------
print('[5] comparing IFN scores by cell type x disease ...')
celltypes = sorted(adata.obs['Manuscript_Identity'].dropna().unique())
rows = []
for ct in celltypes:
    sub = adata[adata.obs['Manuscript_Identity'] == ct]
    copd = sub.obs['Disease'] == 'COPD'
    ctl = sub.obs['Disease'] == 'Control'
    if copd.sum() < 10 or ctl.sum() < 10:
        continue
    for score in ['PTLD_IFN_score', 'IFNa_score', 'IFNg_score']:
        v_copd = sub.obs.loc[copd, score].values.astype(float)
        v_ctl = sub.obs.loc[ctl, score].values.astype(float)
        mean_diff = np.mean(v_copd) - np.mean(v_ctl)
        # Mann-Whitney U (wilcoxon rank-sum)
        from scipy.stats import mannwhitneyu
        try:
            u, p = mannwhitneyu(v_copd, v_ctl, alternative='two-sided')
        except Exception:
            p = np.nan
        rows.append({'celltype': ct, 'score': score.replace('_score',''),
                     'n_COPD': int(copd.sum()), 'n_Control': int(ctl.sum()),
                     'mean_COPD': np.mean(v_copd), 'mean_Control': np.mean(v_ctl),
                     'mean_diff': mean_diff, 'pval': p})

res = pd.DataFrame(rows)
# 全局 BH-FDR 多重校正（Codex 评审修正：所有 celltype × score 一起校正）
from statsmodels.stats.multitest import multipletests
res['padj'] = multipletests(res['pval'].fillna(1), method='fdr_bh')[1]
res.to_csv(os.path.join(OUT, 'GSE136831_IFN_score_by_celltype.csv'), index=False)
print(f'[5] saved results ({len(res)} rows, global BH-FDR)')

# ---------- 6. 伪 bulk：Subject × CellType 聚合，验证方向 ----------
print('[6] pseudo-bulk by subject x celltype ...')
adata.obs['subject_ct'] = adata.obs['Subject_Identity'].astype(str) + '|' + adata.obs['Manuscript_Identity'].astype(str)
# 每个 subject_ct 的平均 IFN score + 疾病
pb = adata.obs.groupby(['subject_ct', 'Subject_Identity', 'Manuscript_Identity', 'Disease']).agg(
    PTLD_IFN=('PTLD_IFN_score','mean'), IFNa=('IFNa_score','mean'), IFNg=('IFNg_score','mean'),
    ncell=('PTLD_IFN_score','size')).reset_index()
pb.to_csv(os.path.join(OUT, 'GSE136831_pseudobulk_subject_celltype.csv'), index=False)

# 伪 bulk 层面：COPD vs Control（每个细胞类型独立 Welch's t 检验 + 95% CI）
from scipy.stats import ttest_ind, t as t_dist
pb_rows = []
for ct in sorted(pb['Manuscript_Identity'].unique()):
    s = pb[pb['Manuscript_Identity'] == ct]
    for score in ['PTLD_IFN','IFNa','IFNg']:
        a = s.loc[s['Disease']=='COPD', score].dropna()
        b = s.loc[s['Disease']=='Control', score].dropna()
        if len(a) < 3 or len(b) < 3:
            continue
        mean_diff = a.mean() - b.mean()
        # Welch's t 检验
        t_stat, p = ttest_ind(a, b, equal_var=False)
        # Welch-Satterthwaite df + 95% CI for mean_diff
        va, vb = a.var(ddof=1), b.var(ddof=1)
        na, nb = len(a), len(b)
        se = np.sqrt(va/na + vb/nb)
        df_welch = (va/na + vb/nb)**2 / ((va/na)**2/(na-1) + (vb/nb)**2/(nb-1))
        t_crit = t_dist.ppf(0.975, df_welch)
        ci_lo = mean_diff - t_crit * se
        ci_hi = mean_diff + t_crit * se
        # Hedges' g（合并 SD 标准化）作为效应量
        sp = np.sqrt(((na-1)*va + (nb-1)*vb) / (na+nb-2))
        if sp > 0:
            d_cohen = mean_diff / sp
            j = 1 - 3/(4*(na+nb)-9)  # 小样本 Hedges 校正
            g = d_cohen * j
        else:
            d_cohen = g = np.nan
        pb_rows.append({'celltype': ct, 'score': score, 'n_subj_COPD': na,
                        'n_subj_Control': nb, 'mean_diff': mean_diff,
                        'ci95_low': ci_lo, 'ci95_high': ci_hi,
                        'cohens_d': d_cohen, 'hedges_g': g,
                        'pval': p})
pbres = pd.DataFrame(pb_rows)
# 全局 BH-FDR（Codex 评审修正：111 次比较一并校正，而非按 score 分组）
from statsmodels.stats.multitest import multipletests
pbres['padj'] = multipletests(pbres['pval'].fillna(1), method='fdr_bh')[1]
pbres['sig'] = np.where(pbres['padj'] < 0.05,
                        np.where(pbres['mean_diff'] > 0, 'up', 'down'),
                        'ns')
pbres.to_csv(os.path.join(OUT, 'GSE136831_pseudobulk_celltype_test.csv'), index=False)
print(f'[6] pseudo-bulk test done ({len(pbres)} rows, global BH-FDR over all {len(pbres)} comparisons, Welch+95%CI)')
# 显式打印 FDR 后显著条目
sig_rows = pbres[pbres['sig'] != 'ns'].sort_values('padj')
if len(sig_rows):
    print(f'[6] FDR<0.05 significant rows: {len(sig_rows)}')
    print(sig_rows[['celltype','score','mean_diff','ci95_low','ci95_high','hedges_g','pval','padj','sig']].to_string(index=False))
else:
    print('[6] no rows reached FDR<0.05 (raw p<0.05 still reported in CSV)')

# ---------- 7. 汇总输出（top 差异细胞类型）----------
print('[7] summary: top cell types with IFN change (cell-level) ...')
top = res.sort_values('pval').head(20)
print(top[['celltype','score','mean_diff','pval','padj','n_COPD','n_Control']].to_string(index=False))
print('\nDONE')
