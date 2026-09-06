# -*- coding: utf-8 -*-
"""
GSE192483 人肺结核肺组织单细胞：donor-aware 配对（H 病灶 vs L 未受累）。
数据：当前 extracted/ 下的全部 .h5（CellRanger filtered_feature_bc_matrix），5 位患者、11 个样本，5 对 H+L 均可配对。
流程：读 h5 -> QC -> normalize/HVG/PCA/UMAP/Leiden -> marker 注释 -> IFN 评分
      -> H vs L（细胞级 MWU + donor-aware 伪bulk配对 Wilcoxon）。
"""
import os, glob
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu, wilcoxon

BASE = os.path.dirname(os.path.abspath(__file__))
EXTRACT = os.path.join(BASE, 'extracted')
OUT = os.path.join(BASE, 'results')
os.makedirs(OUT, exist_ok=True)
sc.settings.verbosity = 1
sc.set_figure_params(dpi=100)
sc.settings.n_jobs = 1
np.random.seed(20260904)

# ---------- 1. 读 h5 ----------
h5files = sorted(glob.glob(os.path.join(EXTRACT, '*.h5')))
print(f'[1] found {len(h5files)} h5 files')
adatas = []
for f in h5files:
    base = os.path.basename(f)
    sample = base.replace('.filtered_feature_bc_matrix.h5', '')   # GSM5747736_SP019H
    label = sample.split('_')[-1]      # SP019H
    patient = label[:-1]               # SP019
    lesion = label[-1]                 # H / L
    adata = sc.read_10x_h5(f, gex_only=True)
    adata.var_names_make_unique()
    adata.obs['sample'] = sample
    adata.obs['patient'] = patient
    adata.obs['lesion'] = lesion
    adatas.append(adata)
    print(f'    {sample}: {adata.n_obs} cells x {adata.n_vars} genes (patient={patient}, lesion={lesion})')

adata = sc.concat(adatas, join='outer', label='batch')
adata.obs_names_make_unique()
print(f'[1] merged: {adata.n_obs} cells x {adata.n_vars} genes')

# ---------- 2. QC ----------
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
before = adata.n_obs
adata = adata[(adata.obs.n_genes_by_counts >= 200) & (adata.obs.pct_counts_mt < 20), :].copy()
print(f'[2] QC: {before} -> {adata.n_obs} cells')

# ---------- 3. normalize + HVG + PCA + UMAP + Leiden（仅对 HVG 缩放，避免稠密化）----------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='sample')
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack')
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.6, flavor='igraph', n_iterations=2, directed=False)
adata.obs['leiden'] = adata_hvg.obs['leiden'].astype(str)
adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
del adata_hvg
print(f'[3] clustering: {adata.obs["leiden"].nunique()} clusters')

# ---------- 4. marker 注释 ----------
markers = {
    'Alveolar_Mac': ['CD68','MARCO','FABP4','MRC1','ITGAX','MERTK'],
    'Monocyte': ['LYZ','CD14','FCGR3A','S100A8','S100A9','VCAN'],
    'cDC': ['FCER1A','CLEC9A','CD1C','BATF3','IRF8'],
    'pDC': ['IL3RA','LILRA4','GZMB','CLEC4C'],
    'T_cell': ['CD3D','CD3E','CD2','IL7R','TRAC'],
    'CD4_T': ['CD4','IL7R'],
    'CD8_T': ['CD8A','CD8B'],
    'NK': ['NKG7','GNLY','KLRD1','NCAM1','KLRB1'],
    'B_cell': ['MS4A1','CD79A','CD79B','CD19','BANK1'],
    'Plasma': ['IGHG1','MZB1','JCHAIN','XBP1'],
    'Mast': ['TPSAB1','TPSB2','CPA3','KIT','MS4A2'],
    'Epithelial': ['EPCAM','KRT18','KRT8'],
    'AT1': ['AGER','HOPX','AQP5','CAV1'],
    'AT2': ['SFTPC','SFTPA1','SFTPB','LAMP3','NAPSA'],
    'Club': ['SCGB1A1','SCGB3A1','SCGB3A2'],
    'Ciliated': ['FOXJ1','TPPP3','PIFO','SNTN'],
    'Basal': ['KRT5','TP63','KRT14','KRT15'],
    'Goblet': ['MUC5AC','MUC5B','TFF3'],
    'Endothelial': ['PECAM1','VWF','CLDN5','CDH5','FLT1'],
    'Fibroblast': ['COL1A1','COL1A2','DCN','PDGFRA','LUM'],
    'Smooth_Muscle': ['ACTA2','MYH11','TAGLN','DES'],
    'Lymphatic': ['PROX1','LYVE1','PDPN'],
}
for name, gs in markers.items():
    gs = [g for g in gs if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=gs, score_name=f'mk_{name}', ctrl_size=50)

mk_scores = [f'mk_{k}' for k in markers if f'mk_{k}' in adata.obs.columns]
cluster_marker = adata.obs.groupby('leiden')[mk_scores].mean().idxmax(axis=1)
cluster_marker = cluster_marker.str.replace('^mk_', '', regex=True)
adata.obs['celltype'] = adata.obs['leiden'].map(cluster_marker)
print('[4] annotation done:')
print(adata.obs.groupby('celltype').size().sort_values(ascending=False))
ann_map = adata.obs.groupby('leiden')['celltype'].first()
ann_map.to_csv(os.path.join(OUT, 'GSE192483_cluster_celltype_map.csv'), header=['celltype'])

# ---------- 5. IFN 模块评分 ----------
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

# ---------- 6. 细胞级 H vs L（每细胞类型）----------
print('[6] cell-level H vs L per celltype ...')
rows = []
for ct in sorted(adata.obs['celltype'].unique()):
    sub = adata[adata.obs['celltype'] == ct]
    h = sub.obs['lesion'] == 'H'
    l = sub.obs['lesion'] == 'L'
    if h.sum() < 10 or l.sum() < 10:
        continue
    for score in ['PTLD_IFN_score','IFNa_score','IFNg_score']:
        vh = sub.obs.loc[h, score].values.astype(float)
        vl = sub.obs.loc[l, score].values.astype(float)
        u, p = mannwhitneyu(vh, vl, alternative='two-sided')
        rows.append({'celltype': ct, 'score': score.replace('_score',''),
                     'n_H': int(h.sum()), 'n_L': int(l.sum()),
                     'mean_H': np.mean(vh), 'mean_L': np.mean(vl),
                     'mean_diff_HmL': np.mean(vh)-np.mean(vl), 'pval': p})
res = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
for score in ['PTLD_IFN','IFNa','IFNg']:
    m = res['score'] == score
    res.loc[m, 'padj'] = multipletests(res.loc[m,'pval'].fillna(1), method='fdr_bh')[1]
res.to_csv(os.path.join(OUT, 'GSE192483_IFN_score_H_vs_L_celllevel.csv'), index=False)
print('[6] done, top by pval:')
print(res.sort_values('pval').head(15)[['celltype','score','mean_diff_HmL','pval','padj','n_H','n_L']].to_string(index=False))

# ---------- 7. donor-aware 伪bulk 配对（H vs L，Wilcoxon signed-rank）----------
print('[7] donor-aware pseudo-bulk paired (H vs L) ...')
# 伪 bulk：patient x celltype 的平均 IFN score
adata.obs['patient_ct'] = adata.obs['patient'].astype(str) + '|' + adata.obs['celltype'].astype(str)
pb = adata.obs.groupby(['patient_ct','patient','celltype','lesion']).agg(
    PTLD_IFN=('PTLD_IFN_score','mean'), IFNa=('IFNa_score','mean'), IFNg=('IFNg_score','mean'),
    ncell=('PTLD_IFN_score','size')).reset_index()

pb_rows = []
for ct in sorted(pb['celltype'].unique()):
    s = pb[pb['celltype'] == ct]
    for score in ['PTLD_IFN','IFNa','IFNg']:
        # 找同时有 H 和 L 的患者（配对）
        piv = s.pivot(index='patient', columns='lesion', values=score).dropna()
        if len(piv) < 3:   # 至少 3 对才做配对检验
            continue
        h = piv['H'].values.astype(float)
        l = piv['L'].values.astype(float)
        try:
            w, p = wilcoxon(h, l)
        except Exception:
            p = np.nan
        pb_rows.append({'celltype': ct, 'score': score, 'n_pairs': len(piv),
                        'mean_H': h.mean(), 'mean_L': l.mean(),
                        'mean_diff_HmL': h.mean()-l.mean(), 'wilcoxon_p': p})
pbres = pd.DataFrame(pb_rows)
if len(pbres):
    for score in ['PTLD_IFN','IFNa','IFNg']:
        m = pbres['score'] == score
        pbres.loc[m, 'padj'] = multipletests(pbres.loc[m,'wilcoxon_p'].fillna(1), method='fdr_bh')[1]
pbres.to_csv(os.path.join(OUT, 'GSE192483_IFN_score_H_vs_L_pseudobulk_paired.csv'), index=False)
print('[7] donor-aware paired done, top by pval:')
if len(pbres):
    print(pbres.sort_values('wilcoxon_p').head(15)[['celltype','score','n_pairs','mean_diff_HmL','wilcoxon_p','padj']].to_string(index=False))

# ---------- 8. 保存 ----------
adata.write(os.path.join(OUT, 'GSE192483_annotated.h5ad'))
print('\nDONE')
