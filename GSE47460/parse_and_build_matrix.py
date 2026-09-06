# -*- coding: utf-8 -*-
"""解析 GSE47460 全部样本 SOFT，构建基因 x 样本表达矩阵（log2 信号强度）。"""
import csv, os, glob
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
SOFT_DIR = os.path.join(BASE, 'soft')

# 1. 读样本表（gsm -> gpl, disease state）
sample_meta = {}
with open(os.path.join(BASE, 'sample_table.tsv'), encoding='utf-8') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        sample_meta[row['gsm']] = row

# 2. 读探针->基因映射（两平台）
def load_map(fn):
    m = {}
    with open(os.path.join(BASE, fn), encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line: continue
            parts = line.split('\t')
            if len(parts) >= 2 and parts[1]:
                m[parts[0]] = parts[1]
    return m
gpl_map = {
    'GPL14550': load_map('gpl14550_id2gene.tsv'),
    'GPL6480': load_map('gpl6480_id2gene.tsv'),
}

# 3. 解析每个样本 SOFT 的表达表
# 结构：gene -> {gsm: value}
gene_vals = defaultdict(dict)   # gene -> gsm -> value
sample_order = []               # gsm 顺序
n_probe_total = 0
n_probe_mapped = 0
missing_files = []
probe_dup_warn = 0

for gsm in sorted(sample_meta.keys()):
    gpl = sample_meta[gsm]['gpl']
    path = os.path.join(SOFT_DIR, f'{gsm}.txt')
    if not os.path.exists(path):
        missing_files.append(gsm)
        continue
    pmap = gpl_map.get(gpl, {})
    # 读表达表
    gene_local = {}   # gene -> value (本样本)
    in_table = False
    header_seen = False
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!sample_table_begin'):
                in_table = True; continue
            if line.startswith('!sample_table_end'):
                in_table = False; break
            if not in_table: continue
            if not header_seen:
                header_seen = True; continue  # 跳过表头 ID_REF VALUE
            parts = line.split('\t')
            if len(parts) < 2: continue
            pid = parts[0].strip('"')
            try:
                val = float(parts[1].strip('"'))
            except ValueError:
                continue
            n_probe_total += 1
            gene = pmap.get(pid, '')
            if not gene:
                continue
            n_probe_mapped += 1
            # 同一样本内多探针映射同一基因：取最大值（与提交者 highest signal 一致）
            if gene in gene_local:
                if val > gene_local[gene]:
                    gene_local[gene] = val
            else:
                gene_local[gene] = val
    for gene, val in gene_local.items():
        gene_vals[gene][gsm] = val
    sample_order.append(gsm)

print(f'samples parsed: {len(sample_order)}')
print(f'missing files: {len(missing_files)}')
if missing_files:
    print('missing:', missing_files[:20])
print(f'total probes (rows): {n_probe_total}')
print(f'probes mapped to gene: {n_probe_mapped}')
print(f'unique genes in matrix: {len(gene_vals)}')

# 4. 输出基因 x 样本矩阵
out_fn = os.path.join(BASE, 'GSE47460_gene_matrix.tsv')
with open(out_fn, 'w', newline='') as out:
    w = csv.writer(out, delimiter='\t')
    header = ['gene'] + sample_order
    w.writerow(header)
    for gene in sorted(gene_vals.keys()):
        row = [gene] + [gene_vals[gene].get(gsm, '') for gsm in sample_order]
        w.writerow(row)
print(f'saved {out_fn}')

# 5. 输出样本元数据（含分组，供 R 使用）
meta_fn = os.path.join(BASE, 'sample_meta_final.tsv')
with open(meta_fn, 'w', newline='') as out:
    fields = ['gsm', 'title', 'gpl', 'disease state', 'gold stage', 'smoker?',
              '%predicted fev1 (pre-bd)', '%predicted fvc (pre-bd)', '%predicted dlco', 'age', 'Sex']
    w = csv.writer(out, delimiter='\t')
    w.writerow(fields)
    for gsm in sample_order:
        m = sample_meta[gsm]
        w.writerow([m.get(f, '') for f in fields])
print(f'saved {meta_fn}')
