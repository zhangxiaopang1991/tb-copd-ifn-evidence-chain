# -*- coding: utf-8 -*-
"""解析 GSE276060 bulk RNA-seq CSV -> counts 矩阵 + 样本分组表"""
import gzip, re, csv, os

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "GSE276060")

# 1) 从 series_matrix 提取 title -> description -> group 映射
title2desc = {}
with gzip.open(f"{BASE}/GSE276060_series_matrix.txt.gz", "rt") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line.startswith("!Sample_"):
            continue
        parts = line.split("\t")
        key = parts[0]
        vals = [v.strip('"') for v in parts[1:]]
        if key == "!Sample_title":
            titles = vals
        elif key == "!Sample_description":
            descs = vals
for t, d in zip(titles, descs):
    title2desc[t] = d

def desc_to_group(d):
    # Case MTB4_18-070 -> MTB ; Case NTM5_18-073-7 -> NTM ; Case Ctl-Lung1_... -> Ctl_Lung ;
    # Case Ctl-Tonsil1_... -> Ctl_Tonsil ; Case Granulo1_... -> Granuloma
    body = d.split("Case ", 1)[-1]
    if body.startswith("MTB"):
        return "MTB"
    if body.startswith("NTM"):
        return "NTM"
    if body.startswith("Ctl-Lung"):
        return "Ctl_Lung"
    if body.startswith("Ctl-Tonsil"):
        return "Ctl_Tonsil"
    if body.startswith("Granulo"):
        return "Granuloma"
    return "Other"

# 2) 解析 CSV header，定位 Total counts 列
with gzip.open(f"{BASE}/GSE276060_ExpressionRNAseq.csv.gz", "rt") as f:
    header = f.readline().rstrip("\n")
cols = header.split(",")
# 去掉首尾引号
cols = [c.strip('"') for c in cols]

meta_cols = cols[:5]  # Sample/Name/Chromosome/Region/Identifier
data_cols = cols[5:]

# 解析每个 data 列 -> (sample, metric)
pat = re.compile(r'^(.+?)_S\d+ \(GE\) - (Total counts|RPKM|TPM|CPM)$')
sample_metric = []
for c in data_cols:
    m = pat.match(c)
    assert m, f"无法解析列名: {c}"
    sample_metric.append((m.group(1), m.group(2)))

# 提取 Total counts 列索引（相对 data_cols）
count_idx = [i for i, (s, mt) in enumerate(sample_metric) if mt == "Total counts"]
count_samples = [sample_metric[i][0] for i in count_idx]
print(f"Total counts 样本数: {len(count_samples)}")

# 3) 逐行读取，抽取 counts
counts = {}  # gene -> list of ints
gene_names = []
with gzip.open(f"{BASE}/GSE276060_ExpressionRNAseq.csv.gz", "rt") as f:
    f.readline()  # skip header
    reader = csv.reader(f)
    for row in reader:
        if len(row) != len(cols):
            continue
        gname = row[1].strip('"')
        if not gname or gname == "Name":
            continue
        vals = []
        for i in count_idx:
            v = row[5 + i]
            v = v.strip('"').strip()
            vals.append(int(float(v)) if v else 0)
        counts.setdefault(gname, []).append(vals)

# 同名基因合并（求和）
merged = {}
for g, rows in counts.items():
    if len(rows) == 1:
        merged[g] = rows[0]
    else:
        # 按位置求和
        n = len(rows[0])
        tot = [0] * n
        for r in rows:
            for k in range(n):
                tot[k] += r[k]
        merged[g] = tot

# 4) 构建样本分组
meta_rows = []
for s in count_samples:
    d = title2desc.get(s, s)
    meta_rows.append([s, desc_to_group(d), d])

# 5) 写出 counts 矩阵 + meta
with open(f"{BASE}/GSE276060_counts.tsv", "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t")
    w.writerow(["gene"] + count_samples)
    for g in sorted(merged.keys()):
        w.writerow([g] + merged[g])

with open(f"{BASE}/GSE276060_meta.tsv", "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t")
    w.writerow(["sample", "group", "description"])
    for r in meta_rows:
        w.writerow(r)

from collections import Counter
print("分组计数:", Counter(m[1] for m in meta_rows))
print("基因数:", len(merged))
print("输出: GSE276060_counts.tsv / GSE276060_meta.tsv")
