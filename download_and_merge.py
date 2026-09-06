import csv, os, urllib.request, concurrent.futures, sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.makedirs("counts", exist_ok=True)

# 1) 解析 SDRF -> 42 个 Source Name + 元数据
rows = list(csv.reader(open("EMTAB17246.sdrf.txt", encoding="utf-8"), delimiter="\t"))
header = rows[0]
idx = {n: i for i, n in enumerate(header)}
def col(r, n):
    return r[idx[n]] if n in idx and idx[n] < len(r) else ""

samples = []
seen = set()
for r in rows[1:]:
    if not r:
        continue
    sn = col(r, "Source Name")
    if not sn or sn in seen:
        continue
    seen.add(sn)
    samples.append({
        "sn": sn,
        "donor": col(r, "Characteristics[individual]"),
        "sex": col(r, "Characteristics[sex]"),
        "age": col(r, "Characteristics[age]"),
        "infect": col(r, "Characteristics[infect]"),
        "stimulus": col(r, "Characteristics[stimulus]"),
    })
print(f"解析到 {len(samples)} 个样本", file=sys.stderr)

# 2) 下载
base = "https://www.ebi.ac.uk/biostudies/files/E-MTAB-17246/"
def dl(s):
    sn = s["sn"]
    p = os.path.join("counts", sn + ".count")
    if os.path.exists(p) and os.path.getsize(p) > 1000:
        return sn, "cached"
    last = None
    for _ in range(3):
        try:
            urllib.request.urlretrieve(base + sn + ".count", p)
            return sn, "ok"
        except Exception as e:
            last = e
    return sn, f"err:{last}"

with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
    for sn, st in ex.map(dl, samples):
        if st not in ("ok", "cached"):
            print(f"  !! {sn} -> {st}", file=sys.stderr)

# 3) 合并
genes = {}   # ENSG -> symbol
mat = {}     # ENSG -> {sn: count}
for s in samples:
    sn = s["sn"]
    p = os.path.join("counts", sn + ".count")
    if not os.path.exists(p):
        continue
    with open(p) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            g, sym, cnt = parts[0], parts[1], parts[2]
            if g.startswith("ENSG"):
                genes[g] = sym
                mat.setdefault(g, {})[sn] = cnt

sns = [s["sn"] for s in samples]
with open("counts_matrix.tsv", "w") as f:
    f.write("gene_id\t" + "\t".join(sns) + "\n")
    for g in sorted(genes):
        row = [mat.get(g, {}).get(sn, "0") for sn in sns]
        f.write(g + "\t" + "\t".join(row) + "\n")

with open("sample_table.tsv", "w") as f:
    f.write("sample\tdonor\tsex\tage\tinfect\tstimulus\n")
    for s in samples:
        f.write(f"{s['sn']}\t{s['donor']}\t{s['sex']}\t{s['age']}\t{s['infect']}\t{s['stimulus']}\n")

print(f"DONE 基因数={len(genes)} 样本数={len(sns)}")
