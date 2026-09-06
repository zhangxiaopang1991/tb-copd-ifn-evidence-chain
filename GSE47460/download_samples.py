# -*- coding: utf-8 -*-
"""批量下载 GSE47460 全部 582 个样本的完整 SOFT（含表达表），并发、断点续传。"""
import csv, os, sys, time
import concurrent.futures
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
SOFT_DIR = os.path.join(BASE, 'soft')
os.makedirs(SOFT_DIR, exist_ok=True)

# 读样本列表
gsms = []
with open(os.path.join(BASE, 'sample_table.tsv'), encoding='utf-8') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        gsms.append(row['gsm'])
print(f'total samples to download: {len(gsms)}', flush=True)

def dl(gsm):
    path = os.path.join(SOFT_DIR, f'{gsm}.txt')
    # 已存在且含表达表结束标记则跳过
    if os.path.exists(path):
        with open(path, 'rb') as f:
            tail = f.read()[-200:]
        if b'sample_table_end' in tail:
            return gsm, 'skip', 0
    url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full'
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        data = urllib.request.urlopen(req, timeout=120).read()
        with open(path, 'wb') as f:
            f.write(data)
        if b'sample_table_end' not in data[-200:]:
            return gsm, 'incomplete', len(data)
        return gsm, 'ok', len(data)
    except Exception as e:
        return gsm, f'err:{e}', 0

t0 = time.time()
done = 0
ok = 0
skip = 0
fail = []
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
    futs = {ex.submit(dl, g): g for g in gsms}
    for fut in concurrent.futures.as_completed(futs):
        gsm, status, size = fut.result()
        done += 1
        if status == 'ok': ok += 1
        elif status == 'skip': skip += 1
        else:
            fail.append((gsm, status))
        if done % 50 == 0 or done == len(gsms):
            el = time.time() - t0
            print(f'[{done}/{len(gsms)}] ok={ok} skip={skip} fail={len(fail)} elapsed={el:.0f}s', flush=True)

print(f'\nDONE ok={ok} skip={skip} fail={len(fail)}', flush=True)
for g, s in fail:
    print('FAIL', g, s, flush=True)
