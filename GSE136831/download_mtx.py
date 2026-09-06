# -*- coding: utf-8 -*-
"""断点续传下载 GSE136831 稀疏矩阵（2G），单连接 + Range，抗中断。"""
import os, time, urllib.request, sys

URL = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE136831&format=file&file=GSE136831_RawCounts_Sparse.mtx.gz'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'GSE136831_RawCounts_Sparse.mtx.gz')

def get_size(url):
    req = urllib.request.Request(url, method='HEAD')
    r = urllib.request.urlopen(req, timeout=60)
    return int(r.headers.get('Content-Length', 0))

total = get_size(URL)
print(f'remote total: {total/1024/1024:.1f} MB', flush=True)

downloaded = os.path.getsize(OUT) if os.path.exists(OUT) else 0
print(f'local existing: {downloaded/1024/1024:.1f} MB', flush=True)

max_retries = 50
while downloaded < total:
    try:
        req = urllib.request.Request(URL, headers={'Range': f'bytes={downloaded}-'})
        r = urllib.request.urlopen(req, timeout=120)
        with open(OUT, 'ab') as f:
            t0 = time.time(); chunk_bytes = 0
            while True:
                chunk = r.read(1024 * 256)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk); chunk_bytes += len(chunk)
                if time.time() - t0 >= 5:
                    speed = chunk_bytes / (time.time() - t0)
                    print(f'  {downloaded/1024/1024:.1f}/{total/1024/1024:.1f} MB  {speed/1024:.0f} KB/s', flush=True)
                    t0 = time.time(); chunk_bytes = 0
        print(f'completed: {downloaded/1024/1024:.1f} MB', flush=True)
        break
    except Exception as e:
        print(f'  interrupted at {downloaded/1024/1024:.1f} MB: {e}  -> retrying', flush=True)
        time.sleep(3)
        max_retries -= 1
        if max_retries <= 0:
            print('too many retries, giving up', flush=True)
            sys.exit(1)

print('=== DONE ===', flush=True)
