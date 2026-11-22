from pathlib import Path
import gzip,hashlib,sys
import numpy as np
import pandas as pd
d=Path(sys.argv[1]);dest=Path(sys.argv[2]);dest.mkdir(exist_ok=True,parents=True)
m=pd.read_csv(d/'sample_meta_final.tsv',sep='\t')
m=m[m['disease state'].isin(['Control','Chronic Obstructive Lung Disease'])]
rows=[]; annotations=[]
for plat,z in m.groupby('gpl'):
    mapping=d/(plat.lower()+'_id2gene.tsv')
    pro=pd.read_csv(mapping,sep='\t',header=None,dtype=str).query("1 == 1") if False else pd.read_csv(mapping,sep='\t',header=None,dtype=str)
    probes=set(pro.loc[pro[1]=='STAT2',0]);assert probes
    official=(d/(plat.lower()+'_full.txt')).read_text(encoding='utf8',errors='replace').splitlines()
    for probe in probes:
        record=next(l for l in official if l.startswith(probe+'\t'));assert 'STAT2' in record.split('\t')
        annotations.append([plat,probe,record[:180],hashlib.sha256(mapping.read_bytes()).hexdigest()])
    with gzip.open(d/f'GSE47460-{plat}_series_matrix.txt.gz','rt',encoding='utf8',errors='replace') as f:
        hits=[line.split('\t',1)[0].strip('"') for line in f if line.split('\t',1)[0].strip('"') in probes]
    for gsm in z.gsm:
        block=(d/'soft'/(gsm+'.txt')).read_text(encoding='utf8',errors='replace').split('!sample_table_begin')[1].split('!sample_table_end')[0]
        vals=[l.split('\t')[1] for l in block.splitlines() if l.split('\t')[0].strip('"') in probes]
        rows.append([plat,gsm,len(vals),np.isfinite(pd.to_numeric(vals,errors='coerce')).sum(),len(hits)])
r=pd.DataFrame(rows,columns=['platform','gsm','matching_probe_rows','finite_values','series_matrix_matching_rows'])
r.to_csv(dest/'STAT2_deposit_check.csv',index=False)
pd.DataFrame(annotations,columns=['platform','probe','annotation_excerpt','mapping_sha256']).to_csv(dest/'STAT2_annotation.csv',index=False)
print(r.groupby('platform').agg({'gsm':'count','matching_probe_rows':'sum','finite_values':'sum','series_matrix_matching_rows':'max'}).to_string())

