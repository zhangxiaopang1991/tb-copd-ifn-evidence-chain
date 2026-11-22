from pathlib import Path
import csv, gzip, re, sys
import fitz
import numpy as np
import pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

out = Path(__file__).resolve().parent
data = Path(sys.argv[1])
with fitz.open(Path(sys.argv[2])) as pdf:
    values = pdf[5].get_text().splitlines()
cases = {}
for i, value in enumerate(values):
    if re.fullmatch(r'Case\d+', value):
        cases[value] = (int(values[i+2]), int(values[i+4]=='Male'), float(values[i+9]))
with gzip.open(data/'GSE8581_series_matrix.txt.gz','rt') as stream:
    for line in stream:
        fields = next(csv.reader([line], delimiter='\t'))
        if not fields: continue
        if fields[0]=='!Sample_title': titles = fields[1:]
        if fields[0]=='!Sample_geo_accession': accessions = fields[1:]
        if fields[0]=='!series_matrix_table_begin': break
    matrix = pd.read_csv(stream,sep='\t',index_col=0,comment='!')
meta = [(g, name.removeprefix('Human_COPD_')) for g,name in zip(accessions,titles) if name.removeprefix('Human_COPD_') in cases]
assert len(meta)==len(cases)==len(set(k for _,k in meta))==15
clinical = pd.DataFrame([cases[k] for _,k in meta],columns=['age','male','fev'])
with gzip.open(data/'GPL570.annot.gz','rt') as stream:
    for line in stream:
        if line.startswith('ID\t'): break
    annot=pd.read_csv(stream,sep='\t',names=line.strip().split('\t'),usecols=['ID','Gene symbol'],comment='!',dtype=str)
x=np.column_stack([np.ones(15),clinical.fev/10,clinical.age,clinical.male])
assert np.linalg.matrix_rank(x)==4
rows=[]
primary=pd.read_csv(out/'Supplement/Supplementary_Table_S12_lung_check.csv')
for gene in primary.gene:
    probes=annot.loc[annot['Gene symbol']==gene,'ID'].tolist()
    m=matrix.loc[probes,[g for g,_ in meta]].astype(float)
    assert np.isfinite(m).all().all() and m.min().min()>0
    m=np.log2(m)
    probe=m.mean(axis=1).idxmax()
    y=m.loc[probe].to_numpy()
    b=np.linalg.lstsq(x,y,rcond=None)[0]
    se=np.sqrt(np.sum((y-x@b)**2)/11*np.linalg.inv(x.T@x)[1,1])
    half=t.ppf(.975,11)*se
    rows.append([gene,probe,15,b[1],b[1]-half,b[1]+half,2*t.sf(abs(b[1]/se),11)])
res=pd.DataFrame(rows,columns=['gene','probe','n','beta_per10','lower95','upper95','p'])
res['q_BH5']=multipletests(res.p,method='fdr_bh')[1]
res['same_direction_GSE47460']=np.sign(res.beta_per10)==np.sign(primary.GSE47460_beta_per10)
res.to_csv(out/'Supplement/Supplementary_Table_S13_GSE8581.csv',index=False)
print(res.to_string(index=False))
