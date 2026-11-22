"""Case-only secondary check; recruiting-site and donor independence remain unresolved."""
from pathlib import Path
import csv, gzip, hashlib
from io import StringIO
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import t as tdist

out=Path(__file__).resolve().parent
import sys
inputs=Path(sys.argv[1])
series=inputs/'GSE76925_series_matrix.txt.gz'
annotation=inputs/'GPL10558.annot.gz'
genes=['IFIT1','IFI44','USP18','BST2','IFITM3']
with gzip.open(annotation,'rt') as f:
    found=False
    for line in f:
        if line.startswith('ID\t'):
            probes=pd.read_csv(StringIO(line+f.read().split('!platform_table_end')[0]),sep='\t',low_memory=False)
            found=True
            break
    assert found
probes=probes.loc[probes['Gene symbol'].isin(genes),['ID','Gene symbol']]
with gzip.open(series,'rt') as f:
    fields={}
    for line in f:
        vals=next(csv.reader([line],delimiter='\t'))
        if not vals:continue
        if vals[0]=='!Sample_geo_accession': ids=vals[1:]
        if vals[0]=='!Sample_characteristics_ch1':
            pairs=[v.split(': ',1) for v in vals[1:]]
            key=pairs[0][0]
            assert all(p[0]==key for p in pairs)
            fields[key]=[p[1] for p in pairs]
        if line.startswith('!series_matrix_table_begin'):break
    x=pd.read_csv(StringIO(f.read().split('!series_matrix_table_end')[0]),sep='\t',index_col=0)
assert x.columns.tolist()==ids and x.index.is_unique and len(ids)==len(set(ids))==151
meta=pd.DataFrame(fields,index=ids)
cases=meta.loc[meta.copd.eq('case')].copy()
assert len(cases)==111 and cases.ID.is_unique
for c in ['age','fev1.pp','packyears']:cases[c]=pd.to_numeric(cases[c],errors='coerce')
cases['male']=cases.Sex.map({'M':1,'F':0})
cases['fev1_per10']=cases['fev1.pp']/10
rows=[]
for gene in genes:
    candidates=probes.loc[(probes['Gene symbol']==gene)&probes.ID.isin(x.index),'ID']
    assert len(candidates)>0
    # Rank using average case expression only, without association statistics.
    means=x.loc[candidates,cases.index].mean(axis=1).sort_index().sort_values(ascending=False,kind='stable')
    for idx,probe in enumerate(means.index):
        for model,cols in [('primary',['fev1_per10','age','male']),('packyears',['fev1_per10','age','male','packyears'])]:
            data=cases[cols].join(x.loc[probe].rename('response')).dropna()
            design=sm.add_constant(data[cols])
            fit=sm.OLS(data.response,design).fit()
            k=1
            a=design.to_numpy(); y=data.response.to_numpy()
            coef=np.linalg.lstsq(a,y,rcond=None)[0]
            se=np.sqrt(np.linalg.inv(a.T@a)[k,k]*np.sum((y-a@coef)**2)/(len(y)-a.shape[1]))
            assert np.allclose([coef[k],se,2*tdist.sf(abs(coef[k]/se),len(y)-a.shape[1])],[fit.params.iloc[k],fit.bse.iloc[k],fit.pvalues.iloc[k]])
            rows.append([gene,probe,idx+1,model,len(data),fit.params.iloc[k],*fit.conf_int().iloc[k],fit.pvalues.iloc[k]])
res=pd.DataFrame(rows,columns=['gene','probe','expression_rank','model','n','beta_per10','lower95','upper95','p'])
for model in ['primary','packyears']:
    mask=(res.model==model)&(res.expression_rank==1)
    assert mask.sum()==5
    res.loc[mask,'q_BH5']=multipletests(res.loc[mask,'p'],method='fdr_bh')[1]
res.to_csv(out/'five_gene_results.csv',index=False)
print('cases',len(cases),'fields',list(fields),'FEV1 summary',cases['fev1.pp'].describe().to_dict())
print('selected',res[res.expression_rank==1].round(5).to_string(index=False))
print('alternative probes',res[(res.expression_rank>1)&res.model.eq('primary')].round(4).to_string(index=False))
print('SHA256',[(p.name,hashlib.sha256(p.read_bytes()).hexdigest()) for p in [series,annotation]])

