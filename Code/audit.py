from pathlib import Path
import re,sys
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import false_discovery_control as bh
out=Path(__file__).resolve().parent/'Evidence';out.mkdir(exist_ok=True)
data=Path(sys.argv[1])
s6=pd.read_csv(sys.argv[2]).set_index('gene')
m=pd.read_csv(data/'sample_meta_final.tsv',sep='\t').set_index('gsm')
matrix=pd.read_csv(data/'GSE47460_gene_matrix.tsv',sep='\t',index_col=0)
m=m[m['disease state'].isin(['Control','Chronic Obstructive Lung Disease'])].copy()
panel=re.findall(r'"([A-Z0-9]+)"',(data/'analyze_gse47460.R').read_text(encoding='utf8').split('core <- c(')[1].split(')')[0])
a=[]
for g in panel:
 for plat,z in m.groupby('gpl'):
    v=matrix.reindex([g],columns=z.index).iloc[0]
    a.append([g,plat,len(z),v.notna().sum(),g in s6.index])
pd.DataFrame(a,columns=['gene','platform','eligible_cases_controls','observed','in_S6']).to_csv(out/'panel_origin_platform.csv',index=False)
m=m[m['disease state']!='Control'].copy()
m['fev']=pd.to_numeric(m['%predicted fev1 (pre-bd)'],errors='coerce')
m['age']=pd.to_numeric(m.age,errors='coerce')
m['smoking']=m['smoker?'].map(lambda s:'Never' if 'Never' in str(s) else 'Current' if 'Current' in str(s) else 'Ever' if pd.notna(s) else None)
m=m.dropna(subset=['fev','age','smoking','Sex','gpl'])
assert len(m)==215
m.join(matrix.loc[s6.index,m.index].T).to_csv(out/'lung_model_public_inputs.csv')
rows=[];miss=[];influence=[]
for g in s6.index:
 m['y']=matrix.loc[g,m.index]
 for k,z in [('all',m),*m.groupby('gpl')]:
    f='y~fev+age+C(Sex)+C(smoking)'+('+C(gpl)' if k=='all' else '')
    model=smf.ols(f,z).fit()
    miss.append([g,k,len(z),z.y.isna().sum()])
    for mode,fit in [('OLS',model),('HC3',smf.ols(f,z).fit(cov_type='HC3',use_t=True))]:
        rows.append([g,k,mode,int(fit.nobs),fit.params.fev,*fit.conf_int().loc['fev'],fit.pvalues.fev])
    if k=='all':
        h=model.get_influence()
        j=model.model.exog_names.index('fev')
        x=model.model.exog
        delta=(x@np.linalg.inv(x.T@x))[:,j]*model.resid.to_numpy()/(1-h.hat_matrix_diag)
        loo=model.params.fev-delta
        influence.append([g,h.cooks_distance[0].max(),loo.min(),loo.max()])
r=pd.DataFrame(rows,columns=['gene','platform','model','n','beta','lower95','upper95','p'])
r['q20']=r.groupby(['platform','model']).p.transform(bh)
r.to_csv(out/'platform_OLS_HC3.csv',index=False)
primary=r.query("platform=='all' and model=='OLS'").set_index('gene')
assert np.allclose(primary.loc[s6.index,'beta'],s6.beta,rtol=0,atol=1e-10)
pd.DataFrame(miss,columns=['gene','platform','eligible','missing']).to_csv(out/'case_missing.csv',index=False)
pd.DataFrame(influence,columns=['gene','max_Cooks_D','min_leave_one_out_beta','max_leave_one_out_beta']).to_csv(out/'influence.csv',index=False)
stats=[]
for k,z in [('all',m),*m.groupby('gpl')]:
 for field in ['age','fev','Sex','smoking']:
    if field in ['age','fev']:
        for v in [.25,.5,.75]:stats.append([k,field,str(v),z[field].quantile(v),len(z)])
    else:
        for value,n in z[field].value_counts().items():stats.append([k,field,value,n,len(z)])
pd.DataFrame(stats,columns=['platform','variable','level_or_quantile','value','n']).to_csv(out/'cohort.csv',index=False)
print(r.query("platform=='all' and q20<.05")[['gene','model','q20']].to_string(index=False))
print('missing',sum(x[-1] for x in miss),'omitted', [i for i in a if i[0] not in s6.index])

