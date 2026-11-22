from pathlib import Path
import re
import shutil
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(__file__).resolve().parent
src=out.parent/'COPD_retarget_draft_2026-11-21'
fig=out/'Figures'
source=out/'Figure_source_data'
fig.mkdir(exist_ok=True)
source.mkdir(exist_ok=True)
s12=pd.read_csv(next((out/'lung15/Supplement').glob('*S12*'))).set_index('gene')
s13=pd.read_csv(next((out/'lung15/Supplement').glob('*S13*'))).set_index('gene')
records=[]
for cohort,table,prefix in [('GSE47460',s12,'GSE47460_'),('GSE76925',s12,'GSE76925_primary_'),('GSE8581',s13,'')]:
    for gene in s12.index:
        row=table.loc[gene]
        records.append([gene,cohort,*[row[prefix+k] for k in ['n','beta_per10','lower95','upper95']]])
data=pd.DataFrame(records,columns=['gene','cohort','n','beta_per10','lower95','upper95'])
assert len(data)==15 and (data.lower95<data.beta_per10).all() and (data.beta_per10<data.upper95).all()
data.to_csv(source/'Figure2_lung_cohorts.csv',index=False)
plt.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':10})
canvas,axes=plt.subplots(1,3,figsize=(8.2,3.65),sharex=True,sharey=True)
for ax,(cohort,group),color in zip(axes,data.groupby('cohort',sort=False),['#125071','#30736b','#75496b']):
    for y,(_,r) in enumerate(group.iterrows()):
        ax.errorbar(r.beta_per10,4-y,xerr=[[r.beta_per10-r.lower95],[r.upper95-r.beta_per10]],fmt='o',color=color,capsize=3,ms=5)
    ax.set_title(f'{cohort}\n(n = {int(group.n.iloc[0])})',fontsize=11)
    ax.axvline(0,color='0.45',lw=.85,ls='--',zorder=0)
    ax.set(xlim=(-.415,.415),ylim=(-.45,4.45),xticks=[-.4,0,.4],yticks=range(5),yticklabels=list(s12.index[::-1]))
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y',length=0)
canvas.subplots_adjust(left=.12,right=.97,top=.78,bottom=.25,wspace=.17)
canvas.supxlabel('Change in log₂ expression per 10 percentage points higher percent-predicted FEV₁',y=.08,fontsize=10)
for ext in ['pdf','png']:canvas.savefig(fig/('Figure2.'+ext),dpi=360)
