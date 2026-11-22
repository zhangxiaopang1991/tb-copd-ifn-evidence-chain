from pathlib import Path
import sys,re
import pandas as pd,numpy as np
data=Path(sys.argv[1]);out=Path(__file__).resolve().parent/'Evidence'
frame=pd.read_csv(out/'lung_model_public_inputs.csv').set_index('gsm')
genes=frame.columns[-20:]
maps={}
for plat in frame.gpl.unique():
    tab=pd.read_csv(data/(plat.lower()+'_id2gene.tsv'),sep='\t',header=None,dtype=str)
    maps[plat]=dict(tab[tab[1].isin(genes)].iloc[:,:2].values)
flags=[]
for gsm,r in frame.iterrows():
    lines=(data/'soft'/(gsm+'.txt')).read_text(encoding='utf8',errors='replace')
    fields=lines.split('!sample_table_begin')[0]
    pack=bool(re.search(r'pack[\s_/-]*years?',fields,re.I))
    values={};mapping=maps[r.gpl]
    for line in lines.split('!sample_table_begin')[1].split('!sample_table_end')[0].splitlines():
        pair=line.split('\t')
        if len(pair)>1 and pair[0].strip('"') in mapping:
            g=mapping[pair[0].strip('"')]
            try:v=float(pair[1].strip('"'))
            except ValueError:continue
            values[g]=max(values.get(g,-np.inf),v)
    assert all(np.isclose(values[g],r[g],rtol=0,atol=1e-9) for g in genes),gsm
    flags.append([gsm,len(values),pack])
pd.DataFrame(flags,columns=['gsm','matched_genes','pack_year_label_in_header']).to_csv(out/'raw_sample_checks.csv',index=False)
print(len(flags),sum(i[1] for i in flags),sum(i[2] for i in flags))

