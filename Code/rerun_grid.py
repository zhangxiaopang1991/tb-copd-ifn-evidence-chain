from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

root = Path(__file__).resolve().parent
data = root / 'panel_inputs'
genes = pd.read_csv(root / 'lung_model_rerun/reference_S6.csv').gene.tolist()
meta = pd.read_csv(data / 'sample_meta_final.tsv', sep='\t').set_index('gsm')
matrix = pd.read_csv(data / 'GSE47460_gene_matrix.tsv', sep='\t', index_col=0)
meta = meta.loc[meta['disease state'] == 'Chronic Obstructive Lung Disease'].copy()
meta['fev'] = pd.to_numeric(meta['%predicted fev1 (pre-bd)'], errors='coerce')
meta['age_num'] = pd.to_numeric(meta.age, errors='coerce')
raw = meta['smoker?'].fillna('').astype(str).str.strip()
meta['smoking3'] = raw.map(lambda x: 'Never' if 'Never' in x else 'Current' if 'Current' in x else 'Ever' if x else None)
meta['smoking2'] = raw.map({'1-Current':0,'2-Ever (>100)':1,'3-Never':1})
meta['ever'] = raw.map({'1-Current':1,'2-Ever (>100)':1,'3-Never':0})
for gene in genes:
    meta[gene] = pd.to_numeric(matrix.loc[gene].reindex(meta.index), errors='coerce')
meta = meta.dropna(subset=['fev', 'age_num', 'smoking3', 'smoking2', 'Sex', 'gpl'] + genes)
records = []
assert meta.smoking3.value_counts().to_dict()=={'Ever':189,'Current':14,'Never':12}
for encoding,term in [('current_vs_other','smoking2'),('never_vs_ever','ever'),('three_category','C(smoking3)')]:
    terms = ['C(gpl)', term, 'age_num', 'C(Sex)']
    for size in range(5):
        for subset in combinations(terms, size):
            for winsor in [False, True]:
                for gene in genes:
                    frame = meta.copy()
                    frame['expression'] = frame[gene]
                    if winsor:
                        for field in ['fev', 'expression']:
                            frame[field] = frame[field].clip(*frame[field].quantile([.01, .99]).tolist())
                    formula = 'expression ~ fev' + ((' + ' + ' + '.join(subset)) if subset else '')
                    fit = smf.ols(formula, data=frame).fit()
                    records.append(dict(encoding=encoding, gene=gene, covariates=' + '.join(subset) or 'none', winsor=winsor, n=int(fit.nobs), beta=fit.params.fev, p=fit.pvalues.fev))
results = pd.DataFrame(records)
rows = []
for (encoding, gene), group in results.groupby(['encoding', 'gene']):
    sign = np.sign(group.beta.median())
    passed = (group.p < .05) & (np.sign(group.beta) == sign)
    rows.append(dict(encoding=encoding, gene=gene, n_specs=len(group), n_pass=int(passed.sum()), fraction=passed.mean(), all_pass=bool(passed.all())))
summary = pd.DataFrame(rows)
results.to_csv(root / 'specification_grid_audit_results.csv', index=False)
summary.to_csv(root / 'specification_grid_audit_summary.csv', index=False)
