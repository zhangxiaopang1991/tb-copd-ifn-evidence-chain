"""Rerun association models from public-derived, RMA-normalised inputs."""
from pathlib import Path
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

BASE = Path(__file__).resolve().parent
GENES = ['IFIT1', 'IFI44', 'USP18', 'BST2', 'IFITM3']
a = pd.read_csv(BASE / 'inputs/GSE42057.csv')
b = pd.read_csv(BASE / 'inputs/GSE54837.csv')
cases = a.loc[a.group.eq('COPD')]
specs = [
    ('GSE42057_PBMC', 'GOLD_1_to_4', cases.loc[cases.finalgold.isin([1, 2, 3, 4])], 'fev1pp', 'C(smokcignow) + age + C(gender) + bmi', 84),
    ('GSE42057_PBMC', 'title_labelled_94', cases, 'fev1pp', 'C(smokcignow) + age + C(gender) + bmi', 94),
]
copd = b.loc[b.subject_type.eq('COPD Subjects')].copy()
assert copd.gold.value_counts().to_dict() == {'Stage II': 68, 'Stage III': 55, 'Stage IV': 13}
copd['gold_num'] = copd.gold.map({'Stage II': 2, 'Stage III': 3, 'Stage IV': 4})
specs.append(('GSE54837_whole_blood', 'original_adjusted_model', copd, 'gold_num', 'age + C(gender)', 136))
rows = []
for cohort, definition, subset, exposure, covariates, expected_n in specs:
    start = len(rows)
    for gene in GENES:
        fit = smf.ols(f'{gene} ~ {exposure} + {covariates}', data=subset).fit()
        assert int(fit.nobs) == expected_n
        ci = fit.conf_int().loc[exposure]
        rows.append(dict(cohort=cohort, sample_definition=definition, gene=gene, n=int(fit.nobs),
                         beta=fit.params[exposure], se=fit.bse[exposure], ci_low=ci.iloc[0], ci_high=ci.iloc[1], p=fit.pvalues[exposure]))
    q = multipletests([r['p'] for r in rows[start:]], method='fdr_bh')[1]
    for row, value in zip(rows[start:], q):
        row['q'] = value
result = pd.DataFrame(rows)
result.to_csv(BASE / 'rerun_results.csv', index=False)
reference = pd.read_csv(BASE / 'reference_S11.csv')
checks = []
for row in rows:
    match = reference.loc[reference.cohort.eq(row['cohort']) & reference.gene.eq(row['gene']) & reference.n.eq(row['n'])]
    assert len(match) == 1, (row['cohort'], row['gene'], row['n'])
    for field in ['beta', 'se', 'ci_low', 'ci_high', 'p', 'q']:
        if pd.isna(match.iloc[0][field]):
            assert field == 'se'
            continue
        difference = abs(row[field] - float(match.iloc[0][field]))
        checks.append(dict(cohort=row['cohort'], gene=row['gene'], n=row['n'], field=field, absolute_difference=difference))
        assert difference < 1e-10, checks[-1]
pd.DataFrame(checks).to_csv(BASE / 'comparison_checks.csv', index=False)
print(f'PASS: 15 models; n=84/94/136; {len(checks)} available numeric comparisons within 1e-10. Missing reference SE not compared.')
