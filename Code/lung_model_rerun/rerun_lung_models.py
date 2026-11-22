"""Reproduce the 20-gene lung OLS model from processed inputs, not raw arrays."""
from pathlib import Path
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

base = Path(__file__).resolve().parent
frame = pd.read_csv(base / 'GSE47460_derived_inputs.csv', index_col='gsm')
reference = pd.read_csv(base / 'reference_S6.csv').set_index('gene')
assert frame.index.is_unique and len(reference) == 20 and reference.index.is_unique
rows = []
for gene in reference.index:
    current = frame.copy()
    current['expression'] = current[gene]
    fit = smf.ols('expression ~ fev + C(gpl) + C(smoking) + age_num + C(Sex)', data=current).fit()
    assert fit.nobs == 215
    ci = fit.conf_int().loc['fev']
    rows.append(dict(gene=gene, n=int(fit.nobs), beta=fit.params['fev'],
                     se=fit.bse['fev'], p=fit.pvalues['fev'], ci_lo=ci.iloc[0], ci_hi=ci.iloc[1]))
results = pd.DataFrame(rows).set_index('gene')
results['q_BH'] = multipletests(results.p, method='fdr_bh')[1]
checks = []
for gene in results.index:
    for field in ['beta', 'se', 'p', 'ci_lo', 'ci_hi', 'q_BH']:
        difference = abs(results.loc[gene, field] - reference.loc[gene, field])
        assert difference < 1e-10, (gene, field, difference)
        checks.append(dict(gene=gene, field=field, absolute_difference=difference))
assert set(results.index[results.q_BH < .05]) == {'IFIT1', 'IFI44', 'USP18', 'BST2', 'IFITM3'}
results.to_csv(base / 'rerun_results.csv')
pd.DataFrame(checks).to_csv(base / 'comparison_checks.csv', index=False)
print('PASS: 20 lung models; n=215; 120 comparisons within 1e-10; five BH-significant genes.')
