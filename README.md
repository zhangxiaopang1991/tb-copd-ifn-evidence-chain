# TB-COPD interferon evidence chain: reproducible analysis package

Scope: this package contains the analysis and figure-generation scripts used for the reported five-dataset reconstruction. It does not redistribute the original public datasets and is not yet a one-command workflow. Reproduction requires downloading the accession-level source data and preserving the input filenames and directory layout expected by each script.

Execution order:
1. analyze_gse276060.R
2. recheck_EMTAB17246.R
3. GSE192483/analyze_gse192483.py
4. GSE47460/analyze_gse47460.R and sensitivity_ifn_covariates.R
5. GSE136831/aggregate_pseudobulk_counts.py, run_pseudobulk_edgeR.R, analyze_cell_composition.py and test_disease_celltype_interaction.py
6. figures/generate_figures.R

Software: R 4.6.1 (edgeR 4.10.1, limma 3.68.4, fgsea 1.38.0); Python 3.12.14 (scanpy 1.12.4, anndata 0.13.3).

Inference uses sample, donor, patient or subject units; cells are not independent replicates.

Key analysis boundaries:
- The five datasets are analysed separately. No pooled cross-platform meta-analysis is performed.
- E-MTAB-17246 uses `none_Air` as the reference level. The Mtb-by-smoke interaction is the difference-in-differences `(Mtb_Smoke - Mtb_Air) - (none_Smoke - none_Air)`.
- GSE192483 cell-level summaries are descriptive; formal lesion comparisons use patient-level H/L pairs.
- GSE47460 disease and FEV1 models use complete cases for included covariates. Remaining expression-matrix missing values after gene filtering are row-mean imputed.
- GSE136831 score, count-based pseudobulk and composition analyses answer different questions and should not be interpreted as interchangeable validation layers.

Frozen submission version: `v1.0.0-submission`.

License: MIT. See `LICENSE`.
