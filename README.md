# Interferon-related gene expression and lung function in COPD

Analysis code and selected aggregate verification tables for the November 2026 manuscript revision. The earlier five-dataset release is at [v1.0.0-submission](https://github.com/zhangxiaopang1991/tb-copd-ifn-evidence-chain/tree/v1.0.0-submission). The current scripts add two external lung series and two blood cohorts, a STAT2 coverage audit, and three smoking codings in the specification grid. See Evidence/Analysis_changes.txt.

## Input sources

- NCBI GEO: GSE47460, GSE76925, GSE8581, GSE42057, GSE54837, GSE276060, GSE192483 and GSE136831.
- EMBL-EBI BioStudies/ArrayExpress: E-MTAB-17246.

The public code archive contains aggregate results and no participant-level input tables. Obtain source records from the original repositories. The computational supplement accompanying the manuscript includes the exact model inputs and the upstream sample-level audit tables. The legacy workflows were written for an existing analysis workspace and may require path configuration. No raw-to-results run for all nine datasets is claimed.

## Entry points

From Code/, run python Figure2.py to regenerate Figure 2 from the bundled S12 and S13 aggregate tables. Matplotlib and pandas are required.

Prepare the GSE47460 gene-level matrix and metadata in panel_inputs/GSE47460_gene_matrix.tsv and panel_inputs/sample_meta_final.tsv, then run python audit.py panel_inputs lung_model_rerun/reference_S6.csv and python rerun_grid.py. The archived setup scripts are in five_dataset_scripts/GSE47460/.

The model verification entry points are lung_model_rerun/rerun_lung_models.py and blood_model_rerun/rerun_blood_models.py. They expect lung_model_rerun/GSE47460_derived_inputs.csv and blood_model_rerun/inputs/GSE42057.csv plus GSE54837.csv, respectively. Use the reference tables, environment files and dependency pins in those folders.

For the additional lung series, run lung111/case_only_check.py with the GSE76925 series matrix and GPL10558 annotation in its GEO input directory. Run lung15/check_GSE8581.py with the GSE8581 series matrix, GPL570 annotation and clinical supplemental PDF from PMC2645534. For the STAT2 coverage audit, check_stat2.py uses GSE47460 platform annotations, probe maps, series matrices and sample SOFT files. The archived contextual workflows are documented in five_dataset_scripts/README.md.

The legacy MIT licence is included with the archived five-dataset scripts.
