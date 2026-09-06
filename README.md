# TB-COPD interferon evidence chain: reproducible analysis package

This package contains the analysis scripts for the five-dataset interferon
evidence-chain study. Raw public data are not included. Each script reads
inputs documented in `dataset_manifest.tsv`.

## Datasets

| Accession | Role |
|---|---|
| GSE276060 | Infection-stage bulk RNA-seq |
| E-MTAB-17246 | Paired airway epithelial smoke/mycobacteria model |
| GSE192483 | Lesion-level single-cell RNA-seq |
| GSE47460 | Chronic-stage bulk microarray with FEV1 metadata |
| GSE136831 | Chronic-stage single-cell RNA-seq |

## Gene sets

`gene_sets.tsv` defines three project-prespecified modules: `PTLD_IFN`,
`IFNa`, and `IFNg`. They are project-defined and are not official MSigDB
Hallmark sets. Some earlier bulk scripts hard-code Hallmark-like IFN-alpha and
IFN-gamma lists; these should be reconciled against `gene_sets.tsv` before
final submission.

## Software

- Python 3.12.14, scanpy 1.12.4, anndata 0.13.3, scipy, statsmodels, pandas, numpy
- R 4.6.1, edgeR 4.10.1, limma 3.68.4, fgsea 1.38.0, Matrix 1.7.5, dplyr, tidyr

## Execution order

1. Download and parse public inputs for each dataset.
2. Run per-dataset analyses.
3. Generate figures and supplementary tables.
4. Run metadata and sample-flow QA audits.

See individual script headers for expected working directory and input paths.

## Statistical units

Formal inference uses sample, donor, patient, or subject units. Individual
cells are not treated as independent replicates.
