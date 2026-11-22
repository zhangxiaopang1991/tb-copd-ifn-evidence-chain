# Five-dataset source snapshot status

This folder is copied from the local `tb-copd-data/submission_repo` and includes the five-dataset scripts, accession manifest and project-defined gene sets. It supports the study's original five-dataset reconstruction: GSE276060, E-MTAB-17246, GSE192483, GSE47460 and GSE136831.

The two peripheral-blood cohorts GSE42057 and GSE54837 are covered separately by `../blood_model_rerun` and the workspace-dependent blood scripts. The local five-dataset source snapshot is not a one-command pipeline: scripts expect accession-derived input files and some working-directory layouts. Raw public data are not copied here.

The gene sets are project-defined exploratory modules, not official MSigDB Hallmark sets. Before public release, reconcile any hard-coded lists in legacy scripts against `gene_sets.tsv`, test the full execution order in a clean environment, and verify redistribution terms for every derived input. No independent biological validation is implied.
