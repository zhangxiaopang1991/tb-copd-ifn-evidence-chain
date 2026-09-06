suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

OUT <- "supplementary"
dir.create(OUT, showWarnings = FALSE)

write_tbl <- function(x, name) {
  write.csv(x, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
}

# S1: dataset and analysis-unit registry.
s1 <- data.frame(
  dataset = c("GSE276060", "E-MTAB-17246", "GSE192483", "GSE47460", "GSE136831"),
  phase = c("infection", "environmental modifier", "lesion-level localization",
            "chronic bulk", "chronic cell-level"),
  assay = c("bulk RNA-seq", "paired bulk RNA-seq", "10x single-cell",
            "two-platform microarray", "10x single-cell"),
  statistical_unit = c("sample", "donor", "patient x celltype pair",
                       "sample/subject", "subject x celltype"),
  primary_role = c("infection-stage anchor", "smoke-amplification evidence",
                   "directional localization", "covariate-sensitive bulk and FEV1",
                   "state and composition decomposition"),
  stringsAsFactors = FALSE
)
write_tbl(s1, "Supplementary_Table_S1_dataset_registry")

# S2: full infection-stage pathway table.
s2 <- read.csv("results/GSE276060_fgsea.csv", check.names = FALSE) |>
  arrange(contrast, padj, pathway)
write_tbl(s2, "Supplementary_Table_S2_GSE276060_fgsea")

# S3: full paired exposure pathway table.
s3 <- read.csv("results/EMTAB17246_fgsea_PTLD_IFN.csv", check.names = FALSE) |>
  arrange(contrast, padj, pathway)
write_tbl(s3, "Supplementary_Table_S3_EMTAB17246_fgsea")

# S4: all available donor-aware lesion comparisons.
s4 <- read.csv("GSE192483/results/GSE192483_IFN_score_H_vs_L_pseudobulk_paired.csv",
               check.names = FALSE) |>
  arrange(padj, celltype, score)
write_tbl(s4, "Supplementary_Table_S4_GSE192483_paired")

# S5: all bulk pathway sensitivity models.
s5 <- read.csv("GSE47460/GSE47460_IFN_covariate_sensitivity_fgsea.csv",
               check.names = FALSE) |>
  arrange(pathway, model)
write_tbl(s5, "Supplementary_Table_S5_GSE47460_pathway_sensitivity")

# S6: all prespecified core-gene FEV1 sensitivity results.
s6 <- read.csv("GSE47460/GSE47460_FEV1_IFN_core_genes.csv",
               check.names = FALSE) |>
  arrange(gene, model)
write_tbl(s6, "Supplementary_Table_S6_GSE47460_FEV1")

# S7: all subject-level score comparisons, including non-significant results.
s7 <- read.csv("GSE136831/results/GSE136831_pseudobulk_celltype_test.csv",
               check.names = FALSE) |>
  arrange(padj, celltype, score)
write_tbl(s7, "Supplementary_Table_S7_GSE136831_score")

# S8: all count-based pathway-by-celltype tests.
s8 <- read.csv("GSE136831/results/GSE136831_edgeR_IFN_summary.csv",
               check.names = FALSE) |>
  arrange(global_padj, celltype, pathway)
write_tbl(s8, "Supplementary_Table_S8_GSE136831_edgeR")

# S9: all composition tests after Multiplet exclusion and zero filling.
s9 <- read.csv("GSE136831/results/GSE136831_celltype_composition_test.csv",
               check.names = FALSE) |>
  arrange(padj, celltype)
write_tbl(s9, "Supplementary_Table_S9_GSE136831_composition")

# S10: exploratory subject-clustered disease-by-celltype interaction.
s10 <- read.csv(
  "GSE136831/results/GSE136831_disease_celltype_interaction_summary.csv",
  check.names = FALSE
) |>
  arrange(padj_BH_across_scores, score)
write_tbl(s10, "Supplementary_Table_S10_GSE136831_interaction")

writeLines(c(
  "Supplementary table generation record",
  "Generated from the current local result tables.",
  "S1 dataset registry",
  "S2 GSE276060 full fgsea results",
  "S3 E-MTAB-17246 full fgsea results",
  "S4 GSE192483 all donor-aware paired comparisons",
  "S5 GSE47460 pathway covariate sensitivity",
  "S6 GSE47460 FEV1 core-gene sensitivity",
  "S7 GSE136831 all subject-level score comparisons",
  "S8 GSE136831 all count-based edgeR pathway-celltype tests",
  "S9 GSE136831 all composition tests",
  "S10 GSE136831 exploratory disease-by-celltype interaction summary",
  "All tables retain non-significant comparisons for transparent review."
), file.path(OUT, "supplementary_table_QA_notes.txt"))

cat("Generated S1-S10 supplementary tables in", OUT, "\n")
