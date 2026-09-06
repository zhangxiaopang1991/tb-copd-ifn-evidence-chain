suppressPackageStartupMessages({ library(limma); library(fgsea) })

args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)][1])
if (!file.exists("GSE47460_gene_matrix.tsv") && !is.na(script_file)) {
  setwd(dirname(script_file))
}

mat <- read.delim("GSE47460_gene_matrix.tsv", row.names = 1, check.names = FALSE)
meta <- read.delim("sample_meta_final.tsv", stringsAsFactors = FALSE, check.names = FALSE)
rownames(meta) <- meta$gsm
meta <- meta[colnames(mat), , drop = FALSE]
keep <- meta$`disease state` %in% c("Control", "Chronic Obstructive Lung Disease")
meta <- meta[keep, , drop = FALSE]
mat <- as.matrix(mat[, rownames(meta), drop = FALSE])
storage.mode(mat) <- "double"
mat[mat == ""] <- NA
mat <- mat[rowMeans(!is.na(mat)) >= 0.80, , drop = FALSE]
for (i in seq_len(nrow(mat))) {
  x <- mat[i, ]
  if (anyNA(x)) x[is.na(x)] <- mean(x, na.rm = TRUE)
  mat[i, ] <- x
}

meta$disease <- factor(ifelse(meta$`disease state` == "Control", "Control", "COPD"),
                       levels = c("Control", "COPD"))
meta$platform <- factor(meta$gpl)
smoking_raw <- trimws(meta$`smoker?`)
smoking_group <- rep(NA_character_, nrow(meta))
known_smoking <- !is.na(smoking_raw) & smoking_raw != ""
smoking_group[known_smoking & grepl("Never", smoking_raw, fixed = TRUE)] <- "Never"
smoking_group[known_smoking & grepl("Current", smoking_raw, fixed = TRUE)] <- "Current"
smoking_group[known_smoking & is.na(smoking_group)] <- "Ever"
meta$smoking <- factor(smoking_group, levels = c("Never", "Ever", "Current"))
meta$age_num <- suppressWarnings(as.numeric(meta$age))
meta$sex_factor <- factor(meta$Sex)

gene_sets <- list(
  PTLD_IFN = c("IFI6","OAS1","OAS2","ISG15","MX1","MX2","IFIT1","IFIT2","IFIT3","STAT1","STAT2","IRF7","IRF9","GBP1","GBP2","GBP5","OASL","IFI44","IFI27","IFI35","IFITM1","IFITM3","RSAD2","USP18","XAF1","BST2","EPSTI1","DDX58","IFIH1","EIF2AK2"),
  IFNa = c("ISG15","MX1","OAS1","OAS2","IFIT1","IFIT2","IFIT3","IFI6","IFI27","IFI44","IFI44L","IFITM1","IFITM2","IFITM3","RSAD2","USP18","XAF1","BST2","EPSTI1","EIF2AK2","GBP1","GBP2","DDX58","IFIH1","IRF7","IRF9","STAT1","STAT2","SAMD9","SAMD9L","PARP9","DTX3L","TRIM22","TRIM25","TNFSF10","ISG20","OASL"),
  IFNg = c("STAT1","STAT2","IRF1","IRF7","IRF9","GBP1","GBP2","GBP5","CXCL9","CXCL10","CXCL11","HLA-A","HLA-B","HLA-C","TAP1","TAP2","B2M","PSMB8","PSMB9","PSME1","PSME2","NLRC5","CIITA","IFITM1","IFITM3","OAS1","OAS2","OASL","MX1","MX2","ISG15","ISG20","IFI35","SOCS1","SOCS3","JAK1","JAK2","CCL2","CCL5","CCL7")
)
gene_sets <- lapply(gene_sets, function(x) unique(x[x %in% rownames(mat)]))

specs <- list(
  base_platform = list(formula = ~ disease + platform, vars = c("disease", "platform")),
  plus_smoking = list(formula = ~ disease + platform + smoking, vars = c("disease", "platform", "smoking")),
  full_covariates = list(formula = ~ disease + platform + smoking + age_num + sex_factor,
                         vars = c("disease", "platform", "smoking", "age_num", "sex_factor"))
)
core <- c("IFI6","OAS1","OAS2","ISG15","MX1","IFIT1","STAT1","IRF9","GBP5","OASL","IFI27","RSAD2","BST2","IFITM1","EIF2AK2")
fg_out <- list(); core_out <- list()

for (model_name in names(specs)) {
  spec <- specs[[model_name]]
  complete <- complete.cases(meta[, spec$vars, drop = FALSE])
  md <- droplevels(meta[complete, , drop = FALSE])
  mm <- mat[, complete, drop = FALSE]
  design <- model.matrix(spec$formula, data = md)
  fit <- eBayes(lmFit(mm, design))
  stats <- fit$t[, "diseaseCOPD"]
  names(stats) <- rownames(mm)
  fg <- fgsea(pathways = gene_sets, stats = stats, nperm = 20000, nproc = 1)
  fg$model <- model_name
  fg$n_samples <- nrow(md)
  fg$n_COPD <- sum(md$disease == "COPD")
  fg$n_Control <- sum(md$disease == "Control")
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ",", FUN.VALUE = character(1))
  fg_out[[model_name]] <- fg[, c("model","pathway","n_samples","n_COPD","n_Control","pval","padj","ES","NES","size","leadingEdge")]

  tt <- topTable(fit, coef = "diseaseCOPD", number = Inf, sort.by = "none")
  for (g in intersect(core, rownames(tt))) {
    core_out[[length(core_out) + 1L]] <- data.frame(model = model_name, gene = g,
      n_samples = nrow(md), n_COPD = sum(md$disease == "COPD"), n_Control = sum(md$disease == "Control"),
      coef = tt[g, "logFC"], pval = tt[g, "P.Value"], padj = tt[g, "adj.P.Val"])
  }
  cat(model_name, ": n=", nrow(md), " (COPD ", sum(md$disease == "COPD"), ", Control ", sum(md$disease == "Control"), ")\n", sep = "")
  print(fg[, c("pathway","NES","padj")])
}

write.csv(do.call(rbind, fg_out), "GSE47460_IFN_covariate_sensitivity_fgsea.csv", row.names = FALSE)
write.csv(do.call(rbind, core_out), "GSE47460_IFN_covariate_sensitivity_core_genes.csv", row.names = FALSE)
write.table(data.frame(variable = c("smoking_missing","COPD_total","Control_total"),
                       value = c(sum(is.na(meta$smoking)), sum(meta$disease == "COPD"), sum(meta$disease == "Control"))),
            "GSE47460_IFN_covariate_sensitivity_metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote GSE47460 covariate sensitivity outputs.\n")
