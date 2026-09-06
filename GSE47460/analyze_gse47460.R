suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
  library(ggplot2)
})
.script_args <- commandArgs(trailingOnly = FALSE)
.script_file <- sub("^--file=", "", .script_args[grep("^--file=", .script_args)][1])
if (!file.exists("GSE47460_gene_matrix.tsv") && !is.na(.script_file)) {
  setwd(dirname(.script_file))
}

# ---------- 读数据 ----------
mat <- read.delim("GSE47460_gene_matrix.tsv", row.names = 1, check.names = FALSE)
meta <- read.delim("sample_meta_final.tsv", stringsAsFactors = FALSE, check.names = FALSE)
rownames(meta) <- meta$gsm
meta <- meta[colnames(mat), ]

cat("矩阵维度:", nrow(mat), "基因 x", ncol(mat), "样本\n")
cat("疾病分布:\n"); print(table(meta$`disease state`))
cat("平台分布:\n"); print(table(meta$gpl))

# ---------- 过滤：仅 COPD vs Control（排除 ILD）----------
keep_s <- meta$`disease state` %in% c("Control", "Chronic Obstructive Lung Disease")
meta <- meta[keep_s, ]
mat <- mat[, rownames(meta)]
meta$disease <- factor(ifelse(meta$`disease state` == "Control", "Control", "COPD"),
                       levels = c("Control", "COPD"))
meta$gpl <- factor(meta$gpl)
cat("\n过滤后:", nrow(meta), "样本 (COPD vs Control)\n")
cat("COPD/Control x 平台交叉表:\n"); print(table(meta$disease, meta$gpl))

# ---------- 表达矩阵数值化 + 缺失过滤 ----------
mat <- as.matrix(mat)
storage.mode(mat) <- "double"
mat[mat == ""] <- NA
cat("\n缺失率: ", round(mean(is.na(mat)) * 100, 2), "%\n")

# 保留在 >= 80% 样本中有值的基因
keep_g <- rowMeans(!is.na(mat)) >= 0.80
mat <- mat[keep_g, ]
cat("保留基因(缺失<20%):", nrow(mat), "\n")

# 剩余 NA 用基因均值填补（芯片数据，少量缺失可接受）
for (i in 1:nrow(mat)) {
  x <- mat[i, ]
  if (any(is.na(x))) x[is.na(x)] <- mean(x, na.rm = TRUE)
  mat[i, ] <- x
}

# ---------- limma：COPD vs Control，gpl 作协变量 ----------
design <- model.matrix(~ 0 + disease + gpl, data = meta)
colnames(design) <- make.names(colnames(design))
fit <- lmFit(mat, design)
cont <- makeContrasts(COPD_vs_Control = diseaseCOPD - diseaseControl, levels = design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)
de <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
de$gene <- rownames(de)
write.csv(de, "GSE47460_DEG_COPD_vs_Control.csv", row.names = FALSE)
cat("\nCOPD vs Control: 显著 DEG 数 (adj.P<0.05):", sum(de$adj.P.Val < 0.05), "\n")
cat("上调:", sum(de$adj.P.Val < 0.05 & de$logFC > 0), " 下调:", sum(de$adj.P.Val < 0.05 & de$logFC < 0), "\n")

# ---------- fgsea ----------
gs_sym <- list(
  PTLD_IFN = c("IFI6","OAS1","OAS2","ISG15","MX1","MX2","IFIT1","IFIT2","IFIT3",
                           "STAT1","STAT2","IRF7","IRF9","GBP1","GBP2","GBP5","OASL","IFI44",
                           "IFI27","IFI35","IFITM1","IFITM3","RSAD2","USP18","XAF1","BST2",
                           "EPSTI1","DDX58","IFIH1","EIF2AK2"),
  IFNa = c("ISG15","MX1","OAS1","OAS2","IFIT1","IFIT2","IFIT3","IFI6","IFI27","IFI44",
                    "IFI44L","IFITM1","IFITM2","IFITM3","RSAD2","USP18","XAF1","BST2","EPSTI1",
                    "EIF2AK2","GBP1","GBP2","DDX58","IFIH1","IRF7","IRF9","STAT1","STAT2",
                    "SAMD9","SAMD9L","PARP9","DTX3L","TRIM22","TRIM25","TNFSF10","ISG20","OASL"),
  IFNg = c("STAT1","STAT2","IRF1","IRF7","IRF9","GBP1","GBP2","GBP5",
                    "CXCL9","CXCL10","CXCL11","HLA-A","HLA-B","HLA-C","TAP1","TAP2","B2M",
                    "PSMB8","PSMB9","PSME1","PSME2","NLRC5","CIITA","IFITM1","IFITM3","OAS1",
                    "OAS2","OASL","MX1","MX2","ISG15","ISG20","IFI35","SOCS1","SOCS3",
                    "JAK1","JAK2","CCL2","CCL5","CCL7")
)
gs <- lapply(gs_sym, function(s) unique(s[s %in% rownames(mat)]))
gs <- gs[sapply(gs, length) >= 3]
cat("\n基因集大小(交集后):\n"); print(sapply(gs, length))

stats <- fit2$t[, 1]
names(stats) <- rownames(mat)
stats <- stats[is.finite(stats)]
fg <- fgsea(pathways = gs, stats = stats, nperm = 20000, nproc = 1)
fg <- fg[order(fg$pval), ]
fg$leadingEdge <- sapply(fg$leadingEdge, function(x) paste(x, collapse = ","))
write.csv(fg, "GSE47460_fgsea_COPD_vs_Control.csv", row.names = FALSE)
cat("\n== fgsea 结果 ==\n")
print(fg[, c("pathway","pval","padj","ES","NES","size")])

# ---------- 点图 ----------
plot_df <- fg
plot_df$signed_logFDR <- sign(plot_df$ES) * (-log10(pmax(plot_df$padj, 1e-10)))
plot_df$pathway <- factor(plot_df$pathway, levels = rev(plot_df$pathway))
p <- ggplot(plot_df, aes(x = signed_logFDR, y = pathway)) +
  geom_point(aes(size = size, color = padj < 0.05), show.legend = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey60"), name = "FDR<0.05") +
  labs(x = "signed -log10(FDR)", y = NULL,
       title = "GSE47460 LTRC lung tissue: COPD vs Control (IFN gene sets)",
       subtitle = "positive = up in COPD; n = COPD 220 / Control 108, platform-adjusted") +
  theme_bw(base_size = 10)
ggsave("GSE47460_fgsea_dotplot.png", p, width = 9, height = 3.5, dpi = 150)

# ---------- FEV1 连续关联（核心 IFN 基因）----------
fev1 <- suppressWarnings(as.numeric(meta$`%predicted fev1 (pre-bd)`))
has_fev1 <- !is.na(fev1) & meta$disease == "COPD"
cat("\nCOPD 样本有 FEV1 记录:", sum(has_fev1), "\n")
if (sum(has_fev1) >= 20) {
  m2 <- meta[has_fev1, ]
  mm <- mat[, rownames(m2)]
  fev <- fev1[has_fev1]
  m2$smoking_raw <- trimws(m2$`smoker?`)
  smoking_group <- rep(NA_character_, nrow(m2))
  known_smoking <- !is.na(m2$smoking_raw) & m2$smoking_raw != ""
  smoking_group[known_smoking & grepl("Never", m2$smoking_raw, fixed = TRUE)] <- "Never"
  smoking_group[known_smoking & grepl("Current", m2$smoking_raw, fixed = TRUE)] <- "Current"
  smoking_group[known_smoking & is.na(smoking_group)] <- "Ever"
  m2$smoking <- factor(smoking_group, levels = c("Never", "Ever", "Current"))
  m2$age_num <- suppressWarnings(as.numeric(m2$age))
  m2$sex_factor <- factor(m2$Sex)
  core <- c("IFI6","OAS1","OAS2","ISG15","MX1","MX2","IFIT1","STAT1","STAT2","IRF7","IRF9",
            "GBP1","GBP2","GBP5","OASL","RSAD2","USP18","IFI44","IFI27","IFITM3","BST2")
  core <- core[core %in% rownames(mm)]
  fev_specs <- list(
    platform = list(formula = ~ fev + gpl, vars = c("fev", "gpl")),
    plus_smoking = list(formula = ~ fev + gpl + smoking, vars = c("fev", "gpl", "smoking")),
    full_covariates = list(formula = ~ fev + gpl + smoking + age_num + sex_factor,
                           vars = c("fev", "gpl", "smoking", "age_num", "sex_factor"))
  )
  fev_out <- list()
  for (model_name in names(fev_specs)) {
    spec <- fev_specs[[model_name]]
    md <- m2
    md$fev <- fev
    complete <- complete.cases(md[, spec$vars, drop = FALSE])
    md <- droplevels(md[complete, , drop = FALSE])
    mm2 <- mm[, complete, drop = FALSE]
    design2 <- model.matrix(spec$formula, data = md)
    fitf <- eBayes(lmFit(mm2, design2))
    fev_res <- topTable(fitf, coef = "fev", number = Inf, sort.by = "none")
    fev_core <- fev_res[core, c("logFC","P.Value","adj.P.Val")]
    fev_core$gene <- rownames(fev_core)
    fev_core$model <- model_name
    fev_core$n_samples <- nrow(md)
    colnames(fev_core)[1:3] <- c("coef_fev1","P.Value","adj.P.Val")
    fev_out[[model_name]] <- fev_core
    cat("\n== FEV1 与核心 IFN 基因关联：", model_name,
        "（COPD 内，正系数=FEV1越高表达越高）==\n", sep = "")
    print(fev_core[order(fev_core$P.Value), c("gene","coef_fev1","P.Value","adj.P.Val")])
  }
  write.csv(do.call(rbind, fev_out), "GSE47460_FEV1_IFN_core_genes.csv", row.names = FALSE)
}

cat("\nDone.\n")
