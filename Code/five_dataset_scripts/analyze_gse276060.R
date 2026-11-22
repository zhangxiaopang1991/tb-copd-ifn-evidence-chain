suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(fgsea)
  library(ggplot2)
})

.script_args <- commandArgs(trailingOnly = FALSE)
.script_file <- sub("^--file=", "", .script_args[grep("^--file=", .script_args)][1])
if (!is.na(.script_file)) setwd(dirname(normalizePath(.script_file)))
dir.create("results", showWarnings = FALSE)

# ---------- 读入 ----------
counts <- read.delim("GSE276060/GSE276060_counts.tsv", row.names = 1, check.names = FALSE)
meta <- read.delim("GSE276060/GSE276060_meta.tsv", stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
counts <- counts[, meta$sample]
stopifnot(all(colnames(counts) == meta$sample))

cat("原始矩阵:", nrow(counts), "基因 x", ncol(counts), "样本\n")
cat("分组计数:\n"); print(table(meta$group))

# 组因子（顺序固定，Ctl_Lung 作参照）
meta$group <- factor(meta$group, levels = c("Ctl_Lung", "Ctl_Tonsil", "MTB", "NTM"))

# ---------- DGEList + 过滤 + voom ----------
dge <- DGEList(counts = counts, samples = meta)
keep <- filterByExpr(dge, group = meta$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
cat("过滤后:", nrow(dge), "基因\n")

design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- make.names(colnames(design))
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)

cont.matrix <- makeContrasts(
  MTB_vs_Ctl   = groupMTB - groupCtl_Lung,
  NTM_vs_Ctl   = groupNTM - groupCtl_Lung,
  MTB_vs_NTM   = groupMTB - groupNTM,
  Myco_vs_Ctl  = (groupMTB + groupNTM)/2 - groupCtl_Lung,
  levels = design
)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# ---------- 输出 DEG 表（含 symbol 列）----------
for (cn in colnames(cont.matrix)) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  write.csv(tt, sprintf("results/GSE276060_DEG_%s.csv", cn), row.names = FALSE)
}

# ---------- 基因集（symbol 级别）----------
gs <- list(
  PTLD_IFN = c("IFI6","OAS1","OAS2",                       # Kenyan 论文 STRING 核心
                           "ISG15","MX1","MX2","IFIT1","IFIT2","IFIT3",
                           "STAT1","STAT2","IRF7","IRF9","GBP1","GBP2","GBP5",
                           "OASL","IFI44","IFI27","IFI35","IFITM1","IFITM3",
                           "RSAD2","USP18","XAF1","BST2","EPSTI1","DDX58","IFIH1","EIF2AK2"),
  IFNa = c("ISG15","MX1","OAS1","OAS2","IFIT1","IFIT2","IFIT3",
                    "IFI6","IFI27","IFI44","IFI44L","IFITM1","IFITM2","IFITM3",
                    "RSAD2","USP18","XAF1","BST2","EPSTI1","EIF2AK2","GBP1","GBP2",
                    "DDX58","IFIH1","IRF7","IRF9","STAT1","STAT2","SAMD9","SAMD9L",
                    "PARP9","DTX3L","TRIM22","TRIM25","TNFSF10","ISG20","OASL"),
  IFNg = c("STAT1","STAT2","IRF1","IRF7","IRF9",
                    "GBP1","GBP2","GBP5","CXCL9","CXCL10","CXCL11",
                    "HLA-A","HLA-B","HLA-C","TAP1","TAP2","B2M",
                    "PSMB8","PSMB9","PSME1","PSME2","NLRC5","CIITA",
                    "IFITM1","IFITM3","OAS1","OAS2","OASL","MX1","MX2","ISG15","ISG20","IFI35",
                    "SOCS1","SOCS3","JAK1","JAK2","CCL2","CCL5","CCL7"),
  IFN_antigen_presentation = c("STAT1","STAT2","IRF1","IRF7","IRF9","TAP1","TAP2","B2M",
                               "HLA-A","HLA-B","HLA-C","HLA-E","NLRC5","CIITA","TAPBP",
                               "GBP1","GBP2","GBP4","GBP5","MX1","MX2","OAS1","OAS2","OASL",
                               "ISG15","ISG20","IFIT1","IFIT2","IFIT3","IFI44","IFI6","IFI27",
                               "IFIH1","DDX58","PSMB8","PSMB9","PSME1","PSME2","JAK1","JAK2","TYK2"),
  Inflammation_cytokines = c("IL1B","IL6","CXCL8","IL18","IL33","TNF",
                             "CCL2","CCL5","CCL20","CXCL1","CXCL10","CXCL11","CSF2","CSF3"),
  Oxidative_stress = c("G6PD","GCLC","GCLM","GSR","TKT","FTL","FTH1","HMOX1","NQO1","TXN","TXNRD1",
                       "SOD1","SOD2","CAT","GPX1","GPX2","PRDX1","NFE2L2","KEAP1","SLC7A11",
                       "GSTM1","ALDH3A1","CYP1B1","CYP4B1")
)

# 交集到矩阵内基因
gs <- lapply(gs, function(syms) unique(intersect(syms, rownames(v))))
gs <- gs[sapply(gs, length) >= 3]
cat("\n基因集大小（交集后）:\n"); print(sapply(gs, length))

# ---------- fgsea ----------
set.seed(20260906)
fgsea_all <- data.frame()
for (cn in colnames(cont.matrix)) {
  stats <- fit2$t[, cn]
  names(stats) <- rownames(fit2)
  stats <- stats[is.finite(stats)]
  res <- fgsea(pathways = gs, stats = stats, nperm = 20000, nproc = 1)
  res <- res[order(res$pval), ]
  res$contrast <- cn
  fgsea_all <- rbind(fgsea_all, res[, c("contrast","pathway","pval","padj","ES","NES","size","leadingEdge")])
}
fgsea_all$leadingEdge <- sapply(fgsea_all$leadingEdge, function(x) paste(x, collapse = ","))
write.csv(fgsea_all, "results/GSE276060_fgsea.csv", row.names = FALSE)

# ---------- 点图 ----------
plot_df <- fgsea_all
plot_df$signed_logFDR <- sign(plot_df$ES) * (-log10(pmax(plot_df$padj, 1e-10)))
plot_df$trend <- plot_df$padj < 0.25

p <- ggplot(plot_df, aes(x = signed_logFDR, y = pathway)) +
  geom_point(aes(size = size, alpha = trend), color = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ contrast, ncol = 4) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35)) +
  labs(x = "signed -log10(FDR) by enrichment direction",
       y = NULL,
       title = "GSE276060 (human lung) pathway enrichment (fgsea)",
       subtitle = "FDR<0.25 highlighted; positive NES = up in numerator") +
  theme_bw(base_size = 9) + theme(legend.position = "top")
ggsave("results/GSE276060_fgsea_dotplot.png", p, width = 12, height = 5.5, dpi = 150)

cat("\n== fgsea 显著/趋势 (FDR<0.25) ==\n")
print(fgsea_all[fgsea_all$padj < 0.25, c("contrast","pathway","NES","pval","padj","size")])

# ---------- 核心 IFN 基因表达（MTB vs Ctl_Lung）----------
core <- c("IFI6","OAS1","OAS2","ISG15","MX1","MX2","IFIT1","IFIT3","STAT1","STAT2",
          "IRF7","IRF9","GBP1","GBP2","GBP5","OASL","IFI44","IFI27","RSAD2","USP18","BST2")
core <- intersect(core, rownames(v))
if (length(core) >= 3) {
  # 只取肺样本
  lung <- meta$sample[meta$group %in% c("Ctl_Lung","MTB","NTM")]
  E <- v$E[core, lung, drop = FALSE]
  ann <- data.frame(group = meta[lung, "group"], row.names = lung)
  # 保存热图矩阵
  write.csv(as.data.frame(E), "results/GSE276060_IFN_genes_expr.csv")
  # 简单热图（无 pheatmap 时用 ggplot 热图）
  library(reshape2)
  Em <- melt(as.matrix(E))
  colnames(Em) <- c("gene","sample","expr")
  Em$group <- meta[as.character(Em$sample), "group"]
  hp <- ggplot(Em, aes(x = sample, y = gene, fill = expr)) +
    geom_tile() + facet_grid(. ~ group, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Core IFN/ISG genes expression (human lung, GSE276060)", x = NULL, y = NULL)
  ggsave("results/GSE276060_IFN_genes_heatmap.png", hp, width = 12, height = 6, dpi = 150)
}

cat("\nDone. 输出见 results/GSE276060_*\n")
