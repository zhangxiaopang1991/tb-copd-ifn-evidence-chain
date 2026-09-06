suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(pheatmap)
})
.script_args <- commandArgs(trailingOnly = FALSE)
.script_file <- sub("^--file=", "", .script_args[grep("^--file=", .script_args)][1])
if (!file.exists("GSE47460_gene_matrix.tsv") && !is.na(.script_file)) {
  setwd(dirname(.script_file))
}

mat <- read.delim("GSE47460_gene_matrix.tsv", row.names=1, check.names=FALSE)
meta <- read.delim("sample_meta_final.tsv", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta$gsm; meta <- meta[colnames(mat),]
keep <- meta$`disease state` %in% c("Control","Chronic Obstructive Lung Disease")
m2 <- meta[keep,]; mm <- mat[, rownames(m2)]
m2$disease <- factor(ifelse(m2$`disease state`=="Control","Control","COPD"), levels=c("Control","COPD"))
m2$gpl <- factor(m2$gpl)
mm <- as.matrix(mm); storage.mode(mm) <- "double"

# ---------- 图1：FEV1 关联森林图（核心 IFN 基因）----------
fev1 <- suppressWarnings(as.numeric(m2$`%predicted fev1 (pre-bd)`))
has <- !is.na(fev1) & m2$disease=="COPD"
mmc <- mm[, has]; fev <- fev1[has]; gplc <- m2$gpl[has]
design2 <- model.matrix(~ fev + gplc)
fitf <- lmFit(mmc, design2); fitf <- eBayes(fitf)
core <- c("IFI6","OAS1","OAS2","ISG15","MX1","IFIT1","STAT1","IRF9","GBP5","OASL",
          "IFI44","IFI27","IFI35","RSAD2","USP18","BST2","IFITM1","IFITM3","EIF2AK2","DDX58")
core <- core[core %in% rownames(mmc)]
res <- topTable(fitf, coef="fev", number=Inf, sort.by="none")
df <- data.frame(gene=core, coef=res[core,"logFC"], p=res[core,"P.Value"], stringsAsFactors=FALSE)
# 用 SE 计算 CI：从 t 和 logFC 反推
df$se <- res[core,"logFC"] / res[core,"t"]
df$lo <- df$coef - 1.96*df$se
df$hi <- df$coef + 1.96*df$se
df$sig <- ifelse(df$p < 0.05, "FDR/P<0.05", "n.s.")
df$gene <- factor(df$gene, levels=df$gene[order(df$coef)])
p1 <- ggplot(df, aes(x=coef, y=gene, color=sig)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50") +
  geom_point(size=2) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.3) +
  scale_color_manual(values=c("FDR/P<0.05"="firebrick","n.s."="grey60")) +
  labs(x="coefficient of FEV1 (%predicted) -> gene expression", y=NULL,
       title="GSE47460 COPD: FEV1 association with IFN genes",
       subtitle="positive coef = better lung function -> higher IFN gene; n(COPD with FEV1)=218") +
  theme_bw(base_size=9) + theme(legend.position="top")
ggsave("GSE47460_FEV1_IFN_forest.png", p1, width=8, height=6, dpi=150)

# ---------- 图2：核心 IFN 基因热图（COPD vs Control）----------
core_hm <- c("IFI6","OAS1","OAS2","ISG15","MX1","MX2","IFIT1","STAT1","IRF7","IRF9",
             "GBP1","GBP5","OASL","IFI44","IFI27","RSAD2","USP18","BST2","IFITM1","IFITM3","EIF2AK2")
core_hm <- core_hm[core_hm %in% rownames(mm)]
# 每个基因按样本均值中心化（行 z-score）
hm <- mm[core_hm, ]
hm <- t(scale(t(hm)))
# 样本列注释
anno <- data.frame(Disease=m2$disease, Platform=m2$gpl, row.names=rownames(m2))
anno_col <- list(Disease=c(Control="#2E86AB", COPD="#D1495B"), Platform=c(GPL14550="#88a825", GPL6480="#7c4dff"))
# 只抽样展示（328 样本太多，按 disease 各取前 60 用于可视化）
set.seed(1)
idx_ctrl <- sample(which(m2$disease=="Control"), min(60, sum(m2$disease=="Control")))
idx_copd <- sample(which(m2$disease=="COPD"), min(60, sum(m2$disease=="COPD")))
idx <- c(idx_ctrl, idx_copd)
hm_sub <- hm[, idx]
anno_sub <- anno[idx, , drop=FALSE]
png("GSE47460_IFN_genes_heatmap.png", width=10, height=6, units="in", res=150)
pheatmap(hm_sub, annotation_col=anno_sub, annotation_colors=anno_col,
         cluster_cols=FALSE, cluster_rows=TRUE, show_colnames=FALSE,
         color=colorRampPalette(c("#2E86AB","white","#D1495B"))(100),
         main="GSE47460 lung tissue: IFN genes (row z-scored), COPD vs Control",
         fontsize_row=7, border_color=NA)
dev.off()

cat("Done. 两张图已保存\n")
