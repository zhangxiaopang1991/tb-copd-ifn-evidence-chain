suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(fgsea)
  library(ggplot2)
})

.script_args <- commandArgs(trailingOnly = FALSE)
.script_file <- sub("^--file=", "", .script_args[grep("^--file=", .script_args)][1])
if (!is.na(.script_file)) setwd(dirname(normalizePath(.script_file)))

# ---------- 复用分析.R 的前半段 ----------
counts <- read.delim("counts_matrix.tsv", row.names = 1, check.names = FALSE)
meta <- read.delim("sample_table.tsv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
one <- read.delim("sample001.count", header = FALSE, col.names = c("gid", "sym", "cnt"), stringsAsFactors = FALSE)
sym_map <- setNames(one$sym, one$gid)
gid_map <- setNames(one$gid, one$sym)

meta$infect_s <- unname(c("none" = "none", "Mycobacterium tuberculosis" = "Mtb", "Mycobacterium avium" = "Mav")[as.character(meta$infect)])
meta$stim_s   <- ifelse(grepl("smoke", meta$stimulus, ignore.case = TRUE), "Smoke", "Air")
meta$group <- factor(paste(meta$infect_s, meta$stim_s, sep = "_"),
                     levels = c("none_Air", "none_Smoke", "Mtb_Air", "Mtb_Smoke", "Mav_Air", "Mav_Smoke"))
meta$donor <- factor(meta$donor)
counts <- counts[, rownames(meta)]

dge <- DGEList(counts = counts)
keep <- filterByExpr(dge, group = meta$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~ donor + group, data = meta)
colnames(design) <- make.names(colnames(design))
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)

cont.matrix <- makeContrasts(
  Smoke_none      = groupnone_Smoke,
  Mtb_Air         = groupMtb_Air,
  Mtb_Smoke       = groupMtb_Smoke - groupnone_Smoke,
  Mtb_Smoke_vs_Air = groupMtb_Smoke - groupMtb_Air,
  Mtb_int         = (groupMtb_Smoke - groupMtb_Air) - groupnone_Smoke,
  Mav_Air         = groupMav_Air,
  Mav_Smoke       = groupMav_Smoke - groupnone_Smoke,
  Mav_int         = (groupMav_Smoke - groupMav_Air) - groupnone_Smoke,
  levels = design
)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# ---------- 自定义通路（symbol -> ENSG） ----------
custom_pathways <- list(
  IFN_antigen_presentation = c("STAT1","STAT2","IRF1","IRF7","IRF9","TAP1","TAP2","B2M",
                               "HLA-A","HLA-B","HLA-C","HLA-E","NLRC5","CIITA","TAPBP",
                               "GBP1","GBP2","GBP4","GBP5","MX1","MX2","OAS1","OAS2","OASL",
                               "ISG15","ISG20","IFIT1","IFIT2","IFIT3","IFI44","IFI6","IFI27",
                               "IFIH1","DDX58","PSMB8","PSMB9","PSME1","PSME2","JAK1","JAK2","TYK2"),
  Innate_sensing = c("TLR2","TLR4","NOD2","MYD88","TBK1","NFKB1","RELA","IRF3","IRF5",
                     "AIM2","NLRP3","CASP1"),
  Inflammation_cytokines = c("IL1B","IL6","CXCL8","IL18","IL33","TNF",
                           "CCL2","CCL5","CCL20","CXCL1","CXCL10","CXCL11","CSF2","CSF3"),
  Mucous_barrier = c("MUC5AC","MUC5B","MUC13","SPDEF","FOXJ1","SCGB1A1","DEFB4A","DEFB1"),
  ECM_remodeling_EMT = c("MMP1","MMP3","MMP7","MMP9","MMP12","SERPINE1","TGFB1","TGFB2",
                         "VIM","SNAI1","SNAI2","ZEB1","ZEB2","CDH1","CDH2","COL13A1","VCAN","TAGLN","CALD1","INHBA"),
  Oxidative_stress = c("G6PD","GCLC","GCLM","GSR","TKT","FTL","FTH1","HMOX1","NQO1","TXN","TXNRD1",
                       "SOD1","SOD2","CAT","GPX1","GPX2","PRDX1","NFE2L2","KEAP1","SLC7A11","GSTM1","ALDH3A1","CYP1B1","CYP4B1")
)

custom_index <- lapply(custom_pathways, function(syms) {
  g <- unique(gid_map[syms])
  g <- g[!is.na(g) & g %in% rownames(v)]
  g
})
custom_index <- custom_index[sapply(custom_index, length) >= 3]
cat("自定义通路基因集大小:\n"); print(sapply(custom_index, length))

# ---------- fgsea across all contrasts ----------
fgsea_all <- data.frame()
for (cn in colnames(cont.matrix)) {
  stats <- fit2$t[, cn]
  names(stats) <- rownames(fit2)
  stats <- stats[is.finite(stats)]
  res <- fgsea(pathways = custom_index, stats = stats, nperm = 10000, nproc = 1)
  res <- res[order(res$pval), ]
  res$contrast <- cn
  fgsea_all <- rbind(fgsea_all, res[, c("contrast","pathway","pval","padj","ES","NES","size","leadingEdge")])
}
fgsea_all$leadingEdge <- sapply(fgsea_all$leadingEdge, function(x) paste(x, collapse = ","))
write.csv(fgsea_all, "results/fgsea_custom_pathways.csv", row.names = FALSE)

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
       title = "Custom pathway enrichment (fgsea) across contrasts",
       subtitle = "FDR<0.25 highlighted; positive NES = up-regulated in numerator") +
  theme_bw(base_size = 9) +
  theme(legend.position = "top")

ggsave("results/fgsea_custom_pathways_dotplot.png", p, width = 13, height = 6, dpi = 150)

cat("\n== fgsea 自定义通路显著/趋势结果 (FDR<0.25) ==\n")
print(fgsea_all[fgsea_all$padj < 0.25, c("contrast","pathway","NES","pval","padj","size")])

cat("\nDone. Output: results/fgsea_custom_pathways.csv, results/fgsea_custom_pathways_dotplot.png\n")
