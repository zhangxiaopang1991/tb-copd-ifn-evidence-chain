suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(fgsea)
  library(ggplot2)
})
.script_args <- commandArgs(trailingOnly = FALSE)
.script_file <- sub("^--file=", "", .script_args[grep("^--file=", .script_args)][1])
if (!is.na(.script_file)) setwd(dirname(normalizePath(.script_file)))

# ---------- 复用 analyze.R 的配对设计 ----------
counts <- read.delim("counts_matrix.tsv", row.names = 1, check.names = FALSE)
meta <- read.delim("sample_table.tsv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
one <- read.delim("sample001.count", header = FALSE, col.names = c("gid","sym","cnt"), stringsAsFactors = FALSE)
sym_map <- setNames(one$sym, one$gid)   # ENSG -> symbol
gid_map <- setNames(one$gid, one$sym)   # symbol -> ENSG

meta$infect_s <- unname(c("none"="none","Mycobacterium tuberculosis"="Mtb","Mycobacterium avium"="Mav")[as.character(meta$infect)])
meta$stim_s   <- ifelse(grepl("smoke", meta$stimulus, ignore.case=TRUE), "Smoke", "Air")
meta$group <- factor(paste(meta$infect_s, meta$stim_s, sep="_"),
                     levels=c("none_Air","none_Smoke","Mtb_Air","Mtb_Smoke","Mav_Air","Mav_Smoke"))
meta$donor <- factor(meta$donor)
counts <- counts[, rownames(meta)]

dge <- DGEList(counts=counts)
keep <- filterByExpr(dge, group=meta$group)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~ donor + group, data=meta)
colnames(design) <- make.names(colnames(design))
v <- voom(dge, design, plot=FALSE)
fit <- lmFit(v, design)
cont.matrix <- makeContrasts(
  Mtb_Air = groupMtb_Air,
  Mtb_Smoke = groupMtb_Smoke - groupnone_Smoke,
  Mtb_int = (groupMtb_Smoke - groupMtb_Air) - groupnone_Smoke,
  Mav_Air = groupMav_Air,
  Mav_Smoke = groupMav_Smoke - groupnone_Smoke,
  Mav_int = (groupMav_Smoke - groupMav_Air) - groupnone_Smoke,
  levels = design
)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# ---------- 新基因集（symbol -> ENSG）----------
gs_sym <- list(
  PTLD_obstructive_IFN = c("IFI6","OAS1","OAS2","ISG15","MX1","MX2","IFIT1","IFIT2","IFIT3",
                           "STAT1","STAT2","IRF7","IRF9","GBP1","GBP2","GBP5","OASL","IFI44",
                           "IFI27","IFI35","IFITM1","IFITM3","RSAD2","USP18","XAF1","BST2",
                           "EPSTI1","DDX58","IFIH1","EIF2AK2"),
  HALLMARK_IFNa = c("ISG15","MX1","OAS1","OAS2","IFIT1","IFIT2","IFIT3","IFI6","IFI27","IFI44",
                    "IFI44L","IFITM1","IFITM2","IFITM3","RSAD2","USP18","XAF1","BST2","EPSTI1",
                    "EIF2AK2","GBP1","GBP2","DDX58","IFIH1","IRF7","IRF9","STAT1","STAT2",
                    "SAMD9","SAMD9L","PARP9","DTX3L","TRIM22","TRIM25","TNFSF10","ISG20","OASL"),
  HALLMARK_IFNg = c("STAT1","STAT2","IRF1","IRF2","IRF7","IRF8","IRF9","GBP1","GBP2","GBP4","GBP5",
                    "CXCL9","CXCL10","CXCL11","HLA-A","HLA-B","HLA-C","TAP1","TAP2","B2M",
                    "PSMB8","PSMB9","PSME1","PSME2","NLRC5","CIITA","IFITM1","IFITM3","OAS1",
                    "OAS2","OASL","MX1","MX2","ISG15","ISG20","IFI35","SOCS1","SOCS3",
                    "JAK1","JAK2","CCL2","CCL5","CCL7")
)
gs <- lapply(gs_sym, function(syms){ g <- unique(gid_map[syms]); g[!is.na(g) & g %in% rownames(v)] })
gs <- gs[sapply(gs, length) >= 3]
cat("基因集大小（交集后）:\n"); print(sapply(gs, length))

# ---------- fgsea ----------
fgsea_all <- data.frame()
for (cn in colnames(cont.matrix)) {
  stats <- fit2$t[, cn]; names(stats) <- rownames(fit2); stats <- stats[is.finite(stats)]
  res <- fgsea(pathways=gs, stats=stats, nperm=20000, nproc=1)
  res <- res[order(res$pval),]; res$contrast <- cn
  fgsea_all <- rbind(fgsea_all, res[, c("contrast","pathway","pval","padj","ES","NES","size","leadingEdge")])
}
fgsea_all$leadingEdge <- sapply(fgsea_all$leadingEdge, function(x) paste(x, collapse=","))
write.csv(fgsea_all, "results/EMTAB17246_fgsea_PTLD_IFN.csv", row.names=FALSE)

plot_df <- fgsea_all
plot_df$signed_logFDR <- sign(plot_df$ES) * (-log10(pmax(plot_df$padj, 1e-10)))
plot_df$trend <- plot_df$padj < 0.25
p <- ggplot(plot_df, aes(x=signed_logFDR, y=pathway)) +
  geom_point(aes(size=size, alpha=trend), color="steelblue") +
  geom_vline(xintercept=0, linetype="dashed", color="grey40") +
  facet_wrap(~contrast, ncol=3) +
  scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.35)) +
  labs(x="signed -log10(FDR)", y=NULL, title="E-MTAB-17246 (epithelium) re-run with PTLD IFN gene sets",
       subtitle="positive NES = up in numerator; FDR<0.25 highlighted") +
  theme_bw(base_size=9) + theme(legend.position="top")
ggsave("results/EMTAB17246_fgsea_PTLD_IFN_dotplot.png", p, width=11, height=4.5, dpi=150)

cat("\n== fgsea 结果 (FDR<0.25) ==\n")
print(fgsea_all[fgsea_all$padj < 0.25, c("contrast","pathway","NES","pval","padj","size")])
cat("\nDone.\n")
