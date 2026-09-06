suppressPackageStartupMessages({ library(edgeR); library(Matrix); library(fgsea) })

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args)) args[1] else "."
setwd(base)
pb_dir <- file.path(base, "pseudobulk_counts")

genes <- read.delim(file.path(pb_dir, "genes.tsv"), stringsAsFactors = FALSE)
gene_names <- make.unique(genes$symbol)

gene_sets <- list(
  PTLD_IFN = c("IFI6","OAS1","OAS2","ISG15","MX1","MX2","IFIT1","IFIT2","IFIT3","STAT1","STAT2","IRF7","IRF9","GBP1","GBP2","GBP5","OASL","IFI44","IFI27","IFI35","IFITM1","IFITM3","RSAD2","USP18","XAF1","BST2","EPSTI1","DDX58","IFIH1","EIF2AK2"),
  IFNa = c("ISG15","MX1","OAS1","OAS2","IFIT1","IFIT2","IFIT3","IFI6","IFI27","IFI44","IFI44L","IFITM1","IFITM2","IFITM3","RSAD2","USP18","XAF1","BST2","EPSTI1","EIF2AK2","GBP1","GBP2","DDX58","IFIH1","IRF7","IRF9","STAT1","STAT2","SAMD9","SAMD9L","PARP9","DTX3L","TRIM22","TRIM25","TNFSF10","ISG20","OASL"),
  IFNg = c("STAT1","STAT2","IRF1","IRF7","IRF9","GBP1","GBP2","GBP5","CXCL9","CXCL10","CXCL11","HLA-A","HLA-B","HLA-C","TAP1","TAP2","B2M","PSMB8","PSMB9","PSME1","PSME2","NLRC5","CIITA","IFITM1","IFITM3","OAS1","OAS2","OASL","MX1","MX2","ISG15","ISG20","IFI35","SOCS1","SOCS3","JAK1","JAK2","CCL2","CCL5","CCL7")
)

files <- list.files(pb_dir, pattern = "^counts_.*\\.mtx\\.gz$", full.names = TRUE)
pathway_out <- list(); de_out <- list(); meta_out <- list()

for (f in files) {
  safe <- sub("^counts_", "", basename(f))
  safe <- sub("\\.mtx\\.gz$", "", safe)
  if (safe == "Multiplet") next
  mfile <- file.path(pb_dir, paste0("metadata_", safe, ".tsv"))
  md <- read.delim(mfile, stringsAsFactors = FALSE, colClasses = "character")
  md$n_cells <- as.numeric(md$n_cells)
  counts <- readMM(gzfile(f))
  counts <- as.matrix(counts)
  rownames(counts) <- gene_names
  colnames(counts) <- md$sample_id
  md <- md[md$n_cells >= 20, , drop = FALSE]
  counts <- counts[, md$sample_id, drop = FALSE]
  md$disease <- factor(md$Disease_Identity, levels = c("Control", "COPD"))
  if (sum(md$disease == "COPD") < 3 || sum(md$disease == "Control") < 3) next
  keep <- filterByExpr(counts, group = md$disease)
  counts <- counts[keep, , drop = FALSE]
  y <- DGEList(counts = counts, group = md$disease)
  y <- calcNormFactors(y)
  design <- model.matrix(~ disease, data = md)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "diseaseCOPD")
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  tt$gene <- rownames(tt)
  tt$celltype <- unique(md$Manuscript_Identity)
  tt$sample_COPD <- sum(md$disease == "COPD")
  tt$sample_Control <- sum(md$disease == "Control")
  de_out[[length(de_out) + 1L]] <- tt[, c("celltype","gene","logFC","logCPM","F","PValue","FDR","sample_COPD","sample_Control")]

  stat <- tt$logFC
  names(stat) <- tt$gene
  stat <- stat[is.finite(stat)]
  gs <- lapply(gene_sets, function(x) intersect(x, names(stat)))
  gs <- gs[lengths(gs) >= 3]
  if (length(gs)) {
    fg <- fgsea(pathways = gs, stats = stat, nproc = 1)
    fg$celltype <- unique(md$Manuscript_Identity)
    fg$sample_COPD <- sum(md$disease == "COPD")
    fg$sample_Control <- sum(md$disease == "Control")
    fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ",", FUN.VALUE = character(1))
    pathway_out[[length(pathway_out) + 1L]] <- fg[, c("celltype","pathway","sample_COPD","sample_Control","pval","padj","ES","NES","size","leadingEdge")]
  }
  meta_out[[length(meta_out) + 1L]] <- data.frame(celltype = unique(md$Manuscript_Identity),
    n_pseudobulk = nrow(md), n_COPD = sum(md$disease == "COPD"), n_Control = sum(md$disease == "Control"),
    median_cells = median(md$n_cells), min_cells = min(md$n_cells))
  cat(unique(md$Manuscript_Identity), ": ", nrow(md), " samples; genes after filter=", nrow(counts), "\n", sep = "")
}

pathway <- do.call(rbind, pathway_out)
pathway$global_padj <- p.adjust(pathway$pval, method = "BH")
de <- do.call(rbind, de_out)
meta_summary <- do.call(rbind, meta_out)
write.csv(pathway, file.path(base, "results", "GSE136831_edgeR_IFN_summary.csv"), row.names = FALSE)
write.csv(de, file.path(base, "results", "GSE136831_edgeR_DEG_all_celltypes.csv"), row.names = FALSE)
write.csv(meta_summary, file.path(base, "results", "GSE136831_edgeR_pseudobulk_metadata.csv"), row.names = FALSE)
cat("Wrote edgeR pseudobulk outputs.\n")
