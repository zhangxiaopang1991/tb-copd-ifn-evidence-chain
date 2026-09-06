suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(svglite)
  library(ragg)
})

OUT <- "figures"
dir.create(OUT, showWarnings = FALSE)

pal <- c(
  dark = "#252525",
  mid = "#6B7280",
  light = "#D9DEE5",
  blue = "#2F6F9F",
  teal = "#2A9D8F",
  red = "#C94C4C",
  orange = "#D98E2B"
)

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.5, colour = pal["mid"]),
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
}

save_pub <- function(plot, stem, width_mm = 183, height_mm = 120) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(file.path(OUT, paste0(stem, ".svg")), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(file.path(OUT, paste0(stem, ".pdf")),
                       width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  grDevices::tiff(file.path(OUT, paste0(stem, ".tiff")),
                  width = w, height = h, units = "in", res = 600,
                  compression = "lzw", type = "cairo")
  print(plot)
  dev.off()
  grDevices::png(file.path(OUT, paste0(stem, ".png")),
                 width = w, height = h, units = "in", res = 300,
                 type = "cairo")
  print(plot)
  dev.off()
}

write_source <- function(df, stem) {
  write.csv(df, file.path(OUT, paste0(stem, "_source.csv")), row.names = FALSE)
}

clean_pathway <- function(x) {
  x |>
    gsub("IFNa", "IFN-alpha", x = _) |>
    gsub("IFNg", "IFN-gamma", x = _) |>
    gsub("IFN_antigen_presentation", "Antigen presentation", x = _) |>
    gsub("_", " ", x = _)
}

# Figure 1: infection-stage enrichment across MTB, NTM and related contrasts.
fg1 <- read.csv("results/GSE276060_fgsea.csv", check.names = FALSE) |>
  filter(pathway %in% c("IFNa", "IFNg",
                        "IFN_antigen_presentation", "PTLD_IFN"),
         contrast %in% c("MTB_vs_Ctl", "NTM_vs_Ctl", "MTB_vs_NTM")) |>
  mutate(
    pathway_label = clean_pathway(pathway),
    contrast = recode(contrast,
                      MTB_vs_Ctl = "MTB vs lung control",
                      NTM_vs_Ctl = "NTM vs lung control",
                      MTB_vs_NTM = "MTB vs NTM"),
    pathway_label = factor(pathway_label,
                           levels = c("IFN-alpha", "IFN-gamma",
                                      "Antigen presentation", "PTLD IFN"))
  )
write_source(fg1, "Figure1_GSE276060")
p1 <- ggplot(fg1, aes(x = NES, y = pathway_label, colour = contrast, size = -log10(pmax(padj, 1e-12)))) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = pal["mid"]) +
  geom_point(alpha = 0.95) +
  scale_colour_manual(values = c("MTB vs lung control" = unname(pal["red"]),
                                 "NTM vs lung control" = unname(pal["blue"]),
                                 "MTB vs NTM" = unname(pal["teal"]))) +
  scale_size_continuous(range = c(1.8, 4.5), name = "-log10(FDR)") +
  labs(x = "Normalized enrichment score", y = NULL,
       title = "Infection-stage interferon activation",
       subtitle = "GSE276060; point size indicates adjusted significance") +
  theme_pub() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2), size = guide_legend(nrow = 2))
save_pub(p1, "Figure1_infection_stage", 183, 92)

# Figure 2: paired smoke-amplification evidence.
fg2 <- read.csv("results/EMTAB17246_fgsea_PTLD_IFN.csv", check.names = FALSE) |>
  filter(pathway %in% c("PTLD_IFN", "IFNa", "IFNg"),
         contrast %in% c("Mtb_Air", "Mtb_Smoke", "Mtb_int")) |>
  mutate(
    pathway_label = factor(clean_pathway(pathway),
                           levels = c("IFN-alpha", "IFN-gamma", "PTLD IFN")),
    contrast = recode(contrast,
                      Mtb_Air = "Mtb + air",
                      Mtb_Smoke = "Mtb + smoke",
                      Mtb_int = "Mtb x smoke interaction")
  )
write_source(fg2, "Figure2_EMTAB17246")
p2 <- ggplot(fg2, aes(x = contrast, y = NES, fill = contrast)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = pal["mid"]) +
  geom_col(width = 0.65, colour = NA) +
  geom_text(aes(label = paste0("FDR=", format.pval(padj, digits = 2, eps = 1e-3))),
            vjust = -0.35, size = 2.4, colour = pal["dark"]) +
  facet_wrap(~ pathway_label, nrow = 1) +
  scale_fill_manual(values = c("Mtb + air" = unname(pal["light"]),
                               "Mtb + smoke" = unname(pal["orange"]),
                               "Mtb x smoke interaction" = unname(pal["red"]))) +
  labs(x = NULL, y = "Normalized enrichment score",
       title = "Smoke exposure amplifies the Mtb-associated interferon response") +
  theme_pub() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
save_pub(p2, "Figure2_smoke_amplification", 183, 92)

# Figure 3: lesion-level directional paired estimates, without significance marks.
les <- read.csv("GSE192483/results/GSE192483_IFN_score_H_vs_L_pseudobulk_paired.csv",
                check.names = FALSE) |>
  filter(score %in% c("PTLD_IFN", "IFNa", "IFNg")) |>
  mutate(
    celltype = factor(celltype),
    module = recode(score, PTLD_IFN = "PTLD IFN",
                    IFNa = "IFN-alpha", IFNg = "IFN-gamma")
  )
les_top <- les |>
  group_by(celltype) |>
  summarise(max_abs = max(abs(mean_diff_HmL), na.rm = TRUE), .groups = "drop") |>
  slice_max(max_abs, n = 12) |>
  pull(celltype)
les <- les |>
  filter(celltype %in% les_top) |>
  mutate(celltype = reorder(celltype, mean_diff_HmL))
write_source(les, "Figure3_GSE192483")
p3 <- ggplot(les, aes(x = mean_diff_HmL, y = celltype, colour = module)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = pal["mid"]) +
  geom_point(size = 1.8, position = position_dodge(width = 0.55)) +
  facet_wrap(~ module, nrow = 1) +
  scale_colour_manual(values = c("PTLD IFN" = unname(pal["blue"]),
                                 "IFN-alpha" = unname(pal["teal"]),
                                 "IFN-gamma" = unname(pal["orange"]))) +
  labs(x = "Lesion minus non-lesion mean score",
       y = NULL,
       title = "Lesion-level single-cell evidence is directional",
       subtitle = "Donor-aware paired comparisons; no FDR-significant result") +
  theme_pub() +
  theme(legend.position = "none")
save_pub(p3, "Figure3_lesion_directional", 183, 120)

# Figure 4: chronic-stage covariate sensitivity and FEV1 association.
bulk4 <- read.csv("GSE47460/GSE47460_IFN_covariate_sensitivity_fgsea.csv",
                  check.names = FALSE) |>
  filter(pathway %in% c("PTLD_IFN", "IFNa")) |>
  mutate(
    model = recode(model, base_platform = "Platform",
                   plus_smoking = "+ smoking",
                   full_covariates = "Full covariates"),
    pathway = recode(pathway, PTLD_IFN = "PTLD IFN", IFNa = "IFN-alpha")
  )
bulk4$model <- factor(bulk4$model,
                      levels = c("Platform", "+ smoking", "Full covariates"))
write_source(bulk4, "Figure4_GSE47460_bulk")
p4a <- ggplot(bulk4, aes(x = model, y = NES, group = pathway, colour = pathway)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = pal["mid"]) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("PTLD IFN" = unname(pal["blue"]),
                                 "IFN-alpha" = unname(pal["teal"]))) +
  labs(x = NULL, y = "NES", title = "Disease contrast is covariate-sensitive") +
  theme_pub() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 25, hjust = 1))

fev4 <- read.csv("GSE47460/GSE47460_FEV1_IFN_core_genes.csv",
                 check.names = FALSE) |>
  filter(gene %in% c("IFIT1", "IFI44", "USP18", "BST2", "ISG15", "STAT1"),
         model %in% c("platform", "plus_smoking", "full_covariates")) |>
  mutate(model = recode(model, platform = "Platform",
                        plus_smoking = "+ smoking",
                        full_covariates = "Full covariates"),
         gene = reorder(gene, coef_fev1))
fev4$model <- factor(fev4$model,
                     levels = c("Platform", "+ smoking", "Full covariates"))
write_source(fev4, "Figure4_GSE47460_FEV1")
p4b <- ggplot(fev4, aes(x = coef_fev1, y = gene, colour = model)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = pal["mid"]) +
  geom_point(size = 1.8, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c("Platform" = unname(pal["blue"]),
                                 "+ smoking" = unname(pal["teal"]),
                                 "Full covariates" = unname(pal["red"]))) +
  labs(x = "FEV1 coefficient", y = NULL,
       title = "Core-gene associations persist") +
  theme_pub() +
  theme(legend.position = "bottom")
p4 <- p4a + p4b + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub(p4, "Figure4_chronic_stage", 183, 112)

# Figure 5: score, count-based state and composition layers.
score5 <- read.csv("GSE136831/results/GSE136831_pseudobulk_celltype_test.csv",
                   check.names = FALSE) |>
  filter(padj < 0.05) |>
  mutate(celltype = reorder(celltype, mean_diff))
write_source(score5, "Figure5_GSE136831_score")
p5a <- ggplot(score5, aes(x = mean_diff, y = celltype, colour = score)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = pal["mid"]) +
  geom_errorbar(aes(xmin = ci95_low, xmax = ci95_high),
                width = 0.25, orientation = "y") +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("PTLD_IFN" = unname(pal["blue"]),
                                 "IFNa" = unname(pal["teal"]),
                                 "IFNg" = unname(pal["orange"]))) +
  labs(x = "COPD minus control score", y = NULL, title = "Subject-level score") +
  theme_pub() +
  theme(legend.position = "none")

edge5 <- read.csv("GSE136831/results/GSE136831_edgeR_IFN_summary.csv",
                  check.names = FALSE) |>
  filter(global_padj < 0.05) |>
  mutate(label = paste(celltype, pathway, sep = " | "),
         label = reorder(label, NES))
write_source(edge5, "Figure5_GSE136831_edgeR")
p5b <- ggplot(edge5, aes(x = NES, y = label, colour = pathway)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = pal["mid"]) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("IFNa" = unname(pal["teal"]),
                                 "IFNg" = unname(pal["orange"]),
                                 "PTLD_IFN" = unname(pal["blue"]))) +
  labs(x = "NES", y = NULL, title = "Count-based differential state") +
  theme_pub() +
  theme(legend.position = "none")

comp5 <- read.csv("GSE136831/results/GSE136831_celltype_composition_test.csv",
                  check.names = FALSE) |>
  filter(padj < 0.05) |>
  mutate(celltype = reorder(celltype, mean_diff_COPD_minus_Control),
         mean_diff_pct = mean_diff_COPD_minus_Control * 100)
write_source(comp5, "Figure5_GSE136831_composition")
p5c <- ggplot(comp5, aes(x = mean_diff_pct, y = celltype)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = pal["mid"]) +
  geom_segment(aes(x = 0, xend = mean_diff_pct, yend = celltype),
               linewidth = 1.2, colour = pal["red"]) +
  geom_point(size = 1.8, colour = pal["red"]) +
  labs(x = "Difference in proportion (percentage points)", y = NULL,
       title = "Cell composition") +
  theme_pub()

p5 <- (p5a | p5b | p5c) + plot_annotation(
  title = "Chronic-stage single-cell results separate state from composition",
  subtitle = "Subject-level inference; Multiplet excluded from composition analysis",
  tag_levels = "a"
) & theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub(p5, "Figure5_single_cell_layers", 270, 125)

writeLines(c(
  "Figure contract",
  "Figure 1: infection-stage pathway activation in GSE276060.",
  "Figure 2: smoke amplification and interaction in E-MTAB-17246.",
  "Figure 3: directional lesion-level paired evidence; no significance marks.",
  "Figure 4: covariate-sensitive disease contrast and FEV1 associations.",
  "Figure 5: score, count-based and composition layers in GSE136831.",
  "All inferential units are sample, donor or subject; cells are not independent replicates.",
  "Exports: SVG, PDF, 600 dpi TIFF, 300 dpi PNG preview.",
  "Source data are saved as *_source.csv beside each figure."
), file.path(OUT, "figure_QA_notes.txt"))

cat("Generated Figure 1-5 and source-data tables in", OUT, "\n")
