# ==============================================================================
# Figure 3 panel C  -- PCA (PC1 vs PC2, PC2 vs PC3) coloured by sample class
# Supplementary panel B1 -- PCA coloured by Experimental Group (condition)
#
# ORIGINAL: 03_qc_exploration/06_pcaPlot_Neha.r, multipanel_pca_combinations()
#           -> 99_historical/pca_plots_legacy/plots/main_figure/pca_combinations_celltype.svg
#           hemisphere-level, 120 acquisitions, PC1 25.4% / PC2 17.7% / PC3 10.2%
# CORRECTED: same prcomp settings (centre + scale, zero-variance rows removed)
#           on the validated animal-level GCT, 48 AnimalID x sample_class units.
#
# The legacy file rendered three combinations side by side (PC1v2, PC2v3, PC1v3);
# the manuscript panel used the first two, so the two-panel layout is reproduced
# here. Point size, palette, theme and axis-label format are the legacy ones.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))
suppressPackageStartupMessages({ library(patchwork); library(digest) })

source(file.path(REPO_ROOT, "R", "analysis_labels.R"))
source(file.path(REPO_ROOT, "R", "protigy_input_utils.R"))
source(file.path(REPO_ROOT, "R", "neha_path_utils.R"))
source(file.path(REPO_ROOT, "R", "pca_animal_level_utils.R"))

STAGE <- ensure_dir(file.path(FULL, "pca"))

# ---- rebuild the canonical animal-level PCA ---------------------------------
gct_sha <- digest::digest(file = GCT_ANIMAL, algo = "sha256")
message("animal-level GCT sha256: ", gct_sha)

parsed    <- validate_protigy_gct_v13(GCT_ANIMAL)
validated <- validate_neha_pca_animal_input(parsed, expected_n = 3L)
prepared  <- prepare_neha_animal_pca(validated$expression_matrix, center = TRUE, scale. = TRUE)

pca  <- prepared$pca
meta <- validated$sample_metadata
stopifnot(identical(rownames(meta), rownames(pca$x)))

varp <- (pca$sdev^2) / sum(pca$sdev^2) * 100
message(sprintf("PC1 %.1f%%  PC2 %.1f%%  PC3 %.1f%%  (n = %d units)", varp[1], varp[2], varp[3], nrow(pca$x)))

# ---- legacy PCA theme (06a_pca_core.r theme_pca_min) ------------------------
theme_pca_min <- function() {
  theme_minimal(base_size = 16, base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#ECECEC", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.title = element_text(color = "#444444", size = 18),
      axis.text  = element_text(color = "#555555", size = 18),
      axis.ticks = element_line(color = "#DDDDDD", linewidth = 0.3),
      axis.line  = element_blank(),
      panel.border = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = "right",
      legend.title = element_text(color = "#444444", size = 12),
      legend.text  = element_text(color = "#666666", size = 13),
      plot.title   = element_text(color = "#222222", face = "bold", size = 16, margin = margin(b = 6)),
      panel.spacing = unit(0.6, "lines"),
      plot.margin = margin(8, 8, 6, 8),
      complete = TRUE
    )
}

scores <- data.frame(
  Sample = rownames(pca$x),
  meta[, c("AnimalID", "condition_code", "condition", "sample_class"), drop = FALSE],
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  stringsAsFactors = FALSE, row.names = NULL
)

# Publication sample-class labels; legacy legend showed the raw 'bg' alias.
scores$sample_class_label <- factor(unname(CLASS_LABELS[scores$sample_class]),
                                    levels = unname(CLASS_LABELS[SAMPLE_CLASSES]))
scores$experimental_group <- factor(scores$condition, levels = CONDITIONS,
                                    labels = unname(COND_LABELS[CONDITIONS]))

pca_scatter <- function(df, xv, yv, xi, yi, title, colour_var, palette, legend_title) {
  ggplot(df, aes(.data[[xv]], .data[[yv]], color = .data[[colour_var]])) +
    geom_point(size = 4, alpha = 0.8) +
    scale_color_manual(values = palette, name = legend_title) +
    theme_pca_min() +
    labs(x = sprintf("%s (%.1f%%)", xv, varp[xi]),
         y = sprintf("%s (%.1f%%)", yv, varp[yi]),
         title = title)
}

# ---- Figure 3 C --------------------------------------------------------------
pal_class <- setNames(PCA_PALETTE, unname(CLASS_LABELS[SAMPLE_CLASSES]))

p12 <- pca_scatter(scores, "PC1", "PC2", 1, 2, "PC1 vs PC2", "sample_class_label", pal_class, "Sample class")
p23 <- pca_scatter(scores, "PC2", "PC3", 2, 3, "PC2 vs PC3", "sample_class_label", pal_class, "Sample class")

fig3c <- p12 + p23 + patchwork::plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

save_panel(fig3c, STAGE, PANELS_F3, "Fig3C_PCA_animal_level.svg", width = 12, height = 6)
save_source_data(
  scores[, c("Sample", "AnimalID", "condition_code", "condition", "sample_class", "PC1", "PC2", "PC3")],
  STAGE, PANELS_F3, "Fig3C_PCA_animal_level_source_data.csv"
)

# also keep the full three-combination version produced by the legacy function,
# so the panel can be swapped for the original 3-up layout if preferred
p13 <- pca_scatter(scores, "PC1", "PC3", 1, 3, "PC1 vs PC3", "sample_class_label", pal_class, "Sample class")
combo3 <- p12 + p23 + p13 + patchwork::plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "right")
ggsave(file.path(STAGE, "pca_combinations_sample_class_animal_level_3up.svg"),
       combo3, width = 18, height = 6, device = svglite::svglite)

# ---- Supplementary B1 : Experimental Group ----------------------------------
pal_group <- setNames(PCA_PALETTE, unname(COND_LABELS[CONDITIONS]))
supB1 <- pca_scatter(scores, "PC1", "PC2", 1, 2, "Experimental Group",
                     "experimental_group", pal_group, "Experimental Group")
save_panel(supB1, STAGE, PANELS_SUP, "SuppB1_PCA_experimental_group_animal_level.svg",
           width = 7.5, height = 6)
save_source_data(
  scores[, c("Sample", "AnimalID", "condition_code", "condition", "sample_class", "PC1", "PC2", "PC3")],
  STAGE, PANELS_SUP, "SuppB1_PCA_experimental_group_animal_level_source_data.csv"
)

# ---- variance table ----------------------------------------------------------
var_tbl <- data.frame(PC = seq_along(varp), percent_variance = varp,
                      cumulative_percent = cumsum(varp))
utils::write.csv(var_tbl, file.path(STAGE, "pca_variance_explained_animal_level.csv"), row.names = FALSE)

audit <- data.frame(
  item = c("input_gct", "input_gct_sha256", "n_units", "n_proteins_after_zero_variance_removal",
           "center", "scale", "PC1_percent", "PC2_percent", "PC3_percent",
           "legacy_PC1_percent", "legacy_PC2_percent", "legacy_PC3_percent",
           "legacy_output", "legacy_sampling_unit"),
  value = c(GCT_ANIMAL, gct_sha, nrow(pca$x), nrow(prepared$matrix),
            TRUE, TRUE, sprintf("%.1f", varp[1]), sprintf("%.1f", varp[2]), sprintf("%.1f", varp[3]),
            "25.4", "17.7", "10.2",
            "99_historical/pca_plots_legacy/plots/main_figure/pca_combinations_celltype.svg",
            "hemisphere/acquisition (120 observations)"),
  stringsAsFactors = FALSE
)
utils::write.csv(audit, file.path(STAGE, "Fig3C_SuppB1_pca_audit.csv"), row.names = FALSE)

message("PCA panels done.")
