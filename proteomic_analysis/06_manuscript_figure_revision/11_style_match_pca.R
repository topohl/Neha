# ==============================================================================
# VISUAL-FIDELITY PASS -- Figure 3 C and Supplementary B
#
# No statistics recomputed. PC scores, variance-explained values and metadata are
# read from the source-data CSVs already written by the first regeneration pass.
# Only the drawing changes.
#
# Fig 3C fixes vs the first pass:
#   legend moved from right to BELOW, horizontal, legend title dropped
#   point size / font sizes / grid colour taken from the historical SVG
#   panel aspect taken from the historical SVG (8 x 6 in per panel, rendered at 50%)
# Supp B fixes vs the first pass:
#   the three PCAs are stacked VERTICALLY in one file, as in the manuscript
#   "Plate" becomes "Collection plate"
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "10_style_match_common.R"))
suppressPackageStartupMessages(library(patchwork))

FIRST_F3  <- file.path(OUT_ROOT, "figure_panels", "Figure3")
FIRST_SUP <- file.path(OUT_ROOT, "figure_panels", "Supplementary_proteomics")

theme_pca_sm <- function() {
  theme_minimal(base_size = PCA_FS_AXIS, base_family = FONT) +
    theme(
      panel.grid.major = element_line(colour = PCA_GRID, linewidth = 0.4 * PT),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.title = element_text(colour = "#444444", size = PCA_FS_AXIS),
      axis.text  = element_text(colour = "#555555", size = PCA_FS_AXIS),
      axis.ticks = element_blank(),
      axis.line  = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(colour = "#222222", face = "bold", size = PCA_FS_TITLE,
                                hjust = 0.5, margin = margin(b = 4)),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_blank(),
      legend.text      = element_text(colour = "#444444", size = PCA_FS_LEGEND),
      legend.key       = element_blank(),
      legend.background = element_blank(),
      plot.margin = margin(4, 6, 4, 4)
    )
}

# ==============================================================================
# FIGURE 3 C
# ==============================================================================
sc <- utils::read.csv(file.path(FIRST_F3, "Fig3C_PCA_animal_level_source_data.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
vr <- utils::read.csv(file.path(FULL, "pca", "pca_variance_explained_animal_level.csv"),
                      stringsAsFactors = FALSE)
varp <- vr$percent_variance
cat(sprintf("PC1 %.1f%%  PC2 %.1f%%  PC3 %.1f%%  (n = %d units)\n", varp[1], varp[2], varp[3], nrow(sc)))

# manuscript legend order: neuropil, cFos, neuron, mCherry
LEG_ORDER <- c("neuropil", "cfos", "neuron", "mcherry")
sc$cls <- factor(sc$sample_class, levels = LEG_ORDER,
                 labels = unname(CLASS_LABELS[LEG_ORDER]))
pal <- setNames(unname(PCA_COLS[LEG_ORDER]), unname(CLASS_LABELS[LEG_ORDER]))

pca_sub <- function(xv, yv, xi, yi, ttl) {
  ggplot(sc, aes(.data[[xv]], .data[[yv]], colour = cls)) +
    geom_point(size = pt_radius_to_size(PCA_PT_R), shape = 16, alpha = PCA_ALPHA) +
    scale_colour_manual(values = pal, drop = FALSE) +
    theme_pca_sm() +
    labs(x = sprintf("%s (%.1f%%)", xv, varp[xi]),
         y = sprintf("%s (%.1f%%)", yv, varp[yi]),
         title = ttl)
}

p12 <- pca_sub("PC1", "PC2", 1, 2, "PC1 vs PC2")
p23 <- pca_sub("PC2", "PC3", 2, 3, "PC2 vs PC3")

fig3c <- p12 + p23 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

save_sm(fig3c, SM_F3, "Fig3C_PCA_animal_level_STYLE_MATCHED",
        width = 2 * PCA_PANEL_W, height = PCA_PANEL_H + 0.35)
copy_source_data(file.path(FIRST_F3, "Fig3C_PCA_animal_level_source_data.csv"),
                 SM_F3, "Fig3C_PCA_animal_level_STYLE_MATCHED_source_data.csv")

# ==============================================================================
# SUPPLEMENTARY B -- three PCAs stacked vertically
# ==============================================================================
# Supplementary B is reproduced at the full legacy canvas (7.5 x 6.2 in) with the
# legacy font sizes (axis 20 px, legend 14 px) rather than at half scale, so the
# axis tick labels sit exactly as they do in the historical base PCA plots.
theme_supb <- function() {
  theme_pca_sm() +
    theme(legend.position = "right", legend.direction = "vertical",
          axis.title = element_text(colour = "#444444", size = SUPB_FS_AX),
          axis.text  = element_text(colour = "#555555", size = SUPB_FS_AX),
          legend.text = element_text(colour = "#444444", size = SUPB_FS_LG),
          plot.title = element_text(hjust = 0, face = "bold", size = SUPB_FS_AX))
}

# B1 : animal level, Experimental Group
b1_src <- utils::read.csv(file.path(FIRST_SUP, "SuppB1_PCA_experimental_group_animal_level_source_data.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
b1_src$grp <- factor(b1_src$condition, levels = CONDITIONS, labels = unname(COND_LABELS[CONDITIONS]))
pal_grp <- setNames(c("#1E90FF", "#61D7FF", "#FFA500", "#FFD700"), unname(COND_LABELS[CONDITIONS]))

pB1 <- ggplot(b1_src, aes(PC1, PC2, colour = grp)) +
  geom_point(size = pt_radius_to_size(SUPB_PT_R), shape = 16, alpha = 0.8) +
  scale_colour_manual(values = pal_grp) +
  theme_supb() +
  labs(x = sprintf("PC1 (%.1f%%)", varp[1]), y = sprintf("PC2 (%.1f%%)", varp[2]),
       title = "Experimental Group")

# B2 / B3 : acquisition level, hemisphere and collection plate
b23 <- utils::read.csv(file.path(FIRST_SUP, "SuppB2_B3_PCA_acquisition_level_QC_source_data.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
qc_var <- utils::read.csv(file.path(FULL, "qc", "SuppB_pca_audit.csv"), stringsAsFactors = FALSE)
gv <- function(k) qc_var$value[qc_var$item == k]
v1 <- as.numeric(gv("acquisition_PC1_percent")); v2 <- as.numeric(gv("acquisition_PC2_percent"))

b23$hemi  <- factor(sub("Hemisphere ", "", b23$Hemisphere), levels = c("1", "2"))
b23$plate <- factor(sub("Plate", "", b23[["Collection plate"]]), levels = c("1", "2"))

pB2 <- ggplot(b23, aes(PC1, PC2, colour = hemi)) +
  geom_point(size = pt_radius_to_size(SUPB_PT_R), shape = 16, alpha = 0.8) +
  scale_colour_manual(values = SUPB_HEMI) +
  theme_supb() +
  labs(x = sprintf("PC1 (%.1f%%)", v1), y = sprintf("PC2 (%.1f%%)", v2), title = "Hemisphere")

pB3 <- ggplot(b23, aes(PC1, PC2, colour = plate)) +
  geom_point(size = pt_radius_to_size(SUPB_PT_R), shape = 16, alpha = 0.8) +
  scale_colour_manual(values = SUPB_PLATE) +
  theme_supb() +
  labs(x = sprintf("PC1 (%.1f%%)", v1), y = sprintf("PC2 (%.1f%%)", v2), title = "Collection plate")

suppB <- pB1 / pB2 / pB3 + plot_layout(ncol = 1)
save_sm(suppB, SM_SUP, "SuppB_PCA_stacked_STYLE_MATCHED",
        width = SUPB_W, height = 3 * SUPB_H)

# individual files too, so panels can be placed separately if preferred
save_sm(pB1, SM_SUP, "SuppB1_PCA_experimental_group_animal_level_STYLE_MATCHED", SUPB_W, SUPB_H)
save_sm(pB2, SM_SUP, "SuppB2_PCA_hemisphere_acquisition_level_QC_STYLE_MATCHED", SUPB_W, SUPB_H)
save_sm(pB3, SM_SUP, "SuppB3_PCA_collection_plate_acquisition_level_QC_STYLE_MATCHED", SUPB_W, SUPB_H)

copy_source_data(file.path(FIRST_SUP, "SuppB1_PCA_experimental_group_animal_level_source_data.csv"),
                 SM_SUP, "SuppB1_PCA_experimental_group_animal_level_STYLE_MATCHED_source_data.csv")
copy_source_data(file.path(FIRST_SUP, "SuppB2_B3_PCA_acquisition_level_QC_source_data.csv"),
                 SM_SUP, "SuppB2_B3_PCA_acquisition_level_QC_STYLE_MATCHED_source_data.csv")

message("PCA panels style-matched.")
