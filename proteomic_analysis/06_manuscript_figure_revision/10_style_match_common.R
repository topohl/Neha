# ==============================================================================
# VISUAL-FIDELITY PASS -- shared constants recovered from the historical outputs
#
# Every value below was read out of a real historical SVG on disk (path given),
# except where marked RECONSTRUCTED, which means it was read off the manuscript
# raster because the historical file differs or does not exist.
#
# NOTHING HERE TOUCHES A STATISTIC.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))

HIST <- file.path(DATA_ROOT, "99_historical")
SM        <- file.path(OUT_ROOT, "figure_panels_style_matched")
SM_F3     <- ensure_dir(file.path(SM, "Figure3"))
SM_SUP    <- ensure_dir(file.path(SM, "Supplementary_proteomics"))
SM_PREV   <- ensure_dir(file.path(SM, "previews"))
SM_MAN    <- ensure_dir(file.path(SM, "manifests"))

PT   <- 0.3527778          # pt -> mm
FONT <- "Arial"
# ggplot geom_point size -> rendered radius in pt is size * .pt / 2 = size * 1.4225
pt_radius_to_size <- function(r_pt) r_pt / 1.4225

save_sm <- function(plot, dir, stem, width, height) {
  svg <- file.path(dir, paste0(stem, ".svg"))
  png <- file.path(SM_PREV, paste0(stem, ".png"))
  ggplot2::ggsave(svg, plot, width = width, height = height, units = "in", device = svglite::svglite)
  ggplot2::ggsave(png, plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
  message("  saved: ", basename(svg), " + preview")
  invisible(c(svg = svg, png = png))
}

copy_source_data <- function(from, dir, to) {
  if (file.exists(from)) {
    file.copy(from, file.path(dir, to), overwrite = TRUE)
  } else {
    warning("source data not found: ", from)
  }
  invisible(file.path(dir, to))
}

# ---- Figure 3 C : 99_historical/pca_plots_legacy/plots/main_figure/pca_combinations_celltype.svg
# canvas 1728 x 432 pt for 3 panels -> 576 x 432 pt (8 x 6 in) per panel, fonts at 2x.
# Reproduced here at 50% nominal scale: 4 x 3 in per panel, fonts halved.
PCA_COLS      <- c(neuropil = "#1E90FF", cfos = "#61D7FF", mcherry = "#FFA500", neuron = "#FFD700")
PCA_PT_R      <- 12.80 / 2          # pt, halved with the canvas
PCA_ALPHA     <- 0.80
PCA_FS_AXIS   <- 36 / 2
PCA_FS_TITLE  <- 32 / 2
PCA_FS_LEGEND <- 28 / 2
PCA_GRID      <- "#DDDDDD"
PCA_PANEL_W   <- 8 / 2
PCA_PANEL_H   <- 6 / 2

# ---- Figure 3 D : .../memory_ensemble/05_Volcano_Plots/Volcano_mcherry2_mcherry4.svg
# canvas 360 x 360 pt (5 x 5 in)
VOLC_UP     <- "#F36E21"
VOLC_DOWN   <- "#465B65"
VOLC_NS     <- "#BABABA"
VOLC_NS_A   <- 0.60
VOLC_PT_R   <- 6.92                 # pt
VOLC_FS_AX  <- 17
VOLC_FS_LAB <- 11.90
VOLC_DASH   <- c(4, 4)
VOLC_IN     <- 5

# ---- Figure 3 E / Supp D : 03_output/qc/Rank-Abundance_Marker_QC.svg
RANK_BG     <- "#D9D9D9"
RANK_COLS   <- c(`Mitochondrial` = "#D35400", `Growth-Related Plasticity` = "#C71585",
                 `Neuron Soma` = "#7F8C8D", `Neuropil/Structure` = "#2980B9")
# category renaming is RECONSTRUCTED from the manuscript legend text; the hex values
# are the ones read out of Rank-Abundance_Marker_QC.svg
RANK_SETS <- list(
  `Mitochondrial`             = c("Mrpl53", "Timm9", "Ndufb7", "Yars2", "Ndufs6", "Mrpl50", "Mrpl4"),
  `Growth-Related Plasticity` = c("Aktip", "Naa10", "Mycbp", "Ikbkb", "Acvr1b", "Strap"),
  `Neuron Soma`               = c("Rps9", "Rpl22", "Rps16", "Brd4", "Rpl35", "H1-1", "H1-2", "H1-3", "Rps18", "Rpl6"),
  `Neuropil/Structure`        = c("Vcan", "Kcna1", "Mog", "Cntnap1", "Cnp")
)

# ---- Figure 3 F / Supp F : .../memory_ensemble/02_Main_Plots/Dotplot_Enrichment_TopGenes_PerComp.svg
# canvas 403.20 x 576.00 pt (5.6 x 8 in). NES ramp identified numerically from the
# embedded legend PNG: colorRampPalette(c("#4C87C6","white","#FAA51A")) fits with
# mean|RGB diff| = 0.0102 versus 0.0798 for the blue-white-red ramp.
DOT_NES_LOW  <- "#4C87C6"
DOT_NES_MID  <- "white"
DOT_NES_HIGH <- "#FAA51A"
DOT_FS       <- 10
DOT_BAR_W    <- 11.52
DOT_BAR_H    <- 57.60
DOT_SIZE_RNG <- c(1.5, 5)
DOT_W        <- 5.6
DOT_H        <- 8

# ---- Figure 3 H : 99_historical/compare_sig_expr/plots/scatter_learning_logfc_vs_cno_logfc_mcherry.svg
# canvas 288 x 432 pt (4 x 6 in)
SCAT_PT      <- "#2C3E50"
SCAT_PT_A    <- 0.30
SCAT_PT_R    <- 4.98
SCAT_FIT     <- "#E74C3C"
SCAT_FIT_LW  <- 4.27
SCAT_ZERO    <- "#CCCCCC"
SCAT_FS_AX   <- 20
SCAT_FS_TIT  <- 22
SCAT_W       <- 4
SCAT_H       <- 4.6   # legacy file is 4 x 6 in; the manuscript panel is nearer square, so the
                      # height is reduced while every recovered mark size is kept unchanged

# ---- Figure 3 I : 99_historical/compare_sig_expr/plots/NES_Heatmap_Sorted_By_Comparison_Group.svg
# canvas 864 x 720 pt (12 x 10 in); pheatmap of t(nes_matrix): terms on x rotated -45,
# two value rows, an annotation row on top coloured by Comparison_Group.
HEAT_LOW     <- "#2d8be9"
HEAT_MID     <- "#F7F7F7"
HEAT_HIGH    <- "#d12f42"
HEAT_ANNOT   <- c(neuropil = "#4c87c6", cfos = "#6ccff6", mcherry = "#faa51a", neuron = "#fdd700")
HEAT_CELL    <- 15          # pt
HEAT_ANN_H   <- 10          # pt
HEAT_FS_TERM <- 8
HEAT_FS_ROW  <- 10
HEAT_BORDER  <- "white"

# ---- Supplementary B : 99_historical/pca_plots_legacy/plots/base/*.svg
# canvas 540 x 446.4 pt (7.5 x 6.2 in), circles r 15.29 pt, fill-opacity 0.80,
# axis text/title 20 px, legend text 14 px. Reproduced at full legacy size.
SUPB_PT_R  <- 15.29
SUPB_W     <- 7.5
SUPB_H     <- 6.2
SUPB_FS_AX <- 20
SUPB_FS_LG <- 14
# Hemisphere uses the historical #1E90FF / #61D7FF. The manuscript recoloured the
# plate panel to red / khaki; those two hex values are RECONSTRUCTED from the
# manuscript raster because the historical SVG still uses the blue pair.
SUPB_HEMI  <- c(`1` = "#1E90FF", `2` = "#61D7FF")
SUPB_PLATE <- c(`1` = "#C0392B", `2` = "#C3BC9F")
