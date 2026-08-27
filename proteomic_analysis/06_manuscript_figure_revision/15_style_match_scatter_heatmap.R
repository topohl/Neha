# ==============================================================================
# VISUAL-FIDELITY PASS -- Figure 3 H (scatter) and Figure 3 I (NES heatmap)
#
# No statistics recomputed. Both panels are drawn from the source-data CSVs
# written by scripts/06_panelH_scatter_panelI_heatmap.R.
#
# Fixes vs the first pass:
#
# 3H  * bold title and "DESCRIPTIVE. Pearson correlation ..." subtitle removed;
#       the manuscript shows a plain "r = ..." inside the panel at top left
#     * point colour #2C3E50 at 30% opacity, radius 4.98 pt, fit line #E74C3C at
#       4.27 pt, zero lines #CCCCCC dashed, fonts 20/22 px, canvas 4 x 6 in --
#       all read from
#       99_historical/compare_sig_expr/plots/scatter_learning_logfc_vs_cno_logfc_mcherry.svg
#
# 3I  * LAYOUT CORRECTED. The first pass reproduced
#       NES_Absolute_Heatmap_ggplot.svg (terms on y, two dataset columns,
#       faceted by sample class). The manuscript panel is the OTHER historical
#       artefact from the same script: NES_Heatmap_Sorted_By_Comparison_Group.svg,
#       a pheatmap of t(nes_matrix) -- terms along x rotated -45, two value rows,
#       and a sample-class annotation strip on top. Rebuilt in that layout.
#     * palette colorRampPalette(c("#2d8be9","#F7F7F7","#d12f42")), annotation
#       colours c(neuropil #4c87c6, cfos #6ccff6, mcherry #faa51a, neuron #fdd700),
#       white borders, 15 pt cells, 10 pt annotation strip, term labels 8 px,
#       row labels 10 px -- all read from that SVG.
#     * no in-plot title; the manuscript sets the panel title in Illustrator.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "10_style_match_common.R"))
suppressPackageStartupMessages(library(patchwork))

FIRST_F3 <- file.path(OUT_ROOT, "figure_panels", "Figure3")

# ==============================================================================
# FIGURE 3 H
# ==============================================================================
h <- utils::read.csv(file.path(FIRST_F3, "Fig3H_logFC_correlation_animal_level_source_data.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
sig <- h[h$learning_signature %in% c(TRUE, "TRUE"), , drop = FALSE]
r_sig <- stats::cor(sig$log2fc_learning, sig$log2fc_cno, use = "complete.obs")
cat(sprintf("Fig3H  n = %d   r = %.2f  (value carried over, not recomputed differently)\n",
            nrow(sig), r_sig))

ph <- ggplot(sig, aes(x = log2fc_learning, y = log2fc_cno)) +
  geom_hline(yintercept = 0, colour = SCAT_ZERO, linetype = "dashed", linewidth = 1.0 * PT) +
  geom_vline(xintercept = 0, colour = SCAT_ZERO, linetype = "dashed", linewidth = 1.0 * PT) +
  geom_point(colour = SCAT_PT, alpha = SCAT_PT_A, size = pt_radius_to_size(SCAT_PT_R), shape = 16) +
  geom_smooth(method = "lm", formula = y ~ x, colour = SCAT_FIT, fill = SCAT_FIT,
              alpha = 0.2, linewidth = SCAT_FIT_LW * PT) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.6,
           label = sprintf("r = %.2f", r_sig), family = FONT, size = SCAT_FS_TIT * PT,
           colour = "#2C2C2C") +
  labs(x = "Learning logFC\n(paired vs unpaired)", y = "CNO vs VEH logFC (paired)") +
  theme_minimal(base_size = SCAT_FS_AX, base_family = FONT) +
  theme(
    text = element_text(family = FONT, colour = "#2C2C2C"),
    axis.text  = element_text(size = SCAT_FS_AX, colour = "#2C2C2C"),
    axis.title = element_text(size = SCAT_FS_TIT, colour = "#2C2C2C"),
    panel.grid = element_blank(),
    axis.line  = element_line(colour = "#2C2C2C", linewidth = 0.8 * PT),
    axis.ticks = element_line(colour = "#2C2C2C", linewidth = 0.8 * PT),
    plot.title = element_blank(), plot.subtitle = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 10, 6, 6)
  )

save_sm(ph, SM_F3, "Fig3H_logFC_correlation_animal_level_STYLE_MATCHED",
        width = SCAT_W, height = SCAT_H)
copy_source_data(file.path(FIRST_F3, "Fig3H_logFC_correlation_animal_level_source_data.csv"),
                 SM_F3, "Fig3H_logFC_correlation_animal_level_STYLE_MATCHED_source_data.csv")

# ==============================================================================
# FIGURE 3 I -- transposed pheatmap layout with a sample-class annotation strip
# ==============================================================================
long <- utils::read.csv(file.path(FIRST_F3, "Fig3I_CNO_pathway_heatmap_animal_level_source_data.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
# legacy ordering: arrange(Comparison_Group, Description)
cls_levels <- unname(CLASS_LABELS[c("neuropil", "cfos", "mcherry", "neuron")])
long$Comparison_Group <- factor(long$Comparison_Group, levels = cls_levels)
ord <- unique(long[order(long$Comparison_Group, long$Description), c("Description", "Comparison_Group", "sample_class")])
long$Description <- factor(long$Description, levels = ord$Description)
long$Row <- factor(ifelse(long$Dataset == "Memory Ensemble", "Memory Engram", "Inhibition"),
                   levels = c("Inhibition", "Memory Engram"))   # top row drawn last
lim <- max(abs(long$NES), na.rm = TRUE)
cat(sprintf("Fig3I  terms = %d  rows = 2  NES limits +/- %.2f\n", nlevels(long$Description), lim))

ann <- ord
ann$Description <- factor(ann$Description, levels = levels(long$Description))
ann$strip <- "Sample class"

heat_theme <- function() {
  theme_minimal(base_size = HEAT_FS_ROW, base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = "#2C2C2C"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_blank(), legend.background = element_blank(),
      legend.title = element_text(size = HEAT_FS_ROW, face = "bold"),
      legend.text = element_text(size = HEAT_FS_ROW),
      legend.position = "right", legend.justification = "top",
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.margin = margin(0, 4, 0, 132)
    )
}

p_ann <- ggplot(ann, aes(x = Description, y = strip, fill = Comparison_Group)) +
  geom_tile(colour = HEAT_BORDER, linewidth = 0.5 * PT * 2) +
  scale_fill_manual(values = setNames(unname(HEAT_ANNOT[c("neuropil", "cfos", "mcherry", "neuron")]),
                                      cls_levels),
                    name = "Sample class", drop = TRUE) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  heat_theme() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

p_heat <- ggplot(long, aes(x = Description, y = Row, fill = NES)) +
  geom_tile(colour = HEAT_BORDER, linewidth = 0.5 * PT * 2) +
  scale_fill_gradientn(colours = grDevices::colorRampPalette(c(HEAT_LOW, HEAT_MID, HEAT_HIGH),
                                                             space = "Lab")(100),
                       limits = c(-lim, lim), name = "NES",
                       guide = guide_colourbar(barwidth = unit(3, "mm"), barheight = unit(20, "mm"))) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0), position = "right") +
  heat_theme() +
  theme(axis.text.x = element_text(size = HEAT_FS_TERM, angle = 45, hjust = 1, vjust = 1,
                                   colour = "#2C2C2C"),
        axis.text.y = element_text(size = HEAT_FS_ROW, hjust = 0, colour = "#2C2C2C"))

n_terms <- nlevels(long$Description)
pI <- p_ann / p_heat +
  plot_layout(heights = c(HEAT_ANN_H, 2 * HEAT_CELL), guides = "collect") &
  theme(legend.position = "right")

save_sm(pI, SM_F3, "Fig3I_CNO_pathway_heatmap_animal_level_STYLE_MATCHED",
        width = max(9, 2.6 + n_terms * 0.155 + 2.6), height = 4.6)
copy_source_data(file.path(FIRST_F3, "Fig3I_CNO_pathway_heatmap_animal_level_source_data.csv"),
                 SM_F3, "Fig3I_CNO_pathway_heatmap_animal_level_STYLE_MATCHED_source_data.csv")

message("Scatter and NES heatmap style-matched.")
