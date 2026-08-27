# ==============================================================================
# VISUAL-FIDELITY PASS -- Figure 3 D, volcano
#
# No statistics recomputed. Every value is read from the source-data CSV already
# written by scripts/02_panelD_volcano.R.
#
# Fixes vs the first pass:
#   palette      the first pass used #CA0020 / #0571B0 / #CCCCCC, taken from the
#                *script* version recovered from git. The historical SVG actually
#                on disk uses #F36E21 (up), #465B65 (down), #BABABA at 60% opacity
#                (n.s.), which is what the manuscript shows. Switched to those.
#   title        the in-plot title is dropped; the manuscript has none.
#   geometry     5 x 5 in canvas, point radius 6.92 pt, axis text 17 px,
#                gene labels 11.90 px, dash pattern 4,4 -- all from the SVG.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "10_style_match_common.R"))
suppressPackageStartupMessages(library(ggrepel))

FIRST_F3 <- file.path(OUT_ROOT, "figure_panels", "Figure3")
src <- utils::read.csv(file.path(FIRST_F3, "Fig3D_mCherry_paired_vs_unpaired_volcano_animal_level_source_data.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
cat("rows:", nrow(src), " FDR up:", sum(src$Significance == "up"),
    " down:", sum(src$Significance == "down"), "\n")
lab <- src[src$labelled_in_panel, , drop = FALSE]
cat("labels:", paste(lab$Gene_Name, collapse = ", "), "\n")

x_limit <- ceiling(max(abs(src$log2fc), 0, na.rm = TRUE)) + 0.5

p <- ggplot(src, aes(x = log2fc, y = -log10(padj))) +
  geom_point(aes(colour = Significance), shape = 16,
             size = pt_radius_to_size(VOLC_PT_R), alpha = VOLC_NS_A) +
  geom_hline(yintercept = -log10(FDR), linetype = "dashed", colour = "#999999",
             linewidth = 1.2 * PT) +
  ggrepel::geom_text_repel(
    data = lab, aes(label = Gene_Name),
    size = VOLC_FS_LAB * PT,            # 11.90 px -> mm, as in the historical SVG
    min.segment.length = 0, max.overlaps = Inf, box.padding = 0.5,
    family = FONT, colour = "black", segment.colour = "black", segment.size = 0.25
  ) +
  scale_colour_manual(values = c(up = VOLC_UP, down = VOLC_DOWN, `n.s.` = VOLC_NS),
                      breaks = c("up", "down", "n.s.")) +
  scale_x_continuous(limits = c(-x_limit, x_limit)) +
  labs(x = expression(log[2] ~ Fold ~ Change), y = expression(-log[10](italic(P)[adj]))) +
  theme_classic(base_size = VOLC_FS_AX, base_family = FONT) +
  theme(
    legend.position = "none",
    axis.text  = element_text(colour = "#2C2C2C", size = VOLC_FS_AX),
    axis.title = element_text(colour = "#2C2C2C", size = VOLC_FS_AX),
    axis.line  = element_line(colour = "#2C2C2C", linewidth = 1.2 * PT),
    axis.ticks = element_line(colour = "#2C2C2C", linewidth = 1.0 * PT),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

save_sm(p, SM_F3, "Fig3D_mCherry_paired_vs_unpaired_volcano_animal_level_STYLE_MATCHED",
        width = VOLC_IN, height = VOLC_IN)
copy_source_data(file.path(FIRST_F3, "Fig3D_mCherry_paired_vs_unpaired_volcano_animal_level_source_data.csv"),
                 SM_F3, "Fig3D_mCherry_paired_vs_unpaired_volcano_animal_level_STYLE_MATCHED_source_data.csv")

message("Volcano style-matched.")
