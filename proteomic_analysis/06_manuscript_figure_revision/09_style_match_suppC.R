# ==============================================================================
# VISUAL-FIDELITY PASS -- Supplementary panel C
# Cell-type identity heatmap (EWCE Baseline, animal level)
#
# NO STATISTICS ARE RECOMPUTED HERE. The values are read verbatim from the
# canonical animal-level EWCE table and are identical to those already written
# by scripts/08_suppC_ewce_suppE_audit.R. Only the plotting style changes.
#
# ------------------------------------------------------------------------------
# AESTHETICS RECOVERED DIRECTLY FROM THE HISTORICAL SVG
#   99_historical/ewce_legacy/Plots/Indiv_Heatmap_Baseline.svg
#
#   canvas            340.16 x 340.16 pt  (4.7244 in square)
#   panel region      x 41.22 -> 301.02, y 14.39 -> 301.31
#   tile              16.04 x 20.21 pt, stroke #FFFFFF width 0.43 pt
#   palette           viridisLite::magma  (identified numerically from the
#                     embedded legend colourbar PNG: mean|RGB diff| = 0.013
#                     vs magma, 0.066 inferno, 0.114 plasma, 0.250 viridis)
#   scale limits      range of the data, NOT centred on zero.
#                     Legend geometry gives [-2.212, 14.891]; the legacy table
#                     range is [-2.181, 14.863].
#   legend            colourbar 8.50 x 42.52 pt at the right, white ticks
#                     (stroke 0.38 pt), title "Z-score" Arial 6 pt BOLD,
#                     tick labels Arial 6 pt, breaks every 5
#   axis text         Arial 6.00 px, y right-anchored, x rotated -45 deg
#                     anchored at the end
#   axis lines        black, stroke-width 0.75 pt, left and bottom only
#   axis ticks        black, 0.75 pt, length 1.74 pt
#   grid              none; panel background white; no panel border
#
# RECONSTRUCTED FROM THE MANUSCRIPT IMAGE (not from the historical SVG)
#   row set           13 rows. The historical SVG has 14 and includes Int12
#                     between Int15 and CA2Pyr2; the manuscript panel does not
#                     show Int12, so it is excluded here.
#   column labels     neuropil_1..4 (the historical SVG says bg_1..bg_4;
#                     the manuscript uses the publication name).
#   title             the historical SVG carries "Significant Subtypes Landscape
#                     (Baseline)" (Arial 8 pt bold). The manuscript panel has no
#                     title, so none is drawn.
#
# STATISTIC VERIFIED EQUIVALENT
#   Inverting the historical tile fills through magma over the recovered legend
#   range reproduces the legacy table's `sd_from_mean` with r = 0.99960 and
#   max abs difference 0.46 (300-step colourbar quantisation). The panel plots
#   EWCE sd_from_mean (Z score), AnalysisType == "Baseline", annotation level 2.
#   The corrected animal-level panel plots the same field, same analysis type,
#   same annotation level, same 16 sample_class x condition targets.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))
suppressPackageStartupMessages({ library(readxl); library(viridisLite) })
source(file.path(REPO_ROOT, "R", "analysis_labels.R"))

SM       <- file.path(OUT_ROOT, "figure_panels_style_matched")
SM_SUP   <- ensure_dir(file.path(SM, "Supplementary_proteomics"))
SM_PREV  <- ensure_dir(file.path(SM, "previews"))
SM_MAN   <- ensure_dir(file.path(SM, "manifests"))

PT <- 0.3527778   # pt -> mm, for ggplot linewidth / unit conversions
FONT <- "Arial"

# ---- recovered constants -----------------------------------------------------
CANVAS_IN      <- 340.16 / 72          # 4.7244 in, square
TILE_STROKE_PT <- 0.43
AXIS_LINE_PT   <- 0.75
TICK_LEN_PT    <- 1.74
BAR_W_PT       <- 8.50
BAR_H_PT       <- 42.52
LEG_TICK_PT    <- 0.38
FS             <- 6                    # Arial 6.00 px

ROWS13 <- c("S1PyrL6", "S1PyrL23", "S1PyrDL", "Oligo6", "Oligo5", "Oligo4",
            "Oligo3", "Oligo2", "Oligo1", "Int15", "CA2Pyr2", "CA1Pyr2", "CA1Pyr1")
COND_CODE <- c(paired_cno = 1L, paired_veh = 2L, unpaired_cno = 3L, unpaired_veh = 4L)
CLASS_ORDER <- c("neuropil", "cfos", "mcherry", "neuron")
COLS16 <- as.vector(t(outer(CLASS_ORDER, 1:4, paste, sep = "_")))

# ---- data: read verbatim, no recomputation ----------------------------------
tbl <- file.path(EWCE_ROOT, "02_Tables_Supplements", "Supplementary_Table_EWCE.xlsx")
pr  <- as.data.frame(readxl::read_excel(tbl, sheet = "Primary_Results"))
b   <- pr[pr$AnalysisType == "Baseline", , drop = FALSE]
stopifnot(length(unique(b$TopN)) == 1L, length(unique(b$AnnotLevel)) == 1L)

h <- b[b$CellType %in% ROWS13, , drop = FALSE]
h$sample_class <- normalize_sample_class(h$Stratum)
h$condition    <- normalize_condition(h$Metric)
h$Column <- paste(h$sample_class, COND_CODE[h$condition], sep = "_")

stopifnot(setequal(h$CellType, ROWS13), setequal(h$Column, COLS16), nrow(h) == 13L * 16L)

h$CellType <- factor(h$CellType, levels = rev(ROWS13))   # rev: first label at top
h$Column   <- factor(h$Column, levels = COLS16)

lims <- range(h$sd_from_mean)
cat(sprintf("corrected Z range: [%.3f, %.3f]   (legacy legend range [-2.212, 14.891])\n",
            lims[1], lims[2]))

# ---- plot --------------------------------------------------------------------
p <- ggplot(h, aes(x = Column, y = CellType, fill = sd_from_mean)) +
  geom_tile(colour = "#FFFFFF", linewidth = TILE_STROKE_PT * PT) +
  scale_fill_gradientn(
    colours = viridisLite::magma(256),
    limits  = lims,
    breaks  = seq(0, floor(lims[2] / 5) * 5, by = 5),
    name    = "Z-score",
    guide   = guide_colourbar(
      barwidth  = unit(BAR_W_PT * PT, "mm"),
      barheight = unit(BAR_H_PT * PT, "mm"),
      ticks.colour = "#FFFFFF",
      ticks.linewidth = LEG_TICK_PT * PT,
      frame.colour = NA
    )
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  # recovered tile geometry: 16.04 pt wide x 20.21 pt tall
  coord_fixed(ratio = 20.21 / 16.04) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = FS, base_family = FONT) +
  theme(
    text             = element_text(family = FONT, colour = "black"),
    axis.text.x      = element_text(size = FS, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y      = element_text(size = FS, colour = "black", hjust = 1),
    axis.line        = element_line(colour = "black", linewidth = AXIS_LINE_PT * PT),
    axis.ticks       = element_line(colour = "black", linewidth = AXIS_LINE_PT * PT),
    axis.ticks.length = unit(TICK_LEN_PT * PT, "mm"),
    legend.title     = element_text(size = FS, face = "bold", family = FONT, colour = "black"),
    legend.text      = element_text(size = FS, family = FONT, colour = "black"),
    legend.position  = "right",
    legend.justification = "center",
    legend.key       = element_blank(),
    legend.background = element_blank(),
    legend.margin    = margin(0, 0, 0, 2, "pt"),
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    plot.title       = element_blank(),
    plot.subtitle    = element_blank(),
    plot.margin      = margin(4, 2, 2, 2, "pt")
  )

svg_path <- file.path(SM_SUP, "SuppC_celltype_identity_heatmap_animal_level_STYLE_MATCHED.svg")
png_path <- file.path(SM_PREV, "SuppC_celltype_identity_heatmap_animal_level_STYLE_MATCHED.png")
ggsave(svg_path, p, width = CANVAS_IN, height = CANVAS_IN, units = "in", device = svglite::svglite)
ggsave(png_path, p, width = CANVAS_IN, height = CANVAS_IN, units = "in", dpi = 400)
message("saved: ", svg_path)
message("saved: ", png_path)

# ---- source data: same values, laid out as plotted --------------------------
src <- h[order(h$Column, match(as.character(h$CellType), ROWS13)),
         c("CellType", "Column", "sample_class", "condition", "Target", "TopN", "AnnotLevel",
           "Direction", "fold_change", "sd_from_mean", "p", "q_target", "q_global", "Significant_Global")]
names(src)[names(src) == "sd_from_mean"] <- "Z_score_sd_from_mean"
utils::write.csv(src, file.path(SM_SUP, "SuppC_celltype_identity_heatmap_animal_level_STYLE_MATCHED_source_data.csv"),
                 row.names = FALSE)

# ---- aesthetic provenance ledger --------------------------------------------
aes_tbl <- data.frame(
  property = c("canvas_size", "aspect", "palette", "scale_type", "scale_limits",
               "tile_stroke_colour", "tile_stroke_width_pt", "axis_font", "axis_font_size_px",
               "x_label_rotation", "x_label_anchor", "axis_line_colour", "axis_line_width_pt",
               "tick_length_pt", "legend_position", "legend_bar_width_pt", "legend_bar_height_pt",
               "legend_tick_colour", "legend_title", "legend_title_face", "legend_breaks",
               "panel_grid", "panel_border", "plot_title", "significance_markers",
               "row_order", "column_labels", "row_set"),
  value = c(sprintf("%.2f x %.2f pt", 340.16, 340.16), "square (1:1)", "viridisLite::magma",
            "sequential continuous, limits = data range (not zero-centred)",
            sprintf("[%.3f, %.3f]", lims[1], lims[2]),
            "#FFFFFF", TILE_STROKE_PT, "Arial", FS, "-45 deg", "end / hjust = 1",
            "black", AXIS_LINE_PT, TICK_LEN_PT, "right", BAR_W_PT, BAR_H_PT, "#FFFFFF",
            "Z-score", "bold", "every 5", "none", "none", "none (omitted)", "none",
            paste(ROWS13, collapse = " > "), paste(COLS16, collapse = ", "),
            "13 cell types (manuscript set)"),
  recovered_from = c(rep("historical SVG", 25),
                     "manuscript image", "manuscript image (SVG says bg_1..bg_4)",
                     "manuscript image (SVG has 14 rows incl. Int12)"),
  stringsAsFactors = FALSE
)
utils::write.csv(aes_tbl, file.path(SM_MAN, "SuppC_recovered_aesthetics.csv"), row.names = FALSE)

message("Supplementary C style-matched. Rows: ", nlevels(h$CellType), "  Columns: ", nlevels(h$Column))
