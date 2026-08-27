# ==============================================================================
# VISUAL-FIDELITY PASS -- Figure 3 F and Supplementary F, GSEA bubble grids
#
# No statistics recomputed. NES, p.adjust, term membership and column membership
# are read verbatim from the source-data CSVs written by
# scripts/04_panelF_suppF_gsea_dotplots.R. Only row/column ORDER (a layout
# operation, the legacy order_dotplot() algorithm) is re-derived here.
#
# Fixes vs the first pass:
#   NES palette   the first pass used colorRampPalette(c("#0571B0","white",
#                 "#CA0020")) taken from the git version of 02_compareGO.r.
#                 The legend colourbar embedded in the historical SVG decodes to
#                 colorRampPalette(c("#4C87C6","white","#FAA51A"))
#                 (mean|RGB diff| 0.0102 vs 0.0798 for blue-white-red), which is
#                 the blue-white-ORANGE ramp the manuscript shows. Switched.
#   title         "Top terms per comparison" removed; the manuscript has no title.
#   geometry      canvas 5.6 x 8 in, text 10 px, colourbar 11.52 x 57.60 pt,
#                 size range 1.5-5 -- all read from the historical SVG.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "10_style_match_common.R"))
suppressPackageStartupMessages(library(stringr))

FIRST_F3  <- file.path(OUT_ROOT, "figure_panels", "Figure3")
FIRST_SUP <- file.path(OUT_ROOT, "figure_panels", "Supplementary_proteomics")

# legacy order_dotplot(): hclust the NES matrix for columns, then rows by
# dominant column / direction / |NES|.  Layout only.
order_dotplot <- function(df, fixed_col_order = NULL) {
  agg <- stats::aggregate(NES ~ Description + Comparison, data = df[!is.na(df$NES), ], FUN = mean)
  mat <- tapply(agg$NES, list(agg$Description, agg$Comparison), function(z) z[1])
  mat[is.na(mat)] <- 0
  if (ncol(mat) < 2 || nrow(mat) < 2) return(list(row_order = rownames(mat), col_order = colnames(mat)))
  # Column order: the manuscript shows the fixed biological order
  # neuropil > cFos > mCherry > neuron, so it is supplied rather than clustered.
  # (The legacy hclust call is kept below for reference.)
  col_order <- if (!is.null(fixed_col_order)) {
    fixed_col_order[fixed_col_order %in% colnames(mat)]
  } else {
    colnames(mat)[stats::hclust(stats::dist(t(mat)))$order]
  }
  mat <- mat[, col_order, drop = FALSE]
  dom_comp <- apply(mat, 1, function(x) colnames(mat)[which.max(abs(x))])
  dom_val  <- apply(mat, 1, function(x) x[which.max(abs(x))])
  od <- data.frame(Term = rownames(mat), Comp = factor(dom_comp, levels = col_order),
                   Val = dom_val, stringsAsFactors = FALSE)
  od$Comp_Index <- as.integer(od$Comp)
  od$Direction  <- ifelse(od$Val > 0, 1, 0)
  od$Abs_Val    <- abs(od$Val)
  od <- od[order(-od$Comp_Index, od$Direction, od$Abs_Val), , drop = FALSE]
  list(row_order = od$Term, col_order = col_order)
}

dot_panel <- function(csv, col_levels) {
  d <- utils::read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
  d <- d[d$Comparison %in% col_levels, , drop = FALSE]
  ord <- order_dotplot(d, fixed_col_order = col_levels)
  d$Comparison  <- factor(as.character(d$Comparison), levels = ord$col_order)
  d$Description <- factor(d$Description, levels = ord$row_order)
  lim <- max(abs(d$NES), na.rm = TRUE)

  p <- ggplot(d, aes(x = Comparison, y = Description, colour = NES, size = -log10(p.adjust))) +
    geom_point(alpha = 0.85, stroke = 0.3) +
    scale_colour_gradientn(
      colours = grDevices::colorRampPalette(c(DOT_NES_LOW, DOT_NES_MID, DOT_NES_HIGH),
                                            space = "Lab")(100),
      name = "NES", limits = c(-lim, lim),
      guide = guide_colourbar(barwidth = unit(DOT_BAR_W * PT, "mm"),
                              barheight = unit(DOT_BAR_H * PT, "mm"), order = 1)
    ) +
    scale_size_continuous(name = expression(-log[10](italic(P)[adj])),
                          range = DOT_SIZE_RNG, limits = c(0, NA)) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)), drop = FALSE) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = DOT_FS, base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = "#2C2C2C"),
      axis.text.x = element_text(size = DOT_FS, colour = "#2C2C2C", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = DOT_FS, colour = "#2C2C2C"),
      axis.line = element_line(colour = "#2C2C2C", linewidth = 0.6 * PT),
      axis.ticks = element_line(colour = "#2C2C2C", linewidth = 0.5 * PT),
      axis.ticks.length = unit(3 * PT, "mm"),
      legend.title = element_text(size = DOT_FS, colour = "#2C2C2C"),
      legend.text = element_text(size = DOT_FS * 0.9, colour = "#2C2C2C"),
      legend.position = "right", legend.justification = "top",
      legend.key = element_blank(), legend.background = element_blank(),
      panel.grid = element_blank(), panel.border = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
  list(plot = p, n_rows = nlevels(d$Description), n_cols = nlevels(d$Comparison))
}

# ---- Figure 3 F --------------------------------------------------------------
f3f <- dot_panel(file.path(FIRST_F3, "Fig3F_memory_engram_GSEA_animal_level_source_data.csv"),
                 unname(CLASS_LABELS[c("neuropil", "cfos", "mcherry", "neuron")]))
cat("Fig3F rows:", f3f$n_rows, " cols:", f3f$n_cols, "\n")
save_sm(f3f$plot, SM_F3, "Fig3F_memory_engram_GSEA_animal_level_STYLE_MATCHED",
        width = DOT_W, height = max(5, 2.5 + f3f$n_rows * 0.25))
copy_source_data(file.path(FIRST_F3, "Fig3F_memory_engram_GSEA_animal_level_source_data.csv"),
                 SM_F3, "Fig3F_memory_engram_GSEA_animal_level_STYLE_MATCHED_source_data.csv")

# ---- Supplementary F ---------------------------------------------------------
sf <- dot_panel(file.path(FIRST_SUP, "SuppF_inhibition_GSEA_animal_level_source_data.csv"),
                unname(CLASS_LABELS[c("neuropil", "cfos", "mcherry")]))
cat("SuppF rows:", sf$n_rows, " cols:", sf$n_cols, "\n")
save_sm(sf$plot, SM_SUP, "SuppF_inhibition_GSEA_animal_level_STYLE_MATCHED",
        width = DOT_W - 0.6, height = max(5, 2.5 + sf$n_rows * 0.25))
copy_source_data(file.path(FIRST_SUP, "SuppF_inhibition_GSEA_animal_level_source_data.csv"),
                 SM_SUP, "SuppF_inhibition_GSEA_animal_level_STYLE_MATCHED_source_data.csv")

message("GSEA dotplots style-matched.")
