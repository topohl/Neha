# ==============================================================================
# Figure 3 panel F      -- GSEA bubble grid, paired_veh vs unpaired_veh
#                          columns neuropil / cFos / mCherry / neuron
# Supplementary panel F -- GSEA bubble grid, paired_cno vs paired_veh
#                          columns neuropil / cFos / mCherry
#
# ORIGINAL: 04_differential_expression_enrichment/02_compareGO.r (pre-2026
#           3140-line version), dotplot block at legacy line 913-969
#           Fig3F  -> 99_historical/compareGO/BP/learning_signature/
#                     memory_ensemble/02_Main_Plots/Dotplot_Enrichment_TopGenes_PerComp.svg
#           SuppF  -> .../effects_inhibition_memory_ensemble/02_Main_Plots/
#                     Dotplot_Enrichment_TopGenes_PerComp.svg
#           input: hemisphere-level GSEA ranked by log2fc
#
# CORRECTED: identical term-selection rule, ordering algorithm, palette, size
#           scale, theme and dimension formula; input replaced by the canonical
#           animal-level GSEA ranked by the MODERATED t STATISTIC
#           (03_output/enrichment/enrichment_t_rank_validation_20260825/).
#
# Selection rule preserved verbatim from the original:
#   significant_only = TRUE (p.adjust < 0.05)
#   best instance per Description per comparison by |NES|
#   per comparison: top 5 by NES among positives + top 5 by NES among negatives
#   union of those Descriptions plotted for EVERY column, significant or not
#   row/column order from order_dotplot(): hclust on the NES matrix for columns,
#   then rows by dominant column, direction and |NES|
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))
suppressPackageStartupMessages({ library(stringr); library(tibble) })

STAGE <- ensure_dir(file.path(FULL, "enrichment"))

# ---- legacy constants --------------------------------------------------------
custom_palette   <- colorRampPalette(c("#0571B0", "white", "#CA0020"), space = "Lab")
target_n_terms   <- 5
min_set_size     <- 10
significant_only <- TRUE

theme_clean <- function(base_size = 9, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "#2C2C2C", family = base_family, lineheight = 1.3),
      axis.line.x = element_line(color = "#2C2C2C", linewidth = 0.6),
      axis.line.y = element_line(color = "#2C2C2C", linewidth = 0.6),
      axis.ticks = element_line(color = "#2C2C2C", linewidth = 0.5),
      axis.ticks.length = unit(3, "pt"),
      axis.text = element_text(color = "#2C2C2C", size = rel(0.95)),
      axis.title = element_text(color = "#2C2C2C", size = rel(1.0)),
      legend.background = element_blank(), legend.box.background = element_blank(),
      legend.key = element_rect(color = NA, fill = NA),
      legend.key.size = unit(10, "pt"), legend.key.height = unit(10, "pt"),
      legend.title = element_text(color = "#2C2C2C", size = rel(0.95)),
      legend.text = element_text(color = "#2C2C2C", size = rel(0.9)),
      legend.position = "right", legend.justification = "top",
      plot.title = element_text(color = "#2C2C2C", size = rel(1.15), face = "bold",
                                hjust = 0, vjust = 1, margin = margin(b = 8)),
      plot.subtitle = element_text(color = "#555555", size = rel(0.95), hjust = 0, margin = margin(b = 5)),
      strip.text = element_text(color = "#2C2C2C", size = rel(0.95)),
      strip.background = element_rect(color = "#EEEEEE", fill = "#EEEEEE", linewidth = 0.4),
      panel.grid = element_blank(), panel.border = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

calc_dims <- function(df_plot) {
  n_cols <- length(unique(as.character(df_plot$Comparison)))
  n_rows <- length(unique(as.character(df_plot$Description)))
  list(w = max(3.35, 3.35 + (n_cols * 0.5)), h = max(5, 2.5 + (n_rows * 0.25)))
}

order_dotplot <- function(df) {
  agg <- stats::aggregate(NES ~ Description + Comparison, data = df[!is.na(df$NES), ], FUN = mean)
  mat <- tapply(agg$NES, list(agg$Description, agg$Comparison), function(z) z[1])
  mat[is.na(mat)] <- 0
  if (ncol(mat) < 2 || nrow(mat) < 2) return(list(row_order = rownames(mat), col_order = colnames(mat)))
  hc_col <- stats::hclust(stats::dist(t(mat)))
  col_order <- colnames(mat)[hc_col$order]
  dom_comp <- apply(mat, 1, function(x) colnames(mat)[which.max(abs(x))])
  dom_val  <- apply(mat, 1, function(x) x[which.max(abs(x))])
  od <- data.frame(Term = rownames(mat), Comp = factor(dom_comp, levels = col_order), Val = dom_val,
                   stringsAsFactors = FALSE)
  od$Comp_Index <- as.integer(od$Comp)
  od$Direction  <- ifelse(od$Val > 0, 1, 0)
  od$Abs_Val    <- abs(od$Val)
  od <- od[order(-od$Comp_Index, od$Direction, od$Abs_Val), , drop = FALSE]
  list(row_order = od$Term, col_order = col_order)
}

# best instance per Description per Comparison, by |NES| (legacy combined_df_best)
best_per_term <- function(df) {
  df <- df[order(df$Comparison, df$Description, -abs(df$NES)), , drop = FALSE]
  df[!duplicated(df[, c("Comparison", "Description")]), , drop = FALSE]
}

select_standard <- function(best) {
  pool <- if (significant_only) best[is.finite(best$p.adjust) & best$p.adjust < FDR, , drop = FALSE] else best
  if (!nrow(pool)) return(pool[0, , drop = FALSE])
  do.call(rbind, lapply(unique(pool$Comparison), function(cp) {
    d  <- pool[pool$Comparison == cp, , drop = FALSE]
    up <- d[d$NES > 0, , drop = FALSE]; up <- head(up[order(-up$NES), , drop = FALSE], target_n_terms)
    dn <- d[d$NES < 0, , drop = FALSE]; dn <- head(dn[order(dn$NES), , drop = FALSE], target_n_terms)
    rbind(up, dn)
  }))
}

build_dotplot <- function(classes, numerator, denominator, plot_title, stub, panel_dir, panel_file) {
  raw <- do.call(rbind, lapply(classes, function(cl) {
    g <- read_gsea(cl, numerator, denominator)
    g$Comparison <- unname(CLASS_LABELS[cl])
    g[, c("Comparison", "sample_class", "canonical_contrast", "ID", "Description",
          "setSize", "NES", "pvalue", "p.adjust", "qvalue", "core_enrichment")]
  }))
  raw$Comparison <- factor(raw$Comparison, levels = unname(CLASS_LABELS[classes]))

  best <- best_per_term(raw)
  sel  <- select_standard(best)

  per_class_fdr <- do.call(rbind, lapply(classes, function(cl) data.frame(
    sample_class = cl,
    column = unname(CLASS_LABELS[cl]),
    n_terms_tested = sum(raw$sample_class == cl),
    n_FDR_significant = sum(raw$sample_class == cl & raw$p.adjust < FDR, na.rm = TRUE),
    n_terms_contributed_to_panel = if (nrow(sel)) sum(sel$sample_class == cl) else 0L,
    stringsAsFactors = FALSE)))
  print(per_class_fdr)

  if (!nrow(sel)) {
    warning("No FDR-significant terms in any column for ", plot_title, "; no panel written.")
    return(list(plot = NULL, data = raw, selected = sel, per_class = per_class_fdr))
  }

  plot_df <- best[best$Description %in% unique(sel$Description), , drop = FALSE]
  ord <- order_dotplot(plot_df)
  plot_df$Comparison  <- factor(as.character(plot_df$Comparison), levels = ord$col_order)
  plot_df$Description <- factor(plot_df$Description, levels = ord$row_order)

  lim <- max(abs(plot_df$NES), na.rm = TRUE)
  p <- ggplot(plot_df, aes(x = Comparison, y = Description, color = NES, size = -log10(p.adjust))) +
    geom_point(alpha = 0.85, stroke = 0.3) +
    scale_color_gradientn(colours = custom_palette(100), name = "NES", limits = c(-lim, lim)) +
    scale_size_continuous(name = expression(-log[10](italic(P)[adj])), range = c(1.5, 5), limits = c(0, NA)) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)), drop = FALSE) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
    labs(title = plot_title, x = NULL, y = NULL) +
    theme_clean(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major.y = element_blank(), legend.position = "right")

  dims <- calc_dims(plot_df)
  save_panel(p, STAGE, panel_dir, panel_file, width = dims$w, height = dims$h)

  src <- plot_df[, c("Comparison", "sample_class", "canonical_contrast", "ID", "Description",
                     "setSize", "NES", "pvalue", "p.adjust", "qvalue")]
  src$fdr_significant <- is.finite(src$p.adjust) & src$p.adjust < FDR
  src$selected_by_column <- src$Description %in% sel$Description
  save_source_data(src, STAGE, panel_dir, sub("\\.svg$", "_source_data.csv", panel_file))

  utils::write.csv(per_class_fdr, file.path(STAGE, paste0(stub, "_per_column_term_counts.csv")), row.names = FALSE)
  list(plot = p, data = raw, selected = sel, per_class = per_class_fdr, plotted = plot_df)
}

# ---- Figure 3 F : Memory engram, paired vs unpaired --------------------------
cat("\n=== Figure 3 F : *_paired_veh over *_unpaired_veh ===\n")
f3f <- build_dotplot(
  classes = c("neuropil", "cfos", "mcherry", "neuron"),
  numerator = "paired_veh", denominator = "unpaired_veh",
  plot_title = "Top terms per comparison",
  stub = "Fig3F_memory_engram_GSEA_animal_level",
  panel_dir = PANELS_F3,
  panel_file = "Fig3F_memory_engram_GSEA_animal_level.svg"
)

# ---- Supplementary F : Inhibition, paired CNO vs paired VEH -----------------
cat("\n=== Supplementary F : *_paired_cno over *_paired_veh (3 columns, as in the manuscript) ===\n")
sf3 <- build_dotplot(
  classes = c("neuropil", "cfos", "mcherry"),
  numerator = "paired_cno", denominator = "paired_veh",
  plot_title = "Top terms per comparison",
  stub = "SuppF_inhibition_GSEA_animal_level",
  panel_dir = PANELS_SUP,
  panel_file = "SuppF_inhibition_GSEA_animal_level.svg"
)

cat("\n=== Supplementary F variant : all four sample classes (stage only) ===\n")
sf4 <- build_dotplot(
  classes = c("neuropil", "cfos", "mcherry", "neuron"),
  numerator = "paired_cno", denominator = "paired_veh",
  plot_title = "Top terms per comparison",
  stub = "SuppF_inhibition_GSEA_animal_level_4class",
  panel_dir = STAGE,
  panel_file = "SuppF_inhibition_GSEA_animal_level_4class.svg"
)

# ---- audit -------------------------------------------------------------------
audit <- rbind(
  data.frame(panel = "Figure3 F", f3f$per_class, stringsAsFactors = FALSE),
  data.frame(panel = "Supplementary F", sf3$per_class, stringsAsFactors = FALSE),
  data.frame(panel = "Supplementary F (4-class variant)", sf4$per_class, stringsAsFactors = FALSE)
)
utils::write.csv(audit, file.path(STAGE, "Fig3F_SuppF_gsea_dotplot_audit.csv"), row.names = FALSE)

notes <- c(
  "Ranking statistic changed from log2fc (legacy) to the moderated limma t statistic (canonical). Term membership therefore changes.",
  "Figure 3 F: the mCherry column contributes ZERO terms because mcherry_paired_veh_over_mcherry_unpaired_veh has 0/3829 FDR-significant GO-BP terms at animal level. The column is retained and its (non-significant) NES values are still drawn, so the absence is visible rather than hidden.",
  "Positive NES = enriched toward the canonical numerator (paired_veh for Fig3F, paired_cno for SuppF).",
  "Selection rule, ordering algorithm, palette, size range, theme and dimension formula are the legacy ones, unchanged."
)
writeLines(notes, file.path(STAGE, "Fig3F_SuppF_notes.txt"))

message("Panels F / Supp F done.")
