# ==============================================================================
# Figure 3 panel F      -- GSEA bubble grid, paired_veh vs unpaired_veh
#                          columns neuropil / cFos / mCherry / neuron
# Supplementary panel F -- GSEA bubble grid, paired_cno vs paired_veh
#                          columns neuropil / cFos / mCherry / neuron (NeuN column
#                          added 2026-11-25; the manuscript panel had three)
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
#           (03_output/enrichment/enrichment_within_compartment_20261125/GO_BP/,
#           fgsea nPermSimple = 100000).
#
# Blank cells: a term is drawn in every column whose GSEA table returns it. A
# cell is blank when the table does not; the audit counts those cells and says
# whether fgsea could not test the set (GO_BP/audits/<comparison>/
# gsea_untested_sets.csv) or the set is missing for another reason.
#
# Selection rule preserved verbatim from the original:
#   significant_only = TRUE (p.adjust < 0.05)
#   best instance per Description per comparison by |NES|
#   per comparison: top 5 by NES among positives + top 5 by NES among negatives
#   union of those Descriptions plotted for EVERY column, significant or not
#   row/column order from order_dotplot(): hclust on the NES matrix for columns,
#   then rows by dominant column, direction and |NES|
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/figure_rebuild_within_compartment_20261125"
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

# gene sets in the size limits that fgsea gave no P value (primary ranking only)
read_untested <- function(sample_class, numerator, denominator) {
  cmp <- sprintf("%s_%s_over_%s_%s", sample_class, numerator, sample_class, denominator)
  f <- file.path(dirname(ENRICH_ROOT), "audits", cmp, "gsea_untested_sets.csv")
  if (!file.exists(f)) stop("Untested-set audit not found: ", f, call. = FALSE)
  u <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  u[u$analysis == "GO_BP_GSEA", , drop = FALSE]
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
  untested <- lapply(stats::setNames(classes, classes), read_untested,
                     numerator = numerator, denominator = denominator)

  if (!nrow(sel)) {
    warning("No FDR-significant terms in any column for ", plot_title, "; no panel written.")
    return(list(plot = NULL, data = raw, selected = sel, per_class = NULL))
  }

  plot_df <- best[best$Description %in% unique(sel$Description), , drop = FALSE]
  # cells left blank: a plotted term that a column's GSEA table does not return
  blank <- do.call(rbind, lapply(classes, function(cl) {
    miss <- setdiff(unique(sel$Description), plot_df$Description[plot_df$sample_class == cl])
    if (!length(miss)) return(NULL)
    data.frame(Comparison = unname(CLASS_LABELS[cl]), sample_class = cl, Description = miss,
               reason = ifelse(miss %in% untested[[cl]]$Description,
                               "gene set not tested by fgsea (gsea_untested_sets.csv)",
                               "not in the GSEA table for another reason"),
               stringsAsFactors = FALSE)
  }))
  if (is.null(blank)) blank <- data.frame(Comparison = character(0), sample_class = character(0),
                                          Description = character(0), reason = character(0))

  per_class_fdr <- do.call(rbind, lapply(classes, function(cl) {
    cells <- plot_df[plot_df$sample_class == cl, , drop = FALSE]
    b <- blank[blank$sample_class == cl, , drop = FALSE]
    data.frame(
      sample_class = cl,
      column = unname(CLASS_LABELS[cl]),
      n_terms_tested = sum(raw$sample_class == cl),
      n_terms_untested = nrow(untested[[cl]]),
      n_FDR_significant = sum(raw$sample_class == cl & raw$p.adjust < FDR, na.rm = TRUE),
      n_terms_contributed_to_panel = sum(sel$sample_class == cl),
      n_terms_plotted = length(unique(sel$Description)),
      n_cells_drawn = nrow(cells),
      n_cells_FDR_significant = sum(is.finite(cells$p.adjust) & cells$p.adjust < FDR),
      n_cells_blank = nrow(b),
      n_cells_blank_untested = sum(grepl("not tested by fgsea", b$reason)),
      n_cells_blank_other = sum(!grepl("not tested by fgsea", b$reason)),
      stringsAsFactors = FALSE)
  }))
  print(per_class_fdr)
  if (nrow(blank)) print(blank, row.names = FALSE)
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

  utils::write.csv(per_class_fdr, fresh(file.path(STAGE, paste0(stub, "_per_column_term_counts.csv"))), row.names = FALSE)
  utils::write.csv(blank, fresh(file.path(STAGE, paste0(stub, "_blank_cells.csv"))), row.names = FALSE)
  list(plot = p, data = raw, selected = sel, per_class = per_class_fdr, plotted = plot_df, blank = blank)
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
# All four compartments (2026-11-25); the earlier panel left out the neuron column.
cat("\n=== Supplementary F : *_paired_cno over *_paired_veh (4 columns) ===\n")
sf4 <- build_dotplot(
  classes = c("neuropil", "cfos", "mcherry", "neuron"),
  numerator = "paired_cno", denominator = "paired_veh",
  plot_title = "Top terms per comparison",
  stub = "SuppF_inhibition_GSEA_animal_level",
  panel_dir = PANELS_SUP,
  panel_file = "SuppF_inhibition_GSEA_animal_level.svg"
)

# ---- audit -------------------------------------------------------------------
audit <- rbind(
  data.frame(panel = "Figure3 F", f3f$per_class, stringsAsFactors = FALSE),
  data.frame(panel = "Supplementary F", sf4$per_class, stringsAsFactors = FALSE)
)
utils::write.csv(audit, fresh(file.path(STAGE, "Fig3F_SuppF_gsea_dotplot_audit.csv")), row.names = FALSE)

count_line <- function(nm, pc) {
  sprintf(paste0("%s: %d terms plotted, %d cells drawn (%d at adjusted P < %.2f), %d blank ",
                 "(%d not tested by fgsea, %d other). Per column: %s."),
          nm, pc$n_terms_plotted[1], sum(pc$n_cells_drawn), sum(pc$n_cells_FDR_significant), FDR,
          sum(pc$n_cells_blank), sum(pc$n_cells_blank_untested), sum(pc$n_cells_blank_other),
          paste(sprintf("%s %d tested, %d untested, %d significant, %d cells, %d blank",
                        pc$column, pc$n_terms_tested, pc$n_terms_untested, pc$n_FDR_significant,
                        pc$n_cells_drawn, pc$n_cells_blank), collapse = "; "))
}
notes <- c(
  "Ranking statistic changed from log2fc (legacy) to the moderated limma t statistic (canonical). Term membership therefore changes.",
  "GSEA input: enrichment_within_compartment_20261125 (GO-BP, fgsea nPermSimple = 100000).",
  count_line("Figure3 F", f3f$per_class),
  count_line("Supplementary F", sf4$per_class),
  "Blank cells are listed with their reason in *_blank_cells.csv; counts per column are in Fig3F_SuppF_gsea_dotplot_audit.csv.",
  "Positive NES = enriched toward the canonical numerator (paired_veh for Fig3F, paired_cno for SuppF).",
  "Selection rule, ordering algorithm, palette, size range, theme and dimension formula are the legacy ones, unchanged."
)
writeLines(notes, fresh(file.path(STAGE, "Fig3F_SuppF_notes.txt")))

message("Panels F / Supp F done.")
