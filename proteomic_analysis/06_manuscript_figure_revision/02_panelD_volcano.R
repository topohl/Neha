# ==============================================================================
# Figure 3 panel D -- Volcano, mCherry paired_veh vs unpaired_veh
#
# ORIGINAL: 04_differential_expression_enrichment/02_compareGO.r (pre-2026
#           3140-line version), volcano block at legacy line 1422-1484
#           -> 99_historical/compareGO/BP/learning_signature/memory_ensemble/
#              05_Volcano_Plots/Volcano_mcherry2_mcherry4.svg
#           input 99_historical/datasets_mapped/learning_signature/
#              memory_ensemble/mcherry2_mcherry4.csv  (hemisphere-level limma)
#
# CORRECTED: identical plot construction, label-selection rule, palette, theme
#           and figure size; input replaced by the canonical animal-level
#           contrast 02_data/animal_level/mapped/forward/
#              mcherry_paired_veh_vs_mcherry_unpaired_veh.csv
#
# Label rule preserved verbatim from the original: top 5 by log2FC among
# FDR<0.05 up, top 5 by log2FC among FDR<0.05 down. Labels are therefore
# re-derived from the ANIMAL-LEVEL padj; none are carried over.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))
suppressPackageStartupMessages({ library(ggrepel); library(stringr) })

STAGE <- ensure_dir(file.path(FULL, "differential"))

# ---- legacy theme_clean (02_compareGO.r lines 65-110), verbatim -------------
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
      plot.subtitle = element_text(color = "#555555", size = rel(0.95), hjust = 0,
                                   margin = margin(b = 5)),
      strip.text = element_text(color = "#2C2C2C", size = rel(0.95)),
      strip.background = element_rect(color = "#EEEEEE", fill = "#EEEEEE", linewidth = 0.4),
      panel.grid = element_blank(), panel.border = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

make_volcano <- function(d, plot_title, out_stub, panel_dir, panel_file = NULL) {
  v <- d[is.finite(d$log2fc) & is.finite(d$padj), , drop = FALSE]
  v$Gene_Name <- ifelse(is.na(v$mapped_gene_symbol) | !nzchar(v$mapped_gene_symbol),
                        v$gene_symbol, v$mapped_gene_symbol)
  v$Significance <- ifelse(v$padj < FDR & v$log2fc > 0, "up",
                    ifelse(v$padj < FDR & v$log2fc < 0, "down", "n.s."))

  up   <- v[v$Significance == "up", , drop = FALSE]
  down <- v[v$Significance == "down", , drop = FALSE]
  top5_up   <- head(up[order(-up$log2fc), , drop = FALSE], 5)
  top5_down <- head(down[order(down$log2fc), , drop = FALSE], 5)
  lab <- rbind(top5_up, top5_down)

  x_limit <- ceiling(max(abs(v$log2fc), 0, na.rm = TRUE)) + 0.5

  p <- ggplot(v, aes(x = log2fc, y = -log10(padj), label = Gene_Name)) +
    geom_point(aes(color = Significance), shape = 16, alpha = 0.65, size = 2) +
    ggrepel::geom_text_repel(data = lab, aes(label = Gene_Name), size = 3,
                             min.segment.length = 0, max.overlaps = Inf, box.padding = 0.5) +
    scale_color_manual(values = c(up = "#CA0020", down = "#0571B0", `n.s.` = "#CCCCCC"),
                       breaks = c("up", "down", "n.s.")) +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    labs(title = plot_title,
         x = expression(log[2] ~ Fold ~ Change),
         y = expression(-log[10](italic(P)[adj]))) +
    theme_clean(base_size = 9) +
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(FDR), linetype = "dashed", color = "#999999",
               linewidth = 0.4, alpha = 0.7)

  fn <- panel_file %||% paste0(out_stub, ".svg")
  save_panel(p, STAGE, panel_dir, fn, width = 3.5, height = 3.5)

  src <- v[, c("gene_symbol", "mapped_gene_symbol", "Gene_Name", "original_protein_id",
               "Description", "uniprot_accession", "log2fc", "t", "pval", "padj",
               "Significance")]
  src$labelled_in_panel <- src$Gene_Name %in% lab$Gene_Name
  save_source_data(src[order(src$padj), ], STAGE, panel_dir,
                   sub("\\.svg$", "_source_data.csv", fn))

  list(plot = p, data = v, labels = lab,
       n = nrow(v), n_up = nrow(up), n_down = nrow(down))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Figure 3 D --------------------------------------------------------------
d <- read_mapped("mcherry", "paired_veh", "unpaired_veh")
res <- make_volcano(
  d,
  # legacy title format "Volcano plot: <comparison>"; the canonical comparison name
  # is too long for the 3.5 in panel, so the manuscript label is used instead
  plot_title = "Volcano plot: mCherry paired vs unpaired",
  out_stub   = "Fig3D_mCherry_paired_vs_unpaired_volcano_animal_level",
  panel_dir  = PANELS_F3
)

cat(sprintf("\nFig3D  n = %d proteins | FDR<%.2f up = %d, down = %d\n", res$n, FDR, res$n_up, res$n_down))
cat("labelled:", paste(res$labels$Gene_Name, collapse = ", "), "\n")

# ---- audit -------------------------------------------------------------------
legacy <- file.path(DATA_ROOT, "99_historical/datasets_mapped/learning_signature/memory_ensemble/mcherry2_mcherry4.csv")
legacy_n_fdr <- NA_integer_
if (file.exists(legacy)) {
  ld <- utils::read.csv(legacy, stringsAsFactors = FALSE, check.names = FALSE)
  pc <- intersect(c("padj", "adj.P.Val", "p.adjust", "FDR"), names(ld))
  if (length(pc)) legacy_n_fdr <- sum(suppressWarnings(as.numeric(ld[[pc[1]]])) < FDR, na.rm = TRUE)
}

audit <- data.frame(
  item = c("panel", "original_output", "original_input", "original_sampling_unit",
           "original_n_FDR_significant", "corrected_input", "corrected_sampling_unit",
           "corrected_n_tested", "corrected_n_FDR_up", "corrected_n_FDR_down",
           "corrected_labels", "label_rule", "caveat"),
  value = c("Figure3 D",
            "99_historical/compareGO/BP/learning_signature/memory_ensemble/05_Volcano_Plots/Volcano_mcherry2_mcherry4.svg",
            "99_historical/datasets_mapped/learning_signature/memory_ensemble/mcherry2_mcherry4.csv",
            "hemisphere/acquisition (n=6 per group)",
            as.character(legacy_n_fdr),
            "02_data/animal_level/mapped/forward/mcherry_paired_veh_vs_mcherry_unpaired_veh.csv",
            "animal (n=3 per group)",
            as.character(res$n), as.character(res$n_up), as.character(res$n_down),
            paste(res$labels$Gene_Name, collapse = "; "),
            "top 5 by log2FC among FDR<0.05 up and top 5 among FDR<0.05 down (legacy rule, re-applied to animal-level padj)",
            paste("The mCherry paired_veh vs unpaired_veh contrast is 100% aliased with collection plate and shows a broad",
                  "one-directional shift (median log2FC 0.333 over all 5327 proteins) with zero FDR-significant canonical",
                  "GO-BP enrichment. Show the volcano as a magnitude statement only; do not label it a learning programme.")),
  stringsAsFactors = FALSE
)
utils::write.csv(audit, file.path(STAGE, "Fig3D_volcano_audit.csv"), row.names = FALSE)

message("Panel D done.")
