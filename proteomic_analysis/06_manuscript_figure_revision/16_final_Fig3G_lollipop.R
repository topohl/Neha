# ==============================================================================
# FIGURE 3 G -- FINAL. mCherry paired_veh vs unpaired_veh, GO-BP lollipop.
# Status: REGENERATED_ANIMAL_LEVEL_NULL_ENRICHMENT
#
# The panel is RETAINED, not retired. It now shows the corrected animal-level
# null: the broad mCherry protein-level difference is NOT accompanied by
# significant GO-BP over-representation.
#
# NO ENRICHMENT IS RERUN. The canonical ORA table is read as-is. It is the full
# tested GO-BP result because the canonical run used pvalue_cutoff = 1 and
# qvalue_cutoff = 1 (audits/run_parameters.csv), so nothing was pre-filtered away.
#
# ------------------------------------------------------------------------------
# HISTORICAL DESIGN, recovered from
#   script  04_differential_expression_enrichment/02_compareGO.r (pre-2026
#           3140-line version), lollipop block at legacy line 2236
#   output  99_historical/compareGO/BP/learning_signature/memory_ensemble/
#           07_Regulated_Protein_GO/GO_Lollipop_mcherry2_mcherry4_Upregulated_BP.svg
#
#   analysis      clusterProfiler::enrichGO (ORA) over the FDR-significant
#                 UPregulated protein list, keyType UNIPROT, ont BP, BH,
#                 minGSSize 10
#   selection     Count >= 10, then arrange(desc(GeneRatio)), slice_head(n = 15)
#                 -- but applied to a table that had ALREADY been filtered to
#                 p < 0.05 / q < 0.2 by enrichGO itself
#   display       y ordered so the largest Gene Ratio sits at the top
#   canvas        432 x 360 pt (6 x 5 in)
#   segment       #BDBDBD, linewidth 1.5
#   points        size = Count, range c(2, 6), stroke 0.64, colour = p.adjust
#   palette       recovered numerically from the embedded legend colourbar PNG:
#                 c("#E64F4F", "#EFA24A", "#DFB74B"), fit 0.0024
#                 (the git version of the script says red-white-blue; the file on
#                 disk and the manuscript both use this red-orange-khaki ramp)
#   axes          x "Gene Ratio", y "GO Term"; term labels 9 px, wrapped at 50
#   title         "Upregulated GO - mcherry2_mcherry4" -- omitted here, the
#                 manuscript sets the panel heading in Illustrator
#
# ------------------------------------------------------------------------------
# SELECTION RULE USED HERE, AND WHY IT DIFFERS
#
# The historical rule ranked by Gene Ratio *within an already significant set*.
# At animal level no term is significant, so ranking by Gene Ratio alone would
# simply pick the largest gene sets and would not be the same operation.
# The faithful analogue is therefore:
#     filter Count >= 10   (historical)
#     select the 15 terms with the lowest corrected adjusted P
#     order the y axis by descending Gene Ratio   (historical display rule)
# This is the "lowest adjusted P" option; it is recorded in the source data and
# in the audit. No term is presented as significant.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "10_style_match_common.R"))
suppressPackageStartupMessages(library(stringr))

CMP  <- "mcherry_paired_veh_over_mcherry_unpaired_veh"
ORA  <- file.path(ENRICH_ROOT, CMP, "ORA_GO_BP_fdr_up.csv")
STAGE <- ensure_dir(file.path(FULL, "enrichment"))

LOLLI_RAMP  <- c("#E64F4F", "#EFA24A", "#DFB74B")
LOLLI_SEG   <- "#BDBDBD"
LOLLI_SEG_LW<- 1.5
LOLLI_SIZE  <- c(2, 6)
LOLLI_W     <- 6
LOLLI_H     <- 5
LOLLI_FS    <- 9
N_TERMS     <- 15
MIN_COUNT   <- 10

theme_clean <- function(base_size = 9, base_family = FONT) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "#2C2C2C", family = base_family, lineheight = 1.3),
      axis.line.x = element_line(color = "#2C2C2C", linewidth = 0.6 * PT),
      axis.line.y = element_line(color = "#2C2C2C", linewidth = 0.6 * PT),
      axis.ticks = element_line(color = "#2C2C2C", linewidth = 0.5 * PT),
      axis.ticks.length = unit(3 * PT, "mm"),
      axis.text = element_text(color = "#2C2C2C", size = rel(0.95)),
      axis.title = element_text(color = "#2C2C2C", size = rel(1.0)),
      legend.background = element_blank(), legend.key = element_rect(color = NA, fill = NA),
      legend.key.size = unit(10 * PT, "mm"),
      legend.title = element_text(color = "#2C2C2C", size = rel(0.95)),
      legend.text = element_text(color = "#2C2C2C", size = rel(0.9)),
      legend.position = "right", legend.justification = "top",
      panel.grid = element_blank(), panel.border = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ---- read the canonical ORA table, unchanged --------------------------------
ora <- utils::read.csv(ORA, stringsAsFactors = FALSE, check.names = FALSE)
ora$Count <- as.numeric(ora$Count)
ora$GeneRatioNum <- vapply(strsplit(as.character(ora$GeneRatio), "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
query_size <- unique(vapply(strsplit(as.character(ora$GeneRatio), "/"),
                            function(x) as.numeric(x[2]), numeric(1)))[1]

n_fdr <- sum(ora$p.adjust < FDR, na.rm = TRUE)
min_padj <- min(ora$p.adjust, na.rm = TRUE)
cat(sprintf("corrected ORA: %d terms tested, %d at FDR<%.2f, min adjusted P = %.4f, query size = %d\n",
            nrow(ora), n_fdr, FDR, min_padj, query_size))
stopifnot(n_fdr == 0L)

pool <- ora[ora$Count >= MIN_COUNT, , drop = FALSE]
sel  <- head(pool[order(pool$p.adjust, pool$pvalue, pool$ID), , drop = FALSE], N_TERMS)
sel  <- sel[order(-sel$GeneRatioNum), , drop = FALSE]           # historical display order
sel$DescriptionWrapped <- stringr::str_wrap(sel$Description, width = 50)
sel$DescriptionWrapped <- factor(sel$DescriptionWrapped, levels = rev(unique(sel$DescriptionWrapped)))

cat("selected terms (top", N_TERMS, "by corrected adjusted P, Count >=", MIN_COUNT, "):\n")
print(sel[, c("Description", "Count", "GeneRatioNum", "p.adjust")], row.names = FALSE)

# ---- plot --------------------------------------------------------------------
build <- function(annotate_null) {
  p <- ggplot(sel, aes(x = GeneRatioNum, y = DescriptionWrapped)) +
    geom_segment(aes(x = 0, xend = GeneRatioNum, y = DescriptionWrapped, yend = DescriptionWrapped),
                 colour = LOLLI_SEG, linewidth = LOLLI_SEG_LW * PT * 2) +
    geom_point(aes(size = Count, colour = p.adjust), alpha = 0.85, stroke = 0.3) +
    scale_colour_gradientn(colours = grDevices::colorRampPalette(LOLLI_RAMP, space = "Lab")(100),
                           name = expression(italic(P)[adj])) +
    scale_size_continuous(range = LOLLI_SIZE, name = "Gene Count") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)), name = "Gene Ratio") +
    labs(y = "GO Term") +
    theme_clean(base_size = LOLLI_FS)
  if (annotate_null) {
    # unobtrusive caption note, not a redesign of the graph itself
    p <- p +
      labs(caption = sprintf("No term reaches FDR < %.2f (minimum adjusted P = %.3f).", FDR, min_padj)) +
      theme(plot.caption = element_text(family = FONT, size = 7, colour = "#6E6E6E",
                                        hjust = 1, margin = margin(t = 6)))
  }
  p
}

save_sm(build(TRUE),  SM_F3, "Fig3G_mCherry_GO_lollipop_animal_level_FINAL", LOLLI_W, LOLLI_H)
save_sm(build(FALSE), SM_F3, "Fig3G_mCherry_GO_lollipop_animal_level_FINAL_no_annotation", LOLLI_W, LOLLI_H)

# ---- source data -------------------------------------------------------------
src <- sel[, c("ID", "Description", "GeneRatio", "GeneRatioNum", "BgRatio", "Count",
               "FoldEnrichment", "pvalue", "p.adjust", "qvalue", "geneID")]
src$plotted <- TRUE
src$n_terms_tested <- nrow(ora)
src$n_FDR_significant <- n_fdr
src$minimum_adjusted_P <- min_padj
src$query_gene_set_size <- query_size
src$selection_rule <- sprintf("Count >= %d, then %d lowest corrected adjusted P; y ordered by descending Gene Ratio (historical display rule)",
                              MIN_COUNT, N_TERMS)
sd_path <- file.path(SM_F3, "Fig3G_source_data.csv")
utils::write.csv(src, sd_path, row.names = FALSE)
utils::write.csv(src, file.path(STAGE, "Fig3G_source_data.csv"), row.names = FALSE)

# ---- audit -------------------------------------------------------------------
audit <- data.frame(
  item = c("panel", "status", "old_analysis", "old_sampling_unit", "old_n_FDR_significant",
           "corrected_analysis", "corrected_sampling_unit", "corrected_input",
           "n_terms_tested", "n_FDR_significant", "minimum_adjusted_P", "query_gene_set_size",
           "enrichment_rerun", "selection_rule", "visual_style", "inferential_caveat",
           "preferred_interpretation"),
  value = c(
    "Figure 3 G", "REGENERATED_ANIMAL_LEVEL_NULL_ENRICHMENT",
    "clusterProfiler::enrichGO ORA over the hemisphere-level FDR-significant upregulated mcherry2_mcherry4 protein list, ont BP, BH, minGSSize 10",
    "hemisphere/acquisition (L and R treated as independent)",
    "15 terms plotted, all below the enrichGO p<0.05 / q<0.2 cutoffs",
    "clusterProfiler::enrichGO ORA over the ANIMAL-LEVEL FDR-significant upregulated protein list; canonical table read as computed, nothing rerun",
    "animal (n = 3 per group)",
    file.path("03_output/enrichment/enrichment_t_rank_validation_20260825/per_comparison", CMP, "ORA_GO_BP_fdr_up.csv"),
    as.character(nrow(ora)), as.character(n_fdr), sprintf("%.4f", min_padj), as.character(query_size),
    "NO - canonical ORA used as-is (run used pvalue_cutoff = 1, qvalue_cutoff = 1, so the table is the full tested result)",
    sprintf("Count >= %d, then the %d terms with the lowest corrected adjusted P; y axis ordered by descending Gene Ratio as in the historical panel", MIN_COUNT, N_TERMS),
    "matches the historical lollipop: 6 x 5 in, #BDBDBD segments, size = Count range 2-6, p.adjust ramp #E64F4F/#EFA24A/#DFB74B recovered from the historical legend colourbar, x 'Gene Ratio', y 'GO Term', no embedded title",
    "NONE of the plotted terms is FDR-significant. Colour encodes adjusted P on the historical ramp, so the legend numbers (not the hues) carry the significance information; the annotated variant states the null explicitly.",
    "The broad mCherry protein-level difference was not accompanied by significant GO-BP over-representation after animal-level correction."),
  stringsAsFactors = FALSE)
utils::write.csv(audit, file.path(SM_F3, "Fig3G_interpretation_audit.csv"), row.names = FALSE)
utils::write.csv(audit, file.path(STAGE, "Fig3G_interpretation_audit.csv"), row.names = FALSE)

message("Figure 3 G FINAL written. n FDR significant = ", n_fdr, ", min adjusted P = ", round(min_padj, 4))
