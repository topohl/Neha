# ==============================================================================
# Figure 3 panel H -- Learning logFC vs CNO logFC scatter, mCherry (DESCRIPTIVE)
# Figure 3 panel I -- NES heatmap, learning vs CNO pathway profiles
#
# ORIGINAL H: 04_differential_expression_enrichment/04_compare_sig_expr.r
#             (pre-2026 360-line version), "Plot F" block at legacy line 340-361
#             -> 99_historical/compare_sig_expr/plots/
#                scatter_learning_logfc_vs_cno_logfc_mcherry.svg   (r = -0.71)
#             x: hemisphere-level limma log2FC for mcherry2_mcherry4, restricted
#                to the FDR-significant learning signature
#             y: RAW group-mean log2 difference (paired CNO minus paired VEH)
#                computed inside the script from the hemisphere-level GCT
#
# ORIGINAL I: 04_differential_expression_enrichment/03_compare_pathways.r
#             (pre-2026 273-line version), ggplot heatmap at legacy line 197-244
#             -> 99_historical/compare_sig_expr/plots/NES_Absolute_Heatmap_ggplot.svg
#             terms FDR<=0.05 in BOTH the learning and the CNO ensemble, per
#             sample class; fill = NES; facet rows = sample class
#
# CORRECTED: both axes / both columns now come from the canonical animal-level
#            moderated-limma contrasts and the canonical moderated-t-ranked GSEA.
#            Plot type, palette, geometry, facet structure and figure size are
#            the originals.
#
# INTERPRETATION CHANGE (applies to BOTH panels, and is the reason they are
# classified REGENERATED_WITH_INTERPRETATION_CHANGE):
#   The two contrasts SHARE paired_veh, and share it on OPPOSITE sides --
#     learning = paired_veh / unpaired_veh   (paired_veh is the numerator)
#     CNO      = paired_cno / paired_veh     (paired_veh is the denominator)
#   Any protein or pathway that is high in paired_veh is therefore pushed up on
#   one axis and down on the other BY CONSTRUCTION, whatever the biology. A
#   negative correlation in panel H, and uniformly opposed NES signs in panel I,
#   are the expected structural consequence of that design and are NOT evidence
#   that chemogenetic inhibition reverses a learning programme.
#   Panel H additionally conditions on learning significance, which adds
#   regression-to-the-mean on the y axis.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))
suppressPackageStartupMessages({ library(tibble) })

STAGE_X <- ensure_dir(file.path(FULL, "cross_compartment"))
STAGE_D <- ensure_dir(file.path(FULL, "differential"))
STAGE_E <- ensure_dir(file.path(FULL, "enrichment"))

# ==============================================================================
# PANEL H
# ==============================================================================
learn <- read_mapped("mcherry", "paired_veh", "unpaired_veh")
cno   <- read_mapped("mcherry", "paired_cno", "paired_veh")
# Structural control. There is no canonical contrast that shares NO arm with the
# learning contrast, so the control instead flips WHICH arm is shared and on which
# side:  learning = paired_veh - unpaired_veh
#        panel H  = paired_cno - paired_veh    (shares paired_veh, OPPOSITE sides)
#        control  = unpaired_cno - unpaired_veh(shares unpaired_veh, SAME side)
# Under pure independent noise of equal variance these give r = -0.5 and r = +0.5
# respectively. If the observed correlations bracket zero in that pattern, the sign
# is being set by the shared arm, not by biology.
ctrl  <- read_mapped("mcherry", "unpaired_cno", "unpaired_veh")

key <- c("uniprot_accession", "original_protein_id", "mapped_gene_symbol")
j <- merge(learn[, c(key, "log2fc", "t", "padj")],
           cno[,   c("uniprot_accession", "log2fc", "t", "padj")],
           by = "uniprot_accession", suffixes = c("_learning", "_cno"))
j <- merge(j, setNames(ctrl[, c("uniprot_accession", "log2fc", "padj")],
                       c("uniprot_accession", "log2fc_unpaired_cno", "padj_unpaired_cno")),
           by = "uniprot_accession")
j$learning_signature <- is.finite(j$padj_learning) & j$padj_learning < FDR

sig <- j[j$learning_signature, , drop = FALSE]

r_sig  <- stats::cor(sig$log2fc_learning, sig$log2fc_cno, use = "complete.obs")
r_all  <- stats::cor(j$log2fc_learning,   j$log2fc_cno,   use = "complete.obs")
r_ctrl_sig <- stats::cor(sig$log2fc_learning, sig$log2fc_unpaired_cno, use = "complete.obs")
r_ctrl_all <- stats::cor(j$log2fc_learning,   j$log2fc_unpaired_cno,   use = "complete.obs")

cat(sprintf("\nPanel H correlations (mCherry, animal level)\n"))
cat(sprintf("  learning signature (n = %5d, legacy panel definition) : r = %+.2f\n", nrow(sig), r_sig))
cat(sprintf("  all measured proteins (n = %5d)                       : r = %+.2f\n", nrow(j), r_all))
cat(sprintf("  STRUCTURAL CONTROL - unpaired_cno vs unpaired_veh, which shares unpaired_veh\n"))
cat(sprintf("  as a common DENOMINATOR (same side) instead of paired_veh on opposite sides:\n"))
cat(sprintf("    on the learning signature                            : r = %+.2f\n", r_ctrl_sig))
cat(sprintf("    on all measured proteins                             : r = %+.2f\n", r_ctrl_all))
cat(sprintf("  pure-noise expectation: -0.5 for a shared arm on opposite sides, +0.5 for the same side\n"))
cat(sprintf("  legacy hemisphere-level value shown in the manuscript  : r = -0.71\n\n"))

scatter_panel <- function(d, xv, yv, xlab, ylab, ttl, rval, n) {
  ggplot(d, aes(x = .data[[xv]], y = .data[[yv]])) +
    geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
    geom_point(color = "#2c3e50", alpha = 0.3, size = 4, shape = 16) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#e74c3c", fill = "#e74c3c", alpha = 0.2) +
    labs(x = xlab, y = ylab, title = ttl,
         subtitle = sprintf("DESCRIPTIVE. Pearson correlation: r = %.2f (n = %d)", rval, n)) +
    # legacy used base_size 22 with short axis labels; reduced here so the longer,
    # unambiguous contrast names fit inside the same panel footprint
    theme_minimal(base_size = 14, base_family = "sans") +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
}

p_h <- scatter_panel(
  sig, "log2fc_learning", "log2fc_cno",
  xlab = "Learning Log2FC (paired_veh / unpaired_veh)",
  ylab = "CNO vs VEH Log2FC (paired)",
  ttl  = "Learning signature vs CNO effect - mCherry",
  rval = r_sig, n = nrow(sig)
)
save_panel(p_h, STAGE_X, PANELS_F3, "Fig3H_logFC_correlation_animal_level.svg", width = 6, height = 7)

src_h <- j[, c("uniprot_accession", "original_protein_id", "mapped_gene_symbol",
               "log2fc_learning", "t_learning", "padj_learning",
               "log2fc_cno", "t_cno", "padj_cno",
               "log2fc_unpaired_cno", "padj_unpaired_cno", "learning_signature")]
save_source_data(src_h[order(src_h$padj_learning), ], STAGE_X, PANELS_F3,
                 "Fig3H_logFC_correlation_animal_level_source_data.csv")

# stage-only companions: the unselected version and the no-shared-reference control
p_all <- scatter_panel(j, "log2fc_learning", "log2fc_cno",
                       "Learning Log2FC (paired_veh / unpaired_veh)",
                       "CNO vs VEH Log2FC (paired)",
                       "All measured proteins (no selection on learning FDR)",
                       r_all, nrow(j))
ggsave(file.path(STAGE_X, "Fig3H_variant_all_proteins.svg"), p_all, width = 6, height = 7, device = svglite::svglite)

p_ctrl <- scatter_panel(sig, "log2fc_learning", "log2fc_unpaired_cno",
                        "Learning Log2FC (paired_veh / unpaired_veh)",
                        "CNO vs VEH Log2FC (unpaired stratum)",
                        "Structural control: shared arm on the SAME side",
                        r_ctrl_sig, nrow(sig))
ggsave(file.path(STAGE_X, "Fig3H_control_shared_arm_same_side.svg"), p_ctrl, width = 6, height = 7, device = svglite::svglite)

corr_tbl <- data.frame(
  comparison = c("learning vs paired-CNO, learning signature (legacy panel definition)",
                 "learning vs paired-CNO, all measured proteins",
                 "learning vs unpaired-CNO structural control, learning signature",
                 "learning vs unpaired-CNO structural control, all measured proteins"),
  shared_arm = c("paired_veh", "paired_veh", "unpaired_veh", "unpaired_veh"),
  shared_arm_side = c("opposite sides (numerator vs denominator)",
                      "opposite sides (numerator vs denominator)",
                      "same side (denominator in both)",
                      "same side (denominator in both)"),
  pure_noise_expected_r = c(-0.5, -0.5, 0.5, 0.5),
  conditioned_on_learning_FDR = c(TRUE, FALSE, TRUE, FALSE),
  n = c(nrow(sig), nrow(j), nrow(sig), nrow(j)),
  pearson_r = c(r_sig, r_all, r_ctrl_sig, r_ctrl_all),
  stringsAsFactors = FALSE
)
utils::write.csv(corr_tbl, file.path(STAGE_X, "Fig3H_correlation_diagnostics.csv"), row.names = FALSE)
print(corr_tbl)

# ==============================================================================
# PANEL I
# ==============================================================================
classes <- c("neuropil", "cfos", "mcherry", "neuron")
per_class <- lapply(classes, function(cl) {
  L <- read_gsea(cl, "paired_veh", "unpaired_veh")
  C <- read_gsea(cl, "paired_cno", "paired_veh")
  Ls <- L[is.finite(L$p.adjust) & L$p.adjust < FDR, , drop = FALSE]
  Cs <- C[is.finite(C$p.adjust) & C$p.adjust < FDR, , drop = FALSE]
  ov <- intersect(Ls$Description, Cs$Description)          # legacy rule: FDR in BOTH
  if (!length(ov)) {
    return(list(rows = NULL, summary = data.frame(
      sample_class = cl, n_FDR_learning = nrow(Ls), n_FDR_cno = nrow(Cs),
      n_overlap = 0L, n_concordant = 0L, n_opposed = 0L, stringsAsFactors = FALSE)))
  }
  m <- merge(Ls[Ls$Description %in% ov, c("ID", "Description", "NES", "p.adjust")],
             Cs[Cs$Description %in% ov, c("Description", "NES", "p.adjust")],
             by = "Description", suffixes = c("_learning", "_cno"))
  m$Comparison_Group <- unname(CLASS_LABELS[cl])
  m$sample_class <- cl
  list(rows = m, summary = data.frame(
    sample_class = cl, n_FDR_learning = nrow(Ls), n_FDR_cno = nrow(Cs),
    n_overlap = nrow(m),
    n_concordant = sum(sign(m$NES_learning) == sign(m$NES_cno)),
    n_opposed    = sum(sign(m$NES_learning) != sign(m$NES_cno)),
    stringsAsFactors = FALSE))
})
i_summary <- do.call(rbind, lapply(per_class, `[[`, "summary"))
print(i_summary)
utils::write.csv(i_summary, file.path(STAGE_E, "Fig3I_overlap_summary.csv"), row.names = FALSE)

rows <- do.call(rbind, Filter(Negate(is.null), lapply(per_class, `[[`, "rows")))
if (is.null(rows) || !nrow(rows)) stop("No FDR-overlapping terms; panel I cannot be drawn.", call. = FALSE)

long <- rbind(
  data.frame(Description = rows$Description, Comparison_Group = rows$Comparison_Group,
             sample_class = rows$sample_class, Dataset = "Memory Ensemble",
             NES = rows$NES_learning, p.adjust = rows$p.adjust_learning, stringsAsFactors = FALSE),
  data.frame(Description = rows$Description, Comparison_Group = rows$Comparison_Group,
             sample_class = rows$sample_class, Dataset = "Effects Inhibition",
             NES = rows$NES_cno, p.adjust = rows$p.adjust_cno, stringsAsFactors = FALSE)
)
long$Dataset <- factor(long$Dataset, levels = c("Memory Ensemble", "Effects Inhibition"))
ord <- rows[order(rows$Comparison_Group, rows$Description), ]
long$Description <- factor(long$Description, levels = rev(unique(ord$Description)))
long$Comparison_Group <- factor(long$Comparison_Group,
                                levels = unname(CLASS_LABELS[classes[classes %in% rows$sample_class]]))

p_i <- ggplot(long, aes(x = Dataset, y = Description, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#2d8be9", mid = "#F7F7F7", high = "#d12f42", midpoint = 0, name = "NES") +
  facet_grid(Comparison_Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
    strip.placement = "outside",
    panel.spacing = unit(0.1, "lines"),
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "#555555"),
    panel.grid = element_blank()
  ) +
  labs(title = "CNO-associated pathway changes in paired animals",
       subtitle = paste("GO-BP terms FDR<0.05 in both contrasts, animal level.",
                        "\nMemory Ensemble = paired_veh / unpaired_veh   |   Effects Inhibition = paired_cno / paired_veh",
                        "\nThe two contrasts share paired_veh on opposite sides, so opposed NES signs are expected by construction."))

h_in <- max(6, 2 + 0.22 * length(unique(long$Description)))
save_panel(p_i, STAGE_E, PANELS_F3, "Fig3I_CNO_pathway_heatmap_animal_level.svg", width = 8, height = h_in)
save_source_data(long, STAGE_E, PANELS_F3, "Fig3I_CNO_pathway_heatmap_animal_level_source_data.csv")
utils::write.csv(rows, file.path(STAGE_E, "Fig3I_overlapping_pathways_animal_level.csv"), row.names = FALSE)

audit <- data.frame(
  item = c("panelH_original_r_hemisphere_level",
           "panelH_corrected_r_learning_signature",
           "panelH_corrected_r_all_proteins",
           "panelH_control_r_shared_arm_same_side_signature",
           "panelH_control_r_shared_arm_same_side_all",
           "panelH_status",
           "panelI_total_overlapping_terms",
           "panelI_concordant_terms",
           "panelI_opposed_terms",
           "panelI_classes_with_zero_overlap",
           "panelI_title_change"),
  value = c("-0.71",
            sprintf("%.2f", r_sig), sprintf("%.2f", r_all),
            sprintf("%.2f", r_ctrl_sig), sprintf("%.2f", r_ctrl_all),
            "descriptive only; shared paired_veh reference and selection on learning FDR both inflate the anti-correlation",
            as.character(sum(i_summary$n_overlap)),
            as.character(sum(i_summary$n_concordant)),
            as.character(sum(i_summary$n_opposed)),
            paste(i_summary$sample_class[i_summary$n_overlap == 0], collapse = "; "),
            "'Inhibition of Memory Engram' -> 'CNO-associated pathway changes in paired animals'"),
  stringsAsFactors = FALSE
)
utils::write.csv(audit, file.path(STAGE_X, "Fig3H_Fig3I_audit.csv"), row.names = FALSE)

message("Panels H and I done.")
