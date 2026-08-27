# ==============================================================================
# FIGURE 3 H and FIGURE 3 I -- FINAL
#
# 3H status: REGENERATED_ANIMAL_LEVEL_DESCRIPTIVE
# 3I status: REGENERATED_ANIMAL_LEVEL_WITH_SHARED_REFERENCE_CAVEAT
#
# Both panels are RETAINED. Nothing is recomputed: values come from the
# source-data CSVs written by scripts/06_..., and the shared-arm control
# correlation is read from the diagnostics table already produced there.
# Only naming, the neutral-title variant and the interpretation audits are new.
#
# HISTORICAL REFERENCES
#   3H  99_historical/compare_sig_expr/plots/scatter_learning_logfc_vs_cno_logfc_mcherry.svg
#       (288 x 432 pt; #2C3E50 @30%, r 4.98 pt; fit #E74C3C 4.27 pt; zero lines
#        #CCCCCC dashed; fonts 20/22 px)  -- historical r = -0.71
#   3I  99_historical/compare_sig_expr/plots/NES_Heatmap_Sorted_By_Comparison_Group.svg
#       (864 x 720 pt; pheatmap of t(nes_matrix): terms on x rotated -45, two
#        value rows, sample-class annotation strip on top; palette
#        c("#2d8be9","#F7F7F7","#d12f42"); annotation colours
#        c(neuropil #4c87c6, cfos #6ccff6, mcherry #faa51a, neuron #fdd700))
#       Term-selection rule RECOVERED from 03_compare_pathways.r: keep terms with
#       p.adjust <= 0.05 in BOTH ensembles within the same sample class.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "10_style_match_common.R"))
suppressPackageStartupMessages(library(patchwork))

FIRST_F3 <- file.path(OUT_ROOT, "figure_panels", "Figure3")
XSTAGE   <- ensure_dir(file.path(FULL, "cross_compartment"))
ESTAGE   <- ensure_dir(file.path(FULL, "enrichment"))

# ==============================================================================
# FIGURE 3 H
# ==============================================================================
h <- utils::read.csv(file.path(FIRST_F3, "Fig3H_logFC_correlation_animal_level_source_data.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
sig <- h[h$learning_signature %in% c(TRUE, "TRUE"), , drop = FALSE]
r_sig <- stats::cor(sig$log2fc_learning, sig$log2fc_cno, use = "complete.obs")

# already-computed diagnostics; NOT rerun
diag <- utils::read.csv(file.path(XSTAGE, "Fig3H_correlation_diagnostics.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
gd <- function(pat) diag$pearson_r[grepl(pat, diag$comparison)][1]
r_all       <- gd("all measured proteins$")
r_ctrl_sig  <- gd("structural control, learning signature")
r_ctrl_all  <- gd("structural control, all measured")
cat(sprintf("3H  n = %d  r = %.2f | all-protein r = %.2f | same-side control r = %+.2f / %+.2f\n",
            nrow(sig), r_sig, r_all, r_ctrl_sig, r_ctrl_all))

pH <- ggplot(sig, aes(x = log2fc_learning, y = log2fc_cno)) +
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
  theme(text = element_text(family = FONT, colour = "#2C2C2C"),
        axis.text  = element_text(size = SCAT_FS_AX, colour = "#2C2C2C"),
        axis.title = element_text(size = SCAT_FS_TIT, colour = "#2C2C2C"),
        panel.grid = element_blank(),
        axis.line  = element_line(colour = "#2C2C2C", linewidth = 0.8 * PT),
        axis.ticks = element_line(colour = "#2C2C2C", linewidth = 0.8 * PT),
        plot.title = element_blank(), plot.subtitle = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 10, 6, 6))

save_sm(pH, SM_F3, "Fig3H_mCherry_logFC_correlation_animal_level_FINAL", SCAT_W, SCAT_H)
utils::write.csv(h, file.path(SM_F3, "Fig3H_source_data.csv"), row.names = FALSE)

audH <- data.frame(
  item = c("panel", "status", "correlation_type", "corrected_r", "n_proteins",
           "x_contrast", "y_contrast", "shared_reference", "shared_arm_sides",
           "pure_noise_expected_r_opposite_sides", "control_contrast",
           "control_shared_arm", "control_shared_arm_sides",
           "control_r_learning_signature", "control_r_all_proteins",
           "r_all_proteins_unselected", "historical_r", "old_sampling_unit",
           "corrected_sampling_unit", "selection_of_points", "control_rerun",
           "inferential_status", "caption_language"),
  value = c(
    "Figure 3 H", "REGENERATED_ANIMAL_LEVEL_DESCRIPTIVE",
    "Pearson product-moment correlation of protein-level moderated-limma log2FC estimates",
    sprintf("%.4f", r_sig), as.character(nrow(sig)),
    "mcherry_paired_veh_over_mcherry_unpaired_veh (animal level)",
    "mcherry_paired_cno_over_mcherry_paired_veh (animal level)",
    "paired_veh",
    "OPPOSITE sides: paired_veh is the numerator of the x contrast and the denominator of the y contrast",
    "-0.5",
    "mcherry_unpaired_cno_over_mcherry_unpaired_veh",
    "unpaired_veh",
    "SAME side: unpaired_veh is the denominator of both",
    sprintf("%+.4f", r_ctrl_sig), sprintf("%+.4f", r_ctrl_all),
    sprintf("%.4f", r_all), "-0.71",
    "hemisphere/acquisition on x (limma) and a raw group-mean difference on y",
    "animal (n = 3 per group), moderated limma log2FC on both axes",
    "points are the 1132 proteins FDR<0.05 in the learning contrast, the historical panel definition; no point was removed",
    "NO - read from full_regenerated/cross_compartment/Fig3H_correlation_diagnostics.csv",
    "DESCRIPTIVE. The negative correlation is structurally expected: the two contrasts share paired_veh on opposite sides, and flipping which arm is shared reverses the sign (same-side control gives r = +0.90). Conditioning the x axis on learning FDR additionally induces regression to the mean on y. This panel cannot on its own establish reversal.",
    "Protein-level effect estimates were negatively correlated across the two contrasts; because both contrasts contain the paired-vehicle condition on opposite sides, this correlation is descriptive and cannot by itself establish reversal of a learning-associated proteomic programme."),
  stringsAsFactors = FALSE)
utils::write.csv(audH, file.path(SM_F3, "Fig3H_interpretation_audit.csv"), row.names = FALSE)
utils::write.csv(audH, file.path(XSTAGE, "Fig3H_interpretation_audit.csv"), row.names = FALSE)

# ==============================================================================
# FIGURE 3 I
# ==============================================================================
long <- utils::read.csv(file.path(FIRST_F3, "Fig3I_CNO_pathway_heatmap_animal_level_source_data.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
cls_levels <- unname(CLASS_LABELS[c("neuropil", "cfos", "mcherry", "neuron")])
long$Comparison_Group <- factor(long$Comparison_Group, levels = cls_levels)
ord <- unique(long[order(long$Comparison_Group, long$Description),
                   c("Description", "Comparison_Group", "sample_class")])
long$Description <- factor(long$Description, levels = ord$Description)
long$Row <- factor(ifelse(long$Dataset == "Memory Ensemble", "Memory Engram", "Inhibition"),
                   levels = c("Inhibition", "Memory Engram"))
lim <- max(abs(long$NES), na.rm = TRUE)

wide <- reshape(long[, c("Description", "Comparison_Group", "Dataset", "NES")],
                idvar = c("Description", "Comparison_Group"), timevar = "Dataset", direction = "wide")
names(wide) <- sub("^NES\\.", "NES_", names(wide))
wide$direction_relation <- ifelse(sign(wide$`NES_Memory Ensemble`) == sign(wide$`NES_Effects Inhibition`),
                                  "concordant", "opposed")
n_terms <- nrow(wide)
n_opp <- sum(wide$direction_relation == "opposed")
cat(sprintf("3I  terms = %d  opposed = %d  concordant = %d  NES limits +/- %.2f\n",
            n_terms, n_opp, n_terms - n_opp, lim))

ann <- ord; ann$Description <- factor(ann$Description, levels = levels(long$Description))
ann$strip <- "Sample class"

heat_theme <- function() {
  theme_minimal(base_size = HEAT_FS_ROW, base_family = FONT) +
    theme(text = element_text(family = FONT, colour = "#2C2C2C"),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          axis.title = element_blank(), axis.ticks = element_blank(),
          legend.key = element_blank(), legend.background = element_blank(),
          legend.title = element_text(size = HEAT_FS_ROW, face = "bold"),
          legend.text = element_text(size = HEAT_FS_ROW),
          legend.position = "right", legend.justification = "top",
          plot.subtitle = element_blank(),
          plot.margin = margin(0, 4, 0, 132))
}

p_ann <- ggplot(ann, aes(x = Description, y = strip, fill = Comparison_Group)) +
  geom_tile(colour = HEAT_BORDER, linewidth = 0.5 * PT * 2) +
  scale_fill_manual(values = setNames(unname(HEAT_ANNOT[c("neuropil", "cfos", "mcherry", "neuron")]),
                                      cls_levels), name = "Sample class", drop = TRUE) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  heat_theme() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                       plot.title = element_blank())

p_heat <- ggplot(long, aes(x = Description, y = Row, fill = NES)) +
  geom_tile(colour = HEAT_BORDER, linewidth = 0.5 * PT * 2) +
  scale_fill_gradientn(colours = grDevices::colorRampPalette(c(HEAT_LOW, HEAT_MID, HEAT_HIGH),
                                                             space = "Lab")(100),
                       limits = c(-lim, lim), name = "NES",
                       guide = guide_colourbar(barwidth = unit(3, "mm"), barheight = unit(20, "mm"))) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0), position = "right") +
  heat_theme() +
  theme(axis.text.x = element_text(size = HEAT_FS_TERM, angle = 45, hjust = 1, vjust = 1, colour = "#2C2C2C"),
        axis.text.y = element_text(size = HEAT_FS_ROW, hjust = 0, colour = "#2C2C2C"),
        plot.title = element_blank())

build_I <- function(title = NULL) {
  top <- if (is.null(title)) p_ann else
    p_ann + ggtitle(title) +
    theme(plot.title = element_text(family = FONT, size = 13, face = "bold", hjust = 0.5,
                                    margin = margin(b = 6)))
  top / p_heat + plot_layout(heights = c(HEAT_ANN_H, 2 * HEAT_CELL), guides = "collect") &
    theme(legend.position = "right")
}

W <- max(9, 2.6 + n_terms * 0.155 + 2.6)
save_sm(build_I("CNO and pairing-associated pathway contrasts"), SM_F3,
        "Fig3I_pathway_heatmap_animal_level_FINAL", W, 4.8)
save_sm(build_I(NULL), SM_F3,
        "Fig3I_pathway_heatmap_animal_level_FINAL_no_title", W, 4.6)

utils::write.csv(long, file.path(SM_F3, "Fig3I_source_data.csv"), row.names = FALSE)
utils::write.csv(wide, file.path(ESTAGE, "Fig3I_term_direction_table.csv"), row.names = FALSE)

ovl <- utils::read.csv(file.path(ESTAGE, "Fig3I_overlap_summary.csv"), stringsAsFactors = FALSE)
audI <- data.frame(
  item = c("panel", "status", "row_1_definition", "row_2_definition", "shared_arm",
           "shared_arm_sides", "ranking_statistic", "enrichment_source",
           "term_selection_rule", "term_selection_rule_provenance",
           "n_overlapping_terms", "n_opposed", "n_concordant",
           "per_class_overlap", "NES_scale_limits", "old_sampling_unit",
           "corrected_sampling_unit", "old_term_count", "enrichment_rerun",
           "visual_style", "warning", "neutral_title"),
  value = c(
    "Figure 3 I", "REGENERATED_ANIMAL_LEVEL_WITH_SHARED_REFERENCE_CAVEAT",
    "Memory Engram = *_paired_veh over *_unpaired_veh (positive NES = enriched toward paired_veh)",
    "Inhibition = *_paired_cno over *_paired_veh (positive NES = enriched toward paired_cno)",
    "paired_veh",
    "OPPOSITE sides: paired_veh is the numerator of row 1 and the denominator of row 2",
    "moderated limma t statistic (canonical); the historical panel ranked by log2FC",
    "03_output/enrichment/enrichment_t_rank_validation_20260825/per_comparison/*/GSEA_GO_BP.csv",
    "terms with p.adjust <= 0.05 in BOTH contrasts within the same sample class",
    "RECOVERED verbatim from 04_differential_expression_enrichment/03_compare_pathways.r (pre-2026 273-line version); not invented, and not tuned to maximise opposition",
    as.character(n_terms), as.character(n_opp), as.character(n_terms - n_opp),
    paste(sprintf("%s: learning %d, cno %d, overlap %d", ovl$sample_class,
                  ovl$n_FDR_learning, ovl$n_FDR_cno, ovl$n_overlap), collapse = "; "),
    sprintf("+/- %.2f", lim),
    "hemisphere/acquisition, GSEA ranked by log2FC", "animal, GSEA ranked by moderated t",
    "32 terms in the historical panel", "NO - canonical GSEA tables used as computed",
    "matches the historical pheatmap: terms on x rotated -45, two value rows, sample-class annotation strip on top, palette c(#2d8be9, #F7F7F7, #d12f42), white borders",
    "Cross-row opposition CANNOT be interpreted as proof that CNO reverses a pairing-associated programme. Because the two rows share paired_veh on opposite sides, any pathway high in paired_veh is pushed up in row 1 and down in row 2 by construction. The uniform opposition observed here is the expected structural consequence of that design.",
    "CNO and pairing-associated pathway contrasts"),
  stringsAsFactors = FALSE)
utils::write.csv(audI, file.path(SM_F3, "Fig3I_interpretation_audit.csv"), row.names = FALSE)
utils::write.csv(audI, file.path(ESTAGE, "Fig3I_interpretation_audit.csv"), row.names = FALSE)

message("Figures 3 H and 3 I FINAL written.")
