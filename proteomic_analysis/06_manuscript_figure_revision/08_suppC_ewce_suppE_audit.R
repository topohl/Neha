# ==============================================================================
# Supplementary panel C -- Cell-type marker / identity heatmap (EWCE Baseline)
# Supplementary panel E -- Cross-compartment bubble plot: CLASSIFICATION AUDIT
#
# ORIGINAL C: SOURCE_NOT_FOUND. The EWCE version that wrote
#             99_historical/ewce_legacy/Plots/Indiv_Heatmap_Baseline.svg
#             (run log Run_Log_20260408_1322.txt) is absent from the repository
#             and from every git commit. On-plot title "Significant Subtypes
#             Landscape (Baseline)", legend "Z-score", rows = Zeisel ctd_allKI
#             annotation-level-2 cell types, columns = bg_1 .. neuron_4.
#             Current successor: 05_celltype_enrichment_EWCE/01_EWCE.r, which
#             still computes AnalysisType == "Baseline" at animal level.
#
# ORIGINAL E: 04_differential_expression_enrichment/02_compareGO.r (pre-2026
#             3140-line version) with ensemble_profiling =
#             "baseline_cell_type_profiling", condition = "CS"
#             -> 99_historical/compareGO/BP/baseline_cell_type_profiling/CS/
#                02_Main_Plots/Dotplot_Enrichment_TopGenes_PerComp.svg
#             comparisons bg2_neuron2, cfos2_neuron2, mcherry2_cfos2,
#             mcherry2_neuron2 -- i.e. CROSS-COMPARTMENT, not one of the 12
#             canonical within-compartment treatment contrasts.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))
suppressPackageStartupMessages({ library(readxl) })
source(file.path(REPO_ROOT, "R", "analysis_labels.R"))
source(file.path(REPO_ROOT, "R", "protigy_input_utils.R"))

STAGE_Q <- ensure_dir(file.path(FULL, "qc"))
STAGE_X <- ensure_dir(file.path(FULL, "cross_compartment"))

# ==============================================================================
# SUPPLEMENTARY C -- EWCE Baseline cell-type landscape, animal level
# ==============================================================================
tbl <- file.path(EWCE_ROOT, "02_Tables_Supplements", "Supplementary_Table_EWCE.xlsx")
pr  <- as.data.frame(readxl::read_excel(tbl, sheet = "Primary_Results"))
b   <- pr[pr$AnalysisType == "Baseline", , drop = FALSE]
stopifnot(nrow(b) > 0, length(unique(b$TopN)) == 1L, length(unique(b$AnnotLevel)) == 1L)
cat(sprintf("EWCE Baseline: %d rows, %d targets, %d cell types, TopN %s, annotation level %s\n",
            nrow(b), length(unique(b$Target)), length(unique(b$CellType)),
            unique(b$TopN), unique(b$AnnotLevel)))

# legacy row rule: keep cell types that are significant somewhere ("Significant Subtypes")
sig_ct <- sort(unique(b$CellType[b$Significant_Global]))
cat("significant cell types:", length(sig_ct), ":", paste(sig_ct, collapse = ", "), "\n")

h <- b[b$CellType %in% sig_ct, , drop = FALSE]
h$sample_class <- normalize_sample_class(h$Stratum)
h$condition    <- normalize_condition(h$Metric)
h$Column <- paste(unname(CLASS_LABELS[h$sample_class]), unname(COND_LABELS[h$condition]), sep = " ")

col_levels <- unlist(lapply(SAMPLE_CLASSES, function(cl)
  paste(unname(CLASS_LABELS[cl]), unname(COND_LABELS[CONDITIONS]), sep = " ")))
h$Column <- factor(h$Column, levels = col_levels[col_levels %in% h$Column])

# order rows the way the legacy panel did: cell types grouped, strongest last
ord <- stats::aggregate(sd_from_mean ~ CellType, data = h, FUN = function(z) max(abs(z)))
h$CellType <- factor(h$CellType, levels = ord$CellType[order(ord$sd_from_mean)])

lim <- max(abs(h$sd_from_mean), na.rm = TRUE)
p_c <- ggplot(h, aes(x = Column, y = CellType, fill = sd_from_mean)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_point(data = h[h$Significant_Global, , drop = FALSE], shape = 8, size = 1.1, colour = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                       limits = c(-lim, lim), name = "Z-score") +
  labs(title = "Significant Subtypes Landscape (Baseline)",
       subtitle = sprintf("EWCE marker over-representation, animal level; annotation level %s, top %s, asterisk = global FDR < %.2f",
                          unique(b$AnnotLevel), unique(b$TopN), FDR),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 6)),
    plot.subtitle = element_text(hjust = 0.5, size = 8, colour = "#555555", margin = margin(b = 8)),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

save_panel(p_c, STAGE_Q, PANELS_SUP, "SuppC_celltype_identity_heatmap_animal_level.svg",
           width = 11, height = 5.5)
save_source_data(
  h[, c("CellType", "Column", "sample_class", "condition", "Target", "TopN", "AnnotLevel",
        "fold_change", "sd_from_mean", "p", "q_target", "q_global", "Significant_Global")],
  STAGE_Q, PANELS_SUP, "SuppC_celltype_identity_heatmap_animal_level_source_data.csv")

legacy_ct <- c("CA1Pyr1", "CA1Pyr2", "CA2Pyr2", "Int12", "Int15",
               "Oligo1", "Oligo2", "Oligo3", "Oligo4", "Oligo5", "Oligo6",
               "S1PyrDL", "S1PyrL23", "S1PyrL6")
audit_c <- data.frame(
  item = c("panel", "status", "original_output", "original_generating_script",
           "original_sampling_unit", "corrected_input", "corrected_sampling_unit",
           "annotation_level", "topN", "n_significant_celltypes_corrected",
           "significant_celltypes_corrected", "significant_celltypes_legacy",
           "retained_from_legacy", "lost_vs_legacy", "gained_vs_legacy", "caption_requirement"),
  value = c("Supplementary C", "REGENERATED_ANIMAL_LEVEL",
            "99_historical/ewce_legacy/Plots/Indiv_Heatmap_Baseline.svg",
            "SOURCE_NOT_FOUND (EWCE run of 2026-04-08; reconstructed from its output)",
            "hemisphere/acquisition-derived per-condition protein lists",
            "03_output/ewce/EWCE_Results_animal_level_validation_20260825/02_Tables_Supplements/Supplementary_Table_EWCE.xlsx (Primary_Results, AnalysisType Baseline)",
            "animal-derived per-condition protein lists",
            as.character(unique(b$AnnotLevel)), as.character(unique(b$TopN)),
            as.character(length(sig_ct)),
            paste(sig_ct, collapse = "; "), paste(legacy_ct, collapse = "; "),
            paste(intersect(sig_ct, legacy_ct), collapse = "; "),
            paste(setdiff(legacy_ct, sig_ct), collapse = "; "),
            paste(setdiff(sig_ct, legacy_ct), collapse = "; "),
            "EWCE tests over-representation of cell-type marker genes in a ranked protein list. It does NOT measure cell abundance, activation or state."),
  stringsAsFactors = FALSE)
utils::write.csv(audit_c, file.path(STAGE_Q, "SuppC_ewce_baseline_audit.csv"), row.names = FALSE)
cat("retained:", paste(intersect(sig_ct, legacy_ct), collapse = ", "), "\n")
cat("lost vs legacy:", paste(setdiff(legacy_ct, sig_ct), collapse = ", "), "\n")

# ==============================================================================
# SUPPLEMENTARY E -- classification audit, no model is fitted here
# ==============================================================================
parsed <- validate_protigy_gct_v13(GCT_ANIMAL)
cm <- parsed$column_metadata
meta <- data.frame(
  Sample = colnames(cm),
  AnimalID = trimws(as.character(cm["AnimalID", ])),
  condition = normalize_condition(cm["condition", ]),
  sample_class = normalize_sample_class(cm["sample_class", ]),
  stringsAsFactors = FALSE)

pairs <- list(
  c("neuropil", "neuron"), c("cfos", "neuron"), c("mcherry", "cfos"), c("mcherry", "neuron")
)
overlap <- do.call(rbind, lapply(pairs, function(p) {
  a <- meta$AnimalID[meta$sample_class == p[1] & meta$condition == "paired_veh"]
  b2 <- meta$AnimalID[meta$sample_class == p[2] & meta$condition == "paired_veh"]
  data.frame(
    legacy_comparison = c(`neuropil|neuron` = "bg2_neuron2", `cfos|neuron` = "cfos2_neuron2",
                          `mcherry|cfos` = "mcherry2_cfos2", `mcherry|neuron` = "mcherry2_neuron2")[[paste(p, collapse = "|")]],
    arm_1 = p[1], arm_2 = p[2], condition = "paired_veh",
    n_animals_arm_1 = length(a), n_animals_arm_2 = length(b2),
    n_shared_animals = length(intersect(a, b2)),
    shared_animals = paste(sort(intersect(a, b2)), collapse = ";"),
    percent_shared = 100 * length(intersect(a, b2)) / max(length(a), length(b2)),
    stringsAsFactors = FALSE)
}))
print(overlap)
utils::write.csv(overlap, file.path(STAGE_X, "SuppE_cross_compartment_animal_overlap.csv"), row.names = FALSE)

scope <- utils::read.csv(file.path(DATA_ROOT, "02_data/animal_level/input_gct/contrast_scope_audit.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
utils::write.csv(scope, file.path(STAGE_X, "SuppE_project_contrast_scope_audit_copy.csv"), row.names = FALSE)

lines <- c(
  "================================================================================",
  "SUPPLEMENTARY PANEL E  --  CLASSIFICATION: REQUIRES_NEW_SECONDARY_MODEL",
  "Bubble plot, 'Cellular Identities Memory Engram'",
  paste("Audit written", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "================================================================================",
  "",
  "1. WHAT STATISTICAL MODEL GENERATED THESE CROSS-COMPARTMENT COMPARISONS",
  "",
  "   Plot     99_historical/compareGO/BP/baseline_cell_type_profiling/CS/",
  "            02_Main_Plots/Dotplot_Enrichment_TopGenes_PerComp.svg",
  "   Script   04_differential_expression_enrichment/02_compareGO.r (pre-2026",
  "            3140-line version) run with ensemble_profiling =",
  "            'baseline_cell_type_profiling' and condition = 'CS'. That script",
  "            only READS enrichment CSVs; it fits no model itself.",
  "   Upstream 99_historical/core_enrichment/BP/baseline_cell_type_profiling/CS/*.csv",
  "            <- clusterProfiler GSEA over",
  "            99_historical/datasets_mapped/baseline_cell_type_profiling/CS/",
  "              bg2_neuron2.csv, cfos2_neuron2.csv, mcherry2_cfos2.csv,",
  "              mcherry2_neuron2.csv",
  "            <- ProTigy ORDINARY TWO-SAMPLE moderated-t on the hemisphere-level",
  "               matrix (02_data/gct/pg.matrix_Two-sample_mod_T_*_n120x5349.gct).",
  "",
  "   So the displayed statistic is a two-independent-groups moderated t test",
  "   between two tissue compartments, taken at condition 2 (paired_veh).",
  "",
  "2. WERE HEMISPHERE OBSERVATIONS TREATED AS INDEPENDENT?",
  "",
  "   Yes. The input matrix has 96 acquisitions = 12 animals x 4 compartments x 2",
  "   hemispheres, and the two-sample test used them as independent replicates.",
  "   That is the same pseudoreplication the animal-level correction removes",
  "   everywhere else.",
  "",
  "   But this panel has a SECOND, independent violation that the 12 canonical",
  "   contrasts do not have:",
  "",
  paste0("     ", sprintf("%-18s %s vs %s : %d and %d animals, %d SHARED (%.0f%%) -- %s",
                          overlap$legacy_comparison, overlap$arm_1, overlap$arm_2,
                          overlap$n_animals_arm_1, overlap$n_animals_arm_2,
                          overlap$n_shared_animals, overlap$percent_shared, overlap$shared_animals)),
  "",
  "   Both arms of every one of these four comparisons are dissected from THE SAME",
  "   ANIMALS. They are paired observations, not two independent groups. An",
  "   unpaired test on fully paired data mis-states the standard error in a way",
  "   that no amount of L/R averaging fixes: averaging hemispheres corrects the",
  "   n, it does not introduce the within-animal blocking factor.",
  "",
  "3. IS THERE AN APPROPRIATE ANIMAL-LEVEL / PAIRED ANALOGUE ALREADY AVAILABLE?",
  "",
  "   No. The canonical animal-level branch contains exactly 12 contrasts, all",
  "   WITHIN a single compartment:",
  "     paired_cno vs paired_veh, paired_veh vs unpaired_veh,",
  "     unpaired_cno vs unpaired_veh   x   {neuropil, cfos, mcherry, neuron}",
  "   No cross-compartment contrast exists at animal level, in",
  "   02_data/animal_level/mapped/ or in",
  "   03_output/enrichment/enrichment_t_rank_validation_20260825/.",
  "",
  "   The project's own scope audit already recorded this decision. From",
  "   02_data/animal_level/input_gct/contrast_scope_audit.csv:",
  "",
  "     candidate_contrast              : cross_sample_class_same_condition",
  "     include_in_primary_contrast_manifest : FALSE",
  "     evidence    : historical baseline profiling/all-pairwise outputs",
  "     disposition : requires paired/animal-aware modeling; do not use ordinary",
  "                   independent two-sample ProTigy",
  "",
  "4. CLASSIFICATION",
  "",
  "     REQUIRES_NEW_SECONDARY_MODEL",
  "",
  "   Not VALID_AS_DESCRIPTIVE: the panel's encoding is inferential throughout --",
  "   dot colour is NES and dot size is -log10(adjusted P). There is no way to",
  "   read it descriptively.",
  "   Not REGENERATABLE_ANIMAL_LEVEL: the required contrast does not exist in the",
  "   canonical animal-level outputs, and producing it means choosing and fitting",
  "   a new model.",
  "   Not RETIRED: compartment identity is a real and reportable property of this",
  "   dataset; the panel is recoverable, it simply needs a model that has not been",
  "   run yet.",
  "",
  "5. WHAT A CORRECT MODEL WOULD LOOK LIKE (NOT RUN -- AWAITING A DECISION)",
  "",
  "   Fit on the animal-level matrix (02_data/animal_level/input_gct/",
  "   neha_protigy_input_animal_level_primary.gct), restricted to condition",
  "   paired_veh, with compartment as the tested factor and AnimalID as a blocking",
  "   factor. Either:",
  "     (a) limma with duplicateCorrelation(block = AnimalID), or",
  "     (b) a paired contrast on within-animal differences (compartment A minus",
  "         compartment B per animal), then a one-sample moderated t.",
  "   n = 3 animals per compartment in the paired_veh stratum, so the paired",
  "   difference has 2 residual degrees of freedom before moderation. Power is",
  "   low and any result must be framed accordingly.",
  "   GSEA would then be re-ranked on the resulting moderated t, matching the",
  "   canonical convention, and the bubble plot rebuilt with the SAME code path",
  "   used for Figure 3 F and Supplementary F.",
  "",
  "6. DECISION REQUESTED BEFORE ANY CODE IS RUN",
  "",
  "   Per the task instruction, no new secondary model has been created. Choose:",
  "     (i)   run the within-animal paired model above and rebuild the panel;",
  "     (ii)  replace the panel with a purely descriptive compartment-identity",
  "           display that carries no p-values -- Supplementary panel C (EWCE",
  "           Baseline) and Supplementary panel D (rank abundance) already serve",
  "           this purpose and are both corrected;",
  "     (iii) retire the panel.",
  "",
  "   The old hemisphere-level statistics must not be reused under any option.",
  "",
  "FILES WRITTEN BY THIS AUDIT",
  "  SuppE_cross_compartment_animal_overlap.csv    animal sharing between the arms",
  "  SuppE_project_contrast_scope_audit_copy.csv   the project's own prior ruling",
  "================================================================================"
)
writeLines(lines, file.path(STAGE_X, "SuppE_REQUIRES_NEW_SECONDARY_MODEL.txt"))
file.copy(file.path(STAGE_X, "SuppE_REQUIRES_NEW_SECONDARY_MODEL.txt"),
          file.path(PANELS_SUP, "SuppE_REQUIRES_NEW_SECONDARY_MODEL.txt"), overwrite = TRUE)

message("Supplementary C regenerated; Supplementary E audit written.")
