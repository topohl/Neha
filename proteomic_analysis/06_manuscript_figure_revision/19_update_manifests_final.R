# ==============================================================================
# Update every manifest after the 2026-08-27 correction that RETAINS
# Figure 3 G / 3 H / 3 I and Supplementary E.
# No statistic is touched; this script only edits status/provenance records and
# moves the two now-superseded retirement notes out of the deliverable folders.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "10_style_match_common.R"))

MAN   <- file.path(OUT_ROOT, "manifests")
SMMAN <- file.path(SM, "manifests")

NEW_STATUS <- c(
  "Figure3|G"                  = "REGENERATED_ANIMAL_LEVEL_NULL_ENRICHMENT",
  "Figure3|H"                  = "REGENERATED_ANIMAL_LEVEL_DESCRIPTIVE",
  "Figure3|I"                  = "REGENERATED_ANIMAL_LEVEL_WITH_SHARED_REFERENCE_CAVEAT",
  "Supplementary_proteomics|E" = "REGENERATED_SECONDARY_PAIRED_ANIMAL_LEVEL")

key_of <- function(fig, pan) paste(fig, pan, sep = "|")

# ---- 1. original_figure_panel_source_map.csv ---------------------------------
p1 <- file.path(MAN, "original_figure_panel_source_map.csv")
x1 <- utils::read.csv(p1, stringsAsFactors = FALSE, check.names = FALSE)
k1 <- key_of(x1$figure, x1$panel)
hit <- k1 %in% names(NEW_STATUS)
x1$regeneration_status[hit] <- unname(NEW_STATUS[k1[hit]])
x1$notes[k1 == "Figure3|G"] <- paste(
  "RETAINED (2026-08-27 correction). The panel is regenerated as the corrected animal-level NULL:",
  "ORA GO-BP over the animal-level FDR-significant up list returns 0 of 3639 terms at FDR<0.05,",
  "minimum adjusted P 0.231. The 15 terms with the lowest corrected adjusted P (Count >= 10) are plotted",
  "in the historical lollipop style, ordered by descending Gene Ratio as in the original, with a caption",
  "stating the null. Enrichment was NOT rerun: the canonical run used pvalue_cutoff = 1 and qvalue_cutoff = 1,",
  "so the stored table is the full tested result. Supersedes the earlier retirement decision.")
x1$notes[k1 == "Supplementary_proteomics|E"] <- paste(
  "RETAINED (2026-08-27 correction). A secondary paired animal-level model was authorised for this panel only.",
  "Historical design recovered and empirically verified: condition paired_veh, numerator = first-named compartment,",
  "gseGO BP over log2FC ranks. Corrected model: limma ~ 0 + sample_class + AnimalID on the 12 paired_veh",
  "animal-level columns (3 animals x 4 compartments), full rank 6/6, residual df 6, eBayes df.total 27.4;",
  "GSEA re-ranked on the paired-model moderated t. Isolated under",
  "full_regenerated/cross_compartment/SuppE_secondary_paired/. Supersedes the earlier BLOCKED decision.")
utils::write.csv(x1, p1, row.names = FALSE)
cat("updated:", basename(p1), "\n")
print(x1[hit, c("figure", "panel", "regeneration_status")], row.names = FALSE)

# ---- 2. regenerated_panel_manifest.csv ---------------------------------------
p2 <- file.path(MAN, "regenerated_panel_manifest.csv")
x2 <- utils::read.csv(p2, stringsAsFactors = FALSE, check.names = FALSE)
k2 <- key_of(x2$figure, x2$panel)
hit2 <- k2 %in% names(NEW_STATUS)
x2$status[hit2] <- unname(NEW_STATUS[k2[hit2]])

setv <- function(df, k, col, val) { df[[col]][key_of(df$figure, df$panel) == k] <- val; df }

x2 <- setv(x2, "Figure3|G", "corrected_source_script", "scripts/16_final_Fig3G_lollipop.R")
x2 <- setv(x2, "Figure3|G", "new_output_path", "figure_panels_style_matched/Figure3/Fig3G_mCherry_GO_lollipop_animal_level_FINAL.svg")
x2 <- setv(x2, "Figure3|G", "source_data_path", "figure_panels_style_matched/Figure3/Fig3G_source_data.csv")
x2 <- setv(x2, "Figure3|G", "visual_layout_preserved", "YES")
x2 <- setv(x2, "Figure3|G", "key_numbers_new", "3639 GO-BP terms tested, 0 at FDR<0.05, minimum adjusted P 0.231; 15 terms plotted (Count >= 10, lowest adjusted P)")
x2 <- setv(x2, "Figure3|G", "illustrator_action", "Drop in. Caption must state that no term reaches FDR<0.05.")

x2 <- setv(x2, "Figure3|H", "corrected_source_script", "scripts/17_final_Fig3H_Fig3I.R")
x2 <- setv(x2, "Figure3|H", "new_output_path", "figure_panels_style_matched/Figure3/Fig3H_mCherry_logFC_correlation_animal_level_FINAL.svg")
x2 <- setv(x2, "Figure3|H", "source_data_path", "figure_panels_style_matched/Figure3/Fig3H_source_data.csv")

x2 <- setv(x2, "Figure3|I", "corrected_source_script", "scripts/17_final_Fig3H_Fig3I.R")
x2 <- setv(x2, "Figure3|I", "new_output_path", "figure_panels_style_matched/Figure3/Fig3I_pathway_heatmap_animal_level_FINAL.svg")
x2 <- setv(x2, "Figure3|I", "source_data_path", "figure_panels_style_matched/Figure3/Fig3I_source_data.csv")

x2 <- setv(x2, "Supplementary_proteomics|E", "corrected_source_script", "scripts/18_suppE_secondary_paired_model.R")
x2 <- setv(x2, "Supplementary_proteomics|E", "corrected_canonical_input",
           "02_data/animal_level/input_gct/neha_protigy_input_animal_level_primary.gct restricted to paired_veh (SECONDARY model, isolated)")
x2 <- setv(x2, "Supplementary_proteomics|E", "new_output_path", "figure_panels_style_matched/Supplementary_proteomics/SuppE_cellular_identities_animal_level_FINAL.svg")
x2 <- setv(x2, "Supplementary_proteomics|E", "source_data_path", "figure_panels_style_matched/Supplementary_proteomics/SuppE_source_data.csv")
x2 <- setv(x2, "Supplementary_proteomics|E", "new_sampling_unit", "animal (3 animals x 4 compartments, within-animal paired)")
x2 <- setv(x2, "Supplementary_proteomics|E", "visual_layout_preserved", "YES")
x2 <- setv(x2, "Supplementary_proteomics|E", "displayed_statistics_changed", "YES")
x2 <- setv(x2, "Supplementary_proteomics|E", "biological_interpretation_changed", "YES")
x2 <- setv(x2, "Supplementary_proteomics|E", "key_numbers_new",
           "limma ~ 0 + sample_class + AnimalID; rank 6/6; residual df 6; eBayes df.total 27.4; FDR proteins 460/0/2574/1077; FDR GO-BP terms 1023/14/1/2; 19 terms plotted")
x2 <- setv(x2, "Supplementary_proteomics|E", "illustrator_action", "Drop in. Caption must say secondary paired cross-compartment analysis, n = 3 animals, not a treatment effect.")
utils::write.csv(x2, p2, row.names = FALSE)
cat("updated:", basename(p2), "\n")
print(x2[hit2, c("figure", "panel", "status")], row.names = FALSE)

# ---- 3. visual_fidelity_audit.csv --------------------------------------------
p3 <- file.path(SMMAN, "visual_fidelity_audit.csv")
x3 <- utils::read.csv(p3, stringsAsFactors = FALSE, check.names = FALSE)
k3 <- key_of(x3$figure, x3$panel)
setv3 <- function(df, k, col, val) { df[[col]][key_of(df$figure, df$panel) == k] <- val; df }

x3 <- setv3(x3, "Figure3|G", "new_input_file", "03_output/enrichment/enrichment_t_rank_validation_20260825/per_comparison/mcherry_paired_veh_over_mcherry_unpaired_veh/ORA_GO_BP_fdr_up.csv")
x3 <- setv3(x3, "Figure3|G", "style_matched_output", "figure_panels_style_matched/Figure3/Fig3G_mCherry_GO_lollipop_animal_level_FINAL.svg")
x3 <- setv3(x3, "Figure3|G", "visual_fidelity_before", "n/a - panel absent in the first pass")
x3 <- setv3(x3, "Figure3|G", "visual_fidelity_after", "EXACT_OR_NEAR_EXACT")
x3 <- setv3(x3, "Figure3|G", "palette_match", "EXACT")
x3 <- setv3(x3, "Figure3|G", "layout_match", "EXACT")
x3 <- setv3(x3, "Figure3|G", "label_match", "EXACT")
x3 <- setv3(x3, "Figure3|G", "legend_match", "EXACT")
x3 <- setv3(x3, "Figure3|G", "remaining_difference", "x-axis Gene Ratio extends to 0.19 rather than 0.02 because the corrected up-list is larger; a caption line states the null")
x3 <- setv3(x3, "Figure3|G", "notes", "RETAINED as the corrected null. 6 x 5 in, #BDBDBD segments, size = Count range 2-6, p.adjust ramp c(#E64F4F,#EFA24A,#DFB74B) recovered numerically from the historical legend colourbar (fit 0.0024), x 'Gene Ratio', y 'GO Term', no embedded title. A _no_annotation variant is also provided.")

x3 <- setv3(x3, "Supplementary_proteomics|E", "new_input_file", "full_regenerated/cross_compartment/SuppE_secondary_paired/enrichment/GSEA_GO_BP_all_contrasts.csv")
x3 <- setv3(x3, "Supplementary_proteomics|E", "style_matched_output", "figure_panels_style_matched/Supplementary_proteomics/SuppE_cellular_identities_animal_level_FINAL.svg")
x3 <- setv3(x3, "Supplementary_proteomics|E", "underlying_statistic_verified", "YES - secondary paired animal-level limma model, GSEA re-ranked on its moderated t; historical design recovered and empirically verified")
x3 <- setv3(x3, "Supplementary_proteomics|E", "visual_fidelity_before", "n/a - panel absent in the first pass")
x3 <- setv3(x3, "Supplementary_proteomics|E", "visual_fidelity_after", "EXACT_OR_NEAR_EXACT")
x3 <- setv3(x3, "Supplementary_proteomics|E", "palette_match", "EXACT")
x3 <- setv3(x3, "Supplementary_proteomics|E", "layout_match", "EXACT")
x3 <- setv3(x3, "Supplementary_proteomics|E", "label_match", "EXACT")
x3 <- setv3(x3, "Supplementary_proteomics|E", "legend_match", "EXACT")
x3 <- setv3(x3, "Supplementary_proteomics|E", "remaining_difference", "term membership differs because the ranking statistic is the paired-model moderated t rather than hemisphere-level log2FC")
x3 <- setv3(x3, "Supplementary_proteomics|E", "notes", "RETAINED via an authorised secondary paired model. Same bubble grammar as the historical dotplot: blue-white-orange NES ramp c(#4C87C6, white, #FAA51A), size = -log10 adjusted P range 1.5-5, terms wrapped at 50, columns in the historical order. A _titled variant carrying 'Cellular Identities Memory Engram' is also provided.")

x3 <- setv3(x3, "Figure3|H", "style_matched_output", "figure_panels_style_matched/Figure3/Fig3H_mCherry_logFC_correlation_animal_level_FINAL.svg")
x3 <- setv3(x3, "Figure3|I", "style_matched_output", "figure_panels_style_matched/Figure3/Fig3I_pathway_heatmap_animal_level_FINAL.svg")
utils::write.csv(x3, p3, row.names = FALSE)
cat("updated:", basename(p3), "\n")

# ---- 4. move the superseded retirement notes out of the deliverable folders ---
banner <- c(
  "################################################################################",
  "SUPERSEDED 2026-08-27.",
  "The panel described below is NO LONGER retired / blocked. It has been retained and",
  "regenerated at animal level. See:",
  "  Figure 3 G   -> figure_panels_style_matched/Figure3/Fig3G_mCherry_GO_lollipop_animal_level_FINAL.svg",
  "  Supplementary E -> figure_panels_style_matched/Supplementary_proteomics/SuppE_cellular_identities_animal_level_FINAL.svg",
  "The evidence recorded below remains factually correct and is kept as the audit trail",
  "for WHY the corrected result is null (3G) and why a secondary paired model was needed (Supp E).",
  "################################################################################",
  "")
moves <- list(
  list(from = file.path(SM, "Figure3", "Fig3G_RETIRED_AFTER_ANIMAL_LEVEL_CORRECTION.txt"),
       to   = file.path(SMMAN, "SUPERSEDED_Fig3G_retirement_note.txt")),
  list(from = file.path(SM, "Supplementary_proteomics", "SuppE_REQUIRES_NEW_SECONDARY_MODEL.txt"),
       to   = file.path(SMMAN, "SUPERSEDED_SuppE_requires_new_model_note.txt")))
for (mv in moves) {
  if (file.exists(mv$from)) {
    writeLines(c(banner, readLines(mv$from, warn = FALSE)), mv$to)
    file.remove(mv$from)
    cat("superseded ->", basename(mv$to), "\n")
  }
}

message("Manifests updated.")
