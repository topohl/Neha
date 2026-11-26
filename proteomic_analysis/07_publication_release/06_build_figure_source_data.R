#!/usr/bin/env Rscript

# Publication release, stage 06 -- figure source data and the panel provenance map.
#
# Produces
#   editor_source_data/figure_source_map.tsv          one row per manuscript panel
#   editor_source_data/figure_source_data/**          the per-panel source-data files,
#                                                     copied byte for byte (G and Supp. F
#                                                     have two each)
#
# Built from the 2026-08-27 revision manifests:
#   manifests/regenerated_panel_manifest.csv          what was done, with the numbers
#   manifests/original_figure_panel_source_map.csv    where the original panel came from
#   06_manuscript_figure_revision/SCRIPT_PROVENANCE.csv  executed vs repository script hashes
#
# Where the records disagree or were later replaced:
#   * FINAL_REPORT.md (11:59) retired Figure 3 G and blocked Supplementary E; the later
#     RETAINED_PANELS_3G_3H_3I_SuppE_REPORT.md (14:13) reversed both, and
#     regenerated_panel_manifest.csv (14:11) records that outcome, so the manifest is used.
#   * 3G / 3H / 3I / Supp E have a first-pass and a final script; both are recorded
#     (generating_script, final_revision_script).
#   * Fig. 2 D, F, G, H, I and Supp. Fig. 4 F (keys Figure3 D/F/G/H/I and
#     Supplementary_proteomics F) were rebuilt on the per-compartment statistics, last on
#     2026-11-25 in 03_output/figure_rebuild_within_compartment_20261125 (GO GSEA at fgsea
#     nPermSimple = 100000, KEGG unchanged). REBUILT_PANELS overrides the manifest for these
#     six; the rebuild has no manifest, so their numbers and notes are read from its
#     source-data and audit files. G is now a GSEA lollipop of GO-BP and GO-CC terms instead
#     of ORA. Changes of 2026-11-25: G rows ordered by adjusted P within each ontology,
#     Supp. Fig. 4 F with all four sample classes, I with one column per term and sample
#     class, D and H axes in SD units. The 2026-11-24 rebuild is superseded.
#
# panel_label and the released file names use the manuscript numbering (Figure 2,
# Supplementary Figure 4: Fig2D_source_data.csv, SuppFig4F_source_data.csv). Internal keys
# (Figure3, Supplementary_proteomics) keep the old names, and
# panel_source_data_original_filename keeps each canonical file name.
#
# The scripts under 06_manuscript_figure_revision/ (including its within_compartment_*
# copies) are snapshots of what ran and are only read here.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_validation.R"))
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_presentation.R"))
release_require("digest")

DATA_ROOT <- release_data_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/06_build_figure_source_data.R"

release_banner("stage 06 -- figure source data")

REVISION_ROOT <- file.path(DATA_ROOT, "03_output", "reviewer_revision_animal_level_20260827")
REBUILD_ROOT <- release_figure_rebuild_root(DATA_ROOT)
REBUILD_REL <- paste(RELEASE_STATISTICS_RUN$figure_rebuild_relative, collapse = "/")
STATS_REL <- paste(RELEASE_STATISTICS_RUN$stats_relative, collapse = "/")
ENRICH_REL <- paste(RELEASE_STATISTICS_RUN$enrichment_relative, collapse = "/")
# generation_pass of the current rebuild in SCRIPT_PROVENANCE.csv; earlier rebuild passes
# ("within_compartment_rebuild", 2026-11-24) are superseded and map no panel.
REBUILD_PASS <- "within_compartment_rebuild_20261125"
REBUILD_STYLE <- "figure_panels_style_matched/Figure3"

paths <- list(
  regenerated = file.path(REVISION_ROOT, "manifests", "regenerated_panel_manifest.csv"),
  original_map = file.path(REVISION_ROOT, "manifests", "original_figure_panel_source_map.csv"),
  script_provenance = file.path(REPO_ROOT, "06_manuscript_figure_revision",
                                "SCRIPT_PROVENANCE.csv"),
  final_report = file.path(REVISION_ROOT, "manifests", "FINAL_REPORT.md"),
  retained_report = file.path(REVISION_ROOT, "manifests",
                              "RETAINED_PANELS_3G_3H_3I_SuppE_REPORT.md"),
  rebuild_record = release_figure_rebuild_record(DATA_ROOT),
  fig_h_audit = file.path(REBUILD_ROOT, REBUILD_STYLE, "Fig3H_interpretation_audit.csv"),
  fig_i_audit = file.path(REBUILD_ROOT, REBUILD_STYLE, "Fig3I_interpretation_audit.csv"),
  enrichment_index = file.path(release_enrichment_root("GO_BP", DATA_ROOT),
                               "indexEnrichmentComparisons.csv"),
  ora_learning_up = file.path(release_enrichment_root("GO_BP", DATA_ROOT), "per_comparison",
                              "mcherry_paired_veh_over_mcherry_unpaired_veh",
                              "ORA_GO_BP_fdr_up.csv"),
  ora_learning_up_cc = file.path(release_enrichment_root("GO_CC", DATA_ROOT), "per_comparison",
                                 "mcherry_paired_veh_over_mcherry_unpaired_veh",
                                 "ORA_GO_CC_fdr_up.csv"),
  contrast_manifest = release_path("metadata", "primary_contrast_manifest.tsv",
                                   create_dir = FALSE)
)
for (nm in names(paths)) release_assert_exists(paths[[nm]], nm)

regenerated <- release_read_csv(paths$regenerated)
original_map <- release_read_csv(paths$original_map)
provenance <- release_read_csv(paths$script_provenance)

release_log("  regenerated panel manifest: ", nrow(regenerated), " panels")

# --------------------------------------------------------------------------------------
# panel -> final revision script, from SCRIPT_PROVENANCE.csv
# --------------------------------------------------------------------------------------
# corrected_source_script_reference holds entries like "Figure3 C; Supplementary_proteomics B1"
# and "Supplementary_proteomics A/B2/B3". Parsed rather than transcribed so the mapping
# always follows the frozen provenance record.

parse_panel_refs <- function(reference) {
  reference <- trimws(as.character(reference))
  if (is.na(reference) || !nzchar(reference)) return(character(0))
  out <- character(0)
  for (chunk in trimws(strsplit(reference, ";", fixed = TRUE)[[1]])) {
    if (!nzchar(chunk)) next
    m <- regmatches(chunk, regexec("^(\\S+)\\s+(.*)$", chunk))[[1]]
    if (length(m) != 3L) next
    figure <- m[[2]]
    for (panel in trimws(strsplit(m[[3]], "/", fixed = TRUE)[[1]])) {
      if (nzchar(panel)) out <- c(out, paste(figure, panel, sep = "|"))
    }
  }
  unique(out)
}

# Script names are given relative to 06_manuscript_figure_revision/, so the 2026-08-27
# scripts read as before (e.g. 17_final_Fig3H_Fig3I.R) and the per-compartment rebuild
# copies as within_compartment_20261125/<script>. The generations are mapped separately:
# a rebuilt panel's script comes only from the current rebuild's rows, and the superseded
# 2026-11-24 rebuild maps no panel.
snapshot_name <- function(p) sub("^06_manuscript_figure_revision/", "", as.character(p))
is_rebuild_row <- as.character(provenance$generation_pass) == REBUILD_PASS
is_superseded_rebuild <- startsWith(as.character(provenance$generation_pass),
                                    "within_compartment_rebuild") & !is_rebuild_row
panel_to_script <- list()
panel_to_rebuild_script <- list()
for (i in seq_len(nrow(provenance))) {
  if (is_superseded_rebuild[i]) next
  refs <- parse_panel_refs(provenance$corrected_source_script_reference[i])
  for (key in refs) {
    if (is_rebuild_row[i]) {
      panel_to_rebuild_script[[key]] <- c(panel_to_rebuild_script[[key]], i)
    } else {
      panel_to_script[[key]] <- c(panel_to_script[[key]],
                                  snapshot_name(provenance$repository_snapshot_path[i]))
    }
  }
}
release_log("  SCRIPT_PROVENANCE.csv maps ", length(panel_to_script),
            " panels to a repository snapshot script and ", length(panel_to_rebuild_script),
            " rebuilt panels to a per-compartment rebuild script")

# Every frozen script must still match the executed copy, otherwise the figure provenance
# recorded in the release would be describing code that is no longer in the repository.
if (!all(toupper(as.character(provenance$sha256_identical)) == "TRUE")) {
  stop("SCRIPT_PROVENANCE.csv reports an executed/repository script hash mismatch; the ",
       "figure provenance cannot be published as exact.", call. = FALSE)
}
repo_snapshot_ok <- vapply(seq_len(nrow(provenance)), function(i) {
  p <- file.path(REPO_ROOT, as.character(provenance$repository_snapshot_path[i]))
  file.exists(p) && identical(release_sha256(p),
                              tolower(as.character(provenance$repository_sha256[i])))
}, logical(1))
if (!all(repo_snapshot_ok)) {
  stop("Frozen revision script(s) no longer hash to the value recorded in ",
       "SCRIPT_PROVENANCE.csv:\n",
       paste("  -", provenance$repository_snapshot_path[!repo_snapshot_ok], collapse = "\n"),
       call. = FALSE)
}
release_log("  all ", nrow(provenance), " frozen revision scripts hash as recorded")

# The rebuild scripts on the shared drive must still be the ones that were snapshotted: a
# later edit there would leave the repository copy describing code that no longer exists.
rebuild_rows <- which(is_rebuild_row)
if (!length(rebuild_rows)) {
  stop("SCRIPT_PROVENANCE.csv has no ", REBUILD_PASS, " rows.", call. = FALSE)
}
rebuild_exec_ok <- vapply(rebuild_rows, function(i) {
  p <- file.path(DATA_ROOT, as.character(provenance$executed_script_path[i]))
  file.exists(p) && identical(release_sha256(p),
                              tolower(as.character(provenance$executed_sha256[i])))
}, logical(1))
if (!all(rebuild_exec_ok)) {
  stop("Per-compartment rebuild script(s) on the shared drive changed since they were ",
       "snapshotted:\n", paste("  -", provenance$executed_script_path[rebuild_rows[!rebuild_exec_ok]],
                               collapse = "\n"), call. = FALSE)
}
release_log("  all ", length(rebuild_rows), " per-compartment rebuild scripts on the shared ",
            "drive hash as snapshotted")

# --------------------------------------------------------------------------------------
# analysis type / statistical unit / primary-secondary per panel
# --------------------------------------------------------------------------------------

PANEL_ANALYSIS <- c(
  "Figure3|A" = "experimental workflow schematic (artwork; no data)",
  "Figure3|B" = "representative immunofluorescence (microscopy; no proteomics data)",
  "Figure3|C" = "PCA of animal-level protein abundance (PC1/PC2 and PC2/PC3)",
  "Figure3|D" = "volcano of animal-level differential proteins, mCherry paired_veh vs unpaired_veh",
  "Figure3|E" = "rank-abundance curves by sample class and condition (descriptive)",
  "Figure3|F" = "GSEA GO-BP bubble grid, paired_veh vs unpaired_veh across sample classes",
  "Figure3|G" = "GSEA lollipop, GO-BP and GO-CC (ranked by the moderated t), mCherry paired_veh vs unpaired_veh",
  "Figure3|H" = "scatter of the learning effect size against the CNO effect size (standardised abundance differences, SD units), mCherry (descriptive)",
  "Figure3|I" = "GSEA NES heatmap, learning and paired-CNO contrasts",
  "Supplementary_proteomics|A" = "protein and peptide identification counts per acquisition (technical QC)",
  "Supplementary_proteomics|B1" = "PCA of animal-level protein abundance coloured by experimental group",
  "Supplementary_proteomics|B2" = "PCA of acquisition-level protein abundance coloured by hemisphere (technical QC)",
  "Supplementary_proteomics|B3" = "PCA of acquisition-level protein abundance coloured by collection plate (technical QC)",
  "Supplementary_proteomics|C" = "EWCE cell-type identity heatmap (Baseline analysis)",
  "Supplementary_proteomics|D" = "rank-abundance curves, neuron and neuropil at unpaired_veh (descriptive)",
  "Supplementary_proteomics|E" = "cross-compartment GSEA bubble plot from a secondary within-animal paired model",
  "Supplementary_proteomics|F" = "GSEA GO-BP bubble grid, paired_cno vs paired_veh across sample classes"
)

PANEL_PRIMARY <- c(
  "Figure3|A" = "not_applicable", "Figure3|B" = "not_applicable",
  "Figure3|C" = "primary", "Figure3|D" = "primary", "Figure3|E" = "primary",
  "Figure3|F" = "primary", "Figure3|G" = "primary", "Figure3|H" = "primary",
  "Figure3|I" = "primary",
  "Supplementary_proteomics|A" = "technical_qc",
  "Supplementary_proteomics|B1" = "primary",
  "Supplementary_proteomics|B2" = "technical_qc",
  "Supplementary_proteomics|B3" = "technical_qc",
  "Supplementary_proteomics|C" = "secondary",
  "Supplementary_proteomics|D" = "primary",
  "Supplementary_proteomics|E" = "secondary",
  "Supplementary_proteomics|F" = "primary"
)

PANEL_RELEASE_TABLE <- c(
  "Figure3|C" = "processed_data/protein_abundance_animal_level.tsv.gz",
  "Figure3|D" = "differential_analysis/primary_differential_proteins.tsv.gz",
  "Figure3|E" = "processed_data/protein_abundance_animal_level.tsv.gz",
  "Figure3|F" = "enrichment/primary_GSEA_GO_BP.tsv.gz",
  "Figure3|G" = "enrichment/primary_GSEA_GO_BP.tsv.gz;enrichment/primary_GSEA_GO_CC.tsv.gz",
  "Figure3|H" = "differential_analysis/primary_differential_proteins.tsv.gz",
  "Figure3|I" = "enrichment/primary_GSEA_GO_BP.tsv.gz",
  "Supplementary_proteomics|A" = "(not in this release: acquisition QC counts, see notes)",
  "Supplementary_proteomics|B1" = "processed_data/protein_abundance_animal_level.tsv.gz",
  "Supplementary_proteomics|B2" = "processed_data/protein_abundance_measurement_level.tsv.gz",
  "Supplementary_proteomics|B3" = "processed_data/protein_abundance_measurement_level.tsv.gz",
  "Supplementary_proteomics|C" = "enrichment/primary_EWCE.tsv.gz",
  "Supplementary_proteomics|D" = "processed_data/protein_abundance_animal_level.tsv.gz",
  "Supplementary_proteomics|E" = "(not in this release: secondary paired model, see notes)",
  "Supplementary_proteomics|F" = "enrichment/primary_GSEA_GO_BP.tsv.gz"
)

PANEL_ROW_FILTER <- c(
  "Figure3|C" = "all 5349 proteins x 48 animal-level units",
  "Figure3|D" = "canonical_comparison == mcherry_paired_veh_over_mcherry_unpaired_veh AND id_mapping_status == mapped",
  "Figure3|E" = "all 5349 proteins; means over the 3 animals per sample_class x condition",
  "Figure3|F" = "contrast_family == paired_veh_vs_unpaired_veh AND analysis_role == canonical AND ontology == GO_BP",
  # the selection rule is appended from the panel's own source data below
  "Figure3|G" = "canonical_comparison == mcherry_paired_veh_over_mcherry_unpaired_veh AND analysis_role == canonical AND ontology in {GO_BP, GO_CC}",
  "Figure3|H" = "canonical_comparison in {mcherry_paired_veh_over_mcherry_unpaired_veh, mcherry_paired_cno_over_mcherry_paired_veh} AND id_mapping_status == mapped AND significant_fdr_0_05 == TRUE in mcherry_paired_veh_over_mcherry_unpaired_veh",
  "Figure3|I" = "contrast_family in {paired_veh_vs_unpaired_veh, paired_cno_vs_paired_veh} AND analysis_role == canonical AND adjusted_p_value <= 0.05 in BOTH within a sample class",
  "Supplementary_proteomics|B1" = "all 5349 proteins x 48 animal-level units",
  "Supplementary_proteomics|B2" = "all 5349 proteins x 96 acquisitions; coloured by sample_metadata hemisphere",
  "Supplementary_proteomics|B3" = "all 5349 proteins x 96 acquisitions; coloured by sample_metadata collection_plate",
  "Supplementary_proteomics|C" = "ewce_analysis_type == Baseline AND is_primary_setting == TRUE",
  "Supplementary_proteomics|D" = "sample_class in {neuron, neuropil} AND condition == unpaired_veh",
  "Supplementary_proteomics|F" = "contrast_family == paired_cno_vs_paired_veh AND analysis_role == canonical AND ontology == GO_BP"
)

# Notes for the six rebuilt panels (D, F, G, H, I, Supp F) are computed from the rebuild's
# files below and added to this map; the pooled-run notes they replace are not kept.
PANEL_NOTE <- c(
  "Supplementary_proteomics|A" = paste("Acquisition-level identification depth. The unit is",
                                       "the acquisition, which is correct for technical QC",
                                       "and is NOT a biological n. Its source workbook lives",
                                       "outside the clusterProfiler tree and is cited in the",
                                       "map rather than re-derived."),
  "Supplementary_proteomics|B2" = paste("Retained at acquisition level by necessity:",
                                        "hemisphere does not exist in the 48-row animal-level",
                                        "metadata. Hemisphere explains 0.11% of PC1 and 0.02%",
                                        "of PC2, neither significant."),
  "Supplementary_proteomics|B3" = paste("Retained at acquisition level by necessity:",
                                        "collection plate does not exist in the animal-level",
                                        "metadata. Collection plate explains 0.01% of PC1 and",
                                        "1.5% of PC2, neither significant. Label reads",
                                        "'Collection plate'; the word 'batch' must not be used."),
  "Supplementary_proteomics|E" = paste("SECONDARY analysis, not one of the 12 primary",
                                       "contrasts and not a treatment effect. Both arms of",
                                       "every cross-compartment comparison come from the SAME",
                                       "three animals, so it was fitted as a within-animal",
                                       "paired model (limma ~ 0 + sample_class + AnimalID;",
                                       "residual df 6). n = 3 animals: low power. An earlier",
                                       "note blocking this panel was superseded."),
  "Supplementary_proteomics|C" = paste("EWCE Baseline analysis over per-condition abundant",
                                       "protein lists; descriptive cell-type identity, not a",
                                       "treatment contrast. 12 of the legacy 14 cell types",
                                       "survive animal-level correction."),
  "Figure3|E" = paste("Descriptive; no statistical test. Numerically indistinguishable from",
                      "the hemisphere-level original (Spearman 0.99999, Pearson on group",
                      "means 1.000), because in a balanced design the mean of animal means",
                      "equals the mean of all hemispheres."),
  "Supplementary_proteomics|D" = paste("Descriptive; no statistical test. Same equivalence as",
                                        "Fig. 2 E.")
)

# --------------------------------------------------------------------------------------
# panels rebuilt on the per-compartment statistics (2026-11-24; again 2026-11-25)
# --------------------------------------------------------------------------------------
# Source data pinned by hash: a later regeneration of the rebuild folder must fail here
# rather than be published unnoticed. A panel may have more than one file (G also publishes
# the terms skipped as redundant, Supp. F the cells left blank because fgsea could not test
# the set); paths and hashes are then ";"-separated, in order. D and H do not use the GO
# enrichment, so their source data are byte-identical to the 2026-11-24 rebuild's.

REBUILT_PANELS <- data.frame(
  key = c("Figure3|D", "Figure3|F", "Figure3|G", "Figure3|H", "Figure3|I",
          "Supplementary_proteomics|F"),
  source_data_path = c(
    "figure_panels/Figure3/Fig3D_mCherry_paired_vs_unpaired_volcano_animal_level_source_data.csv",
    "figure_panels/Figure3/Fig3F_memory_engram_GSEA_animal_level_source_data.csv",
    paste(file.path(REBUILD_STYLE, "Fig3G_source_data.csv"),
          file.path(REBUILD_STYLE, "Fig3G_redundant_terms_skipped.csv"), sep = ";"),
    file.path(REBUILD_STYLE, "Fig3H_source_data.csv"),
    file.path(REBUILD_STYLE, "Fig3I_source_data.csv"),
    paste("figure_panels/Supplementary_proteomics/SuppF_inhibition_GSEA_animal_level_source_data.csv",
          "full_regenerated/enrichment/SuppF_inhibition_GSEA_animal_level_blank_cells.csv",
          sep = ";")),
  source_data_sha256 = c(
    "8033fd31d77b18b57f9bcbf06c1f9d53e383d5e8e9f3537f31bbe2a49bb71b7c",
    "29c5888b43bd11fe322e4397d5511426b6aaa6f44294fc17c7a9965fdb3ed98f",
    paste("e263647698b563f027367304116c6ebdfd49baa87766aabe13527880d1f08758",
          "3dd2c0e0ac90975581d8fff01de90c9eacf6a815fb2f8ea99b9e69fcd8e7b5b5", sep = ";"),
    "f56893b79a6a0e26d0776e3c2e4b5eb4ebd3578e6e7f2fa726155dff6d683bd5",
    "72ac22401f8340ded04f4a01da1a0d649325b1e5931b05ffc6c81d37328d3286",
    paste("5bbd23918ff8c5664d68609e7d925f8bf9d6e24a5454a96314f13273774f96ca",
          "7d737ae3cee2e99cafcf8a55f5b8eb10fcd6acae5ef19f713548bf50c4f59ada", sep = ";")),
  status = c("REBUILT_WITHIN_COMPARTMENT", "REBUILT_WITHIN_COMPARTMENT",
             "REBUILT_WITHIN_COMPARTMENT_GSEA_REPLACES_ORA",
             "REBUILT_WITHIN_COMPARTMENT_DESCRIPTIVE",
             "REBUILT_WITHIN_COMPARTMENT_WITH_SHARED_REFERENCE_CAVEAT",
             "REBUILT_WITHIN_COMPARTMENT"),
  canonical_input = c(
    file.path(STATS_REL, "mapped", "forward", "mcherry_paired_veh_vs_mcherry_unpaired_veh.csv"),
    file.path(ENRICH_REL, "GO_BP", "per_comparison", "*_paired_veh_over_*_unpaired_veh",
              "GSEA_GO_BP.csv"),
    file.path(ENRICH_REL, "{GO_BP,GO_CC}", "per_comparison",
              "mcherry_paired_veh_over_mcherry_unpaired_veh", "GSEA_GO_{BP,CC}.csv"),
    paste(file.path(STATS_REL, "mapped", "forward",
                    "mcherry_paired_veh_vs_mcherry_unpaired_veh.csv"),
          "and mcherry_paired_cno_vs_mcherry_paired_veh.csv"),
    paste(file.path(ENRICH_REL, "GO_BP", "per_comparison"),
          "(learning and paired-CNO contrasts)"),
    file.path(ENRICH_REL, "GO_BP", "per_comparison", "*_paired_cno_over_*_paired_veh",
              "GSEA_GO_BP.csv")),
  # Decided 2026-11-25, relative to the pooled-model release.
  biological_interpretation_changed = c("NO", "YES", "YES", "NO", "YES", "NO"),
  # Title and unit stated here rather than copied from the 2026-08-27 manifest, whose text
  # described the pooled-model panels (and called the H axes "moderated limma log2FC").
  # Axis wording as decided 2026-11-25: D "Effect size (SD)"; H "Learning effect (SD)" against
  # "CNO vs VEH effect, paired (SD)".
  panel_title = c(
    "Volcano: mCherry paired vs unpaired, effect size (SD) against adjusted P",
    "GSEA bubble grid: memory engram, paired vs unpaired",
    "GSEA lollipop, GO-BP and GO-CC: mCherry paired vs unpaired, rows by adjusted P",
    "Learning effect (SD) vs CNO vs VEH effect, paired (SD), mCherry",
    "NES heatmap: learning vs CNO pathway profiles",
    "GSEA bubble grid: inhibition, paired CNO vs paired VEH, four sample classes"),
  statistical_unit = c(
    paste("animal (n = 3 per group); effect size = standardized abundance difference (SD",
          "units); moderated t and BH FDR from the per-sample-class limma fit"),
    "animal; GSEA ranked by the moderated t from the per-sample-class limma fit",
    "animal (n = 3 per group); GSEA ranked by the moderated t from the per-sample-class limma fit",
    paste("animal on both axes; standardized abundance differences (SD units) from the",
          "per-sample-class limma fit"),
    "animal; GSEA ranked by the moderated t from the per-sample-class limma fit",
    "animal; GSEA ranked by the moderated t from the per-sample-class limma fit"),
  stringsAsFactors = FALSE, check.names = FALSE
)

rebuilt_files <- function(key) {
  i <- match(key, REBUILT_PANELS$key)
  list(rel = strsplit(REBUILT_PANELS$source_data_path[i], ";", fixed = TRUE)[[1]],
       sha = strsplit(REBUILT_PANELS$source_data_sha256[i], ";", fixed = TRUE)[[1]])
}
rebuilt_csv <- function(key, n = 1L) {
  release_read_csv(file.path(REBUILD_ROOT, rebuilt_files(key)$rel[[n]]))
}
pin_bad <- character(0)
for (key in REBUILT_PANELS$key) {
  f <- rebuilt_files(key)
  observed <- vapply(file.path(REBUILD_ROOT, f$rel), release_sha256, character(1),
                     USE.NAMES = FALSE)
  if (length(f$rel) != length(f$sha) || !identical(observed, f$sha)) pin_bad <- c(pin_bad, key)
}
if (length(pin_bad)) {
  stop("Rebuilt panel source data missing or no longer matching the pinned hash: ",
       paste(pin_bad, collapse = ", "), call. = FALSE)
}
missing_rebuild_script <- setdiff(REBUILT_PANELS$key, names(panel_to_rebuild_script))
if (length(missing_rebuild_script)) {
  stop("No per-compartment rebuild script recorded in SCRIPT_PROVENANCE.csv for: ",
       paste(missing_rebuild_script, collapse = ", "), call. = FALSE)
}
release_log("  ", nrow(REBUILT_PANELS), " rebuilt panels: source data match their pinned hashes")

# Numbers for the rebuilt panels, read from their files (nothing is recomputed).
contrasts <- release_read_tsv_plain(paths$contrast_manifest)
enrich_index <- release_read_csv(paths$enrichment_index)
fdr_terms_by_class <- function(family) {
  cm <- contrasts[contrasts$contrast_family == family, , drop = FALSE]
  cm <- cm[order(match(cm$sample_class, sample_classes)), , drop = FALSE]
  n <- as.integer(enrich_index$go_gsea_fdr_term_count[
    match(cm$canonical_comparison, enrich_index$canonical_comparison)])
  if (anyNA(n)) stop("GO-BP enrichment index lacks a ", family, " comparison.", call. = FALSE)
  stats::setNames(n, cm$sample_class)
}
by_class_text <- function(x) {
  paste0(release_sample_class_display(names(x)), " ", unname(x), collapse = ", ")
}

# Gene sets fgsea could not test in the blocks a panel draws on, and the permutation settings
# behind every GSEA count.
untested <- release_fgsea_untested(contrasts$canonical_comparison, DATA_ROOT)
permutations <- release_enrichment_permutations(DATA_ROOT)
untested_sentence <- function(families, ontologies = "GO_BP", classes = NULL) {
  blocks <- release_untested_blocks_text(untested, contrasts, families, ontologies, classes)
  if (!nzchar(blocks)) {
    return(paste("fgsea tested every gene set in the blocks this panel uses.",
                 permutations$sentence))
  }
  paste0("fgsea could not test some gene sets in the blocks this panel uses, so they are",
         " absent from its GSEA results rather than non-significant (", blocks, " sets). ",
         permutations$sentence)
}
fmt_r <- function(x) sprintf("%+.2f", as.numeric(x))
audit_value <- function(audit, item) {
  v <- audit$value[audit$item == item]
  if (length(v) != 1L) stop("Audit file has no unique item ", item, call. = FALSE)
  v
}

# Fig. 2 D: volcano of the mapped learning contrast.
d_src <- rebuilt_csv("Figure3|D")
d_up <- sum(d_src$Significance == "up"); d_down <- sum(d_src$Significance == "down")
d_labels <- d_src$mapped_gene_symbol[as.logical(d_src$labelled_in_panel)]
d_key <- paste0(nrow(d_src), " mapped proteins tested, ", d_up + d_down, " FDR<0.05 (",
                d_up, " up / ", d_down, " down); labels ", paste(d_labels, collapse = ", "))

# Fig. 2 F and Supp. Fig. 4 F: GO-BP bubble grids; FDR counts as recorded by the run.
# Classes are listed in the panel's own column order (order of first appearance).
f_src <- rebuilt_csv("Figure3|F")
f_cols <- unique(as.character(f_src$sample_class))
f_counts <- fdr_terms_by_class("paired_veh_vs_unpaired_veh")[f_cols]
f_plot_sig <- tapply(as.logical(f_src$fdr_significant), f_src$sample_class, sum)[f_cols]
f_key <- paste0("FDR<0.05 GO-BP terms per column: ", by_class_text(f_counts), "; ",
                length(unique(f_src$ID)), " distinct terms plotted")
f_note <- paste0("GSEA GO-BP ranked by the moderated t. Of the ", length(unique(f_src$ID)),
                 " plotted terms, FDR < 0.05 per column: ", by_class_text(f_plot_sig), ". ",
                 untested_sentence("paired_veh_vs_unpaired_veh", "GO_BP", f_cols), " ",
                 RELEASE_PLATE_SENTENCE)
s_counts <- fdr_terms_by_class("paired_cno_vs_paired_veh")
s_src <- rebuilt_csv("Supplementary_proteomics|F")
s_shown <- unique(as.character(s_src$sample_class))
s_key <- paste0("FDR<0.05 GO-BP terms per column: ", by_class_text(s_counts[s_shown]),
                if (length(setdiff(names(s_counts), s_shown)))
                  paste0(" (not shown: ", by_class_text(s_counts[setdiff(names(s_counts), s_shown)]),
                         ")") else "",
                "; ", length(unique(s_src$ID)), " distinct terms plotted")
s_blank <- rebuilt_csv("Supplementary_proteomics|F", 2L)
s_note <- paste0("GSEA GO-BP ranked by the moderated t, paired CNO vs paired VEH. ",
                 if (nrow(s_blank)) paste0(
                   nrow(s_blank), " cell", if (nrow(s_blank) > 1L) "s are" else " is",
                   " left blank because fgsea could not test the set in that block (",
                   paste0(release_sample_class_display(s_blank$sample_class), ": ",
                          s_blank$Description, collapse = "; "),
                   "); they are listed in the panel's second source-data file. ") else "",
                 untested_sentence("paired_cno_vs_paired_veh", "GO_BP", s_shown))

# Fig. 2 G: GSEA lollipop of GO-BP and GO-CC terms; it replaced an ORA lollipop.
g_src <- rebuilt_csv("Figure3|G", 1L)
g_skip <- rebuilt_csv("Figure3|G", 2L)
g_ont <- unique(g_src[c("ontology", "n_terms_tested", "n_FDR_significant")])
g_ont$plotted <- as.integer(table(g_src$ontology)[g_ont$ontology])
g_ont$skipped <- vapply(g_ont$ontology, function(o) sum(g_skip$ontology == o), integer(1))
g_rule <- unique(as.character(g_src$selection_rule))
if (nrow(g_ont) != 2L || length(g_rule) != 1L) {
  stop("Fig. 2 G source data no longer holds one GO-BP and one GO-CC block.", call. = FALSE)
}
ora_up <- release_read_csv(paths$ora_learning_up)
ora_up_cc <- release_read_csv(paths$ora_learning_up_cc)
ora_text <- function(d, label) {
  paste0(sum(d$p.adjust < 0.05), " of ", format(nrow(d), big.mark = ","), " ", label,
         " terms at FDR < 0.05 (minimum adjusted P ", signif(min(d$p.adjust), 2), ")")
}
g_key <- paste0("GSEA (moderated t): ",
                paste0("GO-", g_ont$ontology, " ", g_ont$n_FDR_significant, " of ",
                       g_ont$n_terms_tested, " terms at FDR<0.05", collapse = "; "),
                "; plotted ", paste0(g_ont$plotted, " GO-", g_ont$ontology, collapse = " + "),
                " (", sum(g_ont$skipped), " skipped as redundant)")
g_note <- paste0(
  "GSEA ranked by the moderated t. Selection: ", g_rule, ". Rows are ordered by adjusted P",
  " within each ontology. The skipped terms are",
  " published as a second file. Replaces the pooled-run ORA lollipop. The panel shows",
  " GSEA, not over-representation; for comparison, ORA over the FDR-significant proteins",
  " higher in paired_veh finds ", ora_text(ora_up, "GO-BP"), " and ",
  ora_text(ora_up_cc, "GO-CC"), ". Term counts near q = 0.05 move by a few terms with",
  " the GSEA permutation seed. ",
  untested_sentence("paired_veh_vs_unpaired_veh", c("GO_BP", "GO_CC"), "mcherry"), " ",
  RELEASE_PLATE_SENTENCE)
PANEL_ROW_FILTER[["Figure3|G"]] <- paste0(PANEL_ROW_FILTER[["Figure3|G"]], "; ", g_rule)

# Fig. 2 H: descriptive correlation of the learning and CNO effects.
h_audit <- release_read_csv(paths$fig_h_audit)
h_key <- paste0("r = ", fmt_r(audit_value(h_audit, "corrected_r")),
                " on the learning signature (n = ", audit_value(h_audit, "n_proteins"),
                "); r = ", fmt_r(audit_value(h_audit, "r_all_proteins_unselected")),
                " unselected; same-side control r = ",
                fmt_r(audit_value(h_audit, "control_r_learning_signature")))
h_note <- paste0(
  "DESCRIPTIVE ONLY. Pearson r = ", fmt_r(audit_value(h_audit, "corrected_r")), " across the ",
  audit_value(h_audit, "n_proteins"), " proteins at FDR < 0.05 in the learning contrast.",
  " The two contrasts share paired_veh on OPPOSITE sides, so a negative correlation is",
  " expected by construction; the same-side control gives r = ",
  fmt_r(audit_value(h_audit, "control_r_learning_signature")), " on the same protein set.",
  " Conditioning x on learning FDR additionally induces regression to the mean on y. This",
  " panel is not evidence that inhibition reverses learning.")

# Fig. 2 I: NES heatmap of terms significant in both rows. Opposition is counted from the
# source data and must agree with the rebuild's own audit.
i_src <- rebuilt_csv("Figure3|I")
i_audit <- release_read_csv(paths$fig_i_audit)
i_rows <- unique(as.character(i_src$Row))
if (length(i_rows) != 2L) stop("Fig. 2 I source data does not hold two rows.", call. = FALSE)
i_key_cols <- paste(i_src$sample_class, i_src$Description, sep = "\r")
i_nes <- tapply(as.numeric(i_src$NES), list(i_key_cols, as.character(i_src$Row)), identity)
if (!is.numeric(i_nes) || anyNA(i_nes)) {
  stop("Fig. 2 I source data does not give one NES per term and row.", call. = FALSE)
}
i_terms <- data.frame(term = rownames(i_nes),
                      sample_class = sub("\r.*$", "", rownames(i_nes)),
                      opposed = sign(as.numeric(i_nes[, 1])) != sign(as.numeric(i_nes[, 2])),
                      stringsAsFactors = FALSE)
i_per_class <- table(factor(i_terms$sample_class,
                            levels = unique(as.character(i_src$sample_class))))
i_conc_classes <- unique(i_terms$sample_class[!i_terms$opposed])
if (nrow(i_terms) != as.integer(audit_value(i_audit, "n_overlapping_terms")) ||
    sum(i_terms$opposed) != as.integer(audit_value(i_audit, "n_opposed"))) {
  stop("Fig. 2 I source data and its interpretation audit disagree on the term counts.",
       call. = FALSE)
}
i_key <- paste0(nrow(i_terms), " overlapping FDR terms: ",
                by_class_text(stats::setNames(as.integer(i_per_class), names(i_per_class))),
                "; ", sum(i_terms$opposed), " opposed, ", sum(!i_terms$opposed), " concordant",
                if (length(i_conc_classes))
                  paste0(" (", paste(release_sample_class_display(i_conc_classes),
                                     collapse = ", "), ")") else "")
i_note <- paste0(
  "Row 1 (paired_veh / unpaired_veh) and row 2 (paired_cno / paired_veh) share paired_veh on",
  " OPPOSITE sides. Of the ", nrow(i_terms), " terms FDR-significant in both rows within a",
  " sample class, ", sum(i_terms$opposed), " are opposed in sign, the expected structural",
  " consequence of the design",
  if (sum(!i_terms$opposed))
    paste0("; the ", sum(!i_terms$opposed), " concordant term(s) are in ",
           paste(release_sample_class_display(i_conc_classes), collapse = ", ")) else "",
  ". The opposition cannot be read as proof that CNO reverses a pairing-associated programme. ",
  untested_sentence(c("paired_veh_vs_unpaired_veh", "paired_cno_vs_paired_veh"), "GO_BP",
                    unique(as.character(i_src$sample_class))))

REBUILT_PANELS$key_numbers <- c(d_key, f_key, g_key, h_key, i_key, s_key)
PANEL_NOTE[["Figure3|F"]] <- f_note
PANEL_NOTE[["Figure3|G"]] <- g_note
PANEL_NOTE[["Figure3|H"]] <- h_note
PANEL_NOTE[["Figure3|I"]] <- i_note
PANEL_NOTE[["Supplementary_proteomics|F"]] <- s_note
release_log("  rebuilt-panel numbers and notes read from the rebuild's source-data and ",
            "audit files")

# --------------------------------------------------------------------------------------
# assemble the map
# --------------------------------------------------------------------------------------

orig_key <- paste(original_map$figure, trimws(original_map$panel), sep = "|")

rows <- lapply(seq_len(nrow(regenerated)), function(i) {
  r <- regenerated[i, , drop = FALSE]
  figure <- trimws(as.character(r$figure))
  panel <- trimws(as.character(r$panel))
  key <- paste(figure, panel, sep = "|")

  rb <- match(key, REBUILT_PANELS$key)
  rebuilt <- !is.na(rb)
  rebuild_prov <- if (rebuilt) provenance[panel_to_rebuild_script[[key]], , drop = FALSE] else NULL

  final_script <- if (rebuilt) snapshot_name(rebuild_prov$repository_snapshot_path) else
    panel_to_script[[key]]
  final_script <- if (is.null(final_script)) NA_character_ else
    paste(unique(final_script), collapse = ";")

  oi <- match(key, orig_key)
  panel_root <- if (rebuilt) REBUILD_ROOT else REVISION_ROOT
  source_data_rel <- if (rebuilt) rebuilt_files(key)$rel else
    trimws(as.character(r$source_data_path))
  has_source_data <- nzchar(source_data_rel[[1]]) && !source_data_rel[[1]] %in% c("n/a", "NA")
  source_data_abs <- if (has_source_data) file.path(panel_root, source_data_rel) else NA_character_
  # Destination names are flattened to <panel>_source_data.<ext>, named in the manuscript
  # numbering (Fig2D, SuppFig4F). The canonical layout is nested three levels deep with long
  # descriptive filenames, and reproducing it pushed paths past the Windows 260-character
  # limit. The original filename is kept in a column so the correspondence stays exact. A
  # panel's further files keep their own name after the panel slug, with the internal slug
  # replaced (Fig3G_redundant_terms_skipped.csv -> Fig2G_redundant_terms_skipped.csv); all
  # are ";"-separated in order.
  panel_label <- paste(ifelse(figure == "Figure3", "Fig. 2", "Supp. Fig. 4"), panel)
  panel_slug <- gsub("[^A-Za-z0-9]", "", sub("^Supp[.] Fig[.]", "SuppFig", panel_label))
  internal_slug <- paste0(ifelse(figure == "Figure3", "Fig3", "Supp"), panel)
  dest_rel <- if (has_source_data) {
    first <- paste0(panel_slug, "_source_data.", tolower(tools::file_ext(source_data_rel[[1]])))
    extra <- if (length(source_data_rel) > 1L) {
      paste0(panel_slug, "_", sub(paste0("^", internal_slug, "_"), "",
                                  basename(source_data_rel[-1])))
    } else character(0)
    file.path("editor_source_data", "figure_source_data", c(first, extra))
  } else NA_character_
  join <- function(x) paste(x, collapse = ";")

  data.frame(
    figure = figure,
    panel = panel,
    panel_label = panel_label,
    panel_title = if (rebuilt) REBUILT_PANELS$panel_title[rb] else
      trimws(as.character(r$panel_title)),
    analysis = unname(ifelse(key %in% names(PANEL_ANALYSIS), PANEL_ANALYSIS[[key]],
                             NA_character_)),
    source_table = unname(ifelse(key %in% names(PANEL_RELEASE_TABLE),
                                 PANEL_RELEASE_TABLE[[key]], NA_character_)),
    source_rows_filter = unname(ifelse(key %in% names(PANEL_ROW_FILTER),
                                       PANEL_ROW_FILTER[[key]], NA_character_)),
    statistical_unit = if (rebuilt) REBUILT_PANELS$statistical_unit[rb] else
      trimws(as.character(r$new_sampling_unit)),
    original_statistical_unit = trimws(as.character(r$old_sampling_unit)),
    primary_or_secondary = unname(ifelse(key %in% names(PANEL_PRIMARY),
                                         PANEL_PRIMARY[[key]], NA_character_)),
    canonical_source = if (rebuilt) REBUILT_PANELS$canonical_input[rb] else
      trimws(as.character(r$corrected_canonical_input)),
    generating_script = if (is.na(oi)) NA_character_ else
      trimws(as.character(original_map$original_generating_script[oi])),
    corrected_source_script = if (rebuilt) join(rebuild_prov$executed_script_path) else
      trimws(as.character(r$corrected_source_script)),
    final_revision_script = final_script,
    revision_status = if (rebuilt) REBUILT_PANELS$status[rb] else trimws(as.character(r$status)),
    displayed_statistics_changed = if (rebuilt) "YES" else
      trimws(as.character(r$displayed_statistics_changed)),
    biological_interpretation_changed = if (rebuilt)
      REBUILT_PANELS$biological_interpretation_changed[rb] else
        trimws(as.character(r$biological_interpretation_changed)),
    key_numbers_original = trimws(as.character(r$key_numbers_old)),
    key_numbers_revised = if (rebuilt) REBUILT_PANELS$key_numbers[rb] else
      trimws(as.character(r$key_numbers_new)),
    # The superseded pooled-model numbers, kept for the "what changed" comparison.
    key_numbers_pooled_model = if (rebuilt) trimws(as.character(r$key_numbers_new)) else
      NA_character_,
    panel_source_data_release_path = if (has_source_data) join(dest_rel) else NA_character_,
    panel_source_data_original_filename = if (has_source_data) join(basename(source_data_rel)) else
      NA_character_,
    panel_source_data_canonical_path = if (has_source_data) join(source_data_abs) else
      NA_character_,
    panel_source_data_sha256 = if (has_source_data && all(file.exists(source_data_abs)))
      join(vapply(source_data_abs, release_sha256, character(1), USE.NAMES = FALSE)) else
        NA_character_,
    interpretation_note = unname(ifelse(key %in% names(PANEL_NOTE), PANEL_NOTE[[key]],
                                        NA_character_)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
})
figure_map <- do.call(rbind, rows)
rownames(figure_map) <- NULL

REQUIRED_PANELS <- c(
  "Figure3|C", "Figure3|D", "Figure3|E", "Figure3|F", "Figure3|G", "Figure3|H", "Figure3|I",
  "Supplementary_proteomics|A", "Supplementary_proteomics|B1", "Supplementary_proteomics|B2",
  "Supplementary_proteomics|B3", "Supplementary_proteomics|C", "Supplementary_proteomics|D",
  "Supplementary_proteomics|E", "Supplementary_proteomics|F"
)
present <- paste(figure_map$figure, figure_map$panel, sep = "|")
missing_panels <- setdiff(REQUIRED_PANELS, present)
if (length(missing_panels)) {
  stop("figure_source_map is missing required panel(s): ",
       paste(missing_panels, collapse = ", "), call. = FALSE)
}
data_panels <- figure_map[figure_map$primary_or_secondary != "not_applicable", , drop = FALSE]
no_script <- data_panels$panel_label[is.na(data_panels$final_revision_script)]
if (length(no_script)) {
  stop("Data-bearing panel(s) without a final revision script: ",
       paste(no_script, collapse = ", "), call. = FALSE)
}
release_log("  figure map: ", nrow(figure_map), " panels, ", nrow(data_panels),
            " data-bearing")

# The map must carry exactly the pinned source data for every rebuilt panel.
map_key <- paste(figure_map$figure, figure_map$panel, sep = "|")
rb_idx <- match(map_key, REBUILT_PANELS$key)
pin_mismatch <- which(!is.na(rb_idx) &
  (is.na(figure_map$panel_source_data_sha256) |
     figure_map$panel_source_data_sha256 != REBUILT_PANELS$source_data_sha256[rb_idx]))
if (length(pin_mismatch) || sum(!is.na(rb_idx)) != nrow(REBUILT_PANELS)) {
  stop("figure_source_map does not carry the pinned rebuilt source data for: ",
       paste(figure_map$panel_label[pin_mismatch], collapse = ", "), call. = FALSE)
}

# --------------------------------------------------------------------------------------
# copy the per-panel source data verbatim
# --------------------------------------------------------------------------------------

copied <- 0L
missing_source_data <- character(0)
split_paths <- function(x) strsplit(x, ";", fixed = TRUE)[[1]]
for (i in seq_len(nrow(figure_map))) {
  if (is.na(figure_map$panel_source_data_canonical_path[i]) ||
      is.na(figure_map$panel_source_data_release_path[i])) next
  srcs <- split_paths(figure_map$panel_source_data_canonical_path[i])
  rels <- split_paths(figure_map$panel_source_data_release_path[i])
  if (!all(file.exists(srcs))) {
    missing_source_data <- c(missing_source_data,
                             paste0(figure_map$panel_label[i], ": ",
                                    paste(srcs[!file.exists(srcs)], collapse = ";")))
    figure_map$panel_source_data_release_path[i] <- NA_character_
    next
  }
  shas <- split_paths(figure_map$panel_source_data_sha256[i])
  for (k in seq_along(srcs)) {
    dest <- release_path(rels[[k]])
    if (!file.copy(srcs[[k]], dest, overwrite = TRUE, copy.date = FALSE)) {
      stop("Could not copy panel source data: ", srcs[[k]], call. = FALSE)
    }
    if (!identical(release_sha256(dest), shas[[k]])) {
      stop("Copied panel source data does not match its source hash: ", rels[[k]],
           call. = FALSE)
    }
    release_register(rels[[k]],
                     paste0(if (k == 1L) "source data" else
                              paste0("additional source data (", basename(srcs[[k]]), ")"),
                            " for panel ", figure_map$panel_label[i]),
                     srcs[[k]], shas[[k]], STAGE, tolower(tools::file_ext(srcs[[k]])))
    copied <- copied + 1L
  }
}
release_log("  copied ", copied, " panel source-data files byte-for-byte")
if (length(missing_source_data)) {
  release_log("  NOTE: no source-data file on disk for: ",
              paste(missing_source_data, collapse = " | "))
}

# --------------------------------------------------------------------------------------
# write
# --------------------------------------------------------------------------------------

sources <- c(paths$regenerated, paths$original_map, paths$script_provenance,
             paths$retained_report, paths$rebuild_record, paths$fig_h_audit,
             paths$fig_i_audit, paths$enrichment_index, paths$ora_learning_up,
             paths$ora_learning_up_cc)
hashes <- vapply(sources, release_sha256, character(1), USE.NAMES = FALSE)

w <- release_write_table(figure_map,
                         release_path("editor_source_data", "figure_source_map.tsv"))
release_register("editor_source_data/figure_source_map.tsv",
                 "per-panel provenance: analysis, source table, row filter, unit, scripts",
                 sources, hashes, STAGE, "tsv")
release_log("  wrote figure_source_map.tsv (", w$rows, "x", w$cols, ")")
release_log("stage 06 complete")
