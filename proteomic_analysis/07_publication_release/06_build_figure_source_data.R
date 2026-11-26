#!/usr/bin/env Rscript

# Publication release, stage 06 -- figure source data and the panel provenance map.
#
# Produces
#   provenance/figure_source_map.tsv          one row per manuscript panel
#   editor_source_data/figure_source_data/*.csv       reader-facing source data per panel
#                                                     (G and Supp. F have two files each)
#
# The published source data are reader-facing versions of the canonical panel files: plain
# headers (the workbook's display header wherever a column holds the same quantity, as pinned
# in RELEASE_FIGURE_CSV_HEADERS), the display labels of the workbook, no internal columns,
# values carried over as text. The canonical files are checked against pinned hashes and
# named in the map with that hash, but not copied into the release. Supp. Fig. 4 A takes
# only the per-acquisition counts of its source workbook (sheet 1); the workbook's other
# sheet is unrelated data.
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
release_require("digest", "readxl")

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
                                   create_dir = FALSE),
  # Release tables of the earlier stages, for the lookups and cross-checks of the
  # reader-facing source data.
  sample_metadata = release_path("metadata", "sample_metadata.tsv", create_dir = FALSE),
  animal_metadata = release_path("metadata", "animal_level_sample_metadata.tsv",
                                 create_dir = FALSE),
  feature_annotation = release_path("processed_data", "protein_feature_annotation.tsv.gz",
                                    create_dir = FALSE),
  differential = release_path("differential_analysis", "primary_differential_proteins.tsv.gz",
                              create_dir = FALSE),
  gsea_go_bp = release_path("enrichment", "primary_GSEA_GO_BP.tsv.gz", create_dir = FALSE),
  untested_sets = release_path("enrichment", "gsea_untested_gene_sets.tsv", create_dir = FALSE)
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
  "Figure3|H" = "scatter of the learning effect size against the CNO effect size (standardized abundance differences, SD units), mCherry (descriptive)",
  "Figure3|I" = "GSEA NES heatmap, learning and paired-CNO contrasts",
  "Supplementary_proteomics|A" = "protein and peptide identification counts per acquisition (technical QC)",
  "Supplementary_proteomics|B1" = "PCA of animal-level protein abundance colored by experimental group",
  "Supplementary_proteomics|B2" = "PCA of acquisition-level protein abundance colored by hemisphere (technical QC)",
  "Supplementary_proteomics|B3" = "PCA of acquisition-level protein abundance colored by collection plate (technical QC)",
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
  "Supplementary_proteomics|A" = "(acquisition QC counts: in the panel source data only)",
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
  # the rows per curve are appended from the panel's own source data below
  "Figure3|E" = "means over the 3 animals per sample_class x condition",
  "Figure3|F" = "contrast_family == paired_veh_vs_unpaired_veh AND analysis_role == canonical AND ontology == GO_BP",
  # the selection rule is appended from the panel's own source data below
  "Figure3|G" = "canonical_comparison == mcherry_paired_veh_over_mcherry_unpaired_veh AND analysis_role == canonical AND ontology in {GO_BP, GO_CC}",
  "Figure3|H" = "canonical_comparison in {mcherry_paired_veh_over_mcherry_unpaired_veh, mcherry_paired_cno_over_mcherry_paired_veh} AND id_mapping_status == mapped AND significant_fdr_0_05 == TRUE in mcherry_paired_veh_over_mcherry_unpaired_veh",
  "Figure3|I" = "contrast_family in {paired_veh_vs_unpaired_veh, paired_cno_vs_paired_veh} AND analysis_role == canonical AND adjusted_p_value <= 0.05 in BOTH within a sample class",
  "Supplementary_proteomics|B1" = "all 5349 proteins x 48 animal-level units",
  "Supplementary_proteomics|B2" = "all 5349 proteins x 96 acquisitions; colored by sample_metadata hemisphere",
  "Supplementary_proteomics|B3" = "all 5349 proteins x 96 acquisitions; colored by sample_metadata collection_plate",
  "Supplementary_proteomics|C" = "ewce_analysis_type == Baseline AND is_primary_setting == TRUE",
  "Supplementary_proteomics|D" = "sample_class in {neuron, neuropil} AND condition == unpaired_veh; means over the 3 animals",
  "Supplementary_proteomics|F" = "contrast_family == paired_cno_vs_paired_veh AND analysis_role == canonical AND ontology == GO_BP"
)

# The rank-abundance panels plot 2^(mean standardized abundance). The mean is in SD units, so
# the plotted value is a display transform, not an abundance.
RANK_PLOTTED_VALUE_NOTE <- paste(
  "The plotted value is 2 to the power of the mean standardized abundance (SD units): a",
  "display transform, not an abundance.")

# Notes for the six rebuilt panels (D, F, G, H, I, Supp F) are computed from the rebuild's
# files below and added to this map; the pooled-run notes they replace are not kept.
PANEL_NOTE <- c(
  # The sentence on corrected sample classes is added from the published file below.
  "Supplementary_proteomics|A" = paste("Acquisition-level identification depth. The unit is",
                                       "the acquisition, which is correct for technical QC",
                                       "and is NOT a biological n."),
  "Supplementary_proteomics|B2" = paste("Shown at acquisition level because hemisphere is",
                                        "averaged out of the 48 animal-level profiles.",
                                        "Hemisphere explains 0.11% of PC1 and 0.02% of PC2,",
                                        "neither significant."),
  "Supplementary_proteomics|B3" = paste("Acquisition-level technical check. Collection plate",
                                        "explains 0.01% of PC1 and 1.5% of PC2, neither",
                                        "significant.", release_plate_injection_sentence(),
                                        "Label reads 'Collection plate'; the word 'batch' must",
                                        "not be used."),
  "Supplementary_proteomics|E" = paste("SECONDARY analysis, not one of the 12 primary",
                                       "comparisons and not a treatment effect. Both arms of",
                                       "every cross-compartment comparison come from the SAME",
                                       "three animals, so it was fitted as a within-animal",
                                       "paired model (limma ~ 0 + sample_class + AnimalID;",
                                       "residual df 6). n = 3 animals: low power. An earlier",
                                       "note blocking this panel was superseded."),
  "Supplementary_proteomics|C" = paste("EWCE Baseline analysis over per-condition abundant",
                                       "protein lists; descriptive cell-type identity, not a",
                                       "treatment comparison. 12 of the legacy 14 cell types",
                                       "survive animal-level correction."),
  "Figure3|E" = paste("Descriptive; no statistical test.", RANK_PLOTTED_VALUE_NOTE,
                      "Numerically indistinguishable from",
                      "the hemisphere-level original (Spearman 0.99999, Pearson on group",
                      "means 1.000), because in a balanced design the mean of animal means",
                      "equals the mean of all hemispheres."),
  "Supplementary_proteomics|D" = paste("Descriptive; no statistical test.",
                                        RANK_PLOTTED_VALUE_NOTE, "Same equivalence as",
                                        "Fig. 2 E.")
)

# Reader-facing unit of the panels whose manifest unit is written in the revision's terms.
# The manifest's explanation after " - " is kept (the changelog shows it; table 2p does not).
PANEL_UNIT <- c(
  "Figure3|E" = "animal (left/right averaged), then mean of 3 animals",
  "Supplementary_proteomics|B2" = "acquisition (96 runs)",
  "Supplementary_proteomics|B3" = "acquisition (96 runs)",
  "Supplementary_proteomics|C" = "protein list of one compartment x condition (n = 3 animals)",
  "Supplementary_proteomics|D" = "animal (left/right averaged), then mean of 3 animals"
)
panel_unit <- function(key, manifest_unit) {
  if (!key %in% names(PANEL_UNIT)) return(manifest_unit)
  at <- regexpr(" - ", manifest_unit, fixed = TRUE)
  paste0(PANEL_UNIT[[key]], if (at > 0L) substring(manifest_unit, at) else "")
}

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
  # 2026-11-25, relative to the pooled-model release.
  biological_interpretation_changed = c("NO", "YES", "YES", "NO", "YES", "NO"),
  # Title and unit stated here rather than copied from the 2026-08-27 manifest, whose text
  # described the pooled-model panels (and called the H axes "moderated limma log2FC").
  # Axis wording (2026-11-25): D "Effect size (SD)"; H "Learning effect (SD)" against
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

# The canonical source data of the other panels, pinned by hash the same way: the published
# files are parsed from them, so a changed file must stop the build. Keyed by file name.
PANEL_SOURCE_PINS <- c(
  "Fig3C_PCA_animal_level_source_data.csv" =
    "c026c23202c1a63b7e662031a856b5296febc02511129ce107c0f51edbed1a8c",
  "Fig3E_rank_abundance_animal_level_source_data.csv" =
    "597bd3fb3b094a0346a10c1862f45da2ff168722ea02fa1b683ba04f5a848128",
  "SuppA_protein_peptide_counts_source_data.xlsx" =
    "52359117b8b253c41960d00d6a8dc55cd47abbcdbb5289c56fba7ac66e5c5056",
  "SuppB1_PCA_experimental_group_animal_level_source_data.csv" =
    "c026c23202c1a63b7e662031a856b5296febc02511129ce107c0f51edbed1a8c",
  "SuppB2_B3_PCA_acquisition_level_QC_source_data.csv" =
    "6d3421202d3747f7cbe0174359c47ebad7aea2eeab678d0541db15899ef6eb0d",
  "SuppC_celltype_identity_heatmap_animal_level_source_data.csv" =
    "d378bb86c5fc8f773bc5b6303876fc3994898f11dfde0a4f768597010fed1003",
  "SuppD_rank_abundance_animal_level_source_data.csv" =
    "ebdf363bff8a69f28b29545f9efb52de196d7160b9a304c5140814d833371535",
  "SuppE_source_data.csv" =
    "2de96a5209ae7fa89c90c14431e9c5b6b02bc07e8d9ae6b606ca3daab235aed9"
)

# Published name of a panel's further files: <panel>_<name>.csv, in file order.
PANEL_EXTRA_FILE_NAMES <- list(
  "Figure3|G" = "redundant_terms_skipped",
  "Supplementary_proteomics|F" = "blank_cells"
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
# The rule names the ontologies BP and CC; the map and the source data write GO-BP and
# GO-CC.
go_collection_names <- function(x) {
  gsub("(?<!GO-)\\b(BP|CC)\\b", "GO-\\1", x, perl = TRUE)
}
g_rule <- go_collection_names(g_rule)
ora_up <- release_read_csv(paths$ora_learning_up)
ora_up_cc <- release_read_csv(paths$ora_learning_up_cc)
ora_text <- function(d, label) {
  paste0(sum(d$p.adjust < 0.05), " of ", format(nrow(d), big.mark = ","), " ", label,
         " terms at FDR < 0.05 (minimum adjusted P ", release_format_p(min(d$p.adjust)), ")")
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
                " on the proteins at FDR < 0.05 in the Pairing comparison (n = ",
                audit_value(h_audit, "n_proteins"),
                "); r = ", fmt_r(audit_value(h_audit, "r_all_proteins_unselected")),
                " unselected; same-side control r = ",
                fmt_r(audit_value(h_audit, "control_r_learning_signature")))
h_note <- paste0(
  "DESCRIPTIVE ONLY. Pearson r = ", fmt_r(audit_value(h_audit, "corrected_r")), " across the ",
  audit_value(h_audit, "n_proteins"), " proteins at FDR < 0.05 in the Pairing comparison.",
  " The two comparisons share paired_veh on OPPOSITE sides, so a negative correlation is",
  " expected by construction; the same-side control gives r = ",
  fmt_r(audit_value(h_audit, "control_r_learning_signature")), " on the same protein set.",
  " Selecting the proteins on their Pairing FDR (x axis) additionally induces regression to",
  " the mean on y. This panel is not evidence that CNO reverses the pairing effect.")

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
    paste0("; the ", sum(!i_terms$opposed),
           if (sum(!i_terms$opposed) == 1L) " concordant term is in " else
             " concordant terms are in ",
           paste(release_sample_class_display(i_conc_classes), collapse = ", ")) else "",
  ". The opposition cannot be read as proof that CNO reverses a pairing-associated program. ",
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
  # Published files are CSV, named in the manuscript numbering: <panel>_source_data.csv
  # (Fig2D, SuppFig4F), further files <panel>_<name>.csv from PANEL_EXTRA_FILE_NAMES; all are
  # ";"-separated in order. The canonical file names are kept in a column so the
  # correspondence stays exact.
  panel_label <- paste(ifelse(figure == "Figure3", "Fig. 2", "Supp. Fig. 4"), panel)
  panel_slug <- gsub("[^A-Za-z0-9]", "", sub("^Supp[.] Fig[.]", "SuppFig", panel_label))
  dest_rel <- if (has_source_data) {
    extra_names <- PANEL_EXTRA_FILE_NAMES[[key]]
    if (length(extra_names) != length(source_data_rel) - 1L) {
      stop("Panel ", panel_label, " has ", length(source_data_rel) - 1L, " further source ",
           "file(s) but ", length(extra_names), " published name(s).", call. = FALSE)
    }
    file.path("editor_source_data", "figure_source_data",
              c(paste0(panel_slug, "_source_data.csv"),
                if (length(extra_names)) paste0(panel_slug, "_", extra_names, ".csv")))
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
      panel_unit(key, trimws(as.character(r$new_sampling_unit))),
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
# reader-facing panel source data
# --------------------------------------------------------------------------------------
# Each canonical file is checked against its pin, read as text and published as CSV with
# the headers of RELEASE_FIGURE_CSV_HEADERS and the display labels of the workbook. Every
# source column is either published or dropped by name, so a changed file stops the build.
# Values are carried over as text. Added columns are lookups only (sample IDs, hemisphere,
# gene symbols, comparison names), and the comparisons with the release tables are checks,
# not inputs.

# The display labels are the workbook's (R/release_presentation.R); fig_label() stops on a
# value with no label instead of passing it through.
FIG_CLASS_LABEL <- RELEASE_SAMPLE_CLASS_DISPLAY
FIG_CONDITION_LABEL <- RELEASE_CONDITION_DISPLAY
FIG_FAMILY_LABEL <- RELEASE_FAMILY_LABEL
FIG_ONTOLOGY_LABEL <- stats::setNames(RELEASE_ONTOLOGY_DISPLAY[c("GO_BP", "GO_CC")],
                                      c("BP", "CC"))
LEARNING_COMPARISON <- "mcherry_paired_veh_over_mcherry_unpaired_veh"

sample_meta <- release_read_tsv_plain(paths$sample_metadata)
animal_meta <- release_read_tsv_plain(paths$animal_metadata)
read_release_gz <- function(p, ...) {
  utils::read.delim(gzfile(p), sep = "\t", quote = "", stringsAsFactors = FALSE,
                    check.names = FALSE, ...)
}
differential <- read_release_gz(paths$differential)
gsea_go_bp <- read_release_gz(paths$gsea_go_bp)
untested_tbl <- utils::read.delim(paths$untested_sets, sep = "\t", quote = "",
                                  stringsAsFactors = FALSE, check.names = FALSE,
                                  colClasses = "character")
features <- read_release_gz(paths$feature_annotation, colClasses = "character")
has_acc <- !is.na(features$uniprot_accession) & nzchar(features$uniprot_accession)
acc_sym <- unique(data.frame(
  acc = features$uniprot_accession[has_acc],
  sym = ifelse(is.na(features$gene_symbol[has_acc]), "", trimws(features$gene_symbol[has_acc])),
  stringsAsFactors = FALSE))
if (anyDuplicated(acc_sym$acc)) {
  stop("The feature annotation gives more than one gene symbol for an accession.",
       call. = FALSE)
}
symbol_of_accession <- stats::setNames(acc_sym$sym, acc_sym$acc)

fig_check <- function(ok, ...) {
  if (!isTRUE(all(ok))) stop(..., call. = FALSE)
  invisible(TRUE)
}
fig_label <- function(map, x, what) {
  x <- as.character(x)
  out <- unname(map[x])
  if (anyNA(out)) {
    stop("No display label for ", what, ": ", paste(unique(x[is.na(out)]), collapse = ", "),
         call. = FALSE)
  }
  out
}
contrast_of <- function(key) {
  i <- match(key, contrasts$canonical_comparison)
  fig_check(!is.na(i), "Not a primary comparison: ", paste(unique(key[is.na(i)]), collapse = ", "))
  contrasts[i, , drop = FALSE]
}
groups_label <- function(key) {
  r <- contrast_of(key)
  paste(fig_label(FIG_CONDITION_LABEL, r$numerator_condition, "condition"), "vs",
        fig_label(FIG_CONDITION_LABEL, r$denominator_condition, "condition"))
}
comparison_label <- function(key) {
  paste0(fig_label(FIG_CLASS_LABEL, contrast_of(key)$sample_class, "compartment"), ": ",
         groups_label(key))
}
comparison_key <- function(sample_class, family) {
  i <- match(paste(sample_class, family),
             paste(contrasts$sample_class, contrasts$contrast_family))
  fig_check(!is.na(i), "No primary comparison for ",
            paste(unique(paste(sample_class, family)[is.na(i)]), collapse = ", "))
  contrasts$canonical_comparison[i]
}
# Hemisphere as Left/Right where the sample metadata records it from explicit evidence (and,
# when given, the replicate number the panel used agrees); NA otherwise.
certain_hemisphere <- function(sample_ids, replicate_group = NULL) {
  i <- match(sample_ids, sample_meta$sample_id)
  h <- sample_meta$hemisphere[i]
  ok <- !is.na(i) & h %in% c("Left", "Right") &
    sample_meta$hemisphere_evidence_source[i] %in% "explicit_evidence"
  if (!is.null(replicate_group)) {
    rg <- sample_meta$legacy_replicate_group[i]
    ok <- ok & !is.na(rg) & rg == replicate_group
  }
  ifelse(ok, h, NA_character_)
}
# Text values compared with release numbers (the CSVs carry 15 significant digits).
same_numbers <- function(text, reference) {
  x <- suppressWarnings(as.numeric(text))
  length(x) == length(reference) && !anyNA(x) && !anyNA(reference) &&
    all(abs(x - reference) <= 1e-9 * pmax(1, abs(reference)))
}
differential_rows <- function(cmp, protein_group_ids) {
  da <- differential[differential$canonical_comparison == cmp, , drop = FALSE]
  j <- match(protein_group_ids, da$protein_group_id)
  fig_check(!is.na(j), "Panel proteins missing from the release statistics of ", cmp)
  da[j, , drop = FALSE]
}
gsea_go_rows <- function(keys, by, values) {
  g <- gsea_go_bp[gsea_go_bp$canonical_comparison %in% keys, , drop = FALSE]
  table_keys <- paste(g$canonical_comparison, g[[by]], sep = "\r")
  j <- match(paste(keys, values, sep = "\r"), table_keys)
  fig_check(!is.na(j) & !table_keys[j] %in% table_keys[duplicated(table_keys)],
            "Panel gene sets not found once in the release GSEA GO-BP table")
  g[j, , drop = FALSE]
}

read_text_csv <- function(path) {
  d <- utils::read.csv(path, colClasses = "character", check.names = FALSE,
                       stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  names(d) <- sub(paste0("^", intToUtf8(0xFEFF)), "", names(d))
  d
}
# Supp. Fig. 4 A: only the first sheet holds this panel's data.
read_first_sheet <- function(path) {
  d <- as.data.frame(readxl::read_excel(path, sheet = 1L, col_types = "text",
                                        .name_repair = "minimal"),
                     stringsAsFactors = FALSE)
  attr(d, "sheets") <- readxl::excel_sheets(path)
  d
}
expect_columns <- function(d, cols, what) {
  fig_check(!anyDuplicated(names(d)) && setequal(names(d), cols),
            what, " no longer has the expected columns (", paste(cols, collapse = ", "),
            "); found: ", paste(names(d), collapse = ", "))
}
as_table <- function(cols) {
  n <- max(lengths(cols))
  fig_check(lengths(cols) %in% c(1L, n), "Published columns differ in length.")
  out <- data.frame(row.names = seq_len(n))
  for (nm in names(cols)) {
    v <- rep_len(as.character(cols[[nm]]), n)
    out[[nm]] <- ifelse(is.na(v), "", v)
  }
  rownames(out) <- NULL
  out
}
csv_field <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  q <- grepl("[\",\r\n]", x) | grepl("^\\s|\\s$", x)
  x[q] <- paste0("\"", gsub("\"", "\"\"", x[q], fixed = TRUE), "\"")
  x
}
# UTF-8 with a byte order mark, so spreadsheet programs read the labels correctly.
write_reader_csv <- function(df, path) {
  lines <- paste(csv_field(names(df)), collapse = ",")
  if (nrow(df)) lines <- c(lines, do.call(paste, c(unname(lapply(df, csv_field)), sep = ",")))
  parent <- dirname(path)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0xef, 0xbb, 0xbf)), con)
  writeBin(charToRaw(paste0(paste(enc2utf8(lines), collapse = "\n"), "\n")), con)
  invisible(list(rows = nrow(df), cols = ncol(df)))
}

# --- one builder per canonical layout --------------------------------------------------

build_pca_animal <- function(d, spec) {
  expect_columns(d, c("Sample", "AnimalID", "condition_code", "condition", "sample_class",
                      "PC1", "PC2", "PC3"), spec$what)
  i <- match(d$Sample, animal_meta$animal_level_column_name)
  fig_check(!is.na(i) & d$AnimalID == animal_meta$AnimalID[i] &
              d$sample_class == animal_meta$sample_class[i] &
              d$condition == animal_meta$condition[i] &
              d$condition_code == animal_meta$condition_code[i],
            spec$what, " disagrees with the animal-level sample metadata.")
  as_table(list(
    "Animal-level profile" = d$Sample,
    "Animal ID" = d$AnimalID,
    "Compartment" = fig_label(FIG_CLASS_LABEL, d$sample_class, "compartment"),
    "Condition" = fig_label(FIG_CONDITION_LABEL, d$condition, "condition"),
    "PC1" = d$PC1, "PC2" = d$PC2, "PC3" = d$PC3))
}

# gene_symbol holds the accession and Description the source-matrix symbols; Gene_Name is the
# symbol with the accession filled in. Only the mapped symbol is published.
build_volcano <- function(d, spec) {
  expect_columns(d, c("gene_symbol", "mapped_gene_symbol", "Gene_Name", "original_protein_id",
                      "Description", "uniprot_accession", "log2fc", "t", "pval", "padj",
                      "Significance", "labelled_in_panel"), spec$what)
  cmp <- spec$comparison
  r <- contrast_of(cmp)
  fig_check(identical(d$gene_symbol, d$uniprot_accession),
            spec$what, ": gene_symbol no longer repeats uniprot_accession.")
  da <- differential_rows(cmp, d$original_protein_id)
  fig_check(same_numbers(d$log2fc, da$effect_size_sd_units) &&
              same_numbers(d$t, da$moderated_t) && same_numbers(d$pval, da$P.Value) &&
              same_numbers(d$padj, da$adj.P.Val),
            spec$what, " does not hold the release statistics of ", cmp)
  padj <- as.numeric(d$padj); es <- as.numeric(d$log2fc)
  expected <- ifelse(padj < 0.05, ifelse(es > 0, "up", "down"), "n.s.")
  fig_check(d$Significance == expected, spec$what, ": Significance disagrees with padj and sign.")
  direction <- c(up = paste("higher in", fig_label(FIG_CONDITION_LABEL, r$numerator_condition,
                                                  "condition")),
                 down = paste("higher in", fig_label(FIG_CONDITION_LABEL, r$denominator_condition,
                                                    "condition")),
                 n.s. = "not significant")
  as_table(list(
    "Comparison" = comparison_label(cmp),
    "Gene symbol" = d$mapped_gene_symbol,
    "UniProt accession" = d$uniprot_accession,
    "Protein group" = d$original_protein_id,
    "Effect size (SD units)" = d$log2fc,
    "Moderated t" = d$t,
    "P value" = d$pval,
    "Adjusted P value (BH)" = d$padj,
    "Direction" = unname(direction[d$Significance]),
    "Labeled in panel" = d$labelled_in_panel,
    "Comparison key" = cmp))
}

# MeanLog2 is the mean over 3 animals of the per-protein standardized abundance, and the
# panel plots LinearValue = 2^MeanLog2. Only the marker categories the panel shows keep a
# label, under the panel's legend names.
RANK_MARKER_TYPES <- c("None", "General Neuron", "Neuropil/Structure", "Stress Response",
                       "Activation")
build_rank_abundance <- function(d, spec) {
  expect_columns(d, c("Genes", "Condition", "MeanLog2", "LinearValue", "Rank", "n_animals",
                      "animals", "MarkerType"), spec$what)
  fig_check(d$MarkerType %in% RANK_MARKER_TYPES, spec$what, ": unknown MarkerType.")
  category <- unname(spec$categories[d$MarkerType])
  as_table(list(
    "Compartment" = fig_label(FIG_CLASS_LABEL, sub("_.*$", "", d$Condition), "compartment"),
    "Condition" = fig_label(FIG_CONDITION_LABEL,
                            gsub("-", "_", sub("^[^_]*_", "", d$Condition), fixed = TRUE),
                            "condition"),
    "Gene symbol" = d$Genes,
    "Rank" = d$Rank,
    "Mean standardized abundance over animals (SD units)" = d$MeanLog2,
    "Plotted value (2 to the power of the mean; display transform, not an abundance)" =
      d$LinearValue,
    "Animals (n)" = d$n_animals,
    "Animals" = d$animals,
    "Marker category in panel" = ifelse(is.na(category), "", category)))
}

build_gsea_grid <- function(d, spec) {
  expect_columns(d, c("Comparison", "sample_class", "canonical_contrast", "ID", "Description",
                      "setSize", "NES", "pvalue", "p.adjust", "qvalue", "fdr_significant",
                      "selected_by_column"), spec$what)
  key <- d$canonical_contrast
  r <- contrast_of(key)
  fig_check(r$sample_class == d$sample_class & r$contrast_family == spec$family,
            spec$what, ": comparisons outside ", spec$family)
  fig_check(d$selected_by_column == "TRUE", spec$what,
            ": selected_by_column is no longer TRUE on every row.")
  g <- gsea_go_rows(key, "term_id", d$ID)
  fig_check(same_numbers(d$NES, g$NES) && same_numbers(d[["p.adjust"]], g$adjusted_p_value),
            spec$what, " does not hold the release GSEA GO-BP values.")
  as_table(list(
    "Comparison" = comparison_label(key),
    "Term ID" = d$ID,
    "Term" = d$Description,
    "Gene-set size" = d$setSize,
    "NES" = d$NES,
    "P value" = d$pvalue,
    "Adjusted P value (BH)" = d[["p.adjust"]],
    "q value" = d$qvalue,
    "FDR < 0.05" = d$fdr_significant,
    "Comparison key" = key))
}

build_lollipop <- function(d, spec) {
  expect_columns(d, c("ontology", "ID", "Description", "setSize", "leading_edge_proteins",
                      "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue",
                      "core_enrichment", "n_terms_tested", "n_FDR_significant",
                      "selection_rule", "input_file"), spec$what)
  key <- spec$comparison
  fig_check(basename(dirname(gsub("\\\\", "/", d$input_file))) == key,
            spec$what, ": input files are not all from ", key)
  members <- strsplit(d$core_enrichment, "/", fixed = TRUE)
  fig_check(lengths(members) == as.integer(d$leading_edge_proteins),
            spec$what, ": leading_edge_proteins is not the member count.")
  symbols <- vapply(members, function(m) {
    s <- unname(symbol_of_accession[m])
    fig_check(!is.na(s), spec$what, ": leading-edge accession(s) not in the feature annotation.")
    paste(ifelse(nzchar(s), s, m), collapse = "/")
  }, character(1))
  as_table(list(
    "Comparison" = comparison_label(key),
    "Gene-set collection" = fig_label(FIG_ONTOLOGY_LABEL, d$ontology, "ontology"),
    "Term ID" = d$ID,
    "Term" = d$Description,
    "Gene-set size" = d$setSize,
    "Enrichment score" = d$enrichmentScore,
    "NES" = d$NES,
    "P value" = d$pvalue,
    "Adjusted P value (BH)" = d[["p.adjust"]],
    "q value" = d$qvalue,
    "Leading-edge proteins (n)" = d$leading_edge_proteins,
    "Leading-edge genes" = symbols,
    "Leading-edge proteins (UniProt accessions)" = d$core_enrichment,
    "Gene sets reported (n)" = d$n_terms_tested,
    "Gene sets at FDR < 0.05 (n)" = d$n_FDR_significant,
    "Selection rule" = go_collection_names(d$selection_rule),
    "Comparison key" = key))
}

build_redundant_terms <- function(d, spec) {
  expect_columns(d, c("ontology", "ID", "Description", "p.adjust", "max_similarity_to_kept"),
                 spec$what)
  key <- spec$comparison
  as_table(list(
    "Comparison" = comparison_label(key),
    "Gene-set collection" = fig_label(FIG_ONTOLOGY_LABEL, d$ontology, "ontology"),
    "Term ID" = d$ID,
    "Term" = d$Description,
    "Adjusted P value (BH)" = d[["p.adjust"]],
    "Highest leading-edge Jaccard index with a plotted term" = d$max_similarity_to_kept,
    "Comparison key" = key))
}

build_scatter <- function(d, spec) {
  expect_columns(d, c("uniprot_accession", "original_protein_id", "mapped_gene_symbol",
                      "log2fc_learning", "t_learning", "padj_learning", "log2fc_cno", "t_cno",
                      "padj_cno", "log2fc_unpaired_cno", "padj_unpaired_cno",
                      "learning_signature"), spec$what)
  cmp <- spec$comparisons
  ok <- vapply(names(cmp), function(nm) {
    da <- differential_rows(cmp[[nm]], d$original_protein_id)
    same_numbers(d[[paste0("log2fc_", nm)]], da$effect_size_sd_units) &&
      same_numbers(d[[paste0("padj_", nm)]], da$adj.P.Val) &&
      (!paste0("t_", nm) %in% names(d) || same_numbers(d[[paste0("t_", nm)]], da$moderated_t))
  }, logical(1))
  fig_check(ok, spec$what, " does not hold the release statistics of ",
            paste(cmp[!ok], collapse = ", "))
  fig_check(d$learning_signature == ifelse(as.numeric(d$padj_learning) < 0.05, "TRUE", "FALSE"),
            spec$what, ": learning_signature is not FDR < 0.05 in the learning comparison.")
  compartment <- unique(contrast_of(cmp)$sample_class)
  fig_check(length(compartment) == 1L, spec$what, ": comparisons from more than one compartment.")
  g <- vapply(cmp, groups_label, character(1))
  cols <- list("Compartment" = fig_label(FIG_CLASS_LABEL, compartment, "compartment"),
               "Gene symbol" = d$mapped_gene_symbol,
               "UniProt accession" = d$uniprot_accession,
               "Protein group" = d$original_protein_id)
  for (nm in names(cmp)) {
    cols[[paste0("Effect size (SD units): ", g[[nm]])]] <- d[[paste0("log2fc_", nm)]]
    if (paste0("t_", nm) %in% names(d)) {
      cols[[paste0("Moderated t: ", g[[nm]])]] <- d[[paste0("t_", nm)]]
    }
    cols[[paste0("Adjusted P value (BH): ", g[[nm]])]] <- d[[paste0("padj_", nm)]]
  }
  cols[[paste0("Plotted: FDR < 0.05 in ", g[["learning"]])]] <- d$learning_signature
  as_table(cols)
}

# Rows of the heatmap: "Memory Engram" (dataset "Memory Ensemble") is the pairing contrast,
# "Inhibition" ("Effects Inhibition") CNO in paired animals.
HEATMAP_ROW_FAMILY <- c("Memory Engram" = "paired_veh_vs_unpaired_veh",
                        "Inhibition" = "paired_cno_vs_paired_veh")
HEATMAP_DATASET_FAMILY <- c("Memory Ensemble" = "paired_veh_vs_unpaired_veh",
                            "Effects Inhibition" = "paired_cno_vs_paired_veh")
build_heatmap <- function(d, spec) {
  expect_columns(d, c("Description", "Comparison_Group", "sample_class", "Dataset", "NES",
                      "p.adjust", "Column", "Row"), spec$what)
  family <- unname(HEATMAP_ROW_FAMILY[d$Row])
  fig_check(!is.na(family) & family == unname(HEATMAP_DATASET_FAMILY[d$Dataset]),
            spec$what, ": Row and Dataset do not name the same contrast family.")
  key <- comparison_key(d$sample_class, family)
  g <- gsea_go_rows(key, "term_name", d$Description)
  fig_check(same_numbers(d$NES, g$NES) && same_numbers(d[["p.adjust"]], g$adjusted_p_value),
            spec$what, " does not hold the release GSEA GO-BP values.")
  as_table(list(
    "Comparison" = comparison_label(key),
    "Comparison type" = fig_label(FIG_FAMILY_LABEL, family, "comparison type"),
    "Term ID" = g$term_id,
    "Term" = d$Description,
    "NES" = d$NES,
    "Adjusted P value (BH)" = d[["p.adjust"]],
    "Comparison key" = key))
}

build_ewce <- function(d, spec) {
  expect_columns(d, c("CellType", "Column", "sample_class", "condition", "Target", "TopN",
                      "AnnotLevel", "fold_change", "sd_from_mean", "p", "q_target",
                      "q_global", "Significant_Global"), spec$what)
  fig_check(d$Target == paste0(d$sample_class, "__", d$condition),
            spec$what, ": Target is not sample class and condition.")
  as_table(list(
    "Compartment" = fig_label(FIG_CLASS_LABEL, d$sample_class, "compartment"),
    "Condition" = fig_label(FIG_CONDITION_LABEL, d$condition, "condition"),
    "Cell type" = ifelse(d$CellType == "(none)",
                         paste0("unassigned level-", d$AnnotLevel, " cluster"), d$CellType),
    "Proteins in list (n)" = d$TopN,
    "Annotation level" = d$AnnotLevel,
    "Fold change (observed/null)" = d$fold_change,
    "SD from null mean" = d$sd_from_mean,
    "P value" = d$p,
    "Adjusted P value, within list (BH)" = d$q_target,
    "Adjusted P value, global (BH)" = d$q_global,
    "FDR < 0.05, global" = d$Significant_Global))
}

build_pca_acquisition <- function(d, spec) {
  expect_columns(d, c("acquisition", "shortname", "PC1", "PC2", "Hemisphere",
                      "Collection plate", "sample_class", "AnimalID", "condition_code"),
                 spec$what)
  i <- match(d$acquisition, sample_meta$sample_id)
  cls <- release_normalize_annotation_label(d$sample_class)
  fig_check(!is.na(i) & d$shortname == sample_meta$legacy_shortname[i] &
              cls == sample_meta$sample_class[i] & d$AnimalID == sample_meta$AnimalID[i] &
              d$condition_code == sample_meta$condition_code[i] &
              d[["Collection plate"]] == sample_meta$collection_plate[i],
            spec$what, " disagrees with the sample metadata.")
  hemisphere <- certain_hemisphere(d$acquisition, sub("^Hemisphere ", "", d$Hemisphere))
  as_table(list(
    "Sample ID" = d$acquisition,
    "Animal ID" = d$AnimalID,
    "Hemisphere" = ifelse(is.na(hemisphere), d$Hemisphere, hemisphere),
    "Compartment" = fig_label(FIG_CLASS_LABEL, cls, "compartment"),
    "Condition" = fig_label(FIG_CONDITION_LABEL, sample_meta$condition[i], "condition"),
    "Collection plate" = d[["Collection plate"]],
    "PC1" = d$PC1, "PC2" = d$PC2))
}

# ids such as 1_n01_s4<U+2212>a1_1_8782 are the run short names (Plate1_N01_S4-A1_1_8782) in
# lower case, without "Plate", with a Unicode minus. Celltype is the sample class before the
# correction of six acquisitions; both are published.
build_acquisition_counts <- function(d, spec) {
  expect_columns(d, c("Celltype", "Group", "celltype_group", "id", "proteins", "peptides"),
                 spec$what)
  fig_check(d$celltype_group == paste0(d$Celltype, "_", d$Group),
            spec$what, ": celltype_group is not Celltype and Group.")
  run <- tolower(paste0("plate", gsub(intToUtf8(0x2212), "-", d$id, fixed = TRUE)))
  i <- match(run, tolower(sample_meta$legacy_shortname))
  fig_check(!anyNA(i) && !anyDuplicated(i) && length(i) == nrow(sample_meta),
            spec$what, " does not name every acquisition of the sample metadata once.")
  fig_check(d$Group == sample_meta$condition_code[i],
            spec$what, ": Group is not the condition code.")
  before <- release_normalize_annotation_label(d$Celltype)
  fig_check(!is.na(before) & before == sample_meta$original_sample_class[i],
            spec$what, ": Celltype is not the sample class before correction.")
  fig_check(grepl("^[0-9]+$", d$proteins) & grepl("^[0-9]+$", d$peptides),
            spec$what, ": counts are not whole numbers.")
  as_table(list(
    "Sample ID" = sample_meta$sample_id[i],
    "Animal ID" = sample_meta$AnimalID[i],
    "Hemisphere" = certain_hemisphere(sample_meta$sample_id[i]),
    "Compartment" = fig_label(FIG_CLASS_LABEL, sample_meta$sample_class[i], "compartment"),
    "Compartment as recorded" = fig_label(FIG_CLASS_LABEL, before, "compartment"),
    "Condition" = fig_label(FIG_CONDITION_LABEL, sample_meta$condition[i], "condition"),
    "Proteins (n)" = d$proteins,
    "Peptides (n)" = d$peptides))
}

build_cross_compartment <- function(d, spec) {
  expect_columns(d, c("Comparison", "contrast_id", "historical_alias", "ID", "Description",
                      "setSize", "NES", "enrichmentScore", "pvalue", "p.adjust", "qvalue",
                      "n_ranked_genes", "fdr_significant"), spec$what)
  parts <- strsplit(d$Comparison, " vs ", fixed = TRUE)
  fig_check(lengths(parts) == 2L, spec$what, ": Comparison is not '<a> vs <b>'.")
  a <- vapply(parts, `[[`, character(1), 1L)
  b <- vapply(parts, `[[`, character(1), 2L)
  fig_check(d$contrast_id == paste0(a, "_vs_", b), spec$what,
            ": contrast_id does not match Comparison.")
  as_table(list(
    "Comparison" = paste(fig_label(FIG_CLASS_LABEL, a, "compartment"), "vs",
                         fig_label(FIG_CLASS_LABEL, b, "compartment")),
    "Term ID" = d$ID,
    "Term" = d$Description,
    "Gene-set size" = d$setSize,
    "NES" = d$NES,
    "Enrichment score" = d$enrichmentScore,
    "P value" = d$pvalue,
    "Adjusted P value (BH)" = d[["p.adjust"]],
    "q value" = d$qvalue,
    "Ranked proteins (n)" = d$n_ranked_genes,
    "FDR < 0.05" = d$fdr_significant))
}

build_blank_cells <- function(d, spec) {
  expect_columns(d, c("Comparison", "sample_class", "Description", "reason"), spec$what)
  key <- comparison_key(d$sample_class, spec$family)
  u <- untested_tbl[untested_tbl$analysis_role == "canonical" & untested_tbl$ontology == "GO_BP", ,
                    drop = FALSE]
  j <- match(paste(key, d$Description, sep = "\r"),
             paste(u$canonical_comparison, u$term_name, sep = "\r"))
  fig_check(!is.na(j), spec$what, ": blank cells that are not gene sets without a P value.")
  as_table(list(
    "Comparison" = comparison_label(key),
    "Term ID" = u$term_id[j],
    "Term" = d$Description,
    "Reason" = "fgsea gave no P value for this gene set (listed in Supplementary Table 2l)",
    "Comparison key" = key))
}

PANEL_PUBLICATION <- list(
  "Figure3|C" = list(build = list(build_pca_animal)),
  "Figure3|D" = list(build = list(build_volcano), comparison = LEARNING_COMPARISON),
  "Figure3|E" = list(build = list(build_rank_abundance),
                     categories = c("Activation" = "Mitochondrial",
                                    "Stress Response" = "Growth-Related Plasticity")),
  "Figure3|F" = list(build = list(build_gsea_grid), family = "paired_veh_vs_unpaired_veh"),
  "Figure3|G" = list(build = list(build_lollipop, build_redundant_terms),
                     comparison = LEARNING_COMPARISON),
  "Figure3|H" = list(
    build = list(build_scatter),
    comparisons = c(learning = LEARNING_COMPARISON,
                    cno = "mcherry_paired_cno_over_mcherry_paired_veh",
                    unpaired_cno = "mcherry_unpaired_cno_over_mcherry_unpaired_veh")),
  "Figure3|I" = list(build = list(build_heatmap)),
  "Supplementary_proteomics|A" = list(build = list(build_acquisition_counts)),
  "Supplementary_proteomics|B1" = list(build = list(build_pca_animal)),
  "Supplementary_proteomics|B2" = list(build = list(build_pca_acquisition)),
  "Supplementary_proteomics|B3" = list(build = list(build_pca_acquisition)),
  "Supplementary_proteomics|C" = list(build = list(build_ewce)),
  "Supplementary_proteomics|D" = list(build = list(build_rank_abundance),
                                      categories = c("General Neuron" = "Neuron Soma",
                                                     "Neuropil/Structure" = "Neuropil/Structure")),
  "Supplementary_proteomics|E" = list(build = list(build_cross_compartment)),
  "Supplementary_proteomics|F" = list(build = list(build_gsea_grid, build_blank_cells),
                                      family = "paired_cno_vs_paired_veh")
)

published <- 0L
missing_source_data <- character(0)
supp_a_out <- NULL
split_paths <- function(x) strsplit(x, ";", fixed = TRUE)[[1]]
for (i in seq_len(nrow(figure_map))) {
  if (is.na(figure_map$panel_source_data_canonical_path[i]) ||
      is.na(figure_map$panel_source_data_release_path[i])) next
  key <- paste(figure_map$figure[i], figure_map$panel[i], sep = "|")
  label <- figure_map$panel_label[i]
  srcs <- split_paths(figure_map$panel_source_data_canonical_path[i])
  rels <- split_paths(figure_map$panel_source_data_release_path[i])
  if (!all(file.exists(srcs))) {
    missing_source_data <- c(missing_source_data,
                             paste0(label, ": ", paste(srcs[!file.exists(srcs)], collapse = ";")))
    figure_map$panel_source_data_release_path[i] <- NA_character_
    next
  }
  spec <- PANEL_PUBLICATION[[key]]
  if (is.null(spec) || length(spec$build) != length(srcs)) {
    stop("No reader-facing layout for the source data of panel ", label, call. = FALSE)
  }
  shas <- split_paths(figure_map$panel_source_data_sha256[i])
  pins <- if (key %in% REBUILT_PANELS$key) rebuilt_files(key)$sha else
    unname(PANEL_SOURCE_PINS[basename(srcs)])
  if (length(pins) != length(shas) || anyNA(pins) || !identical(pins, shas)) {
    stop("Source data of panel ", label, " no longer match the pinned hash.", call. = FALSE)
  }
  for (k in seq_along(srcs)) {
    spec$what <- paste0("Source data of ", label, " (", basename(srcs[[k]]), ")")
    is_xlsx <- grepl("[.]xlsx$", srcs[[k]], ignore.case = TRUE)
    d <- if (is_xlsx) read_first_sheet(srcs[[k]]) else read_text_csv(srcs[[k]])
    out <- spec$build[[k]](d, spec)
    fig_check(identical(names(out), RELEASE_FIGURE_CSV_HEADERS[[basename(rels[[k]])]]),
              spec$what, ": the published headers are not those of RELEASE_FIGURE_CSV_HEADERS.")
    w <- write_reader_csv(out, release_path(rels[[k]]))
    release_register(rels[[k]],
                     paste0("reader-facing source data for panel ", label, " (from ",
                            basename(srcs[[k]]), ")"),
                     srcs[[k]], shas[[k]], STAGE, "csv")
    published <- published + 1L
    release_log("    ", format(basename(rels[[k]]), width = 34), " ", w$rows, " x ", w$cols,
                if (is_xlsx) paste0("  (sheet 1 of ", length(attr(d, "sheets")),
                                    " of the source workbook)") else "")
    if (key == "Supplementary_proteomics|A") supp_a_out <- out
    # Rank-abundance curves: one row per gene symbol; the count is read from the file.
    if (identical(spec$build[[k]], build_rank_abundance)) {
      per_curve <- unique(as.integer(table(paste(out$Compartment, out$Condition))))
      fig_check(length(per_curve) == 1L && !anyDuplicated(paste(out$Compartment, out$Condition,
                                                                out[["Gene symbol"]])),
                spec$what, ": the curves do not hold one row per gene symbol each.")
      figure_map$source_rows_filter[i] <- paste0(
        figure_map$source_rows_filter[i], "; one row per gene symbol, ",
        format(per_curve, big.mark = ","), " per curve (proteins without a gene symbol are ",
        "not plotted)")
    }
  }
}
release_log("  published ", published, " reader-facing panel source-data files")
if (length(missing_source_data)) {
  release_log("  NOTE: no source-data file on disk for: ",
              paste(missing_source_data, collapse = " | "))
}
# Supp. Fig. 4 A: rows, count ranges and the relabeled acquisitions, read from the published
# file.
a_row <- which(paste(figure_map$figure, figure_map$panel, sep = "|") ==
                 "Supplementary_proteomics|A")
if (length(a_row) == 1L && !is.null(supp_a_out)) {
  count_range <- function(x) {
    paste(format(range(as.integer(x)), big.mark = ",", trim = TRUE), collapse = "-")
  }
  figure_map$source_rows_filter[a_row] <- paste("all", nrow(supp_a_out), "acquisitions")
  figure_map$key_numbers_revised[a_row] <- paste0(
    count_range(supp_a_out[["Proteins (n)"]]), " proteins and ",
    count_range(supp_a_out[["Peptides (n)"]]), " peptides per acquisition")
  relabeled <- supp_a_out$Compartment != supp_a_out[["Compartment as recorded"]]
  if (any(relabeled)) {
    n_words <- c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten",
                 "Eleven", "Twelve")
    n_text <- if (sum(relabeled) <= length(n_words)) n_words[[sum(relabeled)]] else
      as.character(sum(relabeled))
    animals <- sort(unique(supp_a_out[["Animal ID"]][relabeled]))
    animals_text <- if (length(animals) == 1L) animals else
      paste(paste(animals[-length(animals)], collapse = ", "), "and", animals[[length(animals)]])
    figure_map$interpretation_note[a_row] <- paste0(
      figure_map$interpretation_note[a_row], " ", n_text,
      if (sum(relabeled) == 1L) " acquisition of " else " acquisitions of ", animals_text,
      if (sum(relabeled) == 1L) " was" else " were",
      " listed under another compartment before the label correction; the source data",
      " give both labels.")
  }
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
                         release_path("provenance", "figure_source_map.tsv"))
release_register("provenance/figure_source_map.tsv",
                 "per-panel provenance: analysis, source table, row filter, unit, scripts",
                 sources, hashes, STAGE, "tsv")
release_log("  wrote figure_source_map.tsv (", w$rows, "x", w$cols, ")")
release_log("stage 06 complete")
