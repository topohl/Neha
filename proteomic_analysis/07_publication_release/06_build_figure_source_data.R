#!/usr/bin/env Rscript

# Publication release, stage 06 -- figure source data and the panel provenance map.
#
# Produces
#   editor_source_data/figure_source_map.tsv          one row per manuscript panel
#   editor_source_data/figure_source_data/**          the per-panel source-data files,
#                                                     copied byte-for-byte
#
# Built from the finalized revision manifests, which are the last word on every panel:
#   manifests/regenerated_panel_manifest.csv          what was done, with the numbers
#   manifests/original_figure_panel_source_map.csv    where the original panel came from
#   06_manuscript_figure_revision/SCRIPT_PROVENANCE.csv  executed vs repository script hashes
#
# Two supersessions matter and are resolved here rather than left ambiguous:
#   * FINAL_REPORT.md (11:59) retired Figure 3 G and blocked Supplementary E.
#     RETAINED_PANELS_3G_3H_3I_SuppE_REPORT.md (14:13) superseded both decisions, and
#     regenerated_panel_manifest.csv (14:11) records the superseding outcome. The manifest
#     is therefore the authority, and this stage reads it.
#   * Panels 3G / 3H / 3I / Supp E have a later "final" script than their first-pass one.
#     Both are recorded: generating_script (first pass) and final_revision_script.
#
# The frozen scripts under 06_manuscript_figure_revision/ are provenance snapshots and are
# only ever read here.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_validation.R"))
release_require("digest")

DATA_ROOT <- release_data_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/06_build_figure_source_data.R"

release_banner("stage 06 -- figure source data")

REVISION_ROOT <- file.path(DATA_ROOT, "03_output", "reviewer_revision_animal_level_20260827")

paths <- list(
  regenerated = file.path(REVISION_ROOT, "manifests", "regenerated_panel_manifest.csv"),
  original_map = file.path(REVISION_ROOT, "manifests", "original_figure_panel_source_map.csv"),
  script_provenance = file.path(REPO_ROOT, "06_manuscript_figure_revision",
                                "SCRIPT_PROVENANCE.csv"),
  final_report = file.path(REVISION_ROOT, "manifests", "FINAL_REPORT.md"),
  retained_report = file.path(REVISION_ROOT, "manifests",
                              "RETAINED_PANELS_3G_3H_3I_SuppE_REPORT.md")
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
# cannot drift from the frozen provenance record.

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

panel_to_script <- list()
for (i in seq_len(nrow(provenance))) {
  refs <- parse_panel_refs(provenance$corrected_source_script_reference[i])
  for (key in refs) {
    panel_to_script[[key]] <- c(panel_to_script[[key]],
                                basename(as.character(provenance$repository_snapshot_path[i])))
  }
}
release_log("  SCRIPT_PROVENANCE.csv maps ", length(panel_to_script),
            " panels to a repository snapshot script")

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
  "Figure3|G" = "ORA GO-BP lollipop over the mCherry FDR-significant upregulated list",
  "Figure3|H" = "scatter of learning log2FC against CNO log2FC, mCherry (descriptive)",
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
  "Figure3|G" = "enrichment/primary_ORA_GO_BP.tsv.gz",
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
  "Figure3|D" = "canonical_comparison == mcherry_paired_veh_over_mcherry_unpaired_veh",
  "Figure3|E" = "all 5349 proteins; means over the 3 animals per sample_class x condition",
  "Figure3|F" = "contrast_family == paired_veh_vs_unpaired_veh AND analysis_role == canonical AND ontology == GO_BP",
  "Figure3|G" = "canonical_comparison == mcherry_paired_veh_over_mcherry_unpaired_veh AND query_list == fdr_significant_higher_in_numerator",
  "Figure3|H" = "canonical_comparison in {mcherry_paired_veh_over_mcherry_unpaired_veh, mcherry_paired_cno_over_mcherry_paired_veh}",
  "Figure3|I" = "contrast_family in {paired_veh_vs_unpaired_veh, paired_cno_vs_paired_veh} AND analysis_role == canonical AND adjusted_p_value <= 0.05 in BOTH within a sample class",
  "Supplementary_proteomics|B1" = "all 5349 proteins x 48 animal-level units",
  "Supplementary_proteomics|B2" = "all 5349 proteins x 96 acquisitions; coloured by sample_metadata hemisphere",
  "Supplementary_proteomics|B3" = "all 5349 proteins x 96 acquisitions; coloured by sample_metadata collection_plate",
  "Supplementary_proteomics|C" = "ewce_analysis_type == Baseline AND is_primary_setting == TRUE",
  "Supplementary_proteomics|D" = "sample_class in {neuron, neuropil} AND condition == unpaired_veh",
  "Supplementary_proteomics|F" = "contrast_family == paired_cno_vs_paired_veh AND analysis_role == canonical AND ontology == GO_BP"
)

PANEL_NOTE <- c(
  "Figure3|G" = paste("The corrected animal-level ORA returns NOTHING significant:",
                      "3639 GO-BP terms tested, 0 at FDR < 0.05, minimum adjusted P = 0.231.",
                      "The 15 plotted terms are the lowest adjusted P among terms with",
                      "Count >= 10 and are NOT significant. Colour encodes adjusted P on the",
                      "historical ramp, so the legend numbers, not the hues, carry the",
                      "significance information. An earlier note retiring this panel was",
                      "superseded; the panel is retained as a null result."),
  "Figure3|H" = paste("DESCRIPTIVE ONLY. The two contrasts share paired_veh on OPPOSITE",
                      "sides, so a negative correlation is expected by construction. The",
                      "same-side control gives r = +0.90 on the same protein set.",
                      "Conditioning x on learning FDR additionally induces regression to the",
                      "mean on y. This panel is not evidence that inhibition reverses",
                      "learning."),
  "Figure3|I" = paste("Row 1 (paired_veh / unpaired_veh) and row 2 (paired_cno / paired_veh)",
                      "share paired_veh on OPPOSITE sides, so the uniform opposition of all",
                      "38 overlapping terms is the expected structural consequence of the",
                      "design and cannot be read as proof that CNO reverses a",
                      "pairing-associated programme."),
  "Figure3|F" = paste("The mCherry column contributes no FDR-significant GO-BP term.",
                      "That emptiness is a result and is left visible."),
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
                                        "Figure 3 E.")
)

# --------------------------------------------------------------------------------------
# assemble the map
# --------------------------------------------------------------------------------------

orig_key <- paste(original_map$figure, trimws(original_map$panel), sep = "|")

rows <- lapply(seq_len(nrow(regenerated)), function(i) {
  r <- regenerated[i, , drop = FALSE]
  figure <- trimws(as.character(r$figure))
  panel <- trimws(as.character(r$panel))
  key <- paste(figure, panel, sep = "|")

  final_script <- panel_to_script[[key]]
  final_script <- if (is.null(final_script)) NA_character_ else
    paste(unique(final_script), collapse = ";")

  oi <- match(key, orig_key)
  source_data_rel <- trimws(as.character(r$source_data_path))
  has_source_data <- nzchar(source_data_rel) && !source_data_rel %in% c("n/a", "NA")
  source_data_abs <- if (has_source_data) file.path(REVISION_ROOT, source_data_rel) else NA_character_
  # Destination names are flattened to <panel>_source_data.<ext>. The canonical layout is
  # nested three levels deep with long descriptive filenames, and reproducing it pushed
  # paths past the Windows 260-character limit. The original filename is kept in a column
  # so the correspondence stays exact.
  panel_slug <- paste0(ifelse(figure == "Figure3", "Fig3", "Supp"), panel)
  dest_rel <- if (has_source_data) {
    file.path("editor_source_data", "figure_source_data",
              paste0(panel_slug, "_source_data.", tolower(tools::file_ext(source_data_rel))))
  } else NA_character_

  data.frame(
    figure = figure,
    panel = panel,
    panel_label = paste(ifelse(figure == "Figure3", "Fig 3", "Supp"), panel),
    panel_title = trimws(as.character(r$panel_title)),
    analysis = unname(ifelse(key %in% names(PANEL_ANALYSIS), PANEL_ANALYSIS[[key]],
                             NA_character_)),
    source_table = unname(ifelse(key %in% names(PANEL_RELEASE_TABLE),
                                 PANEL_RELEASE_TABLE[[key]], NA_character_)),
    source_rows_filter = unname(ifelse(key %in% names(PANEL_ROW_FILTER),
                                       PANEL_ROW_FILTER[[key]], NA_character_)),
    statistical_unit = trimws(as.character(r$new_sampling_unit)),
    original_statistical_unit = trimws(as.character(r$old_sampling_unit)),
    primary_or_secondary = unname(ifelse(key %in% names(PANEL_PRIMARY),
                                         PANEL_PRIMARY[[key]], NA_character_)),
    canonical_source = trimws(as.character(r$corrected_canonical_input)),
    generating_script = if (is.na(oi)) NA_character_ else
      trimws(as.character(original_map$original_generating_script[oi])),
    corrected_source_script = trimws(as.character(r$corrected_source_script)),
    final_revision_script = final_script,
    revision_status = trimws(as.character(r$status)),
    displayed_statistics_changed = trimws(as.character(r$displayed_statistics_changed)),
    biological_interpretation_changed = trimws(as.character(r$biological_interpretation_changed)),
    key_numbers_original = trimws(as.character(r$key_numbers_old)),
    key_numbers_revised = trimws(as.character(r$key_numbers_new)),
    panel_source_data_release_path = dest_rel,
    panel_source_data_original_filename = if (has_source_data) basename(source_data_rel) else
      NA_character_,
    panel_source_data_canonical_path = source_data_abs,
    panel_source_data_sha256 = if (has_source_data && file.exists(source_data_abs))
      release_sha256(source_data_abs) else NA_character_,
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

# --------------------------------------------------------------------------------------
# copy the per-panel source data verbatim
# --------------------------------------------------------------------------------------

copied <- 0L
missing_source_data <- character(0)
for (i in seq_len(nrow(figure_map))) {
  src <- figure_map$panel_source_data_canonical_path[i]
  rel <- figure_map$panel_source_data_release_path[i]
  if (is.na(src) || is.na(rel)) next
  if (!file.exists(src)) {
    missing_source_data <- c(missing_source_data, paste0(figure_map$panel_label[i], ": ", src))
    figure_map$panel_source_data_release_path[i] <- NA_character_
    next
  }
  dest <- release_path(rel)
  if (!file.copy(src, dest, overwrite = TRUE, copy.date = FALSE)) {
    stop("Could not copy panel source data: ", src, call. = FALSE)
  }
  if (!identical(release_sha256(dest), figure_map$panel_source_data_sha256[i])) {
    stop("Copied panel source data does not match its source hash: ", rel, call. = FALSE)
  }
  release_register(rel, paste0("source data for panel ", figure_map$panel_label[i]),
                   src, figure_map$panel_source_data_sha256[i], STAGE,
                   tolower(tools::file_ext(src)))
  copied <- copied + 1L
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
             paths$retained_report)
hashes <- vapply(sources, release_sha256, character(1), USE.NAMES = FALSE)

w <- release_write_table(figure_map,
                         release_path("editor_source_data", "figure_source_map.tsv"))
release_register("editor_source_data/figure_source_map.tsv",
                 "per-panel provenance: analysis, source table, row filter, unit, scripts",
                 sources, hashes, STAGE, "tsv")
release_log("  wrote figure_source_map.tsv (", w$rows, "x", w$cols, ")")
release_log("stage 06 complete")
