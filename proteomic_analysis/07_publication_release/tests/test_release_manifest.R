#!/usr/bin/env Rscript

# Contract: the release manifest is a completeness claim. Every file in the package is
# listed, every listed hash validates, and the workbook and SDRF meet their own hygiene
# rules.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("07_publication_release", "tests", "test_release_manifest.R")
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/",
                           mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "project_path_utils.R"))
LAYER <- file.path(repo_root, "07_publication_release")
source(file.path(LAYER, "R", "release_utils.R"))
source(file.path(LAYER, "R", "release_validation.R"))

failures <- 0L
expect <- function(condition, message) {
  if (!isTRUE(condition)) {
    failures <<- failures + 1L
    cat("  [FAIL]", message, "\n")
  } else {
    cat("  [ok]  ", message, "\n")
  }
}

OUT_ROOT <- release_output_root()
man_path <- file.path(OUT_ROOT, "provenance", "release_manifest.tsv")
sums_path <- file.path(OUT_ROOT, "provenance", "SHA256SUMS.txt")
if (!file.exists(man_path) || !file.exists(sums_path)) {
  cat("No built publication release at", OUT_ROOT, "\n")
  cat("Skipping publication release manifest contracts.\n")
  quit(save = "no", status = 0L)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  cat("digest not installed; skipping publication release manifest contracts.\n")
  quit(save = "no", status = 0L)
}

rd <- function(p) read.delim(p, sep = "\t", quote = "", stringsAsFactors = FALSE,
                             check.names = FALSE)
man <- rd(man_path)

cat("=== completeness ===\n")
on_disk <- gsub("\\\\", "/", list.files(OUT_ROOT, recursive = TRUE))
should_list <- setdiff(on_disk, c("provenance/release_manifest.tsv",
                                  "provenance/SHA256SUMS.txt"))
expect(length(setdiff(should_list, man$relative_path)) == 0L,
       sprintf("every released file is listed (%d files)", nrow(man)))
expect(length(setdiff(man$relative_path, on_disk)) == 0L,
       "every listed file exists on disk")
expect(!anyDuplicated(man$relative_path), "no duplicate paths in the manifest")

cat("\n=== required manifest fields ===\n")
required_cols <- c("relative_path", "file_role", "format", "rows_if_tabular",
                   "columns_if_tabular", "sha256", "source_artifact",
                   "canonical_source_sha256", "generated_by", "git_commit", "timestamp_utc")
missing_cols <- setdiff(required_cols, names(man))
expect(length(missing_cols) == 0L,
       sprintf("manifest carries every required column%s",
               if (length(missing_cols)) paste0(" (missing: ",
                                                paste(missing_cols, collapse = ", "), ")")
               else ""))
expect(all(nzchar(man$file_role)), "every file has a role")
expect(all(nzchar(man$generated_by)), "every file names its generating script")
expect(all(nzchar(man$sha256)) && all(nchar(man$sha256) == 64L),
       "every file has a 64-character SHA256")
tabular <- man$format %in% c("tsv", "tsv.gz", "csv")
expect(all(!is.na(man$rows_if_tabular[tabular])) &&
         all(!is.na(man$columns_if_tabular[tabular])),
       sprintf("every tabular file reports its shape (%d tabular files)", sum(tabular)))
expect(length(unique(man$git_commit)) == 1L, "one git commit stamped across the manifest")

cat("\n=== hashes validate ===\n")
bad <- vapply(seq_len(nrow(man)), function(i) {
  p <- file.path(OUT_ROOT, man$relative_path[i])
  !file.exists(p) || !identical(release_sha256(p), man$sha256[i])
}, logical(1))
expect(!any(bad),
       sprintf("all %d manifest hashes validate against the files on disk%s", nrow(man),
               if (any(bad)) paste0(" -- first mismatch: ", man$relative_path[which(bad)[1]])
               else ""))

sums <- readLines(sums_path, warn = FALSE)
expect(length(sums) >= nrow(man),
       sprintf("SHA256SUMS covers at least the manifest contents (%d lines)", length(sums)))
expect(any(grepl("provenance/release_manifest.tsv", sums, fixed = TRUE)),
       "SHA256SUMS covers the manifest itself")
sums_paths <- sub("^[0-9a-f]{64}  ", "", sums)
expect(length(setdiff(setdiff(on_disk, "provenance/SHA256SUMS.txt"), sums_paths)) == 0L,
       "SHA256SUMS covers every file except itself")
bad_sums <- vapply(seq_along(sums), function(i) {
  h <- substr(sums[i], 1, 64)
  p <- file.path(OUT_ROOT, sums_paths[i])
  !file.exists(p) || !identical(release_sha256(p), h)
}, logical(1))
expect(!any(bad_sums), "every SHA256SUMS entry validates")

cat("\n=== the expected package structure is present ===\n")
for (f in c("README_DATA.md", "metadata/sample_metadata.tsv",
            "metadata/animal_level_sample_metadata.tsv",
            "metadata/primary_contrast_manifest.tsv",
            "metadata/secondary_analysis_manifest.tsv", "metadata/data_dictionary.tsv",
            "processed_data/protein_abundance_measurement_level.tsv.gz",
            "processed_data/protein_abundance_animal_level.tsv.gz",
            "processed_data/protein_feature_annotation.tsv.gz",
            "differential_analysis/primary_differential_proteins.tsv.gz",
            "enrichment/primary_GSEA_GO_BP.tsv.gz", "enrichment/primary_GSEA_KEGG.tsv.gz",
            "enrichment/primary_ORA_GO_BP.tsv.gz",
            "enrichment/GSEA_log2FC_sensitivity.tsv.gz", "enrichment/primary_EWCE.tsv.gz",
            "editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx",
            "editor_source_data/REVISION_PROTEOMICS_DATA_CHANGELOG.md",
            "editor_source_data/figure_source_map.tsv",
            "pride/sdrf.tsv", "pride/SDRF_MISSING_METADATA.md", "pride/README_PRIDE.md",
            "provenance/data_lineage.tsv", "provenance/software_versions.tsv",
            "provenance/analysis_parameters.tsv",
            "provenance/sessionInfo_release.txt",
            "provenance/UPSTREAM_PREPROCESSING_GAP.md")) {
  expect(f %in% man$relative_path, sprintf("manifest lists %s", f))
}

cat("\n=== editor workbook hygiene ===\n")
inv_path <- file.path(OUT_ROOT, "editor_source_data", "workbook_sheet_inventory.tsv")
if (!file.exists(inv_path)) {
  expect(FALSE, "workbook_sheet_inventory.tsv exists")
} else {
  wb <- rd(inv_path)
  sheet_cols <- stats::setNames(strsplit(wb$columns, ";", fixed = TRUE), wb$sheet)
  problems <- release_check_workbook_hygiene(sheet_cols)
  expect(length(problems) == 0L,
         sprintf("no duplicate or unnamed workbook columns across %d sheets%s", nrow(wb),
                 if (length(problems)) paste0(" -- ", problems[[1]]) else ""))
  expect(all(c("README", "Sample_Metadata", "Animal_Level_Metadata", "Primary_Contrasts",
               "Differential_Proteins", "GSEA_GO_BP", "GSEA_KEGG", "ORA_GO_BP", "EWCE",
               "Secondary_Analyses", "Figure_Source_Map", "Software_Versions") %in%
               wb$sheet),
         "every required sheet is present")
  expect(all(as.integer(wb$n_rows) > 0L), "no sheet is empty")
  expect(all(nchar(wb$sheet) <= 31L), "no sheet name exceeds the Excel limit")
}

cat("\n=== SDRF ===\n")
sdrf_path <- file.path(OUT_ROOT, "pride", "sdrf.tsv")
fs_path <- file.path(OUT_ROOT, "pride", "sdrf_field_status.tsv")
if (!file.exists(sdrf_path) || !file.exists(fs_path)) {
  expect(FALSE, "SDRF and its field-status table exist")
} else {
  sdrf <- rd(sdrf_path); fs <- rd(fs_path)
  expect(nrow(sdrf) == RELEASE_DESIGN_INVARIANTS$n_measurement_records,
         "SDRF has one row per acquisition")
  expect(!anyDuplicated(names(sdrf)), "SDRF has no duplicate column names")
  expect(all(nzchar(names(sdrf))), "SDRF has no unnamed columns")
  expect(!anyDuplicated(sdrf[["comment[data file]"]]),
         "SDRF data file names are unique")
  expect(all(fs$status %in% release_sdrf_field_statuses),
         "every SDRF field carries a sanctioned status")
  req_missing <- fs$sdrf_field[fs$status == "MISSING_REQUIRED_METADATA"]
  expect(all(vapply(intersect(req_missing, names(sdrf)),
                    function(cn) all(sdrf[[cn]] == "not available"), logical(1))),
         sprintf("every unfillable required field is explicitly 'not available' (%d fields)",
                 length(req_missing)))
  pr_path <- file.path(OUT_ROOT, "pride", "pride_readiness.tsv")
  if (file.exists(pr_path)) {
    pr <- rd(pr_path)
    status <- pr$value[pr$key == "pride_status"]
    expect(status %in% release_pride_statuses,
           sprintf("PRIDE status is a sanctioned value (%s)", status))
    expect(!(status == "PRIDE_READY" && length(req_missing) > 0L),
           "PRIDE readiness is not overstated while required fields are missing")
  }
}

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Publication release manifest contracts failed: %d", failures), call. = FALSE)
}
cat("All publication release manifest contracts hold.\n")
