#!/usr/bin/env Rscript

# Checks: every file in the package is in the release manifest, every listed hash matches,
# repository sources are named by paths relative to proteomic_analysis/, and the workbook
# and SDRF meet their own rules.

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
  cat("Skipping publication release manifest checks.\n")
  quit(save = "no", status = 0L)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  cat("digest not installed; skipping publication release manifest checks.\n")
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

cat("\n=== repository files are named by relative paths ===\n")
repo_prefix <- tolower(paste0(sub("/+$", "", repo_root), "/"))
text_files <- on_disk[grepl("[.](tsv|md|txt|csv)$|[.]tsv[.]gz$", on_disk)]
local_hits <- text_files[vapply(text_files, function(f) {
  txt <- gsub("\\\\", "/", readLines(file.path(OUT_ROOT, f), warn = FALSE), useBytes = TRUE)
  txt <- tryCatch(tolower(txt), error = function(e) txt)
  any(grepl(repo_prefix, txt, fixed = TRUE, useBytes = TRUE))
}, logical(1))]
expect(length(local_hits) == 0L,
       sprintf("no release file names the local repository folder (%d text files%s)",
               length(text_files),
               if (length(local_hits)) paste0("; found in: ",
                                              paste(head(local_hits, 3), collapse = "; "))
               else ""))
src <- unique(unlist(strsplit(man$source_artifact[!is.na(man$source_artifact)], ";",
                              fixed = TRUE)))
repo_src <- src[release_is_repo_path(src, repo_root)]
missing_src <- repo_src[!file.exists(release_repo_resolve(repo_src, repo_root))]
expect(length(repo_src) > 0L && length(missing_src) == 0L,
       sprintf("repository sources in the manifest resolve under proteomic_analysis/ (%d%s)",
               length(repo_src),
               if (length(missing_src)) paste0("; not found: ",
                                               paste(head(missing_src, 3), collapse = "; "))
               else ""))

cat("\n=== the expected package structure is present ===\n")
for (f in c("README_DATA.md", "metadata/sample_metadata.tsv",
            "metadata/animal_level_sample_metadata.tsv",
            "metadata/experimenter_metadata.tsv",
            "metadata/sample_preparation_protocol.tsv",
            "metadata/primary_contrast_manifest.tsv",
            "metadata/secondary_analysis_manifest.tsv", "metadata/data_dictionary.tsv",
            "processed_data/protein_abundance_measurement_level.tsv.gz",
            "processed_data/protein_abundance_animal_level.tsv.gz",
            "processed_data/protein_feature_annotation.tsv.gz",
            "differential_analysis/primary_differential_proteins.tsv.gz",
            "differential_analysis/significant_proteins.tsv",
            "enrichment/primary_GSEA_GO_BP.tsv.gz", "enrichment/primary_GSEA_GO_CC.tsv.gz",
            "enrichment/primary_GSEA_KEGG.tsv.gz", "enrichment/primary_ORA_GO_BP.tsv.gz",
            "enrichment/primary_ORA_GO_CC.tsv.gz",
            "enrichment/GSEA_log2FC_sensitivity.tsv.gz", "enrichment/primary_EWCE.tsv.gz",
            "enrichment/gsea_untested_gene_sets.tsv",
            "enrichment/significant_gene_sets.tsv.gz",
            RELEASE_WORKBOOK_REL,
            "editor_source_data/workbook_sheet_inventory.tsv",
            "provenance/REVISION_PROTEOMICS_DATA_CHANGELOG.md",
            "provenance/figure_source_map.tsv",
            "pride/sdrf.tsv", "pride/SDRF_MISSING_METADATA.md", "pride/README_PRIDE.md",
            "pride/SAMPLE_PREPARATION_PROTOCOL.md",
            "provenance/data_lineage.tsv", "provenance/software_versions.tsv",
            "provenance/analysis_parameters.tsv",
            "provenance/sessionInfo_release.txt",
            "provenance/UPSTREAM_PREPROCESSING_GAP.md")) {
  expect(f %in% man$relative_path, sprintf("manifest lists %s", f))
}

cat("\n=== Supplementary Table 2 workbook ===\n")
inv_path <- file.path(OUT_ROOT, "editor_source_data", "workbook_sheet_inventory.tsv")
xlsx_path <- file.path(OUT_ROOT, RELEASE_WORKBOOK_REL)
if (!file.exists(inv_path) || !file.exists(xlsx_path)) {
  expect(FALSE, "the workbook and workbook_sheet_inventory.tsv exist")
} else {
  wb <- utils::read.delim(inv_path, sep = "\t", quote = "", stringsAsFactors = FALSE,
                          check.names = FALSE, colClasses = "character",
                          na.strings = character(0))
  sheet_cols <- stats::setNames(strsplit(wb$columns, ";", fixed = TRUE), wb$sheet)
  problems <- release_check_workbook_hygiene(sheet_cols)
  expect(length(problems) == 0L,
         sprintf("no duplicate or unnamed workbook columns across %d sheets%s", nrow(wb),
                 if (length(problems)) paste0(" -- ", problems[[1]]) else ""))
  display_cols <- stats::setNames(strsplit(wb$display_columns, ";", fixed = TRUE), wb$sheet)
  problems_display <- release_check_workbook_hygiene(display_cols)
  expect(length(problems_display) == 0L,
         sprintf("no duplicate or unnamed DISPLAY headers%s",
                 if (length(problems_display)) paste0(" -- ", problems_display[[1]]) else ""))
  expect(all(lengths(display_cols) == lengths(sheet_cols)),
         "each sheet has one display header per machine column")
  expect(!any(grepl("_", unlist(display_cols, use.names = FALSE), fixed = TRUE)),
         "no display header leaks a snake_case machine name")
  expect(all(nzchar(wb$title)), "every sheet has a title")
  expect(identical(wb$sheet_key, RELEASE_WORKBOOK_LAYOUT$sheet_key),
         "the sheets are README, Column definitions, then 2a to 2q in the expected order")
  expect(identical(wb$supplementary_table, RELEASE_WORKBOOK_LAYOUT$letter),
         sprintf("lettered 2a to 2q; README and Column definitions carry no letter (%s)",
                 paste(wb$supplementary_table[nzchar(wb$supplementary_table)], collapse = " ")))
  expect(!any(c("GSEA_log2FC_sensitivity", "Metadata_Field_Status") %in% wb$sheet_key),
         "the effect-size-ranked GSEA and the metadata field status stay out of the workbook")
  keyed <- vapply(sheet_cols, function(cs) "canonical_comparison" %in% cs, logical(1))
  expect(sum(keyed) > 0L &&
           all(vapply(sheet_cols[keyed], function(cs) {
             cs[[1]] == "comparison_label" && cs[[length(cs)]] == "canonical_comparison"
           }, logical(1))),
         sprintf("comparison-keyed sheets show the comparison first and its key last (%d sheets)",
                 sum(keyed)))
  internal_cols <- grep(RELEASE_WORKBOOK_INTERNAL_COLUMNS, unlist(sheet_cols), value = TRUE)
  expect(length(internal_cols) == 0L,
         sprintf("no contrast, alias, path, hash or curation column in the workbook%s",
                 if (length(internal_cols)) paste0(" -- ", paste(unique(internal_cols),
                                                                 collapse = ", ")) else ""))
  expect(all(as.integer(wb$n_rows) > 0L), "no sheet is empty")
  expect(all(nchar(wb$sheet) <= 31L), "no sheet name exceeds the Excel limit")

  wbx <- release_read_workbook_xml(xlsx_path)
  expect(identical(wbx$tabs, c("Contents", wb$sheet)),
         "the workbook tabs are Contents, then the inventory sheets in order")
  expect(!any(vapply(wbx$sheets, `[[`, logical(1), "inline_strings")),
         "every text cell is a shared string, so the scan below sees all text")
  hits <- release_internal_text_hits(c(wbx$tabs, wbx$shared_strings,
                                       unlist(lapply(wbx$sheets, `[[`, "formulas"))))
  expect(nrow(hits) == 0L,
         sprintf("no cell holds a drive path, hash, script or folder name, 'Figure3' or Fabian%s",
                 if (nrow(hits)) paste0(" -- ", hits$kind[1], ": ", substr(hits$text[1], 1, 80))
                 else ""))
  expect(all(vapply(wbx$sheets, function(s) identical(s$frozen_at, "B5"), logical(1))),
         "every sheet is frozen below the row-4 header and right of column 1")
  british <- release_british_spelling_hits(wbx$shared_strings,
                                           release_workbook_term_names(xlsx_path, wb))
  expect(length(british) == 0L,
         sprintf("no workbook text uses British spelling, GO and KEGG term names of 2f to 2l excepted%s",
                 if (length(british)) paste0(" -- ", substr(british[1], 1, 80)) else ""))
}
us_docs <- file.path(OUT_ROOT, c("README_DATA.md", "metadata/data_dictionary.tsv",
                                 "pride/README_PRIDE.md"))
us_doc_hits <- unlist(lapply(us_docs[file.exists(us_docs)], function(f) {
  release_british_spelling_hits(readLines(f, warn = FALSE, encoding = "UTF-8"))
}))
expect(all(file.exists(us_docs)) && length(us_doc_hits) == 0L,
       sprintf("README_DATA.md, data_dictionary.tsv and README_PRIDE.md use no British spelling%s",
               if (length(us_doc_hits)) paste0(" -- ", substr(us_doc_hits[1], 1, 80)) else ""))

cat("\n=== panel source data ===\n")
fsd_dir <- file.path(OUT_ROOT, "editor_source_data", "figure_source_data")
fsd_files <- list.files(fsd_dir)
expect(length(fsd_files) > 0L && all(grepl("[.]csv$", fsd_files)),
       sprintf("the panel source data are CSV files only (%d files)", length(fsd_files)))
supp_a_path <- file.path(fsd_dir, "SuppFig4A_source_data.csv")
sm_path <- file.path(OUT_ROOT, "metadata", "sample_metadata.tsv")
if (!file.exists(supp_a_path) || !file.exists(sm_path)) {
  expect(FALSE, "SuppFig4A_source_data.csv and the sample metadata exist")
} else {
  supp_a <- release_read_published_csv(supp_a_path)
  sm <- rd(sm_path)
  expect(nrow(supp_a) == RELEASE_DESIGN_INVARIANTS$n_measurement_records &&
           setequal(supp_a[["Sample ID"]], sm$sample_id) &&
           !anyDuplicated(supp_a[["Sample ID"]]),
         "Supp. Fig. 4 A holds the 96 acquisitions of this project, by their table 2a Sample IDs")
  foreign <- unlist(lapply(file.path(fsd_dir, fsd_files), function(p) {
    d <- release_read_published_csv(p)
    if (any(grepl("RFID|Exp9|ActivityIndex|AnimalPos", c(names(d), unlist(d))))) basename(p)
  }))
  expect(length(foreign) == 0L,
         sprintf("no panel source-data file carries data of another project%s",
                 if (length(foreign)) paste0(" -- ", paste(foreign, collapse = ", ")) else ""))
  fsd_header_list <- lapply(stats::setNames(file.path(fsd_dir, fsd_files), fsd_files),
                            function(p) names(release_read_published_csv(p)))
  bad <- fsd_files[!vapply(fsd_files, function(f) {
    identical(fsd_header_list[[f]], RELEASE_FIGURE_CSV_HEADERS[[f]])
  }, logical(1))]
  inv_path <- file.path(OUT_ROOT, "editor_source_data", "workbook_sheet_inventory.tsv")
  wb_display <- if (file.exists(inv_path)) {
    unlist(strsplit(rd(inv_path)$display_columns, ";", fixed = TRUE))
  } else character(0)
  variants <- release_figure_header_variants(unlist(fsd_header_list), wb_display)
  expect(setequal(fsd_files, names(RELEASE_FIGURE_CSV_HEADERS)) && !length(bad) &&
           length(wb_display) > 0L && !length(variants),
         sprintf("every panel source-data file has the pinned headers, with the workbook display header for every workbook quantity%s",
                 if (length(bad) || length(variants))
                   paste0(" -- ", paste(c(bad, variants), collapse = ", ")) else ""))
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
  duplicate_headers <- unique(names(sdrf)[duplicated(names(sdrf))])
  expect(setequal(duplicate_headers,
                  c("comment[cleavage agent details]", "comment[sdrf template]")) &&
           sum(names(sdrf) == "comment[cleavage agent details]") == 2L &&
           sum(names(sdrf) == "comment[sdrf template]") == 2L,
         "repeated SDRF headers are exactly the two enzymes and two templates")
  expect(all(nzchar(names(sdrf))), "SDRF has no unnamed columns")
  expect(!anyDuplicated(sdrf[["comment[data file]"]]),
         "SDRF data file names are unique")
  expect(all(fs$status %in% release_sdrf_field_statuses),
         "every SDRF field carries a sanctioned status")
  req_missing <- fs$sdrf_field[fs$status == "MISSING_REQUIRED_METADATA"]
  expect(all(vapply(intersect(req_missing, names(sdrf)),
                    function(cn) all(sdrf[[cn]] == "not available"), logical(1))),
         sprintf("every unfillable required field is 'not available' (%d fields)",
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
  stop(sprintf("Publication release manifest checks failed: %d", failures), call. = FALSE)
}
cat("All publication release manifest checks hold.\n")
