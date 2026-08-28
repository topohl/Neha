#!/usr/bin/env Rscript

# Data-integrity contracts for the validated animal-level artefacts on shared storage.
#
# Added 2026-08-26 alongside the input/data/output restructure. Folder moves cannot alter
# file contents, but they DO invalidate the absolute paths recorded inside index/manifest
# CSVs -- and those are actively file.exists()- and hash-verified at read time. These
# contracts catch a broken provenance chain immediately instead of at the next pipeline run.
#
# Skips cleanly when the shared drive is unreachable, so it is safe in CI.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("tests", "test_data_integrity.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "project_path_utils.R"))
source(file.path(repo_root, "R", "animal_level_enrichment_utils.R"))

if (!requireNamespace("digest", quietly = TRUE)) {
  cat("digest not installed; skipping data-integrity contracts.\n")
  quit(save = "no", status = 0L)
}

failures <- 0L
expect <- function(condition, message) {
  if (!isTRUE(condition)) {
    failures <<- failures + 1L
    cat("  [FAIL]", message, "\n")
  } else {
    cat("  [ok]  ", message, "\n")
  }
}

option_or_env <- function(option_name, env_name, default) {
  value <- getOption(option_name)
  if (!is.null(value) && nzchar(trimws(as.character(value)))) return(as.character(value))
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(trimws(value))) return(value)
  default
}

data_root <- option_or_env(
  "proteomics.project_data_root", "PROTEOMICS_PROJECT_DATA_ROOT",
  "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
)

PRIMARY_GCT_SHA256 <- "f12cf99e1bfb7c17bbf56bffb6783e924698bce5d5533a8e312bc4bbb733bbb3"
gct_path <- file.path(data_root, "02_data", "animal_level", "input_gct",
                      "neha_protigy_input_animal_level_primary.gct")
mapped_root <- file.path(data_root, "02_data", "animal_level", "mapped")
mapped_index <- file.path(mapped_root, "indexMappedComparisons.csv")

if (!file.exists(gct_path)) {
  cat("Shared data root not reachable; skipping data-integrity contracts.\n")
  cat("  expected:", gct_path, "\n")
  cat("  set PROTEOMICS_PROJECT_DATA_ROOT to run them.\n")
  quit(save = "no", status = 0L)
}

cat("=== validated animal-level GCT ===\n")
observed <- digest::digest(file = gct_path, algo = "sha256")
expect(identical(observed, PRIMARY_GCT_SHA256),
       sprintf("primary GCT SHA256 matches the locked contract (%s)", substr(observed, 1, 16)))

cat("\n=== recorded provenance still resolves and verifies ===\n")
if (!file.exists(mapped_index)) {
  expect(FALSE, paste("mapped index is missing:", mapped_index))
} else {
  ix <- tryCatch(
    read_enrichment_mapped_index(mapped_index, mapped_root, verify_hashes = TRUE),
    error = function(e) e
  )
  if (inherits(ix, "error")) {
    expect(FALSE, paste("mapped index failed hash/path verification:", conditionMessage(ix)))
  } else {
    expect(nrow(ix) == 12L,
           sprintf("mapped index verifies with hashes on and covers 12 comparisons (got %d)", nrow(ix)))
  }
}

cat("\n=== index/manifest CSVs: every recorded path+hash pair agrees ===\n")
idx_files <- list.files(file.path(data_root, c("02_data", "03_output")),
                        pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
idx_files <- idx_files[vapply(idx_files, function(f) {
  hdr <- tryCatch(names(utils::read.csv(f, nrows = 1L, check.names = FALSE)),
                  error = function(e) character())
  any(grepl("sha256", hdr, ignore.case = TRUE)) && any(grepl("path", hdr, ignore.case = TRUE))
}, logical(1))]

checked <- 0L
problems <- character()
for (f in idx_files) {
  d <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d)) next
  for (pc in grep("path$", names(d), ignore.case = TRUE, value = TRUE)) {
    stem <- sub("_?path$", "", pc, ignore.case = TRUE)
    cand <- grep("sha256", names(d), ignore.case = TRUE, value = TRUE)
    sc <- cand[grepl(stem, cand, ignore.case = TRUE)]
    if (!length(sc)) next
    sc <- sc[[1]]
    for (i in seq_len(nrow(d))) {
      p <- as.character(d[[pc]][i]); s <- tolower(as.character(d[[sc]][i]))
      if (is.na(p) || !nzchar(p) || is.na(s) || !nzchar(s)) next
      if (!grepl("clusterProfiler", p, fixed = TRUE)) next
      if (!file.exists(p)) { problems <- c(problems, paste("missing:", p)); next }
      if (dir.exists(p)) next
      checked <- checked + 1L
      if (!identical(digest::digest(file = p, algo = "sha256"), s)) {
        problems <- c(problems, paste("hash mismatch:", p, "recorded in", basename(f), sc))
      }
    }
  }
}
cat("    index/manifest CSVs scanned:", length(idx_files), "| path+hash pairs checked:", checked, "\n")
expect(checked > 0L, "found recorded path+hash pairs to verify")
expect(length(problems) == 0L, "all recorded paths resolve and all recorded hashes match")
if (length(problems)) for (p in unique(problems)) cat("      -", p, "\n")

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Data-integrity contracts failed: %d", failures), call. = FALSE)
}
cat("All data-integrity contracts hold.\n")
