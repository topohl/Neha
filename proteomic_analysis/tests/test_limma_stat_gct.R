#!/usr/bin/env Rscript

# Tests for 01_preprocessing/02b_limma_stat_gct.r.
#   1. pooled mode reproduces the 2026-08-24 ProTigy export field by field;
#   2. within_compartment mode fits one model per sample class (residual df 8), writes a GCT
#      that the split stage accepts, and reproduces the locked primary statistics GCT byte for
#      byte (the GCT embeds no paths or timestamps). That last check is what ties the script's
#      eBayes and BH settings to the published statistics. The model fields of the run manifest
#      must match the canonical run; a different R or limma version is reported, not failed.
# Writes only to a temporary folder. Skips when the shared drive or limma is absent.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("tests", "test_limma_stat_gct.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "project_path_utils.R"))
source(file.path(repo_root, "R", "protigy_stat_gct_utils.R"))

animal_root <- file.path(PROTEOMICS_DATA_ROOT, "02_data", "animal_level")
protigy_export <- file.path(animal_root, "stat_results_for_ssGSEA_neha_proteome.gct")
input_gct <- file.path(animal_root, "input_gct", "neha_protigy_input_animal_level_primary.gct")
if (!file.exists(protigy_export) || !file.exists(input_gct)) {
  cat("Shared data root not reachable; skipping limma stat GCT checks.\n")
  quit(save = "no", status = 0L)
}
if (!requireNamespace("limma", quietly = TRUE) || !requireNamespace("digest", quietly = TRUE)) {
  cat("limma or digest not installed; skipping limma stat GCT checks.\n")
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

run_stage <- function(model, output) {
  rscript <- file.path(R.home("bin"), "Rscript")
  # Outputs stay in this session's temp folder. FORCE is false, so if the output override
  # were ignored the stage would refuse to replace its file on S:.
  stopifnot(project_path_is_within(output, tempdir()), !file.exists(output))
  # system2(env = ) is ignored on Windows, so set the variables for the child process here.
  # FORCE and the input are set here too, so neither is inherited from the outer session.
  vars <- c(PROTEOMICS_LIMMA_VARIANCE_MODEL = model, PROTEOMICS_LIMMA_OUTPUT_GCT = output,
            PROTEOMICS_LIMMA_FORCE = "false", PROTEOMICS_LIMMA_INPUT_GCT = input_gct)
  old <- Sys.getenv(names(vars), unset = NA, names = TRUE)
  do.call(Sys.setenv, as.list(vars))
  on.exit({
    for (n in names(old)) if (is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, stats::setNames(list(old[[n]]), n))
  })
  status <- system2(rscript, shQuote(file.path(repo_root, "01_preprocessing", "02b_limma_stat_gct.r")),
                    stdout = FALSE, stderr = FALSE)
  identical(as.integer(status), 0L) && file.exists(output)
}
read_fields <- function(path) utils::read.delim(path, skip = 2, check.names = FALSE, colClasses = "character")

scratch <- file.path(tempdir(), "limma_stat_gct_test")
dir.create(scratch, showWarnings = FALSE)

cat("pooled mode vs ProTigy export\n")
pooled_path <- file.path(scratch, "pooled.gct")
expect(run_stage("pooled", pooled_path), "pooled run completes")
if (file.exists(pooled_path)) {
  ours <- read_fields(pooled_path)
  theirs <- read_fields(protigy_export)
  expect(identical(readLines(pooled_path, 2L), readLines(protigy_export, 2L)), "GCT version and dimension lines match")
  expect(identical(names(ours), names(theirs)), "field names and order match")
  expect(identical(ours$id, theirs$id), "protein ids match")
  significant <- grep("^significant[.]", names(ours), value = TRUE)
  expect(all(vapply(significant, function(n) identical(ours[[n]], theirs[[n]]), logical(1))), "significance calls match")
  numeric_fields <- setdiff(names(ours), c("id", "Description", significant))
  worst <- max(vapply(numeric_fields, function(n) max(abs(as.numeric(ours[[n]]) - as.numeric(theirs[[n]]))), numeric(1)))
  expect(worst < 1e-8, paste0("all numeric fields agree (max abs difference ", signif(worst, 3), ")"))
}

cat("within_compartment mode\n")
within_path <- file.path(scratch, "within.gct")
expect(run_stage("within_compartment", within_path), "within_compartment run completes")
if (file.exists(within_path)) {
  manifest <- utils::read.csv(sub("[.]gct$", "_run_manifest.csv", within_path), stringsAsFactors = FALSE)
  expect(nrow(manifest) == 12L, "12 comparisons")
  expect(all(manifest$n_columns_in_model == 12L) && all(manifest$df_residual == 8L), "one model per sample class: 12 columns, residual df 8")
  gct <- read_protigy_stat_gct(within_path, strict_primary = TRUE)
  contract <- tryCatch(validate_stat_gct_contract(gct, expected_n_proteins = 5349L), error = function(e) e)
  expect(!inherits(contract, "error"), "split-stage check accepts the GCT")
  ours <- read_fields(within_path)
  theirs <- read_fields(protigy_export)
  logfc <- grep("^logFC[.]", names(ours), value = TRUE)
  expect(all(vapply(logfc, function(n) max(abs(as.numeric(ours[[n]]) - as.numeric(theirs[[n]]))) < 1e-12, logical(1))),
         "effect sizes are the same as in the pooled model")

  # The published statistics: the release locks this GCT by hash (RELEASE_LOCKED_ARTEFACTS).
  release_env <- new.env()
  sys.source(file.path(repo_root, "07_publication_release", "R", "release_utils.R"), envir = release_env)
  sys.source(file.path(repo_root, "07_publication_release", "R", "release_validation.R"), envir = release_env)
  locked_sha <- release_env$RELEASE_LOCKED_ARTEFACTS$stat_gct$sha256
  canonical_gct <- file.path(animal_root, "within_compartment_20261124", "stat_results_within_compartment.gct")
  ours_sha <- project_file_sha256(within_path)
  expect(identical(ours_sha, locked_sha),
         paste0("the script reproduces the locked statistics GCT byte for byte (", substr(ours_sha, 1, 16), ")"))
  if (file.exists(canonical_gct)) {
    expect(identical(project_file_sha256(canonical_gct), locked_sha), "the canonical statistics GCT on the shared drive still has the locked hash")
    canonical_manifest <- utils::read.csv(sub("[.]gct$", "_run_manifest.csv", canonical_gct), stringsAsFactors = FALSE)
    same_cols <- c("canonical_comparison", "n_columns_in_model", "df_residual", "df_prior_median", "s2_prior",
                   "variance_model", "n_fdr_proteins", "input_gct_sha256", "output_gct_sha256")
    expect(identical(manifest[same_cols], canonical_manifest[same_cols]),
           "the run manifest matches the canonical one apart from paths, timestamp and software versions")
    # Another R or limma build is reported, not failed: the byte-for-byte check above already
    # ties the script to the locked statistics.
    versions <- c("r_version", "limma_version")
    ours_versions <- manifest[versions]
    canonical_versions <- canonical_manifest[versions]
    for (col in versions) {
      if (!identical(ours_versions[[col]], canonical_versions[[col]])) {
        cat(sprintf("  [note] %s differs from the canonical run: %s here, %s in the canonical run\n", col,
                    paste(unique(ours_versions[[col]]), collapse = "; "),
                    paste(unique(canonical_versions[[col]]), collapse = "; ")))
      }
    }
  }
}

if (failures > 0L) {
  cat(failures, "limma stat GCT check(s) failed.\n")
  quit(save = "no", status = 1L)
}
cat("All limma stat GCT checks passed.\n")
