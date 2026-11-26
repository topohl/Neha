#!/usr/bin/env Rscript

# Publication release: runs stages 01-12 in order and stops at the first failure. Run
# 13_validate_release.R afterwards with the same output root.
#
# From proteomic_analysis/, build into scratch first:
#   PROTEOMICS_RELEASE_OUTPUT_ROOT=<scratch dir> Rscript 07_publication_release/run_release.R
#
# Without PROTEOMICS_RELEASE_OUTPUT_ROOT the release goes to the shared-drive root
# 03_output/release_within_compartment_20261126b. Any output root on the shared drive also
# needs PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE=true and no uncommitted or untracked changes
# under proteomic_analysis/ (or PROTEOMICS_RELEASE_ALLOW_DIRTY=true). Output roots that
# overlap 01_input, 02_data, 99_historical, the migration backup, an analysis folder or a
# figure rebuild under 03_output or a frozen earlier release (03_output/publication_release,
# 03_output/release_within_compartment_20261125, 03_output/release_within_compartment_20261126)
# are refused.
#
# The locked input hashes are checked and the build registry is reset before stage 01.
# Each stage runs in its own environment and can also be run alone with Rscript.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
LAYER <- file.path(REPO_ROOT, "07_publication_release")
source(file.path(LAYER, "R", "release_validation.R"))

OUT_ROOT <- release_prepare_output_root()

# Remove the previous stage 13 reports (see release_clear_previous_validation_reports()).
cleared_validation_reports <- release_clear_previous_validation_reports(OUT_ROOT)

cat("\n")
cat("================================================================\n")
cat(" Associative Memory Proteomics -- publication release build\n")
cat("================================================================\n")
cat("repository   : ", REPO_ROOT, "\n", sep = "")
cat("git commit   : ", release_git_commit(REPO_ROOT), "\n", sep = "")
cat("data root    : ", release_data_root(), "\n", sep = "")
cat("project root : ", release_project_root(), "\n", sep = "")
cat("release root : ", OUT_ROOT, "\n", sep = "")
cat("statistics   : ", RELEASE_STATISTICS_RUN$id, " (", RELEASE_PRIMARY_MODEL$label, ")\n",
    sep = "")
cat("started      : ", release_timestamp_utc(), "\n", sep = "")
if (length(cleared_validation_reports)) {
  cat("cleared      : previous Stage 13 validation reports (",
      length(cleared_validation_reports), ")\n", sep = "")
}

locked <- release_verify_locked_artefacts()
cat("\nlocked artifacts:\n")
for (i in seq_len(nrow(locked))) {
  cat("  ", ifelse(locked$matches[i], "[ok]  ", "[FAIL]"), " ", locked$artefact[i], " ",
      substr(locked$observed_sha256[i], 1, 16), "\n", sep = "")
}
if (any(!locked$matches)) {
  stop("Locked artifact hash mismatch; refusing to build a release from changed inputs.",
       call. = FALSE)
}

release_registry_reset(OUT_ROOT)

STAGES <- c(
  "01_build_sample_metadata.R",
  "02_build_contrast_manifest.R",
  "03_build_processed_data_exports.R",
  "04_build_differential_results.R",
  "05_build_enrichment_exports.R",
  "06_build_figure_source_data.R",
  "07_build_editor_source_workbook.R",
  "08_build_editor_changelog.R",
  "09_build_pride_sdrf.R",
  "10_build_provenance.R",
  "11_build_readme_and_dictionary.R",
  "12_build_release_manifest.R"
)

results <- data.frame(stage = character(0), status = character(0),
                      seconds = numeric(0), message = character(0),
                      stringsAsFactors = FALSE, check.names = FALSE)

for (stage in STAGES) {
  path <- file.path(LAYER, stage)
  if (!file.exists(path)) {
    stop("Release stage missing: ", path, call. = FALSE)
  }
  started <- Sys.time()
  outcome <- tryCatch({
    env <- new.env(parent = globalenv())
    assign("here", LAYER, envir = env)
    sys.source(path, envir = env, keep.source = FALSE)
    list(status = "PASS", message = "")
  }, error = function(e) list(status = "FAIL", message = conditionMessage(e)))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  results <- rbind(results, data.frame(stage = stage, status = outcome$status,
                                       seconds = round(elapsed, 1),
                                       message = outcome$message,
                                       stringsAsFactors = FALSE, check.names = FALSE))
  if (outcome$status == "FAIL") {
    cat("\n[FAIL] ", stage, ": ", outcome$message, "\n", sep = "")
    break
  }
}

cat("\n================================================================\n")
cat(" build summary\n")
cat("================================================================\n")
for (i in seq_len(nrow(results))) {
  cat(sprintf("  %-6s %-40s %6.1fs\n", results$status[i], results$stage[i],
              results$seconds[i]))
  if (nzchar(results$message[i])) cat("         ", results$message[i], "\n", sep = "")
}
n_fail <- sum(results$status == "FAIL")
n_run <- nrow(results)
cat(sprintf("\n  %d/%d stages PASS, %d FAIL, %d not reached\n",
            n_run - n_fail, length(STAGES), n_fail, length(STAGES) - n_run))
cat("  release root: ", OUT_ROOT, "\n", sep = "")
cat("  finished    : ", release_timestamp_utc(), "\n\n", sep = "")

if (n_fail > 0L) {
  quit(save = "no", status = 1L)
}
cat("Next: Rscript 07_publication_release/13_validate_release.R\n\n")
