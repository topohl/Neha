#!/usr/bin/env Rscript

# Checks: the release layer only reads results. No R file in it fits a model, runs an
# enrichment, bootstraps a cell-type test or computes a PCA. It also checks the output-root
# rules: protected roots are refused, the shared drive needs an opt-in, and an in-place
# rebuild removes only the previous stage 13 reports.
#
# A static scan of the source, so it needs neither the shared drive nor a built release.
# The scanner uses the parse tree, so call names written as strings do not count;
# R/release_validation.R (which holds the list), 13_validate_release.R and this file are
# skipped anyway.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("07_publication_release", "tests", "test_release_no_scientific_recomputation.R")
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

EXEMPT <- c("release_validation.R", "13_validate_release.R",
            "test_release_no_scientific_recomputation.R")

cat("=== publication scripts contain no scientific computation ===\n")

scripts <- list.files(LAYER, pattern = "[.]R$", full.names = TRUE, recursive = TRUE)
builders <- scripts[!basename(scripts) %in% EXEMPT]
expect(length(builders) >= 12L,
       sprintf("found the release builders to scan (%d files)", length(builders)))

hits <- release_scan_prohibited_calls(builders)
expect(nrow(hits) == 0L,
       sprintf("no prohibited computation call in any builder (%d scanned)",
               length(builders)))
if (nrow(hits)) {
  for (i in seq_len(nrow(hits))) {
    cat("      -", basename(hits$file[i]), ":", hits$line[i], " ", hits$call[i], " | ",
        substr(hits$text[i], 1, 90), "\n", sep = "")
  }
}

cat("\n=== the scanner itself is honest ===\n")

tmp <- tempfile(fileext = ".R")
writeLines(c("fit <- lmFit(x, design)", "fit <- eBayes(fit)",
             "res <- gseGO(geneList = g, OntBP = 'BP')", "p <- prcomp(m)",
             "# a comment mentioning eBayes( should not count",
             "label <- 'we did not call enrichGO'"), tmp)
probe <- release_scan_prohibited_calls(tmp)
expect(nrow(probe) >= 4L,
       sprintf("scanner detects real calls in a probe file (%d hits)", nrow(probe)))
expect(!any(grepl("^#", probe$text)),
       "scanner ignores whole-line comments")
unlink(tmp)

tmp2 <- tempfile(fileext = ".R")
writeLines(c("PROHIBITED <- c('lmFit', 'eBayes', 'gseGO')",
             "x <- read_GSEA_table(path)",
             "d <- results[['p.adjust']]",
             "sample_metadata <- build(x)"), tmp2)
probe2 <- release_scan_prohibited_calls(tmp2)
expect(nrow(probe2) == 0L,
       sprintf("scanner does not fire on name strings or lookalike identifiers (%d hits)",
               nrow(probe2)))
unlink(tmp2)

cat("\n=== builders cannot write into a protected canonical root ===\n")

data_root <- release_data_root()
protected <- release_protected_roots(data_root)
expect(length(protected) >= 10L,
       sprintf("protected root list is populated (%d roots)", length(protected)))

for (p in protected) {
  rejected <- inherits(tryCatch(release_validate_output_root(p, data_root),
                                error = function(e) e), "error")
  expect(rejected, sprintf("rejects an output root at %s", sub("^.*/", ".../", p)))
}
for (p in c(file.path(data_root, "02_data", "animal_level", "input_gct"),
            file.path(data_root, "03_output", "enrichment", "anything"),
            data_root,
            file.path(data_root, "03_output"))) {
  rejected <- inherits(tryCatch(release_validate_output_root(p, data_root),
                                error = function(e) e), "error")
  expect(rejected, sprintf("rejects a nested/containing root at %s", sub("^.*/", ".../", p)))
}
accepted <- !inherits(tryCatch(
  release_validate_output_root(file.path(data_root, "03_output",
                                         "release_within_compartment_20261126b"),
                               data_root),
  error = function(e) e), "error")
expect(accepted,
       "accepts the intended release root 03_output/release_within_compartment_20261126b")
saved_out_root <- Sys.getenv("PROTEOMICS_RELEASE_OUTPUT_ROOT", unset = NA)
Sys.unsetenv("PROTEOMICS_RELEASE_OUTPUT_ROOT")
default_root <- release_output_root()
if (!is.na(saved_out_root)) Sys.setenv(PROTEOMICS_RELEASE_OUTPUT_ROOT = saved_out_root)
expect(identical(basename(default_root), "release_within_compartment_20261126b"),
       "the default release root is 03_output/release_within_compartment_20261126b")
frozen_rejected <- inherits(tryCatch(
  release_validate_output_root(file.path(data_root, "03_output", "publication_release"),
                               data_root),
  error = function(e) e), "error")
expect(frozen_rejected, "rejects the frozen pooled-run release root 03_output/publication_release")
first_rejected <- inherits(tryCatch(
  release_validate_output_root(file.path(data_root, "03_output",
                                         "release_within_compartment_20261125"),
                               data_root),
  error = function(e) e), "error")
expect(first_rejected,
       "rejects the frozen first per-compartment release root 03_output/release_within_compartment_20261125")
second_rejected <- inherits(tryCatch(
  release_validate_output_root(file.path(data_root, "03_output",
                                         "release_within_compartment_20261126"),
                               data_root),
  error = function(e) e), "error")
expect(second_rejected,
       "rejects the frozen second per-compartment release root 03_output/release_within_compartment_20261126")
for (rebuild in c("figure_rebuild_within_compartment_20261124",
                  "figure_rebuild_within_compartment_20261125",
                  "figure_rebuild_within_compartment_20261126")) {
  rebuild_rejected <- inherits(tryCatch(
    release_validate_output_root(file.path(data_root, "03_output", rebuild), data_root),
    error = function(e) e), "error")
  expect(rebuild_rejected, paste0("rejects the figure rebuild root 03_output/", rebuild))
}
accepted_scratch <- !inherits(tryCatch(
  release_validate_output_root(tempdir(), data_root), error = function(e) e), "error")
expect(accepted_scratch, "accepts a scratch root")

cat("\n=== writing to the shared drive requires an explicit opt-in ===\n")

# Use a fresh synthetic data root. The canonical publication release may already exist,
# and its existence must not make this write-guard fixture stateful or destructive.
guard_data_root <- tempfile("amp_release_shared_guard_")
shared_root <- file.path(guard_data_root, "03_output",
                         "release_within_compartment_20261126b")
saved <- Sys.getenv("PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE", unset = NA)
Sys.unsetenv("PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE")
blocked <- inherits(tryCatch(release_prepare_output_root(shared_root, guard_data_root),
                             error = function(e) e), "error")
expect(blocked, "a shared-drive release root is refused without PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE")
expect(!dir.exists(shared_root),
       "the refused shared-drive root was not even created as an empty directory")

Sys.setenv(PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE = "true")
saved_dirty <- Sys.getenv("PROTEOMICS_RELEASE_ALLOW_DIRTY", unset = NA)
Sys.unsetenv("PROTEOMICS_RELEASE_ALLOW_DIRTY")
opted_in <- !inherits(tryCatch(release_assert_shared_drive_opt_in(shared_root, guard_data_root,
                                                                  dirty = FALSE),
                               error = function(e) e), "error")
expect(opted_in, "the opt-in permits a deliberate shared-drive release from a clean tree")
dirty_refused <- inherits(tryCatch(release_assert_shared_drive_opt_in(shared_root, guard_data_root,
                                                                      dirty = TRUE),
                                   error = function(e) e), "error")
expect(dirty_refused,
       "a shared-drive release from a dirty tree is refused without PROTEOMICS_RELEASE_ALLOW_DIRTY")
Sys.setenv(PROTEOMICS_RELEASE_ALLOW_DIRTY = "true")
dirty_allowed <- !inherits(tryCatch(release_assert_shared_drive_opt_in(shared_root, guard_data_root,
                                                                       dirty = TRUE),
                                    error = function(e) e), "error")
expect(dirty_allowed, "PROTEOMICS_RELEASE_ALLOW_DIRTY=true permits it")
scratch_dirty_ok <- !inherits(tryCatch(release_assert_shared_drive_opt_in(tempdir(), guard_data_root,
                                                                          dirty = TRUE),
                                       error = function(e) e), "error")
if (is.na(saved_dirty)) {
  Sys.unsetenv("PROTEOMICS_RELEASE_ALLOW_DIRTY")
} else {
  Sys.setenv(PROTEOMICS_RELEASE_ALLOW_DIRTY = saved_dirty)
}
expect(scratch_dirty_ok, "a scratch build from a dirty tree needs no opt-in")
if (is.na(saved)) {
  Sys.unsetenv("PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE")
} else {
  Sys.setenv(PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE = saved)
}
commit <- release_git_commit(repo_root)
dirty_now <- release_git_tree_dirty(repo_root)
expect(identical(commit, "UNKNOWN") ||
         (grepl("^[0-9a-f]{40}(\\+dirty)?$", commit) &&
            identical(grepl("[+]dirty$", commit), isTRUE(dirty_now))),
       sprintf("the recorded commit carries +dirty exactly when the tree differs from it (%s)",
               commit))

scratch_dir <- file.path(tempdir(), "amp_release_guard_probe")
scratch_ok <- !inherits(tryCatch(release_prepare_output_root(scratch_dir, data_root),
                                 error = function(e) e), "error")
expect(scratch_ok, "a scratch root needs no opt-in")
unlink(scratch_dir, recursive = TRUE)

expect(release_is_shared_drive_root(shared_root, guard_data_root) &&
         !release_is_shared_drive_root(tempdir(), guard_data_root),
       "shared-drive detection distinguishes the data root from scratch")

cat("\n=== in-place rebuild clears only prior generated validation reports ===\n")

rebuild_root <- tempfile("amp_release_rebuild_guard_")
dir.create(file.path(rebuild_root, "provenance"), recursive = TRUE)
validation_targets <- file.path(rebuild_root, RELEASE_VALIDATION_REPORT_RELPATHS)
writeLines("old validation", validation_targets[[1]])
writeLines("old validation", validation_targets[[2]])
unrelated <- file.path(rebuild_root, "provenance", "must_remain.txt")
writeLines("unrelated", unrelated)
cleared <- release_clear_previous_validation_reports(rebuild_root)
expect(length(cleared) == 2L && !any(file.exists(validation_targets)),
       "the two Stage 13 reports are cleared before an in-place rebuild")
expect(file.exists(unrelated),
       "an unrelated release file is not removed by validation-report cleanup")
unlink(rebuild_root, recursive = TRUE)

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Publication release recomputation checks failed: %d", failures),
       call. = FALSE)
}
cat("All publication release recomputation checks hold.\n")
