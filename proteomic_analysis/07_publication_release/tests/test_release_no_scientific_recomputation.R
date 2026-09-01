#!/usr/bin/env Rscript

# Contract: the publication release layer READS validated results. It must never fit a
# model, run an enrichment, bootstrap a cell-type test, or recompute a PCA.
#
# This is a static scan of the layer's own source, so it needs neither the shared drive nor
# a built release and always runs -- including in CI.
#
# The scanner's own pattern table names the prohibited calls as strings. Strings do not
# match, because the pattern requires an opening parenthesis after the identifier; but
# R/release_validation.R and 13_validate_release.R are excluded anyway so that adding a
# call name to the table can never be mistaken for making the call.

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
  release_validate_output_root(file.path(data_root, "03_output", "publication_release"),
                               data_root),
  error = function(e) e), "error")
expect(accepted, "accepts the intended release root 03_output/publication_release")
accepted_scratch <- !inherits(tryCatch(
  release_validate_output_root(tempdir(), data_root), error = function(e) e), "error")
expect(accepted_scratch, "accepts a scratch root")

cat("\n=== writing to the shared drive requires an explicit opt-in ===\n")

shared_root <- file.path(data_root, "03_output", "publication_release")
saved <- Sys.getenv("PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE", unset = NA)
Sys.unsetenv("PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE")
blocked <- inherits(tryCatch(release_prepare_output_root(shared_root, data_root),
                             error = function(e) e), "error")
expect(blocked, "a shared-drive release root is refused without PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE")
expect(!dir.exists(shared_root),
       "the refused shared-drive root was not even created as an empty directory")

Sys.setenv(PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE = "true")
opted_in <- !inherits(tryCatch(release_assert_shared_drive_opt_in(shared_root, data_root),
                               error = function(e) e), "error")
expect(opted_in, "the opt-in permits a deliberate shared-drive release")
if (is.na(saved)) {
  Sys.unsetenv("PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE")
} else {
  Sys.setenv(PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE = saved)
}

scratch_dir <- file.path(tempdir(), "amp_release_guard_probe")
scratch_ok <- !inherits(tryCatch(release_prepare_output_root(scratch_dir, data_root),
                                 error = function(e) e), "error")
expect(scratch_ok, "a scratch root needs no opt-in")
unlink(scratch_dir, recursive = TRUE)

expect(release_is_shared_drive_root(shared_root, data_root) &&
         !release_is_shared_drive_root(tempdir(), data_root),
       "shared-drive detection distinguishes the data root from scratch")

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Publication release recomputation contracts failed: %d", failures),
       call. = FALSE)
}
cat("All publication release recomputation contracts hold.\n")
