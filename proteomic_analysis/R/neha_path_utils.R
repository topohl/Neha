# DEPRECATED COMPATIBILITY SHIM -- do not use in new code.
#
# The active implementation moved to R/project_path_utils.R on 2026-08-28, when the collaborator
# name was removed from this project's active branding and identifiers.
#
# This shim exists for exactly one reason. The manuscript figure-revision snapshots in
# 06_manuscript_figure_revision/ are byte-identical, hash-locked provenance (see that folder's
# SCRIPT_PROVENANCE.csv), so they cannot be edited. 01_panelC_and_suppB_pca.R does
#
#     source(file.path(REPO_ROOT, "R", "analysis_labels.R"))
#     source(file.path(REPO_ROOT, "R", "protigy_input_utils.R"))
#     source(file.path(REPO_ROOT, "R", "neha_path_utils.R"))        # <- this file
#     source(file.path(REPO_ROOT, "R", "pca_animal_level_utils.R"))
#     ...
#     validated <- validate_neha_pca_animal_input(parsed, expected_n = 3L)
#     prepared  <- prepare_neha_animal_pca(validated$expression_matrix, ...)
#
# and its runnable twin under 03_output/reviewer_revision_animal_level_20260827/scripts/ reaches
# into this repository for those three names. Removing them outright would silently break a re-run
# of that snapshot. Note the ordering above: this file is sourced BEFORE pca_animal_level_utils.R,
# so it must supply the path primitives that file's guard requires, while the two function
# aliases resolve lazily at call time (by then the real definitions are loaded).
#
# Nothing in the ACTIVE pipeline sources this file. Every active caller, test and workflow uses
# R/project_path_utils.R and the un-prefixed names directly.
#
# REMOVAL CRITERION: delete this file once 06_manuscript_figure_revision/ is either made
# self-contained (the deferred follow-up noted in its README) or formally retired. It carries no
# analysis logic, so deleting it can never change a numerical result.

# Locate the sibling implementation. When a file is sourced, the calling frame holds `ofile`.
local({
  dir <- NULL
  for (i in rev(seq_len(sys.nframe()))) {
    of <- get0("ofile", envir = sys.frame(i), inherits = FALSE)
    if (!is.null(of) && is.character(of) && length(of) == 1L && nzchar(of)) {
      dir <- dirname(normalizePath(of, winslash = "/", mustWork = FALSE))
      break
    }
  }
  if (is.null(dir)) dir <- file.path(normalizePath(getwd(), winslash = "/", mustWork = FALSE), "R")
  target <- file.path(dir, "project_path_utils.R")
  if (!file.exists(target)) {
    stop("Deprecated shim R/neha_path_utils.R cannot locate its replacement R/project_path_utils.R at: ",
         target, call. = FALSE)
  }
  sys.source(target, envir = globalenv())
})

.neha_deprecated <- function(old, new) {
  warning(sprintf(
    "%s() is deprecated and kept only for the frozen 06_manuscript_figure_revision snapshots; use %s().",
    old, new
  ), call. = FALSE)
}

validate_neha_pca_animal_input <- function(...) {
  .neha_deprecated("validate_neha_pca_animal_input", "validate_pca_animal_input")
  if (!exists("validate_pca_animal_input", mode = "function")) {
    stop("R/pca_animal_level_utils.R must be sourced before calling this deprecated alias.", call. = FALSE)
  }
  validate_pca_animal_input(...)
}

prepare_neha_animal_pca <- function(...) {
  .neha_deprecated("prepare_neha_animal_pca", "prepare_animal_pca")
  if (!exists("prepare_animal_pca", mode = "function")) {
    stop("R/pca_animal_level_utils.R must be sourced before calling this deprecated alias.", call. = FALSE)
  }
  prepare_animal_pca(...)
}
