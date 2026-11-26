# Deprecated: use R/project_path_utils.R (renamed 2026-08-28).
#
# Kept only for the frozen scripts in 06_manuscript_figure_revision/, which are hash-locked
# (see its SCRIPT_PROVENANCE.csv). Each 01_panelC_and_suppB_pca.R there, and its executed
# copy under 03_output/, sources this file before R/pca_animal_level_utils.R and calls
# validate_neha_pca_animal_input() and prepare_neha_animal_pca(). So this file loads
# R/project_path_utils.R, which the PCA helpers need, and defines those two names as
# aliases that warn and call the current functions.
#
# No active stage sources it (tests/test_deprecated_path_utils_shim.R checks this). Delete it
# once 06_manuscript_figure_revision/ no longer needs it; it holds no analysis logic.

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
