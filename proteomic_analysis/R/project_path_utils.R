# Shared path primitives for the proteomics pipeline.
#
# Before 2026-08-26 each of R/pca_animal_level_utils.R, R/animal_level_enrichment_utils.R,
# R/mapthatprot_animal_level_utils.R and R/ewce_animal_level_utils.R defined its own
# near-identical normalize / is-within / overlap helpers. They had genuinely divergent
# behaviour, which is a latent bug source rather than harmless duplication:
#
#   * normalisation : PCA/EWCE used bare normalizePath(); enrichment/mapthatprot additionally
#                     applied trimws() + as.character() + path.expand().
#   * case handling : PCA/EWCE always lower-cased both sides (case-insensitive on every OS);
#                     enrichment/mapthatprot only lower-cased on Windows (so they were
#                     case-SENSITIVE on Linux/macOS).
#   * direction     : only enrichment checked containment bidirectionally.
#
# This file resolves those three differences deliberately, in the safer direction each time:
# the fuller normalisation, always-case-insensitive comparison, and bidirectional overlap.
# Comparison is intentionally case-insensitive on all platforms because these paths are
# Windows/SMB paths where casing is not meaningful, and because a missed match here means
# failing to protect a historical data root.
#
# The per-branch validate_*_output_root() functions keep their own names, signatures
# and error messages -- only their internals delegate here.
#
# Sourcing: callers must source this file BEFORE the *_utils.R file that uses it. Every
# entry-point script and test already sources R/analysis_labels.R first, so this line goes
# immediately alongside it.

PROTEOMICS_DATA_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"

project_normalize_path <- function(path, must_work = FALSE) {
  normalizePath(path.expand(trimws(as.character(path))), winslash = "/", mustWork = must_work)
}

project_path_is_within <- function(path, parent) {
  path <- tolower(project_normalize_path(path))
  parent <- sub("/+$", "", tolower(project_normalize_path(parent)))
  identical(path, parent) || startsWith(path, paste0(parent, "/"))
}

project_paths_overlap <- function(left, right) {
  project_path_is_within(left, right) || project_path_is_within(right, left)
}

#' Fail if `output_root` collides with any protected root.
#'
#' Shared engine for the branch-specific validators. `message_prefix` and `name_hit` let each
#' branch keep its established error wording (several tests match on that text).
project_validate_output_root <- function(output_root, protected_roots, message_prefix,
                                      name_hit = FALSE, hit_label = "protected") {
  output_root <- project_normalize_path(output_root)
  protected <- vapply(protected_roots, project_normalize_path, character(1), USE.NAMES = FALSE)
  hit <- vapply(protected, function(p) project_paths_overlap(output_root, p), logical(1))
  if (any(hit)) {
    if (isTRUE(name_hit)) {
      stop(message_prefix, output_root, " (", hit_label, ": ", protected[[which(hit)[[1]]]], ")", call. = FALSE)
    }
    stop(message_prefix, output_root, call. = FALSE)
  }
  output_root
}
