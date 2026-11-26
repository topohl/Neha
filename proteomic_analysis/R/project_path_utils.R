# Path and file-hash helpers shared by the stages and the R/*_utils.R files.
#
# Paths are compared after trimws(), path.expand() and normalizePath(), case-insensitively on
# every OS (they are Windows/SMB paths), and overlap is checked in both directions. The
# per-stage validate_*_output_root() functions keep their own names and error messages and
# use these helpers.
#
# Source this file before any R/*_utils.R file that uses it.

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
#' `message_prefix` and `name_hit` keep each caller's error wording, which tests match.
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

# --------------------------------------------------------------------------------------
# File hashing.
#
# digest::digest(file = ...) calls file.access() first, and on the S: share that call can
# return -1 for files that read fine, so digest() stops with "The specified file is not
# readable". project_file_sha256() reads the bytes itself and returns the same digest.
#
# Still on digest(file = ...): protigy_file_sha256() in R/protigy_stat_gct_utils.R and the
# provenance hashes in 01_preprocessing/02a_prepare_animal_level_protigy_input.r,
# 03_qc_exploration/pca/06a_pca_core.r and 05_celltype_enrichment_EWCE/01_EWCE.r. Switch
# them when those stages are next re-run.
# --------------------------------------------------------------------------------------

#' SHA-256 of a file's bytes, read in binary mode.
#'
#' Stops if the file cannot be opened or read (never returns NA). Returns a lowercase
#' 64-character hex digest.
project_file_sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required to hash file contents.", call. = FALSE)
  }
  path <- as.character(path)
  if (length(path) != 1L || is.na(path) || !nzchar(trimws(path))) {
    stop("project_file_sha256() needs a single non-empty path.", call. = FALSE)
  }
  path <- trimws(path)
  if (dir.exists(path)) {
    stop("Cannot hash a directory: ", path, call. = FALSE)
  }
  con <- tryCatch(suppressWarnings(file(path, open = "rb")), error = function(e) NULL)
  if (is.null(con)) {
    stop("Cannot open file for hashing: ", path, call. = FALSE)
  }
  on.exit(close(con), add = TRUE)
  # Read to EOF in chunks; file.size() can be out of date on SMB.
  chunks <- list()
  repeat {
    chunk <- readBin(con, what = "raw", n = 1048576L)
    if (!length(chunk)) break
    chunks[[length(chunks) + 1L]] <- chunk
  }
  bytes <- if (length(chunks)) unlist(chunks, use.names = FALSE) else raw(0)
  tolower(digest::digest(bytes, algo = "sha256", serialize = FALSE))
}
