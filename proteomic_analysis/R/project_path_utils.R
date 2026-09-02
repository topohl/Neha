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

# --------------------------------------------------------------------------------------
# Shared file-content hashing primitive.
#
# Added 2026-09-02. Four near-identical SHA-256 helpers had grown up independently --
# enrichment_sha256(), mapthatprot_sha256(), protigy_file_sha256() and the publication
# layer's release_sha256() -- and the first three delegated straight to
# digest::digest(file = ...). That is not a safe way to hash a file on this project's
# storage: digest() probes the path with file.access() before opening it, and on the SMB
# share that probe returns -1 -- "No mapping between account names and security IDs was
# done" -- for files that read back perfectly and hash to their locked values. digest()
# then aborts with "The specified file is not readable" on wholly intact data, which is
# how test_data_integrity and test_rank_abundance_animal_level came to fail on unchanged
# canonical artefacts.
#
# So the primitive lives here, alongside the other shared path helpers, and hashes bytes
# it has actually read. The only authority on whether a file is readable is whether the
# connection opens and the bytes come back -- never an ACL probe. Content sensitivity is
# unchanged and deliberately so: these hashes still change if a single byte changes, and
# they still reproduce the locked values recorded across the pipeline.
#
# Callers converted so far: enrichment_sha256(), mapthatprot_sha256(), release_sha256(),
# the rank-abundance stage's input provenance, and the two data-integrity/rank-abundance
# test suites. Still on digest(file = ...) and therefore still SMB-fragile:
# protigy_file_sha256() in R/protigy_stat_gct_utils.R, and the provenance hashes in
# 01_preprocessing/02a_prepare_animal_level_protigy_input.r, 03_qc_exploration/pca/
# 06a_pca_core.r and 05_celltype_enrichment_EWCE/01_EWCE.r. Those were left alone
# deliberately: converting them means editing producing stages for animal aggregation,
# PCA and EWCE, whose outputs are frozen canonical values, and no test currently exercises
# them against the share. Convert them the next time one of those stages is legitimately
# re-run, not as a drive-by.
# --------------------------------------------------------------------------------------

#' SHA-256 of a file's bytes, read in binary mode.
#'
#' Errors -- rather than returning NA -- when the file cannot genuinely be opened and read,
#' so a caller can never mistake an unreadable file for a hashed one. Returns a lowercase
#' 64-character digest. Never calls file.access(); never calls digest(file = ...).
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
  # Read to EOF in chunks rather than trusting a single file.size() call: on SMB the
  # reported size can be stale, and a short read must never quietly hash a partial file.
  chunks <- list()
  repeat {
    chunk <- readBin(con, what = "raw", n = 1048576L)
    if (!length(chunk)) break
    chunks[[length(chunks) + 1L]] <- chunk
  }
  bytes <- if (length(chunks)) unlist(chunks, use.names = FALSE) else raw(0)
  tolower(digest::digest(bytes, algo = "sha256", serialize = FALSE))
}
