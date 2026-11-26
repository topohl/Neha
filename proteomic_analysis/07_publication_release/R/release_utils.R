# Shared helpers for the release layer: repository-relative paths, data and output roots,
# where the primary statistics run lives, output-root checks, reading and writing release
# files, and the build registry.
#
# The release reads analysis outputs and repackages them; it never recomputes a statistic.
# release_validate_output_root() refuses any output root that overlaps a protected input or
# analysis folder, using project_paths_overlap() from R/project_path_utils.R.
#
# After sourcing this file, call release_source_project_helpers(); it loads
# R/analysis_labels.R and R/project_path_utils.R, which these helpers use.

# 2.0.0 (2026-11-25): the primary statistics are the per-compartment limma run.
RELEASE_LAYER_VERSION <- "2.0.0"

# --------------------------------------------------------------------------------------
# repository / helper discovery
# --------------------------------------------------------------------------------------

#' Absolute path to `proteomic_analysis/`, discovered from the running script.
release_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    here <- normalizePath(dirname(sub("^--file=", "", file_arg)), winslash = "/", mustWork = FALSE)
    for (up in c(".", "..", "../..", "../../..")) {
      cand <- normalizePath(file.path(here, up), winslash = "/", mustWork = FALSE)
      if (file.exists(file.path(cand, "R", "analysis_labels.R"))) return(cand)
    }
  }
  cand <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (up in c(".", "..", "../..")) {
    up_cand <- normalizePath(file.path(cand, up), winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(up_cand, "R", "analysis_labels.R"))) return(up_cand)
  }
  stop("Cannot locate proteomic_analysis/: R/analysis_labels.R not found from ", cand, call. = FALSE)
}

#' Source the shared project helpers this layer depends on.
release_source_project_helpers <- function(repo_root = release_repo_root(),
                                          envir = parent.frame()) {
  source(file.path(repo_root, "R", "analysis_labels.R"), local = envir)
  source(file.path(repo_root, "R", "project_path_utils.R"), local = envir)
  invisible(repo_root)
}

#' Write paths inside proteomic_analysis/ relative to that folder.
#'
#' Release tables name repository files this way, so they do not carry the local folder of
#' the computer that built the release. Other paths (the shared drive, paths inside the
#' release, NA) are returned unchanged. A ";"-separated list is handled entry by entry.
release_repo_relative <- function(path, repo_root = release_repo_root()) {
  if (!length(path)) return(character(0))
  root <- tolower(paste0(sub("/+$", "", gsub("\\\\", "/", repo_root)), "/"))
  one <- function(p) {
    if (is.na(p) || !nzchar(p)) return(p)
    parts <- strsplit(p, ";", fixed = TRUE)[[1]]
    norm <- gsub("\\\\", "/", parts)
    inside <- startsWith(tolower(norm), root)
    parts[inside] <- substring(norm[inside], nchar(root) + 1L)
    paste(parts, collapse = ";")
  }
  vapply(as.character(path), one, character(1), USE.NAMES = FALSE)
}

#' Resolve a path as release tables write it: a relative path is read against
#' proteomic_analysis/, an absolute path is returned as it is.
release_repo_resolve <- function(path, repo_root = release_repo_root()) {
  path <- as.character(path)
  relative <- !is.na(path) & nzchar(path) & !grepl("^([A-Za-z]:|[/\\\\])", path)
  path[relative] <- file.path(repo_root, path[relative])
  path
}

#' TRUE for a relative path whose first folder exists in proteomic_analysis/.
#'
#' That is how a repository path written by release_repo_relative() is told apart from a
#' path inside the release (metadata/, provenance/ and so on, which the repository does
#' not have).
release_is_repo_path <- function(path, repo_root = release_repo_root()) {
  path <- gsub("\\\\", "/", as.character(path))
  ok <- !is.na(path) & nzchar(path) & grepl("/", path, fixed = TRUE) &
    !grepl("^([A-Za-z]:|/)", path)
  first <- sub("/.*$", "", path)
  ok & dir.exists(file.path(repo_root, first))
}

# --------------------------------------------------------------------------------------
# option / env override resolution (same idiom as every scientific stage)
# --------------------------------------------------------------------------------------

release_option_or_env <- function(option_name, env_name, default) {
  value <- getOption(option_name)
  if (!is.null(value) && nzchar(trimws(as.character(value)))) return(as.character(value))
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(trimws(value))) return(value)
  default
}

#' Canonical shared-storage data root (read-only for this layer).
release_data_root <- function() {
  release_option_or_env(
    "proteomics.project_data_root", "PROTEOMICS_PROJECT_DATA_ROOT",
    "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
  )
}

#' Enclosing project folder, one level above clusterProfiler.
#'
#' Holds pre-clusterProfiler acquisition-era artefacts the release must cite:
#' sample_annotation.xlsx (original `.d` acquisition paths), pg.matrix_raw.txt (the
#' 5,747-row search output) and protein_count.xlsx (Supplementary A source).
release_project_root <- function() {
  release_option_or_env(
    "proteomics.release_project_root", "PROTEOMICS_RELEASE_PROJECT_ROOT",
    "S:/Lab_Member/Tobi/Experiments/Collabs/Neha"
  )
}

#' Destination for the publication package.
#'
#' Defaults to the real shared-drive release root. Always point this at scratch while
#' testing: PROTEOMICS_RELEASE_OUTPUT_ROOT=<scratch>.
#'
#' The default moved on 2026-11-25 with the primary statistics. The pooled-run release
#' of 2026-09-07 stays frozen in 03_output/publication_release (a protected root).
release_output_root <- function() {
  release_option_or_env(
    "proteomics.release_output_root", "PROTEOMICS_RELEASE_OUTPUT_ROOT",
    file.path(release_data_root(), "03_output",
              "release_within_compartment_20261125")
  )
}

# --------------------------------------------------------------------------------------
# where the primary statistics run lives
# --------------------------------------------------------------------------------------
# Since 2026-11-24 the primary statistics are the per-compartment limma run. This is the
# one place the builders and tests learn where its files are; the model itself is
# described by RELEASE_PRIMARY_MODEL in release_validation.R. The pooled run's folders
# stay on disk for comparison only.
RELEASE_STATISTICS_RUN <- list(
  id = "within_compartment_20261124",
  date = "2026-11-24",
  stats_relative = c("02_data", "animal_level", "within_compartment_20261124"),
  gct = "stat_results_within_compartment.gct",
  run_manifest = "stat_results_within_compartment_run_manifest.csv",
  # GO-BP (GO_BP/) and GO-CC (GO_CC/) GSEA and ORA: the rerun with fgsea nPermSimple = 100000,
  # adopted 2026-11-25. Its ORA tables are byte-identical to the 20261124 run's.
  enrichment_relative = c("03_output", "enrichment", "enrichment_within_compartment_20261125"),
  enrichment_date = "2026-11-25",
  enrichment_runs = c("GO_BP", "GO_CC"),
  # KEGG GSEA was not rerun (its gene sets are downloaded at run time and not version-pinned),
  # so it stays in the GO_BP folder of the 20261124 run, at the fgsea default nPermSimple. That
  # run's GO-BP and GO-CC tables are superseded and are read only as the same-setting
  # reference of the sensitivity analyses.
  kegg_enrichment_relative = c("03_output", "enrichment", "enrichment_within_compartment_20261124"),
  kegg_enrichment_date = "2026-11-24",
  kegg_run = "GO_BP",
  figure_rebuild_relative = c("03_output", "figure_rebuild_within_compartment_20261125"),
  # The record of the rebuild: the script that ran its R scripts, in order.
  figure_rebuild_record = c("scripts", "run_rebuild.sh"),
  # Rebuild of 2026-11-24 on the 20261124 enrichment, superseded by the one above.
  superseded_figure_rebuild_relative = c("03_output", "figure_rebuild_within_compartment_20261124"),
  inferential_relative = c("03_output", "inferential_checks", "within_compartment_20261124"),
  # Superseded pooled run (one fit over all 48 columns), retained for comparison. Read only
  # for the "what changed from the pooled release" tables; its mapped and enrichment indexes
  # must both record this GCT hash before anything is counted from them.
  pooled_stats_relative = c("02_data", "animal_level"),
  pooled_stat_gct = "stat_results_for_ssGSEA_neha_proteome.gct",
  pooled_stat_gct_sha256 = "e1ae20f02e2747cfae3572933f2b23e6c770b92ef6810963a2806afb7adbe2b6",
  pooled_enrichment_relative = c("03_output", "enrichment",
                                 "enrichment_t_rank_validation_20260825"),
  pooled_release_relative = c("03_output", "publication_release"),
  # The GO GSEA rerun with a larger fgsea nPermSimple (enrichment_relative above) was adopted
  # on 2026-11-25. The sensitivity analyses (interaction, leave-one-animal-out, rank
  # robustness) were not rerun: they stay at the fgsea default and are compared with the
  # primary analysis at that same setting.
  enrichment_rerun_status = "adopted"
)

release_rel_path <- function(relative, data_root = release_data_root()) {
  do.call(file.path, c(list(data_root), as.list(relative)))
}

#' Folder holding the primary statistics GCT, its run manifest, split/ and mapped/.
release_stats_root <- function(data_root = release_data_root()) {
  release_rel_path(RELEASE_STATISTICS_RUN$stats_relative, data_root)
}

#' A file next to the primary statistics GCT (the GCT itself or its run manifest).
release_stat_run_path <- function(file, data_root = release_data_root()) {
  file.path(release_stats_root(data_root), file)
}

#' GO enrichment run folder: "GO_BP" or "GO_CC".
release_enrichment_root <- function(run = c("GO_BP", "GO_CC"), data_root = release_data_root()) {
  run <- match.arg(run)
  file.path(release_rel_path(RELEASE_STATISTICS_RUN$enrichment_relative, data_root), run)
}

#' The run folder holding the released KEGG GSEA (per_comparison/*/GSEA_KEGG*.csv).
release_kegg_root <- function(data_root = release_data_root()) {
  file.path(release_rel_path(RELEASE_STATISTICS_RUN$kegg_enrichment_relative, data_root),
            RELEASE_STATISTICS_RUN$kegg_run)
}

release_figure_rebuild_root <- function(data_root = release_data_root()) {
  release_rel_path(RELEASE_STATISTICS_RUN$figure_rebuild_relative, data_root)
}

release_figure_rebuild_record <- function(data_root = release_data_root()) {
  do.call(file.path, c(list(release_figure_rebuild_root(data_root)),
                       as.list(RELEASE_STATISTICS_RUN$figure_rebuild_record)))
}

#' Per-compartment interaction and robustness checks (05/06 in 04_differential_*).
release_inferential_root <- function(data_root = release_data_root()) {
  release_rel_path(RELEASE_STATISTICS_RUN$inferential_relative, data_root)
}

# --------------------------------------------------------------------------------------
# output-root checks
# --------------------------------------------------------------------------------------

#' Roots this layer must never write into.
#'
#' 03_output is not protected as a whole, because the release root lives
#' inside it (03_output/release_within_compartment_20261125). The individual
#' canonical analysis branches under 03_output are protected by name instead, and so is
#' the frozen pooled-run release (03_output/publication_release). Overlap is tested on
#' whole path components, so the new root's name does not collide with the frozen one.
release_protected_roots <- function(data_root = release_data_root()) {
  c(
    file.path(data_root, "01_input"),
    file.path(data_root, "02_data"),
    file.path(data_root, "99_historical"),
    file.path(data_root, "_migration_backup_20260826"),
    file.path(data_root, "03_output", "enrichment"),
    file.path(data_root, "03_output", "ewce"),
    file.path(data_root, "03_output", "pca"),
    file.path(data_root, "03_output", "synthesis"),
    file.path(data_root, "03_output", "inferential_checks"),
    # inside inferential_checks already; named so the guard reports it by name
    release_rel_path(RELEASE_STATISTICS_RUN$inferential_relative, data_root),
    file.path(data_root, "03_output", "qc"),
    file.path(data_root, "03_output", "reviewer_revision_animal_level_20260827"),
    file.path(data_root, "03_output", "manuscript_curation_20260827"),
    release_rel_path(RELEASE_STATISTICS_RUN$figure_rebuild_relative, data_root),
    release_rel_path(RELEASE_STATISTICS_RUN$superseded_figure_rebuild_relative, data_root),
    # pooled-run release of 2026-09-07: superseded, frozen, never rebuilt in place
    release_rel_path(RELEASE_STATISTICS_RUN$pooled_release_relative, data_root)
  )
}

#' Stop unless `output_root` is clear of every protected root.
#'
#' Overlap is checked in both directions and ignoring case (project_paths_overlap()), so an
#' output root that contains a protected root, such as the parent of 02_data, is refused too.
release_validate_output_root <- function(output_root, data_root = release_data_root()) {
  if (is.null(output_root) || !nzchar(trimws(as.character(output_root)))) {
    stop("Publication release output root is empty.", call. = FALSE)
  }
  resolved <- project_normalize_path(output_root)
  protected <- release_protected_roots(data_root)
  hit <- vapply(protected, function(p) project_paths_overlap(resolved, p), logical(1))
  if (any(hit)) {
    stop(
      "Refusing to write the publication release into a protected root: ", resolved,
      " (overlaps: ", project_normalize_path(protected[[which(hit)[[1]]]]), ")",
      call. = FALSE
    )
  }
  resolved
}

#' Is `output_root` on the shared drive, under the canonical data root?
release_is_shared_drive_root <- function(output_root, data_root = release_data_root()) {
  isTRUE(try(project_path_is_within(output_root, data_root), silent = TRUE))
}

#' Stop before writing to the shared drive unless PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE=true.
#' The default release root is on the shared drive, so a bare `Rscript run_release.R` would
#' otherwise write there.
#'
#' A shared-drive build is also refused while the repository has uncommitted or untracked
#' changes, because the recorded commit would not contain the code. Set
#' PROTEOMICS_RELEASE_ALLOW_DIRTY=true to build anyway; the commit is then recorded with a
#' "+dirty" suffix.
release_assert_shared_drive_opt_in <- function(output_root,
                                               data_root = release_data_root(),
                                               dirty = NULL) {
  if (!release_is_shared_drive_root(output_root, data_root)) return(invisible(TRUE))
  allow <- tolower(trimws(release_option_or_env(
    "proteomics.release_allow_shared_drive", "PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE", "")))
  if (!allow %in% c("true", "yes", "1")) {
    stop(
      "Refusing to write the publication release to the shared drive without an explicit ",
      "opt-in.\n  target: ", output_root,
      "\n  Build into scratch first:\n",
      "    PROTEOMICS_RELEASE_OUTPUT_ROOT=<scratch dir> Rscript 07_publication_release/run_release.R\n",
      "  To publish to the shared drive, set ",
      "PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE=true.",
      call. = FALSE)
  }
  if (is.null(dirty)) dirty <- release_git_tree_dirty()
  allow_dirty <- tolower(trimws(release_option_or_env(
    "proteomics.release_allow_dirty", "PROTEOMICS_RELEASE_ALLOW_DIRTY", "")))
  if (isTRUE(dirty) && !allow_dirty %in% c("true", "yes", "1")) {
    stop(
      "Refusing to write the publication release to the shared drive from a repository with ",
      "uncommitted or untracked changes: the recorded commit would not contain the code that ",
      "built it.\n  Commit first, or set PROTEOMICS_RELEASE_ALLOW_DIRTY=true to publish with ",
      "the commit recorded as <hash>+dirty.",
      call. = FALSE)
  }
  invisible(TRUE)
}

#' Create (if needed) and validate the release root, then return it normalised.
#'
#' Both checks run before the directory is created, so a refused target leaves no empty
#' folder.
release_prepare_output_root <- function(output_root = release_output_root(),
                                        data_root = release_data_root()) {
  release_validate_output_root(output_root, data_root)
  release_assert_shared_drive_opt_in(output_root, data_root)
  if (!dir.exists(output_root)) {
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_root)) {
    stop("Could not create publication release root: ", output_root, call. = FALSE)
  }
  release_validate_output_root(output_root, data_root)
}

#' Resolve a path inside the release, validating the root first and creating the parent.
release_path <- function(..., output_root = release_output_root(), create_dir = TRUE) {
  root <- release_validate_output_root(output_root)
  target <- file.path(root, ...)
  if (isTRUE(create_dir)) {
    parent <- dirname(target)
    if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  target
}

# --------------------------------------------------------------------------------------
# reading canonical artefacts
# --------------------------------------------------------------------------------------

release_require <- function(...) {
  pkgs <- c(...)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install required packages for the publication release: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' SHA256 of a file's bytes (project_file_sha256() in R/project_path_utils.R), or NA if the
#' file is absent. The manifest and figure stages use the NA to tell "not released" from
#' "hash mismatch".
release_sha256 <- function(path) {
  release_require("digest")
  if (!file.exists(path)) return(NA_character_)
  if (is.na(file.size(path))) return(NA_character_)
  project_file_sha256(path)
}

release_assert_exists <- function(path, what) {
  if (!file.exists(path)) {
    stop("Canonical input missing (", what, "): ", path, call. = FALSE)
  }
  invisible(path)
}

release_read_csv <- function(path, what = basename(path)) {
  release_assert_exists(path, what)
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                       fileEncoding = "UTF-8-BOM")
  names(d) <- sub("^\ufeff", "", names(d))
  d
}

release_read_tsv_plain <- function(path, what = basename(path)) {
  release_assert_exists(path, what)
  utils::read.delim(path, sep = "\t", quote = "\"", stringsAsFactors = FALSE,
                    check.names = FALSE, colClasses = "character")
}

#' Read a GCT v1.3 file into id / row-descriptor / column-descriptor / matrix parts.
#'
#' Kept separate from R/protigy_input_utils.R so the release layer does not load a GCT
#' writer.
release_read_gct <- function(path, what = basename(path)) {
  release_assert_exists(path, what)
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 3L) stop("GCT too short: ", path, call. = FALSE)
  version <- trimws(lines[[1]])
  dims <- suppressWarnings(as.integer(strsplit(lines[[2]], "\t", fixed = TRUE)[[1]]))
  if (length(dims) < 4L || anyNA(dims[1:4])) {
    stop("GCT dimension line is not four integers: ", path, call. = FALSE)
  }
  n_row <- dims[[1]]; n_col <- dims[[2]]; n_rdesc <- dims[[3]]; n_cdesc <- dims[[4]]

  header <- strsplit(lines[[3]], "\t", fixed = TRUE)[[1]]
  sample_ids <- header[seq.int(2L + n_rdesc, length.out = n_col)]

  cdesc <- NULL
  if (n_cdesc > 0L) {
    cd <- lapply(seq_len(n_cdesc), function(i) strsplit(lines[[3L + i]], "\t", fixed = TRUE)[[1]])
    cdesc_names <- vapply(cd, function(x) x[[1]], character(1))
    cdesc <- as.data.frame(
      stats::setNames(
        lapply(cd, function(x) as.character(x[seq.int(2L + n_rdesc, length.out = n_col)])),
        cdesc_names
      ),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    cdesc$sample_id <- sample_ids
  }

  body <- lines[seq.int(4L + n_cdesc, length.out = n_row)]
  split_body <- strsplit(body, "\t", fixed = TRUE)
  ids <- vapply(split_body, function(x) x[[1]], character(1))
  rdesc <- NULL
  if (n_rdesc > 0L) {
    rdesc_names <- header[seq.int(2L, length.out = n_rdesc)]
    rdesc <- as.data.frame(
      stats::setNames(
        lapply(seq_len(n_rdesc), function(j) {
          vapply(split_body, function(x) as.character(x[[1L + j]]), character(1))
        }),
        rdesc_names
      ),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    rdesc$id <- ids
  }
  mat <- t(vapply(split_body, function(x) {
    suppressWarnings(as.numeric(x[seq.int(2L + n_rdesc, length.out = n_col)]))
  }, numeric(n_col)))
  dimnames(mat) <- list(ids, sample_ids)

  list(version = version,
       dimensions = c(n_row = n_row, n_col = n_col,
                      n_row_descriptors = n_rdesc, n_col_descriptors = n_cdesc),
       ids = ids, sample_ids = sample_ids, rdesc = rdesc, cdesc = cdesc, mat = mat,
       path = path)
}

# --------------------------------------------------------------------------------------
# writing release artefacts
# --------------------------------------------------------------------------------------

#' Strip characters that would corrupt an unquoted TSV.
#'
#' GO term names and the protein descriptions carried in the retained search/quantification
#' output are free text; a stray tab or newline would shift every downstream
#' column. Replaced with a single space rather than dropped, so field content stays
#' readable.
release_sanitize_tsv <- function(df) {
  # Index columns by position so SDRF's standards-sanctioned repeated headers
  # (for example, multiple cleavage agents) are each sanitized independently.
  for (i in seq_along(df)) {
    if (is.character(df[[i]]) || is.factor(df[[i]])) {
      df[[i]] <- gsub("[\t\r\n]+", " ", as.character(df[[i]]))
    }
  }
  df
}

#' Format doubles with full round-trip precision so exports equal their source exactly.
#'
#' write.table() writes 15 significant digits, which does not round-trip every double;
#' 17 always do, so the exports pass the exact-equality checks against the canonical CSVs.
release_format_numeric <- function(df) {
  # Use positional indexing for the same repeated-header reason as above.
  for (i in seq_along(df)) {
    if (is.numeric(df[[i]]) && !is.integer(df[[i]])) {
      df[[i]] <- ifelse(is.na(df[[i]]), NA_character_,
                         vapply(df[[i]], function(x) format(x, digits = 17, trim = TRUE),
                                character(1)))
    } else if (is.logical(df[[i]])) {
      df[[i]] <- ifelse(is.na(df[[i]]), NA_character_,
                         ifelse(df[[i]], "TRUE", "FALSE"))
    }
  }
  df
}

#' Write a release table as TSV, optionally gzipped, and return its shape.
release_write_table <- function(df, path, gzip = grepl("[.]gz$", path)) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  df <- release_sanitize_tsv(df)
  df <- release_format_numeric(df)
  parent <- dirname(path)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  con <- if (isTRUE(gzip)) gzfile(path, open = "wb") else file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE, na = "", eol = "\n")
  invisible(list(path = path, rows = nrow(df), cols = ncol(df), columns = names(df)))
}

release_write_lines <- function(lines, path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeLines(enc2utf8(as.character(lines)), con, useBytes = TRUE)
  invisible(path)
}

# Stage 13 runs after the 12-stage build and writes these two reports. An in-place rebuild
# removes them first, or stage 12 would reject them as unregistered files; any other
# unexpected file still fails the manifest.
RELEASE_VALIDATION_REPORT_RELPATHS <- c(
  "provenance/validation_results.tsv",
  "provenance/VALIDATION_REPORT.md"
)

release_clear_previous_validation_reports <- function(output_root = release_output_root()) {
  root <- release_validate_output_root(output_root)
  targets <- file.path(root, RELEASE_VALIDATION_REPORT_RELPATHS)
  existing <- targets[file.exists(targets)]
  if (length(existing)) unlink(existing, recursive = FALSE, force = FALSE)
  remaining <- targets[file.exists(targets)]
  if (length(remaining)) {
    stop("Could not remove previous generated validation report(s): ",
         paste(remaining, collapse = "; "), call. = FALSE)
  }
  invisible(existing)
}

# --------------------------------------------------------------------------------------
# build registry: each stage records its files; 12_build_release_manifest.R checks them
# --------------------------------------------------------------------------------------

RELEASE_REGISTRY_RELPATH <- "provenance/_build_registry.tsv"

RELEASE_REGISTRY_COLUMNS <- c(
  "relative_path", "file_role", "source_artifact", "canonical_source_sha256",
  "generated_by", "format_hint"
)

release_registry_path <- function(output_root = release_output_root()) {
  release_path(RELEASE_REGISTRY_RELPATH, output_root = output_root)
}

release_registry_empty <- function() {
  as.data.frame(
    stats::setNames(
      replicate(length(RELEASE_REGISTRY_COLUMNS), character(0), simplify = FALSE),
      RELEASE_REGISTRY_COLUMNS),
    stringsAsFactors = FALSE, check.names = FALSE)
}

release_registry_reset <- function(output_root = release_output_root()) {
  p <- release_registry_path(output_root)
  release_write_table(release_registry_empty(), p, gzip = FALSE)
  invisible(p)
}

#' Record one produced release file so the manifest builder can describe it.
#'
#' `source_artifact` and `canonical_source_sha256` accept vectors; they are joined with
#' ";" so a table assembled from twelve canonical CSVs can name all twelve.
release_register <- function(relative_path, file_role, source_artifact,
                             canonical_source_sha256 = NA_character_,
                             generated_by, format_hint = NA_character_,
                             output_root = release_output_root()) {
  p <- release_registry_path(output_root)
  # A source inside this release is recorded by its release-relative path, so the registry
  # and the manifest read the same whichever root the release was built in.
  source_artifact <- gsub("\\\\", "/", as.character(source_artifact))
  root_prefix <- paste0(sub("/+$", "", gsub("\\\\", "/", project_normalize_path(output_root))),
                        "/")
  inside <- !is.na(source_artifact) &
    startsWith(tolower(source_artifact), tolower(root_prefix))
  source_artifact[inside] <- substring(source_artifact[inside], nchar(root_prefix) + 1L)
  # A source inside the repository is recorded relative to proteomic_analysis/.
  repo_root <- tryCatch(release_repo_root(), error = function(e) NULL)
  if (!is.null(repo_root)) source_artifact <- release_repo_relative(source_artifact, repo_root)
  row <- data.frame(
    relative_path = gsub("\\\\", "/", as.character(relative_path)[[1]]),
    file_role = as.character(file_role)[[1]],
    source_artifact = paste(unique(source_artifact), collapse = ";"),
    canonical_source_sha256 = paste(
      unique(stats::na.omit(as.character(canonical_source_sha256))), collapse = ";"),
    generated_by = as.character(generated_by)[[1]],
    format_hint = as.character(format_hint)[[1]],
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!nzchar(row$canonical_source_sha256)) row$canonical_source_sha256 <- NA_character_
  existing <- if (file.exists(p)) release_registry_read(output_root) else release_registry_empty()
  existing <- existing[existing$relative_path != row$relative_path, , drop = FALSE]
  release_write_table(rbind(existing, row), p, gzip = FALSE)
  invisible(row$relative_path)
}

release_registry_read <- function(output_root = release_output_root()) {
  p <- release_registry_path(output_root)
  if (!file.exists(p)) return(release_registry_empty())
  d <- utils::read.delim(p, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
                         colClasses = "character")
  if (!nrow(d)) return(release_registry_empty())
  d
}

# --------------------------------------------------------------------------------------
# small shared conveniences
# --------------------------------------------------------------------------------------

#' TRUE when proteomic_analysis/ has uncommitted or untracked changes, NA if git is absent.
release_git_tree_dirty <- function(repo_root = release_repo_root()) {
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", shQuote(repo_root), "status", "--porcelain",
                                      "--untracked-files=all", "--", "."),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) NULL
  )
  status <- attr(out, "status")
  if (is.null(out) || (!is.null(status) && status != 0L)) return(NA)
  any(nzchar(trimws(out)))
}

#' HEAD of the repository, suffixed "+dirty" when the working tree differs from it.
release_git_commit <- function(repo_root = release_repo_root()) {
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", shQuote(repo_root), "rev-parse", "HEAD"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  if (!length(out) || is.na(out[[1]]) || !nzchar(out[[1]])) return("UNKNOWN")
  commit <- trimws(out[[1]])
  if (isTRUE(release_git_tree_dirty(repo_root))) commit <- paste0(commit, "+dirty")
  commit
}

release_timestamp_utc <- function() {
  format(as.POSIXct(Sys.time()), tz = "UTC", "%Y-%m-%d %H:%M:%S UTC")
}

release_log <- function(...) cat(..., "\n", sep = "")

release_banner <- function(stage) cat("\n=== ", stage, " ===\n", sep = "")

#' Canonical pairing/treatment decomposition of a condition label.
release_split_condition <- function(condition) {
  condition <- as.character(condition)
  list(
    pairing_status = ifelse(grepl("^paired_", condition), "paired",
                            ifelse(grepl("^unpaired_", condition), "unpaired", NA_character_)),
    treatment = ifelse(grepl("_cno$", condition), "cno",
                       ifelse(grepl("_veh$", condition), "veh", NA_character_))
  )
}

#' Historical filename alias for a sample class (`neuropil` was written `bg`).
release_sample_class_alias <- function(sample_class) {
  ifelse(as.character(sample_class) == "neuropil", "bg", as.character(sample_class))
}

#' Normalise the 2024 acquisition-annotation group labels to canonical sample classes.
#'
#' sample_annotation.xlsx uses `Background`, `cFosN`, `mCherryN`, `Neuron` (and `N`).
#' normalize_sample_class() in R/analysis_labels.R returns NA for `cFosN` and `mCherryN`,
#' and a comparison against NA would hide a mismatch, so the extra aliases live here
#' instead of in the shared helper the analysis stages use.
RELEASE_ANNOTATION_LABEL_ALIASES <- c(
  "background" = "neuropil",
  "bg" = "neuropil",
  "cfosn" = "cfos",
  "cfos" = "cfos",
  "mcherryn" = "mcherry",
  "mcherry" = "mcherry",
  "neuron" = "neuron",
  "n" = "neuron"
)

release_normalize_annotation_label <- function(x) {
  key <- tolower(trimws(as.character(x)))
  out <- unname(RELEASE_ANNOTATION_LABEL_ALIASES[key])
  out[is.na(key)] <- NA_character_
  out
}
