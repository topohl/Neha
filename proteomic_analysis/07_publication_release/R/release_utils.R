# Shared primitives for the publication / deposition release layer.
#
# This layer is READ-ONLY with respect to validated analysis data. It reads canonical
# artefacts and repackages them; it never recomputes a statistic and never writes into a
# canonical root. The write-side guard is release_validate_output_root(), which fails
# closed before any builder opens a connection.
#
# Sourcing order: R/analysis_labels.R and R/project_path_utils.R must be sourced first
# (release_source_project_helpers() does that for you). The containment engine here is
# project_paths_overlap() from R/project_path_utils.R -- deliberately reused rather than
# reimplemented, because a divergent normalisation is exactly how a protected root gets
# missed (see the comment block at the top of that file).

RELEASE_LAYER_VERSION <- "1.0.0"

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
release_output_root <- function() {
  release_option_or_env(
    "proteomics.release_output_root", "PROTEOMICS_RELEASE_OUTPUT_ROOT",
    file.path(release_data_root(), "03_output", "publication_release")
  )
}

# --------------------------------------------------------------------------------------
# write-side governance
# --------------------------------------------------------------------------------------

#' Roots this layer must never write into.
#'
#' 03_output is NOT protected wholesale, because the intended real release root lives
#' inside it (03_output/publication_release). The individual canonical analysis branches
#' under 03_output are protected by name instead.
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
    file.path(data_root, "03_output", "qc"),
    file.path(data_root, "03_output", "reviewer_revision_animal_level_20260827"),
    file.path(data_root, "03_output", "manuscript_curation_20260827")
  )
}

#' Fail closed unless `output_root` is clear of every protected root.
#'
#' Containment is checked bidirectionally and case-insensitively via
#' project_paths_overlap(), so an output root that *contains* a protected root is
#' rejected too. Pointing the release at the parent of 02_data must not pass merely
#' because it is not literally 02_data.
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

#' Fail closed unless writing to the shared drive has been explicitly opted into.
#'
#' The default release root is on the shared drive, so a bare `Rscript run_release.R`
#' would populate it. That is the same failure mode as the 2026-08-28 incident recorded in
#' the repository history: a run that was meant for scratch reached the validated tree
#' because nothing required the destination to be stated. Building into scratch stays
#' frictionless; publishing to the shared drive requires saying so.
release_assert_shared_drive_opt_in <- function(output_root,
                                               data_root = release_data_root()) {
  if (!release_is_shared_drive_root(output_root, data_root)) return(invisible(TRUE))
  allow <- tolower(trimws(release_option_or_env(
    "proteomics.release_allow_shared_drive", "PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE", "")))
  if (allow %in% c("true", "yes", "1")) return(invisible(TRUE))
  stop(
    "Refusing to write the publication release to the shared drive without an explicit ",
    "opt-in.\n  target: ", output_root,
    "\n  Build into scratch first:\n",
    "    PROTEOMICS_RELEASE_OUTPUT_ROOT=<scratch dir> Rscript 07_publication_release/run_release.R\n",
    "  To publish to the shared drive deliberately, set ",
    "PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE=true.",
    call. = FALSE)
}

#' Create (if needed) and validate the release root, then return it normalised.
#'
#' Order matters: both guards run BEFORE the directory is created, so a rejected target
#' does not even leave an empty folder behind on the shared drive.
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

#' SHA256 of a file's bytes, or NA if the file is absent.
#'
#' Hashes bytes read here rather than calling digest(file = ...). digest() probes the file
#' with file.access() first, and on this project's SMB share that probe returns -1 --
#' "No mapping between account names and security IDs was done" -- for files that read
#' perfectly well. When the share is in that state every hash in the layer raises "The
#' specified file is not readable" and a whole release build aborts on intact data.
#'
#' Reading the bytes is both more reliable and strictly stronger evidence than an ACL
#' probe: it cannot succeed on a file it could not read. The digest is identical -- verified
#' against both locked GCT hashes, which this function still reproduces exactly.
release_sha256 <- function(path) {
  release_require("digest")
  if (!file.exists(path)) return(NA_character_)
  size <- file.size(path)
  if (is.na(size)) return(NA_character_)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  bytes <- readBin(con, what = "raw", n = size)
  if (length(bytes) != size) {
    stop("Short read while hashing ", path, ": got ", length(bytes), " of ", size,
         " bytes.", call. = FALSE)
  }
  digest::digest(bytes, algo = "sha256", serialize = FALSE)
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
#' Implemented here rather than reused from R/protigy_input_utils.R because that helper
#' belongs to a *producing* stage; the release layer must not pull a GCT writer into
#' scope at all.
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
#' GO term names and DIA-NN protein descriptions are free text; a stray tab or newline
#' would silently shift every downstream column. Replaced with a single space rather
#' than dropped, so field content stays readable.
release_sanitize_tsv <- function(df) {
  for (nm in names(df)) {
    if (is.character(df[[nm]]) || is.factor(df[[nm]])) {
      df[[nm]] <- gsub("[\t\r\n]+", " ", as.character(df[[nm]]))
    }
  }
  df
}

#' Format doubles with full round-trip precision so exports equal their source exactly.
#'
#' write.table() would emit 7 significant digits, which breaks the exact-equality
#' contract the release tests assert against the canonical CSVs. 17 significant digits
#' round-trips an IEEE double.
release_format_numeric <- function(df) {
  for (nm in names(df)) {
    if (is.numeric(df[[nm]]) && !is.integer(df[[nm]])) {
      df[[nm]] <- ifelse(is.na(df[[nm]]), NA_character_,
                         vapply(df[[nm]], function(x) format(x, digits = 17, trim = TRUE),
                                character(1)))
    } else if (is.logical(df[[nm]])) {
      df[[nm]] <- ifelse(is.na(df[[nm]]), NA_character_,
                         ifelse(df[[nm]], "TRUE", "FALSE"))
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

# --------------------------------------------------------------------------------------
# build registry -- the cross-process channel 09_build_release_manifest.R aggregates
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
  row <- data.frame(
    relative_path = gsub("\\\\", "/", as.character(relative_path)[[1]]),
    file_role = as.character(file_role)[[1]],
    source_artifact = paste(unique(as.character(source_artifact)), collapse = ";"),
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

release_git_commit <- function(repo_root = release_repo_root()) {
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", shQuote(repo_root), "rev-parse", "HEAD"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  if (!length(out) || is.na(out[[1]]) || !nzchar(out[[1]])) return("UNKNOWN")
  trimws(out[[1]])
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
#' sample_annotation.xlsx uses `Background`, `cFosN`, `mCherryN`, `Neuron` (and `N` in the
#' hemisphere label). Only two of those four are in sample_class_aliases in
#' R/analysis_labels.R, so normalize_sample_class() returns NA for `cFosN` and `mCherryN`.
#' That matters here: a corroboration check comparing against NA silently passes, which
#' would hide exactly the kind of mismatch the check exists to find. The shared helper is
#' left alone -- it is used by the scientific stages and changing its vocabulary could
#' alter their behaviour -- and the extra aliases live here instead.
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
