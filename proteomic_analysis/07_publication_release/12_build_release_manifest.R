#!/usr/bin/env Rscript

# Publication release, stage 12 -- the release manifest and checksums.
#
# Produces
#   provenance/release_manifest.tsv   one row per released file
#   provenance/SHA256SUMS.txt         checksums for the whole package
#
# The manifest is built by WALKING the release tree, not by trusting the build registry.
# Every file found must be registered, and every registered file must exist; either
# mismatch is a hard error. That is what makes the manifest a completeness claim rather
# than a list of the files someone remembered to add.
#
# Self-reference: release_manifest.tsv lists every file except itself and SHA256SUMS.txt.
# SHA256SUMS.txt covers every file except itself, so the manifest's own integrity is
# checkable from the checksum file.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_validation.R"))
release_require("digest")

OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/12_build_release_manifest.R"
MANIFEST_REL <- "provenance/release_manifest.tsv"
SUMS_REL <- "provenance/SHA256SUMS.txt"

release_banner("stage 12 -- release manifest")

if (!nrow(release_registry_read())) {
  stop("Build registry is empty; run the builder stages first.", call. = FALSE)
}
# The registry is itself a released file: it is the cross-stage record of which builder
# produced what, and it is what this manifest is checked against. Register it so the
# completeness check below covers it rather than tripping over it.
release_register(RELEASE_REGISTRY_RELPATH,
                 "internal build registry: which builder stage produced each release file",
                 "release builder stages", NA_character_, STAGE, "tsv")
registry <- release_registry_read()

all_files <- list.files(OUT_ROOT, recursive = TRUE, all.files = FALSE, full.names = FALSE)
all_files <- gsub("\\\\", "/", all_files)
self_files <- c(MANIFEST_REL, SUMS_REL)
manifest_files <- setdiff(all_files, self_files)

unregistered <- setdiff(manifest_files, registry$relative_path)
if (length(unregistered)) {
  stop("Release files present but not registered by any builder:\n",
       paste("  -", unregistered, collapse = "\n"), call. = FALSE)
}
missing <- setdiff(registry$relative_path, all_files)
if (length(missing)) {
  stop("Registered files that do not exist in the release:\n",
       paste("  -", missing, collapse = "\n"), call. = FALSE)
}
release_log("  ", length(manifest_files), " released files, all registered")

detect_format <- function(rp) {
  low <- tolower(rp)
  if (grepl("[.]tsv[.]gz$", low)) return("tsv.gz")
  if (grepl("[.]tsv$", low)) return("tsv")
  if (grepl("[.]csv$", low)) return("csv")
  if (grepl("[.]xlsx$", low)) return("xlsx")
  if (grepl("[.]md$", low)) return("markdown")
  if (grepl("[.]txt$", low)) return("text")
  sub("^.*[.]", "", low)
}

shape_of <- function(abs_path, fmt) {
  if (!fmt %in% c("tsv", "tsv.gz", "csv")) return(c(rows = NA_integer_, cols = NA_integer_))
  sep <- if (fmt == "csv") "," else "\t"
  con <- if (fmt == "tsv.gz") gzfile(abs_path, "rt") else file(abs_path, "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L, warn = FALSE)
  if (!length(header)) return(c(rows = 0L, cols = 0L))
  ncol <- length(strsplit(header, sep, fixed = TRUE)[[1]])
  nrow <- 0L
  repeat {
    chunk <- readLines(con, n = 50000L, warn = FALSE)
    if (!length(chunk)) break
    nrow <- nrow + length(chunk)
  }
  c(rows = as.integer(nrow), cols = as.integer(ncol))
}

git_commit <- release_git_commit(REPO_ROOT)
stamp <- release_timestamp_utc()

rows <- lapply(sort(manifest_files), function(rp) {
  abs_path <- file.path(OUT_ROOT, rp)
  reg <- registry[registry$relative_path == rp, , drop = FALSE]
  fmt <- detect_format(rp)
  shape <- shape_of(abs_path, fmt)
  data.frame(
    relative_path = rp,
    file_role = reg$file_role[[1]],
    format = fmt,
    rows_if_tabular = unname(shape[["rows"]]),
    columns_if_tabular = unname(shape[["cols"]]),
    bytes = as.integer(file.info(abs_path)$size),
    sha256 = release_sha256(abs_path),
    source_artifact = reg$source_artifact[[1]],
    canonical_source_sha256 = reg$canonical_source_sha256[[1]],
    generated_by = reg$generated_by[[1]],
    git_commit = git_commit,
    timestamp_utc = stamp,
    stringsAsFactors = FALSE, check.names = FALSE)
})
manifest <- do.call(rbind, rows)
rownames(manifest) <- NULL

if (anyDuplicated(manifest$relative_path)) {
  stop("Duplicate paths in the release manifest.", call. = FALSE)
}
if (any(is.na(manifest$sha256))) {
  stop("Could not hash every release file.", call. = FALSE)
}

by_dir <- table(dirname(manifest$relative_path))
for (d in names(by_dir)) {
  release_log("    ", format(d, width = 34), by_dir[[d]], " file(s)")
}
release_log("  total bytes: ", format(sum(manifest$bytes, na.rm = TRUE), big.mark = ","))

release_write_table(manifest, release_path(MANIFEST_REL))

# SHA256SUMS covers the manifest too, so the manifest's own integrity is verifiable.
# The directory is re-listed here, AFTER the manifest has been written: the listing taken
# at the top of this stage predates it, and using that one silently omitted the manifest
# from its own checksum file on a first build.
all_files <- gsub("\\\\", "/", list.files(OUT_ROOT, recursive = TRUE))
sums_targets <- sort(setdiff(all_files, SUMS_REL))
sums_lines <- vapply(sums_targets, function(rp) {
  paste0(release_sha256(file.path(OUT_ROOT, rp)), "  ", rp)
}, character(1), USE.NAMES = FALSE)
release_write_lines(sums_lines, release_path(SUMS_REL))

release_log("  wrote release_manifest.tsv (", nrow(manifest), " files)")
release_log("  wrote SHA256SUMS.txt (", length(sums_lines), " checksums)")
release_log("stage 12 complete")
