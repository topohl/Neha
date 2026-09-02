#!/usr/bin/env Rscript

# Contracts for project_file_sha256(), the shared file-content hashing primitive in
# R/project_path_utils.R.
#
# Added 2026-09-02. The primitive exists because digest::digest(file = ...) decides whether
# a file is readable by probing it with file.access(), and on this project's SMB share that
# probe returns -1 -- "No mapping between account names and security IDs was done" -- for
# files that read back byte-perfect and hash to their locked values. digest() then aborts
# with "The specified file is not readable" on wholly intact data, which is exactly how
# test_data_integrity and test_rank_abundance_animal_level came to fail on unchanged
# canonical artefacts.
#
# These contracts pin the replacement in both directions: it must hash real bytes and agree
# with published SHA-256 digests, and it must still be strictly content sensitive -- a
# single changed byte has to change the result, and an unreadable file has to raise rather
# than yield a value. They run entirely on local temporary files, so they carry no
# shared-drive dependency; the SMB behaviour itself is exercised by the canonical sections
# of the two tests named above.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("tests", "test_project_file_sha256.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "project_path_utils.R"))

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required for the hashing primitive contracts.", call. = FALSE)
}

failures <- 0L
expect <- function(condition, message) {
  if (!isTRUE(condition)) {
    failures <<- failures + 1L
    cat("  [FAIL]", message, "\n")
  } else {
    cat("  [ok]  ", message, "\n")
  }
}

work <- file.path(tempdir(), paste0("project_file_sha256_", as.integer(Sys.time())))
dir.create(work, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

write_bytes <- function(path, bytes) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(bytes, con)
  path
}

# --------------------------------------------------------------------------------------
cat("=== agreement with published SHA-256 digests ===\n")
# Known-answer vectors from FIPS 180-4 / RFC 6234, so correctness is checked against
# published digests rather than against another call into the same library.
known <- list(
  list(id = "empty", label = "empty input", bytes = raw(0),
       sha = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),
  list(id = "abc", label = "the three bytes abc", bytes = charToRaw("abc"),
       sha = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"),
  list(id = "helloworld", label = "the string hello world", bytes = charToRaw("hello world"),
       sha = "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9")
)
for (k in known) {
  p <- write_bytes(file.path(work, paste0("known_", k$id, ".bin")), k$bytes)
  expect(identical(project_file_sha256(p), k$sha),
         sprintf("%s hashes to its published SHA-256", k$label))
}

# Called out separately because the empty file is the one input where a byte-reading
# implementation could plausibly return NA, or error, instead of a digest.
empty <- write_bytes(file.path(work, "empty_again.bin"), raw(0))
expect(identical(project_file_sha256(empty),
                 "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),
       "an empty file hashes rather than erroring or returning NA")

cat("\n=== return shape ===\n")
payload <- write_bytes(file.path(work, "payload.bin"),
                       as.raw(c(0L, 1L, 2L, 255L, 128L, 10L, 13L)))
h <- project_file_sha256(payload)
expect(is.character(h) && length(h) == 1L, "returns a length-1 character vector")
expect(!is.na(h) && nchar(h) == 64L, "returns 64 characters")
expect(grepl("^[0-9a-f]{64}$", h), "returns lowercase hexadecimal only")

cat("\n=== content sensitivity: the primitive must be no weaker than digest(file=) ===\n")
base_bytes <- as.raw(rep(c(1L, 2L, 3L, 4L), length.out = 4096L))
a <- write_bytes(file.path(work, "a.bin"), base_bytes)
flipped <- base_bytes
flipped[[2048L]] <- as.raw(5L)
b <- write_bytes(file.path(work, "b.bin"), flipped)
expect(!identical(project_file_sha256(a), project_file_sha256(b)),
       "a single changed byte changes the digest")

truncated <- write_bytes(file.path(work, "truncated.bin"), base_bytes[seq_len(4095L)])
expect(!identical(project_file_sha256(a), project_file_sha256(truncated)),
       "a truncated file changes the digest")

appended <- write_bytes(file.path(work, "appended.bin"), c(base_bytes, as.raw(0L)))
expect(!identical(project_file_sha256(a), project_file_sha256(appended)),
       "a trailing NUL byte changes the digest")

cat("\n=== path independence: identical bytes hash identically ===\n")
nested <- file.path(work, "nested", "deeper")
dir.create(nested, recursive = TRUE, showWarnings = FALSE)
copy <- file.path(nested, "renamed_copy.dat")
expect(file.copy(a, copy, overwrite = TRUE), "fixture copied to a second path")
expect(identical(project_file_sha256(a), project_file_sha256(copy)),
       "the same bytes at a different path and filename hash the same")

# Bytes are hashed as bytes, so a CRLF/LF difference has to stay visible. This is what
# keeps the primitive honest on a Windows/SMB checkout, where a text-mode read would
# silently make the two agree.
lf <- write_bytes(file.path(work, "lf.txt"), charToRaw("line1\nline2\n"))
crlf <- write_bytes(file.path(work, "crlf.txt"), charToRaw("line1\r\nline2\r\n"))
expect(!identical(project_file_sha256(lf), project_file_sha256(crlf)),
       "LF and CRLF variants of the same text hash differently (binary read, not text)")

cat("\n=== larger-than-one-chunk input ===\n")
# The implementation reads in 1 MiB chunks, so this input spans several: a chunk-boundary
# or reassembly bug cannot hide behind a fixture that fits in a single read.
set.seed(20260902L)
big_bytes <- as.raw(sample.int(256L, 3L * 1024L * 1024L + 12345L, replace = TRUE) - 1L)
big <- write_bytes(file.path(work, "big.bin"), big_bytes)
expect(identical(project_file_sha256(big),
                 tolower(digest::digest(big_bytes, algo = "sha256", serialize = FALSE))),
       "a multi-chunk file hashes to the digest of its full byte vector")
idx <- 2L * 1024L * 1024L
big_flipped <- big_bytes
big_flipped[[idx]] <- as.raw((as.integer(big_bytes[[idx]]) + 1L) %% 256L)
big2 <- write_bytes(file.path(work, "big2.bin"), big_flipped)
expect(!identical(project_file_sha256(big), project_file_sha256(big2)),
       "a byte flipped past the first chunk boundary changes the digest")

cat("\n=== unreadable input errors instead of returning a value ===\n")
errs <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
expect(errs(project_file_sha256(file.path(work, "does_not_exist.bin"))),
       "a nonexistent file raises an error")
expect(errs(project_file_sha256(work)), "a directory raises an error")
expect(errs(project_file_sha256(NA_character_)), "NA raises an error")
expect(errs(project_file_sha256("")), "an empty path raises an error")
expect(errs(project_file_sha256(c(a, b))), "a vector of paths raises an error")

cat("\n=== the implementation takes no forbidden dependency ===\n")
# A source-level contract, because a regression here would be silent: the function would
# keep returning correct digests everywhere except on the share this project runs against.
#
# Matching is done against parsed-and-deparsed code, not raw file text. Every file touched
# by this fix *documents* the digest(file = ...) hazard in a comment, and a grep over raw
# text would match that prose and report a call that is not there. parse() drops comments,
# so what is asserted is the code.
code_of <- function(rel) {
  exprs <- parse(file.path(repo_root, rel), keep.source = FALSE)
  paste(unlist(lapply(exprs, function(e) deparse(e))), collapse = "\n")
}
impl <- code_of(file.path("R", "project_path_utils.R"))
body_txt <- paste(deparse(body(project_file_sha256)), collapse = "\n")
expect(!grepl("file.access", body_txt, fixed = TRUE),
       "project_file_sha256() does not call file.access()")
expect(!grepl("file *= *TRUE", body_txt) && !grepl("digest *\\( *file *=", body_txt),
       "project_file_sha256() does not call digest(file = ...)")
expect(grepl('open *= *"rb"', body_txt),
       "project_file_sha256() opens the file in binary mode")
expect(grepl("readBin", body_txt, fixed = TRUE),
       "project_file_sha256() reads the bytes itself")
expect(!grepl("file.access", impl, fixed = TRUE),
       "R/project_path_utils.R contains no file.access() call at all")

cat("\n=== the hashers that broke on SMB now route through the primitive ===\n")
callers <- list(
  c("R/animal_level_enrichment_utils.R", "enrichment_sha256"),
  c("R/mapthatprot_animal_level_utils.R", "mapthatprot_sha256"),
  c("07_publication_release/R/release_utils.R", "release_sha256"),
  c("tests/test_data_integrity.R", "the data-integrity contracts"),
  c("tests/test_rank_abundance_animal_level.R", "the rank-abundance contracts")
)
for (cl in callers) {
  code <- code_of(cl[[1]])
  expect(grepl("project_file_sha256", code, fixed = TRUE),
         sprintf("%s hashes via project_file_sha256()", cl[[2]]))
  expect(!grepl("digest *\\( *file *=", code) && !grepl("file *= *TRUE", code),
         sprintf("%s no longer calls digest(file = ...)", cl[[2]]))
}

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("project_file_sha256 contracts failed: %d", failures), call. = FALSE)
}
cat("All project_file_sha256 contracts hold.\n")
