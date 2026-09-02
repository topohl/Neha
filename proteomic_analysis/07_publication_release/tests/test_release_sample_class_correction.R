#!/usr/bin/env Rscript

# Contract: the six C46/C47 sample-class corrections are published exactly, both identities
# are retained for all 96 measurements, and nothing about the correction is left as an
# unresolved blocker.
#
# The point of this file is that the six rows cannot drift. Every field of the correction is
# pinned -- the identities, the hemisphere, the animals, the plate positions, both class
# vectors and the direction of the cycle -- and the preserved correction record is verified
# by hash rather than by modification time.
#
# Skips cleanly when no built release is reachable, so it is safe in CI. Point it at one
# with PROTEOMICS_RELEASE_OUTPUT_ROOT.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("07_publication_release", "tests", "test_release_sample_class_correction.R")
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

OUT_ROOT <- release_output_root()
DATA_ROOT <- release_data_root()
sm_path <- file.path(OUT_ROOT, "metadata", "sample_metadata.tsv")
corr_path <- file.path(OUT_ROOT, "metadata", "sample_class_corrections.tsv")
if (!file.exists(sm_path) || !file.exists(corr_path)) {
  cat("No built publication release at", OUT_ROOT, "\n")
  cat("  set PROTEOMICS_RELEASE_OUTPUT_ROOT and run 07_publication_release/run_release.R.\n")
  cat("Skipping publication release sample-class correction contracts.\n")
  quit(save = "no", status = 0L)
}

read_tsv <- function(p) read.delim(p, sep = "\t", quote = "", stringsAsFactors = FALSE,
                                   check.names = FALSE)
sm <- read_tsv(sm_path)
corr <- read_tsv(corr_path)
inv <- RELEASE_DESIGN_INVARIANTS
SIX <- c("N53", "N60", "N67", "N81", "N88", "N95")

# --------------------------------------------------------------------------------------
cat("=== A. exactly six corrected rows ===\n")
# --------------------------------------------------------------------------------------
expect(nrow(corr) == 6L, sprintf("the correction table has exactly 6 rows (got %d)",
                                 nrow(corr)))
expect(nrow(corr) == RELEASE_N_SAMPLE_CLASS_CORRECTED,
       "the correction table matches the release contract count")
expect(sum(as.logical(sm$sample_class_corrected)) == 6L,
       sprintf("exactly 6 of 96 measurements are flagged corrected (got %d)",
               sum(as.logical(sm$sample_class_corrected))))

# --------------------------------------------------------------------------------------
cat("\n=== B. the corrected set is exactly N53 N60 N67 N81 N88 N95 ===\n")
# --------------------------------------------------------------------------------------
expect(setequal(corr$plate_sample_number, SIX),
       paste("corrected set is", paste(sort(corr$plate_sample_number), collapse = " ")))
expect(setequal(sm$plate_sample_number[as.logical(sm$sample_class_corrected)], SIX),
       "the sample metadata flags the same six measurements")
expect(setequal(corr$sample_id, sm$sample_id[as.logical(sm$sample_class_corrected)]),
       "correction table and sample metadata agree on the six acquisition ids")
expect(setequal(corr$AnimalID, c("C46", "C47")),
       "the corrected animals are exactly C46 and C47")
expect(all(table(corr$AnimalID) == 3L), "three corrected measurements per animal")

# --------------------------------------------------------------------------------------
cat("\n=== C. all six are Left hemisphere ===\n")
# --------------------------------------------------------------------------------------
expect(all(corr$hemisphere == "Left"), "every corrected measurement is Left hemisphere")
expect(all(sm$hemisphere[as.logical(sm$sample_class_corrected)] == "Left"),
       "the sample metadata agrees that all six are Left")
expect(!any(as.logical(sm$sample_class_corrected) & sm$hemisphere == "Right"),
       "no Right-hemisphere measurement is flagged corrected")

# --------------------------------------------------------------------------------------
cat("\n=== D. the mappings are exactly the cyclic reassignment ===\n")
# --------------------------------------------------------------------------------------
EXPECTED <- c(N53 = "neuropil->mcherry", N60 = "mcherry->neuron",
              N67 = "neuron->neuropil", N81 = "neuropil->mcherry",
              N88 = "mcherry->neuron", N95 = "neuron->neuropil")
observed <- stats::setNames(
  paste0(corr$original_sample_class, "->", corr$analysis_sample_class),
  corr$plate_sample_number)
for (n in SIX) {
  expect(identical(unname(observed[n]), unname(EXPECTED[n])),
         sprintf("%s: %s (expected %s)", n,
                 ifelse(is.na(observed[n]), "MISSING", observed[n]), EXPECTED[n]))
}
expect(identical(unname(RELEASE_SAMPLE_CLASS_CORRECTION$cycle[corr$original_sample_class]),
                 as.character(corr$analysis_sample_class)),
       "every correction follows neuropil -> mcherry -> neuron -> neuropil")
expect(!any(corr$original_sample_class == corr$analysis_sample_class),
       "no corrected row leaves the class unchanged")
expect(!("cfos" %in% c(corr$original_sample_class, corr$analysis_sample_class)),
       "the cfos-enriched samples are not involved in the correction")

correction_checks <- release_check_sample_class_corrections(corr)
for (i in seq_len(nrow(correction_checks))) {
  expect(correction_checks$pass[i], correction_checks$check[i])
}

# --------------------------------------------------------------------------------------
cat("\n=== E. the other 90 measurements are untouched ===\n")
# --------------------------------------------------------------------------------------
unchanged <- sm[!as.logical(sm$sample_class_corrected), , drop = FALSE]
expect(nrow(unchanged) == 90L, sprintf("90 measurements are not corrected (got %d)",
                                       nrow(unchanged)))
expect(identical(unchanged$original_sample_class, unchanged$analysis_sample_class),
       "original == analysis sample class for all 90 uncorrected measurements")
expect(all(unchanged$sample_class_correction_status == "not_corrected"),
       "the 90 uncorrected rows are marked not_corrected")

# --------------------------------------------------------------------------------------
cat("\n=== both identities are retained, and neither is discarded ===\n")
# --------------------------------------------------------------------------------------
expect(all(c("original_sample_class", "analysis_sample_class", "sample_class_corrected")
           %in% names(sm)),
       "sample_metadata.tsv carries original_sample_class and analysis_sample_class")
expect(nrow(sm) == inv$n_measurement_records &&
         all(nzchar(sm$original_sample_class)) && all(nzchar(sm$analysis_sample_class)),
       "both labels are populated for all 96 measurements")
expect(identical(as.character(sm$analysis_sample_class), as.character(sm$sample_class)),
       "sample_class remains the canonical analysis-time assignment")
expect(all(sm$original_sample_class %in% sample_classes) &&
         all(sm$analysis_sample_class %in% sample_classes),
       "both labels use the canonical sample-class vocabulary")
expect(identical(as.character(sm$sample_class_plate_layout_implied),
                 as.character(sm$original_sample_class)),
       "the plate layout independently reproduces the pre-correction class for all 96")

# --------------------------------------------------------------------------------------
cat("\n=== the correction is a CONFIRMED_INTENTIONAL_SWITCH_CORRECTION ===\n")
# --------------------------------------------------------------------------------------
expect(all(corr$correction_status == "confirmed_intentional_switch_correction"),
       "correction_status is confirmed_intentional_switch_correction on every row")
expect(all(corr$correction_status %in% release_sample_class_correction_statuses),
       "correction_status uses the sanctioned vocabulary")
expect(!any(grepl("unresolved", corr$correction_status, ignore.case = TRUE)),
       "no row is recorded as an unresolved discrepancy")

# --------------------------------------------------------------------------------------
cat("\n=== F. no quantitative value was changed ===\n")
# --------------------------------------------------------------------------------------
expect(all(toupper(as.character(corr$quantitative_values_changed)) == "FALSE"),
       "quantitative_values_changed is FALSE on every corrected row")
expect(isTRUE(RELEASE_SAMPLE_CLASS_CORRECTION$quantitative_values_changed == FALSE) &&
         RELEASE_SAMPLE_CLASS_CORRECTION$max_abs_numeric_difference == 0,
       "the contract records a maximum absolute numeric difference of 0")

# The acquisition identity of the six is unchanged: the run name still yields the raw file.
six_rows <- sm[as.logical(sm$sample_class_corrected), , drop = FALSE]
expect(all(sub("[.]d$", "", six_rows$raw_file_basename) == six_rows$sample_id),
       "the six retain their original acquisition identity")
expect(all(grepl(paste0("_", six_rows$plate_sample_number, "_", collapse = "|"),
                 six_rows$sample_id)),
       "each corrected row's acquisition name still carries its own plate sample number")

# --------------------------------------------------------------------------------------
cat("\n=== the preserved correction record is verified by hash, not mtime ===\n")
# --------------------------------------------------------------------------------------
ref <- release_verify_sample_class_correction_reference(DATA_ROOT, require_match = FALSE)
if (!ref$exists) {
  cat("  [skip]  the correction record is unreachable; hash not checked\n")
} else {
  expect(isTRUE(ref$matches),
         paste("the correction record matches the SHA256 the forensic audit recorded:",
               substr(ref$observed_sha256, 1, 16)))
  expect(all(as.character(corr$correction_reference_sha256) == ref$expected_sha256),
         "every published row cites that same SHA256")
  expect(all(nzchar(corr$correction_reference)),
         "every row names the preserved correction record")
  expect(identical(unique(as.character(corr$correction_reference_mtime)),
                   ref$observed_mtime),
         paste("the recorded mtime matches the file:", ref$observed_mtime))
}
expect(isFALSE(RELEASE_SAMPLE_CLASS_CORRECTION$prose_rationale_exists),
       "the contract records that no prose rationale survives")

# --------------------------------------------------------------------------------------
cat("\n=== UMAP is documented, not treated as ground truth ===\n")
# --------------------------------------------------------------------------------------
if ("umap_nearest_class_suggestion" %in% names(corr)) {
  agree <- corr$umap_nearest_class_suggestion == corr$analysis_sample_class
  expect(sum(agree) == 5L && sum(!agree) == 1L,
         sprintf("the UMAP nearest-centre suggestion agrees for 5 of 6 (got %d)",
                 sum(agree)))
  expect(identical(corr$plate_sample_number[!agree], "N60"),
         "the divergent row is N60, and it is published rather than hidden")
  expect(corr$umap_nearest_class_suggestion[!agree] == "cfos" &&
           corr$analysis_sample_class[!agree] == "neuron",
         "N60: UMAP suggested cfos, the applied correction was neuron")
} else {
  expect(FALSE, "the UMAP suggestion is published for comparison")
}

# --------------------------------------------------------------------------------------
cat("\n=== K. the correction is not a release or PRIDE SKIP ===\n")
# --------------------------------------------------------------------------------------
vr_path <- file.path(OUT_ROOT, "provenance", "validation_results.tsv")
if (!file.exists(vr_path)) {
  cat("  [skip]  no validation report; run 13_validate_release.R\n")
} else {
  vr <- read_tsv(vr_path)
  sample_class_rows <- vr[grepl("sample-class|sample_class", vr$check), , drop = FALSE]
  expect(nrow(sample_class_rows) > 0L, "the validator reports on the sample-class correction")
  expect(!any(sample_class_rows$status == "SKIP"),
         sprintf("no sample-class check is a SKIP (%d row(s) checked)",
                 nrow(sample_class_rows)))
  expect(!any(sample_class_rows$status == "FAIL"), "no sample-class check FAILs")
  expect(any(sample_class_rows$status == "PASS" &
               grepl("RESOLVED", sample_class_rows$check, fixed = TRUE)),
         "the correction is recorded as RESOLVED and PASSing")

  # L. the genuinely missing acquisition metadata must stay unresolved.
  pride_missing <- vr[vr$check == "missing required SDRF fields", , drop = FALSE]
  if (nrow(pride_missing) == 1L && pride_missing$status == "SKIP") {
    expect(!grepl("sample class", pride_missing$detail, fixed = TRUE),
           "the remaining PRIDE SKIP is acquisition metadata, not sample class")
    cat("  [ok]   the missing SDRF acquisition metadata remains truthfully unresolved\n")
  } else if (nrow(pride_missing) == 1L && pride_missing$status == "PASS") {
    cat("  [ok]   all SDRF-required fields are now populated\n")
  } else {
    expect(FALSE, "the missing-SDRF-fields contract is reported exactly once")
  }
}

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Publication release sample-class correction contracts failed: %d", failures),
       call. = FALSE)
}
cat("All publication release sample-class correction contracts hold.\n")
