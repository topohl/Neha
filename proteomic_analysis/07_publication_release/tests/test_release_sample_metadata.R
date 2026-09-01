#!/usr/bin/env Rscript

# Contract: the published sample metadata describes the validated animal-level design and
# nothing else. 96 acquisitions, 48 animal-level units, 12 animals, 16 strata, n = 3.
#
# Skips cleanly when no built release is reachable, so it is safe in CI. Point it at one
# with PROTEOMICS_RELEASE_OUTPUT_ROOT.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("07_publication_release", "tests", "test_release_sample_metadata.R")
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
sm_path <- file.path(OUT_ROOT, "metadata", "sample_metadata.tsv")
al_path <- file.path(OUT_ROOT, "metadata", "animal_level_sample_metadata.tsv")
if (!file.exists(sm_path) || !file.exists(al_path)) {
  cat("No built publication release at", OUT_ROOT, "\n")
  cat("  set PROTEOMICS_RELEASE_OUTPUT_ROOT and run 07_publication_release/run_release.R.\n")
  cat("Skipping publication release sample-metadata contracts.\n")
  quit(save = "no", status = 0L)
}

read_tsv <- function(p) read.delim(p, sep = "\t", quote = "", stringsAsFactors = FALSE,
                                   check.names = FALSE)
sm <- read_tsv(sm_path)
al <- read_tsv(al_path)
inv <- RELEASE_DESIGN_INVARIANTS

cat("=== measurement-level records ===\n")
expect(nrow(sm) == inv$n_measurement_records,
       sprintf("96 measurement-level records (got %d)", nrow(sm)))
expect(!anyDuplicated(sm$sample_id), "sample_id is unique")
expect(all(sm$sample_class %in% sample_classes), "every sample_class is canonical")
expect(all(sm$condition %in% condition_levels), "every condition is canonical")
expect(all(sm$hemisphere %in% c("Left", "Right")), "hemisphere is Left or Right")
expect(all(sm$analysis_unit_primary == "animal"),
       "analysis_unit_primary is animal on every row")
expect(all(as.logical(sm$included_in_animal_level)),
       "every measurement contributed to an animal-level unit")
expect(all(sm$collection_plate %in% c("Plate1", "Plate2")),
       "collection_plate is Plate1 or Plate2")

cat("\n=== raw acquisition identity ===\n")
expect(!anyDuplicated(sm$raw_file_basename), "raw_file_basename is unique")
expect(all(grepl("[.]d$", sm$raw_file_basename)),
       "every raw_file_basename ends in .d")
expect(all(sub("[.]d$", "", sm$raw_file_basename) == sm$sample_id),
       "raw_file_basename minus .d equals sample_id")
expect(all(nzchar(sm$raw_file_original_path)),
       "every row carries the original acquisition path")
expect(all(is.na(sm$instrument_model) | !nzchar(sm$instrument_model)),
       "instrument_model is left blank rather than inferred")
expect(all(sm$instrument_alias_token == "Olive"),
       "the run-name instrument alias is recorded as a token, not as a model")

cat("\n=== animal-level units ===\n")
al$n_hemisphere_measurements <- as.integer(al$n_hemisphere_measurements)
expect(nrow(al) == inv$n_animal_level_units,
       sprintf("48 animal-level units (got %d)", nrow(al)))
expect(length(unique(al$AnimalID)) == inv$n_animals,
       sprintf("12 distinct AnimalIDs (got %d)", length(unique(al$AnimalID))))
expect(length(unique(al$sample_class)) == inv$n_sample_classes, "4 sample classes")
expect(length(unique(al$condition)) == inv$n_conditions, "4 conditions")
expect(nrow(unique(al[c("sample_class", "condition")])) == inv$n_strata,
       "16 sample_class x condition strata")
expect(all(table(al$sample_class, al$condition) == inv$n_animals_per_stratum),
       "3 animals in every stratum")

design_checks <- release_check_animal_level_design(al)
for (i in seq_len(nrow(design_checks))) {
  expect(design_checks$pass[i], design_checks$check[i])
}

cat("\n=== the two tables agree ===\n")
expect(all(sm$animal_level_sample_id %in% al$animal_level_column_name),
       "every measurement maps to a published animal-level unit")
expect(setequal(unique(sm$animal_level_sample_id), al$animal_level_column_name),
       "the animal-level units are exactly those the measurements produce")
per_unit <- table(sm$animal_level_sample_id)
expect(all(per_unit == inv$n_hemispheres_per_unit),
       "exactly 2 measurements feed every animal-level unit")
expect(all(vapply(seq_len(nrow(al)), function(i) {
  ids <- strsplit(al$source_sample_ids[i], ";", fixed = TRUE)[[1]]
  length(ids) == 2L && all(ids %in% sm$sample_id)
}, logical(1))), "source_sample_ids resolve to published measurements")
expect(all(vapply(seq_len(nrow(al)), function(i) {
  sub <- sm[sm$animal_level_sample_id == al$animal_level_column_name[i], , drop = FALSE]
  length(unique(sub$condition)) == 1L && unique(sub$condition) == al$condition[i]
}, logical(1))), "each unit has one condition and it matches its measurements")

cat("\n=== metadata is not invented ===\n")
fp_path <- file.path(OUT_ROOT, "metadata", "metadata_field_provenance.tsv")
if (!file.exists(fp_path)) {
  expect(FALSE, "metadata_field_provenance.tsv exists")
} else {
  fp <- read_tsv(fp_path)
  expect(all(fp$status %in% release_metadata_statuses),
         "every field carries a sanctioned status value")
  expect(any(fp$status == "MISSING_RECOVERABLE" | fp$status == "MISSING_UNKNOWN"),
         "fields that are genuinely unknown are recorded as missing")
  known <- fp[fp$status %in% c("KNOWN_VERIFIED", "KNOWN_BUT_NEEDS_STANDARDIZATION"), ,
              drop = FALSE]
  expect(all(nzchar(known$evidence)), "every KNOWN_* field cites its evidence")
  im <- fp[fp$field == "instrument_model", , drop = FALSE]
  expect(nrow(im) == 1L && im$status == "MISSING_RECOVERABLE",
         "instrument_model is recorded as missing, not guessed")
}

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Publication release sample-metadata contracts failed: %d", failures),
       call. = FALSE)
}
cat("All publication release sample-metadata contracts hold.\n")
