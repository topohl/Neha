#!/usr/bin/env Rscript

# Publication release, stage 01 -- authoritative sample metadata tables.
#
# Produces
#   metadata/sample_metadata.tsv              one row per acquisition/measurement (96)
#   metadata/animal_level_sample_metadata.tsv one row per AnimalID x sample_class (48)
#   metadata/metadata_field_provenance.tsv    per-field status, so nothing is guessed
#   metadata/sample_class_corrections.tsv     the 6 corrected sample-class assignments
#
# Everything is READ from validated artefacts. No value is inferred:
#
#   source_sample_assignment.csv          measurement -> animal assignment used to build
#                                         the locked animal-level GCT
#   collection_plate_provenance_audit.csv  acquisition-run-name tokens already parsed and
#                                         audited (instrument alias, date, LC method,
#                                         collection plate, well, injection index)
#   design_identifiability_audit.csv       independent 48-row animal-level table with
#                                         explicit left/right source samples
#   sample_annotation.xlsx                 the ONLY artefact carrying original `.d`
#                                         acquisition paths; sheet is named report.stats_
#   sample_info.xlsx                       legacy 2024 schema (ReplicateGroup, plate, ...)
#   the locked animal-level GCT             column names, verified against what we build
#
# The 48-row table is DERIVED from the 96-row table and then cross-checked against
# design_identifiability_audit.csv. Two independent constructions agreeing is the point;
# a single read would not catch a stale audit.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_validation.R"))
release_require("readxl", "digest")

DATA_ROOT <- release_data_root()
PROJECT_ROOT <- release_project_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/01_build_sample_metadata.R"

release_banner("stage 01 -- sample metadata")
release_log("release root: ", OUT_ROOT)

# --------------------------------------------------------------------------------------
# canonical inputs
# --------------------------------------------------------------------------------------

INFERENTIAL <- file.path(DATA_ROOT, "03_output", "inferential_checks",
                         "final_inferential_checks_20260826")
INPUT_GCT_DIR <- file.path(DATA_ROOT, "02_data", "animal_level", "input_gct")

paths <- list(
  source_assignment = file.path(INPUT_GCT_DIR, "source_sample_assignment.csv"),
  aggregation_summary = file.path(INPUT_GCT_DIR, "aggregation_summary.csv"),
  plate_provenance = file.path(INFERENTIAL, "collection_plate_provenance_audit.csv"),
  design_identifiability = file.path(INFERENTIAL, "design_identifiability_audit.csv"),
  sample_annotation = file.path(PROJECT_ROOT, "sample_annotation.xlsx"),
  sample_info = file.path(DATA_ROOT, "01_input", "metadata", "sample_info.xlsx"),
  animal_gct = release_locked_artefact_path("animal_level_input_gct", DATA_ROOT)
)
for (nm in names(paths)) release_assert_exists(paths[[nm]], nm)

experimenter_source <- release_experimenter_metadata_source_path(REPO_ROOT)
protocol_source <- release_sample_preparation_source_path(REPO_ROOT)
experimenter_metadata <- release_read_experimenter_metadata(REPO_ROOT)
sample_preparation_protocol <- release_read_sample_preparation_protocol(REPO_ROOT)
release_log("  experimenter metadata: ", nrow(experimenter_metadata),
            " records; protocol: ", nrow(sample_preparation_protocol), " ordered steps")

assignment <- release_read_csv(paths$source_assignment)
plate <- release_read_csv(paths$plate_provenance)
design48 <- release_read_csv(paths$design_identifiability)
agg_summary <- release_read_csv(paths$aggregation_summary)
annotation <- as.data.frame(
  readxl::read_excel(paths$sample_annotation, sheet = 1, .name_repair = "minimal"),
  stringsAsFactors = FALSE, check.names = FALSE)
sample_info <- as.data.frame(
  readxl::read_excel(paths$sample_info, sheet = 1, .name_repair = "minimal"),
  stringsAsFactors = FALSE, check.names = FALSE)
animal_gct <- release_read_gct(paths$animal_gct)

stopifnot(nrow(assignment) == RELEASE_DESIGN_INVARIANTS$n_measurement_records)
stopifnot(nrow(plate) == RELEASE_DESIGN_INVARIANTS$n_measurement_records)
stopifnot(nrow(annotation) == RELEASE_DESIGN_INVARIANTS$n_measurement_records)
stopifnot(nrow(sample_info) == RELEASE_DESIGN_INVARIANTS$n_measurement_records)
stopifnot(nrow(design48) == RELEASE_DESIGN_INVARIANTS$n_animal_level_units)

# --------------------------------------------------------------------------------------
# raw acquisition file identity
# --------------------------------------------------------------------------------------
# sample_annotation.xlsx `Name` holds the original acquisition path, e.g.
#   D:\Proteomics\Fabian\Tobias\Olive_20241217_..._1_8782.d
# Stripping the directory and the `.d` suffix must reproduce the sample id exactly --
# that equality is what licenses publishing a raw_file_basename at all. It is asserted,
# not assumed; if it ever fails the stage stops rather than emitting a plausible name.

annotation_name <- as.character(annotation$Name)
raw_basename_no_ext <- sub("\\.d$", "", sub("^.*[\\\\/]", "", annotation_name), ignore.case = TRUE)
annotation_by_nnumber <- stats::setNames(seq_len(nrow(annotation)), as.character(annotation$id))

info_idx <- annotation_by_nnumber[as.character(sample_info$sampleNumber)]
if (anyNA(info_idx)) {
  stop("sample_annotation.xlsx does not cover every sample_info N-number.", call. = FALSE)
}
if (!identical(raw_basename_no_ext[info_idx], as.character(sample_info$id))) {
  stop("Original `.d` acquisition names do not reproduce the sample ids; ",
       "refusing to publish raw_file_basename.", call. = FALSE)
}
release_log("  raw acquisition identity verified for all ", nrow(sample_info),
            " runs (basename(Name) minus .d == sample_id)")

raw_by_sample <- data.frame(
  sample_id = as.character(sample_info$id),
  raw_file_basename = paste0(raw_basename_no_ext[info_idx], ".d"),
  raw_file_original_path = annotation_name[info_idx],
  raw_file_extension = ".d",
  historical_group_label = as.character(annotation$group_label)[info_idx],
  annotation_plate = as.character(annotation$Plate)[info_idx],
  historical_hemisphere_label = as.character(annotation$label)[info_idx],
  legacy_replicate_group = as.character(sample_info$ReplicateGroup),
  legacy_sample_number = as.character(sample_info$sampleNumber),
  legacy_shortname = as.character(sample_info$shortname),
  legacy_row_name = as.character(sample_info$row.names),
  legacy_group2 = as.character(sample_info$group2),
  stringsAsFactors = FALSE, check.names = FALSE
)

# --------------------------------------------------------------------------------------
# measurement-level table
# --------------------------------------------------------------------------------------

plate_idx <- match(assignment$sample_id, plate$sample_id)
if (anyNA(plate_idx)) stop("collection-plate audit does not cover every source sample.", call. = FALSE)
raw_idx <- match(assignment$sample_id, raw_by_sample$sample_id)
if (anyNA(raw_idx)) stop("acquisition-name table does not cover every source sample.", call. = FALSE)

cond <- release_split_condition(assignment$condition)
acq_date_token <- as.character(plate$date_token)[plate_idx]
acq_date_iso <- ifelse(
  grepl("^[0-9]{8}$", acq_date_token),
  paste(substr(acq_date_token, 1, 4), substr(acq_date_token, 5, 6), substr(acq_date_token, 7, 8),
        sep = "-"),
  NA_character_
)

measurement <- data.frame(
  sample_id = as.character(assignment$sample_id),
  acquisition_run_name = as.character(assignment$sample_id),
  raw_file_basename = raw_by_sample$raw_file_basename[raw_idx],
  raw_file_original_path = raw_by_sample$raw_file_original_path[raw_idx],
  raw_file_extension = raw_by_sample$raw_file_extension[raw_idx],
  AnimalID = as.character(assignment$AnimalID),
  hemisphere = as.character(assignment$hemisphere),
  hemisphere_evidence_source = as.character(assignment$hemisphere_crosscheck_mode),
  historical_hemisphere_label = raw_by_sample$historical_hemisphere_label[raw_idx],
  legacy_replicate_group = raw_by_sample$legacy_replicate_group[raw_idx],
  sample_class = as.character(assignment$sample_class),
  sample_class_historical_alias = release_sample_class_alias(assignment$sample_class),
  sample_class_historical_group_label = raw_by_sample$historical_group_label[raw_idx],
  condition_code = as.character(assignment$condition_code),
  condition = as.character(assignment$condition),
  pairing_status = cond$pairing_status,
  treatment = cond$treatment,
  collection_plate = as.character(plate$plate_token)[plate_idx],
  plate_sample_number = as.character(plate$N_token)[plate_idx],
  injection_index = as.integer(plate$N_number)[plate_idx],
  plate_set_token = as.character(plate$S_token)[plate_idx],
  well_position = as.character(plate$well_position)[plate_idx],
  acquisition_run_trailing_id = as.character(plate$trailing_id)[plate_idx],
  instrument_alias_token = as.character(plate$instrument_token)[plate_idx],
  instrument_model = NA_character_,
  acquisition_date = acq_date_iso,
  acquisition_date_token = acq_date_token,
  lc_and_method_token = as.character(plate$method_token)[plate_idx],
  analysis_unit_primary = "animal",
  included_in_animal_level = !as.logical(assignment$exclude),
  animal_level_sample_id = as.character(assignment$output_column_name),
  measurement_level_matrix_column = as.character(assignment$sample_id),
  source_assignment_status = as.character(assignment$source_assignment_status),
  legacy_sample_number = raw_by_sample$legacy_sample_number[raw_idx],
  legacy_shortname = raw_by_sample$legacy_shortname[raw_idx],
  legacy_group2_hemisphere_coded = raw_by_sample$legacy_group2[raw_idx],
  stringsAsFactors = FALSE, check.names = FALSE
)

measurement <- measurement[order(measurement$sample_class, measurement$condition_code,
                                 measurement$AnimalID, measurement$hemisphere), , drop = FALSE]
rownames(measurement) <- NULL

# --------------------------------------------------------------------------------------
# sample-class correction: both identities are preserved
# --------------------------------------------------------------------------------------
# Six acquisitions carry an analysis-time sample class that differs from every
# pre-correction record. The forensic audit resolved what that is: a deliberate correction
# applied during historical sample-identity QC, not a mislabelling and not an open
# question. Only class METADATA was reassigned -- the acquisition identities and the
# quantitative abundance profiles are the same bytes either way.
#
# This block therefore publishes BOTH identities for all 96 measurements rather than
# choosing one. `sample_class` stays the canonical analysis-time assignment, because that
# is what the locked GCTs and every validated statistic were computed on; changing it
# would desynchronise the published metadata from the results. `original_sample_class`
# carries the pre-correction assignment, so nothing is discarded.
#
# Three independent pre-correction records exist and are compared here:
#
#   1. sample_annotation.xlsx `group_label`  -- the acquisition-era annotation (all 96)
#   2. the autosampler plate layout          -- implied by the well row letter (all 96)
#   3. the retained UMAP correction table    -- the six corrected rows only
#
# The correction is only accepted if all three agree with each other AND the resulting six
# rows equal RELEASE_SAMPLE_CLASS_CORRECTION. Deriving it three ways and then requiring the
# contract to match is the point: a single read would let a stale source pass unnoticed.

CORRECTION <- RELEASE_SAMPLE_CLASS_CORRECTION

# The preserved correction record is cited by hash in the released provenance table, so the
# hash is computed here and required to match what the forensic audit recorded. It is never
# taken on trust from the narrative, and mtime alone is never treated as sufficient.
correction_ref <- release_verify_sample_class_correction_reference(DATA_ROOT)
release_log("  correction reference verified: ", basename(correction_ref$path))
release_log("    sha256 ", correction_ref$observed_sha256)
release_log("    mtime  ", correction_ref$observed_mtime,
            " (audit recorded ", correction_ref$expected_mtime, ")")

# Record 1: the acquisition-era annotation, normalised to the canonical vocabulary.
original_sample_class <- release_normalize_annotation_label(
  raw_by_sample$historical_group_label[match(measurement$sample_id, raw_by_sample$sample_id)])
if (anyNA(original_sample_class)) {
  stop("Could not normalise every sample_annotation.xlsx group_label to a canonical ",
       "sample class: ",
       paste(unique(raw_by_sample$historical_group_label[
         match(measurement$sample_id[is.na(original_sample_class)], raw_by_sample$sample_id)]),
         collapse = ", "),
       ". Refusing to run a comparison that would silently pass on NA.",
       call. = FALSE)
}

# Record 2: the autosampler plate layout. The modal class per well row is computed from the
# ANALYSIS class, so the layout rule is not borrowed from the source it is used to check;
# with only 6 of 96 reassigned the mode still recovers the pre-correction layout.
well_row <- substr(measurement$well_position, 1, 1)
row_modal_class <- vapply(split(measurement$sample_class, well_row), function(x) {
  names(sort(table(x), decreasing = TRUE))[[1]]
}, character(1))
layout_class <- unname(row_modal_class[well_row])

measurement$original_sample_class <- original_sample_class
measurement$analysis_sample_class <- measurement$sample_class
measurement$sample_class_plate_layout_implied <- layout_class
measurement$sample_class_corrected <-
  measurement$original_sample_class != measurement$analysis_sample_class

plate_position <- paste(measurement$plate_set_token, measurement$well_position, sep = "-")
corrected_idx <- which(measurement$sample_class_corrected)

# --------------------------------------------------------------------------------------
# record 3: the retained UMAP correction table
# --------------------------------------------------------------------------------------
# Parsed for the pre-correction class it recorded per run. Its `potential_switch_celltype`
# column is the UMAP nearest-class-centre SUGGESTION and is deliberately NOT used as the
# corrected class: for N60 the suggestion was cFosN while the correction actually applied
# was neuron. UMAP identified the switch; it is not treated as ground truth for what the
# switch was. The applied classes come from the validated analysis assignment, and the
# independent profile-similarity and left/right-pair analyses support them.

correction_table <- utils::read.csv(correction_ref$path, sep = ";", stringsAsFactors = FALSE,
                                    check.names = FALSE)
correction_table <- correction_table[nzchar(trimws(correction_table$outlier)), , drop = FALSE]
ct_n_token <- sub("^.*_(N[0-9]{2})_S[45]-[A-H][0-9].*$", "\\1", correction_table$outlier)
ct_original <- release_normalize_annotation_label(correction_table$outlier_celltype)
ct_suggested <- release_normalize_annotation_label(correction_table$potential_switch_celltype)
names(ct_original) <- ct_n_token
names(ct_suggested) <- ct_n_token

release_log("  correction table lists ", nrow(correction_table),
            " historical outlier acquisition(s); ", length(intersect(ct_n_token,
            CORRECTION$expected$plate_sample_number)), " of them are the corrected six")

# --------------------------------------------------------------------------------------
# assemble the correction provenance table and require it to equal the contract
# --------------------------------------------------------------------------------------

correction_reference_relpath <- paste(c("clusterProfiler", CORRECTION$reference_relative),
                                      collapse = "/")

sample_class_corrections <- data.frame(
  sample_id = measurement$sample_id[corrected_idx],
  plate_sample_number = measurement$plate_sample_number[corrected_idx],
  AnimalID = measurement$AnimalID[corrected_idx],
  hemisphere = measurement$hemisphere[corrected_idx],
  plate_position = plate_position[corrected_idx],
  collection_plate = measurement$collection_plate[corrected_idx],
  condition = measurement$condition[corrected_idx],
  original_sample_class = measurement$original_sample_class[corrected_idx],
  analysis_sample_class = measurement$analysis_sample_class[corrected_idx],
  correction_status = CORRECTION$status,
  correction_method = CORRECTION$method,
  correction_date = CORRECTION$correction_date,
  correction_reference = correction_reference_relpath,
  correction_reference_sha256 = correction_ref$observed_sha256,
  correction_reference_mtime = correction_ref$observed_mtime,
  quantitative_values_changed = FALSE,
  umap_nearest_class_suggestion = unname(ct_suggested[measurement$plate_sample_number[corrected_idx]]),
  original_class_corroborating_records = paste(
    "sample_annotation.xlsx group_label; autosampler plate layout;",
    "umap_outlier_samples_with_switches_CORRECTED.csv outlier_celltype"),
  animal_level_unit_affected = measurement$animal_level_sample_id[corrected_idx],
  notes = paste(
    "Sample-class metadata reassigned during historical sample-identity QC. Acquisition",
    "identity and quantitative protein-abundance profile unchanged. No surviving prose",
    "note records the rationale; the correction table itself is the preserved historical",
    "correction record."),
  stringsAsFactors = FALSE, check.names = FALSE
)
sample_class_corrections <- sample_class_corrections[
  order(match(sample_class_corrections$plate_sample_number,
              CORRECTION$expected$plate_sample_number)), , drop = FALSE]
rownames(sample_class_corrections) <- NULL

correction_checks <- release_check_sample_class_corrections(sample_class_corrections)
if (any(!correction_checks$pass)) {
  stop("The sample-class correction derived from the canonical sources does not match the ",
       "release contract:\n",
       paste("  -", correction_checks$check[!correction_checks$pass], collapse = "\n"),
       "\nRefusing to publish a correction record that disagrees with the forensic audit.",
       call. = FALSE)
}

# The plate layout must independently reproduce the pre-correction class, for all 96.
if (!identical(measurement$sample_class_plate_layout_implied,
               measurement$original_sample_class)) {
  bad <- which(measurement$sample_class_plate_layout_implied !=
                 measurement$original_sample_class)
  stop("The autosampler plate layout and sample_annotation.xlsx disagree about the ",
       "pre-correction sample class for ", length(bad), " measurement(s): ",
       paste(measurement$plate_sample_number[bad], collapse = ", "),
       ". Two pre-correction records that do not agree cannot corroborate each other.",
       call. = FALSE)
}

# The correction table must independently reproduce the pre-correction class for the six.
ct_for_six <- unname(ct_original[sample_class_corrections$plate_sample_number])
if (anyNA(ct_for_six)) {
  stop("The preserved correction table does not cover every corrected measurement: ",
       paste(sample_class_corrections$plate_sample_number[is.na(ct_for_six)], collapse = ", "),
       call. = FALSE)
}
if (!identical(ct_for_six, sample_class_corrections$original_sample_class)) {
  stop("The preserved correction table records a different pre-correction class than ",
       "sample_annotation.xlsx for at least one corrected measurement.", call. = FALSE)
}

# 90 of 96 must be untouched, and no non-corrected row may differ.
n_corrected <- nrow(sample_class_corrections)
n_unchanged <- sum(!measurement$sample_class_corrected)
if (n_corrected != RELEASE_N_SAMPLE_CLASS_CORRECTED ||
    n_unchanged != RELEASE_DESIGN_INVARIANTS$n_measurement_records -
      RELEASE_N_SAMPLE_CLASS_CORRECTED) {
  stop("Expected exactly ", RELEASE_N_SAMPLE_CLASS_CORRECTED, " corrected and ",
       RELEASE_DESIGN_INVARIANTS$n_measurement_records - RELEASE_N_SAMPLE_CLASS_CORRECTED,
       " unchanged sample classes; got ", n_corrected, " and ", n_unchanged, ".",
       call. = FALSE)
}
unchanged <- measurement[!measurement$sample_class_corrected, , drop = FALSE]
if (!identical(unchanged$original_sample_class, unchanged$analysis_sample_class)) {
  stop("A measurement not flagged as corrected nonetheless has differing original and ",
       "analysis sample classes.", call. = FALSE)
}

# Per-row correction provenance. Blank rather than repeated boilerplate on the 90 rows
# where nothing was corrected, so a reader filtering on these fields gets the six.
blank_if_uncorrected <- function(value) {
  ifelse(measurement$sample_class_corrected, value, NA_character_)
}
measurement$sample_class_correction_status <-
  ifelse(measurement$sample_class_corrected, CORRECTION$status, "not_corrected")
measurement$sample_class_correction_method <- blank_if_uncorrected(CORRECTION$method)
measurement$sample_class_correction_provenance <-
  blank_if_uncorrected(correction_reference_relpath)
measurement$sample_class_correction_reference_sha256 <-
  blank_if_uncorrected(correction_ref$observed_sha256)
measurement$sample_class_correction_date <- blank_if_uncorrected(CORRECTION$correction_date)

release_log("  sample class: ", n_unchanged, "/", nrow(measurement),
            " measurements unchanged, ", n_corrected, " corrected (",
            CORRECTION$status, ")")
for (i in seq_len(nrow(sample_class_corrections))) {
  r <- sample_class_corrections[i, , drop = FALSE]
  release_log("      ", r$plate_sample_number, " ", r$plate_position, " ", r$AnimalID,
              " ", r$hemisphere, ": ", r$original_sample_class, " -> ",
              r$analysis_sample_class)
}
release_log("  both original and analysis sample-class labels are retained for all ",
            nrow(measurement), " measurements")

# Cross-checks the design contracts depend on.
stopifnot(nrow(measurement) == RELEASE_DESIGN_INVARIANTS$n_measurement_records)
stopifnot(!anyDuplicated(measurement$sample_id))
stopifnot(!anyDuplicated(measurement$raw_file_basename))
stopifnot(all(measurement$hemisphere %in% c("Left", "Right")))
stopifnot(all(measurement$sample_class %in% sample_classes))
stopifnot(all(measurement$condition %in% condition_levels))
stopifnot(all(measurement$included_in_animal_level))
if (!identical(sort(unique(measurement$animal_level_sample_id)), sort(animal_gct$sample_ids))) {
  stop("animal_level_sample_id does not reproduce the locked GCT column names.", call. = FALSE)
}
# Legacy ReplicateGroup is the numeric hemisphere encoding; confirm it still agrees.
legacy_hemisphere <- ifelse(measurement$legacy_replicate_group == "1", "Left", "Right")
if (!identical(legacy_hemisphere, measurement$hemisphere)) {
  stop("Legacy ReplicateGroup disagrees with the audited hemisphere assignment.", call. = FALSE)
}
release_log("  measurement-level table: ", nrow(measurement), " rows x ", ncol(measurement), " cols")

# --------------------------------------------------------------------------------------
# animal-level table
# --------------------------------------------------------------------------------------

key <- paste(measurement$AnimalID, measurement$sample_class, sep = "\r")
split_by_unit <- split(seq_len(nrow(measurement)), key)

animal_level <- do.call(rbind, lapply(split_by_unit, function(idx) {
  block <- measurement[idx, , drop = FALSE]
  block <- block[order(block$hemisphere), , drop = FALSE]
  if (length(unique(block$condition)) != 1L) {
    stop("An AnimalID x sample_class unit spans more than one condition.", call. = FALSE)
  }
  if (length(unique(block$collection_plate)) != 1L) {
    stop("An AnimalID x sample_class unit spans more than one collection plate.", call. = FALSE)
  }
  data.frame(
    AnimalID = block$AnimalID[[1]],
    sample_class = block$sample_class[[1]],
    sample_class_historical_alias = block$sample_class_historical_alias[[1]],
    condition = block$condition[[1]],
    condition_code = block$condition_code[[1]],
    pairing_status = block$pairing_status[[1]],
    treatment = block$treatment[[1]],
    collection_plate = block$collection_plate[[1]],
    n_hemisphere_measurements = nrow(block),
    hemispheres_present = paste(block$hemisphere, collapse = ";"),
    source_sample_ids = paste(block$sample_id, collapse = ";"),
    source_raw_file_basenames = paste(block$raw_file_basename, collapse = ";"),
    left_source_sample = block$sample_id[block$hemisphere == "Left"][[1]],
    right_source_sample = block$sample_id[block$hemisphere == "Right"][[1]],
    animal_level_column_name = block$animal_level_sample_id[[1]],
    phenotype_within_unit = paste(block$sample_class[[1]], block$condition[[1]], sep = "_"),
    aggregation_policy = "equal_weight_mean_LR_on_existing_imputed_log2_values",
    analysis_unit = "animal",
    stringsAsFactors = FALSE, check.names = FALSE
  )
}))
animal_level <- animal_level[order(animal_level$sample_class, animal_level$condition_code,
                                   animal_level$AnimalID), , drop = FALSE]
rownames(animal_level) <- NULL

# Independent cross-check against the canonical 48-row identifiability audit.
d48 <- design48[order(design48$sample_class, design48$AnimalID), , drop = FALSE]
a48 <- animal_level[order(animal_level$sample_class, animal_level$AnimalID), , drop = FALSE]
for (pair in list(c("AnimalID", "AnimalID"), c("sample_class", "sample_class"),
                  c("condition", "condition"), c("left_source_sample", "left_source_sample"),
                  c("right_source_sample", "right_source_sample"))) {
  if (!identical(as.character(a48[[pair[[1]]]]), as.character(d48[[pair[[2]]]]))) {
    stop("Derived animal-level table disagrees with design_identifiability_audit.csv on ",
         pair[[1]], call. = FALSE)
  }
}
if (!identical(as.character(a48$collection_plate), as.character(d48$Plate))) {
  stop("Derived collection plate disagrees with design_identifiability_audit.csv.", call. = FALSE)
}
release_log("  animal-level table cross-checks against design_identifiability_audit.csv")

design_checks <- release_check_animal_level_design(animal_level)
if (any(!design_checks$pass)) {
  stop("Animal-level design invariants violated:\n",
       paste("  -", design_checks$check[!design_checks$pass], collapse = "\n"), call. = FALSE)
}
# The canonical aggregation summary independently records 3 animals and complete bilateral
# pairing in all 16 strata; require it to still say so.
if (!all(agg_summary$n_unique_animals == 3L) ||
    !all(agg_summary$complete_bilateral_pairs == 3L) ||
    !all(agg_summary$one_sided_pairs == 0L) || !all(agg_summary$missing_pairs == 0L)) {
  stop("Canonical aggregation_summary.csv no longer reports 3/3 complete bilateral pairs.",
       call. = FALSE)
}
release_log("  animal-level table: ", nrow(animal_level), " rows x ", ncol(animal_level), " cols")

# --------------------------------------------------------------------------------------
# per-field provenance -- what is known, what is not, and what would be needed
# --------------------------------------------------------------------------------------

fld <- function(field, status, source, evidence, note = NA_character_) {
  data.frame(field = field, status = status, current_source = source,
             evidence = evidence, note = note,
             stringsAsFactors = FALSE, check.names = FALSE)
}

field_provenance <- rbind(
  fld("sample_id", "KNOWN_VERIFIED", paths$source_assignment,
      "96 unique values; identical to the measurement-level GCT column names and to sample_info.xlsx id",
      "This string is the acquisition run name, not an arbitrary label."),
  fld("raw_file_basename", "KNOWN_VERIFIED", paths$sample_annotation,
      "basename(Name) with the .d suffix removed equals sample_id for all 96 runs (asserted in this stage)",
      "Name values are of the form D:\\Proteomics\\Fabian\\Tobias\\<run>.d"),
  fld("raw_file_original_path", "KNOWN_VERIFIED", paths$sample_annotation,
      "verbatim Name column",
      "Acquisition-machine path; the files themselves are not in the project tree."),
  fld("raw_file_extension", "KNOWN_VERIFIED", paths$sample_annotation,
      "every Name ends in .d",
      paste("`.d` is a vendor-specific acquisition directory format (Bruker).",
            "The instrument MODEL is not recorded anywhere in the project.")),
  fld("AnimalID", "KNOWN_VERIFIED", paths$source_assignment,
      "12 distinct values: C11 C12 C14 C25 C26 C27 C33f C34f C45f C46 C47 C510",
      NA_character_),
  fld("hemisphere", "KNOWN_VERIFIED", paths$sample_annotation,
      "explicit _left/_right suffix on the annotation label; numeric ReplicateGroup 1/2 agrees one-to-one (asserted in this stage)",
      "The 2024 acquisition run names contain no L/R token; hemisphere is NOT inferred from sample_id."),
  fld("sample_class", "KNOWN_VERIFIED", paths$source_assignment,
      paste0("4 classes; historical labels cFosN/Background/mCherryN/neuron normalise ",
             "through R/analysis_labels.R. Analysis-time class agrees with the ",
             "pre-correction records for ", n_unchanged, " of ", nrow(measurement),
             " measurements; the remaining ", n_corrected, " were corrected during ",
             "historical sample-identity QC and are documented row by row in ",
             "metadata/sample_class_corrections.tsv."),
      paste0("Historical alias `bg` == neuropil. This column is the ANALYSIS-TIME class, ",
             "i.e. the class the locked GCTs and every validated statistic were computed ",
             "on. The pre-correction assignment is retained for all 96 measurements in ",
             "`original_sample_class`; neither label is discarded. Marker definitions and ",
             "laser-capture microdissection provenance are in the checked-in experimenter ",
             "metadata record.")),
  fld("original_sample_class", "KNOWN_VERIFIED", paths$sample_annotation,
      paste0("sample_annotation.xlsx group_label normalised to the canonical vocabulary, ",
             "independently reproduced for all 96 measurements by the autosampler ",
             "plate layout and, for the ", n_corrected, " corrected rows, by the retained ",
             "UMAP correction table (asserted in this stage)."),
      "Pre-correction sample-class assignment. Retained so the original identity is not lost."),
  fld("sample_class_corrected", "KNOWN_VERIFIED",
      release_sample_class_correction_reference_path(DATA_ROOT),
      paste0("TRUE for exactly ", n_corrected, " of ", nrow(measurement),
             " measurements: ", paste(sample_class_corrections$plate_sample_number,
                                      collapse = " "),
             ", all left hemisphere, animals ",
             paste(sort(unique(sample_class_corrections$AnimalID)), collapse = " and "),
             ". Correction reference verified by SHA256 ",
             substr(correction_ref$observed_sha256, 1, 16), "..."),
      paste0("Status ", CORRECTION$status, ". Only sample-class metadata were reassigned; ",
             "acquisition identifiers and quantitative protein-abundance profiles were ",
             "unchanged. No surviving prose note records the rationale -- the correction ",
             "table itself is the preserved historical correction record.")),
  fld("condition / condition_code", "KNOWN_VERIFIED", paths$source_assignment,
      "ExpGroup 1..4 maps through condition_code_map to paired_cno/paired_veh/unpaired_cno/unpaired_veh",
      NA_character_),
  fld("collection_plate", "KNOWN_VERIFIED", paths$plate_provenance,
      "Plate1/Plate2 token parsed from the acquisition run name; agrees with sample_info plate",
      paste("Collection plate only. NOT a proteomics preparation, digestion, LC-MS,",
            "acquisition or instrument batch.")),
  fld("well_position / plate_set_token / injection_index", "KNOWN_VERIFIED", paths$plate_provenance,
      "S4/S5 set token, well A1..H7 and running index 1..96 parsed from the acquisition run name",
      NA_character_),
  fld("acquisition_date", "KNOWN_BUT_NEEDS_STANDARDIZATION", paths$plate_provenance,
      "date_token 20241217 in every acquisition run name; single acquisition date for all 96 runs",
      paste("Derived from the run-name token, not from an instrument record.",
            "Confirm against the acquisition log before submission.")),
  fld("instrument_alias_token", "KNOWN_BUT_NEEDS_STANDARDIZATION", paths$plate_provenance,
      "token `Olive` in every acquisition run name",
      paste("A local instrument alias, not a model. Do not publish it as an instrument model.")),
  fld("instrument_model", "MISSING_RECOVERABLE", "NONE",
      "no instrument model string exists anywhere in the project tree",
      "Needed from: the acquisition facility, an instrument method file, or the .d metadata."),
  fld("lc_and_method_token", "KNOWN_BUT_NEEDS_STANDARDIZATION", paths$plate_provenance,
      "token `FCo_Evo2_80SPDzoom_5cmRapid_Tobias` in every acquisition run name",
      paste("Preserved verbatim and deliberately NOT decoded. It is an operator method",
            "label; the acquisition method itself is not documented.")),
  fld("strain or breed", "KNOWN_VERIFIED", experimenter_source,
      release_metadata_value(experimenter_metadata, "strain_or_breed"),
      "Cohort-level metadata; applies to all 12 animals."),
  fld("animal supplier", "KNOWN_VERIFIED", experimenter_source,
      paste(release_metadata_value(experimenter_metadata, "animal_supplier"),
            release_metadata_value(experimenter_metadata, "animal_supplier_location"),
            sep = ", "),
      "Cohort-level provenance; not an AnimalID-level field."),
  fld("sex", "MISSING_RECOVERABLE", experimenter_source,
      release_metadata_value(experimenter_metadata, "cohort_sex_composition"),
      paste("The cohort contained both sexes, but reliable AnimalID-level sex assignments",
            "are unavailable. The SDRF therefore retains `not available` per sample.")),
  fld("age at experiment start", "KNOWN_VERIFIED", experimenter_source,
      release_metadata_value(experimenter_metadata, "age_at_experiment_start"),
      "Cohort-level range at experiment start; not an age-at-collection value."),
  fld("age at tissue collection", "MISSING_RECOVERABLE", experimenter_source,
      release_metadata_value(experimenter_metadata, "age_at_tissue_collection"),
      paste("No established timeline supports conversion from age at experiment start.",
            "The SDRF age field therefore remains `not available`.")),
  fld("organism", "KNOWN_VERIFIED", file.path(DATA_ROOT, "01_input", "references",
                                              "MOUSE_10090_idmapping.dat"),
      "mapping reference is MOUSE_10090_idmapping.dat; 02_id_mapping queries organism_id = 10090; org.Mm.eg.db used throughout; protein entry names end in _MOUSE",
      "Mus musculus, NCBI taxid 10090."),
  fld("organism part / brain region", "KNOWN_VERIFIED", experimenter_source,
      release_metadata_value(experimenter_metadata, "organism_part"),
      paste("Experimenter-supplied CeM provenance. The SDRF uses the supported parent",
            "UBERON term while retaining the exact CeM description in this record.")),
  fld("sample collection method", "KNOWN_VERIFIED", experimenter_source,
      release_metadata_value(experimenter_metadata, "collection_method"),
      "Applies to all four marker-defined sample classes."),
  fld("proteomics labelling", "KNOWN_VERIFIED", experimenter_source,
      release_metadata_value(experimenter_metadata, "labeling_strategy"),
      "Verified as experimenter metadata; not inferred from run counts."),
  fld("sample preparation", "KNOWN_VERIFIED", protocol_source,
      "seven ordered steps; protocol version 2024-11-28",
      paste("Full lysis, heating, Lys-C and trypsin digestion, stop/SpeedVac option,",
            "Evotip cleanup and peptide-resuspension record is released verbatim.")),
  fld("search modification parameters", "MISSING_RECOVERABLE", experimenter_source,
      release_metadata_value(experimenter_metadata, "search_modification_parameters"),
      "CAA/TCEP use is not treated as evidence of search modification settings."),
  fld("proteomics acquisition mode", "MISSING_RECOVERABLE", experimenter_source,
      release_metadata_value(experimenter_metadata, "acquisition_mode"),
      "Needed from the instrument method or acquisition record."),
  fld("technical_replicate", "NOT_APPLICABLE", paths$source_assignment,
      "each acquisition run is one biological hemisphere sample; no run was injected twice",
      "Left/Right are anatomical subsamples of one animal, not technical replicates.")
)
release_assert_metadata_status(field_provenance$status)

# --------------------------------------------------------------------------------------
# write
# --------------------------------------------------------------------------------------

canonical_sources <- c(paths$source_assignment, paths$plate_provenance,
                       paths$design_identifiability, paths$sample_annotation,
                       paths$sample_info, paths$aggregation_summary)
canonical_hashes <- vapply(canonical_sources, release_sha256, character(1), USE.NAMES = FALSE)

w1 <- release_write_table(measurement, release_path("metadata", "sample_metadata.tsv"))
release_register("metadata/sample_metadata.tsv", "sample metadata, one row per acquisition/measurement",
                 canonical_sources, canonical_hashes, STAGE, "tsv")

w2 <- release_write_table(animal_level, release_path("metadata", "animal_level_sample_metadata.tsv"))
release_register("metadata/animal_level_sample_metadata.tsv",
                 "animal-level inferential units, one row per AnimalID x sample_class",
                 canonical_sources, canonical_hashes, STAGE, "tsv")

w3 <- release_write_table(field_provenance,
                          release_path("metadata", "metadata_field_provenance.tsv"))
release_register("metadata/metadata_field_provenance.tsv",
                 "per-field metadata status: verified / needs standardisation / missing",
                 c(canonical_sources, experimenter_source, protocol_source),
                 c(canonical_hashes, release_sha256(experimenter_source),
                   release_sha256(protocol_source)), STAGE, "tsv")

w4 <- release_write_table(sample_class_corrections,
                          release_path("metadata", "sample_class_corrections.tsv"))
release_register("metadata/sample_class_corrections.tsv",
                 paste("the six measurements whose sample-class metadata were corrected",
                       "during historical sample-identity QC; both the original and the",
                       "analysis-time class, with the preserved correction record and its",
                       "SHA256"),
                 c(paths$sample_info, paths$sample_annotation, correction_ref$path),
                 c(release_sha256(paths$sample_info), release_sha256(paths$sample_annotation),
                   correction_ref$observed_sha256),
                 STAGE, "tsv")

w5 <- release_write_table(experimenter_metadata,
                          release_path("metadata", "experimenter_metadata.tsv"))
release_register("metadata/experimenter_metadata.tsv",
                 paste("source-controlled experimenter facts, scope exclusions and",
                       "explicitly unresolved deposition metadata"),
                 experimenter_source, release_sha256(experimenter_source), STAGE, "tsv")

w6 <- release_write_table(sample_preparation_protocol,
                          release_path("metadata", "sample_preparation_protocol.tsv"))
release_register("metadata/sample_preparation_protocol.tsv",
                 "source-controlled 2024-11-28 sample-preparation protocol",
                 protocol_source, release_sha256(protocol_source), STAGE, "tsv")

release_log("  wrote sample_metadata.tsv (", w1$rows, "x", w1$cols, ")")
release_log("  wrote animal_level_sample_metadata.tsv (", w2$rows, "x", w2$cols, ")")
release_log("  wrote metadata_field_provenance.tsv (", w3$rows, "x", w3$cols, ")")
release_log("  wrote sample_class_corrections.tsv (", w4$rows, "x", w4$cols, ")")
release_log("  wrote experimenter_metadata.tsv (", w5$rows, "x", w5$cols, ")")
release_log("  wrote sample_preparation_protocol.tsv (", w6$rows, "x", w6$cols, ")")
release_log("stage 01 complete")
