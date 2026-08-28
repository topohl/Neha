args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("tests", "test_pca_animal_level.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
source(file.path(repo_root, "R", "project_path_utils.R"))
if (!file.exists(file.path(repo_root, "R", "pca_animal_level_utils.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "protigy_input_utils.R"))
source(file.path(repo_root, "R", "pca_animal_level_utils.R"))

expect_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

expect_identical <- function(actual, expected, message) {
  if (!identical(actual, expected)) stop(message, call. = FALSE)
}

expect_equal <- function(actual, expected, message, tolerance = 1e-12) {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance, check.attributes = FALSE))) {
    stop(message, call. = FALSE)
  }
}

expect_error <- function(expression, pattern, message) {
  observed <- tryCatch({
    force(expression)
    NA_character_
  }, error = function(e) conditionMessage(e))
  if (is.na(observed) || !grepl(pattern, observed, ignore.case = TRUE)) {
    stop(message, " Observed: ", ifelse(is.na(observed), "<no error>", observed), call. = FALSE)
  }
}

animal_ids_by_condition <- list(
  paired_cno = c("C11", "C12", "C34f"),
  paired_veh = c("C14", "C33f", "C510"),
  unpaired_cno = c("C25", "C45f", "C46"),
  unpaired_veh = c("C26", "C27", "C47")
)
fixture_metadata <- do.call(rbind, lapply(sample_classes, function(sample_class) {
  do.call(rbind, lapply(names(animal_ids_by_condition), function(condition) {
    ids <- animal_ids_by_condition[[condition]]
    data.frame(
      output_column_name = paste(ids, sample_class, sep = "_"),
      AnimalID = ids,
      condition_code = names(condition_code_map)[match(condition, condition_code_map)],
      condition = condition,
      sample_class = sample_class,
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(fixture_metadata) <- NULL

fixture_matrix <- outer(
  seq_len(8L),
  seq_len(nrow(fixture_metadata)),
  function(protein, sample) protein * 0.3 + sample * 0.01 + ((protein * sample) %% 7) * 0.02
)
fixture_matrix[8, ] <- 5
rownames(fixture_matrix) <- paste0("P", seq_len(nrow(fixture_matrix)), "_MOUSE")
colnames(fixture_matrix) <- fixture_metadata$output_column_name
fixture_descriptions <- paste0("Gene", seq_len(nrow(fixture_matrix)))

fixture_path <- tempfile("animal_pca_", fileext = ".gct")
write_protigy_gct_v13(fixture_matrix, fixture_metadata, fixture_descriptions, fixture_path)
parsed <- validate_protigy_gct_v13(fixture_path, expected_matrix = fixture_matrix)

expect_identical(unname(parsed$dimensions[c("n_rows", "n_samples", "n_row_descriptors", "n_column_metadata")]),
                 c(8L, 48L, 1L, 5L),
                 "Robust GCT v1.3 parsing did not honor declared dimensions and metadata counts.")
expect_identical(rownames(parsed$column_metadata),
                 c("AnimalID", "condition_code", "condition", "sample_class", "phenotypeWithinUnit"),
                 "Robust parser did not retain the five named animal-level metadata rows.")
expect_identical(parsed$protein_ids, rownames(fixture_matrix),
                 "Robust parser changed protein identifiers.")

validated <- validate_pca_animal_input(parsed, expected_n = 3L)
expect_identical(ncol(validated$expression_matrix), 48L,
                 "Animal-level PCA input must contain exactly 48 observations.")
expect_identical(nrow(validated$sample_metadata), 48L,
                 "Animal-level PCA metadata must contain exactly 48 observations.")
expect_true(!anyDuplicated(paste(validated$sample_metadata$AnimalID,
                                 validated$sample_metadata$sample_class,
                                 sep = "\r")),
            "More than one observation entered an AnimalID x sample_class PCA unit.")
expect_true(all(validated$group_counts$n_animals == 3L),
            "Every sample_class x condition must contain exactly three animals.")
expect_identical(rownames(validated$sample_metadata), colnames(validated$expression_matrix),
                 "Animal-level PCA metadata is not aligned to abundance columns.")
expect_identical(validated$post_aggregation_transformations, "none",
                 "Animal-level PCA input validation introduced a transformation.")

prepared <- prepare_animal_pca(validated$expression_matrix, center = TRUE, scale. = TRUE)
expect_identical(prepared$removed_rows$protein_id, "P8_MOUSE",
                 "PCA did not remove exactly the zero-variance fixture protein.")
expect_identical(prepared$removed_rows$reason, "zero_variance_in_scaled_pca",
                 "Zero-variance removal reason is not audited explicitly.")
expect_equal(
  unname(prepared$matrix),
  unname(fixture_matrix[setdiff(rownames(fixture_matrix), "P8_MOUSE"), , drop = FALSE]),
  "PCA input abundances changed beyond zero-variance removal."
)
expect_identical(prepared$center, TRUE, "Primary PCA centering changed.")
expect_identical(prepared$scale, TRUE, "Primary PCA scaling changed.")
expect_identical(prepared$abundance_normalization, "none", "PCA added abundance normalization.")
expect_identical(prepared$abundance_imputation, "none", "PCA added abundance imputation.")
expect_identical(rownames(prepared$pca$x), colnames(prepared$matrix),
                 "PCA scores are not aligned to animal-level observations.")

prepared_repeat <- prepare_animal_pca(validated$expression_matrix, center = TRUE, scale. = TRUE)
expect_equal(prepared_repeat$pca$x, prepared$pca$x, "PCA scores are not deterministic for fixed input.", tolerance = 0)
expect_equal(prepared_repeat$pca$rotation, prepared$pca$rotation,
             "PCA loadings are not deterministic for fixed input.", tolerance = 0)
expect_equal(prepared_repeat$pca$sdev, prepared$pca$sdev,
             "PCA variance is not deterministic for fixed input.", tolerance = 0)

duplicate_fixture <- parsed
duplicate_fixture$column_metadata["AnimalID", 2] <- duplicate_fixture$column_metadata["AnimalID", 1]
expect_error(
  validate_pca_animal_input(duplicate_fixture),
  "duplicated within sample_class",
  "Duplicated AnimalID/sample_class observations must fail closed."
)

hemisphere_fixture <- parsed
hemisphere_fixture$matrix <- cbind(hemisphere_fixture$matrix, hemisphere_fixture$matrix[, 1])
colnames(hemisphere_fixture$matrix)[ncol(hemisphere_fixture$matrix)] <- paste0(colnames(parsed$matrix)[1], "_right")
hemisphere_fixture$column_metadata <- cbind(
  hemisphere_fixture$column_metadata,
  hemisphere_fixture$column_metadata[, 1]
)
colnames(hemisphere_fixture$column_metadata)[ncol(hemisphere_fixture$column_metadata)] <-
  colnames(hemisphere_fixture$matrix)[ncol(hemisphere_fixture$matrix)]
expect_error(
  validate_pca_animal_input(hemisphere_fixture),
  "duplicated within sample_class",
  "Hemisphere rows must not enter biological PCA as independent observations."
)

misaligned_fixture <- parsed
colnames(misaligned_fixture$column_metadata) <- rev(colnames(misaligned_fixture$column_metadata))
expect_error(
  validate_pca_animal_input(misaligned_fixture),
  "not aligned",
  "Misaligned GCT abundance and metadata columns must fail closed."
)

underfilled_fixture <- parsed
underfilled_fixture$matrix <- underfilled_fixture$matrix[, -1, drop = FALSE]
underfilled_fixture$column_metadata <- underfilled_fixture$column_metadata[, -1, drop = FALSE]
expect_error(
  validate_pca_animal_input(underfilled_fixture),
  "exactly 3 animals per sample_class/condition",
  "A sample_class/condition with fewer than three animals must fail closed."
)

nonfinite_fixture <- parsed
nonfinite_fixture$matrix[1, 1] <- NA_real_
expect_error(
  validate_pca_animal_input(nonfinite_fixture),
  "nonfinite values",
  "Nonfinite abundances must fail rather than trigger PCA imputation or filtering."
)

historical_root <- file.path(tempdir(), "historical_pca")
expect_error(
  validate_pca_output_root(file.path(historical_root, "child"), historical_root),
  "cannot overwrite the historical PCA root",
  "Animal-level PCA must reject the historical output root and descendants."
)
isolated_root <- validate_pca_output_root(file.path(tempdir(), "animal_pca"), historical_root)
expect_true(!pca_path_is_within(isolated_root, historical_root),
            "An isolated animal-level PCA output root was rejected or rewritten incorrectly.")

audit <- make_pca_audit(
  validated_input = validated,
  prepared_pca = prepared,
  source_path = fixture_path,
  source_sha256 = paste(rep("a", 64L), collapse = ""),
  output_paths = c(sample_class_plot = file.path(isolated_root, "pca_by_sample_class.svg")),
  execution_status = "success"
)
required_audit_columns <- c(
  "source_gct_path", "source_gct_sha256", "n_protein_rows_loaded", "n_samples_loaded",
  "n_samples_entering_pca", "sample_classes", "conditions", "n_animals",
  "pca_center", "pca_scale", "n_rows_removed_for_pca", "rows_removed_for_pca",
  "row_removal_reason", "pc1_variance_explained", "pc2_variance_explained",
  "output_paths", "execution_status", "error_message"
)
expect_identical(nrow(audit), 16L, "PCA audit must record all 16 sample_class x condition groups.")
expect_true(all(required_audit_columns %in% names(audit)), "PCA audit is missing required fields.")
expect_true(all(audit$n_samples_loaded == 48L & audit$n_samples_entering_pca == 48L),
            "PCA audit does not retain all 48 animal-level observations.")
expect_true(all(audit$n_animals == 3L), "PCA audit does not record n=3 for every group.")
expect_true(all(audit$n_rows_removed_for_pca == 1L & audit$row_removal_reason == "zero_variance_in_scaled_pca"),
            "PCA audit does not record exact zero-variance removal.")
expect_true(all(is.finite(audit$pc1_variance_explained) & is.finite(audit$pc2_variance_explained)),
            "PCA audit does not record PC1/PC2 variance explained.")
failure_audit <- make_pca_audit(
  source_path = fixture_path,
  source_sha256 = paste(rep("a", 64L), collapse = ""),
  output_paths = c(pca_audit = file.path(isolated_root, "animal_level_pca_audit.csv")),
  execution_status = "failed",
  error_message = "fixture failure"
)
expect_identical(nrow(failure_audit), 1L, "A failed PCA input must produce one machine-readable audit row.")
expect_identical(failure_audit$execution_status, "failed", "Failed PCA audit status is incorrect.")
expect_identical(failure_audit$error_message, "fixture failure", "Failed PCA audit did not retain the error message.")

# The PCA workflow was split (2026-08-26) from one monolith into an orchestrator plus ordered
# parts under 03_qc_exploration/pca/. These contract assertions are about the workflow as a
# whole, so read the entry point together with every part rather than just the entry point.
pca_script_files <- c(
  file.path(repo_root, "03_qc_exploration", "06_pcaPlot_animal_level.r"),
  sort(list.files(file.path(repo_root, "03_qc_exploration", "pca"),
                  pattern = "\\.r$", full.names = TRUE))
)
expect_true(length(pca_script_files) >= 2,
            "PCA workflow parts under 03_qc_exploration/pca/ were not found.")
pca_script <- paste(
  unlist(lapply(pca_script_files, function(f) readLines(f, warn = FALSE))),
  collapse = "\n"
)
expect_true(grepl("validate_protigy_gct_v13(gct_file)", pca_script, fixed = TRUE),
            "PCA script does not use the established robust GCT v1.3 parser.")
expect_true(!grepl("all_lines\\[4:12\\]|all_lines\\[13", pca_script),
            "PCA script still parses GCT metadata or proteins by fixed row positions.")
expect_true(grepl("PROTEOMICS_PCA_ANIMAL_LEVEL_INPUT", pca_script, fixed = TRUE) &&
              grepl("PROTEOMICS_PCA_OUTPUT_ROOT", pca_script, fixed = TRUE),
            "PCA script does not expose isolated input/output overrides.")
expect_true(grepl("plot_and_save_group(\"sample_class\"", pca_script, fixed = TRUE) &&
              grepl("plot_and_save_group(\"condition\"", pca_script, fixed = TRUE) &&
              grepl("plot_and_save_group(\"AnimalID\"", pca_script, fixed = TRUE),
            "Primary PCA plots do not expose sample_class, condition, and AnimalID.")
expect_true(!grepl("impute::impute.knn|row_medians|keep_cols <- colMeans\\(is.na", pca_script),
            "Primary biological PCA still performs missingness filtering or imputation.")

unlink(fixture_path)
cat("Animal-level PCA/QC tests passed.\n")
