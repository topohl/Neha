args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("tests", "test_ewce_animal_level.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
source(file.path(repo_root, "R", "neha_path_utils.R"))
if (!file.exists(file.path(repo_root, "R", "ewce_animal_level_utils.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "ewce_animal_level_utils.R"))

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
  error_message <- tryCatch({
    force(expression)
    NA_character_
  }, error = function(e) conditionMessage(e))
  if (is.na(error_message) || !grepl(pattern, error_message, ignore.case = TRUE)) {
    stop(message, " Observed: ", ifelse(is.na(error_message), "<no error>", error_message), call. = FALSE)
  }
}

animal_ids_by_condition <- list(
  paired_cno = paste0("A", 1:3),
  paired_veh = paste0("A", 4:6),
  unpaired_cno = paste0("A", 7:9),
  unpaired_veh = paste0("A", 10:12)
)
fixture_metadata <- do.call(rbind, lapply(sample_classes, function(sample_class) {
  do.call(rbind, lapply(names(animal_ids_by_condition), function(condition) {
    ids <- animal_ids_by_condition[[condition]]
    data.frame(
      Sample = paste(ids, sample_class, sep = "_"),
      AnimalID = ids,
      condition_code = names(condition_code_map)[match(condition, condition_code_map)],
      condition = condition,
      sample_class = sample_class,
      phenotypeWithinUnit = paste(sample_class, condition, sep = "_"),
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(fixture_metadata) <- NULL

fixture_matrix <- outer(
  seq_len(7L),
  seq_len(nrow(fixture_metadata)),
  function(protein, sample) protein + sample / 100
)
rownames(fixture_matrix) <- paste0("protein_", seq_len(nrow(fixture_matrix)))
colnames(fixture_matrix) <- fixture_metadata$Sample
fixture_column_metadata <- rbind(
  AnimalID = fixture_metadata$AnimalID,
  condition_code = fixture_metadata$condition_code,
  condition = fixture_metadata$condition,
  sample_class = fixture_metadata$sample_class,
  phenotypeWithinUnit = fixture_metadata$phenotypeWithinUnit
)
colnames(fixture_column_metadata) <- fixture_metadata$Sample
fixture_gct <- list(matrix = fixture_matrix, column_metadata = fixture_column_metadata)

validated <- validate_neha_ewce_animal_input(fixture_gct, "fixture_animal_level.gct", expected_n = 3L)
expect_identical(ncol(validated$expression_matrix), 48L, "Validated animal-level EWCE input must have 48 animal x sample-class observations.")
expect_identical(nrow(validated$sample_metadata), 48L, "Animal-level metadata row count is incorrect.")
expect_true(!anyDuplicated(paste(validated$sample_metadata$AnimalID, validated$sample_metadata$Stratum, sep = "\r")),
            "More than one sample entered an AnimalID x sample_class unit.")
expect_true(all(validated$count_audit$n_animals == 3L), "Every sample_class x condition must contain exactly three animals.")
expect_identical(validated$transformations_after_aggregation, "none", "Animal-level input contract introduced a transformation.")

manifest <- neha_primary_contrast_manifest()
expect_identical(nrow(manifest), 12L, "EWCE must use exactly the 12 shared Neha primary contrasts.")
expect_identical(
  as.integer(table(factor(manifest$sample_class, levels = sample_classes))),
  rep(3L, 4L),
  "EWCE must use three contrasts for each sample class."
)

specs <- lapply(sample_classes, function(sample_class) {
  prepare_neha_ewce_limma_stratum(
    validated$expression_matrix,
    validated$sample_metadata,
    sample_class,
    expected_n = 3L
  )
})
names(specs) <- sample_classes
for (sample_class in sample_classes) {
  spec <- specs[[sample_class]]
  expect_identical(spec$n_animal_observations, 12L, "Each sample-class limma fit must contain 12 animal observations.")
  expect_true(all(table(spec$sample_metadata$Cond) == 3L), "A limma condition does not contain exactly three animals.")
  expect_true(!anyDuplicated(spec$sample_metadata$AnimalID), "A hemisphere-level duplicate entered an animal-level limma fit.")
  expect_equal(
    unname(spec$expression_matrix),
    unname(validated$expression_matrix[, spec$sample_metadata$Sample, drop = FALSE]),
    "Preparing limma input normalized, filtered, imputed, or otherwise changed abundance values."
  )
  expected_class_manifest <- manifest[manifest$sample_class == sample_class, , drop = FALSE]
  expect_identical(colnames(spec$contrast_matrix), expected_class_manifest$canonical_comparison,
                   "Limma contrast names do not match the shared canonical comparisons.")
  for (i in seq_len(nrow(expected_class_manifest))) {
    contrast <- spec$contrast_matrix[, i]
    expect_identical(unname(contrast[expected_class_manifest$case_condition[[i]]]), 1,
                     "Canonical contrast numerator does not have coefficient +1.")
    expect_identical(unname(contrast[expected_class_manifest$reference_condition[[i]]]), -1,
                     "Canonical contrast denominator does not have coefficient -1.")
    expect_identical(sum(contrast != 0), 2L, "Canonical contrast contains an unexpected design coefficient.")
  }

  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Focused EWCE tests require the pipeline dependency 'limma'.", call. = FALSE)
  }
  fit <- limma::lmFit(spec$expression_matrix, spec$design)
  fit <- limma::eBayes(
    limma::contrasts.fit(fit, spec$contrast_matrix),
    trend = TRUE,
    robust = TRUE
  )
  first_contrast <- spec$contrast_metadata[1, , drop = FALSE]
  modeled <- limma::topTable(fit, coef = first_contrast$canonical_comparison, number = Inf, sort.by = "none")
  expected_logfc <- rowMeans(spec$expression_matrix[, as.character(spec$sample_metadata$Cond) == first_contrast$numerator_condition, drop = FALSE]) -
    rowMeans(spec$expression_matrix[, as.character(spec$sample_metadata$Cond) == first_contrast$denominator_condition, drop = FALSE])
  expect_equal(modeled$logFC, expected_logfc, "Moderated limma logFC does not preserve numerator-minus-denominator direction.")
  expect_true(all(is.finite(modeled$t)), "Moderated limma t statistics are nonfinite in the valid animal-level fixture.")
}

all_contrasts <- unlist(lapply(specs, function(x) colnames(x$contrast_matrix)), use.names = FALSE)
expect_identical(all_contrasts, manifest$canonical_comparison,
                 "Four animal-level limma fits do not produce the 12 manifest comparisons in order.")

duplicate_fixture <- fixture_gct
duplicate_fixture$column_metadata["AnimalID", 2] <- duplicate_fixture$column_metadata["AnimalID", 1]
expect_error(
  validate_neha_ewce_animal_input(duplicate_fixture, "duplicate.gct"),
  "duplicated within sample_class/condition",
  "Duplicate AnimalID within a sample-class condition must fail closed."
)

hemisphere_fixture <- fixture_gct
hemisphere_fixture$matrix <- cbind(hemisphere_fixture$matrix, hemisphere_fixture$matrix[, 1])
colnames(hemisphere_fixture$matrix)[ncol(hemisphere_fixture$matrix)] <- paste0(colnames(fixture_matrix)[1], "_right")
hemisphere_fixture$column_metadata <- cbind(
  hemisphere_fixture$column_metadata,
  hemisphere_fixture$column_metadata[, 1]
)
colnames(hemisphere_fixture$column_metadata)[ncol(hemisphere_fixture$column_metadata)] <- colnames(hemisphere_fixture$matrix)[ncol(hemisphere_fixture$matrix)]
expect_error(
  validate_neha_ewce_animal_input(hemisphere_fixture, "hemisphere_rows.gct"),
  "duplicated within sample_class/condition",
  "Hemisphere rows must not enter limma as independent animal observations."
)

inconsistent_fixture <- fixture_gct
inconsistent_col <- which(
  inconsistent_fixture$column_metadata["sample_class", ] == "neuron" &
    inconsistent_fixture$column_metadata["condition", ] == "paired_cno" &
    inconsistent_fixture$column_metadata["AnimalID", ] == "A3"
)
inconsistent_fixture$column_metadata["AnimalID", inconsistent_col] <- "AX"
expect_error(
  validate_neha_ewce_animal_input(inconsistent_fixture, "inconsistent.gct"),
  "inconsistent AnimalID membership",
  "Inconsistent animal membership across sample classes must fail closed."
)

missing_animal_fixture <- fixture_gct
missing_animal_fixture$column_metadata["AnimalID", 1] <- ""
expect_error(
  validate_neha_ewce_animal_input(missing_animal_fixture, "missing_animal.gct"),
  "AnimalID is missing",
  "Missing AnimalID must fail closed."
)

underfilled_fixture <- fixture_gct
underfilled_fixture$matrix <- underfilled_fixture$matrix[, -1, drop = FALSE]
underfilled_fixture$column_metadata <- underfilled_fixture$column_metadata[, -1, drop = FALSE]
expect_error(
  validate_neha_ewce_animal_input(underfilled_fixture, "underfilled.gct"),
  "exactly 3 animals per sample_class/condition",
  "An animal-level group with fewer than three animals must fail closed."
)

nonfinite_fixture <- fixture_gct
nonfinite_fixture$matrix[1, 1] <- NA_real_
expect_error(
  validate_neha_ewce_animal_input(nonfinite_fixture, "nonfinite.gct"),
  "nonfinite values",
  "Nonfinite animal-level abundance must fail rather than be filtered or imputed."
)

params <- list(
  seed = 42L,
  reps = 10000L,
  annot_levels = c(1L, 2L),
  primary_annot_level = 2L,
  top_n_values = c(100L, 250L, 500L),
  primary_top_n = 250L,
  fdr_alpha = 0.05,
  marker_top_n = 200L,
  conditions = condition_levels,
  reference_condition = reference_condition
)
audits <- do.call(rbind, lapply(specs, function(spec) {
  make_neha_ewce_differential_audit(
    spec$contrast_metadata,
    source_path = "fixture_animal_level.gct",
    output_path = "fixture_output.xlsx",
    n_proteins = spec$n_proteins,
    analysis_params = params,
    execution_status = "success"
  )
}))
expect_identical(nrow(audits), 12L, "Differential audit must contain one row per canonical comparison.")
expect_true(all(audits$numerator_animal_n == 3L & audits$denominator_animal_n == 3L),
            "Differential audit did not record three numerator and denominator animals.")
expect_true(all(audits$execution_status == "success"), "Differential audit execution status is incorrect.")
expect_true(all(audits$post_aggregation_normalization == "none" &
                  audits$post_aggregation_filtering == "none" &
                  audits$post_aggregation_imputation == "none"),
            "Differential audit does not preserve the no-transform contract.")
required_audit_columns <- c(
  "canonical_comparison", "sample_class", "numerator_condition", "denominator_condition",
  "numerator_animal_n", "denominator_animal_n", "animal_ids_used",
  "source_animal_level_input", "n_proteins_entering_differential_analysis",
  "differential_output_path", "execution_status", "error_message",
  "ewce_seed", "ewce_reps", "ewce_annotation_levels", "ewce_primary_top_n", "ewce_fdr_alpha"
)
expect_true(all(required_audit_columns %in% names(audits)),
            "Differential audit is missing one or more required validation fields.")

expect_error(
  validate_neha_ewce_output_root(
    file.path(tempdir(), "historical", "child"),
    file.path(tempdir(), "historical")
  ),
  "cannot overwrite the historical EWCE root",
  "Animal-level EWCE output must reject the historical result root."
)

# Accept path: an isolated root outside the historical tree must be returned unchanged.
# (Previously only the reject path was covered here, unlike tests/test_pca_animal_level.R.)
ewce_isolated_root <- validate_neha_ewce_output_root(
  file.path(tempdir(), "animal_ewce"),
  file.path(tempdir(), "historical")
)
expect_true(grepl("animal_ewce$", ewce_isolated_root),
            "Animal-level EWCE output must accept and return an isolated output root.")
expect_true(!neha_ewce_path_is_within(ewce_isolated_root, file.path(tempdir(), "historical")),
            "Isolated EWCE output root must not be reported as inside the historical root.")

ewce_script <- paste(readLines(file.path(repo_root, "05_celltype_enrichment_EWCE", "01_EWCE.r"), warn = FALSE), collapse = "\n")
expect_true(grepl("validate_neha_ewce_animal_input\\(animal_gct", ewce_script),
            "EWCE does not consume the validated animal-level input contract.")
expect_true(grepl("limma::lmFit\\(limma_input\\$expression_matrix, limma_input\\$design\\)", ewce_script),
            "Limma is not fitted to the validated animal-level matrix and design.")
expect_true(grepl("trend = TRUE", ewce_script, fixed = TRUE) && grepl("robust = TRUE", ewce_script, fixed = TRUE),
            "Existing moderated eBayes settings changed.")
expect_true(grepl("EWCE::bootstrap_enrichment_test", ewce_script, fixed = TRUE),
            "Existing EWCE bootstrap methodology was removed.")
expect_true(!grepl("duplicateCorrelation", ewce_script, fixed = TRUE),
            "Hemisphere pseudoreplication was replaced with duplicateCorrelation rather than animal-level samples.")
expect_true(!grepl("row_medians|rowSums\\(is.finite|x\\[idx\\]", ewce_script),
            "EWCE still filters or imputes abundance values inside the limma path.")

cat("Animal-level EWCE sampling-unit tests passed.\n")
