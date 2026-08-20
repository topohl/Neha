args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("tests", "test_animal_level_proteomics.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "animal_level_proteomics_utils.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "animal_level_proteomics_utils.R"))
source(file.path(repo_root, "R", "protigy_input_utils.R"))

expect_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

expect_identical <- function(actual, expected, message) {
  if (!identical(actual, expected)) {
    stop(
      message, "\nExpected: ", paste(expected, collapse = ", "),
      "\nActual: ", paste(actual, collapse = ", "),
      call. = FALSE
    )
  }
}

expect_error <- function(expr, pattern, message) {
  error <- tryCatch({
    force(expr)
    NULL
  }, error = identity)
  if (is.null(error) || !grepl(pattern, conditionMessage(error), perl = TRUE, ignore.case = TRUE)) {
    stop(message, call. = FALSE)
  }
}

make_metadata <- function() {
  data.frame(
    sample_id = c("run_A1_L_cfos_x", "run_A1_R_cfos_x", "run_A2_L_cfos_x", "run_A2_R_cfos_x", "run_A1_L_neuron_x", "run_A1_R_neuron_x"),
    AnimalID = c("A1", "A1", "A2", "A2", "A1", "A1"),
    ReplicateGroup = c("Left", "Right", "L", "R", "Left", "Right"),
    sample_class = c("cfos", "cfos", "cfos", "cfos", "neuron", "neuron"),
    condition_code = c("1", "1", "2", "2", "1", "1"),
    condition = c("paired_cno", "paired_cno", "paired_veh", "paired_veh", "paired_cno", "paired_cno"),
    exclude = FALSE,
    stringsAsFactors = FALSE
  )
}

metadata <- make_metadata()
expected_units <- unique(metadata[c("AnimalID", "sample_class", "condition_code", "condition")])
plan <- prepare_animal_level_aggregation(metadata, "test", expected_units = expected_units)
values <- matrix(
  c(1, 3, 10, 14, 100, 104, 2, 6, 20, 24, 200, 208),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("P1", "P2"), metadata$sample_id)
)
aggregated <- aggregate_animal_level_matrix(values, rownames(values), plan)

# 1. Exact arithmetic L/R mean.
expect_identical(unname(aggregated$primary["P1", "A1_cfos"]), 2, "Left=1 and Right=3 must aggregate to 2.")

# 2. Different animals remain separate.
expect_identical(unname(aggregated$primary["P1", "A2_cfos"]), 12, "Animals were incorrectly averaged together.")

# 3. Different sample classes remain separate.
expect_identical(unname(aggregated$primary["P1", "A1_neuron"]), 102, "Sample classes were incorrectly averaged together.")

# 4 and 11. An animal cannot span conditions, preventing condition averaging.
multi_condition <- metadata
multi_condition$condition_code[2] <- "2"
multi_condition$condition[2] <- "paired_veh"
expect_error(
  prepare_animal_level_aggregation(multi_condition, "test"),
  "multiple conditions",
  "AnimalID mapping to multiple conditions must fail."
)

# 5. Left-only observations are retained without imputation.
left_only <- metadata[metadata$sample_id != "run_A1_R_cfos_x", , drop = FALSE]
left_plan <- prepare_animal_level_aggregation(left_only, "test", expected_units = expected_units)
left_values <- values[, left_only$sample_id, drop = FALSE]
left_agg <- aggregate_animal_level_matrix(left_values, rownames(values), left_plan)
expect_identical(unname(left_agg$primary["P1", "A1_cfos"]), 1, "Left-only observation was not retained exactly.")
expect_identical(
  left_plan$audit$aggregation_method[left_plan$audit$output_column_name == "A1_cfos"],
  "single_observed_hemisphere_no_imputation",
  "Left-only method was not recorded."
)

# 6. Right-only observations are retained without imputation.
right_only <- metadata[metadata$sample_id != "run_A1_L_cfos_x", , drop = FALSE]
right_plan <- prepare_animal_level_aggregation(right_only, "test", expected_units = expected_units)
right_values <- values[, right_only$sample_id, drop = FALSE]
right_agg <- aggregate_animal_level_matrix(right_values, rownames(values), right_plan)
expect_identical(unname(right_agg$primary["P1", "A1_cfos"]), 3, "Right-only observation was not retained exactly.")

# 7. Missing-both units are represented explicitly and omitted from matrices.
missing_expected <- rbind(
  expected_units,
  data.frame(AnimalID = "A2", sample_class = "neuron", condition_code = "2", condition = "paired_veh")
)
missing_plan <- prepare_animal_level_aggregation(metadata, "test", expected_units = missing_expected)
missing_row <- missing_plan$audit[missing_plan$audit$AnimalID == "A2" & missing_plan$audit$sample_class == "neuron", ]
expect_identical(missing_row$hemisphere_status, "missing_both", "Missing-both unit was not explicit in the audit.")
missing_agg <- aggregate_animal_level_matrix(values, rownames(values), missing_plan)
expect_true(!"A2_neuron" %in% colnames(missing_agg$primary), "Missing-both unit must not be fabricated.")

# 8. Duplicate Left fails.
duplicate_left <- rbind(metadata, transform(metadata[1, ], sample_id = "run_A1_L_cfos_duplicate"))
expect_error(prepare_animal_level_aggregation(duplicate_left, "test"), "duplicate Left", "Duplicate Left must fail.")

# 9. Duplicate Right fails.
duplicate_right <- rbind(metadata, transform(metadata[2, ], sample_id = "run_A1_R_cfos_duplicate"))
expect_error(prepare_animal_level_aggregation(duplicate_right, "test"), "duplicate Right", "Duplicate Right must fail.")

# 10. Metadata/sample-name hemisphere disagreement fails.
disagreement <- metadata
disagreement$ReplicateGroup[1] <- "Right"
expect_error(prepare_animal_level_aggregation(disagreement, "test"), "disagrees", "Hemisphere disagreement must fail.")

# 12. A repeated source ID cannot enter two output units.
duplicate_source <- metadata
duplicate_source$sample_id[5] <- duplicate_source$sample_id[1]
expect_error(prepare_animal_level_aggregation(duplicate_source, "test"), "must be unique", "A source sample cannot enter two outputs.")

# 13. Output column IDs must be unique.
colliding <- data.frame(
  sample_id = c("run_A1_L_class-a_x", "run_A1_R_class-a_x", "run_A1_L_class_a_x", "run_A1_R_class_a_x"),
  AnimalID = "A1",
  ReplicateGroup = rep(c("Left", "Right"), 2),
  sample_class = rep(c("class-a", "class_a"), each = 2),
  condition_code = "1",
  condition = "paired_cno",
  exclude = FALSE,
  stringsAsFactors = FALSE
)
expect_error(prepare_animal_level_aggregation(colliding, "test"), "output column IDs", "Output ID collision must fail.")

# 14. Ordering is deterministic regardless of metadata row order.
set.seed(17)
shuffled <- metadata[sample(seq_len(nrow(metadata))), , drop = FALSE]
shuffled_plan <- prepare_animal_level_aggregation(shuffled, "test", expected_units = expected_units)
expect_identical(shuffled_plan$audit$output_column_name, plan$audit$output_column_name, "Output ordering is not deterministic.")

# 15. Protein identity and order are unchanged.
expect_identical(rownames(aggregated$primary), rownames(values), "Protein IDs/order changed during aggregation.")

# 16. Every numeric value equals the expected L/R mean.
expected_matrix <- cbind(
  A1_cfos = c(P1 = 2, P2 = 4),
  A2_cfos = c(P1 = 12, P2 = 22),
  A1_neuron = c(P1 = 102, P2 = 204)
)
expected_matrix <- expected_matrix[, colnames(aggregated$primary), drop = FALSE]
expect_true(
  isTRUE(all.equal(unname(aggregated$primary), unname(expected_matrix), tolerance = 0, check.attributes = FALSE)),
  "Animal-level numeric values differ from exact expected means."
)

# 17-20. Practical ProTigy GCT v1.3 structure and parser alignment.
gct_path <- tempfile(fileext = ".gct")
descriptions <- c("Protein one", "Protein two")
write_protigy_gct_v13(aggregated$primary, aggregated$primary_audit, descriptions, gct_path)
parsed <- validate_protigy_gct_v13(gct_path, expected_matrix = aggregated$primary)
expect_identical(unname(parsed$dimensions), c(2L, 3L, 1L, 5L), "ProTigy GCT dimensions are incorrect.")
expect_identical(names(parsed$row_descriptors), "Description", "Explicit Description descriptor is missing.")
physical_lines <- readLines(gct_path, warn = FALSE)
metadata_fields <- strsplit(physical_lines[4:8], "\t", fixed = TRUE)
expect_true(all(vapply(metadata_fields, `[[`, character(1), 2) == "na"), "Column metadata filler must be na.")
expect_true(
  identical(rownames(parsed$matrix), rownames(parsed$row_descriptors)),
  "Parser matrix and row descriptor names are misaligned."
)

unlink(gct_path)
message("Animal-level proteomics and ProTigy GCT tests passed (20 focused contracts).")
