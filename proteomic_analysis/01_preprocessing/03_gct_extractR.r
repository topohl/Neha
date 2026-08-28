# ================================================================
# animal-level ProTigy statistical-result GCT extraction
#
# Statistical fields can occur anywhere among the physical GCT fields,
# including the row-descriptor area. This stage extracts existing ProTigy
# results only; it does not fit or rerun any statistical model.
# ================================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  sub("^--file=", "", file_arg)
} else {
  file.path("01_preprocessing", "03_gct_extractR.r")
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "protigy_stat_gct_utils.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "protigy_stat_gct_utils.R"))

option_or_env <- function(option_name, env_name, default) {
  value <- getOption(option_name)
  if (!is.null(value) && nzchar(trimws(as.character(value)))) return(as.character(value))
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(trimws(value))) return(value)
  default
}

input_gct <- option_or_env(
  "proteomics.protigy_stat_gct_input",
  "PROTEOMICS_PROTIGY_STAT_GCT_INPUT",
  "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/02_data/animal_level/stat_results_for_ssGSEA_neha_proteome.gct"
)
output_root <- option_or_env(
  "proteomics.protigy_stat_gct_output_root",
  "PROTEOMICS_PROTIGY_STAT_GCT_OUTPUT_ROOT",
  "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/02_data/animal_level/split"
)
input_gct <- normalizePath(input_gct, winslash = "/", mustWork = TRUE)
output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)

path_is_within <- function(path, parent) {
  path <- tolower(normalizePath(path, winslash = "/", mustWork = FALSE))
  parent <- sub("/+$", "", tolower(normalizePath(parent, winslash = "/", mustWork = FALSE)))
  identical(path, parent) || startsWith(path, paste0(parent, "/"))
}

# Historical hemisphere-level trees were consolidated under 99_historical/ during the
# 2026-08-26 restructure (see CANONICAL_OUTPUTS.md). The former Datasets/raw and
# Datasets/mapped are now 99_historical/datasets_raw and 99_historical/datasets_mapped.
historical_root <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/99_historical"
protected_roots <- c(historical_root, file.path(historical_root, c("datasets_raw", "datasets_mapped")))
if (any(vapply(protected_roots, function(path) path_is_within(output_root, path), logical(1)))) {
  stop(
    "Refusing to write animal-level ProTigy splits under the historical 99_historical tree: ",
    output_root,
    call. = FALSE
  )
}

gct <- read_protigy_stat_gct(input_gct, strict_primary = TRUE)
contract <- validate_stat_gct_contract(gct, expected_n_proteins = 5349L)
expected <- primary_contrast_manifest()

forward_dir <- file.path(output_root, "forward")
reverse_dir <- file.path(output_root, "reverse")
expected_forward_names <- paste0(expected$canonical_contrast, ".csv")
expected_reverse_names <- paste0(
  expected$reference_phenotype,
  "_vs_",
  expected$case_phenotype,
  ".csv"
)

unexpected_existing_csv <- function(directory, expected_names) {
  if (!dir.exists(directory)) return(character())
  existing <- list.files(directory, pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
  setdiff(existing, expected_names)
}
unexpected <- c(
  unexpected_existing_csv(forward_dir, expected_forward_names),
  unexpected_existing_csv(reverse_dir, expected_reverse_names)
)
if (length(unexpected)) {
  stop(
    "Output split directories contain unexpected CSV files; refusing a mixed contract: ",
    paste(unique(unexpected), collapse = ", "),
    call. = FALSE
  )
}

dir.create(forward_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reverse_dir, recursive = TRUE, showWarnings = FALSE)

index_rows <- lapply(seq_len(nrow(expected)), function(i) {
  comparison <- expected$canonical_comparison[[i]]
  forward <- extract_protigy_comparison_table(gct, comparison)
  metric_by_column <- attr(forward, "metric_by_column")
  reverse <- reverse_protigy_metric_frame(forward, metric_by_column)
  forward_path <- normalizePath(
    file.path(forward_dir, expected_forward_names[[i]]),
    winslash = "/",
    mustWork = FALSE
  )
  reverse_path <- normalizePath(
    file.path(reverse_dir, expected_reverse_names[[i]]),
    winslash = "/",
    mustWork = FALSE
  )

  utils::write.csv(forward, forward_path, row.names = FALSE, quote = TRUE, na = "")
  utils::write.csv(reverse, reverse_path, row.names = FALSE, quote = TRUE, na = "")

  forward_check <- utils::read.csv(forward_path, stringsAsFactors = FALSE, check.names = FALSE)
  reverse_check <- utils::read.csv(reverse_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (
    nrow(forward_check) != gct$n_protein_rows_read ||
      nrow(reverse_check) != gct$n_protein_rows_read ||
      !identical(as.character(forward_check$gene_symbol), gct$ids) ||
      !identical(as.character(reverse_check$gene_symbol), gct$ids) ||
      !identical(as.character(forward_check$Description), gct$description) ||
      !identical(as.character(reverse_check$Description), gct$description)
  ) {
    stop("Post-write identity or row-count validation failed for: ", comparison, call. = FALSE)
  }

  fields <- gct$parsed_fields[
    gct$parsed_fields$canonical_comparison == comparison,
    ,
    drop = FALSE
  ]
  data.frame(
    canonical_comparison = comparison,
    canonical_contrast = expected$canonical_contrast[[i]],
    sample_class = expected$sample_class[[i]],
    numerator_condition = expected$case_condition[[i]],
    denominator_condition = expected$reference_condition[[i]],
    numerator_phenotype = expected$case_phenotype[[i]],
    denominator_phenotype = expected$reference_phenotype[[i]],
    historical_comparison_alias = expected$historical_comparison_name[[i]],
    forward_filename = basename(forward_path),
    forward_path = forward_path,
    reverse_filename = basename(reverse_path),
    reverse_path = reverse_path,
    statistical_fields_detected = paste(fields$metric, collapse = ";"),
    source_fields = paste(fields$field, collapse = ";"),
    n_proteins = gct$n_protein_rows_read,
    source_gct_path = gct$path,
    source_gct_sha256 = gct$sha256,
    stringsAsFactors = FALSE
  )
})
index <- do.call(rbind, index_rows)
index_path <- file.path(output_root, "indexComparisons.csv")
utils::write.csv(index, index_path, row.names = FALSE, quote = TRUE, na = "")

forward_files <- sort(list.files(forward_dir, pattern = "\\.csv$", full.names = TRUE), method = "radix")
reverse_files <- sort(list.files(reverse_dir, pattern = "\\.csv$", full.names = TRUE), method = "radix")
if (length(forward_files) != 12L || length(reverse_files) != 12L || nrow(index) != 12L) {
  stop("Post-write split contract failed: expected 12 forward, 12 reverse, and 12 index rows.", call. = FALSE)
}

contract_manifest <- data.frame(
  contract_version = "neha_animal_level_protigy_stat_split_v1",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  status = contract$status,
  source_gct_path = gct$path,
  source_gct_sha256 = gct$sha256,
  declared_n_proteins = unname(gct$dimensions[["nrmat"]]),
  declared_n_matrix_columns = unname(gct$dimensions[["ncmat"]]),
  declared_n_row_descriptors = unname(gct$dimensions[["nrhd"]]),
  declared_n_column_metadata_rows = unname(gct$dimensions[["nchd"]]),
  observed_physical_fields = gct$observed_physical_fields,
  comparison_naming_style = gct$naming_style,
  required_da_metrics = paste(protigy_required_da_metrics(), collapse = ";"),
  detected_statistics = paste(gct$metrics_found, collapse = ";"),
  n_primary_comparisons = nrow(index),
  n_forward_files = length(forward_files),
  n_reverse_files = length(reverse_files),
  n_proteins_per_file = gct$n_protein_rows_read,
  protein_ids_unique = gct$protein_ids_unique,
  description_retained = gct$description_present && gct$description_aligned,
  duplicate_metric_comparison_fields = gct$duplicate_metric_comparison_fields,
  output_root = output_root,
  index_path = normalizePath(index_path, winslash = "/", mustWork = TRUE),
  index_sha256 = protigy_file_sha256(index_path),
  stringsAsFactors = FALSE
)
contract_path <- file.path(output_root, "run_contract_manifest.csv")
utils::write.csv(contract_manifest, contract_path, row.names = FALSE, quote = TRUE, na = "")

message("animal-level ProTigy statistical-result extraction completed successfully.")
message("Source GCT SHA-256: ", gct$sha256)
message("Forward files: ", length(forward_files), " (", gct$n_protein_rows_read, " proteins each)")
message("Reverse files: ", length(reverse_files), " (", gct$n_protein_rows_read, " proteins each)")
message("Index: ", normalizePath(index_path, winslash = "/", mustWork = TRUE))
message("Contract manifest: ", normalizePath(contract_path, winslash = "/", mustWork = TRUE))
