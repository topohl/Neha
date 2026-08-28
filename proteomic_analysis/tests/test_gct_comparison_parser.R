args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("tests", "test_gct_comparison_parser.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "protigy_stat_gct_utils.R"))

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

# Historical comparison-name parsing used by existing demo and legacy consumers.
cases <- c(
  "neuron_3.over.neuron_2" = "neuron_unpaired_cno_vs_neuron_paired_veh",
  "neuron_4.over.neuron_2" = "neuron_unpaired_veh_vs_neuron_paired_veh",
  "neuron_4.over.neuron_3" = "neuron_unpaired_veh_vs_neuron_unpaired_cno",
  "bg_2.over.bg_1" = "neuropil_paired_veh_vs_neuropil_paired_cno",
  "bg_3.over.bg_1" = "neuropil_unpaired_cno_vs_neuropil_paired_cno",
  "bg_4.over.bg_1" = "neuropil_unpaired_veh_vs_neuropil_paired_cno",
  "cfos_1.over.bg_1" = "cfos_paired_cno_vs_neuropil_paired_cno",
  "cfos_2.over.bg_1" = "cfos_paired_veh_vs_neuropil_paired_cno",
  "CA1_so_bg_2.over.CA1_so_bg_1" = "CA1_so_neuropil_paired_veh_vs_CA1_so_neuropil_paired_cno",
  "CA1_so_neuropil_1.over.CA1_so_neuropil_2" = "CA1_so_neuropil_paired_cno_vs_CA1_so_neuropil_paired_veh",
  "mcherry_1.over.mcherry_2" = "mcherry_paired_cno_vs_mcherry_paired_veh"
)
for (key in names(cases)) {
  parsed <- parse_comparison_key(key)
  if (is.null(parsed)) stop("Expected parser success for: ", key, call. = FALSE)
  expect_identical(parsed$name, unname(cases[[key]]), paste("Unexpected parsed name for", key))
}

invalid_cases <- c(
  "neuron.over.neuron_2",
  "neuron_5.over.neuron_2",
  "foo_2.over.neuron_2",
  "cfos_1_vs_bg_1"
)
for (key in invalid_cases) {
  if (!is.null(parse_comparison_key(key))) stop("Expected parser failure for: ", key, call. = FALSE)
}

parsed_codes <- parse_condition_code(c("1", "sample_4", NA, "none"))
expect_identical(parsed_codes, c("1", "4", NA_character_, NA_character_), "parse_condition_code must return one value per input.")
parsed_classes <- parse_sample_class(c("sample_bg_1", "demo_cFos_2", NA, "none"))
expect_identical(parsed_classes, c("neuropil", "cfos", NA_character_, NA_character_), "parse_sample_class must return one value per input.")

# Both ProTigy comparison separators canonicalize to _over_.
stat_fields <- parse_protigy_stat_field(c(
  "logFC.mcherry_paired_cno_over_mcherry_paired_veh",
  "adj.P.Val.mcherry_1.over.mcherry_2"
))
expect_identical(stat_fields$metric, c("logFC", "adj.P.Val"), "Metric prefixes were not parsed correctly.")
expect_identical(
  stat_fields$naming_style,
  c("corrected_underscore_over", "historical_dot_over"),
  "Comparison naming styles were not detected."
)
expect_identical(
  stat_fields$comparison,
  c("mcherry_paired_cno_over_mcherry_paired_veh", "mcherry_1_over_mcherry_2"),
  "Comparison separators were not canonicalized."
)
parsed_styles <- parse_protigy_comparison(stat_fields$comparison)
expect_identical(
  parsed_styles$canonical_comparison,
  rep("mcherry_paired_cno_over_mcherry_paired_veh", 2),
  "Corrected and historical labels did not resolve to one biological comparison."
)

# The shared contrast manifest is exactly the requested 4 classes x 3 contrasts.
manifest <- primary_contrast_manifest()
expect_identical(nrow(manifest), 12L, "The primary manifest must contain exactly 12 contrasts.")
expect_identical(unique(manifest$sample_class), sample_classes, "Primary contrast sample-class order changed.")
validated_manifest <- validate_primary_comparison_contract(manifest$canonical_comparison)
expect_identical(nrow(validated_manifest), 12L, "Exactly 12 primary contrasts were not accepted.")
expect_true(all(validated_manifest$accepted_primary), "A required primary contrast was rejected.")

expect_error(
  validate_protigy_comparisons("cfos_paired_cno_over_neuropil_paired_veh", strict_primary = TRUE),
  "cross_sample_class",
  "Cross-sample-class comparisons must be rejected in strict animal-level mode."
)
expect_error(
  validate_protigy_comparisons("mcherry_paired_cno_over_mcherry_unpaired_cno", strict_primary = TRUE),
  "unsupported_condition_contrast",
  "Unsupported within-class condition comparisons must be rejected."
)

make_stat_gct <- function(path, comparison, duplicate_logfc = FALSE) {
  metrics <- protigy_required_da_metrics()
  fields <- paste0(metrics, ".", comparison)
  values <- c("2.5", "5", "3.5", "0.01", "0.02", "1.2", "TRUE", "2")
  if (isTRUE(duplicate_logfc)) {
    fields <- c(fields, paste0("logFC.", comparison))
    values <- c(values, "2.5")
  }
  matrix_field <- paste0("sign.logP.", comparison)
  header <- c("id", "Description", fields, matrix_field)
  nrhd <- 1L + length(fields)
  metadata <- c("Batch", rep("na", nrhd), "synthetic")
  stopifnot(length(metadata) == length(header))
  lines <- c(
    "#1.3",
    paste(c(2L, 1L, nrhd, 1L), collapse = "\t"),
    paste(header, collapse = "\t"),
    paste(metadata, collapse = "\t"),
    paste(c("A0A000001_MOUSE", "GeneA", values, "4"), collapse = "\t"),
    paste(c("Q00001_MOUSE", "GeneB", replace(values, 7L, "FALSE"), "-4"), collapse = "\t")
  )
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

# Statistical fields in the row-descriptor area, arbitrary descriptor count,
# and a nonzero column-metadata count are parsed from declared dimensions.
synthetic <- tempfile(fileext = ".gct")
make_stat_gct(synthetic, "mcherry_paired_cno_over_mcherry_paired_veh")
gct <- read_protigy_stat_gct(synthetic, strict_primary = FALSE)
expect_identical(unname(gct$dimensions), c(2L, 1L, 9L, 1L), "Synthetic statistical-GCT dimensions were parsed incorrectly.")
expect_identical(gct$n_protein_rows_read, 2L, "Declared protein rows were not preserved.")
expect_true(all(gct$parsed_fields$physical_role[gct$parsed_fields$metric != "sign.logP"] == "row_descriptor"), "Statistical row descriptors were not recognized as physical fields.")
expect_identical(gct$parsed_fields$physical_role[gct$parsed_fields$metric == "sign.logP"], "matrix_column", "Matrix-area statistic was not detected.")

forward <- extract_protigy_comparison_table(gct, "mcherry_paired_cno_over_mcherry_paired_veh")
expect_identical(forward$gene_symbol, c("A0A000001_MOUSE", "Q00001_MOUSE"), "Original GCT protein IDs were not retained.")
expect_identical(forward$Description, c("GeneA", "GeneB"), "Description was not retained and aligned.")
expect_true(is.logical(forward$significant), "significant must remain logical, not numeric.")
expect_identical(forward$significant, c(TRUE, FALSE), "TRUE/FALSE significance values changed.")
expect_identical(forward$log2fc, c(2.5, 2.5), "logFC extraction failed.")
expect_identical(forward$pval, c(0.01, 0.01), "P.Value extraction failed.")
expect_identical(forward$padj, c(0.02, 0.02), "adj.P.Val extraction failed.")
expect_identical(forward$aveExpr, c(5, 5), "AveExpr extraction failed.")
expect_identical(forward$B, c(1.2, 1.2), "B extraction failed.")
expect_identical(forward$logpval, c(2, 2), "Log.P.Value extraction failed.")
expect_identical(forward$sign.logP, c(4, -4), "Known additional ProTigy statistic was not retained.")

reverse <- reverse_protigy_metric_frame(forward, attr(forward, "metric_by_column"))
expect_identical(reverse$log2fc, -forward$log2fc, "Reverse logFC sign is incorrect.")
expect_identical(reverse$t, -forward$t, "Reverse t sign is incorrect.")
for (column in c("pval", "padj", "B", "significant", "aveExpr", "logpval", "sign.logP")) {
  expect_identical(reverse[[column]], forward[[column]], paste("Reverse must not change", column))
}

# A fully historical .over. statistical GCT resolves to the same canonical result.
historical <- tempfile(fileext = ".gct")
make_stat_gct(historical, "mcherry_1.over.mcherry_2")
historical_gct <- read_protigy_stat_gct(historical, strict_primary = FALSE)
expect_identical(
  unique(historical_gct$parsed_fields$canonical_comparison),
  "mcherry_paired_cno_over_mcherry_paired_veh",
  "Historical .over. fields did not canonicalize to the canonical comparison."
)

# Duplicate physical fields are detected after canonicalization and cannot be extracted.
duplicate <- tempfile(fileext = ".gct")
make_stat_gct(duplicate, "mcherry_paired_cno_over_mcherry_paired_veh", duplicate_logfc = TRUE)
duplicate_gct <- read_protigy_stat_gct(duplicate, strict_primary = FALSE)
expect_true(duplicate_gct$duplicate_metric_comparison_fields > 0L, "Duplicate metric/comparison fields were not detected.")
expect_error(
  extract_protigy_comparison_table(duplicate_gct, "mcherry_paired_cno_over_mcherry_paired_veh"),
  "Duplicate metric/comparison",
  "Duplicate metric/comparison fields must prevent extraction."
)

unlink(c(synthetic, historical, duplicate))
message("GCT comparison and statistical-result parser tests passed.")
