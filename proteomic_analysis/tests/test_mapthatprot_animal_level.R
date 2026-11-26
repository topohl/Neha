args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  sub("^--file=", "", file_arg)
} else {
  file.path("tests", "test_mapthatprot_animal_level.R")
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
source(file.path(repo_root, "R", "project_path_utils.R"))
if (!file.exists(file.path(repo_root, "R", "mapthatprot_animal_level_utils.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "mapthatprot_animal_level_utils.R"))

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

make_split_fixture <- function(root) {
  manifest <- primary_contrast_manifest()
  forward_dir <- file.path(root, "forward")
  reverse_dir <- file.path(root, "reverse")
  dir.create(forward_dir, recursive = TRUE)
  dir.create(reverse_dir, recursive = TRUE)

  filenames <- paste0(manifest$canonical_contrast, ".csv")
  make_input <- function(log2fc) {
    data.frame(
      gene_symbol = c("A0A000001_MOUSE", "Q00001_MOUSE"),
      Description = c("GeneA", "GeneB"),
      log2fc = rep(log2fc, 2L),
      aveExpr = c(5, 6),
      t = c(2, -2),
      pval = c(0.01, 0.02),
      padj = c(0.03, 0.04),
      B = c(1, 0),
      significant = c(TRUE, FALSE),
      logpval = c(2, 1.7),
      sign.logP = c(2, -1.7),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  for (filename in filenames) {
    utils::write.csv(make_input(1), file.path(forward_dir, filename), row.names = FALSE)
    utils::write.csv(make_input(-1), file.path(reverse_dir, filename), row.names = FALSE)
  }

  index <- data.frame(
    canonical_comparison = manifest$canonical_comparison,
    canonical_contrast = manifest$canonical_contrast,
    sample_class = manifest$sample_class,
    numerator_condition = manifest$case_condition,
    denominator_condition = manifest$reference_condition,
    historical_comparison_alias = manifest$historical_comparison_name,
    forward_filename = filenames,
    forward_path = normalizePath(file.path(forward_dir, filenames), winslash = "/", mustWork = TRUE),
    reverse_filename = filenames,
    reverse_path = normalizePath(file.path(reverse_dir, filenames), winslash = "/", mustWork = TRUE),
    statistical_fields_detected = "logFC;AveExpr;t;P.Value;adj.P.Val;B;significant;Log.P.Value;sign.logP",
    n_proteins = 2L,
    source_gct_path = "synthetic.gct",
    source_gct_sha256 = paste(rep("a", 64L), collapse = ""),
    stringsAsFactors = FALSE
  )
  utils::write.csv(index, file.path(root, "indexComparisons.csv"), row.names = FALSE)
  invisible(list(index = index, manifest = manifest))
}

fixture_root <- tempfile("mapthatprot_split_")
dir.create(fixture_root)
fixture <- make_split_fixture(fixture_root)

# The index gives 12 forward inputs in manifest order.
indexed <- read_mapthatprot_input_index(fixture_root, direction = "forward", expected_n_proteins = 2L)
expect_identical(nrow(indexed), 12L, "The mapper did not discover exactly 12 indexed comparisons.")
expect_identical(
  indexed$canonical_comparison,
  fixture$manifest$canonical_comparison,
  "Indexed comparisons were not ordered by the shared contrast manifest."
)
expect_true(all(indexed$mapping_direction == "forward"), "Forward selection was not recorded.")
expect_true(all(dirname(indexed$input_path) == normalizePath(file.path(fixture_root, "forward"), winslash = "/")), "Forward discovery escaped split/forward.")
forward_input <- validate_mapthatprot_input_file(indexed$input_path[[1]], expected_n_proteins = 2L)
expect_true(all(forward_input$log2fc == 1), "The indexed forward file was not selected.")
expect_true(is.logical(forward_input$significant), "significant must remain logical at the mapping handoff.")

# Missing indexed files and duplicate comparisons fail before mapping begins.
missing_root <- tempfile("mapthatprot_missing_")
dir.create(missing_root)
missing_fixture <- make_split_fixture(missing_root)
unlink(missing_fixture$index$forward_path[[1]])
expect_error(
  read_mapthatprot_input_index(missing_root, expected_n_proteins = 2L),
  "missing",
  "A missing indexed comparison file must be rejected."
)

duplicate_root <- tempfile("mapthatprot_duplicate_")
dir.create(duplicate_root)
duplicate_fixture <- make_split_fixture(duplicate_root)
duplicate_index_path <- file.path(duplicate_root, "indexComparisons.csv")
duplicate_index <- utils::read.csv(duplicate_index_path, stringsAsFactors = FALSE, check.names = FALSE)
duplicate_index$canonical_comparison[[2]] <- duplicate_index$canonical_comparison[[1]]
utils::write.csv(duplicate_index, duplicate_index_path, row.names = FALSE)
expect_error(
  read_mapthatprot_input_index(duplicate_root, expected_n_proteins = 2L),
  "duplicate",
  "Duplicate indexed comparisons must be rejected."
)

# Mapped and unmapped outputs form a lossless source-row partition and retain
# the original protein ID, Description, and every DA statistic.
input <- data.frame(
  gene_symbol = c("A0A000001_MOUSE", "Q00001_MOUSE", "BAD_ID"),
  Description = c("GeneA", "GeneB", "GeneC"),
  log2fc = c(1, -1, 0.5),
  aveExpr = c(5, 6, 7),
  t = c(2, -2, 1),
  pval = c(0.01, 0.02, 0.5),
  padj = c(0.03, 0.04, 0.7),
  B = c(1, 0, -1),
  significant = c(TRUE, FALSE, FALSE),
  logpval = c(2, 1.7, 0.3),
  sign.logP = c(2, -1.7, 0.3),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
mapping_audit <- data.frame(
  source_row_id = 1:3,
  original_protein_id = input$gene_symbol,
  uniprot_accession = c("A0A000001", "Q00001", NA),
  mapped_gene_symbol = c("GeneA", "GeneB", NA),
  mapping_strategy = c("accept_accession_base", "entry_local_mouse", "ineligible_non_mouse_identifier"),
  mapping_status = c("mapped", "mapped", "unmapped"),
  multi_protein = c(FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
outputs <- build_mapthatprot_output_tables(input, mapping_audit)
expect_identical(nrow(outputs$mapped), 2L, "Mapped output row count is incorrect.")
expect_identical(nrow(outputs$unmapped), 1L, "Unmapped output row count is incorrect.")
expect_identical(outputs$mapped$gene_symbol, c("A0A000001", "Q00001"), "Mapped gene_symbol must remain the UniProt compatibility key.")
expect_identical(outputs$mapped$original_protein_id, input$gene_symbol[1:2], "Original protein IDs were not retained.")
expect_identical(outputs$mapped$Description, input$Description[1:2], "Descriptions were not retained.")
expect_identical(outputs$unmapped$original_protein_id, input$gene_symbol[[3]], "Unmapped original protein ID was lost.")
expect_identical(outputs$unmapped$Description, input$Description[[3]], "Unmapped Description was lost.")
expect_true(all(names(input)[-c(1, 2)] %in% names(outputs$mapped)), "Mapped output dropped DA statistics.")
expect_true(all(names(input)[-c(1, 2)] %in% names(outputs$unmapped)), "Unmapped output dropped DA statistics.")
expect_true(validate_mapthatprot_partition(input, outputs$mapped, outputs$unmapped), "Valid row partition was rejected.")

bad_audit <- mapping_audit
bad_audit$source_row_id[[3]] <- 2L
expect_error(
  build_mapthatprot_output_tables(input, bad_audit),
  "exactly once",
  "Duplicate source rows in the mapping audit must be rejected."
)
bad_unmapped <- outputs$unmapped
bad_unmapped$source_row_id <- outputs$mapped$source_row_id[[1]]
expect_error(
  validate_mapthatprot_partition(input, outputs$mapped, bad_unmapped),
  "row loss or duplication",
  "Mapped/unmapped row loss or duplication must be rejected."
)

# Canonical output cannot target the split inputs or any historical output tree.
defaults <- mapthatprot_default_paths()
for (protected_root in c(defaults$split_root, defaults$historical_roots)) {
  expect_error(
    validate_mapthatprot_output_root(protected_root, defaults$split_root, defaults$historical_roots),
    "protected or historical",
    paste("Exact protected output root was accepted:", protected_root)
  )
  expect_error(
    validate_mapthatprot_output_root(file.path(protected_root, "animal_level_test"), defaults$split_root, defaults$historical_roots),
    "protected or historical",
    paste("Protected output root was accepted:", protected_root)
  )
}
allowed_output <- validate_mapthatprot_output_root(
  file.path(tempdir(), "isolated_mapping_output"),
  defaults$split_root,
  defaults$historical_roots
)
expect_true(grepl("isolated_mapping_output$", allowed_output), "An isolated output root was unexpectedly rejected.")

# Defaults point at the per-compartment run of 2026-11-24, not the pooled split/ and mapped/.
animal_root <- file.path(defaults$project_root, "02_data", "animal_level")
expect_identical(
  defaults$split_root,
  file.path(animal_root, "within_compartment_20261124", "split"),
  "The default split root is not the per-compartment split."
)
expect_identical(
  defaults$output_root,
  file.path(animal_root, "within_compartment_20261124", "mapped"),
  "The default output root is not the per-compartment mapped folder."
)

# Environment overrides and the overwrite refusal, on temp folders. The MapThatProt variables
# are cleared first and restored afterwards.
mapthatprot_env <- c(
  "PROTEOMICS_MAPTHATPROT_SPLIT_ROOT", "PROTEOMICS_MAPTHATPROT_OUTPUT_ROOT",
  "PROTEOMICS_MAPTHATPROT_DIRECTION", "PROTEOMICS_MAPTHATPROT_FORCE"
)
saved_env <- Sys.getenv(mapthatprot_env, unset = NA)
Sys.unsetenv(mapthatprot_env)
output_fixture <- tempfile("mapthatprot_output_")
tryCatch({
  config <- resolve_mapthatprot_config(output_root = output_fixture)
  expect_identical(
    config$split_root,
    normalize_mapthatprot_path(defaults$split_root),
    "Without overrides the default split root was not used."
  )

  Sys.setenv(
    PROTEOMICS_MAPTHATPROT_SPLIT_ROOT = fixture_root,
    PROTEOMICS_MAPTHATPROT_OUTPUT_ROOT = output_fixture
  )
  config <- resolve_mapthatprot_config()
  expect_identical(config$split_root, normalize_mapthatprot_path(fixture_root), "PROTEOMICS_MAPTHATPROT_SPLIT_ROOT was not used.")
  expect_identical(config$output_root, normalize_mapthatprot_path(output_fixture), "PROTEOMICS_MAPTHATPROT_OUTPUT_ROOT was not used.")

  # An output folder that already holds the mapped index is refused and left as it was.
  dir.create(output_fixture)
  forward_index <- file.path(output_fixture, "indexMappedComparisons.csv")
  writeLines("canonical_comparison", forward_index)
  files_before <- list.files(output_fixture, recursive = TRUE, all.files = TRUE, include.dirs = TRUE)
  for (value in c("false", "no", "0")) {
    Sys.setenv(PROTEOMICS_MAPTHATPROT_FORCE = value)
    expect_error(
      resolve_mapthatprot_config(),
      "set PROTEOMICS_MAPTHATPROT_FORCE=true",
      paste("An existing mapped index was not refused with FORCE =", value)
    )
  }
  Sys.unsetenv("PROTEOMICS_MAPTHATPROT_FORCE")
  expect_error(resolve_mapthatprot_config(), basename(output_fixture), "The refusal does not name the output folder.")
  expect_identical(
    list.files(output_fixture, recursive = TRUE, all.files = TRUE, include.dirs = TRUE),
    files_before,
    "The refusal changed the output folder."
  )
  expect_identical(readLines(forward_index), "canonical_comparison", "The refusal changed the existing index.")

  for (value in c("true", "TRUE", "1", "yes")) {
    Sys.setenv(PROTEOMICS_MAPTHATPROT_FORCE = value)
    config <- resolve_mapthatprot_config()
    expect_identical(config$output_root, normalize_mapthatprot_path(output_fixture), paste("FORCE =", value, "did not allow the run."))
  }
  Sys.unsetenv("PROTEOMICS_MAPTHATPROT_FORCE")

  # The check is per direction: a forward index does not block a reverse run; a reverse index does.
  Sys.setenv(PROTEOMICS_MAPTHATPROT_DIRECTION = "reverse")
  config <- resolve_mapthatprot_config()
  expect_identical(config$direction, "reverse", "A reverse run was blocked by the forward index.")
  writeLines("canonical_comparison", file.path(output_fixture, "indexMappedComparisons_reverse.csv"))
  expect_error(
    resolve_mapthatprot_config(),
    "indexMappedComparisons_reverse\\.csv",
    "An existing reverse mapped index was not refused."
  )
  Sys.unsetenv("PROTEOMICS_MAPTHATPROT_DIRECTION")

  # The batch script resolves its configuration (and so the refusal) before it creates or
  # writes any file.
  batch_script <- file.path(repo_root, "02_id_mapping", "01_MapThatProt_batch.r")
  parse_data <- utils::getParseData(parse(batch_script, keep.source = TRUE))
  calls <- parse_data[parse_data$token == "SYMBOL_FUNCTION_CALL", c("line1", "text")]
  resolve_line <- min(calls$line1[calls$text == "resolve_mapthatprot_config"])
  write_calls <- c("dir.create", "write_csv", "write.csv", "writeLines", "saveWorkbook", "file.copy", "sink")
  first_write_line <- min(calls$line1[calls$text %in% write_calls])
  expect_true(
    is.finite(resolve_line) && resolve_line < first_write_line,
    "01_MapThatProt_batch.r writes output before resolve_mapthatprot_config()."
  )

  # Read-only check on the shared drive: the canonical per-compartment mapped index exists,
  # so a run with the defaults is refused.
  Sys.unsetenv(mapthatprot_env)
  default_index <- file.path(defaults$output_root, "indexMappedComparisons.csv")
  if (file.exists(default_index)) {
    expect_error(
      resolve_mapthatprot_config(),
      "set PROTEOMICS_MAPTHATPROT_FORCE=true",
      "A run with the defaults did not refuse the existing per-compartment mapped index."
    )
  } else {
    message("Shared drive not reachable; skipped the default mapped index check.")
  }
}, finally = {
  for (name in mapthatprot_env) {
    if (is.na(saved_env[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, stats::setNames(list(saved_env[[name]]), name))
    }
  }
  unlink(output_fixture, recursive = TRUE)
})

unlink(c(fixture_root, missing_root, duplicate_root), recursive = TRUE)
message("Animal-level MapThatProt handoff tests passed (indexed discovery, identity, row accounting, output protection, default paths and overwrite refusal).")
