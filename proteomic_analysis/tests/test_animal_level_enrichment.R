#!/usr/bin/env Rscript

repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "animal_level_enrichment_utils.R"))

assert_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)
assert_equal <- function(left, right, message) if (!identical(left, right)) stop(message, call. = FALSE)
assert_error <- function(expression, pattern, message) {
  observed <- tryCatch({ force(expression); NULL }, error = conditionMessage)
  if (is.null(observed) || !grepl(pattern, observed, ignore.case = TRUE)) {
    stop(message, if (!is.null(observed)) paste0(" Observed: ", observed) else " No error was raised.", call. = FALSE)
  }
}

manifest <- neha_primary_contrast_manifest()
assert_equal(nrow(manifest), 12L, "Shared Neha manifest must contain exactly 12 comparisons.")
assert_equal(length(unique(manifest$sample_class)), 4L, "Shared Neha manifest must contain four sample classes.")
assert_true(all(table(manifest$sample_class) == 3L), "Each Neha sample class must contain exactly three contrasts.")

fixture_root <- tempfile("neha_enrichment_fixture_")
dir.create(fixture_root, recursive = TRUE)
on.exit(unlink(fixture_root, recursive = TRUE, force = TRUE), add = TRUE)
mapped_root <- file.path(fixture_root, "canonical_mapped")
forward_root <- file.path(mapped_root, "forward")
dir.create(forward_root, recursive = TRUE)

make_mapped <- function(multiplier = 1) data.frame(
  gene_symbol = c("P11111", "O54931", "O54931", "Q22222"),
  original_protein_id = c("ONE_MOUSE", "AKAP2_MOUSE", "PALM2_MOUSE", "TWO_MOUSE"),
  Description = c("One", "Akap2", "Palm2", "Two"),
  mapped_gene_symbol = c("One", "Akap2", "Palm2", "Two"),
  uniprot_accession = c("P11111", "O54931", "O54931", "Q22222"),
  mapping_status = "mapped",
  mapping_strategy = "fixture",
  source_row_id = 1:4,
  multi_protein = FALSE,
  log2fc = c(1.5, 0.5, -2, -0.25) * multiplier,
  aveExpr = c(10, 11, 12, 13),
  t = c(3, 1, -4, -0.5),
  pval = c(0.001, 0.2, 0.002, 0.8),
  padj = c(0.01, 0.3, 0.02, 0.9),
  B = c(2, -1, 3, -4),
  significant = c(TRUE, FALSE, TRUE, FALSE),
  logpval = -log10(c(0.001, 0.2, 0.002, 0.8)),
  sign.logP = c(3, 0.7, -2.7, -0.1),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

paths <- file.path(forward_root, paste0(manifest$canonical_contrast, ".csv"))
for (i in seq_along(paths)) utils::write.csv(make_mapped(i), paths[[i]], row.names = FALSE)
hashes <- vapply(paths, neha_enrichment_sha256, character(1))
index <- data.frame(
  canonical_comparison = manifest$canonical_comparison,
  canonical_contrast = manifest$canonical_contrast,
  sample_class = manifest$sample_class,
  numerator_condition = manifest$case_condition,
  denominator_condition = manifest$reference_condition,
  historical_comparison_alias = manifest$historical_comparison_name,
  mapping_direction = "forward",
  output_mapped_path = paths,
  n_input_proteins = 4L,
  n_output_mapped_rows = 4L,
  n_mapped = 4L,
  n_unmapped = 0L,
  mapped_output_sha256 = hashes,
  source_split_index_path = file.path(fixture_root, "split_index.csv"),
  source_split_index_sha256 = rep("split-index-hash", 12L),
  source_split_sha256 = paste0("split-hash-", seq_len(12L)),
  source_gct_sha256 = rep("gct-hash", 12L),
  mapping_reference_path = file.path(fixture_root, "mapping.csv"),
  mapping_reference_version = "fixture-v1",
  mapping_reference_snapshot_date_utc = "2026-08-24",
  mapping_reference_sha256 = rep("mapping-hash", 12L),
  row_accounting_valid = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
index_path <- file.path(mapped_root, "indexMappedComparisons.csv")
utils::write.csv(index, index_path, row.names = FALSE)

read_index <- read_neha_enrichment_mapped_index(index_path, mapped_root, verify_hashes = TRUE)
assert_equal(nrow(read_index), 12L, "Canonical indexed mapped-path reader lost comparisons.")
assert_equal(read_index$mapped_input_path, unname(vapply(paths, neha_enrichment_normalize_path, character(1), must_work = TRUE)), "Reader did not preserve exact indexed mapped paths.")
assert_equal(read_index$source_gct_sha256, rep("gct-hash", 12L), "Mapped provenance was not propagated.")
assert_equal(read_index$historical_comparison_alias, manifest$historical_comparison_name, "Aliases must be retained as metadata.")

mapped <- read_neha_enrichment_mapped_file(read_index$mapped_input_path[[1]], expected_rows = 4L, expected_sha256 = hashes[[1]])
assert_true(is.logical(mapped$significant), "Mapped significant flags must remain logical.")
collapsed <- collapse_neha_enrichment_accessions(mapped, "log2fc")
assert_equal(collapsed$n_duplicate_rows_collapsed, 1L, "Exactly one duplicate UniProt row should be collapsed.")
assert_equal(collapsed$n_duplicated_accessions, 1L, "Exactly one duplicated UniProt accession should be audited.")
selected <- collapsed$collapsed[collapsed$collapsed$uniprot_accession == "O54931", , drop = FALSE]
assert_equal(selected$original_protein_id, "PALM2_MOUSE", "Largest absolute log2FC representative was not selected.")
assert_equal(selected$log2fc, -2, "Duplicate collapse did not preserve the selected statistic sign.")
assert_true(all(c("AKAP2_MOUSE", "PALM2_MOUSE") %in% collapsed$duplicate_audit$original_protein_id), "Duplicate audit omitted contributing protein IDs.")
assert_equal(sum(collapsed$duplicate_audit$selected_representative), 1L, "Duplicate audit must mark one representative.")

rank_one <- build_neha_gsea_rank(collapsed$collapsed, "log2fc")$rank
reordered <- mapped[c(4, 2, 1, 3), , drop = FALSE]
rank_two <- build_neha_gsea_rank(collapse_neha_enrichment_accessions(reordered, "log2fc")$collapsed, "log2fc")$rank
assert_equal(rank_one, rank_two, "GSEA rank construction must be deterministic under source-row reordering.")
assert_true(!anyDuplicated(names(rank_one)), "GSEA identifiers must be unique after audited collapse.")
assert_true(all(is.finite(rank_one)), "GSEA ranks must not contain NA or infinite values.")
assert_equal(names(rank_one)[[1]], "P11111", "Deterministic GSEA sorting or positive numerator direction is incorrect.")
assert_equal(unname(rank_one[[length(rank_one)]]), -2, "Negative numerator-direction rank was inverted.")

ora <- build_neha_ora_sets(collapsed$collapsed, fdr_threshold = 0.05, top_abs_log2fc = 1)
assert_equal(ora$universe, sort(unique(mapped$uniprot_accession)), "ORA universe must be the measured successfully mapped proteome.")
assert_equal(length(ora$universe), 3L, "Duplicate UniProt accessions inflated the ORA universe.")
assert_equal(length(ora$all_significant), 2L, "FDR-significant ORA list is incorrect.")
assert_true(all(ora$all_significant %in% ora$universe), "ORA foreground escaped the measured universe.")

bad_index <- index
bad_index$output_mapped_path[[2]] <- bad_index$output_mapped_path[[1]]
utils::write.csv(bad_index, file.path(mapped_root, "duplicateIndex.csv"), row.names = FALSE)
assert_error(
  read_neha_enrichment_mapped_index(file.path(mapped_root, "duplicateIndex.csv"), mapped_root, verify_hashes = FALSE),
  "duplicates", "Duplicate indexed mapped paths must be rejected."
)

defaults <- neha_enrichment_default_paths()
assert_error(
  validate_neha_enrichment_mapped_root(defaults$historical_input_roots[[2]], defaults$historical_input_roots),
  "historical", "Canonical mode must reject the historical mapped root."
)
assert_error(
  validate_neha_enrichment_mapped_root(dirname(defaults$historical_input_roots[[2]]), defaults$historical_input_roots),
  "historical", "Canonical mode must reject a mapped root that contains a historical mapped root."
)
assert_error(
  validate_neha_enrichment_output_root(
    defaults$historical_output_roots[[2]], mapped_root, defaults$split_root,
    defaults$historical_input_roots, defaults$historical_output_roots
  ),
  "historical|overlap", "Canonical mode must protect historical output roots."
)
assert_error(
  validate_neha_enrichment_output_root(
    file.path(mapped_root, "enrichment"), mapped_root, defaults$split_root,
    defaults$historical_input_roots, defaults$historical_output_roots
  ),
  "overlap", "Canonical output must not overlap mapped inputs."
)

set.seed(91)
rng_before <- .Random.seed
seeded_one <- with_neha_enrichment_seed(1234L, stats::runif(5))
assert_equal(.Random.seed, rng_before, "Deterministic enrichment seed helper changed caller RNG state.")
seeded_two <- with_neha_enrichment_seed(1234L, stats::runif(5))
assert_equal(seeded_one, seeded_two, "Deterministic enrichment seed helper is not reproducible.")

cat("All animal-level enrichment tests passed.\n")
