#!/usr/bin/env Rscript

# Publication release, stage 04 -- normalised differential-protein results.
#
# Produces
#   differential_analysis/primary_differential_proteins.tsv.gz  12 x 5349 = 64188 rows
#   differential_analysis/primary_differential_summary.tsv      12 rows
#
# NO STATISTIC IS RECOMPUTED. Every numeric value is read from the canonical animal-level
# split tables (02_data/animal_level/split/forward/) and copied through unchanged; the
# stage then re-reads each canonical file and asserts bit-level equality of the exported
# values against it.
#
# The statistical source is split/, not mapped/, on purpose. mapped/ drops the 22 protein
# groups whose identifiers are not mouse (human keratins, porcine trypsin and similar
# contaminants), so exporting from it would make 22 tested proteins per comparison
# disappear from the published record without saying so. Here all 5349 tested rows are
# published and id_mapping_status marks which ones carry no mouse annotation.
#
# Column renaming, and why:
#   split `gene_symbol` -> protein_group_id     it holds UniProt ENTRY NAMES, not symbols
#   split `Description` -> (not published)      it holds gene symbols, not a description;
#                                               the real gene_symbol and protein_description
#                                               come from processed_data/protein_feature_annotation
#   split `log2fc`      -> effect_size_sd_units  (NOT a log2 fold change -- see below)
#                                               standardized abundance difference, SD units
#   split `t`           -> moderated_t
#   split `aveExpr`     -> average_standardized_abundance
#   split `pval`        -> P.Value
#   split `padj`        -> adj.P.Val
#   split `B`           -> B_log_odds

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_validation.R"))
release_require("digest")

DATA_ROOT <- release_data_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/04_build_differential_results.R"

release_banner("stage 04 -- differential protein results")

ANIMAL_ROOT <- file.path(DATA_ROOT, "02_data", "animal_level")
SPLIT_FORWARD <- file.path(ANIMAL_ROOT, "split", "forward")

paths <- list(
  split_index = file.path(ANIMAL_ROOT, "split", "indexComparisons.csv"),
  mapped_index = file.path(ANIMAL_ROOT, "mapped", "indexMappedComparisons.csv"),
  contrast_manifest = release_path("metadata", "primary_contrast_manifest.tsv", create_dir = FALSE),
  feature_annotation = release_path("processed_data", "protein_feature_annotation.tsv.gz",
                                    create_dir = FALSE)
)
for (nm in names(paths)) release_assert_exists(paths[[nm]], nm)

split_index <- release_read_csv(paths$split_index)
mapped_index <- release_read_csv(paths$mapped_index)
contrasts <- release_read_tsv_plain(paths$contrast_manifest)
features <- read.delim(gzfile(paths$feature_annotation), sep = "\t", quote = "",
                       stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")

if (nrow(contrasts) != RELEASE_DESIGN_INVARIANTS$n_primary_contrasts) {
  stop("Contrast manifest does not hold 12 primary contrasts.", call. = FALSE)
}

# Every split table must still hash to what was recorded downstream, otherwise the
# published statistics would not be the ones that were validated. The split index itself
# records only the source GCT hash; the per-table hash lives in the mapped index, which
# consumed these files (`source_split_sha256`).
if (!"source_split_sha256" %in% names(mapped_index)) {
  stop("indexMappedComparisons.csv has no source_split_sha256 column; cannot verify ",
       "the split tables against a recorded hash.", call. = FALSE)
}
split_hashes <- stats::setNames(
  as.character(mapped_index$source_split_sha256[
    match(contrasts$canonical_contrast, mapped_index$canonical_contrast)]),
  contrasts$canonical_contrast)
if (anyNA(split_hashes)) {
  stop("The mapped index does not record a split hash for every primary contrast.",
       call. = FALSE)
}
# The source GCT recorded in the split index must still be the locked statistical GCT.
locked_stat_sha <- RELEASE_LOCKED_ARTEFACTS$protigy_stat_gct$sha256
if (!all(as.character(split_index$source_gct_sha256) == locked_stat_sha)) {
  stop("indexComparisons.csv no longer points every contrast at the locked ProTigy ",
       "statistical GCT.", call. = FALSE)
}

feature_idx_by_id <- stats::setNames(seq_len(nrow(features)), features$protein_group_id)

SPLIT_NUMERIC <- c(log2fc = "effect_size_sd_units", aveExpr = "average_standardized_abundance",
                   t = "moderated_t", pval = "P.Value", padj = "adj.P.Val",
                   B = "B_log_odds", logpval = "neg_log10_P_value",
                   sign.logP = "signed_neg_log10_P_value")

# The abundance matrix these statistics were computed on is standardised separately for
# each protein across the measurement-level dataset: row means are approximately 0 and row
# SDs approximately 1 (measured in stage 03, and corroborated by limma's AveExpr sitting at
# ~0 for every protein). A difference of group means on that scale is a STANDARDIZED
# ABUNDANCE DIFFERENCE expressed in SD units of that scale, not a log2 ratio. The canonical
# column is nevertheless named `log2fc` (ProTigy `logFC`), so publishing it under that name
# would assert a fold-change reading the data do not support. It is renamed here and the
# source name is recorded alongside, so provenance is not lost either.
#
# The wording comes from RELEASE_EFFECT_SIZE in R/release_validation.R rather than being
# written out again here: four builders and the validator all describe this quantity, and a
# local copy is how they drift apart. Note it is NOT called a "standardised mean difference"
# unqualified -- that phrasing invites a Cohen's d reading, i.e. scaling by a pooled
# WITHIN-GROUP SD, which is a different quantity from the across-dataset per-protein SD
# actually used here.
EFFECT_SIZE_DEFINITION <- RELEASE_EFFECT_SIZE$definition
EFFECT_SIZE_SOURCE_COLUMN <- RELEASE_EFFECT_SIZE$source_field_detail
EFFECT_SIZE_UNITS <- RELEASE_EFFECT_SIZE$units

blocks <- vector("list", nrow(contrasts))
verify <- vector("list", nrow(contrasts))

for (i in seq_len(nrow(contrasts))) {
  row <- contrasts[i, , drop = FALSE]
  csv_path <- file.path(SPLIT_FORWARD, paste0(row$canonical_contrast, ".csv"))
  release_assert_exists(csv_path, paste("split table", row$canonical_contrast))

  observed_hash <- release_sha256(csv_path)
  recorded_hash <- unname(split_hashes[[row$canonical_contrast]])
  if (!is.na(recorded_hash) && nzchar(recorded_hash) &&
      !identical(observed_hash, recorded_hash)) {
    stop("Split table hash changed since it was indexed: ", csv_path, call. = FALSE)
  }

  d <- release_read_csv(csv_path)
  if (nrow(d) != RELEASE_DESIGN_INVARIANTS$n_proteins_statistical) {
    stop("Split table ", basename(csv_path), " has ", nrow(d), " rows, expected 5349.",
         call. = FALSE)
  }
  missing_cols <- setdiff(c("gene_symbol", "Description", names(SPLIT_NUMERIC), "significant"),
                          names(d))
  if (length(missing_cols)) {
    stop("Split table ", basename(csv_path), " is missing column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  protein_group_id <- as.character(d$gene_symbol)
  fidx <- feature_idx_by_id[protein_group_id]
  if (anyNA(fidx)) {
    stop("Feature annotation does not cover every protein in ", basename(csv_path),
         call. = FALSE)
  }

  significant <- as.logical(d$significant)
  # Record how `significant` is defined rather than trusting or recomputing it.
  agrees_with_fdr <- identical(significant, as.numeric(d$padj) < 0.05)

  block <- data.frame(
    sample_class = row$sample_class,
    contrast_family = row$contrast_family,
    canonical_comparison = row$canonical_comparison,
    canonical_contrast = row$canonical_contrast,
    historical_comparison_alias = row$historical_comparison_alias,
    numerator_condition = row$numerator_condition,
    denominator_condition = row$denominator_condition,
    n_numerator_animals = as.integer(row$n_numerator_animals),
    n_denominator_animals = as.integer(row$n_denominator_animals),
    analysis_unit = "animal",
    primary_or_secondary = "primary",
    protein_group_id = protein_group_id,
    uniprot_accession = features$uniprot_accession[fidx],
    gene_symbol = features$gene_symbol[fidx],
    protein_description = features$protein_description[fidx],
    id_mapping_status = features$id_mapping_status[fidx],
    source_matrix_gene_symbols = as.character(d$Description),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  for (src in names(SPLIT_NUMERIC)) {
    block[[SPLIT_NUMERIC[[src]]]] <- as.numeric(d[[src]])
  }
  # Short per-row units marker. The full definition would be a ~300-character constant
  # repeated 64,188 times; it lives once in the summary table, the data dictionary and
  # README_DATA.md instead.
  block$effect_size_units <- EFFECT_SIZE_UNITS
  block$effect_size_source_column <- EFFECT_SIZE_SOURCE_COLUMN
  # Bare source token as well as the prose form, so provenance is machine-readable without
  # parsing a sentence. The exported effect value is this field, copied, with no transform.
  block$source_statistic_field <- RELEASE_EFFECT_SIZE$source_field
  block$significant_fdr_0_05 <- significant

  # Mapped-only counts are carried alongside the full-set counts because the manuscript
  # figures were drawn from mapped/ (5327 rows). Without both numbers a reader cannot
  # reconcile, for example, the panel-3D figure count of 1132 with the full tested set.
  is_mapped <- features$id_mapping_status[fidx] == "mapped"

  blocks[[i]] <- block
  verify[[i]] <- list(contrast = row$canonical_contrast, path = csv_path,
                      sha256 = observed_hash, rows = nrow(d),
                      significant_is_padj_lt_0_05 = agrees_with_fdr,
                      n_significant = sum(significant, na.rm = TRUE),
                      n_up = sum(significant & as.numeric(d$log2fc) > 0, na.rm = TRUE),
                      n_down = sum(significant & as.numeric(d$log2fc) < 0, na.rm = TRUE),
                      n_mapped = sum(is_mapped),
                      n_significant_mapped = sum(significant & is_mapped, na.rm = TRUE),
                      n_up_mapped = sum(significant & is_mapped &
                                          as.numeric(d$log2fc) > 0, na.rm = TRUE),
                      n_down_mapped = sum(significant & is_mapped &
                                            as.numeric(d$log2fc) < 0, na.rm = TRUE))
}

differential <- do.call(rbind, blocks)
rownames(differential) <- NULL

expected_rows <- RELEASE_DESIGN_INVARIANTS$n_primary_contrasts *
  RELEASE_DESIGN_INVARIANTS$n_proteins_statistical
if (nrow(differential) != expected_rows) {
  stop("Differential export has ", nrow(differential), " rows, expected ", expected_rows,
       ".", call. = FALSE)
}
if (length(unique(differential$canonical_comparison)) != 12L) {
  stop("Differential export does not cover 12 complete comparison blocks.", call. = FALSE)
}
if (!all(table(differential$canonical_comparison) ==
           RELEASE_DESIGN_INVARIANTS$n_proteins_statistical)) {
  stop("Comparison blocks are not all 5349 rows.", call. = FALSE)
}
if (release_column_is_misleading_gene_symbol(differential$gene_symbol)) {
  stop("Published gene_symbol column contains UniProt identifiers; refusing to write.",
       call. = FALSE)
}
release_log("  assembled ", nrow(differential), " rows across 12 comparisons")

sig_definition_uniform <- all(vapply(verify, function(v) isTRUE(v$significant_is_padj_lt_0_05),
                                     logical(1)))
if (!sig_definition_uniform) {
  release_log("  NOTE: canonical `significant` is not exactly adj.P.Val < 0.05 in every ",
              "table; the canonical flag is published as-is and the definition is recorded.")
}

# --------------------------------------------------------------------------------------
# exact re-verification against the canonical sources
# --------------------------------------------------------------------------------------
# Re-read every split table from disk and require identical doubles. This is what makes
# "we did not recompute anything" checkable rather than asserted.

release_log("  re-verifying exported statistics against the canonical split tables ...")
for (v in verify) {
  d <- release_read_csv(v$path)
  ex <- differential[differential$canonical_contrast == v$contrast, , drop = FALSE]
  if (!identical(as.character(ex$protein_group_id), as.character(d$gene_symbol))) {
    stop("Protein order diverged from the canonical table for ", v$contrast, call. = FALSE)
  }
  for (src in names(SPLIT_NUMERIC)) {
    if (!identical(ex[[SPLIT_NUMERIC[[src]]]], as.numeric(d[[src]]))) {
      stop("Exported ", SPLIT_NUMERIC[[src]], " differs from canonical `", src, "` in ",
           v$contrast, call. = FALSE)
    }
  }
  if (!identical(ex$significant_fdr_0_05, as.logical(d$significant))) {
    stop("Exported significance flag differs from canonical `significant` in ", v$contrast,
         call. = FALSE)
  }
}
release_log("  all 12 comparisons verified identical to their canonical sources")

# --------------------------------------------------------------------------------------
# per-comparison summary (counts of already-flagged rows; no thresholding is applied)
# --------------------------------------------------------------------------------------

summary_tbl <- do.call(rbind, lapply(seq_along(verify), function(i) {
  v <- verify[[i]]
  row <- contrasts[contrasts$canonical_contrast == v$contrast, , drop = FALSE]
  data.frame(
    sample_class = row$sample_class,
    contrast_family = row$contrast_family,
    canonical_comparison = row$canonical_comparison,
    numerator_condition = row$numerator_condition,
    denominator_condition = row$denominator_condition,
    n_numerator_animals = as.integer(row$n_numerator_animals),
    n_denominator_animals = as.integer(row$n_denominator_animals),
    analysis_unit = "animal",
    n_proteins_tested = v$rows,
    n_significant_fdr_0_05 = v$n_significant,
    n_significant_higher_in_numerator = v$n_up,
    n_significant_higher_in_denominator = v$n_down,
    n_proteins_tested_mapped_only = v$n_mapped,
    n_significant_fdr_0_05_mapped_only = v$n_significant_mapped,
    n_significant_higher_in_numerator_mapped_only = v$n_up_mapped,
    n_significant_higher_in_denominator_mapped_only = v$n_down_mapped,
    significance_definition = ifelse(isTRUE(v$significant_is_padj_lt_0_05),
                                     "adj.P.Val < 0.05 (Benjamini-Hochberg)",
                                     "canonical `significant` flag as computed upstream"),
    effect_size_units = EFFECT_SIZE_UNITS,
    effect_size_definition = EFFECT_SIZE_DEFINITION,
    effect_size_source_column = EFFECT_SIZE_SOURCE_COLUMN,
    source_statistic_field = RELEASE_EFFECT_SIZE$source_field,
    canonical_source_path = v$path,
    canonical_source_sha256 = v$sha256,
    stringsAsFactors = FALSE, check.names = FALSE)
}))
summary_tbl <- summary_tbl[order(match(summary_tbl$sample_class, sample_classes),
                                  summary_tbl$contrast_family), , drop = FALSE]
rownames(summary_tbl) <- NULL

for (i in seq_len(nrow(summary_tbl))) {
  release_log("    ", format(summary_tbl$canonical_comparison[i], width = 50), " ",
              summary_tbl$n_significant_fdr_0_05[i], " / ",
              summary_tbl$n_proteins_tested[i], " FDR<0.05")
}

# --------------------------------------------------------------------------------------
# write
# --------------------------------------------------------------------------------------

source_paths <- vapply(verify, function(v) v$path, character(1))
source_hashes <- vapply(verify, function(v) v$sha256, character(1))

w1 <- release_write_table(differential,
                          release_path("differential_analysis",
                                       "primary_differential_proteins.tsv.gz"))
release_register("differential_analysis/primary_differential_proteins.tsv.gz",
                 paste("animal-level differential protein statistics, long format,",
                       "one row per protein x primary comparison"),
                 source_paths, source_hashes, STAGE, "tsv.gz")

w2 <- release_write_table(summary_tbl,
                          release_path("differential_analysis",
                                       "primary_differential_summary.tsv"))
release_register("differential_analysis/primary_differential_summary.tsv",
                 "per-comparison counts of FDR-significant proteins",
                 source_paths, source_hashes, STAGE, "tsv")

release_log("  wrote primary_differential_proteins.tsv.gz (", w1$rows, "x", w1$cols, ")")
release_log("  wrote primary_differential_summary.tsv (", w2$rows, "x", w2$cols, ")")
release_log("stage 04 complete")
