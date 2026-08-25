#!/usr/bin/env Rscript

# Compare protein-level differential statistics for the established learning
# and paired-CNO contrasts. The animal-level mapped branch contains contrast
# statistics, not per-animal expression, so this script deliberately produces
# contrast-statistic concordance outputs only.

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1]]) else "04_differential_expression_enrichment/04_compare_sig_expr.r"
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "R", "analysis_labels.R"))
source(file.path(project_root, "R", "animal_level_enrichment_utils.R"))

config <- resolve_neha_enrichment_config()
enrichment_index_path <- Sys.getenv(
  "NEHA_ENRICHMENT_INDEX",
  unset = file.path(config$output_root, "indexEnrichmentComparisons.csv")
)
enrichment_index_path <- neha_enrichment_normalize_path(enrichment_index_path, must_work = TRUE)
if (!neha_enrichment_path_is_within(enrichment_index_path, config$output_root)) {
  stop("compare_sig_expr requires the canonical animal-level enrichment index.", call. = FALSE)
}
index <- utils::read.csv(enrichment_index_path, stringsAsFactors = FALSE, check.names = FALSE)
validate_neha_enrichment_index(index, require_files = TRUE)
if (any(index$execution_status != "success")) stop("compare_sig_expr requires 12 successful enrichment comparisons.", call. = FALSE)

output_root <- file.path(config$output_root, "compare_sig_expr")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_index_path <- file.path(output_root, "indexCompareSigExpr.csv")
if (file.exists(output_index_path) && !isTRUE(config$force)) {
  stop("compare_sig_expr output already exists; set NEHA_ENRICHMENT_FORCE=true to replace isolated canonical outputs.", call. = FALSE)
}

find_contrast <- function(sample_class, numerator, denominator) {
  rows <- which(
    index$sample_class == sample_class &
      index$numerator_condition == numerator &
      index$denominator_condition == denominator
  )
  if (length(rows) != 1L) stop("Expected exactly one indexed contrast for ", sample_class, ": ", numerator, " vs ", denominator, call. = FALSE)
  index[rows, , drop = FALSE]
}

read_collapsed <- function(row) {
  mapped <- read_neha_enrichment_mapped_file(
    row$mapped_input_path,
    expected_rows = row$n_source_protein_rows,
    expected_sha256 = row$mapped_input_sha256
  )
  collapse_neha_enrichment_accessions(mapped, "log2fc")$collapsed
}

safe_z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  finite <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(finite) >= 2L && stats::sd(x[finite]) > 0) out[finite] <- as.numeric(scale(x[finite]))
  out
}

expected <- neha_primary_contrast_manifest()
records <- lapply(unique(expected$sample_class), function(sample_class) {
  learning_meta <- find_contrast(sample_class, "paired_veh", "unpaired_veh")
  cno_meta <- find_contrast(sample_class, "paired_cno", "paired_veh")
  learning <- read_collapsed(learning_meta)
  cno <- read_collapsed(cno_meta)
  learning <- learning[, c(
    "uniprot_accession", "original_protein_id", "Description", "mapped_gene_symbol",
    "source_row_id", "log2fc", "t", "pval", "padj", "significant"
  ), drop = FALSE]
  cno <- cno[, c(
    "uniprot_accession", "original_protein_id", "Description", "mapped_gene_symbol",
    "source_row_id", "log2fc", "t", "pval", "padj", "significant"
  ), drop = FALSE]
  names(learning)[-1] <- paste0(names(learning)[-1], "_learning")
  names(cno)[-1] <- paste0(names(cno)[-1], "_cno")
  joined <- merge(learning, cno, by = "uniprot_accession", all = TRUE, sort = TRUE)
  joined$sample_class <- sample_class
  joined$learning_canonical_comparison <- learning_meta$canonical_comparison
  joined$cno_canonical_comparison <- cno_meta$canonical_comparison
  joined$learning_historical_alias <- learning_meta$historical_comparison_alias
  joined$cno_historical_alias <- cno_meta$historical_comparison_alias
  joined$z_log2fc_learning <- safe_z(joined$log2fc_learning)
  joined$z_log2fc_cno <- safe_z(joined$log2fc_cno)
  joined$fdr_significant_learning <- is.finite(joined$padj_learning) & joined$padj_learning < config$fdr_threshold
  joined$fdr_significant_cno <- is.finite(joined$padj_cno) & joined$padj_cno < config$fdr_threshold
  joined$direction_relation <- ifelse(
    !is.finite(joined$log2fc_learning) | !is.finite(joined$log2fc_cno), "not_shared",
    ifelse(sign(joined$log2fc_learning) == sign(joined$log2fc_cno), "concordant", "opposed")
  )
  joined <- joined[, c(
    "sample_class", "uniprot_accession",
    "original_protein_id_learning", "Description_learning", "mapped_gene_symbol_learning", "source_row_id_learning",
    "log2fc_learning", "t_learning", "pval_learning", "padj_learning", "significant_learning", "z_log2fc_learning",
    "original_protein_id_cno", "Description_cno", "mapped_gene_symbol_cno", "source_row_id_cno",
    "log2fc_cno", "t_cno", "pval_cno", "padj_cno", "significant_cno", "z_log2fc_cno",
    "fdr_significant_learning", "fdr_significant_cno", "direction_relation",
    "learning_canonical_comparison", "cno_canonical_comparison",
    "learning_historical_alias", "cno_historical_alias"
  ), drop = FALSE]
  all_path <- write_neha_enrichment_csv(joined, file.path(output_root, paste0(sample_class, "_learning_vs_paired_cno_protein_statistics.csv")))
  signature <- joined[joined$fdr_significant_learning, , drop = FALSE]
  signature <- signature[order(signature$padj_learning, -abs(signature$log2fc_learning), signature$uniprot_accession, method = "radix"), , drop = FALSE]
  signature_path <- write_neha_enrichment_csv(signature, file.path(output_root, paste0(sample_class, "_FDR_learning_signature.csv")))

  shared <- is.finite(joined$log2fc_learning) & is.finite(joined$log2fc_cno)
  pearson <- if (sum(shared) >= 3L) stats::cor(joined$log2fc_learning[shared], joined$log2fc_cno[shared], method = "pearson") else NA_real_
  spearman <- if (sum(shared) >= 3L) stats::cor(joined$log2fc_learning[shared], joined$log2fc_cno[shared], method = "spearman") else NA_real_
  signature_shared <- joined$fdr_significant_learning & shared
  signature_pearson <- if (sum(signature_shared) >= 3L) stats::cor(
    joined$log2fc_learning[signature_shared], joined$log2fc_cno[signature_shared], method = "pearson"
  ) else NA_real_

  plot_path <- NA_character_
  if (isTRUE(config$plots_enabled) && any(shared) && requireNamespace("ggplot2", quietly = TRUE)) {
    plot_data <- joined[shared, , drop = FALSE]
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(log2fc_learning, log2fc_cno, colour = fdr_significant_learning)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
      ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
      ggplot2::geom_point(alpha = 0.55, size = 1.25) +
      ggplot2::scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#B2182B")) +
      ggplot2::labs(
        title = paste(sample_class, "protein-level contrast concordance"),
        subtitle = "Positive values indicate higher abundance in each canonical numerator",
        x = "paired VEH vs unpaired VEH log2FC",
        y = "paired CNO vs paired VEH log2FC",
        colour = "Learning FDR < threshold"
      ) + ggplot2::theme_minimal(base_size = 10)
    plot_path <- file.path(output_root, paste0(sample_class, "_learning_vs_paired_cno_log2fc_scatter.png"))
    ggplot2::ggsave(plot_path, plot = plot, width = 8, height = 7, dpi = 300)
    plot_path <- neha_enrichment_normalize_path(plot_path, must_work = TRUE)
  }

  data.frame(
    sample_class = sample_class,
    learning_canonical_comparison = learning_meta$canonical_comparison,
    cno_canonical_comparison = cno_meta$canonical_comparison,
    learning_historical_alias = learning_meta$historical_comparison_alias,
    cno_historical_alias = cno_meta$historical_comparison_alias,
    learning_mapped_input_path = learning_meta$mapped_input_path,
    learning_mapped_input_sha256 = learning_meta$mapped_input_sha256,
    cno_mapped_input_path = cno_meta$mapped_input_path,
    cno_mapped_input_sha256 = cno_meta$mapped_input_sha256,
    shared_unique_uniprot_count = sum(shared),
    learning_fdr_significant_count = nrow(signature),
    cno_fdr_significant_count = sum(joined$fdr_significant_cno),
    concordant_count = sum(joined$direction_relation == "concordant"),
    opposed_count = sum(joined$direction_relation == "opposed"),
    pearson_log2fc = pearson,
    spearman_log2fc = spearman,
    learning_signature_pearson_log2fc = signature_pearson,
    full_statistics_path = all_path,
    learning_signature_path = signature_path,
    scatter_path = plot_path,
    analysis_scope = "contrast_statistics_only;per_animal_expression_not_available_in_mapped_manifest",
    enrichment_index_path = enrichment_index_path,
    run_timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
})

comparison_index <- do.call(rbind, records)
write_neha_enrichment_csv(comparison_index, output_index_path)
message("Canonical compare_sig_expr completed: ", output_index_path)
