#!/usr/bin/env Rscript

# Compare canonical GO GSEA results without rerunning enrichment.

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1]]) else "04_differential_expression_enrichment/02_compareGO.r"
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "R", "analysis_labels.R"))
source(file.path(project_root, "R", "project_path_utils.R"))
source(file.path(project_root, "R", "animal_level_enrichment_utils.R"))

config <- resolve_enrichment_config()
index_path <- Sys.getenv(
  "PROTEOMICS_ENRICHMENT_INDEX",
  unset = file.path(config$output_root, "indexEnrichmentComparisons.csv")
)
index_path <- enrichment_normalize_path(index_path, must_work = TRUE)
if (!enrichment_path_is_within(index_path, config$output_root)) {
  stop("Canonical enrichment index must resolve inside PROTEOMICS_ENRICHMENT_OUTPUT_ROOT.", call. = FALSE)
}
index <- utils::read.csv(index_path, stringsAsFactors = FALSE, check.names = FALSE)
validate_enrichment_index(index, require_files = TRUE)
if (any(index$execution_status != "success")) stop("compareGO requires 12 successful enrichment comparisons.", call. = FALSE)

output_root <- file.path(config$output_root, "compareGO")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_index_path <- file.path(output_root, "indexCompareGO.csv")
if (file.exists(output_index_path) && !isTRUE(config$force)) {
  stop("compareGO output already exists; set PROTEOMICS_ENRICHMENT_FORCE=true to replace isolated canonical outputs.", call. = FALSE)
}

required_go_columns <- c("ID", "Description", "NES", "pvalue", "p.adjust", "qvalue", "core_enrichment")
tables <- lapply(seq_len(nrow(index)), function(i) {
  source_path <- enrichment_normalize_path(index$go_gsea_output[[i]], must_work = TRUE)
  if (!enrichment_path_is_within(source_path, file.path(config$output_root, "per_comparison"))) {
    stop("GO input resolves outside canonical per_comparison outputs: ", source_path, call. = FALSE)
  }
  x <- utils::read.csv(source_path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(required_go_columns, names(x))
  if (length(missing)) stop("GO GSEA output is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!nrow(x)) return(data.frame(
    canonical_comparison = character(), canonical_contrast = character(), sample_class = character(),
    numerator_condition = character(), denominator_condition = character(), historical_comparison_alias = character(),
    ID = character(), Description = character(), NES = numeric(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), core_enrichment = character(), direction = character(), stringsAsFactors = FALSE
  ))
  x$canonical_comparison <- index$canonical_comparison[[i]]
  x$canonical_contrast <- index$canonical_contrast[[i]]
  x$sample_class <- index$sample_class[[i]]
  x$numerator_condition <- index$numerator_condition[[i]]
  x$denominator_condition <- index$denominator_condition[[i]]
  x$historical_comparison_alias <- index$historical_comparison_alias[[i]]
  x$direction <- ifelse(x$NES > 0, "numerator", ifelse(x$NES < 0, "denominator", "zero"))
  x[, c(
    "canonical_comparison", "canonical_contrast", "sample_class", "numerator_condition",
    "denominator_condition", "historical_comparison_alias", "ID", "Description", "NES",
    "pvalue", "p.adjust", "qvalue", "core_enrichment", "direction"
  ), drop = FALSE]
})
combined <- do.call(rbind, tables)
rownames(combined) <- NULL
if (anyDuplicated(paste(combined$canonical_comparison, combined$ID, sep = "\r"))) {
  stop("A GO term occurs more than once within a canonical comparison.", call. = FALSE)
}
combined <- combined[order(combined$canonical_comparison, combined$p.adjust, -abs(combined$NES), combined$ID, method = "radix"), , drop = FALSE]
fdr_threshold <- config$fdr_threshold
significant <- combined[is.finite(combined$p.adjust) & combined$p.adjust < fdr_threshold, , drop = FALSE]

combined_path <- write_enrichment_csv(combined, file.path(output_root, "GO_GSEA_all_comparisons.csv"))
significant_path <- write_enrichment_csv(significant, file.path(output_root, "GO_GSEA_FDR_significant.csv"))

summary_rows <- lapply(seq_len(nrow(index)), function(i) {
  comparison <- index$canonical_comparison[[i]]
  all_i <- combined[combined$canonical_comparison == comparison, , drop = FALSE]
  sig_i <- significant[significant$canonical_comparison == comparison, , drop = FALSE]
  data.frame(
    canonical_comparison = comparison,
    canonical_contrast = index$canonical_contrast[[i]],
    sample_class = index$sample_class[[i]],
    numerator_condition = index$numerator_condition[[i]],
    denominator_condition = index$denominator_condition[[i]],
    historical_comparison_alias = index$historical_comparison_alias[[i]],
    go_term_count = nrow(all_i),
    fdr_significant_term_count = nrow(sig_i),
    numerator_enriched_count = sum(sig_i$NES > 0),
    denominator_enriched_count = sum(sig_i$NES < 0),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
summary_path <- write_enrichment_csv(summary, file.path(output_root, "GO_GSEA_comparison_summary.csv"))

selected_ids <- unique(significant$ID)
if (!length(selected_ids) && nrow(combined)) {
  selected_ids <- unique(head(combined$ID[order(combined$p.adjust, -abs(combined$NES), combined$ID)], 50L))
}
matrix_path <- file.path(output_root, "GO_GSEA_NES_matrix.csv")
plot_path <- NA_character_
if (length(selected_ids)) {
  matrix <- matrix(NA_real_, nrow = length(selected_ids), ncol = nrow(index),
                   dimnames = list(selected_ids, index$canonical_comparison))
  description <- stats::setNames(combined$Description, combined$ID)
  matrix_rows <- combined[combined$ID %in% selected_ids, , drop = FALSE]
  for (i in seq_len(nrow(matrix_rows))) {
    matrix[matrix_rows$ID[[i]], matrix_rows$canonical_comparison[[i]]] <- matrix_rows$NES[[i]]
  }
  matrix_table <- data.frame(ID = rownames(matrix), Description = unname(description[rownames(matrix)]), matrix,
                             check.names = FALSE, stringsAsFactors = FALSE)
  matrix_path <- write_enrichment_csv(matrix_table, matrix_path)
  if (isTRUE(config$plots_enabled) && requireNamespace("ggplot2", quietly = TRUE)) {
    long <- combined[combined$ID %in% selected_ids, c("ID", "canonical_comparison", "NES"), drop = FALSE]
    long$canonical_comparison <- factor(long$canonical_comparison, levels = index$canonical_comparison)
    long$ID <- factor(long$ID, levels = rev(selected_ids))
    plot <- ggplot2::ggplot(long, ggplot2::aes(canonical_comparison, ID, fill = NES)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
      ggplot2::labs(x = NULL, y = NULL, title = "Canonical animal-level GO GSEA", fill = "NES") +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    plot_path <- file.path(output_root, "GO_GSEA_NES_heatmap.png")
    ggplot2::ggsave(plot_path, plot = plot, width = 13, height = max(7, min(20, length(selected_ids) * 0.22)), dpi = 300)
    plot_path <- enrichment_normalize_path(plot_path, must_work = TRUE)
  }
} else {
  matrix_path <- write_enrichment_csv(
    data.frame(ID = character(), Description = character(), stringsAsFactors = FALSE), matrix_path
  )
}

compare_index <- data.frame(
  index[, c(
    "canonical_comparison", "canonical_contrast", "sample_class", "numerator_condition",
    "denominator_condition", "historical_comparison_alias", "mapped_input_path",
    "mapped_input_sha256", "source_split_sha256", "source_gct_sha256", "mapping_reference_sha256"
  ), drop = FALSE],
  enrichment_index_path = index_path,
  source_go_gsea_path = index$go_gsea_output,
  source_go_gsea_sha256 = vapply(index$go_gsea_output, enrichment_sha256, character(1)),
  go_term_count = summary$go_term_count,
  fdr_significant_term_count = summary$fdr_significant_term_count,
  combined_output_path = combined_path,
  significant_output_path = significant_path,
  summary_output_path = summary_path,
  matrix_output_path = matrix_path,
  heatmap_output_path = plot_path,
  run_timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_enrichment_csv(compare_index, output_index_path)
message("Canonical compareGO completed: ", output_index_path)
