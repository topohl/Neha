#!/usr/bin/env Rscript

# Compare the established learning and paired-CNO GO pathway profiles using
# canonical contrast metadata. Historical aliases are carried only as labels.

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1]]) else "04_differential_expression_enrichment/03_compare_pathways.r"
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "R", "analysis_labels.R"))
source(file.path(project_root, "R", "animal_level_enrichment_utils.R"))

config <- resolve_neha_enrichment_config()
compare_go_index_path <- Sys.getenv(
  "NEHA_COMPARE_GO_INDEX",
  unset = file.path(config$output_root, "compareGO", "indexCompareGO.csv")
)
compare_go_index_path <- neha_enrichment_normalize_path(compare_go_index_path, must_work = TRUE)
if (!neha_enrichment_path_is_within(compare_go_index_path, file.path(config$output_root, "compareGO"))) {
  stop("compare_pathways requires the canonical compareGO index.", call. = FALSE)
}
compare_index <- utils::read.csv(compare_go_index_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "canonical_comparison", "canonical_contrast", "sample_class", "numerator_condition",
  "denominator_condition", "historical_comparison_alias", "combined_output_path"
)
missing <- setdiff(required, names(compare_index))
if (length(missing) || nrow(compare_index) != 12L || anyDuplicated(compare_index$canonical_comparison)) {
  stop("compareGO index does not contain the 12 unique canonical comparisons.", call. = FALSE)
}
expected <- neha_primary_contrast_manifest()
if (!setequal(compare_index$canonical_comparison, expected$canonical_comparison)) {
  stop("compareGO index coverage does not match the shared Neha manifest.", call. = FALSE)
}
combined_paths <- unique(compare_index$combined_output_path)
if (length(combined_paths) != 1L) stop("compareGO index must identify one combined canonical GO table.", call. = FALSE)
combined_path <- neha_enrichment_normalize_path(combined_paths[[1]], must_work = TRUE)
if (!neha_enrichment_path_is_within(combined_path, file.path(config$output_root, "compareGO"))) {
  stop("Combined GO table resolves outside the canonical compareGO output.", call. = FALSE)
}
go <- utils::read.csv(combined_path, stringsAsFactors = FALSE, check.names = FALSE)

output_root <- file.path(config$output_root, "compare_pathways")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_index_path <- file.path(output_root, "indexComparePathways.csv")
if (file.exists(output_index_path) && !isTRUE(config$force)) {
  stop("compare_pathways output already exists; set NEHA_ENRICHMENT_FORCE=true to replace isolated canonical outputs.", call. = FALSE)
}

find_contrast <- function(sample_class, numerator, denominator) {
  rows <- which(
    compare_index$sample_class == sample_class &
      compare_index$numerator_condition == numerator &
      compare_index$denominator_condition == denominator
  )
  if (length(rows) != 1L) stop("Expected exactly one canonical contrast for ", sample_class, ": ", numerator, " vs ", denominator, call. = FALSE)
  compare_index[rows, , drop = FALSE]
}

pair_records <- lapply(unique(expected$sample_class), function(sample_class) {
  learning <- find_contrast(sample_class, "paired_veh", "unpaired_veh")
  cno <- find_contrast(sample_class, "paired_cno", "paired_veh")
  learning_go <- go[go$canonical_comparison == learning$canonical_comparison, c("ID", "Description", "NES", "p.adjust"), drop = FALSE]
  cno_go <- go[go$canonical_comparison == cno$canonical_comparison, c("ID", "Description", "NES", "p.adjust"), drop = FALSE]
  names(learning_go)[names(learning_go) %in% c("Description", "NES", "p.adjust")] <- c("Description_learning", "NES_learning", "padj_learning")
  names(cno_go)[names(cno_go) %in% c("Description", "NES", "p.adjust")] <- c("Description_cno", "NES_cno", "padj_cno")
  joined <- merge(learning_go, cno_go, by = "ID", all = TRUE, sort = TRUE)
  joined$Description <- ifelse(!is.na(joined$Description_learning), joined$Description_learning, joined$Description_cno)
  joined$sample_class <- sample_class
  joined$learning_canonical_comparison <- learning$canonical_comparison
  joined$cno_canonical_comparison <- cno$canonical_comparison
  joined$learning_historical_alias <- learning$historical_comparison_alias
  joined$cno_historical_alias <- cno$historical_comparison_alias
  joined$fdr_significant_learning <- is.finite(joined$padj_learning) & joined$padj_learning < config$fdr_threshold
  joined$fdr_significant_cno <- is.finite(joined$padj_cno) & joined$padj_cno < config$fdr_threshold
  joined$direction_relation <- ifelse(
    !is.finite(joined$NES_learning) | !is.finite(joined$NES_cno), "not_shared",
    ifelse(sign(joined$NES_learning) == sign(joined$NES_cno), "concordant", "opposed")
  )
  joined <- joined[, c(
    "sample_class", "ID", "Description", "NES_learning", "padj_learning", "NES_cno", "padj_cno",
    "fdr_significant_learning", "fdr_significant_cno", "direction_relation",
    "learning_canonical_comparison", "cno_canonical_comparison",
    "learning_historical_alias", "cno_historical_alias"
  ), drop = FALSE]
  profile_path <- write_neha_enrichment_csv(joined, file.path(output_root, paste0(sample_class, "_learning_vs_paired_cno_GO_profiles.csv")))
  overlap <- joined[joined$fdr_significant_learning | joined$fdr_significant_cno, , drop = FALSE]
  overlap_path <- write_neha_enrichment_csv(overlap, file.path(output_root, paste0(sample_class, "_FDR_pathway_overlap.csv")))

  shared <- is.finite(joined$NES_learning) & is.finite(joined$NES_cno)
  correlation <- if (sum(shared) >= 3L) stats::cor(joined$NES_learning[shared], joined$NES_cno[shared], method = "spearman") else NA_real_
  plot_path <- NA_character_
  if (isTRUE(config$plots_enabled) && any(shared) && requireNamespace("ggplot2", quietly = TRUE)) {
    plot_data <- joined[shared, , drop = FALSE]
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(NES_learning, NES_cno, colour = direction_relation)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
      ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
      ggplot2::geom_point(alpha = 0.65, size = 1.5) +
      ggplot2::scale_colour_manual(values = c(concordant = "#2166AC", opposed = "#B2182B")) +
      ggplot2::labs(
        title = paste(sample_class, "GO pathway profiles"),
        subtitle = "Canonical paired VEH vs unpaired VEH compared with paired CNO vs paired VEH",
        x = "Learning contrast NES", y = "Paired CNO contrast NES", colour = NULL
      ) + ggplot2::theme_minimal(base_size = 10)
    plot_path <- file.path(output_root, paste0(sample_class, "_learning_vs_paired_cno_GO_scatter.png"))
    ggplot2::ggsave(plot_path, plot = plot, width = 8, height = 7, dpi = 300)
    plot_path <- neha_enrichment_normalize_path(plot_path, must_work = TRUE)
  }
  data.frame(
    sample_class = sample_class,
    learning_canonical_comparison = learning$canonical_comparison,
    cno_canonical_comparison = cno$canonical_comparison,
    learning_historical_alias = learning$historical_comparison_alias,
    cno_historical_alias = cno$historical_comparison_alias,
    shared_GO_term_count = sum(shared),
    fdr_union_term_count = nrow(overlap),
    spearman_NES = correlation,
    paired_profile_path = profile_path,
    fdr_overlap_path = overlap_path,
    scatter_path = plot_path,
    compare_go_index_path = compare_go_index_path,
    source_combined_GO_path = combined_path,
    run_timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
})

pair_index <- do.call(rbind, pair_records)
write_neha_enrichment_csv(pair_index, output_index_path)
message("Canonical compare_pathways completed: ", output_index_path)
