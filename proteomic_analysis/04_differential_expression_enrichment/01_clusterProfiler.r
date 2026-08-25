#!/usr/bin/env Rscript

# Canonical animal-level enrichment for the Neha proteomics study.
#
# This stage consumes only the validated, forward MapThatProt manifest.  It does
# not discover files by name and it cannot fall back to the historical mapped
# or enrichment trees.  Positive log2fc values always mean higher abundance in
# the canonical numerator recorded in indexMappedComparisons.csv.

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1]]) else "04_differential_expression_enrichment/01_clusterProfiler.r"
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

source(file.path(project_root, "R", "analysis_labels.R"))
source(file.path(project_root, "R", "animal_level_enrichment_utils.R"))

required_packages <- c(
  "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "DOSE", "enrichplot",
  "GOSemSim", "ggplot2", "withr", "digest", "fgsea"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Canonical enrichment requires installed packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

config <- resolve_neha_enrichment_config()
mapped_index <- read_neha_enrichment_mapped_index(config$mapped_index, config$mapped_root)

dir.create(config$output_root, recursive = TRUE, showWarnings = FALSE)
per_comparison_root <- file.path(config$output_root, "per_comparison")
audit_root <- file.path(config$output_root, "audits")
dir.create(per_comparison_root, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)

index_path <- file.path(config$output_root, "indexEnrichmentComparisons.csv")
if (file.exists(index_path) && !isTRUE(config$force)) {
  stop(
    "Canonical enrichment index already exists. Refusing to overwrite it. ",
    "Set NEHA_ENRICHMENT_FORCE=true only after reviewing the existing isolated output root.",
    call. = FALSE
  )
}

timestamp_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
package_versions_path <- write_neha_enrichment_csv(
  neha_enrichment_package_versions(),
  file.path(audit_root, "package_database_versions.csv")
)

parameters <- data.frame(
  parameter = c(
    "ontology", "ranking_statistic", "rank_direction", "pvalue_cutoff",
    "qvalue_cutoff", "p_adjust_method", "fdr_threshold", "top_abs_log2fc",
    "min_gs_size", "max_gs_size", "simplify", "simplify_cutoff",
    "gsea_seed_base", "kegg_enabled", "plots_enabled",
    "duplicate_uniprot_rule", "ora_universe"
  ),
  value = c(
    config$ontology, "log2fc", "positive_is_higher_in_canonical_numerator",
    config$pvalue_cutoff, config$qvalue_cutoff, config$p_adjust_method,
    config$fdr_threshold, config$top_abs_log2fc, config$min_gs_size,
    config$max_gs_size, config$simplify, config$simplify_cutoff,
    config$gsea_seed_base, config$kegg_enabled, config$plots_enabled,
    "largest_absolute_log2fc;finite_preferred;ties_by_source_row_id_then_original_protein_id",
    "all_unique_successfully_mapped_measured_uniprot_accessions_per_comparison"
  ),
  stringsAsFactors = FALSE
)
parameters_path <- write_neha_enrichment_csv(parameters, file.path(audit_root, "run_parameters.csv"))

capture_warnings <- function(expression) {
  warnings <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

run_seeded <- function(expression, comparison, analysis) {
  seed <- derive_neha_enrichment_seed(config$gsea_seed_base, comparison, analysis)
  captured <- capture_warnings(with_neha_enrichment_seed(seed, expression))
  captured$seed <- seed
  captured
}

result_count <- function(result) nrow(neha_enrichment_result_table(result))

write_result <- function(result, path) {
  table <- neha_enrichment_result_table(result)
  write_neha_enrichment_csv(table, path)
}

save_dotplot <- function(result, path, title) {
  if (!isTRUE(config$plots_enabled) || result_count(result) == 0L) return(NA_character_)
  plot <- enrichplot::dotplot(result, showCategory = min(20L, result_count(result))) +
    ggplot2::ggtitle(title)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = 10, height = 7, dpi = 300)
  neha_enrichment_normalize_path(path, must_work = TRUE)
}

run_go_gsea <- function(rank, comparison) {
  run_seeded(
    clusterProfiler::gseGO(
      geneList = rank,
      ont = config$ontology,
      keyType = "UNIPROT",
      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
      exponent = 1,
      minGSSize = config$min_gs_size,
      maxGSSize = config$max_gs_size,
      eps = 1e-10,
      pvalueCutoff = config$pvalue_cutoff,
      pAdjustMethod = config$p_adjust_method,
      verbose = FALSE,
      seed = TRUE,
      by = "fgsea"
    ),
    comparison,
    paste0("GO_", config$ontology, "_GSEA")
  )
}

run_go_ora <- function(gene, universe) {
  if (!length(gene)) return(list(value = NULL, warnings = character()))
  capture_warnings(clusterProfiler::enrichGO(
    gene = gene,
    universe = universe,
    OrgDb = org.Mm.eg.db::org.Mm.eg.db,
    keyType = "UNIPROT",
    ont = config$ontology,
    pAdjustMethod = config$p_adjust_method,
    pvalueCutoff = config$pvalue_cutoff,
    qvalueCutoff = config$qvalue_cutoff,
    minGSSize = config$min_gs_size,
    maxGSSize = config$max_gs_size,
    readable = FALSE
  ))
}

build_kegg_rank <- function(collapsed) {
  conversion <- suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db,
    keys = unique(as.character(collapsed$uniprot_accession)),
    keytype = "UNIPROT",
    columns = "ENTREZID"
  ))
  conversion <- conversion[!is.na(conversion$ENTREZID) & nzchar(conversion$ENTREZID), , drop = FALSE]
  conversion <- conversion[order(conversion$UNIPROT, conversion$ENTREZID, method = "radix"), , drop = FALSE]
  conversion <- conversion[!duplicated(conversion$UNIPROT), , drop = FALSE]
  candidate <- merge(
    collapsed[, c("uniprot_accession", "original_protein_id", "source_row_id", "log2fc"), drop = FALSE],
    conversion,
    by.x = "uniprot_accession", by.y = "UNIPROT", all = FALSE, sort = FALSE
  )
  candidate$log2fc <- suppressWarnings(as.numeric(candidate$log2fc))
  candidate <- candidate[is.finite(candidate$log2fc), , drop = FALSE]
  candidate <- candidate[order(
    candidate$ENTREZID, -abs(candidate$log2fc), candidate$source_row_id,
    candidate$uniprot_accession, method = "radix"
  ), , drop = FALSE]
  candidate$selected_representative <- !duplicated(candidate$ENTREZID)
  selected <- candidate[candidate$selected_representative, , drop = FALSE]
  selected <- selected[order(-selected$log2fc, selected$ENTREZID, method = "radix"), , drop = FALSE]
  rank <- selected$log2fc
  names(rank) <- selected$ENTREZID
  if (anyDuplicated(names(rank)) || any(!is.finite(rank))) stop("Invalid KEGG rank after Entrez collapse.", call. = FALSE)
  candidate$selection_rule <- "one_Entrez_per_UniProt_then_largest_absolute_log2fc_per_Entrez;ties_by_source_row_id_then_UniProt"
  list(rank = rank, audit = candidate, conversion = conversion)
}

run_kegg_gsea <- function(rank, comparison, analysis = "KEGG_GSEA") {
  if (!length(rank)) return(list(value = NULL, warnings = character(), seed = NA_integer_))
  run_seeded(
    clusterProfiler::gseKEGG(
      geneList = rank,
      organism = "mmu",
      keyType = "ncbi-geneid",
      exponent = 1,
      minGSSize = config$min_gs_size,
      maxGSSize = config$max_gs_size,
      eps = 1e-10,
      pvalueCutoff = config$pvalue_cutoff,
      pAdjustMethod = config$p_adjust_method,
      verbose = FALSE,
      seed = TRUE,
      by = "fgsea"
    ),
    comparison,
    analysis
  )
}

run_custom_gsea <- function(rank, genes, comparison, analysis, term_name) {
  if (!length(genes)) return(list(value = NULL, warnings = character(), seed = NA_integer_))
  term2gene <- data.frame(term = term_name, gene = unique(genes), stringsAsFactors = FALSE)
  run_seeded(
    clusterProfiler::GSEA(
      geneList = rank,
      TERM2GENE = term2gene,
      exponent = 1,
      minGSSize = 1,
      maxGSSize = 500,
      eps = 1e-10,
      pvalueCutoff = config$pvalue_cutoff,
      pAdjustMethod = config$p_adjust_method,
      verbose = FALSE,
      seed = TRUE,
      by = "fgsea"
    ),
    comparison,
    analysis
  )
}

empty_warning_table <- function() data.frame(
  canonical_comparison = character(), analysis = character(), warning = character(),
  stringsAsFactors = FALSE
)

append_warnings <- function(table, comparison, analysis, warnings) {
  if (!length(warnings)) return(table)
  rbind(table, data.frame(
    canonical_comparison = comparison,
    analysis = analysis,
    warning = as.character(warnings),
    stringsAsFactors = FALSE
  ))
}

write_comparison_manifest <- function(record, path) {
  values <- vapply(record, function(value) {
    if (!length(value) || all(is.na(value))) "" else paste(as.character(value), collapse = ";")
  }, character(1))
  write_neha_enrichment_csv(
    data.frame(field = names(values), value = unname(values), stringsAsFactors = FALSE),
    path
  )
}

process_comparison <- function(index_row) {
  comparison <- as.character(index_row$canonical_comparison)
  comparison_dir <- file.path(per_comparison_root, comparison)
  audit_dir <- file.path(audit_root, comparison)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  started <- Sys.time()
  warning_rows <- empty_warning_table()

  record <- list(
    canonical_comparison = comparison,
    canonical_contrast = as.character(index_row$canonical_contrast),
    sample_class = as.character(index_row$sample_class),
    numerator_condition = as.character(index_row$numerator_condition),
    denominator_condition = as.character(index_row$denominator_condition),
    historical_comparison_alias = as.character(index_row$historical_comparison_alias),
    mapping_direction = as.character(index_row$mapping_direction),
    mapped_input_path = as.character(index_row$mapped_input_path),
    mapped_input_sha256 = as.character(index_row$mapped_output_sha256),
    mapped_index_path = as.character(index_row$mapped_index_path),
    mapped_index_sha256 = as.character(index_row$mapped_index_sha256),
    source_split_index_path = as.character(index_row$source_split_index_path),
    source_split_index_sha256 = as.character(index_row$source_split_index_sha256),
    source_split_sha256 = as.character(index_row$source_split_sha256),
    source_gct_sha256 = as.character(index_row$source_gct_sha256),
    mapping_reference_path = as.character(index_row$mapping_reference_path),
    mapping_reference_version = as.character(index_row$mapping_reference_version),
    mapping_reference_snapshot_date_utc = as.character(index_row$mapping_reference_snapshot_date_utc),
    mapping_reference_sha256 = as.character(index_row$mapping_reference_sha256),
    ranking_statistic = "log2fc",
    rank_direction = "positive_is_higher_in_canonical_numerator",
    duplicate_uniprot_rule = "largest_absolute_log2fc;finite_preferred;ties_by_source_row_id_then_original_protein_id",
    ora_universe_definition = "all_unique_successfully_mapped_measured_uniprot_accessions",
    ontology = config$ontology,
    pvalue_cutoff = config$pvalue_cutoff,
    qvalue_cutoff = config$qvalue_cutoff,
    p_adjust_method = config$p_adjust_method,
    fdr_threshold = config$fdr_threshold,
    top_abs_log2fc_threshold = config$top_abs_log2fc,
    min_gs_size = config$min_gs_size,
    max_gs_size = config$max_gs_size,
    simplify_enabled = config$simplify,
    simplify_cutoff = config$simplify_cutoff,
    execution_status = "failed",
    error_message = NA_character_
  )

  tryCatch({
    mapped <- read_neha_enrichment_mapped_file(
      record$mapped_input_path,
      expected_rows = index_row$n_output_mapped_rows,
      expected_sha256 = record$mapped_input_sha256
    )
    collapsed_info <- collapse_neha_enrichment_accessions(mapped, "log2fc")
    collapsed <- collapsed_info$collapsed
    rank_info <- build_neha_gsea_rank(collapsed, "log2fc")
    ora_sets <- build_neha_ora_sets(collapsed, config$fdr_threshold, config$top_abs_log2fc)

    record$n_source_protein_rows <- nrow(mapped)
    record$n_mapped_protein_rows <- nrow(mapped)
    record$n_unique_mapped_uniprot <- length(ora_sets$universe)
    record$n_duplicate_rows_collapsed <- collapsed_info$n_duplicate_rows_collapsed
    record$n_duplicated_uniprot_accessions <- collapsed_info$n_duplicated_accessions
    record$rank_vector_size <- length(rank_info$rank)
    record$ora_universe_size <- length(ora_sets$universe)
    record$fdr_significant_protein_count <- length(ora_sets$all_significant)
    record$fdr_significant_up_count <- length(ora_sets$up_significant)
    record$fdr_significant_down_count <- length(ora_sets$down_significant)
    record$top_abs_log2fc_count <- length(ora_sets$top_abs_log2fc)

    record$duplicate_uniprot_audit <- write_neha_enrichment_csv(
      collapsed_info$duplicate_audit, file.path(audit_dir, "duplicate_uniprot_audit.csv")
    )
    record$collapsed_mapped_audit <- write_neha_enrichment_csv(
      collapsed, file.path(audit_dir, "mapped_uniprot_collapsed.csv")
    )
    record$gsea_rank_audit <- write_neha_enrichment_csv(
      rank_info$audit, file.path(audit_dir, "gsea_rank_audit.csv")
    )
    record$ora_membership_audit <- write_neha_enrichment_csv(
      data.frame(
        uniprot_accession = ora_sets$universe,
        fdr_significant = ora_sets$universe %in% ora_sets$all_significant,
        fdr_significant_up = ora_sets$universe %in% ora_sets$up_significant,
        fdr_significant_down = ora_sets$universe %in% ora_sets$down_significant,
        top_abs_log2fc = ora_sets$universe %in% ora_sets$top_abs_log2fc,
        stringsAsFactors = FALSE
      ),
      file.path(audit_dir, "ora_measured_universe_audit.csv")
    )

    go_gsea <- run_go_gsea(rank_info$rank, comparison)
    warning_rows <- append_warnings(warning_rows, comparison, paste0("GO_", config$ontology, "_GSEA"), go_gsea$warnings)
    go_gsea_table <- neha_enrichment_result_table(go_gsea$value)
    record$go_gsea_seed <- go_gsea$seed
    record$go_gsea_term_count <- nrow(go_gsea_table)
    record$go_gsea_fdr_term_count <- sum(is.finite(go_gsea_table$p.adjust) & go_gsea_table$p.adjust < config$fdr_threshold)
    record$go_gsea_output <- write_result(go_gsea$value, file.path(comparison_dir, paste0("GSEA_GO_", config$ontology, ".csv")))
    record$go_gsea_plot <- save_dotplot(
      go_gsea$value,
      file.path(comparison_dir, paste0("GSEA_GO_", config$ontology, "_dotplot.png")),
      paste(comparison, "GO", config$ontology, "GSEA")
    )

    if (isTRUE(config$simplify) && nrow(go_gsea_table)) {
      simplified <- capture_warnings(clusterProfiler::simplify(
        go_gsea$value, cutoff = config$simplify_cutoff, by = "p.adjust",
        select_fun = min, measure = "Wang", semData = NULL
      ))
      warning_rows <- append_warnings(warning_rows, comparison, "GO_GSEA_simplify", simplified$warnings)
      record$go_gsea_simplified_term_count <- result_count(simplified$value)
      record$go_gsea_simplified_output <- write_result(
        simplified$value, file.path(comparison_dir, paste0("GSEA_GO_", config$ontology, "_simplified.csv"))
      )
    } else {
      record$go_gsea_simplified_term_count <- NA_integer_
      record$go_gsea_simplified_output <- NA_character_
    }

    ora_definitions <- list(
      fdr_all = ora_sets$all_significant,
      fdr_up = ora_sets$up_significant,
      fdr_down = ora_sets$down_significant,
      top_abs_log2fc = ora_sets$top_abs_log2fc
    )
    for (ora_name in names(ora_definitions)) {
      ora <- run_go_ora(ora_definitions[[ora_name]], ora_sets$universe)
      warning_rows <- append_warnings(warning_rows, comparison, paste0("GO_ORA_", ora_name), ora$warnings)
      record[[paste0("go_ora_", ora_name, "_term_count")]] <- result_count(ora$value)
      record[[paste0("go_ora_", ora_name, "_output")]] <- write_result(
        ora$value, file.path(comparison_dir, paste0("ORA_GO_", config$ontology, "_", ora_name, ".csv"))
      )
      if (isTRUE(config$simplify) && result_count(ora$value)) {
        simplified_ora <- capture_warnings(clusterProfiler::simplify(
          ora$value, cutoff = config$simplify_cutoff, by = "p.adjust",
          select_fun = min, measure = "Wang", semData = NULL
        ))
        warning_rows <- append_warnings(warning_rows, comparison, paste0("GO_ORA_", ora_name, "_simplify"), simplified_ora$warnings)
        record[[paste0("go_ora_", ora_name, "_simplified_term_count")]] <- result_count(simplified_ora$value)
        record[[paste0("go_ora_", ora_name, "_simplified_output")]] <- write_result(
          simplified_ora$value,
          file.path(comparison_dir, paste0("ORA_GO_", config$ontology, "_", ora_name, "_simplified.csv"))
        )
      } else {
        record[[paste0("go_ora_", ora_name, "_simplified_term_count")]] <- NA_integer_
        record[[paste0("go_ora_", ora_name, "_simplified_output")]] <- NA_character_
      }
    }

    if (isTRUE(config$kegg_enabled)) {
      kegg <- build_kegg_rank(collapsed)
      record$kegg_rank_vector_size <- length(kegg$rank)
      record$kegg_mapping_audit <- write_neha_enrichment_csv(
        kegg$audit, file.path(audit_dir, "kegg_entrez_collapse_audit.csv")
      )
      kegg_gsea <- run_kegg_gsea(kegg$rank, comparison)
      warning_rows <- append_warnings(warning_rows, comparison, "KEGG_GSEA", kegg_gsea$warnings)
      record$kegg_gsea_seed <- kegg_gsea$seed
      record$kegg_gsea_term_count <- result_count(kegg_gsea$value)
      record$kegg_gsea_output <- write_result(kegg_gsea$value, file.path(comparison_dir, "GSEA_KEGG.csv"))
      record$kegg_gsea_plot <- save_dotplot(
        kegg_gsea$value, file.path(comparison_dir, "GSEA_KEGG_dotplot.png"),
        paste(comparison, "KEGG GSEA")
      )

      if (length(config$selected_uniprot)) {
        selected_entrez <- unique(kegg$audit$ENTREZID[kegg$audit$uniprot_accession %in% config$selected_uniprot])
        selected_rank <- kegg$rank[names(kegg$rank) %in% selected_entrez]
        selected_rank <- sort(selected_rank, decreasing = TRUE)
        selected_gsea <- run_kegg_gsea(selected_rank, comparison, "KEGG_selected_UniProt_GSEA")
        warning_rows <- append_warnings(warning_rows, comparison, "KEGG_selected_UniProt_GSEA", selected_gsea$warnings)
        record$custom_selected_uniprot_term_count <- result_count(selected_gsea$value)
        record$custom_selected_uniprot_output <- write_result(
          selected_gsea$value, file.path(comparison_dir, "GSEA_KEGG_selected_uniprot.csv")
        )
      } else {
        record$custom_selected_uniprot_term_count <- 0L
        record$custom_selected_uniprot_output <- NA_character_
      }

      if (length(config$path_ids)) {
        if (!requireNamespace("pathview", quietly = TRUE)) {
          stop("NEHA_ENRICHMENT_PATH_IDS was configured but package 'pathview' is unavailable.", call. = FALSE)
        }
        pathview_dir <- file.path(comparison_dir, "pathview")
        dir.create(pathview_dir, recursive = TRUE, showWarnings = FALSE)
        withr::with_dir(pathview_dir, {
          for (path_id in config$path_ids) {
            pathview_result <- capture_warnings(pathview::pathview(
              gene.data = kegg$rank, pathway.id = path_id, species = "mmu",
              gene.idtype = "entrez", out.suffix = comparison, kegg.native = TRUE
            ))
            warning_rows <- append_warnings(
              warning_rows, comparison, paste0("pathview_", path_id), pathview_result$warnings
            )
          }
        })
        record$custom_pathview_output_directory <- neha_enrichment_normalize_path(pathview_dir, must_work = TRUE)
      } else {
        record$custom_pathview_output_directory <- NA_character_
      }
    } else {
      record$kegg_rank_vector_size <- 0L
      record$kegg_gsea_term_count <- 0L
      record$kegg_gsea_output <- NA_character_
      record$kegg_mapping_audit <- NA_character_
      record$custom_selected_uniprot_term_count <- 0L
      record$custom_selected_uniprot_output <- NA_character_
      record$custom_pathview_output_directory <- NA_character_
    }

    if (length(config$nk3r_genes)) {
      nk3r_gsea <- run_custom_gsea(
        rank_info$rank, config$nk3r_genes, comparison,
        "custom_NK3R_GSEA", "NK3R-signalling"
      )
      warning_rows <- append_warnings(warning_rows, comparison, "custom_NK3R_GSEA", nk3r_gsea$warnings)
      record$custom_nk3r_seed <- nk3r_gsea$seed
      record$custom_nk3r_member_count <- sum(names(rank_info$rank) %in% config$nk3r_genes)
      record$custom_nk3r_term_count <- result_count(nk3r_gsea$value)
      record$custom_nk3r_output <- write_result(
        nk3r_gsea$value, file.path(comparison_dir, "custom_NK3R_GSEA.csv")
      )
    } else {
      record$custom_nk3r_output <- NA_character_
      record$custom_nk3r_member_count <- 0L
      record$custom_nk3r_term_count <- 0L
      record$custom_nk3r_seed <- NA_integer_
    }

    record$warning_count <- nrow(warning_rows)
    record$warnings_path <- write_neha_enrichment_csv(warning_rows, file.path(audit_dir, "warnings.csv"))
    record$execution_status <- "success"
  }, error = function(e) {
    record$error_message <<- conditionMessage(e)
    record$warning_count <<- nrow(warning_rows)
    record$warnings_path <<- write_neha_enrichment_csv(warning_rows, file.path(audit_dir, "warnings.csv"))
  })

  record$runtime_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  record$run_timestamp_utc <- timestamp_utc
  record$package_versions_path <- package_versions_path
  record$run_parameters_path <- parameters_path
  manifest_path <- file.path(audit_dir, "comparison_manifest.csv")
  record$comparison_manifest_path <- neha_enrichment_normalize_path(manifest_path)
  write_comparison_manifest(record, manifest_path)
  record$comparison_manifest_path <- neha_enrichment_normalize_path(manifest_path, must_work = TRUE)
  as.data.frame(record, stringsAsFactors = FALSE, check.names = FALSE)
}

message("Canonical mapped index: ", config$mapped_index)
message("Canonical enrichment output: ", config$output_root)
message("Processing exactly ", nrow(mapped_index), " forward comparisons...")

records <- lapply(seq_len(nrow(mapped_index)), function(i) {
  comparison <- mapped_index$canonical_comparison[[i]]
  message("  [", i, "/", nrow(mapped_index), "] ", comparison)
  process_comparison(mapped_index[i, , drop = FALSE])
})
record_columns <- unique(unlist(lapply(records, names), use.names = FALSE))
records <- lapply(records, function(record) {
  missing <- setdiff(record_columns, names(record))
  for (field in missing) record[[field]] <- NA
  record[, record_columns, drop = FALSE]
})
enrichment_index <- do.call(rbind, records)
rownames(enrichment_index) <- NULL
write_neha_enrichment_csv(enrichment_index, index_path)

run_contract <- data.frame(
  field = c(
    "run_timestamp_utc", "mapped_index_path", "mapped_index_sha256", "output_root",
    "expected_primary_comparisons", "observed_primary_comparisons", "successful_comparisons",
    "mapping_direction", "rank_direction", "ora_universe", "historical_outputs_written",
    "package_versions_path", "run_parameters_path"
  ),
  value = c(
    timestamp_utc, config$mapped_index, neha_enrichment_sha256(config$mapped_index), config$output_root,
    12L, nrow(enrichment_index), sum(enrichment_index$execution_status == "success"),
    "forward", "positive_is_higher_in_canonical_numerator",
    "all_unique_successfully_mapped_measured_uniprot_accessions_per_comparison",
    "false", package_versions_path, parameters_path
  ),
  stringsAsFactors = FALSE
)
write_neha_enrichment_csv(run_contract, file.path(audit_root, "run_contract_manifest.csv"))

failed <- enrichment_index$execution_status != "success"
if (any(failed)) {
  stop(
    "Canonical enrichment failed for: ",
    paste(enrichment_index$canonical_comparison[failed], collapse = ", "),
    ". See indexEnrichmentComparisons.csv and per-comparison audits.",
    call. = FALSE
  )
}
validate_neha_enrichment_index(enrichment_index, require_files = TRUE)
message("Canonical animal-level enrichment completed: ", index_path)
