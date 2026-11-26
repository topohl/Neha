#!/usr/bin/env Rscript

# Post-hoc Pairing x CNO interaction under the per-compartment (primary) model.
#
# Question: does the CNO effect differ between paired and unpaired animals? Contrast
#   (paired_cno - paired_veh) - (unpaired_cno - unpaired_veh)
# fitted within each sample class: limma ~0 + condition over that class's 12 animals
# (n = 3 per condition, residual df 8), contrasts.fit, eBayes(robust = TRUE, trend = FALSE),
# Benjamini-Hochberg within the contrast. This is the model of 01_preprocessing/02b_limma_stat_gct.r
# (within_compartment mode); the script checks that its paired-CNO fit reproduces the mapped
# primary statistics. Effects (logFC / log2fc columns) are standardized abundance differences in
# SD units, not log2 ratios. Positive = the CNO effect is more positive in paired animals.
#
# Gene sets: GSEA GO-BP, GO-CC and KEGG ranked by the interaction moderated t, with the settings,
# UniProt collapse, KEGG Entrez collapse and per-job seeds of 01_clusterProfiler.r
# (seed label "<class>_pairing_x_cno_interaction"); the statistics enter GSEA at the 15 significant
# digits the pipeline's text inputs carry. Extra base seeds give the seed range of the GO-BP and
# GO-CC counts.
#
# Replaces posthoc_pairing_x_cno_interaction/ and pairing_cno_interaction_*.csv of the pooled run
# (03_output/inferential_checks/final_inferential_checks_20260826) in the same layout, with GO-CC
# added, so the release stages can read them from the new folder. Also writes the diagnostics used
# by the report of 06_mcherry_learning_robustness.r (compartment-wide shift, plate weight).
#
# Writes only into PROTEOMICS_INFERENTIAL_OUTPUT_ROOT and refuses a non-empty folder unless
# PROTEOMICS_INFERENTIAL_FORCE=true. Run before 06_mcherry_learning_robustness.r.

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1]]) else "04_differential_expression_enrichment/05_pairing_cno_interaction.r"
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

source(file.path(project_root, "R", "analysis_labels.R"))
source(file.path(project_root, "R", "project_path_utils.R"))
source(file.path(project_root, "R", "animal_level_enrichment_utils.R"))
source(file.path(project_root, "R", "inferential_checks_utils.R"))

option_or_env <- function(option_name, env_name, default) {
  value <- getOption(option_name)
  if (!is.null(value) && nzchar(trimws(as.character(value)))) return(as.character(value))
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(trimws(value))) return(value)
  default
}

animal_root <- file.path(PROTEOMICS_DATA_ROOT, "02_data", "animal_level")
input_gct <- project_normalize_path(option_or_env(
  "proteomics.inferential_input_gct", "PROTEOMICS_INFERENTIAL_INPUT_GCT",
  file.path(animal_root, "input_gct", "neha_protigy_input_animal_level_primary.gct")
), must_work = TRUE)
mapped_root <- project_normalize_path(option_or_env(
  "proteomics.inferential_mapped_root", "PROTEOMICS_INFERENTIAL_MAPPED_ROOT",
  file.path(animal_root, "within_compartment_20261124", "mapped")
), must_work = TRUE)
mapped_index_path <- file.path(mapped_root, "indexMappedComparisons.csv")
# Plate per animal (design only; release stage 01 reads the same file).
design_file <- project_normalize_path(option_or_env(
  "proteomics.inferential_design_file", "PROTEOMICS_INFERENTIAL_DESIGN_FILE",
  file.path(PROTEOMICS_DATA_ROOT, "03_output", "inferential_checks", "final_inferential_checks_20260826", "design_identifiability_audit.csv")
), must_work = TRUE)
output_root <- project_normalize_path(option_or_env(
  "proteomics.inferential_output_root", "PROTEOMICS_INFERENTIAL_OUTPUT_ROOT",
  file.path(PROTEOMICS_DATA_ROOT, "03_output", "inferential_checks", "within_compartment_20261124")
))
force <- inferential_flag(option_or_env("proteomics.inferential_force", "PROTEOMICS_INFERENTIAL_FORCE", "false"), "PROTEOMICS_INFERENTIAL_FORCE")
run_gsea <- inferential_flag(option_or_env("proteomics.inferential_gsea", "PROTEOMICS_INFERENTIAL_GSEA", "true"), "PROTEOMICS_INFERENTIAL_GSEA")
workers <- inferential_workers(option_or_env("proteomics.inferential_gsea_workers", "PROTEOMICS_INFERENTIAL_GSEA_WORKERS", "3"), "PROTEOMICS_INFERENTIAL_GSEA_WORKERS")
seed_range_bases <- option_or_env("proteomics.inferential_seed_range_bases", "PROTEOMICS_INFERENTIAL_SEED_RANGE_BASES", "1,2,3,4")
seed_range_bases <- if (tolower(trimws(seed_range_bases)) == "none") integer() else as.integer(trimws(strsplit(seed_range_bases, ",", fixed = TRUE)[[1]]))
if (anyNA(seed_range_bases)) stop("PROTEOMICS_INFERENTIAL_SEED_RANGE_BASES must be whole numbers separated by commas, or none.", call. = FALSE)

required_packages <- c("limma", "digest", if (run_gsea) c("clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "DOSE", "fgsea", "BiocParallel", "withr"))
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages before running this stage: ", paste(missing_packages, collapse = ", "), call. = FALSE)

invisible(validate_inferential_output_root(output_root))
if (!inferential_dir_is_empty(output_root) && !force) {
  stop("Output folder is not empty; set PROTEOMICS_INFERENTIAL_FORCE=true to write into it: ", output_root, call. = FALSE)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- project_normalize_path(output_root, must_work = TRUE)
config <- resolve_enrichment_config(mapped_root = mapped_root, mapped_index = mapped_index_path, output_root = output_root)
if (!identical(config$gsea_rank, "t")) stop("The interaction GSEA ranks by the moderated t (PROTEOMICS_ENRICHMENT_GSEA_RANK must be t).", call. = FALSE)
gsea_settings <- config[c("min_gs_size", "max_gs_size", "pvalue_cutoff", "p_adjust_method", "fdr_threshold")]
timestamp_utc <- format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")

# ---- inputs -------------------------------------------------------------------------------
input <- read_animal_level_input_gct(input_gct)
values <- input$values
annotation <- input$annotation

mapped_index <- read_enrichment_mapped_index(mapped_index_path, mapped_root)  # verifies the 12 mapped-file hashes
mapped_tables <- lapply(seq_len(nrow(mapped_index)), function(i) read_enrichment_mapped_file(mapped_index$mapped_input_path[[i]]))
names(mapped_tables) <- mapped_index$canonical_comparison
mapping_columns <- c("original_protein_id", "uniprot_accession", "source_row_id")
for (k in names(mapped_tables)[-1]) {
  if (!identical(mapped_tables[[k]][, mapping_columns], mapped_tables[[1]][, mapping_columns])) {
    stop("Mapped tables differ in their protein-to-UniProt mapping: ", k, call. = FALSE)
  }
}
mapping <- mapped_tables[[1]]
if (anyDuplicated(mapping$original_protein_id)) stop("Mapped table has duplicated protein ids.", call. = FALSE)
full_matrix <- nrow(values) == as.integer(mapped_index$n_input_proteins[[1]])

design <- utils::read.csv(design_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample_class", "AnimalID", "condition", "Plate") %in% names(design))) {
  stop("Design file needs sample_class, AnimalID, condition and Plate columns: ", design_file, call. = FALSE)
}

# ---- protein level ------------------------------------------------------------------------
posthoc_root <- file.path(output_root, "posthoc_pairing_x_cno_interaction")
fits <- list()
protein_tables <- list()
diagnostics <- list()
plate_rows <- list()
for (cls in sample_classes) {
  a <- annotation[annotation$sample_class == cls, , drop = FALSE]
  fit <- fit_condition_contrast(values[, a$column, drop = FALSE], a$condition, PAIRING_CNO_INTERACTION_CONTRAST)
  expected_weights <- c(paired_cno = 1, paired_veh = -1, unpaired_cno = -1, unpaired_veh = 1)
  if (!identical(fit$weights[names(expected_weights)], expected_weights)) stop("Unexpected interaction contrast weights.", call. = FALSE)
  if (fit$df_residual != 8L) stop("Expected residual df 8 in the per-compartment model; got ", fit$df_residual, call. = FALSE)
  fits[[cls]] <- fit
  tab <- fit$table

  # The paired-CNO contrast of the same fit must reproduce the mapped primary statistics.
  check_diff <- NA_real_
  if (full_matrix) {
    primary <- mapped_tables[[paste0(cls, "_paired_cno_over_", cls, "_paired_veh")]]
    check <- fit_condition_contrast(values[, a$column, drop = FALSE], a$condition, "paired_cno - paired_veh")$table
    check_diff <- max(abs(check[primary$original_protein_id, "t"] - primary$t))
    if (!is.finite(check_diff) || check_diff > 1e-8) {
      stop("The per-compartment fit does not reproduce the mapped primary t for ", cls, " (max diff ", signif(check_diff, 3), ").", call. = FALSE)
    }
  }

  ids <- rownames(tab)
  mapped_row <- match(ids, mapping$original_protein_id)
  protein <- data.frame(
    sample_class = cls,
    original_protein_id = ids,
    Description = unname(input$description[ids]),
    uniprot_accession = mapping$uniprot_accession[mapped_row],
    mapped = !is.na(mapped_row),
    logFC = tab$logFC, AveExpr = tab$AveExpr, t = tab$t, P.Value = tab$P.Value,
    adj.P.Val = tab$adj.P.Val, B = tab$B,
    df_residual = fit$df_residual, df_prior_median = fit$df_prior_median,
    stringsAsFactors = FALSE
  )
  protein <- protein[order(protein$P.Value, method = "radix"), , drop = FALSE]
  rownames(protein) <- NULL
  protein_tables[[cls]] <- protein
  write_enrichment_csv(protein, file.path(posthoc_root, cls, "interaction_protein_table.csv"))

  d <- design[design$sample_class == cls, , drop = FALSE]
  d <- d[match(a$animal, d$AnimalID), , drop = FALSE]
  if (anyNA(d$AnimalID) || !identical(d$condition, a$condition)) stop("Design file does not match the input animals for ", cls, call. = FALSE)
  plate_weight <- contrast_plate_weight(fit$weights, a$condition, d$Plate)
  plate_rows[[cls]] <- do.call(rbind, lapply(names(expected_weights), function(cond) {
    k <- a$condition == cond
    data.frame(
      sample_class = cls, condition = cond, contrast_weight = unname(fit$weights[[cond]]),
      n_animals = sum(k), n_plate1 = sum(d$Plate[k] == "Plate1"), n_plate2 = sum(d$Plate[k] == "Plate2"),
      plate2_fraction = mean(d$Plate[k] == "Plate2"),
      plate1_animals = paste(a$animal[k][d$Plate[k] == "Plate1"], collapse = ";"),
      plate2_animals = paste(a$animal[k][d$Plate[k] == "Plate2"], collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  diagnostics[[cls]] <- data.frame(
    sample_class = cls, n_tested = nrow(tab), n_mapped = sum(!is.na(mapped_row)),
    n_fdr05 = sum(tab$adj.P.Val < 0.05), n_fdr10 = sum(tab$adj.P.Val < 0.10), min_adj_p = min(tab$adj.P.Val),
    df_residual = fit$df_residual, df_prior_median = fit$df_prior_median, s2_prior = fit$s2_prior,
    protein_median_t = stats::median(tab$t), protein_pct_t_positive = 100 * mean(tab$t > 0),
    protein_median_effect = stats::median(tab$logFC),
    plate2_minus_plate1_weight = plate_weight,
    check_max_abs_diff_t_vs_primary = check_diff,
    stringsAsFactors = FALSE
  )
}
plate_composition <- do.call(rbind, plate_rows)
rownames(plate_composition) <- NULL
diagnostics <- do.call(rbind, diagnostics)
rownames(diagnostics) <- NULL

# Release-shaped files: the columns of final_inferential_checks_20260826.
protein_summary <- do.call(rbind, lapply(sample_classes, function(cls) {
  x <- protein_tables[[cls]]
  data.frame(
    sample_class = cls, n_tested = nrow(x), n_nominal_p05 = sum(x$P.Value < 0.05),
    pct_nominal_p05 = round(100 * mean(x$P.Value < 0.05), 2), n_fdr05 = sum(x$adj.P.Val < 0.05),
    n_positive_interaction_fdr = sum(x$adj.P.Val < 0.05 & x$logFC > 0),
    n_negative_interaction_fdr = sum(x$adj.P.Val < 0.05 & x$logFC < 0),
    median_interaction_log2fc = round(stats::median(x$logFC), 4), IQR_interaction_log2fc = round(stats::IQR(x$logFC), 4),
    stringsAsFactors = FALSE
  )
}))
flagged_genes <- list(cfos = c("Chmp2a", "Cntfr"))  # named in the manuscript from the paired-CNO contrast
top_proteins <- do.call(rbind, lapply(sample_classes, function(cls) {
  x <- protein_tables[[cls]]
  x <- x[order(-abs(x$t), method = "radix"), , drop = FALSE]
  top <- x[seq_len(min(10L, nrow(x))), , drop = FALSE]
  top$rank <- "top10_by_abs_t"
  flagged <- x[x$Description %in% flagged_genes[[cls]], , drop = FALSE]
  flagged$rank <- rep("flagged_gene_from_simple_contrast", nrow(flagged))
  y <- rbind(top, flagged)
  data.frame(sample_class = cls, rank = y$rank, original_protein_id = y$original_protein_id, Description = y$Description,
             uniprot_accession = y$uniprot_accession, logFC = y$logFC, t = y$t, P.Value = y$P.Value, adj.P.Val = y$adj.P.Val,
             stringsAsFactors = FALSE)
}))
invisible(write_enrichment_csv(protein_summary, file.path(output_root, "pairing_cno_interaction_protein_summary.csv")))
invisible(write_enrichment_csv(top_proteins, file.path(output_root, "pairing_cno_interaction_top_proteins.csv")))
invisible(write_enrichment_csv(plate_composition, file.path(output_root, "pairing_cno_interaction_plate_composition.csv")))

# ---- gene sets ----------------------------------------------------------------------------
gsea_job_rows <- NULL
if (run_gsea) {
  ontologies <- c("BP", "CC", "KEGG")
  file_of <- c(BP = "GSEA_GO_BP_interaction.csv", CC = "GSEA_GO_CC_interaction.csv", KEGG = "GSEA_KEGG_interaction.csv")
  jobs <- list()
  for (cls in sample_classes) {
    template <- mapped_tables[[paste0(cls, "_paired_cno_over_", cls, "_paired_veh")]]
    collapsed <- collapse_enrichment_accessions(inferential_mapped_stats(template, fits[[cls]]$table), "log2fc")$collapsed
    rank_info <- build_gsea_rank(collapsed, "t", "canonical")
    write_enrichment_csv(rank_info$audit, file.path(posthoc_root, cls, "interaction_gsea_rank_audit.csv"))
    label <- paste0(cls, "_pairing_x_cno_interaction")
    for (ont in ontologies) {
      jobs[[length(jobs) + 1L]] <- inferential_gsea_job_spec(
        paste(cls, ont, config$gsea_seed_base, sep = "|"), cls, ont, collapsed, label, config$gsea_seed_base, gsea_settings,
        output = file.path(posthoc_root, cls, file_of[[ont]]), return_table = TRUE
      )
    }
    for (base in seed_range_bases) for (ont in c("BP", "CC")) {
      jobs[[length(jobs) + 1L]] <- inferential_gsea_job_spec(
        paste(cls, ont, base, sep = "|"), cls, ont, collapsed, label, base, gsea_settings, return_table = FALSE
      )
    }
  }
  message("Running ", length(jobs), " GSEA jobs in ", min(workers, length(jobs)), " R process(es)...")
  results <- run_inferential_gsea_jobs(jobs, workers, project_root)
  gsea_job_rows <- do.call(rbind, lapply(results, `[[`, "summary"))
  rownames(gsea_job_rows) <- NULL
  gsea_job_rows$primary_seed <- gsea_job_rows$seed_base == config$gsea_seed_base
  gsea_job_rows <- gsea_job_rows[order(match(gsea_job_rows$unit, sample_classes), match(gsea_job_rows$ontology, ontologies),
                                       !gsea_job_rows$primary_seed, gsea_job_rows$seed_base), , drop = FALSE]
  counts <- gsea_job_rows[, c("unit", "ontology", "primary_seed", "seed_base", "seed", "rank_size", "n_terms", "n_fdr05",
                              "n_up", "n_down", "pct_NES_positive", "median_NES", "n_warnings", "runtime_s")]
  names(counts)[1] <- "sample_class"
  write_enrichment_csv(counts, file.path(output_root, "pairing_cno_interaction_gsea_counts.csv"))

  table_of <- function(cls, ont) results[[paste(cls, ont, config$gsea_seed_base, sep = "|")]]$table
  themes <- function(x, positive, n = 3) {
    x <- x[is.finite(x$p.adjust) & x$p.adjust < config$fdr_threshold & (if (positive) x$NES > 0 else x$NES < 0), , drop = FALSE]
    x <- x[order(if (positive) -x$NES else x$NES), , drop = FALSE]
    paste(utils::head(x$Description, n), collapse = "; ")
  }
  n_sig <- function(x) length(significant_terms(x, config$fdr_threshold))
  gsea_summary <- do.call(rbind, lapply(sample_classes, function(cls) {
    bp <- table_of(cls, "BP"); cc <- table_of(cls, "CC"); kg <- table_of(cls, "KEGG")
    data.frame(
      sample_class = cls, go_bp_terms_tested = nrow(bp), go_bp_fdr05 = n_sig(bp),
      top_positive_interaction_themes = themes(bp, TRUE), top_negative_interaction_themes = themes(bp, FALSE),
      kegg_terms_tested = nrow(kg), kegg_fdr05 = n_sig(kg),
      go_cc_terms_tested = nrow(cc), go_cc_fdr05 = n_sig(cc),
      top_positive_cc = themes(cc, TRUE), top_negative_cc = themes(cc, FALSE),
      stringsAsFactors = FALSE
    )
  }))
  write_enrichment_csv(gsea_summary, file.path(output_root, "pairing_cno_interaction_gsea_summary.csv"))

  top_terms <- do.call(rbind, lapply(sample_classes, function(cls) do.call(rbind, lapply(ontologies, function(ont) {
    x <- table_of(cls, ont)
    x <- x[is.finite(x$p.adjust) & x$p.adjust < config$fdr_threshold, , drop = FALSE]
    do.call(rbind, lapply(c("positive", "negative"), function(direction) {
      y <- x[if (direction == "positive") x$NES > 0 else x$NES < 0, , drop = FALSE]
      y <- y[order(y$p.adjust, -abs(y$NES), method = "radix"), , drop = FALSE]
      y <- utils::head(y, 5L)
      if (!nrow(y)) return(NULL)
      data.frame(sample_class = cls, ontology = ont, direction = direction, rank = seq_len(nrow(y)), ID = y$ID,
                 Description = y$Description, setSize = y$setSize, NES = y$NES, p.adjust = y$p.adjust, stringsAsFactors = FALSE)
    }))
  }))))
  write_enrichment_csv(top_terms, file.path(output_root, "pairing_cno_interaction_top_terms.csv"))

  for (ont in ontologies) {
    primary <- counts[counts$primary_seed & counts$ontology == ont, , drop = FALSE]
    key <- c(BP = "go_bp", CC = "go_cc", KEGG = "kegg")[[ont]]
    diagnostics[[paste0(key, "_pct_NES_positive")]] <- primary$pct_NES_positive[match(diagnostics$sample_class, primary$sample_class)]
    diagnostics[[paste0(key, "_median_NES")]] <- primary$median_NES[match(diagnostics$sample_class, primary$sample_class)]
  }
}
invisible(write_enrichment_csv(diagnostics, file.path(output_root, "pairing_cno_interaction_diagnostics.csv")))

# ---- run manifest -------------------------------------------------------------------------
manifest_path <- file.path(output_root, "pairing_cno_interaction_run_manifest.csv")
written <- list.files(output_root, recursive = TRUE, full.names = TRUE)
written <- written[!grepl("^robustness/|^interaction_and_robustness_report[.]txt$",
                          substring(project_normalize_path(written), nchar(output_root) + 2L))]
written <- setdiff(written, manifest_path)
run_manifest <- rbind(
  manifest_rows("run", c("script", "run_timestamp_utc", "output_root", "model", "contrast", "gsea_run"),
                c("04_differential_expression_enrichment/05_pairing_cno_interaction.r", timestamp_utc, output_root,
                  "per sample class: limma ~0 + condition (12 animals, residual df 8); eBayes(robust = TRUE, trend = FALSE); BH within contrast",
                  PAIRING_CNO_INTERACTION_CONTRAST, run_gsea)),
  manifest_rows("input", c("input_gct", "mapped_index", "design_file"), c(input_gct, mapped_index_path, design_file),
                c(project_file_sha256(input_gct), unique(mapped_index$mapped_index_sha256), project_file_sha256(design_file))),
  manifest_rows("input", paste0("mapped:", mapped_index$canonical_comparison), mapped_index$mapped_input_path, mapped_index$mapped_output_sha256),
  manifest_rows("parameter", c(names(gsea_settings), "gsea_seed_base", "seed_range_bases", "gsea_workers", "ranking_statistic", "kegg_gene_sets"),
                c(vapply(gsea_settings, as.character, character(1)), config$gsea_seed_base,
                  if (length(seed_range_bases)) paste(seed_range_bases, collapse = ",") else "none", workers, "t",
                  "downloaded from KEGG at run time (not version-pinned)")),
  inferential_software_rows(),
  if (!is.null(gsea_job_rows)) manifest_rows("gsea_job", gsea_job_rows$job_id, gsea_job_rows$seed,
                                             paste0("label=", gsea_job_rows$seed_label, "; analysis=", gsea_job_rows$seed_analysis,
                                                    "; seed_base=", gsea_job_rows$seed_base, "; terms=", gsea_job_rows$n_terms,
                                                    "; fdr05=", gsea_job_rows$n_fdr05)),
  inferential_output_rows(written, output_root)
)
invisible(write_enrichment_csv(run_manifest, manifest_path))

cat("Pairing x CNO interaction (per-compartment model) written to", output_root, "\n")
print(protein_summary, row.names = FALSE)
print(diagnostics[, c("sample_class", "n_fdr05", "n_fdr10", "min_adj_p", "df_prior_median", "protein_median_t", "plate2_minus_plate1_weight")], row.names = FALSE)
if (!is.null(gsea_job_rows)) print(counts[counts$primary_seed, c("sample_class", "ontology", "n_terms", "n_fdr05", "n_up", "n_down")], row.names = FALSE)
