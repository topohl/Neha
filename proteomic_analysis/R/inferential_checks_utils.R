# Shared helpers for the inferential checks on the per-compartment model:
#   04_differential_expression_enrichment/05_pairing_cno_interaction.r
#   04_differential_expression_enrichment/06_mcherry_learning_robustness.r
#
# Model: the primary model of 01_preprocessing/02b_limma_stat_gct.r (within_compartment mode).
# One limma fit per sample class, ~0 + condition over that class's 12 animals (n = 3 per
# condition, residual df 8), contrasts.fit, eBayes(robust = TRUE, trend = FALSE),
# Benjamini-Hochberg within each contrast.
# GSEA: the gseGO / gseKEGG settings of 04_differential_expression_enrichment/01_clusterProfiler.r,
# ranked by the moderated t, one seed per job from derive_enrichment_seed() applied with
# with_enrichment_seed() and seed = FALSE.

if (!exists("project_normalize_path", mode = "function") || !exists("derive_enrichment_seed", mode = "function")) {
  stop("Source R/project_path_utils.R and R/animal_level_enrichment_utils.R before R/inferential_checks_utils.R", call. = FALSE)
}

# Pairing x CNO interaction: the CNO effect in paired animals minus the CNO effect in unpaired
# animals. Positive = the CNO effect is more positive in paired than in unpaired animals.
PAIRING_CNO_INTERACTION_CONTRAST <- "(paired_cno - paired_veh) - (unpaired_cno - unpaired_veh)"
MCHERRY_LEARNING_CONTRAST <- "paired_veh - unpaired_veh"

# ---- output safety -----------------------------------------------------------------------
# The pooled-run inferential folder and the synthesis folder stay as they are; inputs and the
# enrichment runs are never written by these stages.
inferential_protected_roots <- function(data_root = PROTEOMICS_DATA_ROOT) {
  c(
    file.path(data_root, "01_input"),
    file.path(data_root, "02_data"),
    file.path(data_root, "99_historical"),
    file.path(data_root, "03_output", "enrichment"),
    file.path(data_root, "03_output", "synthesis"),
    file.path(data_root, "03_output", "inferential_checks", "final_inferential_checks_20260826")
  )
}

validate_inferential_output_root <- function(output_root, data_root = PROTEOMICS_DATA_ROOT) {
  project_validate_output_root(
    output_root, inferential_protected_roots(data_root),
    "Inferential-check output root overlaps a protected root: ", name_hit = TRUE
  )
}

inferential_dir_is_empty <- function(path) {
  !dir.exists(path) || !length(list.files(path, all.files = TRUE, no.. = TRUE))
}

inferential_flag <- function(value, name) {
  value <- tolower(trimws(as.character(value)))
  if (!value %in% c("true", "false", "1", "0", "yes", "no")) stop(name, " must be true or false.", call. = FALSE)
  value %in% c("true", "1", "yes")
}

inferential_workers <- function(value, name) {
  workers <- suppressWarnings(as.integer(value))
  if (length(workers) != 1L || is.na(workers) || workers < 1L || workers > 3L) {
    stop(name, " must be 1, 2 or 3.", call. = FALSE)
  }
  workers
}

# ---- input and model ----------------------------------------------------------------------
# Parsed as in 01_preprocessing/02b_limma_stat_gct.r, plus the AnimalID row.
read_animal_level_input_gct <- function(path) {
  lines <- readLines(path, warn = FALSE)
  dims <- scan(text = lines[[2]], what = integer(), quiet = TRUE)
  stopifnot(trimws(lines[[1]]) == "#1.3", length(dims) == 4L, dims[[3]] == 1L)
  n_proteins <- dims[[1]]
  n_meta <- dims[[4]]
  split_line <- function(x) strsplit(x, "\t", fixed = TRUE)[[1]]
  header <- split_line(lines[[3]])
  meta_rows <- lapply(lines[3L + seq_len(n_meta)], split_line)
  names(meta_rows) <- vapply(meta_rows, `[[`, character(1), 1L)
  body <- do.call(rbind, lapply(lines[3L + n_meta + seq_len(n_proteins)], split_line))
  stopifnot(ncol(body) == length(header), nrow(body) == n_proteins, !anyDuplicated(body[, 1]))
  sample_cols <- header[-(1:2)]
  column_index <- seq_along(sample_cols) + 2L
  annotation <- data.frame(
    column = sample_cols,
    sample_class = meta_rows$sample_class[column_index],
    condition = meta_rows$condition[column_index],
    animal = meta_rows$AnimalID[column_index],
    stringsAsFactors = FALSE
  )
  stopifnot(
    all(annotation$sample_class %in% sample_classes),
    all(annotation$condition %in% condition_levels),
    all(table(annotation$sample_class, annotation$condition) == 3L),
    !anyNA(annotation$animal)
  )
  values <- matrix(as.numeric(body[, column_index]), nrow = n_proteins, dimnames = list(body[, 1], sample_cols))
  if (anyNA(values)) stop("The input matrix contains missing values; these checks expect the imputed matrix.", call. = FALSE)
  list(values = values, annotation = annotation, description = stats::setNames(body[, 2], body[, 1]))
}

# One per-compartment fit and one contrast. `values` holds that compartment's columns only.
fit_condition_contrast <- function(values, conditions, contrast) {
  groups <- factor(conditions)
  design <- stats::model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  fit <- limma::lmFit(values, design)
  contrast_matrix <- limma::makeContrasts(contrasts = contrast, levels = design)
  eb <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix), robust = TRUE, trend = FALSE)
  tab <- limma::topTable(eb, number = Inf, sort.by = "none", adjust.method = "BH")
  stopifnot(identical(rownames(tab), rownames(values)))
  list(
    table = tab,
    weights = stats::setNames(contrast_matrix[, 1], rownames(contrast_matrix)),
    df_residual = eb$df.residual[[1]],
    df_prior_median = stats::median(eb$df.prior),
    s2_prior = eb$s2.prior[[1]]
  )
}

# Share of an additive plate difference that a contrast keeps: sum over condition groups of
# contrast weight x fraction of the group's animals on `plate_level`.
contrast_plate_weight <- function(weights, conditions, plates, plate_level = "Plate2") {
  fraction <- tapply(plates == plate_level, conditions, mean)
  if (!all(names(fraction) %in% names(weights))) stop("Plate table has conditions outside the contrast.", call. = FALSE)
  sum(weights[names(fraction)] * fraction)
}

# The pipeline's GSEA reads its statistics from text files holding 15 significant digits (stat GCT,
# mapped CSV). Full precision changes the last bits of the ranking statistic and with them some
# permutation p-values, so statistics fitted in memory are rounded the same way before GSEA.
as_text_precision <- function(x) as.numeric(sprintf("%.15g", x))

# Mapped-table rows with the statistics replaced by another contrast's, ready for
# collapse_enrichment_accessions(). Every mapped protein must be in `tab`.
inferential_mapped_stats <- function(mapped, tab) {
  idx <- match(mapped$original_protein_id, rownames(tab))
  if (anyNA(idx)) stop("The fitted table does not cover every mapped protein; GSEA needs the full input matrix.", call. = FALSE)
  out <- mapped[, c("original_protein_id", "Description", "mapped_gene_symbol", "uniprot_accession", "source_row_id"), drop = FALSE]
  out$log2fc <- as_text_precision(tab$logFC[idx])
  out$t <- as_text_precision(tab$t[idx])
  out$pval <- as_text_precision(tab$P.Value[idx])
  out$padj <- as_text_precision(tab$adj.P.Val[idx])
  out
}

# ---- GSEA ---------------------------------------------------------------------------------
# 01_clusterProfiler.r runs the whole enrichment when sourced, so its KEGG rank builder is taken
# from the parsed script rather than copied: these checks always use the pipeline's Entrez collapse.
load_pipeline_kegg_rank_builder <- function(project_root) {
  path <- file.path(project_root, "04_differential_expression_enrichment", "01_clusterProfiler.r")
  expressions <- as.list(parse(path, keep.source = FALSE))
  is_definition <- vapply(expressions, function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("build_kegg_rank"))
  }, logical(1))
  if (sum(is_definition) != 1L) stop("Expected one build_kegg_rank definition in ", path, call. = FALSE)
  eval(expressions[[which(is_definition)]], envir = globalenv())
  invisible(get("build_kegg_rank", envir = globalenv()))
}

# Prepares a session (the main one or a worker) for inferential_gsea_job(). fgsea would otherwise
# start one worker per core for every job; SerialParam keeps each job in its own process and does
# not change the results (the multilevel step passes one seed to every pathway).
inferential_init_gsea_session <- function(project_root) {
  for (file in c("analysis_labels.R", "project_path_utils.R", "animal_level_enrichment_utils.R", "inferential_checks_utils.R")) {
    source(file.path(project_root, "R", file))
  }
  load_pipeline_kegg_rank_builder(project_root)
  for (pkg in c("clusterProfiler", "org.Mm.eg.db", "fgsea", "BiocParallel", "withr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("GSEA requires package ", pkg, call. = FALSE)
  }
  BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  invisible(TRUE)
}

inferential_gsea_job <- function(job) {
  started <- Sys.time()
  s <- job$settings
  rank <- if (identical(job$ontology, "KEGG")) {
    build_kegg_rank(job$collapsed, "t", "canonical")$rank
  } else {
    build_gsea_rank(job$collapsed, "t", "canonical")$rank
  }
  call_gsea <- function() {
    if (identical(job$ontology, "KEGG")) {
      clusterProfiler::gseKEGG(
        geneList = rank, organism = "mmu", keyType = "ncbi-geneid", exponent = 1,
        minGSSize = s$min_gs_size, maxGSSize = s$max_gs_size, eps = 1e-10,
        pvalueCutoff = s$pvalue_cutoff, pAdjustMethod = s$p_adjust_method,
        verbose = FALSE, seed = FALSE, by = "fgsea"
      )
    } else {
      clusterProfiler::gseGO(
        geneList = rank, ont = job$ontology, keyType = "UNIPROT", OrgDb = org.Mm.eg.db::org.Mm.eg.db,
        exponent = 1, minGSSize = s$min_gs_size, maxGSSize = s$max_gs_size, eps = 1e-10,
        pvalueCutoff = s$pvalue_cutoff, pAdjustMethod = s$p_adjust_method,
        verbose = FALSE, seed = FALSE, by = "fgsea"
      )
    }
  }
  warnings <- character()
  result <- withCallingHandlers(
    with_enrichment_seed(job$seed, call_gsea()),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  table <- enrichment_result_table(result)
  if (!is.na(job$output)) write_enrichment_csv(table, job$output)
  significant <- is.finite(table$p.adjust) & table$p.adjust < s$fdr_threshold
  list(
    id = job$id,
    table = if (isTRUE(job$return_table)) table else NULL,
    summary = data.frame(
      job_id = job$id, unit = job$unit, ontology = job$ontology, seed_label = job$label,
      seed_analysis = job$analysis, seed_base = job$seed_base, seed = job$seed,
      rank_size = length(rank), n_terms = nrow(table), n_fdr05 = sum(significant),
      n_up = sum(significant & table$NES > 0), n_down = sum(significant & table$NES < 0),
      pct_NES_positive = if (nrow(table)) 100 * mean(table$NES > 0) else NA_real_,
      median_NES = if (nrow(table)) stats::median(table$NES) else NA_real_,
      n_warnings = length(unique(warnings)),
      runtime_s = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
      output = if (is.na(job$output)) NA_character_ else job$output,
      stringsAsFactors = FALSE
    )
  )
}

inferential_gsea_job_spec <- function(id, unit, ontology, collapsed, label, seed_base, settings,
                                      output = NA_character_, return_table = TRUE) {
  analysis <- if (identical(ontology, "KEGG")) "KEGG_GSEA" else paste0("GO_", ontology, "_GSEA")
  list(
    id = id, unit = unit, ontology = ontology, collapsed = collapsed, label = label,
    analysis = analysis, seed_base = as.integer(seed_base),
    seed = derive_enrichment_seed(seed_base, label, analysis),
    settings = settings, output = output, return_table = return_table
  )
}

# Defined at top level so that only the job itself is sent to a worker.
inferential_safe_gsea_job <- function(job) {
  tryCatch(inferential_gsea_job(job), error = function(e) list(id = job$id, error = conditionMessage(e)))
}

# Runs the jobs in `workers` R processes (1 = this session). Load-balanced, longest jobs first.
run_inferential_gsea_jobs <- function(jobs, workers, project_root) {
  order_key <- match(vapply(jobs, `[[`, character(1), "ontology"), c("BP", "CC", "KEGG"))
  jobs <- jobs[order(order_key)]
  workers <- min(workers, length(jobs))
  if (workers <= 1L) {
    inferential_init_gsea_session(project_root)
    results <- lapply(jobs, inferential_safe_gsea_job)
  } else {
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, inferential_init_gsea_session, project_root)
    results <- parallel::clusterApplyLB(cl, jobs, inferential_safe_gsea_job)
  }
  failed <- vapply(results, function(r) !is.null(r$error), logical(1))
  if (any(failed)) {
    stop("GSEA failed for: ", paste(vapply(results[failed], function(r) paste0(r$id, " (", r$error, ")"), character(1)), collapse = "; "), call. = FALSE)
  }
  stats::setNames(results, vapply(results, `[[`, character(1), "id"))
}

significant_terms <- function(table, threshold = 0.05) {
  table$ID[is.finite(table$p.adjust) & table$p.adjust < threshold]
}

# ---- provenance ---------------------------------------------------------------------------
manifest_rows <- function(section, name, value, detail = NA_character_) {
  data.frame(section = section, name = name, value = as.character(value), detail = as.character(detail), stringsAsFactors = FALSE)
}

inferential_software_rows <- function() {
  versions <- enrichment_package_versions()
  extra <- c("limma", "BiocParallel")
  extra_versions <- vapply(extra, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }, character(1))
  manifest_rows("software", c(versions$component, extra), c(versions$version, extra_versions))
}

inferential_output_rows <- function(paths, root) {
  paths <- sort(paths)
  relative <- substring(project_normalize_path(paths), nchar(project_normalize_path(root)) + 2L)
  manifest_rows("output", relative, vapply(paths, project_file_sha256, character(1)))
}
