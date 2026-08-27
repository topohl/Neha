# Canonical software contracts for the Neha animal-level enrichment branch.
# Biological comparison definitions remain owned by R/analysis_labels.R.

if (!exists("neha_normalize_path", mode = "function")) {
  stop("R/neha_path_utils.R must be sourced before R/animal_level_enrichment_utils.R", call. = FALSE)
}

# Retained names delegate to the shared primitives in R/neha_path_utils.R.
# Behaviour change on non-Windows only: comparison is now case-insensitive everywhere
# (previously case-sensitive off Windows). See the note in R/neha_path_utils.R.
neha_enrichment_normalize_path <- function(path, must_work = FALSE) neha_normalize_path(path, must_work)
neha_enrichment_path_is_within <- function(path, parent) neha_path_is_within(path, parent)
neha_enrichment_paths_overlap <- function(left, right) neha_paths_overlap(left, right)

neha_enrichment_sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for enrichment provenance.", call. = FALSE)
  }
  path <- neha_enrichment_normalize_path(path, must_work = TRUE)
  tolower(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

neha_enrichment_default_paths <- function() {
  # Paths reflect the 2026-08-26 input/data/output restructure; see CANONICAL_OUTPUTS.md.
  project_root <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
  animal_root <- file.path(project_root, "02_data", "animal_level")
  historical_root <- file.path(project_root, "99_historical")
  list(
    project_root = project_root,
    dataset_root = file.path(project_root, "02_data"),
    mapped_root = file.path(animal_root, "mapped"),
    mapped_index = file.path(animal_root, "mapped", "indexMappedComparisons.csv"),
    output_root = file.path(project_root, "03_output", "enrichment"),
    split_root = file.path(animal_root, "split"),
    historical_input_roots = file.path(historical_root, c("datasets_raw", "datasets_mapped")),
    historical_output_roots = c(
      file.path(historical_root, "core_enrichment"),
      historical_root,
      file.path(historical_root, "plots_pairwise")
    )
  )
}

validate_neha_enrichment_mapped_root <- function(mapped_root, historical_input_roots) {
  mapped_root <- neha_enrichment_normalize_path(mapped_root)
  protected <- vapply(historical_input_roots, neha_enrichment_normalize_path, character(1))
  hit <- vapply(protected, function(path) neha_enrichment_paths_overlap(mapped_root, path), logical(1))
  if (any(hit)) {
    stop(
      "Canonical animal-level enrichment cannot read a historical mapped root: ",
      mapped_root,
      call. = FALSE
    )
  }
  mapped_root
}

validate_neha_enrichment_output_root <- function(
    output_root,
    mapped_root,
    split_root,
    historical_input_roots,
    historical_output_roots) {
  output_root <- neha_enrichment_normalize_path(output_root)
  protected <- c(mapped_root, split_root, historical_input_roots, historical_output_roots)
  protected <- vapply(protected, neha_enrichment_normalize_path, character(1))
  hit <- vapply(protected, function(path) neha_enrichment_paths_overlap(output_root, path), logical(1))
  if (any(hit)) {
    stop(
      "Canonical animal-level enrichment output root overlaps a protected input or historical root: ",
      output_root, " (protected: ", protected[[which(hit)[[1]]]], ")",
      call. = FALSE
    )
  }
  output_root
}

neha_enrichment_env_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = if (default) "true" else "false")))
  if (!value %in% c("true", "false", "1", "0", "yes", "no")) {
    stop(name, " must be true or false.", call. = FALSE)
  }
  value %in% c("true", "1", "yes")
}

neha_enrichment_env_number <- function(name, default, lower = -Inf, upper = Inf, integer = FALSE) {
  raw <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.numeric(raw))
  if (length(value) != 1L || is.na(value) || !is.finite(value) || value < lower || value > upper) {
    stop(name, " must be one finite value between ", lower, " and ", upper, ".", call. = FALSE)
  }
  if (isTRUE(integer) && value != as.integer(value)) stop(name, " must be a whole number.", call. = FALSE)
  if (isTRUE(integer)) as.integer(value) else value
}

neha_enrichment_env_vector <- function(name) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(character(0))
  out <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  unique(out[nzchar(out)])
}

neha_enrichment_env_choice <- function(name, choices, default) {
  value <- tolower(trimws(Sys.getenv(name, unset = default)))
  if (length(value) != 1L || !value %in% choices) {
    stop(name, " must be one of: ", paste(choices, collapse = ", "), ".", call. = FALSE)
  }
  value
}

resolve_neha_enrichment_config <- function(
    mapped_root = Sys.getenv("NEHA_ENRICHMENT_MAPPED_ROOT", unset = ""),
    mapped_index = Sys.getenv("NEHA_ENRICHMENT_MAPPED_INDEX", unset = ""),
    output_root = Sys.getenv("NEHA_ENRICHMENT_OUTPUT_ROOT", unset = "")) {
  defaults <- neha_enrichment_default_paths()
  choose <- function(value, default) if (nzchar(trimws(value))) value else default
  mapped_root <- validate_neha_enrichment_mapped_root(
    choose(mapped_root, defaults$mapped_root),
    defaults$historical_input_roots
  )
  mapped_index <- neha_enrichment_normalize_path(
    choose(mapped_index, file.path(mapped_root, "indexMappedComparisons.csv"))
  )
  if (!neha_enrichment_path_is_within(mapped_index, mapped_root)) {
    stop("Canonical mapped index must resolve inside the configured mapped root.", call. = FALSE)
  }
  output_root <- validate_neha_enrichment_output_root(
    choose(output_root, defaults$output_root),
    mapped_root,
    defaults$split_root,
    defaults$historical_input_roots,
    defaults$historical_output_roots
  )

  ontology <- toupper(trimws(Sys.getenv("NEHA_ENRICHMENT_ONTOLOGY", unset = "BP")))
  if (!ontology %in% c("BP", "CC", "MF")) {
    stop("NEHA_ENRICHMENT_ONTOLOGY must be BP, CC, or MF.", call. = FALSE)
  }
  min_gs_size <- neha_enrichment_env_number("NEHA_ENRICHMENT_MIN_GS_SIZE", 10, 1, Inf, integer = TRUE)
  max_gs_size <- neha_enrichment_env_number("NEHA_ENRICHMENT_MAX_GS_SIZE", 800, 1, Inf, integer = TRUE)
  if (min_gs_size > max_gs_size) stop("Minimum gene-set size exceeds maximum gene-set size.", call. = FALSE)
  gsea_rank <- neha_enrichment_env_choice("NEHA_ENRICHMENT_GSEA_RANK", c("t", "log2fc"), "t")

  list(
    mapped_root = mapped_root,
    mapped_index = mapped_index,
    output_root = output_root,
    project_root = neha_enrichment_normalize_path(defaults$project_root),
    dataset_root = neha_enrichment_normalize_path(defaults$dataset_root),
    ontology = ontology,
    pvalue_cutoff = neha_enrichment_env_number("NEHA_ENRICHMENT_PVALUE_CUTOFF", 1, 0, 1),
    qvalue_cutoff = neha_enrichment_env_number("NEHA_ENRICHMENT_QVALUE_CUTOFF", 1, 0, 1),
    fdr_threshold = neha_enrichment_env_number("NEHA_ENRICHMENT_FDR_THRESHOLD", 0.05, 0, 1),
    top_abs_log2fc = neha_enrichment_env_number("NEHA_ENRICHMENT_TOP_ABS_LOG2FC", 1, 0, Inf),
    p_adjust_method = "BH",
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size,
    simplify = neha_enrichment_env_flag("NEHA_ENRICHMENT_SIMPLIFY", FALSE),
    simplify_cutoff = neha_enrichment_env_number("NEHA_ENRICHMENT_SIMPLIFY_CUTOFF", 0.7, 0, 1),
    gsea_rank = gsea_rank,
    gsea_sensitivity_rank = if (identical(gsea_rank, "t")) "log2fc" else NA_character_,
    gsea_seed_base = neha_enrichment_env_number("NEHA_ENRICHMENT_GSEA_SEED_BASE", 20260824, 1, 2147483646, integer = TRUE),
    kegg_enabled = neha_enrichment_env_flag("NEHA_ENRICHMENT_KEGG", TRUE),
    plots_enabled = neha_enrichment_env_flag("NEHA_ENRICHMENT_PLOTS", TRUE),
    force = neha_enrichment_env_flag("NEHA_ENRICHMENT_FORCE", FALSE),
    nk3r_genes = neha_enrichment_env_vector("NEHA_ENRICHMENT_NK3R_GENES"),
    selected_uniprot = neha_enrichment_env_vector("NEHA_ENRICHMENT_SELECTED_UNIPROT"),
    path_ids = neha_enrichment_env_vector("NEHA_ENRICHMENT_PATH_IDS"),
    historical_input_roots = vapply(defaults$historical_input_roots, neha_enrichment_normalize_path, character(1)),
    historical_output_roots = vapply(defaults$historical_output_roots, neha_enrichment_normalize_path, character(1))
  )
}

neha_enrichment_required_mapped_index_columns <- function() {
  c(
    "canonical_comparison", "canonical_contrast", "sample_class",
    "numerator_condition", "denominator_condition", "historical_comparison_alias",
    "mapping_direction", "output_mapped_path", "n_input_proteins", "n_output_mapped_rows",
    "n_mapped", "n_unmapped", "mapped_output_sha256", "source_split_index_path",
    "source_split_index_sha256", "source_split_sha256", "source_gct_sha256",
    "mapping_reference_path", "mapping_reference_version",
    "mapping_reference_snapshot_date_utc", "mapping_reference_sha256", "row_accounting_valid"
  )
}

read_neha_enrichment_mapped_index <- function(
    mapped_index,
    mapped_root,
    expected_manifest = neha_primary_contrast_manifest(),
    verify_hashes = TRUE) {
  mapped_root <- validate_neha_enrichment_mapped_root(
    mapped_root,
    neha_enrichment_default_paths()$historical_input_roots
  )
  mapped_index <- neha_enrichment_normalize_path(mapped_index, must_work = TRUE)
  if (!neha_enrichment_path_is_within(mapped_index, mapped_root)) {
    stop("Canonical mapped index resolves outside the configured mapped root.", call. = FALSE)
  }
  index <- utils::read.csv(mapped_index, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(neha_enrichment_required_mapped_index_columns(), names(index))
  if (length(missing)) stop("Mapped index is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (nrow(index) != 12L || nrow(index) != nrow(expected_manifest)) {
    stop("Canonical animal-level enrichment requires exactly 12 mapped comparisons; found ", nrow(index), ".", call. = FALSE)
  }
  if (any(as.character(index$mapping_direction) != "forward")) {
    stop("Canonical animal-level enrichment requires mapping_direction = forward for every comparison.", call. = FALSE)
  }
  duplicate_fields <- c("canonical_comparison", "canonical_contrast", "output_mapped_path")
  duplicate_fields <- duplicate_fields[vapply(duplicate_fields, function(field) anyDuplicated(index[[field]]) > 0L, logical(1))]
  if (length(duplicate_fields)) stop("Mapped index contains duplicates in: ", paste(duplicate_fields, collapse = ", "), call. = FALSE)

  expected_keys <- expected_manifest$canonical_comparison
  if (!setequal(as.character(index$canonical_comparison), expected_keys)) {
    stop("Mapped index comparison coverage does not match the shared Neha manifest.", call. = FALSE)
  }
  index <- index[match(expected_keys, index$canonical_comparison), , drop = FALSE]
  expected <- data.frame(
    canonical_contrast = expected_manifest$canonical_contrast,
    sample_class = expected_manifest$sample_class,
    numerator_condition = expected_manifest$case_condition,
    denominator_condition = expected_manifest$reference_condition,
    historical_comparison_alias = expected_manifest$historical_comparison_name,
    stringsAsFactors = FALSE
  )
  for (field in names(expected)) {
    if (!identical(as.character(index[[field]]), expected[[field]])) {
      stop("Mapped index ", field, " does not match the shared Neha manifest.", call. = FALSE)
    }
  }
  if (any(!as.logical(index$row_accounting_valid))) stop("Mapped index contains failed upstream row accounting.", call. = FALSE)
  paths <- as.character(index$output_mapped_path)
  if (any(!file.exists(paths))) stop("Mapped index references missing mapped files.", call. = FALSE)
  paths <- vapply(paths, neha_enrichment_normalize_path, character(1), must_work = TRUE)
  forward_root <- file.path(mapped_root, "forward")
  outside <- !vapply(paths, neha_enrichment_path_is_within, logical(1), parent = forward_root)
  if (any(outside)) {
    stop("Indexed mapped path resolves outside canonical mapped/forward: ", paths[[which(outside)[[1]]]], call. = FALSE)
  }
  index$mapped_input_path <- paths
  if (isTRUE(verify_hashes)) {
    observed <- vapply(paths, neha_enrichment_sha256, character(1))
    if (!identical(unname(observed), unname(tolower(as.character(index$mapped_output_sha256))))) {
      stop("Mapped input SHA-256 does not match indexMappedComparisons.csv.", call. = FALSE)
    }
  }
  index$mapped_index_path <- mapped_index
  index$mapped_index_sha256 <- neha_enrichment_sha256(mapped_index)
  rownames(index) <- NULL
  index
}

neha_enrichment_required_mapped_columns <- function() {
  c(
    "gene_symbol", "original_protein_id", "Description", "mapped_gene_symbol",
    "uniprot_accession", "mapping_status", "mapping_strategy", "source_row_id",
    "multi_protein", "log2fc", "aveExpr", "t", "pval", "padj", "B",
    "significant", "logpval"
  )
}

read_neha_enrichment_mapped_file <- function(path, expected_rows = NULL, expected_sha256 = NULL) {
  path <- neha_enrichment_normalize_path(path, must_work = TRUE)
  if (!is.null(expected_sha256) && !identical(neha_enrichment_sha256(path), tolower(as.character(expected_sha256)))) {
    stop("Mapped file SHA-256 mismatch: ", path, call. = FALSE)
  }
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(neha_enrichment_required_mapped_columns(), names(x))
  if (length(missing)) stop("Mapped enrichment input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!is.null(expected_rows) && nrow(x) != as.integer(expected_rows)) {
    stop("Mapped enrichment input row count mismatch: expected ", expected_rows, ", found ", nrow(x), ".", call. = FALSE)
  }
  accessions <- trimws(as.character(x$uniprot_accession))
  if (anyNA(accessions) || any(!nzchar(accessions))) stop("Mapped enrichment input has missing UniProt accessions.", call. = FALSE)
  if (any(as.character(x$mapping_status) != "mapped")) stop("Mapped enrichment input contains non-mapped rows.", call. = FALSE)
  if (!identical(as.character(x$gene_symbol), accessions)) {
    stop("Mapped compatibility gene_symbol is not aligned with uniprot_accession.", call. = FALSE)
  }
  source_rows <- suppressWarnings(as.integer(x$source_row_id))
  if (anyNA(source_rows) || anyDuplicated(source_rows)) stop("Mapped enrichment input source_row_id must be complete and unique.", call. = FALSE)
  numeric_fields <- c("log2fc", "aveExpr", "t", "pval", "padj", "B", "logpval")
  for (field in numeric_fields) x[[field]] <- suppressWarnings(as.numeric(x[[field]]))
  x$source_row_id <- source_rows
  x$uniprot_accession <- accessions
  x
}

collapse_neha_enrichment_accessions <- function(x, ranking_statistic = "log2fc") {
  if (!ranking_statistic %in% names(x)) stop("Ranking statistic is missing: ", ranking_statistic, call. = FALSE)
  rank_value <- suppressWarnings(as.numeric(x[[ranking_statistic]]))
  accession <- trimws(as.character(x$uniprot_accession))
  source_row <- suppressWarnings(as.integer(x$source_row_id))
  if (anyNA(accession) || any(!nzchar(accession)) || anyNA(source_row)) {
    stop("UniProt accession and source_row_id must be complete before duplicate collapse.", call. = FALSE)
  }
  finite <- is.finite(rank_value)
  abs_rank <- ifelse(finite, abs(rank_value), -Inf)
  ordering <- order(
    accession,
    -as.integer(finite),
    -abs_rank,
    source_row,
    as.character(x$original_protein_id),
    method = "radix"
  )
  selected <- !duplicated(accession[ordering])
  selected_rows <- ordering[selected]
  collapsed <- x[selected_rows, , drop = FALSE]
  collapsed <- collapsed[order(collapsed$source_row_id, method = "radix"), , drop = FALSE]
  rownames(collapsed) <- NULL

  duplicate_accessions <- unique(accession[duplicated(accession) | duplicated(accession, fromLast = TRUE)])
  audit_rows <- which(accession %in% duplicate_accessions)
  duplicate_audit <- x[audit_rows, , drop = FALSE]
  if (nrow(duplicate_audit)) {
    selected_source <- stats::setNames(
      source_row[selected_rows],
      accession[selected_rows]
    )
    duplicate_audit$ranking_statistic <- ranking_statistic
    duplicate_audit$ranking_value <- rank_value[audit_rows]
    duplicate_audit$absolute_ranking_value <- abs(rank_value[audit_rows])
    duplicate_audit$selected_representative <- source_row[audit_rows] == unname(selected_source[accession[audit_rows]])
    duplicate_audit$selected_source_row_id <- as.integer(unname(selected_source[accession[audit_rows]]))
    duplicate_audit$selection_rule <- paste0(
      "largest_absolute_", ranking_statistic,
      ";finite_preferred;ties_by_source_row_id_then_original_protein_id"
    )
    duplicate_audit <- duplicate_audit[order(
      duplicate_audit$uniprot_accession,
      -as.integer(duplicate_audit$selected_representative),
      duplicate_audit$source_row_id,
      method = "radix"
    ), , drop = FALSE]
    rownames(duplicate_audit) <- NULL
  } else {
    duplicate_audit <- data.frame(
      x[0, , drop = FALSE],
      ranking_statistic = character(), ranking_value = numeric(),
      absolute_ranking_value = numeric(), selected_representative = logical(),
      selected_source_row_id = integer(), selection_rule = character(),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  if (anyDuplicated(collapsed$uniprot_accession)) stop("Duplicate UniProt collapse failed.", call. = FALSE)
  list(
    collapsed = collapsed,
    duplicate_audit = duplicate_audit,
    n_source_rows = nrow(x),
    n_unique_accessions = nrow(collapsed),
    n_duplicate_rows_collapsed = nrow(x) - nrow(collapsed),
    n_duplicated_accessions = length(duplicate_accessions),
    rule = paste0(
      "largest_absolute_", ranking_statistic,
      ";finite_preferred;ties_by_source_row_id_then_original_protein_id"
    )
  )
}

neha_gsea_tie_diagnostics <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  finite <- values[is.finite(values)]
  n_finite <- length(finite)
  counts <- if (n_finite) table(finite, useNA = "no") else integer()
  tied_counts <- counts[counts > 1L]
  n_unique <- length(counts)
  rows_participating_in_ties <- if (length(tied_counts)) sum(tied_counts) else 0L
  data.frame(
    n_finite = as.integer(n_finite),
    n_unique = as.integer(n_unique),
    redundancy_fraction = if (n_finite) (n_finite - n_unique) / n_finite else NA_real_,
    rows_participating_in_ties = as.integer(rows_participating_in_ties),
    tied_row_fraction = if (n_finite) rows_participating_in_ties / n_finite else NA_real_,
    largest_tie = as.integer(if (length(counts)) max(counts) else 0L),
    stringsAsFactors = FALSE
  )
}

build_neha_gsea_rank <- function(collapsed, ranking_statistic = "t", analysis_role = "canonical") {
  if (!ranking_statistic %in% c("t", "log2fc")) {
    stop("GSEA ranking statistic must be 't' or 'log2fc'.", call. = FALSE)
  }
  if (!analysis_role %in% c("canonical", "sensitivity")) {
    stop("GSEA analysis role must be canonical or sensitivity.", call. = FALSE)
  }
  values <- suppressWarnings(as.numeric(collapsed[[ranking_statistic]]))
  ids <- as.character(collapsed$uniprot_accession)
  included <- is.finite(values) & !is.na(ids) & nzchar(ids)
  if (identical(ranking_statistic, "t") && any(!included)) {
    stop("Moderated t GSEA rank must contain one finite value for every selected UniProt row.", call. = FALSE)
  }
  tie_diagnostics <- neha_gsea_tie_diagnostics(values[included])
  rank_audit <- data.frame(
    uniprot_accession = ids,
    original_protein_id = as.character(collapsed$original_protein_id),
    source_row_id = as.integer(collapsed$source_row_id),
    ranking_statistic = ranking_statistic,
    rank_source_column = ranking_statistic,
    analysis_role = analysis_role,
    ranking_value = values,
    included_in_gsea = included,
    exclusion_reason = ifelse(included, NA_character_, "non_finite_ranking_statistic"),
    rank_direction = "positive_is_higher_in_canonical_numerator",
    stringsAsFactors = FALSE
  )
  values <- values[included]
  ids <- ids[included]
  ordering <- order(-values, ids, method = "radix")
  ranks <- values[ordering]
  names(ranks) <- ids[ordering]
  if (anyDuplicated(names(ranks))) stop("GSEA rank contains duplicate UniProt identifiers.", call. = FALSE)
  if (any(!is.finite(ranks))) stop("GSEA rank contains NA or infinite values.", call. = FALSE)
  if (length(ranks) && !identical(ranks, sort(ranks, decreasing = TRUE))) {
    stop("GSEA rank is not deterministically decreasing.", call. = FALSE)
  }
  list(
    rank = ranks,
    audit = rank_audit,
    statistic = ranking_statistic,
    rank_source_column = ranking_statistic,
    analysis_role = analysis_role,
    tie_diagnostics = tie_diagnostics
  )
}

build_neha_ora_sets <- function(collapsed, fdr_threshold = 0.05, top_abs_log2fc = 1) {
  ids <- as.character(collapsed$uniprot_accession)
  log2fc <- suppressWarnings(as.numeric(collapsed$log2fc))
  padj <- suppressWarnings(as.numeric(collapsed$padj))
  universe <- sort(unique(ids[!is.na(ids) & nzchar(ids)]), method = "radix")
  significant <- is.finite(padj) & padj < fdr_threshold & is.finite(log2fc)
  make_set <- function(mask) sort(unique(ids[mask & ids %in% universe]), method = "radix")
  sets <- list(
    universe = universe,
    all_significant = make_set(significant),
    up_significant = make_set(significant & log2fc > 0),
    down_significant = make_set(significant & log2fc < 0),
    top_abs_log2fc = make_set(is.finite(log2fc) & abs(log2fc) > top_abs_log2fc)
  )
  if (any(vapply(sets[-1], anyDuplicated, integer(1)) > 0L)) stop("ORA foreground contains duplicate accessions.", call. = FALSE)
  if (any(!unlist(sets[-1], use.names = FALSE) %in% universe)) stop("ORA foreground escaped the measured universe.", call. = FALSE)
  sets$fdr_threshold <- fdr_threshold
  sets$top_abs_log2fc_threshold <- top_abs_log2fc
  sets
}

derive_neha_enrichment_seed <- function(base_seed, comparison, analysis_type) {
  modulus <- 2147483646
  key <- paste0(nchar(comparison, type = "bytes"), ":", comparison, "|", analysis_type)
  hash <- as.double(as.integer(base_seed)) %% modulus
  for (code in utf8ToInt(enc2utf8(key))) hash <- (hash * 131 + as.double(code)) %% modulus
  seed <- as.integer(hash)
  if (seed < 1L) seed <- 1L
  seed
}

with_neha_enrichment_seed <- function(seed, expression) {
  if (!requireNamespace("withr", quietly = TRUE)) stop("Package 'withr' is required for deterministic GSEA.", call. = FALSE)
  withr::with_preserve_seed({
    old_kind <- RNGkind()
    tryCatch({
      RNGkind(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
      set.seed(as.integer(seed))
      force(expression)
    }, finally = {
      do.call(RNGkind, as.list(old_kind))
    })
  })
}

neha_enrichment_empty_result <- function() {
  data.frame(
    ID = character(), Description = character(), setSize = integer(),
    enrichmentScore = numeric(), NES = numeric(), pvalue = numeric(),
    p.adjust = numeric(), qvalue = numeric(), rank = integer(),
    leading_edge = character(), core_enrichment = character(),
    stringsAsFactors = FALSE
  )
}

neha_enrichment_result_table <- function(result) {
  if (is.null(result)) return(neha_enrichment_empty_result())
  out <- tryCatch(as.data.frame(result), error = function(e) NULL)
  if (is.null(out)) return(neha_enrichment_empty_result())
  if (!ncol(out) && !nrow(out)) return(neha_enrichment_empty_result())
  out
}

write_neha_enrichment_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop("Failed to verify enrichment CSV: ", path, call. = FALSE)
  }
  neha_enrichment_normalize_path(path, must_work = TRUE)
}

neha_enrichment_package_versions <- function() {
  packages <- c(
    "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "DOSE", "enrichplot",
    "GOSemSim", "GO.db", "KEGGREST", "fgsea", "ggplot2", "digest", "withr"
  )
  data.frame(
    component = c("R", packages),
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      vapply(packages, function(pkg) {
        if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
      }, character(1))
    ),
    available = c(TRUE, vapply(packages, requireNamespace, logical(1), quietly = TRUE)),
    stringsAsFactors = FALSE
  )
}

validate_neha_enrichment_index <- function(index, require_files = TRUE) {
  required <- c(
    "canonical_comparison", "canonical_contrast", "sample_class",
    "numerator_condition", "denominator_condition", "historical_comparison_alias",
    "mapping_direction", "mapped_input_path", "mapped_input_sha256",
    "execution_status", "ranking_statistic", "rank_source_column", "gsea_analysis_role",
    "rank_vector_size", "rank_n_finite", "rank_n_unique", "rank_redundancy_fraction",
    "rank_rows_participating_in_ties", "rank_tied_row_fraction", "rank_largest_tie",
    "ora_universe_size", "gsea_rank_audit", "go_gsea_output", "comparison_manifest_path"
  )
  missing <- setdiff(required, names(index))
  if (length(missing)) stop("Enrichment index is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (nrow(index) != 12L || anyDuplicated(index$canonical_comparison) || anyDuplicated(index$mapped_input_path)) {
    stop("Enrichment index must contain exactly 12 unique comparisons and input paths.", call. = FALSE)
  }
  if (any(index$mapping_direction != "forward")) stop("Enrichment index contains a non-forward comparison.", call. = FALSE)
  if (any(!index$ranking_statistic %in% c("t", "log2fc")) ||
      any(index$rank_source_column != index$ranking_statistic) ||
      any(index$gsea_analysis_role != "canonical")) {
    stop("Enrichment index has an invalid canonical GSEA rank contract.", call. = FALSE)
  }
  t_primary <- index$ranking_statistic == "t"
  if (any(t_primary)) {
    sensitivity_required <- c(
      "sensitivity_ranking_statistic", "sensitivity_rank_source_column",
      "sensitivity_gsea_analysis_role", "sensitivity_rank_vector_size",
      "sensitivity_gsea_rank_audit", "sensitivity_go_gsea_output"
    )
    missing_sensitivity <- setdiff(sensitivity_required, names(index))
    if (length(missing_sensitivity)) {
      stop("t-primary enrichment index is missing log2fc sensitivity fields: ",
           paste(missing_sensitivity, collapse = ", "), call. = FALSE)
    }
    if (any(index$sensitivity_ranking_statistic[t_primary] != "log2fc") ||
        any(index$sensitivity_rank_source_column[t_primary] != "log2fc") ||
        any(index$sensitivity_gsea_analysis_role[t_primary] != "sensitivity")) {
      stop("t-primary enrichment index has an invalid log2fc sensitivity contract.", call. = FALSE)
    }
  }
  expected <- neha_primary_contrast_manifest()
  index_ordered <- index[match(expected$canonical_comparison, index$canonical_comparison), , drop = FALSE]
  contract <- list(
    canonical_comparison = expected$canonical_comparison,
    canonical_contrast = expected$canonical_contrast,
    sample_class = expected$sample_class,
    numerator_condition = expected$case_condition,
    denominator_condition = expected$reference_condition,
    historical_comparison_alias = expected$historical_comparison_name
  )
  for (field in names(contract)) {
    if (!identical(as.character(index_ordered[[field]]), as.character(contract[[field]]))) {
      stop("Enrichment index ", field, " does not match the shared Neha contrast manifest.", call. = FALSE)
    }
  }
  if (isTRUE(require_files)) {
    success <- index$execution_status == "success"
    file_fields <- c("mapped_input_path", "gsea_rank_audit", "go_gsea_output", "comparison_manifest_path")
    if (any(t_primary & success)) {
      file_fields <- c(file_fields, "sensitivity_gsea_rank_audit", "sensitivity_go_gsea_output")
    }
    for (field in file_fields) {
      paths <- as.character(index[[field]][success])
      if (any(!file.exists(paths))) stop("Successful enrichment row references missing ", field, ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}
