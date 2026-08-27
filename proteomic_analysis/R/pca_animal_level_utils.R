# Focused input, preprocessing, and audit contracts for Neha animal-level PCA.

if (!exists("neha_normalize_path", mode = "function")) {
  stop("R/neha_path_utils.R must be sourced before R/pca_animal_level_utils.R", call. = FALSE)
}

# Retained names delegate to the shared primitives in R/neha_path_utils.R.
neha_pca_normalize_path <- function(path, must_work = FALSE) neha_normalize_path(path, must_work)
neha_pca_path_is_within <- function(path, root) neha_path_is_within(path, root)

validate_neha_pca_output_root <- function(output_root, historical_root) {
  # Message text is asserted by tests/test_pca_animal_level.R; keep it verbatim. Note the
  # historical root (not the output root) is named, which is why this does not use the
  # shared neha_validate_output_root() helper.
  output_root <- neha_pca_normalize_path(output_root)
  historical_root <- neha_pca_normalize_path(historical_root)
  if (neha_pca_path_is_within(output_root, historical_root)) {
    stop(
      "Animal-level PCA output cannot overwrite the historical PCA root: ",
      historical_root,
      call. = FALSE
    )
  }
  output_root
}

validate_neha_pca_animal_input <- function(parsed_gct, expected_n = 3L) {
  if (!is.list(parsed_gct) || is.null(parsed_gct$matrix) || is.null(parsed_gct$column_metadata)) {
    stop("Animal-level PCA input must be a parsed ProTigy GCT.", call. = FALSE)
  }

  expression_matrix <- as.matrix(parsed_gct$matrix)
  storage.mode(expression_matrix) <- "numeric"
  if (nrow(expression_matrix) < 2L || ncol(expression_matrix) < 2L) {
    stop("Animal-level PCA expression matrix must contain at least two proteins and samples.", call. = FALSE)
  }
  if (is.null(rownames(expression_matrix)) || anyDuplicated(rownames(expression_matrix)) ||
      any(is.na(rownames(expression_matrix)) | !nzchar(rownames(expression_matrix)))) {
    stop("Animal-level PCA protein identifiers must be complete and unique.", call. = FALSE)
  }
  if (is.null(colnames(expression_matrix)) || anyDuplicated(colnames(expression_matrix)) ||
      any(is.na(colnames(expression_matrix)) | !nzchar(colnames(expression_matrix)))) {
    stop("Animal-level PCA sample identifiers must be complete and unique.", call. = FALSE)
  }
  if (any(!is.finite(expression_matrix))) {
    stop("Animal-level PCA abundance matrix contains nonfinite values; no imputation or missingness filtering is permitted.", call. = FALSE)
  }

  column_metadata <- as.matrix(parsed_gct$column_metadata)
  required_metadata <- c("AnimalID", "condition_code", "condition", "sample_class", "phenotypeWithinUnit")
  missing_metadata <- setdiff(required_metadata, rownames(column_metadata))
  if (length(missing_metadata)) {
    stop("Animal-level PCA GCT metadata is missing: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  }
  if (!identical(colnames(column_metadata), colnames(expression_matrix))) {
    stop("Animal-level PCA expression and metadata sample columns are not aligned.", call. = FALSE)
  }

  metadata <- data.frame(
    Sample = colnames(expression_matrix),
    AnimalID = trimws(as.character(column_metadata["AnimalID", ])),
    condition_code = trimws(as.character(column_metadata["condition_code", ])),
    condition = normalize_condition(column_metadata["condition", ]),
    sample_class = normalize_sample_class(column_metadata["sample_class", ]),
    phenotypeWithinUnit = trimws(as.character(column_metadata["phenotypeWithinUnit", ])),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (any(is.na(metadata$AnimalID) | !nzchar(metadata$AnimalID))) {
    stop("AnimalID is missing from one or more animal-level PCA samples.", call. = FALSE)
  }
  if (any(is.na(metadata$condition) | !metadata$condition %in% condition_levels) ||
      !setequal(unique(metadata$condition), condition_levels)) {
    stop("Animal-level PCA input does not contain exactly the four canonical conditions.", call. = FALSE)
  }
  if (any(is.na(metadata$sample_class) | !metadata$sample_class %in% sample_classes) ||
      !setequal(unique(metadata$sample_class), sample_classes)) {
    stop("Animal-level PCA input does not contain exactly the four canonical sample classes.", call. = FALSE)
  }
  expected_condition <- unname(condition_code_map[metadata$condition_code])
  if (any(is.na(expected_condition)) || !identical(metadata$condition, expected_condition)) {
    stop("Animal-level PCA condition_code and condition metadata disagree.", call. = FALSE)
  }
  expected_phenotype <- paste(metadata$sample_class, metadata$condition, sep = "_")
  if (!identical(metadata$phenotypeWithinUnit, expected_phenotype)) {
    stop("Animal-level PCA phenotypeWithinUnit metadata is inconsistent.", call. = FALSE)
  }

  unit_key <- paste(metadata$AnimalID, metadata$sample_class, sep = "\r")
  if (anyDuplicated(unit_key)) {
    stop("AnimalID is duplicated within sample_class in the animal-level PCA input.", call. = FALSE)
  }
  animal_condition <- unique(metadata[c("AnimalID", "condition_code", "condition")])
  animal_condition_count <- aggregate(condition ~ AnimalID, animal_condition, function(x) length(unique(x)))
  if (any(animal_condition_count$condition != 1L)) {
    stop("AnimalID maps to multiple conditions in the animal-level PCA input.", call. = FALSE)
  }

  group_counts <- aggregate(
    AnimalID ~ sample_class + condition_code + condition,
    metadata,
    function(x) length(unique(x))
  )
  names(group_counts)[names(group_counts) == "AnimalID"] <- "n_animals"
  expected_groups <- expand.grid(
    sample_class = sample_classes,
    condition = condition_levels,
    stringsAsFactors = FALSE
  )
  expected_groups$condition_code <- names(condition_code_map)[match(expected_groups$condition, condition_code_map)]
  group_counts <- merge(
    expected_groups,
    group_counts,
    by = c("sample_class", "condition_code", "condition"),
    all.x = TRUE,
    sort = FALSE
  )
  group_counts$n_animals[is.na(group_counts$n_animals)] <- 0L
  if (nrow(group_counts) != 16L || any(group_counts$n_animals != expected_n)) {
    stop("Animal-level PCA requires exactly 3 animals per sample_class/condition.", call. = FALSE)
  }

  for (condition_name in condition_levels) {
    condition_metadata <- metadata[metadata$condition == condition_name, , drop = FALSE]
    animal_sets <- lapply(sample_classes, function(class_name) {
      sort(condition_metadata$AnimalID[condition_metadata$sample_class == class_name], method = "radix")
    })
    if (!all(vapply(animal_sets[-1], identical, logical(1), animal_sets[[1]]))) {
      stop(
        "Animal-level PCA has inconsistent AnimalID membership across sample classes for condition ",
        condition_name,
        ".",
        call. = FALSE
      )
    }
  }

  expected_samples <- length(sample_classes) * length(condition_levels) * as.integer(expected_n)
  if (nrow(metadata) != expected_samples || ncol(expression_matrix) != expected_samples) {
    stop("Animal-level PCA input must contain exactly 48 AnimalID x sample_class observations.", call. = FALSE)
  }
  rownames(metadata) <- metadata$Sample

  list(
    expression_matrix = expression_matrix,
    sample_metadata = metadata,
    group_counts = group_counts,
    sampling_unit = "AnimalID_x_sample_class",
    post_aggregation_transformations = "none"
  )
}

prepare_neha_animal_pca <- function(expression_matrix, center = TRUE, scale. = TRUE) {
  expression_matrix <- as.matrix(expression_matrix)
  storage.mode(expression_matrix) <- "numeric"
  if (any(!is.finite(expression_matrix))) {
    stop("PCA input contains nonfinite abundance values; no imputation or missingness filtering is permitted.", call. = FALSE)
  }

  row_variance <- apply(expression_matrix, 1L, stats::var)
  remove <- !is.finite(row_variance) | row_variance <= 0
  removed_rows <- data.frame(
    protein_id = rownames(expression_matrix)[remove],
    reason = rep("zero_variance_in_scaled_pca", sum(remove)),
    stringsAsFactors = FALSE
  )
  pca_matrix <- expression_matrix[!remove, , drop = FALSE]
  if (nrow(pca_matrix) < 2L || ncol(pca_matrix) < 2L) {
    stop("PCA matrix is too small after required zero-variance removal.", call. = FALSE)
  }

  pca <- stats::prcomp(t(pca_matrix), center = center, scale. = scale.)
  if (!identical(rownames(pca$x), colnames(pca_matrix)) || any(!is.finite(pca$x))) {
    stop("PCA scores are nonfinite or misaligned with animal-level samples.", call. = FALSE)
  }

  list(
    matrix = pca_matrix,
    pca = pca,
    removed_rows = removed_rows,
    center = isTRUE(center),
    scale = isTRUE(scale.),
    row_removal_policy = "zero_variance_only_required_for_scaled_pca",
    abundance_normalization = "none",
    abundance_imputation = "none"
  )
}

make_neha_pca_audit <- function(
    validated_input = NULL,
    prepared_pca = NULL,
    source_path,
    source_sha256,
    output_paths,
    execution_status,
    error_message = NA_character_) {
  counts <- if (is.null(validated_input)) {
    data.frame(
      sample_class = NA_character_, condition_code = NA_character_, condition = NA_character_,
      n_animals = NA_integer_, stringsAsFactors = FALSE
    )
  } else {
    validated_input$group_counts
  }
  removed <- if (is.null(prepared_pca) || !nrow(prepared_pca$removed_rows)) {
    data.frame(protein_id = character(), reason = character(), stringsAsFactors = FALSE)
  } else {
    prepared_pca$removed_rows
  }
  variance <- if (is.null(prepared_pca)) c(NA_real_, NA_real_) else {
    values <- prepared_pca$pca$sdev^2 / sum(prepared_pca$pca$sdev^2)
    values[seq_len(min(2L, length(values)))]
  }
  if (length(variance) < 2L) variance <- c(variance, rep(NA_real_, 2L - length(variance)))

  counts$source_gct_path <- neha_pca_normalize_path(source_path, must_work = FALSE)
  counts$source_gct_sha256 <- source_sha256
  counts$n_protein_rows_loaded <- if (is.null(validated_input)) NA_integer_ else nrow(validated_input$expression_matrix)
  counts$n_samples_loaded <- if (is.null(validated_input)) NA_integer_ else ncol(validated_input$expression_matrix)
  counts$n_samples_entering_pca <- if (is.null(prepared_pca)) NA_integer_ else ncol(prepared_pca$matrix)
  counts$sample_classes <- paste(sample_classes, collapse = ";")
  counts$conditions <- paste(condition_levels, collapse = ";")
  counts$pca_center <- if (is.null(prepared_pca)) NA else prepared_pca$center
  counts$pca_scale <- if (is.null(prepared_pca)) NA else prepared_pca$scale
  counts$n_rows_removed_for_pca <- nrow(removed)
  counts$rows_removed_for_pca <- if (nrow(removed)) paste(removed$protein_id, collapse = ";") else "none"
  counts$row_removal_reason <- if (nrow(removed)) paste(unique(removed$reason), collapse = ";") else "none"
  counts$n_proteins_entering_pca <- if (is.null(prepared_pca)) NA_integer_ else nrow(prepared_pca$matrix)
  counts$abundance_normalization <- if (is.null(prepared_pca)) NA_character_ else prepared_pca$abundance_normalization
  counts$abundance_imputation <- if (is.null(prepared_pca)) NA_character_ else prepared_pca$abundance_imputation
  counts$pc1_variance_explained <- variance[[1]]
  counts$pc2_variance_explained <- variance[[2]]
  counts$output_paths <- paste(
    paste(names(output_paths), neha_pca_normalize_path(output_paths, must_work = FALSE), sep = "="),
    collapse = ";"
  )
  counts$execution_status <- execution_status
  counts$error_message <- error_message
  counts
}

prepare_neha_pca_variance_treemap_data <- function(variance_explained, max_pcs = 20L) {
  variance_explained <- as.numeric(variance_explained)
  max_pcs <- as.integer(max_pcs)
  if (length(max_pcs) != 1L || is.na(max_pcs) || max_pcs < 1L) {
    stop("max_pcs must be one positive integer.", call. = FALSE)
  }

  n_input <- min(max_pcs, length(variance_explained))
  if (!n_input) {
    result <- data.frame(
      PC = character(), Category = character(), Variance = numeric(),
      stringsAsFactors = FALSE
    )
    attr(result, "n_input_rows") <- 0L
    return(result)
  }

  pc_index <- seq_len(n_input)
  result <- data.frame(
    PC = paste0("PC", pc_index),
    Category = ifelse(
      pc_index <= 5L,
      "Top 5",
      ifelse(pc_index <= 10L, "PC 6-10", "PC 11-20")
    ),
    Variance = variance_explained[pc_index],
    stringsAsFactors = FALSE
  )
  result <- result[is.finite(result$Variance) & result$Variance > 0, , drop = FALSE]
  rownames(result) <- NULL
  attr(result, "n_input_rows") <- n_input
  result
}

prepare_neha_pca_group_centroids <- function(pca_scores, metadata, group_key, n_pcs = 5L) {
  pca_scores <- as.matrix(pca_scores)
  if (is.null(rownames(pca_scores)) || anyDuplicated(rownames(pca_scores))) {
    stop("PCA scores require unique sample row names for optional group plots.", call. = FALSE)
  }
  if (!group_key %in% names(metadata)) {
    stop("Optional group key is absent from PCA metadata: ", group_key, call. = FALSE)
  }
  if (is.null(rownames(metadata)) || !all(rownames(pca_scores) %in% rownames(metadata))) {
    stop("Optional group metadata is not keyed to the PCA score rows.", call. = FALSE)
  }

  n_pcs <- min(as.integer(n_pcs), ncol(pca_scores))
  if (is.na(n_pcs) || n_pcs < 1L) {
    stop("n_pcs must select at least one PCA score column.", call. = FALSE)
  }
  sample_ids <- rownames(pca_scores)
  groups <- trimws(as.character(metadata[sample_ids, group_key]))
  usable <- !is.na(groups) & nzchar(groups) & apply(pca_scores[, seq_len(n_pcs), drop = FALSE], 1L, function(x) all(is.finite(x)))
  if (!any(usable)) {
    result <- data.frame(group = character(), stringsAsFactors = FALSE)
    attr(result, "n_input_rows") <- length(sample_ids)
    return(result)
  }

  centroids <- stats::aggregate(
    pca_scores[usable, seq_len(n_pcs), drop = FALSE],
    by = list(group = groups[usable]),
    FUN = mean
  )
  attr(centroids, "n_input_rows") <- length(sample_ids)
  centroids
}

run_neha_pca_optional_plot <- function(
    plot_name,
    usable_data,
    output_path,
    render_fun,
    n_input_rows = nrow(usable_data)) {
  if (!is.data.frame(usable_data)) {
    stop("Optional PCA plot data must be a data.frame.", call. = FALSE)
  }
  if (!is.function(render_fun)) {
    stop("Optional PCA plot renderer must be a function.", call. = FALSE)
  }

  make_result <- function(status, reason) {
    data.frame(
      visualization = as.character(plot_name),
      execution_status = status,
      reason = reason,
      n_input_rows = as.integer(n_input_rows),
      n_usable_rows = nrow(usable_data),
      output_path = neha_pca_normalize_path(output_path, must_work = FALSE),
      stringsAsFactors = FALSE
    )
  }

  if (!nrow(usable_data)) {
    reason <- "no usable rows after the visualization's existing finite/positive/group filters"
    message("Skipping optional PCA visualization '", plot_name, "': ", reason, ".")
    return(make_result("skipped_no_usable_rows", reason))
  }

  error <- tryCatch({
    render_fun(usable_data, output_path)
    NULL
  }, error = identity)
  if (inherits(error, "error")) {
    reason <- paste0("optional renderer failed: ", conditionMessage(error))
    message("Skipping optional PCA visualization '", plot_name, "': ", reason, ".")
    return(make_result("skipped_render_error", reason))
  }

  make_result("success", NA_character_)
}
