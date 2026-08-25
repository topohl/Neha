# Focused sampling-unit contracts for the Neha animal-level EWCE branch.

neha_ewce_normalize_path <- function(path, must_work = FALSE) {
  normalizePath(path, winslash = "/", mustWork = must_work)
}

neha_ewce_path_is_within <- function(path, root) {
  path <- tolower(neha_ewce_normalize_path(path, must_work = FALSE))
  root <- sub("/+$", "", tolower(neha_ewce_normalize_path(root, must_work = FALSE)))
  identical(path, root) || startsWith(path, paste0(root, "/"))
}

validate_neha_ewce_output_root <- function(output_root, historical_root) {
  output_root <- neha_ewce_normalize_path(output_root, must_work = FALSE)
  historical_root <- neha_ewce_normalize_path(historical_root, must_work = FALSE)
  if (neha_ewce_path_is_within(output_root, historical_root)) {
    stop(
      "Animal-level EWCE output cannot overwrite the historical EWCE root: ",
      historical_root,
      call. = FALSE
    )
  }
  output_root
}

validate_neha_ewce_animal_input <- function(
    parsed_gct,
    source_path,
    expected_n = 3L,
    expected_sample_classes = sample_classes,
    expected_conditions = condition_levels) {
  if (!is.list(parsed_gct) || is.null(parsed_gct$matrix) || is.null(parsed_gct$column_metadata)) {
    stop("Animal-level EWCE input must be a parsed ProTigy GCT.", call. = FALSE)
  }
  expression_matrix <- as.matrix(parsed_gct$matrix)
  storage.mode(expression_matrix) <- "numeric"
  if (!nrow(expression_matrix) || !ncol(expression_matrix)) {
    stop("Animal-level EWCE expression matrix is empty.", call. = FALSE)
  }
  if (is.null(rownames(expression_matrix)) || anyDuplicated(rownames(expression_matrix)) ||
      any(is.na(rownames(expression_matrix)) | !nzchar(rownames(expression_matrix)))) {
    stop("Animal-level EWCE protein identifiers must be complete and unique.", call. = FALSE)
  }
  if (is.null(colnames(expression_matrix)) || anyDuplicated(colnames(expression_matrix)) ||
      any(is.na(colnames(expression_matrix)) | !nzchar(colnames(expression_matrix)))) {
    stop("Animal-level EWCE sample columns must be complete and unique.", call. = FALSE)
  }
  if (any(!is.finite(expression_matrix))) {
    stop("Animal-level EWCE abundance matrix contains nonfinite values; no filtering or imputation is permitted.", call. = FALSE)
  }

  column_metadata <- as.matrix(parsed_gct$column_metadata)
  required_metadata <- c("AnimalID", "condition_code", "condition", "sample_class", "phenotypeWithinUnit")
  missing_metadata <- setdiff(required_metadata, rownames(column_metadata))
  if (length(missing_metadata)) {
    stop("Animal-level EWCE GCT metadata is missing: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  }
  if (!identical(colnames(column_metadata), colnames(expression_matrix))) {
    stop("Animal-level EWCE expression and metadata sample columns are not aligned.", call. = FALSE)
  }

  animal_id <- trimws(as.character(column_metadata["AnimalID", ]))
  condition_code <- trimws(as.character(column_metadata["condition_code", ]))
  condition <- normalize_condition(column_metadata["condition", ])
  sample_class <- normalize_sample_class(column_metadata["sample_class", ])
  phenotype <- trimws(as.character(column_metadata["phenotypeWithinUnit", ]))
  if (any(is.na(animal_id) | !nzchar(animal_id))) {
    stop("AnimalID is missing from one or more animal-level EWCE samples.", call. = FALSE)
  }
  if (any(is.na(condition) | !condition %in% expected_conditions)) {
    stop("Animal-level EWCE input contains an unsupported or missing condition.", call. = FALSE)
  }
  if (any(is.na(sample_class) | !sample_class %in% expected_sample_classes)) {
    stop("Animal-level EWCE input contains an unsupported or missing sample_class.", call. = FALSE)
  }
  expected_condition <- unname(condition_code_map[condition_code])
  if (any(is.na(expected_condition)) || !identical(condition, expected_condition)) {
    stop("Animal-level EWCE condition_code and condition metadata disagree.", call. = FALSE)
  }
  if (!identical(phenotype, paste(sample_class, condition, sep = "_"))) {
    stop("Animal-level EWCE phenotypeWithinUnit metadata is inconsistent.", call. = FALSE)
  }

  sample_metadata <- data.frame(
    Sample = colnames(expression_matrix),
    AnimalID = animal_id,
    Stratum = sample_class,
    Cond = condition,
    condition_code = condition_code,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  unit_key <- paste(sample_metadata$AnimalID, sample_metadata$Stratum, sample_metadata$Cond, sep = "\r")
  if (anyDuplicated(unit_key)) {
    duplicate_rows <- duplicated(unit_key) | duplicated(unit_key, fromLast = TRUE)
    stop(
      "AnimalID is duplicated within sample_class/condition: ",
      paste(unique(unit_key[duplicate_rows]), collapse = "; "),
      call. = FALSE
    )
  }

  animal_condition <- unique(sample_metadata[c("AnimalID", "Cond")])
  condition_count <- aggregate(Cond ~ AnimalID, animal_condition, function(x) length(unique(x)))
  if (any(condition_count$Cond != 1L)) {
    stop("AnimalID maps to multiple conditions in the animal-level EWCE input.", call. = FALSE)
  }

  observed_counts <- aggregate(
    AnimalID ~ Stratum + Cond,
    sample_metadata,
    function(x) length(unique(x))
  )
  names(observed_counts)[names(observed_counts) == "AnimalID"] <- "n_animals"
  expected_grid <- expand.grid(
    Stratum = expected_sample_classes,
    Cond = expected_conditions,
    stringsAsFactors = FALSE
  )
  count_audit <- merge(expected_grid, observed_counts, by = c("Stratum", "Cond"), all.x = TRUE, sort = FALSE)
  count_audit$n_animals[is.na(count_audit$n_animals)] <- 0L
  extra_group <- merge(observed_counts[c("Stratum", "Cond")], expected_grid, by = c("Stratum", "Cond"), all.x = TRUE)
  if (nrow(observed_counts) != nrow(expected_grid) || any(count_audit$n_animals != expected_n) ||
      nrow(extra_group) != nrow(expected_grid)) {
    stop(
      "Animal-level EWCE requires exactly ", expected_n,
      " animals per sample_class/condition.",
      call. = FALSE
    )
  }

  for (condition_name in expected_conditions) {
    condition_meta <- sample_metadata[sample_metadata$Cond == condition_name, , drop = FALSE]
    animal_sets <- lapply(expected_sample_classes, function(class_name) {
      sort(unique(condition_meta$AnimalID[condition_meta$Stratum == class_name]), method = "radix")
    })
    if (!all(vapply(animal_sets[-1], identical, logical(1), animal_sets[[1]]))) {
      stop(
        "Animal-level EWCE has inconsistent AnimalID membership across sample classes for condition ",
        condition_name, ".",
        call. = FALSE
      )
    }
  }

  expected_columns <- length(expected_sample_classes) * length(expected_conditions) * as.integer(expected_n)
  if (ncol(expression_matrix) != expected_columns) {
    stop("Animal-level EWCE sample count is inconsistent with the validated design.", call. = FALSE)
  }
  sample_metadata$Cond <- factor(sample_metadata$Cond, levels = expected_conditions)
  sample_metadata$Stratum <- factor(sample_metadata$Stratum, levels = expected_sample_classes)

  list(
    expression_matrix = expression_matrix,
    sample_metadata = sample_metadata,
    count_audit = count_audit,
    source_path = neha_ewce_normalize_path(source_path, must_work = FALSE),
    sampling_unit = "AnimalID_x_sample_class",
    aggregation_policy = "equal_weight_mean_LR_on_existing_imputed_log2_values",
    transformations_after_aggregation = "none"
  )
}

neha_ewce_contrast_metadata <- function(sample_class, sample_metadata = NULL) {
  manifest <- neha_primary_contrast_manifest()
  manifest <- manifest[manifest$sample_class == sample_class, , drop = FALSE]
  if (nrow(manifest) != 3L) {
    stop("EWCE requires exactly three shared Neha contrasts per sample class.", call. = FALSE)
  }
  out <- data.frame(
    canonical_comparison = manifest$canonical_comparison,
    canonical_contrast = manifest$canonical_contrast,
    sample_class = manifest$sample_class,
    numerator_condition = manifest$case_condition,
    denominator_condition = manifest$reference_condition,
    numerator_animal_n = NA_integer_,
    denominator_animal_n = NA_integer_,
    numerator_animal_ids = NA_character_,
    denominator_animal_ids = NA_character_,
    animal_ids_used = NA_character_,
    stringsAsFactors = FALSE
  )
  if (!is.null(sample_metadata)) {
    for (i in seq_len(nrow(out))) {
      numerator_ids <- sort(unique(as.character(sample_metadata$AnimalID[
        as.character(sample_metadata$Cond) == out$numerator_condition[[i]]
      ])), method = "radix")
      denominator_ids <- sort(unique(as.character(sample_metadata$AnimalID[
        as.character(sample_metadata$Cond) == out$denominator_condition[[i]]
      ])), method = "radix")
      out$numerator_animal_n[[i]] <- length(numerator_ids)
      out$denominator_animal_n[[i]] <- length(denominator_ids)
      out$numerator_animal_ids[[i]] <- paste(numerator_ids, collapse = ";")
      out$denominator_animal_ids[[i]] <- paste(denominator_ids, collapse = ";")
      out$animal_ids_used[[i]] <- paste(sort(unique(c(numerator_ids, denominator_ids)), method = "radix"), collapse = ";")
    }
  }
  out
}

prepare_neha_ewce_limma_stratum <- function(
    expression_matrix,
    sample_metadata,
    sample_class,
    expected_n = 3L) {
  meta <- as.data.frame(sample_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Sample", "AnimalID", "Stratum", "Cond")
  missing <- setdiff(required, names(meta))
  if (length(missing)) stop("EWCE sample metadata is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  meta <- meta[as.character(meta$Stratum) == sample_class, , drop = FALSE]
  if (nrow(meta) != length(condition_levels) * expected_n) {
    stop("Limma requires exactly one animal-level sample for each expected AnimalID x sample_class unit.", call. = FALSE)
  }
  unit_key <- paste(meta$AnimalID, as.character(meta$Stratum), as.character(meta$Cond), sep = "\r")
  if (anyDuplicated(unit_key)) {
    stop("Duplicated AnimalID within sample_class/condition cannot enter limma.", call. = FALSE)
  }
  counts <- table(factor(as.character(meta$Cond), levels = condition_levels))
  if (any(counts != expected_n)) {
    stop("Limma requires exactly three animals per condition and sample_class.", call. = FALSE)
  }
  if (anyDuplicated(meta$Sample) || any(is.na(meta$Sample) | !nzchar(meta$Sample))) {
    stop("Animal-level limma sample identifiers must be complete and unique.", call. = FALSE)
  }

  expression_matrix <- as.matrix(expression_matrix)
  storage.mode(expression_matrix) <- "numeric"
  sample_index <- match(meta$Sample, colnames(expression_matrix))
  if (anyNA(sample_index)) stop("Animal-level limma metadata does not match the expression columns.", call. = FALSE)
  meta <- meta[order(match(as.character(meta$Cond), condition_levels), meta$AnimalID, meta$Sample, method = "radix"), , drop = FALSE]
  x <- expression_matrix[, meta$Sample, drop = FALSE]
  if (any(!is.finite(x))) {
    stop("Animal-level limma input contains nonfinite abundances; no filtering or imputation is permitted.", call. = FALSE)
  }

  meta$Cond <- factor(as.character(meta$Cond), levels = condition_levels)
  rownames(meta) <- meta$Sample
  design <- stats::model.matrix(~ 0 + Cond, data = meta)
  colnames(design) <- sub("^Cond", "", colnames(design))
  if (!identical(colnames(design), condition_levels)) {
    stop("Animal-level limma design does not contain the four canonical conditions in order.", call. = FALSE)
  }
  contrast_metadata <- neha_ewce_contrast_metadata(sample_class, meta)
  contrast_matrix <- matrix(
    0,
    nrow = ncol(design),
    ncol = nrow(contrast_metadata),
    dimnames = list(colnames(design), contrast_metadata$canonical_comparison)
  )
  for (i in seq_len(nrow(contrast_metadata))) {
    contrast_matrix[contrast_metadata$numerator_condition[[i]], i] <- 1
    contrast_matrix[contrast_metadata$denominator_condition[[i]], i] <- -1
  }
  if (any(contrast_metadata$numerator_animal_n != expected_n) ||
      any(contrast_metadata$denominator_animal_n != expected_n)) {
    stop("A canonical EWCE contrast does not have the expected three animals per group.", call. = FALSE)
  }

  list(
    expression_matrix = x,
    sample_metadata = meta,
    design = design,
    contrast_matrix = contrast_matrix,
    contrast_metadata = contrast_metadata,
    n_proteins = nrow(x),
    n_animal_observations = ncol(x),
    sampling_unit = "animal"
  )
}

make_neha_ewce_differential_audit <- function(
    contrast_metadata,
    source_path,
    output_path,
    n_proteins,
    analysis_params,
    execution_status,
    error_message = NA_character_) {
  out <- as.data.frame(contrast_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  out$source_animal_level_input <- neha_ewce_normalize_path(source_path, must_work = FALSE)
  out$n_proteins_entering_differential_analysis <- as.integer(n_proteins)
  out$differential_output_path <- neha_ewce_normalize_path(output_path, must_work = FALSE)
  out$execution_status <- execution_status
  out$error_message <- error_message
  out$sampling_unit <- "animal"
  out$aggregation_policy <- "equal_weight_mean_LR_on_existing_imputed_log2_values"
  out$post_aggregation_normalization <- "none"
  out$post_aggregation_filtering <- "none"
  out$post_aggregation_imputation <- "none"
  out$limma_design <- "~0+condition;one_column_per_AnimalID_x_sample_class"
  out$limma_ebayes <- "trend=TRUE;robust=TRUE"
  out$rank_direction <- "positive_logFC_and_t_are_higher_in_canonical_numerator"
  out$ewce_seed <- analysis_params$seed
  out$ewce_reps <- analysis_params$reps
  out$ewce_annotation_levels <- paste(analysis_params$annot_levels, collapse = ";")
  out$ewce_primary_annotation_level <- analysis_params$primary_annot_level
  out$ewce_top_n_values <- paste(analysis_params$top_n_values, collapse = ";")
  out$ewce_primary_top_n <- analysis_params$primary_top_n
  out$ewce_fdr_alpha <- analysis_params$fdr_alpha
  out$ewce_marker_top_n <- analysis_params$marker_top_n
  out$canonical_conditions <- paste(analysis_params$conditions, collapse = ";")
  out$reference_condition_metadata <- analysis_params$reference_condition
  out
}
