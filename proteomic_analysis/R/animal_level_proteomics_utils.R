canonicalize_hemisphere <- function(x, field = "ReplicateGroup") {
  raw <- trimws(as.character(x))
  key <- tolower(raw)
  out <- rep(NA_character_, length(raw))
  out[key %in% c("left", "l")] <- "Left"
  out[key %in% c("right", "r")] <- "Right"

  invalid <- is.na(out)
  if (any(invalid)) {
    bad <- unique(raw[invalid])
    bad[is.na(bad) | !nzchar(bad)] <- "<missing>"
    stop(
      field, " must contain only Left/L or Right/R; found: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }
  out
}

parse_sample_id_hemisphere <- function(sample_id) {
  sample_id <- as.character(sample_id)
  has_left <- !is.na(sample_id) & grepl("_L_", sample_id, fixed = TRUE)
  has_right <- !is.na(sample_id) & grepl("_R_", sample_id, fixed = TRUE)
  ambiguous <- has_left & has_right
  missing <- !has_left & !has_right

  if (any(ambiguous | missing)) {
    bad <- sample_id[ambiguous | missing]
    reason <- ifelse(ambiguous[ambiguous | missing], "ambiguous", "missing")
    stop(
      "sample_id hemisphere token is missing or ambiguous: ",
      paste(paste0(bad, " [", reason, "]"), collapse = "; "),
      call. = FALSE
    )
  }
  ifelse(has_left, "Left", "Right")
}

normalize_exclude_flag <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  key <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(key))
  out[is.na(key) | key %in% c("", "false", "f", "0", "no", "n")] <- FALSE
  out[key %in% c("true", "t", "1", "yes", "y")] <- TRUE
  if (anyNA(out)) {
    stop(
      "exclude contains unrecognized values: ",
      paste(unique(as.character(x[is.na(out)])), collapse = ", "),
      call. = FALSE
    )
  }
  out
}

sanitize_output_token <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "missing")
}

make_output_column_name <- function(metadata, extra_unit_cols = character()) {
  parts <- c("AnimalID", "sample_class", extra_unit_cols)
  apply(metadata[, parts, drop = FALSE], 1, function(x) {
    paste(sanitize_output_token(x), collapse = "_")
  })
}

validate_animal_level_metadata <- function(
    metadata,
    extra_unit_cols = character(),
    hemisphere_crosscheck = c("sample_id", "explicit_evidence"),
    hemisphere_evidence_col = NULL) {
  hemisphere_crosscheck <- match.arg(hemisphere_crosscheck)
  required <- c(
    "sample_id", "AnimalID", "ReplicateGroup", "sample_class",
    "condition_code", "condition", extra_unit_cols
  )
  missing_cols <- setdiff(required, names(metadata))
  if (length(missing_cols) > 0) {
    stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
  for (field in c("sample_id", "AnimalID", "sample_class", "condition_code", "condition", extra_unit_cols)) {
    value <- trimws(as.character(metadata[[field]]))
    if (any(is.na(value) | !nzchar(value))) {
      stop(field, " is missing for one or more metadata rows.", call. = FALSE)
    }
    metadata[[field]] <- value
  }
  if (anyDuplicated(metadata$sample_id)) {
    stop(
      "source sample_id values must be unique: ",
      paste(unique(metadata$sample_id[duplicated(metadata$sample_id)]), collapse = ", "),
      call. = FALSE
    )
  }

  metadata$hemisphere <- canonicalize_hemisphere(metadata$ReplicateGroup)
  if (hemisphere_crosscheck == "sample_id") {
    evidence <- parse_sample_id_hemisphere(metadata$sample_id)
  } else {
    if (is.null(hemisphere_evidence_col) || !hemisphere_evidence_col %in% names(metadata)) {
      stop("explicit_evidence cross-check requires hemisphere_evidence_col.", call. = FALSE)
    }
    evidence <- canonicalize_hemisphere(
      metadata[[hemisphere_evidence_col]],
      field = hemisphere_evidence_col
    )
  }
  disagreement <- metadata$hemisphere != evidence
  if (any(disagreement)) {
    stop(
      "ReplicateGroup disagrees with hemisphere evidence for sample_id: ",
      paste(metadata$sample_id[disagreement], collapse = ", "),
      call. = FALSE
    )
  }
  metadata$hemisphere_crosscheck <- evidence
  metadata$hemisphere_crosscheck_mode <- hemisphere_crosscheck

  if (!"exclude" %in% names(metadata)) metadata$exclude <- FALSE
  metadata$exclude <- normalize_exclude_flag(metadata$exclude)

  animal_condition <- unique(metadata[c("AnimalID", "condition_code", "condition")])
  condition_counts <- aggregate(
    paste(condition_code, condition, sep = "\r") ~ AnimalID,
    animal_condition,
    function(x) length(unique(x))
  )
  if (any(condition_counts[[2]] != 1L)) {
    stop(
      "AnimalID maps to multiple conditions: ",
      paste(condition_counts$AnimalID[condition_counts[[2]] != 1L], collapse = ", "),
      call. = FALSE
    )
  }
  metadata
}

make_expected_animal_units <- function(metadata, sample_class_levels = NULL, extra_unit_cols = character()) {
  if (length(extra_unit_cols) > 0) {
    stop(
      "Expected units with extra sampling variables must be supplied explicitly; automatic cross-products are unsafe.",
      call. = FALSE
    )
  }
  animals <- unique(metadata[c("AnimalID", "condition_code", "condition")])
  classes <- if (is.null(sample_class_levels)) sort(unique(metadata$sample_class)) else sample_class_levels
  out <- merge(animals, data.frame(sample_class = classes, stringsAsFactors = FALSE), by = NULL)
  out[c("AnimalID", "sample_class", "condition_code", "condition")]
}

prepare_animal_level_aggregation <- function(
    metadata,
    dataset_project,
    extra_unit_cols = character(),
    expected_units = NULL,
    sample_class_levels = NULL,
    hemisphere_crosscheck = c("sample_id", "explicit_evidence"),
    hemisphere_evidence_col = NULL) {
  hemisphere_crosscheck <- match.arg(hemisphere_crosscheck)
  metadata <- validate_animal_level_metadata(
    metadata,
    extra_unit_cols = extra_unit_cols,
    hemisphere_crosscheck = hemisphere_crosscheck,
    hemisphere_evidence_col = hemisphere_evidence_col
  )
  key_cols <- c("AnimalID", "sample_class", extra_unit_cols)

  if (is.null(expected_units)) {
    expected_units <- unique(metadata[c(key_cols, "condition_code", "condition")])
  } else {
    expected_units <- as.data.frame(expected_units, stringsAsFactors = FALSE, check.names = FALSE)
    required_expected <- c(key_cols, "condition_code", "condition")
    missing_expected <- setdiff(required_expected, names(expected_units))
    if (length(missing_expected) > 0) {
      stop("expected_units is missing: ", paste(missing_expected, collapse = ", "), call. = FALSE)
    }
    expected_units <- unique(expected_units[required_expected])
  }

  observed_units <- unique(metadata[c(key_cols, "condition_code", "condition")])
  observed_key <- do.call(paste, c(observed_units[key_cols], sep = "\r"))
  expected_key <- do.call(paste, c(expected_units[key_cols], sep = "\r"))
  unexpected <- setdiff(observed_key, expected_key)
  if (length(unexpected) > 0) {
    stop("One or more source samples do not map to an expected output unit.", call. = FALSE)
  }

  condition_lookup <- unique(metadata[c(key_cols, "condition_code", "condition")])
  condition_key <- do.call(paste, c(condition_lookup[key_cols], sep = "\r"))
  if (anyDuplicated(condition_key)) {
    stop("An analysis unit maps to multiple conditions.", call. = FALSE)
  }
  idx <- match(expected_key, condition_key)
  has_observed <- !is.na(idx)
  mismatch <- has_observed & (
    expected_units$condition_code != condition_lookup$condition_code[idx] |
      expected_units$condition != condition_lookup$condition[idx]
  )
  if (any(mismatch)) stop("expected_units condition metadata disagrees with source metadata.", call. = FALSE)

  class_rank <- if (is.null(sample_class_levels)) {
    match(expected_units$sample_class, sort(unique(expected_units$sample_class)))
  } else {
    match(expected_units$sample_class, sample_class_levels)
  }
  condition_rank <- suppressWarnings(as.numeric(expected_units$condition_code))
  condition_rank[is.na(condition_rank)] <- match(
    expected_units$condition_code[is.na(condition_rank)],
    sort(unique(expected_units$condition_code))
  )
  order_args <- c(list(class_rank, condition_rank, expected_units$AnimalID), expected_units[extra_unit_cols])
  expected_units <- expected_units[do.call(order, order_args), , drop = FALSE]
  rownames(expected_units) <- NULL

  output_names <- make_output_column_name(expected_units, extra_unit_cols)
  if (anyDuplicated(output_names)) {
    stop("Animal-level output column IDs are not unique.", call. = FALSE)
  }

  audit_rows <- lapply(seq_len(nrow(expected_units)), function(i) {
    unit <- expected_units[i, , drop = FALSE]
    match_unit <- rep(TRUE, nrow(metadata))
    for (field in key_cols) match_unit <- match_unit & metadata[[field]] == unit[[field]]
    source <- metadata[match_unit, , drop = FALSE]
    included <- source[!source$exclude, , drop = FALSE]
    left <- included$sample_id[included$hemisphere == "Left"]
    right <- included$sample_id[included$hemisphere == "Right"]

    if (length(left) > 1L) {
      stop("duplicate Left samples remain for ", output_names[i], ": ", paste(left, collapse = ", "), call. = FALSE)
    }
    if (length(right) > 1L) {
      stop("duplicate Right samples remain for ", output_names[i], ": ", paste(right, collapse = ", "), call. = FALSE)
    }

    status <- if (length(left) == 1L && length(right) == 1L) {
      "bilateral_complete"
    } else if (length(left) == 1L) {
      "left_only_observed"
    } else if (length(right) == 1L) {
      "right_only_observed"
    } else {
      "missing_both"
    }
    method <- switch(
      status,
      bilateral_complete = "equal_weight_mean_LR_on_existing_imputed_log2_values",
      left_only_observed = "single_observed_hemisphere_no_imputation",
      right_only_observed = "single_observed_hemisphere_no_imputation",
      missing_both = NA_character_
    )
    inclusion <- switch(
      status,
      bilateral_complete = "included_primary_and_sensitivity",
      left_only_observed = "included_primary_only_incomplete_pair",
      right_only_observed = "included_primary_only_incomplete_pair",
      missing_both = "excluded_missing_both"
    )
    canonical_parts <- c(dataset_project, unlist(unit[key_cols], use.names = FALSE))
    row <- data.frame(
      dataset_project = dataset_project,
      AnimalID = unit$AnimalID,
      condition_code = unit$condition_code,
      condition = unit$condition,
      sample_class = unit$sample_class,
      canonical_analysis_unit = paste(sanitize_output_token(canonical_parts), collapse = "__"),
      left_sample = if (length(left)) left else NA_character_,
      right_sample = if (length(right)) right else NA_character_,
      n_left_source_samples = length(left),
      n_right_source_samples = length(right),
      hemisphere_status = status,
      aggregation_method = method,
      output_column_name = output_names[i],
      inclusion_status = inclusion,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    for (field in extra_unit_cols) row[[field]] <- unit[[field]]
    row
  })
  audit <- do.call(rbind, audit_rows)
  names(audit)[names(audit) == "dataset_project"] <- "dataset/project"

  unit_key <- do.call(paste, c(expected_units[key_cols], sep = "\r"))
  unit_output <- setNames(audit$output_column_name, unit_key)
  metadata_key <- do.call(paste, c(metadata[key_cols], sep = "\r"))
  assignment <- metadata[c(
    "sample_id", "AnimalID", "condition_code", "condition", "sample_class",
    extra_unit_cols, "hemisphere", "exclude", "hemisphere_crosscheck_mode"
  )]
  assignment$output_column_name <- unname(unit_output[metadata_key])
  assignment$source_assignment_status <- ifelse(
    assignment$exclude,
    "excluded_by_metadata",
    "included_in_primary"
  )
  if (any(is.na(assignment$output_column_name))) {
    stop("One or more source samples map to no output unit.", call. = FALSE)
  }
  included_assignment <- assignment[assignment$source_assignment_status == "included_in_primary", , drop = FALSE]
  if (anyDuplicated(included_assignment$sample_id)) {
    stop("A source sample maps to more than one output unit.", call. = FALSE)
  }

  list(metadata = metadata, audit = audit, source_assignment = assignment)
}

aggregate_animal_level_matrix <- function(expression_matrix, protein_ids, aggregation_plan) {
  expression_matrix <- as.matrix(expression_matrix)
  storage.mode(expression_matrix) <- "numeric"
  protein_ids <- as.character(protein_ids)
  if (nrow(expression_matrix) != length(protein_ids)) {
    stop("protein_ids length does not match expression matrix rows.", call. = FALSE)
  }
  if (anyDuplicated(protein_ids) || any(is.na(protein_ids) | !nzchar(protein_ids))) {
    stop("Protein IDs must be complete and unique.", call. = FALSE)
  }
  if (is.null(colnames(expression_matrix)) || anyDuplicated(colnames(expression_matrix))) {
    stop("Expression matrix source sample columns must be named and unique.", call. = FALSE)
  }

  assignment <- aggregation_plan$source_assignment
  included_samples <- assignment$sample_id[assignment$source_assignment_status == "included_in_primary"]
  missing_samples <- setdiff(included_samples, colnames(expression_matrix))
  if (length(missing_samples) > 0) {
    stop("Metadata samples missing from expression matrix: ", paste(missing_samples, collapse = ", "), call. = FALSE)
  }
  extra_samples <- setdiff(colnames(expression_matrix), aggregation_plan$metadata$sample_id)
  if (length(extra_samples) > 0) {
    stop("Expression samples missing from metadata: ", paste(extra_samples, collapse = ", "), call. = FALSE)
  }

  audit <- aggregation_plan$audit
  primary_rows <- audit$hemisphere_status != "missing_both"
  primary_audit <- audit[primary_rows, , drop = FALSE]
  primary <- vapply(seq_len(nrow(primary_audit)), function(i) {
    samples <- stats::na.omit(c(primary_audit$left_sample[i], primary_audit$right_sample[i]))
    if (length(samples) == 2L) {
      rowMeans(expression_matrix[, samples, drop = FALSE])
    } else {
      expression_matrix[, samples, drop = TRUE]
    }
  }, numeric(nrow(expression_matrix)))
  if (is.null(dim(primary))) primary <- matrix(primary, ncol = 1L)
  colnames(primary) <- primary_audit$output_column_name
  rownames(primary) <- protein_ids

  sensitivity_audit <- audit[audit$hemisphere_status == "bilateral_complete", , drop = FALSE]
  sensitivity <- vapply(seq_len(nrow(sensitivity_audit)), function(i) {
    rowMeans(expression_matrix[, c(sensitivity_audit$left_sample[i], sensitivity_audit$right_sample[i]), drop = FALSE])
  }, numeric(nrow(expression_matrix)))
  if (is.null(dim(sensitivity))) sensitivity <- matrix(sensitivity, ncol = nrow(sensitivity_audit))
  colnames(sensitivity) <- sensitivity_audit$output_column_name
  rownames(sensitivity) <- protein_ids

  list(
    primary = primary,
    sensitivity = sensitivity,
    primary_audit = primary_audit,
    sensitivity_audit = sensitivity_audit,
    protein_ids = protein_ids
  )
}

summarize_animal_level_design <- function(aggregation_audit) {
  split_rows <- split(
    aggregation_audit,
    interaction(
      aggregation_audit$sample_class,
      aggregation_audit$condition_code,
      drop = TRUE,
      lex.order = TRUE
    )
  )
  out <- lapply(split_rows, function(x) {
    data.frame(
      sample_class = x$sample_class[1],
      condition_code = x$condition_code[1],
      condition = x$condition[1],
      n_unique_animals = length(unique(x$AnimalID[x$hemisphere_status != "missing_both"])),
      n_left_samples = sum(x$n_left_source_samples),
      n_right_samples = sum(x$n_right_source_samples),
      complete_bilateral_pairs = sum(x$hemisphere_status == "bilateral_complete"),
      one_sided_pairs = sum(x$hemisphere_status %in% c("left_only_observed", "right_only_observed")),
      missing_pairs = sum(x$hemisphere_status == "missing_both"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$sample_class, suppressWarnings(as.numeric(out$condition_code))), , drop = FALSE]
}
