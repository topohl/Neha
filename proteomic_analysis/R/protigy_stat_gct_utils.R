# Generic parsing and validation helpers for Neha ProTigy statistical-result GCTs.
# Dataset-specific class and condition semantics come from R/analysis_labels.R.

protigy_supported_metrics <- function() {
  c(
    "signed.Log.P.Value", "Log.P.Value", "adj.P.Val", "P.Value",
    "RawAveExpr", "RawlogFC", "significant", "sign.logP", "AveExpr",
    "logFC", "t", "B"
  )
}

protigy_required_da_metrics <- function() {
  c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "significant", "Log.P.Value")
}

protigy_metric_output_names <- function() {
  c(
    "adj.P.Val" = "padj",
    "P.Value" = "pval",
    "logFC" = "log2fc",
    "RawlogFC" = "rawlog2fc",
    "Log.P.Value" = "logpval",
    "signed.Log.P.Value" = "signed.logpval",
    "AveExpr" = "aveExpr",
    "RawAveExpr" = "rawAveExpr",
    "t" = "t"
  )
}

canonicalize_protigy_comparison_separator <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  gsub(".over.", "_over_", x, fixed = TRUE)
}

protigy_comparison_naming_style <- function(x) {
  x <- as.character(x)
  out <- rep("unrecognized", length(x))
  out[!is.na(x) & grepl(".over.", x, fixed = TRUE)] <- "historical_dot_over"
  out[!is.na(x) & grepl("_over_", x, fixed = TRUE)] <- "corrected_underscore_over"
  out
}

parse_protigy_stat_field <- function(field, metrics = protigy_supported_metrics()) {
  field <- as.character(field)
  metrics <- metrics[order(nchar(metrics), decreasing = TRUE, method = "radix")]
  metric <- raw_comparison <- rep(NA_character_, length(field))
  naming_style <- rep("unrecognized", length(field))

  for (candidate in metrics) {
    prefix <- paste0(candidate, ".")
    hit <- is.na(metric) & !is.na(field) & startsWith(field, prefix)
    if (!any(hit)) next
    suffix <- substring(field[hit], nchar(prefix) + 1L)
    style <- protigy_comparison_naming_style(suffix)
    parsed <- style != "unrecognized"
    index <- which(hit)
    metric[index[parsed]] <- candidate
    raw_comparison[index[parsed]] <- suffix[parsed]
    naming_style[index[parsed]] <- style[parsed]
  }

  data.frame(
    field = field,
    metric = metric,
    raw_comparison = raw_comparison,
    comparison = canonicalize_protigy_comparison_separator(raw_comparison),
    naming_style = naming_style,
    stringsAsFactors = FALSE
  )
}

parse_neha_protigy_side <- function(side) {
  side <- trimws(as.character(side))
  if (length(side) != 1L || is.na(side) || !nzchar(side)) return(NULL)

  conditions <- condition_levels[order(nchar(condition_levels), decreasing = TRUE, method = "radix")]
  for (condition in conditions) {
    suffix <- paste0("_", condition)
    if (!endsWith(tolower(side), suffix)) next
    class_token <- substring(side, 1L, nchar(side) - nchar(suffix))
    sample_class <- normalize_sample_class(class_token)
    if (length(sample_class) == 1L && !is.na(sample_class)) {
      return(list(
        original = side,
        sample_class = sample_class,
        condition = condition,
        condition_code = unname(names(condition_code_map)[match(condition, condition_code_map)]),
        phenotype = paste(sample_class, condition, sep = "_")
      ))
    }
  }

  match <- regexec("^(.+?)[_]?([1-4])$", side, perl = TRUE)
  parts <- regmatches(side, match)[[1]]
  if (!length(parts)) return(NULL)
  sample_class <- normalize_sample_class(parts[[2]])
  if (length(sample_class) != 1L || is.na(sample_class)) return(NULL)
  condition_code <- parts[[3]]
  condition <- unname(condition_code_map[[condition_code]])
  list(
    original = side,
    sample_class = sample_class,
    condition = condition,
    condition_code = condition_code,
    phenotype = paste(sample_class, condition, sep = "_")
  )
}

parse_neha_protigy_comparison_one <- function(comparison) {
  input <- canonicalize_protigy_comparison_separator(comparison)
  empty <- data.frame(
    input_comparison = input,
    canonical_comparison = NA_character_,
    parsed = FALSE,
    sample_class = NA_character_,
    numerator_sample_class = NA_character_,
    denominator_sample_class = NA_character_,
    numerator_condition = NA_character_,
    denominator_condition = NA_character_,
    same_sample_class = FALSE,
    accepted_primary = FALSE,
    canonical_contrast = NA_character_,
    historical_comparison_name = NA_character_,
    rejection_reason = "unparseable_comparison",
    stringsAsFactors = FALSE
  )
  if (length(input) != 1L || is.na(input)) return(empty)

  sides <- strsplit(input, "_over_", fixed = TRUE)[[1]]
  if (length(sides) != 2L || any(!nzchar(sides))) return(empty)
  numerator <- parse_neha_protigy_side(sides[[1]])
  denominator <- parse_neha_protigy_side(sides[[2]])
  if (is.null(numerator) || is.null(denominator)) return(empty)

  same_class <- identical(numerator$sample_class, denominator$sample_class)
  canonical <- paste(numerator$phenotype, denominator$phenotype, sep = "_over_")
  expected <- neha_primary_contrast_manifest()
  expected_index <- match(canonical, expected$canonical_comparison)
  accepted <- same_class && !is.na(expected_index)
  reason <- if (!same_class) {
    "cross_sample_class"
  } else if (!accepted) {
    "unsupported_condition_contrast"
  } else {
    "accepted_primary_neha_contrast"
  }

  data.frame(
    input_comparison = input,
    canonical_comparison = canonical,
    parsed = TRUE,
    sample_class = if (same_class) numerator$sample_class else NA_character_,
    numerator_sample_class = numerator$sample_class,
    denominator_sample_class = denominator$sample_class,
    numerator_condition = numerator$condition,
    denominator_condition = denominator$condition,
    same_sample_class = same_class,
    accepted_primary = accepted,
    canonical_contrast = if (accepted) expected$canonical_contrast[[expected_index]] else NA_character_,
    historical_comparison_name = if (accepted) expected$historical_comparison_name[[expected_index]] else NA_character_,
    rejection_reason = reason,
    stringsAsFactors = FALSE
  )
}

parse_neha_protigy_comparison <- function(comparison) {
  rows <- lapply(as.character(comparison), parse_neha_protigy_comparison_one)
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

validate_neha_protigy_comparisons <- function(comparisons, strict_primary = FALSE) {
  out <- parse_neha_protigy_comparison(comparisons)
  if (nrow(out)) {
    out$duplicate_canonical_comparison <- duplicated(out$canonical_comparison) |
      duplicated(out$canonical_comparison, fromLast = TRUE)
    rownames(out) <- NULL
  }
  if (isTRUE(strict_primary) && nrow(out) && any(!out$accepted_primary)) {
    rejected <- paste0(
      out$input_comparison[!out$accepted_primary],
      " [", out$rejection_reason[!out$accepted_primary], "]"
    )
    stop(
      "Animal-level ProTigy GCT contains rejected comparison(s): ",
      paste(rejected, collapse = ", "),
      call. = FALSE
    )
  }
  out
}

validate_neha_primary_comparison_contract <- function(comparisons) {
  expected <- neha_primary_contrast_manifest()
  observed <- validate_neha_protigy_comparisons(comparisons, strict_primary = TRUE)
  observed_keys <- observed$canonical_comparison
  failures <- character()
  if (length(observed_keys) != nrow(expected)) failures <- c(failures, "expected exactly 12 comparisons")
  if (anyDuplicated(observed_keys)) failures <- c(failures, "canonical comparisons are duplicated")
  missing <- setdiff(expected$canonical_comparison, observed_keys)
  unexpected <- setdiff(observed_keys, expected$canonical_comparison)
  if (length(missing)) failures <- c(failures, paste0("missing: ", paste(missing, collapse = ", ")))
  if (length(unexpected)) failures <- c(failures, paste0("unexpected: ", paste(unexpected, collapse = ", ")))
  if (length(failures)) {
    stop("Neha primary comparison contract failed: ", paste(failures, collapse = "; "), call. = FALSE)
  }
  observed[match(expected$canonical_comparison, observed$canonical_comparison), , drop = FALSE]
}

split_gct_tab_line <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (endsWith(line, "\t")) c(fields, "") else fields
}

protigy_file_sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required to calculate source GCT SHA-256.", call. = FALSE)
  }
  tolower(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

read_protigy_stat_gct <- function(path, strict_primary = FALSE) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 3L || trimws(sub("^\ufeff", "", lines[[1]])) != "#1.3") {
    stop("Expected a GCT v1.3 file: ", path, call. = FALSE)
  }
  dims <- suppressWarnings(scan(text = lines[[2]], what = integer(), quiet = TRUE))
  if (length(dims) != 4L || anyNA(dims) || any(dims < 0L)) {
    stop("Invalid GCT v1.3 dimension line: ", path, call. = FALSE)
  }
  names(dims) <- c("nrmat", "ncmat", "nrhd", "nchd")
  expected_fields <- 1L + dims[["nrhd"]] + dims[["ncmat"]]

  first_fields <- vapply(lines[-seq_len(2L)], function(line) {
    trimws(sub("\t.*$", "", line))
  }, character(1))
  header_candidates <- which(tolower(first_fields) == "id") + 2L
  if (!length(header_candidates)) stop("Could not find the GCT header row in: ", path, call. = FALSE)
  header_row <- header_candidates[[1]]
  original_names <- split_gct_tab_line(lines[[header_row]])
  if (length(original_names) != expected_fields) {
    stop(
      "GCT physical field count does not match dimensions: expected ", expected_fields,
      ", observed ", length(original_names), ".",
      call. = FALSE
    )
  }
  id_index <- which(tolower(trimws(original_names)) == "id")
  if (length(id_index) != 1L || id_index[[1]] != 1L) {
    stop("GCT header must contain exactly one reserved id field in the first position.", call. = FALSE)
  }

  data_start <- header_row + dims[["nchd"]] + 1L
  data_end <- data_start + dims[["nrmat"]] - 1L
  if (data_end > length(lines)) stop("GCT contains fewer protein rows than declared.", call. = FALSE)
  data_fields <- lapply(lines[seq.int(data_start, data_end)], split_gct_tab_line)
  field_counts <- lengths(data_fields)
  if (any(field_counts != expected_fields)) {
    bad <- which(field_counts != expected_fields)
    stop(
      "GCT protein row field count is invalid at declared protein row(s): ",
      paste(head(bad, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  data_matrix <- do.call(rbind, data_fields)
  internal_names <- make.unique(
    ifelse(nzchar(original_names), original_names, "__blank__"),
    sep = "__duplicate_"
  )
  colnames(data_matrix) <- internal_names
  data <- as.data.frame(data_matrix, stringsAsFactors = FALSE, check.names = FALSE)
  ids <- trimws(as.character(data[[internal_names[[1]]]]))

  description_index <- which(tolower(trimws(original_names)) == "description")
  description <- if (length(description_index) == 1L) {
    as.character(data[[internal_names[[description_index[[1]]]]]])
  } else {
    rep(NA_character_, length(ids))
  }

  fields <- parse_protigy_stat_field(original_names)
  fields$column_index <- seq_along(original_names)
  fields$column_internal <- internal_names
  fields$physical_role <- ifelse(
    fields$column_index <= (1L + dims[["nrhd"]]),
    "row_descriptor",
    "matrix_column"
  )
  parsed_fields <- fields[!is.na(fields$comparison), , drop = FALSE]
  comparisons <- unique(parsed_fields$comparison)
  comparison_validation <- validate_neha_protigy_comparisons(comparisons, strict_primary = strict_primary)
  canonical_match <- match(parsed_fields$comparison, comparison_validation$input_comparison)
  parsed_fields$canonical_comparison <- comparison_validation$canonical_comparison[canonical_match]
  parsed_fields$metric_comparison_key <- paste(
    parsed_fields$metric,
    parsed_fields$canonical_comparison,
    sep = "||"
  )
  parsed_fields$duplicate_metric_comparison <- duplicated(parsed_fields$metric_comparison_key) |
    duplicated(parsed_fields$metric_comparison_key, fromLast = TRUE)
  styles <- sort(unique(parsed_fields$naming_style), method = "radix")

  list(
    path = path,
    sha256 = protigy_file_sha256(path),
    dimensions = dims,
    header_row = header_row,
    data_start = data_start,
    expected_physical_fields = expected_fields,
    observed_physical_fields = length(original_names),
    physical_dimension_match = identical(as.integer(expected_fields), as.integer(length(original_names))),
    data = data,
    original_names = original_names,
    internal_names = internal_names,
    ids = ids,
    description = description,
    fields = fields,
    parsed_fields = parsed_fields,
    comparison_validation = comparison_validation,
    metrics_found = sort(unique(parsed_fields$metric), method = "radix"),
    naming_style = if (!length(styles)) "none" else if (length(styles) == 1L) styles else paste0("mixed:", paste(styles, collapse = "+")),
    protein_ids_unique = !anyNA(ids) && all(nzchar(ids)) && !anyDuplicated(ids),
    description_present = length(description_index) == 1L,
    description_aligned = length(description) == length(ids),
    description_complete = length(description) == length(ids) && all(!is.na(description) & nzchar(trimws(description))),
    duplicate_metric_comparison_fields = sum(parsed_fields$duplicate_metric_comparison),
    n_protein_rows_read = nrow(data)
  )
}

parse_protigy_metric_values <- function(values, metric, field = metric) {
  values <- trimws(as.character(values))
  missing <- is.na(values) | !nzchar(values) | toupper(values) %in% c("NA", "NAN")
  if (identical(metric, "significant")) {
    normalized <- toupper(values)
    invalid <- !missing & !normalized %in% c("TRUE", "FALSE")
    if (any(invalid)) {
      stop("Logical ProTigy field contains values other than TRUE/FALSE: ", field, call. = FALSE)
    }
    out <- rep(NA, length(values))
    out[!missing] <- normalized[!missing] == "TRUE"
    return(out)
  }

  out <- suppressWarnings(as.numeric(values))
  invalid <- !missing & is.na(out)
  if (any(invalid)) stop("Numeric ProTigy field contains non-numeric values: ", field, call. = FALSE)
  out[missing] <- NA_real_
  out
}

protigy_metric_is_signed <- function(metric) {
  as.character(metric) %in% c("logFC", "RawlogFC", "t")
}

reverse_protigy_metric_frame <- function(df, metric_by_column) {
  out <- df
  shared <- intersect(names(metric_by_column), names(out))
  for (column in shared) {
    if (protigy_metric_is_signed(metric_by_column[[column]])) {
      out[[column]] <- -suppressWarnings(as.numeric(out[[column]]))
    }
  }
  out
}

extract_protigy_comparison_table <- function(gct, comparison) {
  comparison <- parse_neha_protigy_comparison_one(comparison)$canonical_comparison[[1]]
  if (is.na(comparison)) stop("Cannot extract an unparseable comparison.", call. = FALSE)
  selected <- gct$parsed_fields[
    gct$parsed_fields$canonical_comparison == comparison,
    ,
    drop = FALSE
  ]
  selected <- selected[order(selected$column_index, method = "radix"), , drop = FALSE]
  if (!nrow(selected)) stop("Comparison not found in GCT: ", comparison, call. = FALSE)
  if (any(selected$duplicate_metric_comparison)) {
    stop("Duplicate metric/comparison fields prevent extraction: ", comparison, call. = FALSE)
  }

  output_map <- protigy_metric_output_names()
  output_names <- vapply(selected$metric, function(metric) {
    if (metric %in% names(output_map)) unname(output_map[[metric]]) else metric
  }, character(1))
  if (anyDuplicated(output_names)) {
    stop("Metric output names collide for comparison: ", comparison, call. = FALSE)
  }
  out <- data.frame(
    # Compatibility contract: gene_symbol is the original GCT id/UniProt-style
    # protein identifier. It is historically misnamed for MapThatProt input.
    gene_symbol = as.character(gct$ids),
    Description = as.character(gct$description),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (i in seq_len(nrow(selected))) {
    out[[output_names[[i]]]] <- parse_protigy_metric_values(
      gct$data[[selected$column_internal[[i]]]],
      selected$metric[[i]],
      selected$field[[i]]
    )
  }
  attr(out, "metric_by_column") <- stats::setNames(selected$metric, output_names)
  out
}

validate_neha_stat_gct_contract <- function(gct, expected_n_proteins = 5349L) {
  failures <- character()
  if (!gct$physical_dimension_match) failures <- c(failures, "physical GCT dimensions do not match")
  if (gct$n_protein_rows_read != unname(gct$dimensions[["nrmat"]])) failures <- c(failures, "declared protein rows were not all read")
  if (gct$n_protein_rows_read != as.integer(expected_n_proteins)) failures <- c(failures, paste0("expected ", expected_n_proteins, " proteins"))
  if (!gct$protein_ids_unique) failures <- c(failures, "protein IDs are missing or duplicated")
  if (!gct$description_present || !gct$description_aligned || !gct$description_complete) {
    failures <- c(failures, "Description is missing, incomplete, or misaligned")
  }
  if (gct$duplicate_metric_comparison_fields > 0L) failures <- c(failures, "metric/comparison fields are duplicated")

  comparison_contract <- tryCatch(
    validate_neha_primary_comparison_contract(gct$comparison_validation$input_comparison),
    error = function(error) {
      failures <<- c(failures, conditionMessage(error))
      NULL
    }
  )
  expected <- neha_primary_contrast_manifest()
  required <- protigy_required_da_metrics()
  required_grid <- expand.grid(
    metric = required,
    canonical_comparison = expected$canonical_comparison,
    stringsAsFactors = FALSE
  )
  present_keys <- gct$parsed_fields$metric_comparison_key
  required_grid$key <- paste(required_grid$metric, required_grid$canonical_comparison, sep = "||")
  counts <- table(factor(present_keys, levels = required_grid$key))
  missing <- required_grid[counts == 0L, c("metric", "canonical_comparison"), drop = FALSE]
  duplicate_required <- required_grid[counts > 1L, c("metric", "canonical_comparison"), drop = FALSE]
  if (nrow(missing)) failures <- c(failures, "required DA metric/comparison fields are missing")
  if (nrow(duplicate_required)) failures <- c(failures, "required DA metric/comparison fields are duplicated")

  value_validation <- list()
  if (!nrow(missing) && !nrow(duplicate_required)) {
    for (i in seq_len(nrow(required_grid))) {
      field <- gct$parsed_fields[gct$parsed_fields$metric_comparison_key == required_grid$key[[i]], , drop = FALSE]
      values <- tryCatch(
        parse_protigy_metric_values(gct$data[[field$column_internal[[1]]]], field$metric[[1]], field$field[[1]]),
        error = function(error) {
          failures <<- c(failures, conditionMessage(error))
          NULL
        }
      )
      value_validation[[i]] <- data.frame(
        metric = required_grid$metric[[i]],
        canonical_comparison = required_grid$canonical_comparison[[i]],
        n_values = if (is.null(values)) 0L else length(values),
        n_missing = if (is.null(values)) NA_integer_ else sum(is.na(values)),
        stringsAsFactors = FALSE
      )
    }
  }
  value_validation <- if (length(value_validation)) do.call(rbind, value_validation) else data.frame()
  if (nrow(value_validation) && any(value_validation$n_missing > 0L, na.rm = TRUE)) {
    failures <- c(failures, "required DA fields contain missing values")
  }
  if (length(failures)) {
    stop("Neha ProTigy statistical-GCT contract failed: ", paste(unique(failures), collapse = "; "), call. = FALSE)
  }

  list(
    status = "PASS",
    expected_n_proteins = as.integer(expected_n_proteins),
    observed_n_proteins = gct$n_protein_rows_read,
    expected_comparison_count = nrow(expected),
    observed_comparison_count = nrow(comparison_contract),
    comparisons = comparison_contract,
    required_metrics = required,
    value_validation = value_validation
  )
}
