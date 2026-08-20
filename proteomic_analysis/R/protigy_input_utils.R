read_gct_v13 <- function(path) {
  if (!file.exists(path)) stop("GCT file does not exist: ", path, call. = FALSE)
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 3L || trimws(lines[1]) != "#1.3") {
    stop("Expected a GCT v1.3 file: ", path, call. = FALSE)
  }
  dims <- scan(text = lines[2], what = integer(), quiet = TRUE)
  if (length(dims) != 4L) {
    stop("GCT v1.3 dimension line must contain four integers.", call. = FALSE)
  }
  n_rows <- dims[1]
  n_samples <- dims[2]
  n_row_descriptors <- dims[3]
  n_column_metadata <- dims[4]
  header <- strsplit(lines[3], "\t", fixed = TRUE)[[1]]
  expected_fields <- 1L + n_row_descriptors + n_samples
  if (length(header) != expected_fields || header[1] != "id") {
    stop("GCT header field count or reserved id field is invalid.", call. = FALSE)
  }

  metadata_lines <- if (n_column_metadata > 0L) lines[seq.int(4L, 3L + n_column_metadata)] else character()
  column_metadata_fillers <- character()
  column_metadata <- if (length(metadata_lines)) {
    fields <- strsplit(metadata_lines, "\t", fixed = TRUE)
    if (any(lengths(fields) != expected_fields)) stop("GCT column metadata field count is invalid.", call. = FALSE)
    mat <- do.call(rbind, fields)
    rownames(mat) <- mat[, 1]
    if (n_row_descriptors > 0L) {
      column_metadata_fillers <- mat[, seq.int(2L, 1L + n_row_descriptors), drop = FALSE]
      colnames(column_metadata_fillers) <- header[seq.int(2L, 1L + n_row_descriptors)]
    }
    mat[, -(seq_len(1L + n_row_descriptors)), drop = FALSE]
  } else {
    matrix(character(), nrow = 0L, ncol = n_samples)
  }
  colnames(column_metadata) <- tail(header, n_samples)

  data_start <- 4L + n_column_metadata
  data_end <- data_start + n_rows - 1L
  if (data_end > length(lines)) stop("GCT contains fewer protein rows than declared.", call. = FALSE)
  data_fields <- strsplit(lines[seq.int(data_start, data_end)], "\t", fixed = TRUE)
  if (any(lengths(data_fields) != expected_fields)) stop("GCT protein row field count is invalid.", call. = FALSE)
  data_mat <- do.call(rbind, data_fields)
  protein_ids <- data_mat[, 1]
  if (anyDuplicated(protein_ids)) stop("GCT protein IDs are duplicated.", call. = FALSE)
  descriptor_names <- if (n_row_descriptors) header[seq.int(2L, 1L + n_row_descriptors)] else character()
  row_descriptors <- if (n_row_descriptors) {
    out <- as.data.frame(data_mat[, seq.int(2L, 1L + n_row_descriptors), drop = FALSE], stringsAsFactors = FALSE)
    names(out) <- descriptor_names
    rownames(out) <- protein_ids
    out
  } else {
    out <- data.frame(row.names = protein_ids)
    out
  }
  values <- data_mat[, seq.int(2L + n_row_descriptors, expected_fields), drop = FALSE]
  storage.mode(values) <- "numeric"
  rownames(values) <- protein_ids
  colnames(values) <- tail(header, n_samples)

  list(
    version = "#1.3",
    dimensions = c(
      n_rows = n_rows,
      n_samples = n_samples,
      n_row_descriptors = n_row_descriptors,
      n_column_metadata = n_column_metadata
    ),
    matrix = values,
    row_descriptors = row_descriptors,
    column_metadata = column_metadata,
    column_metadata_fillers = column_metadata_fillers,
    protein_ids = protein_ids
  )
}

format_gct_number <- function(x) {
  ifelse(is.na(x), "NaN", format(x, digits = 17L, scientific = FALSE, trim = TRUE))
}

write_protigy_gct_v13 <- function(expression_matrix, sample_metadata, descriptions, path) {
  expression_matrix <- as.matrix(expression_matrix)
  storage.mode(expression_matrix) <- "numeric"
  protein_ids <- rownames(expression_matrix)
  sample_ids <- colnames(expression_matrix)
  if (is.null(protein_ids) || anyDuplicated(protein_ids) || any(is.na(protein_ids) | !nzchar(protein_ids))) {
    stop("Expression matrix must have complete, unique protein row names.", call. = FALSE)
  }
  if (is.null(sample_ids) || anyDuplicated(sample_ids) || any(is.na(sample_ids) | !nzchar(sample_ids))) {
    stop("Expression matrix must have complete, unique sample column names.", call. = FALSE)
  }
  descriptions <- as.character(descriptions)
  if (length(descriptions) != nrow(expression_matrix) || any(is.na(descriptions) | !nzchar(descriptions))) {
    stop("Descriptions must be complete and one-to-one with protein rows.", call. = FALSE)
  }
  sample_metadata <- as.data.frame(sample_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  required_meta <- c("output_column_name", "AnimalID", "condition_code", "condition", "sample_class")
  missing_meta <- setdiff(required_meta, names(sample_metadata))
  if (length(missing_meta) > 0) stop("Sample metadata missing: ", paste(missing_meta, collapse = ", "), call. = FALSE)
  if (anyDuplicated(sample_metadata$output_column_name)) stop("Sample metadata output IDs are duplicated.", call. = FALSE)
  idx <- match(sample_ids, sample_metadata$output_column_name)
  if (anyNA(idx)) stop("Sample metadata does not cover every expression column.", call. = FALSE)
  sample_metadata <- sample_metadata[idx, , drop = FALSE]
  phenotype <- paste(sample_metadata$sample_class, sample_metadata$condition, sep = "_")
  meta_values <- list(
    AnimalID = sample_metadata$AnimalID,
    condition_code = sample_metadata$condition_code,
    condition = sample_metadata$condition,
    sample_class = sample_metadata$sample_class,
    phenotypeWithinUnit = phenotype
  )

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines("#1.3", con)
  writeLines(sprintf("%d\t%d\t1\t%d", nrow(expression_matrix), ncol(expression_matrix), length(meta_values)), con)
  writeLines(paste(c("id", "Description", sample_ids), collapse = "\t"), con)
  for (meta_name in names(meta_values)) {
    writeLines(paste(c(meta_name, "na", as.character(meta_values[[meta_name]])), collapse = "\t"), con)
  }
  for (i in seq_len(nrow(expression_matrix))) {
    writeLines(
      paste(c(protein_ids[i], descriptions[i], format_gct_number(expression_matrix[i, ])), collapse = "\t"),
      con
    )
  }
  invisible(path)
}

validate_protigy_gct_v13 <- function(path, expected_matrix = NULL) {
  parsed <- read_gct_v13(path)
  if (unname(parsed$dimensions["n_row_descriptors"]) != 1L) {
    stop("ProTigy GCT must contain exactly one row descriptor.", call. = FALSE)
  }
  if (!identical(names(parsed$row_descriptors), "Description")) {
    stop("ProTigy GCT row descriptor must be named Description.", call. = FALSE)
  }
  required_meta <- c("AnimalID", "condition_code", "condition", "sample_class", "phenotypeWithinUnit")
  if (!identical(rownames(parsed$column_metadata), required_meta)) {
    stop("ProTigy GCT column metadata rows are missing or out of order.", call. = FALSE)
  }
  if (!is.matrix(parsed$column_metadata_fillers) ||
      ncol(parsed$column_metadata_fillers) != 1L ||
      any(parsed$column_metadata_fillers[, 1] != "na")) {
    stop("ProTigy GCT column metadata filler must be na for every row.", call. = FALSE)
  }
  if (!identical(rownames(parsed$matrix), rownames(parsed$row_descriptors))) {
    stop("GCT matrix and row descriptor row names are not aligned.", call. = FALSE)
  }
  if (!is.null(expected_matrix)) {
    expected_matrix <- as.matrix(expected_matrix)
    if (!identical(dim(parsed$matrix), dim(expected_matrix)) ||
        !identical(rownames(parsed$matrix), rownames(expected_matrix)) ||
        !identical(colnames(parsed$matrix), colnames(expected_matrix)) ||
        !isTRUE(all.equal(unname(parsed$matrix), unname(expected_matrix), tolerance = 0, check.attributes = FALSE))) {
      stop("GCT numeric matrix does not exactly match the expected matrix.", call. = FALSE)
    }
  }
  parsed
}
