# Input, output, identity, and provenance contracts for the canonical Neha
# animal-level MapThatProt branch. The protein-resolution cascade remains in
# 02_id_mapping/01_MapThatProt_batch.r.

mapthatprot_is_absolute_path <- function(path) {
  grepl("^([A-Za-z]:[\\\\/]|\\\\\\\\|/)", path)
}

normalize_mapthatprot_path <- function(path, must_work = FALSE) {
  normalizePath(path.expand(trimws(as.character(path))), winslash = "/", mustWork = must_work)
}

mapthatprot_paths_equal <- function(left, right) {
  left <- normalize_mapthatprot_path(left)
  right <- normalize_mapthatprot_path(right)
  if (.Platform$OS.type == "windows") {
    left <- tolower(left)
    right <- tolower(right)
  }
  identical(left, right)
}

mapthatprot_path_is_within <- function(path, parent) {
  path <- normalize_mapthatprot_path(path)
  parent <- sub("/+$", "", normalize_mapthatprot_path(parent))
  if (.Platform$OS.type == "windows") {
    path <- tolower(path)
    parent <- tolower(parent)
  }
  identical(path, parent) || startsWith(path, paste0(parent, "/"))
}

neha_mapthatprot_default_paths <- function() {
  project_root <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
  dataset_root <- file.path(project_root, "Datasets")
  list(
    project_root = project_root,
    split_root = file.path(dataset_root, "data", "protigy_animal_level", "split"),
    output_root = file.path(dataset_root, "data", "protigy_animal_level", "mapped"),
    uniprot_mapping_file = file.path(dataset_root, "MOUSE_10090_idmapping.dat"),
    manual_mapping_file = file.path(dataset_root, "manual_mapping.xlsx"),
    historical_roots = file.path(dataset_root, c("raw", "mapped", "unmapped", "mapping_reports"))
  )
}

validate_neha_mapthatprot_output_root <- function(output_root, split_root, historical_roots) {
  output_root <- normalize_mapthatprot_path(output_root)
  split_root <- normalize_mapthatprot_path(split_root)
  historical_roots <- vapply(historical_roots, normalize_mapthatprot_path, character(1))
  protected <- c(split_root, historical_roots)
  hit <- vapply(protected, function(path) mapthatprot_path_is_within(output_root, path), logical(1))
  if (any(hit)) {
    stop(
      "Canonical animal-level mapping output root is protected or historical: ",
      output_root,
      " (inside ", protected[[which(hit)[[1]]]], ")",
      call. = FALSE
    )
  }
  output_root
}

resolve_neha_mapthatprot_config <- function(
    split_root = Sys.getenv("NEHA_MAPTHATPROT_SPLIT_ROOT", unset = ""),
    output_root = Sys.getenv("NEHA_MAPTHATPROT_OUTPUT_ROOT", unset = ""),
    direction = Sys.getenv("NEHA_MAPTHATPROT_DIRECTION", unset = "forward"),
    uniprot_mapping_file = Sys.getenv("NEHA_MAPTHATPROT_REFERENCE_FILE", unset = ""),
    manual_mapping_file = Sys.getenv("NEHA_MAPTHATPROT_MANUAL_MAPPING_FILE", unset = "")) {
  defaults <- neha_mapthatprot_default_paths()
  choose <- function(value, default) if (nzchar(trimws(value))) value else default
  direction <- tolower(trimws(direction))
  if (!direction %in% c("forward", "reverse")) {
    stop("NEHA_MAPTHATPROT_DIRECTION must be 'forward' or 'reverse'.", call. = FALSE)
  }
  split_root <- normalize_mapthatprot_path(choose(split_root, defaults$split_root))
  output_root <- validate_neha_mapthatprot_output_root(
    choose(output_root, defaults$output_root),
    split_root,
    defaults$historical_roots
  )
  list(
    direction = direction,
    split_root = split_root,
    output_root = output_root,
    project_root = normalize_mapthatprot_path(defaults$project_root),
    uniprot_mapping_file = normalize_mapthatprot_path(choose(uniprot_mapping_file, defaults$uniprot_mapping_file)),
    manual_mapping_file = normalize_mapthatprot_path(choose(manual_mapping_file, defaults$manual_mapping_file)),
    historical_roots = vapply(defaults$historical_roots, normalize_mapthatprot_path, character(1))
  )
}

neha_mapthatprot_required_index_columns <- function() {
  c(
    "canonical_comparison", "canonical_contrast", "sample_class",
    "numerator_condition", "denominator_condition",
    "historical_comparison_alias", "forward_filename", "forward_path",
    "reverse_filename", "reverse_path", "n_proteins", "source_gct_sha256"
  )
}

read_neha_mapthatprot_input_index <- function(
    split_root,
    direction = "forward",
    expected_manifest = neha_primary_contrast_manifest(),
    expected_n_proteins = 5349L) {
  direction <- match.arg(direction, c("forward", "reverse"))
  split_root <- normalize_mapthatprot_path(split_root, must_work = TRUE)
  index_path <- file.path(split_root, "indexComparisons.csv")
  if (!file.exists(index_path)) stop("Missing authoritative split index: ", index_path, call. = FALSE)
  index <- utils::read.csv(index_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- neha_mapthatprot_required_index_columns()
  missing_columns <- setdiff(required, names(index))
  if (length(missing_columns)) {
    stop("Split index is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  if (nrow(index) != nrow(expected_manifest) || nrow(index) != 12L) {
    stop("Canonical animal-level mapping requires exactly 12 indexed comparisons; found ", nrow(index), ".", call. = FALSE)
  }
  duplicate_fields <- c("canonical_comparison", "canonical_contrast", paste0(direction, "_filename"), paste0(direction, "_path"))
  duplicated_fields <- duplicate_fields[vapply(duplicate_fields, function(field) anyDuplicated(index[[field]]) > 0L, logical(1))]
  if (length(duplicated_fields)) {
    stop("Split index contains duplicate values in: ", paste(duplicated_fields, collapse = ", "), call. = FALSE)
  }
  expected_keys <- expected_manifest$canonical_comparison
  if (!setequal(index$canonical_comparison, expected_keys)) {
    stop("Split index comparison coverage does not match the 12 primary Neha contrasts.", call. = FALSE)
  }
  index <- index[match(expected_keys, index$canonical_comparison), , drop = FALSE]
  expected_contrasts <- expected_manifest$canonical_contrast
  if (!identical(index$canonical_contrast, expected_contrasts)) {
    stop("Split index canonical contrast names do not match the shared Neha manifest.", call. = FALSE)
  }
  expected_metadata <- data.frame(
    sample_class = expected_manifest$sample_class,
    numerator_condition = expected_manifest$case_condition,
    denominator_condition = expected_manifest$reference_condition,
    historical_comparison_alias = expected_manifest$historical_comparison_name,
    stringsAsFactors = FALSE
  )
  for (field in names(expected_metadata)) {
    if (!identical(as.character(index[[field]]), expected_metadata[[field]])) {
      stop("Split index ", field, " does not match the shared Neha manifest.", call. = FALSE)
    }
  }
  protein_counts <- suppressWarnings(as.integer(index$n_proteins))
  if (anyNA(protein_counts) || any(protein_counts != as.integer(expected_n_proteins))) {
    stop(
      "Split index must declare ", expected_n_proteins,
      " proteins for every canonical comparison.",
      call. = FALSE
    )
  }

  path_column <- paste0(direction, "_path")
  filename_column <- paste0(direction, "_filename")
  input_paths <- as.character(index[[path_column]])
  input_filenames <- as.character(index[[filename_column]])
  missing_files <- input_paths[!file.exists(input_paths)]
  if (length(missing_files)) {
    stop("Indexed split file(s) are missing: ", paste(missing_files, collapse = ", "), call. = FALSE)
  }
  input_paths <- vapply(input_paths, normalize_mapthatprot_path, character(1), must_work = TRUE)
  expected_input_root <- file.path(split_root, direction)
  outside <- !vapply(input_paths, mapthatprot_path_is_within, logical(1), parent = expected_input_root)
  if (any(outside)) {
    stop("Indexed ", direction, " file resolves outside the authoritative split/", direction, " root: ", input_paths[[which(outside)[[1]]]], call. = FALSE)
  }
  if (!identical(basename(input_paths), input_filenames)) {
    stop("Indexed split filename/path pairs are inconsistent.", call. = FALSE)
  }
  index$mapping_direction <- direction
  index$input_filename <- input_filenames
  index$input_path <- input_paths
  index$output_filename <- input_filenames
  index$split_index_path <- normalize_mapthatprot_path(index_path, must_work = TRUE)
  rownames(index) <- NULL
  index
}

neha_mapthatprot_required_input_columns <- function() {
  c(
    "gene_symbol", "Description", "log2fc", "aveExpr", "t", "pval",
    "padj", "B", "significant", "logpval"
  )
}

validate_neha_mapthatprot_input_file <- function(path, expected_n_proteins = 5349L) {
  path <- normalize_mapthatprot_path(path, must_work = TRUE)
  input <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing_columns <- setdiff(neha_mapthatprot_required_input_columns(), names(input))
  if (length(missing_columns)) {
    stop("Indexed split input is missing required columns: ", paste(missing_columns, collapse = ", "), " [", path, "]", call. = FALSE)
  }
  if (nrow(input) != as.integer(expected_n_proteins)) {
    stop("Indexed split input row count mismatch: expected ", expected_n_proteins, ", found ", nrow(input), " [", path, "]", call. = FALSE)
  }
  ids <- trimws(as.character(input$gene_symbol))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Indexed split input protein identifiers must be complete and unique: ", path, call. = FALSE)
  }
  description <- trimws(as.character(input$Description))
  if (length(description) != length(ids) || anyNA(description) || any(!nzchar(description))) {
    stop("Indexed split input Description must be complete and aligned: ", path, call. = FALSE)
  }
  if (!is.logical(input$significant)) {
    stop("Indexed split input significant column must parse as logical TRUE/FALSE: ", path, call. = FALSE)
  }
  input
}

build_neha_mapthatprot_output_tables <- function(input, mapping_audit) {
  input <- as.data.frame(input, stringsAsFactors = FALSE, check.names = FALSE)
  mapping_audit <- as.data.frame(mapping_audit, stringsAsFactors = FALSE, check.names = FALSE)
  required_audit <- c(
    "source_row_id", "original_protein_id", "uniprot_accession",
    "mapped_gene_symbol", "mapping_strategy", "mapping_status", "multi_protein"
  )
  missing_audit <- setdiff(required_audit, names(mapping_audit))
  if (length(missing_audit)) stop("Mapping audit is missing: ", paste(missing_audit, collapse = ", "), call. = FALSE)
  expected_rows <- seq_len(nrow(input))
  if (
    nrow(mapping_audit) != nrow(input) || anyDuplicated(mapping_audit$source_row_id) ||
      !setequal(as.integer(mapping_audit$source_row_id), expected_rows)
  ) {
    stop("Mapping audit must contain each input source_row_id exactly once.", call. = FALSE)
  }
  mapping_audit <- mapping_audit[match(expected_rows, mapping_audit$source_row_id), , drop = FALSE]
  if (!identical(as.character(mapping_audit$original_protein_id), as.character(input$gene_symbol))) {
    stop("Mapping audit original_protein_id is not aligned with the input identity key.", call. = FALSE)
  }
  if (anyNA(mapping_audit$mapping_status) || any(!mapping_audit$mapping_status %in% c("mapped", "unmapped"))) {
    stop("Mapping audit status must be exactly 'mapped' or 'unmapped'.", call. = FALSE)
  }
  mapped_status <- mapping_audit$mapping_status == "mapped"
  if (any(mapped_status & (is.na(mapping_audit$uniprot_accession) | !nzchar(mapping_audit$uniprot_accession)))) {
    stop("Mapped rows require a non-empty UniProt accession.", call. = FALSE)
  }
  if (any(!mapped_status & !is.na(mapping_audit$uniprot_accession) & nzchar(mapping_audit$uniprot_accession))) {
    stop("Unmapped rows cannot carry a UniProt accession.", call. = FALSE)
  }

  statistic_columns <- setdiff(names(input), c("gene_symbol", "Description"))
  make_output <- function(rows, mapped) {
    data.frame(
      gene_symbol = if (mapped) mapping_audit$uniprot_accession[rows] else rep(NA_character_, length(rows)),
      original_protein_id = as.character(input$gene_symbol[rows]),
      Description = as.character(input$Description[rows]),
      mapped_gene_symbol = as.character(mapping_audit$mapped_gene_symbol[rows]),
      uniprot_accession = as.character(mapping_audit$uniprot_accession[rows]),
      mapping_status = as.character(mapping_audit$mapping_status[rows]),
      mapping_strategy = as.character(mapping_audit$mapping_strategy[rows]),
      source_row_id = as.integer(mapping_audit$source_row_id[rows]),
      multi_protein = as.logical(mapping_audit$multi_protein[rows]),
      input[rows, statistic_columns, drop = FALSE],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  mapped_rows <- which(mapped_status)
  unmapped_rows <- which(!mapped_status)
  mapped <- make_output(mapped_rows, mapped = TRUE)
  unmapped <- make_output(unmapped_rows, mapped = FALSE)
  if (nrow(mapped) + nrow(unmapped) != nrow(input)) {
    stop("Mapped and unmapped outputs do not account for every input row.", call. = FALSE)
  }
  list(mapped = mapped, unmapped = unmapped, audit = mapping_audit)
}

validate_neha_mapthatprot_partition <- function(input, mapped, unmapped) {
  source_rows <- c(mapped$source_row_id, unmapped$source_row_id)
  if (
    length(source_rows) != nrow(input) || anyDuplicated(source_rows) ||
      !setequal(as.integer(source_rows), seq_len(nrow(input)))
  ) {
    stop("Mapped/unmapped output partition has row loss or duplication.", call. = FALSE)
  }
  retained_statistics <- setdiff(names(input), c("gene_symbol", "Description"))
  missing_mapped <- setdiff(retained_statistics, names(mapped))
  missing_unmapped <- setdiff(retained_statistics, names(unmapped))
  if (length(missing_mapped) || length(missing_unmapped)) {
    stop("Mapped/unmapped output dropped source statistical columns.", call. = FALSE)
  }
  TRUE
}

mapthatprot_sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("The digest package is required for SHA-256 provenance.", call. = FALSE)
  tolower(digest::digest(normalize_mapthatprot_path(path, must_work = TRUE), algo = "sha256", file = TRUE, serialize = FALSE))
}

mapthatprot_reference_provenance <- function(
    path,
    version = Sys.getenv(
      "NEHA_MAPTHATPROT_REFERENCE_VERSION",
      unset = "not_encoded_in_local_idmapping_file"
    )) {
  path <- normalize_mapthatprot_path(path, must_work = TRUE)
  info <- file.info(path)
  data.frame(
    mapping_reference_file = basename(path),
    mapping_reference_path = path,
    mapping_reference_version = version,
    mapping_reference_snapshot_date_utc = format(info$mtime, "%Y-%m-%d", tz = "UTC"),
    mapping_reference_sha256 = mapthatprot_sha256(path),
    mapping_reference_bytes = as.numeric(info$size),
    mapping_reference_modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}
