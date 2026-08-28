args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("01_preprocessing", "02a_prepare_animal_level_protigy_input.r")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "animal_level_proteomics_utils.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "animal_level_proteomics_utils.R"))
source(file.path(repo_root, "R", "protigy_input_utils.R"))

required_packages <- c("readxl", "writexl", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Install required packages before running this stage: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

option_or_env <- function(option_name, env_name, default) {
  value <- getOption(option_name)
  if (!is.null(value) && nzchar(as.character(value))) return(as.character(value))
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(value)) return(value)
  default
}

data_root <- option_or_env(
  "proteomics.animal_level_data_root",
  "PROTEOMICS_ANIMAL_LEVEL_DATA_ROOT",
  "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/02_data/gct"
)
project_root <- option_or_env(
  "proteomics.project_data_root",
  "PROTEOMICS_PROJECT_DATA_ROOT",
  "S:/Lab_Member/Tobi/Experiments/Collabs/Neha"
)
quantitative_gct <- option_or_env(
  "proteomics.animal_level_quantitative_gct",
  "PROTEOMICS_ANIMAL_LEVEL_QUANTITATIVE_GCT",
  file.path(data_root, "pg.matrix_filtered_pcaAdjusted_unnormalized.gct")
)
description_xlsx <- option_or_env(
  "proteomics.animal_level_description_xlsx",
  "PROTEOMICS_ANIMAL_LEVEL_DESCRIPTION_XLSX",
  file.path(data_root, "imputed_data.xlsx")
)
sample_info_xlsx <- option_or_env(
  "proteomics.animal_level_sample_info",
  "PROTEOMICS_ANIMAL_LEVEL_SAMPLE_INFO",
  # hand-maintained metadata lives under 01_input/, not with derived data
  "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/01_input/metadata/sample_info.xlsx"
)
sample_annotation_xlsx <- option_or_env(
  "proteomics.animal_level_sample_annotation",
  "PROTEOMICS_ANIMAL_LEVEL_SAMPLE_ANNOTATION",
  file.path(project_root, "sample_annotation.xlsx")
)
output_dir <- option_or_env(
  "proteomics.animal_level_output_dir",
  "PROTEOMICS_ANIMAL_LEVEL_OUTPUT_DIR",
  "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/02_data/animal_level/input_gct"
)

inputs <- c(
  quantitative_gct = quantitative_gct,
  description_xlsx = description_xlsx,
  sample_info_xlsx = sample_info_xlsx,
  sample_annotation_xlsx = sample_annotation_xlsx
)
code_inputs <- c(
  preparation_stage = file.path(repo_root, "01_preprocessing", "02a_prepare_animal_level_protigy_input.r"),
  aggregation_helper = file.path(repo_root, "R", "animal_level_proteomics_utils.R"),
  gct_helper = file.path(repo_root, "R", "protigy_input_utils.R"),
  analysis_labels = file.path(repo_root, "R", "analysis_labels.R")
)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Required historical inputs are unavailable:\n",
    paste(paste0(names(missing_inputs), ": ", missing_inputs), collapse = "\n"),
    call. = FALSE
  )
}

parse_legacy_label_hemisphere <- function(label) {
  label <- trimws(as.character(label))
  left <- !is.na(label) & grepl("_left$", label, ignore.case = TRUE)
  right <- !is.na(label) & grepl("_right$", label, ignore.case = TRUE)
  invalid <- left == right
  if (any(invalid)) {
    stop(
      "Historical annotation label lacks one explicit _left/_right suffix: ",
      paste(unique(label[invalid]), collapse = ", "),
      call. = FALSE
    )
  }
  ifelse(left, "Left", "Right")
}

adapt_verified_historical_metadata <- function(sample_info, sample_annotation, gct_sample_ids) {
  required_info <- c("id", "sampleNumber", "plate", "group2", "AnimalID", "ReplicateGroup", "celltype", "ExpGroup")
  required_annotation <- c("Name", "id", "label", "ReplicateGroup", "AnimalID", "ExpGroup")
  if (length(setdiff(required_info, names(sample_info))) > 0L) {
    stop("sample_info.xlsx does not match the historical 2024 dataset schema.", call. = FALSE)
  }
  if (length(setdiff(required_annotation, names(sample_annotation))) > 0L) {
    stop("sample_annotation.xlsx does not match the historical 2024 dataset schema.", call. = FALSE)
  }
  if (anyDuplicated(sample_info$sampleNumber) || anyDuplicated(sample_annotation$id)) {
    stop("Historical N-number sample identifiers must be unique.", call. = FALSE)
  }
  if (!identical(as.character(sample_info$id), as.character(gct_sample_ids))) {
    stop("sample_info id order does not exactly match the historical quantitative GCT columns.", call. = FALSE)
  }

  annotation_idx <- match(sample_info$sampleNumber, sample_annotation$id)
  if (anyNA(annotation_idx)) stop("sample_annotation.xlsx does not cover every historical ProTigy sample.", call. = FALSE)
  annotation <- sample_annotation[annotation_idx, , drop = FALSE]
  raw_basename <- sub("^.*[\\\\/]", "", as.character(annotation$Name))
  raw_basename <- sub("\\.d$", "", raw_basename, ignore.case = TRUE)
  if (!identical(raw_basename, as.character(sample_info$id))) {
    stop("Original annotation file identity does not match the ProTigy sample IDs.", call. = FALSE)
  }
  if (!identical(as.character(annotation$AnimalID), as.character(sample_info$AnimalID)) ||
      !identical(as.character(annotation$ExpGroup), as.character(sample_info$ExpGroup))) {
    stop("Historical sample_info and original annotation disagree on AnimalID or condition.", call. = FALSE)
  }

  explicit_hemisphere <- parse_legacy_label_hemisphere(annotation$label)
  numeric_group <- suppressWarnings(as.numeric(as.character(sample_info$ReplicateGroup)))
  annotation_numeric_group <- suppressWarnings(as.numeric(as.character(annotation$ReplicateGroup)))
  if (any(is.na(numeric_group) | !numeric_group %in% c(1, 2)) ||
      any(is.na(annotation_numeric_group) | !annotation_numeric_group %in% c(1, 2))) {
    stop("Legacy ReplicateGroup must be exactly numeric 1/2 in both historical metadata files.", call. = FALSE)
  }
  numeric_evidence <- ifelse(numeric_group == 1, "Left", "Right")
  annotation_numeric_evidence <- ifelse(annotation_numeric_group == 1, "Left", "Right")
  if (!identical(explicit_hemisphere, numeric_evidence) ||
      !identical(explicit_hemisphere, annotation_numeric_evidence)) {
    stop("Explicit historical _left/_right labels disagree with legacy ReplicateGroup.", call. = FALSE)
  }

  sample_class <- normalize_sample_class(sample_info$celltype)
  condition_code <- as.character(sample_info$ExpGroup)
  condition <- normalize_condition(condition_code)
  if (anyNA(sample_class) || anyNA(condition)) {
    stop("Historical sample class or condition cannot be mapped to current canonical labels.", call. = FALSE)
  }

  data.frame(
    sample_id = as.character(sample_info$id),
    AnimalID = as.character(sample_info$AnimalID),
    ReplicateGroup = explicit_hemisphere,
    sample_class = sample_class,
    condition_code = condition_code,
    condition = condition,
    exclude = FALSE,
    hemisphere_evidence = explicit_hemisphere,
    historical_annotation_label = as.character(annotation$label),
    historical_ReplicateGroup = as.character(sample_info$ReplicateGroup),
    historical_sample_number = as.character(sample_info$sampleNumber),
    historical_plate = as.character(sample_info$plate),
    sample_id_lr_token_status = "unavailable_in_2024_instrument_sample_id",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

historical <- read_gct_v13(quantitative_gct)
if (unname(historical$dimensions["n_row_descriptors"]) != 0L) {
  stop("Historical quantitative GCT unexpectedly contains row descriptors; audit the source before continuing.", call. = FALSE)
}
description_table <- as.data.frame(
  readxl::read_excel(description_xlsx, .name_repair = "minimal"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!all(c("Genes", "Gene_Group") %in% names(description_table))) {
  stop("Description workbook must contain Genes and Gene_Group.", call. = FALSE)
}
if (!identical(as.character(description_table$Gene_Group), historical$protein_ids)) {
  stop("Description workbook changes the historical protein IDs or order.", call. = FALSE)
}
description_values <- as.matrix(description_table[, -(1:2), drop = FALSE])
storage.mode(description_values) <- "numeric"
if (!identical(dim(description_values), dim(historical$matrix)) ||
    !isTRUE(all.equal(unname(description_values), unname(historical$matrix), tolerance = 0, check.attributes = FALSE))) {
  stop("Description workbook quantitative values do not exactly match the historical ProTigy input.", call. = FALSE)
}
descriptions <- as.character(description_table$Genes)
description_fallback <- is.na(descriptions) | !nzchar(trimws(descriptions))
descriptions[description_fallback] <- historical$protein_ids[description_fallback]

sample_info <- as.data.frame(readxl::read_excel(sample_info_xlsx), stringsAsFactors = FALSE, check.names = FALSE)
sample_annotation <- as.data.frame(readxl::read_excel(sample_annotation_xlsx), stringsAsFactors = FALSE, check.names = FALSE)
metadata <- adapt_verified_historical_metadata(sample_info, sample_annotation, colnames(historical$matrix))

expected_units <- make_expected_animal_units(metadata, sample_class_levels = sample_classes)
aggregation_plan <- prepare_animal_level_aggregation(
  metadata = metadata,
  # LEGACY DATASET IDENTIFIER -- do not rebrand. This value is written verbatim into the
  # validated aggregation_audit.csv as the "dataset/project" column and is embedded in every
  # canonical_analysis_unit key (e.g. "Neha__C11__mcherry"). Changing it would alter a
  # validated output, so it is retained as a frozen dataset label rather than project branding.
  dataset_project = "Neha",
  expected_units = expected_units,
  sample_class_levels = sample_classes,
  hemisphere_crosscheck = "explicit_evidence",
  hemisphere_evidence_col = "hemisphere_evidence"
)
aggregated <- aggregate_animal_level_matrix(
  expression_matrix = historical$matrix,
  protein_ids = historical$protein_ids,
  aggregation_plan = aggregation_plan
)
aggregation_summary <- summarize_animal_level_design(aggregation_plan$audit)

if (any(aggregation_summary$n_unique_animals != 3L)) {
  bad <- aggregation_summary[aggregation_summary$n_unique_animals != 3L, , drop = FALSE]
  stop(
    "Animal-level design does not preserve biological n=3 for every sample class/condition:\n",
    paste(utils::capture.output(print(bad, row.names = FALSE)), collapse = "\n"),
    call. = FALSE
  )
}

primary_metadata <- aggregated$primary_audit
primary_metadata$phenotypeWithinUnit <- paste(primary_metadata$sample_class, primary_metadata$condition, sep = "_")
has_incomplete_pairs <- any(aggregation_plan$audit$hemisphere_status != "bilateral_complete")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
primary_gct <- file.path(output_dir, "neha_protigy_input_animal_level_primary.gct")
primary_xlsx <- file.path(output_dir, "neha_protigy_input_animal_level_primary.xlsx")
sensitivity_gct <- if (has_incomplete_pairs) {
  file.path(output_dir, "neha_protigy_input_animal_level_strict_complete_bilateral.gct")
} else {
  NA_character_
}

write_protigy_gct_v13(aggregated$primary, primary_metadata, descriptions, primary_gct)
primary_gct_validation <- validate_protigy_gct_v13(primary_gct, expected_matrix = aggregated$primary)
if (has_incomplete_pairs) {
  write_protigy_gct_v13(aggregated$sensitivity, aggregated$sensitivity_audit, descriptions, sensitivity_gct)
  sensitivity_gct_validation <- validate_protigy_gct_v13(sensitivity_gct, expected_matrix = aggregated$sensitivity)
} else if (!identical(aggregated$primary, aggregated$sensitivity)) {
  stop("All pairs are complete but primary and sensitivity matrices differ.", call. = FALSE)
}

primary_export <- data.frame(
  id = historical$protein_ids,
  Description = descriptions,
  aggregated$primary,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
writexl::write_xlsx(primary_export, primary_xlsx)

aggregation_audit_path <- file.path(output_dir, "aggregation_audit.csv")
aggregation_summary_path <- file.path(output_dir, "aggregation_summary.csv")
source_assignment_path <- file.path(output_dir, "source_sample_assignment.csv")
feature_audit_path <- file.path(output_dir, "feature_identity_audit.csv")
contrast_manifest_path <- file.path(output_dir, "protigy_contrast_manifest.csv")
contrast_scope_path <- file.path(output_dir, "contrast_scope_audit.csv")
metadata_contract_path <- file.path(output_dir, "metadata_contract_audit.csv")

utils::write.csv(aggregation_plan$audit, aggregation_audit_path, row.names = FALSE, na = "")
utils::write.csv(aggregation_summary, aggregation_summary_path, row.names = FALSE, na = "")
utils::write.csv(aggregation_plan$source_assignment, source_assignment_path, row.names = FALSE, na = "")

feature_audit <- data.frame(
  source_row_index = seq_along(historical$protein_ids),
  source_protein_id = historical$protein_ids,
  output_row_index = seq_along(historical$protein_ids),
  output_protein_id = rownames(aggregated$primary),
  Description = descriptions,
  description_source = ifelse(description_fallback, "protein_id_fallback", "imputed_data.xlsx:Genes"),
  protein_id_unchanged = historical$protein_ids == rownames(aggregated$primary),
  protein_order_unchanged = seq_along(historical$protein_ids) == seq_along(rownames(aggregated$primary)),
  included_primary = TRUE,
  included_sensitivity = historical$protein_ids %in% rownames(aggregated$sensitivity),
  stringsAsFactors = FALSE
)
utils::write.csv(feature_audit, feature_audit_path, row.names = FALSE, na = "")

contrast_manifest <- primary_contrast_manifest()
# Historical hemisphere-level evidence: former Datasets/mapped/ is now
# 99_historical/datasets_mapped/ after the 2026-08-26 restructure.
evidence_subdir <- c(
  paired_cno_vs_paired_veh = "datasets_mapped/effects_chemogenetic_inhibition/paired",
  paired_veh_vs_unpaired_veh = "datasets_mapped/learning_signature/memory_ensemble",
  unpaired_cno_vs_unpaired_veh = "datasets_mapped/effects_chemogenetic_inhibition/unpaired"
)
contrast_manifest$historical_evidence_file <- file.path(
  project_root,
  "clusterProfiler",
  "99_historical",
  unname(evidence_subdir[contrast_manifest$contrast_family]),
  paste0(contrast_manifest$historical_comparison_name, ".csv")
)
contrast_manifest$historical_evidence_exists <- file.exists(contrast_manifest$historical_evidence_file)
contrast_manifest$recommended_model <- "within_sample_class_animal_level_two_sample_ProTigy"
contrast_manifest <- contrast_manifest[c(
  "sample_class", "contrast_family", "case_code", "case_condition", "reference_code",
  "reference_condition", "case_phenotype", "reference_phenotype", "canonical_contrast",
  "historical_comparison_name", "historical_evidence_file", "historical_evidence_exists",
  "recommended_model"
)]
if (any(!contrast_manifest$historical_evidence_exists)) {
  stop("One or more organized historical contrast evidence files are missing.", call. = FALSE)
}
utils::write.csv(contrast_manifest, contrast_manifest_path, row.names = FALSE, na = "")

contrast_scope_audit <- data.frame(
  candidate_contrast = c(
    "paired_cno_vs_paired_veh",
    "paired_veh_vs_unpaired_veh",
    "unpaired_cno_vs_unpaired_veh",
    "paired_cno_vs_unpaired_cno",
    "cross_sample_class_same_condition"
  ),
  include_in_primary_contrast_manifest = c(TRUE, TRUE, TRUE, FALSE, FALSE),
  evidence = c(
    "organized mapped/effects_chemogenetic_inhibition/paired consumers",
    "organized mapped/learning_signature/memory_ensemble consumers",
    "organized mapped/effects_chemogenetic_inhibition/unpaired consumers",
    "present only in generic all-pairwise outputs; absent from organized current downstream consumers",
    "historical baseline profiling/all-pairwise outputs"
  ),
  disposition = c(
    rep("preserve historical within-sample-class orientation", 3),
    "do not add without a verified scientific consumer",
    "requires paired/animal-aware modeling; do not use ordinary independent two-sample ProTigy"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(contrast_scope_audit, contrast_scope_path, row.names = FALSE, na = "")

metadata_contract_audit <- data.frame(
  contract_item = c(
    "hemisphere_primary_source",
    "legacy_ReplicateGroup_crosscheck",
    "sample_id_LR_crosscheck",
    "source_sample_identity",
    "condition_semantics"
  ),
  status = c("verified", "verified", "unavailable_historical_schema", "verified", "verified"),
  evidence = c(
    "sample_annotation.xlsx label has one explicit _left/_right suffix",
    "sample_info.xlsx and sample_annotation.xlsx numeric 1/2 agree one-to-one with explicit labels",
    "2024 instrument IDs contain neither _L_ nor _R_; no inference made from sample_id",
    "sample_info id equals quantitative GCT columns and annotation Name basename for all 96 samples",
    "ExpGroup maps through R/analysis_labels.R to the four canonical conditions"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(metadata_contract_audit, metadata_contract_path, row.names = FALSE, na = "")

audit_xlsx <- file.path(output_dir, "animal_level_handoff_audits.xlsx")
writexl::write_xlsx(
  list(
    Aggregation_Audit = aggregation_plan$audit,
    Aggregation_Summary = aggregation_summary,
    Source_Assignment = aggregation_plan$source_assignment,
    Feature_Identity = feature_audit,
    Contrast_Manifest = contrast_manifest,
    Contrast_Scope = contrast_scope_audit,
    Metadata_Contract = metadata_contract_audit
  ),
  audit_xlsx
)

sha256_file <- function(path) digest::digest(file = path, algo = "sha256")
output_files <- c(
  primary_gct = primary_gct,
  primary_xlsx = primary_xlsx,
  aggregation_audit = aggregation_audit_path,
  aggregation_summary = aggregation_summary_path,
  source_sample_assignment = source_assignment_path,
  feature_identity_audit = feature_audit_path,
  protigy_contrast_manifest = contrast_manifest_path,
  contrast_scope_audit = contrast_scope_path,
  metadata_contract_audit = metadata_contract_path,
  audit_xlsx = audit_xlsx
)
if (has_incomplete_pairs) output_files <- c(output_files, sensitivity_gct = sensitivity_gct)
manifest_files <- rbind(
  data.frame(file_role = "input", artifact = names(inputs), path = unname(inputs), stringsAsFactors = FALSE),
  data.frame(file_role = "code", artifact = names(code_inputs), path = unname(code_inputs), stringsAsFactors = FALSE),
  data.frame(file_role = "output", artifact = names(output_files), path = unname(output_files), stringsAsFactors = FALSE)
)
manifest_files$sha256 <- vapply(manifest_files$path, sha256_file, character(1))
manifest_files$bytes <- file.info(manifest_files$path)$size
manifest_files$run_timestamp_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
manifest_files$aggregation_policy <- "equal_weight_mean_LR_on_existing_imputed_log2_values"
manifest_files$description_policy <- if (any(description_fallback)) {
  "Genes_with_protein_id_fallback_for_missing_descriptions"
} else {
  "Genes_preserved_one_to_one_no_fallback_needed"
}
manifest_files$sensitivity_status <- if (has_incomplete_pairs) {
  "written_because_incomplete_pairs_exist"
} else {
  "not_written_primary_equals_strict_complete_bilateral"
}
manifest_files$sample_id_hemisphere_crosscheck <- "unavailable_in_historical_2024_instrument_ids; explicit_annotation_label_used"
run_manifest_path <- file.path(output_dir, "run_manifest.csv")
utils::write.csv(manifest_files, run_manifest_path, row.names = FALSE, na = "")

cat("Animal-level ProTigy handoff completed.\n")
cat("Quantitative source:", quantitative_gct, "\n")
cat("Metadata sources:", sample_info_xlsx, "and", sample_annotation_xlsx, "\n")
cat("Proteins:", nrow(historical$matrix), "before;", nrow(aggregated$primary), "after\n")
cat("Samples:", ncol(historical$matrix), "hemisphere-level;", ncol(aggregated$primary), "animal-level\n")
cat("Complete pairs:", sum(aggregation_plan$audit$hemisphere_status == "bilateral_complete"), "\n")
cat("One-sided pairs:", sum(aggregation_plan$audit$hemisphere_status %in% c("left_only_observed", "right_only_observed")), "\n")
cat("Missing pairs:", sum(aggregation_plan$audit$hemisphere_status == "missing_both"), "\n")
cat("Primary GCT:", primary_gct, "\n")
cat("Sensitivity GCT:", if (has_incomplete_pairs) sensitivity_gct else "not written (identical to primary)", "\n")
cat("Run manifest:", run_manifest_path, "\n")
