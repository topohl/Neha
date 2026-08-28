#!/usr/bin/env Rscript

# Fail-closed contracts for 03_qc_exploration/02_rank_abundance_by_sample_class.r.
#
# Sections 1-3 are self-contained and always run. They execute the real stage
# against a valid 48-observation animal-level GCT and a syntactically valid but
# biologically wrong 96-observation hemisphere-style GCT. Section 4 performs
# full keyed equivalence against the locked canonical inputs and finalized
# reviewer-revision reference when the shared drive is reachable; that external
# section is reported explicitly as SKIP when its inputs are unavailable.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1) sub("^--file=", "", file_arg) else
  file.path("tests", "test_rank_abundance_animal_level.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "protigy_input_utils.R"))
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required for rank-abundance provenance contracts.", call. = FALSE)
}

stage <- normalizePath(
  file.path(repo_root, "03_qc_exploration", "02_rank_abundance_by_sample_class.r"),
  winslash = "/",
  mustWork = TRUE
)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

failures <- character()
checks <- 0L
external_status <- "NOT_RUN"
ok <- function(label, condition) {
  checks <<- checks + 1L
  if (!isTRUE(condition)) {
    failures <<- c(failures, label)
    cat("  [FAIL] ", label, "\n")
  } else {
    cat("  [ok]   ", label, "\n")
  }
}

strip_r_comments <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parsed <- parse(file = path, keep.source = TRUE)
  parse_data <- getParseData(parsed)
  comments <- parse_data[parse_data$token == "COMMENT", , drop = FALSE]
  if (nrow(comments)) {
    comments <- comments[order(comments$line1, comments$col1, decreasing = TRUE), , drop = FALSE]
    for (i in seq_len(nrow(comments))) {
      line_number <- comments$line1[[i]]
      lines[[line_number]] <- if (comments$col1[[i]] <= 1L) "" else
        substr(lines[[line_number]], 1L, comments$col1[[i]] - 1L)
    }
  }
  paste(lines, collapse = "\n")
}

sha256_file <- function(path) {
  tolower(digest::digest(file = path, algo = "sha256"))
}

restore_environment <- function(previous) {
  for (name in names(previous)) {
    if (is.na(previous[[name]])) Sys.unsetenv(name) else
      do.call(Sys.setenv, stats::setNames(list(previous[[name]]), name))
  }
}

run_stage <- function(input_gct, id_map, output_dir, log_name, scratch_root) {
  env_names <- c(
    "NEHA_RANK_ABUNDANCE_INPUT_GCT",
    "NEHA_RANK_ABUNDANCE_ID_MAP",
    "NEHA_RANK_ABUNDANCE_OUTPUT_DIR"
  )
  previous <- stats::setNames(
    vapply(env_names, function(name) Sys.getenv(name, unset = NA_character_), character(1)),
    env_names
  )
  old_wd <- getwd()
  on.exit({
    setwd(old_wd)
    restore_environment(previous)
  }, add = TRUE)
  Sys.setenv(
    NEHA_RANK_ABUNDANCE_INPUT_GCT = input_gct,
    NEHA_RANK_ABUNDANCE_ID_MAP = id_map,
    NEHA_RANK_ABUNDANCE_OUTPUT_DIR = output_dir
  )
  setwd(scratch_root)
  log_path <- file.path(scratch_root, log_name)
  status <- system2(rscript, c("--vanilla", shQuote(stage)), stdout = log_path, stderr = log_path)
  list(
    status = status,
    log_path = log_path,
    log = if (file.exists(log_path)) paste(readLines(log_path, warn = FALSE), collapse = "\n") else ""
  )
}

read_rank_csv <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

compare_complete_keyed_output <- function(observed_path, reference_path, label) {
  observed <- read_rank_csv(observed_path)
  reference <- read_rank_csv(reference_path)
  required <- c("Genes", "Condition", "MeanLog2", "LinearValue", "Rank", "n_animals", "animals", "MarkerType")
  ok(paste0(label, ": required columns"),
     identical(names(observed), required) && identical(names(reference), required))
  if (!all(required %in% names(observed)) || !all(required %in% names(reference))) return(invisible(FALSE))

  observed_key <- paste(observed$Condition, observed$Genes, sep = "\r")
  reference_key <- paste(reference$Condition, reference$Genes, sep = "\r")
  ok(paste0(label, ": unique Condition x Genes key"),
     !anyDuplicated(observed_key) && !anyDuplicated(reference_key))
  ok(paste0(label, ": complete keyed row set"), setequal(observed_key, reference_key))
  idx <- match(reference_key, observed_key)
  if (anyNA(idx)) return(invisible(FALSE))

  ok(paste0(label, ": gene key exact"), identical(observed$Genes[idx], reference$Genes))
  ok(paste0(label, ": MeanLog2 exact"),
     isTRUE(all.equal(observed$MeanLog2[idx], reference$MeanLog2, tolerance = 0, check.attributes = FALSE)))
  ok(paste0(label, ": LinearValue exact"),
     isTRUE(all.equal(observed$LinearValue[idx], reference$LinearValue, tolerance = 0, check.attributes = FALSE)))
  ok(paste0(label, ": rank exact"), identical(observed$Rank[idx], reference$Rank))
  ok(paste0(label, ": MarkerType exact"), identical(observed$MarkerType[idx], reference$MarkerType))
  ok(paste0(label, ": n_animals exact"), identical(observed$n_animals[idx], reference$n_animals))
  ok(paste0(label, ": animal list exact"), identical(observed$animals[idx], reference$animals))
  invisible(TRUE)
}

cat("\n=== 1. executable source-code contracts (always run) ===\n")

code <- strip_r_comments(stage)
ok("animal-level GCT default appears in executable code",
   grepl("neha_protigy_input_animal_level_primary\\.gct", code))
ok("input override appears in executable code",
   grepl("NEHA_RANK_ABUNDANCE_INPUT_GCT", code, fixed = TRUE))
ok("shared ProTigy GCT validator is called in executable code",
   grepl("validate_protigy_gct_v13", code, fixed = TRUE) && grepl("protigy_input_utils\\.R", code))
ok("does NOT read per-class imputed workbooks", !grepl("pgmatrix_imputed", code, fixed = TRUE))
ok("does NOT point at the historical imputed-workbook directory",
   !grepl('"gct"\\s*,\\s*"imputed"', code) && !grepl("gct/imputed", code, fixed = TRUE))
ok("does NOT use the old directory override",
   !grepl("NEHA_RANK_ABUNDANCE_INPUT_DIR", code, fixed = TRUE))
ok("does NOT reference the historical hemisphere-level GCT",
   !grepl("pg.matrix_filtered_pcaAdjusted_unnormalized.gct", code, fixed = TRUE))
ok("groups using normalized sample_class and condition metadata",
   grepl("normalize_sample_class", code, fixed = TRUE) && grepl("normalize_condition", code, fixed = TRUE))
ok("does NOT renormalize the processed GCT matrix",
   !grepl("scale\\(|normalizeQuantiles|normalizeBetweenArrays|median_center|quantile_norm", code))
ok("does NOT add imputation", !grepl("missForest|\\bknn\\s*\\(", code, ignore.case = TRUE))
ok("linearization remains 2^MeanLog2", grepl("2\\s*\\^\\s*d\\$MeanLog2", code))

for (category in c("General Neuron", "Neuropil/Structure", "Stress Response", "Activation")) {
  ok(paste0("marker category retained: ", category), grepl(category, code, fixed = TRUE))
}
ok("Fig3E group contract retained",
   grepl('"mcherry_paired-veh", "mcherry_unpaired-veh"', code, fixed = TRUE))
ok("SuppD group contract retained",
   grepl('"neuron_unpaired-veh", "neuropil_unpaired-veh"', code, fixed = TRUE))

cat("\n=== 2. self-contained GCT fixtures (always run) ===\n")

scratch_root <- file.path(tempdir(), paste0("rank_abundance_contract_", Sys.getpid()))
if (dir.exists(scratch_root)) unlink(scratch_root, recursive = TRUE)
dir.create(scratch_root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(scratch_root, recursive = TRUE), add = TRUE)

fixture_conditions <- unname(condition_code_map)
fixture_classes <- sample_classes
fixture_metadata <- do.call(rbind, lapply(seq_along(fixture_conditions), function(condition_index) {
  condition_name <- fixture_conditions[[condition_index]]
  animal_ids <- sprintf("A%02d", seq.int((condition_index - 1L) * 3L + 1L, condition_index * 3L))
  do.call(rbind, lapply(fixture_classes, function(class_name) {
    data.frame(
      output_column_name = paste(animal_ids, class_name, sep = "_"),
      AnimalID = animal_ids,
      condition_code = names(condition_code_map)[[condition_index]],
      condition = condition_name,
      sample_class = class_name,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
}))
rownames(fixture_metadata) <- NULL

fixture_protein_ids <- sprintf("FIXTURE_%02d_MOUSE", seq_len(8L))
fixture_genes <- c("Rps9", "Vcan", "Aktip", "Mrpl4", "GeneA", "GeneB", "GeneC", "GeneD")
fixture_matrix <- outer(seq_along(fixture_protein_ids), seq_len(nrow(fixture_metadata)), function(i, j) {
  0.25 * i + 0.01 * j
})
class_effect <- match(fixture_metadata$sample_class, fixture_classes) * 0.03
condition_effect <- match(fixture_metadata$condition, fixture_conditions) * 0.02
fixture_matrix <- sweep(fixture_matrix, 2L, class_effect + condition_effect, "+")
rownames(fixture_matrix) <- fixture_protein_ids
colnames(fixture_matrix) <- fixture_metadata$output_column_name

fixture_map <- file.path(scratch_root, "fixture_id_map.csv")
utils::write.csv(
  data.frame(
    original_protein_id = fixture_protein_ids,
    mapped_gene_symbol = fixture_genes,
    stringsAsFactors = FALSE
  ),
  fixture_map,
  row.names = FALSE
)

valid_gct <- file.path(scratch_root, "valid_animal_level.gct")
write_protigy_gct_v13(
  fixture_matrix,
  fixture_metadata,
  paste("fixture", fixture_genes),
  valid_gct
)

hemisphere_metadata <- fixture_metadata[rep(seq_len(nrow(fixture_metadata)), each = 2L), , drop = FALSE]
hemisphere_side <- rep(c("L", "R"), times = nrow(fixture_metadata))
hemisphere_metadata$output_column_name <- paste(hemisphere_metadata$output_column_name, hemisphere_side, sep = "_")
hemisphere_matrix <- fixture_matrix[, rep(seq_len(ncol(fixture_matrix)), each = 2L), drop = FALSE]
hemisphere_matrix <- sweep(hemisphere_matrix, 2L, rep(c(-0.001, 0.001), times = ncol(fixture_matrix)), "+")
colnames(hemisphere_matrix) <- hemisphere_metadata$output_column_name
wrong_gct <- file.path(scratch_root, "wrong_hemisphere_level_96.gct")
write_protigy_gct_v13(
  hemisphere_matrix,
  hemisphere_metadata,
  paste("fixture", fixture_genes),
  wrong_gct
)

parsed_valid <- validate_protigy_gct_v13(valid_gct, expected_matrix = fixture_matrix)
parsed_wrong <- validate_protigy_gct_v13(wrong_gct, expected_matrix = hemisphere_matrix)
ok("valid fixture is a syntactically valid 48-column ProTigy GCT",
   nrow(parsed_valid$matrix) == 8L && ncol(parsed_valid$matrix) == 48L)
ok("wrong fixture is a syntactically valid 96-column ProTigy GCT",
   nrow(parsed_wrong$matrix) == 8L && ncol(parsed_wrong$matrix) == 96L)

cat("\n=== 3. real-stage fixture execution (always run) ===\n")

wrong_output <- file.path(scratch_root, "wrong_output_must_not_exist")
wrong_run <- run_stage(wrong_gct, fixture_map, wrong_output, "wrong_stage.log", scratch_root)
ok("NEGATIVE INPUT: 96-observation hemisphere-style GCT fails closed", wrong_run$status != 0L)
ok("NEGATIVE INPUT: failure names the animal-level experimental-unit contract",
   grepl("exactly 48 unique animal-level sample columns|hemisphere-level technical", wrong_run$log))
ok("NEGATIVE INPUT: output directory is not created", !dir.exists(wrong_output))

valid_output <- file.path(scratch_root, "valid_output")
valid_run <- run_stage(valid_gct, fixture_map, valid_output, "valid_stage.log", scratch_root)
ok("POSITIVE INPUT: 48-observation animal-level fixture passes", identical(valid_run$status, 0L))
if (valid_run$status != 0L) cat(valid_run$log, "\n")

fixture_provenance_path <- file.path(valid_output, "rank_abundance_run_provenance.csv")
ok("POSITIVE INPUT: provenance manifest is written", file.exists(fixture_provenance_path))
if (file.exists(fixture_provenance_path)) {
  fixture_provenance <- utils::read.csv(fixture_provenance_path, stringsAsFactors = FALSE, check.names = FALSE)
  ok("POSITIVE INPUT: override path was actually consumed",
     identical(unique(fixture_provenance$input_gct_path), normalizePath(valid_gct, winslash = "/", mustWork = TRUE)))
  ok("POSITIVE INPUT: fixture GCT hash is recorded",
     identical(unique(fixture_provenance$input_gct_sha256), sha256_file(valid_gct)))
  ok("POSITIVE INPUT: exactly 48 sample columns are recorded",
     identical(unique(fixture_provenance$gct_sample_columns), 48L))
  ok("POSITIVE INPUT: exactly 12 distinct animals are recorded",
     identical(unique(fixture_provenance$distinct_animal_ids), 12L))
  ok("POSITIVE INPUT: four canonical sample classes are recorded",
     identical(unique(fixture_provenance$sample_classes), paste(sample_classes, collapse = ";")))
  ok("POSITIVE INPUT: four canonical conditions are recorded",
     identical(unique(fixture_provenance$conditions), paste(condition_levels, collapse = ";")))
  ok("POSITIVE INPUT: all 16 strata are present", nrow(fixture_provenance) == 16L)
  ok("POSITIVE INPUT: three observations per stratum",
     all(fixture_provenance$n_observations == 3L))
  ok("POSITIVE INPUT: three distinct animals per stratum",
     all(fixture_provenance$n_distinct_animals == 3L))
  ok("POSITIVE INPUT: mapping override hash is recorded",
     identical(unique(fixture_provenance$mapping_input_sha256), sha256_file(fixture_map)))
}

fixture_ranks <- file.path(valid_output, "processed_protein_ranks_animal_level.csv")
ok("POSITIVE INPUT: complete rank output is written", file.exists(fixture_ranks))
if (file.exists(fixture_ranks)) {
  fixture_rank_data <- read_rank_csv(fixture_ranks)
  ok("POSITIVE INPUT: 16 fixture groups are produced", length(unique(fixture_rank_data$Condition)) == 16L)
  ok("POSITIVE INPUT: every fixture group records three animals", all(fixture_rank_data$n_animals == 3L))
  ok("POSITIVE INPUT: fixture has no duplicated gene key",
     !anyDuplicated(paste(fixture_rank_data$Condition, fixture_rank_data$Genes, sep = "\r")))
}

cat("\n=== 4. locked canonical and full reference equivalence (external) ===\n")

shared_root <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
canonical_gct <- file.path(shared_root, "02_data", "animal_level", "input_gct",
                           "neha_protigy_input_animal_level_primary.gct")
canonical_map <- file.path(shared_root, "02_data", "animal_level", "mapped", "forward",
                           "mcherry_paired_veh_vs_mcherry_unpaired_veh.csv")
revision_root <- file.path(shared_root, "03_output", "reviewer_revision_animal_level_20260827")
reference_all <- file.path(revision_root, "full_regenerated", "qc",
                           "Processed_Protein_Ranks_animal_level.csv")
reference_fig3e <- file.path(revision_root, "figure_panels", "Figure3",
                             "Fig3E_rank_abundance_animal_level_source_data.csv")
reference_suppd <- file.path(revision_root, "figure_panels", "Supplementary_proteomics",
                             "SuppD_rank_abundance_animal_level_source_data.csv")
external_inputs <- c(canonical_gct, canonical_map, reference_all, reference_fig3e, reference_suppd)

if (!all(file.exists(external_inputs))) {
  external_status <- "SKIP: canonical GCT, mapping, or finalized reference unavailable"
  cat("  [SKIP] ", external_status, "\n")
} else {
  external_status <- "PASS"
  locked_gct_sha256 <- "f12cf99e1bfb7c17bbf56bffb6783e924698bce5d5533a8e312bc4bbb733bbb3"
  locked_map_sha256 <- "0e8a2100a34b4c2afe8bcf83949cc5ca7bb0597999edaaf83e970930bcf55a66"
  ok("CANONICAL: locked animal-level GCT SHA-256",
     identical(sha256_file(canonical_gct), locked_gct_sha256))
  ok("CANONICAL: locked mapping-input SHA-256",
     identical(sha256_file(canonical_map), locked_map_sha256))
  canonical_parsed <- validate_protigy_gct_v13(canonical_gct)
  ok("CANONICAL: 5349 protein rows", nrow(canonical_parsed$matrix) == 5349L)
  ok("CANONICAL: 48 sample columns", ncol(canonical_parsed$matrix) == 48L)

  canonical_output <- file.path(scratch_root, "canonical_output")
  canonical_run <- run_stage(
    canonical_gct,
    canonical_map,
    canonical_output,
    "canonical_stage.log",
    scratch_root
  )
  ok("CANONICAL: real stage completes in the isolated output root", identical(canonical_run$status, 0L))
  if (canonical_run$status != 0L) cat(canonical_run$log, "\n")

  canonical_provenance_path <- file.path(canonical_output, "rank_abundance_run_provenance.csv")
  ok("CANONICAL: provenance manifest is written", file.exists(canonical_provenance_path))
  if (file.exists(canonical_provenance_path)) {
    canonical_provenance <- utils::read.csv(canonical_provenance_path, stringsAsFactors = FALSE, check.names = FALSE)
    ok("CANONICAL: provenance carries locked GCT hash",
       identical(unique(canonical_provenance$input_gct_sha256), locked_gct_sha256))
    ok("CANONICAL: provenance carries locked mapping hash",
       identical(unique(canonical_provenance$mapping_input_sha256), locked_map_sha256))
    ok("CANONICAL: provenance records 5349 x 48 input",
       identical(unique(canonical_provenance$gct_protein_rows), 5349L) &&
         identical(unique(canonical_provenance$gct_sample_columns), 48L))
    ok("CANONICAL: provenance records 5310 final features",
       identical(unique(canonical_provenance$final_feature_count), 5310L))
  }

  compare_complete_keyed_output(
    file.path(canonical_output, "processed_protein_ranks_animal_level.csv"),
    reference_all,
    "OUTPUT all 16 groups"
  )
  compare_complete_keyed_output(
    file.path(canonical_output, "Fig3E_rank_abundance_animal_level_source_data.csv"),
    reference_fig3e,
    "OUTPUT Fig3E"
  )
  compare_complete_keyed_output(
    file.path(canonical_output, "SuppD_rank_abundance_animal_level_source_data.csv"),
    reference_suppd,
    "OUTPUT SuppD"
  )
}

cat("\n=== RESULT ===\n")
cat("External canonical/reference section: ", external_status, "\n", sep = "")
if (length(failures)) {
  stop(
    paste(
      c(
        sprintf("Rank-abundance animal-level tests failed (%d of %d checks):", length(failures), checks),
        paste0("  - ", failures)
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}
cat(sprintf("Rank-abundance self-contained and available external contracts passed (%d checks).\n", checks))
