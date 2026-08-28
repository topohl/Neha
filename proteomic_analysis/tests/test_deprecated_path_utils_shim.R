#!/usr/bin/env Rscript

# Contracts for the deprecated R/neha_path_utils.R compatibility shim.
#
# The shim exists only so that the byte-identical, hash-locked figure-revision snapshots in
# 06_manuscript_figure_revision/ (and their runnable twins on the shared drive) keep working
# after the 2026-08-28 rename to R/project_path_utils.R. These checks pin that contract so the
# shim cannot silently rot, and so its eventual removal is a deliberate, visible decision.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("tests", "test_deprecated_path_utils_shim.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

failures <- 0L
expect <- function(condition, message) {
  if (!isTRUE(condition)) {
    failures <<- failures + 1L
    cat("  [FAIL]", message, "\n")
  } else {
    cat("  [ok]  ", message, "\n")
  }
}

shim <- file.path(repo_root, "R", "neha_path_utils.R")
impl <- file.path(repo_root, "R", "project_path_utils.R")

cat("=== files ===\n")
expect(file.exists(impl), "R/project_path_utils.R is the active implementation")
expect(file.exists(shim), "R/neha_path_utils.R deprecated shim is present")

cat("\n=== the shim carries no analysis logic ===\n")
shim_code <- paste(readLines(shim, warn = FALSE), collapse = "\n")
expect(!grepl("normalizePath\\(path", shim_code) && !grepl("prcomp", shim_code, fixed = TRUE),
       "shim defines no path or PCA logic of its own -- it only delegates")

cat("\n=== sourcing the shim supplies the path primitives ===\n")
# Mimic the frozen snapshot's source order: shim first, then the PCA utils.
source(file.path(repo_root, "R", "analysis_labels.R"))
source(shim)
expect(exists("project_normalize_path", mode = "function"),
       "shim loads project_normalize_path() so the PCA utils guard passes")
expect(exists("PROTEOMICS_DATA_ROOT"),
       "shim loads PROTEOMICS_DATA_ROOT")
expect(identical(project_path_is_within("C:/a/b", "C:/a"), TRUE),
       "delegated project_path_is_within() behaves correctly")

source(file.path(repo_root, "R", "protigy_input_utils.R"))
source(file.path(repo_root, "R", "pca_animal_level_utils.R"))

cat("\n=== deprecated aliases work and warn ===\n")
expect(exists("validate_neha_pca_animal_input", mode = "function"),
       "validate_neha_pca_animal_input alias is defined")
expect(exists("prepare_neha_animal_pca", mode = "function"),
       "prepare_neha_animal_pca alias is defined")

# Build a real animal-level GCT fixture and parse it, exactly as test_pca_animal_level.R does:
# 4 sample classes x 4 conditions x 3 animals = 48 units, matching the production contract.
animal_ids_by_condition <- list(
  paired_cno = c("C11", "C12", "C34f"),
  paired_veh = c("C14", "C33f", "C510"),
  unpaired_cno = c("C25", "C45f", "C46"),
  unpaired_veh = c("C26", "C27", "C47")
)
fixture_metadata <- do.call(rbind, lapply(sample_classes, function(sample_class) {
  do.call(rbind, lapply(names(animal_ids_by_condition), function(condition) {
    ids <- animal_ids_by_condition[[condition]]
    data.frame(
      output_column_name = paste(ids, sample_class, sep = "_"),
      AnimalID = ids,
      condition_code = names(condition_code_map)[match(condition, condition_code_map)],
      condition = condition,
      sample_class = sample_class,
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(fixture_metadata) <- NULL
fixture_matrix <- outer(
  seq_len(8L), seq_len(nrow(fixture_metadata)),
  function(protein, sample) protein * 0.3 + sample * 0.01 + ((protein * sample) %% 7) * 0.02
)
fixture_matrix[8, ] <- 5
rownames(fixture_matrix) <- paste0("P", seq_len(nrow(fixture_matrix)), "_MOUSE")
colnames(fixture_matrix) <- fixture_metadata$output_column_name
fixture_path <- tempfile("shim_animal_pca_", fileext = ".gct")
write_protigy_gct_v13(fixture_matrix, fixture_metadata,
                      paste0("Gene", seq_len(nrow(fixture_matrix))), fixture_path)
parsed <- validate_protigy_gct_v13(fixture_path, expected_matrix = fixture_matrix)

new_validated <- validate_pca_animal_input(parsed, expected_n = 3L)
old_warned <- FALSE
old_validated <- withCallingHandlers(
  validate_neha_pca_animal_input(parsed, expected_n = 3L),
  warning = function(w) {
    if (grepl("deprecated", conditionMessage(w))) old_warned <<- TRUE
    invokeRestart("muffleWarning")
  }
)
expect(old_warned, "validate_neha_pca_animal_input() emits a deprecation warning")
expect(identical(old_validated, new_validated),
       "deprecated validator returns exactly what validate_pca_animal_input() returns")

new_prepared <- prepare_animal_pca(new_validated$expression_matrix, center = TRUE, scale. = TRUE)
prep_warned <- FALSE
old_prepared <- withCallingHandlers(
  prepare_neha_animal_pca(new_validated$expression_matrix, center = TRUE, scale. = TRUE),
  warning = function(w) {
    if (grepl("deprecated", conditionMessage(w))) prep_warned <<- TRUE
    invokeRestart("muffleWarning")
  }
)
expect(prep_warned, "prepare_neha_animal_pca() emits a deprecation warning")
expect(identical(old_prepared, new_prepared),
       "deprecated preparer returns exactly what prepare_animal_pca() returns")

cat("\n=== the active pipeline does not depend on the shim ===\n")
active <- list.files(repo_root, recursive = TRUE, full.names = TRUE,
                     pattern = "\\.(r|R|ps1)$")
active <- active[!grepl("/06_manuscript_figure_revision/|/tests/test_deprecated_path_utils_shim\\.R$|/R/neha_path_utils\\.R$",
                        normalizePath(active, winslash = "/", mustWork = FALSE))]
# Look for an actual source() of the shim, not a mere mention of its name: the stale-label
# audit legitimately names the file in its exception list, and the docs describe it.
offenders <- Filter(function(f) {
  any(grepl("source[[:space:]]*\\([^)]*neha_path_utils", readLines(f, warn = FALSE)))
}, active)
expect(length(offenders) == 0L,
       paste0("no active file sources the deprecated shim",
              if (length(offenders)) paste0(" (found: ", paste(basename(offenders), collapse = ", "), ")") else ""))

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Deprecated-shim contracts failed: %d", failures), call. = FALSE)
}
cat("All deprecated-shim contracts hold.\n")
