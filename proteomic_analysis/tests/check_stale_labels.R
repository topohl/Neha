args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1) sub("^--file=", "", file_arg) else file.path("tests", "check_stale_labels.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(repo_root)) repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

all_files <- list.files(repo_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
# Only active pipeline code is checked. /legacy/ and 06_manuscript_figure_revision/ are
# skipped; the latter holds copies of figure scripts that were already run (see its README)
# and is not edited to pass this check.
archived <- "/legacy/|/06_manuscript_figure_revision/"
active_files <- all_files[
  grepl("\\.(r|R|md|Rmd)$", all_files) &
    !grepl(archived, normalizePath(all_files, winslash = "/", mustWork = FALSE)) &
    !grepl("/tests/check_stale_labels\\.R$", normalizePath(all_files, winslash = "/", mustWork = FALSE))
]

read_file <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")
contents <- setNames(lapply(active_files, read_file), normalizePath(active_files, winslash = "/", mustWork = FALSE))

failures <- character()

add_hits <- function(label, pattern, ignore_case = TRUE) {
  hit <- vapply(contents, function(x) grepl(pattern, x, perl = TRUE, ignore.case = ignore_case), logical(1))
  if (any(hit)) {
    failures <<- c(failures, paste0(label, ": ", paste(names(contents)[hit], collapse = ", ")))
  }
}

add_hits("old con/res/sus condition mapping", '"[123]"\\s*=\\s*"(con|res|sus)"|phenotypes\\s*=\\s*c\\([^)]*"(con|res|sus)"|\\b(con|res|sus)_vs_(con|res|sus)\\b')
add_hits("group-code regex missing code 4", "\\[123\\]")
add_hits("old active sample-class parsing labels", "neuron_soma|neuron_neuropil|sample_class[^\\n]*(microglia|celltype_layer)|celltype_layer[^\\n]*sample_class")

# Output names. Prose may mention manuscripts, publications or the nature of a measurement,
# but file and helper names must not carry a journal name, a drafting status or a retired
# theme name. The manuscript/publication rule needs a separator, so prose such as
# "manuscript figure" still passes.
output_label_rules <- list(
  journal_named_output = "(?i)(?<![A-Za-z])nature(?=(?:[_-](?:fig(?:ure)?[0-9]*|panel|plot|supp)|figure|\\.(?:svg|pdf|png|tiff?)))",
  drafting_status_output = "(?i)(?<![A-Za-z])(publication|manuscript)[-_]+(ready|quality)(?:[-_]+figure)?",
  retired_theme_helper = "(?i)(?<![A-Za-z0-9])theme[-_]?nature(?:_qc)?"
)

has_stale_output_label <- function(text) {
  any(vapply(output_label_rules, function(pattern) grepl(pattern, text, perl = TRUE), logical(1)))
}

# Known cases each rule must reject or allow, so a regex change that alters them fails here.
control_cases <- data.frame(
  text = c(
    "Nature_Fig3.svg",
    "NatureFigure.svg",
    "nature_Fig3.svg",
    "publication_ready_theme",
    "publication-quality-figure",
    "manuscript_ready.svg",
    "theme_nature",
    "This script regenerates manuscript panels.",
    "Publication details are documented in README.",
    "Nature of the measurement is descriptive."
  ),
  expected_rejected = c(rep(TRUE, 7L), rep(FALSE, 3L)),
  stringsAsFactors = FALSE
)
control_cases$observed_rejected <- vapply(control_cases$text, has_stale_output_label, logical(1))
misclassified <- control_cases$observed_rejected != control_cases$expected_rejected
if (any(misclassified)) {
  failures <- c(
    failures,
    paste0(
      "stale-label control misclassified: ",
      paste(control_cases$text[misclassified], collapse = " | ")
    )
  )
}

for (rule_name in names(output_label_rules)) {
  add_hits(paste0("stale output identifier (", rule_name, ")"), output_label_rules[[rule_name]])
}

# Collaborator name. Since 2026-08-28 the project is "Associative Memory Proteomics", and the
# collaborator name may appear in active code, file names and prose only as one of the
# exceptions below: the shared-drive path (Collabs/Neha), hashed artefact file names that
# cannot be renamed without breaking recorded hashes, the dataset_project value written into
# the aggregation audit and every canonical_analysis_unit key, the split version token in a
# validated manifest, and the public repository URL with the folder a clone creates.
# R/neha_path_utils.R is the deprecated shim and is exempt (see
# tests/test_deprecated_path_utils_shim.R).
legacy_collaborator_exceptions <- c(
  "Collabs/Neha",
  "Collabs\\Neha",
  "neha_protigy_input_animal_level_strict_complete_bilateral",
  "neha_protigy_input_animal_level_primary",
  "stat_results_for_ssGSEA_neha_proteome",
  "neha_animal_level_protigy_stat_split_v1",
  'dataset_project = "Neha"',
  "Neha__",
  "<Neha>",
  # Old variable name; tests/test_rank_abundance_animal_level.R checks it is no longer used.
  "NEHA_RANK_ABUNDANCE_INPUT_DIR",
  # The deprecated shim's file name, which the docs mention. Whether any active file sources
  # it is checked in tests/test_deprecated_path_utils_shim.R.
  "neha_path_utils.R",
  # Shared-drive files the release cites by path: the ProTigy parameter file for the
  # animal-level run and the 2025-11-07 ProTigy session folder, whose params.txt is the only
  # record of the ProTigy version. The data files keep their names.
  "neha_proteome_parameters.yaml",
  "20251107_pg.matrix_Neha",
  # The folder a clone of the public repository creates (README).
  "cd Neha/proteomic_analysis"
)

# The public repository URL. "Neha" must be the whole repository name, so a longer name such
# as Neha-tools still fails.
repository_url_pattern <- "https://github\\.com/topohl/Neha(?![A-Za-z0-9_-])"

strip_legacy_exceptions <- function(text) {
  text <- gsub(repository_url_pattern, "", text, perl = TRUE)
  for (literal in legacy_collaborator_exceptions) text <- gsub(literal, "", text, fixed = TRUE)
  text
}

# Known cases for the name rule, as for the output-name rules above.
branding_controls <- data.frame(
  text = c(
    "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler",
    "neha_protigy_input_animal_level_primary.gct",
    'dataset_project = "Neha"',
    "# Neha proteomics workflow",
    "NEHA_PCA_OUTPUT_ROOT",
    "validate_neha_pca_animal_input",
    "source(file.path(repo_root, 'R', 'project_path_utils.R'))",
    "neha_proteome_parameters.yaml",
    "20251107_pg.matrix_Neha",
    "neha_publication_release_utils.R",
    "git clone https://github.com/topohl/Neha.git",
    "cd Neha/proteomic_analysis",
    "`https://github.com/topohl/Neha`",
    "https://github.com/topohl/Neha-tools",
    "cd Neha"
  ),
  expected_rejected = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE,
                        FALSE, FALSE, TRUE,
                        FALSE, FALSE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
branding_controls$observed_rejected <- vapply(
  branding_controls$text,
  function(x) grepl("neha", strip_legacy_exceptions(x), ignore.case = TRUE),
  logical(1)
)
branding_misclassified <- branding_controls$observed_rejected != branding_controls$expected_rejected
if (any(branding_misclassified)) {
  failures <- c(
    failures,
    paste0("branding control misclassified: ",
           paste(branding_controls$text[branding_misclassified], collapse = " | "))
  )
}

branding_exempt <- "/R/neha_path_utils\\.R$|/tests/test_deprecated_path_utils_shim\\.R$"
branding_hit <- vapply(names(contents), function(path) {
  if (grepl(branding_exempt, path)) return(FALSE)
  grepl("neha", strip_legacy_exceptions(contents[[path]]), ignore.case = TRUE)
}, logical(1))
if (any(branding_hit)) {
  failures <- c(
    failures,
    paste0("collaborator-name branding in active code: ",
           paste(names(contents)[branding_hit], collapse = ", "))
  )
}

# Active filenames must not carry the collaborator name either (same exemption for the shim).
branding_filenames <- active_files[
  grepl("neha", basename(active_files), ignore.case = TRUE) &
    !grepl("^neha_path_utils\\.R$", basename(active_files))
]
if (length(branding_filenames) > 0) {
  failures <- c(
    failures,
    paste0("collaborator-name branding in active filenames: ",
           paste(basename(branding_filenames), collapse = ", "))
  )
}

if (length(failures) > 0) {
  stop(paste(c("Stale-label audit failed:", failures), collapse = "\n"), call. = FALSE)
}

message(
  "Stale-label controls passed: ", sum(control_cases$expected_rejected), " must-reject, ",
  sum(!control_cases$expected_rejected), " must-allow. ",
  "Branding controls passed: ", sum(branding_controls$expected_rejected), " must-reject, ",
  sum(!branding_controls$expected_rejected), " must-allow. ",
  "Active-code audit passed over ", length(contents), " files."
)
