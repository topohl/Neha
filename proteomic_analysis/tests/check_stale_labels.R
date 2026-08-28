args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1) sub("^--file=", "", file_arg) else file.path("tests", "check_stale_labels.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(repo_root)) repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

all_files <- list.files(repo_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
# Archived, non-runnable snapshots are excluded the same way /legacy/ is: this audit
# checks ACTIVE pipeline code. 06_manuscript_figure_revision/ is a frozen provenance
# snapshot of already-executed code (see its README), not an active stage, and must not
# be edited to satisfy an active-pipeline stale-label lint.
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

# Output-identifier policy. Ordinary prose may discuss manuscripts, publications, or the
# nature of a measurement. Active filenames/helper identifiers must not encode a target
# journal, drafting status, or retired theme name. Separators are intentionally required
# for manuscript/publication identifiers so prose such as "manuscript figure" remains valid.
output_label_rules <- list(
  journal_named_output = "(?i)(?<![A-Za-z])nature(?=(?:[_-](?:fig(?:ure)?[0-9]*|panel|plot|supp)|figure|\\.(?:svg|pdf|png|tiff?)))",
  drafting_status_output = "(?i)(?<![A-Za-z])(publication|manuscript)[-_]+(ready|quality)(?:[-_]+figure)?",
  retired_theme_helper = "(?i)(?<![A-Za-z0-9])theme[-_]?nature(?:_qc)?"
)

has_stale_output_label <- function(text) {
  any(vapply(output_label_rules, function(pattern) grepl(pattern, text, perl = TRUE), logical(1)))
}

# Committed table-driven controls make changes to the output-label regex policy observable.
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

# Collaborator-name branding policy (added 2026-08-28). The project is "Associative Memory
# Proteomics"; the collaborator name must not return to active branding, identifiers, filenames
# or prose. A small, closed set of occurrences is legitimate and is enumerated here:
#
#   * the physical shared-drive storage path (Collabs/Neha) -- a real legacy location
#   * validated, SHA-locked artefact filenames, which cannot be renamed without invalidating
#     the recorded hashes and manifests
#   * dataset_project = "Neha", written verbatim into the validated aggregation audit as the
#     "dataset/project" column and into every canonical_analysis_unit key
#   * the frozen split contract_version token recorded in a validated manifest
#   * the frozen dated audit dump retained as historical material
#
# Anything else is branding creep and fails this audit. R/neha_path_utils.R is exempt because it
# IS the documented deprecated shim (see tests/test_deprecated_path_utils_shim.R), and
# 06_manuscript_figure_revision/ is already excluded above as frozen provenance.
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
  "neha_proteomics_audit_fast.txt",
  # Deliberate stale-name guard: a test asserts the pre-rename override is NOT used.
  "NEHA_RANK_ABUNDANCE_INPUT_DIR",
  # The deprecated shim's own filename, which documentation legitimately refers to. Whether any
  # active file actually SOURCES it is a separate contract, asserted by
  # tests/test_deprecated_path_utils_shim.R.
  "neha_path_utils.R"
)

strip_legacy_exceptions <- function(text) {
  for (literal in legacy_collaborator_exceptions) text <- gsub(literal, "", text, fixed = TRUE)
  text
}

# Committed controls, mirroring the output-label rules above, so changes to this policy are
# observable rather than silent.
branding_controls <- data.frame(
  text = c(
    "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler",
    "neha_protigy_input_animal_level_primary.gct",
    'dataset_project = "Neha"',
    "# Neha proteomics workflow",
    "NEHA_PCA_OUTPUT_ROOT",
    "validate_neha_pca_animal_input",
    "source(file.path(repo_root, 'R', 'project_path_utils.R'))"
  ),
  expected_rejected = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
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
