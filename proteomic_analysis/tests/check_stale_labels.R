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

if (length(failures) > 0) {
  stop(paste(c("Stale-label audit failed:", failures), collapse = "\n"), call. = FALSE)
}

message(
  "Stale-label controls passed: ", sum(control_cases$expected_rejected), " must-reject, ",
  sum(!control_cases$expected_rejected), " must-allow; active-code audit passed."
)
