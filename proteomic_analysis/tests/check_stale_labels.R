args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1) sub("^--file=", "", file_arg) else file.path("tests", "check_stale_labels.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(repo_root)) repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

all_files <- list.files(repo_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
# Archived, non-runnable snapshots are excluded the same way /legacy/ is: this audit
# checks ACTIVE code. 06_manuscript_figure_revision/ is a verbatim provenance copy of
# already-executed code (see its README) and must not be edited to satisfy a lint.
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
# Narrowed 2026-08-27. The original pattern was "\\bNature\\b|publication|manuscript".
# Its purpose is to stop OUTPUT ARTEFACTS being named after a target venue ("Nature_Fig3.svg",
# "publication_ready_theme"), which goes stale the moment the submission target changes.
# Since the repository gained 06_manuscript_figure_revision/, the bare words "manuscript" and
# "publication" became load-bearing domain vocabulary -- the stage exists precisely to regenerate
# manuscript figure panels, and README/CANONICAL_OUTPUTS must be able to say so. The bare words
# produced 15 false positives and 0 true positives ("Nature" appears nowhere in the repository).
# The venue-naming smells below are still caught, as is theme_nature via its own rule.
#
# The journal name is checked case-SENSITIVELY and without \b on the right, because \b treats
# "_" as a word character: the original "\\bNature\\b" never matched "Nature_Fig3.svg", which is
# the exact artefact-naming case the rule exists to catch. Case sensitivity keeps ordinary prose
# ("the nature of the effect") from tripping it.
add_hits("journal-named output artefacts", "(?<![A-Za-z])Nature(?![a-z])", ignore_case = FALSE)
add_hits("drafting-language output names/comments", "publication[- _]?(ready|quality|figure)")
add_hits("old plotting helper names", "theme_nature|theme_nature_qc")

if (length(failures) > 0) {
  stop(paste(c("Stale-label audit failed:", failures), collapse = "\n"), call. = FALSE)
}

message("Stale-label audit passed for active code.")
