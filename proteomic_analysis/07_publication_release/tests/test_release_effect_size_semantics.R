#!/usr/bin/env Rscript

# Contract: the published effect size is named for what it is, is bit-identical to the
# canonical `logFC` it was copied from, and is nowhere described as a log2 fold change.
#
# The central assertion re-reads all 12 canonical split tables and requires every one of
# the 64,188 exported effect values to equal its source `log2fc` (ProTigy `logFC`) exactly.
# That is what makes "renamed, not recomputed" checkable rather than asserted.
#
# Internal provenance is deliberately NOT scrubbed: canonical filenames and the
# `rank_statistic` value keep the historical `log2fc` token, and the tests below require
# that they do, because renaming them would break the link to the canonical run.
#
# Skips cleanly when no built release is reachable, so it is safe in CI.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("07_publication_release", "tests", "test_release_effect_size_semantics.R")
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/",
                           mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "project_path_utils.R"))
LAYER <- file.path(repo_root, "07_publication_release")
source(file.path(LAYER, "R", "release_utils.R"))
source(file.path(LAYER, "R", "release_validation.R"))

failures <- 0L
expect <- function(condition, message) {
  if (!isTRUE(condition)) {
    failures <<- failures + 1L
    cat("  [FAIL]", message, "\n")
  } else {
    cat("  [ok]  ", message, "\n")
  }
}

OUT_ROOT <- release_output_root()
DATA_ROOT <- release_data_root()
da_path <- file.path(OUT_ROOT, "differential_analysis",
                     "primary_differential_proteins.tsv.gz")
sum_path <- file.path(OUT_ROOT, "differential_analysis",
                      "primary_differential_summary.tsv")
if (!file.exists(da_path) || !file.exists(sum_path)) {
  cat("No built publication release at", OUT_ROOT, "\n")
  cat("Skipping publication release effect-size contracts.\n")
  quit(save = "no", status = 0L)
}

read_tsv <- function(p) {
  if (grepl("[.]gz$", p)) {
    read.delim(gzfile(p), sep = "\t", quote = "", stringsAsFactors = FALSE,
               check.names = FALSE)
  } else {
    read.delim(p, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  }
}
da <- read_tsv(da_path)
ds <- read_tsv(sum_path)
inv <- RELEASE_DESIGN_INVARIANTS

# --------------------------------------------------------------------------------------
cat("=== the public field is named correctly ===\n")
# --------------------------------------------------------------------------------------
expect(RELEASE_EFFECT_SIZE$public_field %in% names(da),
       paste("the effect size is published as", RELEASE_EFFECT_SIZE$public_field))
expect(!any(tolower(names(da)) %in% c("logfc", "log2fc", "log2foldchange",
                                      "log2_fold_change", "fc", "foldchange")),
       "no published column is named logFC, log2FC, log2fc or foldchange")
expect(all(c("moderated_t", "P.Value", "adj.P.Val") %in% names(da)),
       "the moderated t, p-value and adjusted p-value are published alongside it")

# --------------------------------------------------------------------------------------
cat("\n=== J. internal provenance is preserved, not scrubbed ===\n")
# --------------------------------------------------------------------------------------
expect("source_statistic_field" %in% names(da) &&
         all(da$source_statistic_field == "logFC"),
       "source_statistic_field records the canonical source field as logFC")
expect("effect_size_source_column" %in% names(da) &&
         all(grepl("logFC", da$effect_size_source_column, fixed = TRUE)),
       "effect_size_source_column names the ProTigy logFC column")
expect(all(grepl("log2fc", da$effect_size_source_column, fixed = TRUE)),
       "it also names the canonical split-table column log2fc")

# --------------------------------------------------------------------------------------
cat("\n=== G. exported effect values equal the canonical logFC exactly ===\n")
# --------------------------------------------------------------------------------------
SPLIT_FORWARD <- file.path(DATA_ROOT, "02_data", "animal_level", "split", "forward")
if (!dir.exists(SPLIT_FORWARD)) {
  cat("  [skip]  canonical split tables unreachable; exact equality not checked\n")
} else {
  bad <- character(0)
  n_checked <- 0L
  for (cc in unique(da$canonical_contrast)) {
    p <- file.path(SPLIT_FORWARD, paste0(cc, ".csv"))
    if (!file.exists(p)) { bad <- c(bad, paste("missing", cc)); next }
    d <- release_read_csv(p)
    ex <- da[da$canonical_contrast == cc, , drop = FALSE]
    if (!identical(as.character(ex$protein_group_id), as.character(d$gene_symbol))) {
      bad <- c(bad, paste("protein order", cc)); next
    }
    src <- as.numeric(d$log2fc)
    got <- as.numeric(ex[[RELEASE_EFFECT_SIZE$public_field]])
    if (!identical(got, src)) bad <- c(bad, paste("effect value", cc))
    n_checked <- n_checked + length(src)
  }
  expect(length(bad) == 0L,
         sprintf("all %d exported effect values are bit-identical to canonical logFC%s",
                 n_checked,
                 if (length(bad)) paste0(" -- ", paste(head(bad, 3), collapse = "; ")) else ""))
  expect(n_checked == inv$n_primary_contrasts * inv$n_proteins_statistical,
         sprintf("all 64188 rows were compared (got %d)", n_checked))
}

# --------------------------------------------------------------------------------------
cat("\n=== H. the release does not call these values log2 fold changes ===\n")
# --------------------------------------------------------------------------------------
expect(all(da$effect_size_units == RELEASE_EFFECT_SIZE$units),
       "the per-row units say per-protein standard deviation")
expect(!any(grepl("fold", da$effect_size_units, ignore.case = TRUE)),
       "the units string contains no fold-change wording")
expect(all(grepl("NOT a log2 fold change", ds$effect_size_definition, fixed = TRUE)),
       "the definition states explicitly that it is not a log2 fold change")
expect(!any(grepl("standardis(ed|ed) mean difference", ds$effect_size_definition)) &&
         !any(grepl("standardised mean difference", ds$effect_size_definition,
                    fixed = TRUE)) &&
         !any(grepl("standardized mean difference", ds$effect_size_definition,
                    fixed = TRUE)),
       "it is not called a standardised mean difference without qualification")

PUBLIC_DOCS <- c("README_DATA.md",
                 "metadata/data_dictionary.tsv",
                 "metadata/primary_contrast_manifest.tsv",
                 "metadata/secondary_analysis_manifest.tsv",
                 "differential_analysis/primary_differential_summary.tsv",
                 "editor_source_data/REVISION_PROTEOMICS_DATA_CHANGELOG.md",
                 "editor_source_data/MANUSCRIPT_TERMINOLOGY_ACTIONS.md",
                 "editor_source_data/figure_source_map.tsv",
                 "pride/README_PRIDE.md",
                 "pride/SDRF_MISSING_METADATA.md",
                 "provenance/UPSTREAM_PREPROCESSING_GAP.md")
total_hits <- 0L
for (doc in PUBLIC_DOCS) {
  p <- file.path(OUT_ROOT, doc)
  if (!file.exists(p)) next
  h <- release_fold_change_mislabels(readLines(p, warn = FALSE))
  if (nrow(h)) {
    total_hits <- total_hits + nrow(h)
    for (i in seq_len(nrow(h))) {
      cat("        ", doc, ":", h$line[i], " ", substr(h$text[i], 1, 110), "\n", sep = "")
    }
  }
}
expect(total_hits == 0L,
       sprintf("no reader-facing document describes the effect size as a fold change (%d hit(s))",
               total_hits))

# --------------------------------------------------------------------------------------
cat("\n=== the methods-safe definition is available to copy ===\n")
# --------------------------------------------------------------------------------------
readme <- file.path(OUT_ROOT, "README_DATA.md")
if (file.exists(readme)) {
  expect(any(grepl(RELEASE_EFFECT_SIZE$methods_sentence, readLines(readme, warn = FALSE),
                   fixed = TRUE)),
         "README_DATA.md carries the methods-safe definition verbatim")
}
expect(!grepl("z-score", RELEASE_EFFECT_SIZE$methods_sentence, fixed = TRUE) &&
         !grepl("z-scored", RELEASE_EFFECT_SIZE$methods_sentence, fixed = TRUE),
       "the methods sentence does not claim z-scoring")

# The standardisation claim must match what the matrix actually shows.
meas_path <- file.path(OUT_ROOT, "processed_data",
                       "protein_abundance_measurement_level.tsv.gz")
if (!file.exists(meas_path)) {
  cat("  [skip]  measurement-level matrix not built; standardisation not re-measured\n")
} else {
  m <- read_tsv(meas_path)
  mat <- as.matrix(m[, setdiff(names(m), "protein_group_id"), drop = FALSE])
  std <- release_standardization_evidence(mat)
  expect(isTRUE(std$approximately_standardized),
         sprintf("the matrix is per-protein standardised (max |row mean| %.2g, max |sd-1| %.2g)",
                 std$max_abs_row_mean, std$max_abs_row_sd_minus_1))
  expect(isFALSE(std$exact_zscore),
         "it is NOT exactly z-scored, which is why the release does not say so")
}

# --------------------------------------------------------------------------------------
cat("\n=== I. canonical GSEA still ranks by the moderated t ===\n")
# --------------------------------------------------------------------------------------
for (f in c("primary_GSEA_GO_BP.tsv.gz", "primary_GSEA_KEGG.tsv.gz")) {
  p <- file.path(OUT_ROOT, "enrichment", f)
  if (!file.exists(p)) { expect(FALSE, paste(f, "exists")); next }
  g <- read_tsv(p)
  expect(all(g$rank_statistic == "moderated_t"),
         paste(f, "is ranked by moderated_t on every row"))
  expect(all(g$analysis_role == "canonical"),
         paste(f, "carries analysis_role = canonical on every row"))
}
sens_path <- file.path(OUT_ROOT, "enrichment", "GSEA_log2FC_sensitivity.tsv.gz")
if (file.exists(sens_path)) {
  s <- read_tsv(sens_path)
  expect(all(s$analysis_role == "sensitivity"),
         "the sensitivity table is entirely analysis_role = sensitivity")
  expect(all(s$rank_statistic == "log2fc"),
         "the sensitivity rank_statistic keeps its canonical log2fc token (provenance)")
  expect(!any(s$rank_statistic == "moderated_t"),
         "no canonical row leaked into the sensitivity table")
}

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Publication release effect-size contracts failed: %d", failures),
       call. = FALSE)
}
cat("All publication release effect-size contracts hold.\n")
