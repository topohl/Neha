#!/usr/bin/env Rscript

# Contract: the published enrichment tables are the canonical ones, canonical and
# sensitivity rows are never pooled, and every row carries enough context to stand alone.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("07_publication_release", "tests", "test_release_enrichment.R")
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
need <- file.path(OUT_ROOT, "enrichment",
                  c("primary_GSEA_GO_BP.tsv.gz", "primary_GSEA_KEGG.tsv.gz",
                    "primary_ORA_GO_BP.tsv.gz", "GSEA_log2FC_sensitivity.tsv.gz",
                    "primary_EWCE.tsv.gz", "enrichment_coverage.tsv"))
if (!all(file.exists(need))) {
  cat("No built publication release at", OUT_ROOT, "\n")
  cat("Skipping publication release enrichment contracts.\n")
  quit(save = "no", status = 0L)
}

rd <- function(p) {
  if (grepl("[.]gz$", p)) {
    read.delim(gzfile(p), sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.delim(p, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  }
}
gsea_go <- rd(need[1]); gsea_kegg <- rd(need[2]); ora <- rd(need[3])
sens <- rd(need[4]); ewce <- rd(need[5]); coverage <- rd(need[6])
inv <- RELEASE_DESIGN_INVARIANTS

cat("=== coverage ===\n")
expect(length(unique(gsea_go$canonical_comparison)) == inv$n_primary_contrasts,
       "GSEA GO-BP covers all 12 comparisons")
expect(length(unique(gsea_kegg$canonical_comparison)) == inv$n_primary_contrasts,
       "GSEA KEGG covers all 12 comparisons")
expect(length(unique(sens$canonical_comparison)) == inv$n_primary_contrasts,
       "the sensitivity table covers all 12 comparisons")
expect(all(as.logical(coverage$agrees_with_run_record)),
       "every exported row count equals the count the canonical run recorded")
expect(any(as.integer(coverage$rows_exported) == 0L),
       "empty analysis blocks are recorded rather than omitted silently")

cat("\n=== canonical ranks by the moderated t; log2FC is sensitivity only ===\n")
expect(all(gsea_go$rank_statistic == "moderated_t"),
       "canonical GSEA GO-BP is ranked by the moderated t")
expect(all(gsea_kegg$rank_statistic == "moderated_t"),
       "canonical GSEA KEGG is ranked by the moderated t")
expect(all(gsea_go$analysis_role == "canonical") &&
         all(gsea_kegg$analysis_role == "canonical"),
       "no sensitivity row appears in a canonical file")
expect(all(sens$rank_statistic == "log2fc") && all(sens$analysis_role == "sensitivity"),
       "the sensitivity file holds only log2FC-ranked sensitivity rows")
expect(!any(gsea_go$analysis_role == "sensitivity") &&
         !any(sens$analysis_role == "canonical"),
       "canonical and sensitivity rows are never pooled in one file")
expect(all(ora$rank_statistic == "not_applicable"),
       "ORA rows declare no ranking statistic")

cat("\n=== every row stands alone ===\n")
context_cols <- c("sample_class", "contrast_family", "canonical_comparison",
                  "numerator_condition", "denominator_condition", "n_numerator_animals",
                  "n_denominator_animals", "analysis_unit", "analysis", "ontology",
                  "term_id", "term_name", "p_value", "adjusted_p_value", "direction",
                  "rank_statistic", "analysis_role")
for (nm in c("gsea_go", "gsea_kegg", "ora", "sens")) {
  d <- get(nm)
  missing <- setdiff(context_cols, names(d))
  expect(length(missing) == 0L,
         sprintf("%s carries every required context column%s", nm,
                 if (length(missing)) paste0(" (missing: ", paste(missing, collapse = ", "), ")")
                 else ""))
  expect(all(nzchar(d$canonical_comparison)) && all(nzchar(d$analysis_role)),
         sprintf("%s has no blank comparison or analysis_role", nm))
  expect(all(d$analysis_unit == "animal"), sprintf("%s is animal-level throughout", nm))
  expect(all(as.integer(d$n_numerator_animals) == 3L),
         sprintf("%s never claims more than 3 animals per arm", nm))
}
expect(all(nzchar(ora$core_enrichment_definition)) &&
         all(nzchar(gsea_go$core_enrichment_definition)),
       "core_enrichment carries an explicit definition per row type")
expect(length(unique(gsea_go$core_enrichment_definition)) == 1L &&
         length(unique(ora$core_enrichment_definition)) == 1L &&
         gsea_go$core_enrichment_definition[1] != ora$core_enrichment_definition[1],
       "GSEA and ORA gene membership are defined differently and say so")

cat("\n=== ORA query lists ===\n")
expect(all(c("fdr_significant_all", "fdr_significant_higher_in_numerator",
             "fdr_significant_higher_in_denominator") %in% ora$query_list),
       "the three FDR-based ORA query lists are present")
expect(any(ora$analysis_role == "alternative_query_list"),
       "the non-FDR query list is labelled as an alternative, not as canonical")
expect(all(is.na(gsea_go$query_list) | !nzchar(gsea_go$query_list)),
       "GSEA rows carry no ORA query list")

cat("\n=== EWCE ===\n")
expect(all(ewce$analysis == "EWCE"), "EWCE rows are labelled EWCE")
expect(all(ewce$ewce_analysis_type %in% c("Baseline", "Differential")),
       "EWCE analysis type is Baseline or Differential")
expect(all(ewce$primary_or_secondary[ewce$ewce_analysis_type == "Baseline"] == "secondary"),
       "EWCE baseline rows are secondary, not treatment contrasts")
expect(all(ewce$analysis_role[ewce$is_primary_setting == "TRUE"] == "canonical"),
       "only the primary top-N / annotation-level setting is canonical")
expect(any(ewce$analysis_role == "sensitivity"),
       "non-primary EWCE settings are retained as sensitivity")
diff_rows <- ewce$ewce_analysis_type == "Differential"
expect(all(ewce$canonical_comparison[diff_rows] %in% gsea_go$canonical_comparison),
       "EWCE differential rows name one of the 12 primary comparisons")

cat("\n=== values equal the canonical enrichment files ===\n")
ENRICH <- file.path(DATA_ROOT, "03_output", "enrichment",
                    "enrichment_t_rank_validation_20260825", "per_comparison")
if (!dir.exists(ENRICH)) {
  cat("  [skip]  canonical enrichment folder unreachable; equality not checked\n")
} else {
  bad <- character(0); n <- 0L
  for (cmp in unique(gsea_go$canonical_comparison)) {
    p <- file.path(ENRICH, cmp, "GSEA_GO_BP.csv")
    if (!file.exists(p)) { bad <- c(bad, paste("missing", cmp)); next }
    d <- release_read_csv(p)
    ex <- gsea_go[gsea_go$canonical_comparison == cmp, , drop = FALSE]
    if (!identical(ex$term_id, as.character(d$ID))) { bad <- c(bad, paste("order", cmp)); next }
    for (pair in list(c("NES", "NES"), c("enrichment_score", "enrichmentScore"),
                      c("p_value", "pvalue"), c("adjusted_p_value", "p.adjust"))) {
      if (!identical(as.numeric(ex[[pair[[1]]]]), as.numeric(d[[pair[[2]]]]))) {
        bad <- c(bad, paste(pair[[1]], cmp))
      }
    }
    n <- n + 1L
  }
  expect(length(bad) == 0L,
         sprintf("GSEA GO-BP values are bit-identical to canonical in %d comparisons%s", n,
                 if (length(bad)) paste0(" -- ", paste(head(bad, 3), collapse = "; ")) else ""))
}

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Publication release enrichment contracts failed: %d", failures), call. = FALSE)
}
cat("All publication release enrichment contracts hold.\n")
