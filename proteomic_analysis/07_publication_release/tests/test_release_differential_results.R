#!/usr/bin/env Rscript

# Contract: the published differential statistics are the canonical ones, unchanged.
#
# The central assertion re-reads every canonical split table and requires bit-level
# equality with the exported values. If that ever fails, the release has recomputed or
# rounded something and must not be published.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("07_publication_release", "tests", "test_release_differential_results.R")
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
da_path <- file.path(OUT_ROOT, "differential_analysis", "primary_differential_proteins.tsv.gz")
sum_path <- file.path(OUT_ROOT, "differential_analysis", "primary_differential_summary.tsv")
if (!file.exists(da_path) || !file.exists(sum_path)) {
  cat("No built publication release at", OUT_ROOT, "\n")
  cat("Skipping publication release differential contracts.\n")
  quit(save = "no", status = 0L)
}

da <- read.delim(gzfile(da_path), sep = "\t", quote = "", stringsAsFactors = FALSE,
                 check.names = FALSE)
ds <- read.delim(sum_path, sep = "\t", quote = "", stringsAsFactors = FALSE,
                 check.names = FALSE)
inv <- RELEASE_DESIGN_INVARIANTS

cat("=== shape ===\n")
expect(length(unique(da$canonical_comparison)) == inv$n_primary_contrasts,
       sprintf("12 comparison blocks (got %d)", length(unique(da$canonical_comparison))))
expect(all(table(da$canonical_comparison) == inv$n_proteins_statistical),
       "every block holds all 5349 tested protein groups")
expect(nrow(da) == inv$n_primary_contrasts * inv$n_proteins_statistical,
       sprintf("64188 rows total (got %d)", nrow(da)))
expect(nrow(ds) == inv$n_primary_contrasts, "summary has one row per comparison")

cat("\n=== the statistical unit is animal everywhere ===\n")
expect(all(da$analysis_unit == "animal"), "analysis_unit is animal on every row")
expect(all(as.integer(da$n_numerator_animals) == 3L),
       "no row claims more than 3 animals in the numerator")
expect(all(as.integer(da$n_denominator_animals) == 3L),
       "no row claims more than 3 animals in the denominator")
expect(all(da$primary_or_secondary == "primary"), "every row is labelled primary")

cat("\n=== field naming is not misleading ===\n")
expect(!release_column_is_misleading_gene_symbol(da$gene_symbol),
       "the public gene_symbol column does not contain UniProt identifiers")
expect(mean(release_uniprot_like(unique(da$protein_group_id))) > 0.9,
       "protein_group_id holds UniProt entry names, as documented")
expect("uniprot_accession" %in% names(da) && "protein_description" %in% names(da),
       "accession and description are published as separate, correctly named fields")
expect(!"Description" %in% names(da),
       "the misleading source column name `Description` is not republished")

cat("\n=== exported statistics equal the canonical split tables ===\n")
SPLIT_FORWARD <- file.path(DATA_ROOT, "02_data", "animal_level", "split", "forward")
if (!dir.exists(SPLIT_FORWARD)) {
  cat("  [skip]  canonical split tables unreachable; equality not checked\n")
} else {
  CHECK <- c(effect_size_sd_units = "log2fc", average_standardized_abundance = "aveExpr", moderated_t = "t",
             P.Value = "pval", adj.P.Val = "padj", B_log_odds = "B",
             neg_log10_P_value = "logpval", signed_neg_log10_P_value = "sign.logP")
  bad <- character(0)
  for (cc in unique(da$canonical_contrast)) {
    p <- file.path(SPLIT_FORWARD, paste0(cc, ".csv"))
    if (!file.exists(p)) { bad <- c(bad, paste("missing", cc)); next }
    d <- release_read_csv(p)
    ex <- da[da$canonical_contrast == cc, , drop = FALSE]
    if (!identical(as.character(ex$protein_group_id), as.character(d$gene_symbol))) {
      bad <- c(bad, paste("order", cc)); next
    }
    for (pub in names(CHECK)) {
      if (!identical(as.numeric(ex[[pub]]), as.numeric(d[[CHECK[[pub]]]]))) {
        bad <- c(bad, paste(pub, cc))
      }
    }
    if (!identical(as.logical(ex$significant_fdr_0_05), as.logical(d$significant))) {
      bad <- c(bad, paste("significance flag", cc))
    }
  }
  expect(length(bad) == 0L,
         sprintf("all 12 comparisons x 8 statistics are bit-identical to canonical%s",
                 if (length(bad)) paste0(" -- ", paste(head(bad, 3), collapse = "; ")) else ""))
}

cat("\n=== summary counts are consistent with the exported rows ===\n")
ok_counts <- vapply(seq_len(nrow(ds)), function(i) {
  sub <- da[da$canonical_comparison == ds$canonical_comparison[i], , drop = FALSE]
  sig <- as.logical(sub$significant_fdr_0_05)
  mapped <- sub$id_mapping_status == "mapped"
  sum(sig, na.rm = TRUE) == as.integer(ds$n_significant_fdr_0_05[i]) &&
    sum(sig & mapped, na.rm = TRUE) == as.integer(ds$n_significant_fdr_0_05_mapped_only[i]) &&
    sum(mapped) == as.integer(ds$n_proteins_tested_mapped_only[i])
}, logical(1))
expect(all(ok_counts), "full-set and mapped-only counts both reconcile")
expect(all(as.integer(ds$n_proteins_tested_mapped_only) == inv$n_proteins_mapped),
       "5327 mapped protein groups in every comparison")
expect(all(as.integer(ds$n_significant_fdr_0_05) >=
             as.integer(ds$n_significant_fdr_0_05_mapped_only)),
       "the mapped-only significant count never exceeds the full-set count")

cat("\n=== every protein group is annotated ===\n")
expect(!any(is.na(da$protein_group_id)) && all(nzchar(da$protein_group_id)),
       "protein_group_id is populated on every row")
expect(all(da$id_mapping_status %in% c("mapped", "unmapped")),
       "id_mapping_status is mapped or unmapped on every row")
expect(sum(da$id_mapping_status == "unmapped") ==
         inv$n_primary_contrasts * (inv$n_proteins_statistical - inv$n_proteins_mapped),
       "22 unmapped protein groups per comparison are retained, not dropped")

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Publication release differential contracts failed: %d", failures),
       call. = FALSE)
}
cat("All publication release differential contracts hold.\n")
