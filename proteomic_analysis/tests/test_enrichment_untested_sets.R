#!/usr/bin/env Rscript

# PROTEOMICS_ENRICHMENT_NPERM_SIMPLE parsing and the audit of gene sets that fgsea returns
# without a P value (01_clusterProfiler.r writes them to audits/<comparison>/gsea_untested_sets.csv).

repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "project_path_utils.R"))
source(file.path(repo_root, "R", "animal_level_enrichment_utils.R"))

assert_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)
assert_equal <- function(left, right, message) if (!identical(left, right)) stop(message, call. = FALSE)
assert_error <- function(expression, pattern, message) {
  observed <- tryCatch({ force(expression); NULL }, error = conditionMessage)
  if (is.null(observed) || !grepl(pattern, observed, ignore.case = TRUE)) {
    stop(message, if (!is.null(observed)) paste0(" Observed: ", observed) else " No error was raised.", call. = FALSE)
  }
}

# ---- option parsing -------------------------------------------------------------------------
fixture_root <- tempfile("untested_sets_fixture_")
mapped_root <- file.path(fixture_root, "mapped")
resolve_fixture <- function() resolve_enrichment_config(
  mapped_root = mapped_root,
  mapped_index = file.path(mapped_root, "indexMappedComparisons.csv"),
  output_root = file.path(fixture_root, "enrichment")
)
old_value <- Sys.getenv("PROTEOMICS_ENRICHMENT_NPERM_SIMPLE", unset = NA_character_)
on.exit({
  if (is.na(old_value)) Sys.unsetenv("PROTEOMICS_ENRICHMENT_NPERM_SIMPLE") else Sys.setenv(PROTEOMICS_ENRICHMENT_NPERM_SIMPLE = old_value)
}, add = TRUE)

Sys.unsetenv("PROTEOMICS_ENRICHMENT_NPERM_SIMPLE")
assert_equal(resolve_fixture()$gsea_n_perm_simple, NA_integer_, "Unset nPermSimple must keep the fgsea default (NA).")
Sys.setenv(PROTEOMICS_ENRICHMENT_NPERM_SIMPLE = " ")
assert_equal(resolve_fixture()$gsea_n_perm_simple, NA_integer_, "Blank nPermSimple must keep the fgsea default (NA).")
Sys.setenv(PROTEOMICS_ENRICHMENT_NPERM_SIMPLE = "100000")
assert_equal(resolve_fixture()$gsea_n_perm_simple, 100000L, "nPermSimple = 100000 was not parsed as an integer.")
for (bad in c("abc", "0", "-5", "1.5", "1e12")) {
  Sys.setenv(PROTEOMICS_ENRICHMENT_NPERM_SIMPLE = bad)
  assert_error(resolve_fixture(), "PROTEOMICS_ENRICHMENT_NPERM_SIMPLE", paste0("Invalid nPermSimple '", bad, "' must stop with an error."))
}
Sys.unsetenv("PROTEOMICS_ENRICHMENT_NPERM_SIMPLE")

# ---- untested-set listing on a hand-made case -----------------------------------------------
gene_sets <- list(
  A = c("g1", "g2", "g3", "g2"),        # 3 unique members in the rank
  B = c("g1", "x9"),                    # 1 member in the rank: below the minimum
  C = c("g1", "g2", "g3", "g4", "g5"),  # 5 members: above the maximum
  D = character(),
  E = c("g6", "g4", "g5")               # 3 members
)
rank_ids <- paste0("g", 1:6)
untested <- enrichment_untested_gsea_sets(gene_sets, rank_ids, tested_ids = "A", min_gs_size = 2, max_gs_size = 4)
assert_equal(untested$ID, "E", "Only the in-limit set missing from the table may be listed.")
assert_equal(untested$setSize, 3L, "setSize must count unique members present in the rank.")
none_tested <- enrichment_untested_gsea_sets(gene_sets, rank_ids, tested_ids = character(), min_gs_size = 2, max_gs_size = 4)
assert_equal(none_tested$ID, c("A", "E"), "Untested sets must be every in-limit set, sorted by ID.")
assert_equal(none_tested$setSize, c(3L, 3L), "Duplicated members must not inflate setSize.")
assert_equal(nrow(enrichment_untested_gsea_sets(gene_sets, rank_ids, c("A", "E"), 2, 4)), 0L, "Fully tested sets must give no rows.")

# ---- warning count ----------------------------------------------------------------------------
fgsea_warning <- paste0(
  "There were 31 pathways for which P-values were not calculated properly due to unbalanced ",
  "(positive and negative) gene-level statistic values. For such pathways pval, padj, NES, log2err ",
  "are set to NA. You can try to increase the value of the argument nPermSimple (for example set it nPermSimple = 10000)"
)
assert_equal(enrichment_unbalanced_set_count(character()), 0L, "No warning must count as zero untested sets.")
assert_equal(
  enrichment_unbalanced_set_count(c("For some of the pathways the P-values were likely overestimated.", fgsea_warning)),
  31L, "The fgsea unbalanced-set count was not read."
)
assert_equal(enrichment_unbalanced_set_count(sub("31", "2", fgsea_warning)), 2L, "A one-digit count was not read.")

# ---- against fgsea itself -----------------------------------------------------------------
# A rank with almost only positive values and a large set at its bottom: fewer than 10 of the
# random sets reach an ES <= 0, so fgsea returns no P value for that set.
if (requireNamespace("fgsea", quietly = TRUE) && requireNamespace("BiocParallel", quietly = TRUE)) {
  stats <- c(seq(5, 0.01, length.out = 980), seq(-0.01, -0.2, length.out = 20))
  names(stats) <- paste0("p", seq_along(stats))
  pathways <- list(
    bottom_300 = names(stats)[701:1000],
    bottom_250 = names(stats)[751:1000],
    top_40 = names(stats)[1:40],
    spread = names(stats)[seq(3, 1000, by = 9)],
    too_small = names(stats)[c(5, 900)]
  )
  set.seed(20261125)
  warnings <- character()
  result <- withCallingHandlers(
    fgsea::fgsea(pathways, stats, minSize = 10, maxSize = 500, BPPARAM = BiocParallel::SerialParam()),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  without_p <- sort(result$pathway[is.na(result$pval)], method = "radix")
  assert_true(length(without_p) >= 1L, "The fixture did not produce an unbalanced set; the check would be empty.")
  # The table DOSE writes keeps only the sets with a P value.
  untested <- enrichment_untested_gsea_sets(pathways, names(stats), result$pathway[!is.na(result$pval)], 10, 500)
  assert_equal(untested$ID, without_p, "The audit must list exactly the sets fgsea left without a P value.")
  assert_equal(untested$setSize, as.integer(result$size[match(untested$ID, result$pathway)]), "Audit setSize must match fgsea's size.")
  assert_equal(nrow(untested), enrichment_unbalanced_set_count(warnings), "Audit rows must match the fgsea warning count.")
} else {
  cat("fgsea or BiocParallel not installed; skipped the fgsea comparison.\n")
}

cat("All untested-set audit tests passed.\n")
