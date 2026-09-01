#!/usr/bin/env Rscript

# Contract: exactly 12 primary contrasts, generated from the canonical helper, and
# secondary analyses kept strictly out of the primary table.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("07_publication_release", "tests", "test_release_contrasts.R")
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
pc_path <- file.path(OUT_ROOT, "metadata", "primary_contrast_manifest.tsv")
sec_path <- file.path(OUT_ROOT, "metadata", "secondary_analysis_manifest.tsv")
if (!file.exists(pc_path) || !file.exists(sec_path)) {
  cat("No built publication release at", OUT_ROOT, "\n")
  cat("Skipping publication release contrast contracts.\n")
  quit(save = "no", status = 0L)
}

read_tsv <- function(p) read.delim(p, sep = "\t", quote = "", stringsAsFactors = FALSE,
                                   check.names = FALSE)
pc <- read_tsv(pc_path)
sec <- read_tsv(sec_path)
inv <- RELEASE_DESIGN_INVARIANTS

cat("=== the primary contrast set ===\n")
expect(nrow(pc) == inv$n_primary_contrasts,
       sprintf("exactly 12 primary contrasts (got %d)", nrow(pc)))
expect(!anyDuplicated(pc$canonical_comparison), "canonical_comparison is unique")
expect(setequal(pc$sample_class, sample_classes), "all four sample classes are covered")
expect(setequal(pc$contrast_family,
                c("paired_cno_vs_paired_veh", "paired_veh_vs_unpaired_veh",
                  "unpaired_cno_vs_unpaired_veh")),
       "the three canonical contrast families are covered")
expect(all(table(pc$sample_class) == 3L),
       "three contrast families per sample class")
expect(all(table(pc$contrast_family) == 4L),
       "four sample classes per contrast family")

cat("\n=== generated from the canonical helper, not re-encoded ===\n")
canonical <- primary_contrast_manifest()
expect(setequal(pc$canonical_comparison, canonical$canonical_comparison),
       "canonical_comparison reproduces primary_contrast_manifest()")
expect(setequal(pc$canonical_contrast, canonical$canonical_contrast),
       "canonical_contrast reproduces primary_contrast_manifest()")
expect(setequal(pc$historical_comparison_alias, canonical$historical_comparison_name),
       "historical aliases reproduce the canonical helper")
i <- match(pc$canonical_comparison, canonical$canonical_comparison)
expect(identical(as.character(pc$numerator_condition),
                 as.character(canonical$case_condition[i])),
       "numerator condition matches the canonical helper")
expect(identical(as.character(pc$denominator_condition),
                 as.character(canonical$reference_condition[i])),
       "denominator condition matches the canonical helper")

cat("\n=== every arm is 3 animals ===\n")
expect(all(as.integer(pc$n_numerator_animals) == inv$n_animals_per_stratum),
       "3 animals in every numerator arm")
expect(all(as.integer(pc$n_denominator_animals) == inv$n_animals_per_stratum),
       "3 animals in every denominator arm")
expect(all(pc$analysis_unit == "animal"), "analysis_unit is animal on every row")
expect(all(pc$primary_or_secondary == "primary"), "every row is labelled primary")
expect(all(pc$gsea_ranking_statistic == "moderated_t"),
       "the recorded GSEA ranking statistic is the moderated t")

cat("\n=== conditions are canonical and correctly oriented ===\n")
expect(all(pc$numerator_condition %in% condition_levels) &&
         all(pc$denominator_condition %in% condition_levels),
       "both arms name canonical conditions")
fam <- pc$contrast_family
expect(all(pc$numerator_condition[fam == "paired_cno_vs_paired_veh"] == "paired_cno") &&
         all(pc$denominator_condition[fam == "paired_cno_vs_paired_veh"] == "paired_veh"),
       "paired_cno_vs_paired_veh is oriented cno over veh")
expect(all(pc$numerator_condition[fam == "paired_veh_vs_unpaired_veh"] == "paired_veh") &&
         all(pc$denominator_condition[fam == "paired_veh_vs_unpaired_veh"] == "unpaired_veh"),
       "paired_veh_vs_unpaired_veh is oriented paired over unpaired")
expect(all(pc$numerator_condition[fam == "unpaired_cno_vs_unpaired_veh"] == "unpaired_cno") &&
         all(pc$denominator_condition[fam == "unpaired_cno_vs_unpaired_veh"] == "unpaired_veh"),
       "unpaired_cno_vs_unpaired_veh is oriented cno over veh")
expect(all(nzchar(pc$interpretation)), "every contrast carries an interpretation")
learning <- pc$interpretation[fam == "paired_veh_vs_unpaired_veh"]
expect(all(grepl("collection.plate", learning)),
       "the learning contrast carries the collection-plate identifiability caveat")

cat("\n=== secondary analyses are separated ===\n")
expect(nrow(sec) > 0L, sprintf("secondary manifest is populated (%d entries)", nrow(sec)))
expect(all(sec$primary_or_secondary == "secondary"),
       "every secondary row is labelled secondary")
expect(all(toupper(sec$included_in_primary_contrast_manifest) == "FALSE"),
       "no secondary analysis claims membership of the primary manifest")
expect(!any(sec$analysis_id %in% pc$canonical_comparison),
       "no secondary analysis_id collides with a primary comparison")
expect(any(grepl("cross_compartment", sec$analysis_id)),
       "the cross-compartment (Supplementary E) analysis is recorded as secondary")
expect(any(grepl("sensitivity", sec$analysis_id)),
       "the log2FC-ranked GSEA sensitivity analysis is recorded as secondary")
expect(all(nzchar(sec$caveat)), "every secondary analysis carries a caveat")

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Publication release contrast contracts failed: %d", failures), call. = FALSE)
}
cat("All publication release contrast contracts hold.\n")
