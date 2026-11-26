#!/usr/bin/env Rscript

# Tests for 04_differential_expression_enrichment/05_pairing_cno_interaction.r, protein level only.
#   1. the interaction contrast weighs paired_cno, paired_veh, unpaired_cno, unpaired_veh as +1, -1, -1, +1;
#      its estimate is the paired CNO effect minus the unpaired CNO effect; adding pairing and CNO main
#      effects changes neither estimate nor t;
#   2. with the study's plate layout the contrast keeps one third of an additive plate difference;
#   3. the output guard rejects the pooled-run inferential folder, the synthesis folder and their parents;
#   4. on a 300-protein subset of the animal-level input (GSEA off), the stage writes the per-class
#      protein tables and the release-shaped summaries with the columns of final_inferential_checks_20260826,
#      fits residual df 8 per class, and refuses a second run into the same folder.
# Writes only to a temporary folder. Part 4 skips cleanly when the shared drive is not reachable.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("tests", "test_pairing_cno_interaction.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "project_path_utils.R"))
source(file.path(repo_root, "R", "animal_level_enrichment_utils.R"))
source(file.path(repo_root, "R", "inferential_checks_utils.R"))

if (!requireNamespace("limma", quietly = TRUE) || !requireNamespace("digest", quietly = TRUE)) {
  cat("limma or digest not installed; skipping Pairing x CNO interaction checks.\n")
  quit(save = "no", status = 0L)
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
group_mean_interaction <- function(values, conditions) {
  m <- vapply(condition_levels, function(cond) rowMeans(values[, conditions == cond, drop = FALSE]), numeric(nrow(values)))
  (m[, "paired_cno"] - m[, "paired_veh"]) - (m[, "unpaired_cno"] - m[, "unpaired_veh"])
}

cat("interaction contrast on synthetic data\n")
set.seed(20261125)
conditions <- rep(condition_levels, each = 3L)
plates <- c("Plate1", "Plate2", "Plate1", rep("Plate1", 3), rep("Plate2", 3), rep("Plate2", 3))
values <- matrix(stats::rnorm(200L * 12L), nrow = 200L, dimnames = list(paste0("p", 1:200), paste0("animal", 1:12)))
paired <- startsWith(conditions, "paired")
cno <- endsWith(conditions, "_cno")
values <- rbind(
  values,
  additive = values["p1", ] + 1.0 * paired + 2.0 * cno,
  interacting = values["p1", ] + 1.5 * (paired & cno),
  plate_shift = values["p1", ] + 0.9 * (plates == "Plate2")
)
fit <- fit_condition_contrast(values, conditions, PAIRING_CNO_INTERACTION_CONTRAST)
tab <- fit$table
expect(identical(fit$weights[condition_levels], c(paired_cno = 1, paired_veh = -1, unpaired_cno = -1, unpaired_veh = 1)),
       "contrast weights are +1 paired_cno, -1 paired_veh, -1 unpaired_cno, +1 unpaired_veh")
expect(fit$df_residual == 8, "12 animals in 4 groups give residual df 8")
expect(max(abs(tab$logFC - group_mean_interaction(values, conditions))) < 1e-12,
       "estimate = (paired_cno - paired_veh) - (unpaired_cno - unpaired_veh) of the group means")
expect(abs(tab["additive", "logFC"] - tab["p1", "logFC"]) < 1e-12 && abs(tab["additive", "t"] - tab["p1", "t"]) < 1e-8,
       "pairing and CNO main effects leave the interaction estimate and t unchanged")
expect(abs(tab["interacting", "logFC"] - tab["p1", "logFC"] - 1.5) < 1e-12,
       "a CNO effect present only in paired animals gives a positive interaction of that size")
expect(isTRUE(all.equal(tab$adj.P.Val, stats::p.adjust(tab$P.Value, "BH"))), "FDR is Benjamini-Hochberg within the contrast")

cat("collection plate\n")
weight <- contrast_plate_weight(fit$weights, conditions, plates)
expect(abs(weight - 1 / 3) < 1e-12, "the contrast keeps one third of an additive Plate2 - Plate1 difference")
expect(abs(tab["plate_shift", "logFC"] - tab["p1", "logFC"] - 0.9 * weight) < 1e-12, "a pure plate shift enters the estimate with that weight")

cat("output guard\n")
guard_error <- function(path) inherits(tryCatch(validate_inferential_output_root(path), error = function(e) e), "error")
old_root <- file.path(PROTEOMICS_DATA_ROOT, "03_output", "inferential_checks", "final_inferential_checks_20260826")
expect(guard_error(old_root), "rejects final_inferential_checks_20260826")
expect(guard_error(file.path(old_root, "posthoc_pairing_x_cno_interaction")), "rejects a folder inside it")
expect(guard_error(file.path(PROTEOMICS_DATA_ROOT, "03_output", "inferential_checks")), "rejects its parent")
expect(guard_error(file.path(PROTEOMICS_DATA_ROOT, "03_output", "synthesis", "animal_level_biological_synthesis_20260826")), "rejects the synthesis folder")
expect(guard_error(file.path(PROTEOMICS_DATA_ROOT, "02_data", "animal_level")), "rejects 02_data")
expect(!guard_error(file.path(PROTEOMICS_DATA_ROOT, "03_output", "inferential_checks", "within_compartment_20261124")), "accepts the per-compartment folder")
expect(!guard_error(file.path(tempdir(), "x")), "accepts a temporary folder")

cat("stage run on a 300-protein subset (no GSEA)\n")
input_gct <- file.path(PROTEOMICS_DATA_ROOT, "02_data", "animal_level", "input_gct", "neha_protigy_input_animal_level_primary.gct")
mapped_index <- file.path(PROTEOMICS_DATA_ROOT, "02_data", "animal_level", "within_compartment_20261124", "mapped", "indexMappedComparisons.csv")
design_file <- file.path(old_root, "design_identifiability_audit.csv")
if (!all(file.exists(c(input_gct, mapped_index, design_file)))) {
  cat("Shared data root not reachable; skipping the stage run.\n")
} else {
  scratch <- file.path(tempdir(), "pairing_cno_interaction_test")
  dir.create(scratch, showWarnings = FALSE)
  lines <- readLines(input_gct, warn = FALSE)
  dims <- scan(text = lines[[2]], what = integer(), quiet = TRUE)
  n_keep <- 300L
  dims[[1]] <- n_keep
  subset_gct <- file.path(scratch, "subset.gct")
  writeLines(c(lines[[1]], paste(dims, collapse = "\t"), lines[seq(3L, 3L + dims[[4]] + n_keep)]), subset_gct)
  output_root <- file.path(scratch, "out")

  run_stage <- function() {
    # Outputs stay in this session's temp folder; FORCE is set to false.
    stopifnot(project_path_is_within(output_root, tempdir()))
    vars <- c(PROTEOMICS_INFERENTIAL_INPUT_GCT = subset_gct, PROTEOMICS_INFERENTIAL_OUTPUT_ROOT = output_root,
              PROTEOMICS_INFERENTIAL_GSEA = "false", PROTEOMICS_INFERENTIAL_FORCE = "false")
    old <- Sys.getenv(names(vars), unset = NA, names = TRUE)
    do.call(Sys.setenv, as.list(vars))  # system2(env = ) is ignored on Windows
    on.exit({
      for (n in names(old)) if (is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, stats::setNames(list(old[[n]]), n))
    })
    log <- tempfile("stage_", tmpdir = scratch, fileext = ".log")
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      shQuote(file.path(repo_root, "04_differential_expression_enrichment", "05_pairing_cno_interaction.r")),
                      stdout = log, stderr = log)
    list(status = as.integer(status), log = log)
  }

  first <- run_stage()
  expect(identical(first$status, 0L), "stage completes with GSEA off")
  if (!identical(first$status, 0L)) cat(utils::tail(readLines(first$log, warn = FALSE), 15), sep = "\n")
  per_class <- file.path(output_root, "posthoc_pairing_x_cno_interaction", sample_classes, "interaction_protein_table.csv")
  summary_path <- file.path(output_root, "pairing_cno_interaction_protein_summary.csv")
  expect(all(file.exists(per_class)), "one interaction protein table per sample class")
  expect(file.exists(summary_path) && file.exists(file.path(output_root, "pairing_cno_interaction_top_proteins.csv")),
         "release-shaped protein summary and top-protein files written")
  expect(!file.exists(file.path(output_root, "pairing_cno_interaction_gsea_summary.csv")), "no GSEA summary without GSEA")
  if (all(file.exists(c(per_class, summary_path)))) {
    summary <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
    expect(identical(names(summary), c("sample_class", "n_tested", "n_nominal_p05", "pct_nominal_p05", "n_fdr05",
                                       "n_positive_interaction_fdr", "n_negative_interaction_fdr",
                                       "median_interaction_log2fc", "IQR_interaction_log2fc")),
           "protein summary has the columns of the 2026-08-26 file")
    top <- utils::read.csv(file.path(output_root, "pairing_cno_interaction_top_proteins.csv"), stringsAsFactors = FALSE)
    expect(identical(names(top), c("sample_class", "rank", "original_protein_id", "Description", "uniprot_accession",
                                   "logFC", "t", "P.Value", "adj.P.Val")),
           "top-protein file has the columns of the 2026-08-26 file")
    expect(identical(summary$sample_class, sample_classes) && all(summary$n_tested == n_keep), "four classes, all subset proteins tested")

    input <- read_animal_level_input_gct(subset_gct)
    worst_effect <- 0
    worst_fdr <- 0
    df_ok <- TRUE
    sorted_ok <- TRUE
    for (cls in sample_classes) {
      p <- utils::read.csv(file.path(output_root, "posthoc_pairing_x_cno_interaction", cls, "interaction_protein_table.csv"), stringsAsFactors = FALSE)
      a <- input$annotation[input$annotation$sample_class == cls, , drop = FALSE]
      expected <- group_mean_interaction(input$values[, a$column, drop = FALSE], a$condition)
      worst_effect <- max(worst_effect, abs(p$logFC - expected[p$original_protein_id]))
      worst_fdr <- max(worst_fdr, abs(p$adj.P.Val - stats::p.adjust(p$P.Value, "BH")))
      df_ok <- df_ok && all(p$df_residual == 8)
      sorted_ok <- sorted_ok && !is.unsorted(p$P.Value)
    }
    expect(worst_effect < 1e-12, paste0("protein effects equal the group-mean interaction (max diff ", signif(worst_effect, 3), ")"))
    expect(worst_fdr < 1e-12, "FDR is BH within each class")
    expect(df_ok, "per-class models have residual df 8")
    expect(sorted_ok, "protein tables are sorted by P value")
    diagnostics <- utils::read.csv(file.path(output_root, "pairing_cno_interaction_diagnostics.csv"), stringsAsFactors = FALSE)
    expect(all(abs(diagnostics$plate2_minus_plate1_weight - 1 / 3) < 1e-12), "plate weight 1/3 in every class")
    manifest <- utils::read.csv(file.path(output_root, "pairing_cno_interaction_run_manifest.csv"), stringsAsFactors = FALSE)
    expect(identical(manifest$detail[manifest$section == "input" & manifest$name == "input_gct"], project_file_sha256(subset_gct)),
           "run manifest records the input hash")

    before <- project_file_sha256(summary_path)
    second <- run_stage()
    expect(!identical(second$status, 0L), "a second run into the non-empty folder is refused without FORCE")
    expect(identical(project_file_sha256(summary_path), before), "the refused run leaves the outputs unchanged")
  }
}

if (failures > 0L) {
  cat(failures, "Pairing x CNO interaction check(s) failed.\n")
  quit(save = "no", status = 1L)
}
cat("All Pairing x CNO interaction checks passed.\n")
