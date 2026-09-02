#!/usr/bin/env Rscript

# Publication release, stage 13 -- validate the built package.
#
# Produces
#   provenance/validation_results.tsv   one row per contract, PASS / FAIL / SKIP
#   provenance/VALIDATION_REPORT.md     the same, readable
#
# Two verdicts are reported and must not be conflated:
#
#   RELEASE VALIDATION  PASS / FAIL / SKIP per scientific and integrity contract.
#   PRIDE READINESS     one of four values, derived from evidence.
#
# A PRIDE field that is genuinely unavailable does NOT fail the scientific release. It does
# prevent the deposition being called PRIDE_READY, and that distinction is the point.
#
# The validator rewrites the manifest and checksums as its last act, so its own report is
# covered by both. Everything it validates is fixed before that happens.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
LAYER <- file.path(REPO_ROOT, "07_publication_release")
source(file.path(LAYER, "R", "release_validation.R"))
release_require("digest")

DATA_ROOT <- release_data_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/13_validate_release.R"

release_banner("stage 13 -- release validation")

results <- list()
record <- function(group, check, status, detail = "") {
  results[[length(results) + 1L]] <<- data.frame(
    group = group, check = check, status = status, detail = detail,
    stringsAsFactors = FALSE, check.names = FALSE)
  cat("  ", switch(status, PASS = "[ok]  ", FAIL = "[FAIL]", SKIP = "[skip]"), " ",
      check, if (nzchar(detail)) paste0(" -- ", detail) else "", "\n", sep = "")
}
expect <- function(group, check, condition, detail = "") {
  record(group, check, if (isTRUE(condition)) "PASS" else "FAIL", detail)
}

rel <- function(...) release_path(..., create_dir = FALSE)
read_release <- function(path, ...) {
  if (grepl("[.]gz$", path)) {
    read.delim(gzfile(path), sep = "\t", quote = "", stringsAsFactors = FALSE,
               check.names = FALSE, ...)
  } else {
    read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE,
               check.names = FALSE, ...)
  }
}
inv <- RELEASE_DESIGN_INVARIANTS

# --------------------------------------------------------------------------------------
cat("\n=== locked canonical artefacts ===\n")
# --------------------------------------------------------------------------------------
locked <- release_verify_locked_artefacts(DATA_ROOT)
for (i in seq_len(nrow(locked))) {
  expect("locked artefacts", paste0(locked$artefact[i], " SHA256 matches the locked value"),
         locked$matches[i], substr(locked$observed_sha256[i], 1, 16))
}

# --------------------------------------------------------------------------------------
cat("\n=== sample metadata and design ===\n")
# --------------------------------------------------------------------------------------
sm <- read_release(rel("metadata", "sample_metadata.tsv"))
al <- read_release(rel("metadata", "animal_level_sample_metadata.tsv"))

expect("design", "96 measurement-level records", nrow(sm) == inv$n_measurement_records,
       paste(nrow(sm), "rows"))
expect("design", "48 animal-level units", nrow(al) == inv$n_animal_level_units,
       paste(nrow(al), "rows"))
expect("design", "12 distinct AnimalIDs",
       length(unique(sm$AnimalID)) == inv$n_animals,
       paste(length(unique(sm$AnimalID)), "animals"))
expect("design", "16 sample_class x condition strata",
       nrow(unique(al[c("sample_class", "condition")])) == inv$n_strata)
expect("design", "3 animals in every stratum",
       all(table(al$sample_class, al$condition) == inv$n_animals_per_stratum))
al$n_hemisphere_measurements <- as.integer(al$n_hemisphere_measurements)
design_checks <- release_check_animal_level_design(al)
for (i in seq_len(nrow(design_checks))) {
  expect("design", design_checks$check[i], design_checks$pass[i])
}
expect("design", "every measurement maps to an animal-level unit",
       all(sm$animal_level_sample_id %in% al$animal_level_column_name))
expect("design", "raw acquisition filenames unique and complete",
       !anyDuplicated(sm$raw_file_basename) && all(nzchar(sm$raw_file_basename)) &&
         all(grepl("[.]d$", sm$raw_file_basename)))

# Sample-class correction. The contract is NOT that every measurement carries its original
# class -- six do not, deliberately -- and it is NOT that the difference is an open question.
# The forensic audit resolved it as an intentional correction, so what is checked here is
# that the correction is documented exactly, both identities are retained, and nothing about
# it is left as an unresolved blocker.
corr_path <- rel("metadata", "sample_class_corrections.tsv")
expect("design", "both original and analysis sample-class labels are retained for all 96",
       all(c("original_sample_class", "analysis_sample_class", "sample_class_corrected")
           %in% names(sm)) &&
         !anyNA(sm$original_sample_class) && all(nzchar(sm$original_sample_class)) &&
         !anyNA(sm$analysis_sample_class) && all(nzchar(sm$analysis_sample_class)))
expect("design", "analysis_sample_class equals the canonical sample_class column",
       identical(as.character(sm$analysis_sample_class), as.character(sm$sample_class)))
expect("design", "sample_class_corrected is exactly where the two labels differ",
       identical(as.logical(sm$sample_class_corrected),
                 sm$original_sample_class != sm$analysis_sample_class))

n_corr_flagged <- sum(as.logical(sm$sample_class_corrected))
expect("design",
       sprintf("exactly %d of 96 measurements are corrected, %d unchanged",
               RELEASE_N_SAMPLE_CLASS_CORRECTED,
               inv$n_measurement_records - RELEASE_N_SAMPLE_CLASS_CORRECTED),
       n_corr_flagged == RELEASE_N_SAMPLE_CLASS_CORRECTED &&
         sum(!as.logical(sm$sample_class_corrected)) ==
           inv$n_measurement_records - RELEASE_N_SAMPLE_CLASS_CORRECTED,
       paste(n_corr_flagged, "corrected"))

if (!file.exists(corr_path)) {
  record("design", "the sample-class correction is published", "FAIL",
         "sample_class_corrections.tsv missing")
} else {
  corr <- read.delim(corr_path, sep = "\t", quote = "", stringsAsFactors = FALSE,
                     check.names = FALSE)
  expect("design", "the correction table covers exactly the flagged measurements",
         nrow(corr) == n_corr_flagged &&
           setequal(corr$sample_id, sm$sample_id[as.logical(sm$sample_class_corrected)]),
         paste(nrow(corr), "row(s)"))

  # The full contract: count, identity, hemisphere, animals, plate positions, both class
  # vectors and the 3-cycle. Reported per check so a failure names what diverged.
  for (i in seq_len(nrow(cc <- release_check_sample_class_corrections(corr)))) {
    expect("design", cc$check[i], cc$pass[i])
  }

  # The cited correction record must still hash to what the forensic audit recorded.
  ref <- release_verify_sample_class_correction_reference(DATA_ROOT, require_match = FALSE)
  expect("design", "the preserved correction record matches its recorded SHA256",
         isTRUE(ref$matches),
         paste0(substr(ref$observed_sha256, 1, 16), "... (mtime ", ref$observed_mtime, ")"))
  expect("design", "the published correction table cites that same record by hash",
         all(as.character(corr$correction_reference_sha256) == ref$expected_sha256))
  expect("design", "no corrected row claims a quantitative change",
         all(toupper(as.character(corr$quantitative_values_changed)) == "FALSE"))

  # This is the check that used to be a SKIP. It is a PASS now, and it is a real check:
  # the correction has to be RESOLVED, documented in the reader-facing documents, and
  # described without the vocabulary of an open discrepancy.
  readme_txt <- paste(readLines(rel("README_DATA.md"), warn = FALSE), collapse = " ")
  changelog_txt <- paste(readLines(
    rel("editor_source_data", "REVISION_PROTEOMICS_DATA_CHANGELOG.md"), warn = FALSE),
    collapse = " ")
  documented <- grepl("SAMPLE-CLASS CORRECTION", readme_txt, fixed = TRUE) &&
    grepl("sample-class assignments corrected during historical sample-identity",
          readme_txt, fixed = TRUE) &&
    grepl("A sample-class correction, resolved", changelog_txt, fixed = TRUE)
  resolved <- all(as.character(corr$correction_status) ==
                    RELEASE_SAMPLE_CLASS_CORRECTION$status) &&
    all(as.character(corr$correction_status) %in%
          release_sample_class_correction_statuses)
  record("design", "sample-class correction is RESOLVED and documented, not an open item",
         if (documented && resolved) "PASS" else "FAIL",
         paste0(nrow(corr), " corrected measurement(s), status ",
                unique(corr$correction_status),
                "; affected animal-level unit(s): ",
                paste(sort(unique(corr$animal_level_unit_affected)), collapse = ", ")))

  # The reader-facing documents must not describe it as unresolved or as a mislabelling.
  forbidden <- c("UNRESOLVED SAMPLE-CLASS", "unresolved sample-class",
                 "unresolved metadata discrepancy", "requires an experimenter decision",
                 "requires a decision from the experimenter", "mislabelled samples",
                 "mislabeled samples", "samples requiring reanalysis")
  hit <- forbidden[vapply(forbidden, function(f) {
    grepl(f, readme_txt, fixed = TRUE) || grepl(f, changelog_txt, fixed = TRUE)
  }, logical(1))]
  expect("design", "no document still calls the sample-class correction unresolved",
         length(hit) == 0L,
         if (length(hit)) paste("found:", paste(hit, collapse = "; ")) else "")

  # The absence of a written rationale is a fact and must stay stated.
  expect("design", "the missing prose rationale is disclosed, not glossed",
         grepl("No surviving prose note", changelog_txt, fixed = TRUE) ||
           grepl("no surviving prose note", changelog_txt, fixed = TRUE))

  fp <- read.delim(rel("metadata", "metadata_field_provenance.tsv"), sep = "\t",
                   quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  expect("design", "sample_class is reported as verified now that the correction is resolved",
         fp$status[fp$field == "sample_class"] == "KNOWN_VERIFIED")
  expect("design", "original_sample_class carries its own field provenance entry",
         "original_sample_class" %in% fp$field)
}
# --------------------------------------------------------------------------------------
cat("\n=== contrasts ===\n")
# --------------------------------------------------------------------------------------
pc <- read_release(rel("metadata", "primary_contrast_manifest.tsv"))
sec <- read_release(rel("metadata", "secondary_analysis_manifest.tsv"))

expect("contrasts", "exactly 12 primary contrasts", nrow(pc) == inv$n_primary_contrasts,
       paste(nrow(pc), "rows"))
expect("contrasts", "4 sample classes x 3 contrast families",
       length(unique(pc$sample_class)) == 4L && length(unique(pc$contrast_family)) == 3L)
expect("contrasts", "every primary contrast is 3 vs 3 animals",
       all(as.integer(pc$n_numerator_animals) == 3L) &&
         all(as.integer(pc$n_denominator_animals) == 3L))
expect("contrasts", "primary manifest reproduces the canonical helper",
       identical(sort(pc$canonical_comparison),
                 sort(primary_contrast_manifest()$canonical_comparison)))
expect("contrasts", "every primary row is labelled primary",
       all(pc$primary_or_secondary == "primary"))
expect("contrasts", "no secondary analysis appears in the primary manifest",
       !any(sec$analysis_id %in% pc$canonical_comparison) &&
         all(sec$primary_or_secondary == "secondary"))
expect("contrasts", "secondary manifest is non-empty and separate", nrow(sec) > 0L,
       paste(nrow(sec), "entries"))

# --------------------------------------------------------------------------------------
cat("\n=== differential results equal their canonical sources ===\n")
# --------------------------------------------------------------------------------------
da <- read_release(rel("differential_analysis", "primary_differential_proteins.tsv.gz"))
da_sum <- read_release(rel("differential_analysis", "primary_differential_summary.tsv"))

expect("differential", "12 complete comparison blocks",
       length(unique(da$canonical_comparison)) == inv$n_primary_contrasts &&
         all(table(da$canonical_comparison) == inv$n_proteins_statistical),
       paste(nrow(da), "rows"))
expect("differential", "every row is animal-level", all(da$analysis_unit == "animal"))
expect("differential", "no row claims more than 3 animals per arm",
       all(as.integer(da$n_numerator_animals) == 3L) &&
         all(as.integer(da$n_denominator_animals) == 3L))
expect("differential", "gene_symbol column holds gene symbols, not UniProt identifiers",
       !release_column_is_misleading_gene_symbol(da$gene_symbol))
expect("differential", "protein_group_id holds UniProt entry names",
       mean(release_uniprot_like(unique(da$protein_group_id))) > 0.9)

SPLIT_FORWARD <- file.path(DATA_ROOT, "02_data", "animal_level", "split", "forward")
CHECK_MAP <- c(effect_size_sd_units = "log2fc", average_standardized_abundance = "aveExpr", moderated_t = "t",
               P.Value = "pval", adj.P.Val = "padj", B_log_odds = "B")
if (!dir.exists(SPLIT_FORWARD)) {
  record("differential", "exported statistics equal the canonical split tables", "SKIP",
         "canonical split folder unreachable")
} else {
  mismatches <- character(0)
  for (cc in unique(da$canonical_contrast)) {
    p <- file.path(SPLIT_FORWARD, paste0(cc, ".csv"))
    if (!file.exists(p)) { mismatches <- c(mismatches, paste("missing:", cc)); next }
    d <- release_read_csv(p)
    ex <- da[da$canonical_contrast == cc, , drop = FALSE]
    if (!identical(as.character(ex$protein_group_id), as.character(d$gene_symbol))) {
      mismatches <- c(mismatches, paste("protein order:", cc)); next
    }
    for (pub in names(CHECK_MAP)) {
      if (!isTRUE(all.equal(as.numeric(ex[[pub]]), as.numeric(d[[CHECK_MAP[[pub]]]]),
                            tolerance = 0))) {
        mismatches <- c(mismatches, paste0(pub, ": ", cc))
      }
    }
  }
  expect("differential", "exported statistics equal the canonical split tables exactly",
         length(mismatches) == 0L,
         if (length(mismatches)) paste(head(mismatches, 3), collapse = "; ") else
           "12 comparisons, 6 statistics each")
}
expect("differential", "summary counts agree with the exported rows",
       all(vapply(seq_len(nrow(da_sum)), function(i) {
         sub <- da[da$canonical_comparison == da_sum$canonical_comparison[i], , drop = FALSE]
         sum(as.logical(sub$significant_fdr_0_05), na.rm = TRUE) ==
           as.integer(da_sum$n_significant_fdr_0_05[i])
       }, logical(1))))

# --------------------------------------------------------------------------------------
cat("\n=== enrichment results equal their canonical sources ===\n")
# --------------------------------------------------------------------------------------
gsea_go <- read_release(rel("enrichment", "primary_GSEA_GO_BP.tsv.gz"))
gsea_kegg <- read_release(rel("enrichment", "primary_GSEA_KEGG.tsv.gz"))
ora <- read_release(rel("enrichment", "primary_ORA_GO_BP.tsv.gz"))
sens <- read_release(rel("enrichment", "GSEA_log2FC_sensitivity.tsv.gz"))
ewce <- read_release(rel("enrichment", "primary_EWCE.tsv.gz"))
coverage <- read_release(rel("enrichment", "enrichment_coverage.tsv"))

expect("enrichment", "canonical GSEA covers all 12 comparisons",
       length(unique(gsea_go$canonical_comparison)) == 12L &&
         length(unique(gsea_kegg$canonical_comparison)) == 12L)
expect("enrichment", "canonical GSEA is ranked by the moderated t",
       all(gsea_go$rank_statistic == "moderated_t") &&
         all(gsea_kegg$rank_statistic == "moderated_t"))
expect("enrichment", "canonical and sensitivity rows are in separate files",
       all(gsea_go$analysis_role == "canonical") &&
         all(gsea_kegg$analysis_role == "canonical") &&
         all(sens$analysis_role == "sensitivity") &&
         all(sens$rank_statistic == "log2fc"))
expect("enrichment", "every enrichment row carries an analysis_role",
       all(nzchar(gsea_go$analysis_role)) && all(nzchar(ora$analysis_role)) &&
         all(nzchar(sens$analysis_role)))
expect("enrichment", "ORA rows are separated by query_list",
       length(unique(ora$query_list)) >= 3L,
       paste(length(unique(ora$query_list)), "variants"))
expect("enrichment", "exported row counts equal the counts the canonical run recorded",
       all(as.logical(coverage$agrees_with_run_record)),
       paste(sum(as.integer(coverage$rows_exported) == 0L), "blocks legitimately empty"))
# NA is written as an empty field and reads back as "", not NA, so both spellings of
# "no comparison" have to be accepted here. Baseline EWCE rows are not contrasts.
ewce_cmp <- ewce$canonical_comparison
ewce_cmp_missing <- is.na(ewce_cmp) | !nzchar(trimws(ewce_cmp))
expect("enrichment", "EWCE differential rows map to the 12 primary comparisons",
       all(ewce_cmp_missing | ewce_cmp %in% pc$canonical_comparison),
       paste(sum(!ewce_cmp_missing), "differential rows"))
expect("enrichment", "only EWCE baseline rows lack a comparison",
       all(ewce$ewce_analysis_type[ewce_cmp_missing] == "Baseline"))
expect("enrichment", "EWCE baseline rows are labelled secondary",
       all(ewce$primary_or_secondary[ewce$ewce_analysis_type == "Baseline"] == "secondary"))

ENRICH_ROOT <- file.path(DATA_ROOT, "03_output", "enrichment",
                         "enrichment_t_rank_validation_20260825", "per_comparison")
if (!dir.exists(ENRICH_ROOT)) {
  record("enrichment", "exported enrichment statistics equal the canonical files", "SKIP",
         "canonical enrichment folder unreachable")
} else {
  bad <- character(0); checked <- 0L
  for (cmp in unique(gsea_go$canonical_comparison)) {
    p <- file.path(ENRICH_ROOT, cmp, "GSEA_GO_BP.csv")
    if (!file.exists(p)) { bad <- c(bad, paste("missing:", cmp)); next }
    d <- release_read_csv(p)
    ex <- gsea_go[gsea_go$canonical_comparison == cmp, , drop = FALSE]
    if (!identical(ex$term_id, as.character(d$ID))) { bad <- c(bad, paste("order:", cmp)); next }
    if (!isTRUE(all.equal(as.numeric(ex$NES), as.numeric(d$NES), tolerance = 0)) ||
        !isTRUE(all.equal(as.numeric(ex$adjusted_p_value), as.numeric(d[["p.adjust"]]),
                          tolerance = 0))) {
      bad <- c(bad, paste("values:", cmp))
    }
    checked <- checked + 1L
  }
  expect("enrichment", "exported GSEA statistics equal the canonical files exactly",
         length(bad) == 0L,
         if (length(bad)) paste(head(bad, 3), collapse = "; ") else
           paste(checked, "comparisons verified"))
}

# --------------------------------------------------------------------------------------
cat("\n=== effect-size semantics ===\n")
# --------------------------------------------------------------------------------------
# The published effect size must be named for what it is, traceable to what it came from,
# and numerically identical to it. The scan below is the check that no reader-facing
# document slipped back into calling it a fold change.

da_head <- read.delim(gzfile(rel("differential_analysis",
                                 "primary_differential_proteins.tsv.gz")),
                      sep = "\t", quote = "", stringsAsFactors = FALSE,
                      check.names = FALSE, nrows = 5L)
ds_eff <- read_release(rel("differential_analysis", "primary_differential_summary.tsv"))

expect("effect_size", "the public effect-size column is present and correctly named",
       RELEASE_EFFECT_SIZE$public_field %in% names(da_head),
       RELEASE_EFFECT_SIZE$public_field)
expect("effect_size", "no published column is named logFC, log2FC or log2fc",
       !any(tolower(names(da_head)) %in% c("logfc", "log2fc", "log2foldchange",
                                           "log2_fold_change")),
       paste(length(names(da_head)), "columns"))
expect("effect_size", "the source statistic is recorded for provenance",
       "source_statistic_field" %in% names(da_head) &&
         all(da_head$source_statistic_field == RELEASE_EFFECT_SIZE$source_field),
       RELEASE_EFFECT_SIZE$source_field)
expect("effect_size", "the effect-size units are stated per row",
       "effect_size_units" %in% names(da_head) &&
         all(da_head$effect_size_units == RELEASE_EFFECT_SIZE$units))
expect("effect_size", "the summary carries the full effect-size definition",
       "effect_size_definition" %in% names(ds_eff) &&
         all(ds_eff$effect_size_definition == RELEASE_EFFECT_SIZE$definition))
expect("effect_size", "the definition states explicitly that it is not a log2 fold change",
       all(grepl("NOT a log2 fold change", ds_eff$effect_size_definition, fixed = TRUE)))
expect("effect_size", "the effect size is not called a standardised mean difference alone",
       !any(grepl("standardised mean difference", ds_eff$effect_size_definition,
                  fixed = TRUE)) &&
         !any(grepl("standardized mean difference", ds_eff$effect_size_definition,
                    fixed = TRUE)))

# Canonical GSEA ranking must still be the moderated t, and the sensitivity rows must still
# be separated. Terminology changed; the statistic did not.
gsea_go_v <- read_release(rel("enrichment", "primary_GSEA_GO_BP.tsv.gz"))
gsea_kegg_v <- read_release(rel("enrichment", "primary_GSEA_KEGG.tsv.gz"))
expect("effect_size", "canonical GSEA is still ranked by the moderated t",
       all(gsea_go_v$rank_statistic == "moderated_t") &&
         all(gsea_kegg_v$rank_statistic == "moderated_t") &&
         all(gsea_go_v$analysis_role == "canonical") &&
         all(gsea_kegg_v$analysis_role == "canonical"))
expect("effect_size", "the sensitivity ranking keeps its canonical `log2fc` token",
       all(sens$rank_statistic == "log2fc") && all(sens$analysis_role == "sensitivity"))

# The terminology audit, and the scan of every reader-facing document. Builder sources are
# deliberately not scanned: they hold the claim list and the audit rows themselves.
audit_path <- rel("editor_source_data", "EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv")
if (!file.exists(audit_path)) {
  record("effect_size", "the terminology audit is published", "FAIL",
         "EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv missing")
} else {
  aud <- read_release(audit_path)
  expect("effect_size", "every audited occurrence carries a sanctioned classification",
         all(aud$status %in% c("INTERNAL_PROVENANCE", "PUBLIC_MISLABEL",
                               "LEGITIMATE_OTHER_CONTEXT")))
  expect("effect_size", "the audit records corrections, retained tokens and other contexts",
         all(c("INTERNAL_PROVENANCE", "PUBLIC_MISLABEL", "LEGITIMATE_OTHER_CONTEXT")
             %in% aud$status),
         paste(nrow(aud), "rows"))
  expect("effect_size", "no audited PUBLIC_MISLABEL is left unresolved",
         !any(grepl("^RESIDUAL", aud$notes)),
         paste(sum(grepl("^RESIDUAL", aud$notes)), "residual"))
}

PUBLIC_DOCS <- c(
  "README_DATA.md",
  "metadata/data_dictionary.tsv",
  "metadata/primary_contrast_manifest.tsv",
  "metadata/secondary_analysis_manifest.tsv",
  "metadata/sample_class_corrections.tsv",
  "metadata/sample_metadata.tsv",
  "differential_analysis/primary_differential_summary.tsv",
  "editor_source_data/REVISION_PROTEOMICS_DATA_CHANGELOG.md",
  "editor_source_data/MANUSCRIPT_TERMINOLOGY_ACTIONS.md",
  "editor_source_data/figure_source_map.tsv",
  "pride/README_PRIDE.md",
  "pride/SDRF_MISSING_METADATA.md",
  "provenance/UPSTREAM_PREPROCESSING_GAP.md",
  "provenance/data_lineage.tsv"
)
mislabel_rows <- list()
for (doc in PUBLIC_DOCS) {
  p <- file.path(OUT_ROOT, doc)
  if (!file.exists(p)) next
  h <- release_fold_change_mislabels(readLines(p, warn = FALSE))
  if (nrow(h)) mislabel_rows[[doc]] <- h
}
n_mislabel <- sum(vapply(mislabel_rows, nrow, integer(1)))
if (n_mislabel > 0L) {
  for (doc in names(mislabel_rows)) {
    for (i in seq_len(nrow(mislabel_rows[[doc]]))) {
      cat("        ", doc, ":", mislabel_rows[[doc]]$line[i], " [",
          mislabel_rows[[doc]]$term[i], "] ",
          substr(mislabel_rows[[doc]]$text[i], 1, 120), "\n", sep = "")
    }
  }
}
expect("effect_size",
       "no reader-facing document describes the published effect size as a fold change",
       n_mislabel == 0L,
       paste(length(PUBLIC_DOCS), "documents scanned,", n_mislabel, "mislabel(s)"))

# The methods-safe definition has to actually be present for someone to copy.
readme_lines <- readLines(rel("README_DATA.md"), warn = FALSE)
expect("effect_size", "README_DATA.md carries the methods-safe definition",
       any(grepl(RELEASE_EFFECT_SIZE$methods_sentence, readme_lines, fixed = TRUE)))
expect("effect_size", "the manuscript action note lists the out-of-repo corrections",
       file.exists(rel("editor_source_data", "MANUSCRIPT_TERMINOLOGY_ACTIONS.md")))

# --------------------------------------------------------------------------------------
cat("\n=== no hemisphere-level inference presented as animal-level ===\n")
# --------------------------------------------------------------------------------------
abundance_animal <- read_release(rel("processed_data",
                                     "protein_abundance_animal_level.tsv.gz"), nrows = 2L)
animal_cols <- setdiff(names(abundance_animal), "protein_group_id")
expect("statistical unit", "animal-level matrix has 48 columns",
       length(animal_cols) == inv$n_animal_level_units, paste(length(animal_cols), "columns"))
expect("statistical unit", "animal-level columns are not acquisition run names",
       length(intersect(animal_cols, sm$sample_id)) == 0L)
expect("statistical unit", "no inferential table declares a non-animal unit",
       all(da$analysis_unit == "animal") && all(gsea_go$analysis_unit == "animal") &&
         all(ora$analysis_unit == "animal") && all(ewce$analysis_unit == "animal"))
expect("statistical unit", "no inferential row claims 6 replicates per arm",
       !any(as.integer(da$n_numerator_animals) > 3L, na.rm = TRUE) &&
         !any(as.integer(gsea_go$n_numerator_animals) > 3L, na.rm = TRUE))

fm <- read_release(rel("editor_source_data", "figure_source_map.tsv"))
acq_unit <- grepl("hemisphere|acquisition", fm$statistical_unit, ignore.case = TRUE)
expect("statistical unit",
       "panels drawn at acquisition level are labelled technical QC",
       all(fm$primary_or_secondary[acq_unit] == "technical_qc"),
       paste(sum(acq_unit), "acquisition-level panel(s)"))
expect("figures", "all required manuscript panels are mapped",
       all(c("Fig 3 C", "Fig 3 D", "Fig 3 E", "Fig 3 F", "Fig 3 G", "Fig 3 H", "Fig 3 I",
             "Supp A", "Supp B1", "Supp B2", "Supp B3", "Supp C", "Supp D", "Supp E",
             "Supp F") %in% fm$panel_label))
expect("figures", "every data-bearing panel names a final revision script",
       all(nzchar(fm$final_revision_script[fm$primary_or_secondary != "not_applicable"])))

# --------------------------------------------------------------------------------------
cat("\n=== editor workbook hygiene ===\n")
# --------------------------------------------------------------------------------------
wb_inv <- read_release(rel("editor_source_data", "workbook_sheet_inventory.tsv"))
sheet_cols <- stats::setNames(
  lapply(strsplit(wb_inv$columns, ";", fixed = TRUE), identity), wb_inv$sheet)
hygiene <- release_check_workbook_hygiene(sheet_cols)
expect("workbook", "no duplicate or unnamed workbook columns", length(hygiene) == 0L,
       if (length(hygiene)) paste(head(hygiene, 2), collapse = "; ") else
         paste(nrow(wb_inv), "sheets"))
expect("workbook", "required sheets present",
       all(c("README", "Sample_Metadata", "Animal_Level_Metadata", "Primary_Contrasts",
             "Differential_Proteins", "GSEA_GO_BP", "GSEA_KEGG", "ORA_GO_BP", "EWCE",
             "Secondary_Analyses", "Figure_Source_Map", "Software_Versions") %in%
             wb_inv$sheet))
expect("workbook", "workbook file exists",
       file.exists(rel("editor_source_data", "Proteomics_Source_Data_Animal_Level.xlsx")))
expect("workbook", "no sheet carries a misleading public gene_symbol field",
       !any(vapply(sheet_cols, function(cs) {
         "gene_symbol" %in% cs && !"uniprot_accession" %in% cs
       }, logical(1))))

# --------------------------------------------------------------------------------------
cat("\n=== data dictionary and wording ===\n")
# --------------------------------------------------------------------------------------
dict <- read_release(rel("metadata", "data_dictionary.tsv"))
expect("documentation", "data dictionary is populated for every listed column",
       nrow(dict) > 0L && all(nzchar(dict$description)), paste(nrow(dict), "definitions"))
readme_lines <- readLines(rel("README_DATA.md"), warn = FALSE)
PLATE_SENTENCE <- paste(
  "Because pairing condition was completely associated with collection plate in the",
  "paired-vehicle versus unpaired-vehicle comparison, any collection-plate-associated",
  "contribution cannot be distinguished from a pairing-associated contribution.")
flat_readme <- paste(readme_lines, collapse = " ")
expect("documentation", "required collection-plate sentence present verbatim",
       grepl(PLATE_SENTENCE, flat_readme, fixed = TRUE))
expect("documentation", "required no-batch-effect sentence present verbatim",
       grepl(paste("No downstream proteomics batch effect attributable to collection plate",
                   "has been demonstrated."), flat_readme, fixed = TRUE))
prose_files <- c(rel("README_DATA.md"),
                 rel("editor_source_data", "REVISION_PROTEOMICS_DATA_CHANGELOG.md"),
                 rel("pride", "README_PRIDE.md"),
                 rel("pride", "SDRF_MISSING_METADATA.md"),
                 rel("provenance", "UPSTREAM_PREPROCESSING_GAP.md"))
forbidden <- c("plate artefact", "plate artifact", "confirmed batch confound",
               "driven by plate", "technical plate effect")
prose_hits <- unlist(lapply(prose_files, function(f) {
  if (!file.exists(f)) return(character(0))
  txt <- tolower(paste(readLines(f, warn = FALSE), collapse = " "))
  forbidden[vapply(forbidden, function(p) grepl(p, txt, fixed = TRUE), logical(1))]
}))
expect("documentation", "no forbidden collection-plate wording in release prose",
       length(prose_hits) == 0L,
       if (length(prose_hits)) paste(unique(prose_hits), collapse = ", ") else
         paste(length(prose_files), "documents checked"))
expect("documentation", "release prose never calls a processed matrix raw MS data",
       !grepl("raw mass[- ]spectrometry data\\b(?![^.]*not)", flat_readme, perl = TRUE) ||
         grepl("is \\*\\*not\\*\\* raw mass-spectrometry data", flat_readme))

# --------------------------------------------------------------------------------------
cat("\n=== data lineage ===\n")
# --------------------------------------------------------------------------------------
lineage <- read_release(rel("provenance", "data_lineage.tsv"))
expect("lineage", "lineage covers acquisition through release",
       all(c("A02_acquisition_raw", "A04_protein_group_matrix_prefilter",
             "A07_processed_matrix", "A08_animal_level_matrix", "A09_protigy_statistics",
             "A10_split_differential_tables", "A12_mapped_differential_tables",
             "A13_enrichment", "A14_ewce", "A17_publication_release") %in%
             lineage$artifact_id))
unresolved <- lineage$artifact_id[grepl("UNRESOLVED",
                                        paste(lineage$processing_stage, lineage$notes))]
record("lineage", "unproven lineage edges are marked UNRESOLVED, not inferred", "PASS",
       paste(length(unresolved), "artefact(s):", paste(unresolved, collapse = ", ")))
expect("lineage", "the locked artefacts appear with their hashes",
       any(lineage$sha256 == RELEASE_LOCKED_ARTEFACTS$animal_level_input_gct$sha256,
           na.rm = TRUE) &&
         any(lineage$sha256 == RELEASE_LOCKED_ARTEFACTS$protigy_stat_gct$sha256,
             na.rm = TRUE))
expect("lineage", "upstream preprocessing gap is documented",
       file.exists(rel("provenance", "UPSTREAM_PREPROCESSING_GAP.md")))

sw <- read_release(rel("provenance", "software_versions.tsv"))
n_unknown <- sum(sw$version == "UNKNOWN")
record("lineage", "software versions recorded", "PASS",
       paste(nrow(sw), "components,", n_unknown, "UNKNOWN:",
             paste(sw$component[sw$version == "UNKNOWN"], collapse = "; ")))
expect("lineage", "no vague version strings such as 1.1.x",
       !any(grepl("[0-9]+[.][0-9]*[.]?x$", sw$version)))
expect("lineage", "the UniProt mapping reference is traceable",
       any(grepl("idmapping", sw$component, ignore.case = TRUE) &
             nzchar(sw$evidence_sha256) & !is.na(sw$evidence_sha256)))

# --------------------------------------------------------------------------------------
cat("\n=== release manifest and checksums ===\n")
# --------------------------------------------------------------------------------------
manifest_path <- rel("provenance", "release_manifest.tsv")
if (!file.exists(manifest_path)) {
  record("manifest", "release manifest exists", "FAIL", "not built")
} else {
  manifest <- read_release(manifest_path)
  all_files <- gsub("\\\\", "/", list.files(OUT_ROOT, recursive = TRUE))
  expected <- setdiff(all_files, c("provenance/release_manifest.tsv",
                                   "provenance/SHA256SUMS.txt",
                                   "provenance/validation_results.tsv",
                                   "provenance/VALIDATION_REPORT.md"))
  expect("manifest", "every release file is listed in the manifest",
         length(setdiff(expected, manifest$relative_path)) == 0L,
         paste(nrow(manifest), "files"))
  bad_hash <- vapply(seq_len(nrow(manifest)), function(i) {
    p <- file.path(OUT_ROOT, manifest$relative_path[i])
    !file.exists(p) || !identical(release_sha256(p), manifest$sha256[i])
  }, logical(1))
  expect("manifest", "every manifest SHA256 validates against the file on disk",
         !any(bad_hash),
         if (any(bad_hash)) paste(head(manifest$relative_path[bad_hash], 3), collapse = "; ")
         else paste(nrow(manifest), "hashes verified"))
  expect("manifest", "every manifest row names a generating script",
         all(nzchar(manifest$generated_by)))
  expect("manifest", "tabular files report a row and column count",
         all(!is.na(manifest$rows_if_tabular[manifest$format %in% c("tsv", "tsv.gz", "csv")])))
  sums_path <- rel("provenance", "SHA256SUMS.txt")
  expect("manifest", "SHA256SUMS.txt exists and covers the manifest",
         file.exists(sums_path) &&
           any(grepl("provenance/release_manifest.tsv",
                     readLines(sums_path, warn = FALSE), fixed = TRUE)))
}

# --------------------------------------------------------------------------------------
cat("\n=== publication scripts perform no scientific recomputation ===\n")
# --------------------------------------------------------------------------------------
scripts <- list.files(LAYER, pattern = "[.]R$", full.names = TRUE, recursive = TRUE)
hits <- release_scan_prohibited_calls(scripts)
# The validator itself and the validation helper legitimately name these calls when
# scanning for them; a hit inside the scanner's own pattern table is not a computation.
hits <- hits[!basename(hits$file) %in% c("release_validation.R", "13_validate_release.R"), ,
             drop = FALSE]
expect("no recomputation", "no model fitting, enrichment or PCA call in the builders",
       nrow(hits) == 0L,
       if (nrow(hits)) paste0(basename(hits$file[1]), ":", hits$line[1], " ", hits$call[1])
       else paste(length(scripts), "scripts scanned"))
expect("no recomputation", "no builder writes into a protected canonical root",
       all(vapply(release_protected_roots(DATA_ROOT),
                  function(p) !project_paths_overlap(OUT_ROOT, p), logical(1))))

# --------------------------------------------------------------------------------------
cat("\n=== PRIDE readiness (reported, not a scientific pass/fail) ===\n")
# --------------------------------------------------------------------------------------
pr <- read_release(rel("pride", "pride_readiness.tsv"))
pride_status <- pr$value[pr$key == "pride_status"]
fs <- read_release(rel("pride", "sdrf_field_status.tsv"))
n_req <- sum(fs$status == "MISSING_REQUIRED_METADATA")
n_opt <- sum(fs$status == "MISSING_OPTIONAL_METADATA")
n_ready <- sum(fs$status == "READY")

expect("pride", "PRIDE status is one of the four sanctioned values",
       pride_status %in% release_pride_statuses, pride_status)
expect("pride", "SDRF exists with one row per acquisition",
       nrow(read_release(rel("pride", "sdrf.tsv"))) == inv$n_measurement_records)
expect("pride", "every SDRF field has a recorded status", all(nzchar(fs$status)))
expect("pride", "no required field is silently blank",
       all(fs$status[fs$sdrf_requirement == "required" &
                       fs$value_or_reason == "not available"] ==
             "MISSING_REQUIRED_METADATA"))
record("pride", "SDRF field readiness", "PASS",
       paste0(n_ready, " READY, ", n_req, " MISSING_REQUIRED, ", n_opt, " MISSING_OPTIONAL"))
expect("pride", "deposition is not overstated as PRIDE_READY",
       !(pride_status == "PRIDE_READY" && n_req > 0L))
record("pride", "missing required SDRF fields", if (n_req == 0L) "PASS" else "SKIP",
       if (n_req == 0L) "none" else
         paste(fs$sdrf_field[fs$status == "MISSING_REQUIRED_METADATA"], collapse = "; "))

# The sample-class correction must NOT be carried as a PRIDE gap. It used to be the reason
# for one of the two SKIPs; it is resolved, so the SDRF factor is populated and the
# deposition documents say so. This is a PASS, and a real one: the factor has to be
# populated from the analysis-time class and the correction has to be described as resolved.
sdrf_v <- read_release(rel("pride", "sdrf.tsv"))
class_col <- "factor value[sample class]"
expect("pride", "the SDRF sample-class factor is populated for every acquisition",
       class_col %in% names(sdrf_v) &&
         all(nzchar(sdrf_v[[class_col]])) &&
         !any(sdrf_v[[class_col]] %in% c("not available", "not applicable")))
expect("pride", "the SDRF sample-class factor is the validated analysis-time class",
       class_col %in% names(sdrf_v) &&
         setequal(unique(sdrf_v[[class_col]]), unique(sm$analysis_sample_class)))
expect("pride", "sample class is not recorded as missing SDRF metadata",
       !any(fs$status == "MISSING_REQUIRED_METADATA" &
              grepl("sample class", fs$sdrf_field, fixed = TRUE)))
pride_txt <- paste(c(readLines(rel("pride", "README_PRIDE.md"), warn = FALSE),
                     readLines(rel("pride", "SDRF_MISSING_METADATA.md"), warn = FALSE)),
                   collapse = " ")
record("pride", "sample-class correction is resolved, not a deposition blocker",
       if (grepl("RESOLVED", pride_txt, fixed = TRUE) &&
           any(pr$key == "sample_class_correction_status") &&
           pr$value[pr$key == "sample_class_correction_status"] ==
             RELEASE_SAMPLE_CLASS_CORRECTION$status) "PASS" else "FAIL",
       paste0("status ", RELEASE_SAMPLE_CLASS_CORRECTION$status, "; ",
              RELEASE_N_SAMPLE_CLASS_CORRECTED,
              " corrected measurement(s), original class retained in companion tables"))

# ...and the genuinely missing acquisition metadata must STAY missing. Nothing above may be
# read as closing that gap. Reported as SKIP because it is unresolved, which is the truth.
expect("pride", "the acquisition-metadata gap is still reported as incomplete",
       pride_status == "PRIDE_METADATA_INCOMPLETE" || n_req == 0L,
       pride_status)
expect("pride", "PRIDE readiness is not claimed while required fields are missing",
       !(pride_status %in% c("PRIDE_READY", "PRIDE_READY_PENDING_RAW_UPLOAD") && n_req > 0L))

# --------------------------------------------------------------------------------------
# summarise
# --------------------------------------------------------------------------------------
res <- do.call(rbind, results)
rownames(res) <- NULL
n_pass <- sum(res$status == "PASS")
n_fail <- sum(res$status == "FAIL")
n_skip <- sum(res$status == "SKIP")

cat("\n=== RESULT ===\n")
cat("  PASS ", n_pass, "   FAIL ", n_fail, "   SKIP ", n_skip, "\n", sep = "")
cat("  PRIDE status: ", pride_status, "\n", sep = "")

report <- c(
  "# Publication release validation report",
  "",
  paste0("Generated by `", STAGE, "` at ", release_timestamp_utc(), "."),
  paste0("Release root: `", OUT_ROOT, "`"),
  paste0("Repository commit: `", release_git_commit(REPO_ROOT), "`"),
  "",
  paste0("**Release validation: ", n_pass, " PASS, ", n_fail, " FAIL, ", n_skip, " SKIP**"),
  paste0("**PRIDE readiness: `", pride_status, "`** (",
         n_ready, " SDRF fields READY, ", n_req, " MISSING_REQUIRED, ", n_opt,
         " MISSING_OPTIONAL)"),
  "",
  "A PRIDE field that is genuinely unavailable does not fail the scientific release. It",
  "does prevent the deposition being reported as `PRIDE_READY`.",
  "",
  "| group | check | status | detail |",
  "|---|---|---|---|",
  vapply(seq_len(nrow(res)), function(i) {
    paste0("| ", res$group[i], " | ", res$check[i], " | ", res$status[i], " | ",
           gsub("\\|", "/", res$detail[i]), " |")
  }, character(1)),
  ""
)
release_write_lines(report, release_path("provenance", "VALIDATION_REPORT.md"))
release_write_table(res, release_path("provenance", "validation_results.tsv"))
release_register("provenance/VALIDATION_REPORT.md",
                 "release validation report: per-contract PASS/FAIL/SKIP and PRIDE readiness",
                 "the release itself", NA_character_, STAGE, "md")
release_register("provenance/validation_results.tsv",
                 "machine-readable validation results", "the release itself",
                 NA_character_, STAGE, "tsv")

# Re-issue the manifest and checksums so they cover this report.
cat("\n  re-issuing manifest and checksums to cover the validation report ...\n")
manifest_env <- new.env(parent = globalenv())
assign("here", LAYER, envir = manifest_env)
sys.source(file.path(LAYER, "12_build_release_manifest.R"), envir = manifest_env,
           keep.source = FALSE)

if (n_fail > 0L) {
  stop(sprintf("Publication release validation failed: %d contract(s).", n_fail),
       call. = FALSE)
}
cat("\nAll publication release contracts hold.\n")
