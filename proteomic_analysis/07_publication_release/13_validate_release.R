#!/usr/bin/env Rscript

# Publication release, stage 13: validate the built package.
#
# Produces
#   provenance/validation_results.tsv   one row per check, PASS / FAIL / SKIP
#   provenance/VALIDATION_REPORT.md     the same as a readable report
#
# Two separate results: PASS / FAIL / SKIP per check, and a PRIDE readiness status (one of
# release_pride_statuses). Missing PRIDE metadata does not fail the release; it only keeps
# the status below PRIDE_READY.
#
# Run after run_release.R with the same PROTEOMICS_RELEASE_OUTPUT_ROOT. The last step
# rebuilds the manifest and checksums (stage 12) so they cover this report; the script then
# stops with an error if any check failed.

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
source(file.path(LAYER, "R", "release_software_versions.R"))
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
cat("\n=== locked canonical artifacts ===\n")
# --------------------------------------------------------------------------------------
locked <- release_verify_locked_artefacts(DATA_ROOT)
for (i in seq_len(nrow(locked))) {
  expect("locked artifacts", paste0(locked$artefact[i], " SHA256 matches the locked value"),
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
plate_layout <- tryCatch(release_plate_layout_sentence(sm), error = function(e) NA_character_)
expect("design", "collection plate and injection order are as the release states them",
       release_plate_injection_blocks_hold(sm) && !is.na(plate_layout),
       paste(release_plate_injection_sentence(), plate_layout))

experimenter_path <- rel("metadata", "experimenter_metadata.tsv")
protocol_path <- rel("metadata", "sample_preparation_protocol.tsv")
experimenter_metadata <- read_release(experimenter_path, na.strings = character(0),
                                      colClasses = "character")
sample_preparation_protocol <- read_release(protocol_path, na.strings = character(0),
                                            colClasses = "character")
source_experimenter <- release_read_experimenter_metadata(REPO_ROOT)
source_protocol <- release_read_sample_preparation_protocol(REPO_ROOT)
expect("metadata provenance", "released experimenter metadata equals the checked-in source",
       identical(experimenter_metadata, source_experimenter),
       paste(nrow(experimenter_metadata), "records"))
expect("metadata provenance", "released sample-preparation protocol equals the checked-in source",
       identical(sample_preparation_protocol, source_protocol),
       paste(nrow(sample_preparation_protocol), "ordered steps"))
expect("metadata provenance", "sample-preparation protocol is the 2024-11-28 version",
       identical(unique(sample_preparation_protocol$protocol_version_date), "2024-11-28"))

# Sample-class correction: six measurements were reassigned (see
# RELEASE_SAMPLE_CLASS_CORRECTION). Check that both labels are kept, the six rows match the
# record, and the documents describe the correction as resolved.
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

  # Count, identity, hemisphere, animals, plate positions, both class vectors and the cycle,
  # one check each so a failure names the field.
  for (i in seq_len(nrow(cc <- release_check_sample_class_corrections(corr)))) {
    expect("design", cc$check[i], cc$pass[i])
  }

  # The cited correction record must still hash to what the audit recorded.
  ref <- release_verify_sample_class_correction_reference(DATA_ROOT, require_match = FALSE)
  expect("design", "the preserved correction record matches its recorded SHA256",
         isTRUE(ref$matches),
         paste0(substr(ref$observed_sha256, 1, 16), "... (mtime ", ref$observed_mtime, ")"))
  expect("design", "the published correction table cites that same record by hash",
         all(as.character(corr$correction_reference_sha256) == ref$expected_sha256))
  expect("design", "no corrected row claims a quantitative change",
         all(toupper(as.character(corr$quantitative_values_changed)) == "FALSE"))

  # The correction must carry the resolved status and be described in README_DATA.md and
  # the changelog.
  readme_txt <- paste(readLines(rel("README_DATA.md"), warn = FALSE), collapse = " ")
  changelog_txt <- paste(readLines(
    rel("provenance", "REVISION_PROTEOMICS_DATA_CHANGELOG.md"), warn = FALSE),
    collapse = " ")
  documented <- grepl("Sample-class correction", readme_txt, fixed = TRUE) &&
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

  # The reader-facing documents must not describe it as unresolved or as a mislabeling.
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

  # The changelog must say that no written rationale survives.
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
expect("contrasts", "every primary row is labeled primary",
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

# The split and mapped tables the release reads must descend from the locked statistics GCT.
for (ix in list(list(label = "split index", path = file.path(release_stats_root(DATA_ROOT),
                                                             "split", "indexComparisons.csv")),
                list(label = "mapped index", path = file.path(release_stats_root(DATA_ROOT),
                                                              "mapped",
                                                              "indexMappedComparisons.csv")))) {
  if (!file.exists(ix$path)) {
    record("differential", paste("the", ix$label, "cites the locked statistics GCT"), "SKIP",
           "canonical index unreachable")
  } else {
    d <- release_read_csv(ix$path)
    expect("differential", paste("the", ix$label, "cites the locked statistics GCT"),
           nrow(d) == inv$n_primary_contrasts &&
             all(as.character(d$source_gct_sha256) == RELEASE_LOCKED_ARTEFACTS$stat_gct$sha256),
           paste(nrow(d), "comparisons"))
  }
}

SPLIT_FORWARD <- file.path(release_stats_root(DATA_ROOT), "split", "forward")
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

# Table 2d: the significant annotated proteins, copied from the protein table by its own
# significance flag and mapping status.
sig_p <- read_release(rel("differential_analysis", "significant_proteins.tsv"))
SIG_PROTEIN_COLUMNS <- c("canonical_comparison", "sample_class", "numerator_condition",
                         "denominator_condition", "gene_symbol", "protein_description",
                         "uniprot_accession", "protein_group_id", "direction",
                         "effect_size_sd_units", "moderated_t", "P.Value", "adj.P.Val")
sig_expected <- da[as.logical(da$significant_fdr_0_05) %in% TRUE &
                     da$id_mapping_status == "mapped", , drop = FALSE]
sig_m <- match(paste(sig_p$canonical_comparison, sig_p$protein_group_id),
               paste(sig_expected$canonical_comparison, sig_expected$protein_group_id))
expect("differential", "significant_proteins.tsv has the expected columns in order",
       identical(names(sig_p), SIG_PROTEIN_COLUMNS), paste(ncol(sig_p), "columns"))
expect("differential",
       "significant_proteins.tsv holds exactly the FDR < 0.05 annotated rows of the protein table",
       nrow(sig_p) == nrow(sig_expected) && !anyNA(sig_m) && !anyDuplicated(sig_m) &&
         nrow(sig_p) == sum(as.integer(da_sum$n_significant_fdr_0_05_mapped_only)),
       paste(nrow(sig_p), "rows;", sum(as.logical(da$significant_fdr_0_05) %in% TRUE &
                                          da$id_mapping_status != "mapped"),
             "unmapped FDR < 0.05 groups stay in the full table only"))
expect("differential", "significant_proteins.tsv values are identical to the protein table",
       !anyNA(sig_m) &&
         all(vapply(setdiff(SIG_PROTEIN_COLUMNS, c("canonical_comparison", "protein_group_id",
                                                   "direction")), function(cn) {
           identical(sig_p[[cn]], sig_expected[[cn]][sig_m])
         }, logical(1))))
expect("differential", "significant_proteins.tsv direction follows the sign of the effect size",
       identical(sig_p$direction,
                 paste0("higher_in_", ifelse(sig_p$effect_size_sd_units > 0,
                                             sig_p$numerator_condition,
                                             sig_p$denominator_condition))))
expect("differential", "significant_proteins.tsv is sorted by comparison, then adjusted P",
       !anyDuplicated(rle(sig_p$canonical_comparison)$values) &&
         all(vapply(split(sig_p$adj.P.Val, sig_p$canonical_comparison),
                    function(p) !is.unsorted(p), logical(1))))

# --------------------------------------------------------------------------------------
cat("\n=== the differential statistics are the per-sample-class model ===\n")
# --------------------------------------------------------------------------------------
# A release built from the superseded pooled run would pass every equality check above
# against its own tables, so the model itself is checked against its run manifest.
stat_run <- tryCatch(release_read_stat_run(DATA_ROOT), error = function(e) e)
if (inherits(stat_run, "error")) {
  record("statistics", "the statistics run manifest is readable and describes the primary model",
         "FAIL", conditionMessage(stat_run))
} else {
  sm_stat <- stat_run$manifest
  expect("statistics", "one manifest row per primary contrast, all the per-sample-class model",
         nrow(sm_stat) == inv$n_primary_contrasts &&
           setequal(sm_stat$canonical_comparison, pc$canonical_comparison) &&
           all(sm_stat$variance_model == RELEASE_PRIMARY_MODEL$variance_model))
  expect("statistics", "each model uses one sample class's 12 animal-level units",
         all(as.integer(sm_stat$n_columns_in_model) == inv$n_animals) &&
           all(as.integer(sm_stat$df_residual) == inv$n_animals - inv$n_conditions),
         paste0("residual df ", paste(unique(sm_stat$df_residual), collapse = ";")))
  expect("statistics", "the run manifest cites the locked input and statistics GCTs",
         identical(stat_run$input_gct_sha256,
                   RELEASE_LOCKED_ARTEFACTS$animal_level_input_gct$sha256) &&
           identical(stat_run$output_gct_sha256, RELEASE_LOCKED_ARTEFACTS$stat_gct$sha256))
  m_fdr <- stats::setNames(as.integer(sm_stat$n_fdr_proteins), sm_stat$canonical_comparison)
  expect("statistics", "released FDR counts equal the counts the statistics run recorded",
         identical(unname(m_fdr[da_sum$canonical_comparison]),
                   as.integer(da_sum$n_significant_fdr_0_05)),
         paste(sum(m_fdr), "FDR-significant protein rows in total"))
  expect("statistics", "every contrast row names the primary model, residual df and limma version",
         all(pc$statistical_model == RELEASE_PRIMARY_MODEL$label) &&
           all(as.integer(pc$residual_df) == stat_run$df_residual) &&
           all(grepl(paste0("limma ", stat_run$limma_version), pc$statistics_software,
                     fixed = TRUE)))
  ave <- release_aveexpr_summary(da)
  expect("statistics", "AveExpr is the sample-class mean; over the four classes it averages ~0",
         ave$n_proteins == inv$n_proteins_statistical && ave$max_abs_mean_over_classes < 0.0015,
         sprintf("max |mean over classes| %.5f; within-class range %.2f to %.2f",
                 ave$max_abs_mean_over_classes, ave$range_within_classes[1],
                 ave$range_within_classes[2]))
}

# --------------------------------------------------------------------------------------
cat("\n=== enrichment results equal their canonical sources ===\n")
# --------------------------------------------------------------------------------------
gsea_go <- read_release(rel("enrichment", "primary_GSEA_GO_BP.tsv.gz"))
gsea_cc <- read_release(rel("enrichment", "primary_GSEA_GO_CC.tsv.gz"))
gsea_kegg <- read_release(rel("enrichment", "primary_GSEA_KEGG.tsv.gz"))
ora <- read_release(rel("enrichment", "primary_ORA_GO_BP.tsv.gz"))
ora_cc <- read_release(rel("enrichment", "primary_ORA_GO_CC.tsv.gz"))
sens <- read_release(rel("enrichment", "GSEA_log2FC_sensitivity.tsv.gz"))
ewce <- read_release(rel("enrichment", "primary_EWCE.tsv.gz"))
coverage <- read_release(rel("enrichment", "enrichment_coverage.tsv"))
untested_rel <- read_release(rel("enrichment", "gsea_untested_gene_sets.tsv"))

expect("enrichment", "canonical GSEA covers all 12 comparisons",
       length(unique(gsea_go$canonical_comparison)) == 12L &&
         length(unique(gsea_cc$canonical_comparison)) == 12L &&
         length(unique(gsea_kegg$canonical_comparison)) == 12L)
expect("enrichment", "canonical GSEA is ranked by the moderated t",
       all(gsea_go$rank_statistic == "moderated_t") &&
         all(gsea_cc$rank_statistic == "moderated_t") &&
         all(gsea_kegg$rank_statistic == "moderated_t"))
expect("enrichment", "canonical and sensitivity rows are in separate files",
       all(gsea_go$analysis_role == "canonical") &&
         all(gsea_cc$analysis_role == "canonical") &&
         all(gsea_kegg$analysis_role == "canonical") &&
         all(sens$analysis_role == "sensitivity") &&
         all(sens$rank_statistic == "log2fc"))
expect("enrichment", "GO-CC is a primary output and part of the sensitivity file",
       all(gsea_cc$ontology == "GO_CC") && all(ora_cc$ontology == "GO_CC") &&
         any(sens$ontology == "GO_CC") && any(coverage$ontology == "GO_CC"),
       paste(nrow(gsea_cc), "GSEA and", nrow(ora_cc), "ORA GO-CC rows"))
expect("enrichment", "every enrichment row carries an analysis_role",
       all(nzchar(gsea_go$analysis_role)) && all(nzchar(gsea_cc$analysis_role)) &&
         all(nzchar(ora$analysis_role)) && all(nzchar(ora_cc$analysis_role)) &&
         all(nzchar(sens$analysis_role)))
expect("enrichment", "ORA rows are separated by query_list",
       length(unique(ora$query_list)) >= 3L && length(unique(ora_cc$query_list)) >= 3L,
       paste(length(unique(ora$query_list)), "variants"))
# The enrichment runs (GO-BP, GO-CC and the KEGG run) must have been computed on the locked
# statistics GCT.
for (run in c(RELEASE_STATISTICS_RUN$enrichment_runs, "KEGG")) {
  ix_path <- file.path(if (run == "KEGG") release_kegg_root(DATA_ROOT) else
    release_enrichment_root(run, DATA_ROOT), "indexEnrichmentComparisons.csv")
  if (!file.exists(ix_path)) {
    record("enrichment", paste(run, "enrichment index cites the locked statistics GCT"), "SKIP",
           "canonical enrichment index unreachable")
  } else {
    ix <- release_read_csv(ix_path)
    expect("enrichment", paste(run, "enrichment index cites the locked statistics GCT"),
           nrow(ix) == inv$n_primary_contrasts &&
             all(as.character(ix$source_gct_sha256) == RELEASE_LOCKED_ARTEFACTS$stat_gct$sha256),
           paste(nrow(ix), "comparisons"))
  }
}
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
expect("enrichment", "EWCE baseline rows are labeled secondary",
       all(ewce$primary_or_secondary[ewce$ewce_analysis_type == "Baseline"] == "secondary"))

check_gsea_equal <- function(tbl, root, file) {
  if (!dir.exists(root)) return(NULL)
  bad <- character(0); checked <- 0L
  for (cmp in unique(tbl$canonical_comparison)) {
    p <- file.path(root, cmp, file)
    if (!file.exists(p)) { bad <- c(bad, paste("missing:", cmp)); next }
    d <- release_read_csv(p)
    ex <- tbl[tbl$canonical_comparison == cmp, , drop = FALSE]
    if (!identical(ex$term_id, as.character(d$ID))) { bad <- c(bad, paste("order:", cmp)); next }
    if (!isTRUE(all.equal(as.numeric(ex$NES), as.numeric(d$NES), tolerance = 0)) ||
        !isTRUE(all.equal(as.numeric(ex$adjusted_p_value), as.numeric(d[["p.adjust"]]),
                          tolerance = 0))) {
      bad <- c(bad, paste("values:", cmp))
    }
    checked <- checked + 1L
  }
  list(bad = bad, checked = checked)
}
for (blk in list(
  list(tbl = gsea_go, run = "GO_BP", file = "GSEA_GO_BP.csv", label = "GO-BP"),
  list(tbl = gsea_cc, run = "GO_CC", file = "GSEA_GO_CC.csv", label = "GO-CC"),
  list(tbl = gsea_kegg, run = "KEGG", file = "GSEA_KEGG.csv", label = "KEGG"))) {
  root <- if (blk$run == "KEGG") release_kegg_root(DATA_ROOT) else
    release_enrichment_root(blk$run, DATA_ROOT)
  r <- check_gsea_equal(blk$tbl, file.path(root, "per_comparison"), blk$file)
  if (is.null(r)) {
    record("enrichment",
           paste("exported", blk$label, "GSEA statistics equal the canonical files"), "SKIP",
           "canonical enrichment folder unreachable")
  } else {
    expect("enrichment",
           paste("exported", blk$label, "GSEA statistics equal the canonical files exactly"),
           length(r$bad) == 0L,
           if (length(r$bad)) paste(head(r$bad, 3), collapse = "; ") else
             paste(r$checked, "comparisons verified"))
  }
}

# Gene members are labeled with their namespace: UniProt accessions for GO, NCBI Gene IDs
# for KEGG (the KEGG run maps proteins to Entrez IDs before testing).
go_tables <- list(gsea_go, gsea_cc, ora, ora_cc, sens[sens$ontology != "KEGG", , drop = FALSE])
kegg_tables <- list(gsea_kegg, sens[sens$ontology == "KEGG", , drop = FALSE])
namespace_ok <- function(tables, type) {
  all(vapply(tables, function(d) {
    all(d$gene_identifier_type == type) &&
      release_members_in_namespace(release_gene_members(d$core_enrichment), type)
  }, logical(1)))
}
expect("enrichment", "KEGG rows are labeled ENTREZID and their members are NCBI Gene IDs",
       nrow(gsea_kegg) > 0L && namespace_ok(kegg_tables, "ENTREZID"),
       paste(nrow(gsea_kegg), "primary and", sum(sens$ontology == "KEGG"), "sensitivity rows"))
expect("enrichment", "GO rows are labeled UNIPROT and their members are UniProt accessions",
       namespace_ok(go_tables, "UNIPROT"),
       paste(sum(vapply(go_tables, nrow, integer(1))), "rows"))
dict_ns <- read_release(rel("metadata", "data_dictionary.tsv"))
dict_ns <- dict_ns$description[dict_ns$column == "gene_identifier_type"]
expect("enrichment", "the data dictionary says KEGG members are NCBI Gene IDs",
       length(dict_ns) > 0L && all(grepl("ENTREZID (NCBI Gene ID) for KEGG", dict_ns, fixed = TRUE)),
       paste(length(dict_ns), "entries"))

# Table 2f: the primary GSEA rows with adjusted P < 0.05, with leading edges as symbols.
sig_sets <- read_release(rel("enrichment", "significant_gene_sets.tsv.gz"))
SIG_SET_COLUMNS <- c("canonical_comparison", "sample_class", "numerator_condition",
                     "denominator_condition", "ontology", "term_id", "term_name", "direction",
                     "NES", "set_size", "n_core_enrichment_genes", "p_value", "adjusted_p_value",
                     "core_enrichment_gene_symbols")
set_cols <- setdiff(SIG_SET_COLUMNS, "core_enrichment_gene_symbols")
sets_expected <- do.call(rbind, lapply(list(gsea_go, gsea_cc, gsea_kegg), function(d) {
  d[d$adjusted_p_value < 0.05, c(set_cols, "analysis_role", "core_enrichment"), drop = FALSE]
}))
sets_m <- match(paste(sig_sets$canonical_comparison, sig_sets$ontology, sig_sets$term_id),
                paste(sets_expected$canonical_comparison, sets_expected$ontology,
                      sets_expected$term_id))
expect("enrichment", "significant_gene_sets.tsv.gz has the expected columns in order",
       identical(names(sig_sets), SIG_SET_COLUMNS), paste(ncol(sig_sets), "columns"))
expect("enrichment",
       "significant_gene_sets.tsv.gz holds exactly the primary GSEA rows with adjusted P < 0.05",
       nrow(sig_sets) == nrow(sets_expected) && !anyNA(sets_m) && !anyDuplicated(sets_m) &&
         all(sets_expected$analysis_role == "canonical"),
       paste0(nrow(sig_sets), " rows (",
              paste(names(table(sig_sets$ontology)), table(sig_sets$ontology), collapse = ", "),
              ")"))
expect("enrichment", "significant_gene_sets.tsv.gz values are identical to the GSEA tables",
       !anyNA(sets_m) &&
         all(vapply(set_cols, function(cn) {
           identical(sig_sets[[cn]], sets_expected[[cn]][sets_m])
         }, logical(1))))
expect("enrichment", "every leading edge has n_core_enrichment_genes members, one per source member",
       !anyNA(sets_m) &&
         identical(lengths(release_gene_members(sig_sets$core_enrichment_gene_symbols)),
                   as.integer(sig_sets$n_core_enrichment_genes)) &&
         identical(lengths(release_gene_members(sig_sets$core_enrichment_gene_symbols)),
                   lengths(release_gene_members(sets_expected$core_enrichment[sets_m]))))
expect("enrichment", "significant gene-set direction follows the sign of NES",
       identical(sig_sets$direction,
                 paste0("enriched_toward_", ifelse(sig_sets$NES > 0, sig_sets$numerator_condition,
                                                   sig_sets$denominator_condition))))

# The gene sets fgsea could not test: the released list equals the runs' own lists, and the
# coverage table counts the same sets.
untested_src <- tryCatch(release_fgsea_untested_sets(pc$canonical_comparison, DATA_ROOT),
                         error = function(e) NULL)
if (is.null(untested_src)) {
  record("enrichment", "the untested gene sets equal the enrichment runs' lists", "SKIP",
         "canonical untested-set audits unreachable")
} else {
  key_src <- sort(paste(untested_src$canonical_comparison,
                        sub("_GSEA.*$", "", untested_src$analysis),
                        ifelse(grepl("_log2fc_sensitivity$", untested_src$analysis),
                               "sensitivity", "canonical"), untested_src$ID))
  key_rel <- sort(paste(untested_rel$canonical_comparison, untested_rel$ontology,
                        untested_rel$analysis_role, untested_rel$term_id))
  expect("enrichment", "the untested gene sets equal the enrichment runs' lists",
         identical(key_src, key_rel),
         paste(sum(untested_rel$analysis_role == "canonical"), "canonical and",
               sum(untested_rel$analysis_role == "sensitivity"), "sensitivity GO sets"))
}
cov_go <- coverage[coverage$analysis == "GSEA" & coverage$ontology %in% c("GO_BP", "GO_CC"), ,
                   drop = FALSE]
expect("enrichment", "the coverage table counts the same untested GO sets",
       all(mapply(function(cmp, ont, role, n) {
         sum(untested_rel$canonical_comparison == cmp & untested_rel$ontology == ont &
               untested_rel$analysis_role == role) == as.integer(n)
       }, cov_go$canonical_comparison, cov_go$ontology, cov_go$analysis_role,
       cov_go$gsea_untested_gene_sets)),
       paste(nrow(cov_go), "GO GSEA blocks"))
expect("enrichment", "every GSEA block records its untested gene sets",
       !anyNA(suppressWarnings(as.integer(coverage$gsea_untested_gene_sets[
         coverage$analysis == "GSEA"]))))
# The adopted enrichment carries no provisional or pending wording, and its settings are read
# from the runs.
perm_docs <- c("README_DATA.md", "provenance/REVISION_PROTEOMICS_DATA_CHANGELOG.md",
               "provenance/MANUSCRIPT_TERMINOLOGY_ACTIONS.md",
               "provenance/figure_source_map.tsv", "enrichment/enrichment_run_parameters.tsv",
               "metadata/secondary_analysis_manifest.tsv")
perm_text <- unlist(lapply(perm_docs[file.exists(rel(perm_docs))], function(f) {
  readLines(rel(f), warn = FALSE, encoding = "UTF-8")
}))
expect("enrichment",
       "the GO enrichment is recorded as adopted and the published documents contain no provisional or pending wording",
       identical(RELEASE_STATISTICS_RUN$enrichment_rerun_status, "adopted") &&
         !any(grepl("provisional|pending", perm_text, ignore.case = TRUE) &
                !grepl("PRIDE_READY_PENDING_RAW_UPLOAD|pending raw upload", perm_text,
                       ignore.case = TRUE)),
       paste(length(perm_docs), "documents scanned"))
perm <- tryCatch(release_enrichment_permutations(DATA_ROOT), error = function(e) NULL)
if (is.null(perm)) {
  record("enrichment", "the permutation statement is published as the runs record it", "SKIP",
         "canonical run parameters unreachable")
} else {
  expect("enrichment", "the permutation statement is published as the runs record it",
         any(grepl(perm$sentence, perm_text, fixed = TRUE)) &&
           any(grepl(paste0("nPermSimple = ", format(perm$go_n_perm_simple, big.mark = ",")),
                     perm_text, fixed = TRUE)),
         perm$setting)
}

# --------------------------------------------------------------------------------------
cat("\n=== effect-size semantics ===\n")
# --------------------------------------------------------------------------------------
# Effect size: the right name, traceable to its source column, identical values, and no
# reader-facing document calling it a fold change.

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
expect("effect_size", "the definition states that it is not a log2 fold change",
       all(grepl("NOT a log2 fold change", ds_eff$effect_size_definition, fixed = TRUE)))
expect("effect_size", "the effect size is not called a standardized mean difference alone",
       !any(grepl("standardi[sz]ed mean difference", ds_eff$effect_size_definition)))

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

# The terminology audit, then a scan of every reader-facing document (not the builder
# sources, which hold the claim list and the audit rows).
audit_path <- rel("provenance", "EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv")
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
  "provenance/REVISION_PROTEOMICS_DATA_CHANGELOG.md",
  "provenance/MANUSCRIPT_TERMINOLOGY_ACTIONS.md",
  "provenance/figure_source_map.tsv",
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

# README_DATA.md must contain the Methods sentence for the effect size.
readme_lines <- readLines(rel("README_DATA.md"), warn = FALSE)
expect("effect_size", "README_DATA.md carries the methods-safe definition",
       any(grepl(RELEASE_EFFECT_SIZE$methods_sentence, readme_lines, fixed = TRUE)))
expect("effect_size", "the manuscript action note lists the out-of-repo corrections",
       file.exists(rel("provenance", "MANUSCRIPT_TERMINOLOGY_ACTIONS.md")))

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
         all(gsea_cc$analysis_unit == "animal") && all(ora$analysis_unit == "animal") &&
         all(ora_cc$analysis_unit == "animal") && all(ewce$analysis_unit == "animal"))
expect("statistical unit", "no inferential row claims 6 replicates per arm",
       !any(as.integer(da$n_numerator_animals) > 3L, na.rm = TRUE) &&
         !any(as.integer(gsea_go$n_numerator_animals) > 3L, na.rm = TRUE))

fm <- read_release(rel("provenance", "figure_source_map.tsv"))
acq_unit <- grepl("hemisphere|acquisition", fm$statistical_unit, ignore.case = TRUE)
expect("statistical unit",
       "panels drawn at acquisition level are labeled technical QC",
       all(fm$primary_or_secondary[acq_unit] == "technical_qc"),
       paste(sum(acq_unit), "acquisition-level panel(s)"))
expect("figures", "all required manuscript panels are mapped (manuscript numbering)",
       all(c(paste("Fig. 2", c("C", "D", "E", "F", "G", "H", "I")),
             paste("Supp. Fig. 4", c("A", "B1", "B2", "B3", "C", "D", "E", "F"))) %in%
             fm$panel_label))
rebuilt <- fm[grepl("^REBUILT_WITHIN_COMPARTMENT", fm$revision_status), , drop = FALSE]
expect("figures", "the panels that use the differential statistics come from the per-sample-class rebuild",
       setequal(rebuilt$panel_label,
                c(paste("Fig. 2", c("D", "F", "G", "H", "I")), "Supp. Fig. 4 F")) &&
         all(grepl(paste(RELEASE_STATISTICS_RUN$figure_rebuild_relative, collapse = "/"),
                   rebuilt$panel_source_data_canonical_path, fixed = TRUE)),
       paste(nrow(rebuilt), "rebuilt panel(s)"))
# The published panel files are reader-facing CSVs written by stage 06 from the canonical
# panel files. The canonical files are not copied: the map names each one with its pinned
# SHA256, and the registry records that file and hash as the source of the published CSV.
split_cell <- function(x) {
  x <- as.character(x)
  if (is.na(x) || !nzchar(x)) character(0) else strsplit(x, ";", fixed = TRUE)[[1]]
}
registry <- release_registry_read(OUT_ROOT)
panel_rows <- which(!is.na(fm$panel_source_data_release_path) &
                      nzchar(fm$panel_source_data_release_path))
panel_published_ok <- vapply(panel_rows, function(i) {
  files <- split_cell(fm$panel_source_data_release_path[i])
  canon <- split_cell(fm$panel_source_data_canonical_path[i])
  hashes <- split_cell(fm$panel_source_data_sha256[i])
  length(files) > 0L && length(files) == length(canon) && length(canon) == length(hashes) &&
    all(grepl("^editor_source_data/figure_source_data/[^/]+[.]csv$", files)) &&
    all(file.exists(file.path(OUT_ROOT, files))) &&
    all(vapply(seq_along(files), function(j) {
      r <- registry[registry$relative_path == files[[j]], , drop = FALSE]
      nrow(r) == 1L && identical(r$canonical_source_sha256, hashes[[j]]) &&
        identical(basename(r$source_artifact), basename(canon[[j]]))
    }, logical(1)))
}, logical(1))
expect("figures",
       "every published panel file exists and is registered with its canonical source and pinned hash",
       length(panel_rows) > 0L && all(panel_published_ok),
       if (all(panel_published_ok)) paste(length(panel_rows), "panels") else
         paste(fm$panel_label[panel_rows][!panel_published_ok], collapse = "; "))
canon_files <- unique(data.frame(
  path = unlist(lapply(fm$panel_source_data_canonical_path[panel_rows], split_cell)),
  sha256 = unlist(lapply(fm$panel_source_data_sha256[panel_rows], split_cell)),
  stringsAsFactors = FALSE))
if (!all(file.exists(canon_files$path))) {
  record("figures", "every canonical panel file still matches the SHA256 the map records",
         "SKIP", "canonical panel files unreachable")
} else {
  canon_ok <- vapply(seq_len(nrow(canon_files)), function(i) {
    identical(release_sha256(canon_files$path[[i]]), canon_files$sha256[[i]])
  }, logical(1))
  expect("figures", "every canonical panel file still matches the SHA256 the map records",
         all(canon_ok),
         if (all(canon_ok)) paste(nrow(canon_files), "canonical files") else
           paste(basename(canon_files$path[!canon_ok]), collapse = "; "))
}
fsd_dir <- rel("editor_source_data", "figure_source_data")
fsd_files <- list.files(fsd_dir)
fsd_named <- basename(unlist(lapply(fm$panel_source_data_release_path[panel_rows], split_cell)))
expect("figures", "the panel source-data folder holds only the mapped CSV files (no .xlsx)",
       length(fsd_files) > 0L && all(grepl("[.]csv$", fsd_files)) &&
         setequal(fsd_files, fsd_named),
       paste(length(fsd_files), "files"))
fsd <- lapply(stats::setNames(file.path(fsd_dir, fsd_files), fsd_files),
              release_read_published_csv)
fsd_headers <- unlist(lapply(fsd, names), use.names = FALSE)
expect("figures", "panel source-data headers are plain (no underscores, log2fc, pval, padj or alias)",
       all(nzchar(fsd_headers)) &&
         !any(grepl("_|log2fc|logfc|pval|padj|alias", fsd_headers, ignore.case = TRUE)),
       paste(length(fsd_headers), "headers in", length(fsd), "files"))
# Quantities that are also in Supplementary Table 2 carry its display headers (README_DATA
# section 12): every file has exactly the headers pinned in RELEASE_FIGURE_CSV_HEADERS, and no
# header is a workbook header written differently (case, spaces, hyphens, punctuation).
fsd_bad <- names(fsd)[!vapply(names(fsd), function(f) {
  identical(names(fsd[[f]]), RELEASE_FIGURE_CSV_HEADERS[[f]])
}, logical(1))]
fsd_wb_inv <- read_release(rel("editor_source_data", "workbook_sheet_inventory.tsv"),
                           colClasses = "character", na.strings = character(0))
fsd_variants <- release_figure_header_variants(
  fsd_headers, unlist(strsplit(fsd_wb_inv$display_columns, ";", fixed = TRUE)))
expect("figures",
       "panel source-data headers are the pinned set, with the workbook display header for every workbook quantity",
       setequal(names(fsd), names(RELEASE_FIGURE_CSV_HEADERS)) && !length(fsd_bad) &&
         !length(fsd_variants),
       if (length(fsd_bad) || length(fsd_variants)) {
         paste(c(fsd_bad, fsd_variants), collapse = "; ")
       } else paste(length(unique(fsd_headers)), "distinct headers in", length(fsd), "files"))
foreign <- "RFID|Exp9|ActivityIndex|AnimalPos"
foreign_hits <- names(fsd)[vapply(fsd, function(d) {
  any(grepl(foreign, c(names(d), unlist(d, use.names = FALSE))))
}, logical(1))]
expect("figures", "no panel source-data file holds data of another project",
       length(foreign_hits) == 0L,
       if (length(foreign_hits)) paste(foreign_hits, collapse = "; ") else
         paste(length(fsd), "files scanned"))
supp_a <- fsd[["SuppFig4A_source_data.csv"]]
expect("figures", "Supp. Fig. 4 A source data are this project's 96 acquisitions, by the Sample IDs of table 2a",
       !is.null(supp_a) && nrow(supp_a) == inv$n_measurement_records &&
         !anyDuplicated(supp_a[["Sample ID"]]) &&
         setequal(supp_a[["Sample ID"]], sm$sample_id) &&
         identical(names(supp_a), RELEASE_FIGURE_CSV_HEADERS[["SuppFig4A_source_data.csv"]]) &&
         all(supp_a[["Animal ID"]] == sm$AnimalID[match(supp_a[["Sample ID"]], sm$sample_id)]),
       if (is.null(supp_a)) "file missing" else
         paste(nrow(supp_a), "rows x", ncol(supp_a), "columns"))
expect("figures", "every data-bearing panel names a final revision script",
       all(nzchar(fm$final_revision_script[fm$primary_or_secondary != "not_applicable"])))
released_panel_files <- basename(unlist(strsplit(
  fm$panel_source_data_release_path[!is.na(fm$panel_source_data_release_path) &
                                      nzchar(fm$panel_source_data_release_path)],
  ";", fixed = TRUE)))
expect("figures", "released panel source data are named in the manuscript numbering",
       length(released_panel_files) > 0L &&
         all(grepl("^(Fig2|SuppFig4)[A-Z][0-9]?_", released_panel_files)),
       paste(length(released_panel_files), "files"))
expect("figures", "no rebuilt panel describes its effect size as moderated or as log2FC",
       !any(grepl("log2FC|moderated limma", c(rebuilt$statistical_unit, rebuilt$panel_title),
                  ignore.case = TRUE)))

# --------------------------------------------------------------------------------------
cat("\n=== Supplementary Table 2 workbook ===\n")
# --------------------------------------------------------------------------------------
# The inventory (written by stage 07) is checked against the expected layout, and the built
# .xlsx against the inventory and the release tables.
wb_inv <- read_release(rel("editor_source_data", "workbook_sheet_inventory.tsv"),
                       colClasses = "character", na.strings = character(0))
sheet_cols <- stats::setNames(strsplit(wb_inv$columns, ";", fixed = TRUE), wb_inv$sheet)
display_cols <- stats::setNames(strsplit(wb_inv$display_columns, ";", fixed = TRUE),
                                wb_inv$sheet)
hygiene <- release_check_workbook_hygiene(sheet_cols)
expect("workbook", "no duplicate or unnamed workbook columns", length(hygiene) == 0L,
       if (length(hygiene)) paste(head(hygiene, 2), collapse = "; ") else
         paste(nrow(wb_inv), "sheets"))
hygiene_display <- release_check_workbook_hygiene(display_cols)
expect("workbook", "no duplicate or unnamed workbook DISPLAY headers",
       length(hygiene_display) == 0L,
       if (length(hygiene_display)) paste(head(hygiene_display, 2), collapse = "; ") else
         "display headers unique and named on every sheet")
expect("workbook", "every sheet has as many display headers as machine columns",
       all(lengths(display_cols) == lengths(sheet_cols)) &&
         all(lengths(display_cols) == as.integer(wb_inv$n_columns)))
expect("workbook", "no display header leaks a snake_case machine name",
       !any(grepl("_", unlist(display_cols, use.names = FALSE), fixed = TRUE)),
       paste(length(unlist(display_cols, use.names = FALSE)), "headers checked"))
layout <- RELEASE_WORKBOOK_LAYOUT
expect("workbook",
       "the sheets are README, Column definitions and sub-tables 2a to 2q, in that order",
       identical(wb_inv$sheet_key, layout$sheet_key) &&
         identical(wb_inv$supplementary_table, layout$letter),
       paste(wb_inv$supplementary_table[nzchar(wb_inv$supplementary_table)], collapse = " "))
expect("workbook", "every sheet has a title", all(nzchar(wb_inv$title)))
xlsx_path <- rel(RELEASE_WORKBOOK_REL)
xlsx_found <- list.files(rel("editor_source_data"), pattern = "[.]xlsx$")
expect("workbook", paste("the workbook is", RELEASE_WORKBOOK_REL, "and no other workbook is released"),
       file.exists(xlsx_path) && identical(xlsx_found, basename(RELEASE_WORKBOOK_REL)),
       paste(xlsx_found, collapse = "; "))
# The revision records and the figure map are in provenance/, not next to the table.
editor_files <- list.files(rel("editor_source_data"))
expect("workbook", "editor_source_data holds only Supplementary Table 2, its inventory and the figure source data",
       setequal(editor_files, c(basename(RELEASE_WORKBOOK_REL), "workbook_sheet_inventory.tsv",
                                "figure_source_data")),
       paste(editor_files, collapse = "; "))
expect("workbook", "no sheet carries a misleading public gene_symbol field",
       !any(vapply(sheet_cols, function(cs) {
         "gene_symbol" %in% cs && !"uniprot_accession" %in% cs
       }, logical(1))))
internal_cols <- unique(grep(RELEASE_WORKBOOK_INTERNAL_COLUMNS, unlist(sheet_cols), value = TRUE))
expect("workbook", "no contrast, alias, path, hash, curation or revision column reaches the workbook",
       length(internal_cols) == 0L, paste(internal_cols, collapse = ", "))
keyed <- vapply(sheet_cols, function(cs) "canonical_comparison" %in% cs, logical(1))
expect("workbook", "comparison-keyed sheets show Comparison first and Comparison key last",
       sum(keyed) > 0L &&
         all(vapply(names(which(keyed)), function(s) {
           cs <- sheet_cols[[s]]; ds <- display_cols[[s]]
           cs[[1]] == "comparison_label" && cs[[length(cs)]] == "canonical_comparison" &&
             ds[[1]] == "Comparison" && ds[[length(ds)]] == "Comparison key"
         }, logical(1))),
       paste(sum(keyed), "sheets"))
inv_rows <- stats::setNames(as.integer(wb_inv$n_rows), wb_inv$sheet_key)
expected_rows <- c(
  Samples = nrow(sm), Animal_Profiles = nrow(al), Comparisons = nrow(pc),
  Significant_Proteins = sum(as.integer(da_sum$n_significant_fdr_0_05_mapped_only)),
  All_Proteins = nrow(da),
  Significant_Gene_Sets = sum(c(gsea_go$adjusted_p_value, gsea_cc$adjusted_p_value,
                                gsea_kegg$adjusted_p_value) < 0.05),
  GSEA_GO_BP = nrow(gsea_go), GSEA_GO_CC = nrow(gsea_cc), GSEA_KEGG = nrow(gsea_kegg),
  ORA_GO_BP = sum(ora$analysis_role == "canonical"),
  ORA_GO_CC = sum(ora_cc$analysis_role == "canonical"),
  GSEA_Untested_Sets = sum(untested_rel$analysis_role == "canonical"), EWCE = nrow(ewce))
row_diff <- names(expected_rows)[inv_rows[names(expected_rows)] != expected_rows]
expect("workbook",
       "sheet row counts equal their release tables (2j, 2k and 2l: primary analysis only)",
       length(row_diff) == 0L,
       if (length(row_diff)) paste(row_diff, collapse = ", ") else
         paste(sum(expected_rows), "rows on", length(expected_rows), "sheets"))

if (!file.exists(xlsx_path)) {
  record("workbook", "the workbook file can be read", "FAIL", "file missing")
} else {
  wbx <- release_read_workbook_xml(xlsx_path)
  record("workbook", "workbook size", "PASS",
         sprintf("%.1f MB, %d tabs", file.info(xlsx_path)$size / 1024^2, length(wbx$tabs)))
  expect("workbook", "tabs are Contents, then the inventory sheets in order",
         identical(wbx$tabs, c("Contents", wb_inv$sheet)), paste(wbx$tabs, collapse = " | "))
  tabs_all <- c("Contents", wb_inv$sheet)
  if (all(tabs_all %in% names(wbx$sheets))) {
    row_cells <- function(tab, r) {
      v <- unname(wbx$sheets[[tab]]$rows[[r]])
      v[!is.na(v)]
    }
    titles <- c(paste0("Supplementary Table ", RELEASE_WORKBOOK_TABLE_NUMBER, " | Contents"),
                paste0("Supplementary Table ", RELEASE_WORKBOOK_TABLE_NUMBER,
                       sub("^[0-9]+", "", wb_inv$supplementary_table), " | ", wb_inv$title))
    bad_title <- tabs_all[!vapply(seq_along(tabs_all), function(i) {
      identical(row_cells(tabs_all[[i]], 1L), titles[[i]])
    }, logical(1))]
    expect("workbook", "row 1 of every sheet reads 'Supplementary Table 2<letter> | <title>'",
           length(bad_title) == 0L, paste(bad_title, collapse = "; "))
    headers <- c(list(c("Table", "Sheet", "Title", "Rows", "Columns", "Description")),
                 unname(display_cols))
    bad_header <- tabs_all[!vapply(seq_along(tabs_all), function(i) {
      identical(row_cells(tabs_all[[i]], RELEASE_WORKBOOK_HEADER_ROW), headers[[i]])
    }, logical(1))]
    expect("workbook", "the header is on row 4 of every sheet and matches the inventory",
           length(bad_header) == 0L, paste(bad_header, collapse = "; "))
    n_data <- c(nrow(wb_inv), as.integer(wb_inv$n_rows))
    bad_rows <- tabs_all[vapply(seq_along(tabs_all), function(i) {
      wbx$sheets[[tabs_all[[i]]]]$last_row - RELEASE_WORKBOOK_HEADER_ROW != n_data[[i]]
    }, logical(1))]
    expect("workbook", "every sheet holds as many data rows as the inventory records",
           length(bad_rows) == 0L, paste(bad_rows, collapse = "; "))
    bad_layout <- tabs_all[!vapply(tabs_all, function(tab) {
      s <- wbx$sheets[[tab]]
      identical(s$frozen_at, "B5") && "A2:H2" %in% s$merges
    }, logical(1))]
    expect("workbook", "panes frozen below the header and right of column 1; caption merged over A2:H2",
           length(bad_layout) == 0L, paste(bad_layout, collapse = "; "))
    # The Sheet column holds the names as text, each with a hyperlink element to its sheet.
    contents_x <- wbx$sheets[["Contents"]]
    link_targets <- sub("^'([^']*)'!A1$", "\\1", contents_x$hyperlinks)
    contents_sheet <- openxlsx::read.xlsx(xlsx_path, sheet = "Contents",
                                          startRow = RELEASE_WORKBOOK_HEADER_ROW,
                                          check.names = FALSE)
    expect("workbook", "the Contents sheet links every sheet name (text, not a formula) to its sheet",
           identical(link_targets, wb_inv$sheet) &&
             identical(as.character(contents_sheet$Sheet), wb_inv$sheet) &&
             identical(contents_x$hyperlink_refs,
                       paste0("B", RELEASE_WORKBOOK_HEADER_ROW + seq_len(nrow(wb_inv)))) &&
             !length(contents_x$formulas),
           paste(length(link_targets), "links"))
  } else {
    record("workbook", "every inventory sheet is in the workbook", "FAIL",
           paste(setdiff(tabs_all, names(wbx$sheets)), collapse = "; "))
  }
  expect("workbook", "every text cell is a shared string, so the text scan covers it",
         !any(vapply(wbx$sheets, `[[`, logical(1), "inline_strings")))
  wb_text <- c(wbx$tabs, wbx$shared_strings,
               unlist(lapply(wbx$sheets, `[[`, "formulas"), use.names = FALSE))
  internal_hits <- release_internal_text_hits(wb_text)
  expect("workbook",
         "no workbook cell holds a drive path (S:/, D:\\, C:/), a hash, a script or folder name, 'Figure3' or the name Fabian",
         nrow(internal_hits) == 0L,
         if (nrow(internal_hits)) paste0(internal_hits$kind[1], ": ",
                                         substr(internal_hits$text[1], 1, 80)) else
           paste(length(wb_text), "strings scanned"))
  wb_term_names <- release_workbook_term_names(xlsx_path, wb_inv)
  british_hits <- release_british_spelling_hits(wbx$shared_strings, wb_term_names)
  expect("workbook",
         "no workbook text uses British spelling (GO and KEGG term names in tables 2f to 2l excepted)",
         length(british_hits) == 0L,
         if (length(british_hits)) paste0(length(british_hits), " string(s), e.g. ",
                                          substr(british_hits[1], 1, 80)) else
           paste(length(wbx$shared_strings), "strings scanned,", length(wb_term_names),
                 "term names excepted"))

  read_sheet <- function(key) {
    i <- match(key, wb_inv$sheet_key)
    d <- openxlsx::read.xlsx(xlsx_path, sheet = wb_inv$sheet[[i]],
                             startRow = RELEASE_WORKBOOK_HEADER_ROW, check.names = FALSE,
                             sep.names = " ")
    if (!identical(names(d), display_cols[[i]])) {
      stop("Sheet ", wb_inv$sheet[[i]], " does not have the inventory headers.", call. = FALSE)
    }
    names(d) <- sheet_cols[[i]]
    d
  }
  same_num <- release_same_numbers

  readme_sheet <- read_sheet("README")
  readme_flat <- paste(readme_sheet$statement, collapse = " ")
  readme_forbidden <- c("submission", "reviewer", "editor", "changelog", "superseded",
                        "canonical", "03_output")
  readme_bad <- readme_forbidden[vapply(readme_forbidden, grepl, logical(1), x = readme_flat,
                                        ignore.case = TRUE)]
  expect("workbook", "the README sheet is written for the reader and ends with the release version",
         length(readme_bad) == 0L &&
           grepl("^Release of [0-9]{4}-[0-9]{2}-[0-9]{2}[.]$",
                 readme_sheet$statement[readme_sheet$item == "Version"]),
         if (length(readme_bad)) paste("found:", paste(readme_bad, collapse = ", ")) else
           paste(nrow(readme_sheet), "items"))

  coldef <- read_sheet("Column_Definitions")
  lettered <- nzchar(wb_inv$supplementary_table)
  shown <- unlist(lapply(which(lettered), function(i) {
    paste(wb_inv$supplementary_table[[i]], display_cols[[i]])
  }), use.names = FALSE)
  expect("workbook", "Column definitions has exactly one defined row per header of 2a to 2q",
         setequal(paste(coldef$table, coldef$column), shown) &&
           !anyDuplicated(paste(coldef$table, coldef$column)) &&
           all(nzchar(coldef$definition) & !is.na(coldef$definition)),
         paste(nrow(coldef), "definitions"))

  cmp_sheet <- read_sheet("Comparisons")
  m <- match(cmp_sheet$canonical_comparison, da_sum$canonical_comparison)
  count_cols <- intersect(c("n_proteins_tested", "n_significant_fdr_0_05",
                            "n_significant_higher_in_numerator",
                            "n_significant_higher_in_denominator",
                            "n_proteins_tested_mapped_only", "n_significant_fdr_0_05_mapped_only",
                            "n_significant_higher_in_numerator_mapped_only",
                            "n_significant_higher_in_denominator_mapped_only"),
                          names(cmp_sheet))
  expect("workbook", "2c counts equal the differential summary",
         nrow(cmp_sheet) == inv$n_primary_contrasts && !anyNA(m) && length(count_cols) == 8L &&
           all(vapply(count_cols, function(cn) {
             identical(as.integer(cmp_sheet[[cn]]), as.integer(da_sum[[cn]][m]))
           }, logical(1))),
         paste(nrow(cmp_sheet), "comparisons,", length(count_cols), "counts"))

  sig_sheet <- read_sheet("Significant_Proteins")
  expected_sig <- da[as.logical(da$significant_fdr_0_05) %in% TRUE &
                       da$id_mapping_status == "mapped", , drop = FALSE]
  m <- match(paste(sig_sheet$canonical_comparison, sig_sheet$protein_group_id),
             paste(expected_sig$canonical_comparison, expected_sig$protein_group_id))
  expect("workbook", "2d holds exactly the FDR < 0.05 annotated rows of the protein table, values unchanged",
         nrow(sig_sheet) == nrow(expected_sig) && !anyNA(m) && !anyDuplicated(m) &&
           same_num(sig_sheet$effect_size_sd_units, expected_sig$effect_size_sd_units[m]) &&
           same_num(sig_sheet$moderated_t, expected_sig$moderated_t[m]) &&
           same_num(sig_sheet$P.Value, expected_sig$P.Value[m]) &&
           same_num(sig_sheet$adj.P.Val, expected_sig$adj.P.Val[m]) &&
           all(grepl("^higher in ", sig_sheet$direction)),
         paste(nrow(sig_sheet), "rows"))

  sets_sheet <- read_sheet("Significant_Gene_Sets")
  gsea_all <- do.call(rbind, lapply(list(gsea_go, gsea_cc, gsea_kegg), function(d) {
    d[d$adjusted_p_value < 0.05, c("canonical_comparison", "term_id", "NES", "p_value",
                                   "adjusted_p_value", "n_core_enrichment_genes")]
  }))
  m <- match(paste(sets_sheet$canonical_comparison, sets_sheet$term_id),
             paste(gsea_all$canonical_comparison, gsea_all$term_id))
  expect("workbook", "2f holds exactly the primary GSEA rows with adjusted P < 0.05, values unchanged",
         nrow(sets_sheet) == nrow(gsea_all) && !anyNA(m) && !anyDuplicated(m) &&
           same_num(sets_sheet$NES, gsea_all$NES[m]) &&
           same_num(sets_sheet$p_value, gsea_all$p_value[m]) &&
           same_num(sets_sheet$adjusted_p_value, gsea_all$adjusted_p_value[m]) &&
           identical(lengths(release_gene_members(sets_sheet$core_enrichment_gene_symbols)),
                     as.integer(gsea_all$n_core_enrichment_genes[m])),
         paste(nrow(sets_sheet), "rows"))
}

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
                 rel("provenance", "REVISION_PROTEOMICS_DATA_CHANGELOG.md"),
                 rel("pride", "README_PRIDE.md"),
                 rel("pride", "SDRF_MISSING_METADATA.md"),
                 rel("provenance", "UPSTREAM_PREPROCESSING_GAP.md"),
                 rel("metadata", "data_dictionary.tsv"),
                 rel("metadata", "metadata_field_provenance.tsv"),
                 rel("pride", "sdrf_field_status.tsv"))
forbidden <- c(tolower(FORBIDDEN_PLATE_WORDING), tolower(RELEASE_PLATE_RETIRED_WORDING))
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
# The reader documents are written in US spelling, like the workbook.
us_docs <- c("README_DATA.md", "metadata/data_dictionary.tsv", "pride/README_PRIDE.md")
us_doc_hits <- unlist(lapply(us_docs, function(f) {
  if (!file.exists(rel(f))) return(paste(f, "missing"))
  h <- release_british_spelling_hits(readLines(rel(f), warn = FALSE, encoding = "UTF-8"))
  if (length(h)) paste0(f, ": ", substr(h, 1, 80))
}))
expect("documentation",
       "README_DATA.md, data_dictionary.tsv and README_PRIDE.md use no British spelling",
       length(us_doc_hits) == 0L,
       if (length(us_doc_hits)) paste(head(us_doc_hits, 2), collapse = "; ") else
         paste(length(us_docs), "documents scanned"))

# --------------------------------------------------------------------------------------
cat("\n=== data lineage ===\n")
# --------------------------------------------------------------------------------------
lineage <- read_release(rel("provenance", "data_lineage.tsv"))
expect("lineage", "lineage covers acquisition through release",
       all(c("A02_acquisition_raw", "A04_protein_group_matrix_prefilter",
             "A07_processed_matrix", "A08_animal_level_matrix", "A09_differential_statistics",
             "A10_split_differential_tables", "A12_mapped_differential_tables",
             "A13_enrichment", "A14_ewce", "A16c_figure_panels_within_compartment",
             "A17_publication_release") %in% lineage$artifact_id))
unresolved <- lineage$artifact_id[grepl("UNRESOLVED",
                                        paste(lineage$processing_stage, lineage$notes))]
record("lineage", "unproven lineage edges are marked UNRESOLVED, not inferred", "PASS",
       paste(length(unresolved), "artifact(s):", paste(unresolved, collapse = ", ")))
expect("lineage", "the locked artifacts appear with their hashes",
       any(lineage$sha256 == RELEASE_LOCKED_ARTEFACTS$animal_level_input_gct$sha256,
           na.rm = TRUE) &&
         any(lineage$sha256 == RELEASE_LOCKED_ARTEFACTS$stat_gct$sha256, na.rm = TRUE))
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
if (!inherits(stat_run, "error")) {
  canon_sw <- sw[sw$applies_to_stage == "canonical animal-level differential statistics", ,
                 drop = FALSE]
  expect("lineage", "the canonical statistics are attributed to limma at the run-recorded version",
         any(canon_sw$component == "limma" & canon_sw$version == stat_run$limma_version) &&
           any(canon_sw$component == "R" & canon_sw$version == stat_run$r_version),
         paste0("limma ", stat_run$limma_version, ", R ", stat_run$r_version))
  expect("lineage", "no ProTigy version is scoped to the canonical statistics",
         !any(grepl("ProTigy", canon_sw$component, ignore.case = TRUE)))
}

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
cat("\n=== repository files are named by relative paths ===\n")
# --------------------------------------------------------------------------------------
# Release tables name repository files relative to proteomic_analysis/, so no file may carry
# the local folder of the computer that built the release.
repo_prefix <- tolower(paste0(sub("/+$", "", gsub("\\\\", "/", REPO_ROOT)), "/"))
text_files <- gsub("\\\\", "/", list.files(OUT_ROOT, pattern = "[.](tsv|md|txt|csv)$|[.]tsv[.]gz$",
                                           recursive = TRUE))
local_hits <- text_files[vapply(text_files, function(f) {
  txt <- gsub("\\\\", "/", readLines(file.path(OUT_ROOT, f), warn = FALSE), useBytes = TRUE)
  txt <- tryCatch(tolower(txt), error = function(e) txt)
  any(grepl(repo_prefix, txt, fixed = TRUE, useBytes = TRUE))
}, logical(1))]
expect("lineage", "no release file names the local repository folder",
       length(local_hits) == 0L,
       if (length(local_hits)) paste(head(local_hits, 3), collapse = "; ") else
         paste(length(text_files), "text files scanned"))
path_cells <- c(lineage$path_or_release_filename,
                read_release(rel("provenance", "analysis_parameters.tsv"))$evidence_path,
                read_release(rel("metadata", "metadata_field_provenance.tsv"))$current_source,
                if (file.exists(manifest_path)) read_release(manifest_path)$source_artifact)
path_refs <- unique(unlist(strsplit(path_cells[!is.na(path_cells)], ";", fixed = TRUE)))
repo_refs <- path_refs[release_is_repo_path(path_refs, REPO_ROOT)]
missing_refs <- repo_refs[!file.exists(release_repo_resolve(repo_refs, REPO_ROOT))]
expect("lineage", "repository paths in the provenance tables resolve under proteomic_analysis/",
       length(repo_refs) > 0L && length(missing_refs) == 0L,
       if (length(missing_refs)) paste("not found:", paste(head(missing_refs, 3), collapse = "; "))
       else paste(length(repo_refs), "repository paths"))

# --------------------------------------------------------------------------------------
cat("\n=== publication scripts perform no scientific recomputation ===\n")
# --------------------------------------------------------------------------------------
scripts <- list.files(LAYER, pattern = "[.]R$", full.names = TRUE, recursive = TRUE)
hits <- release_scan_prohibited_calls(scripts)
# This file and R/release_validation.R name the calls in the scanner's pattern table; skip
# them.
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
sdrf_v <- read_release(rel("pride", "sdrf.tsv"))

# sdrf_field_status.tsv lists two kinds of required gap: SDRF columns that cannot be filled,
# and the raw acquisition files, which deposition needs but which are not an SDRF column.
# They are counted separately; either keeps the status below PRIDE_READY.
fs_is_sdrf_column <- fs$sdrf_field %in% names(sdrf_v)
fs_missing_req <- fs$status == "MISSING_REQUIRED_METADATA"
n_req_sdrf <- sum(fs_missing_req & fs_is_sdrf_column)
n_req_raw <- sum(fs_missing_req & !fs_is_sdrf_column)
n_req_total <- n_req_sdrf + n_req_raw
n_opt <- sum(fs$status == "MISSING_OPTIONAL_METADATA")
n_ready <- sum(fs$status == "READY")

expect("pride", "PRIDE status is one of the four sanctioned values",
       pride_status %in% release_pride_statuses, pride_status)
expect("pride", "SDRF exists with one row per acquisition",
       nrow(sdrf_v) == inv$n_measurement_records)
expect("pride", "every SDRF field has a recorded status", all(nzchar(fs$status)))
expect("pride", "no required field is left blank without a reason",
       all(fs$status[fs$sdrf_requirement == "required" &
                       fs$value_or_reason == "not available"] ==
             "MISSING_REQUIRED_METADATA"))
record("pride", "SDRF field readiness", "PASS",
       paste0(n_ready, " READY, ", n_req_sdrf, " MISSING_REQUIRED, ", n_opt,
              " MISSING_OPTIONAL (SDRF columns only)"))
expect("pride", "deposition is not overstated as PRIDE_READY",
       !(pride_status == "PRIDE_READY" && n_req_total > 0L))
record("pride", "missing required SDRF fields", if (n_req_sdrf == 0L) "PASS" else "SKIP",
       if (n_req_sdrf == 0L) "none" else
         paste(fs$sdrf_field[fs_missing_req & fs_is_sdrf_column], collapse = "; "))
record("pride", "raw acquisition files are a separate required deposition blocker",
       if (n_req_raw == 0L) "PASS" else "SKIP",
       if (n_req_raw == 0L) "none outstanding" else
         paste0(n_req_raw, " non-SDRF required deposition input(s) outstanding: ",
                paste(fs$sdrf_field[fs_missing_req & !fs_is_sdrf_column], collapse = "; ")))

# The sample-class correction is not a PRIDE gap: the SDRF factor must be filled from the
# analysis-time class, and the PRIDE documents must describe the correction as resolved.
class_col <- "factor value[sample class]"
cols_named <- function(df, field) df[, which(names(df) == field), drop = FALSE]
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

expect("pride", "SDRF technology type is explicit",
       all(sdrf_v[["technology type"]] == "proteomic profiling by mass spectrometry"))
expect("pride", "CeM is represented by the supported parent anatomy term",
       all(sdrf_v[["characteristics[organism part]"]] ==
             release_metadata_value(experimenter_metadata, "organism_part", TRUE)))
expect("pride", "marker-defined LCM classes do not claim unsupported cell ontology terms",
       all(sdrf_v[["characteristics[cell type]"]] == "not applicable"))
expected_sampling <- stats::setNames(
  paste0(release_metadata_value(experimenter_metadata, "organism_part"), ": ",
         experimenter_metadata$sdrf_value[
           experimenter_metadata$category == "sample_class_definition"]),
  experimenter_metadata$sample_class[
    experimenter_metadata$category == "sample_class_definition"])
expect("pride", "sampling site preserves each marker-defined LCM category",
       all(sdrf_v[["characteristics[sampling site]"]] ==
             unname(expected_sampling[sdrf_v[[class_col]]])))
expect("pride", "all 96 SDRF rows describe Mus musculus C57BL/6J samples",
       nrow(sdrf_v) == 96L &&
         all(sdrf_v[["characteristics[organism]"]] == "Mus musculus") &&
         all(sdrf_v[["characteristics[strain or breed]"]] == "C57BL/6J"))
expect("pride", "AnimalID-level sex is not invented",
       all(sdrf_v[["characteristics[sex]"]] == "not available"))
expect("pride", "age at experiment start is not relabeled as age at tissue collection",
       all(sdrf_v[["characteristics[age]"]] == "not available") &&
         !any(grepl("7-10", sdrf_v[["characteristics[age]"]], fixed = TRUE)))
expect("pride", "label-free proteomics is encoded with the SDRF term",
       all(sdrf_v[["comment[label]"]] == "label free sample"))
cleavage <- cols_named(sdrf_v, "comment[cleavage agent details]")
expect("pride", "Lys-C and trypsin are two ordered cleavage-agent columns",
       ncol(cleavage) == 2L &&
         all(cleavage[[1]] == "NT=Lys-C;AC=MS:1001309") &&
         all(cleavage[[2]] == "NT=Trypsin;AC=MS:1001251"))
expect("pride", "TCEP and CAA preparation reagents are encoded without inventing search modifications",
       all(sdrf_v[["comment[reduction reagent]"]] ==
             release_metadata_value(experimenter_metadata, "reduction_reagent", TRUE)) &&
         all(sdrf_v[["comment[alkylation reagent]"]] ==
               release_metadata_value(experimenter_metadata, "alkylation_reagent", TRUE)) &&
         all(sdrf_v[["comment[modification parameters]"]] == "not available"))
for (unknown_field in c("comment[instrument]", "comment[proteomics data acquisition method]",
                        "comment[precursor mass tolerance]", "comment[fragment mass tolerance]")) {
  expect("pride", paste(unknown_field, "is recorded as unresolved"),
         all(sdrf_v[[unknown_field]] == "not available"))
}
repeated_headers <- unique(names(sdrf_v)[duplicated(names(sdrf_v))])
expect("pride", "repeated SDRF headers are limited to the verified enzymes and templates",
       setequal(repeated_headers,
                c("comment[cleavage agent details]", "comment[sdrf template]")) &&
         sum(names(sdrf_v) == "comment[cleavage agent details]") == 2L &&
         sum(names(sdrf_v) == "comment[sdrf template]") == 2L)
protocol_txt <- paste(readLines(rel("pride", "SAMPLE_PREPARATION_PROTOCOL.md"), warn = FALSE),
                      collapse = " ")
expect("pride", "the full sample-preparation protocol is published",
       all(vapply(c("5% DDM", "5 mM TCEP", "20 mM CAA", "0.1 M TEAB",
                    "95 degrees C", "4 ng Lys-C per sample", "6 ng trypsin per sample",
                    "TFA to 1%", "SpeedVac", "Evotip", "4.2 uL A buffer"),
                  function(x) grepl(x, protocol_txt, fixed = TRUE), logical(1))))
expect("pride", "remaining required SDRF gaps are exact and evidence-backed",
       setequal(fs$sdrf_field[fs$status == "MISSING_REQUIRED_METADATA" &
                                fs$sdrf_field %in% names(sdrf_v)],
                c("characteristics[developmental stage]", "comment[instrument]",
                  "comment[proteomics data acquisition method]")))
expect("pride", "externally held raw acquisition directories remain a separate blocker",
       any(fs$sdrf_field == "raw acquisition FILES (not an SDRF column)" &
             fs$status == "MISSING_REQUIRED_METADATA"))
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

# The acquisition-metadata gap must still be reported, and no readiness status may hide it.
expect("pride", "the acquisition-metadata gap is still reported as incomplete",
       pride_status == "PRIDE_METADATA_INCOMPLETE" || n_req_sdrf == 0L,
       pride_status)
expect("pride", "PRIDE readiness is not claimed while required SDRF fields are missing",
       !(pride_status %in% c("PRIDE_READY", "PRIDE_READY_PENDING_RAW_UPLOAD") &&
           n_req_sdrf > 0L))
expect("pride", "PRIDE_READY is not claimed while the raw acquisition files are missing",
       !(pride_status == "PRIDE_READY" && n_req_raw > 0L))

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
  # A scratch root is not named: it is a temporary location, never the published copy.
  if (release_is_shared_drive_root(OUT_ROOT, DATA_ROOT)) {
    paste0("Release root: `", OUT_ROOT, "`")
  } else "Release root: a scratch build outside the shared data root (not a published copy).",
  paste0("Repository commit: `", release_git_commit(REPO_ROOT), "`"),
  "",
  paste0("**Release validation: ", n_pass, " PASS, ", n_fail, " FAIL, ", n_skip, " SKIP**"),
  paste0("**PRIDE readiness: `", pride_status, "`** (",
         n_ready, " SDRF fields READY, ", n_req_sdrf, " MISSING_REQUIRED, ", n_opt,
         " MISSING_OPTIONAL)"),
  paste0("Raw acquisition files/directories are a required deposition input but not an SDRF",
         " column, so they are counted separately: ", n_req_raw,
         " outstanding required deposition blocker(s) of that kind."),
  "",
  "A PRIDE field that is unavailable does not fail the scientific release. It",
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
  stop(sprintf("Publication release validation failed: %d check(s).", n_fail),
       call. = FALSE)
}
cat("\nAll publication release checks hold.\n")
