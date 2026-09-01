#!/usr/bin/env Rscript

# Contract: the lineage is honest. Every step from acquisition to release is listed, the
# locked artefacts are hash-verified, unproven edges say UNRESOLVED, missing versions say
# UNKNOWN, and no vague version string survives.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("07_publication_release", "tests", "test_release_lineage.R")
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
lin_path <- file.path(OUT_ROOT, "provenance", "data_lineage.tsv")
sw_path <- file.path(OUT_ROOT, "provenance", "software_versions.tsv")
ap_path <- file.path(OUT_ROOT, "provenance", "analysis_parameters.tsv")
gap_path <- file.path(OUT_ROOT, "provenance", "UPSTREAM_PREPROCESSING_GAP.md")
if (!all(file.exists(c(lin_path, sw_path, ap_path, gap_path)))) {
  cat("No built publication release at", OUT_ROOT, "\n")
  cat("Skipping publication release lineage contracts.\n")
  quit(save = "no", status = 0L)
}

rd <- function(p) read.delim(p, sep = "\t", quote = "", stringsAsFactors = FALSE,
                             check.names = FALSE)
lineage <- rd(lin_path); sw <- rd(sw_path); ap <- rd(ap_path)

cat("=== the locked artefacts still hash as locked ===\n")
if (!file.exists(release_locked_artefact_path("animal_level_input_gct", DATA_ROOT))) {
  cat("  [skip]  shared data root unreachable; locked hashes not verified\n")
} else {
  locked <- release_verify_locked_artefacts(DATA_ROOT)
  for (i in seq_len(nrow(locked))) {
    expect(locked$matches[i],
           sprintf("%s SHA256 matches the locked contract", locked$artefact[i]))
  }
  expect(all(locked$expected_sha256 == locked$observed_sha256),
         "both locked GCT hashes are unchanged")
}

cat("\n=== the lineage covers the whole chain ===\n")
required <- c("A02_acquisition_raw", "A03_search_quantification",
              "A04_protein_group_matrix_prefilter", "A07_processed_matrix",
              "A08_animal_level_matrix", "A09_protigy_statistics",
              "A10_split_differential_tables", "A11_uniprot_mapping_reference",
              "A12_mapped_differential_tables", "A13_enrichment", "A14_ewce", "A15_pca",
              "A16_figure_panels", "A17_publication_release")
missing <- setdiff(required, lineage$artifact_id)
expect(length(missing) == 0L,
       sprintf("every lineage stage is present%s",
               if (length(missing)) paste0(" (missing: ", paste(missing, collapse = ", "), ")")
               else ""))
expect(all(nzchar(lineage$processing_stage)),
       "every artefact records how it was produced")
expect(all(nzchar(lineage$analysis_unit)), "every artefact records its analysis unit")
expect(all(lineage$canonical %in% c("yes", "no")),
       "every artefact is marked canonical or not")

cat("\n=== unproven edges are UNRESOLVED, not invented ===\n")
combined <- paste(lineage$processing_stage, lineage$notes)
expect(any(grepl("UNRESOLVED", combined)),
       "at least one lineage edge is honestly marked UNRESOLVED")
processed <- lineage[lineage$artifact_id == "A07_processed_matrix", , drop = FALSE]
expect(nrow(processed) == 1L && grepl("UNRESOLVED", paste(processed$processing_stage,
                                                          processed$notes)),
       "the step producing the 5,349-protein processed matrix is marked UNRESOLVED")
search_out <- lineage[lineage$artifact_id == "A03_search_quantification", , drop = FALSE]
expect(nrow(search_out) == 1L && grepl("UNRESOLVED", search_out$notes),
       "the absent search-software output is marked UNRESOLVED")
raw <- lineage[lineage$artifact_id == "A02_acquisition_raw", , drop = FALSE]
expect(nrow(raw) == 1L && raw$software_version == "UNKNOWN",
       "the acquisition instrument version is UNKNOWN, not guessed")

cat("\n=== the locked hashes appear in the lineage ===\n")
expect(any(lineage$sha256 == RELEASE_LOCKED_ARTEFACTS$animal_level_input_gct$sha256,
           na.rm = TRUE),
       "the animal-level GCT hash appears in the lineage")
expect(any(lineage$sha256 == RELEASE_LOCKED_ARTEFACTS$protigy_stat_gct$sha256, na.rm = TRUE),
       "the ProTigy statistics GCT hash appears in the lineage")

cat("\n=== software and database versions ===\n")
expect(nrow(sw) > 20L, sprintf("versions recorded for %d components", nrow(sw)))
expect(!any(grepl("[0-9]+[.][0-9]*[.]?x$", sw$version)),
       "no vague version string such as 1.1.x survives")
expect(all(nzchar(sw$version)), "no version field is blank")
unknown <- sw[sw$version == "UNKNOWN", , drop = FALSE]
expect(nrow(unknown) > 0L, "genuinely unknown versions are recorded as UNKNOWN")
expect(all(nzchar(unknown$recorded_by)),
       "every UNKNOWN version says where it would have to come from")
expect(any(grepl("^R$", sw$component)), "the R version is recorded")
for (pkg in c("clusterProfiler", "fgsea", "org.Mm.eg.db", "AnnotationDbi", "limma",
              "EWCE", "ewceData", "DOSE", "GOSemSim", "ggplot2", "digest")) {
  expect(any(sw$component == pkg & grepl("^[0-9]", sw$version)),
         sprintf("%s has an exact version", pkg))
}
expect(any(grepl("idmapping", sw$component, ignore.case = TRUE)),
       "the UniProt mapping reference is recorded as a versioned component")
idm <- sw[grepl("idmapping", sw$component, ignore.case = TRUE), , drop = FALSE]
expect(all(nzchar(idm$evidence_sha256) & !is.na(idm$evidence_sha256)),
       "the mapping reference carries a SHA256")
expect(any(grepl("ProTigy", sw$component)), "ProTigy is recorded")
expect(any(grepl("DIA-NN", sw$component)), "the search software is recorded")

cat("\n=== analysis parameters ===\n")
expect(nrow(ap) > 30L, sprintf("%d analysis parameters recorded", nrow(ap)))
expect(any(ap$parameter == "ranking_statistic" & ap$value == "t"),
       "the canonical ranking statistic is recorded as t")
expect(any(ap$parameter == "ora_universe"),
       "the ORA universe definition is recorded")
expect(any(grepl("seed", ap$parameter)), "seed state is recorded where it exists")
expect(any(ap$status == "MISSING_UNKNOWN"),
       "parameters that were never recorded are marked missing, not filled in")

cat("\n=== the upstream gap is documented ===\n")
gap <- paste(readLines(gap_path, warn = FALSE), collapse = " ")
expect(grepl("01_impute.r", gap, fixed = TRUE),
       "the imputation script is audited by name")
expect(grepl("set.seed", gap, fixed = TRUE),
       "the seed question is addressed explicitly")
expect(grepl("UNRESOLVED", gap, fixed = TRUE),
       "unrecoverable steps are named UNRESOLVED")
expect(grepl("umap_adjusted", gap, fixed = TRUE),
       "the umap_adjusted versus pcaAdjusted naming discrepancy is recorded")
expect(grepl("NOT", gap) && grepl("regenerated", gap),
       "the document states that the matrix was deliberately not regenerated")

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Publication release lineage contracts failed: %d", failures), call. = FALSE)
}
cat("All publication release lineage contracts hold.\n")
