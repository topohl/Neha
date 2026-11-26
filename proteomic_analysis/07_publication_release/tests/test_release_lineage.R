#!/usr/bin/env Rscript

# Checks: every step from acquisition to release is in data_lineage.tsv, the locked
# artefacts match their hashes, unproven steps say UNRESOLVED, missing versions say UNKNOWN,
# no vague version string is left, and repository files are named by paths relative to
# proteomic_analysis/.

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
  cat("Skipping publication release lineage checks.\n")
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
           sprintf("%s SHA256 matches the locked value", locked$artefact[i]))
  }
  expect(all(locked$expected_sha256 == locked$observed_sha256),
         "every locked artefact hash is unchanged")
}

cat("\n=== the lineage covers the whole chain ===\n")
required <- c("A00_experimenter_metadata", "A00_sample_preparation_protocol",
              "A02_acquisition_raw", "A03_search_quantification",
              "A04_protein_group_matrix_prefilter", "A07_processed_matrix",
              "A08_animal_level_matrix", "A09_differential_statistics",
              "A10_split_differential_tables", "A11_uniprot_mapping_reference",
              "A12_mapped_differential_tables", "A13_enrichment", "A14_ewce", "A15_pca",
              "A16_figure_panels", "A16c_figure_panels_within_compartment",
              "A17_publication_release")
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
em_path <- release_experimenter_metadata_source_path(repo_root)
protocol_path <- release_sample_preparation_source_path(repo_root)
expect(any(lineage$artifact_id == "A00_experimenter_metadata" &
             lineage$sha256 == release_sha256(em_path)),
       "experimenter metadata is represented by its checked-in source hash")
expect(any(lineage$artifact_id == "A00_sample_preparation_protocol" &
             lineage$sha256 == release_sha256(protocol_path)),
       "sample-preparation protocol is represented by its checked-in source hash")

cat("\n=== repository files are named by relative paths ===\n")
for (a in list(list(id = "A00_experimenter_metadata", path = em_path),
               list(id = "A00_sample_preparation_protocol", path = protocol_path))) {
  row <- lineage[lineage$artifact_id == a$id, , drop = FALSE]
  p <- if (nrow(row) == 1L) as.character(row$path_or_release_filename) else NA_character_
  expect(identical(p, release_repo_relative(a$path, repo_root)) &&
           identical(release_sha256(release_repo_resolve(p, repo_root)),
                     as.character(row$sha256)),
         sprintf("%s is named relative to proteomic_analysis/ and resolves to its file (%s)",
                 a$id, p))
}
ev <- unique(unlist(strsplit(ap$evidence_path[!is.na(ap$evidence_path)], ";", fixed = TRUE)))
ev_repo <- ev[release_is_repo_path(ev, repo_root)]
expect(length(ev_repo) > 0L && all(file.exists(release_repo_resolve(ev_repo, repo_root))),
       sprintf("repository evidence paths in analysis_parameters.tsv resolve (%d)",
               length(ev_repo)))
repo_prefix <- tolower(paste0(repo_root, "/"))
expect(!any(grepl(repo_prefix, tolower(c(lineage$path_or_release_filename, ap$evidence_path)),
                  fixed = TRUE)),
       "no lineage or parameter path names the local repository folder")
si_path <- file.path(OUT_ROOT, "provenance", "sessionInfo_release.txt")
si_repo <- if (file.exists(si_path)) {
  grep("^Repository:", readLines(si_path, warn = FALSE), value = TRUE)
} else character(0)
expect(length(si_repo) == 1L && !grepl(tolower(repo_root), tolower(si_repo), fixed = TRUE),
       "sessionInfo_release.txt names the repository without a local path")

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
expect(any(lineage$sha256 == RELEASE_LOCKED_ARTEFACTS$stat_gct$sha256, na.rm = TRUE),
       "the limma statistics GCT hash appears in the lineage")

cat("\n=== no unsupported upstream-software identity claim survives ===\n")
# Phrases the PRIDE metadata recovery audit ruled unsupported: DIA-NN use itself is not
# proven from the retained project files, so a bare DIA-NN mention may only appear as an
# qualified recovery example.
unsupported_claims <- c(
  "DIA-NN column schema", "DIA-NN pg_matrix", "DIA-NN output establishes",
  "software identity evidenced by", "DIA-NN First.Protein.Description",
  "expected DIA-NN column", "DIA-NN report / library",
  "DIA-NN report / spectral library", "DIA-NN search output",
  "Quantification software | DIA-NN"
)
claim_files <- c(
  list.files(OUT_ROOT, pattern = "[.]md$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(OUT_ROOT, c("provenance", "pride", "metadata")),
             pattern = "[.]tsv$", full.names = TRUE)
)
claim_files <- unique(claim_files[file.exists(claim_files)])
expect(length(claim_files) > 0L,
       sprintf("release narrative artefacts are present to audit (%d files)",
               length(claim_files)))
for (f in claim_files) {
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  txt_lower <- tolower(txt)
  hits <- unsupported_claims[vapply(tolower(unsupported_claims),
                                    function(b) grepl(b, txt_lower, fixed = TRUE),
                                    logical(1))]
  expect(length(hits) == 0L,
         sprintf("%s makes no unsupported upstream-software identity claim%s",
                 basename(f),
                 if (length(hits)) paste0(" (found: ", paste(hits, collapse = "; "), ")")
                 else ""))
  if (grepl("DIA-?NN", txt, ignore.case = TRUE)) {
    expect(grepl("if DIA-NN was used", txt, fixed = TRUE),
           sprintf("%s mentions DIA-NN only as a qualified recovery example", basename(f)))
  }
}

cat("\n=== software and database versions ===\n")
expect(nrow(sw) > 20L, sprintf("versions recorded for %d components", nrow(sw)))
expect(!any(grepl("[0-9]+[.][0-9]*[.]?x$", sw$version)),
       "no vague version string such as 1.1.x survives")
expect(all(nzchar(sw$version)), "no version field is blank")
unknown <- sw[sw$version == "UNKNOWN", , drop = FALSE]
expect(nrow(unknown) > 0L, "unknown versions are recorded as UNKNOWN")
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
expect(any(grepl("ProTigy", sw$component)), "ProTigy is recorded (superseded runs)")
# The upstream search/quantification software cannot be identified from the retained files
# (no run log, configuration, report table or spectral library; only the `pg.matrix` file
# naming). It must be listed with version UNKNOWN, and no software or lineage row may name
# a specific search engine (the test looks for DIA-NN).
usq <- sw[grepl("search / quantification software", sw$component, fixed = TRUE), ,
          drop = FALSE]
expect(nrow(usq) == 1L,
       "exactly one row covers the upstream search / quantification software")
if (nrow(usq) == 1L) {
  expect(identical(as.character(usq$version), "UNKNOWN"),
         sprintf("the upstream search/quantification software version is UNKNOWN (got %s)",
                 as.character(usq$version)))
  expect(identical(as.character(usq$status), "MISSING_RECOVERABLE"),
         sprintf("that software row is marked MISSING_RECOVERABLE (got %s)",
                 as.character(usq$status)))
  expect(nzchar(usq$recorded_by) && !grepl("DIA-?NN", usq$component, ignore.case = TRUE),
         "the software row names no specific upstream search engine")
}
expect(!any(grepl("DIA-?NN", sw$component, ignore.case = TRUE)),
       "no software component asserts a specific upstream search-software identity")
expect(!any(grepl("DIA-?NN", lineage$software, ignore.case = TRUE)),
       "no lineage row attributes an artefact to a specific upstream search software")

# Since 2026-11-24 the canonical statistics are limma, one model per sample class; its
# versions are read from that run's manifest. The ProTigy version of the superseded pooled
# animal-level run was recovered on 2026-09-02 by cross-source audit (installed 2.4.1 built
# 2026-08-24 12:48 UTC; that run's v2-only YAML parameter export written 13:06 UTC the same
# day). These checks cover both and keep ProTigy scoped to the superseded run.
source(file.path(LAYER, "R", "release_software_versions.R"))
stat_run <- tryCatch(release_read_stat_run(DATA_ROOT), error = function(e) NULL)
canon <- sw[sw$applies_to_stage == "canonical animal-level differential statistics", ,
            drop = FALSE]
expect(any(canon$component == "limma"),
       "the canonical differential statistics are attributed to limma")
expect(!any(grepl("ProTigy", canon$component)),
       "no ProTigy row is scoped to the canonical statistics")
if (!is.null(stat_run)) {
  expect(identical(as.character(canon$version[canon$component == "limma"]),
                   stat_run$limma_version),
         sprintf("the limma version equals the run manifest (%s)", stat_run$limma_version))
  expect(identical(as.character(canon$version[canon$component == "R"]), stat_run$r_version),
         sprintf("the R version equals the run manifest (%s)", stat_run$r_version))
}
prot_animal <- sw[grepl("pooled animal-level statistical GCT", sw$component, fixed = TRUE), ,
                  drop = FALSE]
expect(nrow(prot_animal) == 1L, "exactly one ProTigy row covers the pooled animal-level run")
if (nrow(prot_animal) == 1L) {
  expect(identical(as.character(prot_animal$version), "2.4.1"),
         sprintf("the pooled animal-level ProTigy version is the recovered 2.4.1 (got %s)",
                 as.character(prot_animal$version)))
  expect(identical(as.character(prot_animal$status), "KNOWN_VERIFIED"),
         "the recovered pooled-run ProTigy version is marked KNOWN_VERIFIED")
  expect(grepl("superseded", prot_animal$applies_to_stage, ignore.case = TRUE),
         "the pooled ProTigy run is scoped as superseded")
  expect(!grepl("1.1.8", as.character(prot_animal$version), fixed = TRUE),
         "the pooled animal-level run does NOT claim the 2025 hemisphere-level version")
  expect(grepl("2026-08-24", as.character(prot_animal$notes)) &&
           grepl("DISPROVEN|disproven", as.character(prot_animal$notes)),
         "the notes carry the build/export date evidence and record v1.1.8 as disproven")
  expect(nzchar(as.character(prot_animal$evidence_sha256)) &&
           !is.na(prot_animal$evidence_sha256) &&
           as.character(prot_animal$evidence_sha256) != "NONE",
         "the recovered version cites a hashed evidence artefact")
}
prot_2025 <- sw[grepl("hemisphere-level runs", sw$component, fixed = TRUE), , drop = FALSE]
if (nrow(prot_2025) == 1L) {
  expect(as.character(prot_2025$version) %in% c("1.1.8", "v1.1.8"),
         "the 2025 hemisphere-level ProTigy rows still report v1.1.8")
  expect(!identical(as.character(prot_2025$version),
                    as.character(prot_animal$version)),
         "the 2025 and animal-level ProTigy versions are recorded as different")
}
# The data_lineage.tsv row for the statistics step must give the run's limma version and
# the locked statistics GCT hash.
stat_step <- lineage[lineage$artifact_id == "A09_differential_statistics", , drop = FALSE]
expect(nrow(stat_step) == 1L, "exactly one lineage row is the canonical statistics step")
if (nrow(stat_step) == 1L && !is.null(stat_run)) {
  expect(grepl(paste0("limma ", stat_run$limma_version), stat_step$software_version,
               fixed = TRUE),
         "the lineage statistics step carries the run-recorded limma version")
  expect(identical(as.character(stat_step$sha256), RELEASE_LOCKED_ARTEFACTS$stat_gct$sha256),
         "the lineage statistics step carries the locked statistics GCT hash")
}
pooled_step <- lineage[lineage$artifact_id == "A09p_protigy_statistics_pooled", , drop = FALSE]
expect(nrow(pooled_step) == 1L &&
         identical(as.character(pooled_step$software_version), "2.4.1") &&
         pooled_step$canonical == "no",
       "the superseded pooled ProTigy step is non-canonical and carries its recovered version")

cat("\n=== analysis parameters ===\n")
expect(nrow(ap) > 30L, sprintf("%d analysis parameters recorded", nrow(ap)))
expect(any(ap$parameter == "ranking_statistic" & ap$value == "t"),
       "the canonical ranking statistic is recorded as t")
expect(any(ap$parameter == "ora_universe"),
       "the ORA universe definition is recorded")
expect(any(grepl("seed", ap$parameter)), "seed state is recorded where it exists")
expect(any(ap$parameter == "dose_seed_argument" & ap$value == "FALSE"),
       "the enrichment seed handling (DOSE seed = FALSE) is a documented parameter")
expect(any(ap$parameter == "df_residual" & grepl("^differential statistics", ap$stage)) &&
         any(ap$parameter == "ebayes_trend" & ap$value == "FALSE"),
       "the limma model parameters (residual df, robust, no trend) are recorded")
expect(any(ap$status == "MISSING_UNKNOWN"),
       "parameters that were never recorded are marked missing, not filled in")

cat("\n=== the upstream gap is documented ===\n")
gap <- paste(readLines(gap_path, warn = FALSE), collapse = " ")
expect(grepl("01_impute.r", gap, fixed = TRUE),
       "the imputation script is audited by name")
expect(grepl("set.seed", gap, fixed = TRUE),
       "the seed question is addressed")
expect(grepl("UNRESOLVED", gap, fixed = TRUE),
       "unrecoverable steps are named UNRESOLVED")
expect(grepl("umap_adjusted", gap, fixed = TRUE),
       "the umap_adjusted versus pcaAdjusted naming discrepancy is recorded")
expect(grepl("NOT", gap) && grepl("regenerated", gap),
       "the document states that the matrix was not regenerated")

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Publication release lineage checks failed: %d", failures), call. = FALSE)
}
cat("All publication release lineage checks hold.\n")
