#!/usr/bin/env Rscript

# Publication release, stage 10: lineage, software and database versions, parameters.
#
# Produces
#   provenance/data_lineage.tsv                 raw -> processed -> analysis, step by step
#   provenance/UPSTREAM_PREPROCESSING_GAP.md    what cannot be reproduced, and why
#   provenance/software_versions.tsv            versions, or UNKNOWN with a source
#   provenance/sessionInfo_release.txt          the R session that built this package
#   provenance/analysis_parameters.tsv          every recorded analysis parameter
#
# Steps that cannot be shown from the files are marked UNRESOLVED. The two that matter most:
# the search-software output is not in the project, and the code that made the 5,349-protein
# processed matrix is not in the repository.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_validation.R"))
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_software_versions.R"))
release_require("digest")

DATA_ROOT <- release_data_root()
PROJECT_ROOT <- release_project_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/10_build_provenance.R"
UNRESOLVED <- "UNRESOLVED"
UNKNOWN <- "UNKNOWN"

release_banner("stage 10 -- provenance")

ANIMAL_ROOT <- file.path(DATA_ROOT, "02_data", "animal_level")
# The primary statistics run (per sample class, 2026-11-24) and its enrichment runs.
STAT_ROOT <- release_stats_root(DATA_ROOT)
ENRICH_RUN_ROOT <- release_rel_path(RELEASE_STATISTICS_RUN$enrichment_relative, DATA_ROOT)
ENRICH_ROOT <- release_enrichment_root("GO_BP", DATA_ROOT)      # GO-BP
ENRICH_ROOT_CC <- release_enrichment_root("GO_CC", DATA_ROOT)   # GO-CC
KEGG_ROOT <- release_kegg_root(DATA_ROOT)                       # KEGG (20261124 GO_BP folder)
EWCE_ROOT <- file.path(DATA_ROOT, "03_output", "ewce",
                       "EWCE_Results_animal_level_validation_20260825")
PCA_ROOT <- file.path(DATA_ROOT, "03_output", "pca",
                      "pca_plots_animal_level_validation_20260825_rerun")
REVISION_ROOT <- file.path(DATA_ROOT, "03_output", "reviewer_revision_animal_level_20260827")
REBUILD_ROOT <- release_figure_rebuild_root(DATA_ROOT)
experimenter_source <- release_experimenter_metadata_source_path(REPO_ROOT)
protocol_source <- release_sample_preparation_source_path(REPO_ROOT)
experimenter_metadata <- release_read_experimenter_metadata(REPO_ROOT)
sample_preparation_protocol <- release_read_sample_preparation_protocol(REPO_ROOT)

hash_or <- function(path, fallback = NA_character_) {
  if (!is.na(path) && nzchar(path) && file.exists(path) && !dir.exists(path)) {
    release_sha256(path)
  } else fallback
}

mapped_index <- release_read_csv(file.path(STAT_ROOT, "mapped", "indexMappedComparisons.csv"))
pkg_versions <- release_read_csv(file.path(ENRICH_ROOT, "audits",
                                           "package_database_versions.csv"))
run_parameters <- release_read_csv(file.path(ENRICH_ROOT, "audits", "run_parameters.csv"))
run_parameters_cc <- release_read_csv(file.path(ENRICH_ROOT_CC, "audits", "run_parameters.csv"))
run_parameters_kegg <- release_read_csv(file.path(KEGG_ROOT, "audits", "run_parameters.csv"))
permutations <- release_enrichment_permutations(DATA_ROOT)
stat_run <- release_read_stat_run(DATA_ROOT)
# EWCE's own differential limma fit, as its audit recorded it.
ewce_model <- release_ewce_differential_model(DATA_ROOT)
# Evidence for the enrichment seed handling: the robustness branch reran the mCherry learning
# GSEA with the recorded seeds and DOSE seed = FALSE and compared it with the stored run.
p_seed_checks <- file.path(release_inferential_root(DATA_ROOT), "robustness", "checks.csv")
seed_checks <- release_read_csv(p_seed_checks)
seed_check_rows <- vapply(c("GO-BP", "GO-CC"), function(ont) {
  hit <- grepl(paste0(ont, " GSEA vs stored enrichment"), seed_checks$check, fixed = TRUE)
  if (sum(hit) != 1L || !identical(toupper(as.character(seed_checks$passed[hit])), "TRUE")) {
    stop("robustness/checks.csv lacks a passing ", ont, " rerun-versus-stored check: ",
         p_seed_checks, call. = FALSE)
  }
  which(hit)
}, integer(1))
# The test that ties 01_preprocessing/02b_limma_stat_gct.r to the locked statistics GCT.
p_limma_test <- file.path(REPO_ROOT, "tests", "test_limma_stat_gct.R")
# AveExpr as the per-sample-class model reports it, measured on the released table.
ave <- release_aveexpr_summary(read.delim(
  gzfile(release_path("differential_analysis", "primary_differential_proteins.tsv.gz",
                      create_dir = FALSE)),
  sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE))
protigy_params_2025 <- file.path(DATA_ROOT, "01_input", "raw_proteomics",
                                 "20251107_pg.matrix_Neha", "params.txt")
animal_param_yaml <- file.path(ANIMAL_ROOT, "neha_proteome_parameters.yaml")
ewce_session_path <- file.path(EWCE_ROOT, "03_QC_Mapping_Logs",
                               "reproducibility_session_info.txt")
pca_session_path <- file.path(PCA_ROOT, "tables", "meta", "sessionInfo.txt")

# --------------------------------------------------------------------------------------
# ProTigy version of the superseded pooled animal-level run (2026-08-24)
# --------------------------------------------------------------------------------------
# The version, its evidence and this check are in R/release_software_versions.R; the check
# runs here too so this stage stops on its own if the parameter export changes.
release_assert_protigy_v2_export_signature(animal_param_yaml)

# --------------------------------------------------------------------------------------
# lineage
# --------------------------------------------------------------------------------------

lin <- function(artifact_id, artifact_type, path_or_release_filename, parent_artifact,
                processing_stage, software, software_version, analysis_unit,
                sha256 = NA_character_, canonical = "yes", notes = NA_character_) {
  data.frame(artifact_id = artifact_id, artifact_type = artifact_type,
             path_or_release_filename = release_repo_relative(path_or_release_filename,
                                                              REPO_ROOT),
             parent_artifact = parent_artifact, processing_stage = processing_stage,
             software = software, software_version = software_version,
             analysis_unit = analysis_unit, sha256 = sha256, canonical = canonical,
             notes = notes, stringsAsFactors = FALSE, check.names = FALSE)
}

p_pg_raw <- file.path(PROJECT_ROOT, "pg.matrix_raw.txt")
p_sample_annotation <- file.path(PROJECT_ROOT, "sample_annotation.xlsx")
p_sample_info <- file.path(DATA_ROOT, "01_input", "metadata", "sample_info.xlsx")
p_processed <- file.path(DATA_ROOT, "02_data", "gct",
                         "pg.matrix_filtered_pcaAdjusted_unnormalized.gct")
p_processed_dup <- file.path(DATA_ROOT, "01_input", "raw_proteomics",
                             "20251107_pg.matrix_Neha", "pg.matrix_filtered_pcaAdjusted.gct")
p_animal_gct <- release_locked_artefact_path("animal_level_input_gct", DATA_ROOT)
p_stat_gct <- release_locked_artefact_path("stat_gct", DATA_ROOT)
if (!identical(hash_or(p_stat_gct), stat_run$output_gct_sha256)) {
  stop("The statistics GCT does not hash to the value its run manifest recorded.",
       call. = FALSE)
}
p_stat_manifest <- release_locked_artefact_path("stat_run_manifest", DATA_ROOT)
p_stat_gct_pooled <- file.path(release_rel_path(RELEASE_STATISTICS_RUN$pooled_stats_relative,
                                                DATA_ROOT),
                               RELEASE_STATISTICS_RUN$pooled_stat_gct)
p_pooled_mapped_index <- file.path(release_rel_path(RELEASE_STATISTICS_RUN$pooled_stats_relative,
                                                    DATA_ROOT),
                                   "mapped", "indexMappedComparisons.csv")
p_pooled_enrich_index <- file.path(
  release_rel_path(RELEASE_STATISTICS_RUN$pooled_enrichment_relative, DATA_ROOT),
  "indexEnrichmentComparisons.csv")
if (!identical(hash_or(p_stat_gct_pooled), RELEASE_STATISTICS_RUN$pooled_stat_gct_sha256)) {
  stop("The superseded pooled statistics GCT no longer hashes to its recorded value: ",
       p_stat_gct_pooled, call. = FALSE)
}
p_limma_script <- file.path(REPO_ROOT, RELEASE_PRIMARY_MODEL$script)
p_split_index <- file.path(STAT_ROOT, "split", "indexComparisons.csv")
p_mapped_index <- file.path(STAT_ROOT, "mapped", "indexMappedComparisons.csv")
p_idmapping <- file.path(DATA_ROOT, "01_input", "references", "MOUSE_10090_idmapping.dat")
p_manual_mapping <- file.path(DATA_ROOT, "01_input", "references", "manual_mapping.xlsx")
p_enrich_index <- file.path(ENRICH_ROOT, "indexEnrichmentComparisons.csv")
p_enrich_index_cc <- file.path(ENRICH_ROOT_CC, "indexEnrichmentComparisons.csv")
p_enrich_index_kegg <- file.path(KEGG_ROOT, "indexEnrichmentComparisons.csv")
p_ewce_table <- file.path(EWCE_ROOT, "02_Tables_Supplements", "Supplementary_Table_EWCE.xlsx")
p_pca_variance <- file.path(PCA_ROOT, "tables", "variance", "pca_variance_explained.csv")
p_loom <- file.path(DATA_ROOT, "01_input", "single_cell", "loom", "l1_amygdala.loom")
p_protein_count <- file.path(PROJECT_ROOT, "protein_count.xlsx")

processed_sha <- hash_or(p_processed)
processed_dup_sha <- hash_or(p_processed_dup)
duplicate_note <- if (!is.na(processed_sha) && identical(processed_sha, processed_dup_sha)) {
  paste("Byte-identical to 02_data/gct/pg.matrix_filtered_pcaAdjusted_unnormalized.gct",
        "(same SHA256). The copy sits inside the 2025-11-07 ProTigy session folder as that",
        "session's INPUT, not as its output, so the matrix existed by 2025-11-07 and was",
        "NOT produced by ProTigy. The producing step remains", UNRESOLVED)
} else {
  paste("NOT byte-identical to 02_data/gct/pg.matrix_filtered_pcaAdjusted_unnormalized.gct;",
        "the relationship between the two is", UNRESOLVED)
}

# Candidate numeric ancestor, verified rather than asserted: a Feb-2025 export at the
# project root whose filename encodes the processing steps the numbers actually show.
p_perseus <- file.path(PROJECT_ROOT,
                       "pg.matrix_filtered_70percent-onegroup_imputed_ANOVA_z-scored.txt")

lineage <- rbind(
  lin("A00_experimenter_metadata", "study / cohort / sample-class metadata",
      experimenter_source, NA_character_,
      "experimenter-supplied facts recorded in a release metadata file in the repository",
      "manual curation", "n/a", "study / cohort / sample class",
      hash_or(experimenter_source), "yes",
      paste("Defines CeM origin, marker-defined LCM classes, mouse strain/source, cohort-only",
            "sex and age-at-start facts, label-free status and explicit unresolved fields.",
            "No AnimalID-level sex was invented.")),
  lin("A00_sample_preparation_protocol", "sample-preparation protocol",
      protocol_source, "A00_experimenter_metadata",
      "experimenter-supplied protocol version 2024-11-28 codified as seven ordered steps",
      "manual curation", "n/a", "proteomics sample",
      hash_or(protocol_source), "yes",
      paste("Records 5% DDM / 5 mM TCEP / 20 mM CAA / 0.1 M TEAB lysis, 95 degrees C",
            "for 1 hour, sequential Lys-C and trypsin digestion, TFA stop or recorded",
            "SpeedVac option, Evotip cleanup and 4.2 uL A-buffer resuspension. CAA/TCEP",
            "do not establish search modification parameters.")),
  lin("A01_biological_samples", "biological material",
      "96 marker-defined LCM samples (12 animals x 4 sample classes x 2 hemispheres)",
      NA_character_, "dissection and sample preparation", UNKNOWN, UNKNOWN,
      "hemisphere sample", NA_character_, "yes",
      paste("Experimenter metadata establishes central medial amygdala (CeM), four",
            "marker-defined laser-capture microdissection classes, label-free proteomics",
            "and the versioned sample-preparation protocol. Instrument model, acquisition",
            "mode, search settings, AnimalID-level sex and age at tissue collection remain",
            UNRESOLVED)),
  lin("A02_acquisition_raw", "native acquisition data",
      "96 `.d` acquisition directories, e.g. Olive_20241217_..._8782.d",
      "A01_biological_samples", "LC-MS acquisition (2024-12-17)",
      paste("instrument model", UNKNOWN, "-- run-name alias `Olive`"), UNKNOWN,
      "acquisition run", NA_character_, "yes",
      paste("Filenames are PROVEN for all 96 from sample_annotation.xlsx (original path",
            "D:\\Proteomics\\Fabian\\Tobias\\<run>.d). The FILES are not present anywhere",
            "in the project tree (verified: zero .d/.raw/.wiff/.mzML). Presence:",
            UNRESOLVED)),
  lin("A03_search_quantification", "search / quantification output",
      "run log / report table / spectral library (software UNKNOWN)",
      "A02_acquisition_raw", "peptide and protein identification and quantification",
      "upstream search / quantification software", UNKNOWN, "acquisition run",
      NA_character_, "yes",
      paste("ABSENT from the project tree:", UNRESOLVED, "-- no run log, report table,",
            "spectral library or search configuration was found (no report.tsv,",
            "report.parquet, report.log.txt or .speclib). The retained processed files use",
            "the historical `pg.matrix` naming convention, and sample_annotation.xlsx",
            "carries a sheet named `report.stats_`, but the exact upstream",
            "search/quantification software and configuration could not be recovered from",
            "the retained project files. Software identity and version are both", UNKNOWN)),
  lin("A04_protein_group_matrix_prefilter", "protein-group matrix",
      p_pg_raw, "A03_search_quantification", "protein-group quantification export",
      "upstream search / quantification software", UNKNOWN, "acquisition run",
      hash_or(p_pg_raw), "yes",
      paste("5,747 protein groups x 96 quantitative columns + 7 annotation columns.",
            "This is the file supplied to the journal as",
            "`processed_protein_group_matrix_raw.txt`. It is NOT native raw MS data.",
            "Annotation columns are republished as",
            "processed_data/protein_feature_annotation.tsv.gz.")),
  lin("A05_sample_annotation", "sample annotation",
      p_sample_annotation, "A02_acquisition_raw", "acquisition-to-sample mapping",
      "manual / spreadsheet", UNKNOWN, "acquisition run", hash_or(p_sample_annotation), "yes",
      paste("The only artefact carrying original `.d` acquisition paths. Sheet is named",
            "`report.stats_`. Supplies the explicit _left/_right hemisphere labels.")),
  lin("A06_sample_info", "sample metadata",
      p_sample_info, "A05_sample_annotation", "legacy 2024 metadata schema",
      "manual / spreadsheet", UNKNOWN, "acquisition run", hash_or(p_sample_info), "yes",
      paste("Columns: id, sampleNumber, shortname, plate, group2, AnimalID,",
            "ReplicateGroup, celltype, ExpGroup, celltype_ExpGroup. NOTE: `group2` is",
            "celltype + ReplicateGroup (a HEMISPHERE code), not celltype + condition;",
            "`celltype_ExpGroup` is the condition-coded one.")),
  lin("A06b_sample_class_correction", "sample-class correction record",
      paste(c("clusterProfiler", RELEASE_SAMPLE_CLASS_CORRECTION$reference_relative),
            collapse = "/"),
      "A05_sample_annotation;A06_sample_info",
      paste0("historical sample-identity QC: sample-class metadata reassigned for ",
             RELEASE_N_SAMPLE_CLASS_CORRECTED, " of 96 acquisitions on ",
             RELEASE_SAMPLE_CLASS_CORRECTION$correction_date),
      "UMAP-based sample-identity QC", UNKNOWN,
      "hemisphere sample",
      release_verify_sample_class_correction_reference(DATA_ROOT)$observed_sha256, "yes",
      paste0("METADATA ONLY. Left hemispheres of C46 and C47, cyclically reassigned ",
             "neuropil -> mcherry -> neuron -> neuropil. Acquisition identifiers and ",
             "quantitative abundance values were UNCHANGED (max absolute numeric ",
             "difference 0). Status ", RELEASE_SAMPLE_CLASS_CORRECTION$status,
             ". No surviving prose rationale; this table is the preserved record. Both ",
             "labels are published in metadata/sample_metadata.tsv and ",
             "metadata/sample_class_corrections.tsv.")),
  lin("A07_processed_matrix", "processed protein abundance matrix",
      p_processed, "A04_protein_group_matrix_prefilter",
      paste("filtering, imputation and per-protein standardisation (measured on the values);",
            "log2 transformation, median centring and PCA adjustment are described for this",
            "file but not recorded -- PRODUCING CODE NOT IN REPOSITORY:", UNRESOLVED),
      UNKNOWN, UNKNOWN, "acquisition run", processed_sha, "yes",
      paste("5,349 x 96. The reproducible pipeline STARTS here. See",
            "provenance/UPSTREAM_PREPROCESSING_GAP.md. 01_preprocessing/01_impute.r does",
            "NOT produce this file: it writes per-sample-class .xlsx files, performs no PCA",
            "adjustment, and its input pg.matrix_raw.tsv is absent.")),
  lin("A07b_processed_matrix_duplicate", "processed protein abundance matrix (duplicate)",
      p_processed_dup, "A07_processed_matrix", "no transformation (identical copy)",
      "n/a", "n/a", "acquisition run", processed_dup_sha, "no", duplicate_note),
  lin("A07c_candidate_numeric_ancestor", "candidate numeric ancestor",
      p_perseus, "A04_protein_group_matrix_prefilter",
      paste("filename encodes: 70% one-group valid-value filter, imputation, ANOVA,",
            "z-scoring --", UNRESOLVED, "as a processing record"),
      UNKNOWN, UNKNOWN, "acquisition run", hash_or(p_perseus), "no",
      paste("5,349 x 106 (96 samples + 10 annotation/ANOVA columns), dated Feb 2025.",
            "VERIFIED: identical protein identifiers in identical order, the same 96 sample",
            "columns, and values agreeing with A07 to a maximum absolute difference of",
            "0.005 -- the signature of a 2-decimal serialisation, not of a further",
            "transformation. This makes it a strong candidate ancestor and it is the only",
            "artefact whose name is consistent with what A07's numbers actually show",
            "(per-protein z-scoring). A filename is not a processing record, so the edge",
            "is recorded as a candidate and remains", UNRESOLVED)),
  lin("A08_animal_level_matrix", "animal-level abundance matrix",
      p_animal_gct, "A07_processed_matrix",
      "equal-weight mean of Left and Right within AnimalID x sample_class",
      "R + 01_preprocessing/02a_prepare_animal_level_protigy_input.r",
      paste("R", R.version.string), "animal", hash_or(p_animal_gct), "yes",
      paste("5,349 x 48. LOCKED artefact. No transformation, normalisation, filtering,",
            "imputation or remapping is applied at this step.")),
  lin("A09_differential_statistics", "differential statistics GCT",
      p_stat_gct, "A08_animal_level_matrix", RELEASE_PRIMARY_MODEL$label,
      paste("R + limma +", RELEASE_PRIMARY_MODEL$script),
      paste0("limma ", stat_run$limma_version, "; R ", stat_run$r_version), "animal",
      hash_or(p_stat_gct), "yes",
      paste0("LOCKED artefact, written ", stat_run$run_timestamp_utc, " in the ProTigy GCT",
             " layout so the split, identifier-mapping and enrichment stages run unchanged.",
             " Versions, residual df, variance priors and input/output SHA256 are recorded",
             " per contrast in ", basename(p_stat_manifest), " (also LOCKED). Primary since ",
             RELEASE_PRIMARY_MODEL$date, "; replaces the pooled ProTigy run",
             " (A09p_protigy_statistics_pooled), whose one ",
             RELEASE_DESIGN_INVARIANTS$n_animal_level_units,
             "-column model shared a variance estimate across the four sample classes. logFC",
             " equals the pooled run's to floating-point rounding; AveExpr is now the",
             " sample-class mean; t, P, adjusted P and B differ.")),
  lin("A09p_protigy_statistics_pooled", "differential statistics GCT (superseded)",
      p_stat_gct_pooled, "A08_animal_level_matrix", RELEASE_PRIMARY_MODEL$superseded_label,
      "ProTigy", PROTIGY_ANIMAL_LEVEL_VERSION, "animal", hash_or(p_stat_gct_pooled), "no",
      paste("SUPERSEDED on", RELEASE_PRIMARY_MODEL$date, "as the primary statistics;",
            "retained for comparison and read only for the comparison tables in changelog",
            "section 7b and manuscript actions section 4b, from its mapped tables (index",
            paste0("02_data/animal_level/mapped/indexMappedComparisons.csv, SHA256 ",
                   hash_or(p_pooled_mapped_index), ")"),
            "and its enrichment run (index",
            paste0(paste(RELEASE_STATISTICS_RUN$pooled_enrichment_relative, collapse = "/"),
                   "/indexEnrichmentComparisons.csv, SHA256 ", hash_or(p_pooled_enrich_index),
                   "),"),
            "both of which record this GCT's hash. Parameters in",
            "02_data/animal_level/neha_proteome_parameters.yaml. ProTigy version for this",
            "(2026-08-24) run is", PROTIGY_ANIMAL_LEVEL_VERSION, "-- recovered 2026-09-02 by",
            "cross-source audit and NOT carried over from the 2025 runs, which used v1.1.5",
            "(2025-03-28) and v1.1.8 (2025-11-07, 2025-12-12). v1.1.8 is disproven for this",
            "run: it emits a params.txt header, this run emitted a v2 YAML export.")),
  lin("A10_split_differential_tables", "differential result tables",
      file.path(STAT_ROOT, "split", "forward"), "A09_differential_statistics",
      "split into 12 forward and 12 reverse per-contrast tables",
      "R + 01_preprocessing/03_gct_extractR.r", paste("R", R.version.string), "animal",
      hash_or(p_split_index), "yes",
      paste("5,349 rows per contrast. Index: indexComparisons.csv; run_contract_manifest.csv",
            "records the source GCT SHA256.")),
  lin("A11_uniprot_mapping_reference", "reference database",
      p_idmapping, NA_character_, "external download",
      "UniProt idmapping (MOUSE_10090)",
      as.character(mapped_index$mapping_reference_version[[1]]), "n/a",
      as.character(mapped_index$mapping_reference_sha256[[1]]), "yes",
      paste0("Snapshot date ", mapped_index$mapping_reference_snapshot_date_utc[[1]],
             "; ", mapped_index$mapping_reference_bytes[[1]], " bytes. The file itself does",
             " not encode a UniProt release number.")),
  lin("A11b_manual_mapping", "reference override",
      p_manual_mapping, NA_character_, "hand-curated identifier overrides",
      "manual / spreadsheet", UNKNOWN, "n/a",
      as.character(mapped_index$manual_mapping_sha256[[1]]), "yes",
      paste0(mapped_index$manual_mapping_rows[[1]], " override rows.")),
  lin("A12_mapped_differential_tables", "mapped differential result tables",
      file.path(STAT_ROOT, "mapped", "forward"),
      "A10_split_differential_tables;A11_uniprot_mapping_reference;A11b_manual_mapping",
      "UniProt identifier resolution", "R + 02_id_mapping/01_MapThatProt_batch.r",
      paste("R", R.version.string), "animal", hash_or(p_mapped_index), "yes",
      paste("5,327 of 5,349 mapped; 22 non-mouse contaminant identifiers retained in",
            "mapped/unmapped/. Index: indexMappedComparisons.csv.")),
  lin("A13_enrichment", "GSEA and ORA results",
      ENRICH_RUN_ROOT, "A12_mapped_differential_tables",
      "GSEA ranked by the moderated t (GO-BP, GO-CC, KEGG); ORA over FDR-significant lists",
      "R + clusterProfiler + fgsea",
      paste0("clusterProfiler ",
             pkg_versions$version[pkg_versions$component == "clusterProfiler"],
             "; fgsea ", pkg_versions$version[pkg_versions$component == "fgsea"]),
      "animal", hash_or(p_enrich_index), "yes",
      paste("Runs of 04_differential_expression_enrichment/01_clusterProfiler.r on the",
            "same mapped tables: GO_BP/ (GO-BP) and GO_CC/ (GO-CC) of this folder, and KEGG",
            "from the GO_BP/ folder of",
            paste0(permutations$kegg_relative, ", whose GO tables this run superseded."),
            permutations$sentence, RELEASE_GO_CC_PROVENANCE_NOTE,
            RELEASE_ENRICHMENT_SEED_HANDLING, RELEASE_KEGG_VERSION_NOTE,
            "The effect-size-ranked GSEA in the same runs is a",
            "SENSITIVITY analysis only; its ranking statistic is the standardised-abundance",
            "effect size recorded upstream as log2fc, not a log2 fold change. Indexes:",
            "GO_BP/indexEnrichmentComparisons.csv (hashed here),",
            paste0("GO_CC/indexEnrichmentComparisons.csv (SHA256 ", hash_or(p_enrich_index_cc),
                   ") and the KEGG run's GO_BP/indexEnrichmentComparisons.csv (SHA256 ",
                   hash_or(p_enrich_index_kegg), "); all record the locked statistics GCT",
                   " hash."))),
  lin("A14_ewce", "cell-type enrichment results",
      EWCE_ROOT, "A08_animal_level_matrix",
      "EWCE bootstrap enrichment, 10000 replicates",
      "R + EWCE + ewceData",
      paste0("EWCE 1.18.1; ewceData 1.18.0"), "animal", hash_or(p_ewce_table), "yes",
      paste("Cell-type reference is ewceData::ctd(), NOT the amygdala loom in 01_input.",
            "Primary settings: annotation level 2, top-N 250, FDR 0.05, measured-proteome",
            "background 4859. Run 2026-08-25 from the animal-level matrix alone:",
            ewce_model$sentence, "Its limma settings are recorded in",
            paste0(basename(ewce_model$path), "."))),
  lin("A14b_unused_single_cell_reference", "reference dataset (unused)",
      p_loom, NA_character_, "not used by any canonical analysis", "n/a", "n/a", "n/a",
      NA_character_, "no",
      paste("Present in 01_input/single_cell/ but the canonical EWCE run calls",
            "ewceData::ctd(). Recorded so its presence is not mistaken for evidence about",
            "the dissected brain region. Not hashed: 104 MB.")),
  lin("A15_pca", "PCA results",
      PCA_ROOT, "A08_animal_level_matrix", "principal component analysis",
      "R + prcomp", paste("R", R.version.string), "animal", hash_or(p_pca_variance), "yes",
      "48 animal-level units, 5349 proteins, 0 zero-variance rows removed."),
  lin("A16_figure_panels",
      "manuscript figure panels and source data (panels not using the differential statistics)",
      REVISION_ROOT, "A08_animal_level_matrix;A14_ewce;A15_pca",
      "figure regeneration at animal level (2026-08-27)",
      "R + 06_manuscript_figure_revision/*.R (frozen snapshots)",
      paste("R", R.version.string), "animal",
      hash_or(file.path(REVISION_ROOT, "manifests", "regenerated_panel_manifest.csv")), "yes",
      paste("Source of the panels that do not depend on the differential statistics (PCA,",
            "rank abundance, acquisition QC, EWCE, the secondary within-animal model). Its",
            "Fig. 2 D, F, G, H, I and Supp. Fig. 4 F were drawn from the superseded pooled",
            "statistics and are replaced by A16c. Per-panel source data copied into",
            "editor_source_data/figure_source_data/.")),
  lin("A16c_figure_panels_within_compartment",
      "manuscript figure panels and source data (per-sample-class statistics)",
      REBUILD_ROOT, "A12_mapped_differential_tables;A13_enrichment",
      paste("figure regeneration on the per-sample-class statistics",
            paste0("(", RELEASE_PRIMARY_MODEL$date, "; again 2026-11-25 on the GO enrichment",
                   " at nPermSimple ", format(permutations$go_n_perm_simple, big.mark = ","),
                   ")")),
      paste("R + 06_manuscript_figure_revision/within_compartment_20261125/ (copies of the",
            "rebuild's scripts/, run by run_rebuild.sh; the 2026-08-27 scripts re-pointed at",
            "this run, plus 16b_final_Fig3G_gsea_lollipop.R)"),
      paste("R", R.version.string), "animal", hash_or(release_figure_rebuild_record(DATA_ROOT)),
      "yes",
      paste("Fig. 2 D, F, G, H, I and Supp. Fig. 4 F (internal names Fig3*, SuppF). Fig. 2 G",
            "is now moderated-t-ranked GSEA of GO-BP and GO-CC terms instead of an ORA",
            "lollipop; the ORA results of the same comparison are quoted in",
            "figure_source_map.tsv. Per-panel source data copied into",
            "editor_source_data/figure_source_data/ under manuscript-numbered names and",
            "pinned by SHA256 in figure_source_map.tsv.")),
  lin("A16b_acquisition_qc_counts", "acquisition QC source data",
      p_protein_count, "A03_search_quantification",
      "protein and peptide identification counts per acquisition",
      UNKNOWN, UNKNOWN, "acquisition run", NA_character_, "yes",
      paste("Source of Supp. Fig. 4 A. Lives outside the clusterProfiler tree.",
            "Producing step is", UNRESOLVED, "-- not hashed here; the copied panel source",
            "data in editor_source_data/ carries its hash.")),
  lin("A17_publication_release", "publication / deposition package",
      "this release", paste("A00_experimenter_metadata;A00_sample_preparation_protocol;",
                            "A05_sample_annotation;A06_sample_info;A04_protein_group_matrix_prefilter;",
                            "A07_processed_matrix;A08_animal_level_matrix;A09_differential_statistics;",
                            "A09p_protigy_statistics_pooled;",
                            "A10_split_differential_tables;A12_mapped_differential_tables;",
                            "A13_enrichment;A14_ewce;A16_figure_panels;",
                            "A16c_figure_panels_within_compartment", sep = ""),
      "read-only repackaging; no statistic recomputed",
      "R + 07_publication_release/*.R", paste("R", R.version.string), "animal / acquisition",
      NA_character_, "yes",
      "Per-file hashes are in provenance/release_manifest.tsv and provenance/SHA256SUMS.txt.")
)

n_unresolved <- sum(grepl(UNRESOLVED, paste(lineage$processing_stage, lineage$notes)))
release_log("  lineage: ", nrow(lineage), " artefacts, ", n_unresolved,
            " carrying an UNRESOLVED edge")

# --------------------------------------------------------------------------------------
# upstream preprocessing gap
# --------------------------------------------------------------------------------------

impute_script <- file.path(REPO_ROOT, "01_preprocessing", "01_impute.r")
impute_src <- readLines(impute_script, warn = FALSE)
has_set_seed <- any(grepl("set\\.seed", impute_src))
has_rnorm <- any(grepl("rnorm\\s*\\(", impute_src))
impute_input_default <- grep("pg\\.matrix_raw\\.tsv", impute_src, value = TRUE)
impute_input_path <- file.path(DATA_ROOT, "01_input", "raw_proteomics", "pg.matrix_raw.tsv")

gap <- c(
  "# Upstream preprocessing: what is known and what is not",
  "",
  paste0("Generated by `", STAGE, "`."),
  "",
  "The reproducible part of this project starts at the 5,349-protein processed matrix",
  "",
  paste0("    ", p_processed),
  paste0("    SHA256 ", processed_sha),
  "",
  "Everything downstream of that file is reproducible from this repository. Everything",
  "upstream of it is not, and this document says exactly where the boundary falls rather",
  "than leaving a reader to discover it.",
  "",
  "## What is established",
  "",
  "| Question | Answer | Evidence |",
  "|---|---|---|",
  paste0("| Direct parent input | `pg.matrix_raw.txt`, 5,747 protein groups x 96 quantitative",
         " columns + 7 annotation columns | The 5,349 analysed protein groups are a strict",
         " subset: all 5,349 `id` values match `T: Protein.Names` in `pg.matrix_raw.txt`",
         " (verified at build time by stage 03) |"),
  paste0("| Is that the file the journal received? | Yes | It has exactly the 5,747 rows x 96",
         " measurement columns described for `processed_protein_group_matrix_raw.txt` |"),
  paste0("| Upstream search / quantification software | UNKNOWN | The retained processed",
         " files use the historical `pg.matrix` naming convention (annotation fields",
         " `T: Protein.Group`, `T: Protein.Names`, `T: Genes`,",
         " `T: First.Protein.Description`; `sample_annotation.xlsx` carries a sheet named",
         " `report.stats_`), but the exact upstream search/quantification software and",
         " configuration could not be recovered from the retained project files |"),
  paste0("| Is the processed matrix duplicated? | Yes | `",
         basename(p_processed_dup), "` under `01_input/raw_proteomics/20251107_pg.matrix_Neha/`",
         " is byte-identical (SHA256 `", processed_dup_sha, "`), which dates the processed",
         " matrix to the 2025-11-07 ProTigy-era drop |"),
  paste0("| Filtering threshold | 70% maximum missingness, applied within sample class |",
         " Encoded in the sibling filename",
         " `pg.matrix_filtered_70percent-onegroup_imputed_ANOVA_z-scored.gct` and in the",
         " `max_missing: 70` field of `neha_proteome_parameters.yaml` (the pooled ProTigy",
         " run's parameter export). **The exact",
         " application to THIS file is not independently recorded.** |"),
  paste0("| A 2025 ProTigy session's settings | log transform none, normalisation none,",
         " filter none, NA filter 75, SD filter 10, two-sample moderated t, adj.p 0.05 |",
         " `", protigy_params_2025, "` (2025-11-07 session, ProTigy v1.1.8) |"),
  paste0("| ...but they are not the only ones | a second surviving session",
         " (2025-12-12, also ProTigy v1.1.8) records DIFFERENT settings on the SAME input:",
         " `norm.data Median` and `na.filt.val 100`. Two further hemisphere-level GCTs",
         " (2025-10-30, 2025-12-15) have no surviving params.txt at all. **No parameter set",
         " can therefore be attached to this matrix**, and none is claimed. |"),
  "",
  "## What is NOT established",
  "",
  "| Step | Status | Detail |",
  "|---|---|---|",
  paste0("| Producing script | ", UNRESOLVED, " | No script in this repository produces `",
         basename(p_processed), "`. |"),
  paste0("| Log transformation | ", UNRESOLVED, " | The values are standardised per protein",
         " (see below), so they are not on a log2 scale. Whether a log transformation and",
         " per-sample median centring came before that step is not recorded for this file. |"),
  paste0("| Normalisation | ", UNRESOLVED, " | The filename says `unnormalized`; the",
         " 2025-11-07 hemisphere-level ProTigy session applied `norm.data none`, and the",
         " per-sample-class limma statistics fit the matrix as given. What, if anything, was",
         " applied upstream is not recorded. |"),
  paste0("| Imputation | ", UNRESOLVED, " | No imputation log, seed or record accompanies",
         " the file. |"),
  paste0("| PCA adjustment | ", UNRESOLVED, " | The filename says `pcaAdjusted`. No script,",
         " parameter or log describing which component(s) were removed exists in the",
         " repository or the project tree. |"),
  paste0("| Software and version | ", UNKNOWN, " | Not recorded. |"),
  paste0("| Seed / RNG state | ", UNRESOLVED, " | See below. |"),
  paste0("| Date / version | 2025-11-07 (file mtime and the duplicate's folder name) |",
         " Weak evidence: a modification time, not a provenance record. |"),
  "",
  "## What the numbers say about the processing, when the filenames disagree",
  "",
  "Two things were measured directly from the matrix rather than read off its name:",
  "",
  paste0("1. **It is standardised per protein.** Each of the 5,349 proteins is standardised",
         " separately across the 96 acquisitions: row means are approximately 0 and row",
         " standard deviations approximately 1. So the values are neither log2 abundance nor",
         " raw intensity, and the token `unnormalized` in the filename is not a description",
         " of them. Corroborated independently by the differential statistics: limma's",
         " `AveExpr` is the mean over one sample class's ", stat_run$n_columns_in_model,
         " animal-level units, and its average over the four sample classes stays close to",
         " zero (largest absolute value over all proteins ",
         format(signif(ave$max_abs_mean_over_classes, 3), scientific = FALSE),
         ", checked by 13_validate_release.R), which requires a per-protein centred input.",
         " The values carry only 2 decimals, so on these bytes no row is exactly mean 0 /",
         " sd 1 -- the residual deviations are the size that 2-decimal rounding produces --",
         " and the release therefore does not claim a specific standardisation formula. See",
         " the reporting note below."),
  paste0("2. **A sibling file matches it to within rounding.** `",
         basename(p_perseus), "` at the project root (Feb 2025, 5,349 x 106) has identical",
         " protein identifiers in identical order, the same 96 sample columns, and values",
         " differing by at most 0.005 -- exactly what a 2-decimal serialisation of the same",
         " matrix looks like. Its name encodes a 70% one-group valid-value filter,",
         " imputation, ANOVA and z-scoring, which is consistent with what the numbers show."),
  "",
  "That is suggestive, and it is the best lead available, but a filename is still not a",
  "processing record: it names no software, no parameters, no seed and no operator. The",
  "edge is recorded in the lineage as a CANDIDATE ancestor and remains UNRESOLVED.",
  "",
  paste0("**Consequence for reporting.** The effect size the canonical tables call `log2fc`",
         " (and the statistical GCT stores as `logFC`) is a difference of group means on this",
         " standardised scale -- a ", RELEASE_EFFECT_SIZE$public_term, ", not a log2 ratio."),
  "",
  paste0("The publication release publishes it as `", RELEASE_EFFECT_SIZE$public_field,
         "` and records `source_statistic_field = \"", RELEASE_EFFECT_SIZE$source_field,
         "\"` alongside it. No value was changed; only the label. It is not",
         " called a \"standardised mean difference\" unqualified (that reads as Cohen's d,",
         " which scales by a pooled within-group SD) and not called \"z-scored\" (see above:",
         " the operation is unresolved and the bytes are rounded)."),
  "",
  "Methods-safe wording:",
  "",
  paste0("> ", RELEASE_EFFECT_SIZE$methods_sentence),
  "",
  "## The `umap_adjusted` versus `pcaAdjusted` discrepancy",
  "",
  "The file supplied to the journal was named",
  "`processed_protein_group_matrix_filtered_umap_adjusted.xlsx`; the file in the project",
  "tree is named `pg.matrix_filtered_pcaAdjusted_unnormalized`. Both are 5,349 x 96. The",
  "dimensions match, but **the names disagree about which dimensionality-reduction method",
  "was used to adjust the data**, and nothing in the project resolves the disagreement.",
  "This is recorded as a discrepancy rather than harmonised. Until the producing",
  "workflow is recovered, neither name should be taken as a description of the operation",
  "actually applied.",
  "",
  "## `01_preprocessing/01_impute.r` -- audited, and it is not the producer",
  "",
  "The mission for this audit was to determine whether that script produced the canonical",
  "matrix, whether a seed was recorded, whether its input still exists, and whether it",
  "corresponds to the old journal matrix. Reading the script settles all four:",
  "",
  paste0("1. **Did it produce the canonical matrix? No.** It writes one `.xlsx` PER SAMPLE",
         " CLASS, named `<date>_pgmatrix_imputed_<sample_class>_<n>samples_missing70pct.xlsx`,",
         " into `02_data/gct/`. No such files exist there. It never writes a `.gct`, and it",
         " performs **no PCA adjustment at all** -- so it cannot be the source of a file",
         " named `pcaAdjusted`."),
  paste0("2. **Was a seed recorded? No.** ", ifelse(has_rnorm,
         "The script imputes with `rnorm(n, mean = mean_obs - 1.8 * sd_obs, sd = 0.3 * sd_obs)`",
         "The script does not call rnorm"),
         " and ", ifelse(has_set_seed, "DOES call set.seed", "**contains no `set.seed` call**"),
         ". A repository-wide search finds `set.seed` only in the PCA, EWCE, enrichment and",
         " test code, never in the imputation path. Its output is therefore not",
         " reproducible even if the input were recovered."),
  paste0("3. **Does its input still exist? No.** The declared input is `",
         basename(impute_input_path), "` under `01_input/raw_proteomics/`, which is absent.",
         " The script already fails closed on this with an actionable message, and",
         " `run_pipeline_check.ps1` records the stage as SKIP for the same reason."),
  paste0("4. **Does it correspond to the old journal matrix?** Only by name. `pg.matrix_raw.txt`",
         " (5,747 x 96) exists at the project root and matches the old",
         " `processed_protein_group_matrix_raw.txt`, but it is a `.txt`, not the `.tsv` the",
         " script expects, and no run linking the two was ever recorded."),
  "",
  "### What was NOT done",
  "",
  "A seed was **not** added and the matrix was **not** regenerated. Stochastic imputation",
  "without a recorded seed cannot be reproduced after the fact; re-running it would produce",
  "a numerically different matrix and therefore a new analysis generation, invalidating",
  "every validated downstream artefact and the locked GCT hashes. That is out of scope for",
  "a packaging task.",
  "",
  "## Consequence for the deposition",
  "",
  "- Everything from the 5,349 x 96 processed matrix onward is reproducible and hash-verified.",
  "- The step from the 5,747-row search output to the 5,349-row processed matrix is",
  paste0("  **", UNRESOLVED, "**."),
  "- The step from acquisition to search output is **absent from the project entirely**.",
  "- The processed matrix should therefore be deposited and cited **as a given input**, with",
  "  its SHA256, rather than described as reproducible from the raw data.",
  "",
  "## What would close the gap",
  "",
  "| Missing item | Would resolve |",
  "|---|---|",
  "| The script or notebook that produced the 5,349 x 96 matrix | filtering, log, normalisation, imputation, PCA adjustment, software, seed |",
  paste0("| The original search/quantification software run log or configuration",
         " (e.g. DIA-NN `report.log.txt`, if DIA-NN was used) | search software identity",
         " and version, enzyme, modifications, tolerances, acquisition method |"),
  "| The 96 `.d` acquisition directories | the raw layer, and the instrument model from acquisition metadata |",
  ""
)
release_write_lines(gap, release_path("provenance", "UPSTREAM_PREPROCESSING_GAP.md"))
release_log("  wrote UPSTREAM_PREPROCESSING_GAP.md (01_impute.r: set.seed present = ",
            has_set_seed, ", rnorm present = ", has_rnorm, ")")

# --------------------------------------------------------------------------------------
# software and database versions
# --------------------------------------------------------------------------------------

# Built by the shared helper so the editor workbook's Software_Versions sheet and this
# table are the same rows from the same derivation. See R/release_software_versions.R.
software_versions <- release_build_software_versions(DATA_ROOT, REPO_ROOT)
rownames(software_versions) <- NULL
n_unknown <- sum(software_versions$version == UNKNOWN)
release_log("  software_versions: ", nrow(software_versions), " components, ", n_unknown,
            " UNKNOWN")

# --------------------------------------------------------------------------------------
# analysis parameters
# --------------------------------------------------------------------------------------

ap <- function(stage_name, parameter, value, recorded_by, evidence_path,
               status = "KNOWN_VERIFIED", notes = NA_character_) {
  data.frame(stage = stage_name, parameter = parameter, value = value,
             recorded_by = recorded_by,
             evidence_path = release_repo_relative(evidence_path, REPO_ROOT), status = status,
             notes = notes, stringsAsFactors = FALSE, check.names = FALSE)
}

yaml_lines <- readLines(animal_param_yaml, warn = FALSE)
yaml_tbl <- do.call(rbind, lapply(yaml_lines, function(l) {
  m <- regmatches(l, regexec("^([A-Za-z_][A-Za-z0-9_]*):\\s*(.*)$", l))[[1]]
  if (length(m) != 3L) return(NULL)
  ap("pooled animal-level ProTigy run, 2026-08-24 (superseded)", m[[2]], trimws(m[[3]]),
     "neha_proteome_parameters.yaml", animal_param_yaml)
}))

protigy_2025_tbl <- if (file.exists(protigy_params_2025)) {
  pl <- readLines(protigy_params_2025, warn = FALSE)
  kv <- pl[grepl("\t", pl)]
  do.call(rbind, lapply(kv, function(l) {
    parts <- strsplit(l, "\t", fixed = TRUE)[[1]]
    ap("ProTigy 2025-11-07 (hemisphere-level, superseded)", parts[[1]],
       ifelse(length(parts) > 1L, parts[[2]], ""), "params.txt", protigy_params_2025)
  }))
} else NULL

enrich_tbl <- do.call(rbind, lapply(seq_len(nrow(run_parameters)), function(i) {
  ap("enrichment (GSEA / ORA), GO_BP run (GO-BP)", run_parameters$parameter[i],
     run_parameters$value[i], "run_parameters.csv",
     file.path(ENRICH_ROOT, "audits", "run_parameters.csv"))
}))
enrich_cc_tbl <- do.call(rbind, lapply(seq_len(nrow(run_parameters_cc)), function(i) {
  ap("enrichment (GSEA / ORA), GO_CC run", run_parameters_cc$parameter[i],
     run_parameters_cc$value[i], "run_parameters.csv",
     file.path(ENRICH_ROOT_CC, "audits", "run_parameters.csv"))
}))
enrich_kegg_tbl <- do.call(rbind, lapply(seq_len(nrow(run_parameters_kegg)), function(i) {
  ap("enrichment (GSEA), KEGG run (its GO-BP tables are superseded)",
     run_parameters_kegg$parameter[i], run_parameters_kegg$value[i], "run_parameters.csv",
     file.path(KEGG_ROOT, "audits", "run_parameters.csv"))
}))
# Not recorded by run_parameters.csv. The seed handling is taken from what the run produced:
# a rerun with the recorded seeds and DOSE seed = FALSE reproduces the stored NES. The script
# in the repository is later than the locked runs, so its hash is labelled as such. The
# features it adds are checked against the KEGG run, which predates them; the GO runs record
# them.
p_enrich_script <- file.path(REPO_ROOT, "04_differential_expression_enrichment",
                             "01_clusterProfiler.r")
fmt_e <- function(x) format(signif(as.numeric(x), 2), scientific = TRUE)
later_features <- c(
  if (!"gsea_n_perm_simple" %in% run_parameters_kegg$parameter)
    "sets fgsea nPermSimple",
  if (!"gsea_biocparallel_backend" %in% run_parameters_kegg$parameter)
    "registers a serial BiocParallel backend",
  if (!any(grepl("untested_set_count", names(release_read_csv(p_enrich_index_kegg)), fixed = TRUE)))
    "records the gene sets fgsea could not test")
if (!all(c("gsea_n_perm_simple", "gsea_biocparallel_backend") %in% run_parameters$parameter) ||
    !any(grepl("untested_set_count", names(release_read_csv(p_enrich_index)), fixed = TRUE))) {
  stop("The GO enrichment run does not record its nPermSimple, backend and untested sets.",
       call. = FALSE)
}
seed_tbl <- rbind(
  ap("enrichment (GSEA / ORA)", "dose_seed_argument", "FALSE",
     paste0("rerun check: the mCherry learning GO-BP and GO-CC GSEA, rerun with the recorded",
            " seeds and DOSE seed = FALSE, reproduce the stored FDR sets and NES (max |dNES| ",
            fmt_e(seed_checks$value[seed_check_rows[["GO-BP"]]]), " GO-BP, ",
            fmt_e(seed_checks$value[seed_check_rows[["GO-CC"]]]), " GO-CC)"),
     p_seed_checks, "KNOWN_VERIFIED",
     paste0(RELEASE_ENRICHMENT_SEED_HANDLING, " Evidence SHA256 ", hash_or(p_seed_checks), ".")),
  ap("enrichment (GSEA / ORA)", "enrichment_script_sha256_at_build", hash_or(p_enrich_script),
     "script at build time, later than the locked runs", p_enrich_script, "KNOWN_VERIFIED",
     if (length(later_features)) paste0(
       "Not the exact code of the locked runs. The GO runs record what this version adds;",
       " the KEGG run predates it: this version also ",
       if (length(later_features) > 1L) paste0(
         paste(later_features[-length(later_features)], collapse = ", "), " and ",
         later_features[length(later_features)]) else later_features,
       ", which the KEGG run's audits do not record. It also has no default output root",
       " (PROTEOMICS_ENRICHMENT_OUTPUT_ROOT must name the run folder). The runs are",
       " identified by their locked index hashes.") else
       "The locked runs' audits record every setting this version adds."))

# The per-sample-class limma model behind every released t, P and FDR. The model terms come
# from RELEASE_PRIMARY_MODEL and the script; the numbers from the run manifest.
LIMMA_STAGE <- "differential statistics (limma, one model per sample class)"
LIMMA_TESTED_BY <- function(call) {
  paste0("script source (", call, "); tests/test_limma_stat_gct.R requires the script to",
         " reproduce the locked statistics GCT byte for byte")
}
sm_stat <- stat_run$manifest
limma_tbl <- rbind(
  ap(LIMMA_STAGE, "variance_model", RELEASE_PRIMARY_MODEL$variance_model, "run manifest",
     p_stat_manifest),
  ap(LIMMA_STAGE, "design", "~0 + condition, fitted separately within each sample class",
     "script source", p_limma_script),
  ap(LIMMA_STAGE, "n_columns_in_model", as.character(stat_run$n_columns_in_model),
     "run manifest", p_stat_manifest),
  ap(LIMMA_STAGE, "df_residual", as.character(stat_run$df_residual), "run manifest",
     p_stat_manifest),
  # The run manifest does not record the eBayes settings; the script does, and
  # tests/test_limma_stat_gct.R checks that the script reproduces the locked GCT byte for
  # byte, which ties the script's settings to the published statistics.
  ap(LIMMA_STAGE, "ebayes_robust", "TRUE",
     LIMMA_TESTED_BY("eBayes call"), paste(p_limma_script, p_limma_test, sep = ";")),
  ap(LIMMA_STAGE, "ebayes_trend", "FALSE",
     LIMMA_TESTED_BY("eBayes call, limma default"), paste(p_limma_script, p_limma_test, sep = ";")),
  ap(LIMMA_STAGE, "p_adjust_method", "BH, within each contrast",
     LIMMA_TESTED_BY("topTable call"), paste(p_limma_script, p_limma_test, sep = ";")),
  ap(LIMMA_STAGE, "significance_threshold", "adj.P.Val < 0.05", "script source",
     p_limma_script),
  ap(LIMMA_STAGE, "limma_version", stat_run$limma_version, "run manifest", p_stat_manifest),
  ap(LIMMA_STAGE, "r_version", stat_run$r_version_string, "run manifest", p_stat_manifest),
  do.call(rbind, lapply(unique(sub("_.*$", "", sm_stat$canonical_comparison)), function(cls) {
    r <- sm_stat[startsWith(sm_stat$canonical_comparison, paste0(cls, "_")), , drop = FALSE][1, ]
    rbind(ap(LIMMA_STAGE, paste0("df_prior_median_", cls),
             format(signif(as.numeric(r$df_prior_median), 4)), "run manifest",
             p_stat_manifest),
          ap(LIMMA_STAGE, paste0("s2_prior_", cls),
             format(signif(as.numeric(r$s2_prior), 4)), "run manifest", p_stat_manifest))
  }))
)

ewce_lines <- if (file.exists(ewce_session_path)) readLines(ewce_session_path, warn = FALSE) else character(0)
grab <- function(prefix) {
  hit <- grep(paste0("^", prefix), ewce_lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(paste0("^", prefix, "\\s*:?\\s*"), "", hit[[1]]))
}
ewce_tbl <- do.call(rbind, lapply(
  list(c("bootstrap_replicates", "EWCE reps"), c("top_n_values", "Top-N values"),
       c("annotation_levels", "Annotation levels"), c("primary_top_n", "Primary top-N"),
       c("primary_annotation_level", "Primary annotation level"),
       c("measured_proteome_background_size", "Measured proteome background size"),
       c("aggregation_policy", "Aggregation policy"),
       c("post_aggregation_transformations", "Post-aggregation transformations"),
       c("sampling_unit", "Sampling unit"), c("input_sha256", "Input SHA256")),
  function(x) ap("EWCE", x[[1]], grab(x[[2]]), "reproducibility_session_info.txt",
                 ewce_session_path)))
# EWCE's differential lists come from its own limma fit, not from the released statistics.
ewce_tbl <- rbind(ewce_tbl,
  ap("EWCE", "differential_limma_design", ewce_model$design,
     "animal_level_differential_audit.csv", ewce_model$path, "KNOWN_VERIFIED",
     ewce_model$sentence),
  ap("EWCE", "differential_limma_ebayes", ewce_model$ebayes,
     "animal_level_differential_audit.csv", ewce_model$path, "KNOWN_VERIFIED",
     paste("EWCE's own fit; the released differential statistics use",
           "eBayes(robust = TRUE, trend = FALSE).")),
  ap("EWCE", "differential_n_proteins", as.character(ewce_model$n_proteins),
     "animal_level_differential_audit.csv", ewce_model$path))

impute_tbl <- rbind(
  ap("historical imputation (01_impute.r, NOT the producer of the canonical matrix)",
     "missingness_filter", "<= 70% missing within sample class", "script source",
     impute_script, "KNOWN_BUT_NEEDS_STANDARDIZATION",
     "Read from the script, not from a run record. See UPSTREAM_PREPROCESSING_GAP.md."),
  ap("historical imputation (01_impute.r, NOT the producer of the canonical matrix)",
     "transformation", "log2 then per-sample median centring", "script source",
     impute_script, "KNOWN_BUT_NEEDS_STANDARDIZATION", NA_character_),
  ap("historical imputation (01_impute.r, NOT the producer of the canonical matrix)",
     "imputation", "downshifted normal: mean - 1.8 * sd, sd = 0.3 * sd, per column",
     "script source", impute_script, "KNOWN_BUT_NEEDS_STANDARDIZATION", NA_character_),
  ap("historical imputation (01_impute.r, NOT the producer of the canonical matrix)",
     "seed", ifelse(has_set_seed, "recorded in script", "NONE -- not reproducible"),
     "script source", impute_script,
     ifelse(has_set_seed, "KNOWN_VERIFIED", "MISSING_UNKNOWN"),
     "No set.seed anywhere in the imputation path."),
  ap("upstream processed matrix", "pca_adjustment", UNRESOLVED, "no record exists", "NONE",
     "MISSING_UNKNOWN", "Filename says pcaAdjusted; no script or parameter documents it."),
  ap("upstream processed matrix", "normalisation", UNRESOLVED, "no record exists", "NONE",
     "MISSING_UNKNOWN", "Filename says unnormalized; not independently recorded.")
)

analysis_parameters <- rbind(limma_tbl, yaml_tbl, protigy_2025_tbl, enrich_tbl, enrich_cc_tbl,
                             enrich_kegg_tbl, seed_tbl, ewce_tbl, impute_tbl)
rownames(analysis_parameters) <- NULL
release_log("  analysis_parameters: ", nrow(analysis_parameters), " recorded parameters")

# --------------------------------------------------------------------------------------
# sessionInfo for THIS build
# --------------------------------------------------------------------------------------

si <- c(
  "Publication release build environment",
  paste0("Generated: ", release_timestamp_utc()),
  paste0("Repository: proteomic_analysis/ (repository paths in this release are relative ",
         "to this folder)"),
  paste0("Git commit: ", release_git_commit(REPO_ROOT)),
  paste0("Release layer version: ", RELEASE_LAYER_VERSION),
  "",
  "NOTE: this is the environment that BUILT this package by reading canonical outputs.",
  "It is NOT the environment that computed the science. The versions that produced each",
  "canonical result are in provenance/software_versions.tsv, recorded by those runs.",
  "",
  capture.output(utils::sessionInfo())
)
release_write_lines(si, release_path("provenance", "sessionInfo_release.txt"))

# --------------------------------------------------------------------------------------
# write
# --------------------------------------------------------------------------------------

w1 <- release_write_table(lineage, release_path("provenance", "data_lineage.tsv"))
release_register("provenance/data_lineage.tsv",
                 "acquisition to publication lineage; unproven edges marked UNRESOLVED",
                 c(experimenter_source, protocol_source, p_pg_raw, p_processed,
                   p_animal_gct, p_stat_gct, p_stat_manifest, release_figure_rebuild_record(DATA_ROOT),
                   p_stat_gct_pooled, p_pooled_mapped_index, p_pooled_enrich_index,
                   ewce_model$path),
                 c(hash_or(experimenter_source), hash_or(protocol_source), hash_or(p_pg_raw),
                   processed_sha, hash_or(p_animal_gct), hash_or(p_stat_gct),
                   hash_or(p_stat_manifest), hash_or(release_figure_rebuild_record(DATA_ROOT)),
                   hash_or(p_stat_gct_pooled), hash_or(p_pooled_mapped_index),
                   hash_or(p_pooled_enrich_index), hash_or(ewce_model$path)),
                 STAGE, "tsv")

w2 <- release_write_table(software_versions,
                          release_path("provenance", "software_versions.tsv"))
release_register("provenance/software_versions.tsv",
                 "exact software and reference-database versions, or UNKNOWN with a source",
                 c(file.path(ENRICH_ROOT, "audits", "package_database_versions.csv"),
                   file.path(ENRICH_ROOT_CC, "audits", "package_database_versions.csv"),
                   p_stat_manifest, ewce_session_path, pca_session_path, protigy_params_2025),
                 c(hash_or(file.path(ENRICH_ROOT, "audits", "package_database_versions.csv")),
                   hash_or(file.path(ENRICH_ROOT_CC, "audits", "package_database_versions.csv")),
                   hash_or(p_stat_manifest), hash_or(ewce_session_path),
                   hash_or(pca_session_path), hash_or(protigy_params_2025)), STAGE, "tsv")

w3 <- release_write_table(analysis_parameters,
                          release_path("provenance", "analysis_parameters.tsv"))
release_register("provenance/analysis_parameters.tsv",
                 "every recorded analysis parameter, per stage",
                 c(p_stat_manifest, p_limma_script, p_limma_test, p_seed_checks,
                   p_enrich_script, animal_param_yaml,
                   protigy_params_2025, file.path(ENRICH_ROOT, "audits", "run_parameters.csv"),
                   file.path(ENRICH_ROOT_CC, "audits", "run_parameters.csv"),
                   file.path(KEGG_ROOT, "audits", "run_parameters.csv"),
                   ewce_session_path, ewce_model$path),
                 c(hash_or(p_stat_manifest), hash_or(p_limma_script), hash_or(p_limma_test),
                   hash_or(p_seed_checks), hash_or(p_enrich_script),
                   hash_or(animal_param_yaml), hash_or(protigy_params_2025),
                   hash_or(file.path(ENRICH_ROOT, "audits", "run_parameters.csv")),
                   hash_or(file.path(ENRICH_ROOT_CC, "audits", "run_parameters.csv")),
                   hash_or(file.path(KEGG_ROOT, "audits", "run_parameters.csv")),
                   hash_or(ewce_session_path), hash_or(ewce_model$path)), STAGE, "tsv")

release_register("provenance/UPSTREAM_PREPROCESSING_GAP.md",
                 "what upstream of the 5,349-protein matrix cannot be reproduced, and why",
                 c(impute_script, p_processed), c(hash_or(impute_script), processed_sha),
                 STAGE, "md")
release_register("provenance/sessionInfo_release.txt",
                 "the R environment that built this package (not the science)",
                 "NONE", NA_character_, STAGE, "txt")

release_log("  wrote data_lineage.tsv (", w1$rows, "x", w1$cols, ")")
release_log("  wrote software_versions.tsv (", w2$rows, "x", w2$cols, ")")
release_log("  wrote analysis_parameters.tsv (", w3$rows, "x", w3$cols, ")")
release_log("  wrote sessionInfo_release.txt, UPSTREAM_PREPROCESSING_GAP.md")
release_log("stage 10 complete")
