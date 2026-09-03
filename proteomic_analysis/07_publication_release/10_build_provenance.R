#!/usr/bin/env Rscript

# Publication release, stage 10 -- lineage, software/database versions, parameters.
#
# Produces
#   provenance/data_lineage.tsv                 raw -> processed -> analysis, edge by edge
#   provenance/UPSTREAM_PREPROCESSING_GAP.md    what cannot be reproduced, and why
#   provenance/software_versions.tsv            exact versions, or UNKNOWN with a source
#   provenance/sessionInfo_release.txt          the environment that built THIS package
#   provenance/analysis_parameters.tsv          every recorded analysis parameter
#
# Lineage edges that cannot be proven are marked UNRESOLVED rather than inferred. Two of
# them are load-bearing and are stated as such: the search-software output is absent from
# the project, and the code that produced the 5,349-protein processed matrix is not in the
# repository.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_validation.R"))
release_require("digest")

DATA_ROOT <- release_data_root()
PROJECT_ROOT <- release_project_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/10_build_provenance.R"
UNRESOLVED <- "UNRESOLVED"
UNKNOWN <- "UNKNOWN"

release_banner("stage 10 -- provenance")

ANIMAL_ROOT <- file.path(DATA_ROOT, "02_data", "animal_level")
ENRICH_ROOT <- file.path(DATA_ROOT, "03_output", "enrichment",
                         "enrichment_t_rank_validation_20260825")
EWCE_ROOT <- file.path(DATA_ROOT, "03_output", "ewce",
                       "EWCE_Results_animal_level_validation_20260825")
PCA_ROOT <- file.path(DATA_ROOT, "03_output", "pca",
                      "pca_plots_animal_level_validation_20260825_rerun")
REVISION_ROOT <- file.path(DATA_ROOT, "03_output", "reviewer_revision_animal_level_20260827")
experimenter_source <- release_experimenter_metadata_source_path(REPO_ROOT)
protocol_source <- release_sample_preparation_source_path(REPO_ROOT)
experimenter_metadata <- release_read_experimenter_metadata(REPO_ROOT)
sample_preparation_protocol <- release_read_sample_preparation_protocol(REPO_ROOT)

hash_or <- function(path, fallback = NA_character_) {
  if (!is.na(path) && nzchar(path) && file.exists(path) && !dir.exists(path)) {
    release_sha256(path)
  } else fallback
}

mapped_index <- release_read_csv(file.path(ANIMAL_ROOT, "mapped", "indexMappedComparisons.csv"))
pkg_versions <- release_read_csv(file.path(ENRICH_ROOT, "audits",
                                           "package_database_versions.csv"))
run_parameters <- release_read_csv(file.path(ENRICH_ROOT, "audits", "run_parameters.csv"))
protigy_params_2025 <- file.path(DATA_ROOT, "01_input", "raw_proteomics",
                                 "20251107_pg.matrix_Neha", "params.txt")
animal_param_yaml <- file.path(ANIMAL_ROOT, "neha_proteome_parameters.yaml")
ewce_session_path <- file.path(EWCE_ROOT, "03_QC_Mapping_Logs",
                               "reproducibility_session_info.txt")
pca_session_path <- file.path(PCA_ROOT, "tables", "meta", "sessionInfo.txt")

# --------------------------------------------------------------------------------------
# ProTigy version for the CANONICAL animal-level run (2026-08-24).
#
# Recovered 2026-09-02 by targeted audit. It is recorded as a constant rather than read
# from whatever ProTigy happens to be installed on the machine running this build: the
# question is which version produced a specific run in the past, and a build-time lookup
# would silently answer a different question on a different machine. That is the same trap
# as carrying the 2025 version forward.
#
# The evidence is five strands that agree:
#   1. The only ProTigy present on the analysis machine is 2.4.1, and its installed
#      DESCRIPTION records Packaged 2026-08-24 12:47:48 UTC / Built 2026-08-24 12:48:07 UTC.
#   2. The canonical run's parameter export, neha_proteome_parameters.yaml, was written
#      2026-08-24 13:06:11 UTC -- 18 minutes after that build finished.
#   3. ProTigy v2's exporter (R/tab_export.R) writes paste0(ome, "_parameters.yaml");
#      the canonical file is neha_proteome_parameters.yaml, i.e. the historical project-specific proteome identifier.
#   4. That exporter writes the parameter list minus "gct_file_path". The canonical YAML
#      carries gct_file_name and no gct_file_path, exactly as the code does.
#   5. All 19 schema keys in the canonical YAML are a subset of 2.4.1's
#      setup_parameters/setupDefaults.yaml; the two extras (gct_file_name,
#      annotation_column) are added at runtime by 2.4.1; and the two absent data-filter
#      keys are precisely the ones 2.4.1 sets to NULL when data_filter is None, which is
#      this run's setting.
#
# The negative result matters as much: v1.1.x writes a params.txt whose first lines are
# "## <timestamp>" / "## Protigy (vX.Y.Z)". The canonical run produced no such file, and
# its YAML key set does not exist in the v1.1.x format at all. v1.1.8 is therefore
# disproven for this run, not merely unproven.
#
# Recorded evidence_path stays the canonical YAML because that is the durable, hashed
# artefact on shared storage; the installed-package DESCRIPTION is machine-local and is
# named in the notes instead.
PROTIGY_ANIMAL_LEVEL_VERSION <- "2.4.1"
PROTIGY_ANIMAL_LEVEL_EVIDENCE <- paste(
  "Recovered 2026-09-02 by cross-source audit, not read from the build machine.",
  "The installed Protigy DESCRIPTION (R library, sha256",
  "74ac5f7c35dfeb06e62e472c5a072b136e52ccc8965531e0a5380dc4b65da37d) reports version",
  "2.4.1, Packaged 2026-08-24 12:47:48 UTC and Built 2026-08-24 12:48:07 UTC; this run's",
  "parameter export was written 18 minutes later, at 2026-08-24 13:06:11 UTC. The export",
  "is a v2-only artefact: ProTigy v2's tab_export.R writes <ome>_parameters.yaml with the",
  "parameter list minus gct_file_path, which is exactly this file's name and key set, and",
  "all 19 of its schema keys are a subset of 2.4.1's setupDefaults.yaml. v1.1.8 is",
  "DISPROVEN for this run rather than merely unproven: v1.1.x emits a params.txt with a",
  "'## Protigy (vX.Y.Z)' header and none of these keys. See the recovery audit for the",
  "full evidence table."
)

# Fail closed if the artefact this claim rests on stops looking like a ProTigy v2 export.
# Without this, a future change to the canonical YAML could leave a version claim standing
# on evidence that no longer exists.
if (file.exists(animal_param_yaml)) {
  .pp <- readLines(animal_param_yaml, warn = FALSE)
  .has_v2_keys <- all(vapply(c("gct_file_name:", "annotation_column:", "group_normalization:",
                               "convert_ids_to_gene_symbol:", "id_mapping_species:"),
                             function(k) any(startsWith(.pp, k)), logical(1)))
  # The v1.1.x marker is a literal "## Protigy (vX.Y.Z)" comment line; matched with
  # startsWith on trimmed lines so this check carries no regex escaping of its own.
  .has_v1_header <- any(startsWith(trimws(.pp), "## Protigy (v"))
  .has_gct_file_path <- any(startsWith(.pp, "gct_file_path:"))
  if (!.has_v2_keys || .has_v1_header || .has_gct_file_path) {
    stop("The animal-level ProTigy parameter export no longer carries the ProTigy v2 ",
         "export signature that the recovered version ", PROTIGY_ANIMAL_LEVEL_VERSION,
         " rests on (", animal_param_yaml, "). Re-verify the version before releasing.",
         call. = FALSE)
  }
  rm(.pp, .has_v2_keys, .has_v1_header, .has_gct_file_path)
}

# --------------------------------------------------------------------------------------
# lineage
# --------------------------------------------------------------------------------------

lin <- function(artifact_id, artifact_type, path_or_release_filename, parent_artifact,
                processing_stage, software, software_version, analysis_unit,
                sha256 = NA_character_, canonical = "yes", notes = NA_character_) {
  data.frame(artifact_id = artifact_id, artifact_type = artifact_type,
             path_or_release_filename = path_or_release_filename,
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
p_stat_gct <- release_locked_artefact_path("protigy_stat_gct", DATA_ROOT)
p_split_index <- file.path(ANIMAL_ROOT, "split", "indexComparisons.csv")
p_mapped_index <- file.path(ANIMAL_ROOT, "mapped", "indexMappedComparisons.csv")
p_idmapping <- file.path(DATA_ROOT, "01_input", "references", "MOUSE_10090_idmapping.dat")
p_manual_mapping <- file.path(DATA_ROOT, "01_input", "references", "manual_mapping.xlsx")
p_enrich_index <- file.path(ENRICH_ROOT, "indexEnrichmentComparisons.csv")
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
      "experimenter-supplied facts codified in a source-controlled release contract",
      "manual curation", "n/a", "study / cohort / sample class",
      hash_or(experimenter_source), "yes",
      paste("Defines CeM origin, marker-defined LCM classes, mouse strain/source, cohort-only",
            "sex and age-at-start facts, label-free status, explicit unknowns and the",
            "unrelated rat-experiment exclusion. No AnimalID-level sex was invented.")),
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
      paste("filtering, log2 transformation, median centring, imputation, PCA adjustment --",
            "PRODUCING CODE NOT IN REPOSITORY:", UNRESOLVED),
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
  lin("A09_protigy_statistics", "differential statistics GCT",
      p_stat_gct, "A08_animal_level_matrix",
      "within-sample-class two-sample moderated t", "ProTigy",
      PROTIGY_ANIMAL_LEVEL_VERSION, "animal",
      hash_or(p_stat_gct), "yes",
      paste("LOCKED artefact. Parameters recorded in",
            "02_data/animal_level/neha_proteome_parameters.yaml. ProTigy VERSION for THIS",
            "(2026-08-24) run is", PROTIGY_ANIMAL_LEVEL_VERSION, "-- recovered 2026-09-02",
            "by cross-source audit and NOT carried over from the 2025 runs, which used",
            "v1.1.5 (2025-03-28) and v1.1.8 (2025-11-07, 2025-12-12). v1.1.8 is disproven",
            "for this run: it emits a params.txt header, this run emitted a v2 YAML export.",
            "See provenance/software_versions.tsv for the evidence.")),
  lin("A10_split_differential_tables", "differential result tables",
      file.path(ANIMAL_ROOT, "split", "forward"), "A09_protigy_statistics",
      "split into 12 forward and 12 reverse per-contrast tables",
      "R + 01_preprocessing/03_gct_extractR.r", paste("R", R.version.string), "animal",
      hash_or(p_split_index), "yes",
      "5,349 rows per contrast. Index: indexComparisons.csv."),
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
      file.path(ANIMAL_ROOT, "mapped", "forward"),
      "A10_split_differential_tables;A11_uniprot_mapping_reference;A11b_manual_mapping",
      "UniProt identifier resolution", "R + 02_id_mapping/01_MapThatProt_batch.r",
      paste("R", R.version.string), "animal", hash_or(p_mapped_index), "yes",
      paste("5,327 of 5,349 mapped; 22 non-mouse contaminant identifiers retained in",
            "mapped/unmapped/. Index: indexMappedComparisons.csv.")),
  lin("A13_enrichment", "GSEA and ORA results",
      ENRICH_ROOT, "A12_mapped_differential_tables",
      "GSEA ranked by the moderated t; ORA over FDR-significant lists",
      "R + clusterProfiler + fgsea",
      paste0("clusterProfiler ",
             pkg_versions$version[pkg_versions$component == "clusterProfiler"],
             "; fgsea ", pkg_versions$version[pkg_versions$component == "fgsea"]),
      "animal", hash_or(p_enrich_index), "yes",
      paste("The effect-size-ranked GSEA in the same run is a SENSITIVITY analysis only;",
            "its ranking statistic is the standardised-abundance effect size recorded",
            "upstream as log2fc, not a log2 fold change.",
            "Index: indexEnrichmentComparisons.csv.")),
  lin("A14_ewce", "cell-type enrichment results",
      EWCE_ROOT, "A08_animal_level_matrix;A10_split_differential_tables",
      "EWCE bootstrap enrichment, 10000 replicates",
      "R + EWCE + ewceData",
      paste0("EWCE 1.18.1; ewceData 1.18.0"), "animal", hash_or(p_ewce_table), "yes",
      paste("Cell-type reference is ewceData::ctd(), NOT the amygdala loom in 01_input.",
            "Primary settings: annotation level 2, top-N 250, FDR 0.05, measured-proteome",
            "background 4859.")),
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
  lin("A16_figure_panels", "manuscript figure panels and source data",
      REVISION_ROOT, "A08_animal_level_matrix;A12_mapped_differential_tables;A13_enrichment;A14_ewce",
      "figure regeneration at animal level",
      "R + 06_manuscript_figure_revision/*.R (frozen snapshots)",
      paste("R", R.version.string), "animal",
      hash_or(file.path(REVISION_ROOT, "manifests", "regenerated_panel_manifest.csv")), "yes",
      "17 panels; per-panel source data copied into editor_source_data/figure_source_data/."),
  lin("A16b_acquisition_qc_counts", "acquisition QC source data",
      p_protein_count, "A03_search_quantification",
      "protein and peptide identification counts per acquisition",
      UNKNOWN, UNKNOWN, "acquisition run", NA_character_, "yes",
      paste("Source of Supplementary A. Lives outside the clusterProfiler tree.",
            "Producing step is", UNRESOLVED, "-- not hashed here; the copied panel source",
            "data in editor_source_data/ carries its hash.")),
  lin("A17_publication_release", "publication / deposition package",
      "this release", paste("A00_experimenter_metadata;A00_sample_preparation_protocol;",
                            "A05_sample_annotation;A06_sample_info;A04_protein_group_matrix_prefilter;",
                            "A07_processed_matrix;A08_animal_level_matrix;A10_split_differential_tables;",
                            "A12_mapped_differential_tables;A13_enrichment;A14_ewce;A16_figure_panels",
                            sep = ""),
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
         " `max_missing: 70` field of `neha_proteome_parameters.yaml`. **The exact",
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
  paste0("| Log transformation | ", UNRESOLVED, " | The values are log2-scale and centred",
         " near zero, consistent with log2 plus per-sample median centring, but no record",
         " of the operation applied to this file exists. |"),
  paste0("| Normalisation | ", UNRESOLVED, " | The filename says `unnormalized`; the",
         " downstream ProTigy run applied `norm.data none`. What, if anything, was applied",
         " upstream is not recorded. |"),
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
         " of them. Corroborated independently by the canonical differential tables: limma's",
         " `AveExpr` lies within +/-0.0011 of zero for every protein, which requires a",
         " per-protein centred input. The values carry only 2 decimals, so on these bytes no",
         " row is exactly mean 0 / sd 1 -- the residual deviations are the size that",
         " 2-decimal rounding produces -- and the release therefore does not claim a",
         " specific standardisation formula. See the reporting note below."),
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
         " (and ProTigy stores as `logFC`) is a difference of group means on this",
         " standardised scale -- a ", RELEASE_EFFECT_SIZE$public_term, ", not a log2 ratio."),
  "",
  paste0("The publication release publishes it as `", RELEASE_EFFECT_SIZE$public_field,
         "` and records `source_statistic_field = \"", RELEASE_EFFECT_SIZE$source_field,
         "\"` alongside it. No value was changed; only the label. It is deliberately not",
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
  "This is recorded as a discrepancy rather than silently harmonised. Until the producing",
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
  "### What was deliberately NOT done",
  "",
  "A seed was **not** added and the matrix was **not** regenerated. Stochastic imputation",
  "without a recorded seed cannot be reproduced after the fact; re-running it would produce",
  "a numerically different matrix and therefore a new analysis generation, invalidating",
  "every validated downstream artefact and both locked GCT hashes. That is out of scope for",
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

sv <- function(component, category, version, status, recorded_by, evidence_path,
               applies_to, notes = NA_character_) {
  data.frame(component = component, category = category, version = version, status = status,
             recorded_by = recorded_by, evidence_path = evidence_path,
             evidence_sha256 = hash_or(evidence_path), applies_to_stage = applies_to,
             notes = notes, stringsAsFactors = FALSE, check.names = FALSE)
}

parse_attached <- function(path) {
  if (!file.exists(path)) return(character(0))
  lines <- readLines(path, warn = FALSE)
  start <- grep("^other attached packages", lines)
  if (!length(start)) return(character(0))
  tail_lines <- lines[(start[[1]] + 1L):length(lines)]
  stop_at <- grep("^loaded via a namespace", tail_lines)
  if (length(stop_at)) tail_lines <- tail_lines[seq_len(stop_at[[1]] - 1L)]
  toks <- unlist(strsplit(paste(tail_lines, collapse = " "), "\\s+"))
  unique(toks[grepl("^[A-Za-z][A-Za-z0-9._]*_[0-9]", toks)])
}

enrich_versions <- do.call(rbind, lapply(seq_len(nrow(pkg_versions)), function(i) {
  sv(pkg_versions$component[i],
     ifelse(pkg_versions$component[i] == "R", "language", "R package"),
     pkg_versions$version[i], "KNOWN_VERIFIED",
     "enrichment run audit, 2026-08-25",
     file.path(ENRICH_ROOT, "audits", "package_database_versions.csv"),
     "differential enrichment (GSEA / ORA)")
}))

ewce_toks <- parse_attached(ewce_session_path)
ewce_versions <- if (length(ewce_toks)) {
  do.call(rbind, lapply(ewce_toks, function(tok) {
    sv(sub("_.*$", "", tok), "R package", sub("^[^_]*_", "", tok), "KNOWN_VERIFIED",
       "EWCE run sessionInfo, 2026-08-25", ewce_session_path, "EWCE cell-type enrichment")
  }))
} else NULL

pca_toks <- parse_attached(pca_session_path)
pca_versions <- if (length(pca_toks)) {
  do.call(rbind, lapply(pca_toks, function(tok) {
    sv(sub("_.*$", "", tok), "R package", sub("^[^_]*_", "", tok), "KNOWN_VERIFIED",
       "PCA run sessionInfo, 2026-08-25", pca_session_path, "PCA")
  }))
} else NULL

protigy_version_2025 <- NA_character_
if (file.exists(protigy_params_2025)) {
  hit <- grep("Protigy", readLines(protigy_params_2025, warn = FALSE), value = TRUE,
              ignore.case = TRUE)
  if (length(hit)) {
    m <- regmatches(hit[[1]], regexpr("v[0-9][0-9.]*", hit[[1]]))
    if (length(m)) protigy_version_2025 <- m
  }
}

external <- rbind(
  sv("ProTigy (hemisphere-level runs, 2025-11-07 and 2025-12-12)", "external application",
     ifelse(is.na(protigy_version_2025), UNKNOWN, protigy_version_2025), "KNOWN_VERIFIED",
     "ProTigy params.txt header", protigy_params_2025,
     "hemisphere-level ProTigy statistics (superseded)",
     paste("Both hemisphere-level params.txt headers report this same version. The",
           "2026-09-02 audit also found seven older params.txt files under the project's",
           "protigy/ folder, all reporting v1.1.5 for exploratory runs on 2025-03-28.",
           "The version therefore changed across this project's history",
           "(1.1.5 -> 1.1.8 -> 2.4.1), which is why no version is carried between runs.")),
  sv("ProTigy (animal-level statistical GCT, 2026-08-24)", "external application",
     PROTIGY_ANIMAL_LEVEL_VERSION, "KNOWN_VERIFIED",
     "recovered by cross-source audit, 2026-09-02",
     animal_param_yaml, "canonical animal-level differential statistics",
     PROTIGY_ANIMAL_LEVEL_EVIDENCE),
  sv("upstream search / quantification software", "external application", UNKNOWN,
     "MISSING_RECOVERABLE",
     "no software name and no version string exist in the project tree", "NONE",
     "peptide/protein identification and quantification",
     paste("The retained processed files use the historical `pg.matrix` naming convention,",
           "but the exact upstream search/quantification software and configuration could",
           "not be recovered from the retained project files. Recoverable from the",
           "acquisition facility: the original search/quantification software run log or",
           "configuration (e.g. DIA-NN report.log.txt, if DIA-NN was used).")),
  sv("MS instrument", "instrument", UNKNOWN, "MISSING_RECOVERABLE",
     "no instrument model recorded anywhere", "NONE", "LC-MS acquisition",
     paste("Run-name alias `Olive` is not a model. Acquisition format is `.d`, a",
           "vendor-specific acquisition directory format.")),
  sv("UniProt idmapping (MOUSE_10090)", "reference database",
     as.character(mapped_index$mapping_reference_version[[1]]), "KNOWN_BUT_NEEDS_STANDARDIZATION",
     paste0("mapping run record; snapshot ",
            mapped_index$mapping_reference_snapshot_date_utc[[1]]),
     p_idmapping, "UniProt identifier mapping",
     paste0("SHA256 ", mapped_index$mapping_reference_sha256[[1]], "; ",
            mapped_index$mapping_reference_bytes[[1]], " bytes; modified ",
            mapped_index$mapping_reference_modified_utc[[1]],
            ". The file does not encode a UniProt release number.")),
  sv("manual identifier overrides", "reference override", "n/a", "KNOWN_VERIFIED",
     "mapping run record", p_manual_mapping, "UniProt identifier mapping",
     paste0(mapped_index$manual_mapping_rows[[1]], " rows; SHA256 ",
            mapped_index$manual_mapping_sha256[[1]])),
  sv("ewceData::ctd()", "reference dataset", "ewceData 1.18.0", "KNOWN_VERIFIED",
     "EWCE run sessionInfo and 05_celltype_enrichment_EWCE/01_EWCE.r line 115",
     ewce_session_path, "EWCE cell-type enrichment",
     paste("The canonical EWCE run uses the packaged ewceData CTD. The",
           "l1_amygdala.loom file in 01_input/single_cell/ is NOT used by it.")),
  sv("R (publication release build)", "language",
     paste(R.version$major, R.version$minor, sep = "."), "KNOWN_VERIFIED",
     "this build", "NONE", "publication release layer",
     "The environment that produced THIS package; see sessionInfo_release.txt.")
)

software_versions <- rbind(enrich_versions, ewce_versions, pca_versions, external)
software_versions <- software_versions[
  !duplicated(paste(software_versions$component, software_versions$version,
                    software_versions$applies_to_stage)), , drop = FALSE]
software_versions <- software_versions[order(software_versions$applies_to_stage,
                                              software_versions$component), , drop = FALSE]
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
             recorded_by = recorded_by, evidence_path = evidence_path, status = status,
             notes = notes, stringsAsFactors = FALSE, check.names = FALSE)
}

yaml_lines <- readLines(animal_param_yaml, warn = FALSE)
yaml_tbl <- do.call(rbind, lapply(yaml_lines, function(l) {
  m <- regmatches(l, regexec("^([A-Za-z_][A-Za-z0-9_]*):\\s*(.*)$", l))[[1]]
  if (length(m) != 3L) return(NULL)
  ap("animal-level ProTigy input", m[[2]], trimws(m[[3]]),
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
  ap("enrichment (GSEA / ORA)", run_parameters$parameter[i], run_parameters$value[i],
     "run_parameters.csv", file.path(ENRICH_ROOT, "audits", "run_parameters.csv"))
}))

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

analysis_parameters <- rbind(yaml_tbl, protigy_2025_tbl, enrich_tbl, ewce_tbl, impute_tbl)
rownames(analysis_parameters) <- NULL
release_log("  analysis_parameters: ", nrow(analysis_parameters), " recorded parameters")

# --------------------------------------------------------------------------------------
# sessionInfo for THIS build
# --------------------------------------------------------------------------------------

si <- c(
  "Publication release build environment",
  paste0("Generated: ", release_timestamp_utc()),
  paste0("Repository: ", REPO_ROOT),
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
                   p_animal_gct, p_stat_gct),
                 c(hash_or(experimenter_source), hash_or(protocol_source), hash_or(p_pg_raw),
                   processed_sha, hash_or(p_animal_gct), hash_or(p_stat_gct)), STAGE, "tsv")

w2 <- release_write_table(software_versions,
                          release_path("provenance", "software_versions.tsv"))
release_register("provenance/software_versions.tsv",
                 "exact software and reference-database versions, or UNKNOWN with a source",
                 c(file.path(ENRICH_ROOT, "audits", "package_database_versions.csv"),
                   ewce_session_path, pca_session_path, protigy_params_2025),
                 c(hash_or(file.path(ENRICH_ROOT, "audits", "package_database_versions.csv")),
                   hash_or(ewce_session_path), hash_or(pca_session_path),
                   hash_or(protigy_params_2025)), STAGE, "tsv")

w3 <- release_write_table(analysis_parameters,
                          release_path("provenance", "analysis_parameters.tsv"))
release_register("provenance/analysis_parameters.tsv",
                 "every recorded analysis parameter, per stage",
                 c(animal_param_yaml, protigy_params_2025,
                   file.path(ENRICH_ROOT, "audits", "run_parameters.csv"), ewce_session_path),
                 c(hash_or(animal_param_yaml), hash_or(protigy_params_2025),
                   hash_or(file.path(ENRICH_ROOT, "audits", "run_parameters.csv")),
                   hash_or(ewce_session_path)), STAGE, "tsv")

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
