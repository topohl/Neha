# Shared constants and checks for the release: status vocabularies, the experimenter
# records, locked input hashes, the sample-class correction, the original-submission
# crosswalk, collection-plate and effect-size wording, the primary model, enrichment notes,
# the no-recomputation scan and workbook checks.
#
# Status vocabularies:
#   release_metadata_statuses    what is known about one metadata or SDRF field
#   release_sdrf_field_statuses  READY / MISSING_REQUIRED_METADATA / MISSING_OPTIONAL_METADATA
#   release_pride_statuses       overall PRIDE readiness (release_pride_status())
# Individual release checks report PASS / FAIL / SKIP (13_validate_release.R).
#
# Requires release_utils.R to be sourced first.

# --------------------------------------------------------------------------------------
# metadata field status vocabulary
# --------------------------------------------------------------------------------------

release_metadata_statuses <- c(
  "KNOWN_VERIFIED",
  "KNOWN_BUT_NEEDS_STANDARDIZATION",
  "MISSING_RECOVERABLE",
  "MISSING_UNKNOWN",
  "NOT_APPLICABLE"
)

release_assert_metadata_status <- function(status) {
  bad <- setdiff(unique(as.character(status)), release_metadata_statuses)
  if (length(bad)) {
    stop("Unknown metadata status value(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# --------------------------------------------------------------------------------------
# source-controlled experimenter metadata
# --------------------------------------------------------------------------------------

RELEASE_EXPERIMENTER_METADATA_COLUMNS <- c(
  "metadata_id", "category", "sample_class", "sort_order", "value", "sdrf_value",
  "status", "applies_to", "evidence_source", "evidence_date", "notes"
)

RELEASE_EXPERIMENTER_METADATA_REQUIRED_IDS <- c(
  "dataset_scope", "organism", "organism_part", "collection_method",
  "sample_class_cfos", "sample_class_mcherry", "sample_class_neuron",
  "sample_class_neuropil", "strain_or_breed", "animal_supplier",
  "animal_supplier_location", "cohort_sex_composition", "animalid_sex_assignments",
  "age_at_experiment_start", "age_at_tissue_collection", "labeling_strategy",
  "reduction_reagent", "alkylation_reagent", "cleavage_agent_lys_c",
  "cleavage_agent_trypsin", "instrument_model", "acquisition_mode",
  "search_modification_parameters", "precursor_mass_tolerance",
  "fragment_mass_tolerance", "raw_acquisition_directories"
)

RELEASE_SAMPLE_PREPARATION_COLUMNS <- c(
  "protocol_id", "protocol_version_date", "step_order", "step_name", "details",
  "temperature", "duration", "status", "evidence_source", "evidence_date", "notes"
)

release_experimenter_metadata_source_path <- function(repo_root = release_repo_root()) {
  file.path(repo_root, "07_publication_release", "metadata", "experimenter_metadata.tsv")
}

release_sample_preparation_source_path <- function(repo_root = release_repo_root()) {
  file.path(repo_root, "07_publication_release", "metadata",
            "sample_preparation_protocol.tsv")
}

release_read_checked_tsv <- function(path, required_columns, what) {
  if (!file.exists(path)) stop(what, " is missing: ", path, call. = FALSE)
  out <- utils::read.delim(path, sep = "\t", quote = "\"", stringsAsFactors = FALSE,
                           check.names = FALSE, na.strings = character(0),
                           colClasses = "character")
  missing_columns <- setdiff(required_columns, names(out))
  if (length(missing_columns)) {
    stop(what, " lacks required column(s): ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }
  out
}

release_read_experimenter_metadata <- function(repo_root = release_repo_root()) {
  path <- release_experimenter_metadata_source_path(repo_root)
  out <- release_read_checked_tsv(path, RELEASE_EXPERIMENTER_METADATA_COLUMNS,
                                  "Experimenter metadata record")
  if (any(!nzchar(out$metadata_id)) || anyDuplicated(out$metadata_id)) {
    stop("Experimenter metadata_id values must be non-empty and unique.", call. = FALSE)
  }
  missing_ids <- setdiff(RELEASE_EXPERIMENTER_METADATA_REQUIRED_IDS, out$metadata_id)
  extra_ids <- setdiff(out$metadata_id, RELEASE_EXPERIMENTER_METADATA_REQUIRED_IDS)
  if (length(missing_ids) || length(extra_ids)) {
    stop("Experimenter metadata IDs do not match the expected set. Missing: ",
         paste(missing_ids, collapse = ", "), "; unexpected: ",
         paste(extra_ids, collapse = ", "), call. = FALSE)
  }
  release_assert_metadata_status(out$status)
  classes <- out[out$category == "sample_class_definition", , drop = FALSE]
  if (!setequal(classes$sample_class, c("cfos", "mcherry", "neuron", "neuropil")) ||
      any(!nzchar(classes$sdrf_value)) || anyNA(suppressWarnings(as.integer(classes$sort_order)))) {
    stop("Experimenter metadata must define exactly the four marker-defined sample classes.",
         call. = FALSE)
  }
  out
}

release_read_sample_preparation_protocol <- function(repo_root = release_repo_root()) {
  path <- release_sample_preparation_source_path(repo_root)
  out <- release_read_checked_tsv(path, RELEASE_SAMPLE_PREPARATION_COLUMNS,
                                  "Sample-preparation protocol record")
  ord <- suppressWarnings(as.integer(out$step_order))
  if (nrow(out) != 7L || anyNA(ord) || !identical(ord, seq_len(7L)) ||
      length(unique(out$protocol_id)) != 1L ||
      !identical(unique(out$protocol_version_date), "2024-11-28") ||
      any(out$status != "KNOWN_VERIFIED")) {
    stop("Sample-preparation protocol must contain the seven ordered, verified 2024-11-28 steps.",
         call. = FALSE)
  }
  out
}

release_metadata_row <- function(metadata, metadata_id) {
  hit <- metadata[metadata$metadata_id == metadata_id, , drop = FALSE]
  if (nrow(hit) != 1L) {
    stop("Expected exactly one experimenter metadata row for ", metadata_id, ".",
         call. = FALSE)
  }
  hit
}

release_metadata_value <- function(metadata, metadata_id, standardized = FALSE) {
  hit <- release_metadata_row(metadata, metadata_id)
  value <- if (isTRUE(standardized) && nzchar(hit$sdrf_value[[1]])) {
    hit$sdrf_value[[1]]
  } else {
    hit$value[[1]]
  }
  as.character(value)
}

#' Per-field SDRF readiness, distinguishing required from optional.
release_sdrf_field_statuses <- c(
  "READY",
  "MISSING_REQUIRED_METADATA",
  "MISSING_OPTIONAL_METADATA"
)

# --------------------------------------------------------------------------------------
# PRIDE readiness
# --------------------------------------------------------------------------------------

release_pride_statuses <- c(
  "PRIDE_READY",
  "PRIDE_READY_PENDING_RAW_UPLOAD",
  "PRIDE_METADATA_INCOMPLETE",
  "PRIDE_INPUTS_INCOMPLETE"
)

#' PRIDE readiness status. Conditions are tested worst first, so a better status is
#' returned only when no worse one applies.
#'
#' @param raw_files_identified TRUE only if the acquisition filenames are known for every run.
#' @param raw_files_present TRUE only if the acquisition files themselves are in hand.
#' @param n_missing_required count of SDRF-required fields that cannot be populated.
release_pride_status <- function(raw_files_identified, raw_files_present, n_missing_required) {
  raw_files_identified <- isTRUE(raw_files_identified)
  raw_files_present <- isTRUE(raw_files_present)
  n_missing_required <- as.integer(n_missing_required)

  if (!raw_files_identified) return("PRIDE_INPUTS_INCOMPLETE")
  if (n_missing_required > 0L) return("PRIDE_METADATA_INCOMPLETE")
  if (!raw_files_present) return("PRIDE_READY_PENDING_RAW_UPLOAD")
  "PRIDE_READY"
}

# --------------------------------------------------------------------------------------
# locked-artefact checks
# --------------------------------------------------------------------------------------

#' The files whose bytes the whole animal-level generation depends on.
#'
#' Since 2026-11-24 the statistics are the per-compartment limma run written by
#' 01_preprocessing/02b_limma_stat_gct.r. The pooled ProTigy export
#' (02_data/animal_level/stat_results_for_ssGSEA_neha_proteome.gct, sha256 e1ae20f0...) is
#' superseded and no longer locked; it stays on disk for comparison.
RELEASE_LOCKED_ARTEFACTS <- list(
  animal_level_input_gct = list(
    relative = c("02_data", "animal_level", "input_gct",
                 "neha_protigy_input_animal_level_primary.gct"),
    sha256 = "f12cf99e1bfb7c17bbf56bffb6783e924698bce5d5533a8e312bc4bbb733bbb3",
    role = "animal-level L/R-averaged abundance matrix (5349 x 48)"
  ),
  stat_gct = list(
    relative = c(RELEASE_STATISTICS_RUN$stats_relative, RELEASE_STATISTICS_RUN$gct),
    sha256 = "c7f08a76edb39335e220f9ea16eaae007474c80b229b2b1f9b192528531688fd",
    role = paste("limma moderated-t statistics GCT, one model per sample class",
                 "(ProTigy file layout; written by 01_preprocessing/02b_limma_stat_gct.r)")
  ),
  stat_run_manifest = list(
    relative = c(RELEASE_STATISTICS_RUN$stats_relative, RELEASE_STATISTICS_RUN$run_manifest),
    sha256 = "c2f2871d5a8a442448c4b01b19d2e2338cfd95e0360a1757dc7a7016db4f8029",
    role = paste("run manifest of the statistics GCT: model, residual df, limma and R",
                 "versions, input and output GCT hashes")
  ),
  # The enrichment indexes record each comparison's term counts, seed and source hash, so any
  # rerun changes them; locking them keeps the exports, the figure notes and the rebuilt
  # panels on the same run. GO-BP and GO-CC: the nPermSimple = 100000 rerun (locked
  # 2026-11-25). KEGG: the GO_BP folder of the 20261124 run.
  enrichment_index_go_bp = list(
    relative = c(RELEASE_STATISTICS_RUN$enrichment_relative, "GO_BP",
                 "indexEnrichmentComparisons.csv"),
    sha256 = "8f4a81cd15d088f204c68caa53224cf7534f318a9d08e302338a17cd08d2ad5a",
    role = "index of the GO-BP enrichment run on the locked statistics (fgsea nPermSimple 100000)"
  ),
  enrichment_index_go_cc = list(
    relative = c(RELEASE_STATISTICS_RUN$enrichment_relative, "GO_CC",
                 "indexEnrichmentComparisons.csv"),
    sha256 = "8afb30d15e3f4bc02d445430e2d7f5b2ca34910cff432ba65e4fd364bbc6e451",
    role = "index of the GO-CC enrichment run on the locked statistics (fgsea nPermSimple 100000)"
  ),
  enrichment_index_kegg = list(
    relative = c(RELEASE_STATISTICS_RUN$kegg_enrichment_relative, RELEASE_STATISTICS_RUN$kegg_run,
                 "indexEnrichmentComparisons.csv"),
    sha256 = "7bd4513e393f20e52de367d8b315e689ced2b9a8ee33ed8cafd2d7d464e24875",
    role = paste("index of the enrichment run whose KEGG GSEA is released (fgsea default",
                 "nPermSimple); its GO-BP tables are superseded")
  )
)

release_locked_artefact_path <- function(key, data_root = release_data_root()) {
  spec <- RELEASE_LOCKED_ARTEFACTS[[key]]
  if (is.null(spec)) stop("Unknown locked artefact: ", key, call. = FALSE)
  do.call(file.path, c(list(data_root), as.list(spec$relative)))
}

#' Verify every locked artefact hash. Returns a data frame with one row per artefact.
release_verify_locked_artefacts <- function(data_root = release_data_root()) {
  do.call(rbind, lapply(names(RELEASE_LOCKED_ARTEFACTS), function(key) {
    spec <- RELEASE_LOCKED_ARTEFACTS[[key]]
    path <- release_locked_artefact_path(key, data_root)
    observed <- release_sha256(path)
    data.frame(
      artefact = key,
      path = path,
      role = spec$role,
      expected_sha256 = spec$sha256,
      observed_sha256 = ifelse(is.na(observed), "FILE_ABSENT", observed),
      matches = identical(observed, spec$sha256),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }))
}

# --------------------------------------------------------------------------------------
# sample-class correction
# --------------------------------------------------------------------------------------
#
# Six left-hemisphere acquisitions of C46 and C47 have a sample class that differs from
# every pre-correction record. The audit found this was a correction made during
# sample-identity QC, not a mislabelling. Only the class labels changed; acquisition
# identities and abundance values did not. Both animals follow the same cycle:
#   neuropil -> mcherry,  mcherry -> neuron,  neuron -> neuropil
#
# Stage 01 derives the six rows from the source tables and requires them to equal the record
# below; the other builders and the validator read it from here.

release_sample_class_correction_statuses <- c(
  "confirmed_intentional_switch_correction",
  "unresolved_discrepancy"
)

RELEASE_SAMPLE_CLASS_CORRECTION <- list(
  status = "confirmed_intentional_switch_correction",

  # Date of the applied metadata edit. Evidenced by the 2025-11-07 processed-matrix drop
  # under 01_input/raw_proteomics/20251107_pg.matrix_Neha/, not by a prose note -- see
  # prose_rationale_exists below.
  correction_date = "2025-11-07",

  method = "historical_sample_identity_qc_umap_correction_table",
  method_label = paste(
    "Identified during historical sample-identity QC from UMAP nearest-class-centre",
    "distances, recorded in a retained project correction table, and subsequently supported",
    "by independent proteomic profile-similarity and left/right-pair analyses."),

  # The preserved correction table, under 99_historical (a protected root).
  reference_relative = c("99_historical", "pca_plots_legacy", "tables", "umap",
                         "umap_outlier_samples_with_switches_CORRECTED.csv"),
  reference_sha256 = "76ad0ca3d27d628fedd10a3260dfff1b914fe307849bb935d782eab3e4c08e11",
  reference_mtime = "2025-11-05 15:12:56",

  # The audit compared the matrix before and after the metadata edit: no value changed.
  quantitative_values_changed = FALSE,
  max_abs_numeric_difference = 0,

  # No written rationale survives; the release says so, and the table is the record.
  prose_rationale_exists = FALSE,

  expected = data.frame(
    plate_sample_number   = c("N53", "N60", "N67", "N81", "N88", "N95"),
    AnimalID              = c("C46", "C46", "C46", "C47", "C47", "C47"),
    hemisphere            = rep("Left", 6L),
    plate_position        = c("S5-B6", "S5-C6", "S5-D6", "S5-F6", "S5-G6", "S5-H6"),
    original_sample_class = c("neuropil", "mcherry", "neuron", "neuropil", "mcherry", "neuron"),
    analysis_sample_class = c("mcherry", "neuron", "neuropil", "mcherry", "neuron", "neuropil"),
    stringsAsFactors = FALSE, check.names = FALSE
  ),

  # The cycle, stated independently of the six rows so a typo in either is caught.
  cycle = c(neuropil = "mcherry", mcherry = "neuron", neuron = "neuropil")
)

#' Number of measurements whose sample class was corrected. Exactly 6 of 96.
RELEASE_N_SAMPLE_CLASS_CORRECTED <- nrow(RELEASE_SAMPLE_CLASS_CORRECTION$expected)

release_sample_class_correction_reference_path <- function(data_root = release_data_root()) {
  do.call(file.path, c(list(data_root),
                       as.list(RELEASE_SAMPLE_CLASS_CORRECTION$reference_relative)))
}

#' Hash the preserved correction table and compare it with the recorded SHA256; stops on a
#' mismatch when require_match is TRUE. The mtime is reported but not checked.
release_verify_sample_class_correction_reference <- function(data_root = release_data_root(),
                                                             require_match = TRUE) {
  path <- release_sample_class_correction_reference_path(data_root)
  observed <- release_sha256(path)
  expected <- RELEASE_SAMPLE_CLASS_CORRECTION$reference_sha256
  mtime <- if (file.exists(path)) {
    format(as.POSIXct(file.info(path)$mtime), "%Y-%m-%d %H:%M:%S")
  } else NA_character_
  out <- list(path = path, exists = file.exists(path),
              observed_sha256 = if (is.na(observed)) "FILE_ABSENT" else observed,
              expected_sha256 = expected,
              matches = identical(observed, expected),
              observed_mtime = mtime,
              expected_mtime = RELEASE_SAMPLE_CLASS_CORRECTION$reference_mtime)
  if (isTRUE(require_match) && !out$matches) {
    stop("The preserved sample-class correction table does not match the SHA256 recorded ",
         "by the audit.\n  path:     ", path,
         "\n  expected: ", expected, "\n  observed: ", out$observed_sha256,
         "\nRefusing to publish a correction provenance record citing a file whose bytes ",
         "have changed.", call. = FALSE)
  }
  out
}

#' Compare a correction table with RELEASE_SAMPLE_CLASS_CORRECTION, field by field.
#'
#' Returns a check/pass data frame shaped like release_check_animal_level_design(), so
#' callers report both the same way.
release_check_sample_class_corrections <- function(corrections) {
  exp <- RELEASE_SAMPLE_CLASS_CORRECTION$expected
  cyc <- RELEASE_SAMPLE_CLASS_CORRECTION$cycle
  got <- as.data.frame(corrections, stringsAsFactors = FALSE, check.names = FALSE)
  ord <- if ("plate_sample_number" %in% names(got)) {
    got[order(match(got$plate_sample_number, exp$plate_sample_number)), , drop = FALSE]
  } else got
  same <- function(field) {
    field %in% names(ord) && identical(as.character(ord[[field]]), as.character(exp[[field]]))
  }
  checks <- list(
    list(sprintf("exactly %d corrected sample-class rows", nrow(exp)),
         nrow(ord) == nrow(exp)),
    list("the corrected set is exactly N53 N60 N67 N81 N88 N95",
         "plate_sample_number" %in% names(ord) &&
           setequal(as.character(ord$plate_sample_number), exp$plate_sample_number)),
    list("every corrected measurement is Left hemisphere",
         "hemisphere" %in% names(ord) && all(as.character(ord$hemisphere) == "Left")),
    list("the corrected animals are exactly C46 and C47",
         "AnimalID" %in% names(ord) &&
           setequal(as.character(ord$AnimalID), c("C46", "C47"))),
    list("plate positions match the correction record", same("plate_position")),
    list("original sample classes match the definition", same("original_sample_class")),
    list("analysis sample classes match the definition", same("analysis_sample_class")),
    list("every correction follows the neuropil->mcherry->neuron->neuropil cycle",
         all(c("original_sample_class", "analysis_sample_class") %in% names(ord)) &&
           identical(unname(cyc[as.character(ord$original_sample_class)]),
                     as.character(ord$analysis_sample_class))),
    list("no corrected row leaves the class unchanged",
         all(c("original_sample_class", "analysis_sample_class") %in% names(ord)) &&
           all(as.character(ord$original_sample_class) !=
                 as.character(ord$analysis_sample_class))),
    list("correction status is confirmed_intentional_switch_correction",
         "correction_status" %in% names(ord) &&
           all(as.character(ord$correction_status) ==
                 RELEASE_SAMPLE_CLASS_CORRECTION$status)),
    list("no corrected row claims a quantitative change",
         "quantitative_values_changed" %in% names(ord) &&
           all(toupper(as.character(ord$quantitative_values_changed)) == "FALSE"))
  )
  do.call(rbind, lapply(checks, function(ch) {
    data.frame(check = ch[[1]], pass = isTRUE(ch[[2]]),
               stringsAsFactors = FALSE, check.names = FALSE)
  }))
}

# --------------------------------------------------------------------------------------
# original submission crosswalk
# --------------------------------------------------------------------------------------
#
# The four files of the original submission are not in the project tree under their
# submitted names, so the crosswalk records how each was matched to this package:
#
#   DIRECT_VERIFICATION             the submitted file was read and hashed by this build
#   DIMENSIONS_AND_CONTENT_LINEAGE  matched by shape plus identifier containment
#   ANALYSIS_GENERATION             matched by which analysis generation produced it
#
# Their properties come from the earlier description of the submission, not from reading
# the files. If the package turns up, set PROTEOMICS_RELEASE_OLD_PACKAGE_ROOT to it and the
# build hashes each file it finds.

RELEASE_CROSSWALK_METHODS <- c(
  "DIRECT_VERIFICATION",
  "DIMENSIONS_AND_CONTENT_LINEAGE",
  "ANALYSIS_GENERATION"
)

#' Directory holding the extracted original submission package, if it is reachable.
release_old_package_root <- function() {
  release_option_or_env("proteomics.release_old_package_root",
                        "PROTEOMICS_RELEASE_OLD_PACKAGE_ROOT", "")
}

#' The four submitted files as previously described (shape, nature, replacement). Used to
#' describe them in the changelog; release_old_package_crosswalk() hashes them only when the
#' package is reachable.
RELEASE_OLD_PACKAGE_FILES <- list(
  list(name = "processed_protein_group_matrix_raw.txt",
       n_protein_rows = 5747L, n_measurement_columns = 96L, n_annotation_columns = 7L,
       nature = "processed protein-level data, NOT native raw MS data",
       replacement = paste("processed_data/protein_feature_annotation.tsv.gz (annotation)",
                           "and processed_data/protein_abundance_measurement_level.tsv.gz"),
       superseded = "Partly. Retained as the search-output layer in the lineage."),
  list(name = "processed_protein_group_matrix_filtered_umap_adjusted.xlsx",
       n_protein_rows = 5349L, n_measurement_columns = 96L, n_annotation_columns = NA_integer_,
       nature = "historical measurement/hemisphere-level representation",
       replacement = "processed_data/protein_abundance_measurement_level.tsv.gz",
       superseded = paste("No, but demoted: it is no longer the matrix inference is",
                          "drawn from.")),
  list(name = "GSEA_ORA_all_results.xlsx",
       n_protein_rows = NA_integer_, n_measurement_columns = NA_integer_,
       n_annotation_columns = NA_integer_,
       nature = paste("superseded hemisphere-level analysis generation; did not contain the",
                      "current complete 12-primary-contrast design"),
       replacement = paste("enrichment/*.tsv.gz and",
                           "editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx"),
       superseded = "Yes, entirely."),
  list(name = "README.txt",
       n_protein_rows = NA_integer_, n_measurement_columns = NA_integer_,
       n_annotation_columns = NA_integer_,
       nature = "group codes only",
       replacement = "README_DATA.md, metadata/data_dictionary.tsv, this changelog",
       superseded = "Yes.")
)

# --------------------------------------------------------------------------------------
# collection-plate wording
# --------------------------------------------------------------------------------------
# Required wording, quoted verbatim in README_DATA.md, the changelog, the manuscript action
# list and the figure source map. The design makes a collection-plate contribution
# unidentifiable in one comparison; no batch effect has been demonstrated, and the wording
# must not suggest one.
RELEASE_PLATE_SENTENCE <- paste(
  "Because pairing condition was completely associated with collection plate in the",
  "paired-vehicle versus unpaired-vehicle comparison, any collection-plate-associated",
  "contribution cannot be distinguished from a pairing-associated contribution.")

RELEASE_PLATE_NO_EFFECT_SENTENCE <- paste(
  "No downstream proteomics batch effect attributable to collection plate has been",
  "demonstrated.")

FORBIDDEN_PLATE_WORDING <- c("plate artefact", "plate artifact", "confirmed batch confound",
                             "driven by plate", "technical plate effect")

#' Statement published when the original package is not reachable. The submitted files
#' were never located, so their shapes come from the earlier description of the submission.
RELEASE_OLD_PACKAGE_EXTERNAL_STATEMENT <- paste(
  "The submitted files were not located in the project tree and were not read by this",
  "build. Their properties below are carried from the prior description of the submission",
  "package, not from a direct structural inspection; in particular",
  "GSEA_ORA_all_results.xlsx was never opened, so nothing is claimed about its internal",
  "layout. The local crosswalk rests on dimensions and content lineage.")

#' Resolve how the old-package crosswalk was established, for this build.
#'
#' Returns the method, whether the package was reachable here, and the per-file evidence.
#' When the package is reachable, each file found is hashed and marked DIRECT_VERIFICATION.
release_old_package_crosswalk <- function(root = release_old_package_root()) {
  root <- trimws(as.character(root))
  reachable <- nzchar(root) && dir.exists(root)
  if (!reachable) {
    return(list(
      method = "DIMENSIONS_AND_CONTENT_LINEAGE",
      locally_reachable = FALSE,
      root = if (nzchar(root)) root else NA_character_,
      statement = RELEASE_OLD_PACKAGE_EXTERNAL_STATEMENT,
      files = data.frame(
        file = vapply(RELEASE_OLD_PACKAGE_FILES, function(f) f$name, character(1)),
        found_locally = FALSE,
        sha256 = NA_character_,
        crosswalk_method = c("DIMENSIONS_AND_CONTENT_LINEAGE",
                             "DIMENSIONS_AND_CONTENT_LINEAGE",
                             "ANALYSIS_GENERATION", "ANALYSIS_GENERATION"),
        stringsAsFactors = FALSE, check.names = FALSE)))
  }
  found <- vapply(RELEASE_OLD_PACKAGE_FILES, function(f) {
    file.exists(file.path(root, f$name))
  }, logical(1))
  hashes <- vapply(seq_along(RELEASE_OLD_PACKAGE_FILES), function(i) {
    if (!found[[i]]) return(NA_character_)
    release_sha256(file.path(root, RELEASE_OLD_PACKAGE_FILES[[i]]$name))
  }, character(1))
  list(
    method = if (all(found)) "DIRECT_VERIFICATION" else "DIMENSIONS_AND_CONTENT_LINEAGE",
    locally_reachable = TRUE,
    root = root,
    statement = if (all(found)) {
      paste0("Original submission package read from ", root,
             " and verified file by file against this release.")
    } else {
      paste0("Original submission package partially reachable at ", root, ": ",
             sum(found), " of ", length(found), " files present. ",
             RELEASE_OLD_PACKAGE_EXTERNAL_STATEMENT)
    },
    files = data.frame(
      file = vapply(RELEASE_OLD_PACKAGE_FILES, function(f) f$name, character(1)),
      found_locally = found,
      sha256 = hashes,
      crosswalk_method = ifelse(found, "DIRECT_VERIFICATION",
                                "DIMENSIONS_AND_CONTENT_LINEAGE"),
      stringsAsFactors = FALSE, check.names = FALSE))
}

# --------------------------------------------------------------------------------------
# effect-size semantics
# --------------------------------------------------------------------------------------
#
# The analysis matrix is standardised per protein across the measurement-level dataset, so
# the coefficient stored as `logFC` (limma; `log2fc` in the split tables) is a difference on
# that standardised scale, not a log2 fold change.
#
# Two terms are avoided: "standardised mean difference" without qualification (it reads as
# Cohen's d, which divides by a pooled within-group SD, not the per-protein SD over the whole
# dataset used here), and "z-scored" (the released values have 2 decimals, so row means and
# SDs are 0 and 1 only to rounding, and the producing step is unresolved; see
# release_standardization_evidence()).

RELEASE_EFFECT_SIZE <- list(
  public_field = "effect_size_sd_units",
  public_term = "standardized abundance difference (SD units)",
  public_term_alt = "difference in standardized protein abundance",
  source_field = "logFC",
  source_field_detail = "logFC (statistics GCT, limma) / log2fc (canonical split tables)",
  units = "per-protein standard deviation of the standardised abundance scale",

  definition = paste(
    "Difference of group means (numerator minus denominator) on the protein-standardised",
    "abundance scale -- a standardized abundance difference expressed in SD units of that",
    "scale. It is NOT a log2 fold change, despite being stored as `logFC` in the statistics",
    "GCT (limma's coefficient name) and `log2fc` in the canonical split tables. Positive",
    "means higher in the numerator condition. A value of 2 means two standard deviations on",
    "the standardised scale, not a four-fold change."),

  methods_sentence = paste(
    "Protein abundances in the analysis matrix were standardized separately for each",
    "protein across the measurement-level dataset. Consequently, the coefficient stored by",
    "the differential-analysis software (limma) as `logFC` represents the difference in",
    "standardized protein abundance between groups, expressed on that standardized scale,",
    "rather than a log2 fold change."),

  sensitivity_public_label = "effect-size-ranked sensitivity analysis",
  sensitivity_public_detail = paste(
    "GSEA ranked by the standardised-abundance effect size instead of the moderated t.",
    "Sensitivity analysis only; the canonical ranking statistic is the moderated t."),

  # File names and column values that keep the old token so they still match the canonical
  # run.
  retained_internal_tokens = c("log2fc", "logFC", "GSEA_log2FC_sensitivity",
                               "GSEA_GO_BP_log2fc_sensitivity",
                               "GSEA_GO_CC_log2fc_sensitivity",
                               "GSEA_KEGG_log2fc_sensitivity",
                               "ORA_GO_BP_top_abs_log2fc", "ORA_GO_CC_top_abs_log2fc",
                               "top_absolute_log2fc")
)

#' Phrases that would assert a fold-change reading of the effect size.
RELEASE_FOLD_CHANGE_CLAIMS <- c(
  "log2 fold change", "log2 fold-change", "log2-fold change", "log2-fold-change",
  "log 2 fold change", "log2 ratio", "fold change", "fold-change"
)

#' Words that make a fold-change mention acceptable when they are on the same line: a
#' negation ("is not a log2 fold change"), a source or history marker ("stored as logFC",
#' "historical"), or `ewce`, whose bootstrap statistic is a real fold change.
RELEASE_FOLD_CHANGE_EXEMPT_MARKERS <- c(
  "not a", "not the", "is not", "are not", "was not", "were not", "never",
  "rather than", "instead of", "despite", "no longer", "misleading", "mislabel",
  "incorrect", "incorrectly", "historical", "historically", "legacy", "superseded",
  "pre-correction",
  "source column", "source field", "stored as", "stored by", "named", "called",
  "do not", "does not", "cannot", "would assert", "not supported", "ewce"
)

#' Lines that describe published effect sizes as fold changes with no negation in sight.
#'
#' Returns a data frame of offending lines; empty means the text is clean.
release_fold_change_mislabels <- function(lines,
                                          claims = RELEASE_FOLD_CHANGE_CLAIMS,
                                          exempt = RELEASE_FOLD_CHANGE_EXEMPT_MARKERS) {
  lines <- as.character(lines)
  empty <- data.frame(line = integer(0), term = character(0), text = character(0),
                      stringsAsFactors = FALSE, check.names = FALSE)
  if (!length(lines)) return(empty)
  low <- tolower(lines)
  hits <- list()
  for (i in seq_along(lines)) {
    matched <- claims[vapply(claims, function(cl) grepl(cl, low[[i]], fixed = TRUE),
                             logical(1))]
    if (!length(matched)) next
    if (any(vapply(exempt, function(mk) grepl(mk, low[[i]], fixed = TRUE), logical(1)))) next
    hits[[length(hits) + 1L]] <- data.frame(
      line = i, term = matched[[1]], text = trimws(lines[[i]]),
      stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!length(hits)) return(empty)
  do.call(rbind, hits)
}

#' Row means and SDs (both conventions) of a matrix, and whether its rows are approximately
#' standardised or z-scored to full precision. On the released 2-decimal values
#' `exact_zscore` is FALSE; callers use it to choose their wording.
release_standardization_evidence <- function(mat, tol_exact = 1e-9) {
  row_mean <- rowMeans(mat)
  row_sd_n1 <- apply(mat, 1L, stats::sd)
  row_sd_n <- sqrt(rowMeans((mat - row_mean)^2))
  list(
    n_proteins = nrow(mat),
    n_columns = ncol(mat),
    max_abs_row_mean = max(abs(row_mean)),
    max_abs_row_sd_minus_1 = max(abs(row_sd_n1 - 1)),
    max_abs_row_sd_pop_minus_1 = max(abs(row_sd_n - 1)),
    median_row_sd = stats::median(row_sd_n1),
    approximately_standardized = max(abs(row_mean)) < 1e-2 &&
      max(abs(row_sd_n1 - 1)) < 5e-2,
    exact_zscore = max(abs(row_mean)) < tol_exact &&
      max(abs(row_sd_n1 - 1)) < tol_exact
  )
}

# --------------------------------------------------------------------------------------
# design invariants
# --------------------------------------------------------------------------------------

RELEASE_DESIGN_INVARIANTS <- list(
  n_measurement_records = 96L,
  n_animal_level_units = 48L,
  n_animals = 12L,
  n_sample_classes = 4L,
  n_conditions = 4L,
  n_strata = 16L,
  n_animals_per_stratum = 3L,
  n_primary_contrasts = 12L,
  n_proteins_statistical = 5349L,
  n_proteins_mapped = 5327L,
  n_hemispheres_per_unit = 2L
)

# --------------------------------------------------------------------------------------
# the primary statistical model
# --------------------------------------------------------------------------------------
# One definition of the model behind every released t, P and FDR. Builders and prose read
# the label, the residual df and the script from here; the limma and R versions come from
# the run manifest (release_read_stat_run()), never from the machine running the build.
# release_read_stat_run() asserts the numbers below against that manifest.

RELEASE_PRIMARY_MODEL <- local({
  inv <- RELEASE_DESIGN_INVARIANTS
  n_cols <- inv$n_conditions * inv$n_animals_per_stratum
  df_res <- n_cols - inv$n_conditions
  list(
    id = RELEASE_STATISTICS_RUN$id,
    date = RELEASE_STATISTICS_RUN$date,
    software = "limma",
    script = "01_preprocessing/02b_limma_stat_gct.r",
    variance_model = "within_compartment",
    n_columns_in_model = n_cols,
    residual_df = df_res,
    label = paste0(
      "limma ~0 + condition fitted separately within each sample class (", n_cols,
      " animals, n = ", inv$n_animals_per_stratum, " per condition, residual df ", df_res,
      "); eBayes(robust = TRUE, trend = FALSE); Benjamini-Hochberg within each contrast"),
    methods_sentence = paste0(
      "Differential abundance was tested with limma, fitting one linear model per sample",
      " class (~0 + condition over that class's ", n_cols, " animal-level units, n = ",
      inv$n_animals_per_stratum, " animals per condition, residual df ", df_res,
      "), with robust empirical Bayes moderation without an intensity trend and",
      " Benjamini-Hochberg adjustment within each of the ", inv$n_primary_contrasts,
      " comparisons."),
    superseded_id = "pooled_protigy_20260824",
    superseded_label = paste0(
      "pooled ProTigy export: one limma fit over all ", inv$n_animal_level_units,
      " animal-level columns (", inv$n_strata, " sample_class x condition groups, residual df ",
      inv$n_animal_level_units - inv$n_strata, "); eBayes(robust = TRUE, trend = FALSE)")
  )
})

#' Read the statistics run manifest. Stops if its hash, model, columns, residual df or GCT
#' hashes differ from RELEASE_PRIMARY_MODEL and the locked files.
release_read_stat_run <- function(data_root = release_data_root()) {
  path <- release_stat_run_path(RELEASE_STATISTICS_RUN$run_manifest, data_root)
  if (!file.exists(path)) stop("Statistics run manifest missing: ", path, call. = FALSE)
  if (!identical(release_sha256(path), RELEASE_LOCKED_ARTEFACTS$stat_run_manifest$sha256)) {
    stop("Statistics run manifest no longer matches its locked hash: ", path, call. = FALSE)
  }
  m <- release_read_csv(path)
  one <- function(col) {
    v <- unique(as.character(m[[col]]))
    if (length(v) != 1L || is.na(v) || !nzchar(v)) {
      stop("Run manifest column ", col, " is not one value across contrasts: ", path,
           call. = FALSE)
    }
    v
  }
  pm <- RELEASE_PRIMARY_MODEL
  if (!identical(one("variance_model"), pm$variance_model) ||
      !identical(as.integer(one("n_columns_in_model")), as.integer(pm$n_columns_in_model)) ||
      !identical(as.integer(one("df_residual")), as.integer(pm$residual_df)) ||
      !identical(one("output_gct_sha256"), RELEASE_LOCKED_ARTEFACTS$stat_gct$sha256) ||
      !identical(one("input_gct_sha256"), RELEASE_LOCKED_ARTEFACTS$animal_level_input_gct$sha256)) {
    stop("The statistics run manifest does not describe RELEASE_PRIMARY_MODEL on the ",
         "locked GCTs: ", path, call. = FALSE)
  }
  r_string <- one("r_version")
  list(path = path, manifest = m,
       limma_version = one("limma_version"),
       r_version_string = r_string,
       r_version = sub("^R version ([0-9.]+).*$", "\\1", r_string),
       n_columns_in_model = as.integer(one("n_columns_in_model")),
       df_residual = as.integer(one("df_residual")),
       input_gct_sha256 = one("input_gct_sha256"),
       output_gct_sha256 = one("output_gct_sha256"),
       run_timestamp_utc = one("run_timestamp_utc"))
}

#' AveExpr as the per-compartment model reports it, measured on the released table.
#'
#' Each sample class has its own model, so AveExpr is that class's mean over its 12 units
#' and is not centred at 0 within a class. The mean of the four class values is the mean over
#' all 48 units, which is what corroborates a per-protein centred input. Prose quotes these
#' numbers instead of typing them. AveExpr is identical across a class's three contrasts, so
#' one contrast family gives one value per protein and class.
release_aveexpr_summary <- function(differential) {
  one_family <- differential[differential$contrast_family == "paired_veh_vs_unpaired_veh", ,
                             drop = FALSE]
  ave <- as.numeric(one_family$average_standardized_abundance)
  class_mean <- tapply(ave, one_family$protein_group_id, mean)
  per_class <- tapply(ave, one_family$sample_class, range)
  list(n_proteins = length(class_mean),
       max_abs_mean_over_classes = max(abs(class_mean)),
       range_within_classes = range(ave),
       range_by_class = per_class)
}

# --------------------------------------------------------------------------------------
# enrichment notes carried verbatim by more than one artefact
# --------------------------------------------------------------------------------------

RELEASE_GO_CC_PROVENANCE_NOTE <- paste(
  "GO-CC restored 2026-11-24 (it was reported in the 2025 analysis and dropped in",
  "August 2026).")

RELEASE_ENRICHMENT_SEED_HANDLING <- paste(
  "Since 2026-11-24 every GSEA call passes DOSE seed = FALSE, so the per-comparison seed",
  "derived from gsea_seed_base is the seed actually used. The pooled run used seed = TRUE,",
  "which reseeds every call with 10407 (the RNG-kind code) and ignores the recorded seed.",
  "Term counts near q = 0.05 move by a few terms with the permutation seed.")

RELEASE_KEGG_VERSION_NOTE <- paste(
  "KEGG gene sets are downloaded at run time (KEGGREST) and are not version-pinned.")

RELEASE_EWCE_RELATIVE <- c("03_output", "ewce", "EWCE_Results_animal_level_validation_20260825")

#' The limma fit EWCE ran for its differential protein lists, as its own audit recorded it.
#'
#' EWCE (2026-08-25) fits its own per-sample-class limma model on the proteins present in the
#' cell-type reference, with eBayes(trend = TRUE). Its top-N lists therefore rank by that
#' fit's t, not by the released statistics. Read from the audit, never typed.
release_ewce_differential_model <- function(data_root = release_data_root()) {
  path <- file.path(release_rel_path(RELEASE_EWCE_RELATIVE, data_root), "03_QC_Mapping_Logs",
                    "animal_level_differential_audit.csv")
  a <- release_read_csv(path)
  one <- function(col) {
    v <- unique(as.character(a[[col]]))
    if (length(v) != 1L || is.na(v) || !nzchar(v)) {
      stop("EWCE differential audit column ", col, " is not one value: ", path, call. = FALSE)
    }
    v
  }
  design <- one("limma_design")
  ebayes <- one("limma_ebayes")
  n_proteins <- as.integer(one("n_proteins_entering_differential_analysis"))
  list(path = path, design = design, ebayes = ebayes, n_proteins = n_proteins,
       n_comparisons = nrow(a),
       sentence = paste0(
         "EWCE differential lists come from EWCE's own limma fit per sample class (",
         sub(";.*$", "", design), " over the ", format(n_proteins, big.mark = ","),
         " proteins that map to its cell-type reference; eBayes ", gsub(";", ", ", ebayes),
         "), not from the released statistics, so they rank by that fit's t rather than by",
         " the released moderated t."))
}

#' Count of the fgsea "P values not calculated" sets in one warnings.csv, per GSEA analysis.
release_fgsea_warning_counts <- function(path) {
  if (!file.exists(path)) stop("Enrichment warnings file missing: ", path, call. = FALSE)
  w <- release_read_csv(path)
  hit <- if (nrow(w)) grepl("pathways for which P-values were not calculated properly",
                            w$warning, fixed = TRUE) else logical(0)
  out <- data.frame(analysis = as.character(w$analysis[hit]),
                    n_sets = suppressWarnings(as.integer(sub(
                      "^There were ([0-9]+) pathways.*$", "\\1", w$warning[hit]))),
                    stringsAsFactors = FALSE)
  if (anyNA(out$n_sets) || anyDuplicated(out$analysis)) {
    stop("Could not read the fgsea untested-set counts from ", path, call. = FALSE)
  }
  out
}

#' Gene sets fgsea returned without a P value, one row per set (GO-BP and GO-CC).
#'
#' fgsea leaves P undefined for a set when the ranked statistics are too one-sided for its
#' simple permutations ("unbalanced" statistics), and DOSE drops such rows, so they are absent
#' from the GSEA tables rather than non-significant. The GO runs list them per comparison in
#' audits/<comparison>/gsea_untested_sets.csv, for the canonical and the effect-size-ranked
#' GSEA; each list must agree with the run's warnings.csv and its index. The KEGG run predates
#' that audit, so KEGG has counts only (release_fgsea_untested()).
release_fgsea_untested_sets <- function(comparisons, data_root = release_data_root()) {
  rows <- list()
  for (run in RELEASE_STATISTICS_RUN$enrichment_runs) {
    root <- release_enrichment_root(run, data_root)
    index <- release_read_csv(file.path(root, "indexEnrichmentComparisons.csv"))
    for (cmp in comparisons) {
      p <- file.path(root, "audits", cmp, "gsea_untested_sets.csv")
      if (!file.exists(p)) stop("Untested-set audit missing: ", p, call. = FALSE)
      d <- release_read_csv(p)
      counts <- release_fgsea_warning_counts(file.path(root, "audits", cmp, "warnings.csv"))
      listed <- table(factor(d$analysis, levels = unique(c(d$analysis, counts$analysis))))
      warned <- stats::setNames(counts$n_sets, counts$analysis)
      ix <- index[index$canonical_comparison == cmp, , drop = FALSE]
      recorded <- c(suppressWarnings(as.integer(ix$go_gsea_untested_set_count)),
                    suppressWarnings(as.integer(ix$sensitivity_go_gsea_untested_set_count)))
      canonical <- paste0(run, "_GSEA")
      if (nrow(ix) != 1L || anyNA(recorded) ||
          !setequal(names(listed)[listed > 0], names(warned)) ||
          !all(as.integer(listed[names(warned)]) == unname(warned)) ||
          recorded[[1]] != sum(d$analysis == canonical) ||
          recorded[[2]] != sum(d$analysis == paste0(canonical, "_log2fc_sensitivity"))) {
        stop("The untested-set list, warnings.csv and the index disagree for ", run, " ", cmp,
             call. = FALSE)
      }
      if (nrow(d)) rows[[length(rows) + 1L]] <- d
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else
    data.frame(canonical_comparison = character(0), analysis = character(0), ID = character(0),
               Description = character(0), setSize = integer(0), stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' Gene sets fgsea returned without a P value, counted per comparison and collection.
#'
#' Canonical GSEA only. GO-BP and GO-CC are counted from the GO runs' untested-set lists, KEGG
#' from the KEGG run's audits/<comparison>/warnings.csv.
release_fgsea_untested <- function(comparisons, data_root = release_data_root()) {
  sets <- release_fgsea_untested_sets(comparisons, data_root)
  sets <- sets[sets$analysis %in% paste0(RELEASE_STATISTICS_RUN$enrichment_runs, "_GSEA"), ,
               drop = FALSE]
  go <- if (nrow(sets)) {
    agg <- stats::aggregate(list(n_sets = rep(1L, nrow(sets))),
                            sets[c("canonical_comparison", "analysis")], sum)
    agg$ontology <- sub("_GSEA$", "", agg$analysis)
    agg[c("canonical_comparison", "analysis", "ontology", "n_sets")]
  } else NULL
  kegg <- do.call(rbind, lapply(comparisons, function(cmp) {
    w <- release_fgsea_warning_counts(file.path(release_kegg_root(data_root), "audits", cmp,
                                                "warnings.csv"))
    w <- w[w$analysis == "KEGG_GSEA", , drop = FALSE]
    if (!nrow(w)) return(NULL)
    data.frame(canonical_comparison = cmp, analysis = w$analysis, ontology = "KEGG",
               n_sets = w$n_sets, stringsAsFactors = FALSE)
  }))
  out <- rbind(go, kegg)
  if (is.null(out)) {
    out <- data.frame(canonical_comparison = character(0), analysis = character(0),
                      ontology = character(0), n_sets = integer(0), stringsAsFactors = FALSE)
  }
  out$n_sets <- as.integer(out$n_sets)
  out <- out[order(match(out$canonical_comparison, comparisons), out$analysis), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' How the GSEA permutations were run.
#'
#' GO-BP and GO-CC: the adopted rerun, whose run_parameters.csv records nPermSimple. KEGG and
#' the sensitivity analyses (Pairing x CNO interaction, leave-one-animal-out, t- vs
#' effect-size ranking) record none, so fgsea used its default; the value is quoted only when
#' the installed fgsea is the version those runs recorded. The rank-robustness run must name
#' the KEGG run as its primary reference, which is what makes the comparison same-setting.
#' Returns the facts and the sentences every reader-facing document quotes.
release_enrichment_permutations <- function(data_root = release_data_root()) {
  if (!identical(RELEASE_STATISTICS_RUN$enrichment_rerun_status, "adopted")) {
    stop("RELEASE_STATISTICS_RUN$enrichment_rerun_status must record the decision.",
         call. = FALSE)
  }
  n_perm <- function(root) {
    p <- release_read_csv(file.path(root, "audits", "run_parameters.csv"))
    p$value[p$parameter == "gsea_n_perm_simple"]
  }
  fmt_n <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE)
  go_n <- unique(unlist(lapply(RELEASE_STATISTICS_RUN$enrichment_runs, function(run) {
    n_perm(release_enrichment_root(run, data_root))
  })))
  if (length(go_n) != 1L || is.na(suppressWarnings(as.numeric(go_n)))) {
    stop("The GO-BP and GO-CC enrichment runs do not record one nPermSimple.", call. = FALSE)
  }
  if (length(n_perm(release_kegg_root(data_root)))) {
    stop("The KEGG run records an nPermSimple; the default-setting wording would be wrong.",
         call. = FALSE)
  }
  versions <- release_read_csv(file.path(release_kegg_root(data_root), "audits",
                                         "package_database_versions.csv"))
  fgsea_run <- as.character(versions$version[versions$component == "fgsea"])

  # The sensitivity analyses: no nPermSimple recorded, same fgsea, primary reference = KEGG run.
  inferential <- release_inferential_root(data_root)
  manifests <- list(
    interaction = release_read_csv(file.path(inferential,
                                             "pairing_cno_interaction_run_manifest.csv")),
    robustness = release_read_csv(file.path(inferential, "robustness",
                                            "robustness_run_manifest.csv")))
  for (nm in names(manifests)) {
    m <- manifests[[nm]]
    if (any(grepl("n_?perm", m$name, ignore.case = TRUE)) ||
        !identical(as.character(m$value[m$section == "software" & m$name == "fgsea"]), fgsea_run)) {
      stop("The ", nm, " run does not match the default-setting fgsea runs.", call. = FALSE)
    }
  }
  rob_root <- as.character(manifests$robustness$value[manifests$robustness$section == "parameter" &
                                                        manifests$robustness$name == "enrichment_root"])
  if (length(rob_root) != 1L ||
      !identical(tolower(project_normalize_path(rob_root)),
                 tolower(project_normalize_path(release_rel_path(
                   RELEASE_STATISTICS_RUN$kegg_enrichment_relative, data_root))))) {
    stop("The rank-robustness run does not compare with the default-setting primary run.",
         call. = FALSE)
  }
  default_n <- NA_integer_
  if (length(fgsea_run) == 1L && requireNamespace("fgsea", quietly = TRUE) &&
      identical(as.character(utils::packageVersion("fgsea")), fgsea_run)) {
    default_n <- suppressWarnings(as.integer(eval(formals(fgsea::fgseaMultilevel)$nPermSimple)))
  }
  default_text <- paste0("fgsea's default nPermSimple",
                         if (!is.na(default_n)) paste0(" (", fmt_n(default_n), " in fgsea ",
                                                       fgsea_run, ")") else "")
  go_rel <- paste(RELEASE_STATISTICS_RUN$enrichment_relative, collapse = "/")
  kegg_rel <- paste(RELEASE_STATISTICS_RUN$kegg_enrichment_relative, collapse = "/")
  default_short <- paste0("fgsea default",
                          if (!is.na(default_n)) paste0(" (", fmt_n(default_n), " in fgsea ",
                                                        fgsea_run, ")") else "")
  setting <- paste0("GO-BP and GO-CC ", fmt_n(go_n), "; KEGG ", default_short)
  sensitivities <- paste("Pairing x CNO interaction, leave-one-animal-out and t- versus",
                         "effect-size ranking for mCherry learning")
  same_setting <- paste0("compared with the primary analysis at that same setting (", kegg_rel,
                         "), not with the nPermSimple = ", fmt_n(go_n), " tables")
  # Standalone, for the rows and paragraphs that describe one sensitivity analysis.
  sensitivity_sentence <- paste0("The sensitivity analyses (", sensitivities, ") use ",
                                 default_text, " and are ", same_setting, ".")
  sentence <- paste0(
    "GO-BP and GO-CC GSEA P values come from fgsea with nPermSimple = ", fmt_n(go_n), " (",
    go_rel, "). KEGG GSEA (", kegg_rel, ") and the sensitivity analyses (", sensitivities,
    ") use ", default_text, "; the sensitivity analyses are ", same_setting, ".")
  list(setting = setting, go_n_perm_simple = as.integer(go_n),
       default_n_perm_simple = default_n, fgsea_version = fgsea_run,
       go_relative = go_rel, kegg_relative = kegg_rel,
       sensitivity_sentence = sensitivity_sentence, sentence = sentence)
}

#' Assert the animal-level design invariants against a 48-row animal-level table.
release_check_animal_level_design <- function(animal_level) {
  inv <- RELEASE_DESIGN_INVARIANTS
  checks <- list(
    list("animal-level table has 48 rows",
         nrow(animal_level) == inv$n_animal_level_units),
    list("12 distinct AnimalIDs",
         length(unique(animal_level$AnimalID)) == inv$n_animals),
    list("4 sample classes",
         length(unique(animal_level$sample_class)) == inv$n_sample_classes),
    list("4 conditions",
         length(unique(animal_level$condition)) == inv$n_conditions),
    list("16 sample_class x condition strata",
         nrow(unique(animal_level[c("sample_class", "condition")])) == inv$n_strata),
    list("3 animals in every stratum",
         all(table(animal_level$sample_class, animal_level$condition) == inv$n_animals_per_stratum)),
    list("no duplicated AnimalID x sample_class",
         !anyDuplicated(paste(animal_level$AnimalID, animal_level$sample_class))),
    list("2 hemisphere measurements per animal-level unit",
         all(animal_level$n_hemisphere_measurements == inv$n_hemispheres_per_unit)),
    list("both hemispheres present in every unit",
         all(animal_level$hemispheres_present == "Left;Right"))
  )
  do.call(rbind, lapply(checks, function(ch) {
    data.frame(check = ch[[1]], pass = isTRUE(ch[[2]]),
               stringsAsFactors = FALSE, check.names = FALSE)
  }))
}

# --------------------------------------------------------------------------------------
# prohibited-computation scanning
# --------------------------------------------------------------------------------------

#' Function calls the release layer must not make (model fitting, enrichment, EWCE, PCA,
#' tests and P-value adjustment, random draws). The scan uses the parse tree, so a name
#' inside a string or comment (stage 10 quotes `rnorm(n, mean = ...)`) is not a hit.
RELEASE_PROHIBITED_CALLS <- c(
  "lmFit", "eBayes", "contrasts.fit", "makeContrasts", "model.matrix",
  "duplicateCorrelation", "treat", "topTable", "decideTests",
  "gseGO", "gseKEGG", "enrichGO", "enrichKEGG", "GSEA", "enricher",
  "fgsea", "fgseaMultilevel", "fgseaSimple",
  "bootstrap_enrichment_test", "generate_celltype_data", "ewce_expression_data",
  "prcomp", "princomp", "svd", "irlba",
  "p.adjust", "t.test", "aov", "anova", "lm", "glm", "wilcox.test",
  "rnorm", "sample", "runif"
)

#' Scan R sources for prohibited computation. Returns one row per hit.
#'
#' Uses the R parser (SYMBOL_FUNCTION_CALL tokens). A file that does not parse falls back
#' to a text scan, and its hits are marked as such.
release_scan_prohibited_calls <- function(paths, prohibited = RELEASE_PROHIBITED_CALLS) {
  empty <- data.frame(file = character(0), line = integer(0), call = character(0),
                      text = character(0), detected_by = character(0),
                      stringsAsFactors = FALSE, check.names = FALSE)
  hits <- list()
  for (p in paths) {
    if (!file.exists(p)) next
    lines <- readLines(p, warn = FALSE)
    parsed <- tryCatch(parse(p, keep.source = TRUE), error = function(e) NULL)
    if (!is.null(parsed)) {
      pd <- utils::getParseData(parsed)
      if (!is.null(pd) && nrow(pd)) {
        calls <- pd[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text %in% prohibited, ,
                    drop = FALSE]
        for (i in seq_len(nrow(calls))) {
          ln <- calls$line1[i]
          hits[[length(hits) + 1L]] <- data.frame(
            file = p, line = as.integer(ln), call = calls$text[i],
            text = trimws(if (ln <= length(lines)) lines[[ln]] else ""),
            detected_by = "parse tree",
            stringsAsFactors = FALSE, check.names = FALSE)
        }
        next
      }
    }
    code_mask <- !grepl("^\\s*#", lines)
    for (call_name in prohibited) {
      pattern <- paste0("(^|[^A-Za-z0-9._])", gsub("\\.", "[.]", call_name), "\\s*\\(")
      for (i in which(code_mask & grepl(pattern, lines))) {
        hits[[length(hits) + 1L]] <- data.frame(
          file = p, line = i, call = call_name, text = trimws(lines[[i]]),
          detected_by = "text fallback (file did not parse)",
          stringsAsFactors = FALSE, check.names = FALSE)
      }
    }
  }
  if (!length(hits)) return(empty)
  do.call(rbind, hits)
}

# --------------------------------------------------------------------------------------
# workbook hygiene
# --------------------------------------------------------------------------------------

#' Column-name problems per sheet: no columns, duplicated names, or blank and placeholder
#' names.
release_check_workbook_hygiene <- function(sheet_columns) {
  problems <- character()
  for (sheet in names(sheet_columns)) {
    cols <- as.character(sheet_columns[[sheet]])
    if (!length(cols)) {
      problems <- c(problems, sprintf("sheet '%s' has no columns", sheet))
      next
    }
    if (anyDuplicated(cols)) {
      problems <- c(problems, sprintf("sheet '%s' has duplicate column names: %s", sheet,
                                      paste(unique(cols[duplicated(cols)]), collapse = ", ")))
    }
    unnamed <- cols[is.na(cols) | !nzchar(trimws(cols)) | grepl("^(\\.\\.\\.|V|X)[0-9]+$", cols)]
    if (length(unnamed)) {
      problems <- c(problems, sprintf("sheet '%s' has unnamed/placeholder columns: %s", sheet,
                                      paste(unnamed, collapse = ", ")))
    }
  }
  problems
}

#' TRUE for values that look like UniProt accessions or entry names. Used to check that no
#' published gene_symbol column holds them (in the pipeline's split and mapped tables,
#' `gene_symbol` holds UniProt identifiers).
release_uniprot_like <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(logical(0))
  accession <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", x) |
    grepl("^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", x)
  entry_name <- grepl("_(MOUSE|HUMAN|PIG|BOVIN|RAT)$", x)
  accession | entry_name
}

#' TRUE when a column named as a gene symbol actually contains UniProt identifiers.
release_column_is_misleading_gene_symbol <- function(values, threshold = 0.5) {
  flags <- release_uniprot_like(values)
  if (!length(flags)) return(FALSE)
  mean(flags) >= threshold
}
