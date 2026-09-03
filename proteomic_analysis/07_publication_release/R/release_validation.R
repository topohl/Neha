# Validation vocabulary and shared checks for the publication / deposition release.
#
# Two separate vocabularies live here and must not be conflated:
#
#   * METADATA FIELD STATUS -- what is known about one metadata field
#     (release_metadata_statuses). Applied per SDRF/metadata field.
#   * RELEASE CHECK STATUS  -- PASS / FAIL / SKIP for one validation contract.
#
# and one derived verdict:
#
#   * PRIDE READINESS -- one of four values, chosen from evidence, never optimistic.
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
  "fragment_mass_tolerance", "raw_acquisition_directories",
  "excluded_rat_viral_insertion_experiment"
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
  out <- utils::read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE,
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
    stop("Experimenter metadata ID contract mismatch. Missing: ",
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
  excluded <- out[out$metadata_id == "excluded_rat_viral_insertion_experiment", , drop = FALSE]
  if (nrow(excluded) != 1L || excluded$status != "NOT_APPLICABLE" ||
      !grepl("Sprague-Dawley rat", excluded$value, fixed = TRUE)) {
    stop("The unrelated rat-experiment scope exclusion is absent or malformed.",
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

#' Derive the single PRIDE readiness verdict from evidence.
#'
#' Deliberately ordered worst-first so no optimistic branch can be reached while a
#' pessimistic condition holds.
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
# locked-artefact contracts
# --------------------------------------------------------------------------------------

#' The two GCTs whose bytes the whole animal-level generation depends on.
RELEASE_LOCKED_ARTEFACTS <- list(
  animal_level_input_gct = list(
    relative = c("02_data", "animal_level", "input_gct",
                 "neha_protigy_input_animal_level_primary.gct"),
    sha256 = "f12cf99e1bfb7c17bbf56bffb6783e924698bce5d5533a8e312bc4bbb733bbb3",
    role = "animal-level L/R-averaged abundance matrix (5349 x 48)"
  ),
  protigy_stat_gct = list(
    relative = c("02_data", "animal_level", "stat_results_for_ssGSEA_neha_proteome.gct"),
    sha256 = "e1ae20f02e2747cfae3572933f2b23e6c770b92ef6810963a2806afb7adbe2b6",
    role = "ProTigy two-sample moderated-t statistical results GCT"
  )
)

release_locked_artefact_path <- function(key, data_root = release_data_root()) {
  spec <- RELEASE_LOCKED_ARTEFACTS[[key]]
  if (is.null(spec)) stop("Unknown locked artefact: ", key, call. = FALSE)
  do.call(file.path, c(list(data_root), as.list(spec$relative)))
}

#' Verify both locked GCT hashes. Returns a data frame; never silently passes.
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
# sample-class correction contract
# --------------------------------------------------------------------------------------
#
# Six acquisitions from the LEFT hemispheres of C46 and C47 carry a sample-class assignment
# that differs from every pre-correction record. The forensic audit resolved this: it was a
# deliberate correction applied during historical sample-identity QC, not a mislabelling
# and not an open question. Only class METADATA was reassigned -- acquisition identities and
# quantitative abundance profiles are untouched.
#
# The reassignment is the same 3-cycle for both animals:
#   neuropil -> mcherry,  mcherry -> neuron,  neuron -> neuropil
#
# This constant is the contract, not a description of one. Stage 01 derives the correction
# from the canonical sources independently and then requires the result to equal what is
# written here; the builders and the validator all read this single definition, so the six
# rows cannot drift apart between the metadata table, the provenance table and the prose.

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

  # The preserved historical correction record. READ-ONLY: it lives under 99_historical,
  # which this layer must never write into.
  reference_relative = c("99_historical", "pca_plots_legacy", "tables", "umap",
                         "umap_outlier_samples_with_switches_CORRECTED.csv"),
  reference_sha256 = "76ad0ca3d27d628fedd10a3260dfff1b914fe307849bb935d782eab3e4c08e11",
  reference_mtime = "2025-11-05 15:12:56",

  # Established by the forensic audit: the quantitative matrix was compared across the
  # metadata edit and the maximum absolute numeric difference was exactly 0.
  quantitative_values_changed = FALSE,
  max_abs_numeric_difference = 0,

  # There is NO surviving prose note recording the rationale. That absence is published as
  # a fact rather than papered over; the correction table itself is the preserved record.
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

#' Verify the preserved correction table against the hash the forensic audit recorded.
#'
#' The hash is computed here and REQUIRED to match; it is never taken on trust from the
#' audit narrative. mtime is reported alongside but is deliberately not the check -- a
#' modification time is not an integrity guarantee.
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
         "by the forensic audit.\n  path:     ", path,
         "\n  expected: ", expected, "\n  observed: ", out$observed_sha256,
         "\nRefusing to publish a correction provenance record citing a file whose bytes ",
         "have changed.", call. = FALSE)
  }
  out
}

#' Assert a correction table equals the contract on every field that matters.
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
    list("original sample classes match the contract", same("original_sample_class")),
    list("analysis sample classes match the contract", same("analysis_sample_class")),
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
# old journal package crosswalk
# --------------------------------------------------------------------------------------
#
# The four files of the original journal submission do not exist under their submitted
# names anywhere in the project tree. The crosswalk from them to this package is therefore
# established by evidence, and the STRENGTH of that evidence is reported rather than
# glossed:
#
#   DIRECT_VERIFICATION             the submitted file was read here, byte-for-byte
#   DIMENSIONS_AND_CONTENT_LINEAGE  matched by shape plus identifier containment
#   ANALYSIS_GENERATION             matched by which analysis generation produced it
#
# The original ZIP has been inspected separately (outside this coding environment). That
# inspection is reported as what it is -- an external verification -- and is deliberately
# NOT laundered into a claim that this build read those bytes. If the package is later
# mounted somewhere reachable, point PROTEOMICS_RELEASE_OLD_PACKAGE_ROOT at it and the
# crosswalk upgrades itself to DIRECT_VERIFICATION with real hashes. No local copy is
# fabricated in the meantime.

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

#' Facts about the original submission established outside this environment.
#'
#' These are the properties the external inspection recorded. They are used to describe the
#' submitted files and to CHECK a locally supplied package if one ever appears -- not as a
#' substitute for having read it here.
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

#' Statement to publish when the original package is not reachable from this environment.
RELEASE_OLD_PACKAGE_EXTERNAL_STATEMENT <- paste(
  "Original submission package directly inspected separately; local project-tree crosswalk",
  "remains based on dimensions/content lineage.")

#' Resolve how the old-package crosswalk was established, for this build.
#'
#' Returns the method, whether the package was reachable here, and the per-file evidence.
#' When the package IS reachable each named file is hashed, so the upgrade to
#' DIRECT_VERIFICATION is backed by bytes this build actually read.
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
# The analysis matrix is standardised separately for each protein across the
# measurement-level dataset, so the coefficient ProTigy stores as `logFC` (and the canonical
# split tables carry as `log2fc`) is NOT a log2 fold change. It is a difference on the
# standardised abundance scale.
#
# Two things this vocabulary deliberately does NOT say:
#
#   * not "standardised mean difference" unqualified -- that term invites a Cohen's d
#     reading, i.e. division by a pooled WITHIN-GROUP SD. The scaling here is a per-protein
#     SD taken across the whole measurement dataset before any grouping, which is a
#     different quantity from a pooled within-group SD.
#   * not "z-scored" -- see release_standardization_evidence(). The released matrix is
#     serialised to 2 decimals, no row is exactly mean 0 / sd 1, and the producing operation
#     is UNRESOLVED in the lineage. "Standardised" is what the numbers support.

RELEASE_EFFECT_SIZE <- list(
  public_field = "effect_size_sd_units",
  public_term = "standardized abundance difference (SD units)",
  public_term_alt = "difference in standardized protein abundance",
  source_field = "logFC",
  source_field_detail = "logFC (ProTigy) / log2fc (canonical split tables)",
  units = "per-protein standard deviation of the standardised abundance scale",

  definition = paste(
    "Difference of group means (numerator minus denominator) on the protein-standardised",
    "abundance scale -- a standardized abundance difference expressed in SD units of that",
    "scale. It is NOT a log2 fold change, despite being stored as `logFC` by ProTigy and",
    "`log2fc` in the canonical split tables. Positive means higher in the numerator",
    "condition. A value of 2 means two standard deviations on the standardised scale, not",
    "a four-fold change."),

  methods_sentence = paste(
    "Protein abundances in the analysis matrix were standardized separately for each",
    "protein across the measurement-level dataset. Consequently, the coefficient stored by",
    "the historical analysis software as `logFC` represents the difference in standardized",
    "protein abundance between groups, expressed on that standardized scale, rather than a",
    "log2 fold change."),

  sensitivity_public_label = "effect-size-ranked sensitivity analysis",
  sensitivity_public_detail = paste(
    "GSEA ranked by the standardised-abundance effect size instead of the moderated t.",
    "Sensitivity analysis only; the canonical ranking statistic is the moderated t."),

  # Filenames and internal column VALUES that legitimately keep the historical token.
  # Renaming these would break provenance against the canonical run, so they stay.
  retained_internal_tokens = c("log2fc", "logFC", "GSEA_log2FC_sensitivity",
                               "GSEA_GO_BP_log2fc_sensitivity",
                               "GSEA_KEGG_log2fc_sensitivity",
                               "ORA_GO_BP_top_abs_log2fc", "top_absolute_log2fc")
)

#' Phrases that would assert a fold-change reading of the effect size.
RELEASE_FOLD_CHANGE_CLAIMS <- c(
  "log2 fold change", "log2 fold-change", "log2-fold change", "log2-fold-change",
  "log 2 fold change", "log2 ratio", "fold change", "fold-change"
)

#' Tokens that make a fold-change mention a citation rather than a claim.
#'
#' A public document MUST be able to say "this is not a log2 fold change" and "the source
#' column is named logFC". What it must not do is describe the published values AS log2
#' fold changes. Telling those apart needs the surrounding words, so a mention is accepted
#' only when its line also carries an explicit negation or provenance marker. `ewce` is
#' exempt because EWCE's own bootstrap statistic genuinely IS a fold change and is
#' published under that name.
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

#' Measure the per-protein standardisation of a matrix, without naming the operation.
#'
#' Reports what the numbers show (row means, both SD conventions) and whether an EXACT
#' z-score can be claimed. For this dataset it cannot: the released values carry 2 decimals,
#' so the residual deviations are serialisation rounding, and the producing step is
#' UNRESOLVED in the lineage. Callers read `exact_zscore` to choose wording rather than
#' asserting "z-scored" and hoping.
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

#' Calls that would mean the publication layer recomputed science instead of reading it.
#'
#' Detection is done on the PARSE TREE, not on the text, so a mention inside a
#' documentation string (the provenance stage quotes the historical imputation formula,
#' `rnorm(n, mean = ...)`, verbatim) or inside a comment is not a hit, while an actual
#' invocation is. Regex over source text cannot make that distinction and produced exactly
#' that false positive.
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
#' Uses the R parser: only tokens the parser classifies as a function call count. If a file
#' cannot be parsed, it falls back to a text scan and marks the hit, because an unparseable
#' builder is itself a problem worth surfacing rather than skipping.
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

#' Defects the previous GSEA_ORA_all_results.xlsx had and this release must not repeat.
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

#' Public column names that must never carry UniProt-style values.
#'
#' The pipeline's internal compatibility column `gene_symbol` holds a UniProt entry name
#' in the split tables and a UniProt accession in the mapped tables. Propagating that
#' name into a published table would mislabel an accession as a gene symbol.
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
