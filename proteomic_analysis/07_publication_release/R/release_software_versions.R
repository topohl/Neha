# Single source of truth for the software / database version table.
#
# Two release artefacts report software versions to a reader:
#   provenance/software_versions.tsv                     (stage 10)
#   the Software_Versions sheet of the editor workbook    (stage 07)
#
# They used to derive it independently, and they disagreed. Stage 07 read only the 2025
# hemisphere-level ProTigy params.txt header and published "Protigy (v1.1.8)" with no
# qualification, while stage 10 correctly reported 2.4.1 for the canonical animal-level run
# and confined v1.1.8 to the superseded hemisphere-level runs. An editor reading the
# workbook would have concluded the analysis used ProTigy 1.1.8.
#
# The fix is structural rather than a corrected duplicate: both stages now call
# release_build_software_versions(), so the two artefacts cannot drift apart again.
#
# Requires release_utils.R to be sourced first.

# --------------------------------------------------------------------------------------
# ProTigy version for the CANONICAL animal-level run (2026-08-24)
# --------------------------------------------------------------------------------------
# Recovered 2026-09-02 by targeted audit. Recorded as a constant rather than read from
# whatever ProTigy happens to be installed on the machine running this build: the question
# is which version produced a specific run in the past, and a build-time lookup would
# silently answer a different question on a different machine. That is the same trap as
# carrying the 2025 version forward.
#
# The evidence is five strands that agree:
#   1. The only ProTigy present on the analysis machine is 2.4.1, and its installed
#      DESCRIPTION records Packaged 2026-08-24 12:47:48 UTC / Built 2026-08-24 12:48:07 UTC.
#   2. The canonical run's parameter export, neha_proteome_parameters.yaml, was written
#      2026-08-24 13:06:11 UTC -- 18 minutes after that build finished.
#   3. ProTigy v2's exporter (R/tab_export.R) writes paste0(ome, "_parameters.yaml"); the
#      canonical file is neha_proteome_parameters.yaml, i.e. the historical
#      project-specific proteome identifier.
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

#' Fail closed if the artefact the ProTigy 2.4.1 claim rests on stops looking like a v2 export.
#'
#' Without this, a future change to the canonical YAML could leave a version claim standing
#' on evidence that no longer exists.
release_assert_protigy_v2_export_signature <- function(animal_param_yaml) {
  if (!file.exists(animal_param_yaml)) return(invisible(FALSE))
  pp <- readLines(animal_param_yaml, warn = FALSE)
  has_v2_keys <- all(vapply(c("gct_file_name:", "annotation_column:", "group_normalization:",
                              "convert_ids_to_gene_symbol:", "id_mapping_species:"),
                            function(k) any(startsWith(pp, k)), logical(1)))
  # The v1.1.x marker is a literal "## Protigy (vX.Y.Z)" comment line; matched with
  # startsWith on trimmed lines so this check carries no regex escaping of its own.
  has_v1_header <- any(startsWith(trimws(pp), "## Protigy (v"))
  has_gct_file_path <- any(startsWith(pp, "gct_file_path:"))
  if (!has_v2_keys || has_v1_header || has_gct_file_path) {
    stop("The animal-level ProTigy parameter export no longer carries the ProTigy v2 ",
         "export signature that the recovered version ", PROTIGY_ANIMAL_LEVEL_VERSION,
         " rests on (", animal_param_yaml, "). Re-verify the version before releasing.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# --------------------------------------------------------------------------------------
# the table
# --------------------------------------------------------------------------------------

#' Parse the "other attached packages" block of a captured sessionInfo.
release_parse_attached_packages <- function(path) {
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

#' Build the authoritative software / database version table.
#'
#' Every row records where the version came from, so a reader can tell a version captured
#' by the run itself from one recovered later, and an UNKNOWN from an omission. Rows are
#' scoped by `applies_to_stage`, which is what keeps the superseded 2025 hemisphere-level
#' ProTigy version from reading as the version behind the canonical result.
release_build_software_versions <- function(data_root = release_data_root(),
                                            repo_root = release_repo_root()) {
  UNKNOWN <- "UNKNOWN"
  animal_root <- file.path(data_root, "02_data", "animal_level")
  enrich_root <- file.path(data_root, "03_output", "enrichment",
                           "enrichment_t_rank_validation_20260825")
  ewce_root <- file.path(data_root, "03_output", "ewce",
                         "EWCE_Results_animal_level_validation_20260825")
  pca_root <- file.path(data_root, "03_output", "pca",
                        "pca_plots_animal_level_validation_20260825_rerun")

  protigy_params_2025 <- file.path(data_root, "01_input", "raw_proteomics",
                                   "20251107_pg.matrix_Neha", "params.txt")
  animal_param_yaml <- file.path(animal_root, "neha_proteome_parameters.yaml")
  ewce_session_path <- file.path(ewce_root, "03_QC_Mapping_Logs",
                                 "reproducibility_session_info.txt")
  pca_session_path <- file.path(pca_root, "tables", "meta", "sessionInfo.txt")
  p_idmapping <- file.path(data_root, "01_input", "references", "MOUSE_10090_idmapping.dat")
  p_manual_mapping <- file.path(data_root, "01_input", "references", "manual_mapping.xlsx")

  release_assert_protigy_v2_export_signature(animal_param_yaml)

  mapped_index <- release_read_csv(file.path(animal_root, "mapped",
                                             "indexMappedComparisons.csv"))
  pkg_versions <- release_read_csv(file.path(enrich_root, "audits",
                                             "package_database_versions.csv"))

  hash_or <- function(path, fallback = NA_character_) {
    if (!is.na(path) && nzchar(path) && file.exists(path) && !dir.exists(path)) {
      release_sha256(path)
    } else fallback
  }
  sv <- function(component, category, version, status, recorded_by, evidence_path,
                 applies_to, notes = NA_character_) {
    data.frame(component = component, category = category, version = version,
               status = status, recorded_by = recorded_by, evidence_path = evidence_path,
               evidence_sha256 = hash_or(evidence_path), applies_to_stage = applies_to,
               notes = notes, stringsAsFactors = FALSE, check.names = FALSE)
  }

  enrich_versions <- do.call(rbind, lapply(seq_len(nrow(pkg_versions)), function(i) {
    sv(pkg_versions$component[i],
       ifelse(pkg_versions$component[i] == "R", "language", "R package"),
       pkg_versions$version[i], "KNOWN_VERIFIED",
       "enrichment run audit, 2026-08-25",
       file.path(enrich_root, "audits", "package_database_versions.csv"),
       "differential enrichment (GSEA / ORA)")
  }))

  ewce_toks <- release_parse_attached_packages(ewce_session_path)
  ewce_versions <- if (length(ewce_toks)) {
    do.call(rbind, lapply(ewce_toks, function(tok) {
      sv(sub("_.*$", "", tok), "R package", sub("^[^_]*_", "", tok), "KNOWN_VERIFIED",
         "EWCE run sessionInfo, 2026-08-25", ewce_session_path,
         "EWCE cell-type enrichment")
    }))
  } else NULL

  pca_toks <- release_parse_attached_packages(pca_session_path)
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
       as.character(mapped_index$mapping_reference_version[[1]]),
       "KNOWN_BUT_NEEDS_STANDARDIZATION",
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

  out <- rbind(enrich_versions, ewce_versions, pca_versions, external)
  out <- out[!duplicated(paste(out$component, out$version, out$applies_to_stage)), ,
             drop = FALSE]
  out <- out[order(out$applies_to_stage, out$component), , drop = FALSE]
  rownames(out) <- NULL
  out
}
