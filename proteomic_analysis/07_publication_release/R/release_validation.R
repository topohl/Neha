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
