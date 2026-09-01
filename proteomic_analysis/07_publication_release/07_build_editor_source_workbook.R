#!/usr/bin/env Rscript

# Publication release, stage 07 -- the editor-facing source-data workbook.
#
# Produces
#   editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx
#
# Deliberately NOT like the previous GSEA_ORA_all_results.xlsx. Rules enforced in code:
#   * one table per sheet, starting at A1 -- never independent tables pasted side by side
#   * no duplicate column names on any sheet
#   * no unnamed or placeholder (V1, X1, ...1) columns
#   * every sheet joinable to the others on canonical_comparison / protein_group_id /
#     AnimalID + sample_class
#   * a README sheet that states the statistical unit before anything else
# release_check_workbook_hygiene() asserts the column rules and the build stops on a hit.
#
# Gene-membership columns (GSEA core_enrichment, ORA geneID) are omitted from the workbook
# and kept in the gzipped TSVs. Individual cells reach several thousand characters and the
# column alone would add tens of megabytes to a file an editor has to open. The row counts
# are unchanged, n_core_enrichment_genes is retained, and the README says so explicitly --
# nothing is silently dropped.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_validation.R"))
release_require("digest", "writexl")

DATA_ROOT <- release_data_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/07_build_editor_source_workbook.R"

release_banner("stage 07 -- editor source-data workbook")

rel <- function(...) release_path(..., create_dir = FALSE)
read_release <- function(path) {
  release_assert_exists(path, basename(path))
  if (grepl("[.]gz$", path)) {
    read.delim(gzfile(path), sep = "\t", quote = "", stringsAsFactors = FALSE,
               check.names = FALSE)
  } else {
    read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  }
}

sample_metadata <- read_release(rel("metadata", "sample_metadata.tsv"))
animal_metadata <- read_release(rel("metadata", "animal_level_sample_metadata.tsv"))
primary_contrasts <- read_release(rel("metadata", "primary_contrast_manifest.tsv"))
secondary <- read_release(rel("metadata", "secondary_analysis_manifest.tsv"))
field_provenance <- read_release(rel("metadata", "metadata_field_provenance.tsv"))
differential <- read_release(rel("differential_analysis", "primary_differential_proteins.tsv.gz"))
diff_summary <- read_release(rel("differential_analysis", "primary_differential_summary.tsv"))
gsea_go <- read_release(rel("enrichment", "primary_GSEA_GO_BP.tsv.gz"))
gsea_kegg <- read_release(rel("enrichment", "primary_GSEA_KEGG.tsv.gz"))
ora <- read_release(rel("enrichment", "primary_ORA_GO_BP.tsv.gz"))
sens <- read_release(rel("enrichment", "GSEA_log2FC_sensitivity.tsv.gz"))
ewce <- read_release(rel("enrichment", "primary_EWCE.tsv.gz"))
coverage <- read_release(rel("enrichment", "enrichment_coverage.tsv"))
figure_map <- read_release(rel("editor_source_data", "figure_source_map.tsv"))
run_params <- read_release(rel("enrichment", "enrichment_run_parameters.tsv"))

DROP_GENE_LISTS <- c("core_enrichment", "leading_edge")
drop_cols <- function(d, cols) d[, setdiff(names(d), cols), drop = FALSE]

# --------------------------------------------------------------------------------------
# software / database versions as recorded by the canonical runs
# --------------------------------------------------------------------------------------

pkg_versions_path <- file.path(DATA_ROOT, "03_output", "enrichment",
                               "enrichment_t_rank_validation_20260825", "audits",
                               "package_database_versions.csv")
release_assert_exists(pkg_versions_path, "canonical package_database_versions.csv")
pkg_versions <- release_read_csv(pkg_versions_path)

ewce_session_path <- file.path(DATA_ROOT, "03_output", "ewce",
                               "EWCE_Results_animal_level_validation_20260825",
                               "03_QC_Mapping_Logs", "reproducibility_session_info.txt")
release_assert_exists(ewce_session_path, "canonical EWCE sessionInfo")
ewce_session <- readLines(ewce_session_path, warn = FALSE)
parse_attached <- function(lines) {
  start <- grep("^other attached packages", lines)
  if (!length(start)) return(character(0))
  tail_lines <- lines[(start[[1]] + 1L):length(lines)]
  stop_at <- grep("^loaded via a namespace", tail_lines)
  if (length(stop_at)) tail_lines <- tail_lines[seq_len(stop_at[[1]] - 1L)]
  toks <- unlist(strsplit(paste(tail_lines, collapse = " "), "\\s+"))
  toks <- toks[grepl("^[A-Za-z][A-Za-z0-9._]*_[0-9]", toks)]
  unique(toks)
}
ewce_pkgs <- parse_attached(ewce_session)

protigy_params_path <- file.path(DATA_ROOT, "01_input", "raw_proteomics",
                                 "20251107_pg.matrix_Neha", "params.txt")
protigy_version <- NA_character_
if (file.exists(protigy_params_path)) {
  pl <- readLines(protigy_params_path, warn = FALSE)
  hit <- grep("Protigy", pl, value = TRUE, ignore.case = TRUE)
  if (length(hit)) {
    protigy_version <- trimws(gsub("^##\\s*|\\s*$", "", hit[[1]]))
  }
}

sw <- rbind(
  data.frame(component = "R", version = pkg_versions$version[pkg_versions$component == "R"],
             recorded_by = "enrichment run audit (2026-08-25)",
             evidence_path = pkg_versions_path, stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(component = pkg_versions$component[pkg_versions$component != "R"],
             version = pkg_versions$version[pkg_versions$component != "R"],
             recorded_by = "enrichment run audit (2026-08-25)",
             evidence_path = pkg_versions_path, stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(component = sub("_.*$", "", ewce_pkgs), version = sub("^[^_]*_", "", ewce_pkgs),
             recorded_by = "EWCE run sessionInfo (2026-08-25)",
             evidence_path = ewce_session_path, stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(component = "ProTigy",
             version = ifelse(is.na(protigy_version), "UNKNOWN", protigy_version),
             recorded_by = ifelse(is.na(protigy_version), "NOT RECOVERED",
                                  "ProTigy params.txt header (2025-11-07 run)"),
             evidence_path = protigy_params_path, stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(component = "DIA-NN (search / quantification)", version = "UNKNOWN",
             recorded_by = paste("software identity evidenced by the DIA-NN pg_matrix column",
                                 "schema and a report.stats-derived annotation sheet;",
                                 "no version string exists in the project"),
             evidence_path = "NONE", stringsAsFactors = FALSE, check.names = FALSE)
)
sw <- sw[!duplicated(paste(sw$component, sw$version)), , drop = FALSE]
sw <- sw[order(sw$component), , drop = FALSE]
rownames(sw) <- NULL

# --------------------------------------------------------------------------------------
# README sheet
# --------------------------------------------------------------------------------------

mch <- diff_summary[diff_summary$canonical_comparison ==
                      "mcherry_paired_veh_over_mcherry_unpaired_veh", , drop = FALSE]
n_contrasts_with_hits <- sum(as.integer(diff_summary$n_significant_fdr_0_05) > 0L)

readme <- rbind(
  data.frame(item = "1. STATISTICAL UNIT", statement = paste(
    "The ANIMAL is the independent experimental unit for every inferential analysis in",
    "this package."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "2. HEMISPHERE AVERAGING", statement = paste(
    "Left and right hemisphere measurements were averaged within AnimalID x sample_class",
    "BEFORE any statistical test. 96 acquisition measurements therefore become 48",
    "animal-level units (12 animals x 4 sample classes)."),
    stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "3. GROUP SIZE", statement = paste(
    "n = 3 animals per condition per sample class, in all 16 strata, with complete",
    "left/right pairing."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "4. THE ORIGINAL SUBMISSION IS SUPERSEDED", statement = paste(
    "The originally submitted analysis treated the 96 hemisphere/acquisition measurements",
    "as independent replicates (n = 6 per group), which is pseudoreplication: two",
    "hemispheres of one animal are not independent samples. Every inferential number in",
    "the original submission, and the GSEA_ORA_all_results.xlsx workbook that carried",
    "them, is superseded by this package."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "5. GSEA RANKING STATISTIC", statement = paste(
    "The canonical GSEA ranks proteins by the MODERATED t STATISTIC. The original",
    "submission ranked by log2FC. log2FC-ranked results are retained here only as a",
    "sensitivity analysis, in a separate file, with analysis_role = sensitivity on every",
    "row."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "6. PRIMARY COMPARISONS", statement = paste(
    "There are exactly 12 primary comparisons: 4 sample classes (mcherry, neuropil, cfos,",
    "neuron) x 3 within-sample-class contrast families (paired_cno vs paired_veh;",
    "paired_veh vs unpaired_veh; unpaired_cno vs unpaired_veh). See the",
    "Primary_Contrasts sheet."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "7. SECONDARY ANALYSES ARE SEPARATE", statement = paste(
    "Cross-compartment comparisons (Supplementary E), the log2FC-ranked GSEA sensitivity",
    "analysis, the post-hoc Pairing x CNO interaction and the EWCE baseline analysis are",
    "SECONDARY. They are listed on the Secondary_Analyses sheet and are never mixed into",
    "the primary contrast table."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "8. COLLECTION PLATE", statement = paste(
    "Because pairing condition was completely associated with collection plate in the",
    "paired-vehicle versus unpaired-vehicle comparison, any collection-plate-associated",
    "contribution cannot be distinguished from a pairing-associated contribution.",
    "No downstream proteomics batch effect attributable to collection plate has been",
    "demonstrated. Collection plate is NOT demonstrated proteomics batch metadata: no",
    "preparation, digestion, LC-MS, acquisition or instrument batch information exists",
    "for this dataset."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "9. WHAT CHANGED IN THE RESULTS", statement = paste0(
    "Correcting the statistical unit did not leave every finding intact. Of the 12 primary",
    " comparisons, ", n_contrasts_with_hits, " contain any FDR-significant protein and ",
    12L - n_contrasts_with_hits, " contain none. The mCherry paired-vehicle versus",
    " unpaired-vehicle comparison retains a broad signal (",
    mch$n_significant_fdr_0_05_mapped_only, " of ", mch$n_proteins_tested_mapped_only,
    " mapped proteins at FDR < 0.05), but its origin is unresolved. See",
    " REVISION_PROTEOMICS_DATA_CHANGELOG.md for the panel-by-panel comparison."),
    stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "10. FIELD NAMES", statement = paste(
    "protein_group_id holds UniProt ENTRY NAMES (e.g. A0A0J9YTR2_MOUSE), semicolon-",
    "separated for multi-protein groups. uniprot_accession holds the accession of the",
    "leading protein. gene_symbol holds a resolved gene symbol and never a UniProt",
    "identifier. protein_description holds a real protein description from the search",
    "output. An internal pipeline column named gene_symbol carries UniProt identifiers;",
    "that name is NOT used in this package."),
    stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "11. ENRICHMENT TABLES ARE COMPLETE, NOT PRE-FILTERED", statement = paste(
    "The canonical enrichment run used pvalue_cutoff = 1 and qvalue_cutoff = 1, so the",
    "GSEA and ORA sheets contain every TESTED term, not only significant ones. Filter on",
    "adjusted_p_value. Where an ORA query list was empty (no FDR-significant proteins in",
    "that comparison) no rows exist; enrichment_coverage in the full release records that",
    "explicitly as a result rather than a missing file."),
    stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "12. GENE MEMBERSHIP LISTS", statement = paste(
    "The GSEA core_enrichment / leading_edge and ORA geneID columns are present in the",
    "gzipped TSV release files but omitted from this workbook: single cells reach several",
    "thousand characters. Row counts here are identical to the TSVs and",
    "n_core_enrichment_genes is retained."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "13. WHAT IS NOT HERE", statement = paste(
    "Native raw mass-spectrometry acquisition files. The 96 `.d` acquisition directories",
    "are named per row in Sample_Metadata (raw_file_basename) but are not part of this",
    "package; they are intended for a PRIDE deposition. No file in this package is raw",
    "mass-spectrometry data."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "14. SHEET RELATIONSHIPS", statement = paste(
    "Join on: canonical_comparison (Primary_Contrasts <-> Differential_Proteins <-> GSEA_*",
    "<-> ORA_GO_BP <-> EWCE); protein_group_id (Differential_Proteins <-> the feature",
    "annotation in the full release); AnimalID + sample_class (Animal_Level_Metadata <->",
    "Sample_Metadata via animal_level_sample_id)."),
    stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "15. UNKNOWN METADATA IS MARKED UNKNOWN", statement = paste(
    "Instrument model, digestion enzyme, DIA acquisition method, search-software version,",
    "labelling, fractionation, animal sex and animal age are NOT recorded anywhere in the",
    "project and have deliberately been left blank rather than inferred. See the",
    "Metadata_Field_Status sheet and pride/SDRF_MISSING_METADATA.md."),
    stringsAsFactors = FALSE, check.names = FALSE)
)

# --------------------------------------------------------------------------------------
# assemble sheets
# --------------------------------------------------------------------------------------

sheets <- list(
  README = readme,
  Sample_Metadata = sample_metadata,
  Animal_Level_Metadata = animal_metadata,
  Primary_Contrasts = primary_contrasts,
  Differential_Summary = diff_summary,
  Differential_Proteins = drop_cols(differential, "source_matrix_gene_symbols"),
  GSEA_GO_BP = drop_cols(gsea_go, DROP_GENE_LISTS),
  GSEA_KEGG = drop_cols(gsea_kegg, DROP_GENE_LISTS),
  ORA_GO_BP = drop_cols(ora, DROP_GENE_LISTS),
  GSEA_log2FC_sensitivity = drop_cols(sens, DROP_GENE_LISTS),
  EWCE = ewce,
  Enrichment_Coverage = coverage,
  Enrichment_Parameters = run_params,
  Secondary_Analyses = secondary,
  Figure_Source_Map = figure_map,
  Metadata_Field_Status = field_provenance,
  Software_Versions = sw
)

# Sheet-name length limit is 31 characters; check rather than truncate silently.
too_long <- names(sheets)[nchar(names(sheets)) > 31L]
if (length(too_long)) {
  stop("Sheet name(s) exceed the 31-character Excel limit: ",
       paste(too_long, collapse = ", "), call. = FALSE)
}

hygiene <- release_check_workbook_hygiene(lapply(sheets, names))
if (length(hygiene)) {
  stop("Workbook hygiene violations:\n", paste("  -", hygiene, collapse = "\n"), call. = FALSE)
}

# No published gene_symbol column may contain UniProt identifiers.
for (nm in names(sheets)) {
  d <- sheets[[nm]]
  if ("gene_symbol" %in% names(d) &&
      release_column_is_misleading_gene_symbol(d$gene_symbol)) {
    stop("Sheet ", nm, " has a gene_symbol column holding UniProt identifiers.",
         call. = FALSE)
  }
}

EXCEL_MAX_ROWS <- 1048575L
EXCEL_MAX_CELL <- 32767L
for (nm in names(sheets)) {
  d <- sheets[[nm]]
  if (nrow(d) > EXCEL_MAX_ROWS) {
    stop("Sheet ", nm, " has ", nrow(d), " rows, beyond the Excel limit.", call. = FALSE)
  }
  for (col in names(d)) {
    if (is.character(d[[col]])) {
      longest <- suppressWarnings(max(nchar(d[[col]]), na.rm = TRUE))
      if (is.finite(longest) && longest > EXCEL_MAX_CELL) {
        stop("Sheet ", nm, " column ", col, " has a ", longest,
             "-character cell, beyond the Excel limit.", call. = FALSE)
      }
    }
  }
  release_log("    ", format(nm, width = 26), format(nrow(d), width = 7), " rows x ",
              ncol(d), " cols")
}

xlsx_path <- release_path("editor_source_data", "Proteomics_Source_Data_Animal_Level.xlsx")
writexl::write_xlsx(sheets, path = xlsx_path, col_names = TRUE, format_headers = TRUE)
size_mb <- round(file.info(xlsx_path)$size / 1024 / 1024, 1)
release_log("  wrote Proteomics_Source_Data_Animal_Level.xlsx (", length(sheets),
            " sheets, ", size_mb, " MB)")

release_register("editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx",
                 "editor-facing source-data workbook, one table per sheet",
                 c(rel("metadata", "sample_metadata.tsv"),
                   rel("differential_analysis", "primary_differential_proteins.tsv.gz"),
                   rel("enrichment", "primary_GSEA_GO_BP.tsv.gz"),
                   rel("enrichment", "primary_ORA_GO_BP.tsv.gz"),
                   rel("enrichment", "primary_EWCE.tsv.gz"),
                   pkg_versions_path, ewce_session_path),
                 c(rep(NA_character_, 5), release_sha256(pkg_versions_path),
                   release_sha256(ewce_session_path)),
                 STAGE, "xlsx")

# The workbook sheet inventory, so a validator can check it without opening Excel.
inventory <- data.frame(
  sheet = names(sheets),
  n_rows = vapply(sheets, nrow, integer(1)),
  n_columns = vapply(sheets, ncol, integer(1)),
  columns = vapply(sheets, function(d) paste(names(d), collapse = ";"), character(1)),
  stringsAsFactors = FALSE, check.names = FALSE
)
rownames(inventory) <- NULL
w <- release_write_table(inventory,
                         release_path("editor_source_data", "workbook_sheet_inventory.tsv"))
release_register("editor_source_data/workbook_sheet_inventory.tsv",
                 "sheet-by-sheet row/column inventory of the editor workbook",
                 "editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx",
                 NA_character_, STAGE, "tsv")
release_log("  wrote workbook_sheet_inventory.tsv (", w$rows, "x", w$cols, ")")
release_log("stage 07 complete")
