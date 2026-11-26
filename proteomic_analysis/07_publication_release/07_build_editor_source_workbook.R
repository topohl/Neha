#!/usr/bin/env Rscript

# Publication release, stage 07 -- the editor-facing source-data workbook.
#
# Produces
#   editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx
#
# Workbook rules:
#   * one table per sheet, under its title and caption; never tables side by side
#   * no duplicate, blank or placeholder (V1, X1, ...1) column names
#     (release_check_workbook_hygiene(); the build stops on a hit)
#   * every sheet joins to the others on canonical_comparison / protein_group_id /
#     AnimalID + sample_class
#   * a README sheet that states the statistical unit first
#
# Gene-membership columns (GSEA core_enrichment, ORA geneID) stay in the gzipped TSVs:
# single cells reach thousands of characters and would add tens of megabytes. Row counts are
# unchanged, n_core_enrichment_genes is kept, and the README sheet says so.
#
# The workbook is the human-readable copy: display headers, a Contents sheet with assigned
# supplementary-table numbers, captions, frozen header rows, filters, column widths and
# number formats. Columns that are constant on a sheet move into its caption. The .tsv files
# keep the machine names; R/release_presentation.R maps between the two.

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
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_presentation.R"))
release_require("digest", "openxlsx")

DATA_ROOT <- release_data_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/07_build_editor_source_workbook.R"

release_banner("stage 07 -- editor source-data workbook")

inv <- RELEASE_DESIGN_INVARIANTS

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
experimenter_metadata_wb <- read_release(rel("metadata", "experimenter_metadata.tsv"))
mv_wb <- function(id) release_metadata_value(experimenter_metadata_wb, id)
differential <- read_release(rel("differential_analysis", "primary_differential_proteins.tsv.gz"))
diff_summary <- read_release(rel("differential_analysis", "primary_differential_summary.tsv"))
gsea_go <- read_release(rel("enrichment", "primary_GSEA_GO_BP.tsv.gz"))
gsea_cc <- read_release(rel("enrichment", "primary_GSEA_GO_CC.tsv.gz"))
gsea_kegg <- read_release(rel("enrichment", "primary_GSEA_KEGG.tsv.gz"))
ora <- read_release(rel("enrichment", "primary_ORA_GO_BP.tsv.gz"))
ora_cc <- read_release(rel("enrichment", "primary_ORA_GO_CC.tsv.gz"))
sens <- read_release(rel("enrichment", "GSEA_log2FC_sensitivity.tsv.gz"))
untested_sets <- read_release(rel("enrichment", "gsea_untested_gene_sets.tsv"))
corrections_wb <- read_release(rel("metadata", "sample_class_corrections.tsv"))
n_class_corrected_wb <- nrow(corrections_wb)
ewce <- read_release(rel("enrichment", "primary_EWCE.tsv.gz"))
coverage <- read_release(rel("enrichment", "enrichment_coverage.tsv"))
figure_map <- read_release(rel("editor_source_data", "figure_source_map.tsv"))
run_params <- read_release(rel("enrichment", "enrichment_run_parameters.tsv"))
# The permutation note as stage 05 published it, and EWCE's own differential model.
perm_sentence <- run_params$value[run_params$parameter == "gsea_n_perm_simple (documented)"]
if (length(perm_sentence) != 1L) {
  stop("enrichment_run_parameters.tsv lacks the documented nPermSimple row.", call. = FALSE)
}
ewce_model <- release_ewce_differential_model(DATA_ROOT)

DROP_GENE_LISTS <- c("core_enrichment", "leading_edge")
drop_cols <- function(d, cols) d[, setdiff(names(d), cols), drop = FALSE]
# --------------------------------------------------------------------------------------
# software / database versions
# --------------------------------------------------------------------------------------
# Same rows as provenance/software_versions.tsv (release_build_software_versions()). The
# `status` and `applies_to_stage` columns are kept so superseded versions and UNKNOWN values
# read as such.
sw <- release_build_software_versions(DATA_ROOT, REPO_ROOT)
# The model and its limma / R versions, as the statistics run recorded them.
stat_run <- release_read_stat_run(DATA_ROOT)

# Retained for the registry provenance recorded at the end of this stage.
pkg_versions_path <- file.path(release_enrichment_root("GO_BP", DATA_ROOT), "audits",
                               "package_database_versions.csv")
ewce_session_path <- file.path(DATA_ROOT, "03_output", "ewce",
                               "EWCE_Results_animal_level_validation_20260825",
                               "03_QC_Mapping_Logs", "reproducibility_session_info.txt")
release_assert_exists(pkg_versions_path, "canonical package_database_versions.csv")
release_assert_exists(ewce_session_path, "canonical EWCE sessionInfo")

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
    "submission ranked by the effect size stored as `logFC`. Those",
    "effect-size-ranked results are retained here only as a",
    "sensitivity analysis, in a separate file, with analysis_role = sensitivity on every",
    "row."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "6. PRIMARY COMPARISONS", statement = paste(
    "There are exactly 12 primary comparisons: 4 sample classes (mcherry, neuropil, cfos,",
    "neuron) x 3 within-sample-class contrast families (paired_cno vs paired_veh;",
    "paired_veh vs unpaired_veh; unpaired_cno vs unpaired_veh). See the",
    "Primary_Contrasts sheet."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "7. SECONDARY ANALYSES ARE SEPARATE", statement = paste(
    "Cross-compartment comparisons (Supp. Fig. 4 E), the effect-size-ranked GSEA",
    "sensitivity analysis, the post-hoc Pairing x CNO interaction, the leave-one-animal-out",
    "robustness check of the mCherry learning comparison, the EWCE baseline analysis and",
    "the superseded pooled-model statistics (kept for comparison only) are SECONDARY. They",
    "are listed on the Secondary_Analyses sheet and are never mixed into the primary",
    "contrast table."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "8. COLLECTION PLATE", statement = paste(
    "Because pairing condition was completely associated with collection plate in the",
    "paired-vehicle versus unpaired-vehicle comparison, any collection-plate-associated",
    "contribution cannot be distinguished from a pairing-associated contribution.",
    "No downstream proteomics batch effect attributable to collection plate has been",
    "demonstrated. Collection plate is NOT demonstrated proteomics batch metadata: no",
    "preparation, digestion, LC-MS, acquisition or instrument batch information exists",
    "for this dataset."), stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "9. WHAT CHANGED IN THE RESULTS", statement = paste0(
    "Correcting the statistical unit did not leave every finding intact. Of the ",
    inv$n_primary_contrasts, " primary comparisons, ", n_contrasts_with_hits,
    if (n_contrasts_with_hits == 1L) " contains" else " contain",
    " any FDR-significant protein and ", inv$n_primary_contrasts - n_contrasts_with_hits,
    if (inv$n_primary_contrasts - n_contrasts_with_hits == 1L) " contains" else " contain",
    " none. The mCherry paired-vehicle versus unpaired-vehicle comparison retains a broad",
    " signal (", mch$n_significant_fdr_0_05_mapped_only, " of ",
    mch$n_proteins_tested_mapped_only, " mapped proteins at FDR < 0.05), but its origin is",
    " unresolved. See REVISION_PROTEOMICS_DATA_CHANGELOG.md for the panel-by-panel",
    " comparison."),
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
    "as a result rather than a missing file.", perm_sentence),
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
    "<-> ORA_* <-> EWCE); protein_group_id (Differential_Proteins <-> the feature",
    "annotation in the full release); AnimalID + sample_class (Animal_Level_Metadata <->",
    "Sample_Metadata via animal_level_sample_id).", ewce_model$sentence),
    stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "15. VERIFIED AND UNKNOWN METADATA ARE DISTINGUISHED", statement = paste(
    "Source-controlled experimenter records establish", mv_wb("organism_part"),
    "material,", mv_wb("strain_or_breed"), "mice,", mv_wb("labeling_strategy"),
    "and sequential Lys-C plus trypsin digestion. Instrument model, acquisition mode,",
    "search settings, AnimalID-level sex, developmental stage and age at tissue collection",
    "remain unavailable and have not been inferred. See the Metadata_Field_Status sheet,",
    "metadata/experimenter_metadata.tsv, metadata/sample_preparation_protocol.tsv and",
    "pride/SDRF_MISSING_METADATA.md."),
    stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "16. THE EFFECT SIZE IS NOT A LOG2 FOLD CHANGE", statement = paste(
    "The abundance matrix is standardised separately for each protein, so the coefficient",
    "the statistical GCT stores as `logFC` (limma's name for it) is a standardized",
    "abundance difference in SD units of that scale, not a log2 ratio. A value",
    "of 2 means two standard deviations,",
    "not four-fold. It is published as `effect_size_sd_units`, with",
    "`source_statistic_field = logFC` for provenance, and the exported numbers are",
    "bit-identical to the canonical values. Only the label changed. Figure axes and",
    "manuscript text still say log2FC and must be relabelled by hand; see",
    "editor_source_data/MANUSCRIPT_TERMINOLOGY_ACTIONS.md."),
    stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "17. SIX SAMPLE CLASSES WERE CORRECTED, AND THAT IS RESOLVED", statement = paste0(
    "For ", n_class_corrected_wb, " of the 96 acquisitions -- the left hemispheres of ",
    paste(sort(unique(corrections_wb$AnimalID)), collapse = " and "),
    " -- the sample-class assignment was corrected during historical sample-identity",
    " quality control: a cyclic reassignment among the mCherry-, neuropil- and",
    " neuron-enriched samples. Only metadata were reassigned; acquisition identifiers and",
    " quantitative abundance profiles were unchanged, and the affected animal-level units",
    " were built from the corrected assignments in the first place, so no reanalysis is",
    " implied. Sample_Metadata carries BOTH labels (original_sample_class and",
    " analysis_sample_class) with sample_class_corrected marking the six; `sample_class`",
    " is the analysis-time class. Per-row provenance is in",
    " metadata/sample_class_corrections.tsv. This is NOT an unresolved discrepancy and NOT",
    " a mislabelling."),
    stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(item = "18. STATISTICAL MODEL", statement = paste0(
    RELEASE_PRIMARY_MODEL$methods_sentence, " Software: limma ", stat_run$limma_version,
    " in R ", stat_run$r_version, " (", RELEASE_PRIMARY_MODEL$script, "), as recorded in",
    " the statistics run manifest. Each sample class therefore has its own variance",
    " estimate, and average_standardized_abundance (limma AveExpr) is the mean over that",
    " sample class's ", stat_run$n_columns_in_model, " units. It replaces the earlier",
    " pooled model (", RELEASE_PRIMARY_MODEL$superseded_label, ") used by the release of",
    " 2026-09-07; see section 7b of REVISION_PROTEOMICS_DATA_CHANGELOG.md and the",
    " Software_Versions sheet."),
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
  Software_Versions = sw,
  # GO-CC appended last (S17, S18) so the Contents stays in table-number order and S1-S16
  # keep the numbers already cited.
  GSEA_GO_CC = drop_cols(gsea_cc, DROP_GENE_LISTS),
  ORA_GO_CC = drop_cols(ora_cc, DROP_GENE_LISTS),
  GSEA_Untested_Sets = untested_sets
)


# Sheet-name length limit is 31 characters; stop rather than truncate.
too_long <- names(sheets)[nchar(names(sheets)) > 31L]
if (length(too_long)) {
  stop("Sheet name(s) exceed the 31-character Excel limit: ",
       paste(too_long, collapse = ", "), call. = FALSE)
}

# Every sheet needs a caption; a sheet added without one would publish untitled.
missing_meta <- setdiff(names(sheets), names(RELEASE_SHEET_META))
if (length(missing_meta)) {
  stop("Sheet(s) with no caption in RELEASE_SHEET_META: ",
       paste(missing_meta, collapse = ", "), call. = FALSE)
}

# --------------------------------------------------------------------------------------
# presentation: column order, hoisted constants, display headers
# --------------------------------------------------------------------------------------
# Columns a reader filters or joins on stay even when constant on the current data, so the
# sheet keeps working when the data changes. Everything else that never varies is hoisted
# into the caption instead of repeating down 64,188 rows.

KEEP_WHEN_CONSTANT <- c("sample_class", "condition", "canonical_comparison",
                        "canonical_contrast", "contrast_family", "numerator_condition",
                        "denominator_condition", "analysis_unit", "primary_or_secondary",
                        "analysis_role", "ontology", "analysis", "query_list",
                        "panel_label", "figure", "panel", "component", "field",
                        "parameter", "item", "AnimalID", "sample_id", "term_id",
                        "gene_symbol", "uniprot_accession", "protein_group_id")

presented <- list()
for (nm in names(sheets)) {
  d <- release_order_columns(sheets[[nm]], nm)
  const <- release_sheet_constants(d, keep = KEEP_WHEN_CONSTANT)
  if (length(const$constants)) d <- d[, setdiff(names(d), const$constants), drop = FALSE]
  presented[[nm]] <- list(data = d, constants = const)
}

# Hygiene runs on BOTH namings. The machine names carry the joinability and the
# gene_symbol/uniprot_accession check; the display names are what a reader sees, and a
# duplicate or blank among them would be just as much a defect.
machine_cols <- lapply(presented, function(p) names(p$data))
display_cols <- lapply(presented, function(p) release_display_labels(names(p$data)))
hygiene <- c(release_check_workbook_hygiene(machine_cols),
             release_check_workbook_hygiene(display_cols))
if (length(hygiene)) {
  stop("Workbook hygiene violations:\n", paste("  -", hygiene, collapse = "\n"), call. = FALSE)
}

# Every published column needs a display label in RELEASE_COLUMN_LABELS.
unlabelled <- release_unlabelled_columns(unlist(machine_cols, use.names = FALSE))
if (length(unlabelled)) {
  stop("Column(s) with no display label in RELEASE_COLUMN_LABELS: ",
       paste(unlabelled, collapse = ", "), call. = FALSE)
}

# No published gene_symbol column may contain UniProt identifiers.
for (nm in names(presented)) {
  d <- presented[[nm]]$data
  if ("gene_symbol" %in% names(d) &&
      release_column_is_misleading_gene_symbol(d$gene_symbol)) {
    stop("Sheet ", nm, " has a gene_symbol column holding UniProt identifiers.",
         call. = FALSE)
  }
}

EXCEL_MAX_ROWS <- 1048575L
EXCEL_MAX_CELL <- 32767L
for (nm in names(presented)) {
  d <- presented[[nm]]$data
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
}

# --------------------------------------------------------------------------------------
# build the workbook
# --------------------------------------------------------------------------------------

INK <- "#1F3864"      # header fill
INK_LIGHT <- "#EDF1F8"  # banded / caption fill
RULE <- "#B4C6E7"

style_header <- openxlsx::createStyle(
  textDecoration = "bold", fgFill = INK, fontColour = "#FFFFFF",
  halign = "left", valign = "bottom", wrapText = TRUE,
  border = "bottom", borderColour = INK, borderStyle = "medium")
style_title <- openxlsx::createStyle(textDecoration = "bold", fontSize = 13,
                                     fontColour = INK)
style_caption <- openxlsx::createStyle(fontSize = 10, fontColour = "#333333",
                                       wrapText = TRUE, valign = "top")
style_note <- openxlsx::createStyle(fontSize = 10, fontColour = "#555555",
                                    textDecoration = "italic",
                                    wrapText = TRUE, valign = "top")
style_text <- openxlsx::createStyle(valign = "top")
style_wrap <- openxlsx::createStyle(wrapText = TRUE, valign = "top")

WRAP_ABOVE <- 60L   # free-text columns wider than this get wrapText

wb <- openxlsx::createWorkbook()
openxlsx::modifyBaseFont(wb, fontSize = 10, fontName = "Calibri")

# ---- Contents sheet -------------------------------------------------------------------
# Assigned table numbers, not positional ones, so inserting a sheet cannot renumber a
# table an editor has already cited.
contents <- do.call(rbind, lapply(names(presented), function(nm) {
  meta <- RELEASE_SHEET_META[[nm]]
  cst <- presented[[nm]]$constants
  data.frame(
    Table = ifelse(is.na(meta$number), "--", meta$number),
    Sheet = nm,
    Title = meta$title,
    Rows = nrow(presented[[nm]]$data),
    Columns = ncol(presented[[nm]]$data),
    Description = meta$caption,
    `Applies to the whole table` = if (length(cst$constants)) {
      paste(paste0(release_display_labels(cst$constants), ": ", cst$values),
            collapse = "; ")
    } else "",
    check.names = FALSE, stringsAsFactors = FALSE)
}))

openxlsx::addWorksheet(wb, "Contents", gridLines = FALSE, tabColour = INK)
openxlsx::writeData(wb, "Contents", "Proteomics source data", startRow = 1, startCol = 1)
openxlsx::addStyle(wb, "Contents", style_title, rows = 1, cols = 1)
openxlsx::writeData(
  wb, "Contents",
  paste0("Associative memory proteomics. Animal-level analysis: ",
         inv$n_animals, " animals, ", inv$n_animal_level_units,
         " animal x sample-class units, n = ", inv$n_animals_per_stratum,
         " per condition per sample class, ", inv$n_primary_contrasts,
         " primary contrasts, one limma model per sample class (residual df ",
         stat_run$df_residual, "). The animal is the independent experimental unit."),
  startRow = 2, startCol = 1)
openxlsx::addStyle(wb, "Contents", style_caption, rows = 2, cols = 1)
openxlsx::writeData(
  wb, "Contents",
  paste0("Read the README sheet before any table. Effect sizes are standardized abundance ",
         "differences in SD units, not log2 fold changes. Generated ",
         release_timestamp_utc(), "."),
  startRow = 3, startCol = 1)
openxlsx::addStyle(wb, "Contents", style_note, rows = 3, cols = 1)
openxlsx::writeData(wb, "Contents", contents, startRow = 5, headerStyle = style_header)
openxlsx::addStyle(wb, "Contents", style_wrap,
                   rows = 6:(5 + nrow(contents)), cols = c(3, 6, 7), gridExpand = TRUE)
openxlsx::setColWidths(wb, "Contents",
                       cols = 1:7, widths = c(7, 26, 44, 9, 9, 62, 46))
openxlsx::freezePane(wb, "Contents", firstActiveRow = 6)

# ---- data sheets ----------------------------------------------------------------------
# Layout per sheet: title (row 1), caption (row 2), table-level constants (row 3 when
# present), header row, data. The header row is frozen, so the header stays visible while
# scrolling 64,188 rows.
for (nm in names(presented)) {
  d <- presented[[nm]]$data
  cst <- presented[[nm]]$constants
  meta <- RELEASE_SHEET_META[[nm]]
  headers <- release_display_labels(names(d))

  openxlsx::addWorksheet(wb, nm, gridLines = FALSE)

  title_text <- if (is.na(meta$number)) meta$title else
    paste0("Supplementary Table ", meta$number, ". ", meta$title)
  openxlsx::writeData(wb, nm, title_text, startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, nm, style_title, rows = 1, cols = 1)
  openxlsx::writeData(wb, nm, meta$caption, startRow = 2, startCol = 1)
  openxlsx::addStyle(wb, nm, style_caption, rows = 2, cols = 1)

  header_row <- 3L
  if (length(cst$constants)) {
    openxlsx::writeData(
      wb, nm,
      paste0("Applies to every row: ",
             paste(paste0(release_display_labels(cst$constants), " = ", cst$values),
                   collapse = "; "), "."),
      startRow = 3, startCol = 1)
    openxlsx::addStyle(wb, nm, style_note, rows = 3, cols = 1)
    header_row <- 4L
  }

  names(d) <- headers
  openxlsx::writeData(wb, nm, d, startRow = header_row, headerStyle = style_header,
                      withFilter = nrow(d) > 1L)
  openxlsx::setRowHeights(wb, nm, rows = header_row, heights = 30)
  openxlsx::freezePane(wb, nm, firstActiveRow = header_row + 1L)

  widths <- integer(ncol(d))
  for (j in seq_along(d)) {
    src <- presented[[nm]]$data[[j]]
    widths[[j]] <- release_column_width(src, headers[[j]])
    fmt <- release_numfmt_for(src, names(presented[[nm]]$data)[[j]])
    if (!is.na(fmt)) {
      openxlsx::addStyle(wb, nm, openxlsx::createStyle(numFmt = fmt, valign = "top"),
                         rows = (header_row + 1L):(header_row + nrow(d)), cols = j,
                         gridExpand = TRUE)
    } else if (widths[[j]] >= WRAP_ABOVE) {
      openxlsx::addStyle(wb, nm, style_wrap,
                         rows = (header_row + 1L):(header_row + nrow(d)), cols = j,
                         gridExpand = TRUE)
    } else {
      openxlsx::addStyle(wb, nm, style_text,
                         rows = (header_row + 1L):(header_row + nrow(d)), cols = j,
                         gridExpand = TRUE)
    }
  }
  openxlsx::setColWidths(wb, nm, cols = seq_along(d), widths = widths)

  release_log("    ", format(nm, width = 26),
              format(nrow(d), width = 7), " rows x ", format(ncol(d), width = 2), " cols",
              if (length(cst$constants))
                paste0("  (", length(cst$constants), " constant col(s) hoisted)") else "")
}

xlsx_path <- release_path("editor_source_data", "Proteomics_Source_Data_Animal_Level.xlsx")
openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
size_mb <- round(file.info(xlsx_path)$size / 1024 / 1024, 1)
release_log("  wrote Proteomics_Source_Data_Animal_Level.xlsx (", length(presented) + 1L,
            " sheets incl. Contents, ", size_mb, " MB)")

release_register("editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx",
                 "editor-facing source-data workbook, one table per sheet",
                 c(rel("metadata", "sample_metadata.tsv"),
                   rel("differential_analysis", "primary_differential_proteins.tsv.gz"),
                   rel("enrichment", "primary_GSEA_GO_BP.tsv.gz"),
                   rel("enrichment", "primary_GSEA_GO_CC.tsv.gz"),
                   rel("enrichment", "primary_ORA_GO_BP.tsv.gz"),
                   rel("enrichment", "primary_ORA_GO_CC.tsv.gz"),
                   rel("enrichment", "primary_EWCE.tsv.gz"),
                   pkg_versions_path, ewce_session_path, stat_run$path, ewce_model$path),
                 c(rep(NA_character_, 7), release_sha256(pkg_versions_path),
                   release_sha256(ewce_session_path), release_sha256(stat_run$path),
                   release_sha256(ewce_model$path)),
                 STAGE, "xlsx")

# The workbook sheet inventory, so a validator can check it without opening Excel. Both
# namings are recorded: `columns` is the machine naming the checks key off, and
# `display_columns` is what the reader sees.
inventory <- data.frame(
  sheet = names(presented),
  supplementary_table = vapply(names(presented), function(nm) {
    n <- RELEASE_SHEET_META[[nm]]$number
    ifelse(is.na(n), "", n)
  }, character(1)),
  title = vapply(names(presented), function(nm) RELEASE_SHEET_META[[nm]]$title,
                 character(1)),
  n_rows = vapply(presented, function(p) nrow(p$data), integer(1)),
  n_columns = vapply(presented, function(p) ncol(p$data), integer(1)),
  columns = vapply(presented, function(p) paste(names(p$data), collapse = ";"),
                   character(1)),
  display_columns = vapply(presented, function(p) {
    paste(release_display_labels(names(p$data)), collapse = ";")
  }, character(1)),
  table_level_constants = vapply(presented, function(p) {
    if (!length(p$constants$constants)) return("")
    paste(paste0(p$constants$constants, "=", p$constants$values), collapse = ";")
  }, character(1)),
  stringsAsFactors = FALSE, check.names = FALSE
)
rownames(inventory) <- NULL
w <- release_write_table(inventory,
                         release_path("editor_source_data", "workbook_sheet_inventory.tsv"))
release_register("editor_source_data/workbook_sheet_inventory.tsv",
                 "sheet-by-sheet inventory of the editor workbook: machine and display columns",
                 "editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx",
                 NA_character_, STAGE, "tsv")
release_log("  wrote workbook_sheet_inventory.tsv (", w$rows, "x", w$cols, ")")
release_log("stage 07 complete")
