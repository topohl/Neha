# Presentation layer for the reader-facing workbook.
#
# The release carries every table twice, on purpose, and the two serve different readers:
#
#   *.tsv / *.tsv.gz   machine-readable. Snake_case machine column names, one table per
#                      file, contract-bound: the tests and the data dictionary key off
#                      these names, so they never change for cosmetic reasons.
#   the .xlsx workbook human-readable. This is what an editor or reviewer actually opens,
#                      so it gets display headers, units in the caption rather than
#                      repeated down a column, a table of contents, and formatting.
#
# This file holds only the mapping between the two. It contains no data and no statistics.
#
# Three conventions worth stating, because they are what make the workbook read like a
# supplementary table set rather than a database dump:
#
#   1. A column whose value is IDENTICAL on every row of a sheet is not a column. It is a
#      property of the table, and it belongs in the caption. release_sheet_constants()
#      finds those and the builder hoists them out.
#   2. Columns are ordered scientifically: what the row identifies, then the result, then
#      supporting quantities, then provenance. Machine order is insertion order, which is
#      not the same thing.
#   3. Numeric formats are chosen per column from the values actually present, so a
#      p-value column spanning 1e-30 renders in scientific notation and one spanning
#      0.05-1 renders as a decimal.

# --------------------------------------------------------------------------------------
# display labels
# --------------------------------------------------------------------------------------
# Units belong in the label where the quantity has them. Statistical symbols follow the
# manuscript convention: italic P is not available in a plain header, so "P value" is used
# rather than "p-value" or "pval".

RELEASE_COLUMN_LABELS <- c(
  # --- identity / design ------------------------------------------------------------
  sample_id                            = "Sample ID",
  acquisition_run_name                 = "Acquisition run name",
  raw_file_basename                    = "Raw acquisition file",
  raw_file_original_path               = "Raw acquisition file, original path",
  raw_file_extension                   = "Raw file format",
  AnimalID                             = "Animal ID",
  hemisphere                           = "Hemisphere",
  hemisphere_evidence_source           = "Hemisphere evidence",
  historical_hemisphere_label          = "Hemisphere label (historical)",
  legacy_replicate_group               = "ReplicateGroup (legacy)",
  sample_class                         = "Sample class",
  sample_class_historical_alias        = "Sample class alias (historical)",
  sample_class_historical_group_label  = "Sample class group label (historical)",
  condition_code                       = "Condition code",
  condition                            = "Condition",
  pairing_status                       = "Pairing",
  treatment                            = "Treatment",
  collection_plate                     = "Collection plate",
  plate_sample_number                  = "Plate sample number",
  injection_index                      = "Injection order",
  plate_set_token                      = "Plate set",
  well_position                        = "Well position",
  acquisition_run_trailing_id          = "Run serial",
  instrument_alias_token               = "Instrument alias in run name",
  instrument_model                     = "Instrument model",
  acquisition_date                     = "Acquisition date",
  acquisition_date_token               = "Acquisition date token",
  lc_and_method_token                  = "LC/method token in run name",
  analysis_unit_primary                = "Analysis unit",
  analysis_unit                        = "Analysis unit",
  included_in_animal_level             = "Included in animal-level analysis",
  animal_level_sample_id               = "Animal-level unit",
  animal_level_column_name             = "Animal-level unit",
  measurement_level_matrix_column      = "Measurement-level matrix column",
  source_assignment_status             = "Assignment status",
  legacy_sample_number                 = "Sample number (legacy)",
  legacy_shortname                     = "Short name (legacy)",
  legacy_group2_hemisphere_coded       = "group2, hemisphere-coded (legacy)",
  phenotype_within_unit                = "Phenotype",
  aggregation_policy                   = "Aggregation policy",
  n_hemisphere_measurements            = "Hemisphere measurements (n)",
  hemispheres_present                  = "Hemispheres present",
  source_sample_ids                    = "Source sample IDs",
  source_raw_file_basenames            = "Source raw acquisition files",
  left_source_sample                   = "Left hemisphere sample",
  right_source_sample                  = "Right hemisphere sample",

  # --- sample-class correction --------------------------------------------------------
  original_sample_class                = "Sample class, original record",
  analysis_sample_class                = "Sample class, as analysed",
  sample_class_plate_layout_implied    = "Sample class implied by plate layout",
  sample_class_corrected               = "Sample class corrected",
  sample_class_correction_status       = "Correction status",
  sample_class_correction_method       = "Correction method",
  sample_class_correction_provenance   = "Correction provenance",
  sample_class_correction_reference_sha256 = "Correction record SHA-256",
  sample_class_correction_date         = "Correction date",

  # --- contrasts ----------------------------------------------------------------------
  contrast_family                      = "Contrast family",
  numerator_condition                  = "Numerator condition",
  denominator_condition                = "Denominator condition",
  numerator_code                       = "Numerator code",
  denominator_code                     = "Denominator code",
  numerator_phenotype                  = "Numerator phenotype",
  denominator_phenotype                = "Denominator phenotype",
  canonical_comparison                 = "Comparison",
  canonical_contrast                   = "Contrast",
  historical_comparison_alias          = "Comparison alias (historical)",
  n_numerator_animals                  = "Numerator animals (n)",
  n_denominator_animals                = "Denominator animals (n)",
  primary_or_secondary                 = "Primary or secondary",
  statistical_model                    = "Statistical model",
  gsea_ranking_statistic               = "GSEA ranking statistic",
  interpretation                       = "Interpretation",
  differential_statistics_file         = "Differential statistics file",
  canonical_split_forward_path         = "Canonical differential table",
  canonical_mapped_forward_path        = "Canonical mapped table",
  canonical_enrichment_dir             = "Canonical enrichment folder",

  # --- proteins and differential statistics -------------------------------------------
  protein_group_id                     = "Protein group (UniProt entry names)",
  uniprot_accession                    = "UniProt accession",
  gene_symbol                          = "Gene symbol",
  protein_description                  = "Protein description",
  id_mapping_status                    = "Identifier mapping status",
  effect_size_sd_units                 = "Standardized abundance difference (SD units)",
  average_standardized_abundance       = "Mean standardized abundance",
  moderated_t                          = "Moderated t",
  P.Value                              = "P value",
  adj.P.Val                            = "Adjusted P value (BH)",
  B_log_odds                           = "B (log odds)",
  neg_log10_P_value                    = "-log10 P value",
  signed_neg_log10_P_value             = "Signed -log10 P value",
  significant_fdr_0_05                 = "Significant (FDR < 0.05)",
  effect_size_units                    = "Effect size units",
  effect_size_definition               = "Effect size definition",
  effect_size_source_column            = "Effect size source column",
  source_statistic_field               = "Source statistic field",
  significance_definition              = "Significance definition",
  n_proteins_tested                    = "Proteins tested (n)",
  n_significant_fdr_0_05               = "Significant, FDR < 0.05 (n)",
  n_significant_higher_in_numerator    = "Higher in numerator (n)",
  n_significant_higher_in_denominator  = "Higher in denominator (n)",
  n_proteins_tested_mapped_only        = "Proteins tested, mapped only (n)",
  n_significant_fdr_0_05_mapped_only   = "Significant, mapped only (n)",
  n_significant_higher_in_numerator_mapped_only   = "Higher in numerator, mapped only (n)",
  n_significant_higher_in_denominator_mapped_only = "Higher in denominator, mapped only (n)",
  canonical_source_path                = "Canonical source file",
  canonical_source_sha256              = "Canonical source SHA-256",

  # --- enrichment ---------------------------------------------------------------------
  analysis                             = "Analysis",
  ontology                             = "Ontology",
  query_list                           = "Query list",
  term_id                              = "Term ID",
  term_name                            = "Term name",
  set_size                             = "Gene set size",
  enrichment_score                     = "Enrichment score",
  NES                                  = "NES",
  p_value                              = "P value",
  adjusted_p_value                     = "Adjusted P value (BH)",
  q_value                              = "q value",
  gene_count                           = "Genes in term (n)",
  gene_ratio                           = "Gene ratio",
  background_ratio                     = "Background ratio",
  fold_enrichment                      = "Fold enrichment",
  rich_factor                          = "Rich factor",
  z_score                              = "z score",
  rank_at_max                          = "Rank at maximum",
  core_enrichment_definition           = "Core enrichment definition",
  n_core_enrichment_genes              = "Core enrichment genes (n)",
  direction                            = "Direction",
  rank_statistic                       = "Ranking statistic",
  rank_source_column                   = "Ranking source column",
  analysis_role                        = "Analysis role",
  gene_identifier_type                 = "Gene identifier type",
  ora_universe_definition              = "ORA universe",
  rows_exported                        = "Rows exported (n)",
  term_count_recorded_by_run           = "Terms recorded by run (n)",
  agrees_with_run_record               = "Agrees with run record",
  note                                 = "Note",
  parameter                            = "Parameter",
  value                                = "Value",

  # --- EWCE ---------------------------------------------------------------------------
  ewce_analysis_type                   = "EWCE analysis type",
  protein_list_definition              = "Protein list",
  cell_type                            = "Cell type",
  annotation_level                     = "Annotation level",
  top_n                                = "Top-N proteins",
  n_hits                               = "Hits (n)",
  n_background                         = "Background (n)",
  fold_change                          = "Fold change (observed/null)",
  sd_from_mean                         = "SD from null mean",
  adjusted_p_value_within_target       = "Adjusted P value, within list (BH)",
  adjusted_p_value_global              = "Adjusted P value, global (BH)",
  significant_global_fdr_0_05          = "Significant, global FDR < 0.05",
  ewce_direction                       = "EWCE direction",
  bootstrap_replicates                 = "Bootstrap replicates (n)",
  is_primary_setting                   = "Primary setting",

  # --- secondary analyses -------------------------------------------------------------
  analysis_id                          = "Analysis ID",
  analysis_name                        = "Analysis",
  analysis_type                        = "Analysis type",
  comparisons                          = "Comparisons",
  condition_scope                      = "Condition scope",
  ranking_statistic                    = "Ranking statistic",
  n_per_group                          = "n per group",
  included_in_primary_contrast_manifest = "In primary contrast set",
  canonical_source                     = "Canonical source",
  generating_script                    = "Generating script",
  caveat                               = "Caveat",

  # --- figures ------------------------------------------------------------------------
  figure                               = "Figure",
  panel                                = "Panel",
  panel_label                          = "Panel label",
  panel_title                          = "Panel title",
  source_table                         = "Source table in this release",
  source_rows_filter                   = "Row filter",
  statistical_unit                     = "Statistical unit",
  original_statistical_unit            = "Original statistical unit",
  corrected_source_script              = "Corrected source script",
  final_revision_script                = "Final revision script",
  revision_status                      = "Revision status",
  displayed_statistics_changed         = "Displayed statistics changed",
  biological_interpretation_changed    = "Interpretation changed",
  key_numbers_original                 = "Key numbers, original",
  key_numbers_revised                  = "Key numbers, revised",
  panel_source_data_release_path       = "Panel source data in this release",
  panel_source_data_original_filename  = "Panel source data, original filename",
  panel_source_data_canonical_path     = "Panel source data, canonical path",
  panel_source_data_sha256             = "Panel source data SHA-256",
  interpretation_note                  = "Interpretation note",

  # --- metadata status / software ------------------------------------------------------
  field                                = "Field",
  status                               = "Status",
  current_source                       = "Current source",
  evidence                             = "Evidence",
  component                            = "Component",
  category                             = "Category",
  version                              = "Version",
  recorded_by                          = "Recorded by",
  evidence_path                        = "Evidence path",
  evidence_sha256                      = "Evidence SHA-256",
  applies_to_stage                     = "Applies to",
  notes                                = "Notes",

  # --- README -------------------------------------------------------------------------
  item                                 = "Item",
  statement                            = "Statement"
)

#' Display header for each machine column name.
#'
#' Falls back to a readable transformation rather than the raw name, so a column added
#' upstream without a label here still renders acceptably instead of leaking snake_case.
release_display_labels <- function(cols) {
  cols <- as.character(cols)
  out <- unname(RELEASE_COLUMN_LABELS[cols])
  missing <- is.na(out)
  if (any(missing)) {
    fallback <- gsub("_", " ", cols[missing])
    fallback <- sub("^(.)", "\\U\\1", fallback, perl = TRUE)
    out[missing] <- fallback
  }
  out
}

#' Machine columns that have no label yet -- used by a contract so the map cannot rot.
release_unlabelled_columns <- function(cols) {
  cols <- unique(as.character(cols))
  cols[!cols %in% names(RELEASE_COLUMN_LABELS)]
}

# --------------------------------------------------------------------------------------
# sheet captions
# --------------------------------------------------------------------------------------
# Supplementary-table numbering is assigned here rather than derived from sheet order, so
# inserting a sheet cannot silently renumber a table an editor has already cited.

RELEASE_SHEET_META <- list(
  README = list(
    number = NA_character_, title = "Read me first",
    caption = paste("How to read this workbook: the statistical unit, the sample size, what",
                    "changed relative to the original submission, and the terminology used",
                    "throughout. Read this before any other sheet.")),
  Sample_Metadata = list(
    number = "S1", title = "Acquisition-level sample metadata",
    caption = paste("One row per acquisition. Includes animal, hemisphere, sample class,",
                    "condition, collection plate, acquisition run identity, and the",
                    "sample-class correction record for the six corrected acquisitions.")),
  Animal_Level_Metadata = list(
    number = "S2", title = "Animal-level inferential units",
    caption = paste("One row per animal x sample class. These are the units all inference",
                    "is computed on; each is the mean of its left and right hemisphere",
                    "measurements.")),
  Primary_Contrasts = list(
    number = "S3", title = "Primary contrast definitions",
    caption = paste("The 12 primary within-sample-class comparisons, their orientation,",
                    "group sizes and interpretation.")),
  Differential_Summary = list(
    number = "S4", title = "Differential abundance summary by comparison",
    caption = paste("Counts of tested and significant proteins per primary comparison.",
                    "Counts are given over all tested protein groups and over the mapped",
                    "subset, because the manuscript figures were drawn from the mapped",
                    "subset.")),
  Differential_Proteins = list(
    number = "S5", title = "Differential abundance, protein level",
    caption = paste("One row per protein group per primary comparison. The effect size is a",
                    "standardized abundance difference in SD units, not a log2 fold change.",
                    "Positive values are higher in the numerator condition.")),
  GSEA_GO_BP = list(
    number = "S6", title = "Gene set enrichment analysis, GO biological process",
    caption = paste("Canonical GSEA, ranked by the moderated t statistic. The complete",
                    "tested result is given, not a significant-only subset: filter on the",
                    "adjusted P value.")),
  GSEA_KEGG = list(
    number = "S7", title = "Gene set enrichment analysis, KEGG",
    caption = paste("Canonical GSEA, ranked by the moderated t statistic. Complete tested",
                    "result.")),
  ORA_GO_BP = list(
    number = "S8", title = "Over-representation analysis, GO biological process",
    caption = paste("Over-representation of GO biological process terms across four",
                    "query-list definitions. Complete tested result.")),
  GSEA_log2FC_sensitivity = list(
    number = "S9", title = "Gene set enrichment analysis, effect-size-ranked (sensitivity)",
    caption = paste("SENSITIVITY ANALYSIS ONLY. GSEA ranked by the standardized abundance",
                    "difference instead of the moderated t statistic. Not the canonical",
                    "result and must not be reported as such.")),
  EWCE = list(
    number = "S10", title = "Cell-type enrichment (EWCE)",
    caption = paste("Expression-weighted cell-type enrichment across all top-N and",
                    "annotation-level settings. The primary setting is flagged; other rows",
                    "are sensitivity settings.")),
  Enrichment_Coverage = list(
    number = "S11", title = "Enrichment coverage audit",
    caption = paste("Rows exported per comparison and analysis, cross-checked against the",
                    "term counts the canonical run recorded. Empty query lists are a result",
                    "and are recorded as such, not a missing file.")),
  Enrichment_Parameters = list(
    number = "S12", title = "Enrichment and EWCE parameters",
    caption = "Parameters as recorded by the canonical enrichment and EWCE runs."),
  Secondary_Analyses = list(
    number = "S13", title = "Secondary and excluded analyses",
    caption = paste("Analyses that are deliberately not part of the 12 primary contrasts,",
                    "with the reason and the inferential caveat for each.")),
  Figure_Source_Map = list(
    number = "S14", title = "Figure panel provenance",
    caption = paste("Every proteomics manuscript panel mapped to the table it is drawn",
                    "from, the row filter, the statistical unit, the generating script and",
                    "the interpretation caveat.")),
  Metadata_Field_Status = list(
    number = "S15", title = "Metadata field status",
    caption = paste("Per-field record of what is verified, what needs standardisation and",
                    "what is missing. Missing values are left blank rather than inferred.")),
  Software_Versions = list(
    number = "S16", title = "Software and reference database versions",
    caption = paste("Versions as recorded by the runs themselves, or UNKNOWN with the source",
                    "that would supply them. Scoped by the stage each applies to, so a",
                    "superseded version cannot be read as the version behind the reported",
                    "result."))
)

#' Preferred leading column order per sheet.
#'
#' Only the leading columns are named; anything unnamed keeps its relative order and is
#' appended. So adding a column upstream cannot drop it from the workbook -- it simply
#' lands at the right-hand end until it is given a position here.
RELEASE_SHEET_COLUMN_ORDER <- list(
  Sample_Metadata = c("sample_id", "AnimalID", "hemisphere", "sample_class", "condition",
                      "pairing_status", "treatment", "collection_plate", "well_position",
                      "injection_index", "acquisition_date", "animal_level_sample_id",
                      "raw_file_basename"),
  Animal_Level_Metadata = c("AnimalID", "sample_class", "condition", "pairing_status",
                            "treatment", "collection_plate", "animal_level_column_name",
                            "n_hemisphere_measurements", "hemispheres_present",
                            "left_source_sample", "right_source_sample"),
  Primary_Contrasts = c("canonical_comparison", "sample_class", "contrast_family",
                        "numerator_condition", "denominator_condition",
                        "n_numerator_animals", "n_denominator_animals", "interpretation"),
  Differential_Summary = c("canonical_comparison", "sample_class", "contrast_family",
                           "numerator_condition", "denominator_condition",
                           "n_numerator_animals", "n_denominator_animals",
                           "n_proteins_tested", "n_significant_fdr_0_05",
                           "n_significant_higher_in_numerator",
                           "n_significant_higher_in_denominator",
                           "n_proteins_tested_mapped_only",
                           "n_significant_fdr_0_05_mapped_only",
                           "n_significant_higher_in_numerator_mapped_only",
                           "n_significant_higher_in_denominator_mapped_only"),
  Differential_Proteins = c("canonical_comparison", "sample_class", "contrast_family",
                            "numerator_condition", "denominator_condition",
                            "gene_symbol", "uniprot_accession", "protein_group_id",
                            "protein_description", "effect_size_sd_units", "moderated_t",
                            "P.Value", "adj.P.Val", "significant_fdr_0_05",
                            "average_standardized_abundance", "B_log_odds"),
  GSEA_GO_BP = c("canonical_comparison", "sample_class", "contrast_family", "term_id",
                 "term_name", "NES", "enrichment_score", "p_value", "adjusted_p_value",
                 "q_value", "set_size", "n_core_enrichment_genes", "direction"),
  GSEA_KEGG = c("canonical_comparison", "sample_class", "contrast_family", "term_id",
                "term_name", "NES", "enrichment_score", "p_value", "adjusted_p_value",
                "q_value", "set_size", "n_core_enrichment_genes", "direction"),
  ORA_GO_BP = c("canonical_comparison", "sample_class", "contrast_family", "query_list",
                "term_id", "term_name", "gene_count", "gene_ratio", "background_ratio",
                "fold_enrichment", "p_value", "adjusted_p_value", "q_value", "direction"),
  GSEA_log2FC_sensitivity = c("canonical_comparison", "sample_class", "contrast_family",
                              "ontology", "term_id", "term_name", "NES", "p_value",
                              "adjusted_p_value", "q_value", "set_size", "direction"),
  EWCE = c("sample_class", "ewce_analysis_type", "canonical_comparison",
           "protein_list_definition", "direction", "cell_type", "annotation_level",
           "top_n", "sd_from_mean", "fold_change", "p_value",
           "adjusted_p_value_global", "significant_global_fdr_0_05", "is_primary_setting"),
  Figure_Source_Map = c("panel_label", "figure", "panel", "panel_title", "analysis",
                        "statistical_unit", "primary_or_secondary", "source_table",
                        "source_rows_filter", "key_numbers_revised", "interpretation_note"),
  Software_Versions = c("component", "category", "version", "status", "applies_to_stage",
                        "recorded_by", "notes")
)

#' Reorder a sheet's columns: named leading columns first, the rest in their existing order.
release_order_columns <- function(df, sheet) {
  lead <- RELEASE_SHEET_COLUMN_ORDER[[sheet]]
  if (is.null(lead)) return(df)
  lead <- lead[lead %in% names(df)]
  df[, c(lead, setdiff(names(df), lead)), drop = FALSE]
}

# --------------------------------------------------------------------------------------
# sheet-level constants
# --------------------------------------------------------------------------------------

#' Columns whose value is identical on every row, with that value.
#'
#' A quantity that never varies is a property of the table, not a variable in it. Hoisting
#' these into the caption is the difference between a supplementary table and a database
#' dump -- and on the 64,188-row differential sheet it removes four columns of repetition.
#'
#' `keep` names columns that stay in the sheet even when constant, because a reader filters
#' or joins on them.
release_sheet_constants <- function(df, keep = character(0), min_rows = 2L) {
  if (nrow(df) < min_rows) {
    return(list(constants = character(0), values = character(0)))
  }
  candidate <- setdiff(names(df), keep)
  is_const <- vapply(candidate, function(nm) {
    v <- df[[nm]]
    u <- unique(v)
    length(u) == 1L && !is.na(u[[1]]) && nzchar(trimws(as.character(u[[1]])))
  }, logical(1))
  cols <- candidate[is_const]
  list(constants = cols,
       values = vapply(cols, function(nm) as.character(unique(df[[nm]])[[1]]), character(1)))
}

# --------------------------------------------------------------------------------------
# numeric formats
# --------------------------------------------------------------------------------------

#' Choose an Excel number format for one column from the values actually present.
#'
#' Probability-like columns spanning very small values render in scientific notation;
#' the same column in a table where nothing is small renders as a plain decimal. Integer
#' columns render without a decimal point. Anything non-numeric returns NA and is left
#' as text.
release_numfmt_for <- function(values, column_name = "") {
  if (!is.numeric(values)) return(NA_character_)
  finite <- values[is.finite(values)]
  if (!length(finite)) return(NA_character_)

  if (all(finite == round(finite)) && max(abs(finite)) < 1e15) return("0")

  smallest <- min(abs(finite[finite != 0]), Inf)
  probability_like <- grepl("p_value|P[.]Value|adj[.]P[.]Val|adjusted_p|q_value",
                            column_name)
  if (probability_like) {
    if (is.finite(smallest) && smallest < 1e-4) return("0.00E+00")
    return("0.0000")
  }
  if (is.finite(smallest) && smallest < 1e-4 && max(abs(finite)) < 1e-2) return("0.00E+00")
  if (max(abs(finite)) >= 1000) return("#,##0.00")
  "0.000"
}

# --------------------------------------------------------------------------------------
# display names for controlled-vocabulary values
# --------------------------------------------------------------------------------------
# These are for PROSE and markdown tables only. Data cells keep the controlled vocabulary
# (`mcherry`, `paired_cno`, `cfos_paired_cno_over_cfos_paired_veh`) because that is what
# joins the tables together and what the contracts assert on. A reader-facing sentence or
# table heading, though, should read "cFos: paired-CNO vs paired-VEH" -- the underscored
# machine key in running text is exactly what makes a supplementary section look like a
# database dump. README_DATA.md documents both namings side by side.

RELEASE_SAMPLE_CLASS_DISPLAY <- c(
  mcherry  = "mCherry",
  neuropil = "Neuropil",
  cfos     = "cFos",
  neuron   = "Neuron"
)

RELEASE_CONDITION_DISPLAY <- c(
  paired_cno   = "paired-CNO",
  paired_veh   = "paired-VEH",
  unpaired_cno = "unpaired-CNO",
  unpaired_veh = "unpaired-VEH"
)

release_sample_class_display <- function(x) {
  x <- as.character(x)
  out <- unname(RELEASE_SAMPLE_CLASS_DISPLAY[x])
  ifelse(is.na(out), x, out)
}

release_condition_display <- function(x) {
  x <- as.character(x)
  out <- unname(RELEASE_CONDITION_DISPLAY[x])
  ifelse(is.na(out), x, out)
}

#' Reader-facing label for one comparison, e.g. "cFos: paired-CNO vs paired-VEH".
release_comparison_label <- function(sample_class, numerator_condition,
                                     denominator_condition) {
  paste0(release_sample_class_display(sample_class), ": ",
         release_condition_display(numerator_condition), " vs ",
         release_condition_display(denominator_condition))
}

#' Column display width, bounded so one long free-text field cannot dominate the sheet.
release_column_width <- function(values, header, min_width = 9L, max_width = 52L) {
  n <- nchar(header)
  if (length(values)) {
    sample_values <- utils::head(as.character(values[!is.na(values)]), 400L)
    if (length(sample_values)) {
      n <- max(n, stats::quantile(nchar(sample_values), 0.95, names = FALSE))
    }
  }
  max(min_width, min(max_width, ceiling(n) + 2L))
}
