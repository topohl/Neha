# Display labels, captions, column choice and number formats for Supplementary Table 2
# (editor_source_data/Supplementary_Table_2_proteomics.xlsx).
#
# The release has every table twice:
#   *.tsv / *.tsv.gz     snake_case column names and controlled-vocabulary values, one table
#                        per file; the tests and the data dictionary use these names
#   the .xlsx workbook   display headers, display labels for keys (compartments, conditions,
#                        comparisons, directions), lettered sub-tables 2a-2q, captions
#
# This file only maps between the two and holds no data. Requires release_validation.R to be
# sourced first. Numbers are never changed: number formats change how a value is shown,
# not the value.

# --------------------------------------------------------------------------------------
# display labels for column headers
# --------------------------------------------------------------------------------------
# Units belong in the header where the quantity has them. "P value" rather than "p-value" or
# "pval". Headers never contain an underscore.

RELEASE_COLUMN_LABELS <- c(
  # --- comparisons ----------------------------------------------------------------------
  comparison_label                     = "Comparison",
  canonical_comparison                 = "Comparison key",
  contrast_family                      = "Comparison type",
  numerator_condition                  = "Group 1",
  denominator_condition                = "Group 2",
  n_numerator_animals                  = "Group 1 animals (n)",
  n_denominator_animals                = "Group 2 animals (n)",
  n_proteins_tested                    = "Proteins tested (n)",
  n_significant_fdr_0_05               = "FDR < 0.05 (n)",
  n_significant_higher_in_numerator    = "FDR < 0.05, higher in group 1 (n)",
  n_significant_higher_in_denominator  = "FDR < 0.05, lower in group 1 (n)",
  n_proteins_tested_mapped_only        = "Annotated proteins tested (n)",
  n_significant_fdr_0_05_mapped_only   = "Annotated, FDR < 0.05 (n)",
  n_significant_higher_in_numerator_mapped_only   = "Annotated, FDR < 0.05, higher in group 1 (n)",
  n_significant_higher_in_denominator_mapped_only = "Annotated, FDR < 0.05, lower in group 1 (n)",

  # --- samples and animals ----------------------------------------------------------------
  sample_id                            = "Sample ID",
  raw_file_basename                    = "Raw file",
  AnimalID                             = "Animal ID",
  hemisphere                           = "Hemisphere",
  sample_class                         = "Compartment",
  condition                            = "Condition",
  pairing_status                       = "Pairing",
  treatment                            = "Treatment",
  collection_plate                     = "Collection plate",
  well_position                        = "Well",
  injection_index                      = "Injection order",
  acquisition_date                     = "Acquisition date",
  animal_level_sample_id               = "Animal-level profile",
  original_sample_class                = "Compartment as recorded",
  sample_class_corrected               = "Compartment label corrected",
  animal_level_column_name             = "Animal-level profile",
  left_source_sample                   = "Left hemisphere sample",
  right_source_sample                  = "Right hemisphere sample",
  n_hemisphere_measurements            = "Hemisphere runs averaged (n)",

  # --- proteins ---------------------------------------------------------------------------
  protein_group_id                     = "Protein group",
  uniprot_accession                    = "UniProt accession",
  gene_symbol                          = "Gene symbol",
  protein_description                  = "Protein description",
  id_mapping_status                    = "Annotated",
  effect_size_sd_units                 = "Effect size (SD units)",
  average_standardized_abundance       = "Mean standardized abundance (SD units)",
  moderated_t                          = "Moderated t",
  P.Value                              = "P value",
  adj.P.Val                            = "Adjusted P value (BH)",
  B_log_odds                           = "B (log odds)",
  significant_fdr_0_05                 = "FDR < 0.05",
  direction                            = "Direction",

  # --- gene sets --------------------------------------------------------------------------
  analysis                             = "Analysis",
  ontology                             = "Gene-set collection",
  query_list                           = "Query list",
  term_id                              = "Term ID",
  term_name                            = "Term",
  set_size                             = "Gene-set size",
  enrichment_score                     = "Enrichment score",
  NES                                  = "NES",
  p_value                              = "P value",
  adjusted_p_value                     = "Adjusted P value (BH)",
  q_value                              = "q value",
  gene_count                           = "Query proteins in set (n)",
  gene_ratio                           = "Gene ratio",
  background_ratio                     = "Background ratio",
  fold_enrichment                      = "Fold enrichment",
  rich_factor                          = "Rich factor",
  z_score                              = "z score",
  rank_at_max                          = "Rank at maximum",
  n_core_enrichment_genes              = "Leading-edge proteins (n)",
  core_enrichment_gene_symbols         = "Leading-edge genes",
  rank_statistic                       = "Ranking statistic",
  analysis_role                        = "Analysis role",
  n_perm_simple                        = "fgsea nPermSimple",

  # --- enrichment settings and coverage (table 2m) ----------------------------------------
  gene_set_source                      = "Gene-set source",
  set_size_limits                      = "Gene-set size limits",
  permutations                         = "Permutations",
  seed                                 = "Random seed",
  ora_background                       = "ORA background (n)",
  multiple_testing                     = "Multiple testing",
  rows_exported                        = "Gene sets reported (n)",
  gsea_untested_gene_sets              = "Gene sets without a P value (n)",
  note                                 = "Note",

  # --- EWCE -------------------------------------------------------------------------------
  ewce_analysis_type                   = "EWCE analysis",
  protein_list_definition              = "Protein list",
  top_n                                = "Proteins in list (n)",
  n_hits                               = "Hits (n)",
  annotation_level                     = "Annotation level",
  cell_type                            = "Cell type",
  n_background                         = "Background proteins (n)",
  fold_change                          = "Fold change (observed/null)",
  sd_from_mean                         = "SD from null mean",
  adjusted_p_value_within_target       = "Adjusted P value, within list (BH)",
  adjusted_p_value_global              = "Adjusted P value, global (BH)",
  significant_global_fdr_0_05          = "FDR < 0.05, global",
  ewce_direction                       = "EWCE direction",
  bootstrap_replicates                 = "Bootstrap replicates (n)",
  is_primary_setting                   = "Primary setting",

  # --- secondary analyses (table 2o) ------------------------------------------------------
  analysis_name                        = "Analysis",
  purpose                              = "Purpose",
  model                                = "Model",
  n_per_group                          = "n",
  results_location                     = "Where the results are",
  key_result                           = "Key result",

  # --- figure panels (table 2p) -----------------------------------------------------------
  panel_label                          = "Panel",
  panel_title                          = "Title",
  statistical_unit                     = "Unit of analysis",
  panel_source                         = "Source",
  rows_used                            = "Rows used",
  key_numbers                          = "Key numbers",

  # --- software (table 2q) ----------------------------------------------------------------
  component                            = "Component",
  version                              = "Version",
  used_for                             = "Used for",

  # --- README and Column definitions ------------------------------------------------------
  item                                 = "Item",
  statement                            = "Statement",
  table                                = "Table",
  column                               = "Column",
  definition                           = "Definition",
  units                                = "Units"
)

#' Display header for each machine column name.
#'
#' Falls back to a readable transformation rather than the raw name; stage 07 stops on any
#' shown column without a label (release_unlabelled_columns()), so the fallback is only a
#' guard.
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

#' Machine columns with no display label yet (stage 07 stops if a shown column has none).
release_unlabelled_columns <- function(cols) {
  cols <- unique(as.character(cols))
  cols[!cols %in% names(RELEASE_COLUMN_LABELS)]
}

# --------------------------------------------------------------------------------------
# sheets: letter, tab name, title, caption
# --------------------------------------------------------------------------------------
# The workbook is Supplementary Table 2 as a whole; its sheets are the lettered sub-tables
# 2a-2q, in this order after the Contents sheet. Keys are internal and never shown. A caption
# set to NA is composed by stage 07 from the data (counts that must not go stale).

RELEASE_SUPPLEMENTARY_TABLE <- "2"

RELEASE_SHEET_META <- list(
  README = list(
    letter = NA_character_, tab = "README", title = "Read me",
    caption = "What the tables contain and how to read them."),
  Column_Definitions = list(
    letter = NA_character_, tab = "Column definitions", title = "Column definitions",
    caption = paste("One row per column of tables 2a to 2q: the header as shown, what the",
                    "column holds and its units.")),
  Samples = list(
    letter = "2a", tab = "2a Samples", title = "Samples (acquisitions)",
    caption = paste("One row per LC-MS run: 12 animals x 4 compartments x 2 hemispheres. The",
                    "left and right runs of each animal and compartment were averaged before",
                    "any test (table 2b).")),
  Animal_Profiles = list(
    letter = "2b", tab = "2b Animal profiles", title = "Animal-level profiles",
    caption = paste("One row per animal and compartment: the mean of the left and right",
                    "hemisphere runs. These 48 profiles are the units of every test, with",
                    "n = 3 animals per group in each compartment.")),
  Comparisons = list(
    letter = "2c", tab = "2c Comparisons", title = "Comparisons",
    caption = paste("The 12 comparisons, three per compartment. Group 1 is compared with",
                    "group 2, so a positive effect size means higher in group 1. Counts are",
                    "protein groups at FDR < 0.05 (Benjamini-Hochberg), over all tested",
                    "groups and over the annotated groups (those with a mouse UniProt",
                    "accession).")),
  Significant_Proteins = list(
    letter = "2d", tab = "2d Significant proteins", title = "Significant proteins",
    caption = NA_character_),
  All_Proteins = list(
    letter = "2e", tab = "2e All proteins", title = "All proteins",
    caption = paste("Every protein group tested, in each of the 12 comparisons, sorted by",
                    "comparison and P value. limma, one model per compartment; P values",
                    "adjusted with Benjamini-Hochberg within each comparison. The effect size",
                    "is a difference in standardized abundance in SD units (group 1 minus",
                    "group 2), not a log2 fold change.")),
  Significant_Gene_Sets = list(
    letter = "2f", tab = "2f Significant gene sets", title = "Significant gene sets",
    caption = NA_character_),
  GSEA_GO_BP = list(
    letter = "2g", tab = "2g GSEA GO-BP", title = "GSEA, GO biological process",
    caption = paste("Every GO biological process gene set GSEA tested in each comparison,",
                    "not only the significant ones: filter on the adjusted P value. Proteins",
                    "were ranked by the moderated t. Gene sets fgsea could not test are in",
                    "table 2l. P values at fgsea's lower bound (1E-10) are shown as",
                    "<=1E-10.")),
  GSEA_GO_CC = list(
    letter = "2h", tab = "2h GSEA GO-CC", title = "GSEA, GO cellular component",
    caption = paste("Every GO cellular component gene set GSEA tested in each comparison,",
                    "not only the significant ones: filter on the adjusted P value. Proteins",
                    "were ranked by the moderated t. Gene sets fgsea could not test are in",
                    "table 2l. P values at fgsea's lower bound (1E-10) are shown as",
                    "<=1E-10.")),
  GSEA_KEGG = list(
    letter = "2i", tab = "2i GSEA KEGG", title = "GSEA, KEGG pathways",
    caption = paste("Every KEGG pathway GSEA tested in each comparison, not only the",
                    "significant ones: filter on the adjusted P value. Proteins were ranked",
                    "by the moderated t. The number of pathways fgsea could not test is in",
                    "table 2m. P values at fgsea's lower bound (1E-10) are shown as",
                    "<=1E-10.")),
  ORA_GO_BP = list(
    letter = "2j", tab = "2j ORA GO-BP", title = "ORA, GO biological process",
    caption = NA_character_),
  ORA_GO_CC = list(
    letter = "2k", tab = "2k ORA GO-CC", title = "ORA, GO cellular component",
    caption = NA_character_),
  GSEA_Untested_Sets = list(
    letter = "2l", tab = "2l No P value", title = "Gene sets without a P value",
    caption = paste("GO gene sets within the size limits that fgsea could not give a P value",
                    "in the primary GSEA because the ranking was strongly one-sided, one row",
                    "per comparison and set. They are absent from tables 2g and 2h rather than",
                    "non-significant. KEGG counts are in table 2m.")),
  Enrichment_Settings = list(
    letter = "2m", tab = "2m Enrichment settings", title = "Enrichment settings and coverage",
    caption = paste("Settings of the GSEA and ORA in tables 2f to 2l, and the number of gene",
                    "sets reported for each comparison and collection. The settings and",
                    "coverage of the effect-size-ranked GSEA and of the alternative ORA list",
                    "are in the data release (enrichment/enrichment_run_parameters.tsv,",
                    "enrichment/enrichment_coverage.tsv).")),
  EWCE = list(
    letter = "2n", tab = "2n EWCE", title = "Cell-type enrichment (EWCE)",
    caption = NA_character_),
  Secondary_Analyses = list(
    letter = "2o", tab = "2o Secondary analyses", title = "Secondary analyses",
    caption = paste("Analyses reported in the paper in addition to the 12 comparisons. None of",
                    "them is part of the primary results.")),
  Figure_Panels = list(
    letter = "2p", tab = "2p Figure panels", title = "Figure panels",
    caption = paste("Source of each proteomics panel of Figure 2 and Supplementary Figure 4:",
                    "the table of this workbook and the Source Data file it is drawn from,",
                    "the rows used and the key numbers.")),
  Software = list(
    letter = "2q", tab = "2q Software", title = "Software",
    caption = paste("Software and reference data behind the reported results, with the",
                    "versions the analysis runs recorded."))
)

#' Sheet title as printed in row 1, e.g. "Supplementary Table 2a | Samples (acquisitions)".
release_sheet_title <- function(key) {
  meta <- RELEASE_SHEET_META[[key]]
  paste0("Supplementary Table ", RELEASE_SUPPLEMENTARY_TABLE,
         if (is.na(meta$letter)) "" else substring(meta$letter, 2L), " | ", meta$title)
}

# --------------------------------------------------------------------------------------
# columns per sheet
# --------------------------------------------------------------------------------------
# `show`: the columns of the sheet, in order. `drop`: source columns left out of the
# workbook (they stay in the .tsv files). Stage 07 stops on a source column in neither list,
# so a column added upstream has to be given a place rather than being appended or lost. `keep`:
# columns that stay even when constant or a duplicate (every other constant column moves
# into the line above the header, and a column identical to an earlier one is dropped).
# Columns blank on every row are dropped.

RELEASE_COMMON_KEEP <- c("comparison_label", "canonical_comparison", "sample_class",
                         "AnimalID", "sample_id", "gene_symbol", "uniprot_accession",
                         "protein_group_id", "term_id", "term_name", "query_list",
                         "ontology", "analysis_role", "rank_statistic", "item", "table",
                         "column", "panel_label", "component", "analysis_name")

RELEASE_GSEA_SHOW <- c("comparison_label", "sample_class", "term_id", "term_name",
                       "direction", "NES", "enrichment_score", "set_size",
                       "n_core_enrichment_genes", "rank_at_max", "p_value",
                       "adjusted_p_value", "q_value", "canonical_comparison")
RELEASE_GSEA_DROP <- c("contrast_family", "canonical_contrast", "numerator_condition",
                       "denominator_condition", "n_numerator_animals",
                       "n_denominator_animals", "analysis_unit", "analysis", "ontology",
                       "query_list", "gene_count", "gene_ratio", "background_ratio",
                       "fold_enrichment", "rich_factor", "z_score", "leading_edge",
                       "core_enrichment", "core_enrichment_definition", "rank_statistic",
                       "rank_source_column", "analysis_role", "gene_identifier_type",
                       "ora_universe_definition")
RELEASE_ORA_SHOW <- c("comparison_label", "sample_class", "query_list", "term_id",
                      "term_name", "gene_count", "gene_ratio", "background_ratio",
                      "fold_enrichment", "rich_factor", "z_score", "p_value",
                      "adjusted_p_value", "q_value", "canonical_comparison")
RELEASE_ORA_DROP <- c("contrast_family", "canonical_contrast", "numerator_condition",
                      "denominator_condition", "n_numerator_animals", "n_denominator_animals",
                      "analysis_unit", "analysis", "ontology", "set_size",
                      "enrichment_score", "NES", "rank_at_max", "leading_edge",
                      "core_enrichment", "core_enrichment_definition",
                      "n_core_enrichment_genes", "direction", "rank_statistic",
                      "rank_source_column", "analysis_role", "gene_identifier_type",
                      "ora_universe_definition")

RELEASE_SHEET_COLUMNS <- list(
  README = list(show = c("item", "statement")),
  Column_Definitions = list(show = c("table", "column", "definition", "units")),
  Samples = list(
    show = c("sample_id", "raw_file_basename", "AnimalID", "hemisphere", "sample_class",
             "condition", "pairing_status", "treatment", "collection_plate", "well_position",
             "injection_index", "acquisition_date", "animal_level_sample_id",
             "original_sample_class", "sample_class_corrected"),
    drop = c("acquisition_run_name", "raw_file_original_path", "raw_file_extension",
             "hemisphere_evidence_source", "historical_hemisphere_label",
             "legacy_replicate_group", "sample_class_historical_alias",
             "sample_class_historical_group_label", "condition_code", "plate_sample_number",
             "plate_set_token", "acquisition_run_trailing_id", "instrument_alias_token",
             "instrument_model", "acquisition_date_token", "lc_and_method_token",
             "analysis_unit_primary", "included_in_animal_level",
             "measurement_level_matrix_column", "source_assignment_status",
             "legacy_sample_number", "legacy_shortname", "legacy_group2_hemisphere_coded",
             "analysis_sample_class", "sample_class_plate_layout_implied",
             "sample_class_correction_status", "sample_class_correction_method",
             "sample_class_correction_provenance",
             "sample_class_correction_reference_sha256", "sample_class_correction_date")),
  Animal_Profiles = list(
    show = c("AnimalID", "sample_class", "condition", "pairing_status", "treatment",
             "collection_plate", "animal_level_column_name", "left_source_sample",
             "right_source_sample", "n_hemisphere_measurements"),
    drop = c("sample_class_historical_alias", "condition_code", "hemispheres_present",
             "source_sample_ids", "source_raw_file_basenames", "phenotype_within_unit",
             "aggregation_policy", "analysis_unit")),
  Comparisons = list(
    show = c("comparison_label", "sample_class", "contrast_family", "numerator_condition",
             "denominator_condition", "n_numerator_animals", "n_denominator_animals",
             "n_proteins_tested", "n_significant_fdr_0_05",
             "n_significant_higher_in_numerator", "n_significant_higher_in_denominator",
             "n_proteins_tested_mapped_only", "n_significant_fdr_0_05_mapped_only",
             "n_significant_higher_in_numerator_mapped_only",
             "n_significant_higher_in_denominator_mapped_only", "canonical_comparison"),
    keep = "*"),
  Significant_Proteins = list(
    show = c("comparison_label", "sample_class", "gene_symbol", "protein_description",
             "uniprot_accession", "protein_group_id", "direction", "effect_size_sd_units",
             "moderated_t", "P.Value", "adj.P.Val", "canonical_comparison"),
    drop = c("numerator_condition", "denominator_condition"),
    keep = "direction"),
  All_Proteins = list(
    show = c("comparison_label", "sample_class", "gene_symbol", "protein_description",
             "uniprot_accession", "protein_group_id", "id_mapping_status",
             "effect_size_sd_units", "moderated_t", "P.Value", "adj.P.Val",
             "significant_fdr_0_05", "average_standardized_abundance", "B_log_odds",
             "canonical_comparison"),
    drop = c("contrast_family", "canonical_contrast", "historical_comparison_alias",
             "numerator_condition", "denominator_condition", "n_numerator_animals",
             "n_denominator_animals", "analysis_unit", "primary_or_secondary",
             "source_matrix_gene_symbols", "neg_log10_P_value", "signed_neg_log10_P_value",
             "effect_size_units", "effect_size_source_column", "source_statistic_field")),
  Significant_Gene_Sets = list(
    show = c("comparison_label", "sample_class", "ontology", "term_id", "term_name",
             "direction", "NES", "set_size", "n_core_enrichment_genes", "p_value",
             "adjusted_p_value", "core_enrichment_gene_symbols", "canonical_comparison"),
    drop = c("numerator_condition", "denominator_condition")),
  GSEA_GO_BP = list(show = RELEASE_GSEA_SHOW, drop = RELEASE_GSEA_DROP),
  GSEA_GO_CC = list(show = RELEASE_GSEA_SHOW, drop = RELEASE_GSEA_DROP),
  GSEA_KEGG = list(show = RELEASE_GSEA_SHOW, drop = RELEASE_GSEA_DROP),
  ORA_GO_BP = list(show = RELEASE_ORA_SHOW, drop = RELEASE_ORA_DROP),
  ORA_GO_CC = list(show = RELEASE_ORA_SHOW, drop = RELEASE_ORA_DROP),
  GSEA_Untested_Sets = list(
    show = c("comparison_label", "sample_class", "ontology", "term_id", "term_name",
             "set_size", "n_perm_simple", "canonical_comparison"),
    drop = c("contrast_family", "numerator_condition", "denominator_condition",
             "analysis_unit", "analysis", "rank_statistic", "analysis_role")),
  Enrichment_Settings = list(
    show = c("comparison_label", "sample_class", "analysis", "ontology", "query_list",
             "gene_set_source", "set_size_limits", "rank_statistic", "permutations",
             "seed", "ora_background", "multiple_testing", "rows_exported",
             "gsea_untested_gene_sets", "note", "canonical_comparison"),
    keep = c("analysis", "rows_exported", "gsea_untested_gene_sets")),
  EWCE = list(
    show = c("comparison_label", "sample_class", "ewce_analysis_type",
             "protein_list_definition", "top_n", "n_hits", "annotation_level", "cell_type",
             "fold_change", "sd_from_mean", "p_value", "adjusted_p_value_within_target",
             "adjusted_p_value_global", "significant_global_fdr_0_05", "ewce_direction",
             "is_primary_setting", "n_background", "bootstrap_replicates",
             "canonical_comparison"),
    drop = c("analysis", "contrast_family", "canonical_contrast", "numerator_condition",
             "denominator_condition", "n_numerator_animals", "n_denominator_animals",
             "analysis_unit", "direction", "analysis_role", "primary_or_secondary"),
    keep = c("ewce_analysis_type", "protein_list_definition", "cell_type")),
  Secondary_Analyses = list(
    show = c("analysis_name", "purpose", "model", "n_per_group", "results_location",
             "key_result"),
    keep = "*"),
  Figure_Panels = list(
    show = c("panel_label", "panel_title", "statistical_unit", "panel_source", "rows_used",
             "key_numbers", "note"),
    keep = "*"),
  Software = list(show = c("component", "version", "used_for"), keep = "*")
)

#' Select and order the shown columns of one sheet. Stops on a source column the sheet
#' neither shows nor drops.
release_select_columns <- function(df, key) {
  spec <- RELEASE_SHEET_COLUMNS[[key]]
  if (is.null(spec)) stop("No column specification for sheet ", key, call. = FALSE)
  unplaced <- setdiff(names(df), c(spec$show, spec$drop))
  if (length(unplaced)) {
    stop("Sheet ", key, ": column(s) neither shown nor dropped in RELEASE_SHEET_COLUMNS: ",
         paste(unplaced, collapse = ", "), call. = FALSE)
  }
  absent <- setdiff(spec$show, names(df))
  if (length(absent)) {
    stop("Sheet ", key, ": shown column(s) missing from the source: ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  df[, spec$show, drop = FALSE]
}

#' Columns of a sheet that stay even when constant, blank-free duplicates or repeated.
release_keep_columns <- function(df, key) {
  keep <- RELEASE_SHEET_COLUMNS[[key]]$keep
  if (identical(keep, "*")) return(names(df))
  intersect(names(df), c(RELEASE_COMMON_KEEP, keep))
}

#' Drop columns blank on every row and columns identical to an earlier column, except `keep`.
#' Returns the data and the names dropped.
release_drop_empty_and_duplicate <- function(df, keep = character(0)) {
  blank <- vapply(df, function(v) all(is.na(v) | !nzchar(trimws(as.character(v)))),
                  logical(1))
  dropped <- names(df)[blank & !names(df) %in% keep]
  df <- df[, setdiff(names(df), dropped), drop = FALSE]
  dup <- character(0)
  cols <- names(df)
  for (j in seq_along(cols)) {
    if (cols[[j]] %in% keep || j == 1L) next
    for (i in seq_len(j - 1L)) {
      if (cols[[i]] %in% dup) next
      if (identical(as.character(df[[i]]), as.character(df[[j]]))) {
        dup <- c(dup, cols[[j]])
        break
      }
    }
  }
  list(data = df[, setdiff(names(df), dup), drop = FALSE], dropped = c(dropped, dup))
}

# --------------------------------------------------------------------------------------
# sheet-level constants
# --------------------------------------------------------------------------------------

#' Columns whose value is identical on every row, with that value.
#'
#' The workbook builder moves these into the line above the header. `keep` names columns
#' that stay in the sheet even when constant, because a reader filters or joins on them.
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
# number formats (display only)
# --------------------------------------------------------------------------------------

RELEASE_P_COLUMNS <- c("P.Value", "adj.P.Val", "p_value", "adjusted_p_value", "q_value",
                       "adjusted_p_value_within_target", "adjusted_p_value_global")
RELEASE_THREE_DECIMALS <- c("effect_size_sd_units", "NES", "enrichment_score",
                            "average_standardized_abundance", "sd_from_mean", "fold_change",
                            "fold_enrichment", "rich_factor", "z_score")
RELEASE_TWO_DECIMALS <- c("moderated_t", "B_log_odds")
RELEASE_GSEA_SHEETS <- c("GSEA_GO_BP", "GSEA_GO_CC", "GSEA_KEGG", "Significant_Gene_Sets")

RELEASE_NUMFMT_P <- "[<0.001]0.00E+00;0.000"
# fgsea reports no GSEA P value below 1e-10; such values show as "<=1E-10".
RELEASE_NUMFMT_P_GSEA <- "[<=1E-10]\"<=1E-10\";[<0.001]0.00E+00;0.000"
# 10,000 bootstrap replicates cannot give a P value below 1e-4; 0 shows as "<0.0001".
RELEASE_NUMFMT_P_EWCE <- "[=0]\"<0.0001\";0.0000"

#' Excel number format for one column: P values decimal at or above 0.001 and scientific
#' below; effect sizes 3 decimals; t and B 2 decimals; integers without decimals. NA for
#' non-numeric columns. Values are never changed.
release_numfmt_for <- function(values, column_name = "", sheet = "") {
  if (!is.numeric(values)) return(NA_character_)
  if (column_name %in% RELEASE_P_COLUMNS) {
    if (identical(sheet, "EWCE")) return(RELEASE_NUMFMT_P_EWCE)
    if (column_name == "p_value" && sheet %in% RELEASE_GSEA_SHEETS) return(RELEASE_NUMFMT_P_GSEA)
    return(RELEASE_NUMFMT_P)
  }
  if (column_name %in% RELEASE_THREE_DECIMALS) return("0.000")
  if (column_name %in% RELEASE_TWO_DECIMALS) return("0.00")
  finite <- values[is.finite(values)]
  if (!length(finite)) return(NA_character_)
  if (all(finite == round(finite)) && max(abs(finite)) < 1e15) return("0")
  "0.000"
}

# --------------------------------------------------------------------------------------
# display labels for controlled-vocabulary values
# --------------------------------------------------------------------------------------
# The .tsv files keep the controlled vocabulary (`mcherry`, `paired_cno`,
# `cfos_paired_cno_over_cfos_paired_veh`), which the tables join on. The workbook, the
# figure source data and the prose show these display labels; the comparison key stays in
# the workbook as its last column. Names are in figure order.

RELEASE_SAMPLE_CLASS_DISPLAY <- c(
  neuropil = "neuropil",
  cfos     = "c-Fos+",
  mcherry  = "mCherry+",
  neuron   = "NeuN+"
)

RELEASE_CONDITION_DISPLAY <- c(
  paired_veh   = "paired VEH",
  paired_cno   = "paired CNO",
  unpaired_veh = "unpaired VEH",
  unpaired_cno = "unpaired CNO"
)

RELEASE_TREATMENT_DISPLAY <- c(cno = "CNO", veh = "VEH")

# Contrast families in comparison order: Pairing first, then CNO in paired, then CNO in
# unpaired animals. RELEASE_FAMILY_LABEL is the workbook cell value; RELEASE_FAMILY_DISPLAY
# the short form used inside prose.
RELEASE_FAMILY_LABEL <- c(
  paired_veh_vs_unpaired_veh   = "Pairing (paired VEH vs unpaired VEH)",
  paired_cno_vs_paired_veh     = "CNO in paired animals",
  unpaired_cno_vs_unpaired_veh = "CNO in unpaired animals"
)
RELEASE_FAMILY_DISPLAY <- c(
  paired_veh_vs_unpaired_veh   = "Pairing",
  paired_cno_vs_paired_veh     = "CNO in paired animals",
  unpaired_cno_vs_unpaired_veh = "CNO in unpaired animals"
)

RELEASE_QUERY_LIST_DISPLAY <- c(
  fdr_significant_all                   = "All FDR < 0.05 proteins",
  fdr_significant_higher_in_numerator   = "FDR < 0.05, higher in group 1",
  fdr_significant_higher_in_denominator = "FDR < 0.05, lower in group 1",
  top_absolute_log2fc                   = "|effect size| > 1 SD (no FDR filter)"
)

RELEASE_ANALYSIS_ROLE_DISPLAY <- c(
  canonical              = "primary",
  alternative_query_list = "alternative",
  sensitivity            = "sensitivity"
)

RELEASE_RANK_STATISTIC_DISPLAY <- c(
  moderated_t    = "moderated t",
  t              = "moderated t",
  log2fc         = "effect size (SD units)",
  not_applicable = "none"
)

RELEASE_ONTOLOGY_DISPLAY <- c(GO_BP = "GO-BP", GO_CC = "GO-CC", KEGG = "KEGG")

# Level-1 cell types of the ewceData reference (level-2 names are the reference's cluster
# codes and are shown as they are).
RELEASE_CELL_TYPE_DISPLAY <- c(
  astrocytes_ependymal = "astrocytes/ependymal",
  endothelial_mural    = "endothelial/mural",
  interneurons         = "interneurons",
  microglia            = "microglia",
  oligodendrocytes     = "oligodendrocytes",
  pyramidal_CA1        = "CA1 pyramidal",
  pyramidal_SS         = "somatosensory pyramidal"
)

release_lookup_display <- function(x, map) {
  x <- as.character(x)
  out <- unname(map[x])
  ifelse(is.na(out), x, out)
}

release_sample_class_display <- function(x) release_lookup_display(x, RELEASE_SAMPLE_CLASS_DISPLAY)
release_condition_display <- function(x) release_lookup_display(x, RELEASE_CONDITION_DISPLAY)

#' Reader-facing label for one comparison, e.g. "mCherry+: paired VEH vs unpaired VEH".
release_comparison_label <- function(sample_class, numerator_condition,
                                     denominator_condition) {
  paste0(release_sample_class_display(sample_class), ": ",
         release_condition_display(numerator_condition), " vs ",
         release_condition_display(denominator_condition))
}

#' Sort rank of comparisons: contrast family first (Pairing, CNO in paired, CNO in
#' unpaired), then compartment in figure order. NA for an unknown family or class.
release_comparison_rank <- function(sample_class, contrast_family) {
  f <- match(as.character(contrast_family), names(RELEASE_FAMILY_LABEL))
  s <- match(as.character(sample_class), names(RELEASE_SAMPLE_CLASS_DISPLAY))
  (f - 1L) * length(RELEASE_SAMPLE_CLASS_DISPLAY) + s
}

#' Display direction: "higher in paired VEH" for higher_in_paired_veh,
#' enriched_toward_paired_veh, or up (numerator) / down (denominator).
release_direction_display <- function(direction, numerator_condition, denominator_condition) {
  d <- as.character(direction)
  out <- d
  pre <- grepl("^(higher_in|enriched_toward)_", d)
  out[pre] <- paste("higher in",
                    release_condition_display(sub("^(higher_in|enriched_toward)_", "", d[pre])))
  up <- d %in% "up"
  out[up] <- paste("higher in", release_condition_display(numerator_condition[up]))
  down <- d %in% "down"
  out[down] <- paste("higher in", release_condition_display(denominator_condition[down]))
  out[d %in% "both"] <- "either"
  out[d %in% "abundant"] <- "most abundant"
  out
}

#' Reader-facing list of the GSEA blocks with gene sets fgsea could not test.
#'
#' `untested` is release_fgsea_untested(); `contrasts` the primary contrast manifest. Gives
#' e.g. "Pairing: neuropil 31 GO-BP, c-Fos+ 11 GO-BP; CNO in paired animals: c-Fos+ 114
#' GO-BP", or "" when no block matches.
release_untested_blocks_text <- function(untested, contrasts, families,
                                         ontologies = "GO_BP", classes = NULL) {
  u <- merge(untested, contrasts[c("canonical_comparison", "sample_class", "contrast_family")],
             by = "canonical_comparison")
  keep <- u$contrast_family %in% families & u$ontology %in% ontologies
  if (!is.null(classes)) keep <- keep & u$sample_class %in% classes
  u <- u[keep, , drop = FALSE]
  if (!nrow(u)) return("")
  class_order <- names(RELEASE_SAMPLE_CLASS_DISPLAY)
  parts <- vapply(families[families %in% u$contrast_family], function(f) {
    uf <- u[u$contrast_family == f, , drop = FALSE]
    uf <- uf[order(match(uf$sample_class, class_order), uf$ontology), , drop = FALSE]
    paste0(RELEASE_FAMILY_DISPLAY[[f]], ": ",
           paste0(release_sample_class_display(uf$sample_class), " ", uf$n_sets, " ",
                  sub("_", "-", uf$ontology), collapse = ", "))
  }, character(1))
  paste(parts, collapse = "; ")
}

# --------------------------------------------------------------------------------------
# reader text for the constructed tables (2o, 2p)
# --------------------------------------------------------------------------------------

# Words that do not belong in reader-facing workbook text (captions, README, tables 2o and
# 2p): internal terms, revision history, internal paths. Stage 07 stops on a hit.
RELEASE_READER_FORBIDDEN <- c("canonical", "contract", "audit", "lineage", "superseded",
                              "original submission", "reviewer", "editor", "changelog",
                              "Figure3", "Supplementary_proteomics", "03_output",
                              "99_historical", "S:/", "D:\\", "Fabian", "sha256",
                              "SHA-256", "pooled", "legacy")

# Words that stay in capitals; any other all-capital word of three or more letters is
# lowered (the source notes use capitals for emphasis).
RELEASE_KEEP_CAPS <- c("GSEA", "KEGG", "EWCE", "FDR", "NES", "ORA", "PCA", "CNO", "VEH",
                       "LCM", "CTD", "TRUE", "FALSE", "QC", "BH", "SD", "LC", "MS", "UMAP",
                       "GCT", "DOSE", "RNA", "DNA", "ATP", "GABA", "AND")

#' Translate machine keys in free text into display labels: comparison keys, contrast
#' families, conditions, compartments, column names in row filters, and emphasis capitals.
release_reader_text <- function(x, contrasts) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  keys <- as.character(contrasts$canonical_comparison)
  labels <- release_comparison_label(contrasts$sample_class, contrasts$numerator_condition,
                                     contrasts$denominator_condition)
  for (i in order(-nchar(keys))) x <- gsub(keys[[i]], labels[[i]], x, fixed = TRUE)
  for (f in names(RELEASE_FAMILY_LABEL)) x <- gsub(f, RELEASE_FAMILY_LABEL[[f]], x, fixed = TRUE)
  rules <- c(
    "canonical_comparison == " = "comparison ",
    "canonical_comparison in " = "comparisons ",
    "contrast_family == " = "",
    "contrast_family in " = "",
    "analysis_role == canonical" = "primary analysis",
    "ontology == " = "",
    "ontology in " = "",
    "id_mapping_status == mapped" = "annotated proteins",
    "significant_fdr_0_05 == TRUE in " = "FDR < 0.05 in ",
    "adjusted_p_value" = "adjusted P",
    "ewce_analysis_type == Baseline" = "baseline analysis",
    "is_primary_setting == TRUE" = "primary setting",
    "sample_class in " = "compartments ",
    "condition == " = "condition ",
    "colored by sample_metadata " = "colored by ",
    "collection_plate" = "collection plate",
    "sample_class x condition" = "compartment x condition",
    "per-sample-class" = "per-compartment",
    "sample-class correction" = "compartment-label correction",
    "sample-class" = "compartment",
    "sample_class" = "compartment",
    "sample classes" = "compartments",
    "sample class" = "compartment",
    "AnimalID" = "animal",
    "moderated_t" = "moderated t",
    "log2fc" = "effect size (SD units)",
    "GO_BP" = "GO-BP",
    "GO_CC" = "GO-CC",
    " AND " = "; ",
    "_vs_" = " vs ",
    "mapped proteins" = "annotated proteins",
    " mapped)" = " annotated)",
    " up / " = " higher / ",
    " down)" = " lower)",
    # upstream run records use British spelling; the paper uses US spelling
    "coloured by sample_metadata " = "colored by ",
    "standardisation" = "standardization",
    "standardised" = "standardized",
    "normalisation" = "normalization",
    "normalised" = "normalized",
    "coloured" = "colored",
    "centred" = "centered",
    "analysed" = "analyzed",
    "programme" = "program",
    "labelled" = "labeled",
    "colour" = "color",
    "artefact" = "artifact",
    "judgement" = "judgment",
    "itemised" = "itemized",
    "both n.s." = "neither significant",
    "n.s." = "not significant")
  for (r in names(rules)) x <- gsub(r, rules[[r]], x, fixed = TRUE)
  # notation: "FDR<0.05" -> "FDR < 0.05", "3.2e-06" -> "3.2E-06"
  x <- gsub("FDR ?([<>]=?) ?([0-9])", "FDR \\1 \\2", x, perl = TRUE)
  x <- gsub("(?<![A-Za-z0-9.])([0-9]+(?:[.][0-9]+)?)e([-+][0-9]+)", "\\1E\\2", x, perl = TRUE)
  # {a, b} -> a and b
  x <- gsub("\\{([^},]*), ([^}]*)\\}", "\\1 and \\2", x)
  for (k in names(RELEASE_CONDITION_DISPLAY)) {
    x <- gsub(paste0("\\b", k, "\\b"), RELEASE_CONDITION_DISPLAY[[k]], x, perl = TRUE)
  }
  class_rules <- c("mcherry" = "mCherry+", "mCherry" = "mCherry+", "cfos" = "c-Fos+",
                   "cFos" = "c-Fos+", "neuron" = "NeuN+", "Neuron" = "NeuN+",
                   "Neuropil" = "neuropil")
  for (k in names(class_rules)) {
    x <- gsub(paste0("\\b", k, "\\b(?![-+])"), class_rules[[k]], x, perl = TRUE)
  }
  x <- gsub("\\bcanonical\\b", "primary", x, perl = TRUE)
  # emphasis capitals; a lowered word that opens the text keeps its capital initial
  vapply(x, function(s) {
    m <- gregexpr("\\b[A-Z]{3,}\\b", s, perl = TRUE)
    words <- regmatches(s, m)[[1]]
    if (!length(words)) return(s)
    lowered <- !words %in% RELEASE_KEEP_CAPS
    new <- ifelse(lowered, tolower(words), words)
    if (m[[1]][[1]] == 1L && lowered[[1]]) {
      new[[1]] <- paste0(substring(words[[1]], 1L, 1L), substring(new[[1]], 2L))
    }
    regmatches(s, m) <- list(new)
    s
  }, character(1), USE.NAMES = FALSE)
}

#' Model strings of the secondary-analysis manifest in words. Applied after
#' release_reader_text(); the code syntax of the manifest (eBayes(...), full rank, df.total)
#' is written out.
release_model_words <- function(x) {
  x <- as.character(x)
  x <- gsub("per compartment: limma ~0 + condition (", "limma ~0 + condition per compartment (",
            x, fixed = TRUE)
  x <- gsub("; eBayes(robust = TRUE, trend = FALSE); BH within contrast",
            ", robust empirical Bayes without trend, BH within comparison", x, fixed = TRUE)
  x <- gsub("eBayes(robust = TRUE, trend = FALSE)", "robust empirical Bayes without trend", x,
            fixed = TRUE)
  x <- gsub("BH within contrast", "BH within comparison", x, fixed = TRUE)
  # "limma ~ 0 + a + b; full rank 6/6; residual df 6; eBayes df.total 27.435"
  re <- "limma ~ ?0 [+] ([^;]+); (full )?rank [0-9]+/[0-9]+; residual df ([0-9]+); eBayes df[.]total ([0-9.]+)"
  vapply(x, function(s) {
    m <- regmatches(s, regexec(re, s))[[1]]
    if (!length(m)) return(s)
    sub(m[[1]], paste0("limma ~0 + ", m[[2]], " (residual df ", m[[4]], ", ",
                       sprintf("%.1f", as.numeric(m[[5]])),
                       " after empirical Bayes moderation)"), s, fixed = TRUE)
  }, character(1), USE.NAMES = FALSE)
}

#' snake_case tokens left in reader text, other than file names.
release_residual_machine_tokens <- function(x) {
  x <- as.character(x)
  x <- gsub("[A-Za-z0-9_./-]+[.](csv|tsv|gz|xlsx|md)\\b", "", x)
  hits <- unlist(regmatches(x, gregexpr("\\b[A-Za-z0-9]+_[A-Za-z0-9_]+\\b", x, perl = TRUE)))
  unique(hits)
}

#' Split text into sentences, keeping "Fig. 2" and "Supp. Fig. 4" together.
release_split_sentences <- function(x) {
  for (abbr in c("Supp. Fig.", "Fig.", "e.g.", "i.e.", "n.s.", "vs.")) {
    x <- gsub(abbr, gsub(".", "#", abbr, fixed = TRUE), x, fixed = TRUE)
  }
  parts <- strsplit(x, "(?<=[.])\\s+(?=[A-Za-z(])", perl = TRUE)[[1]]
  parts <- gsub("#", ".", parts, fixed = TRUE)
  parts[nzchar(trimws(parts))]
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

# Fixed widths for free-text columns of the small sheets, which wrap.
RELEASE_WIDTH_OVERRIDES <- list(
  README = c(item = 30, statement = 110),
  Column_Definitions = c(table = 8, column = 36, definition = 90, units = 14),
  Comparisons = c(comparison_label = 40, contrast_family = 36),
  Enrichment_Settings = c(gene_set_source = 48, note = 40, permutations = 30),
  Secondary_Analyses = c(analysis_name = 34, purpose = 40, model = 60, n_per_group = 22,
                         results_location = 40, key_result = 70),
  Figure_Panels = c(panel_label = 16, panel_title = 40, statistical_unit = 36,
                    panel_source = 40, rows_used = 50, key_numbers = 50, note = 70),
  Software = c(component = 40, version = 22, used_for = 50)
)
