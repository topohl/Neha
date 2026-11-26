#!/usr/bin/env Rscript

# Publication release, stage 11 -- README_DATA.md and the data dictionary.
#
# Produces
#   README_DATA.md                 the human entry point to the package
#   metadata/data_dictionary.tsv   every column of every release table, described
#
# The dictionary walks the tables in TABLES and looks up each column in
# TABLE_COLUMN_OVERRIDES, DESCRIPTIONS, WIDE_NOTE (the two matrices) or, for the SDRF,
# pride/sdrf_field_status.tsv; a column with no description stops the build.
#
# README_DATA.md must contain the collection-plate sentence verbatim and none of
# FORBIDDEN_PLATE_WORDING (both defined in R/release_validation.R).

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
release_require("digest")

DATA_ROOT <- release_data_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/11_build_readme_and_dictionary.R"

release_banner("stage 11 -- README_DATA and data dictionary")

# Collection-plate wording, shared with the changelog (R/release_validation.R).
PLATE_SENTENCE <- RELEASE_PLATE_SENTENCE
NO_BATCH_SENTENCE <- RELEASE_PLATE_NO_EFFECT_SENTENCE
FORBIDDEN_WORDING <- FORBIDDEN_PLATE_WORDING

release_check_forbidden_wording <- function(lines) {
  flat <- tolower(paste(lines, collapse = " "))
  hits <- FORBIDDEN_WORDING[vapply(FORBIDDEN_WORDING,
                                   function(p) grepl(tolower(p), flat, fixed = TRUE),
                                   logical(1))]
  if (length(hits)) {
    stop("Forbidden collection-plate wording in generated prose: ",
         paste(hits, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

rel <- function(...) release_path(..., create_dir = FALSE)
read_release <- function(path) {
  if (grepl("[.]gz$", path)) {
    read.delim(gzfile(path), sep = "\t", quote = "", stringsAsFactors = FALSE,
               check.names = FALSE, nrows = 5L)
  } else {
    read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE,
               nrows = 5L)
  }
}

diff_summary <- read.delim(rel("differential_analysis", "primary_differential_summary.tsv"),
                           sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
corrections <- read.delim(rel("metadata", "sample_class_corrections.tsv"), sep = "\t",
                          quote = "", stringsAsFactors = FALSE, check.names = FALSE)
pride_readiness <- read.delim(rel("pride", "pride_readiness.tsv"), sep = "\t", quote = "",
                              stringsAsFactors = FALSE, check.names = FALSE)
pride_status <- pride_readiness$value[pride_readiness$key == "pride_status"]
n_required_missing <- pride_readiness$value[pride_readiness$key == "n_sdrf_required_missing"]
experimenter_metadata <- read.delim(rel("metadata", "experimenter_metadata.tsv"), sep = "\t",
                                    quote = "", stringsAsFactors = FALSE,
                                    check.names = FALSE, na.strings = character(0))
sample_preparation_protocol <- read.delim(
  rel("metadata", "sample_preparation_protocol.tsv"), sep = "\t", quote = "",
  stringsAsFactors = FALSE, check.names = FALSE, na.strings = character(0))
mv <- function(id, standardized = FALSE) {
  release_metadata_value(experimenter_metadata, id, standardized = standardized)
}
# Full reads (read_release above returns 5 rows, enough for the column walk only).
read_full <- function(path) {
  if (grepl("[.]gz$", path)) {
    read.delim(gzfile(path), sep = "\t", quote = "", stringsAsFactors = FALSE,
               check.names = FALSE)
  } else {
    read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  }
}
stat_run <- release_read_stat_run(DATA_ROOT)
ave <- release_aveexpr_summary(
  read_full(rel("differential_analysis", "primary_differential_proteins.tsv.gz")))
gsea_go_all <- read_full(rel("enrichment", "primary_GSEA_GO_BP.tsv.gz"))
gsea_cc_all <- read_full(rel("enrichment", "primary_GSEA_GO_CC.tsv.gz"))
ora_go_all <- read_full(rel("enrichment", "primary_ORA_GO_BP.tsv.gz"))
# Measured on the released measurement-level matrix, as stage 03 measures it.
std_meas <- local({
  m <- read_full(rel("processed_data", "protein_abundance_measurement_level.tsv.gz"))
  release_standardization_evidence(as.matrix(m[, setdiff(names(m), "protein_group_id")]))
})
# Per-protein SDs of the released animal-level matrix (the mean of the two hemispheres).
animal_row_sd <- local({
  m <- read_full(rel("processed_data", "protein_abundance_animal_level.tsv.gz"))
  apply(as.matrix(m[, setdiff(names(m), "protein_group_id")]), 1L, stats::sd)
})

# Leave-one-animal-out check of the mCherry learning comparison, written by the inferential
# branch (04_differential_expression_enrichment/06_mcherry_learning_robustness.r).
LOO_REL <- c(RELEASE_STATISTICS_RUN$inferential_relative, "robustness",
             "mcherry_learning_leave_one_out.csv")
loo <- release_read_csv(release_rel_path(LOO_REL, DATA_ROOT))
loo_drop <- loo[loo$animal_dropped != "none", , drop = FALSE]
loo_full <- loo[loo$animal_dropped == "none", , drop = FALSE]
if (nrow(loo_full) != 1L || !nrow(loo_drop)) {
  stop("Leave-one-animal-out table lacks the full fit or the refits.", call. = FALSE)
}

# Gene sets fgsea could not estimate (unbalanced positive and negative statistics), as the
# enrichment runs recorded them (GO: audits/<comparison>/gsea_untested_sets.csv; KEGG:
# warnings.csv). clusterProfiler drops such sets from its output, so they are absent from
# the GSEA tables, not non-significant.
fgsea_na <- release_fgsea_untested(unique(diff_summary$canonical_comparison), DATA_ROOT)
fgsea_na_mch <- fgsea_na[fgsea_na$canonical_comparison ==
                           "mcherry_paired_veh_over_mcherry_unpaired_veh", , drop = FALSE]
permutations <- release_enrichment_permutations(DATA_ROOT)
contrast_manifest <- read_full(rel("metadata", "primary_contrast_manifest.tsv"))
untested_blocks <- function(families, ontologies = "GO_BP", classes = NULL) {
  x <- release_untested_blocks_text(fgsea_na, contrast_manifest, families, ontologies, classes)
  if (nzchar(x)) paste(x, "sets") else "none"
}
# The sample classes the released source data of Supp. Fig. 4 F show (all four since
# 2026-11-25).
figure_map_full <- read_full(rel("editor_source_data", "figure_source_map.tsv"))
sf_row <- figure_map_full[figure_map_full$figure == "Supplementary_proteomics" &
                            figure_map_full$panel == "F", , drop = FALSE]
sf_classes <- unique(as.character(release_read_csv(
  strsplit(sf_row$panel_source_data_canonical_path, ";", fixed = TRUE)[[1]][[1]])$sample_class))
ora_cc_all <- read_full(rel("enrichment", "primary_ORA_GO_CC.tsv.gz"))
ewce_model <- release_ewce_differential_model(DATA_ROOT)

class_metadata <- experimenter_metadata[
  experimenter_metadata$category == "sample_class_definition", , drop = FALSE]
class_metadata <- class_metadata[order(as.integer(class_metadata$sort_order)), , drop = FALSE]
# README table of sample classes: data value (`mcherry`), display name ("mCherry"),
# historical label and definition.
class_markdown <- vapply(seq_len(nrow(class_metadata)), function(i) {
  historical <- switch(class_metadata$sample_class[i], mcherry = "mCherryN",
                       neuropil = "bg / Background", cfos = "cFosN", neuron = "neuron")
  paste0("| `", class_metadata$sample_class[i], "` | ",
         release_sample_class_display(class_metadata$sample_class[i]), " | `", historical,
         "` | ", class_metadata$value[i], " |")
}, character(1))

# --------------------------------------------------------------------------------------
# column descriptions
# --------------------------------------------------------------------------------------

DESCRIPTIONS <- c(
  # identity / design
  sample_id = "Canonical measurement identifier. This string IS the acquisition run name and is also the column name in the measurement-level abundance matrix.",
  acquisition_run_name = "Same value as sample_id, named so its origin is clear.",
  raw_file_basename = "Original acquisition file name, `<run>.d`. Proven: basename of the sample_annotation.xlsx Name column with the .d suffix removed equals sample_id for all 96 runs.",
  raw_file_original_path = "Verbatim acquisition-machine path recorded at search time, e.g. D:\\Proteomics\\Fabian\\Tobias\\<run>.d.",
  raw_file_extension = "Acquisition file extension. `.d` is a vendor-specific acquisition directory format. It identifies a vendor family, not an instrument model.",
  AnimalID = "Animal identifier. 12 animals: C11 C12 C14 C25 C26 C27 C33f C34f C45f C46 C47 C510.",
  hemisphere = "Left or Right. The anatomical side from which the marker-defined LCM material was obtained.",
  hemisphere_evidence_source = "How hemisphere was established. `explicit_evidence` = the _left/_right suffix on the annotation label, cross-checked against the numeric ReplicateGroup.",
  historical_hemisphere_label = "Legacy annotation label carrying the explicit side, e.g. cFosN_left.",
  legacy_replicate_group = "Legacy numeric hemisphere encoding: 1 = Left, 2 = Right. Retained for traceability.",
  sample_class = "Marker-defined LCM class: mcherry, neuropil, cfos or neuron. THIS IS THE ANALYSIS-TIME ASSIGNMENT -- the class the locked GCTs and every validated statistic were computed on. For 6 of the 96 measurements it differs from the pre-correction assignment, which is retained in original_sample_class; see metadata/sample_class_corrections.tsv.",
  sample_class_historical_alias = "Historical filename alias. `neuropil` was written `bg` (from background) in older filenames; it is a tissue compartment, not a negative control or a blank.",
  sample_class_historical_group_label = "Legacy group label from the annotation workbook, verbatim: cFosN, Background, mCherryN, Neuron.",
  original_sample_class = "PRE-CORRECTION sample class, normalised to the canonical vocabulary from the sample_annotation.xlsx group_label. Independently reproduced for all 96 measurements by the autosampler plate layout, and for the 6 corrected measurements by the retained UMAP correction table. Retained so the original identity is never lost.",
  analysis_sample_class = "ANALYSIS-TIME sample class. Identical to sample_class; named so the pairing with original_sample_class is clear.",
  sample_class_plate_layout_implied = "Sample class implied by the autosampler well row letter (A/E cfos, B/F neuropil, C/G mcherry, D/H neuron). This is a PRE-CORRECTION record and equals original_sample_class for all 96 measurements.",
  sample_class_corrected = "TRUE for the 6 measurements whose sample-class metadata were corrected during historical sample-identity QC, FALSE for the other 90. TRUE does NOT mean the sample was mislabelled or that any measurement changed; only the class metadata were reassigned.",
  sample_class_correction_status = "confirmed_intentional_switch_correction for the 6 corrected measurements, not_corrected for the other 90.",
  sample_class_correction_method = "How the correction was established, for corrected rows: historical_sample_identity_qc_umap_correction_table. Blank where nothing was corrected.",
  sample_class_correction_provenance = "Project-relative path of the preserved historical correction record, for corrected rows. Blank where nothing was corrected.",
  sample_class_correction_reference_sha256 = "SHA256 of the preserved correction record, recomputed and verified at build time against the value the audit recorded. Blank where nothing was corrected.",
  sample_class_correction_date = "Date of the applied metadata edit, for corrected rows. Blank where nothing was corrected.",
  plate_position = "Autosampler position as plate-set token and well, e.g. S5-B6.",
  correction_status = "Verdict on the sample-class correction: confirmed_intentional_switch_correction. A deliberate correction applied during historical sample-identity QC, not an unresolved discrepancy.",
  correction_method = "How the correction was established: historical_sample_identity_qc_umap_correction_table -- identified from UMAP nearest-class-centre distances, recorded in a retained project correction table, and subsequently supported by independent proteomic profile-similarity and left/right-pair analyses.",
  correction_date = "Date of the applied metadata edit.",
  correction_reference = "Project-relative path of the preserved historical correction record. This table is the preserved record; no surviving prose note explains the rationale.",
  correction_reference_sha256 = "SHA256 of the correction record, recomputed at build time and required to match the value the audit recorded.",
  correction_reference_mtime = "Modification time of the correction record. Reported for context; the SHA256, not the mtime, is what the build verifies.",
  quantitative_values_changed = "Always FALSE. The quantitative matrix was compared across the metadata edit and the maximum absolute numeric difference was 0. Only sample-class metadata were reassigned.",
  umap_nearest_class_suggestion = "The class the UMAP correction table recorded as the nearest class centre. Published so the divergence from the applied correction is visible: it agrees for 5 of the 6, and differs for N60 (nearest centre cfos, applied neuron). UMAP identified the switch; it is not treated as ground truth for what the switch was.",
  original_class_corroborating_records = "The independent pre-correction records that agree on original_sample_class for this measurement.",
  animal_level_unit_affected = "The animal-level unit this measurement was averaged into. These units were built from the corrected assignments in the first place, so the validated results already reflect them.",
  condition_code = "Numeric condition code: 1 = paired_cno, 2 = paired_veh, 3 = unpaired_cno, 4 = unpaired_veh.",
  condition = "Condition label: paired_cno, paired_veh, unpaired_cno or unpaired_veh.",
  pairing_status = "paired or unpaired, decomposed from condition.",
  treatment = "cno or veh, decomposed from condition. Tokens are the project's literal condition coding.",
  collection_plate = "Plate1 or Plate2: which plate the sample was collected onto. NOT a proteomics preparation, digestion, LC-MS, acquisition or instrument batch.",
  plate_sample_number = "Sample number token from the acquisition run name, N01..N96.",
  injection_index = "Numeric running index 1..96 parsed from the plate sample number.",
  plate_set_token = "Set token from the acquisition run name (S4 or S5).",
  well_position = "Well position token from the acquisition run name, e.g. A1, H7.",
  acquisition_run_trailing_id = "Trailing numeric identifier in the acquisition run name, e.g. 8782.",
  instrument_alias_token = "Local instrument alias from the acquisition run name (`Olive`). NOT an instrument model.",
  instrument_model = "Instrument model. EMPTY: not recorded anywhere in the project. Not inferred.",
  acquisition_date = "Acquisition date in ISO form, derived from the run-name date token. All 96 runs carry 2024-12-17.",
  acquisition_date_token = "Verbatim date token from the acquisition run name.",
  lc_and_method_token = "Verbatim LC/method token from the acquisition run name. Preserved undecoded; it is an operator label, not a documented acquisition method.",
  analysis_unit_primary = "The unit of inference for this dataset: always `animal`.",
  included_in_animal_level = "TRUE if this measurement contributed to an animal-level unit. All 96 did.",
  animal_level_sample_id = "The animal-level column this measurement was averaged into, e.g. C11_cfos.",
  measurement_level_matrix_column = "Column name of this measurement in the measurement-level abundance matrix.",
  source_assignment_status = "Assignment outcome recorded by the aggregation stage, e.g. included_in_primary.",
  legacy_sample_number = "Legacy sampleNumber field from sample_info.xlsx.",
  legacy_shortname = "Legacy shortname field from sample_info.xlsx.",
  legacy_group2_hemisphere_coded = "Legacy `group2` field. WARNING: its numeric suffix is the ReplicateGroup (hemisphere), NOT the condition code. `celltype_ExpGroup` was the condition-coded field.",
  # animal-level
  n_hemisphere_measurements = "Number of hemisphere measurements averaged into this animal-level unit. Always 2.",
  hemispheres_present = "Semicolon-separated hemispheres contributing to this unit. Always Left;Right.",
  source_sample_ids = "Semicolon-separated acquisition run names averaged into this unit.",
  source_raw_file_basenames = "Semicolon-separated original `.d` acquisition file names for this unit.",
  left_source_sample = "Acquisition run name of the left-hemisphere measurement.",
  right_source_sample = "Acquisition run name of the right-hemisphere measurement.",
  animal_level_column_name = "Column name of this unit in the animal-level abundance matrix.",
  phenotype_within_unit = "sample_class and condition concatenated, e.g. cfos_paired_cno. The statistics model condition separately within each sample class.",
  aggregation_policy = "How the two hemispheres were combined: equal-weight mean of the Left and Right values of the measurement-level matrix, with no further transformation. The value is the policy name recorded when the animal-level matrix was built; the values it averaged are standardised per protein, not log2 abundance.",
  analysis_unit = "The unit of inference for the row: `animal`.",
  # metadata field provenance
  field = "Metadata field being described.",
  status = "Field status: KNOWN_VERIFIED, KNOWN_BUT_NEEDS_STANDARDIZATION, MISSING_RECOVERABLE, MISSING_UNKNOWN or NOT_APPLICABLE.",
  current_source = "File that currently holds the field, or NONE. Files in the code repository are named relative to its proteomic_analysis/ folder.",
  evidence = "How the status was established.",
  note = "Additional caveat or clarification.",
  # contrasts
  contrast_family = "One of paired_cno_vs_paired_veh, paired_veh_vs_unpaired_veh, unpaired_cno_vs_unpaired_veh.",
  numerator_condition = "Condition in the numerator. A positive effect_size_sd_units means higher in this condition.",
  denominator_condition = "Condition in the denominator.",
  numerator_code = "Condition code of the numerator arm.",
  denominator_code = "Condition code of the denominator arm.",
  numerator_phenotype = "sample_class and numerator condition concatenated.",
  denominator_phenotype = "sample_class and denominator condition concatenated.",
  canonical_comparison = "Canonical comparison key, `<numerator_phenotype>_over_<denominator_phenotype>`. The join key across differential, enrichment and EWCE tables.",
  canonical_contrast = "Canonical contrast name, `<numerator_phenotype>_vs_<denominator_phenotype>`. Used for the canonical file names.",
  historical_comparison_alias = "Historical short name for the comparison, e.g. bg2_bg4 for neuropil paired_veh vs unpaired_veh.",
  n_numerator_animals = "Number of independent animals in the numerator arm. Always 3.",
  n_denominator_animals = "Number of independent animals in the denominator arm. Always 3.",
  primary_or_secondary = "Whether the row belongs to the 12 primary contrasts or to a secondary analysis.",
  statistical_model = paste("The model actually fitted for this comparison:",
                            paste0(RELEASE_PRIMARY_MODEL$label, "."),
                            "Contrast = numerator minus denominator. In the secondary manifest, the model of that analysis."),
  residual_df = "Residual degrees of freedom of the fitted model (12 animals minus 4 condition coefficients), before empirical-Bayes moderation.",
  statistics_software = "Software and script that computed the statistics, with the versions the statistics run manifest recorded at run time.",
  gsea_ranking_statistic = "Statistic the canonical GSEA ranked by: the moderated t.",
  interpretation = "What the comparison means, including any identifiability caveat.",
  differential_statistics_file = "Release file holding the per-protein statistics for this comparison.",
  canonical_split_forward_path = "Canonical differential table this comparison was read from.",
  canonical_mapped_forward_path = "Canonical identifier-mapped table for this comparison.",
  canonical_enrichment_dir = "Canonical GO-BP enrichment output folder for this comparison (the GO_BP folder of the GO run).",
  canonical_enrichment_dir_go_cc = "Canonical GO-CC enrichment output folder for this comparison (the GO_CC folder of the GO run).",
  canonical_enrichment_dir_kegg = "Canonical KEGG enrichment output folder for this comparison (the GO_BP folder of the KEGG run, whose GO tables are superseded).",
  # secondary manifest
  analysis_id = "Stable identifier for the secondary or excluded analysis.",
  analysis_name = "Human-readable name of the analysis.",
  analysis_type = "What kind of analysis it is.",
  comparisons = "Comparisons the analysis covers.",
  condition_scope = "Conditions the analysis is restricted to.",
  ranking_statistic = "Statistic used for ranking, where applicable.",
  n_per_group = "Sample size per arm.",
  included_in_primary_contrast_manifest = "Always FALSE here; these rows are excluded from the primary table.",
  canonical_source = "Canonical output location for this analysis.",
  generating_script = "Script that produced the analysis, where one exists.",
  caveat = "The limitation a reader must carry when reading this analysis.",
  # features
  protein_group_id = "Protein-group identifier: UniProt ENTRY NAME(s), semicolon-separated for multi-protein groups, e.g. MFA1A_MOUSE;MFA1B_MOUSE. This is the row identifier of every abundance and statistics table.",
  protein_group_uniprot_accessions = "All UniProt accessions in the protein group, semicolon-separated, from the search output.",
  uniprot_accession = "UniProt accession of the leading (master) protein in the group.",
  gene_symbol = "Resolved gene symbol for the leading protein. Never a UniProt identifier.",
  gene_symbols_in_group = "All gene symbols in the protein group, semicolon-separated, as reported by the search output.",
  protein_description = "Protein description from the First.Protein.Description annotation field in the retained search/quantification output.",
  n_proteins_in_group = "Number of proteins in the group.",
  is_protein_group = "TRUE if the group contains more than one protein.",
  id_mapping_status = "mapped or unmapped. 5327 of 5349 map; the 22 unmapped are non-mouse contaminant identifiers.",
  id_mapping_strategy = "Which mapping strategy resolved the identifier.",
  id_mapping_multi_protein = "Whether the mapping had to choose a leading protein from a group.",
  source_matrix_row_index = "Row index in the canonical differential table, for traceability.",
  source_matrix_description_field = "Verbatim value of the source matrix `Description` field. NOTE: despite its name this field holds GENE SYMBOLS, not a protein description. Retained only for traceability.",
  id_mapping_reference_file = "Reference database used for identifier mapping.",
  id_mapping_reference_version = "Version encoded in the reference file, if any.",
  id_mapping_reference_snapshot_date_utc = "Date the reference database snapshot was taken.",
  id_mapping_reference_sha256 = "SHA256 of the reference database file.",
  # differential
  source_matrix_gene_symbols = "Verbatim `Description` value from the canonical differential table; holds gene symbols, not a description.",
  effect_size_sd_units = "PUBLIC EFFECT SIZE. Difference of group means (numerator minus denominator) on the protein-standardised abundance scale -- a standardized abundance difference expressed in SD units of that scale. NOT a log2 fold change, despite being stored as logFC in the statistical GCT (limma's coefficient name) and log2fc in the canonical split tables. Positive means higher in the numerator condition. A value of 2 means two standard deviations on the standardised scale, not a four-fold change. Deliberately not called a standardised mean difference unqualified, because that phrasing reads as Cohen's d (a pooled within-group SD) rather than the across-dataset per-protein SD actually used. The exported value is bit-identical to the source logFC; no transformation, rescaling or recomputation was applied.",
  effect_size_units = "Units of effect_size_sd_units: per-protein standard deviation of the standardised abundance scale.",
  effect_size_definition = "Full statement of what the effect size is and is not. Carried once per comparison in the summary table rather than repeated on every protein row.",
  source_statistic_field = "Bare name of the canonical source statistic the effect size was copied from: logFC. Machine-readable provenance; the exported effect_size_sd_units equals this source field exactly.",
  effect_size_source_column = "The canonical column the effect size was copied from, named so the published value can be traced: logFC (statistical GCT, limma) / log2fc (canonical split tables).",
  average_standardized_abundance = "limma AveExpr: the protein's mean standardised abundance over the 12 animal-level units of this comparison's sample class (the columns of that class's model), so it is the same for the three comparisons of a sample class. It is not centred at 0 within a sample class; averaged over the four sample classes it is close to 0 for every protein (bound checked by 13_validate_release.R), as a per-protein standardised input implies.",
  moderated_t = "Moderated t statistic. The canonical GSEA ranking statistic.",
  P.Value = "Raw p-value from the moderated t test.",
  adj.P.Val = "Benjamini-Hochberg adjusted p-value.",
  B_log_odds = "Log-odds of differential abundance (limma B statistic).",
  neg_log10_P_value = "-log10 of the raw p-value, as recorded by the canonical analysis.",
  signed_neg_log10_P_value = "-log10 p-value signed by the direction of change. Rounded to 4 decimals in the canonical source.",
  significant_fdr_0_05 = "Canonical significance flag, adj.P.Val < 0.05. Copied, not recomputed.",
  n_proteins_tested = "Protein groups tested in this comparison (all 5349).",
  n_significant_fdr_0_05 = "Protein groups at FDR < 0.05 over all tested rows.",
  n_significant_higher_in_numerator = "Significant protein groups with a positive effect_size_sd_units, i.e. higher in the numerator condition.",
  n_significant_higher_in_denominator = "Significant protein groups with a negative effect_size_sd_units, i.e. higher in the denominator condition.",
  n_proteins_tested_mapped_only = "Protein groups tested that carry a mouse identifier (5327).",
  n_significant_fdr_0_05_mapped_only = "Significant protein groups restricted to mapped identifiers. This is the count the manuscript figures report.",
  n_significant_higher_in_numerator_mapped_only = "Significant, positive effect_size_sd_units, mapped identifiers only.",
  n_significant_higher_in_denominator_mapped_only = "Significant, negative effect_size_sd_units, mapped identifiers only.",
  significance_definition = "How the canonical significance flag is defined.",
  canonical_source_path = "Canonical file the statistics were read from.",
  canonical_source_sha256 = "SHA256 of that canonical file.",
  # enrichment
  analysis = "Analysis type: GSEA, ORA or EWCE.",
  ontology = "Gene-set collection: GO_BP, GO_CC or KEGG.",
  query_list = "For ORA, which query list the terms were tested against.",
  term_id = "Gene-set identifier (GO or KEGG).",
  term_name = "Gene-set name.",
  set_size = "Number of genes from the ranked list in the gene set (GSEA only).",
  enrichment_score = "GSEA enrichment score.",
  NES = "GSEA normalised enrichment score. Positive means enriched toward the numerator condition.",
  p_value = "Raw p-value.",
  adjusted_p_value = "Benjamini-Hochberg adjusted p-value. The canonical run used pvalue_cutoff = 1, so these tables hold ALL tested terms; filter here.",
  q_value = "q-value as reported by clusterProfiler.",
  gene_count = "Number of query genes annotated to the term (ORA only).",
  gene_ratio = "Query genes in term over query genes tested (ORA only).",
  background_ratio = "Background genes in term over background size (ORA only).",
  fold_enrichment = "Fold enrichment (ORA only).",
  rich_factor = "Rich factor (ORA only).",
  z_score = "z-score (ORA only).",
  rank_at_max = "Rank position at which the GSEA enrichment score is maximal.",
  leading_edge = "GSEA leading-edge summary string (tags/list/signal).",
  core_enrichment = "Gene members. For GSEA the leading-edge subset; for ORA the query genes annotated to the term. See core_enrichment_definition.",
  core_enrichment_definition = "States precisely what core_enrichment holds for this row type.",
  n_core_enrichment_genes = "Count of members in core_enrichment.",
  direction = "Direction of the result, expressed relative to the named conditions.",
  rank_statistic = "Statistic the ranked list was ordered by, using the token the canonical run recorded: moderated_t (canonical), log2fc (sensitivity) or not_applicable (ORA). The log2fc token is retained for provenance; the statistic it names is the standardised-abundance effect size, NOT a log2 fold change.",
  rank_source_column = "Column of the differential table used as the ranking statistic.",
  analysis_role = "canonical, sensitivity, or alternative_query_list. Canonical and sensitivity rows must never be pooled.",
  gene_identifier_type = "Identifier namespace of the gene members: UNIPROT.",
  ora_universe_definition = "The background universe used for ORA.",
  # EWCE
  ewce_analysis_type = "Baseline (per-condition abundant protein lists) or Differential (per-contrast up/down lists).",
  protein_list_definition = "Which protein list the cell-type enrichment was computed on.",
  cell_type = "Cell type from the reference cell-type dataset.",
  annotation_level = "Cell-type annotation level (1 or 2). The primary setting is 2.",
  top_n = "Size of the protein list tested. The primary setting is 250.",
  n_hits = "Number of proteins in the tested list.",
  n_background = "Size of the measured-proteome background (4859).",
  fold_change = "EWCE fold change relative to the bootstrap null: the ratio of observed to null mean expression proportion. This is a genuine fold change belonging to EWCE and is unrelated to the differential effect size in effect_size_sd_units.",
  sd_from_mean = "Standard deviations from the bootstrap null mean.",
  adjusted_p_value_within_target = "FDR-adjusted p-value within the target list.",
  adjusted_p_value_global = "FDR-adjusted p-value across all tests. This is the significance criterion used.",
  significant_global_fdr_0_05 = "TRUE if globally FDR-significant.",
  ewce_direction = "Direction of the EWCE effect as reported by the run.",
  bootstrap_replicates = "Number of bootstrap replicates (10000).",
  is_primary_setting = "TRUE for the primary setting (top-N 250, annotation level 2).",
  # coverage / parameters
  rows_exported = "Rows this release exports for the analysis block.",
  term_count_recorded_by_run = "Term count the canonical run recorded for the same block.",
  agrees_with_run_record = "TRUE when the exported row count equals the canonical record.",
  gsea_untested_gene_sets = "GSEA blocks only: gene sets in the size limits that fgsea gave no P value (unbalanced statistics), so they are absent from the block. GO sets are listed one per row in enrichment/gsea_untested_gene_sets.tsv; KEGG is counted from its run's warnings.",
  n_perm_simple = "fgsea nPermSimple of the run that could not test the set.",
  parameter = "Parameter name.",
  value = "Parameter value.",
  stage = "Analysis stage the parameter belongs to.",
  recorded_by = "Artefact that records the value.",
  evidence_path = "Path to the artefact that records the value. Files in the code repository are named relative to its proteomic_analysis/ folder.",
  evidence_sha256 = "SHA256 of that artefact.",
  applies_to_stage = "Analysis stage the version applies to.",
  notes = "Additional caveat or clarification.",
  category = "Kind of component: language, R package, external application, instrument, reference database.",
  component = "Software, database or instrument component.",
  version = "Exact version, or UNKNOWN where no version is recorded.",
  # figures
  figure = "Figure the panel belongs to.",
  panel = "Panel letter.",
  panel_label = "Panel label in the manuscript numbering (Fig. 2 X, Supp. Fig. 4 X). The figure and panel columns keep the internal keys (Figure3, Supplementary_proteomics).",
  panel_title = "Panel title as it appears in the manuscript.",
  source_table = "Release table the panel values can be reproduced from.",
  source_rows_filter = "Row filter to apply to that table to obtain the panel data.",
  statistical_unit = "Unit the revised panel is computed on.",
  original_statistical_unit = "Unit the original panel was computed on.",
  corrected_source_script = "Script that regenerated the panel at animal level.",
  final_revision_script = "Repository snapshot of the script that is the last word for this panel.",
  revision_status = "Outcome of the revision for this panel.",
  displayed_statistics_changed = "Whether the numbers shown on the panel changed.",
  biological_interpretation_changed = "Whether the biological reading of the panel changed. For the panels rebuilt on the per-sample-class statistics (revision_status REBUILT_WITHIN_COMPARTMENT*: Fig. 2 D, F, G, H, I and Supp. Fig. 4 F) the comparison is with the pooled-model release of 2026-09-07; for every other panel it is with the original submission.",
  key_numbers_original = "Key numbers on the original panel.",
  key_numbers_revised = "Key numbers on the revised panel.",
  key_numbers_pooled_model = "For the panels rebuilt on the per-sample-class statistics: the key numbers the same panel showed in the pooled-model release of 2026-09-07 (superseded), kept so the change is visible. Blank for panels that were not rebuilt.",
  panel_source_data_release_path = "Where the panel source data sits in this release.",
  panel_source_data_original_filename = "Original filename of the panel source data in the canonical revision output.",
  panel_source_data_canonical_path = "Canonical path the panel source data was copied from.",
  panel_source_data_sha256 = "SHA256 of the panel source data.",
  interpretation_note = "What a reader must know to read this panel correctly.",
  # workbook
  sheet = "Workbook sheet name.",
  supplementary_table = paste("Assigned supplementary-table number for the sheet (S1, S2,",
                              "...). Assigned, not positional, so inserting a sheet cannot",
                              "renumber a table that has already been cited. Blank for the",
                              "README sheet."),
  title = "Sheet title as printed in the workbook and listed on the Contents sheet.",
  n_rows = "Rows on the sheet.",
  n_columns = "Columns on the sheet.",
  columns = paste("Semicolon-separated MACHINE column names on the sheet -- the naming the",
                  "tests and this dictionary key off, matching the .tsv siblings."),
  display_columns = paste("Semicolon-separated DISPLAY headers actually printed on the",
                          "sheet, in the same order as `columns`. The workbook is the",
                          "reader-facing artefact and carries these; the .tsv files carry",
                          "the machine names."),
  table_level_constants = paste("Columns whose value was identical on every row and were",
                                "therefore hoisted out of the sheet into its caption, as",
                                "`name=value` pairs. They are properties of the table",
                                "rather than variables in it; all of them remain present in",
                                "the corresponding .tsv."),
  # SDRF
  sdrf_field = "SDRF-Proteomics column, or a non-column deposition requirement.",
  sdrf_requirement = "Whether the field is required or optional for an SDRF submission.",
  value_or_reason = "The value used, or a short reason when the field is missing.",
  why_missing = "Why the field could not be populated.",
  likely_source = "The document or party that would supply the value.",
  required_before_submission = "Whether the field must be filled before PRIDE submission.",
  key = "Name of the reported quantity.",
  # lineage
  artifact_id = "Stable identifier for the artefact within the lineage.",
  artifact_type = "What kind of artefact it is.",
  path_or_release_filename = "Where the artefact lives, or its name in this release. Files in the code repository are named relative to its proteomic_analysis/ folder (for example 07_publication_release/metadata/experimenter_metadata.tsv).",
  parent_artifact = "Semicolon-separated artifact_id(s) this artefact was derived from.",
  processing_stage = "The operation that produced this artefact. UNRESOLVED when it cannot be proven.",
  software = "Software that performed the step.",
  software_version = "Version of that software, or UNKNOWN.",
  sha256 = "SHA256 of the artefact, where a single file could be hashed.",
  canonical = "Whether this artefact is the canonical one for its stage.",
  # manifest
  relative_path = "Path of the file relative to the release root.",
  file_role = "What the file is for.",
  format = "File format.",
  rows_if_tabular = "Row count for tabular files.",
  columns_if_tabular = "Column count for tabular files.",
  source_artifact = "Semicolon-separated canonical artefact(s) the file was derived from. A relative path that starts with a release folder (metadata/, provenance/ and so on) names a file in this release; one that starts with a repository folder (for example 07_publication_release/, R/ or tests/) names a file in the proteomic_analysis/ folder of the code repository.",
  generated_by = "Release script that produced the file.",
  git_commit = "Repository commit the release was built from; a +dirty suffix means the working tree differed from that commit.",
  timestamp_utc = "UTC time the file was described.",
  bytes = "File size in bytes.",
  # dictionary itself
  table = "Release file the column belongs to.",
  column = "Column name.",
  description = "What the column holds.",
  # checks
  check = "Name of the validation check.",
  pass = "Whether the check passed."
)

# --------------------------------------------------------------------------------------
# build the dictionary
# --------------------------------------------------------------------------------------

TABLES <- c(
  "metadata/sample_metadata.tsv", "metadata/animal_level_sample_metadata.tsv",
  "metadata/metadata_field_provenance.tsv", "metadata/sample_class_corrections.tsv",
  "metadata/experimenter_metadata.tsv", "metadata/sample_preparation_protocol.tsv",
  "metadata/primary_contrast_manifest.tsv", "metadata/secondary_analysis_manifest.tsv",
  "processed_data/protein_abundance_measurement_level.tsv.gz",
  "processed_data/protein_abundance_animal_level.tsv.gz",
  "processed_data/protein_feature_annotation.tsv.gz",
  "differential_analysis/primary_differential_proteins.tsv.gz",
  "differential_analysis/primary_differential_summary.tsv",
  "enrichment/primary_GSEA_GO_BP.tsv.gz", "enrichment/primary_GSEA_GO_CC.tsv.gz",
  "enrichment/primary_GSEA_KEGG.tsv.gz", "enrichment/primary_ORA_GO_BP.tsv.gz",
  "enrichment/primary_ORA_GO_CC.tsv.gz", "enrichment/GSEA_log2FC_sensitivity.tsv.gz",
  "enrichment/primary_EWCE.tsv.gz", "enrichment/enrichment_coverage.tsv",
  "enrichment/gsea_untested_gene_sets.tsv", "enrichment/enrichment_run_parameters.tsv",
  "editor_source_data/figure_source_map.tsv",
  "editor_source_data/workbook_sheet_inventory.tsv",
  "editor_source_data/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv",
  "pride/sdrf.tsv", "pride/sdrf_field_status.tsv", "pride/pride_readiness.tsv",
  "provenance/data_lineage.tsv", "provenance/software_versions.tsv",
  "provenance/analysis_parameters.tsv"
)

# The two abundance matrices are wide: one column per sample. Describing 144 sample
# columns individually would bury the dictionary, so they get one generic row each.
WIDE <- c("processed_data/protein_abundance_measurement_level.tsv.gz",
          "processed_data/protein_abundance_animal_level.tsv.gz")
WIDE_NOTE <- c(
  "processed_data/protein_abundance_measurement_level.tsv.gz" =
    "One column per acquisition run, named by the acquisition run name (see metadata/sample_metadata.tsv, column sample_id). Values are processed protein abundance, filtered, imputed and standardised per protein: each protein has a mean of about 0 and an SD of about 1 over the 96 acquisitions. Not log2 abundance; the producing step is unresolved (see provenance/UPSTREAM_PREPROCESSING_GAP.md). NOT raw MS data.",
  "processed_data/protein_abundance_animal_level.tsv.gz" =
    "One column per animal-level unit, named `<AnimalID>_<sample_class>` (see metadata/animal_level_sample_metadata.tsv, column animal_level_column_name). Values are the equal-weight mean of the Left and Right measurements on the same per-protein standardised scale: each protein keeps a mean of about 0, and its SD over the 48 units is below 1 because two hemispheres are averaged. Not log2 abundance.")

# SDRF columns are ontology-style names already documented field by field in
# pride/sdrf_field_status.tsv. Sourcing their definitions from that table keeps one
# description per field instead of two that can disagree.
sdrf_status <- read.delim(rel("pride", "sdrf_field_status.tsv"), sep = "\t", quote = "",
                          stringsAsFactors = FALSE, check.names = FALSE)
sdrf_description <- function(cn) {
  hit <- sdrf_status[sdrf_status$sdrf_field == cn, , drop = FALSE]
  if (!nrow(hit)) return(NA_character_)
  parts <- c(paste0("SDRF-Proteomics field (", hit$sdrf_requirement[1], "; ",
                    hit$status[1], "). Value: ", hit$value_or_reason[1], "."))
  if (!is.na(hit$why_missing[1]) && nzchar(hit$why_missing[1])) {
    parts <- c(parts, paste0("Missing because: ", hit$why_missing[1], "."))
  }
  if (!is.na(hit$notes[1]) && nzchar(hit$notes[1])) {
    parts <- c(parts, hit$notes[1])
  }
  paste(parts, collapse = " ")
}

# Per-table descriptions for columns whose meaning depends on the table (for example
# `status` and `notes`); DESCRIPTIONS is keyed by column name alone.
TABLE_COLUMN_OVERRIDES <- list(
  "metadata/experimenter_metadata.tsv" = c(
    metadata_id = "Stable identifier used by release builders to retrieve one experimenter-supplied fact.",
    category = "Fact category: scope, biological source, sample-class definition, animal cohort, proteomics method, unresolved deposition metadata, raw-data availability or scope exclusion.",
    sample_class = "Canonical sample-class token when the row defines one of the four marker-defined LCM classes; otherwise blank.",
    sort_order = "Display or protocol order where order is meaningful; otherwise blank.",
    value = "Experimenter-supplied value preserved in human-readable form.",
    sdrf_value = "Standards-oriented SDRF representation when the fact belongs in the SDRF; blank for cohort-only or scope-only facts.",
    status = "Metadata status under the release vocabulary: verified, recoverable/unknown, or not applicable.",
    applies_to = "Scope to which the fact applies; distinguishes study/cohort facts from per-sample metadata.",
    evidence_source = "Provenance class for the fact, here experimenter-supplied metadata or protocol.",
    evidence_date = "Date the experimenter-supplied fact was recorded in the release metadata file.",
    notes = "Interpretation boundary, exclusion or non-inference rule attached to the fact."
  ),
  "metadata/sample_preparation_protocol.tsv" = c(
    protocol_id = "Stable identifier for the sample-preparation protocol.",
    protocol_version_date = "Date of the experimenter-supplied protocol version.",
    step_order = "One-based order of the preparation step.",
    step_name = "Short operation name for the preparation step.",
    details = "Reagents, amounts, concentrations or action exactly as supplied.",
    temperature = "Temperature recorded for the step, or not applicable.",
    duration = "Duration recorded for the step; `overnight` is preserved verbatim and unavailable exact durations are not invented.",
    status = "Metadata status for the protocol step; all seven supplied steps are KNOWN_VERIFIED.",
    evidence_source = "Provenance class for the step: experimenter-supplied protocol.",
    evidence_date = "Date the protocol was entered into the release metadata file.",
    notes = "Interpretation boundary or unresolved detail for the step."
  ),
  "editor_source_data/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv" = c(
    file = "File or released artefact the occurrence is in.",
    location = "Where in that file: a named constant, section, column or line number.",
    old_term = "The historical fold-change wording that was found at this location.",
    new_term = "What it now reads, or the same token where it was kept for provenance.",
    status = paste("Classification of the occurrence. INTERNAL_PROVENANCE: a canonical",
                   "filename, column value or source-field citation, kept",
                   "because renaming it would break the link to the canonical run.",
                   "PUBLIC_MISLABEL: reader-facing text that described the published effect",
                   "size as a fold change; corrected. LEGITIMATE_OTHER_CONTEXT: a genuine",
                   "fold change belonging to a different statistic (EWCE's bootstrap fold",
                   "change), or a frozen historical artefact that must not be edited."),
    notes = "Why the occurrence was classified that way, and what was done about it."
  )
)

rows <- list()
undocumented <- character(0)
for (t in TABLES) {
  p <- rel(dirname(t), basename(t))
  if (!file.exists(p)) {
    stop("Release table missing when building the dictionary: ", t, call. = FALSE)
  }
  cols <- names(read_release(p))
  cols <- cols[nzchar(trimws(cols))]
  if (t == "pride/sdrf.tsv") {
    for (cn in cols) {
      d <- sdrf_description(cn)
      if (is.na(d)) undocumented <- c(undocumented, paste0(t, " :: ", cn))
      rows[[length(rows) + 1L]] <- data.frame(table = t, column = cn, description = d,
                                              stringsAsFactors = FALSE, check.names = FALSE)
    }
    next
  }
  if (t %in% WIDE) {
    rows[[length(rows) + 1L]] <- data.frame(
      table = t, column = "protein_group_id",
      description = unname(DESCRIPTIONS[["protein_group_id"]]),
      stringsAsFactors = FALSE, check.names = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      table = t, column = "<one column per sample>",
      description = unname(WIDE_NOTE[[t]]),
      stringsAsFactors = FALSE, check.names = FALSE)
    next
  }
  for (cn in cols) {
    ovr <- TABLE_COLUMN_OVERRIDES[[t]]
    d <- if (!is.null(ovr) && cn %in% names(ovr)) unname(ovr[[cn]])
         else if (cn %in% names(DESCRIPTIONS)) unname(DESCRIPTIONS[[cn]])
         else NA_character_
    if (is.na(d)) undocumented <- c(undocumented, paste0(t, " :: ", cn))
    rows[[length(rows) + 1L]] <- data.frame(table = t, column = cn, description = d,
                                            stringsAsFactors = FALSE, check.names = FALSE)
  }
}
if (length(undocumented)) {
  stop("Columns with no data-dictionary entry:\n",
       paste("  -", unique(undocumented), collapse = "\n"), call. = FALSE)
}
dictionary <- do.call(rbind, rows)
rownames(dictionary) <- NULL
release_log("  data dictionary: ", nrow(dictionary), " column definitions over ",
            length(TABLES), " tables")

# --------------------------------------------------------------------------------------
# README_DATA.md
# --------------------------------------------------------------------------------------

n_with_hits <- sum(as.integer(diff_summary$n_significant_fdr_0_05) > 0L)
mch <- diff_summary[diff_summary$canonical_comparison ==
                      "mcherry_paired_veh_over_mcherry_unpaired_veh", , drop = FALSE]

R <- c(
"# Associative Memory Proteomics -- publication data package",
"",
paste0("Animal-level proteomics of four marker-defined laser-capture microdissection classes from ",
       mv("organism_part"), " in a chemogenetic associative-memory"),
"experiment. This package contains the processed data, the differential and enrichment",
"results, the manuscript figure source data, and the metadata needed to interpret and",
"redeposit them.",
"",
"**Read this first:** the independent experimental unit is the **animal**, not the",
"hemisphere. An earlier version of this analysis treated left and right hemispheres as",
"independent replicates. That was pseudoreplication and has been corrected throughout. See",
"`editor_source_data/REVISION_PROTEOMICS_DATA_CHANGELOG.md` for what changed.",
"",
paste0("**Statistics.** Differential abundance is tested with limma, one model per sample",
       " class (residual df ", stat_run$df_residual, "). This package replaces the release of",
       " 2026-09-07, which used one pooled model over all 48 animal-level units; section 7b of",
       " the changelog lists what changed."),
"",
"---",
"",
"## 1. Experimental design",
"",
paste0("12 ", mv("strain_or_breed"), " mice from ", mv("animal_supplier"), ", ",
       mv("animal_supplier_location"), ". All proteomics material came from ",
       mv("organism_part"), " and was obtained by ", mv("collection_method"), "."),
"",
paste0("Animals were ", mv("age_at_experiment_start"), " at experiment start. The cohort contained ",
       mv("cohort_sex_composition"), ". Reliable AnimalID-level sex and age at tissue collection ",
       "are unavailable and are not inferred."),
"",
"From each animal, four marker-defined classes were collected from each hemisphere:",
"",
"| Sample class (data value) | Display name | Historical alias | What it is |",
"|---|---|---|---|",
class_markdown,
"",
"These are sampling classes, not four Cell Ontology cell types. The SDRF therefore uses",
"`characteristics[cell type] = not applicable`, while `characteristics[sampling site]` and",
"`factor value[sample class]` preserve the marker-defined LCM category.",
"",
"Each animal is in exactly one of four conditions:",
"",
"| Code | Condition | Pairing | Treatment |",
"|---|---|---|---|",
"| 1 | `paired_cno` | paired | cno |",
"| 2 | `paired_veh` | paired | veh |",
"| 3 | `unpaired_cno` | unpaired | cno |",
"| 4 | `unpaired_veh` | unpaired | veh |",
"",
"12 animals x 4 sample classes x 2 hemispheres = **96 acquisition measurements**.",
"3 animals per condition, so **n = 3 animals** in every condition x sample class stratum,",
"with complete left/right pairing in all 16 strata.",
"",
"## 2. The statistical unit correction",
"",
"Left and right measurements were averaged with equal weight within each",
"`AnimalID x sample_class` unit, on the existing processed values (standardised per protein,",
"see section 3), **before any",
"statistical test**. 96 measurements become **48 animal-level units**.",
"",
"In a balanced design the mean of animal means equals the mean of all hemispheres, so point",
"estimates barely move. What changes is the inference: standard errors, degrees of freedom",
"and every p-value.",
"",
"## 3. Measurement-level versus animal-level data",
"",
"| File | Shape | Use it for |",
"|---|---|---|",
"| `processed_data/protein_abundance_measurement_level.tsv.gz` | 5,349 x 96 | acquisition-level QC only |",
"| `processed_data/protein_abundance_animal_level.tsv.gz` | 5,349 x 48 | **all inference** |",
"",
"### The value scale, and what the effect size therefore means",
"",
"**These matrices are standardised per protein, not log2 abundance.** Measured directly",
"on the measurement-level matrix: each protein is standardised separately across the 96",
"acquisitions, giving row means of approximately 0 and row standard deviations of",
"approximately 1. The build asserts this and refuses to run if it ever stops holding.",
paste0("The animal-level matrix is the mean of the two hemispheres on that scale: each protein",
       " keeps a mean of approximately 0, and its standard deviation over the 48 units is below",
       " 1 (", sprintf("%.2f", min(animal_row_sd)), " to ", sprintf("%.2f", max(animal_row_sd)),
       ", median ", sprintf("%.2f", stats::median(animal_row_sd)), ") because averaging two",
       " hemispheres reduces the spread."),
"",
"The consequence matters for reading the results. The statistical GCT stores the",
"coefficient as `logFC` (limma's name for it) and the canonical split tables call it",
"`log2fc`, but a difference of group means",
paste0("on a protein-standardised scale is a **", RELEASE_EFFECT_SIZE$public_term, "** --"),
paste0("equivalently, a ", RELEASE_EFFECT_SIZE$public_term_alt, ". It is not a log2 ratio."),
"An effect size of 2 means \"two standard deviations on the standardised scale\", not",
paste0("\"four-fold\". This package therefore publishes the column as `",
       RELEASE_EFFECT_SIZE$public_field, "`, with `effect_size_definition`,"),
paste0("`effect_size_source_column` and `source_statistic_field = \"",
       RELEASE_EFFECT_SIZE$source_field, "\"` on every row, rather than propagating a"),
"fold-change reading the data do not support. **The exported values are bit-identical to",
"the canonical `logFC`** -- the build re-reads every split table and requires exact",
"equality -- so only the label changed.",
"",
"Two phrasings are avoided:",
"",
"- **not** \"standardised mean difference\" on its own, because that invites a Cohen's d",
"  reading, i.e. scaling by a pooled *within-group* standard deviation. The scaling here is",
"  a per-protein standard deviation taken across the whole measurement dataset, before any",
"  grouping -- a different quantity.",
"- **not** \"z-scored\". The released values carry two decimals, so on those bytes no",
paste0("  protein is exactly mean 0 / standard deviation 1; the residual deviations (max |row",
       " mean| ", format(signif(std_meas$max_abs_row_mean, 3), scientific = FALSE),
       ", max |row sd - 1| ",
       format(signif(std_meas$max_abs_row_sd_minus_1, 3), scientific = FALSE),
       ", measured on the measurement-level",
       " matrix) are the size 2-decimal rounding produces. The"),
"  operation that produced the matrix is also unresolved. \"Standardised\" is what the",
"  numbers support.",
"",
paste0("Corroboration: `average_standardized_abundance` (limma `AveExpr`) is the mean over the ",
       stat_run$n_columns_in_model, " animal-level units of one sample class, so within a class",
       " it is not centred (it ranges from ", sprintf("%.2f", ave$range_within_classes[1]),
       " to ", sprintf("%.2f", ave$range_within_classes[2]), "). Averaged over the four sample",
       " classes it stays close to zero (largest absolute value over all ",
       format(ave$n_proteins, big.mark = ","), " proteins ",
       format(signif(ave$max_abs_mean_over_classes, 3), scientific = FALSE),
       "), which is only possible on a per-protein centred input."),
"",
"For methods sections:",
"",
paste0("> ", RELEASE_EFFECT_SIZE$methods_sentence),
"",
"The upstream reason for the standardisation is **not documented** -- see section 15. The",
"nearest thing to a record is a sibling file at the project root named",
"`pg.matrix_filtered_70percent-onegroup_imputed_ANOVA_z-scored.txt`, whose values match",
"this matrix to within rounding (identical protein order, 96 identical sample columns,",
"maximum absolute difference 0.005 -- the signature of a 2-decimal serialisation). Its",
"filename is consistent with what the numbers show, but a filename is not a processing",
"record, and the matrix in use is named `pcaAdjusted_unnormalized`, which is consistent",
"with neither.",
"",
"Do not run inference on the measurement-level matrix. It is provided because two",
"supplementary QC panels legitimately need the acquisition as their unit (hemisphere and",
"collection-plate QC cannot be shown at animal level, because neither variable exists in",
"the 48-row metadata).",
"",
"Neither file is raw mass-spectrometry data. See section 20.",
"",
"## 4. The 12 primary contrasts",
"",
"4 sample classes x 3 within-sample-class contrast families:",
"",
"1. `paired_cno` vs `paired_veh` -- chemogenetic effect in paired animals",
"2. `paired_veh` vs `unpaired_veh` -- pairing-associated (learning) signature",
"3. `unpaired_cno` vs `unpaired_veh` -- chemogenetic effect in unpaired controls",
"",
"All 12 are in `metadata/primary_contrast_manifest.tsv`. A positive",
"`effect_size_sd_units` always means",
"higher in the numerator condition.",
"",
"## 5. Secondary analyses",
"",
"Kept in a separate file, `metadata/secondary_analysis_manifest.tsv`, and never mixed into",
"the primary table: the Supp. Fig. 4 E cross-compartment within-animal paired model, the",
"effect-size-ranked GSEA sensitivity analysis, the post-hoc Pairing x CNO interaction, the",
"leave-one-animal-out robustness check of the mCherry learning comparison, the EWCE",
"baseline analysis, the superseded pooled-model statistics (retained for comparison only)",
"and two candidate contrasts that were not analysed.",
"",
"## 6. Differential statistics",
"",
"`differential_analysis/primary_differential_proteins.tsv.gz` -- one row per protein group",
paste0("per comparison, 12 x 5,349 = 64,188 rows. limma ", stat_run$limma_version, " (R ",
       stat_run$r_version, "), one linear model per sample class (~0 + condition over that",
       " class's ", stat_run$n_columns_in_model, " animal-level units, n = ",
       RELEASE_DESIGN_INVARIANTS$n_animals_per_stratum, " animals per condition, residual df ",
       stat_run$df_residual, "), empirical Bayes moderation (robust, no intensity trend),",
       " Benjamini-Hochberg FDR within each comparison. Each comparison is 3 vs 3 animals.",
       " Model parameters, including the prior variance of each sample class, are in",
       " `provenance/analysis_parameters.tsv`."),
"",
"Two protein counts appear throughout and mean different things:",
"",
"- **5,349** protein groups were tested.",
"- **5,327** carry a mouse identifier and could be annotated and used for enrichment. The",
"  22-row difference is non-mouse contaminant entries (human keratins, porcine trypsin).",
"  The manuscript figures were drawn from the 5,327 set, so both counts are reported side",
"  by side in `primary_differential_summary.tsv`.",
"",
sprintf("Of the 12 primary comparisons, %d %s at least one FDR-significant protein and %d %s none.",
        n_with_hits, if (n_with_hits == 1L) "contains" else "contain",
        12L - n_with_hits, if (12L - n_with_hits == 1L) "contains" else "contain"),
"",
"## 7. GSEA",
"",
"**Canonical GSEA ranks by the moderated t statistic.** The effect-size-ranked GSEA is a",
"sensitivity",
"analysis only and lives in `enrichment/GSEA_log2FC_sensitivity.tsv.gz` with",
"`analysis_role = sensitivity` on every row. Never pool the two.",
"",
"Settings: GO-BP, GO-CC and KEGG, minGSSize 10, maxGSSize 800, BH adjustment,",
"`pvalue_cutoff = 1` and `qvalue_cutoff = 1` -- so the exported tables hold **every tested",
"term** that fgsea could estimate, not a significant-only subset. Filter on",
"`adjusted_p_value` yourself. GSEA P values come from fgsea permutations with one recorded",
"seed per comparison; term counts near FDR 0.05 can move by a few terms with the seed.",
"KEGG gene sets were downloaded at run time and are not version-pinned.",
"",
paste0("**Permutations and untested gene sets.** ", permutations$sentence, " Gene sets",
       " fgsea could not test are absent from the tables rather than non-significant (section",
       " 13); the GO sets are listed in `enrichment/gsea_untested_gene_sets.tsv`. The figure",
       " panels that draw on such blocks: Fig. 2 F (",
       untested_blocks("paired_veh_vs_unpaired_veh"), "), Fig. 2 I (",
       untested_blocks(c("paired_veh_vs_unpaired_veh", "paired_cno_vs_paired_veh")),
       "), Supp. Fig. 4 F (", untested_blocks("paired_cno_vs_paired_veh", "GO_BP", sf_classes),
       ") and Fig. 2 G (",
       untested_blocks("paired_veh_vs_unpaired_veh", c("GO_BP", "GO_CC"), "mcherry"), ")."),
"",
paste("GO-CC is published as `enrichment/primary_GSEA_GO_CC.tsv.gz` and",
      "`enrichment/primary_ORA_GO_CC.tsv.gz`.", RELEASE_GO_CC_PROVENANCE_NOTE),
"",
"## 8. ORA universe",
"",
"All unique successfully mapped measured UniProt accessions, computed per comparison. Four",
"query-list variants are exported and distinguished by the `query_list` column: all",
"FDR-significant proteins, FDR-significant higher in the numerator, FDR-significant higher",
"in the denominator, and an alternative list defined by the absolute effect size (the",
"canonical run named that variant `top_absolute_log2fc`).",
"",
"Where a comparison had no FDR-significant proteins there is no query list and hence no ORA",
"rows. `enrichment/enrichment_coverage.tsv` records every such block, so an",
"absent block reads as a result and not as a missing file.",
"",
"## 9. Identifier mapping",
"",
"UniProt `MOUSE_10090_idmapping.dat`, snapshot 2025-11-03, plus 11 hand-curated overrides.",
"Both are hashed in `provenance/software_versions.tsv` and",
"`processed_data/protein_feature_annotation.tsv.gz`.",
"",
"**Field naming.** The internal pipeline carries a column named `gene_symbol` that actually",
"holds a UniProt identifier -- a UniProt entry name in the split tables and a UniProt",
"accession in the mapped tables. That name is **not** used in this package. Here:",
"",
"| Published column | Holds | Example |",
"|---|---|---|",
"| `protein_group_id` | UniProt entry name(s) | `MFA1A_MOUSE;MFA1B_MOUSE` |",
"| `protein_group_uniprot_accessions` | all accessions in the group | `C0HKD8;C0HKD9` |",
"| `uniprot_accession` | leading accession | `C0HKD8` |",
"| `gene_symbol` | resolved gene symbol | `Mfap1a` |",
"| `gene_symbols_in_group` | all symbols in the group | `Mfap1a;Mfap1b` |",
"| `protein_description` | real protein description | `Microfibrillar-associated protein 1A` |",
"",
"## 10. EWCE",
"",
"`enrichment/primary_EWCE.tsv.gz`. Bootstrap cell-type enrichment, 10,000 replicates,",
"measured-proteome background of 4,859. Primary setting: annotation level 2, top-N 250,",
"global FDR 0.05; other top-N and annotation-level settings are exported with",
"`analysis_role = sensitivity`. The cell-type reference is `ewceData::ctd()`.",
"",
paste(ewce_model$sentence, "Joining EWCE to the differential tables on",
      "`canonical_comparison` matches the comparison, not the protein ranking."),
"",
"Note: an `l1_amygdala.loom` single-cell file exists in the project inputs but is **not**",
"used by the canonical EWCE run, and it is not evidence about which brain region was",
"dissected.",
"",
"## 11. PCA",
"",
"Animal-level PCA over all 5,349 proteins and 48 units, 0 zero-variance rows removed.",
"PC1 29.1%, PC2 20.5%, PC3 12.1%. Acquisition-level PCA is retained only for the two",
"technical QC panels.",
"",
"## 12. Figure source data",
"",
"`editor_source_data/figure_source_map.tsv` maps every manuscript proteomics panel to the",
"release table, row filter, statistical unit and generating script behind it. The per-panel",
"source-data files are copied byte-for-byte into",
"`editor_source_data/figure_source_data/`.",
"",
paste0("Panels and their source-data files are named in the manuscript numbering (Fig. 2,",
       " Supp. Fig. 4; e.g. `Fig2D_source_data.csv`, `SuppFig4F_source_data.csv`); the",
       " canonical file each was copied from is named in",
       " `panel_source_data_original_filename`. Fig. 2 D, F, G, H, I and Supp. Fig. 4 F",
       " were rebuilt on the per-sample-class statistics in `",
       paste(RELEASE_STATISTICS_RUN$figure_rebuild_relative, collapse = "/"),
       "`; the other panels do not use the differential statistics and come from the",
       " 2026-08-27 animal-level revision."),
"",
"## 13. Known design limitations",
"",
if (nrow(corrections) > 0L) c(
sprintf(paste("- **SAMPLE-CLASS CORRECTION, %d of 96 measurements -- RESOLVED.** Six",
              "acquisition-level samples from the left hemispheres of animals %s had",
              "sample-class assignments corrected during historical sample-identity quality",
              "control. Original plate-layout and early metadata records retain the",
              "pre-correction assignments, whereas a preserved UMAP-based correction table",
              "records a cyclic reassignment among mCherry-, neuropil- and neuron-enriched",
              "samples. Acquisition identifiers and quantitative protein-abundance profiles",
              "were unchanged; only sample-class metadata were reassigned. Independent",
              "profile-similarity and left/right-pair analyses supported the corrected",
              "assignments. Both original and analysis-time sample-class labels are",
              "retained in the released metadata."),
        nrow(corrections),
        paste(sort(unique(corrections$AnimalID)), collapse = " and ")),
"",
paste0("  The reassignment is the same 3-cycle for both animals: **neuropil -> mCherry,",
       " mCherry -> neuron, neuron -> neuropil**. The cfos-enriched samples were not",
       " involved. Per-row provenance, including the correction record and its SHA256, is in",
       " `metadata/sample_class_corrections.tsv`; `original_sample_class`,",
       " `analysis_sample_class` and `sample_class_corrected` are on every row of",
       " `metadata/sample_metadata.tsv`. Verdict: `",
       unique(corrections$correction_status), "`."),
"",
paste0("  **The correction table is itself the preserved historical correction record.** No",
       " surviving prose note explains the rationale -- that absence is stated rather than",
       " glossed. It does not make the correction unresolved: nine independent",
       " pre-correction records carry the original assignments, three of which the build",
       " re-derives and requires to agree, and the correction table's hash is recomputed and",
       " checked against the value the audit recorded."),
"",
paste0("  **This is not a mislabelling and requires no reanalysis.** The affected",
       " animal-level units (",
       paste(sort(unique(corrections$animal_level_unit_affected)), collapse = ", "),
       ") were built from the corrected assignments in the first place, so the validated",
       " results already reflect them. No abundance value changed: the quantitative matrix",
       " was compared across the metadata edit and the maximum absolute numeric difference",
       " was 0. Full account in",
       " `editor_source_data/REVISION_PROTEOMICS_DATA_CHANGELOG.md`, section 11.")
) else
"- No sample-class correction is recorded: all 96 measurements carry their original class.",
"- **n = 3 animals per group.** Small, and the corrected analysis reflects that.",
paste0("- **One model per sample class, residual df ", stat_run$df_residual, ".** Each sample",
       " class has its own variance estimate. Results for different sample classes come from",
       " the same 12 animals and are not independent of one another, so differences between",
       " sample classes in what reaches significance are not themselves tested."),
"- **The mCherry learning result is unresolved.** ",
local({
  mc <- "mcherry_paired_veh_over_mcherry_unpaired_veh"
  up <- function(tbl) {
    tbl[tbl$canonical_comparison == mc & tbl$query_list == "fdr_significant_higher_in_numerator",
        , drop = FALSE]
  }
  ou <- up(ora_go_all); oc <- up(ora_cc_all)
  padj <- function(d) as.numeric(d$adjusted_p_value)
  sprintf(paste("  It retains %s FDR-significant proteins of %s mapped (%s higher in",
                "paired_veh, %s higher in unpaired_veh). Moderated-t-ranked GSEA finds %d",
                "GO-BP and %d GO-CC terms at FDR < 0.05. Over-representation among the",
                "proteins higher in paired_veh finds %d of %s GO-BP terms (minimum adjusted P",
                "%.2f) and %d of %s GO-CC terms (minimum adjusted P %s) at FDR < 0.05.",
                "Pairing is completely associated with collection plate in this comparison",
                "and no technical batch metadata exist, so neither a biological nor a",
                "technical explanation is established."),
          mch$n_significant_fdr_0_05_mapped_only, mch$n_proteins_tested_mapped_only,
          mch$n_significant_higher_in_numerator_mapped_only,
          mch$n_significant_higher_in_denominator_mapped_only,
          sum(gsea_go_all$canonical_comparison == mc &
                as.numeric(gsea_go_all$adjusted_p_value) < 0.05, na.rm = TRUE),
          sum(gsea_cc_all$canonical_comparison == mc &
                as.numeric(gsea_cc_all$adjusted_p_value) < 0.05, na.rm = TRUE),
          sum(padj(ou) < 0.05, na.rm = TRUE), format(nrow(ou), big.mark = ","),
          min(padj(ou), na.rm = TRUE),
          sum(padj(oc) < 0.05, na.rm = TRUE), format(nrow(oc), big.mark = ","),
          format(signif(min(padj(oc), na.rm = TRUE), 2)))
}),
local({
  r <- function(col) range(as.numeric(loo_drop[[col]]))
  fdr <- r("n_fdr05"); bp <- r("go_bp_fdr05"); bpk <- r("go_bp_retained_of_full")
  cc <- r("go_cc_fdr05"); cck <- r("go_cc_retained_of_full")
  paste0("- **The mCherry learning counts depend on single animals.** Refitting with one",
         " animal left out (2 vs 3 animals, residual df ",
         paste(unique(as.integer(loo_drop$df_residual)), collapse = "/"),
         ") gives ", fdr[1], " to ", format(fdr[2], big.mark = ","),
         " FDR-significant proteins (full fit ",
         format(as.integer(loo_full$n_fdr05), big.mark = ","), ") and ", bp[1], " to ",
         bp[2], " GO-BP terms, keeping ", bpk[1], " to ", bpk[2], " of the full fit's ",
         as.integer(loo_full$go_bp_fdr05), ". GO-CC is more stable: ", cc[1], " to ", cc[2],
         " terms, keeping ", cck[1], " to ", cck[2], " of ", as.integer(loo_full$go_cc_fdr05),
         ". The effect estimates are stable (correlation with the full fit ",
         sprintf("%.2f", min(as.numeric(loo_drop$r_effect))), " to ",
         sprintf("%.2f", max(as.numeric(loo_drop$r_effect))), "). Source: `",
         paste(LOO_REL, collapse = "/"), "`.")
}),
local({
  can <- fgsea_na
  plural_sets <- function(n) if (n == 1L) " set" else " sets"
  if (!nrow(can)) {
    paste("- **GSEA permutations.**", permutations$sentence,
          "fgsea estimated a P value for every gene set in every block.")
  } else paste0(
         "- **Some gene sets could not be tested by GSEA.**",
         " When a comparison's moderated t values are",
         " strongly one-sided, fgsea cannot estimate a P value for some large gene sets",
         " (\"unbalanced\" statistics), and clusterProfiler drops them, so they are absent",
         " from the GSEA tables rather than non-significant. This affected ",
         nrow(can), " of the ", 3L * nrow(diff_summary), " comparison x collection blocks",
         " (GO-BP, GO-CC, KEGG; ", min(can$n_sets), " to ",
         max(can$n_sets), " sets each); for the mCherry learning",
         " comparison, ", if (nrow(fgsea_na_mch)) paste0(sum(fgsea_na_mch$n_sets), " ",
         paste(gsub("_", "-", sub("_GSEA$", "", fgsea_na_mch$analysis)), collapse = "/"),
         plural_sets(sum(fgsea_na_mch$n_sets))) else
           "none", ". The GO sets are listed in `enrichment/gsea_untested_gene_sets.tsv`;",
         " counts per block are in `enrichment/enrichment_coverage.tsv` and in section 5 of",
         " the changelog. ",
         permutations$sentence,
         local({
           affected <- c(
             "Fig. 2 F" = untested_blocks("paired_veh_vs_unpaired_veh") != "none",
             "Fig. 2 I" = untested_blocks(c("paired_veh_vs_unpaired_veh",
                                            "paired_cno_vs_paired_veh")) != "none",
             "Supp. Fig. 4 F" = untested_blocks("paired_cno_vs_paired_veh", "GO_BP",
                                                sf_classes) != "none")
           if (any(affected)) paste0(" The GO term counts behind ",
                                     paste(names(affected)[affected], collapse = ", "),
                                     " rest on blocks with untested sets (section 7).") else ""
         }))
}),
"- **KEGG gene sets are not version-pinned.** They were downloaded at run time, so a rerun",
"  can change KEGG results for reasons unrelated to the data.",
"- **Collection plate.** In the learning comparison pairing is completely associated with",
"  collection plate, so the two cannot be separated (section 14).",
"- **Fig. 2 H and Fig. 2 I compare contrasts that share an arm.** Both contain `paired_veh`,",
"  on opposite sides, so the design biases them toward opposition. Descriptive only.",
"- **Supp. Fig. 4 E is a secondary paired model** with residual df 6 before moderation.",
"- **No acquisition-level technical QC beyond identification counts.** The QC report for",
"  this acquisition was never produced for this project.",
"",
"## 14. Collection plate",
"",
PLATE_SENTENCE,
"",
NO_BATCH_SENTENCE,
"",
"The `Plate1`/`Plate2` token records **which plate a sample was collected onto**. It is not",
"a proteomics preparation, digestion, LC-MS, acquisition or instrument batch, and no such",
"batch metadata exist for this dataset. 11 of 12 animals follow paired -> Plate1 /",
"unpaired -> Plate2, and in the learning contrast the association is complete in every",
"sample class. At acquisition level, collection plate explains 0.01% of PC1 and 1.5% of",
"PC2, neither significant.",
"",
"## 15. Upstream preprocessing limitation",
"",
"The reproducible pipeline **starts** at the 5,349 x 96 processed matrix. The steps that",
"produced it from the 5,747-row search output -- filtering, log transformation,",
"normalisation, imputation and PCA adjustment -- are **not reproducible from this",
"repository**: no producing script survives, no seed was recorded for the stochastic",
"imputation, and the search-software output is absent from the project. The processed",
"matrix should be treated as a given input, cited by its SHA256.",
"",
"Full account: `provenance/UPSTREAM_PREPROCESSING_GAP.md`.",
"",
"## 16. Provenance and integrity",
"",
"| File | What it gives you |",
"|---|---|",
"| `provenance/data_lineage.tsv` | every artefact from acquisition to release, with unproven edges marked `UNRESOLVED` |",
"| `provenance/software_versions.tsv` | exact versions as recorded by the runs themselves, or `UNKNOWN` with the source that would supply them |",
"| `provenance/analysis_parameters.tsv` | every recorded analysis parameter |",
"| `provenance/release_manifest.tsv` | per-file role, shape, SHA256, source artefact and generating script |",
"| `provenance/SHA256SUMS.txt` | checksums for the whole package |",
"| `provenance/sessionInfo_release.txt` | the environment that built the package |",
"",
paste0(length(RELEASE_LOCKED_ARTEFACTS), " upstream artefacts are locked by hash and",
       " verified at build time:"),
"",
"```",
# Printed from RELEASE_LOCKED_ARTEFACTS so the README shows the current hashes.
vapply(names(RELEASE_LOCKED_ARTEFACTS), function(key) {
  spec <- RELEASE_LOCKED_ARTEFACTS[[key]]
  label <- switch(key, animal_level_input_gct = "animal-level matrix",
                  stat_gct = "limma statistics", stat_run_manifest = "statistics run manifest",
                  enrichment_index_go_bp = "GO-BP enrichment index",
                  enrichment_index_go_cc = "GO-CC enrichment index",
                  enrichment_index_kegg = "KEGG enrichment index",
                  key)
  paste0(formatC(label, width = -30), spec$sha256)
}, character(1)),
"```",
"",
"## 17. How this package maps to PRIDE",
"",
paste0("Current status: **`", pride_status, "`**."),
"",
"The 96 native acquisition files are `.d` directories. Their names are known for all 96 and",
"appear in `metadata/sample_metadata.tsv` (`raw_file_basename`) and `pride/sdrf.tsv`",
"(`comment[data file]`), but **the files themselves are not in this package and were not",
"found anywhere in the project tree**.",
"",
paste0("`pride/sdrf.tsv` follows SDRF-Proteomics conventions as far as the verified metadata",
       " allow. ", n_required_missing, " SDRF-required fields could not be populated from",
       " the project or experimenter records and carry `not available`; each is itemised in",
       " `pride/SDRF_MISSING_METADATA.md` with the document that would supply it."),
"",
paste0("Experimenter metadata resolves ", mv("organism_part"), ", ", mv("strain_or_breed"),
       ", ", mv("labeling_strategy"), " and sequential Lys-C plus trypsin digestion. The full",
       " 2024-11-28 protocol is preserved in `metadata/sample_preparation_protocol.tsv` and",
       " `pride/SAMPLE_PREPARATION_PROTOCOL.md`."),
"",
"No unresolved acquisition/search metadata was inferred. Instrument model, acquisition mode,",
"search-software version, modification parameters, mass tolerances, developmental stage,",
"AnimalID-level sex and age at tissue collection remain unavailable. CAA/TCEP",
"use is not converted into search modification settings.",
"",
"## 18. Package layout",
"",
"```",
"README_DATA.md",
"metadata/         sample and animal-level metadata, contrast manifests, data dictionary",
"processed_data/   abundance matrices and protein feature annotation",
"differential_analysis/  the 12 primary comparisons, long format",
"enrichment/       GSEA, ORA, EWCE, sensitivity, coverage and run parameters",
"editor_source_data/  the source-data workbook, changelog, figure map and panel data",
"pride/            SDRF and deposition notes",
"provenance/       lineage, versions, parameters, manifest and checksums",
"```",
"",
"## 19. Column definitions",
"",
"`metadata/data_dictionary.tsv` describes every column of every table in this package. It",
"is generated by walking the tables, and the build fails if a column has no definition.",
"",
"## 20. What is not here",
"",
"Native raw mass-spectrometry data. In particular, the file supplied with the original",
"submission as `processed_protein_group_matrix_raw.txt` is **not** raw mass-spectrometry",
"data: it is a 5,747-row protein-group matrix from the search software. See",
"`pride/README_PRIDE.md`.",
"",
"---",
"",
paste0("Generated by `", STAGE, "`."),
""
)
release_check_forbidden_wording(R)
if (!any(grepl(PLATE_SENTENCE, paste(R, collapse = " "), fixed = TRUE))) {
  stop("The required collection-plate sentence is not present verbatim in README_DATA.md.",
       call. = FALSE)
}
release_write_lines(R, release_path("README_DATA.md"))
release_log("  wrote README_DATA.md (", length(R), " lines); wording checks pass")

w <- release_write_table(dictionary, release_path("metadata", "data_dictionary.tsv"))
release_register("metadata/data_dictionary.tsv",
                 "definition of every column in every release table",
                 "all release tables", NA_character_, STAGE, "tsv")
release_register("README_DATA.md", "human entry point to the publication data package",
                 "all release tables", NA_character_, STAGE, "md")
release_log("  wrote data_dictionary.tsv (", w$rows, "x", w$cols, ")")
release_log("stage 11 complete")
