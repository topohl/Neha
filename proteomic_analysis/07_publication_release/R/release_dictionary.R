# Column definitions shared by stage 07 (the workbook's Column definitions sheet) and stage 11
# (metadata/data_dictionary.tsv).
#
#   RELEASE_DESCRIPTIONS        one description per machine column of the release .tsv files,
#                               keyed by column name; stage 11 walks the tables with it
#   RELEASE_READER_DEFINITIONS  the reader-facing definition and unit of every column shown in
#                               the Supplementary Table 2 workbook, including the columns that
#                               exist only there (comparison_label, the settings table, ...)
#
# Requires release_validation.R to be sourced first (statistical_model reads
# RELEASE_PRIMARY_MODEL). Holds no data.

RELEASE_DESCRIPTIONS <- c(
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
  original_sample_class = "PRE-CORRECTION sample class, normalized to the canonical vocabulary from the sample_annotation.xlsx group_label. Independently reproduced for all 96 measurements by the autosampler plate layout, and for the 6 corrected measurements by the retained UMAP correction table. Retained so the original identity is never lost.",
  analysis_sample_class = "ANALYSIS-TIME sample class. Identical to sample_class; named so the pairing with original_sample_class is clear.",
  sample_class_plate_layout_implied = "Sample class implied by the autosampler well row letter (A/E cfos, B/F neuropil, C/G mcherry, D/H neuron). This is a PRE-CORRECTION record and equals original_sample_class for all 96 measurements.",
  sample_class_corrected = "TRUE for the 6 measurements whose sample-class metadata were corrected during historical sample-identity QC, FALSE for the other 90. TRUE does NOT mean the sample was mislabeled or that any measurement changed; only the class metadata were reassigned.",
  sample_class_correction_status = "confirmed_intentional_switch_correction for the 6 corrected measurements, not_corrected for the other 90.",
  sample_class_correction_method = "How the correction was established, for corrected rows: historical_sample_identity_qc_umap_correction_table. Blank where nothing was corrected.",
  sample_class_correction_provenance = "Project-relative path of the preserved historical correction record, for corrected rows. Blank where nothing was corrected.",
  sample_class_correction_reference_sha256 = "SHA256 of the preserved correction record, recomputed and verified at build time against the value the audit recorded. Blank where nothing was corrected.",
  sample_class_correction_date = "Date of the applied metadata edit, for corrected rows. Blank where nothing was corrected.",
  plate_position = "Autosampler position as plate-set token and well, e.g. S5-B6.",
  correction_status = "Verdict on the sample-class correction: confirmed_intentional_switch_correction. A deliberate correction applied during historical sample-identity QC, not an unresolved discrepancy.",
  correction_method = "How the correction was established: historical_sample_identity_qc_umap_correction_table -- identified from UMAP nearest-class-center distances, recorded in a retained project correction table, and subsequently supported by independent proteomic profile-similarity and left/right-pair analyses.",
  correction_date = "Date of the applied metadata edit.",
  correction_reference = "Project-relative path of the preserved historical correction record. This table is the preserved record; no surviving prose note explains the rationale.",
  correction_reference_sha256 = "SHA256 of the correction record, recomputed at build time and required to match the value the audit recorded.",
  correction_reference_mtime = "Modification time of the correction record. Reported for context; the SHA256, not the mtime, is what the build verifies.",
  quantitative_values_changed = "Always FALSE. The quantitative matrix was compared across the metadata edit and the maximum absolute numeric difference was 0. Only sample-class metadata were reassigned.",
  umap_nearest_class_suggestion = "The class the UMAP correction table recorded as the nearest class center. Published so the divergence from the applied correction is visible: it agrees for 5 of the 6, and differs for N60 (nearest center cfos, applied neuron). UMAP identified the switch; it is not treated as ground truth for what the switch was.",
  original_class_corroborating_records = "The independent pre-correction records that agree on original_sample_class for this measurement.",
  animal_level_unit_affected = "The animal-level unit this measurement was averaged into. These units were built from the corrected assignments in the first place, so the validated results already reflect them.",
  condition_code = "Numeric condition code: 1 = paired_cno, 2 = paired_veh, 3 = unpaired_cno, 4 = unpaired_veh.",
  condition = "Condition label: paired_cno, paired_veh, unpaired_cno or unpaired_veh.",
  pairing_status = "paired or unpaired, decomposed from condition.",
  treatment = "cno or veh, decomposed from condition. Tokens are the project's literal condition coding.",
  collection_plate = paste("Plate1 or Plate2: the plate the sample was collected onto.", release_plate_injection_sentence()),
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
  aggregation_policy = "How the two hemispheres were combined: equal-weight mean of the Left and Right values of the measurement-level matrix, with no further transformation. The value is the policy name recorded when the animal-level matrix was built; the values it averaged are standardized per protein, not log2 abundance.",
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
  effect_size_sd_units = "PUBLIC EFFECT SIZE. Difference of group means (numerator minus denominator) on the protein-standardized abundance scale -- a standardized abundance difference expressed in SD units of that scale. NOT a log2 fold change, despite being stored as logFC in the statistical GCT (limma's coefficient name) and log2fc in the canonical split tables. Positive means higher in the numerator condition. A value of 2 means two standard deviations on the standardized scale, not a four-fold change. Not called a standardized mean difference unqualified, because that phrasing reads as Cohen's d (a pooled within-group SD) rather than the across-dataset per-protein SD actually used. The exported value is bit-identical to the source logFC; no transformation, rescaling or recomputation was applied.",
  effect_size_units = "Units of effect_size_sd_units: per-protein standard deviation of the standardized abundance scale.",
  effect_size_definition = "Full statement of what the effect size is and is not. Carried once per comparison in the summary table rather than repeated on every protein row.",
  source_statistic_field = "Bare name of the canonical source statistic the effect size was copied from: logFC. Machine-readable provenance; the exported effect_size_sd_units equals this source field exactly.",
  effect_size_source_column = "The canonical column the effect size was copied from, named so the published value can be traced: logFC (statistical GCT, limma) / log2fc (canonical split tables).",
  average_standardized_abundance = "limma AveExpr: the protein's mean standardized abundance over the 12 animal-level units of this comparison's sample class (the columns of that class's model), so it is the same for the three comparisons of a sample class. It is not centered at 0 within a sample class; averaged over the four sample classes it is close to 0 for every protein (bound checked by 13_validate_release.R), as a per-protein standardized input implies.",
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
  NES = "GSEA normalized enrichment score. Positive means enriched toward the numerator condition.",
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
  rank_statistic = "Statistic the ranked list was ordered by, using the token the canonical run recorded: moderated_t (canonical), log2fc (sensitivity) or not_applicable (ORA). The log2fc token is retained for provenance; the statistic it names is the standardized-abundance effect size, NOT a log2 fold change.",
  rank_source_column = "Column of the differential table used as the ranking statistic.",
  analysis_role = "canonical, sensitivity, or alternative_query_list. Canonical and sensitivity rows must never be pooled.",
  gene_identifier_type = paste("Identifier namespace of the gene members in core_enrichment:",
                               "UNIPROT (UniProt accession) for GO gene sets, ENTREZID (NCBI Gene",
                               "ID) for KEGG gene sets."),
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
  gsea_seed = "GSEA blocks only: the random seed of the fgsea permutations for this comparison, as the enrichment run recorded it (one seed per comparison, derived from gsea_seed_base in enrichment/enrichment_run_parameters.tsv). Blank for ORA.",
  gsea_untested_gene_sets = "GSEA blocks only: gene sets in the size limits that fgsea gave no P value (unbalanced statistics), so they are absent from the block. GO sets are listed one per row in enrichment/gsea_untested_gene_sets.tsv; KEGG is counted from its run's warnings.",
  n_perm_simple = "fgsea nPermSimple of the run that could not test the set.",
  parameter = "Parameter name.",
  value = "Parameter value.",
  stage = "Analysis stage the parameter belongs to.",
  recorded_by = "Artifact that records the value.",
  evidence_path = "Path to the artifact that records the value. Files in the code repository are named relative to its proteomic_analysis/ folder.",
  evidence_sha256 = "SHA256 of that artifact.",
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
  panel_source_data_canonical_path = "Canonical file the published panel source data were built from.",
  panel_source_data_sha256 = "SHA256 of the canonical panel source data; the published CSV has its own hash in the release manifest.",
  interpretation_note = "What a reader must know to read this panel correctly.",
  # workbook
  sheet = "Sheet (tab) name as it appears in the workbook.",
  sheet_key = "Stable internal key of the sheet; the builders and checks refer to sheets by this key.",
  supplementary_table = paste("Sub-table letter of the sheet within Supplementary Table 2 (2a,",
                              "2b, ...). Blank for the README and Column definitions sheets."),
  title = "Sheet title as printed in the workbook and listed on the Contents sheet.",
  n_rows = "Rows on the sheet.",
  n_columns = "Columns on the sheet.",
  columns = paste("Semicolon-separated MACHINE column names on the sheet -- the naming the",
                  "tests and this dictionary key off, matching the .tsv siblings."),
  display_columns = paste("Semicolon-separated DISPLAY headers actually printed on the",
                          "sheet, in the same order as `columns`. The workbook is the",
                          "reader-facing artifact and carries these; the .tsv files carry",
                          "the machine names."),
  table_level_constants = paste("Columns whose value was identical on every row and were",
                                "therefore moved into the line above the header, as",
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
  artifact_id = "Stable identifier for the artifact within the lineage.",
  artifact_type = "What kind of artifact it is.",
  path_or_release_filename = "Where the artifact lives, or its name in this release. Files in the code repository are named relative to its proteomic_analysis/ folder (for example 07_publication_release/metadata/experimenter_metadata.tsv).",
  parent_artifact = "Semicolon-separated artifact_id(s) this artifact was derived from.",
  processing_stage = "The operation that produced this artifact. UNRESOLVED when it cannot be proven.",
  software = "Software that performed the step.",
  software_version = "Version of that software, or UNKNOWN.",
  sha256 = "SHA256 of the artifact, where a single file could be hashed.",
  canonical = "Whether this artifact is the canonical one for its stage.",
  # manifest
  relative_path = "Path of the file relative to the release root.",
  file_role = "What the file is for.",
  format = "File format.",
  rows_if_tabular = "Row count for tabular files.",
  columns_if_tabular = "Column count for tabular files.",
  source_artifact = "Semicolon-separated canonical artifact(s) the file was derived from. A relative path that starts with a release folder (metadata/, provenance/ and so on) names a file in this release; one that starts with a repository folder (for example 07_publication_release/, R/ or tests/) names a file in the proteomic_analysis/ folder of the code repository.",
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
  pass = "Whether the check passed.",
  # significant gene sets
  core_enrichment_gene_symbols = paste(
    "GSEA leading-edge members as gene symbols, in the order and with the separator of",
    "core_enrichment. GO members are mapped from UniProt accession and KEGG members from",
    "NCBI Gene ID (through the KEGG run's own Entrez-to-UniProt record for the same",
    "comparison) to gene_symbol of processed_data/protein_feature_annotation.tsv.gz; a member",
    "without a gene symbol keeps its own identifier (UniProt accession for GO, NCBI Gene ID",
    "for KEGG).")
)

# --------------------------------------------------------------------------------------
# reader-facing definitions for the workbook
# --------------------------------------------------------------------------------------
# One entry per machine column shown in the workbook: c(definition, units). Units are blank
# where the quantity has none. RELEASE_READER_DEFINITIONS_BY_SHEET overrides an entry for one
# sheet where the same column means something narrower there. Stage 07 stops if a shown
# column has no entry, so a column cannot reach the workbook undefined.

release_rdef <- function(definition, units = "") c(definition = definition, units = units)

# The comparison name alone, for the sheets with no effect size, NES or direction.
RELEASE_COMPARISON_NAME_DEFINITION <- release_rdef(
  "The comparison, named '<compartment>: <group 1> vs <group 2>'.")

RELEASE_GSEA_DIRECTION_DEFINITION <- release_rdef(paste(
  "Group toward which the gene set is enriched: group 1 for a positive NES, group 2 for a",
  "negative NES."))

RELEASE_READER_DEFINITIONS <- list(
  # comparisons
  comparison_label = release_rdef(paste(
    "The comparison, named '<compartment>: <group 1> vs <group 2>'. Effect sizes, NES and",
    "directions compare group 1 with group 2.")),
  canonical_comparison = release_rdef(paste(
    "Machine-readable comparison key, '<compartment>_<group 1>_over_<compartment>_<group 2>'.",
    "Joins the comparison tables of this workbook and the files of the data release.")),
  contrast_family = release_rdef(paste(
    "Type of comparison: Pairing (paired VEH vs unpaired VEH), CNO in paired animals (paired",
    "CNO vs paired VEH) or CNO in unpaired animals (unpaired CNO vs unpaired VEH).")),
  numerator_condition = release_rdef(paste(
    "Group 1, the first-named group. A positive effect size means higher in group 1.")),
  denominator_condition = release_rdef("Group 2, the reference group."),
  n_numerator_animals = release_rdef("Animals in group 1.", "animals"),
  n_denominator_animals = release_rdef("Animals in group 2.", "animals"),
  n_proteins_tested = release_rdef("Protein groups tested.", "protein groups"),
  n_significant_fdr_0_05 = release_rdef(
    "Protein groups at FDR < 0.05 (Benjamini-Hochberg adjusted P < 0.05).", "protein groups"),
  n_significant_higher_in_numerator = release_rdef(
    "Protein groups at FDR < 0.05 with a positive effect size (higher in group 1).",
    "protein groups"),
  n_significant_higher_in_denominator = release_rdef(
    "Protein groups at FDR < 0.05 with a negative effect size (lower in group 1).",
    "protein groups"),
  n_proteins_tested_mapped_only = release_rdef(paste(
    "Annotated protein groups tested: those with a mouse UniProt accession. The other",
    "groups are non-mouse entries such as trypsin or keratins."), "protein groups"),
  n_significant_fdr_0_05_mapped_only = release_rdef(
    "Annotated protein groups at FDR < 0.05 (the proteins of table 2d).", "protein groups"),
  n_significant_higher_in_numerator_mapped_only = release_rdef(
    "Annotated protein groups at FDR < 0.05, higher in group 1.", "protein groups"),
  n_significant_higher_in_denominator_mapped_only = release_rdef(
    "Annotated protein groups at FDR < 0.05, lower in group 1.", "protein groups"),

  # samples and animals
  sample_id = release_rdef(paste(
    "Identifier of the LC-MS run. Also the column name in the measurement-level abundance",
    "matrix of the data release.")),
  raw_file_basename = release_rdef("Name of the raw acquisition file (a .d directory) of the run."),
  AnimalID = release_rdef("Animal identifier (12 animals)."),
  hemisphere = release_rdef("Hemisphere the material came from: Left or Right."),
  sample_class = release_rdef(paste(
    "Marker-defined compartment collected by laser-capture microdissection: mCherry+",
    "(mCherry-positive material), c-Fos+ (c-Fos-positive material), NeuN+ (NeuN-positive",
    "cells negative for c-Fos and mCherry) or neuropil (areas negative for c-Fos, mCherry",
    "and NeuN).")),
  condition = release_rdef(
    "Group of the animal: paired VEH, paired CNO, unpaired VEH or unpaired CNO."),
  pairing_status = release_rdef("paired or unpaired."),
  treatment = release_rdef("CNO or VEH (vehicle)."),
  collection_plate = release_rdef(paste(
    "Plate the sample was collected onto (Plate1 or Plate2).",
    release_plate_injection_sentence())),
  well_position = release_rdef("Well of the run on the autosampler plate."),
  injection_index = release_rdef("Injection order of the run, 1 to 96."),
  acquisition_date = release_rdef("Date of the LC-MS run."),
  animal_level_sample_id = release_rdef(paste(
    "Animal-level profile the run was averaged into, '<animal>_<compartment>'; joins table",
    "2a with table 2b.")),
  original_sample_class = release_rdef(paste(
    "Compartment in the original sample records. It differs from Compartment for the runs",
    "whose label was corrected.")),
  sample_class_corrected = release_rdef(paste(
    "TRUE for the runs whose compartment label was corrected during sample-identity quality",
    "control. Only the label changed; the measured values did not.")),
  animal_level_column_name = release_rdef(paste(
    "Animal-level profile, '<animal>_<compartment>'. Also the column name in the",
    "animal-level abundance matrix of the data release.")),
  left_source_sample = release_rdef("Sample ID (table 2a) of the left-hemisphere run."),
  right_source_sample = release_rdef("Sample ID (table 2a) of the right-hemisphere run."),
  n_hemisphere_measurements = release_rdef("Hemisphere runs averaged into the profile.", "runs"),

  # proteins
  protein_group_id = release_rdef(paste(
    "Protein group: UniProt entry name(s), separated by ';' when the group holds several",
    "proteins. Joins tables 2d and 2e with the files of the data release.")),
  uniprot_accession = release_rdef("UniProt accession of the leading protein of the group."),
  gene_symbol = release_rdef("Gene symbol of the leading protein. Blank where none is known."),
  protein_description = release_rdef("Protein name of the leading protein."),
  id_mapping_status = release_rdef(paste(
    "TRUE when the group has a mouse UniProt accession (an annotated protein, as counted in",
    "tables 2c and 2d); FALSE for a non-mouse entry such as trypsin or a keratin.")),
  effect_size_sd_units = release_rdef(paste(
    "Difference in standardized abundance, group 1 minus group 2. Abundances were",
    "standardized separately for each protein, so the unit is one standard deviation of that",
    "protein. Not a log2 fold change: 2 means two standard deviations, not four-fold."),
    "SD units"),
  moderated_t = release_rdef(
    "limma moderated t statistic. The ranking statistic of the primary GSEA."),
  P.Value = release_rdef("P value of the moderated t test."),
  adj.P.Val = release_rdef(paste(
    "Benjamini-Hochberg adjusted P value (FDR) within the comparison. Values below 0.05 are",
    "called significant.")),
  significant_fdr_0_05 = release_rdef("TRUE when the adjusted P value is below 0.05."),
  average_standardized_abundance = release_rdef(paste(
    "Mean standardized abundance of the protein over the 12 animal-level profiles of the",
    "compartment (limma AveExpr)."), "SD units"),
  B_log_odds = release_rdef(
    "limma B statistic: log-odds that the protein is differentially abundant.", "log odds"),
  direction = release_rdef(
    "Group in which the protein is higher: group 1 for a positive effect size."),

  # gene sets
  ontology = release_rdef(paste(
    "Gene-set collection: GO-BP (GO biological process), GO-CC (GO cellular component) or",
    "KEGG (KEGG pathways).")),
  term_id = release_rdef("Gene-set identifier: GO term ID or KEGG pathway ID."),
  term_name = release_rdef("Gene-set name."),
  NES = release_rdef(paste(
    "GSEA normalized enrichment score. Positive: the set is enriched among proteins higher in",
    "group 1; negative: among proteins higher in group 2.")),
  enrichment_score = release_rdef("GSEA enrichment score before normalization."),
  set_size = release_rdef("Measured proteins in the gene set.", "proteins"),
  n_core_enrichment_genes = release_rdef(
    "Proteins in the leading-edge (core enrichment) subset of the gene set.", "proteins"),
  core_enrichment_gene_symbols = release_rdef(paste(
    "Leading-edge proteins as gene symbols, separated by '/', in the order GSEA listed them.",
    "A protein without a gene symbol is given by its identifier (UniProt accession for GO,",
    "NCBI Gene ID for KEGG).")),
  rank_at_max = release_rdef(paste(
    "Position in the ranked protein list at which the running enrichment score reaches its",
    "maximum."), "rank"),
  p_value = release_rdef(paste(
    "GSEA P value from fgsea permutations. fgsea does not report values below 1E-10; such",
    "values are shown as <=1E-10.")),
  adjusted_p_value = release_rdef(paste(
    "Benjamini-Hochberg adjusted P value (FDR) within the comparison and gene-set",
    "collection.")),
  q_value = release_rdef("Storey q value as reported by clusterProfiler."),
  query_list = release_rdef(paste(
    "Proteins tested for over-representation: all proteins at FDR < 0.05, those higher in",
    "group 1 or those lower in group 1.")),
  gene_count = release_rdef("Query proteins annotated to the gene set.", "proteins"),
  gene_ratio = release_rdef(paste(
    "Query proteins in the gene set / query proteins with any annotation in the",
    "collection.")),
  background_ratio = release_rdef(paste(
    "Background proteins in the gene set / background proteins with any annotation in the",
    "collection.")),
  fold_enrichment = release_rdef("Gene ratio divided by background ratio."),
  rich_factor = release_rdef(
    "Query proteins in the gene set divided by all background proteins in the gene set."),
  z_score = release_rdef(paste(
    "Difference between the observed and the expected number of query proteins in the gene",
    "set, in standard deviations (clusterProfiler zScore).")),
  rank_statistic = release_rdef(paste(
    "Statistic the proteins were ranked by: moderated t (primary GSEA) or effect size (SD",
    "units; the effect-size-ranked sensitivity GSEA).")),
  analysis_role = release_rdef(paste(
    "primary: the GSEA reported in the paper (tables 2g to 2i); sensitivity: the",
    "effect-size-ranked GSEA, which is in the data release only.")),

  # enrichment settings and coverage (table 2m only)
  analysis = release_rdef(
    "GSEA (gene set enrichment analysis) or ORA (over-representation analysis)."),
  gene_set_source = release_rdef("Gene-set collection, its source and version."),
  set_size_limits = release_rdef(
    "Gene sets with fewer or more measured proteins were not tested.", "proteins"),
  permutations = release_rdef("fgsea permutation setting (nPermSimple)."),
  seed = release_rdef(
    "Random seed of the fgsea permutations for this comparison, as the run recorded it."),
  ora_background = release_rdef(paste(
    "ORA background: measured proteins with a mouse UniProt accession and an annotation in",
    "the collection."), "proteins"),
  multiple_testing = release_rdef("Multiple-testing adjustment of the P values."),
  rows_exported = release_rdef(paste(
    "Gene sets with a result for this comparison and collection: the rows of tables 2g to",
    "2k."), "gene sets"),
  gsea_untested_gene_sets = release_rdef(paste(
    "Gene sets within the size limits that fgsea could not give a P value (strongly",
    "one-sided ranking). They are absent from the GSEA tables; GO sets are listed in table",
    "2l."), "gene sets"),
  note = release_rdef("Note on the row."),

  # EWCE
  ewce_analysis_type = release_rdef(paste(
    "Baseline: the most abundant proteins of one condition; Differential: the proteins with",
    "the largest increase in one group of a comparison, from EWCE's own limma fit.")),
  protein_list_definition = release_rdef("Protein list tested for cell-type enrichment."),
  top_n = release_rdef("Proteins in the list.", "proteins"),
  annotation_level = release_rdef(paste(
    "Cell-type level of the reference: 1, broad classes; 2, subclusters. The primary setting",
    "is level 2.")),
  cell_type = release_rdef(paste(
    "Cell type of the ewceData mouse cortex and hippocampus single-cell reference.",
    "'unassigned level-2 cluster' is the reference's cluster that has no level-2 name.")),
  fold_change = release_rdef(paste(
    "EWCE fold change: mean expression proportion of the cell type in the list divided by",
    "its mean in the bootstrap null lists.")),
  sd_from_mean = release_rdef(
    "Distance of the observed value from the bootstrap null mean, in standard deviations.",
    "SD"),
  adjusted_p_value_within_target = release_rdef(
    "Benjamini-Hochberg adjusted P value within the tests of one protein list."),
  adjusted_p_value_global = release_rdef(paste(
    "Benjamini-Hochberg adjusted P value over all EWCE tests; the significance criterion.",
    "0 is shown as <0.0001.")),
  significant_global_fdr_0_05 = release_rdef(
    "TRUE when the global adjusted P value is below 0.05."),
  ewce_direction = release_rdef(paste(
    "Sign of the cell-type result as EWCE reports it: Positive when the observed value is",
    "above the bootstrap null mean, Negative when below.")),
  is_primary_setting = release_rdef(paste(
    "TRUE for the primary setting (top 250 proteins, annotation level 2); other rows are",
    "sensitivity settings.")),

  # secondary analyses (table 2o)
  analysis_name = release_rdef(paste(
    "Name of the secondary analysis and, where it has one, the figure panel it supports.")),
  purpose = release_rdef("Question the analysis answers."),
  model = release_rdef("Statistical model or method."),
  n_per_group = release_rdef("Sample size.", "animals"),
  results_location = release_rdef(
    "Table of this workbook or file of the data release that holds the results."),
  key_result = release_rdef(paste(
    "Main finding of the analysis, with the counts or ranges it rests on.")),

  # figure panels (table 2p)
  panel_label = release_rdef(paste(
    "Figure and panel letter: Fig. 2 is Figure 2 and Supp. Fig. 4 is Supplementary Figure",
    "4 of the paper.")),
  panel_title = release_rdef("What the panel shows."),
  statistical_unit = release_rdef("Unit each value in the panel is computed on."),
  panel_source = release_rdef(
    "Table of this workbook the panel is drawn from, and its Source Data file(s)."),
  rows_used = release_rdef("Rows of the source table the panel uses."),
  key_numbers = release_rdef(paste(
    "Counts and statistics the panel shows or rests on, such as the number of proteins or",
    "gene sets at FDR < 0.05 and the variance explained by a principal component.")),

  # software (table 2q)
  component = release_rdef("Software, package or reference dataset."),
  version = release_rdef("Version recorded by the analysis run."),
  used_for = release_rdef("Analysis it was used for.")
)

# Per-sheet meaning, keyed by sheet key (RELEASE_SHEET_META in release_presentation.R).
RELEASE_READER_DEFINITIONS_BY_SHEET <- list(
  Significant_Gene_Sets = list(direction = RELEASE_GSEA_DIRECTION_DEFINITION),
  GSEA_GO_BP = list(direction = RELEASE_GSEA_DIRECTION_DEFINITION),
  GSEA_GO_CC = list(direction = RELEASE_GSEA_DIRECTION_DEFINITION),
  GSEA_KEGG = list(direction = RELEASE_GSEA_DIRECTION_DEFINITION),
  ORA_GO_BP = list(
    p_value = release_rdef("Hypergeometric P value of over-representation."),
    adjusted_p_value = release_rdef(paste(
      "Benjamini-Hochberg adjusted P value (FDR) within the comparison, collection and query",
      "list."))),
  ORA_GO_CC = list(
    p_value = release_rdef("Hypergeometric P value of over-representation."),
    adjusted_p_value = release_rdef(paste(
      "Benjamini-Hochberg adjusted P value (FDR) within the comparison, collection and query",
      "list."))),
  GSEA_Untested_Sets = list(comparison_label = RELEASE_COMPARISON_NAME_DEFINITION),
  Enrichment_Settings = list(
    comparison_label = RELEASE_COMPARISON_NAME_DEFINITION,
    rank_statistic = release_rdef("Statistic the proteins were ranked by for GSEA; none for ORA."),
    query_list = release_rdef("ORA only: the proteins tested for over-representation.")),
  EWCE = list(
    p_value = release_rdef(paste(
      "Bootstrap P value (10,000 replicates). 0 means no bootstrap list reached the observed",
      "value and is shown as <0.0001.")),
    sample_class = release_rdef("Compartment the protein list comes from."),
    comparison_label = release_rdef(paste(
      "Comparison the differential protein list comes from, named '<compartment>: <group",
      "1> vs <group 2>'. Blank on Baseline rows, which have no comparison."))),
  Figure_Panels = list(
    note = release_rdef(paste(
      "What a reader needs to read the panel correctly: how the plotted rows were chosen,",
      "caveats of the design and gene sets fgsea could not test.")))
)

#' Reader-facing definition and unit of one workbook column on one sheet, or NULL.
release_reader_definition <- function(column, sheet_key) {
  by_sheet <- RELEASE_READER_DEFINITIONS_BY_SHEET[[sheet_key]]
  if (!is.null(by_sheet) && column %in% names(by_sheet)) return(by_sheet[[column]])
  if (column %in% names(RELEASE_READER_DEFINITIONS)) return(RELEASE_READER_DEFINITIONS[[column]])
  NULL
}
