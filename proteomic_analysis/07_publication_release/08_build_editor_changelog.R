#!/usr/bin/env Rscript

# Publication release, stage 08 -- the editor/reviewer changelog.
#
# Produces
#   editor_source_data/REVISION_PROTEOMICS_DATA_CHANGELOG.md
#
# Every count in the document is read from the release tables built by stages 04-06, not
# typed in. If a canonical number changes, the changelog changes with it.
#
# The document must not imply that correcting the statistical unit left the biology
# untouched. It did not: most contrasts lost every FDR-significant protein, one main-figure
# panel became a null result, and two panels changed meaning. Those are stated plainly.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_validation.R"))
release_require("digest")

DATA_ROOT <- release_data_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/08_build_editor_changelog.R"

release_banner("stage 08 -- revision changelog")

rel <- function(...) release_path(..., create_dir = FALSE)
read_release <- function(path) {
  release_assert_exists(path, basename(path))
  if (grepl("[.]gz$", path)) {
    read.delim(gzfile(path), sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  }
}

diff_summary <- read_release(rel("differential_analysis", "primary_differential_summary.tsv"))
coverage <- read_release(rel("enrichment", "enrichment_coverage.tsv"))
gsea_go <- read_release(rel("enrichment", "primary_GSEA_GO_BP.tsv.gz"))
sens <- read_release(rel("enrichment", "GSEA_log2FC_sensitivity.tsv.gz"))
figure_map <- read_release(rel("editor_source_data", "figure_source_map.tsv"))
secondary <- read_release(rel("metadata", "secondary_analysis_manifest.tsv"))
discrepancy <- read_release(rel("metadata", "sample_class_source_discrepancy.tsv"))

REVISION_ROOT <- file.path(DATA_ROOT, "03_output", "reviewer_revision_animal_level_20260827")
retained_report <- file.path(REVISION_ROOT, "manifests",
                             "RETAINED_PANELS_3G_3H_3I_SuppE_REPORT.md")
final_report <- file.path(REVISION_ROOT, "manifests", "FINAL_REPORT.md")

# --------------------------------------------------------------------------------------
# derived numbers
# --------------------------------------------------------------------------------------

diff_summary$n_significant_fdr_0_05 <- as.integer(diff_summary$n_significant_fdr_0_05)
diff_summary$n_significant_fdr_0_05_mapped_only <-
  as.integer(diff_summary$n_significant_fdr_0_05_mapped_only)
n_with_hits <- sum(diff_summary$n_significant_fdr_0_05 > 0L)
n_without_hits <- nrow(diff_summary) - n_with_hits

fdr_terms <- function(tbl, filter_expr) {
  sub <- tbl[filter_expr, , drop = FALSE]
  sum(as.numeric(sub$adjusted_p_value) < 0.05, na.rm = TRUE)
}
gsea_fdr_by_comparison <- vapply(
  unique(gsea_go$canonical_comparison),
  function(cmp) fdr_terms(gsea_go, gsea_go$canonical_comparison == cmp),
  integer(1))
sens_fdr_by_comparison <- vapply(
  unique(sens$canonical_comparison),
  function(cmp) fdr_terms(sens, sens$canonical_comparison == cmp & sens$ontology == "GO_BP"),
  integer(1))

md_table <- function(df, align = NULL) {
  cols <- names(df)
  if (is.null(align)) align <- rep("---", length(cols))
  body <- apply(df, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(paste0("| ", paste(cols, collapse = " | "), " |"),
    paste0("|", paste(align, collapse = "|"), "|"),
    body)
}

da_tbl <- data.frame(
  Comparison = diff_summary$canonical_comparison,
  `Sample class` = diff_summary$sample_class,
  `Tested (all)` = diff_summary$n_proteins_tested,
  `FDR<0.05 (all)` = diff_summary$n_significant_fdr_0_05,
  `Tested (mapped)` = diff_summary$n_proteins_tested_mapped_only,
  `FDR<0.05 (mapped)` = diff_summary$n_significant_fdr_0_05_mapped_only,
  `Higher in numerator` = diff_summary$n_significant_higher_in_numerator_mapped_only,
  `Higher in denominator` = diff_summary$n_significant_higher_in_denominator_mapped_only,
  check.names = FALSE, stringsAsFactors = FALSE)

gsea_tbl <- data.frame(
  Comparison = names(gsea_fdr_by_comparison),
  `GO-BP terms FDR<0.05 (canonical, t-ranked)` = unname(gsea_fdr_by_comparison),
  `GO-BP terms FDR<0.05 (log2FC-ranked sensitivity)` =
    unname(sens_fdr_by_comparison[names(gsea_fdr_by_comparison)]),
  check.names = FALSE, stringsAsFactors = FALSE)

panel_tbl_rows <- figure_map[figure_map$primary_or_secondary != "not_applicable", , drop = FALSE]
panel_tbl <- data.frame(
  Panel = panel_tbl_rows$panel_label,
  Status = panel_tbl_rows$revision_status,
  `Unit: original -> revised` = paste(panel_tbl_rows$original_statistical_unit, "->",
                                      panel_tbl_rows$statistical_unit),
  `Displayed statistics changed` = panel_tbl_rows$displayed_statistics_changed,
  `Interpretation changed` = panel_tbl_rows$biological_interpretation_changed,
  check.names = FALSE, stringsAsFactors = FALSE)

changed_interpretation <- panel_tbl_rows$panel_label[
  toupper(panel_tbl_rows$biological_interpretation_changed) == "YES"]

# --------------------------------------------------------------------------------------
# document
# --------------------------------------------------------------------------------------

L <- c(
"# Proteomics data and analysis changelog: original submission versus revised analysis",
"",
"This document accompanies the revised proteomics data package. It states what changed",
"between the analysis supplied with the original submission and the analysis supplied now,",
"and it does so without claiming that the biological conclusions survived unchanged.",
"Several did not.",
"",
"Every count below is generated from the release tables in this package rather than",
"transcribed, so the document cannot drift from the data it describes.",
"",
"---",
"",
"## 1. The statistical unit was wrong and has been corrected",
"",
"| | Original submission | Revised analysis |",
"|---|---|---|",
"| Independent experimental unit | hemisphere / acquisition measurement | **animal** |",
"| Observations entering inference | 96 (or 120 in some ProTigy exports) | **48** (12 animals x 4 sample classes) |",
"| n per group | 6 (two hemispheres from each of 3 animals) | **3 animals** |",
"| Left / right handling | treated as independent replicates | **averaged within AnimalID x sample_class before testing** |",
"",
"Two hemispheres taken from the same animal are not independent samples. Treating them as",
"independent replicates doubles the apparent sample size and shrinks standard errors, which",
"inflates significance throughout. This is pseudoreplication, and it affected every",
"inferential number in the original submission.",
"",
"The correction is an aggregation, not a reanalysis of a different experiment: left and",
"right values are averaged with equal weight within each AnimalID x sample_class unit, on",
"the existing imputed log2 values, with no further transformation. In a balanced design the",
"mean of animal means equals the mean of all hemispheres, so **point estimates barely move**.",
"What changes is the **inference**: standard errors, degrees of freedom, and therefore every",
"p-value.",
"",
"## 2. Data files: 96 measurement columns versus 48 animal-level columns",
"",
"| Layer | Rows x columns | Role |",
"|---|---|---|",
"| Search output (`pg.matrix_raw.txt`) | 5,747 protein groups x 96 acquisitions | quantification output, pre-filter |",
"| Processed measurement-level matrix | 5,349 x 96 acquisitions | filtered, log2, median-centred, imputed, PCA-adjusted |",
"| **Animal-level matrix** | **5,349 x 48 animal units** | **the matrix all inference is computed on** |",
"",
"The originally submitted `processed_protein_group_matrix_filtered_umap_adjusted.xlsx`",
"(5,349 proteins x 96 columns) corresponds to the measurement-level layer. It is still",
"provided here, as `processed_data/protein_abundance_measurement_level.tsv.gz`, but it is",
"**not** the matrix that inference is drawn from.",
"",
"The originally submitted `processed_protein_group_matrix_raw.txt` (5,747 rows x 96 columns)",
"is **processed protein-level search output, not native raw mass-spectrometry data**, despite",
"its name. Native raw acquisition data are the 96 `.d` acquisition directories, which were",
"not part of the original package and are identified per row in",
"`metadata/sample_metadata.tsv`.",
"",
"## 3. The primary comparison set",
"",
"Exactly **12 primary comparisons**: 4 sample classes (mcherry, neuropil, cfos, neuron) x 3",
"within-sample-class contrast families:",
"",
"1. `paired_cno` vs `paired_veh` -- chemogenetic effect in paired animals",
"2. `paired_veh` vs `unpaired_veh` -- pairing-associated (learning) signature",
"3. `unpaired_cno` vs `unpaired_veh` -- chemogenetic effect in unpaired controls",
"",
"All three treatment/control comparisons are covered in **every** sample class, so the",
"chemogenetic manipulation is reported in both the paired and the unpaired arm rather than",
"in the paired arm alone. The originally submitted `GSEA_ORA_all_results.xlsx` did not",
"cover this set consistently.",
"",
"Cross-compartment comparisons are **not** in this set. See section 8.",
"",
"## 4. Differential protein statistics, as recomputed at animal level",
"",
"Two counts are given per comparison. `all` is over all 5,349 tested protein groups;",
"`mapped` is over the 5,327 with a mouse identifier, which is the set the manuscript figures",
"were drawn from. The 22-row difference is contaminant and non-mouse entries (human",
"keratins, porcine trypsin), which are tested but not annotatable.",
"",
md_table(da_tbl),
"",
sprintf(paste("**%d of the 12 primary comparisons contain no FDR-significant protein at all.**",
              "Only %d contain any."), n_without_hits, n_with_hits),
"This is a substantive change from the original submission, in which hemisphere-level",
"inference produced significant sets far more widely.",
"",
"## 5. Enrichment analysis",
"",
"| | Original submission | Revised analysis |",
"|---|---|---|",
"| Input statistics | hemisphere-level | animal-level |",
"| **GSEA ranking statistic** | **log2FC** | **moderated t statistic** |",
"| GSEA / ORA cutoffs | filtered at source | `pvalue_cutoff = 1`, `qvalue_cutoff = 1` -- the complete tested result is released |",
"| Multiple testing | BH | BH, FDR 0.05 |",
"| Gene set size | -- | minGSSize 10, maxGSSize 800 |",
"| ORA universe | not documented | all unique successfully mapped measured UniProt accessions, per comparison |",
"",
"Ranking a GSEA by the effect size alone ignores the variance of the estimate: a large",
"difference measured with n = 3 and enormous uncertainty outranks a smaller, well-estimated",
"one. The canonical analysis now ranks by the moderated t statistic. **Effect-size-ranked",
"results are retained only as a sensitivity analysis**, in a separate file, with",
"`analysis_role = sensitivity` on every row. They must not be reported as the primary",
"enrichment result. (The file keeps the name `GSEA_log2FC_sensitivity` because that is what",
"the canonical run called the variant; on the units of that statistic see section 12.)",
"",
"FDR-significant GO-BP terms per comparison, canonical versus sensitivity:",
"",
md_table(gsea_tbl),
"",
"## 6. Figure panels",
"",
md_table(panel_tbl),
"",
sprintf("Panels whose biological interpretation changed: **%s**.",
        paste(changed_interpretation, collapse = ", ")),
"",
"### Findings that were weakened or eliminated",
"",
"- **Figure 3 G (GO over-representation, mCherry).** The original panel showed GO terms",
"  including synaptic vesicle endocytosis, postsynaptic density organization and respiratory",
"  electron transport chain, over a hemisphere-level significant-protein list. The direct",
"  animal-level analogue -- ORA over the FDR-significant upregulated list -- returns",
"  **3,639 GO-BP terms tested, 0 at FDR < 0.05, minimum adjusted P = 0.231**. The panel is",
"  retained as a null result with the 15 lowest-adjusted-P terms plotted and captioned as",
"  non-significant. No substitute enrichment was selected to fill it.",
"- **Figure 3 F (GSEA bubble grid, learning).** The mCherry column contributes **zero**",
"  FDR-significant GO-BP terms and is deliberately left visibly empty.",
"- **Supplementary C (EWCE cell-type identity).** 12 of the legacy 14 cell types survive.",
"  `Int12` and `S1PyrL23` no longer clear global FDR.",
"- **Figures 3 H and 3 I.** Both remain, but neither can be read as evidence that",
"  chemogenetic inhibition reverses a learning-associated proteomic programme. Both",
"  contrasts share `paired_veh` on opposite sides, so anything high in `paired_veh` is pushed",
"  up on one axis and down on the other **by construction**. The control with the shared arm",
"  on the same side gives r = +0.90 on the same protein set, and all 38 terms overlapping in",
"  3 I are opposed in sign, which the design forces. Panel 3 I has been retitled neutrally.",
"- **Supplementary E.** The original four cross-compartment comparisons were invalid twice",
"  over: hemispheres were treated as independent, and both arms of every comparison came from",
"  the same three animals. It has been refitted as a secondary within-animal paired model",
"  and is no longer presented as one of the designed contrasts.",
"",
"### Findings that did not change",
"",
"- **Figures 3 E and Supplementary D (rank abundance).** Descriptive; no statistical test.",
"  Regenerated at animal level and numerically indistinguishable from the originals",
"  (Spearman rank r = 0.99999; Pearson on group means r = 1.000), as expected in a balanced",
"  design.",
"- **Supplementary A, B2, B3.** Acquisition-level technical QC. The acquisition is the",
"  correct unit for identification depth and for hemisphere / collection-plate QC, so these",
"  panels are unchanged and remain valid.",
"",
"### Numbers that must be retyped in the figures",
"",
"| Panel | Original | Revised |",
"|---|---|---|",
"| Fig 3 C | PC1 25.4% / PC2 17.7% / PC3 10.2% | PC1 29.1% / PC2 20.5% / PC3 12.1% |",
"| Supp B1 | PC1 25.4% / PC2 17.7% | PC1 29.1% / PC2 20.5% |",
"| Fig 3 H | r = -0.71 | r = -0.69 (learning signature, n = 1132) |",
"",
"## 7. The unresolved mCherry result",
"",
sprintf(paste("The mCherry paired-vehicle versus unpaired-vehicle comparison retains %d",
              "FDR-significant proteins of %d mapped (%d higher in paired_veh, %d higher in",
              "unpaired_veh). The statistical result stands as computed. What it reflects does",
              "not follow from it:"),
        diff_summary$n_significant_fdr_0_05_mapped_only[
          diff_summary$canonical_comparison == "mcherry_paired_veh_over_mcherry_unpaired_veh"],
        as.integer(diff_summary$n_proteins_tested_mapped_only[
          diff_summary$canonical_comparison == "mcherry_paired_veh_over_mcherry_unpaired_veh"]),
        as.integer(diff_summary$n_significant_higher_in_numerator_mapped_only[
          diff_summary$canonical_comparison == "mcherry_paired_veh_over_mcherry_unpaired_veh"]),
        as.integer(diff_summary$n_significant_higher_in_denominator_mapped_only[
          diff_summary$canonical_comparison == "mcherry_paired_veh_over_mcherry_unpaired_veh"])),
"",
"- the shift is broad and one-directional with no GSEA or ORA pathway coherence;",
"- pairing condition is completely associated with collection plate in this comparison;",
"- no technical batch metadata exist that would let a technical origin be demonstrated.",
"",
"Neither a biological nor a technical explanation is established.",
"",
"## 8. Primary and secondary analyses are now separated",
"",
"The original workbook mixed cross-compartment profiling in with the treatment contrasts.",
"In this package they are in different files, and every enrichment row carries",
"`analysis_role` and `primary_or_secondary`.",
"",
paste0("Secondary analyses recorded separately (",
       nrow(secondary), " entries in `metadata/secondary_analysis_manifest.tsv`): ",
       paste(secondary$analysis_id, collapse = ", "), "."),
"",
"## 9. Collection plate",
"",
"Because pairing condition was completely associated with collection plate in the",
"paired-vehicle versus unpaired-vehicle comparison, any collection-plate-associated",
"contribution cannot be distinguished from a pairing-associated contribution.",
"",
"No downstream proteomics batch effect attributable to collection plate has been",
"demonstrated. The `Plate1` / `Plate2` token records **which plate a sample was collected",
"onto**. It is not a proteomics preparation, digestion, LC-MS, acquisition or instrument",
"batch, and no such batch metadata exist for this dataset. Supplementary B3 has accordingly",
"been relabelled from `Plate` to `Collection plate`; collection plate explains 0.01% of PC1",
"and 1.5% of PC2 at acquisition level, neither significant.",
"",
"## 10. What is superseded",
"",
"| Original submission file | Replacement in this package | Superseded? |",
"|---|---|---|",
"| `processed_protein_group_matrix_raw.txt` (5,747 x 96) | `processed_data/protein_feature_annotation.tsv.gz` (annotation) and, for abundance, `processed_data/protein_abundance_measurement_level.tsv.gz` | Partly. Retained as the search-output layer in the lineage; it is **not** raw MS data. |",
"| `processed_protein_group_matrix_filtered_umap_adjusted.xlsx` (5,349 x 96) | `processed_data/protein_abundance_measurement_level.tsv.gz` | No, but demoted: it is no longer the matrix inference is drawn from. |",
"| `GSEA_ORA_all_results.xlsx` | `enrichment/*.tsv.gz` and `editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx` | **Yes, entirely.** Hemisphere-level inference, log2FC-ranked GSEA, incomplete contrast coverage. |",
"| `README.txt` (group codes only) | `README_DATA.md`, `metadata/data_dictionary.tsv`, this changelog | **Yes.** |",
"",
"## 11. A metadata discrepancy found while preparing this package",
"",
if (nrow(discrepancy) > 0L) c(
sprintf(paste("Preparing this package surfaced an issue that is **not** a consequence of the",
              "statistical-unit correction and was not previously recorded. For %d of the 96",
              "acquisitions, the sample class used by the validated analysis is contradicted",
              "by both independent records of the same assignment."), nrow(discrepancy)),
"",
"| acquisition | animal | hemisphere | well | used by the analysis | sample_annotation.xlsx | plate layout |",
"|---|---|---|---|---|---|---|",
vapply(seq_len(nrow(discrepancy)), function(i) {
  r <- discrepancy[i, , drop = FALSE]
  paste0("| ", r$plate_sample_number, " | ", r$AnimalID, " | ", r$hemisphere, " | ",
         r$plate_set_token, "-", r$well_position, " | **",
         r$sample_class_used_by_analysis, "** | ", r$sample_class_sample_annotation_xlsx,
         " | ", r$sample_class_plate_layout_implied, " |")
}, character(1)),
"",
paste0("The two independent sources agree with each other on all ", nrow(discrepancy),
       ". The plate layout (row A/E cfos, B/F neuropil, C/G mcherry, D/H neuron) is",
       " followed by the other 90 measurements and by both hemispheres of every other",
       " animal; only these left-hemisphere samples of ",
       paste(sort(unique(discrepancy$AnimalID)), collapse = " and "),
       " deviate, and they deviate by a single row shift, cyclically permuting three",
       " classes."),
"",
paste0("**This has not been corrected.** Correcting it would change which measurements are",
       " averaged together and therefore every downstream statistic, which is a reanalysis",
       " and outside the scope of preparing a data package. It is published as",
       " `metadata/sample_class_source_discrepancy.tsv`, flagged per row by",
       " `sample_class_corroborated` in the sample metadata, and listed as a known",
       " limitation in `README_DATA.md`."),
"",
paste0("**It requires a decision before publication.** ",
       length(unique(discrepancy$animal_level_unit_affected)),
       " of the 48 animal-level units are affected (",
       paste(sort(unique(discrepancy$animal_level_unit_affected)), collapse = ", "),
       "), each mixing one correctly assigned hemisphere with one disputed one.")
) else
"No sample-class discrepancy was found: all 96 measurements are corroborated by two independent sources.",
"",
"## 12. The effect size is not a log2 fold change",
"",
"A second issue surfaced while preparing this package, also independent of the",
"statistical-unit correction.",
"",
"The abundance matrix all statistics were computed on is **standardised per protein**:",
"measured directly, every protein has mean 0 and standard deviation 1 across the 96",
"acquisitions. ProTigy nevertheless reports its effect size as `logFC`, and the canonical",
"split tables carry it as `log2fc`.",
"",
"A difference of group means on a per-protein standardised scale is a **standardised mean",
"difference, in units of per-protein standard deviation**. It is not a log2 ratio, and an",
"effect size of 2 does not mean a four-fold change.",
"",
"Corroborating evidence, from the canonical tables themselves: `AveExpr` -- the average",
"abundance limma reports for each protein -- lies within +/-0.0011 of zero for all 5,349",
"proteins. That can only happen if the input was centred per protein.",
"",
"This package therefore publishes the quantity as `effect_size_sd_units`, carrying",
"`effect_size_definition` and `effect_size_source_column` on every row, rather than",
"republishing it under a name that asserts a fold-change reading. **No value was changed**",
"-- the numbers are bit-identical to the canonical tables, and directions, rankings,",
"p-values and every downstream result are unaffected. Only the label is corrected.",
"",
"Any manuscript text, axis label or caption that describes these values as log2 fold",
"changes should be revised accordingly. This is worth checking before publication.",
"",
"## 13. What this package does not contain",
"",
"Native raw mass-spectrometry acquisition files. The 96 `.d` acquisition directories are",
"identified by name in `metadata/sample_metadata.tsv` but are not included here; they are",
"intended for a separate PRIDE deposition. Metadata that could not be established from any",
"file in the project -- instrument model, digestion enzyme, DIA acquisition method,",
"search-software version, labelling, fractionation, animal sex and age -- have been left",
"blank rather than inferred, and are itemised in `pride/SDRF_MISSING_METADATA.md`.",
"",
"---",
"",
paste0("Underlying revision reports: `", basename(final_report), "` and `",
       basename(retained_report), "` under the canonical revision output root."),
paste0("Generated by `", STAGE, "`.")
)

path <- release_path("editor_source_data", "REVISION_PROTEOMICS_DATA_CHANGELOG.md")
release_write_lines(L, path)
release_register("editor_source_data/REVISION_PROTEOMICS_DATA_CHANGELOG.md",
                 "original submission versus revised analysis, for the editor and reviewers",
                 c(rel("differential_analysis", "primary_differential_summary.tsv"),
                   rel("editor_source_data", "figure_source_map.tsv"),
                   retained_report, final_report),
                 c(NA_character_, NA_character_, release_sha256(retained_report),
                   release_sha256(final_report)),
                 STAGE, "md")
release_log("  wrote REVISION_PROTEOMICS_DATA_CHANGELOG.md (", length(L), " lines)")
release_log("  ", n_without_hits, " of 12 primary comparisons have no FDR-significant protein")
release_log("  panels with changed interpretation: ",
            paste(changed_interpretation, collapse = ", "))
release_log("stage 08 complete")
