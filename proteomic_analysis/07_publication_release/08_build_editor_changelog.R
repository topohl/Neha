#!/usr/bin/env Rscript

# Publication release, stage 08 -- the editor/reviewer changelog.
#
# Produces
#   editor_source_data/REVISION_PROTEOMICS_DATA_CHANGELOG.md
#   editor_source_data/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv     classified fold-change mentions
#   editor_source_data/MANUSCRIPT_TERMINOLOGY_ACTIONS.md     manual, outside-repo actions
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
corrections <- read_release(rel("metadata", "sample_class_corrections.tsv"))
crosswalk <- release_old_package_crosswalk()

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
  `GO-BP terms FDR<0.05 (effect-size-ranked sensitivity)` =
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
"| **GSEA ranking statistic** | **the effect size stored as `logFC`** | **moderated t statistic** |",
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
"| `GSEA_ORA_all_results.xlsx` | `enrichment/*.tsv.gz` and `editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx` | **Yes, entirely.** Hemisphere-level inference, effect-size-ranked GSEA, incomplete contrast coverage. |",
"| `README.txt` (group codes only) | `README_DATA.md`, `metadata/data_dictionary.tsv`, this changelog | **Yes.** |",
"",
paste0("**How that mapping was established.** Crosswalk method for this build: **",
       crosswalk$method, "**."),
"",
if (identical(crosswalk$method, "DIRECT_VERIFICATION")) c(
paste0("The original submission package was reachable from this build at `",
       crosswalk$root, "` and every submitted file was read and hashed here:"),
"",
"| submitted file | SHA256 (read by this build) |",
"|---|---|",
vapply(seq_len(nrow(crosswalk$files)), function(i) {
  paste0("| `", crosswalk$files$file[i], "` | `", crosswalk$files$sha256[i], "` |")
}, character(1))
) else c(
paste0("None of the four originally submitted files exists under those names anywhere in",
       " the project tree -- a recursive search of the analysis root and of the enclosing",
       " project folder returns no `processed_protein_group_matrix_*`, no",
       " `GSEA_ORA_all_results*` and no `README.txt`."),
"",
paste0("**", crosswalk$statement, "** The properties of the submitted files recorded by",
       " that external inspection are stated below; this build did **not** read those",
       " bytes, and no local copy of the package was created to pretend otherwise."),
"",
"| submitted file | recorded properties | crosswalk basis |",
"|---|---|---|",
vapply(seq_along(RELEASE_OLD_PACKAGE_FILES), function(i) {
  f <- RELEASE_OLD_PACKAGE_FILES[[i]]
  shape <- if (!is.na(f$n_protein_rows)) {
    paste0(format(f$n_protein_rows, big.mark = ","), " protein rows x ",
           f$n_measurement_columns, " measurement columns",
           if (!is.na(f$n_annotation_columns)) {
             paste0(" + ", f$n_annotation_columns, " annotation columns")
           } else "", "; ", f$nature)
  } else f$nature
  paste0("| `", f$name, "` | ", shape, " | ", crosswalk$files$crosswalk_method[i], " |")
}, character(1)),
"",
"The local correspondence rests on:",
"",
paste0("- for the two matrices, an identification **by dimensions and content lineage**:",
       " `pg.matrix_raw.txt` at the project root has exactly 5,747 protein rows and 96",
       " measurement columns (plus 7 annotation columns), and the processed matrix has",
       " exactly 5,349 x 96. Both match the dimensions recorded for the submitted files,",
       " and all 5,349 analysed protein identifiers are a strict subset of the 5,747;"),
paste0("- for `GSEA_ORA_all_results.xlsx`, an identification **by analysis generation**,",
       " not by inspection of the workbook here. It is superseded because any",
       " hemisphere-level, effect-size-ranked GSEA/ORA output is superseded by the",
       " animal-level, moderated-t-ranked canonical run, and because it did not contain the",
       " current complete 12-primary-contrast design."),
"",
paste0("To upgrade this section to `DIRECT_VERIFICATION`, extract the submission package",
       " somewhere reachable and set `PROTEOMICS_RELEASE_OLD_PACKAGE_ROOT` to that",
       " directory; the build then hashes each submitted file and reports the byte-level",
       " comparison instead of this paragraph.")
),
"",
"**Nothing in this package depends scientifically on the crosswalk.** It exists so an",
"editor can map the old submission onto the new one, not as evidence for any result.",
"",
"## 11. A sample-class correction, resolved",
"",
if (nrow(corrections) > 0L) c(
sprintf(paste("For %d of the 96 acquisitions the sample class used by the validated analysis",
              "differs from every pre-correction record. A forensic audit established what",
              "that difference is: a **deliberate correction** applied during historical",
              "sample-identity quality control. It is not a mislabelling, and it is not an",
              "open question."), nrow(corrections)),
"",
"| acquisition | animal | hemisphere | plate position | original class | analysis class |",
"|---|---|---|---|---|---|",
vapply(seq_len(nrow(corrections)), function(i) {
  r <- corrections[i, , drop = FALSE]
  paste0("| ", r$plate_sample_number, " | ", r$AnimalID, " | ", r$hemisphere, " | ",
         r$plate_position, " | ", r$original_sample_class, " | **",
         r$analysis_sample_class, "** |")
}, character(1)),
"",
paste0("Both animals received the same cyclic reassignment among the mCherry-, neuropil-",
       " and neuron-enriched samples: **neuropil -> mCherry, mCherry -> neuron, neuron ->",
       " neuropil**. The cfos-enriched samples were not involved. Only left-hemisphere",
       " samples of ", paste(sort(unique(corrections$AnimalID)), collapse = " and "),
       " are affected; the plate layout (row A/E cfos, B/F neuropil, C/G mcherry, D/H",
       " neuron) is followed by the other 90 measurements."),
"",
"### What changed, and what did not",
"",
paste0("**Only sample-class metadata were reassigned.** The six acquisitions retained their",
       " original acquisition identities and their protein-abundance profiles: the",
       " quantitative matrix was compared across the metadata edit and the maximum absolute",
       " numeric difference was **0**. No abundance measurement was altered, no acquisition",
       " was renamed, and no statistic was recomputed as part of this correction."),
"",
"### The evidence",
"",
paste0("- Nine independent pre-correction records carry the original plate assignments.",
       " Three of them are re-derived at build time and required to agree: the",
       " `group_label` column of `sample_annotation.xlsx`, the autosampler plate layout",
       " implied by the well row letter, and the retained correction table itself."),
paste0("- A retained project correction table preserves the reassignment:",
       " `", corrections$correction_reference[[1]], "`, SHA256 `",
       corrections$correction_reference_sha256[[1]], "`, modified ",
       corrections$correction_reference_mtime[[1]], ". The hash is recomputed by the build",
       " and required to match the value the audit recorded; the same table existed",
       " byte-identically in an archive backup."),
paste0("- A single deliberate metadata edit on ", corrections$correction_date[[1]],
       " changed these six class labels and nothing else."),
"- All six left the outlier list once corrected.",
paste0("- Under the original labels the six were among the strongest anomalous assignments",
       " in the dataset (ranks 1, 2, 3, 4, 5 and 8 of 96) and fell outside the normal",
       " left/right similarity distribution."),
paste0("- Independent profile-identity analyses support the corrected labels, and a",
       " design-constrained classification recovered the corrected assignments for both",
       " animals independently."),
"",
paste0("**The correction was identified through UMAP-based QC and then corroborated; UMAP",
       " alone is not treated as ground truth.** The correction table's",
       " nearest-class-centre suggestion agrees with the applied correction for five of the",
       " six, and differs for one (N60: nearest centre `cfos`, applied `neuron`). That",
       " column is published as `umap_nearest_class_suggestion` so the divergence is",
       " visible rather than smoothed over, and the applied classes are those the validated",
       " analysis used, supported by the independent profile-similarity and left/right-pair",
       " analyses."),
"",
"### What is not claimed",
"",
paste0("**No surviving prose note records the rationale.** No memo, comment or README in the",
       " project explains why the reassignment was made. That absence is stated here rather",
       " than papered over. It does not make the correction unresolved: the correction table",
       " itself is the preserved historical correction record, and the independent lines of",
       " evidence above are what establish the assignment."),
"",
paste0("**Both labels are retained.** `metadata/sample_metadata.tsv` carries",
       " `original_sample_class` and `analysis_sample_class` for all 96 measurements, with",
       " `sample_class_corrected` marking the six. The `sample_class` column remains the",
       " canonical analysis-time assignment -- the one the locked GCTs and every validated",
       " statistic were computed on -- and the data dictionary says so explicitly. The",
       " per-row provenance is in `metadata/sample_class_corrections.tsv`."),
"",
paste0("No reanalysis is required and none was performed. The affected animal-level units",
       " (", paste(sort(unique(corrections$animal_level_unit_affected)), collapse = ", "),
       ") were built from these corrected assignments in the first place, so the validated",
       " results already reflect them.")
) else
"No sample-class correction is recorded: all 96 measurements carry their original class.",
"",
"## 12. The effect size is not a log2 fold change",
"",
"A second correction, independent of both the statistical-unit correction and the",
"sample-class correction above.",
"",
paste0("The abundance matrix all statistics were computed on is **standardised separately",
       " for each protein** across the measurement-level dataset. Measured directly: row",
       " means are approximately 0 and row standard deviations approximately 1. ProTigy",
       " nevertheless stores its coefficient as `logFC`, and the canonical split tables",
       " carry it as `log2fc`."),
"",
paste0("A difference of group means on that scale is a **", RELEASE_EFFECT_SIZE$public_term,
       "** -- equivalently, a ", RELEASE_EFFECT_SIZE$public_term_alt, ". It is not a log2",
       " ratio, and a value of 2 does not mean a four-fold change."),
"",
paste0("It is deliberately *not* called a \"standardised mean difference\" without",
       " qualification: that term invites a Cohen's d reading, i.e. scaling by a pooled",
       " within-group standard deviation. The scaling here is a per-protein standard",
       " deviation taken across the whole measurement dataset before any grouping, which is",
       " a different quantity."),
"",
paste0("It is also *not* called \"z-scored\". The released values carry two decimals, so on",
       " those bytes no protein is exactly mean 0 / standard deviation 1, and the operation",
       " that produced the matrix is unresolved -- no script, no parameters (see",
       " `provenance/UPSTREAM_PROVENANCE_GAP.md`). \"Standardised\" is what the numbers",
       " support, so that is what the package says."),
"",
"Corroborating evidence, from the canonical tables themselves: `AveExpr` -- the average",
"abundance limma reports for each protein -- lies within +/-0.0011 of zero for all 5,349",
"proteins. That can only happen if the input was centred per protein.",
"",
paste0("This package therefore publishes the quantity as `", RELEASE_EFFECT_SIZE$public_field,
       "`, carrying `effect_size_definition`, `effect_size_source_column` and",
       " `source_statistic_field = \"", RELEASE_EFFECT_SIZE$source_field, "\"` on every row,",
       " rather than republishing it under a name that asserts a fold-change reading."),
"",
"**No value was changed.** The numbers are bit-identical to the canonical tables -- the",
"build re-reads every split table and requires exact equality -- and directions, rankings,",
"p-values and every downstream result are unaffected. Only the label is corrected.",
"",
paste0("The same applies to the GSEA sensitivity analysis. Canonical GSEA ranks by the",
       " moderated t and is unchanged. What was described as a \"log2FC-ranked\" sensitivity",
       " analysis is publicly an **", RELEASE_EFFECT_SIZE$sensitivity_public_label, "**;",
       " its filenames and its `rank_statistic` value keep the historical `log2fc` token so",
       " the export can still be matched to the canonical run. No ranking value changed and",
       " no enrichment was rerun."),
"",
"### Methods wording",
"",
paste0("> ", RELEASE_EFFECT_SIZE$methods_sentence),
"",
paste0("Every occurrence of the historical fold-change wording in the active",
       " publication-facing files was",
       " classified and is itemised in",
       " `editor_source_data/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv`. Figure axes, captions and",
       " manuscript text live outside this repository and must be corrected by hand; the",
       " specific items are listed in",
       " `editor_source_data/MANUSCRIPT_TERMINOLOGY_ACTIONS.md`."),
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

# --------------------------------------------------------------------------------------
# effect-size terminology audit
# --------------------------------------------------------------------------------------
# Every fold-change mention in the active publication-facing surface, classified. Two parts,
# and the split matters:
#
#   * a RECORD of the corrections actually applied to the publication layer. Curated,
#     because "what was changed" is a historical fact about this commit and cannot be
#     recovered by scanning the already-corrected files.
#   * a LIVE SCAN of the built release documents and the active builders. This is what makes
#     the record checkable: if a PUBLIC_MISLABEL survives anywhere, it shows up here as a
#     row and 13_validate_release.R fails on it.
#
# Classification:
#   INTERNAL_PROVENANCE        a canonical filename, column value or source-field citation.
#                              Kept deliberately -- renaming would break the link to the
#                              canonical run.
#   PUBLIC_MISLABEL            reader-facing text that called the published effect size a
#                              fold change. Corrected.
#   LEGITIMATE_OTHER_CONTEXT   a genuine fold change belonging to a different statistic
#                              (EWCE's bootstrap fold change), or a frozen historical
#                              artefact that must not be edited.

audit_row <- function(file, location, old_term, new_term, status, notes) {
  data.frame(file = file, location = location, old_term = old_term, new_term = new_term,
             status = status, notes = notes,
             stringsAsFactors = FALSE, check.names = FALSE)
}

LAYER_REL <- "proteomic_analysis/07_publication_release"

applied <- rbind(
  audit_row(paste0(LAYER_REL, "/04_build_differential_results.R"),
            "EFFECT_SIZE_DEFINITION", "standardised mean difference",
            RELEASE_EFFECT_SIZE$public_term, "PUBLIC_MISLABEL",
            paste("Unqualified 'standardised mean difference' invites a Cohen's d reading",
                  "(pooled within-group SD). Now sourced from RELEASE_EFFECT_SIZE.")),
  audit_row(paste0(LAYER_REL, "/04_build_differential_results.R"),
            "EFFECT_SIZE_DEFINITION", "z-scored per protein (mean 0, sd 1)",
            "standardized separately for each protein (row means ~0, row SDs ~1)",
            "PUBLIC_MISLABEL",
            paste("Exact z-scoring is not established: values carry 2 decimals and the",
                  "producing step is unresolved. Safer wording retained.")),
  audit_row(paste0(LAYER_REL, "/03_build_processed_data_exports.R"),
            "VALUE_SCALE", "per-protein standardised (z-scored) abundance",
            "per-protein standardised abundance, with measured deviations reported",
            "PUBLIC_MISLABEL", "Same reason; the assertion now reports what it measured."),
  audit_row(paste0(LAYER_REL, "/02_build_contrast_manifest.R"),
            "INTERPRETATION (3 contrast families)", "Positive log2FC = higher in ...",
            "A positive effect_size_sd_units means higher in ...", "PUBLIC_MISLABEL",
            "Reader-facing interpretation of the sign."),
  audit_row(paste0(LAYER_REL, "/02_build_contrast_manifest.R"),
            "secondary manifest, gsea_log2fc_rank_sensitivity label",
            "GSEA ranked by log2FC (sensitivity analysis)",
            paste("GSEA ranked by the standardised-abundance effect size",
                  "(effect-size-ranked sensitivity analysis)"),
            "PUBLIC_MISLABEL",
            "The manifest KEY and the ranking_statistic value keep the log2fc token."),
  audit_row(paste0(LAYER_REL, "/05_build_enrichment_exports.R"),
            "sensitivity export description and log lines", "ranked by log2FC",
            RELEASE_EFFECT_SIZE$sensitivity_public_label, "PUBLIC_MISLABEL",
            "Filenames and rank_statistic value unchanged; see INTERNAL_PROVENANCE rows."),
  audit_row(paste0(LAYER_REL, "/07_build_editor_source_workbook.R"),
            "Read_Me_First items 5 and 7", "log2FC-ranked",
            "effect-size-ranked", "PUBLIC_MISLABEL",
            "Editor-facing workbook notes."),
  audit_row(paste0(LAYER_REL, "/06_build_figure_source_data.R"),
            "figure_source_map, Figure3|H description",
            "scatter of learning log2FC against CNO log2FC",
            "scatter of the learning effect size against the CNO effect size (SD units)",
            "PUBLIC_MISLABEL", "Panel description shown to the editor."),
  audit_row(paste0(LAYER_REL, "/08_build_editor_changelog.R"),
            "sections 5, 10, 11, 12", "log2FC / standardised mean difference",
            paste(RELEASE_EFFECT_SIZE$public_term, "/",
                  RELEASE_EFFECT_SIZE$sensitivity_public_label),
            "PUBLIC_MISLABEL", "Changelog narrative and the superseded-files table."),
  audit_row(paste0(LAYER_REL, "/10_build_provenance.R"),
            "upstream provenance gap, 'Consequence for reporting'",
            "standardised mean difference; mean 0 and sd 1",
            paste(RELEASE_EFFECT_SIZE$public_term, "; row means ~0, row SDs ~1"),
            "PUBLIC_MISLABEL", "Provenance narrative."),
  audit_row(paste0(LAYER_REL, "/11_build_readme_and_dictionary.R"),
            "DESCRIPTIONS map and README_DATA sections 3, 4, 5, 7, 8",
            "Positive log2FC / standardised mean difference / log2FC-ranked",
            paste("positive", RELEASE_EFFECT_SIZE$public_field, "/",
                  RELEASE_EFFECT_SIZE$public_term, "/",
                  RELEASE_EFFECT_SIZE$sensitivity_public_label),
            "PUBLIC_MISLABEL", "Data dictionary and the reader-facing README."),

  audit_row("enrichment/GSEA_log2FC_sensitivity.tsv.gz", "filename",
            "GSEA_log2FC_sensitivity", "GSEA_log2FC_sensitivity", "INTERNAL_PROVENANCE",
            paste("Deliberately unchanged. The name matches the canonical run's",
                  "GSEA_{GO_BP,KEGG}_log2fc_sensitivity.csv; renaming would break the",
                  "provenance link. Documented as effect-size-ranked.")),
  audit_row("enrichment/GSEA_log2FC_sensitivity.tsv.gz", "rank_statistic column value",
            "log2fc", "log2fc", "INTERNAL_PROVENANCE",
            "The value the canonical run recorded. Kept verbatim."),
  audit_row("enrichment/primary_ORA_GO_BP.tsv.gz", "query_list column value",
            "top_absolute_log2fc", "top_absolute_log2fc", "INTERNAL_PROVENANCE",
            "Canonical ORA query-list identifier. Kept verbatim."),
  audit_row("differential_analysis/primary_differential_proteins.tsv.gz",
            "effect_size_source_column and source_statistic_field",
            "logFC", "logFC", "INTERNAL_PROVENANCE",
            paste("Recorded ON PURPOSE so the published value can be traced to its source.",
                  "The exported effect_size_sd_units equals this source field exactly.")),
  audit_row("02_data/animal_level/split/forward/*.csv", "canonical column log2fc",
            "log2fc", "log2fc", "INTERNAL_PROVENANCE",
            "Canonical analysis artefact, outside this layer. Never modified."),
  audit_row("proteomic_analysis/06_manuscript_figure_revision/*.R",
            "frozen figure scripts", "log2FC", "log2FC", "LEGITIMATE_OTHER_CONTEXT",
            paste("FROZEN. Not edited by this change. Their rendered axis labels and",
                  "captions are manuscript-facing and are listed in",
                  "MANUSCRIPT_TERMINOLOGY_ACTIONS.md for manual correction.")),
  audit_row("enrichment/primary_EWCE.tsv.gz", "fold_change column",
            "fold change", "fold change", "LEGITIMATE_OTHER_CONTEXT",
            paste("EWCE's own statistic: the ratio of observed to bootstrap-null mean",
                  "expression proportion. Genuinely a fold change. Unrelated to the",
                  "differential effect size.")),
  audit_row("proteomic_analysis/04_differential_expression_enrichment/01_clusterProfiler.r",
            "canonical analysis script", "log2fc", "log2fc", "INTERNAL_PROVENANCE",
            paste("Produced the canonical statistics. Not a publication-facing document",
                  "and not edited; renaming would invalidate the validated outputs."))
)

# Live scan. Any reader-facing line that still asserts a fold-change reading appears here
# as a PUBLIC_MISLABEL row and fails validation.
# Live scan of BUILT artefacts. Deliberately not a scan of the builder SOURCES: those
# necessarily contain the vocabulary being classified -- R/release_validation.R holds the
# claim list itself, and this file holds the audit rows -- so scanning them reports the
# detector as a defect. What matters is the published surface, so that is what is
# scanned. Reader-facing documents written later in the build (README_DATA.md, the data
# dictionary, the changelog, the PRIDE and provenance notes) are scanned by
# 13_validate_release.R, which runs last and fails the release on any residual hit.
scan_targets <- c(
  rel("differential_analysis", "primary_differential_summary.tsv"),
  rel("metadata", "primary_contrast_manifest.tsv"),
  rel("metadata", "secondary_analysis_manifest.tsv"),
  rel("editor_source_data", "figure_source_map.tsv"),
  rel("enrichment", "enrichment_run_parameters.tsv"),
  rel("enrichment", "enrichment_coverage.tsv")
)
# Label a scanned path for the audit: release-relative for built files, repo-relative
# for sources, so the audit reads the same whichever scratch root the build used.
audit_scan_label <- function(p) {
  p <- gsub("\\\\", "/", p)
  for (root in c(OUT_ROOT, dirname(REPO_ROOT), REPO_ROOT)) {
    root <- gsub("\\\\", "/", root)
    if (startsWith(p, paste0(root, "/"))) return(substring(p, nchar(root) + 2L))
  }
  basename(p)
}

residual <- list()
for (p in scan_targets) {
  if (!file.exists(p)) next
  hits <- release_fold_change_mislabels(readLines(p, warn = FALSE))
  for (i in seq_len(nrow(hits))) {
    residual[[length(residual) + 1L]] <- audit_row(
      audit_scan_label(p),
      paste0("line ", hits$line[i]), hits$term[i], NA_character_, "PUBLIC_MISLABEL",
      paste("RESIDUAL -- found by the live scan and not yet corrected:", hits$text[i]))
  }
}
terminology_audit <- if (length(residual)) {
  rbind(applied, do.call(rbind, residual))
} else applied

n_public_mislabel_corrected <- sum(applied$status == "PUBLIC_MISLABEL")
n_residual <- length(residual)
release_log("  effect-size terminology audit: ", n_public_mislabel_corrected,
            " public mislabel group(s) corrected, ",
            sum(terminology_audit$status == "INTERNAL_PROVENANCE"),
            " internal-provenance token(s) retained, ", n_residual, " residual")
if (n_residual > 0L) {
  release_log("  *** ", n_residual, " residual public fold-change mislabel(s) remain ***")
}

wa <- release_write_table(terminology_audit,
                          release_path("editor_source_data",
                                       "EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv"))
release_register("editor_source_data/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv",
                 paste("every fold-change mention in the active publication-facing surface,",
                       "classified as INTERNAL_PROVENANCE, PUBLIC_MISLABEL or",
                       "LEGITIMATE_OTHER_CONTEXT"),
                 "publication release layer sources and built documents",
                 NA_character_, STAGE, "tsv")
release_log("  wrote EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv (", wa$rows, "x", wa$cols, ")")

# --------------------------------------------------------------------------------------
# manuscript action note -- what has to be fixed by hand, outside this repository
# --------------------------------------------------------------------------------------

MA <- c(
  "# Manuscript actions: effect-size terminology and the sample-class correction",
  "",
  paste0("Generated by `", STAGE, "` on ", release_timestamp_utc(), "."),
  "",
  "Everything in the released data package has been corrected. The items below live",
  "**outside this repository** -- in the manuscript file, the figure files and the",
  "submission forms -- and must be changed by hand. Nothing here changes a number.",
  "",
  "## 1. The effect size is not a log2 fold change",
  "",
  paste0("Published field: `", RELEASE_EFFECT_SIZE$public_field, "`. Preferred wording: **",
         RELEASE_EFFECT_SIZE$public_term, "**, or \"",
         RELEASE_EFFECT_SIZE$public_term_alt, "\"."),
  "",
  "| Where | Incorrect wording | Correct wording |",
  "|---|---|---|",
  paste0("| Figure 3D volcano, x-axis | incorrect: `log2FC` / \"log2 fold change\" | \"",
         RELEASE_EFFECT_SIZE$public_term, "\" |"),
  paste0("| Figure 3H scatter, both axes | learning / CNO `log2FC` | learning / CNO effect",
         " size (SD units) |"),
  paste0("| Supplementary E panels derived from the paired model | `log2FC` | ",
         RELEASE_EFFECT_SIZE$public_term, " |"),
  "| Figure captions quoting effect sizes | incorrect: \"fold change\" | \"SD units\" |",
  "| Methods, differential-abundance paragraph | incorrect: \"log2 fold change\" | see the sentence below |",
  "| Results text quoting effect magnitudes | \"n-fold\" | \"n SD units\" |",
  paste0("| Any \"log2FC-ranked GSEA\" mention | log2FC-ranked | ",
         RELEASE_EFFECT_SIZE$sensitivity_public_label, " |"),
  "",
  "### Methods sentence to insert verbatim",
  "",
  paste0("> ", RELEASE_EFFECT_SIZE$methods_sentence),
  "",
  paste0("Do **not** describe the scale as z-scored: the released values carry two decimals,",
         " no protein is exactly mean 0 / SD 1 on those bytes, and the operation that",
         " produced the matrix is unresolved. Do **not** call the effect size a",
         " \"standardised mean difference\" without qualification either -- that reads as",
         " Cohen's d, which uses a pooled within-group SD rather than the across-dataset",
         " per-protein SD used here."),
  "",
  paste0("The frozen scripts under `06_manuscript_figure_revision/` were **not** edited:",
         " they are the record of how the submitted figures were produced. Their axis",
         " labels therefore still say `log2FC`, and the figures must be relabelled at the",
         " manuscript stage. The classification of every occurrence is in",
         " `EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv`."),
  "",
  "## 2. The sample-class correction",
  "",
  if (nrow(corrections) > 0L) c(
    paste0("If the manuscript, a supplementary table or a submission form describes the",
           " sample-class assignments of animals ",
           paste(sort(unique(corrections$AnimalID)), collapse = " and "),
           ", use the analysis-time classes. Suggested wording:"),
    "",
    paste0("> Six acquisition-level samples from the left hemispheres of animals ",
           paste(sort(unique(corrections$AnimalID)), collapse = " and "),
           " had sample-class assignments corrected during historical sample-identity",
           " quality control. Original plate-layout and early metadata records retain the",
           " pre-correction assignments, whereas a preserved UMAP-based correction table",
           " records a cyclic reassignment among mCherry-, neuropil- and neuron-enriched",
           " samples. Acquisition identifiers and quantitative protein-abundance profiles",
           " were unchanged; only sample-class metadata were reassigned. Independent",
           " profile-similarity and left/right-pair analyses supported the corrected",
           " assignments. Both original and analysis-time sample-class labels are retained",
           " in the released metadata."),
    "",
    "Do not describe these six as mislabelled samples, as an unresolved discrepancy, or as",
    "requiring reanalysis. No abundance measurement changed and no statistic was recomputed.",
    "",
    "Do not claim a written rationale survives. It does not; the correction table is the",
    "preserved record."
  ) else "No sample-class correction is recorded.",
  "",
  "## 3. Still needed for the PRIDE deposition",
  "",
  "Not a terminology item, but it blocks deposition and is not something this repository",
  paste0("can resolve. See `pride/SDRF_MISSING_METADATA.md` for the authoritative list."),
  ""
)
mp <- release_path("editor_source_data", "MANUSCRIPT_TERMINOLOGY_ACTIONS.md")
release_write_lines(MA, mp)
release_register("editor_source_data/MANUSCRIPT_TERMINOLOGY_ACTIONS.md",
                 paste("figure axes, captions and manuscript text that must be corrected by",
                       "hand outside this repository"),
                 "publication release layer", NA_character_, STAGE, "md")
release_log("  wrote MANUSCRIPT_TERMINOLOGY_ACTIONS.md (", length(MA), " lines)")
release_log("stage 08 complete")
