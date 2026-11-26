#!/usr/bin/env Rscript

# Publication release, stage 11 -- README_DATA.md and the data dictionary.
#
# Produces
#   README_DATA.md                 the human entry point to the package
#   metadata/data_dictionary.tsv   every column of every release table, described
#
# The dictionary walks the tables in TABLES and looks up each column in
# TABLE_COLUMN_OVERRIDES, DESCRIPTIONS (RELEASE_DESCRIPTIONS in R/release_dictionary.R, shared
# with the workbook), WIDE_NOTE (the two matrices) or, for the SDRF,
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
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_dictionary.R"))
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
figure_map_full <- read_full(rel("provenance", "figure_source_map.tsv"))
sf_row <- figure_map_full[figure_map_full$figure == "Supplementary_proteomics" &
                            figure_map_full$panel == "F", , drop = FALSE]
sf_classes <- unique(as.character(release_read_csv(
  strsplit(sf_row$panel_source_data_canonical_path, ";", fixed = TRUE)[[1]][[1]])$sample_class))
ora_cc_all <- read_full(rel("enrichment", "primary_ORA_GO_CC.tsv.gz"))
ewce_model <- release_ewce_differential_model(DATA_ROOT)

class_metadata <- experimenter_metadata[
  experimenter_metadata$category == "sample_class_definition", , drop = FALSE]
class_metadata <- class_metadata[order(as.integer(class_metadata$sort_order)), , drop = FALSE]
# README table of compartments: data value (`mcherry`), display name ("mCherry+") and
# definition. The definition is display text, so the marker is written c-Fos as in the paper
# (the experimenter records write cFos).
class_markdown <- vapply(seq_len(nrow(class_metadata)), function(i) {
  paste0("| `", class_metadata$sample_class[i], "` | ",
         release_sample_class_display(class_metadata$sample_class[i]), " | ",
         gsub("\\bcFos\\b", "c-Fos", class_metadata$value[i], perl = TRUE), " |")
}, character(1))

# --------------------------------------------------------------------------------------
# column descriptions
# --------------------------------------------------------------------------------------

# One description per machine column, shared with stage 07 (R/release_dictionary.R).
DESCRIPTIONS <- RELEASE_DESCRIPTIONS

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
  "differential_analysis/significant_proteins.tsv",
  "enrichment/significant_gene_sets.tsv.gz",
  "enrichment/primary_GSEA_GO_BP.tsv.gz", "enrichment/primary_GSEA_GO_CC.tsv.gz",
  "enrichment/primary_GSEA_KEGG.tsv.gz", "enrichment/primary_ORA_GO_BP.tsv.gz",
  "enrichment/primary_ORA_GO_CC.tsv.gz", "enrichment/GSEA_log2FC_sensitivity.tsv.gz",
  "enrichment/primary_EWCE.tsv.gz", "enrichment/enrichment_coverage.tsv",
  "enrichment/gsea_untested_gene_sets.tsv", "enrichment/enrichment_run_parameters.tsv",
  "provenance/figure_source_map.tsv",
  "editor_source_data/workbook_sheet_inventory.tsv",
  "provenance/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv",
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
    "One column per acquisition run, named by the acquisition run name (see metadata/sample_metadata.tsv, column sample_id). Values are processed protein abundance, filtered, imputed and standardized per protein: each protein has a mean of about 0 and an SD of about 1 over the 96 acquisitions. Not log2 abundance; the producing step is unresolved (see provenance/UPSTREAM_PREPROCESSING_GAP.md). NOT raw MS data.",
  "processed_data/protein_abundance_animal_level.tsv.gz" =
    "One column per animal-level unit, named `<AnimalID>_<sample_class>` (see metadata/animal_level_sample_metadata.tsv, column animal_level_column_name). Values are the equal-weight mean of the Left and Right measurements on the same per-protein standardized scale: each protein keeps a mean of about 0, and its SD over the 48 units is below 1 because two hemispheres are averaged. Not log2 abundance.")

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
  "differential_analysis/significant_proteins.tsv" = c(
    direction = paste("higher_in_<numerator_condition> when effect_size_sd_units is positive,",
                      "otherwise higher_in_<denominator_condition>.")
  ),
  "enrichment/significant_gene_sets.tsv.gz" = c(
    direction = paste("enriched_toward_<numerator_condition> when NES is positive, otherwise",
                      "enriched_toward_<denominator_condition>; copied from the GSEA table.")
  ),
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
  "provenance/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv" = c(
    file = "File or released artifact the occurrence is in.",
    location = "Where in that file: a named constant, section, column or line number.",
    old_term = "The historical fold-change wording that was found at this location.",
    new_term = "What it now reads, or the same token where it was kept for provenance.",
    status = paste("Classification of the occurrence. INTERNAL_PROVENANCE: a canonical",
                   "filename, column value or source-field citation, kept",
                   "because renaming it would break the link to the canonical run.",
                   "PUBLIC_MISLABEL: reader-facing text that described the published effect",
                   "size as a fold change; corrected. LEGITIMATE_OTHER_CONTEXT: a genuine",
                   "fold change belonging to a different statistic (EWCE's bootstrap fold",
                   "change), or a frozen historical artifact that must not be edited."),
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

# Collection plate: which animals are on which plate, checked against the sample metadata.
plate_layout_text <- release_plate_layout_sentence(read_full(rel("metadata", "sample_metadata.tsv")))

R <- c(
"# Associative memory proteomics: data package",
"",
paste0("Animal-level proteomics of four marker-defined compartments collected by laser-capture",
       " microdissection from ", mv("organism_part"), " in a chemogenetic associative-memory"),
"experiment. This package contains the processed data, the differential and enrichment",
"results, the figure source data and the metadata needed to interpret and redeposit them.",
"Supplementary Table 2 (`editor_source_data/Supplementary_Table_2_proteomics.xlsx`) shows the",
"main tables with display labels; this file describes the complete tables.",
"",
"**The unit of analysis is the animal**, not the hemisphere: the left and right runs of each",
"animal and compartment are averaged before any test (section 2).",
"",
paste0("**Statistics.** Differential abundance is tested with limma, one model per",
       " compartment (residual df ", stat_run$df_residual, ")."),
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
"From each animal, four marker-defined compartments (`sample_class`) were collected from each",
"hemisphere:",
"",
"| Compartment (data value) | Display name | What it is |",
"|---|---|---|",
class_markdown,
"",
"These are sampling compartments, not four Cell Ontology cell types. The SDRF therefore uses",
"`characteristics[cell type] = not applicable`, while `characteristics[sampling site]` and",
"`factor value[sample class]` keep the marker-defined compartment.",
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
"12 animals x 4 compartments x 2 hemispheres = **96 acquisition measurements**.",
"3 animals per condition, so **n = 3 animals** in every condition x compartment stratum,",
"with complete left/right pairing in all 16 strata.",
"",
"## 2. Averaging the hemispheres",
"",
"Left and right measurements were averaged with equal weight within each",
"`AnimalID x sample_class` unit, on the processed values (standardized per protein,",
"see section 3), **before any statistical test**. 96 measurements become",
"**48 animal-level units**, and every test has n = 3 animals per group.",
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
"**These matrices are standardized per protein, not log2 abundance.** Measured directly",
"on the measurement-level matrix: each protein is standardized separately across the 96",
"acquisitions, giving row means of approximately 0 and row standard deviations of",
"approximately 1. The build asserts this and refuses to run if it ever stops holding.",
paste0("The animal-level matrix is the mean of the two hemispheres on that scale: each protein",
       " keeps a mean of approximately 0, and its standard deviation over the 48 units is below",
       " 1 (", sprintf("%.2f", min(animal_row_sd)), " to ", sprintf("%.2f", max(animal_row_sd)),
       ", median ", sprintf("%.2f", stats::median(animal_row_sd)), ") because averaging two",
       " hemispheres reduces the spread."),
"",
"This matters for reading the results. limma stores its coefficient as `logFC`, but a",
"difference of group means",
paste0("on a protein-standardized scale is a **", RELEASE_EFFECT_SIZE$public_term, "** --"),
paste0("equivalently, a ", RELEASE_EFFECT_SIZE$public_term_alt, ". It is not a log2 ratio."),
"An effect size of 2 means \"two standard deviations on the standardized scale\", not",
paste0("\"four-fold\". This package therefore publishes the column as `",
       RELEASE_EFFECT_SIZE$public_field, "`, with `effect_size_definition`,"),
paste0("`effect_size_source_column` and `source_statistic_field = \"",
       RELEASE_EFFECT_SIZE$source_field, "\"` on every row. The values are limma's"),
"coefficients unchanged (the build requires exact equality); only the name differs.",
"",
"Two phrasings are avoided:",
"",
"- **not** \"standardized mean difference\" on its own, because that invites a Cohen's d",
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
"  step that produced the matrix is also not recorded (section 15). \"Standardized\" is what",
"  the numbers support.",
"",
paste0("Corroboration: `average_standardized_abundance` (limma `AveExpr`) is the mean over the ",
       stat_run$n_columns_in_model, " animal-level units of one compartment, so within a",
       " compartment it is not centered (it ranges from ",
       sprintf("%.2f", ave$range_within_classes[1]), " to ",
       sprintf("%.2f", ave$range_within_classes[2]), "). Averaged over the four",
       " compartments it stays close to zero (largest absolute value over all ",
       format(ave$n_proteins, big.mark = ","), " proteins ",
       format(signif(ave$max_abs_mean_over_classes, 3), scientific = FALSE),
       "), which is only possible on a per-protein centered input."),
"",
"For methods sections:",
"",
paste0("> ", RELEASE_EFFECT_SIZE$methods_sentence),
"",
"Do not run inference on the measurement-level matrix. It is provided for the supplementary",
"QC panels that need the acquisition as their unit: hemisphere is averaged out of the",
"animal-level units, and the collection-plate check is made per run.",
"",
"Neither file is raw mass-spectrometry data. See section 20.",
"",
"## 4. The 12 primary comparisons",
"",
"4 compartments x 3 comparisons within each compartment:",
"",
"1. Pairing: `paired_veh` vs `unpaired_veh`",
"2. CNO in paired animals: `paired_cno` vs `paired_veh`",
"3. CNO in unpaired animals: `unpaired_cno` vs `unpaired_veh`",
"",
"All 12 are in `metadata/primary_contrast_manifest.tsv`. A positive",
"`effect_size_sd_units` always means higher in the numerator condition (group 1).",
"",
"## 5. Secondary analyses",
"",
"Kept in a separate file, `metadata/secondary_analysis_manifest.tsv`, and never mixed into",
"the primary table: the Supp. Fig. 4 E cross-compartment within-animal paired model, the",
"effect-size-ranked GSEA sensitivity analysis, the post-hoc Pairing x CNO interaction, the",
"leave-one-animal-out check of the mCherry+ Pairing comparison, the EWCE baseline analysis,",
"a pooled model fitted over all four compartments together (for comparison only) and two",
"candidate contrasts that were not analyzed.",
"",
"## 6. Differential statistics",
"",
"`differential_analysis/primary_differential_proteins.tsv.gz` -- one row per protein group",
paste0("per comparison, 12 x 5,349 = 64,188 rows. limma ", stat_run$limma_version, " (R ",
       stat_run$r_version, "), one linear model per compartment (~0 + condition over that",
       " compartment's ", stat_run$n_columns_in_model, " animal-level units, n = ",
       RELEASE_DESIGN_INVARIANTS$n_animals_per_stratum, " animals per condition, residual df ",
       stat_run$df_residual, "), empirical Bayes moderation (robust, no intensity trend),",
       " Benjamini-Hochberg FDR within each comparison. Each comparison is 3 vs 3 animals.",
       " Model parameters, including the prior variance of each compartment, are in",
       " `provenance/analysis_parameters.tsv`."),
"",
"Two protein counts appear throughout and mean different things:",
"",
"- **5,349** protein groups were tested.",
"- **5,327** carry a mouse identifier and could be annotated and used for enrichment. The",
"  22-row difference is non-mouse entries: human keratins and other skin proteins,",
"  digestion enzymes (porcine trypsin, Lys-C and a bacterial protease) and bovine serum",
"  proteins. The manuscript figures were drawn from the 5,327 set, so both counts are",
"  reported side by side in `primary_differential_summary.tsv`.",
"",
sprintf("Of the 12 primary comparisons, %d %s at least one FDR-significant protein and %d %s none.",
        n_with_hits, if (n_with_hits == 1L) "contains" else "contain",
        12L - n_with_hits, if (12L - n_with_hits == 1L) "contains" else "contain"),
"",
paste0("`differential_analysis/significant_proteins.tsv` lists the ",
       format(sum(as.integer(diff_summary$n_significant_fdr_0_05_mapped_only)), big.mark = ","),
       " annotated proteins at FDR < 0.05 (the rows of `primary_differential_proteins.tsv.gz`",
       " with `significant_fdr_0_05` TRUE and `id_mapping_status` mapped; Supplementary Table",
       " 2d). The ", sum(as.integer(diff_summary$n_significant_fdr_0_05)) -
         sum(as.integer(diff_summary$n_significant_fdr_0_05_mapped_only)),
       " non-mouse groups that also pass FDR < 0.05 are in `primary_differential_proteins.tsv.gz`",
       " only."),
"",
"## 7. GSEA",
"",
"**The primary GSEA ranks by the moderated t statistic.** The effect-size-ranked GSEA is a",
"sensitivity analysis only and is in `enrichment/GSEA_log2FC_sensitivity.tsv.gz` (the file",
"name keeps the analysis run's label), with `analysis_role = sensitivity` on every row. Never",
"pool the two.",
"",
"Settings: GO-BP, GO-CC and KEGG, minGSSize 10, maxGSSize 800, BH adjustment,",
"`pvalue_cutoff = 1` and `qvalue_cutoff = 1` -- so the exported tables hold **every tested",
"term** that fgsea could estimate, not a significant-only subset. Filter on",
"`adjusted_p_value` yourself. GSEA P values come from fgsea permutations with one recorded",
"seed per comparison (`gsea_seed` in `enrichment/enrichment_coverage.tsv`); term counts near",
"FDR 0.05 can move by a few terms with the seed.",
paste0("KEGG gene sets were downloaded at run time (", RELEASE_STATISTICS_RUN$kegg_enrichment_date,
       "); the KEGG release was not recorded."),
"",
paste0("**Permutations and untested gene sets.** ", permutations$reader_sentence, " Gene sets",
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
      "`enrichment/primary_ORA_GO_CC.tsv.gz`."),
"",
paste0("`enrichment/significant_gene_sets.tsv.gz` lists the gene sets at adjusted P < 0.05 in",
       " the three primary GSEA tables (GO-BP, GO-CC, KEGG; ",
       format(nrow(read_full(rel("enrichment", "significant_gene_sets.tsv.gz"))),
              big.mark = ","),
       " rows; Supplementary Table 2f), values unchanged, with the leading-edge members as",
       " gene symbols in `core_enrichment_gene_symbols`. GO members map from UniProt accession",
       " and KEGG members from NCBI Gene ID, through the KEGG run's own Entrez-to-UniProt",
       " record, to the gene symbols of `processed_data/protein_feature_annotation.tsv.gz`."),
"",
"## 8. ORA universe",
"",
"All unique successfully mapped measured UniProt accessions, computed per comparison. Four",
"query-list variants are exported and distinguished by the `query_list` column: all",
"FDR-significant proteins, FDR-significant higher in the numerator, FDR-significant higher",
"in the denominator, and an alternative list defined by the absolute effect size",
"(`top_absolute_log2fc`, the label the analysis run gave it; no FDR filter).",
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
"Protein identifiers in this package:",
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
"`analysis_role = sensitivity`. The cell-type reference is the ewceData mouse cortex and",
"hippocampus single-cell reference (`ewceData::ctd()`; Zeisel et al. 2015).",
"",
paste(ewce_model$sentence, "Joining EWCE to the differential tables on",
      "`canonical_comparison` matches the comparison, not the protein ranking."),
"",
"## 11. PCA",
"",
"Animal-level PCA over all 5,349 proteins and 48 units, 0 zero-variance rows removed.",
"PC1 29.1%, PC2 20.5%, PC3 12.1%. Acquisition-level PCA is used only for the two",
"technical QC panels.",
"",
"## 12. Figure source data",
"",
paste0("`editor_source_data/figure_source_data/` holds the source data of each proteomics",
       " panel of Figure 2 and Supplementary Figure 4 as CSV files named by panel (e.g.",
       " `Fig2D_source_data.csv`, `SuppFig4F_source_data.csv`). Where a column holds the same",
       " quantity as a column of Supplementary Table 2, its header is the display header used",
       " there, and compartments, conditions and comparisons use the same display labels.",
       " Table 2p of Supplementary Table 2 gives, for",
       " each panel, the table it is drawn from, the rows used and the key numbers.",
       " `provenance/figure_source_map.tsv` names the analysis file each panel was built from,",
       " with its SHA256. Supp. Fig. 4 A gives the protein and peptide counts of each",
       " acquisition."),
"",
"## 13. Known design limitations",
"",
if (nrow(corrections) > 0L) c(
sprintf(paste("- **Sample-class correction (%d of 96 measurements).** Six",
              "acquisition-level samples from the left hemispheres of animals %s had",
              "sample-class assignments corrected during historical sample-identity quality",
              "control. The plate layout and the early metadata records keep the",
              "assignments before the correction, and a preserved UMAP-based correction table",
              "records a cyclic reassignment among mCherry-, neuropil- and neuron-enriched",
              "samples. Acquisition identifiers and quantitative protein-abundance profiles",
              "were unchanged; only sample-class metadata were reassigned. Independent",
              "profile-similarity and left/right-pair analyses supported the corrected",
              "assignments. Both the original and the analyzed sample-class labels are",
              "kept in the released metadata."),
        nrow(corrections),
        paste(sort(unique(corrections$AnimalID)), collapse = " and ")),
"",
paste0("  The reassignment is the same 3-cycle for both animals: **neuropil -> mCherry,",
       " mCherry -> neuron, neuron -> neuropil**. The c-Fos+ samples were not",
       " involved. Per-row details, including the correction record and its SHA256, are in",
       " `metadata/sample_class_corrections.tsv`; `original_sample_class`,",
       " `analysis_sample_class` and `sample_class_corrected` are on every row of",
       " `metadata/sample_metadata.tsv`."),
"",
paste0("  The preserved correction table is the record of the correction; no written",
       " rationale survives. Nine independent records made before the correction carry the",
       " original assignments; the build re-derives three of them and requires them to agree,",
       " and checks the correction table against its recorded SHA256."),
"",
paste0("  **No reanalysis is needed.** The affected animal-level units (",
       paste(sort(unique(corrections$animal_level_unit_affected)), collapse = ", "),
       ") were built from the corrected assignments, so the results already reflect them.",
       " No abundance value changed: the quantitative matrix was compared across the",
       " metadata edit and the maximum absolute numeric difference was 0.")
) else
"- No sample-class correction is recorded: all 96 measurements carry their original class.",
paste0("- **Low power:** each comparison is 3 vs 3 animals (residual df ",
       stat_run$df_residual, ")."),
paste0("- **One model per compartment.** Each compartment has its own variance estimate.",
       " Results for different compartments come from the same 12 animals and are not",
       " independent of one another, so differences between compartments in what reaches",
       " significance are not themselves tested."),
"- **The cause of the mCherry+ Pairing result is not established.** ",
local({
  mc <- "mcherry_paired_veh_over_mcherry_unpaired_veh"
  up <- function(tbl) {
    tbl[tbl$canonical_comparison == mc & tbl$query_list == "fdr_significant_higher_in_numerator",
        , drop = FALSE]
  }
  ou <- up(ora_go_all); oc <- up(ora_cc_all)
  padj <- function(d) as.numeric(d$adjusted_p_value)
  paste0(
    sprintf(paste("  It has %s FDR-significant proteins of %s annotated (%s higher in",
                  "%s, %s higher in %s). Moderated-t-ranked GSEA finds %d",
                  "GO-BP and %d GO-CC terms at FDR < 0.05. Over-representation among the",
                  "proteins higher in %s finds %d of %s GO-BP terms (minimum adjusted P",
                  "%.2f) and %d of %s GO-CC terms (minimum adjusted P %s) at FDR < 0.05."),
            release_format_count(mch$n_significant_fdr_0_05_mapped_only),
            release_format_count(mch$n_proteins_tested_mapped_only),
            release_format_count(mch$n_significant_higher_in_numerator_mapped_only),
            release_condition_display("paired_veh"),
            release_format_count(mch$n_significant_higher_in_denominator_mapped_only),
            release_condition_display("unpaired_veh"),
            sum(gsea_go_all$canonical_comparison == mc &
                  as.numeric(gsea_go_all$adjusted_p_value) < 0.05, na.rm = TRUE),
            sum(gsea_cc_all$canonical_comparison == mc &
                  as.numeric(gsea_cc_all$adjusted_p_value) < 0.05, na.rm = TRUE),
            release_condition_display("paired_veh"),
            sum(padj(ou) < 0.05, na.rm = TRUE), format(nrow(ou), big.mark = ","),
            min(padj(ou), na.rm = TRUE),
            sum(padj(oc) < 0.05, na.rm = TRUE), format(nrow(oc), big.mark = ","),
            release_format_p(min(padj(oc), na.rm = TRUE))),
    " Pairing is completely associated with collection plate in this comparison and no",
    " technical batch metadata exist, so neither a biological nor a technical explanation",
    " is established.")
}),
local({
  r <- function(col) range(as.numeric(loo_drop[[col]]))
  fdr <- r("n_fdr05"); bp <- r("go_bp_fdr05"); bpk <- r("go_bp_retained_of_full")
  cc <- r("go_cc_fdr05"); cck <- r("go_cc_retained_of_full")
  mc <- "mcherry_paired_veh_over_mcherry_unpaired_veh"
  paste0("- **The mCherry+ Pairing counts depend on single animals.** Refitting with one",
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
         sprintf("%.2f", max(as.numeric(loo_drop$r_effect))), "). These GSEA counts come",
         " from runs at fgsea's default permutations; the primary GSEA (nPermSimple = ",
         format(permutations$go_n_perm_simple, big.mark = ","), ") gives ",
         sum(gsea_go_all$canonical_comparison == mc &
               as.numeric(gsea_go_all$adjusted_p_value) < 0.05, na.rm = TRUE), " GO-BP and ",
         sum(gsea_cc_all$canonical_comparison == mc &
               as.numeric(gsea_cc_all$adjusted_p_value) < 0.05, na.rm = TRUE),
         " GO-CC terms for this comparison.")
}),
local({
  can <- fgsea_na
  plural_sets <- function(n) if (n == 1L) " set" else " sets"
  if (!nrow(can)) {
    paste("- **GSEA permutations.**", permutations$reader_sentence,
          "fgsea estimated a P value for every gene set in every block.")
  } else paste0(
         "- **Some gene sets could not be tested by GSEA.**",
         " When a comparison's moderated t values are",
         " strongly one-sided, fgsea cannot estimate a P value for some large gene sets",
         " (\"unbalanced\" statistics), and clusterProfiler drops them, so they are absent",
         " from the GSEA tables rather than non-significant. This affected ",
         nrow(can), " of the ", 3L * nrow(diff_summary), " comparison x collection blocks",
         " (GO-BP, GO-CC, KEGG; ", min(can$n_sets), " to ",
         max(can$n_sets), " sets each); for the mCherry+ Pairing",
         " comparison, ", if (nrow(fgsea_na_mch)) paste0(sum(fgsea_na_mch$n_sets), " ",
         paste(gsub("_", "-", sub("_GSEA$", "", fgsea_na_mch$analysis)), collapse = "/"),
         plural_sets(sum(fgsea_na_mch$n_sets))) else
           "none", ". The GO sets are listed in `enrichment/gsea_untested_gene_sets.tsv`;",
         " counts per block are in `enrichment/enrichment_coverage.tsv`. ",
         permutations$reader_sentence,
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
paste0("- **The KEGG release was not recorded.** KEGG gene sets were downloaded at run time (",
       RELEASE_STATISTICS_RUN$kegg_enrichment_date, "), so a rerun can change KEGG results for",
       " reasons unrelated to the data."),
"- **Collection plate.** In the Pairing comparisons pairing is completely associated with",
"  collection plate, so the two cannot be separated (section 14).",
"- **Fig. 2 H and Fig. 2 I set two comparisons that share a group against each other.**",
"  Both contain `paired_veh`, on opposite sides, so the design biases them toward",
"  opposition. Descriptive only.",
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
paste0("`collection_plate` (`Plate1`/`Plate2`) records the plate each sample was collected",
       " onto. ", release_plate_injection_sentence(), " ", plate_layout_text,
       " In the Pairing comparisons the association is complete in every compartment. At",
       " acquisition level, collection plate explains 0.01% of PC1 and 1.5% of PC2, neither",
       " significant."),
"",
"## 15. Upstream preprocessing limitation",
"",
"The reproducible pipeline **starts** at the 5,349 x 96 processed matrix. The steps that",
"produced it from the 5,747-row search output -- filtering, log transformation,",
"normalization, imputation and PCA adjustment -- are **not reproducible from this",
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
"| `provenance/data_lineage.tsv` | every artifact from acquisition to release, with unproven edges marked `UNRESOLVED` |",
"| `provenance/software_versions.tsv` | exact versions as recorded by the runs themselves, or `UNKNOWN` with the source that would supply them |",
"| `provenance/analysis_parameters.tsv` | every recorded analysis parameter |",
"| `provenance/figure_source_map.tsv` | for each figure panel, the analysis file it was built from, with its SHA256 |",
"| `provenance/REVISION_PROTEOMICS_DATA_CHANGELOG.md`, `provenance/MANUSCRIPT_TERMINOLOGY_ACTIONS.md`, `provenance/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv` | records of how this package was revised |",
"| `provenance/release_manifest.tsv` | per-file role, shape, SHA256, source artifact and generating script |",
"| `provenance/SHA256SUMS.txt` | checksums for the whole package |",
"| `provenance/sessionInfo_release.txt` | the environment that built the package |",
"",
paste0(length(RELEASE_LOCKED_ARTEFACTS), " upstream artifacts are locked by hash and",
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
"(`comment[data file]`), but **the files themselves are not in this package**.",
"",
paste0("`pride/sdrf.tsv` follows SDRF-Proteomics conventions as far as the verified metadata",
       " allow. ", n_required_missing, " SDRF-required fields could not be populated from",
       " the project or experimenter records and carry `not available`; each is itemized in",
       " `pride/SDRF_MISSING_METADATA.md` with the document that would supply it."),
"",
paste0("Experimenter metadata resolves ", mv("organism_part"), ", ", mv("strain_or_breed"),
       ", ", mv("labeling_strategy"), " and sequential Lys-C plus trypsin digestion. The full",
       " 2024-11-28 protocol is preserved in `metadata/sample_preparation_protocol.tsv` and",
       " `pride/SAMPLE_PREPARATION_PROTOCOL.md`."),
"",
"No missing acquisition or search metadata was inferred. Instrument model, acquisition mode,",
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
"editor_source_data/  Supplementary Table 2 (Supplementary_Table_2_proteomics.xlsx), its",
"                  sheet inventory, and the figure source data",
"pride/            SDRF and deposition notes",
"provenance/       lineage, versions, parameters, figure map, revision records, manifest",
"                  and checksums",
"```",
"",
"## 19. Column definitions",
"",
"`metadata/data_dictionary.tsv` describes every column of every table in this package. It",
"is generated by walking the tables, and the build fails if a column has no definition.",
"Supplementary Table 2 (`editor_source_data/Supplementary_Table_2_proteomics.xlsx`) defines",
"the columns it shows, under their display headers, on its Column definitions sheet.",
"",
"## 20. What is not here",
"",
"Native raw mass-spectrometry data. A file named `processed_protein_group_matrix_raw.txt`",
"is **not** raw mass-spectrometry data: it is a 5,747-row protein-group matrix from the",
"search software. See `pride/README_PRIDE.md`.",
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
