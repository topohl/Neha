#!/usr/bin/env Rscript

# Publication release, stage 08: the changelog for editors and reviewers.
#
# Produces
#   provenance/REVISION_PROTEOMICS_DATA_CHANGELOG.md
#   provenance/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv     classified fold-change mentions
#   provenance/MANUSCRIPT_TERMINOLOGY_ACTIONS.md     hand edits outside the repo
#
# Counts for the revised analysis are read from the release tables (stages 04-06), the
# figure rebuild's audit and source-data files, the inferential checks and the pooled run's
# tables. A few values from the original submission and earlier audits (PCA variance, EWCE
# cell types, plate QC, UMAP outlier ranks) are typed in. The changelog says which findings
# weakened or changed in the revised analysis.
#
# Panels use the manuscript numbering (Fig. 2 X, Supp. Fig. 4 X), as do the released
# source-data files; internal keys and rebuild file names (Figure3, Fig3G_*) keep the old
# names.

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
figure_map <- read_release(rel("provenance", "figure_source_map.tsv"))
secondary <- read_release(rel("metadata", "secondary_analysis_manifest.tsv"))
corrections <- read_release(rel("metadata", "sample_class_corrections.tsv"))
differential <- read_release(rel("differential_analysis",
                                "primary_differential_proteins.tsv.gz"))
gsea_cc <- read_release(rel("enrichment", "primary_GSEA_GO_CC.tsv.gz"))
gsea_kegg <- read_release(rel("enrichment", "primary_GSEA_KEGG.tsv.gz"))
ora_go <- read_release(rel("enrichment", "primary_ORA_GO_BP.tsv.gz"))
ora_cc <- read_release(rel("enrichment", "primary_ORA_GO_CC.tsv.gz"))
inv <- RELEASE_DESIGN_INVARIANTS
# Stops unless the run manifest describes RELEASE_PRIMARY_MODEL on the locked GCTs.
stat_run <- release_read_stat_run(DATA_ROOT)

# Post hoc Pairing x CNO interaction and leave-one-animal-out check, from the inferential
# checks (not among the 12 primary contrasts; stage 02 checks they used this statistics
# run). Their counts go into the manuscript action list.
INFERENTIAL <- release_inferential_root(DATA_ROOT)
interaction_protein_path <- file.path(INFERENTIAL, "pairing_cno_interaction_protein_summary.csv")
interaction_gsea_path <- file.path(INFERENTIAL, "pairing_cno_interaction_gsea_summary.csv")
interaction_diag_path <- file.path(INFERENTIAL, "pairing_cno_interaction_diagnostics.csv")
loo_path <- file.path(INFERENTIAL, "robustness", "mcherry_learning_leave_one_out.csv")
inferential_report <- file.path(INFERENTIAL, "interaction_and_robustness_report.txt")
for (p in c(interaction_protein_path, interaction_gsea_path, interaction_diag_path, loo_path,
            inferential_report)) {
  release_assert_exists(p, basename(p))
}
interaction_proteins <- release_read_csv(interaction_protein_path)
interaction_gsea <- release_read_csv(interaction_gsea_path)
interaction_diag <- release_read_csv(interaction_diag_path)
loo <- release_read_csv(loo_path)
loo_full <- loo[loo$animal_dropped == "none", , drop = FALSE]
loo_drop <- loo[loo$animal_dropped != "none", , drop = FALSE]
if (nrow(loo_full) != 1L || !nrow(loo_drop)) {
  stop("Leave-one-animal-out table lacks the full fit or the refits: ", loo_path, call. = FALSE)
}

# Experimenter-supplied Methods facts, from the checked-in records the release publishes.
experimenter_metadata <- release_read_experimenter_metadata(REPO_ROOT)
sample_preparation_protocol <- release_read_sample_preparation_protocol(REPO_ROOT)
crosswalk <- release_old_package_crosswalk()

REVISION_ROOT <- file.path(DATA_ROOT, "03_output", "reviewer_revision_animal_level_20260827")
retained_report <- file.path(REVISION_ROOT, "manifests",
                             "RETAINED_PANELS_3G_3H_3I_SuppE_REPORT.md")
final_report <- file.path(REVISION_ROOT, "manifests", "FINAL_REPORT.md")

# Fig. 2 D, F, G, H, I and Supp. Fig. 4 F were rebuilt on the per-sample-class statistics
# (2026-11-25, on the nPermSimple = 100000 GO enrichment). Their numbers are read from the
# rebuild's audit and source-data files; its run script is the record of the rebuild.
REBUILD_ROOT <- release_figure_rebuild_root(DATA_ROOT)
REBUILD_REL <- paste(RELEASE_STATISTICS_RUN$figure_rebuild_relative, collapse = "/")
REBUILD_FIG <- file.path(REBUILD_ROOT, "figure_panels_style_matched", "Figure3")
rebuild_record <- release_figure_rebuild_record(DATA_ROOT)
fig_f_src_path <- file.path(REBUILD_FIG,
  "Fig3F_memory_engram_GSEA_animal_level_STYLE_MATCHED_source_data.csv")
fig_g_src_path <- file.path(REBUILD_FIG, "Fig3G_source_data.csv")
fig_h_audit_path <- file.path(REBUILD_FIG, "Fig3H_interpretation_audit.csv")
fig_i_audit_path <- file.path(REBUILD_FIG, "Fig3I_interpretation_audit.csv")
for (p in c(retained_report, final_report, rebuild_record, fig_f_src_path, fig_g_src_path,
            fig_h_audit_path, fig_i_audit_path)) {
  release_assert_exists(p, basename(p))
}
audit_value <- function(path, item) {
  a <- release_read_csv(path)
  v <- as.character(a$value[a$item == item])
  if (length(v) != 1L || !nzchar(v)) stop("No single '", item, "' in ", path, call. = FALSE)
  v
}
fig_h_r <- as.numeric(audit_value(fig_h_audit_path, "corrected_r"))
fig_h_n <- as.integer(audit_value(fig_h_audit_path, "n_proteins"))
fig_h_control_r <- as.numeric(audit_value(fig_h_audit_path, "control_r_learning_signature"))
fig_i_n <- as.integer(audit_value(fig_i_audit_path, "n_overlapping_terms"))
fig_i_opposed <- as.integer(audit_value(fig_i_audit_path, "n_opposed"))
fig_i_concordant <- as.integer(audit_value(fig_i_audit_path, "n_concordant"))
fig_g_src <- release_read_csv(fig_g_src_path)
fig_f_src <- release_read_csv(fig_f_src_path)
g_bp <- fig_g_src[fig_g_src$ontology == "BP", , drop = FALSE]
g_cc <- fig_g_src[fig_g_src$ontology == "CC", , drop = FALSE]
if (!nrow(g_bp) || !nrow(g_cc) || length(unique(g_bp$n_FDR_significant)) != 1L ||
    length(unique(g_cc$n_FDR_significant)) != 1L) {
  stop("Fig. 2 G source data no longer holds one GO-BP and one GO-CC block.", call. = FALSE)
}
f_terms <- unique(fig_f_src$ID[as.logical(fig_f_src$selected_by_column)])
f_mcherry <- sum(fig_f_src$sample_class == "mcherry" & as.logical(fig_f_src$fdr_significant) &
                   as.logical(fig_f_src$selected_by_column))

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
gsea_cc_fdr_by_comparison <- vapply(
  unique(gsea_cc$canonical_comparison),
  function(cmp) fdr_terms(gsea_cc, gsea_cc$canonical_comparison == cmp), integer(1))
kegg_fdr_by_comparison <- vapply(
  unique(gsea_kegg$canonical_comparison),
  function(cmp) fdr_terms(gsea_kegg, gsea_kegg$canonical_comparison == cmp), integer(1))
term_count <- function(counts, cmp) {
  v <- unname(counts[cmp])
  v[is.na(v)] <- 0L
  as.integer(v)
}
plural <- function(n, one, many) if (n == 1L) one else many

MCH_LEARNING <- "mcherry_paired_veh_over_mcherry_unpaired_veh"
mch <- diff_summary[diff_summary$canonical_comparison == MCH_LEARNING, , drop = FALSE]
ora_up <- ora_go[ora_go$canonical_comparison == MCH_LEARNING &
                   ora_go$query_list == "fdr_significant_higher_in_numerator", , drop = FALSE]
ora_up_padj <- as.numeric(ora_up$adjusted_p_value)
ora_cc_up <- ora_cc[ora_cc$canonical_comparison == MCH_LEARNING &
                      ora_cc$query_list == "fdr_significant_higher_in_numerator", , drop = FALSE]
if (nrow(mch) != 1L || !nrow(ora_up) || !nrow(ora_cc_up)) {
  stop("The mCherry learning comparison or its ORA up-lists are missing from the release.",
       call. = FALSE)
}
ora_cc_up_padj <- as.numeric(ora_cc_up$adjusted_p_value)
ora_cc_up_sig <- sum(ora_cc_up_padj < 0.05, na.rm = TRUE)
# The released name of the Fig. 2 G redundant-terms file (the panel's second file).
g_skip_release <- basename(strsplit(
  figure_map$panel_source_data_release_path[figure_map$figure == "Figure3" &
                                              figure_map$panel == "G"],
  ";", fixed = TRUE)[[1]][[2]])
comparison_labels <- release_comparison_label(diff_summary$sample_class,
                                              diff_summary$numerator_condition,
                                              diff_summary$denominator_condition)
n_prot_mapped <- as.integer(diff_summary$n_significant_fdr_0_05_mapped_only)
n_bp_now <- term_count(gsea_fdr_by_comparison, diff_summary$canonical_comparison)
n_cc_now <- term_count(gsea_cc_fdr_by_comparison, diff_summary$canonical_comparison)

# AveExpr as the per-sample-class model reports it, measured on the released table.
ave <- release_aveexpr_summary(differential)

# Leave-one-animal-out summary of the mCherry learning comparison (read, not refitted).
loo_range <- function(col) range(as.numeric(loo_drop[[col]]))
loo_fdr <- loo_range("n_fdr05")
loo_r_effect_min <- min(as.numeric(loo_drop$r_effect))
loo_bp <- loo_range("go_bp_fdr05"); loo_bp_kept <- loo_range("go_bp_retained_of_full")
loo_cc <- loo_range("go_cc_fdr05"); loo_cc_kept <- loo_range("go_cc_retained_of_full")
loo_df <- unique(as.integer(loo_drop$df_residual))

# Gene sets fgsea could not test, per block, and the permutation settings behind every GSEA
# count (GO at nPermSimple 100000; KEGG and the sensitivity analyses at the fgsea default).
contrast_manifest <- read_release(rel("metadata", "primary_contrast_manifest.tsv"))
untested <- release_fgsea_untested(diff_summary$canonical_comparison, DATA_ROOT)
permutations <- release_enrichment_permutations(DATA_ROOT)
untested_n <- function(cmp, ontology) {
  vapply(cmp, function(x) {
    as.integer(sum(untested$n_sets[untested$canonical_comparison == x &
                                     untested$ontology == ontology]))
  }, integer(1), USE.NAMES = FALSE)
}
untested_blocks <- function(families, ontologies = "GO_BP", classes = NULL) {
  release_untested_blocks_text(untested, contrast_manifest, families, ontologies, classes)
}
PERMUTATION_NOTE <- paste(
  permutations$sentence,
  "The gene sets fgsea could not test are counted per comparison in section 5 of the",
  "changelog and listed in enrichment/gsea_untested_gene_sets.tsv.")

# The superseded pooled model (one fit over all 48 units), read for the "what changed from
# the pooled release" sections. Its folders are not locked, so they are checked first: both
# indexes must name the pooled statistics GCT, each mapped table must match the hash in its
# index, and term counts are compared with the pooled enrichment index.
POOLED_STATS <- release_rel_path(RELEASE_STATISTICS_RUN$pooled_stats_relative, DATA_ROOT)
POOLED_MAPPED <- file.path(POOLED_STATS, "mapped", "forward")
POOLED_ENRICH_ROOT <- release_rel_path(RELEASE_STATISTICS_RUN$pooled_enrichment_relative,
                                       DATA_ROOT)
POOLED_ENRICH <- file.path(POOLED_ENRICH_ROOT, "per_comparison")
pooled_mapped_index_path <- file.path(POOLED_STATS, "mapped", "indexMappedComparisons.csv")
pooled_enrich_index_path <- file.path(POOLED_ENRICH_ROOT, "indexEnrichmentComparisons.csv")
pooled_mapped_index <- release_read_csv(pooled_mapped_index_path)
pooled_enrich_index <- release_read_csv(pooled_enrich_index_path)
for (ix in list(list(d = pooled_mapped_index, p = pooled_mapped_index_path),
                list(d = pooled_enrich_index, p = pooled_enrich_index_path))) {
  if (nrow(ix$d) != inv$n_primary_contrasts ||
      !all(as.character(ix$d$source_gct_sha256) == RELEASE_STATISTICS_RUN$pooled_stat_gct_sha256)) {
    stop("The pooled-model index does not record the pooled statistics GCT (",
         substr(RELEASE_STATISTICS_RUN$pooled_stat_gct_sha256, 1, 16), "): ", ix$p,
         call. = FALSE)
  }
}
pooled_mapped <- function(cmp) {
  p <- file.path(POOLED_MAPPED, paste0(sub("_over_", "_vs_", cmp), ".csv"))
  release_assert_exists(p, basename(p))
  row <- pooled_mapped_index[pooled_mapped_index$canonical_comparison == cmp, , drop = FALSE]
  if (nrow(row) != 1L || !identical(release_sha256(p), as.character(row$mapped_output_sha256))) {
    stop("Pooled mapped table does not hash to the value its index recorded: ", p, call. = FALSE)
  }
  release_read_csv(p)
}
# Row counts are checked against what the pooled run recorded, and FDR counts too where the
# run recorded one (GO-BP; its index has no KEGG FDR column).
POOLED_INDEX_COUNT <- list(GSEA_GO_BP.csv = c(rows = "go_gsea_term_count",
                                              fdr = "go_gsea_fdr_term_count"),
                           GSEA_KEGG.csv = c(rows = "kegg_gsea_term_count"))
pooled_terms <- function(cmp, file) {
  p <- file.path(POOLED_ENRICH, cmp, file)
  release_assert_exists(p, basename(p))
  d <- release_read_csv(p)
  n <- sum(as.numeric(d[["p.adjust"]]) < 0.05, na.rm = TRUE)
  ix_row <- pooled_enrich_index[pooled_enrich_index$canonical_comparison == cmp, , drop = FALSE]
  cols <- POOLED_INDEX_COUNT[[file]]
  ok <- nrow(ix_row) == 1L && identical(as.integer(ix_row[[cols[["rows"]]]]), nrow(d)) &&
    (!"fdr" %in% names(cols) || identical(as.integer(ix_row[[cols[["fdr"]]]]), as.integer(n)))
  if (!isTRUE(ok)) {
    stop("Pooled enrichment table differs from the pooled run's record: ", p, call. = FALSE)
  }
  n
}
# Rows (protein groups) are counted; symbols are kept only to name short lists.
pooled_sig <- lapply(diff_summary$canonical_comparison, function(cmp) {
  d <- pooled_mapped(cmp)
  as.character(d$mapped_gene_symbol[as.logical(d$significant) %in% TRUE])
})
pooled_n <- lengths(pooled_sig)
pooled_vs_now <- data.frame(
  Comparison = comparison_labels,
  `Proteins (mapped), pooled` = pooled_n,
  `Proteins (mapped), now` = n_prot_mapped,
  `GO-BP terms, pooled` = vapply(diff_summary$canonical_comparison, pooled_terms, integer(1),
                                 file = "GSEA_GO_BP.csv", USE.NAMES = FALSE),
  `GO-BP terms, now` = n_bp_now,
  `KEGG terms, pooled` = vapply(diff_summary$canonical_comparison, pooled_terms, integer(1),
                                file = "GSEA_KEGG.csv", USE.NAMES = FALSE),
  `KEGG terms, now` = term_count(kegg_fdr_by_comparison, diff_summary$canonical_comparison),
  check.names = FALSE, stringsAsFactors = FALSE)
lost <- which(pooled_vs_now[["Proteins (mapped), pooled"]] > 0L & n_prot_mapped == 0L)
lost_text <- if (length(lost)) paste(vapply(lost, function(i) {
  n <- pooled_n[i]
  paste0(comparison_labels[i], " (", n,
         if (n <= 10L) paste0(": ", paste(sort(unique(pooled_sig[[i]])), collapse = ", ")) else "",
         ")")
}, character(1)), collapse = "; ") else ""

#' Render a data frame as a GitHub-flavoured markdown table: cells trimmed, numeric
#' columns right-aligned.
md_table <- function(df, align = NULL) {
  cols <- names(df)
  numeric_col <- vapply(df, function(v) is.numeric(v) || is.integer(v), logical(1))
  if (is.null(align)) align <- ifelse(numeric_col, "--:", "---")
  cells <- as.data.frame(lapply(df, function(v) trimws(as.character(v))),
                         stringsAsFactors = FALSE, check.names = FALSE)
  body <- vapply(seq_len(nrow(cells)), function(i) {
    paste0("| ", paste(unlist(cells[i, ], use.names = FALSE), collapse = " | "), " |")
  }, character(1))
  c(paste0("| ", paste(cols, collapse = " | "), " |"),
    paste0("|", paste(align, collapse = "|"), "|"),
    body)
}

# Tables in the document use display labels ("c-Fos+: paired CNO vs paired VEH"); the
# underscored keys stay in the data files.
da_tbl <- data.frame(
  Comparison = release_comparison_label(diff_summary$sample_class,
                                        diff_summary$numerator_condition,
                                        diff_summary$denominator_condition),
  `Proteins tested` = diff_summary$n_proteins_tested,
  `FDR < 0.05` = diff_summary$n_significant_fdr_0_05,
  `Tested, mapped` = diff_summary$n_proteins_tested_mapped_only,
  `FDR < 0.05, mapped` = diff_summary$n_significant_fdr_0_05_mapped_only,
  `Higher in numerator` = diff_summary$n_significant_higher_in_numerator_mapped_only,
  `Higher in denominator` = diff_summary$n_significant_higher_in_denominator_mapped_only,
  check.names = FALSE, stringsAsFactors = FALSE)

gsea_context <- diff_summary[match(names(gsea_fdr_by_comparison),
                                   diff_summary$canonical_comparison), , drop = FALSE]
gsea_tbl <- data.frame(
  Comparison = release_comparison_label(gsea_context$sample_class,
                                        gsea_context$numerator_condition,
                                        gsea_context$denominator_condition),
  `GO-BP, FDR < 0.05 (moderated-t-ranked)` = unname(gsea_fdr_by_comparison),
  `GO-CC, FDR < 0.05 (moderated-t-ranked)` =
    term_count(gsea_cc_fdr_by_comparison, names(gsea_fdr_by_comparison)),
  `KEGG, FDR < 0.05 (moderated-t-ranked)` =
    term_count(kegg_fdr_by_comparison, names(gsea_fdr_by_comparison)),
  `GO-BP, FDR < 0.05 (effect-size-ranked sensitivity)` =
    term_count(sens_fdr_by_comparison, names(gsea_fdr_by_comparison)),
  `Sets fgsea could not test (GO-BP / GO-CC / KEGG)` = paste(
    untested_n(names(gsea_fdr_by_comparison), "GO_BP"),
    untested_n(names(gsea_fdr_by_comparison), "GO_CC"),
    untested_n(names(gsea_fdr_by_comparison), "KEGG"), sep = " / "),
  check.names = FALSE, stringsAsFactors = FALSE)

# "Interpretation changed" has two bases. The panels rebuilt on the per-sample-class
# statistics are judged against the pooled-model release of 2026-09-07; every other panel
# keeps its 2026-08-27 judgment against the original submission. The table says which.
panel_tbl_rows <- figure_map[figure_map$primary_or_secondary != "not_applicable", , drop = FALSE]
rebuilt_panel <- grepl("^REBUILT_WITHIN_COMPARTMENT", panel_tbl_rows$revision_status)
panel_tbl <- data.frame(
  Panel = panel_tbl_rows$panel_label,
  Status = panel_tbl_rows$revision_status,
  `Unit: original -> revised` = paste(panel_tbl_rows$original_statistical_unit, "->",
                                      panel_tbl_rows$statistical_unit),
  `Displayed statistics changed` = panel_tbl_rows$displayed_statistics_changed,
  `Interpretation changed` = panel_tbl_rows$biological_interpretation_changed,
  `Compared with` = ifelse(rebuilt_panel, "pooled release (2026-09-07)",
                           "original submission"),
  check.names = FALSE, stringsAsFactors = FALSE)

changed_interpretation <- panel_tbl_rows$panel_label[
  toupper(panel_tbl_rows$biological_interpretation_changed) == "YES"]
rebuilt_labels <- panel_tbl_rows$panel_label[rebuilt_panel]
fig_label <- function(figure, panel) {
  v <- figure_map$panel_label[figure_map$figure == figure & figure_map$panel == panel]
  if (length(v) != 1L) stop("No single panel label for ", figure, " ", panel, call. = FALSE)
  v
}
LBL <- list(C = fig_label("Figure3", "C"), D = fig_label("Figure3", "D"),
            E = fig_label("Figure3", "E"), F = fig_label("Figure3", "F"),
            G = fig_label("Figure3", "G"), H = fig_label("Figure3", "H"),
            I = fig_label("Figure3", "I"),
            SA = fig_label("Supplementary_proteomics", "A"),
            SB1 = fig_label("Supplementary_proteomics", "B1"),
            SB3 = fig_label("Supplementary_proteomics", "B3"),
            SC = fig_label("Supplementary_proteomics", "C"),
            SD = fig_label("Supplementary_proteomics", "D"),
            SE = fig_label("Supplementary_proteomics", "E"),
            SF = fig_label("Supplementary_proteomics", "F"))
supp_prefix <- sub(" A$", "", LBL$SA)   # e.g. "Supp. Fig. 4"

# The sample classes Supp. Fig. 4 F shows, from its released source data (it omits one).
sf_row <- figure_map[figure_map$figure == "Supplementary_proteomics" & figure_map$panel == "F", ,
                     drop = FALSE]
sf_classes <- unique(as.character(release_read_csv(
  strsplit(sf_row$panel_source_data_canonical_path, ";", fixed = TRUE)[[1]][[1]])$sample_class))
# Untested blocks for a panel, as text; "none" when fgsea tested every set.
panel_untested <- function(...) {
  x <- untested_blocks(...)
  if (nzchar(x)) paste(x, "sets") else "none"
}
# Legend guidance for a GSEA panel whose blocks have untested gene sets.
untested_legend <- function(...) {
  x <- untested_blocks(...)
  if (!nzchar(x)) return("")
  paste0(" Say that fgsea could not test some gene sets in these blocks (", x, " sets), so",
         " they are absent rather than non-significant.")
}

# --------------------------------------------------------------------------------------
# document
# --------------------------------------------------------------------------------------

L <- c(
"# Proteomics data and analysis changelog: original submission versus revised analysis",
"",
"This document accompanies the revised proteomics data package. It states what changed",
"between the analysis supplied with the original submission and the analysis supplied now,",
"including the findings that weakened or no longer hold.",
"",
"Protein and term counts for the revised analysis are read from the release tables in",
"this package.",
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
paste0("| Statistical model | two-sample moderated t (ProTigy), hemisphere level | **limma, one",
       " model per sample class: ", stat_run$n_columns_in_model, " animal-level units each,",
       " residual df ", stat_run$df_residual, "** |"),
"",
"Two hemispheres taken from the same animal are not independent samples. Treating them as",
"independent replicates doubles the apparent sample size and shrinks standard errors, which",
"inflates significance throughout. This is pseudoreplication, and it affected every",
"inferential number in the original submission.",
"",
"The correction is an aggregation, not a reanalysis of a different experiment: left and",
"right values are averaged with equal weight within each AnimalID x sample_class unit, on",
"the existing processed values (standardized per protein, see section 12), with no further",
"transformation. In a balanced design the",
"mean of animal means equals the mean of all hemispheres, so **point estimates barely move**.",
"What changes is the **inference**: standard errors, degrees of freedom, and therefore every",
"p-value.",
"",
"## 2. Data files: 96 measurement columns versus 48 animal-level columns",
"",
"| Layer | Rows x columns | Role |",
"|---|---|---|",
"| Search output (`pg.matrix_raw.txt`) | 5,747 protein groups x 96 acquisitions | quantification output, pre-filter |",
"| Processed measurement-level matrix | 5,349 x 96 acquisitions | filtered, imputed and standardized per protein; not log2 abundance (section 12) |",
"| **Animal-level matrix** | **5,349 x 48 animal units** | **the matrix all inference is computed on**: mean of the two hemispheres on the same per-protein standardized scale |",
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
"were drawn from. The 22-row difference is non-mouse entries (human keratins and other skin",
"proteins, digestion enzymes and bovine serum proteins), which are tested but not",
"annotatable.",
"",
paste0("The statistics come from limma ", stat_run$limma_version, " (R ", stat_run$r_version,
       ", `", RELEASE_PRIMARY_MODEL$script, "`), fitted separately within each sample class:",
       " one linear model ~0 + condition over that class's ", stat_run$n_columns_in_model,
       " animal-level units (n = ", inv$n_animals_per_stratum, " animals per condition,",
       " residual df ", stat_run$df_residual, "), empirical Bayes moderation (robust, no",
       " intensity trend) and Benjamini-Hochberg adjustment within each comparison. The prior",
       " variance differs between sample classes (",
       paste(vapply(unique(diff_summary$sample_class), function(k) {
         r <- stat_run$manifest[startsWith(stat_run$manifest$canonical_comparison,
                                           paste0(k, "_")), , drop = FALSE][1, ]
         sprintf("%s %.2f", release_sample_class_display(k), as.numeric(r$s2_prior))
       }, character(1)), collapse = ", "),
       "), so each class has its own variance estimate rather than sharing one. Section 7b",
       " compares this with the pooled model used by the release of 2026-09-07."),
"",
md_table(da_tbl),
"",
sprintf("**%d of the 12 primary comparisons %s no FDR-significant protein at all.** Only %d %s any.",
        n_without_hits, plural(n_without_hits, "contains", "contain"),
        n_with_hits, plural(n_with_hits, "contains", "contain")),
"This is a substantive change from the original submission, in which hemisphere-level",
"inference produced significant sets far more widely.",
"",
paste0("The mapped proteins at FDR < 0.05 are listed in",
       " `differential_analysis/significant_proteins.tsv` (Supplementary Table 2d); the ",
       sum(diff_summary$n_significant_fdr_0_05) -
         sum(diff_summary$n_significant_fdr_0_05_mapped_only),
       " non-mouse groups that also pass FDR < 0.05 are in",
       " `primary_differential_proteins.tsv.gz` only."),
"",
"## 5. Enrichment analysis",
"",
"| | Original submission | Revised analysis |",
"|---|---|---|",
"| Input statistics | hemisphere-level | animal-level, one limma model per sample class |",
"| **GSEA ranking statistic** | **the effect size stored as `logFC`** | **moderated t statistic** |",
"| GSEA / ORA cutoffs | filtered at source | `pvalue_cutoff = 1`, `qvalue_cutoff = 1` -- the complete tested result is released |",
"| Multiple testing | BH | BH, FDR 0.05 |",
"| Gene set size | -- | minGSSize 10, maxGSSize 800 |",
"| ORA universe | not documented | all unique successfully mapped measured UniProt accessions, per comparison |",
"| Gene-set collections | not documented | GO-BP, GO-CC and KEGG (KEGG downloaded at run time, not version-pinned) |",
"| GSEA permutations | not documented | fgsea, one recorded seed per comparison (DOSE called with `seed = FALSE`) |",
paste0("| fgsea nPermSimple | not documented | ", permutations$setting, " |"),
"",
paste0("GO-CC is a primary output again. ", RELEASE_GO_CC_PROVENANCE_NOTE, " ",
       RELEASE_ENRICHMENT_SEED_HANDLING),
"",
paste0("**Permutations and untested gene sets.** ", permutations$sentence, " When a",
       " comparison's moderated t values are strongly one-sided, fgsea cannot estimate a P",
       " value for some gene sets and DOSE drops them, so they are missing from the GSEA",
       " tables rather than non-significant. The last column of the table below counts them;",
       " the GO sets are listed in `enrichment/gsea_untested_gene_sets.tsv`.",
       " Blocks behind the figure panels: ", LBL$F, " (",
       panel_untested("paired_veh_vs_unpaired_veh"), "); ", LBL$I, " (",
       panel_untested(c("paired_veh_vs_unpaired_veh", "paired_cno_vs_paired_veh")), "); ",
       LBL$SF, " (", panel_untested("paired_cno_vs_paired_veh", "GO_BP", sf_classes), "); ",
       LBL$G, " (",
       panel_untested("paired_veh_vs_unpaired_veh", c("GO_BP", "GO_CC"), "mcherry"), ")."),
"",
"Ranking a GSEA by the effect size alone ignores the variance of the estimate: a large",
"difference measured with n = 3 and enormous uncertainty outranks a smaller, well-estimated",
"one. The canonical analysis now ranks by the moderated t statistic. **Effect-size-ranked",
"results are retained only as a sensitivity analysis**, in a separate file, with",
"`analysis_role = sensitivity` on every row. They must not be reported as the primary",
"enrichment result. (The file keeps the name `GSEA_log2FC_sensitivity` because that is what",
"the canonical run called the variant; on the units of that statistic see section 12.)",
"",
"FDR-significant terms per comparison (GO-BP, GO-CC and KEGG canonical; GO-BP sensitivity):",
"",
md_table(gsea_tbl),
"",
paste("The FDR-significant gene sets of the three canonical GSEA tables are listed, with their",
      "leading edges as gene symbols, in `enrichment/significant_gene_sets.tsv.gz`",
      "(Supplementary Table 2f)."),
"",
paste0("**KEGG identifiers.** KEGG gene sets are keyed by NCBI Gene IDs, so the leading-edge",
       " members (`core_enrichment`) of the KEGG rows in `enrichment/primary_GSEA_KEGG.tsv.gz`",
       " and `enrichment/GSEA_log2FC_sensitivity.tsv.gz` are NCBI Gene IDs, labeled",
       " `gene_identifier_type = ENTREZID`; GO rows hold UniProt accessions (`UNIPROT`). The",
       " first per-compartment release (folder `",
       utils::tail(RELEASE_STATISTICS_RUN$superseded_release_relative, 1),
       "`) labeled the KEGG rows UNIPROT in error; the members themselves are unchanged."),
"",
"## 6. Figure panels",
"",
md_table(panel_tbl),
"",
sprintf(paste("Panels whose biological interpretation changed, each against the comparison",
              "named in the last column: **%s**. %s were rebuilt on the per-sample-class",
              "statistics and are judged against the pooled-model release of 2026-09-07; the",
              "other panels keep their judgment against the original submission."),
        paste(changed_interpretation, collapse = ", "),
        paste(rebuilt_labels, collapse = ", ")),
"",
"### Findings that were weakened, eliminated or changed",
"",
paste0("- **", LBL$G, " (GO enrichment, mCherry learning).** The original panel showed GO",
       " over-representation terms (synaptic vesicle endocytosis, postsynaptic density",
       " organization, respiratory electron transport chain) over a hemisphere-level",
       " significant-protein list. The direct animal-level analogue (ORA over the",
       " FDR-significant proteins higher in paired_veh) finds **",
       sum(ora_up_padj < 0.05, na.rm = TRUE), " of ", format(nrow(ora_up), big.mark = ","),
       " GO-BP terms at FDR < 0.05 (minimum adjusted P = ",
       sprintf("%.2f", min(ora_up_padj, na.rm = TRUE)), ") and ", ora_cc_up_sig, " of ",
       format(nrow(ora_cc_up), big.mark = ","), " GO-CC terms (minimum adjusted P = ",
       release_format_p(min(ora_cc_up_padj, na.rm = TRUE)), ")**. The panel now shows",
       " moderated-t-ranked GSEA for the same comparison, not over-representation: ",
       nrow(g_bp), " GO-BP and ", nrow(g_cc), " GO-CC terms, chosen from the ",
       unique(g_bp$n_FDR_significant), " GO-BP and ", unique(g_cc$n_FDR_significant),
       " GO-CC terms at FDR < 0.05 with redundant terms skipped (listed in `",
       g_skip_release, "`). These are pathway-level results of a comparison whose origin is",
       " unresolved (section 7)."),
paste0("- **", LBL$F, " (GSEA bubble grid, learning).** Ranked by the moderated t at animal",
       " level, the grid shows ", length(f_terms), " terms; the mCherry column contributes ",
       f_mcherry, " FDR-significant GO-BP ", plural(f_mcherry, "term", "terms"), "."),
paste0("- **", LBL$SC, " (EWCE cell-type identity).** 12 of the legacy 14 cell types",
       " survive. `Int12` and `S1PyrL23` no longer clear global FDR."),
paste0("- **", LBL$H, " and ", LBL$I, ".** Both remain, but neither can be read as evidence",
       " that chemogenetic inhibition reverses a learning-associated proteomic program. Both",
       " contrasts share `paired_veh` on opposite sides, so anything high in `paired_veh` is",
       " pushed up on one axis and down on the other **by construction**. ", LBL$H, " gives",
       " r = ", sprintf("%.2f", fig_h_r), " over ", fig_h_n, " proteins; the control with the",
       " shared arm on the same side gives r = ", sprintf("%+.2f", fig_h_control_r), " on the",
       " same protein set. In ", LBL$I, ", ", fig_i_opposed, " of the ", fig_i_n,
       " overlapping terms are opposed in sign (", fig_i_concordant, " concordant), as the",
       " shared arm biases them to be. ", LBL$I, " has been retitled neutrally."),
paste0("- **", LBL$SE, ".** The original four cross-compartment comparisons were invalid",
       " twice over: hemispheres were treated as independent, and both arms of every",
       " comparison came from the same three animals. It has been refitted as a secondary",
       " within-animal paired model and is no longer presented as one of the designed",
       " contrasts."),
"",
"### Findings that did not change",
"",
paste0("- **", LBL$E, " and ", LBL$SD, " (rank abundance).** Descriptive; no statistical",
       " test. Regenerated at animal level and numerically indistinguishable from the",
       " originals (Spearman rank r = 0.99999; Pearson on group means r = 1.000), as expected",
       " in a balanced design."),
paste0("- **", supp_prefix, " A, B2, B3.** Acquisition-level technical QC. The acquisition is",
       " the correct unit for identification depth and for hemisphere / collection-plate QC, so",
       " these panels are unchanged and remain valid."),
"",
"### Numbers that must be retyped in the figures",
"",
"| Panel | Original | Revised |",
"|---|---|---|",
paste0("| ", LBL$C, " | PC1 25.4% / PC2 17.7% / PC3 10.2% | PC1 29.1% / PC2 20.5% / PC3 12.1% |"),
paste0("| ", LBL$SB1, " | PC1 25.4% / PC2 17.7% | PC1 29.1% / PC2 20.5% |"),
paste0("| ", LBL$D, " | hemisphere-level significant set | ",
       mch$n_significant_fdr_0_05_mapped_only, " FDR-significant proteins (",
       mch$n_significant_higher_in_numerator_mapped_only, " higher in paired_veh, ",
       mch$n_significant_higher_in_denominator_mapped_only, " higher in unpaired_veh) |"),
paste0("| ", LBL$H, " | r = -0.71 | r = ", sprintf("%.2f", fig_h_r),
       " (learning signature, n = ", fig_h_n, ") |"),
"",
"## 7. The unresolved mCherry result",
"",
sprintf(paste("The mCherry paired-vehicle versus unpaired-vehicle comparison retains %d",
              "FDR-significant proteins of %d mapped (%d higher in paired_veh, %d higher in",
              "unpaired_veh). The statistical result stands as computed. What it reflects does",
              "not follow from it:"),
        as.integer(mch$n_significant_fdr_0_05_mapped_only),
        as.integer(mch$n_proteins_tested_mapped_only),
        as.integer(mch$n_significant_higher_in_numerator_mapped_only),
        as.integer(mch$n_significant_higher_in_denominator_mapped_only)),
"",
paste0("- the shift is broad and asymmetric, and it is accompanied by pathway-level",
       " enrichment (GSEA: ", term_count(gsea_fdr_by_comparison, MCH_LEARNING), " GO-BP, ",
       term_count(gsea_cc_fdr_by_comparison, MCH_LEARNING), " GO-CC and ",
       term_count(kegg_fdr_by_comparison, MCH_LEARNING), " KEGG terms at FDR < 0.05; ORA of",
       " the higher-in-paired_veh list: ",
       sum(as.numeric(ora_cc_up$adjusted_p_value) < 0.05, na.rm = TRUE), " GO-CC and ",
       sum(ora_up_padj < 0.05, na.rm = TRUE), " GO-BP terms). Pathway enrichment does not by",
       " itself distinguish a biological from a technical origin;"),
paste0("- the protein count depends on single animals. Leaving one animal out (",
       paste(loo_df, collapse = "/"), " residual df) gives between ", loo_fdr[1], " and ",
       format(loo_fdr[2], big.mark = ","), " FDR-significant proteins (full fit ",
       format(as.integer(loo_full$n_fdr05), big.mark = ","), "), and ", loo_bp[1], " to ",
       loo_bp[2], " GO-BP terms keeping ", loo_bp_kept[1], " to ", loo_bp_kept[2], " of the ",
       as.integer(loo_full$go_bp_fdr05), "; GO-CC is more stable (", loo_cc[1], " to ",
       loo_cc[2], " terms keeping ", loo_cc_kept[1], " to ", loo_cc_kept[2], " of ",
       as.integer(loo_full$go_cc_fdr05), "). The effect estimates themselves are stable",
       " (correlation with the full fit at least ", sprintf("%.2f", loo_r_effect_min), "). ",
       permutations$sensitivity_sentence),
"- pairing condition is completely associated with collection plate in this comparison;",
paste("-", release_plate_injection_sentence()),
"- no technical batch effect has been demonstrated.",
"",
"Neither a biological nor a technical explanation is established.",
"",
"## 7b. What changed from the pooled release (2026-09-07)",
"",
paste0("The release of 2026-09-07 used the ", RELEASE_PRIMARY_MODEL$superseded_label,
       ". The four sample classes therefore shared one variance estimate per protein. The",
       " statistics are now",
       " fitted within each sample class (section 4). Effect sizes are unchanged to",
       " floating-point rounding; moderated t, P values, adjusted P values and AveExpr change.",
       " FDR < 0.05 counts, pooled model and now:"),
"",
md_table(pooled_vs_now),
"",
if (length(lost)) paste0("No longer protein-level results: ", lost_text, ".") else
  "No comparison lost all of its FDR-significant proteins.",
paste0("Panels rebuilt on the per-sample-class statistics: ",
       paste(rebuilt_labels, collapse = ", "), " (`", REBUILD_REL, "`)."),
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
paste0("onto**. ", release_plate_injection_sentence(), " ", LBL$SB3, " has been relabeled",
       " from `Plate` to `Collection plate`; collection plate explains 0.01% of PC1 and 1.5%",
       " of PC2 at acquisition level, neither significant."),
"",
"## 10. What is superseded",
"",
"| Original submission file | Replacement in this package | Superseded? |",
"|---|---|---|",
"| `processed_protein_group_matrix_raw.txt` (5,747 x 96) | `processed_data/protein_feature_annotation.tsv.gz` (annotation) and, for abundance, `processed_data/protein_abundance_measurement_level.tsv.gz` | Partly. Retained as the search-output layer in the lineage; it is **not** raw MS data. |",
"| `processed_protein_group_matrix_filtered_umap_adjusted.xlsx` (5,349 x 96) | `processed_data/protein_abundance_measurement_level.tsv.gz` | No, but demoted: it is no longer the matrix inference is drawn from. |",
"| `GSEA_ORA_all_results.xlsx` | `enrichment/*.tsv.gz` and `editor_source_data/Supplementary_Table_2_proteomics.xlsx` (significant gene sets: table 2f) | **Yes, entirely.** Hemisphere-level inference, effect-size-ranked GSEA, incomplete contrast coverage. |",
"| `README.txt` (group codes only) | `README_DATA.md`, `metadata/data_dictionary.tsv`, this changelog | **Yes.** |",
paste0("| Release of 2026-09-07 (`", paste(RELEASE_STATISTICS_RUN$pooled_release_relative,
                                         collapse = "/"),
       "`, pooled model) | this package | **Yes.** Kept frozen for comparison; see",
       " section 7b. |"),
paste0("| First per-compartment release (folder `",
       utils::tail(RELEASE_STATISTICS_RUN$superseded_release_relative, 1),
       "`, deployed 2026-11-26; old workbook layout) | this package: same statistics (every",
       " number identical); workbook renamed",
       " `editor_source_data/Supplementary_Table_2_proteomics.xlsx` with lettered tables 2a-2q",
       " (the effect-size-ranked GSEA, the absolute-effect-size ORA lists, the superseded and",
       " excluded analyses, the metadata field status, the full software list and the",
       " artwork-only panels are now in the data files only); new",
       " files `differential_analysis/significant_proteins.tsv` and",
       " `enrichment/significant_gene_sets.tsv.gz`; KEGG leading-edge members now labeled",
       " ENTREZID (NCBI Gene IDs; previously labeled UNIPROT in error) in",
       " `enrichment/primary_GSEA_KEGG.tsv.gz` and `enrichment/GSEA_log2FC_sensitivity.tsv.gz`;",
       " `enrichment/enrichment_coverage.tsv` adds the `gsea_seed` column; `effect_size_units`,",
       " `effect_size_definition` and other descriptive text in US spelling, including two field",
       " names (`metadata/metadata_field_provenance.tsv` 'proteomics labeling',",
       " `provenance/analysis_parameters.tsv` 'normalization'); figure source data now CSV with",
       " reader headers and display labels (Supp. Fig. 4 F blank cells now",
       " `SuppFig4F_blank_cells.csv`); changelog, terminology files and figure source map moved",
       " to `provenance/`; collection plate described as an injection block | **Yes.** |"),
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
paste0("**", crosswalk$statement, "**"),
"",
"| submitted file | previously described properties | crosswalk basis |",
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
       " and all 5,349 analyzed protein identifiers are a strict subset of the 5,747;"),
paste0("- for `GSEA_ORA_all_results.xlsx`, an identification **by analysis generation**.",
       " The workbook was never located or opened, so no claim is made about its internal",
       " structure. Its supersession does not require one: any hemisphere-level,",
       " effect-size-ranked GSEA/ORA output is superseded by the animal-level,",
       " moderated-t-ranked canonical run, whatever that workbook contained. Its described",
       " contrast coverage is likewise carried from the prior description of the",
       " submission, not from inspection here."),
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
              "differs from every pre-correction record. An audit established what",
              "that difference is: a **deliberate correction** applied during historical",
              "sample-identity quality control. It is not a mislabeling, and it is not an",
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
       " nearest-class-center suggestion agrees with the applied correction for five of the",
       " six, and differs for one (N60: nearest center `cfos`, applied `neuron`). That",
       " column is published as `umap_nearest_class_suggestion` so the difference is",
       " visible. The applied classes are those the validated",
       " analysis used, supported by the independent profile-similarity and left/right-pair",
       " analyses."),
"",
"### What is not claimed",
"",
paste0("**No surviving prose note records the rationale.** No memo, comment or README in the",
       " project explains why the reassignment was made. The correction is still resolved:",
       " the correction table is the preserved record, and the evidence above establishes",
       " the assignment."),
"",
paste0("**Both labels are retained.** `metadata/sample_metadata.tsv` carries",
       " `original_sample_class` and `analysis_sample_class` for all 96 measurements, with",
       " `sample_class_corrected` marking the six. The `sample_class` column remains the",
       " canonical analysis-time assignment (the one the locked GCTs and every validated",
       " statistic were computed on), and the data dictionary says so. The",
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
paste0("The abundance matrix all statistics were computed on is **standardized separately",
       " for each protein** across the measurement-level dataset. Measured directly on the",
       " measurement-level matrix: row means are approximately 0 and row standard deviations",
       " approximately 1. The animal-level matrix is the mean of the two hemispheres on that",
       " scale, so its row means are the same and its row standard deviations are below 1.",
       " The statistical GCT nevertheless stores the coefficient as `logFC` (limma's name for",
       " it), and the canonical split tables carry it as `log2fc`."),
"",
paste0("A difference of group means on that scale is a **", RELEASE_EFFECT_SIZE$public_term,
       "** -- equivalently, a ", RELEASE_EFFECT_SIZE$public_term_alt, ". It is not a log2",
       " ratio, and a value of 2 does not mean a four-fold change."),
"",
paste0("It is *not* called a \"standardized mean difference\" without",
       " qualification: that term invites a Cohen's d reading, i.e. scaling by a pooled",
       " within-group standard deviation. The scaling here is a per-protein standard",
       " deviation taken across the whole measurement dataset before any grouping, which is",
       " a different quantity."),
"",
paste0("It is also *not* called \"z-scored\". The released values carry two decimals, so on",
       " those bytes no protein is exactly mean 0 / standard deviation 1, and the operation",
       " that produced the matrix is unresolved -- no script, no parameters (see",
       " `provenance/UPSTREAM_PREPROCESSING_GAP.md`). \"Standardized\" is what the numbers",
       " support, so that is what the package says."),
"",
paste0("Corroborating evidence, from the canonical tables themselves: `AveExpr`, the",
       " average abundance limma reports for each protein, is the mean over the ",
       stat_run$n_columns_in_model, " animal-level units of one sample class, so it is not",
       " centered within a class (it ranges from ", sprintf("%.2f", ave$range_within_classes[1]),
       " to ", sprintf("%.2f", ave$range_within_classes[2]), "). Averaged over the four",
       " sample classes it stays close to zero: the largest absolute value over all ",
       format(ave$n_proteins, big.mark = ","), " proteins is ",
       format(signif(ave$max_abs_mean_over_classes, 3), scientific = FALSE),
       " (13_validate_release.R checks this). That can only happen if the input was centered",
       " per protein."),
"",
paste0("This package therefore publishes the quantity as `", RELEASE_EFFECT_SIZE$public_field,
       "`, carrying `effect_size_definition`, `effect_size_source_column` and",
       " `source_statistic_field = \"", RELEASE_EFFECT_SIZE$source_field, "\"` on every row,",
       " rather than republishing it under a name that asserts a fold-change reading."),
"",
"**No value was changed.** The numbers are bit-identical to the canonical tables (the",
"build re-reads every split table and requires exact equality), and directions, rankings,",
"p-values and every downstream result are unaffected. Only the label is corrected.",
"",
paste0("The same applies to the GSEA sensitivity analysis. Canonical GSEA ranks by the",
       " moderated t and is unchanged. What was described as a \"log2FC-ranked\" sensitivity",
       " analysis is publicly an **", RELEASE_EFFECT_SIZE$sensitivity_public_label, "**;",
       " its filenames and its `rank_statistic` value keep the historical `log2fc` token so",
       " the export can still be matched to the canonical run. No ranking value changed, and",
       " the relabeling reran no enrichment."),
"",
"### Methods wording",
"",
paste0("> ", RELEASE_EFFECT_SIZE$methods_sentence),
"",
paste0("Every occurrence of the historical fold-change wording in the active",
       " publication-facing files was",
       " classified and is itemized in",
       " `provenance/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv`. Figure axes, captions and",
       " manuscript text live outside this repository and must be corrected by hand; the",
       " specific items are listed in",
       " `provenance/MANUSCRIPT_TERMINOLOGY_ACTIONS.md`."),
"",
"## 13. What this package does not contain",
"",
"Native raw mass-spectrometry acquisition files. The 96 `.d` acquisition directories are",
"identified by name in `metadata/sample_metadata.tsv` but are not included here; they are",
"intended for a separate PRIDE deposition. Experimenter-supplied tissue, animal-cohort,",
"label-free and sample-preparation facts are propagated from checked-in metadata records.",
"Instrument model, acquisition mode, search settings, AnimalID-level sex and age at tissue",
"collection remain unavailable rather than inferred, and are itemized in",
"`pride/SDRF_MISSING_METADATA.md`.",
"",
"---",
"",
paste0("Underlying revision records: `", basename(final_report), "` and `",
       basename(retained_report), "` (2026-08-27; the panels that do not use the differential",
       " statistics) and `", REBUILD_REL, "/",
       paste(RELEASE_STATISTICS_RUN$figure_rebuild_record, collapse = "/"), "` (",
       paste(rebuilt_labels, collapse = ", "), ", rebuilt on the per-sample-class",
       " statistics). Interaction and leave-one-animal-out checks: `",
       paste(RELEASE_STATISTICS_RUN$inferential_relative, collapse = "/"), "/",
       basename(inferential_report), "`."),
paste0("Generated by `", STAGE, "`.")
)

path <- release_path("provenance", "REVISION_PROTEOMICS_DATA_CHANGELOG.md")
release_write_lines(L, path)
# The pooled-model indexes are inputs too: sections 7b and 4b are counted from the tables
# they index, after the checks above.
changelog_inputs <- c(retained_report, final_report, rebuild_record, fig_f_src_path,
                      fig_g_src_path, fig_h_audit_path, fig_i_audit_path, loo_path,
                      inferential_report, stat_run$path, pooled_mapped_index_path,
                      pooled_enrich_index_path,
                      release_locked_artefact_path("enrichment_index_go_bp", DATA_ROOT),
                      release_locked_artefact_path("enrichment_index_go_cc", DATA_ROOT),
                      release_locked_artefact_path("enrichment_index_kegg", DATA_ROOT))
release_register("provenance/REVISION_PROTEOMICS_DATA_CHANGELOG.md",
                 "original submission versus revised analysis, for the editor and reviewers",
                 c(rel("differential_analysis", "primary_differential_summary.tsv"),
                   rel("provenance", "figure_source_map.tsv"), changelog_inputs),
                 c(NA_character_, NA_character_,
                   vapply(changelog_inputs, release_sha256, character(1), USE.NAMES = FALSE)),
                 STAGE, "md")
release_log("  wrote REVISION_PROTEOMICS_DATA_CHANGELOG.md (", length(L), " lines)")
release_log("  ", n_without_hits, " of 12 primary comparisons have no FDR-significant protein")
release_log("  panels with changed interpretation: ",
            paste(changed_interpretation, collapse = ", "))

# --------------------------------------------------------------------------------------
# effect-size terminology audit
# --------------------------------------------------------------------------------------
# Every fold-change mention in the release, classified, in two parts:
#   - `applied`, listed by hand: the wording corrections made in this layer (they can no
#     longer be found by scanning the corrected files), the old tokens that are kept, and
#     the fold changes of other statistics;
#   - a scan of the release tables built so far; each hit becomes a RESIDUAL
#     PUBLIC_MISLABEL row, and 13_validate_release.R fails on it.
#
# Status values:
#   INTERNAL_PROVENANCE        canonical file name, column value or source-field name, kept
#                              so it still matches the canonical run
#   PUBLIC_MISLABEL            reader-facing text that called the effect size a fold change
#   LEGITIMATE_OTHER_CONTEXT   a real fold change of another statistic (EWCE), or a frozen
#                              script that is not edited

audit_row <- function(file, location, old_term, new_term, status, notes) {
  data.frame(file = file, location = location, old_term = old_term, new_term = new_term,
             status = status, notes = notes,
             stringsAsFactors = FALSE, check.names = FALSE)
}

LAYER_REL <- "07_publication_release"

applied <- rbind(
  audit_row(paste0(LAYER_REL, "/04_build_differential_results.R"),
            "EFFECT_SIZE_DEFINITION", "standardised mean difference",
            RELEASE_EFFECT_SIZE$public_term, "PUBLIC_MISLABEL",
            paste("Unqualified 'standardized mean difference' invites a Cohen's d reading",
                  "(pooled within-group SD). Defined once in RELEASE_EFFECT_SIZE.")),
  audit_row(paste0(LAYER_REL, "/04_build_differential_results.R"),
            "EFFECT_SIZE_DEFINITION", "z-scored per protein (mean 0, sd 1)",
            "standardized separately for each protein (row means ~0, row SDs ~1)",
            "PUBLIC_MISLABEL",
            paste("Exact z-scoring is not established: values carry 2 decimals and the",
                  "producing step is unresolved. Wording avoids claiming exact z-scoring.")),
  audit_row(paste0(LAYER_REL, "/03_build_processed_data_exports.R"),
            "VALUE_SCALE", "per-protein standardised (z-scored) abundance",
            "per-protein standardized abundance, with measured deviations reported",
            "PUBLIC_MISLABEL",
            paste("Exact z-scoring is not established; the assertion reports the measured",
                  "row means and SDs.")),
  audit_row(paste0(LAYER_REL, "/02_build_contrast_manifest.R"),
            "INTERPRETATION (3 contrast families)", "Positive log2FC = higher in ...",
            "A positive effect_size_sd_units means higher in ...", "PUBLIC_MISLABEL",
            "Reader-facing interpretation of the sign."),
  audit_row(paste0(LAYER_REL, "/02_build_contrast_manifest.R"),
            "secondary manifest, gsea_log2fc_rank_sensitivity label",
            "GSEA ranked by log2FC (sensitivity analysis)",
            paste("GSEA ranked by the standardized-abundance effect size",
                  "(effect-size-ranked sensitivity analysis)"),
            "PUBLIC_MISLABEL",
            "The manifest KEY and the ranking_statistic value keep the log2fc token."),
  audit_row(paste0(LAYER_REL, "/05_build_enrichment_exports.R"),
            "sensitivity export description and log lines", "ranked by log2FC",
            RELEASE_EFFECT_SIZE$sensitivity_public_label, "PUBLIC_MISLABEL",
            "Filenames and rank_statistic value unchanged; see INTERNAL_PROVENANCE rows."),
  audit_row(paste0(LAYER_REL, "/07_build_editor_source_workbook.R"),
            "README sheet, items 'Effect size' and 'Complete tables and keys'",
            "log2FC-ranked", "effect-size-ranked", "PUBLIC_MISLABEL",
            "Supplementary Table 2 README sheet."),
  audit_row(paste0(LAYER_REL, "/R/release_presentation.R"),
            "workbook display labels (ranking statistic, ORA query list)",
            "log2fc / top_absolute_log2fc", "effect size (SD units) / |effect size| > 1 SD",
            "PUBLIC_MISLABEL",
            paste("Cell values shown in Supplementary Table 2; the .tsv files keep the",
                  "tokens (INTERNAL_PROVENANCE rows below).")),
  audit_row(paste0(LAYER_REL, "/06_build_figure_source_data.R"),
            "figure_source_map, Figure3|H description",
            "scatter of learning log2FC against CNO log2FC",
            "scatter of the learning effect size against the CNO effect size (SD units)",
            "PUBLIC_MISLABEL", "Panel description shown to the editor."),
  audit_row(paste0(LAYER_REL, "/06_build_figure_source_data.R"),
            paste("figure_source_map, Figure3|H panel_title and statistical_unit (also",
                  "Supplementary Table 2p and changelog section 6)"),
            paste("Learning logFC vs CNO logFC scatter, mCherry; animal on both axes, both",
                  "moderated limma log2FC"),
            paste("Learning effect size vs CNO effect size scatter, mCherry; animal on both",
                  "axes; standardized abundance differences (SD units) from the",
                  "per-sample-class limma fit"),
            "PUBLIC_MISLABEL",
            paste("Was copied from the 2026-08-27 regenerated_panel_manifest.csv for the",
                  "rebuilt panel. limma's logFC is not moderated and here is not a log2",
                  "ratio. The six rebuilt panels now take title and unit from stage 06.")),
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
            paste("Unchanged. The name matches the canonical run's",
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
            paste("Recorded so the published value can be traced to its source.",
                  "The exported effect_size_sd_units equals this source field exactly.")),
  audit_row(paste(c(RELEASE_STATISTICS_RUN$stats_relative, "split", "forward", "*.csv"),
                  collapse = "/"),
            "canonical column log2fc", "log2fc", "log2fc", "INTERNAL_PROVENANCE",
            "Canonical analysis artifact, outside this layer. Never modified."),
  audit_row(paste("06_manuscript_figure_revision/*.R,",
                  "06_manuscript_figure_revision/within_compartment_20261124/*.R and",
                  "06_manuscript_figure_revision/within_compartment_20261125/*.R"),
            "frozen figure scripts", "log2FC", "log2FC", "LEGITIMATE_OTHER_CONTEXT",
            paste("Frozen; not edited. The 2026-11-25 rebuild draws the",
                  "released Fig. 2 D and H axes in SD units; the log2FC that remains there is",
                  "in comments, file names and internal column names. The earlier scripts'",
                  "rendered axis labels and captions are listed in",
                  "MANUSCRIPT_TERMINOLOGY_ACTIONS.md for manual correction.")),
  audit_row("enrichment/primary_EWCE.tsv.gz", "fold_change column",
            "fold change", "fold change", "LEGITIMATE_OTHER_CONTEXT",
            paste("EWCE's own statistic: the ratio of observed to bootstrap-null mean",
                  "expression proportion: a real fold change, unrelated to the",
                  "differential effect size.")),
  audit_row("04_differential_expression_enrichment/01_clusterProfiler.r",
            "canonical enrichment script", "log2fc", "log2fc", "INTERNAL_PROVENANCE",
            paste("Produced the canonical enrichment. Not a publication-facing document; its",
                  "log2fc tokens are kept because renaming them would break the link to the",
                  "outputs. The GO runs record its nPermSimple, serial BiocParallel backend",
                  "and untested-set audit; the KEGG run predates them. The script in the",
                  "repository is later than both runs only in its defaults (no default",
                  "output root; per-compartment inputs).")),
  audit_row(RELEASE_PRIMARY_MODEL$script,
            "canonical statistics script", "logFC", "logFC", "INTERNAL_PROVENANCE",
            paste("Writes the per-sample-class limma statistics in the ProTigy GCT layout.",
                  "logFC is limma's coefficient name and is kept so the split, mapping and",
                  "enrichment stages run unchanged; the published field is",
                  RELEASE_EFFECT_SIZE$public_field))
)

# Scan the built tables, not the builder sources (those hold the claim list and the audit
# rows themselves). Documents written later in the build are scanned by
# 13_validate_release.R.
scan_targets <- c(
  rel("differential_analysis", "primary_differential_summary.tsv"),
  rel("metadata", "primary_contrast_manifest.tsv"),
  rel("metadata", "secondary_analysis_manifest.tsv"),
  rel("provenance", "figure_source_map.tsv"),
  rel("enrichment", "enrichment_run_parameters.tsv"),
  rel("enrichment", "enrichment_coverage.tsv")
)
# Label a scanned path for the audit: release-relative for built files, repo-relative
# for sources, so the audit reads the same whichever scratch root the build used.
audit_scan_label <- function(p) {
  p <- gsub("\\\\", "/", p)
  for (root in c(OUT_ROOT, REPO_ROOT)) {
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
                          release_path("provenance",
                                       "EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv"))
release_register("provenance/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv",
                 paste("every fold-change mention in the active publication-facing surface,",
                       "classified as INTERNAL_PROVENANCE, PUBLIC_MISLABEL or",
                       "LEGITIMATE_OTHER_CONTEXT"),
                 "publication release layer sources and built documents",
                 NA_character_, STAGE, "tsv")
release_log("  wrote EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv (", wa$rows, "x", wa$cols, ")")

# --------------------------------------------------------------------------------------
# manuscript action note -- what has to be fixed by hand, outside this repository
# --------------------------------------------------------------------------------------

# Panels kept although their result is structural or unresolved; each needs its limitation
# in the legend. All numbers are read from the release or the rebuild files.
MA_PANEL_CAVEATS <- list(
  list(panel = LBL$F,
       legend = paste0("State that the grid shows GSEA ranked by the moderated t for the",
                       " learning comparison in each sample class.",
                       untested_legend("paired_veh_vs_unpaired_veh"))),
  list(panel = LBL$G,
       legend = paste0("State that the panel shows GSEA ranked by the moderated t for mCherry",
                       " paired-VEH versus unpaired-VEH, and that the ", nrow(g_bp), " GO-BP",
                       " and ", nrow(g_cc), " GO-CC terms plotted were chosen from the ",
                       unique(g_bp$n_FDR_significant), " GO-BP and ",
                       unique(g_cc$n_FDR_significant), " GO-CC terms at FDR < 0.05, redundant",
                       " terms skipped. Do not call it over-representation, because it is",
                       " GSEA. If ORA of the same up-list is mentioned, give both ontologies: ",
                       sum(ora_up_padj < 0.05, na.rm = TRUE), " of ",
                       format(nrow(ora_up), big.mark = ","), " GO-BP terms (minimum adjusted P ",
                       sprintf("%.2f", min(ora_up_padj, na.rm = TRUE)), ") and ",
                       ora_cc_up_sig, " of ", format(nrow(ora_cc_up), big.mark = ","),
                       " GO-CC terms at FDR < 0.05. Say that the",
                       " comparison's origin is unresolved because pairing is completely",
                       " associated with collection plate.",
                       untested_legend("paired_veh_vs_unpaired_veh", c("GO_BP", "GO_CC"),
                                       "mcherry"))),
  list(panel = LBL$H,
       legend = paste0("Descriptive only. Report r as approximately ", sprintf("%.2f", fig_h_r),
                       " (n = ", fig_h_n, " proteins) and state that the two contrasts share",
                       " the paired-vehicle condition on opposite sides, so a negative",
                       " correlation is expected by construction; the control with the shared",
                       " arm on the same side gives r = ", sprintf("%+.2f", fig_h_control_r),
                       " on the same protein set. The correlation must not be interpreted",
                       " causally or as independent evidence of reversal.")),
  list(panel = LBL$I,
       legend = paste0("State that ", fig_i_opposed, " of the ", fig_i_n, " overlapping",
                       " FDR-significant pathways are opposed in direction (", fig_i_concordant,
                       " concordant), and that opposition is the expected structural",
                       " consequence of the two contrasts sharing the paired-vehicle reference",
                       " on opposite sides. Describe the pattern neutrally; it is not evidence",
                       " that CNO reverses a pairing-associated program.",
                       untested_legend(c("paired_veh_vs_unpaired_veh",
                                         "paired_cno_vs_paired_veh")))),
  list(panel = LBL$SF,
       legend = paste0("State that the grid shows GSEA ranked by the moderated t for paired",
                       " CNO versus paired VEH in each sample class shown.",
                       untested_legend("paired_cno_vs_paired_veh", "GO_BP", sf_classes))),
  list(panel = LBL$SE,
       legend = paste("Identify it as a secondary, within-animal paired cross-compartment",
                      "model, not one of the primary contrasts and not a treatment effect.",
                      "State n = 3 animals and that inference is correspondingly limited.")),
  list(panel = paste(supp_prefix, "A, B2, B3"),
       legend = paste("These are acquisition-level technical QC panels. The acquisition is",
                      "the correct unit here and the bars/points are not biological",
                      "replicates; say so, and label the plate axis \"collection plate\"."))
)

# Enrichment software versions, from the same shared builder as provenance/software_versions.
sw_tbl <- release_build_software_versions(DATA_ROOT, REPO_ROOT)
enr <- sw_tbl[sw_tbl$applies_to_stage == "differential enrichment (GSEA / ORA)", , drop = FALSE]
pkg_versions_enrich <- stats::setNames(enr$version, enr$component)
if (!all(c("clusterProfiler", "fgsea") %in% names(pkg_versions_enrich))) {
  stop("The enrichment run audit records no clusterProfiler or fgsea version.", call. = FALSE)
}

# Comparisons by what the release supports, for the checklist and section 4.
cls_protein <- comparison_labels[n_prot_mapped > 0L]
cls_pathway <- comparison_labels[n_prot_mapped == 0L & n_bp_now + n_cc_now > 0L]
cls_null <- comparison_labels[n_prot_mapped == 0L & n_bp_now + n_cc_now == 0L]
or_none <- function(x) if (length(x)) paste(x, collapse = "; ") else "none"

# Established Methods facts, from the checked-in experimenter records.
MA_VERIFIED_METHODS <- {
  mv <- function(id) release_metadata_value(experimenter_metadata, id)
  prot <- sample_preparation_protocol
  c(
    paste0("Organism and cohort: ", mv("organism"), ", ", mv("strain_or_breed"),
           ", supplied by ", mv("animal_supplier"), " (", mv("animal_supplier_location"),
           "); ", mv("cohort_sex_composition"), "; age at experiment start ",
           mv("age_at_experiment_start"), "."),
    paste0("Tissue and capture: ", mv("organism_part"), ", ", mv("collection_method"), "."),
    paste0("Sample classes: mCherry-positive, cFos-positive, NeuN-positive/cFos-negative/",
           "mCherry-negative neurons, and a neuropil compartment negative for all three.",
           " The neuropil class is a tissue compartment, NOT a cell type."),
    paste0("Quantification strategy: ", mv("labeling_strategy"), "."),
    paste0("Sample preparation: ",
           paste(paste0(prot$step_name, " (", prot$details, ")"), collapse = "; "), "."),
    paste0("Statistical software: limma ", stat_run$limma_version, " in R ",
           stat_run$r_version, " (`", RELEASE_PRIMARY_MODEL$script, "`) for the canonical",
           " animal-level differential statistics, one linear model per sample class;",
           " clusterProfiler ", pkg_versions_enrich[["clusterProfiler"]], " with fgsea ",
           pkg_versions_enrich[["fgsea"]], " for enrichment.")
  )
}

MA_UNRESOLVED_METHODS <- c(
  "the exact developmental stage at tissue collection;",
  "the exact age at tissue collection (only age at experiment start is recorded);",
  "per-animal sex assignments (only the cohort composition is recorded);",
  "the MS instrument model -- the run-name token `Olive` is a local alias, not a model;",
  "the acquisition mode;",
  paste("the search / quantification software and its version -- do NOT name a tool.",
        "In particular, do not infer one from `pg.matrix` filenames, table schemas or",
        "historical QC wording;"),
  "search modification parameters, precursor tolerance and fragment tolerance."
)

MA_REVIEWER_CHECKLIST <- c(
  "acknowledges the animal-level correction where the comment concerns sample size, replication or significance;",
  "does not defend hemisphere-level independence anywhere;",
  "states what changed after the corrected analysis, including findings that weakened or became null;",
  paste0("distinguishes protein-level results (", or_none(cls_protein),
         ") from pathway-level-only observations (", or_none(cls_pathway),
         ") and null results (", or_none(cls_null), ");"),
  "quotes no number from the earlier pooled-model release (section 4b);",
  "does not call the effect size a log2 fold change or a fold change;",
  "describes the canonical GSEA as ranked by the moderated t statistic;",
  "uses the collection-plate identifiability wording and claims no demonstrated batch effect;",
  "describes the six sample-class corrections as resolved intentional metadata corrections, not as swaps, discrepancies or open questions;",
  "does not imply the previously submitted files were directly inspected;",
  "quotes the same n, the same significant proteins and the same limitations as the manuscript."
)

MA <- c(
  "# Manuscript revision map and required actions",
  "",
  paste0("Generated by `", STAGE, "` on ", release_timestamp_utc(), "."),
  "",
  "Everything in the released data package has been corrected. The items below live",
  "**outside this repository** -- in the manuscript file, the figure files, the reviewer",
  "response and the submission forms -- and must be changed by hand. Nothing here changes",
  "a number: every count quoted below is read from the released tables at build time.",
  "",
  "The manuscript and the reviewer response are **not** in this repository, so neither has",
  "been read or audited. This document is therefore a revision map to apply to them, not a",
  "record of corrections already made to them. Section 9 lists what has to be supplied",
  "before that audit can actually be done.",
  "",
  "| # | Category | What it covers |",
  "|---|---|---|",
  "| 1 | EFFECT_SIZE_TERMINOLOGY, FIGURE_AXIS | the effect size is not a log2 fold change |",
  "| 2 | SAMPLE_CLASS_CORRECTION | the six corrected acquisitions |",
  "| 3 | INFERENTIAL_UNIT, HEMISPHERE_AVERAGING, SAMPLE_SIZE | animal as the unit, n = 3 |",
  "| 4 | RESULT_SIGNIFICANCE | which claims the animal-level analysis supports; 4b what changed from the pooled release |",
  "| 5 | GSEA_RANKING | canonical ranking is the moderated t |",
  "| 6 | COLLECTION_PLATE | identifiability, not a demonstrated batch effect |",
  "| 7 | FIGURE_LEGEND | panels needing an explicit caveat |",
  "| 8 | METHODS_METADATA, UNSUPPORTED_SOFTWARE, OLD_SUBMISSION_CROSSWALK | what may and may not be stated |",
  "| 9 | REVIEWER_RESPONSE + missing documents | consistency checklist and what is still needed |",
  "",
  "## 1. The effect size is not a log2 fold change",
  "",
  paste0("Published field: `", RELEASE_EFFECT_SIZE$public_field, "`. Preferred wording: **",
         RELEASE_EFFECT_SIZE$public_term, "**, or \"",
         RELEASE_EFFECT_SIZE$public_term_alt, "\"."),
  "",
  "| Where | Incorrect wording | Correct wording |",
  "|---|---|---|",
  paste0("| ", LBL$D, " volcano, x-axis | incorrect: `log2FC` / \"log2 fold change\" | \"",
         RELEASE_EFFECT_SIZE$public_term, "\" |"),
  paste0("| ", LBL$H, " scatter, both axes | learning / CNO `log2FC` | learning / CNO effect",
         " size (SD units) |"),
  paste0("| ", LBL$SE, " panels derived from the paired model | `log2FC` | ",
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
         " \"standardized mean difference\" without qualification either -- that reads as",
         " Cohen's d, which uses a pooled within-group SD rather than the across-dataset",
         " per-protein SD used here."),
  "",
  paste0("The frozen scripts under `06_manuscript_figure_revision/` (and their copies of",
         " the per-sample-class rebuilds under `within_compartment_20261124/` and",
         " `within_compartment_20261125/`) were **not** edited: they are the record of how the",
         " figures were produced. The 2026-11-25 rebuild already labels the axes of ",
         LBL$D, " and ", LBL$H, " in SD units; panels kept from the earlier scripts still say",
         " `log2FC` on their axes and must be relabeled at the manuscript stage. The",
         " classification of every occurrence is in `EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv`."),
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
    "Do not describe these six as mislabeled samples, as an unresolved discrepancy, or as",
    "requiring reanalysis. No abundance measurement changed and no statistic was recomputed.",
    "",
    "Do not claim a written rationale survives. It does not; the correction table is the",
    "preserved record."
  ) else "No sample-class correction is recorded.",
  "",
  "## 3. The inferential unit, hemisphere averaging and sample size",
  "",
  "Any statement implying that hemispheres are independent replicates must be corrected.",
  "The two hemispheres of one animal are repeated measurements contributing to one",
  "animal-level value, not two inferential replicates.",
  "",
  "| Where | Incorrect wording | Correct wording |",
  "|---|---|---|",
  paste0("| Methods, design paragraph | \"n = 6 per group\", \"", inv$n_measurement_records,
         " samples\" as the inferential n | \"n = ", inv$n_animals_per_stratum,
         " animals per condition per sample class\" |"),
  "| Methods, statistics paragraph | hemispheres as replicates | animal as the independent experimental unit; see the sentence below |",
  "| Results, any n | hemisphere-level n | animal-level n |",
  paste0("| Figure legends quoting n | acquisition or hemisphere n | animal n, except the",
         " technical-QC panels (", supp_prefix, " A, B2, B3) where the acquisition IS the",
         " correct unit |"),
  "| Any mention of technical replicates | \"technical replicates\" for L/R | left and right hemisphere measurements of the same animal |",
  "",
  "### Methods sentences to insert verbatim",
  "",
  paste0("> Left and right hemisphere measurements were averaged within each AnimalID x",
         " sample-class combination before any statistical test, reducing ",
         inv$n_measurement_records, " acquisition-level measurements to ",
         inv$n_animal_level_units, " animal-level units (", inv$n_animals,
         " animals x ", inv$n_sample_classes, " sample classes). The animal was treated as",
         " the independent experimental unit throughout, giving n = ",
         inv$n_animals_per_stratum, " animals per condition per sample class across all ",
         inv$n_strata, " strata, with complete left/right pairing."),
  "",
  paste0("> ", RELEASE_PRIMARY_MODEL$methods_sentence, " Significance was assessed at",
         " FDR < 0.05. Software: limma ", stat_run$limma_version, " in R ", stat_run$r_version,
         "."),
  "",
  "## 4. Which result claims the animal-level analysis supports",
  "",
  paste0("Correcting the statistical unit did not leave the findings intact. Of the ",
         inv$n_primary_contrasts, " primary comparisons, ", n_with_hits, " ",
         plural(n_with_hits, "contains", "contain"), " any FDR-significant protein and ",
         n_without_hits, " ", plural(n_without_hits, "contains", "contain"), " none."),
  "",
  "Protein-level and pathway-level outcomes per comparison, generated from the release:",
  "",
  "| Comparison | Proteins, FDR < 0.05 | GO-BP terms, FDR < 0.05 | GO-CC terms, FDR < 0.05 | How to describe it |",
  "|---|--:|--:|--:|---|",
  vapply(seq_len(nrow(diff_summary)), function(i) {
    np <- n_prot_mapped[i]; ng <- n_bp_now[i]; nc <- n_cc_now[i]
    learning <- diff_summary$contrast_family[i] == "paired_veh_vs_unpaired_veh"
    how <- if (np > 100L && ng + nc == 0L) {
      "large protein-level burden with no coherent pathway support; interpret cautiously"
    } else if (np > 100L) {
      paste0("large protein-level burden with pathway-level support, but origin unresolved",
             if (learning) " (pairing completely associated with collection plate)" else "",
             "; interpret cautiously")
    } else if (np > 0L && ng + nc > 0L) {
      "coherent: protein-level and pathway-level agree"
    } else if (np == 0L && ng + nc > 0L) {
      "pathway-level enrichment ONLY -- do not report as significant proteins"
    } else if (np > 0L) {
      "protein-level only; no pathway support"
    } else {
      "null at both levels"
    }
    paste0("| ", comparison_labels[i], " | ", np, " | ", ng, " | ", nc, " | ", how, " |")
  }, character(1)),
  "",
  paste("GO term counts:", PERMUTATION_NOTE),
  "",
  "The specific corrections most likely to be needed:",
  "",
  paste0("- **Protein-level results.** ", length(cls_protein), " of ", inv$n_primary_contrasts,
         " comparisons ", plural(length(cls_protein), "has", "have"),
         " FDR-significant proteins: ",
         if (length(cls_protein)) paste0(cls_protein, " (", n_prot_mapped[n_prot_mapped > 0L],
                                         " mapped)", collapse = "; ") else "none",
         ". Do not name individual proteins as FDR-significant for any other comparison."),
  paste0("- **The mCherry learning result** (", mch$n_significant_fdr_0_05_mapped_only,
         " mapped proteins) has pathway-level support but an unresolved origin, and its",
         " protein count depends on single animals: leaving one animal out gives ",
         loo_fdr[1], " to ", format(loo_fdr[2], big.mark = ","), " FDR-significant proteins",
         " (effect estimates stay correlated with the full fit at r >= ",
         sprintf("%.2f", loo_r_effect_min), "). Report it with that caveat."),
  if (length(lost)) paste0("- **No longer protein-level results.** Under the earlier pooled",
         " model these had FDR-significant proteins and now have none: ", lost_text,
         ". Remove any such claim, and any protein named from them, from the manuscript and",
         " the response.") else character(0),
  paste0("- **Neuropil** has **no** FDR-significant individual proteins in any of its three",
         " comparisons, but pathway-level enrichment (up to ",
         max(n_bp_now[diff_summary$sample_class == "neuropil"]),
         " GO-BP terms). Describe this as pathway-level enrichment. Do **not** describe it",
         " as protein-level FDR significance."),
  {
    ip <- interaction_proteins
    ig <- interaction_gsea[match(ip$sample_class, interaction_gsea$sample_class), , drop = FALSE]
    idg <- interaction_diag[match(ip$sample_class, interaction_diag$sample_class), , drop = FALSE]
    paste0("- **The Pairing x CNO interaction** (post hoc, per-sample-class model, n = ",
           inv$n_animals_per_stratum, " animals per cell) has ",
           sum(as.integer(ip$n_fdr05)), " FDR-significant proteins in any sample class (",
           format(as.integer(ip$n_tested[1]), big.mark = ","), " tested in each), alongside ",
           paste0(vapply(ip$sample_class, release_sample_class_display, character(1)), " ",
                  as.integer(ig$go_bp_fdr05), collapse = ", "),
           " FDR-significant GO-BP terms. Within each sample class the gene-set results",
           " largely follow a shift of the whole interaction contrast (share of all tested",
           " GO-BP terms whose NES has the majority sign: ",
           paste(vapply(seq_len(nrow(idg)), function(j) {
             p <- as.numeric(idg$go_bp_pct_NES_positive[j])
             sprintf("%s %.0f%% %s", release_sample_class_display(idg$sample_class[j]),
                     max(p, 100 - p), if (p >= 50) "positive" else "negative")
           }, character(1)), collapse = ", "),
           "), and the contrast keeps ",
           sprintf("%.0f%%", 100 * as.numeric(idg$plate2_minus_plate1_weight[1])),
           " of any additive collection-plate difference. Describe the pathway observations",
           " as suggestive and post hoc, not as protein-level confirmation of an interaction,",
           " and name no protein as an interaction hit.")
  },
  paste0("- **cFos unpaired CNO versus unpaired vehicle** is essentially null across",
         " analysis layers (", n_prot_mapped[diff_summary$canonical_comparison ==
                                               "cfos_unpaired_cno_over_cfos_unpaired_veh"],
         " proteins, ", term_count(gsea_fdr_by_comparison,
                                   "cfos_unpaired_cno_over_cfos_unpaired_veh"),
         " GO-BP terms). Do not retain stronger historical claims for it."),
  paste0("- **Neuron** comparisons have ",
         sum(n_prot_mapped[diff_summary$sample_class == "neuron"]),
         " FDR-significant proteins in total and ",
         paste(n_bp_now[diff_summary$sample_class == "neuron"], collapse = " / "),
         " GO-BP terms (",
         paste(comparison_labels[diff_summary$sample_class == "neuron"], collapse = " / "),
         "); pathway-level only."),
  "",
  "### 4b. What changed from the pooled release (2026-09-07)",
  "",
  paste0("Text or figures drafted from the earlier release of 2026-09-07 carry numbers that",
         " must be replaced. That release used the ", RELEASE_PRIMARY_MODEL$superseded_label,
         ". FDR < 0.05 counts, pooled model and now:"),
  "",
  md_table(pooled_vs_now),
  "",
  paste("Term counts for now:", PERMUTATION_NOTE),
  "",
  paste0("Also replace: ", LBL$D, " count (now ", mch$n_significant_fdr_0_05_mapped_only,
         " proteins, ", mch$n_significant_higher_in_numerator_mapped_only, " up / ",
         mch$n_significant_higher_in_denominator_mapped_only, " down); ", LBL$H,
         " r and n (now ", sprintf("%.2f", fig_h_r), ", n = ", fig_h_n, "); ", LBL$I,
         " overlap (now ", fig_i_opposed, " of ", fig_i_n, " opposed, ", fig_i_concordant,
         " concordant); ", LBL$F, " (now ", length(f_terms), " terms, mCherry ", f_mcherry,
         "); ", LBL$G, " content (now GSEA, ", nrow(g_bp), " GO-BP + ", nrow(g_cc),
         " GO-CC terms); and any statement that the mCherry learning comparison has no",
         " pathway support."),
  "",
  "## 5. The canonical GSEA ranking statistic",
  "",
  paste0("The canonical GSEA ranks proteins by the **moderated t statistic**. Ranking by the",
         " standardized-abundance effect size exists only as a sensitivity analysis, is",
         " reported separately, and carries `analysis_role = sensitivity` on every row."),
  "",
  "| Where | Incorrect wording | Correct wording |",
  "|---|---|---|",
  paste0("| Methods, enrichment paragraph | incorrect: \"ranked by log2 fold change\" |",
         " \"ranked by the moderated t statistic\" |"),
  paste0("| Any reference to the sensitivity run | incorrect: presented as the primary",
         " result | \"", RELEASE_EFFECT_SIZE$sensitivity_public_label, "\" |"),
  "| GSEA figure legends | incorrect: log2FC-ranked | moderated-t-ranked |",
  "",
  "## 6. Collection plate",
  "",
  "Use the identifiability wording, and only that. Insert verbatim:",
  "",
  paste0("> ", RELEASE_PLATE_SENTENCE),
  "",
  paste0("> ", RELEASE_PLATE_NO_EFFECT_SENTENCE),
  "",
  paste0("Do **not** write any of: ", paste0("\"", FORBIDDEN_PLATE_WORDING, "\"",
                                             collapse = ", "), "."),
  paste0("Collection plate is an animal-level design property. ",
         release_plate_injection_sentence(), " No technical batch effect has been",
         " demonstrated."),
  "",
  "## 7. Figure legends needing an explicit caveat",
  "",
  "These panels are retained, and each needs its limitation stated in the legend rather",
  "than being removed or softened away.",
  "",
  "| Panel | What the legend must say |",
  "|---|---|",
  vapply(MA_PANEL_CAVEATS, function(pc) {
    paste0("| ", pc$panel, " | ", pc$legend, " |")
  }, character(1)),
  "",
  "## 8. Methods metadata, software and the old submission",
  "",
  "Verified and safe to state:",
  "",
  vapply(MA_VERIFIED_METHODS, function(x) paste0("- ", x), character(1)),
  "",
  "Not established -- must not be stated or inferred:",
  "",
  vapply(MA_UNRESOLVED_METHODS, function(x) paste0("- ", x), character(1)),
  "",
  paste0("On software: the canonical animal-level statistics were produced with limma ",
         stat_run$limma_version, " (R ", stat_run$r_version, "), one model per sample class,",
         " by `", RELEASE_PRIMARY_MODEL$script, "`. ProTigy ", PROTIGY_ANIMAL_LEVEL_VERSION,
         " produced the earlier pooled animal-level run (one model over all ",
         inv$n_animal_level_units, " units) that is superseded, and the hemisphere-level 2025",
         " runs used v1.1.8; report neither as the software behind the reported results. No",
         " search or quantification software may be named: the identity and version of the",
         " upstream tool could not be recovered from the retained project files."),
  "",
  paste0("On the previously submitted package: none of the four submitted files was located",
         " or read. Do not describe their internal structure. The crosswalk in",
         " `REVISION_PROTEOMICS_DATA_CHANGELOG.md` section 10 is based on dimensions and",
         " content lineage for the two matrices and on analysis generation for the",
         " enrichment workbook, and the reviewer response must not imply the originals were",
         " inspected."),
  "",
  "## 9. Reviewer response consistency, and what is still missing",
  "",
  "The manuscript and the reviewer response must agree exactly on: the inferential unit,",
  "n, which proteins are significant, which effects are pathway-level only, effect-size",
  "terminology, and the stated limitations. Checklist for each proteomics-related response:",
  "",
  vapply(MA_REVIEWER_CHECKLIST, function(x) paste0("- [ ] ", x), character(1)),
  "",
  "### Documents that must be supplied before a manuscript audit can be performed",
  "",
  "Neither of the following is present in this repository or the shared project tree, so",
  "neither has been audited:",
  "",
  "- the current manuscript file (Results, Methods, Discussion, figure legends);",
  "- the reviewer response / rebuttal document;",
  "- the supplementary figure legend file, if it is separate from the manuscript;",
  paste0("- the figure source files whose axis labels need relabeling (", LBL$D, ", ",
         LBL$H, ", ", LBL$SE, ")."),
  "",
  "## 10. Still needed for the PRIDE deposition",
  "",
  "Not a terminology item, but it blocks deposition and is not something this repository",
  paste0("can resolve. See `pride/SDRF_MISSING_METADATA.md` for the authoritative list."),
  ""
)
mp <- release_path("provenance", "MANUSCRIPT_TERMINOLOGY_ACTIONS.md")
release_write_lines(MA, mp)
actions_inputs <- c(stat_run$path, pooled_mapped_index_path, pooled_enrich_index_path,
                    interaction_protein_path, interaction_gsea_path, interaction_diag_path,
                    loo_path, fig_g_src_path, fig_h_audit_path, fig_i_audit_path)
release_register("provenance/MANUSCRIPT_TERMINOLOGY_ACTIONS.md",
                 paste("figure axes, captions and manuscript text that must be corrected by",
                       "hand outside this repository"),
                 c(rel("differential_analysis", "primary_differential_summary.tsv"),
                   rel("provenance", "figure_source_map.tsv"), actions_inputs),
                 c(NA_character_, NA_character_,
                   vapply(actions_inputs, release_sha256, character(1), USE.NAMES = FALSE)),
                 STAGE, "md")
release_log("  wrote MANUSCRIPT_TERMINOLOGY_ACTIONS.md (", length(MA), " lines)")
release_log("stage 08 complete")
