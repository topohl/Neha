#!/usr/bin/env Rscript

# Publication release, stage 05 -- normalised enrichment and cell-type exports.
#
# Produces
#   enrichment/primary_GSEA_GO_BP.tsv.gz        canonical, ranked by the moderated t
#   enrichment/primary_GSEA_KEGG.tsv.gz         canonical, ranked by the moderated t
#   enrichment/primary_ORA_GO_BP.tsv.gz         all four query-list variants
#   enrichment/GSEA_log2FC_sensitivity.tsv.gz   SENSITIVITY ONLY, GO-BP + KEGG;
#                                               effect-size-ranked, NOT a log2 fold change
#   enrichment/primary_EWCE.tsv.gz              cell-type enrichment, all settings
#   enrichment/enrichment_run_parameters.tsv    the canonical run parameters as recorded
#
# NOTHING IS RECOMPUTED. Rows are read from the canonical enrichment folder and relabelled.
#
# The canonical GSEA ranking statistic is the moderated t. The secondary ranking exists only
# as a sensitivity analysis, is written to a SEPARATE file, and every row in every file
# carries analysis_role so canonical and sensitivity rows can never be silently pooled.
#
# TERMINOLOGY. That sensitivity ranking is publicly an EFFECT-SIZE-RANKED sensitivity
# analysis, not a "log2FC-ranked" one: the statistic it ranks by is the standardised-
# abundance effect size ProTigy stores as `logFC`, which is not a log2 fold change. The
# FILENAMES and the `rank_statistic` VALUE keep the historical `log2fc` token, because that
# is how this export is matched back to the canonical run and renaming them would break the
# provenance chain. The documentation explains what the statistic is; the identifiers stay.
#
# Note on cutoffs: the canonical run used pvalue_cutoff = 1 and qvalue_cutoff = 1, so these
# tables are the COMPLETE tested result, not a significant-only subset. Filter on
# adjusted_p_value yourself; do not assume presence in the file implies significance.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_validation.R"))
release_require("digest", "readxl")

DATA_ROOT <- release_data_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/05_build_enrichment_exports.R"

release_banner("stage 05 -- enrichment exports")

ENRICH_ROOT <- file.path(DATA_ROOT, "03_output", "enrichment",
                         "enrichment_t_rank_validation_20260825")
EWCE_ROOT <- file.path(DATA_ROOT, "03_output", "ewce",
                       "EWCE_Results_animal_level_validation_20260825")

paths <- list(
  enrichment_index = file.path(ENRICH_ROOT, "indexEnrichmentComparisons.csv"),
  run_parameters = file.path(ENRICH_ROOT, "audits", "run_parameters.csv"),
  package_versions = file.path(ENRICH_ROOT, "audits", "package_database_versions.csv"),
  ewce_table = file.path(EWCE_ROOT, "02_Tables_Supplements", "Supplementary_Table_EWCE.xlsx"),
  ewce_session = file.path(EWCE_ROOT, "03_QC_Mapping_Logs", "reproducibility_session_info.txt"),
  contrast_manifest = release_path("metadata", "primary_contrast_manifest.tsv", create_dir = FALSE)
)
for (nm in names(paths)) release_assert_exists(paths[[nm]], nm)

enrich_index <- release_read_csv(paths$enrichment_index)
run_parameters <- release_read_csv(paths$run_parameters)
contrasts <- release_read_tsv_plain(paths$contrast_manifest)

if (nrow(enrich_index) != RELEASE_DESIGN_INVARIANTS$n_primary_contrasts) {
  stop("Enrichment index covers ", nrow(enrich_index), " comparisons, expected 12.",
       call. = FALSE)
}
if (!all(as.character(enrich_index$execution_status) == "success")) {
  stop("Not every canonical enrichment comparison ran to success.", call. = FALSE)
}
# The canonical run must still say it ranked by the moderated t.
if (!all(as.character(enrich_index$ranking_statistic) == "t") ||
    !all(as.character(enrich_index$gsea_analysis_role) == "canonical") ||
    !all(as.character(enrich_index$sensitivity_ranking_statistic) == "log2fc") ||
    !all(as.character(enrich_index$sensitivity_gsea_analysis_role) == "sensitivity")) {
  stop("Canonical enrichment run no longer declares t-ranked canonical / log2fc-ranked ",
       "sensitivity. Refusing to relabel it.", call. = FALSE)
}
release_log("  canonical ranking statistic confirmed: moderated t (sensitivity ranks by ",
            "the standardised-abundance effect size, recorded upstream as log2fc)")

ORA_UNIVERSE <- unique(as.character(enrich_index$ora_universe_definition))
release_log("  ORA universe: ", paste(ORA_UNIVERSE, collapse = "; "))

# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------

contrast_context <- function(canonical_comparison) {
  row <- contrasts[contrasts$canonical_comparison == canonical_comparison, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("No unique primary contrast for comparison ", canonical_comparison, call. = FALSE)
  }
  row
}

read_result_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  d <- tryCatch(release_read_csv(path), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d
}

count_members <- function(x) {
  x <- as.character(x)
  out <- lengths(strsplit(x, "/", fixed = TRUE))
  out[is.na(x) | !nzchar(x)] <- 0L
  as.integer(out)
}

nes_direction <- function(nes, numerator_condition, denominator_condition) {
  ifelse(is.na(nes), NA_character_,
         ifelse(nes > 0, paste0("enriched_toward_", numerator_condition),
                ifelse(nes < 0, paste0("enriched_toward_", denominator_condition), "no_direction")))
}

GSEA_COLUMNS <- c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue",
                  "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")

build_gsea_block <- function(path, comparison, ontology, analysis_role, rank_statistic,
                             rank_source_column) {
  d <- read_result_csv(path)
  if (is.null(d)) return(NULL)
  missing <- setdiff(GSEA_COLUMNS, names(d))
  if (length(missing)) {
    stop("GSEA result ", path, " is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  ctx <- contrast_context(comparison)
  data.frame(
    sample_class = ctx$sample_class,
    contrast_family = ctx$contrast_family,
    canonical_comparison = ctx$canonical_comparison,
    canonical_contrast = ctx$canonical_contrast,
    numerator_condition = ctx$numerator_condition,
    denominator_condition = ctx$denominator_condition,
    n_numerator_animals = as.integer(ctx$n_numerator_animals),
    n_denominator_animals = as.integer(ctx$n_denominator_animals),
    analysis_unit = "animal",
    analysis = "GSEA",
    ontology = ontology,
    query_list = NA_character_,
    term_id = as.character(d$ID),
    term_name = as.character(d$Description),
    set_size = as.integer(d$setSize),
    enrichment_score = as.numeric(d$enrichmentScore),
    NES = as.numeric(d$NES),
    p_value = as.numeric(d$pvalue),
    adjusted_p_value = as.numeric(d[["p.adjust"]]),
    q_value = as.numeric(d$qvalue),
    gene_count = NA_integer_,
    gene_ratio = NA_character_,
    background_ratio = NA_character_,
    fold_enrichment = NA_real_,
    rank_at_max = as.integer(d$rank),
    leading_edge = as.character(d$leading_edge),
    core_enrichment = as.character(d$core_enrichment),
    core_enrichment_definition = "GSEA leading-edge subset (clusterProfiler core_enrichment)",
    n_core_enrichment_genes = count_members(d$core_enrichment),
    direction = nes_direction(as.numeric(d$NES), ctx$numerator_condition,
                              ctx$denominator_condition),
    rank_statistic = rank_statistic,
    rank_source_column = rank_source_column,
    analysis_role = analysis_role,
    gene_identifier_type = "UNIPROT",
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

ORA_COLUMNS <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust",
                 "qvalue", "geneID", "Count")

ORA_VARIANTS <- list(
  list(file = "ORA_GO_BP_fdr_all.csv", query_list = "fdr_significant_all",
       direction = "both", analysis_role = "canonical"),
  list(file = "ORA_GO_BP_fdr_up.csv", query_list = "fdr_significant_higher_in_numerator",
       direction = "up", analysis_role = "canonical"),
  list(file = "ORA_GO_BP_fdr_down.csv", query_list = "fdr_significant_higher_in_denominator",
       direction = "down", analysis_role = "canonical"),
  list(file = "ORA_GO_BP_top_abs_log2fc.csv", query_list = "top_absolute_log2fc",
       direction = "both", analysis_role = "alternative_query_list")
)

build_ora_block <- function(path, comparison, variant) {
  d <- read_result_csv(path)
  if (is.null(d)) return(NULL)
  missing <- setdiff(ORA_COLUMNS, names(d))
  if (length(missing)) {
    stop("ORA result ", path, " is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  ctx <- contrast_context(comparison)
  data.frame(
    sample_class = ctx$sample_class,
    contrast_family = ctx$contrast_family,
    canonical_comparison = ctx$canonical_comparison,
    canonical_contrast = ctx$canonical_contrast,
    numerator_condition = ctx$numerator_condition,
    denominator_condition = ctx$denominator_condition,
    n_numerator_animals = as.integer(ctx$n_numerator_animals),
    n_denominator_animals = as.integer(ctx$n_denominator_animals),
    analysis_unit = "animal",
    analysis = "ORA",
    ontology = "GO_BP",
    query_list = variant$query_list,
    term_id = as.character(d$ID),
    term_name = as.character(d$Description),
    set_size = NA_integer_,
    enrichment_score = NA_real_,
    NES = NA_real_,
    p_value = as.numeric(d$pvalue),
    adjusted_p_value = as.numeric(d[["p.adjust"]]),
    q_value = as.numeric(d$qvalue),
    gene_count = as.integer(d$Count),
    gene_ratio = as.character(d$GeneRatio),
    background_ratio = as.character(d$BgRatio),
    fold_enrichment = if ("FoldEnrichment" %in% names(d)) as.numeric(d$FoldEnrichment) else NA_real_,
    rich_factor = if ("RichFactor" %in% names(d)) as.numeric(d$RichFactor) else NA_real_,
    z_score = if ("zScore" %in% names(d)) as.numeric(d$zScore) else NA_real_,
    rank_at_max = NA_integer_,
    leading_edge = NA_character_,
    core_enrichment = as.character(d$geneID),
    core_enrichment_definition = "ORA query genes annotated to the term (clusterProfiler geneID)",
    n_core_enrichment_genes = count_members(d$geneID),
    direction = variant$direction,
    rank_statistic = "not_applicable",
    rank_source_column = "not_applicable",
    analysis_role = variant$analysis_role,
    gene_identifier_type = "UNIPROT",
    ora_universe_definition = paste(ORA_UNIVERSE, collapse = ";"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

# --------------------------------------------------------------------------------------
# assemble
# --------------------------------------------------------------------------------------

comparisons <- contrasts$canonical_comparison
gsea_go <- list(); gsea_kegg <- list(); ora <- list(); sensitivity <- list()
source_files <- character(0)

for (cmp in comparisons) {
  dir_cmp <- file.path(ENRICH_ROOT, "per_comparison", cmp)
  if (!dir.exists(dir_cmp)) stop("Missing canonical enrichment folder: ", dir_cmp, call. = FALSE)

  p_go <- file.path(dir_cmp, "GSEA_GO_BP.csv")
  p_kegg <- file.path(dir_cmp, "GSEA_KEGG.csv")
  p_go_s <- file.path(dir_cmp, "GSEA_GO_BP_log2fc_sensitivity.csv")
  p_kegg_s <- file.path(dir_cmp, "GSEA_KEGG_log2fc_sensitivity.csv")

  gsea_go[[cmp]] <- build_gsea_block(p_go, cmp, "GO_BP", "canonical", "moderated_t", "t")
  gsea_kegg[[cmp]] <- build_gsea_block(p_kegg, cmp, "KEGG", "canonical", "moderated_t", "t")
  sensitivity[[paste0(cmp, "|GO_BP")]] <-
    build_gsea_block(p_go_s, cmp, "GO_BP", "sensitivity", "log2fc", "log2fc")
  sensitivity[[paste0(cmp, "|KEGG")]] <-
    build_gsea_block(p_kegg_s, cmp, "KEGG", "sensitivity", "log2fc", "log2fc")
  source_files <- c(source_files, p_go, p_kegg, p_go_s, p_kegg_s)

  for (variant in ORA_VARIANTS) {
    p_ora <- file.path(dir_cmp, variant$file)
    ora[[paste0(cmp, "|", variant$query_list)]] <- build_ora_block(p_ora, cmp, variant)
    source_files <- c(source_files, p_ora)
  }
}

bind_blocks <- function(blocks, label) {
  blocks <- blocks[!vapply(blocks, is.null, logical(1))]
  if (!length(blocks)) stop("No rows assembled for ", label, call. = FALSE)
  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}

gsea_go_tbl <- bind_blocks(gsea_go, "GSEA GO-BP")
gsea_kegg_tbl <- bind_blocks(gsea_kegg, "GSEA KEGG")
ora_tbl <- bind_blocks(ora, "ORA GO-BP")
sens_tbl <- bind_blocks(sensitivity, "GSEA effect-size-ranked sensitivity")

for (nm in names(list(gsea_go_tbl = gsea_go_tbl, gsea_kegg_tbl = gsea_kegg_tbl))) {
  tbl <- get(nm)
  if (length(unique(tbl$canonical_comparison)) != 12L) {
    stop(nm, " does not cover all 12 comparisons.", call. = FALSE)
  }
}
if (!all(gsea_go_tbl$analysis_role == "canonical") ||
    !all(gsea_kegg_tbl$analysis_role == "canonical")) {
  stop("A non-canonical row leaked into a canonical GSEA table.", call. = FALSE)
}
if (!all(sens_tbl$analysis_role == "sensitivity") ||
    !all(sens_tbl$rank_statistic == "log2fc")) {
  stop("A canonical row leaked into the sensitivity table.", call. = FALSE)
}
release_log("  GSEA GO-BP:        ", nrow(gsea_go_tbl), " rows")
release_log("  GSEA KEGG:         ", nrow(gsea_kegg_tbl), " rows")
release_log("  ORA GO-BP:         ", nrow(ora_tbl), " rows over ",
            length(unique(ora_tbl$query_list)), " query-list variants")
release_log("  GSEA sensitivity:  ", nrow(sens_tbl),
            " rows (effect-size-ranked, source field logFC; NOT canonical)")

# --------------------------------------------------------------------------------------
# exact re-verification against the canonical enrichment files
# --------------------------------------------------------------------------------------

release_log("  re-verifying exported enrichment values against canonical sources ...")
verify_pairs <- list(
  list(tbl = gsea_go_tbl, file = "GSEA_GO_BP.csv", role = "canonical", ontology = "GO_BP"),
  list(tbl = gsea_kegg_tbl, file = "GSEA_KEGG.csv", role = "canonical", ontology = "KEGG"),
  list(tbl = sens_tbl, file = "GSEA_GO_BP_log2fc_sensitivity.csv", role = "sensitivity",
       ontology = "GO_BP"),
  list(tbl = sens_tbl, file = "GSEA_KEGG_log2fc_sensitivity.csv", role = "sensitivity",
       ontology = "KEGG")
)
GSEA_CHECK <- c(enrichment_score = "enrichmentScore", NES = "NES", p_value = "pvalue",
                adjusted_p_value = "p.adjust", q_value = "qvalue")
n_verified <- 0L
for (vp in verify_pairs) {
  for (cmp in comparisons) {
    d <- read_result_csv(file.path(ENRICH_ROOT, "per_comparison", cmp, vp$file))
    if (is.null(d)) next
    ex <- vp$tbl[vp$tbl$canonical_comparison == cmp & vp$tbl$ontology == vp$ontology &
                   vp$tbl$analysis_role == vp$role, , drop = FALSE]
    if (!identical(ex$term_id, as.character(d$ID))) {
      stop("Term order diverged for ", cmp, " / ", vp$file, call. = FALSE)
    }
    for (pub in names(GSEA_CHECK)) {
      if (!identical(ex[[pub]], as.numeric(d[[GSEA_CHECK[[pub]]]]))) {
        stop("Exported ", pub, " differs from canonical ", GSEA_CHECK[[pub]], " in ", cmp,
             " / ", vp$file, call. = FALSE)
      }
    }
    n_verified <- n_verified + 1L
  }
}
ORA_CHECK <- c(p_value = "pvalue", adjusted_p_value = "p.adjust", q_value = "qvalue",
               gene_count = "Count")
for (variant in ORA_VARIANTS) {
  for (cmp in comparisons) {
    d <- read_result_csv(file.path(ENRICH_ROOT, "per_comparison", cmp, variant$file))
    if (is.null(d)) next
    ex <- ora_tbl[ora_tbl$canonical_comparison == cmp &
                    ora_tbl$query_list == variant$query_list, , drop = FALSE]
    if (!identical(ex$term_id, as.character(d$ID))) {
      stop("Term order diverged for ", cmp, " / ", variant$file, call. = FALSE)
    }
    for (pub in names(ORA_CHECK)) {
      expected <- if (pub == "gene_count") as.integer(d[[ORA_CHECK[[pub]]]]) else
        as.numeric(d[[ORA_CHECK[[pub]]]])
      if (!identical(ex[[pub]], expected)) {
        stop("Exported ", pub, " differs from canonical ", ORA_CHECK[[pub]], " in ", cmp,
             " / ", variant$file, call. = FALSE)
      }
    }
    n_verified <- n_verified + 1L
  }
}
release_log("  verified ", n_verified, " canonical enrichment result files identical")

# --------------------------------------------------------------------------------------
# EWCE
# --------------------------------------------------------------------------------------

ewce_raw <- as.data.frame(
  readxl::read_excel(paths$ewce_table, sheet = "Full_Results", guess_max = 100000,
                     .name_repair = "minimal"),
  stringsAsFactors = FALSE, check.names = FALSE)
release_log("  EWCE Full_Results: ", nrow(ewce_raw), " rows x ", ncol(ewce_raw), " cols")

ewce_session <- readLines(paths$ewce_session, warn = FALSE)
grab <- function(prefix) {
  hit <- grep(paste0("^", prefix), ewce_session, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(paste0("^", prefix, "\\s*:?\\s*"), "", hit[[1]]))
}
ewce_primary_topn <- suppressWarnings(as.integer(grab("Primary top-N")))
ewce_primary_annot <- suppressWarnings(as.integer(grab("Primary annotation level")))
ewce_reps <- grab("EWCE reps")
ewce_background <- grab("Measured proteome background size")
if (is.na(ewce_primary_topn) || is.na(ewce_primary_annot)) {
  stop("Could not read the primary EWCE settings from the canonical session record.",
       call. = FALSE)
}
release_log("  EWCE primary settings: top-N ", ewce_primary_topn, ", annotation level ",
            ewce_primary_annot, ", ", ewce_reps, " bootstrap reps, background ",
            ewce_background)

diff_rows <- as.character(ewce_raw$AnalysisType) == "Differential"
ctx_idx <- match(as.character(ewce_raw$Contrast), contrasts$canonical_comparison)
if (any(diff_rows & is.na(ctx_idx))) {
  stop("An EWCE differential row names a comparison that is not one of the 12 primary ",
       "contrasts.", call. = FALSE)
}

ewce <- data.frame(
  sample_class = as.character(ewce_raw$Stratum),
  analysis = "EWCE",
  ewce_analysis_type = as.character(ewce_raw$AnalysisType),
  contrast_family = ifelse(diff_rows, contrasts$contrast_family[ctx_idx], NA_character_),
  canonical_comparison = ifelse(diff_rows, as.character(ewce_raw$Contrast), NA_character_),
  canonical_contrast = ifelse(diff_rows, contrasts$canonical_contrast[ctx_idx], NA_character_),
  numerator_condition = ifelse(diff_rows, contrasts$numerator_condition[ctx_idx],
                               as.character(ewce_raw$Metric)),
  denominator_condition = ifelse(diff_rows, contrasts$denominator_condition[ctx_idx],
                                 NA_character_),
  n_numerator_animals = ifelse(diff_rows,
                               as.integer(contrasts$n_numerator_animals[ctx_idx]), 3L),
  n_denominator_animals = ifelse(diff_rows,
                                 as.integer(contrasts$n_denominator_animals[ctx_idx]),
                                 NA_integer_),
  analysis_unit = "animal",
  protein_list_definition = as.character(ewce_raw$Metric),
  direction = as.character(ewce_raw$Direction),
  cell_type = as.character(ewce_raw$CellType),
  annotation_level = as.integer(ewce_raw$AnnotLevel),
  top_n = as.integer(ewce_raw$TopN),
  n_hits = as.integer(ewce_raw$N_Hits),
  n_background = as.integer(ewce_raw$N_Background),
  fold_change = as.numeric(ewce_raw$fold_change),
  sd_from_mean = as.numeric(ewce_raw$sd_from_mean),
  p_value = as.numeric(ewce_raw$p),
  adjusted_p_value_within_target = as.numeric(ewce_raw$q_target),
  adjusted_p_value_global = as.numeric(ewce_raw$q_global),
  significant_global_fdr_0_05 = as.logical(ewce_raw$Significant_Global),
  ewce_direction = as.character(ewce_raw$Direction_EWCE),
  bootstrap_replicates = ewce_reps,
  is_primary_setting = as.integer(ewce_raw$TopN) == ewce_primary_topn &
    as.integer(ewce_raw$AnnotLevel) == ewce_primary_annot,
  analysis_role = ifelse(as.integer(ewce_raw$TopN) == ewce_primary_topn &
                           as.integer(ewce_raw$AnnotLevel) == ewce_primary_annot,
                         "canonical", "sensitivity"),
  primary_or_secondary = ifelse(as.character(ewce_raw$AnalysisType) == "Differential",
                                "primary", "secondary"),
  stringsAsFactors = FALSE, check.names = FALSE
)
if (!identical(ewce$p_value, as.numeric(ewce_raw$p)) ||
    !identical(ewce$adjusted_p_value_global, as.numeric(ewce_raw$q_global)) ||
    !identical(ewce$sd_from_mean, as.numeric(ewce_raw$sd_from_mean))) {
  stop("Exported EWCE statistics differ from the canonical table.", call. = FALSE)
}
release_log("  EWCE export: ", nrow(ewce), " rows (", sum(ewce$is_primary_setting),
            " at the primary setting; ", sum(ewce$primary_or_secondary == "primary"),
            " differential)")

# --------------------------------------------------------------------------------------
# write
# --------------------------------------------------------------------------------------

source_files <- unique(source_files[file.exists(source_files)])
source_hashes <- vapply(source_files, release_sha256, character(1), USE.NAMES = FALSE)
index_hash <- release_sha256(paths$enrichment_index)

emit <- function(tbl, relative, role_text, sources, hashes) {
  w <- release_write_table(tbl, release_path(dirname(relative), basename(relative)))
  release_register(relative, role_text, sources, hashes, STAGE, "tsv.gz")
  release_log("  wrote ", basename(relative), " (", w$rows, "x", w$cols, ")")
  invisible(w)
}

canon_go <- grep("GSEA_GO_BP[.]csv$", source_files, value = TRUE)
canon_kegg <- grep("GSEA_KEGG[.]csv$", source_files, value = TRUE)
canon_ora <- grep("ORA_GO_BP", source_files, value = TRUE)
canon_sens <- grep("log2fc_sensitivity[.]csv$", source_files, value = TRUE)
hash_of <- function(x) source_hashes[match(x, source_files)]

emit(gsea_go_tbl, "enrichment/primary_GSEA_GO_BP.tsv.gz",
     "canonical GSEA GO-BP, ranked by the moderated t statistic; complete tested result",
     c(canon_go, paths$enrichment_index), c(hash_of(canon_go), index_hash))
emit(gsea_kegg_tbl, "enrichment/primary_GSEA_KEGG.tsv.gz",
     "canonical GSEA KEGG, ranked by the moderated t statistic; complete tested result",
     c(canon_kegg, paths$enrichment_index), c(hash_of(canon_kegg), index_hash))
emit(ora_tbl, "enrichment/primary_ORA_GO_BP.tsv.gz",
     "GO-BP over-representation analysis over four query-list definitions",
     c(canon_ora, paths$enrichment_index), c(hash_of(canon_ora), index_hash))
emit(sens_tbl, "enrichment/GSEA_log2FC_sensitivity.tsv.gz",
     paste("SENSITIVITY ONLY:", RELEASE_EFFECT_SIZE$sensitivity_public_label,
           "-- GSEA GO-BP and KEGG ranked by the standardised-abundance effect size",
           "(stored as logFC) instead of the moderated t. Not a log2 fold change."),
     c(canon_sens, paths$enrichment_index), c(hash_of(canon_sens), index_hash))
emit(ewce, "enrichment/primary_EWCE.tsv.gz",
     "EWCE cell-type enrichment across all top-N and annotation-level settings",
     paths$ewce_table, release_sha256(paths$ewce_table))

# Coverage. Several ORA query lists are legitimately EMPTY: with no FDR-significant
# proteins in a comparison there is no query list to over-represent, so the canonical run
# produced a zero-row table. Recording that explicitly stops an absent block being read as
# a missing file.
coverage <- do.call(rbind, lapply(comparisons, function(cmp) {
  ix <- enrich_index[as.character(enrich_index$canonical_comparison) == cmp, , drop = FALSE]
  count_from_index <- function(col) {
    if (!col %in% names(ix)) return(NA_integer_)
    suppressWarnings(as.integer(ix[[col]][[1]]))
  }
  rows_in <- function(tbl, ...) {
    sel <- rep(TRUE, nrow(tbl)); filters <- list(...)
    for (nm in names(filters)) sel <- sel & tbl[[nm]] == filters[[nm]]
    sum(sel, na.rm = TRUE)
  }
  rbind(
    data.frame(canonical_comparison = cmp, analysis = "GSEA", ontology = "GO_BP",
               query_list = NA_character_, analysis_role = "canonical",
               rows_exported = rows_in(gsea_go_tbl, canonical_comparison = cmp),
               term_count_recorded_by_run = count_from_index("go_gsea_term_count"),
               stringsAsFactors = FALSE, check.names = FALSE),
    data.frame(canonical_comparison = cmp, analysis = "GSEA", ontology = "KEGG",
               query_list = NA_character_, analysis_role = "canonical",
               rows_exported = rows_in(gsea_kegg_tbl, canonical_comparison = cmp),
               term_count_recorded_by_run = count_from_index("kegg_gsea_term_count"),
               stringsAsFactors = FALSE, check.names = FALSE),
    data.frame(canonical_comparison = cmp, analysis = "GSEA", ontology = "GO_BP",
               query_list = NA_character_, analysis_role = "sensitivity",
               rows_exported = rows_in(sens_tbl, canonical_comparison = cmp,
                                       ontology = "GO_BP"),
               term_count_recorded_by_run = count_from_index("sensitivity_go_gsea_term_count"),
               stringsAsFactors = FALSE, check.names = FALSE),
    data.frame(canonical_comparison = cmp, analysis = "GSEA", ontology = "KEGG",
               query_list = NA_character_, analysis_role = "sensitivity",
               rows_exported = rows_in(sens_tbl, canonical_comparison = cmp,
                                       ontology = "KEGG"),
               term_count_recorded_by_run = count_from_index("sensitivity_kegg_gsea_term_count"),
               stringsAsFactors = FALSE, check.names = FALSE),
    do.call(rbind, lapply(ORA_VARIANTS, function(variant) {
      idx_col <- switch(variant$query_list,
                        fdr_significant_all = "go_ora_fdr_all_term_count",
                        fdr_significant_higher_in_numerator = "go_ora_fdr_up_term_count",
                        fdr_significant_higher_in_denominator = "go_ora_fdr_down_term_count",
                        top_absolute_log2fc = "go_ora_top_abs_log2fc_term_count")
      data.frame(canonical_comparison = cmp, analysis = "ORA", ontology = "GO_BP",
                 query_list = variant$query_list, analysis_role = variant$analysis_role,
                 rows_exported = rows_in(ora_tbl, canonical_comparison = cmp,
                                         query_list = variant$query_list),
                 term_count_recorded_by_run = count_from_index(idx_col),
                 stringsAsFactors = FALSE, check.names = FALSE)
    }))
  )
}))
coverage$agrees_with_run_record <- coverage$rows_exported == coverage$term_count_recorded_by_run
coverage$note <- ifelse(
  coverage$rows_exported == 0L,
  "empty by construction: the canonical run produced no terms for this query list",
  NA_character_)
if (any(!coverage$agrees_with_run_record, na.rm = TRUE)) {
  bad <- coverage[!coverage$agrees_with_run_record, , drop = FALSE]
  stop("Exported enrichment row counts disagree with the counts the canonical run ",
       "recorded:\n",
       paste0("  - ", bad$canonical_comparison, " ", bad$analysis, " ", bad$ontology, " ",
              ifelse(is.na(bad$query_list), "", bad$query_list), ": exported ",
              bad$rows_exported, ", recorded ", bad$term_count_recorded_by_run,
              collapse = "\n"), call. = FALSE)
}
release_log("  coverage table: ", nrow(coverage), " analysis blocks, ",
            sum(coverage$rows_exported == 0L), " legitimately empty; all row counts agree ",
            "with the canonical run record")
w <- release_write_table(coverage, release_path("enrichment", "enrichment_coverage.tsv"))
release_register("enrichment/enrichment_coverage.tsv",
                 paste("per-comparison enrichment row counts, cross-checked against the",
                       "term counts the canonical run recorded"),
                 paths$enrichment_index, index_hash, STAGE, "tsv")
release_log("  wrote enrichment_coverage.tsv (", w$rows, "x", w$cols, ")")

# The canonical run parameters, republished verbatim plus the EWCE settings.
params <- rbind(
  data.frame(analysis = "enrichment (GSEA/ORA)", parameter = as.character(run_parameters$parameter),
             value = as.character(run_parameters$value),
             stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(analysis = "enrichment (GSEA/ORA)", parameter = "ora_universe_definition_per_comparison",
             value = paste(ORA_UNIVERSE, collapse = ";"),
             stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(analysis = "EWCE",
             parameter = c("bootstrap_replicates", "top_n_values", "annotation_levels",
                           "primary_top_n", "primary_annotation_level",
                           "measured_proteome_background_size", "aggregation_policy",
                           "input_sha256"),
             value = c(ewce_reps, grab("Top-N values"), grab("Annotation levels"),
                       as.character(ewce_primary_topn), as.character(ewce_primary_annot),
                       ewce_background, grab("Aggregation policy"), grab("Input SHA256")),
             stringsAsFactors = FALSE, check.names = FALSE)
)
w <- release_write_table(params, release_path("enrichment", "enrichment_run_parameters.tsv"))
release_register("enrichment/enrichment_run_parameters.tsv",
                 "canonical enrichment and EWCE run parameters, as recorded by the runs",
                 c(paths$run_parameters, paths$ewce_session),
                 c(release_sha256(paths$run_parameters), release_sha256(paths$ewce_session)),
                 STAGE, "tsv")
release_log("  wrote enrichment_run_parameters.tsv (", w$rows, "x", w$cols, ")")
release_log("stage 05 complete")
