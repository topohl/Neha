#!/usr/bin/env Rscript

# Publication release, stage 05 -- normalised enrichment and cell-type exports.
#
# Produces
#   enrichment/primary_GSEA_GO_BP.tsv.gz        canonical, ranked by the moderated t
#   enrichment/primary_GSEA_KEGG.tsv.gz         canonical, ranked by the moderated t
#   enrichment/primary_ORA_GO_BP.tsv.gz         all four query-list variants
#   enrichment/primary_GSEA_GO_CC.tsv.gz        canonical, ranked by the moderated t
#   enrichment/primary_ORA_GO_CC.tsv.gz         all four query-list variants
#   enrichment/GSEA_log2FC_sensitivity.tsv.gz   SENSITIVITY ONLY, GO-BP + GO-CC + KEGG;
#                                               effect-size-ranked, NOT a log2 fold change
#   enrichment/primary_EWCE.tsv.gz              cell-type enrichment, all settings
#   enrichment/gsea_untested_gene_sets.tsv      GO gene sets fgsea gave no P value
#   enrichment/enrichment_run_parameters.tsv    the canonical run parameters as recorded
#
# Nothing is recomputed. Rows are read from the enrichment runs on the per-compartment
# statistics and relabelled: GO-BP and GO-CC from enrichment_within_compartment_20261125
# (fgsea nPermSimple = 100000, adopted 2026-11-25), KEGG from the GO_BP folder of
# enrichment_within_compartment_20261124 (not rerun). Every run must name the locked
# statistics GCT as its source. GO-CC has been a primary output again since 2026-11-24.
#
# Canonical GSEA ranks by the moderated t. GSEA ranked by the effect size (stored as
# `logFC`, not a log2 fold change) is a sensitivity analysis in its own file, and every row
# carries analysis_role. Its file names and `rank_statistic` value keep the old `log2fc`
# token so they still match the canonical run.
#
# The runs used pvalue_cutoff = 1 and qvalue_cutoff = 1, so the tables hold every tested
# term; filter on adjusted_p_value.

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

ENRICH_ROOT <- release_enrichment_root("GO_BP", DATA_ROOT)     # GO-BP
ENRICH_ROOT_CC <- release_enrichment_root("GO_CC", DATA_ROOT)  # GO-CC
KEGG_ROOT <- release_kegg_root(DATA_ROOT)                      # KEGG (20261124 GO_BP folder)
EWCE_ROOT <- file.path(DATA_ROOT, "03_output", "ewce",
                       "EWCE_Results_animal_level_validation_20260825")

paths <- list(
  enrichment_index = file.path(ENRICH_ROOT, "indexEnrichmentComparisons.csv"),
  enrichment_index_cc = file.path(ENRICH_ROOT_CC, "indexEnrichmentComparisons.csv"),
  enrichment_index_kegg = file.path(KEGG_ROOT, "indexEnrichmentComparisons.csv"),
  run_parameters = file.path(ENRICH_ROOT, "audits", "run_parameters.csv"),
  run_parameters_cc = file.path(ENRICH_ROOT_CC, "audits", "run_parameters.csv"),
  run_parameters_kegg = file.path(KEGG_ROOT, "audits", "run_parameters.csv"),
  package_versions = file.path(ENRICH_ROOT, "audits", "package_database_versions.csv"),
  package_versions_cc = file.path(ENRICH_ROOT_CC, "audits", "package_database_versions.csv"),
  package_versions_kegg = file.path(KEGG_ROOT, "audits", "package_database_versions.csv"),
  ewce_table = file.path(EWCE_ROOT, "02_Tables_Supplements", "Supplementary_Table_EWCE.xlsx"),
  ewce_session = file.path(EWCE_ROOT, "03_QC_Mapping_Logs", "reproducibility_session_info.txt"),
  contrast_manifest = release_path("metadata", "primary_contrast_manifest.tsv", create_dir = FALSE)
)
for (nm in names(paths)) release_assert_exists(paths[[nm]], nm)

# All three enrichment indexes are locked: a rerun or forced rewrite of any run must stop the
# release here, not be relabelled into it.
for (key in c("enrichment_index_go_bp", "enrichment_index_go_cc", "enrichment_index_kegg")) {
  observed <- release_sha256(release_locked_artefact_path(key, DATA_ROOT))
  if (!identical(observed, RELEASE_LOCKED_ARTEFACTS[[key]]$sha256)) {
    stop("Locked enrichment index changed (", key, "): ",
         release_locked_artefact_path(key, DATA_ROOT), call. = FALSE)
  }
}
enrich_index <- release_read_csv(paths$enrichment_index)
enrich_index_cc <- release_read_csv(paths$enrichment_index_cc)
enrich_index_kegg <- release_read_csv(paths$enrichment_index_kegg)
permutations <- release_enrichment_permutations(DATA_ROOT)
run_parameters <- release_read_csv(paths$run_parameters)
run_parameters_cc <- release_read_csv(paths$run_parameters_cc)
run_parameters_kegg <- release_read_csv(paths$run_parameters_kegg)
contrasts <- release_read_tsv_plain(paths$contrast_manifest)

check_enrichment_index <- function(ix, label, ontology_token) {
  if (nrow(ix) != RELEASE_DESIGN_INVARIANTS$n_primary_contrasts ||
      !setequal(as.character(ix$canonical_comparison), contrasts$canonical_comparison)) {
    stop(label, " enrichment index does not hold the 12 primary comparisons.", call. = FALSE)
  }
  if (!all(as.character(ix$execution_status) == "success")) {
    stop("Not every ", label, " enrichment comparison ran to success.", call. = FALSE)
  }
  if (!all(as.character(ix$ontology) == ontology_token)) {
    stop(label, " enrichment index does not declare ontology ", ontology_token, ".",
         call. = FALSE)
  }
  # The canonical run must still say it ranked by the moderated t.
  if (!all(as.character(ix$ranking_statistic) == "t") ||
      !all(as.character(ix$gsea_analysis_role) == "canonical") ||
      !all(as.character(ix$sensitivity_ranking_statistic) == "log2fc") ||
      !all(as.character(ix$sensitivity_gsea_analysis_role) == "sensitivity")) {
    stop(label, " enrichment run no longer declares t-ranked canonical / log2fc-ranked ",
         "sensitivity. Refusing to relabel it.", call. = FALSE)
  }
  # Enrichment and statistics must come from the same run: the index records the GCT its
  # ranked lists were split from.
  if (!all(as.character(ix$source_gct_sha256) == RELEASE_LOCKED_ARTEFACTS$stat_gct$sha256)) {
    stop(label, " enrichment was not run on the locked statistics GCT.", call. = FALSE)
  }
  invisible(TRUE)
}
check_enrichment_index(enrich_index, "GO-BP", "BP")
check_enrichment_index(enrich_index_cc, "GO-CC", "CC")
check_enrichment_index(enrich_index_kegg, "KEGG", "BP")
# The GO and KEGG runs must share every package they both record (the GO rerun also records
# BiocParallel, which it uses for its serial backend).
pkg_go <- release_read_csv(paths$package_versions)
if (!identical(pkg_go, release_read_csv(paths$package_versions_cc))) {
  stop("The GO-BP and GO-CC enrichment runs recorded different package versions.",
       call. = FALSE)
}
pkg_kegg <- release_read_csv(paths$package_versions_kegg)
if (!all(pkg_kegg$component %in% pkg_go$component) ||
    !identical(as.character(pkg_kegg$version),
               as.character(pkg_go$version[match(pkg_kegg$component, pkg_go$component)]))) {
  stop("The KEGG run and the GO runs recorded different versions of a shared package.",
       call. = FALSE)
}
release_log("  canonical ranking statistic confirmed: moderated t (sensitivity ranks by ",
            "the standardised-abundance effect size, recorded upstream as log2fc)")
release_log("  all three enrichment runs (GO-BP, GO-CC, KEGG) record the locked statistics GCT ",
            substr(RELEASE_LOCKED_ARTEFACTS$stat_gct$sha256, 1, 16))
release_log("  ", permutations$setting)

ORA_UNIVERSE <- unique(as.character(c(enrich_index$ora_universe_definition,
                                      enrich_index_cc$ora_universe_definition)))
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

# The GO-CC index keeps the generic go_* count columns, so one variant list serves both.
ora_variants <- function(ontology) list(
  list(file = sprintf("ORA_%s_fdr_all.csv", ontology), query_list = "fdr_significant_all",
       direction = "both", analysis_role = "canonical",
       index_column = "go_ora_fdr_all_term_count"),
  list(file = sprintf("ORA_%s_fdr_up.csv", ontology),
       query_list = "fdr_significant_higher_in_numerator",
       direction = "up", analysis_role = "canonical",
       index_column = "go_ora_fdr_up_term_count"),
  list(file = sprintf("ORA_%s_fdr_down.csv", ontology),
       query_list = "fdr_significant_higher_in_denominator",
       direction = "down", analysis_role = "canonical",
       index_column = "go_ora_fdr_down_term_count"),
  list(file = sprintf("ORA_%s_top_abs_log2fc.csv", ontology), query_list = "top_absolute_log2fc",
       direction = "both", analysis_role = "alternative_query_list",
       index_column = "go_ora_top_abs_log2fc_term_count")
)
ORA_VARIANTS <- ora_variants("GO_BP")
ORA_VARIANTS_CC <- ora_variants("GO_CC")

build_ora_block <- function(path, comparison, variant, ontology = "GO_BP") {
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
    ontology = ontology,
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
gsea_cc <- list(); ora_cc <- list()
source_files <- character(0)

for (cmp in comparisons) {
  dir_cmp <- file.path(ENRICH_ROOT, "per_comparison", cmp)
  if (!dir.exists(dir_cmp)) stop("Missing canonical enrichment folder: ", dir_cmp, call. = FALSE)

  dir_kegg <- file.path(KEGG_ROOT, "per_comparison", cmp)
  if (!dir.exists(dir_kegg)) stop("Missing canonical KEGG enrichment folder: ", dir_kegg, call. = FALSE)

  p_go <- file.path(dir_cmp, "GSEA_GO_BP.csv")
  p_kegg <- file.path(dir_kegg, "GSEA_KEGG.csv")
  p_go_s <- file.path(dir_cmp, "GSEA_GO_BP_log2fc_sensitivity.csv")
  p_kegg_s <- file.path(dir_kegg, "GSEA_KEGG_log2fc_sensitivity.csv")

  gsea_go[[cmp]] <- build_gsea_block(p_go, cmp, "GO_BP", "canonical", "moderated_t", "t")
  gsea_kegg[[cmp]] <- build_gsea_block(p_kegg, cmp, "KEGG", "canonical", "moderated_t", "t")
  sensitivity[[paste0(cmp, "|GO_BP")]] <-
    build_gsea_block(p_go_s, cmp, "GO_BP", "sensitivity", "log2fc", "log2fc")
  sensitivity[[paste0(cmp, "|KEGG")]] <-
    build_gsea_block(p_kegg_s, cmp, "KEGG", "sensitivity", "log2fc", "log2fc")
  source_files <- c(source_files, p_go, p_kegg, p_go_s, p_kegg_s)

  for (variant in ORA_VARIANTS) {
    p_ora <- file.path(dir_cmp, variant$file)
    ora[[paste0(cmp, "|", variant$query_list)]] <- build_ora_block(p_ora, cmp, variant, "GO_BP")
    source_files <- c(source_files, p_ora)
  }

  dir_cc <- file.path(ENRICH_ROOT_CC, "per_comparison", cmp)
  if (!dir.exists(dir_cc)) stop("Missing canonical GO-CC enrichment folder: ", dir_cc, call. = FALSE)
  p_cc <- file.path(dir_cc, "GSEA_GO_CC.csv")
  p_cc_s <- file.path(dir_cc, "GSEA_GO_CC_log2fc_sensitivity.csv")
  gsea_cc[[cmp]] <- build_gsea_block(p_cc, cmp, "GO_CC", "canonical", "moderated_t", "t")
  sensitivity[[paste0(cmp, "|GO_CC")]] <-
    build_gsea_block(p_cc_s, cmp, "GO_CC", "sensitivity", "log2fc", "log2fc")
  source_files <- c(source_files, p_cc, p_cc_s)
  for (variant in ORA_VARIANTS_CC) {
    p_ora <- file.path(dir_cc, variant$file)
    ora_cc[[paste0(cmp, "|", variant$query_list)]] <- build_ora_block(p_ora, cmp, variant, "GO_CC")
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
gsea_cc_tbl <- bind_blocks(gsea_cc, "GSEA GO-CC")
ora_tbl <- bind_blocks(ora, "ORA GO-BP")
ora_cc_tbl <- bind_blocks(ora_cc, "ORA GO-CC")
sens_tbl <- bind_blocks(sensitivity, "GSEA effect-size-ranked sensitivity")

for (nm in c("gsea_go_tbl", "gsea_kegg_tbl", "gsea_cc_tbl")) {
  tbl <- get(nm)
  if (length(unique(tbl$canonical_comparison)) != 12L) {
    stop(nm, " does not cover all 12 comparisons.", call. = FALSE)
  }
}
if (!all(gsea_go_tbl$analysis_role == "canonical") ||
    !all(gsea_kegg_tbl$analysis_role == "canonical") ||
    !all(gsea_cc_tbl$analysis_role == "canonical")) {
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
release_log("  GSEA GO-CC:        ", nrow(gsea_cc_tbl), " rows")
release_log("  ORA GO-CC:         ", nrow(ora_cc_tbl), " rows over ",
            length(unique(ora_cc_tbl$query_list)), " query-list variants")
release_log("  GSEA sensitivity:  ", nrow(sens_tbl),
            " rows (effect-size-ranked, source field logFC; NOT canonical)")

# --------------------------------------------------------------------------------------
# exact re-verification against the canonical enrichment files
# --------------------------------------------------------------------------------------

release_log("  re-verifying exported enrichment values against canonical sources ...")
verify_pairs <- list(
  list(tbl = gsea_go_tbl, root = ENRICH_ROOT, file = "GSEA_GO_BP.csv", role = "canonical",
       ontology = "GO_BP"),
  list(tbl = gsea_kegg_tbl, root = KEGG_ROOT, file = "GSEA_KEGG.csv", role = "canonical",
       ontology = "KEGG"),
  list(tbl = gsea_cc_tbl, root = ENRICH_ROOT_CC, file = "GSEA_GO_CC.csv", role = "canonical",
       ontology = "GO_CC"),
  list(tbl = sens_tbl, root = ENRICH_ROOT, file = "GSEA_GO_BP_log2fc_sensitivity.csv",
       role = "sensitivity", ontology = "GO_BP"),
  list(tbl = sens_tbl, root = KEGG_ROOT, file = "GSEA_KEGG_log2fc_sensitivity.csv",
       role = "sensitivity", ontology = "KEGG"),
  list(tbl = sens_tbl, root = ENRICH_ROOT_CC, file = "GSEA_GO_CC_log2fc_sensitivity.csv",
       role = "sensitivity", ontology = "GO_CC")
)
GSEA_CHECK <- c(enrichment_score = "enrichmentScore", NES = "NES", p_value = "pvalue",
                adjusted_p_value = "p.adjust", q_value = "qvalue")
n_verified <- 0L
for (vp in verify_pairs) {
  for (cmp in comparisons) {
    d <- read_result_csv(file.path(vp$root, "per_comparison", cmp, vp$file))
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
for (os in list(list(tbl = ora_tbl, root = ENRICH_ROOT, variants = ORA_VARIANTS),
                list(tbl = ora_cc_tbl, root = ENRICH_ROOT_CC, variants = ORA_VARIANTS_CC))) {
  for (variant in os$variants) {
    for (cmp in comparisons) {
      d <- read_result_csv(file.path(os$root, "per_comparison", cmp, variant$file))
      if (is.null(d)) next
      ex <- os$tbl[os$tbl$canonical_comparison == cmp &
                     os$tbl$query_list == variant$query_list, , drop = FALSE]
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
index_cc_hash <- release_sha256(paths$enrichment_index_cc)
index_kegg_hash <- release_sha256(paths$enrichment_index_kegg)

emit <- function(tbl, relative, role_text, sources, hashes) {
  w <- release_write_table(tbl, release_path(dirname(relative), basename(relative)))
  release_register(relative, role_text, sources, hashes, STAGE, "tsv.gz")
  release_log("  wrote ", basename(relative), " (", w$rows, "x", w$cols, ")")
  invisible(w)
}

canon_go <- grep("GSEA_GO_BP[.]csv$", source_files, value = TRUE)
canon_kegg <- grep("GSEA_KEGG[.]csv$", source_files, value = TRUE)
canon_ora <- grep("ORA_GO_BP", source_files, value = TRUE)
canon_cc <- grep("GSEA_GO_CC[.]csv$", source_files, value = TRUE)
canon_ora_cc <- grep("ORA_GO_CC", source_files, value = TRUE)
canon_sens <- grep("log2fc_sensitivity[.]csv$", source_files, value = TRUE)
hash_of <- function(x) source_hashes[match(x, source_files)]

emit(gsea_go_tbl, "enrichment/primary_GSEA_GO_BP.tsv.gz",
     "canonical GSEA GO-BP, ranked by the moderated t statistic; complete tested result",
     c(canon_go, paths$enrichment_index), c(hash_of(canon_go), index_hash))
emit(gsea_kegg_tbl, "enrichment/primary_GSEA_KEGG.tsv.gz",
     "canonical GSEA KEGG, ranked by the moderated t statistic; complete tested result",
     c(canon_kegg, paths$enrichment_index_kegg), c(hash_of(canon_kegg), index_kegg_hash))
emit(ora_tbl, "enrichment/primary_ORA_GO_BP.tsv.gz",
     "GO-BP over-representation analysis over four query-list definitions",
     c(canon_ora, paths$enrichment_index), c(hash_of(canon_ora), index_hash))
emit(gsea_cc_tbl, "enrichment/primary_GSEA_GO_CC.tsv.gz",
     paste("canonical GSEA GO-CC, ranked by the moderated t statistic; complete tested",
           "result.", RELEASE_GO_CC_PROVENANCE_NOTE),
     c(canon_cc, paths$enrichment_index_cc), c(hash_of(canon_cc), index_cc_hash))
emit(ora_cc_tbl, "enrichment/primary_ORA_GO_CC.tsv.gz",
     paste("GO-CC over-representation analysis over four query-list definitions.",
           RELEASE_GO_CC_PROVENANCE_NOTE),
     c(canon_ora_cc, paths$enrichment_index_cc), c(hash_of(canon_ora_cc), index_cc_hash))
emit(sens_tbl, "enrichment/GSEA_log2FC_sensitivity.tsv.gz",
     paste("SENSITIVITY ONLY:", RELEASE_EFFECT_SIZE$sensitivity_public_label,
           "-- GSEA GO-BP, GO-CC and KEGG ranked by the standardised-abundance effect size",
           "(stored as logFC) instead of the moderated t. Not a log2 fold change."),
     c(canon_sens, paths$enrichment_index, paths$enrichment_index_cc,
       paths$enrichment_index_kegg),
     c(hash_of(canon_sens), index_hash, index_cc_hash, index_kegg_hash))
emit(ewce, "enrichment/primary_EWCE.tsv.gz",
     "EWCE cell-type enrichment across all top-N and annotation-level settings",
     paths$ewce_table, release_sha256(paths$ewce_table))

# Gene sets fgsea returned without a P value. GO: one row per set, from the GO runs'
# untested-set lists (checked against warnings.csv and the index). KEGG: counts only, from
# the KEGG run's warnings.csv; that run predates the per-set list.
untested_sets <- release_fgsea_untested_sets(comparisons, DATA_ROOT)
untested_ctx <- contrasts[match(untested_sets$canonical_comparison,
                                contrasts$canonical_comparison), , drop = FALSE]
untested_tbl <- data.frame(
  sample_class = untested_ctx$sample_class,
  contrast_family = untested_ctx$contrast_family,
  canonical_comparison = untested_sets$canonical_comparison,
  numerator_condition = untested_ctx$numerator_condition,
  denominator_condition = untested_ctx$denominator_condition,
  analysis_unit = rep("animal", nrow(untested_sets)),
  analysis = rep("GSEA", nrow(untested_sets)),
  ontology = sub("_GSEA.*$", "", untested_sets$analysis),
  term_id = as.character(untested_sets$ID),
  term_name = as.character(untested_sets$Description),
  set_size = as.integer(untested_sets$setSize),
  rank_statistic = ifelse(grepl("_log2fc_sensitivity$", untested_sets$analysis), "log2fc",
                          "moderated_t"),
  analysis_role = ifelse(grepl("_log2fc_sensitivity$", untested_sets$analysis), "sensitivity",
                         "canonical"),
  n_perm_simple = rep(permutations$go_n_perm_simple, nrow(untested_sets)),
  stringsAsFactors = FALSE, check.names = FALSE)
untested_tbl <- untested_tbl[order(match(untested_tbl$canonical_comparison, comparisons),
                                   untested_tbl$ontology, untested_tbl$analysis_role,
                                   -untested_tbl$set_size, untested_tbl$term_id), , drop = FALSE]
rownames(untested_tbl) <- NULL
# A set fgsea could not test must not also appear in the exported GSEA table it belongs to.
untested_leak <- vapply(seq_len(nrow(untested_tbl)), function(i) {
  r <- untested_tbl[i, ]
  tbl <- if (r$analysis_role == "sensitivity") sens_tbl else
    if (r$ontology == "GO_BP") gsea_go_tbl else gsea_cc_tbl
  any(tbl$canonical_comparison == r$canonical_comparison & tbl$ontology == r$ontology &
        tbl$term_id == r$term_id)
}, logical(1))
if (any(untested_leak)) {
  stop(sum(untested_leak), " gene set(s) listed as untested also appear in a GSEA table.",
       call. = FALSE)
}
kegg_untested <- do.call(rbind, lapply(comparisons, function(cmp) {
  w <- release_fgsea_warning_counts(file.path(KEGG_ROOT, "audits", cmp, "warnings.csv"))
  if (!nrow(w)) return(NULL)
  data.frame(canonical_comparison = cmp, analysis = w$analysis, n_sets = w$n_sets,
             stringsAsFactors = FALSE)
}))
untested_n_block <- function(cmp, analysis) {
  if (startsWith(analysis, "KEGG")) {
    return(as.integer(sum(kegg_untested$n_sets[kegg_untested$canonical_comparison == cmp &
                                                  kegg_untested$analysis == analysis])))
  }
  as.integer(sum(untested_sets$canonical_comparison == cmp & untested_sets$analysis == analysis))
}
w <- release_write_table(untested_tbl, release_path("enrichment", "gsea_untested_gene_sets.tsv"))
untested_sources <- unlist(lapply(RELEASE_STATISTICS_RUN$enrichment_runs, function(run) {
  file.path(release_enrichment_root(run, DATA_ROOT), "audits", comparisons,
            "gsea_untested_sets.csv")
}))
release_register("enrichment/gsea_untested_gene_sets.tsv",
                 paste("GO-BP and GO-CC gene sets fgsea gave no P value (canonical and",
                       "effect-size-ranked GSEA), one row per set; absent from the GSEA tables"),
                 c(untested_sources, paths$enrichment_index, paths$enrichment_index_cc),
                 c(vapply(untested_sources, release_sha256, character(1), USE.NAMES = FALSE),
                   index_hash, index_cc_hash), STAGE, "tsv")
release_log("  wrote gsea_untested_gene_sets.tsv (", w$rows, "x", w$cols, "; ",
            sum(untested_tbl$analysis_role == "canonical"), " canonical: GO-BP ",
            sum(untested_tbl$analysis_role == "canonical" & untested_tbl$ontology == "GO_BP"),
            ", GO-CC ",
            sum(untested_tbl$analysis_role == "canonical" & untested_tbl$ontology == "GO_CC"),
            ")")

# Coverage. Several ORA query lists are legitimately EMPTY: with no FDR-significant
# proteins in a comparison there is no query list to over-represent, so the canonical run
# produced a zero-row table. Recording it stops an absent block being read as
# a missing file.
coverage <- do.call(rbind, lapply(comparisons, function(cmp) {
  ix <- enrich_index[as.character(enrich_index$canonical_comparison) == cmp, , drop = FALSE]
  ixc <- enrich_index_cc[as.character(enrich_index_cc$canonical_comparison) == cmp, ,
                         drop = FALSE]
  ixk <- enrich_index_kegg[as.character(enrich_index_kegg$canonical_comparison) == cmp, ,
                           drop = FALSE]
  count_from <- function(index_row, col) {
    if (!col %in% names(index_row)) return(NA_integer_)
    suppressWarnings(as.integer(index_row[[col]][[1]]))
  }
  rows_in <- function(tbl, ...) {
    sel <- rep(TRUE, nrow(tbl)); filters <- list(...)
    for (nm in names(filters)) sel <- sel & tbl[[nm]] == filters[[nm]]
    sum(sel, na.rm = TRUE)
  }
  block <- function(analysis, ontology, query_list, role, n, recorded, untested = NA_integer_) {
    data.frame(canonical_comparison = cmp, analysis = analysis, ontology = ontology,
               query_list = query_list, analysis_role = role, rows_exported = n,
               term_count_recorded_by_run = recorded, gsea_untested_gene_sets = untested,
               stringsAsFactors = FALSE, check.names = FALSE)
  }
  ora_blocks <- function(tbl, variants, ontology, index_row) {
    do.call(rbind, lapply(variants, function(variant) {
      block("ORA", ontology, variant$query_list, variant$analysis_role,
            rows_in(tbl, canonical_comparison = cmp, query_list = variant$query_list),
            count_from(index_row, variant$index_column))
    }))
  }
  rbind(
    block("GSEA", "GO_BP", NA_character_, "canonical",
          rows_in(gsea_go_tbl, canonical_comparison = cmp),
          count_from(ix, "go_gsea_term_count"), untested_n_block(cmp, "GO_BP_GSEA")),
    # KEGG: the 20261124 run's GO_BP folder.
    block("GSEA", "KEGG", NA_character_, "canonical",
          rows_in(gsea_kegg_tbl, canonical_comparison = cmp),
          count_from(ixk, "kegg_gsea_term_count"), untested_n_block(cmp, "KEGG_GSEA")),
    block("GSEA", "GO_BP", NA_character_, "sensitivity",
          rows_in(sens_tbl, canonical_comparison = cmp, ontology = "GO_BP"),
          count_from(ix, "sensitivity_go_gsea_term_count"),
          untested_n_block(cmp, "GO_BP_GSEA_log2fc_sensitivity")),
    block("GSEA", "KEGG", NA_character_, "sensitivity",
          rows_in(sens_tbl, canonical_comparison = cmp, ontology = "KEGG"),
          count_from(ixk, "sensitivity_kegg_gsea_term_count"),
          untested_n_block(cmp, "KEGG_GSEA_log2fc_sensitivity")),
    ora_blocks(ora_tbl, ORA_VARIANTS, "GO_BP", ix),
    # GO-CC run: its index keeps the generic go_* column names.
    block("GSEA", "GO_CC", NA_character_, "canonical",
          rows_in(gsea_cc_tbl, canonical_comparison = cmp),
          count_from(ixc, "go_gsea_term_count"), untested_n_block(cmp, "GO_CC_GSEA")),
    block("GSEA", "GO_CC", NA_character_, "sensitivity",
          rows_in(sens_tbl, canonical_comparison = cmp, ontology = "GO_CC"),
          count_from(ixc, "sensitivity_go_gsea_term_count"),
          untested_n_block(cmp, "GO_CC_GSEA_log2fc_sensitivity")),
    ora_blocks(ora_cc_tbl, ORA_VARIANTS_CC, "GO_CC", ixc)
  )
}))
if (anyNA(coverage$term_count_recorded_by_run)) {
  stop("An enrichment index records no term count for ",
       sum(is.na(coverage$term_count_recorded_by_run)), " coverage block(s).", call. = FALSE)
}
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
                       "term counts the canonical runs recorded (GO-BP, GO-CC and KEGG), with",
                       "the gene sets fgsea could not test per GSEA block"),
                 c(paths$enrichment_index, paths$enrichment_index_cc,
                   paths$enrichment_index_kegg),
                 c(index_hash, index_cc_hash, index_kegg_hash), STAGE, "tsv")
release_log("  wrote enrichment_coverage.tsv (", w$rows, "x", w$cols, ")")

# The canonical run parameters, republished verbatim plus the EWCE settings. Four rows are
# documented rather than recorded: the runs' run_parameters.csv does not capture the seed
# handling, the permutation settings of all runs together, the KEGG gene-set version or the
# GO-CC history. The KEGG run's own parameters describe its GO-BP tables too, which are
# superseded; only its KEGG rows are released.
params <- rbind(
  data.frame(analysis = "enrichment GO-BP (GSEA/ORA)",
             parameter = as.character(run_parameters$parameter),
             value = as.character(run_parameters$value),
             stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(analysis = "enrichment GO-CC (GSEA/ORA)",
             parameter = as.character(run_parameters_cc$parameter),
             value = as.character(run_parameters_cc$value),
             stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(analysis = "enrichment KEGG (GSEA)",
             parameter = as.character(run_parameters_kegg$parameter),
             value = as.character(run_parameters_kegg$value),
             stringsAsFactors = FALSE, check.names = FALSE),
  data.frame(analysis = "enrichment (GSEA/ORA)",
             parameter = c("ora_universe_definition_per_comparison",
                           "gsea_seed_handling (documented)",
                           "gsea_n_perm_simple (documented)",
                           "kegg_gene_set_version (documented)",
                           "go_cc_status (documented)"),
             value = c(paste(ORA_UNIVERSE, collapse = ";"), RELEASE_ENRICHMENT_SEED_HANDLING,
                       permutations$sentence, RELEASE_KEGG_VERSION_NOTE,
                       RELEASE_GO_CC_PROVENANCE_NOTE),
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
                 paste("canonical enrichment and EWCE run parameters, as recorded by the",
                       "runs, plus four documented rows (seed handling, nPermSimple, KEGG",
                       "version, GO-CC)"),
                 c(paths$run_parameters, paths$run_parameters_cc, paths$run_parameters_kegg,
                   paths$ewce_session),
                 c(release_sha256(paths$run_parameters), release_sha256(paths$run_parameters_cc),
                   release_sha256(paths$run_parameters_kegg), release_sha256(paths$ewce_session)),
                 STAGE, "tsv")
release_log("  wrote enrichment_run_parameters.tsv (", w$rows, "x", w$cols, ")")
release_log("stage 05 complete")
