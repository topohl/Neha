#!/usr/bin/env Rscript

# Publication release, stage 07 -- Supplementary Table 2, the proteomics workbook.
#
# Produces
#   editor_source_data/Supplementary_Table_2_proteomics.xlsx
#   editor_source_data/workbook_sheet_inventory.tsv
#
# The workbook is Supplementary Table 2 as a whole. A Contents sheet (with links) comes
# first, then README, Column definitions and the lettered sub-tables 2a-2q, in the order and
# with the titles of RELEASE_SHEET_META (R/release_presentation.R). Every sheet: title
# (row 1), caption (row 2, merged), the values shared by every row (row 3, when any), header
# (row 4), data from row 5; panes frozen below the header and right of column 1.
#
# Nothing is computed. Rows are copied from the release tables, filtered on values already in
# them (analysis_role, id_mapping_status), sorted, and keys are shown as display labels.
# Comparison-keyed sheets show the display comparison first and the machine key last. The
# effect-size-ranked GSEA, the alternative ORA list, the metadata field status, internal
# paths, hashes, curation and revision fields stay in the release files only; the README
# sheet names the files. Columns blank on every row and exact duplicates are left out.

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
release_require("digest", "openxlsx")

DATA_ROOT <- release_data_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/07_build_editor_source_workbook.R"
XLSX_REL <- "editor_source_data/Supplementary_Table_2_proteomics.xlsx"

release_banner("stage 07 -- Supplementary Table 2 workbook")

inv <- RELEASE_DESIGN_INVARIANTS
fmt_n <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
# Thousands separators for the counts in free text (tables 2o and 2p): whole numbers of 4 to
# 8 digits, and both ends of a numeric range ('3882-5997' -> '3,882-5,997'). Left alone: dates
# (yyyy-mm-dd, yyyymmdd), decimals, identifiers and hashes (a letter, underscore, dot, colon or
# slash next to the digits), a 4-digit year (1900-2099) or year range, numbers of 9 or more
# digits (seeds, identifiers) and numbers already grouped.
group_thousands_one <- function(s) {
  if (is.na(s) || !nzchar(s)) return(s)
  runs <- gregexpr("[0-9]+", s)[[1]]
  if (runs[1] == -1L) return(s)
  st <- as.integer(runs)
  en <- st + attr(runs, "match.length") - 1L
  val <- substring(s, st, en)
  chars <- strsplit(s, "")[[1]]
  ch <- function(i) if (i >= 1L && i <= length(chars)) chars[i] else ""
  is_digit <- function(c) grepl("^[0-9]$", c)
  glued_left <- vapply(st, function(p) {
    b <- ch(p - 1L)
    grepl("^[A-Za-z_.:/]$", b) || (b == "," && is_digit(ch(p - 2L))) ||
      (b == "-" && grepl("^[A-Za-z_.:/]$", ch(p - 2L)))
  }, logical(1))
  glued_right <- vapply(en, function(p) {
    a <- ch(p + 1L)
    grepl("^[A-Za-z_/]$", a) || (a %in% c(".", ",", ":") && grepl("^[0-9A-Za-z]$", ch(p + 2L)))
  }, logical(1))
  # Runs joined by single hyphens form one token: a range (2 runs) or a date/identifier (3+).
  linked <- c(FALSE, st[-1L] == en[-length(en)] + 2L &
                vapply(en[-length(en)], function(p) ch(p + 1L) == "-", logical(1)))
  token <- cumsum(!linked)
  year <- function(v) nchar(v) == 4L && as.numeric(v) >= 1900 && as.numeric(v) <= 2099
  ymd <- function(v) grepl("^(19|20)[0-9]{2}(0[1-9]|1[0-2])(0[1-9]|[12][0-9]|3[01])$", v)
  count_like <- function(v) nchar(v) >= 4L && nchar(v) <= 8L && !ymd(v)
  group <- logical(length(st))
  for (tk in unique(token)) {
    i <- which(token == tk)
    if (glued_left[i[1]] || glued_right[i[length(i)]] || length(i) > 2L) next
    if (length(i) == 2L && as.numeric(val[i[1]]) <= as.numeric(val[i[2]])) {
      if (!(year(val[i[1]]) && year(val[i[2]]))) group[i] <- vapply(val[i], count_like, logical(1))
    } else {
      group[i] <- vapply(val[i], function(v) count_like(v) && !year(v), logical(1))
    }
  }
  for (k in rev(which(group))) {
    s <- paste0(substr(s, 1L, st[k] - 1L), fmt_n(val[k]), substr(s, en[k] + 1L, nchar(s)))
  }
  s
}
group_thousands <- function(x) {
  vapply(as.character(x), group_thousands_one, character(1), USE.NAMES = FALSE)
}
local({
  cases <- c("3882-5997" = "3,882-5,997", "1-1324" = "1-1,324", "2024-12-17" = "2024-12-17",
             "in 2026" = "in 2026", "2024-2026" = "2024-2026", "20241217" = "20241217",
             "seed 946523163" = "seed 946523163", "5349 proteins" = "5,349 proteins",
             "0 of 5349." = "0 of 5,349.", "all 5,349" = "all 5,349", "0.0001" = "0.0001",
             "C46" = "C46", "GO:0007267" = "GO:0007267", "mmu04724" = "mmu04724",
             "S4-A1_1_8782" = "S4-A1_1_8782", "2026-11" = "2026-11")
  got <- group_thousands(names(cases))
  if (!identical(unname(got), unname(cases))) {
    stop("group_thousands(): ", paste0("'", names(cases)[got != cases], "' -> '",
                                       got[got != cases], "'", collapse = "; "))
  }
})

rel <- function(...) release_path(..., create_dir = FALSE)
read_release <- function(path) {
  release_assert_exists(path, basename(path))
  if (grepl("[.]gz$", path)) {
    read.delim(gzfile(path), sep = "\t", quote = "", stringsAsFactors = FALSE,
               check.names = FALSE)
  } else {
    read.delim(path, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  }
}

sample_metadata <- read_release(rel("metadata", "sample_metadata.tsv"))
animal_metadata <- read_release(rel("metadata", "animal_level_sample_metadata.tsv"))
primary_contrasts <- read_release(rel("metadata", "primary_contrast_manifest.tsv"))
secondary <- read_release(rel("metadata", "secondary_analysis_manifest.tsv"))
experimenter_metadata_wb <- read_release(rel("metadata", "experimenter_metadata.tsv"))
mv_wb <- function(id) release_metadata_value(experimenter_metadata_wb, id)
differential <- read_release(rel("differential_analysis", "primary_differential_proteins.tsv.gz"))
diff_summary <- read_release(rel("differential_analysis", "primary_differential_summary.tsv"))
sig_proteins <- read_release(rel("differential_analysis", "significant_proteins.tsv"))
gsea_go <- read_release(rel("enrichment", "primary_GSEA_GO_BP.tsv.gz"))
gsea_cc <- read_release(rel("enrichment", "primary_GSEA_GO_CC.tsv.gz"))
gsea_kegg <- read_release(rel("enrichment", "primary_GSEA_KEGG.tsv.gz"))
ora <- read_release(rel("enrichment", "primary_ORA_GO_BP.tsv.gz"))
ora_cc <- read_release(rel("enrichment", "primary_ORA_GO_CC.tsv.gz"))
sig_sets <- read_release(rel("enrichment", "significant_gene_sets.tsv.gz"))
untested_sets <- read_release(rel("enrichment", "gsea_untested_gene_sets.tsv"))
ewce <- read_release(rel("enrichment", "primary_EWCE.tsv.gz"))
coverage <- read_release(rel("enrichment", "enrichment_coverage.tsv"))
run_params <- read_release(rel("enrichment", "enrichment_run_parameters.tsv"))
figure_map <- read_release(rel("provenance", "figure_source_map.tsv"))
for (col in names(figure_map)[vapply(figure_map, is.character, logical(1))]) {
  figure_map[[col]][is.na(figure_map[[col]])] <- ""
}

# Software versions as the runs recorded them (the rows of provenance/software_versions.tsv).
sw <- release_build_software_versions(DATA_ROOT, REPO_ROOT)
stat_run <- release_read_stat_run(DATA_ROOT)
ewce_model <- release_ewce_differential_model(DATA_ROOT)
permutations <- release_enrichment_permutations(DATA_ROOT)

# Retained for the registry provenance recorded at the end of this stage.
pkg_versions_path <- file.path(release_enrichment_root("GO_BP", DATA_ROOT), "audits",
                               "package_database_versions.csv")
ewce_session_path <- file.path(DATA_ROOT, "03_output", "ewce",
                               "EWCE_Results_animal_level_validation_20260825",
                               "03_QC_Mapping_Logs", "reproducibility_session_info.txt")
release_assert_exists(pkg_versions_path, "canonical package_database_versions.csv")
release_assert_exists(ewce_session_path, "canonical EWCE sessionInfo")

# --------------------------------------------------------------------------------------
# comparisons: display label and sort order
# --------------------------------------------------------------------------------------

cmp <- primary_contrasts[, c("canonical_comparison", "sample_class", "contrast_family",
                             "numerator_condition", "denominator_condition")]
cmp$label <- release_comparison_label(cmp$sample_class, cmp$numerator_condition,
                                      cmp$denominator_condition)
cmp$rank <- release_comparison_rank(cmp$sample_class, cmp$contrast_family)
if (nrow(cmp) != inv$n_primary_contrasts || anyNA(cmp$rank) || anyDuplicated(cmp$rank)) {
  stop("The primary contrast manifest does not give 12 comparisons in known families and ",
       "compartments.", call. = FALSE)
}
cmp_field <- function(keys, field) cmp[[field]][match(as.character(keys), cmp$canonical_comparison)]
cmp_rank <- function(keys) cmp_field(keys, "rank")

#' Add the display comparison; rows without a comparison (EWCE baseline) get "".
add_comparison_label <- function(d, sheet) {
  keys <- as.character(d$canonical_comparison)
  keys[is.na(keys)] <- ""
  unknown <- setdiff(unique(keys[nzchar(keys)]), cmp$canonical_comparison)
  if (length(unknown)) {
    stop("Sheet ", sheet, " has unknown comparison key(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  d$canonical_comparison <- keys
  d$comparison_label <- ifelse(nzchar(keys), cmp_field(keys, "label"), "")
  d
}

#' Display labels for controlled-vocabulary values (values only; numbers untouched).
display_values <- function(d) {
  if ("direction" %in% names(d)) {
    keys <- if ("canonical_comparison" %in% names(d)) d$canonical_comparison else rep("", nrow(d))
    d$direction <- release_direction_display(d$direction, cmp_field(keys, "numerator_condition"),
                                             cmp_field(keys, "denominator_condition"))
  }
  for (col in intersect(c("sample_class", "original_sample_class"), names(d))) {
    d[[col]] <- release_sample_class_display(d[[col]])
  }
  for (col in intersect(c("condition", "numerator_condition", "denominator_condition"),
                        names(d))) {
    d[[col]] <- release_condition_display(d[[col]])
  }
  maps <- list(treatment = RELEASE_TREATMENT_DISPLAY, contrast_family = RELEASE_FAMILY_LABEL,
               query_list = RELEASE_QUERY_LIST_DISPLAY,
               analysis_role = RELEASE_ANALYSIS_ROLE_DISPLAY,
               rank_statistic = RELEASE_RANK_STATISTIC_DISPLAY,
               ontology = RELEASE_ONTOLOGY_DISPLAY)
  for (col in intersect(names(maps), names(d))) {
    d[[col]] <- release_lookup_display(d[[col]], maps[[col]])
  }
  # Shown as "Annotated", an Excel boolean: TRUE when the protein group has a mouse UniProt
  # accession (id_mapping_status mapped), as in the annotated counts of tables 2c and 2d.
  if ("id_mapping_status" %in% names(d)) {
    if (!all(d$id_mapping_status %in% c("mapped", "unmapped"))) {
      stop("id_mapping_status holds values other than mapped and unmapped.", call. = FALSE)
    }
    d$id_mapping_status <- d$id_mapping_status == "mapped"
  }
  d
}

hit_cmps <- diff_summary$canonical_comparison[as.integer(diff_summary$n_significant_fdr_0_05) > 0L]
hit_cmps <- hit_cmps[order(cmp_rank(hit_cmps))]
labels_text <- function(keys) {
  l <- cmp_field(keys, "label")
  if (length(l) <= 2L) paste(l, collapse = " and ") else
    paste0(paste(l[-length(l)], collapse = ", "), " and ", l[length(l)])
}

# --------------------------------------------------------------------------------------
# 2a-2c: samples, animal-level profiles, comparisons
# --------------------------------------------------------------------------------------

# Figure order, as in tables 2c-2n: compartment (neuropil, c-Fos+, mCherry+, NeuN+), then
# group (paired VEH, unpaired VEH, paired CNO, unpaired CNO), animal and hemisphere.
SAMPLE_CONDITION_ORDER <- c("paired_veh", "unpaired_veh", "paired_cno", "unpaired_cno")
sample_order <- function(d) {
  o <- order(match(d$sample_class, names(RELEASE_SAMPLE_CLASS_DISPLAY)),
             match(d$condition, SAMPLE_CONDITION_ORDER), d$AnimalID,
             if ("hemisphere" %in% names(d)) match(d$hemisphere, c("Left", "Right")) else
               rep(0L, nrow(d)))
  if (anyNA(match(d$sample_class, names(RELEASE_SAMPLE_CLASS_DISPLAY))) ||
      anyNA(match(d$condition, SAMPLE_CONDITION_ORDER))) {
    stop("Unknown compartment or condition in the sample metadata.", call. = FALSE)
  }
  d[o, , drop = FALSE]
}
samples <- sample_order(sample_metadata)
animal_profiles <- sample_order(animal_metadata)

shared_cols <- c("sample_class", "contrast_family", "numerator_condition",
                 "denominator_condition", "n_numerator_animals", "n_denominator_animals")
count_cols <- c("n_proteins_tested", "n_significant_fdr_0_05",
                "n_significant_higher_in_numerator", "n_significant_higher_in_denominator",
                "n_proteins_tested_mapped_only", "n_significant_fdr_0_05_mapped_only",
                "n_significant_higher_in_numerator_mapped_only",
                "n_significant_higher_in_denominator_mapped_only")
ds <- diff_summary[match(primary_contrasts$canonical_comparison,
                         diff_summary$canonical_comparison), , drop = FALSE]
if (anyNA(ds$canonical_comparison) ||
    !identical(lapply(primary_contrasts[shared_cols], as.character),
               lapply(ds[shared_cols], as.character))) {
  stop("The contrast manifest and the differential summary disagree on the comparisons.",
       call. = FALSE)
}
comparisons <- cbind(primary_contrasts[, c("canonical_comparison", shared_cols)],
                     ds[, count_cols])
comparisons <- comparisons[order(cmp_rank(comparisons$canonical_comparison)), , drop = FALSE]

# --------------------------------------------------------------------------------------
# 2d, 2e: proteins
# --------------------------------------------------------------------------------------

n_mapped_sig <- sum(as.integer(diff_summary$n_significant_fdr_0_05_mapped_only))
if (nrow(sig_proteins) != n_mapped_sig) {
  stop("significant_proteins.tsv has ", nrow(sig_proteins), " rows; the differential summary ",
       "counts ", n_mapped_sig, " annotated proteins at FDR < 0.05.", call. = FALSE)
}
sig_proteins <- sig_proteins[order(cmp_rank(sig_proteins$canonical_comparison),
                                   sig_proteins$adj.P.Val, sig_proteins$P.Value,
                                   sig_proteins$protein_group_id), , drop = FALSE]

# The non-mouse protein groups that pass FDR < 0.05 are named in the 2d caption.
NONMOUSE_SIGNIFICANT <- c(TRYP_PIG = "porcine trypsin", API_ACHLY = "Lys-C protease",
                          TRFE_BOVIN = "bovine serotransferrin",
                          APOA1_BOVIN = "bovine apolipoprotein A-I",
                          K1C16_HUMAN = "human keratin 16")
nonmouse <- unique(differential$protein_group_id[differential$significant_fdr_0_05 %in% TRUE &
                                                   differential$id_mapping_status != "mapped"])
if (!setequal(nonmouse, names(NONMOUSE_SIGNIFICANT))) {
  stop("The non-mouse protein groups at FDR < 0.05 are now ", paste(nonmouse, collapse = ", "),
       "; update NONMOUSE_SIGNIFICANT and the table 2d caption.", call. = FALSE)
}
sig_cmps <- unique(sig_proteins$canonical_comparison)
caption_2d <- paste0(
  "Annotated proteins at FDR < 0.05 (", fmt_n(nrow(sig_proteins)), ", ",
  if (length(sig_cmps) == 1L) paste0("all in ", cmp_field(sig_cmps, "label")) else
    paste0("in ", length(sig_cmps), " comparisons"),
  "), sorted by comparison and adjusted P value. ", length(nonmouse),
  " non-mouse protein groups (",
  paste(NONMOUSE_SIGNIFICANT[names(NONMOUSE_SIGNIFICANT) %in% nonmouse], collapse = ", "),
  ") also pass FDR < 0.05; they are listed in table 2e but not here. The effect size is a",
  " difference in standardized abundance in SD units (group 1 minus group 2), not a log2",
  " fold change.")

all_proteins <- differential[order(cmp_rank(differential$canonical_comparison),
                                   differential$P.Value, differential$protein_group_id), ,
                             drop = FALSE]

# --------------------------------------------------------------------------------------
# 2f-2l: gene sets
# --------------------------------------------------------------------------------------

ONTOLOGY_ORDER <- names(RELEASE_ONTOLOGY_DISPLAY)
for (nm in c("gsea_go", "gsea_cc", "gsea_kegg")) {
  if (!all(get(nm)$analysis_role == "canonical")) {
    stop(nm, " holds rows that are not the primary GSEA.", call. = FALSE)
  }
}
n_sig_sets_source <- sum(c(gsea_go$adjusted_p_value, gsea_cc$adjusted_p_value,
                           gsea_kegg$adjusted_p_value) < 0.05, na.rm = TRUE)
if (nrow(sig_sets) != n_sig_sets_source) {
  stop("significant_gene_sets.tsv.gz has ", nrow(sig_sets), " rows; the GSEA tables have ",
       n_sig_sets_source, " rows at adjusted P < 0.05.", call. = FALSE)
}
sig_sets <- sig_sets[order(cmp_rank(sig_sets$canonical_comparison), sig_sets$adjusted_p_value,
                           sig_sets$p_value, match(sig_sets$ontology, ONTOLOGY_ORDER),
                           sig_sets$term_id), , drop = FALSE]
n_by_ontology <- vapply(ONTOLOGY_ORDER, function(o) sum(sig_sets$ontology == o), integer(1))
no_sets <- cmp$canonical_comparison[!cmp$canonical_comparison %in% sig_sets$canonical_comparison]
no_sets <- no_sets[order(cmp_rank(no_sets))]
caption_2f <- paste0(
  "Gene sets at FDR < 0.05 in the GSEA of tables 2g to 2i (",
  paste(fmt_n(n_by_ontology), RELEASE_ONTOLOGY_DISPLAY[ONTOLOGY_ORDER], collapse = ", "),
  "), with their leading-edge proteins as gene symbols. ",
  if (length(no_sets)) paste0(labels_text(no_sets),
                              if (length(no_sets) == 1L) " has none. " else " have none. ")
  else "",
  "Proteins were ranked by the moderated t. Sorted by comparison and adjusted P value. P",
  " values at fgsea's lower bound (1E-10) are shown as <=1E-10.")

sort_gsea <- function(d) {
  d[order(cmp_rank(d$canonical_comparison), d$adjusted_p_value, d$p_value, d$term_id), ,
    drop = FALSE]
}
gsea_go <- sort_gsea(gsea_go)
gsea_cc <- sort_gsea(gsea_cc)
gsea_kegg <- sort_gsea(gsea_kegg)

FDR_LISTS <- c("fdr_significant_all", "fdr_significant_higher_in_numerator",
               "fdr_significant_higher_in_denominator")
primary_ora <- function(d, what) {
  d <- d[d$analysis_role == "canonical", , drop = FALSE]
  if (!all(d$query_list %in% FDR_LISTS)) stop(what, ": unexpected primary query list.", call. = FALSE)
  extra <- setdiff(unique(d$canonical_comparison), hit_cmps)
  if (length(extra)) {
    stop(what, " has rows for comparisons without proteins at FDR < 0.05: ",
         paste(extra, collapse = ", "), call. = FALSE)
  }
  d[order(cmp_rank(d$canonical_comparison), match(d$query_list, FDR_LISTS),
          d$adjusted_p_value, d$p_value, d$term_id), , drop = FALSE]
}
ora <- primary_ora(ora, "ORA GO-BP")
ora_cc <- primary_ora(ora_cc, "ORA GO-CC")
ora_caption <- function(collection) {
  paste0("Over-representation of ", collection, " terms among the proteins at FDR < 0.05",
         " (all, higher in group 1, lower in group 1), tested against all measured annotated",
         " proteins (background in table 2m). ",
         if (length(hit_cmps) == 1L) paste0("Only ", cmp_field(hit_cmps, "label"),
                                            " has proteins at FDR < 0.05, so only it has rows. ")
         else "",
         "Every tested term is listed: filter on the adjusted P value.")
}

# Table 2l: the primary GSEA only; the rows of the effect-size-ranked GSEA stay in the file.
untested_sets <- untested_sets[untested_sets$analysis_role == "canonical", , drop = FALSE]
untested_sets <- untested_sets[order(cmp_rank(untested_sets$canonical_comparison),
                                     match(untested_sets$ontology, ONTOLOGY_ORDER),
                                     untested_sets$term_name), , drop = FALSE]

# --------------------------------------------------------------------------------------
# 2m: enrichment settings and coverage
# --------------------------------------------------------------------------------------

PARAM_BLOCK <- c(GO_BP = "enrichment GO-BP (GSEA/ORA)", GO_CC = "enrichment GO-CC (GSEA/ORA)",
                 KEGG = "enrichment KEGG (GSEA)")
param <- function(ontology, name) {
  v <- run_params$value[run_params$analysis == PARAM_BLOCK[[ontology]] &
                          run_params$parameter == name]
  if (length(v) != 1L) {
    stop("enrichment_run_parameters.tsv: no single ", name, " for ", ontology, call. = FALSE)
  }
  v
}
sw_enrich <- sw[sw$applies_to_stage == "differential enrichment (GSEA / ORA)", , drop = FALSE]
swv <- function(component) {
  v <- unique(sw_enrich$version[sw_enrich$component == component])
  if (length(v) != 1L) stop("No single recorded version of ", component, call. = FALSE)
  v
}
GENE_SET_SOURCE <- c(
  GO_BP = paste0("GO biological process (GO.db ", swv("GO.db"), "; mouse annotation ",
                 "org.Mm.eg.db ", swv("org.Mm.eg.db"), ")"),
  GO_CC = paste0("GO cellular component (GO.db ", swv("GO.db"), "; mouse annotation ",
                 "org.Mm.eg.db ", swv("org.Mm.eg.db"), ")"),
  KEGG = paste0("KEGG pathways, Mus musculus (KEGGREST ", swv("KEGGREST"), "), downloaded ",
                RELEASE_STATISTICS_RUN$kegg_enrichment_date, "; KEGG release not recorded"))
ADJUST_METHOD <- c(BH = "Benjamini-Hochberg")
permutation_text <- function(ontology) {
  if (ontology %in% c("GO_BP", "GO_CC")) {
    n <- param(ontology, "gsea_n_perm_simple")
    if (!identical(as.integer(n), permutations$go_n_perm_simple)) {
      stop("The GO nPermSimple of the run parameters and of the run audit differ.", call. = FALSE)
    }
    return(paste0("nPermSimple = ", fmt_n(n)))
  }
  paste0("fgsea default", if (!is.na(permutations$default_n_perm_simple))
    paste0(" (nPermSimple = ", fmt_n(permutations$default_n_perm_simple), ")") else "")
}
ora_tables <- list(GO_BP = ora, GO_CC = ora_cc)
ora_background <- function(key, ontology, query_list) {
  d <- ora_tables[[ontology]]
  r <- d$background_ratio[d$canonical_comparison == key & d$query_list == query_list]
  if (!length(r)) return(NA_integer_)
  u <- unique(sub("^.*/", "", r))
  if (length(u) != 1L) stop("ORA background differs within one block: ", key, call. = FALSE)
  as.integer(u)
}
LIST_COUNT <- c(fdr_significant_all = "n_significant_fdr_0_05_mapped_only",
                fdr_significant_higher_in_numerator = "n_significant_higher_in_numerator_mapped_only",
                fdr_significant_higher_in_denominator = "n_significant_higher_in_denominator_mapped_only")

cov <- coverage[coverage$analysis_role == "canonical", , drop = FALSE]
if (!all(as.logical(cov$agrees_with_run_record) %in% TRUE)) {
  stop("enrichment_coverage.tsv records a block whose rows differ from the run record.",
       call. = FALSE)
}
cov$query_list[is.na(cov$query_list)] <- ""
is_gsea <- cov$analysis == "GSEA"
if (any(is_gsea & is.na(cov$gsea_seed))) {
  stop("enrichment_coverage.tsv lacks the GSEA seed of a block.", call. = FALSE)
}
enrichment_settings <- data.frame(
  canonical_comparison = cov$canonical_comparison,
  sample_class = cmp_field(cov$canonical_comparison, "sample_class"),
  analysis = cov$analysis,
  ontology = cov$ontology,
  query_list = cov$query_list,
  gene_set_source = unname(GENE_SET_SOURCE[cov$ontology]),
  set_size_limits = vapply(cov$ontology, function(o) {
    paste(param(o, "min_gs_size"), "to", param(o, "max_gs_size"))
  }, character(1)),
  rank_statistic = ifelse(is_gsea, vapply(cov$ontology, param, character(1),
                                          name = "ranking_statistic"), "not_applicable"),
  permutations = ifelse(is_gsea, vapply(cov$ontology, permutation_text, character(1)), ""),
  seed = ifelse(is_gsea, as.character(cov$gsea_seed), ""),
  ora_background = vapply(seq_len(nrow(cov)), function(i) {
    if (is_gsea[[i]]) NA_integer_ else
      ora_background(cov$canonical_comparison[[i]], cov$ontology[[i]], cov$query_list[[i]])
  }, integer(1)),
  multiple_testing = unname(ADJUST_METHOD[vapply(cov$ontology, param, character(1),
                                                 name = "p_adjust_method")]),
  rows_exported = as.integer(cov$rows_exported),
  gsea_untested_gene_sets = suppressWarnings(as.integer(cov$gsea_untested_gene_sets)),
  note = vapply(seq_len(nrow(cov)), function(i) {
    if (is_gsea[[i]] || as.integer(cov$rows_exported[[i]]) > 0L) return("")
    n <- as.integer(diff_summary[[LIST_COUNT[[cov$query_list[[i]]]]]][
      match(cov$canonical_comparison[[i]], diff_summary$canonical_comparison)])
    if (n == 0L) "Empty query list: no protein at FDR < 0.05." else
      "No gene set contains a query protein."
  }, character(1)),
  stringsAsFactors = FALSE, check.names = FALSE)
if (anyNA(enrichment_settings$multiple_testing)) {
  stop("Unexpected multiple-testing method in the run parameters.", call. = FALSE)
}
enrichment_settings <- enrichment_settings[order(
  cmp_rank(enrichment_settings$canonical_comparison),
  match(enrichment_settings$analysis, c("GSEA", "ORA")),
  match(enrichment_settings$ontology, ONTOLOGY_ORDER),
  match(enrichment_settings$query_list, c("", FDR_LISTS))), , drop = FALSE]

# --------------------------------------------------------------------------------------
# 2n: EWCE
# --------------------------------------------------------------------------------------

ewce_param <- function(name) {
  v <- run_params$value[run_params$analysis == "EWCE" & run_params$parameter == name]
  if (length(v) != 1L) stop("No single EWCE ", name, " in the run parameters.", call. = FALSE)
  v
}
is_base <- ewce$ewce_analysis_type == "Baseline"
ewce_keys <- ifelse(is.na(ewce$canonical_comparison), "", ewce$canonical_comparison)
list_ok <- ifelse(is_base, ewce$protein_list_definition == ewce$numerator_condition &
                    ewce$direction == "abundant",
                  ewce$protein_list_definition == paste0(ewce_keys, "_", ewce$direction) &
                    ewce$direction %in% c("up", "down"))
if (!all(list_ok %in% TRUE)) {
  stop("EWCE protein lists do not follow <condition> (Baseline) or <comparison>_<up|down>.",
       call. = FALSE)
}
ewce_tbl <- ewce
ewce_tbl$protein_list_definition <- ifelse(
  is_base, paste0("most abundant, ", release_condition_display(ewce$numerator_condition)),
  paste("higher in", release_condition_display(ifelse(
    ewce$direction == "up", cmp_field(ewce_keys, "numerator_condition"),
    cmp_field(ewce_keys, "denominator_condition")))))
level1 <- as.character(ewce$annotation_level) == "1"
if (!all(ewce$cell_type[level1] %in% names(RELEASE_CELL_TYPE_DISPLAY))) {
  stop("Level-1 EWCE cell type(s) without a display label: ",
       paste(setdiff(ewce$cell_type[level1], names(RELEASE_CELL_TYPE_DISPLAY)), collapse = ", "),
       call. = FALSE)
}
ewce_tbl$cell_type <- ifelse(ewce$cell_type == "(none)",
                             paste0("unassigned level-", ewce$annotation_level, " cluster"),
                             ifelse(level1, release_lookup_display(ewce$cell_type,
                                                                   RELEASE_CELL_TYPE_DISPLAY),
                                    ewce$cell_type))
ewce_tbl <- ewce_tbl[order(
  !is_base,
  ifelse(is_base,
         match(ewce$sample_class, names(RELEASE_SAMPLE_CLASS_DISPLAY)) * 10L +
           match(ewce$numerator_condition, names(RELEASE_CONDITION_DISPLAY)),
         cmp_rank(ewce_keys)),
  -as.integer(ewce$is_primary_setting), ewce$top_n, ewce$annotation_level,
  match(ewce$direction, c("abundant", "up", "down")), ewce$adjusted_p_value_global,
  ewce$p_value, ewce$cell_type), , drop = FALSE]
spaced <- function(x) gsub("\\s*([=+])\\s*", " \\1 ", x)
ewce_fit <- release_reader_text(paste0(spaced(sub(";.*$", "", ewce_model$design)), "; eBayes ",
                                       gsub("\\s*;\\s*", ", ", spaced(ewce_model$ebayes))), cmp)
caption_2n <- paste0(
  "Expression-weighted cell-type enrichment (EWCE) of protein lists against the ewceData",
  " mouse cortex and hippocampus single-cell reference (CTD), with ",
  fmt_n(ewce_param("bootstrap_replicates")), " bootstrap replicates against the measured",
  " proteome. Primary setting: annotation level ", ewce_param("primary_annotation_level"),
  ", top ", ewce_param("primary_top_n"), " proteins; other rows are sensitivity settings.",
  " Baseline rows test the most abundant proteins of one condition and have no comparison.",
  " The differential lists come from EWCE's own limma fit per compartment (", ewce_fit,
  ") over the ", fmt_n(ewce_model$n_proteins), " proteins in the reference, not from the",
  " statistics of table 2e. Level-2 cell types are the reference's cluster names (Zeisel et",
  " al. 2015; see the ewceData documentation). P values of 0 (below the bootstrap",
  " resolution) are shown as <0.0001.")

# --------------------------------------------------------------------------------------
# 2o: secondary analyses reported in the paper
# --------------------------------------------------------------------------------------

panel_files <- function(label) {
  p <- figure_map$panel_source_data_release_path[figure_map$panel_label == label]
  if (length(p) != 1L || !nzchar(p)) stop("No Source Data file for ", label, call. = FALSE)
  basename(strsplit(p, ";", fixed = TRUE)[[1]])
}
panel_key_numbers <- function(label) {
  k <- figure_map$key_numbers_revised[figure_map$panel_label == label]
  if (length(k) != 1L) stop("No key numbers for ", label, call. = FALSE)
  k
}
sec_row <- function(id) {
  r <- secondary[secondary$analysis_id == id, , drop = FALSE]
  if (nrow(r) != 1L) stop("Secondary analysis missing from the manifest: ", id, call. = FALSE)
  r
}
sentences_with <- function(text, pattern) {
  s <- release_split_sentences(text)
  s <- s[grepl(pattern, s)]
  if (!length(s)) stop("Expected result text not found: ", pattern, call. = FALSE)
  paste(s, collapse = " ")
}
default_gsea <- paste0("GSEA at fgsea's default permutations",
                       if (!is.na(permutations$default_n_perm_simple))
                         paste0(" (nPermSimple = ", fmt_n(permutations$default_n_perm_simple),
                                ")") else "")

supp_e <- sec_row("supp_E_cross_compartment_paired")
supp_e_cmps <- release_reader_text(strsplit(supp_e$comparisons, ";", fixed = TRUE)[[1]], cmp)
supp_e_key <- panel_key_numbers("Supp. Fig. 4 E")
count_list <- function(pattern) {
  m <- regmatches(supp_e_key, regexec(pattern, supp_e_key))[[1]]
  if (length(m) != 2L) stop("Supp. Fig. 4 E key numbers lack '", pattern, "'.", call. = FALSE)
  v <- strsplit(m[[2]], "/", fixed = TRUE)[[1]]
  if (length(v) != length(supp_e_cmps)) {
    stop("Supp. Fig. 4 E counts do not match its comparisons.", call. = FALSE)
  }
  paste0(supp_e_cmps, " ", v, collapse = ", ")
}
rank_sens <- sec_row("gsea_log2fc_rank_sensitivity")
jaccard <- regmatches(rank_sens$interpretation,
                      regexpr("GO-BP [0-9.]+-[0-9.]+, GO-CC [0-9.]+-[0-9.]+",
                              rank_sens$interpretation))
if (length(jaccard) != 1L) stop("Rank-sensitivity agreement not found.", call. = FALSE)
interaction <- sec_row("posthoc_pairing_x_cno_interaction")
loo <- sec_row("mcherry_learning_leave_one_animal_out")
ewce_base <- sec_row("ewce_baseline_cell_type_identity")
ewce_key <- trimws(strsplit(panel_key_numbers("Supp. Fig. 4 C"), ";", fixed = TRUE)[[1]][[1]])

# The sensitivity analyses ran at fgsea's default permutations; their GSEA counts are set
# against the primary GSEA of the same comparison, counted from tables 2g and 2h.
default_perm_text <- if (!is.na(permutations$default_n_perm_simple))
  paste0("nPermSimple = ", fmt_n(permutations$default_n_perm_simple)) else
    "fgsea's default permutations"
MCH_PAIRING <- "mcherry_paired_veh_over_mcherry_unpaired_veh"
n_fdr_sets <- function(d) {
  sum(d$canonical_comparison == MCH_PAIRING & d$adjusted_p_value < 0.05, na.rm = TRUE)
}
loo_primary_text <- paste0(
  "GSEA counts come from runs at ", default_perm_text, "; the primary GSEA (tables 2f to 2h, ",
  "nPermSimple = ", fmt_n(permutations$go_n_perm_simple), ") gives ", n_fdr_sets(gsea_go),
  " GO-BP and ", n_fdr_sets(gsea_cc), " GO-CC gene sets at FDR < 0.05 for this comparison.")
model_text <- function(x) release_model_words(release_reader_text(x, cmp))

secondary_analyses <- data.frame(
  analysis_name = c(
    "Cross-compartment comparison (Supp. Fig. 4 E)",
    "GSEA ranked by effect size (sensitivity analysis)",
    "Pairing x CNO interaction (post hoc)",
    paste0("Leave one animal out, ", cmp_field(MCH_PAIRING, "label")),
    "EWCE baseline cell-type profile (Supp. Fig. 4 C)"),
  purpose = c(
    "Compare the protein profiles of the compartments within the same animals.",
    "Check whether the GSEA results depend on the ranking statistic.",
    "Ask whether the CNO effect differs between paired and unpaired animals.",
    "Check how much the Pairing result in mCherry+ depends on single animals.",
    "Describe the cell-type identity of each compartment."),
  model = c(
    paste0("Within-animal paired model on the ",
           sub(" only$", "", release_reader_text(supp_e$condition_scope, cmp)),
           " animals: ", model_text(supp_e$statistical_model),
           ". Comparisons: ", paste(supp_e_cmps, collapse = ", "), "."),
    paste0("The primary statistics ranked by ",
           release_reader_text(rank_sens$ranking_statistic, cmp),
           " instead of the moderated t; no new model. The released gene sets use the",
           " permutation settings of the primary GSEA (table 2m); the overlap below compares",
           " runs of both rankings at ", default_perm_text, "."),
    paste0(model_text(interaction$statistical_model), "; ", default_gsea, "."),
    paste0(model_text(loo$statistical_model), "; ", default_gsea, "."),
    paste0(release_reader_text(ewce_base$statistical_model, cmp), "; protein lists by ",
           release_reader_text(ewce_base$ranking_statistic, cmp), ".")),
  n_per_group = release_reader_text(c(supp_e$n_per_group, rank_sens$n_per_group,
                                      interaction$n_per_group, loo$n_per_group,
                                      ewce_base$n_per_group), cmp),
  results_location = c(
    paste("Key result only; Source Data", paste(panel_files("Supp. Fig. 4 E"), collapse = ", "),
          "holds the plotted gene sets"),
    paste("Data release: enrichment/GSEA_log2FC_sensitivity.tsv.gz (the gene sets; the file",
          "name keeps the analysis run's label); the overlap is given in this table only"),
    "This table (key result)",
    "This table (key result)",
    paste0("Table 2n (Baseline rows); Source Data ",
           paste(panel_files("Supp. Fig. 4 C"), collapse = ", "))),
  key_result = c(
    paste0("Proteins at FDR < 0.05: ", count_list("FDR proteins ([0-9/]+)"),
           ". GO-BP gene sets at FDR < 0.05: ", count_list("FDR GO-BP terms ([0-9/]+)"), "."),
    paste0("Overlap (Jaccard index) of the FDR < 0.05 gene sets of the two rankings, over the",
           " comparisons with any: ", jaccard, ". Computed on runs of both rankings at ",
           default_perm_text, ", not on the released tables, so an overlap computed from",
           " those can differ slightly."),
    release_reader_text(sentences_with(interaction$interpretation, "FDR < 0.05"), cmp),
    paste(release_reader_text(paste(release_split_sentences(loo$interpretation),
                                    collapse = " "), cmp), loo_primary_text),
    paste0(release_reader_text(ewce_key, cmp), " at the primary setting.")),
  stringsAsFactors = FALSE, check.names = FALSE)
secondary_analyses$key_result <- group_thousands(secondary_analyses$key_result)

# --------------------------------------------------------------------------------------
# 2p: figure panels
# --------------------------------------------------------------------------------------
# Source tables in this workbook, by letter; anything else is named by its Source Data file.
SOURCE_LETTER <- c(
  "differential_analysis/significant_proteins.tsv" = "2d",
  "differential_analysis/primary_differential_proteins.tsv.gz" = "2e",
  "enrichment/significant_gene_sets.tsv.gz" = "2f",
  "enrichment/primary_GSEA_GO_BP.tsv.gz" = "2g",
  "enrichment/primary_GSEA_GO_CC.tsv.gz" = "2h",
  "enrichment/primary_GSEA_KEGG.tsv.gz" = "2i",
  "enrichment/primary_ORA_GO_BP.tsv.gz" = "2j",
  "enrichment/primary_ORA_GO_CC.tsv.gz" = "2k",
  "enrichment/primary_EWCE.tsv.gz" = "2n")
# Sentences of the source notes that describe the analysis history, the permutation setting
# (stated once, in the README sheet) or internal folders; they are left out of table 2p.
NOTE_DROP <- c("nPermSimple", "03_output", "original", "superseded", "legacy", "earlier",
               "Replaces", "pooled", "must not", "Numerically indistinguishable",
               "Same equivalence", "clusterProfiler tree", "revision",
               "collection-plate-associated")
KEY_DROP <- c("original", "relabeled", "relabelled", "unchanged", "no longer", "retained", "legacy",
              "superseded", "pooled", "identical")

fm <- figure_map[nzchar(figure_map$source_table) | nzchar(figure_map$panel_source_data_release_path), ,
                 drop = FALSE]
if (!setequal(setdiff(figure_map$panel_label, fm$panel_label), c("Fig. 2 A", "Fig. 2 B"))) {
  stop("Only the artwork panels Fig. 2 A and B may lack a source.", call. = FALSE)
}
clean_note <- function(x) {
  s <- release_split_sentences(x)
  s <- s[!vapply(s, function(z) any(vapply(NOTE_DROP, grepl, logical(1), x = z, fixed = TRUE)),
                 logical(1))]
  release_reader_text(paste(s, collapse = " "), cmp)
}
clean_key_numbers <- function(x) {
  x <- gsub("(identical; ", "(", x, fixed = TRUE)
  x <- gsub("relabel(l)?ed Collection plate; explains", "collection plate explains", x)
  parts <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts) &
                   !vapply(parts, function(z) any(vapply(KEY_DROP, grepl, logical(1), x = z,
                                                         fixed = TRUE)), logical(1))]
  release_reader_text(paste(parts, collapse = "; "), cmp)
}
panel_source <- vapply(seq_len(nrow(fm)), function(i) {
  tables <- trimws(strsplit(fm$source_table[[i]], ";", fixed = TRUE)[[1]])
  letters_used <- unname(SOURCE_LETTER[tables[tables %in% names(SOURCE_LETTER)]])
  files <- basename(trimws(strsplit(fm$panel_source_data_release_path[[i]], ";",
                                    fixed = TRUE)[[1]]))
  files <- files[nzchar(files)]
  paste(c(if (length(letters_used)) paste("Table", paste(letters_used, collapse = ", ")),
          if (length(files)) paste("Source Data", paste(files, collapse = ", "))),
        collapse = "; ")
}, character(1))
# Reader titles of the panels. The figure source map keeps the titles of the panel records,
# which use the artwork's own words ("memory engram", "learning").
PANEL_TITLES <- c(
  "Fig. 2 C" = "PCA: PC1 vs PC2 and PC2 vs PC3 by compartment",
  "Fig. 2 D" = "Volcano: mCherry+ paired VEH vs unpaired VEH, effect size (SD units) against adjusted P",
  "Fig. 2 E" = "Rank abundance: mCherry+, paired VEH and unpaired VEH",
  "Fig. 2 F" = "GSEA bubble grid, GO-BP: Pairing (paired VEH vs unpaired VEH) in each compartment",
  "Fig. 2 G" = "GSEA lollipop, GO-BP and GO-CC: mCherry+ paired VEH vs unpaired VEH, rows by adjusted P",
  "Fig. 2 H" = "mCherry+: Pairing effect against the CNO effect in paired animals (both in SD units)",
  "Fig. 2 I" = "NES heatmap: Pairing and CNO in paired animals, gene sets significant in both",
  "Supp. Fig. 4 A" = "Protein and peptide counts per acquisition",
  "Supp. Fig. 4 B1" = "PCA by group",
  "Supp. Fig. 4 B2" = "PCA by hemisphere",
  "Supp. Fig. 4 B3" = "PCA by collection plate",
  "Supp. Fig. 4 C" = "Cell-type identity heatmap (EWCE baseline)",
  "Supp. Fig. 4 D" = "Rank abundance: NeuN+ and neuropil, unpaired VEH",
  "Supp. Fig. 4 E" = "Cross-compartment GSEA bubble plot",
  "Supp. Fig. 4 F" = "GSEA bubble grid, GO-BP: CNO in paired animals in each compartment")
if (!setequal(names(PANEL_TITLES), fm$panel_label)) {
  stop("PANEL_TITLES does not name exactly the panels of table 2p: ",
       paste(setdiff(union(names(PANEL_TITLES), fm$panel_label),
                     intersect(names(PANEL_TITLES), fm$panel_label)), collapse = ", "),
       call. = FALSE)
}
figure_panels <- data.frame(
  panel_label = fm$panel_label,
  panel_title = unname(PANEL_TITLES[fm$panel_label]),
  statistical_unit = release_reader_text(sub(" - .*$", "", fm$statistical_unit), cmp),
  panel_source = panel_source,
  rows_used = release_reader_text(fm$source_rows_filter, cmp),
  key_numbers = vapply(fm$key_numbers_revised, clean_key_numbers, character(1),
                       USE.NAMES = FALSE),
  note = vapply(fm$interpretation_note, clean_note, character(1), USE.NAMES = FALSE),
  stringsAsFactors = FALSE, check.names = FALSE)
# Supp. Fig. 4 E: the counts named by comparison, as in table 2o (the model is in the note).
e_row <- figure_panels$panel_label == "Supp. Fig. 4 E"
e_plotted <- regmatches(supp_e_key, regexec("([0-9]+) terms plotted", supp_e_key))[[1]]
if (sum(e_row) != 1L || length(e_plotted) != 2L) {
  stop("Supp. Fig. 4 E key numbers lack the plotted gene-set count.", call. = FALSE)
}
figure_panels$key_numbers[e_row] <- paste0(
  "Proteins at FDR < 0.05: ", count_list("FDR proteins ([0-9/]+)"),
  "; GO-BP gene sets at FDR < 0.05: ", count_list("FDR GO-BP terms ([0-9/]+)"), "; ",
  e_plotted[[2]], " gene sets plotted")
figure_panels$rows_used <- group_thousands(figure_panels$rows_used)
figure_panels$key_numbers <- group_thousands(figure_panels$key_numbers)

# --------------------------------------------------------------------------------------
# 2q: software behind the reported results
# --------------------------------------------------------------------------------------

# One "Used for" per component; the version is the one recorded by the run of that stage.
SOFTWARE <- list(
  list(stage = "canonical animal-level differential statistics", component = "R",
       used_for = "all analyses"),
  list(stage = "canonical animal-level differential statistics", component = "limma",
       used_for = "differential abundance"),
  list(stage = "differential enrichment (GSEA / ORA)", component = "clusterProfiler",
       used_for = "GSEA and ORA"),
  list(stage = "differential enrichment (GSEA / ORA)", component = "DOSE",
       used_for = "GSEA and ORA"),
  list(stage = "differential enrichment (GSEA / ORA)", component = "fgsea",
       used_for = "GSEA permutation P values"),
  list(stage = "differential enrichment (GSEA / ORA)", component = "GO.db",
       used_for = "GO gene sets (GSEA, ORA)"),
  list(stage = "differential enrichment (GSEA / ORA)", component = "org.Mm.eg.db",
       used_for = "GO gene sets (GSEA, ORA)"),
  list(stage = "differential enrichment (GSEA / ORA)", component = "KEGGREST",
       used_for = "KEGG gene sets (GSEA)"),
  list(stage = "EWCE cell-type enrichment", component = "EWCE",
       used_for = "cell-type enrichment (EWCE)"),
  list(stage = "EWCE cell-type enrichment", component = "ewceData::ctd()",
       used_for = "cell-type enrichment (EWCE)"))
COMPONENT_DISPLAY <- c("ewceData::ctd()" = "ewceData single-cell reference (CTD)")
software <- do.call(rbind, lapply(SOFTWARE, function(s) {
  v <- unique(sw$version[sw$component == s$component & sw$applies_to_stage == s$stage])
  if (length(v) != 1L) stop("No single recorded version of ", s$component, " for ", s$stage,
                            call. = FALSE)
  data.frame(component = s$component, version = v, used_for = s$used_for,
             stringsAsFactors = FALSE)
}))
# R must be the same version in every stage that recorded it.
r_versions <- unique(sw$version[sw$component == "R" &
                                  sw$applies_to_stage %in% vapply(SOFTWARE, `[[`, character(1),
                                                                  "stage")])
if (length(r_versions) != 1L) {
  stop("The analysis stages recorded different R versions: ", paste(r_versions, collapse = ", "),
       call. = FALSE)
}
kegg_row <- data.frame(component = "KEGG pathway gene sets (Mus musculus)",
                       version = paste("downloaded", RELEASE_STATISTICS_RUN$kegg_enrichment_date),
                       used_for = "KEGG gene sets (GSEA)", stringsAsFactors = FALSE)
at <- match("KEGGREST", software$component)
software <- rbind(software[seq_len(at), ], kegg_row, software[-seq_len(at), ])
software$component <- ifelse(software$component %in% names(COMPONENT_DISPLAY),
                             COMPONENT_DISPLAY[software$component], software$component)
rownames(software) <- NULL

# --------------------------------------------------------------------------------------
# README
# --------------------------------------------------------------------------------------

release_date <- format(Sys.time(), "%Y-%m-%d", tz = "UTC")

# Collection plate, as table 2a records it: plate per animal and injection block per plate.
plate_text <- paste(
  release_plate_layout_sentence(sample_metadata), release_plate_injection_sentence(),
  "In the Pairing comparisons, plate and pairing cannot be separated; no technical batch",
  "effect has been shown.")
hits_text <- if (!length(hit_cmps)) "No comparison has a protein at FDR < 0.05." else
  paste0(if (length(hit_cmps) == 1L) "Only " else "", labels_text(hit_cmps),
         if (length(hit_cmps) == 1L) " has" else " have", " proteins at FDR < 0.05 (",
         fmt_n(n_mapped_sig), " annotated proteins; table 2d).")
readme <- data.frame(
  item = c("Design and samples", "Unit of analysis", "Statistical model", "Effect size",
           "Comparisons", "Enrichment and gene sets without a P value", "Collection plate",
           "Complete tables and keys", "Version"),
  statement = c(
    paste0(inv$n_animals, " ", mv_wb("strain_or_breed"), " mice, ", inv$n_animals_per_stratum,
           " per group: paired VEH, paired CNO, unpaired VEH and unpaired CNO. From each",
           " animal, laser-capture microdissection collected four marker-defined compartments",
           " of the ", mv_wb("organism_part"), ": neuropil, c-Fos+, mCherry+ and NeuN+. The",
           " left and right hemisphere of each were measured separately by LC-MS (",
           mv_wb("labeling_strategy"), "): ", inv$n_measurement_records, " runs (table 2a)."),
    paste0("The left and right hemisphere runs of each animal and compartment were averaged",
           " before any test, giving ", inv$n_animal_level_units, " animal-level profiles",
           " (table 2b). The animal is the unit of analysis: n = ", inv$n_animals_per_stratum,
           " animals per group in every compartment."),
    paste0("Differential abundance was tested with limma ", stat_run$limma_version, " in R ",
           stat_run$r_version, ", with one linear model per compartment (~0 + condition over",
           " its ", stat_run$n_columns_in_model, " animal-level profiles; residual df ",
           stat_run$df_residual, ") and robust empirical Bayes moderation without an",
           " intensity trend. P values were adjusted with the Benjamini-Hochberg method",
           " within each comparison; an adjusted P value below 0.05 (FDR < 0.05) is called",
           " significant."),
    paste0("Abundances were standardized separately for each protein across all runs before",
           " analysis. The effect size is therefore a difference in standardized abundance",
           " between the two groups (group 1 minus group 2), in standard deviations of that",
           " protein. It is not a log2 fold change: 2 means two standard deviations, not a",
           " four-fold difference. A positive value means higher in group 1."),
    paste0(inv$n_primary_contrasts, " comparisons, three in each compartment: Pairing (paired",
           " VEH vs unpaired VEH), CNO in paired animals (paired CNO vs paired VEH) and CNO in",
           " unpaired animals (unpaired CNO vs unpaired VEH). Each is named '<compartment>:",
           " <group 1> vs <group 2>' (table 2c). ", hits_text, " Table 2e gives every protein",
           " in every comparison."),
    paste0("GSEA ranked the proteins of each comparison by the moderated t and tested GO",
           " biological process, GO cellular component and KEGG gene sets of ",
           param("GO_BP", "min_gs_size"), " to ", param("GO_BP", "max_gs_size"),
           " proteins (tables 2g to 2i). ORA tested the proteins at FDR < 0.05 for",
           " over-representation of GO terms against all measured annotated proteins",
           " (tables 2j and 2k). These tables list every gene set tested, not only the",
           " significant ones; table 2f lists the gene sets at FDR < 0.05 in the GSEA, with",
           " their leading-edge genes, and table 2m gives the settings. When the ranking of a",
           " comparison is strongly one-sided, fgsea cannot estimate a P value for some gene",
           " sets: these are absent from the GSEA tables rather than non-significant, and are",
           " listed in table 2l (GO) and counted in table 2m (GO and KEGG). GO GSEA used",
           " fgsea with nPermSimple = ", fmt_n(permutations$go_n_perm_simple), "; KEGG GSEA",
           " used the fgsea default",
           if (!is.na(permutations$default_n_perm_simple))
             paste0(" (", fmt_n(permutations$default_n_perm_simple), ")") else "", "."),
    plate_text,
    paste("The complete tables are files of the data release, described in its README_DATA.md",
          "and metadata/data_dictionary.tsv: all columns of every protein and comparison in",
          "differential_analysis/primary_differential_proteins.tsv.gz; the gene members of",
          "every GSEA and ORA row in enrichment/primary_GSEA_GO_BP.tsv.gz,",
          "primary_GSEA_GO_CC.tsv.gz, primary_GSEA_KEGG.tsv.gz, primary_ORA_GO_BP.tsv.gz and",
          "primary_ORA_GO_CC.tsv.gz, which also hold the ORA of the alternative list (proteins",
          "with |effect size| > 1 SD, no FDR filter); the effect-size-ranked GSEA, a",
          "sensitivity analysis, in enrichment/GSEA_log2FC_sensitivity.tsv.gz (the file name",
          "keeps the analysis run's label); and the status",
          "of every metadata field in metadata/metadata_field_provenance.tsv. Tables join on",
          "the comparison key (the last column of each comparison table,",
          "'<compartment>_<group 1>_over_<compartment>_<group 2>', with the compartments",
          "written mcherry, cfos, neuron and neuropil and the groups paired_veh, paired_cno,",
          "unpaired_veh and unpaired_cno), on Protein group (tables 2d and 2e) and on",
          "Animal-level profile (tables 2a and 2b)."),
    paste0("Release of ", release_date, ".")),
  stringsAsFactors = FALSE, check.names = FALSE)

# --------------------------------------------------------------------------------------
# assemble and present the sheets
# --------------------------------------------------------------------------------------

sheet_data <- list(
  README = readme,
  Samples = samples,
  Animal_Profiles = animal_profiles,
  Comparisons = comparisons,
  Significant_Proteins = sig_proteins,
  All_Proteins = all_proteins,
  Significant_Gene_Sets = sig_sets,
  GSEA_GO_BP = gsea_go,
  GSEA_GO_CC = gsea_cc,
  GSEA_KEGG = gsea_kegg,
  ORA_GO_BP = ora,
  ORA_GO_CC = ora_cc,
  GSEA_Untested_Sets = untested_sets,
  Enrichment_Settings = enrichment_settings,
  EWCE = ewce_tbl,
  Secondary_Analyses = secondary_analyses,
  Figure_Panels = figure_panels,
  Software = software)
captions <- vapply(names(RELEASE_SHEET_META), function(k) RELEASE_SHEET_META[[k]]$caption,
                   character(1))
captions[["Significant_Proteins"]] <- caption_2d
captions[["Significant_Gene_Sets"]] <- caption_2f
captions[["ORA_GO_BP"]] <- ora_caption("GO biological process")
captions[["ORA_GO_CC"]] <- ora_caption("GO cellular component")
captions[["EWCE"]] <- caption_2n
if (anyNA(captions)) {
  stop("Sheet(s) without a caption: ", paste(names(captions)[is.na(captions)], collapse = ", "),
       call. = FALSE)
}

present <- function(d, key) {
  if ("canonical_comparison" %in% names(d)) d <- add_comparison_label(d, key)
  d <- display_values(d)
  d <- release_select_columns(d, key)
  keep <- release_keep_columns(d, key)
  pruned <- release_drop_empty_and_duplicate(d, keep)
  d <- pruned$data
  const <- release_sheet_constants(d, keep = keep)
  if (length(const$constants)) d <- d[, setdiff(names(d), const$constants), drop = FALSE]
  rownames(d) <- NULL
  list(data = d, constants = const, dropped = pruned$dropped)
}
presented <- lapply(stats::setNames(names(sheet_data), names(sheet_data)),
                    function(k) present(sheet_data[[k]], k))

# Column definitions: one row per column shown on the lettered sheets.
lettered <- names(RELEASE_SHEET_META)[!is.na(vapply(RELEASE_SHEET_META, `[[`, character(1),
                                                    "letter"))]
undefined <- character(0)
column_definitions <- do.call(rbind, lapply(lettered, function(k) {
  cols <- names(presented[[k]]$data)
  do.call(rbind, lapply(cols, function(cn) {
    def <- release_reader_definition(cn, k)
    if (is.null(def)) {
      undefined <<- c(undefined, paste0(k, " :: ", cn))
      def <- c(definition = NA_character_, units = NA_character_)
    }
    data.frame(table = RELEASE_SHEET_META[[k]]$letter, column = release_display_labels(cn),
               definition = def[["definition"]], units = def[["units"]],
               stringsAsFactors = FALSE, check.names = FALSE)
  }))
}))
if (length(undefined)) {
  stop("Workbook column(s) with no entry in RELEASE_READER_DEFINITIONS:\n",
       paste("  -", undefined, collapse = "\n"), call. = FALSE)
}
presented$Column_Definitions <- present(column_definitions, "Column_Definitions")
SHEETS <- names(RELEASE_SHEET_META)
missing_sheets <- setdiff(SHEETS, names(presented))
if (length(missing_sheets)) {
  stop("Sheet(s) not built: ", paste(missing_sheets, collapse = ", "), call. = FALSE)
}
presented <- presented[SHEETS]

# --------------------------------------------------------------------------------------
# checks before writing
# --------------------------------------------------------------------------------------

tabs <- vapply(SHEETS, function(k) RELEASE_SHEET_META[[k]]$tab, character(1))
if (any(nchar(tabs) > 31L) || anyDuplicated(tolower(tabs)) || any(grepl("[][*?/\\\\:]", tabs))) {
  stop("Sheet names must be unique, at most 31 characters and free of []*?/\\:.", call. = FALSE)
}

# Hygiene on both namings: the machine names carry the joins and the
# gene_symbol/uniprot_accession check; the display names are what a reader sees.
machine_cols <- lapply(presented, function(p) names(p$data))
display_cols <- lapply(presented, function(p) release_display_labels(names(p$data)))
hygiene <- c(release_check_workbook_hygiene(machine_cols),
             release_check_workbook_hygiene(display_cols))
if (length(hygiene)) {
  stop("Workbook hygiene violations:\n", paste("  -", hygiene, collapse = "\n"), call. = FALSE)
}
unlabeled <- release_unlabelled_columns(unlist(machine_cols, use.names = FALSE))
if (length(unlabeled)) {
  stop("Column(s) with no display label in RELEASE_COLUMN_LABELS: ",
       paste(unlabeled, collapse = ", "), call. = FALSE)
}
if (any(grepl("_", unlist(display_cols), fixed = TRUE))) {
  stop("A display header contains an underscore.", call. = FALSE)
}
for (k in SHEETS) {
  d <- presented[[k]]$data
  if ("gene_symbol" %in% names(d) && release_column_is_misleading_gene_symbol(d$gene_symbol)) {
    stop("Sheet ", k, " has a gene_symbol column holding UniProt identifiers.", call. = FALSE)
  }
  if ("canonical_comparison" %in% names(d) &&
      (names(d)[[1]] != "comparison_label" || names(d)[[ncol(d)]] != "canonical_comparison")) {
    stop("Sheet ", k, " must show the comparison first and its key last.", call. = FALSE)
  }
}

EXCEL_MAX_ROWS <- 1048575L
EXCEL_MAX_CELL <- 32767L
for (k in SHEETS) {
  d <- presented[[k]]$data
  if (nrow(d) > EXCEL_MAX_ROWS) stop("Sheet ", k, " exceeds the Excel row limit.", call. = FALSE)
  for (col in names(d)) {
    if (is.character(d[[col]])) {
      longest <- suppressWarnings(max(nchar(d[[col]]), na.rm = TRUE))
      if (is.finite(longest) && longest > EXCEL_MAX_CELL) {
        stop("Sheet ", k, " column ", col, " has a ", longest,
             "-character cell, beyond the Excel limit.", call. = FALSE)
      }
    }
  }
}

constants_line <- function(k) {
  cst <- presented[[k]]$constants
  if (!length(cst$constants)) return("")
  values <- unname(cst$values)
  numeric_value <- grepl("^-?[0-9]+([.][0-9]+)?$", values)
  values[numeric_value] <- vapply(values[numeric_value], fmt_n, character(1),
                                  USE.NAMES = FALSE)
  paste0("Applies to every row: ",
         paste(paste0(release_display_labels(cst$constants), " = ", values),
               collapse = "; "), ".")
}

# Reader-facing text holds no internal term, history, path or hash.
reader_text <- c(captions, vapply(SHEETS, constants_line, character(1)),
                 readme$statement, column_definitions$definition,
                 unlist(secondary_analyses), unlist(figure_panels))
forbidden_hits <- RELEASE_READER_FORBIDDEN[vapply(RELEASE_READER_FORBIDDEN, function(w) {
  any(grepl(tolower(w), tolower(reader_text), fixed = TRUE))
}, logical(1))]
if (length(forbidden_hits)) {
  stop("Reader-facing workbook text contains: ", paste(forbidden_hits, collapse = ", "),
       call. = FALSE)
}
residual <- release_residual_machine_tokens(c(unlist(secondary_analyses), unlist(figure_panels),
                                              captions))
if (length(residual)) {
  stop("Machine names left in tables 2o, 2p or a caption: ", paste(residual, collapse = ", "),
       call. = FALSE)
}
INTERNAL_FIXED <- c("S:/", "D:\\", "Fabian", "03_output", "99_historical", "Figure3",
                    "Supplementary_proteomics")
INTERNAL_REGEX <- c("\\b[0-9a-f]{64}\\b", "\\b[A-Za-z0-9_]+\\.[Rr]\\b")
for (k in SHEETS) {
  d <- presented[[k]]$data
  for (col in names(d)[vapply(d, is.character, logical(1))]) {
    v <- unique(d[[col]])
    hit <- any(vapply(INTERNAL_FIXED, function(p) any(grepl(p, v, fixed = TRUE)), logical(1))) ||
      any(vapply(INTERNAL_REGEX, function(p) any(grepl(p, v, perl = TRUE)), logical(1)))
    if (hit) stop("Sheet ", k, " column ", col, " holds an internal path, name or hash.",
                  call. = FALSE)
  }
}

# --------------------------------------------------------------------------------------
# write the workbook
# --------------------------------------------------------------------------------------

INK <- "#1F3864"
style_header <- openxlsx::createStyle(
  textDecoration = "bold", fgFill = INK, fontColour = "#FFFFFF",
  halign = "left", valign = "bottom", wrapText = TRUE,
  border = "bottom", borderColour = INK, borderStyle = "medium")
style_title <- openxlsx::createStyle(textDecoration = "bold", fontSize = 13, fontColour = INK)
style_caption <- openxlsx::createStyle(fontSize = 10, fontColour = "#333333",
                                       wrapText = TRUE, valign = "top")
style_note <- openxlsx::createStyle(fontSize = 10, fontColour = "#555555",
                                    textDecoration = "italic", wrapText = TRUE, valign = "top")
style_text <- openxlsx::createStyle(valign = "top")
style_wrap <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
style_link <- openxlsx::createStyle(fontColour = "#0563C1", textDecoration = "underline",
                                    valign = "top")

MERGE_COLS <- 8L          # caption cells span the first 8 columns
NO_FILTER_SHEETS <- c("README", "Column_Definitions", "Software")   # prose or short lists
DEFAULT_WIDTH <- 8.43
SMALL_SHEET <- 300L       # sheets up to this many rows wrap long text and get row heights
HEADER_ROW <- 4L

text_height <- function(text, width) {
  lines <- ceiling(nchar(text) * 1.1 / max(width, 1))
  13 * pmax(1L, lines) + 3
}
merged_text <- function(wb, sheet, text, row, style, widths) {
  openxlsx::writeData(wb, sheet, text, startRow = row, startCol = 1)
  openxlsx::mergeCells(wb, sheet, cols = seq_len(MERGE_COLS), rows = row)
  openxlsx::addStyle(wb, sheet, style, rows = row, cols = seq_len(MERGE_COLS), gridExpand = TRUE)
  openxlsx::setRowHeights(wb, sheet, rows = row,
                          heights = text_height(text, sum(widths[seq_len(MERGE_COLS)])))
}
pad_widths <- function(w) c(w, rep(DEFAULT_WIDTH, max(0L, MERGE_COLS - length(w))))

wb <- openxlsx::createWorkbook()
openxlsx::modifyBaseFont(wb, fontSize = 10, fontName = "Calibri")

# ---- Contents -------------------------------------------------------------------------
contents <- data.frame(
  Table = vapply(SHEETS, function(k) {
    l <- RELEASE_SHEET_META[[k]]$letter
    if (is.na(l)) "" else l
  }, character(1)),
  Sheet = unname(tabs),
  Title = vapply(SHEETS, function(k) RELEASE_SHEET_META[[k]]$title, character(1)),
  Rows = vapply(presented, function(p) nrow(p$data), integer(1)),
  Columns = vapply(presented, function(p) ncol(p$data), integer(1)),
  Description = unname(captions[SHEETS]),
  check.names = FALSE, stringsAsFactors = FALSE)
rownames(contents) <- NULL
contents_widths <- c(7, 26, 34, 9, 9, 90)
contents_text <- c(
  paste0("Proteomics of ", inv$n_animals, " mice: ", inv$n_measurement_records,
         " LC-MS runs averaged into ", inv$n_animal_level_units, " animal-level profiles (",
         inv$n_sample_classes, " compartments; n = ", inv$n_animals_per_stratum,
         " animals per group in each), ", inv$n_primary_contrasts,
         " comparisons, one limma model per compartment (residual df ", stat_run$df_residual,
         "). The animal is the unit of analysis."),
  paste("Read the README sheet first. Effect sizes are differences in standardized abundance",
        "in SD units, not log2 fold changes. Select a sheet name to open it."))
openxlsx::addWorksheet(wb, "Contents", gridLines = FALSE, tabColour = INK)
openxlsx::setColWidths(wb, "Contents", cols = seq_along(contents_widths), widths = contents_widths)
openxlsx::writeData(wb, "Contents",
                    paste0("Supplementary Table ", RELEASE_SUPPLEMENTARY_TABLE, " | Contents"),
                    startRow = 1, startCol = 1)
openxlsx::addStyle(wb, "Contents", style_title, rows = 1, cols = 1)
merged_text(wb, "Contents", contents_text[[1]], 2L, style_caption, pad_widths(contents_widths))
merged_text(wb, "Contents", contents_text[[2]], 3L, style_note, pad_widths(contents_widths))
openxlsx::writeData(wb, "Contents", contents, startRow = HEADER_ROW, headerStyle = style_header)
openxlsx::setRowHeights(wb, "Contents", rows = HEADER_ROW, heights = 30)
# The sheet names stay plain text, so viewers that read stored values show them, and each
# carries a hyperlink element to its sheet. openxlsx writes internal links through its
# Hyperlink class (the one it uses to keep the links of a loaded workbook); the saved file
# is checked below.
Hyperlink <- get("Hyperlink", envir = asNamespace("openxlsx"))
contents_ws <- wb$worksheets[[match("Contents", names(wb))]]
contents_ws$hyperlinks <- lapply(seq_len(nrow(contents)), function(i) {
  Hyperlink$new(ref = paste0("B", HEADER_ROW + i), target = NULL,
                location = sprintf("'%s'!A1", contents$Sheet[[i]]),
                display = contents$Sheet[[i]], is_external = FALSE)
})
data_rows <- HEADER_ROW + seq_len(nrow(contents))
openxlsx::addStyle(wb, "Contents", style_text, rows = data_rows, cols = c(1, 4, 5),
                   gridExpand = TRUE)
openxlsx::addStyle(wb, "Contents", style_link, rows = data_rows, cols = 2, gridExpand = TRUE)
openxlsx::addStyle(wb, "Contents", style_wrap, rows = data_rows, cols = c(3, 6),
                   gridExpand = TRUE)
openxlsx::setRowHeights(wb, "Contents", rows = data_rows,
                        heights = text_height(contents$Description, contents_widths[[6]]))
openxlsx::freezePane(wb, "Contents", firstActiveRow = HEADER_ROW + 1L, firstActiveCol = 2)

# ---- data sheets ----------------------------------------------------------------------
for (k in SHEETS) {
  d <- presented[[k]]$data
  tab <- tabs[[k]]
  machine <- names(d)
  headers <- release_display_labels(machine)
  small <- nrow(d) <= SMALL_SHEET
  overrides <- RELEASE_WIDTH_OVERRIDES[[k]]

  widths <- numeric(ncol(d))
  wrap <- logical(ncol(d))
  for (j in seq_along(d)) {
    natural <- release_column_width(d[[j]], headers[[j]], max_width = 1000L)
    widths[[j]] <- if (!is.null(overrides) && machine[[j]] %in% names(overrides)) {
      overrides[[machine[[j]]]]
    } else min(natural, 52L)
    wrap[[j]] <- small && is.character(d[[j]]) && natural > widths[[j]]
  }

  openxlsx::addWorksheet(wb, tab, gridLines = FALSE)
  openxlsx::setColWidths(wb, tab, cols = seq_along(d), widths = widths)
  openxlsx::writeData(wb, tab, release_sheet_title(k), startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, tab, style_title, rows = 1, cols = 1)
  merged_text(wb, tab, captions[[k]], 2L, style_caption, pad_widths(widths))
  line3 <- constants_line(k)
  if (nzchar(line3)) merged_text(wb, tab, line3, 3L, style_note, pad_widths(widths))

  shown <- d
  names(shown) <- headers
  openxlsx::writeData(wb, tab, shown, startRow = HEADER_ROW, headerStyle = style_header,
                      withFilter = nrow(d) > 1L && !k %in% NO_FILTER_SHEETS)
  openxlsx::setRowHeights(wb, tab, rows = HEADER_ROW, heights = 30)
  openxlsx::freezePane(wb, tab, firstActiveRow = HEADER_ROW + 1L, firstActiveCol = 2)

  if (nrow(d)) {
    rows <- HEADER_ROW + seq_len(nrow(d))
    for (j in seq_along(d)) {
      fmt <- release_numfmt_for(d[[j]], machine[[j]], k)
      style <- if (!is.na(fmt)) openxlsx::createStyle(numFmt = fmt, valign = "top") else
        if (wrap[[j]]) style_wrap else style_text
      openxlsx::addStyle(wb, tab, style, rows = rows, cols = j, gridExpand = TRUE)
    }
    if (small && any(wrap)) {
      lines <- do.call(pmax, c(lapply(which(wrap), function(j) {
        ceiling(nchar(ifelse(is.na(d[[j]]), "", d[[j]])) * 1.1 / widths[[j]])
      }), list(1)))
      openxlsx::setRowHeights(wb, tab, rows = rows, heights = 13 * lines + 3)
    }
  }

  release_log("    ", format(tab, width = 26), format(nrow(d), width = 7), " rows x ",
              format(ncol(d), width = 2), " cols",
              if (length(presented[[k]]$constants$constants))
                paste0("  (", length(presented[[k]]$constants$constants),
                       " shared value(s) above the header)") else "",
              if (length(presented[[k]]$dropped))
                paste0("  (blank or duplicate: ", paste(presented[[k]]$dropped,
                                                        collapse = ", "), ")") else "")
}

xlsx_path <- release_path("editor_source_data", basename(XLSX_REL))
openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
saved <- release_read_workbook_xml(xlsx_path)
if (!identical(saved$sheets[["Contents"]]$hyperlinks, sprintf("'%s'!A1", contents$Sheet)) ||
    !identical(unname(saved$sheets[["Contents"]]$rows[[HEADER_ROW]][["B"]]), "Sheet")) {
  stop("The Contents sheet of the saved workbook does not link each sheet name to its sheet.",
       call. = FALSE)
}
size_mb <- round(file.info(xlsx_path)$size / 1024 / 1024, 1)
release_log("  wrote ", basename(XLSX_REL), " (", length(SHEETS) + 1L,
            " sheets incl. Contents, ", size_mb, " MB)")

release_register(XLSX_REL,
                 "Supplementary Table 2: the proteomics workbook, one table per sheet",
                 c(rel("metadata", "sample_metadata.tsv"),
                   rel("metadata", "animal_level_sample_metadata.tsv"),
                   rel("metadata", "primary_contrast_manifest.tsv"),
                   rel("metadata", "secondary_analysis_manifest.tsv"),
                   rel("differential_analysis", "primary_differential_proteins.tsv.gz"),
                   rel("differential_analysis", "primary_differential_summary.tsv"),
                   rel("differential_analysis", "significant_proteins.tsv"),
                   rel("enrichment", "primary_GSEA_GO_BP.tsv.gz"),
                   rel("enrichment", "primary_GSEA_GO_CC.tsv.gz"),
                   rel("enrichment", "primary_GSEA_KEGG.tsv.gz"),
                   rel("enrichment", "primary_ORA_GO_BP.tsv.gz"),
                   rel("enrichment", "primary_ORA_GO_CC.tsv.gz"),
                   rel("enrichment", "significant_gene_sets.tsv.gz"),
                   rel("enrichment", "gsea_untested_gene_sets.tsv"),
                   rel("enrichment", "enrichment_coverage.tsv"),
                   rel("enrichment", "enrichment_run_parameters.tsv"),
                   rel("enrichment", "primary_EWCE.tsv.gz"),
                   rel("provenance", "figure_source_map.tsv"),
                   pkg_versions_path, ewce_session_path, stat_run$path, ewce_model$path),
                 c(rep(NA_character_, 18), release_sha256(pkg_versions_path),
                   release_sha256(ewce_session_path), release_sha256(stat_run$path),
                   release_sha256(ewce_model$path)),
                 STAGE, "xlsx")

# The sheet inventory, so a validator can check the workbook without opening Excel. Both
# namings are recorded: `columns` is the machine naming the checks key off, and
# `display_columns` is what the reader sees.
inventory <- data.frame(
  sheet = unname(tabs),
  sheet_key = SHEETS,
  supplementary_table = unname(contents$Table),
  title = unname(contents$Title),
  n_rows = vapply(presented, function(p) nrow(p$data), integer(1)),
  n_columns = vapply(presented, function(p) ncol(p$data), integer(1)),
  columns = vapply(presented, function(p) paste(names(p$data), collapse = ";"), character(1)),
  display_columns = vapply(presented, function(p) {
    paste(release_display_labels(names(p$data)), collapse = ";")
  }, character(1)),
  table_level_constants = vapply(presented, function(p) {
    if (!length(p$constants$constants)) return("")
    paste(paste0(p$constants$constants, "=", p$constants$values), collapse = ";")
  }, character(1)),
  stringsAsFactors = FALSE, check.names = FALSE)
rownames(inventory) <- NULL
w <- release_write_table(inventory,
                         release_path("editor_source_data", "workbook_sheet_inventory.tsv"))
release_register("editor_source_data/workbook_sheet_inventory.tsv",
                 "sheet-by-sheet inventory of Supplementary Table 2: machine and display columns",
                 XLSX_REL, NA_character_, STAGE, "tsv")
release_log("  wrote workbook_sheet_inventory.tsv (", w$rows, "x", w$cols, ")")
release_log("stage 07 complete")
