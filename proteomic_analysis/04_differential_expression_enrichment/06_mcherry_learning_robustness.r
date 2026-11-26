#!/usr/bin/env Rscript

# Robustness checks under the per-compartment (primary) model, and the report for this stage and
# 05_pairing_cno_interaction.r.
#
# Model: as 01_preprocessing/02b_limma_stat_gct.r (within_compartment mode): one limma fit per
# sample class, ~0 + condition over that class's 12 animals (n = 3 per condition, residual df 8),
# eBayes(robust = TRUE, trend = FALSE), Benjamini-Hochberg within the contrast. Effects are
# standardized abundance differences (SD units).
#
#   1. mCherry learning (paired_veh - unpaired_veh), leave one animal out: each of the six vehicle
#      animals is dropped in turn and the mCherry model refitted (2 vs 3 animals, residual df 7).
#      Records FDR and q < 0.10 counts, correlation of effects and t with the full fit, and GO-BP /
#      GO-CC GSEA (settings and seeds of 01_clusterProfiler.r for this comparison) with the number
#      of the full fit's significant terms kept.
#   2. t-ranked vs effect-size-ranked GSEA: FDR < 0.05 overlap (Jaccard) per contrast for GO-BP and
#      GO-CC, read from the per-compartment enrichment run (no new GSEA).
#   3. Paired-CNO vs unpaired-CNO: GO terms FDR-significant in both arms and how many have
#      opposite sign, from the same enrichment run.
#   4. interaction_and_robustness_report.txt: plain-text summary of these results and of
#      05_pairing_cno_interaction.r, which must have been run into the same output root.
#
# Writes only into <output root>/robustness/ (refused when not empty) and the report file
# (refused when present) unless PROTEOMICS_INFERENTIAL_FORCE=true.

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_file)) sub("^--file=", "", script_file[[1]]) else "04_differential_expression_enrichment/06_mcherry_learning_robustness.r"
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

source(file.path(project_root, "R", "analysis_labels.R"))
source(file.path(project_root, "R", "project_path_utils.R"))
source(file.path(project_root, "R", "animal_level_enrichment_utils.R"))
source(file.path(project_root, "R", "inferential_checks_utils.R"))

option_or_env <- function(option_name, env_name, default) {
  value <- getOption(option_name)
  if (!is.null(value) && nzchar(trimws(as.character(value)))) return(as.character(value))
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(trimws(value))) return(value)
  default
}

animal_root <- file.path(PROTEOMICS_DATA_ROOT, "02_data", "animal_level")
input_gct <- project_normalize_path(option_or_env(
  "proteomics.inferential_input_gct", "PROTEOMICS_INFERENTIAL_INPUT_GCT",
  file.path(animal_root, "input_gct", "neha_protigy_input_animal_level_primary.gct")
), must_work = TRUE)
mapped_root <- project_normalize_path(option_or_env(
  "proteomics.inferential_mapped_root", "PROTEOMICS_INFERENTIAL_MAPPED_ROOT",
  file.path(animal_root, "within_compartment_20261124", "mapped")
), must_work = TRUE)
mapped_index_path <- file.path(mapped_root, "indexMappedComparisons.csv")
enrichment_root <- project_normalize_path(option_or_env(
  "proteomics.inferential_enrichment_root", "PROTEOMICS_INFERENTIAL_ENRICHMENT_ROOT",
  file.path(PROTEOMICS_DATA_ROOT, "03_output", "enrichment", "enrichment_within_compartment_20261124")
), must_work = TRUE)
output_root <- project_normalize_path(option_or_env(
  "proteomics.inferential_output_root", "PROTEOMICS_INFERENTIAL_OUTPUT_ROOT",
  file.path(PROTEOMICS_DATA_ROOT, "03_output", "inferential_checks", "within_compartment_20261124")
))
force <- inferential_flag(option_or_env("proteomics.inferential_force", "PROTEOMICS_INFERENTIAL_FORCE", "false"), "PROTEOMICS_INFERENTIAL_FORCE")
run_gsea <- inferential_flag(option_or_env("proteomics.inferential_gsea", "PROTEOMICS_INFERENTIAL_GSEA", "true"), "PROTEOMICS_INFERENTIAL_GSEA")
workers <- inferential_workers(option_or_env("proteomics.inferential_gsea_workers", "PROTEOMICS_INFERENTIAL_GSEA_WORKERS", "3"), "PROTEOMICS_INFERENTIAL_GSEA_WORKERS")

required_packages <- c("limma", "digest", if (run_gsea) c("clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "DOSE", "fgsea", "BiocParallel", "withr"))
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages before running this stage: ", paste(missing_packages, collapse = ", "), call. = FALSE)

invisible(validate_inferential_output_root(output_root))
if (!dir.exists(output_root)) stop("Run 05_pairing_cno_interaction.r into this output root first: ", output_root, call. = FALSE)
output_root <- project_normalize_path(output_root, must_work = TRUE)
robustness_root <- file.path(output_root, "robustness")
report_path <- file.path(output_root, "interaction_and_robustness_report.txt")
if ((!inferential_dir_is_empty(robustness_root) || file.exists(report_path)) && !force) {
  stop("robustness/ is not empty or the report exists; set PROTEOMICS_INFERENTIAL_FORCE=true to write into them: ", output_root, call. = FALSE)
}
interaction_files <- file.path(output_root, c(
  "pairing_cno_interaction_protein_summary.csv", "pairing_cno_interaction_top_proteins.csv",
  "pairing_cno_interaction_diagnostics.csv", "pairing_cno_interaction_plate_composition.csv",
  "pairing_cno_interaction_gsea_counts.csv", "pairing_cno_interaction_top_terms.csv",
  "pairing_cno_interaction_run_manifest.csv"
))
names(interaction_files) <- sub("^pairing_cno_interaction_|[.]csv$", "", basename(interaction_files))
names(interaction_files) <- sub("[.]csv$", "", names(interaction_files))
if (!all(file.exists(interaction_files))) {
  stop("Run 05_pairing_cno_interaction.r (with GSEA) into this output root first; missing: ",
       paste(basename(interaction_files[!file.exists(interaction_files)]), collapse = ", "), call. = FALSE)
}
dir.create(robustness_root, recursive = TRUE, showWarnings = FALSE)
config <- resolve_enrichment_config(mapped_root = mapped_root, mapped_index = mapped_index_path, output_root = output_root)
gsea_settings <- config[c("min_gs_size", "max_gs_size", "pvalue_cutoff", "p_adjust_method", "fdr_threshold")]
fdr <- config$fdr_threshold
timestamp_utc <- format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")
inputs_read <- list()
note_input <- function(role, path) {
  inputs_read[[length(inputs_read) + 1L]] <<- manifest_rows("input", role, path, project_file_sha256(path))
  path
}
read_input_csv <- function(role, path) utils::read.csv(note_input(role, path), stringsAsFactors = FALSE, check.names = FALSE)

# ---- 1. mCherry learning, leave one animal out --------------------------------------------
learning_comparison <- "mcherry_paired_veh_over_mcherry_unpaired_veh"
input <- read_animal_level_input_gct(note_input("input_gct", input_gct))
mapped_index <- read_enrichment_mapped_index(mapped_index_path, mapped_root, verify_hashes = FALSE)
invisible(note_input("mapped_index", mapped_index_path))
learning_row <- mapped_index[mapped_index$canonical_comparison == learning_comparison, , drop = FALSE]
learning_mapped <- read_enrichment_mapped_file(
  note_input(paste0("mapped:", learning_comparison), learning_row$mapped_input_path),
  expected_rows = learning_row$n_output_mapped_rows, expected_sha256 = learning_row$mapped_output_sha256
)
mch <- input$annotation[input$annotation$sample_class == "mcherry", , drop = FALSE]
fit_learning <- function(keep) fit_condition_contrast(input$values[, mch$column[keep], drop = FALSE], mch$condition[keep], MCHERRY_LEARNING_CONTRAST)
full <- fit_learning(rep(TRUE, nrow(mch)))
learning_t_diff <- max(abs(full$table[learning_mapped$original_protein_id, "t"] - learning_mapped$t))
if (!is.finite(learning_t_diff) || learning_t_diff > 1e-8) {
  stop("The per-compartment fit does not reproduce the mapped primary t for mCherry learning (max diff ", signif(learning_t_diff, 3), ").", call. = FALSE)
}
drop_animals <- mch$animal[mch$condition %in% c("paired_veh", "unpaired_veh")]
loo_fits <- c(list(none = full), stats::setNames(lapply(drop_animals, function(an) fit_learning(mch$animal != an)), drop_animals))
full_tab <- full$table
mapped_ids <- learning_mapped$original_protein_id
loo <- do.call(rbind, lapply(names(loo_fits), function(k) {
  f <- loo_fits[[k]]
  tab <- f$table
  m <- tab[mapped_ids, , drop = FALSE]
  sig <- tab$adj.P.Val < 0.05
  data.frame(
    animal_dropped = k, condition_dropped = if (k == "none") "" else mch$condition[mch$animal == k],
    df_residual = f$df_residual, df_prior_median = f$df_prior_median, s2_prior = f$s2_prior,
    n_fdr05 = sum(sig), n_fdr05_mapped = sum(m$adj.P.Val < 0.05),
    n_up = sum(sig & tab$logFC > 0), n_down = sum(sig & tab$logFC < 0),
    n_retained_of_full = sum(sig & full_tab$adj.P.Val < 0.05),
    n_fdr10 = sum(tab$adj.P.Val < 0.10), n_fdr10_mapped = sum(m$adj.P.Val < 0.10),
    n_nominal_p05_mapped = sum(m$P.Value < 0.05), min_adj_p_mapped = min(m$adj.P.Val),
    r_effect = stats::cor(tab$logFC, full_tab$logFC), r_t = stats::cor(tab$t, full_tab$t),
    median_t_mapped = stats::median(m$t), pct_t_positive_mapped = 100 * mean(m$t > 0),
    stringsAsFactors = FALSE
  )
}))
rownames(loo) <- NULL

loo_gsea <- NULL
full_gsea_check <- NULL
gsea_job_rows <- NULL
if (run_gsea) {
  loo_gsea_root <- file.path(robustness_root, "mcherry_learning_leave_one_out_gsea")
  jobs <- list()
  for (k in names(loo_fits)) {
    collapsed <- collapse_enrichment_accessions(inferential_mapped_stats(learning_mapped, loo_fits[[k]]$table), "log2fc")$collapsed
    unit <- if (k == "none") "full" else paste0("drop_", k)
    for (ont in c("BP", "CC")) {
      jobs[[length(jobs) + 1L]] <- inferential_gsea_job_spec(
        paste(unit, ont, sep = "|"), unit, ont, collapsed, learning_comparison, config$gsea_seed_base, gsea_settings,
        output = file.path(loo_gsea_root, paste0(unit, "_GSEA_GO_", ont, ".csv")), return_table = TRUE
      )
    }
  }
  message("Running ", length(jobs), " GSEA jobs in ", min(workers, length(jobs)), " R process(es)...")
  results <- run_inferential_gsea_jobs(jobs, workers, project_root)
  gsea_job_rows <- do.call(rbind, lapply(results, `[[`, "summary"))
  full_sets <- lapply(c(BP = "BP", CC = "CC"), function(ont) significant_terms(results[[paste("full", ont, sep = "|")]]$table, fdr))
  loo_gsea <- do.call(rbind, lapply(names(loo_fits), function(k) {
    unit <- if (k == "none") "full" else paste0("drop_", k)
    row <- data.frame(animal_dropped = k, stringsAsFactors = FALSE)
    for (ont in c("BP", "CC")) {
      tab <- results[[paste(unit, ont, sep = "|")]]$table
      s <- significant_terms(tab, fdr)
      key <- paste0("go_", tolower(ont))
      row[[paste0(key, "_fdr05")]] <- length(s)
      row[[paste0(key, "_up")]] <- sum(tab$NES[tab$ID %in% s] > 0)
      row[[paste0(key, "_down")]] <- sum(tab$NES[tab$ID %in% s] < 0)
      row[[paste0(key, "_retained_of_full")]] <- length(intersect(s, full_sets[[ont]]))
      row[[paste0(key, "_full_fdr05")]] <- length(full_sets[[ont]])
    }
    row
  }))
  loo <- merge(loo, loo_gsea, by = "animal_dropped", sort = FALSE)
  loo <- loo[match(names(loo_fits), loo$animal_dropped), , drop = FALSE]

  # The full-fit GSEA must reproduce the stored per-compartment enrichment for this comparison.
  full_gsea_check <- do.call(rbind, lapply(c("BP", "CC"), function(ont) {
    stored <- read_input_csv(paste0("enrichment_GO_", ont, ":", learning_comparison),
                             file.path(enrichment_root, paste0("GO_", ont), "per_comparison", learning_comparison, paste0("GSEA_GO_", ont, ".csv")))
    ours <- results[[paste("full", ont, sep = "|")]]$table
    m <- merge(stored[, c("ID", "NES", "pvalue", "core_enrichment")], ours[, c("ID", "NES", "pvalue", "core_enrichment")], by = "ID", suffixes = c("_stored", "_rerun"))
    data.frame(
      check = paste0("mCherry learning full-fit GO-", ont, " GSEA vs stored enrichment"),
      stored_terms = nrow(stored), rerun_terms = nrow(ours), shared_terms = nrow(m),
      stored_fdr05 = length(significant_terms(stored, fdr)), rerun_fdr05 = length(significant_terms(ours, fdr)),
      same_significant_set = setequal(significant_terms(stored, fdr), significant_terms(ours, fdr)),
      max_abs_diff_NES = max(abs(m$NES_stored - m$NES_rerun)), max_abs_diff_pvalue = max(abs(m$pvalue_stored - m$pvalue_rerun)),
      identical_core_enrichment = identical(m$core_enrichment_stored, m$core_enrichment_rerun),
      stringsAsFactors = FALSE
    )
  }))
}
invisible(write_enrichment_csv(loo, file.path(robustness_root, "mcherry_learning_leave_one_out.csv")))

# ---- 2. and 3. from the per-compartment enrichment run ------------------------------------
read_gsea <- function(role, path) {
  path <- project_normalize_path(path, must_work = TRUE)
  if (!project_path_is_within(path, enrichment_root)) stop("Enrichment table outside the configured enrichment root: ", path, call. = FALSE)
  read_input_csv(role, path)
}
indexes <- lapply(c(BP = "BP", CC = "CC"), function(ont) {
  index <- read_input_csv(paste0("enrichment_index_GO_", ont), file.path(enrichment_root, paste0("GO_", ont), "indexEnrichmentComparisons.csv"))
  validate_enrichment_index(index, require_files = TRUE)
  if (any(index$execution_status != "success")) stop("GO-", ont, " enrichment index has failed comparisons.", call. = FALSE)
  index
})
t_tables <- list()
rank_robustness <- do.call(rbind, lapply(c("BP", "CC"), function(ont) {
  index <- indexes[[ont]]
  do.call(rbind, lapply(seq_len(nrow(index)), function(i) {
    cmp <- index$canonical_comparison[[i]]
    t_tab <- read_gsea(paste0("GSEA_GO_", ont, ":", cmp), index$go_gsea_output[[i]])
    l_tab <- read_gsea(paste0("GSEA_GO_", ont, "_log2fc_sensitivity:", cmp), index$sensitivity_go_gsea_output[[i]])
    t_tables[[paste(ont, cmp)]] <<- t_tab
    warnings <- utils::read.csv(index$warnings_path[[i]], stringsAsFactors = FALSE, check.names = FALSE)
    m <- merge(t_tab[, c("ID", "NES")], l_tab[, c("ID", "NES")], by = "ID", suffixes = c("_t", "_l"))
    st <- significant_terms(t_tab, fdr)
    sl <- significant_terms(l_tab, fdr)
    top <- function(x, n) utils::head(x$ID[order(x$p.adjust, -abs(x$NES))], n)
    data.frame(
      ontology = ont, comparison = cmp, seed = index$go_gsea_seed[[i]], n_uniprot = index$rank_vector_size[[i]],
      n_terms_log2fc = nrow(l_tab), n_terms_t = nrow(t_tab), n_common_terms = nrow(m),
      NES_pearson = stats::cor(m$NES_t, m$NES_l), NES_spearman = stats::cor(m$NES_t, m$NES_l, method = "spearman"),
      NES_sign_concordance = mean(sign(m$NES_t) == sign(m$NES_l)),
      fdr05_log2fc = length(sl), fdr05_t = length(st), fdr05_both = length(intersect(st, sl)),
      fdr05_log2fc_only = length(setdiff(sl, st)), fdr05_t_only = length(setdiff(st, sl)),
      fdr05_jaccard = if (length(union(st, sl))) length(intersect(st, sl)) / length(union(st, sl)) else NA_real_,
      top20_overlap = length(intersect(top(t_tab, 20), top(l_tab, 20))),
      top50_overlap = length(intersect(top(t_tab, 50), top(l_tab, 50))),
      log2fc_warning_count = sum(warnings$analysis == paste0("GO_", ont, "_GSEA_log2fc_sensitivity")),
      t_warning_count = sum(warnings$analysis == paste0("GO_", ont, "_GSEA")),
      stringsAsFactors = FALSE
    )
  }))
}))
invisible(write_enrichment_csv(rank_robustness, file.path(robustness_root, "gsea_rank_robustness_summary.csv")))

cno_arms <- do.call(rbind, lapply(sample_classes, function(cls) do.call(rbind, lapply(c("BP", "CC"), function(ont) {
  a <- t_tables[[paste(ont, sprintf("%s_paired_cno_over_%s_paired_veh", cls, cls))]]
  b <- t_tables[[paste(ont, sprintf("%s_unpaired_cno_over_%s_unpaired_veh", cls, cls))]]
  a <- a[a$ID %in% significant_terms(a, fdr), c("ID", "NES")]
  b <- b[b$ID %in% significant_terms(b, fdr), c("ID", "NES")]
  m <- merge(a, b, by = "ID", suffixes = c("_paired", "_unpaired"))
  data.frame(sample_class = cls, ontology = ont, paired_cno_fdr05 = nrow(a), unpaired_cno_fdr05 = nrow(b),
             shared = nrow(m), opposite_sign = sum(sign(m$NES_paired) != sign(m$NES_unpaired)),
             same_sign = sum(sign(m$NES_paired) == sign(m$NES_unpaired)), stringsAsFactors = FALSE)
}))))
invisible(write_enrichment_csv(cno_arms, file.path(robustness_root, "cno_arm_shared_terms.csv")))

checks <- rbind(
  data.frame(check = "mCherry learning full-fit t vs mapped primary table", value = learning_t_diff, passed = learning_t_diff <= 1e-8, stringsAsFactors = FALSE),
  if (!is.null(full_gsea_check)) data.frame(
    check = paste0(full_gsea_check$check, ": same FDR set (", full_gsea_check$rerun_fdr05, " = ", full_gsea_check$stored_fdr05, "), max |dNES|"),
    value = full_gsea_check$max_abs_diff_NES,
    passed = full_gsea_check$same_significant_set & full_gsea_check$max_abs_diff_NES < 1e-10 & full_gsea_check$identical_core_enrichment,
    stringsAsFactors = FALSE
  )
)
invisible(write_enrichment_csv(checks, file.path(robustness_root, "checks.csv")))

# ---- 4. report ----------------------------------------------------------------------------
x <- lapply(interaction_files, function(p) utils::read.csv(note_input(paste0("interaction:", basename(p)), p), stringsAsFactors = FALSE, check.names = FALSE))
label_of <- c(mcherry = "mCherry+", neuropil = "neuropil", cfos = "c-Fos+", neuron = "NeuN+")
f2 <- function(v) sprintf("%.2f", v)
fp <- function(v) trimws(formatC(v, digits = 2, format = "g"))
signed <- function(v, d = 2) sprintf(paste0("%+.", d, "f"), v)
text_table <- function(df) {
  cells <- lapply(df, as.character)
  widths <- pmax(nchar(names(df)), vapply(cells, function(v) max(nchar(v), 0L), integer(1)))
  line <- function(v) sub(" +$", "", paste(sprintf("%-*s", widths, v), collapse = "  "))
  c(line(names(df)), line(vapply(widths, function(w) strrep("-", w), character(1))),
    vapply(seq_len(nrow(df)), function(i) line(vapply(cells, `[[`, character(1), i)), character(1)))
}
ps <- x$protein_summary
dg <- x$diagnostics
gc <- x$gsea_counts
interaction_manifest <- x$run_manifest
interaction_run <- interaction_manifest$value[interaction_manifest$section == "run" & interaction_manifest$name == "run_timestamp_utc"]
seed_bases <- interaction_manifest$value[interaction_manifest$section == "parameter" & interaction_manifest$name == "seed_range_bases"]

protein_rows <- data.frame(
  compartment = label_of[ps$sample_class], tested = ps$n_tested,
  `nominal P<0.05` = sprintf("%d (%.1f%%)", ps$n_nominal_p05, ps$pct_nominal_p05),
  `FDR<0.05` = ps$n_fdr05, `FDR<0.10` = dg$n_fdr10[match(ps$sample_class, dg$sample_class)],
  `min q` = fp(dg$min_adj_p[match(ps$sample_class, dg$sample_class)]),
  `df.prior` = f2(dg$df_prior_median[match(ps$sample_class, dg$sample_class)]),
  `median effect (SD)` = signed(ps$median_interaction_log2fc),
  check.names = FALSE, stringsAsFactors = FALSE
)
protein_tables <- lapply(stats::setNames(sample_classes, sample_classes), function(cls) {
  utils::read.csv(note_input(paste0("interaction:", cls, "/interaction_protein_table.csv"),
                             file.path(output_root, "posthoc_pairing_x_cno_interaction", cls, "interaction_protein_table.csv")),
                  stringsAsFactors = FALSE)
})
named <- list(mcherry = "Exoc3", cfos = c("Chmp2a", "Cntfr"))
named_lines <- unlist(lapply(names(named), function(cls) {
  p <- protein_tables[[cls]]
  p$rank_abs_t <- rank(-abs(p$t), ties.method = "first")
  vapply(named[[cls]], function(g) {
    r <- p[p$Description == g, , drop = FALSE][1, ]
    sprintf("  %s (%s): effect %s SD, t %s, P %s, q %s, rank %d of %d by |t|", g, label_of[[cls]], signed(r$logFC), f2(r$t), fp(r$P.Value), fp(r$adj.P.Val), r$rank_abs_t, nrow(p))
  }, character(1))
}))
top_lines <- vapply(sample_classes, function(cls) {
  p <- protein_tables[[cls]]
  p <- utils::head(p[order(-abs(p$t)), , drop = FALSE], 5L)
  sprintf("  %s: %s", label_of[[cls]], paste(sprintf("%s %s (t %s, q %s)", p$Description, signed(p$logFC), f2(p$t), fp(p$adj.P.Val)), collapse = "; "))
}, character(1))

count_cell <- function(cls, ont) {
  p <- gc[gc$sample_class == cls & gc$ontology == ont & gc$primary_seed, , drop = FALSE]
  all_seeds <- gc$n_fdr05[gc$sample_class == cls & gc$ontology == ont]
  range <- if (length(all_seeds) > 1L) sprintf(" [%d-%d]", min(all_seeds), max(all_seeds)) else ""
  sprintf("%d (%d/%d)%s", p$n_fdr05, p$n_up, p$n_down, range)
}
gsea_rows <- data.frame(
  compartment = label_of[sample_classes],
  `GO-BP` = vapply(sample_classes, count_cell, character(1), ont = "BP"),
  `GO-CC` = vapply(sample_classes, count_cell, character(1), ont = "CC"),
  KEGG = vapply(sample_classes, count_cell, character(1), ont = "KEGG"),
  check.names = FALSE, stringsAsFactors = FALSE
)
shift_rows <- data.frame(
  compartment = label_of[dg$sample_class],
  `% GO-BP NES>0` = sprintf("%.0f", dg$go_bp_pct_NES_positive), `% GO-CC NES>0` = sprintf("%.0f", dg$go_cc_pct_NES_positive),
  `% KEGG NES>0` = sprintf("%.0f", dg$kegg_pct_NES_positive),
  `protein median t` = signed(dg$protein_median_t), `% t>0` = sprintf("%.0f", dg$protein_pct_t_positive),
  `median effect (SD)` = signed(dg$protein_median_effect),
  check.names = FALSE, stringsAsFactors = FALSE
)
tt <- x$top_terms
term_lines <- unlist(lapply(sample_classes, function(cls) {
  vapply(c("positive", "negative"), function(direction) {
    y <- tt[tt$sample_class == cls & tt$direction == direction & tt$rank <= 3, , drop = FALSE]
    y <- y[order(match(y$ontology, c("BP", "CC", "KEGG")), y$rank), , drop = FALSE]
    body <- if (nrow(y)) paste(sprintf("%s %s (NES %s, q %s)", c(BP = "BP", CC = "CC", KEGG = "KEGG")[y$ontology], y$Description, f2(y$NES), fp(y$p.adjust)), collapse = "; ") else "none"
    sprintf("  %s, %s: %s", label_of[[cls]], direction, body)
  }, character(1))
}))

pc <- x$plate_composition
pc1 <- pc[pc$sample_class == sample_classes[[1]], , drop = FALSE]
plate_weight <- unique(round(dg$plate2_minus_plate1_weight, 12))
mixed <- pc1[pc1$plate2_fraction > 0 & pc1$plate2_fraction < 1, , drop = FALSE]
minority <- vapply(seq_len(nrow(mixed)), function(i) {
  if (mixed$plate2_fraction[[i]] < 0.5) sprintf("%s on Plate2 in %s", mixed$plate2_animals[[i]], mixed$condition[[i]]) else sprintf("%s on Plate1 in %s", mixed$plate1_animals[[i]], mixed$condition[[i]])
}, character(1))
animals_or_dash <- function(v) ifelse(is.na(v) | !nzchar(v), "-", gsub(";", ", ", v))
plate_lines <- sprintf("  %-13s weight %+d: Plate1 %d (%s), Plate2 %d (%s)", pc1$condition, as.integer(pc1$contrast_weight),
                       pc1$n_plate1, animals_or_dash(pc1$plate1_animals), pc1$n_plate2, animals_or_dash(pc1$plate2_animals))
weight_text <- if (length(plate_weight) == 1L && abs(plate_weight - 1 / 3) < 1e-9) "one third" else paste(signif(plate_weight, 3), collapse = ", ")
paired_cells <- startsWith(pc1$condition, "paired")
n_aligned <- sum(pc1$n_plate1[paired_cells]) + sum(pc1$n_plate2[!paired_cells])

loo_rows <- data.frame(
  `animal dropped` = ifelse(loo$animal_dropped == "none", "none (full fit)", paste0(loo$animal_dropped, " (", sub("_veh", "", loo$condition_dropped), ")")),
  `FDR<0.05 (all/mapped)` = sprintf("%d / %d", loo$n_fdr05, loo$n_fdr05_mapped),
  `q<0.10 (mapped)` = loo$n_fdr10_mapped, `min q` = fp(loo$min_adj_p_mapped),
  `r(effect)` = sprintf("%.3f", loo$r_effect), `r(t)` = sprintf("%.3f", loo$r_t),
  check.names = FALSE, stringsAsFactors = FALSE
)
if (!is.null(loo_gsea)) {
  loo_rows$`GO-BP (kept of full)` <- sprintf("%d (%d)", loo$go_bp_fdr05, loo$go_bp_retained_of_full)
  loo_rows$`GO-CC (kept of full)` <- sprintf("%d (%d)", loo$go_cc_fdr05, loo$go_cc_retained_of_full)
}
dropped <- loo[loo$animal_dropped != "none", , drop = FALSE]
full_row <- loo[loo$animal_dropped == "none", , drop = FALSE]
worst_effect <- dropped[which.min(dropped$r_effect), , drop = FALSE]
loo_text <- c(
  sprintf("Effect estimates are stable: r(effect) between each leave-one-out fit and the full fit is %s-%s (lowest %.4f, without %s) and r(t) %s-%s.",
          f2(min(dropped$r_effect)), f2(max(dropped$r_effect)), worst_effect$r_effect, worst_effect$animal_dropped, f2(min(dropped$r_t)), f2(max(dropped$r_t))),
  sprintf("The FDR counts depend on single animals: with 2 vs 3 animals (residual df %d) the FDR < 0.05 count ranges from %d to %d proteins (full fit %d) and the q < 0.10 count from %d to %d (full fit %d, mapped).",
          as.integer(min(dropped$df_residual)), min(dropped$n_fdr05), max(dropped$n_fdr05), full_row$n_fdr05,
          min(dropped$n_fdr10_mapped), max(dropped$n_fdr10_mapped), full_row$n_fdr10_mapped),
  if (!is.null(loo_gsea)) sprintf("The GO-BP set also depends on single animals: %d-%d terms, keeping %d-%d of the full fit's %d. GO-CC is more stable: %d-%d terms, keeping %d-%d of %d.",
          min(dropped$go_bp_fdr05), max(dropped$go_bp_fdr05), min(dropped$go_bp_retained_of_full), max(dropped$go_bp_retained_of_full), full_row$go_bp_fdr05,
          min(dropped$go_cc_fdr05), max(dropped$go_cc_fdr05), min(dropped$go_cc_retained_of_full), max(dropped$go_cc_retained_of_full), full_row$go_cc_fdr05)
)

jaccard_wide <- merge(
  rank_robustness[rank_robustness$ontology == "BP", c("comparison", "fdr05_t", "fdr05_log2fc", "fdr05_both", "fdr05_jaccard")],
  rank_robustness[rank_robustness$ontology == "CC", c("comparison", "fdr05_t", "fdr05_log2fc", "fdr05_both", "fdr05_jaccard")],
  by = "comparison", suffixes = c("_BP", "_CC")
)
jaccard_wide <- jaccard_wide[match(primary_contrast_manifest()$canonical_comparison, jaccard_wide$comparison), , drop = FALSE]
fj <- function(v) ifelse(is.na(v), "n/a", sprintf("%.2f", v))
jaccard_rows <- data.frame(
  comparison = sub("_over_.*", "", jaccard_wide$comparison), versus = sub(".*_over_", "", jaccard_wide$comparison),
  `GO-BP t/effect/both` = sprintf("%d/%d/%d", jaccard_wide$fdr05_t_BP, jaccard_wide$fdr05_log2fc_BP, jaccard_wide$fdr05_both_BP),
  `GO-BP Jaccard` = fj(jaccard_wide$fdr05_jaccard_BP),
  `GO-CC t/effect/both` = sprintf("%d/%d/%d", jaccard_wide$fdr05_t_CC, jaccard_wide$fdr05_log2fc_CC, jaccard_wide$fdr05_both_CC),
  `GO-CC Jaccard` = fj(jaccard_wide$fdr05_jaccard_CC),
  check.names = FALSE, stringsAsFactors = FALSE
)
arm_rows <- data.frame(
  compartment = label_of[cno_arms$sample_class], ontology = paste0("GO-", cno_arms$ontology),
  `paired CNO` = cno_arms$paired_cno_fdr05, `unpaired CNO` = cno_arms$unpaired_cno_fdr05,
  `in both` = cno_arms$shared, `opposite sign` = cno_arms$opposite_sign,
  check.names = FALSE, stringsAsFactors = FALSE
)

sig_shift <- dg$sample_class[(dg$go_bp_pct_NES_positive >= 80 | dg$go_bp_pct_NES_positive <= 20)]
report <- c(
  "PAIRING x CNO INTERACTION AND mCHERRY LEARNING ROBUSTNESS, PER-COMPARTMENT MODEL",
  sprintf("Written %s by 04_differential_expression_enrichment/06_mcherry_learning_robustness.r from its own results", timestamp_utc),
  sprintf("and those of 05_pairing_cno_interaction.r (run %s), both in %s.", interaction_run, output_root),
  "",
  "MODEL",
  "- One limma model per sample class: ~0 + condition over that class's 12 animals (n = 3 per condition,",
  "  residual df 8), contrasts.fit, eBayes(robust = TRUE, trend = FALSE), Benjamini-Hochberg within each",
  "  contrast. This is the primary model of 01_preprocessing/02b_limma_stat_gct.r; the refits reproduce the",
  sprintf("  mapped primary t statistics (max difference %s).", fp(max(c(dg$check_max_abs_diff_t_vs_primary, learning_t_diff), na.rm = TRUE))),
  "- Effects are standardized abundance differences in SD units (the logFC / log2fc columns), not log2 ratios.",
  "- Interaction contrast: (paired_cno - paired_veh) - (unpaired_cno - unpaired_veh). Positive = the CNO",
  "  effect is more positive in paired than in unpaired animals. Post hoc, n = 3 animals per cell.",
  sprintf("- GSEA ranked by the moderated t: GO-BP and GO-CC (gseGO), KEGG (gseKEGG); gene sets %d-%d proteins,", gsea_settings$min_gs_size, gsea_settings$max_gs_size),
  sprintf("  BH, FDR < %s. Settings, UniProt/Entrez collapse and per-job seeds as in 01_clusterProfiler.r (base seed", fdr),
  sprintf("  %d, seed = FALSE). Seed ranges in brackets add base seeds %s. KEGG gene sets are downloaded at run time.", config$gsea_seed_base, seed_bases),
  "",
  "1. INTERACTION, PROTEIN LEVEL",
  sprintf("No protein is FDR-significant (q < 0.05) for the Pairing x CNO interaction in any compartment; %s.",
          if (all(dg$n_fdr10 == 0)) "none reaches q < 0.10 either" else sprintf("%d reach q < 0.10", sum(dg$n_fdr10))),
  text_table(protein_rows),
  "Proteins named in the manuscript:",
  named_lines,
  "Top 5 by |t| (effect in SD):",
  top_lines,
  "",
  "2. INTERACTION, GENE SETS (FDR < 0.05; total (up/down), up = CNO effect more positive in paired animals)",
  text_table(gsea_rows),
  "The gene-set results largely follow a compartment-wide shift of the whole interaction contrast rather than a",
  "selective programme: within each compartment most tested terms and most proteins move the same way.",
  text_table(shift_rows),
  if (length(sig_shift)) sprintf("In %s at least 80%% of all tested GO-BP terms have NES of one sign.", paste(label_of[sig_shift], collapse = ", ")),
  "Top terms (up to 3 per ontology and direction):",
  term_lines,
  "",
  "3. COLLECTION PLATE",
  sprintf("Plate per animal (%s), with the interaction contrast weight of each cell:", basename(x$run_manifest$value[x$run_manifest$name == "design_file"])),
  plate_lines,
  sprintf("The interaction contrast keeps %s of any additive Plate2 - Plate1 difference (%s).", weight_text, paste(minority, collapse = "; ")),
  sprintf("Paired animals on Plate1 and unpaired on Plate2 for %d of %d animals, so a Plate x Treatment contribution", n_aligned, sum(pc1$n_animals)),
  "cannot be separated from a Pairing x Treatment one.",
  "",
  "4. mCHERRY LEARNING (paired_veh - unpaired_veh), LEAVE ONE ANIMAL OUT",
  sprintf("Full fit: %d proteins FDR < 0.05 (%d mapped; %d up / %d down), df.prior %s.", full_row$n_fdr05, full_row$n_fdr05_mapped, full_row$n_up, full_row$n_down, f2(full_row$df_prior_median)),
  text_table(loo_rows),
  loo_text,
  "",
  "5. t-RANKED vs EFFECT-SIZE-RANKED GSEA (FDR < 0.05 term sets; Jaccard = both / either)",
  text_table(jaccard_rows),
  "",
  "6. PAIRED CNO vs UNPAIRED CNO: GO terms FDR-significant in both arms",
  text_table(arm_rows),
  sprintf("All %d terms significant in both arms have opposite sign.", sum(cno_arms$shared))[all(cno_arms$opposite_sign == cno_arms$shared)],
  "",
  "7. CHECKS",
  sprintf("  %s: %s (%s)", checks$check, fp(checks$value), ifelse(checks$passed, "passed", "FAILED")),
  "",
  "FILES",
  "  pairing_cno_interaction_*.csv, posthoc_pairing_x_cno_interaction/<class>/  (05_pairing_cno_interaction.r)",
  "  robustness/mcherry_learning_leave_one_out.csv, robustness/mcherry_learning_leave_one_out_gsea/,",
  "  robustness/gsea_rank_robustness_summary.csv, robustness/cno_arm_shared_terms.csv, robustness/checks.csv",
  "  Run manifests with input hashes, software versions and seeds: pairing_cno_interaction_run_manifest.csv,",
  "  robustness/robustness_run_manifest.csv."
)
writeLines(report, report_path, useBytes = TRUE)

# ---- run manifest -------------------------------------------------------------------------
manifest_path <- file.path(robustness_root, "robustness_run_manifest.csv")
written <- c(setdiff(list.files(robustness_root, recursive = TRUE, full.names = TRUE), manifest_path), report_path)
run_manifest <- rbind(
  manifest_rows("run", c("script", "run_timestamp_utc", "output_root", "model", "contrast", "gsea_run"),
                c("04_differential_expression_enrichment/06_mcherry_learning_robustness.r", timestamp_utc, output_root,
                  "per sample class: limma ~0 + condition (12 animals, residual df 8; 7 with one animal left out); eBayes(robust = TRUE, trend = FALSE); BH within contrast",
                  MCHERRY_LEARNING_CONTRAST, run_gsea)),
  do.call(rbind, inputs_read),
  manifest_rows("parameter", c(names(gsea_settings), "gsea_seed_base", "gsea_workers", "ranking_statistic", "enrichment_root"),
                c(vapply(gsea_settings, as.character, character(1)), config$gsea_seed_base, workers, "t", enrichment_root)),
  inferential_software_rows(),
  if (!is.null(gsea_job_rows)) manifest_rows("gsea_job", gsea_job_rows$job_id, gsea_job_rows$seed,
                                             paste0("label=", gsea_job_rows$seed_label, "; analysis=", gsea_job_rows$seed_analysis,
                                                    "; seed_base=", gsea_job_rows$seed_base, "; terms=", gsea_job_rows$n_terms,
                                                    "; fdr05=", gsea_job_rows$n_fdr05)),
  inferential_output_rows(written, output_root)
)
run_manifest <- run_manifest[!duplicated(run_manifest[, c("section", "name")]), , drop = FALSE]
invisible(write_enrichment_csv(run_manifest, manifest_path))

if (!all(checks$passed)) stop("A reproduction check failed; see robustness/checks.csv.", call. = FALSE)
cat(report, sep = "\n")
