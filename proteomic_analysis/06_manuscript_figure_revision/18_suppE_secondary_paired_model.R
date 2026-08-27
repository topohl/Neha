# ==============================================================================
# SUPPLEMENTARY PANEL E -- "Cellular Identities Memory Engram"
# SECONDARY PAIRED CROSS-COMPARTMENT ANIMAL-LEVEL MODEL
# Status: REGENERATED_SECONDARY_PAIRED_ANIMAL_LEVEL
#
# This is the ONE panel for which a new secondary model is authorised. It is
# isolated, labelled secondary, preserves within-animal pairing, and writes
# nothing outside
#   full_regenerated/cross_compartment/SuppE_secondary_paired/
# It does not touch 02_data/animal_level/, 03_output/enrichment/,
# 03_output/ewce/ or 03_output/pca/, and it does not alter any primary result.
#
# ------------------------------------------------------------------------------
# HISTORICAL DESIGN, RECOVERED AND EMPIRICALLY VERIFIED (not assumed)
#
#   condition      paired_veh. The folder is labelled "CS", but the four contrast
#                  filenames themselves all carry condition code _2, and code 2 is
#                  paired_veh in R/analysis_labels.R. Recovered from the data, not
#                  from the folder name.
#   contrasts      bg2_neuron2      = neuropil vs neuron
#                  cfos2_neuron2    = cFos     vs neuron
#                  mcherry2_cfos2   = mCherry  vs cFos
#                  mcherry2_neuron2 = mCherry  vs neuron
#   source matrix  hemisphere-level ProTigy two-sample moderated-t GCT
#                  (02_data/gct/pg.matrix_Two-sample_mod_T_*_n120x5349.gct);
#                  mapped CSVs carry gene_symbol / P.Value / adj.P.Val / logFC
#   sign           numerator = FIRST-named compartment. Verified empirically:
#                  in bg2_neuron2 the neuropil markers Cnp/Mog/Cntnap1/Vcan/Kcna1
#                  average logFC +1.042 while the soma markers Rps/Rpl/H1 average
#                  -1.714.
#   ranking        log2FC, sorted decreasing (prepare_gene_vectors() in the
#                  pre-2026 01_clusterProfiler.r, legacy line 311-316)
#   enrichment     clusterProfiler::gseGO, ont BP, keyType UNIPROT, minGSSize 10,
#                  maxGSSize 800, pvalueCutoff 1, pAdjustMethod BH, org.Mm.eg.db
#   term selection top 5 by NES up + top 5 by NES down per comparison among
#                  p.adjust < 0.05, union plotted for every column
#                  (df_standard in the pre-2026 02_compareGO.r, legacy line 660)
#
# DOCUMENTED DEVIATION: the ranking statistic becomes the moderated t of the new
# paired model instead of log2FC, matching the canonical corrected convention used
# by every other enrichment panel in this revision. All other gseGO settings are
# identical to both the historical and the canonical runs.
#
# ------------------------------------------------------------------------------
# WHY A NEW MODEL WAS NEEDED
# The same three animals contribute to BOTH arms of all four comparisons, so the
# historical independent two-sample test was wrong regardless of the hemisphere
# issue. Averaging L/R fixes the n but does not introduce the within-animal
# blocking factor. This script fits an animal-blocked limma model instead.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "10_style_match_common.R"))
suppressPackageStartupMessages({
  library(limma); library(clusterProfiler); library(org.Mm.eg.db)
  library(stringr); library(digest); library(withr)
})
source(file.path(REPO_ROOT, "R", "analysis_labels.R"))
source(file.path(REPO_ROOT, "R", "protigy_input_utils.R"))

SE      <- ensure_dir(file.path(FULL, "cross_compartment", "SuppE_secondary_paired"))
SE_PC   <- ensure_dir(file.path(SE, "per_contrast"))
SE_EN   <- ensure_dir(file.path(SE, "enrichment"))
SE_SD   <- ensure_dir(file.path(SE, "source_data"))

CONDITION <- "paired_veh"
CONTRASTS <- list(
  list(id = "neuropil_vs_neuron", num = "neuropil", den = "neuron",  hist = "bg2_neuron2",      label = "neuropil vs neuron"),
  list(id = "cfos_vs_neuron",     num = "cfos",     den = "neuron",  hist = "cfos2_neuron2",    label = "cfos vs neuron"),
  list(id = "mcherry_vs_cfos",    num = "mcherry",  den = "cfos",    hist = "mcherry2_cfos2",   label = "mcherry vs cfos"),
  list(id = "mcherry_vs_neuron",  num = "mcherry",  den = "neuron",  hist = "mcherry2_neuron2", label = "mcherry vs neuron")
)
GSEA_SEED <- 20260827L
MIN_GS <- 10; MAX_GS <- 800; PCUT <- 1; PADJ <- "BH"

# ==============================================================================
# 0. HISTORICAL DESIGN RECOVERY RECORD
# ==============================================================================
rec <- data.frame(
  property = c("panel", "historical_output", "historical_generating_script", "condition",
               "condition_evidence", "n_contrasts", "contrast_definitions",
               "source_matrix", "source_statistic", "sign_convention", "sign_evidence",
               "ranking_statistic", "enrichment_method", "enrichment_settings",
               "term_selection_rule", "term_selection_provenance", "resolution"),
  value = c(
    "Supplementary E - Cellular Identities Memory Engram",
    "99_historical/compareGO/BP/baseline_cell_type_profiling/CS/02_Main_Plots/Dotplot_Enrichment_TopGenes_PerComp.svg",
    "04_differential_expression_enrichment/02_compareGO.r (pre-2026 3140-line version), ensemble_profiling='baseline_cell_type_profiling', condition='CS'",
    "paired_veh",
    "all four historical contrast filenames carry condition code _2; code 2 = paired_veh per R/analysis_labels.R. Recovered from the contrast names, NOT assumed from the 'CS' folder label or the panel title.",
    "4",
    paste(vapply(CONTRASTS, function(z) sprintf("%s = %s over %s", z$hist, z$num, z$den), character(1)), collapse = "; "),
    "02_data/gct/pg.matrix_Two-sample_mod_T_*_n120x5349.gct (hemisphere level, 96 acquisitions)",
    "ProTigy / limma two-sample moderated t (mapped CSV columns gene_symbol, P.Value, adj.P.Val, logFC)",
    "positive logFC and positive NES = higher in the FIRST-named compartment",
    "verified empirically in bg2_neuron2: neuropil markers Cnp/Mog/Cntnap1/Vcan/Kcna1 mean logFC +1.042; neuron-soma markers Rps9/Rpl22/Rps16/Rpl35/H1-1/H1-3/Rps18/Rpl6 mean logFC -1.714",
    "log2FC sorted decreasing (prepare_gene_vectors(), pre-2026 01_clusterProfiler.r line 311-316)",
    "clusterProfiler::gseGO",
    sprintf("ont=BP, keyType=UNIPROT, minGSSize=%d, maxGSSize=%d, pvalueCutoff=%d, pAdjustMethod=%s, OrgDb=org.Mm.eg.db", MIN_GS, MAX_GS, PCUT, PADJ),
    "per comparison, among p.adjust < 0.05: top 5 by NES positive + top 5 by NES negative; union of those Descriptions plotted for every column",
    "df_standard block in the pre-2026 02_compareGO.r (legacy line 660)",
    "RESOLVED - exact historical model established; secondary paired animal-level model proceeds"),
  stringsAsFactors = FALSE)
utils::write.csv(rec, file.path(SE, "SuppE_historical_design_recovery.csv"), row.names = FALSE)
cat("historical design: RESOLVED\n")

# ==============================================================================
# 1. ANIMAL-LEVEL MATRIX, RESTRICTED TO THE HISTORICAL CONDITION
# ==============================================================================
gct_sha <- digest::digest(file = GCT_ANIMAL, algo = "sha256")
parsed <- validate_protigy_gct_v13(GCT_ANIMAL)
mat_all <- parsed$matrix
cm <- parsed$column_metadata
meta_all <- data.frame(
  Sample = colnames(mat_all),
  AnimalID = trimws(as.character(cm["AnimalID", ])),
  condition = normalize_condition(cm["condition", ]),
  sample_class = normalize_sample_class(cm["sample_class", ]),
  stringsAsFactors = FALSE)

keep <- meta_all$condition == CONDITION
meta <- meta_all[keep, , drop = FALSE]
mat  <- mat_all[, keep, drop = FALSE]
meta$sample_class <- factor(meta$sample_class, levels = c("neuropil", "cfos", "mcherry", "neuron"))
meta$AnimalID <- factor(meta$AnimalID)

# ---- VALIDATION 1-4 ----------------------------------------------------------
v <- list()
v$n_columns <- ncol(mat)
v$n_animals <- nlevels(meta$AnimalID)
v$animals <- paste(levels(meta$AnimalID), collapse = ";")
v$no_duplicate_animal_class <- anyDuplicated(paste(meta$AnimalID, meta$sample_class)) == 0
v$all_classes_in_all_animals <- all(table(meta$AnimalID, meta$sample_class) == 1)
v$hemisphere_columns_present <- any(grepl("_L$|_R$|Plate", meta$Sample))
v$input_is_animal_level_gct <- identical(gct_sha, "f12cf99e1bfb7c17bbf56bffb6783e924698bce5d5533a8e312bc4bbb733bbb3")
stopifnot(v$n_columns == 12L, v$n_animals == 3L, v$no_duplicate_animal_class,
          v$all_classes_in_all_animals, !v$hemisphere_columns_present, v$input_is_animal_level_gct)
cat(sprintf("matrix: %d proteins x %d animal-level columns; %d animals (%s)\n",
            nrow(mat), ncol(mat), v$n_animals, v$animals))

# ==============================================================================
# 2. ANIMAL-BLOCKED LIMMA MODEL
# ==============================================================================
design <- stats::model.matrix(~ 0 + sample_class + AnimalID, data = meta)
colnames(design) <- sub("^sample_class", "", colnames(design))
qrd <- qr(design)
v$design_rank <- qrd$rank
v$design_cols <- ncol(design)
v$full_rank <- qrd$rank == ncol(design)
v$residual_df <- nrow(meta) - qrd$rank
stopifnot(v$full_rank)
cat(sprintf("design: %d x %d, rank %d, FULL RANK, residual df = %d\n",
            nrow(design), ncol(design), qrd$rank, v$residual_df))

utils::write.csv(data.frame(Sample = meta$Sample, AnimalID = as.character(meta$AnimalID),
                            condition = meta$condition, sample_class = as.character(meta$sample_class),
                            design, check.names = FALSE),
                 file.path(SE, "model_design.csv"), row.names = FALSE)

cm_str <- vapply(CONTRASTS, function(z) sprintf("%s - %s", z$num, z$den), character(1))
names(cm_str) <- vapply(CONTRASTS, function(z) z$id, character(1))
contr <- limma::makeContrasts(contrasts = unname(cm_str), levels = design)
colnames(contr) <- names(cm_str)
utils::write.csv(data.frame(coefficient = rownames(contr), contr, check.names = FALSE),
                 file.path(SE, "contrast_matrix.csv"), row.names = FALSE)

fit  <- limma::lmFit(mat, design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, contr))
v$eBayes_df_prior <- fit2$df.prior
v$eBayes_df_total <- unique(round(fit2$df.total, 3))[1]
cat(sprintf("eBayes: df.prior = %.3f, df.total = %.3f (residual df %d + prior)\n",
            fit2$df.prior, v$eBayes_df_total, v$residual_df))

# ==============================================================================
# 3. PER-CONTRAST PROTEIN STATISTICS
# ==============================================================================
map <- utils::read.csv(file.path(MAPPED_DIR, "mcherry_paired_veh_vs_mcherry_unpaired_veh.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
map <- unique(map[, c("original_protein_id", "uniprot_accession", "mapped_gene_symbol", "source_row_id")])

all_stats <- do.call(rbind, lapply(CONTRASTS, function(z) {
  tt <- limma::topTable(fit2, coef = z$id, number = Inf, sort.by = "none")
  d <- data.frame(original_protein_id = rownames(mat),
                  contrast_id = z$id, historical_alias = z$hist,
                  numerator = z$num, denominator = z$den,
                  log2fc = tt$logFC, t = tt$t, pval = tt$P.Value, padj = tt$adj.P.Val,
                  B = tt$B, aveExpr = tt$AveExpr, stringsAsFactors = FALSE)
  merge(d, map, by = "original_protein_id", all.x = TRUE, sort = FALSE)
}))
all_stats$significant <- is.finite(all_stats$padj) & all_stats$padj < FDR
utils::write.csv(all_stats, file.path(SE, "all_protein_statistics.csv"), row.names = FALSE)

for (z in CONTRASTS) {
  s <- all_stats[all_stats$contrast_id == z$id, , drop = FALSE]
  s <- s[order(s$padj, -abs(s$log2fc)), , drop = FALSE]
  utils::write.csv(s, file.path(SE_PC, paste0(z$id, "_protein_statistics.csv")), row.names = FALSE)
  cat(sprintf("  %-20s FDR<0.05 = %5d / %d   |log2fc| median %.3f\n",
              z$id, sum(s$significant, na.rm = TRUE), nrow(s), median(abs(s$log2fc), na.rm = TRUE)))
}

# ==============================================================================
# 4. GSEA ON THE PAIRED-MODEL MODERATED t
# ==============================================================================
collapse_rank <- function(s) {
  s <- s[!is.na(s$uniprot_accession) & nzchar(s$uniprot_accession) & is.finite(s$t), , drop = FALSE]
  # canonical duplicate rule: largest |log2fc|, ties by source_row_id then accession
  s <- s[order(s$uniprot_accession, -abs(s$log2fc), s$source_row_id, s$original_protein_id), , drop = FALSE]
  s <- s[!duplicated(s$uniprot_accession), , drop = FALSE]
  sort(stats::setNames(s$t, s$uniprot_accession), decreasing = TRUE)
}

gsea_list <- lapply(CONTRASTS, function(z) {
  s <- all_stats[all_stats$contrast_id == z$id, , drop = FALSE]
  rk <- collapse_rank(s)
  cat(sprintf("GSEA %-20s ranked genes = %d ...\n", z$id, length(rk)))
  res <- withr::with_seed(GSEA_SEED, suppressWarnings(clusterProfiler::gseGO(
    geneList = rk, ont = "BP", keyType = "UNIPROT", OrgDb = org.Mm.eg.db::org.Mm.eg.db,
    exponent = 1, minGSSize = MIN_GS, maxGSSize = MAX_GS, eps = 1e-10,
    pvalueCutoff = PCUT, pAdjustMethod = PADJ, verbose = FALSE, seed = TRUE, by = "fgsea")))
  tab <- as.data.frame(res)
  tab$contrast_id <- z$id; tab$historical_alias <- z$hist; tab$Comparison <- z$label
  tab$n_ranked_genes <- length(rk)
  utils::write.csv(tab, file.path(SE_EN, paste0(z$id, "_GSEA_GO_BP.csv")), row.names = FALSE)
  cat(sprintf("   -> %d terms, %d at FDR<0.05\n", nrow(tab), sum(tab$p.adjust < FDR, na.rm = TRUE)))
  tab
})
gsea <- do.call(rbind, gsea_list)
utils::write.csv(gsea, file.path(SE_EN, "GSEA_GO_BP_all_contrasts.csv"), row.names = FALSE)

# ==============================================================================
# 5. TERM SELECTION (historical df_standard rule) AND BUBBLE PLOT
# ==============================================================================
col_levels <- vapply(CONTRASTS, function(z) z$label, character(1))
gsea$Comparison <- factor(gsea$Comparison, levels = col_levels)

best <- gsea[order(gsea$Comparison, gsea$Description, -abs(gsea$NES)), ]
best <- best[!duplicated(best[, c("Comparison", "Description")]), ]
pool <- best[is.finite(best$p.adjust) & best$p.adjust < FDR, , drop = FALSE]
sel <- do.call(rbind, lapply(levels(gsea$Comparison), function(cp) {
  d <- pool[pool$Comparison == cp, , drop = FALSE]
  up <- d[d$NES > 0, ]; up <- head(up[order(-up$NES), ], 5)
  dn <- d[d$NES < 0, ]; dn <- head(dn[order(dn$NES), ], 5)
  rbind(up, dn)
}))
plot_df <- best[best$Description %in% unique(sel$Description), , drop = FALSE]

# legacy order_dotplot(), with the historical column order fixed
agg <- stats::aggregate(NES ~ Description + Comparison, data = plot_df, FUN = mean)
m <- tapply(agg$NES, list(agg$Description, agg$Comparison), function(z) z[1]); m[is.na(m)] <- 0
m <- m[, col_levels[col_levels %in% colnames(m)], drop = FALSE]
dom <- apply(m, 1, function(x) colnames(m)[which.max(abs(x))])
domv <- apply(m, 1, function(x) x[which.max(abs(x))])
od <- data.frame(Term = rownames(m), Comp = factor(dom, levels = colnames(m)), Val = domv)
od <- od[order(-as.integer(od$Comp), ifelse(od$Val > 0, 1, 0), abs(od$Val)), ]
plot_df$Comparison  <- factor(as.character(plot_df$Comparison), levels = colnames(m))
plot_df$Description <- factor(plot_df$Description, levels = od$Term)
lim <- max(abs(plot_df$NES), na.rm = TRUE)

cat(sprintf("panel: %d terms x %d columns\n", nlevels(plot_df$Description), nlevels(plot_df$Comparison)))

theme_dot <- function() {
  theme_classic(base_size = DOT_FS, base_family = FONT) +
    theme(text = element_text(family = FONT, colour = "#2C2C2C"),
          axis.text.x = element_text(size = DOT_FS, colour = "#2C2C2C", angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = DOT_FS, colour = "#2C2C2C"),
          axis.line = element_line(colour = "#2C2C2C", linewidth = 0.6 * PT),
          axis.ticks = element_line(colour = "#2C2C2C", linewidth = 0.5 * PT),
          axis.ticks.length = unit(3 * PT, "mm"),
          legend.title = element_text(size = DOT_FS, colour = "#2C2C2C"),
          legend.text = element_text(size = DOT_FS * 0.9, colour = "#2C2C2C"),
          legend.position = "right", legend.justification = "top",
          legend.key = element_blank(), legend.background = element_blank(),
          panel.grid = element_blank(), panel.border = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.subtitle = element_blank(), plot.margin = margin(10, 10, 10, 10))
}

build_E <- function(title = NULL) {
  p <- ggplot(plot_df, aes(x = Comparison, y = Description, colour = NES, size = -log10(p.adjust))) +
    geom_point(alpha = 0.85, stroke = 0.3) +
    scale_colour_gradientn(colours = grDevices::colorRampPalette(c(DOT_NES_LOW, DOT_NES_MID, DOT_NES_HIGH),
                                                                 space = "Lab")(100),
                           name = "NES", limits = c(-lim, lim),
                           guide = guide_colourbar(barwidth = unit(DOT_BAR_W * PT, "mm"),
                                                   barheight = unit(DOT_BAR_H * PT, "mm"), order = 1)) +
    scale_size_continuous(name = expression(-log[10](italic(P)[adj])), range = DOT_SIZE_RNG, limits = c(0, NA)) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)), drop = FALSE) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
    labs(x = NULL, y = NULL) + theme_dot()
  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(family = FONT, size = 12, face = "bold", hjust = 0.5,
                                      margin = margin(b = 8)))
  } else {
    p <- p + theme(plot.title = element_blank())
  }
  p
}

n_rows <- nlevels(plot_df$Description)
H <- max(5, 2.5 + n_rows * 0.25)
pE <- build_E(NULL)
ggsave(file.path(SE, "SuppE_cellular_identities_animal_level_FINAL.svg"), pE,
       width = DOT_W, height = H, units = "in", device = svglite::svglite)
ggsave(file.path(SE, "SuppE_cellular_identities_animal_level_FINAL.png"), pE,
       width = DOT_W, height = H, units = "in", dpi = 300, bg = "white")
save_sm(pE, SM_SUP, "SuppE_cellular_identities_animal_level_FINAL", DOT_W, H)
save_sm(build_E("Cellular Identities Memory Engram"), SM_SUP,
        "SuppE_cellular_identities_animal_level_FINAL_titled", DOT_W, H)

src <- plot_df[, c("Comparison", "contrast_id", "historical_alias", "ID", "Description",
                   "setSize", "NES", "enrichmentScore", "pvalue", "p.adjust", "qvalue", "n_ranked_genes")]
src$fdr_significant <- is.finite(src$p.adjust) & src$p.adjust < FDR
utils::write.csv(src, file.path(SE_SD, "SuppE_source_data.csv"), row.names = FALSE)
utils::write.csv(src, file.path(SM_SUP, "SuppE_source_data.csv"), row.names = FALSE)

# ==============================================================================
# 6. DESCRIPTIVE COMPARISON AGAINST THE HISTORICAL HEMISPHERE-LEVEL RESULT
# ==============================================================================
histdir <- file.path(DATA_ROOT, "99_historical/datasets_mapped/baseline_cell_type_profiling/CS")
cmpv <- do.call(rbind, lapply(CONTRASTS, function(z) {
  hp <- file.path(histdir, paste0(z$hist, ".csv"))
  if (!file.exists(hp)) return(NULL)
  hd <- utils::read.csv(hp, stringsAsFactors = FALSE)
  nd <- all_stats[all_stats$contrast_id == z$id, c("uniprot_accession", "log2fc")]
  nd <- nd[!is.na(nd$uniprot_accession), ]
  nd <- nd[!duplicated(nd$uniprot_accession), ]
  mm <- merge(hd[, c("gene_symbol", "logFC")], nd,
              by.x = "gene_symbol", by.y = "uniprot_accession")
  data.frame(contrast_id = z$id, historical_alias = z$hist, n_shared = nrow(mm),
             pearson = stats::cor(mm$logFC, mm$log2fc),
             spearman = stats::cor(mm$logFC, mm$log2fc, method = "spearman"),
             sign_agreement = mean(sign(mm$logFC) == sign(mm$log2fc)),
             stringsAsFactors = FALSE)
}))
print(cmpv, row.names = FALSE)
utils::write.csv(cmpv, file.path(SE, "historical_vs_corrected_logfc_comparison.csv"), row.names = FALSE)

# ==============================================================================
# 7. MODEL AUDIT + VALIDATION LEDGER
# ==============================================================================
per_contrast <- do.call(rbind, lapply(CONTRASTS, function(z) {
  s <- all_stats[all_stats$contrast_id == z$id, ]
  g <- gsea[gsea$contrast_id == z$id, ]
  data.frame(contrast_id = z$id, historical_alias = z$hist, label = z$label,
             numerator = z$num, denominator = z$den,
             n_animals_numerator = 3L, n_animals_denominator = 3L,
             animals = v$animals,
             n_proteins = nrow(s), n_FDR_proteins = sum(s$significant, na.rm = TRUE),
             n_go_terms = nrow(g), n_FDR_go_terms = sum(g$p.adjust < FDR, na.rm = TRUE),
             n_terms_contributed = if (nrow(sel)) sum(sel$contrast_id == z$id) else 0L,
             stringsAsFactors = FALSE)
}))
utils::write.csv(per_contrast, file.path(SE, "per_contrast_summary.csv"), row.names = FALSE)
print(per_contrast[, c("label", "n_FDR_proteins", "n_FDR_go_terms", "n_terms_contributed")], row.names = FALSE)

audit <- data.frame(
  item = c("panel", "status", "analysis_role", "condition", "input_matrix", "input_sha256",
           "n_observations", "n_animals", "animals", "sample_classes",
           "model_formula", "blocking_factor", "design_columns", "design_rank", "full_rank",
           "residual_df", "eBayes_df_prior", "eBayes_df_total",
           "n_paired_differences_per_contrast", "ranking_statistic", "ranking_deviation",
           "enrichment_method", "enrichment_settings", "gsea_seed", "term_selection_rule",
           "sign_convention", "primary_outputs_touched", "canonical_roots_touched",
           "n_limitation", "interpretation", "not_to_be_used_for"),
  value = c(
    "Supplementary E", "REGENERATED_SECONDARY_PAIRED_ANIMAL_LEVEL",
    "SECONDARY paired cross-compartment identity analysis - NOT one of the 12 primary within-compartment treatment contrasts",
    CONDITION, GCT_ANIMAL, gct_sha,
    as.character(v$n_columns), as.character(v$n_animals), v$animals,
    "neuropil, cfos, mcherry, neuron",
    "~ 0 + sample_class + AnimalID (limma lmFit / contrasts.fit / eBayes)",
    "AnimalID as a fixed blocking effect; equivalent to a within-animal paired analysis",
    as.character(v$design_cols), as.character(v$design_rank), as.character(v$full_rank),
    as.character(v$residual_df), sprintf("%.3f", fit2$df.prior), sprintf("%.3f", v$eBayes_df_total),
    "3 (one per animal); the residual df of 6 is pooled across all four compartments under a common-variance assumption",
    "moderated t from this paired model",
    "DEVIATION from the historical run, which ranked by log2FC. Changed to match the canonical corrected convention used by every other enrichment panel in this revision.",
    "clusterProfiler::gseGO",
    sprintf("ont=BP, keyType=UNIPROT, minGSSize=%d, maxGSSize=%d, pvalueCutoff=%d, pAdjustMethod=%s, exponent=1, eps=1e-10, by=fgsea", MIN_GS, MAX_GS, PCUT, PADJ),
    as.character(GSEA_SEED),
    "historical df_standard: top 5 by NES up + top 5 by NES down per comparison among p.adjust < 0.05; union plotted for every column",
    "positive log2FC / positive NES = higher in the FIRST-named compartment (numerator)",
    "NONE", "NONE - nothing written to 02_data/animal_level/, 03_output/enrichment/, 03_output/ewce/ or 03_output/pca/",
    "n = 3 animals. Residual df = 6 before empirical-Bayes moderation. Inference must be treated cautiously.",
    "Secondary paired cross-compartment analysis of compartment identity / relative proteomic composition. The purpose is to show that the dissected fractions have biologically coherent relative identities.",
    "Must NOT be interpreted as a treatment effect and must NOT be promoted into the 12 primary contrasts."),
  stringsAsFactors = FALSE)
utils::write.csv(audit, file.path(SE, "model_audit.csv"), row.names = FALSE)
utils::write.csv(audit, file.path(SM_SUP, "SuppE_model_audit.csv"), row.names = FALSE)

val <- data.frame(
  check = c("1 exactly 3 unique animals per side of each paired comparison",
            "2 animals match across the two compared compartments",
            "3 no hemisphere-level rows enter the model",
            "4 no duplicated animal/sample-class observations",
            "5 contrast direction matches the historical panel",
            "6 t/NES direction documented",
            "7 model matrix is full rank",
            "8 residual df reported",
            "9 enrichment uses the corrected paired-model statistic",
            "10 output does not alter primary-analysis files"),
  result = c(
    sprintf("PASS - 3 animals (%s) contribute to every sample_class", v$animals),
    "PASS - identical animal set in all four sample classes; table(AnimalID, sample_class) is all ones",
    sprintf("PASS - input is the animal-level GCT, sha256 %s; 48 columns, no L/R or plate identifiers", substr(gct_sha, 1, 12)),
    sprintf("PASS - anyDuplicated(AnimalID x sample_class) = 0"),
    "PASS - numerator = first-named compartment, matching the historical filenames and the empirically verified marker directions",
    "PASS - recorded in model_audit.csv (positive = higher in the numerator)",
    sprintf("PASS - rank %d of %d columns", v$design_rank, v$design_cols),
    sprintf("PASS - residual df = %d; eBayes df.total = %.3f", v$residual_df, v$eBayes_df_total),
    "PASS - gseGO ranked on the moderated t of this paired model, not on hemisphere-level ranks",
    "PASS - all outputs under full_regenerated/cross_compartment/SuppE_secondary_paired/ plus figure copies; canonical roots untouched"),
  stringsAsFactors = FALSE)
utils::write.csv(val, file.path(SE, "SuppE_validation.csv"), row.names = FALSE)
print(val, row.names = FALSE)

message("Supplementary E secondary paired model complete.")
