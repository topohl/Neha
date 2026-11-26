# ==============================================================================
# Reviewer revision 2026-08-27 -- shared configuration
#
# Purpose: regenerate the EXISTING manuscript proteomics figure panels using the
# corrected ANIMAL-LEVEL analysis. This is a statistics correction, not a figure
# redesign: plot type, palette, axis conventions and labels are preserved from
# the original panels wherever the corrected evidence still supports them.
#
# Panel provenance is recorded in
#   manifests/original_figure_panel_source_map.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(svglite)
})

DATA_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
REPO_ROOT <- "C:/Users/topohl/Documents/GitHub/Neha/proteomic_analysis"

OUT_ROOT   <- file.path(DATA_ROOT, "03_output", "figure_rebuild_within_compartment_20261126")
FULL       <- file.path(OUT_ROOT, "full_regenerated")
PANELS_F3  <- file.path(OUT_ROOT, "figure_panels", "Figure3")
PANELS_SUP <- file.path(OUT_ROOT, "figure_panels", "Supplementary_proteomics")
MANIFESTS  <- file.path(OUT_ROOT, "manifests")

# ---- canonical corrected inputs (CANONICAL_OUTPUTS.md) -----------------------
GCT_ANIMAL   <- file.path(DATA_ROOT, "02_data/animal_level/input_gct/neha_protigy_input_animal_level_primary.gct")
MAPPED_DIR   <- file.path(DATA_ROOT, "02_data/animal_level/within_compartment_20261124/mapped/forward")
# GO GSEA from the nPermSimple = 100000 rerun (2026-11-25). KEGG was not rerun and
# stays in the 2026-11-24 run (KEGG_ROOT).
ENRICH_ROOT  <- file.path(DATA_ROOT, "03_output/enrichment/enrichment_within_compartment_20261125/GO_BP/per_comparison")
KEGG_ROOT    <- file.path(DATA_ROOT, "03_output/enrichment/enrichment_within_compartment_20261124/GO_BP/per_comparison")
EWCE_ROOT    <- file.path(DATA_ROOT, "03_output/ewce/EWCE_Results_animal_level_validation_20260825")
PCA_ROOT     <- file.path(DATA_ROOT, "03_output/pca/pca_plots_animal_level_validation_20260825_rerun")

FDR <- 0.05

# ---- shared labels -----------------------------------------------------------
SAMPLE_CLASSES <- c("neuropil", "cfos", "mcherry", "neuron")
CLASS_LABELS   <- c(neuropil = "Neuropil", cfos = "cFos", mcherry = "mCherry", neuron = "Neuron")

# Legacy 04_compare_sig_expr.r / 03_compare_pathways.r palette, kept verbatim so
# the corrected panels sit next to the original artwork without a colour shift.
CLASS_COLORS <- c(neuropil = "#4c87c6", cfos = "#6ccff6", mcherry = "#faa51a", neuron = "#fdd700")

# Legacy PCA palette (06a_pca_core.r base_hex), first four entries.
PCA_PALETTE <- c("#1E90FF", "#61d7ffff", "#FFA500", "#FFD700")

CONDITIONS <- c("paired_cno", "paired_veh", "unpaired_cno", "unpaired_veh")
COND_LABELS <- c(paired_cno = "paired-cno", paired_veh = "paired-veh",
                 unpaired_cno = "unpaired-cno", unpaired_veh = "unpaired-veh")

# ---- helpers -----------------------------------------------------------------
ensure_dir <- function(p) { dir.create(p, recursive = TRUE, showWarnings = FALSE); p }

# Every output path goes through fresh(): an existing file is never replaced
# unless FIGURE_REBUILD_FORCE=true.
FORCE <- identical(tolower(Sys.getenv("FIGURE_REBUILD_FORCE")), "true")
fresh <- function(p) {
  if (file.exists(p) && !FORCE)
    stop("Refusing to overwrite ", p, " (set FIGURE_REBUILD_FORCE=true to replace it).", call. = FALSE)
  p
}

save_panel <- function(plot, stage_dir, panel_dir, filename, width, height) {
  ensure_dir(stage_dir); ensure_dir(panel_dir)
  stage_path <- fresh(file.path(stage_dir, filename))
  ggplot2::ggsave(stage_path, plot = plot, width = width, height = height,
                  units = "in", device = svglite::svglite)
  panel_path <- file.path(panel_dir, filename)
  if (!identical(normalizePath(stage_path, winslash = "/", mustWork = FALSE),
                 normalizePath(panel_path, winslash = "/", mustWork = FALSE))) {
    file.copy(stage_path, fresh(panel_path), overwrite = FORCE)
  }
  message("  saved: ", stage_path, "  -> ", panel_path)
  # PNG preview for quick visual checking; SVG remains the deliverable
  prev_dir <- ensure_dir(file.path(OUT_ROOT, "previews"))
  prev <- fresh(file.path(prev_dir, sub("\\.svg$", ".png", filename)))
  try(ggplot2::ggsave(prev, plot = plot, width = width, height = height, units = "in", dpi = 150),
      silent = TRUE)
  invisible(c(stage = stage_path, panel = panel_path))
}

save_source_data <- function(df, stage_dir, panel_dir, filename) {
  ensure_dir(stage_dir); ensure_dir(panel_dir)
  stage_path <- fresh(file.path(stage_dir, filename))
  utils::write.csv(df, stage_path, row.names = FALSE)
  panel_path <- file.path(panel_dir, filename)
  if (!identical(normalizePath(stage_path, winslash = "/", mustWork = FALSE),
                 normalizePath(panel_path, winslash = "/", mustWork = FALSE))) {
    file.copy(stage_path, fresh(panel_path), overwrite = FORCE)
  }
  invisible(c(stage = stage_path, panel = panel_path))
}

read_mapped <- function(sample_class, numerator, denominator) {
  f <- file.path(MAPPED_DIR, sprintf("%s_%s_vs_%s_%s.csv", sample_class, numerator, sample_class, denominator))
  if (!file.exists(f)) stop("Canonical mapped contrast not found: ", f, call. = FALSE)
  d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  d$sample_class <- sample_class
  d$canonical_contrast <- sprintf("%s_%s_over_%s_%s", sample_class, numerator, sample_class, denominator)
  d
}

read_gsea <- function(sample_class, numerator, denominator) {
  f <- file.path(ENRICH_ROOT, sprintf("%s_%s_over_%s_%s", sample_class, numerator, sample_class, denominator), "GSEA_GO_BP.csv")
  if (!file.exists(f)) stop("Canonical GSEA table not found: ", f, call. = FALSE)
  d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  d$sample_class <- sample_class
  d$canonical_contrast <- sprintf("%s_%s_over_%s_%s", sample_class, numerator, sample_class, denominator)
  d
}

# Gene sets with identical members are shown once (2026-11-26) -------------------
# Two GO-BP terms are identical when their member proteins among the proteins ranked
# in the GSEA are the same set. Members come from the GO annotation the enrichment
# used (clusterProfiler's GO data for org.Mm.eg.db, ont BP, keyType UNIPROT; a term
# carries the proteins annotated to it or to its descendants), restricted to the
# proteins ranked for the comparison (GO_BP/audits/<comparison>/gsea_rank_audit.csv,
# included_in_gsea). Every term's member count must equal its setSize in the GSEA
# table, which shows that the annotation is the one the enrichment used. Of an
# identical group the term with the shorter name is kept (tie: lower GO ID).
.go_bp <- new.env(parent = emptyenv())
go_bp_members <- function() {
  if (is.null(.go_bp$p2e)) {
    gd <- suppressMessages(clusterProfiler:::get_GO_data(org.Mm.eg.db::org.Mm.eg.db, "BP", "UNIPROT"))
    .go_bp$p2e <- get("PATHID2EXTID", envir = gd)
  }
  .go_bp$p2e
}

gsea_ranked_proteins <- function(comparison) {
  f <- file.path(dirname(ENRICH_ROOT), "audits", comparison, "gsea_rank_audit.csv")
  if (!file.exists(f)) stop("GSEA rank audit not found: ", f, call. = FALSE)
  a <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  unique(a$uniprot_accession[a$included_in_gsea %in% c(TRUE, "TRUE")])
}

# terms: ID, Description and setSize (NA when the term is not in that comparison's
# GSEA table) of one comparison. Returns one row per term that has an identical
# twin, with the term kept for its group; dropped = TRUE for every other member.
identical_member_sets <- function(terms, comparison) {
  terms <- unique(terms[, c("ID", "Description", "setSize")])
  stopifnot(!anyDuplicated(terms$ID))
  ranked <- gsea_ranked_proteins(comparison)
  p2e <- go_bp_members()
  mem <- lapply(stats::setNames(terms$ID, terms$ID), function(id) sort(intersect(p2e[[id]], ranked)))
  n <- unname(lengths(mem))
  bad <- !is.na(terms$setSize) & n != terms$setSize
  if (any(bad)) stop("GO annotation does not reproduce setSize in ", comparison, " for ",
                     paste(terms$ID[bad], collapse = ", "), call. = FALSE)
  key <- unname(vapply(mem, paste, character(1), collapse = ";"))
  twin <- key %in% key[duplicated(key)]
  out <- data.frame(comparison = rep(comparison, sum(twin)), ID = terms$ID[twin],
                    Description = terms$Description[twin], n_members = n[twin],
                    key = key[twin], stringsAsFactors = FALSE)
  out <- out[order(out$key, nchar(out$Description), out$ID), , drop = FALSE]
  if (!nrow(out)) {
    out$kept_ID <- character(0); out$kept_Description <- character(0); out$dropped <- logical(0)
    out$key <- NULL
    return(out)
  }
  out$kept_ID <- stats::ave(out$ID, out$key, FUN = function(z) z[1])
  out$kept_Description <- out$Description[match(out$kept_ID, out$ID)]
  out$dropped <- out$ID != out$kept_ID
  out$key <- NULL
  rownames(out) <- NULL
  out
}

# Provenance ledger ------------------------------------------------------------
.provenance <- new.env(parent = emptyenv())
.provenance$rows <- list()

record_provenance <- function(figure, panel, status, plot_type, corrected_input,
                              output_path, source_data_path, displayed_stats_changed,
                              layout_preserved, interpretation_changed, notes) {
  .provenance$rows[[length(.provenance$rows) + 1L]] <- data.frame(
    figure = figure, panel = panel, status = status, plot_type = plot_type,
    corrected_input = corrected_input, output_path = output_path,
    source_data_path = source_data_path,
    displayed_statistics_changed = displayed_stats_changed,
    visual_layout_preserved = layout_preserved,
    biological_interpretation_changed = interpretation_changed,
    notes = notes,
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

write_provenance <- function(filename) {
  df <- do.call(rbind, .provenance$rows)
  p <- file.path(ensure_dir(MANIFESTS), filename)
  utils::write.csv(df, p, row.names = FALSE)
  message("provenance written: ", p)
  invisible(p)
}
