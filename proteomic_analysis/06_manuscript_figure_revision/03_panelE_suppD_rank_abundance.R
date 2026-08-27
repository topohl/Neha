# ==============================================================================
# Figure 3 panel E        -- Rank-abundance, mCherry paired-veh / unpaired-veh
# Supplementary panel D   -- Rank-abundance, neuron / neuropil unpaired-veh
#
# ORIGINAL: exact generating script SOURCE_NOT_FOUND (absent from the repository
#           and from every git commit). Recovered from its own outputs:
#             plot   03_output/qc/Rank-Abundance_Marker_QC.svg   (2026-04-14)
#             data   03_output/qc/Processed_Protein_Ranks.xlsx   (16 sheets,
#                    one per sample_class x condition, columns Condition / Rank /
#                    Genes / MeanLog2 / LinearValue / MarkerType)
#           Nearest surviving relatives: 03_qc_exploration/02_rank_abundance_by_
#           sample_class.r (current, per sample_class only) and the deleted
#           03_qc_exploration/02_rank_abundance_plot_E9.r (Exp9 ancestor).
#
# ANALYSIS TYPE: DESCRIPTIVE. The curve is a mean-abundance rank ordering with
#           marker labels. No test, no p-value, no inferential n. The only thing
#           the animal-level correction changes is the averaging: the original
#           averaged over hemisphere-level acquisitions within a group, the
#           corrected version averages the L/R-averaged animal-level values.
#
# Everything reconstructed verbatim from the original artefacts:
#   * transformation  MeanLog2 per group -> LinearValue = 2^MeanLog2 -> rank desc
#   * axes            "Protein Rank (ordered by abundance)" /
#                     "Abundance (Intensity, Log10 Scale)"
#   * title/subtitle  "Rank-Abundance Distribution and Marker Validation" /
#                     "Log2-transformed input data re-linearized for Log10 visualization"
#   * marker sets     read back out of Processed_Protein_Ranks.xlsx$MarkerType
#                     plus the two sets refreshed in the 2026-04-14 plot
#   * palette         recovered from the SVG fill attributes
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))
suppressPackageStartupMessages({ library(ggrepel); library(scales); library(readxl) })

source(file.path(REPO_ROOT, "R", "analysis_labels.R"))
source(file.path(REPO_ROOT, "R", "protigy_input_utils.R"))

STAGE <- ensure_dir(file.path(FULL, "qc"))

# ---- marker sets, recovered from the original outputs ------------------------
# Activation      = mitochondrial ribosome / OXPHOS / import  (Processed_Protein_Ranks.xlsx)
# Stress Response = growth / plasticity signalling            (Processed_Protein_Ranks.xlsx)
# General Neuron and Neuropil/Structure as labelled in the 2026-04-14 SVG and
# still carried by 03_qc_exploration/02_rank_abundance_by_sample_class.r
MARKER_SETS <- list(
  `General Neuron`     = c("Rps9", "Rpl22", "Rps16", "Brd4", "Rpl35", "H1-1", "H1-2", "H1-3", "Rps18", "Rpl6"),
  `Neuropil/Structure` = c("Vcan", "Kcna1", "Mog", "Cntnap1", "Cnp"),
  `Stress Response`    = c("Aktip", "Naa10", "Mycbp", "Ikbkb", "Acvr1b", "Strap"),
  `Activation`         = c("Mrpl53", "Timm9", "Ndufb7", "Yars2", "Ndufs6", "Mrpl50", "Mrpl4")
)
# fill colours read directly out of Rank-Abundance_Marker_QC.svg
MARKER_COLORS <- c(
  `General Neuron`     = "#7F8C8D",
  `Activation`         = "#D35400",
  `Stress Response`    = "#C71585",
  `Neuropil/Structure` = "#2980B9",
  `None`               = "#D9D9D9"
)
BACKGROUND_COLOR <- "#D9D9D9"

marker_tbl <- do.call(rbind, lapply(names(MARKER_SETS), function(k)
  data.frame(Genes = MARKER_SETS[[k]], MarkerType = k, stringsAsFactors = FALSE)))

# ---- animal-level abundance -> group means -> ranks --------------------------
parsed <- validate_protigy_gct_v13(GCT_ANIMAL)
mat  <- parsed$matrix
cmet <- parsed$column_metadata
stopifnot(identical(colnames(cmet), colnames(mat)))

sample_class <- normalize_sample_class(cmet["sample_class", ])
condition    <- normalize_condition(cmet["condition", ])
animal_id    <- trimws(as.character(cmet["AnimalID", ]))
grp <- paste(sample_class, gsub("_", "-", condition), sep = "_")

# protein id -> gene symbol, taken from the canonical MapThatProt output
idmap <- utils::read.csv(file.path(MAPPED_DIR, "mcherry_paired_veh_vs_mcherry_unpaired_veh.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
idmap <- unique(idmap[, c("original_protein_id", "mapped_gene_symbol")])
gene_of <- setNames(idmap$mapped_gene_symbol, idmap$original_protein_id)

groups <- sort(unique(grp))
rank_data <- do.call(rbind, lapply(groups, function(g) {
  cols <- which(grp == g)
  mlog <- rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)   # mean over the 3 ANIMALS
  gn   <- unname(gene_of[rownames(mat)])
  keep <- !is.na(gn) & nzchar(gn) & is.finite(mlog)
  d <- data.frame(Condition = g, Genes = gn[keep], MeanLog2 = mlog[keep],
                  stringsAsFactors = FALSE)
  # one row per gene symbol (largest mean, matching the original best-instance rule)
  d <- d[order(-d$MeanLog2), , drop = FALSE]
  d <- d[!duplicated(d$Genes), , drop = FALSE]
  d$LinearValue <- 2 ^ d$MeanLog2
  d <- d[order(-d$LinearValue), , drop = FALSE]
  d$Rank <- seq_len(nrow(d))
  d$n_animals <- length(cols)
  d$animals <- paste(sort(unique(animal_id[cols])), collapse = ";")
  d
}))
rank_data <- merge(rank_data, marker_tbl, by = "Genes", all.x = TRUE, sort = FALSE)
rank_data$MarkerType[is.na(rank_data$MarkerType)] <- "None"
rank_data <- rank_data[order(rank_data$Condition, rank_data$Rank), , drop = FALSE]

cat("groups:", length(groups), " proteins per group:", nrow(rank_data) / length(groups), "\n")

# ---- plotting ----------------------------------------------------------------
theme_rank <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(colour = "#1A1A1A"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.line = element_line(colour = "#1A1A1A", linewidth = 0.35),
      axis.ticks = element_line(colour = "#1A1A1A", linewidth = 0.35),
      axis.text = element_text(colour = "#1A1A1A", size = rel(0.8)),
      axis.title = element_text(colour = "#1A1A1A", size = rel(1.0), face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = rel(0.95)),
      plot.title = element_text(face = "bold", size = rel(1.25)),
      plot.subtitle = element_text(colour = "#555555", size = rel(0.85)),
      legend.position = "right",
      legend.key = element_blank(),
      panel.spacing = unit(1, "lines")
    )
}

rank_plot <- function(d, ncol_facets, show_title = TRUE) {
  lab <- d[d$MarkerType != "None", , drop = FALSE]
  ggplot(d, aes(x = Rank, y = LinearValue)) +
    geom_point(colour = BACKGROUND_COLOR, alpha = 0.6, size = 0.25, shape = 16) +
    geom_line(colour = BACKGROUND_COLOR, linewidth = 0.3) +
    geom_point(data = lab, aes(colour = MarkerType), size = 1.5) +
    ggrepel::geom_label_repel(
      data = lab, aes(label = Genes, fill = MarkerType),
      colour = "white", size = 2.2, fontface = "bold", family = "sans",
      label.padding = unit(0.15, "lines"), label.r = unit(0, "lines"), label.size = 0,
      segment.colour = "#7F7F7F", segment.linewidth = 0.25,
      box.padding = 0.35, point.padding = 0.15, force = 10, max.overlaps = Inf
    ) +
    scale_colour_manual(values = MARKER_COLORS, name = "Biological Category",
                        breaks = names(MARKER_SETS)) +
    scale_fill_manual(values = MARKER_COLORS, guide = "none") +
    # log10 axis. The animal-level values are 2^(log2 ratio) and span barely one
    # decade, so 10^x exponent labels would all be fractional; plain numbers on a
    # log10 axis are used instead. Axis title wording is the original.
    scale_y_log10(expand = expansion(mult = c(0.05, 0.22)),
                  breaks = scales::breaks_log(n = 5, base = 10),
                  labels = scales::label_number(accuracy = 0.1)) +
    scale_x_continuous(expand = expansion(mult = c(0.025, 0.02)), labels = scales::label_comma()) +
    facet_wrap(~ Condition, ncol = ncol_facets, scales = "fixed") +
    labs(
      x = "Protein Rank (ordered by abundance)",
      y = "Abundance (Intensity, Log10 Scale)",
      title = if (show_title) "Rank-Abundance Distribution and Marker Validation" else NULL,
      subtitle = if (show_title) "Log2-transformed input data re-linearized for Log10 visualization" else NULL
    ) +
    theme_rank()
}

# ---- full 16-facet stage output ---------------------------------------------
p_all <- rank_plot(rank_data, ncol_facets = 4)
ggsave(file.path(STAGE, "Rank-Abundance_Marker_QC_animal_level.svg"), p_all,
       width = 16, height = 15, device = svglite::svglite)
utils::write.csv(rank_data, file.path(STAGE, "Processed_Protein_Ranks_animal_level.csv"), row.names = FALSE)

# ---- Figure 3 E : mCherry paired-veh (left) and unpaired-veh (right) --------
e_groups <- c("mcherry_paired-veh", "mcherry_unpaired-veh")
de <- rank_data[rank_data$Condition %in% e_groups, , drop = FALSE]
de$Condition <- factor(de$Condition, levels = e_groups)
p_e <- rank_plot(de, ncol_facets = 2)
save_panel(p_e, STAGE, PANELS_F3, "Fig3E_rank_abundance_animal_level.svg", width = 8.6, height = 4.6)
save_source_data(de, STAGE, PANELS_F3, "Fig3E_rank_abundance_animal_level_source_data.csv")

# ---- Supplementary D : neuron and neuropil, unpaired-veh --------------------
d_groups <- c("neuron_unpaired-veh", "neuropil_unpaired-veh")
dd <- rank_data[rank_data$Condition %in% d_groups, , drop = FALSE]
dd$Condition <- factor(dd$Condition, levels = d_groups)
p_d <- rank_plot(dd, ncol_facets = 2)
save_panel(p_d, STAGE, PANELS_SUP, "SuppD_rank_abundance_animal_level.svg", width = 8.6, height = 4.6)
save_source_data(dd, STAGE, PANELS_SUP, "SuppD_rank_abundance_animal_level_source_data.csv")

# ---- comparison against the original hemisphere-level ranks -----------------
orig_path <- file.path(DATA_ROOT, "03_output/qc/Processed_Protein_Ranks.xlsx")
cmp <- NULL
if (file.exists(orig_path)) {
  cmp <- do.call(rbind, lapply(c(e_groups, d_groups), function(g) {
    o <- as.data.frame(readxl::read_excel(orig_path, sheet = g))
    n <- rank_data[rank_data$Condition == g, c("Genes", "Rank", "MeanLog2")]
    names(n) <- c("Genes", "Rank_animal_level", "MeanLog2_animal_level")
    m <- merge(o[, c("Genes", "Rank", "MeanLog2")], n, by = "Genes")
    data.frame(
      group = g,
      n_shared_genes = nrow(m),
      n_genes_original = nrow(o),
      n_genes_animal_level = nrow(n),
      spearman_rank = stats::cor(m$Rank, m$Rank_animal_level, method = "spearman"),
      pearson_meanlog2 = stats::cor(m$MeanLog2, m$MeanLog2_animal_level),
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(cmp, file.path(STAGE, "rank_abundance_original_vs_animal_level.csv"), row.names = FALSE)
  print(cmp)
}

audit <- data.frame(
  item = c("panel", "original_output", "original_source_data", "original_generating_script",
           "original_sampling_unit", "corrected_input", "corrected_sampling_unit",
           "analysis_type", "n_groups", "n_proteins_per_group", "marker_categories", "note"),
  value = c("Figure3 E and Supplementary D",
            "03_output/qc/Rank-Abundance_Marker_QC.svg",
            "03_output/qc/Processed_Protein_Ranks.xlsx",
            "SOURCE_NOT_FOUND (reconstructed from the outputs above)",
            "hemisphere/acquisition means within sample_class x condition",
            "02_data/animal_level/input_gct/neha_protigy_input_animal_level_primary.gct",
            "animal (L/R averaged first, then mean of 3 animals per group)",
            "descriptive abundance ranking; no inferential statistic is displayed",
            as.character(length(groups)),
            as.character(nrow(rank_data) / length(groups)),
            paste(names(MARKER_SETS), collapse = "; "),
            "Marker sets and palette were recovered from the original artefacts and are preserved verbatim."),
  stringsAsFactors = FALSE
)
utils::write.csv(audit, file.path(STAGE, "Fig3E_SuppD_rank_abundance_audit.csv"), row.names = FALSE)

message("Panel E / Supp D done.")
