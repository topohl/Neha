# ==============================================================================
# Rank-abundance and marker-validation QC, ANIMAL LEVEL
#
# Regenerates the rank-abundance panels used in the corrected proteomics figures:
#   Figure 3E                 mcherry_paired-veh  + mcherry_unpaired-veh
#   Supplementary proteomics D  neuron_unpaired-veh + neuropil_unpaired-veh
#
# ------------------------------------------------------------------------------
# WHY THIS CHANGED (2026-08-28)
#
# The stage used to read four per-sample-class imputed workbooks from
# 02_data/gct/imputed/ and summarise by sample_class alone. That folder is not
# present in the project tree, so the stage was dead -- its output directory had
# never been created -- and run_pipeline_check.ps1 recorded it as SKIP.
#
# It now reads the validated animal-level GCT and summarises by
# sample_class x condition, which is the grouping the manuscript panels use.
# Hemisphere-level observations are never used for a biological group summary.
#
# NOT changed: the ranking definition, the abundance scale, the group-mean
# logic, or any normalisation/filtering/imputation. The GCT values are already
# processed; 2^MeanLog2 is the original re-linearisation for log10 display, not
# a new transformation.
#
# Numerical acceptance reference (verified equal, see the validation audit):
#   03_output/reviewer_revision_animal_level_20260827/full_regenerated/qc/
#     Processed_Protein_Ranks_animal_level.csv          (16 groups x 5310 rows)
#   .../figure_panels/Figure3/Fig3E_rank_abundance_animal_level_source_data.csv
#   .../figure_panels/Supplementary_proteomics/SuppD_rank_abundance_animal_level_source_data.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(scales)
  library(ggrepel)
  library(writexl)
  library(svglite)
})

option_or_env <- function(option_name, env_name, default) {
  value <- getOption(option_name)
  if (!is.null(value) && nzchar(trimws(as.character(value)))) return(as.character(value))
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(trimws(value))) return(value)
  default
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else
  file.path("03_qc_exploration", "02_rank_abundance_by_sample_class.r")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "protigy_input_utils.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "protigy_input_utils.R"))

project_root <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"

# Validated animal-level matrix: 5349 protein groups x 48 animal-level samples,
# 4 sample classes x 4 conditions, n = 3 animals per cell.
input_gct <- option_or_env(
  "neha.rank_abundance_input_gct", "NEHA_RANK_ABUNDANCE_INPUT_GCT",
  file.path(project_root, "02_data", "animal_level", "input_gct",
            "neha_protigy_input_animal_level_primary.gct")
)

# The GCT is keyed by UniProt entry name (e.g. NOVA2_MOUSE); the marker sets and
# the plotted labels are gene symbols. This is the canonical MapThatProt output,
# used ONLY as a deterministic original_protein_id -> mapped_gene_symbol lookup.
# Every one of the 12 mapped contrast files carries the same lookup; the forward
# mcherry learning contrast is named here for definiteness. No fuzzy matching and
# no fallback: an unmapped protein group is dropped, exactly as in the reference.
id_map_path <- option_or_env(
  "neha.rank_abundance_id_map", "NEHA_RANK_ABUNDANCE_ID_MAP",
  file.path(project_root, "02_data", "animal_level", "mapped", "forward",
            "mcherry_paired_veh_vs_mcherry_unpaired_veh.csv")
)

saving_dir <- option_or_env(
  "neha.rank_abundance_output_dir", "NEHA_RANK_ABUNDANCE_OUTPUT_DIR",
  file.path(project_root, "03_output", "qc", "rank_abundance")
)

if (!file.exists(input_gct)) {
  stop(
    "Animal-level GCT not found: ", input_gct,
    "\nSet NEHA_RANK_ABUNDANCE_INPUT_GCT / options(neha.rank_abundance_input_gct=) to its location.",
    call. = FALSE
  )
}
if (!file.exists(id_map_path)) {
  stop(
    "Protein-id to gene-symbol mapping not found: ", id_map_path,
    "\nSet NEHA_RANK_ABUNDANCE_ID_MAP / options(neha.rank_abundance_id_map=) to a canonical",
    "\nMapThatProt output carrying original_protein_id and mapped_gene_symbol.",
    call. = FALSE
  )
}
dir.create(saving_dir, recursive = TRUE, showWarnings = FALSE)

sample_classes <- c("mcherry", "neuropil", "cfos", "neuron")
conditions <- c("paired_cno", "paired_veh", "unpaired_cno", "unpaired_veh")

# Panel definitions, recovered from the validated reviewer-revision source data
# rather than from the Illustrator titles.
panel_groups <- list(
  Fig3E = c("mcherry_paired-veh", "mcherry_unpaired-veh"),
  SuppD = c("neuron_unpaired-veh", "neuropil_unpaired-veh")
)

# Marker categories, derived from the validated reference (MarkerType column of
# Processed_Protein_Ranks_animal_level.csv). Category keys are the reference
# names; the manuscript legend labels are applied at plot time only.
marker_sets <- list(
  `General Neuron`     = c("Brd4", "H1-1", "H1-2", "H1-3", "Rpl22", "Rpl35", "Rpl6", "Rps16", "Rps18", "Rps9"),
  `Neuropil/Structure` = c("Cnp", "Cntnap1", "Kcna1", "Mog", "Vcan"),
  `Stress Response`    = c("Acvr1b", "Aktip", "Ikbkb", "Mycbp", "Naa10", "Strap"),
  `Activation`         = c("Mrpl4", "Mrpl50", "Mrpl53", "Ndufb7", "Ndufs6", "Timm9", "Yars2")
)

# Highlight colours read out of the historical Rank-Abundance_Marker_QC.svg.
marker_colors <- c(
  `General Neuron`     = "#7F8C8D",
  `Neuropil/Structure` = "#2980B9",
  `Stress Response`    = "#C71585",
  `Activation`         = "#D35400"
)
background_color <- "#D9D9D9"

# Manuscript legend wording for the same categories.
marker_display_labels <- c(
  `General Neuron`     = "Neuron Soma",
  `Neuropil/Structure` = "Neuropil/Structure",
  `Stress Response`    = "Growth-Related Plasticity",
  `Activation`         = "Mitochondrial"
)

# Categories shown on each manuscript panel.
panel_categories <- list(
  Fig3E = c("Activation", "Stress Response"),
  SuppD = c("General Neuron", "Neuropil/Structure")
)

sample_class_colors <- c(
  mcherry = "#4C78A8",
  neuropil = "#72B7B2",
  cfos = "#E45756",
  neuron = "#54A24B"
)

# ------------------------------------------------------------------------------
# Input
# ------------------------------------------------------------------------------

read_animal_level_matrix <- function(gct_path) {
  parsed <- validate_protigy_gct_v13(gct_path)
  cm <- parsed$column_metadata
  required <- c("AnimalID", "condition", "sample_class")
  missing <- setdiff(required, rownames(cm))
  if (length(missing)) {
    stop("Animal-level GCT metadata is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  metadata <- data.frame(
    Sample = colnames(parsed$matrix),
    AnimalID = trimws(as.character(cm["AnimalID", ])),
    condition = normalize_condition(cm["condition", ]),
    sample_class = normalize_sample_class(cm["sample_class", ]),
    stringsAsFactors = FALSE
  )
  metadata$group <- paste(metadata$sample_class, gsub("_", "-", metadata$condition), sep = "_")
  list(matrix = parsed$matrix, metadata = metadata)
}

read_id_map <- function(path) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("original_protein_id", "mapped_gene_symbol")
  missing <- setdiff(required, names(d))
  if (length(missing)) {
    stop("Mapping file is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  d <- unique(d[, required, drop = FALSE])
  if (anyDuplicated(d$original_protein_id)) {
    stop("Mapping file assigns more than one gene symbol to a protein id; refusing to guess.", call. = FALSE)
  }
  stats::setNames(d$mapped_gene_symbol, d$original_protein_id)
}

# ------------------------------------------------------------------------------
# Group summaries at the animal level
# ------------------------------------------------------------------------------

make_rank_data <- function(expression_matrix, metadata, gene_of) {
  gene_symbol <- unname(gene_of[rownames(expression_matrix)])
  groups <- sort(unique(metadata$group))
  do.call(rbind, lapply(groups, function(g) {
    cols <- which(metadata$group == g)
    mean_log2 <- rowMeans(expression_matrix[, cols, drop = FALSE], na.rm = TRUE)
    keep <- !is.na(gene_symbol) & nzchar(gene_symbol) & is.finite(mean_log2)
    d <- data.frame(
      Genes = gene_symbol[keep],
      Condition = g,
      MeanLog2 = mean_log2[keep],
      stringsAsFactors = FALSE
    )
    # one row per gene symbol: keep the most abundant protein group for that symbol
    d <- d[order(-d$MeanLog2), , drop = FALSE]
    d <- d[!duplicated(d$Genes), , drop = FALSE]
    d$LinearValue <- 2 ^ d$MeanLog2
    d <- d[order(-d$LinearValue), , drop = FALSE]
    d$Rank <- seq_len(nrow(d))
    d$n_animals <- length(cols)
    d$animals <- paste(sort(unique(metadata$AnimalID[cols])), collapse = ";")
    d
  }))
}

annotate_markers <- function(rank_data) {
  marker_tbl <- purrr::imap_dfr(marker_sets, function(markers, marker_set) {
    tibble::tibble(Genes = markers, MarkerType = marker_set)
  }) %>% dplyr::distinct()
  out <- merge(rank_data, marker_tbl, by = "Genes", all.x = TRUE, sort = FALSE)
  out$MarkerType[is.na(out$MarkerType)] <- "None"
  out[order(out$Condition, out$Rank), , drop = FALSE]
}

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

theme_rank <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "#1A1A1A", family = "sans"),
      axis.line = ggplot2::element_line(colour = "#1A1A1A", linewidth = 0.25),
      axis.ticks = ggplot2::element_line(colour = "#1A1A1A", linewidth = 0.25),
      axis.text = ggplot2::element_text(colour = "#1A1A1A", size = 7),
      axis.title = ggplot2::element_text(colour = "#1A1A1A", size = 9),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.spacing = grid::unit(1.2, "lines"),
      legend.position = c(0.015, 0.015),
      legend.justification = c(0, 0),
      legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.size = grid::unit(3.2, "mm"),
      legend.text = ggplot2::element_text(size = 7, colour = "#1A1A1A"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(4, 6, 4, 4)
    )
}

plot_rank_abundance <- function(plot_data, categories, ncol_facets = 2, x_title = "Protein Rank (ordered by abundance)") {
  labelled <- plot_data[plot_data$MarkerType %in% categories, , drop = FALSE]
  labelled$MarkerType <- factor(labelled$MarkerType, levels = categories)
  pal <- stats::setNames(unname(marker_colors[categories]), categories)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Rank, y = LinearValue)) +
    ggplot2::geom_point(colour = background_color, alpha = 0.6, size = 0.25, shape = 16) +
    ggplot2::geom_line(colour = background_color, linewidth = 0.3) +
    ggplot2::geom_point(data = labelled, ggplot2::aes(colour = MarkerType), size = 1.5) +
    ggrepel::geom_label_repel(
      data = labelled,
      ggplot2::aes(label = Genes, fill = MarkerType),
      colour = "white", size = 2.2, fontface = "bold", family = "sans",
      label.padding = grid::unit(0.15, "lines"), label.r = grid::unit(0, "lines"),
      label.size = 0, segment.colour = "#7F7F7F",
      box.padding = 0.35, point.padding = 0.15, force = 10, max.overlaps = Inf,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = pal, name = NULL, labels = unname(marker_display_labels[categories]),
      guide = ggplot2::guide_legend(override.aes = list(shape = 15, size = 2.6))
    ) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_y_log10(
      expand = ggplot2::expansion(mult = c(0.05, 0.22)),
      breaks = scales::breaks_log(n = 5, base = 10),
      labels = scales::label_number(accuracy = 0.1)
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.025, 0.02)),
      labels = scales::label_comma()
    ) +
    ggplot2::facet_wrap(~ Condition, ncol = ncol_facets, scales = "fixed") +
    ggplot2::labs(x = x_title, y = expression(log[10] ~ Intensity)) +
    theme_rank()
}

save_panel <- function(plot, file_name, width, height) {
  path <- file.path(saving_dir, paste0(file_name, ".svg"))
  ggplot2::ggsave(path, plot, width = width, height = height, units = "in", device = svglite::svglite)
  message("Saved plot: ", path)
  invisible(path)
}

# ==============================================================================
# Run
# ==============================================================================

input <- read_animal_level_matrix(input_gct)
gene_of <- read_id_map(id_map_path)

message(sprintf("Animal-level matrix: %d protein groups x %d samples, %d animals",
                nrow(input$matrix), ncol(input$matrix), length(unique(input$metadata$AnimalID))))

design <- table(input$metadata$sample_class, input$metadata$condition)
if (!all(design == 3L)) {
  stop("Expected exactly 3 animal-level samples per sample_class x condition; got ",
       paste(sort(unique(as.vector(design))), collapse = ", "), call. = FALSE)
}

rank_data <- annotate_markers(make_rank_data(input$matrix, input$metadata, gene_of))
message(sprintf("Rank table: %d groups x %d gene symbols",
                length(unique(rank_data$Condition)), nrow(rank_data) / length(unique(rank_data$Condition))))

utils::write.csv(
  rank_data[, c("Genes", "Condition", "MeanLog2", "LinearValue", "Rank", "n_animals", "animals", "MarkerType")],
  file.path(saving_dir, "processed_protein_ranks_animal_level.csv"),
  row.names = FALSE
)

# ---- one panel per sample_class x condition ----------------------------------
all_categories <- names(marker_sets)
for (g in sort(unique(rank_data$Condition))) {
  d <- rank_data[rank_data$Condition == g, , drop = FALSE]
  save_panel(plot_rank_abundance(d, all_categories, ncol_facets = 1),
             paste0("rank_abundance_", g), width = 4.6, height = 4.6)
}

# ---- manuscript panels -------------------------------------------------------
for (panel in names(panel_groups)) {
  groups <- panel_groups[[panel]]
  d <- rank_data[rank_data$Condition %in% groups, , drop = FALSE]
  d$Condition <- factor(d$Condition, levels = groups)
  x_title <- if (panel == "Fig3E") NULL else "Protein Rank (ordered by abundance)"
  save_panel(plot_rank_abundance(d, panel_categories[[panel]], ncol_facets = 2, x_title = x_title),
             paste0(panel, "_rank_abundance_animal_level"), width = 8.0, height = 4.0)
  utils::write.csv(
    d[, c("Genes", "Condition", "MeanLog2", "LinearValue", "Rank", "n_animals", "animals", "MarkerType")],
    file.path(saving_dir, paste0(panel, "_rank_abundance_animal_level_source_data.csv")),
    row.names = FALSE
  )
}

# ---- marker abundance summary and module scores (retained) -------------------
marker_summary <- rank_data %>%
  dplyr::filter(MarkerType != "None") %>%
  tidyr::separate(Condition, into = c("sample_class", "condition"), sep = "_", remove = FALSE) %>%
  dplyr::group_by(sample_class, Genes, MarkerType) %>%
  dplyr::summarise(
    MeanLog2 = mean(MeanLog2, na.rm = TRUE),
    MedianRank = stats::median(Rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sample_class = factor(sample_class, levels = sample_classes),
    MarkerType = factor(MarkerType, levels = all_categories,
                        labels = unname(marker_display_labels[all_categories]))
  )

marker_summary_plot <- ggplot2::ggplot(marker_summary, ggplot2::aes(x = sample_class, y = MeanLog2)) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, fill = "grey92", color = "black") +
  ggplot2::geom_point(ggplot2::aes(color = sample_class),
                      position = ggplot2::position_jitter(width = 0.12, height = 0), size = 1.6, alpha = 0.9) +
  ggplot2::facet_wrap(~MarkerType, scales = "free_y", nrow = 1) +
  ggplot2::scale_color_manual(values = sample_class_colors) +
  ggplot2::labs(x = NULL, y = "Mean log2 abundance") +
  ggplot2::theme_classic(base_size = 8) +
  ggplot2::theme(
    text = ggplot2::element_text(family = "sans", color = "black"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 7),
    axis.title = ggplot2::element_text(size = 8, face = "bold"),
    axis.text = ggplot2::element_text(size = 6, color = "black"),
    axis.line = ggplot2::element_line(linewidth = 0.3, color = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
    legend.position = "none",
    panel.spacing = grid::unit(1, "lines")
  )

ggplot2::ggsave(file.path(saving_dir, "marker_abundance_summary.svg"), marker_summary_plot,
                width = 180, height = 70, units = "mm", device = svglite::svglite)

# Module scores at the animal level: one score per animal-level sample.
marker_lookup <- purrr::imap_dfr(marker_sets, function(markers, marker_set) {
  tibble::tibble(Genes = markers, MarkerType = marker_set)
})
gene_symbol_all <- unname(gene_of[rownames(input$matrix)])
score_all <- purrr::imap_dfr(marker_sets, function(markers, marker_set) {
  rows <- which(gene_symbol_all %in% markers)
  tibble::tibble(
    Sample = input$metadata$Sample,
    sample_class = input$metadata$sample_class,
    condition = input$metadata$condition,
    Score = colMeans(input$matrix[rows, , drop = FALSE], na.rm = TRUE),
    DetectedMarkers = length(rows),
    MarkerType = marker_set
  )
}) %>%
  dplyr::mutate(
    sample_class = factor(sample_class, levels = sample_classes),
    MarkerType = factor(MarkerType, levels = all_categories,
                        labels = unname(marker_display_labels[all_categories]))
  )

global_stats <- score_all %>%
  dplyr::group_by(MarkerType) %>%
  dplyr::summarise(p_value = stats::kruskal.test(Score ~ sample_class)$p.value, .groups = "drop")

pairwise_stats <- score_all %>%
  dplyr::group_by(MarkerType) %>%
  tidyr::nest() %>%
  dplyr::mutate(pairwise = purrr::map(data, function(df) {
    pw <- stats::pairwise.wilcox.test(df$Score, df$sample_class, p.adjust.method = "holm")
    as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE) %>%
      dplyr::rename(group1 = Var1, group2 = Var2, p_adj = Freq) %>%
      dplyr::filter(!is.na(p_adj))
  })) %>%
  dplyr::select(MarkerType, pairwise) %>%
  tidyr::unnest(pairwise)

module_score_plot <- ggplot2::ggplot(score_all, ggplot2::aes(x = sample_class, y = Score)) +
  ggplot2::geom_boxplot(ggplot2::aes(fill = sample_class), width = 0.48, outlier.shape = NA,
                        linewidth = 0.32, color = "black", alpha = 0.92) +
  ggplot2::geom_point(ggplot2::aes(color = sample_class),
                      position = ggplot2::position_jitter(width = 0.09, height = 0),
                      size = 1.05, alpha = 0.55, stroke = 0) +
  ggplot2::stat_summary(fun = stats::median, geom = "crossbar", width = 0.42, linewidth = 0.28, color = "black") +
  ggplot2::facet_wrap(~MarkerType, scales = "free_y", nrow = 1) +
  ggplot2::scale_fill_manual(values = sample_class_colors) +
  ggplot2::scale_color_manual(values = sample_class_colors) +
  ggplot2::labs(x = NULL, y = "Marker module score\n(mean log2 abundance)") +
  ggplot2::theme_classic(base_size = 7) +
  ggplot2::theme(
    text = ggplot2::element_text(family = "sans", color = "black"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 7.2, margin = ggplot2::margin(b = 6)),
    axis.title.y = ggplot2::element_text(size = 7.2, face = "bold", margin = ggplot2::margin(r = 5)),
    axis.text.x = ggplot2::element_text(size = 6.2, color = "black", angle = 35, hjust = 1, vjust = 1),
    axis.text.y = ggplot2::element_text(size = 6.2, color = "black"),
    axis.line = ggplot2::element_line(linewidth = 0.3, color = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
    axis.ticks.length = grid::unit(1.5, "mm"),
    panel.spacing = grid::unit(1.1, "lines"),
    legend.position = "none",
    plot.margin = ggplot2::margin(4, 4, 4, 4),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA)
  )

ggplot2::ggsave(file.path(saving_dir, "marker_module_scores.svg"), module_score_plot,
                width = 180, height = 75, units = "mm", device = svglite::svglite)

writexl::write_xlsx(
  list(
    marker_abundance_summary = marker_summary,
    marker_module_scores = score_all,
    global_kruskal = global_stats,
    pairwise_wilcox = pairwise_stats,
    rank_marker_rows = rank_data %>%
      dplyr::filter(MarkerType != "None") %>%
      dplyr::arrange(MarkerType, Genes, Condition)
  ),
  path = file.path(saving_dir, "marker_validation_summary.xlsx")
)

message("Animal-level rank-abundance and marker-validation QC completed.")
