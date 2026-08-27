# ==============================================================================
# VISUAL-FIDELITY PASS -- Figure 3 E and Supplementary D, rank abundance
#
# No statistics recomputed. Ranks and mean abundances are read from the
# source-data CSVs written by scripts/03_panelE_suppD_rank_abundance.R.
#
# Fixes vs the first pass:
#   * plot title and subtitle removed (the manuscript panels carry neither)
#   * facet strip labels removed (the manuscript sets those in Illustrator)
#   * legend moved inside the plot at bottom left, legend title removed,
#     keys drawn as plain squares
#   * marker categories split the way the manuscript splits them and renamed to
#     the manuscript wording:
#        Fig 3E  -> Mitochondrial, Growth-Related Plasticity
#        Supp D  -> Neuron Soma, Neuropil/Structure
#     The hex values are the ones recovered from
#     03_output/qc/Rank-Abundance_Marker_QC.svg; only the category NAMES are
#     reconstructed from the manuscript legend text.
#   * y axis title "log10 Intensity"; x axis title only on Supplementary D,
#     matching the manuscript.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "10_style_match_common.R"))
suppressPackageStartupMessages({ library(ggrepel); library(scales) })

FIRST_F3  <- file.path(OUT_ROOT, "figure_panels", "Figure3")
FIRST_SUP <- file.path(OUT_ROOT, "figure_panels", "Supplementary_proteomics")

marker_tbl <- do.call(rbind, lapply(names(RANK_SETS), function(k)
  data.frame(Genes = RANK_SETS[[k]], Category = k, stringsAsFactors = FALSE)))

prep <- function(csv, keep_cats, group_levels) {
  d <- utils::read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
  d$Category <- marker_tbl$Category[match(d$Genes, marker_tbl$Genes)]
  d$Category[!d$Category %in% keep_cats] <- NA_character_
  d$Condition <- factor(d$Condition, levels = group_levels)
  d
}

rank_panel <- function(d, keep_cats, x_title) {
  lab <- d[!is.na(d$Category), , drop = FALSE]
  lab$Category <- factor(lab$Category, levels = keep_cats)
  pal <- RANK_COLS[keep_cats]

  ggplot(d, aes(x = Rank, y = LinearValue)) +
    geom_point(colour = RANK_BG, alpha = 0.6, size = 0.25, shape = 16) +
    geom_line(colour = RANK_BG, linewidth = 0.3) +
    geom_point(data = lab, aes(colour = Category), size = 1.5) +
    ggrepel::geom_label_repel(
      data = lab, aes(label = Genes, fill = Category),
      colour = "white", size = 2.2, fontface = "bold", family = FONT,
      label.padding = unit(0.15, "lines"), label.r = unit(0, "lines"), label.size = 0,
      segment.colour = "#7F7F7F", box.padding = 0.35, point.padding = 0.15,
      force = 10, max.overlaps = Inf, show.legend = FALSE
    ) +
    # square keys inside the panel, as in the manuscript legend
    scale_colour_manual(values = pal, name = NULL,
                        guide = guide_legend(override.aes = list(shape = 15, size = 2.6))) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.22)),
                  breaks = scales::breaks_log(n = 5, base = 10),
                  labels = scales::label_number(accuracy = 0.1)) +
    scale_x_continuous(expand = expansion(mult = c(0.025, 0.02)), labels = scales::label_comma()) +
    facet_wrap(~ Condition, ncol = 2, scales = "fixed") +
    labs(x = x_title, y = expression(log[10] ~ Intensity)) +
    theme_classic(base_size = 9, base_family = FONT) +
    theme(
      text = element_text(colour = "#1A1A1A", family = FONT),
      axis.line  = element_line(colour = "#1A1A1A", linewidth = 0.35 * PT * 2),
      axis.ticks = element_line(colour = "#1A1A1A", linewidth = 0.35 * PT * 2),
      axis.text  = element_text(colour = "#1A1A1A", size = 7),
      axis.title = element_text(colour = "#1A1A1A", size = 9),
      strip.background = element_blank(),
      strip.text = element_blank(),                 # manuscript sets panel titles in Illustrator
      panel.grid = element_blank(),
      panel.spacing = unit(1.2, "lines"),
      legend.position = c(0.015, 0.015),
      legend.justification = c(0, 0),
      legend.direction = "vertical",
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(3.2, "mm"),
      legend.text = element_text(size = 7, colour = "#1A1A1A"),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 6, 4, 4)
    )
}

# ---- Figure 3 E --------------------------------------------------------------
e_cats <- c("Mitochondrial", "Growth-Related Plasticity")
de <- prep(file.path(FIRST_F3, "Fig3E_rank_abundance_animal_level_source_data.csv"),
           e_cats, c("mcherry_paired-veh", "mcherry_unpaired-veh"))
cat("Fig3E labelled markers:", sum(!is.na(de$Category)), "\n")
pe <- rank_panel(de, e_cats, x_title = NULL)
save_sm(pe, SM_F3, "Fig3E_rank_abundance_animal_level_STYLE_MATCHED", width = 8.0, height = 4.0)
copy_source_data(file.path(FIRST_F3, "Fig3E_rank_abundance_animal_level_source_data.csv"),
                 SM_F3, "Fig3E_rank_abundance_animal_level_STYLE_MATCHED_source_data.csv")

# ---- Supplementary D ---------------------------------------------------------
d_cats <- c("Neuron Soma", "Neuropil/Structure")
dd <- prep(file.path(FIRST_SUP, "SuppD_rank_abundance_animal_level_source_data.csv"),
           d_cats, c("neuron_unpaired-veh", "neuropil_unpaired-veh"))
cat("SuppD labelled markers:", sum(!is.na(dd$Category)), "\n")
pd <- rank_panel(dd, d_cats, x_title = "Protein Rank (ordered by abundance)")
save_sm(pd, SM_SUP, "SuppD_rank_abundance_animal_level_STYLE_MATCHED", width = 8.0, height = 4.0)
copy_source_data(file.path(FIRST_SUP, "SuppD_rank_abundance_animal_level_source_data.csv"),
                 SM_SUP, "SuppD_rank_abundance_animal_level_STYLE_MATCHED_source_data.csv")

message("Rank-abundance panels style-matched.")
