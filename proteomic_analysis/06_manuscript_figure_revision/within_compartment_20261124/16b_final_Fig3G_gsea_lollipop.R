# ==============================================================================
# FIGURE 3 G (manuscript Fig. 2G) -- mCherry paired_veh vs unpaired_veh,
# GO GSEA lollipop, per-compartment limma model (2026-11-24).
#
# Replaces the ORA lollipop (16_final_Fig3G_lollipop.R). Against the measured
# proteome no GO-BP term is over-represented among the FDR-significant proteins
# (minimum adjusted P 0.18), so the panel now shows the ranked (GSEA) result,
# which uses every measured protein.
#
# Input   enrichment_within_compartment_20261124/GO_BP and GO_CC,
#         per_comparison/<CMP>/GSEA_GO_BP.csv and GSEA_GO_CC.csv (read as-is).
# Terms   per ontology: adjusted P < 0.05, ranked by adjusted P (then P, then ID);
#         a term is skipped when its leading-edge proteins overlap those of a
#         term already kept by Jaccard >= 0.7 (terms driven by largely the same
#         proteins are shown once; leading-edge Jaccard is enrichplot's term
#         similarity for GSEA results); first 8 BP and 7 CC kept (15 rows, as in
#         the old panel).
# Display x = NES; point size = leading-edge proteins; colour = adjusted P on
#         the old panel's ramp (log scale); rows within each ontology ordered by
#         NES, largest on top. Canvas, segment and type as the old panel.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/figure_rebuild_within_compartment_20261124"
source(file.path(OUT_ROOT, "scripts", "10_style_match_common.R"))

CMP      <- "mcherry_paired_veh_over_mcherry_unpaired_veh"
ENR_RUN  <- dirname(dirname(ENRICH_ROOT))            # .../enrichment_within_compartment_20261124
GSEA     <- c(BP = file.path(ENR_RUN, "GO_BP", "per_comparison", CMP, "GSEA_GO_BP.csv"),
              CC = file.path(ENR_RUN, "GO_CC", "per_comparison", CMP, "GSEA_GO_CC.csv"))
N_KEEP   <- c(BP = 8L, CC = 7L)
SIM_CUT  <- 0.7
ONT_LAB  <- c(BP = "Biological process", CC = "Cellular component")
STAGE    <- ensure_dir(file.path(FULL, "enrichment"))

LOLLI_RAMP   <- c("#E64F4F", "#EFA24A", "#DFB74B")
LOLLI_SEG    <- "#BDBDBD"
LOLLI_SEG_LW <- 1.5
LOLLI_SIZE   <- c(2, 6)
LOLLI_W      <- 6
LOLLI_H      <- 5
LOLLI_FS     <- 9

theme_clean <- function(base_size = 9, base_family = FONT) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "#2C2C2C", family = base_family, lineheight = 1.3),
      axis.line.x = element_line(color = "#2C2C2C", linewidth = 0.6 * PT),
      axis.line.y = element_line(color = "#2C2C2C", linewidth = 0.6 * PT),
      axis.ticks = element_line(color = "#2C2C2C", linewidth = 0.5 * PT),
      axis.ticks.length = unit(3 * PT, "mm"),
      axis.text = element_text(color = "#2C2C2C", size = rel(0.95)),
      axis.title = element_text(color = "#2C2C2C", size = rel(1.0)),
      legend.background = element_blank(), legend.key = element_rect(color = NA, fill = NA),
      legend.key.size = unit(10 * PT, "mm"),
      legend.title = element_text(color = "#2C2C2C", size = rel(0.95)),
      legend.text = element_text(color = "#2C2C2C", size = rel(0.9)),
      legend.position = "right", legend.justification = "top",
      panel.grid = element_blank(), panel.border = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      strip.text.y = element_text(color = "#2C2C2C", size = rel(0.95), angle = 0, hjust = 0),
      panel.spacing.y = unit(4, "mm"),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ---- read and select --------------------------------------------------------
pick <- function(ont) {
  x <- utils::read.csv(GSEA[[ont]], stringsAsFactors = FALSE, check.names = FALSE)
  n_tested <- nrow(x)
  sig <- x[is.finite(x$p.adjust) & x$p.adjust < FDR, , drop = FALSE]
  sig <- sig[order(sig$p.adjust, sig$pvalue, sig$ID), , drop = FALSE]
  le <- stats::setNames(strsplit(sig$core_enrichment, "/", fixed = TRUE), sig$ID)
  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  kept <- character(0); skipped <- data.frame()
  for (id in sig$ID) {
    if (length(kept) >= N_KEEP[[ont]]) break
    if (length(kept)) {
      s <- max(vapply(kept, function(k) jaccard(le[[id]], le[[k]]), numeric(1)))
      if (s >= SIM_CUT) {
        skipped <- rbind(skipped, data.frame(ontology = ont, ID = id,
                                             Description = sig$Description[sig$ID == id],
                                             p.adjust = sig$p.adjust[sig$ID == id],
                                             max_similarity_to_kept = s))
        next
      }
    }
    kept <- c(kept, id)
  }
  sel <- sig[match(kept, sig$ID), , drop = FALSE]
  sel$ontology <- ont
  sel$leading_edge_proteins <- lengths(strsplit(sel$core_enrichment, "/", fixed = TRUE))
  list(sel = sel, n_tested = n_tested, n_sig = nrow(sig), skipped = skipped)
}
picks <- lapply(c(BP = "BP", CC = "CC"), pick)
sel <- do.call(rbind, lapply(picks, `[[`, "sel"))
for (o in names(picks)) cat(sprintf("%s: %d tested, %d at FDR < %.2f, %d plotted, %d skipped as redundant\n",
                                    o, picks[[o]]$n_tested, picks[[o]]$n_sig, FDR,
                                    sum(sel$ontology == o), nrow(picks[[o]]$skipped)))

sel$Ontology <- factor(unname(ONT_LAB[sel$ontology]), levels = unname(ONT_LAB))
sel <- sel[order(sel$Ontology, -sel$NES), , drop = FALSE]
sel$DescriptionWrapped <- stringr::str_wrap(sel$Description, width = 50)
sel$DescriptionWrapped <- factor(sel$DescriptionWrapped, levels = rev(unique(sel$DescriptionWrapped)))
print(sel[, c("ontology", "Description", "setSize", "leading_edge_proteins", "NES", "p.adjust")], row.names = FALSE)

# ---- plot --------------------------------------------------------------------
p <- ggplot(sel, aes(x = NES, y = DescriptionWrapped)) +
  geom_segment(aes(x = 0, xend = NES, y = DescriptionWrapped, yend = DescriptionWrapped),
               colour = LOLLI_SEG, linewidth = LOLLI_SEG_LW * PT * 2) +
  geom_point(aes(size = leading_edge_proteins, colour = p.adjust), alpha = 0.85, stroke = 0) +
  facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradientn(colours = grDevices::colorRampPalette(LOLLI_RAMP, space = "Lab")(100),
                         trans = "log10", name = expression(italic(P)[adj])) +
  scale_size_continuous(range = LOLLI_SIZE, name = "Leading-edge\nproteins") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08)), name = "NES") +
  labs(y = NULL) +
  theme_clean(base_size = LOLLI_FS)

save_sm(p, SM_F3, "Fig3G_mCherry_GO_GSEA_lollipop_within_compartment_FINAL", LOLLI_W, LOLLI_H)

# ---- source data and audit ----------------------------------------------------
src <- sel[, c("ontology", "ID", "Description", "setSize", "leading_edge_proteins", "enrichmentScore",
               "NES", "pvalue", "p.adjust", "qvalue", "core_enrichment")]
src$n_terms_tested <- vapply(src$ontology, function(o) picks[[o]]$n_tested, integer(1))
src$n_FDR_significant <- vapply(src$ontology, function(o) picks[[o]]$n_sig, integer(1))
src$selection_rule <- sprintf(paste0("adjusted P < %.2f, ranked by adjusted P; skip a term whose leading-edge Jaccard ",
                                     "with a kept term is >= %.1f; first %d BP and %d CC"),
                              FDR, SIM_CUT, N_KEEP[["BP"]], N_KEEP[["CC"]])
src$input_file <- unname(GSEA[src$ontology])
utils::write.csv(src, file.path(SM_F3, "Fig3G_source_data.csv"), row.names = FALSE)
utils::write.csv(src, file.path(STAGE, "Fig3G_source_data.csv"), row.names = FALSE)
skipped <- do.call(rbind, lapply(picks, `[[`, "skipped"))
utils::write.csv(if (is.null(skipped)) data.frame() else skipped,
                 file.path(SM_F3, "Fig3G_redundant_terms_skipped.csv"), row.names = FALSE)

message("Figure 3 G (GSEA lollipop) written: ", nrow(sel), " terms.")
