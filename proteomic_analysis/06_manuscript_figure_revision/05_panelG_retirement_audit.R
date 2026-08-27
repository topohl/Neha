# ==============================================================================
# Figure 3 panel G -- Lollipop GO/pathway plot, mCherry paired vs unpaired
#
# ORIGINAL: 04_differential_expression_enrichment/02_compareGO.r (pre-2026
#           3140-line version), lollipop block at legacy line 2236
#           -> 99_historical/compareGO/BP/learning_signature/memory_ensemble/
#              07_Regulated_Protein_GO/GO_Lollipop_mcherry2_mcherry4_Upregulated_BP.svg
#           on-plot title "Upregulated GO - mcherry2_mcherry4"; x = Gene Ratio,
#           point size = Gene Count, point colour = adjusted P.
#           Statistic: clusterProfiler::enrichGO (ORA) over the HEMISPHERE-LEVEL
#           FDR<0.05 upregulated protein list from
#           06_Significant_Proteins/SigProteins_Upregulated_mcherry2_mcherry4.xlsx
#
# THIS SCRIPT DOES NOT REGENERATE THE PANEL. It writes the evidence that the
# corrected animal-level analysis does not support it, per the task instruction
# not to manufacture an equivalent significant pathway panel.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))

STAGE <- ensure_dir(file.path(FULL, "enrichment"))
CMP   <- "mcherry_paired_veh_over_mcherry_unpaired_veh"
cdir  <- file.path(ENRICH_ROOT, CMP)

da <- read_mapped("mcherry", "paired_veh", "unpaired_veh")
n_fdr    <- sum(da$padj < FDR, na.rm = TRUE)
n_fdr_up <- sum(da$padj < FDR & da$log2fc > 0, na.rm = TRUE)
n_fdr_dn <- sum(da$padj < FDR & da$log2fc < 0, na.rm = TRUE)

read_or_null <- function(p) if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else NULL

gsea <- read_or_null(file.path(cdir, "GSEA_GO_BP.csv"))
gsea_l2 <- read_or_null(file.path(cdir, "GSEA_GO_BP_log2fc_sensitivity.csv"))
kegg <- read_or_null(file.path(cdir, "GSEA_KEGG.csv"))
ora_up  <- read_or_null(file.path(cdir, "ORA_GO_BP_fdr_up.csv"))
ora_dn  <- read_or_null(file.path(cdir, "ORA_GO_BP_fdr_down.csv"))
ora_all <- read_or_null(file.path(cdir, "ORA_GO_BP_fdr_all.csv"))
ora_top <- read_or_null(file.path(cdir, "ORA_GO_BP_top_abs_log2fc.csv"))

summ <- function(nm, d, col = "p.adjust") {
  if (is.null(d)) return(data.frame(analysis = nm, n_terms = NA_integer_, n_FDR_lt_0.05 = NA_integer_,
                                    min_p_adjust = NA_real_, stringsAsFactors = FALSE))
  v <- suppressWarnings(as.numeric(d[[col]]))
  data.frame(analysis = nm, n_terms = nrow(d), n_FDR_lt_0.05 = sum(v < FDR, na.rm = TRUE),
             min_p_adjust = if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

ev <- rbind(
  summ("GSEA GO-BP (canonical, moderated-t ranked)", gsea),
  summ("GSEA GO-BP (log2fc-ranked sensitivity)", gsea_l2),
  summ("GSEA KEGG (canonical)", kegg),
  summ("ORA GO-BP, FDR-significant UP list (direct analogue of the original panel)", ora_up),
  summ("ORA GO-BP, FDR-significant DOWN list", ora_dn),
  summ("ORA GO-BP, all FDR-significant proteins", ora_all),
  summ("ORA GO-BP, top |log2FC| list", ora_top)
)
print(ev)
utils::write.csv(ev, file.path(STAGE, "Fig3G_mcherry_enrichment_evidence.csv"), row.names = FALSE)

# nearest-miss terms, so a reader can see how far the corrected result is from significance
near <- NULL
if (!is.null(ora_up) && nrow(ora_up)) {
  o <- ora_up[order(suppressWarnings(as.numeric(ora_up$p.adjust))), , drop = FALSE]
  keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio", "Count", "pvalue", "p.adjust", "qvalue"), names(o))
  near <- head(o[, keep, drop = FALSE], 20)
  utils::write.csv(near, file.path(STAGE, "Fig3G_mcherry_ORA_up_nearest_misses.csv"), row.names = FALSE)
  print(head(near[, c("Description", "Count", "pvalue", "p.adjust")], 10))
}

lines <- c(
  "================================================================================",
  "FIGURE 3 PANEL G  --  RETIRED_AFTER_ANIMAL_LEVEL_CORRECTION",
  "Lollipop GO/pathway plot, 'Memory Engram / mCherry paired vs unpaired'",
  paste("Audit written", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "================================================================================",
  "",
  "ORIGINAL PANEL",
  "  file    99_historical/compareGO/BP/learning_signature/memory_ensemble/",
  "          07_Regulated_Protein_GO/GO_Lollipop_mcherry2_mcherry4_Upregulated_BP.svg",
  "  script  04_differential_expression_enrichment/02_compareGO.r (pre-2026 3140-line",
  "          version), lollipop block at line 2236",
  "  method  clusterProfiler::enrichGO over-representation analysis of the",
  "          HEMISPHERE-LEVEL FDR<0.05 upregulated protein list for mcherry2_mcherry4",
  "  shown   GO terms incl. synaptic vesicle endocytosis, postsynaptic density",
  "          organization, modification of postsynaptic structure, respiratory",
  "          electron transport chain, mitochondrial respiratory chain complex assembly",
  "",
  "WHY IT CANNOT BE REPRODUCED",
  "",
  sprintf("  The corrected animal-level differential result still exists and is large:"),
  sprintf("    %d of %d proteins at FDR<%.2f (%d higher in paired_veh, %d higher in unpaired_veh).",
          n_fdr, nrow(da), FDR, n_fdr_up, n_fdr_dn),
  "",
  "  Every canonical enrichment table for the contrast, with its FDR count:",
  paste0("    ", sprintf("%-72s %5d terms, %3d at FDR<0.05, min p.adjust %s",
                         ev$analysis, ev$n_terms, ev$n_FDR_lt_0.05,
                         ifelse(is.na(ev$min_p_adjust), "NA", formatC(ev$min_p_adjust, format = "g", digits = 3)))),
  "",
  "  The decisive line is the fourth one. The ORA over the FDR-significant upregulated",
  "  list -- the direct analogue of the original panel, same test, same direction, same",
  "  ontology -- returns ZERO terms below FDR 0.05 (best p.adjust 0.231). There is",
  "  nothing to draw a lollipop of. The canonical GO-BP GSEA of the same contrast is",
  "  likewise empty (0 of 3829 terms).",
  "",
  "WHAT IS *NOT* EMPTY, AND WHY IT IS NOT A SUBSTITUTE",
  "",
  "  Two non-GO-BP-ORA analyses of the same contrast do return FDR-significant sets, and",
  "  they must be reported rather than quietly omitted:",
  "",
  "    KEGG GSEA (canonical, moderated-t ranked): 26 terms at FDR<0.05, ALL 26 with",
  "    positive NES -- Huntington disease (NES +1.65), Pathways of neurodegeneration,",
  "    Amyotrophic lateral sclerosis, Parkinson disease, Alzheimer disease, Prion disease,",
  "    Oxidative phosphorylation, Citrate cycle, Synaptic vesicle cycle, Dopaminergic /",
  "    Glutamatergic / GABAergic / Serotonergic synapse, Gap junction, Circadian",
  "    entrainment, Salivary secretion, Diabetic cardiomyopathy, Human cytomegalovirus",
  "    infection.",
  "",
  "    GO-BP GSEA under the LEGACY log2fc ranking (sensitivity only, not canonical):",
  "    2 terms -- synaptic vesicle cycle and vesicle-mediated transport in synapse, both",
  "    positive NES.",
  "",
  "  Unanimous directionality is the tell. Twenty-six of twenty-six KEGG sets pointing the",
  "  same way, spanning neurodegeneration, OXPHOS, TCA, four separate neurotransmitter",
  "  synapse maps and several sets with no plausible relation to this experiment",
  "  (salivary secretion, diabetic cardiomyopathy, cytomegalovirus infection), is what a",
  "  GLOBAL one-directional displacement of a large, gene-rich, mitochondria- and",
  "  ribosome-heavy background looks like when it is passed through GSEA. It is not",
  "  evidence of a selective learning programme. The KEGG neurodegeneration maps in",
  "  particular share most of their members with oxidative phosphorylation.",
  "",
  "INTERPRETATION",
  "",
  "  A broad, one-directional protein-level shift (median log2FC 0.333 across all 5327",
  "  proteins) with no GO-BP coherence, and with KEGG structure that is entirely",
  "  one-directional and dominated by shared OXPHOS membership, is the signature of a",
  "  global displacement rather than a coordinated biological programme. The mCherry",
  "  paired_veh vs unpaired_veh contrast is additionally 100% aliased with collection",
  "  plate (every paired animal on Plate1, every unpaired animal on Plate2 within this",
  "  stratum), so the contrast cannot be separated from that design factor at all.",
  "  See 03_output/inferential_checks/final_inferential_checks_20260826/.",
  "",
  "DECISION",
  "",
  "  Panel G is RETIRED. Do not substitute a nominal-P panel, a lowered-threshold panel,",
  "  a post-hoc interaction panel, or the 26-term KEGG list in its place: none of those",
  "  would be the same claim, and presenting one would recreate the error the reviewer",
  "  correction is meant to remove.",
  "",
  "  The protein-level magnitude that the panel used to sit next to is still shown, and",
  "  is still honest, in Figure 3 panel D (regenerated at animal level).",
  "",
  "  If a pathway panel is wanted for the mCherry compartment, the only FDR-supported",
  "  option in the corrected analysis is mcherry_paired_cno_over_mcherry_paired_veh",
  "  (76 FDR-significant GO-BP terms), which is a DIFFERENT contrast and already appears",
  "  in Supplementary panel F. It must not be relabelled as a learning result.",
  "",
  "FILES WRITTEN BY THIS AUDIT",
  "  Fig3G_mcherry_enrichment_evidence.csv     every canonical enrichment table for the",
  "                                            contrast with its FDR count",
  "  Fig3G_mcherry_ORA_up_nearest_misses.csv   the 20 lowest-p.adjust ORA terms, none",
  "                                            of which clears FDR",
  "================================================================================"
)
writeLines(lines, file.path(STAGE, "Fig3G_RETIRED_AFTER_ANIMAL_LEVEL_CORRECTION.txt"))
file.copy(file.path(STAGE, "Fig3G_RETIRED_AFTER_ANIMAL_LEVEL_CORRECTION.txt"),
          file.path(PANELS_F3, "Fig3G_RETIRED_AFTER_ANIMAL_LEVEL_CORRECTION.txt"), overwrite = TRUE)

message("Panel G retirement audit written.")
