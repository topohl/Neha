# ==============================================================================
# Supplementary panel A  -- protein / peptide counts per acquisition  (TECHNICAL QC)
# Supplementary panel B2 -- PCA coloured by Hemisphere                (TECHNICAL QC)
# Supplementary panel B3 -- PCA coloured by Collection plate          (TECHNICAL QC)
#
# B1 (Experimental Group) is regenerated at ANIMAL level in
# scripts/01_panelC_and_suppB_pca.R. B2 and B3 cannot be: hemisphere and
# collection plate are properties of an ACQUISITION and do not exist as columns
# in the 48-row animal-level metadata. They are therefore retained deliberately
# at acquisition level and must be captioned as technical QC, never as
# biological n.
#
# ORIGINAL A : SOURCE_NOT_FOUND. Output
#              S:/Lab_Member/Tobi/Experiments/Collabs/Neha/Results/Plots/protein_peptide_counts.svg
#              data  S:/Lab_Member/Tobi/Experiments/Collabs/Neha/protein_count.xlsx
#              (03_qc_exploration/01_qc_protein_peptide_plot.r is NOT the generator:
#               its input quicksearch.stats.annotated.xlsx has never existed here.)
# ORIGINAL B : 03_qc_exploration/06_pcaPlot_Neha.r (pre-2026-08-26 monolith)
#              plot_and_save_group("ReplicateGroup", ...) legacy line 317
#              plot_and_save_group("plate", ...)          legacy line 666
#              -> 99_historical/pca_plots_legacy/plots/base/pca_by_replicate_group.svg
#              -> 99_historical/pca_plots_legacy/plots/base/pca_by_plate_precorrect.svg
#
# NAMING: "ReplicateGroup" is the hemisphere replicate (1/2) and "plate" is the
# COLLECTION plate. Collection plate is NOT established as a proteomics batch:
# all 96 acquisitions were run on one instrument on one date. It is relabelled
# "Collection plate" and must not be called a batch effect.
# ==============================================================================

OUT_ROOT <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))
suppressPackageStartupMessages({ library(digest) })

STAGE <- ensure_dir(file.path(FULL, "qc"))

NEHA_ROOT       <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha"
ORIG_SUPPA_SVG  <- file.path(NEHA_ROOT, "Results/Plots/protein_peptide_counts.svg")
ORIG_SUPPA_DATA <- file.path(NEHA_ROOT, "protein_count.xlsx")

HEMI_GCT   <- file.path(DATA_ROOT, "02_data/gct/pg.matrix_filtered_pcaAdjusted_unnormalized.gct")
LEGACY_META<- file.path(DATA_ROOT, "99_historical/pca_plots_legacy/tables/meta/sample_metadata_parsed.csv")

# ==============================================================================
# SUPPLEMENTARY A -- retained as technical QC, copied with provenance
# ==============================================================================
suppA <- data.frame(
  item = c("panel", "status", "original_plot", "original_plot_sha256", "original_plot_mtime",
           "original_source_data", "original_source_data_sha256", "original_generating_script",
           "sampling_unit", "why_acquisition_level_is_correct", "required_relabelling", "caption_requirement"),
  value = c(
    "Supplementary A",
    "VALID_AS_TECHNICAL_QC",
    ORIG_SUPPA_SVG,
    if (file.exists(ORIG_SUPPA_SVG)) digest::digest(file = ORIG_SUPPA_SVG, algo = "sha256") else NA_character_,
    if (file.exists(ORIG_SUPPA_SVG)) format(file.info(ORIG_SUPPA_SVG)$mtime, tz = "UTC", usetz = TRUE) else NA_character_,
    ORIG_SUPPA_DATA,
    if (file.exists(ORIG_SUPPA_DATA)) digest::digest(file = ORIG_SUPPA_DATA, algo = "sha256") else NA_character_,
    "SOURCE_NOT_FOUND - no script in the repository or in any git commit references protein_peptide_counts or protein_count.xlsx",
    "acquisition (one bar per injected sample; both hemispheres present)",
    paste("Identification depth is a property of the acquisition, not of the animal. Averaging L and R would",
          "destroy exactly the information the panel exists to show. No biological contrast is drawn from it,",
          "so the animal-level correction does not apply."),
    "Background -> Neuropil; cFosN -> cFos; mCherryN -> mCherry",
    "The caption must state that bars are acquisitions (technical depth) and that they are not the biological n."),
  stringsAsFactors = FALSE
)
utils::write.csv(suppA, file.path(STAGE, "SuppA_protein_peptide_counts_audit.csv"), row.names = FALSE)

if (file.exists(ORIG_SUPPA_SVG)) {
  file.copy(ORIG_SUPPA_SVG, file.path(STAGE, "SuppA_protein_peptide_counts_TECHNICAL_QC.svg"), overwrite = TRUE)
  file.copy(ORIG_SUPPA_SVG, file.path(PANELS_SUP, "SuppA_protein_peptide_counts_TECHNICAL_QC.svg"), overwrite = TRUE)
  message("Supplementary A retained verbatim (technical QC): ", ORIG_SUPPA_SVG)
} else {
  warning("Original Supplementary A plot not found: ", ORIG_SUPPA_SVG)
}
if (file.exists(ORIG_SUPPA_DATA)) {
  file.copy(ORIG_SUPPA_DATA, file.path(STAGE, "SuppA_protein_peptide_counts_source_data.xlsx"), overwrite = TRUE)
  file.copy(ORIG_SUPPA_DATA, file.path(PANELS_SUP, "SuppA_protein_peptide_counts_source_data.xlsx"), overwrite = TRUE)
}

# ==============================================================================
# SUPPLEMENTARY B2 / B3 -- acquisition-level PCA, hemisphere and collection plate
# ==============================================================================
read_gct_matrix <- function(path) {
  hdr <- utils::read.delim(path, skip = 1, nrows = 1, header = FALSE)
  n_row <- as.integer(hdr[[1]]); n_col <- as.integer(hdr[[2]])
  n_rdesc <- as.integer(hdr[[3]]); n_cdesc <- as.integer(hdr[[4]])
  raw <- utils::read.delim(path, skip = 2, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  samples <- names(raw)[(2 + n_rdesc):ncol(raw)]
  body <- raw[(n_cdesc + 1):nrow(raw), , drop = FALSE]
  mat <- as.matrix(body[, samples, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- body[[1]]
  stopifnot(nrow(mat) == n_row, ncol(mat) == n_col)
  mat
}

mat <- read_gct_matrix(HEMI_GCT)
cat("hemisphere-level matrix:", nrow(mat), "proteins x", ncol(mat), "acquisitions\n")

meta <- utils::read.csv(LEGACY_META, stringsAsFactors = FALSE)
# join the long instrument column names to the legacy shortnames
idx <- vapply(meta$shortname, function(s) {
  hit <- grep(s, colnames(mat), fixed = TRUE)
  if (length(hit) == 1L) hit else NA_integer_
}, integer(1))
stopifnot(!anyNA(idx), !anyDuplicated(idx))
mat <- mat[, idx, drop = FALSE]
rownames(meta) <- colnames(mat)

# same preparation as the animal-level PCA: drop zero-variance rows, centre, scale
keep <- apply(mat, 1, function(x) { v <- stats::var(x, na.rm = TRUE); is.finite(v) && v > 0 })
complete <- stats::complete.cases(mat)
pmat <- mat[keep & complete, , drop = FALSE]
cat("proteins used:", nrow(pmat), "\n")
pca <- stats::prcomp(t(pmat), center = TRUE, scale. = TRUE)
varp <- (pca$sdev^2) / sum(pca$sdev^2) * 100
cat(sprintf("acquisition-level PC1 %.1f%%  PC2 %.1f%%\n", varp[1], varp[2]))

theme_pca_min <- function() {
  theme_minimal(base_size = 16, base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#ECECEC", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.title = element_text(color = "#444444", size = 18),
      axis.text  = element_text(color = "#555555", size = 18),
      axis.ticks = element_line(color = "#DDDDDD", linewidth = 0.3),
      axis.line = element_blank(), panel.border = element_blank(),
      legend.background = element_blank(), legend.key = element_blank(),
      legend.position = "right",
      legend.title = element_text(color = "#444444", size = 12),
      legend.text  = element_text(color = "#666666", size = 13),
      plot.title = element_text(color = "#222222", face = "bold", size = 16, margin = margin(b = 6)),
      plot.subtitle = element_text(color = "#666666", size = 11, margin = margin(b = 6)),
      plot.margin = margin(8, 8, 6, 8), complete = TRUE
    )
}

sc <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                 Hemisphere = factor(paste0("Hemisphere ", meta$ReplicateGroup)),
                 `Collection plate` = factor(paste0("Plate", meta$plate)),
                 sample_class = meta$celltype, AnimalID = meta$AnimalID,
                 condition_code = meta$ExpGroup,
                 check.names = FALSE, stringsAsFactors = FALSE)

qc_pca <- function(colvar, ttl, legend_title) {
  ggplot(sc, aes(PC1, PC2, color = .data[[colvar]])) +
    geom_point(size = 8, shape = 16, alpha = 0.8) +
    scale_color_manual(values = PCA_PALETTE[seq_len(nlevels(sc[[colvar]]))], name = legend_title) +
    theme_pca_min() +
    labs(x = sprintf("PC1 (%.1f%%)", varp[1]), y = sprintf("PC2 (%.1f%%)", varp[2]),
         title = ttl,
         subtitle = "Acquisition-level technical QC (96 acquisitions). Not a biological n.")
}

p_b2 <- qc_pca("Hemisphere", "Hemisphere", "Hemisphere")
save_panel(p_b2, STAGE, PANELS_SUP, "SuppB2_PCA_hemisphere_acquisition_level_QC.svg", width = 7.5, height = 6)

p_b3 <- qc_pca("Collection plate", "Collection plate", "Collection plate")
save_panel(p_b3, STAGE, PANELS_SUP, "SuppB3_PCA_collection_plate_acquisition_level_QC.svg", width = 7.5, height = 6)

src_b <- data.frame(acquisition = colnames(mat), shortname = meta$shortname, sc,
                    check.names = FALSE, stringsAsFactors = FALSE)
save_source_data(src_b, STAGE, PANELS_SUP, "SuppB2_B3_PCA_acquisition_level_QC_source_data.csv")

# how much PC1/PC2 each technical factor explains, so the caption can be quantitative
assoc <- do.call(rbind, lapply(c("Hemisphere", "Collection plate"), function(v) {
  do.call(rbind, lapply(c("PC1", "PC2"), function(pc) {
    fit <- stats::aov(sc[[pc]] ~ sc[[v]])
    s <- summary(fit)[[1]]
    data.frame(factor = v, PC = pc, F_value = s[["F value"]][1], p_value = s[["Pr(>F)"]][1],
               percent_variance_of_PC = 100 * s[["Sum Sq"]][1] / sum(s[["Sum Sq"]]),
               stringsAsFactors = FALSE)
  }))
}))
print(assoc)
utils::write.csv(assoc, file.path(STAGE, "SuppB2_B3_technical_factor_associations.csv"), row.names = FALSE)

audit_b <- data.frame(
  item = c("panel", "B1_status", "B2_status", "B3_status", "B1_representation", "B2_representation",
           "B3_representation", "hemisphere_gct", "hemisphere_gct_sha256", "n_acquisitions",
           "acquisition_PC1_percent", "acquisition_PC2_percent", "legacy_PC1_percent", "legacy_PC2_percent",
           "plate_naming", "plate_counts"),
  value = c("Supplementary B",
            "REGENERATED_ANIMAL_LEVEL (see full_regenerated/pca/SuppB1_PCA_experimental_group_animal_level.svg)",
            "VALID_AS_TECHNICAL_QC (acquisition level retained by necessity)",
            "VALID_AS_TECHNICAL_QC (acquisition level retained by necessity)",
            "animal-level PCA, 48 AnimalID x sample_class units, centred+scaled prcomp",
            "acquisition-level PCA, 96 acquisitions, centred+scaled prcomp",
            "acquisition-level PCA, 96 acquisitions, centred+scaled prcomp",
            HEMI_GCT, digest::digest(file = HEMI_GCT, algo = "sha256"),
            as.character(ncol(mat)), sprintf("%.1f", varp[1]), sprintf("%.1f", varp[2]),
            "25.4", "17.7",
            "'Plate' relabelled 'Collection plate'; NOT a demonstrated proteomics batch (single instrument, single acquisition date)",
            paste(paste0("Plate", names(table(meta$plate)), "=", as.integer(table(meta$plate))), collapse = "; ")),
  stringsAsFactors = FALSE
)
utils::write.csv(audit_b, file.path(STAGE, "SuppB_pca_audit.csv"), row.names = FALSE)

message("Supplementary A / B done.")
