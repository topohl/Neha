# ====================================================================
# PCA core: config, helpers, validated input, base plots, audit
# Part of the PCA workflow split out of the former monolithic
# 03_qc_exploration/06_pcaPlot_Neha.r (2026-08-26). Sourced in order by that
# script, which remains the entry point. Runs at top level and shares the
# globals created by 06a_pca_core.r (mat, meta, pca, output_dir, helpers).
# Defines mat/meta/pca and the shared helpers every other part uses.
# ====================================================================

# ================== PCA analysis and plotting (extended, fixed, organized) ==================
# Set a reliable CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
}

# Load core packages
pacman::p_load(
    data.table, ggplot2, factoextra, reshape2, stats, ggrepel, tools,
    grid, uwot, RColorBrewer, pheatmap, Rtsne, aricode, rospca, irlba,
    pandoc, treemapify, dplyr, digest
)

# Try to load aricode separately (it's only used for clustering metrics)
if (!requireNamespace("aricode", quietly = TRUE)) {
    message("aricode package not available - clustering ARI/NMI metrics will be skipped")
} else {
    library(aricode)
}

set.seed(42)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path("03_qc_exploration", "06_pcaPlot_Neha.r")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "protigy_input_utils.R"))) {
    repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "protigy_input_utils.R"))
source(file.path(repo_root, "R", "neha_path_utils.R"))
source(file.path(repo_root, "R", "pca_animal_level_utils.R"))

# =============== Config =================
option_or_env <- function(option_name, env_name, default) {
    option_value <- getOption(option_name)
    if (!is.null(option_value) && nzchar(trimws(as.character(option_value)))) return(as.character(option_value))
    env_value <- Sys.getenv(env_name, unset = "")
    if (nzchar(trimws(env_value))) return(env_value)
    default
}

gct_file <- option_or_env(
    "neha.pca_animal_level_input",
    "NEHA_PCA_ANIMAL_LEVEL_INPUT",
    "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/02_data/animal_level/input_gct/neha_protigy_input_animal_level_primary.gct"
)
# Legacy hemisphere-level PCA outputs now live under 99_historical/pca_plots_legacy.
historical_output_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/99_historical/pca_plots_legacy"
output_dir <- validate_neha_pca_output_root(
    option_or_env(
        "neha.pca_output_root",
        "NEHA_PCA_OUTPUT_ROOT",
        "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/pca"
    ),
    historical_output_dir
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
pca_audit_path <- file.path(output_dir, "tables", "meta", "animal_level_pca_audit.csv")
optional_visualization_audit_path <- file.path(output_dir, "tables", "meta", "optional_visualization_audit.csv")

# =============== Helpers =================
trim_ws <- function(x){
    if (is.null(x)) return(x)
    x <- as.character(x)
    x <- gsub("[\u00A0\u2007\u202F]", " ", x, perl = TRUE)
    x <- gsub("^\\s+|\\s+$", "", x, perl = TRUE)
    x
}

# robust writer: prefer fwrite, else write.csv
write_dt <- function(df, path){
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    if (requireNamespace("data.table", quietly = TRUE)) {
        data.table::fwrite(df, path)
    } else {
        utils::write.csv(df, path, row.names = FALSE)
    }
    invisible(path)
}

# pheatmap saver that forces supported extensions
save_pheatmap <- function(mat, file, width=8, height=10, scale="row"){
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        message("pheatmap not installed, skipping heatmap: ", file)
        return(invisible(NULL))
    }
    ext <- tools::file_ext(file)
    if (!nzchar(ext) || !(tolower(ext) %in% c("png","pdf","tiff","bmp","jpeg","jpg"))) {
        file <- sub("\\.[A-Za-z0-9]+$", "", file)
        file <- paste0(file, ".png")
    }
    pheatmap::pheatmap(
        mat, show_colnames = FALSE, scale = scale, clustering_method = "complete",
        color = colorRampPalette(c("#2c7fb8","#f7f7f7","#d95f0e"))(101),
        filename = file, width = width, height = height
    )
    invisible(file)
}

# -------- Subfolder helpers ----------
subdir <- function(...){ file.path(output_dir, ...) }
ensure_dir <- function(path){ dir.create(path, showWarnings = FALSE, recursive = TRUE); path }

save_plot <- function(subfolder, filename, plot, width=7.5, height=6.2, dpi=150){
    d <- ensure_dir(subdir(subfolder))
    ggsave(filename = filename, plot = plot, path = d, width = width, height = height, dpi = dpi)
    invisible(file.path(d, filename))
}

save_table <- function(subfolder, filename, df, row.names=FALSE){
    d <- ensure_dir(subdir(subfolder))
    f <- file.path(d, filename)
    if (exists("write_dt")) write_dt(df, f) else utils::write.csv(df, f, row.names = row.names)
    invisible(f)
}

optional_visualization_audit <- data.frame(
    visualization = character(), execution_status = character(), reason = character(),
    n_input_rows = integer(), n_usable_rows = integer(), output_path = character(),
    stringsAsFactors = FALSE
)
record_optional_visualization <- function(result) {
    optional_visualization_audit <<- rbind(optional_visualization_audit, result)
    write_dt(optional_visualization_audit, optional_visualization_audit_path)
    invisible(result)
}

# ================== Build mat/meta with checks ==================
cat("Reading validated animal-level GCT:", gct_file, "\n")
source_sha256 <- if (file.exists(gct_file)) digest::digest(file = gct_file, algo = "sha256") else NA_character_
core_input <- tryCatch({
    parsed_gct <- validate_protigy_gct_v13(gct_file)
    validated_input <- validate_neha_pca_animal_input(parsed_gct, expected_n = 3L)
    prepared_pca <- prepare_neha_animal_pca(
        validated_input$expression_matrix,
        center = TRUE,
        scale. = TRUE
    )
    list(validated = validated_input, prepared = prepared_pca)
}, error = function(e) {
    failure_audit <- make_neha_pca_audit(
        source_path = gct_file,
        source_sha256 = source_sha256,
        output_paths = c(pca_audit = pca_audit_path),
        execution_status = "failed",
        error_message = conditionMessage(e)
    )
    write_dt(failure_audit, pca_audit_path)
    stop(e)
})
mat <- core_input$prepared$matrix
meta <- core_input$validated$sample_metadata
pca <- core_input$prepared$pca

# ================== PCA pipeline ==================
stopifnot(identical(rownames(meta), rownames(pca$x)))

# Derived labels based on actual metadata fields
# AnimalID, condition_code, condition, sample_class, and phenotypeWithinUnit are validated GCT metadata.
meta$sample_class_condition <- paste(meta$sample_class, meta$condition, sep = "_")

# ================== Minimal modern styling ==================
theme_pca_min <- function() {
  theme_minimal(base_size = 16, base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#ECECEC", size = 0.4),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.title = element_text(color = "#444444", size = 18, face = "plain"),
      axis.text  = element_text(color = "#555555", size = 18),
      axis.ticks = element_line(color = "#DDDDDD", size = 0.3),
      axis.line  = element_blank(),
      panel.border = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = "right",
      legend.title = element_text(color = "#444444", size = 12),
      legend.text  = element_text(color = "#666666", size = 13),
      plot.title   = element_text(color = "#222222", face = "bold", size = 16, margin = margin(b = 6)),
      plot.subtitle= element_text(color = "#666666", size = 12, margin = margin(b = 6)),
      plot.caption = element_text(color = "#999999", size = 10),
      strip.text = element_text(color = "#333333", size = 10, face = "bold"),
      panel.spacing = unit(0.6, "lines"),
      plot.margin = margin(8, 8, 6, 8),
      complete = TRUE
    )
}

#base_hex <- c("#F8D247", "#E89369", "#B2B2B2", "#C7D745", "#C6B18B",
#              "#C592C5", "#A89BB0", "#DEC196", "#66C1A4", "#FF6F61",
#              "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251")

#base_hex <- c("#FFB6C1", "#EE82EE", "#FF69B4", "#FF1493", "#C71585",
#              "#E6B0E8", "#DDA0DD", "#ff68c0ff", "#DA70D6", "#BA55D3",
#              "#F4C2C2", "#FFD1DC", "#FFDAB9", "#FFE4E1", "#FFF0F5")

base_hex <- c("#1E90FF", "#61d7ffff", "#FFA500", "#FFD700", "#FF4500",
              "#FF6347", "#FFDAB9", "#FFE4B5", "#F0E68C", "#ADFF2F",
              "#98FB98", "#FF69B4", "#FFB6C1", "#FF1493", "#FF4500")

make_modern_palette <- function(n){
  if (n <= length(base_hex)) return(base_hex[seq_len(n)])
  grDevices::colorRampPalette(base_hex, space = "Lab")(n)
}

build_group <- function(key){
  if (!key %in% names(meta)) stop(sprintf("Key '%s' not in meta.", key))
  v <- trim_ws(as.character(meta[[key]]))
  v[v == ""] <- NA
  factor(v)
}

plot_and_save_group <- function(key, title_prefix, out_file, point_size = 8) {
  grp <- build_group(key)
  keep <- !is.na(grp)
  if (sum(keep) < 2) { message(sprintf("Skipping '%s': <2 non-NA samples.", key)); return(invisible(NULL)) }
  grp2 <- droplevels(grp[keep]); nlev <- nlevels(grp2); pal <- make_modern_palette(nlev)

  X <- pca; X$x <- X$x[keep, , drop = FALSE]

  if (!is.null(pca$sdev) && length(pca$sdev) >= 2) {
    varp <- (pca$sdev^2) / sum(pca$sdev^2)
    lab_x <- sprintf("PC1 (%.1f%%)", varp[1] * 100)
    lab_y <- sprintf("PC2 (%.1f%%)", varp[2] * 100)
  } else {
    lab_x <- "PC1"; lab_y <- "PC2"
  }

  p <- fviz_pca_ind(
    X,
    geom = "point",
    habillage = grp2,
    addEllipses = FALSE,
    ellipse.type = "t",
    ellipse.level = 0.5,
    palette = pal,
    repel = TRUE,
    mean.point = FALSE,
    title = sprintf("%s by %s", title_prefix, key)
  )

  if (length(p$layers) && inherits(p$layers[[1]]$geom, "GeomPoint")) {
    p$layers[[1]] <- ggplot2::geom_point(mapping = p$layers[[1]]$mapping, inherit.aes = TRUE, shape = 16, size = point_size, alpha = 0.8)
  } else {
    p <- p + ggplot2::geom_point(size = point_size, shape = 16, alpha = 0.8)
  }

  p <- p +
    theme_pca_min() +
    theme(
      axis.line.x = element_line(color = "#E0E0E0"),
      axis.line.y = element_line(color = "#E0E0E0")
    ) +
    labs(x = lab_x, y = lab_y, subtitle = NULL, caption = NULL)

  # Organized save
  save_plot("plots/base", out_file, p)
  p
}

# ================== Generate and save base plots (SVG) ==================
primary_output_paths <- c(
  pca_by_sample_class = file.path(output_dir, "plots", "base", "pca_by_sample_class.svg"),
  pca_by_condition = file.path(output_dir, "plots", "base", "pca_by_condition.svg"),
  pca_by_animal_id = file.path(output_dir, "plots", "base", "pca_by_animal_id.svg"),
  pca_by_phenotype = file.path(output_dir, "plots", "base", "pca_by_phenotype_within_unit.svg"),
  pca_scree = file.path(output_dir, "plots", "variance", "pca_scree.svg"),
  pca_cumulative_variance = file.path(output_dir, "plots", "variance", "pca_cumulative_variance.svg"),
  parsed_metadata = file.path(output_dir, "tables", "meta", "sample_metadata_parsed.csv"),
  variance_explained = file.path(output_dir, "tables", "variance", "pca_variance_explained.csv"),
  pca_audit = pca_audit_path
)
tryCatch({
  plot_and_save_group("sample_class", "PCA", "pca_by_sample_class.svg")
  plot_and_save_group("condition", "PCA", "pca_by_condition.svg")
  plot_and_save_group("AnimalID", "PCA", "pca_by_animal_id.svg")
  plot_and_save_group("phenotypeWithinUnit", "PCA", "pca_by_phenotype_within_unit.svg")

  # Export parsed metadata
  save_table("tables/meta", "sample_metadata_parsed.csv", meta, row.names = TRUE)

  # 1) Scree and cumulative variance
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
  df_scree <- data.frame(PC = seq_along(var_explained),
                         Variance = var_explained,
                         Cumulative = cumsum(var_explained))
  save_table("tables/variance", "pca_variance_explained.csv", df_scree)

  p_scree <- ggplot(df_scree, aes(PC, Variance)) +
    geom_col(fill="#6B5B95") +
    geom_point() + geom_line(group=1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_pca_min() + labs(title="PCA Scree", x="Principal Component", y="Variance Explained")
  save_plot("plots/variance", "pca_scree.svg", p_scree)

  p_cum <- ggplot(df_scree, aes(PC, Cumulative)) +
    geom_point(color="#66C1A4") + geom_line(color="#66C1A4") +
    geom_hline(yintercept = 0.8, linetype="dashed", color="#B2B2B2") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,1)) +
    theme_pca_min() + labs(title="Cumulative Variance", x="Principal Component", y="Cumulative Fraction")
  save_plot("plots/variance", "pca_cumulative_variance.svg", p_cum)
}, error = function(e) {
  failure_audit <- make_neha_pca_audit(
    validated_input = core_input$validated,
    prepared_pca = core_input$prepared,
    source_path = gct_file,
    source_sha256 = source_sha256,
    output_paths = primary_output_paths,
    execution_status = "failed",
    error_message = conditionMessage(e)
  )
  write_dt(failure_audit, pca_audit_path)
  stop(e)
})
pca_audit <- make_neha_pca_audit(
  validated_input = core_input$validated,
  prepared_pca = core_input$prepared,
  source_path = gct_file,
  source_sha256 = source_sha256,
  output_paths = primary_output_paths,
  execution_status = "success"
)
write_dt(pca_audit, pca_audit_path)

