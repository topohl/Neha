# ====================================================================
# PC contributions, pairwise PC matrix, variance barplot, distance heatmap, loading magnitude, trajectory, normalization comparison, metadata-PC correlation, 3D, densities
# Part of the PCA workflow split out of the former monolithic
# 03_qc_exploration/06_pcaPlot_animal_level.r (2026-08-26). Sourced in order by that
# script, which remains the entry point. Runs at top level and shares the
# globals created by 06a_pca_core.r (mat, meta, pca, output_dir, helpers).
# ====================================================================

if (!exists("mat") || !exists("meta") || !exists("pca") || !exists("output_dir")) {
  stop("PCA core state missing. Run 03_qc_exploration/06_pcaPlot_animal_level.r, or source pca/06a_pca_core.r first.", call. = FALSE)
}

# ================== Additional Supplementary Figure Plots ==================

# K) PC contribution heatmap: show which PCs each sample loads on
pc_contribution_heatmap <- function(n_pcs = 10) {
    pc_mat <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    pc_mat_scaled <- t(scale(t(pc_mat)))  # Scale per sample
    
    # Use sampleNumber as row names if available
    if ("sampleNumber" %in% names(meta)) {
        rownames(pc_mat_scaled) <- meta[rownames(pc_mat_scaled), "sampleNumber"]
    }
    
    # Annotate samples with metadata
    ann_cols <- intersect(c("sample_class", "condition_code", "ReplicateGroup", "plate"), names(meta))
    if (length(ann_cols) > 0) {
        annotation_row <- meta[rownames(pca$x), ann_cols, drop = FALSE]
        if ("sampleNumber" %in% names(meta)) {
            rownames(annotation_row) <- meta[rownames(pca$x), "sampleNumber"]
        }
        save_pheatmap(pc_mat_scaled, 
                     file.path(ensure_dir(subdir("plots/heatmaps")), "pc_contributions_per_sample.svg"),
                     width = 3, height = 13, scale = "none")
    }
}
pc_contribution_heatmap(10)

# L) Pairwise PC scatterplot matrix for PC1-PC5
pc_pairs_plot <- function(n_pcs = 5) {
    if (!requireNamespace("GGally", quietly = TRUE)) {
        message("GGally not installed; skipping pairs plot")
        return(NULL)
    }
    pc_subset <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    df <- data.frame(pc_subset, sample_class = meta[rownames(pc_subset), "sample_class"])
    pal <- make_modern_palette(nlevels(droplevels(factor(df$sample_class))))
    
    p <- GGally::ggpairs(df, columns = 1:ncol(pc_subset), 
                        aes(color = sample_class, alpha = 0.6),
                        upper = list(continuous = "cor"),
                        lower = list(continuous = "points")) +
        scale_color_manual(values = pal) +
        theme_pca_min()
    
    save_plot("plots/pairs", "pca_pairs_matrix.png", p, width = 12, height = 12)
}
pc_pairs_plot(5)

# M) Variance explained bar plot with cumulative line overlay
variance_barplot <- function() {
    var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
    cum_var <- cumsum(var_exp)
    n_show <- min(20, length(var_exp))
    
    df <- data.frame(
        PC = factor(paste0("PC", 1:n_show), levels = paste0("PC", 1:n_show)),
        Variance = var_exp[1:n_show],
        Cumulative = cum_var[1:n_show]
    )
    
    p <- ggplot(df, aes(x = PC)) +
        geom_col(aes(y = Variance), fill = "#6B5B95", alpha = 0.7) +
        geom_line(aes(y = Cumulative, group = 1), color = "#66C1A4", size = 1) +
        geom_point(aes(y = Cumulative), color = "#66C1A4", size = 2) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                          sec.axis = sec_axis(~., name = "Cumulative Variance")) +
        theme_pca_min() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "PCA Variance Explained", x = "Principal Component", y = "Individual Variance")
    
    save_plot("plots/variance", "pca_variance_barplot_with_cumulative.svg", p, width = 10, height = 6)
}
variance_barplot()

# N) Sample-to-sample distance heatmap (Euclidean distance in PC space)
sample_distance_heatmap <- function(n_pcs = 10) {
    pc_mat <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    dist_mat <- as.matrix(dist(pc_mat, method = "euclidean"))
    
    ann_cols <- intersect(c("sample_class", "condition_code", "plate"), names(meta))
    if (length(ann_cols) > 0) {
        annotation_col <- meta[colnames(dist_mat), ann_cols, drop = FALSE]
        annotation_row <- annotation_col
    }
    
    # Use sampleNumber as row names if available
    if ("sampleNumber" %in% names(meta)) {
        rownames(dist_mat) <- meta[colnames(dist_mat), "sampleNumber"]
    }
    
    save_pheatmap(dist_mat,
                 file.path(ensure_dir(subdir("plots/heatmaps")), "sample_distance_euclidean_pc_space.png"),
                 width = 10, height = 10, scale = "none")
}
sample_distance_heatmap(10)

# O) Loading magnitude plot: show absolute loadings across PCs
loading_magnitude_plot <- function(n_pcs = 5, top_n = 20) {
    rot <- as.data.frame(pca$rotation)
    pc_cols <- paste0("PC", 1:min(n_pcs, ncol(rot)))
    
    # Calculate magnitude across selected PCs
    rot$magnitude <- sqrt(rowSums(rot[, pc_cols]^2))
    rot$protein <- rownames(rot)
    
    top_proteins <- head(rot[order(-rot$magnitude), ], top_n)
    
    # Reshape for plotting
    top_long <- reshape2::melt(top_proteins[, c("protein", pc_cols)], 
                              id.vars = "protein", 
                              variable.name = "PC", 
                              value.name = "loading")
    
    p <- ggplot(top_long, aes(x = reorder(protein, loading), y = loading, fill = PC)) +
        geom_col(position = "dodge") +
        coord_flip() +
        scale_fill_manual(values = make_modern_palette(n_pcs)) +
        theme_pca_min() +
        labs(title = paste("Top", top_n, "Proteins by Loading Magnitude"),
             x = "Protein", y = "Loading Value")
    
    save_plot("plots/loadings", "loading_magnitude_top_proteins.svg", p, width = 8, height = 10)
}
loading_magnitude_plot(5, 20)

# P) PC trajectory plot: connect samples by experimental order if available
pc_trajectory_plot <- function() {
    if (!"sampleNumber" %in% names(meta)) {
        message("sampleNumber not in metadata; skipping trajectory plot")
        return(NULL)
    }
    
    df <- data.frame(pca$x[, 1:2], 
                     sampleNumber = meta[rownames(pca$x), "sampleNumber"],
                     sample_class = meta[rownames(pca$x), "sample_class"])
    df <- df[order(df$sampleNumber), ]
    
    pal <- make_modern_palette(nlevels(droplevels(factor(df$sample_class))))
    
    p <- ggplot(df, aes(PC1, PC2, color = sample_class)) +
        geom_path(color = "gray70", alpha = 0.5, size = 0.5) +
        geom_point(size = 4, alpha = 0.8) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "PCA Trajectory by Sample Order", 
             subtitle = "Lines connect consecutive samples")
    
    save_plot("plots/trajectory", "pca_sample_trajectory.svg", p)
}
pc_trajectory_plot()

# Q) Comparison of different normalization methods
compare_normalizations <- function() {
    if (!exists("mat") || is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
        message("Matrix 'mat' is not defined or is empty; skipping normalization comparison.")
        return(NULL)
    }
    
    # Check if matrix has valid data
    if (!is.matrix(mat) || !is.numeric(mat)) {
        message("Matrix 'mat' is not a numeric matrix; skipping normalization comparison.")
        return(NULL)
    }
    
    # Remove rows/columns with all NAs or infinite values
    mat_clean <- mat[rowSums(is.finite(mat)) > 0, , drop = FALSE]
    mat_clean <- mat_clean[, colSums(is.finite(mat_clean)) > 0, drop = FALSE]
    
    if (nrow(mat_clean) < 2 || ncol(mat_clean) < 2) {
        message("Matrix too small after cleaning; skipping normalization comparison.")
        return(NULL)
    }
    
    # Current data (already processed) - use as-is
    pca_current <- prcomp(t(mat_clean), center = TRUE, scale. = TRUE)
    
    # Z-score normalization (per protein)
    mat_zscore <- t(scale(t(mat_clean), center = TRUE, scale = TRUE))
    pca_zscore <- prcomp(t(mat_zscore), center = TRUE, scale. = TRUE)
    
    # Median-centered normalization
    mat_median <- sweep(mat_clean, 1, apply(mat_clean, 1, median, na.rm = TRUE), "-")
    pca_median <- prcomp(t(mat_median), center = TRUE, scale. = TRUE)
    
    # Plot comparison
    plot_list <- list(
        list(pca = pca_current, title = "Current (PCA-adjusted)", file = "norm_current.svg"),
        list(pca = pca_zscore, title = "Z-score normalized", file = "norm_zscore.svg"),
        list(pca = pca_median, title = "Median-centered", file = "norm_median.svg")
    )
    
    for (item in plot_list) {
        # Get common samples between PCA and metadata
        common_samples <- intersect(rownames(item$pca$x), rownames(meta))
        if (length(common_samples) < 2) {
            message(paste("Skipping", item$title, "- no common samples with metadata"))
            next
        }
        
        df <- data.frame(
            item$pca$x[common_samples, 1:2, drop = FALSE],
            sample_class = meta[common_samples, "sample_class"]
        )
        
        pal <- make_modern_palette(nlevels(droplevels(factor(df$sample_class))))
        
        p <- ggplot(df, aes(PC1, PC2, color = sample_class)) +
            geom_point(size = 4, alpha = 0.8, stroke = 0) +
            scale_color_manual(values = pal) +
            theme_pca_min() +
            labs(title = paste("PCA:", item$title), color = "sample_class")
        
        save_plot("plots/normalization", item$file, p)
    }
    
    message("Normalization comparison completed successfully!")
}
compare_normalizations()

# R) Metadata correlation heatmap with PCs
metadata_pc_correlation <- function(n_pcs = 10) {
    pc_mat <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    
    # Select numeric metadata columns
    numeric_meta <- meta[, sapply(meta, is.numeric), drop = FALSE]
    
    if (ncol(numeric_meta) == 0) {
        message("No numeric metadata found for correlation")
        return(NULL)
    }
    
    # Align samples
    common_samples <- intersect(rownames(pc_mat), rownames(numeric_meta))
    pc_mat <- pc_mat[common_samples, ]
    numeric_meta <- numeric_meta[common_samples, ]
    
    # Calculate correlations
    cor_mat <- cor(pc_mat, numeric_meta, use = "pairwise.complete.obs")
    
    save_pheatmap(t(cor_mat),
                 file.path(ensure_dir(subdir("plots/heatmaps")), "metadata_pc_correlations.png"),
                 width = 8, height = 6, scale = "none")
    
    # Also save as table
    cor_df <- as.data.frame(as.table(cor_mat))
    colnames(cor_df) <- c("PC", "Metadata", "Correlation")
    save_table("tables/correlations", "metadata_pc_correlations.csv", cor_df)
}
metadata_pc_correlation(10)

# S) 3D PCA plot (PC1, PC2, PC3) if plotly available
pca_3d_plot <- function() {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        message("plotly not installed; skipping 3D PCA plot")
        return(NULL)
    }
    
    df <- data.frame(pca$x[, 1:3],
                     sample_class = meta[rownames(pca$x), "sample_class"],
                     sample = rownames(pca$x))
    
    pal <- make_modern_palette(nlevels(droplevels(factor(df$sample_class))))
    
    p <- plotly::plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, 
                        color = ~sample_class, colors = pal,
                        text = ~sample, type = "scatter3d", mode = "markers",
                        marker = list(size = 5, opacity = 0.8)) %>%
        plotly::layout(title = "3D PCA Plot (PC1-PC3)",
                      scene = list(xaxis = list(title = "PC1"),
                                  yaxis = list(title = "PC2"),
                                  zaxis = list(title = "PC3")))
    
    htmlwidgets::saveWidget(p, 
                           file.path(ensure_dir(subdir("plots/3d")), "pca_3d_interactive.html"),
                           selfcontained = FALSE)
}
pca_3d_plot()

# T) Density plots of PC scores by group
pc_density_plots <- function(n_pcs = 3) {
    for (i in 1:min(n_pcs, ncol(pca$x))) {
        pc_name <- paste0("PC", i)
        df <- data.frame(
            score = pca$x[, i],
            sample_class = meta[rownames(pca$x), "sample_class"]
        )
        
        pal <- make_modern_palette(nlevels(droplevels(factor(df$sample_class))))
        
        p <- ggplot(df, aes(x = score, fill = sample_class)) +
            geom_density(alpha = 0.6) +
            scale_fill_manual(values = pal) +
            theme_pca_min() +
            labs(title = paste(pc_name, "Distribution by sample_class"),
                 x = paste(pc_name, "Score"), y = "Density")
        
        save_plot("plots/density", paste0(pc_name, "_density_by_sample_class.svg"), p)
    }
}
pc_density_plots(3)

message("All supplementary plots completed!")

