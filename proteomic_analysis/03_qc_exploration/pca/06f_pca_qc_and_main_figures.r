# ====================================================================
# Sample QC vs PCs, dendrogram, Kaiser scree, PC time series, batch visualization, report-quality panels, main-figure ellipses and loading plots
# Part of the PCA workflow split out of the former monolithic
# 03_qc_exploration/06_pcaPlot_Neha.r (2026-08-26). Sourced in order by that
# script, which remains the entry point. Runs at top level and shares the
# globals created by 06a_pca_core.r (mat, meta, pca, output_dir, helpers).
# ====================================================================

if (!exists("mat") || !exists("meta") || !exists("pca") || !exists("output_dir")) {
  stop("PCA core state missing. Run 03_qc_exploration/06_pcaPlot_Neha.r, or source pca/06a_pca_core.r first.", call. = FALSE)
}

# Y) Sample QC metrics vs PC scores
sample_qc_vs_pcs <- function(n_pcs = 2) {
    # Calculate QC metrics
    qc_df <- data.frame(
        sample = colnames(mat),
        n_detected = colSums(!is.na(mat) & mat > 0),
        mean_intensity = colMeans(mat, na.rm = TRUE),
        median_intensity = apply(mat, 2, median, na.rm = TRUE),
        missing_pct = colMeans(is.na(mat)) * 100
    )
    
    # Merge with PC scores
    qc_df <- cbind(qc_df, pca$x[qc_df$sample, 1:min(n_pcs, ncol(pca$x)), drop = FALSE])
    
    # Plot each QC metric vs PCs
    for (qc_col in c("n_detected", "mean_intensity", "missing_pct")) {
        for (i in 1:min(n_pcs, ncol(pca$x))) {
            pc_name <- paste0("PC", i)
            
            p <- ggplot(qc_df, aes_string(x = pc_name, y = qc_col)) +
                geom_point(alpha = 0.7, size = 3, color = "#6B5B95") +
                geom_smooth(method = "lm", se = TRUE, color = "#66C1A4") +
                theme_pca_min() +
                labs(title = paste(qc_col, "vs", pc_name),
                     x = paste(pc_name, "Score"),
                     y = qc_col)
            
            save_plot("plots/qc", paste0(qc_col, "_vs_", pc_name, ".svg"), p)
        }
    }
    
    # Save QC table
    save_table("tables/qc", "sample_qc_metrics.csv", qc_df)
}
sample_qc_vs_pcs(3)

# Z) Hierarchical clustering dendrogram colored by metadata
sample_dendrogram <- function(n_pcs = 10, group_key = "sample_class") {
    if (!group_key %in% names(meta)) return(NULL)
    
    pc_mat <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    hc <- hclust(dist(pc_mat), method = "complete")
    
    # Color by group
    groups <- factor(meta[rownames(pc_mat), group_key])
    pal <- make_modern_palette(nlevels(droplevels(groups)))
    colors <- pal[as.numeric(droplevels(groups))]
    
    png(file.path(ensure_dir(subdir("plots/dendrogram")), 
                  paste0("dendrogram_by_", group_key, ".png")),
        width = 12, height = 8, units = "in", res = 150)
    plot(hc, hang = -1, cex = 0.6, main = paste("Sample Dendrogram colored by", group_key))
    rect.hclust(hc, k = nlevels(droplevels(groups)), border = pal)
    dev.off()
}
sample_dendrogram(10, "sample_class")

# AA) Scree plot with Kaiser criterion line (eigenvalue = 1)
scree_with_kaiser <- function() {
    eigenvalues <- pca$sdev^2
    df <- data.frame(
        PC = factor(paste0("PC", 1:length(eigenvalues)), 
                   levels = paste0("PC", 1:length(eigenvalues))),
        Eigenvalue = eigenvalues
    )
    df <- head(df, 20)  # Show first 20
    
    p <- ggplot(df, aes(x = PC, y = Eigenvalue)) +
        geom_col(fill = "#6B5B95", alpha = 0.7) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
        annotate("text", x = 2, y = 1.2, label = "Kaiser criterion (λ=1)", 
                color = "red", hjust = 0) +
        theme_pca_min() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Scree Plot with Kaiser Criterion",
             x = "Principal Component", y = "Eigenvalue")
    
    save_plot("plots/variance", "scree_with_kaiser_criterion.svg", p)
}
scree_with_kaiser()

# AB) PC time series (if temporal/sequential metadata available)
pc_timeseries <- function(time_key = "sampleNumber", group_key = "sample_class") {
    if (!time_key %in% names(meta) || !group_key %in% names(meta)) {
        message("Required metadata not available for time series")
        return(NULL)
    }
    
    df <- data.frame(
        pca$x[, 1:2],
        time = meta[rownames(pca$x), time_key],
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$time) & !is.na(df$group), ]
    df <- df[order(df$time), ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(df, aes(x = time, y = PC1, color = group)) +
        geom_line(alpha = 0.6) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "PC1 Over Sample Order",
             x = "Sample Number", y = "PC1 Score", color = group_key)
    
    save_plot("plots/timeseries", "pc1_timeseries.svg", p)
}
pc_timeseries("sampleNumber", "sample_class")

# AC) Batch effect visualization before/after correction
batch_effect_comparison <- function() {
    if (!"plate" %in% names(meta)) return(NULL)
    
    # Before correction
    df_pre <- data.frame(pca$x[, 1:2], plate = meta[rownames(pca$x), "plate"])
    pal <- make_modern_palette(nlevels(droplevels(factor(df_pre$plate))))
    
    p1 <- ggplot(df_pre, aes(PC1, PC2, color = factor(plate))) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "Before Batch Correction", color = "Plate")
    
    save_plot("plots/batch", "batch_effect_before.svg", p1)
}
batch_effect_comparison()

message("All additional supplementary and main figure plots completed!")

# ================== Additional Main Figure Plots ==================

# AD) Side-by-side PCA plots for different groupings (report quality)
combined_pca_panels <- function() {
    # Create 2x2 panel of key groupings
    groupings <- c("sample_class", "condition", "AnimalID", "phenotypeWithinUnit")
    groupings <- intersect(groupings, names(meta))
    
    if (length(groupings) < 2) {
        message("Not enough grouping variables for panel plot")
        return(NULL)
    }
    
    plot_list <- list()
    for (key in groupings) {
        grp <- build_group(key)
        keep <- !is.na(grp)
        if (sum(keep) < 2) next
        
        grp2 <- droplevels(grp[keep])
        df <- data.frame(pca$x[keep, 1:2, drop=FALSE], group = grp2)
        pal <- make_modern_palette(nlevels(grp2))
        
        varp <- (pca$sdev^2)/sum(pca$sdev^2)
        
        p <- ggplot(df, aes(PC1, PC2, color = group)) +
            geom_point(size = 4, alpha = 0.8, stroke = 0) +
            scale_color_manual(values = pal) +
            theme_pca_min() +
            labs(x = sprintf("PC1 (%.1f%%)", 100*varp[1]),
                 y = sprintf("PC2 (%.1f%%)", 100*varp[2]),
                 color = key,
                 title = key)
        
        plot_list[[key]] <- p
    }
    
    if (length(plot_list) > 0 && requireNamespace("gridExtra", quietly = TRUE)) {
        combined <- gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
        ggsave(file.path(ensure_dir(subdir("plots/main_figure")), "pca_combined_panels.svg"),
               combined, width = 14, height = 12, dpi = 300)
    }
}
combined_pca_panels()

# AE) Report-ready variance explained with cumulative overlay
variance_main_figure <- function() {
    var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
    cum_var <- cumsum(var_exp)
    n_show <- min(15, length(var_exp))
    
    df <- data.frame(
        PC = 1:n_show,
        PC_label = factor(paste0("PC", 1:n_show), levels = paste0("PC", 1:n_show)),
        Individual = var_exp[1:n_show] * 100,
        Cumulative = cum_var[1:n_show] * 100
    )
    
    # Dual-axis plot
    p <- ggplot(df, aes(x = PC)) +
        geom_col(aes(y = Individual), fill = "#6B5B95", alpha = 0.8, width = 0.7) +
        geom_line(aes(y = Cumulative), color = "#66C1A4", size = 1.5, group = 1) +
        geom_point(aes(y = Cumulative), color = "#66C1A4", size = 3) +
        geom_hline(yintercept = 80, linetype = "dashed", color = "#E89369", size = 0.8) +
        annotate("text", x = n_show * 0.8, y = 82, label = "80% threshold", 
                color = "#E89369", size = 4) +
        scale_x_continuous(breaks = 1:n_show, labels = paste0("PC", 1:n_show)) +
        scale_y_continuous(limits = c(0, 100)) +
        theme_pca_min() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            panel.grid.major.x = element_blank()
        ) +
        labs(title = "Variance Explained by Principal Components",
             x = "Principal Component",
             y = "Variance Explained (%)")
    
    save_plot("plots/main_figure", "variance_explained_main.svg", p, width = 10, height = 6)
}
variance_main_figure()

# AF) Top contributing proteins heatmap (main figure version)
top_proteins_heatmap_main <- function(n_pcs = 3, top_per_pc = 15) {
    if (is.null(pca) || is.null(pca$rotation) || ncol(pca$rotation) == 0 || 
        is.null(mat) || ncol(mat) == 0 || nrow(mat) == 0) {
        message("PCA object or matrix is not properly defined")
        return(NULL)
    }
    
    rot <- as.data.frame(pca$rotation)
    
    # Get top contributors for each PC
    top_proteins <- c()
    for (i in 1:min(n_pcs, ncol(rot))) {
        pc_col <- paste0("PC", i)
        if (!pc_col %in% colnames(rot)) next
        
        # Top positive
        top_pos <- head(names(sort(rot[[pc_col]], decreasing = TRUE)), top_per_pc)
        # Top negative
        top_neg <- head(names(sort(rot[[pc_col]], decreasing = FALSE)), top_per_pc)
        top_proteins <- c(top_proteins, top_pos, top_neg)
    }
    top_proteins <- unique(top_proteins)
    
    # Filter to proteins that exist in mat
    top_proteins <- intersect(top_proteins, rownames(mat))
    
    if (length(top_proteins) < 4) {
        message("Not enough valid proteins found for heatmap")
        return(NULL)
    }
    
    # Create expression matrix for these proteins
    expr_subset <- mat[top_proteins, , drop = FALSE]
    
    # Remove proteins with all NA or no variation
    valid_rows <- rowSums(!is.na(expr_subset)) > 0 & 
                  apply(expr_subset, 1, function(x) var(x, na.rm = TRUE) > 0)
    expr_subset <- expr_subset[valid_rows, , drop = FALSE]
    
    if (nrow(expr_subset) < 4) {
        message("Not enough valid proteins after filtering")
        return(NULL)
    }
    
    # Order samples by PC1
    sample_order <- order(pca$x[, 1])
    expr_subset <- expr_subset[, sample_order, drop = FALSE]
    
    # Annotation
    ann_col <- data.frame(
        sample_class = meta[colnames(expr_subset), "sample_class"],
        row.names = colnames(expr_subset)
    )
    
    ann_colors <- list(
        sample_class = setNames(make_modern_palette(nlevels(factor(ann_col$sample_class))),
                           levels(factor(ann_col$sample_class)))
    )
    
    pheatmap::pheatmap(
        expr_subset,
        annotation_col = ann_col,
        annotation_colors = ann_colors,
        show_colnames = FALSE,
        scale = "row",
        clustering_method = "complete",
        color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
        filename = file.path(ensure_dir(subdir("plots/main_figure")), 
                           "top_proteins_heatmap_main.png"),
        width = 10,
        height = 12
    )
}
top_proteins_heatmap_main(3, 15)

# AG) PCA with confidence ellipses (main figure version)
pca_with_ellipses_main <- function(group_key = "sample_class", conf_level = 0.95) {
    if (!group_key %in% names(meta)) return(NULL)
    
    grp <- build_group(group_key)
    keep <- !is.na(grp)
    grp2 <- droplevels(grp[keep])
    
    if (nlevels(grp2) < 2) {
        message(paste("Not enough groups for", group_key))
        return(NULL)
    }
    
    df <- data.frame(pca$x[keep, 1:2, drop=FALSE], group = grp2)
    pal <- make_modern_palette(nlevels(grp2))
    
    varp <- (pca$sdev^2)/sum(pca$sdev^2)
    
    p <- ggplot(df, aes(PC1, PC2, color = group, fill = group)) +
        stat_ellipse(type = "t", level = conf_level, geom = "polygon", 
                    alpha = 0.15, show.legend = FALSE) +
        stat_ellipse(type = "t", level = conf_level, size = 1, show.legend = FALSE) +
        geom_point(size = 5, alpha = 0.8, stroke = 0, shape = 16) +
        scale_color_manual(values = pal) +
        scale_fill_manual(values = pal) +
        theme_pca_min() +
        theme(
            legend.position = "right",
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 16, face = "bold")
        ) +
        labs(x = sprintf("PC1 (%.1f%%)", 100*varp[1]),
             y = sprintf("PC2 (%.1f%%)", 100*varp[2]),
             color = group_key,
             fill = group_key,
             title = paste("PCA with", conf_level*100, "% Confidence Ellipses"))
    
    save_plot("plots/main_figure", paste0("pca_ellipses_", group_key, "_main.svg"), p, width = 10, height = 8)
}
pca_with_ellipses_main("sample_class", 0.95)
pca_with_ellipses_main("condition_code", 0.95)

# AH) Loading plot for top contributing proteins (main figure)
loadings_main_figure <- function(pc_x = 1, pc_y = 2, top_n = 25) {
    rot <- as.data.frame(pca$rotation)
    
    # Calculate contribution magnitude
    mag <- sqrt(rot[, pc_x]^2 + rot[, pc_y]^2)
    top_idx <- order(mag, decreasing = TRUE)[1:min(top_n, length(mag))]
    
    df <- data.frame(
        protein = rownames(rot)[top_idx],
        PC1_load = rot[top_idx, pc_x],
        PC2_load = rot[top_idx, pc_y],
        magnitude = mag[top_idx]
    )
    
    # Create arrow plot
    p <- ggplot(df, aes(x = PC1_load, y = PC2_load)) +
        geom_segment(aes(x = 0, y = 0, xend = PC1_load, yend = PC2_load),
                    arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                    color = "#6B5B95", alpha = 0.7, size = 0.8) +
        geom_text_repel(aes(label = protein), size = 3.5, 
                       max.overlaps = 30, box.padding = 0.5) +
        geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
        coord_fixed() +
        theme_pca_min() +
        labs(title = paste("Top", top_n, "Contributing Proteins"),
             x = paste0("PC", pc_x, " Loading"),
             y = paste0("PC", pc_y, " Loading"))
    
    save_plot("plots/main_figure", paste0("loadings_arrows_PC", pc_x, "_PC", pc_y, "_main.svg"), 
             p, width = 10, height = 10)
}
loadings_main_figure(1, 2, 25)

