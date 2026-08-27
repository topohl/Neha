# ====================================================================
# Arc, radial, alluvial, ridgeline, hexbin, loading network, parallel coordinates, variance treemap, contour, chord, small multiples, UMAP hexbin marginals
# Part of the PCA workflow split out of the former monolithic
# 03_qc_exploration/06_pcaPlot_Neha.r (2026-08-26). Sourced in order by that
# script, which remains the entry point. Runs at top level and shares the
# globals created by 06a_pca_core.r (mat, meta, pca, output_dir, helpers).
# Consumes um from 06b; treemap/chord use the optional-visualization audit trail.
# ====================================================================

if (!exists("mat") || !exists("meta") || !exists("pca") || !exists("output_dir")) {
  stop("PCA core state missing. Run 03_qc_exploration/06_pcaPlot_Neha.r, or source pca/06a_pca_core.r first.", call. = FALSE)
}

# AL) Arc diagram showing sample relationships based on PC distance
pc_arc_diagram <- function(group_key = "sample_class", n_connections = 50) {
    if (!requireNamespace("ggraph", quietly = TRUE) || !requireNamespace("igraph", quietly = TRUE)) {
        message("ggraph/igraph not available; skipping arc diagram")
        return(NULL)
    }
    
    # Calculate distances in PC space (first 5 PCs)
    pc_mat <- pca$x[, 1:min(5, ncol(pca$x)), drop = FALSE]
    dist_mat <- as.matrix(dist(pc_mat))
    
    # Get top N closest pairs
    dist_df <- reshape2::melt(dist_mat)
    dist_df <- dist_df[dist_df$Var1 != dist_df$Var2, ]
    dist_df <- dist_df[order(dist_df$value), ]
    dist_df <- head(dist_df, n_connections)
    
    # Create graph
    g <- igraph::graph_from_data_frame(dist_df[, 1:2], directed = FALSE)
    igraph::V(g)$group <- meta[igraph::V(g)$name, group_key]
    
    pal <- make_modern_palette(length(unique(igraph::V(g)$group)))
    
    p <- ggraph::ggraph(g, layout = "linear") +
        ggraph::geom_edge_arc(aes(alpha = ..index..), strength = 0.3, color = "gray70") +
        ggraph::geom_node_point(aes(color = group), size = 5) +
        scale_color_manual(values = pal) +
        theme_void() +
        theme(legend.position = "bottom") +
        labs(title = "Sample Similarity Network (Arc Diagram)", color = group_key)
    
    save_plot("plots/innovative", "pc_arc_diagram.svg", p, width = 14, height = 6)
}
pc_arc_diagram("sample_class", 50)

# AM) Radial/polar PCA plot
radial_pca_plot <- function(group_key = "sample_class") {
    df <- data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    # Convert to polar coordinates
    df$angle <- atan2(df$PC2, df$PC1)
    df$radius <- sqrt(df$PC1^2 + df$PC2^2)
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(df, aes(x = angle, y = radius, color = group)) +
        geom_point(size = 5, alpha = 0.8) +
        coord_polar(theta = "x") +
        scale_color_manual(values = pal) +
        theme_minimal() +
        theme(
            panel.grid.major = element_line(color = "#ECECEC"),
            axis.text.x = element_blank(),
            axis.title = element_blank()
        ) +
        labs(title = "Radial PCA Projection", color = group_key)
    
    save_plot("plots/innovative", "radial_pca.svg", p, width = 8, height = 8)
}
radial_pca_plot("sample_class")

# AN) Alluvial/Sankey diagram showing sample flow across PCA space
alluvial_pca_flow <- function(group_key = "sample_class") {
    if (!requireNamespace("ggalluvial", quietly = TRUE)) {
        message("ggalluvial not available; skipping alluvial plot")
        return(NULL)
    }
    
    # Bin PC1 and PC2 into quartiles
    df <- data.frame(
        sample = rownames(pca$x),
        PC1_bin = cut(pca$x[, 1], breaks = 4, labels = c("Low", "Med-Low", "Med-High", "High")),
        PC2_bin = cut(pca$x[, 2], breaks = 4, labels = c("Low", "Med-Low", "Med-High", "High")),
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    # Count flows
    flow_df <- as.data.frame(table(df$PC1_bin, df$PC2_bin, df$group))
    colnames(flow_df) <- c("PC1_bin", "PC2_bin", "group", "Freq")
    flow_df <- flow_df[flow_df$Freq > 0, ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(flow_df, aes(axis1 = PC1_bin, axis2 = PC2_bin, y = Freq, fill = group)) +
        ggalluvial::geom_alluvium(aes(fill = group), alpha = 0.7, width = 1/12) +
        ggalluvial::geom_flow(aes(fill = group), width = 1/12, alpha = 0.5) +
        scale_fill_manual(values = pal) +
        theme_minimal() +
        labs(title = "Sample Distribution Flow (PC1 → PC2)", 
             x = "PC Bins", 
             y = "Number of Samples",
             fill = group_key)
    
    save_plot("plots/innovative", "alluvial_pca_flow.svg", p, width = 10, height = 6)
}
alluvial_pca_flow("sample_class")

# AO) Ridgeline plot showing PC score distributions across groups
ridgeline_pc_distributions <- function(n_pcs = 3, group_key = "sample_class") {
    if (!requireNamespace("ggridges", quietly = TRUE)) {
        message("ggridges not available; skipping ridgeline plot")
        return(NULL)
    }
    
    # Prepare data for first N PCs
    pc_data_list <- list()
    for (i in 1:min(n_pcs, ncol(pca$x))) {
        pc_name <- paste0("PC", i)
        df <- data.frame(
            PC = pc_name,
            score = pca$x[, i],
            group = factor(meta[rownames(pca$x), group_key])
        )
        pc_data_list[[i]] <- df
    }
    combined_df <- do.call(rbind, pc_data_list)
    combined_df <- combined_df[!is.na(combined_df$group), ]
    
    pal <- make_modern_palette(nlevels(combined_df$group))
    
    p <- ggplot(combined_df, aes(x = score, y = PC, fill = group)) +
        ggridges::geom_density_ridges(alpha = 0.7, scale = 0.9) +
        scale_fill_manual(values = pal) +
        theme_minimal() +
        labs(title = "PC Score Distributions by Group", 
             x = "Score", y = "Principal Component", fill = group_key)
    
    save_plot("plots/innovative", "ridgeline_pc_distributions.svg", p, width = 10, height = 6)
}
ridgeline_pc_distributions(3, "sample_class")

# AP) Hexbin density plot with marginal distributions
hexbin_pca_with_marginals <- function(group_key = "sample_class") {
    if (!requireNamespace("ggExtra", quietly = TRUE)) {
        message("ggExtra not available; skipping marginal plot")
        return(NULL)
    }
    
    df <- data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(df, aes(PC1, PC2, color = group)) +
        geom_point(size = 4, alpha = 0.7) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "PCA with Marginal Densities", color = group_key)
    
    p_marg <- ggExtra::ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
    
    ggsave(file.path(ensure_dir(subdir("plots/innovative")), "pca_marginal_densities.png"),
           p_marg, width = 10, height = 8, dpi = 300)
}
hexbin_pca_with_marginals("sample_class")

# AQ) Network graph of PC loadings (protein-protein similarity)
loading_network_graph <- function(pc_x = 1, pc_y = 2, top_n = 50, cor_threshold = 0.7) {
    if (!requireNamespace("ggraph", quietly = TRUE) || !requireNamespace("igraph", quietly = TRUE)) {
        message("ggraph/igraph not available; skipping network graph")
        return(NULL)
    }
    
    # Get top proteins by loading magnitude
    rot <- as.data.frame(pca$rotation)
    mag <- sqrt(rot[, pc_x]^2 + rot[, pc_y]^2)
    top_proteins <- names(head(sort(mag, decreasing = TRUE), top_n))
    
    # Filter to proteins that exist in mat
    top_proteins <- intersect(top_proteins, rownames(mat))
    
    if (length(top_proteins) < 2) {
        message("Not enough valid proteins found for network graph")
        return(NULL)
    }
    
    # Calculate protein-protein correlations in expression space
    expr_subset <- mat[top_proteins, , drop = FALSE]
    cor_mat <- cor(t(expr_subset), use = "pairwise.complete.obs")
    
    # Create edges for high correlations
    cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
    edges <- which(abs(cor_mat) > cor_threshold, arr.ind = TRUE)
    
    if (nrow(edges) > 0) {
        edge_df <- data.frame(
            from = rownames(cor_mat)[edges[, 1]],
            to = colnames(cor_mat)[edges[, 2]],
            correlation = cor_mat[edges]
        )
        
        g <- igraph::graph_from_data_frame(edge_df, directed = FALSE)
        
        # Node attributes: loading values
        igraph::V(g)$PC1_load <- rot[igraph::V(g)$name, pc_x]
        igraph::V(g)$PC2_load <- rot[igraph::V(g)$name, pc_y]
        
        p <- ggraph::ggraph(g, layout = "fr") +
            ggraph::geom_edge_link(aes(alpha = abs(correlation)), color = "gray70") +
            ggraph::geom_node_point(aes(color = PC1_load, size = abs(PC2_load))) +
            ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
            scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
            theme_void() +
            labs(title = "Protein Co-expression Network (Top Contributors)",
                 color = paste0("PC", pc_x, " Loading"),
                 size = paste0("PC", pc_y, " |Loading|"))
        
        save_plot("plots/innovative", "loading_network_graph.svg", p, width = 12, height = 10)
    }
}
loading_network_graph(1, 2, 50, 0.7)

# AR) Parallel coordinates plot for top PCs
parallel_coordinates_pcs <- function(n_pcs = 5, group_key = "sample_class") {
    if (!requireNamespace("GGally", quietly = TRUE)) {
        message("GGally not available; skipping parallel coordinates")
        return(NULL)
    }
    
    pc_subset <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    df <- data.frame(
        pc_subset,
        group = factor(meta[rownames(pc_subset), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- GGally::ggparcoord(df, columns = 1:n_pcs, groupColumn = "group",
                           scale = "std", alphaLines = 0.6) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "Parallel Coordinates: Sample Profiles Across PCs", 
             color = group_key)
    
    save_plot("plots/innovative", "parallel_coordinates_pcs.svg", p, width = 10, height = 6)
}
parallel_coordinates_pcs(5, "sample_class")

# AS) Sunburst/treemap of variance contribution
variance_sunburst <- function() {
    if (!requireNamespace("treemap", quietly = TRUE) || !requireNamespace("d3r", quietly = TRUE)) {
        message("treemap/d3r not available; trying basic treemap")
        var_exp <- (pca$sdev^2) / sum(pca$sdev^2) * 100
        df <- prepare_neha_pca_variance_treemap_data(var_exp, max_pcs = 20L)
        output_path <- file.path(output_dir, "plots", "innovative", "variance_treemap.svg")
        result <- run_neha_pca_optional_plot(
            plot_name = "variance_treemap",
            usable_data = df,
            output_path = output_path,
            n_input_rows = attr(df, "n_input_rows"),
            render_fun = function(plot_data, output_path) {
                p <- ggplot(plot_data, aes(area = Variance, fill = Category, label = PC)) +
                    treemapify::geom_treemap() +
                    treemapify::geom_treemap_text(color = "white", place = "centre") +
                    scale_fill_manual(values = c("#6B5B95", "#88B04B", "#B2B2B2")) +
                    theme_minimal() +
                    labs(title = "Variance Contribution Treemap")

                save_plot("plots/innovative", basename(output_path), p, width = 10, height = 8)
            }
        )
        record_optional_visualization(result)
        return(NULL)
    }
}
variance_sunburst()

# AT) Contour plot showing density in PC space
contour_density_pca <- function(group_key = "sample_class") {
    df <- data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(df, aes(PC1, PC2, color = group)) +
        geom_density_2d(size = 1) +
        geom_point(alpha = 0.6, size = 3) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "PCA with Density Contours", color = group_key)
    
    save_plot("plots/innovative", "contour_density_pca.svg", p, width = 10, height = 8)
}
contour_density_pca("sample_class")

# AU) Chord diagram for group relationships based on PC overlap
chord_diagram_groups <- function(group_key = "sample_class", pc_threshold = 2) {
    if (!requireNamespace("circlize", quietly = TRUE)) {
        message("circlize not available; skipping chord diagram")
        return(NULL)
    }
    
    # Calculate group centroids from sample IDs, not names(factor), which is NULL.
    centroids <- prepare_neha_pca_group_centroids(pca$x, meta, group_key, n_pcs = 5L)
    output_path <- file.path(output_dir, "plots", "innovative", "chord_diagram_groups.png")
    result <- run_neha_pca_optional_plot(
        plot_name = "chord_diagram_groups",
        usable_data = centroids,
        output_path = output_path,
        n_input_rows = attr(centroids, "n_input_rows"),
        render_fun = function(plot_data, output_path) {
            rownames(plot_data) <- plot_data$group
            plot_data$group <- NULL

            # Calculate pairwise distances between group centroids.
            dist_mat <- as.matrix(dist(plot_data))

            # Convert to flow matrix (inverse of distance).
            flow_mat <- 1 / (1 + dist_mat)
            diag(flow_mat) <- 0

            png(output_path, width = 10, height = 10, units = "in", res = 300)
            on.exit({
                if (grDevices::dev.cur() > 1L) grDevices::dev.off()
                circlize::circos.clear()
            }, add = TRUE)
            circlize::chordDiagram(
                flow_mat,
                grid.col = make_modern_palette(nrow(flow_mat)),
                transparency = 0.5
            )
            title("Group Relationships in PC Space")
        }
    )
    record_optional_visualization(result)
}
chord_diagram_groups("sample_class")

# AV) Small multiples: individual sample trajectories
sample_trajectory_multiples <- function(n_samples = 16, group_key = "sample_class") {
    # Plot all samples, faceted by group
    groups <- factor(meta[rownames(pca$x), group_key])
    
    df <- data.frame(
        sample_id = rownames(pca$x),
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = groups,
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$group), ]
    
    pal <- make_modern_palette(nlevels(droplevels(df$group)))
    
    p <- ggplot(df, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 3, alpha = 0.8) +
        facet_wrap(~ group, ncol = 3) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        theme(strip.text = element_text(size = 10, face = "bold")) +
        labs(title = "Sample Positions in PC Space by Group", color = group_key)
    
    save_plot("plots/innovative", "sample_multiples.svg", p, width = 12, height = 10)
}
sample_trajectory_multiples(16, "sample_class")

message("All innovative visualizations completed!")

# hexbin_umap_with_marginals - UMAP version with marginal distributions
hexbin_umap_with_marginals <- function(group_key = "sample_class") {
    if (is.null(um)) {
        message("UMAP not available; skipping UMAP marginal plot")
        return(NULL)
    }
    if (!requireNamespace("ggExtra", quietly = TRUE)) {
        message("ggExtra not available; skipping marginal plot")
        return(NULL)
    }
    
    # Recreate um_df to ensure consistency with original UMAP plots
    um_df <- data.frame(UMAP1 = um[,1], UMAP2 = um[,2], meta[rownames(um), , drop = FALSE])
    
    # Filter by group_key
    if (!group_key %in% names(um_df)) {
        message(paste("Group key", group_key, "not found in metadata"))
        return(NULL)
    }
    
    grp <- factor(trim_ws(as.character(um_df[[group_key]])))
    um_df$group <- grp
    um_df <- um_df[!is.na(um_df$group), ]
    
    pal <- make_modern_palette(nlevels(droplevels(um_df$group)))
    
    p <- ggplot(um_df, aes(UMAP1, UMAP2, color = group)) +
        geom_point(size = 6, shape = 16, alpha = 0.8, stroke = 0) +
        scale_color_manual(values = pal, na.translate = FALSE) +
        theme_pca_min() +
        theme(legend.position = "bottom") +
        labs(title = paste("UMAP by", group_key, "with Marginal Densities"), color = group_key)
    
    p_marg <- ggExtra::ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE, size = 10)
    
    ggsave(file.path(ensure_dir(subdir("plots/umap")), 
                     paste0("umap_marginal_densities_by_", group_key, ".svg")),
           p_marg, width = 10, height = 8, dpi = 300)
}
hexbin_umap_with_marginals("sample_class")
hexbin_umap_with_marginals("condition")
hexbin_umap_with_marginals("AnimalID")
