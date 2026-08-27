# ====================================================================
# Statistical summary barplots (PCA and UMAP dimensions), PC1 density, multi-panel PC combinations
# Part of the PCA workflow split out of the former monolithic
# 03_qc_exploration/06_pcaPlot_Neha.r (2026-08-26). Sourced in order by that
# script, which remains the entry point. Runs at top level and shares the
# globals created by 06a_pca_core.r (mat, meta, pca, output_dir, helpers).
# Consumes um from 06b.
# ====================================================================

if (!exists("mat") || !exists("meta") || !exists("pca") || !exists("output_dir")) {
  stop("PCA core state missing. Run 03_qc_exploration/06_pcaPlot_Neha.r, or source pca/06a_pca_core.r first.", call. = FALSE)
}

# AI) Statistical summary barplot (ANOVA/Kruskal results)
statistical_summary_barplot <- function(group_key = "sample_class") {
    if (!group_key %in% names(meta)) return(NULL)
    
    # Run tests for multiple PCs
    results <- list()
    for (i in 1:min(10, ncol(pca$x))) {
        pc_name <- paste0("PC", i)
        df <- data.frame(
            score = pca$x[, i],
            group = factor(meta[rownames(pca$x), group_key])
        )
        df <- df[!is.na(df$group), ]
        
        # Check normality
        normality <- by(df$score, df$group, function(x) {
            if (length(x) < 3) return(list(p = NA))
            shapiro.test(x)
        })
        all_normal <- all(sapply(normality, function(x) x$p.value > 0.05), na.rm = TRUE)
        
        # Run appropriate test
        if (all_normal) {
            test_result <- summary(aov(score ~ group, data = df))
            p_val <- test_result[[1]][["Pr(>F)"]][1]
            test_name <- "ANOVA"
        } else {
            test_result <- kruskal.test(score ~ group, data = df)
            p_val <- test_result$p.value
            test_name <- "Kruskal-Wallis"
        }
        
        results[[i]] <- data.frame(
            PC = pc_name,
            P_Value = p_val,
            NegLog10P = -log10(p_val),
            Test = test_name,
            Significant = p_val < 0.05
        )
    }
    
    results_df <- do.call(rbind, results)
    
    p <- ggplot(results_df, aes(x = PC, y = NegLog10P, fill = Significant)) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
                  color = "red", size = 1) +
        annotate("text", x = 2, y = -log10(0.05) + 0.3, 
                label = "p = 0.05", color = "red", size = 4) +
        scale_fill_manual(values = c("FALSE" = "#B2B2B2", "TRUE" = "#6B5B95"),
                         labels = c("ns", "p < 0.05")) +
        theme_pca_min() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("PC Association with", group_key),
             x = "Principal Component",
             y = "-log10(p-value)",
             fill = "Significance")
    
    save_plot("plots/main_figure", paste0("pc_stats_summary_", group_key, ".svg"), 
             p, width = 5, height = 4)
    
    # Save table
    save_table("tables/statistics", 
              paste0("pc_group_associations_", group_key, ".csv"), 
              results_df)
}
statistical_summary_barplot("sample_class")
statistical_summary_barplot("condition_code")

# AI-UMAP) Statistical summary barplot for UMAP dimensions
statistical_summary_barplot_umap <- function(group_key = "sample_class") {
    if (!group_key %in% names(meta)) return(NULL)
    if (is.null(um)) {
        message("UMAP not available; skipping UMAP statistical summary")
        return(NULL)
    }
    
    # Run tests for all UMAP dimensions
    results <- list()
    n_dims <- ncol(um)
    
    for (i in 1:n_dims) {
        dim_name <- paste0("UMAP", i)
        df <- data.frame(
            score = um[, i],
            group = factor(meta[rownames(um), group_key])
        )
        df <- df[!is.na(df$group), ]
        
        # Check normality
        normality <- by(df$score, df$group, function(x) {
            if (length(x) < 3) return(list(p = NA))
            shapiro.test(x)
        })
        all_normal <- all(sapply(normality, function(x) x$p.value > 0.05), na.rm = TRUE)
        
        # Run appropriate test
        if (all_normal) {
            test_result <- summary(aov(score ~ group, data = df))
            p_val <- test_result[[1]][["Pr(>F)"]][1]
            test_name <- "ANOVA"
        } else {
            test_result <- kruskal.test(score ~ group, data = df)
            p_val <- test_result$p.value
            test_name <- "Kruskal-Wallis"
        }
        
        results[[i]] <- data.frame(
            Dimension = dim_name,
            P_Value = p_val,
            NegLog10P = -log10(p_val),
            Test = test_name,
            Significant = p_val < 0.05,
            stringsAsFactors = FALSE
        )
    }
    
    results_df <- do.call(rbind, results)
    
    # Create factor with proper ordering
    results_df$Dimension <- factor(results_df$Dimension, 
                                   levels = paste0("UMAP", 1:n_dims))
    
    p <- ggplot(results_df, aes(x = Dimension, y = NegLog10P, fill = Significant)) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
                  color = "red", size = 1) +
        annotate("text", x = min(2, n_dims/2), y = -log10(0.05) + 0.3, 
                label = "p = 0.05", color = "red", size = 4) +
        scale_fill_manual(values = c("FALSE" = "#B2B2B2", "TRUE" = "#6B5B95"),
                         labels = c("ns", "p < 0.05")) +
        theme_pca_min() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("UMAP Dimension Association with", group_key),
             x = "UMAP Dimension",
             y = "-log10(p-value)",
             fill = "Significance")
    
    save_plot("plots/main_figure", paste0("umap_stats_summary_", group_key, ".svg"), 
             p, width = max(5, n_dims * 0.5), height = 4)
    
    # Save table
    save_table("tables/statistics", 
              paste0("umap_group_associations_", group_key, ".csv"), 
              results_df)
    
    message(paste("Completed UMAP statistical summary for", n_dims, "dimensions by", group_key))
}
statistical_summary_barplot_umap("sample_class")
statistical_summary_barplot_umap("condition_code")

# AJ) Sample distribution along PC1 (density + rug plot)
pc1_distribution_main <- function(group_key = "sample_class") {
    if (!group_key %in% names(meta)) return(NULL)
    
    df <- data.frame(
        PC1 = pca$x[, 1],
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(df, aes(x = PC1, fill = group, color = group)) +
        geom_density(alpha = 0.5, size = 1.2) +
        geom_rug(aes(color = group), alpha = 0.8, size = 1) +
        scale_fill_manual(values = pal) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        theme(legend.position = "top") +
        labs(title = "PC1 Distribution by Group",
             x = "PC1 Score",
             y = "Density",
             fill = group_key,
             color = group_key)
    
    save_plot("plots/main_figure", paste0("pc1_distribution_", group_key, ".svg"), 
             p, width = 10, height = 6)
}
pc1_distribution_main("sample_class")

# AK) Multi-panel: PC1 vs PC2, PC2 vs PC3, PC1 vs PC3
multipanel_pca_combinations <- function(group_key = "sample_class") {
    if (!group_key %in% names(meta)) return(NULL)
    
    grp <- build_group(group_key)
    keep <- !is.na(grp)
    grp2 <- droplevels(grp[keep])
    pal <- make_modern_palette(nlevels(grp2))
    
    varp <- (pca$sdev^2)/sum(pca$sdev^2) * 100
    
    # PC1 vs PC2
    df12 <- data.frame(pca$x[keep, c(1,2), drop=FALSE], group = grp2)
    p12 <- ggplot(df12, aes(PC1, PC2, color = group)) +
        geom_point(size = 4, alpha = 0.8) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(x = sprintf("PC1 (%.1f%%)", varp[1]),
             y = sprintf("PC2 (%.1f%%)", varp[2]),
             title = "PC1 vs PC2")
    
    # PC2 vs PC3
    if (ncol(pca$x) >= 3) {
        df23 <- data.frame(pca$x[keep, c(2,3), drop=FALSE], group = grp2)
        p23 <- ggplot(df23, aes(PC2, PC3, color = group)) +
            geom_point(size = 4, alpha = 0.8) +
            scale_color_manual(values = pal) +
            theme_pca_min() +
            labs(x = sprintf("PC2 (%.1f%%)", varp[2]),
                 y = sprintf("PC3 (%.1f%%)", varp[3]),
                 title = "PC2 vs PC3")
        
        # PC1 vs PC3
        df13 <- data.frame(pca$x[keep, c(1,3), drop=FALSE], group = grp2)
        p13 <- ggplot(df13, aes(PC1, PC3, color = group)) +
            geom_point(size = 4, alpha = 0.8) +
            scale_color_manual(values = pal) +
            theme_pca_min() +
            labs(x = sprintf("PC1 (%.1f%%)", varp[1]),
                 y = sprintf("PC3 (%.1f%%)", varp[3]),
                 title = "PC1 vs PC3")
        
        if (requireNamespace("gridExtra", quietly = TRUE)) {
            combined <- gridExtra::grid.arrange(p12, p23, p13, ncol = 3)
            ggsave(file.path(ensure_dir(subdir("plots/main_figure")), 
                           paste0("pca_combinations_", group_key, ".svg")),
                   combined, width = 18, height = 6, dpi = 300)
        }
    }
}
multipanel_pca_combinations("sample_class")

message("All main figure plots completed!")


# ================== Innovative Main Figure Visualizations ==================

