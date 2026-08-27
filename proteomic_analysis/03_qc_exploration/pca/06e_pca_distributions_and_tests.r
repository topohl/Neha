# ====================================================================
# PC boxplots, violins with statistical tests, loading contribution bars, correlation circles
# Part of the PCA workflow split out of the former monolithic
# 03_qc_exploration/06_pcaPlot_Neha.r (2026-08-26). Sourced in order by that
# script, which remains the entry point. Runs at top level and shares the
# globals created by 06a_pca_core.r (mat, meta, pca, output_dir, helpers).
# ====================================================================

if (!exists("mat") || !exists("meta") || !exists("pca") || !exists("output_dir")) {
  stop("PCA core state missing. Run 03_qc_exploration/06_pcaPlot_Neha.r, or source pca/06a_pca_core.r first.", call. = FALSE)
}

# ================== Additional Main & Supplementary Figure Plots ==================

# U) Boxplots of PC scores by experimental groups
pc_boxplots_by_group <- function(n_pcs = 5, group_key = "condition_code") {
    if (!group_key %in% names(meta)) {
        message(paste(group_key, "not in metadata; skipping PC boxplots"))
        return(NULL)
    }
    
    for (i in 1:min(n_pcs, ncol(pca$x))) {
        pc_name <- paste0("PC", i)
        df <- data.frame(
            score = pca$x[, i],
            group = factor(meta[rownames(pca$x), group_key])
        )
        df <- df[!is.na(df$group), ]
        
        pal <- make_modern_palette(nlevels(df$group))
        
        p <- ggplot(df, aes(x = group, y = score, fill = group)) +
            geom_boxplot(alpha = 0.7, outlier.shape = NA) +
            geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
            scale_fill_manual(values = pal) +
            theme_pca_min() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = paste(pc_name, "Scores by", group_key),
                 x = group_key, y = paste(pc_name, "Score"))
        
        save_plot("plots/boxplots", paste0(pc_name, "_boxplot_by_", group_key, ".svg"), p)
    }
}
pc_boxplots_by_group(5, "condition_code")
pc_boxplots_by_group(5, "sample_class")

# V) Violin plots with overlaid points and statistical tests
pc_violin_plots <- function(n_pcs = 5, group_key = "sample_class") {
    if (!group_key %in% names(meta)) return(NULL)
    
    all_normality <- list()
    all_main_tests <- list()
    all_posthocs <- list()
    
    for (i in 1:min(n_pcs, ncol(pca$x))) {
        pc_name <- paste0("PC", i)
        df <- data.frame(
            score = pca$x[, i],
            group = factor(meta[rownames(pca$x), group_key])
        )
        df <- df[!is.na(df$group), ]
        
        cat("\n=== Processing", pc_name, "for", group_key, "===\n")
        
        pal <- make_modern_palette(nlevels(df$group))
        
        # Create base plot
        p <- ggplot(df, aes(x = group, y = score, fill = group, color = group)) +
            geom_hline(yintercept = 0, linetype = "solid", color = "lightgrey", size = 1) + 
            geom_violin(alpha = 0.4, trim = FALSE, show.legend = FALSE, linetype = "blank") +
            geom_jitter(width = 0.15, alpha = 0.8, size = 4) +
            scale_fill_manual(values = pal) +
            scale_color_manual(values = pal) +
            theme_pca_min() + 
            theme(panel.grid = element_blank(),
                  axis.text.x = element_text(hjust = 0.5),
                  legend.position = "none") +
            labs(title = paste(pc_name, "Distribution by", group_key),
                 x = group_key, y = paste(pc_name, "Score"))
            
        
        if (length(unique(df$group)) > 1) {
            p <- p + stat_summary(fun = median, geom = "crossbar", width = 0.5, 
                                 color = "black", size = 0.4, fatten = 1)
        }
        
        # Statistical testing - normality
        normality_by_group <- by(df$score, df$group, function(x) {
            if (length(x) < 3) return(data.frame(W = NA, p = NA))
            test <- shapiro.test(x)
            data.frame(W = test$statistic, p = test$p.value)
        })
        
        normality_df <- do.call(rbind, lapply(names(normality_by_group), function(g) {
            data.frame(
                PC = pc_name,
                Group = g,
                Test = "Shapiro-Wilk",
                Statistic = normality_by_group[[g]]$W,
                P_Value = normality_by_group[[g]]$p,
                stringsAsFactors = FALSE
            )
        }))
        
        all_normality[[pc_name]] <- normality_df
        
        all_normal <- all(normality_df$P_Value > 0.05, na.rm = TRUE)
        n_groups <- nlevels(df$group)
        
        # Main test
        main_test_df <- NULL
        posthoc_df <- NULL
        
        if (n_groups >= 2) {
            if (all_normal) {
                anova_result <- aov(score ~ group, data = df)
                anova_summary <- summary(anova_result)
                
                main_test_df <- data.frame(
                    PC = pc_name,
                    Grouping = group_key,
                    Test = "ANOVA",
                    N_Groups = n_groups,
                    All_Normal = all_normal,
                    Statistic = anova_summary[[1]][["F value"]][1],
                    P_Value = anova_summary[[1]][["Pr(>F)"]][1],
                    stringsAsFactors = FALSE
                )
                
                if (main_test_df$P_Value < 0.05 && n_groups > 2) {
                    posthoc <- TukeyHSD(anova_result)
                    posthoc_df <- as.data.frame(posthoc$group)
                    posthoc_df$Comparison <- rownames(posthoc_df)
                    posthoc_df$PC <- pc_name
                    posthoc_df$Grouping <- group_key
                    posthoc_df$Test <- "Tukey HSD"
                    posthoc_df <- posthoc_df[, c("PC", "Grouping", "Test", "Comparison", "diff", "lwr", "upr", "p adj")]
                    colnames(posthoc_df) <- c("PC", "Grouping", "Test", "Comparison", "Diff", "Lower_CI", "Upper_CI", "P_Value")
                }
                
            } else {
                kw_result <- kruskal.test(score ~ group, data = df);
                
                main_test_df <- data.frame(
                    PC = pc_name,
                    Grouping = group_key,
                    Test = "Kruskal-Wallis",
                    N_Groups = n_groups,
                    All_Normal = all_normal,
                    Statistic = kw_result$statistic,
                    P_Value = kw_result$p.value,
                    stringsAsFactors = FALSE
                );
                
                if (main_test_df$P_Value < 0.05 && n_groups > 2) {
                    pw_result <- pairwise.wilcox.test(df$score, df$group, 
                                                      p.adjust.method = "bonferroni");
                    
                    pw_matrix <- pw_result$p.value;
                    comparisons <- which(!is.na(pw_matrix), arr.ind = TRUE);
                    
                    if (nrow(comparisons) > 0) {
                        posthoc_df <- data.frame(
                            PC = pc_name,
                            Grouping = group_key,
                            Test = "Pairwise Wilcoxon",
                            Comparison = paste(rownames(pw_matrix)[comparisons[,1]], 
                                             colnames(pw_matrix)[comparisons[,2]], 
                                             sep = " vs "),
                            Diff = NA_real_,
                            Lower_CI = NA_real_,
                            Upper_CI = NA_real_,
                            P_Value = pw_matrix[comparisons],
                            stringsAsFactors = FALSE
                        );
                    }
                }
            }
        }
        
        # Calculate y-axis limits BEFORE adding significance annotations
        y_range <- range(df$score, na.rm = TRUE)
        y_span <- diff(y_range)
        y_min <- y_range[1] - y_span * 0.15;  # Bottom padding
        y_start <- y_range[2] + y_span * 0.15;  # Starting position for brackets
        y_step <- y_span * 0.25;  # Step size between brackets
        
        # Count number of significant comparisons to calculate required space
        n_sig_comparisons <- 0;
        if (!is.null(posthoc_df) && nrow(posthoc_df) > 0) {
            n_sig_comparisons <- sum(posthoc_df$P_Value < 0.05);
        }
        
        # Add significance annotations
        max_y <- y_start;
        if (n_sig_comparisons > 0) {
            posthoc_df$group1 <- sapply(strsplit(posthoc_df$Comparison, " vs "), `[`, 1);
            posthoc_df$replicate_unit <- sapply(strsplit(posthoc_df$Comparison, " vs "), `[`, 2);
            
            sig_comparisons <- posthoc_df[posthoc_df$P_Value < 0.05, ];
            
            for (idx in 1:nrow(sig_comparisons)) {
                y_pos <- y_start + (idx - 1) * y_step;
                max_y <- max(max_y, y_pos);
                
                p_val <- sig_comparisons$P_Value[idx];
                sig_symbol <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "ns";
                
                x1 <- which(levels(df$group) == sig_comparisons$group1[idx]);
                x2 <- which(levels(df$group) == sig_comparisons$replicate_unit[idx]);
                
                p <- p + 
                    annotate("segment", x = x1, xend = x2, y = y_pos, yend = y_pos, color = "black", size = 0.5) +
                    annotate("segment", x = x1, xend = x1, y = y_pos, yend = y_pos - y_span * 0.02, color = "black", size = 0.5) +
                    annotate("segment", x = x2, xend = x2, y = y_pos, yend = y_pos - y_span * 0.02, color = "black", size = 0.5) +
                    annotate("text", x = (x1 + x2) / 2, y = y_pos + y_span * 0.04, label = sig_symbol, size = 5, fontface = "bold");
            }
        }
        
        # Set y-axis limits with generous padding
        # Add extra space above the highest bracket plus text
        y_max <- max_y + y_span * 0.25;  # Generous top padding
        p <- p + coord_cartesian(ylim = c(y_min, y_max), clip = "off");
        
        # Save plot
        save_plot("plots/violin", paste0(pc_name, "_violin_by_", group_key, ".svg"), p, width = 5, height = 5);
        
        if (!is.null(main_test_df)) all_main_tests[[pc_name]] <- main_test_df;
        if (!is.null(posthoc_df)) all_posthocs[[pc_name]] <- posthoc_df;
    }
    
    # Save summary tables
    if (length(all_main_tests) > 0) {
        combined_main_tests <- do.call(rbind, all_main_tests);
        combined_main_tests$Significant <- ifelse(combined_main_tests$P_Value < 0.05, "Yes", "No");
        save_table("tables/statistics", paste0("summary_main_tests_", group_key, ".csv"), combined_main_tests);
    }
    
    if (length(all_posthocs) > 0) {
        combined_posthocs <- do.call(rbind, all_posthocs);
        combined_posthocs$Significant <- ifelse(combined_posthocs$P_Value < 0.05, "Yes", "No");
        save_table("tables/statistics", paste0("summary_posthoc_tests_", group_key, ".csv"), combined_posthocs);
    }
    
    if (length(all_normality) > 0) {
        combined_normality <- do.call(rbind, all_normality);
        save_table("tables/statistics", paste0("summary_normality_tests_", group_key, ".csv"), combined_normality);
    }
}
pc_violin_plots(5, "sample_class")
pc_violin_plots(5, "condition_code")

# W) Loading contribution per PC (top positive and negative)
loading_contribution_bars <- function(n_pcs = 5, top_n = 10) {
    for (i in 1:min(n_pcs, ncol(pca$rotation))) {
        pc_name <- paste0("PC", i)
        loadings <- pca$rotation[, i]
        
        # Get top positive and negative
        top_pos <- head(sort(loadings, decreasing = TRUE), top_n)
        top_neg <- head(sort(loadings, decreasing = FALSE), top_n)
        
        df <- data.frame(
            protein = c(names(top_pos), names(top_neg)),
            loading = c(top_pos, top_neg),
            direction = c(rep("Positive", top_n), rep("Negative", top_n))
        )
        df$protein <- factor(df$protein, levels = df$protein[order(df$loading)])
        
        p <- ggplot(df, aes(x = protein, y = loading, fill = direction)) +
            geom_col() +
            coord_flip() +
            scale_fill_manual(values = c("Positive" = "#FFB5B5", "Negative" = "#B5D4FF")) +
            theme_pca_min() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            labs(title = paste("Top Loadings for", pc_name),
                 x = "Protein", y = "Loading Value", fill = "Direction")
        
        save_plot("plots/loadings", paste0(pc_name, "_top_loadings_bar.svg"), p, width = 6, height = 4)
    }
}
loading_contribution_bars(5, 5)

# X) Correlation circle plot (biplot style, loadings only)
correlation_circle_plot <- function(pc_x = 1, pc_y = 2, top_n = 30) {
    rot <- as.data.frame(pca$rotation)
    
    # Select top contributors to both PCs
    mag <- sqrt(rot[, pc_x]^2 + rot[, pc_y]^2)
    top_idx <- order(mag, decreasing = TRUE)[1:min(top_n, length(mag))]
    
    df <- data.frame(
        protein = rownames(rot)[top_idx],
        PC1 = rot[top_idx, pc_x],
        PC2 = rot[top_idx, pc_y]
    )
    
    # Create circle
    theta <- seq(0, 2*pi, length.out = 100)
    circle <- data.frame(x = cos(theta), y = sin(theta))
    
    p <- ggplot() +
        geom_path(data = circle, aes(x, y), color = "gray70", linetype = "dashed") +
        geom_segment(data = df, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                    arrow = arrow(length = unit(0.2, "cm")), alpha = 0.6, color = "#6B5B95") +
        geom_text_repel(data = df, aes(x = PC1, y = PC2, label = protein),
                       size = 3, max.overlaps = 30) +
        coord_fixed() +
        theme_pca_min() +
        labs(title = paste("Correlation Circle (PC", pc_x, "vs PC", pc_y, ")"),
             x = paste0("PC", pc_x, " Loadings"),
             y = paste0("PC", pc_y, " Loadings"))
    
    save_plot("plots/biplot", paste0("correlation_circle_PC", pc_x, "_PC", pc_y, ".svg"), p)
}
correlation_circle_plot(1, 2, 30)
correlation_circle_plot(2, 3, 30)

