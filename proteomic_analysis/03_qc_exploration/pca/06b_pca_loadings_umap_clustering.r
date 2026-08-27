# ====================================================================
# Loadings, PC-protein correlations, UMAP, outlier detection, clustering, batch check
# Part of the PCA workflow split out of the former monolithic
# 03_qc_exploration/06_pcaPlot_Neha.r (2026-08-26). Sourced in order by that
# script, which remains the entry point. Runs at top level and shares the
# globals created by 06a_pca_core.r (mat, meta, pca, output_dir, helpers).
# Creates rot/top_n, cors_df, um/npc consumed by later parts.
# ====================================================================

if (!exists("mat") || !exists("meta") || !exists("pca") || !exists("output_dir")) {
  stop("PCA core state missing. Run 03_qc_exploration/06_pcaPlot_Neha.r, or source pca/06a_pca_core.r first.", call. = FALSE)
}


# 2) PC ~ metadata ANOVA with effect sizes (sample_class and condition_code), robust and always writes a file
pc_df <- as.data.frame(pca$x)
pc_df$sample <- rownames(pc_df)
pc_meta <- cbind(pc_df[, "sample", drop = FALSE], 
                 pc_df[, grep("^PC", names(pc_df))], 
                 meta[rownames(pc_df), , drop = FALSE])

# Candidate factors
cand_vars <- intersect(c("sample_class","condition_code"), names(pc_meta))
if (length(cand_vars)) {
    pc_meta[cand_vars] <- lapply(pc_meta[cand_vars], function(x) {
        x <- trim_ws(as.character(x)); x[x == ""] <- NA; factor(x)
    })
}

anova_rows <- list()
for (pc_name in colnames(pca$x)) {
    if (!length(cand_vars)) break
    dat <- pc_meta[, c(pc_name, cand_vars), drop = FALSE]
    dat <- dat[complete.cases(dat), , drop = FALSE]
    if (nrow(dat) < 3) next
    vars_ok <- cand_vars[sapply(cand_vars, function(v) nlevels(droplevels(dat[[v]])) >= 2)]
    if (!length(vars_ok)) next
    form <- as.formula(sprintf("%s ~ %s", pc_name, paste(vars_ok, collapse = " + ")))
    fit <- try(aov(form, data = dat), silent = TRUE)
    if (inherits(fit, "try-error")) next
    a <- try(summary(fit)[[1]], silent = TRUE)
    if (inherits(a, "try-error")) next
    ss_total <- sum(a[, "Sum Sq"], na.rm = TRUE)
    for (v in vars_ok) {
        if (v %in% rownames(a)) {
            ss <- a[v, "Sum Sq"]; df <- a[v, "Df"]; ms <- a[v, "Mean Sq"]; f <- a[v, "F value"]; p <- a[v, "Pr(>F)"]
            eta2 <- if (is.finite(ss_total) && ss_total > 0) ss / ss_total else NA_real_
            anova_rows[[length(anova_rows)+1L]] <- data.frame(
                PC = pc_name, term = v, df = as.numeric(df), ss = as.numeric(ss),
                ms = as.numeric(ms), F = as.numeric(f), p = as.numeric(p), eta2 = as.numeric(eta2),
                stringsAsFactors = FALSE
            )
        }
    }
}
anova_tab <- if (length(anova_rows)) do.call(rbind, anova_rows) else
    data.frame(PC=character(), term=character(), df=double(), ss=double(), ms=double(), F=double(), p=double(), eta2=double(), stringsAsFactors = FALSE)

if (nrow(anova_tab)) anova_tab$q <- p.adjust(anova_tab$p, method = "BH") else anova_tab$q <- numeric(0)

# Always write an output file (even if empty)
save_table("tables/associations", "pc_meta_anova.csv", anova_tab)

# 3) Top loadings export for PC1/PC2 (+/-) and enrichment-ready lists
rot <- as.data.frame(pca$rotation)
rot$protein <- rownames(pca$rotation)
top_n <- 50
for (k in c("PC1","PC2")) {
  if (!k %in% colnames(rot)) next
  ord_pos <- rot[order(-rot[[k]]), c("protein", k)]
  ord_neg <- rot[order(rot[[k]]),  c("protein", k)]
  save_table("tables/loadings", sprintf("loadings_%s_toppos.csv", k), head(ord_pos, top_n))
  save_table("tables/loadings", sprintf("loadings_%s_topneg.csv", k), head(ord_neg, top_n))
  ranked <- rot[, c("protein", k)]
  ranked <- ranked[order(-ranked[[k]]), ]
  save_table("tables/loadings", sprintf("loadings_%s_ranked_full.csv", k), ranked)
}

# 4) PC–protein correlations (per protein vs PC1/PC2)
pc_scores <- pca$x[, colnames(pca$x)[1:min(2, ncol(pca$x))], drop = FALSE]
cors <- lapply(colnames(pc_scores), function(pc) {
  s <- pc_scores[, pc]
  ct <- apply(mat, 1, function(v) {
    ok <- is.finite(v) & is.finite(s)
    if (sum(ok) < 3) return(c(r=NA_real_, p=NA_real_))
    test <- suppressWarnings(cor.test(v[ok], s[ok], method = "pearson"))
    c(r = unname(test$estimate), p = unname(test$p.value))
  })
  df <- as.data.frame(t(ct))
  df$protein <- rownames(df)
  df$pc <- pc
  df
})
cors_df <- do.call(rbind, cors)
cors_df$q <- ave(cors_df$p, cors_df$pc, FUN = function(x) p.adjust(x, method="BH"))
save_table("tables/correlations", "protein_pc_correlations.csv", cors_df)

# 5) Signed loading heatmap for PC1/PC2 (top 50 +/−)
heat_dir <- ensure_dir(subdir("plots/heatmaps"))
for (k in c("PC1","PC2")) {
  if (!k %in% colnames(rot)) next
  pos <- head(rot[order(-rot[[k]]), "protein"], top_n)
  neg <- head(rot[order(rot[[k]]),  "protein"], top_n)
  sel <- unique(c(pos, neg))
  sel <- sel[sel %in% rownames(mat)]
  if (length(sel) >= 4) {
    ord_samples <- order(pca$x[,k])
    M <- mat[sel, ord_samples, drop = FALSE]
    save_pheatmap(M, file.path(heat_dir, sprintf("heatmap_%s_topposneg.png", k)), width=8, height=10, scale="row")
  }
}

# 6) UMAP on PCs (first 20 PCs)
run_umap <- function(X, n_neighbors = 15, min_dist = 0.2, n_components = 2, metric = "euclidean"){
  if (!requireNamespace("uwot", quietly = TRUE)) {
    message("uwot not installed; skipping UMAP.")
    return(NULL)
  }
  set.seed(42)
  uwot::umap(X, n_neighbors = n_neighbors, min_dist = min_dist, n_components = n_components, metric = metric, verbose = FALSE)
}
npc <- min(20, ncol(pca$x))
um <- run_umap(pca$x[, 1:npc, drop = FALSE], n_neighbors = 15, min_dist = 0.2)
if (!is.null(um)) {
  um_df <- data.frame(UMAP1 = um[,1], UMAP2 = um[,2], meta[rownames(pca$x), , drop = FALSE])
  plot_umap_group <- function(key, out_file){
    if (!key %in% names(um_df)) return(NULL)
    grp <- factor(trim_ws(as.character(um_df[[key]])))
    pal <- make_modern_palette(nlevels(droplevels(grp)))
    p <- ggplot(um_df, aes(UMAP1, UMAP2, color = grp)) +
      geom_point(shape = 16, alpha = 0.8, size = 6, stroke = 0) +
      scale_color_manual(values = pal, na.translate = FALSE) +
      theme_pca_min() + labs(title = paste("UMAP by", key), color = key)
    save_plot("plots/umap", out_file, p)
    invisible(p)
  }
  plot_umap_group("sample_class", "umap_by_sample_class.svg")
  plot_umap_group("condition", "umap_by_condition.svg")
  plot_umap_group("AnimalID", "umap_by_animal_id.svg")
}
# detect switched samples / outliers based on umap
if (!is.null(um)) {
    # Calculate sample_class group centers in UMAP space using median
    sample_class_centers <- aggregate(um, by = list(sample_class = meta[rownames(um), "sample_class"]), FUN = median)
    rownames(sample_class_centers) <- sample_class_centers$sample_class
    sample_class_centers$sample_class <- NULL
    
    # Calculate distance of each sample to its own sample_class center
    sample_distances <- sapply(rownames(um), function(sample) {
        sample_sample_class <- meta[sample, "sample_class"]
        if (is.na(sample_sample_class)) return(NA_real_)
        
        center <- as.numeric(sample_class_centers[sample_sample_class, ])
        sample_coords <- as.numeric(um[sample, ])
        
        # Euclidean distance to center
        sqrt(sum((sample_coords - center)^2))
    })
    
    # Identify outliers per sample_class (samples far from their own group center)
    outlier_samples <- character(0)
    outlier_info <- list()
    
    for (ct in unique(meta$sample_class)) {
        if (is.na(ct)) next
        
        samples_in_sample_class <- rownames(meta)[meta$sample_class == ct & !is.na(meta$sample_class)]
        samples_in_sample_class <- intersect(samples_in_sample_class, names(sample_distances))
        
        if (length(samples_in_sample_class) < 3) next  # Need at least 3 samples to detect outliers
        
        dists <- sample_distances[samples_in_sample_class]
        
        # Use 95th percentile or 2 MADs as threshold
        threshold <- max(quantile(dists, 0.95, na.rm = TRUE), 
                        median(dists, na.rm = TRUE) + 2 * mad(dists, na.rm = TRUE))
        
        outliers_in_group <- samples_in_sample_class[dists > threshold]
        
        if (length(outliers_in_group) > 0) {
            outlier_samples <- c(outlier_samples, outliers_in_group)
            for (outlier in outliers_in_group) {
                outlier_info[[outlier]] <- data.frame(
                    sample = outlier,
                    sample_class = ct,
                    distance_to_center = dists[outlier],
                    threshold = threshold,
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    
    if (length(outlier_samples) > 0) {
        message("Potential outlier samples detected based on distance to sample_class center:")
        print(outlier_samples)
        
        # Find potential switch partners for each outlier
        switch_partners <- list()
        for (outlier in outlier_samples) {
            outlier_sample_class <- meta[outlier, "sample_class"]
            outlier_coords <- um[outlier, ]
            
            # Find closest sample_class center (excluding own sample_class)
            other_sample_classs <- setdiff(rownames(sample_class_centers), outlier_sample_class)
            if (length(other_sample_classs) > 0) {
                distances_to_centers <- apply(sample_class_centers[other_sample_classs, , drop = FALSE], 1, function(center) {
                    sqrt(sum((outlier_coords - center)^2))
                })
                closest_sample_class <- names(which.min(distances_to_centers))
                
                switch_partners[[outlier]] <- data.frame(
                    outlier = outlier,
                    outlier_sample_class = outlier_sample_class,
                    distance_to_own_center = sample_distances[outlier],
                    potential_switch_sample_class = closest_sample_class,
                    distance_to_switch_center = min(distances_to_centers),
                    stringsAsFactors = FALSE
                )
            }
        }
        
        if (length(switch_partners) > 0) {
            switch_df <- do.call(rbind, switch_partners)
            message("\nPotential sample switch partners:")
            print(switch_df)
            save_table("tables/umap", "umap_outlier_samples_with_switches.csv", switch_df)
        }
        
        # Save detailed outlier info
        outlier_detail_df <- do.call(rbind, outlier_info)
        save_table("tables/umap", "umap_outlier_samples_details.csv", outlier_detail_df)
        
        # Label outliers in UMAP plot
        um_df$outlier <- ifelse(rownames(um_df) %in% outlier_samples, "Outlier", "Normal")
        p <- ggplot(um_df, aes(UMAP1, UMAP2, color = outlier)) +
            geom_point(alpha = 0.8, size = 6) +
            scale_color_manual(values = c("Normal" = "grey", "Outlier" = "red"), na.translate = FALSE) +
            theme_pca_min() + labs(title = "UMAP with Outliers Labeled", color = "Sample Status")
        save_plot("plots/umap", "umap_with_outliers.svg", p)
        
        # Additional plot: show distances to centers
        dist_df <- data.frame(
            sample = names(sample_distances),
            sample_class = meta[names(sample_distances), "sample_class"],
            distance = sample_distances,
            is_outlier = names(sample_distances) %in% outlier_samples,
            stringsAsFactors = FALSE
        )
        p_dist <- ggplot(dist_df, aes(x = sample_class, y = distance, color = is_outlier)) +
            geom_jitter(width = 0.2, alpha = 0.7, size = 3) +
            scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), labels = c("Normal", "Outlier")) +
            theme_pca_min() +
            labs(title = "Distance to sample_class Center (Median)", x = "sample_class", y = "UMAP Distance", color = "Status") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        save_plot("plots/umap", "umap_distance_to_centers.svg", p_dist)
        
    } else {
        message("No outlier samples detected based on distance to sample_class centers.")
    }
}

# 7) Clustering on PCs and ARI/NMI vs metadata
pc_for_cluster <- scale(pca$x[, 1:npc, drop = FALSE])
best_k <- 2:8

if (!requireNamespace("cluster", quietly = TRUE)) {
  message("cluster not installed; silhouette may be skipped.")
}
if (!requireNamespace("aricode", quietly = TRUE)) {
  message("aricode not installed; ARI/NMI will be skipped.")
}

sil_tab <- list()
clu_eval <- list()
for (k in best_k) {
  km <- kmeans(pc_for_cluster, centers = k, nstart = 50, iter.max = 100)
  sil_avg <- NA_real_
  if (requireNamespace("cluster", quietly = TRUE)) {
    s <- cluster::silhouette(km$cluster, dist(pc_for_cluster))
    sil_avg <- mean(s[, "sil_width"], na.rm = TRUE)
  }
  sil_tab[[length(sil_tab)+1]] <- data.frame(k=k, silhouette=sil_avg)

  if (requireNamespace("aricode", quietly = TRUE)) {
    for (lab in c("sample_class", "condition", "AnimalID")) {
      if (!lab %in% names(meta)) next
      ref <- factor(meta[[lab]])
      pred <- factor(km$cluster, levels = sort(unique(km$cluster)))
      ari <- aricode::ARI(ref, pred)
      nmi <- aricode::NMI(ref, pred)
      clu_eval[[length(clu_eval)+1]] <- data.frame(k=k, label=lab, ARI=ari, NMI=nmi, stringsAsFactors = FALSE)
    }
  }
}
sil_df <- if (length(sil_tab)) do.call(rbind, sil_tab) else data.frame()
if (nrow(sil_df)) save_table("tables/clustering", "cluster_silhouette.csv", sil_df)
clu_df <- if (length(clu_eval)) do.call(rbind, clu_eval) else data.frame()
if (nrow(clu_df)) save_table("tables/clustering", "cluster_ARI_NMI.csv", clu_df)

suggest_k <- if (nrow(sil_df) && is.finite(max(sil_df$silhouette, na.rm=TRUE))) sil_df$k[which.max(sil_df$silhouette)] else 3
km_final <- kmeans(pc_for_cluster, centers = suggest_k, nstart = 100, iter.max = 200)
meta$kmeans_cluster <- factor(km_final$cluster)
pc_grp_plot <- plot_and_save_group("kmeans_cluster", "PCA", "pca_by_kmeans_cluster.svg")
if (inherits(pc_grp_plot, "gg")) save_plot("plots/clustering", "pca_by_kmeans_cluster.svg", pc_grp_plot)

# 8) Batch correction check: pre/post plate (or ReplicateGroup) using limma
run_remove_batch <- function(mat_samples_by_feature, meta_df, batch_key = "plate", design_keys = c("sample_class")) {
    if (!requireNamespace("limma", quietly = TRUE)) {
        message("limma not installed; skipping batch correction.")
        return(NULL)
    }
    if (!batch_key %in% names(meta_df)) {
        message("Batch key not in meta: ", batch_key); return(NULL)
    }
    x <- mat_samples_by_feature
    smp <- colnames(x)
    meta_local <- meta_df[smp, , drop = FALSE]
    if (!identical(rownames(meta_local), smp)) stop("Meta/sample alignment failed for batch removal.")
    batch <- factor(meta_local[[batch_key]])
    if (length(batch) != ncol(x)) stop("Batch length != number of samples in matrix.")
    dk <- intersect(design_keys, names(meta_local))
    if (!length(dk)) {
        design <- matrix(1, ncol(x), 1); colnames(design) <- "Intercept"
    } else {
        design <- stats::model.matrix(as.formula(paste0("~ 0 + ", paste(dk, collapse = "+"))), data = meta_local)
    }
    if (nrow(design) != ncol(x)) stop("Row dimension of design must equal number of samples (columns of x).")
    limma::removeBatchEffect(x, batch = batch, design = design)
}

mat_adj <- try(run_remove_batch(mat_samples_by_feature = mat, meta_df = meta, batch_key = "plate", design_keys = c("sample_class")), silent = TRUE)
if (!inherits(mat_adj, "try-error") && !is.null(mat_adj)) {
    pca_prec <- pca
    pca_post <- prcomp(t(mat_adj), center = TRUE, scale. = TRUE)
    pre_plot <- plot_and_save_group("plate", "PCA (pre-correction)", "pca_by_plate_precorrect.svg")
    if (inherits(pre_plot, "gg")) save_plot("plots/batch", "pca_by_plate_precorrect.svg", pre_plot)
    pca <- pca_post
    post_plot <- plot_and_save_group("plate", "PCA (post-correction)", "pca_by_plate_postcorrect.svg")
    if (inherits(post_plot, "gg")) save_plot("plots/batch", "pca_by_plate_postcorrect.svg", post_plot)
    pca <- pca_prec
} else {
    message("Batch correction skipped or failed; check plate variable presence and design alignment.")
}

