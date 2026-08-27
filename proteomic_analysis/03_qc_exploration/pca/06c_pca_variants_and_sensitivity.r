# ====================================================================
# Sparse biplot, sessionInfo, ellipses, leave-one-plate-out, IRLBA, ROBPCA, t-SNE, gap statistic, PC-protein volcanoes
# Part of the PCA workflow split out of the former monolithic
# 03_qc_exploration/06_pcaPlot_Neha.r (2026-08-26). Sourced in order by that
# script, which remains the entry point. Runs at top level and shares the
# globals created by 06a_pca_core.r (mat, meta, pca, output_dir, helpers).
# Consumes cors_df from 06b.
# ====================================================================

if (!exists("mat") || !exists("meta") || !exists("pca") || !exists("output_dir")) {
  stop("PCA core state missing. Run 03_qc_exploration/06_pcaPlot_Neha.r, or source pca/06a_pca_core.r first.", call. = FALSE)
}

# 9) Optional biplot with sparse loadings on PC1/PC2
plot_biplot_sparse <- function(topN = 20, out_file = "pca_biplot_sparse.svg"){
  R <- as.data.frame(pca$rotation)
  R$protein <- rownames(R)
  for (k in c("PC1","PC2")) if (!k %in% names(R)) return(invisible(NULL))
  sel <- unique(c(head(R[order(-R$PC1), "protein"], topN),
                  head(R[order(R$PC1),  "protein"], topN),
                  head(R[order(-R$PC2), "protein"], topN),
                  head(R[order(R$PC2),  "protein"], topN)))
  dd <- data.frame(pca$x[,1:2, drop=FALSE], meta[rownames(pca$x), , drop=FALSE])
  pal <- make_modern_palette(nlevels(factor(dd$sample_class)))
  p <- ggplot(dd, aes(PC1, PC2, color=factor(sample_class))) +
    geom_point(size=3, alpha=0.8) +
    scale_color_manual(values = pal) +
    theme_pca_min() + labs(title="PCA biplot (sparse loadings)", color="sample_class")
  arrows <- R[R$protein %in% sel, c("PC1","PC2","protein")]
  rngx <- diff(range(dd$PC1)); rngy <- diff(range(dd$PC2))
  scale_fac <- 0.5 * min(rngx, rngy)
  arrows$PC1s <- arrows$PC1 * scale_fac
  arrows$PC2s <- arrows$PC2 * scale_fac
  p <- p +
    geom_segment(data=arrows, aes(x=0, y=0, xend=PC1s, yend=PC2s), inherit.aes = FALSE, arrow = arrow(length = unit(0.15,"cm")), color="#444444", alpha=0.7) +
    ggrepel::geom_text_repel(data=arrows, aes(x=PC1s, y=PC2s, label=protein), inherit.aes = FALSE, size=3, max.overlaps = 200)
  save_plot("plots/biplot", out_file, p, width=8, height=6.5)
}
plot_biplot_sparse(20, "pca_biplot_sparse.svg")

# 10) Save session info
ensure_dir(subdir("tables/meta"))
writeLines(c(capture.output(sessionInfo())), con = file.path(subdir("tables/meta"), "sessionInfo.txt"))

# ================== Additional Extensions (organized outputs) ==================


# A) Confidence ellipses and centroids on PCA
add_pca_ellipses <- function(key, fname){
  grp <- build_group(key)
  keep <- !is.na(grp)
  grp2 <- droplevels(grp[keep])
  if (nlevels(grp2) < 2) return(invisible(NULL))
  df <- data.frame(pca$x[keep, 1:2, drop=FALSE], group = grp2, check.names = FALSE)
  varp <- (pca$sdev^2)/sum(pca$sdev^2)
  p <- ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(alpha=0.9, size=3) +
    stat_ellipse(type="t", level=0.95, linewidth=0.6) +
    stat_summary(aes(group=group, color=group), fun = mean, geom = "point", size=4, shape=4, stroke=1.2) +
    scale_color_manual(values = make_modern_palette(nlevels(grp2))) +
    theme_pca_min() +
    labs(title = paste("PCA with ellipses by", key),
         x = sprintf("PC1 (%.1f%%)", 100*varp[1]),
         y = sprintf("PC2 (%.1f%%)", 100*varp[2]),
         color = key)
  save_plot("plots/ellipses", fname, p)
}
add_pca_ellipses("sample_class", "pca_ellipses_sample_class.svg")

# B) Leave-one-batch/plate-out sensitivity (if plate exists)
leave_one_batch_pca <- function(batch_key = "plate"){
  if (!batch_key %in% names(meta)) { message("No batch key ", batch_key); return(invisible(NULL)) }
  res <- list()
  full_rot <- pca$rotation
  if (is.null(full_rot) || nrow(full_rot) == 0) { message("No rotation in base PCA; skipping similarity metrics."); full_rot <- NULL }
  batches <- na.omit(unique(meta[[batch_key]]))
  for (b in batches) {
    keep <- meta[[batch_key]] != b & !is.na(meta[[batch_key]])
    if (sum(keep) < 3) next
    Xsamp <- t(mat[, keep, drop=FALSE])
    pp <- try(prcomp(Xsamp, center=TRUE, scale.=TRUE), silent=TRUE)
    if (inherits(pp, "try-error")) next

    cs <- NA_real_
    if (!is.null(full_rot) && !is.null(pp$rotation)) {
      common <- intersect(rownames(full_rot), rownames(pp$rotation))
      if (length(common) >= 10) {
        F <- full_rot[common, 1:min(5, ncol(full_rot)), drop=FALSE]
        R <- pp$rotation[common, 1:ncol(F), drop=FALSE]
        cs <- sapply(seq_len(ncol(F)), function(i) {
          fi <- F[,i]; ri <- R[,i]
          den <- sqrt(sum(fi^2) * sum(ri^2))
          if (!is.finite(den) || den == 0) return(NA_real_)
          abs(sum(fi*ri) / den)
        })
      } else {
        want <- if (!is.null(full_rot)) min(5, ncol(full_rot)) else 5
        cs <- rep(NA_real_, want)
      }
    }
    var_ex <- (pp$sdev^2)/sum(pp$sdev^2)
    cs_vec <- if (length(cs) < length(var_ex)) c(cs, rep(NA_real_, length(var_ex)-length(cs))) else cs[seq_along(var_ex)]
    res[[length(res)+1]] <- data.frame(batch_left_out = b,
                                       PC = paste0("PC", seq_along(var_ex)),
                                       variance = var_ex,
                                       cos_sim = cs_vec)

    subsamp <- rownames(pp$x)
    subsamp <- intersect(subsamp, rownames(meta))
    df_sc <- data.frame(pp$x[subsamp,1:2, drop=FALSE],
                        sample_class = droplevels(factor(meta[subsamp, "sample_class"])),
                        check.names = FALSE)
    pal <- make_modern_palette(nlevels(droplevels(df_sc$sample_class)))
    p <- ggplot(df_sc, aes(PC1, PC2, color=sample_class)) +
         geom_point(size=3, alpha=0.8) +
         scale_color_manual(values=pal, drop=FALSE, na.translate=FALSE) +
         theme_pca_min() + labs(title=sprintf("PCA (-%s)", b))
    save_plot("plots/leave_one_batch", sprintf("pca_sample_class_minus_%s.svg", b), p)
  }
  if (length(res)) {
    tab <- do.call(rbind, res)
    save_table("leave_one_batch", "pca_leave_one_batch_summary.csv", tab)
  }
}
leave_one_batch_pca("plate")

# C) Robust PCA variants
# Assumes run_irlba_pca() is already defined earlier and returns a prcomp-like list
run_irlba_pca <- function(){
  if (!requireNamespace("irlba", quietly = TRUE)) { message("irlba not installed; skipping IRLBA PCA."); return(NULL) }
  X <- scale(t(mat), center=TRUE, scale=TRUE)
  rownames_X <- rownames(X)
  k <- min(10, ncol(X)-1, nrow(X)-1)
  if (k < 2) { message("Not enough samples/features for IRLBA PCA."); return(NULL) }
  irlba_res <- irlba::irlba(X, nv = k, nu = k)
  pcx <- irlba_res$u %*% diag(irlba_res$d)
  colnames(pcx) <- paste0("PC", seq_len(ncol(pcx)))
  rownames(pcx) <- rownames_X
  list(x = pcx, sdev = irlba_res$d / sqrt(max(1, nrow(X)-1)), rotation = NULL)
}


# C1) IRLBA PCA (fast truncated SVD on scaled data), aligned rownames
pca_irlba <- run_irlba_pca()
if (!is.null(pca_irlba)) {
  if (is.null(rownames(pca_irlba$x))) rownames(pca_irlba$x) <- rownames(t(mat))
  samp <- rownames(pca_irlba$x)
  samp <- intersect(samp, rownames(meta))
  if (length(samp) >= 2) {
    df <- data.frame(pca_irlba$x[samp, 1:2, drop=FALSE],
                     sample_class = droplevels(factor(meta[samp, "sample_class"])),
                     check.names = FALSE)
    pal <- make_modern_palette(nlevels(droplevels(df$sample_class)))
    p <- ggplot(df, aes(PC1, PC2, color = sample_class)) +
      geom_point(size=3, alpha=0.8) +
      scale_color_manual(values=pal, drop=FALSE, na.translate=FALSE) +
      theme_pca_min() + labs(title="IRLBA PCA by sample_class", color="sample_class")
    save_plot("plots/irlba", "irlba_pca_by_sample_class.svg", p)
  } else {
    message("IRLBA: no overlapping samples to plot with metadata.")
  }
}

# C2) ROBPCA (outlier-robust) if available
run_robpca <- function(){
  if (!requireNamespace("rospca", quietly = TRUE)) { message("rospca not installed; skipping ROBPCA."); return(NULL) }
  X <- scale(t(mat), center=TRUE, scale=TRUE)
  rownames_X <- rownames(X)
  rp <- rospca::robpca(X, k = min(10, ncol(X)))
  pcx <- rp$scores
  colnames(pcx) <- paste0("PC", seq_len(ncol(pcx)))
  rownames(pcx) <- rownames_X
  list(x = pcx, sdev = sqrt(rp$eigenvalues), rotation = NULL)
}
pca_ro <- run_robpca()
if (!is.null(pca_ro)) {
  samp <- rownames(pca_ro$x)
  samp <- intersect(samp, rownames(meta))
  df <- data.frame(pca_ro$x[samp,1:2, drop=FALSE],
                   sample_class = droplevels(factor(meta[samp, "sample_class"])),
                   check.names = FALSE)
  pal <- make_modern_palette(nlevels(droplevels(df$sample_class)))
  p <- ggplot(df, aes(PC1, PC2, color = sample_class)) +
       geom_point(size=3, alpha=0.8) +
       scale_color_manual(values=pal, drop=FALSE, na.translate=FALSE) +
       theme_pca_min() + labs(title="ROBPCA by sample_class", color="sample_class")
  save_plot("plots/robpca", "robpca_by_sample_class.svg", p)
}

# D) t-SNE with parameter sweeps (on top PCs), aligned to meta and safe perplexity
run_tsne_and_plot <- function(){
  if (!requireNamespace("Rtsne", quietly = TRUE)) { message("Rtsne not installed; skipping t-SNE."); return(NULL) }
  samp <- rownames(pca$x)
  X <- pca$x[samp, 1:min(20, ncol(pca$x)), drop=FALSE]
  N <- nrow(X)
  if (N < 5) { message("Too few samples for t-SNE."); return(NULL) }
  for (perp in c(10, 30, 50)) {
    p_ok <- max(1, min(perp, floor((N - 1)/3)))
    set.seed(42 + perp)
    ts <- Rtsne::Rtsne(X, perplexity = p_ok, pca = FALSE, theta = 0.5, check_duplicates = FALSE, verbose = FALSE)
    emb <- ts$Y
    rownames(emb) <- samp
    df <- data.frame(tSNE1 = emb[,1], tSNE2 = emb[,2], meta[samp, , drop = FALSE], check.names = FALSE)
    plot_tsne <- function(key, out){
      if (!key %in% names(df)) return(NULL)
      grp <- factor(trim_ws(as.character(df[[key]])))
      pal <- make_modern_palette(nlevels(droplevels(grp)))
      p <- ggplot(df, aes(tSNE1, tSNE2, color = grp)) +
        geom_point(alpha=0.9, size=3) +
        scale_color_manual(values = pal, na.translate = FALSE) +
        theme_pca_min() + labs(title = paste("t-SNE by", key, "(perp", p_ok, ")"), color = key)
      save_plot("plots/tsne", out, p)
    }
    plot_tsne("sample_class", sprintf("tsne_sample_class_perp%d.svg", p_ok))
  }
}
run_tsne_and_plot()

# E) Gap statistic for clustering (PC space)
compute_gap <- function(){
  if (!requireNamespace("cluster", quietly = TRUE)) { message("cluster not installed; skipping gap."); return(NULL) }
  X <- scale(pca$x[, 1:min(10, ncol(pca$x)), drop=FALSE])
  set.seed(42)
  gap <- cluster::clusGap(X, FUN = stats::kmeans, K.max = 10, B = 50, nstart = 25, iter.max = 100)
  gap_df <- as.data.frame(gap$Tab)
  gap_df$k <- seq_len(nrow(gap_df))
  save_table("clustering", "cluster_gap_statistic.csv", gap_df)
  p <- ggplot(gap_df, aes(k, gap)) + geom_line() + geom_point() + theme_pca_min() + labs(title="Gap statistic", x="k", y="Gap")
  save_plot("clustering", "cluster_gap_statistic.svg", p)
}
compute_gap()

# F) PC–protein correlation volcano plots (uses existing correlations if available)
if (exists("cors_df") && nrow(cors_df)) {
  for (pc in unique(cors_df$pc)) {
    df <- cors_df[cors_df$pc == pc, , drop=FALSE]
    df$ml10q <- -log10(pmax(df$q, .Machine$double.xmin))
    p <- ggplot(df, aes(r, ml10q)) +
      geom_point(alpha=0.4, size=1.3, color="#6B5B95") +
      geom_vline(xintercept = c(-0.3, 0.3), linetype="dashed", color="#B2B2B2") +
      geom_hline(yintercept = -log10(0.05), linetype="dashed", color="#B2B2B2") +
      theme_pca_min() + labs(title = paste("Protein–", pc, "correlations"), x="Pearson r", y="-log10(q)")
    save_plot("correlations", sprintf("protein_pc_%s_volcano.svg", pc), p)
  }
}


