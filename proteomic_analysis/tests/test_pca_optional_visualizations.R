args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("tests", "test_pca_optional_visualizations.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "pca_animal_level_utils.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "pca_animal_level_utils.R"))

expect_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

empty_data <- prepare_neha_pca_variance_treemap_data(c(NA_real_, Inf, 0, -1))
empty_renderer_called <- FALSE
primary_audit <- data.frame(execution_status = "success", stringsAsFactors = FALSE)
empty_result <- run_neha_pca_optional_plot(
  plot_name = "variance_treemap",
  usable_data = empty_data,
  output_path = tempfile(fileext = ".svg"),
  n_input_rows = attr(empty_data, "n_input_rows"),
  render_fun = function(plot_data, output_path) {
    empty_renderer_called <<- TRUE
    stop("empty renderer must not run")
  }
)
expect_true(!empty_renderer_called, "Empty treemap data invoked the plotting path.")
expect_true(identical(empty_result$execution_status, "skipped_no_usable_rows"),
            "Empty treemap input was not skipped deterministically.")
expect_true(identical(empty_result$n_input_rows, 4L) && identical(empty_result$n_usable_rows, 0L),
            "Empty treemap skip audit does not record input and usable row counts.")
expect_true(identical(primary_audit$execution_status, "success"),
            "An optional empty treemap changed primary PCA success.")

nonempty_data <- prepare_neha_pca_variance_treemap_data(c(40, 30, 20, 10))
nonempty_renderer_called <- FALSE
nonempty_result <- run_neha_pca_optional_plot(
  plot_name = "variance_treemap",
  usable_data = nonempty_data,
  output_path = tempfile(fileext = ".svg"),
  n_input_rows = attr(nonempty_data, "n_input_rows"),
  render_fun = function(plot_data, output_path) {
    nonempty_renderer_called <<- TRUE
    expect_true(identical(plot_data$Variance, c(40, 30, 20, 10)),
                "Nonempty treemap data changed before plotting.")
  }
)
expect_true(nonempty_renderer_called, "Nonempty treemap input did not follow the existing plotting path.")
expect_true(identical(nonempty_result$execution_status, "success"),
            "Successful nonempty treemap rendering was not audited as success.")

if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("treemapify", quietly = TRUE)) {
  rendered_path <- tempfile(fileext = ".svg")
  rendered_result <- run_neha_pca_optional_plot(
    plot_name = "variance_treemap",
    usable_data = nonempty_data,
    output_path = rendered_path,
    render_fun = function(plot_data, output_path) {
      plot <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(area = Variance, fill = Category, label = PC)
      ) +
        treemapify::geom_treemap() +
        treemapify::geom_treemap_text(color = "white", place = "centre")
      ggplot2::ggsave(output_path, plot, width = 4, height = 3)
    }
  )
  expect_true(identical(rendered_result$execution_status, "success") && file.exists(rendered_path),
              "Nonempty data did not complete the existing treemapify rendering path.")
  unlink(rendered_path)
}

render_failure <- run_neha_pca_optional_plot(
  plot_name = "variance_treemap",
  usable_data = nonempty_data,
  output_path = tempfile(fileext = ".svg"),
  render_fun = function(plot_data, output_path) stop("fixture renderer failure")
)
expect_true(identical(render_failure$execution_status, "skipped_render_error") &&
              grepl("fixture renderer failure", render_failure$reason, fixed = TRUE),
            "A renderer-only failure was not contained and audited.")
expect_true(identical(primary_audit$execution_status, "success"),
            "An optional renderer failure changed primary PCA success.")

sample_ids <- paste0("sample", seq_len(6L))
pca_scores <- matrix(
  seq_len(30L), nrow = 6L,
  dimnames = list(sample_ids, paste0("PC", seq_len(5L)))
)
metadata <- data.frame(
  sample_class = rep(c("mcherry", "neuropil"), each = 3L),
  row.names = sample_ids,
  stringsAsFactors = FALSE
)
centroids <- prepare_neha_pca_group_centroids(pca_scores, metadata, "sample_class", n_pcs = 5L)
expect_true(nrow(centroids) == 2L && setequal(centroids$group, c("mcherry", "neuropil")),
            "Keyed PCA group centroids did not retain the nonempty optional chord path.")

missing_metadata <- metadata
missing_metadata$sample_class <- NA_character_
empty_centroids <- prepare_neha_pca_group_centroids(pca_scores, missing_metadata, "sample_class", n_pcs = 5L)
expect_true(nrow(empty_centroids) == 0L,
            "An optional group plot with no usable metadata did not produce an empty guarded input.")

cat("Optional PCA visualization regression tests passed.\n")
