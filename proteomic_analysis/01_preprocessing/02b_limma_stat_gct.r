# ================================================================
# Animal-level limma statistics, written as a ProTigy-format statistical GCT
#
# Reads the animal-level input GCT (48 columns: 12 animals x 4 sample classes) and writes
# a statistical GCT with the same layout as the ProTigy export, so that
# 03_gct_extractR.r, MapThatProt and the enrichment stages run on it unchanged.
#
# PROTEOMICS_LIMMA_VARIANCE_MODEL:
#   within_compartment (default)  one limma model per sample class: 12 animals,
#       4 conditions, n = 3 per group. The variance of each compartment is estimated
#       from that compartment only. This is the "within sample class" model recorded in
#       input_gct/protigy_contrast_manifest.csv.
#   pooled  one limma model over all 48 columns (16 groups), so the 4 compartments
#       share one variance estimate per protein. Reproduces the 2026-08-24 ProTigy
#       export (stat_results_for_ssGSEA_neha_proteome.gct); kept to check this script.
#
# Both models: lmFit on ~0 + group, contrasts.fit, eBayes(robust = TRUE, trend = FALSE), and
# Benjamini-Hochberg within each comparison ("significant" = adj.P.Val < 0.05).
# logFC is the same under both models. AveExpr is the mean over the columns in the model:
# the compartment's 12 columns (within_compartment) or all 48 (pooled, close to 0).
# ================================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("01_preprocessing", "02b_limma_stat_gct.r")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "project_path_utils.R"))

required_packages <- c("limma", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages before running this stage: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

option_or_env <- function(option_name, env_name, default) {
  value <- getOption(option_name)
  if (!is.null(value) && nzchar(trimws(as.character(value)))) return(as.character(value))
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(trimws(value))) return(value)
  default
}

animal_root <- file.path(PROTEOMICS_DATA_ROOT, "02_data", "animal_level")
protigy_export <- file.path(animal_root, "stat_results_for_ssGSEA_neha_proteome.gct")

variance_model <- option_or_env("proteomics.limma_variance_model", "PROTEOMICS_LIMMA_VARIANCE_MODEL", "within_compartment")
if (!variance_model %in% c("within_compartment", "pooled")) {
  stop("PROTEOMICS_LIMMA_VARIANCE_MODEL must be within_compartment or pooled.", call. = FALSE)
}
input_gct <- project_normalize_path(option_or_env(
  "proteomics.limma_input_gct",
  "PROTEOMICS_LIMMA_INPUT_GCT",
  file.path(animal_root, "input_gct", "neha_protigy_input_animal_level_primary.gct")
), must_work = TRUE)
output_gct <- project_normalize_path(option_or_env(
  "proteomics.limma_output_gct",
  "PROTEOMICS_LIMMA_OUTPUT_GCT",
  file.path(animal_root, "within_compartment_20261124", "stat_results_within_compartment.gct")
))
force <- tolower(option_or_env("proteomics.limma_force", "PROTEOMICS_LIMMA_FORCE", "false")) %in% c("true", "1", "yes")

protected <- c(protigy_export, file.path(PROTEOMICS_DATA_ROOT, "99_historical"), input_gct)
if (any(vapply(protected, function(p) project_paths_overlap(output_gct, p), logical(1)))) {
  stop("Refusing to write over the ProTigy export, the input GCT or 99_historical: ", output_gct, call. = FALSE)
}
if (file.exists(output_gct) && !force) {
  stop("Output exists; set PROTEOMICS_LIMMA_FORCE=true to replace it: ", output_gct, call. = FALSE)
}

# ---- input --------------------------------------------------------------------------
lines <- readLines(input_gct, warn = FALSE)
dims <- scan(text = lines[[2]], what = integer(), quiet = TRUE)
stopifnot(trimws(lines[[1]]) == "#1.3", length(dims) == 4L, dims[[3]] == 1L)
n_proteins <- dims[[1]]
n_meta <- dims[[4]]
split_line <- function(x) strsplit(x, "\t", fixed = TRUE)[[1]]
header <- split_line(lines[[3]])
meta_rows <- lapply(lines[3L + seq_len(n_meta)], split_line)
names(meta_rows) <- vapply(meta_rows, `[[`, character(1), 1L)
body <- do.call(rbind, lapply(lines[3L + n_meta + seq_len(n_proteins)], split_line))
stopifnot(ncol(body) == length(header), nrow(body) == n_proteins, !anyDuplicated(body[, 1]))

sample_cols <- header[-(1:2)]
column_index <- seq_along(sample_cols) + 2L
annotation <- data.frame(
  column = sample_cols,
  sample_class = meta_rows$sample_class[column_index],
  condition = meta_rows$condition[column_index],
  stringsAsFactors = FALSE
)
stopifnot(
  all(annotation$sample_class %in% sample_classes),
  all(annotation$condition %in% condition_levels),
  all(table(annotation$sample_class, annotation$condition) == 3L)
)
values <- matrix(as.numeric(body[, column_index]), nrow = n_proteins, dimnames = list(body[, 1], sample_cols))
if (anyNA(values)) stop("The input matrix contains missing values; this stage expects the imputed matrix.", call. = FALSE)

# ---- models --------------------------------------------------------------------------
manifest <- primary_contrast_manifest()

fit_group_model <- function(columns, groups) {
  groups <- factor(groups)
  design <- stats::model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  list(fit = limma::lmFit(values[, columns, drop = FALSE], design), design = design)
}

contrast_table <- function(model, case, reference) {
  contrast <- limma::makeContrasts(contrasts = paste0(case, " - ", reference), levels = model$design)
  fit <- limma::eBayes(limma::contrasts.fit(model$fit, contrast), robust = TRUE)
  tab <- limma::topTable(fit, number = Inf, sort.by = "none", adjust.method = "BH")
  stopifnot(identical(rownames(tab), rownames(values)))
  list(table = tab, df_prior = stats::median(fit$df.prior), s2_prior = fit$s2.prior, df_residual = fit$df.residual[[1]])
}

results <- list()
model_rows <- list()
if (identical(variance_model, "pooled")) {
  model <- fit_group_model(annotation$column, paste(annotation$sample_class, annotation$condition, sep = "_"))
  for (i in seq_len(nrow(manifest))) {
    r <- contrast_table(model, manifest$case_phenotype[[i]], manifest$reference_phenotype[[i]])
    results[[manifest$canonical_comparison[[i]]]] <- r$table
    model_rows[[i]] <- data.frame(canonical_comparison = manifest$canonical_comparison[[i]], n_columns_in_model = nrow(annotation),
                                  df_residual = r$df_residual, df_prior_median = r$df_prior, s2_prior = r$s2_prior)
  }
} else {
  for (cls in unique(manifest$sample_class)) {
    a <- annotation[annotation$sample_class == cls, , drop = FALSE]
    model <- fit_group_model(a$column, a$condition)
    for (i in which(manifest$sample_class == cls)) {
      r <- contrast_table(model, manifest$case_condition[[i]], manifest$reference_condition[[i]])
      results[[manifest$canonical_comparison[[i]]]] <- r$table
      model_rows[[i]] <- data.frame(canonical_comparison = manifest$canonical_comparison[[i]], n_columns_in_model = nrow(a),
                                    df_residual = r$df_residual, df_prior_median = r$df_prior, s2_prior = r$s2_prior)
    }
  }
}

# ---- output in the ProTigy export layout ---------------------------------------------
fmt_num <- function(x) trimws(formatC(x, digits = 15, format = "g"))
metric_fields <- list()
signed_fields <- list()
for (comparison in manifest$canonical_comparison) {
  tab <- results[[comparison]]
  log_p <- -log10(tab$P.Value)
  block <- data.frame(
    logFC = fmt_num(tab$logFC),
    AveExpr = fmt_num(tab$AveExpr),
    t = fmt_num(tab$t),
    P.Value = fmt_num(tab$P.Value),
    adj.P.Val = fmt_num(tab$adj.P.Val),
    B = fmt_num(tab$B),
    significant = ifelse(tab$adj.P.Val < 0.05, "TRUE", "FALSE"),
    Log.P.Value = fmt_num(log_p),
    stringsAsFactors = FALSE
  )
  names(block) <- paste0(names(block), ".", comparison)
  metric_fields[[comparison]] <- block
  signed_fields[[comparison]] <- stats::setNames(
    data.frame(fmt_num(round(sign(tab$logFC) * log_p, 4)), stringsAsFactors = FALSE),
    paste0("sign.logP.", comparison)
  )
}
out <- cbind(
  data.frame(id = body[, 1], Description = body[, 2], stringsAsFactors = FALSE),
  do.call(cbind, unname(metric_fields)),
  stats::setNames(as.data.frame(body[, column_index], stringsAsFactors = FALSE), sample_cols),
  do.call(cbind, unname(signed_fields))
)
n_data <- nrow(manifest)
n_descriptors <- ncol(out) - 1L - n_data

dir.create(dirname(output_gct), recursive = TRUE, showWarnings = FALSE)
con <- file(output_gct, open = "wb")
writeLines(c("#1.3", paste(n_proteins, n_data, n_descriptors, 0L, sep = "\t"), paste(names(out), collapse = "\t")), con, sep = "\n")
utils::write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, eol = "\n")
close(con)

run_manifest <- do.call(rbind, model_rows)
run_manifest$variance_model <- variance_model
run_manifest$n_fdr_proteins <- vapply(results[run_manifest$canonical_comparison], function(t) sum(t$adj.P.Val < 0.05), integer(1))
run_manifest$input_gct <- input_gct
run_manifest$input_gct_sha256 <- project_file_sha256(input_gct)
run_manifest$output_gct <- output_gct
run_manifest$output_gct_sha256 <- project_file_sha256(output_gct)
run_manifest$limma_version <- as.character(utils::packageVersion("limma"))
run_manifest$r_version <- R.version.string
run_manifest$run_timestamp_utc <- format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")
utils::write.csv(run_manifest, sub("\\.gct$", "_run_manifest.csv", output_gct), row.names = FALSE)

cat("Variance model:", variance_model, "\n")
cat("Wrote:", output_gct, "\n")
print(run_manifest[, c("canonical_comparison", "n_columns_in_model", "df_residual", "df_prior_median", "n_fdr_proteins")], row.names = FALSE)
