# ====================================================================
# Animal-level PCA workflow -- ENTRY POINT
#
# Was a single 2,576-line script; split 2026-08-26 into ordered parts under
# 03_qc_exploration/pca/. This file is a thin orchestrator, so the documented
# invocation is unchanged:
#
#     Rscript 03_qc_exploration/06_pcaPlot_Neha.r
#
# The parts run at top level and share globals (mat, meta, pca, output_dir and the
# plotting/saving helpers) created by 06a_pca_core.r, exactly as when they were one file.
# Part order matters: 06b creates rot/cors_df/um/npc that 06c, 06g and 06h consume.
#
# Behaviour change (deliberate): previously every extension shared ONE tryCatch, so the
# first failure aborted all remaining extensions. Each part is now guarded separately, so
# one failing visualization no longer suppresses the rest. Audit semantics are preserved --
# a failure still writes the failure audit to pca_audit_path (last writer wins) and the
# script still exits non-zero, but only after attempting every part.
# ====================================================================

.neha_script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.neha_script_path <- if (length(.neha_script_arg)) {
  sub("^--file=", "", .neha_script_arg[[1]])
} else {
  file.path("03_qc_exploration", "06_pcaPlot_Neha.r")
}
.neha_pca_dir <- normalizePath(file.path(dirname(.neha_script_path), "pca"), winslash = "/", mustWork = FALSE)
if (!dir.exists(.neha_pca_dir)) {
  .neha_pca_dir <- normalizePath(file.path(getwd(), "03_qc_exploration", "pca"), winslash = "/", mustWork = TRUE)
}

# Core must succeed: it validates the input contract and writes its own audit on failure.
source(file.path(.neha_pca_dir, "06a_pca_core.r"))

.neha_extension_parts <- c(
  "06b_pca_loadings_umap_clustering.r",
  "06c_pca_variants_and_sensitivity.r",
  "06d_pca_supplementary_heatmaps_and_norm.r",
  "06e_pca_distributions_and_tests.r",
  "06f_pca_qc_and_main_figures.r",
  "06g_pca_stats_summaries.r",
  "06h_pca_innovative_visualizations.r"
)

.neha_failed <- character()
for (.neha_part in .neha_extension_parts) {
  .neha_err <- tryCatch({
    source(file.path(.neha_pca_dir, .neha_part))
    NULL
  }, error = identity)
  if (!is.null(.neha_err)) {
    .neha_failed <- c(.neha_failed, sprintf("%s: %s", .neha_part, conditionMessage(.neha_err)))
    message("PCA extension part failed, continuing with the rest: ", .neha_part)
    message("  ", conditionMessage(.neha_err))
    write_dt(
      make_neha_pca_audit(
        validated_input = core_input$validated,
        prepared_pca = core_input$prepared,
        source_path = gct_file,
        source_sha256 = source_sha256,
        output_paths = primary_output_paths,
        execution_status = "failed",
        error_message = paste(.neha_failed, collapse = " | ")
      ),
      pca_audit_path
    )
  }
}

if (length(.neha_failed)) {
  stop(
    "One or more PCA extension parts failed:\n  ",
    paste(.neha_failed, collapse = "\n  "),
    call. = FALSE
  )
}
message("PCA workflow completed: core + ", length(.neha_extension_parts), " extension parts.")
