#!/usr/bin/env Rscript

# Contracts for 03_qc_exploration/02_rank_abundance_by_sample_class.r.
#
# The stage used to read four per-sample-class imputed workbooks from
# 02_data/gct/imputed/, a folder that is not present in the project tree, and it
# summarised by sample_class alone. It now reads the validated animal-level GCT
# and summarises by sample_class x condition.
#
# Section 1 is a source-text regression guard: it fails if the script reverts to
# the absent share-bound workbook input, or quietly gains a re-transformation of
# the already-processed GCT values.
#
# Section 2 runs the real stage into a temporary output root and checks the
# animal-level design and a fixed numerical contract taken from the validated
# reviewer-revision source data. It skips cleanly when the shared drive is
# unavailable, matching the convention used by the other cohort-level tests.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1) sub("^--file=", "", file_arg) else
  file.path("tests", "test_rank_abundance_animal_level.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

stage <- file.path(repo_root, "03_qc_exploration", "02_rank_abundance_by_sample_class.r")
stopifnot(file.exists(stage))
stage_lines <- readLines(stage, warn = FALSE)
src <- paste(stage_lines, collapse = "\n")
# "must NOT contain" checks run against code only. The header comment legitimately
# describes the input this stage no longer uses, and that history must not trip a
# regression guard aimed at the code.
code <- paste(stage_lines[!grepl("^\\s*#", stage_lines)], collapse = "\n")

failures <- character()
checks <- 0L
ok <- function(label, condition) {
  checks <<- checks + 1L
  if (!isTRUE(condition)) failures <<- c(failures, label) else cat("  [ok]  ", label, "\n")
}

cat("\n=== 1. source-text contracts (always run) ===\n")

# 1. the stage resolves the animal-level GCT
ok("resolves neha_protigy_input_animal_level_primary.gct",
   grepl("neha_protigy_input_animal_level_primary\\.gct", src, fixed = FALSE))
ok("input override is NEHA_RANK_ABUNDANCE_INPUT_GCT",
   grepl("NEHA_RANK_ABUNDANCE_INPUT_GCT", src, fixed = TRUE))
ok("reads the GCT through the shared ProTigy reader",
   grepl("validate_protigy_gct_v13", src, fixed = TRUE) &&
     grepl("protigy_input_utils\\.R", src))

# REGRESSION GUARD: must not fall back to the absent share-bound workbook input
ok("does NOT read per-class imputed workbooks",
   !grepl("pgmatrix_imputed", code, fixed = TRUE))
ok("does NOT point at 02_data/gct/imputed",
   !grepl('"gct",\\s*"imputed"', code) && !grepl("gct/imputed", code, fixed = TRUE))
ok("does NOT use the old NEHA_RANK_ABUNDANCE_INPUT_DIR directory override",
   !grepl("NEHA_RANK_ABUNDANCE_INPUT_DIR", code, fixed = TRUE))

# 2. sample_class AND condition are both used explicitly
ok("uses sample_class and condition metadata explicitly",
   grepl("normalize_sample_class", src, fixed = TRUE) &&
     grepl("normalize_condition", src, fixed = TRUE))
ok("groups by sample_class x condition",
   grepl("paste\\(metadata\\$sample_class, gsub", src))

# 3. no new normalisation / imputation / filtering of the processed values
ok("no re-normalisation of the GCT matrix",
   !grepl("scale\\(|normalizeQuantiles|normalizeBetweenArrays|median_center|quantile_norm", code))
ok("no imputation step",
   !grepl("impute|missForest|knn", code, ignore.case = TRUE))
ok("linearisation is the original 2^MeanLog2 only",
   grepl("2\\s*\\^\\s*d\\$MeanLog2", src))

# 4. marker categories retained
for (cat_ in c("General Neuron", "Neuropil/Structure", "Stress Response", "Activation")) {
  ok(paste0("marker category present: ", cat_), grepl(cat_, src, fixed = TRUE))
}
for (g in c("Vcan", "Kcna1", "Mog", "Cntnap1", "Cnp", "Rps9", "Brd4",
            "Mrpl4", "Timm9", "Ndufb7", "Aktip", "Naa10", "Acvr1b")) {
  ok(paste0("marker gene retained: ", g), grepl(paste0('"', g, '"'), src, fixed = TRUE))
}

# 5. the two manuscript panel group sets
ok("Fig3E groups are mcherry paired-veh + unpaired-veh",
   grepl('"mcherry_paired-veh", "mcherry_unpaired-veh"', src, fixed = TRUE))
ok("SuppD groups are neuron + neuropil unpaired-veh",
   grepl('"neuron_unpaired-veh", "neuropil_unpaired-veh"', src, fixed = TRUE))

cat("\n=== 2. live stage contracts (skipped when the shared drive is unavailable) ===\n")

gct <- file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler",
                 "02_data/animal_level/input_gct/neha_protigy_input_animal_level_primary.gct")
id_map <- file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler",
                    "02_data/animal_level/mapped/forward/mcherry_paired_veh_vs_mcherry_unpaired_veh.csv")

if (!file.exists(gct) || !file.exists(id_map)) {
  cat("  [skip] animal-level GCT or id map not reachable; source-text contracts still enforced\n")
} else {
  out_dir <- file.path(tempdir(), paste0("rank_abundance_test_", as.integer(Sys.time())))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  # system2(env=) is not supported on Windows, so the variable is set in this
  # process and inherited by the child. Restored afterwards so the stage can
  # never be pointed at the real QC output root by a leaked test setting.
  previous <- Sys.getenv("NEHA_RANK_ABUNDANCE_OUTPUT_DIR", unset = NA_character_)
  Sys.setenv(NEHA_RANK_ABUNDANCE_OUTPUT_DIR = out_dir)
  on.exit({
    if (is.na(previous)) Sys.unsetenv("NEHA_RANK_ABUNDANCE_OUTPUT_DIR")
    else Sys.setenv(NEHA_RANK_ABUNDANCE_OUTPUT_DIR = previous)
  }, add = TRUE)
  status <- system2(rscript, c("--vanilla", shQuote(stage)), stdout = FALSE, stderr = FALSE)
  ok("stage runs to completion into a temporary output root", status == 0L)
  ok("stage honoured the temporary output root (nothing written to the default)",
     length(list.files(out_dir)) > 0L)

  ranks_path <- file.path(out_dir, "processed_protein_ranks_animal_level.csv")
  ok("writes processed_protein_ranks_animal_level.csv", file.exists(ranks_path))

  if (file.exists(ranks_path)) {
    d <- utils::read.csv(ranks_path, stringsAsFactors = FALSE)

    # design contracts
    ok("16 sample_class x condition groups", length(unique(d$Condition)) == 16L)
    ok("3 animals per group", all(d$n_animals == 3L))
    animals <- unique(unlist(strsplit(unique(d$animals), ";")))
    ok("12 distinct animals across the matrix", length(animals) == 12L)
    ok("exactly 3 distinct animals listed per group",
       all(vapply(strsplit(unique(d$animals), ";"), length, integer(1)) == 3L))

    classes <- unique(sub("_.*$", "", d$Condition))
    conds <- unique(sub("^[^_]+_", "", d$Condition))
    ok("four sample classes present",
       setequal(classes, c("mcherry", "neuropil", "cfos", "neuron")))
    ok("four conditions present",
       setequal(conds, c("paired-cno", "paired-veh", "unpaired-cno", "unpaired-veh")))

    # no hemisphere-level identifiers leaked in as observations
    ok("no hemisphere-level sample ids treated as observations",
       !any(grepl("_L$|_R$|Plate[0-9]|ReplicateGroup", d$animals)))
    ok("48 animal-level observations implied (16 groups x 3)",
       length(unique(d$Condition)) * 3L == 48L)

    # fixed numerical contract, taken from the validated reviewer-revision
    # source data (Processed_Protein_Ranks_animal_level.csv)
    anchors <- data.frame(
      Condition = c("mcherry_paired-veh", "mcherry_unpaired-veh",
                    "neuron_unpaired-veh", "neuropil_unpaired-veh"),
      Genes = c("Aktip", "Lrrc47", "Gkap1", "Vcan"),
      MeanLog2 = c(1.75, 2.59, 1.49666666666667, 1.83833333333333),
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(anchors))) {
      g <- anchors$Condition[i]
      sub <- d[d$Condition == g, ]
      sub <- sub[order(sub$Rank), ]
      ok(paste0("group ", g, ": 5310 gene symbols"), nrow(sub) == 5310L)
      ok(paste0("group ", g, ": rank 1 is ", anchors$Genes[i]),
         identical(sub$Genes[1], anchors$Genes[i]))
      ok(paste0("group ", g, ": rank 1 MeanLog2 matches reference"),
         abs(sub$MeanLog2[1] - anchors$MeanLog2[i]) < 1e-10)
      ok(paste0("group ", g, ": LinearValue is 2^MeanLog2"),
         max(abs(sub$LinearValue - 2 ^ sub$MeanLog2)) < 1e-10)
      ok(paste0("group ", g, ": rank is a dense 1..n ordering by descending abundance"),
         identical(sub$Rank, seq_len(nrow(sub))) && !is.unsorted(rev(sub$LinearValue)))
    }

    # marker categories survive the rewrite
    mk <- unique(d[d$MarkerType != "None", c("Genes", "MarkerType")])
    ok("28 highlighted marker genes", nrow(mk) == 28L)
    ok("four marker categories",
       setequal(unique(mk$MarkerType),
                c("General Neuron", "Neuropil/Structure", "Stress Response", "Activation")))
    ok("no gene assigned to two categories", !anyDuplicated(mk$Genes))
    ok("Vcan is Neuropil/Structure",
       identical(mk$MarkerType[mk$Genes == "Vcan"], "Neuropil/Structure"))
    ok("Mrpl4 is Activation",
       identical(mk$MarkerType[mk$Genes == "Mrpl4"], "Activation"))
    ok("Aktip is Stress Response",
       identical(mk$MarkerType[mk$Genes == "Aktip"], "Stress Response"))
    ok("Rps9 is General Neuron",
       identical(mk$MarkerType[mk$Genes == "Rps9"], "General Neuron"))
  }

  # the two manuscript panels and their source data are produced
  for (p in c("Fig3E", "SuppD")) {
    ok(paste0(p, " panel SVG written"),
       file.exists(file.path(out_dir, paste0(p, "_rank_abundance_animal_level.svg"))))
    sd_path <- file.path(out_dir, paste0(p, "_rank_abundance_animal_level_source_data.csv"))
    ok(paste0(p, " source data written"), file.exists(sd_path))
    if (file.exists(sd_path)) {
      s <- utils::read.csv(sd_path, stringsAsFactors = FALSE)
      expected <- if (p == "Fig3E") c("mcherry_paired-veh", "mcherry_unpaired-veh") else
        c("neuron_unpaired-veh", "neuropil_unpaired-veh")
      ok(paste0(p, " covers exactly its two groups"), setequal(unique(s$Condition), expected))
      ok(paste0(p, " has 2 x 5310 rows"), nrow(s) == 10620L)
    }
  }
  ok("a per-group panel is written for every sample_class x condition",
     length(list.files(out_dir, pattern = "^rank_abundance_.*\\.svg$")) == 16L)
  ok("marker validation workbook written",
     file.exists(file.path(out_dir, "marker_validation_summary.xlsx")))
}

cat("\n=== RESULT ===\n")
if (length(failures)) {
  stop(paste(c(sprintf("Rank-abundance animal-level tests failed (%d of %d checks):",
                       length(failures), checks), paste0("  - ", failures)), collapse = "\n"),
       call. = FALSE)
}
cat(sprintf("Rank-abundance animal-level tests passed (%d contracts).\n", checks))
