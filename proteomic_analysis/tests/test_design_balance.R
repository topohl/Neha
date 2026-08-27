#!/usr/bin/env Rscript

# Design-balance / identifiability checks for the animal-level Neha design.
#
# Motivation: the plate-vs-Pairing confound in this cohort was only found by manual
# provenance archaeology across two separate audits. These checks make the aliasing an
# explicit, machine-checked property of the design rather than something to rediscover.
#
# They are descriptive contracts, not pass/fail judgements about the biology: the confound
# is a real property of the collected cohort and cannot be fixed in software. The tests
# fail only if the design's *identifiability structure* differs from what has been audited
# and documented in CANONICAL_OUTPUTS.md -- i.e. if someone changes the assignment metadata
# without updating the documented caveats.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else file.path("tests", "test_design_balance.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))

failures <- 0L
expect <- function(condition, message) {
  if (!isTRUE(condition)) {
    failures <<- failures + 1L
    cat("  [FAIL]", message, "\n")
  } else {
    cat("  [ok]  ", message, "\n")
  }
}

option_or_env <- function(option_name, env_name, default) {
  value <- getOption(option_name)
  if (!is.null(value) && nzchar(trimws(as.character(value)))) return(as.character(value))
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(trimws(value))) return(value)
  default
}

assignment_path <- option_or_env(
  "neha.source_sample_assignment",
  "NEHA_SOURCE_SAMPLE_ASSIGNMENT",
  "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/02_data/animal_level/input_gct/source_sample_assignment.csv"
)

# ---------------------------------------------------------------------------
# Pure helpers (testable without the shared drive)
# ---------------------------------------------------------------------------

#' Split a canonical condition label into its 2x2 factors.
derive_design_factors <- function(condition) {
  condition <- as.character(condition)
  data.frame(
    condition = condition,
    Pairing = ifelse(grepl("^unpaired", condition), "unpaired", "paired"),
    Treatment = ifelse(grepl("cno$", condition), "cno", "veh"),
    stringsAsFactors = FALSE
  )
}

#' How strongly is `covariate` aliased with `factor_values`?
#'
#' Returns the fraction of observations that follow the single best one-to-one mapping
#' between covariate levels and factor levels. 1.0 means perfect aliasing (the covariate
#' is a relabelling of the factor and the two cannot be separated in any model).
aliasing_fraction <- function(covariate, factor_values) {
  covariate <- as.character(covariate)
  factor_values <- as.character(factor_values)
  keep <- !is.na(covariate) & !is.na(factor_values)
  covariate <- covariate[keep]; factor_values <- factor_values[keep]
  if (!length(covariate)) return(NA_real_)
  tab <- table(covariate, factor_values)
  # best achievable agreement: for each covariate level, credit its most common factor level
  sum(apply(tab, 1L, max)) / sum(tab)
}

#' Is a within-stratum contrast free of covariate variation?
#'
#' TRUE when every observation in the stratum shares one covariate level, meaning an
#' additive covariate effect cannot contribute to a contrast computed inside it.
stratum_is_covariate_free <- function(covariate) {
  length(unique(as.character(covariate[!is.na(covariate)]))) <= 1L
}

cat("=== derive_design_factors ===\n")
d <- derive_design_factors(condition_levels)
expect(identical(d$Pairing, c("paired", "paired", "unpaired", "unpaired")),
       "Pairing derived correctly from the four canonical conditions")
expect(identical(d$Treatment, c("cno", "veh", "cno", "veh")),
       "Treatment derived correctly from the four canonical conditions")

cat("\n=== aliasing_fraction ===\n")
expect(isTRUE(all.equal(aliasing_fraction(c("A","A","B","B"), c("x","x","y","y")), 1)),
       "perfect aliasing detected as 1.0")
expect(isTRUE(all.equal(aliasing_fraction(c("A","B","A","B"), c("x","x","y","y")), 0.5)),
       "balanced covariate detected as 0.5")
expect(isTRUE(all.equal(aliasing_fraction(c("A","A","A","B"), c("x","x","y","y")), 0.75)),
       "partial aliasing quantified correctly")

cat("\n=== stratum_is_covariate_free ===\n")
expect(isTRUE(stratum_is_covariate_free(c("Plate2","Plate2","Plate2"))),
       "single-level stratum reported as covariate-free")
expect(isFALSE(stratum_is_covariate_free(c("Plate1","Plate2","Plate1"))),
       "mixed stratum reported as NOT covariate-free")

# ---------------------------------------------------------------------------
# Real-cohort contracts (skipped when the shared drive is unavailable)
# ---------------------------------------------------------------------------
if (!file.exists(assignment_path)) {
  cat("\n=== cohort checks SKIPPED: assignment metadata not reachable ===\n")
  cat("    expected at:", assignment_path, "\n")
  cat("    set NEHA_SOURCE_SAMPLE_ASSIGNMENT to run the cohort-level contracts.\n")
} else {
  src <- utils::read.csv(assignment_path, stringsAsFactors = FALSE)
  # Plate is embedded in the raw acquisition (.d) filename; it is acquisition-side provenance.
  src$Plate <- sub(".*_(Plate[0-9]+)_.*", "\\1", src$sample_id)
  fac <- derive_design_factors(src$condition)
  src$Pairing <- fac$Pairing
  src$Treatment <- fac$Treatment

  animals <- unique(src[, c("AnimalID", "condition", "Pairing", "Treatment", "Plate")])

  cat("\n=== cohort: Plate is an animal-level property ===\n")
  per_animal <- tapply(animals$Plate, animals$AnimalID, function(x) length(unique(x)))
  expect(all(per_animal == 1L),
         "each AnimalID maps to exactly one Plate (so Plate is animal-level, not per-tissue)")

  cat("\n=== cohort: documented aliasing structure ===\n")
  a_pairing <- aliasing_fraction(animals$Plate, animals$Pairing)
  a_treat <- aliasing_fraction(animals$Plate, animals$Treatment)
  cat("    Plate~Pairing aliasing  :", sprintf("%.3f", a_pairing), "\n")
  cat("    Plate~Treatment aliasing:", sprintf("%.3f", a_treat), "\n")

  expect(a_pairing >= 0.9,
         sprintf("Plate remains near-perfectly aliased with Pairing (%.3f >= 0.9) as documented", a_pairing))
  expect(a_treat < 0.9,
         sprintf("Plate is NOT strongly aliased with Treatment (%.3f < 0.9) as documented", a_treat))

  cat("\n=== cohort: which simple contrasts are plate-protected? ===\n")
  unpaired <- animals[animals$Pairing == "unpaired", ]
  paired <- animals[animals$Pairing == "paired", ]
  expect(isTRUE(stratum_is_covariate_free(unpaired$Plate)),
         "unpaired stratum is single-plate, so unpaired_cno-vs-unpaired_veh is free of an additive plate effect")
  expect(isFALSE(stratum_is_covariate_free(paired$Plate)),
         "paired stratum spans >1 plate, so paired_cno-vs-paired_veh is NOT fully plate-protected")

  cat("\n=== cohort: learning contrast is completely aliased (the headline caveat) ===\n")
  learn <- animals[animals$condition %in% c("paired_veh", "unpaired_veh"), ]
  a_learn <- aliasing_fraction(learn$Plate, learn$condition)
  cat("    Plate~condition aliasing within {paired_veh, unpaired_veh}:", sprintf("%.3f", a_learn), "\n")
  expect(isTRUE(all.equal(a_learn, 1)),
         "paired_veh-vs-unpaired_veh remains 100% aliased with Plate -- do not interpret it as biology alone")
}

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("Design-balance contracts failed: %d", failures), call. = FALSE)
}
cat("All design-balance contracts hold.\n")
