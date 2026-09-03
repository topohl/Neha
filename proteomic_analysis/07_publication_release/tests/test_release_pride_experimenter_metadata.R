#!/usr/bin/env Rscript

# Contract: experimenter-supplied metadata is source-controlled, scoped to this mouse
# proteomics dataset, and propagated to PRIDE without per-sample invention.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else
  file.path("07_publication_release", "tests", "test_release_pride_experimenter_metadata.R")
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/",
                           mustWork = FALSE)
if (!file.exists(file.path(repo_root, "R", "analysis_labels.R"))) {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "project_path_utils.R"))
layer <- file.path(repo_root, "07_publication_release")
source(file.path(layer, "R", "release_utils.R"))
source(file.path(layer, "R", "release_validation.R"))

failures <- 0L
expect <- function(condition, message) {
  if (!isTRUE(condition)) {
    failures <<- failures + 1L
    cat("  [FAIL]", message, "\n")
  } else {
    cat("  [ok]  ", message, "\n")
  }
}

cat("=== checked-in metadata sources ===\n")
metadata <- release_read_experimenter_metadata(repo_root)
protocol <- release_read_sample_preparation_protocol(repo_root)
expect(nrow(metadata) == length(RELEASE_EXPERIMENTER_METADATA_REQUIRED_IDS),
       "every required experimenter metadata record is present exactly once")
expect(nrow(protocol) == 7L && identical(as.integer(protocol$step_order), 1:7),
       "the complete seven-step protocol is ordered")
expect(identical(unique(protocol$protocol_version_date), "2024-11-28"),
       "the protocol version date is 2024-11-28")
expect(setequal(metadata$sample_class[metadata$category == "sample_class_definition"],
                c("cfos", "mcherry", "neuron", "neuropil")),
       "all four marker-defined LCM classes are defined")
expect(release_metadata_value(metadata, "animalid_sex_assignments") == "not available",
       "AnimalID-level sex remains unavailable")
expect(release_metadata_value(metadata, "age_at_experiment_start") == "7-10 weeks" &&
         release_metadata_value(metadata, "age_at_tissue_collection") == "not available",
       "age at experiment start is distinct from age at tissue collection")
expect(release_metadata_value(metadata, "search_modification_parameters") == "not available",
       "CAA/TCEP use does not fill search modification settings")

out_root <- release_output_root()
sdrf_path <- file.path(out_root, "pride", "sdrf.tsv")
if (!file.exists(sdrf_path)) {
  cat("No built publication release at", out_root, "\n")
  cat("Skipping generated PRIDE contracts.\n")
} else {
  rd <- function(p) read.delim(p, sep = "\t", quote = "", stringsAsFactors = FALSE,
                               check.names = FALSE, na.strings = character(0),
                               colClasses = "character")
  sdrf <- rd(sdrf_path)
  released_metadata <- rd(file.path(out_root, "metadata", "experimenter_metadata.tsv"))
  released_protocol <- rd(file.path(out_root, "metadata", "sample_preparation_protocol.tsv"))
  field_status <- rd(file.path(out_root, "pride", "sdrf_field_status.tsv"))

  cat("\n=== generated SDRF and PRIDE artefacts ===\n")
  expect(identical(released_metadata, metadata),
         "released experimenter metadata is byte-semantically identical to source")
  expect(identical(released_protocol, protocol),
         "released preparation protocol is byte-semantically identical to source")
  expect(all(sdrf[["characteristics[organism part]"]] ==
               release_metadata_value(metadata, "organism_part", TRUE)),
         "CeM provenance reaches the SDRF through the supported anatomy representation")
  expect(all(sdrf[["characteristics[cell type]"]] == "not applicable"),
         "marker-defined LCM material does not claim unsupported cell types")
  expect(all(sdrf[["factor value[sample class]"]] %in%
               c("cfos", "mcherry", "neuron", "neuropil")) &&
           all(grepl("laser-capture microdissection",
                     sdrf[["characteristics[sampling site]"]], fixed = TRUE)),
         "sample-class and sampling-site fields preserve the LCM categories")
  expect(nrow(sdrf) == 96L &&
           all(sdrf[["characteristics[organism]"]] == "Mus musculus") &&
           all(sdrf[["characteristics[strain or breed]"]] == "C57BL/6J") &&
           all(sdrf[["characteristics[sex]"]] == "not available") &&
           all(sdrf[["characteristics[age]"]] == "not available"),
         "96 mouse samples carry strain without invented per-sample sex or age at collection")
  expect(all(sdrf[["comment[label]"]] == "label free sample"),
         "label-free proteomics is encoded correctly")
  templates <- sdrf[, which(names(sdrf) == "comment[sdrf template]"), drop = FALSE]
  expect(ncol(templates) == 2L &&
           all(templates[[1]] == "NT=ms-proteomics;VV=v1.1.0") &&
           all(templates[[2]] == "NT=vertebrates;VV=v1.1.0"),
         "the applicable SDRF templates are represented as two repeatable fields")
  cleavage <- sdrf[, which(names(sdrf) == "comment[cleavage agent details]"), drop = FALSE]
  expect(ncol(cleavage) == 2L &&
           all(cleavage[[1]] == "NT=Lys-C;AC=MS:1001309") &&
           all(cleavage[[2]] == "NT=Trypsin;AC=MS:1001251"),
         "Lys-C and trypsin are encoded in two ordered cleavage-agent columns")
  expect(all(sdrf[["comment[modification parameters]"]] == "not available") &&
           all(sdrf[["comment[instrument]"]] == "not available") &&
           all(sdrf[["comment[proteomics data acquisition method]"]] == "not available"),
         "search modifications, instrument and acquisition mode remain unresolved")
  expect(setequal(field_status$sdrf_field[
                    field_status$status == "MISSING_REQUIRED_METADATA" &
                      field_status$sdrf_field %in% names(sdrf)],
                  c("characteristics[developmental stage]", "comment[instrument]",
                    "comment[proteomics data acquisition method]")),
         "the remaining required-field gap is exact")

  protocol_text <- paste(readLines(file.path(out_root, "pride",
                                             "SAMPLE_PREPARATION_PROTOCOL.md"),
                                   warn = FALSE), collapse = " ")
  expected_terms <- c("5% DDM", "5 mM TCEP", "20 mM CAA", "0.1 M TEAB",
                      "95 degrees C", "1 hour", "4 ng Lys-C per sample", "2 hours",
                      "6 ng trypsin per sample", "overnight", "TFA to 1%", "SpeedVac",
                      "Evotip cleanup", "4.2 uL A buffer")
  expect(all(vapply(expected_terms, function(x) grepl(x, protocol_text, fixed = TRUE),
                    logical(1))),
         "the human-readable PRIDE protocol preserves every supplied preparation detail")
}

cat("\n=== RESULT ===\n")
if (failures > 0L) {
  stop(sprintf("PRIDE experimenter-metadata contracts failed: %d", failures),
       call. = FALSE)
}
cat("All PRIDE experimenter-metadata contracts hold.\n")
