#!/usr/bin/env Rscript

# Publication release, stage 09 -- PRIDE / SDRF-Proteomics metadata.
#
# Produces
#   pride/sdrf.tsv                     SDRF-Proteomics sample-and-data-relationship file
#   pride/sdrf_field_status.tsv        per-column READY / MISSING_REQUIRED / MISSING_OPTIONAL
#   pride/SDRF_MISSING_METADATA.md     what is missing, why, and what would supply it
#   pride/README_PRIDE.md              which files are suitable for PRIDE and which are not
#
# RULE: no acquisition metadata is invented. Where a value cannot be proven from a file in
# the project, the SDRF cell carries the SDRF-sanctioned string "not available" and the
# field is recorded as missing, with the specific document that would supply it.
#
# Fields deliberately NOT guessed: instrument model, digestion enzyme, DIA acquisition
# method, search-software version, labelling chemistry, database-search parameters,
# organism part, animal sex, animal age.
#
# One field needs its reasoning stated rather than hidden. comment[fraction identifier] is
# mandatory in SDRF and is set to 1. That encodes the observed one-run-per-biological-sample
# mapping (96 acquisition runs for 96 distinct AnimalID x sample_class x hemisphere
# samples), which is what the value 1 means in SDRF. It is NOT a claim about a fractionation
# protocol, and it is still listed for confirmation in SDRF_MISSING_METADATA.md.

suppressWarnings({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1L) dirname(sub("^--file=", "", file_arg)) else "07_publication_release"
})
source(file.path(here, "R", "release_utils.R"))
REPO_ROOT <- release_repo_root()
release_source_project_helpers(REPO_ROOT)
source(file.path(REPO_ROOT, "07_publication_release", "R", "release_validation.R"))
release_require("digest")

DATA_ROOT <- release_data_root()
PROJECT_ROOT <- release_project_root()
OUT_ROOT <- release_prepare_output_root()
STAGE <- "07_publication_release/09_build_pride_sdrf.R"

release_banner("stage 09 -- PRIDE / SDRF")

NOT_AVAILABLE <- "not available"
NOT_APPLICABLE <- "not applicable"

rel <- function(...) release_path(..., create_dir = FALSE)
sample_metadata <- release_read_tsv_plain(rel("metadata", "sample_metadata.tsv"))
field_provenance <- release_read_tsv_plain(rel("metadata", "metadata_field_provenance.tsv"))
sample_class_corrections <- release_read_tsv_plain(rel("metadata", "sample_class_corrections.tsv"))
n_class_corrected <- nrow(sample_class_corrections)

if (nrow(sample_metadata) != RELEASE_DESIGN_INVARIANTS$n_measurement_records) {
  stop("Sample metadata does not hold 96 measurement records.", call. = FALSE)
}

# --------------------------------------------------------------------------------------
# are the raw acquisition files actually in hand?
# --------------------------------------------------------------------------------------
# The names are proven; the files are a separate question. Point
# PROTEOMICS_RELEASE_RAW_FILE_ROOT at the store holding the `.d` directories to have this
# resolve to TRUE. Absent that, the release reports them as not yet located, which is what
# keeps the PRIDE verdict honest.

raw_root <- release_option_or_env("proteomics.release_raw_file_root",
                                  "PROTEOMICS_RELEASE_RAW_FILE_ROOT", "")
raw_present_count <- 0L
if (nzchar(raw_root) && dir.exists(raw_root)) {
  present <- vapply(sample_metadata$raw_file_basename, function(b) {
    p <- file.path(raw_root, b)
    dir.exists(p) || file.exists(p)
  }, logical(1))
  raw_present_count <- sum(present)
}
raw_files_identified <- all(nzchar(sample_metadata$raw_file_basename)) &&
  !anyNA(sample_metadata$raw_file_basename)
raw_files_present <- raw_present_count == nrow(sample_metadata)
release_log("  raw acquisition filenames identified: ", raw_files_identified,
            " (", sum(nzchar(sample_metadata$raw_file_basename)), "/96)")
release_log("  raw acquisition files located on disk: ", raw_files_present,
            " (", raw_present_count, "/96",
            ifelse(nzchar(raw_root), paste0(" under ", raw_root),
                   "; PROTEOMICS_RELEASE_RAW_FILE_ROOT not set"), ")")

# --------------------------------------------------------------------------------------
# SDRF table
# --------------------------------------------------------------------------------------

sm <- sample_metadata[order(as.integer(sample_metadata$injection_index)), , drop = FALSE]

# Biological replicate index within sample_class x condition (1..3), stable and derived.
stratum <- paste(sm$sample_class, sm$condition, sep = "\r")
animal_by_stratum <- lapply(split(sm$AnimalID, stratum), function(x) sort(unique(x)))
biological_replicate <- vapply(seq_len(nrow(sm)), function(i) {
  match(sm$AnimalID[i], animal_by_stratum[[stratum[i]]])
}, integer(1))
if (max(biological_replicate) != RELEASE_DESIGN_INVARIANTS$n_animals_per_stratum) {
  stop("Biological replicate index does not run 1..3 within every stratum.", call. = FALSE)
}

sdrf <- data.frame(
  `source name` = paste(sm$AnimalID, sm$sample_class, tolower(sm$hemisphere), sep = "_"),
  `characteristics[organism]` = "Mus musculus",
  `characteristics[strain]` = NOT_AVAILABLE,
  `characteristics[organism part]` = NOT_AVAILABLE,
  `characteristics[cell type]` = NOT_AVAILABLE,
  `characteristics[disease]` = NOT_APPLICABLE,
  `characteristics[sex]` = NOT_AVAILABLE,
  `characteristics[age]` = NOT_AVAILABLE,
  `characteristics[individual]` = sm$AnimalID,
  `characteristics[biological replicate]` = biological_replicate,
  `characteristics[sampling site]` = sm$sample_class,
  `characteristics[anatomical side]` = sm$hemisphere,
  `characteristics[collection plate]` = sm$collection_plate,
  `assay name` = paste0("run ", as.integer(sm$injection_index)),
  `comment[data file]` = sm$raw_file_basename,
  `comment[file uri]` = NOT_AVAILABLE,
  `comment[technical replicate]` = 1L,
  `comment[fraction identifier]` = 1L,
  `comment[label]` = NOT_AVAILABLE,
  `comment[instrument]` = NOT_AVAILABLE,
  `comment[cleavage agent details]` = NOT_AVAILABLE,
  `comment[modification parameters]` = NOT_AVAILABLE,
  `comment[precursor mass tolerance]` = NOT_AVAILABLE,
  `comment[fragment mass tolerance]` = NOT_AVAILABLE,
  `comment[proteomics data acquisition method]` = NOT_AVAILABLE,
  `comment[collision energy]` = NOT_AVAILABLE,
  `comment[dissociation method]` = NOT_AVAILABLE,
  `comment[MS2 analyzer type]` = NOT_AVAILABLE,
  `comment[acquisition date]` = sm$acquisition_date,
  `comment[original acquisition run name]` = sm$sample_id,
  `comment[instrument alias in run name]` = sm$instrument_alias_token,
  `comment[lc and method token in run name]` = sm$lc_and_method_token,
  `factor value[pairing]` = sm$pairing_status,
  `factor value[treatment]` = sm$treatment,
  `factor value[sample class]` = sm$sample_class,
  check.names = FALSE, stringsAsFactors = FALSE
)

if (nrow(sdrf) != RELEASE_DESIGN_INVARIANTS$n_measurement_records) {
  stop("SDRF does not have one row per acquisition.", call. = FALSE)
}
if (anyDuplicated(sdrf[["comment[data file]"]])) {
  stop("SDRF comment[data file] is not unique.", call. = FALSE)
}
if (anyDuplicated(sdrf[["source name"]])) {
  stop("SDRF source name is not unique.", call. = FALSE)
}
release_log("  SDRF: ", nrow(sdrf), " rows x ", ncol(sdrf), " columns")

# --------------------------------------------------------------------------------------
# per-field status
# --------------------------------------------------------------------------------------

fs <- function(field, requirement, status, value_or_reason, why, likely_source,
               required_before_submission, notes = NA_character_) {
  data.frame(sdrf_field = field, sdrf_requirement = requirement, status = status,
             value_or_reason = value_or_reason, why_missing = why,
             likely_source = likely_source,
             required_before_submission = required_before_submission, notes = notes,
             stringsAsFactors = FALSE, check.names = FALSE)
}

field_status <- rbind(
  fs("source name", "required", "READY", "AnimalID_sampleclass_hemisphere", NA_character_,
     "derived from validated sample metadata", "no",
     "96 unique biological samples"),
  fs("characteristics[organism]", "required", "READY", "Mus musculus", NA_character_,
     "MOUSE_10090_idmapping.dat; organism_id = 10090 in the mapping stage; org.Mm.eg.db",
     "no", "NCBI taxid 10090"),
  fs("characteristics[strain]", "optional", "MISSING_OPTIONAL_METADATA", NOT_AVAILABLE,
     "no strain field exists in any project metadata file",
     "animal-facility records / manuscript methods", "no", NA_character_),
  fs("characteristics[organism part]", "required", "MISSING_REQUIRED_METADATA", NOT_AVAILABLE,
     "no dissected-region field or statement exists anywhere in the analysis tree",
     "manuscript methods / experimenter", "yes",
     paste("The EWCE cell-type reference is l1_amygdala.loom. That records a choice of",
           "reference, NOT the dissected tissue, and is not sufficient evidence.",
           "Deliberately not guessed.")),
  fs("characteristics[cell type]", "required", "MISSING_REQUIRED_METADATA", NOT_AVAILABLE,
     "the four sample classes are dissected/labelled fractions, not ontology cell types",
     "manuscript methods; then map each fraction to a Cell Ontology term", "yes",
     paste("The fraction label is carried in characteristics[sampling site] and",
           "factor value[sample class] so the information is not lost. `neuropil` in",
           "particular is a tissue compartment and not a cell type.")),
  fs("characteristics[disease]", "required", "READY", NOT_APPLICABLE,
     NA_character_, "experimental design (behavioural/chemogenetic, no disease model)",
     "no", NA_character_),
  fs("characteristics[sex]", "optional", "MISSING_OPTIONAL_METADATA", NOT_AVAILABLE,
     "sample_info.xlsx and sample_annotation.xlsx have no sex column",
     "animal-facility records", "no", NA_character_),
  fs("characteristics[age]", "optional", "MISSING_OPTIONAL_METADATA", NOT_AVAILABLE,
     "sample_info.xlsx and sample_annotation.xlsx have no age column",
     "animal-facility records", "no", NA_character_),
  fs("characteristics[individual]", "optional", "READY", "AnimalID (12 animals)",
     NA_character_, "validated sample metadata", "no", NA_character_),
  fs("characteristics[biological replicate]", "required", "READY", "1..3 within stratum",
     NA_character_, "derived from the validated animal-level design", "no",
     "3 animals per sample_class x condition"),
  fs("characteristics[sampling site]", "optional", "READY",
     "mcherry / neuropil / cfos / neuron", NA_character_, "validated sample metadata",
     "no", "the dissected/labelled fraction"),
  fs("characteristics[anatomical side]", "optional", "READY", "Left / Right", NA_character_,
     "explicit _left/_right labels in sample_annotation.xlsx", "no",
     "non-standard characteristic; retained because hemisphere is a real design variable"),
  fs("characteristics[collection plate]", "optional", "READY", "Plate1 / Plate2",
     NA_character_, "token in the acquisition run name; agrees with sample_info plate",
     "no",
     paste("Collection plate only. NOT a proteomics preparation, digestion, LC-MS,",
           "acquisition or instrument batch.")),
  fs("assay name", "required", "READY", "run 1 .. run 96", NA_character_,
     "injection index parsed from the acquisition run name", "no", NA_character_),
  fs("comment[data file]", "required", "READY", "<run>.d", NA_character_,
     "sample_annotation.xlsx Name column; basename minus .d equals the sample id for all 96",
     "no",
     paste("Names are proven. The FILES are a separate matter -- see README_PRIDE.md.")),
  fs("comment[file uri]", "optional", "MISSING_OPTIONAL_METADATA", NOT_AVAILABLE,
     "no public URI exists until a PRIDE accession is issued",
     "PRIDE, after submission", "no", NA_character_),
  fs("comment[technical replicate]", "required", "READY", "1", NA_character_,
     "each biological sample maps to exactly one acquisition run", "no",
     "Left/Right are anatomical subsamples of one animal, not technical replicates"),
  fs("comment[fraction identifier]", "required", "READY", "1", NA_character_,
     "96 acquisition runs for 96 distinct biological samples: a strict 1:1 mapping", "no",
     paste("The value 1 encodes the observed one-run-per-sample mapping, which is what",
           "SDRF means by fraction 1. It is NOT a statement about a fractionation",
           "protocol, and no fractionation protocol is documented. Confirm with the",
           "acquisition facility before submission.")),
  fs("comment[label]", "required", "MISSING_REQUIRED_METADATA", NOT_AVAILABLE,
     "no labelling chemistry is recorded in any project file",
     "acquisition facility / sample-preparation protocol", "yes",
     paste("The 1:1 sample-to-run mapping is inconsistent with isobaric multiplexing, but",
           "that is an argument from file counts, not a record of the chemistry used.",
           "Deliberately not guessed. If label-free, the SDRF value is",
           "'label free sample'.")),
  fs("comment[instrument]", "required", "MISSING_REQUIRED_METADATA", NOT_AVAILABLE,
     "no instrument model string exists anywhere in the project tree",
     "acquisition facility; instrument method file; or the .d acquisition metadata", "yes",
     paste("The run name contains the local instrument alias 'Olive', which is not a model.",
           "The acquisition format is `.d`, a vendor-specific acquisition directory format",
           "(Bruker). Neither identifies a model. Must be an ontology term from PSI-MS.")),
  fs("comment[cleavage agent details]", "required", "MISSING_REQUIRED_METADATA", NOT_AVAILABLE,
     "no digestion protocol is recorded in any project file",
     "sample-preparation protocol / manuscript methods", "yes",
     paste("A TRYP_PIG (porcine trypsin) entry appears among the identified protein groups,",
           "but standard contaminant databases include trypsin regardless of the enzyme",
           "actually used, so this is NOT evidence. Deliberately not guessed.")),
  fs("comment[modification parameters]", "required", "MISSING_REQUIRED_METADATA", NOT_AVAILABLE,
     "no search-parameter record (run log, config or report) exists in the project tree",
     paste("acquisition facility: the original search/quantification software run log or",
           "configuration (e.g. DIA-NN report.log.txt, if DIA-NN was used)"), "yes",
     NA_character_),
  fs("comment[precursor mass tolerance]", "optional", "MISSING_OPTIONAL_METADATA", NOT_AVAILABLE,
     "no search-parameter record exists in the project tree",
     "acquisition facility: original search/quantification software run log or configuration",
     "no", NA_character_),
  fs("comment[fragment mass tolerance]", "optional", "MISSING_OPTIONAL_METADATA", NOT_AVAILABLE,
     "no search-parameter record exists in the project tree",
     "acquisition facility: original search/quantification software run log or configuration",
     "no", NA_character_),
  fs("comment[proteomics data acquisition method]", "required", "MISSING_REQUIRED_METADATA",
     NOT_AVAILABLE, "no acquisition-method record exists in the project tree",
     "acquisition facility; instrument method file", "yes",
     paste("The retained processed files use the historical pg.matrix naming convention,",
           "but the exact upstream search/quantification software and configuration could",
           "not be recovered from the retained project files, so nothing retained",
           "establishes the acquisition mode used. Deliberately not guessed.")),
  fs("comment[collision energy]", "optional", "MISSING_OPTIONAL_METADATA", NOT_AVAILABLE,
     "no instrument method record exists in the project tree",
     "acquisition facility; instrument method file", "no", NA_character_),
  fs("comment[dissociation method]", "optional", "MISSING_OPTIONAL_METADATA", NOT_AVAILABLE,
     "no instrument method record exists in the project tree",
     "acquisition facility; instrument method file", "no", NA_character_),
  fs("comment[MS2 analyzer type]", "optional", "MISSING_OPTIONAL_METADATA", NOT_AVAILABLE,
     "no instrument method record exists in the project tree",
     "acquisition facility; instrument method file", "no", NA_character_),
  fs("comment[acquisition date]", "optional", "READY", "2024-12-17", NA_character_,
     "date token present in every acquisition run name; single acquisition date",
     "no",
     paste("Derived from the run-name token rather than from an instrument record;",
           "confirm against the acquisition log.")),
  fs("comment[original acquisition run name]", "optional", "READY",
     "verbatim run name", NA_character_, "validated sample metadata", "no",
     "non-standard comment; preserves the original identifier"),
  fs("comment[instrument alias in run name]", "optional", "READY", "Olive", NA_character_,
     "token in every acquisition run name", "no",
     "non-standard comment; a local alias, explicitly NOT an instrument model"),
  fs("comment[lc and method token in run name]", "optional", "READY",
     "FCo_Evo2_80SPDzoom_5cmRapid_Tobias", NA_character_,
     "token in every acquisition run name", "no",
     paste("non-standard comment; preserved verbatim and deliberately NOT decoded into an",
           "LC system, gradient or throughput claim")),
  fs("factor value[pairing]", "required", "READY", "paired / unpaired", NA_character_,
     "validated condition coding", "no", NA_character_),
  fs("factor value[treatment]", "required", "READY", "cno / veh", NA_character_,
     "validated condition coding", "no",
     "literal tokens from the project's condition coding, not expanded to chemical names"),
  fs("factor value[sample class]", "required", "READY",
     "mcherry / neuropil / cfos / neuron", NA_character_, "validated sample metadata",
     "no",
     paste0("The ANALYSIS-TIME sample class, i.e. the class the validated results",
            " were computed on. For ", n_class_corrected, " of 96 acquisitions this",
            " differs from the pre-correction records because the assignment was",
            " corrected during historical sample-identity QC; that correction is",
            " RESOLVED, not missing metadata. The pre-correction class is kept out of",
            " this SDRF factor deliberately -- overloading a factor value with two",
            " identities would misrepresent the design -- and is carried in the",
            " companion tables metadata/sample_metadata.tsv",
            " (original_sample_class) and metadata/sample_class_corrections.tsv.")),
  fs("raw acquisition FILES (not an SDRF column)", "required for deposition",
     ifelse(raw_files_present, "READY", "MISSING_REQUIRED_METADATA"),
     ifelse(raw_files_present, paste0(raw_present_count, "/96 located"),
            paste0(raw_present_count, "/96 located")),
     ifelse(raw_files_present, NA_character_,
            "the 96 .d acquisition directories are not present in the project tree"),
     "acquisition facility instrument store / institutional archive",
     ifelse(raw_files_present, "no", "yes"),
     paste("Set PROTEOMICS_RELEASE_RAW_FILE_ROOT to the directory holding the .d",
           "directories to re-evaluate. Names are proven; only the bytes are absent."))
)
stopifnot(all(field_status$status %in% c(release_sdrf_field_statuses)))

sdrf_columns <- names(sdrf)
covered <- intersect(sdrf_columns, field_status$sdrf_field)
if (length(covered) != length(sdrf_columns)) {
  stop("SDRF column(s) without a recorded status: ",
       paste(setdiff(sdrf_columns, field_status$sdrf_field), collapse = ", "), call. = FALSE)
}

n_missing_required <- sum(field_status$status == "MISSING_REQUIRED_METADATA" &
                            field_status$sdrf_field %in% sdrf_columns)
n_missing_optional <- sum(field_status$status == "MISSING_OPTIONAL_METADATA")
n_ready <- sum(field_status$status == "READY")

pride_status <- release_pride_status(raw_files_identified, raw_files_present,
                                     n_missing_required)
release_log("  SDRF fields: ", n_ready, " READY, ", n_missing_required,
            " MISSING_REQUIRED, ", n_missing_optional, " MISSING_OPTIONAL")
release_log("  PRIDE status: ", pride_status)

# --------------------------------------------------------------------------------------
# documents
# --------------------------------------------------------------------------------------

missing_tbl <- field_status[field_status$status != "READY", , drop = FALSE]
md_row <- function(x) paste0("| ", paste(x, collapse = " | "), " |")

missing_md <- c(
  "# SDRF metadata that could not be established from the project",
  "",
  paste0("Generated by `", STAGE, "`. Status vocabulary: `READY`, ",
         "`MISSING_REQUIRED_METADATA`, `MISSING_OPTIONAL_METADATA`."),
  "",
  paste0("**PRIDE readiness: `", pride_status, "`**"),
  "",
  paste0("- SDRF columns populated from verified metadata: **", n_ready, "**"),
  paste0("- SDRF-required fields that cannot be populated: **", n_missing_required, "**"),
  paste0("- Optional fields that cannot be populated: **", n_missing_optional, "**"),
  paste0("- Raw acquisition filenames identified: **", raw_files_identified,
         "** (96 of 96, proven from `sample_annotation.xlsx`)"),
  paste0("- Raw acquisition files located on disk: **", raw_files_present, "** (",
         raw_present_count, " of 96)"),
  "",
  "None of the values below were inferred. Where a plausible inference existed it is",
  "written out and explicitly rejected, so a later reader can see what was considered.",
  "",
  "| field | status | why_missing | likely_source | required_before_submission | notes |",
  "|---|---|---|---|---|---|",
  vapply(seq_len(nrow(missing_tbl)), function(i) {
    r <- missing_tbl[i, , drop = FALSE]
    md_row(c(paste0("`", r$sdrf_field, "`"), r$status,
             ifelse(is.na(r$why_missing), "", r$why_missing),
             ifelse(is.na(r$likely_source), "", r$likely_source),
             r$required_before_submission,
             ifelse(is.na(r$notes), "", gsub("\\|", "/", r$notes))))
  }, character(1)),
  "",
  "## What one document would fix most of this",
  "",
  "The original search/quantification software run log or configuration (e.g. DIA-NN",
  "`report.log.txt`, if DIA-NN was used) plus the instrument method for the 2024-12-17",
  "acquisition would supply, in one go: the search-software identity and version, the",
  "digestion enzyme, the modification parameters, the mass tolerances, the acquisition",
  "method, and the instrument model. Those six fields are the bulk of what is blocking a",
  "complete SDRF.",
  "",
  "## Explicitly rejected inferences",
  "",
  "| candidate inference | why it was rejected |",
  "|---|---|",
  md_row(c("instrument model from the run-name token `Olive`",
           "a local instrument alias, not a model; no mapping to a PSI-MS term exists in the project")),
  md_row(c("instrument vendor/model from the `.d` acquisition format",
           "the format identifies a vendor family, not a model; SDRF requires a specific instrument term")),
  md_row(c("digestion enzyme = trypsin from the `TRYP_PIG` protein group",
           "standard contaminant databases contain porcine trypsin regardless of the enzyme used")),
  md_row(c("acquisition method = DIA from the retained `pg.matrix` naming convention",
           paste("the historical naming convention does not identify the upstream",
                 "search/quantification software, let alone the acquisition mode used on",
                 "the instrument"))),
  md_row(c("LC system / gradient / throughput from `FCo_Evo2_80SPDzoom_5cmRapid`",
           "an operator method label; decoding it would be reading an acquisition protocol out of a filename")),
  md_row(c("organism part = amygdala from the `l1_amygdala.loom` EWCE reference",
           "the single-cell reference is an analysis choice; it does not record what tissue was dissected")),
  md_row(c("labelling = label free from the 1:1 sample-to-run mapping",
           "argues against isobaric multiplexing but is not a record of the chemistry used")),
  "",
  "## Sample class: RESOLVED, and deliberately not listed above",
  "",
  paste0("The sample-class assignment is **not** missing or unresolved metadata, and it is",
         " not one of the gaps this document lists. For ", n_class_corrected, " of the 96",
         " acquisitions -- the left hemispheres of ",
         paste(sort(unique(sample_class_corrections$AnimalID)), collapse = " and "),
         " -- the assignment was corrected during historical sample-identity quality",
         " control, and that correction is established."),
  "",
  paste0("`factor value[sample class]` and `characteristics[sampling site]` therefore carry",
         " the **analysis-time** class: the class the validated results were computed on.",
         " That is the correct primary factor representation. The pre-correction class is",
         " deliberately NOT packed into the SDRF factor -- overloading a factor value with",
         " two identities would misrepresent the experimental design, and SDRF has no",
         " sanctioned column for a superseded label. It is carried in the companion tables",
         " instead:"),
  "",
  "| where | field |",
  "|---|---|",
  md_row(c("`metadata/sample_metadata.tsv`",
           "`original_sample_class`, `analysis_sample_class`, `sample_class_corrected`")),
  md_row(c("`metadata/sample_class_corrections.tsv`",
           paste0("all ", n_class_corrected, " corrected rows with the preserved correction",
                  " record and its SHA256"))),
  md_row(c("`editor_source_data/REVISION_PROTEOMICS_DATA_CHANGELOG.md`",
           "section 11, the full evidence and what is not claimed")),
  "",
  paste0("Only sample-class metadata were reassigned: acquisition identifiers and",
         " quantitative protein-abundance profiles were unchanged, so no deposited raw file",
         " or abundance value is affected. Status `",
         RELEASE_SAMPLE_CLASS_CORRECTION$status, "`."),
  "",
  "## Fields that ARE established, for contrast",
  "",
  "| field | value | evidence |",
  "|---|---|---|",
  md_row(c("`comment[data file]`", "`<run>.d`, 96 unique",
           "`sample_annotation.xlsx` `Name` column; basename minus `.d` equals the sample id for all 96, asserted at build time")),
  md_row(c("`characteristics[organism]`", "Mus musculus",
           "`MOUSE_10090_idmapping.dat`; `organism_id = 10090` in the mapping stage; `org.Mm.eg.db`; `_MOUSE` entry names")),
  md_row(c("`comment[acquisition date]`", "2024-12-17",
           "date token present in all 96 acquisition run names")),
  md_row(c("`characteristics[individual]`", "12 AnimalIDs",
           "validated animal-level design; 3 animals per sample class per condition")),
  ""
)
release_write_lines(missing_md, release_path("pride", "SDRF_MISSING_METADATA.md"))

readme_pride <- c(
  "# PRIDE deposition notes",
  "",
  paste0("**Current status: `", pride_status, "`**"),
  "",
  "## What is suitable for PRIDE, and what is not",
  "",
  "| Artefact | Suitable for PRIDE as | Notes |",
  "|---|---|---|",
  md_row(c("The 96 `.d` acquisition directories", "**RAW files**",
           paste0("Named in `metadata/sample_metadata.tsv` (`raw_file_basename`) and in ",
                  "`pride/sdrf.tsv` (`comment[data file]`). ",
                  ifelse(raw_files_present, "Located on disk.",
                         "**Not yet located** -- these must be retrieved before submission.")))),
  md_row(c("`pride/sdrf.tsv`", "**SDRF-Proteomics metadata file**",
           "Incomplete: see `SDRF_MISSING_METADATA.md`.")),
  md_row(c(paste("Search / quantification output -- run log, report table and spectral",
                 "library (e.g. DIA-NN `report.log.txt` / `report.tsv` /",
                 "`report.parquet` / `.speclib`, if DIA-NN was used)"),
           "**SEARCH output**",
           "**Not present in the project tree.** Only a downstream protein-group matrix survives.")),
  md_row(c("`pg.matrix_raw.txt` (5,747 x 96)", "SEARCH-derived result file",
           paste("Depositable as a processed result. It is **NOT** native raw",
                 "mass-spectrometry data -- see below."))),
  md_row(c("`processed_data/protein_abundance_measurement_level.tsv.gz`",
           "processed result", "Filtered, log2, median-centred, imputed, PCA-adjusted.")),
  md_row(c("`processed_data/protein_abundance_animal_level.tsv.gz`", "processed result",
           "The matrix all inference is computed on.")),
  md_row(c("`differential_analysis/` and `enrichment/`", "processed result",
           "Downstream statistics; usually supplied to the journal rather than to PRIDE.")),
  "",
  "## The naming trap in the original submission",
  "",
  "The file supplied to the journal as `processed_protein_group_matrix_raw.txt` is",
  "**not native raw mass-spectrometry data**. It is a protein-group matrix produced by the",
  "search/quantification software: 5,747 protein groups x 96 quantitative columns, plus",
  "seven annotation columns (`T: Protein.Group`, `T: Protein.Names`, `T: Genes`,",
  "`T: First.Protein.Description`, and three GO/pathway slim columns). In this project tree",
  "the same file is `pg.matrix_raw.txt`. The word \"raw\" in its name refers to it being",
  "pre-filtering, not to it being an instrument acquisition file.",
  "",
  "A PRIDE submission that offered this file as the raw data would be incomplete. The",
  "genuinely raw data are the 96 `.d` acquisition directories.",
  "",
  "## Layers, so the distinction stays visible",
  "",
  "| # | Layer | Artefact | In this package? |",
  "|---|---|---|---|",
  md_row(c("1", "Native acquisition", "96 `.d` directories",
           ifelse(raw_files_present, "located, not copied here", "**not located**"))),
  md_row(c("2", "Search / quantification output",
           "run log / report table / spectral library (software UNKNOWN)",
           "**absent from the project**")),
  md_row(c("3", "Protein-group matrix, pre-filter", "`pg.matrix_raw.txt` (5,747 x 96)",
           "cited; annotation republished")),
  md_row(c("4", "Processed measurement-level matrix", "5,349 x 96", "yes")),
  md_row(c("5", "Animal-level matrix", "5,349 x 48", "yes")),
  md_row(c("6", "Differential statistics", "12 primary comparisons", "yes")),
  md_row(c("7", "Enrichment / EWCE", "GSEA, ORA, EWCE", "yes")),
  md_row(c("8", "Figure source data", "15 panels", "yes")),
  "",
  "## Sample-class metadata is resolved",
  "",
  paste0("For ", n_class_corrected, " of the 96 acquisitions the sample class was corrected",
         " during historical sample-identity quality control. This is **resolved**: it is",
         " not an open question and not a deposition blocker. The SDRF uses the",
         " analysis-time class as the primary factor; both the original and the",
         " analysis-time labels are retained in `metadata/sample_metadata.tsv` and",
         " `metadata/sample_class_corrections.tsv`. Only metadata were reassigned -- no",
         " acquisition identifier and no abundance value changed. See",
         " `SDRF_MISSING_METADATA.md` for the detail."),
  "",
  "## Before submitting",
  "",
  "1. Retrieve the 96 `.d` acquisition directories from the facility instrument store.",
  paste("2. Obtain the original search/quantification software run log or configuration",
        "(e.g. DIA-NN `report.log.txt`, if DIA-NN was used) and the instrument method for",
        "the 2024-12-17 acquisition;"),
  "   they supply most of the missing SDRF fields in one document.",
  paste0("3. Populate the ", n_missing_required, " SDRF-required fields listed in ",
         "`SDRF_MISSING_METADATA.md`."),
  "4. Confirm the dissected brain region and add it as `characteristics[organism part]`,",
  "   and map each of the four sample classes to a Cell Ontology term for",
  "   `characteristics[cell type]`.",
  "5. Re-run `10_validate_release.R` and confirm the status advances.",
  "",
  "## Status vocabulary",
  "",
  "| status | meaning |",
  "|---|---|",
  md_row(c("`PRIDE_READY`", "raw files in hand and every required SDRF field populated")),
  md_row(c("`PRIDE_READY_PENDING_RAW_UPLOAD`", "metadata complete; raw files still to be uploaded")),
  md_row(c("`PRIDE_METADATA_INCOMPLETE`", "raw files identified, but required SDRF fields are unknown")),
  md_row(c("`PRIDE_INPUTS_INCOMPLETE`", "the raw acquisition files are not even identified")),
  "",
  paste0("This deposition is `", pride_status, "`: the acquisition filenames are fully ",
         "identified, but ", n_missing_required, " SDRF-required fields cannot be ",
         "populated from any file in the project."),
  ""
)
release_write_lines(readme_pride, release_path("pride", "README_PRIDE.md"))

# --------------------------------------------------------------------------------------
# write
# --------------------------------------------------------------------------------------

sm_path <- rel("metadata", "sample_metadata.tsv")
sa_path <- file.path(PROJECT_ROOT, "sample_annotation.xlsx")

w1 <- release_write_table(sdrf, release_path("pride", "sdrf.tsv"))
release_register("pride/sdrf.tsv", "SDRF-Proteomics sample and data relationship file",
                 c(sm_path, sa_path), c(NA_character_, release_sha256(sa_path)), STAGE, "tsv")

w2 <- release_write_table(field_status, release_path("pride", "sdrf_field_status.tsv"))
release_register("pride/sdrf_field_status.tsv",
                 "per-SDRF-field readiness: READY / MISSING_REQUIRED / MISSING_OPTIONAL",
                 sm_path, NA_character_, STAGE, "tsv")

release_register("pride/SDRF_MISSING_METADATA.md",
                 "SDRF fields that could not be established, and what would supply them",
                 sm_path, NA_character_, STAGE, "md")
release_register("pride/README_PRIDE.md",
                 "which artefacts are suitable for PRIDE and which are not",
                 sm_path, NA_character_, STAGE, "md")

status_tbl <- data.frame(
  key = c("pride_status", "raw_files_identified", "raw_files_located",
          "n_raw_files_located", "n_sdrf_fields_ready", "n_sdrf_required_missing",
          "n_sdrf_optional_missing", "sample_class_correction_status",
          "n_sample_class_corrected"),
  value = c(pride_status, as.character(raw_files_identified), as.character(raw_files_present),
            as.character(raw_present_count), as.character(n_ready),
            as.character(n_missing_required), as.character(n_missing_optional),
            RELEASE_SAMPLE_CLASS_CORRECTION$status, as.character(n_class_corrected)),
  stringsAsFactors = FALSE, check.names = FALSE)
release_write_table(status_tbl, release_path("pride", "pride_readiness.tsv"))
release_register("pride/pride_readiness.tsv", "machine-readable PRIDE readiness verdict",
                 sm_path, NA_character_, STAGE, "tsv")

release_log("  wrote sdrf.tsv (", w1$rows, "x", w1$cols, ")")
release_log("  wrote sdrf_field_status.tsv (", w2$rows, "x", w2$cols, ")")
release_log("  wrote SDRF_MISSING_METADATA.md, README_PRIDE.md, pride_readiness.tsv")
release_log("stage 09 complete")
