# ==============================================================================
# OUT OF SCOPE FOR THIS PROJECT -- Exp9_Social-Stress, not this project.
#
# This script reads TPE9_* workbooks from the Exp9_Social-Stress project and
# WRITES two files back into that project's folder:
#     TPE9_samples_males_processed.tsv
#     TPE9_samples_males_long_with_metadata.xlsx
# It reads none of this project's input and produces none of its output. It was inherited when this
# repository was seeded from Exp9 (the same origin as the quicksearch.stats
# default described in CANONICAL_OUTPUTS.md).
#
# It used to sit at 01_preprocessing/05_metadata_create.r, inside the numbered
# active stage sequence, with no path override and no guard. Every precondition
# for an accidental cross-project overwrite was satisfied on the shared drive:
# the Exp9 folder exists, both inputs exist, and TPE9_samples_males_processed.tsv
# already exists, so a direct `Rscript` invocation would have silently replaced
# a real file belonging to another experiment.
#
# It is therefore quarantined here and refuses to run unless you opt in:
#     PROTEOMICS_ALLOW_EXP9=true Rscript 99_out_of_scope/05_metadata_create_EXP9.r
# and the target directory can be redirected away from the live Exp9 folder:
#     PROTEOMICS_EXP9_WORK_DIR=/some/scratch/dir
#
# run_pipeline_check.ps1 already records this stage as SKIP with the reason
# "belongs to a different project (Exp9_Social-Stress)"; the guard below makes
# that true for direct invocation as well.
# ==============================================================================

option_or_env <- function(option_name, env_name, default) {
    value <- getOption(option_name)
    if (!is.null(value) && nzchar(trimws(as.character(value)))) return(as.character(value))
    value <- Sys.getenv(env_name, unset = "")
    if (nzchar(trimws(value))) return(value)
    default
}

allow <- tolower(trimws(option_or_env("proteomics.allow_exp9", "PROTEOMICS_ALLOW_EXP9", "")))
if (!allow %in% c("true", "1", "yes")) {
    stop(
        "Refusing to run: this script belongs to Exp9_Social-Stress, not this project.\n",
        "It writes TPE9_samples_males_processed.tsv and ",
        "TPE9_samples_males_long_with_metadata.xlsx into the Exp9 project folder,\n",
        "overwriting the former if it exists.\n\n",
        "If you really intend to run it, opt in explicitly:\n",
        "  PROTEOMICS_ALLOW_EXP9=true Rscript 99_out_of_scope/05_metadata_create_EXP9.r\n",
        "and consider redirecting the target away from the live Exp9 folder:\n",
        "  PROTEOMICS_EXP9_WORK_DIR=/some/scratch/dir",
        call. = FALSE
    )
}

# read in xlsx files from directory

if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
}
pacman::p_load(readxl, dplyr, tidyr, stringr, purrr, tibble, writexl)

work_direction <- option_or_env(
    "proteomics.exp9_work_dir", "PROTEOMICS_EXP9_WORK_DIR",
    "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
)
message("Exp9 working directory: ", work_direction)
metadata <- readxl::read_xlsx(file.path(work_direction, "TPE9_sample_metadata_males.xlsx"))

files <- list.files(path = work_direction, pattern = "TPE9_samples_males\\.xlsx$", full.names = TRUE)
if (length(files) == 0) {
    stop("No matching .xlsx files found in work_direction: ", work_direction)
}
safe_read <- purrr::possibly(readxl::read_xlsx, otherwise = tibble::tibble())
data <- purrr::map_dfr(files, safe_read) %>%
    dplyr::filter(!is.na(`Protein.Group`))

# save data as filename but with _processed.tsv
output_file <- file.path(work_direction, "TPE9_samples_males_processed.tsv")
# write as tab-separated values, no row names, no quotes, empty strings for NA
write.table(data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# reshape protein table to long and join with metadata, then save
keep_cols <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description")

# identify sample columns (those that are not protein annotation columns)
sample_cols <- setdiff(names(data), keep_cols)

# Prepare a conservative "cleaning" function so column names and metadata IDs can be matched
clean_id <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    # drop common prefixes/suffixes you may have in column names, adjust if needed
    str_remove_all("^Intensity\\.|^X|^Sample[_-]") %>%
    str_remove_all("\\.raw$|\\.mzML$|_1$|_2$") %>%
    str_replace_all("[^A-Za-z0-9_-]+", "_") %>%
    str_to_lower()
}

# add cleaned IDs to metadata so we can match robustly
if (!"sample_id" %in% names(metadata)) {
  stop("metadata must contain a column named 'sample_id' to match against sample columns")
}
metadata <- metadata %>% mutate(sample_id_clean = clean_id(sample_id))

# create mapping from original column names to cleaned sample_id
sample_map <- tibble(original = sample_cols) %>%
  mutate(sample_id_clean = clean_id(original))

# pivot the sample columns to long form using the original column names, then attach cleaned id and metadata
long_proteins <- data %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "original",
    values_to = "intensity"
  ) %>%
  mutate(intensity = as.numeric(intensity)) %>%
  left_join(sample_map, by = "original") %>%             # add cleaned sample id
  left_join(metadata, by = "sample_id_clean") %>%        # add metadata columns
  filter(!is.na(intensity))

writexl::write_xlsx(long_proteins, file.path(work_direction, "TPE9_samples_males_long_with_metadata.xlsx"))