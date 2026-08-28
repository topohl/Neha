# Install and load required packages
packages <- c("readxl", "writexl", "dplyr")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
    install.packages(packages[!installed])
}
lapply(packages, library, character.only = TRUE)
source(file.path("R", "analysis_labels.R"))

# Define the file path (overridable; default follows the 2026-08-26 restructure)
option_or_env <- function(option_name, env_name, default) {
    value <- getOption(option_name)
    if (!is.null(value) && nzchar(trimws(as.character(value)))) return(as.character(value))
    value <- Sys.getenv(env_name, unset = "")
    if (nzchar(trimws(value))) return(value)
    default
}
project_root <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
file_path <- option_or_env(
    "proteomics.format_metadata_input", "PROTEOMICS_FORMAT_METADATA_INPUT",
    file.path(project_root, "01_input", "metadata", "sample_metadata.xlsx")
)
if (!file.exists(file_path)) {
    stop(
        "Sample metadata workbook not found: ", file_path,
        "\nThe former Datasets/sample_metadata/ folder is not present in the current tree.",
        "\nSet PROTEOMICS_FORMAT_METADATA_INPUT / options(proteomics.format_metadata_input=) to its location.",
        call. = FALSE
    )
}

# Get available sheet names
sheet_names <- excel_sheets(file_path)
print(sheet_names)

# Load the sheet named "samples" into a data frame
df <- read_excel(file_path, sheet = "samples")
head(df)

# extract info from column "sample_id".
# A0003 is AnimalID, _L_ or _R_ indicate left or right replicates and go into "ReplicateGroup".
# _13026 is the run number and goes into "run_order"
df <- df %>%
  mutate(
    AnimalID = sub(".*_(A[0-9]+)_.*", "\\1", sample_id),
    ReplicateGroup = ifelse(grepl("_L_", sample_id), "Left", "Right"),
    run_order = sub(".*_(\\d+)\\.d$", "\\1", sample_id),
    raw_sample_token = sub(".*_(?:A[0-9]+)_(?:L|R)_([^_]+(?:_[^_]+)*)_S\\d+_.*", "\\1", sample_id),
    condition_code = parse_condition_code(sample_id),
    condition = normalize_condition(condition_code),
    exclude = ifelse(grepl("CA3_slm", sample_id), TRUE, FALSE),
    plate = ifelse(grepl("20250703", sample_id), "B",
             ifelse(grepl("20250707", sample_id), "C", NA)),
    sample_number = sub(".*_(S\\d+)_.*", "\\1", sample_id),
    sample_location = ifelse(grepl("_[S]\\d+-", sample_id),
                             sub(".*_(S\\d+)-([A-Z]\\d+)_.*", "\\2", sample_id),
                             sub(".*_(S\\d+)_.*", "\\1", sample_id)),
    shortname = paste(plate, sample_number, sample_location, run_order, sep="_"),
    replicate_unit = paste(condition, ReplicateGroup, sep = "_"),
    sample_class = parse_sample_class(sample_id)
  )

# save as excel file using writexl package
output_path <- option_or_env(
    "proteomics.format_metadata_output", "PROTEOMICS_FORMAT_METADATA_OUTPUT",
    file.path(project_root, "01_input", "metadata", "sample_metadata_processed.xlsx")
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_xlsx(df, output_path)
