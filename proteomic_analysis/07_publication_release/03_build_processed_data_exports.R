#!/usr/bin/env Rscript

# Publication release, stage 03 -- processed abundance matrices and feature annotation.
#
# Produces
#   processed_data/protein_abundance_measurement_level.tsv.gz   5349 x 96 acquisitions
#   processed_data/protein_abundance_animal_level.tsv.gz        5349 x 48 animal units
#   processed_data/protein_feature_annotation.tsv.gz            5349 protein groups
#
# NAMING. Nothing here is called "raw". The measurement-level matrix is processed
# protein-level output, standardised per protein (mean 0, sd 1 across acquisitions); the
# genuinely raw acquisition files are the `.d` directories named in
# metadata/sample_metadata.tsv and are NOT part of this package.
#
# FIELD SEMANTICS. The internal pipeline carries two misleading column names, and this
# stage exists partly to stop them reaching a published table:
#
#   split/*.csv     `gene_symbol`  actually holds the UniProt ENTRY NAME(s)  (A0A0J9YTR2_MOUSE)
#                   `Description`  actually holds the GENE SYMBOL(s)         (Ncbp2as2)
#   mapped/*.csv    `gene_symbol`  actually holds the UniProt ACCESSION      (A0A0J9YTR2)
#
# The published names are therefore protein_group_id / uniprot_accession / gene_symbol /
# protein_description, each carrying what its name says. A real protein description does
# exist, but only in the search output (pg.matrix_raw.txt, DIA-NN column
# `T: First.Protein.Description`) -- so it is joined in from there rather than faked from
# the gene symbol.

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
STAGE <- "07_publication_release/03_build_processed_data_exports.R"

release_banner("stage 03 -- processed data exports")

ANIMAL_ROOT <- file.path(DATA_ROOT, "02_data", "animal_level")

paths <- list(
  measurement_gct = file.path(DATA_ROOT, "02_data", "gct",
                              "pg.matrix_filtered_pcaAdjusted_unnormalized.gct"),
  animal_gct = release_locked_artefact_path("animal_level_input_gct", DATA_ROOT),
  search_output = file.path(PROJECT_ROOT, "pg.matrix_raw.txt"),
  mapped_forward_dir = file.path(ANIMAL_ROOT, "mapped", "forward"),
  unmapped_dir = file.path(ANIMAL_ROOT, "mapped", "unmapped"),
  mapped_index = file.path(ANIMAL_ROOT, "mapped", "indexMappedComparisons.csv"),
  sample_metadata = release_path("metadata", "sample_metadata.tsv", create_dir = FALSE),
  animal_metadata = release_path("metadata", "animal_level_sample_metadata.tsv",
                                 create_dir = FALSE)
)
for (nm in setdiff(names(paths), c("mapped_forward_dir", "unmapped_dir"))) {
  release_assert_exists(paths[[nm]], nm)
}

sample_metadata <- release_read_tsv_plain(paths$sample_metadata)
animal_metadata <- release_read_tsv_plain(paths$animal_metadata)
mapped_index <- release_read_csv(paths$mapped_index)

# --------------------------------------------------------------------------------------
# measurement-level abundance matrix
# --------------------------------------------------------------------------------------

meas <- release_read_gct(paths$measurement_gct)
release_log("  measurement-level GCT: ", meas$dimensions[["n_row"]], " x ",
            meas$dimensions[["n_col"]])
if (meas$dimensions[["n_row"]] != RELEASE_DESIGN_INVARIANTS$n_proteins_statistical ||
    meas$dimensions[["n_col"]] != RELEASE_DESIGN_INVARIANTS$n_measurement_records) {
  stop("Measurement-level matrix is not 5349 x 96.", call. = FALSE)
}
if (!setequal(meas$sample_ids, sample_metadata$sample_id)) {
  stop("Measurement-level matrix columns do not match the published sample metadata.",
       call. = FALSE)
}

# What scale are these values actually on? The filename says `pcaAdjusted_unnormalized`,
# which suggests neither, and a package that mislabels its own units is worse than one that
# says it does not know. Measured directly: every protein has mean 0 and sd 1 across the 96
# acquisitions, i.e. the matrix is standardised PER PROTEIN. It is therefore not log2
# abundance, and a difference of group means on it is in units of per-protein SD.
# Asserted here so the description this package publishes cannot drift from the data.
row_mean <- rowMeans(meas$mat)
row_sd <- apply(meas$mat, 1L, stats::sd)
is_row_standardised <- max(abs(row_mean)) < 1e-2 && max(abs(row_sd - 1)) < 5e-2
release_log("  value scale: per-protein mean ", signif(stats::median(row_mean), 3),
            ", per-protein sd ", signif(stats::median(row_sd), 4),
            " -> row-standardised = ", is_row_standardised)
if (!is_row_standardised) {
  stop("The measurement-level matrix is no longer per-protein standardised. The release ",
       "describes these values as standardised; refusing to publish a description the ",
       "data contradict.", call. = FALSE)
}
VALUE_SCALE <- paste(
  "per-protein standardised (z-scored) abundance: each protein has mean 0 and standard",
  "deviation 1 across the 96 acquisitions. NOT log2 abundance and NOT raw intensity.")

meas_order <- match(sample_metadata$sample_id, meas$sample_ids)
measurement_export <- data.frame(
  protein_group_id = meas$ids,
  meas$mat[, meas_order, drop = FALSE],
  stringsAsFactors = FALSE, check.names = FALSE
)
names(measurement_export)[-1] <- sample_metadata$sample_id

# --------------------------------------------------------------------------------------
# animal-level abundance matrix
# --------------------------------------------------------------------------------------

animal <- release_read_gct(paths$animal_gct)
release_log("  animal-level GCT: ", animal$dimensions[["n_row"]], " x ",
            animal$dimensions[["n_col"]])
if (animal$dimensions[["n_row"]] != RELEASE_DESIGN_INVARIANTS$n_proteins_statistical ||
    animal$dimensions[["n_col"]] != RELEASE_DESIGN_INVARIANTS$n_animal_level_units) {
  stop("Animal-level matrix is not 5349 x 48.", call. = FALSE)
}
if (!identical(animal$ids, meas$ids)) {
  stop("Animal-level and measurement-level matrices do not carry the same protein rows ",
       "in the same order.", call. = FALSE)
}
if (!setequal(animal$sample_ids, animal_metadata$animal_level_column_name)) {
  stop("Animal-level matrix columns do not match the published animal-level metadata.",
       call. = FALSE)
}

animal_order <- match(animal_metadata$animal_level_column_name, animal$sample_ids)
animal_export <- data.frame(
  protein_group_id = animal$ids,
  animal$mat[, animal_order, drop = FALSE],
  stringsAsFactors = FALSE, check.names = FALSE
)
names(animal_export)[-1] <- animal_metadata$animal_level_column_name

# --------------------------------------------------------------------------------------
# feature annotation
# --------------------------------------------------------------------------------------
# Three independent sources are joined, each contributing only what it actually holds:
#   the animal-level GCT row descriptor  -> the source matrix "Description" field
#   pg.matrix_raw.txt (search output)    -> accessions, all gene symbols, real description
#   mapped/ + unmapped/ tables           -> the resolved symbol, accession, mapping status

search_out <- release_read_tsv_plain(paths$search_output, "search output pg.matrix_raw.txt")
col_entry_names <- "T: Protein.Names"
col_accessions <- "T: Protein.Group"
col_genes <- "T: Genes"
col_desc <- grep("^T: First[.]Protein", names(search_out), value = TRUE)
if (length(col_desc) != 1L) {
  # The header in the archived search output carries a mangled name
  # ("T: First.Proteinescription"); match it defensively rather than hardcoding the typo.
  col_desc <- grep("^T: First", names(search_out), value = TRUE)
}
for (needed in c(col_entry_names, col_accessions, col_genes)) {
  if (!needed %in% names(search_out)) {
    stop("Search output is missing the expected DIA-NN column: ", needed, call. = FALSE)
  }
}
if (length(col_desc) != 1L) {
  stop("Could not identify the protein-description column in the search output.", call. = FALSE)
}
release_log("  search output: ", nrow(search_out), " protein groups (pre-filter), ",
            ncol(search_out), " columns")

strip_quotes <- function(x) gsub('^"|"$', "", trimws(as.character(x)))
search_entry <- strip_quotes(search_out[[col_entry_names]])
search_idx <- match(meas$ids, search_entry)
if (anyNA(search_idx)) {
  stop("Search output does not cover every analysed protein group (",
       sum(is.na(search_idx)), " unmatched).", call. = FALSE)
}

# Resolved mapping: read one mapped forward table plus its unmapped companion. The
# mapping is a property of the protein, identical across the 12 comparisons; the index
# records identical mapping statistics for all of them, which is asserted below.
if (length(unique(mapped_index$n_output_mapped_rows)) != 1L ||
    length(unique(mapped_index$n_unmapped)) != 1L) {
  stop("Mapping outcome differs between comparisons; a single annotation table would ",
       "be wrong.", call. = FALSE)
}
representative <- mapped_index$canonical_contrast[[1]]
mapped_tbl <- release_read_csv(file.path(paths$mapped_forward_dir,
                                         paste0(representative, ".csv")))
unmapped_tbl <- release_read_csv(file.path(paths$unmapped_dir,
                                           paste0(representative, ".csv")))
mapping_all <- rbind(
  mapped_tbl[c("original_protein_id", "mapped_gene_symbol", "uniprot_accession",
               "mapping_status", "mapping_strategy", "multi_protein", "source_row_id")],
  unmapped_tbl[c("original_protein_id", "mapped_gene_symbol", "uniprot_accession",
                 "mapping_status", "mapping_strategy", "multi_protein", "source_row_id")]
)
if (nrow(mapping_all) != RELEASE_DESIGN_INVARIANTS$n_proteins_statistical) {
  stop("mapped + unmapped rows do not account for all 5349 protein groups (got ",
       nrow(mapping_all), ").", call. = FALSE)
}
map_idx <- match(meas$ids, mapping_all$original_protein_id)
if (anyNA(map_idx)) stop("Mapping tables do not cover every protein group.", call. = FALSE)

# Sanity: the resolved gene symbols must NOT look like UniProt identifiers, and the
# accession column must. If this ever inverts, the misleading internal names have leaked.
resolved_symbol <- as.character(mapping_all$mapped_gene_symbol)[map_idx]
resolved_accession <- as.character(mapping_all$uniprot_accession)[map_idx]
if (release_column_is_misleading_gene_symbol(resolved_symbol)) {
  stop("Resolved gene_symbol column looks like UniProt identifiers; refusing to publish.",
       call. = FALSE)
}
if (!release_column_is_misleading_gene_symbol(resolved_accession)) {
  stop("uniprot_accession column does not look like UniProt accessions.", call. = FALSE)
}

source_description_field <- if (!is.null(animal$rdesc) && "Description" %in% names(animal$rdesc)) {
  as.character(animal$rdesc$Description)
} else {
  rep(NA_character_, length(animal$ids))
}

entry_names <- meas$ids
n_in_group <- lengths(strsplit(entry_names, ";", fixed = TRUE))

feature_annotation <- data.frame(
  protein_group_id = entry_names,
  protein_group_uniprot_accessions = strip_quotes(search_out[[col_accessions]])[search_idx],
  uniprot_accession = resolved_accession,
  gene_symbol = resolved_symbol,
  gene_symbols_in_group = strip_quotes(search_out[[col_genes]])[search_idx],
  protein_description = strip_quotes(search_out[[col_desc]])[search_idx],
  n_proteins_in_group = as.integer(n_in_group),
  is_protein_group = n_in_group > 1L,
  id_mapping_status = as.character(mapping_all$mapping_status)[map_idx],
  id_mapping_strategy = as.character(mapping_all$mapping_strategy)[map_idx],
  id_mapping_multi_protein = as.character(mapping_all$multi_protein)[map_idx],
  source_matrix_row_index = as.integer(mapping_all$source_row_id)[map_idx],
  source_matrix_description_field = source_description_field,
  id_mapping_reference_file = basename(as.character(mapped_index$mapping_reference_file[[1]])),
  id_mapping_reference_version = as.character(mapped_index$mapping_reference_version[[1]]),
  id_mapping_reference_snapshot_date_utc =
    as.character(mapped_index$mapping_reference_snapshot_date_utc[[1]]),
  id_mapping_reference_sha256 = as.character(mapped_index$mapping_reference_sha256[[1]]),
  stringsAsFactors = FALSE, check.names = FALSE
)

if (nrow(feature_annotation) != RELEASE_DESIGN_INVARIANTS$n_proteins_statistical) {
  stop("Feature annotation is not 5349 rows.", call. = FALSE)
}
n_mapped <- sum(feature_annotation$id_mapping_status == "mapped")
if (n_mapped != RELEASE_DESIGN_INVARIANTS$n_proteins_mapped) {
  stop("Feature annotation reports ", n_mapped, " mapped proteins, expected 5327.",
       call. = FALSE)
}
release_log("  feature annotation: 5349 protein groups, ", n_mapped, " mapped, ",
            nrow(feature_annotation) - n_mapped, " unmapped")
release_log("  protein_description populated for ",
            sum(nzchar(feature_annotation$protein_description) &
                  !is.na(feature_annotation$protein_description)), " groups")

# --------------------------------------------------------------------------------------
# write
# --------------------------------------------------------------------------------------

meas_hash <- release_sha256(paths$measurement_gct)
animal_hash <- release_sha256(paths$animal_gct)
search_hash <- release_sha256(paths$search_output)

w1 <- release_write_table(measurement_export,
                          release_path("processed_data",
                                       "protein_abundance_measurement_level.tsv.gz"))
release_register("processed_data/protein_abundance_measurement_level.tsv.gz",
                 paste("processed protein-level abundance, one column per acquisition run;",
                       VALUE_SCALE, "; NOT raw MS data"),
                 paths$measurement_gct, meas_hash, STAGE, "tsv.gz")

w2 <- release_write_table(animal_export,
                          release_path("processed_data",
                                       "protein_abundance_animal_level.tsv.gz"))
release_register("processed_data/protein_abundance_animal_level.tsv.gz",
                 paste("animal-level abundance used for all inference; Left/Right averaged",
                       "within AnimalID x sample_class; ", VALUE_SCALE),
                 paths$animal_gct, animal_hash, STAGE, "tsv.gz")

w3 <- release_write_table(feature_annotation,
                          release_path("processed_data",
                                       "protein_feature_annotation.tsv.gz"))
release_register("processed_data/protein_feature_annotation.tsv.gz",
                 paste("protein-group identifiers, accessions, resolved gene symbols and",
                       "real protein descriptions with mapping provenance"),
                 c(paths$search_output, paths$mapped_index, paths$animal_gct),
                 c(search_hash, release_sha256(paths$mapped_index), animal_hash),
                 STAGE, "tsv.gz")

release_log("  wrote protein_abundance_measurement_level.tsv.gz (", w1$rows, "x", w1$cols, ")")
release_log("  wrote protein_abundance_animal_level.tsv.gz (", w2$rows, "x", w2$cols, ")")
release_log("  wrote protein_feature_annotation.tsv.gz (", w3$rows, "x", w3$cols, ")")
release_log("stage 03 complete")
