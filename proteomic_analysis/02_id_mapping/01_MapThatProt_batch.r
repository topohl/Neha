#' Batch UniProt ID Mapping Script for Proteomics Data
#'
#' This script processes the indexed animal-level ProTigy comparison files
#' in parallel to map protein identifiers, UniProtKB entry names, and aliases to
#' canonical UniProt accessions. It supports
#' multiple mapping strategies including local mapping files, OrgDb annotations,
#' UniProt.ws queries, and manual overrides. The script enforces a _MOUSE filter
#' to ensure species-specific mapping and outputs mapped/unmapped datasets with
#' detailed mapping strategy reports.
#'
#' Key features:
#' - Parallel processing of multiple CSV files using doParallel
#' - Multi-strategy mapping cascade (accession detection, entry name lookup,
#'   gene symbol resolution, family prefix matching, OrgDb, UniProt.ws)
#' - Manual mapping override support via Excel file
#' - Comprehensive mapping statistics and unmapped protein tracking
#' - Species-specific mapping eligibility (_MOUSE suffix enforcement)
#' - Indexed 12-contrast discovery and isolated animal-level output provenance
#' 
#' Dependencies: dplyr, stringr, tidyr, purrr, readr, R.utils, foreach, doParallel, readxl, AnnotationDbi, org.Mm.eg.db, UniProt.ws
#' 
#' @title MapThatProt_batch
#' 
#' @author Tobias Pohl
#' 
#' @date 2025-12-15

cat("====================================================\n")
cat("Starting MapThatProt_batch execution...\n")
cat("====================================================\n")

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) == 1L) {
    sub("^--file=", "", file_arg)
} else {
    file.path("02_id_mapping", "01_MapThatProt_batch.r")
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
source(file.path(repo_root, "R", "project_path_utils.R"))
if (!file.exists(file.path(repo_root, "R", "mapthatprot_animal_level_utils.R"))) {
    repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(repo_root, "R", "analysis_labels.R"))
source(file.path(repo_root, "R", "mapthatprot_animal_level_utils.R"))

# --- Package Management ---
# Automatically install and load necessary CRAN and Bioconductor packages.
cat("Checking and loading required libraries...\n")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, stringr, tidyr, purrr, readr, R.utils, foreach, doParallel, readxl, AnnotationDbi, org.Mm.eg.db, UniProt.ws)

# --- Configuration & Experimental Settings ---
mapping_config <- resolve_mapthatprot_config()
map_direction <- mapping_config$direction
map_reverse <- identical(map_direction, "reverse")
input_manifest <- read_mapthatprot_input_index(
    mapping_config$split_root,
    direction = map_direction
)
csv_files <- input_manifest$input_path
split_index_provenance <- data.frame(
    source_split_index_path = input_manifest$split_index_path[[1]],
    source_split_index_sha256 = mapthatprot_sha256(input_manifest$split_index_path[[1]]),
    stringsAsFactors = FALSE
)
input_manifest$source_split_sha256 <- vapply(csv_files, mapthatprot_sha256, character(1))
input_manifest$validated_input_rows <- vapply(seq_along(csv_files), function(i) {
    expected_rows <- suppressWarnings(as.integer(input_manifest$n_proteins[[i]]))
    if (is.na(expected_rows)) expected_rows <- 5349L
    nrow(validate_mapthatprot_input_file(csv_files[[i]], expected_n_proteins = expected_rows))
}, integer(1))

cat("Mapping branch: canonical animal-level ProTigy split\n")
cat("Mapping direction:", map_direction, "\n")
cat("Authoritative split index:", input_manifest$split_index_path[[1]], "\n")
cat("Indexed comparison files:", length(csv_files), "\n")
cat("Output root:", mapping_config$output_root, "\n")

# Define output directories for mapped datasets, unmapped trackers, and QC info
mapped_dir <- file.path(mapping_config$output_root, map_direction)
info_dir <- file.path(mapping_config$output_root, "info", if (map_reverse) "reverse" else "")
mapped_summary_dir <- file.path(info_dir, "mapped_summaries")
unmapped_dir <- file.path(mapping_config$output_root, "unmapped", if (map_reverse) "reverse" else "")
unmapped_summary_dir <- file.path(info_dir, "unmapped_summaries")
report_dir <- file.path(mapping_config$output_root, "mapping_reports", if (map_reverse) "reverse" else "")
mapped_index_path <- file.path(
    mapping_config$output_root,
    if (map_reverse) "indexMappedComparisons_reverse.csv" else "indexMappedComparisons.csv"
)

# --- Reference Databases ---
# Define path for central UniProt species-specific knowledgebase flatfile
uniprot_mapping_file_path <- mapping_config$uniprot_mapping_file
if (!file.exists(uniprot_mapping_file_path)) {
    stop("Required local mouse UniProt mapping reference is missing: ", uniprot_mapping_file_path, call. = FALSE)
}
mapping_reference_provenance <- mapthatprot_reference_provenance(uniprot_mapping_file_path)

# Parse the UniProt mapping dictionary natively into memory
cat("Parsing UniProt idmapping dictionary into memory... (This may take a moment)\n")
uniprot_mapping <- readr::read_tsv(
    uniprot_mapping_file_path,
    col_names = c("UniProt_Accession", "Type", "Value"),
    col_types = "ccc",
    quote = ""
)

# Extract core UniProtKB canonical accessions to Entry name mapping
entry_name_to_accession <- uniprot_mapping %>%
    filter(Type == "UniProtKB-ID") %>%
    dplyr::select(UniProt_Accession, UniProtKB_ID = Value) %>%
    distinct(UniProtKB_ID, .keep_all = TRUE)

if (nrow(entry_name_to_accession) == 0) stop("No UniProtKB-ID mappings found in mapping file.")

# Canonical accession-to-gene-symbol annotation for explicit mapped identity.
accession_to_gene_symbol <- uniprot_mapping %>%
    dplyr::filter(Type == "Gene_Name") %>%
    dplyr::transmute(
        uniprot_accession = toupper(trimws(UniProt_Accession)),
        mapped_gene_symbol = trimws(Value)
    ) %>%
    dplyr::filter(nzchar(uniprot_accession), nzchar(mapped_gene_symbol)) %>%
    dplyr::distinct(uniprot_accession, .keep_all = TRUE)
accession_gene_symbol_map <- stats::setNames(
    accession_to_gene_symbol$mapped_gene_symbol,
    accession_to_gene_symbol$uniprot_accession
)

# --- Manual Override Configuration ---
# Optional manual curation file. This is crucial for resolving heavily ambiguous 
# protein groups or unannotated gene symbols common in exploratory proteomics.
cat("Checking for manual mapping override file...\n")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")

manual_mapping_path <- mapping_config$manual_mapping_file
manual_override <- TRUE  # TRUE enforces curation over algorithmic mapping

read_manual_xlsx <- function(path) {
    mm <- try(readxl::read_excel(path, sheet = 1), silent = TRUE)
    if (inherits(mm, "try-error") || !is.data.frame(mm) || !nrow(mm)) {
        cat("Notice: Manual mapping Excel unreadable or empty at:", path, "\n")
        return(NULL)
    }
    # Normalize curation column headers
    mm %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), ~ toupper(trimws(as.character(.))))) %>%
        dplyr::rename_with(tolower)
}

manual_mapping <- if (file.exists(manual_mapping_path)) {
    read_manual_xlsx(manual_mapping_path)
} else {
    cat("Notice: Manual mapping Excel not found at:", manual_mapping_path, "\n")
    NULL
}
manual_mapping_provenance <- data.frame(
    manual_mapping_path = normalize_mapthatprot_path(manual_mapping_path),
    manual_mapping_sha256 = if (file.exists(manual_mapping_path)) mapthatprot_sha256(manual_mapping_path) else NA_character_,
    manual_mapping_rows = if (is.null(manual_mapping)) 0L else nrow(manual_mapping),
    stringsAsFactors = FALSE
)

if (!is.null(manual_mapping)) {
    needed <- c("gene_symbol", "mapped_gene_symbol")
    missing_cols <- setdiff(needed, names(manual_mapping))
    if (length(missing_cols)) {
        cat("Warning: Manual mapping is missing expected columns:", paste(missing_cols, collapse = ", "), "\n")
    } else {
        cat("Successfully loaded", nrow(manual_mapping), "rows from manual mapping configuration.\n")
    }
}

# Create only the isolated animal-level output tree, after every input and
# reference contract has passed.
cat("Initializing isolated animal-level mapping output directories...\n")
for (directory in c(
    info_dir, mapped_dir, mapped_summary_dir, unmapped_dir,
    unmapped_summary_dir, report_dir
)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

# --- Core Mapping Function ---
# This function is executed asynchronously for each input proteomics file.
process_file <- function(data_path) {
    # Suppress internal messages for cleaner parallel console, file-level info handles writing
    # Ingest generic CSV arrays, supporting delimiter fallbacks (e.g. European CSVs)
    df_raw <- tryCatch(
        readr::read_csv(data_path, col_names = TRUE, show_col_types = FALSE, trim_ws = TRUE, quote = "\""),
        error = function(e) {
            readr::read_csv2(data_path, col_names = TRUE, show_col_types = FALSE, trim_ws = TRUE)
        }
    )

    required_input <- mapthatprot_required_input_columns()
    missing_input <- setdiff(required_input, names(df_raw))
    if (length(missing_input)) {
        stop("Canonical split input is missing required columns: ", paste(missing_input, collapse = ", "), call. = FALSE)
    }
    source_input <- as.data.frame(df_raw, stringsAsFactors = FALSE, check.names = FALSE)

    # Mark rows that have protein groups (indicated by semicolons)
    df_raw <- df_raw %>%
    dplyr::mutate(
        original_row_id = dplyr::row_number(),
        original_gene_symbol = as.character(gene_symbol),
        original_protein_id = as.character(gene_symbol),
        multi_protein = stringr::str_detect(gene_symbol, ";")
    )

    # --- Bioinformatics String Normalization Functions ---
    # Strip whitespace, illegal characters, and formatting artifacts
    normalize_token <- function(x) {
        x <- as.character(x)
        x <- sub("\\s.*$", "", x)
        x <- gsub("\\|\\|+", "|", x)
        x <- toupper(gsub("\\s+", "", x))
        x <- gsub("\\u00A0", "", x)
        x <- gsub("\\.+", ".", x)
        x <- gsub("__+", "_", x)
        x
    }
    # Decouple isoform numbers (e.g., Q91X11-2 -> Q91X11) and _MOUSE suffixes
    to_base_no_iso_mouse <- function(x) {
        x <- gsub("-\\d+$", "", x)
        gsub("_MOUSE$", "", x)
    }
    # Validate strictly canonical UniProt lengths & formats
    is_uniprot_ac <- function(x) {
        grepl("^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|A0A[0-9A-Z]{7})$", x)
    }
    # Parse likely accessions stringed into fasta headers
    extract_ac <- function(s) {
        s <- as.character(s)
        m <- stringr::str_match(s, "(?i)(?:^|\\|)([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|A0A[0-9A-Z]{7})(?:\\-|\\||$|[^A-Z0-9])")
        out <- ifelse(is.na(m[,2]), NA_character_, toupper(m[,2]))
        gsub("-\\d+$", "", out)
    }
    # Look for Entry Names (e.g. MAPK1_MOUSE)
    extract_entry <- function(s) {
        s <- as.character(s)
        m <- stringr::str_match(s, "(?i)(?:^|\\|)([A-Z0-9]+_MOUSE)(?:\\||$|\\s)")
        out <- ifelse(is.na(m[,2]), NA_character_, m[,2])
        if (all(is.na(out))) {
            m2 <- stringr::str_match(s, "(?i)\\b([A-Z0-9]+_MOUSE)\\b")
            out <- m2[,2]
        }
        toupper(out)
    }
    nz <- function(x) !is.na(x) & nzchar(x)

    # --- Primary Tokenization Pipeline ---
    df_tok_all <- df_raw %>%
        dplyr::mutate(
            gene_symbol_split = stringr::str_split(gene_symbol, ";"),
            n_proteins = purrr::map_int(gene_symbol_split, length),
            kept_symbol = purrr::map_chr(gene_symbol_split, 1), # Keep the leading master protein
            dropped_symbols = purrr::map_chr(
                gene_symbol_split,
                ~ if (length(.) > 1) paste(.[-1], collapse = ";") else NA_character_
            )
        ) %>%
        dplyr::mutate(
            gene_symbol = kept_symbol
        ) %>%
        dplyr::mutate(
            source_row_id = original_row_id,
            token_raw = trimws(as.character(gene_symbol)),
            token_up  = normalize_token(token_raw),
            acc_guess = extract_ac(token_raw),
            entry_guess_up = extract_entry(token_raw)
        ) %>%
        dplyr::mutate(
            # Coerce everything to a stripped logical base identifier
            token_base = dplyr::case_when(
                nz(entry_guess_up) ~ to_base_no_iso_mouse(entry_guess_up),
                TRUE               ~ to_base_no_iso_mouse(token_up)
            ),
            # Stratify identity classes for the resolution cascade mapping
            token_kind = dplyr::case_when(
                nz(acc_guess)                          ~ "AC_GUESS",
                nz(entry_guess_up)                     ~ "ENTRY_GUESS",
                grepl("_MOUSE$", token_up)             ~ "ENTRY",
                is_uniprot_ac(token_base)              ~ "AC",
                TRUE                                   ~ "SYMBOL_OR_ALIAS"
            ),
            mapping_eligible = nzchar(token_up) & grepl("_MOUSE$", token_up)
        )

    # Preserve the historical mouse eligibility rule for the mapping cascade,
    # but keep ineligible rows in the explicit unmapped audit.
    df_tok <- df_tok_all %>% dplyr::filter(mapping_eligible)

    # Extract clean MOUSE reference entities from the local index
    entry_map <- uniprot_mapping %>%
        dplyr::filter(Type == "UniProtKB-ID") %>%
        dplyr::transmute(
            UNIPROT    = toupper(trimws(UniProt_Accession)),
            entry_full = toupper(trimws(Value)),
            entry_base = toupper(gsub("_MOUSE$", "", trimws(Value)))
        ) %>%
        dplyr::filter(grepl("_MOUSE\\s*$", entry_full), nzchar(UNIPROT)) %>%
        dplyr::distinct(entry_base, .keep_all = TRUE)

    # Build primary synonym translation table and resolve review status overlaps (A0A -> SP/TR)
    gene_map <- uniprot_mapping %>%
        dplyr::filter(Type %in% c("Gene_Name", "Gene_Name(synonym)", "Gene_Synonym")) %>%
        dplyr::transmute(
            primaryAccession = toupper(trimws(UniProt_Accession)),
            input            = toupper(trimws(Value))
        ) %>%
        dplyr::filter(nzchar(input), nzchar(primaryAccession)) %>%
        dplyr::mutate(pref = !startsWith(primaryAccession, "A0A")) %>% # Favor non-TREMBL logic
        dplyr::arrange(dplyr::desc(pref), primaryAccession, input) %>%
        dplyr::group_by(input) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>%
        dplyr::select(input, primaryAccession)

    gene_map <- dplyr::bind_rows(
        gene_map,
        entry_map %>% dplyr::transmute(input = entry_base, primaryAccession = UNIPROT)
    ) %>% dplyr::distinct(input, .keep_all = TRUE)

    if (!nrow(entry_map)) stop("entry_map is empty after robust parse")

    # Start tracking mapping states
    resolved <- df_tok %>%
        dplyr::transmute(
            source_row_id, token_raw, token_up, token_base, token_kind,
            acc_guess = toupper(acc_guess),
            entry_guess_up = toupper(entry_guess_up),
            id_class = token_kind,
            Resolved_UNIPROT = NA_character_,
            strategy = NA_character_
        )

    # --- Resolution Cascade / ID Mapping Algorithms ---
    
    # Strategy 1: The token clearly contains a strict UniProt Accession
    idx0 <- which(nz(resolved$acc_guess))
    if (length(idx0)) {
        resolved$Resolved_UNIPROT[idx0] <- resolved$acc_guess[idx0]
        resolved$strategy[idx0] <- "accept_accession_in_token"
    }

    # Strategy 2: Pre-parsed base natively resolves as UniProt Accession
    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    if (length(need)) {
        idx_ac <- need[is_uniprot_ac(resolved$token_base[need])]
        if (length(idx_ac)) {
            resolved$Resolved_UNIPROT[idx_ac] <- resolved$token_base[idx_ac]
            resolved$strategy[idx_ac] <- "accept_accession_base"
        }
    }

    # Strategy 3: Lookup recognized ENTRY names against our index
    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_en_guess <- need[nz(resolved$entry_guess_up[need])]
    if (length(idx_en_guess)) {
        key <- to_base_no_iso_mouse(resolved$entry_guess_up[idx_en_guess])
        hit <- entry_map$UNIPROT[match(key, entry_map$entry_base)]
        ok <- nz(hit)
        if (any(ok)) {
            ii <- idx_en_guess[ok]
            resolved$Resolved_UNIPROT[ii] <- hit[ok]
            resolved$strategy[ii] <- "entry_from_token"
        }
    }

    # Strategy 4: Fallback for explicit ENTRY classes without guess hints
    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_en <- need[resolved$id_class[need] %in% c("ENTRY")]
    if (length(idx_en)) {
        hit <- entry_map$UNIPROT[match(toupper(resolved$token_base[idx_en]), entry_map$entry_base)]
        ok <- nz(hit)
        if (any(ok)) {
            ii <- idx_en[ok]
            resolved$Resolved_UNIPROT[ii] <- hit[ok]
            resolved$strategy[ii] <- "entry_local_mouse"
        }
    }

    # Strategy 5: Translate Symbol/Alias to Primary Accessions locally
    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_sym <- need[resolved$id_class[need] == "SYMBOL_OR_ALIAS"]
    if (length(idx_sym) && nrow(gene_map)) {
        base_need <- toupper(resolved$token_base[idx_sym])
        hit <- gene_map$primaryAccession[match(base_need, gene_map$input)]
        ok <- nz(hit)
        if (any(ok)) {
            ii <- idx_sym[ok]
            resolved$Resolved_UNIPROT[ii] <- hit[ok]
            resolved$strategy[ii] <- "gene_local_mouse"
        }
    }

    # Strategy 6: Re-route unresolved entries via the Gene Name local dictionary
    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_entry_to_gene <- need[resolved$id_class[need] %in% c("ENTRY","ENTRY_GUESS")]
    if (length(idx_entry_to_gene) && nrow(gene_map)) {
        entry_basis <- ifelse(nz(resolved$entry_guess_up[idx_entry_to_gene]),
                                                    to_base_no_iso_mouse(resolved$entry_guess_up[idx_entry_to_gene]),
                                                    resolved$token_base[idx_entry_to_gene])
        gene_keys <- toupper(entry_basis)
        hit <- gene_map$primaryAccession[match(gene_keys, gene_map$input)]
        ok <- nz(hit)
        if (any(ok)) {
            ii <- idx_entry_to_gene[ok]
            resolved$Resolved_UNIPROT[ii] <- hit[ok]
            resolved$strategy[ii] <- "entry_prefix_as_gene_local"
        }
    }

    # Strategy 7: Deduce missing proteins using Entry family prefix inference logic
    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_family <- need[resolved$id_class[need] %in% c("SYMBOL_OR_ALIAS","ENTRY","ENTRY_GUESS")]
    if (length(idx_family) && nrow(entry_map)) {
        bases <- unique(toupper(resolved$token_base[idx_family]))
        pick_idx <- lapply(bases, function(b) {
            ix <- which(grepl(paste0("^", b, "[A-Z0-9]+$"), entry_map$entry_base))
            if (!length(ix)) return(NA_integer_)
            em <- entry_map[ix, , drop = FALSE]
            em$suffix <- sub(paste0("^", b), "", em$entry_base)
            em$score_review <- ifelse(startsWith(em$UNIPROT, "A0A"), 1L, 0L)
            em$score_suffixA <- ifelse(grepl("^A", em$suffix), 0L, 1L)
            ord <- order(em$score_review, em$score_suffixA, em$suffix, em$UNIPROT)
            ix[ord[1]]
        })
        sel <- !is.na(pick_idx)
        if (any(sel)) {
            b_ok <- bases[sel]
            acc_ok <- entry_map$UNIPROT[unlist(pick_idx[sel])]
            map_vec <- stats::setNames(acc_ok, b_ok)
            key_now <- toupper(resolved$token_base[idx_family])
            hit <- unname(map_vec[key_now])
            ok <- nz(hit)
            ii <- idx_family[ok]
            if (length(ii)) {
                resolved$Resolved_UNIPROT[ii] <- hit[ok]
                resolved$strategy[ii] <- "entry_family_prefix_pick"
            }
        }
    }

    # Strategy 8: External DB Integration - Bioconductor OrgDb API (org.Mm.eg.db)
    if (requireNamespace("AnnotationDbi", quietly = TRUE) && requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
        need_idx <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
        need_ids <- toupper(unique(resolved$token_base[need_idx]))
        ids_ent <- unique(need_ids[!is_uniprot_ac(need_ids)])
        if (length(ids_ent)) {
            sel_sym <- try(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = ids_ent, keytype = "SYMBOL", columns = c("MGIID","ENTREZID","UNIPROT","SYMBOL")), silent = TRUE)
            map_sym <- tibble::tibble(input = character(), primaryAccession = character())
            if (!inherits(sel_sym, "try-error") && nrow(sel_sym)) {
                map_sym <- tibble::as_tibble(sel_sym) %>%
                    dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
                    dplyr::group_by(SYMBOL) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>%
                    dplyr::transmute(input = toupper(SYMBOL), primaryAccession = toupper(UNIPROT))
            }
            kt <- try(AnnotationDbi::keytypes(org.Mm.eg.db::org.Mm.eg.db), silent = TRUE)
            map_alias <- tibble::tibble(input = character(), primaryAccession = character())
            if (!inherits(kt, "try-error") && "ALIAS" %in% kt) {
                sel_alias <- try(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = ids_ent, keytype = "ALIAS", columns = c("UNIPROT","ALIAS")), silent = TRUE)
                if (!inherits(sel_alias, "try-error") && nrow(sel_alias)) {
                    map_alias <- tibble::as_tibble(sel_alias) %>%
                        dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
                        dplyr::group_by(ALIAS) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>%
                        dplyr::transmute(input = toupper(ALIAS), primaryAccession = toupper(UNIPROT))
                }
            }
            map_symall <- dplyr::bind_rows(map_sym, map_alias) %>% dplyr::distinct(input, .keep_all = TRUE)
            if (nrow(map_symall)) {
                need_idx_now <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
                base_need <- toupper(resolved$token_base[need_idx_now])
                hit <- map_symall$primaryAccession[match(base_need, map_symall$input)]
                ok <- !is.na(hit) & nzchar(hit)
                ii <- need_idx_now[ok]
                if (length(ii)) {
                    resolved$Resolved_UNIPROT[ii] <- hit[ok]
                    resolved$strategy[ii] <- "orgdb_mgi_symbol_first"
                }
            }
        }

        # Sub-strategy: Route Symbol -> EntrezID -> UniProt to bridge mapping gaps
        need_idx <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
        need_ids <- toupper(unique(resolved$token_base[need_idx]))
        sym_left <- unique(need_ids[grepl("^[A-Z0-9\\-]{2,}$", need_ids)])
        if (length(sym_left)) {
            sym2eg <- try(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = sym_left, keytype = "SYMBOL", columns = c("ENTREZID","SYMBOL")), silent = TRUE)
            eg2up  <- tibble::tibble()
            if (!inherits(sym2eg, "try-error") && nrow(sym2eg)) {
                ekeys <- unique(na.omit(sym2eg$ENTREZID))
                if (length(ekeys)) {
                    egsel <- try(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = ekeys, keytype = "ENTREZID", columns = c("UNIPROT","ENTREZID")), silent = TRUE)
                    if (!inherits(egsel, "try-error") && nrow(egsel)) {
                        eg2up <- tibble::as_tibble(egsel) %>%
                            dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
                            dplyr::group_by(ENTREZID) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup()
                    }
                }
                if (nrow(eg2up)) {
                    map_sym2up <- tibble::as_tibble(sym2eg) %>%
                        dplyr::distinct(SYMBOL, ENTREZID) %>%
                        dplyr::left_join(eg2up, by = "ENTREZID") %>%
                        dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
                        dplyr::transmute(input = toupper(SYMBOL), primaryAccession = toupper(UNIPROT)) %>%
                        dplyr::distinct(input, .keep_all = TRUE)
                    need_idx2 <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
                    base_need <- toupper(resolved$token_base[need_idx2])
                    hit <- map_sym2up$primaryAccession[match(base_need, map_sym2up$input)]
                    ok <- !is.na(hit) & nzchar(hit)
                    ii <- need_idx2[ok]
                    if (length(ii)) {
                        resolved$Resolved_UNIPROT[ii] <- hit[ok]
                        resolved$strategy[ii] <- "orgdb_symbol_entrez_uniprot"
                    }
                }
            }
        }
    }

    # Strategy 9: External DB Integration - Online UniProt REST API Fetching
    need_idx <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
    need_ids <- unique(toupper(resolved$token_base[need_idx]))
    sym_left2 <- unique(need_ids[grepl("^[A-Z0-9\\-]{2,}$", need_ids)])
    if (length(sym_left2) && requireNamespace("UniProt.ws", quietly = TRUE)) {
        batch_vec <- split(sym_left2, ceiling(seq_along(sym_left2) / 50))
        picks <- list()
        for (b in batch_vec) {
            q_list <- lapply(b, function(g) list(organism_id = 10090, gene_primary = g))
            query_once <- function(ql) try(UniProt.ws::queryUniProt(query = ql, fields = c("accession","id","gene_primary","reviewed"), collapse = "OR", n = 10, pageSize = 10), silent = TRUE)
            res_list <- lapply(q_list, function(ql) { out <- query_once(ql); if (inherits(out, "try-error") || !is.data.frame(out)) { Sys.sleep(0.8); out <- query_once(ql) }; out })
            ok <- res_list[!vapply(res_list, inherits, logical(1), "try-error")]
            if (length(ok)) {
                tbl <- dplyr::bind_rows(lapply(ok, tibble::as_tibble))
                if (nrow(tbl)) {
                    tbl <- tbl %>% dplyr::mutate(gene_primary = toupper(.data$gene_primary), accession = toupper(.data$accession))
                    pick <- tbl %>% dplyr::group_by(gene_primary) %>% dplyr::arrange(dplyr::desc(.data$reviewed), accession, .by_group = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>% dplyr::transmute(input = gene_primary, primaryAccession = accession)
                    picks[[length(picks) + 1]] <- pick
                }
            }
        }
        if (length(picks)) {
            map_gene <- dplyr::bind_rows(picks) %>% dplyr::distinct(input, .keep_all = TRUE)
            need_idx3 <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
            base_need <- toupper(resolved$token_base[need_idx3])
            hit <- map_gene$primaryAccession[match(base_need, map_gene$input)]
            ok <- !is.na(hit) & nzchar(hit)
            ii <- need_idx3[ok]
            if (length(ii)) {
                resolved$Resolved_UNIPROT[ii] <- hit[ok]
                resolved$strategy[ii] <- "uniprot_gene_primary_retry"
            }
        }
    }

    # Strategy 10: Late-stage entry resolution explicitly pinging UniProt.ws
    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_entry_ws <- need[resolved$id_class[need] %in% c("ENTRY","ENTRY_GUESS")]
    if (length(idx_entry_ws) && requireNamespace("UniProt.ws", quietly = TRUE)) {
        left_ids <- unique(toupper(ifelse(nz(resolved$entry_guess_up[idx_entry_ws]),
                                                                            resolved$entry_guess_up[idx_entry_ws],
                                                                            paste0(resolved$token_base[idx_entry_ws], "_MOUSE"))))
        left_ids <- gsub("-\\d+$", "", left_ids)
        batch_vec <- split(left_ids, ceiling(seq_along(left_ids) / 50))
        picks <- list()
        for (b in batch_vec) {
            q_list <- lapply(b, function(ua) list(organism_id = 10090, id = ua))
            query_once <- function(ql) try(UniProt.ws::queryUniProt(query = ql, fields = c("accession","id","reviewed"), collapse = "OR", n = 5, pageSize = 5), silent = TRUE)
            res_list <- lapply(q_list, function(ql) { out <- query_once(ql); if (inherits(out, "try-error") || !is.data.frame(out)) { Sys.sleep(0.8); out <- query_once(ql) }; out })
            ok <- res_list[!vapply(res_list, inherits, logical(1), "try-error")]
            if (length(ok)) {
                tbl <- dplyr::bind_rows(lapply(ok, tibble::as_tibble))
                if (nrow(tbl) && all(c("id", "accession") %in% names(tbl))) {
                    tbl <- tbl %>% dplyr::mutate(id = toupper(.data$id), accession = toupper(.data$accession))
                    pick <- tbl %>% dplyr::group_by(id) %>% dplyr::arrange(dplyr::desc(.data$reviewed), accession, .by_group = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>% dplyr::transmute(input = id, primaryAccession = accession)
                    picks[[length(picks) + 1]] <- pick
                }
            }
        }
        if (length(picks)) {
            map_id <- dplyr::bind_rows(picks) %>% dplyr::distinct(input, .keep_all = TRUE)
            keys_now <- toupper(ifelse(nz(resolved$entry_guess_up[idx_entry_ws]), resolved$entry_guess_up[idx_entry_ws], paste0(resolved$token_base[idx_entry_ws], "_MOUSE")))
            hit <- map_id$primaryAccession[match(keys_now, map_id$input)]
            ok <- nz(hit)
            ii <- idx_entry_ws[ok]
            if (length(ii)) {
                resolved$Resolved_UNIPROT[ii] <- hit[ok]
                resolved$strategy[ii] <- "uniprot_id_retry"
            }
        }
    }

    # --- Manual Override Integration ---
    # Apply curation mappings as a priority lookup
    if (!is.null(manual_mapping)) {
        symbol_cols <- intersect(names(manual_mapping), c("original_symbol","symbol","input","gene_symbol","token_raw"))
        mapped_cols <- intersect(names(manual_mapping), c("final_accession","accession","uniprot","uniprot_accession","mapped_gene_symbol"))
        base_cols   <- intersect(names(manual_mapping), c("base_name","token_base","base","symbol_base"))

        if (length(mapped_cols) && (length(symbol_cols) || length(base_cols))) {
            mm_symbol <- if (length(symbol_cols)) symbol_cols[1] else NA
            mm_mapped <- mapped_cols[1]
            mm_base   <- if (length(base_cols)) base_cols[1] else NA

            mm_clean <- manual_mapping
            if (!is.na(mm_symbol)) mm_clean <- mm_clean %>% dplyr::filter(nzchar(.data[[mm_symbol]]))
            if (!is.na(mm_base))   mm_clean <- mm_clean %>% dplyr::filter(nzchar(.data[[mm_base]]))
            mm_clean <- mm_clean %>% dplyr::filter(nzchar(.data[[mm_mapped]]))

            if (nrow(mm_clean)) {
                # Ensure values represent absolute UniProt Accessions
                map_to_acc <- function(vals) {
                    vals <- toupper(trimws(as.character(vals)))
                    out <- ifelse(is_uniprot_ac(vals), vals, NA_character_)
                    need <- is.na(out) | !nzchar(out)
                    if (any(need)) {
                        base <- to_base_no_iso_mouse(vals[need])
                        hit <- entry_map$UNIPROT[match(base, entry_map$entry_base)]
                        ok <- nz(hit)
                        if (any(ok)) out[need][ok] <- hit[ok]
                    }
                    need <- is.na(out) | !nzchar(out)
                    if (any(need) && nrow(gene_map)) {
                        key <- toupper(to_base_no_iso_mouse(vals[need]))
                        hit <- gene_map$primaryAccession[match(key, gene_map$input)]
                        ok <- nz(hit)
                        if (any(ok)) out[need][ok] <- hit[ok]
                    }
                    out
                }

                mm_clean$.__acc__ <- map_to_acc(mm_clean[[mm_mapped]])
                mm_clean <- mm_clean %>% dplyr::filter(!is.na(.__acc__) & nzchar(.__acc__))

                if (nrow(mm_clean)) {
                    if (!is.na(mm_symbol)) {
                        sym_map <- stats::setNames(mm_clean$.__acc__, mm_clean[[mm_symbol]])
                        tgt <- toupper(resolved$token_raw)
                        hit <- toupper(unname(sym_map[tgt]))
                        ok <- nz(hit)
                        idx <- which(ok)
                        if (length(idx)) {
                            if (isTRUE(manual_override)) {
                                resolved$Resolved_UNIPROT[idx] <- hit[idx]
                                resolved$strategy[idx] <- ifelse(is.na(resolved$strategy[idx]), "manual_symbol", paste0(resolved$strategy[idx], "|manual_symbol"))
                            } else {
                                need_idx <- idx[is.na(resolved$Resolved_UNIPROT[idx]) | !nzchar(resolved$Resolved_UNIPROT[idx])]
                                if (length(need_idx)) {
                                    resolved$Resolved_UNIPROT[need_idx] <- hit[need_idx]
                                    resolved$strategy[need_idx] <- "manual_symbol"
                                }
                            }
                        }
                    }
                    if (!is.na(mm_base)) {
                        base_map <- stats::setNames(mm_clean$.__acc__, mm_clean[[mm_base]])
                        tgt <- toupper(resolved$token_base)
                        hit <- toupper(unname(base_map[tgt]))
                        ok <- nz(hit)
                        idx <- which(ok)
                        if (length(idx)) {
                            if (isTRUE(manual_override)) {
                                resolved$Resolved_UNIPROT[idx] <- hit[idx]
                                resolved$strategy[idx] <- ifelse(is.na(resolved$strategy[idx]), "manual_base", paste0(resolved$strategy[idx], "|manual_base"))
                            } else {
                                need_idx <- idx[is.na(resolved$Resolved_UNIPROT[idx]) | !nzchar(resolved$Resolved_UNIPROT[idx])]
                                if (length(need_idx)) {
                                    resolved$Resolved_UNIPROT[need_idx] <- hit[need_idx]
                                    resolved$strategy[need_idx] <- "manual_base"
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    # --- Data Combination & QC Output ---
    # Join by immutable source row rather than identifier text. This prevents
    # duplicated identifiers or protein groups from multiplying/dropping rows.
    eligible_mapping <- resolved %>%
        dplyr::left_join(
            df_tok_all %>% dplyr::select(
                source_row_id, original_protein_id, token_raw, token_base,
                multi_protein
            ),
            by = c("source_row_id", "token_raw", "token_base")
        ) %>%
        dplyr::transmute(
            source_row_id,
            original_protein_id,
            mapping_identifier = token_raw,
            base_name = token_base,
            uniprot_accession = Resolved_UNIPROT,
            mapping_strategy = strategy,
            multi_protein
        )
    ineligible_mapping <- df_tok_all %>%
        dplyr::filter(!mapping_eligible) %>%
        dplyr::transmute(
            source_row_id,
            original_protein_id,
            mapping_identifier = token_raw,
            base_name = token_base,
            uniprot_accession = NA_character_,
            mapping_strategy = dplyr::if_else(
                !nzchar(token_up),
                "ineligible_empty_identifier",
                "ineligible_non_mouse_identifier"
            ),
            multi_protein
        )
    mapping_info <- dplyr::bind_rows(eligible_mapping, ineligible_mapping) %>%
        dplyr::arrange(source_row_id) %>%
        dplyr::mutate(
            mapped_gene_symbol = unname(accession_gene_symbol_map[toupper(uniprot_accession)]),
            mapping_status = dplyr::if_else(
                !is.na(uniprot_accession) & nzchar(uniprot_accession),
                "mapped",
                "unmapped"
            ),
            original_row_id = source_row_id,
            original_symbol = original_protein_id,
            final_accession = uniprot_accession,
            matched_by = mapping_strategy,
            mapped_status = mapping_status
        )

    output_tables <- build_mapthatprot_output_tables(
        source_input,
        mapping_info %>% dplyr::select(
            source_row_id, original_protein_id, uniprot_accession,
            mapped_gene_symbol, mapping_strategy, mapping_status, multi_protein
        )
    )
    df_mapped <- output_tables$mapped
    unmapped_proteins <- output_tables$unmapped
    validate_mapthatprot_partition(source_input, df_mapped, unmapped_proteins)

    # Log biological multiplicity events (e.g. protein grouping)
    multi_protein_log <- df_raw %>%
    dplyr::filter(multi_protein) %>%
    dplyr::transmute(
        source_file = basename(data_path),
        original_row_id,
        original_entry = original_gene_symbol,
        kept_protein = purrr::map_chr(stringr::str_split(original_gene_symbol, ";"), 1),
        dropped_proteins = purrr::map_chr(
            stringr::str_split(original_gene_symbol, ";"),
            ~ if (length(.) > 1) paste(.[-1], collapse = ";") else NA_character_
        )
    )

    # File IO bindings
    manifest_index <- match(basename(data_path), input_manifest$input_filename)
    if (is.na(manifest_index)) stop("Input file is absent from authoritative index: ", data_path, call. = FALSE)
    manifest_row <- input_manifest[manifest_index, , drop = FALSE]
    base <- tools::file_path_sans_ext(manifest_row$output_filename[[1]])
    mapped_file <- file.path(mapped_dir, manifest_row$output_filename[[1]])
    unmapped_file <- file.path(unmapped_dir, manifest_row$output_filename[[1]])
    info_table_file <- file.path(info_dir, paste0(base, "_mapping_info.csv"))
    info_summary_file <- file.path(info_dir, paste0(base, "_info.txt"))
    mapped_summary_file <- file.path(mapped_summary_dir, paste0(base, "_summary.csv"))
    unmapped_summary_file <- file.path(unmapped_summary_dir, paste0(base, "_summary.csv"))

    readr::write_csv(df_mapped, mapped_file)
    readr::write_csv(unmapped_proteins, unmapped_file)
    readr::write_csv(mapping_info, info_table_file)

    # Extract mapped summaries showing conversion strategy metadata
    mapped_summary <- mapping_info %>%
        dplyr::filter(mapping_status == "mapped") %>%
        dplyr::group_by(uniprot_accession) %>%
        dplyr::summarise(
            original_protein_ids = paste(unique(original_protein_id), collapse = "; "),
            mapped_gene_symbols = paste(unique(stats::na.omit(mapped_gene_symbol)), collapse = "; "),
            strategies = paste(unique(stats::na.omit(mapping_strategy)), collapse = "; "),
            .groups = "drop"
        )
    readr::write_csv(mapped_summary, mapped_summary_file)

    # Extract unmapped summary metrics
    unmapped_summary <- unmapped_proteins %>%
        dplyr::group_by(original_protein_id) %>%
        dplyr::summarise(occurrences = dplyr::n(), .groups = "drop")
    readr::write_csv(unmapped_summary, unmapped_summary_file)

    # Compute high-level file summary statistics
    total_in <- nrow(source_input)
    total_mapped <- nrow(df_mapped)
    total_unmapped <- total_in - total_mapped
    total_multi_protein <- sum(mapping_info$multi_protein)
    mapped_accessions <- mapping_info$uniprot_accession[mapping_info$mapping_status == "mapped"]
    duplicate_canonical_mappings <- sum(duplicated(mapped_accessions))
    strategy_counts <- mapping_info %>%
        dplyr::filter(!is.na(mapping_strategy) & nzchar(mapping_strategy)) %>%
        dplyr::count(mapping_strategy, name = "n") %>%
        dplyr::arrange(dplyr::desc(n))
    strategy_summary <- if (nrow(strategy_counts)) {
        paste0(strategy_counts$mapping_strategy, "=", strategy_counts$n, collapse = ";")
    } else {
        "none"
    }

    summary_lines <- c(
        paste0("file: ", basename(data_path)),
        paste0("total_input_rows: ", total_in),
        paste0("mapped: ", total_mapped),
        paste0("unmapped: ", total_unmapped),
        paste0("multi_protein_rows: ", total_multi_protein),
        paste0("duplicate_canonical_mappings: ", duplicate_canonical_mappings),
        paste0("row_accounting: ", total_mapped + total_unmapped, "/", total_in),
        "strategy_counts:"
    )
    if (nrow(strategy_counts) > 0) {
        strategy_str <- paste0("  - ", strategy_counts$mapping_strategy, ": ", strategy_counts$n)
        summary_lines <- c(summary_lines, strategy_str)
    } else {
        summary_lines <- c(summary_lines, "  - none")
    }
    writeLines(summary_lines, info_summary_file)

    file_index <- data.frame(
        canonical_comparison = manifest_row$canonical_comparison[[1]],
        canonical_contrast = manifest_row$canonical_contrast[[1]],
        sample_class = manifest_row$sample_class[[1]],
        numerator_condition = manifest_row$numerator_condition[[1]],
        denominator_condition = manifest_row$denominator_condition[[1]],
        historical_comparison_alias = manifest_row$historical_comparison_alias[[1]],
        mapping_direction = map_direction,
        input_path = normalize_mapthatprot_path(data_path, must_work = TRUE),
        output_mapped_path = normalize_mapthatprot_path(mapped_file, must_work = TRUE),
        output_unmapped_path = normalize_mapthatprot_path(unmapped_file, must_work = TRUE),
        n_input_proteins = total_in,
        n_output_mapped_rows = nrow(df_mapped),
        n_accounted_output_rows = nrow(df_mapped) + nrow(unmapped_proteins),
        n_mapped = total_mapped,
        n_unmapped = total_unmapped,
        mapping_rate = total_mapped / total_in,
        mapping_strategy_summary = strategy_summary,
        n_multi_protein_rows = total_multi_protein,
        multi_protein_handling = "leading master protein mapped; additional group members retained in audit",
        n_duplicate_canonical_mappings = duplicate_canonical_mappings,
        duplicate_canonical_mapping_handling = "retained and counted; no source rows collapsed",
        row_accounting_valid = total_mapped + total_unmapped == total_in,
        source_split_index_path = split_index_provenance$source_split_index_path[[1]],
        source_split_index_sha256 = split_index_provenance$source_split_index_sha256[[1]],
        source_split_sha256 = manifest_row$source_split_sha256[[1]],
        source_gct_sha256 = manifest_row$source_gct_sha256[[1]],
        mapping_reference_file = mapping_reference_provenance$mapping_reference_file[[1]],
        mapping_reference_path = mapping_reference_provenance$mapping_reference_path[[1]],
        mapping_reference_version = mapping_reference_provenance$mapping_reference_version[[1]],
        mapping_reference_snapshot_date_utc = mapping_reference_provenance$mapping_reference_snapshot_date_utc[[1]],
        mapping_reference_sha256 = mapping_reference_provenance$mapping_reference_sha256[[1]],
        mapping_reference_bytes = mapping_reference_provenance$mapping_reference_bytes[[1]],
        mapping_reference_modified_utc = mapping_reference_provenance$mapping_reference_modified_utc[[1]],
        manual_mapping_path = manual_mapping_provenance$manual_mapping_path[[1]],
        manual_mapping_sha256 = manual_mapping_provenance$manual_mapping_sha256[[1]],
        manual_mapping_rows = manual_mapping_provenance$manual_mapping_rows[[1]],
        mapped_output_sha256 = mapthatprot_sha256(mapped_file),
        run_timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        stringsAsFactors = FALSE
    )

    invisible(list(
        mapping_table = mapping_info %>%
            dplyr::mutate(source_file = basename(data_path)),
        unmapped_table = unmapped_proteins %>%
            dplyr::mutate(source_file = basename(data_path)),
        multi_protein_log_table = multi_protein_log,
        file_index = file_index
    ))

}

# -------------------------
# Parallel Execution Engine
# -------------------------
# Distribute file-level mapping across available processor cores
n_files <- length(csv_files)
available_cores <- parallel::detectCores(logical = FALSE)
workers <- max(1, min(available_cores - 1, n_files))

cat("Setting up parallel processing backend using", workers, "out of", available_cores, "available cores...\n")
cl <- parallel::makeCluster(workers)
doParallel::registerDoParallel(cl)

cat("Initiating parallel mapping cascade for all files...\n")
results <- foreach(i = seq_along(csv_files),
                   .packages = c("dplyr", "stringr", "tidyr", "purrr", "readr", "R.utils")) %dopar% {
    process_file(csv_files[i])
}

parallel::stopCluster(cl)
cat("Batch parallel mapping successfully completed for", n_files, "files.\n")

# -------------------------
# Global Mapping Workbooks
# -------------------------

cat("Aggregating overall biology summaries and computing mapping strategy statistics...\n")

# Aggregate outputs back from processing clusters
all_mapping_tables <- purrr::map(results, "mapping_table") %>% dplyr::bind_rows()
all_unmapped_tables <- purrr::map(results, "unmapped_table") %>% dplyr::bind_rows()
all_dropped_proteins <- purrr::map(results, "multi_protein_log_table") %>% dplyr::bind_rows()
mapped_comparison_index <- purrr::map(results, "file_index") %>% dplyr::bind_rows()
if (
    nrow(mapped_comparison_index) != 12L ||
    any(!mapped_comparison_index$row_accounting_valid) ||
    anyDuplicated(mapped_comparison_index$canonical_comparison)
) {
    stop("Mapped comparison index failed the canonical 12-comparison row-accounting contract.", call. = FALSE)
}
mapped_comparison_index <- mapped_comparison_index[
    match(input_manifest$canonical_comparison, mapped_comparison_index$canonical_comparison),
    ,
    drop = FALSE
]
readr::write_csv(mapped_comparison_index, mapped_index_path)
cat("Saved authoritative mapped comparison index to:", mapped_index_path, "\n")


# Consolidate standard mappings
global_mapping_summary <- all_mapping_tables %>%
    dplyr::distinct() %>%
    dplyr::arrange(source_file, original_symbol)

global_mapping_file <- file.path(info_dir, "GLOBAL_mapping_summary.csv")
readr::write_csv(global_mapping_summary, global_mapping_file)

# Backup in summaries dir
global_mapping_summary_file2 <- file.path(mapped_summary_dir, "GLOBAL_mapping_summary.csv")
readr::write_csv(global_mapping_summary, global_mapping_summary_file2)

cat("Saved global mapping summary to:", info_dir, "and", mapped_summary_dir, "\n")

# Consolidate ID unmapped dropouts
global_unmapped_summary <- all_unmapped_tables %>%
    dplyr::group_by(original_protein_id) %>%
    dplyr::summarise(
        occurrences = dplyr::n(),
        files_present = paste(unique(source_file), collapse = "; "),
        .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(occurrences))

global_unmapped_file <- file.path(unmapped_dir, "GLOBAL_unmapped_proteins.csv")
readr::write_csv(global_unmapped_summary, global_unmapped_file)

global_unmapped_summary_file2 <- file.path(unmapped_summary_dir, "GLOBAL_unmapped_proteins.csv")
readr::write_csv(global_unmapped_summary, global_unmapped_summary_file2)

cat("Saved global unmapped tracking to:", unmapped_dir, "and", unmapped_summary_dir, "\n")

# -------------------------
# Strategy Quality Control
# -------------------------

strategy_stats <- global_mapping_summary %>%
    dplyr::filter(mapped_status == "mapped") %>%
    dplyr::count(matched_by, name = "total_mapped") %>%
    dplyr::arrange(dplyr::desc(total_mapped))

strategy_file <- file.path(info_dir, "GLOBAL_strategy_statistics.csv")
readr::write_csv(strategy_stats, strategy_file)
cat("Saved global strategy statistics to:", strategy_file, "\n")

mapping_full <- all_mapping_tables %>%
    dplyr::distinct() %>%
    dplyr::arrange(source_file, original_symbol)

# Produce unique biological entity tracking
protein_summary <- mapping_full %>%
    dplyr::group_by(original_symbol) %>%
    dplyr::summarise(
        mapped_to = paste(unique(na.omit(final_accession)), collapse = "; "),
        strategies = paste(unique(na.omit(matched_by)), collapse = "; "),
        files_present = paste(unique(source_file), collapse = "; "),
        mapped = any(!is.na(final_accession)),
        .groups = "drop"
    )

unmapped_proteins_global <- protein_summary %>%
    dplyr::filter(mapped == FALSE)


# Tally relative efficiency across mapping resolution methods
strategy_stats <- mapping_full %>%
    dplyr::filter(!is.na(final_accession)) %>%
    dplyr::count(matched_by, name = "count") %>%
    dplyr::arrange(desc(count))

strategy_stats_unique <- mapping_full %>%
    dplyr::filter(!is.na(final_accession)) %>%
    dplyr::distinct(matched_by, original_symbol) %>%
    dplyr::count(matched_by, name = "unique_proteins") %>%
    dplyr::arrange(desc(unique_proteins))

# Calculate dataset coverage % 
coverage_stats <- mapping_full %>%
    dplyr::group_by(source_file) %>%
    dplyr::summarise(
        total = dplyr::n(),
        mapped = sum(!is.na(final_accession)),
        coverage = mapped / total,
        .groups = "drop"
    )

# --- Generate Comprehensive Excel QC Document ---
cat("Generating comprehensive Excel QC Report...\n")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

report_file <- file.path(report_dir, "Mapping_QC_Report.xlsx")

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Full_Mapping_Table")
openxlsx::writeData(wb, "Full_Mapping_Table", mapping_full)

openxlsx::addWorksheet(wb, "Protein_Summary")
openxlsx::writeData(wb, "Protein_Summary", protein_summary)

openxlsx::addWorksheet(wb, "Unmapped_Proteins")
openxlsx::writeData(wb, "Unmapped_Proteins", unmapped_proteins_global)

openxlsx::addWorksheet(wb, "Dropped_Proteins")
openxlsx::writeData(wb, "Dropped_Proteins", all_dropped_proteins)

openxlsx::addWorksheet(wb, "Strategy_Stats")
openxlsx::writeData(wb, "Strategy_Stats", strategy_stats)

openxlsx::addWorksheet(wb, "Unique_Strategy_Stats")
openxlsx::writeData(wb, "Unique_Strategy_Stats", strategy_stats_unique)

openxlsx::addWorksheet(wb, "Coverage_Stats")
openxlsx::writeData(wb, "Coverage_Stats", coverage_stats)

openxlsx::addWorksheet(wb, "Mapped_Comparisons")
openxlsx::writeData(wb, "Mapped_Comparisons", mapped_comparison_index)

openxlsx::saveWorkbook(wb, report_file, overwrite = TRUE)

cat("Saved master QC workbook to:", report_file, "\n")
cat("====================================================\n")
cat("MapThatProt_batch execution completed successfully!\n")
cat("====================================================\n")

