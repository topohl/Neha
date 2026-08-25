# Neha Proteomics Analysis Workflow

## Overview

This folder contains the Neha proteomics workflow for preprocessing, quality control, differential analysis, and EWCE enrichment against external reference cell types.

Shared label definitions live in `R/analysis_labels.R`.

Sample classes:

- `mcherry`
- `neuropil`
- `cfos`
- `neuron`

Condition codes:

- `1` = `paired_cno`
- `2` = `paired_veh`
- `3` = `unpaired_cno`
- `4` = `unpaired_veh`

The default reference condition is `paired_veh`.

## System Requirements

The lightweight demo runs on a standard laptop or CI runner.

Tested environments:

- Windows 11 with R 4.5.1
- Ubuntu GitHub Actions runner with R 4.5.1

No non-standard hardware is required for the demo. The full EWCE workflow can benefit from more memory and CPU time because it loads Bioconductor reference data and runs bootstrap enrichment.

## Installation

```bash
git clone https://github.com/topohl/Neha.git
cd Neha/proteomic_analysis
```

For the lightweight demo and terminology checks, base R is sufficient.

For full analysis scripts, install the packages used by the target script. See `requirements_R.md` for tested R and package versions. Common CRAN packages include:

- `dplyr`
- `tidyr`
- `readxl`
- `writexl`
- `ggplot2`
- `patchwork`
- `stringr`
- `openxlsx`

EWCE analysis additionally uses Bioconductor packages such as:

- `limma`
- `EWCE`
- `ewceData`
- `org.Mm.eg.db`
- `AnnotationDbi`

Expected install time for the lightweight demo is under one minute when R is already installed. Full EWCE package installation can take longer depending on network speed and Bioconductor cache state.

## Demo Instructions

Run from this folder:

```bash
Rscript run_demo.R
```

The demo reads:

- `demo/input/demo_sample_metadata.csv`
- `demo/input/demo_pg_matrix.csv`

It validates shared label parsing, metadata matching, and comparison-name parsing, then writes:

- `demo/output/demo_result_table.csv`
- `demo/output/demo_comparison_table.csv`

Expected outputs are stored under:

- `demo/expected_output/expected_demo_result_table.csv`
- `demo/expected_output/expected_demo_comparison_table.csv`

Expected demo runtime is under one minute.

## Checks

Run the active-code terminology audit:

```bash
Rscript tests/check_stale_labels.R
```

Run the comparison parser checks:

```bash
Rscript tests/test_gct_comparison_parser.R
```

Run the animal-level aggregation and ProTigy GCT contract checks:

```bash
Rscript tests/test_animal_level_proteomics.R
```

Run the indexed animal-level MapThatProt handoff checks:

```bash
Rscript tests/test_mapthatprot_animal_level.R
```

Run the canonical animal-level enrichment contract checks:

```bash
Rscript tests/test_animal_level_enrichment.R
```

## Running On Your Own Data

Use the shared label helper in `R/analysis_labels.R` when adding or modifying scripts. Input metadata should contain `sample_id`; if `sample_class` or `condition_code` are absent, active preprocessing scripts infer them from canonical labels where possible.

`01_preprocessing/02_excel_convert.r` can be run against local demo data by default. For full data, set R options or environment variables:

- `neha.metadata_path` or `NEHA_METADATA_PATH`
- `neha.excel_convert_file` or `NEHA_EXCEL_CONVERT_FILE`
- `neha.excel_convert_folder` or `NEHA_EXCEL_CONVERT_FOLDER`
- `neha.excel_convert_output` or `NEHA_EXCEL_CONVERT_OUTPUT`
- `neha.excel_convert_mode` or `NEHA_EXCEL_CONVERT_MODE`

Full-analysis paths on shared storage are still supported when these settings point to those files.

## Full Analysis Notes

Large input files and some generated results are kept outside Git. Active QC scripts remain under `03_qc_exploration/`, while obsolete or superseded scripts are under `legacy/`.

The EWCE workflow in `05_celltype_enrichment_EWCE/01_EWCE.r` uses shared sample-class and condition definitions from `R/analysis_labels.R`.

## Animal-Level ProTigy Handoff

`01_preprocessing/02a_prepare_animal_level_protigy_input.r` creates a separate ProTigy GCT in which valid Left and Right hemisphere values are averaged within `AnimalID × sample_class`. It consumes the existing 5,349-protein processed/imputed matrix and does not transform, normalize, filter, impute, or remap it. The stage writes the corrected matrix, source-sample and feature-identity audits, design summaries, a historically scoped within-class contrast manifest, and SHA-256 provenance under a separate `protigy_input_animal_level` directory. Existing ProTigy inputs and results are not changed.

`01_preprocessing/03_gct_extractR.r` splits the resulting animal-level ProTigy statistical-results GCT into 12 validated forward and 12 derived reverse tables under `Datasets/data/protigy_animal_level/split`. It reads statistical fields from their physical GCT positions rather than treating the file as an expression GCT. In split CSVs, the compatibility column `gene_symbol` retains the original GCT `id` (a UniProt-style protein identifier) because the current MapThatProt stage requires that historical column name; it is not the GCT `Description`. `Description` is preserved separately. Override the defaults with `NEHA_PROTIGY_STAT_GCT_INPUT` and `NEHA_PROTIGY_STAT_GCT_OUTPUT_ROOT`.

`02_id_mapping/01_MapThatProt_batch.r` consumes only the 12 files listed in that split branch's `indexComparisons.csv`; it does not scan the historical `Datasets/raw` tree. Forward mapping is the default, with `NEHA_MAPTHATPROT_DIRECTION=reverse` available only as an explicit override. Outputs are isolated under `Datasets/data/protigy_animal_level/mapped`, including 12 canonical mapped CSVs in `forward/`, per-row mapping/unmapped audits, an `indexMappedComparisons.csv`, reference and source SHA-256 provenance, and the existing mapping-strategy QC report. Override roots with `NEHA_MAPTHATPROT_SPLIT_ROOT` and `NEHA_MAPTHATPROT_OUTPUT_ROOT`. `NEHA_MAPTHATPROT_REFERENCE_FILE` and `NEHA_MAPTHATPROT_MANUAL_MAPPING_FILE` override the mapping references; `NEHA_MAPTHATPROT_REFERENCE_VERSION` can record a known UniProt release label because the local three-column idmapping file does not encode one internally.

In mapped comparison CSVs, the first column `gene_symbol` is retained for current downstream compatibility but remains historically misnamed: it contains the resolved UniProt accession used by clusterProfiler. `original_protein_id` preserves the original GCT `id`, `Description` preserves the source gene-style description, and `mapped_gene_symbol` provides the mapped gene annotation when available. All source DA statistics are carried through unchanged. Unmapped rows are written separately and participate in an exact mapped-plus-unmapped row-accounting audit.

`04_differential_expression_enrichment/01_clusterProfiler.r` consumes the 12 forward files recorded in `Datasets/data/protigy_animal_level/mapped/indexMappedComparisons.csv` and writes an isolated canonical branch under `Datasets/data/protigy_animal_level/enrichment`. It uses `uniprot_accession` explicitly. Repeated accessions are collapsed once, deterministically, by the largest absolute `log2fc` while preserving its sign, with source-row and protein-ID tie breakers and a per-comparison audit. The same selected source row is then used for every rank sensitivity. Canonical GO and KEGG GSEA rank by the ProTigy/limma moderated `t` statistic by default; positive values still mean higher abundance in the indexed numerator. The standard `t` run also writes `log2fc` sensitivity results (`GSEA_GO_BP_log2fc_sensitivity.csv` and `GSEA_KEGG_log2fc_sensitivity.csv`) with the established comparison/analysis seeds. `NEHA_ENRICHMENT_GSEA_RANK=log2fc` provides an explicit legacy-primary run without duplicating that sensitivity output. GO ORA is unchanged: it uses `padj` and signed `log2fc` foreground definitions and all unique successfully mapped, measured UniProt accessions in that comparison as its explicit universe. Rank source, analysis role, tie diagnostics, deterministic seeds, parameters, input hashes, counts, warnings, and output paths are recorded in `indexEnrichmentComparisons.csv` and per-comparison audits.

The downstream `02_compareGO.r`, `03_compare_pathways.r`, and `04_compare_sig_expr.r` scripts consume the canonical indexes instead of historical filenames. Historical comparison aliases remain provenance metadata only. `04_compare_sig_expr.r` now compares animal-level contrast statistics; it does not fabricate per-animal expression from mapped differential-result tables.

Canonical enrichment defaults can be overridden with `NEHA_ENRICHMENT_MAPPED_ROOT`, `NEHA_ENRICHMENT_MAPPED_INDEX`, `NEHA_ENRICHMENT_OUTPUT_ROOT`, and `NEHA_ENRICHMENT_GSEA_RANK`. Existing isolated outputs are protected unless `NEHA_ENRICHMENT_FORCE=true` is set explicitly. Historical `Datasets/raw`, `Datasets/mapped`, `Datasets/core_enrichment`, `Results`, and `Plots` roots are rejected as canonical inputs or outputs.

The 2024 Neha instrument IDs do not contain `_L_`/`_R_`. For that historical dataset only, the stage requires the original annotation's explicit `_left`/`_right` label and independently verifies its one-to-one agreement with both archived `ReplicateGroup` fields and the instrument sample identity. Reusable aggregation helpers otherwise accept only `Left`/`L` and `Right`/`R` and require the sample-ID hemisphere token.

The current EWCE differential branch remains sample-level: it parses `AnimalID`, but `run_limma_stratum()` passes hemisphere-level columns directly to `limma::lmFit()`. This handoff deliberately does not alter EWCE. A follow-up should call the same animal-level aggregation helper before fitting its limma contrasts.

## License

This repository uses the MIT License. See `../LICENSE`.

## Citation

Citation metadata is provided in `../CITATION.cff`.
