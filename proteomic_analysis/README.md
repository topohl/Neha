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

The 2024 Neha instrument IDs do not contain `_L_`/`_R_`. For that historical dataset only, the stage requires the original annotation's explicit `_left`/`_right` label and independently verifies its one-to-one agreement with both archived `ReplicateGroup` fields and the instrument sample identity. Reusable aggregation helpers otherwise accept only `Left`/`L` and `Right`/`R` and require the sample-ID hemisphere token.

The current EWCE differential branch remains sample-level: it parses `AnimalID`, but `run_limma_stratum()` passes hemisphere-level columns directly to `limma::lmFit()`. This handoff deliberately does not alter EWCE. A follow-up should call the same animal-level aggregation helper before fitting its limma contrasts.

## License

This repository uses the MIT License. See `../LICENSE`.

## Citation

Citation metadata is provided in `../CITATION.cff`.
