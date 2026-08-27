# Neha Proteomics Analysis Workflow

## Overview

This folder contains the Neha proteomics workflow for preprocessing, quality control, differential analysis, and EWCE enrichment against external reference cell types.

Shared label definitions live in `R/analysis_labels.R`.

Sample classes:

- `mcherry`
- `neuropil`
- `cfos`
- `neuron`

> **Naming gotcha:** historical filenames and metadata use **`bg`** (from "background") where
> current outputs say **`neuropil`**. So `bg1_bg2.csv`, `bg2_bg4.csv` and `group2` values like
> `bg_1` all refer to the **neuropil tissue compartment** — *not* a negative control, blank, or
> background subtraction. The alias table lives in `R/analysis_labels.R`
> (`sample_class_aliases`), but it is invisible when reading raw filenames.

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

Run the animal-level EWCE sampling-unit checks:

```bash
Rscript tests/test_ewce_animal_level.R
```

Run the downstream animal-level PCA/QC checks:

```bash
Rscript tests/test_pca_animal_level.R
```

### Pipeline smoke test

`run_pipeline_check.ps1` runs the contract tests and then the runnable pipeline stages,
redirecting every producing stage to a scratch directory via its documented env-var override,
so validated outputs are never touched:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1            # Normal (default)
powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1 -Tier Fast # tests only, ~2 min
powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1 -Tier Full # + mapping/enrichment/EWCE, hours
```

It prints PASS/FAIL/SKIP per stage and exits non-zero on any failure. Stages whose inputs are
genuinely absent from the project tree are reported as SKIP with the reason, rather than
failing. `ProTigy` sits in the middle of the chain as an **external tool**, so no script can run
the pipeline truly end to end; the runner covers both sides of that gap using the committed
ProTigy output. The strongest single check it performs is rebuilding the animal-level GCT from
its inputs and confirming it is bit-identical to the locked
`f12cf99e1bfb7c17bbf56bffb6783e924698bce5d5533a8e312bc4bbb733bbb3`.

Run the data-integrity contracts on their own:

```bash
Rscript tests/test_data_integrity.R
```

Run the design-balance / identifiability contracts:

```bash
Rscript tests/test_design_balance.R
```

This asserts the documented confound structure of the cohort — most importantly that
collection/acquisition plate is near-perfectly aliased with Pairing (11/12 animals), and that
the `paired_veh` vs `unpaired_veh` contrast is **100%** aliased with plate. It quantifies which
simple contrasts are plate-protected and which are not, so the confound cannot be rediscovered
by accident. The cohort-level checks skip cleanly when the shared drive is unreachable.

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

Large input files and some generated results are kept outside Git. Active QC scripts remain under `03_qc_exploration/`.

### Data layout on shared storage

Restructured **2026-08-26** under `S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/`:

| Folder | Contents | Pipeline writes here? |
|---|---|---|
| `01_input/` | external / hand-placed: `references/`, `metadata/`, `raw_proteomics/`, `single_cell/` | no |
| `02_data/` | derived intermediates: `gct/`, `animal_level/{input_gct,split,mapped,audits,qc_exports}` | yes |
| `03_output/` | results: `enrichment/`, `ewce/`, `pca/`, `synthesis/`, `inferential_checks/`, `qc/` | yes |
| `99_historical/` | superseded hemisphere-level era, prior generations, archives | no (read-only) |

**See [`CANONICAL_OUTPUTS.md`](CANONICAL_OUTPUTS.md) for which specific folder is authoritative
per stage** (several stages have more than one generation on disk), plus the analysis caveats
carried forward from the validation audits. The migration was folder renames only; a reversible
manifest and the original index files are in `_migration_backup_20260826/`.

The EWCE workflow in `05_celltype_enrichment_EWCE/01_EWCE.r` uses shared sample-class and condition definitions from `R/analysis_labels.R`.

## Animal-Level ProTigy Handoff

`01_preprocessing/02a_prepare_animal_level_protigy_input.r` creates a separate ProTigy GCT in which valid Left and Right hemisphere values are averaged within `AnimalID × sample_class`. It consumes the existing 5,349-protein processed/imputed matrix and does not transform, normalize, filter, impute, or remap it. The stage writes the corrected matrix, source-sample and feature-identity audits, design summaries, a historically scoped within-class contrast manifest, and SHA-256 provenance under a separate `protigy_input_animal_level` directory. Existing ProTigy inputs and results are not changed.

`01_preprocessing/03_gct_extractR.r` splits the resulting animal-level ProTigy statistical-results GCT into 12 validated forward and 12 derived reverse tables under `02_data/animal_level/split`. It reads statistical fields from their physical GCT positions rather than treating the file as an expression GCT. In split CSVs, the compatibility column `gene_symbol` retains the original GCT `id` (a UniProt-style protein identifier) because the current MapThatProt stage requires that historical column name; it is not the GCT `Description`. `Description` is preserved separately. Override the defaults with `NEHA_PROTIGY_STAT_GCT_INPUT` and `NEHA_PROTIGY_STAT_GCT_OUTPUT_ROOT`.

`02_id_mapping/01_MapThatProt_batch.r` consumes only the 12 files listed in that split branch's `indexComparisons.csv`; it does not scan the historical `99_historical/datasets_raw` tree. Forward mapping is the default, with `NEHA_MAPTHATPROT_DIRECTION=reverse` available only as an explicit override. Outputs are isolated under `02_data/animal_level/mapped`, including 12 canonical mapped CSVs in `forward/`, per-row mapping/unmapped audits, an `indexMappedComparisons.csv`, reference and source SHA-256 provenance, and the existing mapping-strategy QC report. Override roots with `NEHA_MAPTHATPROT_SPLIT_ROOT` and `NEHA_MAPTHATPROT_OUTPUT_ROOT`. `NEHA_MAPTHATPROT_REFERENCE_FILE` and `NEHA_MAPTHATPROT_MANUAL_MAPPING_FILE` override the mapping references; `NEHA_MAPTHATPROT_REFERENCE_VERSION` can record a known UniProt release label because the local three-column idmapping file does not encode one internally.

In mapped comparison CSVs, the first column `gene_symbol` is retained for current downstream compatibility but remains historically misnamed: it contains the resolved UniProt accession used by clusterProfiler. `original_protein_id` preserves the original GCT `id`, `Description` preserves the source gene-style description, and `mapped_gene_symbol` provides the mapped gene annotation when available. All source DA statistics are carried through unchanged. Unmapped rows are written separately and participate in an exact mapped-plus-unmapped row-accounting audit.

`04_differential_expression_enrichment/01_clusterProfiler.r` consumes the 12 forward files recorded in `02_data/animal_level/mapped/indexMappedComparisons.csv` and writes an isolated canonical branch under `03_output/enrichment`. It uses `uniprot_accession` explicitly. Repeated accessions are collapsed once, deterministically, by the largest absolute `log2fc` while preserving its sign, with source-row and protein-ID tie breakers and a per-comparison audit. The same selected source row is then used for every rank sensitivity. Canonical GO and KEGG GSEA rank by the ProTigy/limma moderated `t` statistic by default; positive values still mean higher abundance in the indexed numerator. The standard `t` run also writes `log2fc` sensitivity results (`GSEA_GO_BP_log2fc_sensitivity.csv` and `GSEA_KEGG_log2fc_sensitivity.csv`) with the established comparison/analysis seeds. `NEHA_ENRICHMENT_GSEA_RANK=log2fc` provides an explicit legacy-primary run without duplicating that sensitivity output. GO ORA is unchanged: it uses `padj` and signed `log2fc` foreground definitions and all unique successfully mapped, measured UniProt accessions in that comparison as its explicit universe. Rank source, analysis role, tie diagnostics, deterministic seeds, parameters, input hashes, counts, warnings, and output paths are recorded in `indexEnrichmentComparisons.csv` and per-comparison audits.

The downstream `02_compareGO.r`, `03_compare_pathways.r`, and `04_compare_sig_expr.r` scripts consume the canonical indexes instead of historical filenames. Historical comparison aliases remain provenance metadata only. `04_compare_sig_expr.r` now compares animal-level contrast statistics; it does not fabricate per-animal expression from mapped differential-result tables.

Canonical enrichment defaults can be overridden with `NEHA_ENRICHMENT_MAPPED_ROOT`, `NEHA_ENRICHMENT_MAPPED_INDEX`, `NEHA_ENRICHMENT_OUTPUT_ROOT`, and `NEHA_ENRICHMENT_GSEA_RANK`. Existing isolated outputs are protected unless `NEHA_ENRICHMENT_FORCE=true` is set explicitly. Historical roots under `99_historical/` (`datasets_raw`, `datasets_mapped`, `core_enrichment`, `plots_pairwise`) are rejected as canonical inputs or outputs.

The 2024 Neha instrument IDs do not contain `_L_`/`_R_`. For that historical dataset only, the stage requires the original annotation's explicit `_left`/`_right` label and independently verifies its one-to-one agreement with both archived `ReplicateGroup` fields and the instrument sample identity. Reusable aggregation helpers otherwise accept only `Left`/`L` and `Right`/`R` and require the sample-ID hemisphere token.

The EWCE differential branch consumes `protigy_input_animal_level/neha_protigy_input_animal_level_primary.gct`, the validated handoff containing one equal-weight Left/Right mean per `AnimalID × sample_class`. It fits four 12-animal limma models and the 12 forward comparisons from `neha_primary_contrast_manifest()`, requiring exactly three animals per condition before fitting. It does not repeat aggregation or add normalization, filtering, or imputation to the animal-level abundances; the existing EWCE gene-symbol annotation step is otherwise unchanged. Differential sampling-unit provenance and the existing EWCE parameters are written to `03_QC_Mapping_Logs/animal_level_differential_audit.csv`. Override the input and isolated output root with `NEHA_EWCE_ANIMAL_LEVEL_INPUT` and `NEHA_EWCE_OUTPUT_ROOT`; the historical `99_historical/ewce_legacy` root is rejected.

`03_qc_exploration/06_pcaPlot_Neha.r` remains the PCA entry point but is now a thin orchestrator: the former 2,576-line monolith was split (2026-08-26) into ordered parts under `03_qc_exploration/pca/` (`06a_pca_core.r` plus seven extension parts). The parts run at top level and share the globals created by the core, so behaviour is unchanged — a split-vs-validated run produces an identical 140-file output inventory. Part order matters: `06b` creates the `rot`/`cors_df`/`um`/`npc` objects that `06c`, `06g` and `06h` consume. Each part is now guarded separately, so one failing visualization no longer suppresses every later one (previously all extensions shared a single `tryCatch`); failure still writes the failure audit and exits non-zero.

The downstream biological PCA in `03_qc_exploration/06_pcaPlot_Neha.r` reads that same validated primary animal-level GCT with `validate_protigy_gct_v13()`. Each PCA point is one `AnimalID × sample_class` unit; all 48 units and three animals per sample-class/condition are required. The processed/imputed animal-level abundances are not normalized, filtered for missingness, or imputed again. The primary `prcomp()` remains centered and scaled, with only zero-variance protein rows removed because scaled PCA cannot use them; every removal is audited. Animal-level plots include sample class, condition, AnimalID, and phenotype, and a machine-readable audit is written under the isolated `03_output/pca` root. Override paths with `NEHA_PCA_ANIMAL_LEVEL_INPUT` and `NEHA_PCA_OUTPUT_ROOT`; the historical `99_historical/pca_plots_legacy` root is rejected. Acquisition-level protein/peptide QC and rank-abundance scripts remain sample-level because they evaluate technical acquisition and sample-class intensity behavior rather than downstream biological independence.

## License

This repository uses the MIT License. See `../LICENSE`.

## Citation

Citation metadata is provided in `../CITATION.cff`.
