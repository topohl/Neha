# Associative Memory Proteomics — Analysis Workflow

## Overview

This folder contains the Associative Memory Proteomics workflow for preprocessing, quality control, differential analysis, and EWCE enrichment against external reference cell types.

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
git clone https://github.com/topohl/AssociativeMemoryProteomics.git
cd AssociativeMemoryProteomics/proteomic_analysis
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

Run the animal-level rank-abundance checks (also guards against the stage reverting to its
former, absent share-bound input):

```bash
Rscript tests/test_rank_abundance_animal_level.R
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

Run the deprecated compatibility-shim contracts:

```bash
Rscript tests/test_deprecated_path_utils_shim.R
```

This pins `R/neha_path_utils.R`, which is a deprecated shim kept only so the frozen, hash-locked
snapshots in `06_manuscript_figure_revision/` (and their runnable twins on the shared drive) keep
resolving after the 2026-08-28 rename to `R/project_path_utils.R`. No active code sources it, and
the test asserts that.

Run the design-balance / identifiability contracts:

```bash
Rscript tests/test_design_balance.R
```

This asserts the documented association structure of the cohort — most importantly that
collection plate is near-perfectly associated with Pairing (11/12 animals), and that the
`paired_veh` vs `unpaired_veh` contrast is **100%** associated with collection plate.

Because pairing condition was completely associated with collection plate in this comparison,
any collection-plate-associated contribution cannot be distinguished from a pairing-associated
contribution. No downstream proteomics batch effect attributable to collection plate has been
demonstrated — see [`CANONICAL_OUTPUTS.md`](CANONICAL_OUTPUTS.md) for what the `Plate1`/`Plate2`
token does and does not record. The checks also quantify which simple contrasts carry no
collection-plate variation at all, so this design property cannot be rediscovered by accident.
The cohort-level checks skip cleanly when the shared drive is unreachable.

## Running On Your Own Data

Use the shared label helper in `R/analysis_labels.R` when adding or modifying scripts. Input metadata should contain `sample_id`; if `sample_class` or `condition_code` are absent, active preprocessing scripts infer them from canonical labels where possible.

`01_preprocessing/02_excel_convert.r` can be run against local demo data by default. For full data, set R options or environment variables:

- `proteomics.metadata_path` or `PROTEOMICS_METADATA_PATH`
- `proteomics.excel_convert_file` or `PROTEOMICS_EXCEL_CONVERT_FILE`
- `proteomics.excel_convert_folder` or `PROTEOMICS_EXCEL_CONVERT_FOLDER`
- `proteomics.excel_convert_output` or `PROTEOMICS_EXCEL_CONVERT_OUTPUT`
- `proteomics.excel_convert_mode` or `PROTEOMICS_EXCEL_CONVERT_MODE`

Full-analysis paths on shared storage are still supported when these settings point to those files.

## Full Analysis Notes

Large input files and some generated results are kept outside Git. Active QC scripts remain under `03_qc_exploration/`.

### Rank-abundance QC (animal level)

`03_qc_exploration/02_rank_abundance_by_sample_class.r` reads the validated animal-level GCT
(`02_data/animal_level/input_gct/neha_protigy_input_animal_level_primary.gct`) and summarises by
**sample_class × condition** — 16 groups, n = 3 animals each. Hemisphere-level observations are
never used for a biological group summary.

It regenerates the rank-abundance panels used in the corrected proteomics figures:

| panel | groups |
|---|---|
| Figure 3E | `mcherry_paired-veh`, `mcherry_unpaired-veh` |
| Supplementary proteomics D | `neuron_unpaired-veh`, `neuropil_unpaired-veh` |

plus one panel per group, the full `processed_protein_ranks_animal_level.csv`, the compact
`rank_abundance_run_provenance.csv` input/design/mapping manifest, and the marker validation
workbook. Outputs go to `03_output/qc/rank_abundance/`; redirect with
`PROTEOMICS_RANK_ABUNDANCE_OUTPUT_DIR`. Contracts live in
[`tests/test_rank_abundance_animal_level.R`](tests/test_rank_abundance_animal_level.R), which also
always executes valid 48-observation and invalid 96-observation GCT fixtures and, when the shared
reference is available, requires full keyed equality to the finalized animal-level source data.

### Known unrunnable stages

Four stages cannot run against the current project tree. `run_pipeline_check.ps1` records each as
`SKIP` with a reason; the same list is reproduced here so it is discoverable without reading the
PowerShell runner. None of these blocks the animal-level pipeline, which runs end to end.

| Stage | Why it is skipped |
|---|---|
| `01_preprocessing/01_impute.r` | input `pg.matrix_raw.tsv` absent from the project tree |
| `01_preprocessing/04_format_metadata.r` | input `sample_metadata.xlsx` absent from the project tree |
| `03_qc_exploration/01_qc_protein_peptide_plot.r` | input `quicksearch.stats.annotated.xlsx` has never existed for this project — see the caveat in [`CANONICAL_OUTPUTS.md`](CANONICAL_OUTPUTS.md), which also gives the recovery route |
| `99_out_of_scope/05_metadata_create_EXP9.r` | out of scope: belongs to Exp9_Social-Stress, quarantined and guarded by `PROTEOMICS_ALLOW_EXP9` |

Each script stops with an explanatory error naming the missing input and the environment variable
that redirects it, so a skipped stage is never silently a no-op.

### Manuscript figure revision

`06_manuscript_figure_revision/` holds the code that regenerated the manuscript proteomics figures
after the statistical unit was corrected from hemisphere-level to animal-level. It is a
version-control snapshot rather than a runnable copy — the scripts still resolve their paths to
the shared drive. See [`06_manuscript_figure_revision/README.md`](06_manuscript_figure_revision/README.md).

### Publication and data deposition

`07_publication_release/` builds the publication package: editor source data, PRIDE/SDRF
metadata, data lineage, exact software versions and a SHA-256 release manifest. It **reads**
validated canonical outputs and repackages them — no model is fitted, no enrichment is run,
no PCA is recomputed, and a test enforces that by parsing every builder. Its output root is
validated before anything is written and refuses to overlap `01_input/`, `02_data/`,
`99_historical/` or any canonical `03_output` analysis branch.

```bash
PROTEOMICS_RELEASE_OUTPUT_ROOT="$TEMP/amp_release_$(date +%Y%m%d_%H%M%S)" \
  Rscript 07_publication_release/run_release.R
PROTEOMICS_RELEASE_OUTPUT_ROOT="<same root>" \
  Rscript 07_publication_release/13_validate_release.R
```

Two reporting items it surfaced, neither a consequence of the animal-level correction, are
now **resolved in the reporting layer** — no scientific value changed for either:

- **The six C46/C47 sample classes are a confirmed intentional correction**, not an open
  discrepancy. The left hemispheres of C46 and C47 were cyclically reassigned
  `neuropil → mcherry → neuron → neuropil` during historical sample-identity QC. Only
  metadata changed; acquisition identities and quantitative abundance profiles are
  unchanged, and the affected animal-level units were built from the corrected assignments
  in the first place. Both the original and the analysis-time labels are published.
- **The canonical `logFC` is a standardized abundance difference (SD units)**, not a log2
  fold change: the analysis matrix is standardised separately for each protein. Publication
  exports carry it as `effect_size_sd_units` with `source_statistic_field = "logFC"`;
  internal canonical columns and filenames keep the historical token so provenance holds.

See [`07_publication_release/README.md`](07_publication_release/README.md). The PRIDE
deposition remains `PRIDE_METADATA_INCOMPLETE` for genuinely absent acquisition metadata.

### Out-of-scope code

`99_out_of_scope/` holds code that lives here for historical reasons but is **not** part of this
analysis. Nothing in it reads or writes this project's data, nothing is invoked by `run_pipeline_check.ps1`,
and everything in it refuses to run without an explicit opt-in.

Currently one script: `05_metadata_create_EXP9.r`, which belongs to **Exp9_Social-Stress** and
writes back into that project's folder. It was inherited when this repository was seeded from Exp9
and until 2026-08-27 sat at `01_preprocessing/05_metadata_create.r` inside the numbered active
sequence, unguarded. See [`99_out_of_scope/README.md`](99_out_of_scope/README.md).

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

`01_preprocessing/03_gct_extractR.r` splits the resulting animal-level ProTigy statistical-results GCT into 12 validated forward and 12 derived reverse tables under `02_data/animal_level/split`. It reads statistical fields from their physical GCT positions rather than treating the file as an expression GCT. In split CSVs, the compatibility column `gene_symbol` retains the original GCT `id` (a UniProt-style protein identifier) because the current MapThatProt stage requires that historical column name; it is not the GCT `Description`. `Description` is preserved separately. Override the defaults with `PROTEOMICS_PROTIGY_STAT_GCT_INPUT` and `PROTEOMICS_PROTIGY_STAT_GCT_OUTPUT_ROOT`.

`02_id_mapping/01_MapThatProt_batch.r` consumes only the 12 files listed in that split branch's `indexComparisons.csv`; it does not scan the historical `99_historical/datasets_raw` tree. Forward mapping is the default, with `PROTEOMICS_MAPTHATPROT_DIRECTION=reverse` available only as an explicit override. Outputs are isolated under `02_data/animal_level/mapped`, including 12 canonical mapped CSVs in `forward/`, per-row mapping/unmapped audits, an `indexMappedComparisons.csv`, reference and source SHA-256 provenance, and the existing mapping-strategy QC report. Override roots with `PROTEOMICS_MAPTHATPROT_SPLIT_ROOT` and `PROTEOMICS_MAPTHATPROT_OUTPUT_ROOT`. `PROTEOMICS_MAPTHATPROT_REFERENCE_FILE` and `PROTEOMICS_MAPTHATPROT_MANUAL_MAPPING_FILE` override the mapping references; `PROTEOMICS_MAPTHATPROT_REFERENCE_VERSION` can record a known UniProt release label because the local three-column idmapping file does not encode one internally.

In mapped comparison CSVs, the first column `gene_symbol` is retained for current downstream compatibility but remains historically misnamed: it contains the resolved UniProt accession used by clusterProfiler. `original_protein_id` preserves the original GCT `id`, `Description` preserves the source gene-style description, and `mapped_gene_symbol` provides the mapped gene annotation when available. All source DA statistics are carried through unchanged. Unmapped rows are written separately and participate in an exact mapped-plus-unmapped row-accounting audit.

`04_differential_expression_enrichment/01_clusterProfiler.r` consumes the 12 forward files recorded in `02_data/animal_level/mapped/indexMappedComparisons.csv` and writes an isolated canonical branch under `03_output/enrichment`. It uses `uniprot_accession` explicitly. Repeated accessions are collapsed once, deterministically, by the largest absolute `log2fc` while preserving its sign, with source-row and protein-ID tie breakers and a per-comparison audit. The same selected source row is then used for every rank sensitivity. Canonical GO and KEGG GSEA rank by the ProTigy/limma moderated `t` statistic by default; positive values still mean higher abundance in the indexed numerator. The standard `t` run also writes `log2fc` sensitivity results (`GSEA_GO_BP_log2fc_sensitivity.csv` and `GSEA_KEGG_log2fc_sensitivity.csv`) with the established comparison/analysis seeds. `PROTEOMICS_ENRICHMENT_GSEA_RANK=log2fc` provides an explicit legacy-primary run without duplicating that sensitivity output. GO ORA is unchanged: it uses `padj` and signed `log2fc` foreground definitions and all unique successfully mapped, measured UniProt accessions in that comparison as its explicit universe. Rank source, analysis role, tie diagnostics, deterministic seeds, parameters, input hashes, counts, warnings, and output paths are recorded in `indexEnrichmentComparisons.csv` and per-comparison audits. Note on units: `log2fc` here is the ProTigy/limma coefficient on a per-protein standardised abundance scale — a standardized abundance difference in SD units, not a log2 fold change. The internal column name is retained for provenance; publication exports rename it to `effect_size_sd_units`.

The downstream `02_compareGO.r`, `03_compare_pathways.r`, and `04_compare_sig_expr.r` scripts consume the canonical indexes instead of historical filenames. Historical comparison aliases remain provenance metadata only. `04_compare_sig_expr.r` now compares animal-level contrast statistics; it does not fabricate per-animal expression from mapped differential-result tables.

Canonical enrichment defaults can be overridden with `PROTEOMICS_ENRICHMENT_MAPPED_ROOT`, `PROTEOMICS_ENRICHMENT_MAPPED_INDEX`, `PROTEOMICS_ENRICHMENT_OUTPUT_ROOT`, and `PROTEOMICS_ENRICHMENT_GSEA_RANK`. Existing isolated outputs are protected unless `PROTEOMICS_ENRICHMENT_FORCE=true` is set explicitly. Historical roots under `99_historical/` (`datasets_raw`, `datasets_mapped`, `core_enrichment`, `plots_pairwise`) are rejected as canonical inputs or outputs.

The 2024 instrument IDs for this dataset do not contain `_L_`/`_R_`. For that historical dataset only, the stage requires the original annotation's explicit `_left`/`_right` label and independently verifies its one-to-one agreement with both archived `ReplicateGroup` fields and the instrument sample identity. Reusable aggregation helpers otherwise accept only `Left`/`L` and `Right`/`R` and require the sample-ID hemisphere token.

The EWCE differential branch consumes `protigy_input_animal_level/neha_protigy_input_animal_level_primary.gct`, the validated handoff containing one equal-weight Left/Right mean per `AnimalID × sample_class`. It fits four 12-animal limma models and the 12 forward comparisons from `primary_contrast_manifest()`, requiring exactly three animals per condition before fitting. It does not repeat aggregation or add normalization, filtering, or imputation to the animal-level abundances; the existing EWCE gene-symbol annotation step is otherwise unchanged. Differential sampling-unit provenance and the existing EWCE parameters are written to `03_QC_Mapping_Logs/animal_level_differential_audit.csv`. Override the input and isolated output root with `PROTEOMICS_EWCE_ANIMAL_LEVEL_INPUT` and `PROTEOMICS_EWCE_OUTPUT_ROOT`; the historical `99_historical/ewce_legacy` root is rejected.

`03_qc_exploration/06_pcaPlot_animal_level.r` remains the PCA entry point but is now a thin orchestrator: the former 2,576-line monolith was split (2026-08-26) into ordered parts under `03_qc_exploration/pca/` (`06a_pca_core.r` plus seven extension parts). The parts run at top level and share the globals created by the core, so behaviour is unchanged — a split-vs-validated run produces an identical 140-file output inventory. Part order matters: `06b` creates the `rot`/`cors_df`/`um`/`npc` objects that `06c`, `06g` and `06h` consume. Each part is now guarded separately, so one failing visualization no longer suppresses every later one (previously all extensions shared a single `tryCatch`); failure still writes the failure audit and exits non-zero.

The downstream biological PCA in `03_qc_exploration/06_pcaPlot_animal_level.r` reads that same validated primary animal-level GCT with `validate_protigy_gct_v13()`. Each PCA point is one `AnimalID × sample_class` unit; all 48 units and three animals per sample-class/condition are required. The processed/imputed animal-level abundances are not normalized, filtered for missingness, or imputed again. The primary `prcomp()` remains centered and scaled, with only zero-variance protein rows removed because scaled PCA cannot use them; every removal is audited. Animal-level plots include sample class, condition, AnimalID, and phenotype, and a machine-readable audit is written under the isolated `03_output/pca` root. Override paths with `PROTEOMICS_PCA_ANIMAL_LEVEL_INPUT` and `PROTEOMICS_PCA_OUTPUT_ROOT`; the historical `99_historical/pca_plots_legacy` root is rejected. Acquisition-level protein/peptide QC and rank-abundance scripts remain sample-level because they evaluate technical acquisition and sample-class intensity behavior rather than downstream biological independence.

## License

This repository uses the MIT License. See `../LICENSE`.

## Citation

Citation metadata is provided in `../CITATION.cff`.
