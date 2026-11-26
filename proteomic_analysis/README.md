# Associative Memory Proteomics

Analysis code for the Associative Memory Proteomics data set: 12 mice, three per condition, each giving four laser-captured sample classes, with left and right hemispheres averaged per animal. The code covers preprocessing, differential abundance, ID mapping, GO and KEGG enrichment, cell-type enrichment (EWCE), QC and the data release. The data and results are on the lab shared drive (see [Data on the shared drive](#data-on-the-shared-drive)); without it only the demo and the parts of the tests that use small fixtures will run.

Sample classes are `mcherry`, `neuropil`, `cfos` and `neuron`. Condition codes are `1` = `paired_cno`, `2` = `paired_veh`, `3` = `unpaired_cno` and `4` = `unpaired_veh`. Each sample class has three contrasts (`paired_cno` vs `paired_veh`, `paired_veh` vs `unpaired_veh`, `unpaired_cno` vs `unpaired_veh`), 12 in all. All scripts take these labels from `R/analysis_labels.R`.

Older filenames and metadata say `bg` (background) for the neuropil compartment, for example `bg1_bg2.csv` or `group2` = `bg_1`. It is not a blank, a negative control or a background subtraction. The alias table is `sample_class_aliases` in `R/analysis_labels.R`.

## Requirements and installation

Tested with R 4.5.1 on Windows 11 and on the Ubuntu GitHub Actions runner. No special hardware is needed.

```bash
git clone https://github.com/topohl/Neha.git
cd Neha/proteomic_analysis
```

The demo and the tests run in CI need only base R, so with R installed setup takes under a minute. The pipeline stages need CRAN and Bioconductor packages; [`requirements_R.md`](requirements_R.md) lists them per stage with the versions used. Installing the Bioconductor packages (EWCE, clusterProfiler and their annotation data) can take much longer, depending on the network.

## Demo

From `proteomic_analysis/`:

```bash
Rscript run_demo.R
```

`run_demo.R` reads `demo/input/demo_sample_metadata.csv` and `demo/input/demo_pg_matrix.csv`, parses sample classes, condition codes and comparison names with `R/analysis_labels.R`, and writes `demo/output/demo_result_table.csv` and `demo/output/demo_comparison_table.csv`. It stops with an error if either table differs from its copy in `demo/expected_output/`. It takes well under a minute.

### Running on your own data

`01_preprocessing/02_excel_convert.r` adds the sample metadata as header rows to a protein table. Run from `proteomic_analysis/`, it uses the demo files by default. For other files set these R options or environment variables:

| R option | Environment variable | Sets |
|---|---|---|
| `proteomics.metadata_path` | `PROTEOMICS_METADATA_PATH` | sample metadata (CSV or Excel) |
| `proteomics.excel_convert_file` | `PROTEOMICS_EXCEL_CONVERT_FILE` | input table (`excel` mode) |
| `proteomics.excel_convert_folder` | `PROTEOMICS_EXCEL_CONVERT_FOLDER` | folder of `.xlsx`/`.csv` tables (`folder` mode) |
| `proteomics.excel_convert_output` | `PROTEOMICS_EXCEL_CONVERT_OUTPUT` | output folder |
| `proteomics.excel_convert_mode` | `PROTEOMICS_EXCEL_CONVERT_MODE` | `excel` (default) or `folder` |

The metadata needs a `sample_id` column. If `sample_class` or `condition_code` is missing, it is parsed from the sample IDs (`condition_code` from `condition` when that column exists).

## Pipeline

Run each script with `Rscript` from `proteomic_analysis/`. Paths come from environment variables (most stages also accept an R option) and default to the shared drive. Each script sets them near its top; the main variables are in the tables below.

### Main chain

Since 2026-11-24 the primary statistics come from one limma model per sample class (12 animals, n = 3 per condition). The earlier model, run in ProTigy, pooled all four sample classes; `CANONICAL_OUTPUTS.md` gives the reasons for the change and what it changed.

| Step | Script | Default output | Main overrides |
|---|---|---|---|
| 1 | `01_preprocessing/02a_prepare_animal_level_protigy_input.r` | `02_data/animal_level/input_gct/` | `PROTEOMICS_ANIMAL_LEVEL_OUTPUT_DIR` |
| 2 | `01_preprocessing/02b_limma_stat_gct.r` | `02_data/animal_level/within_compartment_20261124/stat_results_within_compartment.gct` | `PROTEOMICS_LIMMA_OUTPUT_GCT` |
| 3 | `01_preprocessing/03_gct_extractR.r` | `02_data/animal_level/within_compartment_20261124/split/` | `PROTEOMICS_PROTIGY_STAT_GCT_INPUT`, `PROTEOMICS_PROTIGY_STAT_GCT_OUTPUT_ROOT` |
| 4 | `02_id_mapping/01_MapThatProt_batch.r` | `02_data/animal_level/within_compartment_20261124/mapped/` | `PROTEOMICS_MAPTHATPROT_SPLIT_ROOT`, `PROTEOMICS_MAPTHATPROT_OUTPUT_ROOT` |
| 5 | `04_differential_expression_enrichment/01_clusterProfiler.r` | none, `PROTEOMICS_ENRICHMENT_OUTPUT_ROOT` must be set | `PROTEOMICS_ENRICHMENT_MAPPED_ROOT` |
| 6 | `04_differential_expression_enrichment/02_compareGO.r`, `03_compare_pathways.r`, `04_compare_sig_expr.r` | inside the enrichment run folder | as step 5 |
| 7 | `04_differential_expression_enrichment/05_pairing_cno_interaction.r`, then `06_mcherry_learning_robustness.r` | `03_output/inferential_checks/within_compartment_20261124/` | `PROTEOMICS_INFERENTIAL_OUTPUT_ROOT` |

1. Averages the valid Left and Right hemisphere values within each `AnimalID × sample_class` (48 columns, 5,349 proteins) without transforming, filtering, imputing or remapping, and writes the matrix with audits and SHA-256 provenance. The 2024 instrument sample IDs carry no `_L_`/`_R_`, so for this data set the hemisphere comes from the `_left`/`_right` label in the original annotation, checked against the archived `ReplicateGroup` fields and the instrument sample identity.
2. Fits `~0 + condition` within each sample class with `eBayes(robust = TRUE, trend = FALSE)` and Benjamini-Hochberg within each contrast, and writes the results in the ProTigy GCT layout with a run manifest. `PROTEOMICS_LIMMA_VARIANCE_MODEL=pooled` fits the old pooled model instead. An existing output is replaced only with `PROTEOMICS_LIMMA_FORCE=true`.
3. Splits the statistics GCT into 12 forward and 12 reverse contrast tables plus `indexComparisons.csv`. An existing split is replaced only with `PROTEOMICS_PROTIGY_STAT_GCT_FORCE=true`.
4. Maps the protein IDs of the 12 forward tables listed in the split's index to UniProt accessions and gene symbols (`PROTEOMICS_MAPTHATPROT_DIRECTION=reverse` maps the reverse tables). It reads `within_compartment_20261124/split/` by default (the pooled run's tables are in `split/` and `mapped/`), and an existing mapped index is replaced only with `PROTEOMICS_MAPTHATPROT_FORCE=true`. `PROTEOMICS_MAPTHATPROT_REFERENCE_FILE` and `PROTEOMICS_MAPTHATPROT_MANUAL_MAPPING_FILE` set the mapping references, and `PROTEOMICS_MAPTHATPROT_REFERENCE_VERSION` records the UniProt release, which the idmapping file does not contain.
5. Runs GO and KEGG GSEA and GO ORA for each contrast, reading `within_compartment_20261124/mapped/` by default. Repeated accessions are collapsed to the row with the largest absolute `log2fc`. GSEA ranks by the moderated `t` (positive means higher in the numerator) and also writes a `log2fc`-ranked sensitivity run; `PROTEOMICS_ENRICHMENT_GSEA_RANK=log2fc` makes `log2fc` the main ranking. The ORA universe is every mapped, measured accession in the contrast. `PROTEOMICS_ENRICHMENT_ONTOLOGY` picks BP (default), CC or MF, `PROTEOMICS_ENRICHMENT_KEGG=false` turns KEGG off, and `PROTEOMICS_ENRICHMENT_NPERM_SIMPLE` sets fgsea `nPermSimple` (unset means the fgsea default, 1000). Gene sets that fgsea leaves without a P value are listed in `audits/<comparison>/gsea_untested_sets.csv`. An existing run folder is reused only with `PROTEOMICS_ENRICHMENT_FORCE=true`. `CANONICAL_OUTPUTS.md` says which enrichment runs are current (GO from `enrichment_within_compartment_20261125/`, KEGG from `enrichment_within_compartment_20261124/`), how they were made and why KEGG is not version-pinned.
6. Compares GO terms, pathways and significant proteins across contrasts. `04_compare_sig_expr.r` compares contrast statistics; it does not rebuild per-animal values.
7. Tests the post-hoc Pairing × CNO interaction within each sample class, then reruns the mCherry learning contrast with each animal left out and compares t- and effect-ranked GSEA. Both refuse to write into a used folder without `PROTEOMICS_INFERENTIAL_FORCE=true`.

The `logFC` and `log2fc` columns are differences in standardized abundance (SD units), because the analysis matrix is scaled per protein; they are not log2 fold changes. The release exports them as `effect_size_sd_units`. One column name is misleading: `gene_symbol` holds the GCT `id` (a UniProt entry name such as `NOVA2_MOUSE`) in the split tables and the UniProt accession in the mapped tables; the gene symbol is in `mapped_gene_symbol`. The old names are kept because downstream code and recorded hashes depend on them.

To rerun the chain without touching the shared drive, send each step to a scratch folder and point the next step at it:

```bash
OUT="$TEMP/amp_chain_$(date +%Y%m%d_%H%M%S)"
PROTEOMICS_LIMMA_OUTPUT_GCT="$OUT/stat_results_within_compartment.gct" \
  Rscript 01_preprocessing/02b_limma_stat_gct.r
PROTEOMICS_PROTIGY_STAT_GCT_INPUT="$OUT/stat_results_within_compartment.gct" \
PROTEOMICS_PROTIGY_STAT_GCT_OUTPUT_ROOT="$OUT/split" \
  Rscript 01_preprocessing/03_gct_extractR.r
PROTEOMICS_MAPTHATPROT_SPLIT_ROOT="$OUT/split" PROTEOMICS_MAPTHATPROT_OUTPUT_ROOT="$OUT/mapped" \
  Rscript 02_id_mapping/01_MapThatProt_batch.r
PROTEOMICS_ENRICHMENT_MAPPED_ROOT="$OUT/mapped" PROTEOMICS_ENRICHMENT_OUTPUT_ROOT="$OUT/enrichment" \
  Rscript 04_differential_expression_enrichment/01_clusterProfiler.r
```

Step 1 and the three branch scripts below have no overwrite check and write straight into their shared-drive defaults, so set their output variable when rerunning them.

### QC and cell-type branches

These read the animal-level GCT from step 1.

| Script | Default output | Main overrides |
|---|---|---|
| `03_qc_exploration/06_pcaPlot_animal_level.r` | `03_output/pca/` | `PROTEOMICS_PCA_ANIMAL_LEVEL_INPUT`, `PROTEOMICS_PCA_OUTPUT_ROOT` |
| `03_qc_exploration/02_rank_abundance_by_sample_class.r` | `03_output/qc/rank_abundance/` | `PROTEOMICS_RANK_ABUNDANCE_INPUT_GCT`, `PROTEOMICS_RANK_ABUNDANCE_ID_MAP`, `PROTEOMICS_RANK_ABUNDANCE_OUTPUT_DIR` |
| `05_celltype_enrichment_EWCE/01_EWCE.r` | `03_output/ewce/` | `PROTEOMICS_EWCE_ANIMAL_LEVEL_INPUT`, `PROTEOMICS_EWCE_OUTPUT_ROOT` |

- PCA: one point per `AnimalID × sample_class` (48), centered and scaled, with zero-variance proteins dropped and listed in the audit. `06_pcaPlot_animal_level.r` runs `pca/06a_pca_core.r` and then the parts `06b` to `06h` in order (`06b` creates objects that `06c`, `06g` and `06h` use). `06a` must succeed; if one of `06b` to `06h` fails, the others still run and the script exits with an error at the end.
- Rank abundance: groups by `sample_class × condition` (16 groups of 3 animals) and writes the panels `Fig3E` (mCherry, paired vs unpaired vehicle) and `SuppD` (neuron and neuropil, unpaired vehicle), a panel for each of the 16 groups, `processed_protein_ranks_animal_level.csv`, `rank_abundance_run_provenance.csv` and a marker workbook. Gene names come from one mapped table of the pooled run, used only as an ID lookup.
- EWCE: fits its own limma model per sample class (`eBayes(trend = TRUE, robust = TRUE)`, so close to but not the same as step 2) for the 12 contrasts and runs 10,000 bootstrap repetitions, which takes hours. The sampling-unit audit is `03_QC_Mapping_Logs/animal_level_differential_audit.csv`.

### Stages that do not run

| Script | Reason |
|---|---|
| `01_preprocessing/01_impute.r` | its input `pg.matrix_raw.tsv` is not in the data folder |
| `01_preprocessing/04_format_metadata.r` | its input `sample_metadata.xlsx` is not in the data folder |
| `03_qc_exploration/01_qc_protein_peptide_plot.r` | `quicksearch.stats.annotated.xlsx` was never produced for this data set, so acquisition-level QC has not been done (see `CANONICAL_OUTPUTS.md`) |

Each stops with an error that names the missing file and the variable that sets its path. Step 1 starts from the processed matrices in `02_data/gct/`, so the main chain does not need them.

`99_out_of_scope/05_metadata_create_EXP9.r` belongs to the Exp9_Social-Stress project and writes into that project's folder. It runs only with `PROTEOMICS_ALLOW_EXP9=true`; see [`99_out_of_scope/README.md`](99_out_of_scope/README.md).

### Figure scripts

`06_manuscript_figure_revision/` holds byte-identical copies of the figure scripts that were run from the shared drive, including the per-compartment rebuilds (`within_compartment_20261126/` is current; `run_rebuild.sh` runs it). They still read and write shared-drive paths, so keep them as a record and do not edit them: `SCRIPT_PROVENANCE.csv` holds the SHA-256 of each copy and of the script that ran. See its [README](06_manuscript_figure_revision/README.md).

## Data on the shared drive

Inputs and outputs are under `S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/`, in this layout since 2026-08-26:

| Folder | Contents | Written by the pipeline |
|---|---|---|
| `01_input/` | hand-placed inputs: `references/` (UniProt idmapping, manual mapping, cell types), `metadata/` (`sample_info.xlsx`), `raw_proteomics/` (raw pg.matrix drops), `single_cell/` (loom) | no |
| `02_data/` | intermediates: `gct/` (processed and imputed matrices), `animal_level/` (`input_gct/`, `split/`, `mapped/`, `within_compartment_20261124/`, `audits/`, `qc_exports/`) | yes |
| `03_output/` | results: `enrichment/`, `ewce/`, `pca/`, `inferential_checks/`, `synthesis/`, `qc/`, figure rebuilds and releases | yes |
| `99_historical/` | the superseded hemisphere-level analysis and older runs | no |

Several stages have more than one run on disk. [`CANONICAL_OUTPUTS.md`](CANONICAL_OUTPUTS.md) says which one is current and lists the caveats that come with the data. The 2026-08-26 reorganization renamed folders and rewrote the paths recorded in index and manifest CSVs; the manifest to undo it and the original files are in `_migration_backup_20260826/`.

## Tests

Run a test with `Rscript tests/<file>` from `proteomic_analysis/` (`test_runner_safety.ps1` with `powershell -ExecutionPolicy Bypass -File tests\test_runner_safety.ps1`). Parts that need the shared drive are skipped when it is not reachable. CI (`.github/workflows/proteomic-analysis-demo.yml`) runs `check_stale_labels.R`, `test_gct_comparison_parser.R`, `test_deprecated_path_utils_shim.R` and the demo on every push and pull request.

| Test | Checks |
|---|---|
| `check_stale_labels.R` | retired labels and naming rules in active code and docs |
| `test_gct_comparison_parser.R` | parsing of comparison names |
| `test_animal_level_proteomics.R` | Left/Right averaging into animal-level units |
| `test_limma_stat_gct.R` | step 2: `pooled` mode reproduces the ProTigy export, the per-compartment mode reproduces the locked statistics GCT (a different R or limma version is reported, not failed) |
| `test_mapthatprot_animal_level.R` | ID mapping reads only the 12 indexed contrasts, uses the per-compartment folders by default and refuses to replace an existing mapped index |
| `test_animal_level_enrichment.R` | contrast manifest, indexed inputs and settings of the enrichment stage |
| `test_enrichment_untested_sets.R` | `PROTEOMICS_ENRICHMENT_NPERM_SIMPLE` and the list of untested gene sets |
| `test_pairing_cno_interaction.R` | the interaction contrast, its output folder checks and a 300-protein run |
| `test_ewce_animal_level.R` | EWCE uses one value per animal and sample class |
| `test_pca_animal_level.R`, `test_pca_optional_visualizations.R` | PCA input checks; optional plots with empty input |
| `test_rank_abundance_animal_level.R` | rank-abundance stage on 48- and 96-column fixtures, and equality with the reference source data |
| `test_data_integrity.R` | SHA-256 of the validated GCT and of every file recorded in the indexes |
| `test_design_balance.R` | collection plate follows pairing in 11 of 12 animals and in every animal of the `paired_veh` vs `unpaired_veh` contrast, so plate and pairing cannot be separated there (see `CANONICAL_OUTPUTS.md`) |
| `test_project_file_sha256.R` | the file-hashing helper |
| `test_deprecated_path_utils_shim.R` | `R/neha_path_utils.R`, kept only for the scripts in `06_manuscript_figure_revision/`, is not sourced by active code |
| `test_runner_safety.ps1` | the safety checks of the pipeline check, on fixtures |
| `07_publication_release/tests/` | the release layer; tests that need a built release skip when none is found |

### Pipeline check

```powershell
powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1             # Normal: tests, GCT rebuild, split, PCA, rank abundance (10-20 min)
powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1 -Tier Fast  # tests only, about 2 min
powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1 -Tier Full  # adds ID mapping, enrichment and EWCE (hours)
```

`run_pipeline_check.ps1` runs every test above and then the stages, with all output in a new folder under `%TEMP%`. It is Windows-only and needs the shared drive. It uses R 4.5.1 from `%LOCALAPPDATA%\Programs\R` if it is installed there, otherwise `Rscript.exe` on `PATH`. Inherited `PROTEOMICS_*` variables are cleared for the run and restored afterwards. Each step is reported as PASS, FAIL, SKIP (missing input, or the old ProTigy step, which is an external program) or ABORT (a check in `run_pipeline_safety.ps1` stopped the stage before R started), and the script exits non-zero on any FAIL or ABORT. The Normal tier also checks that the rebuilt animal-level GCT has SHA-256 `f12cf99e1bfb7c17bbf56bffb6783e924698bce5d5533a8e312bc4bbb733bbb3` and splits the per-compartment statistics; the Full tier maps and enriches the per-compartment tables. Step 2 is not a runner stage: `test_limma_stat_gct.R` rebuilds it.

## Release

`07_publication_release/` builds the data release from finished outputs: Supplementary Table 2 (`editor_source_data/Supplementary_Table_2_proteomics.xlsx`), figure source data, processed and differential tables, enrichment tables, PRIDE/SDRF metadata, data lineage, software versions and a SHA-256 manifest. It fits no model and runs no enrichment or PCA; `07_publication_release/tests/test_release_no_scientific_recomputation.R` checks the builders for such calls. Since 2026-11-24 it releases the per-compartment statistics.

Build into a scratch folder and validate it there:

```bash
PROTEOMICS_RELEASE_OUTPUT_ROOT="$TEMP/amp_release_$(date +%Y%m%d_%H%M%S)" \
  Rscript 07_publication_release/run_release.R
PROTEOMICS_RELEASE_OUTPUT_ROOT="<same root>" \
  Rscript 07_publication_release/13_validate_release.R
```

Without the variable the release goes to `03_output/release_within_compartment_20261126b/` on the shared drive. That needs `PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE=true` and a clean git tree, because the release records the commit that built it. The output folder may not overlap `01_input/`, `02_data/`, `99_historical/`, the analysis folders in `03_output/` or the frozen earlier releases (`03_output/publication_release/`, `03_output/release_within_compartment_20261125/`, `03_output/release_within_compartment_20261126/`). The C46/C47 sample-class correction, the effect-size columns and the PRIDE metadata that is still missing are described in [`07_publication_release/README.md`](07_publication_release/README.md).

## License and citation

MIT License, see [`../LICENSE`](../LICENSE). Citation metadata is in [`../CITATION.cff`](../CITATION.cff).
