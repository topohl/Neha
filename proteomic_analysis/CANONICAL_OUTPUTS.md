# Canonical inputs & outputs (Neha proteomics)

Read-side companion to the write-side guards in `R/*_utils.R`. Those guards stop you
*writing* into a historical root; this file tells you which existing folder to *read*.

Data root: `S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/`

Restructured **2026-08-26** into `01_input` / `02_data` / `03_output` / `99_historical`.
The migration was folder renames only (no file contents changed except recorded paths inside
index/manifest CSVs). Reversible manifest + originals: `_migration_backup_20260826/`.

## Layout

| Folder | Contains | Written by pipeline? |
|---|---|---|
| `01_input/` | external / hand-placed inputs: `references/` (UniProt idmapping, manual mapping, celltypes), `metadata/` (`sample_info.xlsx`), `raw_proteomics/` (raw pg.matrix drops), `single_cell/` (loom) | **No** — never write here |
| `02_data/` | derived intermediates: `gct/` (processed/imputed matrices), `animal_level/` (`input_gct/`, `split/`, `mapped/`, `audits/`, `qc_exports/`) | Yes |
| `03_output/` | analysis outputs: `enrichment/`, `ewce/`, `pca/`, `synthesis/`, `inferential_checks/`, `qc/` | Yes |
| `99_historical/` | superseded hemisphere-level era + prior generations + archives | **No** — read-only by convention |

## Authoritative artefact per stage

| Stage | **Use this** | Superseded / do not use | Why |
|---|---|---|---|
| Animal-level GCT (n=3/group, L/R averaged) | `02_data/animal_level/input_gct/neha_protigy_input_animal_level_primary.gct`<br>SHA256 `f12cf99e…b3bbb3` | any hemisphere-level `pg.matrix_*` in `02_data/gct/` | The validated 48-observation animal-level matrix. Hemisphere-level matrices treat L/R as independent (pseudoreplication). |
| Split stats | `02_data/animal_level/split/` (`indexComparisons.csv`) | `99_historical/datasets_raw/` | 12 canonical forward contrasts. |
| ID mapping | `02_data/animal_level/mapped/` (`indexMappedComparisons.csv`) | `99_historical/datasets_mapped/`, `datasets_unmapped/`, `mapping_reports/` | Current MapThatProt run. |
| **Enrichment (GSEA/ORA)** | **`03_output/enrichment/enrichment_t_rank_validation_20260825/`** | `03_output/enrichment/enrichment_log2fc_rank_legacy/` | **Canonical GSEA ranks by the moderated `t` statistic.** The legacy folder ranked by `log2fc` (it was named just `enrichment/` before the restructure — renamed to make the difference explicit). log2fc-ranked results exist only as a *sensitivity* analysis. |
| EWCE | `03_output/ewce/EWCE_Results_animal_level_validation_20260825/` | `99_historical/ewce_legacy/` | Animal-level differential signatures; primary settings annotation level 2, top-N 250, FDR 0.05. |
| PCA | `03_output/pca/pca_plots_animal_level_validation_20260825_rerun/` | `..._validation_20260825/` (first pass), `99_historical/pca_plots_legacy/` | Three generations exist; the `_rerun` folder is current and also holds `qc_outlier_followup/`. |
| Biological synthesis | `03_output/synthesis/animal_level_biological_synthesis_20260826/` | — | Cross-modal integration of DA + GSEA + ORA + EWCE + PCA over the 12 contrasts. |
| Inferential checks | `03_output/inferential_checks/final_inferential_checks_20260826/` | — | mcherry learning-contrast diagnostics, post-hoc Pairing×CNO interaction, design-identifiability and collection-plate provenance audits. |
| Acquisition QC | `03_output/qc/` | — | Note: `quicksearch.stats.annotated.xlsx` is **absent** from the project tree, so acquisition-level technical QC has never been evaluated. See caveats. |

## Naming gotcha: `bg` means neuropil

Historical filenames and metadata use **`bg`** (from "background") where current outputs say
**`neuropil`**. So `bg1_bg2.csv`, `bg2_bg4.csv`, and `group2` values like `bg_1` all refer to
the **neuropil tissue compartment** — *not* a negative control, blank, or background
subtraction. The mapping lives in `R/analysis_labels.R` (`sample_class_aliases`), but it is
invisible when reading raw filenames, so treat any `bg*` file as neuropil.

Sample classes: `mcherry`, `neuropil` (= `bg`), `cfos`, `neuron`.
Condition codes: `1`=paired_cno, `2`=paired_veh, `3`=unpaired_cno, `4`=unpaired_veh.

## Caveats carried forward from the audits

- **Collection/acquisition plate is confounded with Pairing.** 11 of 12 animals follow
  paired→Plate1 / unpaired→Plate2. The `paired_veh` vs `unpaired_veh` ("learning") contrast is
  *completely* aliased with plate in every sample class. `tests/test_design_balance.R` asserts
  this so it cannot be rediscovered by accident.
- **The mcherry learning contrast (1132 FDR proteins) must not be read as biology** without
  further technical investigation — broad one-directional shift, no GSEA/ORA pathway coherence,
  and complete plate aliasing. See `03_output/inferential_checks/`.
- **The Pairing×CNO interaction is post-hoc and n=3/cell.** It cannot be distinguished from a
  Plate×Treatment technical effect, because Plate and Pairing are aliased.
- **`quicksearch.stats.annotated.xlsx` never existed for this project** (resolved 2026-08-27).
  A filesystem-wide search found that filename only under `Exp9_Social-Stress`, never under
  `Collabs/Neha`. The old hardcoded default —
  `<Neha>/clusterProfiler/Datasets/pg_matrix/raw/quicksearch.stats.annotated.xlsx` — has a tail
  byte-identical to the real Exp9 path
  `<Exp9_Social-Stress>/proteomics/Datasets/pg_matrix/raw/quicksearch.stats.annotated.xlsx`,
  so it appears to be boilerplate inherited when this repo was seeded from that project (the
  same reason `01_preprocessing/05_metadata_create.r` still points at `Exp9_Social-Stress`).
  The Exp9 file has the correct 33 columns but contains **0 of 326** Neha samples — different
  instrument (`Bluto` vs `Olive`) and acquisition date (2025-07-03 vs 2024-12-17) — so it is
  **not** a usable substitute. Consequence: acquisition-level technical QC (identification
  depth, MS1/MS2 signal, mass accuracy, normalisation instability) has never been evaluable
  for this dataset, which is why every audit reported it as unavailable. To enable it, re-export
  the QC report for the Neha acquisition and point `NEHA_QC_QUICKSEARCH_STATS` at it.
- **`03_qc_exploration/01_qc_protein_peptide_plot.r` plots instrument QC metrics, not protein
  abundance.** Its PCA is a PCA *of acquisition QC metrics*. For protein-abundance PCA use
  `03_qc_exploration/06_pcaPlot_Neha.r`.
- `99_historical/datasets_unmapped/` holds header-only CSVs (`gene_symbol` + newline). These are
  legitimate "zero unmapped proteins for this comparison" records, not stray files.

## Overriding defaults

Every stage resolves paths via `getOption()` then an env var then the default above, e.g.
`NEHA_PCA_ANIMAL_LEVEL_INPUT`, `NEHA_PCA_OUTPUT_ROOT`, `NEHA_EWCE_ANIMAL_LEVEL_INPUT`,
`NEHA_EWCE_OUTPUT_ROOT`, `NEHA_ENRICHMENT_MAPPED_ROOT`, `NEHA_ENRICHMENT_OUTPUT_ROOT`,
`NEHA_MAPTHATPROT_SPLIT_ROOT`, `NEHA_MAPTHATPROT_OUTPUT_ROOT`, `NEHA_PROTIGY_STAT_GCT_INPUT`,
`NEHA_ANIMAL_LEVEL_DATA_ROOT`. Before the restructure several defaults pointed at folders that
did not exist and every validated run relied on overrides; the defaults are now correct.
