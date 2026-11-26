# Inputs and outputs (Associative Memory Proteomics)

Which folder to read for each stage. The path checks in `R/*_utils.R` stop scripts writing into
old roots; this file says where the current results are.

Data root: `S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/`

Its folder layout and the 2026-08-26 reorganization are described in `README.md` (Data on the
shared drive).

## Primary statistics: within-compartment model (2026-11-24)

The ProTigy export is one limma model over all 48 animal-level columns (16 groups, residual
df 32), so the four compartments share one variance estimate per protein. The primary
statistics now fit one model per compartment: `~0 + condition` over that compartment's 12
animals (n = 3 per condition, residual df 8), eBayes robust (no trend), BH within each
contrast. Separate models are used because residual variances differ between compartments
(median mCherry 0.33, neuropil 0.35, c-Fos 0.59, neuron 0.36; per-protein max/min ratio median
4.3 against 2.8 expected under equal variances), and the four compartments of one animal are not
independent columns (duplicateCorrelation consensus 0.20). The contrast manifest's label
`within_sample_class_animal_level_two_sample_ProTigy` does not describe what ProTigy ran.

`01_preprocessing/02b_limma_stat_gct.r` writes the statistics in the ProTigy layout; its
pooled mode reproduces the ProTigy export (`tests/test_limma_stat_gct.R`). logFC is unchanged;
t, p, FDR, B and AveExpr (now the compartment mean) change. FDR < 0.05 proteins, pooled →
per compartment: mCherry learning 1137 → 1003 (998 after mapping), mCherry unpaired CNO 37 → 0,
c-Fos learning 6 → 0, c-Fos paired CNO 3 → 0.

| Stage | Per-compartment run |
|---|---|
| Statistics | `02_data/animal_level/within_compartment_20261124/stat_results_within_compartment.gct` |
| Split / ID mapping | `02_data/animal_level/within_compartment_20261124/split/`, `.../mapped/` |
| Enrichment | GO-BP and GO-CC: `03_output/enrichment/enrichment_within_compartment_20261125/{GO_BP,GO_CC}/` (fgsea `nPermSimple = 100000`, adopted 2026-11-25; ORA byte-identical to 20261124). KEGG: `03_output/enrichment/enrichment_within_compartment_20261124/GO_BP/` (fgsea default, not rerun; its GO tables are superseded) |
| Figure panels | `03_output/figure_rebuild_within_compartment_20261125/` (Fig. 2 D, F, G, H, I and Supp. Fig. 4 F, run by its `scripts/run_rebuild.sh`; repository copy in `06_manuscript_figure_revision/within_compartment_20261125/`). The 2026-11-24 rebuild is superseded; its README lists what was not rebuilt |
| Inferential checks | `03_output/inferential_checks/within_compartment_20261124/` (Pairing×CNO interaction per compartment; mCherry learning leave-one-animal-out, t- vs effect-ranked GSEA, CNO-arm opposition; written by `04_differential_expression_enrichment/05_pairing_cno_interaction.r` and `06_mcherry_learning_robustness.r`) |

GSEA term counts near q = 0.05 move by a few terms with the permutation seed (5 base seeds at
`nPermSimple = 100000`: mCherry learning GO-BP 25–26, GO-CC 72–75; neuropil learning GO-BP
108–117; neuropil paired CNO GO-BP 84–93; `enrichment_within_compartment_20261125/seed_range/`).
The pooled enrichment run used DOSE `seed = TRUE`, which reseeds every call with 10407 (the
RNG-kind code) and ignores the recorded seeds, so it can only be reproduced that way. Since
2026-11-24 the script passes `seed = FALSE` and the recorded seed is the one used. KEGG gene
sets are downloaded at run time and have no fixed version.

fgsea gives no P value for "unbalanced" gene sets (large sets against a one-sided ranking), and
DOSE leaves those rows out of its tables. The GO runs were therefore repeated at
`nPermSimple = 100000` (adopted 2026-11-25). Sets still without a P value (0–27 GO-BP and 0–13
GO-CC per contrast, none in mCherry learning) are listed in
`GO_<ONT>/audits/<contrast>/gsea_untested_sets.csv`. KEGG stays at the fgsea default (0–23
pathways per contrast without a P value, counted in the 20261124 run's
`audits/<contrast>/warnings.csv`). The Pairing×CNO interaction, leave-one-animal-out and t- vs
effect-ranked checks and the Supp. Fig. 4E paired model also use the default, and the checks are
compared with the primary run at that same setting (so did an exploratory complete-case check
that is not part of this repository).

The biological synthesis has not been rebuilt and still uses the pooled run. The
design-identifiability and collection-plate audits in `final_inferential_checks_20260826/` do not
depend on the statistics and stay current, as do the figure panels that do not use them (PCA,
rank abundance, the per-acquisition protein and peptide counts of Supp. Fig. 4 A, EWCE, the secondary within-animal model).
The release package is built on this run into
`03_output/release_within_compartment_20261126/` (Supplementary Table 2 layout; not yet deployed
to the shared drive). Until it is, the deployed per-compartment release is
`release_within_compartment_20261125/` (old workbook layout); it and the pooled-run release stay
frozen.

### Rerunning the per-compartment chain

From `proteomic_analysis/` in Git Bash, with the defaults the stages now have. Each stage
refuses to replace existing output unless its `*_FORCE` variable is `true`
(`PROTEOMICS_LIMMA_FORCE`, `PROTEOMICS_PROTIGY_STAT_GCT_FORCE`, `PROTEOMICS_MAPTHATPROT_FORCE`,
`PROTEOMICS_ENRICHMENT_FORCE`, `PROTEOMICS_INFERENTIAL_FORCE`), so point the outputs at new
folders for a rerun (`README.md` shows a scratch rerun).

```bash
E=S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/enrichment
export PROTEOMICS_ENRICHMENT_PLOTS=false   # as recorded in audits/run_parameters.csv
Rscript 01_preprocessing/02b_limma_stat_gct.r   # within_compartment_20261124/stat_results_within_compartment.gct
Rscript 01_preprocessing/03_gct_extractR.r      # within_compartment_20261124/split/
Rscript 02_id_mapping/01_MapThatProt_batch.r    # within_compartment_20261124/mapped/
# 20261124, fgsea default: the released KEGG comes from the BP run
PROTEOMICS_ENRICHMENT_ONTOLOGY=BP PROTEOMICS_ENRICHMENT_OUTPUT_ROOT=$E/enrichment_within_compartment_20261124/GO_BP \
  Rscript 04_differential_expression_enrichment/01_clusterProfiler.r
PROTEOMICS_ENRICHMENT_ONTOLOGY=CC PROTEOMICS_ENRICHMENT_KEGG=false \
  PROTEOMICS_ENRICHMENT_OUTPUT_ROOT=$E/enrichment_within_compartment_20261124/GO_CC \
  Rscript 04_differential_expression_enrichment/01_clusterProfiler.r
# 20261125: the released GO-BP and GO-CC
for ONT in BP CC; do
  PROTEOMICS_ENRICHMENT_ONTOLOGY=$ONT PROTEOMICS_ENRICHMENT_KEGG=false PROTEOMICS_ENRICHMENT_NPERM_SIMPLE=100000 \
    PROTEOMICS_ENRICHMENT_OUTPUT_ROOT=$E/enrichment_within_compartment_20261125/GO_$ONT \
    Rscript 04_differential_expression_enrichment/01_clusterProfiler.r
done
Rscript 04_differential_expression_enrichment/05_pairing_cno_interaction.r
Rscript 04_differential_expression_enrichment/06_mcherry_learning_robustness.r
```

The 20261125 run as it was executed, with the 02-04 comparison scripts and the seed-range runs,
is `run_rerun.sh` and `seed_range/` in that enrichment folder on the shared drive. The figure
rebuild is `scripts/run_rebuild.sh` in its output folder, and the release is built as described
in `07_publication_release/README.md`.

## Folder per stage (pooled run, 2026-08)

For split statistics, ID mapping, enrichment and the learning-contrast and interaction checks,
the per-compartment folders above replace these rows; the rows here stay as the pooled-run
record.

| Stage | Use this | Superseded | Notes |
|---|---|---|---|
| Animal-level GCT (n=3/group, L/R averaged) | `02_data/animal_level/input_gct/neha_protigy_input_animal_level_primary.gct`<br>SHA256 `f12cf99e…33bbb3` | any hemisphere-level `pg.matrix_*` in `02_data/gct/` | The validated 48-observation animal-level matrix. Hemisphere-level matrices treat L/R as independent (pseudoreplication). |
| Split stats | `02_data/animal_level/split/` (`indexComparisons.csv`) | `99_historical/datasets_raw/` | 12 canonical forward contrasts. |
| ID mapping | `02_data/animal_level/mapped/` (`indexMappedComparisons.csv`) | `99_historical/datasets_mapped/`, `datasets_unmapped/`, `mapping_reports/` | MapThatProt output of the pooled run. |
| Enrichment (GSEA/ORA) | `03_output/enrichment/enrichment_t_rank_validation_20260825/` | `03_output/enrichment/enrichment_log2fc_rank_legacy/` | GSEA ranks by the moderated `t` statistic. The legacy folder (called `enrichment/` before the restructure) ranked by `log2fc`, which is now only a sensitivity analysis. Folder and column names keep the old `log2fc` token (what it measures: `README.md`, Pipeline). |
| EWCE | `03_output/ewce/EWCE_Results_animal_level_validation_20260825/` | `99_historical/ewce_legacy/` | Animal-level differential signatures; primary settings annotation level 2, top-N 250, FDR 0.05. |
| PCA | `03_output/pca/pca_plots_animal_level_validation_20260825_rerun/` | `..._validation_20260825/` (first pass), `99_historical/pca_plots_legacy/` | Three generations exist; the `_rerun` folder is current and also holds `qc_outlier_followup/`. |
| Biological synthesis | `03_output/synthesis/animal_level_biological_synthesis_20260826/` | none | Combines DA, GSEA, ORA, EWCE and PCA over the 12 contrasts; not rebuilt on the per-compartment run. |
| Inferential checks | `03_output/inferential_checks/final_inferential_checks_20260826/` | — | mcherry learning-contrast diagnostics, post-hoc Pairing×CNO interaction, design-identifiability and collection-plate provenance audits. |
| Acquisition QC | `03_output/qc/` | — | Acquisition-level QC has not been run because `quicksearch.stats.annotated.xlsx` does not exist for this project (see caveats). |

## What the `Plate1`/`Plate2` token records

The `Plate1`/`Plate2` token in the sample identifier is the collection plate: the plate the
sample was collected onto. No technical batch metadata exist for this dataset. When a contrast
is described below as associated with collection plate, that is a statement about the design,
not an observed batch effect; wording such as "confirmed plate confound", "driven by plate" or
"technical plate effect" is not supported.

## Caveats

- **Collection plate follows Pairing.** 11 of 12 animals are paired→Plate1 / unpaired→Plate2,
  and in the `paired_veh` vs `unpaired_veh` (learning) contrast the association is complete in
  every sample class, so there a plate contribution cannot be separated from a pairing
  contribution. No proteomics batch effect from collection plate has been shown.
  `tests/test_design_balance.R` checks these counts.
- **mCherry learning contrast.** Pooled run: 1132 mapped FDR-significant proteins and no GO-BP
  pathway at q < 0.05. Per-compartment run: 998 mapped FDR proteins and coherent
  synaptic GSEA terms (GO-BP 26, GO-CC 72 at `nPermSimple = 100000`; 25 and 73 at the fgsea
  default; nearly all up). With one animal left out, effect estimates stay close (r 0.95–0.98) but the FDR count
  ranges from 1 to 1,324 and GO-BP keeps 5–25 of the 25 terms the full fit has at the fgsea
  default (GO-CC 33–72 of 73):
  `03_output/inferential_checks/within_compartment_20261124/`.
- **The Pairing×CNO interaction is post-hoc and n=3/cell.** Because collection plate and Pairing
  are associated, a Plate×Treatment contribution cannot be distinguished from a Pairing×Treatment
  one. Neither has been shown. In the per-compartment run no protein reaches FDR < 0.05 in
  any compartment; the gene-set results mostly follow a shift of the whole compartment.
- **`quicksearch.stats.annotated.xlsx` does not exist for this project** (checked 2026-08-27;
  the filename occurs only under `Exp9_Social-Stress`). The old default path
  `<Neha>/clusterProfiler/Datasets/pg_matrix/raw/quicksearch.stats.annotated.xlsx` came with the
  Exp9 code this repository was started from. The Exp9 file has the same 33 columns, but none of its
  326 rows is one of this project's 96 acquisitions (instrument `Bluto` vs `Olive`, acquired 2025-07-03 vs 2024-12-17), so it
  cannot stand in. Acquisition-level QC (MS1/MS2 signal, mass accuracy,
  normalization instability) has therefore never been run. To run it, re-export the QC report
  for this acquisition and point `PROTEOMICS_QC_QUICKSEARCH_STATS` at it.
- **`03_qc_exploration/01_qc_protein_peptide_plot.r` plots instrument QC metrics, not protein
  abundance.** Its PCA is a PCA *of acquisition QC metrics*. For protein-abundance PCA use
  `03_qc_exploration/06_pcaPlot_animal_level.r`.
- `99_historical/datasets_unmapped/` holds header-only CSVs (`gene_symbol` + newline). These are
  legitimate "zero unmapped proteins for this comparison" records, not stray files.

## Overriding defaults

Stages take paths from `getOption()`, then an environment variable, then the default
(enrichment and MapThatProt read only the environment variable). `README.md` lists the main
variables per stage. Not listed there: `PROTEOMICS_ANIMAL_LEVEL_DATA_ROOT` (step 1),
`PROTEOMICS_LIMMA_INPUT_GCT` (step 2), `PROTEOMICS_ENRICHMENT_PLOTS`, and for 05/06
`PROTEOMICS_INFERENTIAL_INPUT_GCT`, `_MAPPED_ROOT`, `_GSEA`, `_GSEA_WORKERS`, `_DESIGN_FILE` and
`_SEED_RANGE_BASES` (05 only) and `_ENRICHMENT_ROOT` (06 only).
