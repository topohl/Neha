# Publication release layer

Builds the release package from the validated analysis: Supplementary Table 2, figure source
data, PRIDE metadata, provenance and checksums.

This layer computes no statistics. It reads the analysis outputs and repackages them: no model
is fitted, no enrichment is run, no PCA is recomputed and no P value is adjusted.
[`tests/test_release_no_scientific_recomputation.R`](tests/test_release_no_scientific_recomputation.R)
parses every builder with the R parser (`getParseData()`) and fails on calls such as `lmFit`,
`eBayes`, `gseGO`, `enrichGO`, `prcomp`, `p.adjust` or the EWCE bootstrap functions; a name
inside a string or comment does not count.

It is kept apart from the numbered analysis stages (`01_preprocessing/` to
`06_manuscript_figure_revision/`) so that formatting can change without touching code that
produced a number.

## Primary statistics

Since 2026-11-24 the released statistics are the per-compartment limma run written by
`01_preprocessing/02b_limma_stat_gct.r` (`~0 + condition` within each sample class, n = 3 per
condition; the model and the reasons for it are in `../CANONICAL_OUTPUTS.md`). Two definitions
describe the run, and builders, tests and documents read them from there:

- `RELEASE_STATISTICS_RUN` (`R/release_utils.R`) says where the run's files are: the
  statistics GCT and run manifest under `02_data/animal_level/within_compartment_20261124/`
  with its `split/` and `mapped/` tables, the GO-BP and GO-CC enrichment
  `03_output/enrichment/enrichment_within_compartment_20261125/{GO_BP,GO_CC}`, the KEGG
  enrichment in the `GO_BP/` folder of `enrichment_within_compartment_20261124`, the figure
  rebuild `03_output/figure_rebuild_within_compartment_20261125` and the per-compartment
  inferential checks. `release_enrichment_permutations()` reads the permutation settings from
  the runs and returns the sentence the documents quote.
- `RELEASE_PRIMARY_MODEL` (`R/release_validation.R`) says what was fitted. The limma and R
  versions come from the run manifest through `release_read_stat_run()`, which stops
  unless the manifest matches the model and the locked hashes.

The statistics GCT (`stat_gct`) and its run manifest are locked by SHA256 in
`RELEASE_LOCKED_ARTEFACTS`, next to the animal-level input matrix and the GO-BP, GO-CC and
KEGG enrichment indexes. Stages 04 and 05 check that the split, mapped and enrichment indexes
record the locked GCT hash.

**What changed from the pooled release (2026-09-07).** That release used the ProTigy 2.4.1
export, one model over all 48 animal-level columns. Effect sizes agree to floating-point
rounding; t, P, adjusted P, B and AveExpr (now the sample-class mean) differ.
GO-CC is a primary output again (tables 2h and 2k of the workbook); the GO gene sets fgsea
could not test are table 2l (`enrichment/gsea_untested_gene_sets.tsv`). Fig. 2
D, F, G, H, I and Supp. Fig. 4 F come from the per-compartment figure rebuild, the other
panels from the 2026-08-27 revision. The pooled release stays frozen in
`03_output/publication_release`, and the pooled statistics are listed in the secondary
manifest as `superseded_pooled_model_statistics`. The changelog (section 7b) and the
manuscript actions (section 4b) give the numbers, computed from both runs' tables.

## Running it

Always build into scratch first:

```bash
PROTEOMICS_RELEASE_OUTPUT_ROOT="$TEMP/amp_release_$(date +%Y%m%d_%H%M%S)" \
  Rscript 07_publication_release/run_release.R

PROTEOMICS_RELEASE_OUTPUT_ROOT="<that same root>" \
  Rscript 07_publication_release/13_validate_release.R
```

With no override the release root defaults to
`S:/…/clusterProfiler/03_output/release_within_compartment_20261126`. The earlier
releases (`03_output/publication_release`, pooled model; `03_output/release_within_compartment_20261125`,
old workbook layout) are frozen and protected. To publish,
commit first, then rebuild and validate on the shared drive with
`PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE=true`; never copy a scratch build there.

Every stage is also runnable on its own with plain `Rscript`, in the order below.

## Where it may write

`release_validate_output_root()` runs before anything is written and refuses an output root
that overlaps any of these, in either direction:

```
01_input/            02_data/             99_historical/
_migration_backup_20260826/
03_output/{enrichment,ewce,pca,synthesis,inferential_checks,qc,
           reviewer_revision_animal_level_20260827,manuscript_curation_20260827,
           figure_rebuild_within_compartment_20261124, figure_rebuild_within_compartment_20261125,
           inferential_checks/within_compartment_20261124,
           publication_release,                    (the frozen pooled-model release)
           release_within_compartment_20261125}    (the frozen first per-compartment release)
```

`03_output/` is not protected as a whole because the release root lives inside it; the
analysis folders under it are listed by name. Overlap is tested with `project_paths_overlap()`
from `R/project_path_utils.R`, the helper the analysis stages also use.

Any release root under the data root is also refused unless
`PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE=true` is set, because the default root is on the shared
drive and a bare `Rscript run_release.R` would otherwise write there. Both checks run before the
folder is created, so a refused run leaves no empty folder behind.

The rest of the shared drive is read-only to this layer, including
`S:/Lab_Member/Tobi/Experiments/Collabs/Neha/` one level above `clusterProfiler/`, which
holds `sample_annotation.xlsx`, `pg.matrix_raw.txt` and `protein_count.xlsx`.

## Machine and display names

Each table ships twice:

| Artifact | Names | For |
|---|---|---|
| `*.tsv`, `*.tsv.gz` | machine (`adj.P.Val`, `canonical_comparison`, `mcherry`) | scripts and the data dictionary |
| `Supplementary_Table_2_proteomics.xlsx` | display (`Adjusted P value (BH)`, `Comparison`, `mCherry+`) | readers of the paper |

The workbook is Supplementary Table 2: a Contents sheet, a README, Column definitions and the
sub-tables 2a–2q (samples, animal profiles, comparisons, significant proteins, all proteins,
significant gene sets, GSEA, ORA, gene sets without a P value, enrichment settings, EWCE,
secondary analyses, figure panels, software). It holds what a reader needs; the complete
tables (the sensitivity GSEA, the alternative ORA list, gene members as identifiers, the
metadata field status) stay in the release files that the README sheet names.

`R/release_presentation.R` holds the name mapping and the workbook layout rules (no data, no
statistics); `R/release_dictionary.R` holds the column definitions shared by the data
dictionary and the Column definitions sheet:

- every sheet has its title on row 1, the caption on row 2 and the header on row 4; row 3
  ("Applies to every row") holds the values that are the same on every row, either from a
  constant column that the sheet does not show (the `.tsv` keeps it) or, for 2m, from the run
  parameters; the inventory lists them under `table_level_constants`;
- cells show display labels (`mCherry+: paired VEH vs unpaired VEH`); the machine key is
  kept as the last column, `Comparison key`, so the sheets can still be joined to the `.tsv`
  files;
- each sheet lists the columns it shows and the ones it drops, and stage 07 stops on a column
  that is in neither list or has no definition;
- number formats follow the values in each column: P values below 0.001 in scientific
  notation, the rest as decimals; GSEA P values at fgsea's lower bound show as `<=1E-10` and
  EWCE P values of 0 as `<0.0001`.

The sub-table letters are set in `RELEASE_SHEET_META` and pinned for the checks in
`RELEASE_WORKBOOK_LAYOUT` (`R/release_validation.R`), so they cannot change by accident.
Display headers must be unique and non-blank, like the machine names, and contain no
underscores.

Prose and Markdown tables use display names ("c-Fos+: paired CNO vs paired VEH", not
`cfos_paired_cno_over_cfos_paired_veh`); the `.tsv` files keep the controlled vocabulary
because the tables are joined on it. `README_DATA.md` lists both.

## Stages

| Stage | Produces |
|---|---|
| `01_build_sample_metadata.R` | `metadata/sample_metadata.tsv` (96), `animal_level_sample_metadata.tsv` (48), `metadata_field_provenance.tsv`, `sample_class_corrections.tsv` (6), checked-in experimenter metadata and protocol records |
| `02_build_contrast_manifest.R` | `metadata/primary_contrast_manifest.tsv` (12), `secondary_analysis_manifest.tsv` |
| `03_build_processed_data_exports.R` | `processed_data/` abundance matrices and feature annotation |
| `04_build_differential_results.R` | `differential_analysis/` long-format statistics, summary and `significant_proteins.tsv` |
| `05_build_enrichment_exports.R` | `enrichment/` GSEA, ORA, EWCE, sensitivity, coverage, parameters, `significant_gene_sets.tsv.gz` and `gsea_untested_gene_sets.tsv` |
| `06_build_figure_source_data.R` | `editor_source_data/figure_source_data/` (one or two CSVs per panel, with plain headers; the canonical panel files are checked against pinned hashes, not copied) and `provenance/figure_source_map.tsv` |
| `07_build_editor_source_workbook.R` | `editor_source_data/Supplementary_Table_2_proteomics.xlsx`, `workbook_sheet_inventory.tsv` |
| `08_build_editor_changelog.R` | `provenance/REVISION_PROTEOMICS_DATA_CHANGELOG.md`, `EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv`, `MANUSCRIPT_TERMINOLOGY_ACTIONS.md` |
| `09_build_pride_sdrf.R` | `pride/sdrf.tsv`, field status, sample-preparation protocol, missing-metadata and deposition notes |
| `10_build_provenance.R` | `provenance/` lineage, upstream gap, versions, parameters, sessionInfo |
| `11_build_readme_and_dictionary.R` | `README_DATA.md`, `metadata/data_dictionary.tsv` |
| `12_build_release_manifest.R` | `provenance/release_manifest.tsv`, `SHA256SUMS.txt` |
| `13_validate_release.R` | `provenance/VALIDATION_REPORT.md`, `validation_results.tsv` |

`run_release.R` runs 01–12 in order, each in its own environment, and stops at the first
failure.

The data dictionary is built from the tables actually produced:
`11_build_readme_and_dictionary.R` looks up each column in a description map and stops if a
column has no entry, so a new field cannot go undocumented.

## Design notes

### Repository paths

Repository files are named relative to `proteomic_analysis/`. Lineage, parameter,
field-provenance and manifest tables write repository sources as relative paths (for example
`07_publication_release/metadata/experimenter_metadata.tsv`), so no local folder is released;
shared-drive data paths are written in full. `release_repo_relative()` in `R/release_utils.R`
does this, and stage 13 and the tests check that every such path resolves in the repository.

### Copied statistics

Statistics are copied, then checked. Stages 04 and 05 re-read every source CSV after
building their export and require the values to be identical. `release_format_numeric()` writes
17 significant digits so a double survives the round trip (`write.table()` writes 15, which does not always round-trip).

### Column names

Internal column names are not carried over. The pipeline's `gene_symbol` column holds a
UniProt entry name in `split/` and a UniProt accession in `mapped/`, and its `Description`
column holds gene symbols. The release uses `protein_group_id`, `uniprot_accession`,
`gene_symbol` and `protein_description`, each holding what its name says. Protein descriptions
are joined from the search output (`First.Protein.Description` in `pg.matrix_raw.txt`; 5,347 of
5,349 filled). `release_column_is_misleading_gene_symbol()` stops the build if a released
`gene_symbol` column fills with accessions.

### Source of the differential statistics

Differential statistics come from `split/`, not `mapped/`. `mapped/` drops the 22 non-mouse
contaminant identifiers, so exporting from it would lose 22 tested proteins per comparison. All
5,349 tested rows are released, and `id_mapping_status` marks those without mouse annotation.
Significance counts are given for the full set and for the mapped set, which the figures used.

### Empty results

Empty results are recorded. Many enrichment blocks are empty (with no FDR-significant
proteins there is no ORA query list). `enrichment/enrichment_coverage.tsv` lists every block and
checks its row count against the count the analysis run recorded, so an empty block is not
mistaken for a missing file.

### Metadata

Metadata are not inferred. Facts come from the checked-in records in
`07_publication_release/metadata/`; what is still missing is listed under the open item below.
`pride/SDRF_MISSING_METADATA.md` also lists the inferences that were considered and rejected,
with the reasons.

### PRIDE status

PRIDE status is computed. `release_pride_status()` tests the worst condition first, so the
status cannot read better than the inputs allow. A missing acquisition field does not fail the
release, but it keeps the deposition from being called `PRIDE_READY`.

## Overrides

| Variable | Default |
|---|---|
| `PROTEOMICS_RELEASE_OUTPUT_ROOT` | `<data root>/03_output/release_within_compartment_20261126` |
| `PROTEOMICS_PROJECT_DATA_ROOT` | `S:/…/Collabs/Neha/clusterProfiler` |
| `PROTEOMICS_RELEASE_PROJECT_ROOT` | `S:/…/Collabs/Neha` |
| `PROTEOMICS_RELEASE_RAW_FILE_ROOT` | unset; set it to the store holding the 96 `.d` directories to re-evaluate PRIDE readiness |
| `PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE` | unset; must be `true` to write anywhere under the shared data root |
| `PROTEOMICS_RELEASE_ALLOW_DIRTY` | unset; a shared-drive release is refused while the repository has uncommitted or untracked changes, and `true` writes it anyway with the commit recorded as `<hash>+dirty` |
| `PROTEOMICS_RELEASE_OLD_PACKAGE_ROOT` | unset; point it at the extracted original submission package to upgrade the crosswalk from `DIMENSIONS_AND_CONTENT_LINEAGE` to `DIRECT_VERIFICATION` with real hashes |

Each can also be set as an R option (`getOption()`), as in the analysis stages.

## Tests

```bash
Rscript 07_publication_release/tests/test_release_no_scientific_recomputation.R  # no release needed
Rscript 07_publication_release/tests/test_release_sample_metadata.R
Rscript 07_publication_release/tests/test_release_contrasts.R
Rscript 07_publication_release/tests/test_release_differential_results.R
Rscript 07_publication_release/tests/test_release_enrichment.R
Rscript 07_publication_release/tests/test_release_lineage.R
Rscript 07_publication_release/tests/test_release_manifest.R
Rscript 07_publication_release/tests/test_release_sample_class_correction.R
Rscript 07_publication_release/tests/test_release_effect_size_semantics.R
Rscript 07_publication_release/tests/test_release_pride_experimenter_metadata.R
```

All but the first skip when no built release is reachable (`test_release_pride_experimenter_metadata.R` still checks the checked-in metadata). Point them at a build with
`PROTEOMICS_RELEASE_OUTPUT_ROOT`.

`../run_pipeline_check.ps1` runs them at every tier. Without an override they check the release
at the default root on the shared drive. They skip until a release is deployed there; after that,
a FAIL can mean the deployed release is older than the code rather than a fault in the code.

## Two reporting corrections

**The six C46/C47 sample classes are a recorded correction, not an open discrepancy.** For 6 of
the 96 acquisitions the sample class used in the analysis differs from every earlier record. An
audit found that this was an intentional correction made during sample-identity QC: the left
hemispheres of C46 and C47 were reassigned in a cycle, `neuropil -> mcherry -> neuron ->
neuropil`. Only the labels changed. The acquisitions and abundance profiles are the same, and
the animal-level units were built from the corrected labels, so nothing needs to be reanalyzed.

`RELEASE_SAMPLE_CLASS_CORRECTION` in [`R/release_validation.R`](R/release_validation.R) defines
the correction. Stage 01 re-derives it from three earlier records (the `group_label` column of
`sample_annotation.xlsx`, the autosampler plate layout and the retained UMAP correction table),
requires all three to agree, checks the correction table against the SHA256 the audit recorded,
and requires the result to match the definition. Both labels are released for all 96
measurements (`original_sample_class`, `analysis_sample_class`, `sample_class_corrected`), with
per-row provenance in `metadata/sample_class_corrections.tsv`.

No written note of the reason survives; the correction table is the record, and the release says
so. UMAP is not treated as ground truth: its nearest-class-center suggestion agrees with the
applied class for five of the six and differs for N60 (suggested `cfos`, applied `neuron`), so
the suggestion is released in `umap_nearest_class_suggestion` next to the applied class.

**The effect size is not a log2 fold change.** The analysis matrix is standardized separately
for each protein across the measurement-level dataset, so the coefficient the statistics GCT
stores as `logFC` (limma's name) is a standardized abundance difference in SD units, not a log2
ratio. `RELEASE_EFFECT_SIZE` holds the definition that the builders and the validator read.

Two phrasings are avoided: "standardized mean difference" on its own, which reads as Cohen's d
(a pooled within-group SD, not the per-protein SD across the dataset used here); and "z-scored",
because the released values have two decimals, no row is exactly mean 0 / SD 1, and the step
that produced the scaling is UNRESOLVED in the lineage. `release_standardization_evidence()`
measures the matrix and returns `exact_zscore`, and the builders choose their wording from that.

Internal names are kept so the link to the analysis run holds:
`enrichment/GSEA_log2FC_sensitivity.tsv.gz` keeps its filename, `rank_statistic` keeps the value
`log2fc`, and every differential row has `source_statistic_field = "logFC"`. The exported
`effect_size_sd_units` is identical to the source `logFC` for all 64,188 rows, checked at build
time and again by
[`tests/test_release_effect_size_semantics.R`](tests/test_release_effect_size_semantics.R).
Every fold-change mention in the released documents is classified in
`provenance/EFFECT_SIZE_TERMINOLOGY_AUDIT.tsv`; what has to be fixed by hand outside this
repository (figure axes, captions, manuscript text) is listed in
`provenance/MANUSCRIPT_TERMINOLOGY_ACTIONS.md`. That scan reads the built documents, not
the builder sources, which necessarily contain the terms it classifies.

## Open item

The PRIDE deposition is `PRIDE_METADATA_INCOMPLETE` because metadata are missing, not because of
formatting. Instrument model, acquisition mode and developmental stage are required and
unresolved; search-software version, modification parameters, mass tolerances, AnimalID-level
sex and age at tissue collection are also unknown. The 96 `.d` acquisition files are known by
name but held elsewhere. `pride/SDRF_MISSING_METADATA.md` lists what is missing and what would
supply it. The checked-in records in `metadata/` resolve CeM, the C57BL/6J cohort, the
marker-defined LCM classes, label-free status and the 2024-11-28 Lys-C/trypsin preparation
protocol.
