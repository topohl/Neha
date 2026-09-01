# Publication release layer

Builds the publication / deposition package from the validated analysis: editor source
data, PRIDE metadata, provenance and checksums.

**This layer never computes science.** It reads canonical outputs and repackages them. No
model is fitted, no enrichment is run, no PCA is recomputed, no p-value is adjusted. A test
enforces that by parsing every builder and rejecting any call to `lmFit`, `eBayes`,
`gseGO`, `enrichGO`, `prcomp`, `p.adjust`, EWCE bootstrap functions and similar — see
[`tests/test_release_no_scientific_recomputation.R`](tests/test_release_no_scientific_recomputation.R).

It is kept separate from the numbered scientific stages (`01_preprocessing/` ...
`06_manuscript_figure_revision/`) so that publication-facing formatting can change without
touching anything that produced a number.

---

## Running it

Always build into scratch first:

```bash
PROTEOMICS_RELEASE_OUTPUT_ROOT="$TEMP/amp_release_$(date +%Y%m%d_%H%M%S)" \
  Rscript 07_publication_release/run_release.R

PROTEOMICS_RELEASE_OUTPUT_ROOT="<that same root>" \
  Rscript 07_publication_release/13_validate_release.R
```

With no override the release root defaults to
`S:/…/clusterProfiler/03_output/publication_release`.

Every stage is also runnable on its own with plain `Rscript`, in the order below.

## Write-side governance

`release_validate_output_root()` runs before anything is written and rejects an output root
that overlaps, in either direction:

```
01_input/            02_data/             99_historical/
_migration_backup_20260826/
03_output/{enrichment,ewce,pca,synthesis,inferential_checks,qc,
           reviewer_revision_animal_level_20260827,manuscript_curation_20260827}
```

`03_output/` is not protected wholesale, because the intended release root lives inside it;
the canonical analysis branches under it are protected by name. Containment is checked
bidirectionally through `project_paths_overlap()` from `R/project_path_utils.R` — reused
rather than reimplemented, because a divergent normalisation is how a protected root gets
missed.

A second gate covers the rest of the shared drive. The *default* release root is on the
shared drive, so a bare `Rscript run_release.R` would populate it — the same failure mode
as the 2026-08-28 incident in the repository history, where a run meant for scratch reached
the validated tree because nothing required the destination to be stated. Any release root
under the data root is therefore refused unless `PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE=true`
is set. Both gates run *before* the directory is created, so a refused target does not even
leave an empty folder behind.

The shared drive is otherwise **read-only** to this layer, including
`S:/Lab_Member/Tobi/Experiments/Collabs/Neha/` one level above `clusterProfiler/`, which
holds `sample_annotation.xlsx`, `pg.matrix_raw.txt` and `protein_count.xlsx`.

## Stages

| Stage | Produces | Mission phase |
|---|---|---|
| `01_build_sample_metadata.R` | `metadata/sample_metadata.tsv` (96), `animal_level_sample_metadata.tsv` (48), `metadata_field_provenance.tsv`, `sample_class_source_discrepancy.tsv` | 3 |
| `02_build_contrast_manifest.R` | `metadata/primary_contrast_manifest.tsv` (12), `secondary_analysis_manifest.tsv` | 4 |
| `03_build_processed_data_exports.R` | `processed_data/` abundance matrices and feature annotation | 5 |
| `04_build_differential_results.R` | `differential_analysis/` long-format statistics and summary | 6 |
| `05_build_enrichment_exports.R` | `enrichment/` GSEA, ORA, EWCE, sensitivity, coverage, parameters | 7 |
| `06_build_figure_source_data.R` | `editor_source_data/figure_source_map.tsv` and the copied panel data | 16 |
| `07_build_editor_source_workbook.R` | `editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx` | 8 |
| `08_build_editor_changelog.R` | `editor_source_data/REVISION_PROTEOMICS_DATA_CHANGELOG.md` | 9 |
| `09_build_pride_sdrf.R` | `pride/sdrf.tsv`, field status, missing-metadata and deposition notes | 10, 19 |
| `10_build_provenance.R` | `provenance/` lineage, upstream gap, versions, parameters, sessionInfo | 11, 12, 13 |
| `11_build_readme_and_dictionary.R` | `README_DATA.md`, `metadata/data_dictionary.tsv` | 15 |
| `12_build_release_manifest.R` | `provenance/release_manifest.tsv`, `SHA256SUMS.txt` | 14 |
| `13_validate_release.R` | `provenance/VALIDATION_REPORT.md`, `validation_results.tsv` | 17, 18, 19 |

`run_release.R` runs 01–12 in order, each in its own environment, and stops at the first
failure.

There is no `schemas/` directory: the column contract is generated rather than declared.
`11_build_readme_and_dictionary.R` walks the tables that were actually produced and looks
each column up in a description map, and the build **fails** if any column has no entry —
so the dictionary cannot fall behind a table that gained a field, which a hand-maintained
schema file would.

## Design decisions worth knowing

**Statistics are copied, then re-verified.** Stages 04 and 05 re-read every canonical CSV
after assembling their export and require bit-level equality. `release_format_numeric()`
writes 17 significant digits so a double round-trips exactly; `write.table()` defaults to 7
and would have broken that contract silently.

**Misleading internal field names are not propagated.** The pipeline's `gene_symbol` column
holds a UniProt entry name in `split/` and a UniProt accession in `mapped/`, and its
`Description` column holds gene symbols. The published names are `protein_group_id`,
`uniprot_accession`, `gene_symbol` and `protein_description`, each carrying what its name
says. A real protein description exists only in the search output (`pg.matrix_raw.txt`,
DIA-NN `First.Protein.Description`) and is joined in from there — 5,347 of 5,349 populated.
`release_column_is_misleading_gene_symbol()` fails the build if a published `gene_symbol`
column ever fills with accessions.

**Differential statistics come from `split/`, not `mapped/`.** `mapped/` drops the 22
non-mouse contaminant identifiers, so exporting from it would make 22 tested proteins per
comparison vanish without saying so. All 5,349 tested rows are published and
`id_mapping_status` marks which carry no mouse annotation. Both full-set and mapped-only
significance counts are reported, because the figures were drawn from the mapped set.

**Empty results are recorded as results.** 26 of the 96 enrichment analysis blocks are
legitimately empty — with no FDR-significant proteins there is no ORA query list.
`enrichment/enrichment_coverage.tsv` records every block and cross-checks its row count
against the count the canonical run recorded, so an absent block cannot be read as a missing
file.

**Nothing is guessed.** Instrument model, digestion enzyme, DIA acquisition method,
search-software version, labelling chemistry, search parameters, organism part, animal sex
and age are all absent from the project and are published as missing.
`pride/SDRF_MISSING_METADATA.md` also writes out the inferences that were *considered and
rejected*, so a later reader can see the reasoning rather than re-derive it.

**PRIDE readiness is derived, not asserted.** `release_pride_status()` is ordered
worst-first, so no optimistic branch is reachable while a pessimistic condition holds. A
missing acquisition field does not fail the scientific release; it does prevent the
deposition being called `PRIDE_READY`.

**The prohibited-call scan parses.** Regex over source text flagged a documentation string
in stage 10 that quotes the historical imputation formula verbatim. The scan now uses
`getParseData()` and counts only tokens the parser classifies as function calls.

## Overrides

| Variable | Default |
|---|---|
| `PROTEOMICS_RELEASE_OUTPUT_ROOT` | `<data root>/03_output/publication_release` |
| `PROTEOMICS_PROJECT_DATA_ROOT` | `S:/…/Collabs/Neha/clusterProfiler` |
| `PROTEOMICS_RELEASE_PROJECT_ROOT` | `S:/…/Collabs/Neha` |
| `PROTEOMICS_RELEASE_RAW_FILE_ROOT` | unset — set it to the store holding the 96 `.d` directories to re-evaluate PRIDE readiness |
| `PROTEOMICS_RELEASE_ALLOW_SHARED_DRIVE` | unset — required to be `true` before the release may be written anywhere under the shared data root |

Each also has a `getOption()` equivalent, matching the idiom every scientific stage uses.

## Tests

```bash
Rscript 07_publication_release/tests/test_release_no_scientific_recomputation.R  # no release needed
Rscript 07_publication_release/tests/test_release_sample_metadata.R
Rscript 07_publication_release/tests/test_release_contrasts.R
Rscript 07_publication_release/tests/test_release_differential_results.R
Rscript 07_publication_release/tests/test_release_enrichment.R
Rscript 07_publication_release/tests/test_release_lineage.R
Rscript 07_publication_release/tests/test_release_manifest.R
```

All but the first skip cleanly when no built release is reachable, so they are safe in CI.
Point them at a build with `PROTEOMICS_RELEASE_OUTPUT_ROOT`.

## Open item carried by the package

Preparing this package surfaced a metadata discrepancy that is **not** a consequence of the
statistical-unit correction: for 6 of the 96 acquisitions the sample class used by the
validated analysis is contradicted by both independent records of the same assignment
(`sample_annotation.xlsx`, and the autosampler plate layout the other 90 follow). It has
**not** been corrected — that would change the animal-level aggregation and every
downstream statistic — and is published as
`metadata/sample_class_source_discrepancy.tsv`, flagged per row by
`sample_class_corroborated`, and stated in `README_DATA.md` and the editor changelog. It
needs a decision from the experimenter before publication.
