# Software Checklist

This document maps checklist items to repository files for the Associative Memory Proteomics analysis workflow.

| Requirement | Repository location |
| --- | --- |
| Source code | `proteomic_analysis/01_preprocessing/`, `proteomic_analysis/02_id_mapping/`, `proteomic_analysis/03_qc_exploration/`, `proteomic_analysis/04_differential_expression_enrichment/`, `proteomic_analysis/05_celltype_enrichment_EWCE/`, `proteomic_analysis/06_manuscript_figure_revision/` (figure scripts as run), `proteomic_analysis/07_publication_release/`, `proteomic_analysis/R/` (shared helpers) |
| Demo dataset | `proteomic_analysis/demo/input/` |
| README system requirements | `proteomic_analysis/README.md` |
| Dependencies and versions | `proteomic_analysis/requirements_R.md` |
| Tested OS and R versions | `proteomic_analysis/README.md`, `proteomic_analysis/requirements_R.md`; CI workflow under `.github/workflows/proteomic-analysis-demo.yml` |
| Non-standard hardware | `proteomic_analysis/README.md` notes that the demo requires none |
| Installation instructions | `proteomic_analysis/README.md` |
| Install time | `proteomic_analysis/README.md` |
| Demo instructions | `proteomic_analysis/README.md`; runnable script is `proteomic_analysis/run_demo.R` |
| Expected output | `proteomic_analysis/demo/expected_output/` |
| Demo runtime | `proteomic_analysis/README.md` |
| Running on own data | `proteomic_analysis/README.md`; configurable script is `proteomic_analysis/01_preprocessing/02_excel_convert.r` |
| Optional full reproduction instructions | `proteomic_analysis/README.md` and stage folders under `proteomic_analysis/` |
| License | `LICENSE` |
| Citation | `CITATION.cff` |
| Code repository link | `https://github.com/topohl/Neha` |

The demo needs neither the lab shared drive nor the enrichment packages.

## Data release

Built by `proteomic_analysis/07_publication_release/` (see its
[README](../07_publication_release/README.md)). Paths in the last column are relative to the
release folder.

| Requirement | Repository location | Release artifact |
| --- | --- | --- |
| Release code | `proteomic_analysis/07_publication_release/` | |
| Processed data, differential results, enrichment tables | `07_publication_release/03_build_processed_data_exports.R`, `04_build_differential_results.R`, `05_build_enrichment_exports.R` | `processed_data/`, `differential_analysis/`, `enrichment/` |
| Supplementary Table 2 | `07_publication_release/07_build_editor_source_workbook.R` | `editor_source_data/Supplementary_Table_2_proteomics.xlsx` |
| Revision changelog | `07_publication_release/08_build_editor_changelog.R` | `provenance/REVISION_PROTEOMICS_DATA_CHANGELOG.md` |
| Figure source data and panel provenance | `07_publication_release/06_build_figure_source_data.R` | `provenance/figure_source_map.tsv`, `editor_source_data/figure_source_data/` |
| PRIDE metadata | `07_publication_release/09_build_pride_sdrf.R` | `pride/README_PRIDE.md`, `pride/pride_readiness.tsv` |
| Experimenter metadata | `07_publication_release/metadata/*.tsv`, `07_publication_release/01_build_sample_metadata.R` | `metadata/experimenter_metadata.tsv`, `metadata/sample_preparation_protocol.tsv` |
| SDRF and preparation protocol | `07_publication_release/09_build_pride_sdrf.R` | `pride/sdrf.tsv`, `pride/sdrf_field_status.tsv`, `pride/SDRF_MISSING_METADATA.md`, `pride/SAMPLE_PREPARATION_PROTOCOL.md` |
| Data lineage | `07_publication_release/10_build_provenance.R` | `provenance/data_lineage.tsv`, `provenance/UPSTREAM_PREPROCESSING_GAP.md` |
| Exact software and database versions | `07_publication_release/10_build_provenance.R` | `provenance/software_versions.tsv`, `provenance/analysis_parameters.tsv`, `provenance/sessionInfo_release.txt` |
| Release SHA manifest | `07_publication_release/12_build_release_manifest.R` | `provenance/release_manifest.tsv`, `provenance/SHA256SUMS.txt` |
| Data dictionary | `07_publication_release/11_build_readme_and_dictionary.R` | `README_DATA.md`, `metadata/data_dictionary.tsv` |
| Release validation | `07_publication_release/13_validate_release.R`, `07_publication_release/tests/` | `provenance/VALIDATION_REPORT.md`, `provenance/validation_results.tsv` |

`provenance/software_versions.tsv` records exact versions per run: limma and R from the
statistics run manifest, the enrichment packages from the enrichment run audit, and the EWCE
and PCA packages from their `sessionInfo`. The search software and the MS instrument are
`UNKNOWN`, each with the place the information could be recovered from. `requirements_R.md`
lists the same versions by stage.
