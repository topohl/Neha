# Software Checklist

This document maps checklist items to repository files for the Associative Memory Proteomics analysis workflow.

| Requirement | Repository location |
| --- | --- |
| Source code | `proteomic_analysis/01_preprocessing/`, `proteomic_analysis/02_id_mapping/`, `proteomic_analysis/03_qc_exploration/`, `proteomic_analysis/04_differential_expression_enrichment/`, `proteomic_analysis/05_celltype_enrichment_EWCE/`, `proteomic_analysis/R/analysis_labels.R` |
| Demo dataset | `proteomic_analysis/demo/input/` |
| README system requirements | `proteomic_analysis/README.md` |
| Dependencies and versions | `proteomic_analysis/requirements_R.md`; exact local versions can be captured with `sessionInfo()` |
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
| Code repository link | `https://github.com/topohl/AssociativeMemoryProteomics` |

The lightweight demo is designed to run without private shared-drive paths or heavy enrichment analysis.

## Publication and data deposition

Built by `proteomic_analysis/07_publication_release/` (see its
[README](../07_publication_release/README.md)). Every item below is produced into a release
root that is validated against the protected canonical roots before anything is written.

| Requirement | Repository location | Release artefact |
| --- | --- | --- |
| Publication release code | `proteomic_analysis/07_publication_release/` | — |
| Editor source data | `07_publication_release/07_build_editor_source_workbook.R` | `editor_source_data/Proteomics_Source_Data_Animal_Level.xlsx` |
| Revision changelog for editors | `07_publication_release/08_build_editor_changelog.R` | `editor_source_data/REVISION_PROTEOMICS_DATA_CHANGELOG.md` |
| Figure source data and panel provenance | `07_publication_release/06_build_figure_source_data.R` | `editor_source_data/figure_source_map.tsv`, `editor_source_data/figure_source_data/` |
| PRIDE metadata | `07_publication_release/09_build_pride_sdrf.R` | `pride/README_PRIDE.md`, `pride/pride_readiness.tsv` |
| Canonical experimenter metadata | `07_publication_release/metadata/*.tsv`, `07_publication_release/01_build_sample_metadata.R` | `metadata/experimenter_metadata.tsv`, `metadata/sample_preparation_protocol.tsv` |
| SDRF and preparation protocol | `07_publication_release/09_build_pride_sdrf.R` | `pride/sdrf.tsv`, `pride/sdrf_field_status.tsv`, `pride/SDRF_MISSING_METADATA.md`, `pride/SAMPLE_PREPARATION_PROTOCOL.md` |
| Data lineage | `07_publication_release/10_build_provenance.R` | `provenance/data_lineage.tsv`, `provenance/UPSTREAM_PREPROCESSING_GAP.md` |
| Exact software and database versions | `07_publication_release/10_build_provenance.R` | `provenance/software_versions.tsv`, `provenance/analysis_parameters.tsv`, `provenance/sessionInfo_release.txt` |
| Release SHA manifest | `07_publication_release/12_build_release_manifest.R` | `provenance/release_manifest.tsv`, `provenance/SHA256SUMS.txt` |
| Data dictionary | `07_publication_release/11_build_readme_and_dictionary.R` | `README_DATA.md`, `metadata/data_dictionary.tsv` |
| Release validation | `07_publication_release/13_validate_release.R`, `07_publication_release/tests/` | `provenance/VALIDATION_REPORT.md`, `provenance/validation_results.tsv` |

Exact package versions are recorded per canonical run rather than as version families:
`provenance/software_versions.tsv` takes them from the enrichment run audit, the EWCE
`sessionInfo`, and the PCA `sessionInfo`, and records `UNKNOWN` with a named recovery source
where no version exists (ProTigy for the animal-level run, the search software, the
instrument). This supersedes the `1.1.x`-style families in `requirements_R.md` for
publication reporting.
