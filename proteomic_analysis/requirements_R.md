# R Requirements

## Tested environment

- Windows 11 with R 4.5.1
- Ubuntu GitHub Actions runner with R 4.5.1

`run_demo.R` and the tests run in CI need only base R.

## Packages by stage

Versions are the ones recorded on the shared drive: the statistics run manifest (2026-11-24), the package audit of the GO enrichment run `enrichment_within_compartment_20261125`, and the `sessionInfo()` of the EWCE and PCA runs (2026-08-25). Where a package appears in more than one record the version is the same. Packages without a version were not recorded. The data release collects these versions in `provenance/software_versions.tsv`.

| Stage | Packages |
| --- | --- |
| `01_preprocessing/02_excel_convert.r` | readxl, dplyr, writexl |
| `01_preprocessing/02a_prepare_animal_level_protigy_input.r` | readxl, writexl, digest 0.6.39 |
| `01_preprocessing/02b_limma_stat_gct.r` | Bioconductor: limma 3.66.0; CRAN: digest 0.6.39 |
| `01_preprocessing/03_gct_extractR.r` | digest 0.6.39 |
| `02_id_mapping/01_MapThatProt_batch.r` | pacman, dplyr, stringr, tidyr, purrr, readr, tibble, openxlsx, R.utils, foreach, doParallel, readxl; Bioconductor: AnnotationDbi 1.72.0, org.Mm.eg.db 3.22.0, UniProt.ws |
| `04_differential_expression_enrichment/` | Bioconductor: clusterProfiler 4.18.4, DOSE 4.4.0, enrichplot 1.30.5, fgsea 1.36.2, GOSemSim 2.36.0, GO.db 3.22.0, KEGGREST 1.50.0, BiocParallel 1.44.0, org.Mm.eg.db 3.22.0, AnnotationDbi 1.72.0; CRAN: ggplot2 4.0.2, withr 3.0.2, digest 0.6.39. `05_pairing_cno_interaction.r` and `06_mcherry_learning_robustness.r` also need limma. |
| `03_qc_exploration/06_pcaPlot_animal_level.r` | pacman, data.table 1.18.2.1, ggplot2 4.0.2, factoextra 2.0.0, reshape2 1.4.5, ggrepel 0.9.8, uwot 0.2.4, RColorBrewer 1.1-3, pheatmap 1.0.13, Rtsne 0.17, aricode 1.1.0, rospca 1.1.1, irlba 2.3.7, pandoc 0.2.0, treemapify 2.6.0, dplyr 1.2.1, digest 0.6.39. Some optional plots use igraph, ggraph, GGally, plotly, circlize, ggalluvial, ggExtra, gridExtra, cluster, treemap or d3r and are skipped when these are not installed. |
| `03_qc_exploration/02_rank_abundance_by_sample_class.r` | readxl, tidyverse, scales, ggrepel, writexl, svglite, digest |
| `05_celltype_enrichment_EWCE/01_EWCE.r` | Bioconductor: EWCE 1.18.1, ewceData 1.18.0, limma 3.66.0, org.Mm.eg.db 3.22.0, AnnotationDbi 1.72.0; CRAN: dplyr 1.2.1, tibble 3.3.1, tidyr 1.3.2, ggplot2 4.0.2, pheatmap 1.0.13, svglite 2.2.2, ggridges 0.5.7, ggrepel 0.9.8, ggsci 5.0.0, viridis 0.6.5, openxlsx 4.2.8.1, stringr 1.6.0, patchwork 1.3.2, future 1.70.0, future.apply 1.20.2, digest 0.6.39, BiocManager 1.30.27 |
| `07_publication_release/` | digest, readxl, openxlsx; fgsea 1.36.2 if installed (used only to quote its default `nPermSimple` in the documents; without it the count is left out) |

Most stages stop with the names of any missing packages. `01_EWCE.r`, `01_preprocessing/04_format_metadata.r` and the ID-mapping and PCA stages (`pacman::p_load()` and `install.packages()`) install missing packages instead.
