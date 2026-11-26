# Rank-abundance stage and the frozen source map (2026-08-28)

The shared-drive `manifests/original_figure_panel_source_map.csv` is frozen and was written
before the rank-abundance stage was reworked. It says that the nearest active script,
`03_qc_exploration/02_rank_abundance_by_sample_class.r`, draws one curve per `sample_class`
rather than per `sample_class × condition`. That is no longer true. The current script:

- reads the validated 48-column animal-level GCT;
- requires exactly one observation per `AnimalID × sample_class × condition`;
- summarizes the four sample classes across the four conditions: 16 groups of three animals;
- regenerates the Figure 3E and Supplementary D rank-abundance source data, which match the
  final reviewer-revision reference exactly when rows are joined on their key columns.

The source map is left as it is; this note records the difference. `SCRIPT_PROVENANCE.csv` maps
the executed revision scripts to their copies in this folder.
