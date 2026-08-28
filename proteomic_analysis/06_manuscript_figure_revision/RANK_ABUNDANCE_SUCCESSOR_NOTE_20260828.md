# Rank-abundance active-successor note — 2026-08-28

The frozen shared-drive `manifests/original_figure_panel_source_map.csv` was written before the
active rank-abundance stage was hardened in C1. Its description says that the nearest active
successor, `03_qc_exploration/02_rank_abundance_by_sample_class.r`, produces one curve per
`sample_class` rather than per `sample_class × condition`.

That historical description now predates the active implementation. The current successor:

- reads the validated 48-column animal-level GCT;
- requires one unique `AnimalID × sample_class × condition` observation per biological unit;
- summarizes the four canonical sample classes across the four canonical conditions, producing
  16 groups with three distinct animals per group; and
- regenerates the Figure 3E and Supplementary D rank-abundance source data with exact keyed
  equality to the finalized reviewer-revision reference.

The original shared-drive source map remains unchanged because it is a finalized historical
artifact. This note is the non-destructive crosswalk from that frozen description to the current
active successor. See `SCRIPT_PROVENANCE.csv` for the independent mapping between the executed
revision scripts and their byte-identical repository snapshots.
