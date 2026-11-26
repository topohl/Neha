# Manuscript figure revision (animal-level correction, 2026-08-27)

The scripts that regenerated the manuscript proteomics figures after the statistical unit was
changed from hemisphere to animal: 20 scripts, about 3,600 lines. `within_compartment_20261124/`,
`within_compartment_20261125/` and `within_compartment_20261126/` hold the later per-compartment
rebuilds (see below).

Script and panel names use the 2026-08 figure numbering: Fig 3 here is Fig. 2 in the current
manuscript, and Supp A–F are Supp. Fig. 4 A–F.

## A snapshot, not a runnable copy

These files are exact copies of the scripts that ran on the shared drive, kept here so the code
is under version control. They are not self-contained. Each script hardcodes

```r
OUT_ROOT <- "S:/.../clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))
```

so a copy run from this folder still sources `00_common.R` from the shared drive and writes
there. The runnable copy is the one in `03_output/reviewer_revision_animal_level_20260827/scripts/`.
Making these scripts self-contained (`option_or_env()` paths, relative `source()`, a
`run_all.R` driver and a `run_pipeline_check.ps1` tier) has not been done; until then, do not
edit or run them from here.

The `corrected_source_script` cells in the final manifests on the shared drive name the scripts
that ran in `03_output/reviewer_revision_animal_level_20260827/scripts/`; they are not rewritten
to point at this repository.

[`SCRIPT_PROVENANCE.csv`](SCRIPT_PROVENANCE.csv) pairs every executed script, including those of
the per-compartment rebuilds, with its copy here and records both SHA-256 values; all 63 pairs
match. The local `.gitattributes` turns off line-ending conversion for the `.R` and `.sh` files
so the hashes also hold on Windows checkouts.

[`RANK_ABUNDANCE_SUCCESSOR_NOTE_20260828.md`](RANK_ABUNDANCE_SUCCESSOR_NOTE_20260828.md) records
that the current rank-abundance stage groups by `sample_class × condition`, which the frozen
source map on the shared drive does not say.

## Old names in these snapshots

On 2026-08-28 the collaborator name was removed from the active code. The snapshots were not
edited, so their hashes still hold and they still use the old names:

| name in these snapshots | name in the active code |
|---|---|
| `R/neha_path_utils.R` | `R/project_path_utils.R` |
| `validate_neha_pca_animal_input()` | `validate_pca_animal_input()` |
| `prepare_neha_animal_pca()` | `prepare_animal_pca()` |
| `03_qc_exploration/06_pcaPlot_Neha.r` | `03_qc_exploration/06_pcaPlot_animal_level.r` |
| `NEHA_*` environment variables | `PROTEOMICS_*` |

`01_panelC_and_suppB_pca.R` loads the first three from this repository, so `R/neha_path_utils.R`
is kept as a deprecated shim that calls the new functions and warns; the shared-drive copy of
the script therefore still runs. `tests/test_deprecated_path_utils_shim.R` checks the shim and
that no active file sources it. The shim holds no analysis code and can be deleted once this
folder is made self-contained or retired.

Data filenames such as `neha_protigy_input_animal_level_primary.gct` were not renamed: they are
hash-locked validated files, so the name is kept as a dataset identifier.

## Generations

The scripts were written in three passes in 2026-08, followed by the two per-compartment rebuilds
of 2026-11 (see Per-compartment rebuilds below). Later passes supersede earlier ones for the panels they touch; all are kept because the
manifests cite specific script names as provenance.

### 1. First regeneration pass (`00`–`08`)

Established panel provenance and produced the first corrected animal-level panels.

| script | panels |
|---|---|
| `00_common.R` | shared paths, palettes, labels, save helpers |
| `01_panelC_and_suppB_pca.R` | Fig 3C, Supp B1 |
| `02_panelD_volcano.R` | Fig 3D |
| `03_panelE_suppD_rank_abundance.R` | Fig 3E, Supp D |
| `04_panelF_suppF_gsea_dotplots.R` | Fig 3F, Supp F |
| `05_panelG_retirement_audit.R` | Fig 3G (superseded; wrote the retirement note) |
| `06_panelH_scatter_panelI_heatmap.R` | Fig 3H, Fig 3I |
| `07_suppA_suppB_qc.R` | Supp A, Supp B2/B3 |
| `08_suppC_ewce_suppE_audit.R` | Supp C; Supp E (superseded; wrote the note that withheld Supp E) |

### 2. Visual-fidelity pass (`09`–`15`)

Restyled the corrected panels to match the earlier manuscript figures. Palettes were decoded from
the legend colour-bar PNGs embedded in the old SVGs; geometry and fonts were read from the SVG
attributes.

| script | panels |
|---|---|
| `09_style_match_suppC.R` | Supp C |
| `10_style_match_common.R` | recovered constants, each annotated with the file it came from |
| `11_style_match_pca.R` | Fig 3C, Supp B |
| `12_style_match_volcano.R` | Fig 3D |
| `13_style_match_rank_abundance.R` | Fig 3E, Supp D |
| `14_style_match_dotplots.R` | Fig 3F, Supp F |
| `15_style_match_scatter_heatmap.R` | Fig 3H, Fig 3I |

### 3. Retained-panel pass (`16`–`19`)

Reversed the earlier decisions to retire Fig 3G and withhold Supp E; both are kept as corrected
animal-level panels.

| script | panels |
|---|---|
| `16_final_Fig3G_lollipop.R` | Fig 3G: corrected enrichment (null result), retained |
| `17_final_Fig3H_Fig3I.R` | Fig 3H, Fig 3I: final exports and interpretation audits |
| `18_suppE_secondary_paired_model.R` | Supp E: the secondary paired model (below) |
| `19_update_manifests_final.R` | manifest status updates |

## The one new statistical model

`18_suppE_secondary_paired_model.R` is the only script here that fits a new model. Everything else
reads already-computed canonical statistics.

```
limma:  ~ 0 + sample_class + AnimalID
data:   animal-level GCT restricted to paired_veh -> 12 columns = 3 animals x 4 compartments
rank:   6 of 6 columns (full rank)   residual df: 6   eBayes df.total: 27.4
```

It is a secondary paired cross-compartment identity analysis, not a treatment effect, and is
not one of the 12 primary within-compartment contrasts. Its outputs are isolated under
`03_output/reviewer_revision_animal_level_20260827/full_regenerated/cross_compartment/SuppE_secondary_paired/`
and it writes nothing into `02_data/animal_level/`, `03_output/enrichment/`, `03_output/ewce/` or
`03_output/pca/`. All ten of its validation checks are recorded in `SuppE_validation.csv`.

## Per-compartment rebuilds (2026-11-24 to 2026-11-26)

`within_compartment_20261124/`, `within_compartment_20261125/` and `within_compartment_20261126/`
are copies of the `scripts/` folders of the matching `03_output/figure_rebuild_within_compartment_*/`
folders on the shared drive: the scripts above,
re-pointed at the rebuild folder and the per-compartment statistics, plus
`16b_final_Fig3G_gsea_lollipop.R`, which draws panel G as a GSEA lollipop in place of
`16_final_Fig3G_lollipop.R`. The 2026-11-26 rebuild is current; it differs from 2026-11-25 only
in showing gene sets with identical members once (Fig. 2 I: 66 columns instead of 70; Supp.
Fig. 4 F: 38 terms instead of 39; `identical_member_sets()` in `00_common.R`). The 2026-11-24 and
2026-11-25 rebuilds are superseded. The outputs are listed in `../CANONICAL_OUTPUTS.md`.

Only part of each folder ran:

| | scripts |
|---|---|
| run, in this order | `02`, `04`, `06`, `12`, `14`, `15`, `16b`, `17` |
| sourced by those | `00_common.R`, `10_style_match_common.R` |
| copied, not run | `01`, `03`, `05`, `07`, `08`, `09`, `11`, `13`, `16`, `18`, `19` |

Do not run the copied-only scripts. Their panels do not use the changed statistics, `16` is
replaced by `16b`, `19` would rewrite panel manifests, and all of them would write into the
rebuild folder on the shared drive. `within_compartment_20261126/run_rebuild.sh` runs the set
above. The `notes` column of `SCRIPT_PROVENANCE.csv` records, per script, whether it ran and
what changed. The copies keep the paths they ran with: `00_common.R` sets `REPO_ROOT` to the
local checkout, and `run_rebuild.sh` calls a fixed Windows Rscript path and `powershell.exe`.

## Reports and manifests

On the shared drive, under `03_output/reviewer_revision_animal_level_20260827/`:

- `manifests/original_figure_panel_source_map.csv`: where each historical panel came from
- `manifests/regenerated_panel_manifest.csv`: per-panel status and numbers
- `manifests/FINAL_REPORT.md`: the statistical correction
- `manifests/RETAINED_PANELS_3G_3H_3I_SuppE_REPORT.md`: the retained panels
- `figure_panels_style_matched/manifests/visual_fidelity_audit.csv`: style match per panel
