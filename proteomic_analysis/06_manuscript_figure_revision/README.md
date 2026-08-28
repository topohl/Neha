# Manuscript figure revision — animal-level correction (2026-08-27)

The code that regenerated the manuscript proteomics figures after the statistical unit was
corrected from **hemisphere-level to animal-level**. 20 scripts, ~3,600 lines.

## ⚠️ This is a version-control snapshot, not a runnable copy

These files are a **verbatim** copy of what was executed. They are stored here because they were
previously the only unversioned analysis code in the project — every other stage lives in this
repository while its outputs live on the shared drive, and this set had ended up entirely on the
shared drive.

They are **not yet self-contained**. Each script hardcodes

```r
OUT_ROOT <- "S:/.../clusterProfiler/03_output/reviewer_revision_animal_level_20260827"
source(file.path(OUT_ROOT, "scripts", "00_common.R"))
```

so running a copy from this folder still sources `00_common.R` from the shared drive and still
writes there. **The runnable copy remains the one under `03_output/reviewer_revision_animal_level_20260827/scripts/`.**
Converting these to the repository's `option_or_env()` convention with relative `source()`, adding
a `run_all.R` driver and a `run_pipeline_check.ps1` tier is deferred follow-up work; until that is
done, treat this folder as read-only provenance.

The `corrected_source_script` cells in the finalized manifests intentionally remain historical
execution provenance: they identify the scripts under
`03_output/reviewer_revision_animal_level_20260827/scripts/` that actually ran. Those finalized
shared-drive manifests are not rewritten to point at Git.

[`SCRIPT_PROVENANCE.csv`](SCRIPT_PROVENANCE.csv) links every executed script to its durable
`06_manuscript_figure_revision/` repository snapshot and records both SHA-256 values. All 20 pairs
are byte-identical. The local `.gitattributes` keeps the R snapshots out of automatic line-ending
conversion so that this byte-level contract is preserved on Windows checkouts.

The dated [`RANK_ABUNDANCE_SUCCESSOR_NOTE_20260828.md`](RANK_ABUNDANCE_SUCCESSOR_NOTE_20260828.md)
documents, without changing the frozen source map, that the current active rank-abundance
successor now groups by `sample_class × condition`. Neither crosswalk makes this snapshot runnable
or part of the active pipeline.

## Name crosswalk after the 2026-08-28 de-branding

On 2026-08-28 the collaborator name was removed from the active workflow's branding and
identifiers. These snapshots were **not** edited — their SHA-256 values above still hold — so they
still reference the pre-rename names. The mapping is:

| name inside these snapshots | active name today |
|---|---|
| `R/neha_path_utils.R` | `R/project_path_utils.R` |
| `validate_neha_pca_animal_input()` | `validate_pca_animal_input()` |
| `prepare_neha_animal_pca()` | `prepare_animal_pca()` |
| `03_qc_exploration/06_pcaPlot_Neha.r` | `03_qc_exploration/06_pcaPlot_animal_level.r` |
| `NEHA_*` environment variables | `PROTEOMICS_*` |

The first three are what `01_panelC_and_suppB_pca.R` actually resolves against this repository, so
`R/neha_path_utils.R` was **kept as a deprecated shim** that delegates to the new implementation
and emits a deprecation warning. Re-running the shared-drive twin of these scripts therefore still
works. `tests/test_deprecated_path_utils_shim.R` pins that contract and also asserts that no
active file depends on the shim. The shim can be deleted once this folder is made self-contained
or formally retired; it holds no analysis logic, so removing it cannot change a result.

Data filenames such as `neha_protigy_input_animal_level_primary.gct` were deliberately **not**
renamed: they are SHA-locked validated artefacts, so the name is now treated as a legacy dataset
identifier rather than branding.

## Generations

The scripts were written in three passes. Later passes supersede earlier ones for the panels they
touch; all three are kept because the manifests cite specific script names as provenance.

### 1. First regeneration pass — `00`–`08`

Established panel provenance and produced the first corrected animal-level panels.

| script | panels |
|---|---|
| `00_common.R` | shared paths, palettes, labels, save helpers |
| `01_panelC_and_suppB_pca.R` | Fig 3C, Supp B1 |
| `02_panelD_volcano.R` | Fig 3D |
| `03_panelE_suppD_rank_abundance.R` | Fig 3E, Supp D |
| `04_panelF_suppF_gsea_dotplots.R` | Fig 3F, Supp F |
| `05_panelG_retirement_audit.R` | Fig 3G — **superseded**, produced the retirement note |
| `06_panelH_scatter_panelI_heatmap.R` | Fig 3H, Fig 3I |
| `07_suppA_suppB_qc.R` | Supp A, Supp B2/B3 |
| `08_suppC_ewce_suppE_audit.R` | Supp C; Supp E — **superseded**, produced the blocked note |

### 2. Visual-fidelity pass — `09`–`15`

Restyled the corrected panels to match the historical manuscript appearance. Aesthetics were
recovered numerically from the historical SVGs (palettes decoded from embedded legend colourbar
PNGs, geometry and fonts read from the SVG attributes) rather than approximated.

| script | panels |
|---|---|
| `09_style_match_suppC.R` | Supp C |
| `10_style_match_common.R` | recovered constants, each annotated with the file it came from |
| `11_style_match_pca.R` | Fig 3C, Supp B |
| `12_style_match_volcano.R` | Fig 3D |
| `13_style_match_rank_abundance.R` | Fig 3E, Supp D |
| `14_style_match_dotplots.R` | Fig 3F, Supp F |
| `15_style_match_scatter_heatmap.R` | Fig 3H, Fig 3I |

### 3. Retained-panel pass — `16`–`19`

Reversed the earlier decisions to retire Fig 3G and block Supp E; both are now retained as
corrected animal-level equivalents.

| script | panels |
|---|---|
| `16_final_Fig3G_lollipop.R` | Fig 3G — corrected **null** enrichment, retained |
| `17_final_Fig3H_Fig3I.R` | Fig 3H, Fig 3I — FINAL exports + interpretation audits |
| `18_suppE_secondary_paired_model.R` | Supp E — **the one authorised secondary model** |
| `19_update_manifests_final.R` | manifest status updates |

## The one new statistical model

`18_suppE_secondary_paired_model.R` is the only script here that fits a new model. Everything else
reads already-computed canonical statistics.

```
limma:  ~ 0 + sample_class + AnimalID
data:   animal-level GCT restricted to paired_veh -> 12 columns = 3 animals x 4 compartments
rank:   6 of 6 columns (full rank)   residual df: 6   eBayes df.total: 27.4
```

It is a **secondary paired cross-compartment identity analysis**, not a treatment effect, and is
not one of the 12 primary within-compartment contrasts. Its outputs are isolated under
`03_output/reviewer_revision_animal_level_20260827/full_regenerated/cross_compartment/SuppE_secondary_paired/`
and it writes nothing into `02_data/animal_level/`, `03_output/enrichment/`, `03_output/ewce/` or
`03_output/pca/`. All ten of its validation checks are recorded in `SuppE_validation.csv`.

## Reports and manifests

On the shared drive, under `03_output/reviewer_revision_animal_level_20260827/`:

- `manifests/original_figure_panel_source_map.csv` — where each historical panel came from
- `manifests/regenerated_panel_manifest.csv` — per-panel status and numbers
- `manifests/FINAL_REPORT.md` — the statistical correction
- `manifests/RETAINED_PANELS_3G_3H_3I_SuppE_REPORT.md` — the retained panels
- `figure_panels_style_matched/manifests/visual_fidelity_audit.csv` — style fidelity per panel
