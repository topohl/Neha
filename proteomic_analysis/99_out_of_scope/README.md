# Out of scope

Code that lives in this repository for historical reasons but does **not** belong to the Neha
proteomics analysis. Nothing here is part of the numbered pipeline, nothing here is invoked by
`run_pipeline_check.ps1`, and nothing here reads or writes Neha data.

Everything in this folder must refuse to run unless it is explicitly opted into.

## `05_metadata_create_EXP9.r`

Belongs to **Exp9_Social-Stress**, not Neha. It reads `TPE9_*` workbooks from the Exp9 project
folder and writes two files back into it:

- `TPE9_samples_males_processed.tsv`
- `TPE9_samples_males_long_with_metadata.xlsx`

It was inherited when this repository was seeded from Exp9 — the same origin as the
`quicksearch.stats.annotated.xlsx` default described in `CANONICAL_OUTPUTS.md`.

### Why it moved

It previously sat at `01_preprocessing/05_metadata_create.r`, inside the numbered active stage
sequence, and it was the only active-stage script with no `option_or_env()` override, no write
guard, and no sourcing of `R/*_utils.R`. On the shared drive every precondition for an accidental
cross-project overwrite was satisfied:

| precondition | state |
|---|---|
| `Exp9_Social-Stress/proteomics/` exists | yes |
| `TPE9_sample_metadata_males.xlsx` (read) | exists |
| `TPE9_samples_males.xlsx` (glob) | exists |
| `TPE9_samples_males_processed.tsv` (write) | **already exists → would be overwritten** |

`run_pipeline_check.ps1` already recorded the stage as `SKIP` with the reason *"belongs to a
different project (Exp9_Social-Stress)"*, so the runner never executed it. The gap was direct
invocation: `Rscript 01_preprocessing/05_metadata_create.r` — the pattern the README documents for
running a single stage — ran to completion with no warning.

### Running it deliberately

```bash
NEHA_ALLOW_EXP9=true Rscript 99_out_of_scope/05_metadata_create_EXP9.r
```

The target directory can be redirected away from the live Exp9 folder, which is strongly
recommended for any test run:

```bash
NEHA_ALLOW_EXP9=true NEHA_EXP9_WORK_DIR=/some/scratch/dir \
  Rscript 99_out_of_scope/05_metadata_create_EXP9.r
```

Both settings are also available as R options: `neha.allow_exp9`, `neha.exp9_work_dir`.

The analysis body of the script is unchanged from the version that sat in `01_preprocessing/`;
only the banner, the guard and the two path overrides were added.
