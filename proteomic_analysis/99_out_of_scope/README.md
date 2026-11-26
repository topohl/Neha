# Out of scope

Code kept here for historical reasons that does not belong to this proteomics analysis.
Nothing here is part of the numbered pipeline or run by `run_pipeline_check.ps1`, and
nothing here reads or writes this project's data. Each script refuses to run without an opt-in.

## `05_metadata_create_EXP9.r`

Belongs to Exp9_Social-Stress. It reads `TPE9_*` workbooks from the Exp9 project folder
(`Exp9_Social-Stress/proteomics/`) and writes two files back into it:

- `TPE9_samples_males_processed.tsv` (already exists there and would be overwritten)
- `TPE9_samples_males_long_with_metadata.xlsx`

It came with the Exp9 code this repository was started from, like the
`quicksearch.stats.annotated.xlsx` default described in `CANONICAL_OUTPUTS.md`. Until 2026-08-27
it sat at `01_preprocessing/05_metadata_create.r` with no guard, so running it directly would
have overwritten the Exp9 file. It was moved here, and an opt-in check and the
`PROTEOMICS_EXP9_WORK_DIR` override were added; the rest of the script is unchanged.

To run it, preferably with the output redirected to scratch:

```bash
PROTEOMICS_ALLOW_EXP9=true PROTEOMICS_EXP9_WORK_DIR=/some/scratch/dir \
  Rscript 99_out_of_scope/05_metadata_create_EXP9.r
```

Without `PROTEOMICS_EXP9_WORK_DIR` it writes to the live Exp9 folder. The R options
`proteomics.allow_exp9` and `proteomics.exp9_work_dir` do the same as the two variables.
