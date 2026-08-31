<#
    Associative Memory Proteomics pipeline smoke test.

    Verifies the pipeline still works after the 2026-08-26 folder restructure, WITHOUT
    touching any validated output: every producing stage is redirected to a scratch root
    via its documented env-var override.

    Tiers:
      -Tier Fast    tests + data-integrity gates only            (~2 min)
      -Tier Normal  Fast + animal-level GCT rebuild + split + PCA (~10-20 min)   [default]
      -Tier Full    Normal + ID mapping + enrichment + EWCE       (hours; EWCE is 10k reps)

    Usage:
      powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1
      powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1 -Tier Fast
      powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1 -Tier Full

    NOTE ON COVERAGE: the chain is
        02a_prepare_animal_level -> [ProTigy, EXTERNAL TOOL] -> 03_gct_extractR -> ...
    ProTigy is a separate GUI/tool, so no script can run the chain truly end to end. This
    runner covers both sides of that gap using the existing committed ProTigy output.
#>

param(
    [ValidateSet("Fast", "Normal", "Full")]
    [string]$Tier = "Normal"
)

$ErrorActionPreference = "Continue"

$repo    = $PSScriptRoot
$rscript = "C:\Users\topohl\AppData\Local\Programs\R\R-4.5.1\bin\x64\Rscript.exe"
$dataRoot = "S:\Lab_Member\Tobi\Experiments\Collabs\Neha\clusterProfiler"
$scratch = Join-Path $env:TEMP ("proteomics_pipeline_check_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Force -Path $scratch | Out-Null

# --------------------------------------------------------------------------------
# Fail-closed safety gates (added 2026-08-31 after the 2026-08-28 16:36 incident).
# See run_pipeline_safety.ps1 for the full incident write-up. In short: the runner was
# loaded from one commit, the tree was checked out to another mid-run, the stage scripts
# then expected the pre-rename NEHA_* variables, and 37 files were written into the
# validated shared-drive tree instead of scratch. Three gates now stand in front of every
# producing stage: repository identity, override contract, and scratch-root containment.
# --------------------------------------------------------------------------------
. (Join-Path $PSScriptRoot "run_pipeline_safety.ps1")

# Pin repository identity at startup. Every producing stage re-checks against this.
$startupRepoState = Get-GitRepoState -RepoRoot $repo
Write-Host ""
Write-Host "repository pinned at startup:" -ForegroundColor Cyan
Write-Host ("  root   : {0}" -f $startupRepoState.Root)
Write-Host ("  HEAD   : {0}" -f $startupRepoState.Head)
Write-Host ("  branch : {0}" -f $startupRepoState.Branch)

# --------------------------------------------------------------------------------
# Neutralise inherited PROTEOMICS_* / NEHA_* overrides.
# Every stage resolves paths as option -> env var -> default. A stale variable in the
# calling shell is therefore silently honoured by the child Rscript processes, which makes
# results depend on session history instead of the repo. The pre-rename NEHA_* names are
# swept too: a leftover legacy variable could still steer an older stage script at a
# canonical path. Cleared up front, restored on exit, and reported.
# --------------------------------------------------------------------------------
$inheritedOverride = Clear-InheritedOverrides -Prefixes @("PROTEOMICS_", "NEHA_")
if ($inheritedOverride.Count -gt 0) {
    Write-Host ""
    Write-Host "WARNING: cleared inherited override(s) from your shell for this run:" -ForegroundColor Yellow
    foreach ($k in $inheritedOverride.Keys) {
        Write-Host ("  {0} = {1}" -f $k, $inheritedOverride[$k]) -ForegroundColor Yellow
    }
    Write-Host "  (restored when this script exits; set them deliberately if you meant to override)" -ForegroundColor Yellow
    Write-Host ""
}

$results = New-Object System.Collections.ArrayList
function Record($name, $status, $detail, $seconds) {
    [void]$results.Add([pscustomobject]@{
        Stage = $name; Status = $status; Seconds = [math]::Round($seconds, 1); Detail = $detail
    })
    $color = "White"
    if ($status -eq "PASS") { $color = "Green" }
    if ($status -eq "FAIL") { $color = "Red" }
    if ($status -eq "SKIP") { $color = "Yellow" }
    if ($status -eq "ABORT") { $color = "Magenta" }
    Write-Host ("  [{0,-5}] {1,-38} {2,7:N1}s  {3}" -f $status, $name, $seconds, $detail) -ForegroundColor $color
}

$script:safetyAborts = 0

function Invoke-Stage($name, $scriptRelPath, $envMap, $requiredInputs) {
    foreach ($p in $requiredInputs) {
        if (-not (Test-Path -LiteralPath $p)) {
            Record $name "SKIP" ("missing input: " + (Split-Path $p -Leaf)) 0
            return
        }
    }

    # ---- fail-closed safety gates: all of these run BEFORE Rscript is invoked ----------
    # Any violation records ABORT and returns without touching R, so a drifted tree or an
    # unrecognised override can never reach a producing stage.
    try {
        # A. repository identity has not changed since startup
        $currentRepoState = Get-GitRepoState -RepoRoot $repo
        [void](Assert-NoCheckoutDrift -StartupState $startupRepoState -CurrentState $currentRepoState -StageName $name)

        # A. the stage script still resolves where we expect, inside that repository
        $stageFull = Assert-StageScriptResolves -RepoRoot $repo -StageRelPath $scriptRelPath -StageName $name

        # B. the code about to run actually recognises every override we are handing it
        [void](Assert-OverrideContract -RepoRoot $repo -StageFullPath $stageFull -EnvNames @($envMap.Keys) -StageName $name)

        # C. every destination we direct lands inside this run's scratch root, never in a
        #    protected canonical root
        [void](Assert-SafeOutputTargets -EnvMap $envMap -ScratchRoot $scratch -StageName $name)
    } catch {
        $script:safetyAborts++
        Record $name "ABORT" $_.Exception.Message 0
        return
    }

    foreach ($k in $envMap.Keys) { Set-Item -Path ("Env:" + $k) -Value $envMap[$k] }

    # D. nothing else in the environment can redirect this stage
    try {
        [void](Assert-NoStaleOverridesPresent -EnvMap $envMap -StageName $name)
    } catch {
        foreach ($k in $envMap.Keys) { Remove-Item -Path ("Env:" + $k) -ErrorAction SilentlyContinue }
        $script:safetyAborts++
        Record $name "ABORT" $_.Exception.Message 0
        return
    }

    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    $log = Join-Path $scratch ((Split-Path $scriptRelPath -Leaf) + ".log")
    & $rscript $stageFull *> $log
    $code = $LASTEXITCODE
    $sw.Stop()
    foreach ($k in $envMap.Keys) { Remove-Item -Path ("Env:" + $k) -ErrorAction SilentlyContinue }
    if ($code -eq 0) {
        Record $name "PASS" "ok" $sw.Elapsed.TotalSeconds
    } else {
        $tail = (Get-Content -LiteralPath $log -Tail 3 -ErrorAction SilentlyContinue) -join " | "
        Record $name "FAIL" $tail $sw.Elapsed.TotalSeconds
    }
}

Write-Host ""
Write-Host "Proteomics pipeline smoke test  (tier: $Tier)" -ForegroundColor Cyan
Write-Host "scratch output: $scratch"
Write-Host ""

# ---------------------------------------------------------------- 1. contracts
Write-Host "1. Contract / unit tests" -ForegroundColor Cyan
# The tests are read-only, but a tree that has already drifted makes every result below
# meaningless, so check here too rather than only in front of the producing stages.
try {
    [void](Assert-NoCheckoutDrift -StartupState $startupRepoState `
                                  -CurrentState (Get-GitRepoState -RepoRoot $repo) `
                                  -StageName "contract tests")
} catch {
    Write-Host $_.Exception.Message -ForegroundColor Magenta
    exit 1
}
# The safety net gates every producing stage below, so verify it before relying on it.
# Fixture-only: touches no shared drive and invokes no R stage.
$safetyTest = Join-Path $repo "tests\test_runner_safety.ps1"
if (Test-Path -LiteralPath $safetyTest) {
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    $log = Join-Path $scratch "test_runner_safety.log"
    & powershell -NoProfile -ExecutionPolicy Bypass -File $safetyTest *> $log
    $code = $LASTEXITCODE
    $sw.Stop()
    if ($code -eq 0) {
        Record "test_runner_safety (ps1)" "PASS" "ok" $sw.Elapsed.TotalSeconds
    } else {
        $tail = (Get-Content -LiteralPath $log -Tail 3 -ErrorAction SilentlyContinue) -join " | "
        Record "test_runner_safety (ps1)" "FAIL" $tail $sw.Elapsed.TotalSeconds
    }
} else {
    Record "test_runner_safety (ps1)" "FAIL" "safety contracts missing from tests/" 0
}

Get-ChildItem -LiteralPath (Join-Path $repo "tests") -Filter "*.R" | ForEach-Object {
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    $log = Join-Path $scratch ($_.BaseName + ".log")
    & $rscript $_.FullName *> $log
    $code = $LASTEXITCODE
    $sw.Stop()
    if ($code -eq 0) {
        Record $_.BaseName "PASS" "ok" $sw.Elapsed.TotalSeconds
    } else {
        $tail = (Get-Content -LiteralPath $log -Tail 3 -ErrorAction SilentlyContinue) -join " | "
        Record $_.BaseName "FAIL" $tail $sw.Elapsed.TotalSeconds
    }
}

# Data integrity is covered by tests/test_data_integrity.R, which the loop above already ran
# (validated GCT SHA256 + every recorded path/hash pair + the production index reader).

if ($Tier -eq "Fast") {
    Write-Host ""
    Write-Host "Tier=Fast: stopping before producing stages." -ForegroundColor Yellow
}

# ------------------------------------------------------- 3. producing stages
if ($Tier -ne "Fast") {
    Write-Host ""
    Write-Host "3. Pipeline stages (all output redirected to scratch)" -ForegroundColor Cyan

    # 3a. rebuild the animal-level GCT and check it reproduces the locked hash
    $gctOut = Join-Path $scratch "animal_level_gct"
    Invoke-Stage "02a_prepare_animal_level_gct" "01_preprocessing\02a_prepare_animal_level_protigy_input.r" `
        @{ PROTEOMICS_ANIMAL_LEVEL_OUTPUT_DIR = $gctOut } `
        @( (Join-Path $dataRoot "02_data\gct\pg.matrix_filtered_pcaAdjusted_unnormalized.gct"),
           (Join-Path $dataRoot "02_data\gct\imputed_data.xlsx"),
           (Join-Path $dataRoot "01_input\metadata\sample_info.xlsx") )

    $rebuilt = Join-Path $gctOut "neha_protigy_input_animal_level_primary.gct"
    if (Test-Path -LiteralPath $rebuilt) {
        $h = (Get-FileHash -LiteralPath $rebuilt -Algorithm SHA256).Hash.ToLower()
        if ($h -eq "f12cf99e1bfb7c17bbf56bffb6783e924698bce5d5533a8e312bc4bbb733bbb3") {
            Record "02a_reproduces_locked_hash" "PASS" "bit-identical to validated GCT" 0
        } else {
            Record "02a_reproduces_locked_hash" "FAIL" ("hash differs: " + $h.Substring(0,16) + "...") 0
        }
    }

    # 3b. ProTigy is external -- flag the gap explicitly
    Record "ProTigy_stat_results" "SKIP" "EXTERNAL TOOL - cannot be scripted; using committed output" 0

    # 3c. split the committed ProTigy stat GCT
    $splitOut = Join-Path $scratch "split"
    Invoke-Stage "03_gct_extractR (split)" "01_preprocessing\03_gct_extractR.r" `
        @{ PROTEOMICS_PROTIGY_STAT_GCT_OUTPUT_ROOT = $splitOut } `
        @( (Join-Path $dataRoot "02_data\animal_level\stat_results_for_ssGSEA_neha_proteome.gct") )

    # 3d. animal-level PCA (independent branch off the validated GCT)
    Invoke-Stage "06_pcaPlot_animal_level (PCA workflow)" "03_qc_exploration\06_pcaPlot_animal_level.r" `
        @{ PROTEOMICS_PCA_OUTPUT_ROOT = (Join-Path $scratch "pca") } `
        @( (Join-Path $dataRoot "02_data\animal_level\input_gct\neha_protigy_input_animal_level_primary.gct") )

    # 3e. animal-level rank-abundance QC (regenerates the Fig 3E / Supp D panels)
    Invoke-Stage "02_rank_abundance (animal level)" "03_qc_exploration\02_rank_abundance_by_sample_class.r" `
        @{ PROTEOMICS_RANK_ABUNDANCE_OUTPUT_DIR = (Join-Path $scratch "rank_abundance") } `
        @( (Join-Path $dataRoot "02_data\animal_level\input_gct\neha_protigy_input_animal_level_primary.gct"),
           (Join-Path $dataRoot "02_data\animal_level\mapped\forward\mcherry_paired_veh_vs_mcherry_unpaired_veh.csv") )

    # 3f. stages whose inputs are genuinely absent from the project tree
    Record "01_impute"               "SKIP" "input pg.matrix_raw.tsv absent from project tree" 0
    Record "04_format_metadata"      "SKIP" "input sample_metadata.xlsx absent from project tree" 0
    Record "01_qc_protein_peptide"   "SKIP" "input quicksearch.stats.annotated.xlsx absent (known gap)" 0
    Record "05_metadata_create_EXP9" "SKIP" "out of scope: Exp9_Social-Stress; quarantined in 99_out_of_scope/, guarded by PROTEOMICS_ALLOW_EXP9" 0
}

# ------------------------------------------------------------ 4. slow stages
if ($Tier -eq "Full") {
    Write-Host ""
    Write-Host "4. Slow stages (ID mapping, enrichment, EWCE)" -ForegroundColor Cyan

    $mapOut = Join-Path $scratch "mapped"
    Invoke-Stage "01_MapThatProt_batch" "02_id_mapping\01_MapThatProt_batch.r" `
        @{ PROTEOMICS_MAPTHATPROT_SPLIT_ROOT  = (Join-Path $dataRoot "02_data\animal_level\split");
           PROTEOMICS_MAPTHATPROT_OUTPUT_ROOT = $mapOut } `
        @( (Join-Path $dataRoot "01_input\references\MOUSE_10090_idmapping.dat"),
           (Join-Path $dataRoot "02_data\animal_level\split\indexComparisons.csv") )

    Invoke-Stage "01_clusterProfiler (GSEA/ORA)" "04_differential_expression_enrichment\01_clusterProfiler.r" `
        @{ PROTEOMICS_ENRICHMENT_MAPPED_ROOT  = (Join-Path $dataRoot "02_data\animal_level\mapped");
           PROTEOMICS_ENRICHMENT_MAPPED_INDEX = (Join-Path $dataRoot "02_data\animal_level\mapped\indexMappedComparisons.csv");
           PROTEOMICS_ENRICHMENT_OUTPUT_ROOT  = (Join-Path $scratch "enrichment") } `
        @( (Join-Path $dataRoot "02_data\animal_level\mapped\indexMappedComparisons.csv") )

    Invoke-Stage "01_EWCE (10k bootstrap reps)" "05_celltype_enrichment_EWCE\01_EWCE.r" `
        @{ PROTEOMICS_EWCE_OUTPUT_ROOT = (Join-Path $scratch "ewce") } `
        @( (Join-Path $dataRoot "02_data\animal_level\input_gct\neha_protigy_input_animal_level_primary.gct") )
} elseif ($Tier -eq "Normal") {
    Write-Host ""
    Write-Host "Tier=Normal: skipped ID mapping / enrichment / EWCE. Use -Tier Full for those." -ForegroundColor Yellow
}

# ----------------------------------------------------------------- summary
Write-Host ""
Write-Host "===== SUMMARY =====" -ForegroundColor Cyan
$results | Format-Table -AutoSize
# NOTE: wrap in @() -- in PowerShell 5.1 a single-item pipeline result is a scalar whose
# .Count is empty, which would silently under-report failures.
$pass  = @($results | Where-Object { $_.Status -eq "PASS" }).Count
$fail  = @($results | Where-Object { $_.Status -eq "FAIL" }).Count
$skip  = @($results | Where-Object { $_.Status -eq "SKIP" }).Count
$abort = @($results | Where-Object { $_.Status -eq "ABORT" }).Count
Write-Host ("PASS {0}   FAIL {1}   SKIP {2}   ABORT {3}" -f $pass, $fail, $skip, $abort)
if ($abort -gt 0) {
    Write-Host ""
    Write-Host "SAFETY ABORT(S): a fail-closed gate stopped a stage BEFORE R was invoked." -ForegroundColor Magenta
    Write-Host "No producing stage ran for those entries, so nothing was written anywhere." -ForegroundColor Magenta
}
$results | Export-Csv -LiteralPath (Join-Path $scratch "summary.csv") -NoTypeInformation
Write-Host "logs + summary: $scratch"

# restore any PROTEOMICS_* overrides we cleared at startup
foreach ($k in $inheritedOverride.Keys) { Set-Item -Path ("Env:" + $k) -Value $inheritedOverride[$k] }
if ($inheritedOverride.Count -gt 0) {
    Write-Host ("Restored {0} inherited PROTEOMICS_* override(s) in your shell." -f $inheritedOverride.Count) -ForegroundColor Yellow
}

if ($fail -gt 0 -or $abort -gt 0) {
    Write-Host "FAILURES OR SAFETY ABORTS PRESENT - inspect the output above." -ForegroundColor Red
    exit 1
}
Write-Host "No failures. Validated outputs were not modified." -ForegroundColor Green
exit 0
