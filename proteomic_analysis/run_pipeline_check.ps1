<#
    Smoke test for the Associative Memory Proteomics pipeline. Runs the tests, then the
    producing stages with every output redirected to a scratch folder under %TEMP% through
    each stage's env-var override. Validated outputs on S: are only read.

    Tiers:
      -Tier Fast    tests only, including the data-integrity test        (~2 min)
      -Tier Normal  Fast + GCT rebuild + split + PCA + rank abundance   (~10-20 min) [default]
      -Tier Full    Normal + ID mapping + enrichment + EWCE             (hours; EWCE is 10k reps)

    Usage:
      powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1
      powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1 -Tier Fast
      powershell -ExecutionPolicy Bypass -File .\run_pipeline_check.ps1 -Tier Full

    Coverage: the chain is
        02a_prepare_animal_level -> 02b_limma_stat_gct -> 03_gct_extractR -> ...
    01_preprocessing/02b_limma_stat_gct.r produces the primary per-compartment statistics. It
    is not run as a stage here: tests/test_limma_stat_gct.R rebuilds its output in a temporary
    folder and checks it against the locked hash. ProTigy, an external tool, made only the
    superseded pooled statistics and is listed as SKIP. Stage 3c and the Full tier read the
    per-compartment files on S:.
#>

param(
    [ValidateSet("Fast", "Normal", "Full")]
    [string]$Tier = "Normal"
)

$ErrorActionPreference = "Continue"

$repo    = $PSScriptRoot
# Rscript: the R 4.5.1 install under the current user's LOCALAPPDATA, else Rscript.exe on PATH.
$rscript = ""
if (-not [string]::IsNullOrEmpty($env:LOCALAPPDATA)) {
    $rscript = Join-Path $env:LOCALAPPDATA "Programs\R\R-4.5.1\bin\x64\Rscript.exe"
}
if ($rscript -eq "" -or -not (Test-Path -LiteralPath $rscript)) {
    $onPath = Get-Command Rscript.exe -CommandType Application -ErrorAction SilentlyContinue | Select-Object -First 1
    if ($null -eq $onPath) {
        Write-Host "Rscript.exe not found: install R 4.5.1 under %LOCALAPPDATA%\Programs\R or put Rscript.exe on PATH." -ForegroundColor Red
        exit 1
    }
    $rscript = $onPath.Path
}
$dataRoot = "S:\Lab_Member\Tobi\Experiments\Collabs\Neha\clusterProfiler"
$scratch = Join-Path $env:TEMP ("proteomics_pipeline_check_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Force -Path $scratch | Out-Null

# Safety checks from run_pipeline_safety.ps1 (its header says why), run before every
# producing stage: same git checkout as at startup, every override read by the stage code,
# every output path inside the scratch folder.
. (Join-Path $PSScriptRoot "run_pipeline_safety.ps1")

# Record the checkout at startup; each producing stage compares against it.
$startupRepoState = Get-GitRepoState -RepoRoot $repo
Write-Host ""
Write-Host "repository pinned at startup:" -ForegroundColor Cyan
Write-Host ("  root   : {0}" -f $startupRepoState.Root)
Write-Host ("  HEAD   : {0}" -f $startupRepoState.Head)
Write-Host ("  branch : {0}" -f $startupRepoState.Branch)

# Clear inherited PROTEOMICS_* and NEHA_* (pre-rename) variables for this run; they are
# listed here and restored at the end. Stages read option -> env var -> default, so a variable
# left in the calling shell would change where the child Rscript processes read and write.
$inheritedOverride = Clear-InheritedOverrides -Prefixes @("PROTEOMICS_", "NEHA_")
if ($inheritedOverride.Count -gt 0) {
    Write-Host ""
    Write-Host "WARNING: cleared inherited override(s) from your shell for this run:" -ForegroundColor Yellow
    foreach ($k in $inheritedOverride.Keys) {
        Write-Host ("  {0} = {1}" -f $k, $inheritedOverride[$k]) -ForegroundColor Yellow
    }
    Write-Host "  (restored when this script exits; set them again if you meant to override)" -ForegroundColor Yellow
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

    # Safety checks, all before Rscript starts. A failed check records ABORT and skips the stage.
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

# 1. Tests
Write-Host "1. Check / unit tests" -ForegroundColor Cyan
# Stop if the checkout has already changed since startup.
try {
    [void](Assert-NoCheckoutDrift -StartupState $startupRepoState `
                                  -CurrentState (Get-GitRepoState -RepoRoot $repo) `
                                  -StageName "tests")
} catch {
    Write-Host $_.Exception.Message -ForegroundColor Magenta
    exit 1
}
# Test the safety checks first. This does not touch the shared drive or run any R stage.
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
    Record "test_runner_safety (ps1)" "FAIL" "safety checks missing from tests/" 0
}

# Every *.R file in tests/ and 07_publication_release/tests/. Release tests that need a
# built release skip when none is reachable.
$testDirs = @(
    (Join-Path $repo "tests"),
    (Join-Path $repo "07_publication_release\tests")
)
foreach ($testDir in $testDirs) {
    if (-not (Test-Path -LiteralPath $testDir)) {
        Record ("tests: " + (Split-Path $testDir -Leaf)) "FAIL" "test directory missing: $testDir" 0
        continue
    }
    Get-ChildItem -LiteralPath $testDir -Filter "*.R" | ForEach-Object {
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
}

# The data-integrity checks are in tests/test_data_integrity.R, run in the loop above.

if ($Tier -eq "Fast") {
    Write-Host ""
    Write-Host "Tier=Fast: stopping before producing stages." -ForegroundColor Yellow
}

# 3. Producing stages
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

    # 3b. ProTigy, an external tool, made only the superseded pooled statistics. The primary
    #     statistics come from 02b_limma_stat_gct.r, rebuilt by tests/test_limma_stat_gct.R above.
    Record "ProTigy_stat_results" "SKIP" "external tool, superseded pooled branch; 02b is checked by test_limma_stat_gct" 0

    # 3c. split the primary per-compartment statistics GCT (read only; output to scratch)
    $splitOut = Join-Path $scratch "split"
    $statGct = Join-Path $dataRoot "02_data\animal_level\within_compartment_20261124\stat_results_within_compartment.gct"
    Invoke-Stage "03_gct_extractR (split)" "01_preprocessing\03_gct_extractR.r" `
        @{ PROTEOMICS_PROTIGY_STAT_GCT_INPUT       = $statGct;
           PROTEOMICS_PROTIGY_STAT_GCT_OUTPUT_ROOT = $splitOut } `
        @( $statGct )

    # 3d. animal-level PCA (independent branch off the validated GCT)
    Invoke-Stage "06_pcaPlot_animal_level (PCA workflow)" "03_qc_exploration\06_pcaPlot_animal_level.r" `
        @{ PROTEOMICS_PCA_OUTPUT_ROOT = (Join-Path $scratch "pca") } `
        @( (Join-Path $dataRoot "02_data\animal_level\input_gct\neha_protigy_input_animal_level_primary.gct") )

    # 3e. animal-level rank-abundance QC (regenerates the Fig 3E / Supp D panels)
    Invoke-Stage "02_rank_abundance (animal level)" "03_qc_exploration\02_rank_abundance_by_sample_class.r" `
        @{ PROTEOMICS_RANK_ABUNDANCE_OUTPUT_DIR = (Join-Path $scratch "rank_abundance") } `
        @( (Join-Path $dataRoot "02_data\animal_level\input_gct\neha_protigy_input_animal_level_primary.gct"),
           (Join-Path $dataRoot "02_data\animal_level\mapped\forward\mcherry_paired_veh_vs_mcherry_unpaired_veh.csv") )

    # 3f. stages whose inputs are not in the project tree
    Record "01_impute"               "SKIP" "input pg.matrix_raw.tsv absent from project tree" 0
    Record "04_format_metadata"      "SKIP" "input sample_metadata.xlsx absent from project tree" 0
    Record "01_qc_protein_peptide"   "SKIP" "input quicksearch.stats.annotated.xlsx absent (known gap)" 0
    Record "05_metadata_create_EXP9" "SKIP" "out of scope: Exp9_Social-Stress; quarantined in 99_out_of_scope/, guarded by PROTEOMICS_ALLOW_EXP9" 0
}

# 4. Slow stages
if ($Tier -eq "Full") {
    Write-Host ""
    Write-Host "4. Slow stages (ID mapping, enrichment, EWCE)" -ForegroundColor Cyan

    # Mapping and enrichment read the per-compartment split/ and mapped/ folders on S:
    # (read only); their outputs go to scratch.
    $mapOut = Join-Path $scratch "mapped"
    Invoke-Stage "01_MapThatProt_batch" "02_id_mapping\01_MapThatProt_batch.r" `
        @{ PROTEOMICS_MAPTHATPROT_SPLIT_ROOT  = (Join-Path $dataRoot "02_data\animal_level\within_compartment_20261124\split");
           PROTEOMICS_MAPTHATPROT_OUTPUT_ROOT = $mapOut } `
        @( (Join-Path $dataRoot "01_input\references\MOUSE_10090_idmapping.dat"),
           (Join-Path $dataRoot "02_data\animal_level\within_compartment_20261124\split\indexComparisons.csv") )

    Invoke-Stage "01_clusterProfiler (GSEA/ORA)" "04_differential_expression_enrichment\01_clusterProfiler.r" `
        @{ PROTEOMICS_ENRICHMENT_MAPPED_ROOT  = (Join-Path $dataRoot "02_data\animal_level\within_compartment_20261124\mapped");
           PROTEOMICS_ENRICHMENT_MAPPED_INDEX = (Join-Path $dataRoot "02_data\animal_level\within_compartment_20261124\mapped\indexMappedComparisons.csv");
           PROTEOMICS_ENRICHMENT_OUTPUT_ROOT  = (Join-Path $scratch "enrichment") } `
        @( (Join-Path $dataRoot "02_data\animal_level\within_compartment_20261124\mapped\indexMappedComparisons.csv") )

    Invoke-Stage "01_EWCE (10k bootstrap reps)" "05_celltype_enrichment_EWCE\01_EWCE.r" `
        @{ PROTEOMICS_EWCE_OUTPUT_ROOT = (Join-Path $scratch "ewce") } `
        @( (Join-Path $dataRoot "02_data\animal_level\input_gct\neha_protigy_input_animal_level_primary.gct") )
} elseif ($Tier -eq "Normal") {
    Write-Host ""
    Write-Host "Tier=Normal: skipped ID mapping / enrichment / EWCE. Use -Tier Full for those." -ForegroundColor Yellow
}

# Summary
Write-Host ""
Write-Host "===== SUMMARY =====" -ForegroundColor Cyan
$results | Format-Table -AutoSize
# @() keeps .Count right when only one row matches (PowerShell 5.1 returns a scalar).
$pass  = @($results | Where-Object { $_.Status -eq "PASS" }).Count
$fail  = @($results | Where-Object { $_.Status -eq "FAIL" }).Count
$skip  = @($results | Where-Object { $_.Status -eq "SKIP" }).Count
$abort = @($results | Where-Object { $_.Status -eq "ABORT" }).Count
Write-Host ("PASS {0}   FAIL {1}   SKIP {2}   ABORT {3}" -f $pass, $fail, $skip, $abort)
if ($abort -gt 0) {
    Write-Host ""
    Write-Host "SAFETY ABORT(S): a safety check stopped a stage before R started." -ForegroundColor Magenta
    Write-Host "No producing stage ran for those entries, so nothing was written anywhere." -ForegroundColor Magenta
}
$results | Export-Csv -LiteralPath (Join-Path $scratch "summary.csv") -NoTypeInformation
Write-Host "logs + summary: $scratch"

# Restore the PROTEOMICS_* and NEHA_* variables cleared at startup.
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
