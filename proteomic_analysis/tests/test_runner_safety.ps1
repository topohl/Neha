<#
    Contracts for the fail-closed safety gates in run_pipeline_safety.ps1.

    These cover the exact 2026-08-28 incident class: a runner loaded from one commit, a tree
    checked out to another mid-run, stage scripts expecting the pre-rename NEHA_* variables,
    and producing stages silently writing to canonical shared-drive defaults.

    Everything here runs against temporary fixtures. NOTHING touches the real shared drive:
    the protected-root checks are string/containment assertions against synthetic paths, and
    no R stage is ever invoked.
#>

$ErrorActionPreference = "Stop"
. (Join-Path (Split-Path $PSScriptRoot -Parent) "run_pipeline_safety.ps1")

$script:pass = 0
$script:fail = 0

function Test-Case([string]$name, [scriptblock]$body) {
    try {
        & $body
        $script:pass++
        Write-Host ("  [ok]    " + $name) -ForegroundColor Green
    } catch {
        $script:fail++
        Write-Host ("  [FAIL]  " + $name) -ForegroundColor Red
        Write-Host ("            " + $_.Exception.Message) -ForegroundColor Red
    }
}

function Assert-True($cond, $msg) { if (-not $cond) { throw $msg } }

function Assert-Throws([scriptblock]$body, [string]$mustContain, [string]$msg) {
    $threw = $false
    $seen = ""
    try { & $body } catch { $threw = $true; $seen = $_.Exception.Message }
    if (-not $threw) { throw ($msg + " (nothing was thrown)") }
    if ($mustContain -ne "" -and $seen -notlike ("*" + $mustContain + "*")) {
        throw ($msg + " (message did not mention '" + $mustContain + "'): " + $seen)
    }
}

# --------------------------------------------------------------------------- fixtures
$tmp = Join-Path $env:TEMP ("runner_safety_test_" + [guid]::NewGuid().ToString("N").Substring(0, 8))
$scratch = Join-Path $tmp "scratch"
$fakeRepo = Join-Path $tmp "repo"
New-Item -ItemType Directory -Force -Path $scratch | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $fakeRepo "stages") | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $fakeRepo "R") | Out-Null

# A stage fixture that correctly consumes the modern override, via a sourced helper --
# mirroring how the real PCA and enrichment stages resolve their roots.
$modernStage = Join-Path $fakeRepo "stages\modern_stage.r"
@'
# fixture: resolves its output root through a sourced helper
source(file.path(repo_root, "R", "modern_helper.R"))
'@ | Set-Content -LiteralPath $modernStage -Encoding utf8
@'
out <- option_or_env("proteomics.demo_output_root", "PROTEOMICS_DEMO_OUTPUT_ROOT", default_root)
'@ | Set-Content -LiteralPath (Join-Path $fakeRepo "R\modern_helper.R") -Encoding utf8

# A stage fixture frozen at the pre-rename generation: it only knows NEHA_*.
# This is precisely the file the incident ran.
$legacyStage = Join-Path $fakeRepo "stages\legacy_stage.r"
@'
out <- option_or_env("neha.demo_output_root", "NEHA_DEMO_OUTPUT_ROOT", default_root)
'@ | Set-Content -LiteralPath $legacyStage -Encoding utf8

# A stage fixture that names no override at all.
$muteStage = Join-Path $fakeRepo "stages\mute_stage.r"
'out <- default_root' | Set-Content -LiteralPath $muteStage -Encoding utf8

Write-Host ""
Write-Host "=== 1. normal current-tree execution ===" -ForegroundColor Cyan

Test-Case "modern stage: override contract recognised (via sourced helper)" {
    [void](Assert-OverrideContract -RepoRoot $fakeRepo -StageFullPath $modernStage `
            -EnvNames @("PROTEOMICS_DEMO_OUTPUT_ROOT") -StageName "modern")
}

Test-Case "scratch destination accepted" {
    [void](Assert-SafeOutputTarget -EnvName "PROTEOMICS_DEMO_OUTPUT_ROOT" `
            -Value (Join-Path $scratch "demo") -ScratchRoot $scratch -StageName "modern")
}

Test-Case "identical repo fingerprints report no drift" {
    $s = @{ Root = "c:/x"; Head = "abc123"; Branch = "main" }
    [void](Assert-NoCheckoutDrift -StartupState $s -CurrentState $s.Clone() -StageName "modern")
}

Test-Case "env map with only its own overrides passes the stale-variable gate" {
    $env:PROTEOMICS_DEMO_OUTPUT_ROOT = (Join-Path $scratch "demo")
    try {
        [void](Assert-NoStaleOverridesPresent -EnvMap @{ PROTEOMICS_DEMO_OUTPUT_ROOT = "x" } -StageName "modern")
    } finally {
        Remove-Item Env:PROTEOMICS_DEMO_OUTPUT_ROOT -ErrorAction SilentlyContinue
    }
}

Write-Host ""
Write-Host "=== 2. override-contract mismatch (the incident's root cause) ===" -ForegroundColor Cyan

Test-Case "runner sets PROTEOMICS_*, stage only knows NEHA_* -> ABORT" {
    Assert-Throws {
        [void](Assert-OverrideContract -RepoRoot $fakeRepo -StageFullPath $legacyStage `
                -EnvNames @("PROTEOMICS_DEMO_OUTPUT_ROOT") -StageName "legacy")
    } "OVERRIDE CONTRACT NOT SATISFIED" "a legacy stage must not be accepted"
}

Test-Case "the abort names the pre-rename variable the stage actually expects" {
    Assert-Throws {
        [void](Assert-OverrideContract -RepoRoot $fakeRepo -StageFullPath $legacyStage `
                -EnvNames @("PROTEOMICS_DEMO_OUTPUT_ROOT") -StageName "legacy")
    } "NEHA_DEMO_OUTPUT_ROOT" "the diagnostic should point at the legacy name"
}

Test-Case "stage naming no override at all -> ABORT" {
    Assert-Throws {
        [void](Assert-OverrideContract -RepoRoot $fakeRepo -StageFullPath $muteStage `
                -EnvNames @("PROTEOMICS_DEMO_OUTPUT_ROOT") -StageName "mute")
    } "OVERRIDE CONTRACT NOT SATISFIED" "a stage with no override must not be accepted"
}

Test-Case "no output directory is created by a contract abort" {
    $target = Join-Path $scratch "never_created"
    try {
        [void](Assert-OverrideContract -RepoRoot $fakeRepo -StageFullPath $legacyStage `
                -EnvNames @("PROTEOMICS_DEMO_OUTPUT_ROOT") -StageName "legacy")
    } catch { }
    Assert-True (-not (Test-Path -LiteralPath $target)) "abort must not create the output directory"
}

Write-Host ""
Write-Host "=== 3. checkout / HEAD drift ===" -ForegroundColor Cyan

Test-Case "HEAD change between startup and stage -> ABORT" {
    $startup = @{ Root = "c:/x"; Head = "aaaaaaaaaaaa"; Branch = "rebrand" }
    $current = @{ Root = "c:/x"; Head = "bbbbbbbbbbbb"; Branch = "main" }
    Assert-Throws {
        [void](Assert-NoCheckoutDrift -StartupState $startup -CurrentState $current -StageName "02a")
    } "REPOSITORY CHECKOUT DRIFT DETECTED" "a HEAD change must abort"
}

Test-Case "the drift abort reports both commits and both branches" {
    $startup = @{ Root = "c:/x"; Head = "aaaaaaaaaaaa"; Branch = "rebrand" }
    $current = @{ Root = "c:/x"; Head = "bbbbbbbbbbbb"; Branch = "main" }
    $seen = ""
    try { [void](Assert-NoCheckoutDrift -StartupState $startup -CurrentState $current -StageName "02a") } catch { $seen = $_.Exception.Message }
    foreach ($needle in @("aaaaaaaaaaaa", "bbbbbbbbbbbb", "rebrand", "main")) {
        Assert-True ($seen -like ("*" + $needle + "*")) ("drift message should mention " + $needle)
    }
}

Test-Case "repository root change -> ABORT" {
    Assert-Throws {
        [void](Assert-NoCheckoutDrift -StartupState @{ Root = "c:/a"; Head = "h"; Branch = "b" } `
                -CurrentState @{ Root = "c:/other"; Head = "h"; Branch = "b" } -StageName "02a")
    } "REPOSITORY CHECKOUT DRIFT DETECTED" "a root change must abort"
}

Test-Case "renamed-away stage script -> ABORT (the PCA symptom in the incident)" {
    Assert-Throws {
        [void](Assert-StageScriptResolves -RepoRoot $fakeRepo -StageRelPath "stages\vanished_stage.r" -StageName "pca")
    } "stage script not found" "a missing stage script must abort"
}

Test-Case "real repository fingerprint is readable and stable across two reads" {
    $repoRoot = Split-Path $PSScriptRoot -Parent
    $a = Get-GitRepoState -RepoRoot $repoRoot
    $b = Get-GitRepoState -RepoRoot $repoRoot
    Assert-True ($a.Head -eq $b.Head) "two consecutive HEAD reads disagreed"
    Assert-True ($a.Head.Length -ge 40) "HEAD does not look like a full SHA"
    [void](Assert-NoCheckoutDrift -StartupState $a -CurrentState $b -StageName "self")
}

Write-Host ""
Write-Host "=== 4. unsafe output root ===" -ForegroundColor Cyan

Test-Case "destination outside the scratch root -> ABORT" {
    Assert-Throws {
        [void](Assert-SafeOutputTarget -EnvName "PROTEOMICS_DEMO_OUTPUT_ROOT" `
                -Value (Join-Path $tmp "outside_scratch") -ScratchRoot $scratch -StageName "demo")
    } "OUTSIDE this run's scratch root" "a non-scratch destination must abort"
}

# Synthetic protected-root strings only -- these are never opened, only compared.
Test-Case "destination inside the protected 02_data root -> ABORT" {
    Assert-Throws {
        [void](Assert-SafeOutputTarget -EnvName "PROTEOMICS_ANIMAL_LEVEL_OUTPUT_DIR" `
                -Value "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/02_data/animal_level/input_gct" `
                -ScratchRoot $scratch -StageName "02a")
    } "PROTECTED canonical root" "the exact incident destination must be rejected"
}

Test-Case "destination inside the protected 03_output root -> ABORT" {
    Assert-Throws {
        [void](Assert-SafeOutputTarget -EnvName "PROTEOMICS_PCA_OUTPUT_ROOT" `
                -Value "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/03_output/pca" `
                -ScratchRoot $scratch -StageName "pca")
    } "PROTECTED canonical root" "03_output must be rejected"
}

Test-Case "destination inside the protected 99_historical root -> ABORT" {
    Assert-Throws {
        [void](Assert-SafeOutputTarget -EnvName "PROTEOMICS_ENRICHMENT_OUTPUT_ROOT" `
                -Value "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/99_historical/whatever" `
                -ScratchRoot $scratch -StageName "enrichment")
    } "PROTECTED canonical root" "99_historical must be rejected"
}

Test-Case "protected root is rejected in UNC spelling too" {
    Assert-Throws {
        [void](Assert-SafeOutputTarget -EnvName "PROTEOMICS_ANIMAL_LEVEL_OUTPUT_DIR" `
                -Value "\\mdc-berlin.net\fs\AG_Hoernberg\Lab_Member\Tobi\Experiments\Collabs\Neha\clusterProfiler\02_data\animal_level\split" `
                -ScratchRoot $scratch -StageName "03_split")
    } "PROTECTED canonical root" "the UNC spelling of the protected root must also be rejected"
}

Test-Case "protected root is rejected with backslashes and mixed case" {
    Assert-Throws {
        [void](Assert-SafeOutputTarget -EnvName "PROTEOMICS_PROTIGY_STAT_GCT_OUTPUT_ROOT" `
                -Value "s:\LAB_MEMBER\Tobi\Experiments\Collabs\NEHA\clusterProfiler\02_DATA\animal_level\split" `
                -ScratchRoot $scratch -StageName "03_split")
    } "PROTECTED canonical root" "separator/case variation must not bypass the protected-root check"
}

Test-Case "empty override value -> ABORT" {
    Assert-Throws {
        [void](Assert-SafeOutputTarget -EnvName "PROTEOMICS_DEMO_OUTPUT_ROOT" -Value "" `
                -ScratchRoot $scratch -StageName "demo")
    } "is empty" "an empty override must abort"
}

Test-Case "env-map sweep checks only destination keys, leaving read-only inputs alone" {
    # Input overrides legitimately point at the shared drive; only destinations are constrained.
    [void](Assert-SafeOutputTargets -EnvMap @{
        PROTEOMICS_MAPTHATPROT_SPLIT_ROOT  = "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/02_data/animal_level/split"
        PROTEOMICS_MAPTHATPROT_OUTPUT_ROOT = (Join-Path $scratch "mapped")
    } -ScratchRoot $scratch -StageName "mapping")
}

Test-Case "env-map sweep still catches a canonical destination among safe inputs" {
    Assert-Throws {
        [void](Assert-SafeOutputTargets -EnvMap @{
            PROTEOMICS_MAPTHATPROT_SPLIT_ROOT  = "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/02_data/animal_level/split"
            PROTEOMICS_MAPTHATPROT_OUTPUT_ROOT = "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/02_data/animal_level/mapped"
        } -ScratchRoot $scratch -StageName "mapping")
    } "PROTECTED canonical root" "a canonical destination must be caught even alongside valid inputs"
}

Write-Host ""
Write-Host "=== 5. environment hygiene ===" -ForegroundColor Cyan

Test-Case "a stale NEHA_* variable in the parent shell is swept and restored" {
    $env:NEHA_STALE_PROBE = "canonical-path"
    $saved = Clear-InheritedOverrides -Prefixes @("PROTEOMICS_", "NEHA_")
    try {
        Assert-True ($saved.ContainsKey("NEHA_STALE_PROBE")) "legacy variable was not swept"
        Assert-True ($null -eq $env:NEHA_STALE_PROBE) "legacy variable still set after sweep"
        [void](Restore-InheritedOverrides -Saved $saved)
        Assert-True ($env:NEHA_STALE_PROBE -eq "canonical-path") "legacy variable was not restored"
    } finally {
        Remove-Item Env:NEHA_STALE_PROBE -ErrorAction SilentlyContinue
    }
}

Test-Case "an unexpected override present at invoke time -> ABORT" {
    $env:PROTEOMICS_UNEXPECTED_PROBE = "somewhere-else"
    try {
        Assert-Throws {
            [void](Assert-NoStaleOverridesPresent -EnvMap @{ PROTEOMICS_DEMO_OUTPUT_ROOT = "x" } -StageName "demo")
        } "unexpected override variable" "an unexpected override must abort"
    } finally {
        Remove-Item Env:PROTEOMICS_UNEXPECTED_PROBE -ErrorAction SilentlyContinue
    }
}

Write-Host ""
Write-Host "=== 6. the real runner wires the gates in ===" -ForegroundColor Cyan

Test-Case "run_pipeline_check.ps1 dot-sources the safety library and gates Invoke-Stage" {
    $runner = Get-Content -LiteralPath (Join-Path (Split-Path $PSScriptRoot -Parent) "run_pipeline_check.ps1") -Raw
    foreach ($needle in @(
        "run_pipeline_safety.ps1",
        "Get-GitRepoState",
        "Assert-NoCheckoutDrift",
        "Assert-StageScriptResolves",
        "Assert-OverrideContract",
        "Assert-SafeOutputTargets",
        "Assert-NoStaleOverridesPresent")) {
        Assert-True ($runner -like ("*" + $needle + "*")) ("runner does not reference " + $needle)
    }
    Assert-True ($runner -like "*ABORT*") "runner does not record an ABORT status"
}

Test-Case "every real producing stage satisfies its own override contract right now" {
    $repoRoot = Split-Path $PSScriptRoot -Parent
    $stages = @{
        "01_preprocessing\02a_prepare_animal_level_protigy_input.r"        = @("PROTEOMICS_ANIMAL_LEVEL_OUTPUT_DIR")
        "01_preprocessing\03_gct_extractR.r"                               = @("PROTEOMICS_PROTIGY_STAT_GCT_OUTPUT_ROOT")
        "03_qc_exploration\06_pcaPlot_animal_level.r"                      = @("PROTEOMICS_PCA_OUTPUT_ROOT")
        "03_qc_exploration\02_rank_abundance_by_sample_class.r"            = @("PROTEOMICS_RANK_ABUNDANCE_OUTPUT_DIR")
        "02_id_mapping\01_MapThatProt_batch.r"                             = @("PROTEOMICS_MAPTHATPROT_SPLIT_ROOT", "PROTEOMICS_MAPTHATPROT_OUTPUT_ROOT")
        "04_differential_expression_enrichment\01_clusterProfiler.r"       = @("PROTEOMICS_ENRICHMENT_MAPPED_ROOT", "PROTEOMICS_ENRICHMENT_MAPPED_INDEX", "PROTEOMICS_ENRICHMENT_OUTPUT_ROOT")
        "05_celltype_enrichment_EWCE\01_EWCE.r"                            = @("PROTEOMICS_EWCE_OUTPUT_ROOT")
    }
    foreach ($rel in $stages.Keys) {
        $full = Assert-StageScriptResolves -RepoRoot $repoRoot -StageRelPath $rel -StageName $rel
        [void](Assert-OverrideContract -RepoRoot $repoRoot -StageFullPath $full -EnvNames $stages[$rel] -StageName $rel)
    }
}

# --------------------------------------------------------------------------- teardown
Remove-Item -Recurse -Force -LiteralPath $tmp -ErrorAction SilentlyContinue

Write-Host ""
Write-Host "=== RESULT ===" -ForegroundColor Cyan
Write-Host ("passed: {0}   failed: {1}" -f $script:pass, $script:fail)
if ($script:fail -gt 0) {
    Write-Host "Runner-safety contracts FAILED." -ForegroundColor Red
    exit 1
}
Write-Host "All runner-safety contracts hold." -ForegroundColor Green
exit 0
