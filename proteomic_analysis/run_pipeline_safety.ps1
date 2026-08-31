<#
    Fail-closed safety helpers for run_pipeline_check.ps1.

    WHY THIS EXISTS (incident 2026-08-28 16:36 UTC+2)
    -------------------------------------------------
    A smoke run was started from the rebrand branch. PowerShell had already parsed and loaded
    the runner, which sets PROTEOMICS_* overrides. Mid-run the working tree was checked out
    back to `main`, whose stage scripts still read the pre-rename NEHA_* variables. The
    overrides were therefore silently ignored, every stage fell back to its hardcoded default,
    and two producing stages wrote 37 files into the VALIDATED shared-drive tree
    (02_data/animal_level/{input_gct,split}) instead of the scratch root.

    Nothing detected it at the time: the runner had no notion of repository identity, no check
    that a stage actually understands the override it is being handed, and no check that the
    resolved destination is inside scratch. Each of those is now a hard gate.

    These are separate, individually testable functions on purpose -- see
    tests/test_runner_safety.ps1. They THROW on violation and return $true when satisfied, so
    a caller cannot accidentally treat a failure as a pass.

    PowerShell 5.1 compatible: no ternary, no null-coalescing, no inline if-expressions.
#>

# NOTE: deliberately no Set-StrictMode here. This file is dot-sourced into the runner, and
# strict mode would leak into the caller's scope and change unrelated behaviour. Every
# function below validates its own inputs instead (Mandatory parameters + explicit checks).

# The canonical roots a smoke run must never write into, as forward-slash lowercase path
# fragments. Matching is segment-aligned CONTAINMENT rather than prefix matching, so the same
# fragment catches every spelling of the same location:
#   S:\...\Collabs\Neha\clusterProfiler\02_data\...              (mapped drive, backslashes)
#   \\mdc-berlin.net\fs\AG_Hoernberg\...\clusterProfiler\02_data (UNC, share prefix in front)
# A prefix-only check missed the UNC form, which is exactly the kind of near-miss that lets a
# protected write through.
$script:ProtectedRootFragments = @(
    'collabs/neha/clusterprofiler/02_data',
    'collabs/neha/clusterprofiler/03_output',
    'collabs/neha/clusterprofiler/99_historical'
)

function Get-ProtectedRootFragments {
    return $script:ProtectedRootFragments
}

function ConvertTo-ComparablePath {
    <# Normalise a path for containment comparison: absolute, forward slashes, lowercase,
       no trailing separator. Works on paths that do not exist yet (scratch subdirectories). #>
    param([Parameter(Mandatory = $true)][AllowEmptyString()][string]$Path)
    if ([string]::IsNullOrWhiteSpace($Path)) { return '' }
    $p = $Path.Trim()
    try { $p = [System.IO.Path]::GetFullPath($p) } catch { }
    $p = $p.Replace('\', '/').ToLowerInvariant()
    while ($p.Length -gt 1 -and $p.EndsWith('/')) { $p = $p.Substring(0, $p.Length - 1) }
    return $p
}

function Test-PathWithin {
    <# Is $Path inside $Parent (or equal to it)? #>
    param(
        [Parameter(Mandatory = $true)][AllowEmptyString()][string]$Path,
        [Parameter(Mandatory = $true)][AllowEmptyString()][string]$Parent
    )
    $a = ConvertTo-ComparablePath $Path
    $b = ConvertTo-ComparablePath $Parent
    if ($a -eq '' -or $b -eq '') { return $false }
    if ($a -eq $b) { return $true }
    return $a.StartsWith($b + '/')
}

# ------------------------------------------------------------------ A. repository identity

function Get-GitRepoState {
    <# Fingerprint the checked-out tree. Fails closed: if the repository identity cannot be
       established we refuse to run producing stages at all, because the whole point is to
       detect the tree changing underneath us. #>
    param([Parameter(Mandatory = $true)][string]$RepoRoot)

    $git = Get-Command git -ErrorAction SilentlyContinue
    if ($null -eq $git) {
        throw "SAFETY ABORT: 'git' is not available, so repository checkout drift cannot be detected. Refusing to run producing stages."
    }
    $root = (& git -C $RepoRoot rev-parse --show-toplevel 2>$null)
    if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($root)) {
        throw "SAFETY ABORT: '$RepoRoot' is not inside a git repository, so checkout drift cannot be detected. Refusing to run producing stages."
    }
    $head = (& git -C $RepoRoot rev-parse HEAD 2>$null)
    if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($head)) {
        throw "SAFETY ABORT: could not read HEAD for '$RepoRoot'. Refusing to run producing stages."
    }
    $branch = (& git -C $RepoRoot rev-parse --abbrev-ref HEAD 2>$null)
    if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($branch)) { $branch = '(detached)' }

    return @{
        Root   = (ConvertTo-ComparablePath $root.Trim())
        Head   = $head.Trim()
        Branch = $branch.Trim()
    }
}

function Assert-NoCheckoutDrift {
    <# Pure comparison of two repo fingerprints. This is the gate that would have stopped the
       2026-08-28 incident: the runner was loaded on one commit and the tree became another. #>
    param(
        [Parameter(Mandatory = $true)][hashtable]$StartupState,
        [Parameter(Mandatory = $true)][hashtable]$CurrentState,
        [Parameter(Mandatory = $true)][string]$StageName
    )
    if ($StartupState.Root -ne $CurrentState.Root) {
        throw ("SAFETY ABORT before stage '{0}': REPOSITORY CHECKOUT DRIFT DETECTED DURING THE SMOKE RUN. " -f $StageName) +
              ("Repository root changed from '{0}' to '{1}'. " -f $StartupState.Root, $CurrentState.Root) +
              "The loaded runner may no longer match the scripts on disk; refusing to invoke R."
    }
    if ($StartupState.Head -ne $CurrentState.Head) {
        throw ("SAFETY ABORT before stage '{0}': REPOSITORY CHECKOUT DRIFT DETECTED DURING THE SMOKE RUN. " -f $StageName) +
              ("HEAD changed from {0} ({1}) to {2} ({3}). " -f $StartupState.Head, $StartupState.Branch, $CurrentState.Head, $CurrentState.Branch) +
              "The stage scripts on disk are no longer the ones this run started with, so the " +
              "override contract cannot be trusted and outputs could land outside scratch. Refusing to invoke R."
    }
    return $true
}

function Assert-StageScriptResolves {
    <# The stage file must exist and live inside the repository we fingerprinted. #>
    param(
        [Parameter(Mandatory = $true)][string]$RepoRoot,
        [Parameter(Mandatory = $true)][string]$StageRelPath,
        [Parameter(Mandatory = $true)][string]$StageName
    )
    $full = Join-Path $RepoRoot $StageRelPath
    if (-not (Test-Path -LiteralPath $full -PathType Leaf)) {
        throw ("SAFETY ABORT before stage '{0}': stage script not found at '{1}'. " -f $StageName, $full) +
              "This is the signature of a renamed or checked-out-away script; refusing to invoke R."
    }
    if (-not (Test-PathWithin $full $RepoRoot)) {
        throw ("SAFETY ABORT before stage '{0}': stage script '{1}' resolves outside the repository root '{2}'." -f $StageName, $full, $RepoRoot)
    }
    return $full
}

# ------------------------------------------------------- B. override-contract validation

function Get-StageContractSearchSet {
    <# The files that together define a stage's override contract: the stage script plus the
       files it source()s, followed two levels deep.

       Following source() matters: the PCA orchestrator does not name PROTEOMICS_PCA_OUTPUT_ROOT
       itself (pca/06a_pca_core.r does), and the mapping/enrichment stages resolve their roots
       through helpers in R/. A check that only read the stage file would pass vacuously. #>
    param(
        [Parameter(Mandatory = $true)][string]$RepoRoot,
        [Parameter(Mandatory = $true)][string]$StageFullPath,
        [int]$Depth = 2
    )
    $set = New-Object System.Collections.Generic.List[string]
    $set.Add($StageFullPath)
    $seen = @{}
    $seen[(ConvertTo-ComparablePath $StageFullPath)] = $true

    # Index every R file in the repo once, by lowercase basename, excluding frozen provenance.
    $index = @{}
    Get-ChildItem -LiteralPath $RepoRoot -Recurse -File -Include *.R, *.r -ErrorAction SilentlyContinue |
        Where-Object { $_.FullName -notmatch '06_manuscript_figure_revision' } |
        ForEach-Object {
            $k = $_.Name.ToLowerInvariant()
            if (-not $index.ContainsKey($k)) { $index[$k] = New-Object System.Collections.Generic.List[string] }
            $index[$k].Add($_.FullName)
        }

    $frontier = @($StageFullPath)
    for ($d = 0; $d -lt $Depth; $d++) {
        $next = New-Object System.Collections.Generic.List[string]
        foreach ($f in $frontier) {
            if (-not (Test-Path -LiteralPath $f -PathType Leaf)) { continue }
            $text = Get-Content -LiteralPath $f -Raw -ErrorAction SilentlyContinue
            if ($null -eq $text) { continue }
            # Any quoted *.R / *.r literal. Catches both source(file.path(..,"x.R")) and the
            # orchestrator's list of part filenames.
            foreach ($m in [regex]::Matches($text, '"([^"\\/:*?<>|]+\.[Rr])"')) {
                $base = $m.Groups[1].Value.ToLowerInvariant()
                if (-not $index.ContainsKey($base)) { continue }
                foreach ($cand in $index[$base]) {
                    $key = ConvertTo-ComparablePath $cand
                    if ($seen.ContainsKey($key)) { continue }
                    $seen[$key] = $true
                    $set.Add($cand)
                    $next.Add($cand)
                }
            }
        }
        if ($next.Count -eq 0) { break }
        $frontier = $next.ToArray()
    }
    return $set.ToArray()
}

function Test-OverrideRecognized {
    <# Does any file in the contract search set literally reference $EnvName? #>
    param(
        [Parameter(Mandatory = $true)][string[]]$SearchSet,
        [Parameter(Mandatory = $true)][string]$EnvName
    )
    foreach ($f in $SearchSet) {
        if (-not (Test-Path -LiteralPath $f -PathType Leaf)) { continue }
        $text = Get-Content -LiteralPath $f -Raw -ErrorAction SilentlyContinue
        if ($null -eq $text) { continue }
        if ($text.Contains($EnvName)) { return $true }
    }
    return $false
}

function Assert-OverrideContract {
    <# Every override the runner is about to set must actually be consumed by the code that is
       about to run. This is the direct guard against the incident's root cause: runner on
       PROTEOMICS_*, stage file on NEHA_*. #>
    param(
        [Parameter(Mandatory = $true)][string]$RepoRoot,
        [Parameter(Mandatory = $true)][string]$StageFullPath,
        [Parameter(Mandatory = $true)][string[]]$EnvNames,
        [Parameter(Mandatory = $true)][string]$StageName
    )
    if ($EnvNames.Count -eq 0) { return $true }
    $searchSet = Get-StageContractSearchSet -RepoRoot $RepoRoot -StageFullPath $StageFullPath
    $missing = New-Object System.Collections.Generic.List[string]
    foreach ($n in $EnvNames) {
        if (-not (Test-OverrideRecognized -SearchSet $searchSet -EnvName $n)) { $missing.Add($n) }
    }
    if ($missing.Count -gt 0) {
        $legacyHint = ''
        foreach ($n in $missing) {
            $legacy = $n -replace '^PROTEOMICS_', 'NEHA_'
            if ($legacy -ne $n -and (Test-OverrideRecognized -SearchSet $searchSet -EnvName $legacy)) {
                $legacyHint += (" Stage code still expects the pre-rename '{0}'." -f $legacy)
            }
        }
        throw ("SAFETY ABORT before stage '{0}': OVERRIDE CONTRACT NOT SATISFIED. " -f $StageName) +
              ("The runner sets {0}, but nothing in the stage's code recognises: {1}." -f ($EnvNames -join ', '), ($missing -join ', ')) +
              $legacyHint +
              " Without a recognised override the stage would silently fall back to its canonical" +
              " shared-drive default. Refusing to invoke R."
    }
    return $true
}

# ------------------------------------------- C. scratch-root / protected-root enforcement

function Assert-SafeOutputTarget {
    <# The destination the runner hands a stage must be inside this run's scratch root, and
       must never be inside a protected canonical root. #>
    param(
        [Parameter(Mandatory = $true)][string]$EnvName,
        [Parameter(Mandatory = $true)][AllowEmptyString()][string]$Value,
        [Parameter(Mandatory = $true)][string]$ScratchRoot,
        [Parameter(Mandatory = $true)][string]$StageName,
        [string[]]$ProtectedFragments = $null
    )
    if ($null -eq $ProtectedFragments) { $ProtectedFragments = Get-ProtectedRootFragments }

    $norm = ConvertTo-ComparablePath $Value
    if ($norm -eq '') {
        throw ("SAFETY ABORT before stage '{0}': override {1} is empty; refusing to invoke R." -f $StageName, $EnvName)
    }
    # Segment-aligned containment: wrap both sides in '/' so a fragment can only match whole
    # path segments, anywhere in the path (drive-letter prefix or UNC share prefix alike).
    $padded = '/' + $norm.Trim('/') + '/'
    foreach ($frag in $ProtectedFragments) {
        if ($padded.Contains('/' + $frag.Trim('/') + '/')) {
            throw ("SAFETY ABORT before stage '{0}': override {1} resolves to '{2}', which is inside the PROTECTED canonical root '{3}'. " -f $StageName, $EnvName, $norm, $frag) +
                  "A smoke run must never write there. Refusing to invoke R."
        }
    }
    if (-not (Test-PathWithin $Value $ScratchRoot)) {
        throw ("SAFETY ABORT before stage '{0}': override {1} resolves to '{2}', which is OUTSIDE this run's scratch root '{3}'. " -f $StageName, $EnvName, $norm, (ConvertTo-ComparablePath $ScratchRoot)) +
              "Refusing to invoke R."
    }
    return $true
}

function Assert-SafeOutputTargets {
    <# Apply Assert-SafeOutputTarget to the output-directing entries of an env map. Only keys
       that name a destination are checked; input-pointing overrides legitimately reference the
       shared drive read-only. #>
    param(
        [Parameter(Mandatory = $true)][hashtable]$EnvMap,
        [Parameter(Mandatory = $true)][string]$ScratchRoot,
        [Parameter(Mandatory = $true)][string]$StageName,
        [string[]]$ProtectedFragments = $null
    )
    foreach ($k in $EnvMap.Keys) {
        if ($k -match '(OUTPUT|OUT_DIR|OUTPUT_DIR|OUTPUT_ROOT|WORK_DIR)') {
            [void](Assert-SafeOutputTarget -EnvName $k -Value ([string]$EnvMap[$k]) -ScratchRoot $ScratchRoot -StageName $StageName -ProtectedFragments $ProtectedFragments)
        }
    }
    return $true
}

# ------------------------------------------------------------------ D. environment hygiene

function Clear-InheritedOverrides {
    <# Remove every environment variable matching the given prefixes and return what was
       removed so the caller can restore it. Both the current PROTEOMICS_* names and the
       pre-rename NEHA_* names are swept: a stale legacy variable in the parent shell could
       still steer an older stage script to a canonical path. #>
    param([string[]]$Prefixes = @('PROTEOMICS_', 'NEHA_'))
    $saved = @{}
    foreach ($prefix in $Prefixes) {
        Get-ChildItem Env: | Where-Object { $_.Name -like ($prefix + '*') } | ForEach-Object {
            $saved[$_.Name] = $_.Value
            Remove-Item -Path ('Env:' + $_.Name) -ErrorAction SilentlyContinue
        }
    }
    return $saved
}

function Restore-InheritedOverrides {
    param([Parameter(Mandatory = $true)][hashtable]$Saved)
    foreach ($k in $Saved.Keys) { Set-Item -Path ('Env:' + $k) -Value $Saved[$k] }
    return $Saved.Count
}

function Assert-NoStaleOverridesPresent {
    <# Belt-and-braces: at the moment of invoking a stage, the only PROTEOMICS_*/NEHA_*
       variables in the environment must be the ones this stage deliberately set. #>
    param(
        [Parameter(Mandatory = $true)][hashtable]$EnvMap,
        [Parameter(Mandatory = $true)][string]$StageName,
        [string[]]$Prefixes = @('PROTEOMICS_', 'NEHA_')
    )
    $unexpected = New-Object System.Collections.Generic.List[string]
    foreach ($prefix in $Prefixes) {
        Get-ChildItem Env: | Where-Object { $_.Name -like ($prefix + '*') } | ForEach-Object {
            if (-not $EnvMap.ContainsKey($_.Name)) { $unexpected.Add($_.Name) }
        }
    }
    if ($unexpected.Count -gt 0) {
        throw ("SAFETY ABORT before stage '{0}': unexpected override variable(s) present in the environment: {1}. " -f $StageName, (($unexpected | Sort-Object -Unique) -join ', ')) +
              "These could redirect the stage away from scratch. Refusing to invoke R."
    }
    return $true
}
