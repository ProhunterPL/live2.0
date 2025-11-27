# Add Formamide Extended to Phase 2B Analysis
# ===========================================
#
# Run this script after formamide_extended runs complete to add them to analysis.
# The script will automatically detect available runs and include formamide in the complete analysis.
#
# Usage:
#   .\scripts\add_formamide_analysis.ps1

$ErrorActionPreference = "Stop"

$INPUT_DIR = "results/phase2b_additional"
$OUTPUT_DIR = "paper/results_data"

Write-Host "======================================================================"
Write-Host "ADDING FORMAMIDE EXTENDED TO PHASE 2B ANALYSIS"
Write-Host "======================================================================"
Write-Host ""

# Check if formamide directory exists
$formamideDir = Join-Path $INPUT_DIR "formamide_extended"
if (-not (Test-Path $formamideDir)) {
    Write-Host "ERROR: Formamide directory not found: $formamideDir" -ForegroundColor Red
    Write-Host "Make sure formamide_extended runs have completed." -ForegroundColor Yellow
    exit 1
}

# Count completed runs
$runDirs = Get-ChildItem -Path $formamideDir -Directory -Filter "run_*"
$completed = 0
foreach ($runDir in $runDirs) {
    $resultsFile = Join-Path $runDir.FullName "results.json"
    if (Test-Path $resultsFile) {
        $completed++
    }
}

Write-Host "Found $completed completed runs in formamide_extended" -ForegroundColor Green
Write-Host ""

if ($completed -eq 0) {
    Write-Host "WARNING: No completed runs found. Make sure simulations have finished." -ForegroundColor Yellow
    exit 1
}

# Run analysis
Write-Host "Running complete Phase 2B analysis (will include formamide_extended)..."
Write-Host ""

python scripts/analyze_phase2b_complete.py `
    --input $INPUT_DIR `
    --output $OUTPUT_DIR

if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "ERROR: Analysis failed!" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "======================================================================"
Write-Host "FORMAMIDE EXTENDED ADDED TO ANALYSIS!" -ForegroundColor Green
Write-Host "Results updated in: $OUTPUT_DIR"
Write-Host "======================================================================"
Write-Host ""

