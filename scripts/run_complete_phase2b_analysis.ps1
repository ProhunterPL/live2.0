# Complete Phase 2B Analysis Workflow (PowerShell)
# ================================================
#
# This script runs the complete analysis pipeline for Phase 2B results:
# 1. Fix molecules (extract from snapshots if catalog is empty)
# 2. Batch analysis
# 3. Complete analysis (autocatalysis, complexity)
#
# Usage:
#   .\scripts\run_complete_phase2b_analysis.ps1
#   .\scripts\run_complete_phase2b_analysis.ps1 -SkipFix  # Skip molecule fixing
#

param(
    [switch]$SkipFix = $false
)

$ErrorActionPreference = "Stop"

$BASE_DIR = "results/phase2b_additional"
$ANALYSIS_DIR = "analysis/phase2b_complete"
$PAPER_DIR = "paper/results_data"

Write-Host "======================================================================"
Write-Host "PHASE 2B COMPLETE ANALYSIS WORKFLOW"
Write-Host "======================================================================"
Write-Host ""

# Step 1: Fix molecules (if needed)
if (-not $SkipFix) {
    Write-Host "[STEP 1/3] Fixing molecules from snapshots..."
    Write-Host "------------------------------------------------------"
    python scripts/fix_run1_molecules.py --all --base-dir $BASE_DIR
    Write-Host ""
    
    # Verify
    Write-Host "Verifying extraction..."
    python -c "import json; from pathlib import Path; base = Path('$BASE_DIR'); count = sum(1 for f in base.rglob('results.json') if len(json.load(open(f)).get('molecules_detected', [])) > 0); print(f'Runs with molecules: {count}')"
    Write-Host ""
} else {
    Write-Host "[STEP 1/3] SKIPPED (-SkipFix flag)"
    Write-Host ""
}

# Step 2: Batch analysis
Write-Host "[STEP 2/3] Running batch analysis..."
Write-Host "------------------------------------------------------"
python scripts/analyze_phase2_batch.py `
    --input $BASE_DIR `
    --output $ANALYSIS_DIR `
    --recursive `
    --use-matcher
Write-Host ""

# Step 3: Complete analysis
Write-Host "[STEP 3/3] Running complete analysis..."
Write-Host "------------------------------------------------------"
python scripts/analyze_phase2b_complete.py `
    --input $BASE_DIR `
    --output $PAPER_DIR
Write-Host ""

Write-Host "======================================================================"
Write-Host "ANALYSIS COMPLETE!"
Write-Host "======================================================================"
Write-Host ""
Write-Host "Results saved to:"
Write-Host "  - Batch analysis: $ANALYSIS_DIR"
Write-Host "  - Complete analysis: $PAPER_DIR"
Write-Host ""
Write-Host "Next steps:"
Write-Host "  1. Review analysis reports"
Write-Host "  2. Generate figures (if scripts exist)"
Write-Host "  3. Generate tables (if scripts exist)"
Write-Host "  4. Fill LaTeX placeholders"
Write-Host ""

