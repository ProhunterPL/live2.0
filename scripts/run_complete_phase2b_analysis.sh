#!/bin/bash
# Complete Phase 2B Analysis Workflow
# ===================================
#
# This script runs the complete analysis pipeline for Phase 2B results:
# 1. Fix molecules (extract from snapshots if catalog is empty)
# 2. Batch analysis
# 3. Complete analysis (autocatalysis, complexity)
#
# Usage:
#   bash scripts/run_complete_phase2b_analysis.sh
#   bash scripts/run_complete_phase2b_analysis.sh --skip-fix  # Skip molecule fixing
#

set -e  # Exit on error

BASE_DIR="results/phase2b_additional"
ANALYSIS_DIR="analysis/phase2b_complete"
PAPER_DIR="paper/results_data"

SKIP_FIX=false
if [[ "$1" == "--skip-fix" ]]; then
    SKIP_FIX=true
fi

echo "======================================================================"
echo "PHASE 2B COMPLETE ANALYSIS WORKFLOW"
echo "======================================================================"
echo ""

# Step 1: Fix molecules (if needed)
if [ "$SKIP_FIX" = false ]; then
    echo "[STEP 1/3] Fixing molecules from snapshots..."
    echo "------------------------------------------------------"
    python scripts/fix_run1_molecules.py --all --base-dir "$BASE_DIR"
    echo ""
    
    # Verify
    echo "Verifying extraction..."
    python -c "
import json
from pathlib import Path
base = Path('$BASE_DIR')
count = 0
for results_file in base.rglob('results.json'):
    try:
        with open(results_file) as f:
            d = json.load(f)
            mols = len(d.get('molecules_detected', []))
            if mols > 0:
                count += 1
    except:
        pass
print(f'Runs with molecules: {count}')
"
    echo ""
else
    echo "[STEP 1/3] SKIPPED (--skip-fix flag)"
    echo ""
fi

# Step 2: Batch analysis
echo "[STEP 2/3] Running batch analysis..."
echo "------------------------------------------------------"
python scripts/analyze_phase2_batch.py \
    --input "$BASE_DIR" \
    --output "$ANALYSIS_DIR" \
    --recursive \
    --use-matcher
echo ""

# Step 3: Complete analysis
echo "[STEP 3/3] Running complete analysis..."
echo "------------------------------------------------------"
python scripts/analyze_phase2b_complete.py \
    --input "$BASE_DIR" \
    --output "$PAPER_DIR"
echo ""

echo "======================================================================"
echo "ANALYSIS COMPLETE!"
echo "======================================================================"
echo ""
echo "Results saved to:"
echo "  - Batch analysis: $ANALYSIS_DIR"
echo "  - Complete analysis: $PAPER_DIR"
echo ""
echo "Next steps:"
echo "  1. Review analysis reports"
echo "  2. Generate figures (if scripts exist)"
echo "  3. Generate tables (if scripts exist)"
echo "  4. Fill LaTeX placeholders"
echo ""

