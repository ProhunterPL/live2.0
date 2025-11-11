#!/bin/bash
# Complete fix for Phase 2B issues
# ==================================
# 1. Kill stuck hydrothermal simulations
# 2. Limit parallel simulations to 4
# 3. Apply hydrothermal fix (already done in config)
# 4. Optionally restart stuck simulations

set -e

echo "================================================================================"
echo "🔧 PHASE 2B COMPLETE FIX"
echo "================================================================================"
echo ""
echo "This script will:"
echo "  1. Kill stuck hydrothermal_extended simulations (runs 1-10)"
echo "  2. Limit parallel simulations to max 4"
echo "  3. Verify hydrothermal fix is applied"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 1
fi

cd ~/live2.0 || exit 1

# Make scripts executable
chmod +x aws_test/scripts/*.sh

echo ""
echo "================================================================================"
echo "STEP 1: Kill Stuck Hydrothermal Simulations"
echo "================================================================================"
echo ""

bash aws_test/scripts/kill_stuck_hydrothermal.sh

echo ""
echo "================================================================================"
echo "STEP 2: Limit Parallel Simulations to 4"
echo "================================================================================"
echo ""

bash aws_test/scripts/limit_parallel_simulations.sh 4

echo ""
echo "================================================================================"
echo "STEP 3: Verify Hydrothermal Fix"
echo "================================================================================"
echo ""

CONFIG_FILE="aws_test/configs/phase2_hydrothermal_extended_SUPER_FAST.yaml"
if grep -q "cluster_check_interval: 999999999" "$CONFIG_FILE"; then
    echo "✅ Hydrothermal fix applied: cluster_check_interval = 999999999 (disabled)"
else
    echo "❌ Hydrothermal fix NOT applied!"
    echo "   Expected: cluster_check_interval: 999999999"
    echo "   Please check: $CONFIG_FILE"
fi

echo ""
echo "================================================================================"
echo "STEP 4: Current Status"
echo "================================================================================"
echo ""

RUNNING=$(ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l | tr -d ' ')
echo "📊 Currently running simulations: $RUNNING"

if [ "$RUNNING" -gt 4 ]; then
    echo "⚠️  WARNING: $RUNNING simulations running (should be <= 4)"
    echo "   Run: bash aws_test/scripts/limit_parallel_simulations.sh 4"
else
    echo "✅ Within limit ($RUNNING <= 4)"
fi

echo ""
echo "================================================================================"
echo "✅ FIX COMPLETE"
echo "================================================================================"
echo ""
echo "Next steps:"
echo "  1. Monitor progress: python3 aws_test/scripts/check_actual_progress.py"
echo "  2. If needed, restart stuck simulations:"
echo "     bash aws_test/scripts/restart_from_checkpoint.sh ~/live2.0/results/phase2b_additional hydrothermal_extended"
echo "  3. Ensure max 4 parallel: bash aws_test/scripts/limit_parallel_simulations.sh 4"
echo ""

