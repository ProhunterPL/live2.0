#!/bin/bash
# ONE-COMMAND DEPLOYMENT: Fix cluster deadlock and restart Phase 2B
# ==================================================================

set -e

echo "================================================================================"
echo "🚀 PHASE 2B CLUSTER DEADLOCK - EMERGENCY FIX DEPLOYMENT"
echo "================================================================================"
echo ""
echo "This script will:"
echo "  1. Kill stuck simulations (keep run_9)"
echo "  2. Apply code hotfix to disable cluster detection"
echo "  3. Wait for run_9 to complete or get stuck"
echo "  4. Restart 9 new simulations with safer config"
echo ""
echo "⚠️  WARNING: This will kill 7 stuck processes and modify Python code"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled. See docs/aws_test/CLUSTER_FIX_INSTRUCTIONS.md for manual steps."
    exit 1
fi

echo ""
echo "================================================================================"
echo "STEP 1: Kill Stuck Simulations"
echo "================================================================================"
echo ""

cd ~/live2.0

# Make scripts executable
chmod +x aws_test/scripts/kill_stuck_simulations.sh
chmod +x aws_test/scripts/apply_cluster_fix.sh
chmod +x aws_test/scripts/restart_phase2b_safe.sh

# Kill stuck processes
bash aws_test/scripts/kill_stuck_simulations.sh

echo ""
echo "Remaining processes:"
ps aux | grep "run_phase2_full.py" | grep -v grep || echo "  (none - all killed)"
echo ""

echo "================================================================================"
echo "STEP 2: Apply Code Hotfix"
echo "================================================================================"
echo ""

bash aws_test/scripts/apply_cluster_fix.sh

echo ""
echo "================================================================================"
echo "STEP 3: Monitor run_9"
echo "================================================================================"
echo ""

# Check if run_9 is still running
if ps aux | grep "run_phase2_full.py" | grep "run_9" | grep -v grep > /dev/null; then
    echo "✅ run_9 is still running. We'll monitor it..."
    echo ""
    echo "Current progress:"
    python3 ~/live2.0/aws_test/scripts/check_actual_progress.py | grep -A 20 "run_9"
    echo ""
    echo "⏳ Let's wait for run_9 to reach 100,000 steps or complete..."
    echo "   (This could take 1-2 hours depending on current step)"
    echo ""
    echo "You can monitor with:"
    echo "  watch -n 300 'python3 ~/live2.0/aws_test/scripts/check_actual_progress.py'"
    echo ""
    echo "Press ENTER when ready to start new simulations (or Ctrl+C to exit)"
    read -r
else
    echo "ℹ️  run_9 is not running (may have completed or been killed)"
fi

echo ""
echo "================================================================================"
echo "STEP 4: Start New Simulations (SAFER Config)"
echo "================================================================================"
echo ""

echo "Starting 9 new simulations with cluster detection disabled..."
echo ""

bash aws_test/scripts/restart_phase2b_safe.sh 9 110

echo ""
echo "================================================================================"
echo "✅ DEPLOYMENT COMPLETE"
echo "================================================================================"
echo ""
echo "What was done:"
echo "  ✅ Killed 7 stuck simulations"
echo "  ✅ Applied code hotfix for cluster detection"
echo "  ✅ Started 9 new simulations (runs 10-18)"
echo ""
echo "Current status:"
echo "  - run_1: Completed"
echo "  - run_9: Monitored/completed"
echo "  - runs 10-18: Running with safer config"
echo ""
echo "Monitor progress:"
echo "  python3 ~/live2.0/aws_test/scripts/check_actual_progress.py"
echo ""
echo "Check logs:"
echo "  tail -f ~/live2.0/logs/phase2b_safe/run_10.log"
echo ""
echo "Expected completion: ~60-90 minutes"
echo "================================================================================"

