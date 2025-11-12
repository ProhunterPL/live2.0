#!/bin/bash
# Kill Stuck Phase 2B Simulations - Emergency Fix
# =================================================
# This script kills simulations that are deadlocked in cluster detection
# 
# Based on progress report from 2025-11-12 08:16:53:
# - All 10 hydrothermal runs stuck at step 1000 (41+ hours)
# - miller_urey run_9 stuck at 97K (41+ hours)
# - Remaining miller_urey runs (2-8, 10-18) are ACTIVE - DO NOT KILL

set -e

echo "================================================================================"
echo "🚨 Kill Stuck Phase 2B Simulations"
echo "================================================================================"
echo ""
echo "This script will kill:"
echo "  - All 10 hydrothermal_extended runs (stuck at step ~1000)"
echo "  - miller_urey_extended run_9 (stuck at step ~97K)"
echo ""
echo "Will KEEP running (actively progressing):"
echo "  - miller_urey runs 1-8, 10-18"
echo ""
echo "⚠️  WARNING: This will kill processes. Make sure you're on the correct machine!"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 1
fi

echo ""
echo "================================================================================"
echo "Step 1: Kill All Hydrothermal Runs (10 processes)"
echo "================================================================================"
echo ""

# Find and kill all hydrothermal processes
HYDRO_PIDS=$(ps aux | grep "run_phase2_full.py" | grep "hydrothermal_extended" | grep -v grep | awk '{print $2}')

if [ -z "$HYDRO_PIDS" ]; then
    echo "✅ No hydrothermal processes found (already killed or not running)"
else
    HYDRO_COUNT=$(echo "$HYDRO_PIDS" | wc -l)
    echo "🔍 Found $HYDRO_COUNT hydrothermal processes:"
    ps aux | grep "run_phase2_full.py" | grep "hydrothermal_extended" | grep -v grep
    echo ""
    echo "🚫 Killing hydrothermal processes..."
    echo "$HYDRO_PIDS" | xargs kill -9 2>/dev/null || true
    echo "✅ Killed $HYDRO_COUNT hydrothermal processes"
fi

echo ""
echo "================================================================================"
echo "Step 2: Kill miller_urey run_9 (Stuck)"
echo "================================================================================"
echo ""

# Find and kill miller_urey run_9
RUN9_PID=$(ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended/run_9" | grep -v grep | awk '{print $2}')

if [ -z "$RUN9_PID" ]; then
    echo "✅ miller_urey run_9 not found (already killed or not running)"
else
    echo "🔍 Found miller_urey run_9 process:"
    ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended/run_9" | grep -v grep
    echo ""
    echo "🚫 Killing run_9..."
    kill -9 $RUN9_PID 2>/dev/null || true
    echo "✅ Killed run_9"
fi

echo ""
echo "================================================================================"
echo "✅ CLEANUP COMPLETE"
echo "================================================================================"
echo ""

# Show remaining processes
REMAINING=$(ps aux | grep "run_phase2_full.py" | grep -v grep)

if [ -z "$REMAINING" ]; then
    echo "⚠️  No simulations running (all killed)"
else
    REMAINING_COUNT=$(echo "$REMAINING" | wc -l)
    echo "✅ $REMAINING_COUNT simulations still running (should be ~16 miller_urey):"
    echo ""
    echo "$REMAINING"
fi

echo ""
echo "================================================================================"
echo "Next Steps:"
echo "================================================================================"
echo ""
echo "1. Verify cluster detection fix is deployed:"
echo "   grep 'cluster_interval = getattr' ~/live2.0/backend/sim/core/stepper.py"
echo ""
echo "2. Check config files have cluster_check_interval: 999999999:"
echo "   grep 'cluster_check_interval' ~/live2.0/aws_test/configs/*SUPER_FAST.yaml"
echo ""
echo "3. Monitor remaining simulations:"
echo "   python3 ~/live2.0/aws_test/scripts/check_actual_progress.py"
echo ""
echo "4. Wait for remaining 16 miller_urey runs to complete (~12-36 hours)"
echo ""
echo "5. Optionally restart hydrothermal (now with fix):"
echo "   python3 ~/live2.0/aws_test/scripts/run_phase2b_additional.py --scenario hydrothermal_extended"
echo ""
echo "================================================================================"

