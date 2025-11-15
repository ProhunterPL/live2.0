#!/bin/bash
# Quick monitor for Phase 2B progress
# ====================================

PROJECT_ROOT="${PROJECT_ROOT:-$HOME/live2.0}"
RESULTS_DIR="$PROJECT_ROOT/results/phase2b_additional/miller_urey_extended"

echo "================================================================================"
echo "📊 PHASE 2B QUICK STATUS - $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================================"
echo ""

# Count running
RUNNING=$(ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | wc -l | tr -d ' ')
echo "🏃 Running: $RUNNING simulations"
echo ""

# Count completed
COMPLETED=0
for i in {1..18}; do
    if [ -f "$RESULTS_DIR/run_${i}/results.json" ]; then
        COMPLETED=$((COMPLETED + 1))
    fi
done
echo "✅ Completed: $COMPLETED / 18"
echo ""

# Show running processes with progress
if [ "$RUNNING" -gt 0 ]; then
    echo "📋 Running processes:"
    for i in {1..18}; do
        if ps aux | grep "run_phase2_full.py" | grep "run_${i}" | grep -v grep > /dev/null; then
            PID=$(ps aux | grep "run_phase2_full.py" | grep "run_${i}" | grep -v grep | awk '{print $2}')
            
            # Get latest step from log
            LOG_FILE="$RESULTS_DIR/run_${i}/simulation.log"
            if [ -f "$LOG_FILE" ]; then
                LAST_STEP=$(grep -oP "Step \K[0-9]+" "$LOG_FILE" | tail -1)
                if [ -n "$LAST_STEP" ]; then
                    PROGRESS_PERCENT=$((LAST_STEP * 100 / 500000))
                    echo "   run_${i} (PID $PID): Step ${LAST_STEP}/500K (${PROGRESS_PERCENT}%)"
                else
                    echo "   run_${i} (PID $PID): Starting..."
                fi
            else
                echo "   run_${i} (PID $PID): No log yet"
            fi
        fi
    done
fi

echo ""
echo "================================================================================"
echo "🎯 Progress: $COMPLETED completed + $RUNNING running = $((COMPLETED + RUNNING)) / 18"
echo "================================================================================"

