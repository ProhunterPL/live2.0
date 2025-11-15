#!/bin/bash
# Start Phase 2B simulations NOW - No questions asked
# ====================================================
# This script will start 5 simulations immediately without asking

set -e

PROJECT_ROOT="$HOME/live2.0"
cd "$PROJECT_ROOT"

RESULTS_DIR="$PROJECT_ROOT/results/phase2b_additional/miller_urey_extended"
CONFIG_FILE="$PROJECT_ROOT/aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml"
SCRIPT_FILE="$PROJECT_ROOT/scripts/run_phase2_full.py"

echo "================================================================================"
echo "🚀 STARTING PHASE 2B SIMULATIONS NOW"
echo "================================================================================"
echo ""

# Check current status
RUNNING=$(ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | wc -l | tr -d ' ')
echo "Currently running: $RUNNING simulations"
echo ""

# Priority runs to start
PRIORITY_RUNS=(3 4 5 2 10)
MAX_TO_START=5
STARTED=0

for run_id in "${PRIORITY_RUNS[@]}"; do
    # Check if already completed
    if [ -f "$RESULTS_DIR/run_${run_id}/results.json" ]; then
        echo "⏭️  run_${run_id} already completed, skipping"
        continue
    fi
    
    # Check if already running
    if ps aux | grep "run_phase2_full.py" | grep "run_${run_id}" | grep -v grep > /dev/null 2>&1; then
        echo "⏭️  run_${run_id} already running, skipping"
        continue
    fi
    
    # Calculate seed (run_id + 99)
    seed=$((run_id + 99))
    
    echo "🚀 Starting run_${run_id} (seed ${seed})..."
    
    OUTPUT_DIR="$RESULTS_DIR/run_${run_id}"
    mkdir -p "$OUTPUT_DIR"
    
    # Start simulation
    nohup python3 "$SCRIPT_FILE" \
        --config "$CONFIG_FILE" \
        --output "$OUTPUT_DIR" \
        --seed "$seed" \
        --steps 500000 \
        --force-cpu \
        >> "$OUTPUT_DIR/simulation.log" 2>&1 &
    
    PID=$!
    echo "   ✅ Started with PID: $PID"
    
    ((STARTED++))
    
    # Don't exceed max parallel
    if [ "$STARTED" -ge "$MAX_TO_START" ]; then
        break
    fi
    
    # Small delay between starts
    sleep 2
done

echo ""
echo "================================================================================"
echo "✅ STARTED $STARTED SIMULATIONS"
echo "================================================================================"
echo ""

# Wait a moment and verify
sleep 5
NEW_RUNNING=$(ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | wc -l | tr -d ' ')

echo "📊 Final Status:"
echo "  Running: $NEW_RUNNING simulations"
echo ""

if [ "$NEW_RUNNING" -gt 0 ]; then
    echo "✅ Simulations are running!"
    echo ""
    echo "Processes:"
    ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | awk '{print "  PID " $2 ": run_" $NF}'
    echo ""
    echo "Monitor with:"
    echo "  ps aux | grep run_phase2_full.py | grep -v grep"
    echo "  python3 aws_test/scripts/check_actual_progress.py"
    echo "  tail -f $RESULTS_DIR/run_3/simulation.log"
else
    echo "❌ No simulations running! Check logs:"
    echo "  tail -50 $RESULTS_DIR/run_3/simulation.log"
    echo "  tail -50 $RESULTS_DIR/run_4/simulation.log"
fi

echo ""

