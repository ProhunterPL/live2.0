#!/bin/bash
# Quick Start Phase 2B - Start 5 simulations immediately
# =========================================================

set -e

PROJECT_ROOT="$HOME/live2.0"
cd "$PROJECT_ROOT"

RESULTS_DIR="$PROJECT_ROOT/results/phase2b_additional/miller_urey_extended"
CONFIG_FILE="$PROJECT_ROOT/aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml"
SCRIPT_FILE="$PROJECT_ROOT/scripts/run_phase2_full.py"

echo "================================================================================"
echo "🚀 QUICK START - Phase 2B Simulations"
echo "================================================================================"
echo ""

# Check how many are already running
RUNNING=$(ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | wc -l | tr -d ' ')
echo "Currently running: $RUNNING simulations"
echo ""

MAX_TO_START=5
STARTED=0

# Priority order: runs that need to be restarted
PRIORITY_RUNS=(3 4 5 2 10)

for run_id in "${PRIORITY_RUNS[@]}"; do
    # Check if already completed
    if [ -f "$RESULTS_DIR/run_${run_id}/results.json" ]; then
        echo "⏭️  run_${run_id} already completed, skipping"
        continue
    fi
    
    # Check if already running
    if ps aux | grep "run_phase2_full.py" | grep "run_${run_id}" | grep -v grep > /dev/null; then
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

echo "📊 Status:"
echo "  Running: $NEW_RUNNING simulations"
echo ""
echo "Monitor with:"
echo "  ps aux | grep run_phase2_full.py | grep -v grep"
echo "  python3 aws_test/scripts/check_actual_progress.py"
echo ""
echo "Check logs:"
echo "  tail -f $RESULTS_DIR/run_3/simulation.log"
echo ""

