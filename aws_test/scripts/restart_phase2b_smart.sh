#!/bin/bash
# Smart restart of Phase 2B miller_urey simulations
# ===================================================
# - Restarts incomplete runs from scratch (checkpoints don't work yet)
# - Max 5 parallel simulations
# - Priority: run_5 (95% done), then others by progress

set -e

PROJECT_ROOT="${PROJECT_ROOT:-$HOME/live2.0}"
RESULTS_DIR="$PROJECT_ROOT/results/phase2b_additional/miller_urey_extended"
CONFIG_FILE="$PROJECT_ROOT/aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml"
MAX_PARALLEL=5

cd "$PROJECT_ROOT" || exit 1

echo "================================================================================"
echo "🚀 SMART RESTART - PHASE 2B MILLER_UREY SIMULATIONS"
echo "================================================================================"
echo ""
echo "Settings:"
echo "  Max parallel: $MAX_PARALLEL"
echo "  Results dir: $RESULTS_DIR"
echo "  Config: $CONFIG_FILE"
echo ""

# Check config exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ Config file not found: $CONFIG_FILE"
    exit 1
fi

# Check current running simulations
RUNNING=$(ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | wc -l | tr -d ' ')
echo "📊 Currently running: $RUNNING miller_urey simulations"
echo ""

# List of runs to restart in priority order
# Priority: run_5 (95% done), run_3 (12%), run_10-18 (33%), run_4 (9%), run_9 (19%)
RUNS_TO_CHECK=(5 3 10 11 12 13 14 15 16 17 18 4 9)

STARTED=0
SKIPPED=0

for run_id in "${RUNS_TO_CHECK[@]}"; do
    RUN_DIR="$RESULTS_DIR/run_${run_id}"
    
    # Check if already completed
    if [ -f "$RUN_DIR/results.json" ]; then
        echo "⏭️  run_${run_id}: Already completed"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    # Check if already running
    if ps aux | grep "run_phase2_full.py" | grep "run_${run_id}" | grep -v grep > /dev/null; then
        echo "⏭️  run_${run_id}: Already running"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    # Check if we've reached the limit
    CURRENT_RUNNING=$(ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | wc -l | tr -d ' ')
    if [ "$CURRENT_RUNNING" -ge "$MAX_PARALLEL" ]; then
        echo "⏸️  Reached max parallel limit ($MAX_PARALLEL), stopping"
        break
    fi
    
    # Calculate seed (100 + run_id - 1)
    seed=$((100 + run_id - 1))
    
    # Show progress info if available
    PROGRESS_INFO=""
    if [ -f "$RUN_DIR/simulation.log" ]; then
        LAST_STEP=$(grep -oP "Step \K[0-9]+" "$RUN_DIR/simulation.log" | tail -1)
        if [ -n "$LAST_STEP" ]; then
            PROGRESS_PERCENT=$((LAST_STEP * 100 / 500000))
            PROGRESS_INFO=" (was at ${PROGRESS_PERCENT}%, restarting from 0)"
        fi
    fi
    
    echo "🚀 Starting run_${run_id} (seed ${seed})${PROGRESS_INFO}..."
    
    # Create/ensure output directory exists
    mkdir -p "$RUN_DIR"
    
    # Backup old log if exists
    if [ -f "$RUN_DIR/simulation.log" ]; then
        TIMESTAMP=$(date +%Y%m%d_%H%M%S)
        mv "$RUN_DIR/simulation.log" "$RUN_DIR/simulation.log.backup_${TIMESTAMP}"
    fi
    
    # Start simulation in background
    nohup python3 "$PROJECT_ROOT/scripts/run_phase2_full.py" \
        --config "$CONFIG_FILE" \
        --output "$RUN_DIR" \
        --seed "$seed" \
        --steps 500000 \
        --force-cpu \
        > "$RUN_DIR/simulation.log" 2>&1 &
    
    PID=$!
    echo "   ✅ Started with PID $PID"
    
    STARTED=$((STARTED + 1))
    
    # Small delay to avoid overwhelming the system
    sleep 3
done

echo ""
echo "================================================================================"
echo "📊 SUMMARY"
echo "================================================================================"
echo "  Started: $STARTED simulations"
echo "  Skipped: $SKIPPED simulations (completed or running)"
echo ""

# Show current status
TOTAL_RUNNING=$(ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | wc -l | tr -d ' ')
echo "  Total running now: $TOTAL_RUNNING simulations"
echo ""

# Show which runs are running
if [ "$TOTAL_RUNNING" -gt 0 ]; then
    echo "📋 Currently running:"
    ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | while read line; do
        PID=$(echo "$line" | awk '{print $2}')
        RUN=$(echo "$line" | grep -oP "run_\K[0-9]+")
        echo "   - run_${RUN} (PID $PID)"
    done
    echo ""
fi

# List incomplete runs
INCOMPLETE=0
echo "📝 Still incomplete:"
for i in {1..18}; do
    if [ ! -f "$RESULTS_DIR/run_${i}/results.json" ]; then
        if ! ps aux | grep "run_phase2_full.py" | grep "run_${i}" | grep -v grep > /dev/null; then
            echo "   - run_${i}"
            INCOMPLETE=$((INCOMPLETE + 1))
        fi
    fi
done

if [ "$INCOMPLETE" -eq 0 ]; then
    echo "   (none - all runs completed or running!)"
fi

echo ""
echo "================================================================================"
echo "✅ DONE"
echo "================================================================================"
echo ""
echo "Monitor progress with:"
echo "  python3 aws_test/scripts/check_actual_progress.py"
echo "  ps aux | grep run_phase2_full.py | grep -v grep"
echo "  tail -f $RESULTS_DIR/run_5/simulation.log"
echo ""
echo "To start more (if < $MAX_PARALLEL running):"
echo "  bash $0"
echo ""

