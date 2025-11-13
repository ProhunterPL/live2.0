#!/bin/bash
# Kill Only Stuck Processes - Safe version that preserves active runs
# ===================================================================

echo "================================================================================"
echo "🔍 IDENTIFYING STUCK PROCESSES"
echo "================================================================================"
echo ""

RESULTS_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended"
STUCK_PIDS=()

# Check each known PID
declare -A PROCESS_MAP
PROCESS_MAP[73085]="run_6"
PROCESS_MAP[76917]="run_5"
PROCESS_MAP[79713]="run_7"
PROCESS_MAP[79914]="run_8"
PROCESS_MAP[80149]="run_3"

for PID in 73085 76917 79713 79914 80149; do
    RUN=${PROCESS_MAP[$PID]}
    
    # Check if process exists
    if ! ps -p $PID > /dev/null 2>&1; then
        continue
    fi
    
    STATE=$(ps -p $PID -o state= 2>/dev/null | tr -d ' ')
    CPU=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
    LOG_FILE="$RESULTS_DIR/$RUN/simulation.log"
    
    if [ ! -f "$LOG_FILE" ]; then
        echo "⚠️  PID $PID ($RUN): No log file - skipping"
        continue
    fi
    
    # Get log age
    LAST_MOD=$(stat -c %Y "$LOG_FILE" 2>/dev/null)
    NOW=$(date +%s)
    AGE_MIN=$(( (NOW - LAST_MOD) / 60 ))
    AGE_HOURS=$((AGE_MIN / 60))
    
    # Get last step
    LAST_STEP=$(tail -100 "$LOG_FILE" 2>/dev/null | grep -E "Step [0-9]+" | tail -1 | grep -oE "Step [0-9]+" | head -1)
    
    # Criteria for stuck: Sleeping + High CPU + Old log (>60min)
    if [ "$STATE" = "R" ]; then
        echo "✅ PID $PID ($RUN): ACTIVE - State R, CPU ${CPU}%, log ${AGE_MIN}min old"
    elif [ "$STATE" = "S" ] && [ "${CPU%.*}" -gt 100 ] && [ "$AGE_MIN" -gt 60 ]; then
        echo "🚨 PID $PID ($RUN): STUCK - State S, CPU ${CPU}%, log ${AGE_HOURS}h ${AGE_MIN}min old, $LAST_STEP"
        STUCK_PIDS+=($PID)
    elif [ "$STATE" = "S" ] && [ "${CPU%.*}" -gt 100 ] && [ "$AGE_MIN" -le 60 ]; then
        echo "✅ PID $PID ($RUN): ACTIVE (sleeping) - State S, CPU ${CPU}%, log ${AGE_MIN}min old"
    else
        echo "❓ PID $PID ($RUN): UNKNOWN - State $STATE, CPU ${CPU}%, log ${AGE_MIN}min old"
    fi
done

echo ""
echo "================================================================================"
echo "📊 SUMMARY"
echo "================================================================================"
echo ""

if [ ${#STUCK_PIDS[@]} -eq 0 ]; then
    echo "✅ No stuck processes found. All processes appear to be active."
    exit 0
fi

echo "Found ${#STUCK_PIDS[@]} stuck process(es):"
for PID in "${STUCK_PIDS[@]}"; do
    RUN=${PROCESS_MAP[$PID]}
    echo "  - PID $PID ($RUN)"
done

echo ""
echo "⚠️  WARNING: This will kill ${#STUCK_PIDS[@]} stuck process(es)."
echo "   Active processes will be preserved."
read -p "Continue? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 1
fi

echo ""
echo "================================================================================"
echo "🚫 KILLING STUCK PROCESSES"
echo "================================================================================"
echo ""

KILLED=0
for PID in "${STUCK_PIDS[@]}"; do
    RUN=${PROCESS_MAP[$PID]}
    echo "Killing PID $PID ($RUN)..."
    kill -9 $PID 2>/dev/null || true
    KILLED=$((KILLED + 1))
done

echo ""
echo "✅ Killed $KILLED stuck process(es)"
echo ""

REMAINING=$(ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l | tr -d ' ')
echo "Remaining processes: $REMAINING"
echo ""

if [ "$REMAINING" -gt 0 ]; then
    echo "Active processes:"
    ps aux | grep "run_phase2_full.py" | grep -v grep | awk '{printf "  PID %s: %s (CPU: %s%%, State: %s)\n", $2, $NF, $3, $8}'
fi

echo ""
echo "================================================================================"
echo "💡 NEXT STEPS"
echo "================================================================================"
echo ""
echo "1. Verify remaining processes are active:"
echo "   bash aws_test/scripts/analyze_current_situation.sh"
echo ""
echo "2. Restart killed runs (maintaining ≤5 total processes):"
echo "   bash aws_test/scripts/auto_queue_restart.sh"
echo ""
echo "3. Monitor progress:"
echo "   python3 aws_test/scripts/check_actual_progress.py"
echo ""

