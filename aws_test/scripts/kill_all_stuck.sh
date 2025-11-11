#!/bin/bash
# Kill ALL stuck simulations (hydrothermal + miller_urey)
# Based on log file age threshold
# ========================================================

STUCK_THRESHOLD_MIN=${1:-60}  # Default: 60 minutes
RESULTS_DIR="${2:-$HOME/live2.0/results/phase2b_additional}"

echo "🔍 Checking for stuck simulations (threshold: ${STUCK_THRESHOLD_MIN} min)..."
echo "   Results dir: $RESULTS_DIR"
echo ""

# Get all simulation PIDs
PIDS=$(ps aux | grep "run_phase2_full.py" | grep -v grep | awk '{print $2}')

if [ -z "$PIDS" ]; then
    echo "✅ No simulation processes found"
    exit 0
fi

KILLED=0
KEPT=0
STUCK_RUNS=()

for PID in $PIDS; do
    # Get the full command line
    CMD=$(ps -p $PID -o args= 2>/dev/null)
    
    if [ -z "$CMD" ]; then
        continue
    fi
    
    # Extract output directory
    OUTPUT_DIR=$(echo "$CMD" | grep -o "\--output [^ ]*" | awk '{print $2}')
    
    if [ -z "$OUTPUT_DIR" ]; then
        # Try to extract from results directory pattern
        OUTPUT_DIR=$(echo "$CMD" | grep -o "results/[^ ]*" | head -1)
    fi
    
    if [ -z "$OUTPUT_DIR" ] || [ ! -d "$OUTPUT_DIR" ]; then
        # Try to find from PID's working directory
        OUTPUT_DIR=$(pwdx $PID 2>/dev/null | awk '{print $2}' || echo "")
        if [ -z "$OUTPUT_DIR" ] || [ ! -d "$OUTPUT_DIR" ]; then
            echo "⚠️  PID $PID: Cannot determine output directory, checking process state..."
            STATE=$(ps -p $PID -o state= 2>/dev/null)
            if [ "$STATE" = "D" ] || [ "$STATE" = "Z" ]; then
                echo "🚫 Killing PID $PID - process in bad state ($STATE)"
                kill -9 $PID 2>/dev/null
                KILLED=$((KILLED + 1))
            else
                echo "✅ Keeping PID $PID - cannot verify but process state OK"
                KEPT=$((KEPT + 1))
            fi
            continue
        fi
    fi
    
    LOG_FILE="$OUTPUT_DIR/simulation.log"
    RUN_NAME=$(basename "$OUTPUT_DIR")
    SCENARIO=$(basename $(dirname "$OUTPUT_DIR"))
    
    if [ -f "$LOG_FILE" ]; then
        # Get minutes since last modification
        if command -v stat >/dev/null 2>&1; then
            # Linux stat
            if stat -c %Y "$LOG_FILE" >/dev/null 2>&1; then
                AGE_MIN=$(echo "($(date +%s) - $(stat -c %Y "$LOG_FILE")) / 60" | bc)
            else
                # macOS stat
                AGE_MIN=$(echo "($(date +%s) - $(stat -f %m "$LOG_FILE")) / 60" | bc)
            fi
        else
            AGE_MIN=999999
        fi
        
        if [ -z "$AGE_MIN" ] || [ "$AGE_MIN" -gt "$STUCK_THRESHOLD_MIN" ]; then
            echo "🚫 Killing $SCENARIO/$RUN_NAME (PID $PID) - stuck for ${AGE_MIN} min"
            kill -9 $PID 2>/dev/null
            KILLED=$((KILLED + 1))
            STUCK_RUNS+=("$SCENARIO/$RUN_NAME")
        else
            echo "✅ Keeping $SCENARIO/$RUN_NAME (PID $PID) - log updated $AGE_MIN min ago"
            KEPT=$((KEPT + 1))
        fi
    else
        # No log file - might be very new or broken
        echo "⚠️  PID $PID ($SCENARIO/$RUN_NAME): No log file found, checking process state..."
        STATE=$(ps -p $PID -o state= 2>/dev/null)
        if [ "$STATE" = "D" ] || [ "$STATE" = "Z" ]; then
            echo "🚫 Killing $SCENARIO/$RUN_NAME (PID $PID) - process in bad state ($STATE)"
            kill -9 $PID 2>/dev/null
            KILLED=$((KILLED + 1))
            STUCK_RUNS+=("$SCENARIO/$RUN_NAME")
        else
            echo "✅ Keeping $SCENARIO/$RUN_NAME (PID $PID) - process state OK"
            KEPT=$((KEPT + 1))
        fi
    fi
done

echo ""
echo "📊 Summary:"
echo "   Killed: $KILLED stuck processes"
echo "   Kept: $KEPT active processes"
echo ""

if [ ${#STUCK_RUNS[@]} -gt 0 ]; then
    echo "📋 Stuck runs that were killed:"
    for run in "${STUCK_RUNS[@]}"; do
        echo "   - $run"
    done
    echo ""
    echo "💡 To restart these runs:"
    echo "   bash aws_test/scripts/restart_from_checkpoint.sh $RESULTS_DIR <scenario>"
    echo ""
fi

echo "Remaining processes:"
REMAINING=$(ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l | tr -d ' ')
echo "   Count: $REMAINING"
if [ "$REMAINING" -gt 0 ]; then
    ps aux | grep "run_phase2_full.py" | grep -v grep
fi

