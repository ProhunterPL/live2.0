#!/bin/bash
# Kill stuck hydrothermal_extended simulations (runs 1-10)
# =========================================================

echo "🔍 Checking for stuck hydrothermal_extended simulations..."
echo ""

RESULTS_DIR="$HOME/live2.0/results/phase2b_additional"
STUCK_THRESHOLD_MIN=60  # Consider stuck if no activity for 60 minutes

# Get all simulation PIDs
PIDS=$(ps aux | grep "run_phase2_full.py" | grep -v grep | awk '{print $2}')

if [ -z "$PIDS" ]; then
    echo "❌ No simulation processes found"
    exit 0
fi

KILLED=0
KEPT=0

for PID in $PIDS; do
    # Get the full command line
    CMD=$(ps -p $PID -o args=)
    
    # Check if this is a hydrothermal_extended run
    if [[ "$CMD" != *"hydrothermal_extended"* ]]; then
        continue
    fi
    
    # Extract run number and output directory
    RUN_MATCH=$(echo "$CMD" | grep -o "run_[0-9]*")
    OUTPUT_DIR=$(echo "$CMD" | grep -o "\--output [^ ]*" | awk '{print $2}')
    
    if [ -z "$OUTPUT_DIR" ]; then
        # Try to extract from results directory pattern
        OUTPUT_DIR=$(echo "$CMD" | grep -o "results/[^ ]*" | head -1)
    fi
    
    if [ -z "$OUTPUT_DIR" ] || [ ! -d "$OUTPUT_DIR" ]; then
        echo "⚠️  PID $PID: Cannot determine output directory, skipping"
        continue
    fi
    
    LOG_FILE="$OUTPUT_DIR/simulation.log"
    
    if [ -f "$LOG_FILE" ]; then
        # Get minutes since last modification
        if command -v stat >/dev/null 2>&1; then
            # Linux stat
            AGE_MIN=$(echo "($(date +%s) - $(stat -c %Y "$LOG_FILE")) / 60" | bc)
        else
            # macOS stat
            AGE_MIN=$(echo "($(date +%s) - $(stat -f %m "$LOG_FILE")) / 60" | bc)
        fi
        
        if [ -z "$AGE_MIN" ] || [ "$AGE_MIN" -gt "$STUCK_THRESHOLD_MIN" ]; then
            echo "🚫 Killing $RUN_MATCH (PID $PID) - stuck for ${AGE_MIN} minutes"
            kill -9 $PID 2>/dev/null
            KILLED=$((KILLED + 1))
        else
            echo "✅ Keeping $RUN_MATCH (PID $PID) - log updated $AGE_MIN min ago"
            KEPT=$((KEPT + 1))
        fi
    else
        # No log file - might be very new or broken
        echo "⚠️  PID $PID ($RUN_MATCH): No log file found, checking process state..."
        STATE=$(ps -p $PID -o state= 2>/dev/null)
        if [ "$STATE" = "D" ] || [ "$STATE" = "Z" ]; then
            echo "🚫 Killing $RUN_MATCH (PID $PID) - process in bad state ($STATE)"
            kill -9 $PID 2>/dev/null
            KILLED=$((KILLED + 1))
        else
            echo "✅ Keeping $RUN_MATCH (PID $PID) - process state OK"
            KEPT=$((KEPT + 1))
        fi
    fi
done

echo ""
echo "📊 Summary:"
echo "   Killed: $KILLED stuck hydrothermal processes"
echo "   Kept: $KEPT active hydrothermal processes"
echo ""
echo "Remaining processes:"
ps aux | grep "run_phase2_full.py" | grep -v grep || echo "   (none)"

