#!/bin/bash
# Kill stuck simulations (keep only run_9 which is working)
# =========================================================

echo "🔍 Checking for stuck simulations..."
echo ""

# Get all simulation PIDs
PIDS=$(ps aux | grep "run_phase2_full.py" | grep -v grep | awk '{print $2}')

for PID in $PIDS; do
    # Get the run number
    CMD=$(ps -p $PID -o args= | grep -o "run_[0-9]*" | head -1)
    
    # Skip run_9 (it's working)
    if [[ "$CMD" == "run_9" ]]; then
        echo "✅ Keeping $CMD (PID $PID) - it's progressing normally"
        continue
    fi
    
    # Check last log modification time for this run
    OUTPUT_DIR=$(ps -p $PID -o args= | grep -o "\--output [^ ]*" | awk '{print $2}')
    LOG_FILE="$OUTPUT_DIR/simulation.log"
    
    if [ -f "$LOG_FILE" ]; then
        # Get minutes since last modification
        AGE_MIN=$(echo "($(date +%s) - $(stat -c %Y "$LOG_FILE")) / 60" | bc)
        
        if [ $AGE_MIN -gt 60 ]; then
            echo "🚫 Killing $CMD (PID $PID) - stuck for $AGE_MIN minutes"
            kill -9 $PID
        else
            echo "✅ Keeping $CMD (PID $PID) - log updated $AGE_MIN min ago"
        fi
    else
        echo "❓ Unknown state for PID $PID - check manually"
    fi
done

echo ""
echo "Done! Check remaining processes:"
ps aux | grep "run_phase2_full.py" | grep -v grep

