#!/bin/bash
# Identify which processes correspond to which runs
# ==================================================

echo "🔍 Identifying running simulation processes..."
echo ""

# Get all simulation PIDs
PIDS=$(ps aux | grep "run_phase2_full.py" | grep -v grep | awk '{print $2}')

if [ -z "$PIDS" ]; then
    echo "❌ No simulation processes found"
    exit 0
fi

echo "Found $(echo "$PIDS" | wc -l | tr -d ' ') processes:"
echo ""

for PID in $PIDS; do
    # Get the full command line
    CMD=$(ps -p $PID -o args= 2>/dev/null)
    
    if [ -z "$CMD" ]; then
        continue
    fi
    
    # Extract info
    OUTPUT_DIR=$(echo "$CMD" | grep -o "\--output [^ ]*" | awk '{print $2}')
    if [ -z "$OUTPUT_DIR" ]; then
        OUTPUT_DIR=$(echo "$CMD" | grep -o "results/[^ ]*" | head -1)
    fi
    
    RUN_NAME="unknown"
    SCENARIO="unknown"
    
    if [ -n "$OUTPUT_DIR" ] && [ -d "$OUTPUT_DIR" ]; then
        RUN_NAME=$(basename "$OUTPUT_DIR")
        SCENARIO=$(basename $(dirname "$OUTPUT_DIR"))
        
        # Check log age
        LOG_FILE="$OUTPUT_DIR/simulation.log"
        if [ -f "$LOG_FILE" ]; then
            if command -v stat >/dev/null 2>&1; then
                if stat -c %Y "$LOG_FILE" >/dev/null 2>&1; then
                    AGE_MIN=$(echo "($(date +%s) - $(stat -c %Y "$LOG_FILE")) / 60" | bc)
                else
                    AGE_MIN=$(echo "($(date +%s) - $(stat -f %m "$LOG_FILE")) / 60" | bc)
                fi
            else
                AGE_MIN="unknown"
            fi
        else
            AGE_MIN="no log"
        fi
    else
        AGE_MIN="unknown dir"
    fi
    
    # Get process state
    STATE=$(ps -p $PID -o state= 2>/dev/null)
    CPU=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
    MEM=$(ps -p $PID -o %mem= 2>/dev/null | tr -d ' ')
    
    STATUS=""
    if [ "$AGE_MIN" != "unknown" ] && [ "$AGE_MIN" != "no log" ] && [ "$AGE_MIN" != "unknown dir" ]; then
        if [ "$AGE_MIN" -gt 60 ]; then
            STATUS="⚠️  STUCK"
        else
            STATUS="✅ Active"
        fi
    else
        STATUS="❓ Unknown"
    fi
    
    echo "PID $PID: $SCENARIO/$RUN_NAME"
    echo "   State: $STATE | CPU: ${CPU}% | Mem: ${MEM}%"
    echo "   Log age: $AGE_MIN min | Status: $STATUS"
    echo "   Output: $OUTPUT_DIR"
    echo ""
done

