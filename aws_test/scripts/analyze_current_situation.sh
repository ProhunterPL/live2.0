#!/bin/bash
# Analyze current situation - Check which runs are actually stuck
# ================================================================

echo "=================================================================================="
echo "🔍 CURRENT SITUATION ANALYSIS"
echo "=================================================================================="
echo ""

# Process mapping from check_process_details.sh output
declare -A PROCESS_MAP
PROCESS_MAP[73085]="run_6"
PROCESS_MAP[76917]="run_5"
PROCESS_MAP[79713]="run_7"
PROCESS_MAP[79914]="run_8"
PROCESS_MAP[80149]="run_3"

RESULTS_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended"

echo "📊 Process Status:"
echo ""

STUCK_COUNT=0
ACTIVE_COUNT=0

for PID in 73085 76917 79713 79914 80149; do
    RUN=${PROCESS_MAP[$PID]}
    
    # Get process state
    STATE=$(ps -p $PID -o state= 2>/dev/null | tr -d ' ')
    CPU=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
    
    # Check log file
    LOG_FILE="$RESULTS_DIR/$RUN/simulation.log"
    
    if [ -f "$LOG_FILE" ]; then
        # Get last step
        LAST_STEP=$(tail -100 "$LOG_FILE" 2>/dev/null | grep -E "Step [0-9]+" | tail -1 | grep -oE "Step [0-9]+" | head -1)
        
        # Get last modification time
        LAST_MOD=$(stat -c %Y "$LOG_FILE" 2>/dev/null)
        NOW=$(date +%s)
        AGE_HOURS=$(( (NOW - LAST_MOD) / 3600 ))
        AGE_MIN=$(( (NOW - LAST_MOD) / 60 ))
        
        # Determine status
        if [ "$STATE" = "S" ] && [ "${CPU%.*}" -gt 100 ]; then
            if [ "$AGE_MIN" -gt 60 ]; then
                STATUS="🚨 STUCK"
                STUCK_COUNT=$((STUCK_COUNT + 1))
            else
                STATUS="⚠️  SUSPICIOUS"
            fi
        elif [ "$STATE" = "R" ]; then
            STATUS="✅ ACTIVE"
            ACTIVE_COUNT=$((ACTIVE_COUNT + 1))
        else
            STATUS="❓ UNKNOWN"
        fi
        
        echo "PID $PID ($RUN): $STATUS"
        echo "  State: $STATE | CPU: ${CPU}%"
        echo "  Last log: $LAST_STEP"
        if [ "$AGE_HOURS" -ge 1 ]; then
            echo "  Age: ${AGE_HOURS}h ${AGE_MIN}min ago"
        else
            echo "  Age: ${AGE_MIN}min ago"
        fi
        
        # Check if stuck at 160K (deadlock pattern)
        if echo "$LAST_STEP" | grep -qE "160000|160,000"; then
            echo "  ⚠️  STUCK AT 160K - Cluster detection deadlock!"
        fi
    else
        echo "PID $PID ($RUN): ❌ NO LOG FILE"
    fi
    echo ""
done

echo "=================================================================================="
echo "📊 SUMMARY"
echo "=================================================================================="
echo ""
echo "Total processes: 5"
echo "Stuck processes: $STUCK_COUNT"
echo "Active processes: $ACTIVE_COUNT"
echo ""

if [ "$STUCK_COUNT" -gt 0 ]; then
    echo "🚨 ACTION REQUIRED:"
    echo ""
    echo "These processes are deadlocked (Ssl state + high CPU + old logs):"
    for PID in 73085 76917 79713 79914 80149; do
        RUN=${PROCESS_MAP[$PID]}
        STATE=$(ps -p $PID -o state= 2>/dev/null | tr -d ' ')
        CPU=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
        LOG_FILE="$RESULTS_DIR/$RUN/simulation.log"
        
        if [ -f "$LOG_FILE" ]; then
            LAST_MOD=$(stat -c %Y "$LOG_FILE" 2>/dev/null)
            NOW=$(date +%s)
            AGE_MIN=$(( (NOW - LAST_MOD) / 60 ))
            
            if [ "$STATE" = "S" ] && [ "${CPU%.*}" -gt 100 ] && [ "$AGE_MIN" -gt 60 ]; then
                echo "  - PID $PID ($RUN): CPU ${CPU}%, log ${AGE_MIN}min old"
            fi
        fi
    done
    echo ""
    echo "Kill command:"
    echo "  kill -9 73085 76917 79713 79914 80149"
    echo ""
    echo "Or use smart kill script:"
    echo "  bash aws_test/scripts/kill_stuck_smart.sh"
    echo ""
    echo "Then restart maintaining ≤5 processes:"
    echo "  bash aws_test/scripts/auto_queue_restart.sh"
else
    echo "✅ All processes appear to be active!"
fi

echo "=================================================================================="

