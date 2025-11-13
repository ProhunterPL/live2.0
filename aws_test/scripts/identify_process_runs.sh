#!/bin/bash
# Identify which run each process is handling
# ===========================================

echo "=================================================================================="
echo "🔍 PROCESS TO RUN MAPPING"
echo "=================================================================================="
echo ""

ps aux | grep 'run_phase2_full.py' | grep -v grep | while read line; do
    PID=$(echo "$line" | awk '{print $2}')
    STATE=$(echo "$line" | awk '{print $8}')
    CPU=$(echo "$line" | awk '{print $3}')
    
    # Extract run directory from command
    RUN_MATCH=$(echo "$line" | grep -o "run_[0-9]*" | head -1)
    SCENARIO=$(echo "$line" | grep -o "[a-z_]*_extended" | head -1)
    
    # Get last log entry for this run
    if [ -n "$RUN_MATCH" ] && [ -n "$SCENARIO" ]; then
        LOG_FILE="$HOME/live2.0/results/phase2b_additional/$SCENARIO/$RUN_MATCH/simulation.log"
        if [ -f "$LOG_FILE" ]; then
            LAST_STEP=$(tail -100 "$LOG_FILE" | grep -E "Step [0-9]+" | tail -1 | grep -oE "Step [0-9]+" | head -1)
            LAST_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
            if [ -z "$LAST_STEP" ]; then
                LAST_STEP="(no step info)"
            fi
        else
            LAST_STEP="(no log)"
            LAST_TIME="(no log)"
        fi
    else
        LAST_STEP="(unknown run)"
        LAST_TIME="(unknown)"
    fi
    
    STATE_DESC=""
    case "$STATE" in
        R*) STATE_DESC="✅ Running" ;;
        S*) STATE_DESC="😴 Sleeping" ;;
        D*) STATE_DESC="⚠️  I/O Wait" ;;
        Z*) STATE_DESC="💀 Zombie" ;;
        *) STATE_DESC="❓ Unknown" ;;
    esac
    
    echo "PID $PID: $STATE_DESC (CPU: ${CPU}%)"
    echo "   Run: $SCENARIO/$RUN_MATCH"
    echo "   Last: $LAST_STEP @ $LAST_TIME"
    echo ""
done

echo "=================================================================================="
echo "💡 RECOMMENDATION"
echo "=================================================================================="
echo ""
echo "Processes in 'Ssl' state with high CPU (>100%) are likely DEADLOCKED."
echo "Kill them and restart with fixed code."
echo ""
echo "To kill stuck processes:"
echo "  bash aws_test/scripts/kill_stuck_phase2b.sh"
echo ""
echo "Then restart with:"
echo "  bash aws_test/scripts/auto_queue_restart.sh"
echo ""

