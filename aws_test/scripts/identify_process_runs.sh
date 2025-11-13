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
    
    # Extract output directory from --output argument
    # Command format: python3 ... --output /path/to/results/scenario/run_X ...
    OUTPUT_ARG=$(echo "$line" | sed -n 's/.*--output[= ]\([^ ]*\).*/\1/p')
    
    if [ -z "$OUTPUT_ARG" ]; then
        # Try alternative format: --output=/path/to/...
        OUTPUT_ARG=$(echo "$line" | sed -n 's/.*--output=\([^ ]*\).*/\1/p')
    fi
    
    # Extract scenario and run from output path
    # Expected format: .../results/phase2b_additional/scenario/run_X
    if [ -n "$OUTPUT_ARG" ]; then
        SCENARIO=$(echo "$OUTPUT_ARG" | sed -n 's/.*\/\([^/]*_extended\)\/run_[0-9]*.*/\1/p')
        RUN_MATCH=$(echo "$OUTPUT_ARG" | sed -n 's/.*\/run_\([0-9]*\).*/\1/p')
        if [ -n "$RUN_MATCH" ]; then
            RUN_MATCH="run_$RUN_MATCH"
        fi
    fi
    
    # Fallback: try to extract from full command line
    if [ -z "$SCENARIO" ] || [ -z "$RUN_MATCH" ]; then
        RUN_MATCH=$(echo "$line" | grep -oE "run_[0-9]+" | head -1)
        SCENARIO=$(echo "$line" | grep -oE "[a-z_]+_extended" | head -1)
    fi
    
    # Get last log entry for this run
    if [ -n "$RUN_MATCH" ] && [ -n "$SCENARIO" ]; then
        LOG_FILE="$HOME/live2.0/results/phase2b_additional/$SCENARIO/$RUN_MATCH/simulation.log"
        if [ -f "$LOG_FILE" ]; then
            LAST_STEP=$(tail -100 "$LOG_FILE" 2>/dev/null | grep -E "Step [0-9]+" | tail -1 | grep -oE "Step [0-9]+" | head -1)
            LAST_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
            if [ -z "$LAST_STEP" ]; then
                LAST_STEP="(no step info)"
            fi
        else
            LAST_STEP="(no log file)"
            LAST_TIME="(no log file)"
        fi
    else
        # Try to get from output directory directly
        if [ -n "$OUTPUT_ARG" ]; then
            LOG_FILE="$OUTPUT_ARG/simulation.log"
            if [ -f "$LOG_FILE" ]; then
                LAST_STEP=$(tail -100 "$LOG_FILE" 2>/dev/null | grep -E "Step [0-9]+" | tail -1 | grep -oE "Step [0-9]+" | head -1)
                LAST_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
                if [ -z "$LAST_STEP" ]; then
                    LAST_STEP="(no step info)"
                fi
                # Extract run name from path
                RUN_MATCH=$(basename "$OUTPUT_ARG")
                SCENARIO=$(echo "$OUTPUT_ARG" | sed -n 's/.*\/\([^/]*_extended\)\/.*/\1/p')
            else
                LAST_STEP="(no log)"
                LAST_TIME="(no log)"
            fi
        else
            LAST_STEP="(unknown run)"
            LAST_TIME="(unknown)"
        fi
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
    if [ -n "$SCENARIO" ] && [ -n "$RUN_MATCH" ]; then
        echo "   Run: $SCENARIO/$RUN_MATCH"
    else
        echo "   Run: (could not identify)"
        echo "   Command: $(echo "$line" | awk '{for(i=11;i<=NF;i++) printf "%s ", $i; print ""}')"
    fi
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

