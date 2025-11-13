#!/bin/bash
# Check detailed process information including full command line
# ============================================================

echo "=================================================================================="
echo "🔍 DETAILED PROCESS INFORMATION"
echo "=================================================================================="
echo ""

ps aux | grep 'run_phase2_full.py' | grep -v grep | while read line; do
    PID=$(echo "$line" | awk '{print $2}')
    STATE=$(echo "$line" | awk '{print $8}')
    CPU=$(echo "$line" | awk '{print $3}')
    MEM=$(echo "$line" | awk '{print $4}')
    TIME=$(echo "$line" | awk '{print $10}')
    
    echo "PID: $PID"
    echo "  State: $STATE | CPU: ${CPU}% | Memory: ${MEM}% | Time: $TIME"
    echo "  Full command:"
    
    # Get full command line from /proc/PID/cmdline
    if [ -f "/proc/$PID/cmdline" ]; then
        CMDLINE=$(cat "/proc/$PID/cmdline" | tr '\0' ' ')
        echo "    $CMDLINE"
        
        # Extract key information using sed instead of grep
        OUTPUT_DIR=$(echo "$CMDLINE" | sed -n 's/.*--output[= ]\([^ ]*\).*/\1/p')
        CONFIG=$(echo "$CMDLINE" | sed -n 's/.*--config[= ]\([^ ]*\).*/\1/p')
        SEED=$(echo "$CMDLINE" | sed -n 's/.*--seed[= ]\([^ ]*\).*/\1/p')
        
        echo "  Extracted:"
        [ -n "$OUTPUT_DIR" ] && echo "    Output: $OUTPUT_DIR"
        [ -n "$CONFIG" ] && echo "    Config: $CONFIG"
        [ -n "$SEED" ] && echo "    Seed: $SEED"
        
        # Extract run number and scenario from output directory
        if [ -n "$OUTPUT_DIR" ]; then
            RUN_NUM=$(echo "$OUTPUT_DIR" | grep -oE "run_[0-9]+" | head -1)
            SCENARIO=$(echo "$OUTPUT_DIR" | grep -oE "[a-z_]+_extended" | head -1)
            
            if [ -n "$RUN_NUM" ] && [ -n "$SCENARIO" ]; then
                echo "    Run: $SCENARIO/$RUN_NUM"
                
                # Check log file
                LOG_FILE="$OUTPUT_DIR/simulation.log"
                if [ -f "$LOG_FILE" ]; then
                    LAST_STEP=$(tail -100 "$LOG_FILE" 2>/dev/null | grep -E "Step [0-9]+" | tail -1 | grep -oE "Step [0-9]+" | head -1)
                    LAST_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
                    AGE_HOURS=$(($(date +%s) - $(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)) / 3600)
                    
                    echo "    Log: $LAST_STEP @ $LAST_TIME (${AGE_HOURS}h ago)"
                    
                    if [ "$AGE_HOURS" -ge 24 ]; then
                        echo "    ⚠️  STUCK: Log is ${AGE_HOURS}h old"
                    fi
                else
                    echo "    ⚠️  No log file found"
                fi
            fi
        fi
    else
        echo "    (Could not read /proc/$PID/cmdline)"
    fi
    
    echo ""
done

echo "=================================================================================="

