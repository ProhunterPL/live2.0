#!/bin/bash
# Check Auto Queue Status
# ========================

PROJECT_ROOT="$HOME/live2.0"
QUEUE_FILE="$PROJECT_ROOT/logs/phase2b_queue.txt"
RESULTS_DIR="$PROJECT_ROOT/results/phase2b_additional/miller_urey_extended"
LOG_FILE="$PROJECT_ROOT/logs/phase2b_auto_restart.log"

echo "=================================================================================="
echo "📊 AUTO QUEUE STATUS"
echo "=================================================================================="
echo ""

# Check currently running
RUNNING=$(ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | wc -l | tr -d ' ')
echo "Currently running: $RUNNING / 5 (max)"
echo ""

if [ "$RUNNING" -gt 0 ]; then
    echo "Active processes:"
    ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | while read line; do
        PID=$(echo "$line" | awk '{print $2}')
        CPU=$(echo "$line" | awk '{print $3}')
        STATE=$(echo "$line" | awk '{print $8}')
        
        # Extract run from command
        OUTPUT_DIR=$(cat "/proc/$PID/cmdline" 2>/dev/null | tr '\0' ' ' | sed -n 's/.*--output[= ]\([^ ]*\).*/\1/p')
        if [ -n "$OUTPUT_DIR" ]; then
            RUN=$(basename "$OUTPUT_DIR")
        else
            RUN="unknown"
        fi
        
        echo "  PID $PID: $RUN (CPU: ${CPU}%, State: $STATE)"
    done
fi

echo ""
echo "=================================================================================="
echo "📋 QUEUE STATUS"
echo "=================================================================================="
echo ""

if [ -f "$QUEUE_FILE" ]; then
    QUEUE_SIZE=$(wc -l < "$QUEUE_FILE" | tr -d ' ')
    echo "Queue size: $QUEUE_SIZE runs waiting"
    echo ""
    
    if [ "$QUEUE_SIZE" -gt 0 ]; then
        echo "Next runs in queue:"
        head -5 "$QUEUE_FILE" | while read item; do
            if [ -n "$item" ]; then
                RUN_ID=$(echo "$item" | cut -d: -f1)
                SEED=$(echo "$item" | cut -d: -f2)
                
                # Check if completed
                RESULTS_FILE="$RESULTS_DIR/run_${RUN_ID}/results.json"
                if [ -f "$RESULTS_FILE" ]; then
                    STATUS="✅ COMPLETED"
                else
                    # Check if running
                    if ps aux | grep "run_phase2_full.py" | grep "run_${RUN_ID}" | grep -v grep > /dev/null; then
                        STATUS="🏃 RUNNING"
                    else
                        STATUS="⏳ WAITING"
                    fi
                fi
                
                echo "  run_${RUN_ID} (seed ${SEED}): $STATUS"
            fi
        done
    else
        echo "✅ Queue is empty - all runs completed or running"
    fi
else
    echo "⚠️  Queue file not found: $QUEUE_FILE"
fi

echo ""
echo "=================================================================================="
echo "📝 RECENT AUTO RESTART LOG"
echo "=================================================================================="
echo ""

if [ -f "$LOG_FILE" ]; then
    echo "Last 20 lines:"
    tail -20 "$LOG_FILE"
else
    echo "⚠️  Log file not found: $LOG_FILE"
fi

echo ""
echo "=================================================================================="
echo "💡 RECOMMENDATIONS"
echo "=================================================================================="
echo ""

AVAILABLE_SLOTS=$((5 - RUNNING))

if [ "$AVAILABLE_SLOTS" -gt 0 ]; then
    echo "✅ $AVAILABLE_SLOTS slot(s) available"
    echo ""
    echo "The auto_queue_restart.sh script should automatically start new runs."
    echo "If it doesn't, check:"
    echo "  1. Is auto_queue_restart.sh running? (ps aux | grep auto_queue_restart)"
    echo "  2. Check log: tail -f $LOG_FILE"
    echo "  3. Manually trigger: bash $PROJECT_ROOT/aws_test/scripts/auto_queue_restart.sh"
else
    echo "⚠️  At maximum capacity (5/5 processes)"
    echo "   Wait for a run to complete, then auto_queue will start next one."
fi

echo ""

