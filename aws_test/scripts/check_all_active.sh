#!/bin/bash
# Check all active simulation processes

echo "=================================================================================="
echo "🔍 ALL ACTIVE SIMULATION PROCESSES"
echo "=================================================================================="
echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

RESULTS_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended"

# Get all running processes
ps aux | grep 'run_phase2_full.py' | grep -v grep | while read line; do
    PID=$(echo "$line" | awk '{print $2}')
    CPU=$(echo "$line" | awk '{print $3}')
    MEM=$(echo "$line" | awk '{print $4}')
    
    # Extract run number from command line
    CMDLINE=$(cat "/proc/$PID/cmdline" 2>/dev/null | tr '\0' ' ')
    RUN_NUM=$(echo "$CMDLINE" | sed -n 's/.*run_\([0-9]*\).*/\1/p')
    
    if [ -z "$RUN_NUM" ]; then
        echo "⚠️  PID $PID: Could not identify run"
        continue
    fi
    
    RUN_DIR="$RESULTS_DIR/run_$RUN_NUM"
    
    echo "🔍 PID $PID → run_$RUN_NUM"
    echo "   CPU: ${CPU}% | Memory: ${MEM}%"
    
    # Check if completed
    if [ -f "$RUN_DIR/results.json" ]; then
        echo "   ✅ COMPLETED - Kill this process!"
        continue
    fi
    
    # Check log
    LOG_FILE="$RUN_DIR/simulation.log"
    RESTART_LOG="$RUN_DIR/simulation_restart.log"
    
    [ ! -f "$LOG_FILE" ] && LOG_FILE="$RESTART_LOG"
    
    if [ -f "$LOG_FILE" ]; then
        LAST_STEP=$(grep -o "Step [0-9,]*" "$LOG_FILE" 2>/dev/null | tail -1 | tr -d ',' | grep -o "[0-9]*")
        LOG_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
        LOG_MTIME=$(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)
        CURRENT_TIME=$(date +%s)
        LOG_AGE_SECONDS=$((CURRENT_TIME - LOG_MTIME))
        LOG_AGE_MINUTES=$((LOG_AGE_SECONDS / 60))
        
        if [ -n "$LAST_STEP" ]; then
            PROGRESS=$((LAST_STEP * 100 / 500000))
            echo "   📊 Step: $LAST_STEP/500K ($PROGRESS%)"
        fi
        echo "   ⏰ Log: $LOG_TIME (${LOG_AGE_MINUTES}min ago)"
        
        # Show last log entry
        LAST_LINE=$(tail -1 "$LOG_FILE" 2>/dev/null)
        if [ -n "$LAST_LINE" ]; then
            echo "   📝 Last: $(echo "$LAST_LINE" | cut -c1-80)..."
        fi
    else
        echo "   ⚠️  No log file yet"
    fi
    
    echo ""
done

# Count processes
ACTIVE_COUNT=$(ps aux | grep 'run_phase2_full.py' | grep -v grep | wc -l)
echo "=================================================================================="
echo "📊 SUMMARY: $ACTIVE_COUNT active processes"
echo "=================================================================================="

