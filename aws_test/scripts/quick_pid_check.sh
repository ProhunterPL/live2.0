#!/bin/bash
# Quick check: Map PIDs to runs and show status

PIDS=(126286 146987 147202 171230 171445)
RESULTS_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended"

echo "=================================================================================="
echo "🔍 QUICK PID → RUN MAPPING"
echo "=================================================================================="
echo ""

for PID in "${PIDS[@]}"; do
    if [ ! -f "/proc/$PID/cmdline" ]; then
        echo "❌ PID $PID: Process not found"
        continue
    fi
    
    CMDLINE=$(cat "/proc/$PID/cmdline" | tr '\0' ' ')
    
    # Extract run number from output directory
    RUN_NUM=$(echo "$CMDLINE" | sed -n 's/.*run_\([0-9]*\).*/\1/p')
    
    if [ -z "$RUN_NUM" ]; then
        echo "⚠️  PID $PID: Could not identify run"
        echo "   Command: $CMDLINE"
        continue
    fi
    
    RUN_DIR="$RESULTS_DIR/run_$RUN_NUM"
    CPU=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
    
    echo "🔍 PID $PID → run_$RUN_NUM (CPU: ${CPU}%)"
    
    # Check if completed
    if [ -f "$RUN_DIR/results.json" ]; then
        echo "   ✅ COMPLETED - KILL this process!"
        echo "   💡 Command: kill $PID"
        continue
    fi
    
    # Check log
    LOG_FILE="$RUN_DIR/simulation.log"
    [ ! -f "$LOG_FILE" ] && LOG_FILE="$RUN_DIR/simulation_restart.log"
    
    if [ -f "$LOG_FILE" ]; then
        LAST_STEP=$(grep -o "Step [0-9,]*" "$LOG_FILE" 2>/dev/null | tail -1 | tr -d ',' | grep -o "[0-9]*")
        LOG_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
        LOG_MTIME=$(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)
        CURRENT_TIME=$(date +%s)
        LOG_AGE_SECONDS=$((CURRENT_TIME - LOG_MTIME))
        LOG_AGE_HOURS=$((LOG_AGE_SECONDS / 3600))
        
        if [ -n "$LAST_STEP" ]; then
            PROGRESS=$((LAST_STEP * 100 / 500000))
            echo "   📊 Step: $LAST_STEP/500K ($PROGRESS%)"
        fi
        echo "   ⏰ Log: $LOG_TIME (${LOG_AGE_HOURS}h ago)"
        
        # Check if stuck
        if [ "$LOG_AGE_HOURS" -ge 24 ] && (( $(echo "$CPU < 100" | bc -l) )); then
            echo "   ⚠️  STUCK - Old log + low CPU - KILL and restart"
            echo "   💡 Command: kill $PID"
        elif [ "$LOG_AGE_HOURS" -ge 24 ] && (( $(echo "$CPU > 100" | bc -l) )); then
            echo "   💡 Old log but high CPU - may be log buffering"
            if [ -n "$LAST_STEP" ] && [ "$LAST_STEP" -ge 490000 ]; then
                echo "   🎯 Very close to completion - wait a bit"
            else
                echo "   ⚠️  Consider checking if stuck"
            fi
        else
            echo "   ✅ Appears OK"
        fi
    else
        echo "   ⚠️  No log file"
    fi
    echo ""
done

echo "=================================================================================="

