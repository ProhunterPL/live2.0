#!/bin/bash
# Kill specific PIDs after verification

PIDS=(126286 146987 147202 171230 171445)

echo "=================================================================================="
echo "🔍 VERIFY PIDs BEFORE KILLING"
echo "=================================================================================="

RESULTS_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended"

TO_KILL=()
TO_KEEP=()

for PID in "${PIDS[@]}"; do
    if [ ! -f "/proc/$PID/cmdline" ]; then
        echo "❌ PID $PID: Process not found (already dead?)"
        continue
    fi
    
    CMDLINE=$(cat "/proc/$PID/cmdline" | tr '\0' ' ')
    OUTPUT_DIR=$(echo "$CMDLINE" | grep -oP '--output\s+\K[^\s]+' || echo "")
    RUN_NUM=$(echo "$OUTPUT_DIR" | grep -oP 'run_\K[0-9]+' || echo "")
    
    if [ -z "$RUN_NUM" ]; then
        echo "⚠️  PID $PID: Could not identify run - SKIPPING"
        TO_KEEP+=($PID)
        continue
    fi
    
    RUN_DIR="$RESULTS_DIR/run_$RUN_NUM"
    CPU=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
    
    echo ""
    echo "🔍 PID $PID → run_$RUN_NUM (CPU: ${CPU}%)"
    
    # Check if completed
    if [ -f "$RUN_DIR/results.json" ]; then
        echo "   ✅ COMPLETED - Safe to kill"
        TO_KILL+=($PID)
        continue
    fi
    
    # Check log
    LOG_FILE="$RUN_DIR/simulation.log"
    [ ! -f "$LOG_FILE" ] && LOG_FILE="$RUN_DIR/simulation_restart.log"
    
    if [ -f "$LOG_FILE" ]; then
        LOG_AGE_HOURS=$(($(date +%s) - $(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)) / 3600)
        LAST_STEP=$(grep -o "Step [0-9,]*" "$LOG_FILE" 2>/dev/null | tail -1 | tr -d ',' | grep -o "[0-9]*")
        
        # Stuck criteria: old log (>24h) + low CPU (<100%)
        if [ "$LOG_AGE_HOURS" -ge 24 ] && (( $(echo "$CPU < 100" | bc -l) )); then
            if [ -n "$LAST_STEP" ]; then
                PROGRESS=$((LAST_STEP * 100 / 500000))
                echo "   ⚠️  STUCK - Step $LAST_STEP/500K ($PROGRESS%), log ${LOG_AGE_HOURS}h old, CPU ${CPU}%"
                echo "   💡 Safe to kill and restart"
                TO_KILL+=($PID)
            else
                echo "   ⚠️  STUCK - Log ${LOG_AGE_HOURS}h old, CPU ${CPU}%"
                echo "   💡 Safe to kill and restart"
                TO_KILL+=($PID)
            fi
        elif [ "$LOG_AGE_HOURS" -ge 48 ]; then
            echo "   ⚠️  VERY OLD LOG (${LOG_AGE_HOURS}h) - Likely stuck"
            echo "   💡 Safe to kill and restart"
            TO_KILL+=($PID)
        else
            echo "   ✅ Appears to be running normally - KEEP"
            TO_KEEP+=($PID)
        fi
    else
        echo "   ⚠️  No log file - Unknown state - SKIPPING"
        TO_KEEP+=($PID)
    fi
done

echo ""
echo "=================================================================================="
echo "📊 SUMMARY"
echo "=================================================================================="
echo "🔴 To KILL (${#TO_KILL[@]} processes):"
for PID in "${TO_KILL[@]}"; do
    echo "   - PID $PID"
done

echo ""
echo "🟢 To KEEP (${#TO_KEEP[@]} processes):"
for PID in "${TO_KEEP[@]}"; do
    echo "   - PID $PID"
done

echo ""
echo "=================================================================================="

if [ ${#TO_KILL[@]} -eq 0 ]; then
    echo "✅ No processes to kill - all appear to be running normally"
    exit 0
fi

echo "⚠️  About to kill ${#TO_KILL[@]} process(es)"
echo ""
read -p "Continue? (yes/no): " CONFIRM

if [ "$CONFIRM" != "yes" ]; then
    echo "❌ Cancelled"
    exit 1
fi

echo ""
echo "🔴 Killing processes..."

for PID in "${TO_KILL[@]}"; do
    echo "   Killing PID $PID..."
    kill $PID 2>/dev/null
    
    # Wait a moment
    sleep 1
    
    # Check if still alive
    if ps -p $PID > /dev/null 2>&1; then
        echo "   ⚠️  Process still alive - sending SIGKILL..."
        kill -9 $PID 2>/dev/null
        sleep 1
    fi
    
    if ! ps -p $PID > /dev/null 2>&1; then
        echo "   ✅ PID $PID killed"
    else
        echo "   ❌ Failed to kill PID $PID"
    fi
done

echo ""
echo "=================================================================================="
echo "✅ Done"
echo "=================================================================================="
echo ""
echo "💡 Next steps:"
echo "   1. Restart killed runs: bash ~/live2.0/aws_test/scripts/restart_stuck_runs.sh"
echo "   2. Monitor progress: python3 ~/live2.0/aws_test/scripts/check_real_progress.py"
echo ""

