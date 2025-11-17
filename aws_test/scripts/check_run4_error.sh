#!/bin/bash
# Check why run_4 failed to start

RUN4_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended/run_4"

echo "=================================================================================="
echo "🔍 CHECKING RUN_4 STARTUP ERROR"
echo "=================================================================================="
echo ""

# Check if process is running
PID=$(ps aux | grep "run_4" | grep run_phase2_full | grep -v grep | awk '{print $2}')

if [ -n "$PID" ]; then
    echo "✅ Process is RUNNING (PID: $PID)"
    CPU=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
    echo "   CPU: ${CPU}%"
else
    echo "❌ Process NOT running"
fi

echo ""

# Check restart log
RESTART_LOG="$RUN4_DIR/simulation_restart.log"

if [ -f "$RESTART_LOG" ]; then
    echo "📄 simulation_restart.log:"
    echo "--------------------------------------------------------------------------------"
    echo "Last 50 lines:"
    tail -50 "$RESTART_LOG"
    echo ""
    
    # Check for errors
    echo "🔍 Checking for errors:"
    ERROR_COUNT=$(grep -i "error\|exception\|failed\|traceback\|crash" "$RESTART_LOG" 2>/dev/null | wc -l)
    if [ "$ERROR_COUNT" -gt 0 ]; then
        echo "   ⚠️  Found $ERROR_COUNT error messages:"
        grep -i "error\|exception\|failed\|traceback\|crash" "$RESTART_LOG" 2>/dev/null | tail -10 | sed 's/^/   /'
    else
        echo "   ✅ No obvious errors found"
    fi
    
    # Check if it started successfully
    if grep -q "Starting simulation\|Step 0\|Step 1" "$RESTART_LOG" 2>/dev/null; then
        LAST_STEP=$(grep -o "Step [0-9,]*" "$RESTART_LOG" 2>/dev/null | tail -1 | tr -d ',' | grep -o "[0-9]*")
        if [ -n "$LAST_STEP" ]; then
            echo ""
            echo "   📊 Last step logged: $LAST_STEP"
        fi
    else
        echo ""
        echo "   ⚠️  No simulation start messages found"
    fi
else
    echo "❌ No simulation_restart.log found"
    echo "   This suggests the process never started or crashed immediately"
fi

echo ""
echo "=================================================================================="
echo "💡 RECOMMENDATION"
echo "=================================================================================="

if [ -z "$PID" ]; then
    echo "Run_4 process is not running. Possible causes:"
    echo "1. Immediate crash (check log above for errors)"
    echo "2. Permission issues"
    echo "3. Missing dependencies"
    echo "4. Configuration error"
    echo ""
    echo "Try restarting manually:"
    echo "cd ~/live2.0"
    echo "python3 scripts/run_phase2_full.py \\"
    echo "    --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml \\"
    echo "    --output results/phase2b_additional/miller_urey_extended/run_4 \\"
    echo "    --seed 103 \\"
    echo "    --steps 500000 \\"
    echo "    --force-cpu"
fi

echo "=================================================================================="

