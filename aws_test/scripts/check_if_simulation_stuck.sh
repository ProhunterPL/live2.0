#!/bin/bash
# Check if simulation is stuck or near completion

RUN_DIR="$1"
if [ -z "$RUN_DIR" ]; then
    echo "Usage: $0 <run_directory>"
    exit 1
fi

echo "=================================================================================="
echo "🔍 CHECKING SIMULATION STATUS: $(basename $RUN_DIR)"
echo "=================================================================================="
echo ""

LOG_FILE="$RUN_DIR/simulation.log"
RESULTS_FILE="$RUN_DIR/results.json"

# Check if completed
if [ -f "$RESULTS_FILE" ]; then
    echo "✅ Simulation COMPLETED - results.json exists"
    exit 0
fi

# Get process info
PID=$(ps aux | grep "$RUN_DIR" | grep run_phase2_full | grep -v grep | awk '{print $2}')
if [ -z "$PID" ]; then
    echo "❌ Process not found - simulation may have crashed"
    exit 1
fi

echo "🔄 Process is RUNNING (PID: $PID)"
echo ""

# Check last step from log
if [ -f "$LOG_FILE" ]; then
    LAST_STEP=$(grep -o "Step [0-9]* completed" "$LOG_FILE" 2>/dev/null | tail -1 | grep -o "[0-9]*")
    LAST_LOG_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d' ' -f1,2 | cut -d'.' -f1)
    
    if [ -n "$LAST_STEP" ]; then
        PROGRESS=$((LAST_STEP * 100 / 500000))
        REMAINING=$((500000 - LAST_STEP))
        
        echo "📊 Last logged step: $LAST_STEP/500,000 ($PROGRESS%)"
        echo "⏰ Last log update: $LAST_LOG_TIME"
        echo "📉 Remaining steps: $REMAINING"
        echo ""
        
        # Estimate time remaining (assuming ~8 steps/sec)
        if [ "$REMAINING" -gt 0 ]; then
            EST_SECONDS=$((REMAINING / 8))
            EST_HOURS=$((EST_SECONDS / 3600))
            EST_MINUTES=$((EST_SECONDS % 3600 / 60))
            echo "⏱️  Estimated time remaining: ~${EST_HOURS}h ${EST_MINUTES}m (at 8 steps/sec)"
        else
            echo "🎯 Simulation should be complete or very close!"
        fi
        echo ""
    fi
fi

# Check process CPU usage
CPU_USAGE=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
if [ -n "$CPU_USAGE" ]; then
    echo "💻 CPU usage: ${CPU_USAGE}%"
    if (( $(echo "$CPU_USAGE > 100" | bc -l) )); then
        echo "   ✅ Process is actively computing (using multiple cores)"
    elif (( $(echo "$CPU_USAGE > 10" | bc -l) )); then
        echo "   ⚠️  Process is using CPU but may be in I/O wait"
    else
        echo "   ⚠️  Process is using very little CPU - may be stuck or waiting"
    fi
    echo ""
fi

# Check checkpoints
if [ -d "$RUN_DIR/checkpoints" ]; then
    LATEST_CHECKPOINT=$(ls -t "$RUN_DIR/checkpoints"/checkpoint_*.json 2>/dev/null | head -1)
    if [ -n "$LATEST_CHECKPOINT" ]; then
        CHECKPOINT_STEP=$(basename "$LATEST_CHECKPOINT" | grep -o "[0-9]*")
        CHECKPOINT_TIME=$(stat -c %y "$LATEST_CHECKPOINT" 2>/dev/null | cut -d' ' -f1,2 | cut -d'.' -f1)
        echo "📁 Latest checkpoint: Step $CHECKPOINT_STEP at $CHECKPOINT_TIME"
        echo ""
    fi
fi

# Check if log mentions completion or errors
if [ -f "$LOG_FILE" ]; then
    echo "📋 Checking for completion messages or errors:"
    echo "--------------------------------------------------------------------------------"
    
    # Check for success messages
    SUCCESS_MSGS=$(grep -i "success\|completed\|finished" "$LOG_FILE" 2>/dev/null | tail -3)
    if [ -n "$SUCCESS_MSGS" ]; then
        echo "✅ Found completion messages:"
        echo "$SUCCESS_MSGS" | sed 's/^/   /'
    else
        echo "   ⚠️  No completion messages found"
    fi
    
    # Check for errors
    ERRORS=$(grep -i "error\|exception\|failed\|crash" "$LOG_FILE" 2>/dev/null | tail -5)
    if [ -n "$ERRORS" ]; then
        echo ""
        echo "❌ Found errors:"
        echo "$ERRORS" | sed 's/^/   /'
    fi
    echo ""
fi

# Recommendations
echo "=================================================================================="
echo "💡 RECOMMENDATIONS"
echo "=================================================================================="

if [ -n "$LAST_STEP" ] && [ "$LAST_STEP" -ge 490000 ]; then
    echo "🎯 Simulation is very close to completion (${LAST_STEP}/500000)"
    echo "   - Check again in 30-60 minutes"
    echo "   - Process may be in final steps or saving results"
elif [ -n "$LAST_STEP" ] && [ "$LAST_STEP" -ge 450000 ]; then
    echo "📊 Simulation is near completion (${LAST_STEP}/500000)"
    echo "   - Should complete within 1-2 hours"
    echo "   - Check periodically for results.json"
else
    echo "⏳ Simulation is still running"
    echo "   - Logs are buffered (old code)"
    echo "   - Process is active (CPU usage confirms)"
    echo "   - Check periodically for results.json"
fi

echo ""

