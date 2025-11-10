#!/bin/bash
# Check if a simulation run completed successfully

RUN_DIR="$1"
if [ -z "$RUN_DIR" ]; then
    echo "Usage: $0 <run_directory>"
    exit 1
fi

echo "=================================================================================="
echo "🔍 CHECKING COMPLETION STATUS: $(basename $RUN_DIR)"
echo "=================================================================================="
echo ""

RESULTS_FILE="$RUN_DIR/results.json"
SUMMARY_FILE="$RUN_DIR/summary.txt"
LOG_FILE="$RUN_DIR/simulation.log"

# Check for results files
if [ -f "$RESULTS_FILE" ]; then
    echo "✅ results.json EXISTS - Simulation completed successfully!"
    echo "   File: $RESULTS_FILE"
    echo "   Size: $(ls -lh "$RESULTS_FILE" | awk '{print $5}')"
    echo "   Modified: $(stat -c %y "$RESULTS_FILE" | cut -d' ' -f1,2 | cut -d'.' -f1)"
    echo ""
    
    # Check summary
    if [ -f "$SUMMARY_FILE" ]; then
        echo "✅ summary.txt EXISTS"
        echo "   File: $SUMMARY_FILE"
        echo "   Size: $(ls -lh "$SUMMARY_FILE" | awk '{print $5}')"
        echo ""
    fi
    
    exit 0
fi

if [ -f "$SUMMARY_FILE" ]; then
    echo "✅ summary.txt EXISTS - Simulation completed!"
    echo "   File: $SUMMARY_FILE"
    echo "   Size: $(ls -lh "$SUMMARY_FILE" | awk '{print $5}')"
    echo "   Modified: $(stat -c %y "$SUMMARY_FILE" | cut -d' ' -f1,2 | cut -d'.' -f1)"
    echo ""
    echo "⚠️  Note: results.json not found, but summary.txt exists"
    echo "   This may indicate partial completion or different output format"
    echo ""
    exit 0
fi

# Check log for completion messages
if [ -f "$LOG_FILE" ]; then
    echo "📋 Checking log file for completion status..."
    echo "--------------------------------------------------------------------------------"
    
    # Get last step
    LAST_STEP=$(grep -o "Step [0-9]* completed" "$LOG_FILE" 2>/dev/null | tail -1 | grep -o "[0-9]*")
    LAST_LOG_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d' ' -f1,2 | cut -d'.' -f1)
    
    if [ -n "$LAST_STEP" ]; then
        PROGRESS=$((LAST_STEP * 100 / 500000))
        echo "📊 Last logged step: $LAST_STEP/500,000 ($PROGRESS%)"
        echo "⏰ Last log update: $LAST_LOG_TIME"
        echo ""
    fi
    
    # Check for completion messages
    echo "🔍 Searching for completion messages:"
    COMPLETION_MSGS=$(grep -i "success\|completed.*simulation\|finished\|phase 2.*complete" "$LOG_FILE" 2>/dev/null | tail -5)
    if [ -n "$COMPLETION_MSGS" ]; then
        echo "✅ Found completion messages:"
        echo "$COMPLETION_MSGS" | sed 's/^/   /'
    else
        echo "   ⚠️  No explicit completion messages found"
    fi
    echo ""
    
    # Check for errors
    ERRORS=$(grep -i "error\|exception\|failed\|crash\|traceback" "$LOG_FILE" 2>/dev/null | tail -10)
    if [ -n "$ERRORS" ]; then
        echo "❌ Found errors in log:"
        echo "$ERRORS" | sed 's/^/   /'
        echo ""
    fi
    
    # Check last 10 lines
    echo "📋 Last 10 lines of log:"
    echo "--------------------------------------------------------------------------------"
    tail -10 "$LOG_FILE" 2>/dev/null | sed 's/^/   /'
    echo ""
fi

# Check if process is still running
PID=$(ps aux | grep "$RUN_DIR" | grep run_phase2_full | grep -v grep | awk '{print $2}')
if [ -n "$PID" ]; then
    echo "🔄 Process is still RUNNING (PID: $PID)"
    echo "   Simulation may still be in progress"
else
    echo "⏸️  Process is NOT running"
    echo "   Simulation has stopped"
    
    if [ -z "$LAST_STEP" ] || [ "$LAST_STEP" -lt 500000 ]; then
        echo "   ⚠️  Simulation stopped before reaching 500,000 steps"
        echo "   Possible causes:"
        echo "      - Timeout (old 6h limit)"
        echo "      - Error/crash"
        echo "      - Manual termination"
    fi
fi

echo ""
echo "=================================================================================="
echo "💡 RECOMMENDATIONS"
echo "=================================================================================="

if [ -f "$RESULTS_FILE" ] || [ -f "$SUMMARY_FILE" ]; then
    echo "✅ Simulation completed successfully!"
    echo "   - Results are available"
    echo "   - You can proceed with analysis"
elif [ -n "$LAST_STEP" ] && [ "$LAST_STEP" -ge 490000 ]; then
    echo "🎯 Simulation was very close to completion (${LAST_STEP}/500000)"
    echo "   - May have completed but failed to save results.json"
    echo "   - Check for any error messages in log"
    echo "   - Consider re-running from checkpoint if available"
elif [ -n "$LAST_STEP" ] && [ "$LAST_STEP" -ge 400000 ]; then
    echo "📊 Simulation made significant progress (${LAST_STEP}/500000)"
    echo "   - May have been stopped by timeout"
    echo "   - Check for checkpoints to resume"
    echo "   - Consider re-running with longer timeout"
else
    echo "⚠️  Simulation stopped early (${LAST_STEP:-0}/500000)"
    echo "   - Check logs for errors"
    echo "   - May need to restart"
fi

echo ""

