#!/bin/bash
# Quick status check for analysis process
# Usage: bash aws_test/scripts/check_analysis_status.sh

cd ~/live2.0

echo "=== Analysis Status Check ==="
echo "Time: $(date)"
echo ""

# Check if process is running
if ps -p 754197 > /dev/null 2>&1; then
    echo "Process: RUNNING (PID 754197)"
    ps -p 754197 -o etime,pcpu,pmem,cmd | tail -1
else
    echo "Process: FINISHED or NOT FOUND"
    echo "Checking for any analyze processes..."
    ps aux | grep '[a]nalyze_phase2b_complete' || echo "No analysis processes running"
fi

echo ""
echo "--- Last Activity (from log) ---"
tail -1 logs/analyze_phase2b_restart.log 2>/dev/null | head -c 100 || echo "Log file not found"

echo ""
echo "--- File Status ---"
if [ -f paper/results_data/hydrothermal_extended_analysis.json ]; then
    SIZE=$(stat -c%s paper/results_data/hydrothermal_extended_analysis.json 2>/dev/null || echo 0)
    echo "File size: $SIZE bytes"
    echo "Last modified: $(stat -c%y paper/results_data/hydrothermal_extended_analysis.json 2>/dev/null | cut -d. -f1)"
else
    echo "Analysis file not found yet"
fi

echo ""
echo "--- Progress Check ---"
if [ -f logs/analyze_phase2b_restart.log ]; then
    PROCESSED=$(grep -c "Processing run_" logs/analyze_phase2b_restart.log 2>/dev/null || echo 0)
    LAST_RUN=$(grep "Processing run_" logs/analyze_phase2b_restart.log 2>/dev/null | tail -1 | grep -o "run_[0-9]*" || echo "unknown")
    echo "Runs processed: $PROCESSED"
    echo "Last run: $LAST_RUN"
    
    # Check if stuck (same run for >30 minutes)
    if [ -f logs/analyze_phase2b_restart.log ]; then
        LAST_TIME=$(grep "Processing run_" logs/analyze_phase2b_restart.log 2>/dev/null | tail -1 | awk '{print $1, $2}')
        if [ -n "$LAST_TIME" ]; then
            echo "Last activity: $LAST_TIME"
        fi
    fi
fi

echo ""
echo "=== End of Status ==="

