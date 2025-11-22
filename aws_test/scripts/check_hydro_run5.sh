#!/bin/bash
# Quick diagnostic for hydrothermal run_5
# Checks if it's stuck, completed, or still running

RUN_DIR="$HOME/live2.0/results/phase2b_additional/hydrothermal_extended/run_5"
LOG_FILE="$RUN_DIR/simulation.log"
RESULTS_FILE="$RUN_DIR/results.json"

echo "=================================================================================="
echo "🔍 HYDROTHERMAL RUN_5 DIAGNOSTIC"
echo "=================================================================================="
echo ""

# Check if completed
if [ -f "$RESULTS_FILE" ]; then
    echo "✅ Simulation COMPLETED - results.json exists"
    echo ""
    echo "📊 Results summary:"
    python3 -c "import json; d=json.load(open('$RESULTS_FILE')); print(f\"Final step: {d.get('final_step', 'N/A')}\"); print(f\"Total time: {d.get('total_time_seconds', 0)/3600:.2f} hours\")" 2>/dev/null || echo "   (Could not parse results.json)"
    exit 0
fi

# Find process
PID=$(ps aux | grep "run_phase2_full.py" | grep "hydrothermal_extended" | grep "run_5" | grep -v grep | awk '{print $2}')

if [ -z "$PID" ]; then
    echo "❌ Process not found - simulation may have crashed or completed"
    echo ""
    if [ -f "$LOG_FILE" ]; then
        echo "📋 Last log entries:"
        tail -10 "$LOG_FILE" | grep -E "Step|completed|ERROR|Exception" || tail -5 "$LOG_FILE"
    fi
    exit 1
fi

echo "🔄 Process is RUNNING (PID: $PID)"
echo ""

# Get process details
PROC_INFO=$(ps -p $PID -o pid,state,%cpu,%mem,etime,cmd --no-headers 2>/dev/null)
if [ -n "$PROC_INFO" ]; then
    STATE=$(echo "$PROC_INFO" | awk '{print $2}')
    CPU=$(echo "$PROC_INFO" | awk '{print $3}')
    MEM=$(echo "$PROC_INFO" | awk '{print $4}')
    ETIME=$(echo "$PROC_INFO" | awk '{print $5}')
    
    echo "📊 Process Details:"
    echo "   State: $STATE"
    echo "   CPU: ${CPU}%"
    echo "   Memory: ${MEM}%"
    echo "   Running time: $ETIME"
    echo ""
    
    # Check if stuck (sleeping state + high CPU = deadlock)
    if [[ "$STATE" == *"S"* ]] && (( $(echo "$CPU > 100" | bc -l 2>/dev/null || echo 0) )); then
        echo "🚨 WARNING: Process in SLEEPING state but using high CPU (${CPU}%)"
        echo "   This indicates a DEADLOCK (likely cluster detection issue)"
        echo ""
    elif [[ "$STATE" == *"R"* ]]; then
        echo "✅ Process is actively RUNNING (state: R)"
        echo ""
    fi
fi

# Check log file
if [ -f "$LOG_FILE" ]; then
    # Get last step
    LAST_STEP=$(grep -oE "Step [0-9]+ completed" "$LOG_FILE" 2>/dev/null | tail -1 | grep -oE "[0-9]+")
    
    # Get last log modification time
    LAST_LOG_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d' ' -f1,2 | cut -d'.' -f1)
    LOG_AGE_SECONDS=$(( $(date +%s) - $(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0) ))
    LOG_AGE_MINUTES=$((LOG_AGE_SECONDS / 60))
    LOG_AGE_HOURS=$((LOG_AGE_MINUTES / 60))
    
    echo "📝 Log File Status:"
    echo "   Last modified: $LAST_LOG_TIME"
    echo "   Age: ${LOG_AGE_HOURS}h ${LOG_AGE_MINUTES}m"
    echo ""
    
    if [ -n "$LAST_STEP" ]; then
        PROGRESS=$((LAST_STEP * 100 / 500000))
        REMAINING=$((500000 - LAST_STEP))
        
        echo "📊 Progress:"
        echo "   Step: $LAST_STEP / 500,000 ($PROGRESS%)"
        echo "   Remaining: $REMAINING steps"
        echo ""
        
        # Check if should be complete
        if [ "$LAST_STEP" -ge 490000 ]; then
            echo "🎯 Simulation is VERY CLOSE to completion!"
            echo "   Should finish soon (within 1-2 hours)"
            echo ""
            
            # Check if stuck near completion
            if [ "$LOG_AGE_HOURS" -gt 2 ]; then
                echo "⚠️  WARNING: Log is ${LOG_AGE_HOURS}h old but simulation should be done"
                echo "   Process may be stuck in final steps"
                echo ""
            fi
        elif [ "$LAST_STEP" -ge 450000 ]; then
            echo "📈 Simulation is near completion (${LAST_STEP}/500000)"
            echo "   Should complete within 2-3 hours"
            echo ""
        fi
        
        # Estimate time remaining
        if [ "$REMAINING" -gt 0 ]; then
            # Estimate based on typical speed (5-8 steps/sec for CPU)
            EST_SECONDS=$((REMAINING / 6))
            EST_HOURS=$((EST_SECONDS / 3600))
            EST_MINUTES=$((EST_SECONDS % 3600 / 60))
            echo "⏱️  Estimated time remaining: ~${EST_HOURS}h ${EST_MINUTES}m (at 6 steps/sec)"
            echo ""
        fi
    else
        echo "⚠️  Could not find step information in log"
        echo ""
    fi
    
    # Check for errors
    ERRORS=$(grep -iE "error|exception|failed|crash|deadlock" "$LOG_FILE" 2>/dev/null | tail -5)
    if [ -n "$ERRORS" ]; then
        echo "❌ Found errors in log:"
        echo "$ERRORS" | sed 's/^/   /'
        echo ""
    fi
    
    # Show last few log lines
    echo "📋 Last 5 log entries:"
    tail -5 "$LOG_FILE" | sed 's/^/   /'
    echo ""
else
    echo "⚠️  No log file found at: $LOG_FILE"
    echo ""
fi

# Check snapshots
if [ -d "$RUN_DIR/snapshots" ]; then
    SNAPSHOT_COUNT=$(ls -1 "$RUN_DIR/snapshots"/*.json 2>/dev/null | wc -l)
    if [ "$SNAPSHOT_COUNT" -gt 0 ]; then
        LATEST_SNAPSHOT=$(ls -t "$RUN_DIR/snapshots"/*.json 2>/dev/null | head -1)
        SNAPSHOT_TIME=$(stat -c %y "$LATEST_SNAPSHOT" 2>/dev/null | cut -d' ' -f1,2 | cut -d'.' -f1)
        echo "📸 Snapshots: $SNAPSHOT_COUNT files"
        echo "   Latest: $(basename $LATEST_SNAPSHOT) at $SNAPSHOT_TIME"
        echo ""
    fi
fi

# Recommendations
echo "=================================================================================="
echo "💡 RECOMMENDATIONS"
echo "=================================================================================="

if [ -n "$LAST_STEP" ] && [ "$LAST_STEP" -ge 490000 ]; then
    if [ "$LOG_AGE_HOURS" -gt 2 ]; then
        echo "🚨 Simulation should be complete but log is old"
        echo "   Action: Check if process is stuck in final steps"
        echo "   Command: kill -9 $PID (if confirmed stuck)"
    else
        echo "✅ Simulation is near completion - wait a bit longer"
    fi
elif [ -n "$LAST_STEP" ] && [ "$LOG_AGE_HOURS" -gt 4 ]; then
    if [[ "$STATE" == *"S"* ]] && (( $(echo "$CPU > 100" | bc -l 2>/dev/null || echo 0) )); then
        echo "🚨 STUCK - Deadlock detected (sleeping + high CPU)"
        echo "   Action: Kill and restart"
        echo "   Command: kill -9 $PID"
        echo "   Then: bash ~/live2.0/aws_test/scripts/restart_from_checkpoint.sh ~/live2.0/results/phase2b_additional hydrothermal_extended"
    else
        echo "⚠️  Log is old (${LOG_AGE_HOURS}h) but process state unclear"
        echo "   Check CPU usage and process state above"
    fi
elif [ -n "$LAST_STEP" ]; then
    echo "✅ Simulation appears to be running normally"
    echo "   Log buffering is normal - check again in 1-2 hours"
else
    echo "❓ Status unclear - check log file manually"
fi

echo ""

