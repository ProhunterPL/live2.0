#!/bin/bash
# Check status of all runs - identify stuck vs running vs completed

echo "=================================================================================="
echo "🔍 COMPREHENSIVE STATUS CHECK: All Phase 2B Runs"
echo "=================================================================================="
echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

RESULTS_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended"

# Counters
COMPLETED=0
RUNNING=0
STUCK=0
NOT_STARTED=0

echo "📊 Status by Run:"
echo "--------------------------------------------------------------------------------"

for run_num in {1..18}; do
    RUN_DIR="$RESULTS_DIR/run_$run_num"
    
    # Check if completed
    if [ -f "$RUN_DIR/results.json" ]; then
        echo "✅ run_$run_num: COMPLETED"
        COMPLETED=$((COMPLETED + 1))
        continue
    fi
    
    # Check if process is running
    PID=$(ps aux | grep "run_$run_num" | grep run_phase2_full | grep -v grep | awk '{print $2}')
    
    if [ -n "$PID" ]; then
        CPU=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
        STATE=$(ps -p $PID -o stat= 2>/dev/null | tr -d ' ')
        
        # Find log file
        LOG_FILE="$RUN_DIR/simulation.log"
        [ ! -f "$LOG_FILE" ] && LOG_FILE="$RUN_DIR/simulation_restart.log"
        
        if [ -f "$LOG_FILE" ]; then
            LAST_STEP=$(grep -o "Step [0-9,]*" "$LOG_FILE" 2>/dev/null | tail -1 | tr -d ',' | grep -o "[0-9]*")
            LOG_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
            LOG_AGE_HOURS=$(($(date +%s) - $(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)) / 3600)
            
            if [ -n "$LAST_STEP" ]; then
                PROGRESS=$((LAST_STEP * 100 / 500000))
                
                # Determine if stuck
                if [ "$LOG_AGE_HOURS" -ge 24 ] && (( $(echo "$CPU < 50" | bc -l) )); then
                    echo "⚠️  run_$run_num: STUCK (PID: $PID, CPU: ${CPU}%, Step: $LAST_STEP/500K ($PROGRESS%), Log: ${LOG_AGE_HOURS}h old)"
                    STUCK=$((STUCK + 1))
                elif [ "$LOG_AGE_HOURS" -ge 48 ]; then
                    echo "⚠️  run_$run_num: LIKELY STUCK (PID: $PID, CPU: ${CPU}%, Step: $LAST_STEP/500K ($PROGRESS%), Log: ${LOG_AGE_HOURS}h old)"
                    STUCK=$((STUCK + 1))
                else
                    echo "🔄 run_$run_num: RUNNING (PID: $PID, CPU: ${CPU}%, Step: $LAST_STEP/500K ($PROGRESS%), Log: ${LOG_AGE_HOURS}h old)"
                    RUNNING=$((RUNNING + 1))
                fi
            else
                echo "🔄 run_$run_num: RUNNING (PID: $PID, CPU: ${CPU}%, no step info in log)"
                RUNNING=$((RUNNING + 1))
            fi
        else
            echo "🔄 run_$run_num: RUNNING (PID: $PID, CPU: ${CPU}%, no log file)"
            RUNNING=$((RUNNING + 1))
        fi
    else
        # No process - check if there's a log file
        LOG_FILE="$RUN_DIR/simulation.log"
        [ ! -f "$LOG_FILE" ] && LOG_FILE="$RUN_DIR/simulation_restart.log"
        
        if [ -f "$LOG_FILE" ]; then
            LAST_STEP=$(grep -o "Step [0-9,]*" "$LOG_FILE" 2>/dev/null | tail -1 | tr -d ',' | grep -o "[0-9]*")
            LOG_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
            LOG_AGE_HOURS=$(($(date +%s) - $(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)) / 3600)
            
            if [ -n "$LAST_STEP" ]; then
                PROGRESS=$((LAST_STEP * 100 / 500000))
                echo "❌ run_$run_num: STOPPED (Step: $LAST_STEP/500K ($PROGRESS%), Log: ${LOG_AGE_HOURS}h old) - needs restart"
                STUCK=$((STUCK + 1))
            else
                echo "❓ run_$run_num: UNKNOWN (log exists but no step info)"
                NOT_STARTED=$((NOT_STARTED + 1))
            fi
        else
            echo "❓ run_$run_num: NOT STARTED"
            NOT_STARTED=$((NOT_STARTED + 1))
        fi
    fi
done

echo ""
echo "=================================================================================="
echo "📊 SUMMARY"
echo "=================================================================================="
echo "✅ Completed: $COMPLETED runs"
echo "🔄 Running: $RUNNING processes"
echo "⚠️  Stuck/Stopped: $STUCK runs"
echo "❓ Not started: $NOT_STARTED runs"
echo ""

# List stuck runs that need restart
if [ "$STUCK" -gt 0 ]; then
    echo "⚠️  Runs that need restart:"
    for run_num in {1..18}; do
        RUN_DIR="$RESULTS_DIR/run_$run_num"
        
        # Skip if completed
        [ -f "$RUN_DIR/results.json" ] && continue
        
        # Check if process exists
        PID=$(ps aux | grep "run_$run_num" | grep run_phase2_full | grep -v grep | awk '{print $2}')
        [ -n "$PID" ] && continue
        
        # Check log age
        LOG_FILE="$RUN_DIR/simulation.log"
        [ ! -f "$LOG_FILE" ] && LOG_FILE="$RUN_DIR/simulation_restart.log"
        
        if [ -f "$LOG_FILE" ]; then
            LOG_AGE_HOURS=$(($(date +%s) - $(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)) / 3600)
            if [ "$LOG_AGE_HOURS" -ge 24 ]; then
                LAST_STEP=$(grep -o "Step [0-9,]*" "$LOG_FILE" 2>/dev/null | tail -1 | tr -d ',' | grep -o "[0-9]*")
                if [ -n "$LAST_STEP" ]; then
                    PROGRESS=$((LAST_STEP * 100 / 500000))
                    echo "   - run_$run_num: Step $LAST_STEP/500K ($PROGRESS%), log ${LOG_AGE_HOURS}h old"
                fi
            fi
        fi
    done
    echo ""
fi

# Show active processes
echo "🔄 Currently running processes:"
ps aux | grep run_phase2_full | grep -v grep | while read line; do
    PID=$(echo "$line" | awk '{print $2}')
    CPU=$(echo "$line" | awk '{print $3}')
    RUN_MATCH=$(echo "$line" | grep -oE "run_[0-9]+")
    echo "   - $RUN_MATCH: PID $PID, CPU ${CPU}%"
done

echo ""
echo "=================================================================================="

