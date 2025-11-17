#!/bin/bash
# Detailed check for Run 12 and Run 4

echo "=================================================================================="
echo "🔍 DETAILED CHECK: Run 12 & Run 4"
echo "=================================================================================="
echo ""

# Check Run 12
RUN12_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended/run_12"
echo "📊 Run 12 (PID: 126286):"
echo "--------------------------------------------------------------------------------"

# Check process details
PID12=$(ps aux | grep "run_12" | grep run_phase2_full | grep -v grep | awk '{print $2}')
if [ -n "$PID12" ]; then
    CPU12=$(ps -p $PID12 -o %cpu= 2>/dev/null | tr -d ' ')
    MEM12=$(ps -p $PID12 -o %mem= 2>/dev/null | tr -d ' ')
    ETIME12=$(ps -p $PID12 -o etime= 2>/dev/null | tr -d ' ')
    STATE12=$(ps -p $PID12 -o stat= 2>/dev/null | tr -d ' ')
    
    echo "✅ Process is RUNNING"
    echo "   PID: $PID12"
    echo "   CPU: ${CPU12}%"
    echo "   Memory: ${MEM12}%"
    echo "   Elapsed time: $ETIME12"
    echo "   State: $STATE12"
    
    # Check if state indicates stuck
    if [[ "$STATE12" == *"D"* ]]; then
        echo "   ⚠️  WARNING: Process in uninterruptible sleep (I/O wait)"
    elif [[ "$STATE12" == *"R"* ]]; then
        echo "   ✅ Process is actively running"
    fi
    
    # Check log files
    LOG_FILE="$RUN12_DIR/simulation.log"
    RESTART_LOG="$RUN12_DIR/simulation_restart.log"
    
    if [ -f "$LOG_FILE" ]; then
        echo ""
        echo "📄 simulation.log:"
        LAST_STEP12=$(grep -o "Step [0-9]* completed" "$LOG_FILE" 2>/dev/null | tail -1 | grep -o "[0-9]*")
        LOG_SIZE12=$(ls -lh "$LOG_FILE" | awk '{print $5}')
        LOG_TIME12=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
        LOG_MTIME=$(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)
        CURRENT_TIME=$(date +%s)
        LOG_AGE_SECONDS=$((CURRENT_TIME - LOG_MTIME))
        LOG_AGE_HOURS12=$((LOG_AGE_SECONDS / 3600))
        
        if [ -n "$LAST_STEP12" ]; then
            PROGRESS12=$((LAST_STEP12 * 100 / 500000))
            echo "   Last step: $LAST_STEP12/500,000 ($PROGRESS12%)"
        fi
        echo "   Size: $LOG_SIZE12"
        echo "   Last modified: $LOG_TIME12 (${LOG_AGE_HOURS12}h ago)"
        
        # Check last 5 log lines
        echo ""
        echo "   Last 5 log entries:"
        tail -5 "$LOG_FILE" 2>/dev/null | sed 's/^/   /'
        
        # Check for errors
        echo ""
        echo "   Checking for errors:"
        ERROR_COUNT=$(grep -i "error\|exception\|failed\|crash" "$LOG_FILE" 2>/dev/null | wc -l)
        if [ "$ERROR_COUNT" -gt 0 ]; then
            echo "   ⚠️  Found $ERROR_COUNT error/exception messages"
            grep -i "error\|exception\|failed\|crash" "$LOG_FILE" 2>/dev/null | tail -3 | sed 's/^/      /'
        else
            echo "   ✅ No errors found"
        fi
        
        # Check for completion messages
        COMPLETION=$(grep -i "success\|completed\|finished\|simulation complete" "$LOG_FILE" 2>/dev/null | tail -3)
        if [ -n "$COMPLETION" ]; then
            echo ""
            echo "   🎯 Found completion messages:"
            echo "$COMPLETION" | sed 's/^/      /'
        fi
    fi
    
    if [ -f "$RESTART_LOG" ]; then
        echo ""
        echo "📄 simulation_restart.log:"
        RESTART_SIZE=$(ls -lh "$RESTART_LOG" | awk '{print $5}')
        RESTART_TIME=$(stat -c %y "$RESTART_LOG" 2>/dev/null | cut -d'.' -f1)
        echo "   Size: $RESTART_SIZE"
        echo "   Last modified: $RESTART_TIME"
        
        LAST_STEP_RESTART=$(grep -o "Step [0-9]* completed" "$RESTART_LOG" 2>/dev/null | tail -1 | grep -o "[0-9]*")
        if [ -n "$LAST_STEP_RESTART" ]; then
            PROGRESS_RESTART=$((LAST_STEP_RESTART * 100 / 500000))
            echo "   Last step: $LAST_STEP_RESTART/500,000 ($PROGRESS_RESTART%)"
        fi
        
        echo ""
        echo "   Last 5 log entries:"
        tail -5 "$RESTART_LOG" 2>/dev/null | sed 's/^/   /'
    fi
    
    # Check snapshots/checkpoints
    echo ""
    echo "📸 Snapshots/Checkpoints:"
    if [ -d "$RUN12_DIR/snapshots" ]; then
        SNAPSHOT_COUNT=$(ls "$RUN12_DIR/snapshots"/*.json 2>/dev/null | wc -l)
        LATEST_SNAPSHOT=$(ls -t "$RUN12_DIR/snapshots"/*.json 2>/dev/null | head -1)
        if [ -n "$LATEST_SNAPSHOT" ]; then
            SNAPSHOT_TIME=$(stat -c %y "$LATEST_SNAPSHOT" 2>/dev/null | cut -d'.' -f1)
            echo "   Snapshots: $SNAPSHOT_COUNT files"
            echo "   Latest: $(basename $LATEST_SNAPSHOT) @ $SNAPSHOT_TIME"
        else
            echo "   Snapshots: 0 files"
        fi
    else
        echo "   Snapshots: directory not found"
    fi
    
    if [ -d "$RUN12_DIR/checkpoints" ]; then
        CHECKPOINT_COUNT=$(ls "$RUN12_DIR/checkpoints"/*.json 2>/dev/null | wc -l)
        LATEST_CHECKPOINT=$(ls -t "$RUN12_DIR/checkpoints"/*.json 2>/dev/null | head -1)
        if [ -n "$LATEST_CHECKPOINT" ]; then
            CHECKPOINT_TIME=$(stat -c %y "$LATEST_CHECKPOINT" 2>/dev/null | cut -d'.' -f1)
            echo "   Checkpoints: $CHECKPOINT_COUNT files"
            echo "   Latest: $(basename $LATEST_CHECKPOINT) @ $CHECKPOINT_TIME"
        else
            echo "   Checkpoints: 0 files"
        fi
    else
        echo "   Checkpoints: directory not found"
    fi
    
    # Analysis
    echo ""
    echo "🔍 Analysis:"
    if [ -n "$LAST_STEP12" ] && [ -n "$LOG_AGE_HOURS12" ] && [ "$LAST_STEP12" -ge 490000 ]; then
        echo "   🎯 Run 12 is VERY CLOSE to completion (${LAST_STEP12}/500,000)"
        echo "   💡 Process is still running - may be in final steps or extracting results"
        echo "   💡 High CPU (${CPU12}%) suggests active work"
        echo "   💡 Old log (${LOG_AGE_HOURS12}h) suggests log buffering"
        echo "   💡 Recommendation: Wait 1-2 hours, then check for results.json"
    elif [ -n "$LAST_STEP12" ] && [ "$LAST_STEP12" -ge 450000 ]; then
        echo "   📊 Run 12 is near completion (${LAST_STEP12}/500,000)"
        echo "   💡 Should complete within 2-3 hours"
    elif [ "$LOG_AGE_HOURS12" -ge 24 ] && (( $(echo "$CPU12 < 50" | bc -l) )); then
        echo "   ⚠️  WARNING: Run 12 may be STUCK"
        echo "   💡 Old log (${LOG_AGE_HOURS12}h) + low CPU (${CPU12}%)"
        echo "   💡 Recommendation: Check process state, consider killing if stuck"
    else
        echo "   ✅ Run 12 appears to be progressing normally"
        echo "   💡 High CPU (${CPU12}%) indicates active computation"
    fi
else
    echo "❌ Process not found"
fi

echo ""
echo "=================================================================================="

# Check Run 4
RUN4_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended/run_4"
echo "📊 Run 4:"
echo "--------------------------------------------------------------------------------"

PID4=$(ps aux | grep "run_4" | grep run_phase2_full | grep -v grep | awk '{print $2}')
if [ -n "$PID4" ]; then
    echo "✅ Process is RUNNING (PID: $PID4)"
    CPU4=$(ps -p $PID4 -o %cpu= 2>/dev/null | tr -d ' ')
    echo "   CPU: ${CPU4}%"
else
    echo "❌ Process NOT running"
    
    # Check last known state
    LOG_FILE="$RUN4_DIR/simulation.log"
    RESTART_LOG="$RUN4_DIR/simulation_restart.log"
    
    if [ -f "$LOG_FILE" ] || [ -f "$RESTART_LOG" ]; then
        [ ! -f "$LOG_FILE" ] && LOG_FILE="$RESTART_LOG"
        
        LAST_STEP4=$(grep -o "Step [0-9]* completed" "$LOG_FILE" 2>/dev/null | tail -1 | grep -o "[0-9]*")
        LOG_TIME4=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
        LOG_MTIME=$(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)
        CURRENT_TIME=$(date +%s)
        LOG_AGE_SECONDS=$((CURRENT_TIME - LOG_MTIME))
        LOG_AGE_HOURS4=$((LOG_AGE_SECONDS / 3600))
        
        if [ -n "$LAST_STEP4" ]; then
            PROGRESS4=$((LAST_STEP4 * 100 / 500000))
            echo ""
            echo "📄 Last known state:"
            echo "   Last step: $LAST_STEP4/500,000 ($PROGRESS4%)"
            echo "   Last log update: $LOG_TIME4 (${LOG_AGE_HOURS4}h ago)"
            
            echo ""
            echo "   Last 5 log entries:"
            tail -5 "$LOG_FILE" 2>/dev/null | sed 's/^/   /'
            
            echo ""
            echo "💡 Recommendation:"
            if [ -f "$RUN4_DIR/results.json" ]; then
                echo "   ✅ Run 4 is COMPLETED (results.json exists)"
            else
                echo "   🔄 Restart Run 4 from beginning"
                echo "   📝 Command:"
                echo "      cd ~/live2.0"
                echo "      setsid nohup python3 scripts/run_phase2_full.py \\"
                echo "          --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml \\"
                echo "          --output results/phase2b_additional/miller_urey_extended/run_4 \\"
                echo "          --seed 103 \\"
                echo "          --steps 500000 \\"
                echo "          --force-cpu \\"
                echo "          >> results/phase2b_additional/miller_urey_extended/run_4/simulation_restart.log 2>&1 < /dev/null &"
            fi
        fi
    else
        echo "   ⚠️  No log file found - run may not have started"
    fi
fi

echo ""
echo "=================================================================================="

