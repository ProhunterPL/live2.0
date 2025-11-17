#!/bin/bash
# Identify which runs correspond to which PIDs and check if they're stuck

echo "=================================================================================="
echo "🔍 IDENTIFY PIDs → RUNS & CHECK STATUS"
echo "=================================================================================="
echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# PIDs from user
PIDS=(126286 146987 147202 171230 171445)

RESULTS_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended"

echo "📊 PID → Run Mapping:"
echo "--------------------------------------------------------------------------------"

for PID in "${PIDS[@]}"; do
    # Get full command line
    if [ ! -f "/proc/$PID/cmdline" ]; then
        echo "❌ PID $PID: Process not found"
        continue
    fi
    
    CMDLINE=$(cat "/proc/$PID/cmdline" | tr '\0' ' ')
    
    # Extract output directory
    OUTPUT_DIR=$(echo "$CMDLINE" | grep -oP '--output\s+\K[^\s]+' || echo "")
    
    if [ -z "$OUTPUT_DIR" ]; then
        echo "⚠️  PID $PID: Could not extract output directory"
        continue
    fi
    
    # Extract run number
    RUN_NUM=$(echo "$OUTPUT_DIR" | grep -oP 'run_\K[0-9]+' || echo "")
    
    if [ -z "$RUN_NUM" ]; then
        echo "⚠️  PID $PID: Could not extract run number from: $OUTPUT_DIR"
        continue
    fi
    
    RUN_DIR="$RESULTS_DIR/run_$RUN_NUM"
    
    # Get process info
    CPU=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
    MEM=$(ps -p $PID -o %mem= 2>/dev/null | tr -d ' ')
    ETIME=$(ps -p $PID -o etime= 2>/dev/null | tr -d ' ')
    STATE=$(ps -p $PID -o stat= 2>/dev/null | tr -d ' ')
    
    echo ""
    echo "🔍 PID $PID → run_$RUN_NUM"
    echo "   CPU: ${CPU}% | Memory: ${MEM}% | Elapsed: $ETIME | State: $STATE"
    
    # Check if completed
    if [ -f "$RUN_DIR/results.json" ]; then
        RESULTS_TIME=$(stat -c %y "$RUN_DIR/results.json" 2>/dev/null | cut -d'.' -f1)
        echo "   ✅ COMPLETED - results.json exists (created: $RESULTS_TIME)"
        echo "   💡 Recommendation: KILL this process (it's completed but still running)"
        continue
    fi
    
    # Check log files
    LOG_FILE="$RUN_DIR/simulation.log"
    RESTART_LOG="$RUN_DIR/simulation_restart.log"
    
    [ ! -f "$LOG_FILE" ] && LOG_FILE="$RESTART_LOG"
    
    if [ -f "$LOG_FILE" ]; then
        LAST_STEP=$(grep -o "Step [0-9,]*" "$LOG_FILE" 2>/dev/null | tail -1 | tr -d ',' | grep -o "[0-9]*")
        LOG_TIME=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
        LOG_AGE_HOURS=$(($(date +%s) - $(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)) / 3600)
        
        if [ -n "$LAST_STEP" ]; then
            PROGRESS=$((LAST_STEP * 100 / 500000))
            echo "   📊 Last step: $LAST_STEP/500,000 ($PROGRESS%)"
            echo "   ⏰ Log time: $LOG_TIME (${LOG_AGE_HOURS}h ago)"
        else
            echo "   ⏰ Log time: $LOG_TIME (${LOG_AGE_HOURS}h ago)"
        fi
        
        # Check if stuck
        if [ "$LOG_AGE_HOURS" -ge 24 ] && (( $(echo "$CPU < 100" | bc -l) )); then
            echo "   ⚠️  WARNING: STUCK - Old log (${LOG_AGE_HOURS}h) + low CPU (${CPU}%)"
            echo "   💡 Recommendation: KILL this process and restart"
        elif [ "$LOG_AGE_HOURS" -ge 48 ]; then
            echo "   ⚠️  WARNING: LIKELY STUCK - Very old log (${LOG_AGE_HOURS}h)"
            echo "   💡 Recommendation: KILL this process and restart"
        elif [ "$LOG_AGE_HOURS" -ge 24 ] && (( $(echo "$CPU > 100" | bc -l) )); then
            echo "   💡 Note: Old log but high CPU - may be log buffering (normal)"
            echo "   💡 Check if near completion (should be close to 500K)"
            if [ -n "$LAST_STEP" ] && [ "$LAST_STEP" -ge 490000 ]; then
                echo "   🎯 Very close to completion - wait a bit longer"
            fi
        else
            echo "   ✅ Appears to be running normally"
        fi
    else
        echo "   ⚠️  No log file found"
    fi
    
    # Show last few log lines
    if [ -f "$LOG_FILE" ]; then
        echo ""
        echo "   📋 Last 3 log entries:"
        tail -3 "$LOG_FILE" 2>/dev/null | sed 's/^/      /'
    fi
done

echo ""
echo "=================================================================================="
echo "💡 RECOMMENDATIONS"
echo "=================================================================================="

# Check for completed runs that are still running
echo "🔍 Checking for completed runs still running:"
for PID in "${PIDS[@]}"; do
    if [ ! -f "/proc/$PID/cmdline" ]; then
        continue
    fi
    
    CMDLINE=$(cat "/proc/$PID/cmdline" | tr '\0' ' ')
    OUTPUT_DIR=$(echo "$CMDLINE" | grep -oP '--output\s+\K[^\s]+' || echo "")
    RUN_NUM=$(echo "$OUTPUT_DIR" | grep -oP 'run_\K[0-9]+' || echo "")
    
    if [ -n "$RUN_NUM" ]; then
        RUN_DIR="$RESULTS_DIR/run_$RUN_NUM"
        if [ -f "$RUN_DIR/results.json" ]; then
            echo "   ⚠️  PID $PID (run_$RUN_NUM): COMPLETED but still running - KILL: kill $PID"
        fi
    fi
done

echo ""
echo "🔍 Checking for stuck runs (old log + low CPU):"
for PID in "${PIDS[@]}"; do
    if [ ! -f "/proc/$PID/cmdline" ]; then
        continue
    fi
    
    CMDLINE=$(cat "/proc/$PID/cmdline" | tr '\0' ' ')
    OUTPUT_DIR=$(echo "$CMDLINE" | grep -oP '--output\s+\K[^\s]+' || echo "")
    RUN_NUM=$(echo "$OUTPUT_DIR" | grep -oP 'run_\K[0-9]+' || echo "")
    
    if [ -n "$RUN_NUM" ]; then
        RUN_DIR="$RESULTS_DIR/run_$RUN_NUM"
        
        # Skip if completed
        [ -f "$RUN_DIR/results.json" ] && continue
        
        LOG_FILE="$RUN_DIR/simulation.log"
        [ ! -f "$LOG_FILE" ] && LOG_FILE="$RUN_DIR/simulation_restart.log"
        
        if [ -f "$LOG_FILE" ]; then
            LOG_AGE_HOURS=$(($(date +%s) - $(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)) / 3600)
            CPU=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
            
            if [ "$LOG_AGE_HOURS" -ge 24 ] && (( $(echo "$CPU < 100" | bc -l) )); then
                LAST_STEP=$(grep -o "Step [0-9,]*" "$LOG_FILE" 2>/dev/null | tail -1 | tr -d ',' | grep -o "[0-9]*")
                if [ -n "$LAST_STEP" ]; then
                    PROGRESS=$((LAST_STEP * 100 / 500000))
                    echo "   ⚠️  PID $PID (run_$RUN_NUM): STUCK - Step $LAST_STEP/500K ($PROGRESS%), log ${LOG_AGE_HOURS}h old, CPU ${CPU}%"
                    echo "      KILL: kill $PID"
                fi
            fi
        fi
    fi
done

echo ""
echo "=================================================================================="

