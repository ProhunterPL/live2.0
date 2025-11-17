#!/bin/bash
# Restart stuck runs (runs without processes or with very old logs)

echo "=================================================================================="
echo "🔄 RESTART STUCK RUNS"
echo "=================================================================================="
echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

RESULTS_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended"
CONFIG_FILE="$HOME/live2.0/aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml"

# Seed mapping (from run_phase2b_additional.py)
declare -A SEEDS=(
    [1]=100 [2]=101 [3]=102 [4]=103 [5]=104 [6]=105 [7]=106 [8]=107
    [9]=108 [10]=109 [11]=110 [12]=111 [13]=112 [14]=113 [15]=114
    [16]=115 [17]=116 [18]=117
)

RESTARTED=0
SKIPPED=0

# Check each run
for run_num in {1..18}; do
    RUN_DIR="$RESULTS_DIR/run_$run_num"
    
    # Skip if completed
    if [ -f "$RUN_DIR/results.json" ]; then
        echo "⏭️  run_$run_num: COMPLETED - skipping"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    # Check if process exists
    PID=$(ps aux | grep "run_$run_num" | grep run_phase2_full | grep -v grep | awk '{print $2}')
    
    if [ -n "$PID" ]; then
        # Process exists - check if it's stuck
        CPU=$(ps -p $PID -o %cpu= 2>/dev/null | tr -d ' ')
        LOG_FILE="$RUN_DIR/simulation.log"
        [ ! -f "$LOG_FILE" ] && LOG_FILE="$RUN_DIR/simulation_restart.log"
        
        if [ -f "$LOG_FILE" ]; then
            LOG_AGE_HOURS=$(($(date +%s) - $(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)) / 3600)
            
            # If log is >48h old and CPU is low, consider it stuck
            if [ "$LOG_AGE_HOURS" -ge 48 ] && (( $(echo "$CPU < 50" | bc -l) )); then
                echo "⚠️  run_$run_num: STUCK (PID: $PID, CPU: ${CPU}%, log: ${LOG_AGE_HOURS}h old)"
                echo "   💡 Kill this process first: kill $PID"
                echo "   💡 Then restart manually"
                continue
            else
                echo "⏭️  run_$run_num: RUNNING (PID: $PID, CPU: ${CPU}%) - skipping"
                SKIPPED=$((SKIPPED + 1))
                continue
            fi
        else
            echo "⏭️  run_$run_num: RUNNING (PID: $PID) - skipping"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi
    fi
    
    # No process - check if log exists and is old
    LOG_FILE="$RUN_DIR/simulation.log"
    [ ! -f "$LOG_FILE" ] && LOG_FILE="$RUN_DIR/simulation_restart.log"
    
    if [ -f "$LOG_FILE" ]; then
        LOG_AGE_HOURS=$(($(date +%s) - $(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)) / 3600)
        LAST_STEP=$(grep -o "Step [0-9,]*" "$LOG_FILE" 2>/dev/null | tail -1 | tr -d ',' | grep -o "[0-9]*")
        
        # Only restart if log is >24h old
        if [ "$LOG_AGE_HOURS" -ge 24 ]; then
            SEED=${SEEDS[$run_num]}
            
            if [ -z "$SEED" ]; then
                echo "❌ run_$run_num: Unknown seed - skipping"
                continue
            fi
            
            echo ""
            echo "🚀 Restarting run_$run_num (seed: $SEED, last step: $LAST_STEP, log age: ${LOG_AGE_HOURS}h)"
            
            cd "$HOME/live2.0"
            
            # Start simulation
            setsid nohup python3 scripts/run_phase2_full.py \
                --config "$CONFIG_FILE" \
                --output "$RUN_DIR" \
                --seed "$SEED" \
                --steps 500000 \
                --force-cpu \
                >> "$RUN_DIR/simulation_restart.log" 2>&1 < /dev/null &
            
            NEW_PID=$!
            echo "   ✅ Started with PID: $NEW_PID"
            RESTARTED=$((RESTARTED + 1))
            
            # Wait a moment to verify it started
            sleep 2
            if ! ps -p $NEW_PID > /dev/null 2>&1; then
                echo "   ⚠️  WARNING: Process may have crashed immediately"
            fi
        else
            echo "⏭️  run_$run_num: Log is recent (${LOG_AGE_HOURS}h old) - skipping"
            SKIPPED=$((SKIPPED + 1))
        fi
    else
        echo "⏭️  run_$run_num: No log file - may not have started yet - skipping"
        SKIPPED=$((SKIPPED + 1))
    fi
done

echo ""
echo "=================================================================================="
echo "📊 SUMMARY"
echo "=================================================================================="
echo "🔄 Restarted: $RESTARTED runs"
echo "⏭️  Skipped: $SKIPPED runs"
echo ""

# Count active processes
ACTIVE=$(ps aux | grep run_phase2_full | grep -v grep | wc -l)
echo "🔄 Currently active processes: $ACTIVE"
echo ""

if [ "$RESTARTED" -gt 0 ]; then
    echo "💡 Restarted runs will appear in logs in a few minutes"
    echo "💡 Monitor with: python3 ~/live2.0/aws_test/scripts/check_real_progress.py"
fi

echo "=================================================================================="

