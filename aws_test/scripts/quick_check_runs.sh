#!/bin/bash
# Quick check for specific runs - check if completed or stuck

echo "=================================================================================="
echo "🔍 QUICK CHECK: Run 12 & Run 4"
echo "=================================================================================="
echo ""

# Check Run 12
RUN12_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended/run_12"
echo "📊 Run 12:"
echo "--------------------------------------------------------------------------------"

if [ -f "$RUN12_DIR/results.json" ]; then
    echo "✅ COMPLETED - results.json exists"
    RESULTS_SIZE=$(ls -lh "$RUN12_DIR/results.json" | awk '{print $5}')
    RESULTS_TIME=$(stat -c %y "$RUN12_DIR/results.json" 2>/dev/null | cut -d'.' -f1)
    echo "   📄 Size: $RESULTS_SIZE"
    echo "   ⏰ Created: $RESULTS_TIME"
else
    echo "⏳ NOT completed (no results.json)"
    
    # Check if process exists
    PID12=$(ps aux | grep "run_12" | grep run_phase2_full | grep -v grep | awk '{print $2}')
    if [ -n "$PID12" ]; then
        CPU12=$(ps -p $PID12 -o %cpu= 2>/dev/null | tr -d ' ')
        echo "   🔄 Process running (PID: $PID12, CPU: ${CPU12}%)"
        
        # Check last log entry
        if [ -f "$RUN12_DIR/simulation.log" ]; then
            LAST_STEP12=$(grep -o "Step [0-9]* completed" "$RUN12_DIR/simulation.log" 2>/dev/null | tail -1 | grep -o "[0-9]*")
            if [ -n "$LAST_STEP12" ]; then
                PROGRESS12=$((LAST_STEP12 * 100 / 500000))
                echo "   📊 Last logged step: $LAST_STEP12/500,000 ($PROGRESS12%)"
            fi
        fi
        
        # Check if process is stuck (low CPU + old log)
        if [ -n "$CPU12" ] && (( $(echo "$CPU12 < 50" | bc -l) )); then
            LOG_AGE12=$(($(date +%s) - $(stat -c %Y "$RUN12_DIR/simulation.log" 2>/dev/null || echo 0)) / 3600)
            if [ "$LOG_AGE12" -ge 1 ]; then
                echo "   ⚠️  WARNING: Low CPU (${CPU12}%) + old log (${LOG_AGE12}h) - may be stuck"
            fi
        fi
    else
        echo "   ❌ No process found - may have crashed or completed"
    fi
fi

echo ""

# Check Run 4
RUN4_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended/run_4"
echo "📊 Run 4:"
echo "--------------------------------------------------------------------------------"

if [ -f "$RUN4_DIR/results.json" ]; then
    echo "✅ COMPLETED - results.json exists"
    RESULTS_SIZE=$(ls -lh "$RUN4_DIR/results.json" | awk '{print $5}')
    RESULTS_TIME=$(stat -c %y "$RUN4_DIR/results.json" 2>/dev/null | cut -d'.' -f1)
    echo "   📄 Size: $RESULTS_SIZE"
    echo "   ⏰ Created: $RESULTS_TIME"
else
    echo "⏳ NOT completed (no results.json)"
    
    # Check if process exists
    PID4=$(ps aux | grep "run_4" | grep run_phase2_full | grep -v grep | awk '{print $2}')
    if [ -n "$PID4" ]; then
        CPU4=$(ps -p $PID4 -o %cpu= 2>/dev/null | tr -d ' ')
        echo "   🔄 Process running (PID: $PID4, CPU: ${CPU4}%)"
        
        # Check last log entry
        if [ -f "$RUN4_DIR/simulation.log" ] || [ -f "$RUN4_DIR/simulation_restart.log" ]; then
            LOG_FILE="$RUN4_DIR/simulation.log"
            [ ! -f "$LOG_FILE" ] && LOG_FILE="$RUN4_DIR/simulation_restart.log"
            
            LAST_STEP4=$(grep -o "Step [0-9]* completed" "$LOG_FILE" 2>/dev/null | tail -1 | grep -o "[0-9]*")
            if [ -n "$LAST_STEP4" ]; then
                PROGRESS4=$((LAST_STEP4 * 100 / 500000))
                LOG_TIME4=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
                echo "   📊 Last logged step: $LAST_STEP4/500,000 ($PROGRESS4%)"
                echo "   ⏰ Last log update: $LOG_TIME4"
            fi
        fi
        
        # Check if process is stuck (low CPU + old log)
        if [ -n "$CPU4" ] && (( $(echo "$CPU4 < 50" | bc -l) )); then
            LOG_FILE="$RUN4_DIR/simulation.log"
            [ ! -f "$LOG_FILE" ] && LOG_FILE="$RUN4_DIR/simulation_restart.log"
            LOG_AGE4=$(($(date +%s) - $(stat -c %Y "$LOG_FILE" 2>/dev/null || echo 0)) / 3600)
            if [ "$LOG_AGE4" -ge 1 ]; then
                echo "   ⚠️  WARNING: Low CPU (${CPU4}%) + old log (${LOG_AGE4}h) - may be stuck"
            fi
        fi
    else
        echo "   ❌ No process found - may have crashed or needs restart"
        
        # Check last known state
        if [ -f "$RUN4_DIR/simulation.log" ] || [ -f "$RUN4_DIR/simulation_restart.log" ]; then
            LOG_FILE="$RUN4_DIR/simulation.log"
            [ ! -f "$LOG_FILE" ] && LOG_FILE="$RUN4_DIR/simulation_restart.log"
            
            LAST_STEP4=$(grep -o "Step [0-9]* completed" "$LOG_FILE" 2>/dev/null | tail -1 | grep -o "[0-9]*")
            LOG_TIME4=$(stat -c %y "$LOG_FILE" 2>/dev/null | cut -d'.' -f1)
            if [ -n "$LAST_STEP4" ]; then
                PROGRESS4=$((LAST_STEP4 * 100 / 500000))
                echo "   📊 Last known step: $LAST_STEP4/500,000 ($PROGRESS4%)"
                echo "   ⏰ Last log update: $LOG_TIME4"
                echo "   💡 Recommendation: Restart this run"
            fi
        fi
    fi
fi

echo ""
echo "=================================================================================="
echo "📊 SUMMARY"
echo "=================================================================================="

# Count active processes
ACTIVE=$(ps aux | grep run_phase2_full | grep -v grep | wc -l)
echo "🔄 Active simulations: $ACTIVE"

# Count completed
COMPLETED=$(find "$HOME/live2.0/results/phase2b_additional/miller_urey_extended" -name "results.json" 2>/dev/null | wc -l)
echo "✅ Completed runs: $COMPLETED"

echo ""

