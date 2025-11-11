#!/bin/bash
# Automatically restart stuck simulations when slots become available
# =================================================================

RESULTS_DIR="${1:-$HOME/live2.0/results/phase2b_additional}"
MAX_PARALLEL=${2:-4}
CHECK_INTERVAL=${3:-300}  # Check every 5 minutes

echo "🔄 Auto-restart monitor"
echo "   Results dir: $RESULTS_DIR"
echo "   Max parallel: $MAX_PARALLEL"
echo "   Check interval: $CHECK_INTERVAL seconds"
echo ""
echo "This will monitor and restart stuck simulations when slots become available."
echo "Press Ctrl+C to stop."
echo ""

cd ~/live2.0 || exit 1

# Make scripts executable
chmod +x aws_test/scripts/*.sh

ITERATION=0

while true; do
    ITERATION=$((ITERATION + 1))
    echo ""
    echo "================================================================================"
    echo "Check #$ITERATION - $(date '+%Y-%m-%d %H:%M:%S')"
    echo "================================================================================"
    
    # Count currently running
    CURRENT_RUNNING=$(ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l | tr -d ' ')
    AVAILABLE_SLOTS=$((MAX_PARALLEL - CURRENT_RUNNING))
    
    echo "📊 Currently running: $CURRENT_RUNNING / $MAX_PARALLEL"
    echo "   Available slots: $AVAILABLE_SLOTS"
    
    if [ "$AVAILABLE_SLOTS" -le 0 ]; then
        echo "⏳ No slots available, waiting..."
        sleep $CHECK_INTERVAL
        continue
    fi
    
    echo ""
    echo "✅ Slots available! Checking for stuck simulations..."
    
    # Check hydrothermal runs
    HYDRO_RESTARTED=0
    for i in {1..10}; do
        if [ "$AVAILABLE_SLOTS" -le 0 ]; then
            break
        fi
        
        RUN_DIR="$RESULTS_DIR/hydrothermal_extended/run_$i"
        if [ ! -d "$RUN_DIR" ]; then
            continue
        fi
        
        # Skip if completed
        if [ -f "$RUN_DIR/results.json" ]; then
            continue
        fi
        
        # Skip if process is running
        PID=$(ps aux | grep "run_phase2_full.py" | grep "run_$i" | grep "hydrothermal" | grep -v grep | awk '{print $2}')
        if [ -n "$PID" ]; then
            continue
        fi
        
        # Check if log is old (>60 min)
        LOG_FILE="$RUN_DIR/simulation.log"
        if [ -f "$LOG_FILE" ]; then
            if stat -c %Y "$LOG_FILE" >/dev/null 2>&1; then
                AGE_MIN=$(echo "($(date +%s) - $(stat -c %Y "$LOG_FILE")) / 60" | bc)
            else
                AGE_MIN=$(echo "($(date +%s) - $(stat -f %m "$LOG_FILE")) / 60" | bc)
            fi
            
            if [ -n "$AGE_MIN" ] && [ "$AGE_MIN" -gt 60 ]; then
                echo "🔄 Restarting hydrothermal_extended/run_$i (stuck for ${AGE_MIN} min)"
                bash aws_test/scripts/restart_from_checkpoint.sh "$RESULTS_DIR" hydrothermal_extended 2>&1 | grep -A 5 "run_$i" || true
                HYDRO_RESTARTED=$((HYDRO_RESTARTED + 1))
                AVAILABLE_SLOTS=$((AVAILABLE_SLOTS - 1))
                sleep 5  # Small delay between restarts
            fi
        else
            # No log file - restart it
            echo "🔄 Restarting hydrothermal_extended/run_$i (no log file)"
            bash aws_test/scripts/restart_from_checkpoint.sh "$RESULTS_DIR" hydrothermal_extended 2>&1 | grep -A 5 "run_$i" || true
            HYDRO_RESTARTED=$((HYDRO_RESTARTED + 1))
            AVAILABLE_SLOTS=$((AVAILABLE_SLOTS - 1))
            sleep 5
        fi
    done
    
    # Check miller_urey runs (2-4, 9)
    MILLER_RESTARTED=0
    for i in 2 3 4 9; do
        if [ "$AVAILABLE_SLOTS" -le 0 ]; then
            break
        fi
        
        RUN_DIR="$RESULTS_DIR/miller_urey_extended/run_$i"
        if [ ! -d "$RUN_DIR" ]; then
            continue
        fi
        
        # Skip if completed
        if [ -f "$RUN_DIR/results.json" ]; then
            continue
        fi
        
        # Skip if process is running
        PID=$(ps aux | grep "run_phase2_full.py" | grep "run_$i" | grep "miller_urey" | grep -v grep | awk '{print $2}')
        if [ -n "$PID" ]; then
            continue
        fi
        
        # Check if log is old (>60 min)
        LOG_FILE="$RUN_DIR/simulation.log"
        if [ -f "$LOG_FILE" ]; then
            if stat -c %Y "$LOG_FILE" >/dev/null 2>&1; then
                AGE_MIN=$(echo "($(date +%s) - $(stat -c %Y "$LOG_FILE")) / 60" | bc)
            else
                AGE_MIN=$(echo "($(date +%s) - $(stat -f %m "$LOG_FILE")) / 60" | bc)
            fi
            
            if [ -n "$AGE_MIN" ] && [ "$AGE_MIN" -gt 60 ]; then
                echo "🔄 Restarting miller_urey_extended/run_$i (stuck for ${AGE_MIN} min)"
                bash aws_test/scripts/restart_from_checkpoint.sh "$RESULTS_DIR" miller_urey_extended 2>&1 | grep -A 5 "run_$i" || true
                MILLER_RESTARTED=$((MILLER_RESTARTED + 1))
                AVAILABLE_SLOTS=$((AVAILABLE_SLOTS - 1))
                sleep 5
            fi
        else
            # No log file - restart it
            echo "🔄 Restarting miller_urey_extended/run_$i (no log file)"
            bash aws_test/scripts/restart_from_checkpoint.sh "$RESULTS_DIR" miller_urey_extended 2>&1 | grep -A 5 "run_$i" || true
            MILLER_RESTARTED=$((MILLER_RESTARTED + 1))
            AVAILABLE_SLOTS=$((AVAILABLE_SLOTS - 1))
            sleep 5
        fi
    done
    
    if [ "$HYDRO_RESTARTED" -eq 0 ] && [ "$MILLER_RESTARTED" -eq 0 ]; then
        echo "✅ No stuck simulations to restart at this time"
    else
        echo ""
        echo "📊 Restarted: $HYDRO_RESTARTED hydrothermal + $MILLER_RESTARTED miller_urey = $((HYDRO_RESTARTED + MILLER_RESTARTED)) total"
    fi
    
    echo ""
    echo "⏳ Waiting $CHECK_INTERVAL seconds before next check..."
    sleep $CHECK_INTERVAL
done

