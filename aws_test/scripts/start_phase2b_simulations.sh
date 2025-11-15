#!/bin/bash
# Start Phase 2B simulations - Complete Solution
# ================================================
# This script will:
# 1. Check current status
# 2. Start auto_queue_restart if not running
# 3. Or start simulations manually if needed

set -e

PROJECT_ROOT="$HOME/live2.0"
cd "$PROJECT_ROOT"

echo "================================================================================"
echo "🚀 STARTING PHASE 2B SIMULATIONS"
echo "================================================================================"
echo ""

# Check current status
echo "📊 Checking current status..."
RUNNING=$(ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | wc -l | tr -d ' ')
AUTO_RESTART=$(ps aux | grep "auto_queue_restart.sh" | grep -v grep | wc -l | tr -d ' ')

echo "  Currently running simulations: $RUNNING"
echo "  Auto-restart script running: $AUTO_RESTART"
echo ""

# Check completed runs
COMPLETED=0
RESULTS_DIR="$PROJECT_ROOT/results/phase2b_additional/miller_urey_extended"
for i in {1..18}; do
    if [ -f "$RESULTS_DIR/run_${i}/results.json" ]; then
        ((COMPLETED++))
    fi
done
echo "  Completed runs: $COMPLETED / 18"
echo ""

# Option 1: Start auto_queue_restart if not running
if [ "$AUTO_RESTART" -eq 0 ]; then
    echo "🔄 Auto-restart script is NOT running"
    echo ""
    
    # Default to option 2 (manual start) if no input
    OPTION="${1:-2}"
    
    if [ "$OPTION" = "1" ] || [ "$OPTION" = "auto" ]; then
        echo ""
        echo "Starting auto_queue_restart..."
        
        # Make sure script is executable
        chmod +x aws_test/scripts/auto_queue_restart.sh
        
        # Start in background
        nohup bash aws_test/scripts/auto_queue_restart.sh > logs/auto_restart_main.log 2>&1 &
        AUTO_PID=$!
        
        echo "✅ Started auto_queue_restart with PID: $AUTO_PID"
        echo ""
        echo "Monitor with:"
        echo "  tail -f logs/auto_restart_main.log"
        echo "  tail -f logs/phase2b_auto_restart.log"
        echo "  python3 aws_test/scripts/check_actual_progress.py"
        echo ""
        
        # Wait a moment and check if it started simulations
        sleep 5
        NEW_RUNNING=$(ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | wc -l | tr -d ' ')
        if [ "$NEW_RUNNING" -gt "$RUNNING" ]; then
            echo "✅ Simulations started! ($RUNNING -> $NEW_RUNNING running)"
        else
            echo "⏳ Waiting for simulations to start (check logs above)"
        fi
        
    elif [ "$OPTION" = "2" ]; then
        echo ""
        echo "Starting simulations manually..."
        
        # Start first 5 runs that are not completed
        STARTED=0
        for run_id in 2 3 4 5 10; do
            if [ ! -f "$RESULTS_DIR/run_${run_id}/results.json" ]; then
                # Get seed (run_id + 100)
                seed=$((run_id + 99))
                
                echo "  Starting run_${run_id} (seed ${seed})..."
                
                bash aws_test/scripts/start_simulation_background.sh $run_id $seed
                
                ((STARTED++))
                sleep 2
                
                # Don't exceed MAX_PARALLEL
                if [ "$STARTED" -ge 5 ]; then
                    break
                fi
            fi
        done
        
        echo ""
        echo "✅ Started $STARTED simulations"
        echo ""
        echo "Monitor with:"
        echo "  ps aux | grep run_phase2_full.py | grep -v grep"
        echo "  python3 aws_test/scripts/check_actual_progress.py"
        
    else
        echo "Invalid option. Exiting."
        exit 1
    fi
    
else
    echo "✅ Auto-restart script is already running"
    echo ""
    echo "Current status:"
    echo "  Running: $RUNNING simulations"
    echo "  Completed: $COMPLETED / 18"
    echo ""
    echo "Monitor with:"
    echo "  tail -f logs/phase2b_auto_restart.log"
    echo "  python3 aws_test/scripts/check_actual_progress.py"
fi

echo ""
echo "================================================================================"
echo "✅ DONE"
echo "================================================================================"

