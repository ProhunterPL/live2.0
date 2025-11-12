#!/bin/bash
# Simple Progress Monitor - Phase 2B Miller-Urey
# Checks progress every 5 minutes without restarting

while true; do
    echo "============================================================"
    echo "=== $(date '+%Y-%m-%d %H:%M:%S') ==="
    echo "============================================================"
    
    RUNNING=$(ps aux | grep 'run_phase2_full.py' | grep -v grep | wc -l)
    echo "🏃 Running: $RUNNING processes"
    echo ""
    
    for run in 1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18; do
        RUN_DIR="results/phase2b_additional/miller_urey_extended/run_$run"
        
        if [ -f "$RUN_DIR/results.json" ]; then
            echo "  ✅ run_$run: COMPLETED"
        elif [ -f "$RUN_DIR/simulation_restart.log" ]; then
            # Try restart log first
            LAST=$(grep "Step.*000/500,000" "$RUN_DIR/simulation_restart.log" 2>/dev/null | tail -1)
            if [ -n "$LAST" ]; then
                echo "  🏃 run_$run: $LAST"
            else
                echo "  🔄 run_$run: Starting..."
            fi
        elif [ -f "$RUN_DIR/simulation.log" ]; then
            # Fallback to main log
            LAST=$(grep "Step.*000/500,000" "$RUN_DIR/simulation.log" 2>/dev/null | tail -1)
            if [ -n "$LAST" ]; then
                echo "  🏃 run_$run: $LAST (from old log)"
            else
                echo "  ⏸️  run_$run: No progress yet"
            fi
        else
            echo "  ⏸️  run_$run: Not started"
        fi
    done
    
    echo ""
    echo "⏰ Next check in 5 minutes..."
    echo ""
    sleep 300
done

