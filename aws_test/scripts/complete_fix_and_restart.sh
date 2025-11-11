#!/bin/bash
# Complete fix: Kill stuck + Restart with proper config
# =====================================================

set -e

RESULTS_DIR="${1:-$HOME/live2.0/results/phase2b_additional}"
MAX_PARALLEL=4

echo "================================================================================"
echo "🔧 COMPLETE FIX AND RESTART"
echo "================================================================================"
echo ""
echo "This will:"
echo "  1. Identify all stuck simulations"
echo "  2. Kill stuck processes (if any running)"
echo "  3. Limit parallel simulations to $MAX_PARALLEL"
echo "  4. Restart stuck hydrothermal_extended runs (with fixed config)"
echo "  5. Restart stuck miller_urey_extended runs"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 1
fi

cd ~/live2.0 || exit 1

# Make scripts executable
chmod +x aws_test/scripts/*.sh

echo ""
echo "================================================================================"
echo "STEP 1: Identify Running Processes"
echo "================================================================================"
echo ""

bash aws_test/scripts/identify_running_processes.sh

echo ""
echo "================================================================================"
echo "STEP 2: Kill All Stuck Simulations (>60 min no activity)"
echo "================================================================================"
echo ""

bash aws_test/scripts/kill_all_stuck.sh 60 "$RESULTS_DIR"

echo ""
echo "================================================================================"
echo "STEP 3: Limit Parallel Simulations to $MAX_PARALLEL"
echo "================================================================================"
echo ""

bash aws_test/scripts/limit_parallel_simulations.sh $MAX_PARALLEL

echo ""
echo "================================================================================"
echo "STEP 4: Restart Stuck Hydrothermal Runs"
echo "================================================================================"
echo ""

# Count how many hydrothermal runs need restart
HYDRO_NEED_RESTART=0
for i in {1..10}; do
    RUN_DIR="$RESULTS_DIR/hydrothermal_extended/run_$i"
    if [ -d "$RUN_DIR" ] && [ ! -f "$RUN_DIR/results.json" ]; then
        # Check if process is running
        PID=$(ps aux | grep "run_phase2_full.py" | grep "run_$i" | grep "hydrothermal" | grep -v grep | awk '{print $2}')
        if [ -z "$PID" ]; then
            HYDRO_NEED_RESTART=$((HYDRO_NEED_RESTART + 1))
        fi
    fi
done

echo "Found $HYDRO_NEED_RESTART hydrothermal runs that need restart"
echo ""

if [ "$HYDRO_NEED_RESTART" -gt 0 ]; then
    # Check how many are currently running
    CURRENT_RUNNING=$(ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l | tr -d ' ')
    AVAILABLE_SLOTS=$((MAX_PARALLEL - CURRENT_RUNNING))
    
    echo "Currently running: $CURRENT_RUNNING"
    echo "Available slots: $AVAILABLE_SLOTS"
    echo ""
    
    if [ "$AVAILABLE_SLOTS" -le 0 ]; then
        echo "⚠️  No available slots. Waiting for some simulations to complete..."
        echo "   You can restart hydrothermal runs manually later:"
        echo "   bash aws_test/scripts/restart_from_checkpoint.sh $RESULTS_DIR hydrothermal_extended"
    else
        echo "🔄 Restarting up to $AVAILABLE_SLOTS hydrothermal runs..."
        bash aws_test/scripts/restart_from_checkpoint.sh "$RESULTS_DIR" hydrothermal_extended
    fi
else
    echo "✅ No hydrothermal runs need restart"
fi

echo ""
echo "================================================================================"
echo "STEP 5: Restart Stuck Miller-Urey Runs"
echo "================================================================================"
echo ""

# Count how many miller_urey runs need restart (runs 2-4, 9)
MILLER_NEED_RESTART=0
for i in 2 3 4 9; do
    RUN_DIR="$RESULTS_DIR/miller_urey_extended/run_$i"
    if [ -d "$RUN_DIR" ] && [ ! -f "$RUN_DIR/results.json" ]; then
        # Check if process is running
        PID=$(ps aux | grep "run_phase2_full.py" | grep "run_$i" | grep "miller_urey" | grep -v grep | awk '{print $2}')
        if [ -z "$PID" ]; then
            MILLER_NEED_RESTART=$((MILLER_NEED_RESTART + 1))
        fi
    fi
done

echo "Found $MILLER_NEED_RESTART miller_urey runs that need restart (runs 2-4, 9)"
echo ""

if [ "$MILLER_NEED_RESTART" -gt 0 ]; then
    # Check how many are currently running
    CURRENT_RUNNING=$(ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l | tr -d ' ')
    AVAILABLE_SLOTS=$((MAX_PARALLEL - CURRENT_RUNNING))
    
    echo "Currently running: $CURRENT_RUNNING"
    echo "Available slots: $AVAILABLE_SLOTS"
    echo ""
    
    if [ "$AVAILABLE_SLOTS" -le 0 ]; then
        echo "⚠️  No available slots. Waiting for some simulations to complete..."
        echo "   You can restart miller_urey runs manually later:"
        echo "   bash aws_test/scripts/restart_from_checkpoint.sh $RESULTS_DIR miller_urey_extended"
    else
        echo "🔄 Restarting up to $AVAILABLE_SLOTS miller_urey runs..."
        # Note: restart_from_checkpoint.sh will skip runs that are already running
        bash aws_test/scripts/restart_from_checkpoint.sh "$RESULTS_DIR" miller_urey_extended
    fi
else
    echo "✅ No miller_urey runs need restart"
fi

echo ""
echo "================================================================================"
echo "✅ FIX COMPLETE"
echo "================================================================================"
echo ""
echo "Current status:"
CURRENT_RUNNING=$(ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l | tr -d ' ')
echo "   Running simulations: $CURRENT_RUNNING / $MAX_PARALLEL"
echo ""
echo "Monitor progress:"
echo "   python3 aws_test/scripts/check_actual_progress.py"
echo ""
echo "If you need to restart more simulations later:"
echo "   bash aws_test/scripts/restart_from_checkpoint.sh $RESULTS_DIR <scenario>"
echo ""

