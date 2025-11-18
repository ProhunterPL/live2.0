#!/bin/bash
# Quick restart script for run_4 with limited CPU threads

RUN_NUM=4
SEED=103
CPU_THREADS=12

RESULTS_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended"
CONFIG_FILE="$HOME/live2.0/aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml"
RUN_DIR="$RESULTS_DIR/run_$RUN_NUM"

echo "=================================================================================="
echo "🔄 RESTARTING RUN_4 WITH LIMITED CPU THREADS"
echo "=================================================================================="
echo ""
echo "📊 Configuration:"
echo "   Run: run_$RUN_NUM"
echo "   Seed: $SEED"
echo "   CPU threads: $CPU_THREADS"
echo "   Last step: 277500/500000 (55%)"
echo ""

# Check if run is already running and kill it
PID=$(ps aux | grep "run_$RUN_NUM" | grep run_phase2_full | grep -v grep | awk '{print $2}')
if [ -n "$PID" ]; then
    echo "⚠️  Found existing process (PID: $PID) - killing..."
    kill $PID 2>/dev/null
    sleep 2
    if ps -p $PID > /dev/null 2>&1; then
        echo "   Force killing..."
        kill -9 $PID 2>/dev/null
        sleep 1
    fi
    echo "   ✅ Process killed"
    echo ""
fi

# Check if completed
if [ -f "$RUN_DIR/results.json" ]; then
    echo "⚠️  Run 4 is already COMPLETED (results.json exists)"
    echo "   Skipping restart"
    exit 0
fi

# Create output directory
mkdir -p "$RUN_DIR"

echo "🚀 Starting run_4..."
echo ""

cd "$HOME/live2.0"

# Start simulation
setsid nohup python3 scripts/run_phase2_full.py \
    --config "$CONFIG_FILE" \
    --output "$RUN_DIR" \
    --seed "$SEED" \
    --steps 500000 \
    --force-cpu \
    --cpu-threads "$CPU_THREADS" \
    >> "$RUN_DIR/simulation_restart.log" 2>&1 < /dev/null &

NEW_PID=$!
echo "✅ Started with PID: $NEW_PID"
echo ""

# Wait a moment to verify it started
sleep 3
if ps -p $NEW_PID > /dev/null 2>&1; then
    CPU_USAGE=$(ps -p $NEW_PID -o %cpu= 2>/dev/null | tr -d ' ')
    THREADS=$(ps -p $NEW_PID -o nlwp= 2>/dev/null | tr -d ' ')
    echo "✅ Process is running!"
    echo "   CPU: ${CPU_USAGE}%"
    echo "   Threads: ${THREADS}"
    echo ""
    echo "📝 Monitor with:"
    echo "   tail -f $RUN_DIR/simulation_restart.log"
    echo "   ps aux | grep $NEW_PID"
    echo ""
    echo "📊 Check progress:"
    echo "   python3 ~/live2.0/aws_test/scripts/check_actual_progress.py | grep run_4"
else
    echo "❌ Process may have crashed - check log:"
    echo "   tail -50 $RUN_DIR/simulation_restart.log"
    exit 1
fi

echo "=================================================================================="

