#!/bin/bash
# Restart a specific run with limited CPU threads

if [ $# -lt 2 ]; then
    echo "Usage: $0 <run_number> <seed> [cpu_threads]"
    echo "Example: $0 4 103 12"
    echo ""
    echo "If cpu_threads not specified, will use 12 (optimal for 5 parallel runs)"
    exit 1
fi

RUN_NUM=$1
SEED=$2
CPU_THREADS=${3:-12}  # Default to 12 threads

RESULTS_DIR="$HOME/live2.0/results/phase2b_additional/miller_urey_extended"
CONFIG_FILE="$HOME/live2.0/aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml"
RUN_DIR="$RESULTS_DIR/run_$RUN_NUM"

echo "=================================================================================="
echo "🔄 RESTART RUN_$RUN_NUM WITH LIMITED CPU THREADS"
echo "=================================================================================="
echo ""
echo "📊 Configuration:"
echo "   Run: run_$RUN_NUM"
echo "   Seed: $SEED"
echo "   CPU threads: $CPU_THREADS"
echo "   Config: $CONFIG_FILE"
echo ""

# Check if run is already running
PID=$(ps aux | grep "run_$RUN_NUM" | grep run_phase2_full | grep -v grep | awk '{print $2}')
if [ -n "$PID" ]; then
    echo "⚠️  Run $RUN_NUM is already running (PID: $PID)"
    read -p "Kill and restart? (yes/no): " CONFIRM
    if [ "$CONFIRM" != "yes" ]; then
        echo "❌ Cancelled"
        exit 1
    fi
    echo "   Killing PID $PID..."
    kill $PID 2>/dev/null
    sleep 2
    if ps -p $PID > /dev/null 2>&1; then
        kill -9 $PID 2>/dev/null
    fi
    echo "   ✅ Process killed"
    echo ""
fi

# Check if completed
if [ -f "$RUN_DIR/results.json" ]; then
    echo "⚠️  Run $RUN_NUM is already COMPLETED (results.json exists)"
    read -p "Restart anyway? (yes/no): " CONFIRM
    if [ "$CONFIRM" != "yes" ]; then
        echo "❌ Cancelled"
        exit 1
    fi
    echo ""
fi

# Create output directory
mkdir -p "$RUN_DIR"

echo "🚀 Starting run_$RUN_NUM..."
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
    echo "✅ Process is running (CPU: ${CPU_USAGE}%)"
    echo ""
    echo "📝 Monitor with:"
    echo "   tail -f $RUN_DIR/simulation_restart.log"
    echo "   ps aux | grep $NEW_PID"
else
    echo "❌ Process may have crashed - check log:"
    echo "   tail -50 $RUN_DIR/simulation_restart.log"
fi

echo ""
echo "=================================================================================="

