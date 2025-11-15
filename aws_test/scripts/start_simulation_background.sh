#!/bin/bash
# Start simulation in background with proper logging
# ==================================================
# Usage: bash start_simulation_background.sh <run_id> <seed>

set -e

if [ $# -lt 2 ]; then
    echo "Usage: $0 <run_id> <seed>"
    echo "Example: $0 3 103"
    exit 1
fi

RUN_ID=$1
SEED=$2

PROJECT_ROOT="$HOME/live2.0"
cd "$PROJECT_ROOT"

CONFIG_FILE="$PROJECT_ROOT/aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml"
OUTPUT_DIR="$PROJECT_ROOT/results/phase2b_additional/miller_urey_extended/run_${RUN_ID}"
SCRIPT_FILE="$PROJECT_ROOT/scripts/run_phase2_full.py"

mkdir -p "$OUTPUT_DIR"

echo "================================================================================"
echo "🚀 Starting Simulation in Background: run_${RUN_ID}"
echo "================================================================================"
echo ""
echo "Config: $CONFIG_FILE"
echo "Output: $OUTPUT_DIR"
echo "Seed: $SEED"
echo ""

# Validate files
if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

if [ ! -f "$SCRIPT_FILE" ]; then
    echo "❌ ERROR: Script not found: $SCRIPT_FILE"
    exit 1
fi

# Start in background
echo "Starting process..."
nohup python3 "$SCRIPT_FILE" \
    --config "$CONFIG_FILE" \
    --output "$OUTPUT_DIR" \
    --seed "$SEED" \
    --steps 500000 \
    --force-cpu \
    >> "$OUTPUT_DIR/simulation.log" 2>&1 &

PID=$!
echo "✅ Started with PID: $PID"
echo ""
echo "Monitor with:"
echo "  tail -f $OUTPUT_DIR/simulation.log"
echo "  ps aux | grep $PID"
echo ""
echo "Check status:"
echo "  python3 aws_test/scripts/check_actual_progress.py"
echo ""

# Wait a moment and check if still running
sleep 3
if ps -p $PID > /dev/null 2>&1; then
    echo "✅ Process is running!"
else
    echo "❌ Process died immediately!"
    echo "📋 Check log:"
    tail -50 "$OUTPUT_DIR/simulation.log"
    exit 1
fi

