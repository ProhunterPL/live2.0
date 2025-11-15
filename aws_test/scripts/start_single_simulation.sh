#!/bin/bash
# Start a single simulation with full debugging
# ============================================
# Usage: bash start_single_simulation.sh <run_id> <seed>

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
LOG_FILE="$OUTPUT_DIR/startup.log"

mkdir -p "$OUTPUT_DIR"

echo "================================================================================"
echo "🚀 Starting Simulation: run_${RUN_ID}"
echo "================================================================================"
echo ""
echo "Config: $CONFIG_FILE"
echo "Output: $OUTPUT_DIR"
echo "Seed: $SEED"
echo "Steps: 500000"
echo ""

# Check if config exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# Check if script exists
SCRIPT_FILE="$PROJECT_ROOT/scripts/run_phase2_full.py"
if [ ! -f "$SCRIPT_FILE" ]; then
    echo "❌ ERROR: Script not found: $SCRIPT_FILE"
    exit 1
fi

# Test Python first
echo "Testing Python..."
python3 --version || {
    echo "❌ ERROR: Python3 not found!"
    exit 1
}

# Test imports
echo "Testing imports..."
python3 << 'PYEOF'
import sys
sys.path.insert(0, '/home/ubuntu/live2.0')
try:
    import taichi as ti
    from backend.sim.config import SimulationConfig
    print("✅ Imports OK")
except Exception as e:
    print(f"❌ Import failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
PYEOF

if [ $? -ne 0 ]; then
    echo "❌ ERROR: Python imports failed!"
    exit 1
fi

echo ""
echo "Starting simulation..."
echo ""

# Start simulation - run in foreground first to see errors
python3 "$SCRIPT_FILE" \
    --config "$CONFIG_FILE" \
    --output "$OUTPUT_DIR" \
    --seed "$SEED" \
    --steps 500000 \
    --force-cpu \
    2>&1 | tee "$LOG_FILE"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "✅ Simulation completed successfully!"
else
    echo ""
    echo "❌ Simulation failed with exit code: $EXIT_CODE"
    echo "📋 Check log: $LOG_FILE"
    echo "📋 Check simulation log: $OUTPUT_DIR/simulation.log"
fi

exit $EXIT_CODE

