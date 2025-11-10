#!/bin/bash
# Restart Phase 2B with SAFER configuration (no cluster detection)
# =================================================================

set -e

echo "================================================================================"
echo "🔄 RESTARTING PHASE 2B WITH SAFER CONFIGURATION"
echo "================================================================================"
echo ""

# Configuration
RESULTS_DIR="$HOME/live2.0/results/phase2b_additional"
CONFIG_FILE="$HOME/live2.0/aws_test/configs/phase2_miller_urey_extended_SAFER.yaml"
LOG_DIR="$HOME/live2.0/logs/phase2b_safe"

# Create log directory
mkdir -p "$LOG_DIR"

# Check how many runs to do
if [ -z "$1" ]; then
    echo "Usage: $0 <number_of_runs> [starting_seed]"
    echo ""
    echo "Example: $0 10 110    # Run 10 simulations starting from seed 110"
    echo ""
    exit 1
fi

NUM_RUNS=$1
START_SEED=${2:-110}  # Default starting seed = 110

echo "Configuration:"
echo "  Config file: $CONFIG_FILE"
echo "  Number of runs: $NUM_RUNS"
echo "  Starting seed: $START_SEED"
echo "  Results dir: $RESULTS_DIR"
echo "  Log dir: $LOG_DIR"
echo ""

if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

echo "⚠️  This will start $NUM_RUNS new simulations in the background."
read -p "Continue? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 1
fi

echo ""
echo "🚀 Starting simulations..."
echo ""

# Track PIDs
PIDS=()

for ((i=0; i<NUM_RUNS; i++)); do
    RUN_NUM=$((10 + i))  # Start from run_10
    SEED=$((START_SEED + i))
    OUTPUT_DIR="$RESULTS_DIR/miller_urey_extended/run_$RUN_NUM"
    LOG_FILE="$LOG_DIR/run_${RUN_NUM}.log"
    
    echo "Starting run_$RUN_NUM (seed=$SEED)..."
    
    # Run with -u for unbuffered output and nohup for persistence
    nohup python3 -u $HOME/live2.0/scripts/run_phase2_full.py \
        --config "$CONFIG_FILE" \
        --output "$OUTPUT_DIR" \
        --seed $SEED \
        --steps 500000 \
        --force-cpu \
        > "$LOG_FILE" 2>&1 &
    
    PID=$!
    PIDS+=($PID)
    
    echo "  ✅ Started (PID=$PID, log=$LOG_FILE)"
    
    # Small delay to avoid race conditions
    sleep 2
done

echo ""
echo "================================================================================"
echo "✅ ALL SIMULATIONS STARTED"
echo "================================================================================"
echo ""
echo "Running PIDs: ${PIDS[@]}"
echo ""
echo "Monitor progress with:"
echo "  python3 ~/live2.0/aws_test/scripts/check_actual_progress.py"
echo ""
echo "Check logs:"
echo "  tail -f $LOG_DIR/run_10.log"
echo ""
echo "Kill all if needed:"
echo "  kill ${PIDS[@]}"
echo ""

