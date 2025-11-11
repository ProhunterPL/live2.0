#!/bin/bash
# Restart stuck simulations from latest checkpoint/snapshot
# ========================================================

RESULTS_DIR="${1:-$HOME/live2.0/results/phase2b_additional}"
SCENARIO="${2:-hydrothermal_extended}"

echo "🔄 Restarting $SCENARIO simulations from checkpoints..."
echo "   Results dir: $RESULTS_DIR"
echo ""

SCENARIO_DIR="$RESULTS_DIR/$SCENARIO"

if [ ! -d "$SCENARIO_DIR" ]; then
    echo "❌ Scenario directory not found: $SCENARIO_DIR"
    exit 1
fi

# Find all run directories
RUN_DIRS=$(find "$SCENARIO_DIR" -maxdepth 1 -type d -name "run_*" | sort)

if [ -z "$RUN_DIRS" ]; then
    echo "❌ No run directories found in $SCENARIO_DIR"
    exit 1
fi

RESTARTED=0
SKIPPED=0

for RUN_DIR in $RUN_DIRS; do
    RUN_NAME=$(basename "$RUN_DIR")
    
    # Check if already completed
    if [ -f "$RUN_DIR/results.json" ]; then
        echo "⏭️  Skipping $RUN_NAME - already completed"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    # Check if process is still running
    PID=$(ps aux | grep "run_phase2_full.py" | grep "$RUN_NAME" | grep -v grep | awk '{print $2}')
    if [ -n "$PID" ]; then
        echo "⏭️  Skipping $RUN_NAME - process still running (PID $PID)"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    # Find latest checkpoint
    CHECKPOINT_DIR="$RUN_DIR/checkpoints"
    LATEST_CHECKPOINT=""
    LATEST_STEP=0
    
    if [ -d "$CHECKPOINT_DIR" ]; then
        # Find latest checkpoint JSON file
        for CP in "$CHECKPOINT_DIR"/checkpoint_*.json; do
            if [ -f "$CP" ]; then
                # Extract step number from filename (checkpoint_00050000.json -> 50000)
                STEP=$(basename "$CP" | sed 's/checkpoint_\([0-9]*\)\.json/\1/' | sed 's/^0*//')
                if [ -z "$STEP" ]; then
                    STEP=0
                fi
                if [ "$STEP" -gt "$LATEST_STEP" ]; then
                    LATEST_STEP=$STEP
                    LATEST_CHECKPOINT="$CP"
                fi
            fi
        done
    fi
    
    # Find latest snapshot as fallback
    SNAPSHOT_DIR="$RUN_DIR/snapshots"
    LATEST_SNAPSHOT=""
    if [ -z "$LATEST_CHECKPOINT" ] && [ -d "$SNAPSHOT_DIR" ]; then
        # Find latest snapshot JSON file
        for SNAP in "$SNAPSHOT_DIR"/snapshot_*.json; do
            if [ -f "$SNAP" ]; then
                # Extract step number
                STEP=$(basename "$SNAP" | sed 's/snapshot_\([0-9]*\)\.json/\1/' | sed 's/^0*//')
                if [ -z "$STEP" ]; then
                    STEP=0
                fi
                if [ "$STEP" -gt "$LATEST_STEP" ]; then
                    LATEST_STEP=$STEP
                    LATEST_SNAPSHOT="$SNAP"
                fi
            fi
        done
    fi
    
    if [ "$LATEST_STEP" -eq 0 ]; then
        echo "⚠️  $RUN_NAME: No checkpoint/snapshot found, will restart from beginning"
        LATEST_STEP=0
    else
        echo "📦 $RUN_NAME: Found checkpoint at step $LATEST_STEP"
    fi
    
    # Extract config file from log or use default
    CONFIG_FILE=""
    if [ -f "$RUN_DIR/simulation.log" ]; then
        CONFIG_MATCH=$(grep -o "\--config [^ ]*" "$RUN_DIR/simulation.log" | head -1 | awk '{print $2}')
        if [ -n "$CONFIG_MATCH" ] && [ -f "$CONFIG_MATCH" ]; then
            CONFIG_FILE="$CONFIG_MATCH"
        fi
    fi
    
    # Default config if not found
    if [ -z "$CONFIG_FILE" ]; then
        CONFIG_FILE="$HOME/live2.0/aws_test/configs/phase2_${SCENARIO}_SUPER_FAST.yaml"
        if [ ! -f "$CONFIG_FILE" ]; then
            CONFIG_FILE="$HOME/live2.0/aws_test/configs/phase2_${SCENARIO}.yaml"
        fi
    fi
    
    if [ ! -f "$CONFIG_FILE" ]; then
        echo "❌ $RUN_NAME: Config file not found, skipping"
        continue
    fi
    
    # Extract seed from log or use default
    SEED=$(grep -o "\--seed [0-9]*" "$RUN_DIR/simulation.log" 2>/dev/null | head -1 | awk '{print $2}')
    if [ -z "$SEED" ]; then
        # Extract from run number
        RUN_NUM=$(echo "$RUN_NAME" | sed 's/run_//')
        SEED=$((110 + RUN_NUM - 1))  # Default seed for hydrothermal
    fi
    
    # Calculate remaining steps (assuming 500K total)
    MAX_STEPS=500000
    REMAINING_STEPS=$((MAX_STEPS - LATEST_STEP))
    
    if [ "$REMAINING_STEPS" -le 0 ]; then
        echo "✅ $RUN_NAME: Already completed (step $LATEST_STEP >= $MAX_STEPS)"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    echo "🚀 Restarting $RUN_NAME from step $LATEST_STEP (remaining: $REMAINING_STEPS steps)"
    echo "   Config: $CONFIG_FILE"
    echo "   Seed: $SEED"
    echo ""
    
    # Note: run_phase2_full.py doesn't support --checkpoint yet
    # For now, we'll restart from beginning but log the checkpoint info
    # TODO: Add checkpoint loading to run_phase2_full.py
    
    cd "$HOME/live2.0" || exit 1
    
    # Run in background with nohup
    nohup python3 scripts/run_phase2_full.py \
        --config "$CONFIG_FILE" \
        --output "$RUN_DIR" \
        --seed "$SEED" \
        --steps "$MAX_STEPS" \
        --force-cpu \
        >> "$RUN_DIR/simulation.log" 2>&1 &
    
    NEW_PID=$!
    echo "   ✅ Started with PID $NEW_PID"
    echo "   📝 Note: Restarting from beginning (checkpoint loading not yet implemented)"
    echo "   💡 Checkpoint at step $LATEST_STEP is available but not used"
    echo ""
    
    RESTARTED=$((RESTARTED + 1))
    
    # Small delay to avoid overwhelming the system
    sleep 2
done

echo ""
echo "📊 Summary:"
echo "   Restarted: $RESTARTED simulations"
echo "   Skipped: $SKIPPED simulations (completed or running)"
echo ""
echo "⚠️  NOTE: Checkpoint loading is not yet implemented in run_phase2_full.py"
echo "   Simulations will restart from beginning, but checkpoints are preserved"
echo ""

