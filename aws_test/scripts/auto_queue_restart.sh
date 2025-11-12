#!/bin/bash
# Auto Queue Restart - Phase 2B Miller-Urey
# ==========================================
# 
# Monitors running simulations and automatically starts next batch
# when current batch completes. Runs 4 simulations in parallel max.
#
# Usage:
#   bash auto_queue_restart.sh
#
# This will:
# 1. Monitor runs 5-8 (currently running)
# 2. When they complete, start runs 2-4 + run 10
# 3. When those complete, start runs 11-14
# 4. Continue until all runs complete

set -e

PROJECT_ROOT="$HOME/live2.0"
RESULTS_DIR="$PROJECT_ROOT/results/phase2b_additional/miller_urey_extended"
CONFIG_FILE="$PROJECT_ROOT/aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml"
MAX_PARALLEL=5
MONITOR_INTERVAL=300  # Check every 5 minutes

cd "$PROJECT_ROOT"

echo "================================================================================"
echo "🔄 Auto Queue Restart - Miller-Urey Extended"
echo "================================================================================"
echo ""
echo "Settings:"
echo "  Max parallel: $MAX_PARALLEL"
echo "  Monitor interval: ${MONITOR_INTERVAL}s (5 min)"
echo "  Results dir: $RESULTS_DIR"
echo "  Config: $CONFIG_FILE"
echo ""

# Define queue of runs to execute
# Format: run_id:seed
# Runs 5-8 are already running (seeds 104-107)
# Queue contains: 2-4 (seeds 101-103) and 10-18 (seeds 109-117)
QUEUE=(
    "2:101"
    "3:102"
    "4:103"
    "10:109"
    "11:110"
    "12:111"
    "13:112"
    "14:113"
    "15:114"
    "16:115"
    "17:116"
    "18:117"
)

QUEUE_FILE="$PROJECT_ROOT/logs/phase2b_queue.txt"
LOG_FILE="$PROJECT_ROOT/logs/phase2b_auto_restart.log"
mkdir -p "$PROJECT_ROOT/logs"

# Initialize queue file if doesn't exist
if [ ! -f "$QUEUE_FILE" ]; then
    printf "%s\n" "${QUEUE[@]}" > "$QUEUE_FILE"
    echo "📝 Queue initialized with ${#QUEUE[@]} runs"
fi

# Function to check if a run is completed
is_completed() {
    local run_id=$1
    local results_file="$RESULTS_DIR/run_${run_id}/results.json"
    [ -f "$results_file" ]
}

# Function to check if a run is currently running
is_running() {
    local run_id=$1
    ps aux | grep "run_phase2_full.py" | grep "run_${run_id}" | grep -v grep > /dev/null
}

# Function to count currently running simulations
count_running() {
    ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | wc -l
}

# Function to start a simulation
start_simulation() {
    local run_id=$1
    local seed=$2
    local output_dir="$RESULTS_DIR/run_${run_id}"
    
    mkdir -p "$output_dir"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 🚀 Starting run_${run_id} (seed ${seed})..." | tee -a "$LOG_FILE"
    
    nohup python3 "$PROJECT_ROOT/scripts/run_phase2_full.py" \
        --config "$CONFIG_FILE" \
        --output "$output_dir" \
        --seed "$seed" \
        --steps 500000 \
        --force-cpu \
        >> "$output_dir/simulation.log" 2>&1 &
    
    local pid=$!
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]    ✅ Started with PID $pid" | tee -a "$LOG_FILE"
    
    # Small delay to avoid overwhelming system
    sleep 3
}

# Function to get next items from queue
get_next_from_queue() {
    local count=$1
    local items=()
    
    if [ ! -f "$QUEUE_FILE" ]; then
        echo ""
        return
    fi
    
    local i=0
    while IFS= read -r line && [ $i -lt $count ]; do
        if [ -n "$line" ]; then
            items+=("$line")
            ((i++))
        fi
    done < "$QUEUE_FILE"
    
    printf "%s\n" "${items[@]}"
}

# Function to remove items from queue
remove_from_queue() {
    local count=$1
    
    if [ ! -f "$QUEUE_FILE" ]; then
        return
    fi
    
    # Keep all lines except first $count
    tail -n +$((count + 1)) "$QUEUE_FILE" > "${QUEUE_FILE}.tmp"
    mv "${QUEUE_FILE}.tmp" "$QUEUE_FILE"
}

# Function to get queue size
get_queue_size() {
    if [ ! -f "$QUEUE_FILE" ] || [ ! -s "$QUEUE_FILE" ]; then
        echo "0"
    else
        grep -c "^" "$QUEUE_FILE"
    fi
}

echo "================================================================================"
echo "📊 Initial Status Check"
echo "================================================================================"
echo ""

# Check currently running
CURRENT_RUNNING=$(count_running)
echo "Currently running: $CURRENT_RUNNING simulations"

if [ "$CURRENT_RUNNING" -gt 0 ]; then
    echo ""
    echo "Running processes:"
    ps aux | grep "run_phase2_full.py" | grep "miller_urey_extended" | grep -v grep | awk '{print "  - PID " $2 ": " $NF}'
fi

# Check completed
COMPLETED=0
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18; do
    if is_completed $i; then
        ((COMPLETED++))
    fi
done

QUEUE_SIZE=$(get_queue_size)

echo ""
echo "Progress:"
echo "  ✅ Completed: $COMPLETED / 18"
echo "  🏃 Running: $CURRENT_RUNNING"
echo "  ⏳ In Queue: $QUEUE_SIZE"
echo ""

if [ "$QUEUE_SIZE" -eq 0 ]; then
    echo "✅ Queue is empty. Monitoring current runs until completion..."
    echo ""
fi

echo "================================================================================"
echo "🔍 Monitoring Loop Started"
echo "================================================================================"
echo ""
echo "Press Ctrl+C to stop monitoring"
echo ""

# Main monitoring loop
ITERATION=1
while true; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Iteration $ITERATION - Checking status..." | tee -a "$LOG_FILE"
    
    # Count currently running
    RUNNING=$(count_running)
    QUEUE_SIZE=$(get_queue_size)
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Running: $RUNNING, Queue: $QUEUE_SIZE" | tee -a "$LOG_FILE"
    
    # If we have capacity and queue has items
    if [ "$RUNNING" -lt "$MAX_PARALLEL" ] && [ "$QUEUE_SIZE" -gt 0 ]; then
        CAPACITY=$((MAX_PARALLEL - RUNNING))
        echo "[$(date '+%Y-%m-%d %H:%M:%S')]   🆓 Capacity available: $CAPACITY slots" | tee -a "$LOG_FILE"
        
        # Get next items from queue
        NEXT_ITEMS=$(get_next_from_queue $CAPACITY)
        
        if [ -n "$NEXT_ITEMS" ]; then
            STARTED=0
            while IFS= read -r item; do
                if [ -n "$item" ]; then
                    RUN_ID=$(echo "$item" | cut -d: -f1)
                    SEED=$(echo "$item" | cut -d: -f2)
                    
                    # Check if already completed or running
                    if is_completed "$RUN_ID"; then
                        echo "[$(date '+%Y-%m-%d %H:%M:%S')]   ⏭️  run_${RUN_ID} already completed, skipping" | tee -a "$LOG_FILE"
                    elif is_running "$RUN_ID"; then
                        echo "[$(date '+%Y-%m-%d %H:%M:%S')]   ⏭️  run_${RUN_ID} already running, skipping" | tee -a "$LOG_FILE"
                    else
                        start_simulation "$RUN_ID" "$SEED"
                        ((STARTED++))
                    fi
                fi
            done <<< "$NEXT_ITEMS"
            
            # Remove started items from queue
            if [ "$STARTED" -gt 0 ]; then
                remove_from_queue $STARTED
                echo "[$(date '+%Y-%m-%d %H:%M:%S')]   ✅ Started $STARTED new simulations" | tee -a "$LOG_FILE"
            fi
        fi
    fi
    
    # Check if all done
    COMPLETED=0
    for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18; do
        if is_completed $i; then
            ((COMPLETED++))
        fi
    done
    
    RUNNING=$(count_running)
    QUEUE_SIZE=$(get_queue_size)
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]   📊 Progress: $COMPLETED completed, $RUNNING running, $QUEUE_SIZE queued" | tee -a "$LOG_FILE"
    
    # Exit condition: all completed and none running
    if [ "$COMPLETED" -eq 18 ] && [ "$RUNNING" -eq 0 ]; then
        echo ""
        echo "================================================================================"
        echo "🎉 ALL SIMULATIONS COMPLETED!"
        echo "================================================================================"
        echo "" | tee -a "$LOG_FILE"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✅ All 18 runs completed successfully" | tee -a "$LOG_FILE"
        echo ""
        echo "Next steps:"
        echo "  1. Analyze results: python3 aws_test/scripts/analyze_additional_results.py"
        echo "  2. Check status: python3 aws_test/scripts/check_actual_progress.py"
        echo ""
        exit 0
    fi
    
    # Wait before next check
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]   ⏰ Sleeping for ${MONITOR_INTERVAL}s..." | tee -a "$LOG_FILE"
    echo ""
    sleep $MONITOR_INTERVAL
    ((ITERATION++))
done

