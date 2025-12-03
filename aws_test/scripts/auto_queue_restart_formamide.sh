#!/bin/bash
# Auto Queue & Restart System - FORMAMIDE
# ===========================================
# 
# Based on successful hydrothermal auto_queue_restart_hydro.sh
# Manages formamide simulations with automatic restart
# 
# Features:
# - Max 4 parallel simulations
# - Auto-restart on completion
# - Progress monitoring
# - Stuck detection
# 
# Usage:
#   nohup bash aws_test/scripts/auto_queue_restart_formamide.sh > logs/auto_restart_formamide.log 2>&1 &

set -e

# ============================================================================
# CONFIGURATION
# ============================================================================

# Base configuration
BASE_DIR="$HOME/live2.0"
RESULTS_BASE="results/phase2b_additional/formamide_extended"
CONFIG="aws_test/configs/phase2_formamide_extended_SUPER_FAST.yaml"
MAX_PARALLEL=4
CHECK_INTERVAL=300  # 5 minutes

# Run configuration - 8 runs (run_1 to run_8)
TOTAL_RUNS=8
RUNS=($(seq 1 $TOTAL_RUNS))
SEEDS=(200 201 202 203 204 205 206 207)  # run_1=200, run_2=201, etc.

# Log files
MAIN_LOG="$BASE_DIR/logs/auto_restart_formamide_main.log"
PROGRESS_LOG="$BASE_DIR/logs/auto_restart_formamide_progress.log"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MAIN_LOG"
}

log_progress() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$PROGRESS_LOG"
}

get_running_count() {
    ps aux | grep "run_phase2_full.py" | grep "formamide" | grep -v grep | wc -l
}

is_run_completed() {
    local run_num=$1
    local run_dir="$BASE_DIR/$RESULTS_BASE/run_$run_num"
    
    if [ -f "$run_dir/results.json" ]; then
        return 0  # Completed
    else
        return 1  # Not completed
    fi
}

is_run_in_progress() {
    local run_num=$1
    ps aux | grep "run_phase2_full.py" | grep "formamide" | grep "run_$run_num" | grep -v grep > /dev/null
    return $?
}

start_run() {
    local run_num=$1
    local seed=$2
    local run_dir="$BASE_DIR/$RESULTS_BASE/run_$run_num"
    local log_file="$BASE_DIR/logs/formamide_run_${run_num}.log"
    
    log "[START] Starting run_$run_num (seed=$seed)"
    
    cd "$BASE_DIR"
    mkdir -p "$run_dir"
    mkdir -p "$(dirname "$log_file")"
    
    # Use 16 threads per run (optimal for 4 parallel runs on 64-core system)
    CPU_THREADS=16
    
    nohup python3 "$BASE_DIR/scripts/run_phase2_full.py" \
        --config "$BASE_DIR/$CONFIG" \
        --output "$run_dir" \
        --seed "$seed" \
        --steps 500000 \
        --force-cpu \
        --cpu-threads "$CPU_THREADS" \
        >> "$log_file" 2>&1 &
    
    local pid=$!
    log "[START] ✅ Started run_$run_num with PID $pid (seed=$seed)"
    log_progress "run_$run_num:STARTED:$(date '+%Y-%m-%d %H:%M:%S'):PID=$pid"
    
    # Small delay to avoid overwhelming system
    sleep 5
}

get_next_run_to_start() {
    for run_num in "${RUNS[@]}"; do
        if ! is_run_completed "$run_num" && ! is_run_in_progress "$run_num"; then
            echo "$run_num"
            return 0
        fi
    done
    echo ""
    return 1
}

count_completed() {
    local count=0
    for run_num in "${RUNS[@]}"; do
        if is_run_completed "$run_num"; then
            ((count++))
        fi
    done
    echo "$count"
}

check_stuck_processes() {
    # Check for processes that haven't updated log in last hour
    local stuck_count=0
    for run_num in "${RUNS[@]}"; do
        if is_run_in_progress "$run_num"; then
            local log_file="$BASE_DIR/logs/formamide_run_${run_num}.log"
            if [ -f "$log_file" ]; then
                local last_update=$(stat -c %Y "$log_file" 2>/dev/null || echo 0)
                local now=$(date +%s)
                local age=$((now - last_update))
                
                # If log hasn't updated in 2 hours, consider stuck
                if [ $age -gt 7200 ]; then
                    log "[WARNING] run_$run_num appears stuck (log age: ${age}s)"
                    ((stuck_count++))
                fi
            fi
        fi
    done
    
    if [ $stuck_count -gt 0 ]; then
        log "[WARNING] Found $stuck_count potentially stuck processes"
    fi
}

print_status() {
    local running=$(get_running_count)
    local completed=$(count_completed)
    local pending=$((TOTAL_RUNS - completed - running))
    
    log "=========================================="
    log "FORMAMIDE QUEUE STATUS"
    log "=========================================="
    log "Running: $running / $MAX_PARALLEL"
    log "Completed: $completed / $TOTAL_RUNS"
    log "Pending: $pending"
    log "=========================================="
    
    # Show which runs are in each state
    log "Completed runs:"
    for run_num in "${RUNS[@]}"; do
        if is_run_completed "$run_num"; then
            log "  ✅ run_$run_num"
        fi
    done
    
    log "Running runs:"
    for run_num in "${RUNS[@]}"; do
        if is_run_in_progress "$run_num"; then
            log "  🔄 run_$run_num"
        fi
    done
    
    log "Pending runs:"
    for run_num in "${RUNS[@]}"; do
        if ! is_run_completed "$run_num" && ! is_run_in_progress "$run_num"; then
            log "  ⏳ run_$run_num"
        fi
    done
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

cd "$BASE_DIR"

# Create logs directory
mkdir -p "$BASE_DIR/logs"
mkdir -p "$BASE_DIR/$RESULTS_BASE"

log "=========================================="
log "FORMAMIDE AUTO-RESTART SYSTEM"
log "=========================================="
log "Configuration:"
log "  Base dir: $BASE_DIR"
log "  Results: $RESULTS_BASE"
log "  Config: $CONFIG"
log "  Total runs: $TOTAL_RUNS"
log "  Max parallel: $MAX_PARALLEL"
log "  Check interval: ${CHECK_INTERVAL}s"
log "=========================================="

# Initial status
print_status

# Main loop
ITERATION=1
while true; do
    log "----------------------------------------"
    log "[CHECK] Iteration $ITERATION - Checking queue status..."
    
    # Count running processes
    running=$(get_running_count)
    log "[INFO] Currently running: $running / $MAX_PARALLEL"
    
    # Start new runs if slots available
    if [ $running -lt $MAX_PARALLEL ]; then
        slots_available=$((MAX_PARALLEL - running))
        log "[INFO] Slots available: $slots_available"
        
        for ((i=1; i<=slots_available; i++)); do
            next_run=$(get_next_run_to_start)
            
            if [ -n "$next_run" ]; then
                # Get seed for this run (array is 0-indexed)
                seed_index=$((next_run - 1))
                seed=${SEEDS[$seed_index]}
                
                log "[QUEUE] Starting next run: run_$next_run (seed=$seed)"
                start_run "$next_run" "$seed"
                sleep 5  # Small delay between starts
            else
                log "[INFO] No more runs to start"
                break
            fi
        done
    else
        log "[INFO] All slots full, waiting..."
    fi
    
    # Check for stuck processes
    check_stuck_processes
    
    # Print status
    print_status
    
    # Check if all completed
    completed=$(count_completed)
    if [ $completed -eq $TOTAL_RUNS ]; then
        log "=========================================="
        log "✅ ALL FORMAMIDE RUNS COMPLETED!"
        log "=========================================="
        log "Total runs: $TOTAL_RUNS"
        log "Completed: $completed"
        log "=========================================="
        break
    fi
    
    # Wait before next check
    log "[WAIT] Waiting ${CHECK_INTERVAL}s before next check..."
    sleep "$CHECK_INTERVAL"
    ((ITERATION++))
done

log "[DONE] Formamide queue system finished"

