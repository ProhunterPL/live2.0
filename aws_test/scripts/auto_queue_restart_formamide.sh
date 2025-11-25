#!/bin/bash
# Auto Queue & Restart System - FORMAMIDE
<<<<<<< HEAD
# =======================================
=======
# ===========================================
>>>>>>> 841e48d (koniec testów hydro, próba wycigniecia danych autokatalizy)
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
<<<<<<< HEAD
CONFIG="aws_test/configs/phase2_formamide_AWS_OPTIMIZED.yaml"
MAX_PARALLEL=4
CHECK_INTERVAL=300  # 5 minutes

# Run configuration - 8 runs as requested
TOTAL_RUNS=8
RUNS=($(seq 1 $TOTAL_RUNS))
SEEDS=(100 101 102 103 104 105 106 107)  # run_1=100, run_2=101, etc.
=======
CONFIG="aws_test/configs/phase2_formamide_extended_SUPER_FAST.yaml"
MAX_PARALLEL=4
CHECK_INTERVAL=300  # 5 minutes

# Run configuration - 8 runs (run_1 to run_8)
TOTAL_RUNS=8
RUNS=($(seq 1 $TOTAL_RUNS))
SEEDS=(200 201 202 203 204 205 206 207)  # run_1=200, run_2=201, etc.
>>>>>>> 841e48d (koniec testów hydro, próba wycigniecia danych autokatalizy)

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
<<<<<<< HEAD
    
    nohup python3 scripts/run_phase2_full.py \
        --config "$CONFIG" \
        --output "$run_dir" \
        --steps 500000 \
        --seed "$seed" \
        --force-cpu \
        > "$log_file" 2>&1 &
    
    local pid=$!
    log "[PID] run_$run_num started with PID=$pid"
    echo "$pid" > "$BASE_DIR/logs/formamide_run_${run_num}.pid"
=======
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
>>>>>>> 841e48d (koniec testów hydro, próba wycigniecia danych autokatalizy)
}

get_next_run_to_start() {
    for run_num in "${RUNS[@]}"; do
        if ! is_run_completed "$run_num" && ! is_run_in_progress "$run_num"; then
            echo "$run_num"
            return 0
        fi
    done
    echo ""
<<<<<<< HEAD
=======
    return 1
>>>>>>> 841e48d (koniec testów hydro, próba wycigniecia danych autokatalizy)
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

<<<<<<< HEAD
count_in_progress() {
    local count=0
    for run_num in "${RUNS[@]}"; do
        if is_run_in_progress "$run_num"; then
            ((count++))
        fi
    done
    echo "$count"
}

print_status() {
    log "=========================================="
    log "FORMAMIDE QUEUE STATUS"
    log "=========================================="
    log "Time: $(date '+%Y-%m-%d %H:%M:%S')"
    
    local completed=$(count_completed)
    local in_progress=$(count_in_progress)
    local remaining=$((TOTAL_RUNS - completed - in_progress))
    
    log "Completed: $completed / $TOTAL_RUNS"
    log "In Progress: $in_progress"
    log "Remaining: $remaining"
    log "=========================================="
    
    # Show which runs are in progress
    if [ $in_progress -gt 0 ]; then
        log "Currently running:"
        for run_num in "${RUNS[@]}"; do
            if is_run_in_progress "$run_num"; then
                log "  - run_$run_num"
            fi
        done
    fi
    
    log "=========================================="
}

check_stuck_processes() {
    # Check for processes that haven't written to log in 1 hour
    for pid_file in "$BASE_DIR"/logs/formamide_run_*.pid; do
        if [ -f "$pid_file" ]; then
            pid=$(cat "$pid_file")
            run_num=$(basename "$pid_file" | sed 's/formamide_run_//' | sed 's/.pid//')
            log_file="$BASE_DIR/logs/formamide_run_${run_num}.log"
            
            if [ -f "$log_file" ]; then
                # Check if log was modified in last hour
                log_age=$(( $(date +%s) - $(stat -c %Y "$log_file" 2>/dev/null || echo 0) ))
                if [ $log_age -gt 3600 ]; then
                    log "[WARNING] run_$run_num appears stuck (no log update for ${log_age}s)"
                    log "[ACTION] Consider killing PID=$pid and restarting"
                    # Optionally auto-kill: kill -9 $pid
=======
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
>>>>>>> 841e48d (koniec testów hydro, próba wycigniecia danych autokatalizy)
                fi
            fi
        fi
    done
<<<<<<< HEAD
}

# ============================================================================
# MAIN LOOP
# ============================================================================

=======
    
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

>>>>>>> 841e48d (koniec testów hydro, próba wycigniecia danych autokatalizy)
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

<<<<<<< HEAD
# Create logs directory
mkdir -p "$BASE_DIR/logs"

=======
>>>>>>> 841e48d (koniec testów hydro, próba wycigniecia danych autokatalizy)
# Initial status
print_status

# Main loop
<<<<<<< HEAD
while true; do
    log "----------------------------------------"
    log "[CHECK] Checking queue status..."
=======
ITERATION=1
while true; do
    log "----------------------------------------"
    log "[CHECK] Iteration $ITERATION - Checking queue status..."
>>>>>>> 841e48d (koniec testów hydro, próba wycigniecia danych autokatalizy)
    
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
<<<<<<< HEAD
        log "🎉 ALL RUNS COMPLETED!"
        log "=========================================="
        log "Total completed: $completed / $TOTAL_RUNS"
        log "Results location: $BASE_DIR/$RESULTS_BASE/"
=======
        log "✅ ALL FORMAMIDE RUNS COMPLETED!"
        log "=========================================="
        log "Total runs: $TOTAL_RUNS"
        log "Completed: $completed"
>>>>>>> 841e48d (koniec testów hydro, próba wycigniecia danych autokatalizy)
        log "=========================================="
        break
    fi
    
    # Wait before next check
<<<<<<< HEAD
    log "[WAIT] Sleeping for ${CHECK_INTERVAL}s before next check..."
    sleep $CHECK_INTERVAL
done

log "=========================================="
log "Auto-restart system finished"
log "=========================================="
=======
    log "[WAIT] Waiting ${CHECK_INTERVAL}s before next check..."
    sleep "$CHECK_INTERVAL"
    ((ITERATION++))
done

log "[DONE] Formamide queue system finished"
>>>>>>> 841e48d (koniec testów hydro, próba wycigniecia danych autokatalizy)

