#!/bin/bash
# Direct Progress Checker for Phase 2B
# Works by reading simulation logs directly

echo "================================================================================"
echo "🔍 PHASE 2B PROGRESS - DIRECT CHECK"
echo "================================================================================"
echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Function to get last step from log
get_last_step() {
    local logfile="$1"
    grep "Step [0-9]* completed" "$logfile" 2>/dev/null | tail -1 | grep -oP 'Step \K[0-9]+' || echo "0"
}

# Function to get last log time
get_last_time() {
    local logfile="$1"
    tail -1 "$logfile" 2>/dev/null | grep -oP '^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}'
}

# Function to calculate age in minutes
calculate_age() {
    local log_time="$1"
    if [ -n "$log_time" ]; then
        local last_epoch=$(date -d "$log_time" +%s 2>/dev/null || echo 0)
        local now_epoch=$(date +%s)
        echo $(( ($now_epoch - $last_epoch) / 60 ))
    else
        echo "?"
    fi
}

# Base directory
BASE_DIR="$HOME/live2.0/results/phase2b_additional"

# Check each scenario
for scenario in miller_urey_extended hydrothermal_extended formamide_extended; do
    echo "📊 $scenario"
    echo "----------------------------------------"
    
    scenario_dir="$BASE_DIR/$scenario"
    
    if [ ! -d "$scenario_dir" ]; then
        echo "  ⏳ Scenario not started yet"
        echo ""
        continue
    fi
    
    for run_id in {1..10}; do
        logfile="$scenario_dir/run_$run_id/simulation.log"
        
        if [ -f "$logfile" ]; then
            last_step=$(get_last_step "$logfile")
            last_time=$(get_last_time "$logfile")
            age_min=$(calculate_age "$last_time")
            
            if [ "$last_step" -gt 0 ]; then
                pct=$(echo "scale=1; $last_step / 5000" | bc 2>/dev/null || echo "?")
                
                # Determine status based on age
                if [ "$age_min" != "?" ] && [ "$age_min" -lt 5 ]; then
                    status="🔄 RUNNING"
                elif [ "$age_min" != "?" ] && [ "$age_min" -lt 60 ]; then
                    status="⏸️ PAUSED "
                else
                    status="⏹️ STOPPED"
                fi
                
                printf "  %-10s run_%d: Step %7s / 500,000 (%5s%%) - %4s min ago\n" \
                    "$status" "$run_id" "$last_step" "$pct" "$age_min"
            else
                echo "  🔄 STARTING run_$run_id: Initializing..."
            fi
        else
            echo "  ⏳ WAITING  run_$run_id: Not started yet"
        fi
    done
    
    echo ""
done

# Show process count
echo "================================================================================"
echo "💻 SYSTEM STATUS"
echo "================================================================================"
process_count=$(ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l)
echo "Active processes: $process_count"

if [ "$process_count" -gt 0 ]; then
    echo ""
    echo "Running simulations:"
    ps aux | grep "run_phase2_full.py" | grep -v grep | awk '{print "  PID " $2 ": CPU " $3 "%, MEM " $4 "%, TIME " $10}' | head -10
fi

echo ""
echo "================================================================================"

