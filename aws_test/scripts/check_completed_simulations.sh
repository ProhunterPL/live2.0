#!/bin/bash
# Quick check for completed simulations

echo "=================================================================================="
echo "🔍 CHECKING FOR COMPLETED SIMULATIONS"
echo "=================================================================================="
echo ""

RESULTS_DIR="${1:-results/phase2b_additional}"

for scenario in miller_urey_extended hydrothermal_extended formamide_extended; do
    scenario_dir="$RESULTS_DIR/$scenario"
    if [ ! -d "$scenario_dir" ]; then
        continue
    fi
    
    echo "📁 $scenario:"
    echo "--------------------------------------------------------------------------------"
    
    for run_dir in "$scenario_dir"/run_*; do
        if [ ! -d "$run_dir" ]; then
            continue
        fi
        
        run_id=$(basename "$run_dir")
        results_file="$run_dir/results.json"
        summary_file="$run_dir/summary.txt"
        
        if [ -f "$results_file" ] || [ -f "$summary_file" ]; then
            echo "  ✅ $run_id: COMPLETED"
            if [ -f "$results_file" ]; then
                echo "     📄 results.json exists"
            fi
            if [ -f "$summary_file" ]; then
                echo "     📄 summary.txt exists"
            fi
        else
            # Check if process is still running
            log_file="$run_dir/simulation.log"
            if [ -f "$log_file" ]; then
                last_step=$(grep -o "Step [0-9]* completed" "$log_file" | tail -1 | grep -o "[0-9]*")
                if [ -n "$last_step" ]; then
                    progress=$((last_step * 100 / 500000))
                    echo "  🔄 $run_id: Step $last_step/500000 ($progress%)"
                else
                    echo "  🔄 $run_id: Running (no step info)"
                fi
            else
                echo "  ❓ $run_id: No log file"
            fi
        fi
    done
    echo ""
done

echo "=================================================================================="
echo "📊 SUMMARY"
echo "=================================================================================="

completed=$(find "$RESULTS_DIR" -name "results.json" -o -name "summary.txt" | wc -l)
echo "✅ Completed simulations: $completed"
echo ""

