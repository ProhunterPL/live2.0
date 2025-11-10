#!/bin/bash
# Check detailed status of simulations

echo "=================================================================================="
echo "🔍 DETAILED SIMULATION STATUS CHECK"
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
        log_file="$run_dir/simulation.log"
        
        # Check if completed
        if [ -f "$results_file" ] || [ -f "$summary_file" ]; then
            echo "  ✅ $run_id: COMPLETED"
            if [ -f "$results_file" ]; then
                echo "     📄 results.json: $(ls -lh "$results_file" | awk '{print $5, $6, $7, $8}')"
            fi
            if [ -f "$summary_file" ]; then
                echo "     📄 summary.txt: $(ls -lh "$summary_file" | awk '{print $5, $6, $7, $8}')"
            fi
        elif [ -f "$log_file" ]; then
            # Check last step
            last_step=$(grep -o "Step [0-9]* completed" "$log_file" 2>/dev/null | tail -1 | grep -o "[0-9]*")
            last_log_time=$(stat -c %y "$log_file" 2>/dev/null | cut -d' ' -f1,2 | cut -d'.' -f1)
            
            if [ -n "$last_step" ]; then
                progress=$((last_step * 100 / 500000))
                echo "  🔄 $run_id: Step $last_step/500000 ($progress%)"
                echo "     📝 Last log update: $last_log_time"
                
                # Check if log mentions completion
                if grep -q "SUCCESS\|completed\|finished" "$log_file" 2>/dev/null | tail -5; then
                    echo "     ⚠️  Log mentions completion - check for results.json"
                fi
                
                # Check last few lines
                echo "     📋 Last 3 log lines:"
                tail -3 "$log_file" 2>/dev/null | sed 's/^/        /'
            else
                echo "  🔄 $run_id: Running (no step info in log)"
                echo "     📝 Last log update: $last_log_time"
            fi
        else
            echo "  ❓ $run_id: No log file (not started or failed)"
        fi
    done
    echo ""
done

echo "=================================================================================="
echo "📊 SUMMARY"
echo "=================================================================================="

completed=$(find "$RESULTS_DIR" -name "results.json" 2>/dev/null | wc -l)
running=$(ps aux | grep run_phase2_full | grep -v grep | wc -l)

echo "✅ Completed: $completed"
echo "🔄 Running processes: $running"
echo ""

if [ "$completed" -gt 0 ]; then
    echo "📄 Completed simulations:"
    find "$RESULTS_DIR" -name "results.json" 2>/dev/null | while read f; do
        dir=$(dirname "$f")
        run_id=$(basename "$dir")
        scenario=$(basename $(dirname "$dir"))
        echo "   ✅ $scenario/$run_id"
    done
fi

echo ""

