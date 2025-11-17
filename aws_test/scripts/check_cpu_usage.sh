#!/bin/bash
# Check CPU usage by all simulation processes

echo "=================================================================================="
echo "🔍 CPU USAGE ANALYSIS"
echo "=================================================================================="
echo ""

# Get total CPU cores
TOTAL_CORES=$(nproc)
echo "📊 System Info:"
echo "   Total CPU cores: $TOTAL_CORES"
echo ""

# Get all simulation processes
echo "📊 Active Simulation Processes:"
echo "--------------------------------------------------------------------------------"

TOTAL_CPU=0
PROCESS_COUNT=0

ps aux | grep 'run_phase2_full.py' | grep -v grep | while read line; do
    PID=$(echo "$line" | awk '{print $2}')
    CPU=$(echo "$line" | awk '{print $3}')
    
    # Extract run number
    CMDLINE=$(cat "/proc/$PID/cmdline" 2>/dev/null | tr '\0' ' ')
    RUN_NUM=$(echo "$CMDLINE" | sed -n 's/.*run_\([0-9]*\).*/\1/p')
    
    if [ -n "$RUN_NUM" ]; then
        echo "   PID $PID (run_$RUN_NUM): ${CPU}% CPU"
        PROCESS_COUNT=$((PROCESS_COUNT + 1))
        # Note: CPU is already a percentage, so we can't just add them
    fi
done

echo ""
echo "📊 CPU Usage Summary:"
echo "--------------------------------------------------------------------------------"

# Get actual CPU usage from top
echo "Current CPU usage (from top):"
top -b -n 1 | head -5

echo ""
echo "CPU usage by process (sorted):"
ps aux | grep 'run_phase2_full.py' | grep -v grep | sort -k3 -rn | head -10 | awk '{printf "   PID %s: %s%% CPU\n", $2, $3}'

echo ""
echo "💡 Analysis:"
echo "   - Each run tries to use ALL $TOTAL_CORES threads"
echo "   - With 5 parallel runs: 5 × $TOTAL_CORES = $((5 * TOTAL_CORES)) logical threads"
echo "   - This causes oversubscription and context switching overhead"
echo ""
echo "💡 Recommendation:"
echo "   - Limit each run to ~$((TOTAL_CORES / 5)) threads per run"
echo "   - Or use: $((TOTAL_CORES / 5))-13 threads per run for 5 parallel runs"
echo "   - This will reduce context switching and improve overall performance"

echo ""
echo "=================================================================================="

