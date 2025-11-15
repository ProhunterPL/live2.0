#!/bin/bash
# Debug what processes are actually running
# =========================================

PROJECT_ROOT="${PROJECT_ROOT:-$HOME/live2.0}"
RESULTS_DIR="$PROJECT_ROOT/results/phase2b_additional/miller_urey_extended"

echo "================================================================================"
echo "🔍 PROCESS DEBUG - $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================================"
echo ""

echo "1️⃣ RAW PROCESS LIST (all run_phase2_full.py):"
echo "--------------------------------------------------------------------------------"
ps aux | grep 'run_phase2_full.py' | grep -v grep | while read line; do
    PID=$(echo "$line" | awk '{print $2}')
    CMD=$(echo "$line" | grep -oP '\-\-output.*' | head -c 100)
    echo "PID $PID: $CMD"
done
echo ""

echo "2️⃣ MILLER_UREY PROCESSES ONLY:"
echo "--------------------------------------------------------------------------------"
ps aux | grep 'run_phase2_full.py' | grep 'miller_urey_extended' | grep -v grep | while read line; do
    PID=$(echo "$line" | awk '{print $2}')
    OUTPUT=$(echo "$line" | grep -oP 'miller_urey_extended/run_\d+' | head -1)
    if [ -n "$OUTPUT" ]; then
        RUN_ID=$(echo "$OUTPUT" | grep -oP 'run_\K\d+')
        echo "PID $PID -> run_${RUN_ID}"
    else
        echo "PID $PID -> (couldn't parse run ID)"
    fi
done
MILLER_COUNT=$(ps aux | grep 'run_phase2_full.py' | grep 'miller_urey_extended' | grep -v grep | wc -l)
echo ""
echo "Total miller_urey processes: $MILLER_COUNT"
echo ""

echo "3️⃣ COMPLETED RUNS (have results.json):"
echo "--------------------------------------------------------------------------------"
for i in {1..18}; do
    if [ -f "$RESULTS_DIR/run_${i}/results.json" ]; then
        SIZE=$(du -h "$RESULTS_DIR/run_${i}/results.json" | awk '{print $1}')
        MTIME=$(stat -c '%y' "$RESULTS_DIR/run_${i}/results.json" 2>/dev/null | cut -d'.' -f1)
        echo "✅ run_${i}: $SIZE, completed at $MTIME"
    fi
done
echo ""

echo "4️⃣ ACTIVE LOG FILES (modified in last 5 minutes):"
echo "--------------------------------------------------------------------------------"
FIVE_MIN_AGO=$(date -d '5 minutes ago' +%s)
for i in {1..18}; do
    LOG="$RESULTS_DIR/run_${i}/simulation.log"
    if [ -f "$LOG" ]; then
        MTIME=$(stat -c '%Y' "$LOG")
        if [ "$MTIME" -gt "$FIVE_MIN_AGO" ]; then
            SIZE=$(du -h "$LOG" | awk '{print $1}')
            LAST_STEP=$(grep -oP "Step \K[0-9]+" "$LOG" | tail -1)
            AGO=$(($(date +%s) - MTIME))
            echo "🔄 run_${i}: Step $LAST_STEP (${SIZE}, ${AGO}s ago)"
        fi
    fi
done
echo ""

echo "5️⃣ SUMMARY:"
echo "--------------------------------------------------------------------------------"
COMPLETED=$(ls "$RESULTS_DIR"/run_*/results.json 2>/dev/null | wc -l)
echo "Completed: $COMPLETED / 18"
echo "Running processes: $MILLER_COUNT"
echo "Expected incomplete: $((18 - COMPLETED))"
echo ""

if [ "$MILLER_COUNT" -gt 5 ]; then
    echo "⚠️  WARNING: More than 5 processes running!"
elif [ "$MILLER_COUNT" -lt 5 ] && [ "$COMPLETED" -lt 18 ]; then
    echo "💡 Tip: Can start more processes (current: $MILLER_COUNT, max: 5)"
elif [ "$COMPLETED" -eq 18 ]; then
    echo "🎉 ALL COMPLETE!"
else
    echo "✅ Looking good! ($MILLER_COUNT running, max 5)"
fi

echo ""
echo "================================================================================"

