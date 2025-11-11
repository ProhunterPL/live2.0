#!/bin/bash
# Limit parallel simulations to max 4
# Kills excess processes if more than MAX_PARALLEL are running
# ==============================================================

MAX_PARALLEL=${1:-4}  # Default to 4, can be overridden

echo "🔍 Checking parallel simulation count..."
echo "   Max allowed: $MAX_PARALLEL"
echo ""

# Get all simulation PIDs
PIDS=$(ps aux | grep "run_phase2_full.py" | grep -v grep | awk '{print $2}')

if [ -z "$PIDS" ]; then
    echo "✅ No simulations running"
    exit 0
fi

COUNT=$(echo "$PIDS" | wc -l | tr -d ' ')
echo "📊 Currently running: $COUNT simulations"

if [ "$COUNT" -le "$MAX_PARALLEL" ]; then
    echo "✅ Within limit ($COUNT <= $MAX_PARALLEL)"
    exit 0
fi

echo "⚠️  Exceeds limit ($COUNT > $MAX_PARALLEL)"
echo "   Killing excess processes..."
echo ""

# Sort PIDs by CPU time (kill oldest/lowest CPU first)
# Get PIDs with their CPU time
PID_CPU=$(for PID in $PIDS; do
    CPU=$(ps -p $PID -o etime= 2>/dev/null | awk '{print $1}')
    echo "$PID $CPU"
done | sort -k2 -n | head -n -$MAX_PARALLEL | awk '{print $1}')

KILLED=0
for PID in $PID_CPU; do
    CMD=$(ps -p $PID -o args= 2>/dev/null)
    RUN_MATCH=$(echo "$CMD" | grep -o "run_[0-9]*" | head -1)
    SCENARIO=$(echo "$CMD" | grep -o "[a-z_]*_extended" | head -1)
    
    echo "🚫 Killing PID $PID ($SCENARIO/$RUN_MATCH)"
    kill -9 $PID 2>/dev/null
    KILLED=$((KILLED + 1))
done

echo ""
echo "✅ Killed $KILLED excess processes"
echo ""
echo "Remaining processes ($MAX_PARALLEL max):"
ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l | xargs echo "   Count:"

