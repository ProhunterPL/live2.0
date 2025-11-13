#!/bin/bash
# Smart Kill Stuck Simulations - Identifies and kills only deadlocked processes
# ============================================================================
# This script identifies processes in "Ssl" state with high CPU and old logs,
# then kills only those confirmed stuck processes.

set -e

PROJECT_ROOT="$HOME/live2.0"
RESULTS_DIR="$PROJECT_ROOT/results/phase2b_additional"
STUCK_THRESHOLD_HOURS=24  # Consider stuck if no log activity for 24+ hours

echo "================================================================================"
echo "🔍 SMART STUCK PROCESS DETECTION"
echo "================================================================================"
echo ""

# Function to check if a run is stuck based on log age
is_stuck_by_log() {
    local run_path=$1
    local log_file="$run_path/simulation.log"
    
    if [ ! -f "$log_file" ]; then
        return 1  # No log = can't determine
    fi
    
    # Get last modification time
    local last_mod=$(stat -c %Y "$log_file" 2>/dev/null || echo "0")
    local now=$(date +%s)
    local age_hours=$(( (now - last_mod) / 3600 ))
    
    # Check if log is old (24+ hours) or stuck at 160K
    if [ "$age_hours" -ge "$STUCK_THRESHOLD_HOURS" ]; then
        # Check if stuck at 160K (deadlock pattern)
        local last_step=$(tail -100 "$log_file" 2>/dev/null | grep -E "Step [0-9]+" | tail -1 | grep -oE "Step [0-9]+" | head -1)
        if echo "$last_step" | grep -qE "160000|160,000"; then
            return 0  # Stuck at 160K
        fi
        # Or if very old (48+ hours)
        if [ "$age_hours" -ge 48 ]; then
            return 0  # Stuck
        fi
    fi
    
    return 1  # Not stuck
}

echo "📊 Analyzing processes..."
echo ""

# Create temp file for stuck PIDs
STUCK_FILE=$(mktemp)
trap "rm -f $STUCK_FILE" EXIT

# Analyze each process
ps aux | grep "run_phase2_full.py" | grep -v grep | while read line; do
    PID=$(echo "$line" | awk '{print $2}')
    STATE=$(echo "$line" | awk '{print $8}')
    CPU=$(echo "$line" | awk '{print $3}')
    
    # Extract run info
    RUN_MATCH=$(echo "$line" | grep -o "run_[0-9]*" | head -1)
    SCENARIO=$(echo "$line" | grep -o "[a-z_]*_extended" | head -1)
    
    if [ -z "$RUN_MATCH" ] || [ -z "$SCENARIO" ]; then
        continue
    fi
    
    RUN_PATH="$RESULTS_DIR/$SCENARIO/$RUN_MATCH"
    
    # Check if process is in suspicious state
    if echo "$STATE" | grep -q "^S"; then
        # Sleeping state - check CPU
        CPU_VAL=$(echo "$CPU" | sed 's/%//' | cut -d'.' -f1)
        if [ "$CPU_VAL" -gt 100 ] 2>/dev/null; then
            # High CPU but sleeping = likely deadlock
            if is_stuck_by_log "$RUN_PATH"; then
                echo "🚨 STUCK: PID $PID ($SCENARIO/$RUN_MATCH) - Ssl, CPU ${CPU}%, old log"
                echo "$PID" >> "$STUCK_FILE"
            else
                echo "✅ ACTIVE: PID $PID ($SCENARIO/$RUN_MATCH) - Recent activity"
            fi
        else
            echo "✅ ACTIVE: PID $PID ($SCENARIO/$RUN_MATCH) - Low CPU sleep"
        fi
    elif echo "$STATE" | grep -q "^R"; then
        # Running state - definitely active
        echo "✅ ACTIVE: PID $PID ($SCENARIO/$RUN_MATCH) - Running"
    else
        echo "⚠️  UNKNOWN: PID $PID ($SCENARIO/$RUN_MATCH) - State: $STATE"
    fi
done

echo ""
echo "================================================================================"
echo "📊 SUMMARY"
echo "================================================================================"
echo ""

STUCK_COUNT=$(wc -l < "$STUCK_FILE" 2>/dev/null || echo "0")
TOTAL=$(ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l | tr -d ' ')

echo "Total processes: $TOTAL"
echo "Stuck processes: $STUCK_COUNT"
echo "Active processes: $((TOTAL - STUCK_COUNT))"
echo ""

if [ "$STUCK_COUNT" -eq 0 ]; then
    echo "✅ No stuck processes found. All simulations are active."
    exit 0
fi

echo "================================================================================"
echo "🚨 STUCK PROCESSES TO KILL"
echo "================================================================================"
echo ""

# List stuck processes
while read PID; do
    if [ -n "$PID" ]; then
        LINE=$(ps aux | grep "^[^ ]* *$PID " | grep -v grep)
        if [ -n "$LINE" ]; then
            RUN_MATCH=$(echo "$LINE" | grep -o "run_[0-9]*" | head -1)
            SCENARIO=$(echo "$LINE" | grep -o "[a-z_]*_extended" | head -1)
            CPU=$(echo "$LINE" | awk '{print $3}')
            STATE=$(echo "$LINE" | awk '{print $8}')
            echo "  PID $PID: $SCENARIO/$RUN_MATCH (CPU: ${CPU}%, State: $STATE)"
        fi
    fi
done < "$STUCK_FILE"

echo ""
echo "⚠️  WARNING: This will kill $STUCK_COUNT stuck processes."
read -p "Continue? (y/n) " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 1
fi

echo ""
echo "================================================================================"
echo "🚫 KILLING STUCK PROCESSES"
echo "================================================================================"
echo ""

KILLED=0
while read PID; do
    if [ -n "$PID" ]; then
        LINE=$(ps aux | grep "^[^ ]* *$PID " | grep -v grep)
        if [ -n "$LINE" ]; then
            RUN_MATCH=$(echo "$LINE" | grep -o "run_[0-9]*" | head -1)
            SCENARIO=$(echo "$LINE" | grep -o "[a-z_]*_extended" | head -1)
            echo "🚫 Killing PID $PID ($SCENARIO/$RUN_MATCH)..."
            kill -9 "$PID" 2>/dev/null || true
            KILLED=$((KILLED + 1))
        fi
    fi
done < "$STUCK_FILE"

echo ""
echo "✅ Killed $KILLED stuck processes"
echo ""

echo "================================================================================"
echo "✅ CLEANUP COMPLETE"
echo "================================================================================"
echo ""

REMAINING=$(ps aux | grep "run_phase2_full.py" | grep -v grep | wc -l | tr -d ' ')
echo "Remaining processes: $REMAINING"
echo ""
if [ "$REMAINING" -gt 0 ]; then
    echo "Active processes:"
    ps aux | grep "run_phase2_full.py" | grep -v grep | awk '{printf "  PID %s: %s (CPU: %s%%, State: %s)\n", $2, $NF, $3, $8}'
fi

echo ""
echo "================================================================================"
echo "💡 NEXT STEPS"
echo "================================================================================"
echo ""
echo "1. Verify fix is applied:"
echo "   python3 ~/live2.0/aws_test/scripts/diagnose_phase2b_issues.py"
echo ""
echo "2. Restart killed runs (maintaining ≤5 total processes):"
echo "   bash ~/live2.0/aws_test/scripts/auto_queue_restart.sh"
echo ""
echo "3. Monitor progress:"
echo "   python3 ~/live2.0/aws_test/scripts/check_actual_progress.py"
echo ""
