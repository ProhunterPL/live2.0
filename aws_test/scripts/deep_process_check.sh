#!/bin/bash
# Deep Process Diagnostics - Check what simulations are ACTUALLY doing
# =====================================================================

echo "================================================================================"
echo "🔬 DEEP PROCESS DIAGNOSTICS"
echo "================================================================================"
echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Find all Python simulation processes
PIDS=$(ps aux | grep "run_phase2_full.py" | grep -v grep | awk '{print $2}')

if [ -z "$PIDS" ]; then
    echo "❌ No simulation processes found!"
    exit 1
fi

echo "Found processes: $PIDS"
echo ""

for PID in $PIDS; do
    echo "================================================================================"
    echo "🔍 Process PID: $PID"
    echo "================================================================================"
    
    # Get command line
    CMD=$(ps -p $PID -o args= | head -c 200)
    echo "Command: $CMD"
    echo ""
    
    # Get resource usage
    echo "📊 Resource Usage:"
    ps -p $PID -o pid,ppid,etime,%cpu,%mem,stat,cmd | tail -1
    echo ""
    
    # Check what files are open
    echo "📁 Open files (recent writes):"
    lsof -p $PID 2>/dev/null | grep -E "\.log|\.json|\.npz" | tail -10
    echo ""
    
    # Check system calls (sample for 3 seconds)
    echo "🔬 System calls (sampling for 3 seconds)..."
    timeout 3 strace -c -p $PID 2>&1 | tail -20 || echo "  (strace may require sudo)"
    echo ""
    
    # Check what the process is doing RIGHT NOW (stack trace)
    echo "🎯 Current stack trace (what is it doing?):"
    sudo gdb -batch -ex "attach $PID" -ex "thread apply all bt" -ex "detach" -ex "quit" 2>/dev/null | head -50 || echo "  (gdb may require sudo or not installed)"
    echo ""
    
    # Check CPU affinity
    echo "💻 CPU Affinity:"
    taskset -p $PID 2>/dev/null || echo "  (taskset not available)"
    echo ""
    
    # Check thread count
    echo "🧵 Thread count:"
    ps -Lp $PID | wc -l
    echo ""
    
    echo "================================================================================"
    echo ""
done

# Check disk I/O
echo "================================================================================"
echo "💾 DISK I/O ACTIVITY"
echo "================================================================================"
iostat -x 2 2 | tail -20 2>/dev/null || echo "(iostat not available)"
echo ""

# Check memory pressure
echo "================================================================================"
echo "💭 MEMORY PRESSURE"
echo "================================================================================"
free -h
echo ""
vmstat 1 3 | tail -3
echo ""

# Check if any processes are in D state (uninterruptible sleep - often I/O wait)
echo "================================================================================"
echo "⏳ PROCESSES IN UNINTERRUPTIBLE SLEEP (D state - I/O wait)"
echo "================================================================================"
D_STATE=$(ps aux | grep "run_phase2_full.py" | grep " D" | wc -l)
if [ $D_STATE -gt 0 ]; then
    echo "⚠️  WARNING: $D_STATE simulation process(es) in D state (I/O wait)"
    ps aux | grep "run_phase2_full.py" | grep " D"
else
    echo "✅ No processes in D state"
fi
echo ""

# Check load average
echo "================================================================================"
echo "📈 SYSTEM LOAD"
echo "================================================================================"
uptime
echo ""

# Summary
echo "================================================================================"
echo "💡 DIAGNOSTIC SUMMARY"
echo "================================================================================"
echo "If processes show:"
echo "  - High CPU (>100%) + R state → Actually computing ✅"
echo "  - High CPU + D state → Stuck in I/O operation ⚠️"
echo "  - Low CPU + S state → Sleeping/waiting 😴"
echo "  - Many write() syscalls → Making progress ✅"
echo "  - Many futex() syscalls → Thread synchronization (maybe slow) ⚠️"
echo ""
echo "Check system calls above to see what processes are spending time on."
echo "================================================================================"

