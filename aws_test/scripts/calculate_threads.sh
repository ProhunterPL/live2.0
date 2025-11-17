#!/bin/bash
# Calculate optimal CPU threads per run based on active runs and total cores

TOTAL_CORES=$(nproc)
ACTIVE_RUNS=$(ps aux | grep 'run_phase2_full.py' | grep -v grep | wc -l)

echo "=================================================================================="
echo "🔢 CPU THREADS CALCULATOR"
echo "=================================================================================="
echo ""
echo "📊 System Info:"
echo "   Total CPU cores: $TOTAL_CORES"
echo "   Active simulation runs: $ACTIVE_RUNS"
echo ""

if [ "$ACTIVE_RUNS" -eq 0 ]; then
    echo "⚠️  No active runs found"
    echo "   Recommendation: Use all cores (no --cpu-threads needed)"
    exit 0
fi

# Calculate threads per run
# Leave some cores for system (4-8 cores)
SYSTEM_RESERVE=4
AVAILABLE_CORES=$((TOTAL_CORES - SYSTEM_RESERVE))
THREADS_PER_RUN=$((AVAILABLE_CORES / ACTIVE_RUNS))

# Round down to nearest multiple of 2 (better for CPU)
THREADS_PER_RUN=$((THREADS_PER_RUN - (THREADS_PER_RUN % 2)))

# Minimum 4 threads per run
if [ "$THREADS_PER_RUN" -lt 4 ]; then
    THREADS_PER_RUN=4
fi

# Maximum 16 threads per run (to avoid oversubscription)
if [ "$THREADS_PER_RUN" -gt 16 ]; then
    THREADS_PER_RUN=16
fi

TOTAL_THREADS=$((ACTIVE_RUNS * THREADS_PER_RUN))

echo "💡 Recommendation:"
echo "   Threads per run: $THREADS_PER_RUN"
echo "   Total threads used: $TOTAL_THREADS / $TOTAL_CORES"
echo "   System reserve: $SYSTEM_RESERVE cores"
echo ""

echo "📝 Usage:"
echo "   Add to your command: --cpu-threads $THREADS_PER_RUN"
echo ""

echo "📊 Breakdown:"
echo "   $ACTIVE_RUNS runs × $THREADS_PER_RUN threads = $TOTAL_THREADS threads"
echo "   + $SYSTEM_RESERVE cores reserved for system = $TOTAL_CORES total"
echo ""

if [ "$TOTAL_THREADS" -gt "$TOTAL_CORES" ]; then
    echo "⚠️  WARNING: Total threads ($TOTAL_THREADS) exceeds available cores ($TOTAL_CORES)"
    echo "   Consider reducing number of parallel runs"
else
    echo "✅ Configuration looks good - no oversubscription"
fi

echo ""
echo "=================================================================================="

