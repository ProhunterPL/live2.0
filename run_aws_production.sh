#!/bin/bash
# Run Production Simulations on AWS
# Usage: bash run_aws_production.sh [parallel_jobs]

set -e

# Detect CPU count
NCPUS=$(nproc)
DEFAULT_PARALLEL=$((NCPUS / 4))  # 4 cores per job
PARALLEL=${1:-$DEFAULT_PARALLEL}

echo "================================"
echo "AWS Production Run - Phase 2"
echo "================================"
echo "CPUs detected: $NCPUS"
echo "Parallel jobs: $PARALLEL"
echo "================================"
echo ""

cd live2.0

# Create session log
SESSION_DIR="results/aws_session_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$SESSION_DIR"

echo "Session directory: $SESSION_DIR"
echo ""

# Option 1: Use master orchestrator (RECOMMENDED)
echo "Starting master orchestrator..."
python3 scripts/phase2_master_1M.py \
    --mode full \
    --scenarios all \
    --max-parallel $PARALLEL \
    --output-dir results/phase2_aws \
    2>&1 | tee "$SESSION_DIR/master.log"

echo ""
echo "================================"
echo "✅ Production run complete!"
echo "================================"
echo ""
echo "Results saved to: results/phase2_aws/"
echo ""
echo "To download results to local machine:"
echo "  scp -r -i your-key.pem ubuntu@<IP>:~/live2.0/results/phase2_aws ./results/"
echo ""

