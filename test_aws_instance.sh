#!/bin/bash
# Test AWS Instance Performance
# Usage: bash test_aws_instance.sh

set -e

echo "================================"
echo "Testing AWS Instance Performance"
echo "================================"

cd live2.0

# Run quick test
echo "Running test simulation (1000 steps)..."
python3 scripts/run_phase2_full.py \
    --config configs/phase2_quick_test.yaml \
    --output results/aws_test \
    --steps 1000 \
    --seed 42

echo ""
echo "================================"
echo "✅ Test complete!"
echo "================================"
echo ""
echo "Check results:"
echo "  cat results/aws_test/summary.txt"
echo ""
echo "Expected performance: 4-6 steps/s"
echo "If performance is good, proceed with production runs."
echo ""

