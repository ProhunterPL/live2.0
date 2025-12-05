#!/bin/bash
# Quick script to build reaction networks on AWS
# Usage: bash scripts/aws_build_reaction_networks.sh

set -e

echo "=========================================="
echo "BUILDING REACTION NETWORKS ON AWS"
echo "=========================================="

cd ~/live2.0

# Pull latest code
echo "[1/3] Pulling latest code..."
git pull origin main || echo "Warning: git pull failed (continuing...)"

# Run batch script with parallel processing
echo "[2/3] Building reaction networks (4 parallel workers)..."
python3 scripts/build_reaction_networks_batch.py \
    --scenario hydrothermal_extended \
    --base-dir results/phase2b_additional \
    --parallel 4

# Verify results
echo "[3/3] Verifying results..."
SUCCESS_COUNT=$(find results/phase2b_additional/hydrothermal_extended/run_*/reaction_network.json 2>/dev/null | wc -l)
echo "Generated reaction networks: $SUCCESS_COUNT/17"

if [ "$SUCCESS_COUNT" -eq 17 ]; then
    echo ""
    echo "=========================================="
    echo "[OK] ALL REACTION NETWORKS GENERATED!"
    echo "=========================================="
    echo ""
    echo "Next steps:"
    echo "1. Run autocatalysis detection:"
    echo "   python3 scripts/analyze_phase2b_complete.py \\"
    echo "       --input results/phase2b_additional \\"
    echo "       --output paper/results_data"
    echo ""
    echo "2. Or download results locally (optional):"
    echo "   scp -r -i <your-key.pem> ubuntu@<AWS_IP>:~/live2.0/results/phase2b_additional/hydrothermal_extended/run_*/reaction_network.json \\"
    echo "       results/phase2b_additional/hydrothermal_extended/"
    exit 0
else
    echo ""
    echo "=========================================="
    echo "[WARNING] Some runs failed ($SUCCESS_COUNT/17)"
    echo "=========================================="
    exit 1
fi

