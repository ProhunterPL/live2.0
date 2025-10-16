#!/bin/bash
# AWS Results Pipeline - Complete Automation
# ==========================================
# One-command pipeline: Download → Analyze → Report
#
# Usage: bash scripts/aws_pipeline.sh <aws-host> <ssh-key>

set -e

AWS_HOST="$1"
SSH_KEY="$2"

if [ -z "$AWS_HOST" ] || [ -z "$SSH_KEY" ]; then
    echo "Usage: bash scripts/aws_pipeline.sh <aws-host> <ssh-key>"
    echo ""
    echo "Example:"
    echo "  bash scripts/aws_pipeline.sh 54.123.45.67 ~/.ssh/aws_key.pem"
    exit 1
fi

echo "================================"
echo "AWS RESULTS PIPELINE"
echo "================================"
echo "Host: $AWS_HOST"
echo "Key: $SSH_KEY"
echo ""

# 1. Download results
echo "STEP 1/3: Downloading results from AWS..."
echo "-----------------------------------"
python scripts/aws_results_downloader.py \
    --host "$AWS_HOST" \
    --key "$SSH_KEY" \
    --local-base ./results/aws_batch

if [ $? -ne 0 ]; then
    echo "❌ Download failed!"
    exit 1
fi

echo ""
echo "✅ Download complete!"
echo ""

# 2. Analyze results
echo "STEP 2/3: Analyzing results..."
echo "-----------------------------------"
python scripts/aws_results_analyzer.py \
    --input ./results/aws_batch

if [ $? -ne 0 ]; then
    echo "❌ Analysis failed!"
    exit 1
fi

echo ""
echo "✅ Analysis complete!"
echo ""

# 3. Generate final report
echo "STEP 3/3: Generating publication-ready report..."
echo "-----------------------------------"

ANALYSIS_DIR="./results/aws_batch/analysis"

if [ -f "$ANALYSIS_DIR/batch_analysis_report.json" ]; then
    echo "📊 Analysis Report:"
    cat "$ANALYSIS_DIR/batch_analysis_report.json" | python -m json.tool | head -30
    echo ""
    
    echo "📁 Output files:"
    echo "  - Analysis report: $ANALYSIS_DIR/batch_analysis_report.json"
    echo "  - Figures: $ANALYSIS_DIR/figures/"
    echo "  - Comparison: $ANALYSIS_DIR/comparison/"
    echo ""
    
    # Count results
    MILLER=$(find $ANALYSIS_DIR/miller_urey -name "*_molecules.json" 2>/dev/null | wc -l)
    HYDRO=$(find $ANALYSIS_DIR/hydrothermal -name "*_molecules.json" 2>/dev/null | wc -l)
    FORM=$(find $ANALYSIS_DIR/formamide -name "*_molecules.json" 2>/dev/null | wc -l)
    TOTAL=$((MILLER + HYDRO + FORM))
    
    echo "📈 Results summary:"
    echo "  - Miller-Urey: $MILLER runs"
    echo "  - Hydrothermal: $HYDRO runs"
    echo "  - Formamide: $FORM runs"
    echo "  - TOTAL: $TOTAL runs analyzed"
    echo ""
fi

echo "================================"
echo "✅ PIPELINE COMPLETE!"
echo "================================"
echo ""
echo "Next steps:"
echo "  1. Review figures in: $ANALYSIS_DIR/figures/"
echo "  2. Check analysis report: $ANALYSIS_DIR/batch_analysis_report.json"
echo "  3. Generate paper figures: python scripts/generate_all_figures.py"
echo ""

