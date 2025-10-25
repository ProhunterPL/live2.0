#!/bin/bash
# AWS Phase 2B Runner Script
# Uruchamia Phase 2B na instancji AWS

set -e  # Exit on error

echo "========================================"
echo "AWS PHASE 2B RUNNER"
echo "========================================"
echo "Starting Phase 2B additional runs..."
echo ""

# Check if we're in the right directory
if [ ! -f "run_phase2b_master.py" ]; then
    echo "❌ Error: run_phase2b_master.py not found"
    echo "Please run this script from aws_test/ directory"
    exit 1
fi

# Check system resources
echo "🔍 Checking system resources..."
echo "CPU cores: $(nproc)"
echo "Memory: $(free -h | grep Mem | awk '{print $2}')"
echo "Disk space: $(df -h . | tail -1 | awk '{print $4}')"
echo ""

# Check if Python dependencies are installed
echo "🐍 Checking Python dependencies..."
python3 -c "import taichi as ti; print(f'✅ Taichi: {ti.__version__}')" || {
    echo "❌ Taichi not installed. Running setup..."
    pip3 install -r ../requirements.txt
}

echo ""

# Create necessary directories
echo "📁 Creating directories..."
mkdir -p results/phase2b_additional/logs
mkdir -p results/phase2b_additional/formamide_debug/logs
mkdir -p results/phase2b_additional/miller_urey_extended
mkdir -p results/phase2b_additional/hydrothermal_extended
mkdir -p results/phase2b_additional/formamide_extended

echo "✅ Directories created"
echo ""

# Start Phase 2B
echo "🚀 Starting Phase 2B..."
echo "This will run:"
echo "  - Debug formamide (9 tests)"
echo "  - 30 additional simulations (500K steps each)"
echo "  - Analysis and reporting"
echo ""
echo "Estimated time: 3-4 days"
echo "Estimated cost: $180-240 (c6i.16xlarge)"
echo ""

# Ask for confirmation
read -p "Continue? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "❌ Cancelled by user"
    exit 1
fi

echo ""
echo "🎯 Starting Phase 2B Master Script..."
echo "========================================"

# Run Phase 2B master script
python3 run_phase2b_master.py --mode all

echo ""
echo "========================================"
echo "🎉 PHASE 2B COMPLETED!"
echo "========================================"
echo ""
echo "📊 Results available in:"
echo "  - results/phase2b_additional/phase2b_summary_report.md"
echo "  - results/phase2b_additional/phase2b_analysis_report.md"
echo "  - results/phase2b_additional/formamide_debug/formamide_debug_report.md"
echo ""
echo "📈 Next steps:"
echo "  1. Review analysis reports"
echo "  2. Generate publication figures"
echo "  3. Proceed to Phase 3 (Paper Writing)"
echo ""
echo "💰 Total cost: ~$180-240"
echo "⏱️  Total time: ~3-4 days"
echo ""
echo "✅ Phase 2B complete!"
