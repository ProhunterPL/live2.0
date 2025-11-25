#!/bin/bash
# Fix Taichi Version for RTX 5070 Compatibility
# Taichi 1.7+ has memory allocation bugs with new GPUs

echo "================================"
echo "Taichi Version Fix"
echo "================================"
echo ""

echo "Problem: Taichi 1.7+ has GPU memory allocation issues"
echo "Solution: Downgrade to stable Taichi 1.6.0"
echo ""

# Check current version
echo "Current Taichi version:"
python -c "import taichi as ti; print(f'  {ti.__version__}')" 2>/dev/null || echo "  Not installed"

echo ""
echo "Uninstalling current Taichi..."
pip uninstall taichi -y

echo ""
echo "Installing Taichi 1.6.0 (stable)..."
pip install taichi==1.6.0

echo ""
echo "Verifying installation..."
python -c "import taichi as ti; print(f'Installed: Taichi {ti.__version__}')"

echo ""
echo "================================"
echo "Complete!"
echo "================================"
echo ""
echo "Now you can run:"
echo "  ./run_hybrid_test.ps1  (or python tests/benchmark_hybrid.py)"
echo ""

