#!/bin/bash
# Diagnose why simulations won't start
# ====================================

set -e

PROJECT_ROOT="$HOME/live2.0"
cd "$PROJECT_ROOT"

echo "================================================================================"
echo "🔍 DIAGNOSING STARTUP ISSUE"
echo "================================================================================"
echo ""

# 1. Check Python
echo "1️⃣ Checking Python..."
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version 2>&1)
    echo "   ✅ Python3 found: $PYTHON_VERSION"
    PYTHON_PATH=$(which python3)
    echo "   📍 Path: $PYTHON_PATH"
else
    echo "   ❌ Python3 not found!"
    exit 1
fi
echo ""

# 2. Check if script exists
echo "2️⃣ Checking script file..."
SCRIPT_PATH="$PROJECT_ROOT/scripts/run_phase2_full.py"
if [ -f "$SCRIPT_PATH" ]; then
    echo "   ✅ Script exists: $SCRIPT_PATH"
    if [ -r "$SCRIPT_PATH" ]; then
        echo "   ✅ Script is readable"
    else
        echo "   ❌ Script is NOT readable!"
    fi
else
    echo "   ❌ Script NOT found: $SCRIPT_PATH"
    exit 1
fi
echo ""

# 3. Check config file
echo "3️⃣ Checking config file..."
CONFIG_PATH="$PROJECT_ROOT/aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml"
if [ -f "$CONFIG_PATH" ]; then
    echo "   ✅ Config exists: $CONFIG_PATH"
    if [ -r "$CONFIG_PATH" ]; then
        echo "   ✅ Config is readable"
    else
        echo "   ❌ Config is NOT readable!"
    fi
else
    echo "   ❌ Config NOT found: $CONFIG_PATH"
    echo "   📋 Available configs:"
    ls -la "$PROJECT_ROOT/aws_test/configs/"*.yaml 2>/dev/null | head -5 || echo "      (none found)"
fi
echo ""

# 4. Check Python imports
echo "4️⃣ Testing Python imports..."
python3 << 'PYEOF'
import sys
print(f"   Python path: {sys.executable}")
print(f"   Python version: {sys.version}")

try:
    import taichi as ti
    print(f"   ✅ Taichi imported: {ti.__version__}")
except ImportError as e:
    print(f"   ❌ Taichi import failed: {e}")
    sys.exit(1)

try:
    import numpy as np
    print(f"   ✅ NumPy imported: {np.__version__}")
except ImportError as e:
    print(f"   ❌ NumPy import failed: {e}")
    sys.exit(1)

try:
    sys.path.insert(0, '/home/ubuntu/live2.0')
    from backend.sim.config import SimulationConfig
    print(f"   ✅ SimulationConfig imported")
except ImportError as e:
    print(f"   ❌ SimulationConfig import failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("   ✅ All imports successful")
PYEOF

if [ $? -ne 0 ]; then
    echo "   ❌ Python imports failed!"
    exit 1
fi
echo ""

# 5. Test script syntax
echo "5️⃣ Testing script syntax..."
if python3 -m py_compile "$SCRIPT_PATH" 2>&1; then
    echo "   ✅ Script syntax is valid"
else
    echo "   ❌ Script syntax error!"
    exit 1
fi
echo ""

# 6. Test dry-run (just parse arguments)
echo "6️⃣ Testing script argument parsing..."
python3 "$SCRIPT_PATH" --help > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "   ✅ Script can parse arguments"
else
    echo "   ❌ Script argument parsing failed!"
    echo "   📋 Trying with verbose output:"
    python3 "$SCRIPT_PATH" --help 2>&1 | head -20
fi
echo ""

# 7. Check output directory permissions
echo "7️⃣ Checking output directory permissions..."
OUTPUT_DIR="$PROJECT_ROOT/results/phase2b_additional/miller_urey_extended/run_3"
mkdir -p "$OUTPUT_DIR"
if [ -w "$OUTPUT_DIR" ]; then
    echo "   ✅ Output directory is writable: $OUTPUT_DIR"
else
    echo "   ❌ Output directory is NOT writable: $OUTPUT_DIR"
fi
echo ""

# 8. Test minimal run (100 steps)
echo "8️⃣ Testing minimal run (100 steps)..."
TEST_OUTPUT="$PROJECT_ROOT/results/test_startup"
mkdir -p "$TEST_OUTPUT"

echo "   🚀 Attempting to start test simulation..."
python3 "$SCRIPT_PATH" \
    --config "$CONFIG_PATH" \
    --output "$TEST_OUTPUT" \
    --seed 999 \
    --steps 100 \
    --force-cpu \
    > "$TEST_OUTPUT/startup_test.log" 2>&1 &
    
TEST_PID=$!
echo "   📍 Started with PID: $TEST_PID"

# Wait a bit to see if it starts
sleep 5

if ps -p $TEST_PID > /dev/null 2>&1; then
    echo "   ✅ Process is running!"
    echo "   📋 Checking log..."
    if [ -f "$TEST_OUTPUT/simulation.log" ]; then
        echo "   ✅ Log file created"
        tail -10 "$TEST_OUTPUT/simulation.log" 2>/dev/null || echo "      (log empty or not readable)"
    else
        echo "   ⚠️  Log file not created yet"
    fi
    
    # Kill test process
    echo "   🛑 Killing test process..."
    kill $TEST_PID 2>/dev/null || true
    sleep 2
    kill -9 $TEST_PID 2>/dev/null || true
else
    echo "   ❌ Process died immediately!"
    echo "   📋 Error output:"
    tail -30 "$TEST_OUTPUT/startup_test.log" 2>/dev/null || echo "      (no log file)"
fi
echo ""

# 9. Check for common issues
echo "9️⃣ Checking for common issues..."
echo "   📋 Current directory: $(pwd)"
echo "   📋 User: $(whoami)"
echo "   📋 Home: $HOME"
echo "   📋 Disk space:"
df -h "$PROJECT_ROOT" | tail -1
echo ""

# 10. Check if screen/nohup work
echo "🔟 Testing screen/nohup..."
if command -v screen &> /dev/null; then
    echo "   ✅ screen is available"
else
    echo "   ⚠️  screen not found (install: sudo apt-get install screen)"
fi

if command -v nohup &> /dev/null; then
    echo "   ✅ nohup is available"
else
    echo "   ❌ nohup not found!"
fi
echo ""

echo "================================================================================"
echo "✅ DIAGNOSIS COMPLETE"
echo "================================================================================"
echo ""
echo "If all checks passed, try running manually:"
echo ""
echo "  cd ~/live2.0"
echo "  python3 scripts/run_phase2_full.py \\"
echo "    --config aws_test/configs/phase2_miller_urey_extended_SUPER_FAST.yaml \\"
echo "    --output results/phase2b_additional/miller_urey_extended/run_3 \\"
echo "    --seed 103 \\"
echo "    --steps 500000 \\"
echo "    --force-cpu"
echo ""

