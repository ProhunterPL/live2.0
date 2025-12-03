#!/bin/bash
# Quick Start - Formamide Queue on AWS
# ========================================
# 
# Simple launcher for formamide auto-restart system
# Run this to start 8 formamide runs (4 parallel max)

set -e

cd ~/live2.0

echo "=========================================="
echo "🧪 FORMAMIDE QUEUE - AWS LAUNCHER"
echo "=========================================="
echo ""

# Check if already running
if pgrep -f "auto_queue_restart_formamide.sh" > /dev/null; then
    echo "[ERROR] Formamide queue already running!"
    echo ""
    echo "Current processes:"
    ps aux | grep "auto_queue_restart_formamide.sh" | grep -v grep
    echo ""
    echo "To stop: pkill -f auto_queue_restart_formamide.sh"
    exit 1
fi

# Check if other simulations still running
other_running=$(ps aux | grep "run_phase2_full.py" | grep -v "formamide" | grep -v grep | wc -l)
if [ $other_running -gt 0 ]; then
    echo "[WARNING] Other simulations still running ($other_running processes)"
    echo ""
    read -p "Continue anyway? (y/N): " confirm
    if [ "$confirm" != "y" ] && [ "$confirm" != "Y" ]; then
        echo "Cancelled."
        exit 0
    fi
fi

# System info
echo ""
echo "📊 System Info:"
echo "   CPU Cores: $(nproc)"
echo "   Memory: $(free -h | grep Mem | awk '{print $2}')"
echo "   Disk Free: $(df -h ~ | tail -1 | awk '{print $4}')"
echo ""

# Config info
echo "⚙️  Configuration:"
echo "   Config: aws_test/configs/phase2_formamide_extended_SUPER_FAST.yaml"
echo "   Total runs: 8 (run_1 to run_8)"
echo "   Max parallel: 4"
echo "   Steps per run: 500,000"
echo "   Expected time per run: 6-8 hours"
echo "   Total expected time: ~12-16 hours"
echo ""

# Confirm
echo "=========================================="
read -p "🚀 Start formamide queue? (y/N): " confirm
echo ""

if [ "$confirm" != "y" ] && [ "$confirm" != "Y" ]; then
    echo "Cancelled."
    exit 0
fi

# Create logs directory
mkdir -p logs

# Start the queue
echo "Starting auto-restart system..."
nohup bash aws_test/scripts/auto_queue_restart_formamide.sh > logs/auto_restart_formamide_launcher.log 2>&1 &

echo ""
echo "=========================================="
echo "✅ FORMAMIDE QUEUE STARTED!"
echo "=========================================="
echo ""
echo "Monitor with:"
echo "  tail -f logs/auto_restart_formamide_main.log"
echo ""
echo "Check status:"
echo "  ps aux | grep formamide | grep run_phase2_full"
echo ""
echo "Stop queue:"
echo "  pkill -f auto_queue_restart_formamide.sh"
echo ""

