#!/bin/bash
# Quick Start - Formamide Queue on AWS
# ====================================
# 
# Simple launcher for formamide auto-restart system
# Run this after hydrothermal runs complete

set -e

cd ~/live2.0

echo "=========================================="
echo "🧪 FORMAMIDE QUEUE - AWS LAUNCHER"
echo "=========================================="
echo ""

# Check if already running
if pgrep -f "auto_queue_restart_formamide.sh" > /dev/null; then
    echo "❌ ERROR: Formamide queue already running!"
    echo ""
    echo "Current processes:"
    ps aux | grep "auto_queue_restart_formamide.sh" | grep -v grep
    echo ""
    echo "To stop: pkill -f auto_queue_restart_formamide.sh"
    exit 1
fi

# Check if hydrothermal still running
hydro_running=$(ps aux | grep "run_phase2_full.py" | grep "hydrothermal" | grep -v grep | wc -l)
if [ $hydro_running -gt 0 ]; then
    echo "⚠️  WARNING: Hydrothermal simulations still running ($hydro_running processes)"
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
echo "   Config: aws_test/configs/phase2_formamide_AWS_OPTIMIZED.yaml"
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

sleep 2

# Check if started
if pgrep -f "auto_queue_restart_formamide.sh" > /dev/null; then
    echo ""
    echo "✅ SUCCESS! Formamide queue started"
    echo ""
    echo "=========================================="
    echo "📋 MONITORING COMMANDS"
    echo "=========================================="
    echo ""
    echo "1. Check progress:"
    echo "   python3 ~/live2.0/aws_test/scripts/check_actual_progress.py"
    echo ""
    echo "2. View main log:"
    echo "   tail -f logs/auto_restart_formamide_main.log"
    echo ""
    echo "3. View progress log:"
    echo "   tail -f logs/auto_restart_formamide_progress.log"
    echo ""
    echo "4. View specific run log (e.g., run_1):"
    echo "   tail -f logs/formamide_run_1.log"
    echo ""
    echo "5. Check running processes:"
    echo "   ps aux | grep run_phase2_full | grep formamide"
    echo ""
    echo "6. Stop the queue:"
    echo "   pkill -f auto_queue_restart_formamide.sh"
    echo "   # Then kill individual simulations:"
    echo "   pkill -f 'run_phase2_full.py.*formamide'"
    echo ""
    echo "=========================================="
    echo ""
    echo "🎯 Queue will manage 8 runs automatically"
    echo "   Expected completion: ~12-16 hours from now"
    echo ""
    echo "Good luck! 🚀"
    echo ""
else
    echo ""
    echo "❌ ERROR: Failed to start queue"
    echo "Check logs/auto_restart_formamide_launcher.log for details"
    exit 1
fi

