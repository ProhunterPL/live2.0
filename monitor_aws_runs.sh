#!/bin/bash
# Monitor AWS Runs
# Usage: bash monitor_aws_runs.sh

echo "================================"
echo "AWS Simulation Monitoring"
echo "================================"
echo ""

cd live2.0/results

# Count completed runs
COMPLETED=$(find . -name "summary.txt" -type f 2>/dev/null | wc -l)
echo "✅ Completed simulations: $COMPLETED"

# Show recent activity
echo ""
echo "Recent log activity:"
echo "-------------------"
if [ -f "phase2_aws/master.log" ]; then
    tail -20 phase2_aws/master.log
else
    echo "No master log found yet."
fi

echo ""
echo "================================"
echo "System Resources:"
echo "================================"
echo ""

# CPU usage
echo "CPU Usage:"
top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print "  Used: " 100 - $1"%"}'

# Memory usage
echo ""
echo "Memory Usage:"
free -h | grep Mem | awk '{print "  Used: "$3" / "$2" ("$3/$2*100"%)"}'

# Disk usage
echo ""
echo "Disk Usage:"
df -h . | tail -1 | awk '{print "  Used: "$3" / "$2" ("$5")"}'

echo ""
echo "For live monitoring, run: htop"
echo "For detailed logs, run: tail -f results/phase2_aws/master.log"
echo ""

