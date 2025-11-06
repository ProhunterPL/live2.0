#!/usr/bin/env python3
"""
Check why simulations stopped - detailed diagnostics
"""
import re
import json
from pathlib import Path
from datetime import datetime
import os
import time

results_dir = Path("results/phase2b_additional")

print("=" * 80)
print("🔍 DETAILED DIAGNOSTICS - Why Simulations Stopped")
print("=" * 80)
print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Check phase2b_results.json for failed runs
results_json = results_dir / "phase2b_results.json"
if results_json.exists():
    try:
        with open(results_json, 'r') as f:
            data = json.load(f)
        
        print("📊 Failed Runs Details:")
        print("-" * 80)
        for scenario_name, scenario_data in data.get('scenarios', {}).items():
            for run in scenario_data.get('runs', []):
                if run.get('status') != 'success':
                    print(f"\n{scenario_name} run_{run.get('run_id')}:")
                    print(f"  Status: {run.get('status', 'unknown')}")
                    print(f"  Duration: {run.get('duration', 0):.1f}s")
                    if 'error' in run:
                        error = run['error']
                        if len(error) > 200:
                            error = error[:200] + "..."
                        print(f"  Error: {error}")
    except Exception as e:
        print(f"Error reading JSON: {e}")

# Check simulation logs for errors
print("\n" + "=" * 80)
print("📝 Checking Logs for Errors:")
print("=" * 80)

scenarios = ["miller_urey_extended", "hydrothermal_extended", "formamide_extended"]

for scenario in scenarios:
    scenario_dir = results_dir / scenario
    if not scenario_dir.exists():
        continue
    
    runs = sorted([d for d in scenario_dir.iterdir() if d.is_dir() and d.name.startswith('run_')])
    
    for run_dir in runs:
        log_file = run_dir / "simulation.log"
        if not log_file.exists():
            continue
        
        print(f"\n{scenario}/{run_dir.name}:")
        print("-" * 80)
        
        try:
            with open(log_file, 'r') as f:
                lines = f.readlines()
            
            # Get last 50 lines
            last_lines = lines[-50:]
            
            # Check for errors
            errors = []
            for i, line in enumerate(last_lines):
                if any(keyword in line.upper() for keyword in ['ERROR', 'CRITICAL', 'EXCEPTION', 'TRACEBACK', 'FAILED', 'CRASHED', 'KILLED', 'OOM']):
                    errors.append((len(lines) - 50 + i, line.strip()))
            
            if errors:
                print(f"  ⚠️  Found {len(errors)} error(s) in last 50 lines:")
                for line_num, error_line in errors[-5:]:  # Show last 5 errors
                    print(f"    Line {line_num}: {error_line[:100]}")
            else:
                print("  ✅ No obvious errors in last 50 lines")
            
            # Show last few lines
            print(f"\n  Last 5 lines of log:")
            for line in last_lines[-5:]:
                print(f"    {line.strip()}")
            
            # Check if log ends abruptly
            last_line = lines[-1] if lines else ""
            if 'completed' not in last_line.lower() and 'error' not in last_line.lower():
                print(f"\n  ⚠️  Log ends abruptly - process may have been killed")
            
        except Exception as e:
            print(f"  ❌ Error reading log: {e}")

# Check master runner log
print("\n" + "=" * 80)
print("📋 Master Runner Log (last 30 lines):")
print("=" * 80)

master_log = results_dir / "logs" / "phase2b_runner.log"
if master_log.exists():
    try:
        with open(master_log, 'r') as f:
            lines = f.readlines()
        for line in lines[-30:]:
            print(line.rstrip())
    except Exception as e:
        print(f"Error reading master log: {e}")
else:
    print("Master log not found")

print("\n" + "=" * 80)

