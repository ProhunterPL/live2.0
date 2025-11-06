#!/usr/bin/env python3
"""
Check if simulations were restarted - compare current vs previous state
"""
import re
import json
from pathlib import Path
from datetime import datetime
import os
import time

results_dir = Path("results/phase2b_additional")

print("=" * 80)
print("🔍 CHECKING IF SIMULATIONS WERE RESTARTED")
print("=" * 80)
print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Check simulation logs for both runs
scenarios = ["miller_urey_extended"]

for scenario in scenarios:
    scenario_dir = results_dir / scenario
    if not scenario_dir.exists():
        continue
    
    runs = sorted([d for d in scenario_dir.iterdir() if d.is_dir() and d.name.startswith('run_')])
    
    for run_dir in runs:
        log_file = run_dir / "simulation.log"
        if not log_file.exists():
            continue
        
        print(f"{scenario}/{run_dir.name}:")
        print("-" * 80)
        
        try:
            # Get file size and modification time
            file_size = os.path.getsize(log_file)
            mtime = os.path.getmtime(log_file)
            age_minutes = (time.time() - mtime) / 60
            
            print(f"  Log file size: {file_size:,} bytes")
            print(f"  Last modified: {age_minutes:.1f} minutes ago")
            print(f"  Modified time: {datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')}")
            
            # Read all lines to check for multiple starts
            with open(log_file, 'r') as f:
                lines = f.readlines()
            
            print(f"  Total log lines: {len(lines):,}")
            
            # Find all step numbers
            steps = []
            for i, line in enumerate(lines):
                match = re.search(r'Step (\d+) completed', line)
                if match:
                    step_num = int(match.group(1))
                    steps.append((i, step_num))
            
            if len(steps) > 0:
                print(f"  Total step entries: {len(steps)}")
                print(f"  First step: {steps[0][1]:,} (line {steps[0][0]})")
                print(f"  Last step: {steps[-1][1]:,} (line {steps[-1][0]})")
                
                # Check if steps reset (went backwards)
                if len(steps) > 1:
                    max_step = max(s[1] for s in steps)
                    last_step = steps[-1][1]
                    
                    if max_step > last_step:
                        print(f"  ⚠️  WARNING: Steps went backwards!")
                        print(f"      Max step seen: {max_step:,}")
                        print(f"      Last step: {last_step:,}")
                        print(f"      This suggests log was overwritten or simulation restarted")
                    
                    # Check for step 0 or very low steps after high steps
                    high_steps = [s for s in steps if s[1] > 100000]
                    
                    if high_steps:
                        low_steps_after_high = [s for s in steps if s[1] < 20000 and s[0] > high_steps[-1][0]]
                        
                        if low_steps_after_high:
                            print(f"  ⚠️  WARNING: Found low steps after high steps - simulation was restarted!")
                            print(f"      Last high step: {high_steps[-1][1]:,} at line {high_steps[-1][0]}")
                            print(f"      First low step after: {low_steps_after_high[0]}")
            
            # Check first few lines for start indicators
            print(f"\n  First 5 lines:")
            for i, line in enumerate(lines[:5]):
                print(f"    {i+1}: {line.strip()[:80]}")
            
            # Check last few lines
            print(f"\n  Last 5 lines:")
            for i, line in enumerate(lines[-5:]):
                line_num = len(lines) - 5 + i + 1
                print(f"    {line_num}: {line.strip()[:80]}")
            
        except Exception as e:
            print(f"  ❌ Error reading log: {e}")
        
        print()

# Check if there are multiple log files or backups
print("=" * 80)
print("📁 Checking for backup log files:")
print("=" * 80)

for scenario in scenarios:
    scenario_dir = results_dir / scenario
    if not scenario_dir.exists():
        continue
    
    runs = sorted([d for d in scenario_dir.iterdir() if d.is_dir() and d.name.startswith('run_')])
    
    for run_dir in runs:
        log_files = list(run_dir.glob("simulation.log*"))
        if len(log_files) > 1:
            print(f"{scenario}/{run_dir.name}: Found {len(log_files)} log files")
            for log_file in log_files:
                size = os.path.getsize(log_file)
                mtime = os.path.getmtime(log_file)
                print(f"  - {log_file.name}: {size:,} bytes, modified {datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')}")

print("\n" + "=" * 80)

