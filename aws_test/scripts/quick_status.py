#!/usr/bin/env python3
"""
Quick Phase 2B Status Check - Can be run directly on AWS
"""
import re
from pathlib import Path
from datetime import datetime

results_dir = Path("results/phase2b_additional")
target_steps = 500000

print("=" * 80)
print("🔍 PHASE 2B STATUS CHECK")
print("=" * 80)
print(f"📅 Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Check each scenario
scenarios = ["miller_urey_extended", "hydrothermal_extended", "formamide_extended"]

for scenario in scenarios:
    scenario_dir = results_dir / scenario
    if not scenario_dir.exists():
        continue
    
    print(f"📁 {scenario.replace('_', ' ').title()}:")
    print("-" * 80)
    
    runs = sorted([d for d in scenario_dir.iterdir() if d.is_dir() and d.name.startswith('run_')])
    
    for run_dir in runs:
        run_id = run_dir.name
        log_file = run_dir / "simulation.log"
        results_file = run_dir / "results.json"
        
        if results_file.exists():
            print(f"  ✅ {run_id}: COMPLETED (results.json exists)")
            continue
        
        if log_file.exists():
            # Read last few lines
            try:
                with open(log_file, 'r') as f:
                    lines = f.readlines()
                    # Find last step
                    last_step = None
                    for line in reversed(lines):
                        match = re.search(r'Step (\d+) completed', line)
                        if match:
                            last_step = int(match.group(1))
                            break
                    
                    if last_step:
                        progress = (last_step / target_steps) * 100
                        # Get last line timestamp
                        last_line = lines[-1] if lines else ""
                        print(f"  🔄 {run_id}: Step {last_step:,}/{target_steps:,} ({progress:.1f}%)")
                        print(f"      Last log line: {last_line.strip()[:80]}")
                    else:
                        print(f"  ❓ {run_id}: Log exists but no steps found")
            except Exception as e:
                print(f"  ❌ {run_id}: Error reading log: {e}")
        else:
            print(f"  ❌ {run_id}: No log file")
    
    print()

# Check for results.json
results_count = len(list(results_dir.rglob("results.json")))
print(f"📄 results.json files: {results_count}")
print("=" * 80)

