#!/usr/bin/env python3
"""
Quick diagnostic script - run directly on AWS
"""
import re
import json
import time
from pathlib import Path
from datetime import datetime
import os

results_dir = Path("results/phase2b_additional")
target_steps = 500000

print("=" * 80)
print("🔍 PHASE 2B DIAGNOSTICS")
print("=" * 80)
print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Check for results.json
results_files = list(results_dir.rglob("results.json"))
print(f"📄 results.json files: {len(results_files)}")
if results_files:
    print("   Found in:")
    for f in results_files[:5]:
        print(f"   - {f.relative_to(results_dir)}")

# Check phase2b_results.json
results_json = results_dir / "phase2b_results.json"
if results_json.exists():
    try:
        with open(results_json, 'r') as f:
            data = json.load(f)
        print(f"\n📊 Status from phase2b_results.json:")
        print(f"   Completed: {data.get('completed_runs', 0)}/{data.get('total_runs', 0)}")
        print(f"   Failed: {data.get('failed_runs', 0)}/{data.get('total_runs', 0)}")
    except Exception as e:
        print(f"   Error reading JSON: {e}")

# Check simulation logs
print(f"\n📝 Simulation Logs:")
scenarios = ["miller_urey_extended", "hydrothermal_extended", "formamide_extended"]

for scenario in scenarios:
    scenario_dir = results_dir / scenario
    if not scenario_dir.exists():
        continue
    
    runs = sorted([d for d in scenario_dir.iterdir() if d.is_dir() and d.name.startswith('run_')])
    
    for run_dir in runs:
        log_file = run_dir / "simulation.log"
        results_file = run_dir / "results.json"
        
        if results_file.exists():
            print(f"  ✅ {scenario}/{run_dir.name}: COMPLETED")
            continue
        
        if log_file.exists():
            try:
                # Get file modification time
                mtime = os.path.getmtime(log_file)
                age_minutes = (time.time() - mtime) / 60
                
                # Read last step
                with open(log_file, 'r') as f:
                    lines = f.readlines()
                
                last_step = None
                for line in reversed(lines):
                    match = re.search(r'Step (\d+) completed', line)
                    if match:
                        last_step = int(match.group(1))
                        break
                
                # Check for errors
                errors = [l for l in lines[-50:] if 'ERROR' in l or 'CRITICAL' in l or 'Exception' in l]
                
                if last_step:
                    progress = (last_step / target_steps) * 100
                    status = "🔄 RUNNING" if age_minutes < 10 else "⏸️ STOPPED"
                    print(f"  {status} {scenario}/{run_dir.name}: Step {last_step:,}/{target_steps:,} ({progress:.1f}%)")
                    print(f"      Last update: {age_minutes:.1f} minutes ago")
                    if errors:
                        print(f"      ⚠️  Found {len(errors)} errors in last 50 lines")
                        print(f"      Last error: {errors[-1].strip()[:80]}")
                else:
                    print(f"  ❓ {scenario}/{run_dir.name}: Log exists but no steps found")
                    
            except Exception as e:
                print(f"  ❌ {scenario}/{run_dir.name}: Error reading log: {e}")
        else:
            print(f"  ❌ {scenario}/{run_dir.name}: No log file")

print("\n" + "=" * 80)

