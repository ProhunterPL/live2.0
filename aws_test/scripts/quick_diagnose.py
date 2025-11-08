#!/usr/bin/env python3
"""
Quick diagnostic script - run directly on AWS
Enhanced version that checks CPU usage, not just log timestamps
"""
import re
import json
import time
import subprocess
from pathlib import Path
from datetime import datetime
import os

results_dir = Path("results/phase2b_additional")
target_steps = 500000

print("=" * 80)
print("🔍 PHASE 2B DIAGNOSTICS (Enhanced)")
print("=" * 80)
print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Get running processes
def get_running_processes():
    """Get CPU usage and info for running simulation processes"""
    processes = {}
    try:
        result = subprocess.run(
            ['ps', 'aux'],
            capture_output=True,
            text=True,
            check=True
        )
        
        for line in result.stdout.split('\n'):
            if 'run_phase2_full.py' in line and 'grep' not in line:
                parts = line.split()
                if len(parts) >= 11:
                    pid = parts[1]
                    cpu_percent = float(parts[2])
                    mem_percent = float(parts[3])
                    cpu_time = parts[9]
                    
                    # Extract run number from command line
                    run_match = re.search(r'run_(\d+)', line)
                    scenario_match = re.search(r'(miller_urey|hydrothermal|formamide)_extended', line)
                    
                    if run_match and scenario_match:
                        run_num = int(run_match.group(1))
                        scenario = scenario_match.group(0)
                        key = f"{scenario}/run_{run_num}"
                        
                        processes[key] = {
                            'pid': pid,
                            'cpu_percent': cpu_percent,
                            'mem_percent': mem_percent,
                            'cpu_time': cpu_time,
                            'is_active': cpu_percent > 100  # >100% means actively computing
                        }
    except Exception as e:
        print(f"⚠️  Warning: Could not check processes: {e}\n")
    
    return processes

# Get process info
running_processes = get_running_processes()

# Show process summary
if running_processes:
    print("💻 ACTIVE PROCESSES:")
    total_cpu = sum(p['cpu_percent'] for p in running_processes.values())
    active_count = sum(1 for p in running_processes.values() if p['is_active'])
    print(f"   Total processes: {len(running_processes)}")
    print(f"   Actively computing: {active_count} (CPU > 100%)")
    print(f"   Total CPU usage: {total_cpu:.0f}%")
    print()

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

# Check simulation logs and correlate with processes
print(f"\n📝 SIMULATION STATUS:")
scenarios = ["miller_urey_extended", "hydrothermal_extended", "formamide_extended"]

for scenario in scenarios:
    scenario_dir = results_dir / scenario
    if not scenario_dir.exists():
        continue
    
    runs = sorted([d for d in scenario_dir.iterdir() if d.is_dir() and d.name.startswith('run_')])
    
    for run_dir in runs:
        log_file = run_dir / "simulation.log"
        results_file = run_dir / "results.json"
        run_key = f"{scenario}/{run_dir.name}"
        
        # Check if process is running
        process_info = running_processes.get(run_key)
        
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
                    
                    # Determine status based on PROCESS, not log age
                    if process_info and process_info['is_active']:
                        status = "🔄 RUNNING"
                        status_note = f"CPU: {process_info['cpu_percent']:.0f}%, TIME: {process_info['cpu_time']}"
                    elif process_info and not process_info['is_active']:
                        status = "⏸️ PAUSED"
                        status_note = f"CPU: {process_info['cpu_percent']:.0f}% (idle)"
                    elif age_minutes < 10:
                        status = "🔄 RUNNING"
                        status_note = "Recent log update"
                    else:
                        status = "⏸️ STOPPED"
                        status_note = f"No process found, log {age_minutes:.0f}min old"
                    
                    print(f"  {status} {scenario}/{run_dir.name}: Step {last_step:,}/{target_steps:,} ({progress:.1f}%)")
                    print(f"      Status: {status_note}")
                    
                    # Add note about log buffering
                    if process_info and process_info['is_active'] and age_minutes > 10:
                        print(f"      ℹ️  Log buffered ({age_minutes:.0f}min old) but process actively computing")
                    
                    if errors:
                        print(f"      ⚠️  Found {len(errors)} errors in last 50 lines")
                        print(f"      Last error: {errors[-1].strip()[:80]}")
                else:
                    if process_info:
                        print(f"  🔄 STARTING {scenario}/{run_dir.name}: Initializing (CPU: {process_info['cpu_percent']:.0f}%)")
                    else:
                        print(f"  ❓ {scenario}/{run_dir.name}: Log exists but no steps found")
                    
            except Exception as e:
                print(f"  ❌ {scenario}/{run_dir.name}: Error reading log: {e}")
        else:
            if process_info:
                print(f"  🔄 STARTING {scenario}/{run_dir.name}: Process running (CPU: {process_info['cpu_percent']:.0f}%) but no log yet")
            else:
                print(f"  ⏳ WAITING {scenario}/{run_dir.name}: Not started yet")

print("\n" + "=" * 80)
print("ℹ️  NOTE: Status based on CPU usage (>100% = actively computing)")
print("   Log timestamps may be delayed due to buffering - this is normal!")
print("=" * 80)

