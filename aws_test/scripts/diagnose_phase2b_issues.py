#!/usr/bin/env python3
"""
Diagnose Phase 2B Issues - Comprehensive Analysis
==================================================
Identifies stuck simulations, deadlocks, and configuration problems.
"""

import subprocess
import time
from pathlib import Path
from datetime import datetime
import json
import re

def run_command(cmd):
    """Run shell command and return output"""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=10)
        return result.stdout.strip()
    except Exception as e:
        return f"ERROR: {e}"

def check_cluster_fix_applied():
    """Check if cluster detection fix is applied in stepper.py"""
    print("="*80)
    print("🔍 CHECKING CLUSTER DETECTION FIX")
    print("="*80)
    
    stepper_path = Path.home() / "live2.0" / "backend" / "sim" / "core" / "stepper.py"
    
    if not stepper_path.exists():
        print(f"❌ stepper.py not found at {stepper_path}")
        return False
    
    content = stepper_path.read_text()
    
    # Check for the fix pattern
    if "cluster_check_interval" in content and "getattr" in content:
        print("✅ Cluster detection fix appears to be applied")
        
        # Extract the relevant section
        lines = content.split('\n')
        for i, line in enumerate(lines):
            if "cluster_check_interval" in line or "Update clusters" in line:
                start = max(0, i-2)
                end = min(len(lines), i+5)
                print("\n📄 Relevant code section:")
                for j in range(start, end):
                    marker = ">>>" if j == i else "   "
                    print(f"{marker} {j+1:4d}: {lines[j]}")
                break
        return True
    else:
        print("❌ Cluster detection fix NOT found in stepper.py")
        print("   Expected: cluster_check_interval with getattr()")
        return False

def check_config_cluster_setting():
    """Check if config files have cluster_check_interval disabled"""
    print("\n" + "="*80)
    print("🔍 CHECKING CONFIG FILES")
    print("="*80)
    
    config_dir = Path.home() / "live2.0" / "aws_test" / "configs"
    
    if not config_dir.exists():
        print(f"❌ Config directory not found: {config_dir}")
        return
    
    config_files = list(config_dir.glob("*SUPER_FAST.yaml"))
    
    if not config_files:
        print("⚠️  No SUPER_FAST config files found")
        return
    
    for config_file in config_files:
        print(f"\n📄 {config_file.name}:")
        content = config_file.read_text()
        
        # Check for cluster_check_interval
        if "cluster_check_interval" in content:
            # Extract the value
            match = re.search(r'cluster_check_interval:\s*(\d+)', content)
            if match:
                value = int(match.group(1))
                if value >= 999999999:
                    print(f"   ✅ Disabled (value: {value})")
                else:
                    print(f"   ⚠️  Enabled (value: {value}) - may cause deadlock!")
            else:
                print(f"   ❓ Found but couldn't parse value")
        else:
            print(f"   ❌ cluster_check_interval NOT SET - will use default!")

def analyze_stuck_runs():
    """Analyze which runs are stuck and why"""
    print("\n" + "="*80)
    print("🔍 ANALYZING STUCK RUNS")
    print("="*80)
    
    results_dir = Path.home() / "live2.0" / "results" / "phase2b_additional" / "miller_urey_extended"
    
    if not results_dir.exists():
        print(f"❌ Results directory not found: {results_dir}")
        return
    
    stuck_runs = []
    active_runs = []
    completed_runs = []
    
    for run_dir in sorted(results_dir.glob("run_*")):
        run_id = run_dir.name
        log_file = run_dir / "simulation.log"
        results_file = run_dir / "results.json"
        
        if results_file.exists():
            completed_runs.append(run_id)
            continue
        
        if not log_file.exists():
            continue
        
        # Get last log entry with step info
        last_lines = run_command(f"tail -50 {log_file} | grep -E 'Step [0-9]+' | tail -1")
        
        if not last_lines:
            continue
        
        # Extract step number
        match = re.search(r'Step (\d+)', last_lines)
        if not match:
            continue
        
        step = int(match.group(1))
        
        # Get file modification time
        stat = log_file.stat()
        age_hours = (time.time() - stat.st_mtime) / 3600
        
        # Check if process is running
        is_running = run_command(f"ps aux | grep 'run_phase2_full.py' | grep '{run_id}' | grep -v grep")
        
        # Determine status
        if age_hours > 24:  # No activity for 24+ hours
            stuck_runs.append({
                'run_id': run_id,
                'step': step,
                'age_hours': age_hours,
                'is_running': bool(is_running),
                'last_log': last_lines[:100]
            })
        elif is_running:
            active_runs.append({
                'run_id': run_id,
                'step': step,
                'age_hours': age_hours,
                'last_log': last_lines[:100]
            })
    
    print(f"\n✅ Completed: {len(completed_runs)} runs")
    print(f"   {', '.join(sorted(completed_runs))}")
    
    print(f"\n🏃 Active: {len(active_runs)} runs")
    for run in active_runs:
        print(f"   {run['run_id']}: Step {run['step']:,} ({run['age_hours']:.1f}h ago)")
    
    print(f"\n❌ Stuck: {len(stuck_runs)} runs")
    for run in stuck_runs:
        status = "🔄 Process running" if run['is_running'] else "💀 Process dead"
        print(f"   {run['run_id']}: Step {run['step']:,} ({run['age_hours']:.1f}h ago) - {status}")
        if run['step'] == 160000:
            print(f"      ⚠️  STUCK AT 160K - Likely cluster detection deadlock!")
    
    return stuck_runs, active_runs, completed_runs

def check_process_states():
    """Check detailed process states"""
    print("\n" + "="*80)
    print("🔍 PROCESS STATE ANALYSIS")
    print("="*80)
    
    MAX_PARALLEL = 5  # Recommended limit for stability
    
    output = run_command("ps aux | grep 'run_phase2_full.py' | grep -v grep")
    
    if not output:
        print("❌ No simulation processes found!")
        return
    
    processes = []
    for line in output.split('\n'):
        if not line:
            continue
        parts = line.split()
        if len(parts) < 11:
            continue
        
        pid = parts[1]
        cpu = parts[2]
        mem = parts[3]
        stat = parts[7]
        time_running = parts[9]
        
        # Extract run info
        run_match = [p for p in parts if 'run_' in p]
        run_name = run_match[0].split('/')[-1] if run_match else "unknown"
        
        processes.append({
            'pid': pid,
            'run': run_name,
            'cpu': cpu,
            'mem': mem,
            'stat': stat,
            'time': time_running
        })
    
    print(f"\n📊 Found {len(processes)} processes (max recommended: {MAX_PARALLEL}):")
    
    if len(processes) > MAX_PARALLEL:
        print(f"   ⚠️  WARNING: {len(processes)} processes running (exceeds limit of {MAX_PARALLEL})")
        print(f"   💡 Recommendation: Kill excess processes to maintain stability")
    elif len(processes) == MAX_PARALLEL:
        print(f"   ✅ At maximum capacity ({MAX_PARALLEL} processes)")
    else:
        print(f"   ✅ Within limit ({len(processes)}/{MAX_PARALLEL} processes)")
    
    for p in processes:
        state_desc = {
            'R': 'Running',
            'S': 'Sleeping',
            'D': 'I/O Wait',
            'Z': 'Zombie',
            'T': 'Stopped'
        }
        state = p['stat'][0] if p['stat'] else '?'
        desc = state_desc.get(state, 'Unknown')
        
        print(f"   PID {p['pid']:6s} ({p['run']:15s}): {state} ({desc}) | CPU: {p['cpu']:>6s}% | Time: {p['time']}")
        
        if state == 'S' and float(p['cpu'].replace('%', '')) > 100:
            print(f"      ⚠️  High CPU but sleeping - possible deadlock!")

def main():
    print("="*80)
    print("🔬 PHASE 2B DIAGNOSTIC ANALYSIS")
    print("="*80)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Check if fix is applied
    fix_applied = check_cluster_fix_applied()
    
    # Check config files
    check_config_cluster_setting()
    
    # Analyze stuck runs
    stuck, active, completed = analyze_stuck_runs()
    
    # Check process states
    check_process_states()
    
    # Summary and recommendations
    print("\n" + "="*80)
    print("💡 RECOMMENDATIONS")
    print("="*80)
    
    MAX_PARALLEL = 5  # Recommended limit for stability
    
    # Check current process count
    process_output = run_command("ps aux | grep 'run_phase2_full.py' | grep -v grep")
    if process_output and process_output != "ERROR":
        process_lines = [l for l in process_output.split('\n') if l.strip() and 'grep' not in l]
        process_count = len(process_lines)
    else:
        process_count = 0
    
    if process_count > MAX_PARALLEL:
        print(f"\n⚠️  PARALLEL LIMIT EXCEEDED: {process_count} processes running (max: {MAX_PARALLEL})")
        print(f"   Action: Run 'bash aws_test/scripts/limit_parallel_simulations.sh {MAX_PARALLEL}'")
        print(f"   This will kill excess processes to maintain stability")
    
    if not fix_applied:
        print("\n🚨 CRITICAL: Cluster detection fix not applied!")
        print("   Action: Apply the fix from CLUSTER_DEADLOCK_FIX.md")
    
    if stuck:
        print(f"\n🚨 Found {len(stuck)} stuck runs:")
        stuck_at_160k = [r for r in stuck if r['step'] == 160000]
        if stuck_at_160k:
            print(f"   - {len(stuck_at_160k)} runs stuck at 160K steps (cluster deadlock)")
            print("   Action: Kill these processes and restart with fixed code")
        
        running_stuck = [r for r in stuck if r['is_running']]
        if running_stuck:
            print(f"   - {len(running_stuck)} stuck processes still running (wasting CPU)")
            print("   Action: Kill these processes immediately")
            print(f"   Note: After killing, ensure total processes <= {MAX_PARALLEL}")
    
    if active:
        print(f"\n✅ {len(active)} runs are actively progressing")
        print("   Action: Monitor these - they should complete successfully")
        if len(active) < MAX_PARALLEL:
            print(f"   💡 Capacity available: {MAX_PARALLEL - len(active)} slots for new runs")
    
    if completed:
        print(f"\n✅ {len(completed)} runs completed successfully")
        print("   Action: Extract molecules from these runs")
    
    print("\n" + "="*80)

if __name__ == "__main__":
    main()

