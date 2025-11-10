#!/usr/bin/env python3
"""
Check ACTUAL Progress - Verify if simulations are really progressing
====================================================================

Checks:
1. File modification times (should be changing if writing snapshots)
2. Log file sizes (should be growing if logs being written)
3. Directory sizes (should grow with snapshots/checkpoints)
4. Process state (should be R = running, not D = I/O wait)
"""

import subprocess
import time
from pathlib import Path
from datetime import datetime
import json

def run_command(cmd):
    """Run shell command and return output"""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=10)
        return result.stdout.strip()
    except Exception as e:
        return f"ERROR: {e}"

def check_process_state():
    """Check if processes are actively running or stuck"""
    print("="*80)
    print("🔍 PROCESS STATE CHECK")
    print("="*80)
    
    output = run_command("ps aux | grep 'run_phase2_full.py' | grep -v grep")
    if not output:
        print("❌ No simulation processes found!")
        return
    
    for line in output.split('\n'):
        if not line:
            continue
        parts = line.split()
        if len(parts) < 11:
            continue
        
        pid = parts[1]
        cpu = parts[2]
        mem = parts[3]
        stat = parts[7]  # Process state
        time_running = parts[9]
        
        # Extract run number from command
        run_match = [p for p in parts if 'run_' in p]
        run_name = run_match[0].split('/')[-1] if run_match else "unknown"
        
        print(f"\n📊 PID {pid} ({run_name}):")
        print(f"   State: {stat} | CPU: {cpu}% | Memory: {mem}% | Running: {time_running}")
        
        # Interpret state
        if 'R' in stat:
            print(f"   ✅ Process is actively running")
        elif 'D' in stat:
            print(f"   ⚠️  Process in uninterruptible sleep (I/O wait)")
        elif 'S' in stat:
            print(f"   😴 Process sleeping (waiting for something)")
        elif 'Z' in stat:
            print(f"   💀 Process is zombie (dead but not reaped)")
        else:
            print(f"   ❓ Unknown state: {stat}")
        
        # Check thread count
        thread_count = run_command(f"ps -Lp {pid} | wc -l")
        print(f"   🧵 Threads: {thread_count}")

def check_file_activity(results_dir):
    """Check if files are being modified"""
    print("\n" + "="*80)
    print("📁 FILE ACTIVITY CHECK")
    print("="*80)
    
    results_dir = Path(results_dir).expanduser()
    
    for run_dir in sorted(results_dir.glob("*/run_*")):
        run_name = f"{run_dir.parent.name}/{run_dir.name}"
        
        # Check various directories
        log_file = run_dir / "simulation.log"
        snapshot_dir = run_dir / "snapshots"
        checkpoint_dir = run_dir / "checkpoints"
        results_file = run_dir / "results.json"
        
        print(f"\n🔍 {run_name}:")
        
        # Log file
        if log_file.exists():
            stat = log_file.stat()
            age_min = (time.time() - stat.st_mtime) / 60
            print(f"   📄 simulation.log: {stat.st_size:,} bytes, "
                  f"modified {age_min:.1f} min ago")
        else:
            print(f"   ❌ No simulation.log found")
        
        # Snapshots
        if snapshot_dir.exists():
            snapshots = list(snapshot_dir.glob("*.json")) + list(snapshot_dir.glob("*.npz"))
            if snapshots:
                latest = max(snapshots, key=lambda p: p.stat().st_mtime)
                age_min = (time.time() - latest.stat().st_mtime) / 60
                print(f"   📸 Snapshots: {len(snapshots)} files, "
                      f"latest {age_min:.1f} min ago")
            else:
                print(f"   📸 Snapshot directory exists but EMPTY")
        else:
            print(f"   ❌ No snapshots directory")
        
        # Checkpoints
        if checkpoint_dir.exists():
            checkpoints = list(checkpoint_dir.glob("*.json")) + list(checkpoint_dir.glob("*.npz"))
            if checkpoints:
                latest = max(checkpoints, key=lambda p: p.stat().st_mtime)
                age_min = (time.time() - latest.stat().st_mtime) / 60
                print(f"   💾 Checkpoints: {len(checkpoints)} files, "
                      f"latest {age_min:.1f} min ago")
            else:
                print(f"   💾 Checkpoint directory exists but EMPTY")
        else:
            print(f"   ❌ No checkpoints directory")
        
        # Results
        if results_file.exists():
            print(f"   ✅ COMPLETED (results.json exists)")
        else:
            print(f"   ⏳ In progress (no results.json yet)")

def check_directory_sizes(results_dir):
    """Check if directory sizes are growing"""
    print("\n" + "="*80)
    print("💾 DIRECTORY SIZE CHECK")
    print("="*80)
    
    results_dir = Path(results_dir).expanduser()
    
    for run_dir in sorted(results_dir.glob("*/run_*")):
        run_name = f"{run_dir.parent.name}/{run_dir.name}"
        
        # Get directory size
        size_output = run_command(f"du -sh {run_dir}")
        size = size_output.split()[0] if size_output else "unknown"
        
        print(f"   {run_name}: {size}")

def check_last_log_entry(results_dir):
    """Check what the last log entry says"""
    print("\n" + "="*80)
    print("📝 LAST LOG ENTRIES")
    print("="*80)
    
    results_dir = Path(results_dir).expanduser()
    
    for run_dir in sorted(results_dir.glob("*/run_*")):
        run_name = f"{run_dir.parent.name}/{run_dir.name}"
        log_file = run_dir / "simulation.log"
        
        if log_file.exists():
            # Get last few lines with step info
            last_lines = run_command(f"tail -100 {log_file} | grep -E 'Step [0-9]' | tail -3")
            if last_lines:
                print(f"\n🔍 {run_name}:")
                for line in last_lines.split('\n'):
                    if line:
                        print(f"   {line}")
            else:
                # Just get last line
                last_line = run_command(f"tail -1 {log_file}")
                print(f"\n🔍 {run_name}:")
                print(f"   {last_line}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Check actual simulation progress")
    parser.add_argument('--results-dir', default='~/live2.0/results/phase2b_additional',
                       help='Results directory')
    
    args = parser.parse_args()
    
    print("="*80)
    print("🔬 ACTUAL PROGRESS CHECK")
    print("="*80)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    check_process_state()
    check_file_activity(args.results_dir)
    check_directory_sizes(args.results_dir)
    check_last_log_entry(args.results_dir)
    
    print("\n" + "="*80)
    print("💡 INTERPRETATION")
    print("="*80)
    print("If you see:")
    print("  ✅ Recent file modifications → Progress is happening")
    print("  ⚠️  Old file modifications + high CPU → Log buffering (normal)")
    print("  ❌ No snapshots/checkpoints + old logs → STUCK or not progressing")
    print("  💀 Zombie processes → Simulation crashed")
    print("="*80)

if __name__ == "__main__":
    main()

