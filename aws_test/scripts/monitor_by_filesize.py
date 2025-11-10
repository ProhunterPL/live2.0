#!/usr/bin/env python3
"""
Monitor Phase2B Progress by File Size Changes
==============================================

Since logs are buffered, monitor checkpoint/snapshot file creation
and simulation.log file size changes to detect progress.
"""

import subprocess
import time
from pathlib import Path
from datetime import datetime
import json

def get_file_info(filepath):
    """Get file size and modification time"""
    try:
        stat = filepath.stat()
        return {
            'size': stat.st_size,
            'mtime': stat.st_mtime
        }
    except FileNotFoundError:
        return None

def monitor_progress(results_dir, interval=300):
    """
    Monitor progress by tracking file changes
    
    Args:
        results_dir: Path to results directory
        interval: Check interval in seconds (default: 5 minutes)
    """
    results_dir = Path(results_dir)
    cache_file = Path.home() / '.phase2b_filesize_cache.json'
    
    # Load previous state
    if cache_file.exists():
        with open(cache_file) as f:
            prev_state = json.load(f)
    else:
        prev_state = {}
    
    current_state = {}
    changes_detected = []
    
    print("=" * 80)
    print(f"🔍 MONITORING FILE-BASED PROGRESS")
    print("=" * 80)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Check all simulation runs
    for run_dir in sorted(results_dir.glob("*/run_*")):
        run_name = f"{run_dir.parent.name}/{run_dir.name}"
        
        # Check log file size
        log_file = run_dir / "simulation.log"
        log_info = get_file_info(log_file)
        
        # Check for checkpoint files
        checkpoints = list(run_dir.glob("checkpoint_*.npz"))
        snapshots = list(run_dir.glob("snapshot_*.npz"))
        
        # Check for results
        results_file = run_dir / "results.json"
        has_results = results_file.exists()
        
        current_state[run_name] = {
            'log_size': log_info['size'] if log_info else 0,
            'log_mtime': log_info['mtime'] if log_info else 0,
            'checkpoints': len(checkpoints),
            'snapshots': len(snapshots),
            'completed': has_results
        }
        
        # Compare with previous state
        if run_name in prev_state:
            prev = prev_state[run_name]
            curr = current_state[run_name]
            
            # Detect changes
            log_grew = curr['log_size'] > prev['log_size']
            new_checkpoints = curr['checkpoints'] > prev['checkpoints']
            new_snapshots = curr['snapshots'] > prev['snapshots']
            newly_completed = curr['completed'] and not prev['completed']
            
            if newly_completed:
                print(f"✅ {run_name}: COMPLETED!")
                changes_detected.append(run_name)
            elif log_grew or new_checkpoints or new_snapshots:
                print(f"🔄 {run_name}: ACTIVE")
                print(f"   Log: {prev['log_size']:,} → {curr['log_size']:,} bytes "
                      f"({curr['log_size'] - prev['log_size']:+,})")
                if new_checkpoints:
                    print(f"   Checkpoints: {prev['checkpoints']} → {curr['checkpoints']} "
                          f"(+{curr['checkpoints'] - prev['checkpoints']})")
                if new_snapshots:
                    print(f"   Snapshots: {prev['snapshots']} → {curr['snapshots']} "
                          f"(+{curr['snapshots'] - prev['snapshots']})")
                changes_detected.append(run_name)
            elif curr['completed']:
                print(f"✅ {run_name}: Already completed")
            else:
                print(f"⏸️  {run_name}: NO CHANGES")
                print(f"   Log: {curr['log_size']:,} bytes (unchanged)")
                print(f"   Checkpoints: {curr['checkpoints']}, Snapshots: {curr['snapshots']}")
        else:
            # New simulation
            status = "COMPLETED" if current_state[run_name]['completed'] else "NEW"
            print(f"🆕 {run_name}: {status}")
            print(f"   Log: {current_state[run_name]['log_size']:,} bytes")
            print(f"   Checkpoints: {current_state[run_name]['checkpoints']}, "
                  f"Snapshots: {current_state[run_name]['snapshots']}")
        
        print()
    
    # Save current state
    with open(cache_file, 'w') as f:
        json.dump(current_state, f, indent=2)
    
    print("=" * 80)
    if changes_detected:
        print(f"✅ {len(changes_detected)} simulations showing progress!")
    else:
        print("⚠️  No changes detected since last check")
        print("   This could mean:")
        print("   1. Simulations just started (need more time)")
        print("   2. Log buffering (normal - progress still happening)")
        print("   3. All simulations completed")
    print("=" * 80)
    print(f"\n💡 Run this script again in {interval//60} minutes to check for changes")
    print(f"   Previous state saved to: {cache_file}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Monitor Phase2B progress by file changes")
    parser.add_argument('--results-dir', default='~/live2.0/results/phase2b_additional',
                       help='Results directory (default: ~/live2.0/results/phase2b_additional)')
    parser.add_argument('--interval', type=int, default=300,
                       help='Suggested check interval in seconds (default: 300 = 5 min)')
    
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir).expanduser()
    monitor_progress(results_dir, args.interval)

