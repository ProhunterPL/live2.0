#!/usr/bin/env python3
"""
Phase 2B Progress Checker
==========================

Parses simulation.log files to show actual step progress and ETA.
"""

import re
import os
import sys
import time
from pathlib import Path
from datetime import datetime, timedelta
import argparse

def parse_last_step(log_file):
    """Parse the last step number from simulation.log"""
    if not log_file.exists():
        return None
    
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
            # Look for step lines in reverse order
            for line in reversed(lines):
                # Match: "Step 95000 completed in 102.7ms"
                match = re.search(r'Step (\d+) completed', line)
                if match:
                    return int(match.group(1))
    except Exception as e:
        return None
    
    return None

def parse_step_times(log_file, num_samples=5):
    """Parse recent step times to calculate average speed using real timestamps"""
    if not log_file.exists():
        return None
    
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
            step_entries = []
            
            # Look for step completion times with timestamps
            for line in reversed(lines):
                # Match: "2025-11-06 15:01:11,097 - backend.sim.core.stepper - INFO - Step 95000 completed in 108.3ms"
                match = re.search(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}),\d+.*Step (\d+) completed', line)
                if match:
                    timestamp_str = match.group(1)
                    step_num = int(match.group(2))
                    try:
                        timestamp = datetime.strptime(timestamp_str, '%Y-%m-%d %H:%M:%S')
                        step_entries.append((step_num, timestamp))
                        if len(step_entries) >= num_samples:
                            break
                    except:
                        continue
            
            if len(step_entries) < 2:
                return None
            
            # Calculate real time between steps (steps are logged every 10K steps)
            # Use timestamps to get actual elapsed time
            step_entries.sort(key=lambda x: x[0])  # Sort by step number
            
            total_steps = step_entries[-1][0] - step_entries[0][0]
            total_time = (step_entries[-1][1] - step_entries[0][1]).total_seconds()
            
            if total_time <= 0 or total_steps <= 0:
                return None
            
            steps_per_sec = total_steps / total_time
            
            return {
                'steps_per_sec': steps_per_sec,
                'total_steps': total_steps,
                'total_time_sec': total_time
            }
    except Exception as e:
        return None
    
    return None

def check_progress(results_dir="results/phase2b_additional", target_steps=500000):
    """Check progress of all Phase 2B simulations"""
    results_dir = Path(results_dir).expanduser()
    
    print("=" * 80)
    print("🔍 PHASE 2B PROGRESS CHECK")
    print("=" * 80)
    print(f"📅 Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"🎯 Target: {target_steps:,} steps per simulation")
    print()
    
    scenarios = ["miller_urey_extended", "hydrothermal_extended", "formamide_extended"]
    
    total_completed = 0
    total_running = 0
    total_failed = 0
    
    for scenario in scenarios:
        scenario_dir = results_dir / scenario
        if not scenario_dir.exists():
            print(f"📁 {scenario}: Not started")
            continue
        
        print(f"📁 {scenario.replace('_', ' ').title()}:")
        print("-" * 80)
        
        runs = sorted([d for d in scenario_dir.iterdir() if d.is_dir() and d.name.startswith('run_')],
                     key=lambda x: int(x.name.split('_')[1]) if x.name.split('_')[1].isdigit() else 0)
        
        for run_dir in runs:
            run_id = run_dir.name
            log_file = run_dir / "simulation.log"
            results_file = run_dir / "results.json"
            summary_file = run_dir / "summary.txt"
            
            # Check if completed
            if results_file.exists() and summary_file.exists():
                print(f"  ✅ {run_id}: COMPLETED")
                total_completed += 1
                continue
            
            # Check if running
            if log_file.exists():
                last_step = parse_last_step(log_file)
                if last_step is not None:
                    progress_pct = (last_step / target_steps) * 100
                    step_info = parse_step_times(log_file)
                    
                    # Check if log was updated recently (within last 10 minutes)
                    log_age_seconds = time.time() - log_file.stat().st_mtime
                    is_active = log_age_seconds < 600
                    
                    status_icon = "🔄" if is_active else "⏸️"
                    
                    print(f"  {status_icon} {run_id}: Step {last_step:,}/{target_steps:,} ({progress_pct:.1f}%)")
                    
                    if step_info:
                        steps_per_sec = step_info['steps_per_sec']
                        remaining_steps = target_steps - last_step
                        
                        if steps_per_sec > 0 and is_active:
                            eta_seconds = remaining_steps / steps_per_sec
                            eta_time = datetime.now() + timedelta(seconds=eta_seconds)
                            print(f"      Speed: {steps_per_sec:.2f} steps/sec")
                            print(f"      ETA: {eta_time.strftime('%Y-%m-%d %H:%M:%S')} ({eta_seconds/3600:.1f} hours)")
                        else:
                            print(f"      Speed: {steps_per_sec:.2f} steps/sec (may be paused)")
                    else:
                        # Fallback: estimate from log file age
                        if is_active:
                            log_age_seconds = time.time() - log_file.stat().st_mtime
                            if log_age_seconds > 0:
                                # Rough estimate: assume 8-12 steps/sec average
                                estimated_speed = 10.0  # conservative estimate
                                remaining_steps = target_steps - last_step
                                eta_seconds = remaining_steps / estimated_speed
                                eta_time = datetime.now() + timedelta(seconds=eta_seconds)
                                print(f"      Speed: ~{estimated_speed:.1f} steps/sec (estimated)")
                                print(f"      ETA: {eta_time.strftime('%Y-%m-%d %H:%M:%S')} ({eta_seconds/3600:.1f} hours)")
                    
                    if is_active:
                        total_running += 1
                    else:
                        print(f"      ⚠️  Last activity: {log_age_seconds/60:.1f} minutes ago")
                        total_failed += 1
                else:
                    print(f"  ❓ {run_id}: Log exists but no steps found")
                    total_failed += 1
            else:
                print(f"  ❌ {run_id}: No log file (not started or failed)")
                total_failed += 1
        
        print()
    
    print("=" * 80)
    print("📊 SUMMARY")
    print("=" * 80)
    print(f"✅ Completed: {total_completed}/30")
    print(f"🔄 Running: {total_running}/30")
    print(f"❌ Failed/Not Started: {total_failed}/30")
    print()
    
    # Check for results.json files
    results_json_count = len(list(results_dir.rglob("results.json")))
    print(f"📄 results.json files found: {results_json_count}")
    print()
    
    if results_json_count == 0:
        print("ℹ️  Note: results.json files are created only after simulations complete all steps.")
        print("   This is normal - simulations are still running!")
    
    print("=" * 80)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check Phase 2B simulation progress")
    parser.add_argument("--results-dir", default="results/phase2b_additional",
                       help="Results directory")
    parser.add_argument("--target-steps", type=int, default=500000,
                       help="Target number of steps per simulation")
    parser.add_argument("--watch", action="store_true",
                       help="Watch mode - update every 60 seconds")
    
    args = parser.parse_args()
    
    if args.watch:
        try:
            while True:
                os.system('clear' if os.name == 'posix' else 'cls')
                check_progress(args.results_dir, args.target_steps)
                print("\n⏱️  Refreshing in 60 seconds... (Ctrl+C to stop)")
                time.sleep(60)
        except KeyboardInterrupt:
            print("\n👋 Stopped")
    else:
        check_progress(args.results_dir, args.target_steps)

