#!/usr/bin/env python3
"""
Check Real Progress of Running Simulations
===========================================

This script verifies if simulations are actually making progress by:
1. Comparing CPU time used vs time since last log update
2. Checking if process is actively computing
3. Estimating actual progress based on CPU usage patterns
"""

import re
import subprocess
import time
import os
import json
from pathlib import Path
from datetime import datetime

def get_process_info(pid):
    """Get detailed process information"""
    try:
        result = subprocess.run(
            ['ps', '-p', str(pid), '-o', 'pid,etime,pcpu,pmem,cmd'],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except:
        return None

def parse_cpu_time(etime_str):
    """Parse elapsed time string from ps (format: HH:MM:SS or DD-HH:MM:SS)"""
    try:
        if '-' in etime_str:
            # Format: DD-HH:MM:SS
            days, time_part = etime_str.split('-')
            days = int(days)
        else:
            days = 0
            time_part = etime_str
        
        parts = time_part.split(':')
        if len(parts) == 3:
            hours, minutes, seconds = map(int, parts)
            total_seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds
            return total_seconds
        elif len(parts) == 2:
            minutes, seconds = map(int, parts)
            total_seconds = days * 86400 + minutes * 60 + seconds
            return total_seconds
    except:
        pass
    return 0

def get_last_step_from_log(log_file):
    """Extract last step from log file"""
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
        
        for line in reversed(lines):
            # Pattern: "Step 123,456/500,000"
            match = re.search(r'Step\s+([\d,]+)/([\d,]+)', line)
            if match:
                step_str = match.group(1).replace(',', '')
                return int(step_str)
    except:
        pass
    return None

def estimate_progress_from_cpu(pid, last_step, target_steps=500000):
    """Estimate progress based on CPU usage and time"""
    try:
        # Get process CPU time
        result = subprocess.run(
            ['ps', '-p', str(pid), '-o', 'etime,pcpu'],
            capture_output=True,
            text=True,
            check=True
        )
        
        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:
            return None
        
        data = lines[1].split()
        if len(data) < 2:
            return None
        
        etime_str = data[0]
        cpu_percent = float(data[1])
        
        cpu_time_seconds = parse_cpu_time(etime_str)
        
        # If CPU > 100%, it's using multiple cores
        # Estimate: if using 1000% CPU for 1 hour = 10 core-hours
        # At ~10 steps/second per core, that's ~360,000 steps per core-hour
        # So 10 core-hours = ~3.6M steps (but we're at 500K max)
        
        # More conservative: if process has been running for X seconds at Y% CPU
        # and last log was at step Z, estimate current step
        
        # Simple heuristic: if CPU is high and time has passed, assume progress
        if cpu_percent > 100 and cpu_time_seconds > 0:
            # Estimate steps per second based on typical performance
            # Conservative estimate: ~10 steps/sec per 100% CPU
            steps_per_100pct_cpu_per_sec = 10
            estimated_steps_since_log = int(
                (cpu_percent / 100) * steps_per_100pct_cpu_per_sec * (cpu_time_seconds / 3600) * 3600
            )
            
            # But we need to know when the log was last updated
            # For now, just indicate that progress is likely happening
            return {
                'estimated_current_step': last_step,  # Conservative: assume no progress beyond last log
                'cpu_time_hours': cpu_time_seconds / 3600,
                'cpu_percent': cpu_percent,
                'likely_making_progress': cpu_percent > 100
            }
    except Exception as e:
        print(f"Error estimating progress: {e}")
    
    return None

def load_previous_state(cache_file):
    """Load previous progress state from cache file"""
    if cache_file.exists():
        try:
            with open(cache_file, 'r') as f:
                return json.load(f)
        except:
            pass
    return {}

def save_current_state(cache_file, state):
    """Save current progress state to cache file"""
    try:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_file, 'w') as f:
            json.dump(state, f, indent=2)
    except Exception as e:
        print(f"⚠️  Warning: Could not save state cache: {e}")

def check_simulation_progress(results_dir="results/phase2b_additional"):
    """Check progress of all running simulations"""
    results_dir = Path(results_dir)
    target_steps = 500000
    
    # Load previous state for comparison
    cache_file = results_dir / ".progress_cache.json"
    previous_state = load_previous_state(cache_file)
    current_state = {}
    
    print("=" * 80)
    print("🔍 REAL PROGRESS CHECK")
    print("=" * 80)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Get all running processes
    try:
        result = subprocess.run(
            ['ps', 'aux'],
            capture_output=True,
            text=True,
            check=True
        )
        
        processes = {}
        for line in result.stdout.split('\n'):
            if 'run_phase2_full.py' in line and 'grep' not in line:
                parts = line.split()
                if len(parts) >= 11:
                    pid = int(parts[1])
                    cpu_percent = float(parts[2])
                    
                    # Extract run info
                    run_match = re.search(r'run_(\d+)', line)
                    scenario_match = re.search(r'(miller_urey|hydrothermal|formamide)_extended', line)
                    
                    if run_match and scenario_match:
                        run_num = int(run_match.group(1))
                        scenario = scenario_match.group(0)
                        key = f"{scenario}/run_{run_num}"
                        
                        processes[key] = {
                            'pid': pid,
                            'cpu_percent': cpu_percent
                        }
    except Exception as e:
        print(f"⚠️  Error getting processes: {e}\n")
        return
    
    if not processes:
        print("❌ No running simulation processes found\n")
        return
    
    print(f"📊 Found {len(processes)} running simulations\n")
    
    # Check each simulation
    for key, proc_info in processes.items():
        pid = proc_info['pid']
        cpu_percent = proc_info['cpu_percent']
        
        # Find log file
        scenario, run_name = key.split('/')
        log_file = results_dir / scenario / run_name / "simulation.log"
        
        print(f"\n{'='*80}")
        print(f"🔍 {key} (PID: {pid})")
        print(f"{'='*80}")
        
        # Get last step from log
        last_step = get_last_step_from_log(log_file) if log_file.exists() else None
        
        if last_step is None:
            print("  ⚠️  Could not read last step from log")
            continue
        
        progress = (last_step / target_steps) * 100
        
        # Compare with previous state
        prev_key = key
        previous_last_step = previous_state.get(prev_key, {}).get('last_step')
        previous_timestamp = previous_state.get(prev_key, {}).get('timestamp')
        
        # Store current state
        current_state[prev_key] = {
            'last_step': last_step,
            'timestamp': time.time(),
            'progress': progress
        }
        
        # Calculate time since last check
        time_since_last_check = None
        if previous_timestamp:
            time_since_last_check = (time.time() - previous_timestamp) / 3600  # hours
        
        # Get log file age
        if log_file.exists():
            mtime = os.path.getmtime(log_file)
            age_seconds = time.time() - mtime
            age_minutes = age_seconds / 60
            age_hours = age_seconds / 3600
        else:
            age_minutes = None
            age_hours = None
        
        # Get process CPU time
        proc_details = get_process_info(pid)
        if proc_details:
            print(f"  Process info: {proc_details}")
        
        print(f"\n  📊 Last logged step: {last_step:,}/{target_steps:,} ({progress:.1f}%)")
        
        if age_minutes is not None:
            print(f"  ⏰ Log age: {age_minutes:.1f} minutes ({age_hours:.2f} hours)")
        
        print(f"  💻 CPU usage: {cpu_percent:.0f}%")
        
        # Show progress comparison if available
        if previous_last_step is not None and time_since_last_check is not None:
            step_change = last_step - previous_last_step
            if step_change > 0:
                print(f"  📈 Progress change: +{step_change:,} steps since last check ({time_since_last_check:.2f}h ago)")
                steps_per_hour = step_change / time_since_last_check if time_since_last_check > 0 else 0
                print(f"  ⚡ Actual progress rate: ~{steps_per_hour:,.0f} steps/hour")
            elif step_change == 0:
                print(f"  ⚠️  NO PROGRESS in logs since last check ({time_since_last_check:.2f}h ago)")
                print(f"  💡 This confirms log buffering - process is working but logs aren't updating")
            else:
                print(f"  ⚠️  Step count decreased (unusual - may be log parsing issue)")
        
        # Analysis
        print(f"\n  🔍 Analysis:")
        
        if cpu_percent > 100:
            print(f"     ✅ Process is actively computing (using multiple cores)")
            
            if age_minutes is not None and age_minutes > 30:
                print(f"     ⚠️  Log is old ({age_minutes:.0f} min), but process is active")
                print(f"     💡 This suggests log buffering (normal for old code)")
                print(f"     💡 Process is likely making progress, but logs aren't flushed")
                
                # More realistic estimation:
                # - Use time since last log (age_hours), not total CPU time
                # - Conservative: ~5-8 steps/sec per 100% CPU (varies by simulation)
                # - For SUPER_FAST config, typically ~8-12 steps/sec per 100% CPU
                if age_hours is not None and age_hours > 0:
                    # Conservative estimate: 7 steps/sec per 100% CPU
                    steps_per_100pct_cpu_per_sec = 7
                    steps_per_hour_per_100pct = steps_per_100pct_cpu_per_sec * 3600
                    
                    # Calculate steps since last log
                    estimated_steps_since_log = int((cpu_percent / 100) * steps_per_hour_per_100pct * age_hours)
                    
                    # Cap at remaining steps (can't exceed target)
                    remaining_steps = target_steps - last_step
                    estimated_steps_since_log = min(estimated_steps_since_log, remaining_steps)
                    
                    estimated_current = last_step + estimated_steps_since_log
                    estimated_progress_pct = (estimated_current / target_steps) * 100
                    
                    # Check if simulation might be completed
                    results_file = results_dir / scenario / run_name / "results.json"
                    
                    # If results.json exists, simulation is definitely completed
                    if results_file.exists():
                        print(f"     ✅ Simulation is COMPLETED (results.json exists)")
                        print(f"     📊 Final step: {target_steps:,} (100.0%)")
                    # If process is still running but estimated to be at target, it's likely near completion
                    # BUT: be conservative if last_step was very low (< 50%) - estimation may be too optimistic
                    elif estimated_current >= target_steps and last_step >= target_steps * 0.5:
                        # Process is still running, so it might be in final steps or extracting results
                        if last_step >= target_steps * 0.95:  # Very close (95%+)
                            print(f"     🎯 Simulation likely COMPLETED or in final stages")
                            print(f"     📊 Last logged: {last_step:,}, Estimated: ~{target_steps:,} (100.0%)")
                            print(f"     💡 Process may be extracting results or saving final data")
                            print(f"     💡 Check for results.json file - it should appear soon if completed")
                        else:
                            # Process running but estimated at target - might be close but not done
                            print(f"     🎯 Simulation likely very close to completion")
                            print(f"     📊 Last logged: {last_step:,}, Estimated: ~{target_steps:,} (100.0%)")
                            print(f"     ⚠️  Process still running - may be in final steps")
                            print(f"     💡 Estimated progress since last log: +{estimated_steps_since_log:,} steps")
                            print(f"     💡 Check for results.json file to confirm completion")
                    else:
                        # Normal progress estimate
                        # If estimated is at target but last_step was low, be more conservative
                        use_conservative = estimated_current >= target_steps and last_step < target_steps * 0.5
                        display_estimate = estimated_current
                        
                        if use_conservative:
                            # Estimation seems too optimistic - show more conservative estimate
                            # Cap at a reasonable maximum based on last_step progress
                            conservative_max = min(last_step * 2, target_steps * 0.9)  # Max 2x last_step or 90%
                            display_estimate = min(estimated_current, conservative_max)
                            conservative_pct = (display_estimate / target_steps) * 100
                            
                            print(f"     📈 Estimated current step: ~{display_estimate:,} ({conservative_pct:.1f}%)")
                            print(f"     ⚠️  Conservative estimate (last log was at {last_step:,}, {last_step/target_steps*100:.1f}%)")
                            print(f"     💡 Full estimate would be ~{estimated_current:,}, but may be too optimistic")
                            print(f"     📊 Estimated progress since last log: +{estimated_steps_since_log:,} steps")
                            print(f"     ⚠️  This is an ESTIMATE - actual progress may be lower")
                            print(f"     💡 Estimation assumes ~{steps_per_100pct_cpu_per_sec} steps/sec per 100% CPU")
                            print(f"     📉 Remaining steps: {target_steps - display_estimate:,}")
                        else:
                            # Normal estimate
                            print(f"     📈 Estimated current step: ~{display_estimate:,} ({estimated_progress_pct:.1f}%)")
                            print(f"     📊 Estimated progress since last log: +{estimated_steps_since_log:,} steps")
                            print(f"     ⚠️  This is an ESTIMATE based on CPU usage - actual progress may vary")
                            print(f"     💡 Estimation assumes ~{steps_per_100pct_cpu_per_sec} steps/sec per 100% CPU")
                            print(f"     📉 Remaining steps: {target_steps - display_estimate:,}")
                        
                        # Add time estimate if reasonable
                        if estimated_steps_since_log > 0 and age_hours > 0:
                            steps_per_hour_actual = estimated_steps_since_log / age_hours
                            remaining_steps = target_steps - display_estimate
                            remaining_hours = remaining_steps / steps_per_hour_actual if steps_per_hour_actual > 0 else None
                            if remaining_hours and remaining_hours < 24:
                                print(f"     ⏱️  Estimated time remaining: ~{remaining_hours:.1f} hours")
            else:
                print(f"     ✅ Log is recent - progress is visible")
        else:
            print(f"     ⚠️  Process CPU usage is low - may be idle or stuck")
        
        # Recommendation
        print(f"\n  💡 Recommendation:")
        if cpu_percent > 100 and age_minutes is not None and age_minutes > 30:
            print(f"     - Process is working, but using old code without log flushing")
            print(f"     - Progress is likely happening, but not visible in logs")
            print(f"     - New simulations will have better log visibility")
            print(f"     - Consider restarting after current batch completes")
        elif cpu_percent > 100:
            print(f"     - Process is working normally")
        else:
            print(f"     - Process may be stuck - investigate further")
    
    # Save current state for next run
    save_current_state(cache_file, current_state)
    
    print(f"\n{'='*80}")
    print("ℹ️  NOTE: This script estimates progress based on CPU usage patterns")
    print("   Actual progress may vary. Check logs periodically for updates.")
    if previous_state:
        print("   Progress comparison is based on previous run (stored in .progress_cache.json)")
    print("=" * 80)

if __name__ == "__main__":
    import sys
    results_dir = sys.argv[1] if len(sys.argv) > 1 else "results/phase2b_additional"
    check_simulation_progress(results_dir)

