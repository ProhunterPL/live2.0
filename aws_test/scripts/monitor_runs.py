#!/usr/bin/env python3
"""
Phase 2B Monitoring Tool
========================

Real-time monitoring tool for Phase 2B additional runs.
Shows progress, performance metrics, and system status.
"""

import os
import sys
import json
import time
import psutil
import argparse
from pathlib import Path
from datetime import datetime, timedelta
import logging

class Phase2BMonitor:
    def __init__(self, results_dir="results/phase2b_additional"):
        self.results_dir = Path(results_dir)
        self.start_time = datetime.now()
        
        # Setup logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
    def get_system_status(self):
        """Get current system status"""
        cpu_percent = psutil.cpu_percent(interval=1)
        memory = psutil.virtual_memory()
        disk = psutil.disk_usage('/')
        
        return {
            "cpu_percent": cpu_percent,
            "memory_percent": memory.percent,
            "memory_available_gb": memory.available / (1024**3),
            "disk_percent": disk.percent,
            "disk_free_gb": disk.free / (1024**3)
        }
    
    def get_simulation_progress(self):
        """Get simulation progress"""
        progress = {
            "total_runs": 30,
            "completed_runs": 0,
            "failed_runs": 0,
            "running_runs": 0,
            "scenarios": {}
        }
        
        scenarios = ["miller_urey_extended", "hydrothermal_extended", "formamide_extended"]
        
        for scenario in scenarios:
            scenario_dir = self.results_dir / scenario
            if not scenario_dir.exists():
                progress["scenarios"][scenario] = {
                    "completed": 0,
                    "failed": 0,
                    "running": 0,
                    "total": 10
                }
                continue
            
            completed = 0
            failed = 0
            running = 0
            
            for run_dir in scenario_dir.iterdir():
                if not run_dir.is_dir():
                    continue
                
                # Check if run is complete
                summary_file = run_dir / "summary.txt"
                results_file = run_dir / "results.json"
                
                if summary_file.exists() and results_file.exists():
                    completed += 1
                    progress["completed_runs"] += 1
                else:
                    # Check if simulation is still running
                    log_file = run_dir / "simulation.log"
                    if log_file.exists():
                        # Check if log was updated recently (within last 5 minutes)
                        if time.time() - log_file.stat().st_mtime < 300:
                            running += 1
                            progress["running_runs"] += 1
                        else:
                            failed += 1
                            progress["failed_runs"] += 1
                    else:
                        failed += 1
                        progress["failed_runs"] += 1
            
            progress["scenarios"][scenario] = {
                "completed": completed,
                "failed": failed,
                "running": running,
                "total": 10
            }
        
        return progress
    
    def get_recent_activity(self):
        """Get recent activity from logs"""
        activity = []
        
        # Check for recent log files
        log_dir = self.results_dir / "logs"
        if log_dir.exists():
            for log_file in log_dir.glob("*.log"):
                if time.time() - log_file.stat().st_mtime < 3600:  # Last hour
                    try:
                        with open(log_file, 'r') as f:
                            lines = f.readlines()
                            # Get last few lines
                            recent_lines = lines[-10:] if len(lines) > 10 else lines
                            activity.extend(recent_lines)
                    except:
                        pass
        
        return activity[-20:]  # Last 20 lines
    
    def display_status(self):
        """Display current status"""
        os.system('clear' if os.name == 'posix' else 'cls')
        
        print("=" * 80)
        print("🔍 PHASE 2B MONITORING DASHBOARD")
        print("=" * 80)
        print(f"📅 Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"⏱️  Uptime: {datetime.now() - self.start_time}")
        print()
        
        # System status
        system = self.get_system_status()
        print("💻 SYSTEM STATUS")
        print("-" * 40)
        print(f"CPU Usage: {system['cpu_percent']:.1f}%")
        print(f"Memory: {system['memory_percent']:.1f}% ({system['memory_available_gb']:.1f} GB available)")
        print(f"Disk: {system['disk_percent']:.1f}% ({system['disk_free_gb']:.1f} GB free)")
        print()
        
        # Simulation progress
        progress = self.get_simulation_progress()
        print("📊 SIMULATION PROGRESS")
        print("-" * 40)
        print(f"Total Runs: {progress['completed_runs']}/{progress['total_runs']} completed")
        print(f"Failed: {progress['failed_runs']}")
        print(f"Running: {progress['running_runs']}")
        print(f"Success Rate: {progress['completed_runs']/(progress['completed_runs']+progress['failed_runs'])*100:.1f}%" if (progress['completed_runs']+progress['failed_runs']) > 0 else "N/A")
        print()
        
        # Scenario breakdown
        print("🎯 SCENARIO BREAKDOWN")
        print("-" * 40)
        for scenario, data in progress["scenarios"].items():
            print(f"{scenario.replace('_', ' ').title()}:")
            print(f"  ✅ Completed: {data['completed']}/{data['total']}")
            print(f"  ❌ Failed: {data['failed']}")
            print(f"  🔄 Running: {data['running']}")
            print()
        
        # Recent activity
        activity = self.get_recent_activity()
        if activity:
            print("📝 RECENT ACTIVITY")
            print("-" * 40)
            for line in activity[-10:]:  # Last 10 lines
                print(line.strip())
            print()
        
        print("=" * 80)
        print("Press Ctrl+C to exit")
        print("=" * 80)
    
    def monitor_continuous(self, interval=30):
        """Monitor continuously with specified interval"""
        try:
            while True:
                self.display_status()
                time.sleep(interval)
        except KeyboardInterrupt:
            print("\n👋 Monitoring stopped")
    
    def generate_report(self):
        """Generate monitoring report"""
        progress = self.get_simulation_progress()
        system = self.get_system_status()
        
        report = {
            "timestamp": datetime.now().isoformat(),
            "uptime_hours": (datetime.now() - self.start_time).total_seconds() / 3600,
            "system_status": system,
            "simulation_progress": progress,
            "estimated_completion": self.estimate_completion_time(progress)
        }
        
        # Save report
        report_file = self.results_dir / "monitoring_report.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"📄 Monitoring report saved: {report_file}")
        return report
    
    def estimate_completion_time(self, progress):
        """Estimate completion time"""
        if progress["running_runs"] == 0:
            return "No active simulations"
        
        # Rough estimate: 2-4 hours per simulation
        avg_time_per_run = 3 * 3600  # 3 hours in seconds
        remaining_runs = progress["total_runs"] - progress["completed_runs"] - progress["failed_runs"]
        
        if remaining_runs <= 0:
            return "All simulations completed"
        
        estimated_seconds = remaining_runs * avg_time_per_run
        estimated_time = datetime.now() + timedelta(seconds=estimated_seconds)
        
        return {
            "remaining_runs": remaining_runs,
            "estimated_completion": estimated_time.isoformat(),
            "estimated_hours": estimated_seconds / 3600
        }

def main():
    parser = argparse.ArgumentParser(description="Phase 2B Monitoring Tool")
    parser.add_argument("--results-dir", default="results/phase2b_additional",
                       help="Results directory to monitor")
    parser.add_argument("--interval", type=int, default=30,
                       help="Monitoring interval in seconds")
    parser.add_argument("--report-only", action="store_true",
                       help="Generate report and exit")
    
    args = parser.parse_args()
    
    monitor = Phase2BMonitor(args.results_dir)
    
    if args.report_only:
        monitor.generate_report()
    else:
        monitor.monitor_continuous(args.interval)

if __name__ == "__main__":
    main()
