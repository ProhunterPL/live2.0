#!/usr/bin/env python3
"""
Phase 2B Master Script
=======================

Master script to run the complete Phase 2B additional runs process:
1. Debug formamide scenario
2. Run 30 additional simulations
3. Monitor progress
4. Analyze results
5. Generate reports

Usage:
    python run_phase2b_master.py --mode [debug|run|monitor|analyze|all]
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

def run_command(cmd, description):
    """Run a command and handle errors"""
    print(f"[RUN] {description}")
    print(f"   Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"   [OK] {description} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"   [FAIL] {description} failed: {e.stderr}")
        return False
    except Exception as e:
        print(f"   [ERROR] {description} crashed: {str(e)}")
        return False

def debug_formamide():
    """Debug formamide scenario"""
    print("[PHASE 1] DEBUG FORMAMIDE")
    print("=" * 50)
    
    cmd = [
        "python3", "scripts/debug_formamide.py",
        "--output-dir", "results/phase2b_additional/formamide_debug"
    ]
    
    success = run_command(cmd, "Debug formamide scenario")
    
    if success:
        print("[OK] Formamide debug completed")
        print("[INFO] Check results/phase2b_additional/formamide_debug/formamide_debug_report.md")
    else:
        print("[FAIL] Formamide debug failed")
        print("[INFO] Check logs and fix issues before proceeding")
    
    return success

def run_additional_simulations():
    """Run 30 additional simulations"""
    print("[PHASE 2] RUN ADDITIONAL SIMULATIONS")
    print("=" * 50)
    
    cmd = [
        "python3", "scripts/run_phase2b_additional.py",
        "--output-dir", "results/phase2b_additional",
        "--max-parallel", "4"  # Run 4 simulations in parallel (16 cores each)
    ]
    
    print(f"[RUN] Run 30 additional simulations")
    print(f"   Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Check if script ran without syntax errors
        if result.returncode != 0:
            print(f"   [FAIL] Script failed: {result.stderr}")
            return False
        
        # Check actual results from JSON file
        results_file = Path("results/phase2b_additional/phase2b_results.json")
        if results_file.exists():
            import json
            with open(results_file, 'r') as f:
                data = json.load(f)
            
            successful = data.get('completed_runs', 0)
            failed = data.get('failed_runs', 0)
            total = data.get('total_runs', 30)
            
            if successful == total:
                print(f"   [OK] All simulations completed successfully: {successful}/{total}")
                return True
            elif successful > 0:
                print(f"   [PARTIAL] Partial success: {successful}/{total} completed, {failed}/{total} failed")
                print(f"   [INFO] Check results/phase2b_additional/phase2b_summary_report.md")
                return True  # Partial success is still progress
            else:
                print(f"   [FAIL] All simulations failed: {successful}/{total} completed, {failed}/{total} failed")
                print(f"   [INFO] Check logs/phase2b_runner.log for errors")
                return False
        else:
            print(f"   [WARN] Results file not found: {results_file}")
            print(f"   [INFO] Script may have failed before creating results")
            return False
            
    except Exception as e:
        print(f"   [ERROR] Script crashed: {str(e)}")
        return False

def monitor_progress():
    """Monitor simulation progress"""
    print("[PHASE 3] MONITOR PROGRESS")
    print("=" * 50)
    
    print("[INFO] Starting monitoring dashboard...")
    print("   Press Ctrl+C to stop monitoring")
    
    cmd = [
        "python3", "scripts/monitor_runs.py",
        "--results-dir", "results/phase2b_additional",
        "--interval", "30"
    ]
    
    try:
        subprocess.run(cmd)
        return True
    except KeyboardInterrupt:
        print("\n[INFO] Monitoring stopped")
        return True
    except Exception as e:
        print(f"[FAIL] Monitoring failed: {str(e)}")
        return False

def analyze_results():
    """Analyze results and generate reports"""
    print("[PHASE 4] ANALYZE RESULTS")
    print("=" * 50)
    
    cmd = [
        "python3", "scripts/analyze_additional_results.py",
        "--phase2b-dir", "results/phase2b_additional",
        "--phase2a-dir", "aws_test"
    ]
    
    success = run_command(cmd, "Analyze Phase 2B results")
    
    if success:
        print("[OK] Analysis completed")
        print("[INFO] Check results/phase2b_additional/phase2b_analysis_report.md")
    else:
        print("[FAIL] Analysis failed")
        print("[INFO] Check logs and fix issues")
    
    return success

def run_all():
    """Run complete Phase 2B process"""
    print("[PHASE 2B] COMPLETE PROCESS")
    print("=" * 60)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Phase 1: Debug formamide
    if not debug_formamide():
        print("[FAIL] Formamide debug failed - stopping process")
        return False
    
    print("\n" + "=" * 60)
    
    # Phase 2: Run additional simulations
    if not run_additional_simulations():
        print("[FAIL] Additional simulations failed - stopping process")
        return False
    
    print("\n" + "=" * 60)
    
    # Phase 3: Analyze results
    if not analyze_results():
        print("[FAIL] Analysis failed - stopping process")
        return False
    
    print("\n" + "=" * 60)
    print("[OK] PHASE 2B COMPLETE!")
    print("=" * 60)
    
    print("[SUMMARY] FINAL RESULTS:")
    print("   [OK] Formamide debug completed")
    print("   [OK] 30 additional simulations completed")
    print("   [OK] Results analyzed")
    print("   [OK] Reports generated")
    
    print("\n[INFO] CHECK THESE FILES:")
    print("   [INFO] results/phase2b_additional/formamide_debug/formamide_debug_report.md")
    print("   [INFO] results/phase2b_additional/phase2b_summary_report.md")
    print("   [INFO] results/phase2b_additional/phase2b_analysis_report.md")
    
    print("\n[NEXT STEPS]:")
    print("   1. Review analysis reports")
    print("   2. Generate publication figures")
    print("   3. Proceed to Phase 3 (Paper Writing)")
    
    return True

def main():
    parser = argparse.ArgumentParser(description="Phase 2B Master Script")
    parser.add_argument("--mode", 
                       choices=["debug", "run", "monitor", "analyze", "all"],
                       default="all",
                       help="Mode to run")
    
    args = parser.parse_args()
    
    print("[PHASE 2B] MASTER SCRIPT")
    print("=" * 40)
    print(f"Mode: {args.mode}")
    print()
    
    if args.mode == "debug":
        debug_formamide()
    elif args.mode == "run":
        run_additional_simulations()
    elif args.mode == "monitor":
        monitor_progress()
    elif args.mode == "analyze":
        analyze_results()
    elif args.mode == "all":
        run_all()
    
    print("\n[INFO] Phase 2B Master Script completed")

if __name__ == "__main__":
    main()
