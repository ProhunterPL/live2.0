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
    print(f"🚀 {description}")
    print(f"   Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"   ✅ {description} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"   ❌ {description} failed: {e.stderr}")
        return False
    except Exception as e:
        print(f"   💥 {description} crashed: {str(e)}")
        return False

def debug_formamide():
    """Debug formamide scenario"""
    print("🔍 PHASE 1: DEBUG FORMAMIDE")
    print("=" * 50)
    
    cmd = [
        "python", "scripts/debug_formamide.py",
        "--output-dir", "results/phase2b_additional/formamide_debug"
    ]
    
    success = run_command(cmd, "Debug formamide scenario")
    
    if success:
        print("✅ Formamide debug completed")
        print("📄 Check results/phase2b_additional/formamide_debug/formamide_debug_report.md")
    else:
        print("❌ Formamide debug failed")
        print("🔧 Check logs and fix issues before proceeding")
    
    return success

def run_additional_simulations():
    """Run 30 additional simulations"""
    print("🚀 PHASE 2: RUN ADDITIONAL SIMULATIONS")
    print("=" * 50)
    
    cmd = [
        "python", "scripts/run_phase2b_additional.py",
        "--output-dir", "results/phase2b_additional"
    ]
    
    success = run_command(cmd, "Run 30 additional simulations")
    
    if success:
        print("✅ Additional simulations completed")
        print("📄 Check results/phase2b_additional/phase2b_summary_report.md")
    else:
        print("❌ Additional simulations failed")
        print("🔧 Check logs and retry failed runs")
    
    return success

def monitor_progress():
    """Monitor simulation progress"""
    print("📊 PHASE 3: MONITOR PROGRESS")
    print("=" * 50)
    
    print("🔍 Starting monitoring dashboard...")
    print("   Press Ctrl+C to stop monitoring")
    
    cmd = [
        "python", "scripts/monitor_runs.py",
        "--results-dir", "results/phase2b_additional",
        "--interval", "30"
    ]
    
    try:
        subprocess.run(cmd)
        return True
    except KeyboardInterrupt:
        print("\n👋 Monitoring stopped")
        return True
    except Exception as e:
        print(f"❌ Monitoring failed: {str(e)}")
        return False

def analyze_results():
    """Analyze results and generate reports"""
    print("📈 PHASE 4: ANALYZE RESULTS")
    print("=" * 50)
    
    cmd = [
        "python", "scripts/analyze_additional_results.py",
        "--phase2b-dir", "results/phase2b_additional",
        "--phase2a-dir", "aws_test"
    ]
    
    success = run_command(cmd, "Analyze Phase 2B results")
    
    if success:
        print("✅ Analysis completed")
        print("📄 Check results/phase2b_additional/phase2b_analysis_report.md")
    else:
        print("❌ Analysis failed")
        print("🔧 Check logs and fix issues")
    
    return success

def run_all():
    """Run complete Phase 2B process"""
    print("🎯 PHASE 2B: COMPLETE PROCESS")
    print("=" * 60)
    print(f"📅 Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Phase 1: Debug formamide
    if not debug_formamide():
        print("❌ Formamide debug failed - stopping process")
        return False
    
    print("\n" + "=" * 60)
    
    # Phase 2: Run additional simulations
    if not run_additional_simulations():
        print("❌ Additional simulations failed - stopping process")
        return False
    
    print("\n" + "=" * 60)
    
    # Phase 3: Analyze results
    if not analyze_results():
        print("❌ Analysis failed - stopping process")
        return False
    
    print("\n" + "=" * 60)
    print("🎉 PHASE 2B COMPLETE!")
    print("=" * 60)
    
    print("📊 FINAL RESULTS:")
    print("   ✅ Formamide debug completed")
    print("   ✅ 30 additional simulations completed")
    print("   ✅ Results analyzed")
    print("   ✅ Reports generated")
    
    print("\n📄 CHECK THESE FILES:")
    print("   📄 results/phase2b_additional/formamide_debug/formamide_debug_report.md")
    print("   📄 results/phase2b_additional/phase2b_summary_report.md")
    print("   📄 results/phase2b_additional/phase2b_analysis_report.md")
    
    print("\n🎯 NEXT STEPS:")
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
    
    print("🎯 PHASE 2B MASTER SCRIPT")
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
    
    print("\n👋 Phase 2B Master Script completed")

if __name__ == "__main__":
    main()
