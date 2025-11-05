#!/usr/bin/env python3
"""
AWS Phase 2B Results Downloader
===============================

Downloads Phase 2B results from AWS instance and analyzes them.
"""

import os
import sys
import json
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

def download_results(host, key_path, local_dir="results/phase2b_local"):
    """Download Phase 2B results from AWS instance"""
    
    print(f"📥 Downloading Phase 2B results from {host}")
    print(f"📁 Local directory: {local_dir}")
    
    # Create local directory
    local_path = Path(local_dir)
    local_path.mkdir(parents=True, exist_ok=True)
    
    # Download results
    remote_path = "~/live2.0/aws_test/results/phase2b_additional"
    
    cmd = [
        "scp", "-r", "-i", key_path,
        f"ubuntu@{host}:{remote_path}/*",
        str(local_path)
    ]
    
    print(f"🔗 Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✅ Download completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Download failed: {e.stderr}")
        return False

def check_remote_status(host, key_path):
    """Check status of Phase 2B runs on AWS"""
    
    print(f"🔍 Checking Phase 2B status on {host}")
    
    # Check multiple indicators of completion
    check_cmd = [
        "ssh", "-i", key_path,
        f"ubuntu@{host}",
        """cd ~/live2.0/aws_test/results/phase2b_additional && \
        (find . -name 'summary.txt' | wc -l; \
         find . -name 'results.json' | wc -l; \
         find . -name '*.md' -type f | wc -l; \
         ls -d */run_*/ 2>/dev/null | wc -l)"""
    ]
    
    try:
        result = subprocess.run(check_cmd, check=True, capture_output=True, text=True)
        lines = result.stdout.strip().split('\n')
        
        summary_count = int(lines[0]) if len(lines) > 0 else 0
        results_count = int(lines[1]) if len(lines) > 1 else 0
        reports_count = int(lines[2]) if len(lines) > 2 else 0
        run_dirs = int(lines[3]) if len(lines) > 3 else 0
        
        completed_runs = max(summary_count, results_count, run_dirs)
        
        print(f"📊 Completed runs (summary.txt): {summary_count}")
        print(f"📊 Completed runs (results.json): {results_count}")
        print(f"📊 Run directories: {run_dirs}")
        print(f"📊 Reports (MD files): {reports_count}")
        
        # Check if still running
        running_cmd = [
            "ssh", "-i", key_path,
            f"ubuntu@{host}",
            "ps aux | grep -E 'run_phase2b|run_phase2_full' | grep -v grep | wc -l"
        ]
        
        result = subprocess.run(running_cmd, check=True, capture_output=True, text=True)
        running_processes = int(result.stdout.strip())
        print(f"🔄 Running processes: {running_processes}")
        
        # Check for report files (indicates completion)
        if reports_count >= 3:  # At least 3 report files
            print("✅ Phase 2B reports found - looks completed!")
            return "completed"
        
        if running_processes > 0:
            print("⏳ Phase 2B still running...")
            return "running"
        elif completed_runs >= 30:
            print("✅ Phase 2B completed!")
            return "completed"
        elif completed_runs > 0:
            print(f"⚠️ Phase 2B partially completed ({completed_runs}/30)")
            return "partial"
        else:
            print("⚠️ Phase 2B may have issues or not started yet")
            return "error"
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Status check failed: {e.stderr}")
        print(f"Debug output: {e.stdout}")
        return "error"

def analyze_downloaded_results(local_dir):
    """Analyze downloaded results"""
    
    print(f"📊 Analyzing downloaded results in {local_dir}")
    
    local_path = Path(local_dir)
    if not local_path.exists():
        print(f"❌ Local directory {local_dir} not found")
        return False
    
    # Run analysis - use aws_test/scripts path
    script_path = Path(__file__).parent.parent / "scripts" / "analyze_additional_results.py"
    if not script_path.exists():
        # Try alternative path
        script_path = Path(__file__).parent.parent.parent / "aws_test" / "scripts" / "analyze_additional_results.py"
    
    if not script_path.exists():
        print(f"⚠️ Analysis script not found at {script_path}")
        print("📄 Skipping analysis - check results manually")
        return True  # Don't fail, just skip analysis
    
    analysis_cmd = [
        "python", str(script_path),
        "--phase2b-dir", str(local_path),
        "--phase2a-dir", str(Path(__file__).parent.parent.parent / "aws_test")
    ]
    
    try:
        result = subprocess.run(analysis_cmd, check=True, capture_output=True, text=True)
        print("✅ Analysis completed")
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Analysis failed: {e.stderr}")
        print("📄 Check results manually - analysis script had issues")
        return True  # Don't fail download because of analysis issues

def main():
    parser = argparse.ArgumentParser(description="AWS Phase 2B Results Downloader")
    parser.add_argument("--host", required=True, help="AWS instance IP or hostname")
    parser.add_argument("--key", required=True, help="Path to SSH private key")
    parser.add_argument("--local-dir", default="results/phase2b_local",
                       help="Local directory for downloaded results")
    parser.add_argument("--status-only", action="store_true",
                       help="Check status only, don't download")
    parser.add_argument("--analyze-only", action="store_true",
                       help="Analyze existing local results only")
    
    args = parser.parse_args()
    
    print("🔍 AWS Phase 2B Results Downloader")
    print("=" * 50)
    
    if args.analyze_only:
        # Analyze existing local results
        analyze_downloaded_results(args.local_dir)
        return
    
    # Check remote status
    status = check_remote_status(args.host, args.key)
    
    if args.status_only:
        return
    
    if status == "running":
        print("⏳ Phase 2B still running. Download partial results? (y/N)")
        response = input().strip().lower()
        if response != 'y':
            print("❌ Download cancelled")
            return
    
    # Download results
    if download_results(args.host, args.key, args.local_dir):
        # Analyze downloaded results
        analyze_downloaded_results(args.local_dir)
        
        print("\n🎉 Download and analysis complete!")
        print(f"📁 Results in: {args.local_dir}")
        print("📄 Check analysis reports for details")

if __name__ == "__main__":
    main()
