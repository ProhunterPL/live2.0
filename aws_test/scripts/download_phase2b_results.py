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
    
    # Check if Phase 2B is running
    check_cmd = [
        "ssh", "-i", key_path,
        f"ubuntu@{host}",
        "cd ~/live2.0/aws_test && find results/phase2b_additional -name 'summary.txt' | wc -l"
    ]
    
    try:
        result = subprocess.run(check_cmd, check=True, capture_output=True, text=True)
        completed_runs = int(result.stdout.strip())
        print(f"📊 Completed runs: {completed_runs}")
        
        # Check if still running
        running_cmd = [
            "ssh", "-i", key_path,
            f"ubuntu@{host}",
            "ps aux | grep run_phase2b | grep -v grep | wc -l"
        ]
        
        result = subprocess.run(running_cmd, check=True, capture_output=True, text=True)
        running_processes = int(result.stdout.strip())
        print(f"🔄 Running processes: {running_processes}")
        
        if running_processes > 0:
            print("⏳ Phase 2B still running...")
            return "running"
        elif completed_runs >= 30:
            print("✅ Phase 2B completed!")
            return "completed"
        else:
            print("⚠️ Phase 2B may have issues")
            return "partial"
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Status check failed: {e.stderr}")
        return "error"

def analyze_downloaded_results(local_dir):
    """Analyze downloaded results"""
    
    print(f"📊 Analyzing downloaded results in {local_dir}")
    
    local_path = Path(local_dir)
    if not local_path.exists():
        print(f"❌ Local directory {local_dir} not found")
        return False
    
    # Run analysis
    analysis_cmd = [
        "python", "scripts/analyze_additional_results.py",
        "--phase2b-dir", str(local_path),
        "--phase2a-dir", "aws_test"
    ]
    
    try:
        result = subprocess.run(analysis_cmd, check=True, capture_output=True, text=True)
        print("✅ Analysis completed")
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Analysis failed: {e.stderr}")
        return False

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
