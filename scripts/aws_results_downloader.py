#!/usr/bin/env python3
"""
AWS Results Downloader & Organizer
==================================

Automatically download and organize simulation results from AWS instance.
Handles batch downloading, verification, and organization.

Usage:
    python scripts/aws_results_downloader.py --host <aws-ip> --key <ssh-key>
    
Author: Live 2.0 Team
Date: October 2025
"""

import argparse
import subprocess
import json
import os
from pathlib import Path
from datetime import datetime
import sys

class AWSResultsDownloader:
    """Download and organize results from AWS"""
    
    def __init__(self, aws_host, ssh_key, remote_base='~/live2.0/results', 
                 local_base='./results/aws_batch'):
        self.aws_host = aws_host
        self.ssh_key = ssh_key
        self.remote_base = remote_base
        self.local_base = Path(local_base)
        self.local_base.mkdir(parents=True, exist_ok=True)
        
        self.download_log = self.local_base / 'download_log.json'
        self.load_log()
    
    def load_log(self):
        """Load download history"""
        if self.download_log.exists():
            with open(self.download_log) as f:
                self.log = json.load(f)
        else:
            self.log = {
                'downloads': [],
                'last_update': None,
                'total_files': 0,
                'total_size_mb': 0
            }
    
    def save_log(self):
        """Save download history"""
        with open(self.download_log, 'w') as f:
            json.dump(self.log, f, indent=2)
    
    def check_connection(self):
        """Test SSH connection to AWS"""
        print("🔍 Testing AWS connection...")
        cmd = f"ssh -i {self.ssh_key} ubuntu@{self.aws_host} 'echo OK'"
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, 
                                  text=True, timeout=10)
            if result.returncode == 0 and 'OK' in result.stdout:
                print("✅ Connection successful!")
                return True
            else:
                print(f"❌ Connection failed: {result.stderr}")
                return False
        except Exception as e:
            print(f"❌ Connection error: {e}")
            return False
    
    def list_remote_results(self):
        """List all completed simulations on AWS"""
        print("📋 Listing remote results...")
        
        # Count completed runs (have summary.txt)
        cmd = f"""ssh -i {self.ssh_key} ubuntu@{self.aws_host} '
        cd {self.remote_base} && 
        find . -name "summary.txt" -type f | 
        sed "s|/summary.txt||" | 
        sort
        '"""
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, 
                                  text=True, timeout=30)
            if result.returncode == 0:
                paths = [p.strip() for p in result.stdout.strip().split('\n') if p.strip()]
                print(f"✅ Found {len(paths)} completed simulations")
                return paths
            else:
                print(f"❌ Error listing results: {result.stderr}")
                return []
        except Exception as e:
            print(f"❌ Error: {e}")
            return []
    
    def get_simulation_info(self, remote_path):
        """Get info about a simulation"""
        cmd = f"""ssh -i {self.ssh_key} ubuntu@{self.aws_host} '
        cd {self.remote_base}/{remote_path} && 
        if [ -f summary.txt ]; then cat summary.txt; fi
        '"""
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, 
                                  text=True, timeout=10)
            if result.returncode == 0:
                return result.stdout
            return None
        except:
            return None
    
    def download_simulation(self, remote_path):
        """Download a single simulation"""
        # Parse scenario and run number
        parts = remote_path.strip('./').split('/')
        if len(parts) < 2:
            print(f"⚠️  Invalid path format: {remote_path}")
            return False
        
        scenario = parts[0]
        run_id = parts[1]
        
        local_dir = self.local_base / scenario / run_id
        local_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"📥 Downloading {scenario}/{run_id}...")
        
        # Use rsync for efficient transfer (only new/changed files)
        cmd = f"""rsync -avz --progress \
        -e "ssh -i {self.ssh_key}" \
        ubuntu@{self.aws_host}:{self.remote_base}/{remote_path}/ \
        {local_dir}/
        """
        
        try:
            result = subprocess.run(cmd, shell=True, timeout=600)
            if result.returncode == 0:
                print(f"✅ Downloaded {scenario}/{run_id}")
                
                # Verify download
                if (local_dir / 'summary.txt').exists():
                    # Log successful download
                    self.log['downloads'].append({
                        'scenario': scenario,
                        'run_id': run_id,
                        'path': str(local_dir),
                        'timestamp': datetime.now().isoformat()
                    })
                    self.log['total_files'] += 1
                    self.save_log()
                    return True
                else:
                    print(f"⚠️  Download incomplete (no summary.txt)")
                    return False
            else:
                print(f"❌ Download failed")
                return False
        except subprocess.TimeoutExpired:
            print(f"⏰ Download timed out (may be too large)")
            return False
        except Exception as e:
            print(f"❌ Error: {e}")
            return False
    
    def download_all(self, skip_existing=True):
        """Download all completed simulations"""
        remote_paths = self.list_remote_results()
        
        if not remote_paths:
            print("No results to download.")
            return
        
        downloaded = []
        skipped = []
        failed = []
        
        for i, path in enumerate(remote_paths, 1):
            print(f"\n[{i}/{len(remote_paths)}] Processing {path}")
            
            # Check if already downloaded
            parts = path.strip('./').split('/')
            if len(parts) >= 2:
                local_check = self.local_base / parts[0] / parts[1] / 'summary.txt'
                if skip_existing and local_check.exists():
                    print(f"⏭️  Already downloaded, skipping")
                    skipped.append(path)
                    continue
            
            # Download
            if self.download_simulation(path):
                downloaded.append(path)
            else:
                failed.append(path)
        
        # Summary
        print("\n" + "="*60)
        print("📊 DOWNLOAD SUMMARY")
        print("="*60)
        print(f"✅ Downloaded: {len(downloaded)}")
        print(f"⏭️  Skipped (existing): {len(skipped)}")
        print(f"❌ Failed: {len(failed)}")
        print(f"📁 Total local: {self.log['total_files']}")
        
        if failed:
            print("\n❌ Failed downloads:")
            for path in failed:
                print(f"  - {path}")
        
        self.log['last_update'] = datetime.now().isoformat()
        self.save_log()
        
        print(f"\n💾 Log saved to: {self.download_log}")
    
    def verify_downloads(self):
        """Verify integrity of downloaded files"""
        print("🔍 Verifying downloaded simulations...")
        
        scenarios = ['miller_urey', 'hydrothermal', 'formamide']
        valid = []
        invalid = []
        
        for scenario in scenarios:
            scenario_dir = self.local_base / scenario
            if not scenario_dir.exists():
                continue
            
            for run_dir in sorted(scenario_dir.iterdir()):
                if not run_dir.is_dir():
                    continue
                
                # Check required files
                required = ['summary.txt', 'config.yaml']
                has_all = all((run_dir / f).exists() for f in required)
                
                if has_all:
                    valid.append(str(run_dir.relative_to(self.local_base)))
                else:
                    invalid.append(str(run_dir.relative_to(self.local_base)))
        
        print(f"\n✅ Valid simulations: {len(valid)}")
        print(f"❌ Invalid simulations: {len(invalid)}")
        
        if invalid:
            print("\n⚠️  Invalid simulations (missing files):")
            for path in invalid:
                print(f"  - {path}")
        
        return valid, invalid
    
    def get_download_status(self):
        """Get current download status"""
        valid, invalid = self.verify_downloads()
        
        status = {
            'last_update': self.log.get('last_update'),
            'total_downloaded': self.log['total_files'],
            'valid_simulations': len(valid),
            'invalid_simulations': len(invalid),
            'by_scenario': {}
        }
        
        # Count by scenario
        for scenario in ['miller_urey', 'hydrothermal', 'formamide']:
            scenario_dir = self.local_base / scenario
            if scenario_dir.exists():
                count = len([d for d in scenario_dir.iterdir() 
                           if d.is_dir() and (d / 'summary.txt').exists()])
                status['by_scenario'][scenario] = count
            else:
                status['by_scenario'][scenario] = 0
        
        return status


def main():
    parser = argparse.ArgumentParser(
        description='Download simulation results from AWS',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download all new results
  python scripts/aws_results_downloader.py --host 54.123.45.67 --key ~/.ssh/aws_key.pem
  
  # Force re-download everything
  python scripts/aws_results_downloader.py --host 54.123.45.67 --key ~/.ssh/aws_key.pem --force
  
  # Just check status
  python scripts/aws_results_downloader.py --host 54.123.45.67 --key ~/.ssh/aws_key.pem --status-only
        """
    )
    
    parser.add_argument('--host', required=True, help='AWS instance IP or hostname')
    parser.add_argument('--key', required=True, help='Path to SSH private key')
    parser.add_argument('--remote-base', default='~/live2.0/results',
                       help='Remote results directory (default: ~/live2.0/results)')
    parser.add_argument('--local-base', default='./results/aws_batch',
                       help='Local download directory (default: ./results/aws_batch)')
    parser.add_argument('--force', action='store_true',
                       help='Re-download existing files')
    parser.add_argument('--status-only', action='store_true',
                       help='Only show download status, don\'t download')
    
    args = parser.parse_args()
    
    # Create downloader
    downloader = AWSResultsDownloader(
        aws_host=args.host,
        ssh_key=args.key,
        remote_base=args.remote_base,
        local_base=args.local_base
    )
    
    # Check connection
    if not downloader.check_connection():
        print("\n❌ Cannot connect to AWS. Check your credentials and host.")
        sys.exit(1)
    
    # Status only?
    if args.status_only:
        status = downloader.get_download_status()
        print("\n" + "="*60)
        print("📊 DOWNLOAD STATUS")
        print("="*60)
        print(f"Last update: {status['last_update'] or 'Never'}")
        print(f"Total downloaded: {status['total_downloaded']}")
        print(f"Valid simulations: {status['valid_simulations']}")
        print(f"Invalid simulations: {status['invalid_simulations']}")
        print("\nBy scenario:")
        for scenario, count in status['by_scenario'].items():
            print(f"  {scenario}: {count}")
        sys.exit(0)
    
    # Download all
    downloader.download_all(skip_existing=not args.force)
    
    # Verify
    print("\n" + "="*60)
    downloader.verify_downloads()


if __name__ == '__main__':
    main()

