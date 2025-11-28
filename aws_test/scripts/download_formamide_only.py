#!/usr/bin/env python3
"""
Download only formamide_extended results from AWS
"""
import os
import sys
import subprocess
from pathlib import Path

def fix_key_permissions(key_path):
    """Try to fix key permissions on Windows"""
    if sys.platform == 'win32':
        # Try using icacls with current user
        try:
            import getpass
            username = getpass.getuser()
            cmd = ['icacls', key_path, '/inheritance:r', f'/grant:r', f'{username}:(R)']
            subprocess.run(cmd, check=False, capture_output=True)
        except:
            pass

def download_formamide(host, key_path, local_dir="results/phase2b_additional"):
    """Download formamide_extended from AWS"""
    
    print(f"📥 Downloading formamide_extended from {host}")
    
    # Fix permissions
    fix_key_permissions(key_path)
    
    # Create local directory
    local_path = Path(local_dir) / "formamide_extended"
    local_path.mkdir(parents=True, exist_ok=True)
    
    # Download using scp
    remote_path = "~/live2.0/aws_test/results/phase2b_additional/formamide_extended"
    
    cmd = [
        "scp", "-r", "-i", key_path,
        "-o", "StrictHostKeyChecking=no",
        f"ubuntu@{host}:{remote_path}/*",
        str(local_path)
    ]
    
    print(f"🔗 Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✅ Download completed successfully")
        print(f"📁 Results saved to: {local_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Download failed")
        print(f"Error: {e.stderr}")
        print(f"\n💡 Try fixing key permissions manually:")
        print(f"   icacls \"{key_path}\" /inheritance:r /grant:r \"%USERNAME%:(R)\"")
        return False

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="63.178.224.65")
    parser.add_argument("--key", default=r"D:\OneDrive\Pulpit\live_aws_credentials\key-do-live.pem")
    parser.add_argument("--local-dir", default="results/phase2b_additional")
    args = parser.parse_args()
    
    download_formamide(args.host, args.key, args.local_dir)

